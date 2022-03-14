#include "cuckoo.h"
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <unistd.h>

#ifndef CUCKOO_MALLOC
#define CUCKOO_MALLOC malloc
#define CUCKOO_CALLOC calloc
#define CUCKOO_REALLOC realloc
#define CUCKOO_FREE free
#endif

// int globalCuckooHash64Bit;


void *rotate(void *arg){
    CuckooFilter *filter = (CuckooFilter *)arg;
    sleep(filter->ttl*2);

    while (true)
    {
        pthread_rwlock_wrlock(&lock);

        time_t now = time(0);
        int del_count = 0;
        for(int i=0; i < filter->numFilters; i++){
            if (filter->filters[i].expire_time < now){
                del_count +=1;
                CUCKOO_FREE(filter->filters[i].data);
            }else{
                break;
            }
        }
        if(del_count > 0){
            uint16_t new_len = filter->numFilters - del_count;
            if (new_len == 0){
                CUCKOO_FREE(filter->filters);
                    SubCF *filtersArray = CUCKOO_REALLOC(filter->filters, sizeof(*filtersArray));
                    SubCF *currentFilter = filtersArray;
                    currentFilter->bucketSize = filter->bucketSize;
                    currentFilter->numBuckets = filter->numBuckets;
                    currentFilter->expire_time = time(0) + filter->ttl;
                    currentFilter->data =
                    CUCKOO_CALLOC((size_t)currentFilter->numBuckets * filter->bucketSize, sizeof(CuckooBucket));
                    filter->numFilters = 1;
                    filter->filters = filtersArray;
            }else {
                SubCF *filtersArray = CUCKOO_REALLOC(filter->filters, sizeof(*filtersArray) * (new_len));
                for(int i=del_count; i < filter->numFilters; i++){
                    filtersArray[i-del_count].data = filter->filters[i].data;
                    filtersArray[i-del_count].bucketSize = filter->filters[i].bucketSize;
                    filtersArray[i-del_count].expire_time = filter->filters[i].expire_time;
                    filtersArray[i-del_count].numBuckets = filter->filters[i].numBuckets;
                }
                CUCKOO_FREE(filter->filters);
                filter->filters = filtersArray;
                filter->numFilters = new_len;
            }
        }

        pthread_rwlock_unlock(&lock);

        sleep(filter->ttl);
        sleep(10);
    }

  return NULL;
}

static int CuckooFilter_Grow(CuckooFilter *filter);

static int isPower2(uint64_t num) { return (num & (num - 1)) == 0 && num != 0; }

static uint64_t getNextN2(uint64_t n) {
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n |= n >> 32;
    n++;
    return n;
}

int CuckooFilter_Init(CuckooFilter *filter, uint64_t capacity, uint16_t bucketSize,
                      uint16_t maxIterations, uint16_t ttl) {
    memset(filter, 0, sizeof(*filter));
    filter->ttl = getNextN2(ttl);
    filter->bucketSize = bucketSize;
    filter->maxIterations = maxIterations;
    filter->numBuckets = getNextN2(capacity / bucketSize);
    if (filter->numBuckets == 0) {
        filter->numBuckets = 1;
    }
    assert(isPower2(filter->numBuckets));

    if (CuckooFilter_Grow(filter) != 0) {
        return -1; // LCOV_EXCL_LINE memory failure
    }
    pthread_rwlock_init(&lock, NULL);

	pthread_t thread_rotate;
	pthread_create(&thread_rotate, NULL, rotate, filter);
    return 0;
}

void CuckooFilter_Free(CuckooFilter *filter) {
    pthread_rwlock_wrlock(&lock);
    for (uint16_t ii = 0; ii < filter->numFilters; ++ii) {
        CUCKOO_FREE(filter->filters[ii].data);
    }
    CUCKOO_FREE(filter->filters);
    pthread_rwlock_unlock(&lock);
    pthread_rwlock_destroy(&lock);
}

static int CuckooFilter_Grow(CuckooFilter *filter) {
    pthread_rwlock_wrlock(&lock);
    SubCF *filtersArray =
        CUCKOO_REALLOC(filter->filters, sizeof(*filtersArray) * (filter->numFilters + 1));

    if (!filtersArray) {
        pthread_rwlock_unlock(&lock);
        return -1; // LCOV_EXCL_LINE memory failure
    }

    SubCF *currentFilter = filtersArray + filter->numFilters;
    currentFilter->bucketSize = filter->bucketSize;
    currentFilter->numBuckets = filter->numBuckets;
    currentFilter->expire_time = time(0) + filter->ttl;
    currentFilter->data =
        CUCKOO_CALLOC((size_t)currentFilter->numBuckets * filter->bucketSize, sizeof(CuckooBucket));
    if (!currentFilter->data) {
        pthread_rwlock_unlock(&lock);
        return -1; // LCOV_EXCL_LINE memory failure
    }

    filter->numFilters++;
    filter->filters = filtersArray;
    pthread_rwlock_unlock(&lock);
    return 0;
}

typedef struct {
    CuckooHash h1;
    CuckooHash h2;
    CuckooFingerprint fp;
} LookupParams;

static CuckooHash getAltHash(CuckooFingerprint fp, CuckooHash index) {
    return ((CuckooHash)(index ^ ((CuckooHash)fp * 0x5bd1e995)));
}

static void getLookupParams(CuckooHash hash, LookupParams *params) {
    params->fp = hash % 255 + 1;

    params->h1 = hash;
    params->h2 = getAltHash(params->fp, params->h1);
    // assert(getAltHash(params->fp, params->h2, numBuckets) == params->h1);
}

static uint32_t SubCF_GetIndex(const SubCF *subCF, CuckooHash hash) {
    return (hash % subCF->numBuckets) * subCF->bucketSize;
}

static uint8_t *Bucket_Find(CuckooBucket bucket, uint16_t bucketSize, CuckooFingerprint fp) {
    for (uint16_t ii = 0; ii < bucketSize; ++ii) {
        if (bucket[ii] == fp) {
            return bucket + ii;
        }
    }
    return NULL;
}

static int Filter_Find(const SubCF *filter, const LookupParams *params) {
    uint8_t bucketSize = filter->bucketSize;
    uint64_t loc1 = SubCF_GetIndex(filter, params->h1);
    uint64_t loc2 = SubCF_GetIndex(filter, params->h2);
    return Bucket_Find(&filter->data[loc1], bucketSize, params->fp) != NULL ||
           Bucket_Find(&filter->data[loc2], bucketSize, params->fp) != NULL;
}

static int Bucket_Delete(CuckooBucket bucket, uint16_t bucketSize, CuckooFingerprint fp) {
    for (uint16_t ii = 0; ii < bucketSize; ii++) {
        if (bucket[ii] == fp) {
            bucket[ii] = CUCKOO_NULLFP;
            return 1;
        }
    }
    return 0;
}

static int Filter_Delete(const SubCF *filter, const LookupParams *params) {
    uint8_t bucketSize = filter->bucketSize;
    uint64_t loc1 = SubCF_GetIndex(filter, params->h1);
    uint64_t loc2 = SubCF_GetIndex(filter, params->h2);
    return Bucket_Delete(&filter->data[loc1], bucketSize, params->fp) ||
           Bucket_Delete(&filter->data[loc2], bucketSize, params->fp);
}

static int CuckooFilter_CheckFP(const CuckooFilter *filter, const LookupParams *params) {
    pthread_rwlock_rdlock (&lock);
    for (uint16_t ii = 0; ii < filter->numFilters; ++ii) {
        if (Filter_Find(&filter->filters[ii], params)) {
            pthread_rwlock_unlock(&lock);
            return 1;
        }
    }
    pthread_rwlock_unlock(&lock);
    return 0;
}

int CuckooFilter_Check(const CuckooFilter *filter, CuckooHash hash) {
    LookupParams params;
    getLookupParams(hash, &params);
    return CuckooFilter_CheckFP(filter, &params);
}

static uint16_t bucketCount(const CuckooBucket bucket, uint16_t bucketSize, CuckooFingerprint fp) {
    uint16_t ret = 0;
    for (uint16_t ii = 0; ii < bucketSize; ++ii) {
        if (bucket[ii] == fp) {
            ret++;
        }
    }
    return ret;
}

static uint64_t subFilterCount(const SubCF *filter, const LookupParams *params) {
    uint8_t bucketSize = filter->bucketSize;
    uint64_t loc1 = SubCF_GetIndex(filter, params->h1);
    uint64_t loc2 = SubCF_GetIndex(filter, params->h2);

    return bucketCount(&filter->data[loc1], bucketSize, params->fp) +
           bucketCount(&filter->data[loc2], bucketSize, params->fp);
}

static uint8_t *Bucket_FindAvailable(CuckooBucket bucket, uint16_t bucketSize) {
    for (uint16_t ii = 0; ii < bucketSize; ++ii) {
        if (bucket[ii] == CUCKOO_NULLFP) {
            return &bucket[ii];
        }
    }
    return NULL;
}

static uint8_t *Filter_FindAvailable(SubCF *filter, const LookupParams *params) {
    uint8_t *slot;
    uint8_t bucketSize = filter->bucketSize;
    uint64_t loc1 = SubCF_GetIndex(filter, params->h1);
    uint64_t loc2 = SubCF_GetIndex(filter, params->h2);
    if ((slot = Bucket_FindAvailable(&filter->data[loc1], bucketSize)) ||
        (slot = Bucket_FindAvailable(&filter->data[loc2], bucketSize))) {
        return slot;
    }
    return NULL;
}

static CuckooInsertStatus Filter_KOInsert(CuckooFilter *filter, SubCF *curFilter,
                                          const LookupParams *params);

static CuckooInsertStatus CuckooFilter_InsertFP(CuckooFilter *filter, const LookupParams *params) {
    pthread_rwlock_wrlock(&lock);
    for (uint16_t ii = filter->numFilters; ii > 0; --ii) {
        uint8_t *slot = Filter_FindAvailable(&filter->filters[ii - 1], params);
        if (slot) {
            *slot = params->fp;
            filter->numItems++;
            pthread_rwlock_unlock(&lock);
            return CuckooInsert_Inserted;
        }
    }

    // No space. Time to evict!
    CuckooInsertStatus status =
        Filter_KOInsert(filter, &filter->filters[filter->numFilters - 1], params);
    if (status == CuckooInsert_Inserted) {
        filter->numItems++;
        pthread_rwlock_unlock(&lock);
        return CuckooInsert_Inserted;
    }

    pthread_rwlock_unlock(&lock);

    if (CuckooFilter_Grow(filter) != 0) {
        return CuckooInsert_MemAllocFailed;
    }

    // Try to insert the filter again
    return CuckooFilter_InsertFP(filter, params);
}

CuckooInsertStatus CuckooFilter_Insert(CuckooFilter *filter, CuckooHash hash) {
    LookupParams params;
    getLookupParams(hash, &params);
    if (CuckooFilter_CheckFP(filter, &params)) {
        return CuckooInsert_Exists;
    }
    return CuckooFilter_InsertFP(filter, &params);
}

static void swapFPs(uint8_t *a, uint8_t *b) {
    uint8_t temp = *a;
    *a = *b;
    *b = temp;
}

static CuckooInsertStatus Filter_KOInsert(CuckooFilter *filter, SubCF *curFilter,
                                          const LookupParams *params) {
    uint16_t maxIterations = filter->maxIterations;
    uint32_t numBuckets = curFilter->numBuckets;
    uint16_t bucketSize = filter->bucketSize;
    CuckooFingerprint fp = params->fp;

    uint16_t counter = 0;
    uint32_t victimIx = 0;
    uint32_t ii = params->h1 % numBuckets;

    while (counter++ < maxIterations) {
        uint8_t *bucket = &curFilter->data[ii * bucketSize];
        swapFPs(bucket + victimIx, &fp);
        ii = getAltHash(fp, ii) % numBuckets;
        // Insert the new item in potentially the same bucket
        uint8_t *empty = Bucket_FindAvailable(&curFilter->data[ii * bucketSize], bucketSize);
        if (empty) {
            // printf("Found slot. Bucket[%lu], Pos=%lu\n", ii, empty - curFilter[ii]);
            // printf("Old FP Value: %d\n", *empty);
            // printf("Setting FP: %p\n", empty);
            *empty = fp;
            return CuckooInsert_Inserted;
        }
        victimIx = (victimIx + 1) % bucketSize;
    }

    // If we weren't able to insert, we roll back and try to insert new element in new filter
    counter = 0;
    while (counter++ < maxIterations) {
        victimIx = (victimIx + bucketSize - 1) % bucketSize;
        ii = getAltHash(fp, ii) % numBuckets;
        uint8_t *bucket = &curFilter->data[ii * bucketSize];
        swapFPs(bucket + victimIx, &fp);
    }

    return CuckooInsert_NoSpace;
}

#define RELOC_EMPTY 0
#define RELOC_OK 1
#define RELOC_FAIL -1

/**
 * Attempt to move a slot from one bucket to another filter
 */
static int relocateSlot(CuckooFilter *cf, CuckooBucket bucket, uint16_t filterIx, uint64_t bucketIx,
                        uint16_t slotIx) {
    LookupParams params = {0};
    if ((params.fp = bucket[slotIx]) == CUCKOO_NULLFP) {
        // Nothing in this slot.
        return RELOC_EMPTY;
    }

    // Because We try to insert in sub filter with less or equal number of
    // buckets, our current fingerprint is sufficient
    params.h1 = bucketIx;
    params.h2 = getAltHash(params.fp, bucketIx);

    // Look at all the prior filters and attempt to find a home
    for (uint16_t ii = 0; ii < filterIx; ++ii) {
        uint8_t *slot = Filter_FindAvailable(&cf->filters[ii], &params);
        if (slot) {
            *slot = params.fp;
            bucket[slotIx] = CUCKOO_NULLFP;
            return RELOC_OK;
        }
    }
    return RELOC_FAIL;
}

/**
 * Attempt to strip a single filter moving it down a slot
 */
static uint64_t CuckooFilter_CompactSingle(CuckooFilter *cf, uint16_t filterIx) {
    SubCF *currentFilter = &cf->filters[filterIx];
    MyCuckooBucket *filter = currentFilter->data;
    int rv = RELOC_OK;

    for (uint64_t bucketIx = 0; bucketIx < currentFilter->numBuckets; ++bucketIx) {
        for (uint16_t slotIx = 0; slotIx < currentFilter->bucketSize; ++slotIx) {
            int status = relocateSlot(cf, &filter[bucketIx * currentFilter->bucketSize], filterIx,
                                      bucketIx, slotIx);
            if (status == RELOC_FAIL) {
                rv = RELOC_FAIL;
            }
        }
    }
    // we free a filter only if it the latest one
    if (rv == RELOC_OK && filterIx == cf->numFilters - 1) {
        CUCKOO_FREE(filter);
        cf->numFilters--;
    }
    return rv;
}

/**
 * Attempt to move elements to older filters. If latest filter is emptied, it is freed.
 * `bool` determines whether to continue iteration on other filters once a filter cannot
 * be freed and therefore following filter cannot be freed either.
 */
void CuckooFilter_Compact(CuckooFilter *cf, bool cont) {
    for (uint64_t ii = cf->numFilters; ii > 1; --ii) {
        if (CuckooFilter_CompactSingle(cf, ii - 1) == RELOC_FAIL && !cont) {
            // if compacting failed, stop as lower filters cannot be freed.
            break;
        }
    }
    cf->numDeletes = 0;
}

/* CF.DEBUG uses another function
void CuckooFilter_GetInfo(const CuckooFilter *cf, CuckooHash hash, CuckooKey *out) {
    LookupParams params;
    getLookupParams(hash, cf->numBuckets, &params);
    out->fp = params.fp;
    out->h1 = params.h1;
    out->h2 = params.h2;
    assert(getAltHash(params.fp, out->h1, cf->numBuckets) == out->h2);
    assert(getAltHash(params.fp, out->h2, cf->numBuckets) == out->h1);
}*/
