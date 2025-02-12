#ifndef CUCKOO_H
#define CUCKOO_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include "murmurhash2.h"

// Defines whether 32bit or 64bit hash function will be used.
// It will be deprecated in RedisBloom 3.0.

#define CUCKOO_BKTSIZE 2
#define CUCKOO_NULLFP 0
// extern int globalCuckooHash64Bit;

typedef uint8_t CuckooFingerprint;
typedef uint64_t CuckooHash;
typedef uint8_t CuckooBucket[1];
typedef uint8_t MyCuckooBucket;

typedef struct {
    uint32_t numBuckets;
    uint8_t bucketSize;
    time_t expire_time;
    MyCuckooBucket *data;
} SubCF;

typedef struct {
    uint64_t numBuckets;
    uint64_t numItems;
    uint64_t numDeletes;
    uint16_t numFilters;
    uint16_t bucketSize;
    uint16_t maxIterations;
    uint16_t ttl;
    SubCF *filters;
    pthread_rwlock_t lock;
    pthread_t thread_rotate;
} CuckooFilter;


#define CUCKOO_GEN_HASH(s, n) MurmurHash64A_Bloom(s, n, 0)

/*
#define CUCKOO_GEN_HASH(s, n)                       \
            globalCuckooHash64Bit == 1 ?            \
                MurmurHash64A_Bloom(s, n, 0) :      \
                murmurhash2(s, n, 0)
*/
typedef struct {
    uint64_t i1;
    uint64_t i2;
    CuckooFingerprint fp;
} CuckooKey;

typedef enum {
    CuckooInsert_Inserted = 1,
    CuckooInsert_Exists = 0,
    CuckooInsert_NoSpace = -1,
    CuckooInsert_MemAllocFailed = -2
} CuckooInsertStatus;

int CuckooFilter_Init(CuckooFilter *filter, uint64_t capacity, uint16_t bucketSize,
                      uint16_t maxIterations, uint16_t ttl);
void CuckooFilter_Free(CuckooFilter *filter);
CuckooInsertStatus CuckooFilter_Insert(CuckooFilter *filter, CuckooHash hash);
int CuckooFilter_Check(CuckooFilter *filter, CuckooHash hash);
void CuckooFilter_Compact(CuckooFilter *filter, bool cont);
void CuckooFilter_GetInfo(const CuckooFilter *cf, CuckooHash hash, CuckooKey *out);
#endif