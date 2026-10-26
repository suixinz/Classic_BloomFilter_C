
#include "bloom_filter.h"


#include <math.h>       //log() ceil() pow() exp()
#include <string.h>     //memset()

#ifndef BITS_PER_BYTE

#define BITS_PER_BYTE 8

#endif

struct Bloom_Filter{

    unsigned char* bit_array;
    size_t bit_size;

    size_t hash_count;

    size_t ele_count;

};

/*

n-->​	Expected number of elements​ to be inserted into the filter
p-->​	Target false positive rate​ (acceptable probability of false positives, between 0 and 1)
m​-->	Size of the bit array​ (total number of bits)
k​-->	Number of hash functions​ used per element

relation: n, p ==> m, k

m = - n * ln(p) / (ln(2))^2;
k = - ln(p) / ln(2);

*/

//==================== Memory Management ====================

static const double ln_2 = 0.693147;

static inline size_t align_up(size_t n){
    return (n + BITS_PER_BYTE - 1) & (~(size_t)(BITS_PER_BYTE - 1));
}

static inline size_t max(size_t a, size_t b){
    return a > b ? a : b;
}

Bloom_Filter* Bloom_Filter_Create(size_t n, double p){

    Bloom_Filter* bf = NULL;

    if(n == 0 || p <= 0.0 || p >= 1.0){
        goto cleanup;
    }

    bf = (Bloom_Filter*)calloc(1, sizeof(Bloom_Filter));

    if(bf == NULL){
        goto cleanup;
    }

    const double ln_p = log(p);

    const size_t m = (size_t)ceil(- (double)n * ln_p / (ln_2 * ln_2));

    /*
    * WARNING: Applying unary minus to a size_t value does NOT produce
    * a negative number due to unsigned integer wrap-around.
    *
    * For example:
    *   size_t n = 10000;
    *   -n  →  SIZE_MAX - 10000 + 1  (a huge positive value on 64-bit)
    *
    * Always cast to double first before negation:
    *   double result = -(double)n;   // Correct: -10000.0
    * 2026.10.26
    */

    const size_t k = max((size_t)ceil(- ln_p / ln_2), 1);

    bf->bit_size = align_up(m);
    bf->bit_array = (unsigned char*)calloc(bf->bit_size / BITS_PER_BYTE, sizeof(unsigned char));

    if(bf->bit_array == NULL){
        goto cleanup;
    }

    bf->hash_count = k;
    bf->ele_count = 0;

    return bf;


    cleanup:

        free(bf);
        return NULL;

}

void Bloom_Filter_Clear(Bloom_Filter* bf){

    if(bf == NULL){
        return;
    }

    memset(bf->bit_array, 0, bf->bit_size / BITS_PER_BYTE);

    bf->ele_count = 0;

}

void Bloom_Filter_Destroy(Bloom_Filter* bf){

    if(bf == NULL){
        return;
    }

    free(bf->bit_array);
    free(bf);

}


//==================== Operating Functions ====================

static size_t BKDRHash(const char *str){

    size_t seed = 131;
    size_t hash = 0;

    while(*str){
        hash = hash * seed + (unsigned char)(*str++);
    }

    return hash;

}

static size_t DJBHash(const char *str){

    size_t hash = 5381;

    while(*str){
        hash = ((hash << 5) + hash) + (unsigned char)(*str++);
    }

    return hash;

}


void Element_Mark(Bloom_Filter* bf, const char* str){

    const size_t bkdr_hash = BKDRHash(str);
    const size_t djb_hash  = DJBHash(str);

    for(size_t i = 0; i < bf->hash_count; i++){

        size_t pos_i = (bkdr_hash + i * djb_hash) % bf->bit_size;

        bf->bit_array[pos_i / BITS_PER_BYTE] |= (1 << (pos_i % BITS_PER_BYTE));

    }

    bf->ele_count++;

}

bool Element_Query(Bloom_Filter* bf, const char* str){

    const size_t bkdr_hash = BKDRHash(str);
    const size_t djb_hash  = DJBHash(str);

    for(size_t i = 0; i < bf->hash_count; i++){

        size_t pos_i = (bkdr_hash + i * djb_hash) % bf->bit_size;

        if(bf->bit_array[pos_i / BITS_PER_BYTE] & (1 << (pos_i % BITS_PER_BYTE))){
            continue;
        }

        return false;

    }

    return true;

}

double bloom_filter_current_false_positive_rate(Bloom_Filter *bf){

    // Formula: p = (1 - e^(-kn/m))^k
    // k: number of hash functions
    // n: number of inserted elements
    // m: size of bit array

    double k = bf->hash_count;
    double n = bf->ele_count;
    double m = bf->bit_size;
    
    double exponent = -k * n / m;
    double p = pow(1 - exp(exponent), k);
    return p;

}