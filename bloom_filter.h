#pragma once

#include <stdbool.h>    //bool
#include <stdlib.h>     //calloc() free() size_t


typedef struct Bloom_Filter Bloom_Filter;

//==================== Memory Management ====================

Bloom_Filter* Bloom_Filter_Create(size_t n, double p);

void Bloom_Filter_Clear(Bloom_Filter* bf);

void Bloom_Filter_Destroy(Bloom_Filter* bf);


//==================== Operating Functions ====================

void Element_Mark(Bloom_Filter* bf, const char* str);

bool Element_Query(Bloom_Filter* bf, const char* str);

double bloom_filter_current_false_positive_rate(Bloom_Filter *bf);