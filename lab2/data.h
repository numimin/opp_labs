#ifndef LAB1_DATA_H
#define LAB1_DATA_H

#include <stddef.h>

#define NX 100
#define NY 40
#define N (NX * NY)
#define E 0.000001

typedef struct {
    double* mat;
    double* b;
    double* x;
} proc_data;

void fill_pdata(proc_data* pdata, size_t n);
void free_pdata(proc_data* this);

#endif // !LAB1_DATA_H
