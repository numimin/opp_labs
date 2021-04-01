#ifndef LAB1_DATA_H
#define LAB1_DATA_H

#include <stddef.h>

#define NX 100
#define NY 40
#define N (NX * NY)
#define E 0.000001

typedef struct {
    int proc_rank;
    int proc_count;

    size_t local_size;
    size_t global_size;
    size_t data_per_proc;
} proc_meta;

typedef struct {
    double* mat;
    double* b;
    double* x;
} proc_data;

size_t get_pidx(const proc_meta* this, int rank);
size_t get_spidx(const proc_meta* this);
void free_pdata(proc_data* this);
double fill_mat(size_t i, size_t j);
size_t get_lsize(const proc_meta* this, int rank);

#endif // !LAB1_DATA_H
