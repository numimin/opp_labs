#include "data.h"

#include <math.h>
#include <malloc.h>
#include <stdlib.h>

#include "matrix.h"

size_t get_lsize(const proc_meta* this, int rank) {
    return (rank < this->proc_count - 1) 
            ? this->data_per_proc
            : this->global_size - this->data_per_proc * (this->proc_count - 1);
}

size_t get_pidx(const proc_meta* this, int rank) {
    return this->data_per_proc * rank;
}

size_t get_spidx(const proc_meta* this) {
    return get_pidx(this, this->proc_rank);
}

void free_pdata(proc_data* this) {
    free(this->b);
    free(this->x);
    free(this->mat);
}

double fill_mat(size_t i, size_t j) {
    return (i == j) ? (-4) 
                    : (abs(i - j) == 1) ? 1
                    : (abs(i - j) == NX) ? 1
                    : 0;
}

void fill_pdata(proc_data* pdata, const proc_meta* pmeta) {
    const size_t n = pmeta->local_size;
    const size_t m = pmeta->global_size;

    pdata->b = malloc(n * sizeof *pdata->b);
    pdata->x = malloc(n * sizeof *pdata->x);
    for (size_t i = 0; i < n; ++i) {
        pdata->x[i] = sin(2 * M_PI * (i + get_spidx(pmeta)) / m);
    }

    pdata->mat = malloc(n * m * sizeof *pdata->mat);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            pdata->mat[n * i + j] = fill_mat(i, j + get_spidx(pmeta));
        }
    }

    mul_mv_s(pdata->b, pdata->mat, pdata->x, pmeta);

    for (size_t i = 0; i < n; ++i) {
        pdata->x[i] = 0;
    }
}
