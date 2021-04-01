#include "data.h"

#include <math.h>
#include <malloc.h>
#include <stdlib.h>

#include "matrix.h"

double fill_mat(size_t i, size_t j) {
    return (i == j) ? (-4) 
                    : (abs(i - j) == 1) ? 1
                    : (abs(i - j) == NX) ? 1
                    : 0;
}

void fill_pdata(proc_data* pdata, size_t n) {
    pdata->b = malloc(n * sizeof *pdata->b);
    pdata->x = malloc(n * sizeof *pdata->x);
    for (size_t i = 0; i < n; ++i) {
        pdata->x[i] = sin(2 * M_PI * i / n);
    }

    pdata->mat = malloc(n * n * sizeof *pdata->mat);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            pdata->mat[n * i + j] = fill_mat(i, j);
        }
    }

    mul_mv_s(pdata->b, pdata->mat, pdata->x, n);

    for (size_t i = 0; i < n; ++i) {
        pdata->x[i] = 0;
    }
}

void free_pdata(proc_data* this) {
    free(this->b);
    free(this->x);
    free(this->mat);
}
