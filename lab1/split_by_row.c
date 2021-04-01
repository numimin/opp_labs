#include "main.h"

#include <mpi.h>
#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double* par_mul_mv(double* y, const double* mat_pt, const double* x_pt, const proc_meta* pmeta) {
    double* x = malloc(pmeta->data_per_proc * pmeta->proc_count * sizeof *x);
    MPI_Allgather(x_pt, pmeta->data_per_proc, MPI_DOUBLE, x, pmeta->data_per_proc, MPI_DOUBLE, MPI_COMM_WORLD);
    mul_mv_nm(y, mat_pt, x, pmeta->global_size, pmeta->local_size);
    free(x);
    return y;
}

void fill_pdata(proc_data* pdata, const proc_meta* pmeta) {
    const size_t n = pmeta->local_size;
    const size_t m = pmeta->global_size;

    pdata->b = malloc(n * sizeof *pdata->b);
    pdata->x = malloc(pmeta->data_per_proc * sizeof *pdata->x);
    for (size_t i = 0; i < n; ++i) {
        pdata->x[i] = sin(2 * M_PI * (i + get_spidx(pmeta)) / m);
    }

    pdata->mat = malloc(n * m * sizeof *pdata->mat);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            pdata->mat[m * i + j] = fill_mat(i + get_spidx(pmeta), j);
        }
    }

    par_mul_mv(pdata->b, pdata->mat, pdata->x, pmeta);

    for (size_t i = 0; i < n; ++i) {
        pdata->x[i] = 0;
    }
}
