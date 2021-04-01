#include "main.h"

#include <mpi.h>
#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double* par_mul_mv(double* y, const double* mat_pt, const double* x_pt, const proc_meta* pmeta) {
    double* glob_res = malloc(pmeta->global_size * sizeof *glob_res);
    mul_mv_nm(glob_res, mat_pt, x_pt, pmeta->local_size, pmeta->global_size);

    for (size_t rank = 0; rank < pmeta->proc_count; ++rank) {
        MPI_Reduce(&glob_res[get_pidx(pmeta, rank)], y,
                   get_lsize(pmeta, rank), MPI_DOUBLE, MPI_SUM, 
                   rank, MPI_COMM_WORLD);
    }

    free(glob_res);
    return y;
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

    par_mul_mv(pdata->b, pdata->mat, pdata->x, pmeta);

    for (size_t i = 0; i < n; ++i) {
        pdata->x[i] = 0;
    }
}
