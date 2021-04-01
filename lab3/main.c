#include <mpi.h>
#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "data.h"

//a[n x m]
//b[m x k]
//out[n x k]
double* mul_mm(const double* a, const double* b, double* out, size_t n, size_t m, size_t k) {
    for (size_t i = 0; i < n; ++i) {
        double* const out_row = &out[i * k];
        memset(out_row, 0, sizeof *out_row * k);

        for (size_t j = 0; j < m; ++j) {
            const double factor = a[i * m + j];
            const double* const b_row = &b[j * k];
            for (size_t p = 0; p < k; ++p) {
                out_row[p] += factor * b_row[p];
            }
        }
    }
}

void fill_mat(double* mat, size_t n, size_t m) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            mat[m * i + j] = 123.0 + (i - 50.0) * (j + 123.0) * j;
        }
    }
}

void print_mat(const double* mat, size_t n, size_t m) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            printf("%f ", mat[i * m + j]);
        }
        printf("\n");
    }
}

void scatter_matrices(proc_meta* pmeta, proc_data* pdata, int tag) {
    double* global_a = malloc(get_n(pmeta) * get_m(pmeta) * sizeof *global_a);
    double* global_b = malloc(get_m(pmeta) * get_k(pmeta) * sizeof *global_b);

    fill_mat(global_a, get_n(pmeta), get_m(pmeta));
    fill_mat(global_b, get_m(pmeta), get_k(pmeta));

    for (size_t row = 1; row < get_row_count(pmeta); ++row) {
        const int rank = get_cart_rank(pmeta, row, 0);

        MPI_Send(&global_a[get_row_idx(pmeta, rank) * get_m(pmeta)], 
            get_row_size(pmeta, row) * get_m(pmeta), MPI_DOUBLE, 
            rank, tag, get_glob_comm(pmeta));
    }

    for (size_t col = 1; col < get_col_count(pmeta); ++col) {
        const int rank = get_cart_rank(pmeta, 0, col);

        MPI_Send(&global_b[get_col_idx(pmeta, rank)], 
            1, get_b_type(pmeta, col), 
            rank, tag, get_glob_comm(pmeta));
    }

    for (size_t i = 0; i < pdata->mat_n; ++i) {
        memcpy(&pdata->mat_a[i * pdata->mat_m], &global_a[i * get_m(pmeta)], 
            pdata->mat_m * sizeof *pdata->mat_a);
    }

    for (size_t i = 0; i < pdata->mat_m; ++i) {
        memcpy(&pdata->mat_b[i * pdata->mat_k], &global_b[i * get_k(pmeta)], 
            pdata->mat_k * sizeof *pdata->mat_b);
    }

    free(global_a);
    free(global_b);
}

void gather_matrices(proc_meta* pmeta, proc_data* pdata, int tag) {
    if (get_scol(pmeta) == 0) {
        MPI_Recv(pdata->mat_a, pdata->mat_n * pdata->mat_m, MPI_DOUBLE, 
            get_root(pmeta), tag, get_glob_comm(pmeta), MPI_STATUS_IGNORE);
    } 

    if (get_srow(pmeta) == 0) {
        MPI_Recv(pdata->mat_b, pdata->mat_m * pdata->mat_k, MPI_DOUBLE, 
            get_root(pmeta), tag, get_glob_comm(pmeta), MPI_STATUS_IGNORE);
    }
}

void gather_pdata(proc_data* pdata, proc_meta* pmeta) {
    const int tag = 123;

    if (is_root(pmeta)) {
        scatter_matrices(pmeta, pdata, tag);
    } else {
        gather_matrices(pmeta, pdata, tag);
    }

    MPI_Bcast(pdata->mat_a, pdata->mat_n * pdata->mat_m, MPI_DOUBLE,
        get_row_rank(pmeta, 0), get_srow_comm(pmeta));

    MPI_Bcast(pdata->mat_b, pdata->mat_m * pdata->mat_k, MPI_DOUBLE,
        get_col_rank(pmeta, 0), get_scol_comm(pmeta));
}

void par_mul_mm(proc_meta* pmeta) {
    proc_data pdata;
    init_pdata(&pdata, pmeta);
    gather_pdata(&pdata, pmeta);
    mul_mm(pdata.mat_a, pdata.mat_b, pdata.mat_out, pdata.mat_n, pdata.mat_m, pdata.mat_k);

    const int tag = 234;
    if (is_root(pmeta)) {
        double* out = malloc(get_n(pmeta) * get_k(pmeta) * sizeof *out);

        for (size_t i = 0; i < pdata.mat_n; ++i) {
            memcpy(&out[i * get_k(pmeta)], &pdata.mat_out[i * pdata.mat_k], 
                pdata.mat_k * sizeof *out);
        }
        free_pdata(&pdata);

        for (size_t rank = 0; rank < get_proc_count(pmeta); ++rank) {
            if (rank == get_cart_rank(pmeta, 0, 0)) continue;

            MPI_Recv(&out[get_k(pmeta) * get_row_idx(pmeta, rank) + get_col_idx(pmeta, rank)],
                1, get_out_type(pmeta, rank),
                rank, tag, get_glob_comm(pmeta), MPI_STATUS_IGNORE);
        }

        free(out);
        return;
    }

    MPI_Send(pdata.mat_out, 
        pdata.mat_n * pdata.mat_k, MPI_DOUBLE, 
        get_cart_rank(pmeta, 0, 0), tag, get_glob_comm(pmeta));
    free_pdata(&pdata);
}

void test_mul(size_t row_count, size_t col_count) {
    proc_meta pmeta;
    init_pmeta(&pmeta, row_count, col_count);

    const double start = MPI_Wtime();
    par_mul_mm(&pmeta);
    const double end = MPI_Wtime();

    if (get_rank(&pmeta) == 0) {
        printf("Cores %d [rows: %d | cols: %d]: %f\n", 
            get_proc_count(&pmeta), row_count, col_count, end - start);
    }

    free_pmeta(&pmeta);
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int proc_count;
    MPI_Comm_size(MPI_COMM_WORLD, &proc_count);

    for (size_t cols = 1; cols <= proc_count; ++cols) {
        if (proc_count % cols) continue;

        const size_t rows = proc_count / cols;
        test_mul(rows, cols);
    }


    MPI_Finalize();
    return 0;
}
