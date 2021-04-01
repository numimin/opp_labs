#include "data.h"

#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

size_t get_data_per_proc(size_t global_size, int proc_count) {
    return global_size / proc_count + (global_size % proc_count ? 1 : 0);
}

size_t get_last_cut_size(size_t global_size, int proc_count) {
    const size_t data_per_proc = get_data_per_proc(global_size, proc_count);
    return global_size - data_per_proc * (proc_count - 1);
}

size_t get_cut_size(size_t global_size, int proc_count, int rank) {
    return (rank < proc_count - 1) 
            ? get_data_per_proc(global_size, proc_count)
            : get_last_cut_size(global_size, proc_count);
}

size_t get_msize(const proc_meta* this) {
    return this->global_m;
}

size_t get_nsize(const proc_meta* this, int row) {
    return get_cut_size(this->global_n, this->dimensions[0], row);
}

size_t get_snsize(const proc_meta* this) {
    return get_nsize(this, this->coords[0]);
}

size_t get_ksize(const proc_meta* this, int col) {
    return get_cut_size(this->global_k, this->dimensions[1], col);
}

size_t get_sksize(const proc_meta* this) {
    return get_ksize(this, this->coords[1]);
}

void free_pdata(proc_data* this) {
    free(this->mat_a);
    free(this->mat_b);
    free(this->mat_out);
}   

int init_pmeta(proc_meta* this, size_t row_cnt, size_t col_cnt) {
    MPI_Comm_size(MPI_COMM_WORLD, &this->proc_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &this->proc_rank);

    if (col_cnt * row_cnt != this->proc_count) {
        fprintf(stderr, "Expected %d processes for [%d x %d], found %d\n", 
            col_cnt * row_cnt, row_cnt, col_cnt, this->proc_count);
        return -1;
    }

    this->dimensions[0] = row_cnt;
    this->dimensions[1] = col_cnt;

    this->global_n = N;
    this->global_m = M;
    this->global_k = K;

    this->row_per_proc = get_data_per_proc(this->global_n, row_cnt);
    this->col_per_proc = get_data_per_proc(this->global_k, col_cnt);
    this->last_row_size = get_last_cut_size(this->global_n, row_cnt);
    this->last_col_size = get_last_cut_size(this->global_k, col_cnt);

    MPI_Type_vector(this->global_m, this->col_per_proc, this->global_k, MPI_DOUBLE, &this->b_type);
    MPI_Type_commit(&this->b_type);

    MPI_Type_vector(this->global_m, this->last_col_size, this->global_k, MPI_DOUBLE, &this->last_b_type);
    MPI_Type_commit(&this->last_b_type);

    MPI_Type_vector(this->row_per_proc, this->col_per_proc, this->global_k, MPI_DOUBLE, &this->out_type);
    MPI_Type_commit(&this->out_type);

    MPI_Type_vector(this->row_per_proc, this->last_col_size, this->global_k, MPI_DOUBLE, &this->last_col_out_type);
    MPI_Type_commit(&this->last_col_out_type);

    MPI_Type_vector(this->last_row_size, this->col_per_proc, this->global_k, MPI_DOUBLE, &this->last_row_out_type);
    MPI_Type_commit(&this->last_row_out_type);

    MPI_Type_vector(this->last_row_size, this->last_col_size, this->global_k, MPI_DOUBLE, &this->last_out_type);
    MPI_Type_commit(&this->last_out_type);

    int periodic[2] = {};
    MPI_Cart_create(MPI_COMM_WORLD, 2, this->dimensions, periodic, 1, &this->global_comm);

    MPI_Cart_coords(this->global_comm, this->proc_rank, 2, this->coords);
    MPI_Cart_rank(this->global_comm, this->coords, &this->proc_rank);
    
    int row_comm_dims[] = {0, 1};
    MPI_Cart_sub(this->global_comm, row_comm_dims, &this->row_comm);

    int col_comm_dims[] = {1, 0};
    MPI_Cart_sub(this->global_comm, col_comm_dims, &this->col_comm);

    return 0;
}

void init_pdata(proc_data* pdata, const proc_meta* pmeta) {
    pdata->mat_n = get_snsize(pmeta);
    pdata->mat_m = get_msize(pmeta);
    pdata->mat_k = get_sksize(pmeta);

    pdata->mat_a = malloc(pdata->mat_n * pdata->mat_m * sizeof *pdata->mat_a);
    pdata->mat_b = malloc(pdata->mat_m * pdata->mat_k * sizeof *pdata->mat_b);
    pdata->mat_out = malloc(pdata->mat_n * pdata->mat_k * sizeof *pdata->mat_out);
}

void debug(const proc_meta* this) {
    fprintf(stderr, "Coords: [%d, %d], Dims: [%d, %d]\n", get_srow(this), get_scol(this),
        get_row_count(this), get_col_count(this));

    fprintf(stderr, "[n x m x k] == [%d x %d x %d]\n", this->global_n, this->global_m, this->global_k);

    fprintf(stderr, "Per process: [%d, %d]\n", this->row_per_proc, this->col_per_proc);
    fprintf(stderr, "Last process: [%d, %d]\n", this->last_row_size, this->last_col_size);
}

int get_cart_rank(const proc_meta* this, int row, int col) {
    int rank;
    int coords[] = {row, col};
    MPI_Cart_rank(this->global_comm, coords, &rank);
    return rank;
}

int get_row_rank(const proc_meta* this, int col) {
    int rank;
    MPI_Cart_rank(this->row_comm, &col, &rank);
    return rank;
}

int get_col_rank(const proc_meta* this, int row) {
    int rank;
    MPI_Cart_rank(this->col_comm, &row, &rank);
    return rank;
}

size_t get_row_size(const proc_meta* this, int row) {
    if (row + 1 == get_row_count(this)) {
        return this->last_row_size;
    } else {
        return this->row_per_proc;
    }
}

size_t get_col_size(const proc_meta* this, int col) {
    if (col + 1 == get_col_count(this)) {
        return this->last_col_size;
    } else {
        return this->col_per_proc;
    }
}

int get_row(const proc_meta* this, int rank) {
    int coords[2];
    MPI_Cart_coords(this->global_comm, rank, 2, coords);
    return coords[0];
}

int get_col(const proc_meta* this, int rank) {
    int coords[2];
    MPI_Cart_coords(this->global_comm, rank, 2, coords);
    return coords[1];
}

int get_srow(const proc_meta* this) {
    return this->coords[0];
}

int get_scol(const proc_meta* this) {
    return this->coords[1];
}

int get_row_count(const proc_meta* this) {
    return this->dimensions[0];
}

int get_col_count(const proc_meta* this) {
    return this->dimensions[1];
}

MPI_Datatype get_b_type(const proc_meta* this, int col) {
    if (col + 1 == get_col_count(this)) {
        return this->last_b_type;
    } else {
        return this->b_type;
    }
}

MPI_Datatype get_out_type(const proc_meta* this, int rank) {
    const int col = get_col(this, rank);
    const int row = get_row(this, rank);

    if (col + 1 == get_col_count(this)) {
        if (row + 1 == get_row_count(this)) {
            return this->last_out_type;
        }

        return this->last_col_out_type;
    }

    if (row + 1 == get_row_count(this)) {
        return this->last_row_out_type;
    }

    return this->out_type;
}

int free_pmeta(proc_meta* this) {
    MPI_Comm_free(&this->global_comm);
    MPI_Comm_free(&this->row_comm);
    MPI_Comm_free(&this->col_comm);

    MPI_Type_free(&this->b_type);
    MPI_Type_free(&this->last_b_type);

    MPI_Type_free(&this->last_out_type);
    MPI_Type_free(&this->last_row_out_type);
    MPI_Type_free(&this->last_col_out_type);
    MPI_Type_free(&this->out_type);
}

int is_root(const proc_meta* this) {
    return this->proc_rank == get_cart_rank(this, 0, 0);
}

int get_row_idx(const proc_meta* this, int rank) {
    return this->row_per_proc * get_row(this, rank);
}

int get_col_idx(const proc_meta* this, int rank) {
    return this->col_per_proc * get_col(this, rank);
}

size_t get_n(const proc_meta* this) {
    return this->global_n;
}

size_t get_m(const proc_meta* this) {
    return this->global_m;
}

size_t get_k(const proc_meta* this) {
    return this->global_k;
}

int get_proc_count(const proc_meta* this) {
    return this->proc_count;
}

int get_rank(const proc_meta* this) {
    return this->proc_rank;
}

MPI_Comm get_glob_comm(const proc_meta* this) {
    return this->global_comm;
}

MPI_Comm get_srow_comm(const proc_meta* this) {
    return this->row_comm;
}

MPI_Comm get_scol_comm(const proc_meta* this) {
    return this->col_comm;
}

int get_root(const proc_meta* this) {
    return get_cart_rank(this, 0, 0);
}
