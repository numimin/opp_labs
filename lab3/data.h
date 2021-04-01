#ifndef LAB1_DATA_H
#define LAB1_DATA_H

#include <stddef.h>
#include <mpi.h>

#define N 1000
#define M 12000
#define K 3000

typedef struct {
    int proc_rank;
    int proc_count;

    MPI_Comm global_comm;

    MPI_Comm row_comm;
    MPI_Comm col_comm;

    MPI_Datatype b_type;
    MPI_Datatype last_b_type;

    MPI_Datatype out_type;
    MPI_Datatype last_row_out_type;
    MPI_Datatype last_col_out_type;
    MPI_Datatype last_out_type;

    int coords[2]; //(row, col)
    int dimensions[2]; 

    size_t row_per_proc;
    size_t col_per_proc;

    size_t last_row_size;
    size_t last_col_size;

    size_t global_n;
    size_t global_m;
    size_t global_k;
} proc_meta;

typedef struct {
    double* mat_a;
    double* mat_b;
    double* mat_out;

    size_t mat_n;
    size_t mat_m;
    size_t mat_k;
} proc_data;

size_t get_row_size(const proc_meta* this, int row);
size_t get_col_size(const proc_meta* this, int col);

int get_row(const proc_meta* this, int rank);
int get_col(const proc_meta* this, int rank);
int get_srow(const proc_meta* this);
int get_scol(const proc_meta* this);

int get_row_idx(const proc_meta* this, int rank);
int get_col_idx(const proc_meta* this, int rank);

int get_row_count(const proc_meta* this);
int get_col_count(const proc_meta* this);

size_t get_n(const proc_meta* this);
size_t get_m(const proc_meta* this);
size_t get_k(const proc_meta* this);

int is_root(const proc_meta* this);
int get_root(const proc_meta* this);

int get_proc_count(const proc_meta* this);
int get_rank(const proc_meta* this);
int get_cart_rank(const proc_meta* this, int row, int col);
int get_row_rank(const proc_meta* this, int col);
int get_col_rank(const proc_meta* this, int row);

MPI_Datatype get_b_type(const proc_meta* this, int col);
MPI_Datatype get_out_type(const proc_meta* this, int rank);

MPI_Comm get_glob_comm(const proc_meta* this);
MPI_Comm get_srow_comm(const proc_meta* this);
MPI_Comm get_scol_comm(const proc_meta* this);

void free_pdata(proc_data* this);
int init_pmeta(proc_meta* this, size_t row_cnt, size_t col_cnt);
int free_pmeta(proc_meta* this);
void init_pdata(proc_data* pdata, const proc_meta* pmeta);

void debug(const proc_meta* this);

#endif // !LAB1_DATA_H
