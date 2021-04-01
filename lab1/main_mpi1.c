#include <mpi.h>
#include <malloc.h>
#include <string.h>
#include <math.h>

#include "matrix.h"


/*unsigned ulog2(unsigned x) {
    if (x == 0) return 0;

    unsigned power = 0;
    while (x >>= 1) ++power;
    return power;
}*/

/*double* par_sum_pyramid(double* res_vec, const double* local_vec, 
                        size_t vec_size, size_t height, int index) {
    
    const size_t proc_count = 1 << height;
    
    memcpy(res_vec, local_vec, vec_size * sizeof *res_vec);

    double* buf_vec = malloc(vec_size * sizeof *buf_vec);
    const int pyramid_tag = 123;

    for (size_t i = 0; i < height; ++i) {
        const size_t offset = 1 << i;
        const size_t group_count = offset;

        if ((index / group_count) % 2 == 0) {
            MPI_Sendrecv(res_vec, vec_size, MPI_DOUBLE, (index + offset) % proc_count, pyramid_tag, 
                         buf_vec, vec_size, MPI_DOUBLE, (index + offset) % proc_count, pyramid_tag,
                         MPI_COMM_WORLD, NULL);
        } else {
            MPI_Sendrecv(res_vec, vec_size, MPI_DOUBLE, (index - offset) % proc_count, pyramid_tag, 
                         buf_vec, vec_size, MPI_DOUBLE, (index - offset) % proc_count, pyramid_tag,
                         MPI_COMM_WORLD, NULL);
        }
        for (size_t i = 0; i < vec_size; ++i) {
            res_vec[i] += buf_vec[i];
        }
    }
    free(buf_vec);

    return res_vec;
}*/

/*double dot_product(const double* a, const double* b, size_t n) {
    double res = 0;
    for (size_t i = 0; i < n; ++i) {
        res += a[i] * b[i];
    }
    return res;
}

double* mul_mv_nm(double* y, const double* mat, const double* x, size_t n, size_t m) {
    for (size_t i = 0; i < m; ++i) {
        y[i] = dot_product(&mat[m * i], n);
    }
    return y;
}

typedef struct {
    int proc_rank;
    int proc_count;

    size_t pidx;
    size_t data_per_proc;
    size_t local_size;
    size_t global_size;
} proc_data;

double* par_mul_mv(double* y_pt, const double* mat_pt, const double* x_pt, const proc_data* pdata) {
    double* res = malloc(pdata->data_per_proc * pdata->proc_count * sizeof *res);

    mul_mv_nm(res, mat_pt, x_pt, local_size, global_size);

    if (pdata->proc_count != 1) {
        double* buf_vec = malloc(pdata->data_per_proc * pdata->proc_count * sizeof *buf_vec);
        const int mul_tag = 134;

        const size_t height = ulog2(pdata->proc_count);
        for (size_t i = 0; i < height; ++i) {
            const size_t offset = 1 << i;
            const size_t group_count = offset;

            const size_t group_id = pdata->proc_rank % group_count;
            const size_t proc_per_group = pdata->proc_count / group_count;



            const size_t ingroup_index = index / group_count;
            if (ingroup_index % 2 == 0) {
                MPI_Sendrecv(res_vec, vec_size, MPI_DOUBLE, (pdata->proc_rank + offset) % pdata->proc_count, mul_tag, 
                            buf_vec, vec_size, MPI_DOUBLE, (pdata->proc_rank + offset) % pdata->proc_count, mul_tag,
                            MPI_COMM_WORLD, NULL);
            } else {
                MPI_Sendrecv(res_vec, vec_size, MPI_DOUBLE, (pdata->proc_rank - offset) % pdata->proc_count, mul_tag, 
                            buf_vec, vec_size, MPI_DOUBLE, (pdata->proc_rank - offset) % pdata->proc_count, mul_tag,
                            MPI_COMM_WORLD, NULL);
            }

            for (size_t i = 0; i < vec_size; ++i) {
                res[i] += buf_vec[i];
            }
        }
        free(buf_vec);
    }

    memcpy(y_pt, &res[pdata->pidx], pdata->local_size);
    free(res);
    return y_pt;
}*/

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

void free_pdata(proc_data* this) {
    free(this->b);
    free(this->x);
    free(this->mat);
}

size_t get_pidx(const proc_meta* this, int rank) {
    return this->data_per_proc * rank;
}

size_t get_spidx(const proc_meta* this) {
    return get_pidx(this, this->proc_rank);
}

size_t get_lsize(const proc_meta* this, int rank) {
    return (rank < this->proc_count - 1) 
            ? this->data_per_proc
            : this->global_size - this->data_per_proc * (this->proc_count - 1);
}

void solve_linear(proc_data* pdata, const proc_meta* pmeta, double e);

void* cut(const void* ptr, size_t idx, size_t size, size_t byte_size) {
    return memcpy(
        malloc(size * byte_size),
        &ptr[idx],
        size * byte_size
    );
}

void make_cuts(const proc_data* global_data, 
                 proc_data* local_data,
                 const proc_meta* pmeta) {

    local_data->mat = malloc((pmeta->local_size * pmeta->global_size) * 
                             sizeof *local_data->mat);
    for (size_t i = 0; i < global_size; ++i) {
        memcpy(
            &local_data->mat[i * pmeta->local_size], 
            &global_data->mat[i * pmeta->global_size + get_spidx(pmeta)],
            pmeta->local_size * sizeof *local_data->mat
        );
    }
    local_data->b = cut(global_data->b, get_spidx(pmeta), 
                        pmeta->local_size, sizeof *local_data->b);
    local_data->x = cut(global_data->x, get_spidx(pmeta), 
                        pmeta->local_size, sizeof *local_data->x);
}

void print_vec(const double* vec, size_t n) {
    printf("[ ");
    for (size_t i = 0; i < n; ++i) {
        printf("%f%s ", vec[i], (i == n - 1) ? "" : ",");
    }
    printf("]\n");
}

void init_pmeta(proc_meta* this, size_t global_size) {
    MPI_Comm_size(MPI_COMM_WORLD, &this->proc_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &this->proc_rank);

    this->data_per_proc = global_size / this->proc_count + (global_size % this->proc_count ? 1 : 0) ;
    this->global_size = global_size;
    this->local_size = get_lsize(this, this->proc_rank);
}

int main(int argc, char* argv[]) {
    const size_t global_size/*= init global_size*/;
    proc_data global_data;
    
    global_data.b = malloc(global_size * sizeof *global_data.b);
    //init b
    global_data.x = malloc(global_size * sizeof *global_data.x);
    //init x
    global_data.mat = malloc((global_size * global_size) *
                          sizeof *global_data.mat);
    //init mat
    const double epsilon/*= init epsilon*/;

    MPI_Init(&argc, &argv);
    proc_meta pmeta;
    init_pmeta(&pmeta, global_size);

    proc_data local_data;
    make_cuts(&global_data, &local_data, &pmeta);
    free_pdata(&global_data);

    solve_linear(&local_data, &pmeta, epsilon);

    double* final_x = NULL;
    if (pmeta.proc_rank == 0) {
        final_x = malloc(pmeta.data_per_proc * pmeta.proc_count * sizeof *final_x);
    }

    if (get_lsize(pmeta) < pmeta.data_per_proc) {
        local_data.x = realloc(local_data.x, pmeta.data_per_proc * sizeof *local_data.x);
    }

    MPI_Gather(local_data.x, pmeta.data_per_proc, MPI_DOUBLE, final_x, pmeta.data_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (pmeta.proc_rank == 0) {
        print_vec(final_x, pmeta.global_size);
    }

    free_pdata(&local_data);
    MPI_Finalize();
    return 0;
}

double par_length(const double* y_pt, const proc_meta* pmeta) {
    const double dp_pt = dot_product(y_pt, y_pt, pmeta->local_size);
    if (pmeta->proc_count == 1) return sqrt(dp_pt);

    double dp;
    MPI_Allreduce(&dp_pt, &dp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sqrt(dp);
}

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

//dot_product(ay, y, local_size) / dot_product(ay, ay, local_size)
double calc_tau(const double* ay, const double* y, const proc_meta* pmeta) {
    double nom_denom_pt[2];
    nom_denom_pt[0] = dot_product(ay, y, pmeta->local_size);
    nom_denom_pt[1] = dot_product(ay, ay, pmeta->local_size);

    if (pmeta->proc_count == 1) return nom_denom_pt[0] / nom_denom_pt[1];

    double nom_denom[2];
    MPI_Allreduce(nom_denom_pt, nom_denom, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return nom_denom[0] / nom_denom[1];
}

//Ax_(n+1) = Ax_n - t_n * Ay_n
//mat_pt = <local_size x global_size>;
void solve_linear(proc_data* pdata, const proc_meta* pmeta, double e) {
    const double b_len = par_length(pdata->b, pmeta);
    e *= b_len;

    double* ax_pt = par_mul_mv(malloc(pmeta->local_size * sizeof *ax_pt), 
                               pdata->mat, pdata->x, pmeta);

    double* y_pt = sub(ax_pt, pdata->b, pmeta->local_size);

    double* ay_pt = malloc(pmeta->local_size * sizeof *ay_pt);

    while (par_length(y_pt, pmeta) >= e) {
        par_mul_mv(ay_pt, pdata->mat, y_pt, pmeta); //ay = mat * y
        const double t = calc_tau(ay_pt, y_pt, pmeta);

        saxpy(x_pt, y_pt, -t, pmeta->local_size);
        saxpy(ax_pt, ay_pt, -t, pmeta->local_size);
        sub_s(
            memcpy(y_pt, ax_pt, pmeta->local_size * sizeof *y_pt), 
            pdata->b, 
            pmeta->local_size
        ); //y = ax - b
    }

    free(ax_pt);
    free(y_pt);
    free(ay_pt);
}
