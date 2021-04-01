#include <mpi.h>
#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "main.h"

void init_pmeta(proc_meta* this, size_t global_size) {
    MPI_Comm_size(MPI_COMM_WORLD, &this->proc_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &this->proc_rank);

    this->data_per_proc = global_size / this->proc_count + (global_size % this->proc_count ? 1 : 0) ;
    this->global_size = global_size;
    this->local_size = get_lsize(this, this->proc_rank);
}

double par_length(const double* y_pt, const proc_meta* pmeta) {
    const double dp_pt = dot_product(y_pt, y_pt, pmeta->local_size);

    double dp;
    MPI_Allreduce(&dp_pt, &dp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sqrt(dp);
}

//dot_product(ay, y, local_size) / dot_product(ay, ay, local_size)
double calc_tau(const double* ay, const double* y, const proc_meta* pmeta) {
    double nom_denom_pt[2];
    nom_denom_pt[0] = dot_product(ay, y, pmeta->local_size);
    nom_denom_pt[1] = dot_product(ay, ay, pmeta->local_size);

    double nom_denom[2];
    MPI_Allreduce(nom_denom_pt, nom_denom, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return nom_denom[0] / nom_denom[1];
}

void par_solve_linear(proc_data* pdata, const proc_meta* pmeta, double e) {
    const double b_len = par_length(pdata->b, pmeta);
    e *= b_len;

    double* ax_pt = par_mul_mv(malloc(pmeta->local_size * sizeof *ax_pt), 
                               pdata->mat, pdata->x, pmeta);

    double* y_pt = sub(ax_pt, pdata->b, pmeta->local_size);

    double* ay_pt = malloc(pmeta->local_size * sizeof *ay_pt);

    while (par_length(y_pt, pmeta) >= e) {
        par_mul_mv(ay_pt, pdata->mat, y_pt, pmeta); //ay = mat * y
        const double t = calc_tau(ay_pt, y_pt, pmeta);

        saxpy(pdata->x, y_pt, -t, pmeta->local_size);
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

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    const double epsilon = E;

    proc_meta pmeta;
    init_pmeta(&pmeta, N);

    proc_data pdata;
    fill_pdata(&pdata, &pmeta);

    const double start = MPI_Wtime();
    par_solve_linear(&pdata, &pmeta, epsilon);

    double* final_x = NULL;
    if (pmeta.proc_rank == 0) {
        final_x = malloc(pmeta.data_per_proc * pmeta.proc_count * sizeof *final_x);
    }

    if (get_lsize(&pmeta, pmeta.proc_rank) < pmeta.data_per_proc) {
        pdata.x = realloc(pdata.x, pmeta.data_per_proc * sizeof *pdata.x);
    }

    MPI_Gather(pdata.x, pmeta.data_per_proc, MPI_DOUBLE, final_x, pmeta.data_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (pmeta.proc_rank == 0) {
        free(final_x);
    }
    free_pdata(&pdata);
    const double end = MPI_Wtime();
    printf("Cores %d: %f\n", pmeta.proc_count, end - start);

    MPI_Finalize();
    return 0;
}
