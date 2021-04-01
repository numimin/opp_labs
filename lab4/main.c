#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include "problem_data.h"
#include "matrix_3d.h"
#include "local_problem.h"
#include "problem_data.h"
#include "proc_meta.h"

typedef struct {
    ProcMeta meta;
    ProblemData task;
    LocalProblemData local_task;
    LocalProblem problem;
} ProcData;

int init_procdata(ProcData* this, int dimensions[DIMS]) {
    pd_init(&this->task);
    if (pm_init(&this->meta, &this->local_task, &this->task, dimensions) != EXIT_SUCCESS) return EXIT_FAILURE;
    init_local_problem(&this->problem, &this->local_task);

    return EXIT_SUCCESS;
}

bool local_solved( Matrix3D* old,  Matrix3D* new, double epsilon) {
    for (int x = 0; x < mt_x_len(old); ++x) {
        for (int y = 0; y < mt_y_len(old); ++y) {
            for (int z = 0; z < mt_z_len(old); ++z) {
                if (fabs(mt_get(new, x, y, z) - mt_get(old, x, y, z)) 
                    >= epsilon) return false;
            }
        }
    }

    return true;
}

bool global_solved( ProcMeta* meta,  Matrix3D* old,  Matrix3D* new, double epsilon) {
     bool solved_here = local_solved(old, new, epsilon);
    bool* buf = NULL;
    if (pm_is_root(meta)) {
        buf = malloc(pm_count(meta) * sizeof *buf);
    } 

    MPI_Gather(&solved_here, 1, MPI_C_BOOL, 
               buf, pm_count(meta), MPI_C_BOOL, 
               pm_root(meta), pm_comm(meta));

    bool solved_all = true;
    if (pm_is_root(meta)) {
        for (int i = 0; i < pm_count(meta); ++i) {
            if (!buf[i]) {
                solved_all = false;
                break;
            }
        }
        free(buf);
    }

    MPI_Bcast(&solved_all, 1, MPI_C_BOOL, pm_root(meta), pm_comm(meta));
    return solved_all;
}

void isend_plane( ProcData* data, int plane_index, MPI_Request out_requests[NB_COUNT], int tag) {
     double* plane = get_out_plane_c(&data->problem, plane_index);
    MPI_Isend(plane, 1, pm_plane(&data->meta, plane_index),
        pm_neighbour(&data->meta, plane_index), tag, pm_comm(&data->meta), 
        &out_requests[plane_index]);
}

void irecv_plane(ProcData* data, int plane_index, MPI_Request in_requests[NB_COUNT], int tag) {
     Plane* plane = get_in_plane(&data->problem, plane_index);
    MPI_Irecv(plane->data, pl_x_len(plane) * pl_y_len(plane), MPI_DOUBLE,
        pm_neighbour(&data->meta, plane_index), tag, pm_comm(&data->meta), 
        &in_requests[plane_index]);
}

void isend_neighbour_planes( ProcData* data, MPI_Request out_requests[NB_COUNT], int tag) {
    for (int i = 0; i < NB_COUNT; ++i) {
        isend_plane(data, i, out_requests, tag);
    }
}

void solve_equation(ProcData* data, double epsilon) {
    const int tag = 123;

    MPI_Request out_requests[NB_COUNT];
    isend_neighbour_planes(data, out_requests, tag);

    ProblemSolver solver;
    init_solver(&solver, &data->local_task, &data->task);

    while (!global_solved(
        &data->meta,
        get_old(&data->problem),
        get_new(&data->problem), epsilon)) {

        MPI_Request in_requests[NB_COUNT];
        for (int i = 0; i < NB_COUNT; ++i) {
            irecv_plane(data, i, in_requests, tag);
        }

        ps_swap_cubes(&data->problem);
        ps_iterate(&solver, &data->problem);

        int planes_processed = 0;
        while (planes_processed < NB_COUNT) {
            int plane_index;
            MPI_Waitany(NB_COUNT, in_requests, &plane_index, MPI_STATUS_IGNORE);
            in_requests[plane_index] = MPI_REQUEST_NULL;

            ps_finish_plane(&solver, &data->problem, plane_index);
        }
        ps_multiply_edges(get_new(&data->problem), solver.factor);

        MPI_Waitall(NB_COUNT, out_requests, MPI_STATUSES_IGNORE);
        isend_neighbour_planes(data, out_requests, tag);
    }
}

void test_solution(size_t x, size_t y, size_t z) {
    ProcData data;
    int dimensions[DIMS] = {x, y, z};
    init_procdata(&data, dimensions);
    
    const double start = MPI_Wtime();
    solve_equation(&data, EPSILON);
    const double end = MPI_Wtime();

    if (pm_is_root(&data.meta)) {
        printf("Cores %d [%d x %d x %d]: %f\n", 
            pm_count(&data.meta), x, y, z, end - start);
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int proc_count;
    MPI_Comm_size(MPI_COMM_WORLD, &proc_count);

    for (size_t x = 1; x <= proc_count; ++x) {
        if (proc_count % x) continue;

        for (size_t y = 1; y <= (proc_count / x); ++y) {
            if ((proc_count / x) % y) continue;

            const size_t z = (proc_count / x) / y;
            test_solution(x, y, z);
        }
    }


    MPI_Finalize();
    return EXIT_SUCCESS;
}
