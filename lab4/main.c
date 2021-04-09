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
    LocalProblem problem;
} ProcData;

int init_procdata(ProcData* this) {
    pd_init(&this->task);
    pm_init(&this->meta, &this->task);
    lp_init(&this->problem, &this->task, &this->meta);

    return EXIT_SUCCESS;
}

void free_procdata(ProcData* this) {
    pm_free(&this->meta);
    lp_free(&this->problem);
}

bool local_solved(const Matrix3D* old, const Matrix3D* new, double epsilon) {
    double delta;

    for (int x = 1; x < mt_x_len(old) - 1; ++x) {
        for (int y = 1; y < mt_y_len(old) - 1; ++y) {
            for (int z = 1; z < mt_z_len(old) - 1; ++z) {
                delta = fabs(mt_get(new, x, y, z) - mt_get(old, x, y, z));
                if (delta >= epsilon) return false;
            }
        }
    }

    return true;
}

bool global_solved(const ProcMeta* meta, const Matrix3D* old, const Matrix3D* new, double epsilon) {
    const bool solved_here = local_solved(old, new, epsilon);
    bool* buf = NULL;
    if (pm_is_root(meta)) {
        buf = malloc(pm_count(meta) * sizeof *buf);
    } 

    MPI_Gather(&solved_here, 1, MPI_C_BOOL, 
               buf, 1, MPI_C_BOOL, 
               pm_root(meta), pm_comm(meta));

    bool solved_all;
    if (pm_is_root(meta)) {
        solved_all = true;
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

void isend_plane(const ProcData* data, int plane_index, MPI_Request out_requests[NB_COUNT], int tag) {
    const int neighbour = pm_neighbour(&data->meta, plane_index);
    if (neighbour == NO_NEIGHBOUR) {
        out_requests[plane_index] = MPI_REQUEST_NULL;
        return;
    }

    const double* plane = lp_out_plane_c(&data->problem, plane_index);
    
    MPI_Isend(plane, pm_disc_len_y(&data->meta) * pm_disc_len_z(&data->meta), MPI_DOUBLE, 
        neighbour, tag, pm_comm(&data->meta), 
        &out_requests[plane_index]);
}

void irecv_plane(ProcData* data, int plane_index, MPI_Request in_requests[NB_COUNT], int tag) {
    const int neighbour = pm_neighbour(&data->meta, plane_index);
    if (neighbour == NO_NEIGHBOUR) {
        in_requests[plane_index] = MPI_REQUEST_NULL;
        return;
    }

    double* plane = lp_in_plane(&data->problem, plane_index);

    MPI_Irecv(plane, pm_disc_len_y(&data->meta) * pm_disc_len_z(&data->meta), MPI_DOUBLE,
        neighbour, tag, pm_comm(&data->meta), 
        &in_requests[plane_index]);
}

void isend_neighbour_planes(const ProcData* data, MPI_Request out_requests[NB_COUNT], int tag) {
    for (int i = 0; i < NB_COUNT; ++i) {
        isend_plane(data, i, out_requests, tag);
    }
}

void irecv_neighbour_planes(ProcData* data, MPI_Request in_requests[NB_COUNT], int tag) {
    for (int i = 0; i < NB_COUNT; ++i) {
        irecv_plane(data, i, in_requests, tag);
    }
}

void solve_equation(ProcData* data, double epsilon) {
    const int tag = 123;
    MPI_Request out_requests[NB_COUNT];
    for (int i = 0; i < NB_COUNT; ++i) {
        out_requests[i] = MPI_REQUEST_NULL;
    }
    MPI_Request in_requests[NB_COUNT];

    int iterations = 0;
    do {
        ++iterations;

        MPI_Waitall(NB_COUNT, out_requests, MPI_STATUSES_IGNORE);
        isend_neighbour_planes(data, out_requests, tag);

        lp_swap_cubes(&data->problem);

        irecv_neighbour_planes(data, in_requests, tag);
        lp_iterate(&data->problem);

        for (int i = 0; i < NB_COUNT; ++i) {
            if (in_requests[i] != MPI_REQUEST_NULL) continue;
            lp_finish_plane(&data->problem, i);
        }

        for (;;) {
            int plane_index;
            MPI_Waitany(NB_COUNT, in_requests, &plane_index, MPI_STATUS_IGNORE);
            if (plane_index == MPI_UNDEFINED) break;

            in_requests[plane_index] = MPI_REQUEST_NULL;

            lp_finish_plane(&data->problem, plane_index);
        }
        
    } while (!global_solved(
        &data->meta,
        lp_old_mat(&data->problem),
        lp_new_mat(&data->problem), epsilon)
    );

    if (pm_is_root(&data->meta)) {
        fprintf(stderr, "Iterations: %d\n", iterations);
    }
}

double get_precision(const double* solution, const ProblemData* task) {
    double max_delta = 0.0;

    const int len[DIMS] = {pd_disc_x(task), pd_disc_y(task), pd_disc_z(task)};

    for (int x = 0; x < len[X]; ++x) {
        for (int y = 0; y < len[Y]; ++y) {
            for (int z = 0; z < len[Z]; ++z) {
                const double phi = pd_phi(task, x, y, z);
                const double solved = solution[len[Y] * len[Z] * x + len[Z] * y + z];

                const double delta = fabs(solved - phi) / phi;
                if (delta > max_delta) {
                    max_delta = delta;
                }
            }
        }
    }

    return max_delta;
}

double gather_solution(const ProcData* data) {
    double* solution = NULL;
    int* recvcounts = NULL;
    int* displs = NULL;

    const Matrix3D* local_solution = lp_new_mat_c(&data->problem);
    const int plane_size = mt_y_len(local_solution) * mt_z_len(local_solution);
    
    if (pm_is_root(&data->meta)) {
        solution = malloc(pd_disc_x(&data->task) * 
            pd_disc_y(&data->task) * pd_disc_z(&data->task) * sizeof *solution);
        recvcounts = malloc(pm_count(&data->meta) * sizeof *recvcounts);
        displs = malloc(pm_count(&data->meta) * sizeof *displs);

        for (int i = 0; i < pm_count(&data->meta) - 1; ++i) {
            recvcounts[i] = data->meta.data_per_proc;
        }
        recvcounts[pm_count(&data->meta) - 1] = data->meta.last_data_count;

        ++recvcounts[0];
        ++recvcounts[pm_count(&data->meta) - 1];

        for (int i = 0; i < pm_count(&data->meta); ++i) {
            recvcounts[i] *= plane_size;
        }

        int curr_displ = 0;
        for (int i = 0; i < pm_count(&data->meta); ++i) {
            displs[i] = curr_displ;
            curr_displ += recvcounts[i];
        }
    }

    const bool left = pm_neighbour(&data->meta, NB_LEFT) == NO_NEIGHBOUR;
    const bool right = pm_neighbour(&data->meta, NB_RIGHT) == NO_NEIGHBOUR;

    int plane_count = pm_disc_len_x(&data->meta);
    if (left) ++plane_count;
    if (right) ++plane_count;

    MPI_Gatherv(mt_get_ref(local_solution, (left) ? 0 : 1, 0, 0), 
        plane_count * plane_size, MPI_DOUBLE, 
        solution, recvcounts, displs,
        MPI_DOUBLE, pm_root(&data->meta), pm_comm(&data->meta));

    if (pm_is_root(&data->meta)) {
        const double precision = get_precision(solution, &data->task);

        free(solution);
        free(displs);
        free(recvcounts);

        return precision;
    }
    return 0;
}

void test_solution() {
    ProcData data;
    init_procdata(&data);

    const double start = MPI_Wtime();
    solve_equation(&data, EPSILON);
    const double precision = gather_solution(&data);
    const double end = MPI_Wtime();

    if (pm_is_root(&data.meta)) {
        printf("Cores %d: %f; Precision: %.10f\n", 
            pm_count(&data.meta), end - start, precision);
    }

    free_procdata(&data);
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    test_solution();
    MPI_Finalize();

    return EXIT_SUCCESS;
}
