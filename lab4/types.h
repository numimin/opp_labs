#ifndef LAB4_TYPES_H
#define LAB4_TYPES_H

#include <stddef.h>
#include <stdbool.h>
#include <mpi.h>

#define NB_LEFT 0 
#define NB_RIGHT 1 

#define NB_COUNT 2

#define DIMS 3
#define X 0 
#define Y 1 
#define Z 2 

typedef struct {
    int proc_count;
    int proc_rank;
    int proc_coord;

    int data_per_proc;
    int last_data_count;

    MPI_Comm global_comm;

    int neighbours[NB_COUNT];

    int task_dims[DIMS];
    int task_coords[DIMS];
} ProcMeta;

typedef struct {
    int dimensions[3];
    double* data;
} Matrix3D;

typedef struct {
    double factor;
    double h_2[DIMS];    
    Matrix3D rho;
} Constants;

typedef struct {
    Matrix3D slice[2]; //old and new
    int slice_dims[DIMS];
    int inner_dims[DIMS];
    int inner_offset[DIMS];

    int plane_x_coords[NB_COUNT];
    int plane_neighbour_x_coords[NB_COUNT];

    int old;
    int new;

    Constants constants;
} LocalProblem;

typedef double (*func_r3) (double x, double y, double z);

typedef struct {
    double area_offset[DIMS];
    double global_area[DIMS];

    double a;
    func_r3 phi;
    func_r3 rho;

    int discrete_dimensions[DIMS];
    double local_area[DIMS];
} ProblemData;

#endif // !LAB4_TYPES_H