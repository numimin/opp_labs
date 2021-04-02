#ifndef LAB4_TYPES_H
#define LAB4_TYPES_H

#include <stddef.h>
#include <stdbool.h>
#include <mpi.h>

#include "functional.h"

#define NB_LEFT 0 //(X, F)
#define NB_RIGHT 1 //(X, T)
#define NB_FRONT 2 //(Y, F)
#define NB_BACK 3 //(X, T)
#define NB_UP 4 //(Z, F)
#define NB_DOWN 5 //(Z, T)

#define NB_COUNT 6

#define DIMS 3
#define X 0 
#define Y 1 
#define Z 2 

typedef struct {
    int proc_count;
    int proc_rank;

    MPI_Comm global_comm;

    int comm_dims[DIMS]; //(x, y, z)
    int comm_coords[DIMS];

    MPI_Datatype plane_types[DIMS];
    int neighbours[NB_COUNT];

    size_t task_dims[DIMS];
    size_t task_coords[DIMS];
} ProcMeta;

typedef struct {
    size_t dimensions[3];
    double* data;
} Matrix3D;

typedef struct {
    size_t dimensions[2];
    double* data;
} Plane;

typedef struct {
    double factor;
    double h_2[DIMS];    
    Matrix3D rho;
} Constants;

typedef struct {
    Matrix3D slice[2]; //old and new
    Plane neighbours[NB_COUNT];

    size_t old;
    size_t new;

    Constants constants;

    PlaneIterator plane_mul_iterators[NB_COUNT];
    PlaneIterator plane_iterators[NB_COUNT];
    PlaneIterator plane_edged_iterators[NB_COUNT];
    size_t plane_neighbor_indices[NB_COUNT];
} LocalProblem;

typedef double (*func_r3) (double x, double y, double z);

typedef struct {
    double area_offset[DIMS];
    double global_area[DIMS];

    double a;
    func_r3 phi;
    func_r3 rho;

    size_t discrete_dimensions[DIMS];
    double local_area[DIMS];
} ProblemData;

#endif // !LAB4_TYPES_H