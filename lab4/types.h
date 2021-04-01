#ifndef LAB4_TYPES_H
#define LAB4_TYPES_H

#include <stddef.h>
#include <stdbool.h>
#include <mpi.h>

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
    double dimensions[DIMS];
} Vector;

typedef struct {
    int dimensions[3];
    double* data;
} Matrix3D;

typedef struct {
    int coords[2];
    int dimensions[2];
    int excluded_coord;
    bool right;
    double* data;
} Plane;

typedef struct {
    Matrix3D slice[2]; //old and new
    Plane neighbours[NB_COUNT];

    int old;
    int new;
} LocalProblem;

typedef struct {
    double factor;
    double h_2[DIMS];    
    Matrix3D rho;
} ProblemSolver;

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

/*typedef struct {
    Vector offset;
    int dimensions[DIMS];
    int coords[DIMS];
} LocalProblemData;*/

#endif // !LAB4_TYPES_H