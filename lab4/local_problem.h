#ifndef LAB4_LOCAL_PROBLEM_H
#define LAB4_LOCAL_PROBLEM_H

#include "types.h"

#include "matrix_3d.h"
#include "problem_data.h"

void ps_iterate( ProblemSolver* solver, LocalProblem* problem);
void ps_swap_cubes(LocalProblem* this);
void ps_finish_plane( ProblemSolver* solver, LocalProblem* this, int plane_index);
void ps_multiply_edges(Matrix3D* new, double factor);
void init_solver(ProblemSolver* this,  LocalProblemData* local_data,  ProblemData* data);

Matrix3D* get_old(LocalProblem* this);
Matrix3D* get_new(LocalProblem* this);

 Matrix3D* get_old_c( LocalProblem* this);
 Matrix3D* get_new_c( LocalProblem* this);

Plane* get_in_plane(LocalProblem* this, int index);
 Plane* get_in_plane_c( LocalProblem* this, int index);
double* get_out_plane(LocalProblem* this, int index);
 double* get_out_plane_c( LocalProblem* this, int index);

void init_local_problem(LocalProblem* this,  LocalProblemData* local_task);

#endif // !LAB4_LOCAL_PROBLEM_H