#ifndef LAB4_LOCAL_PROBLEM_H
#define LAB4_LOCAL_PROBLEM_H

#include "types.h"

#include "matrix_3d.h"
#include "problem_data.h"

void lp_init(LocalProblem* this, const ProblemData* data, const ProcMeta* meta);

void lp_iterate(LocalProblem* this);
void lp_swap_cubes(LocalProblem* this);
void lp_finish_plane(LocalProblem* this, size_t plane_index);
void lp_multiply_edges(LocalProblem* this);

Matrix3D* lp_old_mat(LocalProblem* this);
Matrix3D* lp_new_mat(LocalProblem* this);

Plane* lp_in_plane(LocalProblem* this, size_t index);
double* lp_out_plane(LocalProblem* this, size_t index);

#endif // !LAB4_LOCAL_PROBLEM_H