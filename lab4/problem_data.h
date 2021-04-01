#ifndef LAB4_PROBLEM_DATA_H
#define LAB4_PROBLEM_DATA_H

#include "types.h"
#include "local_problem.h"

#define X_0 (-1)
#define Y_0 (-1)
#define Z_0 (-1)

#define D_X 2
#define D_Y 2
#define D_Z 2

#define A 100000

#define EPSILON 0.00000001

#define N_X 1000 
#define N_Y 1000 
#define N_Z 1000 

double phi(double x, double y, double z);
double rho(double x, double y, double z);

double pd_h(const ProblemData* this, size_t coord);
double pd_h_x(const ProblemData* this);
double pd_h_y(const ProblemData* this);
double pd_h_z(const ProblemData* this);

double pd_area(const ProblemData* this, size_t coord, size_t index);
double pd_area_x(const ProblemData* this, size_t x);
double pd_area_y(const ProblemData* this, size_t y);
double pd_area_z(const ProblemData* this, size_t z);

size_t pd_disc(const ProblemData* this, size_t coord);
size_t pd_disc_x(const ProblemData* this);
size_t pd_disc_y(const ProblemData* this);
size_t pd_disc_z(const ProblemData* this);

double pd_a(const ProblemData* this);
double pd_factor(const ProblemData* this);

double pd_func(const ProblemData* this, const size_t discrete_point[DIMS], func_r3 func);
double pd_rho(const ProblemData* this, const size_t discrete_point[DIMS]);
double pd_phi(const ProblemData* this, const size_t discrete_point[DIMS]);

void pd_init(ProblemData* this);

#endif // !LAB4_PROBLEM_DATA_H