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

#define N_X 300
#define N_Y 300
#define N_Z 300

double phi(double x, double y, double z);
double rho(double x, double y, double z);

double pd_h(const ProblemData* this, int coord);
double pd_h_x(const ProblemData* this);
double pd_h_y(const ProblemData* this);
double pd_h_z(const ProblemData* this);

double pd_area(const ProblemData* this, int coord, int index);
double pd_area_x(const ProblemData* this, int x);
double pd_area_y(const ProblemData* this, int y);
double pd_area_z(const ProblemData* this, int z);

int pd_disc(const ProblemData* this, int coord);
int pd_disc_x(const ProblemData* this);
int pd_disc_y(const ProblemData* this);
int pd_disc_z(const ProblemData* this);

double pd_a(const ProblemData* this);
double pd_factor(const ProblemData* this);

double pd_func(const ProblemData* this, const int discrete_point[DIMS], func_r3 func);
double pd_rho_ar(const ProblemData* this, const int discrete_point[DIMS]);
double pd_rho(const ProblemData* this, int x, int y, int z);
double pd_phi_ar(const ProblemData* this, const int discrete_point[DIMS]);
double pd_phi(const ProblemData* this, int x, int y, int z);

void pd_init(ProblemData* this);

#endif // !LAB4_PROBLEM_DATA_H