#ifndef  LAB1_MAIN_H
#define  LAB1_MAIN_H

#include "matrix.h"
#include "data.h"

void init_pmeta(proc_meta* this, size_t global_size);
void fill_pdata(proc_data* pdata, const proc_meta* pmeta);

double par_length(const double* y_pt, const proc_meta* pmeta);
double calc_tau(const double* ay, const double* y, const proc_meta* pmeta);
double* par_mul_mv(double* y, const double* mat_pt, const double* x_pt, const proc_meta* pmeta);
void par_solve_linear(proc_data* pdata, const proc_meta* pmeta, double e);

#endif // ! LAB1_MAIN_H