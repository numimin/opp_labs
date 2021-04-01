#ifndef LAB_1_MATRIX_H
#define LAB_1_MATRIX_H

#include <stddef.h>

double* copy_of(const double* x, size_t n);

double* mul_mv(const double* mat, const double* x, size_t n);
double* mul_mv_s(double* y, const double* mat, const double* x, size_t n);
double* mul_mv_nm(double* y, const double* mat, const double* x, size_t n, size_t m);

double dot_product(const double* a, const double* b, size_t n);
double length(const double* a, size_t n);

double* sub(const double* a, const double* b, size_t n);
double* sub_s(double* a, const double* b, size_t n);

double* saxpy(double* y, const double* x, double alpha, size_t n);

#endif