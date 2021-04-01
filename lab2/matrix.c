#include "matrix.h"

#include <malloc.h>
#include <string.h>
#include <math.h>

double length(const double* a, size_t n) {
    return sqrt(dot_product(a, a, n));
}

double* copy_of(const double* x, size_t n) {
    double* res = malloc(n * sizeof *res);
    return memcpy(res, x, n * sizeof *res);
}

double* sub(const double* a, const double* b, size_t n) {
    return sub_s(copy_of(a, n), b, n);
}

double* saxpy(double* y, const double* x, double alpha, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        y[i] += alpha * x[i];
    }
    return y;
}

double* sub_s(double* a, const double* b, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        a[i] -= b[i];
    }
    return a;
}

double dot_product(const double* a, const double* b, size_t n) {
    double res = 0;
    for (size_t i = 0; i < n; ++i) {
        res += a[i] * b[i];
    }
    return res;
}

double* mul_mv(const double* mat, const double* x, size_t n) {
    return mul_mv_s(malloc(n * sizeof *x), mat, x, n);
}

double* mul_mv_s(double* y, const double* mat, const double* x, size_t n) {
    return mul_mv_nm(y, mat, x, n, n);
}

double* mul_mv_nm(double* y, const double* mat, const double* x, size_t n, size_t m) {
    for (size_t i = 0; i < m; ++i) {
        y[i] = dot_product(&mat[n * i], x, n);
    }
    return y;
}
