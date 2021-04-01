#include "matrix.h"

#include <malloc.h>
#include <string.h>
#include <math.h>

#include <omp.h>

double* par_memcpy(double* dest, const double* src, size_t n) {
    const size_t num_threads = omp_get_num_threads();
    const size_t thread_num = omp_get_thread_num();

    const size_t data_per_thread = n / num_threads + (n % num_threads ? 1 : 0);
    const size_t chunk_size = (thread_num < num_threads - 1) 
                            ? data_per_thread
                            : n - data_per_thread * (num_threads - 1);

    const size_t index = data_per_thread * thread_num;
    memcpy(&dest[index], &src[index], chunk_size * sizeof(double));
    #pragma omp barrier
    return dest;
}

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
    #pragma omp for schedule(runtime)
    for (size_t i = 0; i < n; ++i) {
        y[i] += alpha * x[i];
    }
    return y;
}

double* sub_s(double* a, const double* b, size_t n) {
    #pragma omp for schedule(runtime)
    for (size_t i = 0; i < n; ++i) {
        a[i] -= b[i];
    }
    return a;
}

double dot_product(const double* a, const double* b, size_t n) {
    static double res;
    
    #pragma omp barrier
    res = 0.0;
    #pragma omp barrier
    #pragma omp for reduction(+:res) schedule(runtime)
    for (size_t i = 0; i < n; ++i) {
        res += a[i] * b[i];
    }
    return res;
}

double np_dot_product(const double* a, const double* b, size_t n) {
    double res = 0.0;
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
    #pragma omp for schedule(runtime)
    for (size_t i = 0; i < m; ++i) {
        y[i] = np_dot_product(&mat[n * i], x, n);
    }
    return y;
}
