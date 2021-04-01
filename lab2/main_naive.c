#include <malloc.h>
#include <string.h>

#include <omp.h>

#include "data.h"
#include "matrix.h"

void solve_linear(const double* mat, const double* b, double* x, size_t n, double e);

int main(int argc, char* argv[]) {
    const double epsilon = E;
    const size_t n = N;

    proc_data pdata;
    fill_pdata(&pdata, n);

    const double start = omp_get_wtime();
    solve_linear(pdata.mat, pdata.b, pdata.x, n, epsilon);
    free_pdata(&pdata);
    const double end = omp_get_wtime();

    printf("Time: %f\n", end - start);
    return 0;
}

//Ax_(n+1) = Ax_n - t_n * Ay_n
void solve_linear(const double* mat, const double* b, double* x, size_t n, double e) {
    const double b_len = length(b, n);
    e *= b_len;

    double* ax = mul_mv(mat, x, n);
    double* y = sub(ax, b, n);
    double* ay = malloc(n * sizeof *ay);

    int i = 0;
    while (length(y, n) >= e) {
        ++i;
        mul_mv_s(ay, mat, y, n); //ay = mat * y

        const double t = dot_product(ay, y, n) / dot_product(ay, ay, n);
        saxpy(x, y, -t, n);

        saxpy(ax, ay, -t, n);

        sub_s(memcpy(y, ax, n * sizeof *y), b, n); //y = ax - b
    }
    printf("Iterations: %d\n", i);

    free(ax);
    free(y);
    free(ay);
}
