#include <malloc.h>
#include <string.h>

#include <unistd.h>
#include <sys/times.h>

#include "data.h"
#include "matrix.h"

void solve_linear(const double* mat, const double* b, double* x, size_t n, double e);

void fill_pdata(proc_data* pdata, const proc_meta* pmeta) {
    const size_t n = pmeta->local_size;
    const size_t m = pmeta->global_size;

    pdata->b = malloc(n * sizeof *pdata->b);
    pdata->x = malloc(n * sizeof *pdata->x);
    for (size_t i = 0; i < n; ++i) {
        pdata->x[i] = sin(2 * M_PI * (i + get_spidx(pmeta)) / m);
    }

    pdata->mat = malloc(n * m * sizeof *pdata->mat);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            pdata->mat[n * i + j] = fill_mat(i, j + get_spidx(pmeta));
        }
    }

    mul_mv_s(pdata->b, pdata->mat, pdata->x, n);

    for (size_t i = 0; i < n; ++i) {
        pdata->x[i] = 0;
    }
}

int main(int argc, char* argv[]) {
    const double epsilon = E;
    proc_meta pmeta = {.data_per_proc=N, .global_size=N, .local_size=N, .proc_count=1, .proc_rank=0};
    proc_data pdata;
    fill_pdata(&pdata, &pmeta);

    struct tms start, end;
    const long clocks_per_sec = sysconf(_SC_CLK_TCK);

    times(&start);
    solve_linear(pdata.mat, pdata.b, pdata.x, pmeta.local_size, epsilon);
    free_pdata(&pdata);
    times(&end);

    const long clocks = end.tms_utime - start.tms_utime;
    printf("Cores: 1: %f\n", (double) clocks / clocks_per_sec);
    return 0;
}

//Ax_(n+1) = Ax_n - t_n * Ay_n
void solve_linear(const double* mat, const double* b, double* x, size_t n, double e) {
    double* ax = mul_mv(mat, x, n);
    double* y = sub(ax, b, n);
    double* ay = malloc(n * sizeof *ay);

    int i = 0;
    while ((length(y, n) / length(b, n)) >= e) {
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
