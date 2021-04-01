#include "data.h"

#include <math.h>
#include <malloc.h>
#include <stdlib.h>

#include "matrix.h"

size_t get_lsize(const proc_meta* this, int rank) {
    return (rank < this->proc_count - 1) 
            ? this->data_per_proc
            : this->global_size - this->data_per_proc * (this->proc_count - 1);
}

size_t get_pidx(const proc_meta* this, int rank) {
    return this->data_per_proc * rank;
}

size_t get_spidx(const proc_meta* this) {
    return get_pidx(this, this->proc_rank);
}

void free_pdata(proc_data* this) {
    free(this->b);
    free(this->x);
    free(this->mat);
}

double fill_mat(size_t i, size_t j) {
    return (i == j) ? (-4) 
                    : (abs(i - j) == 1) ? 1
                    : (abs(i - j) == NX) ? 1
                    : 0;
}
