#include "proc_meta.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int get_data_per_proc(int global_size, int proc_count) {
    return global_size / proc_count + (global_size % proc_count ? 1 : 0);
}

int get_last_cut_size(int global_size, int proc_count) {
     int data_per_proc = get_data_per_proc(global_size, proc_count);
    return global_size - data_per_proc * (proc_count - 1);
}

int get_cut_size(int global_size, int proc_count, int rank) {
    return (rank < proc_count - 1) 
            ? get_data_per_proc(global_size, proc_count)
            : get_last_cut_size(global_size, proc_count);
}

int get_cart_rank(const ProcMeta* this, int coord) {
    int rank;
    MPI_Cart_rank(pm_comm(this), &coord, &rank);
    return rank;
}

void determine_neighbours(ProcMeta* this) {
    const int coords[NB_COUNT] = {pm_coord(this) - 1, pm_coord(this) + 1};

    for (int i = 0; i < NB_COUNT; ++i) {
        if (coords[i] < 0 || coords[i] >= pm_count(this)) {
            this->neighbours[i] = NO_NEIGHBOUR;
        } else {
            this->neighbours[i] = get_cart_rank(this, coords[i]);
        }
    }
}

void cube_without_edges(const ProblemData* data, int dimensions[DIMS]) {
    for (int coord = 0; coord < DIMS; ++coord) {
        dimensions[coord] = pd_disc(data, coord) - 2;
    }
}

void calculate_cuts(ProcMeta* this, const int dimensions[DIMS]) {
    memcpy(this->task_dims, dimensions, DIMS * sizeof *this->task_dims);
    this->task_dims[X] = get_cut_size(dimensions[X], pm_count(this), pm_coord(this));

    memset(this->task_coords, 0, DIMS * sizeof *this->task_dims);
    this->task_coords[X] = pm_coord(this) * get_data_per_proc(dimensions[X], pm_count(this));
}

void pm_init(ProcMeta* this, const ProblemData* data) {
    MPI_Comm_size(MPI_COMM_WORLD, &this->proc_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &this->proc_rank);

    const int periodic = 0;
    MPI_Cart_create(MPI_COMM_WORLD, 1, &this->proc_count, &periodic, 1, &this->global_comm);
    MPI_Cart_coords(this->global_comm, this->proc_rank, 1, &this->proc_coord);
    MPI_Cart_rank(this->global_comm, &this->proc_coord, &this->proc_rank);

    determine_neighbours(this);

    int dims_without_edges[DIMS];
    cube_without_edges(data, dims_without_edges);
    calculate_cuts(this, dims_without_edges);

    this->data_per_proc = get_data_per_proc(dims_without_edges[X], pm_count(this));
    this->last_data_count = get_last_cut_size(dims_without_edges[X], pm_count(this));
}

void pm_free(ProcMeta* this) {
    MPI_Comm_free(&this->global_comm);
}

int pm_neighbour(const ProcMeta* this, int index) {
    return this->neighbours[index];
}

MPI_Comm pm_comm(const ProcMeta* this) {
    return this->global_comm;
}

int pm_coord(const ProcMeta* this) {
    return this->proc_coord;
}

int pm_root(const ProcMeta* this) {
    return get_cart_rank(this, 0);
}

int pm_rank(const ProcMeta* this) {
    return this->proc_rank;
}

bool pm_is_root(const ProcMeta* this) {
    return (pm_rank(this) == pm_root(this));
}

int pm_count(const ProcMeta* this) {
    return this->proc_count;
}

int pm_disc(const ProcMeta* this, int coord, int local_coord) {
    return this->task_coords[coord] + 1 + local_coord;
}

int pm_disc_x(const ProcMeta* this, int local_x) {
    return pm_disc(this, X, local_x);
}

int pm_disc_y(const ProcMeta* this, int local_y) {
    return pm_disc(this, Y, local_y);
}

int pm_disc_z(const ProcMeta* this, int local_z) {
    return pm_disc(this, Z, local_z);
}

int pm_disc_len(const ProcMeta* this, int coord) {
    return this->task_dims[coord];
}

int pm_disc_len_x(const ProcMeta* this) {
    return pm_disc_len(this, X);
}

int pm_disc_len_y(const ProcMeta* this) {
    return pm_disc_len(this, Y);
}

int pm_disc_len_z(const ProcMeta* this) {
    return pm_disc_len(this, Z);
}
