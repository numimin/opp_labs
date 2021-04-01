#include "proc_meta.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

size_t get_data_per_proc(size_t global_size, int proc_count) {
    return global_size / proc_count + (global_size % proc_count ? 1 : 0);
}

size_t get_last_cut_size(size_t global_size, int proc_count) {
     size_t data_per_proc = get_data_per_proc(global_size, proc_count);
    return global_size - data_per_proc * (proc_count - 1);
}

size_t get_cut_size(size_t global_size, int proc_count, int rank) {
    return (rank < proc_count - 1) 
            ? get_data_per_proc(global_size, proc_count)
            : get_last_cut_size(global_size, proc_count);
}

int get_cart_rank(const ProcMeta* this, const int coords[DIMS]) {
    int rank;
    MPI_Cart_rank(pm_comm(this), coords, &rank);
    return rank;
}

void determine_neighbours(ProcMeta* this) {
    const int coords[NB_COUNT][DIMS] = {
        {pm_x(this) + 1, pm_y(this), pm_z(this)},
        {pm_x(this) - 1, pm_y(this), pm_z(this)},
        {pm_x(this), pm_y(this) + 1, pm_z(this)},
        {pm_x(this), pm_y(this) - 1, pm_z(this)},
        {pm_x(this), pm_y(this), pm_z(this) + 1},
        {pm_x(this), pm_y(this), pm_z(this) - 1},
    };

    int changed_coord = 0;
    for (int i = 0; i < NB_COUNT; ++i) {
        if (coords[i][changed_coord] < 0 || coords[i][changed_coord] > pm_len(this, changed_coord)) {
            this->neighbours[i] = NO_NEIGHBOUR;
        } else {
            this->neighbours[i] = get_cart_rank(this, coords[i]);
        }

        if (!(i % 2)) ++changed_coord;
    }
}

void cube_without_edges(const ProblemData* data, size_t offset[DIMS], size_t dimensions[DIMS]) {
    for (size_t coord = 0; coord < DIMS; ++coord) {
        dimensions[coord] = pd_disc(data, coord);
        offset[coord] = 1;
    }
}

void calculate_cuts(ProcMeta* this, const size_t offset[DIMS], const size_t dimensions[DIMS]) {
    for (size_t coord = 0; coord < DIMS; ++coord) {
        this->task_dims[coord] = get_cut_size(dimensions[coord], pm_len(this, coord), pm_coord(this, coord));
        this->task_coords[coord] = offset[coord] + pm_coord(this, coord) * get_data_per_proc(dimensions[coord], pm_len(this, coord));
    }
}

void create_plane_types(ProcMeta* this, const size_t dimensions[DIMS]) {
    for (size_t coord = 0; coord < DIMS; ++coord) {
        size_t sizes[DIMS];
        memcpy(sizes, this->task_dims, DIMS * sizeof *sizes);
        sizes[coord] = 1;

        MPI_Type_create_subarray(DIMS, dimensions, 
            sizes, this->task_coords, MPI_ORDER_C, 
            MPI_DOUBLE, &this->plane_types[coord]);
    }
}

int pm_init(ProcMeta* this, const ProblemData* data, const int dimensions[DIMS]) {
    MPI_Comm_size(MPI_COMM_WORLD, &this->proc_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &this->proc_rank);

    if (dimensions[X] * dimensions[Y] * dimensions[Z] != this->proc_count) {
        fprintf(stderr, "Expected %d processes for [%d x %d x %d], found %d\n", 
            dimensions[X] * dimensions[Y] * dimensions[Z],
            dimensions[X], dimensions[Y], dimensions[Z], this->proc_count);
        return EXIT_FAILURE;
    }

    memcpy(this->comm_dims, dimensions, DIMS * sizeof *this->comm_dims);

    int periodic[DIMS] = {};
    MPI_Cart_create(MPI_COMM_WORLD, DIMS, this->comm_dims, periodic, 1, &this->global_comm);

    MPI_Cart_coords(this->global_comm, this->proc_rank, DIMS, this->comm_coords);
    MPI_Cart_rank(this->global_comm, this->comm_coords, &this->proc_rank);

    determine_neighbours(this);

    size_t cube_offset[DIMS];
    size_t dims_without_edges[DIMS];
    cube_without_edges(data, cube_offset, dims_without_edges);

    calculate_cuts(this, cube_offset, dims_without_edges);
    create_plane_types(this, dims_without_edges);

    return EXIT_SUCCESS;
}

int pm_neighbour(const ProcMeta* this, size_t index) {
    return this->neighbours[index];
}

MPI_Comm pm_comm(const ProcMeta* this) {
    return this->global_comm;
}

MPI_Datatype pm_plane(const ProcMeta* this, size_t plane_index) {
    return this->plane_types[s_plane_orth_coord[plane_index]];
}

int pm_coord(const ProcMeta* this, size_t coord) {
    return this->comm_coords[coord];
}

int pm_x(const ProcMeta* this) {
    return pm_coord(this, X);
}

int pm_y(const ProcMeta* this) {
    return pm_coord(this, Y);
}

int pm_z(const ProcMeta* this) {
    return pm_coord(this, Z);
}

int pm_len(const ProcMeta* this, size_t coord) {
    return this->comm_dims[coord];
}

int pm_len_x(const ProcMeta* this) {
    return pm_len(this, X);
}

int pm_len_y(const ProcMeta* this) {
    return pm_len(this, Y);
}

int pm_len_z(const ProcMeta* this) {
    return pm_len(this, Z);
}

int pm_root(const ProcMeta* this) {
    static  int coords[DIMS] = {};
    return get_cart_rank(this, coords);
}

int pm_rank(const ProcMeta* this) {
    return this->proc_rank;
}

bool pm_is_root(const ProcMeta* this) {
    return (pm_rank(this) == pm_root(this));
}

size_t pm_count(const ProcMeta* this) {
    return this->proc_count;
}
