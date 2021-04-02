#include "functional.h"

#include <string.h>

void fn_iterate(const size_t order[DIMS], const Extent extents[DIMS], iterator func, void* data) {
    size_t vec[DIMS]; //(X, Y, Z)
    for (vec[order[0]] = extents[0].start; vec[order[0]] < extents[0].end; ++vec[order[0]]) {
        for (vec[order[1]] = extents[1].start; vec[order[1]] < extents[1].end; ++vec[order[1]]) {
            for (vec[order[2]] = extents[2].start; vec[order[2]] < extents[2].end; ++vec[order[2]]) {
                func(data, vec[X], vec[Y], vec[Z]);
            }
        }
    }
}

typedef struct {
    iterator_n3 func;
} FuncWrapper; 

void func_wrapper(FuncWrapper* wrapper, size_t x, size_t y, size_t z) {
    wrapper->func(x, y, z);
};

void fn_iterate_n3(const size_t order[DIMS], const Extent extents[DIMS], iterator_n3 func) {
    FuncWrapper wrapper = {.func=func};
    fn_iterate(order, extents, &wrapper, NULL);
}

void pi_init(PlaneIterator* this, size_t ort_coord, size_t ort_index, const size_t other_coords[2], const Extent plane_extents[2]) {
    this->index.ort_coord = ort_coord;
    memcpy(this->index.plane_coords, other_coords, 2 * sizeof *this->index.plane_coords);

    this->order[X] = ort_coord;
    this->order[Y] = other_coords[X];
    this->order[Z] = other_coords[Y];

    for (size_t i = 0; i < 2; ++i) {
        memcpy(&this->extents[i + 1], &plane_extents[i], sizeof(Extent));
    }

    this->extents[0].start = ort_index;
    this->extents[0].end = ort_index + 1;
}

typedef struct {
    void* data;
    pi_iterator func;
    PlaneIndex* index;
} PiWrapper;

void pi_wrapper(PiWrapper* wrapper, size_t x, size_t y, size_t z) {
    wrapper->index->dims[X] = x;
    wrapper->index->dims[Y] = y;
    wrapper->index->dims[Z] = z;
    wrapper->func(wrapper->index, wrapper->data);
}

void pi_iterate(PlaneIterator* this, pi_iterator func, void* data) {
    PiWrapper wrapper = {.data=data, .func=func, .index=&this->index};
    fn_iterate(this->order, this->extents, pi_wrapper, &wrapper);
}

size_t* pi_set(PlaneIndex* this, size_t index) {
    return &this->dims[index];
}

size_t* pi_set_ort(PlaneIndex* this) {
    return pi_set(this, pi_ort_idx(this));
}

size_t* pi_set_coord(PlaneIndex* this, size_t coord) {
    return pi_set(this, pi_coord_idx(this, coord));
}

size_t* pi_set_x(PlaneIndex* this) {
    return pi_set_coord(this, X);
}

size_t* pi_set_y(PlaneIndex* this) {
    return pi_set_coord(this, Y);
}

size_t pi_get(PlaneIndex* this, size_t index) {
    return this->dims[index];
}

size_t pi_get_ort(const PlaneIndex* this) {
    return pi_get(this, pi_ort_idx(this));
}

size_t pi_get_coord(const PlaneIndex* this, size_t coord) {
    return pi_get(this, pi_coord_idx(this, coord));
}

size_t pi_get_x(const PlaneIndex* this) {
    return pi_get_coord(this, X);
}

size_t pi_get_y(const PlaneIndex* this) {
    return pi_get_coord(this, Y);
}

size_t pi_ort_idx(const PlaneIndex* this) {
    return this->ort_coord;
}

size_t pi_coord_idx(const PlaneIndex* this, size_t coord) {
    return this->plane_coords[coord];
}

size_t pi_x_idx(const PlaneIndex* this) {
    return pi_coord_idx(this, X);
}

size_t pi_y_idx(const PlaneIndex* this) {
    return pi_coord_idx(this, Y);
}
