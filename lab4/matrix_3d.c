#include "matrix_3d.h"

#include <stdlib.h>
#include <string.h>

//coordinate ordering in planes orthogonal to X, Y, Z respectively:
const static size_t s_plane_coords[3][2] = {{Y, Z}, {Z, X}, {X, Y}};

const static size_t s_plane_indices[3][2] = { {NB_LEFT, NB_RIGHT}, {NB_BACK, NB_FRONT}, {NB_DOWN, NB_UP}};
const static size_t s_plane_orth_coord[NB_COUNT] = {X, X, Y, Y, Z, Z};
const static bool s_plane_right[NB_COUNT] = {false, true, false, true, false, true};

size_t get_index(const Matrix3D* this, size_t x, size_t y, size_t z) {
    return mt_y_len(this) * mt_z_len(this) * x 
         + mt_z_len(this) * y 
         + z;
}

double mt_get(const Matrix3D* this, size_t x, size_t y, size_t z) {
    return mt_get_by_i(this, get_index(this, x, y, z));
}

double mt_get_ar(const Matrix3D* this, const size_t dimensions[DIMS]) {
    return mt_get(this, dimensions[X], dimensions[Y], dimensions[Z]);
}

const double* mt_get_ref(const Matrix3D* this, size_t x, size_t y, size_t z) {
    return mt_get_ref_by_i(this, get_index(this, x, y, z));
}

const double* mt_get_ref_ar(const Matrix3D* this, const size_t dimensions[DIMS]) {
    return mt_get_ref(this, dimensions[X], dimensions[Y], dimensions[Z]);
}

double* mt_set(Matrix3D* this, size_t x, size_t y, size_t z) {
    return mt_set_by_i(this, get_index(this, x, y, z));
}

double* mt_set_ar(Matrix3D* this, const size_t dimensions[DIMS]) {
    return mt_set(this, dimensions[X], dimensions[Y], dimensions[Z]);
}

double mt_get_by_i(const Matrix3D* this, size_t index) {
    return this->data[index];
}

const double* mt_get_ref_by_i(const Matrix3D* this, size_t index) {
    return &this->data[index];
}

double* mt_set_by_i(Matrix3D* this, size_t index) {
    return &this->data[index];
}

void mt_init(Matrix3D* this, const size_t dimensions[DIMS]) {
    memcpy(this->dimensions, dimensions, DIMS * sizeof *this->dimensions);
    this->data = calloc(mt_x_len(this) * mt_y_len(this) * mt_z_len(this),  sizeof *this->data);
}

size_t mt_len(const Matrix3D* this, size_t coord) {
    return this->dimensions[coord];
}

size_t mt_x_len(const Matrix3D* this) {
    return mt_len(this, X);
}

size_t mt_y_len(const Matrix3D* this) {
    return mt_len(this, Y);
}

size_t mt_z_len(const Matrix3D* this) {
    return mt_len(this, Z);
}

double pl_get(const Plane* this, size_t x, size_t y) {
    return this->data[x * pl_y_len(this) + y];
}

const double* pl_get_ref(const Plane* this, size_t x, size_t y) {
    return &this->data[x * pl_y_len(this) + y];
}

double* pl_set(Plane* this, size_t x, size_t y) {
    return &this->data[x * pl_y_len(this) + y];
}

size_t pl_len(const Plane* this, size_t coord) {
    return this->dimensions[coord];
}

size_t pl_x_len(const Plane* this) {
    return pl_len(this, X);
}

size_t pl_y_len(const Plane* this) {
    return pl_len(this, Y);
}

void pl_init(Plane* this, const size_t dimensions[2]) {
    for (size_t i = 0; i < 2; ++i) {
        this->dimensions[i] = dimensions[i];
    }

    this->data = calloc(pl_x_len(this) * pl_y_len(this), sizeof *this->data);
}

void pl_free(Plane* this) {
    free(this->data);
}

void mt_free(Matrix3D* this) {
    free(this->data);
}
