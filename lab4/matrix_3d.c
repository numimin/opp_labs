#include "matrix_3d.h"

#include <stdlib.h>
#include <string.h>

int get_index(const Matrix3D* this, int x, int y, int z) {
    return mt_y_len(this) * mt_z_len(this) * x 
         + mt_z_len(this) * y 
         + z;
}

double mt_get(const Matrix3D* this, int x, int y, int z) {
    return mt_get_by_i(this, get_index(this, x, y, z));
}

double mt_get_ar(const Matrix3D* this, const int dimensions[DIMS]) {
    return mt_get(this, dimensions[X], dimensions[Y], dimensions[Z]);
}

const double* mt_get_ref(const Matrix3D* this, int x, int y, int z) {
    return mt_get_ref_by_i(this, get_index(this, x, y, z));
}

const double* mt_get_ref_ar(const Matrix3D* this, const int dimensions[DIMS]) {
    return mt_get_ref(this, dimensions[X], dimensions[Y], dimensions[Z]);
}

double* mt_set(Matrix3D* this, int x, int y, int z) {
    return mt_set_by_i(this, get_index(this, x, y, z));
}

double* mt_set_ar(Matrix3D* this, const int dimensions[DIMS]) {
    return mt_set(this, dimensions[X], dimensions[Y], dimensions[Z]);
}

double mt_get_by_i(const Matrix3D* this, int index) {
    return this->data[index];
}

const double* mt_get_ref_by_i(const Matrix3D* this, int index) {
    return &this->data[index];
}

double* mt_set_by_i(Matrix3D* this, int index) {
    return &this->data[index];
}

void mt_init(Matrix3D* this, const int dimensions[DIMS]) {
    memcpy(this->dimensions, dimensions, DIMS * sizeof *this->dimensions);
    this->data = calloc(mt_x_len(this) * mt_y_len(this) * mt_z_len(this),  sizeof *this->data);
}

int mt_len(const Matrix3D* this, int coord) {
    return this->dimensions[coord];
}

int mt_x_len(const Matrix3D* this) {
    return mt_len(this, X);
}

int mt_y_len(const Matrix3D* this) {
    return mt_len(this, Y);
}

int mt_z_len(const Matrix3D* this) {
    return mt_len(this, Z);
}

void mt_free(Matrix3D* this) {
    free(this->data);
}
