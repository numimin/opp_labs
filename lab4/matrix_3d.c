#include "matrix_3d.h"

#include <stdlib.h>
#include <string.h>

//coordinate ordering in planes orthogonal to X, Y, Z respectively:
const static size_t s_plane_coords[3][2] = {{Y, Z}, {Z, X}, {X, Y}};

const static size_t s_plane_indices[3][2] = { {NB_LEFT, NB_RIGHT}, {NB_BACK, NB_FRONT}, {NB_DOWN, NB_UP}};
const static size_t s_plane_orth_coord[NB_COUNT] = {X, X, Y, Y, Z, Z};
const static bool s_plane_right[NB_COUNT] = {false, true, false, true, false, true};

int get_index( Matrix3D* this, int x, int y, int z) {
    return this->dimensions[Y] * this->dimensions[Z] * x + this->dimensions[Z] * y + z;
}

double ts_get( Matrix3D* this, int x, int y, int z) {
    return this->data[get_index(this, x, y, z)];
}

double ts_get_ar( Matrix3D* this,  int dimensions[DIMS]) {
    return ts_get(this, dimensions[X], dimensions[Y], dimensions[Z]);
}

double* ts_set(Matrix3D* this, int x, int y, int z) {
    return &this->data[get_index(this, x, y, z)];
}

 double* ts_set_c( Matrix3D* this, int x, int y, int z) {
    return &this->data[get_index(this, x, y, z)];
}

double* ts_set_ar(Matrix3D* this,  int dimensions[DIMS]) {
    return ts_set(this, dimensions[X], dimensions[Y], dimensions[Z]);
}

 double* ts_set_ar_c( Matrix3D* this,  int dimensions[DIMS]) {
    return ts_set_c(this, dimensions[X], dimensions[Y], dimensions[Z]);
}


void init_3d_matrix(Matrix3D* this,  int dimensions[DIMS]) {
    memcpy(this->dimensions, dimensions, DIMS * sizeof *this->dimensions);
    this->data = calloc(this->dimensions[X] * this->dimensions[Y] * this->dimensions[Z],  sizeof *this->data);
}

int ts_coord_len( Matrix3D* this, int coord) {
    return this->dimensions[coord];
}

int ts_x_len( Matrix3D* this) {
    return ts_coord_len(this, X);
}

int ts_y_len( Matrix3D* this) {
    return ts_coord_len(this, Y);
}

int ts_z_len( Matrix3D* this) {
    return ts_coord_len(this, Z);
}

double pl_get( Plane* this, int x, int y) {
    return this->data[x * this->dimensions[Y] + y];
}

int pl_coord_len( Plane* this, int coord) {
    return this->dimensions[coord];
}

int pl_x_len( Plane* this) {
    return pl_coord_len(this, X);
}

int pl_y_len( Plane* this) {
    return pl_coord_len(this, Y);
}

double vc_x( Vector* this) {
    return this->dimensions[X];
}

double vc_y( Vector* this) {
    return this->dimensions[Y];
}

double vc_z( Vector* this) {
    return this->dimensions[Z];
}

void init_plane(Plane* this,  int dimensions[DIMS], int excluded_coord,  int coords[2], bool right) {
    for (int i = 0; i < 2; ++i) {
        this->coords[i] = coords[i];
        this->dimensions[i] = dimensions[coords[i]];
    }

    this->right = right;
    this->excluded_coord = excluded_coord;

    this->data = calloc(pl_x_len(this) * pl_y_len(this), sizeof *this->data);
}

void free_plane(Plane* this) {
    free(this->data);
}

void free_3d_matrix(Matrix3D* this) {
    free(this->data);
}
