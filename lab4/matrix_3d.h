#ifndef LAB4_TENSOR_H
#define LAB4_TENSOR_H

#include "types.h"

extern size_t s_plane_coords[DIMS][2];

extern size_t s_plane_indices[DIMS][2];
extern size_t s_plane_orth_coord[NB_COUNT];
extern bool s_plane_right[NB_COUNT];

void mt_init(Matrix3D* this, const size_t dimensions[DIMS]);

size_t mt_len(const Matrix3D* this, size_t coord);
size_t mt_x_len(const Matrix3D* this);
size_t mt_y_len(const Matrix3D* this);
size_t mt_z_len(const Matrix3D* this);

double mt_get(const Matrix3D* this, size_t x, size_t y, size_t z);
const double* mt_get_ref(const Matrix3D* this, size_t x, size_t y, size_t z);
double* mt_set(Matrix3D* this, size_t x, size_t y, size_t z);

double mt_get_ar(const Matrix3D* this, const size_t dimensions[DIMS]);
const double* mt_get_ref_ar(const Matrix3D* this, const size_t dimensions[DIMS]);
double* mt_set_ar(Matrix3D* this, const size_t dimensions[DIMS]);

void pl_init(Plane* this, const size_t dimensions[DIMS], size_t excluded_coord, const size_t coords[2], bool right);

double pl_get(const Plane* this, size_t x, size_t y);
const double* pl_get_ref(const Plane* this, size_t x, size_t y);
double* pl_set(Plane* this, size_t x, size_t y);

size_t pl_len(const Plane* this, size_t coord);
size_t pl_x_len(const Plane* this);
size_t pl_y_len(const Plane* this);

void pl_free(Plane* this);
void mt_free(Matrix3D* this);

#endif // !LAB4_TENSOR_H