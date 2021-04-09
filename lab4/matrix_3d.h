#ifndef LAB4_TENSOR_H
#define LAB4_TENSOR_H

#include "types.h"

void mt_init(Matrix3D* this, const int dimensions[DIMS]);

int mt_len(const Matrix3D* this, int coord);
int mt_x_len(const Matrix3D* this);
int mt_y_len(const Matrix3D* this);
int mt_z_len(const Matrix3D* this);

double mt_get(const Matrix3D* this, int x, int y, int z);
const double* mt_get_ref(const Matrix3D* this, int x, int y, int z);
double* mt_set(Matrix3D* this, int x, int y, int z);

double mt_get_ar(const Matrix3D* this, const int dimensions[DIMS]);
const double* mt_get_ref_ar(const Matrix3D* this, const int dimensions[DIMS]);
double* mt_set_ar(Matrix3D* this, const int dimensions[DIMS]);

double mt_get_by_i(const Matrix3D* this, int index);
const double* mt_get_ref_by_i(const Matrix3D* this, int index);
double* mt_set_by_i(Matrix3D* this, int index);

void mt_free(Matrix3D* this);

#endif // !LAB4_TENSOR_H