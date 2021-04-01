#ifndef LAB4_TENSOR_H
#define LAB4_TENSOR_H

#include "types.h"

extern size_t s_plane_coords[DIMS][2];

extern size_t s_plane_indices[DIMS][2];
extern size_t s_plane_orth_coord[NB_COUNT];
extern bool s_plane_right[NB_COUNT];

void init_3d_matrix(Matrix3D* this,  int dimensions[DIMS]);

int ts_coord_len( Matrix3D* this, int coord);
int ts_x_len( Matrix3D* this);
int ts_y_len( Matrix3D* this);
int ts_z_len( Matrix3D* this);

double ts_get( Matrix3D* this, int x, int y, int z);
double* ts_set(Matrix3D* this, int x, int y, int z);

 double* ts_set_c( Matrix3D* this, int x, int y, int z);

double ts_get_ar( Matrix3D* this,  int dimensions[DIMS]);
double* ts_set_ar(Matrix3D* this,  int dimensions[DIMS]);

 double* ts_set_ar_c( Matrix3D* this,  int dimensions[DIMS]);

void init_plane(Plane* this,  int dimensions[DIMS], int excluded_coord,  int coords[2], bool right);

double pl_get( Plane* this, int x, int y);

int pl_coord_len( Plane* this, int coord);
int pl_x_len( Plane* this);
int pl_y_len( Plane* this);

double vc_x( Vector* this);
double vc_y( Vector* this);
double vc_z( Vector* this);

void free_plane(Plane* this);
void free_3d_matrix(Matrix3D* this);

#endif // !LAB4_TENSOR_H