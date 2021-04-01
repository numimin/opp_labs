#include "problem_data.h"

double pd_h(const ProblemData* this, size_t coord) {
    return this->local_area[coord];
}

double pd_h_x(const ProblemData* this) {
    return pd_h(this, X);
}

double pd_h_y(const ProblemData* this) {
    return pd_h(this, Y);
}

double pd_h_z(const ProblemData* this) {
    return pd_h(this, Z);
}

double pd_a(const ProblemData* this) {
    return this->a;
}

double pd_factor(const ProblemData* this) {
    return 1 / (
        2 / (pd_h_x(this) * pd_h_x(this)) +
        2 / (pd_h_y(this) * pd_h_y(this)) +
        2 / (pd_h_z(this) * pd_h_z(this)) +
        pd_a(this)
    );
}

double pd_area(const ProblemData* this, size_t coord, size_t index) {
    return this->area_offset[coord] + index * pd_h(this, coord);
}

double pd_area_x(const ProblemData* this, size_t x) {
    return pd_area(this, X, x);
}

double pd_area_y(const ProblemData* this, size_t y) {
    return pd_area(this, Y, y);
}

double pd_area_z(const ProblemData* this, size_t z) {
    return pd_area(this, Z, z);
}

double pd_func(const ProblemData* this, const size_t discrete_point[DIMS], func_r3 func) {
    const double x = pd_area_x(this, discrete_point[X]);
    const double y = pd_area_y(this, discrete_point[Y]);
    const double z = pd_area_z(this, discrete_point[Z]);

    return func(x, y, z);
}

double pd_rho(const ProblemData* this, const size_t discrete_point[DIMS]) {
    return pd_func(this, discrete_point, this->rho);
}

double pd_phi(const ProblemData* this, const size_t discrete_point[DIMS]) {
    return pd_func(this, discrete_point, this->phi);
}

void count_local_area(double local_area[DIMS], double global_area[DIMS], int dimensions[DIMS]) {
    for (int i = 0; i < DIMS; ++i) {
        if (dimensions[i] == 0) 
            local_area[i] = 0;
        local_area[i] = global_area[i] / 
                             (dimensions[i] - 1);
    }
}

void pd_init(ProblemData* this) {
    this->area_offset[X] = X_0;
    this->area_offset[Y] = Y_0;
    this->area_offset[Z] = Z_0;

    this->global_area[X] = D_X;
    this->global_area[Y] = D_Y;
    this->global_area[Z] = D_Z;

    this->a = A;

    this->phi = phi;
    this->rho = rho;

    this->discrete_dimensions[X] = N_X;
    this->discrete_dimensions[Y] = N_Y;
    this->discrete_dimensions[Z] = N_Z;

    count_local_area(this->local_area, this->global_area, this->discrete_dimensions);
}

double phi(double x, double y, double z) {
    return x * x + y * y + z * z;
}

double rho(double x, double y, double z) {
    return 6 - A * phi(x, y, z);
}

size_t pd_disc(const ProblemData* this, size_t coord) {
    return this->discrete_dimensions[coord];
}

size_t pd_disc_x(const ProblemData* this) {
    return pd_disc(this, X);
}

size_t pd_disc_y(const ProblemData* this) {
    return pd_disc(this, Y);
}

size_t pd_disc_z(const ProblemData* this) {
    return pd_disc(this, Z);
}
