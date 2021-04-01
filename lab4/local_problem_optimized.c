#include "local_problem.h"

void coords_without(int excluded, int coords[2]) {
    int coord_idx = 0;
    for (int i = 0; i < 3; ++i) {
        if (excluded != i) {
            coords[coord_idx++] = i;
        }
    }
}

double get_h_2( ProblemSolver* this, int coord) {
    return this->h_2[coord];
}

double get_h_x_2( ProblemSolver* this) {
    return get_h_2(this, X);
}

double get_h_y_2( ProblemSolver* this) {
    return get_h_2(this, Y);
}

double get_h_z_2( ProblemSolver* this) {
    return get_h_2(this, Z);
}

//iterates everything that knows but doesn't add what doesn't know and multiply by factor;
void iterate_plane( ProblemSolver* solver,  Matrix3D* old, Matrix3D* new, int coord, bool right) {
    int coords[2];
    coords_without(coord, coords);

     int coord_index = right ? mt_len(old, coord) : 0;
     int neighbor_index = right ? mt_len(old, coord) - 1 : 1;

    int dims[DIMS] = {};
    dims[coord] = coord_index;

    for (dims[coords[0]] = 0; dims[coords[0]] < mt_len(old, coords[0]); ++dims[coords[0]]) {
        for (dims[coords[1]] = 0; dims[coords[1]] < mt_len(old, coords[1]); ++dims[coords[1]]) {
            dims[coord] = neighbor_index;
            double res = mt_get_ar(old, dims) / get_h_2(solver, coord);
            dims[coord] = coord_index;
            *mt_set_ar(new, dims) = res;
        }
    }
    
    for (dims[coords[0]] = 1; dims[coords[0]] < mt_len(old, coords[0]) - 1; ++dims[coords[0]]) {
        for (dims[coords[1]] = 1; dims[coords[1]] < mt_len(old, coords[1]) - 1; ++dims[coords[1]]) {
            for (int i = 0; i < 2; ++i) {
                double res = 0;

                dims[coords[i]] += 1;
                res += mt_get_ar(old, dims);

                dims[coords[i]] -= 2;
                res += mt_get_ar(old, dims) / get_h_2(solver, coords[i]);

                dims[coords[i]] += 1;
                res /= get_h_2(solver, coords[i]);
                *mt_set_ar(new, dims) += res;
            }

            *mt_set_ar(new, dims) -= mt_get_ar(&solver->rho, dims);
        }
    }
}

void ps_swap_cubes(LocalProblem* this) {
    int tmp = this->old;
    this->old = this->new;
    this->new = tmp;
}

//add evetything that the Plane can provide, doesn't multiply
void finish_iteration( ProblemSolver* solver,  Plane* old, Matrix3D* new) {
    int dims[DIMS] = {};
    dims[old->excluded_coord] = old->right ? mt_len(new, old->excluded_coord) : 0;
    for (dims[old->coords[X]] = 0; dims[old->coords[X]] < mt_len(new, old->coords[X]); ++dims[old->coords[X]]) {
        for (dims[old->coords[Y]] = 0; dims[old->coords[Y]] < mt_len(new, old->coords[Y]); ++dims[old->coords[Y]]) {
            *mt_set_ar(new, dims) += pl_get(old, dims[old->coords[X]], dims[old->coords[Y]]) / get_h_2(solver, old->excluded_coord);
        }
    }
}

void multiply_plane(Matrix3D* new, double factor, int coord, bool right, int coords[2], bool shrink[2]) {
    int start[2];
    int end[2];
    for (int i = 0; i < 2; ++i) {
        start[i] = shrink[i] ? 1 : 0;
        end[i] = shrink[i] ? mt_len(new, coords[i]) - 1 : mt_len(new, coords[i]);
    }

    int dims[DIMS] = {};
    dims[coord] = right ? mt_len(new, coord) : 0;
    for (dims[coords[0]] = start[0]; dims[coords[0]] < end[0]; ++dims[coords[0]]) {
        for (dims[coords[1]] = start[1]; dims[coords[1]] < end[1]; ++dims[coords[1]]) {
            *mt_set_ar(new, dims) *= factor;
        }
    }
}

void ps_multiply_edges(Matrix3D* new, double factor) {
     static bool shrink[3][2] = {{false, false}, {false, true}, {true, true}};

    for (int coord = 0; coord < DIMS; ++coord) {
        for (bool right = true; right; right = !right) {
            multiply_plane(new, factor, coord, right, s_plane_coords[coord], shrink[coord]);
        }
    }
}

void ps_finish_plane( ProblemSolver* solver, LocalProblem* this, int plane_index) {
    finish_iteration(solver, get_in_plane_c(this, plane_index), get_new(this));
}

//iterates over everything that exists
void iterate( ProblemSolver* this,  Matrix3D* old, Matrix3D* new) {
    for (int x = 1; x < mt_x_len(old) - 1; ++x) {
        for (int y = 1; y < mt_y_len(old) - 1; ++y) {
            for (int z = 1; z < mt_z_len(old) - 1; ++z) {
                *mt_set(new, x, y, z) = this->factor * (
                    (mt_get(old, x + 1, y, z) + mt_get(old, x - 1, y, z)) / get_h_x_2(this) +
                    (mt_get(old, x, y + 1, z) + mt_get(old, x, y - 1, z)) / get_h_y_2(this) +
                    (mt_get(old, x, y, z + 1) + mt_get(old, x, y, z - 1)) / get_h_z_2(this) -
                     mt_get(&this->rho, x, y, z));
            }
        }
    }

    for (int coord = 0; coord < DIMS; ++coord) {
        for (bool right = true; right; right = !right) {
            iterate_plane(this, old, new, coord, right);
        }
    }
}

void ps_iterate( ProblemSolver* solver, LocalProblem* problem) {
    iterate(solver, get_old(problem), get_new(problem));
}

void init_solver(ProblemSolver* this,  LocalProblemData* local_task,  ProblemData* data) {
    this->factor = pd_factor(data);
    mt_init(&this->rho, local_task->dimensions);
    get_rho(data, local_task, &this->rho);

    this->h_2[X] = pd_h_x(data) * pd_h_x(data);
    this->h_2[Y] = pd_h_y(data) * pd_h_y(data);
    this->h_2[Z] = pd_h_z(data) * pd_h_z(data);
}

Matrix3D* get_old(LocalProblem* this) {
    return &this->slice[this->old];
}

 Matrix3D* get_old_c( LocalProblem* this) {
    return &this->slice[this->old];
}

Matrix3D* get_new(LocalProblem* this) {
    return &this->slice[this->new];
}

 Matrix3D* get_new_c( LocalProblem* this) {
    return &this->slice[this->new];
}

Plane* get_in_plane(LocalProblem* this, int index) {
    return &this->neighbours[index];
}

 Plane* get_in_plane_c( LocalProblem* this, int index) {
    return &this->neighbours[index];
}

double* get_out_plane(LocalProblem* this, int index) {
     int coord = s_plane_orth_coord[index];
     int coord_index = s_plane_right[index] ? mt_len(get_old(this), coord) : 0;

    int dims[DIMS] = {};
    dims[coord] = coord_index;
    return mt_set_ar(get_new(this), dims);
}

 double* get_out_plane_c( LocalProblem* this, int index) {
     int coord = s_plane_orth_coord[index];
     int coord_index = s_plane_right[index] ? mt_len(get_old(this), coord) : 0;

    int dims[DIMS] = {};
    dims[coord] = coord_index;
    return ts_set_ar_c(get_new(this), dims);
}

void init_local_problem(LocalProblem* this,  LocalProblemData* local_task) {
    for (int i = 0; i < 2; ++i) {
        mt_init(&this->slice[i], local_task->dimensions);
    }

    this->old = 0;
    this->new = 1;

    for (int coord = 0; coord < 3; ++coord) {
        for (bool right = true; right; right = !right) {
            pl_init(&this->neighbours[s_plane_indices[coord][right]], 
                local_task->dimensions, coord, s_plane_coords[coord], right);
        }
    }
}
