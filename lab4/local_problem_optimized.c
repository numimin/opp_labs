#include "local_problem.h"

#include "functional.h"
#include "proc_meta.h"

double get_h_2(const LocalProblem* this, size_t coord) {
    return this->constants.h_2[coord];
}

double get_h_x_2(const LocalProblem* this) {
    return get_h_2(this, X);
}

double get_h_y_2(const LocalProblem* this) {
    return get_h_2(this, Y);
}

double get_h_z_2(const LocalProblem* this) {
    return get_h_2(this, Z);
}

const Matrix3D* get_rho(const LocalProblem* this) {
    return &this->constants.rho;
}

double get_factor(const LocalProblem* this) {
    return this->constants.factor;
}

void lp_swap_cubes(LocalProblem* this) {
    size_t tmp = this->old;
    this->old = this->new;
    this->new = tmp;
}

typedef struct {
    Matrix3D* matrix;
    double factor;
} PlaneMulData;

void plane_mul(PlaneIndex* index, LocalProblem* problem) {
    *mt_set_ar(lp_new_mat(problem), index->dims) *= get_factor(problem);
}

void lp_multiply_edges(LocalProblem* this) {
    for (size_t i = 0; i < NB_COUNT; ++i) {
        pi_iterate(&this->plane_mul_iterators[i], plane_mul, this);
    }
}

typedef struct {
    LocalProblem* problem;
    Plane* plane;
} IteratorData;

void plane_finish(PlaneIndex* index, IteratorData* data) {
    *mt_set_ar(lp_new_mat(data->problem), index->dims) += 
        pl_get(data->plane, pi_get_x(index), pi_get_y(index)) / 
        data->problem->constants.h_2[pi_get_ort(index)];
}

//add evetything that the Plane can provide, doesn't multiply
void lp_finish_plane(LocalProblem* this, size_t plane_index) {
    Plane* plane = lp_in_plane(this, plane_index);
    IteratorData data = {.plane=plane, .problem=this};
    pi_iterate(&this->plane_iterators[plane_index], plane_finish, &data);
}

void calculate_next_phi(LocalProblem* this, size_t x, size_t y, size_t z) {
    const Matrix3D* old = lp_old_mat(this);
    Matrix3D* new = lp_new_mat(this);

    *mt_set(new, x, y, z) = this->constants.factor * (
        (mt_get(old, x + 1, y, z) + mt_get(old, x - 1, y, z)) / get_h_x_2(this) +
        (mt_get(old, x, y + 1, z) + mt_get(old, x, y - 1, z)) / get_h_y_2(this) +
        (mt_get(old, x, y, z + 1) + mt_get(old, x, y, z - 1)) / get_h_z_2(this) -
            mt_get(get_rho(this), x, y, z));
}

typedef struct {
    LocalProblem* problem;
    size_t neighbor_index;
} NeighborAddData;

void plane_neighbor_add(PlaneIndex* index, NeighborAddData* data) {
    const size_t coord_index = pi_get_ort(index);

    *pi_set_ort(index) = data->neighbor_index;
    const double res = mt_get_ar(lp_old_mat(data), index->dims) / get_h_2(data->problem, coord_index);

    *pi_set_ort(index) = coord_index;
    *mt_set_ar(lp_new_mat(data->problem), index->dims) = res;
}

void plane_inner_add(PlaneIndex* index, LocalProblem* problem) {
    for (size_t coord = 0; coord < 2; ++coord) {
        double res = 0;

        *pi_set_coord(index, coord) += 1;
        res += mt_get_ar(lp_old_mat(problem), index->dims);

        *pi_set_coord(index, coord) -= 2;
        res += mt_get_ar(lp_old_mat(problem), index->dims);
        res /= get_h_2(problem, pi_coord_idx(index, coord));

        *pi_set_coord(index, coord) += 1;
        *mt_set_ar(lp_new_mat(problem), index->dims) += res;
    }

    *mt_set_ar(lp_new_mat(problem), index->dims) -= mt_get_ar(get_rho(problem), index->dims);
}

//iterates over everything that exists
void lp_iterate(LocalProblem* this) {
    const size_t order[DIMS] = {X, Y, Z};
    const Extent extents[DIMS] = {
        {.start=1, .end=mt_x_len(lp_old_mat(this)) - 1},
        {.start=1, .end=mt_y_len(lp_old_mat(this)) - 1},
        {.start=1, .end=mt_z_len(lp_old_mat(this)) - 1},
    };

    fn_iterate_n3(order, extents, calculate_next_phi);

    for (size_t i = 0; i < NB_COUNT; ++i) {
        NeighborAddData data = {.neighbor_index=this->plane_neighbor_indices[i], .problem=this};
        pi_iterate(&this->plane_iterators[i], plane_neighbor_add, &data);
        pi_iterate(&this->plane_edged_iterators[i], plane_inner_add, this);
    }
}

typedef struct {
    Matrix3D* rho;
    const ProblemData* data;
    const ProcMeta* meta;
} RhoInitData;

void fill_rho(RhoInitData* data, size_t x, size_t y, size_t z) {
    size_t discrete_point[DIMS];
    discrete_point[X] = pm_disc_x(data->meta, x);
    discrete_point[Y] = pm_disc_y(data->meta, y);
    discrete_point[Z] = pm_disc_z(data->meta, z);
    *mt_set(data, x, y, z) = pd_rho(data->data, discrete_point);
}

void init_constants(Constants* this, const ProblemData* data, const ProcMeta* meta) {
    this->factor = pd_factor(data);

    this->h_2[X] = pd_h_x(data) * pd_h_x(data);
    this->h_2[Y] = pd_h_y(data) * pd_h_y(data);
    this->h_2[Z] = pd_h_z(data) * pd_h_z(data);

    mt_init(&this->rho, meta->task_dims);

    const size_t order[DIMS] = {X, Y, Z};
    const Extent extents[DIMS] = {
        {.start=0, .end=mt_x_len(&this->rho)},
        {.start=0, .end=mt_y_len(&this->rho)},
        {.start=0, .end=mt_z_len(&this->rho)},
    };

    RhoInitData rho_data = {.rho=&this->rho, .data=data, .meta=meta};
    fn_iterate(order, extents, fill_rho, &rho_data);
}

void lp_init(LocalProblem* this, const ProblemData* data, const ProcMeta* meta) {
    init_constants(&this->constants, data, meta);

    for (int i = 0; i < 2; ++i) {
        mt_init(&this->slice[i], meta->task_dims[i]);
    }

    this->old = 0;
    this->new = 1;

    const static bool shrink[3][2] = {{false, false}, {false, true}, {true, true}};

    for (size_t coord = 0; coord < DIMS; ++coord) {
        Extent plane_extents[2];
        Extent edged_plane_extents[2];
        Extent mul_plane_extents[2];

        for (size_t i = 0; i < 2; ++i) {
            plane_extents[i].start = 0;
            plane_extents[i].end = meta->task_dims[s_plane_coords[coord][i]];

            edged_plane_extents[i].start = 1;
            edged_plane_extents[i].start = meta->task_dims[s_plane_coords[coord][i]] - 1;

            mul_plane_extents[i].start = shrink[coord][i] ? edged_plane_extents[i].start : plane_extents[i].start;
            mul_plane_extents[i].end = shrink[coord][i] ? edged_plane_extents[i].end : plane_extents[i].end;
        }

        const size_t plane_dims[2] = {plane_extents[0].end, plane_extents[1].end};

        for (bool right = true; right; right = !right) {
            const size_t plane_index = s_plane_indices[coord][right];
            pl_init(&this->neighbours[plane_index], plane_dims);

            const size_t coord_index = right ? meta->task_dims[coord] - 1 : 0;
            pi_init(&this->plane_iterators[plane_index], 
                coord, coord_index,
                s_plane_coords[coord], plane_extents);
            pi_init(&this->plane_edged_iterators[plane_index], 
                coord, coord_index,
                s_plane_coords[coord], edged_plane_extents);
            pi_init(&this->plane_mul_iterators[plane_index], 
                coord, coord_index,
                s_plane_coords[coord], edged_plane_extents);
        }
    }
}

Matrix3D* lp_old_mat(LocalProblem* this) {
    return &this->slice[this->old];
}

Matrix3D* lp_new_mat(LocalProblem* this) {
    return &this->slice[this->new];
}

Plane* lp_in_plane(LocalProblem* this, size_t index) {
    return &this->neighbours[index];
}

double* lp_out_plane(LocalProblem* this, size_t index) {
     int coord = s_plane_orth_coord[index];
     int coord_index = s_plane_right[index] ? mt_len(lp_old_mat(this), coord) : 0;

    int dims[DIMS] = {};
    dims[coord] = coord_index;
    return mt_set_ar(lp_new_mat(this), dims);
}
