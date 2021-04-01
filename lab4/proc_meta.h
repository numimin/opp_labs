#ifndef LAB4_PROC_META_H
#define LAB4_PROC_META_H

#include "types.h"
#include "local_problem.h"

#define NO_NEIGHBOUR (-1)

int pm_init(ProcMeta* this, const ProblemData* data, const int dimensions[DIMS]);

int pm_neighbour(const ProcMeta* this, size_t index);
MPI_Comm pm_comm(const ProcMeta* this);

MPI_Datatype pm_plane(const ProcMeta* this, size_t plane_index);

int pm_coord(const ProcMeta* this, size_t coord);
int pm_x(const ProcMeta* this);
int pm_y(const ProcMeta* this);
int pm_z(const ProcMeta* this);

int pm_len(const ProcMeta* this, size_t coord);
int pm_len_x(const ProcMeta* this);
int pm_len_y(const ProcMeta* this);
int pm_len_z(const ProcMeta* this);

size_t pm_disc(const ProcMeta* this, size_t coord, size_t local_coord);
size_t pm_disc_x(const ProcMeta* this, size_t local_x);
size_t pm_disc_y(const ProcMeta* this, size_t local_y);
size_t pm_disc_z(const ProcMeta* this, size_t local_z);

size_t pm_disc_len(const ProcMeta* this, size_t coord);
size_t pm_disc_len_x(const ProcMeta* this);
size_t pm_disc_len_y(const ProcMeta* this);
size_t pm_disc_len_z(const ProcMeta* this);

int pm_root(const ProcMeta* this);
int pm_rank(const ProcMeta* this);
bool pm_is_root(const ProcMeta* this);
size_t pm_count(const ProcMeta* this);

#endif // !LAB4_PROC_META_H