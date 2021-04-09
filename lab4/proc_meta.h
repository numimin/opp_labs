#ifndef LAB4_PROC_META_H
#define LAB4_PROC_META_H

#include "types.h"
#include "local_problem.h"

#define NO_NEIGHBOUR (-1)

void pm_init(ProcMeta* this, const ProblemData* data);
void pm_free(ProcMeta* this);

int pm_count(const ProcMeta* this);
int pm_coord(const ProcMeta* this);
int pm_rank(const ProcMeta* this);
int pm_neighbour(const ProcMeta* this, int index);

MPI_Comm pm_comm(const ProcMeta* this);

int pm_disc(const ProcMeta* this, int coord, int local_coord);
int pm_disc_x(const ProcMeta* this, int local_x);
int pm_disc_y(const ProcMeta* this, int local_y);
int pm_disc_z(const ProcMeta* this, int local_z);

int pm_disc_len(const ProcMeta* this, int coord);
int pm_disc_len_x(const ProcMeta* this);
int pm_disc_len_y(const ProcMeta* this);
int pm_disc_len_z(const ProcMeta* this);

int pm_root(const ProcMeta* this);
bool pm_is_root(const ProcMeta* this);

#endif // !LAB4_PROC_META_H