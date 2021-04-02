#ifndef LAB4_FUNCTIONAL_H
#define LAB4_FUNCTIONAL_H

#include <stddef.h>

#include "types.h"

#define DIMS 3

typedef void (*iterator) (void* data, size_t x, size_t y, size_t z);
typedef void (*iterator_n3) (size_t x, size_t y, size_t z);

typedef void (*pi_iterator) (PlaneIndex* index, void* data);

typedef struct {
    size_t start;
    size_t end;
} Extent;

typedef struct {
    size_t ort_coord;
    size_t plane_coords[2];
    size_t dims[DIMS];
} PlaneIndex;

size_t* pi_set_ort(PlaneIndex* this);
size_t* pi_set_coord(PlaneIndex* this, size_t coord);
size_t* pi_set_x(PlaneIndex* this);
size_t* pi_set_y(PlaneIndex* this);

size_t pi_get_ort(const PlaneIndex* this);
size_t pi_get_coord(const PlaneIndex* this, size_t coord);
size_t pi_get_x(const PlaneIndex* this);
size_t pi_get_y(const PlaneIndex* this);

size_t pi_ort_idx(const PlaneIndex* this);
size_t pi_coord_idx(const PlaneIndex* this, size_t coord);
size_t pi_x_idx(const PlaneIndex* this);
size_t pi_y_idx(const PlaneIndex* this);

typedef struct {
    PlaneIndex index;
    Extent extents[DIMS];
    size_t order[DIMS];
} PlaneIterator;

void pi_init(PlaneIterator* this, size_t ort_coord, size_t ort_index, const size_t other_coords[2], const Extent plane_extents[2]);
void pi_iterate(PlaneIterator* this, pi_iterator func, void* data);

void fn_iterate(const size_t order[DIMS], const Extent extents[DIMS], iterator func, void* data);
void fn_iterate_n3(const size_t order[DIMS], const Extent extents[DIMS], iterator_n3 func);

#endif // !LAB4_FUNCTIONAL_H