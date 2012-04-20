#ifndef _MORPH_H_
#define _MORPH_H_

SEXP morph_R (SEXP x, SEXP kernel, SEXP value, SEXP value_not, SEXP n_neighbours, SEXP n_neighbours_not, SEXP is_brush);

int is_compatible_value (const SEXP x, const R_len_t index, const SEXP include_list, const SEXP exclude_list, const int is_integer);

int is_compatible_neighbourhood (const SEXP x, const int *x_dims, const int n_dims, const int *neighbourhood_dims, const R_len_t index, const SEXP include_list, const SEXP exclude_list, const int is_integer);

void apply_kernel (const SEXP x, SEXP y, const int *x_dims, const int n_dims, const int *x_loc, const SEXP kernel, const int *kernel_dims, const int is_integer, const int is_brush);

#endif
