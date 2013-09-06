#ifndef _MORPH_H_
#define _MORPH_H_

typedef union {
    int *i;
    double *d;
} int_or_double_ptr;

SEXP morph_R (SEXP x, SEXP kernel, SEXP value, SEXP value_not, SEXP n_neighbours, SEXP n_neighbours_not, SEXP operator, SEXP merge);

int is_compatible_value (const int_or_double_ptr x_p, const R_len_t index, const SEXP include_list, const SEXP exclude_list, const int is_integer);

int is_compatible_neighbourhood (const int_or_double_ptr x_p, const int *x_dims, const int n_dims, const int neighbourhood_len, const int *neighbourhood_matrix_locs, const R_len_t index, const SEXP include_list, const SEXP exclude_list, const int is_integer, int *temp);

void apply_kernel (const int_or_double_ptr x_p, int_or_double_ptr y_p, const int *x_dims, const int n_dims, const int *x_loc, const int_or_double_ptr kernel_p, const R_len_t kernel_len, const int *kernel_matrix_locs, const double kernel_sum, const int is_integer, const char *operator, const char *merge, int *temp, double *values);

#endif
