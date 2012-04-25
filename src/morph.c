#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

#include "morph.h"
#include "index.h"

SEXP morph_R (SEXP x, SEXP kernel, SEXP value, SEXP value_not, SEXP n_neighbours, SEXP n_neighbours_not, SEXP is_brush, SEXP is_eraser)
{
    R_len_t i, ii;
    SEXP y;
    int j;
    double kernel_sum = 0.0;
    
    SEXP dims = getAttrib(x, R_DimSymbol);
    int n_dims = (int) LENGTH(dims);
    int kernel_width = INTEGER(getAttrib(kernel, R_DimSymbol))[0];
    int kernel_centre = (kernel_width - 1) / 2;
    int integer_x = (int) IS_INTEGER(x);
    R_len_t len = LENGTH(x);
    R_len_t kernel_len = LENGTH(kernel);
    R_len_t neighbourhood_len = ((R_len_t) R_pow_di(3.0, n_dims)) - 1;
    
    int *x_dims = (int *) R_alloc(n_dims, sizeof(int));
    for (j=0; j<n_dims; j++)
        x_dims[j] = INTEGER(dims)[j];
    
    int *kernel_dims = (int *) R_alloc(n_dims, sizeof(int));
    for (j=0; j<n_dims; j++)
        kernel_dims[j] = kernel_width;
    
    int *neighbourhood_dims = (int *) R_alloc(n_dims, sizeof(int));
    for (j=0; j<n_dims; j++)
        neighbourhood_dims[j] = 3;
    
    int *temp = (int *) R_alloc(2*n_dims, sizeof(int));
    
    int_or_double_ptr kernel_p;
    if (integer_x)
        kernel_p.i = INTEGER(kernel);
    else
        kernel_p.d = REAL(kernel);
    
    int *kernel_matrix_locs = (int *) R_alloc(kernel_len*n_dims, sizeof(int));
    for (i=0; i<kernel_len; i++)
    {
        vector_to_matrix_loc(i, kernel_dims, n_dims, temp);
        for (j=0; j<n_dims; j++)
            kernel_matrix_locs[i + (j*kernel_len)] = temp[j] - kernel_centre;
        
        if (integer_x)
            kernel_sum += (double) kernel_p.i[i];
        else
            kernel_sum += kernel_p.d[i];
    }
    
    int *neighbourhood_matrix_locs = (int *) R_alloc(neighbourhood_len*n_dims, sizeof(int));
    for (i=0; i<neighbourhood_len; i++)
    {
        // The neighbourhood does not include the centre point itself
        ii = (i >= neighbourhood_len/2) ? i + 1 : i;
        
        vector_to_matrix_loc(ii, neighbourhood_dims, n_dims, temp);
        for (j=0; j<n_dims; j++)
            neighbourhood_matrix_locs[i + (j*neighbourhood_len)] = temp[j] - 1;
    }
    
    int *loc = (int *) R_alloc(n_dims, sizeof(int));
    
    if (integer_x)
    {
        PROTECT(y = NEW_INTEGER(len));
        memcpy(INTEGER(y), INTEGER(x), ((size_t) len)*sizeof(int));
    }
    else
    {
        PROTECT(y = NEW_NUMERIC(len));
        memcpy(REAL(y), REAL(x), ((size_t) len)*sizeof(double));
    }
    
    for (i=0; i<len; i++)
    {
        if (is_compatible_value(x, i, value, value_not, integer_x) && is_compatible_neighbourhood(x, x_dims, n_dims, neighbourhood_len, neighbourhood_matrix_locs, i, n_neighbours, n_neighbours_not, integer_x, temp))
        {
            vector_to_matrix_loc((size_t) i, x_dims, n_dims, loc);
            apply_kernel(x, y, x_dims, n_dims, loc, kernel, kernel_len, kernel_matrix_locs, kernel_sum, integer_x, INTEGER(is_brush)[0], INTEGER(is_eraser)[0], temp);
        }
    }
    
    UNPROTECT(1);
    return y;
}

int is_compatible_value (const SEXP x, const R_len_t index, const SEXP include_list, const SEXP exclude_list, const int is_integer)
{
    R_len_t i;
    
    if (LENGTH(include_list) > 0)
    {
        for (i=0; i<LENGTH(include_list); i++)
        {
            if (is_integer && INTEGER(x)[index] == (int) REAL(include_list)[i])
                return TRUE;
            else if (!is_integer && REAL(x)[index] == REAL(include_list)[i])
                return TRUE;
        }
        return FALSE;
    }
    
    if (LENGTH(exclude_list) > 0)
    {
        for (i=0; i<LENGTH(exclude_list); i++)
        {
            if (is_integer && INTEGER(x)[index] == (int) REAL(exclude_list)[i])
                return FALSE;
            else if (!is_integer && REAL(x)[index] == REAL(exclude_list)[i])
                return FALSE;
        }
        return TRUE;
    }
    
    return TRUE;
}

int is_compatible_neighbourhood (const SEXP x, const int *x_dims, const int n_dims, const int neighbourhood_len, const int *neighbourhood_matrix_locs, const R_len_t index, const SEXP include_list, const SEXP exclude_list, const int is_integer, int *temp)
{
    if (LENGTH(include_list) == 0 && LENGTH(exclude_list) == 0)
        return TRUE;
    
    size_t i, l;
    int j, n_neighbours = 0;
    
    int_or_double_ptr x_p;
    if (is_integer)
        x_p.i = INTEGER(x);
    else
        x_p.d = REAL(x);
    
    vector_to_matrix_loc((size_t) index, x_dims, n_dims, temp+3);
    
    for (i=0; i<neighbourhood_len; i++)
    {
        for (j=0; j<n_dims; j++)
            temp[j] = temp[j+3] + neighbourhood_matrix_locs[i + (j*neighbourhood_len)];
        
        if (loc_in_bounds(temp, x_dims, n_dims))
        {
            matrix_to_vector_loc(temp, x_dims, n_dims, &l);
        
            if (is_integer && x_p.i[l] != 0)
                n_neighbours++;
            else if (!is_integer && x_p.d[l] != 0.0)
                n_neighbours++;
        }
    }
    
    if (LENGTH(include_list) > 0)
    {
        for (i=0; i<LENGTH(include_list); i++)
        {
            if (n_neighbours == INTEGER(include_list)[i])
                return TRUE;
        }
        return FALSE;
    }
    
    if (LENGTH(exclude_list) > 0)
    {
        for (i=0; i<LENGTH(exclude_list); i++)
        {
            if (n_neighbours == INTEGER(exclude_list)[i])
                return FALSE;
        }
        return TRUE;
    }
    
    return TRUE;
}

void apply_kernel (const SEXP x, SEXP y, const int *x_dims, const int n_dims, const int *x_loc, const SEXP kernel, const R_len_t kernel_len, const int *kernel_matrix_locs, const double kernel_sum, const int is_integer, const int is_brush, const int is_eraser, int *temp)
{
    size_t i, l;
    int j;
    
    int_or_double_ptr x_p, y_p, kernel_p;
    if (is_integer)
    {
        x_p.i = INTEGER(x);
        y_p.i = INTEGER(y);
        kernel_p.i = INTEGER(kernel);
    }
    else
    {
        x_p.d = REAL(x);
        y_p.d = REAL(y);
        kernel_p.d = REAL(kernel);
    }
    
    double value = 0, visited_kernel_sum = 0;
    
    for (i=0; i<kernel_len; i++)
    {
        for (j=0; j<n_dims; j++)
            temp[j] = x_loc[j] + kernel_matrix_locs[i + (j*kernel_len)];
        
        if (loc_in_bounds(temp, x_dims, n_dims))
        {
            matrix_to_vector_loc(temp, x_dims, n_dims, &l);
            
            if (is_brush)
            {
                if (is_integer)
                {
                    if (kernel_p.i[i] == NA_INTEGER)
                        y_p.i[l] = x_p.i[l];
                    else if (is_eraser)
                        y_p.i[l] = kernel_p.i[i] == 0 ? x_p.i[l] : 0;
                    else
                        y_p.i[l] = kernel_p.i[i];
                }
                else
                {
                    if (ISNA(kernel_p.d[i]))
                        y_p.d[l] = x_p.d[l];
                    else if (is_eraser)
                        y_p.d[l] = kernel_p.d[i] == 0.0 ? x_p.d[l] : 0.0;
                    else
                        y_p.d[l] = kernel_p.d[i];
                }
            }
            else
            {
                if (is_integer)
                {
                    value += (double) kernel_p.i[i] * x_p.i[l];
                    visited_kernel_sum += (double) kernel_p.i[i];
                }
                else
                {
                    value += kernel_p.d[i] * x_p.d[l];
                    visited_kernel_sum += kernel_p.d[i];
                }
            }
        }
    }
    
    if (!is_brush)
    {
        matrix_to_vector_loc(x_loc, x_dims, n_dims, &l);
        if (is_integer)
            y_p.i[l] = (int) round(value * kernel_sum / visited_kernel_sum);
        else
            y_p.d[l] = value * kernel_sum / visited_kernel_sum;
    }
}
