#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

#include "morph.h"
#include "index.h"

SEXP morph_R (SEXP x, SEXP kernel, SEXP value, SEXP value_not, SEXP n_neighbours, SEXP n_neighbours_not, SEXP is_brush, SEXP is_eraser)
{
    R_len_t i;
    int j, n;
    SEXP y;
    
    SEXP dims = getAttrib(x, R_DimSymbol);
    int n_dims = (int) LENGTH(dims);
    int kernel_width = INTEGER(getAttrib(kernel, R_DimSymbol))[0];
    int kernel_centre = (kernel_width - 1) / 2;
    int integer_x = (int) IS_INTEGER(x);
    R_len_t len = LENGTH(x);
    R_len_t kernel_len = LENGTH(kernel);
    R_len_t neighbourhood_len = (R_len_t) R_pow_di(3.0, n_dims);
    
    int *x_dims = (int *) R_alloc(n_dims, sizeof(int));
    for (j=0; j<n_dims; j++)
        x_dims[j] = INTEGER(dims)[j];
    
    int *kernel_dims = (int *) R_alloc(n_dims, sizeof(int));
    for (j=0; j<n_dims; j++)
        kernel_dims[j] = kernel_width;
    
    int *neighbourhood_dims = (int *) R_alloc(n_dims, sizeof(int));
    for (j=0; j<n_dims; j++)
        neighbourhood_dims[j] = 3;
    
    int *temp = (int *) R_alloc(n_dims, sizeof(int));
    
    int *kernel_matrix_locs = (int *) R_alloc(kernel_len*n_dims, sizeof(int));
    for (i=0; i<kernel_len; i++)
    {
        vector_to_matrix_loc(i, kernel_dims, n_dims, temp);
        for (j=0; j<n_dims; j++)
            kernel_matrix_locs[i + (j*kernel_len)] = temp[j] - kernel_centre;
    }
    
    int *neighbourhood_matrix_locs = (int *) R_alloc(neighbourhood_len*n_dims, sizeof(int));
    for (i=0; i<neighbourhood_len; i++)
    {
        vector_to_matrix_loc(i, neighbourhood_dims, n_dims, temp);
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
        if (is_compatible_value(x, i, value, value_not, integer_x) && is_compatible_neighbourhood(x, x_dims, n_dims, neighbourhood_len, neighbourhood_matrix_locs, i, n_neighbours, n_neighbours_not, integer_x))
        {
            vector_to_matrix_loc((size_t) i, x_dims, n_dims, loc);
            apply_kernel(x, y, x_dims, n_dims, loc, kernel, kernel_len, kernel_matrix_locs, integer_x, INTEGER(is_brush)[0], INTEGER(is_eraser)[0]);
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

int is_compatible_neighbourhood (const SEXP x, const int *x_dims, const int n_dims, const int neighbourhood_len, const int *neighbourhood_matrix_locs, const R_len_t index, const SEXP include_list, const SEXP exclude_list, const int is_integer)
{
    if (LENGTH(include_list) == 0 && LENGTH(exclude_list) == 0)
        return TRUE;
    
    size_t i, l;
    int *x_loc = (int *) R_alloc(n_dims, sizeof(int));
    int *neighbourhood_loc = (int *) R_alloc(n_dims, sizeof(int));
    int j, n_neighbours = 0;
    
    vector_to_matrix_loc((size_t) index, x_dims, n_dims, x_loc);
    
    for (i=0; i<neighbourhood_len; i++)
    {
        for (j=0; j<n_dims; j++)
            neighbourhood_loc[j] = x_loc[j] + neighbourhood_matrix_locs[i + (j*neighbourhood_len)];
        
        if (loc_in_bounds(neighbourhood_loc, x_dims, n_dims))
        {
            matrix_to_vector_loc(neighbourhood_loc, x_dims, n_dims, &l);
        
            if (is_integer && INTEGER(x)[l] != 0)
                n_neighbours++;
            else if (!is_integer && REAL(x)[l] != 0.0)
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

void apply_kernel (const SEXP x, SEXP y, const int *x_dims, const int n_dims, const int *x_loc, const SEXP kernel, const R_len_t kernel_len, const int *kernel_matrix_locs, const int is_integer, const int is_brush, const int is_eraser)
{
    size_t i, l;
    int *current_loc = (int *) R_alloc(n_dims, sizeof(int));
    int j;
    
    double value = 0;
    
    for (i=0; i<kernel_len; i++)
    {
        for (j=0; j<n_dims; j++)
            current_loc[j] = x_loc[j] + kernel_matrix_locs[i + (j*kernel_len)];
        
        if (loc_in_bounds(current_loc, x_dims, n_dims))
        {
            matrix_to_vector_loc(current_loc, x_dims, n_dims, &l);
            
            if (is_brush)
            {
                if (is_integer)
                {
                    if (INTEGER(kernel)[i] == NA_INTEGER)
                        INTEGER(y)[l] = INTEGER(x)[l];
                    else if (is_eraser)
                        INTEGER(y)[l] = INTEGER(kernel)[i] == 0 ? INTEGER(x)[l] : 0;
                    else
                        INTEGER(y)[l] = INTEGER(kernel)[i];
                }
                else
                {
                    if (ISNA(REAL(kernel)[i]))
                        REAL(y)[l] = REAL(x)[l];
                    else if (is_eraser)
                        REAL(y)[l] = REAL(kernel)[i] == 0.0 ? REAL(x)[l] : 0.0;
                    else
                        REAL(y)[l] = REAL(kernel)[i];
                }
            }
            else
            {
                if (is_integer)
                    value += (double) INTEGER(kernel)[i] * INTEGER(x)[l];
                else
                    value += REAL(kernel)[i] * REAL(x)[l];
            }
        }
    }
    
    if (!is_brush)
    {
        matrix_to_vector_loc(x_loc, x_dims, n_dims, &l);
        if (is_integer)
            INTEGER(y)[l] = (int) round(value);
        else
            REAL(y)[l] = value;
    }
}
