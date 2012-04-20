#include <R.h>
#include <Rdefines.h>
// #include <Rinternals.h>
#include <Rmath.h>

#include "morph.h"
#include "index.h"

SEXP morph_R (SEXP x, SEXP kernel, SEXP value, SEXP value_not, SEXP n_neighbours, SEXP n_neighbours_not, SEXP is_brush)
{
    R_len_t i;
    int j, n;
    SEXP y;
    
    SEXP dims = getAttrib(x, R_DimSymbol);
    int n_dims = (int) LENGTH(dims);
    int kernel_width = INTEGER(getAttrib(kernel, R_DimSymbol))[0];
    int integer_x = (int) IS_INTEGER(x);
    R_len_t len = LENGTH(x);
    
    int *x_dims = (int *) R_alloc(n_dims, sizeof(int));
    for (j=0; j<n_dims; j++)
        x_dims[j] = INTEGER(dims)[j];
    
    int *kernel_dims = (int *) R_alloc(n_dims, sizeof(int));
    for (j=0; j<n_dims; j++)
        kernel_dims[j] = kernel_width;
    
    int *neighbourhood_dims = (int *) R_alloc(n_dims, sizeof(int));
    for (j=0; j<n_dims; j++)
        neighbourhood_dims[j] = 3;
    
    int *loc = (int *) R_alloc(n_dims, sizeof(int));
    
    if (integer_x)
        PROTECT(y = NEW_INTEGER(len));
    else
        PROTECT(y = NEW_NUMERIC(len));
    
    for (i=0; i<len; i++)
    {
        if (is_compatible_value(x, i, value, value_not, integer_x) && is_compatible_neighbourhood(x, x_dims, n_dims, neighbourhood_dims, i, n_neighbours, n_neighbours_not, integer_x))
        {
            vector_to_matrix_loc((size_t) i, x_dims, n_dims, loc);
            apply_kernel(x, y, x_dims, n_dims, loc, kernel, kernel_dims, integer_x, INTEGER(is_brush)[0]);
        }
    }
    
    UNPROTECT(1);
    return y;
}

int is_compatible_value (const SEXP x, const R_len_t index, const SEXP include_list, const SEXP exclude_list, const int is_integer)
{
    R_len_t i;
    
    if (!isNull(include_list))
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
    
    if (!isNull(exclude_list))
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

int is_compatible_neighbourhood (const SEXP x, const int *x_dims, const int n_dims, const int *neighbourhood_dims, const R_len_t index, const SEXP include_list, const SEXP exclude_list, const int is_integer)
{
    if (isNull(include_list) && isNull(exclude_list))
        return TRUE;
    
    size_t i, l, neighbourhood_length = (size_t) R_pow_di(3.0, n_dims);
    int *x_loc = (int *) R_alloc(n_dims, sizeof(int));
    int *neighbourhood_loc = (int *) R_alloc(n_dims, sizeof(int));
    int j, n_neighbours = 0;
    
    vector_to_matrix_loc((size_t) index, x_dims, n_dims, x_loc);
    
    for (i=0; i<neighbourhood_length; i++)
    {
        vector_to_matrix_loc(i, neighbourhood_dims, n_dims, neighbourhood_loc);
        for (j=0; j<n_dims; j++)
            neighbourhood_loc[j] += x_loc[j] - 1;
        
        if (loc_in_bounds(neighbourhood_loc, x_dims, n_dims))
        {
            matrix_to_vector_loc(neighbourhood_loc, x_dims, n_dims, &l);
        
            if (is_integer && INTEGER(x)[l] != 0)
                n_neighbours++;
            else if (!is_integer && REAL(x)[l] != 0.0)
                n_neighbours++;
        }
    }
    
    if (!isNull(include_list))
    {
        for (i=0; i<LENGTH(include_list); i++)
        {
            if (n_neighbours == (int) REAL(include_list)[i])
                return TRUE;
        }
        return FALSE;
    }
    
    if (!isNull(exclude_list))
    {
        for (i=0; i<LENGTH(exclude_list); i++)
        {
            if (n_neighbours == (int) REAL(exclude_list)[i])
                return FALSE;
        }
        return TRUE;
    }
    
    return TRUE;
}

void apply_kernel (const SEXP x, SEXP y, const int *x_dims, const int n_dims, const int *x_loc, SEXP kernel, const int *kernel_dims, const int is_integer, const int is_brush)
{
    int kernel_centre = (kernel_dims[0] - 1) / 2;
    size_t i, l, kernel_length = n_dims * kernel_dims[0];
    int *current_loc = (int *) R_alloc(n_dims, sizeof(int));
    int j;
    
    for (i=0; i<kernel_length; i++)
    {
        vector_to_matrix_loc(i, kernel_dims, n_dims, current_loc);
        for (j=0; j<n_dims; j++)
            current_loc[j] += x_loc[j] - kernel_centre;
        
        if (loc_in_bounds(current_loc, x_dims, n_dims))
        {
            matrix_to_vector_loc(current_loc, x_dims, n_dims, &l);
            
            if (is_brush)
            {
                if (is_integer)
                    INTEGER(y)[l] = INTEGER(kernel)[i];
                else
                    REAL(y)[l] = REAL(kernel)[i];
            }
            else
            {
                if (is_integer)
                    INTEGER(y)[l] = INTEGER(kernel)[i] * INTEGER(x)[l];
                else
                    REAL(y)[l] = REAL(kernel)[i] * REAL(x)[l];
            }
        }
        
    }
}
