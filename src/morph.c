#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

#include "morph.h"
#include "index.h"

SEXP morph_R (SEXP x, SEXP kernel, SEXP value, SEXP value_not, SEXP n_neighbours, SEXP n_neighbours_not, SEXP operator, SEXP merge)
{
    R_len_t i, ii;
    SEXP y;
    int j;
    double kernel_sum = 0.0;
    
    int *x_dims = INTEGER(getAttrib(x, R_DimSymbol));
    int *kernel_dims = INTEGER(getAttrib(kernel, R_DimSymbol));
    int n_dims = (int) LENGTH(getAttrib(x, R_DimSymbol));
    int integer_x = (int) IS_INTEGER(x);
    R_len_t len = LENGTH(x);
    R_len_t kernel_len = LENGTH(kernel);
    R_len_t neighbourhood_len = ((R_len_t) R_pow_di(3.0, n_dims)) - 1;
    
    int *kernel_centre = (int *) R_alloc(n_dims, sizeof(int));
    for (j=0; j<n_dims; j++)
        kernel_centre[j] = (kernel_dims[j] - 1) / 2;
    
    int *neighbourhood_dims = (int *) R_alloc(n_dims, sizeof(int));
    for (j=0; j<n_dims; j++)
        neighbourhood_dims[j] = 3;
    
    int *temp = (int *) R_alloc(2*n_dims, sizeof(int));
    double *values = (double *) R_alloc(kernel_len, sizeof(double));
    
    int_or_double_ptr x_p, y_p, kernel_p;
    if (integer_x)
    {
        x_p.i = INTEGER(x);
        kernel_p.i = INTEGER(kernel);
    }
    else
    {
        x_p.d = REAL(x);
        kernel_p.d = REAL(kernel);
    }
    
    int *kernel_matrix_locs = (int *) R_alloc(kernel_len*n_dims, sizeof(int));
    for (i=0; i<kernel_len; i++)
    {
        vector_to_matrix_loc(i, kernel_dims, n_dims, temp);
        for (j=0; j<n_dims; j++)
            kernel_matrix_locs[i + (j*kernel_len)] = temp[j] - kernel_centre[j];
        
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
        y_p.i = INTEGER(y);
        memcpy(y_p.i, x_p.i, ((size_t) len)*sizeof(int));
    }
    else
    {
        PROTECT(y = NEW_NUMERIC(len));
        y_p.d = REAL(y);
        memcpy(y_p.d, x_p.d, ((size_t) len)*sizeof(double));
    }
    
    for (i=0; i<len; i++)
    {
        if (is_compatible_value(x_p, i, value, value_not, integer_x) && is_compatible_neighbourhood(x_p, x_dims, n_dims, neighbourhood_len, neighbourhood_matrix_locs, i, n_neighbours, n_neighbours_not, integer_x, temp))
        {
            vector_to_matrix_loc((size_t) i, x_dims, n_dims, loc);
            apply_kernel(x_p, y_p, x_dims, n_dims, loc, kernel_p, kernel_len, kernel_matrix_locs, kernel_sum, integer_x, CHAR(STRING_ELT(operator,0)), CHAR(STRING_ELT(merge,0)), temp, values);
        }
    }
    
    UNPROTECT(1);
    return y;
}

int is_compatible_value (const int_or_double_ptr x_p, const R_len_t index, const SEXP include_list, const SEXP exclude_list, const int is_integer)
{
    R_len_t i;
    
    if (LENGTH(include_list) > 0)
    {
        for (i=0; i<LENGTH(include_list); i++)
        {
            if (is_integer && x_p.i[index] == (int) REAL(include_list)[i])
                return TRUE;
            else if (!is_integer && x_p.d[index] == REAL(include_list)[i])
                return TRUE;
        }
        return FALSE;
    }
    
    if (LENGTH(exclude_list) > 0)
    {
        for (i=0; i<LENGTH(exclude_list); i++)
        {
            if (is_integer && x_p.i[index] == (int) REAL(exclude_list)[i])
                return FALSE;
            else if (!is_integer && x_p.d[index] == REAL(exclude_list)[i])
                return FALSE;
        }
        return TRUE;
    }
    
    return TRUE;
}

int is_compatible_neighbourhood (const int_or_double_ptr x_p, const int *x_dims, const int n_dims, const int neighbourhood_len, const int *neighbourhood_matrix_locs, const R_len_t index, const SEXP include_list, const SEXP exclude_list, const int is_integer, int *temp)
{
    if (LENGTH(include_list) == 0 && LENGTH(exclude_list) == 0)
        return TRUE;
    
    size_t i, l;
    int j, n_neighbours = 0;
    
    vector_to_matrix_loc((size_t) index, x_dims, n_dims, temp+n_dims);
    
    for (i=0; i<neighbourhood_len; i++)
    {
        for (j=0; j<n_dims; j++)
            temp[j] = temp[j+n_dims] + neighbourhood_matrix_locs[i + (j*neighbourhood_len)];
        
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

void apply_kernel (const int_or_double_ptr x_p, int_or_double_ptr y_p, const int *x_dims, const int n_dims, const int *x_loc, const int_or_double_ptr kernel_p, const R_len_t kernel_len, const int *kernel_matrix_locs, const double kernel_sum, const int is_integer, const char *operator, const char *merge, int *temp, double *values)
{
    R_len_t i;
    int j, values_present, kernel_valid, element;
    size_t l;
    double final_value;
    
    values_present = 0;
    for (i=0; i<kernel_len; i++)
    {
        if (is_integer && !ISNA(kernel_p.i[i]))
            kernel_valid = TRUE;
        else if (!is_integer && !ISNA(kernel_p.d[i]))
            kernel_valid = TRUE;
        else
            kernel_valid = FALSE;
            
        if (kernel_valid)
        {
            for (j=0; j<n_dims; j++)
                temp[j] = x_loc[j] + kernel_matrix_locs[i + (j*kernel_len)];
        
            if (loc_in_bounds(temp, x_dims, n_dims))
            {
                matrix_to_vector_loc(temp, x_dims, n_dims, &l);
            
                if (is_integer)
                {
                    if (strcmp(operator,"+") == 0)
                        values[i] = (double) x_p.i[l] + kernel_p.i[i];
                    else if (strcmp(operator,"-") == 0)
                        values[i] = (double) x_p.i[l] - kernel_p.i[i];
                    else if (strcmp(operator,"*") == 0)
                        values[i] = (double) x_p.i[l] * kernel_p.i[i];
                    else if (strcmp(operator,"i") == 0)
                        values[i] = (double) (kernel_p.i[i] != 0 ? x_p.i[l] : NA_REAL);
                    else if (strcmp(operator,"1") == 0)
                        values[i] = (kernel_p.i[i] != 0 ? 1 : 0);
                    else if (strcmp(operator,"0") == 0)
                        values[i] = (kernel_p.i[i] != 0 ? 0 : 1);
                }
                else
                {
                    if (strcmp(operator,"+") == 0)
                        values[i] = x_p.d[l] + kernel_p.d[i];
                    else if (strcmp(operator,"-") == 0)
                        values[i] = x_p.d[l] - kernel_p.d[i];
                    else if (strcmp(operator,"*") == 0)
                        values[i] = x_p.d[l] * kernel_p.d[i];
                    else if (strcmp(operator,"i") == 0)
                        values[i] = (kernel_p.d[i] != 0.0 ? x_p.d[l] : NA_REAL);
                    else if (strcmp(operator,"1") == 0)
                        values[i] = (kernel_p.i[i] != 0.0 ? 1.0 : 0.0);
                    else if (strcmp(operator,"0") == 0)
                        values[i] = (kernel_p.i[i] != 0.0 ? 0.0 : 1.0);
                }
            }
        }
        else
            values[i] = NA_REAL;
        
        if (!ISNA(values[i]))
            values_present++;
    }
    
    if (values_present == 0)
        final_value = NA_REAL;
    else if (kernel_len == 1)
        final_value = values[0];
    else if (strcmp(merge,"sum") == 0 || strcmp(merge,"mean") == 0)
    {
        final_value = 0.0;
        for (i=0; i<kernel_len; i++)
            final_value += (ISNA(values[i]) ? 0.0 : values[i]);
        
        if (strcmp(merge,"mean") == 0)
            final_value /= (double) values_present;
    }
    else if (strcmp(merge,"min") == 0)
    {
        // Partial sort: ensure that element 0 is in the right place
        rPsort(values, (int) kernel_len, 0);
        final_value = values[0];
    }
    else if (strcmp(merge,"max") == 0)
    {
        element = values_present - 1;
        rPsort(values, (int) kernel_len, element);
        final_value = values[element];
    }
    else if (strcmp(merge,"median") == 0)
    {
        element = (int) floor((double) values_present / 2.0);
        rPsort(values, (int) kernel_len, element);
        final_value = values[element];
    }
    
    matrix_to_vector_loc(x_loc, x_dims, n_dims, &l);
    if (is_integer)
        y_p.i[l] = (ISNA(final_value) ? NA_INTEGER : ((int) round(final_value)));
    else
        y_p.d[l] = final_value;
}
