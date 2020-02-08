#include <Rcpp.h>

#include "Distancer.h"

double initialTransform (const double &x) { return x == 0.0 ? R_PosInf : 0.0; }

inline double intersectionPoint (Array<double>::Iterator &it, const int &loc, const int &vertex)
{
    if (!R_FINITE(it[loc]) || !R_FINITE(it[vertex]))
        return R_PosInf;
    else
    {
        // This is the solution (for x) to the equation
        //   y_l + (x_l - x)^2 = y_v + (x_v - x)^2,
        // where the iterator provides the mapping from x to y
        return ((it[loc] + loc*loc) - (it[vertex] + vertex*vertex)) / (2 * (loc - vertex));
    }
}

Array<double> * Distancer::run ()
{
    // Transform the source array so that distances are zero within the region and infinite elsewhere
    Array<double> *result = new Array<double>(*original);
    std::transform(original->begin(), original->end(), result->begin(), initialTransform);
    
    const std::vector<int> &dims = original->getDimensions();
    int nDims = original->getDimensionality();
    
    for (int i=0; i<nDims; i++)
    {
        for (size_t j=0; j<result->countLines(i); j++)
        {
            std::vector<int> vertices(dims[i]);
            std::vector<double> intersections(dims[i]);
            
            vertices[0] = 0;
            intersections[0] = R_NegInf;
            intersections[1] = R_PosInf;
            
            int k = 0;
            Array<double>::Iterator it = result->beginLine(j,i);
            for (int l=1; l<dims[i]; l++)
            {
                double s = intersectionPoint(it, l, vertices[k]);
                while (s <= intersections[k])
                    s = intersectionPoint(it, l, vertices[--k]);
                vertices[++k] = l;
                intersections[k] = s;
                intersections[k+1] = R_PosInf;
            }
            
            k = 0;
            for (int l=0; l<dims[i]; l++)
            {
                double q = static_cast<double>(l);
                while (intersections[k+1] < q)
                    k++;
                double dx = q - vertices[k];
                it[l] = it[vertices[k]] + dx * dx;
            }
        }
    }
    
    return result;
}
