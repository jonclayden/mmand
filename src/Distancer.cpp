#include <Rcpp.h>

#include "Parallel.h"
#include "Distancer.h"


double initialTransform (const double &x) { return x == 0.0 ? R_PosInf : 0.0; }

inline double intersectionPoint (const std::vector<double> &vec, const int &loc, const int &vertex, const double &sqPixdim)
{
    // This is the solution (for x) to the equation
    //   y_l + p^2 * (x_l - x)^2 = y_v + p^2 * (x_v - x)^2,
    // where the first argument provides the mapping from x to y
    return (vec[loc] - vec[vertex] + sqPixdim * (loc*loc - vertex*vertex)) / (2 * sqPixdim * (loc - vertex));
}

Array<double> * Distancer::run ()
{
    // Transform the source array so that distances are zero within the region
    // and infinite elsewhere
    Array<double> *result = new Array<double>(*original);
    std::transform(original->begin(), original->end(), result->begin(), initialTransform);
    
    const std::vector<int> &dims = original->getDimensions();
    const int nDims = original->getDimensionality();
    const std::vector<double> &pixdims = original->getPixelDimensions();
    
    // This form is separable, so we apply it in one direction at a time
    for (int i=0; i<nDims; i++)
    {
        const double sqPixdim = usePixdim ? pixdims[i] * pixdims[i] : 1.0;
        
        // Lines are independent, so can be processed in parallel
        PARALLEL_LOOP_START(j, result->countLines(i))
            // The vertices are the minima of a series of parabolas. The
            // intersections are the locations where they cross
            std::vector<int> vertices;
            std::vector<double> intersections;
            
            // At least two parabolas are needed for a "real" intersection to
            // occur. Parabolas k-1 and k intersect at intersections[k]. The
            // first value here is an "off the left end" extreme
            intersections.push_back(R_NegInf);
            
            // Get an iterator to access the line of data, and take a copy
            Array<double>::Iterator it = result->beginLine(j,i);
            const std::vector<double> line(result->beginLine(j,i), result->endLine(j,i));
            for (int l=0; l<dims[i]; l++)
            {
                // Don't place a parabola if the transformed data is infinite
                if (!R_FINITE(line[l]))
                    continue;
                
                // If at least one other parabola has been placed, find the
                // relevant intersection with the new parabola
                if (!vertices.empty())
                {
                    // If the intersection with the most recently placed
                    // parabola is to the "left" of its intersection with its
                    // predecessor, the new one replaces the previous one
                    // (and so on, back through the chain)
                    double s = intersectionPoint(line, l, vertices.back(), sqPixdim);
                    while (s <= intersections.back())
                    {
                        vertices.pop_back();
                        intersections.pop_back();
                        s = intersectionPoint(line, l, vertices.back(), sqPixdim);
                    }
                    intersections.push_back(s);
                }
                
                // Place the new parabola, centred at l
                vertices.push_back(l);
            }
            
            // Add an "off the right end" extreme value for use below
            intersections.push_back(R_PosInf);
            
            // If no parabolas have been placed, there's nothing else to do
            if (vertices.empty())
                continue;
            
            // Step back over the data, replacing each element with the value
            // of the lowest parabola at that location
            for (int k=0, l=0; l<dims[i]; l++)
            {
                // The relevant parabola is the last one whose intersection
                // point we haven't yet passed
                const double q = static_cast<double>(l);
                while (intersections[k+1] < q)
                    k++;
                double dx = q - vertices[k];
                if (usePixdim)
                    dx *= pixdims[i];
                it[l] = line[vertices[k]] + dx * dx;
            }
        PARALLEL_LOOP_END
    }
    
    // Take the square-root of each value to get Euclidean distance. The
    // function pointer cast is needed to resolve the overload on sqrt()
    std::transform(result->begin(), result->end(), result->begin(), (double(*)(double)) ::sqrt);
    return result;
}
