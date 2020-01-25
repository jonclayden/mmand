#include <Rcpp.h>

#include "Distancer.h"

const std::vector<double> & Distancer::run ()
{
    const std::vector<int> &dims = original->getDimensions();
    int nDims = original->getDimensionality();
    const size_t nSamples = original->size();
    
    for (int i=0; i<nDims; i++)
    {
        for (size_t j=0; j<original->countLines(i); j++)
        {
            std::vector<int> vertices;
            std::vector<double> intersections;
            
            vertices.push_back(0);
            intersections.push_back(R_NegInf);
            intersections.push_back(R_PosInf);
            
            int k = 0;
            Array<double>::Iterator it = original->beginLine(j,i);
            for (int l=0; l<dims[i]; l++)
            {
                intersections.push_back(R_PosInf);
                double q = static_cast<double>(l);
                int &p = vertices[k];
                double s = ((it[l] + l*l) - (it[p] + p*p)) / (2 * (l - p));
                while (s <= intersections[k])
                {
                    p = vertices[--k];
                    s = ((it[l] + l*l) - (it[p] + p*p)) / (2 * (l - p));
                }
                vertices[++k] = l;
                intersections[k] = s;
                intersections.push_back(R_PosInf);
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
    
    return original->getData();
}
