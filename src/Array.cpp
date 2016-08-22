#include <RcppEigen.h>

#include "Array.h"

template <typename DataType>
Neighbourhood Array<DataType>::getNeighbourhood () const
{
    return this->getNeighbourhood(dims);
}

template <typename DataType>
Neighbourhood Array<DataType>::getNeighbourhood (const int width) const
{
    std::vector<int> widths(nDims, width);
    return this->getNeighbourhood(widths);
}

template <typename DataType>
Neighbourhood Array<DataType>::getNeighbourhood (const std::vector<int> &widths) const
{
    Neighbourhood neighbourhood;
    
    neighbourhood.widths = widths;
    std::vector<int> extremes(nDims);
    for (int i=0; i<nDims; i++)
    {
        if (neighbourhood.widths[i] % 2 == 0)
            neighbourhood.widths[i]++;
        extremes[i] = (neighbourhood.widths[i] - 1) / 2;
    }
    
    neighbourhood.size = 1;
    std::vector<size_t> steps(nDims+1);
    steps[0] = 1;
    for (int i=0; i<nDims; i++)
    {
        neighbourhood.size *= neighbourhood.widths[i];
        steps[i+1] = steps[i] * dims[i];
    }
    
    neighbourhood.locs.resize(neighbourhood.size, nDims);
    neighbourhood.offsets.resize(neighbourhood.size);
    
    for (int j=0; j<neighbourhood.size; j++)
    {
        if (j==0)
        {
            for (int i=0; i<nDims; i++)
                neighbourhood.locs(j,i) = -extremes[i];
        }
        else
        {
            neighbourhood.locs(j,0) = neighbourhood.locs(j-1,0) + 1;
            for (int i=0; i<nDims; i++)
            {
                if (neighbourhood.locs(j,i) > extremes[i])
                {
                    neighbourhood.locs(j,i) = -extremes[i];
                    neighbourhood.locs(j,i+1) = neighbourhood.locs(j-1,i+1) + 1;
                }
                else if (i < (nDims-1))
                    neighbourhood.locs(j,i+1) = neighbourhood.locs(j-1,i+1);
            }
        }
        
        neighbourhood.offsets[j] = 0;
        for (int i=0; i<nDims; i++)
            neighbourhood.offsets[j] += neighbourhood.locs(j,i) * steps[i];
    }
    
    return neighbourhood;
}

template <typename DataType>
void Array<DataType>::flattenIndex (const std::vector<int> &loc, size_t &result) const
{
    // Dimensionalities 1-4 are most common so treat them as special cases for speed
    switch (nDims)
    {
        case 1:
        result = loc[0];
        break;
        
        case 2:
        result = loc[0] + dims[0] * loc[1];
        break;
        
        case 3:
        result = loc[0] + dims[0] * (loc[1] + dims[1] * loc[2]);
        break;
        
        case 4:
        result = loc[0] + dims[0] * (loc[1] + dims[1] * (loc[2] + dims[2] * loc[3]));
        
        default:
        {
            result = loc[0];
            for (int i=1; i<nDims; i++)
                result += loc[i] * strides[i];
        }
    }
}

template <typename DataType>
void Array<DataType>::expandIndex (const size_t &loc, std::vector<int> &result) const
{
    result[0] = loc % dims[0];
    for (int i=1; i<nDims; i++)
        result[i] = (loc / strides[i]) % dims[i];
}

// Tell the compiler that we're going to need these specialisations (otherwise
// it won't generate the relevant code and we'll get a linker error)
template class Array<double>;
