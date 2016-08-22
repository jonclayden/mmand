#include <RcppEigen.h>

#include "Resampler.h"

template <class InputIterator, class OutputIterator>
void Resampler::interpolate (InputIterator begin, InputIterator end, const std::vector<double> &locs, OutputIterator result)
{
    const ptrdiff_t len = end - begin;
    std::vector<double> data(len);
    
    if (presharpen)
    {
        InputIterator input = begin;
        std::vector<double> coefs(len);
        for (int i=0; ; i++)
        {
            coefs[i] = (c / (b - (i==0 ? 0 : a*coefs[i-1])));
            data[i] = (*input - (i==0 ? 0 : a*data[i-1])) / (b - (i==0 ? 0 : a*coefs[i-1]));
            if (++input == end)
                break;
        }
        
        const int len = data.size();
        for (int i=(len-1); i>=0; i--)
            data[i] -= (i==(len-1) ? 0 : coefs[i]*data[i+1]);
    }
    else
        std::copy(begin, end, data.begin());
    
    for (int j=0; j<locs.size(); j++, ++result)
    {
        // NB: This will only work for kernels of width up to 4
        int base = static_cast<int>(floor(locs[j])) - (kernelWidth>2 ? 1 : 0);
        double value = 0.0;
        for (int k=base; k<base+kernelWidth; k++)
        {
            if (k >= 0 && k < len)
                value += data[k] * kernel->evaluate(static_cast<double>(k) - locs[j]);
        }
        
        *result = value;
    }
}

// double Resampler::samplePoint (const std::vector<int> &base, const Eigen::VectorXd &offset, const int dim)
// {
//     if (dim == 0)
//     {
//         Array<double>::ConstIterator it = original->beginLine(base, 0);
//         return interpolate(it, it+kernelWidth, offset[0]);
//     }
//     else
//     {
//         std::vector<double> elements(kernelWidth);
//         for (int i=0; i<kernelWidth; i++)
//         {
//             std::vector<int> temp = base;
//             temp[dim] += i;
//             elements[i] = samplePoint(temp, offset, dim-1);
//         }
//         return interpolate(elements.begin(), elements.end(), offset[dim]);
//     }
// }
//
// std::vector<double> & Resampler::run (const Eigen::MatrixXd &locations)
// {
//     const int nDims = original->getDimensionality();
//     samples.resize(locations.rows());
//
//     for (size_t i=0; i<samples.size(); i++)
//     {
//         std::vector<int> base(nDims);
//         Eigen::VectorXd offset(nDims);
//         for (int j=0; j<nDims; j++)
//         {
//             // NB: This will only work for kernels of width up to 4
//             base[j] = static_cast<int>(floor(locations(i,j))) - (kernelWidth>2 ? 1 : 0);
//             offset[j] = locations(i,j) - floor(locations(i,j)) + (kernelWidth>2 ? 1.0 : 0.0);
//         }
//         samples[i] = samplePoint(base, offset, nDims-1);
//     }
//
//     return samples;
// }

const std::vector<double> & Resampler::run (const std::vector<dbl_vector> &locations)
{
    const int nDims = locations.size();
    int_vector dims = original->getDimensions();
    
    working = const_cast<Array<double>*>(original);
    for (int i=0; i<nDims; i++)
    {
        dims[i] = locations[i].size();
        Array<double> *result = new Array<double>(dims, NA_REAL);
        for (size_t j=0; j<working->countLines(i); j++)
            interpolate(working->beginLine(j,i), working->endLine(j,i), locations[i], result->beginLine(j,i));
        
        if (working != original)
            delete working;
        working = result;
    }
    
    return working->getData();
}
