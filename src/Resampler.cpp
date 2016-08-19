#include <RcppEigen.h>

#include "Resampler.h"

template <class IteratorType>
double Resampler::interpolate (IteratorType begin, IteratorType end, const double loc)
{
    std::vector<double> data;
    IteratorType it;
    
    if (presharpen)
    {
        int i;
        std::vector<double> coefs;
        for (i=0, it=begin; ; i++)
        {
            coefs.push_back(c / (b - (i==0 ? 0 : a*coefs[i-1])));
            data.push_back((*it - (i==0 ? 0 : a*data[i-1])) / (b - (i==0 ? 0 : a*coefs[i-1])));
            if (++it == end)
                break;
        }
        
        const int len = data.size();
        for (i=(len-1); i>=0; i--)
            data[i] -= (i==(len-1) ? 0 : coefs[i]*data[i+1]);
    }
    else
    {
        data.resize(end - begin);
        std::copy(begin, end, data.begin());
    }
    
    double result = 0.0;
    for (int i=0; i<data.size(); i++)
        result += data[i] * kernel->evaluate(static_cast<double>(i) - loc);
    
    return result;
}

double Resampler::samplePoint (const std::vector<int> &base, const Eigen::VectorXd &offset, const int dim)
{
    if (dim == 0)
    {
        Array<double>::ConstIterator it = original->beginLine(base, 0);
        return interpolate(it, it+kernelWidth, offset[0]);
    }
    else
    {
        std::vector<double> elements(kernelWidth);
        for (int i=0; i<kernelWidth; i++)
        {
            std::vector<int> temp = base;
            temp[dim] += i;
            elements[i] = samplePoint(temp, offset, dim-1);
        }
        return interpolate(elements.begin(), elements.end(), offset[dim]);
    }
}

std::vector<double> & Resampler::run (const Eigen::MatrixXd &locations)
{
    const int nDims = original->getDimensionality();
    samples.resize(locations.rows());
    
    for (size_t i=0; i<samples.size(); i++)
    {
        std::vector<int> base(nDims);
        Eigen::VectorXd offset(nDims);
        for (int j=0; j<nDims; j++)
        {
            // NB: This will only work for kernels of width up to 4
            base[j] = static_cast<int>(floor(locations(i,j))) - (kernelWidth>2 ? 1 : 0);
            offset[j] = locations(i,j) - floor(locations(i,j)) + (kernelWidth>2 ? 1.0 : 0.0);
        }
        samples[i] = samplePoint(base, offset, nDims-1);
    }
    
    return samples;
}
