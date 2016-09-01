#include <RcppEigen.h>

#include "Resampler.h"

template <class InputIterator, class OutputIterator>
void Resampler::presharpen (InputIterator begin, InputIterator end, OutputIterator result)
{
    const ptrdiff_t len = end - begin;
    std::vector<double> coefs(len, 0.0);
    
    *result = *begin;
    for (int i=1; i<(len-1); i++)
    {
        coefs[i] = c / (b - a*coefs[i-1]);
        const double temp = a * (*result);
        *(++result) = (*(++begin) - temp) / (b - a*coefs[i-1]);
    }
    *(++result) = *(++begin);
    
    for (int i=(len-1); i>0; i--)
    {
        const double temp = coefs[i-1] * (*result);
        *(--result) -= temp;
    }
}

void Resampler::presharpen ()
{
    delete working;
    working = new Array<double>(*original);
    
    if (toPresharpen)
    {
        for (int i=0; i<working->getDimensionality(); i++)
        {
            for (size_t j=0; j<working->countLines(i); j++)
                presharpen(working->beginLine(j,i), working->endLine(j,i), working->beginLine(j,i));
        }
    }
}

template <class OutputIterator>
void Resampler::interpolate (Interpolant data, const std::vector<double> &locs, OutputIterator result)
{
    const int baseOffset = std::max(0, kernelWidth/2 - 1);
    
    for (int j=0; j<locs.size(); j++, ++result)
    {
        const int base = static_cast<int>(floor(locs[j])) - baseOffset;
        double value = 0.0;
        double kernelSum = 0.0;
        for (int k=base; k<base+kernelWidth; k++)
        {
            const double kernelValue = kernel->evaluate(static_cast<double>(k) - locs[j]);
            if (kernelValue != 0.0)
            {
                value += data(k) * kernelValue;
                kernelSum += kernelValue;
            }
        }
        
        // if (kernelSum != 1.0 && kernelSum != 0.0)
        //     value /= kernelSum;;
        *result = value;
    }
}

template <class InputIterator, class OutputIterator>
void Resampler::interpolate (InputIterator begin, InputIterator end, const std::vector<double> &locs, OutputIterator result)
{
    const ptrdiff_t len = end - begin;
    std::vector<double> data(begin, end);
    
    for (int j=0; j<locs.size(); j++, ++result)
    {
        // NB: This will only work for kernels of width up to 4
        const int base = static_cast<int>(floor(locs[j])) - (kernelWidth>2 ? 1 : 0);
        double value = 0.0;
        double kernelSum = 0.0;
        for (int k=base; k<base+kernelWidth; k++)
        {
            if (k >= 0 && k < len)
            {
                const double kernelValue = kernel->evaluate(static_cast<double>(k) - locs[j]);
                if (kernelValue != 0.0)
                {
                    value += data[k] * kernelValue;
                    kernelSum += kernelValue;
                }
            }
        }
        
        // if (kernelSum != 1.0 && kernelSum != 0.0)
        //     value /= kernelSum;;
        *result = value;
    }
}

template <class InputIterator>
double Resampler::interpolate (InputIterator begin, InputIterator end, const double &loc)
{
    const ptrdiff_t len = end - begin;
    std::vector<double> data(begin, end);
    
    double value = 0.0;
    double kernelSum = 0.0;
    for (int k=0; k<kernelWidth; k++)
    {
        if (k >= 0 && k < len)
        {
            const double kernelValue = kernel->evaluate(static_cast<double>(k) - loc);
            if (kernelValue != 0.0)
            {
                value += data[k] * kernelValue;
                kernelSum += kernelValue;
            }
        }
    }
    
    // if (kernelSum != 1.0 && kernelSum != 0.0)
    //     value /= kernelSum;
    return value;
}

double Resampler::samplePoint (const std::vector<int> &base, const std::vector<double> &offset, const int dim)
{
    if (dim == 0)
    {
        Array<double>::Iterator start = working->beginLine(base, 0);
        Array<double>::Iterator end = working->endLine(base, 0);
        if (end > start+kernelWidth)
            return interpolate(start, start+kernelWidth, offset[0]);
        else
            return interpolate(start, end, offset[0]);
    }
    else
    {
        std::vector<double> elements;
        const std::vector<int> &dims = working->getDimensions();
        for (int i=0; i<kernelWidth; i++)
        {
            std::vector<int> temp = base;
            temp[dim] += i;
            if (temp[dim] < dims[dim])
                elements.push_back(samplePoint(temp, offset, dim-1));
        }
        return interpolate(elements.begin(), elements.end(), offset[dim]);
    }
}

const std::vector<double> & Resampler::run (const Eigen::MatrixXd &locations)
{
    const int_vector &dims = original->getDimensions();
    const int nDims = locations.cols();
    const int nSamples = locations.rows();
    working = new Array<double>(*original);
    
    if (toPresharpen)
    {
        for (int i=0; i<nDims; i++)
        {
            for (size_t j=0; j<working->countLines(i); j++)
                presharpen(working->beginLine(j,i), working->endLine(j,i), working->beginLine(j,i));
        }
    }
    
    samples.resize(nSamples);
    int_vector base(nDims);
    dbl_vector offset(nDims);
    const int baseOffset = std::max(0, kernelWidth/2 - 1);
    for (size_t k=0; k<nSamples; k++)
    {
        for (int i=0; i<nDims; i++)
        {
            base[i] = static_cast<int>(floor(locations(k,i))) - baseOffset;
            offset[i] = locations(k,i) - static_cast<double>(base[i]);
            if (base[i] < 0)
            {
                offset[i] += static_cast<double>(base[i]);
                base[i] = 0;
            }
            else if (base[i] >= dims[i])
            {
                offset[i] += static_cast<double>(base[i]) - dims[i] + 1.0;
                base[i] = dims[i] - 1;
            }
        }
        samples[k] = samplePoint(base, offset, nDims-1);
    }
    
    return samples;
}

const std::vector<double> & Resampler::run (const std::vector<dbl_vector> &locations)
{
    const int nDims = locations.size();
    int_vector dims = original->getDimensions();
    
    presharpen();

    for (int i=0; i<nDims; i++)
    {
        dims[i] = locations[i].size();
        Array<double> *result = new Array<double>(dims, NA_REAL);
        for (size_t j=0; j<working->countLines(i); j++)
        {
            Interpolant interpolant(working->beginLine(j,i), working->endLine(j,i));
            interpolate(interpolant, locations[i], result->beginLine(j,i));
        }

        delete working;
        working = result;
    }
    
    return working->getData();
}
