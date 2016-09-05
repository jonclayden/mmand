#ifndef _RESAMPLER_H_
#define _RESAMPLER_H_

#include "Array.h"
#include "Kernel.h"

typedef std::vector<double> dbl_vector;
typedef std::vector<int>    int_vector;

// Base class representing a line of presharpened data
// Has a value at index -1 and len; zero beyond these extended limits
class Interpolant
{
protected:
    size_t len;
    double prestart, postend;
    
public:
    virtual double operator() (ptrdiff_t i) const { return 0.0; }
    
    const size_t & length () const { return len; }
};

// Cached version, which stores the target data for repeated access
class CachedInterpolant : public Interpolant
{
private:
    dbl_vector data;
    
public:
    template <class IteratorType>
    CachedInterpolant (IteratorType start, IteratorType end)
        : data(start,end)
    {
        len = data.size();
        prestart = 2*data[0] - data[1];
        postend = 2*data[len-1] - data[len-2];
    }
    
    double operator() (ptrdiff_t i) const
    {
        if (i > -1 && i < len)
            return data[i];
        else if (i == -1)
            return prestart;
        else if (i == len)
            return postend;
        else
            return 0.0;
    }
};

// Uncached version, which is faster for single-pass access
template <class IteratorType>
class UncachedInterpolant : public Interpolant
{
private:
    IteratorType start, end;
    
public:
    UncachedInterpolant (IteratorType start, IteratorType end)
        : start(start), end(end)
    {
        len = end - start;
        prestart = 2*(*start) - (*(++start));
        --end;
        postend = 2*(*end) - (*(--end));
    }
    
    double operator() (ptrdiff_t i) const
    {
        if (i > -1 && i < len)
            return *(start + i);
        else if (i == -1)
            return prestart;
        else if (i == len)
            return postend;
        else
            return 0.0;
    }
};

// Main class responsible for resampling
class Resampler
{
protected:
    const Array<double> *original;
    Array<double> *working;
    
    Kernel *kernel;
    int kernelWidth, baseOffset;
    double a, b, c;
    bool toPresharpen;
    
    dbl_vector samples;
    
    template <class InputIterator, class OutputIterator>
    void presharpen (InputIterator begin, InputIterator end, OutputIterator result);
    
    void presharpen ();
    
    // For some reason, the (virtual) call operator doesn't function as expected
    // if these methods take an Interpolant rather than a derived class
    template <class InputIterator>
    double interpolate (const UncachedInterpolant<InputIterator> data, const double &loc);
    
    template <class OutputIterator>
    void interpolate (const CachedInterpolant data, const std::vector<double> &locs, OutputIterator result);
    
    double samplePoint (const std::vector<int> &base, const std::vector<double> &offset, const int dim);
    
public:
    Resampler ()
        : original(NULL), working(NULL), kernel(NULL) {}
    
    Resampler (const Array<double> *original, Kernel * const kernel)
        : original(original), working(NULL), kernel(kernel)
    {
        kernelWidth = static_cast<int>(floor(2.0 * kernel->getSupportMax()));
        if (kernelWidth > 4)
            throw std::runtime_error("Kernel widths of greater than 4 are not supported");
        baseOffset = std::max(0, kernelWidth/2 - 1);
        
        a = kernel->evaluate(-1.0);
        b = kernel->evaluate(0.0);
        c = kernel->evaluate(1.0);
        toPresharpen = (a != 0.0 || b != 1.0 || c != 0.0);
    }
    
    ~Resampler ()
    {
        delete original;
        delete working;
        delete kernel;
    }
    
    const std::vector<double> & run (const Rcpp::NumericMatrix &locations);
    
    const std::vector<double> & run (const std::vector<dbl_vector> &locations);
};

#endif
