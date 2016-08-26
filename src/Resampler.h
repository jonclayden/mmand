#ifndef _RESAMPLER_H_
#define _RESAMPLER_H_

#include "Array.h"
#include "Kernel.h"

typedef std::vector<double> dbl_vector;
typedef std::vector<int>    int_vector;

template <class IteratorType>
class Interpolant
{
private:
    IteratorType start, end;
    size_t len;
    double prestart, postend;
    
public:
    Interpolant (IteratorType start, IteratorType end)
        : start(start), end(end)
    {
        if (end - start < 0)
            throw std::runtime_error("End point must be after start point");
        len = size_t(end - start);
        prestart = 2*(*start) - *(start+1);
        postend = 2*(*(end-1)) - *(end-2);
    }
    
    const double operator[] (size_t i) const
    {
        if (i == -1)
            return prestart;
        else if (i == len)
            return postend;
        else if (i < -1 || i > len)
            return 0.0;
        else
            return *(start + ptrdiff_t(i));
    }
    
    const size_t & length () const { return len; }
};

class Resampler
{
protected:
    const Array<double> *original;
    Array<double> *working;
    
    Kernel *kernel;
    int kernelWidth;
    double a, b, c;
    bool toPresharpen;
    
    dbl_vector samples;
    
    template <class InputIterator, class OutputIterator>
    void presharpen (InputIterator begin, InputIterator end, OutputIterator result);
    
    template <class InputIterator, class OutputIterator>
    void interpolate (Interpolant<InputIterator> data, const std::vector<double> &locs, OutputIterator result);
    
    template <class InputIterator, class OutputIterator>
    void interpolate (InputIterator begin, InputIterator end, const std::vector<double> &locs, OutputIterator result);
    
    template <class InputIterator>
    double interpolate (InputIterator begin, InputIterator end, const double &loc);
    
    double samplePoint (const std::vector<int> &base, const std::vector<double> &offset, const int dim);
    
public:
    Resampler ()
        : original(NULL), working(NULL), kernel(NULL) {}
    
    Resampler (const Array<double> *original, Kernel * const kernel)
        : original(original), working(NULL), kernel(kernel)
    {
        kernelWidth = static_cast<int>(floor(2.0 * kernel->getSupportMax()));
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
    
    const std::vector<double> & run (const Eigen::MatrixXd &locations);
    
    const std::vector<double> & run (const std::vector<dbl_vector> &locations);
};

#endif
