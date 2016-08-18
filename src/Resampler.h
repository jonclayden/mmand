#ifndef _RESAMPLER_H_
#define _RESAMPLER_H_

#include "Array.h"
#include "Kernel.h"

typedef std::vector<double> dbl_vector;
typedef std::vector<int>    int_vector;

class Resampler
{
protected:
    const Array<double> *original;
    Kernel *kernel;
    double a, b, c;
    bool presharpen;
    
public:
    Resampler () {}
    
    Resampler (const Array<double> *original, Kernel * const kernel)
        : original(original), kernel(kernel)
    {
        a = kernel->evaluate(-1.0);
        b = kernel->evaluate(0.0);
        c = kernel->evaluate(1.0);
        presharpen = (a != 0.0 || b != 1.0 || c != 0.0);
    }
    
    ~Resampler ()
    {
        delete original;
        delete kernel;
    }
    
    double interpolate (Array<double>::ConstIterator begin, Array<double>::ConstIterator end, const double loc);
};

#endif
