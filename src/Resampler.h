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
    int kernelWidth;
    double a, b, c;
    bool presharpen;
    
    dbl_vector samples;
    
    template <class IteratorType>
    double interpolate (IteratorType begin, IteratorType end, const double loc);
    
    double samplePoint (const std::vector<int> &base, const Eigen::VectorXd &offset, const int dim);
    
public:
    Resampler () {}
    
    Resampler (const Array<double> *original, Kernel * const kernel)
        : original(original), kernel(kernel)
    {
        kernelWidth = static_cast<int>(floor(2.0 * kernel->getSupportMax()));
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
    
    std::vector<double> & run (const Eigen::MatrixXd &locations);
};

#endif
