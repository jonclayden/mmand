#ifndef _RESAMPLER_H_
#define _RESAMPLER_H_

#include "Array.h"
#include "Kernel.h"

typedef std::vector<double> dbl_vector;
typedef std::vector<int>    int_vector;

class SamplingScheme
{
public:
    SamplingScheme () {}
    
    virtual const double at (int sample, int dim) const { return NA_REAL; }
    
    virtual int getNDims () const { return NA_INTEGER; }
    
    virtual int getNSamples () const { return NA_INTEGER; }
};

class GeneralSamplingScheme: public SamplingScheme
{
private:
    Eigen::MatrixXd locations;
    
public:
    GeneralSamplingScheme (const Eigen::MatrixXd &locations)
        : locations(locations) {}
    
    const double at (int sample, int dim) const { return locations(sample,dim); }
    
    int getNDims () const { return locations.cols(); }
    
    int getNSamples () const { return locations.rows(); }
};

class GriddedSamplingScheme: public SamplingScheme
{
private:
    std::vector<dbl_vector> locations;
    std::vector<size_t> steps;
    int_vector dims;
    int nSamples;
    
public:
    GriddedSamplingScheme (const std::vector<dbl_vector> &locations)
        : locations(locations)
    {
        int nDims = locations.size();
        dims = int_vector(nDims);
        steps = std::vector<size_t>(nDims+1);
        steps[0] = 1;
        nSamples = 1;
        
        for (int i=0; i<nDims; i++)
        {
            dims[i] = locations[i].size();
            nSamples *= dims[i];
            steps[i+1] = steps[i] * dims[i];
        }
    }
    
    const double at (int sample, int dim) const
    {
        int index = (sample / steps[dim]) % dims[dim];
        return locations[dim][index];
    }
    
    int getNDims () const { return locations.size(); }
    
    int getNSamples () const { return this->nSamples; }
};

class Resampler
{
private:
    Array<double> *original;
    Kernel *kernel;
    SamplingScheme *sampler;
    dbl_vector samples;
    
public:
    Resampler () {}
    
    Resampler (Array<double> * const original, Kernel * const kernel)
        : original(original), kernel(kernel) {}
    
    ~Resampler ()
    {
        delete original;
        delete kernel;
        delete sampler;
    }
    
    void setSamplingScheme (SamplingScheme * const sampler) { this->sampler = sampler; }
    
    dbl_vector & run ();
};

#endif
