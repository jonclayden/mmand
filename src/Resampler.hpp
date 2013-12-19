#ifndef _RESAMPLER_HPP_
#define _RESAMPLER_HPP_

#include "Array.hpp"
#include "Kernel.hpp"

typedef std::vector<double> dbl_vector;
typedef std::vector<int>    int_vector;

class SamplingScheme
{
public:
    SamplingScheme () {}
    
    virtual const double at (int sample, int dim) { return NA_REAL; }
    
    virtual int getNDims () { return NA_INTEGER; }
    
    virtual int getNSamples () { return NA_INTEGER; }
};

class GeneralSamplingScheme: public SamplingScheme
{
private:
    arma::mat locations;
    
public:
    GeneralSamplingScheme (const arma::mat &locations)
        : locations(locations) {}
    
    const double at (int sample, int dim) { return locations(sample,dim); }
    
    int getNDims () { return locations.n_cols; }
    
    int getNSamples () { return locations.n_rows; }
};

class GriddedSamplingScheme: public SamplingScheme
{
private:
    std::vector<dbl_vector> locations;
    std::vector<long> steps;
    int_vector dims;
    int nSamples;
    
public:
    GriddedSamplingScheme (const std::vector<dbl_vector> &locations)
        : locations(locations)
    {
        int nDims = locations.size();
        dims = int_vector(nDims);
        steps = std::vector<long>(nDims+1);
        steps[0] = 1;
        nSamples = 1;
        
        for (int i=0; i<nDims; i++)
        {
            dims[i] = locations[i].size();
            nSamples *= dims[i];
            steps[i+1] = steps[i] * dims[i];
        }
    }
    
    const double at (int sample, int dim)
    {
        int index = (sample / steps[dim]) % dims[dim];
        return locations[dim][index];
    }
    
    int getNDims () { return locations.size(); }
    
    int getNSamples () { return this->nSamples; }
};

class Resampler
{
private:
    Array *original;
    Kernel *kernel;
    SamplingScheme *sampler;
    dbl_vector samples;
    
public:
    Resampler () {}
    
    Resampler (Array * const original, Kernel * const kernel)
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
