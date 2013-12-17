#ifndef _RESAMPLER_HPP_
#define _RESAMPLER_HPP_

#include "Array.hpp"
#include "Kernel.hpp"

typedef std::vector<double> dbl_vector;
typedef std::vector<int>    int_vector;

class Resampler
{
private:
    Array *original;
    Kernel *kernel;
    
    std::vector<int_vector> samplingLocations;
    std::vector<double> samples;
    
    int nDims;
    size_type nSamples;
    
public:
    Resampler () {}
    
    Resampler (Array * const original, Kernel * const kernel)
        : original(original), kernel(kernel) {}
    
    ~Resampler ()
    {
        delete original;
        delete kernel;
    }
    
    void setSamplingLocations (const std::vector<int_vector> &samplingLocations)
    {
        this->samplingLocations = samplingLocations;
        nDims = samplingLocations.size();
        
        nSamples = 1;
        for (int i=0; i<nDims; i++)
            nSamples *= samplingLocations[i].size();
        
        samples.resize(nSamples);
    }
    
    std::vector<double> & run ();
};

#endif
