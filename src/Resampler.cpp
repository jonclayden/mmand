#include <RcppArmadillo.h>

#include "Resampler.hpp"

std::vector<double> & Resampler::run ()
{
    // Precalculate index deltas for the neighbourhood constituting the support region for the kernel
    int neighbourhoodWidth = static_cast<int>(floor(2.0*kernel->getSupportMax()));
    Neighbourhood neighbourhood = original->getNeighbourhood(neighbourhoodWidth);
    
    std::vector<int> &dims = original->getDims();
    int nDims = sampler->getNDims();
    int nSamples = sampler->getNSamples();
    samples.resize(nSamples);
    
    std::vector<int> nearestNeighbour(nDims);
    std::vector<double> nearestNeighbourOffset(nDims);
    long nearestNeighbourIndex;
    
    for (int i=0; i<nSamples; i++)
    {
        for (int j=0; j<nDims; j++)
        {
            nearestNeighbour[j] = static_cast<int>(round(sampler->at(i,j)));
            nearestNeighbourOffset[j] = static_cast<double>(nearestNeighbour[j]) - sampler->at(i,j);
        }
        original->flattenIndex(nearestNeighbour, nearestNeighbourIndex);
        
        samples[i] = 0.0;
        double kernelTotal = 0.0;
        for (int k=0; k<neighbourhood.size; k++)
        {
            double kernelValue = 1.0;
            for (int j=0; j<nDims; j++)
            {
                double delta = nearestNeighbourOffset[j] + static_cast<double>(neighbourhood.locs(k,j));
                int currentDimIndex = nearestNeighbour[j] + neighbourhood.locs(k,j);
                if (currentDimIndex < 0 || currentDimIndex >= dims[j])
                {
                    kernelValue = 0.0;
                    break;
                }
                else
                    kernelValue *= kernel->evaluate(delta);
            }
            
            if (kernelValue != 0.0)
            {
                long currentIndex = nearestNeighbourIndex + neighbourhood.offsets[k];
                samples[i] += kernelValue * original->at(currentIndex);
                kernelTotal += kernelValue;
            }
        }
        
        if (kernelTotal != 1.0)
            samples[i] /= kernelTotal;
    }
    
    return samples;
}
