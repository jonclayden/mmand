#include <RcppArmadillo.h>

#include "Resampler.hpp"

std::vector<double> & Resampler::run ()
{
    // Precalculate index deltas for the neighbourhood constituting the support region for the kernel
    int neighbourhoodWidth = static_cast<int>(floor(2.0*kernel->getSupportMax()));
    Neighbourhood neighbourhood = original->getNeighbourhood(neighbourhoodWidth);
    
    int nDims = sampler->getNDims();
    int nSamples = sampler->getNSamples();
    samples.resize(nSamples);
    
    std::vector<int> nearestNeighbour(nDims);
    std::vector<double> nearestNeighbourOffset(nDims);
    size_type nearestNeighbourIndex;
    
    for (int i=0; i<nSamples; i++)
    {
        for (int j=0; j<nDims; j++)
        {
            nearestNeighbour[j] = static_cast<int>(round(sampler->at(i,j)));
            nearestNeighbourOffset[j] = static_cast<double>(nearestNeighbour[j]) - sampler->at(i,j);
        }
        original->flattenIndex(nearestNeighbour, nearestNeighbourIndex);
        
        samples[i] = 0.0;
        for (int k=0; k<neighbourhood.size; k++)
        {
            size_type currentIndex = nearestNeighbourIndex + neighbourhood.offsets[k];
            if (currentIndex >= 0 && currentIndex < original->size())
            {
                double sampleContribution = original->at(currentIndex);
                for (int j=0; j<nDims; j++)
                {
                    double delta = nearestNeighbourOffset[j] + static_cast<double>(neighbourhood.locs(k,j));
                    sampleContribution *= kernel->evaluate(delta);
                }
                samples[i] += sampleContribution;
            }
        }
    }
    
    return samples;
}
