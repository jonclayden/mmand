#include <RcppArmadillo.h>

#include "Morpher.hpp"

bool Morpher::meetsRestrictions (const long n)
{
    double value = original->at(n);
    
    if (includedValues.size() > 0)
    {
        bool found = false;
        for (int i=0; i<includedValues.size(); i++)
        {
            if (value == includedValues[i])
                found = true;
        }
        
        if (!found)
            return false;
    }
    else if (excludedValues.size() > 0)
    {
        for (int i=0; i<excludedValues.size(); i++)
        {
            if (value == includedValues[i])
                return false;
        }
    }
    
    if (includedNeighbourhoods.size() > 0 || excludedNeighbourhoods.size() > 0)
    {
        int nDims = original->getNDims();
        int_vector currentLoc(nDims);
        original->expandIndex(n, currentLoc);
        const std::vector<int> &dims = original->getDims();
        
        int nNeighbours = 0;
        int neighbourhoodCentre = (immediateNeighbourhood.size - 1) / 2;
        for (int k=0; k<immediateNeighbourhood.size; k++)
        {
            if (k == neighbourhoodCentre)
                continue;
            
            bool validLoc = true;
            for (int j=0; j<nDims; j++)
            {
                int currentDimIndex = currentLoc[j] + immediateNeighbourhood.locs(k,j);
                if (currentDimIndex < 0 || currentDimIndex >= dims[j])
                    validLoc = false;
            }
            
            if (validLoc && original->at(n+immediateNeighbourhood.offsets[k]) != 0.0)
                nNeighbours++;
        }
        
        if (includedNeighbourhoods.size() > 0)
        {
            bool found = false;
            for (int i=0; i<includedNeighbourhoods.size(); i++)
            {
                if (nNeighbours == includedNeighbourhoods[i])
                    found = true;
            }
        
            if (!found)
                return false;
        }
        else if (excludedNeighbourhoods.size() > 0)
        {
            for (int i=0; i<excludedNeighbourhoods.size(); i++)
            {
                if (nNeighbours == excludedNeighbourhoods[i])
                    return false;
            }
        }
    }
    
    return true;
}

std::vector<double> & Morpher::run ()
{
    Array * kernelArray = kernel->getArray();
    const Neighbourhood &kernelNeighbourhood = kernelArray->getNeighbourhood();
    const Neighbourhood &sourceNeighbourhood = original->getNeighbourhood(kernelArray->getDims());
    
    const std::vector<int> &dims = original->getDims();
    int nDims = original->getNDims();
    long nSamples = original->size();
    samples.resize(nSamples);
    
    double kernelSum = 0.0;
    for (int k=0; k<kernelNeighbourhood.size; k++)
        kernelSum += kernelArray->at(k);
    
    int_vector currentLoc(nDims);
    for (long i=0; i<nSamples; i++)
    {
        if (!meetsRestrictions(i))
        {
            samples[i] = original->at(i);
            continue;
        }
        
        samples[i] = 0.0;
        double visitedKernelSum = 0.0;
        
        original->expandIndex(i, currentLoc);
        
        for (int k=0; k<sourceNeighbourhood.size; k++)
        {
            bool validLoc = true;
            for (int j=0; j<nDims; j++)
            {
                int currentDimIndex = currentLoc[j] + sourceNeighbourhood.locs(k,j);
                if (currentDimIndex < 0 || currentDimIndex >= dims[j])
                    validLoc = false;
            }
            
            if (validLoc)
            {
                if (kernel->isBrush && kernelArray->at(k) != NA_REAL)
                {
                    if (kernel->isEraser && kernelArray->at(k) != 0.0)
                        samples[i+sourceNeighbourhood.offsets[k]] = 0.0;
                    else if (!kernel->isEraser)
                        samples[i+sourceNeighbourhood.offsets[k]] = kernelArray->at(k);
                }
                else if (!kernel->isBrush)
                {
                    samples[i] += kernelArray->at(k) * original->at(i+sourceNeighbourhood.offsets[k]);
                    visitedKernelSum += kernelArray->at(k);
                }
            }
        }
        
        if (!kernel->isBrush)
            samples[i] *= (visitedKernelSum == 0.0) ? 1.0 : (kernelSum / visitedKernelSum);
    }
    
    return samples;
}
