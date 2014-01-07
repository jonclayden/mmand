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
    const long neighbourhoodSize = kernelNeighbourhood.size;
    
    const int_vector &dims = original->getDims();
    int nDims = original->getNDims();
    long nSamples = original->size();
    samples.resize(nSamples);
    
    double kernelSum = 0.0;
    for (long k=0; k<neighbourhoodSize; k++)
        kernelSum += kernelArray->at(k);
    
    int_vector currentLoc(nDims);
    for (long i=0; i<nSamples; i++)
    {
        if (!meetsRestrictions(i))
        {
            samples[i] = original->at(i);
            continue;
        }
        
        original->expandIndex(i, currentLoc);
        
        dbl_vector values;
        for (long k=0; k<neighbourhoodSize; k++)
        {
            bool validLoc = true;
            for (int j=0; j<nDims; j++)
            {
                int currentDimIndex = currentLoc[j] + sourceNeighbourhood.locs(k,j);
                if (currentDimIndex < 0 || currentDimIndex >= dims[j])
                    validLoc = false;
            }
            
            if (validLoc && !R_IsNA(kernelArray->at(k)))
            {
                switch (elementOp)
                {
                    case PlusOp:
                    values.push_back(kernelArray->at(k) + original->at(i+sourceNeighbourhood.offsets[k]));
                    break;
                    
                    case MinusOp:
                    values.push_back(kernelArray->at(k) - original->at(i+sourceNeighbourhood.offsets[k]));
                    break;
                    
                    case MultiplyOp:
                    values.push_back(kernelArray->at(k) * original->at(i+sourceNeighbourhood.offsets[k]));
                    break;
                    
                    case IdentityOp:
                    if (kernelArray->at(k) != 0.0)
                        values.push_back(original->at(i+sourceNeighbourhood.offsets[k]));
                    break;
                    
                    case OneOp:
                    if (kernelArray->at(k) != 0.0)
                        values.push_back(1.0);
                    break;
                    
                    case ZeroOp:
                    if (kernelArray->at(k) != 0.0)
                        values.push_back(0.0);
                    break;
                }
            }
        }
        
        if (values.size() == 0)
            samples[i] = NA_REAL;
        else if (neighbourhoodSize == 1)
            samples[i] = values[0];
        else
        {
            switch (mergeOp)
            {
                case SumOp:
                {
                    double sum = 0.0;
                    for (int l=0; l<values.size(); l++)
                        sum += values[l];
                    samples[i] = sum;
                    break;
                }
                
                case MinOp:
                samples[i] = *(std::min_element(values.begin(), values.end()));
                break;
                
                case MaxOp:
                samples[i] = *(std::max_element(values.begin(), values.end()));
                break;
                
                case MeanOp:
                {
                    double sum = 0.0;
                    for (int l=0; l<values.size(); l++)
                        sum += values[l];
                    samples[i] = sum / static_cast<double>(values.size());
                    break;
                }
                
                case MedianOp:
                {
                    int middleIndex = values.size() / 2;
                    std::partial_sort(values.begin(), values.begin()+middleIndex, values.end());
                    if (values.size() % 2 == 0)
                        samples[i] = (values[middleIndex-1] + values[middleIndex]) / 2.0;
                    else
                        samples[i] = values[middleIndex];
                    break;
                }
            }
        }
    }
    
    return samples;
}
