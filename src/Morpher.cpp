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
            if (value == excludedValues[i])
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

void Morpher::resetValues ()
{
    values.clear();
    if (mergeOp == MinOp)
        values.push_back(R_PosInf);
    else if (mergeOp == MaxOp)
        values.push_back(R_NegInf);
}

void Morpher::accumulateValue (double value)
{
    if (R_IsNA(value))
        return;
    
    if (mergeOp == MinOp && value < values[0])
        values[0] = value;
    else if (mergeOp == MaxOp && value > values[0])
        values[0] = value;
    else if (mergeOp != MinOp && mergeOp != MaxOp)
        values.push_back(value);
}

double Morpher::mergeValues ()
{
    if (values.size() == 0)
        return NA_REAL;
    else if (values.size() == 1)
        return values[0];
    else
    {
        switch (mergeOp)
        {
            case SumOp:
            {
                double sum = 0.0;
                for (int l=0; l<values.size(); l++)
                    sum += values[l];
                return sum;
            }
            
            case MeanOp:
            {
                double sum = 0.0;
                for (int l=0; l<values.size(); l++)
                    sum += values[l];
                return (sum / static_cast<double>(values.size()));
            }
            
            case MedianOp:
            {
                int middleIndex = values.size() / 2;
                std::partial_sort(values.begin(), values.begin()+middleIndex, values.end());
                if (values.size() % 2 == 0)
                    return ((values[middleIndex-1] + values[middleIndex]) / 2.0);
                else
                    return values[middleIndex];
            }
        }
    }
    
    return NA_REAL;
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
        
        resetValues();
        original->expandIndex(i, currentLoc);
        
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
                    accumulateValue(kernelArray->at(k) + original->at(i+sourceNeighbourhood.offsets[k]));
                    break;
                    
                    case MinusOp:
                    accumulateValue(kernelArray->at(k) - original->at(i+sourceNeighbourhood.offsets[k]));
                    break;
                    
                    case MultiplyOp:
                    accumulateValue(kernelArray->at(k) * original->at(i+sourceNeighbourhood.offsets[k]));
                    break;
                    
                    case IdentityOp:
                    if (kernelArray->at(k) != 0.0)
                        accumulateValue(original->at(i+sourceNeighbourhood.offsets[k]));
                    break;
                    
                    case OneOp:
                    if (kernelArray->at(k) != 0.0)
                        accumulateValue(1.0);
                    break;
                    
                    case ZeroOp:
                    if (kernelArray->at(k) != 0.0)
                        accumulateValue(0.0);
                    break;
                }
            }
        }
        
        samples[i] = mergeValues();
    }
    
    return samples;
}
