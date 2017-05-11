#ifndef _MORPHER_H_
#define _MORPHER_H_

#include "Array.h"
#include "Kernel.h"

typedef std::vector<double> dbl_vector;
typedef std::vector<int>    int_vector;

enum ElementOp { PlusOp, MinusOp, MultiplyOp, IdentityOp, OneOp, ZeroOp, EqualOp };
enum MergeOp { SumOp, MinOp, MaxOp, MeanOp, MedianOp, AllOp, AnyOp };

class Morpher
{
private:
    Array<double> *original;
    DiscreteKernel *kernel;
    
    ElementOp elementOp;
    MergeOp mergeOp;
    
    Neighbourhood immediateNeighbourhood;
    int_vector currentLoc;
    
    dbl_vector includedValues, excludedValues;
    int_vector includedNeighbourhoods, excludedNeighbourhoods;
    
    bool renormalise;
    
    dbl_vector values;
    dbl_vector samples;
    
    bool meetsRestrictions (const size_t n);
    
    void resetValues ();
    void accumulateValue (double value);
    double mergeValues ();
    
public:
    Morpher (Array<double> * const original, DiscreteKernel * const kernel, const ElementOp elementOp, const MergeOp mergeOp)
        : original(original), kernel(kernel), elementOp(elementOp), mergeOp(mergeOp), renormalise(true)
    {
        this->immediateNeighbourhood = original->getNeighbourhood(3);
    }
    
    ~Morpher ()
    {
        delete original;
        delete kernel;
    }
    
    void setValidNeighbourhoods (const int_vector &include, const int_vector &exclude)
    {
        this->includedNeighbourhoods = include;
        this->excludedNeighbourhoods = exclude;
    }
    
    void setValidValues (const dbl_vector &include, const dbl_vector &exclude)
    {
        this->includedValues = include;
        this->excludedValues = exclude;
    }
    
    void shouldRenormalise (const bool renormalise)
    {
        this->renormalise = renormalise;
    }
    
    dbl_vector & run ();
};

#endif
