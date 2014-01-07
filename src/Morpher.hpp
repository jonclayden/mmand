#ifndef _MORPHER_HPP_
#define _MORPHER_HPP_

#include "Array.hpp"
#include "Kernel.hpp"

typedef std::vector<double> dbl_vector;
typedef std::vector<int>    int_vector;

enum ElementOp { PlusOp, MinusOp, MultiplyOp, IdentityOp, OneOp, ZeroOp };
enum MergeOp { SumOp, MinOp, MaxOp, MeanOp, MedianOp };

class Morpher
{
private:
    Array *original;
    DiscreteKernel *kernel;
    
    ElementOp elementOp;
    MergeOp mergeOp;
    
    Neighbourhood immediateNeighbourhood;
    
    dbl_vector includedValues, excludedValues;
    int_vector includedNeighbourhoods, excludedNeighbourhoods;
    
    dbl_vector samples;
    
    bool meetsRestrictions (const long n);
    
public:
    Morpher (Array * const original, DiscreteKernel * const kernel, const ElementOp elementOp, const MergeOp mergeOp)
        : original(original), kernel(kernel), elementOp(elementOp), mergeOp(mergeOp)
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
    
    dbl_vector & run ();
};

#endif
