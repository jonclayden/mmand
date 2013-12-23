#ifndef _MORPHER_HPP_
#define _MORPHER_HPP_

#include "Array.hpp"
#include "Kernel.hpp"

typedef std::vector<double> dbl_vector;
typedef std::vector<int>    int_vector;

class Morpher
{
private:
    Array *original;
    DiscreteKernel *kernel;
    
    Neighbourhood immediateNeighbourhood;
    
    dbl_vector includedValues, excludedValues;
    int_vector includedNeighbourhoods, excludedNeighbourhoods;
    
    dbl_vector samples;
    
    bool meetsRestrictions (const long n);
    
public:
    Morpher (Array * const original, DiscreteKernel * const kernel)
        : original(original), kernel(kernel)
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
