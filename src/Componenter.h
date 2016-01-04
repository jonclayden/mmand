#ifndef _COMPONENTER_H_
#define _COMPONENTER_H_

#include "Array.h"
#include "Kernel.h"

#include "lemon/smart_graph.h"

class Componenter
{
private:
    Array<double> *original;
    DiscreteKernel *kernel;
    
    lemon::SmartGraph connections;
    std::vector<int> labels;
    
public:
    Componenter (Array<double> * const original, DiscreteKernel * const kernel)
        : original(original), kernel(kernel) {}
    
    ~Componenter ()
    {
        delete original;
        delete kernel;
    }
    
    std::vector<int> & run ();
};

#endif
