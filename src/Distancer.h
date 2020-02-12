#ifndef _DISTANCER_H_
#define _DISTANCER_H_

#include "Array.h"

class Distancer
{
private:
    Array<double> *original;
    bool usePixdim;
    
public:
    Distancer (Array<double> * const original, const bool usePixdim)
        : original(original), usePixdim(usePixdim) {}
    
    ~Distancer ()
    {
        delete original;
    }
    
    Array<double> * run ();
};

#endif
