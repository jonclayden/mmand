#ifndef _DISTANCER_H_
#define _DISTANCER_H_

#include "Array.h"

class Distancer
{
private:
    Array<double> *original;
    
public:
    Distancer (Array<double> * const original)
        : original(original) {}
    
    ~Distancer ()
    {
        delete original;
    }
    
    Array<double> * run ();
};

#endif
