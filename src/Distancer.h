#ifndef _DISTANCER_H_
#define _DISTANCER_H_

#include "Array.h"

typedef std::vector<double> dbl_vector;

class Distancer
{
private:
    Array<double> *original;
    
public:
    Distancer (Array<double> * const original)
        : original(original) {}
    
    const dbl_vector & run ();
};

#endif