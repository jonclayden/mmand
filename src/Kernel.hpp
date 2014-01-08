#ifndef _KERNEL_HPP_
#define _KERNEL_HPP_

#include "Array.hpp"

// Base class for all kernels
// Defines support but will always evaluate to zero
class Kernel
{
protected:
    double supportMin, supportMax;
    
public:
    Kernel ()
        : supportMin(0.0), supportMax(0.0) {}
    Kernel (const double supportMin, const double supportMax)
        : supportMin(supportMin), supportMax(supportMax) {}
    
    virtual double evaluate (const double x) const { return 0.0; }
    
    double getSupportMin () const { return supportMin; }
    
    double getSupportMax () const { return supportMax; }
    
    bool isWithinSupport (const double x) const
    {
        double absX = fabs(x);
        return (absX >= supportMin && absX <= supportMax);
    }
};

// Discrete kernel
// This kernel is defined only at integral locations in a grid
// It is used when the image dimensions aren't changing
class DiscreteKernel : public Kernel
{
protected:
    Array *values;
    
public:
    DiscreteKernel (Array * const values)
        : values(values)
    {
        supportMin = 0.0;
        supportMax = 0.0;
        
        const std::vector<int> &dims = values->getDims();
        for (std::vector<int>::const_iterator i = dims.begin(); i != dims.end(); i++)
        {
            double currentSupportMax = floor(static_cast<double>(*i) / 2.0);
            if (currentSupportMax > supportMax)
                supportMax = currentSupportMax;
        }
    }
    
    ~DiscreteKernel ()
    {
        delete values;
    }
    
    Array * getArray () const { return values; }
};

// General polynomial kernel
// Evaluates to a polynomial function of location, within the support region
class PolynomialKernel : public Kernel
{
private:
    double term (const double x, const int i) const;
    
protected:
    int degree;
    arma::vec coefficients;
    
public:
    PolynomialKernel (const arma::vec &coefficients, const double supportMin, const double supportMax)
        : Kernel(supportMin,supportMax), coefficients(coefficients)
    {
        this->degree = coefficients.n_elem - 1;
    }
    
    double evaluate (const double x) const;
};

// Composite kernel
// A kernel comprised of several other kernels with different supports
class CompositeKernel : public Kernel
{
protected:
    std::vector<Kernel*> kernels;
    
public:
    CompositeKernel (const std::vector<Kernel*> &kernels)
        : kernels(kernels)
    {
        supportMin = R_PosInf;
        supportMax = R_NegInf;
        
        for (std::vector<Kernel*>::iterator i = this->kernels.begin(); i != this->kernels.end(); i++)
        {
            if ((*i)->getSupportMin() < supportMin)
                supportMin = (*i)->getSupportMin();
            if ((*i)->getSupportMax() > supportMax)
                supportMax = (*i)->getSupportMax();
        }
    }
    
    ~CompositeKernel ()
    {
        for (std::vector<Kernel*>::iterator i = this->kernels.begin(); i != this->kernels.end(); i++)
            delete *i;
    }
    
    double evaluate (const double x) const;
};

// Kernel generator
// A container for static functions to generate special case kernels easily
class KernelGenerator
{
public:
    static PolynomialKernel * box ();
    static PolynomialKernel * triangle ();
    static CompositeKernel * mitchellNetravali (const double B, const double C);
};

#endif
