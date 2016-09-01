#ifndef _KERNEL_H_
#define _KERNEL_H_

#include "Array.h"

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
    Array<double> *values;
    
public:
    DiscreteKernel (Array<double> * const values)
        : values(values)
    {
        supportMin = 0.0;
        supportMax = 0.0;
        
        const std::vector<int> &dims = values->getDimensions();
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
    
    Array<double> * getArray () const { return values; }
};

template <int N>
struct PolynomialEvaluator
{
private:
    const Eigen::VectorXd *coefficients;
    const PolynomialEvaluator<N-1> child;

public:
    PolynomialEvaluator (const Eigen::VectorXd *coefficients)
        : coefficients(coefficients), child(coefficients) {}

    double operator() (const double x) const { return (*coefficients)[N] + x * child(x); }
};

template <>
struct PolynomialEvaluator<0>
{
private:
    const Eigen::VectorXd *coefficients;

public:
    PolynomialEvaluator (const Eigen::VectorXd *coefficients)
        : coefficients(coefficients) {}

    double operator() (const double x) const { return (*coefficients)[0]; }
};

// General polynomial kernel
// Evaluates to a polynomial function of location, within the support region
template <int Degree>
class PolynomialKernel : public Kernel
{
protected:
    Eigen::VectorXd coefficients;
    PolynomialEvaluator<Degree> evaluator;
    
public:
    PolynomialKernel (const Eigen::VectorXd &coefficients, const double supportMin, const double supportMax)
        : Kernel(supportMin,supportMax), coefficients(coefficients), evaluator(&this->coefficients)
    {
        std::reverse(this->coefficients.data(), this->coefficients.data()+this->coefficients.size());
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
    static PolynomialKernel<0> * box ();
    static PolynomialKernel<1> * triangle ();
    static CompositeKernel * mitchellNetravali (const double B, const double C);
};

#endif
