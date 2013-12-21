#ifndef _KERNEL_HPP_
#define _KERNEL_HPP_

// Base class for all kernels
// Defines support but will always evaluate to zero
class Kernel
{
protected:
    double supportMin, supportMax;
    bool usingCache;
    std::map<double,double> cache;
    
    std::pair<bool,double> checkCache (const double key)
    {
        if (!usingCache || cache.count(key) == 0)
            return std::pair<bool,double>(false, 0.0);
        else
            return std::pair<bool,double>(true, cache[key]);
    }
    
    void updateCache (const double key, const double value)
    {
        if (usingCache && cache.count(key) == 0)
            cache.insert(std::pair<double,double>(key,value));
    }
    
public:
    Kernel ()
        : supportMin(0.0), supportMax(0.0), usingCache(false) {}
    Kernel (const bool usingCache)
        : supportMin(0.0), supportMax(0.0), usingCache(usingCache) {}
    Kernel (const double supportMin, const double supportMax)
        : supportMin(supportMin), supportMax(supportMax), usingCache(false) {}
    Kernel (const double supportMin, const double supportMax, const bool usingCache)
        : supportMin(supportMin), supportMax(supportMax), usingCache(usingCache) {}
    
    virtual double evaluate (const double x) { return 0.0; }
    
    double getSupportMin () { return supportMin; }
    
    double getSupportMax () { return supportMax; }
    
    bool isWithinSupport (const double x)
    {
        double absX = fabs(x);
        return (absX >= supportMin && absX <= supportMax);
    }
};

// General polynomial kernel
// Evaluates to a polynomial function of location, within the support region
class PolynomialKernel : public Kernel
{
private:
    double term (const double x, const int i);
    
protected:
    int degree;
    arma::vec coefficients;
    
public:
    PolynomialKernel (const arma::vec &coefficients, const double supportMin, const double supportMax)
        : Kernel(supportMin,supportMax), coefficients(coefficients)
    {
        this->degree = coefficients.n_elem - 1;
    }
    
    double evaluate (const double x);
};

// Composite kernel
// A kernel comprised of several other kernels with different supports
class CompositeKernel : public Kernel
{
protected:
    std::vector<Kernel*> kernels;
    
public:
    CompositeKernel (const std::vector<Kernel*> &kernels)
        : Kernel(true), kernels(kernels)
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
    
    double evaluate (const double x);
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
