#include <RcppEigen.h>

#include "Kernel.h"

template <int Degree>
double PolynomialKernel<Degree>::evaluate (const double x) const
{
    if (!isWithinSupport(x))
        return 0.0;
    else
        return evaluator(fabs(x));
}

double CompositeKernel::evaluate (const double x) const
{
    if (kernels.size() == 0)
        return 0.0;
    else if (!isWithinSupport(x))
        return 0.0;
    else
    {
        for (std::vector<Kernel*>::const_iterator i = kernels.begin(); i != kernels.end(); i++)
        {
            if ((*i)->isWithinSupport(x))
                return (*i)->evaluate(x);
        }
        
        return 0.0;
    }
}

// Box kernel: constant value of 1.0, support of 0.5
// Used for nearest-neighbour sampling
PolynomialKernel<0> * KernelGenerator::box ()
{
    Eigen::VectorXd coefficients(1);
    coefficients[0] = 1.0;
    return new PolynomialKernel<0>(coefficients, 0.0, 0.5);
}

// Triangle kernel: linear slope downwards from 0 to 1
// Used for linear interpolation
PolynomialKernel<1> * KernelGenerator::triangle ()
{
    Eigen::VectorXd coefficients(2);
    coefficients[0] = 1.0;
    coefficients[1] = -1.0;
    return new PolynomialKernel<1>(coefficients, 0.0, 1.0);
}

// Mitchell-Netravali family of cubic kernels
CompositeKernel * KernelGenerator::mitchellNetravali (const double B, const double C)
{
    Eigen::VectorXd coefficients1 = Eigen::VectorXd::Zero(4);
    coefficients1[0] = 1.0 - B/3.0;
    coefficients1[2] = -3.0 + 2.0*B + C;
    coefficients1[3] = 2.0 - 1.5*B - C;
    PolynomialKernel<3> *kernel1 = new PolynomialKernel<3>(coefficients1, 0.0, 1.0);
    
    Eigen::VectorXd coefficients2 = Eigen::VectorXd::Zero(4);
    coefficients2[0] = 4.0*B/3.0 + 4.0*C;
    coefficients2[1] = -2.0*B - 8.0*C;
    coefficients2[2] = B + 5.0*C;
    coefficients2[3] = -B/6.0 - C;
    PolynomialKernel<3> *kernel2 = new PolynomialKernel<3>(coefficients2, 1.0, 2.0);
    
    std::vector<Kernel*> kernels;
    kernels.push_back(kernel1);
    kernels.push_back(kernel2);
    return new CompositeKernel(kernels);
}
