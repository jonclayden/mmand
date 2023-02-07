#ifndef KERNEL_ARRAY_H_
#define KERNEL_ARRAY_H_

#include "Array.h"

#include <functional>

enum struct ElementOp { Plus, Minus, Multiply, Identity, One, Zero, Equal };
enum struct MergeOp { Sum, Min, Max, Mean, Median, All, Any };

// This kernel is defined only at integral locations in a grid.
// It is used when the image dimensions aren't changing
class KernelArray
{
protected:
    std::vector<int> spaceDims;
    std::vector<double> values;
    std::vector<ptrdiff_t> offsets;
    
    std::function<double(double&,double&)> combine;
    
public:
    KernelArray (Array<double> * const kernel, const Neighbourhood &neighbourhood, ElementOp op, MergeOp merge)
    {
        if (kernel->size() != neighbourhood.size)
            throw Rcpp::exception("Kernel array and target neighbourhood should match in size");
        
        // Where the kernel is NA, the element's contribution can be ignored.
        // For certain operations this is true for another kernel value too
        double identity = NA_REAL;
        if (op == ElementOp::Identity || op == ElementOp::One || op == ElementOp::Zero)
            identity = 0.0;
        else if (op == ElementOp::Multiply && (merge == MergeOp::Sum || merge == MergeOp::Any))
            identity = 0.0;
        
        for (size_t i=0; i<neighbourhood.size; i++)
        {
            if (!R_IsNA(kernel->at(i)) && kernel->at(i) != identity)
            {
                values.push_back(kernel->at(i));
                offsets.push_back(neighbourhood.offsets[i]);
            }
        }
        
        switch (op)
        {
            case ElementOp::Plus:       combine = std::plus<double>();                                      break;
            case ElementOp::Minus:      combine = std::minus<double>();                                     break;
            case ElementOp::Multiply:   combine = std::multiplies<double>();                                break;
            case ElementOp::Identity:   combine = [](double &x, double &k) { return x; };                   break;
            case ElementOp::One:        combine = [](double &x, double &k) { return 1.0; };                 break;
            case ElementOp::Zero:       combine = [](double &x, double &k) { return 0.0; };                 break;
            case ElementOp::Equal:      combine = [](double &x, double &k) { return x==k ? 1.0 : 0.0; };    break;
        }
    }
};

#endif
