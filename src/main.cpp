#include <Rcpp.h>

#include "Componenter.h"
#include "Resampler.h"
#include "Morpher.h"

using namespace Rcpp;
using namespace std;

Array<double> * arrayFromData (SEXP data_)
{
    NumericVector data(data_);
    
    int_vector dim;
    if (data.hasAttribute("dim"))
        dim = as<int_vector>(data.attr("dim"));
    else
    {
        dim = int_vector(1);
        dim[0] = data.length();
    }
        
    Array<double> *array = new Array<double>(dim, as<dbl_vector>(data));
    return array;
}

Kernel * kernelFromElements (SEXP kernel_)
{
    List kernelElements(kernel_);
    string kernelName = as<string>(kernelElements["name"]);
    Kernel *kernel = NULL;
    
    if (kernelName.compare("box") == 0)
        kernel = KernelGenerator::box();
    else if (kernelName.compare("triangle") == 0)
        kernel = KernelGenerator::triangle();
    else if (kernelName.compare("mitchell-netravali") == 0)
        kernel = KernelGenerator::mitchellNetravali(as<double>(kernelElements["B"]), as<double>(kernelElements["C"]));
    else if (kernelName.compare("lanczos") == 0)
        kernel = KernelGenerator::lanczos();
    
    return kernel;
}

RcppExport SEXP is_binary (SEXP data_)
{
BEGIN_RCPP
    NumericVector data(data_);
    bool isBinary = true;
    double nonzeroValue = NA_REAL;
    
    for (int i=0; i<data.length(); i++)
    {
        if (data[i] != 0.0)
        {
            if (R_IsNA(nonzeroValue))
                nonzeroValue = data[i];
            else if (nonzeroValue != data[i])
            {
                isBinary = false;
                break;
            }
        }
    }
    
    RObject result = wrap(isBinary);
    if (isBinary)
        result.attr("value") = nonzeroValue;
    return result;
END_RCPP
}

RcppExport SEXP is_symmetric (SEXP data_)
{
BEGIN_RCPP
    NumericVector data(data_);
    bool isSymmetric = true;
    
    for (int i=0; i<data.length()/2; i++)
    {
        if (!R_IsNA(data[i]) && data[i] != data[data.length()-i-1])
        {
            isSymmetric = false;
            break;
        }
    }
    
    return wrap(isSymmetric);
END_RCPP
}

RcppExport SEXP get_neighbourhood (SEXP data_, SEXP width_)
{
BEGIN_RCPP
    Array<double> *array = arrayFromData(data_);
    Neighbourhood neighbourhood = array->getNeighbourhood(as<int_vector>(width_));
    
    delete array;
    
    const int size = static_cast<int>(neighbourhood.size);
    std::vector<double> offsets(size);
    for (int i=0; i<size; i++)
        offsets[i] = static_cast<double>(neighbourhood.offsets[i]);
    
    return List::create(Named("widths")=neighbourhood.widths, Named("size")=size, Named("locs")=neighbourhood.locs, Named("offsets")=offsets);
END_RCPP
}

RcppExport SEXP sample_kernel (SEXP kernel_, SEXP values_)
{
BEGIN_RCPP
    Kernel *kernel = kernelFromElements(kernel_);
    NumericVector values(values_);
    NumericVector result(values.length());
    
    for (int i=0; i<values.length(); i++)
        result[i] = kernel->evaluate(values[i]);
    
    delete kernel;
    
    return result;
END_RCPP
}

RcppExport SEXP resample (SEXP data_, SEXP kernel_, SEXP samplingScheme_)
{
BEGIN_RCPP
    Array<double> *array = arrayFromData(data_);
    Kernel *kernel = kernelFromElements(kernel_);
    Resampler resampler(array, kernel);

    List samplingScheme(samplingScheme_);
    string schemeType = as<string>(samplingScheme["type"]);
    
    if (schemeType.compare("general") == 0)
    {
        NumericMatrix points = samplingScheme["points"];
        const dbl_vector &samples = resampler.run(points);
        return wrap(samples);
    }
    else if (schemeType.compare("grid") == 0)
    {
        List points = samplingScheme["points"];
        vector<dbl_vector> samplingVector(points.length());
        for (int i=0; i<points.length(); i++)
            samplingVector[i] = as<dbl_vector>(points[i]);
        const dbl_vector &samples = resampler.run(samplingVector);
        return wrap(samples);
    }
    else
        throw std::runtime_error("Scheme type unsupported");
END_RCPP
}

RcppExport SEXP morph (SEXP data_, SEXP kernel_, SEXP elementOp_, SEXP mergeOp_, SEXP restrictions_)
{
BEGIN_RCPP
    Array<double> *array = arrayFromData(data_);
    
    Array<double> *kernelArray = arrayFromData(kernel_);
    DiscreteKernel *kernel = new DiscreteKernel(kernelArray);
    
    const string elementOpString = as<string>(elementOp_);
    ElementOp elementOp;
    if (elementOpString.compare("+") == 0)
        elementOp = PlusOp;
    else if (elementOpString.compare("-") == 0)
        elementOp = MinusOp;
    else if (elementOpString.compare("*") == 0)
        elementOp = MultiplyOp;
    else if (elementOpString.compare("i") == 0)
        elementOp = IdentityOp;
    else if (elementOpString.compare("1") == 0)
        elementOp = OneOp;
    else if (elementOpString.compare("0") == 0)
        elementOp = ZeroOp;
    else if (elementOpString.compare("==") == 0)
        elementOp = EqualOp;
    else
        throw new runtime_error("Unsupported element operation specified");
    
    const string mergeOpString = as<string>(mergeOp_);
    MergeOp mergeOp;
    if (mergeOpString.compare("sum") == 0)
        mergeOp = SumOp;
    else if (mergeOpString.compare("min") == 0)
        mergeOp = MinOp;
    else if (mergeOpString.compare("max") == 0)
        mergeOp = MaxOp;
    else if (mergeOpString.compare("mean") == 0)
        mergeOp = MeanOp;
    else if (mergeOpString.compare("median") == 0)
        mergeOp = MedianOp;
    else if (mergeOpString.compare("all") == 0)
        mergeOp = AllOp;
    else if (mergeOpString.compare("any") == 0)
        mergeOp = AnyOp;
    else
        throw new runtime_error("Unsupported merge operation specified");
    
    Morpher morpher(array, kernel, elementOp, mergeOp);
    
    List restrictions(restrictions_);
    morpher.setValidNeighbourhoods(as<int_vector>(restrictions["nNeighbours"]), as<int_vector>(restrictions["nNeighboursNot"]));
    morpher.setValidValues(as<dbl_vector>(restrictions["value"]), as<dbl_vector>(restrictions["valueNot"]));
    vector<double> &samples = morpher.run();
    return wrap(samples);
END_RCPP
}

RcppExport SEXP connected_components (SEXP data_, SEXP kernel_)
{
BEGIN_RCPP
    Array<double> *array = arrayFromData(data_);
    
    Array<double> *kernelArray = arrayFromData(kernel_);
    DiscreteKernel *kernel = new DiscreteKernel(kernelArray);
    
    Componenter componenter(array, kernel);
    vector<int> &labels = componenter.run();
    return wrap(labels);
END_RCPP
}
