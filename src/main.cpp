#include <RcppArmadillo.h>

#include "Resampler.hpp"
#include "Morpher.hpp"

using namespace Rcpp;
using namespace std;

Array * arrayFromData (SEXP data_)
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
        
    Array *array = new Array(as<dbl_vector>(data), dim);
    return array;
}

Kernel * kernelFromElements (SEXP kernel_)
{
    List kernelElements(kernel_);
    string kernelName = as<string>(kernelElements["name"]);
    Kernel *kernel;
    
    if (kernelName.compare("box") == 0)
        kernel = KernelGenerator::box();
    else if (kernelName.compare("triangle") == 0)
        kernel = KernelGenerator::triangle();
    else if (kernelName.compare("mitchell-netravali") == 0)
        kernel = KernelGenerator::mitchellNetravali(as<double>(kernelElements["B"]), as<double>(kernelElements["C"]));
    else if (kernelName.compare("discrete") == 0)
    {
        Array *kernelArray = arrayFromData(kernelElements["values"]);
        kernel = new DiscreteKernel(kernelArray, as<bool>(kernelElements["brush"]), as<bool>(kernelElements["eraser"]));
    }
    
    return kernel;
}

RcppExport SEXP get_neighbourhood (SEXP data_, SEXP width_)
{
BEGIN_RCPP
    Array *array = arrayFromData(data_);
    Neighbourhood neighbourhood = array->getNeighbourhood(as<int>(width_));
    
    delete array;
    
    return List::create(Named("widths")=neighbourhood.widths, Named("size")=neighbourhood.size, Named("locs")=neighbourhood.locs, Named("offsets")=neighbourhood.offsets);
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
    Array *array = arrayFromData(data_);
    Kernel *kernel = kernelFromElements(kernel_);
    Resampler resampler(array, kernel);
    
    List samplingScheme(samplingScheme_);
    string schemeType = as<string>(samplingScheme["type"]);
    SamplingScheme *sampler;
    
    if (schemeType.compare("general") == 0)
        sampler = new GeneralSamplingScheme(as<arma::mat>(samplingScheme["points"]));
    else if (schemeType.compare("grid") == 0)
    {
        List points = samplingScheme["points"];
        vector<dbl_vector> samplingVector(points.length());
        for (int i=0; i<points.length(); i++)
            samplingVector[i] = as<dbl_vector>(points[i]);
        sampler = new GriddedSamplingScheme(samplingVector);
    }
    
    resampler.setSamplingScheme(sampler);
    vector<double> &samples = resampler.run();
    return wrap(samples);
END_RCPP
}

RcppExport SEXP morph (SEXP data_, SEXP kernel_, SEXP restrictions_)
{
BEGIN_RCPP
    Array *array = arrayFromData(data_);
    
    List kernelElements(kernel_);
    Array *kernelArray = arrayFromData(kernelElements["values"]);
    DiscreteKernel *kernel = new DiscreteKernel(kernelArray, as<bool>(kernelElements["brush"]), as<bool>(kernelElements["eraser"]));
    
    Morpher morpher(array, kernel);
    
    List restrictions(restrictions_);
    morpher.setValidNeighbourhoods(as<int_vector>(restrictions["nNeighbours"]), as<int_vector>(restrictions["nNeighboursNot"]));
    morpher.setValidValues(as<dbl_vector>(restrictions["value"]), as<dbl_vector>(restrictions["valueNot"]));
    vector<double> &samples = morpher.run();
    return wrap(samples);
END_RCPP
}
