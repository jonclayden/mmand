#include <RcppArmadillo.h>

#include "Resampler.hpp"

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
    return (array);
}

RcppExport SEXP get_neighbourhood (SEXP data_, SEXP width_)
{
BEGIN_RCPP
    Array *array = arrayFromData(data_);
    
    Neighbourhood neighbourhood = array->getNeighbourhood(as<int>(width_));
    
    return List::create(Named("width")=neighbourhood.width, Named("size")=neighbourhood.size, Named("locs")=neighbourhood.locs, Named("offsets")=neighbourhood.offsets);
END_RCPP
}

RcppExport SEXP resample (SEXP data_, SEXP kernel_, SEXP samplingLocations_)
{
BEGIN_RCPP
    Array *array = arrayFromData(data_);
    
    List kernelElements(kernel_);
    string kernelName = as<string>(kernelElements["name"]);
    Kernel *kernel;
    
    if (kernelName.compare("box") == 0)
        kernel = KernelGenerator::box();
    else if (kernelName.compare("triangle") == 0)
        kernel = KernelGenerator::triangle();
    else if (kernelName.compare("mitchell-netravali") == 0)
        kernel = KernelGenerator::mitchellNetravali(as<double>(kernelElements["B"]), as<double>(kernelElements["C"]));
    
    Resampler resampler(array, kernel);
    
    List samplingLocations(samplingLocations_);
    vector<int_vector> samplingVector(samplingLocations.length());
    for (int i=0; i<samplingLocations.length(); i++)
        samplingVector[i] = as<int_vector>(samplingLocations[i]);
    resampler.setSamplingLocations(samplingVector);
    
    vector<double> &samples = resampler.run();
    return wrap(samples);
END_RCPP
}
