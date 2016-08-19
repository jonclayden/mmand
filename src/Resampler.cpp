#include <RcppEigen.h>

#include "Resampler.h"

double Resampler::interpolate (Array<double>::ConstIterator begin, Array<double>::ConstIterator end, const double loc)
{
    std::vector<double> data;
    Array<double>::ConstIterator it;
    
    if (presharpen)
    {
        int i;
        std::vector<double> coefs;
        for (i=0, it=begin; ; i++)
        {
            coefs.push_back(c / (b - (i==0 ? 0 : a*coefs[i-1])));
            data.push_back((*it - (i==0 ? 0 : a*data[i-1])) / (b - (i==0 ? 0 : a*coefs[i-1])));
            if (++it > end)
                break;
        }
        
        const int len = data.size();
        for (i=(len-1); i>=0; i--)
            data[i] -= (i==(len-1) ? 0 : coefs[i]*data[i+1]);
    }
    else
        std::copy(begin, end+1, data.begin());
    
    double result = 0.0;
    for (int i=0; i<data.size(); i++)
        result += data[i] * kernel->evaluate(static_cast<double>(i) - loc);
    
    return result;
}
