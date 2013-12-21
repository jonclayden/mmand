#ifndef _ARRAY_HPP_
#define _ARRAY_HPP_

struct Neighbourhood
{
    long size;
    std::vector<int> widths;
    arma::Mat<int> locs;
    std::vector<long> offsets;
};

class Array
{
private:
    std::vector<double> data;
    std::vector<int> dims;
    int nDims;

public:
    Array () {}
    
    Array (const std::vector<double> &data, const std::vector<int> &dims)
        : data(data), dims(dims)
    {
        nDims = dims.size();
    }
    
    const double & at (long n) { return data.at(n); }
    
    long size () { return data.size(); }
    
    std::vector<int> & getDims () { return dims; }
    
    Neighbourhood getNeighbourhood ();
    
    Neighbourhood getNeighbourhood (const int width);
    
    Neighbourhood getNeighbourhood (const std::vector<int> &widths);
    
    void flattenIndex (const std::vector<int> &loc, long &result);
    
    void expandIndex (const long &loc, std::vector<int> &result);
};

#endif
