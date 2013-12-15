#ifndef _ARRAY_HPP_
#define _ARRAY_HPP_

typedef std::vector<double>::size_type size_type;

struct Neighbourhood
{
    int width;
    int size;
    arma::imat locs;
    std::vector<size_type> offsets;
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
    
    // double & at (size_type n) { data.at(n); }
    
    const double & at (size_type n) { return data.at(n); }
    
    size_type size () { return data.size(); }
    
    Neighbourhood getNeighbourhood (const int width);
    
    void flattenIndex (const std::vector<int> &loc, size_type &result);
    
    void expandIndex (const size_type &loc, std::vector<int> &result);
};

#endif
