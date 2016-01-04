#ifndef _ARRAY_H_
#define _ARRAY_H_

#include <RcppEigen.h>

struct Neighbourhood
{
    size_t size;
    std::vector<int> widths;
    Eigen::ArrayXXi locs;
    std::vector<ptrdiff_t> offsets;
};

template <typename DataType> class Array
{
protected:
    std::vector<DataType> data;
    std::vector<int> dims;
    int nDims;

public:
    typedef typename std::vector<DataType>::const_iterator const_iterator;
    typedef typename std::vector<DataType>::iterator iterator;
    typedef typename std::vector<DataType>::const_reference const_reference;
    typedef typename std::vector<DataType>::reference reference;
    
    Array () { nDims = 0; }
    
    Array (const std::vector<int> &dims, const DataType &value)
        : dims(dims)
    {
        nDims = dims.size();
        
        size_t length = 1;
        for (int i=0; i<nDims; i++)
            length *= dims[i];
        
        data = std::vector<DataType>(length, value);
    }
    
    Array (const std::vector<int> &dims, const std::vector<DataType> &data)
        : data(data), dims(dims)
    {
        nDims = dims.size();
    }
    
    size_t size () const { return data.size(); }
    bool empty () const { return (data.size() == 0); }
    
    void fill (const DataType &value) { data.assign(data.size(), value); }
    
    const_iterator begin () const { return data.begin(); }
    iterator begin () { return data.begin(); }
    const_iterator end () const { return data.end(); }
    iterator end () { return data.end(); }
    
    const_reference at (const size_t n) const { return data.at(n); }
    const_reference at (const std::vector<int> &loc) const
    {
        size_t n;
        flattenIndex(loc, n);
        return data.at(n);
    }
    
    reference operator[] (const size_t n) { return data[n]; }
    reference operator[] (const std::vector<int> &loc)
    {
        size_t n;
        flattenIndex(loc, n);
        return data[n];
    }
    
    const_reference operator[] (const size_t n) const { return data[n]; }
    const_reference operator[] (const std::vector<int> &loc) const
    {
        size_t n;
        flattenIndex(loc, n);
        return data[n];
    }
    
    const std::vector<DataType> & getData () const { return data; }
    const std::vector<int> & getDimensions () const { return dims; }
    int getDimensionality () const { return nDims; }
    
    Neighbourhood getNeighbourhood () const;
    Neighbourhood getNeighbourhood (const int width) const;
    Neighbourhood getNeighbourhood (const std::vector<int> &widths) const;
    
    void flattenIndex (const std::vector<int> &loc, size_t &result) const;
    void expandIndex (const size_t &loc, std::vector<int> &result) const;
};

#endif
