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
    std::vector<size_t> strides;
    
    template <typename ValueType>
    class IteratorType : public std::iterator<std::forward_iterator_tag, ValueType>
    {
    private:
        ValueType *ptr;
        size_t step;
        
    public:
        IteratorType ()
            : ptr(NULL), step(0) {}
        IteratorType (ValueType *p, const size_t step = 1)
            : ptr(p), step(step) {}
        IteratorType (const IteratorType<ValueType> &other)
            : ptr(other.ptr), step(other.step) {}
        
        IteratorType<ValueType> & operator++ () { ptr += step; return *this; }
        IteratorType<ValueType> operator+ (ptrdiff_t n) { return IteratorType<ValueType>(ptr + n*step, step); }
        
        bool operator== (const IteratorType<ValueType> &other) { return (ptr==other.ptr && step==other.step); }
        bool operator!= (const IteratorType<ValueType> &other) { return (ptr!=other.ptr || step!=other.step); }
        bool operator> (const IteratorType<ValueType> &other) { return (ptr > other.ptr); }
        
        ValueType & operator* () { return *ptr; }
        
        bool valid () { return (ptr != NULL && step > 0); }
    };
    
    void calculateStrides ()
    {
        strides = std::vector<size_t>(nDims+1);
        strides[0] = 1;
        for (int i=0; i<nDims; i++)
            strides[i+1] = strides[i] * size_t(dims[i]);
    }

public:
    typedef IteratorType<const DataType> ConstIterator;
    typedef IteratorType<DataType> Iterator;
    typedef const DataType & ConstReference;
    typedef DataType & Reference;
    
    Array () { nDims = 0; }
    
    Array (const std::vector<int> &dims, const DataType &value)
        : dims(dims)
    {
        nDims = dims.size();
        calculateStrides();
        
        size_t length = 1;
        for (int i=0; i<nDims; i++)
            length *= dims[i];
        
        data = std::vector<DataType>(length, value);
    }
    
    Array (const std::vector<int> &dims, const std::vector<DataType> &data)
        : data(data), dims(dims)
    {
        nDims = dims.size();
        calculateStrides();
    }
    
    size_t size () const { return data.size(); }
    bool empty () const { return (data.size() == 0); }
    
    void fill (const DataType &value) { data.assign(data.size(), value); }
    
    ConstIterator begin () const { return ConstIterator(&data.front()); }
    Iterator begin () { return Iterator(&data.front()); }
    ConstIterator end () const { return ConstIterator(&data.back() + 1); }
    Iterator end () { return Iterator(&data.back() + 1); }
    
    ConstIterator beginLine (const std::vector<int> &begin, const int dim) const { return ConstIterator(&at(begin), strides[dim]); }
    Iterator beginLine (const std::vector<int> &begin, const int dim) { return Iterator(&at(begin), strides[dim]); }
    ConstIterator endLine (const std::vector<int> &begin, const int dim) const
    {
        const size_t offset = (dims[dim] - begin[dim]) * strides[dim];
        return ConstIterator(&at(begin) + offset, strides[dim]);
    }
    Iterator endLine (const std::vector<int> &begin, const int dim)
    {
        const size_t offset = (dims[dim] - begin[dim]) * strides[dim];
        return Iterator(&at(begin) + offset, strides[dim]);
    }
    
    ConstReference at (const size_t n) const { return data.at(n); }
    ConstReference at (const std::vector<int> &loc) const { return data.at(flattenIndex(loc)); }
    
    Reference at (const size_t n) { return data.at(n); }
    Reference at (const std::vector<int> &loc) { return data.at(flattenIndex(loc)); }
    
    ConstReference operator[] (const size_t n) const { return data[n]; }
    ConstReference operator[] (const std::vector<int> &loc) const { return data[flattenIndex(loc)]; }
    
    Reference operator[] (const size_t n) { return data[n]; }
    Reference operator[] (const std::vector<int> &loc) { return data[flattenIndex(loc)]; }
    
    const std::vector<DataType> & getData () const { return data; }
    const std::vector<int> & getDimensions () const { return dims; }
    int getDimensionality () const { return nDims; }
    
    Neighbourhood getNeighbourhood () const;
    Neighbourhood getNeighbourhood (const int width) const;
    Neighbourhood getNeighbourhood (const std::vector<int> &widths) const;
    
    void flattenIndex (const std::vector<int> &loc, size_t &result) const;
    void expandIndex (const size_t &loc, std::vector<int> &result) const;
    
    size_t flattenIndex (const std::vector<int> &loc) const
    {
        size_t n;
        flattenIndex(loc, n);
        return n;
    }
};

#endif
