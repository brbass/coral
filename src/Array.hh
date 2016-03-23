#ifndef Array_hh
#define Array_hh

#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "Check.hh"

using std::ostream;
using std::string;
using std::vector;

template <class T> class Array
{
private:

    int number_of_dimensions_;
    int total_size_;
    string description_;
    vector<int> size_;
    vector<T> data_;

public:

    Array(string const description = "");
    Array(vector<int> const size,
          string const description = "");
    Array(vector<T> const &data, 
          vector<int> const size,
          string const description = "");
    
    void resize(vector<int> const size);
    void assign(vector<int> const size, T const value);
    void set_description(string const description);

    int get_total_size(vector<int> const size);

    T &operator()(int const index);
    T &operator()(vector<int> const subscript);
    template<class U> friend ostream &operator<<(ostream &out, Array<U> &array);
    
    int subscript_to_index(vector<int> const subscript);
    vector<int> index_to_subscript(int const index);
        
    int size() const
    {
        return size_;
    }
    int dimension(int const dimension)
    {
        Check(dimension < number_of_dimensions_, "number of dimensions");
        
        return size_[dimension];
    }
    int number_of_dimensions()
    {
        return number_of_dimensions_;
    }
    string description() const
    {
        return description_;
    }
    vector<int> dimensions() const
    {
        return size_;
    }
    vector<T> data() const
    {
        return data_;
    }
};

#endif
