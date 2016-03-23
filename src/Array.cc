#include "Array.hh"

#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using namespace std;

template<class T> Array<T>::
Array(string const description)
{
    vector<int> size(0);
    resize(size);
    description_ = description;
}

template<class T> Array<T>::
Array(vector<int> const size,
      string const description)
{
    resize(size);
    description_ = description;
}

template<class T> Array<T>::
Array(vector<T> const &data,
      vector<int> const size,
      string const description)
{
    resize(size);
    description_ = description;
    data_ = data;
}

template<class T> void Array<T>::
resize(vector<int> const size)
{
    number_of_dimensions_ = size.size();
    total_size_ = get_total_size(size);
    size_ = size;
    data_.resize(total_size_);
}

template<class T> void Array<T>::
assign(vector<int> const size, T const value)
{
    number_of_dimensions_ = size.size();
    total_size_ = get_total_size(size);
    size_ = size;
    data_.assign(total_size_, value);
}

template<class T> void Array<T>::
set_description(string const description)
{
    description_ = description;
}

template<class T> int Array<T>::
get_total_size(vector<int> const size)
{
    if (size.size() == 0)
    {
        return 0;
    }
    else
    {
        int product = 1;
        for (int i = 0; i < size.size(); ++i)
        {
            product *= size[i];
        }
        
        return product;
    }
}

template<class T> int Array<T>::
subscript_to_index(vector<int> const subscript)
{
    Check(subscript.size() == number_of_dimensions_, "subscript size");
    
    int sum = 0;
    
    for (int i = 0; i < number_of_dimensions_; ++i)
    {
        Check(subscript[i] < size_[i], "subscript size");
        
        sum = subscript[i] + size_[i] * sum;
    }
 
    return sum;
}

template<class T> vector<int> Array<T>::
index_to_subscript(int const index)
{
    Check(index < total_size_, "index");

    vector<int> subscript(number_of_dimensions_);

    int product = 1;
    for (int i = 1; i < number_of_dimensions_; ++i)
    {
        product *= size_[i];
    }
    
    int sum = index;
    for (int i = 0; i < number_of_dimensions_ - 1; ++i)
    {
        subscript[i] = floor(static_cast<double>(sum) / product);
        
        sum -= product * subscript[i];
        product /= size_[i];
    }
    
    subscript[number_of_dimensions_ - 1] = floor(static_cast<double>(sum) / product);
    
    return subscript;
}

template<class T> T &Array<T>::
operator()(int const index)
{
    Check(index < total_size_, "index");
    
    return data_[index];
}

template<class T> T &Array<T>::
operator()(vector<int> const subscript)
{
    int index = subscript_to_index(subscript);

    return data_[index];
}

template <class U>
ostream &operator<<(ostream &out, Array<U> &array)
{
    out << array.description_ << endl;
    
    for (int i = 0; i < array.total_size_; ++i)
    {
        vector<int> subscript = array.index_to_subscript(i);
        
        for (unsigned j = 0; j < array.number_of_dimensions_; ++j)
        {
            out << subscript[j] << "\t";
        }
        out << array.data_[i] << endl;
    }

    return out;
}
