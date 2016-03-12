#include "Matrix_Operations.hh"

#include <iostream>
#include <vector>

#include "Check.hh"

using namespace std;

Matrix_Operations::
Matrix_Operations()
{
}

// dot product a . b
void Matrix_Operations::
dot(vector<double> const &a_data,
    vector<double> const &b_data,
    double &x_data,
    unsigned number_of_elements)
{
    Check(a_data.size() == number_of_elements, "a size");
    Check(b_data.size() == number_of_elements, "b size");
    
    double sum = 0;

    for (unsigned i=0; i<number_of_elements; ++i)
    {
        sum += a_data[i]*b_data[i];
    }

    x_data = sum;
}


// cross product A x B
void Matrix_Operations::
cross(vector<double> const &a_data,
      vector<double> const &b_data,
      vector<double> &x_data,
      unsigned number_of_elements)
{
    Check(number_of_elements == 3, "size must be 3 for cross product");
    Check(a_data.size() == number_of_elements, "a size");
    Check(b_data.size() == number_of_elements, "b size");
    Check(x_data.size() == number_of_elements, "x size");
    
    x_data[0] = a_data[1]*b_data[2] - a_data[2]*b_data[1];
    x_data[1] = a_data[2]*b_data[0] - a_data[0]*b_data[2];
    x_data[2] = a_data[0]*b_data[1] - a_data[1]*b_data[0];
}

// multiplies matrices A (n x m) and B (m x p)
void Matrix_Operations::
multiply(vector<double> const &a_data,
         vector<double> const &b_data,
         vector<double> &x_data,
         unsigned n,
         unsigned m,
         unsigned p)
{
    Check(a_data.size() == n * m, "a size");
    Check(b_data.size() == m * p, "b size");
    Check(x_data.size() == n * p, "x size");
    
    for (unsigned i=0; i<n; ++i)
    {
        for (unsigned j=0; j<p; ++j)
        {
            double sum = 0;
            
            for (unsigned k=0; k<m; ++k)
            {
                unsigned index1 = k + m*i;
                unsigned index2 = j + p*k;
                
                sum += a_data[index1]*b_data[index2];
            }
            
            unsigned index = j + p*i;

            x_data[index] = sum;
        }
    }
}

