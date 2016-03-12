#ifndef Linear_Algebra_hh
#define Linear_Algebra_hh

#include "GSL_Linear_Algebra.hh"
#include "Trilinos_Linear_Algebra.hh"

#include <vector>

using std::vector;

class Linear_Algebra
{
private:
    
    unsigned size_;
    
    GSL_Linear_Algebra gsl_solver;
    Trilinos_Linear_Algebra trilinos_solver;
    
public:
    
    Linear_Algebra(unsigned size);
    
    void solve(vector<double> &a_data, 
               vector<double> &b_data,
               vector<double> &x_data);
    
    unsigned size()
    {
        return size_;
    }
};

#endif
