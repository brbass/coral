#include "Linear_Algebra.hh"

#include <iostream>
#include <vector>

#include "Check.hh"

using namespace std;

Linear_Algebra::
Linear_Algebra(unsigned size):
    size_(size)
{
}

// solves linear system Ax=b
void Linear_Algebra::
solve(vector<double> &a_data,
      vector<double> &b_data,
      vector<double> &x_data)
{
    Check(a_data.size() == size()*size(), "A size");
    Check(b_data.size() == size(), "b size");
    Check(x_data.size() == size(), "x size");

    if (size() < 280) // fastest for small problems
    {
        gsl_solver.lu_solve(a_data,
                            b_data,
                            x_data,
                            size());
    }
    else // fastest for large problems
    {
        trilinos_solver.epetra_solve(a_data,
                                     b_data,
                                     x_data,
                                     size());
    }
}

