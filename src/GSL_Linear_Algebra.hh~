#ifndef Linear_Algebra_hh
#define Linear_Algebra_hh

#include <vector>

using std::vector;

class Linear_Algebra
{
private:
    unsigned size_;

public:
    Linear_Algebra(unsigned size);
    
    void solve(vector<double> &a_data, 
               vector<double> &b_data,
               vector<double> &x_data);

    void gsl_solve(vector<double> &a_data, 
                   vector<double> &b_data,
                   vector<double> &x_data);

    void boost_solve(vector<double> &a_data, 
                     vector<double> &b_data,
                     vector<double> &x_data);

    void epetra_solve(vector<double> &a_data,
                      vector<double> &b_data,
                      vector<double> &x_data);

    void amesos_solve(vector<double> &a_data,
                      vector<double> &b_data,
                      vector<double> &x_data);

    void aztec_solve(vector<double> &a_data,
                     vector<double> &b_data,
                     vector<double> &x_data);
    
    void dot(vector<double> const &a_data,
             vector<double> const &b_data,
             double &x_data);
    void cross(vector<double> const &a_data,
               vector<double> const &b_data,
               vector<double> &x_data);

    void multiply(vector<double> const &a_data,
                  vector<double> const &b_data,
                  vector<double> &x_data);
    void multiply_vector_vector(vector<double> const &a_data,
                                vector<double> const &b_data,
                                vector<double> &x_data);
    void multiply_matrix_vector(vector<double> const &a_data,
                                vector<double> const &b_data,
                                vector<double> &x_data);
    void multiply_matrix_matrix(vector<double> const &a_data,
                                vector<double> const &b_data,
                                vector<double> &x_data);
    
    unsigned size()
    {
        return size_;
    }
};

#endif
