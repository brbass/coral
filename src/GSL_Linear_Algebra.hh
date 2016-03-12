#ifndef GSL_Linear_Algebra_hh
#define GSL_Linear_Algebra_hh

#include <vector>

using std::vector;

class GSL_Linear_Algebra
{
private:
    
public:

    GSL_Linear_Algebra();
    
    void lu_solve(vector<double> &a_data,
                  vector<double> &b_data,
                  vector<double> &x_data,
                  unsigned number_of_rows);
  
    void qr_solve(vector<double> &a_data,
                  vector<double> &b_data,
                  vector<double> &x_data,
                  unsigned number_of_rows);

    void qr_lssolve(vector<double> &a_data,
                    vector<double> &b_data,
                    vector<double> &x_data,
                    unsigned number_of_rows,
                    unsigned number_of_columns);

};

#endif
