#ifndef Ordinates_hh
#define Ordinates_hh

#include <string>
#include <vector>

#include "Gauss_Legendre.hh"

using std::string;
using std::vector;
    
class Ordinates
{
private:

    unsigned number_of_ordinates_;

    vector<double> ordinates_;
    vector<double> weights_;
    // vector<double> alpha_;
    // vector<double> alpha_half_;
        
public:
    
    Ordinates(unsigned number_of_ordinates);

    void check_class_invariants() const;
    void write_data(string folder) const;
    
    unsigned number_of_ordinates() const
    {
        return number_of_ordinates_;
    }

    double ordinates(unsigned const o) const
    {
        return ordinates_[o];
    }

    double weights(unsigned const o) const
    {
        return weights_[o];
    }

    // double alpha(unsigned o)
    // {
    //     return alpha_[o];
    // }

    // double alpha_half(unsigned o)
    // {
    //     return alpha_half_[o];
    // }
};

#endif
