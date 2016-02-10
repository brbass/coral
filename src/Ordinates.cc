#include "Ordinates.hh"

#include <string>

#include "Check.hh"
#include "Gauss_Legendre.hh"
#include "Parser.hh"

using namespace std;

Ordinates::
Ordinates(unsigned number_of_ordinates):
    number_of_ordinates_(number_of_ordinates)
{
    gauss_legendre_vec(number_of_ordinates_, ordinates_, weights_);

    // alpha_.resize(number_of_ordinates_, 0);
    // alpha_half_.resize(number_of_ordinates_, 0);
        
    // alpha_half_[0] = -ordinates_[0] * weights_[0];
    // alpha_[0] = alpha_half_[0];
        
    // for (unsigned o = 1; o < number_of_ordinates_; ++o)
    // {
    //     alpha_half_[o] = alpha_half_[o-1] - ordinates_[o] * weights_[o];
    //     alpha_[o] = alpha_half_[o] + alpha_half_[o-1];
    // }

    check_class_invariants();
}

void Ordinates::
check_class_invariants() const
{
    Check(ordinates_.size() == number_of_ordinates_, "ordinates size");
    Check(weights_.size() == number_of_ordinates_, "weights size");
}

void Ordinates::
write_data(string folder) const
{
    Parser parser(folder);
    
    parser.write_data(number_of_ordinates_, "ordinates_number_of_ordinates");
    parser.write_data(ordinates_, "ordinates_ordinates");
    parser.write_data(weights_, "ordinates_weights");
}
