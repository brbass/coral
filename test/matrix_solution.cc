#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "mpi.h"

#include "Check.hh"
#include "GSL_Linear_Algebra.hh"
#include "Matrix_Operations.hh"
#include "Random_Number_Generator.hh"
#include "Timer.hh"
#include "Trilinos_Linear_Algebra.hh"

using namespace std;

GSL_Linear_Algebra gsl_solver;
Matrix_Operations matrix_operator;
Trilinos_Linear_Algebra trilinos_solver;
Random_Number_Generator rng(0, 1);

void matrix_solution(vector<double> const &mat,
                     vector<double> const &lhs,
                     unsigned number_of_elements,
                     unsigned number_of_solves = 1,
                     bool print_results = false,
                     bool print_debug = false,
                     bool print_timing = true)
{
    Check(mat.size() == number_of_elements*number_of_elements, "mat size");
    Check(lhs.size() == number_of_elements, "lhs size");

    // calculate rhs
    vector<double> rhs(number_of_elements);
    
    matrix_operator.multiply(mat, 
                             lhs, 
                             rhs,
                             number_of_elements,
                             number_of_elements,
                             1);
    
    enum Methods
    {
        AMESOS,
        AZTEC,
        BOOST,
        EPETRA,
        GSL_LU,
        GSL_QR,
        GSL_QR_LS
    };
    
    vector<Methods> methods = {AMESOS, EPETRA, GSL_LU, GSL_QR, GSL_QR_LS};
    unsigned num_methods = methods.size();
    
    vector<string> method_description(num_methods);
    vector<double> time(num_methods, 0);
    vector<double> error(num_methods, 0);
    vector<vector<double>> lhs_store;
    
    // solve problems

    Timer timer;

    for (unsigned i=0; i<num_methods; ++i)
    {
        timer.start();
        
        vector<double> mat_temp(number_of_elements * number_of_elements);
        vector<double> rhs_temp(number_of_elements);
        vector<double> lhs_temp(number_of_elements, 0);
        
        for (unsigned j = 0; j < number_of_solves; ++j)
        {
            mat_temp = mat;
            rhs_temp = rhs;
            lhs_temp.assign(number_of_elements, 0);
            
            switch (methods[i])
            {
            case AMESOS:
                trilinos_solver.amesos_dense_solve(mat_temp, 
                                                   rhs_temp, 
                                                   lhs_temp,
                                                   number_of_elements);
                method_description[i] = "Amesos";
                break;

            case AZTEC:
                trilinos_solver.aztec_dense_solve(mat_temp, 
                                                  rhs_temp, 
                                                  lhs_temp,
                                                  number_of_elements);
                method_description[i] = "Aztec";
                break;

            case EPETRA:
                trilinos_solver.epetra_solve(mat_temp, 
                                             rhs_temp, 
                                             lhs_temp,
                                             number_of_elements);
                method_description[i] = "Epetra";
                break;

            case GSL_LU:
                gsl_solver.lu_solve(mat_temp, 
                                    rhs_temp, 
                                    lhs_temp,
                                    number_of_elements);
                method_description[i] = "GSL LU";
                break;
            case GSL_QR:
                gsl_solver.qr_solve(mat_temp, 
                                    rhs_temp, 
                                    lhs_temp,
                                    number_of_elements);
                method_description[i] = "GSL QR";
                break;
            case GSL_QR_LS:
                gsl_solver.qr_lssolve(mat_temp, 
                                      rhs_temp, 
                                      lhs_temp,
                                      number_of_elements,
                                      number_of_elements);
                method_description[i] = "GSL QR LS";
                break;
            }
        }
        timer.stop();
        time[i] = timer.time() / number_of_solves;

        double sum = 0;
        for (unsigned j=0; j<number_of_elements; ++j)
        {
            sum += pow(lhs_temp[j] - lhs[j], 2);
        }
        error[i] = sum / number_of_elements;
        
        if (print_results)
        {
            lhs_store.push_back(lhs_temp);
        }
    }

    // print results

    unsigned w = 16;
    
    if (print_timing)
    {
        cout << setw(w) << "Method";
        cout << setw(w) << "Timing";
        cout << setw(w) << "Mean Sq. Error";
        cout << endl;

        for (unsigned i=0; i<num_methods; ++i)
        {
            cout << setw(w) << method_description[i];
            cout << setw(w) << time[i];
            cout << setw(w) << error[i];
            cout << endl;
        }
        cout << endl;
    }

    if (print_results)
    {
        unsigned number_to_print = 10;
        unsigned print_every = ceil(1.*number_of_elements / number_to_print);
        
        cout << setw(w) << "cell";
        for (unsigned i=0; i<num_methods; ++i)
        {
            cout << setw(w) << method_description[i];
        }
        cout << endl;

        for (unsigned i=0; i<number_of_elements; ++i)
        {
            if (i % print_every == 0)
            {
                string cell = to_string(i) + " / " + to_string(number_of_elements);
                cout << setw(w) << cell;
                for (unsigned j=0; j<num_methods; ++j)
                {
                    cout << setw(w) << lhs_store[j][i];
                }
                cout << endl;
            }
        }
    }
    
    if (print_debug)
    {
        cerr << "debug printing not yet implemented" << endl;
    }
}

void random_matrix_solution(unsigned number_of_elements,
                            unsigned number_of_solves = 1,
                            bool print_results = false,
                            bool print_debug = false,
                            bool print_timing = true)
{
    // initialize data

    vector<double> mat(rng.random_double_vector(number_of_elements*number_of_elements));
    vector<double> lhs(rng.random_double_vector(number_of_elements));

    matrix_solution(mat, 
                    lhs, 
                    number_of_elements,
                    number_of_solves,
                    print_results,
                    print_debug,
                    print_timing);
}

void ill_matrix_solution(unsigned number_of_elements,
                         unsigned number_of_solves = 1,
                         bool print_results = false,
                         bool print_debug = false,
                         bool print_timing = true)
{
    vector<double> mat(number_of_elements*number_of_elements);
    vector<double> lhs(rng.random_double_vector(number_of_elements));
    
    // use Hilbert matrix for test

    for (unsigned i = 0; i < number_of_elements; ++i)
    {
        for (unsigned j = 0; j < number_of_elements; ++j)
        {
            mat[j + number_of_elements * i] = 1./(j+i+1); 
            // mat[j + number_of_elements * i] = 10000 + 1./(j+i+1); 
        }
    }
    
    matrix_solution(mat, 
                    lhs, 
                    number_of_elements,
                    number_of_solves,
                    print_results,
                    print_debug,
                    print_timing);
    
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    // parse input

    if (argc < 2)
    {
        cerr << "usage: matrix_solution [number_of_elements number_of_solves=1 print_results=false print_debug=false]" << endl;
        return 1;
    }
    
    unsigned number_of_elements = 2;
    unsigned number_of_solves = 1;
    bool print_results = false;
    bool print_debug = false;
    bool print_timing = true;

    for (int i=1; i<argc; ++i)
    {
        switch(i)
        {
        case 1:
            number_of_elements = atoi(argv[1]);
            break;
        case 2:
            number_of_solves = atoi(argv[2]);
            break;
        case 3:
            print_results = atoi(argv[3]);
            break;
        case 4:
            print_debug = atoi(argv[4]);
            break;
        }
    }

    cout << "random" << endl;
    random_matrix_solution(number_of_elements,
                           number_of_solves,
                           print_results,
                           print_debug,
                           print_timing);

    cout << "ill-conditioned" << endl;
    ill_matrix_solution(number_of_elements,
                        number_of_solves,
                        print_results,
                        print_debug,
                        print_timing);


    MPI_Finalize();
}
