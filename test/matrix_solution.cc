#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "mpi.h"

#include "Check.hh"
#include "Linear_Algebra.hh"
#include "Random_Number_Generator.hh"
#include "Timer.hh"

using namespace std;

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    // parse input

    if (argc < 2)
    {
        cerr << "usage: matrix_solution [number_of_elements print_results=false print_debug=false]" << endl;
        return 1;
    }
    
    unsigned number_of_elements = 2;
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
            print_results = atoi(argv[2]);
            break;
        case 3:
            print_debug = atoi(argv[3]);
            break;
        }
    }

    // initialize data

    Random_Number_Generator random(0, 1);
    
    vector<double> mat(random.random_double_vector(number_of_elements*number_of_elements));
    vector<double> lhs(random.random_double_vector(number_of_elements));
    
    Linear_Algebra solver(number_of_elements);
    vector<double> rhs(number_of_elements, 0);
    
    solver.multiply(mat, lhs, rhs);
    
    enum Methods
    {
        AMESOS,
        AZTEC,
        BOOST,
        EPETRA,
        GSL
    };
    
    vector<Methods> methods = {AMESOS, AZTEC, EPETRA, GSL};
    unsigned num_methods = methods.size();
    
    vector<string> method_description(num_methods);
    vector<double> time(num_methods, 0);
    vector<double> error(num_methods, 0);
    vector<vector<double>> lhs_store;
    
    // solve problems

    Timer timer;

    for (unsigned i=0; i<num_methods; ++i)
    {
        vector<double> mat_temp(mat);
        vector<double> rhs_temp(rhs);
        vector<double> lhs_temp(number_of_elements, 0);
        
        timer.start();

        switch (methods[i])
        {
        case AMESOS:
            solver.amesos_solve(mat_temp, rhs_temp, lhs_temp);
            method_description[i] = "Amesos";
            break;

        case AZTEC:
            solver.aztec_solve(mat_temp, rhs_temp, lhs_temp);
            method_description[i] = "Aztec";
            break;

        case BOOST:
            solver.boost_solve(mat_temp, rhs_temp, lhs_temp);
            method_description[i] = "Boost";
            break;

        case EPETRA:
            solver.epetra_solve(mat_temp, rhs_temp, lhs_temp);
            method_description[i] = "Epetra";
            break;

        case GSL:
            solver.gsl_solve(mat_temp, rhs_temp, lhs_temp);
            method_description[i] = "GSL";
            break;
        }
        
        timer.stop();
        time[i] = timer.time();

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

    MPI_Finalize();
}
