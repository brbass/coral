#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "Array.hh"
#include "Check.hh"
#include "Random_Number_Generator.hh"
#include "Timer.hh"

using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cerr << "usage: array [first_dimension second_dimension ...]" << endl;
        return 1;
    }
    int number_of_dimensions = argc - 1;
    
    Random_Number_Generator rng(-100, 100);

    vector<int> size(number_of_dimensions);
    int total_size = 1;
    for (int i = 0; i < number_of_dimensions; ++i)
    {
        size[i] = atoi(argv[i+1]);
        total_size *= size[i];
    }
    vector<double> data(rng.random_double_vector(total_size));
    
    string description = "hi";

    Array<double> array(data,
                        size,
                        description);
    
    cout << "subscript check" << endl;
    unsigned w = 16;
    // cout << setw(w) << "original" << setw(w) << "check value";
    // for (int d = 0; d < number_of_dimensions; ++d)
    // {
    //     cout << setw(w) << "dimension " + to_string(d);
    // }
    // cout << endl;
    // for (int i = 0; i < total_size; ++i)
    // {
    //     vector<int> subscript = array.index_to_subscript(i);
    //     int new_index = array.subscript_to_index(subscript);
    //     cout << setw(w) << i << setw(w) << new_index;
    //     for (int d = 0; d < number_of_dimensions; ++d)
    //     {
    //         cout << setw(w) << subscript[d];
    //     }
    //     cout << endl;
    // }
    for (int i = 0; i < total_size; ++i)
    {
        Check(i == array.subscript_to_index(array.index_to_subscript(i)), "");
    }
    cout << "passed" << endl << endl;
    
    cout << "time to multiply two arrays" << endl;
    vector<vector<int> > indices;
    for (unsigned i = 0; i < total_size; ++i)
    {
        indices.push_back(array.index_to_subscript(i));
    }
    Check(indices.size() == total_size, "");

    Timer timer;

    unsigned num_iterations = ceil(1e6 / total_size);

    timer.start();
    for (unsigned t = 0; t < num_iterations; ++t)
    {
        for (unsigned i = 0; i < total_size; ++i)
        {
            data[i] = 0.99*data[i];
        }
    }
    timer.stop();
    cout << "vector: " << timer.time() << endl;

    timer.start();
    for (unsigned t = 0; t < num_iterations; ++t)
    {
        for (unsigned i = 0; i < total_size; ++i)
        {
            array(indices[i]) = 0.99*array(indices[i]);
        }
    }
    timer.stop();
    cout << "array, subscript: " << timer.time() << endl;

    timer.start();
    for (unsigned t = 0; t < num_iterations; ++t)
    {
        for (unsigned i = 0; i < total_size; ++i)
        {
            array(i) = 0.99*array(i);
        }
    }
    timer.stop();
    cout << "array, linear: " << timer.time() << endl;
}
