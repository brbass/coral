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
    int size1 = 20;
    int size2 = 19;
    int size3 = 18;
    
    int number_of_dimensions = 3;
    
    Random_Number_Generator rng(-100, 100);
    
    vector<int> size = {size1, size2, size3};
    int total_size = size1 * size2 * size3;

    vector<double> data(rng.random_double_vector(total_size));
    
    string description = "hi";
    
    Array<double> array(data,
                        size,
                        description);
    
    // check subscripts

    cout << "subscript check" << endl;
    unsigned w = 16;
    for (int i = 0; i < total_size; ++i)
    {
        Check(i == array.subscript_to_index(array.index_to_subscript(i)), "");
    }
    cout << "passed" << endl << endl;

    // multiply arrays
    
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
    cout << "vector, standard: " << timer.time() << endl;

    timer.start();
    for (unsigned t = 0; t < num_iterations; ++t)
    {
        for (unsigned i = 0; i < size1; ++i)
        {
            for (unsigned j = 0; j < size2; ++j)
            {
                for (unsigned k = 0; k < size3; ++k)
                {
                    unsigned l = k + size3 * (j + size2 * i);
                    
                    data[l] = 0.99*data[l];
                }
            }
        }
    }
    timer.stop();
    cout << "vector, linear: " << timer.time() << endl;

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

    timer.start();
    for (unsigned t = 0; t < num_iterations; ++t)
    {
        for (unsigned i = 0; i < size1; ++i)
        {
            for (unsigned j = 0; j < size2; ++j)
            {
                for (unsigned k = 0; k < size3; ++k)
                {
                    array(i, j, k) = 0.99*array(i, j, k);
                }
            }
        }
    }
    timer.stop();
    cout << "array, multiple: " << timer.time() << endl;
    cout << endl;
}
