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

    Array<double> array1(description);
    Array<double> array2(size,
                         description);
    Array<double> array3(data,
                         size,
                         description);
}
