#ifndef Parser_hh
#define Parser_hh

#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "Check.hh"

using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;

class Parser
{
private:

    string folder_;

public:

    Parser(string &folder);

    template<class T> void parse_data(T &data, string data_filename)
    {
        string data_path = get_filepath(data_filename);
        
        ifstream data_file(data_path.c_str());
        
        if (data_file.is_open())
        {
            data_file >> data;
        }
        else
        {
            Check(false, "unable to load " + data_filename);
        }
    }

    template<class T> void parse_data(vector<T> &data, string data_filename)
    {
        data.resize(0);
        
        string data_path = get_filepath(data_filename);
        
        ifstream data_file(data_path.c_str());
        
        if (data_file.is_open())
        {
            unsigned temp;
            
            while (data_file >> temp)
            {
                data.push_back(temp);
            }
        }
        else
        {
            Check(false, "unable to load " + data_filename);
        }
    }

    template<class T> void write_data(T const &data, string data_filename)
    {
        string data_path = get_filepath(data_filename);
        
        ofstream data_file(data_path.c_str());
    
        if (data_file.is_open())
        {
            data_file << data;
        }
        else
        {
            Check(false, "unable to save to " + data_filename);
        }
    }

    template<class T> void write_data(vector<T> const &data, string data_filename)
    {
        string data_path = get_filepath(data_filename);
    
        ofstream data_file(data_path.c_str());
        
        if (data_file.is_open())
        {
            for (unsigned i = 0; i < data.size(); ++i)
            {
                data_file << data[i] << std::endl;
            }
        }
        else
        {
            Check(false, "unable to save to " + data_filename);
        }
    }

    void set_folder(string folder);

    string get_filepath(string const &value_filename);
};

#endif
