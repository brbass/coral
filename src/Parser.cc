#include "Parser.hh"

#include <cstdlib>
#include <cstring>
#include <string>

using namespace std;

Parser::
Parser(string &folder)
{
    folder_ = folder;
        
    system(("mkdir -p " + folder_).c_str());
}
    
void Parser::
set_folder(string folder)
{
    folder_ = folder;
        
    system(("mkdir -p " + folder_).c_str());
}

string Parser::
get_filepath(string const &data_filename)
{
    return folder_ + "/" + data_filename;
}
