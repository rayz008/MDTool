#include <fstream>
#include <sstream>
#include <algorithm>

#include "settings.h"


Settings::Settings(const char* filename)
{
    ReadSettings(filename);
}

void Settings::ReadSettings(const char* filename)
{
    std::ifstream file(filename);
    if(file.is_open())
    {
        std::string line;
        while(getline(file, line))
        {
            std::istringstream iss;
            auto SkipArgs = [&]()
            {
                getline(file, line);
                iss.clear();
                iss.str(line);
                iss.ignore(256, ' ');
                iss.ignore(256, ' ');
            };
            if(line == "[Files]"){
                // ADD FUNCTION THAT DOES THIS
                SkipArgs();
                iss >> infile;
                SkipArgs();
                iss >> outfile;
                SkipArgs();
                iss >> increfile;
            }
            else if (line == "[System]"){
                SkipArgs();
                iss >> atom1;
                SkipArgs();
                iss >> atom2;
                getline(file, line);
                box = line;
                if(line.size() != 0){
                    NVT = true;
                }
            }
            else if (line == "[Settings]"){
                SkipArgs();
                iss >> r_min;
                SkipArgs();
                iss >> r_max;
                SkipArgs();
                iss >> bins;
                SkipArgs();
                iss >> increment;
            }
        }
        file.close();
    } else {
        throw std::logic_error("Unexpected file format for settings.");
    }


    return;
}
