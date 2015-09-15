//
//  IOFunctions.h
//  Convolutions
//
//  Created by Stephan Meesters on 18/08/15.
//
//

#ifndef Convolutions_IOFunctions_h
#define Convolutions_IOFunctions_h

#include <fstream>
#include <string>

inline bool file_exists (const std::string& name)
{
    std::ifstream f(name.c_str());
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }
}

#endif
