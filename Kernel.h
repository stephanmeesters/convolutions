//
//  Kernel.h
//  Convolutions
//
//  Created by Stephan Meesters on 14/08/15.
//
//

#ifndef __Convolutions__Kernel__
#define __Convolutions__Kernel__

#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "MathFunctions.h"
#include "IOFunctions.h"

#define PI 3.14159265358979323

class Kernel
{
private:
    
    double D33;
    double D44;
    double t;
    
    void DefineOrientations();
    
    double k2(double* x, double* y, double* r, double* v);
    double kernel(double* c);
    double* coordinateMap(double x, double y, double z, double beta, double gamma);

    double GetMaximumValue();
    void EstimateKernelSize();
    
public:
    
    Kernel(double D33, double D44, double t);
    ~Kernel();
    
    double* Create(std::string fileName);
    
    double* orientations; // 162 * x,y,z
    
    double* lut;
    
    int xsize;
    int ysize;
    int zsize;
};

#endif /* defined(__Convolutions__Kernel__) */
