//
//  MathFunctions.h
//  Convolutions
//
//  Created by Stephan Meesters on 14/08/15.
//
//

#ifndef Convolutions_MathFunctions_h
#define Convolutions_MathFunctions_h

#include <math.h>
#include <complex>

inline double cot(double i) { return(1 / tan(i)); }
inline double acosT(double i) { return(1.57079632679 - i - i*i*i/6 - 3*i*i*i*i*i/40); } // O(i)^6
inline double asinT(double i) { return(i + i*i*i/6 + 3*i*i*i*i*i/40); } // O(i)^6

inline double* Difference(double* vec, double* vec2)
{
    double* dp = (double*) malloc(3*sizeof(double));
    dp[0] = vec[0] - vec2[0];
    dp[1] = vec[1] - vec2[1];
    dp[2] = vec[2] - vec2[2];
    return dp;
}

inline double* Flip(double* vec)
{
    double* dp = (double*) malloc(3*sizeof(double));
    dp[0] = -vec[0];
    dp[1] = -vec[1];
    dp[2] = -vec[2];
    return dp;
}

inline double* HalvedDifference(double* vec, double* vec2)
{
    double* dp = (double*) malloc(3*sizeof(double));
    dp[0] = (vec[0] - vec2[0])/2.0;
    dp[1] = (vec[1] - vec2[1])/2.0;
    dp[2] = (vec[2] - vec2[2])/2.0;
    return dp;
}

inline double* Normalize(double* vec)
{
    double length = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    double* nvec = (double*) malloc(3*sizeof(double));
    nvec[0] = vec[0]/length;
    nvec[1] = vec[1]/length;
    nvec[2] = vec[2]/length;
    return nvec;
}

inline double* Cross(double* vec, double* vec2)
{
    double* c = (double*) malloc(3*sizeof(double));
    c[0] = vec[1]*vec2[2] - vec[2]*vec2[1];
    c[1] = vec[2]*vec2[0] - vec[0]*vec2[2];
    c[2] = vec[0]*vec2[1] - vec[1]*vec2[0];
    return c;
}

inline double Norm(double* vec)
{
    double output = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
    free(vec);
    return output;
}

//void PrintVector3(double* vec)
//{
//    if(vec == NULL)
//        return;
//    printf("%f, %f, %f \n", vec[0], vec[1], vec[2]);
//}
//
//void PrintVector2(double* vec)
//{
//    if(vec == NULL)
//        return;
//    printf("{%f, %f} \n", vec[0], vec[1]);
//}
//
//void PrintMatrix3x3(double* mat)
//{
//    printf("{%f, %f, %f \n", mat[0], mat[1], mat[2]);
//    printf(" %f, %f, %f \n", mat[3], mat[4], mat[5]);
//    printf(" %f, %f, %f} \n", mat[6], mat[7], mat[8]);
//}


inline double* EulerAngles(double* input)
{
    
    
    double x = input[0];
    double y = input[1];
    double z = input[2];
    double* output = (double*)malloc(2*sizeof(double));
    
    
    // based on EulerAngles from HARDIAlgorithms
//    if(x*x < 1e-6 && y*y < 1e-6 && (z-1)*(z-1) < 1e-6)
//    {
//        output[0] = 0;
//        output[1] = 0;
//    }
//    else if(x*x < 1e-6 && y*y < 1e-6 && (z+1)*(z+1) < 1e-6)
//    {
//        output[0] = -3.14159265358;
//        output[1] = 0;
//    }
//    else if((x - 1)*(x - 1) < 1e-10)
//    {
//        output[0] = 1.57079632679;
//        output[1] = 0;
//    }
//    else if((x + 1)*(x + 1) < 1e-10)
//    {
//        output[0] = -1.57079632679;
//        output[1] = 0;
//    }
//    else
//    {
//        output[0] = acos(sqrt(y*y + z*z)*(z>=0?1:-1))*(x>0?1:-1);
//        output[1] = -asin(y/sqrt(y*y + z*z)*(z>=0?1:-1));
//    }
    
    
    //EulerAnglesStandard
    if(x*x < 10e-6 && y*y < 10e-6 && (z-1)*(z-1) < 10e-6) // handle the case (0,0,1)
    {
        output[0] = 0;
        output[1] = 0;
    }
    else if(x*x < 10e-6 && y*y < 10e-6 && (z+1)*(z+1) < 10e-6) // handle the case (0,0,-1)
    {
        output[0] = M_PI;
        output[1] = 0;
    }
    else                                                    // all other cases
    {
        output[0] = acos(z);
        std::complex<double> complexxy (x,y);
        output[1] = std::arg(complexxy);
    }
    
    return output;
}

inline double* R(double* input)
{
    double beta = input[0];
    double gamma = input[1];
    double* output = (double*)malloc(9*sizeof(double));
    
    double cb = cos(beta);
    double sb = sin(beta);
    double cg = cos(gamma);
    double sg = sin(gamma);
    
    // old R
    //output[0] = cb;
    //output[1] = 0;
    //output[2] = sb;
    //output[3] = sb*sg;
    //output[4] = cg;
    //output[5] = -cb*sg;
    //output[6] = -cg*sb;
    //output[7] = sg;
    //output[8] = cb*cg;

    // Rstandard
    output[0] = cb*cg;
    output[1] = -sg;
    output[2] = cg*sb;
    output[3] = cb*sg;
    output[4] = cg;
    output[5] = sb*sg;
    output[6] = -sb;
    output[7] = 0;
    output[8] = cb;

    
    free(input);
    
    return output;
}

inline double* Transpose3x3(double *src)
{
    double* dst = (double*)malloc(sizeof(double)*9);
    for(int n = 0; n<9; n++)
    {
        int i = n/3;
        int j = n%3;
        dst[n] = src[3*j + i];
    }
    
    free(src);
    
    return dst;
}

inline double* Multiply(double* mat, double* vec)
{
    double* dst = (double*)malloc(sizeof(double)*3);
    dst[0] = mat[0]*vec[0] + mat[1]*vec[1]+mat[2]*vec[2];
    dst[1] = mat[3]*vec[0] + mat[4]*vec[1]+mat[5]*vec[2];
    dst[2] = mat[6]*vec[0] + mat[7]*vec[1]+mat[8]*vec[2];
    
    //free(mat);
    
    return dst;
}

inline double* Subtract3(double* vec, double* vec2)
{
    double* dst = (double*)malloc(sizeof(double)*3);
    dst[0] = vec[0]-vec2[0];
    dst[1] = vec[1]-vec2[1];
    dst[2] = vec[2]-vec2[2];
    return dst;
}

#endif
