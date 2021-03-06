//
//  Kernel.cpp
//  Convolutions
//
//  Created by Stephan Meesters on 14/08/15.
//
//

#include "Kernel.h"
#include <stdlib.h>
#include <cstring>

Kernel::Kernel(double D33, double D44, double t)
{
    // save kernel parameters
    this->D33 = D33;
    this->D44 = D44;
    this->t = t;
    
    // define orientations from tesselation
    DefineOrientations();
    
    // check if the kernel was already computed before with these settings
    char kernelName[50];
    sprintf(kernelName, "kernel_d33@%4.2f_d44@%4.2f_t@%4.2f.dat",D33,D44,t);
    bool kernelFileExists = file_exists(kernelName);
    
    if(kernelFileExists)
    {
        printf("Loading existing kernel lookup table\n");
        
        std::ifstream file(kernelName, std::ifstream::in | std::ifstream::binary);
        if(!file.is_open())
        {
            printf("File:%s not found.",kernelName);
            return;
        }
        
        // read data
        int numOrientations = 0, N = 0, M = 0;
        file.read(reinterpret_cast<char*>(&numOrientations),sizeof(int));
        file.read(reinterpret_cast<char*>(&N),sizeof(int));
        file.read(reinterpret_cast<char*>(&M),sizeof(int));
        
        if(numOrientations == 0 || N == 0 || M == 0)
        {
            printf("lookup table corrupted!\n");
            return;
        }
        printf("the kernel is %dx%dx%d\n",M,M,N);
        this->xsize = M;
        this->ysize = M;
        this->zsize = N;
        
        double* lookupTable = (double*)malloc(numOrientations*numOrientations*N*M*M*sizeof(double));
        file.read(reinterpret_cast<char*>(lookupTable),numOrientations*numOrientations*N*M*M*sizeof(double));
        file.close();
        
        this->lut = lookupTable;
    }
    
    // if not, create
    else
    {
        printf("Creating new kernel lookup table\n");
        this->lut = Create(kernelName);
    }
}

Kernel::~Kernel()
{
    free(orientations);
    free(this->lut);
}

// --------------------------- //

double* Kernel::Create(std::string fileName)
{
    double kernelMax = GetMaximumValue();
    EstimateKernelSize();

    int M = this->xsize;
    int N = this->zsize;
int numOrientations = 162;
    
    int minx = -(M-1)/2;
    int maxx = (M-1)/2+1;
    int miny = minx;
    int maxy = maxx;
    int minz = -(N-1)/2;
    int maxz = (N-1)/2+1;
    
    
    
    // loop over all orientations r which creates a lookup table of 486*N*M*M scalars
    double* lookupTable = (double*)malloc(numOrientations*numOrientations*N*M*M*sizeof(double));
    #pragma omp parallel for
    for (int angv = 0; angv<numOrientations; angv++)
    {
	double x[3] = {0,0,0};
        double y[3] = {0,0,0};
        double r[3] = {0,0,1};
        double v[3] = {0,0,1};

        v[0] = orientations[angv*3 + 0];
        v[1] = orientations[angv*3 + 1];
        v[2] = orientations[angv*3 + 2];
        
        //printf("v: (%f,%f,%f)\n",v[0],v[1],v[2]);
        
        for(int angr = 0; angr<numOrientations; angr++)
        {
            r[0] = orientations[angr*3 + 0];
            r[1] = orientations[angr*3 + 1];
            r[2] = orientations[angr*3 + 2];
            
            //printf("r: (%f,%f,%f)\n",r[0],r[1],r[2]);
            
            for (int xp = minx; xp<maxx; xp++)
            {
                for(int yp = miny; yp<maxy; yp++)
                {
                    for(int zp = minz; zp<maxz; zp++)
                    {
                        x[0] = xp;
                        x[1] = yp;
                        x[2] = zp;
                        
                        lookupTable[(xp-minx) + (yp-miny)*M + (zp-minz)*M*M + angr*M*M*N + angv*numOrientations*M*M*N] = k2(x,y,r,v)/kernelMax;
                        
                        //printf("y: (%d,%d,%d): %.10f at POS:%d\n",(int)x[0],(int)x[1],(int)x[2],k2(x,y,r,v)/kernelMax,(xp-minx) + (yp-miny)*M + (zp-minz)*M*M + angr*M*M*N + angv*numOrientations*M*M*N);
                        
                        //if(isnan(lookupTable[angv * numOrientations*M*M*N + angr*M*M*N + xp*M*N + yp*N + zp]))
                        //    printf("(%d:%f-%f-%f,%d:%f-%f-%f,%d,%d,%d) %f\n",angv,v[0],v[1],v[2],angr,r[0],r[1],r[2],xp,yp,zp,k2(x,y,r,v));
                        
                        //printf("%f\n",lookupTable[angv * numOrientations*M*M*N + angr*M*M*N + xp*M*N + yp*N + zp]);
                    }
                }
            }
        }
    }
    printf("lookup table complete.");
    
    // export the lookup table
    std::ofstream myfile;
    myfile.open (fileName.c_str(), std::ios::out | std::ios::binary);
    if ( !myfile.is_open() )
        return NULL;
    myfile.write(reinterpret_cast<const char*>(&numOrientations), std::streamsize(sizeof(int))); // num orientations
    myfile.write(reinterpret_cast<const char*>(&N), std::streamsize(sizeof(int)));
    myfile.write(reinterpret_cast<const char*>(&M), std::streamsize(sizeof(int)));
    myfile.write(reinterpret_cast<const char*>(lookupTable), std::streamsize(numOrientations*numOrientations*N*M*M*sizeof(double)));
    myfile.close();
    
    return lookupTable;
    
//    
//    double x[3] = {0,0,0};
//    double y[3] = {0,0,0};
//    double r[3] = {0,0,1};
//    double* v = &orientations[50];
//    
//    double test = k2(x,y,r,v);
//    printf("kernel output: %f \n",test);
//    
//    double* test2 = coordinateMap(0, 0, 0, 0.5, 0.22);
//    for(int i = 0; i<6; i++)
//        printf("coords output: %f \n",test2[i]);
//    
//    double test3 = kernel(test2);
//    printf("kernel output: %f \n",test3);
//    
//    double vv[3] = {0.2,0.4,0.6};
//    double* test4 = R(EulerAngles(vv));
//    for(int i = 0; i<9; i++)
//        printf("R output: %f \n",test4[i]);
//    
//    double* test5 = EulerAngles(vv);
//    for(int i = 0; i<2; i++)
//        printf("euler output: %f \n",test5[i]);
}

double Kernel::GetMaximumValue()
{
    double x[3] = {0,0,0};
    double y[3] = {0,0,0};
    double r[3] = {0,0,1};
    double v[3] = {0,0,1};
    
    // evaluate at origin
    double kernelMax = k2(x, y, r, v);
    printf("max kernel val: %f \n",kernelMax);

    return kernelMax;
}

void Kernel::EstimateKernelSize()
{
    double x[3] = {0,0,0};
    double y[3] = {0,0,0};
    double r[3] = {0,0,1};
    double v[3] = {0,0,1};
    
    // evaluate at origin
    double kernelMax = k2(x, y, r, v);
    printf("max kernel val: %f \n",kernelMax);
    
    // determine a good kernel size NxMxM
    // move along z-axis
    double i = 0.0;
    while(true)
    {
        i += 0.1;
        x[2] = i;
        double kval = k2(x,y,r,v)/kernelMax;
        //printf("i:%f val: %f \n",i,kval);
        if(kval < 0.1)
            break;
    }
    int N = ceil(i)*2;
    if(N%2 == 0)
        N -= 1;
    //printf("--------\n");
    
    // move along x-axis
    i = 0.0;
    x[2] = 0;
    while(true)
    {
        i += 0.1;
        x[0] = i;
        double kval = k2(x,y,r,v)/kernelMax;
        //printf("i:%f val: %f \n",i,kval);
        if(kval < 0.1)
            break;
    }
    int M = ceil(i)*2;
    if(M%2 == 0)
        M -= 1;
    
    N = std::max(N,M);
    M = std::max(N,M);
    
    printf("the kernel is %dx%dx%d\n",M,M,N);
    this->xsize = M;
    this->ysize = M;
    this->zsize = N;
}

double Kernel::k2(double* x, double* y, double* r, double* v)
{
    double* a = Subtract3(x,y);
    double* transm = Transpose3x3(R(EulerAngles(v)));
    double* arg1 = Multiply(transm,a);
    double* arg2p = Multiply(transm,r);
    double* arg2 = EulerAngles(arg2p);
    
    if(isnan(arg2[0]))
    {
        printf("nan! (%f %f %f) ---> (%f %f) \n",arg2p[0],arg2p[1],arg2p[2],arg2[0],arg2[1]);

    }
    
    double* c = coordinateMap(arg1[0], arg1[1], arg1[2], arg2[0], arg2[1]);
    double kernelval = kernel(c);

    
    //kernelval *= sin(eulerv[0]); // "usual" surface measure for sphere
    
    free(arg1);
    free(a);
    free(arg2);
    free(arg2p);
    free(c);
    free(transm);
    
    return kernelval;

}

inline double* Kernel::coordinateMap(double x, double y, double z, double beta, double gamma)
{
    // new coordinate chart of Jorg
    
    //printf("coords: %f %f %f %f %f \n",x,y,z,beta,gamma);
    double* c = (double*)malloc(6*sizeof(double));
    
    if(beta == 0)
    {
        c[0] = x;
        c[1] = y;
        c[2] = z;
        c[3] = c[4] = c[5] = 0;
    }
    else
    {
        double q = fabs(beta);
        double cg = cos(gamma);
        double sg = sin(gamma);
        double cotq2 = cot(q/2);
        
        c[0] = -0.5*z*beta*cg + x*(1-(beta*beta*cg*cg*(1 - 0.5*q*cotq2))/(q*q)) - (y*beta*beta*cg*(1-0.5*q*cotq2)*sg)/(q*q);
        
        c[1] = -0.5*z*beta*sg - (x*beta*beta*cg*(1-0.5*q*cotq2)*sg)/(q*q) + y*(1-(beta*beta*(1-0.5*q*cotq2)*sg*sg)/(q*q));
        
        c[2] = 0.5*x*beta*cg + 0.5*y*beta*sg + z*(1+((1-0.5*q*cotq2)*(-beta*beta*cg*cg - beta*beta*sg*sg))/(q*q));
        
        c[3] = beta * (-sg);
        
        c[4] = beta * cg;
        
        c[5] = 0;
    }
    
    return c;
}

inline double Kernel::kernel(double* c)
{
    return 1/(8*sqrt(2))*sqrt(PI)*t*sqrt(t*D33)*sqrt(D33*D44) *
            1/(16*PI*PI*D33*D33*D44*D44*t*t*t*t) *
            exp(-sqrt( (c[0]*c[0] + c[1]*c[1])/(D33*D44) + (c[2]*c[2]/D33 + (c[3]*c[3]+c[4]*c[4])/D44)*(c[2]*c[2]/D33 + (c[3]*c[3]+c[4]*c[4])/D44) + c[5]*c[5]/D44)/(4*t));
}

void Kernel::DefineOrientations()
{
    double orientationsList[486] = {0.0,0.0,1.,
        -0.20318274393026847,0.14762090442216314,0.9679487802288661,
        -0.42532540417602,0.3090169943749474,0.8506508083520399,
        -0.6095482317908053,0.4428627132664893,0.6575131712132815,
        -0.7236067977499789,0.5257311121191336,0.4472135954999579,
        0.07760890225389615,0.2388556408050596,0.9679487802288659,
        -0.1381966011250105,0.42532540417601994,0.8944271909999157,
        -0.36180339887498947,0.587785252292473,0.7236067977499789,
        -0.531939329536909,0.6817183540715489,0.5022953667054891,
        0.16245984811645317,0.5,0.8506508083520399,
        -0.052786404500042065,0.6881909602355867,0.7236067977499789,
        -0.2628655560595668,0.8090169943749473,0.5257311121191336,
        0.23282670676168846,0.7165669224151787,0.6575131712132815,
        0.02964396283142001,0.8641878268373419,0.5022953667054891,
        0.27639320225002106,0.8506508083520399,0.4472135954999579,
        -0.20318274393026847,-0.14762090442216314,0.9679487802288661,
        -0.42532540417602,-0.3090169943749474,0.8506508083520399,
        -0.6095482317908053,-0.4428627132664893,0.6575131712132815,
        -0.7236067977499789,-0.5257311121191336,0.4472135954999579,
        -0.4472135954999579,0.,0.8944271909999157,
        -0.6708203932499369,-0.16245984811645314,0.7236067977499789,
        -0.8127309757210738,-0.2952418088443262,0.5022953667054891,
        -0.6708203932499369,0.16245984811645314,0.7236067977499789,
        -0.85065080835204,0.,0.5257311121191336,
        -0.8127309757210738,0.2952418088443262,0.5022953667054891,
        0.07760890225389615,-0.2388556408050596,0.9679487802288659,
        0.16245984811645317,-0.5,0.8506508083520399,
        0.23282670676168846,-0.7165669224151787,0.6575131712132815,
        0.27639320225002106,-0.8506508083520399,0.4472135954999579,
        -0.1381966011250105,-0.42532540417601994,0.8944271909999157,
        -0.052786404500042065,-0.6881909602355867,0.7236067977499789,
        0.02964396283142001,-0.8641878268373419,0.5022953667054891,
        -0.36180339887498947,-0.587785252292473,0.7236067977499789,
        -0.2628655560595668,-0.8090169943749473,0.5257311121191336,
        -0.531939329536909,-0.6817183540715489,0.5022953667054891,
        0.2511476833527446,0.,0.9679487802288661,
        0.5257311121191336,0.,0.8506508083520399,
        0.7534430500582336,0.,0.6575131712132815,
        0.8944271909999159,0.,0.4472135954999579,
        0.36180339887498947,-0.2628655560595668,0.8944271909999157,
        0.6381966011250104,-0.2628655560595668,0.7236067977499789,
        0.8310519523121299,-0.23885564080505955,0.5022953667054891,
        0.44721359549995787,-0.5257311121191336,0.7236067977499789,
        0.6881909602355868,-0.5,0.5257311121191336,
        0.48397439011443294,-0.7165669224151788,0.5022953667054891,
        0.36180339887498947,0.2628655560595668,0.8944271909999157,
        0.44721359549995787,0.5257311121191336,0.7236067977499789,
        0.48397439011443294,0.7165669224151788,0.5022953667054891,
        0.6381966011250104,0.2628655560595668,0.7236067977499789,
        0.6881909602355868,0.5,0.5257311121191336,
        0.8310519523121299,0.23885564080505955,0.5022953667054891,
        0.7236067977499789,-0.5257311121191336,-0.4472135954999579,
        0.6095482317908053,-0.4428627132664893,-0.6575131712132815,
        0.42532540417602,-0.3090169943749474,-0.8506508083520399,
        0.20318274393026847,-0.14762090442216314,-0.9679487802288661,
        0.,0.,-1.,
        0.531939329536909,-0.6817183540715489,-0.5022953667054891,
        0.36180339887498947,-0.587785252292473,-0.7236067977499789,
        0.1381966011250105,-0.42532540417601994,-0.8944271909999157,
        -0.07760890225389615,-0.2388556408050596,-0.9679487802288659,
        0.2628655560595668,-0.8090169943749473,-0.5257311121191336,
        0.052786404500042065,-0.6881909602355867,-0.7236067977499789,
        -0.16245984811645317,-0.5,-0.8506508083520399,
        -0.02964396283142001,-0.8641878268373419,-0.5022953667054891,
        -0.23282670676168846,-0.7165669224151787,-0.6575131712132815,
        -0.27639320225002106,-0.8506508083520399,-0.4472135954999579,
        0.7236067977499789,0.5257311121191336,-0.4472135954999579,
        0.6095482317908053,0.4428627132664893,-0.6575131712132815,
        0.42532540417602,0.3090169943749474,-0.8506508083520399,
        0.20318274393026847,0.14762090442216314,-0.9679487802288661,
        0.8127309757210738,0.2952418088443262,-0.5022953667054891,
        0.6708203932499369,0.16245984811645314,-0.7236067977499789,
        0.4472135954999579,0.,-0.8944271909999157,
        0.85065080835204,0.,-0.5257311121191336,
        0.6708203932499369,-0.16245984811645314,-0.7236067977499789,
        0.8127309757210738,-0.2952418088443262,-0.5022953667054891,
        -0.27639320225002106,0.8506508083520399,-0.4472135954999579,
        -0.23282670676168846,0.7165669224151787,-0.6575131712132815,
        -0.16245984811645317,0.5,-0.8506508083520399,
        -0.07760890225389615,0.2388556408050596,-0.9679487802288659,
        -0.02964396283142001,0.8641878268373419,-0.5022953667054891,
        0.052786404500042065,0.6881909602355867,-0.7236067977499789,
        0.1381966011250105,0.42532540417601994,-0.8944271909999157,
        0.2628655560595668,0.8090169943749473,-0.5257311121191336,
        0.36180339887498947,0.587785252292473,-0.7236067977499789,
        0.531939329536909,0.6817183540715489,-0.5022953667054891,
        -0.8944271909999159,0.,-0.4472135954999579,
        -0.7534430500582336,0.,-0.6575131712132815,
        -0.5257311121191336,0.,-0.8506508083520399,
        -0.2511476833527446,0.,-0.9679487802288661,
        -0.8310519523121299,0.23885564080505955,-0.5022953667054891,
        -0.6381966011250104,0.2628655560595668,-0.7236067977499789,
        -0.36180339887498947,0.2628655560595668,-0.8944271909999157,
        -0.6881909602355868,0.5,-0.5257311121191336,
        -0.44721359549995787,0.5257311121191336,-0.7236067977499789,
        -0.48397439011443294,0.7165669224151788,-0.5022953667054891,
        -0.48397439011443294,-0.7165669224151788,-0.5022953667054891,
        -0.44721359549995787,-0.5257311121191336,-0.7236067977499789,
        -0.36180339887498947,-0.2628655560595668,-0.8944271909999157,
        -0.6881909602355868,-0.5,-0.5257311121191336,
        -0.6381966011250104,-0.2628655560595668,-0.7236067977499789,
        -0.8310519523121299,-0.23885564080505955,-0.5022953667054891,
        0.1552178045077923,0.9554225632202383,0.25114768335274457,
        -0.1381966011250105,0.9510565162951535,0.276393202250021,
        -0.4472135954999579,0.8506508083520399,0.276393202250021,
        -0.6871571340447014,0.6817183540715489,0.25114768335274457,
        0.,1.,0.,
        -0.3090169943749474,0.9510565162951535,0.,
        -0.5877852522924731,0.8090169943749473,0.,
        -0.1552178045077923,0.9554225632202383,-0.25114768335274457,
        -0.4360094506919568,0.8641878268373419,-0.25114768335274457,
        -0.8606959151435498,0.4428627132664894,0.2511476833527446,
        -0.9472135954999581,0.1624598481164532,0.2763932022500211,
        -0.9472135954999581,-0.1624598481164532,0.2763932022500211,
        -0.8606959151435498,-0.4428627132664894,0.2511476833527446,
        -0.9510565162951535,0.3090169943749474,0.,
        -0.9999999999999999,0.,0.,
        -0.9510565162951535,-0.3090169943749474,0.,
        -0.9566257939885021,0.1476209044221631,-0.25114768335274457,
        -0.9566257939885021,-0.1476209044221631,-0.25114768335274457,
        -0.6871571340447014,-0.6817183540715489,0.25114768335274457,
        -0.4472135954999579,-0.8506508083520399,0.276393202250021,
        -0.1381966011250105,-0.9510565162951535,0.276393202250021,
        0.1552178045077923,-0.9554225632202383,0.25114768335274457,
        -0.5877852522924731,-0.8090169943749473,0.,
        -0.3090169943749474,-0.9510565162951535,0.,
        0.,-1.,0.,
        -0.4360094506919568,-0.8641878268373419,-0.25114768335274457,
        -0.1552178045077923,-0.9554225632202383,-0.25114768335274457,
        0.4360094506919568,-0.8641878268373419,0.25114768335274457,
        0.6708203932499369,-0.6881909602355867,0.276393202250021,
        0.8618033988749894,-0.42532540417601994,0.276393202250021,
        0.9566257939885021,-0.1476209044221631,0.25114768335274457,
        0.5877852522924731,-0.8090169943749473,0.,
        0.8090169943749475,-0.587785252292473,0.,
        0.9510565162951535,-0.3090169943749474,0.,
        0.6871571340447014,-0.6817183540715489,-0.25114768335274457,
        0.8606959151435498,-0.4428627132664894,-0.2511476833527446,
        0.9566257939885021,0.1476209044221631,0.25114768335274457,
        0.8618033988749894,0.42532540417601994,0.276393202250021,
        0.6708203932499369,0.6881909602355867,0.276393202250021,
        0.4360094506919568,0.8641878268373419,0.25114768335274457,
        0.9510565162951535,0.3090169943749474,0.,
        0.8090169943749475,0.587785252292473,0.,
        0.5877852522924731,0.8090169943749473,0.,
        0.8606959151435498,0.4428627132664894,-0.2511476833527446,
        0.6871571340447014,0.6817183540715489,-0.25114768335274457,
        0.4472135954999579,-0.8506508083520399,-0.276393202250021,
        0.1381966011250105,-0.9510565162951535,-0.276393202250021,
        0.3090169943749474,-0.9510565162951535,0.,
        0.9472135954999581,0.1624598481164532,-0.2763932022500211,
        0.9472135954999581,-0.1624598481164532,-0.2763932022500211,
        1.,0.,0.,
        0.1381966011250105,0.9510565162951535,-0.276393202250021,
        0.4472135954999579,0.8506508083520399,-0.276393202250021,
        0.3090169943749474,0.9510565162951535,0.,
        -0.8618033988749894,0.42532540417601994,-0.276393202250021,
        -0.6708203932499369,0.6881909602355867,-0.276393202250021,
        -0.8090169943749475,0.587785252292473,0.,
        -0.6708203932499369,-0.6881909602355867,-0.276393202250021,
        -0.8618033988749894,-0.42532540417601994,-0.276393202250021,
        -0.8090169943749475,-0.587785252292473,0.
    };
    orientations = (double*)malloc(486*sizeof(double));
    std::memcpy(orientations,&orientationsList[0],486*sizeof(double)); // copy to heap
    
    //orientations = orientationsList;
}


