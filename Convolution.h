//
//  Convolution.h
//  Convolutions
//
//  Created by Stephan Meesters on 18/08/15.
//
//

#ifndef __Convolutions__Convolution__
#define __Convolutions__Convolution__

#include "Kernel.h"

/** Includes -- NIfTI */
#include "nifti1.h"
#include "nifti1_io.h"

#include <cstring>
#include <limits>

template <class DType> 
class Convolution
{
private:
    
    Kernel* kernel;
    std::string path_dsf;
    std::string path_output;
    int numIterations;
    bool padding;
    bool saveIterations;
    nifti_image* DSF;
    double* dataD;
    double* dataDEnhanced;
    int nvox;
    
    void NormalizeToDatatypeRange();
    void ExportData(std::string outputName);
    
public:
    
    Convolution(Kernel* kernel,
                std::string path_dsf,
                std::string path_output,
                int iterations,
                bool padding,
                bool saveIterations);
    ~Convolution();
};


/* FUNCTION DECLARATIONS */

template <class DType>
Convolution<DType>::Convolution(Kernel* kernel,
                                std::string path_dsf,
                                std::string path_output,
                                int numIterations,
                                bool padding,
                                bool saveIterations)
{
    this->kernel = kernel;
    this->path_dsf = path_dsf;
    this->path_output = path_output;
    this->numIterations = numIterations;
    this->padding = padding;
    this->saveIterations = saveIterations;
    
    // read input data
    this->DSF = nifti_image_read(path_dsf.c_str(), 1);
    int nx = DSF->nx;
    int ny = DSF->ny;
    int nz = DSF->nz;
    int norient = DSF->nu;
    nvox = DSF->nvox;
    DType* data = (DType*)DSF->data;
    
    // convert data to double for precision
    dataD = (double*)malloc(nvox*sizeof(double));
    std::copy(data,data+nvox,dataD);
    
    // define kernel dimensions
    int kernelx = kernel->xsize;
    int kernely = kernel->ysize;
    int kernelz = kernel->zsize;
    int kernelxh = kernel->xsize/2;
    int kernelyh = kernel->ysize/2;
    int kernelzh = kernel->zsize/2;
    
    // create an empty dataset to store the output
    dataDEnhanced = (double*)calloc(nvox,sizeof(double)); // zero initialize
    
    // convolution including zero-padding for the boundaries
    if(padding)
    {
        printf("applying zero padding.\n");
        
        // loop over FOD cx,cy,cz
        for(int iterations = 0; iterations<numIterations; iterations++)
        {
#pragma omp parallel for
            for (int corient = 0; corient<norient; corient++)
            {
                for (int cy = 0; cy<ny-1; cy++)
                {
                    for (int cz = 0; cz<nz-1; cz++)
                    {
                        for(int cx = 0; cx<nx-1; cx++)
                        {
                            // loop over kernel x,y,z
                            double output = 0.0;
                            for (int x = cx-kernelxh; x<=cx+kernelxh; x++)
                            {
                                for (int y = cy-kernelyh; y<=cy+kernelyh; y++)
                                {
                                    for (int z = cz-kernelzh; z<=cz+kernelzh; z++)
                                    {
                                        for (int orient = 0; orient < norient; orient++)
                                        {
                                            // skip iteration if the data is outside bounds
                                            if(cy < 0 || cy >= ny ||
                                               cx < 0 || cx >= nx ||
                                               cz < 0 || cz >= nz )
                                                continue;
                                            
                                            double U = dataD[x + y*nx + z*nx*ny + orient*nx*ny*nz];
                                            // ez centered around corient
                                            double pt = kernel->lut[
                                                                    (x-(cx-kernelxh)) +
                                                                    (y-(cy-kernelyh))*kernelx +
                                                                    (z-(cz-kernelzh))*kernelx*kernely +
                                                                    orient*kernelx*kernely*kernelz +
                                                                    corient*norient*kernelx*kernely*kernelz
                                                                    ];
                                            output += pt*U;
                                        }
                                    }
                                }
                            }
                            dataDEnhanced[corient * nx*ny*nz + cz*nx*ny + cy*nx + cx] = output;
                        }
                    }
                }
            }
            if(numIterations > 1)
                std::memcpy(dataD, dataDEnhanced, nvox*sizeof(double));
        }
    }
    
    // convolution without padding (does not compute the boundaries)
    else
    {
        // loop over FOD cx,cy,cz
        for(int iterations = 0; iterations<numIterations; iterations++)
        {
#pragma omp parallel for
            for (int corient = 0; corient<norient; corient++)
            {
                for (int cy = kernelyh+1; cy<ny-kernelyh-1; cy++)
                {
                    for (int cz = kernelzh+1; cz<nz-kernelzh-1; cz++)
                    {
                        for(int cx = kernelxh+1; cx<nx-kernelxh-1; cx++)
                        {
                            // loop over kernel x,y,z
                            double output = 0.0;
                            for (int x = cx-kernelxh; x<=cx+kernelxh; x++)
                            {
                                for (int y = cy-kernelyh; y<=cy+kernelyh; y++)
                                {
                                    for (int z = cz-kernelzh; z<=cz+kernelzh; z++)
                                    {
                                        for (int orient = 0; orient < norient; orient++)
                                        {
                                            double U = dataD[x + y*nx + z*nx*ny + orient*nx*ny*nz];
                                            // ez centered around corient
                                            double pt = kernel->lut[
                                                                    (x-(cx-kernelxh)) +
                                                                    (y-(cy-kernelyh))*kernelx +
                                                                    (z-(cz-kernelzh))*kernelx*kernely +
                                                                    orient*kernelx*kernely*kernelz +
                                                                    corient*norient*kernelx*kernely*kernelz
                                                                    ];
                                            output += pt*U;
                                        }
                                    }
                                }
                            }
                            dataDEnhanced[corient * nx*ny*nz + cz*nx*ny + cy*nx + cx] = output;
                        }
                    }
                }
            }
            if(numIterations > 1)
            {
                // normalize within range of datatype. numerical inprecision?
                NormalizeToDatatypeRange();
                if(saveIterations)
                {
                    std::string trunc = path_output.substr(0,path_output.size()-4);
                    char outName[100];
                    sprintf(outName, "%s_iteration%02d.nii",trunc.c_str(),iterations+1);
                    ExportData(outName);
                }
                std::memcpy(dataD, dataDEnhanced, DSF->nvox*sizeof(double));
            }
        }
    }
    
    // output the endresult
    ExportData(path_output);
}
template <class DType>
Convolution<DType>::~Convolution()
{
    this->kernel = NULL;
    free(dataDEnhanced);
    free(DSF);
    free(dataD);
}

template <class DType>
void Convolution<DType>::NormalizeToDatatypeRange()
{
    // check datatype bounds
    double maxval = 0;
    for(int i = 0; i<nvox; i++)
    {
        if(dataDEnhanced[i] > maxval)
            maxval = dataDEnhanced[i];
    }
    
    // bring back into datatype range, if necessary
    if(maxval > std::numeric_limits<DType>::max())
    {
        DType maxd = std::numeric_limits<DType>::max();
        for(int i = 0; i<nvox; i++)
        {
            dataDEnhanced[i] = dataDEnhanced[i] / maxval * maxd;
        }
    }
}

template <class DType>
void Convolution<DType>::ExportData(std::string outputName)
{
    // convert back to target datatype
    DType* dataSEnhanced = (DType*)malloc(DSF->nvox*sizeof(DType));
    std::copy(dataDEnhanced,dataDEnhanced+DSF->nvox,dataSEnhanced);
    
    // export end result NIfTI
    nifti_image* DSF_output = nifti_image_read(path_dsf.c_str(), 1); // load the header of the original DSF
    DSF_output->data = dataSEnhanced;
    DSF_output->fname = (char*)outputName.c_str();
    nifti_image_write(DSF_output);
    
    free(dataSEnhanced);
    free(DSF_output);
}







#endif /* defined(__Convolutions__Convolution__) */
