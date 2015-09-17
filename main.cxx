//
//  main.cpp
//  Convolutions
//
//  Created by Stephan Meesters on 14/08/15.
//  Copyright (c) 2015 Stephan Meesters. All rights reserved.
//

#include <iostream>

/** Includes - TCLAP */
#include "tclap/CmdLine.h"

/** Includes -- project */
#include "Kernel.h"
#include "Convolution.h"

/** Includes -- NIfTI */
#include "nifti1.h"
#include "nifti1_io.h"

/** Includes OpenMP */
#ifdef __linux__
#include "omp.h"
#endif

int main(int argc, const char * argv[])
{
    try
    {
        // Command line arguments setup
        TCLAP::CmdLine cmd("Command description message", ' ', "0.1");
        TCLAP::ValueArg<std::string> path_dsf("",
                                              "data",
                                              "Path to the .nii file containg the DSF",
                                              true,
                                              "",
                                              "Path to .nii file");
        TCLAP::ValueArg<std::string> path_output("",
                                              "output",
                                              "Path to the output .vtk file",
                                              true,
                                              "",
                                              "Path to output .vtk file");
        
        TCLAP::ValueArg<float> param_d33("",
                                              "d33",
                                              "Value of D33 parameter. Default value: 1.0",
                                              false,
                                              1.0,
                                              "Value of D33 parameter");
        TCLAP::ValueArg<float> param_d44("",
                                         "d44",
                                         "Value of D44 parameter. Default value: 0.04",
                                         false,
                                         0.04,
                                         "Value of D44 parameter");
        TCLAP::ValueArg<float> param_t("",
                                         "t",
                                         "Value of diffusion time parameter. Default value: 1.4",
                                         false,
                                         2,
                                         "Value of diffusion time parameter");
        
        TCLAP::ValueArg<bool> param_padding("",
                                       "padding",
                                       "Apply padding. Default value: false",
                                       false,
                                       false,
                                       "Apply padding");
        
        TCLAP::ValueArg<int> param_iterations("",
                                            "iterations",
                                            "Number of convolution iterations",
                                            false,
                                            1,
                                            "Iterations");
        
        TCLAP::ValueArg<int> param_saveiterations("",
                                              "save-iterations",
                                              "Save the result after each iteration",
                                              false,
                                              1,
                                              "Save-Iterations");
        
        TCLAP::ValueArg<int> param_numthreads("",
                                                  "numthreads",
                                                  "Number of threads for parallel processing. Default=max",
                                                  false,
                                                  0,
                                                  "Num Threads");
        

        TCLAP::SwitchArg param_verbose("","verbose","Verbose mode",false);
        
        cmd.add( param_verbose );
        cmd.add( param_t );
        cmd.add( param_d44 );
        cmd.add( param_d33 );
        cmd.add( path_dsf );
        cmd.add( path_output);
        cmd.add( param_iterations );
        cmd.add( param_padding );
        cmd.add( param_saveiterations );
        cmd.add( param_numthreads );

        // Parse the args.
        cmd.parse( argc, argv );
        
        // Configure OpenMP multithreading (only available with GCC in Linux)
#ifdef __linux__
        if(param_numthreads.getValue() != 0)
        {
            omp_set_dynamic(0);     // Explicitly disable dynamic teams
            omp_set_num_threads(param_numthreads.getValue()); // Nr of threads for all consecutive parallel regions
        }
#endif
        
        // Initialize kernel, calculate it or load it
        Kernel* kernel = new Kernel(param_d33.getValue(),
                                    param_d44.getValue(),
                                    param_t.getValue());
        
        // Check the datatype -- read header
        nifti_image* DSF = nifti_image_read(path_dsf.getValue().c_str(), 0);
        
        // Perform convolution
        switch (DSF->datatype)
        {
            case DT_UNSIGNED_CHAR:   // "Byte"
            {
                Convolution<char>* convolution = new Convolution<char>(kernel,
                                                                         path_dsf.getValue(),
                                                                         path_output.getValue(),
                                                                         param_iterations.getValue(),
                                                                         param_padding.getValue(),
                                                                         param_saveiterations.getValue());
                free(convolution);
                break;
            }
                
            case DT_SIGNED_SHORT:   // "Integer16"
            {
                Convolution<short>* convolution = new Convolution<short>(kernel,
                                                           path_dsf.getValue(),
                                                           path_output.getValue(),
                                                           param_iterations.getValue(),
                                                           param_padding.getValue(),
                                                           param_saveiterations.getValue());
                free(convolution);
                break;
            }
                
            case DT_SIGNED_INT:     // "Integer32"
            {
                Convolution<int>* convolution = new Convolution<int>(kernel,
                                                                         path_dsf.getValue(),
                                                                         path_output.getValue(),
                                                                         param_iterations.getValue(),
                                                                         param_padding.getValue(),
                                                                         param_saveiterations.getValue());
                free(convolution);
                break;
            }
                
            case DT_FLOAT:          // "Real32"
            {
                Convolution<float>* convolution = new Convolution<float>(kernel,
                                                                         path_dsf.getValue(),
                                                                         path_output.getValue(),
                                                                         param_iterations.getValue(),
                                                                         param_padding.getValue(),
                                                                         param_saveiterations.getValue());
                free(convolution);
                break;
            }
                
            case DT_DOUBLE:         // "Real64"
            {
                Convolution<double>* convolution = new Convolution<double>(kernel,
                                                                         path_dsf.getValue(),
                                                                         path_output.getValue(),
                                                                         param_iterations.getValue(),
                                                                         param_padding.getValue(),
                                                                         param_saveiterations.getValue());
                free(convolution);
                break;
            }
                
            case DT_UINT16:         // "UnsignedInteger16"
            {
                Convolution<unsigned int>* convolution = new Convolution<unsigned int>(kernel,
                                                                         path_dsf.getValue(),
                                                                         path_output.getValue(),
                                                                         param_iterations.getValue(),
                                                                         param_padding.getValue(),
                                                                         param_saveiterations.getValue());
                free(convolution);
                break;
            }
                
//            case DT_UINT32:         // "UnsignedInteger32" -- not supported by GCC 4.4.7 used on the BigMath2
//            {
//                Convolution<uint32_t>* convolution = new Convolution<uint32_t>(kernel,
//                                                                         path_dsf.getValue(),
//                                                                         path_output.getValue(),
//                                                                         param_iterations.g.etValue(),
//                                                                         param_padding.getValue(),
//                                                                         param_saveiterations.getValue());
//                free(convolution);
//                break;
//            }
                
            default:
            {
                printf("This datatype (code:%d) is not supported. Aborting.\n",DSF->datatype);
                break;
            }
        }
        
        // Clean up
        free(kernel);       
        
        printf("end\n");
        
        
    }
    catch (TCLAP::ArgException &e)  // catch any exceptions
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }
    
    
    return 0;
}












