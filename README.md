# Convolutions

A C++ implementation of shift-twist convolutions for diffusion data.

USAGE: 

   ./Convolutions  [--save-iterations <Save-Iterations>] [--padding <Apply
                   padding>] [--iterations <Iterations>] --output <Path to
                   output .vtk file> --data <Path to .nii file> [--d33
                   <Value of D33 parameter>] [--d44 <Value of D44
                   parameter>] [--t <Value of diffusion time parameter>]
                   [--verbose] [--] [--version] [-h]


Where: 

   --save-iterations <Save-Iterations>
     Save the result after each iteration

   --padding <Apply padding>
     Apply padding. Default value: false

   --iterations <Iterations>
     Number of convolution iterations

   --output <Path to output .vtk file>
     (required)  Path to the output .vtk file

   --data <Path to .nii file>
     (required)  Path to the .nii file containg the DSF

   --d33 <Value of D33 parameter>
     Value of D33 parameter. Default value: 1.0

   --d44 <Value of D44 parameter>
     Value of D44 parameter. Default value: 0.04

   --t <Value of diffusion time parameter>
     Value of diffusion time parameter. Default value: 1.4

   --verbose
     Verbose mode

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.
