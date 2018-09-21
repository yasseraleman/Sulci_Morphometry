/*********************************************************************
 * Demo.cpp
 *
 * This file shows the basics of setting up a mex file to work with
 * Matlab.  This example shows how to use 2D matricies.  This may
 * 
 * Keep in mind:
 * <> Use 0-based indexing as always in C or C++
 * <> Indexing is column-based as in Matlab (not row-based as in C)
 * <> Use linear indexing.  [x*dimy+y] instead of [x][y]
 *
 * For more information, see my site: www.shawnlankton.com
 * by: Shawn Lankton
 *
 ********************************************************************/
#include <matrix.h>
#include <mex.h>   

/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

//declare variables
    mxArray *a_in_m, *c_out_m;
    const mwSize *dims;
    double *Trip, *Neigh;
    int dimx, dimy, numdims;
    int i,z, indexs, index2, t, dimtemp;
    double temp;
    
//associate inputs
    a_in_m = mxDuplicateArray(prhs[0]);
    
//figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    numdims = mxGetNumberOfDimensions(prhs[0]);
    dimx = (int)dims[0]; dimy = (int)dims[1];
    
//associate outputs
    c_out_m = plhs[0] = mxCreateDoubleMatrix(dimx,dimy-1,mxREAL);
    
    //c_out_m = plhs[0] = mxCreateNumericMatrix(dimx, dimy, mxINT32_CLASS, mxREAL);
    
//associate pointers
    Trip = mxGetPr(a_in_m);
    Neigh = mxGetPr(c_out_m);
    
//do something
    //mexPrintf("element = %u  == >> %u\n",dimx,dimy);
    for(i=0;i<dimx;i++)                             // Rows Loop
    {
        dimtemp = int(Trip[int(dimx)+int(i)]);      // Number of elements different from zero inside the Trip row
        for(z=0;z<dimtemp;z++)                      // For each column with an elemnt different from zero
        {
            index2 = dimx*(2+z)+i; // Index in the Trip matrix
            temp = Trip[index2];   // Trip Value
            t = int(temp)-1;       // Float to integer
            indexs = int(Neigh[t])*dimx+t; // Index in the Neighbors Matrix
             Neigh[int(indexs)+int(dimx)]=index2+1;  // Store the Trip Index in the Neighbors Matrix
            Neigh[t] = Neigh[t] + 1;                 // Increase the number of Neighbors
        } 
    }
    
    return;
}
