#include "mex.h"
#include <stdio.h>
#include <math.h>

void Vert_Neib(double * Tri, double * faces, double Nv, double Nf)
{
    double posc;
    double ind;
    double posc2;
    
    for (int i = 0; i < Nv; i++)
    {
        Tri[i] = (double)i+1;
        Tri[(int)(i+Nv)] = 0;
    }
               
    for (int j = 0; j < Nf; j++)
        {
            ind = faces[(int)(j)];
            posc = ind+Nv-1;
            Tri[(int)(posc)] = Tri[(int)(posc)]+1;
            posc2 = posc+Tri[(int)(posc)]*Nv;
            Tri[(int)(posc2)] = j+1;
            
            ind = faces[(int)(j+Nf)];
            posc = ind+Nv-1;
            Tri[(int)(posc)] = Tri[(int)(posc)]+1;
            posc2 = posc+Tri[(int)(posc)]*Nv;
            Tri[(int)(posc2)] = j+1;
            
            ind = faces[(int)(j+2*Nf)];
            posc = ind+Nv-1;
            Tri[(int)(posc)] = Tri[(int)(posc)]+1;
            posc2 = posc+Tri[(int)(posc)]*Nv;
            Tri[(int)(posc2)] = j+1;
        }
        
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    double * Nv;
    double * Nf;
    double * faces;
    double * Tri;
    
    // Arguments:
    faces     = mxGetPr(prhs[0]);
    Nv        = mxGetPr(prhs[1]);
    Nf        = mxGetPr(prhs[2]);
    
    // Return Matrices:
    plhs[0] = mxCreateDoubleMatrix((int)Nv[0], 100, mxREAL); //crea la matriz inicializada en cero
    Tri = mxGetPr(plhs[0]);
    
    Vert_Neib(Tri, faces, Nv[0], Nf[0]);
}