#include "materialPar.h"
#include <mpi.h>
#include <stdio.h>

// initialize the local material (NxLocal x Ny) to have basic data, set alpha value, figure out 
// padding and starting index rows
int initMaterialLoc(materialLoc *aMaterial, unsigned int Nx, unsigned int NyTotal, unsigned int nPadRows, float dx, float dy, float alpha){
    // check rank and number of processes
    int rank, nProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // information about how many columns and spacing between each column in grid
    aMaterial->Nx = Nx;
    aMaterial->dx = dx;
    
    // information about how many rows overall, in this unpadded subset (NxLocal), how much padding
    aMaterial->NyTotal = NyTotal; // total rows in the entire distributed material grid
    aMaterial->dy = dy; // physical spacing between rows (float)
    
    // ===================STUDENT CODE START HERE==============================================

    // standard number of rows in any subset (integer division) for NyLocal

    // number of rows in global material in lower ranks for startYId

    // if not evenly divisible by number of processes, add one extra row to some processes
      // if a late enough rank to not have an extra row add rows for earlier processes with one extra row each

      // if rank is less than the remainder of NyTotal divided by number of processes total

        // add an extra row for this process

        // every rank before this one had one extra row, so start index is pushed back
       
    // ===================STUDENT CODE END HERE===============================================
    aMaterial->nPadRows = nPadRows;
    aMaterial->NyPadded = nPadRows + aMaterial->NyLocal + nPadRows; // nPadRows rows, then NyLocal rows, then nPadRows in the local padded subarray
    
    // check that the number of padding rows isn't bigger than the number of local rows
    if(nPadRows > aMaterial->NyLocal){
        printf("WARNING: nPadRows must be <= NyLocal");
        return 1;
    }
    
    // alpha parameter that governs how quickly heat spreads out
    aMaterial->alpha = alpha;
    
    return 0;    
};
