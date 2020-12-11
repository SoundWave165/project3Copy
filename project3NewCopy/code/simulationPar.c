#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "simulationPar.h"
#include "materialPar.h"
#include "checkPtPar.h"
#include <mpi.h>

// Calculate the maximum stable time step allowed following CFL condition
float calcMaxTimeStepLoc(simLoc *thisSimLoc){
	// CFL condition for 2D heat equation is dt <= min(dx^2,dy^2)/(4*alpha)
	// get properties of the material and discretization
	float alpha = (thisSimLoc->thisMaterialLoc)->alpha;
	float dx = (thisSimLoc->thisMaterialLoc)->dx;
	float dy = (thisSimLoc->thisMaterialLoc)->dy;

	// figure out the min step size squared
	float minStepSq = dx*dx;
	if(dy < dx) minStepSq = dy*dy;
	float dtMax = minStepSq/(4*alpha);

	return dtMax;
};

// Initialize the loccal initial state (temperature matrix at T = 0, copied from valsForInitStateGlobal) 
// and current state of the system (matrix of temperature values).
// valsForInitStateGlobal will be thisMaterial.Nx x thisMaterial.NyTotal, and only the values relevant 
// to this process's local subset of the material should be copied in. 
// Note: initializing the simulation does not also initialize the material. Do that separately.
int initSimLoc(simLoc *thisSimLoc, float timeStep, float *valsForInitStateGlobal, float bdryVal, materialLoc *thisMaterialLoc){
	// return flag, 0 if okay, 1 if not
	int flag = 0;

	// tie this simulation to a material
	thisSimLoc->thisMaterialLoc = thisMaterialLoc;	

	// calculate maximum allowed stable time step (CFL condition)
	thisSimLoc->dtMax = calcMaxTimeStepLoc(thisSimLoc); // maximum stable time step
	thisSimLoc->dt = timeStep;
	if(timeStep >= thisSimLoc->dtMax){  // if too large a time step is requested, send out a warning
		printf("WARNING: In initSimLoc(), requested time step exceeds stability limit. Unphysical behavior is likely. \n");
		flag = 1;
	} 

	// start out at time 0
	thisSimLoc->currentTimeIdx = 0;
	// calculate starting index for where to start copying this subarray from the global initial state array
	int startIdx = (thisSimLoc->thisMaterialLoc)->startYId * (thisSimLoc->thisMaterialLoc)->Nx;
	// create space for initial state (unpadded size) and fill in values
	int nx = (thisSimLoc->thisMaterialLoc)->Nx;
	int ny = (thisSimLoc->thisMaterialLoc)->NyLocal;
	int nPts = nx * ny;
	thisSimLoc->initStateLoc = malloc(nPts*sizeof(float));
	int i;
	for(i=0; i<nPts; ++i) thisSimLoc->initStateLoc[i] = valsForInitStateGlobal[startIdx + i];
	
	// check rank and number of processes
	int rank, nProcs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	
	// add boundary conditions on 0th and last columns (just in case valsForInitStateGlobal didn't follow the bdryVal)
	thisSimLoc->bdryVal = bdryVal;
	for(i=0; i<nPts; i+=nx) thisSimLoc->initStateLoc[i] = bdryVal; // column 0
	for(i=nx-1; i<nPts; i+=nx) thisSimLoc->initStateLoc[i] = bdryVal; // column nx-1
	// ===================================BEGIN STUDENT CODE==============================
	// if 0th or last rank, fill boundary conditions in 0th or last row

        // row 0 of 0th process
	
	// last row of last process
		

	// create padded state array for current state and fill in with initial state values
	
        // number of points of padding on at the start of local padded array
	
        // total number of points including padding on both sides
	
	// allocate the prior array
	// fill with values
	
	// create padded state array for current state
	
	// ===============================END OF STUDENT CODE==================================

	if((thisSimLoc->priorStateLoc == NULL) || (thisSimLoc->currentStateLoc == NULL)) flag = 1;
	// return a 0 if all was ok, but a 1 if there were issues
	return flag;
};

// Share ghost region information (must be done before each step of the simulation). Send first
// row of unpadded part of local state to previous process (except rank 0), and send last row
// of unpadded part of local state to next process (except last rank). 
// Need to get these ghost regions filled in into the priorStateLoc (so they can be used for next computation).
int exchangeGhostRegions(simLoc *thisSimLoc){
    
    // ====================BEGIN STUDENT CODE===========================================
    
    // check rank and number of processes
    
    // calculate previous and next rank 
     // note, will ignore previous for rank == 0
     // note will ignore next for last rank, size-1
    
    // pointer to location to start filling in data received from previous rank (start of padded array)
    
    // spot to start sending data to previous rank (start of unpadded array, # of padding points inside padded array)
    
    // spot to start sending data to next rank (# padding points inside the end of the unpadded array, 2* # padding points indside end of padded  array)
    
    // spot to start filling in data from next rank (# padding points inside end of padded array)
    
    // calculate number of messages to send/recv and set up statuses and requests

    // counter to track message ID

    // interactions with previous rank if rank > 0
      
    // receive from previous (up) rank 
      
    // send to previous (up) rank
    
    // otherwise fill appropriate padding with the boundary value

    // interactions with next rank if not the last rank
      
    // receive from next (down) rank
      
    // send to next (down)rank
    
    // otherwise fill appropriate padding with the boundary values

    // make sure all the requests are done
    
    // ====================END STUDENT CODE===========================================

    return 0;
};

// Share ghost regions, then move the simulation forward by one time step
int oneStepLoc(simLoc *thisSimLoc){
    int flag = 0;
	// grab the prior state and current (i.e. to update) state
	float *newStateLoc = thisSimLoc->currentStateLoc;
	float *priorStateLoc = thisSimLoc->priorStateLoc; 

	// if you run into problems return a 1
	if((newStateLoc == NULL) || (priorStateLoc == NULL)){ 
		printf("WARNING: null pointer for state encountered in oneStep() \n");
		return 1;
	}
	
	// shuffle around ghost region information into the priorStateLoc padded regions
	flag += exchangeGhostRegions(thisSimLoc);

	// grab the dimensions and material properties
	int nCols  = (thisSimLoc->thisMaterialLoc)->Nx;
	int nRowsTotal  = (thisSimLoc->thisMaterialLoc)->NyPadded;
	int nRowsUnpadded = (thisSimLoc->thisMaterialLoc)->NyLocal;
	int nPadRows = (thisSimLoc->thisMaterialLoc)->nPadRows;
	int nRowsGlobal = (thisSimLoc->thisMaterialLoc)->NyTotal;
	float dx = (thisSimLoc->thisMaterialLoc)->dx;
	float dy = (thisSimLoc->thisMaterialLoc)->dy;
	float alpha = (thisSimLoc->thisMaterialLoc)->alpha;

	// ============================BEGIN STUDENT CODE==============================
	// go one entry at a time filling in newStateLoc based on the values in priorStateLoc

		    // calculate the row index of this row in the global array
		    
		    // index of current location in padded subarray		    
			
			// case for all points not on the boundaries is to use the 5 pt stencil to update

				// calculate d^2/dx^2 term
				  // index of point to left in padded subarray
				  // index of point to right in padded subarray

				// calculate d^2/dy^2 term
				  // index of point on row above in padded subarray
				  // index of point on row below in padded subarray

            
             // case for all points on the boundaries is to fill with boundary value

	// =======================END STUDENT CODE========================================

	// Now that the new state is all updated and the prior state is no longer needed,
	// copy the values in the unpadded part of newState into priorState, and move onto the next time step.
	thisSimLoc->currentTimeIdx = thisSimLoc->currentTimeIdx + 1; // have completed a time step
	for(row=nPadRows; row<nRowsUnpadded+nPadRows; ++row){
		for(col=0; col<nCols; ++col){
			int idx = (row*nCols) + col;
			priorStateLoc[idx] = newStateLoc[idx];
		}
	}

	// since no problems were found earlier, return a 0
	return 0;

};

// Simulate nSteps time steps and record snapshots of the whole temperature field
// every stepsPerCheckPt time steps. runSimLoc does do the initialization of the checkPtTimeLoc
// struct automatically at the beginning of the simulation. 
// Note: running the simulation doesn't also initialize the sim or the material. Do them separately.
int runSimLoc(simLoc *thisSimLoc, int nSteps, int stepsPerCheckPt, checkPtTimeLoc *theseTimesLoc){
	int flag = 0; // return flag, 0 if no problem, but nonzero if there's a problem

	// do the initialization of the checkpointing struct
	int nSnaps = calcNSnapsLoc(nSteps,stepsPerCheckPt);
	int checkPtInitFlag = initCheckPtTimeLoc(theseTimesLoc, thisSimLoc->thisMaterialLoc, thisSimLoc, nSnaps);
	if(checkPtInitFlag){
		printf("WARNING: issue initializing checkpoint in runSim \n");
		flag =  checkPtInitFlag;
	}
	
	// check rank
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// run through the steps
	int step;
	recordSnapLoc(theseTimesLoc); // always record 0th  time step's prior state
	for(step=1; step<nSteps; ++step){
	    // share ghost regions and update simulation
		int stepFlag = oneStepLoc(thisSimLoc);
		if(stepFlag){
			printf("WARNING: issue in simulation at %d time step on rank %d \n",step, rank);
			flag = stepFlag;
		}
		if(step % stepsPerCheckPt == 0) recordSnapLoc(theseTimesLoc); // this checks if a checkpoint needs to be recorded, and records it if needed
	}
	return flag;

};

// deallocate memory associated with currentStateLoc, initStateLoc, priorStateLoc
int cleanupSimLoc(simLoc *thisSimLoc){
	free(thisSimLoc->initStateLoc);
	thisSimLoc->initStateLoc = NULL;
	free(thisSimLoc->priorStateLoc);
	thisSimLoc->priorStateLoc = NULL;
	free(thisSimLoc->currentStateLoc);
	thisSimLoc->currentStateLoc = NULL;
	return 0;
};
