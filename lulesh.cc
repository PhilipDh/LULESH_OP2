/*
  This is a Version 2.0 MPI + OpenMP implementation of LULESH

                 Copyright (c) 2010-2013.
      Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.
                  LLNL-CODE-461231
                All rights reserved.

This file is part of LULESH, Version 2.0.
Please also read this link -- http://www.opensource.org/licenses/index.php

//////////////
DIFFERENCES BETWEEN THIS VERSION (2.x) AND EARLIER VERSIONS:
* Addition of regions to make work more representative of multi-material codes
* Default size of each domain is 30^3 (27000 elem) instead of 45^3. This is
  more representative of our actual working set sizes
* Single source distribution supports pure serial, pure OpenMP, MPI-only, 
  and MPI+OpenMP
* Addition of ability to visualize the mesh using VisIt 
  https://wci.llnl.gov/codes/visit/download.html
* Various command line options (see ./lulesh2.0 -h)
 -q              : quiet mode - suppress stdout
 -i <iterations> : number of cycles to run
 -s <size>       : length of cube mesh along side
 -r <numregions> : Number of distinct regions (def: 11)
 -b <balance>    : Load balance between regions of a domain (def: 1)
 -c <cost>       : Extra cost of more expensive regions (def: 1)
 -f <filepieces> : Number of file parts for viz output (def: np/9)
 -p              : Print out progress
 -v              : Output viz file (requires compiling with -DVIZ_MESH
 -h              : This message

 printf("Usage: %s [opts]\n", execname);
      printf(" where [opts] is one or more of:\n");
      printf(" -q              : quiet mode - suppress all stdout\n");
      printf(" -i <iterations> : number of cycles to run\n");
      printf(" -s <size>       : length of cube mesh along side\n");
      printf(" -r <numregions> : Number of distinct regions (def: 11)\n");
      printf(" -b <balance>    : Load balance between regions of a domain (def: 1)\n");
      printf(" -c <cost>       : Extra cost of more expensive regions (def: 1)\n");
      printf(" -f <numfiles>   : Number of files to split viz dump into (def: (np+10)/9)\n");
      printf(" -p              : Print out progress\n");
      printf(" -v              : Output viz file (requires compiling with -DVIZ_MESH\n");
      printf(" -h              : This message\n");
      printf("\n\n");

*Notable changes in LULESH 2.0

* Split functionality into different files
lulesh.cc - where most (all?) of the timed functionality lies
lulesh-comm.cc - MPI functionality
lulesh-init.cc - Setup code
lulesh-viz.cc  - Support for visualization option
lulesh-util.cc - Non-timed functions
*
* The concept of "regions" was added, although every region is the same ideal
*    gas material, and the same sedov blast wave problem is still the only
*    problem its hardcoded to solve.
* Regions allow two things important to making this proxy app more representative:
*   Four of the LULESH routines are now performed on a region-by-region basis,
*     making the memory access patterns non-unit stride
*   Artificial load imbalances can be easily introduced that could impact
*     parallelization strategies.  
* The load balance flag changes region assignment.  Region number is raised to
*   the power entered for assignment probability.  Most likely regions changes
*   with MPI process id.
* The cost flag raises the cost of ~45% of the regions to evaluate EOS by the
*   entered multiple. The cost of 5% is 10x the entered multiple.
* MPI and OpenMP were added, and coalesced into a single version of the source
*   that can support serial builds, MPI-only, OpenMP-only, and MPI+OpenMP
* Added support to write plot files using "poor mans parallel I/O" when linked
*   with the silo library, which in turn can be read by VisIt.
* Enabled variable timestep calculation by default (courant condition), which
*   results in an additional reduction.
* Default domain (mesh) size reduced from 45^3 to 30^3
* Command line options to allow numerous test cases without needing to recompile
* Performance optimizations and code cleanup beyond LULESH 1.0
* Added a "Figure of Merit" calculation (elements solved per microsecond) and
*   output in support of using LULESH 2.0 for the 2017 CORAL procurement
*
* Possible Differences in Final Release (other changes possible)
*
* High Level mesh structure to allow data structure transformations
* Different default parameters
* Minor code performance changes and cleanup

TODO in future versions
* Add reader for (truly) unstructured meshes, probably serial only
* CMake based build system

//////////////

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the disclaimer below.

   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the disclaimer (as noted below)
     in the documentation and/or other materials provided with the
     distribution.

   * Neither the name of the LLNS/LLNL nor the names of its contributors
     may be used to endorse or promote products derived from this software
     without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Additional BSD Notice

1. This notice is required to be provided under our contract with the U.S.
   Department of Energy (DOE). This work was produced at Lawrence Livermore
   National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.

2. Neither the United States Government nor Lawrence Livermore National
   Security, LLC nor any of their employees, makes any warranty, express
   or implied, or assumes any liability or responsibility for the accuracy,
   completeness, or usefulness of any information, apparatus, product, or
   process disclosed, or represents that its use would not infringe
   privately-owned rights.

3. Also, reference herein to any specific commercial products, process, or
   services by trade name, trademark, manufacturer or otherwise does not
   necessarily constitute or imply its endorsement, recommendation, or
   favoring by the United States Government or Lawrence Livermore National
   Security, LLC. The views and opinions of authors expressed herein do not
   necessarily state or reflect those of the United States Government or
   Lawrence Livermore National Security, LLC, and shall not be used for
   advertising or product endorsement purposes.

*/

#include <climits>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <unistd.h>
#include <iomanip>

#include <op_seq.h>

#if _OPENMP
# include <omp.h>
#endif


#include "lulesh.h"
const Real_t  m_hgcoef = Real_t(3.0);            // hourglass control 
// Cutoffs (treat as constants)
const Real_t  m_e_cut = Real_t(1.0e-7);             // energy tolerance 
const Real_t  m_p_cut = Real_t(1.0e-7) ;             // pressure tolerance 
const Real_t  m_q_cut = Real_t(1.0e-7) ;             // q tolerance 
const Real_t  m_v_cut = Real_t(1.0e-10) ;             // relative volume tolerance 
const Real_t  m_u_cut = Real_t(1.0e-7) ;             // velocity tolerance 

// Other constants (usually setable, but hardcoded in this proxy app)


const Real_t  m_ss4o3 = (Real_t(4.0)/Real_t(3.0));
const Real_t  m_qstop = Real_t(1.0e+12);             // excessive q indicator 
const Real_t  m_monoq_max_slope = Real_t(1.0);
const Real_t  m_monoq_limiter_mult = Real_t(2.0);
const Real_t  m_qlc_monoq = Real_t(0.5);         // linear term coef for q 
const Real_t  m_qqc_monoq = (Real_t(2.0)/Real_t(3.0));         // quadratic term coef for q 
const Real_t  m_qqc = Real_t(2.0);
const Real_t  m_eosvmax = Real_t(1.0e+9);
const Real_t  m_eosvmin = Real_t(1.0e-9);
const Real_t  m_pmin = Real_t(0.);              // pressure floor 
const Real_t  m_emin = Real_t(-1.0e+15);              // energy floor 
const Real_t  m_dvovmax = Real_t(0.1);           // maximum allowable volume change 
const Real_t  m_refdens = Real_t(1.0);           // reference density 

#include "CalcForceForNodes.h"
#include "CalcAccelerationForNodes.h"
#include "CalcVelocityForNodes.h"
#include "CalcPositionForNodes.h"
#include "CalcKinematicsForElems.h"
#include "CalcLagrangeElements.h"
#include "CalcMonotonicQGradientsForElems.h"
#include "UpdateVolumesForElems.h"
#include "NoExcessiveArtificialViscosity.h"
// Define arrays and constants
// Should maybe be moved to a header file to avoid clutter or to the main function to not be global
// Element-centered

#define DEBUGGING 0

   // Region information
Int_t    m_numReg ;
Int_t    m_cost; //imbalance cost
Index_t *m_regElemSize ;   // Size of region sets
Index_t *m_regNumList ;    // Region number per domain element
Index_t **m_regElemlist ;  // region indexset 

Index_t* nodelist;

Index_t* lxim; /* element connectivity across each face */
Index_t* lxip;
Index_t* letam;
Index_t* letap;
Index_t* lzetam;
Index_t* lzetap;

Int_t* elemBC;            /* symmetry/free-surface flags for each elem face */

Real_t* dxx ;  /* principal strains -- temporary */
Real_t* dyy ;
Real_t* dzz ;

Real_t* delv_xi ;    /* velocity gradient -- temporary */
Real_t* delv_eta ;
Real_t* delv_zeta ;

Real_t* delx_xi ;    /* coordinate gradient -- temporary */
Real_t* delx_eta ;
Real_t* delx_zeta ;

Real_t* e;
Real_t* p;

Real_t* q;
Real_t* ql;
Real_t* qq;

Real_t* v;     /* relative volume */        
Real_t* vnew;  /* reference volume */
Real_t* volo;  /* new relative volume -- temporary */
Real_t* delv;  /* m_vnew - m_v */
Real_t* m_vdov;  /* volume derivative over volume */

Real_t* arealg; /* characteristic length of an element */

Real_t* ss;

Real_t* elemMass;

/* Node-centered */
Real_t* x;                  /* coordinates */
Real_t* y;
Real_t* z;

Real_t* xd;                 /* velocities */
Real_t* yd;
Real_t* zd;

Real_t* xdd;              /* accelerations */
Real_t* ydd;
Real_t* zdd;

Real_t* m_fx;                 /* forces */
Real_t* m_fy;
Real_t* m_fz;

Real_t* nodalMass;          /* mass */

Index_t* symmX;              /* symmetry plane nodesets */
Index_t* symmY;
Index_t* symmZ;

/* Other Constants */

// Variables to keep track of timestep, simulation time, and cycle
Real_t  m_dtcourant ;         // courant constraint 
Real_t  m_dthydro ;           // volume change constraint 
Int_t   m_cycle ;             // iteration count for simulation 
Real_t  m_dtfixed ;           // fixed time increment 
Real_t  m_time ;              // current time 
Real_t  m_deltatime ;         // variable time increment 
Real_t  m_deltatimemultlb ;
Real_t  m_deltatimemultub ;
Real_t  m_dtmax ;             // maximum allowable time increment 
Real_t  m_stoptime ;          // end time for simulation 

Index_t sizeX ;
Index_t sizeY ;
Index_t sizeZ ;
Index_t m_numElem ;
Index_t m_numNode ;


// OMP hack 
Index_t* m_nodeElemStart ;
Index_t* m_nodeElemCornerList ;

// Used in setup
Index_t rowMin, rowMax;
Index_t colMin, colMax;
Index_t planeMin, planeMax ;

// OP2 Variables
op_set nodes;
op_set elems;
op_set symmetry;
op_set temp_vols;

op_map p_nodelist;
op_map p_symmX;
op_map p_symmY;
op_map p_symmZ;

op_dat p_fx;
op_dat p_fy;
op_dat p_fz;

op_dat p_xd;
op_dat p_yd;
op_dat p_zd;

op_dat p_xdd;
op_dat p_ydd;
op_dat p_zdd;

op_dat p_x;
op_dat p_y;
op_dat p_z;

op_dat p_p;
op_dat p_q;
op_dat p_v;
op_dat p_volo;
op_dat p_arealg;
op_dat p_delv;
op_dat p_vnew;
op_dat p_vdov;

op_dat p_ss;
op_dat p_elemMass;
op_dat p_nodalMass;

op_dat p_delv_xi ;    /* velocity gradient -- temporary */
op_dat p_delv_eta ;
op_dat p_delv_zeta ;

op_dat p_delx_xi ;    /* coordinate gradient -- temporary */
op_dat p_delx_eta ;
op_dat p_delx_zeta ;

//Temporary OP2 Vars
Real_t *sigxx; 
Real_t *sigyy; 
Real_t *sigzz; 
Real_t *determ;
op_dat p_sigxx;
op_dat p_determ;

Real_t* dvdx;
Real_t* dvdy;
Real_t* dvdz;
Real_t* x8n;
Real_t* y8n;
Real_t* z8n;
op_dat p_dvdx;
op_dat p_dvdy;
op_dat p_dvdz;
op_dat p_x8n;
op_dat p_y8n;
op_dat p_z8n;


op_dat p_dxx;
op_dat p_dyy;
op_dat p_dzz;

static inline 
void AllocStrains(Int_t numElem){
   dxx = (Real_t*) malloc(numElem * sizeof(Real_t));
   dyy = (Real_t*) malloc(numElem * sizeof(Real_t));
   dzz = (Real_t*) malloc(numElem * sizeof(Real_t));
}

static inline
void AllocGradients(Int_t numElem, Int_t allElem){
   delx_xi = (Real_t*) malloc(numElem * sizeof(Real_t));
   delx_eta = (Real_t*) malloc(numElem * sizeof(Real_t));
   delx_zeta = (Real_t*) malloc(numElem * sizeof(Real_t));

   delv_xi = (Real_t*) malloc(allElem * sizeof(Real_t));
   delv_eta = (Real_t*) malloc(allElem * sizeof(Real_t));
   delv_zeta = (Real_t*) malloc(allElem * sizeof(Real_t));
}

static inline
void DeallocStrains(){
   Release(&dzz);
   Release(&dyy);
   Release(&dxx);
}

static inline
void DeallocGradients(){
   Release(&delx_zeta);
   Release(&delx_eta) ;
   Release(&delx_xi)  ;

   Release(&delv_zeta);
   Release(&delv_eta) ;
   Release(&delv_xi)  ;
}

static inline
void allocateElems(Int_t edgeElems){
   Int_t numElem = edgeElems * edgeElems * edgeElems;

   nodelist = (Index_t*) malloc(numElem*8 * sizeof(Index_t));

   lxim = (Index_t*) malloc(numElem * sizeof(Index_t));
   lxip = (Index_t*) malloc(numElem * sizeof(Index_t));
   letam = (Index_t*) malloc(numElem * sizeof(Index_t));
   letap = (Index_t*) malloc(numElem * sizeof(Index_t));
   lzetam = (Index_t*) malloc(numElem * sizeof(Index_t));
   lzetap = (Index_t*) malloc(numElem * sizeof(Index_t));

   elemBC = (Int_t*) malloc(numElem * sizeof(Int_t));

   e = (Real_t*) malloc(numElem * sizeof(Real_t));
   p = (Real_t*) malloc(numElem * sizeof(Real_t));

   q = (Real_t*) malloc(numElem * sizeof(Real_t));
   ql = (Real_t*) malloc(numElem * sizeof(Real_t));
   qq = (Real_t*) malloc(numElem * sizeof(Real_t));

   v = (Real_t*) malloc(numElem * sizeof(Real_t));

   volo = (Real_t*) malloc(numElem * sizeof(Real_t));
   delv = (Real_t*) malloc(numElem * sizeof(Real_t));
   m_vdov = (Real_t*) malloc(numElem * sizeof(Real_t));

   arealg = (Real_t*) malloc(numElem * sizeof(Real_t));

   ss = (Real_t*) malloc(numElem * sizeof(Real_t));

   elemMass = (Real_t*) malloc(numElem * sizeof(Real_t));

   vnew = (Real_t*) malloc(numElem * sizeof(Real_t));

   sigxx = (Real_t*) malloc(numElem* 3 * sizeof(Real_t));
   determ = (Real_t*) malloc(numElem * sizeof(Real_t));

   dvdx = (Real_t*) malloc(numElem * 8 * sizeof(Real_t));
   dvdy = (Real_t*) malloc(numElem * 8 * sizeof(Real_t));
   dvdz = (Real_t*) malloc(numElem * 8 * sizeof(Real_t));
   x8n = (Real_t*) malloc(numElem * 8 * sizeof(Real_t));
   y8n = (Real_t*) malloc(numElem * 8 * sizeof(Real_t));
   z8n = (Real_t*) malloc(numElem * 8 * sizeof(Real_t));

   dxx = (Real_t*) malloc(numElem * sizeof(Real_t));
   dyy = (Real_t*) malloc(numElem * sizeof(Real_t));
   dzz = (Real_t*) malloc(numElem * sizeof(Real_t));

   delx_xi = (Real_t*) malloc(numElem * sizeof(Real_t));
   delx_eta = (Real_t*) malloc(numElem * sizeof(Real_t));
   delx_zeta = (Real_t*) malloc(numElem * sizeof(Real_t));

   delv_xi = (Real_t*) malloc(numElem * sizeof(Real_t));
   delv_eta = (Real_t*) malloc(numElem * sizeof(Real_t));
   delv_zeta = (Real_t*) malloc(numElem * sizeof(Real_t));
}

static inline
void allocateNodes(Index_t colLoc, Index_t rowLoc, Index_t planeLoc,
                  Int_t edgeNodes){

   Int_t numNode = edgeNodes * edgeNodes * edgeNodes;

   x = (Real_t*) malloc(numNode * sizeof(Real_t)); // Coordinates
   y = (Real_t*) malloc(numNode * sizeof(Real_t));
   z = (Real_t*) malloc(numNode * sizeof(Real_t));

   xd = (Real_t*) malloc(numNode * sizeof(Real_t)); // Velocities
   yd = (Real_t*) malloc(numNode * sizeof(Real_t));
   zd = (Real_t*) malloc(numNode * sizeof(Real_t));

   xdd = (Real_t*) malloc(numNode * sizeof(Real_t)); //Accelerations
   // xdd = (Real_t*) malloc(numNode * 3 * sizeof(Real_t)); //Accelerations
   ydd = (Real_t*) malloc(numNode * sizeof(Real_t));
   zdd = (Real_t*) malloc(numNode * sizeof(Real_t));

   m_fx = (Real_t*) malloc(numNode * sizeof(Real_t)); // Forces
   // m_fx = (Real_t*) malloc(numNode * 3 * sizeof(Real_t)); // Forces
   m_fy = (Real_t*) malloc(numNode * sizeof(Real_t));
   m_fz = (Real_t*) malloc(numNode * sizeof(Real_t));


   nodalMass = (Real_t*) malloc(numNode * sizeof(Real_t));  // mass

        // Boundary nodesets
   if (colLoc == 0)
      symmX = (Index_t*) malloc(edgeNodes*edgeNodes*sizeof(Index_t));
   if (rowLoc == 0)
      symmY = (Index_t*) malloc(edgeNodes*edgeNodes*sizeof(Index_t));
   if (planeLoc == 0)
      symmZ = (Index_t*) malloc(edgeNodes*edgeNodes*sizeof(Index_t));
}

static inline
void initialise(Index_t colLoc,
               Index_t rowLoc, Index_t planeLoc,
               Index_t nx, Int_t tp, Int_t nr, Int_t balance, Int_t cost){
   Index_t edgeElems = nx ;
   Index_t edgeNodes = edgeElems+1 ;

   sizeX = edgeElems ;
   sizeY = edgeElems ;
   sizeZ = edgeElems ;
   m_numElem = edgeElems * edgeElems * edgeElems;

   m_numNode = edgeNodes * edgeNodes * edgeNodes;

   m_regNumList = (Index_t*) malloc(m_numElem * sizeof(Index_t));

   //! Setup Comm Buffer, Should not be necessary in final app
   rowMin = (rowLoc == 0)        ? 0 : 1;
   rowMax = (rowLoc == tp-1)     ? 0 : 1;
   colMin = (colLoc == 0)        ? 0 : 1;
   colMax = (colLoc == tp-1)     ? 0 : 1;
   planeMin = (planeLoc == 0)    ? 0 : 1;
   planeMax = (planeLoc == tp-1) ? 0 : 1;

   for(Index_t i=0; i<m_numElem;++i){
      p[i] = Real_t(0.0);
      e[i] = Real_t(0.0);
      q[i] = Real_t(0.0);
      ss[i] = Real_t(0.0);
   }

   for(Index_t i=0; i<m_numElem;++i){
      v[i] = Real_t(1.0);
   }

   for (Index_t i = 0; i<m_numNode;++i){
      xd[i] = Real_t(0.0);
      yd[i] = Real_t(0.0);
      zd[i] = Real_t(0.0);
   }

   for (Index_t i = 0; i<m_numNode;++i){
   // for (Index_t i = 0; i<m_numNode*3;++i){
      xdd[i] = Real_t(0.0);
      ydd[i] = Real_t(0.0);
      zdd[i] = Real_t(0.0);
   }

   for (Index_t i=0; i<m_numNode; ++i) {
      nodalMass[i] = Real_t(0.0) ;
   }

   //!Build Mesh function
   Index_t meshEdgeElems = tp * nx;
   Index_t nidx = 0;
   Real_t tz = Real_t(1.125)*Real_t(planeLoc*nx)/Real_t(meshEdgeElems) ;
   for (Index_t plane=0; plane<edgeNodes; ++plane) {
      Real_t ty = Real_t(1.125)*Real_t(rowLoc*nx)/Real_t(meshEdgeElems) ;
      for (Index_t row=0; row<edgeNodes; ++row) {
         Real_t tx = Real_t(1.125)*Real_t(colLoc*nx)/Real_t(meshEdgeElems) ;
         for (Index_t col=0; col<edgeNodes; ++col) {
         x[nidx] = tx ;
         y[nidx] = ty ;
         z[nidx] = tz ;
         ++nidx ;
         // tx += ds ; // may accumulate roundoff... 
         tx = Real_t(1.125)*Real_t(colLoc*nx+col+1)/Real_t(meshEdgeElems) ;
         }
         // ty += ds ;  // may accumulate roundoff... 
         ty = Real_t(1.125)*Real_t(rowLoc*nx+row+1)/Real_t(meshEdgeElems) ;
      }
      // tz += ds ;  // may accumulate roundoff... 
      tz = Real_t(1.125)*Real_t(planeLoc*nx+plane+1)/Real_t(meshEdgeElems) ;
   }

     // embed hexehedral elements in nodal point lattice 
   Index_t zidx = 0 ;
   nidx = 0 ;
   for (Index_t plane=0; plane<edgeElems; ++plane) {
      for (Index_t row=0; row<edgeElems; ++row) {
         for (Index_t col=0; col<edgeElems; ++col) {
         Index_t *localNode = &nodelist[zidx*Index_t(8)] ; //!TODO Check original implementation for reference

         localNode[0] = nidx                                       ;
         localNode[1] = nidx                                   + 1 ;
         localNode[2] = nidx                       + edgeNodes + 1 ;
         localNode[3] = nidx                       + edgeNodes     ;
         localNode[4] = nidx + edgeNodes*edgeNodes                 ;
         localNode[5] = nidx + edgeNodes*edgeNodes             + 1 ;
         localNode[6] = nidx + edgeNodes*edgeNodes + edgeNodes + 1 ;
         localNode[7] = nidx + edgeNodes*edgeNodes + edgeNodes     ;
         ++zidx ;
         ++nidx ;
         }
         ++nidx ;
      }
      nidx += edgeNodes ;
   }
   //! End Build Mesh Function

   //! Start Create Region Sets
   srand(0);
   Index_t myRank = 0;

   m_numReg = nr;
   m_regElemSize = (Index_t*) malloc(m_numReg * sizeof(Index_t));
   m_regElemlist = (Index_t**) malloc(m_numReg * sizeof(Index_t));
   Index_t nextIndex = 0;
   //if we only have one region just fill it
   // Fill out the regNumList with material numbers, which are always
   // the region index plus one 
   if(m_numReg == 1){
      while(nextIndex < m_numElem){
         m_regNumList[nextIndex] = 1;
         nextIndex++;
      }
      m_regElemSize[0] = 0;
   } else {//If we have more than one region distribute the elements.
      Int_t regionNum;
      Int_t regionVar;
      Int_t lastReg = -1;
      Int_t binSize;
      Index_t elements;
      Index_t runto = 0;
      Int_t costDenominator = 0;
      Int_t* regBinEnd = (Int_t*) malloc(m_numReg * sizeof(Int_t));
      //Determine the relative weights of all the regions.  This is based off the -b flag.  Balance is the value passed into b.  
      for(Index_t i=0; i<m_numReg;++i){
         m_regElemSize[i] = 0;
         costDenominator += pow((i+1), balance);//Total sum of all regions weights
         regBinEnd[i] = costDenominator;  //Chance of hitting a given region is (regBinEnd[i] - regBinEdn[i-1])/costDenominator
      }
      //Until all elements are assigned
      while (nextIndex < m_numElem) {
         //pick the region
         regionVar = rand() % costDenominator;
         Index_t i = 0;
         while(regionVar >= regBinEnd[i]) i++;
         //rotate the regions based on MPI rank.  Rotation is Rank % m_numRegions this makes each domain have a different region with 
         //the highest representation
         regionNum = ((i+myRank)% m_numReg) +1;
         while(regionNum == lastReg){
            regionVar = rand() % costDenominator;
            i = 0;
            while(regionVar >= regBinEnd[i]) i++;
            regionNum = ((i + myRank) % m_numReg) + 1;
         }
         //Pick the bin size of the region and determine the number of elements.
         binSize = rand() % 1000;
         if(binSize < 773) {
	         elements = rand() % 15 + 1;
	      }
	      else if(binSize < 937) {
	         elements = rand() % 16 + 16;
	      }
	      else if(binSize < 970) {
	         elements = rand() % 32 + 32;
	      }
	      else if(binSize < 974) {
	         elements = rand() % 64 + 64;
	      } 
	      else if(binSize < 978) {
	         elements = rand() % 128 + 128;
	      }
	      else if(binSize < 981) {
	         elements = rand() % 256 + 256;
	      }
	      else
	         elements = rand() % 1537 + 512;
         runto = elements + nextIndex;
	      //Store the elements.  If we hit the end before we run out of elements then just stop.
         while (nextIndex < runto && nextIndex < m_numElem) {
	         m_regNumList[nextIndex] = regionNum;
	         nextIndex++;
	      }
         lastReg = regionNum;
      }
      delete [] regBinEnd; 
   }
      // Convert m_regNumList to region index sets
   // First, count size of each region 
   for (Index_t i=0 ; i<m_numElem ; ++i) {
      int r = m_regNumList[i]-1; // region index == regnum-1
      m_regElemSize[r]++;
   }
   // Second, allocate each region index set
   for (Index_t i=0 ; i<m_numReg ; ++i) {
      m_regElemlist[i] = (Index_t*) malloc(m_regElemSize[i]*sizeof(Index_t));
      m_regElemSize[i] = 0;
   }
   // Third, fill index sets
   for (Index_t i=0 ; i<m_numElem ; ++i) {
      Index_t r = m_regNumList[i]-1;       // region index == regnum-1
      Index_t regndx = m_regElemSize[r]++; // Note increment
      m_regElemlist[r][regndx] = i;
   }

   //! End Create Region Sets


   //! Setup Symmetry Planes Function
   nidx = 0 ;
   for (Index_t i=0; i<edgeNodes; ++i) {
      Index_t planeInc = i*edgeNodes*edgeNodes ;
      Index_t rowInc   = i*edgeNodes ;
      for (Index_t j=0; j<edgeNodes; ++j) {
         if (planeLoc == 0) {
            symmZ[nidx] = rowInc   + j ;
         }
         if (rowLoc == 0) {
            symmY[nidx] = planeInc + j ;
         }
         if (colLoc == 0) {
            symmX[nidx] = planeInc + j*edgeNodes ;
         }
         ++nidx ;
      }
   }
   //! End Setup Symmetry Plnes Function

   //! Setup Elem Connectivity Function
   lxim[0] = 0 ;
   for (Index_t i=1; i<m_numElem; ++i) {
      lxim[i]   = i-1 ;
      lxip[i-1] = i ;
   }
   lxip[m_numElem-1] = m_numElem-1 ;

   for (Index_t i=0; i<edgeElems; ++i) {
      letam[i] = i ; 
      letap[m_numElem-edgeElems+i] = m_numElem-edgeElems+i ;
   }
   for (Index_t i=edgeElems; i<m_numElem; ++i) {
      letam[i] = i-edgeElems ;
      letap[i-edgeElems] = i ;
   }

   for (Index_t i=0; i<edgeElems*edgeElems; ++i) {
      lzetam[i] = i ;
      lzetap[m_numElem-edgeElems*edgeElems+i] = m_numElem-edgeElems*edgeElems+i ;
   }
   for (Index_t i=edgeElems*edgeElems; i<m_numElem; ++i) {
      lzetam[i] = i - edgeElems*edgeElems ;
      lzetap[i-edgeElems*edgeElems] = i ;
   }
   //! End Selem Connectivity Function

   //! Setup Boundary Conditions Function
   // set up boundary condition information
   Index_t ghostIdx[6] ;  // offsets to ghost locations
   for (Index_t i=0; i<m_numElem; ++i) {
      elemBC[i] = Int_t(0) ;
   }
     for (Index_t i=0; i<6; ++i) {
    ghostIdx[i] = INT_MIN ;
  }

  Int_t pidx = m_numElem ;
  if (planeMin != 0) {
    ghostIdx[0] = pidx ;
    pidx += sizeX*sizeY ;
  }

  if (planeMax != 0) {
    ghostIdx[1] = pidx ;
    pidx += sizeX*sizeY ;
  }

  if (rowMin != 0) {
    ghostIdx[2] = pidx ;
    pidx += sizeX*sizeZ ;
  }

  if (rowMax != 0) {
    ghostIdx[3] = pidx ;
    pidx += sizeX*sizeZ ;
  }

  if (colMin != 0) {
    ghostIdx[4] = pidx ;
    pidx += sizeY*sizeZ ;
  }

  if (colMax != 0) {
    ghostIdx[5] = pidx ;
  }
  // symmetry plane or free surface BCs 
  for (Index_t i=0; i<edgeElems; ++i) {

   Index_t planeInc = i*edgeElems*edgeElems ;
   Index_t rowInc   = i*edgeElems ;
   for (Index_t j=0; j<edgeElems; ++j) {
      if (planeLoc == 0) {
	      elemBC[rowInc+j] |= ZETA_M_SYMM ;
      }
      else {
	      elemBC[rowInc+j] |= ZETA_M_COMM ;
	      lzetam[rowInc+j] = ghostIdx[0] + rowInc + j ;
      }
      if (planeLoc == tp-1) {
	      elemBC[rowInc+j+m_numElem-edgeElems*edgeElems] |= ZETA_P_FREE;
      }
      else {
	      elemBC[rowInc+j+m_numElem-edgeElems*edgeElems] |= ZETA_P_COMM ;
	      lzetap[rowInc+j+m_numElem-edgeElems*edgeElems] = ghostIdx[1] + rowInc + j ;
      }
      if (rowLoc == 0) {
	      elemBC[planeInc+j] |= ETA_M_SYMM ;
      }
      else {
	      elemBC[planeInc+j] |= ETA_M_COMM ;
	      letam[planeInc+j] = ghostIdx[2] + rowInc + j ;
      }
      if (rowLoc == tp-1) {
	      elemBC[planeInc+j+edgeElems*edgeElems-edgeElems] |= ETA_P_FREE ;
      }
      else {
	      elemBC[planeInc+j+edgeElems*edgeElems-edgeElems] |=  ETA_P_COMM ;
	      letap[planeInc+j+edgeElems*edgeElems-edgeElems] = ghostIdx[3] +  rowInc + j ;
      }
      if (colLoc == 0) {
	      elemBC[planeInc+j*edgeElems] |= XI_M_SYMM ;
      }
      else {
	      elemBC[planeInc+j*edgeElems] |= XI_M_COMM ;
	      lxim[planeInc+j*edgeElems] = ghostIdx[4] + rowInc + j ;
      }

      if (colLoc == tp-1) {
	      elemBC[planeInc+j*edgeElems+edgeElems-1] |= XI_P_FREE ;
      }
      else {
	      elemBC[planeInc+j*edgeElems+edgeElems-1] |= XI_P_COMM ;
	      lxip[planeInc+j*edgeElems+edgeElems-1] = ghostIdx[5] + rowInc + j ;
      }
    }
   }
   //! End SetupBC Function

   // Setup defaults

   // These can be changed (requires recompile) if you want to run
   // with a fixed timestep, or to a different end time, but it's
   // probably easier/better to just run a fixed number of timesteps
   // using the -i flag in 2.x

   m_dtfixed = Real_t(-1.0e-6) ; // Negative means use courant condition
   m_stoptime  = Real_t(1.0e-2); // *Real_t(edgeElems*tp/45.0) ;

   // Initial conditions
   m_deltatimemultlb = Real_t(1.1) ;
   m_deltatimemultub = Real_t(1.2) ;
   m_dtcourant = Real_t(1.0e+20) ;
   m_dthydro   = Real_t(1.0e+20) ;
   m_dtmax     = Real_t(1.0e-2) ;
   m_time    = Real_t(0.) ;
   m_cycle   = Int_t(0) ;

   // initialize field data 
   for (Index_t i=0; i<m_numElem; ++i) {
      Real_t x_local[8], y_local[8], z_local[8] ;
      Index_t *elemToNode = &nodelist[i*Index_t(8)] ;
      for( Index_t lnode=0 ; lnode<8 ; ++lnode )
      {
        Index_t gnode = elemToNode[lnode];
        x_local[lnode] = x[gnode];
        y_local[lnode] = y[gnode];
        z_local[lnode] = z[gnode];
      }

      // volume calculations
      Real_t volume = CalcElemVolume(x_local, y_local, z_local );
      volo[i] = volume ;
      elemMass[i] = volume ;
      for (Index_t j=0; j<8; ++j) {
         Index_t idx = elemToNode[j] ;
         nodalMass[idx] += volume / Real_t(8.0) ;
      }
   }

   // deposit initial energy
   // An energy of 3.948746e+7 is correct for a problem with
   // 45 zones along a side - we need to scale it
   const Real_t ebase = Real_t(3.948746e+7);
   Real_t scale = (nx * tp)/Real_t(45.0);
   Real_t einit = ebase*scale*scale*scale;
   if (rowLoc + colLoc + planeLoc == 0) {
      // Dump into the first zone (which we know is in the corner)
      // of the domain that sits at the origin
      e[0] = einit;
   }
   m_deltatime = (Real_t(.5)*cbrt(volo[0]))/sqrt(Real_t(2.0)*einit);

}

/* Work Routines */

static inline
void TimeIncrement(Domain& domain)
{
   // Real_t targetdt = domain.stoptime() - domain.time() ;
   Real_t targetdt = m_stoptime - m_time ;

   // if ((domain.dtfixed() <= Real_t(0.0)) && (domain.cycle() != Int_t(0))) {
   if ((m_dtfixed <= Real_t(0.0)) && (m_cycle != Int_t(0))) {
   // if ((dtfixed <= Real_t(0.0)) && (cycle != Int_t(0))) {
      Real_t ratio ;
      // Real_t olddt = domain.deltatime() ;
      Real_t olddt = m_deltatime ;

      /* This will require a reduction in parallel */
      Real_t gnewdt = Real_t(1.0e+20) ;
      Real_t newdt ;
      if (m_dtcourant < gnewdt) {
         gnewdt = m_dtcourant / Real_t(2.0) ;
      }
      if (m_dthydro < gnewdt) {
         gnewdt = m_dthydro * Real_t(2.0) / Real_t(3.0) ;
      }
#if USE_MPI      
      MPI_Allreduce(&gnewdt, &newdt, 1,
                    ((sizeof(Real_t) == 4) ? MPI_FLOAT : MPI_DOUBLE),
                    MPI_MIN, MPI_COMM_WORLD) ;
#else
      //op_par_loop(timeIncrem, "TimeIncrement", nodes, 
      //             op_arg_dat(&gnewdt, -1, OP_ID, 1, "double", OP_READ),
      //             op_arg_gbl(&newdt,   1,           "double", OP_MIN))
      newdt = gnewdt;
#endif
      
      ratio = newdt / olddt ;
      if (ratio >= Real_t(1.0)) {
         // if (ratio < domain.deltatimemultlb()) {
         if (ratio < m_deltatimemultlb) {
            newdt = olddt ;
         }
         // else if (ratio > domain.deltatimemultub()) {
         //    newdt = olddt*domain.deltatimemultub() ;
         // }
         else if (ratio > m_deltatimemultub) {
            newdt = olddt*m_deltatimemultub ;
         }
      }

      // if (newdt > domain.dtmax()) {
      //    newdt = domain.dtmax() ;
      // }
      if (newdt > m_dtmax) {
         newdt = m_dtmax ;
      }
      // domain.deltatime() = newdt ;
      m_deltatime = newdt ;
   }

   /* TRY TO PREVENT VERY SMALL SCALING ON THE NEXT CYCLE */
   // if ((targetdt > domain.deltatime()) &&
   //     (targetdt < (Real_t(4.0) * domain.deltatime() / Real_t(3.0))) ) {
   //    targetdt = Real_t(2.0) * domain.deltatime() / Real_t(3.0) ;
   // }

   // if (targetdt < domain.deltatime()) {
   //    domain.deltatime() = targetdt ;
   // }
   if ((targetdt > m_deltatime) &&
       (targetdt < (Real_t(4.0) * m_deltatime / Real_t(3.0))) ) {
      targetdt = Real_t(2.0) * m_deltatime / Real_t(3.0) ;
   }

   if (targetdt < m_deltatime) {
      m_deltatime = targetdt ;
   }

   // domain.time() += domain.deltatime() ;
   m_time += m_deltatime ;

   // ++domain.cycle() ;
   ++m_cycle ;
}

/******************************************/

static inline
void CollectDomainNodesToElemNodes(Domain &domain,
                                   const Index_t* elemToNode,
                                   Real_t elemX[8],
                                   Real_t elemY[8],
                                   Real_t elemZ[8])
{
   Index_t nd0i = elemToNode[0] ;
   Index_t nd1i = elemToNode[1] ;
   Index_t nd2i = elemToNode[2] ;
   Index_t nd3i = elemToNode[3] ;
   Index_t nd4i = elemToNode[4] ;
   Index_t nd5i = elemToNode[5] ;
   Index_t nd6i = elemToNode[6] ;
   Index_t nd7i = elemToNode[7] ;

   // elemX[0] = domain.x(nd0i);
   // elemX[1] = domain.x(nd1i);
   // elemX[2] = domain.x(nd2i);
   // elemX[3] = domain.x(nd3i);
   // elemX[4] = domain.x(nd4i);
   // elemX[5] = domain.x(nd5i);
   // elemX[6] = domain.x(nd6i);
   // elemX[7] = domain.x(nd7i);
   elemX[0] = x[nd0i];
   elemX[1] = x[nd1i];
   elemX[2] = x[nd2i];
   elemX[3] = x[nd3i];
   elemX[4] = x[nd4i];
   elemX[5] = x[nd5i];
   elemX[6] = x[nd6i];
   elemX[7] = x[nd7i];
   // std::cout << x[nd7i];

   // elemY[0] = domain.y(nd0i);
   // elemY[1] = domain.y(nd1i);
   // elemY[2] = domain.y(nd2i);
   // elemY[3] = domain.y(nd3i);
   // elemY[4] = domain.y(nd4i);
   // elemY[5] = domain.y(nd5i);
   // elemY[6] = domain.y(nd6i);
   // elemY[7] = domain.y(nd7i);
   elemY[0] = y[nd0i];
   elemY[1] = y[nd1i];
   elemY[2] = y[nd2i];
   elemY[3] = y[nd3i];
   elemY[4] = y[nd4i];
   elemY[5] = y[nd5i];
   elemY[6] = y[nd6i];
   elemY[7] = y[nd7i];

   // elemZ[0] = domain.z(nd0i);
   // elemZ[1] = domain.z(nd1i);
   // elemZ[2] = domain.z(nd2i);
   // elemZ[3] = domain.z(nd3i);
   // elemZ[4] = domain.z(nd4i);
   // elemZ[5] = domain.z(nd5i);
   // elemZ[6] = domain.z(nd6i);
   // elemZ[7] = domain.z(nd7i);
   elemZ[0] = z[nd0i];
   elemZ[1] = z[nd1i];
   elemZ[2] = z[nd2i];
   elemZ[3] = z[nd3i];
   elemZ[4] = z[nd4i];
   elemZ[5] = z[nd5i];
   elemZ[6] = z[nd6i];
   elemZ[7] = z[nd7i];

}

/******************************************/

static inline
void InitStressTermsForElems(Domain &domain,
                             Real_t *sigxx, Real_t *sigyy, Real_t *sigzz,
                             Index_t numElem)
{
   //
   // pull in the stresses appropriate to the hydro integration
   op_par_loop(initStressTerms, "initStressTerms", elems,
               op_arg_dat(p_sigxx, -1, OP_ID, 3, "double", OP_WRITE),
               op_arg_dat(p_q, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_p, -1, OP_ID, 1, "double", OP_READ));
// #pragma omp parallel for firstprivate(numElem)
//    for (Index_t i = 0 ; i < numElem ; ++i){
//       // sigxx[i] = sigyy[i] = sigzz[i] =  - domain.p(i) - domain.q(i) ;
//       sigxx[i] = sigyy[i] = sigzz[i] =  - p[i] - q[i] ;
//    }
}

/******************************************/

static inline
void CalcElemShapeFunctionDerivatives( Real_t const x[],
                                       Real_t const y[],
                                       Real_t const z[],
                                       Real_t b[][8],
                                       Real_t* const volume )
{
  const Real_t x0 = x[0] ;   const Real_t x1 = x[1] ;
  const Real_t x2 = x[2] ;   const Real_t x3 = x[3] ;
  const Real_t x4 = x[4] ;   const Real_t x5 = x[5] ;
  const Real_t x6 = x[6] ;   const Real_t x7 = x[7] ;

  const Real_t y0 = y[0] ;   const Real_t y1 = y[1] ;
  const Real_t y2 = y[2] ;   const Real_t y3 = y[3] ;
  const Real_t y4 = y[4] ;   const Real_t y5 = y[5] ;
  const Real_t y6 = y[6] ;   const Real_t y7 = y[7] ;

  const Real_t z0 = z[0] ;   const Real_t z1 = z[1] ;
  const Real_t z2 = z[2] ;   const Real_t z3 = z[3] ;
  const Real_t z4 = z[4] ;   const Real_t z5 = z[5] ;
  const Real_t z6 = z[6] ;   const Real_t z7 = z[7] ;

  Real_t fjxxi, fjxet, fjxze;
  Real_t fjyxi, fjyet, fjyze;
  Real_t fjzxi, fjzet, fjzze;
  Real_t cjxxi, cjxet, cjxze;
  Real_t cjyxi, cjyet, cjyze;
  Real_t cjzxi, cjzet, cjzze;

  fjxxi = Real_t(.125) * ( (x6-x0) + (x5-x3) - (x7-x1) - (x4-x2) );
  fjxet = Real_t(.125) * ( (x6-x0) - (x5-x3) + (x7-x1) - (x4-x2) );
  fjxze = Real_t(.125) * ( (x6-x0) + (x5-x3) + (x7-x1) + (x4-x2) );

  fjyxi = Real_t(.125) * ( (y6-y0) + (y5-y3) - (y7-y1) - (y4-y2) );
  fjyet = Real_t(.125) * ( (y6-y0) - (y5-y3) + (y7-y1) - (y4-y2) );
  fjyze = Real_t(.125) * ( (y6-y0) + (y5-y3) + (y7-y1) + (y4-y2) );

  fjzxi = Real_t(.125) * ( (z6-z0) + (z5-z3) - (z7-z1) - (z4-z2) );
  fjzet = Real_t(.125) * ( (z6-z0) - (z5-z3) + (z7-z1) - (z4-z2) );
  fjzze = Real_t(.125) * ( (z6-z0) + (z5-z3) + (z7-z1) + (z4-z2) );

  /* compute cofactors */
  cjxxi =    (fjyet * fjzze) - (fjzet * fjyze);
  cjxet =  - (fjyxi * fjzze) + (fjzxi * fjyze);
  cjxze =    (fjyxi * fjzet) - (fjzxi * fjyet);

  cjyxi =  - (fjxet * fjzze) + (fjzet * fjxze);
  cjyet =    (fjxxi * fjzze) - (fjzxi * fjxze);
  cjyze =  - (fjxxi * fjzet) + (fjzxi * fjxet);

  cjzxi =    (fjxet * fjyze) - (fjyet * fjxze);
  cjzet =  - (fjxxi * fjyze) + (fjyxi * fjxze);
  cjzze =    (fjxxi * fjyet) - (fjyxi * fjxet);

  /* calculate partials :
     this need only be done for l = 0,1,2,3   since , by symmetry ,
     (6,7,4,5) = - (0,1,2,3) .
  */
  b[0][0] =   -  cjxxi  -  cjxet  -  cjxze;
  b[0][1] =      cjxxi  -  cjxet  -  cjxze;
  b[0][2] =      cjxxi  +  cjxet  -  cjxze;
  b[0][3] =   -  cjxxi  +  cjxet  -  cjxze;
  b[0][4] = -b[0][2];
  b[0][5] = -b[0][3];
  b[0][6] = -b[0][0];
  b[0][7] = -b[0][1];

  b[1][0] =   -  cjyxi  -  cjyet  -  cjyze;
  b[1][1] =      cjyxi  -  cjyet  -  cjyze;
  b[1][2] =      cjyxi  +  cjyet  -  cjyze;
  b[1][3] =   -  cjyxi  +  cjyet  -  cjyze;
  b[1][4] = -b[1][2];
  b[1][5] = -b[1][3];
  b[1][6] = -b[1][0];
  b[1][7] = -b[1][1];

  b[2][0] =   -  cjzxi  -  cjzet  -  cjzze;
  b[2][1] =      cjzxi  -  cjzet  -  cjzze;
  b[2][2] =      cjzxi  +  cjzet  -  cjzze;
  b[2][3] =   -  cjzxi  +  cjzet  -  cjzze;
  b[2][4] = -b[2][2];
  b[2][5] = -b[2][3];
  b[2][6] = -b[2][0];
  b[2][7] = -b[2][1];

  /* calculate jacobian determinant (volume) */
  *volume = Real_t(8.) * ( fjxet * cjxet + fjyet * cjyet + fjzet * cjzet);
//   std::cout << "Volume: " << *volume << "\n";
}

/******************************************/

static inline
void SumElemFaceNormal(Real_t *normalX0, Real_t *normalY0, Real_t *normalZ0,
                       Real_t *normalX1, Real_t *normalY1, Real_t *normalZ1,
                       Real_t *normalX2, Real_t *normalY2, Real_t *normalZ2,
                       Real_t *normalX3, Real_t *normalY3, Real_t *normalZ3,
                       const Real_t x0, const Real_t y0, const Real_t z0,
                       const Real_t x1, const Real_t y1, const Real_t z1,
                       const Real_t x2, const Real_t y2, const Real_t z2,
                       const Real_t x3, const Real_t y3, const Real_t z3)
{
   Real_t bisectX0 = Real_t(0.5) * (x3 + x2 - x1 - x0);
   Real_t bisectY0 = Real_t(0.5) * (y3 + y2 - y1 - y0);
   Real_t bisectZ0 = Real_t(0.5) * (z3 + z2 - z1 - z0);
   Real_t bisectX1 = Real_t(0.5) * (x2 + x1 - x3 - x0);
   Real_t bisectY1 = Real_t(0.5) * (y2 + y1 - y3 - y0);
   Real_t bisectZ1 = Real_t(0.5) * (z2 + z1 - z3 - z0);
   Real_t areaX = Real_t(0.25) * (bisectY0 * bisectZ1 - bisectZ0 * bisectY1);
   Real_t areaY = Real_t(0.25) * (bisectZ0 * bisectX1 - bisectX0 * bisectZ1);
   Real_t areaZ = Real_t(0.25) * (bisectX0 * bisectY1 - bisectY0 * bisectX1);

   *normalX0 += areaX;
   *normalX1 += areaX;
   *normalX2 += areaX;
   *normalX3 += areaX;

   *normalY0 += areaY;
   *normalY1 += areaY;
   *normalY2 += areaY;
   *normalY3 += areaY;

   *normalZ0 += areaZ;
   *normalZ1 += areaZ;
   *normalZ2 += areaZ;
   *normalZ3 += areaZ;
}

/******************************************/

static inline
void CalcElemNodeNormals(Real_t pfx[8],
                         Real_t pfy[8],
                         Real_t pfz[8],
                         const Real_t x[8],
                         const Real_t y[8],
                         const Real_t z[8])
{
   for (Index_t i = 0 ; i < 8 ; ++i) {
      pfx[i] = Real_t(0.0);
      pfy[i] = Real_t(0.0);
      pfz[i] = Real_t(0.0);
   }
   /* evaluate face one: nodes 0, 1, 2, 3 */
   SumElemFaceNormal(&pfx[0], &pfy[0], &pfz[0],
                  &pfx[1], &pfy[1], &pfz[1],
                  &pfx[2], &pfy[2], &pfz[2],
                  &pfx[3], &pfy[3], &pfz[3],
                  x[0], y[0], z[0], x[1], y[1], z[1],
                  x[2], y[2], z[2], x[3], y[3], z[3]);
   /* evaluate face two: nodes 0, 4, 5, 1 */
   SumElemFaceNormal(&pfx[0], &pfy[0], &pfz[0],
                  &pfx[4], &pfy[4], &pfz[4],
                  &pfx[5], &pfy[5], &pfz[5],
                  &pfx[1], &pfy[1], &pfz[1],
                  x[0], y[0], z[0], x[4], y[4], z[4],
                  x[5], y[5], z[5], x[1], y[1], z[1]);
   /* evaluate face three: nodes 1, 5, 6, 2 */
   SumElemFaceNormal(&pfx[1], &pfy[1], &pfz[1],
                  &pfx[5], &pfy[5], &pfz[5],
                  &pfx[6], &pfy[6], &pfz[6],
                  &pfx[2], &pfy[2], &pfz[2],
                  x[1], y[1], z[1], x[5], y[5], z[5],
                  x[6], y[6], z[6], x[2], y[2], z[2]);
   /* evaluate face four: nodes 2, 6, 7, 3 */
   SumElemFaceNormal(&pfx[2], &pfy[2], &pfz[2],
                  &pfx[6], &pfy[6], &pfz[6],
                  &pfx[7], &pfy[7], &pfz[7],
                  &pfx[3], &pfy[3], &pfz[3],
                  x[2], y[2], z[2], x[6], y[6], z[6],
                  x[7], y[7], z[7], x[3], y[3], z[3]);
   /* evaluate face five: nodes 3, 7, 4, 0 */
   SumElemFaceNormal(&pfx[3], &pfy[3], &pfz[3],
                  &pfx[7], &pfy[7], &pfz[7],
                  &pfx[4], &pfy[4], &pfz[4],
                  &pfx[0], &pfy[0], &pfz[0],
                  x[3], y[3], z[3], x[7], y[7], z[7],
                  x[4], y[4], z[4], x[0], y[0], z[0]);
   /* evaluate face six: nodes 4, 7, 6, 5 */
   SumElemFaceNormal(&pfx[4], &pfy[4], &pfz[4],
                  &pfx[7], &pfy[7], &pfz[7],
                  &pfx[6], &pfy[6], &pfz[6],
                  &pfx[5], &pfy[5], &pfz[5],
                  x[4], y[4], z[4], x[7], y[7], z[7],
                  x[6], y[6], z[6], x[5], y[5], z[5]);
}

/******************************************/

static inline
void SumElemStressesToNodeForces( const Real_t B[][8],
                                  const Real_t stress_xx,
                                  const Real_t stress_yy,
                                  const Real_t stress_zz,
                                  Real_t fx[], Real_t fy[], Real_t fz[] )
{
   for(Index_t i = 0; i < 8; i++) {
      fx[i] = -( stress_xx * B[0][i] );
      fy[i] = -( stress_yy * B[1][i]  );
      fz[i] = -( stress_zz * B[2][i] );
   }
}

/******************************************/

static inline
void IntegrateStressForElems( Domain &domain,
                              Real_t *sigxx, Real_t *sigyy, Real_t *sigzz,
                              Real_t *determ, Index_t numElem, Index_t numNode)
{
#if _OPENMP
   Index_t numthreads = omp_get_max_threads();
#else
   Index_t numthreads = 1;
#endif

   Index_t numElem8 = numElem * 8 ;
   Real_t *fx_elem;
   Real_t *fy_elem;
   Real_t *fz_elem;
   Real_t fx_local[8] ;
   Real_t fy_local[8] ;
   Real_t fz_local[8] ;


  if (numthreads > 1) {
     fx_elem = Allocate<Real_t>(numElem8) ;
     fy_elem = Allocate<Real_t>(numElem8) ;
     fz_elem = Allocate<Real_t>(numElem8) ;
  }
  // loop over all elements
  op_par_loop(IntegrateStressForElemsLoop, "IntegrateStressForElemsLoop", elems,
               op_arg_dat(p_x, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_x, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_x, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_y, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_y, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_y, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_z, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_z, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_z, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_fx, 0, p_nodelist, 1, "double", OP_WRITE), op_arg_dat(p_fx, 1, p_nodelist, 1, "double", OP_WRITE), op_arg_dat(p_fx, 2, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fx, 3, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fx, 4, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fx, 5, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fx, 6, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fx, 7, p_nodelist, 1, "double", OP_WRITE),
               op_arg_dat(p_fy, 0, p_nodelist, 1, "double", OP_WRITE), op_arg_dat(p_fy, 1, p_nodelist, 1, "double", OP_WRITE), op_arg_dat(p_fy, 2, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fy, 3, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fy, 4, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fy, 5, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fy, 6, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fy, 7, p_nodelist, 1, "double", OP_WRITE),
               op_arg_dat(p_fz, 0, p_nodelist, 1, "double", OP_WRITE), op_arg_dat(p_fz, 1, p_nodelist, 1, "double", OP_WRITE), op_arg_dat(p_fz, 2, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fz, 3, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fz, 4, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fz, 5, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fz, 6, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fz, 7, p_nodelist, 1, "double", OP_WRITE),
               op_arg_dat(p_determ, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_sigxx, -1, OP_ID, 3, "double", OP_READ)
               );
      

// #pragma omp parallel for firstprivate(numElem)
//   for( Index_t k=0 ; k<numElem ; ++k )
//   {
//    //  const Index_t* const elemToNode = domain.nodelist(k);
//     const Index_t* const elemToNode = &nodelist[k*Index_t(8)];
//     Real_t B[3][8] ;// shape function derivatives
//     Real_t x_local[8] ;
//     Real_t y_local[8] ;
//     Real_t z_local[8] ;

//     // get nodal coordinates from global arrays and copy into local arrays.
//     CollectDomainNodesToElemNodes(domain, elemToNode, x_local, y_local, z_local);

//     // Volume calculation involves extra work for numerical consistency
//     CalcElemShapeFunctionDerivatives(x_local, y_local, z_local,
//                                          B, &determ[k]);

//     CalcElemNodeNormals( B[0] , B[1], B[2],
//                           x_local, y_local, z_local );

//     if (numthreads > 1) {
//        // Eliminate thread writing conflicts at the nodes by giving
//        // each element its own copy to write to
//        SumElemStressesToNodeForces( B, sigxx[k], sigyy[k], sigzz[k],
//                                     &fx_elem[k*8],
//                                     &fy_elem[k*8],
//                                     &fz_elem[k*8] ) ;
//     }
//     else {
//        SumElemStressesToNodeForces( B, sigxx[k], sigyy[k], sigzz[k],
//                                     fx_local, fy_local, fz_local ) ;

//        // copy nodal force contributions to global force arrray.
//        for( Index_t lnode=0 ; lnode<8 ; ++lnode ) {
//           Index_t gnode = elemToNode[lnode];
//          //  domain.fx(gnode) += fx_local[lnode];
//          //  domain.fy(gnode) += fy_local[lnode];
//          //  domain.fz(gnode) += fz_local[lnode];
//           m_fx[gnode] += fx_local[lnode];
//           m_fy[gnode] += fy_local[lnode];
//           m_fz[gnode] += fz_local[lnode];
//        }
//       // std::cout<< "X:" << fx_local[0] << ",Y: " << fy_local[0] << ",Z: " << fz_local[0] << "\n";
//     }
//   }

  if (numthreads > 1) {
     // If threaded, then we need to copy the data out of the temporary
     // arrays used above into the final forces field
#pragma omp parallel for firstprivate(numNode)
     for( Index_t gnode=0 ; gnode<numNode ; ++gnode )
     {
        Index_t count = domain.nodeElemCount(gnode) ;
        Index_t *cornerList = domain.nodeElemCornerList(gnode) ;
        Real_t fx_tmp = Real_t(0.0) ;
        Real_t fy_tmp = Real_t(0.0) ;
        Real_t fz_tmp = Real_t(0.0) ;
        for (Index_t i=0 ; i < count ; ++i) {
           Index_t ielem = cornerList[i] ;
           fx_tmp += fx_elem[ielem] ;
           fy_tmp += fy_elem[ielem] ;
           fz_tmp += fz_elem[ielem] ;
        }
      //   domain.fx(gnode) = fx_tmp ;
      //   domain.fy(gnode) = fy_tmp ;
      //   domain.fz(gnode) = fz_tmp ;
        m_fx[gnode] = fx_tmp ;
        m_fy[gnode] = fy_tmp ;
        m_fz[gnode] = fz_tmp ;
     }
     Release(&fz_elem) ;
     Release(&fy_elem) ;
     Release(&fx_elem) ;
  }
}

/******************************************/

static inline
void VoluDer(const Real_t x0, const Real_t x1, const Real_t x2,
             const Real_t x3, const Real_t x4, const Real_t x5,
             const Real_t y0, const Real_t y1, const Real_t y2,
             const Real_t y3, const Real_t y4, const Real_t y5,
             const Real_t z0, const Real_t z1, const Real_t z2,
             const Real_t z3, const Real_t z4, const Real_t z5,
             Real_t* dvdx, Real_t* dvdy, Real_t* dvdz)
{
   const Real_t twelfth = Real_t(1.0) / Real_t(12.0) ;

   *dvdx =
      (y1 + y2) * (z0 + z1) - (y0 + y1) * (z1 + z2) +
      (y0 + y4) * (z3 + z4) - (y3 + y4) * (z0 + z4) -
      (y2 + y5) * (z3 + z5) + (y3 + y5) * (z2 + z5);
   *dvdy =
      - (x1 + x2) * (z0 + z1) + (x0 + x1) * (z1 + z2) -
      (x0 + x4) * (z3 + z4) + (x3 + x4) * (z0 + z4) +
      (x2 + x5) * (z3 + z5) - (x3 + x5) * (z2 + z5);

   *dvdz =
      - (y1 + y2) * (x0 + x1) + (y0 + y1) * (x1 + x2) -
      (y0 + y4) * (x3 + x4) + (y3 + y4) * (x0 + x4) +
      (y2 + y5) * (x3 + x5) - (y3 + y5) * (x2 + x5);

   *dvdx *= twelfth;
   *dvdy *= twelfth;
   *dvdz *= twelfth;
}

/******************************************/

static inline
void CalcElemVolumeDerivative(Real_t dvdx[8],
                              Real_t dvdy[8],
                              Real_t dvdz[8],
                              const Real_t x[8],
                              const Real_t y[8],
                              const Real_t z[8])
{
   VoluDer(x[1], x[2], x[3], x[4], x[5], x[7],
           y[1], y[2], y[3], y[4], y[5], y[7],
           z[1], z[2], z[3], z[4], z[5], z[7],
           &dvdx[0], &dvdy[0], &dvdz[0]);
   VoluDer(x[0], x[1], x[2], x[7], x[4], x[6],
           y[0], y[1], y[2], y[7], y[4], y[6],
           z[0], z[1], z[2], z[7], z[4], z[6],
           &dvdx[3], &dvdy[3], &dvdz[3]);
   VoluDer(x[3], x[0], x[1], x[6], x[7], x[5],
           y[3], y[0], y[1], y[6], y[7], y[5],
           z[3], z[0], z[1], z[6], z[7], z[5],
           &dvdx[2], &dvdy[2], &dvdz[2]);
   VoluDer(x[2], x[3], x[0], x[5], x[6], x[4],
           y[2], y[3], y[0], y[5], y[6], y[4],
           z[2], z[3], z[0], z[5], z[6], z[4],
           &dvdx[1], &dvdy[1], &dvdz[1]);
   VoluDer(x[7], x[6], x[5], x[0], x[3], x[1],
           y[7], y[6], y[5], y[0], y[3], y[1],
           z[7], z[6], z[5], z[0], z[3], z[1],
           &dvdx[4], &dvdy[4], &dvdz[4]);
   VoluDer(x[4], x[7], x[6], x[1], x[0], x[2],
           y[4], y[7], y[6], y[1], y[0], y[2],
           z[4], z[7], z[6], z[1], z[0], z[2],
           &dvdx[5], &dvdy[5], &dvdz[5]);
   VoluDer(x[5], x[4], x[7], x[2], x[1], x[3],
           y[5], y[4], y[7], y[2], y[1], y[3],
           z[5], z[4], z[7], z[2], z[1], z[3],
           &dvdx[6], &dvdy[6], &dvdz[6]);
   VoluDer(x[6], x[5], x[4], x[3], x[2], x[0],
           y[6], y[5], y[4], y[3], y[2], y[0],
           z[6], z[5], z[4], z[3], z[2], z[0],
           &dvdx[7], &dvdy[7], &dvdz[7]);
}

/******************************************/

static inline
void CalcElemFBHourglassForce(Real_t *xd, Real_t *yd, Real_t *zd,  Real_t hourgam[][4],
                              Real_t coefficient,
                              Real_t *hgfx, Real_t *hgfy, Real_t *hgfz )
{
   Real_t hxx[4];
   for(Index_t i = 0; i < 4; i++) {
      hxx[i] = hourgam[0][i] * xd[0] + hourgam[1][i] * xd[1] +
               hourgam[2][i] * xd[2] + hourgam[3][i] * xd[3] +
               hourgam[4][i] * xd[4] + hourgam[5][i] * xd[5] +
               hourgam[6][i] * xd[6] + hourgam[7][i] * xd[7];
   }
   for(Index_t i = 0; i < 8; i++) {
      hgfx[i] = coefficient *
                (hourgam[i][0] * hxx[0] + hourgam[i][1] * hxx[1] +
                 hourgam[i][2] * hxx[2] + hourgam[i][3] * hxx[3]);
   }
   for(Index_t i = 0; i < 4; i++) {
      hxx[i] = hourgam[0][i] * yd[0] + hourgam[1][i] * yd[1] +
               hourgam[2][i] * yd[2] + hourgam[3][i] * yd[3] +
               hourgam[4][i] * yd[4] + hourgam[5][i] * yd[5] +
               hourgam[6][i] * yd[6] + hourgam[7][i] * yd[7];
   }
   for(Index_t i = 0; i < 8; i++) {
      hgfy[i] = coefficient *
                (hourgam[i][0] * hxx[0] + hourgam[i][1] * hxx[1] +
                 hourgam[i][2] * hxx[2] + hourgam[i][3] * hxx[3]);
   }
   for(Index_t i = 0; i < 4; i++) {
      hxx[i] = hourgam[0][i] * zd[0] + hourgam[1][i] * zd[1] +
               hourgam[2][i] * zd[2] + hourgam[3][i] * zd[3] +
               hourgam[4][i] * zd[4] + hourgam[5][i] * zd[5] +
               hourgam[6][i] * zd[6] + hourgam[7][i] * zd[7];
   }
   for(Index_t i = 0; i < 8; i++) {
      hgfz[i] = coefficient *
                (hourgam[i][0] * hxx[0] + hourgam[i][1] * hxx[1] +
                 hourgam[i][2] * hxx[2] + hourgam[i][3] * hxx[3]);
   }
}

/******************************************/

static inline
void CalcFBHourglassForceForElems( Domain &domain,
                                   Real_t *determ,
                                   Real_t *x8n, Real_t *y8n, Real_t *z8n,
                                   Real_t *dvdx, Real_t *dvdy, Real_t *dvdz,
                                   Real_t hourg, Index_t numElem,
                                   Index_t numNode)
{

#if _OPENMP
   Index_t numthreads = omp_get_max_threads();
#else
   Index_t numthreads = 1;
#endif
   /*************************************************
    *
    *     FUNCTION: Calculates the Flanagan-Belytschko anti-hourglass
    *               force.
    *
    *************************************************/
  
   Index_t numElem8 = numElem * 8 ;

   Real_t *fx_elem; 
   Real_t *fy_elem; 
   Real_t *fz_elem; 

   if(numthreads > 1) {
      fx_elem = Allocate<Real_t>(numElem8) ;
      fy_elem = Allocate<Real_t>(numElem8) ;
      fz_elem = Allocate<Real_t>(numElem8) ;
   }

   // Real_t  gamma[4][8];

   // gamma[0][0] = Real_t( 1.);
   // gamma[0][1] = Real_t( 1.);
   // gamma[0][2] = Real_t(-1.);
   // gamma[0][3] = Real_t(-1.);
   // gamma[0][4] = Real_t(-1.);
   // gamma[0][5] = Real_t(-1.);
   // gamma[0][6] = Real_t( 1.);
   // gamma[0][7] = Real_t( 1.);
   // gamma[1][0] = Real_t( 1.);
   // gamma[1][1] = Real_t(-1.);
   // gamma[1][2] = Real_t(-1.);
   // gamma[1][3] = Real_t( 1.);
   // gamma[1][4] = Real_t(-1.);
   // gamma[1][5] = Real_t( 1.);
   // gamma[1][6] = Real_t( 1.);
   // gamma[1][7] = Real_t(-1.);
   // gamma[2][0] = Real_t( 1.);
   // gamma[2][1] = Real_t(-1.);
   // gamma[2][2] = Real_t( 1.);
   // gamma[2][3] = Real_t(-1.);
   // gamma[2][4] = Real_t( 1.);
   // gamma[2][5] = Real_t(-1.);
   // gamma[2][6] = Real_t( 1.);
   // gamma[2][7] = Real_t(-1.);
   // gamma[3][0] = Real_t(-1.);
   // gamma[3][1] = Real_t( 1.);
   // gamma[3][2] = Real_t(-1.);
   // gamma[3][3] = Real_t( 1.);
   // gamma[3][4] = Real_t( 1.);
   // gamma[3][5] = Real_t(-1.);
   // gamma[3][6] = Real_t( 1.);
   // gamma[3][7] = Real_t(-1.);

/*************************************************/
/*    compute the hourglass modes */

   op_par_loop(FBHourglassForceForElems, "CalcFBHourglassForceForElems", elems,
               op_arg_dat(p_xd, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_xd, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_xd, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_yd, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_yd, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_yd, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_zd, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_zd, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_zd, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_fx, 0, p_nodelist, 1, "double", OP_WRITE), op_arg_dat(p_fx, 1, p_nodelist, 1, "double", OP_WRITE), op_arg_dat(p_fx, 2, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fx, 3, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fx, 4, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fx, 5, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fx, 6, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fx, 7, p_nodelist, 1, "double", OP_WRITE),
               op_arg_dat(p_fy, 0, p_nodelist, 1, "double", OP_WRITE), op_arg_dat(p_fy, 1, p_nodelist, 1, "double", OP_WRITE), op_arg_dat(p_fy, 2, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fy, 3, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fy, 4, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fy, 5, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fy, 6, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fy, 7, p_nodelist, 1, "double", OP_WRITE),
               op_arg_dat(p_fz, 0, p_nodelist, 1, "double", OP_WRITE), op_arg_dat(p_fz, 1, p_nodelist, 1, "double", OP_WRITE), op_arg_dat(p_fz, 2, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fz, 3, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fz, 4, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fz, 5, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fz, 6, p_nodelist, 1, "double", OP_WRITE),op_arg_dat(p_fz, 7, p_nodelist, 1, "double", OP_WRITE),
               op_arg_dat(p_dvdx, -1, OP_ID, 8, "double", OP_READ), 
               op_arg_dat(p_dvdy, -1, OP_ID, 8, "double", OP_READ), 
               op_arg_dat(p_dvdz, -1, OP_ID, 8, "double", OP_READ), 
               op_arg_dat(p_x8n, -1, OP_ID, 8, "double", OP_READ), 
               op_arg_dat(p_y8n, -1, OP_ID, 8, "double", OP_READ), 
               op_arg_dat(p_z8n, -1, OP_ID, 8, "double", OP_READ),
               op_arg_dat(p_determ, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_ss, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_elemMass, -1, OP_ID, 1, "double", OP_READ)
   );
// #pragma omp parallel for firstprivate(numElem, hourg)
//    for(Index_t i2=0;i2<numElem;++i2){
//       Real_t *fx_local, *fy_local, *fz_local ;
//       Real_t hgfx[8], hgfy[8], hgfz[8] ;
//       Real_t coefficient;

//       Real_t hourgam[8][4];
//       Real_t xd1[8], yd1[8], zd1[8] ;

//       // const Index_t *elemToNode = domain.nodelist(i2);
//       const Index_t *elemToNode = &nodelist[i2*Index_t(8)];
//       Index_t i3=8*i2;
//       Real_t volinv=Real_t(1.0)/determ[i2];
//       Real_t ss1, mass1, volume13 ;
//       for(Index_t i1=0;i1<4;++i1){

//          Real_t hourmodx =
//             x8n[i3] * gamma[i1][0] + x8n[i3+1] * gamma[i1][1] +
//             x8n[i3+2] * gamma[i1][2] + x8n[i3+3] * gamma[i1][3] +
//             x8n[i3+4] * gamma[i1][4] + x8n[i3+5] * gamma[i1][5] +
//             x8n[i3+6] * gamma[i1][6] + x8n[i3+7] * gamma[i1][7];

//          Real_t hourmody =
//             y8n[i3] * gamma[i1][0] + y8n[i3+1] * gamma[i1][1] +
//             y8n[i3+2] * gamma[i1][2] + y8n[i3+3] * gamma[i1][3] +
//             y8n[i3+4] * gamma[i1][4] + y8n[i3+5] * gamma[i1][5] +
//             y8n[i3+6] * gamma[i1][6] + y8n[i3+7] * gamma[i1][7];

//          Real_t hourmodz =
//             z8n[i3] * gamma[i1][0] + z8n[i3+1] * gamma[i1][1] +
//             z8n[i3+2] * gamma[i1][2] + z8n[i3+3] * gamma[i1][3] +
//             z8n[i3+4] * gamma[i1][4] + z8n[i3+5] * gamma[i1][5] +
//             z8n[i3+6] * gamma[i1][6] + z8n[i3+7] * gamma[i1][7];

//          hourgam[0][i1] = gamma[i1][0] -  volinv*(dvdx[i3  ] * hourmodx +
//                                                   dvdy[i3  ] * hourmody +
//                                                   dvdz[i3  ] * hourmodz );

//          hourgam[1][i1] = gamma[i1][1] -  volinv*(dvdx[i3+1] * hourmodx +
//                                                   dvdy[i3+1] * hourmody +
//                                                   dvdz[i3+1] * hourmodz );

//          hourgam[2][i1] = gamma[i1][2] -  volinv*(dvdx[i3+2] * hourmodx +
//                                                   dvdy[i3+2] * hourmody +
//                                                   dvdz[i3+2] * hourmodz );

//          hourgam[3][i1] = gamma[i1][3] -  volinv*(dvdx[i3+3] * hourmodx +
//                                                   dvdy[i3+3] * hourmody +
//                                                   dvdz[i3+3] * hourmodz );

//          hourgam[4][i1] = gamma[i1][4] -  volinv*(dvdx[i3+4] * hourmodx +
//                                                   dvdy[i3+4] * hourmody +
//                                                   dvdz[i3+4] * hourmodz );

//          hourgam[5][i1] = gamma[i1][5] -  volinv*(dvdx[i3+5] * hourmodx +
//                                                   dvdy[i3+5] * hourmody +
//                                                   dvdz[i3+5] * hourmodz );

//          hourgam[6][i1] = gamma[i1][6] -  volinv*(dvdx[i3+6] * hourmodx +
//                                                   dvdy[i3+6] * hourmody +
//                                                   dvdz[i3+6] * hourmodz );

//          hourgam[7][i1] = gamma[i1][7] -  volinv*(dvdx[i3+7] * hourmodx +
//                                                   dvdy[i3+7] * hourmody +
//                                                   dvdz[i3+7] * hourmodz );

//       }

//       /* compute forces */
//       /* store forces into h arrays (force arrays) */

//       // ss1=domain.ss(i2);
//       ss1=ss[i2];

//       // mass1=domain.elemMass(i2);
//       mass1=elemMass[i2];
//       volume13=CBRT(determ[i2]);

//       Index_t n0si2 = elemToNode[0];
//       Index_t n1si2 = elemToNode[1];
//       Index_t n2si2 = elemToNode[2];
//       Index_t n3si2 = elemToNode[3];
//       Index_t n4si2 = elemToNode[4];
//       Index_t n5si2 = elemToNode[5];
//       Index_t n6si2 = elemToNode[6];
//       Index_t n7si2 = elemToNode[7];

//       xd1[0] = xd[n0si2];
//       xd1[1] = xd[n1si2];
//       xd1[2] = xd[n2si2];
//       xd1[3] = xd[n3si2];
//       xd1[4] = xd[n4si2];
//       xd1[5] = xd[n5si2];
//       xd1[6] = xd[n6si2];
//       xd1[7] = xd[n7si2];

//       yd1[0] = yd[n0si2];
//       yd1[1] = yd[n1si2];
//       yd1[2] = yd[n2si2];
//       yd1[3] = yd[n3si2];
//       yd1[4] = yd[n4si2];
//       yd1[5] = yd[n5si2];
//       yd1[6] = yd[n6si2];
//       yd1[7] = yd[n7si2];

//       zd1[0] = zd[n0si2];
//       zd1[1] = zd[n1si2];
//       zd1[2] = zd[n2si2];
//       zd1[3] = zd[n3si2];
//       zd1[4] = zd[n4si2];
//       zd1[5] = zd[n5si2];
//       zd1[6] = zd[n6si2];
//       zd1[7] = zd[n7si2];

//       coefficient = - hourg * Real_t(0.01) * ss1 * mass1 / volume13;

//       CalcElemFBHourglassForce(xd1,yd1,zd1,
//                       hourgam,
//                       coefficient, hgfx, hgfy, hgfz);

//       // With the threaded version, we write into local arrays per elem
//       // so we don't have to worry about race conditions
//       if (numthreads > 1) {
//          fx_local = &fx_elem[i3] ;
//          fx_local[0] = hgfx[0];
//          fx_local[1] = hgfx[1];
//          fx_local[2] = hgfx[2];
//          fx_local[3] = hgfx[3];
//          fx_local[4] = hgfx[4];
//          fx_local[5] = hgfx[5];
//          fx_local[6] = hgfx[6];
//          fx_local[7] = hgfx[7];

//          fy_local = &fy_elem[i3] ;
//          fy_local[0] = hgfy[0];
//          fy_local[1] = hgfy[1];
//          fy_local[2] = hgfy[2];
//          fy_local[3] = hgfy[3];
//          fy_local[4] = hgfy[4];
//          fy_local[5] = hgfy[5];
//          fy_local[6] = hgfy[6];
//          fy_local[7] = hgfy[7];

//          fz_local = &fz_elem[i3] ;
//          fz_local[0] = hgfz[0];
//          fz_local[1] = hgfz[1];
//          fz_local[2] = hgfz[2];
//          fz_local[3] = hgfz[3];
//          fz_local[4] = hgfz[4];
//          fz_local[5] = hgfz[5];
//          fz_local[6] = hgfz[6];
//          fz_local[7] = hgfz[7];
//       }
//       else {
//          // domain.fx(n0si2) += hgfx[0];
//          // domain.fy(n0si2) += hgfy[0];
//          // domain.fz(n0si2) += hgfz[0];

//          // domain.fx(n1si2) += hgfx[1];
//          // domain.fy(n1si2) += hgfy[1];
//          // domain.fz(n1si2) += hgfz[1];

//          // domain.fx(n2si2) += hgfx[2];
//          // domain.fy(n2si2) += hgfy[2];
//          // domain.fz(n2si2) += hgfz[2];

//          // domain.fx(n3si2) += hgfx[3];
//          // domain.fy(n3si2) += hgfy[3];
//          // domain.fz(n3si2) += hgfz[3];

//          // domain.fx(n4si2) += hgfx[4];
//          // domain.fy(n4si2) += hgfy[4];
//          // domain.fz(n4si2) += hgfz[4];

//          // domain.fx(n5si2) += hgfx[5];
//          // domain.fy(n5si2) += hgfy[5];
//          // domain.fz(n5si2) += hgfz[5];

//          // domain.fx(n6si2) += hgfx[6];
//          // domain.fy(n6si2) += hgfy[6];
//          // domain.fz(n6si2) += hgfz[6];

//          // domain.fx(n7si2) += hgfx[7];
//          // domain.fy(n7si2) += hgfy[7];
//          // domain.fz(n7si2) += hgfz[7];
//          m_fx[n0si2] += hgfx[0];
//          m_fy[n0si2] += hgfy[0];
//          m_fz[n0si2] += hgfz[0];

//          m_fx[n1si2] += hgfx[1];
//          m_fy[n1si2] += hgfy[1];
//          m_fz[n1si2] += hgfz[1];

//          m_fx[n2si2] += hgfx[2];
//          m_fy[n2si2] += hgfy[2];
//          m_fz[n2si2] += hgfz[2];

//          m_fx[n3si2] += hgfx[3];
//          m_fy[n3si2] += hgfy[3];
//          m_fz[n3si2] += hgfz[3];

//          m_fx[n4si2] += hgfx[4];
//          m_fy[n4si2] += hgfy[4];
//          m_fz[n4si2] += hgfz[4];

//          m_fx[n5si2] += hgfx[5];
//          m_fy[n5si2] += hgfy[5];
//          m_fz[n5si2] += hgfz[5];

//          m_fx[n6si2] += hgfx[6];
//          m_fy[n6si2] += hgfy[6];
//          m_fz[n6si2] += hgfz[6];

//          m_fx[n7si2] += hgfx[7];
//          m_fy[n7si2] += hgfy[7];
//          m_fz[n7si2] += hgfz[7];
//       }
//    }

//    if (numthreads > 1) {
//      // Collect the data from the local arrays into the final force arrays
// #pragma omp parallel for firstprivate(numNode)
//       for( Index_t gnode=0 ; gnode<numNode ; ++gnode )
//       {
//          Index_t count = domain.nodeElemCount(gnode) ;
//          Index_t *cornerList = domain.nodeElemCornerList(gnode) ;
//          Real_t fx_tmp = Real_t(0.0) ;
//          Real_t fy_tmp = Real_t(0.0) ;
//          Real_t fz_tmp = Real_t(0.0) ;
//          for (Index_t i=0 ; i < count ; ++i) {
//             Index_t ielem = cornerList[i] ;
//             fx_tmp += fx_elem[ielem] ;
//             fy_tmp += fy_elem[ielem] ;
//             fz_tmp += fz_elem[ielem] ;
//          }
//          // domain.fx(gnode) += fx_tmp ;
//          // domain.fy(gnode) += fy_tmp ;
//          // domain.fz(gnode) += fz_tmp ;
//          m_fx[gnode] += fx_tmp ;
//          m_fy[gnode] += fy_tmp ;
//          m_fz[gnode] += fz_tmp ;
//       }
//       Release(&fz_elem) ;
//       Release(&fy_elem) ;
//       Release(&fx_elem) ;
//    }
}

/******************************************/

static inline
void CalcHourglassControlForElems(Domain& domain,
                                  Real_t determ[], Real_t hgcoef)
{
   // Index_t numElem = domain.numElem() ;
   Index_t numElem = m_numElem ;
   Index_t numElem8 = numElem * 8 ;
   // Real_t *dvdx = Allocate<Real_t>(numElem8) ;
   // Real_t *dvdy = Allocate<Real_t>(numElem8) ;
   // Real_t *dvdz = Allocate<Real_t>(numElem8) ;
   // Real_t *x8n  = Allocate<Real_t>(numElem8) ;
   // Real_t *y8n  = Allocate<Real_t>(numElem8) ;
   // Real_t *z8n  = Allocate<Real_t>(numElem8) ;

   // op_dat p_dvdx = op_decl_dat_temp(elems, 8, "double",dvdx, "dvdx");
   // op_dat p_dvdy = op_decl_dat_temp(elems, 8, "double",dvdy, "dvdy");
   // op_dat p_dvdz = op_decl_dat_temp(elems, 8, "double",dvdz, "dvdz");
   // op_dat p_x8n = op_decl_dat_temp(elems, 8, "double",x8n, "x8n");
   // op_dat p_y8n = op_decl_dat_temp(elems, 8, "double",y8n, "y8n");
   // op_dat p_z8n = op_decl_dat_temp(elems, 8, "double",z8n, "z8n");

   op_par_loop(CalcVolumeDerivatives, "CalcVolumeDerivatives", elems,
               op_arg_dat(p_x, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_x, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_x, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_y, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_y, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_y, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_z, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_z, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_z, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_dvdx, -1, OP_ID, 8, "double", OP_RW), 
               op_arg_dat(p_dvdy, -1, OP_ID, 8, "double", OP_RW), 
               op_arg_dat(p_dvdz, -1, OP_ID, 8, "double", OP_RW), 
               op_arg_dat(p_x8n, -1, OP_ID, 8, "double", OP_WRITE), 
               op_arg_dat(p_y8n, -1, OP_ID, 8, "double", OP_WRITE), 
               op_arg_dat(p_z8n, -1, OP_ID, 8, "double", OP_WRITE), 
               op_arg_dat(p_v, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_determ, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_volo, -1, OP_ID, 1, "double", OP_READ)
               );
//    /* start loop over elements */
// #pragma omp parallel for firstprivate(numElem)
//    for (Index_t i=0 ; i<numElem ; ++i){
//       Real_t  x1[8],  y1[8],  z1[8] ;
//       Real_t pfx[8], pfy[8], pfz[8] ;

//       // Index_t* elemToNode = domain.nodelist(i);
//       Index_t* elemToNode = &nodelist[i*Index_t(8)];
//       CollectDomainNodesToElemNodes(domain, elemToNode, x1, y1, z1);

//       CalcElemVolumeDerivative(pfx, pfy, pfz, x1, y1, z1);

//       /* load into temporary storage for FB Hour Glass control */
//       for(Index_t ii=0;ii<8;++ii){
//          Index_t jj=8*i+ii;

//          dvdx[jj] = pfx[ii];
//          dvdy[jj] = pfy[ii];
//          dvdz[jj] = pfz[ii];
//          x8n[jj]  = x1[ii];
//          y8n[jj]  = y1[ii];
//          z8n[jj]  = z1[ii];
//       }

//       // determ[i] = domain.volo(i) * domain.v(i);
//       determ[i] = volo[i] * v[i];

//       /* Do a check for negative volumes */
//       // if ( domain.v(i) <= Real_t(0.0) ) {
//       if ( v[i] <= Real_t(0.0) ) {
// #if USE_MPI         
//          MPI_Abort(MPI_COMM_WORLD, VolumeError) ;
// #else
//          exit(VolumeError);
// #endif
//       }
//    }

   if ( hgcoef > Real_t(0.) ) {
      CalcFBHourglassForceForElems( domain,
                                    determ, x8n, y8n, z8n, dvdx, dvdy, dvdz,
                                    // hgcoef, numElem, domain.numNode()) ;
                                    hgcoef, numElem, m_numNode) ;
   }

   // op_free_dat_temp(p_dvdx);
   // op_free_dat_temp(p_dvdy);
   // op_free_dat_temp(p_dvdz);
   // op_free_dat_temp(p_x8n);
   // op_free_dat_temp(p_y8n);
   // op_free_dat_temp(p_z8n);
   // Release(&z8n) ;
   // Release(&y8n) ;
   // Release(&x8n) ;
   // Release(&dvdz) ;
   // Release(&dvdy) ;
   // Release(&dvdx) ;

   return ;
}

/******************************************/

static inline
void CalcVolumeForceForElems(Domain& domain)
{
   // Index_t numElem = domain.numElem() ;
   Index_t numElem = m_numElem ;
   if (numElem != 0) {
      // Real_t  hgcoef = domain.hgcoef() ;
      Real_t  hgcoef = m_hgcoef ;
      // Real_t *sigxx  = Allocate<Real_t>(numElem) ;
      // Real_t *sigyy  = Allocate<Real_t>(numElem) ;
      // Real_t *sigzz  = Allocate<Real_t>(numElem) ;
      // Real_t *determ = Allocate<Real_t>(numElem) ;

      // Real_t *sigxx  = Allocate<Real_t>(3*numElem) ;
      //op_decl_dat_temp(elems, 3, "double", sigxx, "sigxx");
      //op_decl_dat_temp(elems, 1, "double", determ, "determ");

      /* Sum contributions to total stress tensor */
      InitStressTermsForElems(domain, sigxx, sigyy, sigzz, numElem);
      #if DEBUGGING
      std::cout << "Done initting stress terms\n";
      #endif

      // call elemlib stress integration loop to produce nodal forces from
      // material stresses.
      IntegrateStressForElems( domain,
                               sigxx, sigyy, sigzz, determ, numElem,
                              //  domain.numNode()) ;
                               m_numNode) ;

      // check for negative element volume
      op_par_loop(CheckForNegativeElementVolume, "CheckForNegativeElementVolume", elems,
                  op_arg_dat(p_determ, -1, OP_ID, 1, "double", OP_READ));
// #pragma omp parallel for firstprivate(numElem)
//       for ( Index_t k=0 ; k<numElem ; ++k ) {      
//          if (determ[k] <= Real_t(0.0)) {
// #if USE_MPI            
//             MPI_Abort(MPI_COMM_WORLD, VolumeError) ;
// #else
//          std::cout<<"After volume check\n";
//             exit(VolumeError);
// #endif
//          }
//       }

      CalcHourglassControlForElems(domain, determ, hgcoef) ;

      //op_free_dat_temp(determ);
      //op_free_dat_temp(sigxx);
      // Release(&determ) ;
      // Release(&sigzz) ;
      // Release(&sigyy) ;
      // Release(&sigxx) ;
   }
}

/******************************************/

static inline void CalcForceForNodes(Domain& domain)
{
//   Index_t numNode = domain.numNode() ;
  Index_t numNode = m_numNode ;

#if USE_MPI  
  CommRecv(domain, MSG_COMM_SBN, 3,
           domain.sizeX() + 1, domain.sizeY() + 1, domain.sizeZ() + 1,
           true, false) ;
#endif  

   op_par_loop(setForceToZero, "setForceToZero", nodes,
               op_arg_dat(p_fx, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_fy, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_fz, -1, OP_ID, 1, "double", OP_WRITE));
   #if DEBUGGING
   std::cout << "Done Setting force to zero\n";
   #endif
// #pragma omp parallel for firstprivate(numNode)
//   for (Index_t i=0; i<numNode; ++i) {
//    //   domain.fx(i) = Real_t(0.0) ;
//    //   domain.fy(i) = Real_t(0.0) ;
//    //   domain.fz(i) = Real_t(0.0) ;
//      m_fx[i] = Real_t(0.0) ;
//      m_fy[i] = Real_t(0.0) ;
//      m_fz[i] = Real_t(0.0) ;
//   }

  /* Calcforce calls partial, force, hourq */
  CalcVolumeForceForElems(domain) ;

#if USE_MPI  
  Domain_member fieldData[3] ;
  fieldData[0] = &Domain::fx ;
  fieldData[1] = &Domain::fy ;
  fieldData[2] = &Domain::fz ;
  
  CommSend(domain, MSG_COMM_SBN, 3, fieldData,
           domain.sizeX() + 1, domain.sizeY() + 1, domain.sizeZ() +  1,
           true, false) ;
  CommSBN(domain, 3, fieldData) ;
#endif  
}

/******************************************/

static inline
void CalcAccelerationForNodes(Domain &domain, Index_t numNode)
{   
   op_par_loop(CalcAccelForNodes, "CalcAccelForNodes", nodes,
               op_arg_dat(p_xdd, -1, OP_ID, 1, "double", OP_WRITE), op_arg_dat(p_ydd, -1, OP_ID, 1, "double", OP_WRITE), op_arg_dat(p_zdd, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_fx, -1, OP_ID, 1, "double", OP_READ), op_arg_dat(p_fy, -1, OP_ID, 1, "double", OP_READ), op_arg_dat(p_fz, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_nodalMass, -1, OP_ID, 1, "double", OP_READ)
               );
// #pragma omp parallel for firstprivate(numNode)
//    for (Index_t i = 0; i < numNode; ++i) {
//       // domain.xdd(i) = domain.fx(i) / domain.nodalMass(i);
//       // domain.ydd(i) = domain.fy(i) / domain.nodalMass(i);
//       // domain.zdd(i) = domain.fz(i) / domain.nodalMass(i);
//       xdd[i] = m_fx[i] / nodalMass[i];
//       ydd[i] = m_fy[i] / nodalMass[i];
//       zdd[i] = m_fz[i] / nodalMass[i];
   // }
}

/******************************************/

static inline
void ApplyAccelerationBoundaryConditionsForNodes(Domain& domain)
{
   // Index_t size = domain.sizeX();
   Index_t size = sizeX;
   Index_t numNodeBC = (size+1)*(size+1) ;

#pragma omp parallel
   {
      if (!domain.symmXempty() != 0) {
#pragma omp for nowait firstprivate(numNodeBC)
         for(Index_t i=0 ; i<numNodeBC ; ++i)
            // xdd[domain.symmX(i)] = Real_t(0.0) ;
            xdd[symmX[i]] = Real_t(0.0) ;
      }

      if (!domain.symmYempty() != 0) {
#pragma omp for nowait firstprivate(numNodeBC)
         for(Index_t i=0 ; i<numNodeBC ; ++i)
            // ydd[domain.symmY(i)] = Real_t(0.0) ;
            ydd[symmY[i]] = Real_t(0.0) ;
      }

      if (!domain.symmZempty() != 0) {
#pragma omp for nowait firstprivate(numNodeBC)
         for(Index_t i=0 ; i<numNodeBC ; ++i)
            // zdd[domain.symmZ(i)] = Real_t(0.0) ;
            zdd[symmZ[i]] = Real_t(0.0) ;
      }
   }
}

/******************************************/

static inline
void CalcVelocityForNodes(Domain &domain,const Real_t dt, const Real_t u_cut,
                          Index_t numNode)
{
   op_par_loop(CalcVeloForNodes, "CalcVeloForNodes", nodes,
               op_arg_dat(p_xd, -1, OP_ID, 1, "double", OP_WRITE), op_arg_dat(p_yd, -1, OP_ID, 1, "double", OP_WRITE), op_arg_dat(p_zd, -1, OP_ID, 1, "double", OP_WRITE), 
               op_arg_dat(p_xdd, -1, OP_ID, 1, "double", OP_READ), op_arg_dat(p_ydd, -1, OP_ID, 1, "double", OP_READ), op_arg_dat(p_zdd, -1, OP_ID, 1, "double", OP_READ),
               op_arg_gbl(&dt, 1, "double", OP_READ) 
               );
// #pragma omp parallel for firstprivate(numNode)
//    for ( Index_t i = 0 ; i < numNode ; ++i )
//    {
//      Real_t xdtmp, ydtmp, zdtmp ;

//    //   xdtmp = domain.xd(i) + domain.xdd(i) * dt ;
//      xdtmp = xd[i] + xdd[i] * dt ;
//      if( FABS(xdtmp) < u_cut ) xdtmp = Real_t(0.0);
//    //   domain.xd(i) = xdtmp ;
//      xd[i] = xdtmp ;

//    //   ydtmp = domain.yd(i) + domain.ydd(i) * dt ;
//      ydtmp = yd[i] + ydd[i] * dt ;
//      if( FABS(ydtmp) < u_cut ) ydtmp = Real_t(0.0);
//    //   domain.yd(i) = ydtmp ;
//      yd[i] = ydtmp ;

//    //   zdtmp = domain.zd(i) + domain.zdd(i) * dt ;
//      zdtmp = zd[i] + zdd[i] * dt ;
//      if( FABS(zdtmp) < u_cut ) zdtmp = Real_t(0.0);
//    //   domain.zd(i) = zdtmp ;
//      zd[i] = zdtmp ;
//    }
}

/******************************************/

static inline
void CalcPositionForNodes(Domain &domain, const Real_t dt, Index_t numNode)
{

   op_par_loop(CalcPosForNodes, "CalcPosForNodes", nodes,
               op_arg_dat(p_x, -1, OP_ID, 1, "double", OP_WRITE), op_arg_dat(p_y, -1, OP_ID, 1, "double", OP_WRITE), op_arg_dat(p_z, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_xd, -1, OP_ID, 1, "double", OP_READ), op_arg_dat(p_yd, -1, OP_ID, 1, "double", OP_READ), op_arg_dat(p_zd, -1, OP_ID, 1, "double", OP_READ),
               op_arg_gbl(&dt, 1, "double", OP_READ)
               );
// #pragma omp parallel for firstprivate(numNode)
//    for ( Index_t i = 0 ; i < numNode ; ++i )
//    {
//    //   domain.x(i) += domain.xd(i) * dt ;
//      x[i] += xd[i] * dt ;
//    //   domain.y(i) += domain.yd(i) * dt ;
//      y[i] += yd[i] * dt ;
//    //   domain.z(i) += domain.zd(i) * dt ;
//      z[i] += zd[i] * dt ;
//    }
}

/******************************************/

static inline
void LagrangeNodal(Domain& domain)
{
#ifdef SEDOV_SYNC_POS_VEL_EARLY
   Domain_member fieldData[6] ;
#endif

   // const Real_t delt = domain.deltatime() ;
   const Real_t delt = m_deltatime ;
   // Real_t u_cut = domain.u_cut() ;
   Real_t u_cut = m_u_cut ;

  /* time of boundary condition evaluation is beginning of step for force and
   * acceleration boundary conditions. */
  CalcForceForNodes(domain);

#if USE_MPI  
#ifdef SEDOV_SYNC_POS_VEL_EARLY
   CommRecv(domain, MSG_SYNC_POS_VEL, 6,
            domain.sizeX() + 1, domain.sizeY() + 1, domain.sizeZ() + 1,
            false, false) ;
#endif
#endif
   
   // CalcAccelerationForNodes(domain, domain.numNode());
   CalcAccelerationForNodes(domain, m_numNode);
   
   ApplyAccelerationBoundaryConditionsForNodes(domain);

   // CalcVelocityForNodes( domain, delt, u_cut, domain.numNode()) ;
   CalcVelocityForNodes( domain, delt, u_cut, m_numNode) ;

   CalcPositionForNodes( domain, delt, m_numNode );
#if USE_MPI
#ifdef SEDOV_SYNC_POS_VEL_EARLY
  fieldData[0] = &Domain::x ;
  fieldData[1] = &Domain::y ;
  fieldData[2] = &Domain::z ;
  fieldData[3] = &Domain::xd ;
  fieldData[4] = &Domain::yd ;
  fieldData[5] = &Domain::zd ;

   CommSend(domain, MSG_SYNC_POS_VEL, 6, fieldData,
            domain.sizeX() + 1, domain.sizeY() + 1, domain.sizeZ() + 1,
            false, false) ;
   CommSyncPosVel(domain) ;
#endif
#endif
   
  return;
}

/******************************************/

static inline
Real_t CalcElemVolume( const Real_t x0, const Real_t x1,
               const Real_t x2, const Real_t x3,
               const Real_t x4, const Real_t x5,
               const Real_t x6, const Real_t x7,
               const Real_t y0, const Real_t y1,
               const Real_t y2, const Real_t y3,
               const Real_t y4, const Real_t y5,
               const Real_t y6, const Real_t y7,
               const Real_t z0, const Real_t z1,
               const Real_t z2, const Real_t z3,
               const Real_t z4, const Real_t z5,
               const Real_t z6, const Real_t z7 )
{
  Real_t twelveth = Real_t(1.0)/Real_t(12.0);

  Real_t dx61 = x6 - x1;
  Real_t dy61 = y6 - y1;
  Real_t dz61 = z6 - z1;

  Real_t dx70 = x7 - x0;
  Real_t dy70 = y7 - y0;
  Real_t dz70 = z7 - z0;

  Real_t dx63 = x6 - x3;
  Real_t dy63 = y6 - y3;
  Real_t dz63 = z6 - z3;

  Real_t dx20 = x2 - x0;
  Real_t dy20 = y2 - y0;
  Real_t dz20 = z2 - z0;

  Real_t dx50 = x5 - x0;
  Real_t dy50 = y5 - y0;
  Real_t dz50 = z5 - z0;

  Real_t dx64 = x6 - x4;
  Real_t dy64 = y6 - y4;
  Real_t dz64 = z6 - z4;

  Real_t dx31 = x3 - x1;
  Real_t dy31 = y3 - y1;
  Real_t dz31 = z3 - z1;

  Real_t dx72 = x7 - x2;
  Real_t dy72 = y7 - y2;
  Real_t dz72 = z7 - z2;

  Real_t dx43 = x4 - x3;
  Real_t dy43 = y4 - y3;
  Real_t dz43 = z4 - z3;

  Real_t dx57 = x5 - x7;
  Real_t dy57 = y5 - y7;
  Real_t dz57 = z5 - z7;

  Real_t dx14 = x1 - x4;
  Real_t dy14 = y1 - y4;
  Real_t dz14 = z1 - z4;

  Real_t dx25 = x2 - x5;
  Real_t dy25 = y2 - y5;
  Real_t dz25 = z2 - z5;

#define TRIPLE_PRODUCT(x1, y1, z1, x2, y2, z2, x3, y3, z3) \
   ((x1)*((y2)*(z3) - (z2)*(y3)) + (x2)*((z1)*(y3) - (y1)*(z3)) + (x3)*((y1)*(z2) - (z1)*(y2)))

  Real_t volume =
    TRIPLE_PRODUCT(dx31 + dx72, dx63, dx20,
       dy31 + dy72, dy63, dy20,
       dz31 + dz72, dz63, dz20) +
    TRIPLE_PRODUCT(dx43 + dx57, dx64, dx70,
       dy43 + dy57, dy64, dy70,
       dz43 + dz57, dz64, dz70) +
    TRIPLE_PRODUCT(dx14 + dx25, dx61, dx50,
       dy14 + dy25, dy61, dy50,
       dz14 + dz25, dz61, dz50);

#undef TRIPLE_PRODUCT

  volume *= twelveth;

  return volume ;
}

/******************************************/

//inline
Real_t CalcElemVolume( const Real_t x[8], const Real_t y[8], const Real_t z[8] )
{
return CalcElemVolume( x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7],
                       y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7],
                       z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7]);
}

/******************************************/

static inline
Real_t AreaFace( const Real_t x0, const Real_t x1,
                 const Real_t x2, const Real_t x3,
                 const Real_t y0, const Real_t y1,
                 const Real_t y2, const Real_t y3,
                 const Real_t z0, const Real_t z1,
                 const Real_t z2, const Real_t z3)
{
   Real_t fx = (x2 - x0) - (x3 - x1);
   Real_t fy = (y2 - y0) - (y3 - y1);
   Real_t fz = (z2 - z0) - (z3 - z1);
   Real_t gx = (x2 - x0) + (x3 - x1);
   Real_t gy = (y2 - y0) + (y3 - y1);
   Real_t gz = (z2 - z0) + (z3 - z1);
   Real_t area =
      (fx * fx + fy * fy + fz * fz) *
      (gx * gx + gy * gy + gz * gz) -
      (fx * gx + fy * gy + fz * gz) *
      (fx * gx + fy * gy + fz * gz);
   return area ;
}

/******************************************/

static inline
Real_t CalcElemCharacteristicLength( const Real_t x[8],
                                     const Real_t y[8],
                                     const Real_t z[8],
                                     const Real_t volume)
{
   Real_t a, charLength = Real_t(0.0);

   a = AreaFace(x[0],x[1],x[2],x[3],
                y[0],y[1],y[2],y[3],
                z[0],z[1],z[2],z[3]) ;
   charLength = std::max(a,charLength) ;

   a = AreaFace(x[4],x[5],x[6],x[7],
                y[4],y[5],y[6],y[7],
                z[4],z[5],z[6],z[7]) ;
   charLength = std::max(a,charLength) ;

   a = AreaFace(x[0],x[1],x[5],x[4],
                y[0],y[1],y[5],y[4],
                z[0],z[1],z[5],z[4]) ;
   charLength = std::max(a,charLength) ;

   a = AreaFace(x[1],x[2],x[6],x[5],
                y[1],y[2],y[6],y[5],
                z[1],z[2],z[6],z[5]) ;
   charLength = std::max(a,charLength) ;

   a = AreaFace(x[2],x[3],x[7],x[6],
                y[2],y[3],y[7],y[6],
                z[2],z[3],z[7],z[6]) ;
   charLength = std::max(a,charLength) ;

   a = AreaFace(x[3],x[0],x[4],x[7],
                y[3],y[0],y[4],y[7],
                z[3],z[0],z[4],z[7]) ;
   charLength = std::max(a,charLength) ;

   charLength = Real_t(4.0) * volume / SQRT(charLength);

   return charLength;
}

/******************************************/

static inline
void CalcElemVelocityGradient( const Real_t* const xvel,
                                const Real_t* const yvel,
                                const Real_t* const zvel,
                                const Real_t b[][8],
                                const Real_t detJ,
                                Real_t* const d )
{
  const Real_t inv_detJ = Real_t(1.0) / detJ ;
  Real_t dyddx, dxddy, dzddx, dxddz, dzddy, dyddz;
  const Real_t* const pfx = b[0];
  const Real_t* const pfy = b[1];
  const Real_t* const pfz = b[2];

  d[0] = inv_detJ * ( pfx[0] * (xvel[0]-xvel[6])
                     + pfx[1] * (xvel[1]-xvel[7])
                     + pfx[2] * (xvel[2]-xvel[4])
                     + pfx[3] * (xvel[3]-xvel[5]) );

  d[1] = inv_detJ * ( pfy[0] * (yvel[0]-yvel[6])
                     + pfy[1] * (yvel[1]-yvel[7])
                     + pfy[2] * (yvel[2]-yvel[4])
                     + pfy[3] * (yvel[3]-yvel[5]) );

  d[2] = inv_detJ * ( pfz[0] * (zvel[0]-zvel[6])
                     + pfz[1] * (zvel[1]-zvel[7])
                     + pfz[2] * (zvel[2]-zvel[4])
                     + pfz[3] * (zvel[3]-zvel[5]) );

  dyddx  = inv_detJ * ( pfx[0] * (yvel[0]-yvel[6])
                      + pfx[1] * (yvel[1]-yvel[7])
                      + pfx[2] * (yvel[2]-yvel[4])
                      + pfx[3] * (yvel[3]-yvel[5]) );

  dxddy  = inv_detJ * ( pfy[0] * (xvel[0]-xvel[6])
                      + pfy[1] * (xvel[1]-xvel[7])
                      + pfy[2] * (xvel[2]-xvel[4])
                      + pfy[3] * (xvel[3]-xvel[5]) );

  dzddx  = inv_detJ * ( pfx[0] * (zvel[0]-zvel[6])
                      + pfx[1] * (zvel[1]-zvel[7])
                      + pfx[2] * (zvel[2]-zvel[4])
                      + pfx[3] * (zvel[3]-zvel[5]) );

  dxddz  = inv_detJ * ( pfz[0] * (xvel[0]-xvel[6])
                      + pfz[1] * (xvel[1]-xvel[7])
                      + pfz[2] * (xvel[2]-xvel[4])
                      + pfz[3] * (xvel[3]-xvel[5]) );

  dzddy  = inv_detJ * ( pfy[0] * (zvel[0]-zvel[6])
                      + pfy[1] * (zvel[1]-zvel[7])
                      + pfy[2] * (zvel[2]-zvel[4])
                      + pfy[3] * (zvel[3]-zvel[5]) );

  dyddz  = inv_detJ * ( pfz[0] * (yvel[0]-yvel[6])
                      + pfz[1] * (yvel[1]-yvel[7])
                      + pfz[2] * (yvel[2]-yvel[4])
                      + pfz[3] * (yvel[3]-yvel[5]) );
  d[5]  = Real_t( .5) * ( dxddy + dyddx );
  d[4]  = Real_t( .5) * ( dxddz + dzddx );
  d[3]  = Real_t( .5) * ( dzddy + dyddz );
}

/******************************************/

//static inline
void CalcKinematicsForElems( Domain &domain,
                             Real_t deltaTime, Index_t numElem, op_dat p_dxx, op_dat p_dyy, op_dat p_dzz )
{

   op_par_loop(CalcKinematicsForElem, "CalcKinematicsForElem", elems,
               op_arg_dat(p_x, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_x, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_x, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_y, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_y, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_y, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_z, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_z, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_z, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_xd, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_xd, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_xd, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_yd, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_yd, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_yd, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_zd, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_zd, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_zd, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_dxx, -1, OP_ID, 1, "double", OP_WRITE), op_arg_dat(p_dyy, -1, OP_ID, 1, "double", OP_WRITE), op_arg_dat(p_dzz, -1, OP_ID, 1, "double", OP_WRITE), 
               op_arg_dat(p_vnew, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_volo, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_delv, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_v, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_arealg, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_gbl(&deltaTime, 1, "double", OP_READ)
   );

// //   // loop over all elements
// #pragma omp parallel for firstprivate(numElem, deltaTime)
//   for( Index_t k=0 ; k<numElem ; ++k )
//   {
//     Real_t B[3][8] ; /** shape function derivatives */
//     Real_t D[6] ;
//     Real_t x_local[8] ;
//     Real_t y_local[8] ;
//     Real_t z_local[8] ;
//     Real_t xd_local[8] ;
//     Real_t yd_local[8] ;
//     Real_t zd_local[8] ;
//     Real_t detJ = Real_t(0.0) ;

//     Real_t volume ;
//     Real_t relativeVolume ;
//    //  const Index_t* const elemToNode = domain.nodelist(k) ;
//     const Index_t* const elemToNode = &nodelist[k*Index_t(8)] ;

//     // get nodal coordinates from global arrays and copy into local arrays.
//     CollectDomainNodesToElemNodes(domain, elemToNode, x_local, y_local, z_local);

//     // volume calculations
//     volume = CalcElemVolume(x_local, y_local, z_local );
//     relativeVolume = volume / volo[k] ;
//    //  domain.vnew(k) = relativeVolume ;
//     vnew[k] = relativeVolume ;
//    //  domain.delv(k) = relativeVolume - domain.v(k) ;
//     delv[k] = relativeVolume - v[k] ;

//     // set characteristic length
//    //  domain.arealg(k) = CalcElemCharacteristicLength(x_local, y_local, z_local,
//    //                                           volume);
//     arealg[k] = CalcElemCharacteristicLength(x_local, y_local, z_local,
//                                              volume);

//     // get nodal velocities from global array and copy into local arrays.
//     for( Index_t lnode=0 ; lnode<8 ; ++lnode )
//     {
//       Index_t gnode = elemToNode[lnode];
//       // xd_local[lnode] = domain.xd(gnode);
//       xd_local[lnode] = xd[gnode];
//       // yd_local[lnode] = domain.yd(gnode);
//       yd_local[lnode] = yd[gnode];
//       // zd_local[lnode] = domain.zd(gnode);
//       zd_local[lnode] = zd[gnode];
//     }

//     Real_t dt2 = Real_t(0.5) * deltaTime;
//     for ( Index_t j=0 ; j<8 ; ++j )
//     {
//        x_local[j] -= dt2 * xd_local[j];
//        y_local[j] -= dt2 * yd_local[j];
//        z_local[j] -= dt2 * zd_local[j];
//     }

//     CalcElemShapeFunctionDerivatives( x_local, y_local, z_local,
//                                       B, &detJ );

//     CalcElemVelocityGradient( xd_local, yd_local, zd_local,
//                                B, detJ, D );

//     // put velocity gradient quantities into their global arrays.
//    //  domain.dxx(k) = D[0];
//    //  domain.dyy(k) = D[1];
//    //  domain.dzz(k) = D[2];
//     dxx[k] = D[0];
//     dyy[k] = D[1];
//     dzz[k] = D[2];
//   }
}

/******************************************/

static inline
void CalcLagrangeElements(Domain& domain)
{
   // Index_t numElem = domain.numElem() ;
   Index_t numElem = m_numElem ;
   if (numElem > 0) {
      // const Real_t deltatime = domain.deltatime() ;
      const Real_t deltatime = m_deltatime ;

      // domain.AllocateStrains(numElem);
      // AllocStrains(numElem);
      // op_dat p_dxx = op_decl_dat_temp(elems, 1, "double", dxx, "p_dxx");
      // op_dat p_dyy = op_decl_dat_temp(elems, 1, "double", dyy, "p_dyy");
      // op_dat p_dzz = op_decl_dat_temp(elems, 1, "double", dzz, "p_dzz");

      CalcKinematicsForElems(domain, deltatime, numElem, p_dxx, p_dyy, p_dzz) ;

      op_par_loop(CalcLagrangeElemRemaining, "CalcLagrangeElemRemaining", elems,
                  op_arg_dat(p_dxx, -1, OP_ID, 1, "double", OP_RW), op_arg_dat(p_dyy, -1, OP_ID, 1, "double", OP_RW), op_arg_dat(p_dzz, -1, OP_ID, 1, "double", OP_RW), 
                  op_arg_dat(p_vdov, -1, OP_ID, 1, "double", OP_WRITE),
                  op_arg_dat(p_vnew, -1, OP_ID, 1, "double", OP_READ)
                  );
      // element loop to do some stuff not included in the elemlib function.
// #pragma omp parallel for firstprivate(numElem)
//       for ( Index_t k=0 ; k<numElem ; ++k )
//       {
//          // calc strain rate and apply as constraint (only done in FB element)
//          // Real_t vdov = domain.dxx(k) + domain.dyy(k) + domain.dzz(k) ;
//          Real_t vdov = dxx[k] + dyy[k] + dzz[k] ;
//          Real_t vdovthird = vdov/Real_t(3.0) ;

//          // make the rate of deformation tensor deviatoric
//          // domain.vdov(k) = vdov ;
//          m_vdov[k] = vdov ;
//          // domain.dxx(k) -= vdovthird ;
//          // domain.dyy(k) -= vdovthird ;
//          // domain.dzz(k) -= vdovthird ;

//          dxx[k] -= vdovthird ;
//          dyy[k] -= vdovthird ;
//          dzz[k] -= vdovthird ;

//         // See if any volumes are negative, and take appropriate action.
//          // if (domain.vnew(k) <= Real_t(0.0))
//          if (vnew[k] <= Real_t(0.0))
//         {
// #if USE_MPI           
//            MPI_Abort(MPI_COMM_WORLD, VolumeError) ;
// #else
//            exit(VolumeError);
// #endif
//         }
//       }
      // domain.DeallocateStrains();
      // op_free_dat_temp(p_dxx);
      // op_free_dat_temp(p_dyy);
      // op_free_dat_temp(p_dzz);
      // DeallocStrains();
   }
}

/******************************************/

static inline
void CalcMonotonicQGradientsForElems(Domain& domain)
{
   // Index_t numElem = domain.numElem();
   Index_t numElem = m_numElem;

   op_par_loop(CalcMonotonicQGradientsForElem, "CalcMonotonicQGradientsForElem", elems,
               op_arg_dat(p_x, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_x, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_x, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_y, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_y, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_y, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_z, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_z, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_z, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_xd, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_xd, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_xd, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_yd, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_yd, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_yd, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_zd, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_zd, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_zd, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_volo, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_vnew, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_delx_zeta, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_delv_zeta, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_delv_xi, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_delx_xi, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_delx_eta, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_delv_eta, -1, OP_ID, 1, "double", OP_WRITE)
               );
// #pragma omp parallel for firstprivate(numElem)
//    for (Index_t i = 0 ; i < numElem ; ++i ) {
//       const Real_t ptiny = Real_t(1.e-36) ;
//       Real_t ax,ay,az ;
//       Real_t dxv,dyv,dzv ;

//       // const Index_t *elemToNode = domain.nodelist(i);
//       const Index_t *elemToNode = &nodelist[i*Index_t(8)];
//       Index_t n0 = elemToNode[0] ;
//       Index_t n1 = elemToNode[1] ;
//       Index_t n2 = elemToNode[2] ;
//       Index_t n3 = elemToNode[3] ;
//       Index_t n4 = elemToNode[4] ;
//       Index_t n5 = elemToNode[5] ;
//       Index_t n6 = elemToNode[6] ;
//       Index_t n7 = elemToNode[7] ;

//       // Real_t x0 = domain.x(n0) ;
//       // Real_t x1 = domain.x(n1) ;
//       // Real_t x2 = domain.x(n2) ;
//       // Real_t x3 = domain.x(n3) ;
//       // Real_t x4 = domain.x(n4) ;
//       // Real_t x5 = domain.x(n5) ;
//       // Real_t x6 = domain.x(n6) ;
//       // Real_t x7 = domain.x(n7) ;
//       Real_t x0 = x[n0] ;
//       Real_t x1 = x[n1] ;
//       Real_t x2 = x[n2] ;
//       Real_t x3 = x[n3] ;
//       Real_t x4 = x[n4] ;
//       Real_t x5 = x[n5] ;
//       Real_t x6 = x[n6] ;
//       Real_t x7 = x[n7] ;

//       // Real_t y0 = domain.y(n0) ;
//       // Real_t y1 = domain.y(n1) ;
//       // Real_t y2 = domain.y(n2) ;
//       // Real_t y3 = domain.y(n3) ;
//       // Real_t y4 = domain.y(n4) ;
//       // Real_t y5 = domain.y(n5) ;
//       // Real_t y6 = domain.y(n6) ;
//       // Real_t y7 = domain.y(n7) ;
//       Real_t y0 = y[n0] ;
//       Real_t y1 = y[n1] ;
//       Real_t y2 = y[n2] ;
//       Real_t y3 = y[n3] ;
//       Real_t y4 = y[n4] ;
//       Real_t y5 = y[n5] ;
//       Real_t y6 = y[n6] ;
//       Real_t y7 = y[n7] ;

//       // Real_t z0 = domain.z(n0) ;
//       // Real_t z1 = domain.z(n1) ;
//       // Real_t z2 = domain.z(n2) ;
//       // Real_t z3 = domain.z(n3) ;
//       // Real_t z4 = domain.z(n4) ;
//       // Real_t z5 = domain.z(n5) ;
//       // Real_t z6 = domain.z(n6) ;
//       // Real_t z7 = domain.z(n7) ;
//       Real_t z0 = z[n0] ;
//       Real_t z1 = z[n1] ;
//       Real_t z2 = z[n2] ;
//       Real_t z3 = z[n3] ;
//       Real_t z4 = z[n4] ;
//       Real_t z5 = z[n5] ;
//       Real_t z6 = z[n6] ;
//       Real_t z7 = z[n7] ;

//       Real_t xv0 = xd[n0] ;
//       Real_t xv1 = xd[n1] ;
//       Real_t xv2 = xd[n2] ;
//       Real_t xv3 = xd[n3] ;
//       Real_t xv4 = xd[n4] ;
//       Real_t xv5 = xd[n5] ;
//       Real_t xv6 = xd[n6] ;
//       Real_t xv7 = xd[n7] ;

//       Real_t yv0 = yd[n0] ;
//       Real_t yv1 = yd[n1] ;
//       Real_t yv2 = yd[n2] ;
//       Real_t yv3 = yd[n3] ;
//       Real_t yv4 = yd[n4] ;
//       Real_t yv5 = yd[n5] ;
//       Real_t yv6 = yd[n6] ;
//       Real_t yv7 = yd[n7] ;

//       Real_t zv0 = zd[n0] ;
//       Real_t zv1 = zd[n1] ;
//       Real_t zv2 = zd[n2] ;
//       Real_t zv3 = zd[n3] ;
//       Real_t zv4 = zd[n4] ;
//       Real_t zv5 = zd[n5] ;
//       Real_t zv6 = zd[n6] ;
//       Real_t zv7 = zd[n7] ;

//       // Real_t vol = domain.volo(i)*domain.vnew(i) ;
//       Real_t vol = volo[i]*vnew[i] ;
//       Real_t norm = Real_t(1.0) / ( vol + ptiny ) ;

//       Real_t dxj = Real_t(-0.25)*((x0+x1+x5+x4) - (x3+x2+x6+x7)) ;
//       Real_t dyj = Real_t(-0.25)*((y0+y1+y5+y4) - (y3+y2+y6+y7)) ;
//       Real_t dzj = Real_t(-0.25)*((z0+z1+z5+z4) - (z3+z2+z6+z7)) ;

//       Real_t dxi = Real_t( 0.25)*((x1+x2+x6+x5) - (x0+x3+x7+x4)) ;
//       Real_t dyi = Real_t( 0.25)*((y1+y2+y6+y5) - (y0+y3+y7+y4)) ;
//       Real_t dzi = Real_t( 0.25)*((z1+z2+z6+z5) - (z0+z3+z7+z4)) ;

//       Real_t dxk = Real_t( 0.25)*((x4+x5+x6+x7) - (x0+x1+x2+x3)) ;
//       Real_t dyk = Real_t( 0.25)*((y4+y5+y6+y7) - (y0+y1+y2+y3)) ;
//       Real_t dzk = Real_t( 0.25)*((z4+z5+z6+z7) - (z0+z1+z2+z3)) ;

//       /* find delvk and delxk ( i cross j ) */

//       ax = dyi*dzj - dzi*dyj ;
//       ay = dzi*dxj - dxi*dzj ;
//       az = dxi*dyj - dyi*dxj ;

//       // domain.delx_zeta(i) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;
//       delx_zeta[i] = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

//       ax *= norm ;
//       ay *= norm ;
//       az *= norm ;

//       dxv = Real_t(0.25)*((xv4+xv5+xv6+xv7) - (xv0+xv1+xv2+xv3)) ;
//       dyv = Real_t(0.25)*((yv4+yv5+yv6+yv7) - (yv0+yv1+yv2+yv3)) ;
//       dzv = Real_t(0.25)*((zv4+zv5+zv6+zv7) - (zv0+zv1+zv2+zv3)) ;

//       // domain.delv_zeta(i) = ax*dxv + ay*dyv + az*dzv ;
//       delv_zeta[i] = ax*dxv + ay*dyv + az*dzv ;

//       /* find delxi and delvi ( j cross k ) */

//       ax = dyj*dzk - dzj*dyk ;
//       ay = dzj*dxk - dxj*dzk ;
//       az = dxj*dyk - dyj*dxk ;

//       // domain.delx_xi(i) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;
//       delx_xi[i] = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

//       ax *= norm ;
//       ay *= norm ;
//       az *= norm ;

//       dxv = Real_t(0.25)*((xv1+xv2+xv6+xv5) - (xv0+xv3+xv7+xv4)) ;
//       dyv = Real_t(0.25)*((yv1+yv2+yv6+yv5) - (yv0+yv3+yv7+yv4)) ;
//       dzv = Real_t(0.25)*((zv1+zv2+zv6+zv5) - (zv0+zv3+zv7+zv4)) ;

//       // domain.delv_xi(i) = ax*dxv + ay*dyv + az*dzv ;
//       delv_xi[i] = ax*dxv + ay*dyv + az*dzv ;

//       /* find delxj and delvj ( k cross i ) */

//       ax = dyk*dzi - dzk*dyi ;
//       ay = dzk*dxi - dxk*dzi ;
//       az = dxk*dyi - dyk*dxi ;

//       // domain.delx_eta(i) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;
//       delx_eta[i] = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

//       ax *= norm ;
//       ay *= norm ;
//       az *= norm ;

//       dxv = Real_t(-0.25)*((xv0+xv1+xv5+xv4) - (xv3+xv2+xv6+xv7)) ;
//       dyv = Real_t(-0.25)*((yv0+yv1+yv5+yv4) - (yv3+yv2+yv6+yv7)) ;
//       dzv = Real_t(-0.25)*((zv0+zv1+zv5+zv4) - (zv3+zv2+zv6+zv7)) ;

//       // domain.delv_eta(i) = ax*dxv + ay*dyv + az*dzv ;
//       delv_eta[i] = ax*dxv + ay*dyv + az*dzv ;
//    }
}

/******************************************/

static inline
void CalcMonotonicQRegionForElems(Domain &domain, Int_t r,
                                  Real_t ptiny)
{
   // Real_t monoq_limiter_mult = domain.monoq_limiter_mult();
   // Real_t monoq_max_slope = domain.monoq_max_slope();
   // Real_t qlc_monoq = domain.qlc_monoq();
   // Real_t qqc_monoq = domain.qqc_monoq();
   Real_t monoq_limiter_mult = m_monoq_limiter_mult;
   Real_t monoq_max_slope =m_monoq_max_slope;
   Real_t qlc_monoq = m_qlc_monoq;
   Real_t qqc_monoq = m_qqc_monoq;

#pragma omp parallel for firstprivate(qlc_monoq, qqc_monoq, monoq_limiter_mult, monoq_max_slope, ptiny)
   // for ( Index_t i = 0 ; i < domain.regElemSize(r); ++i ) {
   for ( Index_t i = 0 ; i < m_regElemSize[r]; ++i ) {
      // Index_t ielem = domain.regElemlist(r,i);
      Index_t ielem = m_regElemlist[r][i];
      Real_t qlin, qquad ;
      Real_t phixi, phieta, phizeta ;
      // Int_t bcMask = domain.elemBC(ielem) ;
      Int_t bcMask = elemBC[ielem] ;
      Real_t delvm = 0.0, delvp =0.0;

      /*  phixi     */
      // Real_t norm = Real_t(1.) / (domain.delv_xi(ielem)+ ptiny ) ;
      Real_t norm = Real_t(1.) / (delv_xi[ielem]+ ptiny ) ;

      switch (bcMask & XI_M) {
         case XI_M_COMM: /* needs comm data */
         case 0:         delvm = delv_xi[lxim[ielem]]; break ;
         case XI_M_SYMM: delvm = delv_xi[ielem] ;       break ;
         case XI_M_FREE: delvm = Real_t(0.0) ;      break ;
         default:          fprintf(stderr, "Error in switch at %s line %d\n",
                                   __FILE__, __LINE__);
            delvm = 0; /* ERROR - but quiets the compiler */
            break;
      }
      switch (bcMask & XI_P) {
         case XI_P_COMM: /* needs comm data */
         case 0:         delvp = delv_xi[lxip[ielem]] ; break ;
         case XI_P_SYMM: delvp = delv_xi[ielem] ;       break ;
         case XI_P_FREE: delvp = Real_t(0.0) ;      break ;
         default:          fprintf(stderr, "Error in switch at %s line %d\n",
                                   __FILE__, __LINE__);
            delvp = 0; /* ERROR - but quiets the compiler */
            break;
      }

      delvm = delvm * norm ;
      delvp = delvp * norm ;

      phixi = Real_t(.5) * ( delvm + delvp ) ;

      delvm *= monoq_limiter_mult ;
      delvp *= monoq_limiter_mult ;

      if ( delvm < phixi ) phixi = delvm ;
      if ( delvp < phixi ) phixi = delvp ;
      if ( phixi < Real_t(0.)) phixi = Real_t(0.) ;
      if ( phixi > monoq_max_slope) phixi = monoq_max_slope;


      /*  phieta     */
      // norm = Real_t(1.) / ( domain.delv_eta(ielem) + ptiny ) ;
      norm = Real_t(1.) / ( delv_eta[ielem] + ptiny ) ;

      switch (bcMask & ETA_M) {
         case ETA_M_COMM: /* needs comm data */
         case 0:          delvm = delv_eta[letam[ielem]] ; break ;
         case ETA_M_SYMM: delvm = delv_eta[ielem] ;        break ;
         case ETA_M_FREE: delvm = Real_t(0.0) ;        break ;
         default:          fprintf(stderr, "Error in switch at %s line %d\n",
                                   __FILE__, __LINE__);
            delvm = 0; /* ERROR - but quiets the compiler */
            break;
      }
      switch (bcMask & ETA_P) {
         case ETA_P_COMM: /* needs comm data */
         case 0:          delvp = delv_eta[letap[ielem]] ; break ;
         case ETA_P_SYMM: delvp = delv_eta[ielem] ;        break ;
         case ETA_P_FREE: delvp = Real_t(0.0) ;        break ;
         default:          fprintf(stderr, "Error in switch at %s line %d\n",
                                   __FILE__, __LINE__);
            delvp = 0; /* ERROR - but quiets the compiler */
            break;
      }

      delvm = delvm * norm ;
      delvp = delvp * norm ;

      phieta = Real_t(.5) * ( delvm + delvp ) ;

      delvm *= monoq_limiter_mult ;
      delvp *= monoq_limiter_mult ;

      if ( delvm  < phieta ) phieta = delvm ;
      if ( delvp  < phieta ) phieta = delvp ;
      if ( phieta < Real_t(0.)) phieta = Real_t(0.) ;
      if ( phieta > monoq_max_slope)  phieta = monoq_max_slope;

      /*  phizeta     */
      // norm = Real_t(1.) / ( domain.delv_zeta(ielem) + ptiny ) ;
      norm = Real_t(1.) / ( delv_zeta[ielem] + ptiny ) ;

      switch (bcMask & ZETA_M) {
         case ZETA_M_COMM: /* needs comm data */
         // case 0:           delvm = domain.delv_zeta(domain.lzetam(ielem)) ; break ;
         case 0:           delvm = delv_zeta[lzetam[ielem]] ; break ;
         case ZETA_M_SYMM: delvm = delv_zeta[ielem] ;         break ;
         case ZETA_M_FREE: delvm = Real_t(0.0) ;          break ;
         default:          fprintf(stderr, "Error in switch at %s line %d\n",
                                   __FILE__, __LINE__);
            delvm = 0; /* ERROR - but quiets the compiler */
            break;
      }
      switch (bcMask & ZETA_P) {
         case ZETA_P_COMM: /* needs comm data */
         // case 0:           delvp = domain.delv_zeta(domain.lzetap(ielem)) ; break ;
         case 0:           delvp = delv_zeta[lzetap[ielem]] ; break ;
         case ZETA_P_SYMM: delvp = delv_zeta[ielem] ;         break ;
         case ZETA_P_FREE: delvp = Real_t(0.0) ;          break ;
         default:          fprintf(stderr, "Error in switch at %s line %d\n",
                                   __FILE__, __LINE__);
            delvp = 0; /* ERROR - but quiets the compiler */
            break;
      }

      delvm = delvm * norm ;
      delvp = delvp * norm ;

      phizeta = Real_t(.5) * ( delvm + delvp ) ;

      delvm *= monoq_limiter_mult ;
      delvp *= monoq_limiter_mult ;

      if ( delvm   < phizeta ) phizeta = delvm ;
      if ( delvp   < phizeta ) phizeta = delvp ;
      if ( phizeta < Real_t(0.)) phizeta = Real_t(0.);
      if ( phizeta > monoq_max_slope  ) phizeta = monoq_max_slope;

      /* Remove length scale */

      // if ( domain.vdov(ielem) > Real_t(0.) )  {
      if ( m_vdov[ielem] > Real_t(0.) )  {
         qlin  = Real_t(0.) ;
         qquad = Real_t(0.) ;
      }
      else {
         // Real_t delvxxi   = domain.delv_xi(ielem)   * domain.delx_xi(ielem)   ;
         // Real_t delvxeta  = domain.delv_eta(ielem)  * domain.delx_eta(ielem)  ;
         // Real_t delvxzeta = domain.delv_zeta(ielem) * domain.delx_zeta(ielem) ;
         Real_t delvxxi   = delv_xi[ielem]   * delx_xi[ielem]   ;
         Real_t delvxeta  = delv_eta[ielem]  * delx_eta[ielem]  ;
         Real_t delvxzeta = delv_zeta[ielem] * delx_zeta[ielem] ;

         if ( delvxxi   > Real_t(0.) ) delvxxi   = Real_t(0.) ;
         if ( delvxeta  > Real_t(0.) ) delvxeta  = Real_t(0.) ;
         if ( delvxzeta > Real_t(0.) ) delvxzeta = Real_t(0.) ;

         // Real_t rho = domain.elemMass(ielem) / (domain.volo(ielem) * domain.vnew(ielem)) ;
         Real_t rho = elemMass[ielem] / (volo[ielem] * vnew[ielem]) ;

         qlin = -qlc_monoq * rho *
            (  delvxxi   * (Real_t(1.) - phixi) +
               delvxeta  * (Real_t(1.) - phieta) +
               delvxzeta * (Real_t(1.) - phizeta)  ) ;

         qquad = qqc_monoq * rho *
            (  delvxxi*delvxxi     * (Real_t(1.) - phixi*phixi) +
               delvxeta*delvxeta   * (Real_t(1.) - phieta*phieta) +
               delvxzeta*delvxzeta * (Real_t(1.) - phizeta*phizeta)  ) ;
      }

      // domain.qq(ielem) = qquad ;
      // domain.ql(ielem) = qlin  ;
      qq[ielem] = qquad ;
      ql[ielem] = qlin  ;
   }
}

/******************************************/

static inline
void CalcMonotonicQForElems(Domain& domain)
{  
   //
   // initialize parameters
   // 
   const Real_t ptiny = Real_t(1.e-36) ;

   //
   // calculate the monotonic q for all regions
   //
   // for (Index_t r=0 ; r<domain.numReg() ; ++r) {
   for (Index_t r=0 ; r<m_numReg ; ++r) {
      // if (domain.regElemSize(r) > 0) {
      if (m_regElemSize[r] > 0) {
         CalcMonotonicQRegionForElems(domain, r, ptiny) ;
      }
   }
}

/******************************************/

static inline
void CalcQForElems(Domain& domain)
{
   //
   // MONOTONIC Q option
   //

   // Index_t numElem = domain.numElem() ;
   Index_t numElem = m_numElem ;

   if (numElem != 0) {
      Int_t allElem = numElem +  /* local elem */
            2*sizeX*sizeY + /* plane ghosts */
            2*sizeX*sizeZ + /* row ghosts */
            2*sizeY*sizeZ ; /* col ghosts */
            // 2*domain.sizeX()*domain.sizeY() + /* plane ghosts */
            // 2*domain.sizeX()*domain.sizeZ() + /* row ghosts */
            // 2*domain.sizeY()*domain.sizeZ() ; /* col ghosts */

      // domain.AllocateGradients(numElem, allElem);
      // AllocGradients(numElem, allElem);

#if USE_MPI      
      CommRecv(domain, MSG_MONOQ, 3,
               domain.sizeX(), domain.sizeY(), domain.sizeZ(),
               true, true) ;
#endif      

      /* Calculate velocity gradients */
      CalcMonotonicQGradientsForElems(domain);

#if USE_MPI      
      Domain_member fieldData[3] ;
      
      /* Transfer veloctiy gradients in the first order elements */
      /* problem->commElements->Transfer(CommElements::monoQ) ; */

      fieldData[0] = &Domain::delv_xi ;
      fieldData[1] = &Domain::delv_eta ;
      fieldData[2] = &Domain::delv_zeta ;

      CommSend(domain, MSG_MONOQ, 3, fieldData,
               domain.sizeX(), domain.sizeY(), domain.sizeZ(),
               true, true) ;

      CommMonoQ(domain) ;
#endif      
      CalcMonotonicQForElems(domain);
      // Free up memory
      // domain.DeallocateGradients();
      // DeallocGradients();

      op_par_loop(NoExcessiveArtificialViscosity, "NoExcessiveArtificialViscosity", elems,
                  op_arg_dat(p_q, -1, OP_ID, 1, "double", OP_READ));
      /* Don't allow excessive artificial viscosity */
//       Index_t idx = -1; 
//       for (Index_t i=0; i<numElem; ++i) {
//          // if ( domain.q(i) > domain.qstop() ) {
//          if ( q[i] > m_qstop ) {
//             idx = i ;
//             break ;
//          }
//       }

//       if(idx >= 0) {
// #if USE_MPI         
//          MPI_Abort(MPI_COMM_WORLD, QStopError) ;
// #else
//          exit(QStopError);
// #endif
//       }
   }
}

/******************************************/

static inline
void CalcPressureForElems(Real_t* p_new, Real_t* bvc,
                          Real_t* pbvc, Real_t* e_old,
                          Real_t* compression, Real_t *vnewc,
                          Real_t pmin,
                          Real_t p_cut, Real_t eosvmax,
                          Index_t length, Index_t *regElemList)
{
#pragma omp parallel for firstprivate(length)
   for (Index_t i = 0; i < length ; ++i) {
      Real_t c1s = Real_t(2.0)/Real_t(3.0) ;
      bvc[i] = c1s * (compression[i] + Real_t(1.));
      pbvc[i] = c1s;
   }

#pragma omp parallel for firstprivate(length, pmin, p_cut, eosvmax)
   for (Index_t i = 0 ; i < length ; ++i){
      Index_t ielem = regElemList[i];
      
      p_new[i] = bvc[i] * e_old[i] ;

      if    (FABS(p_new[i]) <  p_cut   )
         p_new[i] = Real_t(0.0) ;

      if    ( vnewc[ielem] >= eosvmax ) /* impossible condition here? */
         p_new[i] = Real_t(0.0) ;

      if    (p_new[i]       <  pmin)
         p_new[i]   = pmin ;
   }
}

/******************************************/

static inline
void CalcEnergyForElems(Real_t* p_new, Real_t* e_new, Real_t* q_new,
                        Real_t* bvc, Real_t* pbvc,
                        Real_t* p_old, Real_t* e_old, Real_t* q_old,
                        Real_t* compression, Real_t* compHalfStep,
                        Real_t* vnewc, Real_t* work, Real_t* delvc, Real_t pmin,
                        Real_t p_cut, Real_t  e_cut, Real_t q_cut, Real_t emin,
                        Real_t* qq_old, Real_t* ql_old,
                        Real_t rho0,
                        Real_t eosvmax,
                        Index_t length, Index_t *regElemList)
{
   Real_t *pHalfStep = Allocate<Real_t>(length) ;

#pragma omp parallel for firstprivate(length, emin)
   for (Index_t i = 0 ; i < length ; ++i) {
      e_new[i] = e_old[i] - Real_t(0.5) * delvc[i] * (p_old[i] + q_old[i])
         + Real_t(0.5) * work[i];

      if (e_new[i]  < emin ) {
         e_new[i] = emin ;
      }
   }

   CalcPressureForElems(pHalfStep, bvc, pbvc, e_new, compHalfStep, vnewc,
                        pmin, p_cut, eosvmax, length, regElemList);

#pragma omp parallel for firstprivate(length, rho0)
   for (Index_t i = 0 ; i < length ; ++i) {
      Real_t vhalf = Real_t(1.) / (Real_t(1.) + compHalfStep[i]) ;

      if ( delvc[i] > Real_t(0.) ) {
         q_new[i] /* = qq_old[i] = ql_old[i] */ = Real_t(0.) ;
      }
      else {
         Real_t ssc = ( pbvc[i] * e_new[i]
                 + vhalf * vhalf * bvc[i] * pHalfStep[i] ) / rho0 ;

         if ( ssc <= Real_t(.1111111e-36) ) {
            ssc = Real_t(.3333333e-18) ;
         } else {
            ssc = SQRT(ssc) ;
         }

         q_new[i] = (ssc*ql_old[i] + qq_old[i]) ;
      }

      e_new[i] = e_new[i] + Real_t(0.5) * delvc[i]
         * (  Real_t(3.0)*(p_old[i]     + q_old[i])
              - Real_t(4.0)*(pHalfStep[i] + q_new[i])) ;
   }

#pragma omp parallel for firstprivate(length, emin, e_cut)
   for (Index_t i = 0 ; i < length ; ++i) {

      e_new[i] += Real_t(0.5) * work[i];

      if (FABS(e_new[i]) < e_cut) {
         e_new[i] = Real_t(0.)  ;
      }
      if (     e_new[i]  < emin ) {
         e_new[i] = emin ;
      }
   }

   CalcPressureForElems(p_new, bvc, pbvc, e_new, compression, vnewc,
                        pmin, p_cut, eosvmax, length, regElemList);

#pragma omp parallel for firstprivate(length, rho0, emin, e_cut)
   for (Index_t i = 0 ; i < length ; ++i){
      const Real_t sixth = Real_t(1.0) / Real_t(6.0) ;
      Index_t ielem = regElemList[i];
      Real_t q_tilde ;

      if (delvc[i] > Real_t(0.)) {
         q_tilde = Real_t(0.) ;
      }
      else {
         Real_t ssc = ( pbvc[i] * e_new[i]
                 + vnewc[ielem] * vnewc[ielem] * bvc[i] * p_new[i] ) / rho0 ;

         if ( ssc <= Real_t(.1111111e-36) ) {
            ssc = Real_t(.3333333e-18) ;
         } else {
            ssc = SQRT(ssc) ;
         }

         q_tilde = (ssc*ql_old[i] + qq_old[i]) ;
      }

      e_new[i] = e_new[i] - (  Real_t(7.0)*(p_old[i]     + q_old[i])
                               - Real_t(8.0)*(pHalfStep[i] + q_new[i])
                               + (p_new[i] + q_tilde)) * delvc[i]*sixth ;

      if (FABS(e_new[i]) < e_cut) {
         e_new[i] = Real_t(0.)  ;
      }
      if (     e_new[i]  < emin ) {
         e_new[i] = emin ;
      }
   }

   CalcPressureForElems(p_new, bvc, pbvc, e_new, compression, vnewc,
                        pmin, p_cut, eosvmax, length, regElemList);

#pragma omp parallel for firstprivate(length, rho0, q_cut)
   for (Index_t i = 0 ; i < length ; ++i){
      Index_t ielem = regElemList[i];

      if ( delvc[i] <= Real_t(0.) ) {
         Real_t ssc = ( pbvc[i] * e_new[i]
                 + vnewc[ielem] * vnewc[ielem] * bvc[i] * p_new[i] ) / rho0 ;

         if ( ssc <= Real_t(.1111111e-36) ) {
            ssc = Real_t(.3333333e-18) ;
         } else {
            ssc = SQRT(ssc) ;
         }

         q_new[i] = (ssc*ql_old[i] + qq_old[i]) ;

         if (FABS(q_new[i]) < q_cut) q_new[i] = Real_t(0.) ;
      }
   }

   Release(&pHalfStep) ;

   return ;
}

/******************************************/

static inline
void CalcSoundSpeedForElems(Domain &domain,
                            Real_t *vnewc, Real_t rho0, Real_t *enewc,
                            Real_t *pnewc, Real_t *pbvc,
                            Real_t *bvc, Real_t ss4o3,
                            Index_t len, Index_t *regElemList)
{
#pragma omp parallel for firstprivate(rho0, ss4o3)
   for (Index_t i = 0; i < len ; ++i) {
      Index_t ielem = regElemList[i];
      Real_t ssTmp = (pbvc[i] * enewc[i] + vnewc[ielem] * vnewc[ielem] *
                 bvc[i] * pnewc[i]) / rho0;
      if (ssTmp <= Real_t(.1111111e-36)) {
         ssTmp = Real_t(.3333333e-18);
      }
      else {
         ssTmp = SQRT(ssTmp);
      }
      // domain.ss(ielem) = ssTmp ;
      ss[ielem] = ssTmp ;

   }
}

/******************************************/

static inline
void EvalEOSForElems(Domain& domain, Real_t *vnewc,
                     Int_t numElemReg, Index_t *regElemList, Int_t rep)
{
   // Real_t  e_cut = domain.e_cut() ;
   // Real_t  p_cut = domain.p_cut() ;
   // Real_t  ss4o3 = domain.ss4o3() ;
   // Real_t  q_cut = domain.q_cut() ;

   // Real_t eosvmax = domain.eosvmax() ;
   // Real_t eosvmin = domain.eosvmin() ;
   // Real_t pmin    = domain.pmin() ;
   // Real_t emin    = domain.emin() ;
   // Real_t rho0    = domain.refdens() ;
   Real_t  e_cut = m_e_cut ;
   Real_t  p_cut = m_p_cut ;
   Real_t  ss4o3 = m_ss4o3 ;
   Real_t  q_cut = m_q_cut ;

   Real_t eosvmax = m_eosvmax ;
   Real_t eosvmin = m_eosvmin ;
   Real_t pmin    = m_pmin ;
   Real_t emin    = m_emin ;
   Real_t rho0    = m_refdens ;

   // These temporaries will be of different size for 
   // each call (due to different sized region element
   // lists)
   Real_t *e_old = Allocate<Real_t>(numElemReg) ;
   Real_t *delvc = Allocate<Real_t>(numElemReg) ;
   Real_t *p_old = Allocate<Real_t>(numElemReg) ;
   Real_t *q_old = Allocate<Real_t>(numElemReg) ;
   Real_t *compression = Allocate<Real_t>(numElemReg) ;
   Real_t *compHalfStep = Allocate<Real_t>(numElemReg) ;
   Real_t *qq_old = Allocate<Real_t>(numElemReg) ;
   Real_t *ql_old = Allocate<Real_t>(numElemReg) ;
   Real_t *work = Allocate<Real_t>(numElemReg) ;
   Real_t *p_new = Allocate<Real_t>(numElemReg) ;
   Real_t *e_new = Allocate<Real_t>(numElemReg) ;
   Real_t *q_new = Allocate<Real_t>(numElemReg) ;
   Real_t *bvc = Allocate<Real_t>(numElemReg) ;
   Real_t *pbvc = Allocate<Real_t>(numElemReg) ;
 
   //loop to add load imbalance based on region number 
   for(Int_t j = 0; j < rep; j++) {
      /* compress data, minimal set */
#pragma omp parallel
      {
#pragma omp for nowait firstprivate(numElemReg)
         for (Index_t i=0; i<numElemReg; ++i) {
            Index_t ielem = regElemList[i];
            // e_old[i] = domain.e(ielem) ;
            // delvc[i] = domain.delv(ielem) ;
            e_old[i] = e[ielem] ;
            delvc[i] = delv[ielem] ;
            // p_old[i] = domain.p(ielem) ;
            // q_old[i] = domain.q(ielem) ;
            
            p_old[i] = p[ielem] ;
            q_old[i] = q[ielem] ;

            // qq_old[i] = domain.qq(ielem) ;
            // ql_old[i] = domain.ql(ielem) ;

            qq_old[i] = qq[ielem] ;
            ql_old[i] = ql[ielem] ;
         }

#pragma omp for firstprivate(numElemReg)
         for (Index_t i = 0; i < numElemReg ; ++i) {
            Index_t ielem = regElemList[i];
            Real_t vchalf ;
            compression[i] = Real_t(1.) / vnewc[ielem] - Real_t(1.);
            vchalf = vnewc[ielem] - delvc[i] * Real_t(.5);
            compHalfStep[i] = Real_t(1.) / vchalf - Real_t(1.);
         }

      /* Check for v > eosvmax or v < eosvmin */
         if ( eosvmin != Real_t(0.) ) {
#pragma omp for nowait firstprivate(numElemReg, eosvmin)
            for(Index_t i=0 ; i<numElemReg ; ++i) {
               Index_t ielem = regElemList[i];
               if (vnewc[ielem] <= eosvmin) { /* impossible due to calling func? */
                  compHalfStep[i] = compression[i] ;
               }
            }
         }
         if ( eosvmax != Real_t(0.) ) {
#pragma omp for nowait firstprivate(numElemReg, eosvmax)
            for(Index_t i=0 ; i<numElemReg ; ++i) {
               Index_t ielem = regElemList[i];
               if (vnewc[ielem] >= eosvmax) { /* impossible due to calling func? */
                  p_old[i]        = Real_t(0.) ;
                  compression[i]  = Real_t(0.) ;
                  compHalfStep[i] = Real_t(0.) ;
               }
            }
         }

#pragma omp for nowait firstprivate(numElemReg)
         for (Index_t i = 0 ; i < numElemReg ; ++i) {
            work[i] = Real_t(0.) ; 
         }
      }
      CalcEnergyForElems(p_new, e_new, q_new, bvc, pbvc,
                         p_old, e_old,  q_old, compression, compHalfStep,
                         vnewc, work,  delvc, pmin,
                         p_cut, e_cut, q_cut, emin,
                         qq_old, ql_old, rho0, eosvmax,
                         numElemReg, regElemList);
   }

#pragma omp parallel for firstprivate(numElemReg)
   for (Index_t i=0; i<numElemReg; ++i) {
      Index_t ielem = regElemList[i];
      // domain.p(ielem) = p_new[i] ;
      // domain.e(ielem) = e_new[i] ;
      // domain.q(ielem) = q_new[i] ;
      p[ielem] = p_new[i];
      e[ielem] = e_new[i];
      q[ielem] = q_new[i];
   }

   CalcSoundSpeedForElems(domain,
                          vnewc, rho0, e_new, p_new,
                          pbvc, bvc, ss4o3,
                          numElemReg, regElemList) ;

   Release(&pbvc) ;
   Release(&bvc) ;
   Release(&q_new) ;
   Release(&e_new) ;
   Release(&p_new) ;
   Release(&work) ;
   Release(&ql_old) ;
   Release(&qq_old) ;
   Release(&compHalfStep) ;
   Release(&compression) ;
   Release(&q_old) ;
   Release(&p_old) ;
   Release(&delvc) ;
   Release(&e_old) ;
}

/******************************************/

static inline
void ApplyMaterialPropertiesForElems(Domain& domain)
{
   // Index_t numElem = domain.numElem() ;
   Index_t numElem = m_numElem ;

  if (numElem != 0) {
    /* Expose all of the variables needed for material evaluation */
   //  Real_t eosvmin = domain.eosvmin() ;
   //  Real_t eosvmax = domain.eosvmax() ;
    Real_t eosvmin = m_eosvmin ;
    Real_t eosvmax = m_eosvmax ;
    Real_t *vnewc = Allocate<Real_t>(numElem) ;

#pragma omp parallel
    {
#pragma omp for firstprivate(numElem)
       for(Index_t i=0 ; i<numElem ; ++i) {
         //  vnewc[i] = domain.vnew(i) ;
          vnewc[i] = vnew[i] ;
       }

       // Bound the updated relative volumes with eosvmin/max
       if (eosvmin != Real_t(0.)) {
#pragma omp for nowait firstprivate(numElem)
          for(Index_t i=0 ; i<numElem ; ++i) {
             if (vnewc[i] < eosvmin)
                vnewc[i] = eosvmin ;
          }
       }

       if (eosvmax != Real_t(0.)) {
#pragma omp for nowait firstprivate(numElem)
          for(Index_t i=0 ; i<numElem ; ++i) {
             if (vnewc[i] > eosvmax)
                vnewc[i] = eosvmax ;
          }
       }

       // This check may not make perfect sense in LULESH, but
       // it's representative of something in the full code -
       // just leave it in, please
#pragma omp for nowait firstprivate(numElem)
       for (Index_t i=0; i<numElem; ++i) {
         //  Real_t vc = domain.v(i) ;
          Real_t vc = v[i] ;
          if (eosvmin != Real_t(0.)) {
             if (vc < eosvmin)
                vc = eosvmin ;
          }
          if (eosvmax != Real_t(0.)) {
             if (vc > eosvmax)
                vc = eosvmax ;
          }
          if (vc <= 0.) {
#if USE_MPI
             MPI_Abort(MPI_COMM_WORLD, VolumeError) ;
#else
             exit(VolumeError);
#endif
          }
       }
    }

   //  for (Int_t r=0 ; r<domain.numReg() ; r++) {
    for (Int_t r=0 ; r<m_numReg ; r++) {
      //  Index_t numElemReg = domain.regElemSize(r);
       Index_t numElemReg = m_regElemSize[r];
      //  Index_t *regElemList = domain.regElemlist(r);
       Index_t *regElemList = m_regElemlist[r];
       Int_t rep;
       //Determine load imbalance for this region
       //round down the number with lowest cost
      //  if(r < domain.numReg()/2)
       if(r < m_numReg/2)
	 rep = 1;
       //you don't get an expensive region unless you at least have 5 regions
      //  else if(r < (domain.numReg() - (domain.numReg()+15)/20))
       else if(r < (m_numReg - (m_numReg+15)/20))
         // rep = 1 + domain.cost();
         rep = 1 + m_cost;
       //very expensive regions
       else
	//  rep = 10 * (1+ domain.cost());
	 rep = 10 * (1+ m_cost);
       EvalEOSForElems(domain, vnewc, numElemReg, regElemList, rep);
    }

    Release(&vnewc) ;
  }
}

/******************************************/

static inline
void UpdateVolumesForElems(Domain &domain,
                           Real_t v_cut, Index_t length)
{
   if (length != 0) {
      op_par_loop(updateVolumesForElem, "UpdateVolumesForElems", elems,
               op_arg_dat(p_vnew, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_v, -1, OP_ID, 1, "double", OP_WRITE)
               );
// #pragma omp parallel for firstprivate(length, v_cut)
//       for(Index_t i=0 ; i<length ; ++i) {
//          // Real_t tmpV = domain.vnew(i) ;
//          Real_t tmpV = vnew[i] ;

//          if ( FABS(tmpV - Real_t(1.0)) < v_cut )
//             tmpV = Real_t(1.0) ;

//          // domain.v(i) = tmpV ;
//          v[i] = tmpV ;
//       }
   }

   return ;
}

/******************************************/

static inline
void LagrangeElements(Domain& domain, Index_t numElem)
{
  CalcLagrangeElements(domain) ;

  /* Calculate Q.  (Monotonic q option requires communication) */
  CalcQForElems(domain) ;

  ApplyMaterialPropertiesForElems(domain) ;

//   UpdateVolumesForElems(domain, 
//                         domain.v_cut(), numElem) ;
  UpdateVolumesForElems(domain, 
                        m_v_cut, numElem) ;
}

/******************************************/

static inline
void CalcCourantConstraintForElems(Domain &domain, Index_t length,
                                   Index_t *regElemlist,
                                   Real_t qqc, Real_t& dtcourant)
{
#if _OPENMP
   const Index_t threads = omp_get_max_threads();
   Index_t courant_elem_per_thread[threads];
   Real_t dtcourant_per_thread[threads];
#else
   Index_t threads = 1;
   Index_t courant_elem_per_thread[1];
   Real_t  dtcourant_per_thread[1];
#endif

#pragma omp parallel firstprivate(length, qqc)
   {
      Real_t   qqc2 = Real_t(64.0) * qqc * qqc ;
      Real_t   dtcourant_tmp = dtcourant;
      Index_t  courant_elem  = -1 ;

#if _OPENMP
      Index_t thread_num = omp_get_thread_num();
#else
      Index_t thread_num = 0;
#endif      

#pragma omp for 
      for (Index_t i = 0 ; i < length ; ++i) {
         Index_t indx = regElemlist[i] ;
         // Real_t dtf = domain.ss(indx) * domain.ss(indx) ;
         Real_t dtf = ss[indx] * ss[indx] ;

         // if ( domain.vdov(indx) < Real_t(0.) ) {
         if ( m_vdov[indx] < Real_t(0.) ) {
            dtf = dtf
               //  + qqc2 * domain.arealg(indx) * domain.arealg(indx)
               //  * domain.vdov(indx) * domain.vdov(indx) ;
                + qqc2 * arealg[indx] * arealg[indx]
                * m_vdov[indx] * m_vdov[indx] ;
         }

         dtf = SQRT(dtf) ;
         dtf = arealg[indx] / dtf ;

         // if (domain.vdov(indx) != Real_t(0.)) {
         if (m_vdov[indx] != Real_t(0.)) {
            if ( dtf < dtcourant_tmp ) {
               dtcourant_tmp = dtf ;
               courant_elem  = indx ;
            }
         }
      }

      dtcourant_per_thread[thread_num]    = dtcourant_tmp ;
      courant_elem_per_thread[thread_num] = courant_elem ;
   }

   for (Index_t i = 1; i < threads; ++i) {
      if (dtcourant_per_thread[i] < dtcourant_per_thread[0] ) {
         dtcourant_per_thread[0]    = dtcourant_per_thread[i];
         courant_elem_per_thread[0] = courant_elem_per_thread[i];
      }
   }

   if (courant_elem_per_thread[0] != -1) {
      dtcourant = dtcourant_per_thread[0] ;
   }

   return ;

}

/******************************************/

static inline
void CalcHydroConstraintForElems(Domain &domain, Index_t length,
                                 Index_t *regElemlist, Real_t dvovmax, Real_t& dthydro)
{
#if _OPENMP
   const Index_t threads = omp_get_max_threads();
   Index_t hydro_elem_per_thread[threads];
   Real_t dthydro_per_thread[threads];
#else
   Index_t threads = 1;
   Index_t hydro_elem_per_thread[1];
   Real_t  dthydro_per_thread[1];
#endif

#pragma omp parallel firstprivate(length, dvovmax)
   {
      Real_t dthydro_tmp = dthydro ;
      Index_t hydro_elem = -1 ;

#if _OPENMP
      Index_t thread_num = omp_get_thread_num();
#else      
      Index_t thread_num = 0;
#endif      

#pragma omp for
      for (Index_t i = 0 ; i < length ; ++i) {
         Index_t indx = regElemlist[i] ;

         // if (domain.vdov(indx) != Real_t(0.)) {
         //    Real_t dtdvov = dvovmax / (FABS(domain.vdov(indx))+Real_t(1.e-20)) ;
         if (m_vdov[indx] != Real_t(0.)) {
            Real_t dtdvov = dvovmax / (FABS(m_vdov[indx])+Real_t(1.e-20)) ;

            if ( dthydro_tmp > dtdvov ) {
                  dthydro_tmp = dtdvov ;
                  hydro_elem = indx ;
            }
         }
      }

      dthydro_per_thread[thread_num]    = dthydro_tmp ;
      hydro_elem_per_thread[thread_num] = hydro_elem ;
   }

   for (Index_t i = 1; i < threads; ++i) {
      if(dthydro_per_thread[i] < dthydro_per_thread[0]) {
         dthydro_per_thread[0]    = dthydro_per_thread[i];
         hydro_elem_per_thread[0] =  hydro_elem_per_thread[i];
      }
   }

   if (hydro_elem_per_thread[0] != -1) {
      dthydro =  dthydro_per_thread[0] ;
   }

   return ;
}

/******************************************/

static inline
void CalcTimeConstraintsForElems(Domain& domain) {

   // Initialize conditions to a very large value
   // domain.dtcourant() = 1.0e+20;
   // domain.dthydro() = 1.0e+20;
   m_dtcourant = 1.0e+20;
   m_dthydro = 1.0e+20;

   // for (Index_t r=0 ; r < domain.numReg() ; ++r) {
   for (Index_t r=0 ; r < m_numReg ; ++r) {
      /* evaluate time constraint */
      // CalcCourantConstraintForElems(domain, domain.regElemSize(r),
                                    // domain.regElemlist(r),
                                    // domain.qqc(),
                                    // domain.dtcourant()) ;
      CalcCourantConstraintForElems(domain, m_regElemSize[r],
                                    m_regElemlist[r],
                                    m_qqc,
                                    m_dtcourant) ;

      /* check hydro constraint */
      // CalcHydroConstraintForElems(domain, domain.regElemSize(r),
                                 //  domain.regElemlist(r),
                                 //  domain.dvovmax(),
      CalcHydroConstraintForElems(domain, m_regElemSize[r],
                                  m_regElemlist[r],
                                  m_dvovmax,
                                  m_dthydro) ;
   }
}

/******************************************/

static inline
void LagrangeLeapFrog(Domain& domain)
{
#ifdef SEDOV_SYNC_POS_VEL_LATE
   Domain_member fieldData[6] ;
#endif

   /* calculate nodal forces, accelerations, velocities, positions, with
    * applied boundary conditions and slide surface considerations */
   LagrangeNodal(domain);


#ifdef SEDOV_SYNC_POS_VEL_LATE
#endif

   /* calculate element quantities (i.e. velocity gradient & q), and update
    * material states */
   // LagrangeElements(domain, domain.numElem());
   LagrangeElements(domain, m_numElem);

#if USE_MPI   
#ifdef SEDOV_SYNC_POS_VEL_LATE
   CommRecv(domain, MSG_SYNC_POS_VEL, 6,
            domain.sizeX() + 1, domain.sizeY() + 1, domain.sizeZ() + 1,
            false, false) ;

   fieldData[0] = &Domain::x ;
   fieldData[1] = &Domain::y ;
   fieldData[2] = &Domain::z ;
   fieldData[3] = &Domain::xd ;
   fieldData[4] = &Domain::yd ;
   fieldData[5] = &Domain::zd ;
   
   CommSend(domain, MSG_SYNC_POS_VEL, 6, fieldData,
            domain.sizeX() + 1, domain.sizeY() + 1, domain.sizeZ() + 1,
            false, false) ;
#endif
#endif   

   CalcTimeConstraintsForElems(domain);

#if USE_MPI   
#ifdef SEDOV_SYNC_POS_VEL_LATE
   CommSyncPosVel(domain) ;
#endif
#endif   
}


/******************************************/

int main(int argc, char *argv[])
{
   op_init(argc, argv, 2);
   Domain *locDom ;

   int numRanks ;
   int myRank ;
   struct cmdLineOpts opts;


#if USE_MPI   
   Domain_member fieldData ;
   
#ifdef _OPENMP
   int thread_support;

   MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &thread_support);
   if (thread_support==MPI_THREAD_SINGLE)
    {
        fprintf(stderr,"The MPI implementation has no support for threading\n");
        MPI_Finalize();
        exit(1);
    }
#else
   MPI_Init(&argc, &argv);
#endif
    
   MPI_Comm_size(MPI_COMM_WORLD, &numRanks) ;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;
#else
   numRanks = 1;
   myRank = 0;
#endif   

   /* Set defaults that can be overridden by command line opts */
   opts.its = 9999999; // Iterations
   opts.nx  = 30; //Size
   opts.numReg = 11;
   opts.numFiles = (int)(numRanks+10)/9;
   opts.showProg = 0;
   opts.quiet = 0;
   opts.viz = 0;
   opts.balance = 1;
   opts.cost = 1;

   ParseCommandLineOptions(argc, argv, myRank, &opts);
   std::cout << "Testing\n";

   if ((myRank == 0) && (opts.quiet == 0)) {
      std::cout << "Running problem size " << opts.nx << "^3 per domain until completion\n";
      std::cout << "Num processors: "      << numRanks << "\n";
#if _OPENMP
      std::cout << "Num threads: " << omp_get_max_threads() << "\n";
#endif
      std::cout << "Total number of elements: " << ((Int8_t)numRanks*opts.nx*opts.nx*opts.nx) << " \n\n";
      std::cout << "To run other sizes, use -s <integer>.\n";
      std::cout << "To run a fixed number of iterations, use -i <integer>.\n";
      std::cout << "To run a more or less balanced region set, use -b <integer>.\n";
      std::cout << "To change the relative costs of regions, use -c <integer>.\n";
      std::cout << "To print out progress, use -p\n";
      std::cout << "To write an output file for VisIt, use -v\n";
      std::cout << "See help (-h) for more options\n\n";
   }

   // Set up the mesh and decompose. Assumes regular cubes for now
   Int_t col, row, plane, side;
   InitMeshDecomp(numRanks, myRank, &col, &row, &plane, &side);

   // Build the main data structure and initialize it
   locDom = new Domain(numRanks, col, row, plane, opts.nx,
                       side, opts.numReg, opts.balance, opts.cost) ;
   
   int edgeElems = opts.nx;
   int edgeNodes = edgeElems+1;

   allocateElems(edgeElems);
   allocateNodes(0,0,0,edgeNodes);
   //! The 0 zeroes are for the row, col and plane locations, might have to be removed/changed
   initialise(0,0,0,opts.nx,1,opts.numReg,opts.balance, opts.cost);

   nodes = op_decl_set(m_numNode, "nodes");
   elems = op_decl_set(m_numElem, "elems");
   symmetry = op_decl_set(edgeNodes*edgeNodes, "symmetries");
   temp_vols = op_decl_set(m_numElem*8, "tempVols");

   p_nodelist = op_decl_map(elems, nodes, 8, nodelist, "nodelist");
   p_symmX = op_decl_map(nodes, symmetry, 1, symmX, "symmX");
   p_symmY = op_decl_map(nodes, symmetry, 1, symmY, "symmY");
   p_symmZ = op_decl_map(nodes, symmetry, 1, symmZ, "symmZ");

   //Node Centred
   p_x = op_decl_dat(nodes, 1, "double", x, "p_x");
   p_y = op_decl_dat(nodes, 1, "double", y, "p_y");
   p_z = op_decl_dat(nodes, 1, "double", z, "p_z");
   p_xd = op_decl_dat(nodes, 1, "double", xd, "p_xd");
   p_yd = op_decl_dat(nodes, 1, "double", yd, "p_yd");
   p_zd = op_decl_dat(nodes, 1, "double", zd, "p_zd");
   p_xdd = op_decl_dat(nodes, 1, "double", xdd, "p_xdd");
   p_ydd = op_decl_dat(nodes, 1, "double", ydd, "p_ydd");
   p_zdd = op_decl_dat(nodes, 1, "double", zdd, "p_zdd");
   p_fx = op_decl_dat(nodes, 1, "double", m_fx, "p_fx");
   p_fy = op_decl_dat(nodes, 1, "double", m_fy, "p_fy");
   p_fz = op_decl_dat(nodes, 1, "double", m_fz, "p_fz");

   p_nodalMass = op_decl_dat(nodes, 1, "double", nodalMass, "p_nodalMass");
   // //Elem Centred
   op_dat p_e = op_decl_dat(elems, 1, "double", e, "p_e");
   p_p = op_decl_dat(elems, 1, "double", p, "p_p");
   p_q = op_decl_dat(elems, 1, "double", q, "p_q");
   op_dat p_ql = op_decl_dat(elems, 1, "double", ql, "p_ql");
   op_dat p_qq = op_decl_dat(elems, 1, "double", qq, "p_qq");
   p_v = op_decl_dat(elems, 1, "double", v, "p_v");
   p_volo = op_decl_dat(elems, 1, "double", volo, "p_volo");
   p_delv = op_decl_dat(elems, 1, "double", delv, "p_delv");
   p_vdov = op_decl_dat(elems, 1, "double", m_vdov, "p_vdov");
   p_arealg = op_decl_dat(elems, 1, "double", arealg, "p_arealg");

   p_dxx = op_decl_dat(elems, 1, "double", dxx, "p_dxx");
   p_dyy = op_decl_dat(elems, 1, "double", dyy, "p_dyy");
   p_dzz = op_decl_dat(elems, 1, "double", dzz, "p_dzz");
   
   p_ss = op_decl_dat(elems, 1, "double", ss, "p_ss");
   p_elemMass = op_decl_dat(elems, 1, "double", elemMass, "p_elemMass");
   p_vnew = op_decl_dat(elems, 1, "double", vnew, "p_vnew");

   //Temporary
   p_sigxx = op_decl_dat(elems, 3, "double", sigxx, "p_sigxx");
   p_determ = op_decl_dat(elems, 1, "double", determ, "p_determ");

   p_dvdx = op_decl_dat(elems, 8, "double",dvdx, "dvdx");
   p_dvdy = op_decl_dat(elems, 8, "double",dvdy, "dvdy");
   p_dvdz = op_decl_dat(elems, 8, "double",dvdz, "dvdz");
   p_x8n = op_decl_dat(elems, 8, "double",x8n, "x8n");
   p_y8n = op_decl_dat(elems, 8, "double",y8n, "y8n");
   p_z8n = op_decl_dat(elems, 8, "double",z8n, "z8n");

   p_delv_xi = op_decl_dat(elems, 1, "double", delv_xi, "p_delv_xi"); 
   p_delv_eta = op_decl_dat(elems, 1, "double", delv_eta, "p_delv_eta"); 
   p_delv_zeta = op_decl_dat(elems, 1, "double", delv_zeta, "p_delv_zeta"); 

   p_delx_xi = op_decl_dat(elems, 1, "double", delx_xi, "p_delx_xi"); 
   p_delx_eta = op_decl_dat(elems, 1, "double", delx_eta, "p_delx_eta"); 
   p_delx_zeta = op_decl_dat(elems, 1, "double", delx_zeta, "p_delx_zeta"); 

   // //Declare Constants
   op_decl_const(1, "double", &m_e_cut);
   op_decl_const(1, "double", &m_p_cut);
   op_decl_const(1, "double", &m_q_cut);
   op_decl_const(1, "double", &m_v_cut);
   op_decl_const(1, "double", &m_u_cut);

   op_decl_const(1, "double", &m_hgcoef);
   op_decl_const(1, "double", &m_ss4o3);
   op_decl_const(1, "double", &m_qstop);
   op_decl_const(1, "double", &m_monoq_max_slope);
   op_decl_const(1, "double", &m_monoq_limiter_mult);
   op_decl_const(1, "double", &m_qlc_monoq);
   op_decl_const(1, "double", &m_qqc_monoq);
   op_decl_const(1, "double", &m_qqc);
   op_decl_const(1, "double", &m_eosvmax);
   op_decl_const(1, "double", &m_eosvmin);
   op_decl_const(1, "double", &m_pmin);
   op_decl_const(1, "double", &m_emin);
   op_decl_const(1, "double", &m_dvovmax);
   op_decl_const(1, "double", &m_refdens);
   op_diagnostic_output();

   // op_par_loop(test_e, "test_e", elems, 
   //             op_arg_dat(p_e, -1, OP_ID, 1, "double", OP_WRITE));
   
   // Real_t* x_t = (Real_t*) malloc(m_numNode*sizeof(Real_t));
   // op_fetch_data(p_x, x_t);
   // std::cout << std::scientific << std::setprecision(6);
   // // std::cout << std::setw(12) << locDom.e(ElemId) << "\n";
   // std::cout << std::setw(12) << x[3] << "\n";
   // std::cout << std::setw(12) << x_t[3] << "\n";


#if USE_MPI   
   fieldData = &Domain::nodalMass ;

   // Initial domain boundary communication 
   CommRecv(*locDom, MSG_COMM_SBN, 1,
            locDom->sizeX() + 1, locDom->sizeY() + 1, locDom->sizeZ() + 1,
            true, false) ;
   CommSend(*locDom, MSG_COMM_SBN, 1, &fieldData,
            locDom->sizeX() + 1, locDom->sizeY() + 1, locDom->sizeZ() +  1,
            true, false) ;
   CommSBN(*locDom, 1, &fieldData) ;

   // End initialization
   MPI_Barrier(MPI_COMM_WORLD);
#endif   
   
   // BEGIN timestep to solution */
#if USE_MPI   
   double start = MPI_Wtime();
#else
   timeval start;
   gettimeofday(&start, NULL) ;
#endif
//debug to see region sizes
//   for(Int_t i = 0; i < locDom->numReg(); i++)
//      std::cout << "region" << i + 1<< "size" << locDom->regElemSize(i) <<std::endl;

   // !Main Loop
   // while((locDom->time() < locDom->stoptime()) && (locDom->cycle() < opts.its)) {
   while((m_time < m_stoptime) && (m_cycle < opts.its)) {
   // while((time < stoptime) && (cycle < opts.its)) {
      TimeIncrement(*locDom) ;
      LagrangeLeapFrog(*locDom) ;

      if ((opts.showProg != 0) && (opts.quiet == 0) && (myRank == 0)) {
         // std::cout << "cycle = " << locDom->cycle()       << ", "
         std::cout << "cycle = " << m_cycle       << ", "
                   << std::scientific
                  //  << "time = " << double(locDom->time()) << ", "
                  //  << "dt="     << double(locDom->deltatime()) << "\n";
                   << "time = " << double(m_time) << ", "
                   << "dt="     << double(m_deltatime) << "\n";
         std::cout.unsetf(std::ios_base::floatfield);
      }
   }

   // Use reduced max elapsed time
   double elapsed_time;
#if USE_MPI   
   elapsed_time = MPI_Wtime() - start;
#else
   timeval end;
   gettimeofday(&end, NULL) ;
   elapsed_time = (double)(end.tv_sec - start.tv_sec) + ((double)(end.tv_usec - start.tv_usec))/1000000 ;
#endif
   double elapsed_timeG;
#if USE_MPI   
   MPI_Reduce(&elapsed_time, &elapsed_timeG, 1, MPI_DOUBLE,
              MPI_MAX, 0, MPI_COMM_WORLD);
#else
   elapsed_timeG = elapsed_time;
#endif

   // Write out final viz file */
   if (opts.viz) {
      DumpToVisit(*locDom, opts.numFiles, myRank, numRanks) ;
   }
   
   if ((myRank == 0) && (opts.quiet == 0)) {
      VerifyAndWriteFinalOutput(elapsed_timeG, *locDom,m_cycle, opts.nx, numRanks, e);
   }

   delete locDom;
   op_exit(); 

#if USE_MPI
   MPI_Finalize() ;
#endif

   return 0 ;
}
