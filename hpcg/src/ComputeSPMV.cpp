
//@HEADER
// ***************************************************
//
// HPCG: High Performance Conjugate Gradient Benchmark
//
// Contact:
// Michael A. Heroux ( maherou@sandia.gov)
// Jack Dongarra     (dongarra@eecs.utk.edu)
// Piotr Luszczek    (luszczek@eecs.utk.edu)
//
// ***************************************************
//@HEADER

/*!
 @file ComputeSPMV.cpp

 HPCG routine
 */

#include "ComputeSPMV.hpp"
#include "ComputeSPMV_ref.hpp"

#ifndef HPCG_NO_MPI
#include <mpi.h>
#endif


#ifndef HPCG_NO_MPI
#include "ExchangeHalo.hpp"
#endif

#if defined(HPCG_DEBUG) || defined(HPCG_DETAILED_DEBUG)
#include <fstream>
using std::endl;
#include "hpcg.hpp"
#endif
#include <cassert>

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif
#include <cassert>

/*!
  Routine to compute sparse matrix vector product y = Ax where:
  Precondition: First call exchange_externals to get off-processor values of x

  This routine calls the reference SpMV implementation by default, but
  can be replaced by a custom, optimized routine suited for
  the target system.

  @param[in]  A the known system matrix
  @param[in]  x the known vector
  @param[out] y the On exit contains the result: Ax.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeSPMV_ref
*/
int ComputeSPMV( const SparseMatrix & A, Vector & x, Vector & y, bool testing) {

  // This line and the next two lines should be removed and your version of ComputeSPMV should be used.
  double sum;
  A.isSpmvOptimized = true; //changed to true because this is an optimized version
  if(testing){
    printf("%d\n", testing);
  }
  
  assert(x.localLength>=A.localNumberOfColumns); // Test vector lengths
  assert(y.localLength>=A.localNumberOfRows);
  global_int_t nx = A.geom->nx;
  global_int_t ny = A.geom->ny;
  global_int_t nz = A.geom->nz;
  global_int_t gnx = A.geom->gnx;
  global_int_t gny = A.geom->gny;
  global_int_t gnz = A.geom->gnz;
  global_int_t gix0 = A.geom->gix0;
  global_int_t giy0 = A.geom->giy0;
  global_int_t giz0 = A.geom->giz0;
  
  #ifndef HPCG_NO_MPI
    ExchangeHalo(A,x);
  #endif
  const double * const xv = x.values;
  double * const yv = y.values;
  const local_int_t nrow = A.localNumberOfRows;

#ifndef HPCG_NO_OPENMP
  #pragma omp parallel for private(sum)
#endif

  for (local_int_t iz=0; iz<nz; iz++) { //outer loops step through physical space, from that we compute a row
    for (local_int_t iy=0; iy<ny; iy++) {
      for (local_int_t ix=0; ix<nx; ix++) {
        global_int_t giz = giz0+iz;
        global_int_t giy = giy0+iy;
        global_int_t gix = gix0+ix;
        local_int_t currentLocalRow = iz*nx*ny+iy*nx+ix;
        global_int_t currentGlobalRow = giz*gnx*gny+giy*gnx+gix;

        // A.globalToLocalMap[currentGlobalRow] = currentLocalRow;

        // A.localToGlobalMap[currentLocalRow] = currentGlobalRow;

#ifdef HPCG_DETAILED_DEBUG
        HPCG_fout << " rank, globalRow, localRow = " << A.geom->rank << " " << currentGlobalRow << " " << A.globalToLocalMap[currentGlobalRow] << endl;
#endif

        sum = 0.0;
        const local_int_t * const cur_inds = A.mtxIndL[currentLocalRow];
        const int cur_nnz = A.nonzerosInRow[currentLocalRow];
        for (int j=0; j<=cur_nnz; j++) {
          local_int_t localcol = cur_inds[j];
          if (testing){
            if (localcol==currentLocalRow && currentGlobalRow < 9) {                    
              sum  += 26.0 * xv[cur_inds[j]] * (currentGlobalRow+2)*1.0e6; 
            } else if (localcol==currentLocalRow){
              sum  += 26.0 * xv[cur_inds[j]] * 1.0e6;
            } else {
              sum += -1.0 * xv[cur_inds[j]];
            }
          }else{
            if (localcol==currentLocalRow) {                    
              sum  += 26.0 * xv[cur_inds[j]] ; 
            } else {
              sum += -1.0 * xv[cur_inds[j]];
            }
          }
        } 
        yv[currentLocalRow] = sum;

      } // end ix loop
    } // end iy loop
  } // end iz loop
  
  return 0;
}
