
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
 @file ComputeSPMV_ref.cpp

 HPCG routine
 */

#include "ComputeSPMV_ref.hpp"

#ifndef HPCG_NO_MPI
#include "ExchangeHalo.hpp"
#endif

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif
#include <cassert>

/*!
  Routine to compute matrix vector product y = Ax where:
  Precondition: First call exchange_externals to get off-processor values of x

  This is the reference SPMV implementation.  It CANNOT be modified for the
  purposes of this benchmark.

  @param[in]  A the known system matrix
  @param[in]  x the known vector
  @param[out] y the On exit contains the result: Ax.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeSPMV
*/
int ComputeSPMV_ref( const SparseMatrix & A, Vector & x, Vector & y) {

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
  #pragma omp parallel for
#endif
  for (local_int_t iz=0; iz<nz; iz++) {
    global_int_t giz = giz0+iz;
    for (local_int_t iy=0; iy<ny; iy++) {
      global_int_t giy = giy0+iy;
      for (local_int_t ix=0; ix<nx; ix++) {
        global_int_t gix = gix0+ix;
        local_int_t currentLocalRow = iz*nx*ny+iy*nx+ix;
        global_int_t currentGlobalRow = giz*gnx*gny+giy*gnx+gix;
#ifndef HPCG_NO_OPENMP
// C++ std::map is not threadsafe for writing
        #pragma omp critical
#endif
        A.globalToLocalMap[currentGlobalRow] = currentLocalRow;

        A.localToGlobalMap[currentLocalRow] = currentGlobalRow;
#ifdef HPCG_DETAILED_DEBUG
        HPCG_fout << " rank, globalRow, localRow = " << A.geom->rank << " " << currentGlobalRow << " " << A.globalToLocalMap[currentGlobalRow] << endl;
#endif

        double * matrixValues = new double[];

        char numberOfNonzerosInRow = 0;
        double * currentValuePointer = matrixValues[currentLocalRow]; // Pointer to current value in current row
        global_int_t * currentIndexPointerG = mtxIndG[currentLocalRow]; // Pointer to current index in current row
        for (int sz=-1; sz<=1; sz++) {
          if (giz+sz>-1 && giz+sz<gnz) {
            for (int sy=-1; sy<=1; sy++) {
              if (giy+sy>-1 && giy+sy<gny) {
                for (int sx=-1; sx<=1; sx++) {
                  if (gix+sx>-1 && gix+sx<gnx) {
                    global_int_t curcol = currentGlobalRow+sz*gnx*gny+sy*gnx+sx;
                    if (curcol==currentGlobalRow) {
                    // matrixDiagonal[currentLocalRow] = currentValuePointer;
                      *currentValuePointer++ = 26.0;
                    } else {
                      *currentValuePointer++ = -1.0;
                    }
                    *currentIndexPointerG++ = curcol;
                    numberOfNonzerosInRow++;
                  } // end x bounds test
                } // end sx loop
              } // end y bounds test
            } // end sy loop
          } // end z bounds test
        } // end sz loop
        nonzerosInRow[currentLocalRow] = numberOfNonzerosInRow;
#ifndef HPCG_NO_OPENMP
        #pragma omp critical
#endif
        localNumberOfNonzeros += numberOfNonzerosInRow; // Protect this with an atomic
        if (b!=0)      bv[currentLocalRow] = 26.0 - ((double) (numberOfNonzerosInRow-1));
        if (x!=0)      xv[currentLocalRow] = 0.0;
        if (xexact!=0) xexactv[currentLocalRow] = 1.0;

        double sum = 0.0;
        const double * const cur_vals = A.matrixValues[i];
        const local_int_t * const cur_inds = A.mtxIndL[i];
        const int cur_nnz = A.nonzerosInRow[i];

        for (int j=0; j< cur_nnz; j++)
          sum += cur_vals[j]*xv[cur_inds[j]];
        yv[i] = sum;

      } // end ix loop
    } // end iy loop
  } // end iz loop
  
  return 0;
}
