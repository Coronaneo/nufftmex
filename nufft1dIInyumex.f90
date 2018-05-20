#include "fintrf.h"

! standard header for mex functions
subroutine mexfunction(nlhs, plhs, nrhs, prhs)

  implicit double precision (a-h,o-z)

  ! number of input arguments, number of output arguments
  integer      :: nlhs, nrhs  
  ! pointer to inputs and outputs
  mwPointer    :: plhs(*), prhs(*) 
  ! get some of the matlab mex functions
  mwPointer    :: mxGetPr, mxGetPi, mxCreateDoubleMatrix 
  ! define a size integer so that we can get its type
  mwSize       :: nj

  ! input variable
  real*8 eps
  integer iflag,ms
  real*8, allocatable ::  xj(:)
  complex*16, allocatable ::  fk(:)

  ! output variable
  complex*16, allocatable ::  cj(:)

  ! aux variable
  integer ier

  nj     = mxGetM(prhs(1))


  ! copy the right-hand side argument in matlab to dnu and t
  allocate(cj(nj),xj(nj),fk(ms))

  call mxCopyPtrToReal8(mxGetPr(prhs(1)), xj, nj)
  call mxCopyPtrToInteger4(mxGetPr(prhs(2)),iflag,1)
  call mxCopyPtrToReal8(mxGetPr(prhs(3)), eps, 1)
  call mxCopyPtrToInteger4(mxGetPr(prhs(4)),ms,1)
  call mxCopyPtrToComplex16(mxGetPr(prhs(5)),mxGetPi(prhs(5)), fk, ms)

  ! initialize the left-hand side of the function
  plhs(1) = mxCreateDoubleMatrix(nj, 1, 0)


  call nufft1d1f90(nj,xj,cj,iflag,eps,ms,fk,ier)


  ! copy the output to the left-hand side
  call mxCopyComplex16ToPtr(cj, mxGetPr(plhs(1)),mxGetPi(plhs(1)),nj) 

  deallocate(cj,xj,fk)

end subroutine mexfunction