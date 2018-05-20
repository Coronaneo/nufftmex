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
  mwSize       :: m, n

  ! input variable
  real*8 eps
  real*8, allocatable ::  xj(:,:)
  complex*16, allocatable ::  cj(:,:)
  complex*16, allocatable ::  v(:)

  ! output variable
  complex*16, allocatable ::  fk(:,:)

  ! aux variable
  integer ier
  integer j

  m     = mxGetM(prhs(1))
  n     = mxGetN(prhs(1))

  ! copy the right-hand side argument in matlab to dnu and t
  allocate(cj(m,n),xj(m,n),fk(m,n),v(m))

  call mxCopyPtrToReal8(mxGetPr(prhs(1)), xj, m*n)
  call mxCopyPtrToComplex16(mxGetPr(prhs(2)),mxGetPi(prhs(2)), cj, m*n)
  call mxCopyPtrToComplex16(mxGetPr(prhs(3)),mxGetPi(prhs(3)), v, m)
  call mxCopyPtrToReal8(mxGetPr(prhs(4)), eps, 1)

  ! initialize the left-hand side of the function
  plhs(1) = mxCreateDoubleMatrix(m, n, 0)

  do j=1,n
     call nufft1d1f90(m,xj(:,j),cj(:,j)*v,1,eps,m,fk(:,j),ier)
  enddo

  ! copy the output to the left-hand side
  call mxCopyComplex16ToPtr(fk, mxGetPr(plhs(1)),mxGetPi(plhs(1)),m*n) 

  deallocate(cj,xj,fk,v)

end subroutine mexfunction
