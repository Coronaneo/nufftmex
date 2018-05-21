format long
nj=10;
ms=10;
iflag=-1
eps=1e-12
xj=(pi*[0:nj-1]/nj)'
fk=[0:ms-1]'
[cj] = nufft1dIInyumex(xj,eps,fk)
nufftc = nufftII(nj*xj/2/pi,iflag,ms,60,eps);
fftc = nufftc(fk)
