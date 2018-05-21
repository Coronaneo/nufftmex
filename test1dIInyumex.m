nj=10;
ms=10;
iflag=-1
eps=1e-12
xj=rand(nj,1);
fk=randn(ms,1)+i*rand(ms,1)
[cj] = nufft1dIInyumex(xj,iflag,eps,ms,fk)
