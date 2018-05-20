N = 10;
A = randn(N,N); f = randn(N,N);
v=randn(N);
y = nufft1d1multi_mex(A,f,v,1e-5);



