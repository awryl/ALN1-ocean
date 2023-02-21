!ifort -c -O -fPIC -vec-report0 m_subspace_iter.f90
mex mex_subspace_iter.F90 subspace_iter.f90 m_subspace_iter.f90  FC="ifort" FFLAGS="-O -fPIC -vec-report0" CC="" LD="ifort" FLIBS="-L/applications/matlabr2011a/bin/glnxa64 -lmx -lmex -lmat -lm -L/mnt/n7fs/ens/tp_guivarch/opt/OpenBLAS/ -lopenblas";

n=256;
m=64;
p=1;
[q r]=qr(rand(n,n));
a=q*diag(n:-1:1)*q';
v0=rand(n,m);
percentage=0.2;
maxit=5000;
eps=1e-12;

for version=0:2
  tic;
  [vout,w,res_ev,it_ev,it,n_ev]=mex_subspace_iter(version,a,p,v0,percentage,maxit,eps);
  toc;
  it
  n_ev
end


