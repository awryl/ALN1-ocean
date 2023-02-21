#include "fintrf.h"     
subroutine mexfunction(nlhs,plhs,nrhs,prhs)
! matlab interface:
!    [vout, w, res_ev, it_ev, it, n_ev]
!  = mex_subspace_iter(version, a, p, v, percentage, maxit, eps)
!
! input:
! version      0,1,2 depending on the version to be called
! a            input matrix
! p            number of products to perform before reorthogonalizing
! v            initial guess for the subspace
! percentage   fraction of the trace to the recover
! maxit        maximum number of iterations
! eps          tolerance for the stopping criterion
!
! ouput:
! vout         computed eigenvectors
! w            computed eigenvalues
! res_ev       residual for each eigenvector (except for version 0)
! it_ev        number of iterations for each eigenvector
! it           total number of iterations
! n_ev         number of converged eigenpairs

  implicit none
 
  mwpointer :: mxgetpr, mxcreatedoublematrix 
  mwsize    :: mxgetn, mxgetm 

  integer   :: nlhs
  integer   :: nrhs
  mwpointer :: plhs(*)
  mwpointer :: prhs(*)

  integer   :: version, p, maxit, it, n_ev, lwork, ierr
  integer, allocatable, dimension(:) :: it_ev
  double precision, allocatable, dimension(:):: work
  mwsize    :: n, m
  mwpointer :: a, v, vout, w, res_ev, eps, percentage
  mwpointer :: d_version, d_p, d_it, d_n_ev, d_it_ev, d_maxit

  if(nrhs .ne. 7) then
    call mexerrmsgidandtxt ('matlab:mysubmex:ninput','6 inputs required.')
  endif
  if(nlhs .gt. 6) then
    call mexerrmsgidandtxt ('matlab:mysubmex:noutput','too many output arguments.')
  endif     

  d_version   = mxgetpr(prhs(1))
  a           = mxgetpr(prhs(2))       
  n           = mxgetm(prhs(2))
  d_p         = mxgetpr(prhs(3)) 
  v          = mxgetpr(prhs(4)) 
  m           = mxgetn(prhs(4)) 
  percentage  = mxgetpr(prhs(5))  
  d_maxit     = mxgetpr(prhs(6)) 
  eps         = mxgetpr(prhs(7))

  call convert_to_int(%val(d_version), version)
  call convert_to_int(%val(d_p), p)
  call convert_to_int(%val(d_maxit), maxit)

  plhs(1)   = mxcreatedoublematrix(n, m, 0)
  vout      = mxgetpr(plhs(1))
  plhs(2)   = mxcreatedoublematrix(m, 1, 0)
  w         = mxgetpr(plhs(2))
  plhs(3)   = mxcreatedoublematrix(m, 1, 0)
  res_ev    = mxgetpr(plhs(3))
  plhs(4)   = mxcreatedoublematrix(m, 1, 0)
  d_it_ev   = mxgetpr(plhs(4))
  plhs(5)   = mxcreatedoublematrix(1, 1, 0)
  d_it      = mxgetpr(plhs(5))
  plhs(6)   = mxcreatedoublematrix(1, 1, 0)
  d_n_ev    = mxgetpr(plhs(6))

  lwork = m*m +5*m + n*n
  allocate(work(lwork))
  allocate(it_ev(m))

  call subspace_iter(version, %val(a), p, %val(v), %val(w), n, m,    &
                     %val(percentage), maxit, %val(eps), work, lwork, &
                     %val(res_ev), it_ev, it, n_ev, ierr)

  if(ierr.ne. 0) then
    write(*,*)'Error in subspace_iter; ierr=',ierr
  end if

  call convert_to_double(it, %val(d_it), 1)
  call convert_to_double(n_ev, %val(d_n_ev), 1)
  call convert_to_double(it_ev, %val(d_it_ev), m)
  call copy_array(%val(v), %val(vout), n, m)

  deallocate(it_ev)
  deallocate(work)

end subroutine mexfunction

subroutine convert_to_int(d_i, i)

  integer :: i
  double precision :: d_i

  i = int(d_i)

end subroutine convert_to_int

subroutine convert_to_double(i, d_i, n)
  
  integer :: i(n)
  double precision :: d_i(n)
  
  d_i = i

end subroutine convert_to_double

subroutine copy_array(v, w, n, m)

  integer :: n, m
  double precision :: v(n,m), w(n,m)

  w = v

end subroutine copy_array

