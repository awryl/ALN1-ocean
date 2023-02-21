! This file is provided as part of the "projet long" for the Algebre Lineaire Numerique course
! at ENSEEIHT
! Date: 04/04/2012
! Authors: P. Amestoy, P. Berger, A. Buttari, Y. Diouane, S. Gratton, F.H. Rouet
!
! This file contains a subroutine that wraps the three routines related to three
! different version of the subspace iteration method.
! Note that not all of the arguments are used for the call to the first, simpler
! routine.
! 
subroutine subspace_iter(version, a, p, v, w, n, m, percentage, maxit, eps, work,    &
     lwork, res_ev, it_ev, it, n_ev, ierr)
  use m_subspace_iter
  implicit none

  !! the version to be called : 0, 1 or 2
  integer,          intent(in)                            :: version
  !! the subspace dimensions
  integer,          intent(in)                            :: n, m
  !! the number of products to perform before reorthogonalizing
  integer,          intent(in)                            :: p
  !! the target matrix 
  double precision, dimension(n,n),intent(in)             :: a
  !! maximum # of iteration
  integer,          intent(in)                            :: maxit
  !! the number of dominant eigenvectors to compute
  double precision, intent(in)                            :: percentage
  !! the tolerance for the stopping criterion
  double precision, intent(in)                            :: eps
  !! length of the workspace
  integer                                                 :: lwork
  !! the workspace
  double precision, dimension(lwork)                      :: work
  !! the starting subspace. The computed eigenvectors will be
  !! returned in this array
  double precision, dimension(n,m),intent(inout)          :: v
  !! the residuals for each eigenvector        
  double precision, dimension(m),intent(out)              :: res_ev
  !! the n_ev dominant eigenvalues                        
  double precision, dimension(m),intent(out)              :: w
  !! the number of iteration  to converge for each eigenvector
  integer,          dimension(m),intent(out)              :: it_ev
  !! the global number iteration to converge
  integer,          intent(out)                           :: it
  !! the nummber of converged eigenvectors
  integer,          intent(out)                           :: n_ev
  !! a flag for signaling errors
  integer,intent(out)                                     :: ierr


  select case(version)
  case(0)
     ! Basic subspace iteration. Note that in this case not all the arguments provided
     ! in input to this routine are passed. Moreover, for this basic version one global
     ! residual exists (and not one residual for each eigenvector). 
     call subspace_iter_v0(a, v, w, n, m, maxit, eps, work, lwork, res_ev(1), it, ierr)
     ! Fill the unused argument with meaningful stuff
     res_ev(2:m) = -1
     it_ev       = it
     n_ev        = m
  case(1)
     ! V1 is the subspace iteration method with Rayleigh-Ritz.
     call subspace_iter_v1(a, v, w, n, m, percentage, maxit, eps, work, lwork, res_ev, it_ev, it, n_ev, ierr)
  case(2)
     ! V2 is the subspace iteration method with Rayleigh-Ritz plus acceleration methods.
     ! V1 can be accelerated in two ways
     ! 1) by performing a block of p products y=A*v before reorthogonlizing
     ! 2) by defating each time some eigenvalue has converged
     call subspace_iter_v2(a, p, v, w, n, m, percentage, maxit, eps, work, lwork, res_ev, it_ev, it, n_ev, ierr)
  case default
     write(*,'("Version is not valid")')

  end select


  return

end subroutine subspace_iter
