! This file is provided as part of the "projet long" for the Algebre Lineaire Numerique course
! at ENSEEIHT
! Date: 04/04/2012
! Authors: P. Amestoy, P. Berger, A. Buttari, Y. Diouane, S. Gratton, F.H. Rouet
!
! This file contains the driver for testing the developed methods
! 
program main
  use m_subspace_iter
  implicit none
  integer                       :: ierr, n, p, maxit, lwork, i, it, n_ev, disp, m, v
  double precision, allocatable :: a(:,:), v0(:,:), work(:), e_vect(:,:)
  double precision, allocatable :: w(:), res_ev(:), e_vals(:) 
  integer, allocatable          :: it_ev(:)
  integer                       :: t_start, t_end, t_rate, len
  double precision              :: time, eps, res=0d0, percentage

  !! parameter initializations 
  ierr=0

  maxit=1e4
  eps=1d-13
  percentage=10d-3

  !! specify the dimensions n, p and the number of the dominante eigen vectors needed
  write(*,'("Version?")')
  read(*,*) v
  write(*,'("p (block-product size)?")')
  read(*,*) p
  write(*,'("n (dimension of the matrix a)?")')
  read(*,*) n
  write(*,'("m ( dimension of the invariant subspace )? ")')
  read(*,*) m
  write(*,'("Percentage (the percentage of the trace needed ) ")')
  read(*,*) percentage
  write(*,'("Do you want to display the eigenvalues?(1 :yes or 0 : no) ")')
  read(*,*) disp

  lwork= m*m + 5*m + n*n
  allocate(w(m), res_ev(m), it_ev(m), &
       & a(n,n), work(lwork), v0(n,m), e_vect(n,m), e_vals(n), stat=ierr)
  if(ierr .ne.0) then
     write(*,'("Can not allocate data")')
     stop
  endif


  !! initialize the matrix A as semi-positive symmetric.
  !! Its eigenvalues are n, n-1, ..., 1 
  a = 0.d0
  do i=1,n
     e_vals(i) = real(n-i+1)
     a(i,i)    = e_vals(i)
  end do

  e_vect = 0.d0
  do i=1, m
     e_vect(i,i)=1.d0
  end do

  call random_number(v0)
  call dgeqrf(n, m, v0, n, w, work, lwork, ierr)
  call dormqr('l', 'n', n, n, m, v0, n, w, a, n, work, lwork, ierr)
  call dormqr('r', 't', n, n, m, v0, n, w, a, n, work, lwork, ierr)
  call dormqr('l', 'n', n, m, m, v0, n, w, e_vect, n, work, lwork, ierr)

  !! initialize the subspace v0
  call random_number(v0)
  w = 0.d0

  call system_clock(t_start)
  call subspace_iter(v, a, p, v0, w, n, m, percentage, maxit, eps, work, lwork, res_ev, it_ev, it, n_ev, ierr)
  call system_clock(t_end, t_rate)
  time = real(t_end-t_start)/real(t_rate)
  write(*,'(" ")')
  write(*,'(" ")')
  write(*,'("Done. ")')

  if(ierr .gt. 0)then
     write(*,'("Error detected in subspace iteration ",i4)')ierr
  else
     write(*,'("Elapsed time                           : ",f10.2)')time
     write(*,'("Number of computed eigenpairs          : ",i5)')n_ev
     write(*,'("Forward error on computed eigenvalues  : ",es10.2)')maxval(abs(w(1:n_ev))-abs(e_vals(1:n_ev))) / &
          & maxval(abs(e_vals(1:n_ev)))
     write(*,'("Forward error on computed eigenvectors : ",es10.2)')maxval(abs(v0(:,1:n_ev))-abs(e_vect(:,1:n_ev))) / &
          & maxval(abs(e_vect(:,1:n_ev)))
     if( disp.eq.1) then
        write(*,'("Residuals       : ")',advance='no')
        do i=1,n_ev
           write(*,'(es8.1,x)',advance='no')res_ev(i)
        end do
        write(*,'(" ")')
     end if
     if( disp.eq.1) then
        write(*,'("Eigenvalues     : ")',advance='no')
        do i=1, n_ev
           write(*,'(f8.2,x)',advance='no')w(i)
        end do
        write(*,'(" ")')
     end if
     write(*,'("# of iterations : ",i6)')it
     if( disp.eq.1) then
        write(*,'("Iter per eig    : ")',advance='no')
        do i=1, n_ev
           write(*,'(i8,x)',advance='no')it_ev(i)
        end do
        write(*,'(" ")')
     end if
  end if

  stop

end program main
