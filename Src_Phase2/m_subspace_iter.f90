! This file is provided as part of the "projet long" for the Algebre Lineaire Numerique course
! at ENSEEIHT
! Date: 04/04/2012
! Authors: P. Amestoy, P. Berger, A. Buttari, Y. Diouane, S. Gratton, F.H. Rouet
!
! This file holds a FORTRAN module that includes the implementation 
! of the three variants of the subspace
! iteration method.
! All routines are empty, only the interface with related 
! declaration is provided
! 
module m_subspace_iter
  implicit none
  private
  public :: subspace_iter_v0, subspace_iter_v1, subspace_iter_v2
contains


! Renvoie le produit scalaire d deux vecteurs de taille dim
function prod_scal(v1, v2, dim)
	implicit none
	integer, intent(in)								:: dim
	double precision, dimension(dim), intent(in)	:: v1, v2
	double precision								:: prod_scal
	integer											:: i
	! Initialisation
	prod_scal = 0.0
	! calcul de v1.v2
	do i=1,dim
		prod_scal = prod_scal + v1(i)*v2(i)
	end do
	return
end function prod_scal


! Renvoie la norme 2 d'un vecteur de R^n
function norme2 (a, n)
	integer, intent(in)							:: n
	double precision, dimension(n), intent(in)	:: a
	double precision							:: norme2
	integer										:: i
	! Initialisation
	norme2 = 0
	! Calcul de la norme2
	do i=1,n
		norme2 = norme2 + (a(i))**2		
	end do
	norme2 = sqrt(norme2)
	return
end function norme2


! Renvoie la norme 2 d'une matrice de taille (m,n)
function norme2mat (a, m, n)
	integer, intent(in)								:: m, n
	double precision, dimension(m,n), intent(in)	:: a
	double precision								:: norme2mat
	integer											:: i, j
	! Initialisation
	norme2mat = 0
	! Calcul de la norme2 de la matrice
	do i=1,n
		do j=1,n
			norme2mat = norme2mat + (a(i,j))**2	
		end do	
	end do
	norme2mat = sqrt(norme2mat)
	return
end function norme2mat


! Orthonormalisation de V de taille (n,m)
subroutine ortho(v, n, m)
	implicit none
	integer, intent(in)								:: n, m
	double precision, dimension(n,m), intent(inout)	:: v
	double precision						     	:: p
	integer											:: i, j, k

	do j=1,m
		! On projette Vj sur Vi
		do i=1,(j-1)
			p = prod_scal(v(:,i),v(:,j),n)
			do k=1,n 
				v(k,j)=v(k,j) - p*v(k,i)
			end do
	   	end do
		! On rend unitaire les vecteurs
		p = norme2(v(:,j),n)
		do k=1,n
			v(k,j)=v(k,j)/p
		end do
	end do
	return
end subroutine ortho

! Retourne un vecteur
subroutine reverse(w,n)
	implicit none
	integer, intent(in)								:: n
	double precision, dimension(n), intent(inout)	:: w
	double precision, dimension(n)					:: wbis
	integer											:: i
	
	do i=1,n
		wbis(i) = w(n-i+1)
	end do
	w = wbis
	return
end subroutine reverse
	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Simple, basic version
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine subspace_iter_v0(a, v, w, n , m, maxit, eps, work,  &
lwork, res, it, ierr)
    implicit none
    !! the subspace dimensions
    integer,          intent(in)						:: n, m
    !! the target matrix                              
    double precision, dimension(n, n), intent(in)		:: a
    !! maximum # of iteration                         
    integer,          intent(in)						:: maxit
    !! the tolerance for the stopping criterion       
    double precision, intent(in)						:: eps
    !! the length of the workspace                    
    integer												:: lwork
    !! the workspace                                  
    double precision, dimension(lwork)					:: work
    !! the starting subspace. The computed eigenvectors will be
    !! returned in this array
    double precision, dimension(n, m), intent(inout)	:: v
    !! the m dominant eigenvalues
    double precision, dimension(m), intent(out)			:: w
    !! the returned residual                          
    double precision, intent(out)						:: res
    !! the number of iteration to converge            
    integer,          intent(out)						:: it
    !! a flag for signaling errors                    
    integer,          intent(out)						:: ierr
                 	                                     
	double precision, dimension(n,m)					:: y
    double precision, dimension(m,m)			     	:: h
    double precision					     			:: p
    integer							     				:: i, j, k

    ierr = 0

	! Si m>=n, il y a une erreur
    if(m.gt.n)then
       ierr=1
       return
    end if
	
	! On initialise Y avec m vecteurs linéairement indép.
	call random_number(y)
	
	res = eps+1
	it=0

 	do while ((res >= eps).and.(it <= maxit)) 
		! On orthonormalise Y
		call ortho(y, n, m)
		! V:=Y
		v=y
		! Y := AV
		y=matmul(a, v)
		! On forme le quotient de Rayleigh
		h=matmul(transpose(v), matmul(a, v))
		! On incrémente le nombre d'itérations
		it=it+1
		! On met à jour le résidu
		res=(norme2mat(matmul(a,v)-matmul(v,h),n,m)/norme2mat(a,n,n)) 
	end do 
	
	! On stocke les valeurs propres dans w
	do i=1,m
		w(i)=h(i,i)
	end do
	
    return
  end subroutine subspace_iter_v0

   


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! V1 version of the subspace iteration with Rayleigh-Ritz.
! In this case the convergence can be checked for each eigenvector separately
! and the method can be stopped when the convereged eigenvectors capture
! as much as the trace of A as requested by the "percentage" argument
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine subspace_iter_v1(a, v, w, n, m, percentage, maxit, eps, work,    &
       lwork, res_ev, it_ev, it, n_ev, ierr)
    implicit none
    !! the subspace dimensions
    integer,          intent(in)							:: n, m
    !! the traget matrix 
    double precision, dimension(n,n), intent(in)			:: a
    !! maximum # of iteration
    integer,          intent(in)							:: maxit
    !! the number of the dominant eigenvectors to compute
    double precision, intent(in)							:: percentage
    !! the tolerance for the stopping criterion
    double precision, intent(in)							:: eps
    !! length of the workspace
    integer													:: lwork
    !! the workspace
    double precision, dimension(lwork)						:: work
    !! the starting subspace. The computed eigenvectors will be
    !! returned in this array
    double precision, dimension(n,m), intent(inout)			:: v
    !! the residuals for each each eigenvector        
    double precision, dimension(m),   intent(out)			:: res_ev
    !! the n_ev dominant eigenvalues                        
    double precision, dimension(m),   intent(out)			:: w
    !! the number of iteration to converge for each eigenvector
    integer,          dimension(m),   intent(out)			:: it_ev
    !! the global number iteration to converge
    integer,          intent(out)							:: it
    !! the number of converged eigenvectors
    integer,          intent(out)							:: n_ev
    !! a flag for signaling errors
    integer,          intent(out)							:: ierr

    double precision										:: percentReached, aux
	double precision, dimension(m,m)						:: h
    double precision										:: trace 
	double precision										:: p
	integer													:: i, info
	
    ierr = 0
	
    ! Si on n'a pas 0<=percentage<=1, il y a une erreur
    if((percentage.gt.1d0)  .or. (percentage.lt.0d0)) then
       ierr=1
       return
    end if
	
	! Si m>=n, il y a une erreur
    if(m.gt.n) then  
       ierr=1
       return
    end if
	
	! On initialise V avec m vecteurs linéairement indép.
	call random_number(v)

	! On calcule la trace de A
    trace = 0
    do i=1,n
		trace = trace + a(i,i)
    end do
		
	! Initialisation des variables
	it = 0
	percentReached =0
	n_ev = 0

	do while ((percentReached <= percentage).and.(it<=maxit).and.(n_ev<m))
		! V:=AV
		v = matmul(a, v)
		
		! On orthonormalise V
		call ortho(v,n,m)
		
		! On calcule la projection de Rayleigh-Ritz 
		h=matmul(transpose(v),matmul(a,v))
		call dsyev('V','U', m, h, m, w,work,lwork,info)
		
		! V:=VH
		v=matmul(v,h)
		
		aux = norme2(matmul(a,v(:,m-n_ev))-w(m-n_ev)*v(:,m-n_ev),n)
		
		! Sauvegarde des vecteurs propres
		do while (aux < trace*eps) 
			it_ev(n_ev+1) = it                	
			res_ev(n_ev+1) = aux
			percentReached = percentReached + w(m-n_ev)/trace 			
			n_ev = n_ev + 1 
			if (n_ev.gt.m) then
			   aux = trace*eps+1 
			else
				aux = norme2(matmul(a,v(:,m-n_ev))-w(m-n_ev)*v(:,m-n_ev),n)
			end if	
		end do 
		
		! On incrémente le nombre d'itérations
		it = it + 1 
	end do
	
	call reverse(w,m)
	
	do i=1,n
		call reverse(v(i,:),m)
	end do
	
    return

  end subroutine subspace_iter_v1




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! V2 is the subspace iteration method with Rayleigh-Ritz plus acceleration methods.
! V1 can be accelerated in two ways
! 1) by performing a block of p products y=A*v before reorthogonlizing
! 2) by defating each time some eigenvalue has converged
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine subspace_iter_v2(a, p, v, w, n, m, percentage, maxit, eps, work,    &
       lwork, res_ev, it_ev, it, n_ev, ierr)
    implicit none
    !! the subspace dimensions
    integer,          intent(in)							:: n, m
    !! the iteration number ( the power)
    integer,          intent(in)							:: p
    !! the regarded matrix 
    double precision, dimension(n,n),intent(in)				:: a
    !! maximum of iteration
    integer,          intent(in)							:: maxit
    !! the number of the dominant eigen vector to compute
    double precision, intent(in)							:: percentage
    !! a tolerance
    double precision, intent(in)							:: eps
    !! length of the spacework
    integer													:: lwork
    !! the spacework
    double precision, dimension(lwork)						:: work
    !! the starting subspace
    double precision, dimension(n,m),intent(inout)			:: v
    !! the starting subspace. The computed eigenvectors will be
    !! returned in this array
    double precision, dimension(m),intent(out)				:: res_ev
    !! the n_ev dominant eigenvalues                        
    double precision, dimension(m),intent(out)				:: w
    !! the number of iteration  to converge for each eigenvector
    integer,          dimension(m),intent(out)				:: it_ev
    !! the number iteration to converge
    integer,          intent(out)							:: it
    !! the nummber of the converged eigenvectors
    integer,          intent(out)							:: n_ev
    !! a flag
    integer,intent(out)										:: ierr

	double precision										:: percentReached, aux
	double precision, dimension(m,m)						:: h
    double precision										:: trace
	integer													:: i, j, k, info
	
    ierr = 0
    ! Si on n'a pas 0<=percentage<=1, il y a une erreur
    if(percentage .gt.1d0  .or.percentage .lt.0d0 )then
       ierr=1
       return
    end if
	
	! Si m>=n, il y a une erreur
    if(m.gt.n)then
       ierr=1
       return
    end if
	
	! Si p<=0, il y a une erreur
    if( p .le.0 )then
       ierr=1
       return
    end if

	! On initialise V avec m vecteurs linéairement indép.
    call random_number(v)

	! On calcule la trace de A
    trace = 0
    do i=1,n
		trace = trace + a(i,i)
    end do
	
	! Initialisation des variables
	it = 0 
	percentReached = 0
	n_ev = 0

	do while ((percentReached <= percentage).and.(it<=maxit).and.(n_ev<m))
	
		! Multiplication par blocs et déflation
		do k =1,p		
			do j = m-n_ev, m
				v(:,j) = matmul(a,v(:,j))
			end do
		end do
		
		! V:=AV
		v = matmul(a, v)
		
		! On orthonormalise V
		call ortho(v,n,m)
		
		! On calcule la projection de Rayleigh-Ritz 
		h=matmul(transpose(v),matmul(a,v))
		call dsyev('V','U', m, h, m, w,work,lwork,info)
		
		! V:=VH
		v=matmul(v,h)
		
		aux = norme2(matmul(a,v(:,m-n_ev))-w(m-n_ev)*v(:,m-n_ev),n)
		
		! Sauvegarde des vecteurs propres
		do while (aux < trace*eps) 
			it_ev(n_ev+1) = it                	
			res_ev(n_ev+1) = aux
			percentReached = percentReached + w(m-n_ev)/trace 			
			n_ev = n_ev + 1 
			if (n_ev.gt.m) then
			   aux = trace*eps+1 
			else
				aux = norme2(matmul(a,v(:,m-n_ev))-w(m-n_ev)*v(:,m-n_ev),n)
			end if	
		end do 
		
		! On incrémente le nombre d'itérations
		it = it + 1 
	end do
	
	call reverse(w,m)
	
	do i=1,n
		call reverse(v(i,:),m)
	end do

    return
  end subroutine subspace_iter_v2


  
end module m_subspace_iter



