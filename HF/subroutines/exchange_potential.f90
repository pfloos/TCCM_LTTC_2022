subroutine exchange_potential(nBas,P,ERI,K)

! Compute Coulomb matrix

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: mu,nu,la,si

! Output variables

  double precision,intent(out)  :: K(nBas,nBas)

  K = 0d0
  do mu=1,nBas
    do nu=1,nBas
      do la=1,nBas
        do si=1,nBas
          K(mu,nu) = K(mu,nu) - 0.5d0*P(la,si)*ERI(mu,la,si,nu)
        enddo
      enddo
    enddo
  enddo


end subroutine exchange_potential
