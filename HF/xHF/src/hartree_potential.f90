subroutine hartree_potential(nBas,P,ERI,J)

! Compute Coulomb matrix

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: mu,nu,la,si

! Output variables

  double precision,intent(out)  :: J(nBas,nBas)

  J = 0d0


end subroutine hartree_potential
