subroutine lda_exchange_energy(nGrid,weight,rho,Ex)

! Compute LDA exchange energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Local variables

  integer                       :: iG
  double precision              :: alpha
  double precision              :: r,ra,rb
  double precision,parameter    :: threshold = 1d-15

! Output variables

  double precision              :: Ex

! Cx coefficient for Slater LDA exchange

  alpha = - (1d0/2d0)**(1d0/3d0)*(3d0/2d0)*(3d0/(4d0*pi))**(1d0/3d0)

! Compute LDA exchange energy

  Ex = 0d0

  do iG=1,nGrid

    r = max(0d0,rho(iG))

    if(r > threshold) then

      Ex = Ex  + weight(iG)*alpha*r**(4d0/3d0)

    end if

  end do

end subroutine lda_exchange_energy
