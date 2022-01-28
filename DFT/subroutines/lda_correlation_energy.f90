subroutine lda_correlation_energy(nGrid,weight,rho,Ec)

! Compute LDA correlation energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Local variables

  integer                       :: iG
  double precision              :: a,d
  double precision              :: r,ra,rb,epsc
  double precision,parameter    :: threshold = 1d-15

! Output variables

  double precision              :: Ec


! Coefficients for Wigner's LDA correlation
  
  a = 0.04918d0
  d = 0.349d0

! Compute LDA correlation energy
  
  Ec = 0d0
  
  do iG=1,nGrid
  
    ra = max(0d0,rho(iG)/2d0)
    rb = max(0d0,rho(iG)/2d0)

    if(ra > threshold .and. rb > threshold) then

      r = ra + rb

      epsc = ra*rb/(r + d*r**(2d0/3d0))

      Ec = Ec + weight(iG)*epsc

    end if

  end do

  Ec = -4d0*a*Ec

end subroutine lda_correlation_energy
