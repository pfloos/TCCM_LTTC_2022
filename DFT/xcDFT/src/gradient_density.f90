subroutine gradient_density(nGrid,nBas,P,AO,dAO,drho)

! Calculate gradient of the one-electron density

  implicit none
  include 'parameters.h'

! Input variables

  double precision,parameter    :: thresh = 1d-15

  integer,intent(in)            :: nGrid
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: dAO(3,nBas,nGrid)

! Local variables

  integer                       :: ixyz,iG,mu,nu
  double precision,external     :: trace_matrix

! Output variables

  double precision,intent(out)  :: drho(3,nGrid)

  drho(:,:) = 0d0

end subroutine gradient_density
