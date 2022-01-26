function electron_number(nGrid,w,rho) result(nEl)

! Compute the number of electrons via integration of the one-electron density

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: w(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Output variables

  double precision              :: nEl

  nEl = 0d0

end function electron_number
