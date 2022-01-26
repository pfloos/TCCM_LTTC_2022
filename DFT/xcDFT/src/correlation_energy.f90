function correlation_energy(rung,nGrid,weight,rho,drho) result(Ec)

! Compute the correlation energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: rung
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(3,nGrid)

! Local variables

  double precision              :: EcLDA,EcGGA
  double precision              :: aC
  double precision              :: Ec

! Output variables

! Memory allocation

  Ec    = 0d0
  EcLDA = 0d0
  EcGGA = 0d0

  select case (rung)

!   Hartree calculation
    case(0) 

      Ec = 0d0

!   LDA functionals
    case(1) 

!     call lda_correlation_energy(nGrid,weight,rho,EcLDA)

      Ec = EcLDA

!   GGA functionals
    case(2) 

!     call gga_correlation_energy(nGrid,weight,rho,drho,EcGGA)

      Ec = EcGGA

!   Hybrid functionals
    case(4) 

      aC = 0.81d0

!     call lda_correlation_energy(nGrid,weight,rho,EcLDA)
!     call gga_correlation_energy(nGrid,weight,rho,drho,EcGGA)

      Ec = EcLDA + aC*(EcGGA - EcLDA) 

!   Hartree-Fock calculation
    case(666) 

      Ec = 0d0

  end select
 
end function correlation_energy
