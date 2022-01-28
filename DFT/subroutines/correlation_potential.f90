subroutine correlation_potential(rung,nGrid,weight,nBas,AO,dAO,rho,drho,Fc)

! Compute the exchange potential 

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: rung
  integer,intent(in)            :: nGrid
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: dAO(3,nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(3,nGrid)

! Local variables

  double precision,allocatable  :: FcLDA(:,:),FcGGA(:,:)
  double precision              :: aC

! Output variables

  double precision,intent(out)  :: Fc(nBas,nBas)

! Memory allocation

  allocate(FcLDA(nBas,nBas),FcGGA(nBas,nBas))

  FcLDA(:,:) = 0d0
  FcGGA(:,:) = 0d0

  select case (rung)

!   Hartree calculation
    case(0) 

      Fc(:,:) = 0d0

!   LDA functionals
    case(1) 

      call lda_correlation_potential(nGrid,weight,nBas,AO,rho,FcLDA)

      Fc(:,:) = FcLDA(:,:)

!   GGA functionals
    case(2) 

!     call gga_correlation_potential(nGrid,weight,nBas,AO,dAO,rho,drho,FcGGA)
      Fc(:,:) = 0d0

      Fc(:,:) = FcGGA(:,:)

!   Hybrid functionals
    case(4) 

      aC = 0.81d0

      call lda_correlation_potential(nGrid,weight,nBas,AO,rho,FcLDA)
!     call gga_correlation_potential(nGrid,weight,nBas,AO,dAO,rho,drho,FcGGA)

      Fc(:,:) = FcLDA(:,:) + aC*(FcGGA(:,:) - FcLDA(:,:))  

!   Hartree-Fock calculation
    case(666) 

      Fc(:,:) = 0d0

  end select

end subroutine correlation_potential
