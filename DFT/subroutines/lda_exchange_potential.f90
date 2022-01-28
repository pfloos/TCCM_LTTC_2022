subroutine lda_exchange_potential(nGrid,weight,nBas,AO,rho,Fx)

! Compute LDA exchange potential

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Local variables

  integer                       :: mu,nu,iG
  double precision              :: alpha
  double precision              :: r,dFxdr
  double precision,parameter    :: threshold = 1d-15

! Output variables

  double precision,intent(out)  :: Fx(nBas,nBas)

! Cx coefficient for Slater LDA exchange

  alpha = - (1d0/2d0)**(1d0/3d0)*(3d0/2d0)*(3d0/(4d0*pi))**(1d0/3d0)

! Compute LDA exchange matrix in the AO basis

  Fx(:,:) = 0d0
  do mu=1,nBas
    do nu=1,nBas
      do iG=1,nGrid

        r = max(0d0,rho(iG))

        if(r > threshold) then 

          dFxdr = 4d0/3d0*alpha*r**(1d0/3d0)

          Fx(mu,nu) = Fx(mu,nu) + weight(iG)*dFxdr*AO(mu,iG)*AO(nu,iG)

         end if

      end do
    end do
  end do

end subroutine lda_exchange_potential
