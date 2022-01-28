subroutine lda_correlation_potential(nGrid,weight,nBas,AO,rho,Fc)

! Compute LDA correlation potential

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
  double precision              :: a,d
  double precision              :: r,ra,rb,ec,dFcdr
  double precision,parameter    :: threshold = 1d-15

! Output variables

  double precision,intent(out)  :: Fc(nBas,nBas)

! Coefficients for Wigner's LDA correlation

  a = 0.04918d0
  d = 0.349d0

! Compute LDA correlation matrix in the AO basis

  Fc(:,:) = 0d0

  do mu=1,nBas
    do nu=1,nBas
      do iG=1,nGrid

        ra = max(0d0,rho(iG)/2d0)
        rb = max(0d0,rho(iG)/2d0)
        r  = ra + rb
      
        if(r > threshold) then

          ec    = ra*rb/(r + d*r**(2d0/3d0))
          dFcdr = ec*(d/(3d0*r**(4d0/3d0)*(1d0 + d*r**(-1d0/3d0))) - 1d0/r + 1d0/ra)
  
          Fc(mu,nu) = Fc(mu,nu) + weight(iG)*AO(mu,iG)*AO(nu,iG)*dFcdr

        end if

      end do
    end do
  end do

  Fc(:,:) = -4d0*a*Fc(:,:)

end subroutine lda_correlation_potential
