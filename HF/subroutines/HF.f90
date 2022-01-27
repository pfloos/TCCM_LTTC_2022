subroutine HF(nBas,nO,c,Hc,ERI)

! Compute the HF energy in the MO basis

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nO
  double precision,intent(in)   :: c(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: i,j
  integer                       :: mu,nu
  double precision,allocatable  :: h(:)
  double precision              :: EHF

! Output variables

! Memory allocation
  
  allocate(h(nO))

! AO to MO transformation for Hc
  
  h(:) = 0d0

  do i=1,nO
    do mu=1,nBas
      do nu=1,nBas
        h(i) = h(i) + c(mu,i)*Hc(mu,nu)*c(nu,i)
      end do
    end do
  end do

! Compute HF energy

  EHF = 0d0

  do i=1,nO
    EHF = EHF + 2d0*h(i)
    do j=1,nO
      EHF = EHF + 2d0*ERI(i,j,i,j) - ERI(i,j,j,i)
    end do
  end do

  write(*,'(A30,F15.10)') 'HF energy',EHF

end subroutine HF
