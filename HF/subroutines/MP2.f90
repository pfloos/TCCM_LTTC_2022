subroutine MP2(nBas,nO,nV,e,ERI)

! Compute the HF energy in the MO basis

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: i,j
  integer                       :: a,b
  double precision              :: eps
  double precision              :: EcMP2

! Output variables

! Compute MP2 correlation energy

  EcMP2 = 0d0

  do i=1,nO
    do j=1,nO
      do a=nO+1,nBas
        do b=nO+1,nBas

          eps = e(i) + e(j) - e(a) - e(b)
          EcMP2 = EcMP2 + ERI(i,j,a,b)*(2d0*ERI(i,j,a,b) - ERI(i,j,b,a))/eps

        end do
      end do
    end do
  end do

  write(*,'(A30,F15.10)') 'MP2 correlation energy',EcMP2

end subroutine MP2 
