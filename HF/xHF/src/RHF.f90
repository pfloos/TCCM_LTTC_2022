subroutine RHF(nBas,nO,S,T,V,Hc,ERI,X,ENuc,EHF)

! Perform a restricted Hartree-Fock calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas

  integer,intent(in)            :: nO
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas) 
  double precision,intent(in)   :: X(nBas,nBas) 
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ENuc

! Local variables

  integer,parameter             :: maxSCF = 64
  double precision,parameter    :: thresh = 1d-5
  integer                       :: nSCF
  double precision              :: Conv
  double precision              :: Gap
  double precision              :: ET,EV,EJ
  double precision              :: EK
  double precision,allocatable  :: e(:)
  double precision,allocatable  :: c(:,:),cp(:,:)
  double precision,allocatable  :: P(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: F(:,:),Fp(:,:)
  double precision,allocatable  :: error(:,:)
  double precision,external     :: trace_matrix

! Output variables

  double precision,intent(out)  :: EHF

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|      Restricted Hartree-Fock calculation     |'
  write(*,*)'************************************************'
  write(*,*)

! Memory allocation

  allocate(e(nBas),c(nBas,nBas),cp(nBas,nBas),P(nBas,nBas),      &
           J(nBas,nBas),K(nBas,nBas),F(nBas,nBas),Fp(nBas,nBas), &
           error(nBas,nBas))

! Guess coefficients and eigenvalues

  F(:,:) = Hc(:,:)

! Initialization

  nSCF = 0
  Conv = 1d0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------

  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'| RHF calculation                                  |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','HF energy','|','Conv','|','HL Gap','|'
  write(*,*)'----------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

!   Increment 

    nSCF = nSCF + 1

! *** Exercise 2.1 *** !
!   Transform for the Fock matrix F in the orthogonal basis 
! ****************** !

! *** Exercise 2.2 *** !
!   Diagonalize F' to get MO coefficients (eigenvectors in the orthogonal basis) c' and MO energies (eigenvalues) e
! ****************** !

! *** Exercise 2.3 *** !
!   Back-transform the MO coefficients c in the original non-orthogonal basis
! ****************** !

! *** Exercise 2.4 *** !
!   Compute the density matrix P
! ****************** !

! *** Exercise 2.5 *** !
!   Compute the Hartree potential J
!   call hartree_potential(nBas,P,ERI,J)
! ****************** !

! *** Exercise 2.6 *** !
!   Compute the exchange potential K
!   call exchange_potential(nBas,P,ERI,K)
! ****************** !

!   Build Fock operator
  
    F(:,:) = Hc(:,:) + J(:,:) + K(:,:)

! *** Exercise 2.7 *** !
!   Compute the error vector and extract the convergence criterion
! ****************** !

!------------------------------------------------------------------------
!   Compute HF energy
!------------------------------------------------------------------------

! *** Exercise 3.1 *** !
!   Compute the kinetic energy
    ET = 0d0
! ****************** !

! *** Exercise 3.2 *** !
!   Compute the potential energy
    EV = 0d0
! ****************** !

! *** Exercise 3.3 *** !
!   Compute the Hartree energy
    EJ = 0d0
! ****************** !

! *** Exercise 3.4 *** !
!   Compute the exchange energy
    EK = 0d0
! ****************** !

!   Total HF energy

    EHF = ET + EV + EJ + EK

! *** Exercise 3.5 *** !
!   Compute HOMO-LUMO gap
    Gap = 0d0
! ****************** !


!   Dump results

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
      '|',nSCF,'|',EHF+ENuc,'|',Conv,'|',Gap,'|'
 
  enddo
  write(*,*)'----------------------------------------------------'
!------------------------------------------------------------------------
! End of SCF loop
!------------------------------------------------------------------------

! Did it actually converge?

  if(nSCF == maxSCF) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    stop

  endif

! Compute final HF energy

  call print_RHF(nBas,nO,e,C,ENuc,ET,EV,EJ,EK,EHF)

end subroutine RHF
