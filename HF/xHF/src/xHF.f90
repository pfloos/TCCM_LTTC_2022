program xHF

! Perform a Hartree-Fock calculation

  include 'parameters.h'

  integer                       :: nAt,nBas,nEl,nO,nV
  double precision              :: ENuc,EHF

  double precision,allocatable  :: ZNuc(:)
  double precision,allocatable  :: rAt(:,:)

  integer                       :: nShell
  integer,allocatable           :: TotAngMomShell(:)
  integer,allocatable           :: KShell(:)
  double precision,allocatable  :: CenterShell(:,:)
  double precision,allocatable  :: DShell(:,:)
  double precision,allocatable  :: ExpShell(:,:)

  double precision,allocatable  :: S(:,:)
  double precision,allocatable  :: T(:,:)
  double precision,allocatable  :: V(:,:)
  double precision,allocatable  :: Hc(:,:)
  double precision,allocatable  :: X(:,:)
  double precision,allocatable  :: ERI(:,:,:,:)

  double precision              :: start_HF,end_HF,t_HF

! Hello World

  write(*,*)
  write(*,*) '*******************************'
  write(*,*) '* TCCM winter school 2020: HF *'
  write(*,*) '*******************************'
  write(*,*)

!------------------------------------------------------------------------
! Read input information
!------------------------------------------------------------------------

! Read number of atoms, number of electrons of the system
! nO   = number of occupied orbitals
! nV   = number of virtual orbitals (see below)
! nBas = number of basis functions (see below)
!      = nO + nV

  call read_molecule(nAt,nEl,nO)
  allocate(ZNuc(nAt),rAt(nAt,3))

! Read geometry

  call read_geometry(nAt,ZNuc,rAt,ENuc)

  allocate(CenterShell(maxShell,3),TotAngMomShell(maxShell),KShell(maxShell), &
           DShell(maxShell,maxK),ExpShell(maxShell,maxK))

!------------------------------------------------------------------------
! Read basis set information
!------------------------------------------------------------------------

  call read_basis(nAt,rAt,nBas,nO,nV,nShell,TotAngMomShell,CenterShell,KShell,DShell,ExpShell)

!------------------------------------------------------------------------
! Read one- and two-electron integrals
!------------------------------------------------------------------------

! Memory allocation for one- and two-electron integrals

  allocate(S(nBas,nBas),T(nBas,nBas),V(nBas,nBas),Hc(nBas,nBas),X(nBas,nBas), &
           ERI(nBas,nBas,nBas,nBas))

! Read integrals

  call read_integrals(nBas,S,T,V,Hc,ERI)  

! Orthogonalization X = S^(-1/2)

! *** Exercise 1 *** !
! Compute the orthogonalization matrix X = S^(-1/2)
! call orthogonalization_matrix(nBas,S,X)
! ****************** !

!------------------------------------------------------------------------
! Compute restricted HF energy
!------------------------------------------------------------------------

    call cpu_time(start_HF)
    call RHF(nBas,nO,S,T,V,Hc,ERI,X,ENuc,EHF)
    call cpu_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for HF = ',t_HF,' seconds'
    write(*,*)

!------------------------------------------------------------------------
! End of xHF
!------------------------------------------------------------------------
end program xHF
