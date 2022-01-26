subroutine RKS(rung,nGrid,weight,nBas,AO,dAO,nO,S,T,V,Hc,ERI,X,ENuc,EKS)

! Perform a restricted Kohn-Sham calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: rung
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: dAO(3,nBas,nGrid)

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
  integer,parameter             :: n_diis = 1
  integer                       :: nSCF
  double precision              :: Conv
  double precision              :: ET,EV,EJ
  double precision              :: Ex
  double precision              :: Ec
  double precision,allocatable  :: e(:)
  double precision,allocatable  :: c(:,:),cp(:,:)
  double precision,allocatable  :: P(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: F(:,:),Fp(:,:)
  double precision,allocatable  :: Fx(:,:)
  double precision,allocatable  :: FxHF(:,:)
  double precision,allocatable  :: Fc(:,:)
  double precision,allocatable  :: error(:,:)
  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: F_diis(:,:)
  double precision,external     :: trace_matrix
  double precision,external     :: exchange_energy
  double precision,external     :: correlation_energy
  double precision,external     :: electron_number

  double precision,allocatable  :: rho(:)
  double precision,allocatable  :: drho(:,:)
  double precision              :: nEl

! Output variables

  double precision,intent(out)  :: EKS

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|      Restricted Kohn-Sham calculation        |'
  write(*,*)'************************************************'
  write(*,*)

!------------------------------------------------------------------------
! Rung of Jacob's ladder
!------------------------------------------------------------------------

  call select_rung(rung)

! Memory allocation

  allocate(e(nBas),c(nBas,nBas),cp(nBas,nBas),P(nBas,nBas),                       &
           J(nBas,nBas),F(nBas,nBas),Fp(nBas,nBas),Fx(nBas,nBas),FxHF(nBas,nBas), &
           Fc(nBas,nBas),error(nBas,nBas),rho(nGrid),drho(3,nGrid),             &
           error_diis(nBas*nBas,n_diis),F_diis(nBas*nBas,n_diis))

! Guess coefficients and eigenvalues

  F(:,:) = Hc(:,:)

! Initialization

  nSCF = 0
  Conv = 1d0
  nEl  = 0d0

  Ex   = 0d0
  Ec   = 0d0

  Fx(:,:)   = 0d0
  FxHF(:,:) = 0d0
  Fc(:,:)   = 0d0

  F_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------

  write(*,*)
  write(*,*)'------------------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') & 
            '|','#','|','E(KS)','|','Ex(KS)','|','Ec(KS)','|','Conv','|','nEl','|'
  write(*,*)'------------------------------------------------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

!   Increment 

    nSCF = nSCF + 1

!  Transform Fock matrix in orthogonal basis

    Fp = matmul(transpose(X),matmul(F,X))

!  Diagonalize Fock matrix to get eigenvectors and eigenvalues

    cp(:,:) = Fp(:,:)
    call diagonalize_matrix(nBas,cp,e)

!   Back-transform eigenvectors in non-orthogonal basis

    c = matmul(X,cp)

!   Compute density matrix 

    P(:,:)  = 2d0*matmul(c(:,1:nO),transpose(c(:,1:nO)))

!   Compute one-electron density and its gradient if necessary

!   call density(nGrid,nBas,P,AO,rho)
!   if(rung > 1 .and. rung < 666) call gradient_density(nGrid,nBas,P,AO,dAO,drho)

!   Build Coulomb repulsion
    
    call hartree_coulomb(nBas,P,ERI,J)

!   Compute exchange potential

    call exchange_potential(rung,nGrid,weight,nBas,P,ERI,AO,dAO,rho,drho,Fx,FxHF)

!   Compute correlation potential

!   call correlation_potential(rung,nGrid,weight,nBas,AO,dAO,rho,drho,Fc)

!   Build Fock operator
  
    F(:,:) = Hc(:,:) + J(:,:) + Fx(:,:) + Fc(:,:)

!   Check convergence 

    error = matmul(F,matmul(P,S)) - matmul(matmul(S,P),F)
    Conv  = maxval(abs(error))

!   DIIS extrapolation

    call DIIS_extrapolation(nBas*nBas,min(n_diis,nSCF),error_diis,F_diis,error,F)

!------------------------------------------------------------------------
!   Compute KS energy
!------------------------------------------------------------------------

!  Kinetic energy

    ET = trace_matrix(nBas,matmul(P,T))

!  Potential energy

    EV = trace_matrix(nBas,matmul(P,V))

!  Coulomb energy

    EJ = 0.5d0*trace_matrix(nBas,matmul(P,J))

!   Exchange energy

    Ex = exchange_energy(rung,nGrid,weight,nBas,P,FxHF,rho,drho)

!   Correlation energy

!   Ec = correlation_energy(rung,nGrid,weight,rho,drho,Ec)

    EKS = ET + EV + EJ + Ex + Ec

!   Check the grid accuracy by computing the number of electrons 

!   nEl = electron_number(nGrid,weight,rho)

!   Dump results

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') & 
      '|',nSCF,'|',EKS+ENuc,'|',Ex,'|',Ec,'|',Conv,'|',nEl,'|'
 
  enddo
  write(*,*)'------------------------------------------------------------------------------------------'
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

! Compute final KS energy

  call print_RKS(nBas,nO,e,C,ENuc,ET,EV,EJ,Ex,Ec,EKS)

end subroutine RKS
