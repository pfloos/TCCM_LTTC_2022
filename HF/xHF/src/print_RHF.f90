subroutine print_RHF(nBas,nO,e,C,ENuc,ET,EV,EJ,Ex,EHF)

! Print one- and two-electron energies and other stuff for RHF calculation

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas,nO
  double precision,intent(in)        :: e(nBas),c(nBas,nBas),ENuc,ET,EV,EJ,Ex,EHF

  integer                            :: HOMO,LUMO
  double precision                   :: Gap

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = e(LUMO) - e(HOMO)

! Dump results


  write(*,*)
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A32)')           ' Summary              '
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A32,1X,F16.10)') ' One-electron energy  ',ET + EV
  write(*,'(A32,1X,F16.10)') ' Kinetic      energy  ',ET
  write(*,'(A32,1X,F16.10)') ' Potential    energy  ',EV
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A32,1X,F16.10)') ' Two-electron energy  ',EJ + Ex
  write(*,'(A32,1X,F16.10)') ' Coulomb      energy  ',EJ
  write(*,'(A32,1X,F16.10)') ' Exchange     energy  ',Ex
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A32,1X,F16.10)') ' Electronic   energy  ',EHF
  write(*,'(A32,1X,F16.10)') ' Nuclear   repulsion  ',ENuc
  write(*,'(A32,1X,F16.10)') ' Hartree-Fock energy  ',EHF + ENuc
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A36,F13.6)')     ' HF HOMO      energy (eV):',e(HOMO)*autoeV
  write(*,'(A36,F13.6)')     ' HF LUMO      energy (eV):',e(LUMO)*autoev
  write(*,'(A36,F13.6)')     ' HF HOMO-LUMO gap    (eV):',Gap*autoev
  write(*,'(A50)')           '---------------------------------------'
  write(*,*)

! Print results

  write(*,'(A50)') '---------------------------------------'
  write(*,'(A50)') 'Hartree-Fock orbital coefficients      '
  write(*,'(A50)') '---------------------------------------'
  call matout(nBas,nBas,C)
  write(*,*)
  write(*,'(A50)') '---------------------------------------'
  write(*,'(A50)') ' Hartree-Fock orbital energies         '
  write(*,'(A50)') '---------------------------------------'
  call matout(nBas,1,e)
  write(*,*)

end subroutine print_RHF


