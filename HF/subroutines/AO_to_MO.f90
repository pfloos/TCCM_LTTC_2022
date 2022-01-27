subroutine AO_to_MO(nBas,c,ERI_AO,ERI_MO)

! Expression of bi-electronic integrals in the MO basis set

  implicit none
!  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: c(nBas,nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: mu,nu,la,si
  integer                       :: p,q,r,s
  double precision,allocatable  :: scr(:,:,:,:)

! Output variables

  double precision,intent(out)  :: ERI_MO(nBas,nBas,nBas,nBas)

! Memory allocation

  allocate(scr(nBas,nBas,nBas,nBas))
  
!--------------------------------------
! AO to MO transformation starts here !
!--------------------------------------

! Transform 4th index

  scr(:,:,:,:) = 0d0

  do mu=1,nBas
     do nu=1,nBas
       do la=1,nBas
         do si=1,nBas

           do s=1,nBas

             scr(mu,nu,la,s) = scr(mu,nu,la,s) + c(si,s)*ERI_AO(mu,nu,la,si)

           enddo    

         enddo
      enddo
    enddo
  enddo

! Transform 3rd index

  ERI_MO(:,:,:,:) = 0d0

  do mu=1,nBas
    do nu=1,nBas
      do la=1,nBas

        do r=1,nBas
          do s=1,nBas

            ERI_MO(mu,nu,r,s) = ERI_MO(mu,nu,r,s) + c(la,r)*scr(mu,nu,la,s)

          enddo
        enddo    

      enddo
    enddo
  enddo

! Transform 2nd index

  scr(:,:,:,:) = 0d0

  do mu=1,nBas
    do nu=1,nBas

      do q=1,nBas
        do r=1,nBas
          do s=1,nBas

            scr(mu,q,r,s) = scr(mu,q,r,s) + c(nu,q)*ERI_MO(mu,nu,r,s)

          enddo
        enddo
      enddo    

    enddo
  enddo

! Transform 1st (and last) index

  ERI_MO(:,:,:,:) = 0d0

  do mu=1,nBas

    do p=1,nBas
      do q=1,nBas
        do r=1,nBas
          do s=1,nBas

            ERI_MO(p,q,r,s) = ERI_MO(p,q,r,s) + c(mu,p)*scr(mu,q,r,s)

          enddo
        enddo
      enddo
    enddo

  enddo

end subroutine AO_to_MO
