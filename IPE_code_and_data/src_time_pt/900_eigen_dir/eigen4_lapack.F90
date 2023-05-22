!!**********************************************************************************************!!
SUBROUTINE EV_diag3(NBASIS,H,S,root,vec)
!!**********************************************************************************************!!
implicit none
integer, intent(in) :: NBASIS
double precision, intent(in) :: H(NBASIS,NBASIS)
double precision, intent(in) :: S(NBASIS,NBASIS)
double precision, intent(out) :: root(NBASIS)
double precision, intent(out) :: vec(NBASIS,NBASIS)

double precision :: A(NBASIS,NBASIS)
double precision :: B(NBASIS,NBASIS)

!!--------------------!!
!! Lapack variables
!!--------------------!!
integer :: INFO
integer,parameter :: ITYPE = 1
integer :: LWORK
integer :: iostat
double precision,allocatable :: WORK(:)
character(1),parameter  :: JOBZ = 'V'
character(1),parameter  :: UPLO = 'U'

!! FOR HP
	A = H
	B = S

! Get work array size
    allocate(WORK(1))
    LWORK = -1
	call DSYGV(ITYPE,JOBZ,UPLO,NBASIS,A,NBASIS,B,NBASIS,root,WORK,LWORK,INFO)
    LWORK = int(WORK(1))
    deallocate(WORK)

!   write(*,*) 
!   write(*,*) 'LWORK FOR EV_diag3 :',LWORK
!   if(iostat /= 0) then
!       write(*,*) 'IOSTAT : ',iostat
!       STOP 'ALLOCATION FAILED FOR WORK IN EV_diag3'
!   end if

    allocate(WORK(LWORK))
	call DSYGV(ITYPE,JOBZ,UPLO,NBASIS,A,NBASIS,B,NBASIS,root,WORK,LWORK,INFO)
	vec = A
    call eigenSign(NBASIS,vec)

!   write(*,*) 'INFO FOR EV_diag3 :',INFO
!   write(*,*) 

END SUBROUTINE EV_diag3

!!**********************************************************************************************!!
SUBROUTINE EV_diag4(NBASIS,H,root,vec)
!!**********************************************************************************************!!
implicit none
integer, intent(in) :: NBASIS
double precision, intent(in) :: H(NBASIS,NBASIS)
double precision, intent(out) :: root(NBASIS)
double precision, intent(out) :: vec(NBASIS,NBASIS)



integer :: i,j,k

!-----------------------!
! Lapack variables
!-----------------------!
integer :: INFO
integer :: LWORK
integer :: iostat
double precision,allocatable :: WORK(:)
character(1),parameter  :: JOBZ = 'V'
character(1),parameter  :: UPLO = 'U'

!! FOR HP 
    vec = H

! Get work array size
    allocate(WORK(1))
    LWORK = -1
	call DSYEV(JOBZ,UPLO,NBASIS,vec,NBASIS,root,WORK,LWORK,INFO)
    LWORK = int(WORK(1))
    deallocate(WORK)


    write(*,*) 
    write(*,*) 'LWORK FOR EV_diag4 :',LWORK
    allocate(WORK(LWORK),stat=iostat)

    if(iostat /= 0) then
        write(*,*) 'IOSTAT : ',iostat
        STOP 'ALLOCATION FAILED FOR WORK IN Ev_diag4'
    end if

	call DSYEV(JOBZ,UPLO,NBASIS,vec,NBASIS,root,WORK,LWORK,INFO)
    call eigenSign(NBASIS,vec)

    write(*,*) 'INFO FOR EV_diag4 :',INFO
    write(*,*) 

return
END SUBROUTINE EV_diag4



!!*********************************************!!
SUBROUTINE EV_diag5(NBASIS,H,S,root,vec)
!!*********************************************!!
! Solve generalized eigenvalue problem using UPPER PACKED storage mode
! for H and S
!
implicit none
integer, intent(in) :: NBASIS
double precision, intent(inout) :: H(NBASIS*(NBASIS+1)/2)
double precision, intent(inout) :: S(NBASIS*(NBASIS+1)/2)
double precision, intent(out) :: root(NBASIS)
double precision, intent(out) :: vec(NBASIS,NBASIS)

!!--------------------!!
!! Lapack variables
!!--------------------!!
integer :: INFO
integer :: LWORK
integer,parameter :: ITYPE = 1
double precision :: WORK(3*NBASIS)
character(1),parameter  :: JOBZ = 'V'
character(1),parameter  :: UPLO = 'U'
!
! Solves generalized eigenvalue problem
!
! Calling LAPACK subroutine
    call DSPGV(ITYPE,JOBZ,UPLO,NBASIS,H,S,root,vec,NBASIS,WORK,INFO)
    call eigenSign(NBASIS,vec)
END SUBROUTINE EV_diag5


!!*********************************************!!
SUBROUTINE EV_diag6(NB,H,root,vec)
!!*********************************************!!
implicit none
integer, intent(in) :: NB
double precision, intent(in) :: H(NB*(NB+1)/2)
double precision, intent(out) :: root(NB)
double precision, intent(out) :: vec(NB,NB)
double precision  :: T(NB*(NB+1)/2)

integer :: INFO
character(1) :: JOBZ
character(1) :: UPLO
double precision :: WORK(3*NB)

    JOBZ  = 'V'
    UPLO  = 'U'

    T(:) = H(:)
    call DSPEV(JOBZ,UPLO,NB,T,root,vec,NB,WORK,INFO)
    call eigenSign(NB,vec)

END SUBROUTINE EV_diag6










