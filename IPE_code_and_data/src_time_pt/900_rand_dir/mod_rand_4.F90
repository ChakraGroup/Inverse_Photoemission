!******************************************!
MODULE mod_rand
!******************************************!
implicit none
private
public :: init_rand_with_sysclock
public :: init_rand_with_seed
public :: rand_num
public :: rand_num_inf_inf
public :: box_muller
public :: gen_box_muller
public :: rand_sign
public :: irand_sign
public :: rand_unit_sphere

integer :: nseed

   interface rand_num
      module procedure rand_num_1
      module procedure rand_num_2
      module procedure irand_num_3
   end interface rand_num

CONTAINS

!*********************************************!
SUBROUTINE init_rand_with_sysclock()
!*********************************************!
implicit none
integer :: iseed
call system_clock(iseed)
write(*,*)'-----------------------------------'
write(*,*) "seed_for_random_stream= ",iseed
write(*,*) "tag2019_seed_for_random_stream= ",iseed
write(*,*)'-----------------------------------'
call init_rand_with_seed(iseed)
end SUBROUTINE init_rand_with_sysclock



!*********************************************!
SUBROUTINE init_rand_with_seed(iseed)
!*********************************************!
implicit none
integer,intent(in) :: iseed
integer :: nseed 
integer,allocatable ::  imat(:)

   nseed = 0
   call random_seed(size=nseed)
   allocate(imat(nseed))
   imat(:) = iseed
   call random_seed(put=imat)
   deallocate(imat)

END SUBROUTINE init_rand_with_seed


!*********************************************!
FUNCTION rand_num_1() result(ans)
!*********************************************!
implicit none
double precision :: ans
double precision :: a(1)
   
   call random_number(a)
   ans = a(1)

END FUNCTION rand_num_1

!*********************************************!
FUNCTION rand_num_2(xmin,xmax) result(ans)
!*********************************************!
implicit none
double precision,intent(in) :: xmin
double precision,intent(in) :: xmax
double precision :: ans

   ans = xmin + ((xmax-xmin)*rand_num_1())

END FUNCTION rand_num_2

!*********************************************!
FUNCTION irand_num_3(imin,imax) result(ians)
!*********************************************!
implicit none
integer,intent(in) :: imin
integer,intent(in) :: imax
integer :: ians
   ians = imin + nint((imax-imin)*rand_num_1())
END FUNCTION irand_num_3

!*********************************************!
FUNCTION rand_num_inf_inf() result(ans)
!*********************************************!
implicit none
double precision :: x
double precision :: y
double precision :: z
double precision :: ans

   x = rand_num_1()
   y = rand_num_1()
   z = y/(1.0d0-y)
   ans = z   
   if(x < 0.5) then
     ans = -z
   endif  

END FUNCTION rand_num_inf_inf

FUNCTION rand_sign() result(ans)
implicit none
double precision :: x
double precision :: ans
   x = rand_num_1()
   ans = 1.0d0
   if( x < 0.50d0 ) ans = -1.0d0
END FUNCTION rand_sign

FUNCTION irand_sign() result(ans)
implicit none
double precision :: x
integer :: ans
   x = rand_num_1()
   ans = 1
   if( x < 0.50d0 ) ans = -1
END FUNCTION irand_sign

! theta = [0,pi]
! phi = [0,2pi]
!*********************************************!
SUBROUTINE rand_unit_sphere(theta,phi)
!*********************************************!
implicit none
double precision,intent(out) :: theta,phi
double precision :: u,v
   u = rand_num_1()
   v = rand_num_1()
   phi   = 2.0d0 * dacos(-1.0d0) * u
   theta = dacos((2.0d0*v)-1.0d0)
END SUBROUTINE rand_unit_sphere

!*********************************************!
FUNCTION box_muller() result(ans) 
!*********************************************!
implicit none
double precision,parameter :: PI = 3.14159265358979323846264338327950d0
double precision :: ans 
double precision :: x1
double precision :: x2
   ans = 0.0d0
   x1  = rand_num()
   x2  = rand_num()
   ans = sqrt(-2.0*log(x1)) * cos(2.0*PI*x2)

END FUNCTION box_muller

!! Samples Exp[-alp*(x-A)^2]
SUBROUTINE gen_box_muller(alp,A,x,wt)
implicit none
double precision,intent(in) :: alp
double precision,intent(in) :: A
double precision,intent(out) :: x
double precision,intent(out) :: wt
double precision,parameter :: PI = 3.14159265358979323846264338327950d0
  wt = sqrt(pi/alp)
  x  = (box_muller()/sqrt(2.0e0*alp))+A
END SUBROUTINE gen_box_muller





END MODULE mod_rand

! program test
!  use mod_rand,only: init_rand_with_seed,rand_num
!  implicit none
!  integer :: i
!  integer :: iseed
!  integer :: imin,imax
!  double precision :: xmin,xmax
! 
! 
!     write(*,*) 'enter seed'
!     read(*,*) iseed
!     call init_rand_with_seed(iseed)
! 
!     imin = 10
!     imax = 100
! 
!     xmin = 1.0d0
!     xmax = 100.0d0
! 
!     call init_rand_with_seed(iseed)
!  do i=1,10
!     write(*,*) rand_num()
!  enddo
!  do i=11,20
!     write(*,*) rand_num()
!  enddo
!      write(*,*) rand_num(imin,imax)
!     write(*,*) rand_num(xmin,xmax)
!    
!
! end program test
