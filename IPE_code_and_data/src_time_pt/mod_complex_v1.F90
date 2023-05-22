module mod_complex
use mod_param,only: wp
implicit none
private
integer,parameter :: cxn = 2
public :: cxn
public :: cx_pack_ab_format
public :: cx_pack_exp_format
public :: cx_unpack_ab_format
public :: cx_conj
public :: cx_sqmod


contains
! z = a + ib
!************************************************************************!
subroutine cx_pack_ab_format(a,b,cx)
!************************************************************************!
implicit none
real(wp),intent(in) :: a
real(wp),intent(in) :: b
real(wp),intent(out) :: cx(cxn,cxn)
    cx(1,:) = (/a,b/)
    cx(2,:) = (/-b,a/)
end subroutine cx_pack_ab_format

! z = r*exp[theta]
!************************************************************************!
subroutine cx_pack_exp_format(r,theta,cx)
!************************************************************************!
implicit none
real(wp),intent(in) :: r
real(wp),intent(in) :: theta
real(wp),intent(out) :: cx(cxn,cxn)
real(wp) :: a,b
    a = r*cos(theta)
    b = r*sin(theta)
    call cx_pack_ab_format(a,b,cx)
end subroutine cx_pack_exp_format

!************************************************************************!
subroutine cx_unpack_ab_format(cx,a,b)
!************************************************************************!
implicit none
real(wp),intent(in) :: cx(cxn,cxn)
real(wp),intent(out) :: a
real(wp),intent(out) :: b
    a = cx(1,1)
    b = cx(1,2)
end subroutine cx_unpack_ab_format

!************************************************************************!
subroutine cx_conj(cx,cx_star)
!************************************************************************!
implicit none
real(wp),intent(in) :: cx(cxn,cxn)
real(wp),intent(out) :: cx_star(cxn,cxn)
real(wp) :: a,b
    call cx_unpack_ab_format(cx,a,b)
    call cx_pack_ab_format(a,-b,cx_star)
end subroutine cx_conj

!************************************************************************!
subroutine cx_sqmod(cx,sqmod)
!************************************************************************!
implicit none
real(wp),intent(in) :: cx(cxn,cxn)
real(wp),intent(out) :: sqmod
real(wp) :: a,b
    call cx_unpack_ab_format(cx,a,b)
    sqmod = (a*a)+(b*b)
end subroutine cx_sqmod

end module mod_complex















