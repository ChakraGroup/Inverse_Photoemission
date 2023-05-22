module mod_param
implicit none
integer,parameter :: wp = 8
real(wp),parameter :: underflow_tol = 1.0e-100_wp
real(wp),parameter :: overflow_tol  = 1.0e100_wp

contains

!***********************************************!
subroutine check_for_underflow(x)
!***********************************************!
implicit none
real(wp),intent(inout) :: x 
    if(abs(x) < underflow_tol) x = 0.0e0_wp
end subroutine check_for_underflow

end module mod_param

