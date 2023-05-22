module mod_cano_ortho
implicit none
private
integer,parameter :: wp=8
public :: calc_cano_ortho

contains
! X.T S X = 1
! NOTE: X.T X != 1
!!************************************************!!
subroutine calc_cano_ortho(nb,eigen_tol,smat,X)
!!************************************************!!
implicit none
integer,intent(in) :: nb
real(wp),intent(in) :: eigen_tol
real(wp),intent(in) :: smat(nb,nb)
real(wp),intent(out) :: X(nb,nb)
real(wp) :: root(nb)
real(wp) :: vec(nb,nb)
real(wp) :: diag(nb,nb)
integer :: i

    diag = 0.0e0_wp
    call EV_diag4(nb,smat,root,vec)
    do i=1,nb
        if(root(i) < eigen_tol) then
            write(*,*) "***WARNING: SMAT NOT positive definite in calc_cano_ortho***" 
            write(*,*) "***setting sqrt[1/s(i)] = 0***" 
            write(*,*) "i,s(i),tol ",i,root(i),eigen_tol
            diag(i,i) = 0.0e0_wp
        else
            diag(i,i) = 1.0e0_wp/sqrt(root(i))
        end if
    end do
    X = matmul(vec,diag)
end subroutine calc_cano_ortho

end module mod_cano_ortho

