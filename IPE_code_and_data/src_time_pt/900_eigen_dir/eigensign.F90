SUBROUTINE eigenSign(N,V)
implicit none
integer,intent(in) :: N
double precision,intent(inout) :: V(N,N)
double precision :: U(N,N)

!
! This subroutine imposes the sign convention
! to all the eigenvectors
!
! The component of the vector with the highest
! magnitude will be defined as positive
!

integer :: i,j

    U = dabs(V)
    do i=1,N
        j = maxloc(U(1:N,i),dim=1)
        if(V(j,i) < 0.0d0) V(1:N,i) = -V(1:N,i)
    end do


END SUBROUTINE eigenSign


