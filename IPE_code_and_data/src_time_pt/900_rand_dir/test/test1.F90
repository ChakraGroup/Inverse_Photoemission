program test1
use mod_param,only: wp
use mod_rand
use mod_rand_util,only: shuffle_vec
implicit none
integer :: iseed

integer,parameter :: n=5
integer :: i
real(wp) :: v(n)

call init_rand_with_sysclock()
write(*,*) rand_num()

    do i=1,n
        v(i) = real(i,wp)
    end do
    write(*,*) int(v)
    do i=1,10
        call shuffle_vec(n,v)
        write(*,*) int(v)
    end do

end program test1

