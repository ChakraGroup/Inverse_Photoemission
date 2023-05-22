implicit none
integer :: i(5)
integer :: A(5)
    i    = 0
    i(1) = 1
    i(2) = 2
    A = 1
    A(i) = 100
    write(*,*) A
end


