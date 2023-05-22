program test1
use mod_param,only: wp
use mod_rand
use mod_rand_util,only: shuffle_vec
implicit none
integer :: iseed

integer,parameter :: n=5
integer :: i
real(wp) :: v(n)
real(wp) :: wt(n)
integer :: npt(n)
integer :: npt_new(n)
integer :: tgt_npt
integer :: sum_npt

    tgt_npt = 1000
    call init_rand_with_sysclock()

    do i=1,n
        wt(i) = rand_num()
    end do
    wt = wt/sum(wt)
    call calc_npt_using_wt(5,tgt_npt,wt,npt)

    do i=1,n
        write(*,*) "--wt--",wt(i),npt(i)
    end do

    write(*,*) "sum_wt= ",sum(wt)
    write(*,*) "sum_npt= ",sum(npt)
    write(*,*) "tgt_npt= ",tgt_npt

end program test1

subroutine calc_npt_using_wt(nbin,tgt_npt,wt,npt)
use mod_param,only: wp
implicit none
integer,intent(in) :: nbin
integer,intent(in) :: tgt_npt
real(wp),intent(in) :: wt(nbin)
integer,intent(out) :: npt(nbin)
integer :: this_idx
integer :: buff(1)
real(wp) :: sum_npt
    npt(:)  = int(wt(:)*tgt_npt)
    sum_npt = sum(npt)
    if(sum_npt < tgt_npt) then
        buff     = maxloc(wt)
        this_idx = buff(1)
        npt(this_idx) = npt(this_idx) + (tgt_npt-sum_npt)
    else if (sum_npt > tgt_npt) then
        buff     = minloc(wt)
        this_idx = buff(1)
        npt(this_idx) = npt(this_idx) - (sum_npt-tgt_npt)
    end if
end subroutine calc_npt_using_wt






