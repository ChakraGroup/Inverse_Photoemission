module mod_rand_util
use mod_param,only: wp
use mod_rand
implicit none
private
public :: shuffle_vec
public :: calc_randmat
public :: calc_npt_per_bin

contains

!***********************************************!
subroutine shuffle_vec(n,v)
!***********************************************!
implicit none
integer,intent(in) :: n
real(8),intent(inout) :: v(n)
integer :: i,j
real(8) :: t
    do i=1,n
        j = rand_num(i,n)
        t = v(i)
        v(i) = v(j)
        v(j) = t
    end do
end subroutine shuffle_vec

!***********************************************!
subroutine rand_sample_from_1d(nbin,npt_tot,npt_per_bin, & 
                               bin_start,rand_vec)
!***********************************************!
implicit none
integer,intent(in) :: nbin
integer,intent(in) :: npt_tot
integer,intent(in) :: npt_per_bin(nbin)
real(wp),intent(in) :: bin_start(nbin+1)
real(wp),intent(out) :: rand_vec(npt_tot)

!--local--!
integer :: ibin
real(wp) :: xmin
real(wp) :: xmax
real(wp) :: r
integer  :: ipt
integer  :: jpt
integer  :: jpt_max
    if(sum(npt_per_bin) /= npt_tot) then
        write(*,*) "***ERROR sum(npt_per_bin) != npt_tot***"
        write(*,*) "npt_tot= ",npt_tot
        write(*,*) "sum(npt_per_bin)= ",sum(npt_per_bin)
        do ibin=1,nbin
            write(*,*) ibin,npt_per_bin(ibin)
        end do
        STOP "***FATAL ERROR in rand_sample_from_1d***"
    end if
    ipt = 0
    bin_loop: do ibin=1,nbin
        xmin = bin_start(ibin)
        xmax = bin_start(ibin+1)
        jpt = 0
        jpt_max = npt_per_bin(ibin)
        points_loop: do 
            r = real(rand_num(),wp)
            jpt = jpt + 1
            if(jpt > jpt_max) exit points_loop
            ipt = ipt + 1
            rand_vec(ipt) = xmin + ((xmax-xmin)*r)
            !--antithetic pair--!
            jpt = jpt + 1
            if(jpt > jpt_max) exit points_loop
            ipt = ipt + 1
            rand_vec(ipt) = xmin + ((xmax-xmin)*(1.0e0-r))
        end do points_loop
    end do bin_loop 
end subroutine rand_sample_from_1d


!ZZZ
!***********************************************!
subroutine calc_randmat(ndim,nbin,npt_tot,npt_per_bin_dim, & 
                        bin_start,randmat)
!***********************************************!
implicit none
integer,intent(in) :: ndim
integer,intent(in) :: nbin
integer,intent(in) :: npt_tot
integer,intent(in) :: npt_per_bin_dim(ndim,nbin)
real(wp),intent(in) :: bin_start(nbin+1)
real(wp),intent(out) :: randmat(npt_tot,ndim)
!--local--!
integer :: idim
    do idim=1,ndim
        call rand_sample_from_1d(nbin,npt_tot,npt_per_bin_dim(idim,:), & 
                                 bin_start,randmat(:,idim))
    end do
end subroutine calc_randmat




!***********************************************!
subroutine calc_npt_per_bin(nbin,target_npt,base_npt, & 
                            bin_start,wt,npt_per_bin)
!***********************************************!
implicit none
integer,intent(in) :: nbin
integer,intent(in) :: target_npt
integer,intent(in) :: base_npt
real(wp),intent(in) :: bin_start(nbin+1)
real(wp),intent(in) :: wt(nbin)
integer,intent(out) :: npt_per_bin(nbin)
!--local--!
integer :: ibin
integer :: rem_npt
integer :: curr_npt
integer :: diff
integer :: imax
integer :: ibuff(1)
!--check if wt are valid--!
    if( abs(sum(wt)-1.0e0_wp) > 1.0e-8_wp ) then
        write(*,*) "***ERROR: wt don't add up to 1***"
        write(*,*) "sum_wt= ",sum(wt)
        do ibin=1,nbin
            write(*,*) ibin,wt(ibin)
        end do
        STOP "***FATAL ERROR***"
    end if

!--first: all bins gets base number of points--!
    do ibin=1,nbin
        npt_per_bin(ibin) = base_npt 
    end do
!--calc remaining points--!
    rem_npt = target_npt - (base_npt*nbin)
!--assign remaining pt proportional to wt--!
    if(rem_npt > 0) then
        do ibin=1,nbin
            npt_per_bin(ibin) = npt_per_bin(ibin) + int(wt(ibin)*rem_npt)
        end do
    else
        write(*,*) "***WARNING: rem_npt is less than zero***"
    end if
!--check if we have the target number of pts--!
    curr_npt = sum(npt_per_bin)
!--if we have less points,then we give more pt to the max wt--!
!--if we have more points,then we take away pt randomly--!
    if(curr_npt < target_npt) then
        diff = target_npt - curr_npt
        ibuff = maxloc(wt)
        imax  = ibuff(1)
        npt_per_bin(imax) = npt_per_bin(imax) + diff
    else if(curr_npt > target_npt) then
        STOP "--PLEASE IMPLEMENT curr_npt > target_npt--"
    end if

    if(sum(npt_per_bin) /= target_npt) then
        write(*,*) "***ERROR: sum(npt_per_bin) != target_npt***"
        write(*,*) "target_npt= ",target_npt
        do ibin=1,nbin
            write(*,*) ibin,npt_per_bin(ibin)
        end do
        STOP "***FATAL ERROR***"
    end if

!   write(*,*) "--from calc_npt_per_bin--"
!   do ibin=1,nbin
!       write(*,*) ibin,npt_per_bin(ibin)
!   end do
end subroutine calc_npt_per_bin




end module mod_rand_util





