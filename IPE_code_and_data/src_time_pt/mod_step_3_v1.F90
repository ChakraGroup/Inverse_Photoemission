module mod_step_3
use mod_param,only: wp
use mod_complex,only: cxn,cx_pack_ab_format,cx_pack_exp_format,cx_sqmod
implicit none
private
real(wp),parameter :: zero = 0.0e0_wp
real(wp),parameter :: one  = 1.0e0_wp
real(wp),parameter :: ev_to_hartree = 0.03674930e0_wp

public :: driver_1

contains

!****************************************************************!
subroutine calc_optical_cross_sec_v1(nbas,omega_mo,omega_ext,wmat,eta, &
                             init_state_idx,cross_sec)
!****************************************************************!
implicit none
integer,intent(in)  :: nbas
real(wp),intent(in) :: omega_mo(nbas,nbas)
real(wp),intent(in) :: omega_ext
real(wp),intent(in) :: wmat(nbas,nbas)
real(wp),intent(in) :: eta
integer,intent(in) :: init_state_idx
real(wp),intent(out) :: cross_sec
!--local--!
real(wp) :: xsum
real(wp) :: omega_ik
real(wp) :: wik
real(wp) :: deno
integer :: k
!-------------------------------------------------------------------!
    xsum = 0.0e0_wp
    do k=1,nbas
        omega_ik = omega_mo(init_state_idx,k) 
        wik      = wmat(init_state_idx,k)
        if(omega_ik > 0.0e0_wp) then
            deno = ((omega_ik-omega_ext)**2) + (eta**2)
            xsum = xsum + ((wik**2)/deno)
        end if
    end do
    cross_sec = xsum
!-------------------------------------------------------------------!
end subroutine calc_optical_cross_sec_v1



!****************************************************************!
subroutine core_v3(nbas,num_omega)
!****************************************************************!
implicit none
integer,intent(in) :: nbas
integer,intent(in) :: num_omega
!--system info--!
integer :: t_nbas
real(wp) :: moeng(nbas)
real(wp) :: omega_mo(nbas,nbas)
real(wp) :: wmat(nbas,nbas)
integer  :: i,j
integer  :: t_i,t_j
!--PT--!
real(wp) :: tmin
real(wp) :: tmax
real(wp) :: tstep
real(wp) :: tcurr
real(wp) :: cx_bold(nbas,cxn,cxn)
real(wp) :: cx_bnew(nbas,cxn,cxn)
real(wp) :: omega_ext
real(wp) :: eta
real(wp) :: prob(nbas)
integer :: init_state_idx
real(wp) :: omega_min,omega_max,omega_step
real(wp) :: cross_sec
real(wp) :: csmat(num_omega,nbas)
integer :: iom
real(wp) :: omega_ext_list(num_omega)
integer  :: target_i,target_j
real(wp) :: omega_delta_ij
!--ecap--!
real(wp) :: ecap_prob(nbas)
!--diff--!
real(wp) :: diff,min_diff
integer  :: this_i,this_j
!------------------------------------------------------------------!
    open(unit=100,file='inp_dipole_info.txt',action='read')
        read(100,*) !title
        read(100,*) t_nbas
        if(t_nbas /= nbas) STOP "***ERROR:nbas in core_v3 is different from nbas in file inp_dipole_info.txt***"
        read(100,*) !title
        do i=1,nbas
        do j=i,nbas
            read(100,*) t_i,t_j,moeng(i),moeng(j),wmat(i,j)
            wmat(j,i) = wmat(i,j)
        end do
        end do
    close(100)
!------------------------------------------------------------------!
    open(unit=100,file='inp_ecap_info.txt',action='read')
        read(100,*) !title
        read(100,*) t_nbas
        if(t_nbas /= nbas) STOP "***ERROR:nbas in core_v3 is different from nbas in file inp_dipole_info.txt***"
        read(100,*) !title
        do i=1,nbas
            read(100,*) t_i,moeng(i),ecap_prob(i)
        end do
    close(100)
!------------------------------------------------------------------!
    open(unit=100,file='inp_step_3.txt',action='read') 
        read(100,*) !title
        read(100,*) target_i,target_j,omega_delta_ij
    close(100)
!------------------------------------------------------------------!
!--normalize ecap_prob--!
    ecap_prob(:) = ecap_prob(:)/sum(ecap_prob)
!------------------------------------------------------------------!
    write(*,*) 
    write(*,*) "--normalized capture prob--" 
    do i=1,nbas
        write(*,*) "idx_moeng_prob ",i,moeng(i),ecap_prob(i)
    end do
!------------------------------------------------------------------!
    do i=1,nbas
        wmat(i,i) = 0.0e0_wp
    end do
!------------------------------------------------------------------!
! this is only for de-excitation
! do not allow excitation. 
    do i=1,nbas
    do j=1,nbas
        omega_mo(i,j) = moeng(i)-moeng(j)
        if(omega_mo(i,j) <= 0.0e0_wp) then
            omega_mo(i,j) = 0.0e0_wp
            wmat(i,j)     = 0.0e0_wp
        end if
    end do
    end do
!------------------------------------------------------------------!
    if(omega_mo(target_i,target_j) > 0.0e0_wp) then
       omega_mo(target_i,target_j) = omega_mo(target_i,target_j) + omega_delta_ij
    end if
    if(omega_mo(target_j,target_i) > 0.0e0_wp) then
       omega_mo(target_j,target_i) = omega_mo(target_j,target_i) + omega_delta_ij
    end if
!------------------------------------------------------------------!
    do i=1,nbas
    do j=1,nbas
        if(omega_mo(i,j) <= 0.0e0_wp) then
            omega_mo(i,j) = 0.0e0_wp
            wmat(i,j)     = 0.0e0_wp
        end if
    end do
    end do
!------------------------------------------------------------------!
    omega_ext = 0.0e0_wp
    eta       = 1.0e-6_wp
!------------------------------------------------------------------!
    omega_min  =  0.0e0_wp   * ev_to_hartree
    omega_max  =  6.0e0_wp   * ev_to_hartree
    omega_step =  0.10e0_wp  * ev_to_hartree 

    omega_loop1: do iom=1,num_omega
        omega_ext      = omega_min + ((iom-1)*(omega_max-omega_min)/(num_omega-1))
        omega_ext_list(iom) = omega_ext
        init_state_loop: do i=1,nbas
            init_state_idx = i
            call calc_optical_cross_sec_v1(nbas,omega_mo,omega_ext,wmat,eta, &
                           init_state_idx,cross_sec)
            csmat(iom,i) = cross_sec * ecap_prob(i)
        end do init_state_loop
    end do omega_loop1
!------------------------------------------------------------------!
!--normalize--!
    csmat(:,:) = csmat(:,:)/sum(csmat)
!------------------------------------------------------------------!
    open(unit=100,file='out_ips_spectra.txt',action='write')
        write(100,*) "OMEGA NORMALIZED_IPS_PROB"
        write(*,*) 
        do iom=1,num_omega
            omega_ext = omega_ext_list(iom)
            write(*,*) "--omega_tot_cross_sec--", omega_ext,sum(csmat(iom,:))
            write(100,*) omega_ext,sum(csmat(iom,:))
        end do
    close(100)

    open(unit=100,file='out_tot_mo_cross_sec.txt',action='write')
        write(100,*) "mo_idx,E-ELUMO,total_ips_cross_sec"
        write(*,*) 
        do i=1,nbas
            write(*,*) "--state_tot_cross_sec--", i,moeng(i),sum(csmat(:,i))
            write(100,*)  i,moeng(i)-moeng(1),sum(csmat(:,i))
        end do
    close(100)
!------------------------------------------------------------------!
    write(*,*)
    write(*,*) "--mapping omega_ext to state_info--"
    do iom=1,num_omega
        omega_ext = omega_ext_list(iom)
        this_i   = 0
        this_j   = 0
        min_diff = 100.0e0_wp
        do i=1,nbas
        do j=1,nbas
            diff = abs(omega_ext-omega_mo(i,j))
            if(diff < min_diff) then
                this_i = i
                this_j = j
                min_diff =  diff
            end if
        end do 
        end do 
        write(*,*) "omega_ext-omega0--",omega_ext,omega_mo(this_i,this_j),this_i,this_j
    end do
!------------------------------------------------------------------!
end subroutine core_v3















!****************************************************************!
subroutine driver_1
!****************************************************************!
implicit none
integer :: nbas
integer :: num_omega
!-------------------------------------------------------------------!
    open(unit=100,file='inp_dipole_info.txt',action='read')
        read(100,*) !title
        read(100,*) nbas
    close(100)
!-------------------------------------------------------------------!
    num_omega = 100
    call core_v3(nbas,num_omega)
!-------------------------------------------------------------------!
end subroutine driver_1









end module mod_step_3

