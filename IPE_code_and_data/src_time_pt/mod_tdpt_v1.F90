module mod_tdpt
use mod_param,only: wp
use mod_complex,only: cxn,cx_pack_ab_format,cx_pack_exp_format,cx_sqmod
implicit none
private
real(wp),parameter :: zero = 0.0e0_wp
real(wp),parameter :: one  = 1.0e0_wp
real(wp),parameter :: ev_to_hartree = 0.03674930e0_wp

public :: driver_1

contains

! v(t) = v * cos(wt)
!****************************************************************!
subroutine calc_single_state_bval(tcurr,omega_mo,omega_ext,vext_coupl_elem,eta,cx_b)
!****************************************************************!
implicit none
real(wp),intent(in) :: tcurr
real(wp),intent(in) :: omega_mo
real(wp),intent(in) :: omega_ext
real(wp),intent(in) :: vext_coupl_elem
real(wp),intent(in) :: eta
real(wp),intent(out) :: cx_b(cxn,cxn)
!--local--!
real(wp) :: theta
real(wp) :: deno
real(wp) :: cx_one(cxn,cxn)
real(wp) :: cx_pure_img(cxn,cxn)
real(wp) :: cx_exp_plus(cxn,cxn)
real(wp) :: cx_term_plus(cxn,cxn)
real(wp) :: cx_exp_minus(cxn,cxn)
real(wp) :: cx_term_minus(cxn,cxn)
real(wp) :: cx_term_vext(cxn,cxn)
!-------------------------------------------------------------------!
    call cx_pack_ab_format(one,zero,cx_one)
!-------------------------------------------------------------------!
! first: w(if) + w_ext
    theta = tcurr*(omega_mo + omega_ext)
    deno  = omega_mo + omega_ext + eta
    call cx_pack_exp_format(one,theta,cx_exp_plus)
    cx_term_plus = cx_one - cx_exp_plus
    cx_term_plus = cx_term_plus * (1.0e0_wp/deno)
    !write(*,*) "plus_term= ",cx_term_plus(1,:)
!-------------------------------------------------------------------!
!second: w(if)-w_ext
    theta = tcurr*(omega_mo - omega_ext)
    deno  = omega_mo - omega_ext + eta
    call cx_pack_exp_format(one,theta,cx_exp_minus)
    cx_term_minus = cx_one - cx_exp_minus
    cx_term_minus = cx_term_minus * (1.0e0_wp/deno)
    !write(*,*) "minus_term= ",cx_term_minus(1,:)
!-------------------------------------------------------------------!
! v(t) = v * cos(wt)
    cx_term_vext = (cx_term_plus + cx_term_minus)*(vext_coupl_elem*0.50e0_wp) 
    !write(*,*) "vext_term= ",cx_term_vext(1,:)
!-------------------------------------------------------------------!
    call cx_pack_ab_format(zero,one,cx_pure_img)
    cx_b = -1.0e0_wp * matmul(cx_pure_img,cx_term_vext)
!-------------------------------------------------------------------!
end subroutine calc_single_state_bval

!****************************************************************!
subroutine calc_multi_state_bvec(tcurr, &
                                 nbas,omega_mo,omega_ext,wmat,eta, &
                                 cx_bold,cx_bnew,prob)
!****************************************************************!
implicit none
real(wp),intent(in) :: tcurr
integer,intent(in)  :: nbas
real(wp),intent(in) :: omega_mo(nbas,nbas)
real(wp),intent(in) :: omega_ext
real(wp),intent(in) :: wmat(nbas,nbas)
real(wp),intent(in) :: eta
real(wp),intent(in) :: cx_bold(nbas,cxn,cxn)
real(wp),intent(out) :: cx_bnew(nbas,cxn,cxn)
real(wp),intent(out) :: prob(nbas)
!--local--!
integer :: i,j
real(wp) :: cx_ij(cxn,cxn)
real(wp) :: norm_prob
!-------------------------------------------------------------------!
    do i=1,nbas
        cx_bnew(i,:,:) = 1.0e0_wp
        do j=1,nbas
            call calc_single_state_bval(tcurr,omega_mo(i,j),omega_ext,wmat(i,j),eta,cx_ij)
            cx_bnew(i,:,:) = cx_bnew(i,:,:) + matmul(cx_ij,cx_bold(j,:,:))
        end do 
    end do
!-------------------------------------------------------------------!
!--normalize--!
    do i=1,nbas
        call cx_sqmod(cx_bnew(i,:,:),prob(i))
    end do
    norm_prob = sum(prob)
    cx_bnew   = cx_bnew/sqrt(norm_prob)
    prob      = prob/norm_prob
!-------------------------------------------------------------------!
!-------------------------------------------------------------------!
end subroutine calc_multi_state_bvec


!make sure that: omega(init_state_idx,final_state) > 0 
! all transitions with omega < 0 will be ignored
!****************************************************************!
subroutine calc_transition_prob_v1(tcurr, &
                                nbas,omega_mo,omega_ext,wmat,eta, &
                                init_state_idx,prob)
!****************************************************************!
implicit none
real(wp),intent(in) :: tcurr
integer,intent(in)  :: nbas
real(wp),intent(in) :: omega_mo(nbas,nbas)
real(wp),intent(in) :: omega_ext
real(wp),intent(in) :: wmat(nbas,nbas)
real(wp),intent(in) :: eta
integer,intent(in) :: init_state_idx
real(wp),intent(out) :: prob(nbas)
!--local--!
integer :: i,j
real(wp) :: cx_bvec(nbas,cxn,cxn)
real(wp) :: norm_prob
!-------------------------------------------------------------------!
    cx_bvec = 0.0e0_wp
    call cx_pack_ab_format(one,zero,cx_bvec(init_state_idx,:,:))
!-------------------------------------------------------------------!
    do j=1,nbas
        if(omega_mo(init_state_idx,j) > 0.0e0_wp) then
            call calc_single_state_bval(tcurr,omega_mo(init_state_idx,j), & 
                                    omega_ext,wmat(init_state_idx,j),eta,cx_bvec(j,:,:))
        end if
    end do 
!-------------------------------------------------------------------!
!--normalize--!
    prob(:) = 0.0e0_wp
    do i=1,nbas
        call cx_sqmod(cx_bvec(i,:,:),prob(i))
    end do
    norm_prob = sum(prob)
    if(norm_prob == 0.0e0_wp) STOP "***ERROR: total_prob is zero***"
    cx_bvec   = cx_bvec/sqrt(norm_prob)
    prob      = prob/norm_prob
!-------------------------------------------------------------------!
end subroutine calc_transition_prob_v1

!****************************************************************!
subroutine calc_transition_prob_v2(tcurr, &
                                nbas,omega_mo,omega_ext,wmat,eta, &
                                init_state_idx,prob)
!****************************************************************!
implicit none
real(wp),intent(in) :: tcurr
integer,intent(in)  :: nbas
real(wp),intent(in) :: omega_mo(nbas,nbas)
real(wp),intent(in) :: omega_ext
real(wp),intent(in) :: wmat(nbas,nbas)
real(wp),intent(in) :: eta
integer,intent(in) :: init_state_idx
real(wp),intent(out) :: prob(nbas)
!--local--!
integer :: i,j
real(wp) :: cx_bvec(nbas,cxn,cxn)
real(wp) :: norm_prob
real(wp) :: x,y,z
!-------------------------------------------------------------------!
    prob(:) = 0.0e0_wp
    prob(init_state_idx) = 1.0e0_wp
    do j=1,nbas
        if(omega_mo(init_state_idx,j) > 0.0e0_wp) then
            x = omega_mo(init_state_idx,j)-omega_ext
            if(abs(x) < eta .and. x >= 0.0e0_wp) x = eta
            if(abs(x) < eta .and. x < 0.0e0_wp) x  = -eta
            y = x*tcurr*0.50e0_wp
            z = (sin(y)/y)**2
            prob(j) = (wmat(init_state_idx,j)**2) * y * (tcurr**2)
        end if
    end do 
!-------------------------------------------------------------------!
!--normalize--!
    norm_prob = sum(prob)
    if(norm_prob == 0.0e0_wp) STOP "***ERROR: total_prob is zero***"
    prob = prob/norm_prob
!-------------------------------------------------------------------!
end subroutine calc_transition_prob_v2

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
    cross_sec = omega_ext*xsum
!-------------------------------------------------------------------!
end subroutine calc_optical_cross_sec_v1



!****************************************************************!
subroutine core_v1
!****************************************************************!
implicit none
!--bval--!
real(wp) :: omega_mo
real(wp) :: omega_ext
real(wp) :: vext_coupl_elem
real(wp) :: eta
real(wp) :: cx_b(cxn,cxn)
real(wp) :: prob_b
!--PT--!
real(wp) :: tmin
real(wp) :: tmax
real(wp) :: tstep
real(wp) :: tcurr
!------------------------------------------------------------------!
    omega_mo  = 1.0e0_wp
    omega_ext = 0.0e0_wp
    vext_coupl_elem = 0.10e0_wp
    eta = 1.0e-6_wp
    cx_b = 0.0e0_wp
!------------------------------------------------------------------!
    tmin = 0.0e0_wp
    tmax = 10.0e0_wp
    tstep = 0.10e0_wp
    tcurr = tmin
    timeloop1: do
        if(tcurr > tmax) exit timeloop1
        call calc_single_state_bval(tcurr,omega_mo,omega_ext,vext_coupl_elem,eta,cx_b)
        call cx_sqmod(cx_b,prob_b)
        write(*,*) "--t_prob--",tcurr,prob_b
        tcurr = tcurr + tstep
    end do timeloop1
!------------------------------------------------------------------!
end subroutine core_v1




!****************************************************************!
subroutine core_v2(nbas,num_omega)
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
!------------------------------------------------------------------!
    open(unit=100,file='inp_dipole_info.txt',action='read')
        read(100,*) !title
        read(100,*) t_nbas
        if(t_nbas /= nbas) STOP "***ERROR:nbas in core_v2 is different from nbas in file inp_dipole_info.txt***"
        read(100,*) !title
        do i=1,nbas
        do j=i,nbas
            read(100,*) t_i,t_j,moeng(i),moeng(j),wmat(i,j)
            wmat(j,i) = wmat(i,j)
        end do
        end do
    close(100)
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
    init_state_idx = 12
    omega_ext = 0.0e0_wp
    eta       = 1.0e-6_wp
    cx_bold(:,:,:) = 0.0e0_wp
    call cx_pack_ab_format(one,zero,cx_bold(init_state_idx,:,:))
!------------------------------------------------------------------!
    omega_min  =  3.0e0_wp * ev_to_hartree
    omega_max  =  6.0e0_wp * ev_to_hartree
    omega_step =  0.10e0_wp * ev_to_hartree 

    omega_loop1: do iom=1,num_omega
        omega_ext      = omega_min + ((iom-1)*(omega_max-omega_min)/(num_omega-1))
        omega_ext_list(iom) = omega_ext
        init_state_loop: do i=1,nbas
            init_state_idx = i
            call calc_optical_cross_sec_v1(nbas,omega_mo,omega_ext,wmat,eta, &
                           init_state_idx,cross_sec)
            csmat(iom,i) = cross_sec
        end do init_state_loop
    end do omega_loop1

    write(*,*) 
    do iom=1,num_omega
        omega_ext = omega_ext_list(iom)
        write(*,*) "--omega_tot_cross_sec--", omega_ext,sum(csmat(iom,:))
    end do

    write(*,*) 
    do i=1,nbas
        write(*,*) "--state_tot_cross_sec--", i,moeng(i),sum(csmat(:,i))
    end do
!------------------------------------------------------------------!
    open(unit=100,file='out_csmat.txt',action='write')
        write(100,*) "num_omega nmo"
        write(100,*) num_omega,nbas
        write(100,*) "list_of_omega_ext"
        write(100,*) omega_ext_list(:)
        write(100,*) "cross_sec_mat(iomega_ext,imo)"
        write(100,*) csmat(:,:)
    close(100)
!------------------------------------------------------------------!
end subroutine core_v2

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
!--ecap--!
real(wp) :: ecap_prob(nbas)
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
!--normalize ecap_prob--!
    ecap_prob(:) = ecap_prob(:)/sum(ecap_prob)
!------------------------------------------------------------------!
    write(*,*) 
    write(*,*) "--normalized capture prob--" 
    do i=1,nbas
        write(*,*) "idx_moeng_prob ",i,moeng(i),ecap_prob(i)
    end do
!------------------------------------------------------------------!
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
    init_state_idx = 12
    omega_ext = 0.0e0_wp
    eta       = 1.0e-6_wp
    cx_bold(:,:,:) = 0.0e0_wp
    call cx_pack_ab_format(one,zero,cx_bold(init_state_idx,:,:))
!------------------------------------------------------------------!
    omega_min  =  3.0e0_wp * ev_to_hartree
    omega_max  =  6.0e0_wp * ev_to_hartree
    omega_step =  0.10e0_wp * ev_to_hartree 

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
    !call core_v1
    num_omega = 100
    !call core_v2(nbas,num_omega)
    call core_v3(nbas,num_omega)
!-------------------------------------------------------------------!
end subroutine driver_1









end module mod_tdpt

