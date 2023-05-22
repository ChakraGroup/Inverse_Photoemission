module mod_dipole
use class_zipped
use mod_rand,only: rand_num
use mod_param,only: wp
implicit none
private
integer,parameter :: n3d = 3
real(wp),parameter :: PI = 3.141592653589793238460e0_wp
public :: driver_direct_int

contains

!********************************************************************!
subroutine driver_direct_int(ZP_h,ZP_e)
!********************************************************************!
implicit none
type(zip_param),intent(in) :: ZP_h
type(zip_param),intent(in) :: ZP_e

integer,parameter :: nloop  = 10
!set nuse_h to number of occupied states and nuse_e to number of virtual states! 
integer,parameter :: nuse_h = 1
integer,parameter :: nuse_e = 26

integer,parameter :: eV_to_har = 0.01837465
integer :: i1,i2,i3,ig1(n3d),ig2(n3d),iloop
real(wp) :: rv1(n3d)
real(wp) :: rv2(n3d)
real(wp) :: moval_h_ig1(ZP_h%nmo)
real(wp) :: moval_h_ig2(ZP_h%nmo)
real(wp) :: moval_e_ig1(ZP_e%nmo)
real(wp) :: moval_e_ig2(ZP_e%nmo)
integer :: imo,jmo,amo,bmo,idim
real(wp) :: npt
real(wp) :: vol3d
real(wp) :: dipole_vec(n3d)
real(wp) :: norm_h_term
real(wp) :: norm_h_sq(nuse_h)
real(wp) :: norm_sq_int_h(nuse_h)
real(wp) :: norm_int_h(nuse_h)
real(wp) :: norm_e_term
real(wp) :: norm_e_sq(nuse_e)
real(wp) :: norm_sq_int_e(nuse_e)
real(wp) :: norm_int_e(nuse_e)

real(wp) :: dip_matA(nuse_e,nuse_e,n3d)
real(wp) :: dip_matB(nuse_e,nuse_e,n3d)
real(wp) :: dip_matA_int(nuse_e,nuse_e,n3d)
real(wp) :: dip_matB_int(nuse_e,nuse_e,n3d)
real(wp) :: dip_matA_ratio(nuse_e,nuse_e,n3d)
real(wp) :: dip_matB_ratio(nuse_e,nuse_e,n3d)
real(wp) :: dipole_sq_mod(nuse_e,nuse_e)
real(wp) :: omega(nuse_e,nuse_e)
real(wp) :: moeng_h(nuse_h)
real(wp) :: moeng_e(nuse_e)
real(wp) :: coeff(nuse_e,nuse_e)
real(wp) :: osc_str(nuse_e,nuse_e)

!----------------------------------------------------------------------
    vol3d  = product(ZP_h%rmax-ZP_h%rmin)

!----------------------------------------------------------------------
    dip_matA = 0.0e0_wp
    dipole_vec = 0.0e0_wp
    dip_matB = 0.0e0_wp
!----------------------------------------------------------------------
    npt = 0.0e0_wp
    norm_h_sq = 0.0e0_wp
    norm_e_sq = 0.0e0_wp
    do i1=1,100
    do i2=1,100
    do i3=1,100
        ig1(:) = (/i1,i2,i3/)
        call grid_idx_to_rvec(ZP_h,ig1,rv1)
        !call unzip_movec_v3(ZP_h,ig1,moval_h_ig1)
        !call unzip_movec_v3(ZP_e,ig1,moval_e_ig1)
        npt = npt + 1
        do imo=1,nuse_h
            call unzip_moval_v3(ZP_h,ig1,imo,moval_h_ig1(imo))
            norm_h_term = moval_h_ig1(imo) * moval_h_ig1(imo)
            norm_h_sq(imo) = norm_h_sq(imo) + norm_h_term
        end do
        do amo=1,nuse_e
            call unzip_moval_v3(ZP_e,ig1,amo,moval_e_ig1(amo))
            norm_e_term = moval_e_ig1(amo) * moval_e_ig1(amo)
            norm_e_sq(amo) = norm_e_sq(amo) + norm_e_term
        end do
            
        do bmo=1,nuse_e
        do amo=1,nuse_e
        do idim = 1,n3d
            dipole_vec(idim) = moval_e_ig1(bmo)*moval_e_ig1(amo)*rv1(idim)
            dip_matA(bmo,amo,idim) = dip_matA(bmo,amo,idim) + dipole_vec(idim)
            !write(*,*) dip_matA(imo,amo,idim)
        end do
        end do
        end do
        
        do bmo=1,nuse_e
        do amo=1,nuse_e
        do idim = 1,n3d
            dipole_vec(idim) = moval_e_ig1(bmo)*moval_e_ig1(amo)*rv1(idim)
            dip_matB(bmo,amo,idim) = dip_matB(bmo,amo,idim) + dipole_vec(idim)
            !write(*,*) smatB(imo,amo,idim)
        end do
        end do
        end do
    end do
    end do
    end do
!---------------------------end of grid loop---------------------------!   
    norm_sq_int_h = norm_h_sq*(vol3d) / npt
    norm_sq_int_e = norm_e_sq * vol3d / npt
    !do jmo = 1,nuse_h
        !write(*,*) jmo,norm_sq_int_h(jmo)
    !end do
    !do amo = 1,nuse_e
        !write(*,*) amo,norm_sq_int_e(amo)
    !end do
    norm_int_h = sqrt(norm_sq_int_h)
    norm_int_e = sqrt(norm_sq_int_e)
    do jmo = 1,nuse_h
        norm_int_h(jmo) = sqrt(norm_sq_int_h(jmo))
        !write(*,*) jmo,norm_int_h(jmo)
    end do
    do amo = 1,nuse_e
        norm_int_e(amo) = sqrt(norm_sq_int_e(amo))
        !write(*,*) amo,norm_int_e(amo)
    end do
    !STOP "==after norm=="
    
    do bmo=1,nuse_e
    do amo=1,nuse_e
    do idim=1,n3d
    dip_matA_int(bmo,amo,idim) = dip_matA(bmo,amo,idim)*vol3d*(1.0e0_wp/real(npt,wp))
    dip_matA_ratio(bmo,amo,idim) = dip_matA_int(bmo,amo,idim)*(1.0e0_wp/((norm_int_e(amo))*(norm_int_e(amo))))
    !write(*,*)jmo,amo,dip_matA_int(jmo,amo,idim),dip_matA_ratio(jmo,amo,idim)
    end do 
    end do
    end do
    
    do bmo=1,nuse_e
    do amo=1,nuse_e
    do idim=1,n3d
    dip_matB_int(bmo,amo,idim) = dip_matB(bmo,amo,idim)*vol3d*(1.0e0_wp/real(npt,wp))
    dip_matB_ratio(bmo,amo,idim) = dip_matB_int(bmo,amo,idim)*(1.0e0_wp/((norm_int_e(amo))*(norm_int_e(amo))))
    !write(*,*)jmo,amo,dip_matB_int(jmo,amo,idim),dip_matB_ratio(jmo,amo,idim)
    end do 
    end do
    end do
    omega = 0.0e0_wp
    dipole_sq_mod = 0.0e0_wp
    open(unit=13,file='dipole_info.txt',action ='write')
!writing number of virtual orbitals used!
    write(13,*)"NMO"
    write(13,*) nuse_e
    write(13,*)"initial_idx final_idx initial_moeng final_moeng dipole_sq_mod"
    do bmo=1,nuse_e !bmo loop
    do amo=1,nuse_e !amo loop   
        do idim=1,n3d
        dipole_sq_mod(bmo,amo) = dipole_sq_mod(bmo,amo) + (dip_matA_ratio(bmo,amo,idim)*dip_matB_ratio(bmo,amo,idim))
        !write(*,*)jmo,amo,dip_matA_ratio(jmo,amo,idim)*dip_matB_ratio(jmo,amo,idim),dipole_sq_mod(jmo,amo)
        end do
        call get_moeng(ZP_h,nuse_h,moeng_h)
        call get_moeng(ZP_e,nuse_e,moeng_e)
        omega(bmo,amo) = (moeng_e(amo)-moeng_e(bmo))
!coefficient for oscillator strenght on line 224!
        coeff(bmo,amo) = (4.0e0/3.0e0)*PI*omega(bmo,amo)
!calculating oscillator strenght in case we want that info!
        osc_str(bmo,amo) = dipole_sq_mod(bmo,amo)*coeff(bmo,amo)
! we restrict to transitions with energies less than or = to 5 ev!
        if (amo >= bmo .AND. omega(amo,bmo) <= 0.183746518e0_wp) THEN
        !writing bmo and amo indices and energies and square modulus of dipole moment vectors!
        write(13,*)bmo,amo,moeng_e(bmo),moeng_e(amo),dipole_sq_mod(bmo,amo)
        end if
        !write(13,*)jmo,amo,moeng_h,moeng_e,omega,osc_str
    end do !end amo loop
    end do !end bmo loop
    close(13)


end subroutine driver_direct_int

subroutine grid_idx_to_rvec(ZP,ig,rvec)
!********************************************************************!
!here we find the position vectors corresponding to MO indices!
implicit none
type(zip_param),intent(in) :: ZP
integer,intent(in) :: ig(ZP%ndim)
real(wp),intent(out) :: rvec(ZP%ndim) 
integer :: idim
integer :: idx

    do idim=1,ZP%ndim
        idx = ig(idim)
        rvec(idim) = ZP%rmin(idim) + ((idx-1)*(ZP%rmax(idim)-ZP%rmin(idim))/(ZP%npt1d(idim)-1))
    end do

end subroutine grid_idx_to_rvec

end module mod_dipole



