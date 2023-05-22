module mod_direct
use class_zipped
use mod_rand,only: rand_num
use mod_param,only: wp
implicit none
private
integer,parameter :: n3d = 3
real(wp),parameter :: PI = 3.141592653589793238460e0_wp
public :: get_unit_vec
public :: driver_direct_int

contains

subroutine get_unit_vec(infile,nunit,ndim,unit_vec)
implicit none
integer,parameter :: wp = 8
character(*),intent(in) :: infile
integer,intent(out) :: nunit                        !number_of_unit_vec
integer,intent(out) :: ndim                         !dim of space
real(wp),allocatable,intent(inout) :: unit_vec(:,:) !nunit x ndim 

integer :: i
real(wp) :: d
    open(unit=100,file=infile,action='read')
        read(100,*) !title
        read(100,*) nunit,ndim
        !write(*,*) "num_unit_vec_ndim= ",nunit,ndim
        if(nunit <= 0) then
            STOP "***I am not going to allocate an empty array***"
        end if
        if(allocated(unit_vec)) deallocate(unit_vec)
        allocate(unit_vec(nunit,ndim))
        unit_vec = 0.0e0_wp
        !write(*,*) "---list_of_normalized_unit_vectors--"
        do i=1,nunit
            read(100,*) unit_vec(i,:)
            d = dot_product(unit_vec(i,:),unit_vec(i,:))
            if( d > 0.0e0_wp) then
                unit_vec(i,:) = unit_vec(i,:)/sqrt(d)
            else
                STOP "***I am not going to use ZERO vector***"
            end if
            !write(*,*) unit_vec(i,:)
        end do
    close(100)
end subroutine get_unit_vec

!********************************************************************!
subroutine driver_direct_int(ZP_h,ZP_e)
!********************************************************************!
implicit none
type(zip_param),intent(in) :: ZP_h
type(zip_param),intent(in) :: ZP_e
integer,parameter :: nloop  = 10
integer,parameter :: nuse_h = 216
integer,parameter :: nuse_e = 456
integer,parameter :: homo_idx = 1
integer,parameter :: lumo_idx = 1
integer,parameter :: ione = 1 
integer,parameter :: ihomo = 1 
integer,parameter :: alumo = 1 
integer,parameter :: ihun = 100 
real(wp),parameter :: one = 1.0e0_wp
real(wp),parameter :: gamma = 0.360e0_wp
real(wp),parameter :: mass_e = 9.1093837015e-31_wp
real(wp),parameter :: h_bar = 1.054571817e-34_wp
real(wp),parameter :: charge_e = -1.602176634e-19_wp
integer,parameter :: nr = 4
integer,parameter :: eV_to_har = 0.01837465
integer :: i1,i2,i3,ig1(n3d),ig2(n3d),iloop,i
real(wp) :: rv1(n3d)
real(wp) :: rv2(n3d)
real(wp) :: moval_h_ig1(ZP_h%nmo)
real(wp) :: moval_h_ig2(ZP_h%nmo)
real(wp) :: moval_e_ig1(ZP_e%nmo)
real(wp) :: moval_e_ig2(ZP_e%nmo)
integer :: imo,jmo,amo,bmo,idim,kmo,cmo,ndim
real(wp) :: npt
real(wp) :: ov_h(nuse_h)
real(wp) :: ov_e(nuse_h)
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
real(wp) :: energy_j(nuse_h)
real(wp) :: energy_a(nuse_e)
real(wp) :: rmin
real(wp) :: rmax

real(wp) :: smatA(nuse_e,n3d)
real(wp) :: smatB(nuse_e,n3d)
real(wp) :: smatA_int(nuse_e,n3d)
real(wp) :: smatB_int(nuse_e,n3d)
real(wp) :: smatA_ratio(nuse_e,n3d)
real(wp) :: smatB_ratio(nuse_e,n3d)
real(wp) :: dipole_sq_mod(nuse_e)
real(wp) :: omega(nuse_e)
real(wp) :: constant
real(wp) :: moeng_h(nuse_h)
real(wp) :: moeng_e(nuse_e)
real(wp) :: moeng(nuse_e)
real(wp) :: coeff(nuse_e)
real(wp) :: osc_str(nuse_e)
integer :: iunit                        !number_of_unit_vec
integer :: nunit                        !number_of_unit_vec
integer :: nmo_e                  !number_of_unit_vec


character(*),parameter :: infile = 'unit_vec.txt'
real(wp),allocatable :: unit_vec(:,:) !nunit x ndim 
real(wp),allocatable :: smat_e_sin(:,:) !nuse_e x nunit 
real(wp),allocatable :: smat_e_cos(:,:) !nuse_e x nunit 
real(wp) :: k_vec(n3d) 
real(wp),allocatable :: Re_ratio(:,:) !nuse_e x nunit 
real(wp),allocatable :: Im_ratio(:,:) !nuse_e x nunit 
real(wp),allocatable :: sq_mod_overlap(:,:)
                    


real(wp) :: r12
real(wp) :: v12
real(wp) :: norm_tot

real(wp) :: Ediff(nuse_h)
real(wp) :: Ediff_e(nuse_e)
real(wp) :: theta
real(wp) :: st
real(wp) :: k_mag(nuse_e)
real(wp) :: ct
real(wp) :: moe_st(nuse_e)
real(wp) :: moe_ct(nuse_e)
real(wp) :: moeng_homo(homo_idx)
real(wp) :: moe_energy(nuse_e)

!----------------------------------------------------------------------


call get_unit_vec(infile,nunit,ndim,unit_vec)

    if(allocated(smat_e_sin)) deallocate(smat_e_sin)
    allocate(smat_e_sin(nuse_e,nunit))
    
    if(allocated(smat_e_cos)) deallocate(smat_e_cos)
    allocate(smat_e_cos(nuse_e,nunit))
    
    
    if(allocated(Re_ratio)) deallocate(Re_ratio)
    allocate(Re_ratio(nuse_e,nunit))
    
    if(allocated(Im_ratio)) deallocate(Im_ratio)
    allocate(Im_ratio(nuse_e,nunit))
    
    if(allocated(sq_mod_overlap)) deallocate(sq_mod_overlap)
    allocate(sq_mod_overlap(nuse_e,nunit))

!---------------------------------------------------------------------
    vol3d  = product(ZP_h%rmax-ZP_h%rmin)

    call get_moeng(ZP_e,nuse_e,moeng_e)
    
    !call get_moeng(ZP_e,nuse_e,moeng_e)
    call get_moeng(ZP_h,nuse_h,moeng_h)
    
   ! write(*,*) "getting moeng_e_array"
   ! write(*,*) moeng_e
   ! write(*,*) "getting moeng_h_array"
   ! write(*,*) moeng_h(1)
    !call get_moeng(ZP_e,nuse_e,moeng_e)

!----------------------------------------------------------------------
!----------------------------------------------------------------------
        npt = 0.0e0_wp
        norm_h_sq = 0.0e0_wp
        smat_e_cos = 0.0e0_wp
        smat_e_sin = 0.0e0_wp
        norm_h_sq = 0.0e0_wp
        norm_e_sq = 0.0e0_wp
        moe_st = 0.0e0_wp
        moe_ct = 0.0e0_wp
        do i1=1,100
        do i2=1,100
        do i3=1,100
            npt = npt + 1.0e0_wp
            ig1(:) = (/i1,i2,i3/)
            call grid_idx_to_rvec(ZP_h,ig1,rv1)
            do amo=1,nuse_e
                call unzip_moval_v3(ZP_e,ig1,amo,moval_e_ig1(amo))
                k_mag(amo) = sqrt(2*abs(moeng_e(amo)-moeng_h(homo_idx)))
                do iunit=1,nunit
                    k_vec(:) = unit_vec(iunit,:)*k_mag(amo)
                    k_vec(:) = unit_vec(iunit,:)
                    theta= dot_product(k_vec(:),rv1(:))
                    !write(*,*) theta
                    ct = cos(theta)
                    st = sin(theta)
                        norm_e_term = moval_e_ig1(amo) * moval_e_ig1(amo)
                        norm_e_sq(amo) = norm_e_sq(amo) + norm_e_term
                        moe_st(amo) =  moval_e_ig1(amo) * st
                        !write(*,*) moe_st(amo)
                        smat_e_sin(amo,iunit) = smat_e_sin(amo,iunit) + moe_st(amo)
                        moe_ct(amo) = moval_e_ig1(amo) * ct
                        smat_e_cos(amo,iunit) = smat_e_cos(amo,iunit) + moe_ct(amo) 
                end do
            end do
        end do
        end do
        end do
        !do amo=1,nuse_e
         !   write(*,*) amo,moeng_e(amo),moeng_h(1),k_mag(amo)
        !end do
        
        open(unit=11,file='overlap_info.txt',action ='write')
        write(11,*) "amo,iunit,moeng,k_mag,Im_int,Re_int,sq_mod"
        write(11,*) nuse_e,nunit
        do amo=1,nuse_e
            do iunit=1,nunit
               Im_ratio(amo,iunit) = smat_e_sin(amo,iunit)/norm_e_sq(amo)
               Re_ratio(amo,iunit) = smat_e_cos(amo,iunit)/norm_e_sq(amo)
               sq_mod_overlap(amo,iunit) = ((Im_ratio(amo,iunit))**2) + ((Re_ratio(amo,iunit))**2)
               write(11,*) amo,iunit,moeng_e(amo),k_mag(amo),Im_ratio(amo,iunit),Re_ratio(amo,iunit),sq_mod_overlap(amo,iunit)
            end do
        end do
        close(11)
        
        
        !do amo=1,nuse_e
         !   call unzip_moval_v3(ZP_e,ig1,amo,moval_e_ig1(amo))
          !  do iunit=1,nunit    
           !     call get_unit_vec(infile,nunit,ndim,unit_vec)
            !    write(*,*)"iunit,amo,Im,Re",iunit,amo,smat_e_sin(amo,iunit),smat_e_cos(amo,iunit)
           ! end do
        !end do
        


!---------------------------end of grid loop---------------------------!   

    
    if(allocated(smat_e_sin)) deallocate(smat_e_sin)
    if(allocated(smat_e_cos)) deallocate(smat_e_cos)
end subroutine driver_direct_int

subroutine grid_idx_to_rvec(ZP,ig,rvec)
!********************************************************************!
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

end module mod_direct


