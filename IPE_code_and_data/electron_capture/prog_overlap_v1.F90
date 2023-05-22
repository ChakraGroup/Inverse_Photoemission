!********************************************************************!
program direct
!********************************************************************!
use mod_param,only: wp,check_for_underflow
use mod_rand,only: init_rand_with_sysclock
use class_zipped,only: zip_param,            & 
                       read_zipped_grid_v2,  &
                       check_norm
use mod_direct,only:get_unit_vec,driver_direct_int
implicit none
character(100),parameter :: file_zipped_hole = "zipped_grid_moh.ufm" 
character(100),parameter :: file_zipped_elec = "zipped_grid_moe.ufm"
type(zip_param) :: ZP_h
type(zip_param) :: ZP_e
character(*),parameter :: infile = "unit_vec.txt"
integer :: nunit !number_of_unit_vec
integer :: ndim  !dim of space
real(wp),allocatable :: unit_vec(:,:) !nunit x ndim 
    call init_rand_with_sysclock()

    call read_zipped_grid_v2(file_zipped_hole,ZP_h)
    !call check_norm(ZP_h)

    call read_zipped_grid_v2(file_zipped_elec,ZP_e)     
    !call check_norm(ZP_e)
    
    call get_unit_vec(infile,nunit,ndim,unit_vec)
    !write(*,*) "--in_main_code--"
    !write(*,*) "unit_vec: ",unit_vec(:,:)

    !if(allocated(unit_vec)) deallocate(unit_vec)

    call driver_direct_int(ZP_h,ZP_e)

end program direct


