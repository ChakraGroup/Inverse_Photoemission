!********************************************************************!
program prog_1
!********************************************************************!
use mod_param,only: wp,check_for_underflow
use mod_rand,only: init_rand_with_sysclock
use mod_tdpt,only: driver_1
implicit none

real(wp) :: tstart,tend
character(100) :: start_date
character(100) :: end_date
!--------------------------------------------------------------------!
     call fdate(start_date)
     write(*,*) "------------------------------------------------------"
     write(*,*) "START_DATE: ",trim(adjustL(start_date))
     write(*,*) "------------------------------------------------------"
     write(*,*)
     call cpu_time(tstart)
!--------------------------------------------------------------------!
    call init_rand_with_sysclock()
!--------------------------------------------------------------------!
    call driver_1
!--------------------------------------------------------------------!
    call cpu_time(tend)
    write(*,*) 
    write(*,*) "tot_calc_time= ",tend-tstart
    call fdate(end_date)
    write(*,*) "------------------------------------------------------"
    write(*,*) "END_DATE: ",trim(adjustL(end_date))
    write(*,*) "------------------------------------------------------"
!--------------------------------------------------------------------!
end program prog_1
