C23456789
      subroutine init_rand_seed() 
      implicit none
      integer*4 i
      call system_clock(i)
      i = i + 632
      call srand(i)
      end

      subroutine int_rand_num(imin,imax,i)
      implicit none
      integer imin,imax,i
      i = imin + nint((imax-imin)*rand())
      end 

