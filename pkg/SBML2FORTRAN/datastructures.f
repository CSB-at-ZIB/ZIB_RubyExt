      program datastructures
c
      implicit none
c
      integer        j, n
c
      character      a(9)*11
c
c
      n = 3
c
      call get_str_array( a )
c
c
      write( * , 100 ) ( a(j), j = 1,n )
c
100   format( 5('"',A,'", ') ) 
c
      end program datastructures
