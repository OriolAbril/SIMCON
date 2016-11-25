      subroutine init_random_seed()
      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))

      CALL SYSTEM_CLOCK(COUNT=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)

      DEALLOCATE(seed)
      end
      
      program montecarlo
      integer :: i,j,k, n=10 ,Nhit
      real*8, parameter :: Aexacta=13.5916
      integer, parameter :: m=100, l=6
      real*8 :: ran(2),ran2(2),integral,funcio,mitja,desvest
      
      open(1,file="integrals_monte.dat",status="replace")
      open(2,file="estadistica_monte.dat",status="replace")
      call init_random_seed()
      do  k=1,l
         n=n*10
         write(1,*) n
         mitja=0
         desvest=0
         do i=1,m
            Nhit=0
            do j=1,n
               call  random_number(ran)
               ran2(2)=ran(2)*11
               ran2(1)=ran(1)
               funcio=1.0/ran2(1)
               if (funcio>ran2(2)) then
                  Nhit=Nhit+1
               end if
            end do
            integral=44.0*Nhit/n
            mitja=mitja+integral
            desvest=desvest+(integral-Aexacta)**2
            write(1,*) integral
!            print *, integral
         end do
         mitja=mitja/m
         desvest=sqrt(desvest/m)
         write(2,"(I10,5x,F10.6,5x,F10.6)") n,mitja, desvest
      end do
      

      end program montecarlo
