      program montecarlo
      integer :: i,k, n=0
      real(8) :: meane,desveste,en,press,meanp,desvestp
      meane=0.0d0
      meanp=0.0d0
      desveste=0.0d0
      desvestp=0.0d0
      open(1,file="energy-press.dat",status="old")
      open(2,file='stats.dat',status='new')
      do
         read(1,*,end=1) i,en,press
         n = n+1
         meane=meane+en
         meanp=meanp+press
      end do
 1    rewind(1)
      meanp=meanp/(n*1.d0)
      meane=meane/(n*1.d0)
      do k=1,n
         read(1,*) i,en,press
         desveste=desveste+(en-meane)**2
         desvestp=desvestp+(press-meanp)**2
      end do
      desvestp=sqrt(desvestp/(n*1.d0))
      desveste=sqrt(desveste/(n*1.d0))
      write(2,*) n,mitja, desvest
      

      end program montecarlo
