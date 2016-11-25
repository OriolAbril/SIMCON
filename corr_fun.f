      program corrfun
      integer :: i,j, n=0, columns=2
      real*8, allocatable :: data(:,:),S(:),w(:)
      real*8, parameter :: dt=0.005, dw=0.05
      integer, parameter :: unit = 1,m=8000 
      n = 0
      allocate(S(m))
      allocate(w(m))
      open(unit,file="corr.dat",status="old",action="read")
      do
         read(unit,*,end=1)
         n = n+1
      end do
 1    rewind(unit)
      allocate(data(n,columns))
      do i = 1, n
         read(unit,*) data(i,:)
!         print *, data(i,:) 
      end do
      close(unit)

      open(unit,file="integral.dat",status="replace")
      do i=1,m
         w(i)=(i-1)*dw
         S(i)=data(1,2)*dcos(data(1,1)*w(i))+data(n,2)*dcos(data(n,1)
     &*w(i))
         do j=2,n-1
            if (mod(j,2)==0) THEN
               S(i)=S(i)+2*data(j,2)*dcos(data(j,1)*w(i))
            else
               S(i)=S(i)+4*data(j,2)*dcos(data(j,1)*w(i))
            end if
         end do
         S(i)=(dt/3.0)*S(i)
         write(unit,"(F10.6,5x,ES18.8E3)") w(i),S(i)/S(1)
      end do
      

      end program corrfun
