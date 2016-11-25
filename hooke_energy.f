      program hookeenergy
      real(8), parameter :: k=2.0, m=0.2 
      integer :: i
      open(1,file="leapfrog0p01.dat",status="old")
      open(2,file="verlet0.01.dat",status="old")
      open(3,file="velocityverlet0.01.dat",status="old")
      open(4,file="eulermod0.01.dat",status="old")
      open(5,file="euler0.01.dat",status="old")
      open(11,file="en_leapfrog0p01.dat",status="replace")
      open(12,file="en_verlet0.01.dat",status="replace")
      open(13,file="en_velocityverlet0.01.dat",status="replace")
      open(14,file="en_eulermod0.01.dat",status="replace")
      open(15,file="en_euler0.01.dat",status="replace")
      do i=1,5
         call energy(i,10+i,k,m)
         close(i)
         close(i+10)
      end do

      CONTAINS
      subroutine energy(inunit,outunit,k,m)
      real(8) :: llegit(3),escrit(4)
      integer, intent(IN) :: inunit,outunit
      real(8), intent(in) :: k,m
      do
         read(inunit,*,end=1) llegit(:)
         escrit(4)=0.5*(m*llegit(3)**2+k*(llegit(2)**2))
         escrit(1:3)=llegit(:)
         write(outunit,*) escrit(:)
      end do
 1    rewind(inunit)
      end subroutine energy

      function hooke(t,x,param)
      Implicit none
      real(8), dimension(:), intent(in) :: param
      real(8), dimension(:), intent(in) :: x
      real(8), intent(in) :: t
      real(8), dimension(size(x)) :: hooke
      hooke=x*t
      hooke=-(param(1)*x)/param(2)
      end function hooke

      end program hookeenergy
