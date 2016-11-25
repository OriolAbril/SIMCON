
      
      program euler
      integer :: i !,j, n=0, columns=2
      !real(8), allocatable :: y(:,:),temps(:)
      real(8), parameter :: h=0.1, t0=2, tf=12, k=5.3
      integer, parameter :: unit = 1,m=80 
      real(8) :: x(2)=[1,2]
      !real(8) :: F
      open(unit,file="proves.dat",status="replace")
      call eulere(hooke,x,t0,tf,h,k,temps,y)
      do i=1,m
         write(unit,"(I5, 5x, F10.4)") i, hooke(i*x,k)
      end do
      
      CONTAINS
      real(8) function hooke(x,k)
      Implicit none
      real(8), intent(in) :: x(2),k
      hooke=-(k*x)
      end function hooke

      subroutine eulere(x0,t0,tf,h,param,temps,y)
      real(8), dimension(:), intent(in) :: x0,param
      real(8), intent(in) :: t0,tf,h
      
      print *, fun(x0,param)
      end subroutine eulere
      end program euler
