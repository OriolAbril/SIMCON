      program harmonic
      use integradors
      real(8), parameter :: t0=2, tf=12
      real(8),parameter :: h=0.01
      real(8), dimension(2), parameter :: p=[2.0,0.2] 
      real(8),dimension(1) :: x=[0.05],v=[0.0]
      open(1,file="leapfrog0p01.dat",status="replace")
      call leapfrog(hooke,x,v,t0,tf,h,p,1)
      open(2,file="verlet0.01.dat",status="replace")
      call verlet(hooke,x,v,t0,tf,h,p,2)
      open(3,file="velocityverlet0.01.dat",status="replace")
      call velocityverlet(hooke,x,v,t0,tf,h,p,3)
      open(4,file="eulermod0.01.dat",status="replace")
      call eulermod(hooke,x,v,t0,tf,h,p,4)
      open(5,file="euler0.01.dat",status="replace")
      call eulere(hooke,x,v,t0,tf,h,p,5)
      close(1)
      close(2)
      close(3)
      close(4)
      close(5)
      CONTAINS
      function hooke(t,x,param)
      Implicit none
      real(8), dimension(:), intent(in) :: param
      real(8), dimension(:), intent(in) :: x
      real(8), intent(in) :: t
      real(8), dimension(size(x)) :: hooke
      hooke=x*t
      hooke=-(param(1)*x)/param(2)
      end function hooke

      end program harmonic
