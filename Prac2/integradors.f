!     Modul amb subrutines per integrar trajectories a partir de particules
      module integradors
      contains
!     Per cridar-les cal haver obert el fitxer i passar la unit
!     Els parametres son:
!     fun: la funcio que retorna l'acceleracio. la funcio ha de tenir 3 inputs, temps (real), x i param que han de ser arrays
!     (per molt que tinguin dimensio 1)
!     x0,v0,param: array de condicions inicials (posicio i velocitat) i de parametres de la funcio
!     t0,tf,h: temps inicial, final i pas de temps respectivament (reals)
!     unit: unitat on s'ha obert el fitxer (integer)

!     Euler explicit
      subroutine eulere(fun,x0,v0,t0,tf,h,param,unit)
      integer, intent(in) :: unit
      real(8), dimension(:), intent(in) :: x0,v0,param
      real(8), intent(in) :: t0,tf,h
      real(8), dimension(size(x0)) :: x,v,a
      real(8) :: t
      integer :: N,i
      interface
         function fun(t,x,p)
         real(8), dimension(:), intent(in) :: x
         real(8), dimension(:), intent(in) :: p
         real(8), intent(in) :: t
         real(8), dimension(size(x)) :: fun
         end function
      end interface
      t=t0
      x=x0
      v=v0
      N=idnint((tf-t0)/h-0.9)
      do i=1,N
         write(unit,*) t,x,v
         a=fun(t,x,param)
         x=x+h*v+h**2/2.0*a
         v=v+h*a
         t=t+h
      end do
      write(unit,*) t,x,v
      end subroutine eulere

!     Euler modificat
      subroutine eulermod(fun,x0,v0,t0,tf,h,param,unit)
      integer, intent(in) :: unit
      real(8), dimension(:), intent(in) :: x0,v0,param
      real(8), intent(in) :: t0,tf,h
      real(8), dimension(size(x0)) :: x,x2,v,a,a1
      real(8) :: t
      integer :: N,i
      interface
         function fun(t,x,p)
         real(8), dimension(:), intent(in) :: x
         real(8), dimension(:), intent(in) :: p
         real(8), intent(in) :: t
         real(8), dimension(size(x)) :: fun
         end function
      end interface
      t=t0
      x=x0
      v=v0
      N=idnint((tf-t0)/h-0.9)
      do i=1,N
         write(unit,*) t,x,v
         a=fun(t,x,param)
         t=t+h
         x2=x+h*v+(h**2)/2.0*a
         a1=(a+fun(t,x2,param))/2.0
         x=x+h*v+h**2/2.0*a1
         v=v+h*a1
      end do
      write(unit,*) t,x,v
      end subroutine eulermod

!     Original Verlet
      subroutine verlet(fun,x0,v0,t0,tf,h,param,unit)
      integer, intent(in) :: unit
      real(8), dimension(:), intent(in) :: x0,v0,param
      real(8), intent(in) :: t0,tf,h
      real(8), dimension(size(x0)) :: x,v,a,x2,xaux
      real(8) :: t
      integer :: N,i
      interface
         function fun(t,x,p)
         real(8), dimension(:), intent(in) :: x
         real(8), dimension(:), intent(in) :: p
         real(8), intent(in) :: t
         real(8), dimension(size(x)) :: fun
         end function
      end interface
      t=t0
      x=x0
      xaux=x0
      v=v0
      write(unit,*) t,x,v
      a=fun(t,x,param)
      t=t+h
      x2=x+h*v+(h**2)/2.0*a
      N=idnint((tf-t0)/h-0.9)
      do i=1,N
         a=fun(t,x2,param)
         xaux=2*x2-x+(h**2)*a
         v=(xaux-x)/(2*h)
         x=x2
         x2=xaux
         write(unit,*) t,x,v
         t=t+h
      end do
      end subroutine verlet

!     Velocity Verlet
      subroutine velocityverlet(fun,x0,v0,t0,tf,h,param,unit)
      integer, intent(in) :: unit
      real(8), dimension(:), intent(in) :: x0,v0,param
      real(8), intent(in) :: t0,tf,h
      real(8), dimension(size(x0)) :: x,v,a,a2
      real(8) :: t
      integer :: N,i
      interface
         function fun(t,x,p)
         real(8), dimension(:), intent(in) :: x
         real(8), dimension(:), intent(in) :: p
         real(8), intent(in) :: t
         real(8), dimension(size(x)) :: fun
         end function
      end interface
      t=t0
      x=x0
      v=v0
      a=fun(t,x,param)
      write(unit,*) t,x,v
      N=idnint((tf-t0)/h-0.9)
      do i=1,N
         x=x+v*h+a*h**2/2.0
         t=t+h
         a2=fun(t,x,param)
         v=v+(a+a2)/2.0*h
         a=a2
         write(unit,*) t,x,v
      end do
      end subroutine velocityverlet

!     Leapfrog Verlet
      subroutine leapfrog(fun,x0,v0,t0,tf,h,param,unit)
      integer, intent(in) :: unit
      real(8), dimension(:), intent(in) :: x0,v0,param
      real(8), intent(in) :: t0,tf,h
      real(8), dimension(size(x0)) :: x,v,a,x2,v12,v22
      real(8) :: t
      integer :: N,i
      interface
         function fun(t,x,p)
         real(8), dimension(:), intent(in) :: x
         real(8), dimension(:), intent(in) :: p
         real(8), intent(in) :: t
         real(8), dimension(size(x)) :: fun
         end function
      end interface
      t=t0!0
      x=x0!0
      v=v0!0
      write(unit,*) t,x,v
      a=fun(t,x,param)!0
      t=t+h!1
      x2=x+h*v+(h**2)/2.0*a!1
      v12=(x2-x)/h!0.5
      N=idnint((tf-t0)/h-0.9)
      do i=1,N
         a=fun(t,x2,param)!2
         v22=v12+a*h!1.5
         x=x2!1
         x2=x2+v22*h!2
         v=(v22+v12)/(2)!1
         v12=v22!1.5
         write(unit,*) t,x,v
         t=t+h!2
      end do
      end subroutine leapfrog
      
      end module integradors
