*****************************************************************************
*          Molecular Dynamics code to simulate at NVE collectivity
*          a system of N atoms of Ar inside a cubic box, interacting
*          through a truncated Lennard-Jones force field at 2.5 sigma.
*          Leap-frog Verlet integration algorithm.
*
*          Students must consider several algorithms of numerical
*          integration and compare their efficiencies
*****************************************************************************
      PROGRAM leapfroglj
      use modfrog
c     1. Defining variables and dimensions
      implicit none
      real(8), dimension(3,1000) :: r,vinf,accel
      integer, parameter :: nhis=500,it0=50,t0max=200
      integer :: nconf,natoms,nf,is,i,l
      real(8) :: mass, sigma,epsil,deltat,rc,pi,boxlength,vol,ro,uvel
      real(8) :: delg,epot,press,ecin,temp,etot,pressio,rc3,ptail
      real(8) :: tpress,ptot
      real(8), allocatable :: g(:),vacf(:),r2t(:)
      integer, allocatable :: ntime(:)
      real(8), dimension(1000,3,t0max) :: x0,vx0
      integer :: t0
      integer, dimension(t0max) :: time0



c     2. Reading data and computing related quantities

      open(1,file='leap-lj.data',status='old')
      read(1,*) nconf
      read(1,*) natoms
      read(1,*) mass
      read(1,*) sigma,epsil
      read(1,*) deltat
      close(1)

      allocate(vacf(nconf))
      allocate(r2t(nconf))
      allocate(ntime(nconf))
      nf = 3*natoms-3           !number of degrees of freedom
      rc = 2.5d0                !sigma = range of the potential in reduced units
      pi = 4*datan(1.d0)
c     3. Reading initial configuration (positions, velocities) in A and A/ps

      open(2,file='leap-conf.data',status='old')
      do is = 1,natoms
         read(2,*) (r(l,is),l=1,3)
         read(2,*) (vinf(l,is),l=1,3)
      end do
      read(2,*) boxlength
      close(2)

c     Opening files to write results

      open(3,file='energy-leap.dat',status='replace')
      open(4,file='temp-leap.dat',status='replace')
      open(2,file='press-leap.dat',status='replace')

c     4. Change to reduced units
      print *,natoms/boxlength**3
      call reduced(natoms,r,vinf,boxlength,deltat,epsil,sigma,
     & mass,uvel)
      vol=boxlength**3
      ro=natoms/vol
      print *,ro
!     Initialize g(r)
      allocate(g(nhis))
      call grinit(nhis,g,boxlength,delg)
c     5. Start the loop to generate new configurations
      do i = 1,nconf
         call forces(natoms,r,boxlength,accel,rc,epot,press,g,delg,
     &        nhis,boxlength)
         call velpos(natoms,vinf,accel,deltat,r,nf,ecin,temp,
     &        boxlength,i,nconf,it0,t0,t0max,time0,x0,vx0,vacf,
     &        r2t,ntime)         
         etot = ecin + epot
         pressio = press/(3.0*vol)
         rc3=rc**3
         ptail = (16.d0/3.d0)*ro**2*pi*(2.d0/3.d0/rc3**3-1.0/rc3)
         pressio = pressio + ptail
         !print *, pressio
         tpress = ro * temp
         ptot = pressio + tpress
!         write(*,*) temp,etot
         write(3,*) i*deltat, etot
         write(4,*) i*deltat, temp
         write(2,*) i*deltat, ptot, pressio, tpress
      end do
      close(2)
      close(3)
      close(4)

c     6. Saving last configuration in A and A/ps

      open(5,file='gdr.dat',status='replace')
      call grend(5,nhis,delg,g,natoms,ro,nconf)

      open(12,file='vacf-r2t.dat')
      do is=1,nconf
         write(12,*) is*deltat,vacf(is),r2t(is)
      end do
      close(12)
      open(11,file='newconf.data',status='unknown')
      do is = 1,natoms
         write(11,*) (r(l,is)*sigma,l=1,3)
         write(11,*) (vinf(l,is)*uvel,l=1,3)
      end do
      write(11,*) boxlength*sigma
      close(11)
      
      end program leapfroglj
