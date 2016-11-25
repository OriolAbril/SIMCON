*****************************************************************************
*          Molecular Dynamics code to simulate at NVE collectivity
*          a system of N atoms of Ar inside a cubic box, interacting
*          through a truncated Lennard-Jones force field at 2.5 sigma.
*          Leap-frog Verlet integration algorithm.
*
*          Students must consider several algorithms of numerical
*          integration and compare their efficiencies
*****************************************************************************
      PROGRAM montelj
      use utilsmonte
c     1. Defining variables and dimensions
      implicit none
      real(8), allocatable :: r(:,:)
      integer, parameter :: nhis=500
      integer :: ncycle,natoms,nf,is,i,l,nconf=0,nsample,move=0
      real(8) :: sigma,epsil,delta,rc,pi,boxlength,vol,ro,epot,etail
      real(8) :: delg,press,pressio,ptail,energy,beta,temp
      real(8), dimension(nhis) :: g



c     2. Reading data and computing related quantities

      open(1,file='lj-monte.data',status='old')
      read(1,*) ncycle,nsample
      read(1,*) natoms
      read(1,*) temp
      read(1,*) sigma,epsil
      read(1,*) delta
      read(1,*) ro
      close(1)
      print *,ro
      nf = 3*natoms-3           !number of degrees of freedom
      rc = 4.d0                !sigma = range of the potential in reduced units
      pi = 4*datan(1.d0)
      beta=1.0d0/temp
c     3. Reading initial configuration (positions, velocities) in A and A/ps
      allocate(r(3,natoms))
      open(2,file='conf-monte.data',status='old')
      do is = 1,natoms
         read(2,*) (r(l,is),l=1,3)
      end do
      read(2,*) boxlength
      close(2)
      
c     Initiate random seed
      call init_random_seed()
      
c     Opening files to write results

      open(3,file='energy-press.dat',status='replace')

c     4. Change to reduced units
      call reduced(natoms,r,boxlength,delta,sigma)
      vol=boxlength**3
      ro=natoms/vol
      print *,ro
      etail=8.0d0*pi/3.0d0*natoms*ro*(1.0d0/(3.0d0*rc**9)-1.0d0/rc**3)
      ptail=16.0d0*pi**2/3.0d0*ro**2*(3.0d0/(3.0d0*rc**9)-1.0d0/rc**3)
      call potenergy(r,energy,etail,natoms,rc,boxlength)

!     Initialize g(r)
      call grinit(nhis,g,boxlength,delg)

c     5. Start the loop to generate new configurations
      do i = 1,ncycle
         call mcmove(natoms,r,boxlength,rc,delta,beta,energy,etail,move)
         if (mod(i,nsample).eq.0) then
            call mcsample(natoms,r,boxlength,rc,epot,press,g,
     &           delg,nhis)
            epot=epot+etail
            pressio = press/(3.0*vol)
            pressio = pressio + ptail
            nconf=nconf+1
            print *,i
            print *,epot
            print *,move
            write(3,*) i,epot, pressio
         end if
      end do
      close(3)
      write(*,*) move*1.0d0/(ncycle*1.0d0)


c     6. Saving last configuration in A and A/ps

      open(5,file='gdr.dat',status='replace')
      call grend(5,nhis,delg,g,natoms,ro,nconf)
      open(11,file='newconf_mc.data',status='unknown')
      
      do is = 1,natoms
         write(11,*) (r(l,is)*sigma,l=1,3)
      end do
      write(11,*) boxlength*sigma
      close(11)
      
      end program montelj
