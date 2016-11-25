      module modfrog
      
      CONTAINS
*********************************************************
*********************************************************
c              subroutine reduced
*********************************************************
*********************************************************

      subroutine reduced(natoms,r,vinf,boxlength,deltat,epsil,
     & sigma,mass,uvel)
      implicit double precision(a-h,o-z)
      double precision mass
      dimension r(3,1000),vinf(3,1000)

      rgas = 8.314472673d0      !J/(mol*K)
      utime = sigma*dsqrt(mass/epsil)*dsqrt(10.d0/rgas)
      uvel = sigma/utime        !unit of velocity, expressed in A/ps

      boxlength = boxlength/sigma
      deltat = deltat/utime
      do is = 1,natoms
         do l = 1,3
            r(l,is) = r(l,is)/sigma
            vinf(l,is) = vinf(l,is)/uvel
         end do
      end do

      end subroutine reduced

*********************************************************
*********************************************************
c              subroutine forces
*********************************************************
*********************************************************

      subroutine forces(natoms,r,boxlength,accel,rc,epot,press,g,
     &     delg,nhis,gbox)
      implicit none
      real(8), dimension(:) :: g
      real(8), dimension(3,1000) :: r,accel
      integer, intent(in) :: nhis,natoms
      real(8), intent(in) :: boxlength,rc,delg,gbox
      real(8), intent(inout) :: press,epot
      integer :: is,l,js
      real(8) :: pot,p

      do is = 1,natoms
         do l = 1,3
            accel(l,is) = 0.d0  !sets accelerations to 0
         end do
      end do
      epot = 0.d0
      press = 0.d0

c      atom-atom interactions

      do is = 1,natoms-1
         do js = is+1,natoms
            call lj(is,js,r,boxlength,accel,rc,pot,p,g,delg,nhis,gbox)
            epot = epot + pot
            press = press + p
         end do
      end do


      end subroutine forces

*********************************************************
*********************************************************
c              subroutine Lennard-Jones
*********************************************************
*********************************************************

      subroutine lj(is,js,r,boxlength,accel,rc,pot,p,g,delg,nhis,gbox)
      implicit none
      real(8), dimension(3,1000) :: r,accel
      real(8), dimension(:), intent(inout) :: g
      real(8), dimension(3) :: rij
      integer, intent(in) :: nhis,is,js
      real(8),intent(in) :: boxlength,rc,delg,gbox
      real(8), intent(out) :: pot,p
      real(8) :: rr2,rijl,rr,ynvrr2,ynvrr6,ynvrr12,forcedist
      integer :: l
      rr2 = 0.d0
      p = 0.d0
      pot = 0.d0
      do l = 1,3
         rijl = r(l,js) - r(l,is)
         rij(l) = rijl - boxlength*dnint(rijl/boxlength)
         rr2 = rr2 + rij(l)*rij(l)
      end do

      rr = dsqrt(rr2)
      call griterforce(rr,gbox,g,delg,nhis)
      if (rr.lt.rc) then
         ynvrr2 = 1.d0/rr2
         ynvrr6 = ynvrr2*ynvrr2*ynvrr2
         ynvrr12 = ynvrr6*ynvrr6
         forcedist = 24.d0*(2.d0*ynvrr12-ynvrr6)*ynvrr2
         pot = 4.d0*(ynvrr12-ynvrr6)
         p=rr2*forcedist
         do l = 1,3
            accel(l,is) = accel(l,is) - forcedist*rij(l)
            accel(l,js) = accel(l,js) + forcedist*rij(l)
         end do
      end if


      end subroutine lj

*********************************************************
*********************************************************
c              subroutine velpos
*********************************************************
*********************************************************

c       Calculating velocity at instants t+delta/2 and t
c       and getting temporal position at instant t + deltat.

      subroutine velpos(natoms,vinf,accel,deltat,r,nf,ecin,temp,
     & boxlength,i,nconf,it0,t0,t0max,time0,x0,vx0,vacf,r2t,ntime)
      implicit none
      real(8), dimension(3,1000), intent(inout) :: r,vinf
      real(8), dimension(3,1000), intent(in) :: accel
      integer, intent(in) :: it0,t0max,nconf,natoms,i,nf
      real(8), intent(in) :: deltat,boxlength
      real(8),intent(inout) :: ecin,temp
      real(8), dimension(nconf),intent(inout) :: vacf(:),r2t(:)
      real(8), dimension(1000,3,t0max) :: x0,vx0
      real(8), dimension(3,1000) :: rsup
      integer, intent(inout) :: t0
      integer, intent(inout) :: ntime(nconf),time0(t0max)
      integer :: is,l
      real(8) :: v,v2,vsup
      
      ecin = 0.d0
      do is = 1,natoms
         v2 = 0.d0
         do l = 1,3
            vsup = vinf(l,is) + accel(l,is)*deltat
            rsup(l,is) = r(l,is) + vsup*deltat
            v = (vsup+vinf(l,is))/2.d0
            v2 = v2 + v*v
            vinf(l,is) = vsup
         end do
         ecin = ecin + 0.5d0*v2
      end do
      call diffusion(i,natoms,nconf,it0,t0,t0max,time0,x0,vx0,r,
     & vinf,vacf,r2t,ntime)
      temp = 2.d0*ecin/dfloat(nf)
      
ccccc       Applying  periodic boundary conditions
      r=rsup
      Do is=1,natoms
         do l=1,3
            if (r(l,is).lt.0) r(l,is) = r(l,is) + boxlength
            if (r(l,is).gt.boxlength) r(l,is) = r(l,is) - boxlength
         end do
      end do
      end subroutine velpos

*********************************
!     Radial distribution function
*********************************
      
      subroutine grinit(nhis,g,gbox,delg)
      integer, intent(in) :: nhis
      real(8), dimension(nhis) :: g
      real(8), intent(in) :: gbox
      real(8),intent(out) :: delg
      integer :: i
      delg=gbox/(2.0d0*(nhis+1)) !The nhis+1 is in order to allow
!     the distances that fulfill (dist.lt.gbox/2.0d0) to have ig = int(dist/delg)=nhis
!     if it were nhis, the last value of g(r) would be zero instead of the actual value
      do i=1,nhis
         g(i)=0
      end do
      end subroutine grinit

      subroutine griter(delg,g,r,npart,box)
      implicit none
      integer, intent(in) :: npart,box
      real(8), dimension(:) :: g
      real(8), intent(in) :: r(3,1000)
      real(8), intent(in) :: delg
      real(8) :: rr2,dist,rijl,rij(3)
      integer :: ig,i,j,l
      do i=1,npart-1
         do j=i+1,npart
            rr2 = 0.d0
            do l = 1,3
               rijl = r(l,j) - r(l,i)
               rij(l) = rijl - box*dnint(rijl/box)
               rr2 = rr2 + rij(l)*rij(l)
            end do
            dist = sqrt(rr2)
            if (dist.lt.box/2) then
               ig = int(dist/delg)
               g(ig) = g(ig) + 2
            end if
         end do
      end do
      end subroutine griter

      subroutine griterforce(dist,gbox,g,delg,nhis)
      integer, intent(in) :: nhis
      real(8), intent(inout) :: g(nhis)
      real(8), intent(in) :: delg,gbox,dist
      if (dist.lt.gbox/2.0d0) then
         ig = int(dist/delg) 
         g(ig) = g(ig) + 2
      end if
      end subroutine griterforce
      
      subroutine grend(unit,nhis,delg,g,npart,ro,ngr)
      implicit none
      integer, intent(in) :: npart,nhis,ngr,unit
      real(8), intent(inout),dimension(:) :: g
      real(8), intent(in) :: delg,ro
      real(8) :: r,vb,nid,pi
      integer :: i
      pi=4*datan(1.d0)
      do i=1,nhis
         r=delg*(i-0.5)
         vb = (i**3-(i-1)**3)*delg**3
         nid = (4.0/3)*pi*vb*ro
         g(i) =g(i)/(ngr*npart*nid)
         write(unit,*) r,g(i)
      end do
      end subroutine grend

*********************************
!     Subroutines to integrate
*********************************
      subroutine simpson(vec,dt,I)
      implicit none
      real(8), intent(in), dimension(:):: vec
      real(8), intent(in) :: dt
      real(8), intent(out) :: I
      integer :: j,n
      n=size(vec)
      I=vec(1)+vec(n)
      do j=2,n-1
         if (mod(j,2)==0) THEN
            I=I+2*vec(j)
         else
            I=I+4*vec(j)
         end if
      end do
      I=(dt/3.0)*I
      end subroutine simpson

      subroutine trapezoidal(vec,dt,I)
      implicit none
      real(8), intent(in), dimension(:):: vec
      real(8), intent(in) :: dt
      real(8), intent(out) :: I
      integer :: j,n
      n=size(vec)
      I=vec(1)+vec(n)
      do j=2,n-1
         I=I+2*vec(j)
      end do
      I=(dt/2.0)*I
      end subroutine trapezoidal

************************************
!     Calculation of vacf and msd
************************************
      subroutine diffusion(ntel,npart,nconf,it0,t0,t0max,time0,
     &     x0,vx0,x,vx,vacf,r2t,ntime)
      integer, intent(in) :: npart,nconf,it0,t0max,ntel
      integer, intent(inout) :: t0,time0(t0max)
      integer, intent(inout),dimension(nconf) :: ntime
      real(8), intent(in), dimension(1000,3) :: x,vx
      real(8), intent(inout), dimension(npart,3,t0max) :: x0,vx0
      real(8), intent(inout), dimension(nconf) :: vacf,r2t
      integer :: tt0,i,t,delt,j,intime
      
      if (mod(ntel,it0).eq.0 .OR. ntel==1) then
         t0 = t0 + 1
         tt0 = mod(t0-1,t0max) + 1
         time0(tt0) = ntel
         do i = 1,npart
            x0(i,:,tt0) = x(i,:)
            vx0(i,:,tt0) = vx(i,:)
         end do
      endif
      
      do t=1,min(t0,t0max)
         delt = ntel-time0(t)+1
         if (delt.lt.nconf) then
            intime=ntime(delt)
            ntime(delt) = intime + 1
            do i = 1,npart
               do j=1,3
                  vacf(delt) = vacf(delt) + vx(i,j)*vx0(i,j,t)
                  r2t(delt) = r2t(delt) + (x(i,j)-x0(i,j,t))**2
               end do
            end do
         endif
      end do
      end subroutine diffusion
      
      
      end module modfrog
