      module utilsmonte
      
      CONTAINS
*********************************************************
*********************************************************
c              subroutine reduced
*********************************************************
*********************************************************

      subroutine reduced(natoms,r,boxlength,delta,
     &     sigma)
      implicit none 
      integer, intent(in) :: natoms
      real(8), dimension(3,natoms), intent(inout) :: r
      real(8), intent(in) :: sigma
      real(8), intent(inout) :: boxlength,delta
      integer :: is,l

      boxlength = boxlength/sigma
      delta = delta/sigma

      do is = 1,natoms
         do l = 1,3
            r(l,is) = r(l,is)/sigma
         end do
      end do

      end subroutine reduced

      subroutine initmonte(unit,natoms,ro,sigma)
!     initialization of a cubic cell. ro has to be in the same units of sigma
      implicit none
      integer, intent(in) :: unit, natoms
      real(8), intent(in) :: ro,sigma
      real(8) :: vol, box,space,pos(3)
      integer :: i
      vol = natoms/ro
      box = vol**(1.0d0/3.0d0)*sigma
      space = box/natoms**(1.0d0/3.0d0)
      do i=1,3
         pos(i)=0.0d0
      end do
      do i=1,natoms
         write(unit,*) pos
         pos(1)=pos(1)+space
         if (pos(1).gt.box) then
            pos(1)=0
            pos(2)=pos(2)+space
            if (pos(2).gt.box) then
               pos(2)=0
               pos(3)=pos(3)+space
            end if
         end if
      end do
      write(unit,*) box
      end subroutine initmonte
      
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
      
      subroutine mcsample(natoms,r,boxlength,rc,epot,press,g,
     &     delg,nhis)
      implicit none
      integer, intent(in) :: nhis,natoms
      real(8), dimension(nhis) :: g
      real(8), dimension(3,natoms) :: r
      real(8), intent(in) :: boxlength,rc,delg
      real(8), intent(out) :: press,epot
      integer :: is,js
      real(8) :: pot,p,ep,pr
      ep = 0.d0
      pr = 0.d0

c      atom-atom interactions

      do is = 1,natoms-1
         do js = is+1,natoms
            call lj(natoms,is,js,r,boxlength,rc,pot,p,g,delg,nhis)
            ep = ep + pot
            pr = pr + p
         end do
      end do
      epot=ep
      press=pr

      end subroutine mcsample

      subroutine mcmove(natoms,r,box,rc,delta,beta,energy,etail,move)
      implicit none
      integer, intent(in) :: natoms
      real(8), dimension(3,natoms), intent(inout) :: r
      integer, intent(inout) :: move
      real(8), intent(inout) :: energy
      real(8), intent(in) :: box,rc,delta,beta,etail
      real(8), dimension(3,natoms) :: rtry
      real(8) :: ran,U,accept
      integer :: atom,l
      call random_number(ran)
      rtry=r
      atom=nint(ran*natoms+0.5)
      do l=1,3
         call random_number(ran)
         rtry(l,atom)=r(l,atom)+delta*(ran-0.5d0)
         if (rtry(l,atom).lt.0) rtry(l,atom) = rtry(l,atom) + box
         if (rtry(l,atom).gt.box) rtry(l,atom) = rtry(l,atom) - box
      end do

      call potenergy(rtry,U,etail,natoms,rc,box)
      accept=min(1.0d0,dexp(-beta*(U-energy)))
      call random_number(ran)
      if (ran.lt.accept) then
         do l=1,3
            r(l,atom)=rtry(l,atom)
         end do
         energy=U
         move=move+1
      end if 
      end subroutine mcmove

      subroutine potenergy(r,epot,etail,natoms,rc,box)
      implicit none
      integer, intent(in) :: natoms
      real(8), dimension(3,natoms), intent(inout) :: r
      real(8), intent(out) :: epot
      real(8), intent(in) :: rc,box,etail
      real(8) :: pot
      integer :: is,js
      epot = 0.d0
      
c     atom-atom interactions

      do is = 1,natoms-1
         do js = is+1,natoms
            call ljenergy(natoms,is,js,r,box,rc,pot)
            epot = epot + pot
         end do
      end do
      epot=epot+etail
      end subroutine potenergy

      subroutine potatom(atom,rold,rnew,eold,enew,etail,natoms,rc,box)
      implicit none
      integer, intent(in) :: natoms,atom
      real(8), dimension(3,natoms), intent(in) :: rold,rnew
      real(8), intent(out) :: enew
      real(8), intent(in) :: rc,box,etail,eold
      real(8) :: pot,ediffo,ediffn
      integer :: is
      ediffo = 0.d0
      ediffn = 0.d0
      
c     atom-atom interactions

      do is = 1,atom-1
         call ljenergy(natoms,is,ATOM,rold,box,rc,pot)
         ediffo = ediffo + pot
      end do
      do is = atom+1,natoms
         call ljenergy(natoms,is,atom,rold,box,rc,pot)
         ediffo = ediffo + pot
      end do

      do is = 1,atom-1
         call ljenergy(natoms,is,ATOM,rnew,box,rc,pot)
         ediffn = ediffn + pot
      end do
      do is = atom+1,natoms
         call ljenergy(natoms,is,atom,rnew,box,rc,pot)
         ediffn = ediffn + pot
      end do
      
      
      enew=eold+etail-ediffo+ediffn
      end subroutine potatom
      
      subroutine ljenergy(natoms,is,js,r,boxlength,rc,pot)
      implicit none
      integer, intent(in) :: is,js,natoms
      real(8), dimension(3,natoms) :: r
      real(8), dimension(3) :: rij
      real(8),intent(in) :: boxlength,rc
      real(8), intent(out) :: pot
      real(8) :: rr2,rijl,rr,ynvrr2,ynvrr6,ynvrr12
      integer :: l
      rr2 = 0.d0
      pot = 0.d0
      do l = 1,3
         rijl = r(l,js) - r(l,is)
         rij(l) = rijl - boxlength*dnint(rijl/boxlength)
         rr2 = rr2 + rij(l)*rij(l)
      end do

      rr = dsqrt(rr2)
      if (rr.lt.rc) then
         ynvrr2 = 1.d0/rr2
         ynvrr6 = ynvrr2*ynvrr2*ynvrr2
         ynvrr12 = ynvrr6*ynvrr6
         pot = 4.d0*(ynvrr12-ynvrr6)
      end if
      end subroutine ljenergy
      
*********************************************************
*********************************************************
c              subroutine Lennard-Jones
*********************************************************
*********************************************************

      subroutine lj(natoms,is,js,r,boxlength,rc,pot,p,g,delg,nhis)
      implicit none
      integer, intent(in) :: nhis,is,js,natoms
      real(8), dimension(3,natoms) :: r
      real(8), dimension(nhis), intent(inout) :: g
      real(8), dimension(3) :: rij
      real(8),intent(in) :: boxlength,rc,delg
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
      call griterforce(rr,g,delg,nhis)
      if (rr.lt.rc) then
         ynvrr2 = 1.d0/rr2
         ynvrr6 = ynvrr2*ynvrr2*ynvrr2
         ynvrr12 = ynvrr6*ynvrr6
         forcedist = 24.d0*(2.d0*ynvrr12-ynvrr6)*ynvrr2
         pot = 4.d0*(ynvrr12-ynvrr6)
         p=rr2*forcedist
      end if


      end subroutine lj


*********************************
!     Radial distribution function
*********************************
      
      subroutine grinit(nhis,g,gbox,delg)
      integer, intent(in) :: nhis
      real(8), dimension(nhis) :: g
      real(8), intent(in) :: gbox
      real(8),intent(out) :: delg
      integer :: i
      delg=gbox/(2.0d0*(nhis+1))                   !The nhis+1 is in order to allow
!     the distances that fulfill (dist.lt.gbox/2.0d0) to have ig = int(dist/delg)=nhis
!     if it were nhis, the last value of g(r) would be zero instead of the actual value
      do i=1,nhis
         g(i)=0
      end do
      end subroutine grinit

      subroutine griterforce(dist,g,delg,nhis)
      integer, intent(in) :: nhis
      real(8), intent(inout) :: g(nhis)
      real(8), intent(in) :: delg,dist
      ig = int(dist/delg) 
      if (ig.lt.(nhis+1)) then
         g(ig) = g(ig) + 2
      end if
      end subroutine griterforce
      
      subroutine grend(unit,nhis,delg,g,npart,ro,ngr)
      implicit none
      integer, intent(in) :: npart,nhis,ngr,unit
      real(8), intent(inout),dimension(nhis) :: g
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
      
      end module utilsmonte
