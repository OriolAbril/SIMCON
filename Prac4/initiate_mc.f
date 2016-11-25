      PROGRAM initiatemc
      use utilsmonte
      
!     1. Defining variables and dimensions
      implicit none
      integer, parameter :: unit=1
      integer :: ncycle,natoms
      real(8) :: temp,sigma,epsil,delta,ro


!     2. Reading data and computing related quantities

      open(1,file='lj-monte.data',status='old')
      read(1,*) ncycle
      read(1,*) natoms
      read(1,*) temp
      read(1,*) sigma,epsil
      read(1,*) delta
      read(1,*) ro
      close(1)

      
      open(unit,file='conf-monte.data',status='replace')
      call initmonte(unit,natoms,ro,sigma)
      close(unit)

      end program initiatemc
