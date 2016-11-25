!    Primer programa amb fortran que serveix per calcular pi 
      function leibniz(cpi,k)
      Implicit none
      Integer k
      double precision cpi,leibniz
      cpi=cpi+4.0*((-1.0)**k)/(2.0*dble(k)+1.0)
      return
      end function leibniz

      Program leib
      integer :: i,j,k=0
      double precision :: pi,cpi,leibniz
      pi=datan(1.0d0)*4.0
      print *, pi
      do while (abs(cpi-pi).gt.1e-9) 
      cpi=leibniz(cpi,k)
      k=k+1
      end do
      
      print *, k
      print *, cpi
      end program leib
