      PROGRAM suma
      real :: a=22.34
      integer :: b=5
      open(1,file="suma.dat",status="replace")
!     Escribim amb el format predeterminat
      write(1,*) (a+b)
!     Escribim amb format enginyeering
      write(1,"(EN19.8E3)") (a+b)
!     Escribim amb format cientific
      write(1,"(ES19.8E4)") (a+b)
!     Escribim amb el format decimal
      write(1,"(F9.6)") (a+b)
      close(1)
      
      END PROGRAM suma
