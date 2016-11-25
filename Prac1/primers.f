        PROGRAM PRIMERS
        
        DIMENSION PRIMER(1000)        
        REAL*8 PIV
        OPEN (5, FILE = 'primers.dat')

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        N = 2
        L = 2
        PRIMER(1) = 1
        PRIMER(2) = 2
10      N = N + 1
        DO I = 2, N-1
           IF (MOD(DFLOAT(N),DFLOAT(I)).EQ.0) THEN
              GOTO 10  
           ENDIF
        ENDDO
        L = L + 1 
        IF (L.EQ.1000) THEN 
           GOTO 20 
        ENDIF
        PRIMER(L) = N
        GOTO 10
20      CONTINUE
        WRITE(5,*) (PRIMER(K),K=1,L-1)
        STOP
        END  
