C                                                                 
C                                                                   
C                                                                 
               PROGRAM WATER                                     
C                                                                  
C         + + + + + + + + + + + + + + + + + + + + + + + + + + +   
C         +                                                   +  
C         +                                                   + 
C         +      PROGRAMA PRINCIPAL DE DINAMICA MOLECULAR     +
C         +               LIQUIDS POLIATOMICS                 +    
C         +           PREPARAT PER A 216 MOLECULES            +   
C         +                                                   +  
C         +       Criteri geometric per                       +
C         +         definir els ponts d'hidrogen              +
C         +                                                   +    
C         +       CORRECCIO DE VELOCITATS PELS GRAUS          +    
C         +          DE LLIBERTAT INTERNS I DEL CM            +     
C         +                 PER SEPARAT                       +    
C         +                                                   +   
C         +     ALGORISME D'INTEGRACIO : LEAP-FROG VERLET     +  
C         +         AMB ACOBLAMENT A UN "BANY TERMIC"         + 
C         +     Berendsen et al. J.Chem.Phys.81,3684(1984)    +  
C         +                                                   +    
C         +        TRACTANT LES FORCES DE LLARG ABAST         +     
C         +               AMB LA SUMA D'EWALD                 +    
C         +                                                   +   
C         + + + + + + + + + + + + + + + + + + + + + + + + + + +  
C                                                               
C                                                              
C                                                             
        IMPLICIT REAL*8 (A-H,O-Z)                            
        REAL*8 KB,MASSA                              
        COMMON/CONF/R(3,3,250),NPASS                      
        COMMON/LEAPF/VA(3,3,216)                         
        COMMON/DADES/TREF,PREF,MASSA(5),NATS,N          
        COMMON/BANY/KB,TAU,NF,NFCM,NFINT               
        COMMON/NTP/NTP                                
        DIMENSION PTOTAL(3)                          
c	
        open(5,file='input',status='old')
	open(60,file='result.dat')
	open(11,file='temp.dat')   
	open(17,file='data.in')   
	open(27,file='data.dat')
	open(1,file='gdroo.dat')
        open(2,file='gdroh.dat')
	open(3,file='gdrhh.dat')
        open(10,file='favo.dat')
        open(20,file='favh.dat')
        open(31,file='r2o.dat')
        open(32,file='r2h.dat')
c
        READ(5,*) NPL,NP11,NTP,KINI,KFCOR    
        READ(5,*) TREF,PREF,N,NATS,(MASSA(J),J=1,NATS)          
C                                                              
        WRITE(60,111)                                          
        WRITE(60,222) NPL,NP11,NTP,KINI,KFCOR  
        WRITE(60,333) TREF,PREF,N,NATS,(MASSA(J),J=1,NATS)        
C                                                               
        INPL = 0                                               
        INP11= 0                                              
        NPASS = 0                                        
        EPSLOW = 0.D0                                  
C
        IF (KINI.EQ.0) THEN                         
           CALL INICIA(COSTAT)                     
        ELSE                                      
           DO I=1,N                              
              DO J=1,NATS                                      
                 READ(17,*) (VA(L,J,I),L=1,3)                    
                 READ(17,*) (R(L,J,I),L=1,3)                      
              END DO                                                
           END DO                                                  
           READ(17,*) COSTAT                                      
           WRITE(60,1200)                                         
           READ(5,*)                                            
           READ(5,*)                                           
        END IF                                                
c
        dens = (N*180.D0)/(6.023D0*COSTAT*COSTAT*COSTAT) 
c
        write(60,*) ' DENSITAT =  ',  dens
C                                                            
        COSINV = 1.0D0/COSTAT                               
C                                                          
        CALL PREVIA(COSTAT,RF,RF2,ROCUT2,DT,UEMOL) 
C                                                       
        IF (KFCOR.EQ.0) THEN                           
           CALL IFNCOR(N,NATS,DT,RF,NTP)              
        END IF                                       
C                                                   
C                                                  
100     CONTINUE                                  
C                                                
C             *     *    PART ITERATIVA DE CALCUL    *     *       
C                                                                 
C                                                                 
        CALL FASFOR(EPFAST)                                         
C                                                                  
        CALL SLOFOR(COSTAT,COSINV,RF2,ROCUT2,EPSLOW,UEMOL)   
C                                                             
        CALL NCONFI(DT,COSTAT,ECIN,ECIN1,ECIN2)            
        NPASS = NPASS + 1                                          
C                                                                   
        INP11 = INP11 + 1                                          
        IF (INP11.EQ.NP11) THEN                                   
           TEMP = 2.0D0*ECIN/(DFLOAT(NF)*KB)                     
           TEMP1 = 2.0D0*ECIN1/(DFLOAT(NFCM)*KB)                
           TEMP2 = 2.0D0*ECIN2/(DFLOAT(NFINT)*KB)              
           ECIN = (ECIN/DFLOAT(N))/UEMOL                      
           EPOT = ((EPSLOW+EPFAST)/DFLOAT(N))/UEMOL          
           ETOTAL = ECIN + EPOT                             
           WRITE(11,*) TEMP1,TEMP2,TEMP                    
           WRITE(11,*) TEMP,ECIN,EPOT,ETOTAL              
           INP11 = 0                                     
           INPL = INPL + 1                              
           IF (INPL.EQ.NPL) THEN                       
              CALL CNTMOM(PTOTAL)                     
              WRITE(60,2000) NPASS,TEMP,EPOT,ETOTAL   
              WRITE(60,2010) PTOTAL                  
              INPL = 0                             
           END IF                                 
        END IF                                   
C                                               
        IF (NPASS.EQ.NTP) GO TO 150                              
        GO TO 100                                               
C                                                                 
150     CONTINUE                                              
C                                                           
        IF (KFCOR.EQ.0) THEN                               
           CALL GFNCOR(N,NATS,COSTAT,NPASS,KIG)           
        END IF                                           
C                                                       
        DO I=1,N                                       
           DO J=1,NATS                                
              WRITE(27,*) (VA(L,J,I),L=1,3)          
              WRITE(27,*) (R(L,J,I),L=1,3)          
           END DO                                  
        END DO                                    
        WRITE(27,*) COSTAT                       
C                                               
        close(5)
	close(60)
	close(11)
	close(17)
	close(27)
	close(1)
	close(2)
	close(3)
        close(10)
        close(20)
        close(31)
        close(32)
C                                        
111     FORMAT(///16X,'D I N A M I C A    M O L E C U L A R'///7X, 
     &  'A L G O R I S M E     L E A P - F R O G    V E R L E T'//2X,
     &  'A M B   A C O B L A M E N T',                              
     &  '   A   U N   " B A N Y   T E R M I C "'                   
     &  //////18X,'*  *  *    D A D E S    *  *  *'/)             
222     FORMAT(//3X,'NPL',5X,'NP11',5X,'NTP',   
     &  4X,'KINI',3X,'KFCOR'/                        
     &  ' ',I5,4X,I4,4X,I5,3X,I6,3X,I7)           
333     FORMAT(/3X,'TEMP REF :',D12.4,' K',7X,'PRESS REF :',        
     &  D11.3,' atm'//6X,'No. MOLECULES :',I4,5X,                  
     &  'No. ATOMS PER MOLECULA :',I3//3X,                        
     &  'MASSES ATOMIQUES :',5F9.4/)                             
1200    FORMAT(///' ','JA ESTA LLEGIDA LA CONFIGURACIO INICIALITZADORA')
2000    FORMAT(' ','No. PASSOS :',I7/' ','TEMP :',D12.5,' K',   
     &  3X,'EPOT :',D13.6,3X,'ETOTAL :',D13.6,' Kcal/mol')     
2010    FORMAT(8X,'QUANT MOVIMENT ',3(3X,D13.6))              
        STOP                                                       
        END                                                       
C                                                                
C******************************************************************
C                                                                   
        SUBROUTINE PREVIA(COSTAT,RF,RF2,ROCUT2,DT,UEMOL)        
        IMPLICIT REAL*8(A-H,O-Z)                                 
        REAL*8 KB,KC,NAV,MASSA,KANGL,KRANGL,KRR              
        COMMON/DADES/TREF,PREF,MASSA(5),NATS,N                     
        COMMON/SPC/AO,CO,CARR(5),ROO2                             
        COMMON/EWALD/AEW,DOSQPI,DOSPIL,H2MAX,V2C                 
        COMMON/MASINV/YMASSA(5)                                 
        COMMON/TAURAH/DOH,ALPHA,KANGL,KRANGL,KRR,DISOH,DISHH,  
     #                YMASS1,YMASS2,YMASS3                    
        COMMON/BANY/KB,TAU,NF,NFCM,NFINT                    
        COMMON/BOND/UBONDMAX,ROHCUTOFF,PSIMAX
C                                                        
        READ(5,*) RF,AEWL,H2MAX                                    
        READ(5,*) AO,CO,(CARR(J),J=1,NATS)                       
        READ(5,*) DOH,ALPHA,KANGL,KRANGL,KRR,DISOH,DISHH        
        READ(5,*) DT,TAU,ROO                             
        READ(5,*) ROCUTOFF,UBONDMAX,ROHCUTOFF,PSIMAX         
C                                                          
        UMA = 1.6605655D-27                               
        UT = 1.0D-12                                     
        UL = 1.0D-10                                    
        CAL = 4.184D0                                  
        NAV = 6.02204D23                              
        KB = 1.38062D-23                             
        E = 1.6021D-19                              
        KC = 8.988D9                               
C                                                 
        UEA = UMA*((UL*UL)/(UT*UT))              
        UEMOL = 1.0D3*CAL/UEA/NAV                                   
        UFOR = UMA*(UL/(UT*UT))*1.0D3*1.0D2/1.0D-3                 
C                                                                 
        WRITE(60,1000) AO,CO,(CARR(J),J=1,NATS)                   
        WRITE(60,1010) DOH,ALPHA,KANGL,KRANGL,KRR,DISOH,DISHH    
C                                                              
        RF = RF*COSTAT                                        
        WRITE(60,1100) RF,AEWL,H2MAX                          
        WRITE(60,1030) ROO                                   
        WRITE(60,1101) ROCUTOFF, UBONDMAX         
C                                                       
        AO = AO*UEMOL                                  
        CO = CO*UEMOL                                 
        DO J=1,NATS                                                
           CARR(J) = CARR(J)*E*DSQRT(KC)/DSQRT(UL)/DSQRT(UEA)      
        END DO                                                    
        DOH = DOH/UFOR                                              
        KANGL = KANGL/UFOR                                         
        KRANGL = KRANGL/UFOR                                      
        KRR = KRR/UFOR                                           
        DT = DT/UT                                              
        TAU = TAU/UT                                           
        WRITE(60,1200) DT,TAU                            
        KB = KB/UEA                                          
        ROO2 = ROO*ROO                                      
        RF2 = RF*RF                                        
        ROCUT2 = ROCUTOFF*ROCUTOFF                        
C                                                        
        CALL HVECTS(N,NATS,AEWL,COSTAT)                 
C                                                      
        DO J=1,NATS                                   
           YMASSA(J) = 1.0D0/MASSA(J)                
        END DO                                      
        YMASS1 = YMASSA(1)                         
        YMASS2 = YMASSA(2)                        
        YMASS3 = YMASSA(3)                       
        NF = 3*NATS*N - 3                       
        NFCM = 3*N - 3                         
        NFINT = NF - NFCM                     
C                                            
1000    FORMAT(//16X,'PARAMETRES DEL POTENCIAL   :::'//4X,       
     &  'AO :',D11.4,' Kcal*A12/mol',5X,'CO :',D11.4,' Kcal*A6/mol'
     &  /10X,'CARREGUES :',3D13.3,' e'/)                           
1010    FORMAT (/15X,'DOH :',D10.3,' mdyn A',5X,'alpha :',D12.4,' A-1'
     &  //5X,'K angle :',D12.4,' mdyn A-1 ',                       
     &  5X,'K r-angle :',D12.4,' mdyn A-1',                         
     &  5X,'K r-r :',D12.4,' mdyn A-1'//                           
     &  5X,'dist.  O-H  H-H :',2D16.6, 'A'/)                      
1030    FORMAT(//5X,                                             
     &  'RADI DE LA 1a. CAPA D''HIDRATACIO :',D11.4,' A'//)     
1100    FORMAT(/' ','RF :',F9.4,' A',5X,'ALFA*L:',F6.2,5X,     
     &              'H2 MAX :',F6.2)                                
1101    FORMAT(/' ','ROCUTOFF:',F8.3,' A',4X,    
     &              'UBONDMAX :',F8.3,' KJ/mol'//)    
1032    FORMAT(/' ','DIST. AL QUADR. O-H BONDED:',F8.4,' A2',3X,   
     &  'DIST. AL QUADR. O-H NO BONDED:',F8.4,' A2'//               
     &  'DIST. AL QUADR. O-H BRIDGE:',F8.4,' A2',3X,                
     &  'DIST. AL QUADR. O-O BONDED:',F8.4,' A2')                  
1200    FORMAT(//' ','PAS D''INTEGRACIO :', D10.3,' ps',          
     &  5X,'PARAMETRE D''ACOBLAMENT AL BANY :',D10.3,' ps'/)    
        RETURN                                                 
        END                                                   
C                                                            
C*****************************************************************
C                                                                  
        SUBROUTINE HVECTS(N,NATS,AEWL,COSTAT)                      
        IMPLICIT REAL*8(A-H,O-Z)                                  
        INTEGER HMU                                           
        COMMON/EWALD/AEW,DOSQPI,DOSPIL,H2MAX,V2C                   
        COMMON/PAIRDAT/V2Cw
        COMMON/HH/HMU(3,300),ADH2(300),NDH                        
        COMMON/SPC/AO,CO,CARR(5),ROO2                            
C                                                               
        AEW = AEWL/COSTAT                                      
        PI = 4.0D0*ATAN(1.0D0)                                
        DOSQPI = 2.0D0/DSQRT(PI)                             
        DOSPIL = 2.0D0*PI/COSTAT                            
        PI2 = PI*PI                                        
        COST2 = COSTAT*COSTAT                             
        AEWL2 = AEW*AEW*COST2                            
C                                                       
        NARREL=INT(DSQRT(H2MAX))                       
        NOMBRE=NARREL*2+1                             
        IH=0                                         
        DO I=1,NOMBRE                              
           DO J=1,NOMBRE                          
              DO K=1,NARREL+1                    
                 I1=NARREL+1-I                             
                 J1=NARREL+1-J                                 
                 K1=K-1                                  
                 H2=I1**2+J1**2+K1**2                               
                 IF(H2.EQ.0) GO TO 10                              
                 IF(K1.EQ.0) THEN                                 
                    IF(I1.LT.0) GO TO 10                         
                    IF(I1.EQ.0.AND.J1.LT.0) GO TO 10            
                 END IF                                        
                 IF (H2.LE.H2MAX) THEN                        
                    IH = IH + 1                              
                    HMU(1,IH) = I1                                
                    HMU(2,IH) = J1                               
                    HMU(3,IH) = K1                               
                    ADH2(IH) = 4.0D0/COST2*DEXP(-PI2*H2/AEWL2)/H2  
                 END IF                                  
10               CONTINUE                                           
              END DO                                               
           END DO                                                 
        END DO                                                   
C                                                               
        NDH = IH                                               
C                                                             
        WRITE(60,2000) NDH                                    
C                                                           
        V2C = 0.D0                                         
        DO I=1,N                                          
           DO J=1,NATS                                   
              V2C = V2C - CARR(J)*CARR(J)               
           END DO                                      
        END DO                                        
        V2C = V2C*AEW/DSQRT(PI)                      
c
c    Valor per a 1 unica aigua
c
        V2Cw = 0.0
        DO J=1,NATS
           V2Cw = V2Cw - CARR(J)*CARR(J)
        END DO
c
        V2Cw = V2Cw*AEW/DSQRT(PI)
c
2000    FORMAT(/10X,'No. VECTORS H :',I5/)                         
        RETURN                                                    
        END                                                     
C                                                                   
C******************************************************************
C                                                               
        SUBROUTINE INICIA(COSTAT)                                 
        IMPLICIT REAL*8(A-H,O-Z)                              
        REAL*8 KB,MASSA,MTOTAL                                      
        COMMON/CONF/R(3,3,250),NPASS                               
        COMMON/LEAPF/VA(3,3,216)                                  
        COMMON/DADES/TREF,PREF,MASSA(5),NATS,N                   
        DIMENSION RAUX(3,700),PTOTAL(3)                         
C                                                              
        UMA = 1.6605655D-27                                   
        UL = 1.0D-10                                         
        UT = 1.0D-12                                        
        UEA = UMA*((UL*UL)/(UT*UT))                        
        KB = 1.380662D-23/UEA                             
        URO = UMA/(UL*UL*UL)/1.0D3                       
C                                                       
        DO I=1,N                                       
           DO J=1,NATS                                
              DO L=1,3                               
                 VA(L,J,I) = 0.D0                   
                 R(L,J,I) = 0.D0                   
              END DO                              
           END DO                                
        END DO                                  
        READ(5,*) RO,NCOST,KDIST,XIN,IX,NBOND  
        WRITE(60,111)                          
C                                            
C          &&&&         CONFIGURACIO INICIAL         &&&&        
C                                                               
C                                                             
C           ***     GENERACIO DE LES VELOCITATS   ***               
C                                                                  
        YR   = DFLOAT(IX)                                         
        VMV  = 0.D0                                              
        MTOTAL = 0.D0                                           
        DO J=1,NATS                                            
           MTOTAL = MTOTAL + MASSA(J)                         
        END DO                                               
        NF = (3*NATS-NBOND)*N - 3                           
        ECIN = 0.D0                                        
        DO I=1,N                                          
           DO L=1,3                                      
              SV=DSQRT(DFLOAT(NF)*KB*TREF/(3.0D0*MTOTAL*DFLOAT(N)))
              CALL GAUSS(YR,SV,VMV,VCM)                         
              DO J=1,NATS                                           
                 VA(L,J,I) = VCM*DSQRT(MTOTAL/DFLOAT(NATS)/MASSA(J)) 
                 ECIN = ECIN + .5D0*(MASSA(J)*VA(L,J,I)*VA(L,J,I))  
              END DO                                               
           END DO                                                 
        END DO                                                   
        TEMP  = 2.0D0*ECIN/(DFLOAT(NF)*KB)                      
        FCV2   = DSQRT(TREF/TEMP)                              
        DO L=1,3                                              
           PTOTAL(L) = 0.D0                                  
        END DO                                              
        DO I=1,N                                           
           DO J=1,NATS                                    
              DO L=1,3                                   
                 VA(L,J,I) = VA(L,J,I)*FCV2             
                 PTOTAL(L) = PTOTAL(L) + MASSA(J)*VA(L,J,I)     
              END DO                                         
           END DO                                                   
        END DO                                                     
        DO L=1,3                                                  
           PTOTAL(L) = PTOTAL(L)/DFLOAT(N)                       
        END DO                                                  
        DO I=1,N                                               
           DO J=1,NATS                                        
              DO L=1,3                                       
                 VA(L,J,I) = VA(L,J,I) - PTOTAL(L)/MTOTAL         
              END DO                                             
           END DO                                               
        END DO                                                 
        ECIN = 0.D0                                           
        DO L=1,3                                             
           PTOTAL(L) = 0.D0                                 
        END DO                                             
        DO I=1,N                                          
           DO J=1,NATS                                   
              DO L=1,3                                  
                 ECIN = ECIN + .5D0*MASSA(J)*VA(L,J,I)*VA(L,J,I) 
                 PTOTAL(L) = PTOTAL(L) + MASSA(J)*VA(L,J,I)         
              END DO                                               
           END DO                                                 
        END DO                                                   
        DO L=1,3                                                
           PTOTAL(L) = PTOTAL(L)/DFLOAT(N)                     
        END DO                                                
        TEMP  = 2.0D0*ECIN/(DFLOAT(NF)*KB)                   
        WRITE(60,1000) FCV2,TEMP,PTOTAL                         
C                                                                   
C                                                                  
C      ***     GENERACIO SIMETRICA DE LES POSICIONS     ***       
C                                                                
        COSTAT = (DFLOAT(N)*MTOTAL/(RO/URO))**(1.0D0/3.0D0)     
        JJ = 1                                                 
        XCOST = COSTAT/DFLOAT(NCOST)                          
        DO I=1,NCOST                                         
           DO J=1,NCOST                                     
              DO K=1,NCOST                                 
                 RAUX(1,JJ) = DFLOAT(I-1)*XCOST           
                 RAUX(2,JJ) = DFLOAT(J-1)*XCOST          
                 RAUX(3,JJ) = DFLOAT(K-1)*XCOST         
                 JJ = JJ + 1                           
              END DO                                  
           END DO                                    
        END DO                                      
C                                                                 
        NAUX = NCOST*NCOST*NCOST                               
        IF (KDIST.EQ.0) THEN                                        
           WRITE (6,112) RO,IX,NCOST,KDIST,XIN                     
C                                                                 
C            *   MODIFICACIO DE L'ESTRUCTURA SIMETRICA   *       
C                                                               
           XIN = XIN*XCOST                                     
           YR = DFLOAT(IX)                                    
           DO I=1,NAUX                                       
              DO J=1,3                                      
                 YR = DMOD(1.96290D4*YR+1.0D0,2.684354560D8)
                 YFL = SNGL(YR*3.725290D-9)                
                 X = -1.D0                                
                 IF (YFL.GT.0.5D0) X=.1D01               
                 YR = DMOD(1.96290D4*YR+1.0D0,2.684354560D8)    
                 YFL = SNGL(YR*3.725290D-9)                         
                 RAUX(J,I) = RAUX(J,I) + X*YFL*XIN                 
              END DO                                              
           END DO                                                
             ELSE                                               
           WRITE(60,113) RO,IX,NCOST,KDIST                      
        END IF                                                
        WRITE(60,2000) COSTAT                                 
C                                                                  
        IF (NAUX.NE.N) THEN                                       
           DO I=1,N                                             
 1202         YR = DMOD(1.96290D4*YR+1.0D0,2.684354560D8)           
              YFL = SNGL(YR*3.725290D-9)                           
              J = NAUX*YFL + 1                                    
              IF (RAUX(1,J).EQ.0.D0) GO TO 1202                  
              DO L=1,3                                          
                 R(L,1,I) = RAUX(L,J)                          
              END DO                                          
              RAUX(1,J) = 0.D0                               
           END DO                                           
               ELSE                                        
           DO I=1,N                                       
              DO L=1,3                                   
                 R(L,1,I) = RAUX(L,I)                   
              END DO                                   
           END DO                                     
        END IF                                       
C                                                   
        CALL ATOMS(N,NATS,COSTAT)                  
C                                                
111     FORMAT(/////5X,'*  *  *    ',           
     &  'C O N F I G U R A C I O    I N I C I A L     *  * *')      
112     FORMAT(//10X,'DENSITAT',5X,'IX',5X,'NCOST',5X,'KDIST',5X,'XIN'/
     &  8X,D11.4,3X,I3,6X,I3,8X,I2,2X,D9.2/)                       
113     FORMAT(//10X,'DENSITAT',5X,'IX',5X,'NCOST',5X,'KDIST'/     
     &  8X,D11.4,3X,I3,6X,I3,8X,I2/)                                
1000    FORMAT(//' ','F.C.VEL :',D10.3,3X,'TEMP :',D12.5,3X,       
     &  'QUANT. MOV. :',3D12.4)                                   
2000    FORMAT(/15X,'& &       COSTAT :',D14.6,' A      &  &'/)  
4000    FORMAT(////7X,'VELOCITAT',14X,'POSICIO')                
4100    FORMAT(/5X,'M O L E C U L A',I4)                       
4110    FORMAT(/5X,'ATOM ',I4/)                               
4200    FORMAT(2(4X,D16.9))                                  
        RETURN                                              
        END                                                
C                                                         
C*******************************************************************
C                                                                   
C       Aquesta rutina genera variables aleatories amb una distribucio 
C                                                                  
C       p(X) = exp(-((X-mu)*(X-mu)/(2*sigma*sigma))/(2*pi*sigma*sigma) 
C                                                                  
C                                                                 
        SUBROUTINE GAUSS(SEED,SIGMA,MU,X)                           
        IMPLICIT REAL*8(A-H,O-Z)                                   
        REAL*8 MU                                                 
        A = 0.D0                                                 
        DO 1010 I=1,12                                          
           SEED = DMOD(1.96290D4*SEED + 1.0D0,2.684354560D8)   
           YFL = SNGL(SEED*3.725290D-9)                       
           A = A + YFL                                       
 1010   CONTINUE                                            
        X = (A-.6D01)*SIGMA + MU                           
        RETURN                                            
        END                                              
C                                                       
C*******************************************************************
C                                                                  
        SUBROUTINE ATOMS(N,NATS1,COSTAT)                           
        IMPLICIT REAL*8(A-H,O-Z)                                 
        COMMON/CONF/R(3,3,250),NPASS                                
        DIMENSION AN(3)                                            
C                                                                 
        READ(5,*) DOH,AHOH                                       
        WRITE(60,1000) DOH,AHOH                                  
        PI = 4.D0*DATAN(1.0D0)                                 
        ANGLE = (AHOH*PI/180.0D0)/2.0D0                       
        XO = - 2.0D0*DOH*DCOS(ANGLE)/18.0D0                  
        YO = 0.0D0                                          
        ZO = 0.0D0                                         
        XH1 = XO + DOH*DCOS(ANGLE)                        
        YH1 = DOH*DSIN(ANGLE)                            
        ZH1 = 0.0D0                                     
        XH2 = XH1                                      
        YH2 = - YH1                                   
        ZH2= 0.0D0                                                  
C                                                                   
        DO I=1,N                                                   
           XCM= R(1,1,I)                                          
           YCM= R(2,1,I)                                            
           ZCM= R(3,1,I)                                           
           YR = 0.D0                                              
           YR = DMOD(1.96290D4*YR + 1.0D0 , 2.684354560D8)       
           YFL   = SNGL(YR*3.725290D-9)                         
           AN(1) = PI*YFL                                      
           YR = DMOD(1.96290D4*YR + 1.0D0 , 2.684354560D8)    
           YFL   = SNGL(YR*3.725290D-9)                      
           AN(2) = 2.0D0*PI*YFL                             
           YR = DMOD(1.96290D4*YR + 1.0D0 , 2.684354560D8) 
           YFL   = SNGL(YR*3.725290D-9)                   
           AN(3) = 2.0D0*PI*YFL                          
           A11 = COS(AN(3))*COS(AN(2))-COS(AN(1))*SIN(AN(2))*SIN(AN(3)) 
           A12 =-SIN(AN(3))*COS(AN(2))-COS(AN(1))*SIN(AN(2))*COS(AN(3)) 
           A21 = COS(AN(3))*SIN(AN(2))+COS(AN(1))*COS(AN(2))*SIN(AN(3)) 
           A22 =-SIN(AN(3))*SIN(AN(2))+COS(AN(1))*COS(AN(2))*COS(AN(3)) 
           A31 = SIN(AN(1))*SIN(AN(3))                              
           A32 = SIN(AN(1))*COS(AN(3))                              
           R(1,1,I) = XCM + XO*A11                                  
           R(2,1,I) = YCM + XO*A21                                 
           R(3,1,I) = ZCM + XO*A31                                 
           R(1,2,I) = XCM + XH1*A11 + YH1*A12                      
           R(2,2,I) = YCM + XH1*A21 + YH1*A22                     
           R(3,2,I) = ZCM + XH1*A31 + YH1*A32                       
           R(1,3,I) = XCM + XH2*A11 + YH2*A12                      
           R(2,3,I) = YCM + XH2*A21 + YH2*A22                     
           R(3,3,I) = ZCM + XH2*A31 + YH2*A32                    
        END DO                                                  
C                                                              
        DO I=1,N                                             
           DO L=1,3                                         
              IF (R(L,1,I).LT.0.D0) THEN                       
                 DO J=1,NATS1                                  
                    R(L,J,I) = R(L,J,I) + COSTAT                   
                 END DO                                             
              END IF                                                
           END DO                                                  
        END DO                                                      
1000    FORMAT(/' ','DISTANCIA  O-H :',D13.6,' A',                  
     &          5X,' ANGLE  H-O-H  :',D13.6,' deg ')                
        RETURN                                                      
        END                                                         
C                                                                   
C******************************************************************
C                                                                   
        SUBROUTINE CNTMOM(PTOTAL)                                  
        IMPLICIT REAL*8(A-H,O-Z)                                   
        REAL*8 MASSA                                          
        COMMON/VELS/V(3,3,216)                                      
        COMMON/DADES/TREF,PREF,MASSA(5),NATS,N                     
        DIMENSION PTOTAL(3)                                         
C                                                                   
        DO L=1,3                                                    
           PTOTAL(L) = 0.D0                                         
        END DO                                                      
C                                                                   
        DO I=1,N                                                   
           DO J=1,NATS                                              
              DO L=1,3                                              
                 PTOTAL(L) = PTOTAL(L)+MASSA(J)*V(L,J,I)            
              END DO                                                
           END DO                                                   
        END DO                                                      
        DO L=1,3                                                  
           PTOTAL(L) = PTOTAL(L)/DFLOAT(N)                          
        END DO                                                      
C                                                                   
        RETURN                                                      
        END                                                        
C                                                                   
C****************************************************************   
C                                                                   
        SUBROUTINE SLOFOR(COSTAT,COSINV,RF2,ROCUT2,EPSLOW,UEMOL)    
        IMPLICIT REAL*8 (A-H,O-Z)                                  
        REAL*8 MASSA
        INTEGER HMU                                                 
        COMMON/CONF/R(3,3,250),NPASS                              
        COMMON/ACCELS/ASLOW(3,3,216),AFAST(3,3,216)                 
        COMMON/DADES/TREF,PREF,MASSA(5),NATS,N                      
        COMMON/MASINV/YMASSA(5)                                     
        COMMON/SPC/AO,CO,CARR(5),ROO2                              
        COMMON/EWALD/AEW,DOSQPI,DOSPIL,H2MAX,V2C                    
        COMMON/PAIRDAT/V2Cw
        COMMON/HH/HMU(3,300),ADH2(300),NDH                          
        COMMON/HIDRAT/NMOLH(216),IMOLH(75,216),ISOLV(216,216)       
        COMMON/SELF/RRELW(3,3,250),DISTW(3,250)                     
        COMMON/BOND/UBONDMAX,ROHCUTOFF,PSIMAX 
        COMMON/HBOS/NOXINT(4000,5)  
        DIMENSION USIN(3,216),UCOS(3,216)
        DIMENSION FTOT(3,3,216),RJJ(3,5,5),R2(5,5),XF(5,5),         
     &            TMS(3),RREL(3),RREL12(3),RREL13(3),RREL23(3)      
C                                                                   
C                                                                   
        EPSLOW = 0.D0                                               
        KOX = 0
        KOH = 0                                                    
c
        DO I=1,N
           DO J=1,NATS                                              
              DO L=1,3                                             
                 FTOT(L,J,I) = 0.D0                       
              END DO                                       
           END DO                                               
        END DO                                                   
        DO IA=1,N                                               
           NMOLH(IA) = 0                                           
           DO I=1,25                                               
              IMOLH(I,IA) = 0                                     
           END DO
        END DO                                                      
        DO IA=1,N
           DO IB=1,N                                               
              ISOLV(IA,IB) = 0                                  
           END DO                                              
        END DO                                               
C                                                           
        DO IA=1,N-1                                            
           DO IB=IA+1,N                                       
              DIST2 = 0.D0                                      
C                                                      
              DO L=1,3                                            
                 RREL(L) = R(L,1,IA)-R(L,1,IB)                     
                 TMS(L) = COSTAT*ANINT(RREL(L)*COSINV)              
                 RREL(L) = RREL(L) - TMS(L)                        
                 DIST2 = DIST2 + RREL(L)*RREL(L)                    
              END DO                                               
              IF (DIST2.LT.RF2) THEN                               
                 DO JB=1,NATS                                      
                    DO JA=1,NATS                                    
                       R2(JA,JB) = 0.D0                             
                       DO L=1,3                                     
                          RJJ(L,JA,JB)=R(L,JA,IA)-R(L,JB,IB)       
     &                                           -TMS(L)          
                          R2(JA,JB) = R2(JA,JB) +                 
     &                                RJJ(L,JA,JB)*RJJ(L,JA,JB)    
                       END DO                                      
                    END DO                                         
                 END DO                                             
                 CALL FORCES(IA,IB,NATS,R2,XF,YEP)
                 DO JB=1,NATS
                    DO JA=1,NATS
                       DO L=1,3
                          FJJ = XF(JA,JB)*RJJ(L,JA,JB)
                          FTOT(L,JA,IA) = FTOT(L,JA,IA) + FJJ
                          FTOT(L,JB,IB) = FTOT(L,JB,IB) - FJJ
                       END DO
                    END DO
                 END DO
                 EPSLOW = EPSLOW + YEP
              END IF
C
C     PRIMER DRAFT D'OXIGENS POSSIBLES CANDIDATS A O---H- - -O
C        (ESTAN DINS D'UNA ESFERA DE RADI ROCUTOFF)
C
                 IF (DIST2.LT.ROCUT2) THEN
                    KOX = KOX + 1
                    NOXINT(KOX,1) = IA
                    NOXINT(KOX,2) = IB
                 END IF
c
           END DO                                                  
        END DO                                                      
C
C       MIREM SI HI HA PONTS D'HIDROGEN : criteri geometric !!!
C
        IF (KOX.GE.1) THEN
           CALL HBONDS(KOX,N,NATS,COSTAT,COSINV)
        END IF
C                                                                       
C       SELF-ENERGY TERM                                                
C                                                                       
        YEP = 0.D0                                                      
        DO I=1,N                                                        
           R12 = DISTW(1,I)                                             
           R13 = DISTW(2,I)                                             
           R23 = DISTW(3,I)                                             
           DO L=1,3                                                     
              RREL12(L) = RRELW(L,1,I)                                  
              RREL13(L) = RRELW(L,2,I)                                  
              RREL23(L) = RRELW(L,3,I)                                  
           END DO                                                       
           AEWR1 = AEW*R12                                              
           AEWR2 = AEW*R13                                              
           AEWR3 = AEW*R23                                             
	   CALL TRANEF(AEWR1,U1)
           CALL TRANEF(AEWR2,U2)
       	   CALL TRANEF(AEWR3,U3)
           XF12 =  CARR(1)*CARR(2)*                                     
     $            (- U1/(R12*R12) +                    
     $             DOSQPI*AEW/R12*DEXP(- AEWR1*AEWR1) )                 
           XF13 =  CARR(1)*CARR(3)*                                     
     $            (- U2/(R13*R13) +                    
     $             DOSQPI*AEW/R13*DEXP(- AEWR2*AEWR2) )                 
           XF23 =  CARR(2)*CARR(3)*                                     
     $            (- U3/(R23*R23) +                    
     $             DOSQPI*AEW/R23*DEXP(- AEWR3*AEWR3) )                 
           DO L=1,3                                                     
              FTOT(L,1,I) = FTOT(L,1,I)                                 
     $                      + XF12*RREL12(L) + XF13*RREL13(L)          
              FTOT(L,2,I) = FTOT(L,2,I)                               
     $                      - XF12*RREL12(L) + XF23*RREL23(L)        
              FTOT(L,3,I) = FTOT(L,3,I)                                 
     $                      - XF13*RREL13(L) - XF23*RREL23(L)           
           END DO                                                       
           YEP = YEP - CARR(1)*CARR(2)*U1/R12         
     $               - CARR(1)*CARR(3)*U2/R13          
     $               - CARR(2)*CARR(3)*U3/R23          
        END DO                                                          
        EPSLOW = EPSLOW + YEP                                          
C                                                                       
        YEP = 0.D0
c
        DO IIH=1,NDH
           ADH2IIH=ADH2(IIH)
           HMUX=HMU(1,IIH)
           HMUY=HMU(2,IIH)
           HMUZ=HMU(3,IIH)
           SUMSIN = 0.0D0
           SUMCOS = 0.0D0
           DO I=1,N
              DO J=1,NATS
                 THETA = HMUX*R(1,J,I)
                 THETA = THETA + HMUY*R(2,J,I)
                 THETA = THETA + HMUZ*R(3,J,I)
                 THETA=DOSPIL*THETA
                 USIN(J,I)=CARR(J)*SIN(THETA)
                 UCOS(J,I)=CARR(J)*COS(THETA)
                 SUMSIN = SUMSIN + USIN(J,I)
                 SUMCOS = SUMCOS + UCOS(J,I)
              END DO
           END DO
           DO I=1,N
              DO J=1,NATS
                 AA=(USIN(J,I)*SUMCOS-UCOS(J,I)*SUMSIN)*ADH2IIH
                 FTOT(1,J,I) = FTOT(1,J,I) + AA*HMUX
                 FTOT(2,J,I) = FTOT(2,J,I) + AA*HMUY
                 FTOT(3,J,I) = FTOT(3,J,I) + AA*HMUZ
              END DO
           END DO
           YEP = YEP + ADH2IIH*(SUMCOS*SUMCOS + SUMSIN*SUMSIN)
        END DO
c
        DO I=1,N
           DO J=1,NATS
              DO L=1,3
                 ASLOW(L,J,I) = FTOT(L,J,I)*YMASSA(J)
              END DO
           END DO
        END DO
c
        PI = 4.0D0*ATAN(1.0D0)                                          
        YEP = (.25D0*COSTAT/PI)*YEP                                     
        EPSLOW = EPSLOW + YEP + V2C                                     
        RETURN                                                         
        END                                                           
C
C*************************************************************************
C
	SUBROUTINE TRANEF(X,U)
	IMPLICIT REAL*8 (A-H,O-Z)
	INTEGER    J,NCFC,NCFD
	DIMENSION C(18), D(17)
C
	DATA NCFC, NCFD / 18, 17 /
	DATA SQRTPI / 1.7724538509055160 /
	DATA C
     x  / 1.9449071068178803,4.20186582324414D-2,-1.86866103976769D-2,
     x  5.1281061839107D-3, -1.0683107461726D-3, 1.744737872522D-4,
     x  -2.15642065714D-5, 1.7282657974D-6, -2.00479241D-8,
     x  -1.64782105D-8, 2.0008475D-9, 2.57716D-11, -3.06343D-11,
     x  1.9158D-12, 3.703D-13, -5.43D-14, -4.0D-15, 1.2D-15 /
C
        DATA D
     x  / 1.4831105640848036,-3.010710733865950D-1,6.89948306898316D-2,
     x  -1.39162712647222D-2, 2.4207995224335D-3, -3.658639685849D-4,
     x  4.86209844323D-5, -5.7492565580D-6, 6.113243578D-7,
     x  -5.89910153D-8, 5.2070091D-9, -4.232976D-10, 3.18811D-11,
     x  -2.2361D-12, 1.467D-13, -9.0D-15, 5.0D-16 /
C  
        XV = ABS(X)
        IF(XV .GE. 6.25) THEN
        U = 1.
        IF(X.LT.0) U = -1.
	RETURN
	END IF
C
        IF(XV .GT. 2.0) THEN
           X2 = 2. - 20./(3.+XV)
           BJP2 = 0.
           BJP1 = C(NCFC)
	   J = NCFC - 1
           DO WHILE (J.GT.1)
           BJ = X2*BJP1 - BJP2 + C(J)
           BJP2 = BJP1
           BJP1 = BJ
	   J = J - 1
           END DO
C
           BJ = X2*BJP1 - BJP2 + C(J)
	   X2 = 0.5 * (BJ-BJP2) / XV * EXP (-X*X) / SQRTPI
	   U = 1. - X2
           IF(X .LT. 0.0) U = -U
	   RETURN
C
        ELSE
C
	   X2 = X*X - 2.
           BJP2 = 0.
           BJP1 = D(NCFD)
           J = NCFD - 1
           DO WHILE (J.GT.1)
           BJ = X2*BJP1 - BJP2 + D(J)
           BJP2 = BJP1
	   BJP1 = BJ
           J = J - 1
           END DO
C
           BJ = X2*BJP1 - BJP2 + D(J)
           U = 0.5 * (BJ-BJP2) * X
        END IF
C
        RETURN
        END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C                                                                    
        SUBROUTINE HBONDS(KOX,N,NATS,COSTAT,COSINV)
        IMPLICIT REAL*8 (A-H,O-Z)
        COMMON/CONF/R(3,3,250),NPASS
        COMMON/BOND/UBONDMAX,ROHCUTOFF,PSIMAX
        COMMON/HBOS/NOXINT(4000,5)
        COMMON/HBONS/NOHBOND(2000,5),NOHNBOND(2000,5),KOH,KNOH
        COMMON/HBONS2/NHBCOUNT(5000,7,2),IND(5000)
        COMMON/BANG/PSI(20000)
        COMMON/PERCENT/rat(200000)
        COMMON/HBFAV/KHB0,NHNBS(5000,2),KHB1,NHBS(5000,2)
        COMMON/LINBIF/nhblin(10000,2),nhbbif(10000,2)
        DIMENSION IO1(2), IO2(2)
        DIMENSION ROHPB(3), ROHMOL(3), ROOPB(3)
C
        PI = 4.0D0*DATAN(1.0D0)
        KOH = 0
C
        DO JBO = 1, KOX
           IO1(1) = NOXINT(JBO,1)
           IO2(1) = NOXINT(JBO,2)
           IO1(2) = NOXINT(JBO,2)
           IO2(2) = NOXINT(JBO,1)
C
           DO K = 1,2
              IOA = IO1(K)
              IOB = IO2(K)
              DO J = 2, 3
                 ROHPB2 = 0.D0
                 ROOPB2 = 0.D0
                 ROHMOL2 = 0.D0
                 ROOPEROH = 0.D0
                 DO L = 1,3
                    ROHPB(L) = R(L,1,IOA) - R(L,J,IOB)
                    ROHMOL(L) = R(L,J,IOB) - R(L,1,IOB)
                    ROOPB(L) = R(L,1,IOA) - R(L,1,IOB)
C
                    ROHPB(L)=ROHPB(L)-COSTAT*ANINT(ROHPB(L)/COSTAT)
                    ROHMOL(L)=ROHMOL(L)-COSTAT*ANINT(ROHMOL(L)/COSTAT)
                    ROOPB(L)=ROOPB(L)-COSTAT*ANINT(ROOPB(L)/COSTAT)
C
                    ROHPB2 = ROHPB2 + ROHPB(L)*ROHPB(L)
                    ROOPB2 = ROOPB2 + ROOPB(L)*ROOPB(L)
                    ROHMOL2 = ROHMOL2 + ROHMOL(L)*ROHMOL(L)
                 END DO
                 DISROHPB = DSQRT(ROHPB2)
                 DISROO = DSQRT(ROOPB2)
                 DISROHM = DSQRT(ROHMOL2)
                 PMOD = (ROOPB2+ROHMOL2-ROHPB2)/(2.D0*DISROO*DISROHM)
                 PSI2 = DACOS(PMOD)
                 PSIDEG = (PSI2*DFLOAT(180)/PI)
C
                 IF (DISROHPB.LE.ROHCUTOFF) THEN
                   IF (PSIDEG.LE.PSIMAX) THEN
C
C     KOH ES L'INDEX QUE ENS COMPTA ELS HB QUE HI HA.
C     NOHBOND(KOH,1) DIU QUE L'HB NUM. KOH CORRESPON A
C     LA MOLECULA IOB; NOHBOND(KOH,2)
C     ENS DIU AMB QUINA ALTRA MOLECULA (O) FORMA PONT.
C     NOHBOND(KOH,3) indica quin dels H's de NOHBOND(KOH,1)
C     forma el HB.
C
                        KOH = KOH + 1
                        PSI(KOH) = PSIDEG
                        NOHBOND(KOH,1) = IOB
                        NOHBOND(KOH,2) = IOA
                        NOHBOND(KOH,3) = J
                   END IF
                 END IF
C
              END DO
           END DO
        END DO
C
        DO II = 1,KOH
           DO IJ = 1,KOH
              IF (NOHBOND(IJ,1).GT.NOHBOND(II,1)) THEN
                 VARIABLE1 = NOHBOND(II,1)
                 NOHBOND(II,1) = NOHBOND(IJ,1)
                 NOHBOND(IJ,1) = VARIABLE1
                 VARIABLE2 = NOHBOND(II,2)
                 NOHBOND(II,2) = NOHBOND(IJ,2)
                 NOHBOND(IJ,2) = VARIABLE2
                 VARIABLE3 = NOHBOND(II,3)
                 NOHBOND(II,3) = NOHBOND(IJ,3)
                 NOHBOND(IJ,3) = VARIABLE3
              END IF
           END DO
        END DO
C
C       ANEM A VEURE QUINS HIDROGENS DELS POSSIBLES RESULTA QUE
C       FORMEN PONT I ELS HEM DE DESCARTAR PELS NO-ENLLASATS.
C       TAMBE PRENEM NOTA DELS QUE NO FORMEN PONT AMB SEGURETAT.
C
        KNOH = 0
C
        DO JOX = 1,N
           DO JH = 2,NATS
              JNH = 1
98            IF (JOX.EQ.NOHBOND(JNH,1)) THEN
                 IF (JH.EQ.NOHBOND(JNH,3)) THEN
                    GO TO 89
                 END IF
              ENDIF
79            IF (JOX.NE.NOHBOND(JNH,1)) THEN
                 IF (JNH.EQ.KOH) THEN
                    KNOH = KNOH + 1
                    NOHNBOND(KNOH,1) = JOX
                    NOHNBOND(KNOH,2) = JH
                    GO TO 89
                 END IF
                 JNH = JNH + 1
                 GO TO 98
              END IF
89            CONTINUE
           END DO
        END DO
C
C       SEPAREM ARA ELS OXIGENS PER CATEGORIES:
C       ELS QUE EN FORMEN 0,1,2,3,4,5,6
C       IND(JK) INDICA QUANTS PONTS FORMA L'OXIGEN JK.
C
        DO JK = 1, N
           IND(JK) = 0
           DO LK = 1, KOH
              IF (JK.EQ.NOHBOND(LK,1)) THEN
                 IND(JK) = IND(JK) + 1
                 NHBCOUNT(JK,IND(JK),1) = JK
                 NHBCOUNT(JK,IND(JK),2) = NOHBOND(LK,2)
              ENDIF
              IF (JK.EQ.NOHBOND(LK,2)) THEN
                 IND(JK) = IND(JK) + 1
                 NHBCOUNT(JK,IND(JK),1) = JK
                 NHBCOUNT(JK,IND(JK),2) = NOHBOND(LK,1)
              ENDIF
           ENDDO
        ENDDO
C
C       ABUNDANCIA DE PONTS D'HIDROGEN
c
        KOB0 = 0
        KOB1 = 0
        KOB2 = 0
        KOB3 = 0
        KOB4 = 0
        KOB5 = 0
        KOB6 = 0
c
        DO JK = 1,N
           IF (IND(JK).EQ.0) THEN
              KOB0 = KOB0 + 1
           END IF
           IF (IND(JK).EQ.1) THEN
              KOB1 = KOB1 + 1
           END IF
           IF (IND(JK).EQ.2) THEN
              KOB2 = KOB2 + 1
           END IF
           IF (IND(JK).EQ.3) THEN
              KOB3 = KOB3 + 1
           END IF
           IF (IND(JK).EQ.4) THEN
              KOB4 = KOB4 + 1
           END IF
           IF (IND(JK).EQ.5) THEN
              KOB5 = KOB5 + 1
           END IF
           IF (IND(JK).EQ.6) THEN
              KOB6 = KOB6 + 1
           END IF
        END DO
C
c        WRITE(55,'(7i5)') KOB0,KOB1,KOB2,KOB3,KOB4,KOB5,KOB6
        RETURN
        END 
C                                                                   
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c     calculation of the pair potential energy between molecules IA and IB
c
        subroutine pair(IA,IB,NATS,R2,UEMOL,vt)
        IMPLICIT REAL*8 (A-H,O-Z) 
        COMMON/PAIRDAT/V2Cw
        COMMON/CONF/R(3,3,250),NPASS
        COMMON/SPC/AO,CO,CARR(5),ROO2
        COMMON/EWALD/AEW,DOSQPI,DOSPIL,H2MAX,V2C
        COMMON/HH/HMU(3,300),ADH2(300),NDH
        COMMON/SELF/RRELW(3,3,250),DISTW(3,250)
        DIMENSION USIN(3),UCOS(3),R2(5,5)
c
c     One water 'IA' interacts with 'IB' 
c
        INTEGER HMU
c
        vt = 0.D0                    
c
c     We calculate the LJ + Coulomb pot. energy for the 
c      'IA'-'IB' interaction...for particles at distance < RF
c
        DO JA=1,NATS
           DO JB=1,NATS
              RAB2 = R2(JA,JB)
              RAB = DSQRT(RAB2)
              AEWRAB = AEW*RAB
              RAB3 = RAB*RAB2
              CALL TRANEF(AEWRAB,U)
              ERFC = 1.0D0 - U
              UAB = CARR(JA)*CARR(JB)*ERFC/RAB
              vt = vt + UAB
           END DO
        END DO
c
        coul = vt
c
        X = 1.0D0/(R2(1,1)*R2(1,1)*R2(1,1))
        Y = X/R2(1,1)
        alj  = X*(AO*X-CO)
        vt = vt + alj 
C                                                         
C       SELF-ENERGY TERM                                 
C                                                       
        YEP = 0.D0                                     
c        
c     1  Water...only
c
        R12 = DISTW(1,IA)
        R13 = DISTW(2,IA)
        R23 = DISTW(3,IA)
        AEWR1 = AEW*R12                                          
        AEWR2 = AEW*R13                                         
        AEWR3 = AEW*R23                                        
        call tranef(AEWR1,EAW1)
        call tranef(AEWR2,EAW2)
        call tranef(AEWR3,EAW3)
        YEP = YEP - CARR(1)*CARR(2)*EAW1/R12
     $            - CARR(1)*CARR(3)*EAW2/R13
     $            - CARR(2)*CARR(3)*EAW3/R23
c
        self = YEP 
c
        vt = vt + YEP                                 
c
c      Single water k-part...1 water only...
c
        YEP = 0.0
c 
        DO IIH=1,NDH
           ADH2IIH=ADH2(IIH)
           HMUX=HMU(1,IIH)
           HMUY=HMU(2,IIH)
           HMUZ=HMU(3,IIH)
           SUMSIN = 0.0D0
           SUMCOS = 0.0D0
c
           do J = 1, NATS
              THETA = HMUX*R(1,J,IA)
              THETA = THETA + HMUY*R(2,J,IA)
              THETA = THETA + HMUZ*R(3,J,IA)
              THETA=DOSPIL*THETA
              USIN(J)=CARR(J)*SIN(THETA)
              UCOS(J)=CARR(J)*COS(THETA)
              SUMSIN = SUMSIN + USIN(J)
              SUMCOS = SUMCOS + UCOS(J)
           enddo
c
           YEP = YEP + ADH2IIH*(SUMCOS*SUMCOS+SUMSIN*SUMSIN)
c
        enddo
c
        akp = YEP
c
        PI = 4.0D0*ATAN(1.0D0)                                   
        YEP = (.25D0*COSTAT/PI)*YEP                             
        vt = vt + YEP + V2Cw                            
c
        CAL = 4.184D0
c
        vt = vt*CAL/UEMOL
c
        RETURN
        END
C                                                  
C                    + + + + + + + + + + + + + + + + +                 
C                    + + + +        SPC0       + + + +                
C                    + + + + + + + + + + + + + + + + +               
C                                                                   
C                                                                  
        SUBROUTINE FORCES(IA,IB,NATS,R2,XF,YEPOT)           
        IMPLICIT REAL*8 (A-H,O-Z)                                
        COMMON/SPC/AO,CO,CARR(5),ROO2                           
        COMMON/EWALD/AEW,DOSQPI,DOSPIL,H2MAX,V2C               
        COMMON/DD0/DR(5),GDR(5,500),VOL(5,500),IZ(5),NG,NMDR(5,500),    
     #             NZ(5)                                               
        COMMON/HIDRAT/NMOLH(216),IMOLH(75,216),ISOLV(216,216)         
        DIMENSION R2(5,5),XF(5,5)                                    
C                                                                   
        YEPOT = 0.D0                                               
C                                                                 
        DO JA=1,NATS                                             
           DO JB=1,NATS                                                 
              RAB2 = R2(JA,JB)                                         
              RAB = DSQRT(RAB2)                                       
              AEWRAB = AEW*RAB                                       
              RAB3 = RAB*RAB2                                       
	      CALL TRANEF(AEWRAB,U)    
   	      ERFC = 1.0D0 - U
              UAB = CARR(JA)*CARR(JB)*ERFC/RAB                    
              XF(JA,JB) = CARR(JA)*CARR(JB)*( ERFC +                   
     #                  DOSQPI*AEWRAB*DEXP(-AEWRAB*AEWRAB) )/RAB3     
              YEPOT = YEPOT + UAB                                    
              IF (JA.EQ.1.AND.JB.EQ.1) THEN                         
C                DR(1) = 0.05                                      
                 II = (RAB/DR(1)) + 1                             
                 NMDR(1,II) = NMDR(1,II) + 2               
                    ELSE                                        
                 IF (JA.EQ.1.OR.JB.EQ.1) THEN                  
C                   DR(2) = 0.05                              
                    II = (RAB/DR(2)) + 1                     
                    NMDR(2,II) = NMDR(2,II) + 1       
                      ELSE                                 
C                   DR(3) = 0.05                          
                    II = (RAB/DR(3)) + 1                 
                    NMDR(3,II) = NMDR(3,II) + 2   
                 END IF                                
              END IF                                  
           END DO                                    
        END DO                                      
        X = 1.0D0/(R2(1,1)*R2(1,1)*R2(1,1))                             
        Y = X/R2(1,1)                                                  
        XF(1,1) = XF(1,1) + 6.0D0*(2.0D0*AO*X-CO)*Y                   
        YEPOT = YEPOT + X*(AO*X-CO)                                  
C                                                                   
        IF (R2(1,1).LT.ROO2) THEN                                  
           NMOLH(IA) = NMOLH(IA) + 1                              
           IMOLH(NMOLH(IA),IA) = IB                              
           NMOLH(IB) = NMOLH(IB) + 1                            
           IMOLH(NMOLH(IB),IB) = IA                            
           ISOLV(IA,IB) = 1                                   
           ISOLV(IB,IA) = 1                                  
        END IF                                                          
C                                                                      
        RETURN                                                        
        END                                                          
C                                                                   
C                + + + + + + + + + + + + + +                       
C                  ++      TAURAH0       ++                       
C                + + + + + + + + + + + + + +                     
C                                                               
        SUBROUTINE FASFOR(EPFAST)                              
        IMPLICIT REAL*8 (A-H,O-Z)                             
        REAL*8 MASSA,KANGL,KRANGL,KRR                        
        COMMON/CONF/R(3,3,250),NPASS                        
        COMMON/ACCELS/ASLOW(3,3,216),AFAST(3,3,216)        
        COMMON/DADES/TREF,PREF,MASSA(5),NATS,N            
        COMMON/TAURAH/DOH,ALPHA,KANGL,KRANGL,KRR,DISOH,DISHH,          
     #                YMASS1,YMASS2,YMASS3                            
        COMMON/SELF/RRELW(3,3,250),DISTW(3,250)                      
        COMMON/DD0/DR(5),GDR(5,500),VOL(5,500),IZ(5),NG,NMDR(5,500),
     #             NZ(5)                                           
        COMMON/DD45/DANGLE,IANMIN,IANMAX,PHOH(400),               
     #          VMHOH,SHOH2,VMOH,SOH2                                    
        DIMENSION RREL12(3),RREL13(3),RREL23(3)                 
c
C                                                              
        EPFAST = 0.D0                                         
        DO I=1,N                                          
           R1 = 0.D0                                     
           R2 = 0.D0                                    
           R3 = 0.D0                                       
           DO L=1,3                                    
              RREL12(L) = R(L,1,I) - R(L,2,I)         
              RREL13(L) = R(L,1,I) - R(L,3,I)        
              RREL23(L) = R(L,2,I) - R(L,3,I)       
              R1 = R1 + RREL12(L)*RREL12(L)        
              R2 = R2 + RREL13(L)*RREL13(L)    
              R3 = R3 + RREL23(L)*RREL23(L)   
           END DO                                
           R1 = DSQRT(R1)                    
           R2 = DSQRT(R2)                   
           R3 = DSQRT(R3)                                               
           DISTW(1,I) = R1                                             
           DISTW(2,I) = R2                                            
           DISTW(3,I) = R3                                           
           II = (R1/DR(2)) + 1                                      
           NMDR(2,II) = NMDR(2,II) + 1                             
           II = (R2/DR(2)) + 1                                    
           NMDR(2,II) = NMDR(2,II) + 1                           
           II = (R3/DR(3)) + 1                                  
           NMDR(3,II) = NMDR(3,II) + 2                         
c
           DO L=1,3                                           
              RREL12(L) = RREL12(L)/R1                       
              RREL13(L) = RREL13(L)/R2                      
              RREL23(L) = RREL23(L)/R3                     
              RRELW(L,1,I) = RREL12(L)                    
              RRELW(L,2,I) = RREL13(L)                   
              RRELW(L,3,I) = RREL23(L)                  
           END DO                                      
c
           COSHOH = 0.D0                              
           DO L=1,3                                  
              COSHOH = COSHOH + RREL12(L)*RREL13(L) 
           END DO                       
           ANGLE = DACOS(COSHOH)                  
           IAN = (ANGLE/DANGLE) + 1              
           PHOH(IAN) = PHOH(IAN) + 1            
           VMHOH = VMHOH + ANGLE               
           SHOH2 = SHOH2 + ANGLE*ANGLE        
           VMOH = VMOH + R1 + R2             
           SOH2 = SOH2 + R1*R1 + R2*R2      
           R1 = R1 - DISOH                 
           R2 = R2 - DISOH                
           R3 = R3 - DISHH      
c         
           EXP1 = DEXP(-ALPHA*R1)       
           EXP2 = DEXP(-ALPHA*R2)      
           XF12 = - 2.0D0*DOH*ALPHA*EXP1*(1.0D0-EXP1)                   
     #            - KRANGL*R3 - KRR*R2                                 
           XF13 = - 2.0D0*DOH*ALPHA*EXP2*(1.0D0-EXP2)                 
     #            - KRANGL*R3 - KRR*R1                               
           XF23 = - KANGL*R3 - KRANGL*(R1+R2)                       
c
           DO L=1,3                                                     
              AFAST(L,1,I) = (XF12*RREL12(L)+XF13*RREL13(L))*YMASS1    
              AFAST(L,2,I) =(-XF12*RREL12(L)+XF23*RREL23(L))*YMASS2   
              AFAST(L,3,I) =(-XF13*RREL13(L)-XF23*RREL23(L))*YMASS3  
           END DO                                                   
           EPFAST = EPFAST + DOH*( (1.0D0-EXP1)*(1.0D0-EXP1) +     
     &                            (1.0D0-EXP2)*(1.0D0-EXP2 ) )    
     &             + .5D0*KANGL*R3*R3 + KRANGL*R3*(R1+R2) + KRR*(R1*R2) 
        END DO                                                         
c
        RETURN                                                        
        END                                                          
C                                                                   
C                     + + + + + + + + + + + + + + + + +            
C                     + + + +       H2ONC       + + + +           
C                     + + + + + + + + + + + + + + + + +          
C                                                               
C                                                              
        SUBROUTINE NCONFI(DT,COSTAT,ECIN,ECIN1,ECIN2)         
        IMPLICIT REAL*8(A-H,O-Z)                             
        REAL*8 KB,MASSA,MTOTAL,LBDA1,LBDA2                  
        COMMON/CONF/R(3,3,250),NPASS                       
        COMMON/ACCELS/ASLOW(3,3,216),AFAST(3,3,216)       
        COMMON/LEAPF/VA(3,3,216)                         
        COMMON/SPC/AO,CO,CARR(5),ROO2                   
        COMMON/VELS/V(3,3,216)                         
        COMMON/DADES/TREF,PREF,MASSA(5),NATS,N        
        COMMON/BANY/KB,TAU,NF,NFCM,NFINT             
        COMMON/HIDRAT/NMOLH(216),IMOLH(75,216),ISOLV(216,216)         
        COMMON/HBONS/NOHBOND(2000,5),NOHNBOND(2000,5),KOH,KNOH
        COMMON/HBONS2/NHBCOUNT(5000,7,2),IND(5000)
        COMMON/HBFAV/KHB0,NHNBS(5000,2),KHB1,NHBS(5000,2) 
        COMMON/DD10/FAV(3,10000),V20(3),NPA10(100),NPTS10,K10,
     &       N10MAX,K10BIS,NA10,JA10,IA10,N10,I10,N10F,I10BIS
        COMMON/DD8/FDQM(5,10000),NPA8(60),NPTS8,K8,N8MAX,K8BIS,NA8,
     &             JA8,IA8,N8,I8,N8F,I8BIS  
        DIMENSION VS(3,3,216),R0(60,3,3,1728),TR(3,1728)   
        DIMENSION V0(60,3,3,216)
        DIMENSION VACM(3),VSCM(3),VCM(3)                          
        DIMENSION RJJ(3,5,5),R2(5,5),TMS(3),RREL(3)        
C                                                                
        MTOTAL = 0.D0                                           
        DO J=1,NATS                                            
           MTOTAL = MTOTAL + MASSA(J)                         
        END DO                                               
        YMTOT = 1.D0/MTOTAL                                 
C                                                          
        ECIN1 = 0.D0                                      
        ECIN2 = 0.D0                                     
        DO I=1,N                                       
           DO L = 1,3                                 
              VACM(L) = 0.D0                         
              DO J=1,NATS                                           
                VACM(L) = VACM(L) + MASSA(J)*VA(L,J,I)       
              END DO                                               
              VACM(L) = VACM(L)*YMTOT                            
              ECIN1 = ECIN1 + .5D0*MTOTAL*VACM(L)*VACM(L)    
           END DO                                           
           DO J=1,NATS                                     
              DO L = 1,3                                  
                 ECIN2 = ECIN2 + .5D0*MASSA(J)*(VA(L,J,I)-VACM(L))     
     #                                        *(VA(L,J,I)-VACM(L))    
              END DO                                                 
           END DO                                                   
        END DO                                                     
        TEMP1 = 2.0D0*ECIN1/(DFLOAT(NFCM)*KB)                     
        LBDA1 = DSQRT( 1.0D0 + (DT/TAU)*(TREF/TEMP1-1.0D0) )     
        TEMP2 = 2.0D0*ECIN2/(DFLOAT(NFINT)*KB)                  
        LBDA2 = DSQRT( 1.0D0 + (DT/TAU)*(TREF/TEMP2-1.0D0) )   
C                                                             
        DO I=1,N                                             
           DO L=1,3                                         
              VSCM(L) = 0.D0                               
           END DO                                         
           DO J=1,NATS                                   
              DO L=1,3                                  
                 VS(L,J,I) = VA(L,J,I) +               
     #                      (ASLOW(L,J,I)+AFAST(L,J,I))*DT             
              END DO                                                  
           END DO                                                    
           DO L=1,3                                                 
              DO J=1,NATS                                      
                 VSCM(L) = VSCM(L) + MASSA(J)*VS(L,J,I)            
              END DO                                              
              VSCM(L) = VSCM(L)*YMTOT                           
           END DO                                             
           DO J=1,NATS                                       
              DO L=1,3                                    
                 VS(L,J,I) = LBDA1*VSCM(L)               
     #                     + LBDA2*(VS(L,J,I)-VSCM(L))  
              END DO                                   
           END DO                                 
        END DO                                   
C                                                    
        ECIN = 0.D0                                 
        ECINCM = 0.D0                              
        DO I=1,N                                                        
           DO L=1,3                                                    
              VCM(L) = 0.D0                                           
           END DO                                                    
           DO J=1,NATS                                              
              DO L=1,3                                             
                 V(L,J,I) = .5D0*( VS(L,J,I)+VA(L,J,I) )         
                 ECIN = ECIN + .5D0*MASSA(J)*V(L,J,I)*V(L,J,I)  
                 R(L,J,I) = R(L,J,I) + VS(L,J,I)*DT            
                 VA(L,J,I) = VS(L,J,I)                        
              END DO                                         
           END DO                                         
           DO L=1,3                                         
              DO J=1,NATS                                  
                 VCM(L) = VCM(L) + MASSA(J)*V(L,J,I)                   
              END DO                                                  
              VCM(L) = VCM(L)*YMTOT                                  
              ECINCM = ECINCM + .5D0*MTOTAL*VCM(L)*VCM(L)           
           END DO                                                  
           DO L=1,3                                               
              IF (I8BIS.EQ.0) TR(L,I) = 0.D0   
              ZAB = R(L,1,I) - COSTAT                            
              IF ((R(L,1,I)*ZAB).GT.0.D0) THEN                  
                 IF (R(L,1,I).LT.0.D0) THEN                  
                    DO J=1,NATS                                
                       R(L,J,I) = R(L,J,I) + COSTAT           
                    END DO                                           
                    TR(L,I) = TR(L,I) + COSTAT   
                  ELSE                                              
                    DO J=1,NATS                                    
                       R(L,J,I) = R(L,J,I) - COSTAT               
                    END DO                                       
                    TR(L,I) = TR(L,I) - COSTAT
                  END IF                                        
              END IF                                           
           END DO                                             
        END DO                                               
C                                                             
        I10BIS = I10BIS + 1                                  
        IF (I10BIS.EQ.K10BIS) THEN                          
C                                                          
C              & & & &        F.A.V.      & & & &         
C                                                        
           I10 = I10 + 1                                
           IF (I10.EQ.K10) THEN                        
              N10 = N10 + 1                           
              DO I=1,N                               
                 DO J=1,NATS                        
                    DO L=1,3                                   
                       V0(N10,L,J,I) = V(L,J,I)               
                       V20(J) = V20(J) + V(L,J,I)*V(L,J,I)   
                    END DO                                  
                 END DO                                    
              END DO                                      
              IA10 = IA10 + 1                            
              I10 = 0                                   
           END IF                                      
           IF (JA10.LT.N10) THEN                      
              JB10 = JA10 + 1                        
              DO L10=JB10,N10                       
                 NPA10(L10) = NPA10(L10) + 1       
                 IJ = NPA10(L10)                  
                 IF (IJ.EQ.NPTS10) JA10 = L10    
                 DO I=1,N                       
                    DO J=1,NATS                
                       DO L=1,3                                        
                          FAV(J,IJ)=FAV(J,IJ)+V0(L10,L,J,I)*V(L,J,I)
                       END DO                                        
                    END DO                                          
                 END DO                                            
              END DO                                              
           END IF                                                
C                                                               
           IF (N10.EQ.N10MAX) THEN                                    
              JB10 = JA10 + 1
              DO L10=JB10,N10MAX                                       
                 NPA10(L10-JA10) = NPA10(L10)                         
                 DO I=1,N                                            
                    DO J=1,NATS                                     
                       DO L=1,3                                    
                          V0(L10-JA10,L,J,I) = V0(L10,L,J,I)      
                       END DO                                    
                    END DO                                      
                 END DO                                        
              END DO                                          
              N10 = N10MAX - JA10                            
              DO L10=N10+1,N10MAX                           
                 NPA10(L10) = 0                            
              END DO                                      
              N10F = N10F + JA10                         
              JA10 = 0                                
           END IF                                       
           I10BIS = 0                                  
        END IF                                                   
C
C       ELS INDEXS 8 CORRESPONEN A LA DESVIACIO QUADRATICA MITJANA
C
        I8BIS = I8BIS + 1
        IF (I8BIS.EQ.K8BIS) THEN
C
C              & & & &        F.D.Q.M.       & & & &
C
           I8 = I8 + 1
           IF (I8.EQ.K8) THEN
CCCCCCCC      IF (IA8.LT.NA8) THEN
                 N8 = N8 + 1
                 DO I=1,N
                    DO J=1,NATS
                       DO L=1,3
                          R0(N8,L,J,I) = R(L,J,I)
                       END DO
                    END DO
                 END DO
                 IA8 = IA8 + 1
                 I8 = 0
CCCCCCCC      END IF
           END IF
           IF (JA8.LT.N8) THEN
              JB8 = JA8 + 1
              DO L8=JB8,N8
                 NPA8(L8) = NPA8(L8) + 1
                 IJ = NPA8(L8)
                 IF (IJ.EQ.NPTS8) JA8 = L8
                 DO I=1,N
                    DO J=1,NATS
                       DO L=1,3
                        IF (IJ.NE.1) R0(L8,L,J,I)=R0(L8,L,J,I)+TR(L,I)
                        FDQM(J,IJ) = FDQM(J,IJ)+(R(L,J,I)-R0(L8,L,J,I))
     #                            *(R(L,J,I)-R0(L8,L,J,I))
                       END DO
                    END DO
                 END DO
              END DO
           END IF
           IF (N8.EQ.N8MAX) THEN
              JB8 = JA8 + 1
              DO L8=JB8,N8MAX
                 NPA8(L8-JA8) = NPA8(L8)
                 DO I=1,N
                    DO J=1,NATS
                       DO L=1,3 
                          R0(L8-JA8,L,J,I) = R0(L8,L,J,I)
                       END DO
                    END DO
                 END DO
              END DO
              N8 = N8MAX - JA8
              DO L8=N8+1,N8MAX
                 NPA8(L8) = 0
              END DO
              N8F = N8F + JA8
             JA8 = 0
           END IF
           I8BIS = 0
        END IF
C
        RETURN                                                 
        END                                                   
c
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C                                                                       
C
C                       * * * * * * * * * * * * *
C                          *     FUNCOR      *
C                       * * * * * * * * * * * * *
C
        SUBROUTINE IFNCOR(N,NATS,DT,RF,NPASS)
        IMPLICIT REAL*8(A-H,O-Z)
        COMMON/FUNCOR/KGRS,KFAV,KDQM,KGEOM
C
        READ(5,*) KGDRS,KFAV,KDQM,KGEOM
        IF (KGDRS.EQ.0) THEN
           CALL IGDRS(RF)
        END IF
        IF (KFAV.EQ.0) THEN
           CALL IFAV(NATS,NPASS)
        END IF
        IF (KDQM.EQ.0) THEN
           CALL IDQM(NATS,NPASS)
        END IF
        IF (KGEOM.EQ.0) THEN
           CALL IGEOM
        END IF
        RETURN
        END
c
C***********************************************************************
C
        SUBROUTINE GFNCOR(N,NATS,COSTAT,NPASS,KIG)
        IMPLICIT REAL*8(A-H,O-Z)
        COMMON/FUNCOR/KGDRS,KFAV,KDQM,KGEOM
c
        IF (KGDRS.EQ.0) THEN
           CALL GGDRS(N,COSTAT,NPASS)
        END IF
        IF (KFAV.EQ.0) THEN
           CALL GFAV(N,NATS)
        END IF
        IF (KDQM.EQ.0) THEN
           CALL GDQM(N,NATS,KIG)
        END IF
        IF (KGEOM.EQ.0) THEN
           CALL GGEOM(N,NPASS)
        END IF
C
        RETURN
        END
C
C***********************************************************************
C                                                               
C                         + + + + + + + + + +
C                          ++    FDQM    ++
C                         + + + + + + + + + +
C
C              Ha de ser   N8MAX*K8   mes gran que NPTS8
C
        SUBROUTINE IDQM(NATS,NTP)
        IMPLICIT REAL*8 (A-H,O-Z)
        COMMON/DD8/FDQM(5,10000),NPA8(60),NPTS8,K8,N8MAX,K8BIS,
     #              NA8,JA8,IA8,N8,I8,N8F,I8BIS
        READ(5,*) NPTS8,K8,N8MAX,K8BIS
        WRITE(60,5000) NPTS8,K8,N8MAX,K8BIS
        NA8 = 0
        DO IJ=1,NPTS8
           DO J=1,NATS
              FDQM(J,IJ) = 0.D0
           END DO
        END DO
        DO L8=1,N8MAX
           NPA8(L8) = 0
        END DO
C
        IA8 = 0
        JA8 = 0
        N8 = 0
        I8 = K8 - 1
        I8BIS = K8BIS - 1
        IF (NTP.GE.(NPTS8*K8BIS)) THEN
            NA8 = ((NTP-NPTS8*K8BIS)/(K8*K8BIS)) + 1
        END IF
        WRITE(60,5100) NA8
        N8F = 0
C
C C
5000    FORMAT(///10X,'******      ',
     #  'DESPLACAMENT QUADRATIC MITJA     ******'
     @  //' ','N PTS :',I6,6X,'K8 :',I5,6X,'N ORIGENS MAXIM :',I4,
     #  6X,'K8BIS :',I5)
5100    FORMAT(//21X,'PODREM CALCULAR :',I5,'  F.D.Q.M. COMPLETES  ')
        RETURN
        END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
        SUBROUTINE GDQM(N,NATS,KIG)
        IMPLICIT REAL*8(A-H,O-Z)
        COMMON/DD8/FDQM(5,10000),NPA8(60),NPTS8,K8,N8MAX,K8BIS,
     #              NA8,JA8,IA8,N8,I8,N8F,I8BIS
C
CCCCCCCCC
        KIG = 1
CCCCCCCCC
C       IF (NA8.NE.0) THEN
           IF (KIG.EQ.0) THEN
              DO IJ=1,NPTS8
                 DO J=1,NATS
                    FDQM(J,IJ) = FDQM(J,IJ)/DFLOAT(N*NA8)
                 END DO
              END DO
                ELSE
              MA8 = 0.D0
              DO IJ=1,NPTS8
                 MA8 = 0
                 DO L8=1,N8
                    IF (NPA8(L8).GE.IJ) MA8 = MA8 + 1
                 END DO
                 MA8 = MA8 + N8F
                 IF (MA8.NE.0) THEN
                    DO J=1,NATS
                       FDQM(J,IJ)=FDQM(J,IJ)/DFLOAT(N*MA8)
                    END DO
                 END IF
              END DO
           END IF
C
           WRITE(31,*) NPTS8
           WRITE(31,*) ( FDQM(1,IJ),IJ=1,NPTS8 )
           WRITE(31,*) N8F,N8,(NPA8(L8),L8=1,N8)
C
C
           WRITE(60,1000)
           IF (NPTS8.LT.100) THEN
              WRITE(60,1100)( FDQM(1,IJ),IJ=1,NPTS8 )
           ELSE
              WRITE(60,1100)( FDQM(1,IJ),IJ=1,100 )
           END IF
C
C             &&       pels hidrogens de l'aigua      &&
C
           DO IJ=1,NPTS8
              FDQM(2,IJ) = ( FDQM(2,IJ) + FDQM(3,IJ) )/2.0D0
           END DO  
C
C
           WRITE(32,*) NPTS8
           WRITE(32,*) ( FDQM(2,IJ),IJ=1,NPTS8 )
           WRITE(32,*) N8F,N8,(NPA8(L8),L8=1,N8)
C
C
           WRITE(60,1001)
           IF (NPTS8.LT.100) THEN
              WRITE(60,1100)( FDQM(2,IJ), IJ=1,NPTS8 )
           ELSE
              WRITE(60,1100)( FDQM(2,IJ), IJ=1,100 )
           END IF
C
1000    FORMAT(///25X, '**     ',
     &  'F . D . Q . M .    D E L S    O X I G E N S    **'//)
1001    FORMAT(///25X, '**     ',
     &  'F . D . Q . M .    D E L S    H I D R O G E N S    **'//)
1100    FORMAT(100(/' ',10F11.5))
        RETURN
        END
C
C***********************************************************************
C                                                                      
C                       + + + + + + + + + + + +               
C                          ++     GDRS    ++                 
C                       + + + + + + + + + + + +             
C                                                         
        SUBROUTINE IGDRS(RF)                             
        IMPLICIT REAL*8(A-H,O-Z)                                        
        COMMON/DD0/DR(5),GDR(5,500),VOL(5,500),IZ(5),NG,NMDR(5,500),   
     #             NZ(5)                                              
C                                                                    
        READ(5,*) NZM,NG,(DR(JJ),IZ(JJ),JJ=1,NG)                    
        PI = 4.0D0*DATAN(1.0D0)                                    
        WRITE(60,1000)                                             
        WRITE(60,1100) NG,(DR(JJ),JJ=1,NG)                        
C                                                               
        DO JJ=1,NG                                             
           NZ(JJ) = (RF/DR(JJ))                               
           DO I=IZ(JJ),NZ(JJ)                                
              GDR(JJ,I) = 0.D0                              
              NMDR(JJ,I) = 0                               
              RJ = DFLOAT(I-1)*DR(JJ)                     
              VOL(JJ,I) = (4.0D0/3.0D0)*PI*              
     #                   ( (RJ+DR(JJ))**3.0D0 - RJ**3.0D0)              
           END DO                                                      
        END DO                                                        
1000    FORMAT(///' ',10X,'******        ',                          
     @  'FUNCIONS DE DISTRIBUCIO RADIAL        *******'/)           
1100    FORMAT(/15X,'No. g.d.r. :',I4,5X,'AMPLADES :',5D10.2 /)    
1300    FORMAT(15X,                                               
     #         '&   LES TAULES DE VOLUMS JA ESTAN DISSENYADES   &'//)   
C                                                                      
        RETURN                                                        
        END                                                          
C                                                                   
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                                                                      
        SUBROUTINE GGDRS(N,COSTAT,NPASS)                              
        IMPLICIT REAL*8(A-H,O-Z)                                     
        COMMON/DD0/DR(5),GDR(5,500),VOL(5,500),IZ(5),NG,NMDR(5,500),
     #             NZ(5)                                           
        DIMENSION COEF(5),RMIN(5),NA(5),IV(5)                     
C                                                                
        VOLUM = COSTAT*COSTAT*COSTAT                            
C                                                              
C          &&        nomes per a l'aigua       &&             
C                                                            
        COEF(1) = VOLUM/DFLOAT(N*N)                         
        COEF(2) = VOLUM/DFLOAT(2*N*N)                      
        COEF(3) = VOLUM/DFLOAT(2*N*2*N)                   
C                                                        
        DO JJ=1,NG                                      
           IV(JJ) = 0                                  
           DO I=IZ(JJ),NZ(JJ)                         
              IF (NMDR(JJ,I).NE.0) THEN              
                 IF (IV(JJ).EQ.0) THEN              
                    NA(JJ) = I                                          
                    RMIN(JJ) = DFLOAT(I-1)*DR(JJ) + DR(JJ)/2.0D0       
                    IV(JJ) = IV(JJ) + 1                               
                 END IF                                              
                 GDR(JJ,I) = COEF(JJ)*                              
     #                       (NMDR(JJ,I)/DFLOAT(NPASS))/VOL(JJ,I)  
              END IF                                              
           END DO                                                
        END DO                                                  
C                                                              
        WRITE(60,1101)                                                   
        R1 = DFLOAT(IZ(1)-1)*DR(1) + DR(1)/2.0D0                       
        WRITE(60,1000) DR(1),NZ(1)                                     
        WRITE(60,2000) IZ(1),R1,NA(1),RMIN(1)                         
        WRITE(60,3000) ( GDR(1,I),I=IZ(1),NZ(1) )                    
        WRITE(1,*) DR(1),NA(1),NZ(1),RMIN(1),NPASS                 
        WRITE(1,*) (GDR(1,I),I=NA(1),NZ(1))                       
        WRITE(60,1102)                                            
        R1 = DFLOAT(IZ(2)-1)*DR(2) + DR(2)/2.0D0                
        WRITE(60,1000) DR(2),NZ(2)                              
        WRITE(60,2000) IZ(2),R1,NA(2),RMIN(2)                  
        WRITE(60,3000) ( GDR(2,I),I=IZ(2),NZ(2) )             
        WRITE(2,*) DR(2),NA(2),NZ(2),RMIN(2),NPASS          
        WRITE(2,*) (GDR(2,I),I=NA(2),NZ(2))                
        WRITE(60,1103)                                     
        R1 = DFLOAT(IZ(3)-1)*DR(3) + DR(3)/2.0D0         
        WRITE(60,1000) DR(3),NZ(3)                       
        WRITE(60,2000) IZ(3),R1,NA(3),RMIN(3)           
        WRITE(60,3000) ( GDR(3,I),I=IZ(3),NZ(3) )      
        WRITE(3,*) DR(3),NA(3),NZ(3),RMIN(3),NPASS   
        WRITE(3,*) (GDR(3,I),I=NA(3),NZ(3))         
C                                                  
C                                                 
1101  FORMAT(//25X,                                                     
     #'**       G. D . R.    O X I G E N - O X I G E N     ** ')       
1102  FORMAT(//25X, '**    ',                                         
     #'G.D.R.  O X I G E N - H I D R O G E N    **')                 
1103  FORMAT(//25X, '**    ',                                       
     #'G.D.R.    H I D R O G EN - H I D R O G E N      **')        
1000  FORMAT(//35X,'AMPLADA DIVISIONS :',D9.2,5X,'N. ZONES :',I4) 
2000  FORMAT(16X,'N. 1.A ZONA :',I4,5X,'R  1.A ZONA  :',F8.4,10X,
     *'NA : ',I4,5X,'R MIN =',F8.4/)                            
2100  FORMAT(32X,'N.  1.A ZONA  :',I4,10X,'R 1.A ZONA :',F8.4/)
3000  FORMAT(100(/' ',10F11.5))                               
      RETURN                                                 
      END                                                   
c
C                                                
C                         + + + + + + + + + +   
C                           ++    FAV    ++    
C                         + + + + + + + + + + 
C                                            
C              Ha de ser   N10MAX*K10   mes gran que NPTS10       
C                                                                
        SUBROUTINE IFAV(NATS,NTP)                               
        IMPLICIT REAL*8 (A-H,O-Z)                                    
        COMMON/DD10/FAV(3,10000),V20(3),NPA10(100),NPTS10,K10,
     #       N10MAX,K10BIS,NA10,JA10,IA10,N10,I10,N10F,I10BIS      
c
        READ(5,*) NPTS10,K10,N10MAX,K10BIS                        
        WRITE(60,5000) NPTS10,K10,N10MAX,K10BIS                   
c
        NA10 = 0                                                
        DO IJ=1,NPTS10                                         
           DO J=1,NATS                                       
              FAV(J,IJ) = 0                              
           END DO                                          
        END DO                                            
        DO L10=1,N10MAX                                  
           NPA10(L10) = 0                               
        END DO                                    
C                                                  
        IA10 = 0                                    
        JA10 = 0                                     
        N10 = 0                                       
        I10 = K10 - 1                                  
        I10BIS = K10BIS - 1                              
        IF (NTP.GE.(NPTS10*K10BIS)) THEN                
           NA10 = ((NTP-NPTS10*K10BIS)/(K10*K10BIS)) + 1  
        END IF                                              
        WRITE(60,5100) NA10                                   
        N10F = 0                                              
C                                                              
5000    FORMAT(///10X,'******      ',                            
     #  'FUNCIO D''AUTOCORRELACIO DE VELOCITATS O i H    ******'      
     @  //' ','N PTS :',I6,6X,'K10 :',I5,6X,'N ORIGENS MAXIM :',I4,
     #  6X,'K10BIS :',I4)                                           
5100    FORMAT(//21X,'PODREM CALCULAR :',I5,'  F.A.V. COMPLETES  ')
        RETURN                                                       
        END                                                           
C                                                                       
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                                                                   
        SUBROUTINE GFAV(N,NATS)                                  
        IMPLICIT REAL*8(A-H,O-Z)                                      
        COMMON/DD10/FAV(3,10000),V20(3),NPA10(100),NPTS10,K10,
     #       N10MAX,K10BIS,NA10,JA10,IA10,N10,I10,N10F,I10BIS
C                                                                  
        MA10 = 0.D0                                         
        DO L10=1,N10                                         
           IF (NPA10(L10).GE.1) MA10 = MA10 + 1               
        END DO                           
        MA10 = MA10 + N10F                
        IF (MA10.NE.0) THEN                
           DO J=1,NATS                      
              V20(J) = V20(J)/DFLOAT(N*MA10) 
           END DO                             
        END IF                                 
        DO IJ=1,NPTS10                          
           MA10 = 0                              
           DO L10=1,N10                           
              IF (NPA10(L10).GE.IJ) MA10 = MA10 + 1
           END DO                                          
           MA10 = MA10 + N10F                               
           IF (MA10.NE.0) THEN                               
              DO J=1,NATS                                     
                if (V20(J).ne.0.0) then 
                  FAV(J,IJ)=(FAV(J,IJ)/DFLOAT(N*MA10))/V20(J)   
                endif
              END DO                                            
           END IF                                       
        END DO                                           
c
           WRITE(10,*) NPTS10,V20(1)                                
           WRITE(10,*) ( FAV(1,IJ),IJ=1,NPTS10 )                     
c           WRITE(10,*) N10F,N10,(NPA10(L10),L10=1,N10)                
C                                                                     
           WRITE(60,1000) V20(1)                                        
           IF (NPTS10.LT.100) THEN                                      
              WRITE(60,1100)( FAV(1,IJ),IJ=1,NPTS10 )    
                ELSE                                     
              WRITE(60,1100)( FAV(1,IJ),IJ=1,100 )         
           END IF                                          
C                                                           
C             &&       nomes per a un io en aigua      &&    
C                                                             
           DO IJ=1,NPTS10                                      
              FAV(2,IJ) = ( FAV(2,IJ) + FAV(3,IJ) )/2.0D0       
           END DO                                        
           V20(2) = ( V20(2) + V20(3) )/2.0D0             
C                                                          
           WRITE(20,*) NPTS10,V20(2)                    
           WRITE(20,*) ( FAV(2,IJ),IJ=1,NPTS10 )         
c           WRITE(20,*) N10F,N10,(NPA10(L10),L10=1,N10)    
C                                                 
           WRITE(60,1001) V20(2)                    
           IF (NPTS10.LT.100) THEN                  
              WRITE(60,1100)( FAV(2,IJ), IJ=1,NPTS10 )
                ELSE                                 
              WRITE(60,1100)( FAV(2,IJ), IJ=1,100 )
           END IF                                     
C                                                       
1000    FORMAT(///25X, '**     ',                        
     &  'F . A . V.    D E L S    O X I G E N S    **'//  
     #  35X,'V2(0) = ',D15.6,' A2/s2'/)                    
1001    FORMAT(///25X, '**     ',                           
     &  'F . A . V.    D E L S    H I D R O G E N S    **'// 
     #  35X,'V2(0) = ',D15.6,' A2/s2'/)                       
1100    FORMAT(100(/' ',10F11.5))                              
        RETURN                                                  
        END
C                                                          
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
C                                                                    
C                           + + + + + + + + + + + + + +             
C                            ++      GEOMETRIA     ++              
C                           + + + + + + + + + + + + + +           
C                                                                
        SUBROUTINE IGEOM                                        
        IMPLICIT REAL*8(A-H,O-Z)                                     
        COMMON/DD45/DANGLE,IANMIN,IANMAX,PHOH(400),VMHOH,SHOH2,     
     #                                             VMOH,SOH2       
C                                                                 
        READ(5,*) DANGLE,IANMIN,IANMAX                           
        WRITE(60,1000) DANGLE,IANMIN,IANMAX                      
        PI = 4.0D0*DATAN(1.0D0)                                
        DANGLE = DANGLE*PI/DFLOAT(180)                        
        VMHOH = 0.D0                                         
        SHOH2 = 0.D0                                        
        VMOH = 0.D0                                        
        SOH2 = 0.D0                                       
        DO IAN=IANMIN,IANMAX                             
           PHOH(IAN) = 0.D0                            
        END DO                                        
1000    FORMAT(///10X,'******      ',                            
     #  'GEOMETRIA DE LES MOLECULES       ******'//             
     #  ' ','INTERVAL ANGLE HOH :',D10.2,' deg',5X,            
     #  'INDEX PRIMERA I DARRERA ZONES :',2I6/)               
        RETURN                                               
        END                                                 
C
C********************************************************************
C                                                         
        SUBROUTINE GGEOM(N,NPASS)                                       
        IMPLICIT REAL*8(A-H,O-Z)                                       
        COMMON/DD45/DANGLE,IANMIN,IANMAX,PHOH(400),VMHOH,SHOH2,       
     #                                             VMOH,SOH2         
C                                                                   
        PI = 4.0D0*DATAN(1.0D0)                                    
        DANGLE = DANGLE*DFLOAT(180)/PI                            
        NHOH = N*NPASS                                           
        NOH = 2*N*NPASS                                         
        VMHOH = (VMHOH/DFLOAT(NHOH))*DFLOAT(180)/PI            
        SHOH2 = SHOH2*(DFLOAT(180)/PI)**2.0D0                 
        SIGHOH = DSQRT((SHOH2 - DFLOAT(NHOH)*VMHOH*VMHOH)/   
     #                  DFLOAT(NHOH-1))                     
        VMOH = VMOH/DFLOAT(NOH)                            
        SIGOH = DSQRT((SOH2 - DFLOAT(NOH)*VMOH*VMOH)/     
     #                  DFLOAT(NOH-1))                   
        DO IAN=IANMIN,IANMAX                            
           PHOH(IAN) = PHOH(IAN)/DFLOAT(NHOH)          
        END DO                                        
C                                                                   
c        WRITE(15,*) DANGLE,IANMIN,IANMAX                           
c        WRITE(15,*) (PHOH(IAN),IAN=IANMIN,IANMAX)                 
C                                                                
        ANGMIN = DFLOAT(IANMIN-1)*DANGLE + DANGLE/2.0D0         
        ANGMAX = DFLOAT(IANMAX-1)*DANGLE + DANGLE/2.0D0        
        WRITE(60,1000) ANGMIN,ANGMAX                           
        WRITE(60,1010) (PHOH(IAN),IAN=IANMIN,IANMAX)          
        WRITE(60,1100) VMHOH,SIGHOH                          
        WRITE(60,1200) VMOH,SIGOH                           
1000    FORMAT(///10X,'*****     ',                                     
     #  'DISTRIBUCIO DE PROBABILITATS DE L''ANGLE H-O-H      *****'//  
     #  15X,'ANGLE MINIM :',D13.6,' deg',                             
     #   5X,'ANGLE MAXIM :',D13.6,' deg'/)                           
1010    FORMAT(100(/' ',10F11.5))                                   
1100    FORMAT(//15X,'VALOR MIG DE L''ANGLE H-O-H  :',D13.6,' deg'//    
     #  20X,'desv. standard :',D13.6)                                  
1200    FORMAT(//15X,'VALOR MIG DE LA DISTANCIA O-H  :',D13.6,' A'//  
     #  20X,'desv. standard :',D13.6)                                
        RETURN                                                      
        END                                                        
