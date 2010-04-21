C=====================================================================
C     Subroutine for finding Pd (Pressure Difference)
C     and gradient of cost function in the shape optimization of an airfoil
C     by solving 2-d Euler equations and corresponding adjoint equation 
C     Copyright: ZHANG Zhengke,TSL,NUS, Singapore, March 2003
C     Modified : Lum Kai Yew,TSL,NUS, Singapore, May 2003
C     
C     input: 
C     MACHIN: Mach number
C     AOAIN:  angle of attack
C     ALPHIN(ND)
C     output:
C     PDOUT:  1/2 Squared pressure differences
C     GOUT(ND): gradient of the cost function
C     RSDEU: residual of Euler equation at its last iteration
C     RSDADJ: residual of the adjoint equation at its last iteration
C==========================================================================
      SUBROUTINE PDGRADIN(MACHIN,AOAIN,ALPHIN,PDOUT,GOUT,RSDEU,RSDADJ)
      PARAMETER (ML=4,MI=161,MJ=91,MM=161)  
      PARAMETER (ND=24)	  
      REAL MACHIN,AOAIN
      DIMENSION ALPHIN(ND),GOUT(ND)
C     ------ variables,arrays for 2-d Euler eq. ----------      
      COMMON /UVV/ U(ML,MI,MJ),V(ML,MI,MJ)
      COMMON /TSP/ DTL(MI,MJ),SM(MI,MJ),PP(MI,MJ)
      COMMON /FQQ/ FDS(ML,MI,MJ),Q(ML,MI,MJ)
      COMMON /DDD/ D(ML,MI,MJ),AA(MM),BB(MM),CC(ML,MM)
      COMMON /DPP/ DP(MI,MJ),DPM(MI,MJ)
      COMMON /RAD/ RADI(MI,MJ),RADJ(MI,MJ)

      COMMON /GRD/ XX(MI,MJ),YY(MI,MJ)
      COMMON /GSFI/ SIX(MI,MJ),SIY(MI,MJ)
      COMMON /GSFJ/ SJX(MI,MJ),SJY(MI,MJ)
      COMMON /RXI/ XXI(MI,MJ),YXI(MI,MJ)
      COMMON /RET/ XET(MI,MJ),YET(MI,MJ)
      COMMON /FBNV/ UNX(MI),UNY(MI)
      COMMON /BDOWD/DSTRBD(MI)
      COMMON /CPBD/ CPBD(MI),XBD(MI),YBD(MI)
      COMMON /VOL/ VOL(MI,MJ)       
      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1
      
      COMMON /SWT/ INI,ITM,NST,NRSD,NFSV
      COMMON /E22/ E2I,E2J
      COMMON /E44/ E4I,E4J
      COMMON /RSD/ RSDM,ECRI,TSTA,COUR4,EPS,ITL
      COMMON /FAR/ GM,RMF,VXF,VYF,HF,AL
      COMMON /MAX/ IM,JM

C     ---- variables,arrays for 2-d adjoint eq. only---------
      COMMON /PSI/ PS0(ML,MI,MJ),PSI(ML,MI,MJ)
      COMMON /AT/ A1T(ML,ML,MI,MJ),A2T(ML,ML,MI,MJ)
      COMMON /PD/ PD(MI)
      COMMON /ADJFB/ TBT(ML,ML,MI)

C     ---- variables,arrays for gradient etc.---------
      COMMON /NEWG/ XNEW(MI,MJ),YNEW(MI,MJ)
      COMMON /DSV/ALPHA(ND),ALPHA0(ND)
      COMMON /DLI/DLI
	  COMMON /G/G(ND)
C     ------------------------------------------------
      COMMON /COST/COST
C     ------------------------------------------------
C     WRITE(*,*)
C     WRITE(*,*) '-----------------------------------'
C     WRITE(*,*) 'Begin Euler and adjoint solvers:'

      DO I=1,ND
         ALPHA0(I)=ALPHIN(I)
         ALPHA(I)=ALPHA0(I)
      END DO

C     ----- read in original parameters for Euler solver -----
      include 'eu.org'
C     ----- replace freestream conditions with arguments -----
      RMF=MACHIN
      AL=AOAIN

C     ----------------------------------- 
      AMPLF1=10.0
      AMPLF2=0.0      
      PI=3.141592653589793
      ALP=AL*PI/180.
C     ----- read in the grid data of initial shape -----
      OPEN(UNIT=3,NAME='afoil.grd',TYPE='OLD')
      DO J=1,NJ
         DO I=1,NI
            READ(3,*) XX(I,J),YY(I,J)
         END DO
	  END DO
      CLOSE(3)

C     ----- read in the desired(target) pressure -----
      OPEN(UNIT=2,NAME='pd.org',TYPE='OLD')
      DO I=1,NI-1
         READ(2,*) PD(I)
	  END DO
      CLOSE(2)

	  CALL NEWGRD
	  DO J=1,NJ
         DO I=1,NI
            XX(I,J)=XNEW(I,J)
            YY(I,J)=YNEW(I,J)
         END DO
      END DO
      
      CALL EUL2D
      RSDEU=RSDM
      ITLEU=ITL

      CALL COST_I 
      PDOUT = COST

      IF (RSDEU.GT.ECRI) THEN
         RSDADJ=1.0
         DO I=1,ND
            GOUT(I) = 0.0
         END DO
C     WRITE(*,*) '  Euler solve did not converged!'
C     WRITE(*,*) '  Skipping adjoint solver ...'
         GOTO 10
      END IF

      CALL ADJT2D
      RSDADJ=RSDM
      ITLADJ=ITL

	  CALL GRAD

      DO I=1,ND
         GOUT(I)=G(I)
      END DO

C 10   WRITE(*,20) PDOUT, RSDEU, RSDADJ
C 20   FORMAT(/3X,'Cost = ',F12.8,/3X,
C     1     'EULER residu   = ',F12.10,/3X,
C     2     'Adjoint residu = ',F12.10)
C     WRITE(*,*) 'End solvers'
C     WRITE(*,*) '-----------------------------------'
C     WRITE(*,*)

 10   RETURN
      END

C------------------Main sub-program -------------------------------
      SUBROUTINE EUL2D
C------------------------------------------------------------------
C     Numerical solution of 2-d Euler equations using finite volume	
C     method for compressible flows 
C-------------------------Program Eu2d.f -------------------------
C     Copyright: ZHANG Zhengke,TSL,NUS, Singapore, Jul.29, 2002
C     ----------  PROGRAM FLOW CHART ----------
C     Main sub-prog.-----------------------------
C     |          ( call Sub. ORGNL  )        
C     |          ( call Sub. FIELD  )        
C     Sub.ORGNL------------------------------ 
C     |     
C     |          ( call Sub. XYS   )
C     |
C     Sub. FIELD-----------------------------
C     |          ( call Sub. STEP   )
C     |          ( call Sub. RADIUS  )
C     |          ( call Sub. ARTVSC )
C     |          ( call Sub. FLUX   )
C     |          ( call Sub. EXPLIC )
C     |          ( call Sub. FLUX   )
C     |          ( call Sub. EXPLIC ) 
C     |          ( call Sub. FLUX   )
C     |          ( call Sub. SMOOTH )
C     |
C     Sub. XYS------------------------------
C     | 
C     Sub. STEP------------------------------
C     |
C     Sub. FLUX------------------------------ 
C     
C     Sub. BOUND-----------------------------
C     |
C     Sub. EXPLIC----------------------------
C     |
C     Sub. SMOOTH----------------------------
C     |
C     Sub. ARTVSC----------------------------
C     
C-------------------------------------------------------------------
      PARAMETER (ML=4,MI=161,MJ=91,MM=161)        
      COMMON /UVV/ U(ML,MI,MJ),V(ML,MI,MJ)
      COMMON /TSP/ DTL(MI,MJ),SM(MI,MJ),PP(MI,MJ)
      COMMON /FQQ/ FDS(ML,MI,MJ),Q(ML,MI,MJ)
      COMMON /DDD/ D(ML,MI,MJ),AA(MM),BB(MM),CC(ML,MM)
      COMMON /DPP/ DP(MI,MJ),DPM(MI,MJ)
      COMMON /RAD/ RADI(MI,MJ),RADJ(MI,MJ)

      COMMON /GRD/ XX(MI,MJ),YY(MI,MJ)
      COMMON /GSFI/ SIX(MI,MJ),SIY(MI,MJ)
      COMMON /GSFJ/ SJX(MI,MJ),SJY(MI,MJ)
      COMMON /RXI/ XXI(MI,MJ),YXI(MI,MJ)
      COMMON /RET/ XET(MI,MJ),YET(MI,MJ)
      COMMON /FBNV/ UNX(MI),UNY(MI)
      COMMON /BDOWD/DSTRBD(MI)
      COMMON /CPBD/ CPBD(MI),XBD(MI),YBD(MI)
      COMMON /VOL/ VOL(MI,MJ)       
      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1
      
      COMMON /SWT/ INI,ITM,NST,NRSD,NFSV
      COMMON /E22/ E2I,E2J
      COMMON /E44/ E4I,E4J
      COMMON /RSD/ RSDM,ECRI,TSTA,COUR4,EPS,ITL
      COMMON /FAR/ GM,RMF,VXF,VYF,HF,AL
      COMMON /MAX/ IM,JM


      CALL ORGNL
      CALL FIELD

      RETURN
      END
C------------------------------------------------------------------
      SUBROUTINE ORGNL
      PARAMETER (ML=4,MI=161,MJ=91,MM=161)              
      COMMON /UVV/ U(ML,MI,MJ),V(ML,MI,MJ)
      COMMON /TSP/ DTL(MI,MJ),SM(MI,MJ),PP(MI,MJ)
      COMMON /FQQ/ FDS(ML,MI,MJ),Q(ML,MI,MJ)
      COMMON /DDD/ D(ML,MI,MJ),AA(MM),BB(MM),CC(ML,MM)
      COMMON /DPP/ DP(MI,MJ),DPM(MI,MJ)
      COMMON /RAD/ RADI(MI,MJ),RADJ(MI,MJ)

      COMMON /GRD/ XX(MI,MJ),YY(MI,MJ)
      COMMON /GSFI/ SIX(MI,MJ),SIY(MI,MJ)
      COMMON /GSFJ/ SJX(MI,MJ),SJY(MI,MJ)
      COMMON /RXI/ XXI(MI,MJ),YXI(MI,MJ)
      COMMON /RET/ XET(MI,MJ),YET(MI,MJ)
      COMMON /FBNV/ UNX(MI),UNY(MI)
      COMMON /BDOWD/DSTRBD(MI)
      COMMON /CPBD/ CPBD(MI),XBD(MI),YBD(MI)              
      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1
      
      COMMON /SWT/ INI,ITM,NST,NRSD,NFSV
      COMMON /E22/ E2I,E2J
      COMMON /E44/ E4I,E4J
      COMMON /RSD/ RSDM,ECRI,TSTA,COUR4,EPS,ITL
      COMMON /FAR/ GM,RMF,VXF,VYF,HF,AL
      COMMON /MAX/ IM,JM


C     OPEN(UNIT=1,NAME='eu.org',TYPE='OLD')
C     READ(1,*) NI,NJ,IHF
C     READ(1,*) INI,ITM,NST,NRSD,NFSV
C     READ(1,*) E2I,E2J
C     READ(1,*) E4I,E4J
C     READ(1,*) COUR4,EPS
C     READ(1,*) RMF,AL,ECRI
C     CLOSE(1)

C     NI=129
C     NJ=65
C     IHF=65
C     INI=0
C     ITM=9000
C     NST=1
C     NRSD=500
C     NFSV=500
C     E2I=0.9
C     E2J=0.5
C     E4I=0.06
C     E4J=0.06
C     COUR4=0.7
C     EPS=0.333
C     RMF=0.8
C     AL=1.0
C     ECRI=0.00001

      NI1=NI-1
      NJ1=NJ-1


      GM=1.4
      RMF2=RMF**2
      PI=3.141592653589793
      ALP=AL*PI/180.
      VXF=COS(ALP)
      VYF=SIN(ALP)
      PPF=1./GM/RMF2
      ENF=PPF/(GM-1.)+0.5
      HF=ENF+PPF

      CALL XYS

      IF(INI.EQ.1) THEN
         OPEN(UNIT=2,NAME='eu.fld',TYPE='OLD')
         READ(2,*) (((U(L,I,J),L=1,ML)
     1        ,I=1,NI1),J=1,NJ1)
         CLOSE(2)
         DO J=1,NJ1
            DO I=1,NI1
               VV=U(2,I,J)**2+U(3,I,J)**2
               PP(I,J)=(GM-1.)*(U(4,I,J)-0.5*VV/U(1,I,J))
            END DO
	     END DO
      ELSE
         DO J=1,NJ1
            DO I=1,NI1
               U(1,I,J)=1.0
               U(2,I,J)=VXF
               U(3,I,J)=VYF
               U(4,I,J)=ENF
               PP(I,J)=PPF
            END DO
	     END DO
      END IF

      DO J=1,NJ1
         DO I=1,NI1
            DO L=1,ML
               V(L,I,J)=U(L,I,J)
            END DO
         END DO
	  END DO

      NM=MAX0(NI,NJ)
      BB(1)=1.+2.*EPS
      AA(1)=EPS/BB(1)
      DO I=2,NM
         BB(I)=BB(1)-EPS*EPS/BB(I-1)
         AA(I)=EPS/BB(I)
      END DO

      RETURN
      END
C-----------------------------------------------------------------
      SUBROUTINE FIELD
      PARAMETER (ML=4,MI=161,MJ=91,MM=161)                
      COMMON /UVV/ U(ML,MI,MJ),V(ML,MI,MJ)
      COMMON /TSP/ DTL(MI,MJ),SM(MI,MJ),PP(MI,MJ)
      
      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1

      COMMON /SWT/ INI,ITM,NST,NRSD,NFSV
      COMMON /RSD/ RSDM,ECRI,TSTA,COUR4,EPS,ITL
      COMMON /MAX/ IM,JM

C     WRITE(*,10)
 10   FORMAT(/3X,'EULER',5X,'ITE',9X,'RSDM',11X,'IM',3X,'JM')

      ITE=0
 20   CONTINUE
      CALL RADIUS
      K=MOD(ITE,NST)
      IF(K.EQ.0) CALL STEP
      CALL ARTVSC


      TSTA=1./4.
      CALL FLUX
      CALL EXPLIC
      TSTA=1./2.
      CALL FLUX
      CALL EXPLIC
      TSTA=1.
      CALL FLUX
      CALL SMOOTH

      ITE=ITE+1
      KMD=MOD(ITE,NFSV)
      IF(KMD.NE.0) GOTO 30
C     WRITE(*,23) ITE,RSDM,IM,JM
 23   FORMAT(11X,I5,4X,F12.8,5X,3I5,2X,'(Field saved)')
      GO TO 45

 30   CONTINUE
      KK=MOD(ITE,NRSD)
      IF(KK.NE.0) GOTO 45
C     WRITE(*,40) ITE,RSDM,IM,JM
 40   FORMAT(11X,I5,4X,F12.8,5X,3I5)

 45   CONTINUE

      IF(RSDM.LT.ECRI) GOTO 50
      IF(ITE.GE.ITM) GOTO 70
      GOTO 20

 50   CONTINUE
C     WRITE(*,60) ITE,ECRI
 60   FORMAT(6X,'ITE =',I5,5X,'RSDM <'
     1     ,1X,'ECRI =',F15.9,5X,'STOP')
      GOTO 90

 70   CONTINUE
C     WRITE(*,80) ITM,RSDM
 80   FORMAT(6X,'ITE >',1X,'ITEmax ='
     1     ,I5,5X,'RSDM =',F15.9,5X,'STOP')

 90   CONTINUE

      ITL=ITE
      RETURN
      END
C------------------------------------------------------------------
      SUBROUTINE XYS
      PARAMETER (ML=4,MI=161,MJ=91,MM=161)
      COMMON /TSP/ DTL(MI,MJ),SM(MI,MJ),PP(MI,MJ)
      
      COMMON /GRD/ XX(MI,MJ),YY(MI,MJ)
      COMMON /GSFI/ SIX(MI,MJ),SIY(MI,MJ)
      COMMON /GSFJ/ SJX(MI,MJ),SJY(MI,MJ)
      COMMON /RXI/ XXI(MI,MJ),YXI(MI,MJ)
      COMMON /RET/ XET(MI,MJ),YET(MI,MJ)
      COMMON /FBNV/ UNX(MI),UNY(MI)
      COMMON /BDOWD/DSTRBD(MI)
      COMMON /CPBD/ CPBD(MI),XBD(MI),YBD(MI) 
      COMMON /VOL/ VOL(MI,MJ)  

      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1


      DO I=1,NI1
         I1=I+1
         XBD(I)=(XX(I,1)+XX(I1,1))/2.
         YBD(I)=(YY(I,1)+YY(I1,1))/2.
	  END DO

      DO I=1,NI
         DO J=1,NJ1
            J1=J+1
            RETX=XX(I,J1)-XX(I,J)
            RETY=YY(I,J1)-YY(I,J)
            SIX(I,J)=RETY
            SIY(I,J)=-RETX
         END DO
      END DO

      DO J=1,NJ
         DO I=1,NI1
            I1=I+1
            RXIX=XX(I1,J)-XX(I,J)
            RXIY=YY(I1,J)-YY(I,J)
            SJX(I,J)=-RXIY
            SJY(I,J)=RXIX
         END DO
      END DO

      DO I=1,NI1
         ARE=SQRT(SJX(I,NJ)**2+SJY(I,NJ)**2)
         IF(ARE.NE.0.) THEN
            UNX(I)=SJX(I,NJ)/ARE
            UNY(I)=SJY(I,NJ)/ARE
         ELSE 
            UNX(I)=0.
            UNY(I)=0.
         END IF
      END DO

      J=2
      J1=J-1
      J2=J+1
      DO I=1,NI1
         I1=I+1
         X1=(XX(I,J1)+XX(I1,J1))/2.
         Y1=(YY(I,J1)+YY(I1,J1))/2.
         X2=(XX(I,J)+XX(I1,J))/2.
         Y2=(YY(I,J)+YY(I1,J))/2.
         DST1=SQRT((X2-X1)**2+(Y2-Y1)**2)
         X3=(XX(I,J2)+XX(I1,J2))/2.
         Y3=(YY(I,J2)+YY(I1,J2))/2.
         DST2=SQRT((X3-X2)**2+(Y3-Y2)**2)
         DSTRBD(I)=DST1/(DST1+DST2)
      END DO
      
      DO I=1,NI1
         I1=I+1
         DO J=1,NJ1
            J1=J+1
            SS1=SQRT(SIX(I,J)**2+SIY(I,J)**2)   
            SS2=SQRT(SJX(I,J)**2+SJY(I,J)**2)
            SS3=SQRT(SIX(I1,J)**2+SIY(I1,J)**2)
            SS4=SQRT(SJX(I,J1)**2+SJY(I,J1)**2)
            SM(I,J)=AMAX1(SS1,SS2,SS3,SS4)
         END DO
      END DO

      DO J=1,NJ1
         J1=J+1
         DO I=1,NI1
            I1=I+1
            XXI(I,J)=(XX(I1,J)+XX(I1,J1))/2.-(XX(I,J)+XX(I,J1))/2.
            YXI(I,J)=(YY(I1,J)+YY(I1,J1))/2.-(YY(I,J)+YY(I,J1))/2. 
            XET(I,J)=(XX(I,J1)+XX(I1,J1))/2.-(XX(I,J)+XX(I1,J))/2.
            YET(I,J)=(YY(I,J1)+YY(I1,J1))/2.-(YY(I,J)+YY(I1,J))/2.
         END DO
      END DO

      DO J=1,NJ1
         J1=J+1
         DO I=1,NI1
            I1=I+1
            R1X=XX(I1,J1)-XX(I,J)
            R1Y=YY(I1,J1)-YY(I,J)
            R2X=XX(I,J1)-XX(I1,J)
            R2Y=YY(I,J1)-YY(I1,J)
            VOL(I,J)=ABS(R1X*R2Y-R1Y*R2X)/2.
         END DO
      END DO

      RETURN
      END
C------------------------------------------------------------------
      SUBROUTINE  RADIUS
      PARAMETER (ML=4,MI=161,MJ=91,MM=161)        
      COMMON /UVV/ U(ML,MI,MJ),V(ML,MI,MJ)
      COMMON /TSP/ DTL(MI,MJ),SM(MI,MJ),PP(MI,MJ)
      COMMON /RAD/ RADI(MI,MJ),RADJ(MI,MJ)

      COMMON /RXI/ XXI(MI,MJ),YXI(MI,MJ)
      COMMON /RET/ XET(MI,MJ),YET(MI,MJ)     
      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1
      COMMON /FAR/ GM,RMF,VXF,VYF,HF,AL


      DO J=1,NJ1
         DO I=1,NI1
            VX=V(2,I,J)/V(1,I,J)
            VY=V(3,I,J)/V(1,I,J)
            SPSD=SQRT(GM*PP(I,J)/V(1,I,J))
            RADI(I,J)=ABS(VX*YET(I,J)-VY*XET(I,J))
     &           +SPSD*SQRT(XET(I,J)**2+YET(I,J)**2)
            RADJ(I,J)=ABS(VY*XXI(I,J)-VX*YXI(I,J))
     &           +SPSD*SQRT(XXI(I,J)**2+YXI(I,J)**2)
         END DO
      END DO

      RETURN
      END
C-----------------------------------------------------------------
      SUBROUTINE STEP
      PARAMETER (ML=4,MI=161,MJ=91,MM=161)               
      COMMON /UVV/ U(ML,MI,MJ),V(ML,MI,MJ)
      COMMON /TSP/ DTL(MI,MJ),SM(MI,MJ),PP(MI,MJ)
      COMMON /RAD/ RADI(MI,MJ),RADJ(MI,MJ)

      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1

      COMMON /FAR/ GM,RMF,VXF,VYF,HF,AL
      COMMON /RSD/ RSDM,ECRI,TSTA,COUR4,EPS,ITL

      DO J=1,NJ1
         DO I=1,NI1
            RR=V(1,I,J)
            VX=V(2,I,J)/RR
            VY=V(3,I,J)/RR
            AA=SQRT(GM*PP(I,J)/RR)
            QQ=SQRT(VX**2+VY**2)
            DTL(I,J)=COUR4/AMAX1(QQ,AA)/SM(I,J)
c---- DTL(I,J)=COUR4/(RADI(I,J)+RADJ(I,J))
         END DO
	  END DO

      RETURN
      END
C------------------------------------------------------------------
      SUBROUTINE FLUX
      PARAMETER (ML=4,MI=161,MJ=91,MM=161)               
      COMMON /UVV/ U(ML,MI,MJ),V(ML,MI,MJ)
      COMMON /TSP/ DTL(MI,MJ),SM(MI,MJ),PP(MI,MJ)
      COMMON /FQQ/ FDS(ML,MI,MJ),Q(ML,MI,MJ)

      COMMON /GSFI/ SIX(MI,MJ),SIY(MI,MJ)
      COMMON /GSFJ/ SJX(MI,MJ),SJY(MI,MJ)
      COMMON /FBNV/ UNX(MI),UNY(MI)

      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1
      
      COMMON /FAR/ GM,RMF,VXF,VYF,HF,AL
      
      CALL BOUND(1)
      DO 10 I=2,NI1
         I1=I-1
         DO 10 J=1,NJ1
            P2=(PP(I,J)+PP(I1,J))/2.
            W1=(V(1,I,J)+V(1,I1,J))/2.
            W2=(V(2,I,J)+V(2,I1,J))/2.
            W3=(V(3,I,J)+V(3,I1,J))/2.
            W4=(V(4,I,J)+V(4,I1,J))/2.
            FDS(1,I,J)=SIX(I,J)*W2+SIY(I,J)*W3
            FDS(2,I,J)=FDS(1,I,J)*W2/W1+P2*SIX(I,J)
            FDS(3,I,J)=FDS(1,I,J)*W3/W1+P2*SIY(I,J)
            FDS(4,I,J)=FDS(1,I,J)*((W4+P2)/W1+HF)/2.
 10      CONTINUE

         DO 20 I=1,NI1
            I1=I+1
            DO 20 J=1,NJ1
               DO 20 L=1,ML
                  Q(L,I,J)=FDS(L,I1,J)-FDS(L,I,J)
 20            CONTINUE

               
               CALL BOUND(2)
               DO 60 J=2,NJ1
                  J1=J-1
                  DO 40 I=1,NI1
                     P2=(PP(I,J)+PP(I,J1))/2.
                     W1=(V(1,I,J)+V(1,I,J1))/2.
                     W2=(V(2,I,J)+V(2,I,J1))/2.
                     W3=(V(3,I,J)+V(3,I,J1))/2.
                     W4=(V(4,I,J)+V(4,I,J1))/2.
                     FDS(1,I,J)=SJX(I,J)*W2+SJY(I,J)*W3
                     FDS(2,I,J)=FDS(1,I,J)*W2/W1+P2*SJX(I,J)
                     FDS(3,I,J)=FDS(1,I,J)*W3/W1+P2*SJY(I,J)
                     FDS(4,I,J)=FDS(1,I,J)*((W4+P2)/W1+HF)/2.
 40               CONTINUE
 60            CONTINUE

               DO 80 J=1,NJ1
                  J1=J+1
                  DO 80 I=1,NI1
                     DO 80 L=1,ML
                        Q(L,I,J)=Q(L,I,J)+FDS(L,I,J1)-FDS(L,I,J)
 80                  CONTINUE

                     RETURN
                     END
C------------------------------------------------------------------
      SUBROUTINE BOUND(II)
      PARAMETER (ML=4,MI=161,MJ=91,MM=161)               
      COMMON /UVV/ U(ML,MI,MJ),V(ML,MI,MJ)
      COMMON /TSP/ DTL(MI,MJ),SM(MI,MJ),PP(MI,MJ)
      COMMON /FQQ/ FDS(ML,MI,MJ),Q(ML,MI,MJ)

      COMMON /GSFI/ SIX(MI,MJ),SIY(MI,MJ)
      COMMON /GSFJ/ SJX(MI,MJ),SJY(MI,MJ)
      COMMON /FBNV/ UNX(MI),UNY(MI)
      COMMON /BDOWD/DSTRBD(MI)

      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1 
      
      COMMON /FAR/ GM,RMF,VXF,VYF,HF,AL 	          

      IF(II.EQ.2) GOTO 400

      DO J=1,NJ1 
         PP2=(PP(1,J)+PP(NI1,J))/2.
         W1=(V(1,1,J)+V(1,NI1,J))/2.
         W2=(V(2,1,J)+V(2,NI1,J))/2.
         W3=(V(3,1,J)+V(3,NI1,J))/2.
         W4=(V(4,1,J)+V(4,NI1,J))/2.
         FDS(1,1,J)=SIX(1,J)*W2+SIY(1,J)*W3
         FDS(2,1,J)=FDS(1,1,J)*W2/W1+PP2*SIX(1,J)
         FDS(3,1,J)=FDS(1,1,J)*W3/W1+PP2*SIY(1,J)
         FDS(4,1,J)=FDS(1,1,J)*((W4+PP2)/W1+HF)/2.
         FDS(1,NI,J)=FDS(1,1,J)
         FDS(2,NI,J)=FDS(2,1,J)
         FDS(3,NI,J)=FDS(3,1,J)
         FDS(4,NI,J)=FDS(4,1,J)
      END DO
      RETURN

 400  CONTINUE
      DO 460 I=1,NI1
         RR=V(1,I,NJ1)
         VX=V(2,I,NJ1)/RR
         VY=V(3,I,NJ1)/RR
         EN=V(4,I,NJ1)
         P2=PP(I,NJ1)
         AE=SQRT(GM*P2/RR)
         QNE=VX*UNX(I)+VY*UNY(I)
C     ---- RME=SQRT(VX**2+VY**2)/AE

         AF=1./RMF
         QNF=VXF*UNX(I)+VYF*UNY(I)

         QNSB=0.5*(QNE+QNF)+(AE-AF)/(GM-1.0)
         ASB=(GM-1.)*(QNE-QNF)/4.0+0.5*(AE+AF)

         IF(ABS(QNSB).GE.ASB) GOTO 420

         IF(QNSB.GT.0.) GOTO 410
         SON=1./GM/RMF**2
         QTX=VXF-QNF*UNX(I)
         QTY=VYF-QNF*UNY(I)
         VXFB=QTX+QNSB*UNX(I)
         VYFB=QTY+QNSB*UNY(I)
         RRFB=(ASB**2/SON/GM)**(1./(GM-1.))
         PPFB=ASB**2*RRFB/GM
         GOTO 440
 410     CONTINUE
         SON=P2/RR**GM
         QTX=VX-QNE*UNX(I)
         QTY=VY-QNE*UNY(I)
         VXFB=QTX+QNSB*UNX(I)
         VYFB=QTY+QNSB*UNY(I)
         RRFB=(ASB**2/SON/GM)**(1./(GM-1.))
         PPFB=ASB**2*RRFB/GM
         GOTO 440

 420     CONTINUE
         QNSPI=QNF
         ASPI=AF
         QNSPO=QNE
         ASPO=AE
         IF (QNSPI.GT.0..AND.QNSPO.GT.0.) GO TO 430
         VXFB=VXF
         VYFB=VYF
         RRFB=1.0
         PPFB=1./GM/RMF**2
         GOTO 440
 430     CONTINUE
         VXFB=VX
         VYFB=VY
         RRFB=RR
         PPFB=P2


 440     CONTINUE
         VV=RRFB*(VXFB**2+VYFB**2)
         W1=RRFB
         W2=RRFB*VXFB
         W3=RRFB*VYFB
         W4=PPFB/(GM-1.)+0.5*VV
         PP2=PPFB
         FDS(1,I,NJ)=SJX(I,NJ)*W2+SJY(I,NJ)*W3
         FDS(2,I,NJ)=FDS(1,I,NJ)*W2/W1+PP2*SJX(I,NJ)
         FDS(3,I,NJ)=FDS(1,I,NJ)*W3/W1+PP2*SJY(I,NJ)
         FDS(4,I,NJ)=FDS(1,I,NJ)*((W4+PP2)/W1+HF)/2.
 460  CONTINUE
      
      DO 480 I=1,NI1
         AA=1.+DSTRBD(I)
         BB=-DSTRBD(I)
         PPB=AA*PP(I,1)+BB*PP(I,2)
         FDS(1,I,1)=0.0
         FDS(2,I,1)=PPB*SJX(I,1)
         FDS(3,I,1)=PPB*SJY(I,1)
         FDS(4,I,1)=0.0
 480  CONTINUE
      RETURN
      END
C------------------------------------------------------------------
      SUBROUTINE EXPLIC
      PARAMETER (ML=4,MI=161,MJ=91,MM=161)               
      COMMON /UVV/ U(ML,MI,MJ),V(ML,MI,MJ)
      COMMON /TSP/ DTL(MI,MJ),SM(MI,MJ),PP(MI,MJ)
      COMMON /FQQ/ FDS(ML,MI,MJ),Q(ML,MI,MJ)
      COMMON /DDD/ D(ML,MI,MJ),AA(MM),BB(MM),CC(ML,MM)

      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1

      COMMON /RSD/ RSDM,ECRI,TSTA,COUR4,EPS,ITL
      COMMON /FAR/ GM,RMF,VXF,VYF,HF,AL       

      DO 20 J=1,NJ1
         DO 20 I=1,NI1
            DO 10 L=1,ML
               V(L,I,J)=U(L,I,J)-TSTA*DTL(I,J)*(Q(L,I,J)
     1              -D(L,I,J))
 10         CONTINUE
            VV=V(2,I,J)**2+V(3,I,J)**2
            PP(I,J)=(GM-1.)*(V(4,I,J)-0.5*VV/V(1,I,J))
 20      CONTINUE

         RETURN
         END
C-----------------------------------------------------------------
      SUBROUTINE SMOOTH
      PARAMETER (ML=4,MI=161,MJ=91,MM=161)               
      COMMON /UVV/ U(ML,MI,MJ),V(ML,MI,MJ)
      COMMON /TSP/ DTL(MI,MJ),SM(MI,MJ),PP(MI,MJ)
      COMMON /FQQ/ FDS(ML,MI,MJ),Q(ML,MI,MJ)
      COMMON /DDD/ D(ML,MI,MJ),AA(MM),BB(MM),CC(ML,MM)

      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1

      COMMON /RSD/ RSDM,ECRI,TSTA,COUR4,EPS,ITL
      COMMON /FAR/ GM,RMF,VXF,VYF,HF,AL       
      COMMON /MAX/ IM,JM

      RSDM=0.0

      DO 30 J=1,NJ1
         DO 20 I=1,NI1
            DO 10 L=1,ML
               Q(L,I,J)=-TSTA*DTL(I,J)*(Q(L,I,J)-D(L,I,J))
 10         CONTINUE
            RSD=ABS(Q(4,I,J)/V(4,I,J))
            IF(RSD.LT.RSDM) GOTO 20
            RSDM=RSD
            IM=I
            JM=J
 20      CONTINUE
 30   CONTINUE

      EPS1=EPS+1.0
C     ---------------- I direction ----------------------------------
      DO 110 J=1,NJ1
         DO 50 L=1,ML
            CC(L,1)=EPS1*Q(L,1,J)
            Q(L,NI1,J)=EPS1*Q(L,NI1,J)
 50      CONTINUE
         DO 70 I=2,NI1
            I1=I-1
            DO 60 L=1,ML
               CC(L,I)=Q(L,I,J)+AA(I1)*CC(L,I1)
 60         CONTINUE
 70      CONTINUE
         DO 80 L=1,ML
            Q(L,NI1,J)=CC(L,NI1)/BB(NI1)
 80      CONTINUE
         DO 100 I1=1,NI1-1
            I=NI1-I1
            I2=I+1
            DO 90 L=1,ML
               Q(L,I,J)=CC(L,I)/BB(I)+AA(I)*Q(L,I2,J)
 90         CONTINUE
 100     CONTINUE
 110  CONTINUE

C     ----------------- J direction --------------------------------
      DO 270 I=1,NI1
         DO 210 L=1,ML
            CC(L,1)=EPS1*Q(L,I,1)
            Q(L,I,NJ1)=EPS1*Q(L,I,NJ1)
 210     CONTINUE
         DO 230 J=2,NJ1
            J1=J-1
            DO 220 L=1,ML
               CC(L,J)=Q(L,I,J)+AA(J1)*CC(L,J1)
 220        CONTINUE
 230     CONTINUE
         DO 240 L=1,ML
            Q(L,I,NJ1)=CC(L,NJ1)/BB(NJ1)
 240     CONTINUE
         DO 260 J1=1,NJ1-1
            J=NJ1-J1
            J2=J+1
            DO 250 L=1,ML
               Q(L,I,J)=CC(L,J)/BB(J)+AA(J)*Q(L,I,J2)
 250        CONTINUE
 260     CONTINUE
 270  CONTINUE

C     -----------------------------------------------------------

      DO 810 J=1,NJ1
         DO 810 I=1,NI1
            DO 800 L=1,ML
               V(L,I,J)=U(L,I,J)+Q(L,I,J)
               U(L,I,J)=V(L,I,J)
 800        CONTINUE
            VV=V(2,I,J)**2+V(3,I,J)**2
            PP(I,J)=(GM-1.)*(V(4,I,J)-0.5*VV/V(1,I,J))
 810     CONTINUE

         RETURN
         END

C-------------------------------------------------------------
      SUBROUTINE  ARTVSC
      PARAMETER (ML=4,MI=161,MJ=91,MM=161)        
      COMMON /UVV/ U(ML,MI,MJ),V(ML,MI,MJ)
      COMMON /TSP/ DTL(MI,MJ),SM(MI,MJ),PP(MI,MJ)
      COMMON /FQQ/ FDS(ML,MI,MJ),Q(ML,MI,MJ)
      COMMON /DDD/ D(ML,MI,MJ),AA(MM),BB(MM),CC(ML,MM)
      COMMON /DPP/ DP(MI,MJ),DPM(MI,MJ)
      COMMON /RAD/ RADI(MI,MJ),RADJ(MI,MJ)

      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1

      COMMON /E22/ E2I,E2J
      COMMON /E44/ E4I,E4J
      COMMON /FAR/ GM,RMF,VXF,VYF,HF,AL
      
      NI2=NI1-1
      NJ2=NJ1-1
      NI3=NI2-1
      NJ3=NJ2-1
C     omg=1.0 for transonic flow
C     omg=0.5 for supersonic or hypersonic
	  OMG=0.5

C     ---------- i direction ----------------------------
      DO J=1,NJ1
         DO I=2,NI2
            DP(I,J)=ABS(PP(I+1,J)-2.*PP(I,J)+PP(I-1,J))/
     1           ((1.-OMG)*(ABS(PP(I+1,J)-PP(I,J))
     2           +ABS(PP(I,J)-PP(I-1,J)))
     3           +OMG*(PP(I+1,J)+2.*PP(I,J)+PP(I-1,J)))
         END DO
	  END DO

      DO J=1,NJ1
         DO I=3,NI3
            DPM(I,J)=AMAX1(DP(I-1,J),DP(I,J),DP(I+1,J))
         END DO
         DPM(2,J)=DPM(3,J)
         DPM(1,J)=DPM(2,J)
         DPM(NI2,J)=DPM(NI3,J)
         DPM(NI1,J)=DPM(NI2,J)
	  END DO

      DO 30 L=1,ML
         DO 30 J=1,NJ1
            D(L,1,J)=0.0
            D(L,NI1,J)=0.0
 30      CONTINUE

         DO 80 I=1,NI2
            I1=I-1
            I2=I+1
            I3=I+2
            DO 60 J=1,NJ1
               RADN=(RADI(I,J)+RADI(I+1,J))/2.
               EP2=E2I*DPM(I,J)
               DO 40 L=1,ML
                  FDS(L,I,J)=RADN*EP2*(V(L,I2,J)-V(L,I,J))
 40            CONTINUE
               IF(I.EQ.1.OR.I.EQ.NI2) GOTO 60
               EP4=AMAX1(0.0,E4I-EP2)
               DO 50 L=1,ML
                  FDS(L,I,J)=FDS(L,I,J)-RADN*EP4*
     1                 (V(L,I3,J)-3.*V(L,I2,J)
     2                 +3.*V(L,I,J)-V(L,I1,J))
 50            CONTINUE
 60         CONTINUE
 80      CONTINUE

         DO 90 I=2,NI2
            I1=I-1
            DO 90 J=1,NJ1
               DO 90 L=1,ML
                  D(L,I,J)=FDS(L,I,J)-FDS(L,I1,J)
 90            CONTINUE

C     ---------- j direction ----------------------------
               DO I=1,NI1
                  DO J=2,NJ2
                     DP(I,J)=ABS(PP(I,J+1)-2.*PP(I,J)+PP(I,J-1))/
     1                    ((1.-OMG)*(ABS(PP(I,J+1)-PP(I,J))
     2                    +ABS(PP(I,J)-PP(I,J-1)))
     3                    +OMG*(PP(I,J+1)+2.*PP(I,J)+PP(I,J-1)))
                  END DO
               END DO

               DO I=1,NI1
                  DO J=3,NJ3
                     DPM(I,J)=AMAX1(DP(I,J-1),DP(I,J),DP(I,J+1))
                  END DO
                  DPM(I,2)=DPM(I,3)
                  DPM(I,1)=DPM(I,2)
                  DPM(I,NJ2)=DPM(I,NJ3)
                  DPM(I,NJ1)=DPM(I,NJ2)
               END DO

               DO 180 J=1,NJ2
                  J1=J-1
                  J2=J+1
                  J3=J+2
                  DO 160 I=1,NI1
                     RADN=(RADJ(I,J)+RADJ(I,J+1))/2.
                     EP2=E2J*DPM(I,J)
                     DO 140 L=1,ML
                        FDS(L,I,J)=RADN*EP2*(V(L,I,J2)-V(L,I,J))
 140                 CONTINUE
                     IF(J.EQ.1.OR.J.EQ.NJ2) GOTO 160
                     EP4=AMAX1(0.0,E4J-EP2)
                     DO 150 L=1,ML
                        FDS(L,I,J)=FDS(L,I,J)-RADN*EP4*
     1                       (V(L,I,J3)-3.*V(L,I,J2)
     2                       +3.*V(L,I,J)-V(L,I,J1))
 150                 CONTINUE
 160              CONTINUE
 180           CONTINUE

               DO 220 J=2,NJ2
                  J1=J-1
                  DO 210 I=1,NI1
                     DO 200 L=1,ML
                        D(L,I,J)=D(L,I,J)+FDS(L,I,J)-FDS(L,I,J1)
 200                 CONTINUE
 210              CONTINUE
 220           CONTINUE

               RETURN
               END
C--------------------------------------------------------------------
C------the end 2-d Euler solver ------------------------------------
C====================================================================
C-------------------------------------------------------------------
C-------------------------------------------------------------------
      SUBROUTINE COST_I
      PARAMETER (ML=4,MI=161,MJ=91,MM=161) 
      COMMON /TSP/ DTL(MI,MJ),SM(MI,MJ),PP(MI,MJ)

      COMMON /BDOWD/DSTRBD(MI)

      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1 

      COMMON /PD/ PD(MI)

      COMMON /COST/COST

	  COST=0.0
      DO I=1,NI1
         AA=1.+DSTRBD(I)
         BB=-DSTRBD(I)
         PPB=AA*PP(I,1)+BB*PP(I,2)
         COST=COST+(PPB-PD(I))*(PPB-PD(I))
      END DO
      COST=COST/2.0

      RETURN
      END
C------------------------------------------------------------------
C     
C===================================================================
C-----------------Main sub-program --------------------------------
      SUBROUTINE ADJT2D
C------------------------------------------------------------------
C     Solution of 2-d adjoint equation by discretizing the integration 
C     form of the equation directly in physical spacec using finite 
C     volume method.  The adjoint equation corresponds to the 
C     ompressible-flow Euler equations.      
C-----------------------Program adjt2d.f -------------------------
C     Copyright: ZHANG Zhengke,TSL,NUS, Singapore, Aug., 2002
C     ----------  PROGRAM FLOW CHART ----------
C     Main Sub-prog.-----------------------------
C     |          ( call Sub. ADJ_ORG  )        
C     |          ( call Sub. ADJ_FIELD  )        
C     |          ( call Sub. ADJ_OUTPUT ) 
C     Sub.ADJ_ORG------------------------------ 
C     |     
C     |          ( call Sub. A1TA2T )
C     |          ( call Sub. RADIUS )
C     |          ( call Sub. STEP )
C     |
C     Sub. ADJ_FIELD-----------------------------
C     |          ( call Sub. ADJ_ARTV )
C     |          ( call Sub. ADJ_FLUX )
C     |          ( call Sub. ADJ_EXPLIC )
C     |          ( call Sub. ADJ_FLUX )
C     |          ( call Sub. ADJ_EXPLIC ) 
C     |          ( call Sub. ADJ_FLUX )
C     |          ( call Sub. ADJ_SMTH )
C     Sub. ADJ_OUTPUT----------------------------
C     |          
C     Sub. A1TA2T
C     |          
C     Sub. ADJ_FLUX------------------------------ 
C     |
C     Sub. ADJ_BOUND-----------------------------
C     |
C     Sub. ADJ_EXPLIC----------------------------
C     |
C     Sub. ADJ_SMTH----------------------------
C     |
C     Sub. ADJ_ARTV----------------------------
C     
C--------------------------------------------------------------------
      PARAMETER (ML=4,MI=161,MJ=91,MM=161)
      PARAMETER (ND=24)
      COMMON /UVV/ U(ML,MI,MJ),V(ML,MI,MJ)
      COMMON /TSP/ DTL(MI,MJ),SM(MI,MJ),PP(MI,MJ)
      COMMON /FQQ/ FDS(ML,MI,MJ),Q(ML,MI,MJ)
      COMMON /DDD/ D(ML,MI,MJ),AA(MM),BB(MM),CC(ML,MM)
      COMMON /DPP/ DP(MI,MJ),DPM(MI,MJ)
      COMMON /RAD/ RADI(MI,MJ),RADJ(MI,MJ)

      COMMON /GRD/ XX(MI,MJ),YY(MI,MJ)
      COMMON /GSFI/ SIX(MI,MJ),SIY(MI,MJ)
      COMMON /GSFJ/ SJX(MI,MJ),SJY(MI,MJ)
      COMMON /RXI/ XXI(MI,MJ),YXI(MI,MJ)
      COMMON /RET/ XET(MI,MJ),YET(MI,MJ)
      COMMON /FBNV/ UNX(MI),UNY(MI)
      COMMON /BDOWD/DSTRBD(MI)
      COMMON /CPBD/ CPBD(MI),XBD(MI),YBD(MI)
      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1
      
      COMMON /SWT/ INI,ITM,NST,NRSD,NFSV
      COMMON /E22/ E2I,E2J
      COMMON /E44/ E4I,E4J
      COMMON /RSD/ RSDM,ECRI,TSTA,COUR4,EPS,ITL
      COMMON /FAR/ GM,RMF,VXF,VYF,HF,AL
      COMMON /MAX/ IM,JM

      COMMON /PSI/ PS0(ML,MI,MJ),PSI(ML,MI,MJ)
      COMMON /AT/ A1T(ML,ML,MI,MJ),A2T(ML,ML,MI,MJ)
      COMMON /PD/ PD(MI)
      COMMON /ADJFB/ TBT(ML,ML,MI)


      CALL ADJ_ORG
      CALL ADJ_FIELD

      RETURN
      END
C------------------------------------------------------------------
      SUBROUTINE ADJ_ORG
      PARAMETER (ML=4,MI=161,MJ=91,MM=161)
      COMMON /UVV/ U(ML,MI,MJ),V(ML,MI,MJ)
      COMMON /TSP/ DTL(MI,MJ),SM(MI,MJ),PP(MI,MJ)
      COMMON /FQQ/ FDS(ML,MI,MJ),Q(ML,MI,MJ)
      COMMON /DDD/ D(ML,MI,MJ),AA(MM),BB(MM),CC(ML,MM)
      COMMON /DPP/ DP(MI,MJ),DPM(MI,MJ)
      COMMON /RAD/ RADI(MI,MJ),RADJ(MI,MJ) 

      COMMON /GRD/ XX(MI,MJ),YY(MI,MJ)
      COMMON /GSFI/ SIX(MI,MJ),SIY(MI,MJ)
      COMMON /GSFJ/ SJX(MI,MJ),SJY(MI,MJ)
      COMMON /RXI/ XXI(MI,MJ),YXI(MI,MJ)
      COMMON /RET/ XET(MI,MJ),YET(MI,MJ)
      COMMON /FBNV/ UNX(MI),UNY(MI)
      COMMON /BDOWD/DSTRBD(MI)
      COMMON /CPBD/ CPBD(MI),XBD(MI),YBD(MI)
      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1
      
      COMMON /SWT/ INI,ITM,NST,NRSD,NFSV
      COMMON /E22/ E2I,E2J
      COMMON /E44/ E4I,E4J
      COMMON /RSD/ RSDM,ECRI,TSTA,COUR4,EPS,ITL
      COMMON /FAR/ GM,RMF,VXF,VYF,HF,AL
      COMMON /MAX/ IM,JM

      COMMON /PSI/ PS0(ML,MI,MJ),PSI(ML,MI,MJ)
      COMMON /AT/ A1T(ML,ML,MI,MJ),A2T(ML,ML,MI,MJ)
      COMMON /PD/ PD(MI)
      COMMON /ADJFB/ TBT(ML,ML,MI)


C     OPEN(UNIT=1,NAME='adj.org',TYPE='OLD')
C     READ(1,*) NI,NJ,IHF
C     READ(1,*) INI,ITM,NST,NRSD,NFSV
C     READ(1,*) E2I,E2J
C     READ(1,*) E4I,E4J
C     READ(1,*) COUR4,EPS
C     READ(1,*) ECRI
C     CLOSE(1)

      include 'adj.org'

      IF(INI.EQ.1) THEN
         OPEN(UNIT=2,NAME='adj.fld',TYPE='OLD')
         READ(2,*) (((PS0(L,I,J),L=1,ML)
     1        ,I=1,NI1),J=1,NJ1)
         CLOSE(2)
      ELSE
         DO J=1,NJ1
            DO I=1,NI1
               PS0(1,I,J)=0.0
               PS0(2,I,J)=0.0
               PS0(3,I,J)=0.0
               PS0(4,I,J)=0.0
            END DO
         END DO
      END IF

      DO J=1,NJ1
         DO I=1,NI1
            DO L=1,ML
               PSI(L,I,J)=PS0(L,I,J)
            END DO
         END DO
	  END DO


      CALL A1TA2T
      CALL RADIUS
      CALL STEP

      NM=MAX0(NI,NJ)
      BB(1)=1.+2.*EPS
      AA(1)=EPS/BB(1)
      DO I=2,NM
         BB(I)=BB(1)-EPS*EPS/BB(I-1)
         AA(I)=EPS/BB(I)
      END DO

      RETURN
      END
C-----------------------------------------------------------------
      SUBROUTINE ADJ_FIELD
      PARAMETER (ML=4,MI=161,MJ=91,MM=161) 
      COMMON /PSI/ PS0(ML,MI,MJ),PSI(ML,MI,MJ) 
      
      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1

      COMMON /SWT/ INI,ITM,NST,NRSD,NFSV
      COMMON /RSD/ RSDM,ECRI,TSTA,COUR4,EPS,ITL
      COMMON /MAX/ IM,JM

C     WRITE(*,10)
 10   FORMAT(/3X,'ADJOINT',3X,'ITE',9X,'RSDM',11X,'IM',3X,'JM')

      ITE=0
 20   CONTINUE

      CALL ADJ_ARTV

      TSTA=1./4.
      CALL ADJ_FLUX
      CALL ADJ_EXPLIC
      TSTA=1./2.
      CALL ADJ_FLUX
      CALL ADJ_EXPLIC
      TSTA=1.
      CALL ADJ_FLUX
      CALL ADJ_SMTH



      ITE=ITE+1
      KMD=MOD(ITE,NFSV)
      IF(KMD.NE.0) GOTO 30
C     WRITE(*,23) ITE,RSDM,IM,JM
 23   FORMAT(11X,I5,4X,F12.8,5X,3I5,2X,'(Field saved)')
      GO TO 45

 30   CONTINUE
      KK=MOD(ITE,NRSD)
      IF(KK.NE.0) GOTO 45
C     WRITE(*,40) ITE,RSDM,IM,JM
 40   FORMAT(11X,I5,4X,F12.8,5X,3I5)

 45   CONTINUE

      IF(RSDM.LT.ECRI) GOTO 50
      IF(ITE.GE.ITM) GOTO 70
      GOTO 20

 50   CONTINUE
C     WRITE(*,60) ITE,ECRI
 60   FORMAT(6X,'ITE =',I5,5X,'RSDM <'
     1     ,1X,'ECRI =',F15.9,5X,'STOP')
      GOTO 90

 70   CONTINUE
C     WRITE(*,80) ITM,RSDM
 80   FORMAT(6X,'ITE >',1X,'ITEmax ='
     1     ,I5,5X,'RSDM =',F15.9,5X,'STOP')

 90   CONTINUE

      ITL=ITE
      RETURN
      END
C------------------------------------------------------------------
      SUBROUTINE A1TA2T
      PARAMETER (ML=4,MI=161,MJ=91,MM=161)       
      COMMON /UVV/ U(ML,MI,MJ),V(ML,MI,MJ)
      COMMON /AT/ A1T(ML,ML,MI,MJ),A2T(ML,ML,MI,MJ)
      COMMON /RXI/ XXI(MI,MJ),YXI(MI,MJ)
      COMMON /RET/ XET(MI,MJ),YET(MI,MJ)

      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1

      COMMON /FAR/ GM,RMF,VXF,VYF,HF,AL
      COMMON /ADJFB/ TBT(ML,ML,MI)

      DIMENSION C2T(ML,ML),CT(ML,ML),BETA(ML),BTA(ML,ML)
      DIMENSION WR(ML),WI(ML),SCALE(ML),INT(ML), IS(ML),JS(ML)
      DIMENSION TR(ML,ML),TRI(ML,ML),PROD(ML,ML)

C     OPEN(UNIT=9,NAME='adjfar_eigen_invs.dat')

      DO J=1,NJ1
         DO I=1,NI1
            VX=V(2,I,J)/V(1,I,J)
            VY=V(3,I,J)/V(1,I,J)
            EN=V(4,I,J)/V(1,I,J)
            A1T(1,1,I,J)=0.0
            A1T(1,2,I,J)=0.5*((GM-3.)*VX*VX+(GM-1.)*VY*VY)
            A1T(1,3,I,J)=-VX*VY
            A1T(1,4,I,J)=-GM*EN*VX+(GM-1.)*VX*(VX*VX+VY*VY)
            A1T(2,1,I,J)=1.0
            A1T(2,2,I,J)=(3.-GM)*VX		  
            A1T(2,3,I,J)=VY		  
            A1T(2,4,I,J)=GM*EN-0.5*(GM-1.)*(3.*VX*VX+VY*VY)
            A1T(3,1,I,J)=0.0
            A1T(3,2,I,J)=-(GM-1.)*VY
            A1T(3,3,I,J)=VX
            A1T(3,4,I,J)=-(GM-1.)*VX*VY
            A1T(4,1,I,J)=0.0
            A1T(4,2,I,J)=GM-1.
            A1T(4,3,I,J)=0.
            A1T(4,4,I,J)=GM*VX

            A2T(1,1,I,J)=0.0
            A2T(1,2,I,J)=-VX*VY
            A2T(1,3,I,J)=0.5*((GM-1.)*VX*VX+(GM-3.)*VY*VY)
            A2T(1,4,I,J)=-GM*EN*VY+(GM-1.)*(VX*VX+VY*VY)*VY
            A2T(2,1,I,J)=0.0
            A2T(2,2,I,J)=VY		  
            A2T(2,3,I,J)=-(GM-1.)*VX		  
            A2T(2,4,I,J)=-(GM-1.)*VX*VY
            A2T(3,1,I,J)=1.0
            A2T(3,2,I,J)=VX
            A2T(3,3,I,J)=(3.-GM)*VY
            A2T(3,4,I,J)=GM*EN-0.5*(GM-1.)*(VX*VX+3.*VY*VY)
            A2T(4,1,I,J)=0.0
            A2T(4,2,I,J)=0.0
            A2T(4,3,I,J)=GM-1.
            A2T(4,4,I,J)=GM*VY
         END DO
	  END DO

      DO I=1,NI1
         DO L1=1,ML
            DO L2=1,ML
               C2T(L1,L2)=-YXI(I,NJ1)*A1T(L1,L2,I,NJ1)
     &              +XXI(I,NJ1)*A2T(L1,L2,I,NJ1)
               CT(L1,L2)=C2T(L1,L2)
            END DO
         END DO

C     -------- find eigenvalues and eigenvectors of C2T ------
C     -------- and transformation matrix TR too (each   ------
C     -------- column of TR is an eigenvector of C2T)   ------
         NM=ML
         N=ML
	     CALL BALANC(NM,N,CT,LOW,IGH,SCALE)
	     CALL ELMHES(NM,N,LOW,IGH,CT,INT)
         CALL ELTRAN(NM,N,LOW,IGH,CT,INT,TR)
         CALL HQR2(NM,N,LOW,IGH,CT,WR,WI,TR,IERR)
C     ------ if all eigenvalues are determined within 
C     ------ 30 iterations, IERR is set to zero;
C     ------ if more 30 iterations are required to 
C     ------ determine an eigenvalue, the subroutine
C     ------ terminates with IERR set to the index of
C     ------ the eigenvalue for which the failure occurs.
	     IF(IERR .EQ. 0) THEN
            M=N
            CALL BALBAK(NM,N,LOW,IGH,SCALE,M,TR)
         END IF

C     -------- find inverse of TR ------------
	     DO L2=1,ML
            DO L1=1,ML
               TRI(L1,L2)=TR(L1,L2)
            END DO
	     END DO
         CALL BRINV(TRI,N,INDC,IS,JS)
C     ------- if the returned INDC=0, it indicates that the---
C     ------- original matrix is ill-conditioned and there ---
C     ------- doesn't exist an inverse for the original matrix
C     ------------------------------

         DO L=1,ML
            IF(WR(L) .GT. 0.0) BETA(L)=0.0
            IF(WR(L) .LT. 0.0) BETA(L)=1.0
         END DO

         DO L1=1,ML
            DO L2=1,ML
               BTA(L1,L2)=0.0
            END DO
            BTA(L1,L1)=BETA(L1)
         END DO

C     ----- form coefficient matrix TBT for the I-th cell---
         DO L2=1,ML
            DO L1=1,ML
               PRD=0.0
               DO L3=1,ML
                  PRD=PRD+TR(L1,L3)*BTA(L3,L2);
               END DO
               PROD(L1,L2)=PRD
            END DO
         END DO

         DO L2=1,ML
            DO L1=1,ML
               PRD=0.0
               DO L3=1,ML
                  PRD=PRD+PROD(L1,L3)*TRI(L3,L2);
               END DO
               TBT(L1,L2,I)=PRD
            END DO
         END DO

      END DO
C     CLOSE(9)  

      RETURN 
	  END 

C-----------------------------------------------------------------
      SUBROUTINE ADJ_STEP
      PARAMETER (ML=4,MI=161,MJ=91,MM=161)               
      COMMON /UVV/ U(ML,MI,MJ),V(ML,MI,MJ)
      COMMON /TSP/ DTL(MI,MJ),SM(MI,MJ),PP(MI,MJ)
      COMMON /RAD/ RADI(MI,MJ),RADJ(MI,MJ)
      COMMON /VOL/ VOL(MI,MJ)
      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1

      COMMON /FAR/ GM,RMF,VXF,VYF,HF,AL
      COMMON /RSD/ RSDM,ECRI,TSTA,COUR4,EPS,ITL


      DO J=1,NJ1
         DO I=1,NI1
            RR=V(1,I,J)
            VX=V(2,I,J)/RR
            VY=V(3,I,J)/RR
            AA=SQRT(GM*PP(I,J)/RR)
            QQ=SQRT(VX**2+VY**2)                      
            DTL(I,J)=COUR4/AMAX1(QQ,AA)/SM(I,J)
         END DO
	  END DO

      RETURN
      END
C------------------------------------------------------------------
      SUBROUTINE ADJ_FLUX
      PARAMETER (ML=4,MI=161,MJ=91,MM=161)  
      COMMON /FQQ/ FDS(ML,MI,MJ),Q(ML,MI,MJ)
      COMMON /PSI/ PS0(ML,MI,MJ),PSI(ML,MI,MJ)
      COMMON /AT/ A1T(ML,ML,MI,MJ),A2T(ML,ML,MI,MJ)

      COMMON /GSFI/ SIX(MI,MJ),SIY(MI,MJ)
      COMMON /GSFJ/ SJX(MI,MJ),SJY(MI,MJ)
      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1

      
      CALL ADJ_BOUND(1)
      DO I=2,NI1
         I1=I-1
         DO J=1,NJ1
            DO L=1,ML
               FDS(L,I,J)=(PSI(L,I,J)+PSI(L,I1,J))/2.
            END DO
         END DO
      END DO

      DO I=1,NI1
         I1=I+1
         DO J=1,NJ1
            PS1SX=FDS(1,I1,J)*SIX(I1,J)-FDS(1,I,J)*SIX(I,J)
            PS2SX=FDS(2,I1,J)*SIX(I1,J)-FDS(2,I,J)*SIX(I,J)
            PS3SX=FDS(3,I1,J)*SIX(I1,J)-FDS(3,I,J)*SIX(I,J)
            PS4SX=FDS(4,I1,J)*SIX(I1,J)-FDS(4,I,J)*SIX(I,J)

            PS1SY=FDS(1,I1,J)*SIY(I1,J)-FDS(1,I,J)*SIY(I,J)
            PS2SY=FDS(2,I1,J)*SIY(I1,J)-FDS(2,I,J)*SIY(I,J)
            PS3SY=FDS(3,I1,J)*SIY(I1,J)-FDS(3,I,J)*SIY(I,J)
            PS4SY=FDS(4,I1,J)*SIY(I1,J)-FDS(4,I,J)*SIY(I,J)

            DO L1=1,ML
               Q(L1,I,J)=A1T(L1,1,I,J)*PS1SX
     1              +A1T(L1,2,I,J)*PS2SX
     2              +A1T(L1,3,I,J)*PS3SX
     3              +A1T(L1,4,I,J)*PS4SX
     4              +A2T(L1,1,I,J)*PS1SY
     5              +A2T(L1,2,I,J)*PS2SY
     6              +A2T(L1,3,I,J)*PS3SY
     7              +A2T(L1,4,I,J)*PS4SY
            END DO
         END DO
	  END DO

      CALL ADJ_BOUND(2)
      DO J=2,NJ1
         J1=J-1
         DO I=1,NI1
            DO L=1,ML
               FDS(L,I,J)=(PSI(L,I,J)+PSI(L,I,J1))/2. 
            END DO 
         END DO
      END DO

      DO J=1,NJ1
         J1=J+1
         DO I=1,NI1
            PS1SX=FDS(1,I,J1)*SJX(I,J1)-FDS(1,I,J)*SJX(I,J)
            PS2SX=FDS(2,I,J1)*SJX(I,J1)-FDS(2,I,J)*SJX(I,J)
            PS3SX=FDS(3,I,J1)*SJX(I,J1)-FDS(3,I,J)*SJX(I,J)
            PS4SX=FDS(4,I,J1)*SJX(I,J1)-FDS(4,I,J)*SJX(I,J)

            PS1SY=FDS(1,I,J1)*SJY(I,J1)-FDS(1,I,J)*SJY(I,J)
            PS2SY=FDS(2,I,J1)*SJY(I,J1)-FDS(2,I,J)*SJY(I,J)
            PS3SY=FDS(3,I,J1)*SJY(I,J1)-FDS(3,I,J)*SJY(I,J)
            PS4SY=FDS(4,I,J1)*SJY(I,J1)-FDS(4,I,J)*SJY(I,J)

            DO L1=1,ML
               Q(L1,I,J)=Q(L1,I,J)
     1              +A1T(L1,1,I,J)*PS1SX
     2              +A1T(L1,2,I,J)*PS2SX
     3              +A1T(L1,3,I,J)*PS3SX
     4              +A1T(L1,4,I,J)*PS4SX
     5              +A2T(L1,1,I,J)*PS1SY
     6              +A2T(L1,2,I,J)*PS2SY
     7              +A2T(L1,3,I,J)*PS3SY
     8              +A2T(L1,4,I,J)*PS4SY
            END DO
         END DO
	  END DO

      RETURN
      END
C------------------------------------------------------------------
      SUBROUTINE ADJ_BOUND(II)
      PARAMETER (ML=4,MI=161,MJ=91,MM=161)               
      COMMON /UVV/ U(ML,MI,MJ),V(ML,MI,MJ)
      COMMON /TSP/ DTL(MI,MJ),SM(MI,MJ),PP(MI,MJ)
      COMMON /FQQ/ FDS(ML,MI,MJ),Q(ML,MI,MJ)
      COMMON /PSI/ PS0(ML,MI,MJ),PSI(ML,MI,MJ)

      COMMON /GRD/ XX(MI,MJ),YY(MI,MJ)
      COMMON /BDOWD/DSTRBD(MI)

      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1 
      
      COMMON /FAR/ GM,RMF,VXF,VYF,HF,AL

      COMMON /PD/ PD(MI)
      
      COMMON /ADJFB/ TBT(ML,ML,MI)
	  
	  IF(II.EQ.2) GO TO 400 
      DO J=1,NJ1
         DO L=1,ML 
            FDS(L,1,J)=(PSI(L,1,J)+PSI(L,NI1,J))/2.
            FDS(L,NI,J)=FDS(L,1,J)
         END DO
      END DO
      RETURN


 400  CONTINUE
      DO I=1,NI1
         AA=1.+DSTRBD(I)
         BB=-DSTRBD(I)
         PPB=AA*PP(I,1)+BB*PP(I,2)
         XXX=XX(I+1,1)-XX(I,1)
         YYX=YY(I+1,1)-YY(I,1)
         DSSQ=XXX**2+YYX**2
         DS=SQRT(DSSQ)
         FDS(1,I,1)=PSI(1,I,1)
         FDS(2,I,1)=(-YYX/DSSQ)*(PPB-PD(I))
     1        +(XXX**2/DSSQ)*PSI(2,I,1)
     2        +(XXX*YYX/DSSQ)*PSI(3,I,1)
         FDS(3,I,1)=(XXX/DSSQ)*(PPB-PD(I))
     1        +(XXX*YYX/DSSQ)*PSI(2,I,1)
     2        +(YYX**2/DSSQ)*PSI(3,I,1)
         FDS(4,I,1)=PSI(4,I,1)
      END DO

C     ------- far-field boundary ----------------------------
      DO I=1,NI1        
         DO L=1,ML
            FDS(L,I,NJ)=0.0
            DO L2=1,ML
               FDS(L,I,NJ)=FDS(L,I,NJ)+TBT(L,L2,I)*PSI(L2,I,NJ1)
            END DO
         END DO
      END DO

      RETURN
      END
C------------------------------------------------------------------
      SUBROUTINE ADJ_EXPLIC
      PARAMETER (ML=4,MI=161,MJ=91,MM=161) 
      COMMON /TSP/ DTL(MI,MJ),SM(MI,MJ),PP(MI,MJ)
      COMMON /FQQ/ FDS(ML,MI,MJ),Q(ML,MI,MJ)
      COMMON /DDD/ D(ML,MI,MJ),AA(MM),BB(MM),CC(ML,MM)
      COMMON /PSI/ PS0(ML,MI,MJ),PSI(ML,MI,MJ)

      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1

      COMMON /RSD/ RSDM,ECRI,TSTA,COUR4,EPS,ITL  

      DO 20 J=1,NJ1
         DO 20 I=1,NI1
            DO 10 L=1,ML
               PSI(L,I,J)=PS0(L,I,J)+TSTA*DTL(I,J)*(Q(L,I,J)
     1              +D(L,I,J))
 10         CONTINUE
 20      CONTINUE

         RETURN
         END
C-----------------------------------------------------------------
      SUBROUTINE ADJ_SMTH
      PARAMETER (ML=4,MI=161,MJ=91,MM=161) 
      COMMON /TSP/ DTL(MI,MJ),SM(MI,MJ),PP(MI,MJ)
      COMMON /FQQ/ FDS(ML,MI,MJ),Q(ML,MI,MJ)
      COMMON /DDD/ D(ML,MI,MJ),AA(MM),BB(MM),CC(ML,MM)
      COMMON /PSI/ PS0(ML,MI,MJ),PSI(ML,MI,MJ)

      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1

      COMMON /RSD/ RSDM,ECRI,TSTA,COUR4,EPS,ITL  
      COMMON /MAX/ IM,JM


      RSDM=0.0
      DO 30 J=1,NJ1
         DO 20 I=1,NI1
            DO 10 L=1,ML
               Q(L,I,J)=TSTA*DTL(I,J)*(Q(L,I,J)+D(L,I,J))
               RSD=ABS(Q(L,I,J))
               IF(RSD.LT.RSDM) GOTO 10
               RSDM=RSD
               IM=I
               JM=J
 10         CONTINUE
 20      CONTINUE
 30   CONTINUE

      EPS1=EPS+1.0
C     ---------------- I direction ----------------------------------
      DO 110 J=1,NJ1
         DO 50 L=1,ML
            CC(L,1)=EPS1*Q(L,1,J)
            Q(L,NI1,J)=EPS1*Q(L,NI1,J)
 50      CONTINUE
         DO 70 I=2,NI1
            I1=I-1
            DO 60 L=1,ML
               CC(L,I)=Q(L,I,J)+AA(I1)*CC(L,I1)
 60         CONTINUE
 70      CONTINUE
         DO 80 L=1,ML
            Q(L,NI1,J)=CC(L,NI1)/BB(NI1)
 80      CONTINUE
         DO 100 I1=1,NI1-1
            I=NI1-I1
            I2=I+1
            DO 90 L=1,ML
               Q(L,I,J)=CC(L,I)/BB(I)+AA(I)*Q(L,I2,J)
 90         CONTINUE
 100     CONTINUE
 110  CONTINUE

C     ----------------- J direction --------------------------------
      DO 270 I=1,NI1
         DO 210 L=1,ML
            CC(L,1)=EPS1*Q(L,I,1)
            Q(L,I,NJ1)=EPS1*Q(L,I,NJ1)
 210     CONTINUE
         DO 230 J=2,NJ1
            J1=J-1
            DO 220 L=1,ML
               CC(L,J)=Q(L,I,J)+AA(J1)*CC(L,J1)
 220        CONTINUE
 230     CONTINUE
         DO 240 L=1,ML
            Q(L,I,NJ1)=CC(L,NJ1)/BB(NJ1)
 240     CONTINUE
         DO 260 J1=1,NJ1-1
            J=NJ1-J1
            J2=J+1
            DO 250 L=1,ML
               Q(L,I,J)=CC(L,J)/BB(J)+AA(J)*Q(L,I,J2)
 250        CONTINUE
 260     CONTINUE
 270  CONTINUE

C     -----------------------------------------------------------

      DO 810 J=1,NJ1
         DO 810 I=1,NI1
            DO 800 L=1,ML
               PSI(L,I,J)=PS0(L,I,J)+Q(L,I,J)
               PS0(L,I,J)=PSI(L,I,J)
 800        CONTINUE
 810     CONTINUE

         RETURN
         END

C-------------------------------------------------------------
      SUBROUTINE  ADJ_ARTV
      PARAMETER (ML=4,MI=161,MJ=91,MM=161)
      COMMON /TSP/ DTL(MI,MJ),SM(MI,MJ),PP(MI,MJ)
      COMMON /FQQ/ FDS(ML,MI,MJ),Q(ML,MI,MJ)
      COMMON /DDD/ D(ML,MI,MJ),AA(MM),BB(MM),CC(ML,MM)
      COMMON /DPP/ DP(MI,MJ),DPM(MI,MJ)
      COMMON /RAD/ RADI(MI,MJ),RADJ(MI,MJ)

      COMMON /PSI/ PS0(ML,MI,MJ),PSI(ML,MI,MJ)

      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1

      COMMON /E22/ E2I,E2J
      COMMON /E44/ E4I,E4J
      COMMON /FAR/ GM,RMF,VXF,VYF,HF,AL
      
      NI2=NI1-1
      NJ2=NJ1-1
      NI3=NI2-1
      NJ3=NJ2-1
C     omg=1.0 for transonic flow
C     omg=0.5 for supersonic or hypersonic
	  OMG=0.5

C     ---------- i direction ----------------------------
      DO J=1,NJ1
         DO I=2,NI2
            DP(I,J)=ABS(PP(I+1,J)-2.*PP(I,J)+PP(I-1,J))/
     1           ((1.-OMG)*(ABS(PP(I+1,J)-PP(I,J))
     2           +ABS(PP(I,J)-PP(I-1,J)))
     3           +OMG*(PP(I+1,J)+2.*PP(I,J)+PP(I-1,J)))
         END DO
	  END DO

      DO J=1,NJ1
         DO I=3,NI3
            DPM(I,J)=AMAX1(DP(I-1,J),DP(I,J),DP(I+1,J))
         END DO
         DPM(2,J)=DPM(3,J)
         DPM(1,J)=DPM(2,J)
         DPM(NI2,J)=DPM(NI3,J)
         DPM(NI1,J)=DPM(NI2,J)
	  END DO

      DO 30 L=1,ML
         DO 30 J=1,NJ1
            D(L,1,J)=0.0
            D(L,NI1,J)=0.0
 30      CONTINUE

         DO 80 I=1,NI2
            I1=I-1
            I2=I+1
            I3=I+2
            DO 60 J=1,NJ1
               RADN=(RADI(I,J)+RADI(I+1,J))/2.
               EP2=E2I*DPM(I,J)
               DO 40 L=1,ML
                  FDS(L,I,J)=RADN*EP2*(PSI(L,I2,J)-PSI(L,I,J))
 40            CONTINUE
               IF(I.EQ.1.OR.I.EQ.NI2) GOTO 60
               EP4=AMAX1(0.0,E4I-EP2)
               DO 50 L=1,ML
                  FDS(L,I,J)=FDS(L,I,J)-RADN*EP4*
     1                 (PSI(L,I3,J)-3.*PSI(L,I2,J)
     2                 +3.*PSI(L,I,J)-PSI(L,I1,J))
 50            CONTINUE
 60         CONTINUE
 80      CONTINUE

         DO 90 I=2,NI2
            I1=I-1
            DO 90 J=1,NJ1
               DO 90 L=1,ML
                  D(L,I,J)=FDS(L,I,J)-FDS(L,I1,J)
 90            CONTINUE

C     ---------- j direction ----------------------------
               DO I=1,NI1
                  DO J=2,NJ2
                     DP(I,J)=ABS(PP(I,J+1)-2.*PP(I,J)+PP(I,J-1))/
     1                    ((1.-OMG)*(ABS(PP(I,J+1)-PP(I,J))
     2                    +ABS(PP(I,J)-PP(I,J-1)))
     3                    +OMG*(PP(I,J+1)+2.*PP(I,J)+PP(I,J-1)))
                  END DO
               END DO

               DO I=1,NI1
                  DO J=3,NJ3
                     DPM(I,J)=AMAX1(DP(I,J-1),DP(I,J),DP(I,J+1))
                  END DO
                  DPM(I,2)=DPM(I,3)
                  DPM(I,1)=DPM(I,2)
                  DPM(I,NJ2)=DPM(I,NJ3)
                  DPM(I,NJ1)=DPM(I,NJ2)
               END DO

               DO 180 J=1,NJ2
                  J1=J-1
                  J2=J+1
                  J3=J+2
                  DO 160 I=1,NI1
                     RADN=(RADJ(I,J)+RADJ(I,J+1))/2.
                     EP2=E2J*DPM(I,J)
                     DO 140 L=1,ML
                        FDS(L,I,J)=RADN*EP2*(PSI(L,I,J2)-PSI(L,I,J))
 140                 CONTINUE
                     IF(J.EQ.1.OR.J.EQ.NJ2) GOTO 160
                     EP4=AMAX1(0.0,E4J-EP2)
                     DO 150 L=1,ML
                        FDS(L,I,J)=FDS(L,I,J)-RADN*EP4*
     1                       (PSI(L,I,J3)-3.*PSI(L,I,J2)
     2                       +3.*PSI(L,I,J)-PSI(L,I,J1))
 150                 CONTINUE
 160              CONTINUE
 180           CONTINUE

               DO 220 J=2,NJ2
                  J1=J-1
                  DO 210 I=1,NI1
                     DO 200 L=1,ML
                        D(L,I,J)=D(L,I,J)+FDS(L,I,J)-FDS(L,I,J1)
 200                 CONTINUE
 210              CONTINUE
 220           CONTINUE

               RETURN
               END
C--------end of adjoint equation solution subroutines -------------
C===================================================================


C-------------------------------------------------------------------
C     The following five subroutines (i.e.,BALANC,ELMHES,ELTRAN,HQR2,
C     BALBAK) are typed from the book written by B.T.Smith, J.M.Boyle,
C     J.J.Dongarra, B.S.Garbow, Y.Ikebe,V.C.Klema, and C.B.Moler:   
C     Matrix Eigensystem Routines--EISPACK Guide, Lecture Notes in 
C     Computer Science, Vol.6 (Second Edition),
C     edited by G.Goos and J.Hartmanis, Springger-Verlag, 1976
C     The routines was supposed to calculate all eigenvalues and 
C     corresponding eigenvectors of a real general matrix system.  
C     The typer and debugger: Zhang Zhengke, TSL, NUS, Nov.2002
C     -----------------------------------------------------------------
C     
C     First, an original real nonsymmetric matrix A is balanced first  
C     in subroutine "BALANC".On input, A is a matrix of order N to be
C     balanced. On output, A cntains the transformed matrix.NM(input):
C     row dimension	of A. N(input):must be not greater than NM. LOW,IGH
C     are output indicating the the boundary indices for the balanced 
C     matrix here and will be input in the subsequent subroutines. 
C     Scale:real output one-dimensional array of dimension at least N 
C     containing the information about the similarity transformation.
C     ******The subroutine BALBAK should be used to retrieve ***** 
C     ******the eigenvectors of the original matrix after the **** 
C     ******eigenvectors of the transformed matrix have been  ****
C     ******determined.									     *****
C 	  ************************************************************
C     
C     Second,the balanced matrix A is reduced to an upper Hessenberg
C     matrix by stablized elementary similarity transformations in
C     subroutine ELMHES. NM,N,LOW,IGH are all input here. A: on input,
C     A is a balanced matrix of order N to be reduced to upper 
C     Hessenberg form, on output, A contains the upper Hessenberg 
C     matrix as well as the multipliers used in the reduction. INT:
C     an integer one-dimensional array of dimension at least IGH 
C     identifying the rows and columns interchanged during the 
C     reduction.
C     
C     Third, call subroutine ELTRAN to accumulate the stabilized 
C     elementary similarity transformations used in the reduction 
C     of a real general matrix to upper Hessenberg form by ELMHES.
C     Here, NM,N,LOW,IGH,A,INT are all input. A has row dimension NM
C     and column dimension at least IGH.   Z is a real output 2-d
C     array with row dimension NM and column dimension at least N. 
C     It contains the transformation matrix produced in the reduction
C     by ELMHES to the upper Hessenberg form. 
C     
C     Fourth, calls HQR2 to compute the eigenvalues and eigenvectors
C     of a real upper Hessenberg matrix using the QR method.The 
C     eigenvectors of a real general matrix can also be computed if 
C     ELMHES and ELTRAN have been used to reduce this general matrix 
C     to Hessenberg form. NM,N,LOW,IGH,H are input. H: on input, it 
C     contains the upper Hessenberg matrix. In fact, H here can be 
C     taken as the matrix A out from ELTRAN. WR,WI: real output 1-d
C     array of dimension at least N containing the real and imaginary 
C     parts, respectively, of the eigenvalues of the Hessenberg matrix.
C     The eigenvalues are unordered except that complext conjugate 
C     pairs of eigenvalues appear consecutively with the eigenvalue 
C     having the positive imaginary part first. Z: if the eigenvectors
C     of a real general matrix are desired, Z contains the 
C     transformation matrix produced in ELTRAN which reduced the 
C     general matrrix to Hessenberg form, and on output, contains
C     the real and imaginary parts of the eigenvectors of this real
C     general matrix. If the j-th eigenvalue is real, the j-th column
C     of Z contains its eigenvector. If the j-th eigenvalue is complex
C     with positive imaginary part, the j-th and (j+1)-th columns of Z
C     contain the real and imaginary parts of its eigenvector.  
C     IERR: integer output. If all the eigenvalues are determined 
C     within 30 iterations, IERR is set zero. If more than 30 
C     iterations are required to determine an eigenvalue, this 
C     subroutine terminates with IERR set to the index of the 
C     eigenvalue for which the failure occurs. The eigenvalues in the 
C     WR and WI arrays should be correct for indices IERR+1, IERR+2,
C     ..., N, but no eigenvectors are computed.
C     
C     Finally, call BALBAK. BALBAK back transforms the eigenvectors of 
C     that real matrix transformed by BALANC. NM,N,LOW,IGH,SCALE are 
C     all input obtained in the above subroutines. M:integer input 
C     variable set to the number of columns of Z to be back transformed.
C     Z: 2-d array with row dimension NM and column dimension at least 
C     M. On input the first M columns of Z contain the real and 
C     imaginary parts of the eigenvectors to be transformed. On output,
C     these M columns of Z contain the real and imaginary parts of the
C     transformed eigenvectors. This Z should be the Z from HQR2.
C-------------------------------------------------------------------
      SUBROUTINE BALANC(NM,N,A,LOW,IGH,SCALE)
      INTEGER I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC
      REAL A(NM,N),SCALE(N)
      REAL C,F,G,R,S,B2,RADIX
      REAL ABS
      LOGICAL NOCONV
C     ----RADIX is a machine dependent parameter specifying the 
C     ----base of the machine floating point representation ---
      RADIX=2.

      B2=RADIX*RADIX
      K=1
      L=N
      GO TO 100
C     ------in-line procedure for row and column exchange-------
 20   SCALE(M)=J	 
      IF(J.EQ.M) GO TO 50
      
      DO I=1,L
         F=A(I,J)
         A(I,J)=A(I,M)
         A(I,M)=F
      END DO

      DO I=K,N
         F=A(J,I)
         A(J,I)=A(M,I)
         A(M,I)=F
      END DO

 50   GO TO (80,130), IEXC
C     ------search for rows isolating an eigenvalue and push them 
C     ------down ------------------------------------------------
 80   IF(L.EQ. 1) GO TO 280
      L=L-1
C     ------- for J=L step -1 until 1 DO ------------------------
 100  DO 120 JJ=1,L
         J=L+1-JJ

         DO 110 I=1,L
            IF(I.EQ.J) GO TO 110
            IF(A(J,I) .NE. 0.0) GO TO 120
 110     CONTINUE
         
         M=L
         IEXC=1
         GO TO 20
 120  CONTINUE
      GO TO 140
C     --------- search for columns isolating an eigenvalue 
C     --------- and push them left -----------------------
 130  K=K+1
 140  DO 170 J=K,L
         DO 150 I=K,L
            IF(I.EQ.J) GO TO 150
            IF(A(I,J) .NE. 0.0) GO TO 170
 150     CONTINUE  
         M=K
         IEXC=2
         GO TO 20
 170  CONTINUE
C     --------- now balance the submatrix in rows K to L -----------
      DO I=K,L
         SCALE(I)=1.0
      END DO
C     --------- iterative loop for norm reduction -----------
 190  NOCONV=.FALSE.

      DO 270 I=K,L
         C=0.0
         R=0.0
         DO 200 J=K,L
            IF(J.EQ.I) GO TO 200
            C=C+ABS(A(J,I))
            R=R+ABS(A(I,J))
 200     CONTINUE
C     ------- guard against zero C or R due to underflow ----------
         IF(C.EQ. 0.0 .OR. R.EQ. 0.0) GO TO 270
         G=R/RADIX
         F=1.0
         S=C+R
 210     IF(C .GE. G) GO TO 220
         F=F*RADIX
         C=C*B2
         GO TO 210
 220     G=R*RADIX
 230     IF(C .LT. G) GO TO 240
         F=F/RADIX
         C=C/B2
         GO TO 230
C     ------------- now balance ---------------
 240     IF((C+R)/F .GE. 0.95*S) GO TO 270
         G=1.0/F
         SCALE(I)=SCALE(I)*F
         NOCONV=.TRUE.

         DO J=K, N
            A(I,J)=A(I,J)*G
         END DO

         DO J=1,L
            A(J,I)=A(J,I)*F
         END DO
         
 270  CONTINUE

      IF(NOCONV) GO TO 190

 280  LOW=K
      IGH=L

      RETURN
      END
C--------------------------------------------------------------------
      SUBROUTINE ELMHES(NM,N,LOW,IGH,A,INT)
      INTEGER I,J,M,N,LA,NM,IGH,KP1,LOW,MM1,MP1
      REAL A(NM,N)
      REAL X,Y
      REAL ABS
      INTEGER INT(IGH)
      
      LA=IGH-1
      KP1=LOW+1
      IF(LA.LT.KP1) GO TO 200

      DO 180 M=KP1, LA
         MM1=M-1
         X=0.0
         I=M
         DO 100 J=M, IGH
            IF(ABS(A(J,MM1)) .LE. ABS(X)) GO TO 100
            X=A(J,MM1)
            I=J	   
 100     CONTINUE
         
         INT(M)=I
         IF(I.EQ.M) GO TO 130
C     ---------- interchange rows and columns of A ------------
         DO 110 J=MM1, N
            Y=A(I,J)
            A(I,J)=A(M,J)
            A(M,J)=Y	   
 110     CONTINUE

         DO 120 J=1, IGH
            Y=A(J,I)
            A(J,I)=A(J,M)
            A(J,M)=Y	   
 120     CONTINUE
C     --------- end interchange ------------
 130     IF(X .EQ. 0.0) GO TO 180
         MP1=M+1

         DO 160 I=MP1, IGH
            Y=A(I,MM1)
            IF(Y .EQ. 0.0) GO TO 160
            Y=Y/X      
            A(I,MM1)=Y

            DO J=M, N
               A(I,J)=A(I,J)-Y*A(M,J)
            END DO

            DO J=1, IGH
               A(J,M)=A(J,M)+Y*A(J,I)
            END DO

 160     CONTINUE

 180  CONTINUE

 200  RETURN
      END
C--------------------------------------------------------------------
      SUBROUTINE ELTRAN(NM,N,LOW,IGH,A,INT,Z)
      INTEGER I,J,N,KL,MM,MP,NM,IGH,LOW,MP1
C     REAL A(NM,IGH),Z(NM,N)
      REAL A(NM,N),Z(NM,N)
      INTEGER INT(IGH)

C     -------- initialize Z to identity matrix -------------
      DO 80 I=1, N
         
         DO J=1, N
            Z(I,J)=0.0
         END DO

         Z(I,I)=1.0
 80   CONTINUE
      
      KL=IGH-LOW-1
      IF(KL .LT. 1) GO TO 200
C     -------- for MP=IGH-1 step -1 until LOW+1 DO --------------
      DO 140 MM=1, KL
         MP=IGH-MM
         MP1=MP+1

         DO I=MP1, IGH
            Z(I,MP)=A(I,MP-1)  
         END DO

         I=INT(MP)
         IF(I .EQ. MP) GO TO 140

         DO J=MP, IGH
            Z(MP,J)=Z(I,J)
            Z(I,J)=0.0
         END DO

         Z(I,MP)=1.0
 140  CONTINUE

 200  RETURN
      END
C--------------------------------------------------------------------
      SUBROUTINE HQR2(NM,N,LOW,IGH,H,WR,WI,Z,IERR)
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,MM,NA,NM,NN,IGH,ITS,LOW,
     &     MP2,ENM2,IERR
      REAL H(NM,N),WR(N),WI(N),Z(NM,N)
      REAL P,Q,R,S,T,W,X,Y,RA,SA,VI,VR,ZZ,NORM,MACHEP
      REAL SQRT,ABS,SIGN
      INTEGER MIN0
      LOGICAL NOTLAS
      COMPLEX Z3
C---- COMPLEX CMPLX
      REAL REAL,AIMAG
      
C     ------- MACHEP is a machine dependent parameter specifying
C     ------- the relative precision of floating point arithmetic	 
      
      MACHEP=2.**(-26)	 
	  
      IERR=0
      NORM=0.0
      K=1
C     ---store roots isolated by balanc and compute matrix norm---
      DO 50 I=1,N

         DO J=K,N
            NORM=NORM+ABS(H(I,J))
         END DO

         K=I
         IF(I.GE.LOW .AND. I.LE.IGH) GO TO 50
         WR(I)=H(I,I)
         WI(I)=0.0
 50   CONTINUE
      
      EN=IGH
      T=0.0
C     ------- search for next eigenvalues -----------
 60   IF(EN. LT. LOW) GO TO 340
      ITS=0
      NA=EN-1
      ENM2=NA-1
C     ------- look for single small sub-diagonal element 
C     ------- for L=EN step -1 until 1 DO  -------------
 70   DO LL=LOW, EN
         L=EN+LOW-LL
         IF(L. EQ. LOW) GO TO 100
         S=ABS(H(L-1,L-1))+ABS(H(L,L)) 
         IF(S .EQ. 0.0) S=NORM
         IF(ABS(H(L,L-1)) .LE. MACHEP*S) GO TO 100
      END DO
C     -------- form shift ---------------
 100  X=H(EN,EN)
      IF(L.EQ.EN) GO TO 270
      Y=H(NA,NA)
      W=H(EN,NA)*H(NA,EN)
      IF(L .EQ. NA) GO TO 280
      IF(ITS .EQ. 30) GO TO 1000
      IF(ITS .NE. 10 .AND. ITS .NE. 20) GO TO 130
C     --------- form exceptional shift ----------
      T=T+X

      DO I=LOW, EN
         H(I,I)=H(I,I)-X
      END DO

      S=ABS(H(EN,NA))+ABS(H(NA,ENM2))
      X=0.75*S
      Y=X
      W=-0.4375*S*S
 130  ITS=ITS+1
C     ----looking for two consecutive small subdiagonal elements 
C     ----for M=EN-2 step -1 until L DO  -----------------------
      DO 140 MM=L, ENM2
         M=ENM2+L-MM
         ZZ=H(M,M)
         R=X-ZZ
         S=Y-ZZ
         P=(R*S-W)/H(M+1,M)+H(M,M+1)
         Q=H(M+1,M+1)-ZZ-R-S
         R=H(M+2,M+1)
         S=ABS(P)+ABS(Q)+ABS(R)
         P=P/S
         Q=Q/S
         R=R/S
         IF(M.EQ.L) GO TO 150
         IF(ABS(H(M,M-1))*(ABS(Q)+ABS(R)) .LE. MACHEP*ABS(P) 
     &        *(ABS(H(M-1,M-1))+ABS(ZZ)+ABS(H(M+1,M+1)))) GO TO 150
 140  CONTINUE

 150  MP2=M+2

      DO 160 I=MP2, EN
         H(I,I-2)=0.0
         IF(I.EQ.MP2) GO TO 160
         H(I,I-3)=0.0
 160  CONTINUE
C     ----------- double QR step involving rows L to EN 
C     ----------- and columns M to EN -----------------

      DO 260 K=M, NA
         NOTLAS=K.NE.NA
         IF(K.EQ.M) GO TO 170
         P=H(K,K-1) 
         Q=H(K+1,K-1)
         R=0.0
         IF(NOTLAS) R=H(K+2,K-1)
         X=ABS(P)+ABS(Q)+ABS(R)
         IF(X .EQ. 0.0) GO TO 260
         P=P/X		 
         Q=Q/X
         R=R/X
 170     S=SIGN(SQRT(P*P+Q*Q+R*R),P)
         IF(K.EQ.M) GO TO 180
         H(K,K-1)=-S*X
         GO TO 190
 180     IF(L.NE.M) H(K,K-1)=-H(K,K-1)
 190     P=P+S
         X=P/S
         Y=Q/S
         ZZ=R/S
         Q=Q/P
         R=R/P
C     -------- row modification ----------
         DO 210 J=K,N
            P=H(K,J)+Q*H(K+1,J)
            IF(.NOT. NOTLAS) GO TO 200
            P=P+R*H(K+2,J)
            H(K+2,J)=H(K+2,J)-P*ZZ
 200        H(K+1,J)=H(K+1,J)-P*Y
            H(K,J)=H(K,J)-P*X
 210     CONTINUE

         J=MIN0(EN,K+3)
C     --------- column modification ---------
         DO 230 I=1,J
            P=X*H(I,K)+Y*H(I,K+1)
            IF(.NOT. NOTLAS) GO TO 220
            P=P+ZZ*H(I,K+2)
            H(I,K+2)=H(I,K+2)-P*R
 220        H(I,K+1)=H(I,K+1)-P*Q
            H(I,K)=H(I,K)-P
 230     CONTINUE
C     --------- accumulate transformations ---------------
         DO 250 I=LOW, IGH
            P=X*Z(I,K)+Y*Z(I,K+1)
            IF(.NOT. NOTLAS) GO TO 240
            P=P+ZZ*Z(I,K+2)
            Z(I,K+2)=Z(I,K+2)-P*R
 240        Z(I,K+1)=Z(I,K+1)-P*Q
            Z(I,K)=Z(I,K)-P
 250     CONTINUE

 260  CONTINUE

      GO TO 70
C     ----------- one root found --------------------------
 270  H(EN,EN)=X+T
      WR(EN)=H(EN,EN)
      WI(EN)=0.0
      EN=NA
      GO TO 60
C     ----------- two roots found -------------------------
 280  P=(Y-X)/2.0
      Q=P*P+W
      ZZ=SQRT(ABS(Q))
      H(EN,EN)=X+T
      X=H(EN,EN)
      H(NA,NA)=Y+T
      IF(Q .LT. 0.0) GO TO 320
C     --------- real pair ---------------
      ZZ=P+SIGN(ZZ,P)
      WR(NA)=X+ZZ
      WR(EN)=WR(NA)
      IF(ZZ .NE. 0.0) WR(EN)=X-W/ZZ
      WI(NA)=0.0
      WI(EN)=0.0
      X=H(EN,NA)
      S=ABS(X)+ABS(ZZ)
      P=X/S
      Q=ZZ/S
      R=SQRT(P*P+Q*Q)
      P=P/R
      Q=Q/R
C     ----------- row modification -----------------
      DO 290 J=NA, N
         ZZ=H(NA,J)
         H(NA,J)=Q*ZZ+P*H(EN,J)
         H(EN,J)=Q*H(EN,J)-P*ZZ
 290  CONTINUE
C     ----------- column modification ----------------
      DO 300 I=1, EN
         ZZ=H(I,NA)
         H(I,NA)=Q*ZZ+P*H(I,EN)
         H(I,EN)=Q*H(I,EN)-P*ZZ
 300  CONTINUE
C     ----------- accumulate transformations ----------------
      DO 310 I=LOW, IGH
         ZZ=Z(I,NA)
         Z(I,NA)=Q*ZZ+P*Z(I,EN)
         Z(I,EN)=Q*Z(I,EN)-P*ZZ
 310  CONTINUE

      GO TO 330

C     ----------- complex pair -------------------
 320  WR(NA)=X+P
      WR(EN)=X+P
      WI(NA)=ZZ
      WI(EN)=-ZZ
 330  EN=ENM2
      GO TO 60
C     ---------- all roots found.-------------------------------------
C     ---- backsubstitute to find vecrors of upper traiangular form---
 340  IF(NORM .EQ. 0.0) GO TO 1001
C     --------- for EN=N step -1 until 1 DO  ---------------------
      DO 800 NN=1, N
         EN=N+1-NN 
         P=WR(EN)
         Q=WI(EN)
         NA=EN-1
         IF(Q) 710, 600, 800
C     ----------- real vector --------------
 600     M=EN
         H(EN,EN)=1.0
         IF(NA .EQ. 0) GO TO 800
C     ----------- for I=EN-1 step -1 until 1 DO ----------------------
         DO 700 II=1, NA
            I=EN-II
            W=H(I,I)-P
            R=H(I,EN)
            IF(M.GT.NA) GO TO 620
C     
            DO J=M, NA
               R=R+H(I,J)*H(J,EN)
            END DO
            
 620        IF(WI(I) .GE. 0.0) GO TO 630
            ZZ=W
            S=R
            GO TO 700
 630        M=I
            IF(WI(I) .NE. 0.0) GO TO 640       
            T=W
            IF(W .EQ. 0.0) T=MACHEP*NORM
            H(I,EN)=-R/T
            GO TO 700
C     ---------- solve real equations ----------
 640        X=H(I,I+1)
            Y=H(I+1,I)
            Q=(WR(I)-P)*(WR(I)-P)+WI(I)*WI(I)
            T=(X*S-ZZ*R)/Q
            H(I,EN)=T
            IF(ABS(X) .LE. ABS(ZZ)) GO TO 650
            H(I+1,EN)=(-R-W*T)/X
            GO TO 700
 650        H(I+1,EN)=(-S-Y*T)/ZZ
 700     CONTINUE
C     ---------- end real vector ---------------
         GO TO 800
C     ---------- complex vector ----------------
 710     M=NA
C     ---------- last vector component chosen imaginary so that
C     ---------- eigenvector matrix is triangular -------------
         IF(ABS(H(EN,NA)) .LE. ABS(H(NA,EN))) GO TO 720
         H(NA,NA)=Q/H(EN,NA)
         H(NA,EN)=-(H(EN,EN)-P)/H(EN,NA)
         GO TO 730
 720     Z3=CMPLX(0.0,-H(NA,EN))/CMPLX(H(NA,NA)-P,Q)
         H(NA,NA)=REAL(Z3)
         H(NA,EN)=AIMAG(Z3)
 730     H(EN,NA)=0.0
         H(EN,EN)=1.0
         ENM2=NA-1
         IF(ENM2 .EQ. 0) GO TO 800
C     --------- for I=EN-2 step -1 until 1 DO -----------------
         DO 790 II=1, ENM2
            I=NA-II
            W=H(I,I)-P
            RA=0.0
            SA=H(I,EN)

            DO J=M, NA
               RA=RA+H(I,J)*H(J,NA)
               SA=SA+H(I,J)*H(J,EN)
            END DO

            IF(WI(I) .GE. 0.0) GO TO 770
            ZZ=W
            R=RA
            S=SA
            GO TO 790
 770        M=I
            IF(WI(I) .NE. 0.0) GO TO 780
            Z3=CMPLX(-RA,-SA)/CMPLX(W,Q)
            H(I,NA)=REAL(Z3)
            H(I,EN)=AIMAG(Z3)
            GO TO 790
C     ---------- solve complex equations -----------------
 780        X=H(I,I+1)             
            Y=H(I+1,I)
            VR=(WR(I)-P)*(WR(I)-P)+WI(I)*WI(I)-Q*Q
            VI=(WR(I)-P)*2.0*Q
            IF(VR .EQ. 0.0 .AND. VI .EQ. 0.0) VR=MACHEP*NORM
     &           *(ABS(W)+ABS(Q)+ABS(X)+ABS(Y)+ABS(ZZ))
            Z3=CMPLX(X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA)/CMPLX(VR,VI)
            H(I,NA)=REAL(Z3)
            H(I,EN)=AIMAG(Z3)
            IF(ABS(X) .LE. ABS(ZZ)+ABS(Q)) GO TO 785 
            H(I+1,NA)=(-RA-W*H(I,NA)+Q*H(I,EN))/X
            H(I+1,EN)=(-SA-W*H(I,EN)-Q*H(I,NA))/X
            GO TO 790
 785        Z3=CMPLX(-R-Y*H(I,NA),-S-Y*H(I,EN))/CMPLX(ZZ,Q)
            H(I+1,NA)=REAL(Z3)
            H(I+1,EN)=AIMAG(Z3)
 790     CONTINUE
C     ------------ end complex vector -------------
 800  CONTINUE
C     ------------ end back substitution -------------------------
C     ------------ vectors of isolated roots ---------------------
      DO 840 I=1, N
         IF(I.GE.LOW .AND. I.LE.IGH) GO TO 840
         
         DO J=I, N
            Z(I,J)=H(I,J)
         END DO

 840  CONTINUE

C     ------ multiply by transformation matrix to give ------------
C     ------ vectors of original full matrix           ------------
C     ------ for J=N step -1 until LOW DO  ------------------------
      DO 880 JJ=LOW, N
         J=N+LOW-JJ
         M=MIN0(J,IGH)

         DO 880 I=LOW, IGH
            ZZ=0.0
            
            DO K=LOW, M
               ZZ=ZZ+Z(I,K)*H(K,J)
            END DO

            Z(I,J)=ZZ
 880     CONTINUE

         GO TO 1001
C     ---- set error --no convergence to an eigenvalue ------------
C     ---- after 30 iterations                         ------------
 1000    IERR=EN
 1001    RETURN
         END
C--------------------------------------------------------------------
      SUBROUTINE BALBAK(NM,N,LOW,IGH,SCALE,M,Z)
      INTEGER I,J,K,M,N,II,NM,IGH,LOW
      REAL SCALE(N),Z(NM,N)
      REAL S

      IF(M .EQ. 0) GO TO 200
      IF(IGH .EQ. LOW) GO TO 120

      DO 110 I=LOW, IGH
         S=SCALE(I)
C     --- left hand eigenvectors are back transformed if the ------ 
C     --- foregoing statement is replaced by	S=1.0/SCALE(I)	------
         DO J=1, M
            Z(I,J)=Z(I,J)*S
         END DO

 110  CONTINUE
C     ------ for I=LOW-1 step -1 until 1, -------------
C     ------     IGH+1 step 1 until N DO  -------------
 120  DO 140 II=1, N
         I=II
         IF(I.GE.LOW .AND. I.LE.IGH) GO TO 140
         IF(I.LT.LOW) I=LOW-II
         K=SCALE(I)
         IF(K.EQ.I) GO TO 140

         DO 130 J=1, M
            S=Z(I,J)
            Z(I,J)=Z(K,J)
            Z(K,J)=S
 130     CONTINUE

 140  CONTINUE
      
 200  RETURN
      END
C---- end of the five subroutines for finding eigenvalues and-------
C---- corresponding eigenvectors of a general real matrix ----------
C------------------------------------------------------------------
C     
C---------------------------------------------------------------------
C     The following subroutines BRINV and BCINV are used to find the 
C     inverse of a real matrix and a complext matrix, respectively 
C     by Gauss-Jordan elimination method.
C     They are copied directly from the disk acompanying a book titled
C     "A Collection of Commonly Used FORTRAN Algorithms" (in Chinese)
C     written by Xu Shiliang and published by Tsinghua University Press.
C     The book and the disk are supplied by Dr. Wang Zhengyi of TSL.
C---------------------------------------------------------------------
      SUBROUTINE BRINV(A,N,L,IS,JS)
      DIMENSION A(N,N),IS(N),JS(N)
C---- DOUBLE PRECISION A,T,D
C--------------------------------------------------------------------
C     A: real N*N array,the original matrix to be inversed on input,
C     will be the inverse matrix on output returned.
C     N: integer input vwriable, order of A.
C     L: integer output vwriable, L=0 means that original matrix A is 
C     ill-conditioned and failure to find the inverse. 
C     IS,JS:integer 1-d array with dimension N, working array.
C--------------------------------------------------------------------
      L=1
      DO 100 K=1,N
         D=0.0
         DO 10 I=K,N
            DO 10 J=K,N
               IF (ABS(A(I,J)).GT.D) THEN
                  D=ABS(A(I,J))
                  IS(K)=I
                  JS(K)=J
               END IF
 10         CONTINUE
            IF (D+1.0.EQ.1.0) THEN
               L=0
               RETURN
            END IF
 20         FORMAT(1X,'ERR**NOT INV')
            DO 30 J=1,N
               T=A(K,J)
               A(K,J)=A(IS(K),J)
               A(IS(K),J)=T
 30         CONTINUE
            DO 40 I=1,N
               T=A(I,K)
               A(I,K)=A(I,JS(K))
               A(I,JS(K))=T
 40         CONTINUE
            A(K,K)=1/A(K,K)
            DO 50 J=1,N
               IF (J.NE.K) THEN
                  A(K,J)=A(K,J)*A(K,K)
               END IF
 50         CONTINUE
            DO 70 I=1,N
               IF (I.NE.K) THEN
                  DO 60 J=1,N
                     IF (J.NE.K) THEN
                        A(I,J)=A(I,J)-A(I,K)*A(K,J)
                     END IF
 60               CONTINUE
               END IF
 70         CONTINUE
            DO 80 I=1,N
               IF (I.NE.K) THEN
                  A(I,K)=-A(I,K)*A(K,K)
               END IF
 80         CONTINUE
 100     CONTINUE
         DO 130 K=N,1,-1
            DO 110 J=1,N
               T=A(K,J)
               A(K,J)=A(JS(K),J)
               A(JS(K),J)=T
 110        CONTINUE
            DO 120 I=1,N
               T=A(I,K)
               A(I,K)=A(I,IS(K))
               A(I,IS(K))=T
 120        CONTINUE
 130     CONTINUE
         RETURN
         END
C--------------------------------------------------------------------
      SUBROUTINE BCINV(AR,AI,N,L,IS,JS)
      DIMENSION AR(N,N),AI(N,N),IS(N),JS(N)
C---- DOUBLE PRECISION AR,AI,D,P,T,Q,S,B
C--------------------------------------------------------------------
C     AR,AI: real N*N array, real and imaginary parts of the the original 
C     complex matrix to be inversed on input,
C     will be the real and imaginary parts of inverse matrix 
C     on output returned.
C     N: integer input vwriable, order of the matrices.
C     L: integer output vwriable, L=0 means that original matrix  is 
C     ill-conditioned and failure to find the inverse. 
C     IS,JS:integer 1-d array with dimension N, working array.
C--------------------------------------------------------------------
      L=1
      DO 100 K=1,N
         D=0.0
         DO 10 I=K,N
            DO 10 J=K,N
               P=AR(I,J)*AR(I,J)+AI(I,J)*AI(I,J)
               IF (P.GT.D) THEN
                  D=P
                  IS(K)=I
                  JS(K)=J
               END IF
 10         CONTINUE
            IF (D+1.0.EQ.1.0) THEN
               L=0
               RETURN
            END IF
 20         FORMAT(1X,'ERR**NOT INV')
            DO 30 J=1,N
               T=AR(K,J)
               AR(K,J)=AR(IS(K),J)
               AR(IS(K),J)=T
               T=AI(K,J)
               AI(K,J)=AI(IS(K),J)
               AI(IS(K),J)=T
 30         CONTINUE
            DO 40 I=1,N
               T=AR(I,K)
               AR(I,K)=AR(I,JS(K))
               AR(I,JS(K))=T
               T=AI(I,K)
               AI(I,K)=AI(I,JS(K))
               AI(I,JS(K))=T
 40         CONTINUE
            AR(K,K)=AR(K,K)/D
            AI(K,K)=-AI(K,K)/D
            DO 50 J=1,N
               IF (J.NE.K) THEN
                  P=AR(K,J)*AR(K,K)
                  Q=AI(K,J)*AI(K,K)
                  S=(AR(K,J)+AI(K,J))*(AR(K,K)+AI(K,K))
                  AR(K,J)=P-Q
                  AI(K,J)=S-P-Q
               END IF
 50         CONTINUE
            DO 70 I=1,N
               IF (I.NE.K) THEN
                  DO 60 J=1,N
                     IF (J.NE.K) THEN
                        P=AR(K,J)*AR(I,K)
                        Q=AI(K,J)*AI(I,K)
                        S=(AR(K,J)+AI(K,J))*(AR(I,K)+AI(I,K))
                        T=P-Q
                        B=S-P-Q
                        AR(I,J)=AR(I,J)-T
                        AI(I,J)=AI(I,J)-B
                     END IF
 60               CONTINUE
               END IF
 70         CONTINUE
            DO 80 I=1,N
               IF (I.NE.K) THEN
                  P=AR(I,K)*AR(K,K)
                  Q=AI(I,K)*AI(K,K)
                  S=(AR(I,K)+AI(I,K))*(AR(K,K)+AI(K,K))
                  AR(I,K)=Q-P
                  AI(I,K)=P+Q-S
               END IF
 80         CONTINUE
 100     CONTINUE
         DO 130 K=N,1,-1
            DO 110 J=1,N
               T=AR(K,J)
               AR(K,J)=AR(JS(K),J)
               AR(JS(K),J)=T
               T=AI(K,J)
               AI(K,J)=AI(JS(K),J)
               AI(JS(K),J)=T
 110        CONTINUE
            DO 120 I=1,N
               T=AR(I,K)
               AR(I,K)=AR(I,IS(K))
               AR(I,IS(K))=T
               T=AI(I,K)
               AI(I,K)=AI(I,IS(K))
               AI(I,IS(K))=T
 120        CONTINUE
 130     CONTINUE
         RETURN
         END
C-------end of subroutines for finding inverse of a matrix --------
C---- also end of all subroutines for solving adjoint equation -----
C===================================================================
C===================================================================


C===================================================================
C===================================================================
C-------------------------------------------------------------------
C-------------------------------------------------------------------
C-------------------------------------------------------------------
      SUBROUTINE GRAD
      PARAMETER (ML=4,MI=161,MJ=91,MM=161)
      PARAMETER (ND=24)
C     ------------------------------------------------------
      COMMON /UVV/ U(ML,MI,MJ),V(ML,MI,MJ)
      COMMON /TSP/ DTL(MI,MJ),SM(MI,MJ),PP(MI,MJ)
      COMMON /GRD/ XX(MI,MJ),YY(MI,MJ)
      COMMON /RXI/ XXI(MI,MJ),YXI(MI,MJ)
      COMMON /RET/ XET(MI,MJ),YET(MI,MJ)
      COMMON /BDOWD/DSTRBD(MI)
      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1
C     ---------------------------------------- 
      COMMON /PSI/ PS0(ML,MI,MJ),PSI(ML,MI,MJ)
      COMMON /PD/ PD(MI)
      COMMON /ADJFB/ TBT(ML,ML,MI)
C     ----------------------------------------	       	  
      COMMON /NEWG/ XNEW(MI,MJ),YNEW(MI,MJ)
      COMMON /DSV/ALPHA(ND),ALPHA0(ND)
      COMMON /DLI/DLI
	  COMMON /G/G(ND)
C     ----------------------------------------------------------

      DLTA=0.000001
	  DO 100 K=1,ND

         DO I=1,ND
            IF(I.EQ.K) THEN
	           ALPHA(I)=ALPHA0(I)+DLTA
            ELSE 
	           ALPHA(I)=ALPHA0(I)
            END IF
         END DO

         CALL NEWGRD

         CALL DLTI

         G(K)=DLI/DLTA

 100  CONTINUE
      ALPHA(ND)=ALPHA0(ND)

c$$$	  GMAX=0.0
c$$$      DO I=1,ND
c$$$         IF(ABS(G(I)) .GT. GMAX) GMAX=ABS(G(I))
c$$$      END DO
c$$$
c$$$      DO I=1,ND
c$$$         G(I)=G(I)/GMAX
c$$$      END DO
      
	  RETURN
      END
C------------------------------------------------------------------
      SUBROUTINE NEWGRD
      PARAMETER (ML=4,MI=161,MJ=91,MM=161)
      PARAMETER (ND=24)
C     -----------------------------------------------------------
      COMMON /GRD/ XX(MI,MJ),YY(MI,MJ)
      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1
      
      COMMON /NEWG/ XNEW(MI,MJ),YNEW(MI,MJ)
C     -----------------------------------------------------------
      COMMON /DSV/ALPHA(ND),ALPHA0(ND)
      DIMENSION ARC(MM)
C     -----------------------------------------------------------

      CHORD=ABS(XX(IHF,1)-XX(1,1))
      XLE=0.0

C     ------ upper half body surface ----------------
      NDHF=ND/2
      DO I=1,IHF
         XI=FLOAT(I-1)/FLOAT(IHF-1)
         XNEW(IHF+I-1,1)=XX(IHF+I-1,1)
         XC=(XNEW(IHF+I-1,1)-XLE)/CHORD
         YNEW(IHF+I-1,1)=YUS(ALPHA,XC)
      END DO
C     ------ lower half body surface ----------------
      DO I=1,IHF
         XI=FLOAT(I-1)/FLOAT(IHF-1)
         XNEW(I,1)=XX(I,1)
         XC=(XNEW(I,1)-XLE)/CHORD
         YNEW(I,1)=YLS(ALPHA(NDHF+1),XC)
      END DO
C     -----------------------------------------------
C     -----------------------------------------------
      DO I=1,NI

         ARC(NJ)=0.0
         DO JINV=1,NJ
            J=NJ-JINV+1
            ARC(J)=ARC(J+1)+SQRT((XX(I,J+1)-XX(I,J))**2
     &           +(YY(I,J+1)-YY(I,J))**2)
         END DO
         ARCTT=ARC(1)

         DO J=1,NJ
            ARC(J)=ARC(J)/ARCTT
         END DO
         
         DO J=1,NJ
            XNEW(I,J)=XX(I,J)+ARC(J)*(XNEW(I,1)-XX(I,1))
            YNEW(I,J)=YY(I,J)+ARC(J)*(YNEW(I,1)-YY(I,1))
         END DO

      END DO

      RETURN
      END

C-------------------------------------------------------------------
      SUBROUTINE DLTI
      PARAMETER (ML=4,MI=161,MJ=91,MM=161)
C     ------------------------------------------------------
      COMMON /UVV/ U(ML,MI,MJ),V(ML,MI,MJ)
      COMMON /TSP/ DTL(MI,MJ),SM(MI,MJ),PP(MI,MJ)
      COMMON /GRD/ XX(MI,MJ),YY(MI,MJ)
      COMMON /RXI/ XXI(MI,MJ),YXI(MI,MJ)
      COMMON /RET/ XET(MI,MJ),YET(MI,MJ)
      COMMON /BDOWD/DSTRBD(MI)
      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1
C     ---------------------------------------- 
      COMMON /PSI/ PS0(ML,MI,MJ),PSI(ML,MI,MJ)
      COMMON /PD/ PD(MI)
      COMMON /ADJFB/ TBT(ML,ML,MI)
C     ----------------------------------------	       	  
      COMMON /NEWG/ XNEW(MI,MJ),YNEW(MI,MJ)
      COMMON /DLI/DLI
      DIMENSION DLTS(ML,ML,MI,MJ),DS21SF(MI),DS22SF(MI)
      DIMENSION DS1F(ML,MI,MJ),DS2F(ML,MI,MJ)
      DIMENSION PSBD(ML,MI),PSFB(ML,MI)
      DIMENSION PSXI(ML,MI,MJ),PSET(ML,MI,MJ)  
C     ----------------------------------------------------------
C     ------------- dSij in interior field -----------------
      DO J=1,NJ1
         J1=J+1
         DO I=1,NI1
            I1=I+1
            XXINEW=
     &           (XNEW(I1,J)+XNEW(I1,J1))/2.-(XNEW(I,J)+XNEW(I,J1))/2.
            YXINEW=
     &           (YNEW(I1,J)+YNEW(I1,J1))/2.-(YNEW(I,J)+YNEW(I,J1))/2.
            XETNEW=
     &           (XNEW(I,J1)+XNEW(I1,J1))/2.-(XNEW(I,J)+XNEW(I1,J))/2.
            YETNEW=
     &           (YNEW(I,J1)+YNEW(I1,J1))/2.-(YNEW(I,J)+YNEW(I1,J))/2.
            S11NEW=YETNEW
            S12NEW=-XETNEW
            S21NEW=-YXINEW
            S22NEW=XXINEW

            S11=YET(I,J)
            S12=-XET(I,J)
            S21=-YXI(I,J)
            S22=XXI(I,J)
            
            DLTS(1,1,I,J)=S11NEW-S11       
            DLTS(1,2,I,J)=S12NEW-S12
            DLTS(2,1,I,J)=S21NEW-S21       
            DLTS(2,2,I,J)=S22NEW-S22
         END DO
      END DO

C     ----------- dS21, dS22 on airfoil surface -------------
      DO I=1,NI1
         XXIO=XX(I+1,1)-XX(I,1)
         YXIO=YY(I+1,1)-YY(I,1)
         XXIN=XNEW(I+1,1)-XNEW(I,1)
         YXIN=YNEW(I+1,1)-YNEW(I,1)
         DS21SF(I)=(-YXIN)-(-YXIO)
         DS22SF(I)=(XXIN)-(XXIO)
      END DO

C     -------- dS11*f1+dS12*f2 in interior field -------------
      DO J=1,NJ1
         DO I=1,NI1
            DS1F(1,I,J)=DLTS(1,1,I,J)*V(2,I,J)
     &           +DLTS(1,2,I,J)*V(3,I,J)
            DS1F(2,I,J)=DLTS(1,1,I,J)
     1           *(V(2,I,J)*V(2,I,J)/V(1,I,J)+PP(I,J))
     2           +DLTS(1,2,I,J)*V(3,I,J)*V(2,I,J)/V(1,I,J)
            DS1F(3,I,J)=DLTS(1,1,I,J)*V(2,I,J)*V(3,I,J)/V(1,I,J)
     1           +DLTS(1,2,I,J)
     2           *(V(3,I,J)*V(3,I,J)/V(1,I,J)+PP(I,J))
            DS1F(4,I,J)=DLTS(1,1,I,J)*(V(2,I,J)/V(1,I,J))
     1           *(V(4,I,J)+PP(I,J))
     2           +DLTS(1,2,I,J)*(V(3,I,J)/V(1,I,J))
     3           *(V(4,I,J)+PP(I,J))
         END DO
      END DO

C     --------- dS21*f1+dS22*f2 in interior field --------------
      DO J=1,NJ1
         DO I=1,NI1
            DS2F(1,I,J)=DLTS(2,1,I,J)*V(2,I,J)
     &           +DLTS(2,2,I,J)*V(3,I,J)
            DS2F(2,I,J)=DLTS(2,1,I,J)
     1           *(V(2,I,J)*V(2,I,J)/V(1,I,J)+PP(I,J))
     2           +DLTS(2,2,I,J)*V(3,I,J)*V(2,I,J)/V(1,I,J)
            DS2F(3,I,J)=DLTS(2,1,I,J)*V(2,I,J)*V(3,I,J)/V(1,I,J)
     1           +DLTS(2,2,I,J)
     2           *(V(3,I,J)*V(3,I,J)/V(1,I,J)+PP(I,J))
            DS2F(4,I,J)=DLTS(2,1,I,J)*(V(2,I,J)/V(1,I,J))
     1           *(V(4,I,J)+PP(I,J))
     2           +DLTS(2,2,I,J)*(V(3,I,J)/V(1,I,J))
     3           *(V(4,I,J)+PP(I,J))
         END DO
      END DO

C     ------- PSI on airfoil surface -------------
      DO I=1,NI1
         AA=1.+DSTRBD(I)
         BB=-DSTRBD(I)
         PPB=AA*PP(I,1)+BB*PP(I,2)
         XXX=XX(I+1,1)-XX(I,1)
         YYX=YY(I+1,1)-YY(I,1)
         DSSQ=XXX**2+YYX**2
         PSBD(1,I)=PSI(1,I,1)
         PSBD(2,I)=(-YYX/DSSQ)*(PPB-PD(I))
     1        +(XXX**2/DSSQ)*PSI(2,I,1)
     2        +(XXX*YYX/DSSQ)*PSI(3,I,1)
         PSBD(3,I)=(XXX/DSSQ)*(PPB-PD(I))
     1        +(XXX*YYX/DSSQ)*PSI(2,I,1)
     2        +(YYX**2/DSSQ)*PSI(3,I,1)
         PSBD(4,I)=PSI(4,I,1)
      END DO

C     ------- PSI on far-field boundary -------------
      DO I=1,NI1        
         DO L=1,ML
            PSFB(L,I)=0.0
            DO L2=1,ML
               PSFB(L,I)=PSFB(L,I)+TBT(L,L2,I)*PSI(L2,I,NJ1)
            END DO
         END DO
      END DO

C     ------- dPSI/dXI,dPSI/dET in interior field ---------
      DO L=1,ML
         DO J=1,NJ1
            DO I=2,NI1-1
               PSXI(L,I,J)=(PSI(L,I+1,J)-PSI(L,I-1,J))/2.0
            END DO
            PSXI(L,1,J)=(PSI(L,2,J)-PSI(L,NI1,J))/2.0
            PSXI(L,NI1,J)=(PSI(L,1,J)-PSI(L,NI1-1,J))/2.0
         END DO
      END DO

      DO L=1,ML
         DO I=1,NI1
            DO J=2,NJ1-1
               PSET(L,I,J)=(PSI(L,I,J+1)-PSI(L,I,J-1))/2.0
            END DO
            PSET(L,I,1)=(PSI(L,I,1)+PSI(L,I,2))/2.0-PSBD(L,I)
            PSET(L,I,NJ1)=PSFB(L,I)-(PSI(L,I,NJ1-1)+PSI(L,I,NJ1))/2.0
         END DO
      END DO

C     ----------------------------------------------------
      DLI=0.0
      DO J=1,NJ1
         DO I=1,NI1
            DLI=DLI+PSXI(1,I,J)*DS1F(1,I,J)
     1           +PSXI(2,I,J)*DS1F(2,I,J)
     2           +PSXI(3,I,J)*DS1F(3,I,J)
     3           +PSXI(4,I,J)*DS1F(4,I,J)
     4           +PSET(1,I,J)*DS2F(1,I,J)
     5           +PSET(2,I,J)*DS2F(2,I,J)
     6           +PSET(3,I,J)*DS2F(3,I,J)
     7           +PSET(4,I,J)*DS2F(4,I,J)
         END DO
      END DO

      DO I=1,NI1
         AA=1.+DSTRBD(I)
         BB=-DSTRBD(I)
         PPB=AA*PP(I,1)+BB*PP(I,2)
         DLI=DLI+(PSBD(2,I)*DS21SF(I)+PSBD(3,I)*DS22SF(I))*PPB
      END DO

      DLI=-DLI 

	  RETURN
      END

C     --------------------------------------------------------------
      FUNCTION YUS(A,XC)
      PARAMETER (NDHF=12)
      DIMENSION A(NDHF),X(NDHF),P(NDHF),F(NDHF)
      PI=ATAN(1.0)*4.0

      X(1) = 0.03
      DO I=1,NDHF
         X(I)=I/(NDHF+1.)
      END DO

      DO I=1,NDHF
         P(I)=ALOG(0.5)/ALOG(X(I))
      END DO

      F(1)=XC**0.25*(1.-XC)/EXP(20.*XC)
      DO I=2,NDHF-2
         F(I)=(SIN(PI*XC**P(I)))**4
      END DO
      F(NDHF-1)=SIN(PI*XC**P(NDHF-1))
      F(NDHF)=SIN(PI*XC**P(NDHF))

      YUS=THCKAF(XC)
      DO I=1,NDHF
         YUS=YUS+A(I)*F(I)
      END DO

      RETURN
      END
C     ---------------------------------------------------------------
      FUNCTION YLS(A,XC)
      PARAMETER (NDHF=12)
      DIMENSION A(NDHF),X(NDHF),P(NDHF),F(NDHF)
      PI=ATAN(1.0)*4.0

      X(1)=0.03
      DO I=2,NDHF
         X(I)=I/(NDHF+1.)
      END DO

      DO I=1,NDHF
         P(I)=ALOG(0.5)/ALOG(X(I))
      END DO

      F(1)=XC**0.25*(1.-XC)/EXP(20.*XC)
      DO I=2,NDHF-2
         F(I)=(SIN(PI*XC**P(I)))**4
      END DO
      F(NDHF-1)=SIN(PI*XC**P(NDHF-1))
      F(NDHF)=SIN(PI*XC**P(NDHF))
 
      YLS=-THCKAF(XC)
      DO I=1,NDHF
         YLS=YLS+A(I)*F(I)
      END DO

      RETURN
      END
C     --------------------------------------------------------------
      FUNCTION THCKAF(XC)
      a1=1.4779155
      a2=-0.624424
      a3=-1.727016
      a4=1.384087
      a5=-0.510563
      t=0.15
      IF (XC.LT.0..OR.XC.GT.1.) THCKAF=0.
      xc2=xc*xc
      xc3=xc2*xc
      xc4=xc3*xc
      xc5=xc4*xc
      THCKAF=t*(a1*sqrt(xc)+a2*xc+a3*xc2+a4*xc3+a5*xc4)
      RETURN
      END
C     ---------------------------------------------------------------
C================================================================
      SUBROUTINE SHAPE(ALPHIN,DNAME)
      PARAMETER (ND=24)
      DIMENSION ALPHIN(ND)
      CHARACTER DNAME*9
     
      COMMON /DSV/ALPHA(ND),ALPHA0(ND)     

      DIMENSION XSF(201),YSF(201)

      DO I=1,ND
         ALPHA0(I)=ALPHIN(I)
         ALPHA(I)=ALPHA0(I)
      END DO

C     ------ plot design shape --------------------------
C     --- upper half surface ---------
      NDHF=ND/2
      DO I=1,101
         XI=FLOAT(I-1)/FLOAT(101-1)
         XC=0.+1.0*PLYNOM(0.,1.,0.01,0.02,XI)
         XSF(101+I-1)=XC
         YSF(101+I-1)=YUS(ALPHA,XC)
      END DO
C     --- lower half body surface --------
      DO I=1,101
         XI=FLOAT(I-1)/FLOAT(101-1)
         XSF(I)=XSF(201-I+1)
         XC=XSF(I)
         YSF(I)=YLS(ALPHA(NDHF+1),XC)
      END DO
C     -----------------------------------------------
      DO I=1,9
         IF (DNAME(I:I).EQ.CHAR(0)) GOTO 10
      END DO

 10   OPEN (9,FILE=DNAME(:I-1) // ".plt")
      WRITE (9,*) 'VARIABLES=','   X','    Y'
      WRITE (9,*) '  ZONE  ','I=', 201,  '    F=POINT' 
      DO I=1,201 
         WRITE (9,'(2X,2(2X,F19.10))') XSF(I),YSF(I)
      END DO
      CLOSE(9)

      RETURN
      END
C     ------------------------------------------------------------
      FUNCTION PLYNOM(A,B,C,D,ZT)
      A0=A
      A1=C
      A2=3.*(B-A)-(2.*C+D)
      A3=(C+D)-2.*(B-A)
      PLYNOM=A0+A1*ZT+A2*ZT**2+A3*ZT**3
      RETURN
      END
C-------------------------------------------------------------
C-------------------------------------------------------------
      SUBROUTINE EUL2DPLT(ALPHIN,DNAME,RMFIN,ALIN)
C     SUBROUTINE EUL2DPLT(ALPHIN,DNAME)
C------------------------------------------------------------------
C     Numerical solution of 2-d Euler equations using finite volume	
C     method for compressible flows 
C     Modified after SOUBROUTINE RESULT (optm2d04.f)
C     Euler solution cp and field output to files <DNAME>-*.plt
C------------------------------------------------------------------
C     Copyright: ZHANG Zhengke,TSL,NUS, Singapore, Jul.29, 2002
C     Modified : LUM Kai Yew, Apr.12, 20003
C-------------------------------------------------------------------
      PARAMETER (ML=4,MI=161,MJ=91,MM=161)        
      PARAMETER (ND=24)	  
      DIMENSION ALPHIN(ND),GOUT(ND)
      CHARACTER DNAME*9
      REAL RMFIN,ALIN
C     ------ variables,arrays for 2-d Euler eq. ----------      
      COMMON /UVV/ U(ML,MI,MJ),V(ML,MI,MJ)
      COMMON /TSP/ DTL(MI,MJ),SM(MI,MJ),PP(MI,MJ)
      COMMON /FQQ/ FDS(ML,MI,MJ),Q(ML,MI,MJ)
      COMMON /DPP/ DP(MI,MJ),DPM(MI,MJ)
      COMMON /RAD/ RADI(MI,MJ),RADJ(MI,MJ)

      COMMON /GRD/ XX(MI,MJ),YY(MI,MJ)
      COMMON /NEWG/ XNEW(MI,MJ),YNEW(MI,MJ)
      COMMON /GSFI/ SIX(MI,MJ),SIY(MI,MJ)
      COMMON /GSFJ/ SJX(MI,MJ),SJY(MI,MJ)
      COMMON /RXI/ XXI(MI,MJ),YXI(MI,MJ)
      COMMON /RET/ XET(MI,MJ),YET(MI,MJ)
      COMMON /FBNV/ UNX(MI),UNY(MI)
      COMMON /BDOWD/DSTRBD(MI)
      COMMON /CPBD/ CPBD(MI),XBD(MI),YBD(MI)
      COMMON /VOL/ VOL(MI,MJ)       
      COMMON /NGD/ NI,NJ,IHF
      COMMON /NCL/ NI1,NJ1
      
      COMMON /SWT/ INI,ITM,NST,NRSD,NFSV
      COMMON /E22/ E2I,E2J
      COMMON /E44/ E4I,E4J
      COMMON /RSD/ RSDM,ECRI,TSTA,COUR4,EPS,ITL
      COMMON /FAR/ GM,RMF,VXF,VYF,HF,AL
      COMMON /MAX/ IM,JM

C     ---- variables,arrays for gradient etc.---------
      COMMON /DSV/ ALPHA(ND),ALPHA0(ND)
      CHARACTER PDSTR*4
C     ------------------------------------------------
      REAL PBOUT(MI)
C     ------------------------------------------------
      INTEGER INAME
      PI=3.141592653589793

      DO I=1,ND
         ALPHA0(I)=ALPHIN(I)
         ALPHA(I)=ALPHA0(I)
      END DO

      include 'eu.org'
      RMF=RMFIN
      AL=ALIN

      ECRI=0.00001
      ALP=AL*PI/180.

      WRITE(PDSTR,'(2I2.2)') INT(10*RMF+0.5), INT(AL)
C     WRITE(*,*) 'PDSTR=',PDSTR
      
C     ----- read in the grid data of initial shape -----
C     ----- and generate new grid of current shape -----
      OPEN(UNIT=3,NAME='afoil.grd',TYPE='OLD')
      DO J=1,NJ
         DO I=1,NI
            READ(3,*) XX(I,J),YY(I,J)
         END DO
	  END DO
      CLOSE(3)
	  CALL NEWGRD
	  DO J=1,NJ
         DO I=1,NI
            XX(I,J)=XNEW(I,J)
            YY(I,J)=YNEW(I,J)
         END DO
      END DO
  
      CALL EUL2D

C     ---- remove trailing blanks from DNAME ---
      DO INAME=1,9
         IF (DNAME(INAME:INAME).EQ.CHAR(0)) GOTO 10
      END DO
 10   CONTINUE

      CHD=1.0
	  SSB=CHD*1.0
	  XAC=0.0    
C     -------------contribution of pressure--------------------------
      FNX=0.0
      FNY=0.0
      PITCH=0.

      PF=1./1.4/RMF**2

      DO I=1,NI1
         AA=1.+DSTRBD(I)
         BB=-DSTRBD(I)
         PB=AA*PP(I,1)+BB*PP(I,2)
         FNX=FNX-PB*SJX(I,1)
         FNY=FNY-PB*SJY(I,1)
         PITCH=PITCH+(XAC-XBD(I))*(-PB*SJY(I,1))
         CPBD(I)=2.*(PB-PF)
         PBOUT(I)=PB
      END DO
      CN=2.*FNY/SSB
      CA=2.*FNX/SSB
      CM=2.*PITCH/(SSB*CHD)
      CL=CN*COS(ALP)-CA*SIN(ALP)
      CD=CN*SIN(ALP)+CA*COS(ALP)
 
C     ---------------------------------------------------------------
C     output wing surface pressure coefficient distribution
	  XLE=0.0
      OPEN(UNIT=7,NAME=DNAME(:INAME-1)//'-cpu.plt')
      WRITE(7,*) ' VARIABLES=   ','X , ','"-CP"'
      WRITE(7,*) 'ZONE', '  I=',IHF-1, '   F=POINT'    
      DO I=IHF,NI1
         XC=(XBD(I)-XLE)/CHD
         WRITE(7,*) XC,-CPBD(I)
      END DO
      CLOSE(7)

      OPEN(UNIT=7,NAME=DNAME(:INAME-1)//'-cpl.plt')
      WRITE(7,*) '  VARIABLES=   ','X , ','"-CP"'
      WRITE(7,*) 'ZONE', '  I=',IHF-1, '   F=POINT'          
      DO I=1,IHF-1
         XC=(XBD(I)-XLE)/CHD
         WRITE(7,*) XC,-CPBD(I)
	  END DO
      CLOSE(7)

      
      OPEN(UNIT=7,NAME='pd'//PDSTR//DNAME(:INAME-1)//'.org')
      DO I=1,NI1
         WRITE(7,*) PBOUT(I)
      END DO
      CLOSE(7)

C     --------------------------------------------------------

      DO J=1,NJ1
         J1=J+1
         DO I=1,NI1
            I1=I+1
            XG1=XX(I,J)
            XG2=XX(I1,J)
            XG3=XX(I,J1)
            XG4=XX(I1,J1)
            FDS(1,I,J)=(XG1+XG2+XG3+XG4)/4.
            YG1=YY(I,J)
            YG2=YY(I1,J)
            YG3=YY(I,J1)
            YG4=YY(I1,J1)
            FDS(2,I,J)=(YG1+YG2+YG3+YG4)/4. 
         END DO
	  END DO

      DO J=1,NJ1
         DO I=1,NI1
            VX=U(2,I,J)/U(1,I,J)
            VY=U(3,I,J)/U(1,I,J)
            VEL=SQRT(VX**2+VY**2)
            SND=SQRT(GM*PP(I,J)/U(1,I,J))
            RMCH=VEL/SND
            CFFP=2.*(PP(I,J)-PF)
            Q(1,I,J)=VX
            Q(2,I,J)=VY
            Q(3,I,J)=RMCH
            Q(4,I,J)=CFFP
         END DO
	  END DO

      DO J=1,NJ1
         FDS(1,NI,J)=FDS(1,1,J)
         FDS(2,NI,J)=FDS(2,1,J)
         Q(1,NI,J)=Q(1,1,J)
         Q(2,NI,J)=Q(2,1,J)
         Q(3,NI,J)=Q(3,1,J)
         Q(4,NI,J)=Q(4,1,J)
	  END DO


      OPEN(UNIT=2,NAME=DNAME(:INAME-1)//'-fld.plt')
      REWIND(2)
      WRITE(2,*) '  VARIABLES= ','X, ','Y, ', 
     &     'U, ','V, ','Ma, ','P, ','Rho'
      WRITE(2,*) 'ZONE  ','I=',NI,'  J=',NJ1,'  F=POINT'
      DO J=1,NJ1
         DO I=1,NI
            WRITE(2,20) FDS(1,I,J),FDS(2,I,J),
     1           Q(1,I,J),Q(2,I,J),Q(3,I,J),
     2           Q(4,I,J),U(1,I,J)
 20         FORMAT(7(1X,F15.7))
         END DO
      END DO
      CLOSE(2)

      WRITE(*,'(3X,"CD =",F12.8,3X,"CL =",F12.8,3X,"CL/CD =",F12.8)')
     1     CD, CL, CL/CD

      RETURN
      END
C     ---------------------------------------------------------------
