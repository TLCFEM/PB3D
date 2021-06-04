C 
C                        *** PB3D (V1.2) ***
C
C            *   DYNAMIC SIMULATION OF WOOD BUILDINGS   *
C
C              checked work in 2020 with model input file "UBC3S.IN"
C                    UPDATED ON 15TH oct 2020
C  
C            Timber Engineering & Applied Mechanics Group
C                   Department of Wood Science
C                  University of British Columbia
C                  Vancouver, BC, V6T 1Z4  CANADA
C
C
C---------------------------------------------------------------------
C
C     * SHEAR WALLS REPRESENTD BY PSEUDO-NAIL WALL ANALOG MODEL
C     * MAX. NUMBER OF NODES                     = 1000
C     * MAX. NUMBER OF ELEMENTS                  = 1000
C     * MAX. NUMBER OF BOUNDARY CONDITIONS       = 1000 
C     * MAX. NUMBER OF SHEAR WALLS               = 100
C     * MAX. NUMBER OF SHEAR WALL TYPES          = 10
C     * MAX. NUMBER OF FRAME MEMBER TYPES        = 10


      MODULE COMMON_FRAME
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:,:)::NN
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:)::MM,MT
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:)::ATR
      REAL*8,ALLOCATABLE,DIMENSION(:)::DELTAF

      REAL*8,ALLOCATABLE,DIMENSION(:)::EM,GM,BM,HM,
     1XA,XB,XC
      INTEGER,ALLOCATABLE,DIMENSION(:,:)::IEF
      INTEGER,ALLOCATABLE,DIMENSION(:)::MTF,MATID,IETYP,IWALLTMP
      END MODULE COMMON_FRAME

      MODULE COMMON_A
      REAL*8, ALLOCATABLE, DIMENSION(:)::A,A0,AC,AT,AF,AMASS
      END MODULE COMMON_A

      MODULE COMMON_XV
      REAL*8,ALLOCATABLE,DIMENSION(:)::X,XP,VP,ACP,PSI,RMUP,RMU,RC2,RA,
     1R0,PSIU,RFIX
      END MODULE COMMON_XV


C     **************START OF THE MAIN PROGRAM ***********

      USE COMMON_FRAME
      USE COMMON_A
      USE COMMON_XV
      USE DFLIB

      IMPLICIT REAL*8 (A-H, O-Z)

      CHARACTER*100 TITLE,DYNA,ACCELX,ACCELY,STCFOR,STCDIS,
     1FRCFILE,DSPFILE,MSTIFF,CHKFILE,OUTFILE,OUTFILE_A
      CHARACTER*8 XLABEL,YLABEL,ZLABLE
      REAL*8  N0F,N1F,M0F,M1F,M2F,L0F,L1F,L2F,K0F,K1F        !SHAPE FUNCTIONS
      REAL*8  N0SW,N1SW,M0SW,M1SW,M2SW
      REAL*8  N0, N1, M0, M1, M2
      REAL*8  time_begin, time_end
      INTEGER DUMMYINT,ERROR,NUMARG
      CHARACTER*512 FNAME


      COMMON /WALL/                           !SHEAR WALLS
     1N0(100,2,6,9),N1(100,2,6,9),
     1M0(100,2,6,9),M1(100,2,6,9),M2(100,2,6,9),
     1X0SW(100,63),D0PSW(100,20,9),D0NSW(100,20,9),
     1P0PSW(100,20,9),P0NSW(100,20,9),
     1W0PSW(100,20,9),W0NSW(100,20,9),
	1S0SW(100,20,9,16),EPS0SW(100,20,9,16),
     1X0SW2(100,63),D0PSW2(100,20,9),D0NSW2(100,20,9),
     1P0PSW2(100,20,9),P0NSW2(100,20,9),
     1W0PSW2(100,20,9),W0NSW2(100,20,9),
	1S0SW2(100,20,9,16),EPS0SW2(100,20,9,16),
     1FSW0(100),DFSW0(100),
	1FSW2(100),DFSW2(100),DSPSW2(100),
     1FSW(100) ,DFSW(100) ,DSPSW(100),DELTA(100,2),
	1IESW(100,2),IWTYPE(100),IWORIEN(100),IWDMG(100),IWDMG2(100),
     1NELSW,NELSW_1

      COMMON/GAUSSI/                                         !Gauss quadrature
     1EGXF(16),HGXF(16),EGYF(16),HGYF(16),EGZF(16),HGZF(16),
     1EGXW(16),HGXW(16),EGYW(16),HGYW(16),
     1NGXF,NGYF,NGZF,NGXW,NGYW
      COMMON/FSHAPE/                                         !Frame shapes
     1N0F(12,5),N1F(12,5),M0F(12,5),M1F(12,5),M2F(12,5),
     1L0F(12,5),L1F(12,5),L2F(12,5),K0F(12,5),K1F(12,5)
      COMMON /PNAIL/                                         !Pseudo nail and embedment
     1XLSW(10),XDSW(10),ESW(10),SYSW(10),TTC(10),WDMGC(10),
     1Q0SWP(10),Q1SWP(10),Q2SWP(10),Q3SWP(10),Q4SWP(10),
     1XKSWP(10),DMAXSWP(10),PMAXSWP(10),SDFP(10),
     1Q0SWN(10),Q1SWN(10),Q2SWN(10),Q3SWN(10),Q4SWN(10),
     1XKSWN(10),DMAXSWN(10),PMAXSWN(10),SDFN(10)
      COMMON/SUPPORTS/NB(1000),MNB(1000),IB(1000,6),PL(1000,6) !Boundary conditions
      COMMON/SUPL/ NSUPL, ISUPL(1000)                    !Nodes to calculate reaction forces
      COMMON/CONV/ PSI1, PSI0, PSU0, PSU, PSI20, PSI2, FSWT
      COMMON/PLOT/ NEP,NJOTOT,NODMASS(1000),FNODMAS(1000)
      COMMON/GACC/ AXX,AYY,GA

C
      DIMENSION ETA3G(3),H3G(3),ETA4G(4),H4G(4),ETA5G(5),H5G(5),
     1ETA7G(7),H7G(7),ETA9G(9),H9G(9),
     1ETA3L(3),H3L(3),ETA4L(4),H4L(4),ETA5L(5),H5L(5),
     1ETA7L(7),H7L(7),ETA9L(9),H9L(9)
      DIMENSION QW(2,12),QPAR(100,12),XWE(12),PSIW(12),CSW(12,12),
     1N0SW(6,9),N1SW(6,9),M0SW(6,9),M1SW(6,9),M2SW(6,9)
      DIMENSION FACTOR(5000),AX(5000),AY(5000)
      DIMENSION PEAK(2),ADX(10000),ADY(10000),TIM(10000)
      DIMENSION IEO(100), KEO(3)
      DIMENSION XE(12), XTEMP(12), PSIE(12)
      DIMENSION NODEL(10),NODELF(10),IDIS(10),DISPLT(10),
     1PFX(10),PFY(10),PX(10),PY(10)
      DIMENSION CF(12,12)
      DIMENSION NSTEP(3)
      DIMENSION DT00(3)
      DIMENSION XCE(36)
      DIMENSION REA(6)
      DIMENSION NB0(1000),MNB0(1000),IB0(1000,6),PL0(1000,6)
      DIMENSION DISEP0(1000),IEP(1000)
      DIMENSION XLABEL(1000),YLABEL(1000),ZLABEL(1000)
      DIMENSION IDAMWALL(100), IWALLCOUNT(100)
      DIMENSION FYSW(500,100)


C--------------------------------------------------------------------
C     ***** STARTS DATA INPUT ******
C--------------------------------------------------------------------

      INTEGER(4) i2            !MAX WINDOW SIZE
      TYPE (qwinfo) winss
      winss.type = QWIN$MAX
      i2 = setwsizeqq(0, winss)

      WRITE (*,1000)
1000  FORMAT(//
     1'    *************************************************'/
     1'    *                                               *'/
     1'    *                PB3D (V1.2)                    *'/
     1'    *                                               *'/
     1'    *      DYNAMIC SIMULATION OF TIMBER BUILDINGS   *'/
     1'    *                                               *'/
     1'    *                                               *'/
     1'    *          Department of Wood Science           *'/
     1'    *        University of British Columbia         *'/
     1'    *             Vancouver, Canada                 *'/
     1'    *                                               *'/
     1'    *************************************************'/)

      call GET_COMMAND_ARGUMENT(1,DYNA,DUMMYINT,ERROR)
      if (ERROR.GT.0)THEN
         WRITE (*,1041)
1041     FORMAT (1X,'INPUT BUILDING MODEL FILE:' )
         READ (*,*) DYNA
      endif

      OPEN (1,FILE=DYNA,  STATUS='UNKNOWN')      !INPUT DYNAMIC FILE
      OPEN (2,FILE='PB3D-PLOT.OUT', STATUS='UNKNOWN')   !TO PLOT RESULTS (e.g.,deformations)
      OPEN (3,FILE='PB3D-CHECK.OUT', STATUS='UNKNOWN')     !TO CHECK CONVERGENCE
      OPEN (5,FILE='PB3D-WALLS.OUT',STATUS='UNKNOWN')


      READ (1,1001) TITLE
1001  FORMAT(A)
      WRITE (3,1001) TITLE
      WRITE (*,1002) TITLE
1002  FORMAT(/' PROBLEM TITLE:  ',A)
      WRITE (3,1003)
1003  FORMAT('***************************************'/)
C
      READ (1,*) NGXF,NGYF,NGZF           !No of Gaussian integration points for frame members
      CALL GAUSS (NGXF, EGXF, HGXF, IERR) !(Max. NGXF = 9, Max. NGYF = 16, Max. NGZF = 16)
      CALL GAUSS (NGYF, EGYF, HGYF, IERR)
      CALL GAUSS (NGZF, EGZF, HGZF, IERR)
C
      READ (1,*) NGXW,NGYW                !No of Gaussian integration points for 'pseudo nail' in shear wall model
      NGXW=5
      NGYW=9
      CALL GAUSS (NGXW, EGXW, HGXW, IERR) !(Max. NGXW = 5, Max. NGYW = 16)
      CALL GAUSS (NGYW, EGYW, HGYW, IERR)
C
C    * NUNITS:   ENTER 1 IF UNITS ARE IMPERIAL (Lb,in)
C                ENTER 2 IF UNITS ARE METRIC (kN, m)
C                ENTER 3 IF UNITS ARE METRIC (KN, mm)
C
      READ (1,*) NUNITS

      GA = 9810               ! gravity acceleration (in m/s/s)
C
C    * NDYN = 0  STATIC CASE (Load or displacement factor history given in file DISP)               
C                NLOAD = 0 for displacement control
C                NLOAD = 1 for load control
C    * NDYN = 1  EARTHQUAKE EXCITATION, (Accelerograms in files:
C                ACCELX for horizontal accelerations,
C                ACCELY for vertical accelerations.
C                Enter NONE for the file name if no acceleration is applied in the corresponding direction)               
C                DT0 = Time Step in accelerograms (Common to all,will be checked if not common to all).                     
C                DAMP = Damping ratio in %.  
C
      READ (1,*) NDYN
      IF (NDYN.EQ.0) THEN
         READ(1,*) NLOAD
         IF (NLOAD.EQ.1) OPEN(UNIT=4,FILE=FRCFILE,STATUS='OLD')
         IF (NLOAD.EQ.0) OPEN(UNIT=4,FILE=DSPFILE,STATUS='OLD')
      END IF
      IF (NDYN.EQ.1) THEN
         READ (1,*) IACCELX
         READ (1,*) IACCELY
         READ (1,*) DAMP,PERIOD
         READ (1,*) FAA
         IF (IACCELX.EQ.1) READ (1,1001) ACCELX        !ONLY FILE NAME (DO NOT KEY SPACE)
         IF (IACCELY.EQ.1)  READ (1,1001) ACCELY        !ONLY FILE NAME (DO NOT KEY SPACE)
      END IF
C
C    *  NUMBER OF APPLIED LOADS KEPT FIXED DURING THE HISTORY (Static or Dynamic)
C  
      NLFIX=0
C
C   *** FRAME DATA:
C        * NJOTOT = Total number of nodes in the structure
C        * NETOT  = Total number of elements
C        * NMEMF  = Total number of frame members
C        * NMATF  = Total number of different types of framing members (including diagonal braces)
C        * NWALL  = Total number of shear walls    
C        * NWTYP  = Total number of shear wall types  
C      
C
      READ (1,*) NJOTOT,NMEMBM,NMEMBR,NWALL,NMATF,NWTYP
      NETOT=NMEMBM+NMEMBR+NWALL
      NMEMF=NMEMBM+NMEMBR

      OPEN (7,FILE='PB3D-DSP.OUT',   STATUS='UNKNOWN')
      IF (NDYN.EQ.1) THEN
         OPEN (9,FILE='PB3D-ACC.OUT',    STATUS='UNKNOWN')
      ENDIF

      ALLOCATE(NN(NMEMF,5,6,6,12),MT(NMEMF,5,12,12))
      ALLOCATE(ATR(NMEMF,12,12))
      ALLOCATE(DELTAF(NMEMF))
      ALLOCATE(IEF(NMEMF,2),MTF(NMEMF))
      ALLOCATE(EM(NMEMF),GM(NMEMF),BM(NMEMF),
     1HM(NMEMF),MATID(NMEMF))
      ALLOCATE(XA(NJOTOT),XB(NJOTOT),XC(NJOTOT))
      ALLOCATE(IETYP(NETOT))
      ALLOCATE(IWALLTMP(NWTYP))


      IF (NJOTOT.GT.1000) THEN
         WRITE(*,1005)
1005     FORMAT(/' MAX. TOTAL NUMBER OF NODES EXCEEDS 1000!'/)
         GO TO 6000          !GO TO THE END OF THE PROGRAM
      END IF

      IF (NETOT.GT.1000) THEN
         WRITE(*,1006)
1006     FORMAT(/' MAX. TOTAL NUMBER OF ELEMENTS EXCEEDS 1000!'/)
         GO TO 6000          ! GO TO THE END OF THE PROGRAM
      END IF

      IF (NWALL.GT.100) THEN
         WRITE(*,1007)
1007     FORMAT(/' MAX. TOTAL NUMBER OF SHEAR WALLS EXCEEDS 100!'/)
         GO TO 6000          !GO TO THE END OF THE PROGRAM
      END IF

C
C    * Enter global coordinates of nodes
C

      WRITE (2,*) NJOTOT
      DO I=1,NJOTOT                 !Global coordinates (X,Y,Z) for frame nodes.
         READ(1,*) INODE,XA(I),XB(I),XC(I)
         WRITE(2,*) INODE,XA(I),XB(I),XC(I)
      END DO

C
C   * Enter properties for each of frame member types:
C     MATID = Material identifier: 1 for wood, 2 for steel.
C     EM = Modulus of elasticity
C     GM = Shear modulus (converted by program to torsional rigidity)
C     BM = Width of cross-section, or length perpendicular to the plane of the frame
C     HM = Depth of cross-section, or length in the plane of the frame
C
C
      DO  I=1,NMATF
         READ(1,*) MATID(I),EM(I),GM(I),BM(I),HM(I)
         IF (BM(I).LE.HM(I)) THEN
            RAT = 0.3333*(1.0-EXP(-0.55*HM(I)/BM(I)))
            GM(I) = GM(I)*RAT*HM(I)*BM(I)**3         !converted to torsional rigidity
         END IF
         IF (HM(I).LT.BM(I)) THEN
            RAT = 0.3333*(1.0-EXP(-0.55*BM(I)/HM(I)))
            GM(I) = GM(I)*RAT*BM(I)*HM(I)**3         !converted to torsional rigidity
         END IF
      END DO

C              
C    * Enter Pseudo nail wall data:
C    
      NELSW=20              !No of beam elements of pseudo nail
      NELSW_1=18             !No. of beam elements in timber member
      NEQSW=(NELSW+1)*3     !No of total DOFs
      FY0=0.0               !Initialize shear force in the wall
      DO I=1,NWTYP
         READ(1,*) XLSW(I),XDSW(I),ESW(I),SYSW(I),TTC(I),WDMGC(I)  !SUPER NAIL L, D, E, YIELD STRENGTH, WALL ULT. DRIFT
         READ(1,*) Q0SWP(I),Q1SWP(I),Q3SWP(I),XKSWP(I),DMAXSWP(I),
     1   SDFP(I) !EMBEDMENT PARAMETERS
         READ(1,*) Q0SWN(I),Q1SWN(I),Q3SWN(I),XKSWN(I),DMAXSWN(I),
     1   SDFN(I)
         IWALLCOUNT(I)=0                     !A COUNTER TO CHECK WALL HYSTERESIS CONSISTENCY
         Q2SWP(I)=0.8
         Q2SWN(I)=0.8
         PMAXSWP(I)=(Q0SWP(I)+Q1SWP(I)*DMAXSWP(I))
     1   *(1.0-DEXP(-XKSWP(I)*DMAXSWP(I)/Q0SWP(I)))
         PMAXSWN(I)=(Q0SWN(I)+Q1SWN(I)*DMAXSWN(I))
     1   *(1.0-DEXP(-XKSWN(I)*DMAXSWN(I)/Q0SWN(I)))
         Q4SWP(I)=DLOG(Q2SWP(I))/(((Q3SWP(I)-1.0)*DMAXSWP(I))**2)

         Q4SWN(I)=DLOG(Q2SWN(I))/(((Q3SWN(I)-1.0)*DMAXSWN(I))**2)
      ENDDO



C    * Enter framing member connectivity (beams + braces)
C      For each member, calculate transformation matrix and appropriate vectors and matrices
C      IELEMF = Element ID
C      IETYP  = Element type (1: 3D beam; 2: 3D bar; 3: Pseudo-nail wall )
C      IEF(I,1), iEF(I,2) = 1st node ID and 2nd node ID
C      MTF = Frame member type (cross section ID)


      WRITE(2,*) NMEMF
      DO I=1,NMEMF
         READ (1,*) IELEMF,IETYP(I),IEF(I,1),IEF(I,2),MTF(I)
         WRITE(2,*) IEF(I,1),IEF(I,2),IETYP(I)
         DX = XA(IEF(I,2)) - XA(IEF(I,1))
         DY = XB(IEF(I,2)) - XB(IEF(I,1))
         DZ = XC(IEF(I,2)) - XC(IEF(I,1))
         DELTAF(I) = SQRT(DX**2+DY**2+DZ**2)
         CALL SHAPESF(I,NGXF)


         IF(ABS(DX).GT.1.0.OR.ABS(DY).GT.1.0) THEN    !NON VERTICAL MEMBERS
            COSL = DX/DELTAF(I)
            COSM = DY/DELTAF(I)
            COSN = DZ/DELTAF(I)
            DLM=SQRT(COSL**2+COSM**2)
            DO II = 1, 12
               DO J = 1, 12
                  ATR(I,II,J) = 0.0
               END DO
            END DO
            ATR(I,1,1) = COSL
            ATR(I,1,2) = COSM
            ATR(I,1,3) = COSN
            ATR(I,2,1) = -COSM/DLM
            ATR(I,2,2) = COSL/DLM
            ATR(I,2,3) = 0.0
            ATR(I,3,4) = -COSN*COSL/DLM
            ATR(I,3,5) = -COSN*COSM/DLM
            ATR(I,3,6) = DLM
            ATR(I,4,1) = -COSN*COSL/DLM
            ATR(I,4,2) = -COSN*COSM/DLM
            ATR(I,4,3) = DLM
            ATR(I,5,4) = COSM/DLM
            ATR(I,5,5) = -COSL/DLM
            ATR(I,5,6) = 0.0
            ATR(I,6,4) = COSL
            ATR(I,6,5) = COSM
            ATR(I,6,6) = COSN
            DO II = 1, 6
               DO J = 1, 6
                  ATR(I,II+6, J+6) = ATR(I,II,J)
               END DO
            END DO

         ENDIF


         IF(ABS(DX).LT.1.0.AND.ABS(DY).LT.1.0) THEN   ! vertical members
            DO II = 1, 12
               DO J = 1, 12
                  ATR(I,II,J) = 0.0
               END DO
            END DO
            ATR(I,1,3) = 1.0
            ATR(I,2,2) = 1.0
            ATR(I,3,4) = -1.0
            ATR(I,4,1) = -1.0
            ATR(I,5,5) = -1.0
            ATR(I,6,6) = 1.0
            DO II = 1, 6
               DO J = 1, 6
                  ATR(I,II+6, J+6) = ATR(I,II,J)
               END DO
            END DO
         END IF


         DO 48 IGX = 1, NGXF               !obtain relevant vectors
            DO 45 II = 1, 12
               DO 45 J = 1, 12
                  MT(I,IGX,II,J) = 0.0
                  DO 42 K = 1, 12
                     DO 42 L = 1, 12
                        MT(I,IGX,II,J) = MT(I,IGX,II,J) +
     1                  ATR(I,K,II)*K1F(K,IGX)*K1F(L,IGX)*ATR(I,L,J)
 42               CONTINUE
 45         CONTINUE
 48      CONTINUE
         DO 58 IGX = 1, NGXF
            DO 56 II = 1, 12
               A1 = 0.0
               A2 = 0.0
               A3 = 0.0
               DO K = 1,12
                  A1 = A1 + ATR(I,K,II)*N1F(K,IGX)
                  A2 = A2 + ATR(I,K,II)*M2F(K,IGX)
                  A3 = A3 + ATR(I,K,II)*L2F(K,IGX)
               END DO
               DO 52 IGY = 1, NGYF
                  Y = HM(MTF(I))*EGYF(IGY)/2.0
                  DO 51 IGZ = 1, NGZF
                     Z = BM(MTF(I))*EGZF(IGZ)/2.0
                     IF (IETYP(I).EQ.1) NN(I,IGX,IGY,IGZ,II)=
     1               A1-Y*A2-Z*A3         !FOR ROOF/FLOOR BEAMS
                     IF (IETYP(I).EQ.2) NN(I,IGX,IGY,IGZ,II)=A1  !FOR BRACING MEMBERS
                     IF (ABS(DZ).GE.2000)NN(I,IGX,IGY,IGZ,II)=
     1               A1-Y*A2-Z*A3        !FOR WALL POSTS
 51               CONTINUE
 52            CONTINUE
 56         CONTINUE
 58      CONTINUE

      END DO
C
C    *Enter shear wall connectivity, wall type ID and wall orientation
C     (1: X direction, 2: Y direction)
C     Calculate initial stiffness if needed (USE a version of SHYST).

      WRITE(3,1009)
1009  FORMAT(/' INDIVIDUAL SHEAR WALL INITIAL STIFFNESS:')

      WRITE(2,*) NWALL
      DO I=1,NWALL

         IDAMWALL(I)=0        !STATUS INDEX FOR SHEAR WALLS (0=Okay, 1=FAILURE)
         READ (1,*) IWALL, IESW(I,1), IESW(I,2), IWTYPE(I), IWORIEN(I)
         WRITE(2,*) IESW(I,1), IESW(I,2),IWORIEN(I)
         IWTYP=IWTYPE(I)
         DELTA(I,1)=(XLSW(IWTYP)-TTC(IWTYP))/NELSW_1   !ELEMENT LENGTH ON NAIL POINTY SIDE
         DELTA(I,2)=TTC(IWTYP)/(NELSW-NELSW_1)         !ELEMENT LENGTH ON NAIL HEAD SIDE
         DO 180 JL=1,2                                  !TWO EMBEDMENT LAYERS
            CALL SHAPES (N0SW,N1SW,M0SW,M1SW,M2SW,NGXW,EGXW,DELTA(I,JL))
            DO 180 II=1,6
               DO 180 J=1,NGXW
                  N1(I,JL,II,J)=N1SW(II,J)
                  M0(I,JL,II,J)=M0SW(II,J)
                  M1(I,JL,II,J)=M1SW(II,J)
                  M2(I,JL,II,J)=M2SW(II,J)
180      CONTINUE
         DO I1 = 1, 2
            DO J1 = 1, 12
               QW(I1,J1) = 0.0
            END DO
         ENDDO
         QW(1,1) = 1.0
         QW(1,7) = -1.0
         QW(2,2) = 1.0
         QW(2,8) = -1.0
         DO J=1,12
            IF (IWORIEN(I).EQ.1) QPAR(I,J)=QW(1,J)
            IF (IWORIEN(I).EQ.2) QPAR(I,J)=QW(2,J)
         END DO

         DO IX=1,NEQSW
            X0SW2(I,IX)=0.0
         END DO
         DO IX=1,NELSW
            DO IY=1,NGXW
               D0PSW2(I,IX,IY)=0.0
               D0NSW2(I,IX,IY)=0.0
               P0PSW2(I,IX,IY)=0.0
               P0NSW2(I,IX,IY)=0.0
               W0PSW2(I,IX,IY)=0.0
               W0NSW2(I,IX,IY)=0.0
            END DO
         END DO
         DO IX=1,NELSW
            DO IY=1,NGXW
               DO IZ=1,NGYW
                  S0SW2(I,IX,IY,IZ)  =0.0
                  EPS0SW2(I,IX,IY,IZ)=0.0
               END DO
            END DO
         END DO

         CALL SHYSTW (I,IWTYP,0.0D0,1.0D0,FY0,0,IERRORS)
         FSW0(I)=FY0
         DFSW0(I)=FY0/1.0D0
         WRITE(3,1010) I, DFSW0(I)
1010     FORMAT(' WALL No.=',I3,' INITIAL STIFFENESS=',E12.4)

         DO IX=1,NEQSW
            X0SW2(I,IX)=0.0
         END DO
         DO IX=1,NELSW
            DO IY=1,NGXW
               D0PSW2(I,IX,IY)=0.0
               D0NSW2(I,IX,IY)=0.0
               P0PSW2(I,IX,IY)=0.0
               P0NSW2(I,IX,IY)=0.0
               W0PSW2(I,IX,IY)=0.0
               W0NSW2(I,IX,IY)=0.0
            END DO
         END DO
         DO IX=1,NELSW
            DO IY=1,NGXW
               DO IZ=1,NGYW
                  S0SW2(I,IX,IY,IZ)  =0.0
                  EPS0SW2(I,IX,IY,IZ)=0.0
               END DO
            END DO
         END DO

      END DO



C
C    * Allocate arrays size for disp, load and stiffness matrix.
C            
      NEQ=NJOTOT*6
      NA=NEQ*(NEQ+1)/2
      ALLOCATE (X(NEQ),XP(NEQ),VP(NEQ),ACP(NEQ),PSI(NEQ),RMUP(NEQ),
     1RMU(NEQ),RC2(NEQ),RA(NEQ),R0(NEQ),PSIU(NEQ),RFIX(NEQ))
      ALLOCATE (A(NA),A0(NA),AC(NA),AT(NA),AF(NA),AMASS(NA))

C
C
C  * FOR DYNAMICS: 
C                  Enter lumped mass at nodes, calculate the total applied mass FMASST,     
C                  Estimate frequency based on previously calculated stiffness 
C                  (using a ramp displacement)    
C 

      IF (NDYN.EQ.1) THEN
         FMASST = 0.0
         READ (1,*) NMASS
         DO  I=1,NMASS
            READ (1,*) IMASS,NODMASS(I),FMASS
            FNODMAS(I)=FMASS*(1E-6)
            FMASST = FMASST + FNODMAS(I)     !CALCULATE TOTAL MASS
         END DO

C
C        * Plotting the input model 
C
         CALL DRAWMODEL (NMEMF,NWALL,NMASS)


!
!         WRITE(*,1011)
!1011     FORMAT(//' ENTER STRUCTURAL FUNDAMENTAL PERIOD (SEC): ')
!         READ (*,*) PERIOD
         PI=4.0*ATAN(1.0)
         OMEGA  = 2.0*PI/PERIOD
         FREQ   = OMEGA/(2.0*PI)
         PERIOD = 1.0/FREQ
         DAMP_M = 3.0/2.0*(DAMP/100.0)*OMEGA  ! ALPHA (MASS-PROPORTIONAL DAMPING COEFFICIENT)
         DAMP_K = 1.0/2.0*(DAMP/100.0)/OMEGA  ! BETA (STIFFNESS-PROPORTIONAL DAMPING COEFFICIENT)
         WRITE(*,1012) FREQ, PERIOD, DAMP
         WRITE(3,1012) FREQ, PERIOD, DAMP
1012     FORMAT(/' FUNDAMENTAL FREQUENCY = ',F12.3,'(Hz)'/
     1   ' FUNDAMENTAL PERIOD = ',F12.3,' (seconds)'/
     1   ' DAMPING RATIO = ',F12.1,' % '/ )

         DO I = 1, NEQ   !FIX GRAVITY LOADS
            RFIX(I) = 0.0
         END DO

         DO I = 1, NMASS
            II = (NODMASS(I)-1)*6
            RFIX(II+3)=RFIX(II+3)-FNODMAS(I)*GA
         END DO

      END IF
C

      IF (NDYN.EQ.0) THEN
         MEST=1
      END IF
C
C
C  * SUPPORT CONDITIONS:
C     * Enter:
C     * NSUPP0 = Number of total support conditions
C     * NB0  = Frame node where specified
C     * MNB0 = Number of support conditions specified at NB0
C     * IB0  = Array with MNB0 support condition Codes:
C                         1 : X-DIRECTION (Global)
C                         2 : Y-DIRECTION
C                         3 : Z-DIRECTION
C                         4 : ROTATION ABOUT X-AXIS
C                         5 : ROTATION ABOUT Y-AXIS
C                         6 : ROTATION ABOUT Z-AXIS
C     * MTB = 0 If all support displacements and rotations are ZERO
C     * MTB = 1 If at least one of the support conditions is NOT ZERO
C     * IF MTB = 1, PL0 is the array of non-zero specified values.
C
      READ(1,*) NSUPP0
      IF (NSUPP0.GT.1000) WRITE (*,1013)
1013  FORMAT(/' SUPPORT CONDITIONS EXCEED THE NUMBER ALLOWED (1000)'/)

      NSUPP = NSUPP0                  !APPLY FIXED NODES ON BOTTOM
      DO 220 KBC = 1, NSUPP
         READ (1,*) ISUPPORT,NB0(KBC),MNB0(KBC),MTB,(IB0(KBC,I),I=1,6)
         NB(KBC) = NB0(KBC)
         MNB(KBC) = MNB0(KBC)
         K = 0
         DO I = 1, 6
            IF (IB0(KBC,I).NE.0) THEN
               K = K + 1
               IB(KBC,K) = IB0(KBC,I)
            END IF
         END DO
         IF (MTB.NE.0) GO TO 210
         DO I=1,MNB0(KBC)
            PL0(KBC,I) = 0.0
            PL(KBC,I) = 0.0
         END DO
         GO TO 220
 210     READ (1,*) (PL0(KBC,I),I=1,6)
         K = 0
         DO I = 1, 6
            IF(PL0(KBC,I).NE.0) THEN
               K = K + 1
               PL(KBC,K) = PL0(KBC,I)
            END IF
         END DO
 220  CONTINUE

C
C     Static Case:
C     Enter basic displacement or loading pattern, to be amplified
C     by the factors in the history DISP. Displacements are enforced
C     at NODL frame nodes in the X or Y directions.
C     Loads may be enforced at NODL frame nodes (in X or Y directions)
C
 230  IF (NDYN.EQ.0) THEN
         IF (NLOAD.EQ.0) THEN          !DISPLACEMENT ON NODE
            READ (1,*) NODL
            DO 250 I=1,NODL
               READ (1,*) NODEL(I),IDIS(I),DISPLT(I)  !NODE NO., DIRECTION ID,MAG. FACTOR
 250        CONTINUE
         END IF

         IF (NLOAD.EQ.1) THEN          !FORCE LOADS ON NODE
            READ (1,*) NODL
            DO 260 I=1,NODL
               READ (1,*) NODEL(I),PX(I),PY(I)  !NODE NO, FORCE X, FORCE Y.
 260        CONTINUE
         END IF

         READ (4,*) NSTEPS
         IF (NSTEPS.GT.5000) THEN
            WRITE (*,1015)
1015        FORMAT(/' DEMAND FILE HAS MORE THAN 5000 STEPS!'/)
            CLOSE (4)
            STOP
         END IF
         DO K = 1, NSTEPS
            READ (4,*) I, FACTOR(I)
         END DO
         CLOSE (4)
      END IF
C
C     Dynamics (Earthquake) case:
C     Reads accelerograms, in proper units of acceleration.
C     Number of steps NSTEPS must be the same in all files.
C
      IF (NDYN.EQ.1) THEN
         IPEAK = 0
         DO I = 1, 2
            NSTEP(I) = 0
            DT00(I) = 0.0
         END DO
         DIFAM = 0.0
         IF (IACCELX.EQ.1) THEN
            OPEN (UNIT=4,FILE=ACCELX,STATUS='UNKNOWN')
            READ (4,*) NSTEP(1),DT00(1)
            IF (NSTEP(1).GT.5000) THEN
               WRITE (*,1015)
               CLOSE (4)
            END IF
            DO I = 1, NSTEP(1)
               READ (4,*) TIM(I), AX(I)
            END DO
            DO I = 1, NSTEP(1)
               AX(I)=AX(I)*FAA
            END DO
            PEAK(1) = 0.0
            DACXMAX = 0.0
            DO I = 1, NSTEP(1)
               IF (ABS(AX(I)).GT.PEAK(1)) THEN
                  PEAK(1) = DABS(AX(I))
                  IPE = I
               END IF
               IF (I.GT.1) THEN
                  DIFA = DABS(AX(I)-AX(I-1))
                  IF (DIFA.GT.DACXMAX) DACXMAX = DIFA
               END IF
            END DO
            WRITE(*,1016) PEAK(1)
1016        FORMAT(/' PEAK X-ACCELERATION =',E12.4,' (g)')
            IF (IPE.GT.IPEAK) IPEAK = IPE
            IF (DACXMAX.GT.DIFAM) DIFAM = DACXMAX
            CLOSE (4)
         END IF

         IF (IACCELY.EQ.1) THEN
            OPEN (UNIT=4,FILE=ACCELY,STATUS='UNKNOWN')
            READ (4,*) NSTEP(2),DT00(2)
            IF (NSTEP(2).GT.5000) THEN
               WRITE (*,1015)
               CLOSE (4)
               STOP
            END IF
            DO I = 1, NSTEP(2)
               READ(4,*) TIM(I), AY(I)
            END DO
            DO I = 1, NSTEP(2)
               AY(I)=AY(I)*FAA
            END DO
            PEAK(2) = 0.0
            DACYMAX = 0.0
            DO I = 1, NSTEP(2)
               IF (DABS(AY(I)).GT.PEAK(2)) THEN
                  PEAK(2) = DABS(AY(I))
                  IPE = I
               END IF
               IF (I.GT.1) THEN
                  DIFA = DABS(AY(I)-AY(I-1))
                  IF (DIFA.GT.DACYMAX) DACYMAX = DIFA
               END IF
            END DO
            WRITE(*,1017) PEAK(2)
1017        FORMAT(' PEAK Y-ACCELERATION =',E12.4,' (g)'/)
            IF (IPE.GT.IPEAK) IPEAK = IPE
            IF (DACYMAX.GT.DIFAM) DIFAM = DACYMAX
            CLOSE (4)
         END IF

         DO 281 I = 1, 2
            IF (NSTEP(I).EQ.0) GO TO 281
            NSTEPS = NSTEP(I)
            DT0 = DT00(I)
            GO TO 282
 281     CONTINUE
 282     DO 283 I = 1, 2
            IF (NSTEP(I).EQ.0) GO TO 283
            IF (NSTEP(I).NE.NSTEPS.OR.DT00(I).NE.DT0) THEN
               WRITE(*,1018)
1018           FORMAT(/' ACCELEROGRAMS DO NOT SHARE EITHER THE'/
     1         ' NUMBER OF STEPS OR THE TIME STEP!'/)
            END IF
 283     CONTINUE
         IF (IACCELX.EQ.0) THEN
            DO I = 1, NSTEPS
               AX(I) = 0.0
            END DO
         END IF
         IF (IACCELY.EQ.0) THEN
            DO I = 1, NSTEPS
               AY(I) = 0.0
            END DO
         END IF

         DO I = 1, NSTEPS
            AX(I) = AX(I)*GA
            AY(I) = AY(I)*GA
         END DO

         WRITE(*,1019) DIFAM
1019     FORMAT(' MAX ACC STEP SIZE = ',E12.4,' (g)')
!         WRITE(*,1020)
!1020     FORMAT(' ENTER ALLOWABLE MAX ACC STEP SIZE (g)')
!         READ(*,*) ASTEP
         ASTEP=2D-1
         ASTEP = ASTEP*GA

4001     M=0
         DO 4401 I=1,NSTEPS
            NX=0
            NY=0
            IF (I.EQ.1) THEN
               IF (IACCELX.EQ.1) NX=DABS(AX(1))/ASTEP
               IF (IACCELY.EQ.1) NY=DABS(AY(1))/ASTEP
            ELSE
               IF (IACCELX.EQ.1) NX=(DABS(AX(I)-AX(I-1)))/ASTEP
               IF (IACCELY.EQ.1) NY=(DABS(AY(I)-AY(I-1)))/ASTEP
            END IF
            N=NX
            IF (NY.GT.N) N=NY
            N=N+1
            DO 4301 J=1,N
               M=M+1
               IF (M.GT.10000) THEN
                  WRITE (*,1021)
1021              FORMAT(/' RESULTING NUMBER OF POINTS IN THE'/
     1            ' MODIFIED ACCELEROGRAMS EXCEEDS 10000!'/
     1            ' INCREASE THE LARGEST ALLOWED ACCELERATION STEP'/)
                  WRITE (*,1022)
1022              FORMAT(/' ENTER NEW LARGEST ALLOWED ACC STEP (g)'/)
                  READ (*,*) ASTEP
                  ASTEP = ASTEP*GA
                  GO TO 4001
               END IF
               TIM(M)=(I-1)*DT0+(J*DT0)/N
               IF (I.EQ.1) THEN
                  IF (IACCELX.EQ.1) ADX(M)=(AX(1)*J)/N
                  IF (IACCELY.EQ.1) ADY(M)=(AY(1)*J)/N
               ELSE
                  IF (IACCELX.EQ.1) ADX(M)=AX(I-1)+((AX(I)-AX(I-1))*J)/N
                  IF (IACCELY.EQ.1) ADY(M)=AY(I-1)+((AY(I)-AY(I-1))*J)/N
               END IF
 4301       CONTINUE
 4401    CONTINUE
         WRITE (*,1023) NSTEPS,M
1023     FORMAT(/' NUMBER OF STEPS IN ORIGINAL ACCELEROGRAMS =',I6/
     1   ' NUMBER OF STEPS IN MODIFIED ACCELEROGRAMS =',I6/)
!         WRITE (*,1024)
!1024     FORMAT(' IF MODIFIED NUMBER OF STEPS IS OK, ENTER 1'/
!     1   ' IF NOT, ENTER 0')
!         READ (*,*) NANS
!         IF (NANS.EQ.1) GO TO 4701
!         WRITE (*,1022)
!         READ (*,*) ASTEP
!         ASTEP=ASTEP*GA
!         GO TO 4001
 4701    NSTEPS = M
         IF (IACCELX.EQ.1) THEN
            OPEN (UNIT=11,FILE='AXNEW.ACC',STATUS='UNKNOWN')
            WRITE (11,1025) NSTEPS
1025        FORMAT(I5)
            DO 5001 I=1,NSTEPS
               WRITE (11,1026) TIM(I),ADX(I)
1026           FORMAT(F12.3,F12.5)
 5001       CONTINUE
            CLOSE (11,STATUS='DELETE')
         END IF
         IF (IACCELY.EQ.1) THEN
            OPEN (UNIT=12,FILE='AYNEW.ACC',STATUS='UNKNOWN')
            WRITE (12,1025) NSTEPS
            DO 5101 I=1,NSTEPS
               WRITE (12,1026) TIM(I),ADY(I)
 5101       CONTINUE
            CLOSE (12,STATUS='DELETE')
         END IF
         IF (IACCELX.EQ.0) THEN
            DO I = 1, NSTEPS
               ADX(I) = 0.0
            END DO
         END IF
         IF (IACCELY.EQ.0) THEN
            DO I = 1, NSTEPS
               ADY(I) = 0.0
            END DO
         END IF


      END IF

C
c
C  * OUTPUT DATA:
C     NEO = Number of nodes in the output file
C     IEO = Node numbers included in the output
C     KEO = Indicator for displacement direction:
C           KEO(1) = 1 to include displacement for global X direction,0 to skip                   
C           KEO(2) = 1 to include displacement for global Y direction,0 to skip                   
C           KEO(3) = 1 to include displacement for global Z direction,0 to skip                   
C     IEP  = Node number for which data histories are stored
C            for plotting, also used to compute approximate frequency in dynamics    
C
      READ (1,*) NSUPL       ! NO OF SUPPORT NODES TO CALCULATE REACTION FORCES
      READ (1,*) (ISUPL(I),I=1,NSUPL)

      WRITE (2,*) NSTEPS

      NEO=1
      IEO(1)=1
      KEO(1)=1
      KEO(2)=1
      KEO(3)=1

      NEP=NJOTOT               !NO OF TOTAL NODES IN OUTPUT FILE
C
C    * HEADERS FOR DISP OUTPUT FILE
C
      DO I = 1, NEP
         IEP(I)=I
         DISEP0(I)=0.0
         XLABEL(I)=' UX(mm)'          !CHARACTERS X AND Y FOR HEADERS IN OUTPUT FILE
         YLABEL(I)=' UY(mm)'
         ZLABEL(I)=' NODE'
      END DO


      WRITE(7,1060) (ZLABEL(I),I,I=1,NEP)           !HEADERS FOR THE OUTPUT FILE
1060  FORMAT(' *TIME HISTORY OF DISP* ',1000(A,I4))
      WRITE(7,1065)  (XLABEL(I),YLABEL(I),I=1,NEP)
1065  FORMAT(' NS/TIME RXT(kN) RYT(kN) RZT(kN) ',1000A8)
      WRITE(7,1075)0.0,0.0,0.0,(DISEP0(I),DISEP0(I),I=1,NEP)
1075  FORMAT(1000F5.1)


C
C    * HEADERS FOR ACC OUTPUT FILE
C

      IF (NDYN.EQ.1) THEN
         DO I = 1, NEP
            IEP(I)=I
            DISEP0(I)=0.0
            XLABEL(I)=' AccX(g)'          !CHARACTERS X AND Y FOR HEADERS IN OUTPUT FILE
            YLABEL(I)=' AccY(g)'
         END DO
         WRITE(9,1066) (ZLABEL(I),I,I=1,NEP)        !HEADERS FOR THE OUTPUT FILE
1066     FORMAT(' *TIME HISTORY OF ACC* ',1000(A,I4))
         WRITE(9,1065)  (XLABEL(I),YLABEL(I),I=1,NEP)
         WRITE(9,1075)0.0,0.0,0.0,(DISEP0(I),DISEP0(I),I=1,NEP)
      ENDIF

C
C    * HEADERS FOR FOR WALL LOAD-DRIFT OUTPUT FILE 
C		

      DO I = 1, NWALL
         IEP(I)=I
         XLABEL(I)='mm'          !CHARACTERS X AND Y FOR HEADERS IN OUTPUT FILE
         YLABEL(I)='KN'
         ZLABEL(I)='DMG'
      END DO
      WRITE (5,1061) (IEP(I),IEP(I),IEP(I),I=1,NWALL)
1061  FORMAT(' WALL_NO:',1000I4 )
      WRITE (5,1062) (XLABEL(I),YLABEL(I),ZLABEL(I),I=1,NWALL)
1062  FORMAT(' NS/TIME ',1000A10)


C
C  * ENTER GLOBAL TOLERANCES:
C       TOLGFR = RELATIVE CONVERGENCE FOR PSI VECTOR, OR OUT-OF-BALANCE VECTOR (e.g., 1.0E-02)                       
C       TOLGFA = ABSOLUTE CONVERGENCE FOR PSI VECTOR, OR OUT-OF-BALANCE VECTOR (e.g., 1.0E-05)                       
C       TOLGXR = RELATIVE CONVERGENCE FOR VECTOR INCREMENTIN X (e.g., 1.0E-05)                     
C       TOLGXA = ABSOLUTE CONVERGENCE FOR VECTOR INCREMENT IN X (e.g., 1.0E-02)                        
C       TOLGUR = RELATIVE CONVERGENCE IN ENERGY (e.g., 1.0E-05)                       
C       TOLGUA = ABSOLUTE CONVERGENCE IN ENERGY (e.g., 1.0E-02)
C     NOTE: THE ABSOLUTE TOLERANCES ENTERED THROUGH THE GRAPHICS
C           INTERFACE ARE OVERWRITTEN AS BELOW
C           Tolerance in force is 1 N ( or 1.0D-03 kN)
C           Tolerance in displacement is 0.5 mm
C           Accordingly, tolerance in energy is 0.5D-03
C           These tolerances applied to the other cases of NUNITS
C

      READ(1,*) TOLGFR,TOLGXR,TOLGUR
      READ(1,*) TOLGFA,TOLGXA,TOLGUA
      CLOSE (1)
C
C   ****** END OF DATA INPUT ******
C------------------------------------------------------------------
C------------------------------------------------------------------
C------------------------------------------------------------------
C------------------------------------------------------------------
C------------------------------------------------------------------ 
C------------------------------------------------------------------
C------------------------------------------------------------------
C------------------------------------------------------------------
C    ****** START SOLUTION *******
C
C

      CALL CPU_TIME(time_begin)

      IF (NDYN.EQ.0.AND.NLOAD.EQ.0) WRITE(*,1028) NSTEPS
      IF (NDYN.EQ.0.AND.NLOAD.NE.0) WRITE(*,1029) NSTEPS
      IF (NDYN.NE.0) WRITE(*,1030) NSTEPS
1028  FORMAT(/' ***STARTS SOLUTION:'/
     1' LOAD HISTORY HAS ',I6,' STEPS (Statics,',
     1' Displacement Control)')
1029  FORMAT(/' ***STARTS SOLUTION:'/
     1' LOAD HISTORY HAS ',I6,' STEPS (Statics,',
     1' Load Control)')
1030  FORMAT(/' ***STARTS SOLUTION:'/
     1' LOAD HISTORY HAS ',I6,' STEPS (Earthquake)')

C
C    * Calculate mass matrix in the case of dynamics
C  
      IF (NDYN.EQ.1) THEN
         DO I = 1, NA
            AMASS(I) = 0.0
            A(I)=0.0
         END DO
         DO IM=1,NMASS
            I1 = (NODMASS(IM)-1)*6 + 1
            I2 = (NODMASS(IM)-1)*6 + 2
            I3 = (NODMASS(IM)-1)*6 + 3
            IJ1 = ((2*NEQ-I1)*(I1-1))/2 + I1
            IJ2 = ((2*NEQ-I2)*(I2-1))/2 + I2
            IJ3 = ((2*NEQ-I2)*(I2-1))/2 + I3
            AMASS(IJ1) = AMASS(IJ1) + FNODMAS(IM)
            AMASS(IJ2) = AMASS(IJ2) + FNODMAS(IM)
!	      AMASS(IJ3) = AMASS(IJ3) + FNODMAS(IM)
         END DO
      END IF
C
C *  For static case, if loads are specified, build basic load
C    vector {R0}, to be multiplied later by the load history factor
C
      IF (NDYN.EQ.0) THEN
         DO I=1,NEQ
            R0(I)=0.0
         END DO
         IF (NLOAD.EQ.1) THEN
            DO I=1,NODL
               II=(NODEL(I)-1)*6
               R0(II+1)=R0(II+1)+PX(I)
               R0(II+2)=R0(II+2)+PY(I)
            END DO
         END IF
      END IF
C
      IF (NDYN.EQ.0) THEN
         FACT0=0.0
      ELSE IF (NDYN.EQ.1) THEN
         AXX0=0.0
         AYY0=0.0
      END IF


C
C    * ESTABLISH LINEAR ELASTIC STIFFNESS MATRIX FOR BEAM/BRACE MEMBERS (1: for beam, 2: for brace)
C
      DO I = 1, NA
         A(I)=0.0
         AF(I) = 0.0
      END DO

      DO JF = 1, NMEMF
         DELE = DELTAF(JF)
         BB = BM(MTF(JF))
         HH = HM(MTF(JF))
         E = EM(MTF(JF))
         MATD = MATID(MTF(JF))
         IETYPE=IETYP(JF)
         DO I = 1, 12
            DO J = 1, 12
               CF(I,J) = 0.0D0
            END DO
         END DO

         IF (IETYPE.EQ.1) THEN

            DO IGX = 1, NGXF

               DO I = 1, 12
                  DO J = 1, 12
                     CF(I,J) = CF(I,J) +
     1               GM(MTF(JF))*MT(JF,IGX,I,J)*HGXF(IGX)*
     1               DELE/2.0
                  END DO
               END DO

               DO  IGY = 1, NGYF
                  DO  IGZ = 1, NGZF
                     DO I = 1, 12
                        DO J = 1, 12
                           CF(I,J) = CF(I,J)+E*NN(JF,IGX,IGY,IGZ,I)*
     1                     NN(JF,IGX,IGY,IGZ,J)*HGXF(IGX)*HGYF(IGY)
     1                     *HGZF(IGZ)
     1                     *DELE*BB*HH/8.0
                        END DO
                     END DO
                  ENDDO
               ENDDO

            ENDDO
         END IF

         IF (IETYPE.EQ.2) THEN
            DO  IGX = 1, NGXF
               DO IGY = 1, NGYF
                  DO IGZ = 1, NGZF
                     DO I = 1, 12
                        DO J = 1, 12
                           CF(I,J) = CF(I,J)+E*NN(JF,IGX,IGY,IGZ,I)*
     1                     NN(JF,IGX,IGY,IGZ,J)*HGXF(IGX)*HGYF(IGY)*
     1                     HGZF(IGZ)
     1                     *DELE*BB*HH/8.0
                        END DO
                     END DO

                  ENDDO
               ENDDO
            ENDDO
         END IF

         DO I = 1, 12
            DO J = 1, 12
               IF(ABS(CF(I,J)).LE.1E-5) CF(I,J)=0.0
            END DO
         END DO

         CALL GLOBMF (JF, NEQ, CF)

      END DO


      DO I = 1, NA
         AF(I)=A(I)            !CONSTANT FRAME/BRACE STIFFNESS MATRIX
      END DO

C
C    * Zeroes out initial stiffness matrix (for dynamic case)
C
      DO I = 1, NA
         A(I)=0.0
         A(I)=A(I)+AF(I)
         AT(I) = 0.0
      END DO
C
C    * Zeroes out initial displacement, velocity and acceleration vectors
C
      DO I=1,NEQ
         XP(I)=0.0
         IF (NDYN.EQ.1) THEN
            VP(I)=0.0
            ACP(I)=0.0
         END IF
      END DO
C
C    * Zeroes out initial vectors for frame member properties
C
      DO JF = 1, NMEMF
         MATD = MATID(MTF(JF))
      END DO
C
C    * Zeroes out initial vectors for pseudo nail walls 
C
      DO KC=1,NWALL
         DSPSW(KC)= 0.0          !WALL DRIFT
         FSW(KC)  = 0.0          !LATERAL FORCE
         DFSW(KC) = DFSW0(KC)    !LATERAL STIFFNESS
         DSPSW2(KC)= 0.0          !WALL DRIFT
         FSW2(KC)  = 0.0          !LATERAL FORCE
         DFSW2(KC) = DFSW0(KC)    !LATERAL STIFFNESS
         IWDMG(KC)= 0
         IWDMG2(KC)= 0            !WALL DAMAGE INDICATOR 0=UNDAMAGED
         DO I=1,NEQSW
            X0SW(KC,I)=0.0
         END DO
         DO IX=1,NELSW
            DO I=1,NGXW
               D0PSW(KC,IX,I)=0.0
               D0NSW(KC,IX,I)=0.0
               P0PSW(KC,IX,I)=0.0
               P0NSW(KC,IX,I)=0.0
               W0PSW(KC,IX,I)=0.0
               W0NSW(KC,IX,I)=0.0
               DO J=1,NGYW
                  S0SW(KC,IX,I,J)=0.0D0
                  EPS0SW(KC,IX,I,J)=0.0D0
               END DO
            END DO
         END DO
      END DO
C               
C
C
C    **** STARTS ITERATIONS OVER SPECIFIED STEPS ****
C
C

      DO 5000 NS=1,NSTEPS                    ! INTERATION OVER STEPS


         IF (NDYN.EQ.0) THEN
            FACT=FACTOR(NS)
         END IF
         IF (NDYN.EQ.1) THEN
            AXX=ADX(NS)
            AYY=ADY(NS)
            IF (NS.EQ.1) DT1 = TIM(1)
            IF (NS.GT.1) DT1 = TIM(NS)-TIM(NS-1)
            DT = DT1
         END IF
C
C **** Begins constructing the vector PSI and the matrix A
C
         NDET=0
         NDETT=0
         NITER=0

501      DO I = 1, NEQ
            X(I) = XP(I)
         END DO

502      IF (NDYN.EQ.1) THEN
            DO I = 1, NA
               AC(I)= DAMP_M*AMASS(I)+DAMP_K*AT(I)  ! Rayleigh damping matrix
               A0(I)= 4.0*AMASS(I)/(DT**2)+2.0*AC(I)/DT
            END DO
         END IF

         IF (NDYN.EQ.0) THEN   !Zeroes out stiffness matrix and nodal force vector
            DO I = 1, NA
               A(I) = AF(I)
            END DO
            DO I = 1, NEQ
               PSI(I) = 0.0D0
               RC2(I) = 0.0D0
            END DO
         END IF

         IF (NDYN.EQ.1) THEN   !Zeroes out stiffness matrix and nodal force vector
            DO I = 1, NA
               A(I) = AF(I)
            END DO
            DO I = 1, NEQ
               PSI(I) = 0.0D0
               RC2(I) = 0.0D0
            END DO
         END IF


C
C    * Enter the contribution from the frame members (1: for beam, 2: for brace)
C
         DO JF = 1, NMEMF

            DELE = DELTAF(JF)
            BB = BM(MTF(JF))
            HH = HM(MTF(JF))
            E = EM(MTF(JF))
            MATD = MATID(MTF(JF))
            IETYPE=IETYP(JF)
            DO I = 1, 12
               PSIE(I) = 0.0
            END DO
            DO I = 1, 6
               II = (IEF(JF,1)-1)*6 + I
               JJ = (IEF(JF,2)-1)*6 + I
               XE(I) = X(II)
               XE(I+6) = X(JJ)
            END DO

            IF (IETYPE.EQ.1) THEN
               DO 332 IGX = 1, NGXF
                  DO I = 1, 12
                     XMT = 0.0
                     DO K = 1, 12
                        XMT =XMT + MT(JF,IGX,I,K)*XE(K)      ! TORSION component
                     END DO
                     PSIE(I) = PSIE(I)+GM(MTF(JF))*XMT*HGXF(IGX)*DELE/2.0
                  END DO

                  DO 331 IGY = 1, NGYF
                     DO 330 IGZ = 1, NGZF
                        EPS = 0.0
                        DO I = 1, 12
                           EPS = EPS + NN(JF,IGX,IGY,IGZ,I)*XE(I)      ! AXIAL LOADING
                        END DO
                        S = E * EPS         !Normal stress
                        DO I = 1, 12
                           PSIE(I) = PSIE(I) + S*NN(JF,IGX,IGY,IGZ,I)*
     1                     HGXF(IGX)*HGYF(IGY)*HGZF(IGZ)*DELE*BB*HH/8.0
                        END DO
330                  CONTINUE
331               CONTINUE

332            CONTINUE
            END IF

            IF (IETYPE.EQ.2) THEN
               DO 333 IGX = 1, NGXF
                  DO 334 IGY = 1, NGYF
                     DO 335 IGZ = 1, NGZF
                        EPS = 0.0
                        DO I = 1, 12
                           EPS = EPS + NN(JF,IGX,IGY,IGZ,I)*XE(I)
                        END DO
                        S = E * EPS
                        DO I = 1, 12
                           PSIE(I) = PSIE(I) + S*NN(JF,IGX,IGY,IGZ,I)*
     1                     HGXF(IGX)*HGYF(IGY)*HGZF(IGZ)*DELE*BB*HH/8.0
                        END DO
335                  CONTINUE
334               CONTINUE
333            CONTINUE

            END IF

            CALL GLOBRF (JF, PSIE)

         END DO
C
C
C    * ENTER CONTRIBUTIONS OF SHEAR WALLS
C
C
485      DO KC=1,NWALL
            NI=IESW(KC,1)           ! 1ST NODE CONNECTING SHEAR SPRING
            NJ=IESW(KC,2)           ! 2RD NODE CONNECTING SHEAR SPRING
            IWTYP=IWTYPE(KC)         ! WALL TYPE
            IWORIN=IWORIEN(KC)

            IF (IWDMG(KC).EQ.1) THEN     !IF SHEAR WALL FAILS
               FW=0.0
               DFW=0.0
               GOTO 370
            ENDIF

            DO I =1,12
               PSIW(I)=0.0
               DO J=1,12
                  CSW(I,J)=0.0
               END DO
            END DO
            IF (IWORIN.EQ.1) DSPW=X((NI-1)*6+1)-X((NJ-1)*6+1)
            IF (IWORIN.EQ.2) DSPW=X((NI-1)*6+2)-X((NJ-1)*6+2)

C
C   * Calculate the lateral hysteretic forces, starting from the last converged values.
C  	 

            DO IX=1,NEQSW
               X0SW2(KC,IX)=X0SW(KC,IX)
            END DO
            DO IX=1,NELSW
               DO IY=1,NGXW
                  D0PSW2(KC,IX,IY)=D0PSW(KC,IX,IY)
                  D0NSW2(KC,IX,IY)=D0NSW(KC,IX,IY)
                  P0PSW2(KC,IX,IY)=P0PSW(KC,IX,IY)
                  P0NSW2(KC,IX,IY)=P0NSW(KC,IX,IY)
                  W0PSW2(KC,IX,IY)=W0PSW(KC,IX,IY)
                  W0NSW2(KC,IX,IY)=W0NSW(KC,IX,IY)
               END DO
            END DO
            DO IX=1,NELSW
               DO IY=1,NGXW
                  DO IZ=1,NGYW
                     S0SW2(KC,IX,IY,IZ)  =S0SW(KC,IX,IY,IZ)
                     EPS0SW2(KC,IX,IY,IZ)=EPS0SW(KC,IX,IY,IZ)
                  END DO
               END DO
            END DO

            IWDMG2(KC)=IWDMG(KC)   !ASSIGN DAMAGE INDICATOR (1=failure; 0=working)

            IF (NITER.EQ.0) THEN
               DFW=DFSW(KC)
               FW =FSW(KC)
               GOTO 370
            END IF

C
            IF (NITER.GE.1) THEN

               IF (DABS(DSPW).GE.WDMGC(IWTYP).and.
     1         DABS(DSPSW(KC)).LT.WDMGC(IWTYP)) THEN     !RUN AN EXTRA LOOP IF WALL FAILURE CRITERIA IS EXCEEDED
                  IWDMG2(KC)=1
                  FW=0.0
                  DFW=0.0
                  GOTO 370
               ELSE
                  CALL SHYSTW(KC,IWTYP,DSPSW(KC),DSPW,FW,IWDMG(KC),
     1            IERRORS)
               ENDIF

               DENOM  = DSPW-DSPSW(KC)
               XNUMEN = FW-FSW(KC)
               DFW    = XNUMEN/DENOM
               TOLX  = 1E-6
               IF(ABS(DENOM).LT.TOLX) THEN
                  DSPW= DSPSW(KC)
                  FW  = FSW(KC)
                  DFW = DFSW(KC)
               END IF

            END IF
C
370      DSPSW2(KC) = DSPW
         FSW2(KC) = FW
         DFSW2(KC) = DFW

380      DO I = 1, 12
            PSIW(I) = PSIW(I) + FW*QPAR(KC,I)
            DO J = 1, 12
               CSW(I,J)=CSW(I,J)+DFW*QPAR(KC,I)*QPAR(KC,J)
            END DO
         END DO
         DO J = 1, 6
            PSI((NI-1)*6+J) = PSI((NI-1)*6+J) + PSIW(J)
            PSI((NJ-1)*6+J) = PSI((NJ-1)*6+J) + PSIW(J+6)
            RC2((NI-1)*6+J) = RC2((NI-1)*6+J) + PSIW(J)
            RC2((NJ-1)*6+J) = RC2((NJ-1)*6+J) + PSIW(J+6)
         END DO
         CALL GLOBSW(NI,NJ,NEQ,CSW)

         END DO


c    *  For each step, set tangent stiffness matrix for Rayleigh damping (stiffness proportional)

         IF (NDYN.EQ.1) THEN
            DO I = 1, NA
               AT(I) = A(I)
            END DO
         ENDIF

C
C    *  construct the out-of-balance force vector
C

390      IF (NDYN.EQ.0) THEN
            DO I=1,NEQ
               PSI(I) = PSI(I)-R0(I)*FACT
               RC2(I) = RC2(I)
            END DO
         END IF

         IF (NDYN.EQ.1) THEN
            CALL GMULT (NEQ,1.0D0,1.0D0,DT,1)
            CALL GMULT (NEQ,FDYN1,FDYN2,DT,2)
            DO I = 1, NEQ
               RA(I) = 0.0
            END DO
            DO I=1,NMASS
               NOD=NODMASS(I)
               XMA=FNODMAS(I)
               CALL MASSNOD (NOD, XMA, AXX, AYY)    !Inertial forces
            END DO
            DO I=1,NEQ
               PSI(I) = PSI(I)+ RMU(I)+RA(I)-RMUP(I)-RFIX(I) !Unbalanced force vector
               RC2(I) = RC2(I)+ RMU(I)+RA(I)-RMUP(I)-RFIX(I)
            END DO
         END IF

         IF (NDYN.EQ.1) THEN
            DO I = 1, NA
               A(I) = A(I) + A0(I)
            END DO
         END IF
C
C
C *  Move PSI to right-hand side of the equation
C
         DO I = 1, NEQ
            PSI(I) = -PSI(I)
         END DO
C
C
C * Enforces support conditions and prescribed displacements for
C   static case (NDYN=0, NLOAD=0).
C
487      IF (NDYN.EQ.0.AND.NLOAD.EQ.0) THEN
            DO II = 1, NODL
               IF (IDIS(II).EQ.1) THEN
                  IBC = (NODEL(II)-1)*6 + 1
                  PLL = DISPLT(II)*FACT
                  CALL BC1(IBC,PLL,NEQ,1)
                  CALL BC2(IBC, NEQ)
               END IF
               IF (IDIS(II).EQ.2) THEN
                  IBC = (NODEL(II)-1)*6 + 2
                  PLL = DISPLT(II)*FACT
                  CALL BC1(IBC,PLL,NEQ,1)
                  CALL BC2(IBC, NEQ)
               END IF
            END DO
         END IF
C
         DO II=1,NJOTOT
            DO I=1,MNB(II)
               PLL = PL(II,I)
               IBC=(NB(II)-1)*6+IB(II,I)
               CALL BC1 (IBC, PLL, NEQ,1)
               CALL BC2 (IBC, NEQ)
            END DO
         END DO
C
         DO II=1,NJOTOT
            DO I=1,MNB(II)
               PLL = PL(II,I)
               IBC=(NB(II)-1)*6+IB(II,I)
               CALL BC1(IBC,PLL,NEQ,2)
            END DO
         END DO
C
         IF (NDYN.EQ.0.AND.NLOAD.EQ.0) THEN
            DO II=1,NODL
               IF (IDIS(II).EQ.1) THEN
                  IBC=(NODEL(II)-1)*6+ 1
                  PLL = DISPLT(II)*FACT
                  CALL BC1(IBC,PLL,NEQ,2)
               END IF
               IF (IDIS(II).EQ.2) THEN
                  IBC = (NODEL(II)-1)*6 + 2
                  PLL = DISPLT(II)*FACT
                  CALL BC1(IBC,PLL,NEQ,2)
               END IF
            END DO
         END IF

C
C  * Decomposes the matrix A
C   
         CALL DECOMP (NEQ, NEQ, IERROR)
         IF (IERROR.EQ.1)THEN
            WRITE(*,*) '****** MATRIX DECOMPOSITION FAILS ! ******'
         END IF
         IF (IERROR.EQ.1) GO TO 400
C
C * Solves the system of equations
C
         DO I = 1, NEQ
            PSIU(I) = PSI(I)
         END DO
C	  
         CALL SOLV (NEQ, NEQ)
C	  
C
C * For NITER > 0, Compute norm of the out-of-balance vector for the current {x}  
C   and the norm of the increment in {x}
C
         PSI1 = 0.0
         DO I = 1, NEQ
            PSI1 = PSI1 + PSIU(I)**2
         END DO
         PSI1 = DSQRT(PSI1)

         PSI2 = 0.0
         DO I = 1, NEQ
            PSI2 = PSI2 + PSI(I)**2
         END DO
         PSI2 = DSQRT(PSI2)

         PSU = PSI1*PSI2

         IF (NDET.EQ.0.AND.NITER.EQ.0) THEN
            PSI0=PSI1
            PSI20=PSI2
            PSU0=PSU
         ENDIF
C
C    * Checks convergence, and updates vector {x} if convergence not satisfied.       
C
         WRITE(*,1034) NS,NDETT,NDET,NITER,PSI1,PSI2,PSU
1034     FORMAT(/' STEP=',I5,' INTERMEDIATE CONVERGENCES=',I3/
     1   ' CURRENT INTERMEDIATE SPLITS= ',I2,' ITER= ',I2/
     1   ' PSI1=',E12.5,' PSI2=',E12.5,' PSU=',E12.5)
C  	 
         KC1 = 0
         KC2 = 0
         KC3 = 0
         IF (PSI1.LT.TOLGFA) KC1 = 1
         IF (PSI2.LT.TOLGXA) KC2 = 1
         IF (PSU.LT.TOLGUA) KC3=1
         IF (KC1.EQ.1.OR.KC2.EQ.1.OR.KC3.EQ.1) GO TO 410
C
         DO I=1,NEQ
            X(I)=X(I)+PSI(I)
         END DO
         NITER = NITER + 1
         IF (NITER.GT.10) GO TO 400   ! max. equilibrium iterations = 10
         GO TO 502
C
C *  Splits step, reduces the right-hand side (max. = 10)
C
400      NDET=NDET+1
         NITER=0
         IF (NDET.LE.10) THEN
            IF (NDYN.EQ.0) FACT=(FACT+FACT0)/2.0
            IF (NDYN.EQ.1) THEN
               DT = DT/2.0
               AXX=(AXX+AXX0)/2.0
               AYY=(AYY+AYY0)/2.0
            END IF
            GO TO 501

         END IF

         IF (NDET.GT.10) THEN
            WRITE(*,1035) NS, NDET-1
1035        FORMAT(' AT STEP ',I5,':'/
     1      ' NO SOLUTION HAS BEEN ACHIEVED AFTER ',I2,' SPLITS'//
     1      ' *** PROGRAM STOPS ***'/)
            GO TO 6000
         END IF
C
C * Convergence achieved (intermediate or final step).
C    For dynamics, updates displacements, velocities and accelerations.
C    In general, updates connector and member hysteresis parameters.
C    Go to next step.
C
C
410      IF (NDYN.EQ.1) THEN
            DO I=1,NEQ
               ACP(I)=(X(I)-VP(I)*DT-XP(I))*4.0/(DT**2)-ACP(I)
               VP(I)=(X(I)-XP(I))*2.0D0/DT-VP(I)
            END DO
         END IF

         DO I=1,NEQ
            XP(I)=X(I)
         END DO

C 
C  *** UPDATING SHEAR WALLS 
C
         IF (NDYN.EQ.1) THEN
            DO KC=1,NWALL
               NI=IESW(KC,1)           ! 1ST NODE CONNECTING SHEAR SPRING
               NJ=IESW(KC,2)           ! 2RD NODE CONNECTING SHEAR SPRING
               IWTYP=IWTYPE(KC)         ! WALL TYPE
               IWORIN=IWORIEN(KC)
               IF (IWORIN.EQ.1) DSPW=ABS(X((NI-1)*6+1)-X((NJ-1)*6+1))
               IF (IWORIN.EQ.2) DSPW=ABS(X((NI-1)*6+2)-X((NJ-1)*6+2))
               IF (DSPW.GT.WDMGC(IWTYP)) THEN
                  IDAMWALL(KC)=1
               ENDIF
            ENDDO
         ENDIF

         FSWT=0.0
         DO KC=1,NWALL
            DSPSW(KC)= DSPSW2(KC)
            FSW(KC)  = FSW2(KC)
            DFSW(KC) = DFSW2(KC)
            IWDMG(KC)= IWDMG2(KC)
            FSWT=FSWT+FSW(KC)           !MINGHAO
            DO I=1,NEQSW
               X0SW(KC,I)= X0SW2(KC,I)
            END DO
            DO IX = 1, NELSW
               DO IY = 1, NGXW
                  D0PSW(KC,IX,IY)=D0PSW2(KC,IX,IY)
                  D0NSW(KC,IX,IY)=D0NSW2(KC,IX,IY)
                  P0PSW(KC,IX,IY)=P0PSW2(KC,IX,IY)
                  P0NSW(KC,IX,IY)=P0NSW2(KC,IX,IY)
                  W0PSW(KC,IX,IY)=W0PSW2(KC,IX,IY)
                  W0NSW(KC,IX,IY)=W0NSW2(KC,IX,IY)
               END DO
            END DO
            DO IX=1,NELSW
               DO I=1,NGXW
                  DO J=1,NGYW
                     S0SW(KC,IX,I,J)= S0SW2(KC,IX,I,J)
                     EPS0SW(KC,IX,I,J)= EPS0SW2(KC,IX,I,J)
                  END DO
               END DO
            END DO
         END DO

C
C
         IF (NDYN.EQ.0) THEN

            IF (FACT.NE.FACTOR(NS)) THEN    ! If intermediate convergence (STATIC case):
               NDET=0
               NDETT = NDETT + 1
               NITER=0
               GO TO 465
            END IF
            IF (FACT.EQ.FACTOR(NS)) GO TO 460  !Go to 460 if final result has been achieved for the step

         ENDIF
C

         IF (NDYN.EQ.1) THEN
            IF ((AXX.NE.ADX(NS)).OR.(AYY.NE.ADY(NS)))  THEN    !INTERMEDIATE CONVERGENCE ACHIEVED
               DT1 = DT1 - DT
               DT = DT1
               TCONV = TIM(NS) - DT
               NDET=0
               NDETT = NDETT + 1
               NITER=0
               GO TO 465
            END IF
         ENDIF
C
C * When final convergence within the step has been achieved:
C
C
 460     WRITE (*,1037) NS,NDETT, NITER, KC1,KC2,KC3
1037     FORMAT(/' **** FINISHED STEP =',I5,' :'/
     1   ' CONVERGED AFTER',I3,' INTERMEDIATE CONVERGENCES'/
     1   ' ITERATIONS DURING LAST INTERMEDIATE CONVERGENCE=',I3/
	1   ' KCF=',I2,' KCX=',I2,' KCU=',I2//)

C
C * Calls OUTPUT
C
C * RE is the vector used to compute the nodal applied loads, including
C   the reactions at the supports
C
 465     IF (NDYN.EQ.0) THEN
            IF (FACT.NE.FACTOR(NS)) THEN
               FACT0=FACT
               FACT=FACTOR(NS)
               CALL OUT (NDYN,NEO,IEO,KEO,IEP,NS,TCONV,NSUPP)
            END IF
            IF (FACT.EQ.FACTOR(NS)) THEN
               CALL OUT (NDYN,NEO,IEO,KEO,IEP,NS,TCONV,NSUPP)
               FACT0=FACT
            END IF
         END IF

         IF (NDYN.EQ.1) THEN
            IF ((AXX.NE.ADX(NS)).OR.(AYY.NE.ADY(NS))) THEN
C            CALL OUT (NDYN,NEO,IEO,KEO,IEP,NS,TCONV,NSUPP)
                  AXX0=AXX
                  AXX=ADX(NS)
                  AYY0=AYY
                  AYY=ADY(NS)
               END IF
               IF ((AXX.EQ.ADX(NS)).AND.(AYY.EQ.ADY(NS))) THEN
                  TCONV = TIM(NS)
                  CALL OUT (NDYN,NEO,IEO,KEO,IEP,NS,TCONV,NSUPP)
               END IF
            END IF

C
C   **** CHECK EACH SHEAR WALL DRIFT AND LOAD AT EACH CONVERGED LOAD STEP
C  

            WRITE(5,480) TIM(NS),
     1      (DSPSW(KC),FSW(KC),IWDMG(KC),KC=1,NWALL)
480         FORMAT(F10.3,100(2F12.4,I5))



5000  CONTINUE      ! Returns for the next step

      CALL CPU_TIME(time_end)
      WRITE (3,110) (time_end-time_begin)/60.0/4.0    !FOR QUAD-CORE CPU
110   FORMAT (//'TIME OF CPU OPERATION WAS ', F9.2, '  MINUTES.')

      CLOSE (2)
6000  CLOSE (3)
      CLOSE (5)
      CLOSE (8)
      CLOSE (9)

      END

C
C
C  ****** END OF MAIN PROGRAM ********
C
C




C
C
C
C
C--------------------------------------------------------------------
C
      SUBROUTINE MASSNOD (NOD, XMA, AX, AY)
C     (Checked OK!)
C
C     *CONTRIBUTIONS FROM A MASS XMA AT NODE NOD TO THE
C      RIGHT HAND SIDE, FOR THE GLOBAL ACCELERATION
C      VECTOR AX,AY
C
      USE COMMON_XV
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION ACC(2)
C
      ACC(1)=AX
      ACC(2)=AY
      DO I=1,2
         II=6*(NOD-1)+I
         RA(II)=RA(II)+XMA*ACC(I)
      END DO
      RETURN
      END
C-----------------------------------------------------------------
C
      SUBROUTINE OUT (NDYN,NEO,IEO,KEO,IEP,NS,T,NSUPP)
C     (Checked OK!)
C
C *  Output Subroutine ,information is written into unit #3,
C    plotting files at a specific node are written into unit #8.
C
      USE COMMON_XV
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      COMMON/SUPPORTS/NB(1000),MNB(1000),IB(1000,6),PL(1000,6)
C
      COMMON/SUPL/ NSUPL, ISUPL(1000)
C
      COMMON/CONV/ PSI1, PSI0, PSU0, PSU, PSI20, PSI2, FSWT
C
      COMMON/PLOT/ NEP,NJOTOT,NODMASS(1000),FNODMAS(1000)
      COMMON/GACC/ AXX,AYY,GA
C
      DIMENSION IEO(100), KEO(3)
      DIMENSION REA(6),IEP(1000),UX(1000),UY(1000),UZ(1000),
     1ACX(1000),ACY(1000)
C
      IF (NDYN.EQ.0) WRITE (3,10) NS
 10   FORMAT(/' STEP=',I4)
      IF (NDYN.EQ.1) WRITE (3,20) NS,T
 20   FORMAT(/' STEP=',I4,'  TIME=',F10.4)
      WRITE(3,25) PSI0, PSI1,PSI20,PSI2,PSU0,PSU
 25   FORMAT(/' CONVERGENCE DATA:'/
     1' PSI(initial) =',E11.3,' PSI(converged) =',E11.3/
	1' PSI2(initial)=',E11.3,' PSI2(converged)=',E11.3/
	1' PSU(initial) =',E11.3,' PSU(converged) =',E11.3)
      WRITE (3,30)
 30   FORMAT(/' NODE',6X,' DX',12X,' DY',12X,' DZ'/)
      DO I=1,NEO
         IJ=(IEO(I)-1)*6
         WX = KEO(1)*X(IJ+1)
         WY = KEO(2)*X(IJ+2)
         WZ = KEO(3)*X(IJ+3)
         IF (DABS(WX).LT.1.0E-10) WX = 0.0D0
         IF (DABS(WY).LT.1.0E-10) WY = 0.0D0
         IF (DABS(WZ).LT.1.0E-10) WZ = 0.0D0
         WRITE (3,35) IEO(I),WX,WY,WZ
 35      FORMAT(I4,3E12.4)
      END DO
      WRITE (3,40)
 40   FORMAT(/' SUPPORT'/
     1'  NODE',3X,' RX',8X,' RY',8X,' RZ',8X/)
C
C **  REA(I) are the reactions at support node I
C
C  * RXT, RYT, RZT are the summation of support force reactions,
C    in the X, Y and Z directions.
C  * RXT, RYT, RZT are given in absolute value.
C
      RXT = 0.0D0
      RYT = 0.0D0
      RZT = 0.0D0
      DO II = 1,NSUPP
         DO I=1,6
            IBC=(NB(II)-1)*6+I
            REA(I)=RC2(IBC)
            IF (DABS(REA(I)).LT.1.0E-10) REA(I) = 0.0D0
         END DO
         DO I = 1, NSUPL
            IF (ISUPL(I).EQ.NB(II)) THEN
               RXT = RXT + REA(1)
               RYT = RYT + REA(2)
               RZT = RZT + REA(3)
               WRITE (3,50) NB(II), REA(1),REA(2),REA(3)
 50            FORMAT(I4,3E12.4)
            END IF
         END DO
      END DO
      IF (DABS(RXT).LT.1.0E-10) RXT = 0.0D0
      IF (DABS(RYT).LT.1.0E-10) RYT = 0.0D0
      IF (DABS(RZT).LT.1.0E-10) RZT = 0.0D0
      DO I= 1,NEP
         IJ=(IEP(I)-1)*6
         UX(I) = KEO(1)*X(IJ+1)
         UY(I) = KEO(2)*X(IJ+2)
         IF (DABS(UX(I)).LT.1.0E-6) UX(I) = 0.0D0
         IF (DABS(UY(I)).LT.1.0E-6) UY(I) = 0.0D0
         ACX(I) = KEO(1)*ACP(IJ+1)
         ACY(I) = KEO(2)*ACP(IJ+2)
         IF (DABS(ACX(I)).LT.1.0E-6) ACX(I) = 0.0D0
         IF (DABS(ACY(I)).LT.1.0E-6) ACY(I) = 0.0D0
      END DO

      IF (NDYN.EQ.0) THEN
         WRITE(7,60) NS,RXT,RYT,RZT,(UX(I),UY(I),I=1,NEP)
60       FORMAT(I5,3F12.4,500(2F10.3))
      END IF


      IF (NDYN.EQ.1) THEN
C	   WRITE(7,70) T, RXT,RYT,RZT,FSWT,(UX(I),UY(I),I=1,NEP)  !DISP RESPONSE
            WRITE(7,70) T, RXT,RYT,RZT,(UX(I),UY(I),I=1,NEP)  !DISP RESPONSE
70          FORMAT(F8.3,3F12.4,500(2F10.3))

            WRITE(9,75)T,RXT,RYT,RZT,
     1      (ACX(I)/GA+AXX/GA,ACY(I)/GA+AYY/GA,I=1,NEP) !ACC RESPONSE
75          FORMAT(F8.3,3F12.4,500(2E15.4))

         END IF

         WRITE (2,80) NS, (UX(I),UY(I),UZ(I),I=1,NJOTOT)    !WRITING DEFORMATION DATA TO ResPlot file.
80       FORMAT (I4, 1000(3F10.2))
         RETURN
      END


C-----------------------------------------------------------------
C
      SUBROUTINE DECOMP (N, LHB, IERROR)
C     (Checked OK!)
C
C  * NEXT SUBROUTINES (DECOMP AND SOLV) ARE FOR SOLUTION OF SYSTEM 
C     OF EQUATIONS,USING CHOLESKY FOR BANDED,SYMMETRIC,
C     POSITIVE DEFINITE MATRIX*
C
      USE COMMON_A
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
C    * A IS STORED COLUMNWISE*
C
      IERROR=0
      KB=LHB-1
      IF (A(1).LE.0.0D0) THEN
         IERROR=1
         RETURN
      END IF
      XTMP=A(1)
      A(1)=DSQRT(XTMP)
      IF (N.EQ.1) RETURN
      DO 10 I=2,LHB
 10   A(I)=A(I)/A(1)
      DO 60 J=2,N
         J1=J-1
         IJD=(2*N-J)*(J-1)/2+J
         SUM=A(IJD)
         KO=1
         IF (J.GT.LHB) KO=J-KB
         DO 20 K=KO,J1
            JK=(2*N-K)*(K-1)/2+J
 20      SUM=SUM-A(JK)*A(JK)
         IF (SUM.LE.0.0D0) THEN
            IERROR=1
            RETURN
         END IF
         A(IJD)=DSQRT(SUM)
         DO 50 I=1,KB
            II=J+I
            IF (II.GT.N) GO TO 60
            KO=1
            IF (II.GT.LHB) KO=II-KB
            SUM=A(IJD+I)
            IF (I.EQ.KB) GO TO 40
            DO 30 K=KO,J1
               JK=(2*N-K)*(K-1)/2+J
               IK=(2*N-K)*(K-1)/2+II
 30         SUM=SUM-A(IK)*A(JK)
 40         A(IJD+I)=SUM/A(IJD)
 50      CONTINUE
 60   CONTINUE
      RETURN
      END
C
C  ---------------------------------------------------------------------
C
      SUBROUTINE SOLV (N, LHB)
C     (Checked OK!)
C
      USE COMMON_A
      USE COMMON_XV
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
C    * FORWARD SUBSTITUTION*
C
      KB=LHB-1
      PSI(1)=PSI(1)/A(1)
      IF (N.EQ.1) GO TO 30
      DO 20 I=2,N
         I1=I-1
         KO=1
         IF (I.GT.LHB) KO=I-KB
         SUM=PSI(I)
         II=(2*N-I)*(I-1)/2+I
         DO 10 K=KO,I1
            IK=(2*N-K)*(K-1)/2+I
 10      SUM=SUM-A(IK)*PSI(K)
         PSI(I)=SUM/A(II)
 20   CONTINUE
C
C    * BACKWARD SUBSTITUTION*
C
 30   N1=N-1
      LB=N*(N+1)/2
      PSI(N)=PSI(N)/A(LB)
      IF (N.EQ.1) RETURN
      DO 50 I=1,N1
         I1=N-I+1
         NI=N-I
         KO=N
         IF (I.GT.KB) KO=NI+KB
         SUM=PSI(NI)
         II=(2*N-NI)*(NI-1)/2+NI
         DO 40 K=I1,KO
            IK=(2*N-NI)*(NI-1)/2+K
 40      SUM=SUM-A(IK)*PSI(K)
         PSI(NI)=SUM/A(II)
 50   CONTINUE
      RETURN
      END

C---------------------------------------------------------------
C
      SUBROUTINE GLOBMF (NE, NEQ, CF)
C     (Checked OK!)
C
C  * Subroutine to store the frame member contribution
C    into the global matrix.
C
      USE COMMON_FRAME
      USE COMMON_A

      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION CF(12,12)
C
      DO 30 I=1,2
         I1=6*(I-1)
         II=6*(IEF(NE,I)-1)
         DO 20 J=1,2
            IF (IEF(NE,J).GT.IEF(NE,I)) GO TO 20
            J1=6*(J-1)
            JJ=6*(IEF(NE,J)-1)
            NJ=6
            DO 10 K=1,6
               IF (I.EQ.J) NJ=K
               IL=I1+K
               IIE=II+K
               DO 10 L=1,NJ
                  JE=J1+L
                  JJE=JJ+L
                  IJ=((2*NEQ-JJE)*(JJE-1))/2+IIE
                  A(IJ)=A(IJ)+CF(IL,JE)
 10         CONTINUE
 20      CONTINUE
 30   CONTINUE
      RETURN
      END

C---------------------------------------------------------------
C
      SUBROUTINE GLOBSW (II1,II2,NEQ,CF)
C     (Checked OK!)
C
C  * Subroutine to store the matrix for pseudo nail shear wall
C    into the global matrix. 
C
      USE COMMON_A
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION CF(12,12), IEC(2)

      IEC(1) = II1
      IEC(2) = II2
      DO 30 I=1,2
         I1=6*(I-1)
         II=6*(IEC(I)-1)
         DO 20 J=1,2
            IF (IEC(J).GT.IEC(I)) GO TO 20
            J1=6*(J-1)
            JJ=6*(IEC(J)-1)
            NJ=6
            DO 10 K=1,6
               IF (I.EQ.J) NJ=K
               IL=I1+K
               IIE=II+K
               DO 10 L=1,NJ
                  JE=J1+L
                  JJE=JJ+L
                  IJ=((2*NEQ-JJE)*(JJE-1))/2+IIE
                  A(IJ)=A(IJ)+CF(IL,JE)
 10         CONTINUE
 20      CONTINUE
 30   CONTINUE
      RETURN
      END
C---------------------------------------------------------------
C
      SUBROUTINE GLOBRF (NE, PSIE)
C     (Checked OK!)
C
C  * Subroutine to store the local PSIE vector for a frame member
C    into the global right hand side PSI.
C
      USE COMMON_FRAME
      USE COMMON_A
      USE COMMON_XV
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PSIE(12)

      DO 20 I=1,2
         I1=6*(I-1)
         II=6*(IEF(NE,I)-1)
         DO 10 K=1,6
            IL=I1+K
            IIE=II+K
            PSI(IIE)=PSI(IIE)+PSIE(IL)
            RC2(IIE)=RC2(IIE)+PSIE(IL)
 10      CONTINUE
 20   CONTINUE
      RETURN
      END
C
C---------------------------------------------------------------
C
      SUBROUTINE GMULT (NEQ, FDYN1, FDYN2,DT,KCASE)
C     (Checked OK!)
C
C  * Multiplication subroutine
C    KCASE = 1 : to multiply the matrix A0 times the current disp vector X, (Kt(BAR)-Ki)*X
C    KCASE = 2 : to multiply the matrix [M] and [C] with converged vector XP, VP AND AP.
C                [M]*(AP+4/DT*VP+4/DT/DT*XP)+[C]*(AP+4/DT*VP+4/DT/DT*XP)  
C                See p622 Table 15.3.2 of [Dynamics of Structures, 2nd Ed., A. K. Chopra]
C
      USE COMMON_A
      USE COMMON_XV
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DO I=1,NEQ
         IF (KCASE.EQ.1) RMU(I)=0.0D0
         IF (KCASE.EQ.2) RMUP(I)=0.0D0
         DO J=1,I
            IJ=((2*NEQ-J)*(J-1))/2+I
            IF (KCASE.EQ.1) THEN
               RMU(I) = RMU(I) + A0(IJ)*X(J)
            END IF
            IF (KCASE.EQ.2) THEN
               RMUP(I)=RMUP(I)+
     1         AMASS(IJ)*(4.0/DT/DT*XP(J)+4.0/DT*VP(J)+ACP(J))
     1         +AC(IJ)*(2.0/DT*XP(J)+VP(J))
            END IF
         END DO
         DO J=I+1,NEQ
            IJ=((2*NEQ-I)*(I-1))/2+J
            IF (KCASE.EQ.1) THEN
               RMU(I)=RMU(I)+A0(IJ)*X(J)
            END IF
            IF (KCASE.EQ.2) THEN
               RMUP(I)=RMUP(I)+
     1         AMASS(IJ)*(4.0/DT/DT*XP(J)+4.0/DT*VP(J)+ACP(J))
     1         +AC(IJ)*(2.0/DT*XP(J)+VP(J))
            END IF
         END DO
      END DO
      RETURN
      END
C-----------------------------------------------------------------
C
      SUBROUTINE BC1 (IBC, BCONF, NEQ, NPASS)
C     (Checked OK!)
C
C  * Subroutine to introduce boundary condition BCONF in degree of
C    freedom (row) IBC. This subroutine modifies the
C    right hand side PSI. NPASS=1 modifies the right hand side
C    except for the IBC row. NPASS = 2 modifies the right hand
C    side at row IBC.
C
      USE COMMON_A
      USE COMMON_XV

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

C
      IF (NPASS.EQ.2) GO TO 50
      IF (IBC.EQ.1) GO TO 20
      DO 10 J=1,IBC-1
         IJ=((2*NEQ-J)*(J-1))/2+IBC
         PSI(J)=PSI(J)-A(IJ)*(BCONF-X(IBC))
 10   CONTINUE
 20   IF (IBC.EQ.NEQ) GO TO 35
      DO 30 J=IBC+1,NEQ
         IJ=((2*NEQ-IBC)*(IBC-1))/2+J
         PSI(J)=PSI(J)-A(IJ)*(BCONF-X(IBC))
 30   CONTINUE
 35   RETURN
 50   IJ=((2*NEQ-IBC)*(IBC-1))/2+IBC
      PSI(IBC)=BCONF-X(IBC)
      RETURN
      END
C-----------------------------------------------------------------
C
      SUBROUTINE BC2 (IBC, NEQ)
C     (Checked OK!)
C
C  * Subroutine to introduce boundary condition BCONF in degree of
C    freedom (row) IBC. This subroutine modifies the matrix A.
C
      USE COMMON_A
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

C
      IF (IBC.EQ.1) GO TO 20
      DO 10 J=1,IBC-1
         IJ=((2*NEQ-J)*(J-1))/2+IBC
         A(IJ) = 0.0D0
 10   CONTINUE
 20   IF (IBC.EQ.NEQ) GO TO 35
      DO 30 J=IBC+1,NEQ
         IJ=((2*NEQ-IBC)*(IBC-1))/2+J
         A(IJ) = 0.0D0
 30   CONTINUE
 35   IJ=((2*NEQ-IBC)*(IBC-1))/2+IBC
      A(IJ) = 1.0D0
      RETURN
      END
C
C    -------------------------------------------------------------
C
      SUBROUTINE SHAPESF(JF,NG)
C
C  ** Obtains shape functions for the frame member elements **
C
      USE COMMON_FRAME

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      REAL*8 M0F, M1F, M2F, N0F, N1F, L0F, L1F, L2F, K0F, K1F
C
      COMMON/GAUSSI/                                         !Gauss quadrature
     1EGXF(16),HGXF(16),EGYF(16),HGYF(16),EGZF(16),HGZF(16),
     1EGXW(16),HGXW(16),EGYW(16),HGYW(16),
     1NGXF,NGYF,NGZF,NGXW,NGYW
C
      COMMON/FSHAPE/                                         !Frame shapes
     1N0F(12,5),N1F(12,5),M0F(12,5),M1F(12,5),M2F(12,5),
     1L0F(12,5),L1F(12,5),L2F(12,5),K0F(12,5),K1F(12,5)
C
      DELTA = DELTAF(JF)
      DO 10 I=1,NG
      X=(1.0 + EGXF(I))*DELTA/2.0D0
      X1 = X/DELTA
      X2 = (X/DELTA)**2
      X3 = (X/DELTA)**3
C
C     Shapes for W(X)  (Bending displacement, rotation around z-axis)
C
      M0F(1,I)= 0.0D0
      M0F(2,I)= 1.0 - 3.0D0*X2 + 2.0D0*X3
      M0F(3,I)= X*(1.0 - 2.0D0*X1 + X2)
      M0F(4,I)= 0.0D0
      M0F(5,I)= 0.0D0
      M0F(6,I)= 0.0D0
      M0F(7,I)= 0.0D0
      M0F(8,I)= 3.0D0*X2 - 2.0D0*X3
      M0F(9,I)= X*(X2-X1)
      M0F(10,I)=0.0D0
      M0F(11,I)=0.0D0
      M0F(12,I)=0.0D0
C
C    Shapes for W'(X)
C
      M1F(1,I)= 0.0D0
      M1F(2,I)= -6.0D0*X1/DELTA + 6.0D0*X2/DELTA
      M1F(3,I)= 1.0 - 4.0*X1 + 3.0D0*X2
      M1F(4,I)= 0.0D0
      M1F(5,I)= 0.0D0
      M1F(6,I)= 0.0D0
      M1F(7,I)= 0.0D0
      M1F(8,I)= 6.0D0*X1/DELTA - 6.0D0*X2/DELTA
      M1F(9,I)= 3.0D0*X2 - 2.0D0*X1
      M1F(10,I)= 0.0D0
      M1F(11,I)=0.0D0
      M1F(12,I)=0.0D0
C
C     Shapes for W"(X)
C
      M2F(1,I)= 0.0D0
      M2F(2,I)= -6.0D0/(DELTA**2) + 12.0D0*X1/(DELTA**2)
      M2F(3,I)= -4.0/DELTA + 6.0D0*X1/DELTA
      M2F(4,I)= 0.0D0
      M2F(5,I)= 0.0D0
      M2F(6,I)= 0.0D0
      M2F(7,I)= 0.0D0
      M2F(8,I)= 6.0D0/(DELTA**2) - 12.0D0*X1/(DELTA**2)
      M2F(9,I)= 6.0D0*X1/DELTA - 2.0D0/DELTA
      M2F(10,I)= 0.0D0
      M2F(11,I)=0.0D0
      M2F(12,I)=0.0D0
C
C     Shapes for V(X)  (Bending displacement, rotation around y-axis)
C
      L0F(1,I)= 0.0D0
      L0F(2,I)= 0.0D0
      L0F(3,I)= 0.0D0
      L0F(4,I)= 1.0 - 3.0D0*X2 + 2.0D0*X3
      L0F(5,I)= X*(1.0 - 2.0D0*X1 + X2)
      L0F(6,I)= 0.0D0
      L0F(7,I)= 0.0D0
      L0F(8,I)=0.0D0
      L0F(9,I)=0.0D0
      L0F(10,I)= 3.0D0*X2 - 2.0D0*X3
      L0F(11,I)= X*(X2-X1)
      L0F(12,I)= 0.0D0
C
C    Shapes for V'(X)
C
      L1F(1,I)= 0.0D0
      L1F(2,I)= 0.0D0
      L1F(3,I)= 0.0D0
      L1F(4,I)= -6.0D0*X1/DELTA + 6.0D0*X2/DELTA
      L1F(5,I)= 1.0 - 4.0*X1 + 3.0D0*X2
      L1F(6,I)= 0.0D0
      L1F(7,I)= 0.0D0
      L1F(8,I)= 0.0D0
      L1F(9,I)=0.0D0
      L1F(10,I)= 6.0D0*X1/DELTA - 6.0D0*X2/DELTA
      L1F(11,I)= 3.0D0*X2 - 2.0D0*X1
      L1F(12,I)= 0.0D0
C
C     Shapes for V"(X)
C
      L2F(1,I)= 0.0D0
      L2F(2,I)= 0.0D0
      L2F(3,I)= 0.0D0
      L2F(4,I)= -6.0D0/(DELTA**2) + 12.0D0*X1/(DELTA**2)
      L2F(5,I)= -4.0/DELTA + 6.0D0*X1/DELTA
      L2F(6,I)= 0.0D0
      L2F(7,I)= 0.0D0
      L2F(8,I)= 0.0D0
      L2F(9,I)=0.0D0
      L2F(10,I)= 6.0D0/(DELTA**2) - 12.0D0*X1/(DELTA**2)
      L2F(11,I)= 6.0D0*X1/DELTA - 2.0D0/DELTA
      L2F(12,I)= 0.0D0
C
C     Shapes for U(X)   (Axial displacement)
C
      N0F(1,I)= 1.0 - X1
      N0F(2,I)= 0.0D0
      N0F(3,I)= 0.0D0
      N0F(4,I)= 0.0D0
      N0F(5,I)= 0.0D0
      N0F(6,I)= 0.0D0
      N0F(7,I)= X1
      N0F(8,I)= 0.0D0
      N0F(9,I)= 0.0D0
      N0F(10,I)= 0.0D0
      N0F(11,I)= 0.0D0
      N0F(12,I)= 0.0D0
C
C     SHAPES FOR U'(X)
C
      N1F(1,I)= -1.0D0/DELTA
      N1F(2,I)= 0.0D0
      N1F(3,I)= 0.0D0
      N1F(4,I)= 0.0D0
      N1F(5,I)= 0.0D0
      N1F(6,I)= 0.0D0
      N1F(7,I)= 1.0D0/DELTA
      N1F(8,I)= 0.0D0
      N1F(9,I)= 0.0D0
      N1F(10,I)= 0.0D0
      N1F(11,I)= 0.0D0
      N1F(12,I)= 0.0D0
C
C
C     Shapes for PHI(X) (Torsional rotation)
C
      K0F(1,I)= 0.0D0
      K0F(2,I)= 0.0D0
      K0F(3,I)= 0.0D0
      K0F(4,I)= 0.0D0
      K0F(5,I)= 0.0D0
      K0F(6,I)= 1.0D0 - X1
      K0F(7,I)= 0.0D0
      K0F(8,I)= 0.0D0
      K0F(9,I)= 0.0D0
      K0F(10,I)= 0.0D0
      K0F(11,I)= 0.0D0
      K0F(12,I)= X1

C
C     SHAPES FOR PHI'(X)
C
      K1F(1,I)= 0.0D0
      K1F(2,I)= 0.0D0
      K1F(3,I)= 0.0D0
      K1F(4,I)= 0.0D0
      K1F(5,I)= 0.0D0
      K1F(6,I)= -1.0D0/DELTA
      K1F(7,I)= 0.0D0
      K1F(8,I)= 0.0D0
      K1F(9,I)= 0.0D0
      K1F(10,I)= 0.0D0
      K1F(11,I)= 0.0D0
      K1F(12,I)= 1.0/DELTA
C
 10   CONTINUE
      RETURN
      END


C   
C    ---------------------------------------------------------
C
      SUBROUTINE SOLVS (N, LHB, A, B)
C     (Checked OK!)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C     REAL*8 A, B, TEMP, SUM
      DIMENSION A(2688), B(168)
C
C     * FORWARD SUBSTITUTION *
C
      KB=LHB-1
      TEMP=A(1)
      B(1)=B(1)/TEMP
      DO 20 I=2,N
         I1=I-1
         KO=1
         IF (I.GT.LHB) KO=I-KB
         SUM=B(I)
         II=LHB*I-KB
         DO 10 K=KO,I1
            IK=KB*K+I-KB
            TEMP=A(IK)
            SUM=SUM-TEMP*B(K)
 10      CONTINUE
         B(I)=SUM/A(II)
 20   CONTINUE
C
C     * BACKWARD SUBSTITUTION *
C
      N1=N-1
      LB=LHB*N-KB
      TEMP=A(LB)
      B(N)=B(N)/TEMP
      DO 40 I=1,N1
         I1=N-I+1
         NI=N-I
         KO=N
         IF (I.GT.KB) KO=NI+KB
         SUM=B(NI)
         II=LHB*NI-KB
         DO 30 K=I1,KO
            IK=KB*NI+K-KB
            TEMP=A(IK)
            SUM=SUM-TEMP*B(K)
 30      CONTINUE
         B(NI)=SUM/A(II)
 40   CONTINUE
      RETURN
      END
C------------------------------------------------------------
      SUBROUTINE GAUSS (N, E, H, IERR)
C      (Checked OK!)
C
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION E(16), H(16)

      M=(N-2)*(N-3)*(N-4)*(N-5)*(N-6)*(N-7)*(N-8)
      M=M*(N-9)*(N-10)*(N-11)*(N-12)*(N-15)*(N-16)
      IF (M.NE.0) GO TO 260
      IERR=0
      IF (N.EQ.16) GO TO 220
      IF (N.EQ.15) GO TO 200
      IF (N.EQ.12) GO TO 180
      IF (N.EQ.11) GO TO 160
      IF (N.EQ.10) GO TO 140
      IF (N.EQ.9) GO TO 120
      IF (N.EQ.8) GO TO 100
      IF (N.EQ.7) GO TO 80
      IF (N.EQ.6) GO TO 60
      IF (N.EQ.5) GO TO 40
      IF (N.EQ.4) GO TO 20
      IF (N.EQ.3) GO TO 10
      E(1)=0.577350269189626D0
      E(2)=-E(1)
      H(1)=1.0
      H(2)=H(1)
      RETURN
 10   E(1)=0.774596669241483D0
      E(2)=0.0D0
      E(3)=-E(1)
      H(1)=0.555555555555556D0
      H(2)=0.888888888888889D0
      H(3)=H(1)
      RETURN
 20   E(1)=0.861136311594053D0
      E(2)=0.339981043584856D0
      H(1)=0.347854845137454D0
      H(2)=0.652145154862546D0
      DO 30 I=1,2
         E(5-I)=-E(I)
 30   H(5-I)=H(I)
      RETURN
 40   E(1)=0.906179845938664D0
      E(2)=0.538469310105683D0
      E(3)=0.0D0
      H(1)=0.236926885056189D0
      H(2)=0.478628670499366D0
      H(3)=0.568888888888889D0
      DO 50 I=1,2
         E(6-I)=-E(I)
 50   H(6-I)=H(I)
      RETURN
 60   E(1)=0.932469514203152D0
      E(2)=0.661209386466265D0
      E(3)=0.238619186083197D0
      H(1)=0.171324492379170D0
      H(2)=0.360761573048139D0
      H(3)=0.467913934572691D0
      DO 70 I=1,3
         E(7-I)=-E(I)
 70   H(7-I)=H(I)
      RETURN
 80   E(1)=0.949107912342759D0
      E(2)=0.741531185599394D0
      E(3)=0.405845151377397D0
      E(4)=0.0D0
      H(1)=0.129484966168870D0
      H(2)=0.279705391489277D0
      H(3)=0.381830050505119D0
      H(4)=0.417959183673469D0
      DO 90 I=1,3
         E(8-I)=-E(I)
 90   H(8-I)=H(I)
      RETURN
 100  E(1)=0.960289856497536D0
      E(2)=0.796666477413627D0
      E(3)=0.525532409916329D0
      E(4)=0.183434642495650D0
      H(1)=0.101228536290376D0
      H(2)=0.222381034453374D0
      H(3)=0.313706645877887D0
      H(4)=0.362683783378362D0
      DO 110 I=1,4
         E(9-I)=-E(I)
 110  H(9-I)=H(I)
      RETURN
 120  E(1)=0.968160239507626D0
      E(2)=0.836031107326636D0
      E(3)=0.613371432700590D0
      E(4)=0.324253423403809D0
      E(5)=0.0D0
      H(1)=0.081274388361574D0
      H(2)=0.180648160694857D0
      H(3)=0.260610696402935D0
      H(4)=0.312347077040003D0
      H(5)=0.330239355001260D0
      DO 130 I=1,4
         E(10-I)=-E(I)
 130  H(10-I)=H(I)
      RETURN
 140  E(1)=0.973906528517172D0
      E(2)=0.865063366688985D0
      E(3)=0.679409568299024D0
      E(4)=0.433395394129247D0
      E(5)=0.148874338981631D0
      H(1)=0.066671344308688D0
      H(2)=0.149451349150581D0
      H(3)=0.219086362515982D0
      H(4)=0.269266719309996D0
      H(5)=0.295524224714753D0
      DO 150 I=1,5
         E(11-I)=-E(I)
 150  H(11-I)=H(I)
      RETURN
 160  E(1)=0.978228658146057D0
      E(2)=0.887062599768095D0
      E(3)=0.730152005574049D0
      E(4)=0.519096129206812D0
      E(5)=0.269543155952345D0
      E(6)=0.0D0
      H(1)=0.055668567116174D0
      H(2)=0.125580369464905D0
      H(3)=0.186290210927734D0
      H(4)=0.233193764591990D0
      H(5)=0.262804544510247D0
      H(6)=0.272925086777901D0
      DO 170 I=1,5
         E(12-I)=-E(I)
 170  H(12-I)=H(I)
      RETURN
 180  E(1)=0.981560634246719D0
      E(2)=0.904117256370475D0
      E(3)=0.769902674194305D0
      E(4)=0.587317954286617D0
      E(5)=0.367831498998180D0
      E(6)=0.125233408511469D0
      H(1)=0.047175336386512D0
      H(2)=0.106939325995318D0
      H(3)=0.160078328543346D0
      H(4)=0.203167426723066D0
      H(5)=0.233492536538355D0
      H(6)=0.249147045813403D0
      DO 190 I=1,6
         E(13-I)=-E(I)
 190  H(13-I)=H(I)
      RETURN
 200  E(1)=0.987992518020485D0
      E(2)=0.937273392400706D0
      E(3)=0.848206583410427D0
      E(4)=0.724417731360170D0
      E(5)=0.570972172608539D0
      E(6)=0.394151347077563D0
      E(7)=0.201194093997435D0
      E(8)=0.0D0
      H(1)=0.030753241996117D0
      H(2)=0.070366047488108D0
      H(3)=0.107159220467172D0
      H(4)=0.139570677926154D0
      H(5)=0.166269205816994D0
      H(6)=0.186161000015562D0
      H(7)=0.198431485327112D0
      H(8)=0.202578241925561D0
      DO 210 I=1,7
         E(16-I)=-E(I)
 210  H(16-I)=H(I)
      RETURN
 220  E(1)=0.989400934991650D0
      E(2)=0.944575023073233D0
      E(3)=0.865631202387832D0
      E(4)=0.755404408355003D0
      E(5)=0.617876244402644D0
      E(6)=0.458016777657227D0
      E(7)=0.281603550779259D0
      E(8)=0.095012509837637D0
      H(1)=0.027152459411754D0
      H(2)=0.062253523938648D0
      H(3)=0.095158511682493D0
      H(4)=0.124628971255534D0
      H(5)=0.149595988816577D0
      H(6)=0.169156519395003D0
      H(7)=0.182603415044924D0
      H(8)=0.189450610455068D0
      DO 230 I=1,8
         E(17-I)=-E(I)
 230  H(17-I)=H(I)
      RETURN
 260  WRITE (*,270)
 270  FORMAT(' WRONG CHOICE FOR GAUSS INTEGRATION POINTS IN SHYST'/)
      IERR=1
      RETURN
      END
C
C
C    *----------------------------------------------------------------------------
C
      SUBROUTINE SHYSTW(KC,IMAT,A00,A11,FY0,IDMG,IERRORS)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      REAL*8 M0, M1, M2, N0, N1

      COMMON/GAUSSI/                                         !Gauss quadrature
     1EGXF(16),HGXF(16),EGYF(16),HGYF(16),EGZF(16),HGZF(16),
     1EGXW(16),HGXW(16),EGYW(16),HGYW(16),
     1NGXF,NGYF,NGZF,NGXW,NGYW

      COMMON /WALL/                           !SHEAR WALLS
     1N0(100,2,6,9),N1(100,2,6,9),
     1M0(100,2,6,9),M1(100,2,6,9),M2(100,2,6,9),
     1X0SW(100,63),D0PSW(100,20,9),D0NSW(100,20,9),
     1P0PSW(100,20,9),P0NSW(100,20,9),
     1W0PSW(100,20,9),W0NSW(100,20,9),
	1S0SW(100,20,9,16),EPS0SW(100,20,9,16),
     1X0SW2(100,63),D0PSW2(100,20,9),D0NSW2(100,20,9),
     1P0PSW2(100,20,9),P0NSW2(100,20,9),
     1W0PSW2(100,20,9),W0NSW2(100,20,9),
	1S0SW2(100,20,9,16),EPS0SW2(100,20,9,16),
     1FSW0(100),DFSW0(100),
	1FSW2(100),DFSW2(100),DSPSW2(100),
     1FSW(100) ,DFSW(100) ,DSPSW(100),DELTA(100,2),
	1IESW(100,2),IWTYPE(100),IWORIEN(100),IWDMG(100),IWDMG2(100),
	1NELSW,NELSW_1

      COMMON /PNAIL/                            !Pseudo nail and embedment
     1XLSW(10),XDSW(10),ESW(10),SYSW(10),TTC(10),WDMGC(10),
     1Q0SWP(10),Q1SWP(10),Q2SWP(10),Q3SWP(10),Q4SWP(10),
     1XKSWP(10),DMAXSWP(10),PMAXSWP(10),SDFP(10),
     1Q0SWN(10),Q1SWN(10),Q2SWN(10),Q3SWN(10),Q4SWN(10),
     1XKSWN(10),DMAXSWN(10),PMAXSWN(10),SDFN(10)

      DIMENSION S0(20,9,16),EPS0(20,9,16),DSDE0(20,9,16),
     1X0(63),X0PR(63)
      DIMENSION XE(6),CE(6,6),PSIE(6)
      DIMENSION C(630),PSI(63),RM(63),RMPR(63)
      DIMENSION NBC(21),KBC(21),IBC(21,3)
      DIMENSION EPS3(20,9)
      DIMENSION D0P(20,9),D0N(20,9),P0(20,9),DPDW0(20,9),SW0(20,9)
      DIMENSION P0P(20,9),P0N(20,9),W0P(20,9),W0N(20,9)
      DIMENSION B(16)
      DIMENSION ITYP(10),ILAY(10),IDISP(50),IBEG(10),IFIN(10)

C
C    -----------------------------------------------------------------
C
      DF=XDSW(IMAT)      !Nail diameter  (mm)
      XLF=XLSW(IMAT)     !Nail length    (mm)
      EF=ESW(IMAT)       !Steel MOE      (kN/mm/mm)
      SYF=SYSW(IMAT)     !Yield strength (kN/mm/mm)
      TC=TTC(IMAT)           !Panel thickness(mm)
      IERRORS = 0
      NELEM=20        !Total number of element
      NELEM1=18      !Number of elements on nail pointy side
      NELEM2=NELEM-NELEM1 !Number of elements on nail head side
      NLAYER=2            ! NUMBER OF WOODEN LAYERS (TIMBER + OSB)
      ITYP(1)=1           ! WOOD MATERIAL TYPE
      ITYP(2)=2

      IBEG(1)=1
      IFIN(1)=NELEM1
      IBEG(2)=NELEM1+1
      IFIN(2)=NELEM


      NNBC=1               !B.C. CONDITIONS
      NBC(1)=1
      KBC(1)=1
      IBC(1,1)=3

C
C * Sets up indicators of moving layers. Always layer 1 (frame) moving relative to other layers
C  
      DO J=1,NELEM1
         IDISP(J)=0
      ENDDO

      DO J=NELEM1+1,NELEM
         IDISP(J)=1
      ENDDO

      DO I=1,NGYW
         B(I)=DF*DSQRT(1.0D0-EGYW(I)**2)
      END DO
C    ENTER TOLERANCES FOR HYSTERESIS SUBROUTINE SHYST:
C                 TOLHFA = ABSOLUTE CONVERGENCE FOR PSI VECTOR, OR
C                          OUT-OF-BALANCE VECTOR (e.g., 1.0E-06)
C                 TOLHXA = ABSOLUTE CONVERGENCE FOR VECTOR INCREMENT
C                          IN X (e.g., 1.0E-06)

C
C    ---------------------------------------------------------
C
C  *  MAIN PROGRAM, SIZE OF VECTORS
C
      NEQ=(NELEM+1)*3
      LHB=6
      NA=NEQ*LHB
C
C    * STARTS SOLUTION FOR THE STEP FROM A0 TO AA
      A0=A00
      AA=A11
      NDE=0
600   NITER=0
      KNITER=0
C
C    * SELECT A STARTING DELFECTED SHAPE, BETWEEN A0 AND AA.
C
      DO I=1,NEQ
         X0(I)=X0SW2(KC,I)
      END DO

      DO I=1,NELEM+1
         KAA=(I-1)*3+1
         X0(KAA)=X0(KAA)+(AA-A0)/2.0D0
      END DO


700   DO IX=1,NELEM
         DO IY=1,NGXW
            D0P(IX,IY)=D0PSW2(KC,IX,IY)
            D0N(IX,IY)=D0NSW2(KC,IX,IY)
            P0P(IX,IY)=P0PSW2(KC,IX,IY)
            P0N(IX,IY)=P0NSW2(KC,IX,IY)
            W0P(IX,IY)=W0PSW2(KC,IX,IY)
            W0N(IX,IY)=W0NSW2(KC,IX,IY)
         END DO
      END DO
      DO IX=1,NELEM
         DO IY=1,NGXW
            DO IZ=1,NGYW
               S0(IX,IY,IZ)  =S0SW2(KC,IX,IY,IZ)
               EPS0(IX,IY,IZ)=EPS0SW2(KC,IX,IY,IZ)
            END DO
         END DO
      END DO


C
C    * OBTAIN REACTION LOADS AND STRESSES FOR CURRENT VECTOR X{0}
C     
      Q00P=Q0SWP(IMAT)         !TIMBER LAYER
      Q10P=Q1SWP(IMAT)
      Q40P=Q4SWP(IMAT)
      DMAX0P=DMAXSWP(IMAT)
      XK0P=XKSWP(IMAT)
      PMAX0P=PMAXSWP(IMAT)
      SDF0P=SDFP(IMAT)

      Q00N=Q0SWN(IMAT)        !SHEATHING LAYER
      Q10N=Q1SWN(IMAT)
      Q40N=Q4SWN(IMAT)
      DMAX0N=DMAXSWN(IMAT)
      XK0N=XKSWN(IMAT)
      PMAX0N=PMAXSWN(IMAT)
      SDF0N=SDFN(IMAT)

      IF (IDMG.EQ.2) THEN
         SDF0P=0.0
         SDF0N=0.0
      ENDIF
      IF (IDMG.EQ.1) THEN
         SDF0P=1.0
         SDF0N=1.0
      ENDIF


      DO 800 JL=1,NLAYER

C     **MINGHAO
         IF (JL.EQ.1) THEN
            Q00=Q00P
            Q10=Q10P
            Q40=Q40P
            XK0=XK0P
            DMAX0=DMAX0P
            SDF0=SDF0P
            PMAX0=PMAX0P
         ENDIF
         IF (JL.EQ.2) THEN
            Q00=Q00N
            Q10=Q10N
            Q40=Q40N
            XK0=XK0N
            DMAX0=DMAX0N
            SDF0=SDF0N
            PMAX0=PMAX0N
         ENDIF
C    **MINGHAO

      DO 800 IE=IBEG(JL),IFIN(JL)
         DO  I=1,6
            JE=(IE-1)*3+I
            XE(I)=X0(JE)
         END DO
         DO 790 IX=1,NGXW
            W=0.0
            EPS1=0.0
            EPS2=0.0
            EPS3(IE,IX)=0.0
            DO K=1,6
               W=W+M0(KC,JL,K,IX)*XE(K)
               EPS1=EPS1+N1(KC,JL,K,IX)*XE(K)
               EPS2=EPS2+M2(KC,JL,K,IX)*XE(K)
            END DO
            IF(IDISP(IE).EQ.0) WT=W
            IF(IDISP(IE).EQ.1) WT=W-AA
            IF (WT.GE.0.0) D0=D0PSW2(KC,IE,IX)
            IF (WT.LT.0.0) D0=D0NSW2(KC,IE,IX)
            IF (WT.GE.0.0) SW0(IE,IX)=1.0
            IF (WT.LT.0.0) SW0(IE,IX)=-1.0
            IF (WT.GE.0.0) PP0=P0PSW2(KC,IE,IX)
            IF (WT.LT.0.0) PP0=P0NSW2(KC,IE,IX)
            IF (WT.GE.0.0) WW0=W0PSW2(KC,IE,IX)
            IF (WT.LT.0.0) WW0=W0NSW2(KC,IE,IX)
            AW=DABS(WT)
            CALL PSUP (D0,WW0,PP0,AW,Q00,Q10,Q40,XK0,
     1      DMAX0,SDF0,PMAX0,P,DPDW)
            IF (WT.GE.0.0) D0P(IE,IX)=D0
            IF (WT.LT.0.0) D0N(IE,IX)=D0
            IF (WT.GE.0.0) P0P(IE,IX)=P
            IF (WT.LT.0.0) P0N(IE,IX)=P
            IF (WT.GE.0.0) W0P(IE,IX)=AW
            IF (WT.LT.0.0) W0N(IE,IX)=AW
            P0(IE,IX)=P
            DPDW0(IE,IX)=DPDW
            DO IY=1,NGYW
               Y=EGYW(IY)*DF/2.0
               EPS=EPS1-Y*EPS2
               ST0=S0SW2(KC,IE,IX,IY)
               EP0=EPS0SW2(KC,IE,IX,IY)
               CALL STRESS (ST0, EP0, EPS, EF, SYF, S, DSDEPS)
               S0(IE,IX,IY)=S
               EPS0(IE,IX,IY)=EPS
               DSDE0(IE,IX,IY)=DSDEPS
            END DO
 790     CONTINUE
 800     CONTINUE

!C
!  *  FOR EACH ELEMENT, CONSTRUCT TANGENT STIFFNESS MATRIX AND
!     RIGHT HAND SIDE, AND ADD TO GLOBALS
!C
         DO  I=1,NA
            C(I)=0.0D0
         END DO
C
      DO 900 JL=1,NLAYER
         DO 900 IE=IBEG(JL),IFIN(JL)
            DO 860 I=1,6
               DO 850 J=1,I
                  CE(I,J)=0.0
                  DO 840 IX=1,NGXW
                     DPDW=DPDW0(IE,IX)
                     CE(I,J)=CE(I,J)+HGXW(IX)*DPDW*M0(KC,JL,I,IX)
     1               *M0(KC,JL,J,IX)*DELTA(KC,JL)/2.0
                     FACI=N1(KC,JL,I,IX)
                     FACJ=N1(KC,JL,J,IX)
                     DO 830 IY=1,NGYW
                        Y=EGYW(IY)*DF/2.0
                        DSDEPS=DSDE0(IE,IX,IY)
                        S=S0(IE,IX,IY)
                        FAC1=FACJ-Y*M2(KC,JL,J,IX)
                        FAC2=FACI-Y*M2(KC,JL,I,IX)
                        CE(I,J)=CE(I,J)+HGXW(IX)*HGYW(IY)*(DSDEPS*FAC1*FAC2
     1                  *B(IY))*(DELTA(KC,JL)*DF)/4.0
 830                 CONTINUE
 840              CONTINUE
 850           CONTINUE
 860        CONTINUE
            DO 880 I=1,6
               II=(IE-1)*3+I
               DO 870 J=1,I
                  JJ=(IE-1)*3+J
                  IJ=(JJ-1)*(LHB-1)+II
 870           C(IJ)=C(IJ)+CE(I,J)
 880        CONTINUE
 900  CONTINUE
C

      DO I=1,NEQ
         PSI(I)=0.0
      END DO
      DO 960 JL=1,NLAYER
         DO 960 IE=IBEG(JL),IFIN(JL)
            DO 940 I=1,6
               PSIE(I)=0.0D0
               DO 930 IX=1,NGXW
                  SW=SW0(IE,IX)
                  P=P0(IE,IX)
                  PSIE(I)=PSIE(I)+HGXW(IX)*P*SW*M0(KC,JL,I,IX)
     1            *DELTA(KC,JL)/2.0D0
                  FACI=N1(KC,JL,I,IX)
                  DO 920 IY=1,NGYW
                     Y=EGYW(IY)*DF/2.0D0
                     S=S0(IE,IX,IY)
                     FAC2=FACI-Y*M2(KC,JL,I,IX)
                     PSIE(I)=PSIE(I)+HGXW(IX)*HGYW(IY)*S*FAC2*B(IY)*
     1               (DELTA(KC,JL)*DF)/4.0D0
 920              CONTINUE
 930           CONTINUE
 940        CONTINUE
            DO 950 I=1,6
               II=(IE-1)*3+I
               PSI(II)=PSI(II)+PSIE(I)
 950        CONTINUE
 960  CONTINUE
C
      DO I=1,NEQ
         RM(I)=-PSI(I)
      END DO

C
C    * INTRODUCE B.C.
C
      DO 1090 IB=1,NNBC
         DO 1080 J=1,KBC(IB)
            K=(NBC(IB)-1)*3+IBC(IB,J)
            IF (K.EQ.1) GO TO 1050
            K1=K-LHB+1
            K2=K-1
            IF (K1.LE.0) K1=1
            DO 1040 I=K1,K2
               II=(I-1)*(LHB-1)+K
               RM(I)=RM(I)+C(II)*X0(K)
 1040       C(II)=0.0
 1050       IF (K.EQ.NEQ) GO TO 1070
            K1=K+1
            K2=K+LHB-1
            IF (K2.GT.NEQ) K2=NEQ
            DO 1060 I=K1,K2
               II=(K-1)*(LHB-1)+I
               RM(I)=RM(I)+C(II)*X0(K)
 1060       C(II)=0.0
 1070       KK=(K-1)*(LHB-1)+K
            C(KK)=1.0
            RM(K)=-X0(K)
 1080    CONTINUE
 1090 CONTINUE


      IF (KNITER.NE.0) GO TO 1133
      IF (KNITER.EQ.0) THEN
         DO I = 1, NEQ
            RMPR(I) = RM(I)
            X0PR(I) = X0(I)
         END DO
      END IF


1113  CALL DECOMPSW (NEQ, LHB, C, IERROR)
      IF (IERROR.EQ.1) GO TO 1170
C
C    * Use decomposed [C] and solve (solution stored in {RM}.)
C
      CALL SOLVSW (NEQ, LHB, C, RM)

C
C    * Update the vector {X0}
C      
      DO I=1,NEQ
         X0(I)=X0(I)+RM(I)
      ENDDO
      NITER = NITER + 1
      KNITER=1
      IF (NITER.EQ.50) GOTO 1170
      GOTO 700


C    * CHECK CONVERGENCE
C
C  *  Checks for convergence     
1133  KCW1 = 0
      KCW2 = 0
      PSI2 = 0.0D0
      DX2 = 0.0D0

      DO I=1,NEQ
         PSI2=PSI2+(RM(I)-RMPR(I))**2
      ENDDO
      PSI2 = DSQRT(PSI2)

      DO IE=1,NELEM
         JA=(IE-1)*3+1
         DX2=DX2+(X0PR(JA)-X0(JA))**2
      ENDDO
      DX2 = DSQRT(DX2)

      TOLHFA = 1.0D-06
      TOLHXA = 1.0D-06
      IF (PSI2.LE.TOLHFA) KCW1=1
      IF (DX2.LE.TOLHXA)  KCW2=1
      IF (KCW1.EQ.1.AND.KCW2.EQ.1) GO TO 1190
      DO I = 1, NEQ
         X0PR(I) = X0(I)
         RMPR(I) = RM(I)
      END DO
      GO TO 1113


1170  AA=(A0+AA)/2.0D0
      NDE=NDE+1
      IF (NDE.GT.10) THEN
         IERRORS = 1
         RETURN
      ENDIF
      GOTO 600

1190  IF (AA.EQ.A11) GO TO 1260
C
C *  Solution is finished at AA not equal to A1
C
      DO I=1,NEQ
         X0SW2(KC,I)=X0(I)
      END DO
      DO IE=1,NELEM
         DO IX=1,NGXW
            D0PSW2(KC,IE,IX)=D0P(IE,IX)
            D0NSW2(KC,IE,IX)=D0N(IE,IX)
            P0PSW2(KC,IE,IX)=P0P(IE,IX)
            P0NSW2(KC,IE,IX)=P0N(IE,IX)
            W0PSW2(KC,IE,IX)=W0P(IE,IX)
            W0NSW2(KC,IE,IX)=W0N(IE,IX)
         END DO
      END DO

      DO IE=1,NELEM
         DO IX=1,NGXW
            DO IY=1,NGYW
               S0SW2(KC,IE,IX,IY)=S0(IE,IX,IY)
               EPS0SW2(KC,IE,IX,IY)=EPS0(IE,IX,IY)
            END DO
         END DO
      END DO
      A0=AA
      AA=A11
      NDE=0
      GOTO 600

C
C  *  Solution is finished at AA equal to A1
C  *  Stores the current nodal displacement vector, 
C     residual wood deformations and pseudo nail stress and strain.
C
1260  A0=AA
      DO I=1,NEQ
         X0SW2(KC,I)=X0(I)
      END DO
      DO IE=1,NELEM
         DO IX=1,NGXW
            D0PSW2(KC,IE,IX)=D0P(IE,IX)
            D0NSW2(KC,IE,IX)=D0N(IE,IX)
            P0PSW2(KC,IE,IX)=P0P(IE,IX)
            P0NSW2(KC,IE,IX)=P0N(IE,IX)
            W0PSW2(KC,IE,IX)=W0P(IE,IX)
            W0NSW2(KC,IE,IX)=W0N(IE,IX)
            DO  IY=1,NGYW
               S0SW2(KC,IE,IX,IY)=S0(IE,IX,IY)
               EPS0SW2(KC,IE,IX,IY)=EPS0(IE,IX,IY)
            END DO
         END DO
      END DO
C      
C  *  Computes the force FY0 corresponding to current {x}.
C
      FY0=0.0
      DO JL=1,1
         DO  IE=IBEG(JL),IFIN(JL)
            DO  IX=1,NGXW
               IF(IDISP(IE).EQ.0) FY0=FY0+HGXW(IX)*P0(IE,IX)
     1         *SW0(IE,IX)*DELTA(KC,JL)/2.0
            END DO
         END DO
      ENDDO

1300  RETURN
      END


C    * -------------------------------------------------------------

      SUBROUTINE SHAPES (N0, N1, M0, M1, M2, NG, EG, DELTA)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      REAL*8 M0, M1, M2, N0, N1
      DIMENSION EG(9), M0(6,9), M1(6,9), M2(6,9)
      DIMENSION N0(6,9), N1(6,9)
      DO 10 I=1,NG
         X=EG(I)
         X2=X**2
         X3=X**3

C     SHAPES FOR W(X)
C 
         M0(1,I)=(2.0-3.0*X+X3)/4.0
         M0(2,I)=(1.0-X-X2+X3)*(DELTA/8.0)
         M0(3,I)=0.0D0
         M0(4,I)=(2.0+3.0*X-X3)/4.0
         M0(5,I)=(-1.0-X+X2+X3)*(DELTA/8.0)
         M0(6,I)=0.0D0

C 
C    SHAPES FOR PHI(X)
C 
         M1(1,I)=(-3.0+3.0*X2)*2.0/(4.0*DELTA)
         M1(2,I)=(-1.0-2.0*X+3.0*X2)/4.0
         M1(3,I)=0.0D0
         M1(4,I)=(3.0-3.0*X2)*2.0/(4.0*DELTA)
         M1(5,I)=(-1.0+2.0*X+3.0*X2)/4.0
         M1(6,I)=0.0D0

C     SHAPES FOR PHI'(X)
C 
         M2(1,I)=6.0*X/(DELTA**2)
         M2(2,I)=(-2.0+6.0*X)/(2.0*DELTA)
         M2(3,I)=0.0D0
         M2(4,I)=-6.0*X/(DELTA**2)
         M2(5,I)=(2.0+6.0*X)/(2.0*DELTA)
         M2(6,I)=0.0D0

C 
C     SHAPES FOR U(X)
C 
         N0(1,I)=0.0D0
         N0(2,I)=0.0D0
         N0(3,I)=(2.0D0-3.0D0*X+X3)/4.0D0
         N0(4,I)=0.0D0
         N0(5,I)=0.0D0
         N0(6,I)=(2.0D0+3.0D0*X-X3)/4.0D0
C 
C     SHAPES FOR U'(X)
C 
         N1(1,I)=0.0D0
         N1(2,I)=0.0D0
         N1(3,I)=(-3.0D0+3.0D0*X2)/(2.0D0*DELTA)
         N1(4,I)=0.0D0
         N1(5,I)=0.0D0
         N1(6,I)=(3.0D0-3.0D0*X2)/(2.0D0*DELTA)
 10   CONTINUE
      RETURN
      END

C 
C    * -------------------------------------------------------------   
C    
      SUBROUTINE PSUP (D0,W0,P0,W,Q0,Q1,Q4,XK,DMAX,SDF,PMAX,P,DPDW)

      IMPLICIT REAL*8(A-H,O-Z)

      DY=Q0/(XK-Q1)
      ALF=1.0D0
      IF (D0.GT.DY) THEN
         ALF=(DY/D0)**SDF
      ENDIF

      IF (W.LT.D0) THEN
         P=0.0D0
         DPDW=0.0D0
      ENDIF


      IF ((W.GE.D0).AND.(W.LE.W0))  THEN        !UNLOADING
         P=P0+XK*(W-W0)
         DPDW=XK
         IF (P.LT.0.0D0) THEN
            P=0.0D0
            DPDW=0.0D0
         ENDIF
      ENDIF

      IF ((W.GT.D0).AND.(W.GT.W0)) THEN         !LOADING OR RELOADING
         P=P0+ALF*XK*(W-W0)
         DPDW=ALF*XK
         D0=W-P/XK
         IF (W.LE.DMAX) THEN
            EX=DEXP(-XK*W/Q0)
            PUB=(Q0+Q1*W)*(1.0-EX)
            IF (P.GE.PUB) THEN
               P=PUB
               DPDW=Q1*(1.0-EX)+(Q0+Q1*W)*XK*EX/Q0
               D0=W-P/XK
            ENDIF
         ENDIF

         IF (W.GT.DMAX) THEN
            PUB=PMAX*DEXP(Q4*(W-DMAX)**2)
            IF (P.GE.PUB) THEN
               P=PUB
               DPDW=P*2.0D0*Q4*(W-DMAX)
               D0=W-P/XK
            ENDIF
         ENDIF

      ENDIF

      END

C
C    * Subroutine for the hysteresis loop of pseudo nail
C     (Checked OK!)
      SUBROUTINE STRESS (ST0, EP0, EPS, E, SY, S, DSDEPS)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      S=ST0+E*(EPS-EP0)
      DSDEPS=E
      S1=SY + 0.001D0*E*EPS
      S2=-SY + 0.001D0*E*EPS
      IF (S.GE.S1) DSDEPS=0.001D0*E
      IF (S.GE.S1) S=S1
      IF (S.LE.S2) DSDEPS=0.001D0*E
      IF (S.LE.S2) S=S2
      RETURN
      END


C  
C    * ---------------------------------------------------------
C
      SUBROUTINE DECOMPSW (N, LHB, A, IERROR)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(630)
      IERROR=0
      KB=LHB-1
      TEMP=A(1)
      IF (TEMP.LE.0.0D0) IERROR=1
      IF (IERROR.EQ.1) RETURN
      TEMP=DSQRT(TEMP)
      A(1)=TEMP
      DO 10 I=2,LHB
         A(I)=A(I)/TEMP
 10   CONTINUE
      DO 60 J=2,N
         J1=J-1
         IJD=LHB*J-KB
         SUM=A(IJD)
         KO=1
         IF (J.GT.LHB) KO=J-KB
         DO 20 K=KO,J1
            JK=KB*K+J-KB
            TEMP=A(JK)
            SUM=SUM-(TEMP**2)
 20      CONTINUE
         IF (SUM.LE.0.0D0) IERROR=1
         IF (IERROR.EQ.1) RETURN
         A(IJD)=DSQRT(SUM)
         DO 50 I=1,KB
            II=J+I
            KO=1
            IF (II.GT.LHB) KO=II-KB
            SUM=A(IJD+I)
            IF (I.EQ.KB) GO TO 40
            DO 30 K=KO,J1
               JK=KB*K+J-KB
               IK=KB*K+II-KB
               TEMP=A(JK)
               SUM=SUM-A(IK)*TEMP
 30         CONTINUE
 40         A(IJD+I)=SUM/A(IJD)
 50      CONTINUE
 60   CONTINUE
      RETURN
      END
C
C    ---------------------------------------------------------
C
      SUBROUTINE SOLVSW (N, LHB, A, B)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(630), B(63)

C     * FORWARD SUBSTITUTION *

      KB=LHB-1
      TEMP=A(1)
      B(1)=B(1)/TEMP
      DO 20 I=2,N
         I1=I-1
         KO=1
         IF (I.GT.LHB) KO=I-KB
         SUM=B(I)
         II=LHB*I-KB
         DO 10 K=KO,I1
            IK=KB*K+I-KB
            TEMP=A(IK)
            SUM=SUM-TEMP*B(K)
 10      CONTINUE
         B(I)=SUM/A(II)
 20   CONTINUE

C     * BACKWARD SUBSTITUTION *

      N1=N-1
      LB=LHB*N-KB
      TEMP=A(LB)
      B(N)=B(N)/TEMP
      DO 40 I=1,N1
         I1=N-I+1
         NI=N-I
         KO=N
         IF (I.GT.KB) KO=NI+KB
         SUM=B(NI)
         II=LHB*NI-KB
         DO 30 K=I1,KO
            IK=KB*NI+K-KB
            TEMP=A(IK)
            SUM=SUM-TEMP*B(K)
 30      CONTINUE
         B(NI)=SUM/A(II)
 40   CONTINUE
      RETURN
      END
C
C    * ----------------------------------------------------------------------
C    * IDRAW=0: DRAW UNDEFORMED STRUCTURE
C    * IDRAW=1: DRAW DEFORMED STRUCTURE
C    * IDRAW=2: DRAW UNDEFORMED AND DEFORMED FLOOR/ROOF DIAPHRAGM (PLAN VIEW)
C    * ----------------------------------------------------------------------

      SUBROUTINE DRAWMODEL (NF, NW, NMAS)
      USE DFLIB
      USE DFNLS
      USE COMMON_FRAME

      IMPLICIT REAL*8 (A-H, O-Z)
      REAL*8 N0,N1,M0,M1,M2

      COMMON /WALL/                           !SHEAR WALLS
     1N0(100,2,6,9),N1(100,2,6,9),
     1M0(100,2,6,9),M1(100,2,6,9),M2(100,2,6,9),
     1X0SW(100,63),D0PSW(100,20,9),D0NSW(100,20,9),
     1P0PSW(100,20,9),P0NSW(100,20,9),
     1W0PSW(100,20,9),W0NSW(100,20,9),
	1S0SW(100,20,9,16),EPS0SW(100,20,9,16),
     1X0SW2(100,63),D0PSW2(100,20,9),D0NSW2(100,20,9),
     1P0PSW2(100,20,9),P0NSW2(100,20,9),
     1W0PSW2(100,20,9),W0NSW2(100,20,9),
	1S0SW2(100,20,9,16),EPS0SW2(100,20,9,16),
     1FSW0(100),DFSW0(100),
	1FSW2(100),DFSW2(100),DSPSW2(100),
     1FSW(100) ,DFSW(100) ,DSPSW(100),DELTA(100,2),
	1IESW(100,2),IWTYPE(100),IWORIEN(100),IWDMG(100),IWDMG2(100),
	1NELSW,NELSW_1

      COMMON/PLOT/ NEP,NJOTOT,NODMASS(1000),FNODMAS(1000)

      DIMENSION AAT(3,3)
      DIMENSION XBL(21),YBL(21),ZBL(21)  !BASE GRID LINES
      TYPE (wxycoord) xy

      INTEGER(4) i4,BKCOLOR
      LOGICAL(4) result
      TYPE (qwinfo) qw
      TYPE (qwinfo) winfo


C    *   Maximize display windows
      OPEN(20, file='user', title='Plotting Model Configuration')
      qw.type = QWIN$MAX
      i4 = setwsizeqq(20, qw)

      write (20,100)
100   format (/,' ** Model Schematics **')


      ISCALE=25         !SCALE FACTOR 1:ISCALE

      ZMAX=0.0          !obtain building height
      DO I=1,NF
         IF (XC(IEF(I,1)).GE.ZMAX) ZMAX=XC(IEF(I,1))
         IF (XC(IEF(I,2)).GE.ZMAX) ZMAX=XC(IEF(I,2))
      ENDDO

      IF (ZMAX.LE.4000.0) iyoff=450     !ajustment of position of drawing on the screen
      IF (ZMAX.GT.4000.0.AND.ZMAX.LE.7000.0) iyoff=550
      IF (ZMAX.GT.7000.0) iyoff=650
      ixoff=150

      rotx=-75*3.14159265/180.0                 !X AXIS ROTATION -80 DEGREE
      roty=0*3.14159265/180.0                   !Y AXIS ROTATION 0 DEGREE
      rotz=-30*3.14159265/180.0                 !Z AXIS ROTATION -30 DEGREE

      sz=sin(rotz)
      cz=cos(rotz)
      sx=sin(rotx)
      cx=cos(rotx)
      sy=sin(roty)
      cy=cos(roty)

      AAT(1,1)=cy*cz
      AAT(1,2)=cy*sz
      AAT(1,3)=-sy
      AAT(2,1)=SX*SY*CZ-CX*SZ
      AAT(2,2)=SX*SY*SZ+CX*CZ
      AAT(2,3)=CY*SX
      AAT(3,1)=CX*SY*CZ+SX*SZ
      AAT(3,2)=CX*SY*SZ-SX*CZ
      AAT(3,3)=CY*CX


C    **  DRAW BASE GRID LINES

      DO I=1,21
         XBL(I)=-1200.0+(I-1)*1200.0
         YBL(I)=-1200.0+(I-1)*1200.0
         ZBL(I)=0.0
      ENDDO

      DO I=1,21
         XBB1=XBL(1)*AAT(1,1)-YBL(I)*AAT(1,2)+ZBL(I)*AAT(1,3)
         YBB1=XBL(1)*AAT(2,1)-YBL(I)*AAT(2,2)+ZBL(I)*AAT(2,3)
         ZBB1=XBL(1)*AAT(3,1)-YBL(I)*AAT(3,2)+ZBL(I)*AAT(3,3)

         XBB2=XBL(21)*AAT(1,1)-YBL(I)*AAT(1,2)+ZBL(I)*AAT(1,3)
         YBB2=XBL(21)*AAT(2,1)-YBL(I)*AAT(2,2)+ZBL(I)*AAT(2,3)
         ZBB2=XBL(21)*AAT(3,1)-YBL(I)*AAT(3,2)+ZBL(I)*AAT(3,3)

         XBB3=XBL(I)*AAT(1,1)-YBL(1)*AAT(1,2)+ZBL(I)*AAT(1,3)
         YBB3=XBL(I)*AAT(2,1)-YBL(1)*AAT(2,2)+ZBL(I)*AAT(2,3)
         ZBB3=XBL(I)*AAT(3,1)-YBL(1)*AAT(3,2)+ZBL(I)*AAT(3,3)

         XBB4=XBL(I)*AAT(1,1)-YBL(21)*AAT(1,2)+ZBL(I)*AAT(1,3)
         YBB4=XBL(I)*AAT(2,1)-YBL(21)*AAT(2,2)+ZBL(I)*AAT(2,3)
         ZBB4=XBL(I)*AAT(3,1)-YBL(21)*AAT(3,2)+ZBL(I)*AAT(3,3)

         RESULT = SETCOLORRGB(#778877)
c         CALL SETLINESTYLE( INT2( #AA3C ))
         CALL SETLINESTYLE( INT2( #AAAA ))
         CALL MOVETO_W (XBB1/ISCALE+ixoff, YBB1/ISCALE+iyoff,xy)
         status = LINETO_W (XBB2/ISCALE+ixoff,YBB2/ISCALE+iyoff)

         CALL MOVETO_W (XBB3/ISCALE+ixoff, YBB3/ISCALE+iyoff,xy)
         status = LINETO_W (XBB4/ISCALE+ixoff,YBB4/ISCALE+iyoff)
      ENDDO


C    **  !DRAW FRAMING MEMBERS	

      CALL SETLINESTYLE(#FFFF)
      DO I=1,NF
         XAT1=XA(IEF(I,1))*AAT(1,1)-XB(IEF(I,1))*AAT(1,2)
     1   +XC(IEF(I,1))*AAT(1,3)
         XBT1=XA(IEF(I,1))*AAT(2,1)-XB(IEF(I,1))*AAT(2,2)
     1   +XC(IEF(I,1))*AAT(2,3)
         XCT1=XA(IEF(I,1))*AAT(3,1)-XB(IEF(I,1))*AAT(3,2)
     1   +XC(IEF(I,1))*AAT(3,3)

         XAT2=XA(IEF(I,2))*AAT(1,1)-XB(IEF(I,2))*AAT(1,2)
     1   +XC(IEF(I,2))*AAT(1,3)
         XBT2=XA(IEF(I,2))*AAT(2,1)-XB(IEF(I,2))*AAT(2,2)
     1   +XC(IEF(I,2))*AAT(2,3)
         XCT2=XA(IEF(I,2))*AAT(3,1)-XB(IEF(I,2))*AAT(3,2)
     1   +XC(IEF(I,2))*AAT(3,3)

         WPL=ABS(XC(IEF(I,2))-XC(IEF(I,1)))

         IF (IETYP(I).EQ.2) THEN
            RESULT = SETCOLORRGB(#998877) ! light grey #998877
            CALL MOVETO_W ((XAT1/ISCALE+ixoff), (XBT1/ISCALE+iyoff),xy)
            status = LINETO_W ((XAT2/ISCALE+ixoff), (XBT2/ISCALE+iyoff))
         ENDIF

         IF (IETYP(I).EQ.1) THEN
            RESULT = SETCOLORRGB(#0FF0FF) ! white #FFFFFF; red #0000FF
            XDELTA=SQRT((XAT1-XAT2)**2+(XBT1-XBT2)**2)
            IF (XDELTA.GT.0.0) THEN
               XSIN=(XBT2-XBT1)/XDELTA
               XCOS=(XAT2-XAT1)/XDELTA
               IF (WPL.GT.2000.0) THEN
                  BWD=20.0
               ELSE
                  BWD=40.0
               ENDIF

               XT1=XAT1+BWD*XSIN
               YT1=XBT1-BWD*XCOS
               XT2=XAT1-BWD*XSIN
               YT2=XBT1+BWD*XCOS
               XT3=XAT2-BWD*XSIN
               YT3=XBT2+BWD*XCOS
               XT4=XAT2+BWD*XSIN
               YT4=XBT2-BWD*XCOS

               CALL MOVETO_W ((XT1/ISCALE+ixoff),(YT1/ISCALE+iyoff),xy)
               status = LINETO_W ((XT2/ISCALE+ixoff),(YT2/ISCALE+iyoff))
               CALL MOVETO_W ((XT2/ISCALE+ixoff),(YT2/ISCALE+iyoff),xy)
               status = LINETO_W ((XT3/ISCALE+ixoff),(YT3/ISCALE+iyoff))
               CALL MOVETO_W ((XT3/ISCALE+ixoff),(YT3/ISCALE+iyoff),xy)
               status = LINETO_W ((XT4/ISCALE+ixoff),(YT4/ISCALE+iyoff))
               CALL MOVETO_W ((XT4/ISCALE+ixoff),(YT4/ISCALE+iyoff),xy)
               status = LINETO_W ((XT1/ISCALE+ixoff),(YT1/ISCALE+iyoff))
            ENDIF
         ENDIF

      END DO

C    **  !DRAW SHEAR WALLS (DASHED DIAGONAL LINES)      

      DO I=1,NW
         XAT1=XA(IESW(I,1))*AAT(1,1)-XB(IESW(I,1))*AAT(1,2)
     1   +XC(IESW(I,1))*AAT(1,3)
         XBT1=XA(IESW(I,1))*AAT(2,1)-XB(IESW(I,1))*AAT(2,2)
     1   +XC(IESW(I,1))*AAT(2,3)
         XCT1=XA(IESW(I,1))*AAT(3,1)-XB(IESW(I,1))*AAT(3,2)
     1   +XC(IESW(I,1))*AAT(3,3)

         XAT2=XA(IESW(I,2))*AAT(1,1)-XB(IESW(I,2))*AAT(1,2)
     1   +XC(IESW(I,2))*AAT(1,3)
         XBT2=XA(IESW(I,2))*AAT(2,1)-XB(IESW(I,2))*AAT(2,2)
     1   +XC(IESW(I,2))*AAT(2,3)
         XCT2=XA(IESW(I,2))*AAT(3,1)-XB(IESW(I,2))*AAT(3,2)
     1   +XC(IESW(I,2))*AAT(3,3)
         RESULT = SETCOLORRGB(#00FF00) !SET GREEN COLOR #00FF00  RED #0000FF
         CALL SETLINESTYLE( INT2( #AA3C ))
         CALL MOVETO_W ((XAT1/ISCALE+ixoff), (XBT1/ISCALE+iyoff),xy)
         status = LINETO_W ((XAT2/ISCALE+ixoff), (XBT2/ISCALE+iyoff))
      END DO


C    **  !DRAW MASSES (DASHED DIAGONAL LINES)        
      DO I=1,NMAS

         IF (FNODMAS(I).GE.10E-6) THEN
            XAT1=XA(NODMASS(I))*AAT(1,1)-XB(NODMASS(I))*AAT(1,2)
     1      +XC(NODMASS(I))*AAT(1,3)
            XBT1=XA(NODMASS(I))*AAT(2,1)-XB(NODMASS(I))*AAT(2,2)
     1      +XC(NODMASS(I))*AAT(2,3)
            XCT1=XA(NODMASS(I))*AAT(3,1)-XB(NODMASS(I))*AAT(3,2)
     1      +XC(NODMASS(I))*AAT(3,3)

            wxtmp1=XAT1/ISCALE+ixoff-6.0
            wytmp1=XBT1/ISCALE+iyoff-6.0
            wxtmp2=XAT1/ISCALE+ixoff+6.0
            wytmp2=XBT1/ISCALE+iyoff+6.0
            status = SETCOLORRGB(#FFCC99)
            RESULT =
     1      ELLIPSE_w( $GFILLINTERIOR,wxtmp1,wytmp1,wxtmp2,wytmp2)

         ENDIF

      END DO

      END subroutine