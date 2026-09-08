c     MUSC 2 Immediate postop
c     Ethan Kung  keo@ucsd.edu

c     Created by Mahdi Esmaily Moghadam 12-01-2010
c     Please report any problem to mesmaily@ucsd.edu, memt63@yahoo.com

c This subroutine initializes the parameters, you need to read the
c comments and specify them carefuly
c--------------------------------------------------------------------
c This is an example for RCR boundary condition with parameters
c Rd, R, and C which are distal, and proximal resistance and capacitor.

      SUBROUTINE INITIALIZE(nTimeStep)
      USE COM
      IMPLICIT NONE
      INTENT(OUT) nTimeStep

      LOGICAl ierr
      INTEGER i, nTimeStep
      REAL(KIND=8), ALLOCATABLE :: tZeroX(:)
c
c********************************************************************
c For instance if pressure in 3D solver is in cgs and here mmHg
c pConv=1334=1.334D3, also if flowrate in 3D solver is mm^3/s and
c here is mL/s qConv=1000=1D3. In the case both solver are using same
c unites you can set these two conversion coefficients to 1D0
      pConv = 1.334D3
      qConv = 1D0

c Only when all the surfaces of you model are coupled with NeumannSrfs
c you may set this to .TRUE.
      pCorr = .FALSE.
      qCorr = .FALSE.

c********************************************************************
c Block for your inputs

!      nbr = 24    !SET NUMBER OF PULMONARY BRANCHES#####
!      ninf = 4    !SET NUMBER OF INFERIOR BRANCHES#####
!      nX=19+nbr+ninf+2 !num of variables (19+nbr+12) ########

c These two value should match to that of defined in solver.inp
      nDirichletSrfs = 0
      nNeumannSrfs   = 8
c Number of unknowns that you need inside your lumped parameter network
      nUnknowns      = 37
c Number of time step between N and N+alpha, for higher values this code
c would be more accurate and more costly
      nTimeStep = 100

c Number of parameters to be printed in AllData file (the first
c nUnknowns columns are the "X" values)
      nXprint = 22

c--------------------------------------------------------------------
c You don't need to change this part                                 !
      ALLOCATE (tZeroX(nUnknowns), srfToXdPtr(nDirichletSrfs))       !
      ALLOCATE (srfToXPtr(nNeumannSrfs))                             !
      tZeroX = 0D0
c--------------------------------------------------------------------

      INCLUDE "initial_values_final.f"

c--------------------------------------------------------------------
c You don't need to change this part                                 !
      INQUIRE (FILE='InitialData', EXIST=ierr)                       !
      IF (.NOT.ierr) THEN                                            !
c         PRINT *, 'Initializing unknowns in LPM'                     !
         OPEN (1, FILE='InitialData',STATUS='NEW',FORM='UNFORMATTED')!
         WRITE (1) 0D0                                               !
         DO i=1, nUnknowns                                           !
            WRITE (1) tZeroX(i)                                      !
         END DO                                                      !
         CLOSE(1)                                                    !
      END IF                                                         !
c--------------------------------------------------------------------

c Surface to X pointer: this defines which Unknown corresponds to which
c suface in "List of Neumann Surfaces" inside solver.inp
c For example, if you have "List of Neumann Surfaces= 2 8 4 ...."
c and you set "srfToXPtr = (/5,3,9,.../)"
C this means X(5) corresponds to surface 2, X(3) <-> surface 8,
c and X(9) <-> surface 4
c Also, Q(1) corresponds to the first element in srfToXPtr, in this
c example, Q(1)<=>X(5), Q(2)<=>X(3)
      srfToXPtr  = (/32, 33, 34, 35, 36, 37,14, 5/)
C      srfToXdPtr = (/30, 31/)

      END SUBROUTINE INITIALIZE

c####################################################################
c Here you should find the f_i=dx_i/dt, based on the following parameters:
c  current x_i:                   x(i)
c  Current time:                  t
c  Flowrates from 3D code:        Q(i)
c  Pressure from Dirichlet faces: P(i)

      SUBROUTINE FINDF(t, x, f, Q, P)
      USE COM
      IMPLICIT NONE
      INTENT(IN) t, Q
      INTENT(OUT) f

c      INTEGER, PARAMETER :: nFaces= ninf + nbr + 1 = 29   !SET NUMBER OF INFERIOR BRANCHES#####

      REAL(KIND=8) t, x(nUnknowns), f(nUnknowns), Q(nNeumannSrfs),
     2   P(nDirichletSrfs)

      REAL(KIND=8) Qao, Psa, Psv, Tc, Tm, Tsad, Tr,
     2    Tmr, Tsr, Aa, Av, Ar, dAr, Pith, dPith,
     3    Tve, Pbar, elas, elase, exR,
     4    Aacon, Tmcon, Tcond, Tcon, Qsvc, Qinf,
     5    Vol, CdC, DVsa, DVsv, fsize, Psvp,
     6    Psva, fS, fD, fCAo, dt, Tma, Tmv, Qcoll_l, Qcoll_r,
     7    bloodvol, H, Rfaca, Qbolus, sigmoid, Hpump,
     8    Rfac, Qub, Qvub, Qlegc, Qlegv, Qabivc, Qthivc,
     9    Qla, Qlv, Qka, Qkv, Qia, Qiv, Qvpl, Qvpr,Qcoll

      COMPLEX(KIND=8) fdP(20), fP(20), fe(20), fee(30)
      INTEGER k

      REAL(KIND=8), parameter :: C(68) = (/
     &   0.00,      0.000,     0.0000,     0.0000, 
     &   11.56,     -0.190,     0.0000,     0.0000,
     &   0.30,      4.140,     0.4938,     0.0196, 
     &  -25.42,      4.270,     1.5450,     0.1466,
     &   16.86,     -5.275,     1.4200,     0.3568,
     &  -3.29,      1.045,     0.8912,     0.4848, 
     &   1.95,     -0.200,     0.9969,     0.6061, 
     &  -4.31,      0.545,     1.0400,     0.7314, 
     &  -0.11,     -1.070,     0.9738,     0.8615, 
     & -30.56,     -1.095,     0.7031,     0.9663, 
     &  23.76,    -12.560,    -1.0037,     0.9774, 
     &  32.02,     -3.630,    -3.0281,     0.7021, 
     & -15.49,      8.365,    -2.4356,     0.3294, 
     &   0.35,      2.575,    -1.0681,     0.1254, 
     &  -7.67,      2.700,    -0.4088,     0.0328, 
     &   3.04,     -0.190,    -0.0950,     0.0089, 
     &   0.00,      0.000,     0.0000,     0.0000 /)

      INCLUDE "parameters_final.f"

!     Time periods
      Tc   = 60.0/HR

      fsize = 1.0    ! volume scale factor for heart size
      fS = 0.797       ! systolic pressure factor
      fD = 1.0       ! diastolic pressure factor
      fCAo = 1.21     ! atrial arch compliance factor

!     Atrium
      Tma = MODULO(t+tsas-t1,Tc) ! atrium time
      IF (Tma .LT. Tsas) THEN
         Aa = (1.0 - COS(2.0*pi*Tma/Tsas))/2.0
      ELSE
         Aa = 0.0
      END IF
      DVsa = (x(27) - Vsa0)/fsize
      Psa = fS*Aa*DVsa/Ccsa + fD*csa*(EXP(dsa*DVsa) - 1D0)

! Ventricle.  normalized elastance, these functions peak at 0.3/T
      Tmv = t + 0.125 * Tsvs
      Tmv = MODULO(Tmv,Tc)   ! time within cardio cycle
      Tmv = Tmv/Tsvs     ! change to systole time scale
      IF (Tmv .LT. 1.875)  THEN    ! spline section for systole
         k = floor(Tmv/0.125)+1  ! segment, can be zero for FP error case Tmv<0
         dt = Tmv-(k-1)*0.125    ! time within spline segment
         k = 4*k;   ! row in 2-D spline coefficient array
         Av = ((C(k+1)*dt+C(k+2))*dt+C(k+3))*dt+C(k+4)   ! calculate spline
      ELSE  ! in extended diastole section
         Av = 0.0
      END IF

      DVsv=(x(29) - Vsv0)/fsize        ! size adjustment
      Psvp = csv*(EXP(dsv*DVsv)-1.0)   ! passive curve
      Psva = (a*DVsv+b)*DVsv           ! active curve (parabola)
      Psv = fS*Av*Psva + fD*Psvp    ! theoretical ventricle pressure

!     Respiration
      Tr = 6D1/fr
      Tmr = MOD(t,Tr)
      Tsr = Trr*Tr  !don't need this if prescribing fourier waveform

      IF (Tmr .LT. Tsr) THEN
         Ar  = 5D-1*(1D0 - COS(2D0*pi*Tmr/Tsr))
         dAr = (pi/Tsr)*SIN(2D0*pi*Tmr/Tsr)
      ELSE
         Ar  = 0D0
         dAr = 0D0
      END IF
      Pith  = APith*Ar + P0ith
      dPith = APith*dAr
      Pith = 0
      dPith = 0
      

!!!!!! Now connect Pith to the heart pressures !!!
      Psv = Psv + Pith
      Psa = Psa + Pith

!     This is the only case that x is directly manipulated
!     Atrial-Ventricular valve flow
      IF (x(28) .LT. 0D0) x(28) = 0D0 

!     Aortic flow with effective pressure
      IF (   Psv .GT. x(1)) THEN
          Qao = (SQRT(Rmyo*Rmyo + 4D0*Kao*( Psv - x(1))) - Rmyo)/2D0/Kao
      ELSE
       Qao = 0D0
      END IF
      

      Rfac = 1D0   

      !     resistive flows  
      Qub = (x(3) - x(4))/(Rub*Rfac)  ! upper body (capillaries)
      Qvub = (x(4) - x(5))/Rubv       ! upper body veins
      Qsvc = -Q(8)  ! SVC
      Qlegc = (x(11) - x(12))/(Rlegc*Rfac)
      Qlegv = (abs(x(12) - x(13)) + x(12) - x(13))/Rlegv/2D0
      Qabivc = (x(13) - x(14))/Rabivc

      Qthivc = -Q(7)         ! thoracic IVC flow with Dirichlet inertance
      
      Qla=(x(7) - x(15))/Rla     ! liver
      Qlv=(x(15) - x(14))/Rlv
      Qia=(x(9) - x(16))/Ria     ! intestines
      Qiv=(x(16) - x(15))/Riv    
      Qka=(x(7) - x(17))/Rka     ! kidneys
      Qkv=(x(17) - x(14))/Rkv
      Qvpr=(x(22) - Psa)/Rvpr  ! RPV flow
      Qvpl=(x(26) - Psa)/Rvpl  ! LPV flow

      Qcoll_r = coll_split * (x(7) - x(22)) / Rcoll
      Qcoll_l = (1-coll_split) * (x(7) - x(26)) / Rcoll
      Qcoll = Qcoll_r + Qcoll_l
      
      !********************************************************************
      
      !     The main body of equations 
      f(1) = (Qao - x(2) - x(6))/(fCAo*Cao) + dPith    ! aorta arch pressure
      f(2)  = (x(1) - x(2)*Ruba - x(3))/Luba    ! arteries upper body flow
      f(3)  = (x(2) - Qub)/Caub                 ! arteries upper body pressure
      f(4)  = (Qub - Qvub)/Cvub                 ! veins upper body pressure
      f(5)  = (Qvub - Qsvc)/Csvc + dPith        ! SVC pressure
      f(6) = (x(1) - x(6)*Rthao - x(7))/Lthao   ! thoracic Aorta flow
      f(7) = (x(6) - Qcoll- x(8) -Qla - Qka)/Cthao + dPith   ! thAo pressure
      f(8) = (x(7) - x(8)*Rabao - x(9))/Labao   ! abdominal Aorta flow
      f(9) = (x(8) - Qia - x(10))/Cabao         ! abdAo pressure
      f(10) = (x(9) - x(10)*Rlega - x(11))/Llega   ! leg arteries flow
      f(11) = (x(10) - Qlegc)/Clega             ! leg arteries pressure
      f(12) = (Qlegc - Qlegv)/Clegv             ! leg capillaries pressure
      f(13) = (Qlegv - Qabivc)/Cabivc           ! abdominal IVC pressure
      f(14) = (Qabivc + Qkv + Qlv - Qthivc)/Cthivc + dPith  ! thoracic IVC Press
      f(15) = (Qiv + Qla - Qlv)/Cl              ! liver pressure
      f(16) = (Qia - Qiv)/Ci                    ! intestines pressure
      f(17) = (Qka - Qkv)/Ck                    ! kidneys pressure

      f(18) = (Qsvc + Q(1) - x(19) - x(23))/Ctcpc + dPith   ! Fontan pressure

      f(33) = (Q(2) - x(19)*0.48)/Crpa                     ! pressure node coupled to RPA surface 1
      f(34) = (Q(3) - x(19)*0.13)/Crpa                     ! pressure node coupled to RPA surface 2
      f(35) = (Q(4) - x(19)*0.15)/Crpa                     ! pressure node coupled to RPA surface 3
      f(36) = (Q(5) - x(19)*0.17)/Crpa                     ! pressure node coupled to RPA surface 4
      f(37) = (Q(6) - x(19)*0.07)/Crpa                     ! pressure node coupled to RPA surface 5

      f(19) = (SUM(x(33:37))/5 - x(19)*Rapolr - x(20))/Lapolr   ! RPA flow
      f(20) = (x(19) - x(21))/Clungr + dPith          ! right lung pressure
      f(21) = (x(20) - x(21)*Rlungr - x(22))/Llungr   ! right lung flow
      f(22) = (x(21) + Qcoll_r- Qvpr)/Cvpr + dPith             ! RPV pressure


      f(32) = (Q(1) - x(23))/Clpa                     ! pressure node coupled to LPA surface
      f(23) = (x(32) - x(23)*Rapoll - x(24))/Lapoll   ! LPA flow
      f(24) = (x(23) - x(25))/Clungl + dPith          ! left lung pressure
      f(25) = (x(24) - x(25)*Rlungl - x(26))/Llungl   ! left lung flow
      f(26) = (x(25) + Qcoll_l- Qvpl)/Cvpl + dPith             ! LPV pressure

      f(27) = Qvpr + Qvpl - x(28)   ! Atrium Volume
      
      ! A-V valve flow acceleration
      IF (x(28).EQ.0D0 .AND. Psa.LT.(Psv-Rmyo*Qao)) THEN
          f(28) = 0D0
      ELSE
          f(28) = (Psa - Psv - Kav*x(28)*x(28))/Lav
      END IF

      f(29) = x(28) - Qao  ! ventricle volume
!      f(30) = -(x(14) - x(30)*Rthivc - P(1))/Llc  ! ThIVC flow w/Dirichlet L

!      IF (x(30) .LE. 0D0 .AND. (x(14) - x(30)*Rthivc).LT. P(1)) THEN
!            f(44) = 0D0
!      END IF

!      f(31) = -(x(5) - x(31)*Rsvc - P(2))/Luc  ! SVC flow w/Dirichlet L



!     Assign the additional parameters to be printed
      Xprint(1) = t ! 38
      Xprint(2) = Aa !39
      Xprint(3) = Psa
      Xprint(4) = Qcoll
      Xprint(5) = Psvp
      Xprint(6) = Psva
      Xprint(7) = Psv
      Xprint(8) = Qao
      Xprint(9) = Qub
      Xprint(10) = Qvub
      Xprint(11) = Qsvc
      Xprint(12) = Qlegc
      Xprint(13) = Qlegv
      Xprint(14) = Qabivc
      Xprint(15) = Qla
      Xprint(16) = Qlv
      Xprint(17) = Qvpr
      Xprint(18) = Qvpl
      Xprint(19) = Q(7)
      Xprint(20) = Q(8)
      Xprint(21) = Q(1)
      Xprint(22) = SUM(Q(2:6))



      RETURN
      END SUBROUTINE FINDF


