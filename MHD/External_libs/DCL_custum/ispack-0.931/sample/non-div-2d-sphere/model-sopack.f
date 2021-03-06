************************************************************************
* ISPACK FORTRAN SUBROUTINE LIBRARY FOR SCIENTIFIC COMPUTING
* Copyright (C) 1998--2011 Keiichi Ishioka <ishioka@gfd-dennou.org>
*
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
* 
* You should have received a copy of the GNU Lesser General Public
* License along with this library; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
* 02110-1301 USA.
************************************************************************
***********************************************************************
*     MODEL TO STUDY SATURATION OF SHEAR INSTABILITY         1999/03/24
*            LEGENDRE POLYNOMIALS FLOWS AND NEW TYPE DISTURBANCE
*            BY RUNGE-KUTTA METHOD
***********************************************************************
      PROGRAM MODEL

      IMPLICIT REAL*8(A-H,O-Z)

      CALL SAMAIN(3,-0.07D0, 11,'output-sop.dat')

      END
************************************************************************
*     TEST MODEL FOR NEW ISPACK
************************************************************************
      SUBROUTINE SAMAIN(IPAT,APSI,IU,CF)

*     IPAT: 基本流をP_nモードとして与える場合のn
*     APSI: P_nモードの振幅
*     IU: 出力ファイルの装置番号
*     CF: 出力ファイル名
 
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER CF*(*)
      PARAMETER(PI=3.1415926535897932385D0)
      PARAMETER(MM=21)
      PARAMETER(LM=(MM+1)*(MM+1))
      PARAMETER(NB=23476)
      DIMENSION AVT(LM),W(LM*3)
      EXTERNAL SBDAVT,SBDISS

*---- OPEN SUBROUTINE PACKAGE AND INITIALIZE VARIABLES -----------------

      LEV=10                    !高階粘性項のラプラシアンの階数
      ITM=30                    !時間発展するステップ数
      NDV=20                    !Runge-Kuttaでのステップ分割数
      DT=1D0                    !ファイル出力の時間間隔

      DNU=10D0/(1D0*(MM*(MM+1)-2))**LEV
*      DNU=0

*/ ガウシアン型の擾乱のパラメター /*

      DA=1D-3     !擾乱の振幅
      DB=100D0    !擾乱の幅の逆数のルートに相当するパラメター
      X0=0        !擾乱の中心位置のλ座標
      Y0=PI/4     !擾乱の中心位置のφ座標

      CALL SBOPEN(LEV,DNU,NDV,DT)
      CALL SBINIT(AVT,IPAT,APSI,DA,DB,X0,Y0)

*---- OPEN FHPACK ------------------------------------------------------

      CALL FHUOPN(IU,CF,'W',NB)

*---- TIME EVOLUTION BY RUNGE-KUTTA METHOD -----------------------------

      I=0
      TIM=0
      CALL FEPUTS(IU,LM,AVT)
      CALL SBCHCK(I,AVT)

      DO I=1,ITM
        CALL TDRKNU(LM,NDV,DT,TIM,AVT,W,SBDISS,SBDAVT)
        CALL FEPUTS(IU,LM,AVT)
        CALL SBCHCK(I,AVT)
      END DO

*---- CLOSE FHPACK -----------------------------------------------------

      CALL FHUCLS(IU)

      END
************************************************************************
*     OPEN SUBROUTINE PACKAGE
************************************************************************
      SUBROUTINE SBOPEN(LEV,DNU,NDV,DT)

      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(MM=21,JM=32,IM=64,LM=(MM+1)*(MM+1))
      PARAMETER(ID=IM+1,JD=JM+1)
      PARAMETER(IW=IM+5,JW=JD)
      PARAMETER(PI=3.1415926535897932385D0)
      PARAMETER(SQRT3=1.7320508075688772935D0)
      PARAMETER(OMG=1)

      DIMENSION PSI((MM+2)*(MM+2)),AVT(LM),DAVT(LM)
      DIMENSION DN(LM),RN(LM),DRN(LM),IPM(LM),PM(LM)
      DIMENSION XG(IM),YG(JM)
      DIMENSION G(ID,JD)

      DIMENSION IT(5),T(IM*2),Y(JM/2,4)
      DIMENSION IP(((MM+2)/2+MM+2)*2),P(((MM+2)/2+MM+2)*JM)
      DIMENSION R(((MM+2)/2*2+3)*((MM+1)/2+1))
      DIMENSION IP2(2*((MM+2)/2+MM+2)*2),P2(2*((MM+2)/2+MM+2)*JM)
      DIMENSION R2(2*((MM+2)/2*2+3)*((MM+1)/2+1))
      DIMENSION IA((MM+2)*(MM+2)*4),A((MM+2)*(MM+2)*6)
      DIMENSION Q(2*((MM+2)/2+MM+2)*JM),WW(2*JW*IW),WS(2*JW*IW)
      DIMENSION IB((MM+1)*(MM+1),3),B((MM+1)*(MM+1),3)
      SAVE

      CALL SNINIT(MM+1,IM,JM,IT,T,Y,IP,P,R,IA,A)
      CALL SNKINI(MM+1,JM,2,IP,P,R,IP2,P2,R2)
      CALL SOINIT(MM,IB,B)

      DO L=1,LM
        CALL SNL2NM(L,N,M)
        CALL SNNM2L(N,-M,LD)
        IPM(L)=LD
        PM(L)=M
      END DO

      DO J=1,JM/2
        YG(JM/2+J  )= ASIN(Y(J,1))
        YG(JM/2-J+1)=-ASIN(Y(J,1))
      END DO

      DO J=1,IM
        XG(J)=2*PI*(J-1)/IM
      END DO

      DTDNU=DNU*DT/(2*NDV)

      RN(1)=0
      DRN(1)=1
      DN(1)=0
      DO N=1,MM
        DO M=-N,N
          CALL SNNM2L(N,M,L)
          DN(L)=-N*(N+1)
          RN(L)=-1D0/(N*(N+1))
          DRN(L)=EXP(-DTDNU*(1D0*N*(N+1)-2)**LEV)
        END DO
      END DO

      RETURN
*-----------------------------------------------------------------------
*     CALCULATION OF d(AVT)/dt
*-----------------------------------------------------------------------
      ENTRY SBDAVT(TIM,AVT,DAVT)

      DO L=1,LM
        PSI(L)=RN(L)*AVT(L)
      END DO
      PSI(3)=PSI(3)+OMG/SQRT3

      CALL SONDNL(MM,IM,ID,JM,JD,
     &    PSI,DAVT,IT,T,Y,IP2,P2,R2,IA,A,IB,B,Q,WS,WW)

      DO L=1,LM
        DAVT(L)=-2*PI*(DAVT(L)-2*PM(L)*PSI(IPM(L)))
      END DO

      RETURN
*-----------------------------------------------------------------------
*     CALCULATION OF THE EFFECT OF VISCOSITY
*-----------------------------------------------------------------------
      ENTRY SBDISS(TIM,DTIM,AVT)

      DO L=1,LM
        AVT(L)=DRN(L)*AVT(L)
      END DO

      RETURN
*-----------------------------------------------------------------------
*     CHECK ENERGY AND ENSTROPHY CONSERVATION
*-----------------------------------------------------------------------
      ENTRY SBCHCK(I,AVT)

      DO L=1,LM
        PSI(L)=RN(L)*AVT(L)
      END DO
      PSI(3)=PSI(3)+OMG/SQRT3

      ENE=0
      ENS=0
      DO L=2,LM
        ENE=ENE-DN(L)*PSI(L)*PSI(L)
        ENS=ENS+DN(L)*DN(L)*PSI(L)*PSI(L)
      END DO
      ENE=ENE/2
      ENS=ENS/2

      WRITE(6,'(I5,2F15.10)') I,ENE,ENS

      RETURN
*-----------------------------------------------------------------------
*     INITIALIZATION OF AVT
*-----------------------------------------------------------------------
      ENTRY SBINIT(AVT,IPAT,APSI,DA,DB,X0,Y0)

*     / BASIC FIELD /

      N=IPAT
      CALL BSSET0(LM,AVT)
      CALL SNNM2L(N,0,L)
      AVT(L)=-N*(N+1)*APSI
      AVT(3)=AVT(3)+2*OMG/SQRT3

*     / DISTURBANCE FIELD /

      CALL SCINID(JM,IM,JD,ID,G,XG,YG,DA,DB,X0,Y0)
      CALL SNTG2S(MM+1,IM,ID,JM,JD,1,
     &  G,PSI,IT,T,Y,IP,P,R,IA,A,Q,WS,WW,0,0)

*     / BASIC FIELD + DISTURBANCE FIELD /

      DO L=1,LM
        AVT(L)=AVT(L)+PSI(L)
      END DO

      END
************************************************************************
*     CALCULATE DISTURBANCE FIELD
************************************************************************
      SUBROUTINE SCINID(JM,IM,JD,ID,G,XG,YG,DA,DB,X0,Y0)

      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(PI=3.1415926535897932385D0)
      DIMENSION G(ID,JD)
      DIMENSION YG(JM),XG(IM)

      CB=COS(PI/2-Y0)
      SB=SIN(PI/2-Y0)

      DO I=1,IM
        CA=COS(XG(I)-X0)
        DO J=1,JM
          CC=COS(PI/2-YG(J))
          SC=SIN(PI/2-YG(J))
          CR=CB*CC+SB*SC*CA
          G(I,J)=DA*(EXP(DB*(CR-1))-(1-EXP(-2*DB))/(2*DB))
        END DO
      END DO

      END
