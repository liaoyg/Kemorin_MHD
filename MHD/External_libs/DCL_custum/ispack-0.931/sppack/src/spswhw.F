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
************************************************************************
*   ISPACK FORTRAN SUBROUTINE LIBRARY FOR SCIENTIFIC COMPUTING         *
*   Copyright (C) 1999 Keiichi Ishioka                                 *
*                                                                      *
*   This library is free software; you can redistribute it and/or      *
*   modify it under the terms of the GNU Library General Public        *
*   License as published by the Free Software Foundation; either       *
*   version 2 of the License, or (at your option) any later version.   *
*                                                                      *
*   This library is distributed in the hope that it will be useful,    *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of     *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  *
*   Library General Public License for more details.                   *
*                                                                      *
*   You should have received a copy of the GNU Library General Public  *
*   License along with this library; if not, write to the Free         *
*   Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. *
************************************************************************
************************************************************************
*     浅水方程式の非線形項の計算
*     (角運動量を保存する高階粘性項も含む)
*     (若干高速化を図った版)                                  2000/08/16
************************************************************************
      SUBROUTINE SPSWHW(MM,IM,ID,JM,JD,OMEGA,BARPHI,DNU,ALPHA,
     & LEV,AVT,DIV,PHI,DAVT,DDIV,DPHI,RN,IRM,IT,T,Y,IP,P,R,IA,A,Q,WS,WW)

      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AVT((MM+1)*(MM+1))
      DIMENSION DIV((MM+1)*(MM+1))
      DIMENSION PHI((MM+1)*(MM+1))
      DIMENSION DAVT((MM+1)*(MM+1))
      DIMENSION DDIV((MM+1)*(MM+1))
      DIMENSION DPHI((MM+1)*(MM+1))
      DIMENSION RN((MM+1)*(MM+1),2)
      DIMENSION IRM((MM+1)*(MM+1),2)      
      DIMENSION IT(5),T(IM*2)
      DIMENSION Y(JM*2)
      DIMENSION IP(((MM+1)/2+MM+1)*2)
      DIMENSION P(((MM+1)/2+MM+1)*JM)
      DIMENSION R(((MM+1)/2*2+3)*(MM/2+1))
      DIMENSION IA((MM+1)*(MM+1)*4)
      DIMENSION A((MM+1)*(MM+1)*6)
      DIMENSION Q(((MM+1)/2+MM+1)*JM)
      DIMENSION WS(*)      
      DIMENSION WW(*)

      LMD=((MM+1)/2*2+3)*(MM/2+2)*2
      LMD2=JD*((MM+1)/2+MM+1)*2
      MAXDIM=MAX(ID*JD,LMD,LMD2)

*/ 実際の操作 */
      CALL SPSWHS(MAXDIM,MM,IM,ID,JM,JD,OMEGA,BARPHI,DNU,ALPHA,
     & LEV,AVT,DIV,PHI,DAVT,DDIV,DPHI,RN,IRM,IT,T,Y,IP,P,R,IA,A,Q,WS,WW)

      END
************************************************************************
      SUBROUTINE SPSWHS(MAXDIM,MM,IM,ID,JM,JD,OMEGA,BARPHI,DNU,ALPHA,
     & LEV,AVT,DIV,PHI,DAVT,DDIV,DPHI,RN,IRM,IT,T,Y,IP,P,R,IA,A,Q,WS,WW)

      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(SQRT3=1.7320508075688772935D0)
      DIMENSION AVT((MM+1)*(MM+1))
      DIMENSION DIV((MM+1)*(MM+1))
      DIMENSION PHI((MM+1)*(MM+1))
      DIMENSION DAVT((MM+1)*(MM+1))
      DIMENSION DDIV((MM+1)*(MM+1))
      DIMENSION DPHI((MM+1)*(MM+1))
      DIMENSION RN((MM+1)*(MM+1),2)
      DIMENSION IRM((MM+1)*(MM+1),2)
      DIMENSION IT(5),T(IM*2)
      DIMENSION Y(JM*2)
      DIMENSION IP(((MM+1)/2+MM+1)*2)
      DIMENSION P(((MM+1)/2+MM+1)*JM)
      DIMENSION R(((MM+1)/2*2+3)*(MM/2+1))
      DIMENSION IA((MM+1)*(MM+1),4)
      DIMENSION A((MM+1)*(MM+1),6)
      DIMENSION Q(((MM+1)/2+MM+1)*JM)
      DIMENSION WW(MAXDIM,11)
      DIMENSION WS(MAXDIM,11)

      LM=(MM+1)*(MM+1)
      LMD=((MM+1)/2*2+3)*(MM/2+2)*2
      DNUD=DNU/BARPHI

*/ ψとχの計算(一時的に, ψとしてWS(*,3)を, χとしてWS(*,4)を使う)
*/ 静止系から見ることにする.
      
      CALL BSSET0(LMD,WS(1,3))
      CALL BSSET0(LMD,WS(1,4))      
      DO L=1,LM
        WS(L,3)=RN(L,2)*AVT(L)
        WS(L,4)=RN(L,2)*DIV(L)
      END DO

*/ ｕの計算(WS(*,1))

      CALL BSSET0(LMD,WS(1,7))
      CALL BSSET0(LMD,WS(1,8))
      CALL BSSET0(LMD,WS(1,1))
      DO L=1,LM
        WS(IA(L,2),1)=A(L,2)*WS(L,4)
        WS(IA(L,3),7)=-A(L,3)*WS(L,3)
        WS(IA(L,4),8)=-A(L,4)*WS(L,3)
      END DO
      DO L=1,LMD
        WS(L,1)=WS(L,1)+WS(L,7)+WS(L,8)
      END DO

*/ ｖの計算(WS(*,2))

      CALL BSSET0(LMD,WS(1,7))
      CALL BSSET0(LMD,WS(1,8))
      CALL BSSET0(LMD,WS(1,2))
      DO L=1,LM
        WS(IA(L,2),2)=A(L,2)*WS(L,3)
        WS(IA(L,3),7)=A(L,3)*WS(L,4)
        WS(IA(L,4),8)=A(L,4)*WS(L,4)
      END DO
      DO L=1,LMD
        WS(L,2)=WS(L,2)+WS(L,7)+WS(L,8)
      END DO

*/ u^, v^ の計算の準備

      DO I=1,LEV
        DO L=1,LM
          WS(L,3)=-(RN(L,1)+2)*WS(L,3)
          WS(L,4)=-(RN(L,1)+2)*WS(L,4)
        END DO
      END DO

*/ u^の計算(WS(*,5))

      CALL BSSET0(LMD,WS(1,7))
      CALL BSSET0(LMD,WS(1,8))
      CALL BSSET0(LMD,WS(1,5))
      DO L=1,LM
        WS(IA(L,2),5)=A(L,2)*WS(L,4)
        WS(IA(L,3),7)=-A(L,3)*WS(L,3)
        WS(IA(L,4),8)=-A(L,4)*WS(L,3)
      END DO
      DO L=1,LMD
        WS(L,5)=WS(L,5)+WS(L,7)+WS(L,8)
      END DO

*/ v^の計算(WS(*,6))

      CALL BSSET0(LMD,WS(1,7))
      CALL BSSET0(LMD,WS(1,8))
      CALL BSSET0(LMD,WS(1,6))
      DO L=1,LM
        WS(IA(L,2),6)=A(L,2)*WS(L,3)
        WS(IA(L,3),7)=A(L,3)*WS(L,4)
        WS(IA(L,4),8)=A(L,4)*WS(L,4)
      END DO
      DO L=1,LMD
        WS(L,6)=WS(L,6)+WS(L,7)+WS(L,8)
      END DO

*/ ∂Φ/∂λ(WS(*,3)), cosφ∂Φ/∂φ(WS(*,4))の計算

      CALL BSSET0(LMD,WS(1,3))
      CALL BSSET0(LMD,WS(1,4))
      DO L=1,LM
        WS(IA(L,3),3)=A(L,3)*PHI(L)
        WS(IA(L,4),4)=A(L,4)*PHI(L)
      END DO
      DO L=1,LMD
        WS(L,4)=WS(L,4)+WS(L,3)
      END DO

      CALL BSSET0(LMD,WS(1,3))
      DO L=1,LM
        WS(IA(L,2),3)=A(L,2)*PHI(L)
      END DO

*/ ｑ(WS(*,7))とΦ(WS(*,8))と△Φ(WS(*,9))とq^(WS(*,10))とD^(WS(*,11))
*/          の計算の準備

      DO L=1,LM
        WS(L,7)=AVT(L)
        WS(L,8)=DIV(L)
      END DO
      
      DO I=1,LEV
        DO L=1,LM
          WS(L,7)=-(RN(L,1)+2)*WS(L,7)
          WS(L,8)=-(RN(L,1)+2)*WS(L,8)
        END DO
      END DO

      CALL BSSET0(LMD,WS(1,10))
      CALL BSSET0(LMD,WS(1,11))
      DO L=1,LM
        WS(IA(L,1),10)=A(L,1)*WS(L,7)
        WS(IA(L,1),11)=A(L,1)*WS(L,8)
      END DO
      
      CALL BSSET0(LMD,WS(1,7))
      CALL BSSET0(LMD,WS(1,8))
      CALL BSSET0(LMD,WS(1,9))      
      DO L=1,LM
        WS(IA(L,1),7)=A(L,1)*AVT(L)
        WS(IA(L,1),8)=A(L,1)*PHI(L)
        WS(IA(L,1),9)=A(L,1)*RN(L,1)*PHI(L)
      END DO

*/  ルジャンドル変換
      
      CALL SPLSGB(MM,JM,1,
     &  WS(1,1),WS(1,2),WS(1,3),WS(1,4),WS(1,5),WS(1,6),WS(1,7),WS(1,8),
     &  WS(1,9),WS(1,10),WS(1,11),
     &  WW(1,1),WW(1,2),WW(1,3),WW(1,4),WW(1,5),WW(1,6),WW(1,7),WW(1,8),
     &  WW(1,9),WW(1,10),WW(1,11),
     &  Y,P,R,Q)

*/ スペクトル→グリッド
      DO IV=1,6
*/      パリティ変換
        CALL SNPS2G(MM,JM,JD,1,WW(1,IV),WS(1,IV),IP,Y,1)
*/      フーリエ変換
        CALL SNFS2G(MM,IM,JD,1,WS(1,IV),WW(1,IV),IT,T)
*/      添字の並べ替え
        CALL SNGS2G(IM,ID,JD,1,WW(1,IV),WS(1,IV))
      END DO

      DO IV=7,11
*/      パリティ変換
        CALL SNPS2G(MM,JM,JD,1,WW(1,IV),WS(1,IV),IP,Y,0)
*/      フーリエ変換
        CALL SNFS2G(MM,IM,JD,1,WS(1,IV),WW(1,IV),IT,T)
*/      添字の並べ替え
        CALL SNGS2G(IM,ID,JD,1,WW(1,IV),WS(1,IV))
      END DO

*/ 非線形項の計算

      DO IJ=1,ID*JD
        U=WS(IJ,1)
        V=WS(IJ,2)
        PHID=WS(IJ,8)-BARPHI
        WS(IJ,1)=WS(IJ,7)*U
     &    +DNUD*( 2*WS(IJ,9)*WS(IJ,6)+WS(IJ,10)*WS(IJ,3)
     &    +ALPHA*WS(IJ,11)*WS(IJ,4))
        WS(IJ,2)=WS(IJ,7)*V
     &    +DNUD*(-2*WS(IJ,9)*WS(IJ,5)+WS(IJ,10)*WS(IJ,4)
     &    -ALPHA*WS(IJ,11)*WS(IJ,3))
        WS(IJ,3)=PHID*U
        WS(IJ,4)=PHID*V
        WS(IJ,5)=PHID*WS(IJ,5)
        WS(IJ,6)=PHID*WS(IJ,6)
        WS(IJ,7)=PHID*WS(IJ,10)        
        WS(IJ,8)=(U*U+V*V)*0.5D0+DNUD*ALPHA*PHID*WS(IJ,11)
      END DO

*/ グリッド→スペクトル(IV=5,6の成分だけは SNPG2Sの変換が異なるので注意)
      
      DO IV=1,6
*/      添字の並べ替え
        CALL SNGG2S(IM,ID,JD,1,WS(1,IV),WW(1,IV))
*/      フーリエ変換
        CALL SNFG2S(MM,IM,JD,1,WW(1,IV),WS(1,IV),IT,T)
*/      パリティ変換
        CALL SNPG2S(MM,JM,JD,1,WS(1,IV),WW(1,IV),IP,Y,1)
      END DO
 
      DO IV=7,8
*/      添字の並べ替え
        CALL SNGG2S(IM,ID,JD,1,WS(1,IV),WW(1,IV))
*/      フーリエ変換
        CALL SNFG2S(MM,IM,JD,1,WW(1,IV),WS(1,IV),IT,T)
*/      パリティ変換
        CALL SNPG2S(MM,JM,JD,1,WS(1,IV),WW(1,IV),IP,Y,0)
      END DO

*/      ルジャンドル変換      

      CALL SPLGS8(MM,JM,1,
     &  WW(1,1),WW(1,2),WW(1,3),WW(1,4),WW(1,5),WW(1,6),WW(1,7),WW(1,8),
     &  WS(1,1),WS(1,2),WS(1,3),WS(1,4),WS(1,5),WS(1,6),WS(1,7),WS(1,8),
     &  Y,P,R,Q)
 
      DO L=1,LM
        DAVT(L)=A(L,2)*WS(IA(L,2),1)
     &    +A(L,3)*WS(IA(L,3),2)+A(L,4)*WS(IA(L,4),2)
     &    -DNUD*RN(L,1)*A(L,1)*WS(IA(L,1),7)
     &    +DNUD*2*(RN(L,1)+1)*(-A(L,2)*WS(IA(L,2),6)
     &    +A(L,3)*WS(IA(L,3),5)+A(L,4)*WS(IA(L,4),5))
        DDIV(L)=-A(L,2)*WS(IA(L,2),2)
     &    +A(L,3)*WS(IA(L,3),1)+A(L,4)*WS(IA(L,4),1)
     &    -RN(L,1)*A(L,1)*WS(IA(L,1),8)
     &    +DNUD*2*(RN(L,1)+1)*(-A(L,2)*WS(IA(L,2),5)
     &    -A(L,3)*WS(IA(L,3),6)-A(L,4)*WS(IA(L,4),6))
        DPHI(L)=A(L,2)*WS(IA(L,2),3)
     &    +A(L,3)*WS(IA(L,3),4)+A(L,4)*WS(IA(L,4),4)
      END DO

*/ 回転系に戻るために以下の変換が必要
      
      DO L=1,LM
        DAVT(L)=DAVT(L)-OMEGA*IRM(L,2)*AVT(IRM(L,1))
        DDIV(L)=DDIV(L)-OMEGA*IRM(L,2)*DIV(IRM(L,1))
        DPHI(L)=DPHI(L)-OMEGA*IRM(L,2)*PHI(IRM(L,1))
      END DO
      DDIV(7)=DDIV(7)-OMEGA*OMEGA*2/SQRT(5D0)

      END
***********************************************************************
      SUBROUTINE SPLSGB(MM,JM,KM,
     &  S1,S2,S3,S4,S5,S6,S7,S8,S9,SA,SB,
     &  W1,W2,W3,W4,W5,W6,W7,W8,W9,WA,WB,
     &  Y,P,R,Q)
*-----------------------------------------------------------------------
* 以下にある ifdef によるループの分解は, 今一つ賢くない SX4 のコンパイラ
* でも外側ループのアンローリングができるようにするためのもので, VPでは,
* 不要(むしろ有害)である.
*-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION S1(KM*((MM+1)/2*2+3),0:MM/2+1,2)
      DIMENSION S2(KM*((MM+1)/2*2+3),0:MM/2+1,2)
      DIMENSION S3(KM*((MM+1)/2*2+3),0:MM/2+1,2)
      DIMENSION S4(KM*((MM+1)/2*2+3),0:MM/2+1,2)
      DIMENSION S5(KM*((MM+1)/2*2+3),0:MM/2+1,2)
      DIMENSION S6(KM*((MM+1)/2*2+3),0:MM/2+1,2)
      DIMENSION S7(KM*((MM+1)/2*2+3),0:MM/2+1,2)
      DIMENSION S8(KM*((MM+1)/2*2+3),0:MM/2+1,2)
      DIMENSION S9(KM*((MM+1)/2*2+3),0:MM/2+1,2)
      DIMENSION SA(KM*((MM+1)/2*2+3),0:MM/2+1,2)
      DIMENSION SB(KM*((MM+1)/2*2+3),0:MM/2+1,2)      
      DIMENSION W1(KM*((MM+1)/2+MM+1),2,JM/2,2)
      DIMENSION W2(KM*((MM+1)/2+MM+1),2,JM/2,2)
      DIMENSION W3(KM*((MM+1)/2+MM+1),2,JM/2,2)
      DIMENSION W4(KM*((MM+1)/2+MM+1),2,JM/2,2)
      DIMENSION W5(KM*((MM+1)/2+MM+1),2,JM/2,2)
      DIMENSION W6(KM*((MM+1)/2+MM+1),2,JM/2,2)
      DIMENSION W7(KM*((MM+1)/2+MM+1),2,JM/2,2)
      DIMENSION W8(KM*((MM+1)/2+MM+1),2,JM/2,2)                  
      DIMENSION W9(KM*((MM+1)/2+MM+1),2,JM/2,2)
      DIMENSION WA(KM*((MM+1)/2+MM+1),2,JM/2,2)
      DIMENSION WB(KM*((MM+1)/2+MM+1),2,JM/2,2)                  
      DIMENSION Y(JM/2)
      DIMENSION P(KM*((MM+1)/2+MM+1)*JM)
      DIMENSION R(KM*((MM+1)/2*2+3),0:MM/2)
      DIMENSION Q(KM*((MM+1)/2+MM+1),2,JM/2)

      MMP=(MM+1)/2
      MMD=MM/2
      JH=JM/2

      CALL BSCOPY(KM*((MM+1)/2+MM+1)*JM,P,Q)
      CALL BSSET0(KM*(MM+1+MMP)*2*JM,W1)
      CALL BSSET0(KM*(MM+1+MMP)*2*JM,W2)
      CALL BSSET0(KM*(MM+1+MMP)*2*JM,W3)      
      CALL BSSET0(KM*(MM+1+MMP)*2*JM,W4)
      CALL BSSET0(KM*(MM+1+MMP)*2*JM,W5)
      CALL BSSET0(KM*(MM+1+MMP)*2*JM,W6)
      CALL BSSET0(KM*(MM+1+MMP)*2*JM,W7)
      CALL BSSET0(KM*(MM+1+MMP)*2*JM,W8)
      CALL BSSET0(KM*(MM+1+MMP)*2*JM,W9)
      CALL BSSET0(KM*(MM+1+MMP)*2*JM,WA)
      CALL BSSET0(KM*(MM+1+MMP)*2*JM,WB)
      
      L=0
      DO J=1,JH
        DO K=1,KM*(2*MMP+2)
          M=K+L*KM
          W1(M,1,J,1)=W1(M,1,J,1)+S1(K+KM,L,1)*Q(M,1,J)
          W1(M,2,J,1)=W1(M,2,J,1)+S1(K+KM,L,2)*Q(M,1,J)
          W2(M,1,J,1)=W2(M,1,J,1)+S2(K+KM,L,1)*Q(M,1,J)
          W2(M,2,J,1)=W2(M,2,J,1)+S2(K+KM,L,2)*Q(M,1,J)
          W3(M,1,J,1)=W3(M,1,J,1)+S3(K+KM,L,1)*Q(M,1,J)
          W3(M,2,J,1)=W3(M,2,J,1)+S3(K+KM,L,2)*Q(M,1,J)
          W4(M,1,J,1)=W4(M,1,J,1)+S4(K+KM,L,1)*Q(M,1,J)
          W4(M,2,J,1)=W4(M,2,J,1)+S4(K+KM,L,2)*Q(M,1,J)
          W5(M,1,J,1)=W5(M,1,J,1)+S5(K+KM,L,1)*Q(M,1,J)
          W5(M,2,J,1)=W5(M,2,J,1)+S5(K+KM,L,2)*Q(M,1,J)
          W6(M,1,J,1)=W6(M,1,J,1)+S6(K+KM,L,1)*Q(M,1,J)
          W6(M,2,J,1)=W6(M,2,J,1)+S6(K+KM,L,2)*Q(M,1,J)
          W7(M,1,J,1)=W7(M,1,J,1)+S7(K+KM,L,1)*Q(M,1,J)
          W7(M,2,J,1)=W7(M,2,J,1)+S7(K+KM,L,2)*Q(M,1,J)
          W8(M,1,J,1)=W8(M,1,J,1)+S8(K+KM,L,1)*Q(M,1,J)
          W8(M,2,J,1)=W8(M,2,J,1)+S8(K+KM,L,2)*Q(M,1,J)
          W9(M,1,J,1)=W9(M,1,J,1)+S9(K+KM,L,1)*Q(M,1,J)
          W9(M,2,J,1)=W9(M,2,J,1)+S9(K+KM,L,2)*Q(M,1,J)
          WA(M,1,J,1)=WA(M,1,J,1)+SA(K+KM,L,1)*Q(M,1,J)
          WA(M,2,J,1)=WA(M,2,J,1)+SA(K+KM,L,2)*Q(M,1,J)
          WB(M,1,J,1)=WB(M,1,J,1)+SB(K+KM,L,1)*Q(M,1,J)
          WB(M,2,J,1)=WB(M,2,J,1)+SB(K+KM,L,2)*Q(M,1,J)
          W1(M,1,J,2)=W1(M,1,J,2)+S1(K,L+1,1)*Q(M,2,J)
          W1(M,2,J,2)=W1(M,2,J,2)+S1(K,L+1,2)*Q(M,2,J)
          W2(M,1,J,2)=W2(M,1,J,2)+S2(K,L+1,1)*Q(M,2,J)
          W2(M,2,J,2)=W2(M,2,J,2)+S2(K,L+1,2)*Q(M,2,J)
          W3(M,1,J,2)=W3(M,1,J,2)+S3(K,L+1,1)*Q(M,2,J)
          W3(M,2,J,2)=W3(M,2,J,2)+S3(K,L+1,2)*Q(M,2,J)
          W4(M,1,J,2)=W4(M,1,J,2)+S4(K,L+1,1)*Q(M,2,J)
          W4(M,2,J,2)=W4(M,2,J,2)+S4(K,L+1,2)*Q(M,2,J)
          W5(M,1,J,2)=W5(M,1,J,2)+S5(K,L+1,1)*Q(M,2,J)
          W5(M,2,J,2)=W5(M,2,J,2)+S5(K,L+1,2)*Q(M,2,J)
          W6(M,1,J,2)=W6(M,1,J,2)+S6(K,L+1,1)*Q(M,2,J)
          W6(M,2,J,2)=W6(M,2,J,2)+S6(K,L+1,2)*Q(M,2,J)
          W7(M,1,J,2)=W7(M,1,J,2)+S7(K,L+1,1)*Q(M,2,J)
          W7(M,2,J,2)=W7(M,2,J,2)+S7(K,L+1,2)*Q(M,2,J)
          W8(M,1,J,2)=W8(M,1,J,2)+S8(K,L+1,1)*Q(M,2,J)
          W8(M,2,J,2)=W8(M,2,J,2)+S8(K,L+1,2)*Q(M,2,J)
          W9(M,1,J,2)=W9(M,1,J,2)+S9(K,L+1,1)*Q(M,2,J)
          W9(M,2,J,2)=W9(M,2,J,2)+S9(K,L+1,2)*Q(M,2,J)
          WA(M,1,J,2)=WA(M,1,J,2)+SA(K,L+1,1)*Q(M,2,J)
          WA(M,2,J,2)=WA(M,2,J,2)+SA(K,L+1,2)*Q(M,2,J)
          WB(M,1,J,2)=WB(M,1,J,2)+SB(K,L+1,1)*Q(M,2,J)
          WB(M,2,J,2)=WB(M,2,J,2)+SB(K,L+1,2)*Q(M,2,J)
        END DO
      END DO
      DO L=2,MMD-1,2
        DO J=1,JH
          DO K=1,KM*(2*MMP+2)
            M=K+L*KM
            Q(M,1,J)=Q(M,1,J)+Y(J)*R(K+KM,L-1)*Q(M,2,J)
            W1(M,1,J,1)=W1(M,1,J,1)+S1(K+KM,L,1)*Q(M,1,J)
            W1(M,2,J,1)=W1(M,2,J,1)+S1(K+KM,L,2)*Q(M,1,J)
            W2(M,1,J,1)=W2(M,1,J,1)+S2(K+KM,L,1)*Q(M,1,J)
            W2(M,2,J,1)=W2(M,2,J,1)+S2(K+KM,L,2)*Q(M,1,J)
            W3(M,1,J,1)=W3(M,1,J,1)+S3(K+KM,L,1)*Q(M,1,J)
            W3(M,2,J,1)=W3(M,2,J,1)+S3(K+KM,L,2)*Q(M,1,J)
            W4(M,1,J,1)=W4(M,1,J,1)+S4(K+KM,L,1)*Q(M,1,J)
            W4(M,2,J,1)=W4(M,2,J,1)+S4(K+KM,L,2)*Q(M,1,J)
            W5(M,1,J,1)=W5(M,1,J,1)+S5(K+KM,L,1)*Q(M,1,J)
            W5(M,2,J,1)=W5(M,2,J,1)+S5(K+KM,L,2)*Q(M,1,J)
            W6(M,1,J,1)=W6(M,1,J,1)+S6(K+KM,L,1)*Q(M,1,J)
            W6(M,2,J,1)=W6(M,2,J,1)+S6(K+KM,L,2)*Q(M,1,J)
            W7(M,1,J,1)=W7(M,1,J,1)+S7(K+KM,L,1)*Q(M,1,J)
            W7(M,2,J,1)=W7(M,2,J,1)+S7(K+KM,L,2)*Q(M,1,J)
            W8(M,1,J,1)=W8(M,1,J,1)+S8(K+KM,L,1)*Q(M,1,J)
            W8(M,2,J,1)=W8(M,2,J,1)+S8(K+KM,L,2)*Q(M,1,J)
            W9(M,1,J,1)=W9(M,1,J,1)+S9(K+KM,L,1)*Q(M,1,J)
            W9(M,2,J,1)=W9(M,2,J,1)+S9(K+KM,L,2)*Q(M,1,J)
            WA(M,1,J,1)=WA(M,1,J,1)+SA(K+KM,L,1)*Q(M,1,J)
            WA(M,2,J,1)=WA(M,2,J,1)+SA(K+KM,L,2)*Q(M,1,J)
            WB(M,1,J,1)=WB(M,1,J,1)+SB(K+KM,L,1)*Q(M,1,J)
            WB(M,2,J,1)=WB(M,2,J,1)+SB(K+KM,L,2)*Q(M,1,J)
#ifdef SPLIT
          END DO
        END DO
        DO J=1,JH
          DO K=1,KM*(2*MMP+2)
            M=K+L*KM
#endif
            Q(M,2,J)=Q(M,2,J)+Y(J)*R(K,L)*Q(M,1,J)
            W1(M,1,J,2)=W1(M,1,J,2)+S1(K,L+1,1)*Q(M,2,J)
            W1(M,2,J,2)=W1(M,2,J,2)+S1(K,L+1,2)*Q(M,2,J)
            W2(M,1,J,2)=W2(M,1,J,2)+S2(K,L+1,1)*Q(M,2,J)
            W2(M,2,J,2)=W2(M,2,J,2)+S2(K,L+1,2)*Q(M,2,J)
            W3(M,1,J,2)=W3(M,1,J,2)+S3(K,L+1,1)*Q(M,2,J)
            W3(M,2,J,2)=W3(M,2,J,2)+S3(K,L+1,2)*Q(M,2,J)
            W4(M,1,J,2)=W4(M,1,J,2)+S4(K,L+1,1)*Q(M,2,J)
            W4(M,2,J,2)=W4(M,2,J,2)+S4(K,L+1,2)*Q(M,2,J)
            W5(M,1,J,2)=W5(M,1,J,2)+S5(K,L+1,1)*Q(M,2,J)
            W5(M,2,J,2)=W5(M,2,J,2)+S5(K,L+1,2)*Q(M,2,J)
            W6(M,1,J,2)=W6(M,1,J,2)+S6(K,L+1,1)*Q(M,2,J)
            W6(M,2,J,2)=W6(M,2,J,2)+S6(K,L+1,2)*Q(M,2,J)
            W7(M,1,J,2)=W7(M,1,J,2)+S7(K,L+1,1)*Q(M,2,J)
            W7(M,2,J,2)=W7(M,2,J,2)+S7(K,L+1,2)*Q(M,2,J)
            W8(M,1,J,2)=W8(M,1,J,2)+S8(K,L+1,1)*Q(M,2,J)
            W8(M,2,J,2)=W8(M,2,J,2)+S8(K,L+1,2)*Q(M,2,J)
            W9(M,1,J,2)=W9(M,1,J,2)+S9(K,L+1,1)*Q(M,2,J)
            W9(M,2,J,2)=W9(M,2,J,2)+S9(K,L+1,2)*Q(M,2,J)
            WA(M,1,J,2)=WA(M,1,J,2)+SA(K,L+1,1)*Q(M,2,J)
            WA(M,2,J,2)=WA(M,2,J,2)+SA(K,L+1,2)*Q(M,2,J)
            WB(M,1,J,2)=WB(M,1,J,2)+SB(K,L+1,1)*Q(M,2,J)
            WB(M,2,J,2)=WB(M,2,J,2)+SB(K,L+1,2)*Q(M,2,J)
          END DO
        END DO
      END DO
      IF(MOD(MMD,2).EQ.0) THEN
        L=MMD
        DO J=1,JH
          DO K=1,KM*(2*MMP+1)
            M=K+L*KM
            Q(M,1,J)=Q(M,1,J)+Y(J)*R(K+KM,L-1)*Q(M,2,J)
            W1(M,1,J,1)=W1(M,1,J,1)+S1(K+KM,L,1)*Q(M,1,J)
            W1(M,2,J,1)=W1(M,2,J,1)+S1(K+KM,L,2)*Q(M,1,J)
            W2(M,1,J,1)=W2(M,1,J,1)+S2(K+KM,L,1)*Q(M,1,J)
            W2(M,2,J,1)=W2(M,2,J,1)+S2(K+KM,L,2)*Q(M,1,J)
            W3(M,1,J,1)=W3(M,1,J,1)+S3(K+KM,L,1)*Q(M,1,J)
            W3(M,2,J,1)=W3(M,2,J,1)+S3(K+KM,L,2)*Q(M,1,J)
            W4(M,1,J,1)=W4(M,1,J,1)+S4(K+KM,L,1)*Q(M,1,J)
            W4(M,2,J,1)=W4(M,2,J,1)+S4(K+KM,L,2)*Q(M,1,J)
            W5(M,1,J,1)=W5(M,1,J,1)+S5(K+KM,L,1)*Q(M,1,J)
            W5(M,2,J,1)=W5(M,2,J,1)+S5(K+KM,L,2)*Q(M,1,J)
            W6(M,1,J,1)=W6(M,1,J,1)+S6(K+KM,L,1)*Q(M,1,J)
            W6(M,2,J,1)=W6(M,2,J,1)+S6(K+KM,L,2)*Q(M,1,J)
            W7(M,1,J,1)=W7(M,1,J,1)+S7(K+KM,L,1)*Q(M,1,J)
            W7(M,2,J,1)=W7(M,2,J,1)+S7(K+KM,L,2)*Q(M,1,J)
            W8(M,1,J,1)=W8(M,1,J,1)+S8(K+KM,L,1)*Q(M,1,J)
            W8(M,2,J,1)=W8(M,2,J,1)+S8(K+KM,L,2)*Q(M,1,J)
            W9(M,1,J,1)=W9(M,1,J,1)+S9(K+KM,L,1)*Q(M,1,J)
            W9(M,2,J,1)=W9(M,2,J,1)+S9(K+KM,L,2)*Q(M,1,J)
            WA(M,1,J,1)=WA(M,1,J,1)+SA(K+KM,L,1)*Q(M,1,J)
            WA(M,2,J,1)=WA(M,2,J,1)+SA(K+KM,L,2)*Q(M,1,J)
            WB(M,1,J,1)=WB(M,1,J,1)+SB(K+KM,L,1)*Q(M,1,J)
            WB(M,2,J,1)=WB(M,2,J,1)+SB(K+KM,L,2)*Q(M,1,J)
#ifdef SPLIT
          END DO
        END DO
        DO J=1,JH
          DO K=1,KM*(2*MMP+2)
            M=K+L*KM
#endif
            Q(M,2,J)=Q(M,2,J)+Y(J)*R(K,L)*Q(M,1,J)
            W1(M,1,J,2)=W1(M,1,J,2)+S1(K,L+1,1)*Q(M,2,J)
            W1(M,2,J,2)=W1(M,2,J,2)+S1(K,L+1,2)*Q(M,2,J)
            W2(M,1,J,2)=W2(M,1,J,2)+S2(K,L+1,1)*Q(M,2,J)
            W2(M,2,J,2)=W2(M,2,J,2)+S2(K,L+1,2)*Q(M,2,J)
            W3(M,1,J,2)=W3(M,1,J,2)+S3(K,L+1,1)*Q(M,2,J)
            W3(M,2,J,2)=W3(M,2,J,2)+S3(K,L+1,2)*Q(M,2,J)
            W4(M,1,J,2)=W4(M,1,J,2)+S4(K,L+1,1)*Q(M,2,J)
            W4(M,2,J,2)=W4(M,2,J,2)+S4(K,L+1,2)*Q(M,2,J)
            W5(M,1,J,2)=W5(M,1,J,2)+S5(K,L+1,1)*Q(M,2,J)
            W5(M,2,J,2)=W5(M,2,J,2)+S5(K,L+1,2)*Q(M,2,J)
            W6(M,1,J,2)=W6(M,1,J,2)+S6(K,L+1,1)*Q(M,2,J)
            W6(M,2,J,2)=W6(M,2,J,2)+S6(K,L+1,2)*Q(M,2,J)
            W7(M,1,J,2)=W7(M,1,J,2)+S7(K,L+1,1)*Q(M,2,J)
            W7(M,2,J,2)=W7(M,2,J,2)+S7(K,L+1,2)*Q(M,2,J)
            W8(M,1,J,2)=W8(M,1,J,2)+S8(K,L+1,1)*Q(M,2,J)
            W8(M,2,J,2)=W8(M,2,J,2)+S8(K,L+1,2)*Q(M,2,J)
            W9(M,1,J,2)=W9(M,1,J,2)+S9(K,L+1,1)*Q(M,2,J)
            W9(M,2,J,2)=W9(M,2,J,2)+S9(K,L+1,2)*Q(M,2,J)
            WA(M,1,J,2)=WA(M,1,J,2)+SA(K,L+1,1)*Q(M,2,J)
            WA(M,2,J,2)=WA(M,2,J,2)+SA(K,L+1,2)*Q(M,2,J)
            WB(M,1,J,2)=WB(M,1,J,2)+SB(K,L+1,1)*Q(M,2,J)
            WB(M,2,J,2)=WB(M,2,J,2)+SB(K,L+1,2)*Q(M,2,J)
          END DO
        END DO
      ELSE
        L=MMD+1
        DO J=1,JH
          DO K=1,KM*(2*MMP+1)
            M=K+L*KM-KM
            Q(M,1,J)=Q(M,1,J)+Y(J)*R(K,L-1)*Q(M,2,J)
            W1(M,1,J,1)=W1(M,1,J,1)+S1(K,L,1)*Q(M,1,J)
            W1(M,2,J,1)=W1(M,2,J,1)+S1(K,L,2)*Q(M,1,J)
            W2(M,1,J,1)=W2(M,1,J,1)+S2(K,L,1)*Q(M,1,J)
            W2(M,2,J,1)=W2(M,2,J,1)+S2(K,L,2)*Q(M,1,J)
            W3(M,1,J,1)=W3(M,1,J,1)+S3(K,L,1)*Q(M,1,J)
            W3(M,2,J,1)=W3(M,2,J,1)+S3(K,L,2)*Q(M,1,J)
            W4(M,1,J,1)=W4(M,1,J,1)+S4(K,L,1)*Q(M,1,J)
            W4(M,2,J,1)=W4(M,2,J,1)+S4(K,L,2)*Q(M,1,J)
            W5(M,1,J,1)=W5(M,1,J,1)+S5(K,L,1)*Q(M,1,J)
            W5(M,2,J,1)=W5(M,2,J,1)+S5(K,L,2)*Q(M,1,J)
            W6(M,1,J,1)=W6(M,1,J,1)+S6(K,L,1)*Q(M,1,J)
            W6(M,2,J,1)=W6(M,2,J,1)+S6(K,L,2)*Q(M,1,J)
            W7(M,1,J,1)=W7(M,1,J,1)+S7(K,L,1)*Q(M,1,J)
            W7(M,2,J,1)=W7(M,2,J,1)+S7(K,L,2)*Q(M,1,J)
            W8(M,1,J,1)=W8(M,1,J,1)+S8(K,L,1)*Q(M,1,J)
            W8(M,2,J,1)=W8(M,2,J,1)+S8(K,L,2)*Q(M,1,J)
            W9(M,1,J,1)=W9(M,1,J,1)+S9(K,L,1)*Q(M,1,J)
            W9(M,2,J,1)=W9(M,2,J,1)+S9(K,L,2)*Q(M,1,J)
            WA(M,1,J,1)=WA(M,1,J,1)+SA(K,L,1)*Q(M,1,J)
            WA(M,2,J,1)=WA(M,2,J,1)+SA(K,L,2)*Q(M,1,J)
            WB(M,1,J,1)=WB(M,1,J,1)+SB(K,L,1)*Q(M,1,J)
            WB(M,2,J,1)=WB(M,2,J,1)+SB(K,L,2)*Q(M,1,J)
          END DO
        END DO
      END IF

      END
************************************************************************
      SUBROUTINE SPLGS8(MM,JM,KM,
     &  W1,W2,W3,W4,W5,W6,W7,W8,S1,S2,S3,S4,S5,S6,S7,S8,Y,P,R,Q)
*-----------------------------------------------------------------------
* 以下にある ifdef によるループの分解は, 今一つ賢くない SX4 のコンパイラ
* でも外側ループのアンローリングができるようにするためのもので, VPでは,
* 不要(むしろ有害)である.
*-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION W1(KM*((MM+1)/2+MM+1),2,JM/2,2)
      DIMENSION W2(KM*((MM+1)/2+MM+1),2,JM/2,2)
      DIMENSION W3(KM*((MM+1)/2+MM+1),2,JM/2,2)
      DIMENSION W4(KM*((MM+1)/2+MM+1),2,JM/2,2)
      DIMENSION W5(KM*((MM+1)/2+MM+1),2,JM/2,2)
      DIMENSION W6(KM*((MM+1)/2+MM+1),2,JM/2,2)
      DIMENSION W7(KM*((MM+1)/2+MM+1),2,JM/2,2)
      DIMENSION W8(KM*((MM+1)/2+MM+1),2,JM/2,2)
      DIMENSION S1(KM*((MM+1)/2*2+3),0:MM/2+1,2)
      DIMENSION S2(KM*((MM+1)/2*2+3),0:MM/2+1,2)
      DIMENSION S3(KM*((MM+1)/2*2+3),0:MM/2+1,2)
      DIMENSION S4(KM*((MM+1)/2*2+3),0:MM/2+1,2)
      DIMENSION S5(KM*((MM+1)/2*2+3),0:MM/2+1,2)
      DIMENSION S6(KM*((MM+1)/2*2+3),0:MM/2+1,2)      
      DIMENSION S7(KM*((MM+1)/2*2+3),0:MM/2+1,2)
      DIMENSION S8(KM*((MM+1)/2*2+3),0:MM/2+1,2)      
      DIMENSION Y(JM/2)
      DIMENSION P(KM*((MM+1)/2+MM+1),2,JM/2)
      DIMENSION R(KM*((MM+1)/2*2+3),0:MM/2)
      DIMENSION Q(KM*((MM+1)/2+MM+1),2,JM/2)

      MMP=(MM+1)/2
      MMD=MM/2
      JH=JM/2

      CALL BSSET0(KM*((MM+1)/2*2+3)*(MM/2+2)*2,S1)
      CALL BSSET0(KM*((MM+1)/2*2+3)*(MM/2+2)*2,S2)
      CALL BSSET0(KM*((MM+1)/2*2+3)*(MM/2+2)*2,S3)
      CALL BSSET0(KM*((MM+1)/2*2+3)*(MM/2+2)*2,S4)
      CALL BSSET0(KM*((MM+1)/2*2+3)*(MM/2+2)*2,S5)
      CALL BSSET0(KM*((MM+1)/2*2+3)*(MM/2+2)*2,S6)
      CALL BSSET0(KM*((MM+1)/2*2+3)*(MM/2+2)*2,S7)
      CALL BSSET0(KM*((MM+1)/2*2+3)*(MM/2+2)*2,S8)
      CALL BSCOPY(KM*((MM+1)/2+MM+1)*JM,P,Q)

      L=0
      DO J=1,JH
        DO K=1,KM*(2*MMP+2)
          M=K+L*KM
          S1(K+KM,L,1)=S1(K+KM,L,1)+W1(M,1,J,1)*Q(M,1,J)
          S1(K+KM,L,2)=S1(K+KM,L,2)+W1(M,2,J,1)*Q(M,1,J)
          S2(K+KM,L,1)=S2(K+KM,L,1)+W2(M,1,J,1)*Q(M,1,J)
          S2(K+KM,L,2)=S2(K+KM,L,2)+W2(M,2,J,1)*Q(M,1,J)
          S3(K+KM,L,1)=S3(K+KM,L,1)+W3(M,1,J,1)*Q(M,1,J)
          S3(K+KM,L,2)=S3(K+KM,L,2)+W3(M,2,J,1)*Q(M,1,J)
          S4(K+KM,L,1)=S4(K+KM,L,1)+W4(M,1,J,1)*Q(M,1,J)
          S4(K+KM,L,2)=S4(K+KM,L,2)+W4(M,2,J,1)*Q(M,1,J)
          S5(K+KM,L,1)=S5(K+KM,L,1)+W5(M,1,J,1)*Q(M,1,J)
          S5(K+KM,L,2)=S5(K+KM,L,2)+W5(M,2,J,1)*Q(M,1,J)
          S6(K+KM,L,1)=S6(K+KM,L,1)+W6(M,1,J,1)*Q(M,1,J)
          S6(K+KM,L,2)=S6(K+KM,L,2)+W6(M,2,J,1)*Q(M,1,J)
          S7(K+KM,L,1)=S7(K+KM,L,1)+W7(M,1,J,1)*Q(M,1,J)
          S7(K+KM,L,2)=S7(K+KM,L,2)+W7(M,2,J,1)*Q(M,1,J)
          S8(K+KM,L,1)=S8(K+KM,L,1)+W8(M,1,J,1)*Q(M,1,J)
          S8(K+KM,L,2)=S8(K+KM,L,2)+W8(M,2,J,1)*Q(M,1,J)
          S1(K,L+1,1)=S1(K,L+1,1)+W1(M,1,J,2)*Q(M,2,J)
          S1(K,L+1,2)=S1(K,L+1,2)+W1(M,2,J,2)*Q(M,2,J)
          S2(K,L+1,1)=S2(K,L+1,1)+W2(M,1,J,2)*Q(M,2,J)
          S2(K,L+1,2)=S2(K,L+1,2)+W2(M,2,J,2)*Q(M,2,J)
          S3(K,L+1,1)=S3(K,L+1,1)+W3(M,1,J,2)*Q(M,2,J)
          S3(K,L+1,2)=S3(K,L+1,2)+W3(M,2,J,2)*Q(M,2,J)
          S4(K,L+1,1)=S4(K,L+1,1)+W4(M,1,J,2)*Q(M,2,J)
          S4(K,L+1,2)=S4(K,L+1,2)+W4(M,2,J,2)*Q(M,2,J)
          S5(K,L+1,1)=S5(K,L+1,1)+W5(M,1,J,2)*Q(M,2,J)
          S5(K,L+1,2)=S5(K,L+1,2)+W5(M,2,J,2)*Q(M,2,J)
          S6(K,L+1,1)=S6(K,L+1,1)+W6(M,1,J,2)*Q(M,2,J)
          S6(K,L+1,2)=S6(K,L+1,2)+W6(M,2,J,2)*Q(M,2,J)
          S7(K,L+1,1)=S7(K,L+1,1)+W7(M,1,J,2)*Q(M,2,J)
          S7(K,L+1,2)=S7(K,L+1,2)+W7(M,2,J,2)*Q(M,2,J)
          S8(K,L+1,1)=S8(K,L+1,1)+W8(M,1,J,2)*Q(M,2,J)
          S8(K,L+1,2)=S8(K,L+1,2)+W8(M,2,J,2)*Q(M,2,J)
        END DO
      END DO
      DO L=2,MMD-1,2
        DO J=1,JH
          DO K=1,KM*(2*MMP+2)
            M=K+L*KM
            Q(M,1,J)=Q(M,1,J)+Y(J)*R(K+KM,L-1)*Q(M,2,J)
            S1(K+KM,L,1)=S1(K+KM,L,1)+W1(M,1,J,1)*Q(M,1,J)
            S1(K+KM,L,2)=S1(K+KM,L,2)+W1(M,2,J,1)*Q(M,1,J)
            S2(K+KM,L,1)=S2(K+KM,L,1)+W2(M,1,J,1)*Q(M,1,J)
            S2(K+KM,L,2)=S2(K+KM,L,2)+W2(M,2,J,1)*Q(M,1,J)
            S3(K+KM,L,1)=S3(K+KM,L,1)+W3(M,1,J,1)*Q(M,1,J)
            S3(K+KM,L,2)=S3(K+KM,L,2)+W3(M,2,J,1)*Q(M,1,J)
            S4(K+KM,L,1)=S4(K+KM,L,1)+W4(M,1,J,1)*Q(M,1,J)
            S4(K+KM,L,2)=S4(K+KM,L,2)+W4(M,2,J,1)*Q(M,1,J)
            S5(K+KM,L,1)=S5(K+KM,L,1)+W5(M,1,J,1)*Q(M,1,J)
            S5(K+KM,L,2)=S5(K+KM,L,2)+W5(M,2,J,1)*Q(M,1,J)
            S6(K+KM,L,1)=S6(K+KM,L,1)+W6(M,1,J,1)*Q(M,1,J)
            S6(K+KM,L,2)=S6(K+KM,L,2)+W6(M,2,J,1)*Q(M,1,J)
            S7(K+KM,L,1)=S7(K+KM,L,1)+W7(M,1,J,1)*Q(M,1,J)
            S7(K+KM,L,2)=S7(K+KM,L,2)+W7(M,2,J,1)*Q(M,1,J)
            S8(K+KM,L,1)=S8(K+KM,L,1)+W8(M,1,J,1)*Q(M,1,J)
            S8(K+KM,L,2)=S8(K+KM,L,2)+W8(M,2,J,1)*Q(M,1,J)
#ifdef SPLIT
          END DO
        END DO
        DO J=1,JH
          DO K=1,KM*(2*MMP+2)
            M=K+L*KM
#endif
            Q(M,2,J)=Q(M,2,J)+Y(J)*R(K,L)*Q(M,1,J)
            S1(K,L+1,1)=S1(K,L+1,1)+W1(M,1,J,2)*Q(M,2,J)
            S1(K,L+1,2)=S1(K,L+1,2)+W1(M,2,J,2)*Q(M,2,J)
            S2(K,L+1,1)=S2(K,L+1,1)+W2(M,1,J,2)*Q(M,2,J)
            S2(K,L+1,2)=S2(K,L+1,2)+W2(M,2,J,2)*Q(M,2,J)
            S3(K,L+1,1)=S3(K,L+1,1)+W3(M,1,J,2)*Q(M,2,J)
            S3(K,L+1,2)=S3(K,L+1,2)+W3(M,2,J,2)*Q(M,2,J)
            S4(K,L+1,1)=S4(K,L+1,1)+W4(M,1,J,2)*Q(M,2,J)
            S4(K,L+1,2)=S4(K,L+1,2)+W4(M,2,J,2)*Q(M,2,J)
            S5(K,L+1,1)=S5(K,L+1,1)+W5(M,1,J,2)*Q(M,2,J)
            S5(K,L+1,2)=S5(K,L+1,2)+W5(M,2,J,2)*Q(M,2,J)
            S6(K,L+1,1)=S6(K,L+1,1)+W6(M,1,J,2)*Q(M,2,J)
            S6(K,L+1,2)=S6(K,L+1,2)+W6(M,2,J,2)*Q(M,2,J)
            S7(K,L+1,1)=S7(K,L+1,1)+W7(M,1,J,2)*Q(M,2,J)
            S7(K,L+1,2)=S7(K,L+1,2)+W7(M,2,J,2)*Q(M,2,J)
            S8(K,L+1,1)=S8(K,L+1,1)+W8(M,1,J,2)*Q(M,2,J)
            S8(K,L+1,2)=S8(K,L+1,2)+W8(M,2,J,2)*Q(M,2,J)
          END DO
        END DO
      END DO
      IF(MOD(MMD,2).EQ.0) THEN
        L=MMD
        DO J=1,JH
          DO K=1,KM*(2*MMP+1)
            M=K+L*KM
            Q(M,1,J)=Q(M,1,J)+Y(J)*R(K+KM,L-1)*Q(M,2,J)
            S1(K+KM,L,1)=S1(K+KM,L,1)+W1(M,1,J,1)*Q(M,1,J)
            S1(K+KM,L,2)=S1(K+KM,L,2)+W1(M,2,J,1)*Q(M,1,J)
            S2(K+KM,L,1)=S2(K+KM,L,1)+W2(M,1,J,1)*Q(M,1,J)
            S2(K+KM,L,2)=S2(K+KM,L,2)+W2(M,2,J,1)*Q(M,1,J)
            S3(K+KM,L,1)=S3(K+KM,L,1)+W3(M,1,J,1)*Q(M,1,J)
            S3(K+KM,L,2)=S3(K+KM,L,2)+W3(M,2,J,1)*Q(M,1,J)
            S4(K+KM,L,1)=S4(K+KM,L,1)+W4(M,1,J,1)*Q(M,1,J)
            S4(K+KM,L,2)=S4(K+KM,L,2)+W4(M,2,J,1)*Q(M,1,J)
            S5(K+KM,L,1)=S5(K+KM,L,1)+W5(M,1,J,1)*Q(M,1,J)
            S5(K+KM,L,2)=S5(K+KM,L,2)+W5(M,2,J,1)*Q(M,1,J)
            S6(K+KM,L,1)=S6(K+KM,L,1)+W6(M,1,J,1)*Q(M,1,J)
            S6(K+KM,L,2)=S6(K+KM,L,2)+W6(M,2,J,1)*Q(M,1,J)
            S7(K+KM,L,1)=S7(K+KM,L,1)+W7(M,1,J,1)*Q(M,1,J)
            S7(K+KM,L,2)=S7(K+KM,L,2)+W7(M,2,J,1)*Q(M,1,J)
            S8(K+KM,L,1)=S8(K+KM,L,1)+W8(M,1,J,1)*Q(M,1,J)
            S8(K+KM,L,2)=S8(K+KM,L,2)+W8(M,2,J,1)*Q(M,1,J)
#ifdef SPLIT
          END DO
        END DO
        DO J=1,JH
          DO K=1,KM*(2*MMP+2)
            M=K+L*KM
#endif
            Q(M,2,J)=Q(M,2,J)+Y(J)*R(K,L)*Q(M,1,J)
            S1(K,L+1,1)=S1(K,L+1,1)+W1(M,1,J,2)*Q(M,2,J)
            S1(K,L+1,2)=S1(K,L+1,2)+W1(M,2,J,2)*Q(M,2,J)
            S2(K,L+1,1)=S2(K,L+1,1)+W2(M,1,J,2)*Q(M,2,J)
            S2(K,L+1,2)=S2(K,L+1,2)+W2(M,2,J,2)*Q(M,2,J)
            S3(K,L+1,1)=S3(K,L+1,1)+W3(M,1,J,2)*Q(M,2,J)
            S3(K,L+1,2)=S3(K,L+1,2)+W3(M,2,J,2)*Q(M,2,J)
            S4(K,L+1,1)=S4(K,L+1,1)+W4(M,1,J,2)*Q(M,2,J)
            S4(K,L+1,2)=S4(K,L+1,2)+W4(M,2,J,2)*Q(M,2,J)
            S5(K,L+1,1)=S5(K,L+1,1)+W5(M,1,J,2)*Q(M,2,J)
            S5(K,L+1,2)=S5(K,L+1,2)+W5(M,2,J,2)*Q(M,2,J)
            S6(K,L+1,1)=S6(K,L+1,1)+W6(M,1,J,2)*Q(M,2,J)
            S6(K,L+1,2)=S6(K,L+1,2)+W6(M,2,J,2)*Q(M,2,J)
            S7(K,L+1,1)=S7(K,L+1,1)+W7(M,1,J,2)*Q(M,2,J)
            S7(K,L+1,2)=S7(K,L+1,2)+W7(M,2,J,2)*Q(M,2,J)
            S8(K,L+1,1)=S8(K,L+1,1)+W8(M,1,J,2)*Q(M,2,J)
            S8(K,L+1,2)=S8(K,L+1,2)+W8(M,2,J,2)*Q(M,2,J)
          END DO
        END DO
      ELSE
        L=MMD+1
        DO J=1,JH
          DO K=1,KM*(2*MMP+1)
            M=K+L*KM-KM
            Q(M,1,J)=Q(M,1,J)+Y(J)*R(K,L-1)*Q(M,2,J)
            S1(K,L,1)=S1(K,L,1)+W1(M,1,J,1)*Q(M,1,J)
            S1(K,L,2)=S1(K,L,2)+W1(M,2,J,1)*Q(M,1,J)
            S2(K,L,1)=S2(K,L,1)+W2(M,1,J,1)*Q(M,1,J)
            S2(K,L,2)=S2(K,L,2)+W2(M,2,J,1)*Q(M,1,J)
            S3(K,L,1)=S3(K,L,1)+W3(M,1,J,1)*Q(M,1,J)
            S3(K,L,2)=S3(K,L,2)+W3(M,2,J,1)*Q(M,1,J)
            S4(K,L,1)=S4(K,L,1)+W4(M,1,J,1)*Q(M,1,J)
            S4(K,L,2)=S4(K,L,2)+W4(M,2,J,1)*Q(M,1,J)
            S5(K,L,1)=S5(K,L,1)+W5(M,1,J,1)*Q(M,1,J)
            S5(K,L,2)=S5(K,L,2)+W5(M,2,J,1)*Q(M,1,J)
            S6(K,L,1)=S6(K,L,1)+W6(M,1,J,1)*Q(M,1,J)
            S6(K,L,2)=S6(K,L,2)+W6(M,2,J,1)*Q(M,1,J)
            S7(K,L,1)=S7(K,L,1)+W7(M,1,J,1)*Q(M,1,J)
            S7(K,L,2)=S7(K,L,2)+W7(M,2,J,1)*Q(M,1,J)
            S8(K,L,1)=S8(K,L,1)+W8(M,1,J,1)*Q(M,1,J)
            S8(K,L,2)=S8(K,L,2)+W8(M,2,J,1)*Q(M,1,J)
          END DO
        END DO
      END IF

      END
