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
*     TRANSFORM GRID TO SPECTRA                               2009/08/15
*     (OpenMP 並列化版)
*-----------------------------------------------------------------------
*     W(IM*(JM+1)), Q(JM/2*7*NT), WS((NN+1)*2*NT) (NTは最大スレッド数)
************************************************************************
      SUBROUTINE SJTGOS(MM,NM,NN,IM,JM,S,G,IT,T,P,Q,R,WS,WG,W,IPOW)

      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION S((2*NN+1-MM)*MM+NN+1),G(0:IM-1,JM)
      DIMENSION IT(2,2),T(IM*3,2)      
      DIMENSION P(JM/2,MM+4),Q(JM/2*7,0:MM),R(*)
      DIMENSION WS((NN+1)*2,0:*),WG(0:IM+1,JM),W(0:*)
!$    INTEGER omp_get_thread_num

      IP=0      

!$omp parallel do
      DO J=1,JM
        CALL FJRRUN(G(0,J),WG(0,J),W(IM*(J-1)),T(1,2),IT(1,2))
      END DO
!$omp end parallel do

      F=1D0/IM      

!$omp parallel do
      DO M=0,MM
        IF(M.EQ.0) THEN
          DO J=1,JM
            W(J-1)=WG(0,J)*F
          END DO
        ELSE
          DO J=1,JM        
            W((JM+1)*(2*M-1)+2*(J-1))=WG(2*M,J)*F
            W((JM+1)*(2*M-1)+2*(J-1)+1)=WG(2*M+1,J)*F
          END DO
        END IF
      END DO
!$omp end parallel do            

!$omp parallel do private(NS,M,IP)
      DO MD=0,MM/2
!$      IP=omp_get_thread_num()        
        M=MD
        IF(M.EQ.0) THEN
          CALL LJTGZS(NM,NN,JM,S,W(0),P,Q(1,IP),R,WS(1,IP),IPOW)
        ELSE
          NS=1+NN+1+(M-1)*(NN+NN+2-M)
          CALL LJTGWS(NM,NN,JM,M,S(NS),W((JM+1)*(2*M-1)),P,Q(1,IP),R,
     &      WS(1,IP),IPOW)
        END IF
        M=MM-MD
        IF(2*M.NE.MM) THEN
          NS=1+NN+1+(M-1)*(NN+NN+2-M)
          CALL LJTGWS(NM,NN,JM,M,S(NS),W((JM+1)*(2*M-1)),P,Q(1,IP),R,
     &      WS(1,IP),IPOW)
        END IF
      END DO
!$omp end parallel do      

      END
