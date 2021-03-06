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
*     CALCULATE ENERGY SPECTRUM FOR 3D EULER EQ.              2002/05/02
************************************************************************
      SUBROUTINE P3EMPT(NM,MM,LM,KMAX,Z,ES,W)

      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mpif.h'
      DIMENSION Z(-NM:NM,-MM:MM,2,0:*)
      DIMENSION ES(KMAX),W(KMAX)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,IP,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NP,IERR)

      LP=LM/NP+1
      LS=LP*IP
      LE=MIN(LP*(IP+1)-1,LM)
      IF(LE.GE.LS) THEN
        LC=LE-LS+1
      ELSE
        LC=0
        LS=0
        LE=0
      END IF

      LT=2*LC-1+LS

      CALL BSSET0(KMAX,W)

      IF(LC.GT.0) THEN
        DO L=MAX(LS,1),LE
          DO M=-MM,MM
            DO N=-NM,NM
              K=SQRT(1D0*(L*L+M*M+N*N))+0.5D0
              IF(K.LE.KMAX) THEN
                W(K)=W(K)+((L*L+M*M)*Z(N,M,1,L-LS)*Z(N,M,1,L-LS)
     &            +(L*L+N*N)*Z(N,M,2,L-LS)*Z(N,M,2,L-LS)
     &            +2*M*N*Z(N,M,1,L-LS)*Z(N,M,2,L-LS))
     &            /(1D0*(L*L+M*M+N*N)*L*L)
                W(K)=W(K)+((L*L+M*M)*Z(N,M,1,LT-L)*Z(N,M,1,LT-L)
     &            +(L*L+N*N)*Z(N,M,2,LT-L)*Z(N,M,2,LT-L)
     &            +2*M*N*Z(N,M,1,LT-L)*Z(N,M,2,LT-L))
     &            /(1D0*(L*L+M*M+N*N)*L*L)
              END IF              
            END DO
          END DO
        END DO

        IF(LS.EQ.0) THEN
          L=0
          DO M=-MM,-1
            DO N=-NM,NM
              K=SQRT(1D0*(M*M+N*N))+0.5D0
              IF(K.LE.KMAX) THEN
                W(K)=W(K)+(M*M*Z(N,M,2,0)*Z(N,M,2,0)
     &            +(M*M+N*N)*Z(N,M,1,0)*Z(N,M,1,0))
     &            /(1D0*M*M*(M*M+N*N))
              END IF
            END DO
          END DO
          DO M=1,MM
            DO N=-NM,NM
              K=SQRT(1D0*(M*M+N*N))+0.5D0
              IF(K.LE.KMAX) THEN
                W(K)=W(K)+(M*M*Z(N,M,2,0)*Z(N,M,2,0)
     &            +(M*M+N*N)*Z(N,M,1,0)*Z(N,M,1,0))
     &            /(1D0*M*M*(M*M+N*N))
              END IF
            END DO
          END DO
      
          L=0
          M=0
          DO N=-NM,-1
            K=-N
            IF(K.LE.KMAX) THEN
              W(K)=W(K)
     &          +(Z(N,M,1,0)*Z(N,M,1,0)+Z(N,M,2,0)*Z(N,M,2,0))/(N*N)
            END IF
          END DO
          DO N=1,NM
            K=N
            IF(K.LE.KMAX) THEN
              W(K)=W(K)
     &          +(Z(N,M,1,0)*Z(N,M,1,0)+Z(N,M,2,0)*Z(N,M,2,0))/(N*N)
            END IF
          END DO
        END IF

      END IF

      CALL BSSET0(KMAX,ES)      

      CALL MPI_ALLREDUCE(W,ES,KMAX,
     &  MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)      

      END
