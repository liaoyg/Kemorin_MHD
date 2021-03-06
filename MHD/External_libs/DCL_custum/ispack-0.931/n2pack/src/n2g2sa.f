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
*     TRANSFORM GRID TO SPECTRA                                 96/10/26
************************************************************************
      SUBROUTINE N2G2SA(LM,KM,JM,IM,G,S,W,ITJ,TJ,ITI,TI)

      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION G(0:JM-1,2,0:IM/2-1)
      DIMENSION S(-LM:LM,-KM:KM)
      DIMENSION W(JM*IM)
      DIMENSION ITJ(5),TJ(JM*2),ITI(5),TI(IM*2)

      CALL FTTRUF(JM,IM,G,W,ITI,TI)

      CALL BSSET0(JM,G(0,2,0))

      DO K=0,KM
        CALL FTTZUF(1,JM,G(0,1,K),W,ITJ,TJ)
      END DO

      DO K=1,KM
        DO L=1,LM
          S( L, K)=G(L,1,K)
          S(-L,-K)=G(L,2,K)
          S(-L, K)=G(JM-L,1,K)
          S( L,-K)=G(JM-L,2,K)
        END DO
      END DO

      DO K=1,KM
        S(0, K)=G(0,1,K)
        S(0,-K)=G(0,2,K)
      END DO

      DO L=1,LM
        S( L,0)=G(L,1,0)
        S(-L,0)=G(L,2,0)
      END DO

      S(0,0)=G(0,1,0)

      END
