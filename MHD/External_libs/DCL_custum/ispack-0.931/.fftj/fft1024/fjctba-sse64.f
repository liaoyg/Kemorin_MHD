************************************************************************
* FTTJ:  An FFT library
* Copyright (C) 2008 Keiichi Ishioka <ishioka@gfd-dennou.org>
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
*     initialize table for fft1024 (for sse64) (backward)
*-----------------------------------------------------------------------
      SUBROUTINE FJCTBA(XT)

      IMPLICIT REAL*8(A-H,O-Y),COMPLEX*16(Z)
      PARAMETER(N=1024)
      DIMENSION XT(2,2,0:1023)

      PI=4*ATAN(1D0)

      DO J=0,3
        DO I=0,15
          XT(1,1,I+J*16)=COS(2*PI*I*J/64)
          XT(2,1,I+J*16)=COS(2*PI*I*J/64)
          XT(1,2,I+J*16)=-SIN(2*PI*I*J/64)
          XT(2,2,I+J*16)=SIN(2*PI*I*J/64)        
        END DO
      END DO

      DO J=0,63
        DO K=1,15
          XT(1,1,K+J*15+63)=COS(2*PI*J*K/1024)
          XT(2,1,K+J*15+63)=COS(2*PI*J*K/1024)
          XT(1,2,K+J*15+63)=-SIN(2*PI*J*K/1024)
          XT(2,2,K+J*15+63)=SIN(2*PI*J*K/1024)        
        END DO
      END DO

      END
