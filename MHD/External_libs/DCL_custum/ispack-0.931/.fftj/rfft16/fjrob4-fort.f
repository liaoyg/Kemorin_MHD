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
* rfft 16 out-of-place backward
*-----------------------------------------------------------------------
      SUBROUTINE FJROB4(Z,ZD,ZDD,ZT)

      IMPLICIT REAL*8(A-H,O-Y),COMPLEX*16(Z)
      PARAMETER(C1P8=0.7071067811865475244D0) ! COS(2*PI/8)
      PARAMETER(ZT1=(-0.3826834323650897717D0,0.9238795325112867561D0))
      PARAMETER(ZT2=(-0.7071067811865475244D0,0.7071067811865475244D0))
      PARAMETER(ZT3=(-0.9238795325112867561D0,0.3826834323650897717D0))
      PARAMETER(ZI=(0,1))
      DIMENSION Z(0:7),ZD(0:7),ZDD(*)
      DIMENSION ZT(*)
      
*-----------------------------------------------------------------------
* rfft のための前処理
*-----------------------------------------------------------------------

      ZD(0)=(1+ZI)*CONJG(Z(0))
      
      ZTMP1=(Z(1)-CONJG(Z(7)))*ZT1
      ZTMP0=Z(1)+CONJG(Z(7))
      ZD(1)=ZTMP0+ZTMP1
      ZD(7)=CONJG(ZTMP0-ZTMP1)

      ZTMP1=(Z(2)-CONJG(Z(6)))*ZT2
      ZTMP0=Z(2)+CONJG(Z(6))
      ZD(2)=ZTMP0+ZTMP1
      ZD(6)=CONJG(ZTMP0-ZTMP1)

      ZTMP1=(Z(3)-CONJG(Z(5)))*ZT3
      ZTMP0=Z(3)+CONJG(Z(5))
      ZD(3)=ZTMP0+ZTMP1
      ZD(5)=CONJG(ZTMP0-ZTMP1)

      ZD(4)=2*CONJG(Z(4))

*-----------------------------------------------------------------------
* fft 8 in-place backward
*-----------------------------------------------------------------------
      Z0=ZD(1)
      Z1=ZD(5)
      Z2=ZD(3)
      Z3=ZD(7)
      
      Z7=Z0
      Z0=Z0-Z1    ! Z(1)-Z(5)
      Z1=Z1+Z7    ! Z(1)+Z(5)

      Z7=Z2
      Z2=Z2-Z3    ! Z(3)-Z(7)
      Z3=Z3+Z7    ! Z(3)+Z(7)

      Z7=Z1
      Z1=Z1-Z3    ! (Z(1)+Z(5))-(Z(3)+Z(7))
      Z3=Z3+Z7    ! (Z(1)+Z(5))+(Z(3)+Z(7))

      Z7=Z0
      Z0=Z0-Z2    ! (Z(1)-Z(5))-(Z(3)-Z(7))
      Z2=Z2+Z7    ! (Z(1)-Z(5))+(Z(3)-Z(7))

      Z1=Z1*ZI      ! i*((Z(1)+Z(5))-(Z(3)+Z(7)))
      Z2=Z2*ZI*C1P8 ! i*C1P8((Z(1)-Z(5))+(Z(3)-Z(7)))
      Z0=Z0*C1P8    ! C1P8*((Z(1)-Z(5))-(Z(3)-Z(7)))

      Z4=ZD(2)
      Z5=ZD(6)

      Z7=Z4
      Z4=Z4-Z5 ! Z(2)-Z(6)
      Z5=Z5+Z7 ! Z(2)+Z(6)

      Z4=Z4*ZI ! i*(Z(2)-Z(6))

      Z7=Z2
      Z2=Z2-Z4 ! -i*(Z(2)-Z(6))+i*C1P8((Z(1)-Z(5))+(Z(3)-Z(7)))
      Z4=Z4+Z7 !  i*(Z(2)-Z(6))+i*C1P8((Z(1)-Z(5))+(Z(3)-Z(7)))

      ZD(2)=Z3 ! 一時退避(Z3空き)

      Z6=ZD(0)
      Z7=ZD(4)
      Z3=Z6
      Z6=Z6-Z7 ! Z(0)-Z(4)
      Z7=Z7+Z3 ! Z(0)+Z(4)

      Z3=Z6
      Z6=Z6-Z0 ! (Z(0)-Z(4))-C1P8*((Z(1)-Z(5))-(Z(3)-Z(7)))
      Z0=Z0+Z3 ! (Z(0)-Z(4))+C1P8*((Z(1)-Z(5))-(Z(3)-Z(7)))

      Z3=Z6
      Z6=Z6-Z2
      Z2=Z2+Z3
      ZD(3)=Z2
      ZD(5)=Z6

      Z2=Z0
      Z0=Z0-Z4
      Z4=Z4+Z2
      ZD(7)=Z0
      ZD(1)=Z4


      Z2=Z7
      Z7=Z7-Z5 ! (Z(0)+Z(4))-(Z(2)+Z(6))
      Z5=Z5+Z2 ! (Z(0)+Z(4))+(Z(2)+Z(6))

      Z3=ZD(2) ! Z3の復元

      Z2=Z5
      Z5=Z5-Z3 ! (Z(0)+Z(4))+(Z(2)+Z(6))-((Z(1)+Z(5))+(Z(3)+Z(7)))
      Z3=Z3+Z2 ! (Z(0)+Z(4))+(Z(2)+Z(6))+((Z(1)+Z(5))+(Z(3)+Z(7)))
      ZD(4)=Z5
      ZD(0)=Z3

      Z2=Z7
      Z7=Z7-Z1 ! (Z(0)+Z(4))-(Z(2)+Z(6))-i((Z(1)+Z(5))-(Z(3)+Z(7)))
      Z1=Z1+Z2 ! (Z(0)+Z(4))-(Z(2)+Z(6))+i((Z(1)+Z(5))-(Z(3)+Z(7)))
      ZD(2)=Z1
      ZD(6)=Z7

      END
