########################################################################
# FTTJ:  An FFT library
# Copyright (C) 2008 Keiichi Ishioka <ishioka@gfd-dennou.org>
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301 USA.
########################################################################
# rfft 16 out-of-place forward
.text
.globl fjrof4_
fjrof4_:
       movl    4(%esp), %eax
       movl    8(%esp), %ecx

#-------------------------------------------------------------
# fft 8 in-place backward
#.text
#.globl fjcib3_
#fjcib3_:
#      movl	4(%esp), %eax

      movaps    (%eax), %xmm6 # Z(0)       
      movaps  64(%eax), %xmm7 # Z(4)
      movaps  16(%eax), %xmm0 # Z(1)
      movaps  80(%eax), %xmm1 # Z(5)              
      movaps  48(%eax), %xmm2 # Z(3)       
      movaps 112(%eax), %xmm3 # Z(7)
       
      movaps %xmm0,%xmm4
      subpd %xmm1,%xmm0
      addpd %xmm4,%xmm1
       
      movaps %xmm2,%xmm4
      subpd %xmm3,%xmm2
      addpd %xmm4,%xmm3
       
      movaps %xmm1,%xmm4
      subpd %xmm3,%xmm1
      addpd %xmm4,%xmm3
       
      movaps %xmm0,%xmm4
      subpd %xmm2,%xmm0
      addpd %xmm4,%xmm2
       
      movaps %xmm6,%xmm4
      subpd %xmm7,%xmm6
      addpd %xmm4,%xmm7
       
      movaps 32(%eax), %xmm4 # Z(2)
      movaps 96(%eax), %xmm5 # Z(6)              
      movaps %xmm6, (%ecx)

      movaps %xmm4,%xmm6
      subpd %xmm5,%xmm4
      addpd %xmm6,%xmm5
       
      movaps .CI,%xmm6

      shufpd $0x1,%xmm2,%xmm2
      mulpd .C2D,%xmm2
      shufpd $0x1,%xmm4,%xmm4
      xorpd %xmm6,%xmm4
      shufpd $0x1,%xmm1,%xmm1
      xorpd %xmm6,%xmm1
      mulpd .C2,%xmm0              

      movaps %xmm2,%xmm6
      subpd %xmm4,%xmm2
      addpd %xmm6,%xmm4
       
      movaps %xmm7,%xmm6
      subpd %xmm5,%xmm7
      addpd %xmm6,%xmm5
       
      movaps %xmm5,%xmm6
      subpd %xmm3,%xmm5
      addpd %xmm6,%xmm3

      movaps (%ecx),%xmm6
      movaps %xmm5,64(%ecx)
      movaps %xmm3,  (%ecx)

      movaps %xmm7,%xmm5
      subpd %xmm1,%xmm7
      addpd %xmm5,%xmm1

      movaps %xmm7, 96(%ecx)
      movaps %xmm1, 32(%ecx)
       
      movaps %xmm6,%xmm3
      subpd %xmm0,%xmm6
      addpd %xmm3,%xmm0
       
      movaps %xmm6,%xmm3
      subpd %xmm2,%xmm6
      addpd %xmm3,%xmm2
       
      movaps %xmm6, 80(%ecx)       
      movaps %xmm2, 48(%ecx)

      movaps %xmm0,%xmm2
      subpd %xmm4,%xmm0
      addpd %xmm2,%xmm4
       
      movaps %xmm0,112(%ecx)
      movaps %xmm4, 16(%ecx)
#-------------------------------------------------------------       
#  rfft のための後処理

      movaps .CID,%xmm7
      
      movaps (%ecx),%xmm0
      movaps %xmm0,%xmm1            
      xorpd %xmm7,%xmm0      
      shufpd $1,%xmm1,%xmm1      
      addpd %xmm1,%xmm0
      movapd %xmm0,(%ecx)
      
      movaps  16(%ecx),%xmm0
      movaps 112(%ecx),%xmm1      
      xorpd %xmm7,%xmm0
      movaps %xmm1,%xmm2      
      subpd %xmm0,%xmm1
      
      movaps %xmm1,%xmm3
      shufpd $1,%xmm1,%xmm1
      mulpd .A1,%xmm3
      mulpd .B1,%xmm1
      addpd %xmm3,%xmm1
      
      addpd %xmm1,%xmm0
      subpd %xmm1,%xmm2
      xorpd %xmm7,%xmm2
      movaps %xmm0, 16(%ecx)      
      movaps %xmm2,112(%ecx)
      
      movaps 32(%ecx),%xmm0
      movaps 96(%ecx),%xmm1      
      xorpd %xmm7,%xmm0
      movaps %xmm1,%xmm2      
      subpd %xmm0,%xmm1
      
      movaps %xmm1,%xmm3
      shufpd $1,%xmm1,%xmm1
      mulpd .A2,%xmm3
      mulpd .B2,%xmm1
      addpd %xmm3,%xmm1
      
      addpd %xmm1,%xmm0
      subpd %xmm1,%xmm2
      xorpd %xmm7,%xmm2
      movaps %xmm0,32(%ecx)      
      movaps %xmm2,96(%ecx)

      movaps 48(%ecx),%xmm0
      movaps 80(%ecx),%xmm1      
      xorpd %xmm7,%xmm0
      movaps %xmm1,%xmm2      
      subpd %xmm0,%xmm1
      
      movaps %xmm1,%xmm3
      shufpd $1,%xmm1,%xmm1
      mulpd .A3,%xmm3
      mulpd .B3,%xmm1
      addpd %xmm3,%xmm1
      
      addpd %xmm1,%xmm0
      subpd %xmm1,%xmm2
      xorpd %xmm7,%xmm2
      movaps %xmm0,48(%ecx)      
      movaps %xmm2,80(%ecx)

      movaps 64(%ecx),%xmm0
      xorpd %xmm7,%xmm0
      movaps %xmm0,64(%ecx)

#--------------------------------------------------       
.LE:

       ret

#--------------------------------------------------
.section .rodata
.align 16
.CID:
      .long 0x0, 0x0, 0x0, 0x80000000
.align 16
.CI:
      .long 0x0, 0x80000000, 0x0, 0x0
.align 16
.C2:
      .long 0x667f3bcd, 0x3fe6a09e,0x667f3bcd, 0x3fe6a09e
.align 16
.C2D:
      .long 0x667f3bcd, 0xbfe6a09e,0x667f3bcd, 0x3fe6a09e
.align 16
.A1:
      .long 0xaca8ab4e,0x3fd3c10e,0xaca8ab4e,0x3fd3c10e
.align 16
.B1:
      .long 0xcf328d46,0x3fdd906b,0xcf328d46,0xbfdd906b
.align 16
.A2:
      .long 0x33018866,0x3fc2bec3,0x33018866,0x3fc2bec3
.align 16
.B2:
      .long 0x667f3bcd,0x3fd6a09e,0x667f3bcd,0xbfd6a09e
.align 16
.A3:
      .long 0x866b95d0,0x3fa37ca1,0x866b95d0,0x3fa37ca1
.align 16
.B3:
      .long 0xa6aea963,0x3fc87de2,0xa6aea963,0xbfc87de2
