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
# fft 4 in-place forward
.text
.globl fjcif2_
fjcif2_:
      movl 4(%esp),%eax
      movaps .CI,%xmm7

      movaps   (%eax),%xmm0
      movaps 16(%eax),%xmm1
      movaps 32(%eax),%xmm2
      movaps 48(%eax),%xmm3       

      movaps %xmm0,%xmm4
      subpd  %xmm2,%xmm0
      addpd  %xmm4,%xmm2

      movaps %xmm1,%xmm4
      subpd  %xmm3,%xmm1
      addpd  %xmm4,%xmm3

      xorpd  %xmm7,%xmm1
      shufpd $0x1,%xmm1,%xmm1
       
      movaps %xmm2,%xmm4
      addpd  %xmm3,%xmm2
      subpd  %xmm3,%xmm4

      movaps %xmm0,%xmm5
      addpd  %xmm1,%xmm0
      subpd  %xmm1,%xmm5
       
      movaps %xmm2,  (%eax)
      movaps %xmm0,16(%eax)       
      movaps %xmm4,32(%eax)
      movaps %xmm5,48(%eax)
      ret

.align 16
.CI:
      .long  0x0, 0x80000000, 0x0, 0x0
