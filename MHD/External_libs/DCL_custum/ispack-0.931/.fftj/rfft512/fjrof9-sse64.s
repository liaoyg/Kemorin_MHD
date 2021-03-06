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
# rfft 512 out-of-place forward
.text
.globl fjrof9_
fjrof9_:
       movq %rsi,%r11
       movq %rcx,%r10
#-------------------------------------------------------------              
# fft 256 out-of-place backward (16x16 type)
#.text
#.globl fjcob8_
#fjcob8_:

#------------------------------------
# fftj16 のループ(前処理の方)
#------------------------------------

       movq %rdi,%r8
       addq $256,%r8
       
       movaps  256(%rdi), %xmm0
       movaps 1280(%rdi), %xmm2       
       movaps 2304(%rdi), %xmm1
       movaps 3328(%rdi), %xmm3

       movaps %xmm0,%xmm7
       subpd %xmm1,%xmm0
       addpd %xmm7,%xmm1

       movaps %xmm2,%xmm7
       subpd %xmm3,%xmm2
       addpd %xmm7,%xmm3

       movaps %xmm1,%xmm7
       subpd %xmm3,%xmm1
       addpd %xmm7,%xmm3

       movaps %xmm3,16(%rsi)
       
#-------------

       movaps  768(%rdi), %xmm5
       movaps 1792(%rdi), %xmm6       
       movaps 2816(%rdi), %xmm3
       movaps 3840(%rdi), %xmm4

       movaps %xmm5,%xmm7
       subpd %xmm3,%xmm5
       addpd %xmm7,%xmm3

       movaps %xmm6,%xmm7
       subpd %xmm4,%xmm6
       addpd %xmm7,%xmm4

       movaps %xmm3,%xmm7
       subpd %xmm4,%xmm3
       addpd %xmm7,%xmm4

#-------------

       movaps .C1,%xmm8
       movaps .C3,%xmm9       

       movaps %xmm0,%xmm7
       subpd %xmm6,%xmm0
       addpd %xmm7,%xmm6

       movaps %xmm5,%xmm7
       subpd %xmm2,%xmm5
       addpd %xmm7,%xmm2

       movaps %xmm8,%xmm7
       movaps %xmm9,%xmm10
       
       mulpd %xmm5,%xmm7
       mulpd %xmm0,%xmm10       
       subpd %xmm7,%xmm10
       
       mulpd %xmm8,%xmm0
       mulpd %xmm9,%xmm5
       addpd %xmm5,%xmm0
       
       movaps %xmm8,%xmm7
       
       mulpd %xmm6,%xmm7
       mulpd %xmm9,%xmm6
       
       mulpd %xmm2,%xmm9
       subpd %xmm9,%xmm7

       mulpd %xmm8,%xmm2
       addpd %xmm2,%xmm6

#-----------------------------------

       movaps     (%rdi), %xmm11
       movaps 2048(%rdi), %xmm9
       movaps 1024(%rdi), %xmm2
       movaps 3072(%rdi), %xmm8
       
       movaps %xmm11,%xmm14
       subpd %xmm9,%xmm11       
       addpd %xmm14,%xmm9

       movaps %xmm2,%xmm14
       subpd %xmm8,%xmm2
       addpd %xmm14,%xmm8
       
       movaps %xmm9,%xmm14
       subpd %xmm8,%xmm9
       addpd %xmm14,%xmm8

       movaps %xmm8,(%rsi)
       
#-----------------------------------

       movaps  512(%rdi), %xmm5
       movaps 1536(%rdi), %xmm13       
       movaps 2560(%rdi), %xmm12
       movaps 3584(%rdi), %xmm8

       movaps %xmm5,%xmm14
       subpd %xmm12,%xmm5
       addpd %xmm14,%xmm12
       
       movaps %xmm13,%xmm14
       subpd %xmm8,%xmm13
       addpd %xmm14,%xmm8
       
       movaps %xmm12,%xmm14
       subpd %xmm8,%xmm12
       addpd %xmm14,%xmm8
       
       movaps %xmm5,%xmm14
       subpd %xmm13,%xmm5
       addpd %xmm14,%xmm13
       
#-----------------------------------

       movaps .C2,%xmm15
       mulpd %xmm15,%xmm5
       mulpd %xmm15,%xmm13       
       movaps %xmm11,%xmm14
       subpd %xmm5,%xmm11
       addpd %xmm14,%xmm5

       movaps %xmm2,%xmm14
       subpd %xmm13,%xmm2
       addpd %xmm14,%xmm13
       
       movaps %xmm11,64(%rsi)
       movaps %xmm2,192(%rsi)
       
#-----------------------------------

       mulpd %xmm15,%xmm1
       mulpd %xmm15,%xmm3
       
       movaps .CI,%xmm15

       movaps %xmm1,%xmm11
       subpd %xmm3,%xmm1
       addpd %xmm11,%xmm3

       shufpd $0x1,%xmm12,%xmm12
       xorpd %xmm15,%xmm12
       
       shufpd $0x1,%xmm3,%xmm3
       xorpd %xmm15,%xmm3

       movaps %xmm9,%xmm11
       subpd %xmm1,%xmm9
       addpd %xmm11,%xmm1

       movaps %xmm3,%xmm11
       subpd %xmm12,%xmm3
       addpd %xmm11,%xmm12
       
       movaps %xmm1,%xmm11
       subpd %xmm12,%xmm1
       addpd %xmm11,%xmm12

       movaps %xmm9,%xmm11
       subpd %xmm3,%xmm9
       addpd %xmm11,%xmm3
       
       movaps %xmm12,32(%rsi)
       movaps %xmm3,96(%rsi)
       movaps %xmm9,160(%rsi)
       movaps %xmm1,224(%rsi)       
       
#-----------------------------------       
       
       movaps 16(%rsi), %xmm1

       movaps %xmm5,%xmm11
       subpd %xmm0,%xmm5
       addpd %xmm11,%xmm0
       
       movaps %xmm6,%xmm11
       subpd %xmm13,%xmm6
       addpd %xmm11,%xmm13
       shufpd $0x1,%xmm6,%xmm6
       xorpd %xmm15,%xmm6
       shufpd $0x1,%xmm13,%xmm13
       xorpd %xmm15,%xmm13

       movaps %xmm0,%xmm11
       subpd %xmm13,%xmm0
       addpd %xmm11,%xmm13
       
       movaps %xmm5,%xmm11
       subpd %xmm6,%xmm5
       addpd %xmm11,%xmm6
       
       movaps %xmm13,16(%rsi)              
       movaps %xmm6,112(%rsi)
       movaps %xmm5,144(%rsi)
       movaps %xmm0,240(%rsi)       
       
#-----------------------------------

       movaps (%rsi), %xmm6
       movaps 64(%rsi), %xmm0
       movaps 192(%rsi), %xmm13

       movaps %xmm6,%xmm11
       subpd %xmm8,%xmm6
       addpd %xmm11,%xmm8

       movaps %xmm1,%xmm11
       subpd %xmm4,%xmm1
       addpd %xmm11,%xmm4

       shufpd $0x1,%xmm1,%xmm1
       xorpd %xmm15,%xmm1

       movaps %xmm8,%xmm11
       subpd %xmm4,%xmm8
       addpd %xmm11,%xmm4
       
       movaps %xmm6,%xmm11
       subpd %xmm1,%xmm6
       addpd %xmm11,%xmm1
       
       movaps %xmm4,(%rsi)
       movaps %xmm1,64(%rsi)
       movaps %xmm8,128(%rsi)       
       movaps %xmm6,192(%rsi)
       
#-----------------------------------       
       
       movaps %xmm0,%xmm11
       subpd %xmm10,%xmm0
       addpd %xmm11,%xmm10

       movaps %xmm7,%xmm11
       subpd %xmm13,%xmm7
       addpd %xmm11,%xmm13

       shufpd $0x1,%xmm7,%xmm7
       xorpd %xmm15,%xmm7
       shufpd $0x1,%xmm13,%xmm13
       xorpd %xmm15,%xmm13
       
       movaps %xmm10,%xmm11
       subpd %xmm7,%xmm10
       addpd %xmm11,%xmm7
       
       movaps %xmm0,%xmm11
       subpd %xmm13,%xmm0
       addpd %xmm11,%xmm13
       
       movaps %xmm7,48(%rsi)
       movaps %xmm13,80(%rsi)
       movaps %xmm0,176(%rsi)
       movaps %xmm10,208(%rsi)              

       addq $16,%rdi
       addq $256,%rsi
       addq $512,%rcx       
#---------------------------       
.L0:
       movaps  256(%rdi), %xmm0
       movaps 1280(%rdi), %xmm2       
       movaps 2304(%rdi), %xmm1
       movaps 3328(%rdi), %xmm3

       movaps %xmm0,%xmm7
       subpd %xmm1,%xmm0
       addpd %xmm7,%xmm1

       movaps %xmm2,%xmm7
       subpd %xmm3,%xmm2
       addpd %xmm7,%xmm3

       movaps %xmm1,%xmm7
       subpd %xmm3,%xmm1
       addpd %xmm7,%xmm3

       movaps %xmm3,16(%rsi)
       
#-------------

       movaps  768(%rdi), %xmm5
       movaps 1792(%rdi), %xmm6       
       movaps 2816(%rdi), %xmm3
       movaps 3840(%rdi), %xmm4

       movaps %xmm5,%xmm7
       subpd %xmm3,%xmm5
       addpd %xmm7,%xmm3

       movaps %xmm6,%xmm7
       subpd %xmm4,%xmm6
       addpd %xmm7,%xmm4

       movaps %xmm3,%xmm7
       subpd %xmm4,%xmm3
       addpd %xmm7,%xmm4

#-------------

       movaps .C1,%xmm8
       movaps .C3,%xmm9       

       movaps %xmm0,%xmm7
       subpd %xmm6,%xmm0
       addpd %xmm7,%xmm6

       movaps %xmm5,%xmm7
       subpd %xmm2,%xmm5
       addpd %xmm7,%xmm2

       movaps %xmm8,%xmm7
       movaps %xmm9,%xmm10
       
       mulpd %xmm5,%xmm7
       mulpd %xmm0,%xmm10       
       subpd %xmm7,%xmm10
       
       mulpd %xmm8,%xmm0
       mulpd %xmm9,%xmm5
       addpd %xmm5,%xmm0
       
       movaps %xmm8,%xmm7
       
       mulpd %xmm6,%xmm7
       mulpd %xmm9,%xmm6
       
       mulpd %xmm2,%xmm9
       subpd %xmm9,%xmm7

       mulpd %xmm8,%xmm2
       addpd %xmm2,%xmm6

#-----------------------------------

       movaps     (%rdi), %xmm11
       movaps 2048(%rdi), %xmm9
       movaps 1024(%rdi), %xmm2
       movaps 3072(%rdi), %xmm8
       
       movaps %xmm11,%xmm14
       subpd %xmm9,%xmm11       
       addpd %xmm14,%xmm9

       movaps %xmm2,%xmm14
       subpd %xmm8,%xmm2
       addpd %xmm14,%xmm8
       
       movaps %xmm9,%xmm14
       subpd %xmm8,%xmm9
       addpd %xmm14,%xmm8

       movaps %xmm8,(%rsi)
       
#-----------------------------------

       movaps  512(%rdi), %xmm5
       movaps 1536(%rdi), %xmm13       
       movaps 2560(%rdi), %xmm12
       movaps 3584(%rdi), %xmm8

       movaps %xmm5,%xmm14
       subpd %xmm12,%xmm5
       addpd %xmm14,%xmm12
       
       movaps %xmm13,%xmm14
       subpd %xmm8,%xmm13
       addpd %xmm14,%xmm8
       
       movaps %xmm12,%xmm14
       subpd %xmm8,%xmm12
       addpd %xmm14,%xmm8
       
       movaps %xmm5,%xmm14
       subpd %xmm13,%xmm5
       addpd %xmm14,%xmm13
       
#-----------------------------------

       movaps .C2,%xmm15
       mulpd %xmm15,%xmm5
       mulpd %xmm15,%xmm13       
       movaps %xmm11,%xmm14
       subpd %xmm5,%xmm11
       addpd %xmm14,%xmm5

       movaps %xmm2,%xmm14
       subpd %xmm13,%xmm2
       addpd %xmm14,%xmm13
       
       movaps %xmm11,64(%rsi)
       movaps %xmm2,192(%rsi)
       
#-----------------------------------

       mulpd %xmm15,%xmm1
       mulpd %xmm15,%xmm3
       
       movaps .CI,%xmm15

       movaps %xmm1,%xmm11
       subpd %xmm3,%xmm1
       addpd %xmm11,%xmm3

       shufpd $0x1,%xmm12,%xmm12
       xorpd %xmm15,%xmm12
       
       shufpd $0x1,%xmm3,%xmm3
       xorpd %xmm15,%xmm3

       movaps %xmm9,%xmm11
       subpd %xmm1,%xmm9
       addpd %xmm11,%xmm1

       movaps %xmm3,%xmm11
       subpd %xmm12,%xmm3
       addpd %xmm11,%xmm12
       
       movaps %xmm1,%xmm11
       subpd %xmm12,%xmm1
       addpd %xmm11,%xmm12

       movaps %xmm9,%xmm11
       subpd %xmm3,%xmm9
       addpd %xmm11,%xmm3
       
###
       movaps %xmm12,%xmm11
       shufpd $1,%xmm12,%xmm12
       mulpd  64(%rcx),%xmm11
       mulpd  80(%rcx),%xmm12
       addpd %xmm11,%xmm12
###       
###
       movaps %xmm3,%xmm11
       shufpd $1,%xmm3,%xmm3
       mulpd 192(%rcx),%xmm11
       mulpd 208(%rcx),%xmm3
       addpd %xmm11,%xmm3
###       
###
       movaps %xmm9,%xmm11
       shufpd $1,%xmm9,%xmm9
       mulpd 320(%rcx),%xmm11
       mulpd 336(%rcx),%xmm9
       addpd %xmm11,%xmm9
###       
###
       movaps %xmm1,%xmm11
       shufpd $1,%xmm1,%xmm1
       mulpd 448(%rcx),%xmm11
       mulpd 464(%rcx),%xmm1
       addpd %xmm11,%xmm1
###       

       movaps %xmm12,32(%rsi)
       movaps %xmm3,96(%rsi)
       movaps %xmm9,160(%rsi)
       movaps %xmm1,224(%rsi)       
       
#-----------------------------------       
       
       movaps 16(%rsi), %xmm1

       movaps %xmm5,%xmm11
       subpd %xmm0,%xmm5
       addpd %xmm11,%xmm0
       
       movaps %xmm6,%xmm11
       subpd %xmm13,%xmm6
       addpd %xmm11,%xmm13
       shufpd $0x1,%xmm6,%xmm6
       xorpd %xmm15,%xmm6
       shufpd $0x1,%xmm13,%xmm13
       xorpd %xmm15,%xmm13

       movaps %xmm0,%xmm11
       subpd %xmm13,%xmm0
       addpd %xmm11,%xmm13
       
       movaps %xmm5,%xmm11
       subpd %xmm6,%xmm5
       addpd %xmm11,%xmm6
       
###
       movaps %xmm13,%xmm11
       shufpd $1,%xmm13,%xmm13
       mulpd  32(%rcx),%xmm11
       mulpd  48(%rcx),%xmm13
       addpd %xmm11,%xmm13
###       
###
       movaps %xmm6,%xmm11
       shufpd $1,%xmm6,%xmm6
       mulpd 224(%rcx),%xmm11
       mulpd 240(%rcx),%xmm6
       addpd %xmm11,%xmm6
###       
###
       movaps %xmm5,%xmm11
       shufpd $1,%xmm5,%xmm5
       mulpd 288(%rcx),%xmm11
       mulpd 304(%rcx),%xmm5
       addpd %xmm11,%xmm5
###       
###
       movaps %xmm0,%xmm11
       shufpd $1,%xmm0,%xmm0
       mulpd 480(%rcx),%xmm11
       mulpd 496(%rcx),%xmm0
       addpd %xmm11,%xmm0
###       

       movaps %xmm13,16(%rsi)              
       movaps %xmm6,112(%rsi)
       movaps %xmm5,144(%rsi)
       movaps %xmm0,240(%rsi)       
       
#-----------------------------------

       movaps (%rsi), %xmm6
       movaps 64(%rsi), %xmm0
       movaps 192(%rsi), %xmm13

       movaps %xmm6,%xmm11
       subpd %xmm8,%xmm6
       addpd %xmm11,%xmm8

       movaps %xmm1,%xmm11
       subpd %xmm4,%xmm1
       addpd %xmm11,%xmm4

       shufpd $0x1,%xmm1,%xmm1
       xorpd %xmm15,%xmm1

       movaps %xmm8,%xmm11
       subpd %xmm4,%xmm8
       addpd %xmm11,%xmm4
       
       movaps %xmm6,%xmm11
       subpd %xmm1,%xmm6
       addpd %xmm11,%xmm1
       
###
       movaps %xmm1,%xmm11
       shufpd $1,%xmm1,%xmm1
       mulpd 128(%rcx),%xmm11
       mulpd 144(%rcx),%xmm1
       addpd %xmm11,%xmm1
###       
###
       movaps %xmm8,%xmm11
       shufpd $1,%xmm8,%xmm8
       mulpd 256(%rcx),%xmm11
       mulpd 272(%rcx),%xmm8
       addpd %xmm11,%xmm8
###       
###
       movaps %xmm6,%xmm11
       shufpd $1,%xmm6,%xmm6
       mulpd 384(%rcx),%xmm11
       mulpd 400(%rcx),%xmm6
       addpd %xmm11,%xmm6
###       

       movaps %xmm4,(%rsi)
       movaps %xmm1,64(%rsi)
       movaps %xmm8,128(%rsi)       
       movaps %xmm6,192(%rsi)
       
#-----------------------------------       
       
       movaps %xmm0,%xmm11
       subpd %xmm10,%xmm0
       addpd %xmm11,%xmm10

       movaps %xmm7,%xmm11
       subpd %xmm13,%xmm7
       addpd %xmm11,%xmm13

       shufpd $0x1,%xmm7,%xmm7
       xorpd %xmm15,%xmm7
       shufpd $0x1,%xmm13,%xmm13
       xorpd %xmm15,%xmm13
       
       movaps %xmm10,%xmm11
       subpd %xmm7,%xmm10
       addpd %xmm11,%xmm7
       
       movaps %xmm0,%xmm11
       subpd %xmm13,%xmm0
       addpd %xmm11,%xmm13
       
###
       movaps %xmm7,%xmm11
       shufpd $1,%xmm7,%xmm7
       mulpd  96(%rcx),%xmm11
       mulpd 112(%rcx),%xmm7
       addpd %xmm11,%xmm7
###       
###
       movaps %xmm13,%xmm11
       shufpd $1,%xmm13,%xmm13
       mulpd 160(%rcx),%xmm11
       mulpd 176(%rcx),%xmm13
       addpd %xmm11,%xmm13
###       
###
       movaps %xmm0,%xmm11
       shufpd $1,%xmm0,%xmm0
       mulpd 352(%rcx),%xmm11
       mulpd 368(%rcx),%xmm0
       addpd %xmm11,%xmm0
###       
###
       movaps %xmm10,%xmm11
       shufpd $1,%xmm10,%xmm10
       mulpd 416(%rcx),%xmm11
       mulpd 432(%rcx),%xmm10
       addpd %xmm11,%xmm10
###       

       movaps %xmm7,48(%rsi)
       movaps %xmm13,80(%rsi)
       movaps %xmm0,176(%rsi)
       movaps %xmm10,208(%rsi)              

########################################
       addq $16,%rdi
       addq $256,%rsi
       addq $512,%rcx       
       cmpq %rdi,%r8
       jne .L0

########################       

#       subl $256,%eax
#       subl $4096,%ecx
       subq $256,%rdi
       subq $4096,%rsi
       
#-------------------------------------------------------------------------       
#------------------------------------
# fftj16 のループ(後処理)
#------------------------------------
       movq %rsi,%r8
       addq $256,%r8

.align 16       
.L1:
       movaps  256(%rsi), %xmm0
       movaps 1280(%rsi), %xmm2       
       movaps 2304(%rsi), %xmm1
       movaps 3328(%rsi), %xmm3

       movaps %xmm0,%xmm7
       subpd %xmm1,%xmm0
       addpd %xmm7,%xmm1

       movaps %xmm2,%xmm7
       subpd %xmm3,%xmm2
       addpd %xmm7,%xmm3

       movaps %xmm1,%xmm7
       subpd %xmm3,%xmm1
       addpd %xmm7,%xmm3

       movaps %xmm3,256(%rsi)
       
#-------------

       movaps  768(%rsi), %xmm5
       movaps 1792(%rsi), %xmm6       
       movaps 2816(%rsi), %xmm3
       movaps 3840(%rsi), %xmm4

       movaps %xmm5,%xmm7
       subpd %xmm3,%xmm5
       addpd %xmm7,%xmm3

       movaps %xmm6,%xmm7
       subpd %xmm4,%xmm6
       addpd %xmm7,%xmm4

       movaps %xmm3,%xmm7
       subpd %xmm4,%xmm3
       addpd %xmm7,%xmm4

#-------------

       movaps .C1,%xmm8
       movaps .C3,%xmm9       

       movaps %xmm0,%xmm7
       subpd %xmm6,%xmm0
       addpd %xmm7,%xmm6

       movaps %xmm5,%xmm7
       subpd %xmm2,%xmm5
       addpd %xmm7,%xmm2

       movaps %xmm8,%xmm7
       movaps %xmm9,%xmm10
       
       mulpd %xmm5,%xmm7
       mulpd %xmm0,%xmm10       
       subpd %xmm7,%xmm10
       
       mulpd %xmm8,%xmm0
       mulpd %xmm9,%xmm5
       addpd %xmm5,%xmm0
       
       movaps %xmm8,%xmm7
       
       mulpd %xmm6,%xmm7
       mulpd %xmm9,%xmm6
       
       mulpd %xmm2,%xmm9
       subpd %xmm9,%xmm7

       mulpd %xmm8,%xmm2
       addpd %xmm2,%xmm6

#-----------------------------------

       movaps     (%rsi), %xmm11
       movaps 2048(%rsi), %xmm9
       movaps 1024(%rsi), %xmm2
       movaps 3072(%rsi), %xmm8
       
       movaps %xmm11,%xmm14
       subpd %xmm9,%xmm11       
       addpd %xmm14,%xmm9

       movaps %xmm2,%xmm14
       subpd %xmm8,%xmm2
       addpd %xmm14,%xmm8
       
       movaps %xmm9,%xmm14
       subpd %xmm8,%xmm9
       addpd %xmm14,%xmm8

       movaps %xmm8,(%rsi)
       
#-----------------------------------

       movaps  512(%rsi), %xmm5
       movaps 1536(%rsi), %xmm13       
       movaps 2560(%rsi), %xmm12
       movaps 3584(%rsi), %xmm8

       movaps %xmm5,%xmm14
       subpd %xmm12,%xmm5
       addpd %xmm14,%xmm12
       
       movaps %xmm13,%xmm14
       subpd %xmm8,%xmm13
       addpd %xmm14,%xmm8
       
       movaps %xmm12,%xmm14
       subpd %xmm8,%xmm12
       addpd %xmm14,%xmm8
       
       movaps %xmm5,%xmm14
       subpd %xmm13,%xmm5
       addpd %xmm14,%xmm13
       
#-----------------------------------

       movaps .C2,%xmm15
       mulpd %xmm15,%xmm5
       mulpd %xmm15,%xmm13       
       movaps %xmm11,%xmm14
       subpd %xmm5,%xmm11
       addpd %xmm14,%xmm5

       movaps %xmm2,%xmm14
       subpd %xmm13,%xmm2
       addpd %xmm14,%xmm13
       
       movaps %xmm11,1024(%rsi)
       movaps %xmm2,3072(%rsi)
       
#-----------------------------------

       mulpd %xmm15,%xmm1
       mulpd %xmm15,%xmm3
       
       movaps .CI,%xmm15

       movaps %xmm1,%xmm11
       subpd %xmm3,%xmm1
       addpd %xmm11,%xmm3

       shufpd $0x1,%xmm12,%xmm12
       xorpd %xmm15,%xmm12
       
       shufpd $0x1,%xmm3,%xmm3
       xorpd %xmm15,%xmm3

       movaps %xmm9,%xmm11
       subpd %xmm1,%xmm9
       addpd %xmm11,%xmm1

       movaps %xmm3,%xmm11
       subpd %xmm12,%xmm3
       addpd %xmm11,%xmm12
       
       movaps %xmm1,%xmm11
       subpd %xmm12,%xmm1
       addpd %xmm11,%xmm12

       movaps %xmm9,%xmm11
       subpd %xmm3,%xmm9
       addpd %xmm11,%xmm3
       
       movaps %xmm12,512(%rsi)
       movaps %xmm3,1536(%rsi)
       movaps %xmm9,2560(%rsi)
       movaps %xmm1,3584(%rsi)       
       
#-----------------------------------       
       
       movaps 256(%rsi), %xmm1

       movaps %xmm5,%xmm11
       subpd %xmm0,%xmm5
       addpd %xmm11,%xmm0
       
       movaps %xmm6,%xmm11
       subpd %xmm13,%xmm6
       addpd %xmm11,%xmm13
       shufpd $0x1,%xmm6,%xmm6
       xorpd %xmm15,%xmm6
       shufpd $0x1,%xmm13,%xmm13
       xorpd %xmm15,%xmm13

       movaps %xmm0,%xmm11
       subpd %xmm13,%xmm0
       addpd %xmm11,%xmm13
       
       movaps %xmm5,%xmm11
       subpd %xmm6,%xmm5
       addpd %xmm11,%xmm6
       
       movaps %xmm13,256(%rsi)
       movaps %xmm6,1792(%rsi)
       movaps %xmm5,2304(%rsi)
       movaps %xmm0,3840(%rsi)       
       
#-----------------------------------

       movaps (%rsi), %xmm6
       movaps 1024(%rsi), %xmm0
       movaps 3072(%rsi), %xmm13

       movaps %xmm6,%xmm11
       subpd %xmm8,%xmm6
       addpd %xmm11,%xmm8

       movaps %xmm1,%xmm11
       subpd %xmm4,%xmm1
       addpd %xmm11,%xmm4

       shufpd $0x1,%xmm1,%xmm1
       xorpd %xmm15,%xmm1

       movaps %xmm8,%xmm11
       subpd %xmm4,%xmm8
       addpd %xmm11,%xmm4
       
       movaps %xmm6,%xmm11
       subpd %xmm1,%xmm6
       addpd %xmm11,%xmm1
       
       movaps %xmm4,(%rsi)
       movaps %xmm1,1024(%rsi)
       movaps %xmm8,2048(%rsi)       
       movaps %xmm6,3072(%rsi)
       
#-----------------------------------       
       
       movaps %xmm0,%xmm11
       subpd %xmm10,%xmm0
       addpd %xmm11,%xmm10

       movaps %xmm7,%xmm11
       subpd %xmm13,%xmm7
       addpd %xmm11,%xmm13

       shufpd $0x1,%xmm7,%xmm7
       xorpd %xmm15,%xmm7
       shufpd $0x1,%xmm13,%xmm13
       xorpd %xmm15,%xmm13
       
       movaps %xmm10,%xmm11
       subpd %xmm7,%xmm10
       addpd %xmm11,%xmm7
       
       movaps %xmm0,%xmm11
       subpd %xmm13,%xmm0
       addpd %xmm11,%xmm13
       
       movaps %xmm7,768(%rsi)
       movaps %xmm13,1280(%rsi)
       movaps %xmm0,2816(%rsi)
       movaps %xmm10,3328(%rsi)              

       addq $16,%rsi
       cmpq %rsi,%r8
       jne .L1
#-------------------------------------------------------------       
#  rfft のための後処理

      addq $8192,%r10
      movq $4096,%r9  # %r9 に (N/2)*16 が入る
 
      addq %r11,%r9
      movaps .CID,%xmm7
      
      movaps (%r11),%xmm0
      movaps %xmm0,%xmm1            
      xorpd %xmm7,%xmm0      
      shufpd $1,%xmm1,%xmm1      
      addpd %xmm1,%xmm0
      movapd %xmm0,(%r11)
      
      movaps  16(%r11),%xmm0
      movaps -16(%r9),%xmm1      
      xorpd %xmm7,%xmm0
      movaps %xmm1,%xmm2      
      subpd %xmm0,%xmm1
      
      movaps %xmm1,%xmm3
      shufpd $1,%xmm1,%xmm1
      mulpd 32(%r10),%xmm3
      mulpd 48(%r10),%xmm1
      addpd %xmm3,%xmm1
      
      addpd %xmm1,%xmm0
      subpd %xmm1,%xmm2
      xorpd %xmm7,%xmm2
      movaps %xmm0,16(%r11)      
      movaps %xmm2,-16(%r9)
      
      addq $32,%r11
      subq $32,%r9
      addq $64,%r10      

.align 16
.L1R:
      movaps (%r11),%xmm0
      movaps (%r9),%xmm1      
      xorpd %xmm7,%xmm0
      movaps %xmm1,%xmm2      
      subpd %xmm0,%xmm1
      
      movaps %xmm1,%xmm3
      shufpd $1,%xmm1,%xmm1
      mulpd (%r10),%xmm3
      mulpd 16(%r10),%xmm1
      addpd %xmm3,%xmm1
      
      addpd %xmm1,%xmm0
      subpd %xmm1,%xmm2
      xorpd %xmm7,%xmm2
      movaps %xmm0,(%r11)      
      movaps %xmm2,(%r9)

      movaps 16(%r11),%xmm0
      movaps -16(%r9),%xmm1      
      xorpd %xmm7,%xmm0
      movaps %xmm1,%xmm2      
      subpd %xmm0,%xmm1
      
      movaps %xmm1,%xmm3
      shufpd $1,%xmm1,%xmm1
      mulpd 32(%r10),%xmm3
      mulpd 48(%r10),%xmm1
      addpd %xmm3,%xmm1
      
      addpd %xmm1,%xmm0
      subpd %xmm1,%xmm2
      xorpd %xmm7,%xmm2
      movaps %xmm0,16(%r11)      
      movaps %xmm2,-16(%r9)

      addq $32,%r11
      subq $32,%r9
      addq $64,%r10      
      cmpq %r11,%r9
      jne .L1R

      movaps (%r11),%xmm0
      xorpd %xmm7,%xmm0
      movaps %xmm0,(%r11)

#-------------------------------------------------------------
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
.C1:
      .long 0xcf328d46,0x3fed906b,0xcf328d46,0x3fed906b
.align 16
.C2:
      .long 0x667f3bcd,0x3fe6a09e,0x667f3bcd,0x3fe6a09e
.align 16
.C3:
      .long 0xa6aea963,0x3fd87de2,0xa6aea963,0x3fd87de2
