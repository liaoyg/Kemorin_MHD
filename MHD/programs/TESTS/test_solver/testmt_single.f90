      program testmt
!***************************************
!*
!*   test program for solve matrix
!*
!***************************************
!*
!***************************************
!*
!*        |  1  3  0                |
!*        | -4  2  6  0             |
!*        |  0  1  4  2  0          |
!*    A = |     0  3 -2 -9  0       |
!*        |        0  1  3 -6  0    |
!*        |           0  2  6  1  0 |
!*        |              0  4 -1  3 |
!*        |                 0 -1  5 |
!*
!*        |  1 -4  0                |
!*        |  3  2  1  0             |
!*     t  |  0  6  4  3  0          |
!*    A = |     0  2 -2  1  0       |
!*        |        0  9  3  2  0    |
!*        |           0 -6  6  4  0 |
!*        |              0  1 -1 -1 |
!*        |                 0  3  5 |
!*
!***************************************
!*
!*
!*        |   3 |            | -18 |             | -6 |
!*        |  42 |            | -10 |             |  3 |
!*        |  13 |            |  27 |             |  2 |
!*    b = | -32 |       c =  |   6 |        x =  |  1 |
!*        |  31 |            | - 3 |             |  4 |
!*        | - 5 |            | -22 |             | -3 |
!*        |   1 |            | -14 |             |  5 |
!*        |  25 |            |  45 |             |  6 |
!*
!*
!***************************************
!*
!**********************************************************************
!*
!*       define of a(i,j)
!*
!*               | a(1,1)  a(1,2)  ........  a(1,n)  |
!*               | a(2,1)  a(2,2)  ........  a(2,n)  |
!*          A =  |   .       .      a(i,j)      .    |
!*               | a(n,1)  a(n,2)  ........  a(n,n)  |
!*
!*       subroutine
!*
!      call ludcmp(a,N,NP,ip,ep)      : LU factrization
!      call lubksb(a,N,NP,ip,x)     :  solve
!
!**********************************************************************
!*
!
!**********************************************************************
!     Kemo's 3-band solver
!
!       a(i,i) =   dg(i)
!       a(i+1,i) = al(i)
!       a(i,i+1) = au(i)
!
!*               | dg(1)  au(1)  ........     0         0     |
!*               | al(1)  dg(2)  ........     .         .     |
!*               |   0    al(2)  ........     .         .     |
!*    a(i,j)  =  |   .       0   ........     0         .     |
!*               | ...... al(k-1)  dg(k)     au(k) .......... |
!*               |   .       .   ........  au(N-2)     0      |
!*               |   .       .   ........  dg(N-1)  au(N-1)   |
!*               |   0       0   ........  al(N-1)  dg(N)     |
!
!      SUBROUTINE ludcmp_3band(n, dg, al, au)
!      SUBROUTINE lubksb_3band(n, dg, al, au, x)
!
!**********************************************************************
!*
!*
!**********************************************************************
!*
!*     band matrix for IMSL subroutine 
!*
!*               | a(2,1)  a(1,2)  ........     0         0     |
!*               | a(3,1)  a(2,2)  ........     .         .     |
!*               |   0     a(3,2)  ........     .         .     |
!*    a(i,j)  =  |   .       0     ........     0         .     |
!*               | ...... a(3,k-1)  a(2,k)  a(1,k+1) .......... |
!*               |   .       .     ........  a(1,N-2)     0     |
!*               |   .       .     ........  a(2,N-2)  a(1,N-1) |
!*               |   0       0     ........  a(3,N-2)  a(2,N-1) |
!*
!*    a(1,1) = a(3,N-1) = a(m,j) = 0.0    { m > 4 }
!*
!*       subroutine
!*
!*    dlftrb ( N-1 ,a ,N-1 ,1 ,1 ,a ,N-1 ,ip )       : LU factrization
!*    dlfsrb ( N-1 ,a ,N-1 ,1 ,1 ,ip ,b ,1 ,x )      :  solve
!*
!**********************************************************************
!*
!**********************************************************************
!*
!*     band matrix for SX3r super computer
!*
!*               | a(2,1)  a(3,1)  ........     0         0     |
!*               | a(1,2)  a(2,2)  ........     .         .     |
!*               |   0     a(1,3)  ........     .         .     |
!*    a(i,j)  =  |   .       0     ........     0         .     |
!*               | ........ a(1,k)  a(2,k)  a(3,k) ............ |
!*               |   .       .     ........  a(3,N-3)     0     |
!*               |   .       .     ........  a(2,N-2)  a(3,N-2) |
!*               |   0       0     ........  a(1,N-1)  a(2,N-1) |
!*
!*    a(1,1) = a(3,N-1) = a(m,j) = 0.0    { m > 4 }
!*
!*       subroutine
!*
!*    dbbdlu ( a ,N-1 ,N-1 ,1 ,1 ,ip ,ic )       : LU factrization
!*    dbbdls ( a ,N-1 ,N-1 ,1 ,1 ,b ,ip ,ic )    :  solve
!*
!**********************************************************************
!*
      use m_precision
      use m_ludcmp
      use m_ludcmp_3band
      use lubksb_357band
!
      implicit none
!*
!* ------  define  --------------
!*
      integer(kind = kint), parameter :: ismp = 4, nvect = 8, ncp = 8
      real(kind = kreal) :: a(ncp,ncp) ,band_a(3,ncp), band_lu(5,ncp)
      real(kind = kreal) :: b(ncp), x(ncp)
      real(kind = kreal) :: ep
!*
      integer(kind = kint) :: ip(ncp)
      integer(kind = kint) :: i, ierr, j, k
!*
       data a/ 1.0 ,-4.0,0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0                   &
     &        ,3.0 ,2.0 ,1.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0                   &
     &        ,0.0 ,6.0 ,4.0 ,3.0 ,0.0 ,0.0 ,0.0 ,0.0                   &
     &        ,0.0 ,0.0 ,2.0 ,-2.0,1.0 ,0.0 ,0.0 ,0.0                   &
     &        ,0.0 ,0.0 ,0.0 ,-9.0,3.0 ,2.0 ,0.0 ,0.0                   &
     &        ,0.0 ,0.0 ,0.0 ,0.0 ,-6.0,6.0 ,4.0 ,0.0                   &
     &        ,0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,1.0 ,-1.0,-1.0                  &
     &        ,0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,3.0 ,5.0 /
       data b/ 3.0 , 42.0 , 13.0 ,-32.0 , 31.0 ,-5.0 , 1.0 , 25.0/
!*
!*  -------    set band matrix ab ----------
!*
      do 9 i = 1 ,ncp
        write(6,600) (a(i,j) ,j=1,ncp)
 600    format(1p8e10.1)
   9  continue
!*
!*
!*
      band_a = 0.0d0
      band_a(2,1) = a(1,1)
      band_a(3,1) = a(2,1)
      do j = 2, ncp-1
        band_a(1,j) = a(j-1,j)
        band_a(2,j) = a(j,  j)
        band_a(3,j) = a(j+1,j)
      end do
      band_a(1,ncp) = a(ncp-1,ncp)
      band_a(2,ncp) = a(ncp,  ncp)
!*
!*  --------  solve vector by 3band -------------
!*
!
      x(1:ncp) = b(1:ncp)
!
      call ludcmp_3band(ncp, band_a, ip, ierr, band_lu, ep)
!
      write(6,*) 'alu by 3band'
      do i = 1, 5
        write(6,600) (band_lu(i,j) ,j=1,ncp)
      end do
!
      call lubksb_3band(ncp,band_lu,ip,x)
!
      write(6,*) 'RHS and soultion by band matrix',ierr
      do 20 k = 1 ,ncp
        write(6,*) ip(k), b(k), x(k)
  20  continue
!*
!*  --------  solve vector by full -------------
!*
!
      x(1:ncp) = b(1:ncp)
!
      call ludcmp(a,ncp,ncp,ip,ep)
!
      write(6,*) 'alu by full'
      do 8 i = 1 ,ncp
        write(6,600) (a(i,j) ,j=1,ncp)
   8  continue
!
      call lubksb(a,ncp,ncp,ip,x)
!
      write(6,*) 'RHS and soultion by full matrix'
      do 21 k = 1 ,ncp
        write(6,*) b(k) ,x(k), ip(k)
  21  continue
!*
      stop
      end
