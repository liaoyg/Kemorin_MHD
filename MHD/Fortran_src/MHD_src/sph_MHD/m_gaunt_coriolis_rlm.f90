!>@file   m_gaunt_coriolis_rlm.f90
!!@brief  module m_gaunt_coriolis_rlm
!!
!!@author H. Matsui
!!@date Programmed in 1994
!@n     Modified in Dec., 2013
!
!>@brief Adams-Gaunt integrals for Coriolis term
!!       and coefficients for Coriolis term on f(r,l,m)
!!
!!@verbatim
!!      subroutine alloacte_gaunt_coriolis_rlm(jmax_rlm)
!!      subroutine alloc_coriolis_coef_tri_rlm(jmax_rlm)
!!      subroutine dealloacte_gaunt_coriolis_rlm
!!      subroutine dealloc_coriolis_coef_tri_rlm
!!
!!      integer function find_local_sph_rlm_address(jmax_rlm,           &
!!     &       idx_gl_1d_rlm_j, j_gl)
!!
!!*************************************************
!!*
!!*  Rotation of the Coriolos term
!!*     (wss) = wss(jc,1,j3)*w*dyb/r**2
!!*            + wss(jc,2,j3)*dw*yb/r**2
!!*
!!*     (wts) = wts(j3)*w*yb/r**2
!!*
!!*     (wst) = wst(1,j3)*( dw*dyb/r**2 + w*d2yb/r**2 - 2*w*dyb/r**3 )
!!*            + wst(2,j3)*( d2w/r**2 - 2*dw/r**3 )*yb
!!*
!!*     (wtt) = wtt(jc,1,j3)*dw*yb/r**2
!!*            + wtt(jc,2,j3)*w*( dyb/r**2 - 2*yb/r**3 )
!!*
!!*   Divergence of the Coriolis term
!!*     (wsd) = wsd(jc,1,j3)*w*wsb/r**4
!!*            + wsd(jc,2,j3)*dw*dwsb/r**2
!!*     (wtd) = wtd(j3)*dw*dwtb/r**2
!!*
!!*  Radial componenet of the Coriolis term
!!*     (wsr) = wsr(jc,1,j3)*dw*dusb/r**2
!!*     (wtr) = wtr(j3)*dw*wtb/r**2
!!
!!*************************************************
!!*
!!*************************************************
!!*
!!*     wss(jc,1,j3) = sw_rlm(jc,1,j3)
!!*     wss(jc,2,j3) = sw_rlm(jc,2,j3)
!!*     wts(jc,j3)   = sw_rlm(jc,3,j3)
!!*     wst(jc,1,j3) = tw_rlm(jc,1,j3)
!!*     wst(jc,2,j3) = tw_rlm(jc,2,j3)
!!*     wtt(jc,1,j3) = tw_rlm(jc,3,j3)
!!*     wtt(jc,2,j3) = tw_rlm(jc,4,j3)
!!*
!!*     wsd(jc,1,j3) = sd_rlm(jc,1,j3)
!!*     wsd(jc,2,j3) = sd_rlm(jc,2,j3)
!!*     wtd(jc,j3)   = td_rlm(jc,j3)
!!*
!!*     wsr(jc,j3) =   sr_rlm(jc,j3)
!!*     wtr(jc,j3) =   tr_rlm(jc,j3)
!!*
!!*************************************************
!!*
!!*************************************************************
!!*
!!*      gk_cor(j3,jk,2) : gaunt integral for Coriolis term
!!*              jgi_cor_rlm(j3,1) : index for gi_cor_rlm(j3,1)
!!*              jgi_cor_rlm(j3,2) : index for gi_cor_rlm(j3,2)
!!*      ei_cor_rlm(j3,1) : elsasser integral for Coriolis term
!!*              jei_cor_rlm(j3,1) : index for ei_cor_rlm(j3,1)
!!*
!!*************************************************************
!!*
!!*******************************************************************
!!*                                                                 *
!!*  Adams - Gaunt integrals                                        *
!!*                                                                 *
!!*  (m1,m2,m3)  //  (m1)   (m2)   (m3)                             *
!!* Ki         = || Y    * Y    * Y     sin(theta) d(theta)d(phi)   *
!!*  (l1,l2,l3)  //  (l1)   (l2)   (l3)                             *
!!*                                                                 *
!!*                            (m2)        (m3)                     *
!!*  (m1,m2,m3)  //  (m1)    dy(l2)      dy(l3)                     *
!!* Li         = || Y    *[ -------- * ----------                   *
!!*  (l1,l2,l3)  //  (l1)   d(theta)    d(phi)                      *
!!*                                                                 *
!!*                    (m2)        (m3)                             *
!!*                  dy(l2)      dy(l3)                             *
!!*              - ---------- * -------- ] d(theta)d(phi)           *
!!*                  d(phi)     d(theta)                            *
!!*                                                                 *
!!*  where                                                          *
!!*                   (m)   (m)  | sin(m*phi) |                     *
!!*                  Y   = P   * |   1        |                     *
!!*                   (l)   (l)  | cos(m*phi) |                     *
!!*                                                                 *
!!*                         (m)     2(l-m)!     1                   *
!!*                        P   = [----------]**--- * P(l,m)         *
!!*                         (l)      (l+m)!     2                   *
!!*                         (0)                                     *
!!*                        P   =           P(l,0)                   *
!!*                         (l)                                     *
!!*                                                                 *
!!*******************************************************************
!!@endverbatim
!!
      module m_gaunt_coriolis_rlm
!*
      use m_precision
      use m_constants
!
      implicit none
!
!>      Local address for Coriolis term using Gaunt integral
      integer(kind = kint), allocatable :: jgi_cor_rlm(:,:)
!>      Local address for Coriolis term using Elsasser integral
      integer(kind = kint), allocatable :: jei_cor_rlm(:,:)
!
!>      Gaunt integral for Coriolis term
      real(kind = kreal), allocatable :: gi_cor_rlm(:,:)
!>      Elsasser integral for Coriolis term
      real(kind = kreal), allocatable :: ei_cor_rlm(:,:)
!
!
!>      Coefficients for curl of Coriolis force for poloidal vorticity
      real(kind = kreal), allocatable :: sw_rlm(:,:,:)
!>      Coefficients for curl of Coriolis force for toroidal vorticity
      real(kind = kreal), allocatable :: tw_rlm(:,:,:)
!>      Coefficients for divergence of Coriolis force 
!!       for poloidal vorticity by poloidal velocity
      real(kind = kreal), allocatable :: sd_rlm(:,:,:)
!>      Coefficients for divergence of Coriolis force 
!!       for poloidal vorticity by Toroidal velocity
      real(kind = kreal), allocatable :: td_rlm(:,:)
!>      Coefficients for radial compoonent of Coriolis force 
!!       for poloidal vorticity by poloidal velocity
      real(kind = kreal), allocatable :: sr_rlm(:,:)
!>      Coefficients for radial compoennt of Coriolis force 
!!       for poloidal vorticity by Toroidal velocity
      real(kind = kreal), allocatable :: tr_rlm(:,:)
!
!*   ------------------------------------------------------------------
!*
      contains
!*
!*   ------------------------------------------------------------------
!
      subroutine alloacte_gaunt_coriolis_rlm(jmax_rlm)
!
      integer(kind = kint), intent(in) :: jmax_rlm
!
!
      allocate( jgi_cor_rlm(jmax_rlm,2) )
      allocate( jei_cor_rlm(jmax_rlm,1) )
      allocate( gi_cor_rlm(jmax_rlm,2) )
      allocate( ei_cor_rlm(jmax_rlm,1) )
!
      jgi_cor_rlm = 0
      jei_cor_rlm = 0
      gi_cor_rlm = 0.0d0
      ei_cor_rlm = 0.0d0
!
      end subroutine alloacte_gaunt_coriolis_rlm
!
!-----------------------------------------------------------------------
!
      subroutine alloc_coriolis_coef_tri_rlm(jmax_rlm)
!
      integer(kind = kint), intent(in) :: jmax_rlm
!
!
      allocate( sw_rlm(2,3,jmax_rlm) )
      allocate( tw_rlm(2,4,jmax_rlm) )
!
      allocate( sd_rlm(2,2,jmax_rlm) )
      allocate( td_rlm(2,jmax_rlm) )
!
      allocate( sr_rlm(2,jmax_rlm) )
      allocate( tr_rlm(2,jmax_rlm) )
!
      sw_rlm = zero
      tw_rlm = zero
!
      sd_rlm = zero
      td_rlm = zero
!
      sr_rlm = zero
      tr_rlm = zero
!
      end subroutine alloc_coriolis_coef_tri_rlm
!
! ----------------------------------------------------------------------
!
      subroutine dealloacte_gaunt_coriolis_rlm
!
!
      deallocate(jgi_cor_rlm, gi_cor_rlm)
      deallocate(jei_cor_rlm, ei_cor_rlm)
!
      end subroutine dealloacte_gaunt_coriolis_rlm
!
!-----------------------------------------------------------------------
!
      subroutine dealloc_coriolis_coef_tri_rlm
!
      deallocate(sw_rlm, tw_rlm)
      deallocate(sd_rlm, td_rlm)
      deallocate(sr_rlm, tr_rlm)
!
      end subroutine dealloc_coriolis_coef_tri_rlm
!
! ----------------------------------------------------------------------
!
      integer function find_local_sph_rlm_address(jmax_rlm,             &
     &       idx_gl_1d_rlm_j, j_gl)
!
      integer(kind = kint), intent(in) :: jmax_rlm
      integer(kind = kint), intent(in) :: idx_gl_1d_rlm_j(jmax_rlm,3)
      integer(kind = kint), intent(in) :: j_gl
!
      integer(kind = kint) :: j, l, m
!
!
      l = int( aint(sqrt(dble(j_gl))) )
      m = j_gl - l*(l+1)
!
      find_local_sph_rlm_address = 0
      do j = 1, jmax_rlm
        if (   idx_gl_1d_rlm_j(j,2) .eq. l                              &
     &   .and. idx_gl_1d_rlm_j(j,3) .eq. m) then
          find_local_sph_rlm_address = j
          return
        end if
      end do
!
      end function find_local_sph_rlm_address
!
!-----------------------------------------------------------------------
!
      end module m_gaunt_coriolis_rlm
