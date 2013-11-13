!>@file   cal_sph_exp_nod_icb_ins.f90
!!@brief  module cal_sph_exp_nod_icb_ins
!!
!!@author H. Matsui
!!@date Programmed in Jan., 2010
!!@date Modified in   Nov., 2013
!
!>@brief  Set insulated magnetic boundary condition for ICB
!!
!!@verbatim
!!      subroutine cal_sph_nod_icb_ins_b_and_j(jmax, kr_in,             &
!!     &           fdm2_fix_fld_ICB, fdm2_fix_dr_ICB, is_fld, is_rot)
!!      subroutine cal_sph_nod_icb_ins_mag2(jmax, kr_in, is_fld)
!!
!!      subroutine cal_sph_nod_icb_ins_vp_rot2(jmax, kr_in,             &
!!     &          is_fld, is_rot)
!!      subroutine cal_sph_nod_icb_ins_rot2(jmax, kr_in,                &
!!     &          fdm2_fix_fld_ICB, fdm2_fix_dr_ICB, is_fld, is_rot)
!!      subroutine cal_sph_nod_icb_ins_diffuse2(jmax, kr_in,            &
!!     &          fdm2_fix_fld_ICB, fdm2_fix_dr_ICB, coef_d,            &
!!     &          is_fld, is_diffuse)
!!@endverbatim
!!
!!@n @param jmax  Number of modes for spherical harmonics @f$L*(L+2)@f$
!!@n @param kr_in       Radial ID for inner boundary
!!@n @param fdm2_fix_fld_ICB(0:2,3)
!!!        Matrix to evaluate radial derivative at ICB with fiexed field
!!@n @param fdm2_fix_dr_ICB(-1:1,3)
!!         Matrix to evaluate field at ICB with fiexed radial derivative
!!
!!@n @param coef_d  Coefficient for diffusion term
!!@n @param is_fld       Field address of input field
!!@n @param is_rot       Field address for curl of field
!!@n @param is_diffuse   Field address for diffusion of field
!
      module cal_sph_exp_nod_icb_ins
!
      use m_precision
!
      use m_constants
      use m_spheric_parameter
      use m_schmidt_poly_on_rtm
      use m_sph_spectr_data
      use m_fdm_coefs
!
      implicit none
!
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine cal_sph_nod_icb_ins_b_and_j(jmax, kr_in,               &
     &           fdm2_fix_fld_ICB, fdm2_fix_dr_ICB, is_fld, is_rot)
!
      use m_coef_fdm_fixed_ICB
!
      integer(kind = kint), intent(in) :: jmax, kr_in
      integer(kind = kint), intent(in) :: is_fld, is_rot
      real(kind = kreal), intent(in) :: fdm2_fix_fld_ICB(0:2,3)
      real(kind = kreal), intent(in) :: fdm2_fix_dr_ICB(-1:1,3)
!
      real(kind = kreal) :: d1s_dr1, d2s_dr2, d1t_dr1
      integer(kind = kint) :: j, inod, i_p1, i_p2
!
!
!$omp parallel do private(inod,i_p1,i_p2,d1s_dr1,d2s_dr2,d1t_dr1)
      do j = 1, jmax
        inod = j + (kr_in-1) * jmax
        i_p1 = inod + jmax
        i_p2 = i_p1 + jmax
!
        d1s_dr1 =  (g_sph_rj(j,1)+one) * ar_1d_rj(kr_in,1)              &
     &            * d_rj(inod,is_fld)
        d2s_dr2 =  ( (g_sph_rj(j,1)+one) * ar_1d_rj(kr_in,1)            &
     &          *  fdm2_fix_dr_ICB(-1,3)                                &
     &           + fdm2_fix_dr_ICB( 0,3) ) * d_rj(inod,is_fld  )        &
     &           + fdm2_fix_dr_ICB( 1,3) * d_rj(i_p1,is_fld  )
        d1t_dr1 =  fdm2_fix_fld_ICB( 1,2) * d_rj(i_p1,is_fld+2)         &
     &           + fdm2_fix_fld_ICB( 2,2) * d_rj(i_p2,is_fld+2)
!
        d_rj(inod,is_fld+1) = d1s_dr1
        d_rj(inod,is_fld+2) = zero
        d_rj(inod,is_rot  ) = zero
        d_rj(inod,is_rot+1) = d1t_dr1
        d_rj(inod,is_rot+2) = - ( d2s_dr2 - g_sph_rj(j,3)               &
     &                   *ar_1d_rj(kr_in,2)*d_rj(inod,is_fld  ) )
      end do
!$omp end parallel do
!
      end subroutine cal_sph_nod_icb_ins_b_and_j
!
! -----------------------------------------------------------------------
!
      subroutine cal_sph_nod_icb_ins_mag2(jmax, kr_in, is_fld)
!
      integer(kind = kint), intent(in) :: jmax, kr_in
      integer(kind = kint), intent(in) :: is_fld
!
      real(kind = kreal) :: d1s_dr1
      integer(kind = kint) :: j, inod
!
!
!$omp parallel do private(inod,d1s_dr1)
      do j = 1, jmax
        inod = j + (kr_in-1) * jmax
        d1s_dr1 =  (g_sph_rj(j,1)+one) * ar_1d_rj(kr_in,1)              &
     &            * d_rj(inod,is_fld)
!
        d_rj(inod,is_fld+1) = d1s_dr1
        d_rj(inod,is_fld+2) = zero
      end do
!$omp end parallel do
!
      end subroutine cal_sph_nod_icb_ins_mag2
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine cal_sph_nod_icb_ins_vp_rot2(jmax, kr_in,               &
     &          is_fld, is_rot)
!
      integer(kind = kint), intent(in) :: jmax, kr_in
      integer(kind = kint), intent(in) :: is_fld, is_rot
!
      real(kind = kreal) :: d1t_dr1
      integer(kind = kint) :: j, inod
!
!
!$omp parallel do private(inod,d1t_dr1)
      do j = 1, jmax
        inod = j + (kr_in-1) * jmax
        d1t_dr1 =  (g_sph_rj(j,1)+one) * ar_1d_rj(kr_in,1)              &
     &            * d_rj(inod,is_fld+2)
!
        d_rj(inod,is_rot  ) = d_rj(inod,is_fld+2)
        d_rj(inod,is_rot+1) = d1t_dr1
        d_rj(inod,is_rot+2) = zero
      end do
!$omp end parallel do
!
      end subroutine cal_sph_nod_icb_ins_vp_rot2
!
! -----------------------------------------------------------------------
!
      subroutine cal_sph_nod_icb_ins_rot2(jmax, kr_in,                  &
     &          fdm2_fix_fld_ICB, fdm2_fix_dr_ICB, is_fld, is_rot)
!
      use m_coef_fdm_fixed_ICB
!
      integer(kind = kint), intent(in) :: jmax, kr_in
      integer(kind = kint), intent(in) :: is_fld, is_rot
      real(kind = kreal), intent(in) :: fdm2_fix_fld_ICB(0:2,3)
      real(kind = kreal), intent(in) :: fdm2_fix_dr_ICB(-1:1,3)
!
      real(kind = kreal) :: d2s_dr2, d1t_dr1
      integer(kind = kint) :: j, inod, i_p1, i_p2
!
!
!$omp parallel do private(inod,i_p1,i_p2,d2s_dr2,d1t_dr1)
      do j = 1, jmax
        inod = j + (kr_in-1) * jmax
        i_p1 = inod + jmax
        i_p2 = i_p1 + jmax
!
        d2s_dr2 =  ( (g_sph_rj(j,1)+one) * ar_1d_rj(kr_in,1)            &
     &          *  fdm2_fix_dr_ICB(-1,3)                                &
     &           + fdm2_fix_dr_ICB( 0,3) ) * d_rj(inod,is_fld  )        &
     &           + fdm2_fix_dr_ICB( 1,3) * d_rj(i_p1,is_fld  )
        d1t_dr1 =  fdm2_fix_fld_ICB( 0,2) * d_rj(inod,is_fld+2)         &
     &           + fdm2_fix_fld_ICB( 1,2) * d_rj(i_p1,is_fld+2)         &
     &           + fdm2_fix_fld_ICB( 2,2) * d_rj(i_p2,is_fld+2)
!
        d_rj(inod,is_rot  ) = d_rj(inod,is_fld+2)
        d_rj(inod,is_rot+1) = d1t_dr1
        d_rj(inod,is_rot+2) = - ( d2s_dr2                               &
     &    - g_sph_rj(j,3)*ar_1d_rj(kr_in,2)*d_rj(inod,is_fld  ) )
      end do
!$omp end parallel do
!
      end subroutine cal_sph_nod_icb_ins_rot2
!
! -----------------------------------------------------------------------
!
      subroutine cal_sph_nod_icb_ins_diffuse2(jmax, kr_in,              &
     &          fdm2_fix_fld_ICB, fdm2_fix_dr_ICB, coef_d,              &
     &          is_fld, is_diffuse)
!
      use m_coef_fdm_fixed_ICB
!
      integer(kind = kint), intent(in) :: jmax, kr_in
      integer(kind = kint), intent(in) :: is_fld, is_diffuse
      real(kind = kreal), intent(in) :: coef_d
      real(kind = kreal), intent(in) :: fdm2_fix_fld_ICB(0:2,3)
      real(kind = kreal), intent(in) :: fdm2_fix_dr_ICB(-1:1,3)
!
      real(kind = kreal) :: d2s_dr2,d2t_dr2
      integer(kind = kint) :: j, inod,i_p1,i_p2
!
!$omp parallel do private(inod,i_p1,i_p2,d2s_dr2,d2t_dr2)
      do j = 1, jmax
        inod = j + (kr_in-1) * jmax
        i_p1 = inod + jmax
        i_p2 = i_p1 + jmax
!
        d2s_dr2 =  ( (g_sph_rj(j,1)+one) * ar_1d_rj(kr_in,1)            &
     &           * fdm2_fix_dr_ICB(-1,3)                                &
     &           + fdm2_fix_dr_ICB( 0,3) ) * d_rj(inod,is_fld  )        &
     &           + fdm2_fix_dr_ICB( 1,3) * d_rj(i_p1,is_fld  )
        d2t_dr2 =  fdm2_fix_fld_ICB( 0,3) * d_rj(inod,is_fld+2)         &
     &           + fdm2_fix_fld_ICB( 1,3) * d_rj(i_p1,is_fld+2)         &
     &           + fdm2_fix_fld_ICB( 2,3) * d_rj(i_p2,is_fld+2)
!
        d_rj(inod,is_diffuse  ) = coef_d * (d2s_dr2                     &
     &    - g_sph_rj(j,3)*ar_1d_rj(kr_in,2)*d_rj(inod,is_fld  ) )
        d_rj(inod,is_diffuse+2) = coef_d * (d2t_dr2                     &
     &    - g_sph_rj(j,3)*ar_1d_rj(kr_in,2)*d_rj(inod,is_fld+2) )
      end do
!$omp end parallel do
!
      end subroutine cal_sph_nod_icb_ins_diffuse2
!
! -----------------------------------------------------------------------
!
      end module cal_sph_exp_nod_icb_ins
