!>@file   cal_sph_exp_fixed_scalar.f90
!!@brief  module cal_sph_exp_fixed_scalar
!!
!!@author H. Matsui
!!@date Programmed in Jan., 2010
!
!>@brief  Set fixed scalar boundarry for explicit method
!!
!!@verbatim
!!      subroutine dsdr_sph_fix_scalar_in_2(jmax, kr_in, r_ICB,         &
!!     &          fdm2_fix_fld_ICB, fix_ICB, is_fld, is_grd)
!!      subroutine dsdr_sph_lm0_fix_scalar_in_2(idx_rj_degree_zero,     &
!!     &          jmax, kr_in, r_ICB, fdm2_fix_fld_ICB, fix_ICB,        &
!!     &          is_fld, is_grd)
!!      subroutine cal_dsdr_sph_fix_scalar_in_2(jmax, kr_in, r_ICB,     &
!!     &          fdm2_fix_fld_ICB, fix_ICB, is_fld, is_grd)
!!      subroutine cal_dsdr_sph_no_bc_in_2(jmax, kr_in,                 &
!!     &          fdm2_fix_fld_ICB, is_fld, is_grd)
!!
!!      subroutine dsdr_sph_fix_scalar_out_2(jmax, kr_out, r_CMB,       &
!!     &          fdm2_fix_fld_CMB, fix_CMB, is_fld, is_grd)
!!      subroutine dsdr_sph_lm0_fix_scalar_out_2(idx_rj_degree_zero,    &
!!     &          jmax, kr_out, r_CMB, fdm2_fix_fld_CMB, fix_CMB,       &
!!     &          is_fld, is_grd)
!!      subroutine cal_sph_out_fix_scalar_diffuse2(jmax, kr_out, r_CMB, &
!!     &          fdm2_fix_fld_CMB, fix_CMB, coef_d, is_fld, is_diffuse)
!!      subroutine cal_dsdr_sph_no_bc_out_2(jmax, kr_out,               &
!!     &          fdm2_fix_fld_CMB, is_fld, is_grd)
!!
!!      subroutine cal_sph_div_flux_4_fix_in(jmax, kr_in, r_ICB,        &
!!     &          fdm2_fix_fld_ICB, fix_ICB, is_fld, is_div)
!!      subroutine cal_sph_div_flux_4_fix_out(jmax, kr_out, r_CMB,      &
!!     &          fdm2_fix_fld_CMB, fix_CMB, is_fld, is_div)
!!
!!      subroutine dsdr_sph_fixed_ctr_2(jmax, r_CTR1,                   &
!!     &          fdm2_fix_fld_center, is_fld, is_grd)
!!      subroutine dsdr_sph_lm0_fixed_ctr_2(idx_rj_degree_zero, jmax,   &
!!     &          r_CTR1, fdm2_fix_fld_center, fix_CTR, is_fld, is_grd)
!!      subroutine cal_sph_div_flux_4_fix_ctr(jmax, r_CTR1, fix_CTR,    &
!!     &          fdm2_fix_fld_center, is_fld, is_div)
!!      subroutine cal_sph_fixed_center_diffuse2(jmax, r_CTR1,          &
!!     &          fdm2_fix_fld_center, fix_CTR, coef_d,                 &
!!     &          is_fld, is_diffuse)
!!@endverbatim
!!
!!@n @param idx_rj_degree_zero    Local address for degree 0
!!@n @param jmax         Number of local spherical harmonics mode
!!@n @param kr_in       Radial ID for inner boundary
!!@n @param kr_out        Radial ID for outer boundary
!!@n @param r_ICB(0:2 )   Radius at ICB
!!@n @param r_CMB(0:2)    Radius at CMB
!!@n @param r_CTR1(0:2)   Radius at innermost point
!!@n @param fdm2_fix_fld_ICB(0:2,3)
!!         Matrix to evaluate radial derivative at ICB with fiexed field
!!@n @param fdm2_fix_fld_CMB(0:2,3)
!!         Matrix to evaluate radial derivative at CMB with fiexed field
!!@n @param fdm2_fix_fld_center(-1:1,3)
!!         Matrix to evaluate radial derivative
!!         for center with fiexed field
!!
!!@n @param fix_ICB(jmax) Spectr data for fixed fixed scalar at ICB
!!@n @param fix_CMB(jmax) Spectr data for fixed fixed scalar at CMB
!!@n @param fix_CTR(jmax) Spectr data for fixed fixed scalar at center
!!@n @param coef_d        Coefficient for diffusion term
!!
!!@n @param is_fld       Field address of input field
!!@n @param is_grd       Field address of radial gradient of field
!!@n @param is_diffuse   Field address for diffusion of field
!
      module cal_sph_exp_fixed_scalar
!
      use m_precision
!
      use m_constants
      use m_schmidt_poly_on_rtm
      use m_sph_spectr_data
!
      implicit none
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine dsdr_sph_fix_scalar_in_2(jmax, kr_in, r_ICB,           &
     &          fdm2_fix_fld_ICB, fix_ICB, is_fld, is_grd)
!
      integer(kind = kint), intent(in) :: is_fld, is_grd
      integer(kind = kint), intent(in) :: jmax, kr_in
      real(kind = kreal), intent(in) :: r_ICB(0:2)
      real(kind = kreal), intent(in) :: fix_ICB(jmax)
      real(kind = kreal), intent(in) :: fdm2_fix_fld_ICB(0:2,3)
!
      real(kind = kreal) :: d1t_dr1
      integer(kind = kint) :: inod, i_p1, i_p2, j
!
!
!$omp parallel do private(inod,i_p1,i_p2,d1t_dr1)
      do j = 1, jmax
        inod = j + (kr_in-1) * jmax
        i_p1 = inod + jmax
        i_p2 = i_p1 + jmax
!
        d1t_dr1 =  fdm2_fix_fld_ICB( 0,2) * fix_ICB(j)                  &
     &           + fdm2_fix_fld_ICB( 1,2) * d_rj(i_p1,is_fld)           &
     &           + fdm2_fix_fld_ICB( 2,2) * d_rj(i_p2,is_fld)
!
        d_rj(inod,is_fld  ) = fix_ICB(j)
        d_rj(inod,is_grd  ) = d1t_dr1 * g_sph_rj(j,13) * r_ICB(0)**2
        d_rj(inod,is_grd+1) = fix_ICB(j)
        d_rj(inod,is_grd+2) = zero
      end do
!$omp end parallel do
!
      end subroutine dsdr_sph_fix_scalar_in_2
!
! -----------------------------------------------------------------------
!
      subroutine dsdr_sph_lm0_fix_scalar_in_2(idx_rj_degree_zero,       &
     &          jmax, kr_in, r_ICB, fdm2_fix_fld_ICB, fix_ICB,          &
     &          is_fld, is_grd)
!
      integer(kind = kint), intent(in) :: idx_rj_degree_zero
      integer(kind = kint), intent(in) :: is_fld, is_grd
      integer(kind = kint), intent(in) :: jmax, kr_in
      real(kind = kreal), intent(in) :: r_ICB(0:2)
      real(kind = kreal), intent(in) :: fix_ICB(jmax)
      real(kind = kreal), intent(in) :: fdm2_fix_fld_ICB(0:2,3)
!
      real(kind = kreal) :: d1t_dr1
      integer(kind = kint) :: inod, i_p1, i_p2
!
!
      if(idx_rj_degree_zero .eq. 0) return
        inod = idx_rj_degree_zero + (kr_in-1) * jmax
        i_p1 = inod + jmax
        i_p2 = i_p1 + jmax
        d1t_dr1 = fdm2_fix_fld_ICB( 0,2) * fix_ICB(idx_rj_degree_zero)  &
     &          + fdm2_fix_fld_ICB( 1,2) * d_rj(i_p1,is_fld)            &
     &          + fdm2_fix_fld_ICB( 2,2) * d_rj(i_p2,is_fld)
!
        d_rj(inod,is_grd  ) = d1t_dr1 * r_ICB(0)**2
!
      end subroutine dsdr_sph_lm0_fix_scalar_in_2
!
! -----------------------------------------------------------------------
!
      subroutine cal_sph_fix_scalar_in_diffuse2(jmax, kr_in, r_ICB,     &
     &          fdm2_fix_fld_ICB, fix_ICB, coef_d, is_fld, is_diffuse)
!
      integer(kind = kint), intent(in) :: is_fld, is_diffuse
      integer(kind = kint), intent(in) :: jmax, kr_in
      real(kind = kreal), intent(in) :: fix_ICB(jmax)
      real(kind = kreal), intent(in) :: r_ICB(0:2)
      real(kind = kreal), intent(in) :: fdm2_fix_fld_ICB(0:2,3)
      real(kind = kreal), intent(in) :: coef_d
!
      real(kind = kreal) :: d1t_dr1, d2t_dr2
      integer(kind = kint) :: inod, i_p1, i_p2, j
!
!
!$omp parallel do private(inod,i_p1,i_p2,d1t_dr1,d2t_dr2)
      do j = 1, jmax
        inod = j + (kr_in-1) * jmax
        i_p1 = inod + jmax
        i_p2 = i_p1 + jmax
!
        d1t_dr1 =  fdm2_fix_fld_ICB( 0,2) * fix_ICB(j)                  &
     &           + fdm2_fix_fld_ICB( 1,2) * d_rj(i_p1,is_fld)           &
     &           + fdm2_fix_fld_ICB( 2,2) * d_rj(i_p2,is_fld)
        d2t_dr2 =  fdm2_fix_fld_ICB( 0,3) * fix_ICB(j)                  &
     &           + fdm2_fix_fld_ICB( 1,3) * d_rj(i_p1,is_fld)           &
     &           + fdm2_fix_fld_ICB( 2,3) * d_rj(i_p2,is_fld)
!
        d_rj(inod,is_fld) = fix_ICB(j)
        d_rj(inod,is_diffuse)                                           &
     &         = coef_d * (d2t_dr2 + two*r_ICB(1)*d1t_dr1               &
     &          - g_sph_rj(j,3)*r_ICB(2) * d_rj(inod,is_fld) )
!
      end do
!$omp end parallel do
!
      end subroutine cal_sph_fix_scalar_in_diffuse2
!
! -----------------------------------------------------------------------
!
      subroutine cal_dsdr_sph_no_bc_in_2(jmax, kr_in,                   &
     &          fdm2_fix_fld_ICB, is_fld, is_grd)
!
      integer(kind = kint), intent(in) :: jmax, kr_in
      integer(kind = kint), intent(in) :: is_fld, is_grd
      real(kind = kreal), intent(in) :: fdm2_fix_fld_ICB(0:2,3)
!
      integer(kind = kint) :: inod, i_p1, i_p2, j
!
!
!$omp parallel do private(inod,i_p1,i_p2)
      do j = 1, jmax
        inod = j + (kr_in-1) * jmax
        i_p1 = inod + jmax
        i_p2 = i_p1 + jmax
!
        d_rj(inod,is_grd) =  fdm2_fix_fld_ICB(0,2) * d_rj(inod,is_fld)  &
     &                     + fdm2_fix_fld_ICB(1,2) * d_rj(i_p1,is_fld)  &
     &                     + fdm2_fix_fld_ICB(2,2) * d_rj(i_p2,is_fld)
      end do
!$omp end parallel do
!
      end subroutine cal_dsdr_sph_no_bc_in_2
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine dsdr_sph_fix_scalar_out_2(jmax, kr_out, r_CMB,         &
     &          fdm2_fix_fld_CMB, fix_CMB, is_fld, is_grd)
!
      integer(kind = kint), intent(in) :: is_fld, is_grd
      integer(kind = kint), intent(in) :: jmax, kr_out
      real(kind = kreal), intent(in) :: r_CMB(0:2)
      real(kind = kreal), intent(in) :: fix_CMB(jmax)
      real(kind = kreal), intent(in) :: fdm2_fix_fld_CMB(0:2,3)
!
      real(kind = kreal) :: d1t_dr1
      integer(kind = kint) :: inod, i_n1, i_n2, j
!
!
!$omp parallel do private(inod,i_n1,i_n2,d1t_dr1)
      do j = 1, jmax
        inod = j + (kr_out-1) * jmax
        i_n1 = inod - jmax
        i_n2 = i_n1 - jmax
!
        d1t_dr1 =  fdm2_fix_fld_CMB(2,2) * d_rj(i_n2,is_fld)            &
     &           + fdm2_fix_fld_CMB(1,2) * d_rj(i_n1,is_fld)            &
     &           + fdm2_fix_fld_CMB(0,2) * fix_CMB(j)
!
        d_rj(inod,is_fld  ) = fix_CMB(j)
        d_rj(inod,is_grd  ) = d1t_dr1 * g_sph_rj(j,13) * r_CMB(0)**2
        d_rj(inod,is_grd+1) = fix_CMB(j)
        d_rj(inod,is_grd+2) = zero
      end do
!$omp end parallel do
!
      end subroutine dsdr_sph_fix_scalar_out_2
!
! -----------------------------------------------------------------------
!
      subroutine dsdr_sph_lm0_fix_scalar_out_2(idx_rj_degree_zero,      &
     &          jmax, kr_out, r_CMB, fdm2_fix_fld_CMB, fix_CMB,         &
     &          is_fld, is_grd)
!
      integer(kind = kint), intent(in) :: idx_rj_degree_zero
      integer(kind = kint), intent(in) :: is_fld, is_grd
      integer(kind = kint), intent(in) :: jmax, kr_out
      real(kind = kreal), intent(in) :: r_CMB(0:2)
      real(kind = kreal), intent(in) :: fix_CMB(jmax)
      real(kind = kreal), intent(in) :: fdm2_fix_fld_CMB(0:2,3)
!
      real(kind = kreal) :: d1t_dr1
      integer(kind = kint) :: inod, i_n1, i_n2
!
!
      inod = idx_rj_degree_zero + (kr_out-1) * jmax
      i_n1 = inod - jmax
      i_n2 = i_n1 - jmax
!
      d1t_dr1 = fdm2_fix_fld_CMB(2,2) * d_rj(i_n2,is_fld)               &
     &          + fdm2_fix_fld_CMB(1,2) * d_rj(i_n1,is_fld)             &
     &          + fdm2_fix_fld_CMB(0,2) * fix_CMB(idx_rj_degree_zero)
!
      d_rj(inod,is_grd  ) = d1t_dr1 * r_CMB(0)**2
!
      end subroutine dsdr_sph_lm0_fix_scalar_out_2
!
! -----------------------------------------------------------------------
!
      subroutine cal_sph_out_fix_scalar_diffuse2(jmax, kr_out, r_CMB,   &
     &          fdm2_fix_fld_CMB, fix_CMB, coef_d, is_fld, is_diffuse)
!
      integer(kind = kint), intent(in) :: jmax, kr_out
      integer(kind = kint), intent(in) :: is_fld, is_diffuse
      real(kind = kreal), intent(in) :: coef_d
      real(kind = kreal), intent(in) :: fix_CMB(jmax)
      real(kind = kreal), intent(in) :: r_CMB(0:2)
      real(kind = kreal), intent(in) :: fdm2_fix_fld_CMB(0:2,3)
!
      real(kind = kreal) :: d1t_dr1, d2t_dr2
      integer(kind = kint) :: inod, i_n1, i_n2, j
!
!
!$omp parallel do private(inod,i_n1,i_n2,d2t_dr2,d1t_dr1)
      do j = 1, jmax
        inod = j + (kr_out-1) * jmax
        i_n1 = inod - jmax
        i_n2 = i_n1 - jmax
!
        d1t_dr1 =  fdm2_fix_fld_CMB(2,2) * d_rj(i_n2,is_fld)            &
     &           + fdm2_fix_fld_CMB(1,2) * d_rj(i_n1,is_fld)            &
     &           + fdm2_fix_fld_CMB(0,2) * fix_CMB(j)
        d2t_dr2 =  fdm2_fix_fld_CMB(2,3) * d_rj(i_n2,is_fld)            &
     &           + fdm2_fix_fld_CMB(1,3) * d_rj(i_n1,is_fld)            &
     &           + fdm2_fix_fld_CMB(0,3) * fix_CMB(j)
!
        d_rj(inod,is_fld) = fix_CMB(j)
        d_rj(inod,is_diffuse)                                           &
     &          =  coef_d * (d2t_dr2 + two*r_CMB(1) * d1t_dr1           &
     &            - g_sph_rj(j,3)*r_CMB(2) * d_rj(inod,is_fld))
!
      end do
!$omp end parallel do
!
      end subroutine cal_sph_out_fix_scalar_diffuse2
!
! -----------------------------------------------------------------------
!
      subroutine cal_dsdr_sph_no_bc_out_2(jmax, kr_out,                 &
     &          fdm2_fix_fld_CMB, is_fld, is_grd)
!
      integer(kind = kint), intent(in) :: jmax, kr_out
      integer(kind = kint), intent(in) :: is_fld, is_grd
      real(kind = kreal), intent(in) :: fdm2_fix_fld_CMB(0:2,3)
!
      integer(kind = kint) :: inod, i_n1, i_n2, j
!
!
!$omp parallel do private(inod,i_n1,i_n2)
      do j = 1, jmax
        inod = j + (kr_out-1) * jmax
        i_n1 = inod - jmax
        i_n2 = i_n1 - jmax
!
        d_rj(inod,is_grd) = fdm2_fix_fld_CMB(2,2) * d_rj(i_n2,is_fld) &
     &                    + fdm2_fix_fld_CMB(1,2) * d_rj(i_n1,is_fld) &
     &                    + fdm2_fix_fld_CMB(0,2) * d_rj(inod,is_fld)
      end do
!$omp end parallel do
!
      end subroutine cal_dsdr_sph_no_bc_out_2
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine cal_sph_div_flux_4_fix_in(jmax, kr_in, r_ICB,          &
     &          fdm2_fix_fld_ICB, fix_ICB, is_fld, is_div)
!
      integer(kind = kint), intent(in) :: is_fld, is_div
      integer(kind = kint), intent(in) :: jmax, kr_in
      real(kind = kreal), intent(in) :: fix_ICB(jmax)
      real(kind = kreal), intent(in) :: r_ICB(0:2)
      real(kind = kreal), intent(in) :: fdm2_fix_fld_ICB(0:2,3)
!
      real(kind = kreal) :: d1s_dr1
      integer(kind = kint) :: inod, i_p1, i_p2, j
!
!
!$omp parallel do private(inod,i_p1,i_p2,j,d1s_dr1)
!cdir nodep
      do j = 1, jmax
        inod = j + (kr_in-1) * jmax
        i_p1 = inod + jmax
        i_p2 = i_p1 + jmax
!
        d1s_dr1 =  fdm2_fix_fld_ICB( 0,2) * fix_ICB(j)                  &
     &           + fdm2_fix_fld_ICB( 1,2) * d_rj(i_p1,is_fld)           &
     &           + fdm2_fix_fld_ICB( 2,2) * d_rj(i_p2,is_fld)
!
        d_rj(inod,is_div) =  (d1s_dr1 - d_rj(inod,is_fld+1) )           &
     &                     * max(g_sph_rj(j,3),half) * r_ICB(2)
      end do
!$omp end parallel do
!
      end subroutine cal_sph_div_flux_4_fix_in
!
! -----------------------------------------------------------------------
!
      subroutine cal_sph_div_flux_4_fix_out(jmax, kr_out, r_CMB,        &
     &          fdm2_fix_fld_CMB, fix_CMB, is_fld, is_div)
!
      integer(kind = kint), intent(in) :: is_fld, is_div
      integer(kind = kint), intent(in) :: jmax, kr_out
      real(kind = kreal), intent(in) :: fix_CMB(jmax)
      real(kind = kreal), intent(in) :: r_CMB(0:2)
      real(kind = kreal), intent(in) :: fdm2_fix_fld_CMB(0:2,3)
!
      real(kind = kreal) :: d1s_dr1
      integer(kind = kint) :: inod, i_n1, i_n2, j
!
!
!$omp parallel do private(inod,i_n1,i_n2,j,d1s_dr1)
!cdir nodep
      do j = 1, jmax
        inod = j + (kr_out-1) * jmax
        i_n1 = inod - jmax
        i_n2 = i_n1 - jmax
!
        d1s_dr1 =  fdm2_fix_fld_CMB(2,2) * d_rj(i_n2,is_fld)            &
     &           + fdm2_fix_fld_CMB(1,2) * d_rj(i_n1,is_fld)            &
     &           + fdm2_fix_fld_CMB(0,2) * fix_CMB(j)
!
        d_rj(inod,is_div) =  (d1s_dr1 - d_rj(inod,is_fld+1) )           &
     &                     * max(g_sph_rj(j,3),half) * r_CMB(2)
      end do
!$omp end parallel do
!
      end subroutine cal_sph_div_flux_4_fix_out
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine dsdr_sph_fixed_ctr_2(jmax, r_CTR1,                     &
     &          fdm2_fix_fld_center, is_fld, is_grd)
!
      integer(kind = kint), intent(in) :: jmax
      integer(kind = kint), intent(in) :: is_fld, is_grd
      real(kind = kreal), intent(in) :: r_CTR1(0:2)
      real(kind = kreal), intent(in) :: fdm2_fix_fld_center(-1:1,3)
!
      integer(kind = kint) :: inod, i_p1
      real(kind = kreal) :: d1sdr
!
!
!$omp parallel do private(inod,i_p1,d1sdr)
      do inod = 1, jmax
        i_p1 = inod + jmax
!
        d1sdr =  fdm2_fix_fld_center( 0,2) * d_rj(inod,is_fld)          &
     &         + fdm2_fix_fld_center( 1,2) * d_rj(i_p1,is_fld)
!
        d_rj(inod,is_grd  ) = d1sdr * g_sph_rj(inod,13) * r_CTR1(0)**2
        d_rj(inod,is_grd+1) = d_rj(inod,is_fld)
        d_rj(inod,is_grd+2) = zero
      end do
!$omp end parallel do
!
      end subroutine dsdr_sph_fixed_ctr_2
!
! -----------------------------------------------------------------------
!
      subroutine dsdr_sph_lm0_fixed_ctr_2(idx_rj_degree_zero, jmax,     &
     &          r_CTR1, fdm2_fix_fld_center, fix_CTR, is_fld, is_grd)
!
      integer(kind = kint), intent(in) :: idx_rj_degree_zero
      integer(kind = kint), intent(in) :: jmax
      integer(kind = kint), intent(in) :: is_fld, is_grd
      real(kind = kreal), intent(in) :: fix_CTR(jmax)
      real(kind = kreal), intent(in) :: r_CTR1(0:2)
      real(kind = kreal), intent(in) :: fdm2_fix_fld_center(-1:1,3)
!
      integer(kind = kint) :: inod, i_p1
      real(kind = kreal) :: d1sdr
!
!
      if(idx_rj_degree_zero .eq. 0) return
!
      inod = idx_rj_degree_zero
      i_p1 = inod + jmax
!
      d1sdr =  fdm2_fix_fld_center(-1,2) * fix_CTR(inod)                &
     &         + fdm2_fix_fld_center( 0,2) * d_rj(inod,is_fld)          &
     &         + fdm2_fix_fld_center( 1,2) * d_rj(i_p1,is_fld)
!
      d_rj(inod,is_grd) = d1sdr * r_CTR1(0)**2
!
      end subroutine dsdr_sph_lm0_fixed_ctr_2
!
! -----------------------------------------------------------------------
!
      subroutine cal_sph_div_flux_4_fix_ctr(jmax, r_CTR1, fix_CTR,      &
     &          fdm2_fix_fld_center, is_fld, is_div)
!
      integer(kind = kint), intent(in) :: jmax
      integer(kind = kint), intent(in) :: is_fld, is_div
      real(kind = kreal), intent(in) :: r_CTR1(0:2)
      real(kind = kreal), intent(in) :: fix_CTR(jmax)
      real(kind = kreal), intent(in) :: fdm2_fix_fld_center(-1:1,3)
!
      real(kind = kreal) :: d1s_dr1
      integer(kind = kint) :: i_p1, j
!
!
!$omp parallel do private(i_p1,j,d1s_dr1)
!cdir nodep
      do j = 1, jmax
        i_p1 = j + jmax
!
        d1s_dr1 =  fdm2_fix_fld_center(-1,2) * fix_CTR(j)               &
     &           + fdm2_fix_fld_center( 0,2) * d_rj(j,is_fld)           &
     &           + fdm2_fix_fld_center( 1,2) * d_rj(i_p1,is_fld)
!
        d_rj(j,is_div) =  (d1s_dr1 - d_rj(j,is_fld+1) )                 &
     &                   * max(g_sph_rj(j,3),half) * r_CTR1(2)
      end do
!$omp end parallel do
!
      end subroutine cal_sph_div_flux_4_fix_ctr
!
! -----------------------------------------------------------------------
!
      subroutine cal_sph_fixed_center_diffuse2(jmax, r_CTR1,            &
     &          fdm2_fix_fld_center, fix_CTR, coef_d,                   &
     &          is_fld, is_diffuse)
!
      integer(kind = kint), intent(in) :: is_fld, is_diffuse
      integer(kind = kint), intent(in) :: jmax
      real(kind = kreal), intent(in) :: fix_CTR(jmax)
      real(kind = kreal), intent(in) :: r_CTR1(0:2)
      real(kind = kreal), intent(in) :: fdm2_fix_fld_center(-1:1,3)
      real(kind = kreal), intent(in) :: coef_d
!
      real(kind = kreal) :: d1t_dr1, d2t_dr2
      integer(kind = kint) :: i_p1, i_p2, j
!
!
!$omp parallel do private(i_p1,i_p2,d1t_dr1,d2t_dr2)
      do j = 1, jmax
        i_p1 = j +    jmax
        i_p2 = i_p1 + jmax
!
        d1t_dr1 =  fdm2_fix_fld_center(-1,2) * fix_CTR(j)               &
     &           + fdm2_fix_fld_center( 0,2) * d_rj(i_p1,is_fld)        &
     &           + fdm2_fix_fld_center( 1,2) * d_rj(i_p2,is_fld)
        d2t_dr2 =  fdm2_fix_fld_center(-1,3) * fix_CTR(j)               &
     &           + fdm2_fix_fld_center( 0,3) * d_rj(i_p1,is_fld)        &
     &           + fdm2_fix_fld_center( 1,3) * d_rj(i_p2,is_fld)
!
        d_rj(j,is_diffuse)                                              &
     &         = coef_d * (d2t_dr2 + two*r_CTR1(1)*d1t_dr1              &
     &          - g_sph_rj(j,3)*r_CTR1(2) * d_rj(j,is_fld) )
!
      end do
!$omp end parallel do
!
      end subroutine cal_sph_fixed_center_diffuse2
!
! -----------------------------------------------------------------------
!
      end module cal_sph_exp_fixed_scalar
