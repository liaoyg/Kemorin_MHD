!> @file  cal_inner_core_rotation.f90
!!      module cal_inner_core_rotation
!!
!! @author  H. Matsui
!! @date Programmed in Nov., 2012
!
!> @brief Evaluate torques for inner core rotation
!!
!!@verbatim
!!      subroutine set_inner_core_rotation(kr_in, rj_fld)
!!        type(phys_data), intent(inout) :: rj_fld
!!      subroutine set_icore_viscous_matrix(kr_in, fdm1_fix_fld_ICB)
!!      subroutine cal_icore_viscous_drag_explicit(kr_in,               &
!!     &          fdm1_fix_fld_ICB, coef_d, it_velo, it_viscous,        &
!!     &          ntot_phys_rj, d_rj)
!!      subroutine copy_icore_rot_to_tor_coriolis                       &
!!     &         (kr_in, ntot_phys_rj, d_rj)
!!      subroutine inner_core_coriolis_rj                               &
!!     &         (kr_in, idx_rj_degree_one, nnod_rj, nri, jmax,         &
!!     &          radius_1d_rj_r, ntot_phys_rj, d_rj)
!!      subroutine int_icore_toroidal_lorentz(kr_in, rj_fld)
!!        type(phys_data), intent(inout) :: rj_fld
!!@endverbatim
!!
!!@n @param coef_d  Coefficient for diffusion term
!!@n @param it_velo       Field address for toroidal velocity
!!@n @param it_viscous    Field address for toroidal viscous dissipation
!
      module cal_inner_core_rotation
!
      use m_precision
!
      use m_constants
      use m_sph_phys_address
!
      implicit  none
!
      private :: int_icore_tor_lorentz_l1, cal_icore_viscous_drag_l1
      private :: set_rotate_icb_vt_sph_mat, set_inner_core_rot_l1
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine set_inner_core_rotation(kr_in, rj_fld)
!
      use m_spheric_parameter
      use t_phys_data
!
      integer(kind = kint), intent(in) :: kr_in
      type(phys_data), intent(inout) :: rj_fld
!
!
      call set_inner_core_rot_l1(sph_rj1%idx_rj_degree_one(-1),         &
     &    kr_in, sph_rj1%nidx_rj(1), sph_rj1%nidx_rj(2),                &
     &    sph_rj1%radius_1d_rj_r, sph_rj1%ar_1d_rj,                     &
     &    rj_fld%n_point, rj_fld%ntot_phys, rj_fld%d_fld)
      call set_inner_core_rot_l1(sph_rj1%idx_rj_degree_one( 0),         &
     &    kr_in, sph_rj1%nidx_rj(1), sph_rj1%nidx_rj(2),                &
     &    sph_rj1%radius_1d_rj_r, sph_rj1%ar_1d_rj,                     &
     &    rj_fld%n_point, rj_fld%ntot_phys, rj_fld%d_fld)
      call set_inner_core_rot_l1(sph_rj1%idx_rj_degree_one( 1),         &
     &    kr_in, sph_rj1%nidx_rj(1), sph_rj1%nidx_rj(2),                &
     &    sph_rj1%radius_1d_rj_r, sph_rj1%ar_1d_rj,                     &
     &    rj_fld%n_point, rj_fld%ntot_phys, rj_fld%d_fld)
!
      end subroutine set_inner_core_rotation
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine set_icore_viscous_matrix(kr_in, fdm1_fix_fld_ICB)
!
      use m_t_int_parameter
      use m_spheric_parameter
      use m_physical_property
!
      integer(kind = kint), intent(in) :: kr_in
      real(kind = kreal), intent(in) :: fdm1_fix_fld_ICB(0:1,2)
!
!
      call set_rotate_icb_vt_sph_mat(sph_rj1%idx_rj_degree_one(-1),     &
     &    kr_in, sph_rj1%nidx_rj(1), sph_rj1%ar_1d_rj,                  &
     &    fdm1_fix_fld_ICB, coef_imp_v, coef_d_velo)
      call set_rotate_icb_vt_sph_mat(sph_rj1%idx_rj_degree_one( 0),     &
     &    kr_in, sph_rj1%nidx_rj(1), sph_rj1%ar_1d_rj,                  &
     &    fdm1_fix_fld_ICB, coef_imp_v, coef_d_velo)
      call set_rotate_icb_vt_sph_mat(sph_rj1%idx_rj_degree_one( 1),     &
     &    kr_in, sph_rj1%nidx_rj(1), sph_rj1%ar_1d_rj,                  &
     &    fdm1_fix_fld_ICB, coef_imp_v, coef_d_velo)
!!
      end subroutine set_icore_viscous_matrix
!
! ----------------------------------------------------------------------
!
      subroutine cal_icore_viscous_drag_explicit(kr_in,                 &
     &          fdm1_fix_fld_ICB, coef_d, it_velo, it_viscous,          &
     &          ntot_phys_rj, d_rj)
!
      use m_spheric_parameter
!
      integer(kind = kint), intent(in) :: kr_in
      integer(kind = kint), intent(in) :: it_velo, it_viscous
      real(kind = kreal), intent(in) :: coef_d
      real(kind = kreal), intent(in) :: fdm1_fix_fld_ICB(0:1,2)
      integer(kind = kint), intent(in) :: ntot_phys_rj
!
      real(kind = kreal), intent(inout) :: d_rj(nnod_rj,ntot_phys_rj)
!
!
      call cal_icore_viscous_drag_l1(sph_rj1%idx_rj_degree_one(-1),     &
     &    kr_in, fdm1_fix_fld_ICB, coef_d, it_velo, it_viscous,         &
     &    sph_rj1%nidx_rj(1), sph_rj1%nidx_rj(2), sph_rj1%ar_1d_rj,     &
     &    nnod_rj, ntot_phys_rj, d_rj)
      call cal_icore_viscous_drag_l1(sph_rj1%idx_rj_degree_one( 0),     &
     &    kr_in, fdm1_fix_fld_ICB, coef_d, it_velo, it_viscous,         &
     &    sph_rj1%nidx_rj(1), sph_rj1%nidx_rj(2), sph_rj1%ar_1d_rj,     &
     &    nnod_rj, ntot_phys_rj, d_rj)
      call cal_icore_viscous_drag_l1(sph_rj1%idx_rj_degree_one( 1),     &
     &    kr_in, fdm1_fix_fld_ICB, coef_d, it_velo, it_viscous,         &
     &    sph_rj1%nidx_rj(1), sph_rj1%nidx_rj(2), sph_rj1%ar_1d_rj,     &
     &    nnod_rj, ntot_phys_rj, d_rj)
!
      end subroutine cal_icore_viscous_drag_explicit
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine copy_icore_rot_to_tor_coriolis                         &
     &         (kr_in, idx_rj_degree_one, nnod_rj, jmax,                &
     &          ntot_phys_rj, d_rj)
!
      integer(kind = kint), intent(in) :: idx_rj_degree_one(-1:1)
      integer(kind = kint), intent(in) :: nnod_rj, jmax
      integer(kind = kint), intent(in) :: kr_in
      integer(kind = kint), intent(in) :: ntot_phys_rj
!
      real(kind = kreal), intent(inout) :: d_rj(nnod_rj,ntot_phys_rj)
!
      integer(kind = kint) :: m, i1
!
!
      do m = -1, 1
        if(idx_rj_degree_one(m) .gt. 0) then
          i1 = idx_rj_degree_one(m) + (kr_in-1)*jmax
          d_rj(i1,itor%i_coriolis) = d_rj(i1,ipol%i_rot_Coriolis)
        end if
      end do
!
      end subroutine copy_icore_rot_to_tor_coriolis
!
! ----------------------------------------------------------------------
!
      subroutine inner_core_coriolis_rj                                 &
     &         (kr_in, idx_rj_degree_one, nnod_rj, nri, jmax,           &
     &          radius_1d_rj_r, ntot_phys_rj, d_rj)
!
      use m_physical_property
      use m_poloidal_rotation
!
      integer(kind = kint), intent(in) :: idx_rj_degree_one(-1:1)
      integer(kind = kint), intent(in) :: nnod_rj, nri, jmax
      integer(kind = kint), intent(in) :: kr_in
      integer(kind = kint), intent(in) :: ntot_phys_rj
      real(kind= kreal), intent(in) :: radius_1d_rj_r(nri)
!
      real(kind = kreal), intent(inout) :: d_rj(nnod_rj,ntot_phys_rj)
!
      integer(kind = kint) :: i11s, i10c, i11c
!
!
      if(idx_rj_degree_one( 1) .le. 0) return
!
!
      i11s = idx_rj_degree_one(-1) + (kr_in-1)*jmax
      i10c = idx_rj_degree_one( 0) + (kr_in-1)*jmax
      i11c = idx_rj_degree_one( 1) + (kr_in-1)*jmax
!
      d_rj(i11s,ipol%i_rot_Coriolis)                                    &
     &       =  omega_rj(kr_in,0,2)*d_rj(i11c,ipol%i_vort)              &
     &        - omega_rj(kr_in,0,3)*d_rj(i10c,ipol%i_vort)
      d_rj(i11c,ipol%i_rot_Coriolis)                                    &
     &       =  omega_rj(kr_in,0,1)*d_rj(i10c,ipol%i_vort)              &
     &        - omega_rj(kr_in,0,2)*d_rj(i11s,ipol%i_vort)
      d_rj(i10c,ipol%i_rot_Coriolis)                                    &
     &       =  omega_rj(kr_in,0,3)*d_rj(i11s,ipol%i_vort)              &
     &        - omega_rj(kr_in,0,1)*d_rj(i11c,ipol%i_vort)
!
      d_rj(i11s,ipol%i_rot_Coriolis)                                    &
     &       = -two*coef_cor*radius_1d_rj_r(kr_in)                      &
     &        * d_rj(i11s,ipol%i_rot_Coriolis)
      d_rj(i11c,ipol%i_rot_Coriolis)                                    &
     &       = -two*coef_cor*radius_1d_rj_r(kr_in)                      &
     &        * d_rj(i11c,ipol%i_rot_Coriolis)
      d_rj(i10c,ipol%i_rot_Coriolis)                                    &
     &       = -two*coef_cor*radius_1d_rj_r(kr_in)                      &
     &        * d_rj(i10c,ipol%i_rot_Coriolis)
!
      end subroutine inner_core_coriolis_rj
!
! ----------------------------------------------------------------------
!
      subroutine int_icore_toroidal_lorentz(kr_in, rj_fld)
!
      use m_spheric_parameter
      use t_phys_data
!
      integer(kind = kint), intent(in) :: kr_in
      type(phys_data), intent(inout) :: rj_fld
!
!
      call int_icore_tor_lorentz_l1(sph_rj1%idx_rj_degree_one(-1),      &
     &    kr_in, sph_rj1%nidx_rj(1), sph_rj1%nidx_rj(2),                &
     &    sph_rj1%radius_1d_rj_r, sph_rj1%ar_1d_rj,                     &
     &    rj_fld%n_point, rj_fld%ntot_phys, rj_fld%d_fld)
      call int_icore_tor_lorentz_l1(sph_rj1%idx_rj_degree_one( 0),      &
     &    kr_in, sph_rj1%nidx_rj(1), sph_rj1%nidx_rj(2),                &
     &    sph_rj1%radius_1d_rj_r, sph_rj1%ar_1d_rj,                     &
     &    rj_fld%n_point, rj_fld%ntot_phys, rj_fld%d_fld)
      call int_icore_tor_lorentz_l1(sph_rj1%idx_rj_degree_one( 1),      &
     &    kr_in, sph_rj1%nidx_rj(1), sph_rj1%nidx_rj(2),                &
     &    sph_rj1%radius_1d_rj_r, sph_rj1%ar_1d_rj,                     &
     &    rj_fld%n_point, rj_fld%ntot_phys, rj_fld%d_fld)
!
      end subroutine int_icore_toroidal_lorentz
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine set_rotate_icb_vt_sph_mat(idx_rj_l0, kr_in,            &
     &          nri, ar_1d_rj, fdm1_fix_fld_ICB, coef_imp, coef_d)
!
      use m_t_int_parameter
      use m_schmidt_poly_on_rtm
      use m_radial_matrices_sph
      use m_fdm_coefs
!
      integer(kind = kint), intent(in) :: nri
      integer(kind = kint), intent(in) :: kr_in, idx_rj_l0
      real(kind = kreal), intent(in) :: fdm1_fix_fld_ICB(0:1,2)
      real(kind = kreal), intent(in) :: coef_imp, coef_d
      real(kind= kreal), intent(in) :: ar_1d_rj(nri,3)
!
!
      if(idx_rj_l0 .le. 0) return
!
!       vt_evo_mat(3,kr_in-1,idx_rj_l0) = zero
        vt_evo_mat(2,kr_in,  idx_rj_l0)                                 &
     &          = one - coef_imp*dt*coef_d * five                       &
     &           * (fdm1_fix_fld_ICB(0,2) - two*ar_1d_rj(kr_in,1) )     &
     &           * ar_1d_rj(kr_in,1)
        vt_evo_mat(1,kr_in+1,idx_rj_l0)                                 &
     &          = - coef_imp*dt*coef_d * five * ar_1d_rj(kr_in,1)       &
     &             * fdm1_fix_fld_ICB(1,2)
!
      end subroutine set_rotate_icb_vt_sph_mat
!
! -----------------------------------------------------------------------
!
      subroutine set_inner_core_rot_l1(idx_rj_l0, kr_in,                &
     &          nri, jmax, radius_1d_rj_r, ar_1d_rj,                    &
     &          nnod_rj, ntot_phys_rj, d_rj)
!
      integer(kind = kint), intent(in) :: nri, jmax
      integer(kind = kint), intent(in) :: kr_in, idx_rj_l0
      integer(kind = kint), intent(in) :: nnod_rj, ntot_phys_rj
      real(kind= kreal), intent(in) :: radius_1d_rj_r(nri)
      real(kind= kreal), intent(in) :: ar_1d_rj(nri,3)
!
      real(kind = kreal), intent(inout) :: d_rj(nnod_rj,ntot_phys_rj)
!
      integer(kind = kint) :: i10c, i10c_ri
      integer(kind = kint) :: k
      real(kind = kreal) :: ratio
!
!
      if(idx_rj_l0 .le. 0) return
!
      i10c_ri = idx_rj_l0 + (kr_in-1)*jmax
!$omp parallel do private(k,i10c,ratio)
      do k = 1, kr_in-1
        i10c = idx_rj_l0 + (k-1)*jmax
!
        ratio = radius_1d_rj_r(k)*radius_1d_rj_r(k) * ar_1d_rj(kr_in,2)
!
        d_rj(i10c,itor%i_velo) =   ratio * d_rj(i10c_ri,itor%i_velo)
        d_rj(i10c,ipol%i_vort) =   ratio * d_rj(i10c_ri,ipol%i_vort)
        d_rj(i10c,ipol%i_vort+1) = two *   d_rj(i10c_ri,ipol%i_vort)    &
     &                            * radius_1d_rj_r(k)*ar_1d_rj(kr_in,2)
      end do
!$omp end parallel do
!
      i10c = idx_rj_l0 + (kr_in-1)*jmax
      d_rj(i10c,ipol%i_vort+1) = two *   d_rj(i10c_ri,ipol%i_vort)      &
     &                          * ar_1d_rj(kr_in,1)
!
      end subroutine set_inner_core_rot_l1
!
! ----------------------------------------------------------------------
!
      subroutine cal_icore_viscous_drag_l1(idx_rj_l0, kr_in,            &
     &          fdm1_fix_fld_ICB, coef_d, it_velo, it_viscous,          &
     &          nri, jmax, ar_1d_rj, nnod_rj, ntot_phys_rj, d_rj)
!
      use m_fdm_coefs
!
      integer(kind = kint), intent(in) :: nnod_rj, nri, jmax
      real(kind = kreal), intent(in) :: coef_d
      real(kind = kreal), intent(in) :: fdm1_fix_fld_ICB(0:1,2)
      integer(kind = kint), intent(in) :: kr_in, idx_rj_l0
      integer(kind = kint), intent(in) :: it_velo, it_viscous
      integer(kind = kint), intent(in) :: ntot_phys_rj
      real(kind= kreal), intent(in) :: ar_1d_rj(nri,3)
!
      real(kind = kreal), intent(inout) :: d_rj(nnod_rj,ntot_phys_rj)
!
      integer(kind = kint) ::  i10c_ri, i10c_r1
      real(kind = kreal) :: mat_1, mat_0
!
!
      if(idx_rj_l0 .le. 0) return
!
      i10c_ri = idx_rj_l0 + (kr_in-1)*jmax
      i10c_r1 = idx_rj_l0 +  kr_in * jmax
!
      mat_0 = fdm1_fix_fld_ICB(0,2) - two*ar_1d_rj(kr_in,1)
      mat_1 = fdm1_fix_fld_ICB(1,2)
!
      d_rj(i10c_ri,it_viscous)                                          &
     &                   =  five  * coef_d * ar_1d_rj(kr_in,1)          &
     &                          * (mat_0 * d_rj(i10c_ri,it_velo)        &
     &                           + mat_1 * d_rj(i10c_r1,it_velo))
!
      end subroutine cal_icore_viscous_drag_l1
!
! ----------------------------------------------------------------------
!
      subroutine int_icore_tor_lorentz_l1(idx_rj_l0, kr_in,             &
     &          nri, jmax, radius_1d_rj_r, ar_1d_rj,                    &
     &          nnod_rj, ntot_phys_rj, d_rj)
!
      integer(kind = kint), intent(in) :: nri, jmax
      integer(kind = kint), intent(in) :: kr_in, idx_rj_l0
      integer(kind = kint), intent(in) :: nnod_rj, ntot_phys_rj
      real(kind= kreal), intent(in) :: radius_1d_rj_r(nri)
      real(kind= kreal), intent(in) :: ar_1d_rj(nri,3)
!
      real(kind = kreal), intent(inout) :: d_rj(nnod_rj,ntot_phys_rj)
!
      integer(kind = kint) :: k
      integer(kind = kint) :: i10c_i, i10c_o
      real(kind = kreal) :: sk_10c
!
!
      if(idx_rj_l0 .le. 0) return
!
      i10c_o = idx_rj_l0
      sk_10c = d_rj(i10c_o,itor%i_lorentz) * radius_1d_rj_r(1)**3
!
!$omp parallel do reduction(+:sk_10c) private(i10c_i,i10c_o)
      do k = 1, kr_in-1
        i10c_i = idx_rj_l0 + (k-1)*jmax
        i10c_o = idx_rj_l0 + (k  )*jmax
!
        sk_10c = sk_10c                                                 &
     &        + (d_rj(i10c_i,itor%i_lorentz) * radius_1d_rj_r(k  )**2   &
     &         + d_rj(i10c_o,itor%i_lorentz) * radius_1d_rj_r(k+1)**2)  &
     &        * (radius_1d_rj_r(k+1) - radius_1d_rj_r(k))
      end do
!$omp end parallel do
!
      i10c_o = idx_rj_l0 + (kr_in-1)*jmax
      d_rj(i10c_o,itor%i_lorentz) = half * five * sk_10c                &
     &                           * ar_1d_rj(kr_in,1)**3
      d_rj(i10c_o,ipol%i_rot_Lorentz) = d_rj(i10c_o,itor%i_lorentz)
!
      end subroutine int_icore_tor_lorentz_l1
!
! ----------------------------------------------------------------------
!
      end module cal_inner_core_rotation
