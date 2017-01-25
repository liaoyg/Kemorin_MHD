!>@file   cal_nonlinear_sph_MHD.f90
!!@brief  module cal_nonlinear_sph_MHD
!!
!!@author H. Matsui
!!@date Programmed in Oct., 2009
!
!>@brief  Evaluate nonlinear terms in spherical coordinate grid
!!
!!@verbatim
!!      subroutine nonlinear_terms_in_rtp                               &
!!     &         (sph_rtp, fl_prop, cd_prop, ht_prop, cp_prop1,         &
!!     &          b_trns, f_trns, ncomp_rj_2_rtp, ncomp_rtp_2_rj,       &
!!     &          fld_rtp, frc_rtp)
!!        type(sph_rtp_grid), intent(in) :: sph_rtp
!!        type(fluid_property), intent(in) :: fl_prop
!!        type(conductive_property), intent(in) :: cd_prop
!!        type(scalar_property), intent(in) :: ht_prop, cp_prop
!!        type(phys_address), intent(in) :: b_trns, f_trns
!!      subroutine add_ref_advect_sph_MHD(sph_rj,                       &
!!     &          ht_prop, cp_prop, ref_param_T, ref_param_C,           &
!!     &          leg, reftemp_rj, refcomp_rj, ipol, rj_fld)
!!        type(sph_rj_grid), intent(in) ::  sph_rj
!!        type(scalar_property), intent(in) :: ht_prop
!!        type(scalar_property), intent(in) :: cp_prop
!!        type(reference_scalar_param), intent(in) :: ref_param_T
!!        type(reference_scalar_param), intent(in) :: ref_param_C
!!        type(legendre_4_sph_trans), intent(in) :: leg
!!        type(phys_address), intent(in) :: ipol
!!        type(phys_data), intent(inout) :: rj_fld
!!      subroutine add_reference_advect_sph                             &
!!     &         (kr_in, kr_out, nidx_rj, ar_1d_rj, g_sph_rj,           &
!!     &          coef_advect, is_h_advect, is_velo,                    &
!!     &          nnod_rj, ntot_phys_rj, reftemp_rj, d_rj)
!!@endverbatim
!!
!!@n @param kr_in       Radial ID for inner boundary
!!@n @param kr_out      Radial ID for outer boundary
!
      module cal_nonlinear_sph_MHD
!
      use m_precision
!
      use m_machine_parameter
      use m_constants
!
      use t_physical_property
      use t_reference_scalar_param
      use t_spheric_rj_data
      use t_spheric_rtp_data
      use t_phys_address
      use t_phys_data
      use t_schmidt_poly_on_rtm
!
      implicit none
!
      private :: add_reference_advect_sph
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine nonlinear_terms_in_rtp                                 &
     &         (sph_rtp, fl_prop, cd_prop, ht_prop, cp_prop,            &
     &          b_trns, f_trns, ncomp_rj_2_rtp, ncomp_rtp_2_rj,         &
     &          fld_rtp, frc_rtp)
!
      use m_control_parameter
      use const_wz_coriolis_rtp
      use cal_products_smp
!
      type(sph_rtp_grid), intent(in) :: sph_rtp
      type(fluid_property), intent(in) :: fl_prop
      type(conductive_property), intent(in) :: cd_prop
      type(scalar_property), intent(in) :: ht_prop, cp_prop
      type(phys_address), intent(in) :: b_trns, f_trns
      integer(kind = kint), intent(in) :: ncomp_rj_2_rtp
      integer(kind = kint), intent(in) :: ncomp_rtp_2_rj
      real(kind = kreal), intent(in)                                    &
     &                   :: fld_rtp(sph_rtp%nnod_rtp,ncomp_rj_2_rtp)
      real(kind = kreal), intent(inout)                                 &
     &                   :: frc_rtp(sph_rtp%nnod_rtp,ncomp_rtp_2_rj)
!
!
!$omp parallel
      if( (f_trns%i_m_advect*evo_velo%iflag_scheme) .gt. 0) then
        call cal_cross_prod_w_coef_smp                                  &
     &     (sph_rtp%nnod_rtp, fl_prop%coef_velo,                        &
     &      fld_rtp(1,b_trns%i_vort), fld_rtp(1,b_trns%i_velo),         &
     &      frc_rtp(1,f_trns%i_m_advect) )
      end if
!
      if( (f_trns%i_lorentz*iflag_4_lorentz) .gt. 0) then
        call cal_cross_prod_w_coef_smp                                  &
     &     (sph_rtp%nnod_rtp, fl_prop%coef_lor,                         &
     &      fld_rtp(1,b_trns%i_current), fld_rtp(1,b_trns%i_magne),     &
     &      frc_rtp(1,f_trns%i_lorentz) )
      end if
!
!
!
      if( (f_trns%i_vp_induct * evo_magne%iflag_scheme) .gt. 0) then
        call cal_cross_prod_w_coef_smp                                  &
     &     (sph_rtp%nnod_rtp, cd_prop%coef_induct,                      &
     &      fld_rtp(1,b_trns%i_velo), fld_rtp(1,b_trns%i_magne),        &
     &      frc_rtp(1,f_trns%i_vp_induct) )
      end if
!
!
      if( (f_trns%i_h_flux * evo_temp%iflag_scheme) .gt. 0) then
        call cal_vec_scalar_prod_w_coef_smp                             &
     &     (sph_rtp%nnod_rtp, ht_prop%coef_advect,                      &
     &      fld_rtp(1,b_trns%i_velo), fld_rtp(1,b_trns%i_temp),         &
     &      frc_rtp(1,f_trns%i_h_flux) )
      end if
!
      if( (f_trns%i_c_flux * evo_comp%iflag_scheme) .gt. 0) then
        call cal_vec_scalar_prod_w_coef_smp                             &
     &     (sph_rtp%nnod_rtp, cp_prop%coef_advect,                      &
     &      fld_rtp(1,b_trns%i_velo), fld_rtp(1,b_trns%i_light),        &
     &      frc_rtp(1,f_trns%i_c_flux) )
      end if
!
!      if( (f_trns%i_Coriolis*iflag_4_coriolis) .gt. 0) then
!        call cal_wz_coriolis_rtp                                       &
!     &     (sph_rtp%nnod_rtp, sph_rtp%nidx_rtp, fl_prop%coef_cor,      &
!     &      fld_rtp(1,b_trns%i_velo), frc_rtp(1,f_trns%i_Coriolis))
!      end if
!$omp end parallel
!
      end subroutine nonlinear_terms_in_rtp
!
!-----------------------------------------------------------------------
!
      subroutine add_ref_advect_sph_MHD(sph_rj,                         &
     &          ht_prop, cp_prop, ref_param_T, ref_param_C,             &
     &          leg, reftemp_rj, refcomp_rj, ipol, rj_fld)
!
      use m_boundary_params_sph_MHD
!
      type(sph_rj_grid), intent(in) ::  sph_rj
      type(scalar_property), intent(in) :: ht_prop
      type(scalar_property), intent(in) :: cp_prop
      type(reference_scalar_param), intent(in) :: ref_param_T
      type(reference_scalar_param), intent(in) :: ref_param_C
      type(legendre_4_sph_trans), intent(in) :: leg
      type(phys_address), intent(in) :: ipol
!
      real(kind = kreal), intent(in)                                    &
     &                   :: reftemp_rj(sph_rj%nidx_rj(1),0:1)
      real(kind = kreal), intent(in)                                    &
     &                   :: refcomp_rj(sph_rj%nidx_rj(1),0:1)
!
      type(phys_data), intent(inout) :: rj_fld
!
!
!   ----  Lead advection of reference field
      if (ref_param_T%iflag_reference .eq. id_sphere_ref_temp) then
        call add_reference_advect_sph(sph_bc_T%kr_in, sph_bc_T%kr_out,  &
     &      sph_rj%nidx_rj, sph_rj%ar_1d_rj, leg%g_sph_rj,              &
     &      ht_prop%coef_advect, ipol%i_h_advect, ipol%i_velo,          &
     &      rj_fld%n_point, rj_fld%ntot_phys, reftemp_rj, rj_fld%d_fld)
      end if
      if (ref_param_C%iflag_reference .eq. id_sphere_ref_temp) then
        call add_reference_advect_sph(sph_bc_C%kr_in, sph_bc_T%kr_out,  &
     &      sph_rj%nidx_rj, sph_rj%ar_1d_rj, leg%g_sph_rj,              &
     &      cp_prop%coef_advect, ipol%i_c_advect, ipol%i_velo,          &
     &      rj_fld%n_point, rj_fld%ntot_phys, refcomp_rj, rj_fld%d_fld)
      end if
!
      end subroutine add_ref_advect_sph_MHD
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine add_reference_advect_sph                               &
     &         (kr_in, kr_out, nidx_rj, ar_1d_rj, g_sph_rj,             &
     &          coef_advect, is_h_advect, is_velo,                      &
     &          nnod_rj, ntot_phys_rj, reftemp_rj, d_rj)
!
      integer(kind = kint), intent(in) :: kr_in, kr_out
      integer(kind = kint), intent(in) :: nidx_rj(2)
      integer(kind = kint), intent(in) :: is_h_advect, is_velo
      integer(kind = kint), intent(in) :: nnod_rj, ntot_phys_rj
      real(kind = kreal), intent(in) :: coef_advect
      real(kind = kreal), intent(in) :: g_sph_rj(nidx_rj(2),13)
      real(kind = kreal), intent(in) :: ar_1d_rj(nidx_rj(1),3)
      real(kind = kreal), intent(in) :: reftemp_rj(nidx_rj(1),0:1)
      real(kind = kreal), intent(inout) :: d_rj(nnod_rj,ntot_phys_rj)
!
      integer(kind= kint) :: ist, ied, inod, j, k
!
!
      ist = (kr_in-1)*nidx_rj(2) + 1
      ied = kr_out * nidx_rj(2)
!$omp parallel do private (inod,j,k)
      do inod = ist, ied
        j = mod((inod-1),nidx_rj(2)) + 1
        k = 1 + (inod- j) / nidx_rj(2)
!
        d_rj(inod,is_h_advect) = d_rj(inod,is_h_advect)                 &
     &                    + coef_advect * g_sph_rj(j,3) * ar_1d_rj(k,2) &
     &                     * reftemp_rj(k,1) * d_rj(inod,is_velo)
      end do
!$omp end parallel do
!
      end subroutine add_reference_advect_sph
!
!-----------------------------------------------------------------------
!
      end module cal_nonlinear_sph_MHD
