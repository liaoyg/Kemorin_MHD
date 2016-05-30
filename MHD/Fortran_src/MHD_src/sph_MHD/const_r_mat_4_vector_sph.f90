!>@file   const_r_mat_4_vector_sph.f90
!!@brief  module const_r_mat_4_vector_sph
!!
!!@date  Programmed by H.Matsui on Apr., 2009
!
!>@brief Construct matrix for time evolution of vector fields
!!
!!@verbatim
!!      subroutine const_radial_mat_vort_2step(sph_rj, band_vs_poisson, &
!!     &          band_vp_evo, band_vt_evo, band_wt_evo)
!!      subroutine const_radial_mat_4_magne_sph                         &
!!     &         (sph_rj, band_bp_evo, band_bt_evo)
!!        type(sph_rj_grid), intent(in) :: sph_rj
!!@endverbatim
!
      module const_r_mat_4_vector_sph
!
      use m_precision
      use calypso_mpi
!
      use m_constants
      use m_machine_parameter
      use m_t_int_parameter
      use m_physical_property
      use m_ludcmp_3band
!
      use t_spheric_rj_data
      use t_sph_matrices
!
      use set_radial_mat_sph
!
      implicit none
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine const_radial_mat_vort_2step(sph_rj, band_vs_poisson,   &
     &          band_vp_evo, band_vt_evo, band_wt_evo)
!
      use m_boundary_params_sph_MHD
      use m_coef_fdm_to_center
      use m_coef_fdm_free_ICB
      use m_coef_fdm_free_CMB
      use m_ludcmp_band
      use set_sph_scalar_mat_bc
      use cal_inner_core_rotation
      use center_sph_matrices
      use mat_product_3band_mul
      use check_sph_radial_mat
!
      type(sph_rj_grid), intent(in) :: sph_rj
      type(band_matrices_type), intent(inout) :: band_vp_evo
      type(band_matrices_type), intent(inout) :: band_vt_evo
      type(band_matrices_type), intent(inout) :: band_wt_evo
      type(band_matrices_type), intent(inout) :: band_vs_poisson
!
!      integer(kind = kint) :: j
      real(kind = kreal) :: coef_dvt
!
!
      call alloc_band_mat_sph(ifive, sph_rj, band_vp_evo)
      call alloc_band_mat_sph(ithree, sph_rj, band_vt_evo)
      call alloc_band_mat_sph(ithree, sph_rj, band_wt_evo)
      call alloc_band_mat_sph(ithree, sph_rj, band_vs_poisson)
!
      call set_unit_on_diag(band_vp_evo)
      call set_unit_on_diag(band_vt_evo)
      call set_unit_on_diag(band_wt_evo)
!
      if(coef_d_velo .eq. zero) then
        coef_dvt = one
        call set_unit_mat_4_poisson                                     &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2),                       &
     &      sph_bc_U%kr_in, sph_bc_U%kr_out, band_vt_evo%mat)
        call set_unit_mat_4_poisson                                     &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2),                       &
     &      sph_bc_U%kr_in, sph_bc_U%kr_out, band_wt_evo%mat)
      else
        coef_dvt = coef_imp_v * coef_d_velo * dt
        call set_unit_mat_4_time_evo                                    &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2), band_vt_evo%mat)
        call set_unit_mat_4_time_evo                                    &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2), band_wt_evo%mat)
      end if
!
      call set_unit_mat_4_poisson                                       &
     &   (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2),                         &
     &    sph_bc_U%kr_in, sph_bc_U%kr_out, band_vs_poisson%mat)
!
      call add_vector_poisson_mat_sph                                   &
     &   (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2), sph_rj%ar_1d_rj,        &
     &    sph_bc_U%kr_in, sph_bc_U%kr_out, coef_dvt, band_vt_evo%mat)
      call add_vector_poisson_mat_sph                                   &
     &   (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2), sph_rj%ar_1d_rj,        &
     &    sph_bc_U%kr_in, sph_bc_U%kr_out, coef_dvt, band_wt_evo%mat)
      call add_vector_poisson_mat_sph                                   &
     &   (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2), sph_rj%ar_1d_rj,        &
     &    sph_bc_U%kr_in, sph_bc_U%kr_out, one, band_vs_poisson%mat)
!
!   Boundary condition for ICB
!
      if(sph_bc_U%iflag_icb .eq. iflag_sph_fill_center) then
        call add_vector_poisson_mat_center                              &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2), sph_bc_U%r_ICB,       &
     &      fdm2_fix_fld_ctr1, coef_dvt, band_vt_evo%mat)
        call add_vector_poisson_mat_center                              &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2), sph_bc_U%r_ICB,       &
     &      fdm2_fix_fld_ctr1, coef_dvt, band_wt_evo%mat)
        call add_vector_poisson_mat_center                              &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2), sph_bc_U%r_ICB,       &
     &      fdm2_fix_fld_ctr1, one, band_vs_poisson%mat)
      else
        call add_fix_flux_icb_poisson_mat                               &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2),                       &
     &      sph_bc_U%kr_in, sph_bc_U%r_ICB, sph_bc_U%fdm2_fix_dr_ICB,   &
     &      coef_dvt, band_wt_evo%mat)
!
        if(sph_bc_U%iflag_icb .eq. iflag_free_slip) then
          call add_fix_flux_icb_poisson_mat                             &
     &       (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2),                     &
     &        sph_bc_U%kr_in, sph_bc_U%r_ICB, fdm2_free_vt_ICB,         &
     &        coef_dvt, band_vt_evo%mat)
          call add_fix_flux_icb_poisson_mat                             &
     &       (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2),                     &
     &        sph_bc_U%kr_in, sph_bc_U%r_ICB, fdm2_free_vp_ICB,         &
     &        one, band_vs_poisson%mat)
        else
          call set_fix_fld_icb_poisson_mat                              &
     &       (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2),                     &
     &        sph_bc_U%kr_in, band_vt_evo%mat)
          call add_fix_flux_icb_poisson_mat                             &
     &       (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2),                     &
     &        sph_bc_U%kr_in, sph_bc_U%r_ICB, sph_bc_U%fdm2_fix_dr_ICB, &
     &        one, band_vs_poisson%mat)
        end if
      end if
!
!   Overwrite rotation for inner core
!
      if(sph_bc_U%iflag_icb .eq. iflag_rotatable_ic) then
        call set_icore_viscous_matrix                                   &
     &     (sph_bc_U%kr_in, sph_bc_U%fdm1_fix_fld_ICB,                  &
     &      sph_rj, band_vt_evo)
      end if
!
!   Boundary condition for CMB
!
      call add_fix_flux_cmb_poisson_mat                                 &
     &   (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2),                         &
     &    sph_bc_U%kr_out, sph_bc_U%r_CMB, sph_bc_U%fdm2_fix_dr_CMB,    &
     &    coef_dvt, band_wt_evo%mat)
!
      if(sph_bc_U%iflag_cmb .eq. iflag_free_slip) then
        call add_fix_flux_cmb_poisson_mat                               &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2),                       &
     &      sph_bc_U%kr_out, sph_bc_U%r_CMB, fdm2_free_vt_CMB,          &
     &      coef_dvt, band_vt_evo%mat)
        call add_fix_flux_cmb_poisson_mat                               &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2),                       &
     &      sph_bc_U%kr_out, sph_bc_U%r_CMB, fdm2_free_vp_CMB,          &
     &      one, band_vs_poisson%mat)
      else
        call set_fix_fld_cmb_poisson_mat                                &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2),                       &
     &      sph_bc_U%kr_out, band_vt_evo%mat)
        call add_fix_flux_cmb_poisson_mat                               &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2),                       &
     &      sph_bc_U%kr_out, sph_bc_U%r_CMB, sph_bc_U%fdm2_fix_dr_CMB,  &
     &      one, band_vs_poisson%mat)
      end if
!
!
      call cal_mat_product_3band_mul                                    &
     &   (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2),                         &
     &    sph_bc_U%kr_in, sph_bc_U%kr_out, band_wt_evo%mat,             &
     &    band_vs_poisson%mat, band_vp_evo%mat)
!
      call ludcmp_5band_mul_t                                           &
     &   (np_smp, sph_rj%istack_rj_j_smp, band_vp_evo)
      call ludcmp_3band_mul_t                                           &
     &   (np_smp, sph_rj%istack_rj_j_smp, band_vt_evo)
      call ludcmp_3band_mul_t                                           &
     &   (np_smp, sph_rj%istack_rj_j_smp, band_wt_evo)
      call ludcmp_3band_mul_t                                           &
     &   (np_smp, sph_rj%istack_rj_j_smp, band_vs_poisson)
!
      if(i_debug .eq. iflag_full_msg) then
        write(50+my_rank,'(a)') 'evolution matrix for toroidal velocity'
        call check_radial_band_mat(my_rank, sph_rj, band_vt_evo)
!
        write(50+my_rank,'(a)') 'evolution matrix for toroidal vorticity'
        call check_radial_band_mat(my_rank, sph_rj, band_wt_evo)
!
        write(50+my_rank,'(a)') 'evolution matrix for poloidal velocity'
        call check_radial_band_mat(my_rank, sph_rj, band_vp_evo)
      end if
!
!      do j = 1, sph_rj%nidx_rj(2)
!        do k = 1, sph_rj%nidx_rj(1)
!          band_vp_evo%det(j)                                           &
!     &                = band_vp_evo%det(j) * band_vp_evo%lu(5,k,j)
!          band_vt_evo%det(j)                                           &
!     &                = band_vt_evo%det(j) * band_vt_evo%lu(3,k,j)
!        end do
!        write(my_rank+60,*) 'det vp', j,                               &
!     &                       band_vp_evo%det(j), band_vt_evo%det(j)
!      end do
!
      end subroutine const_radial_mat_vort_2step
!
! -----------------------------------------------------------------------
!
      subroutine const_radial_mat_4_magne_sph                           &
     &         (sph_rj, band_bp_evo, band_bt_evo)
!
      use m_boundary_params_sph_MHD
      use m_coef_fdm_to_center
      use set_sph_scalar_mat_bc
      use set_sph_magne_mat_bc
      use center_sph_matrices
      use check_sph_radial_mat
!
      type(sph_rj_grid), intent(in) :: sph_rj
      type(band_matrices_type), intent(inout) :: band_bp_evo
      type(band_matrices_type), intent(inout) :: band_bt_evo
!
      real(kind = kreal) :: coef_dbt
!
!
      call alloc_band_mat_sph(ithree, sph_rj, band_bp_evo)
      call alloc_band_mat_sph(ithree, sph_rj, band_bt_evo)
!
      call set_unit_on_diag(band_bp_evo)
      call set_unit_on_diag(band_bt_evo)
!
      if(coef_d_magne .eq. zero) then
        coef_dbt = one
        call set_unit_mat_4_poisson                                     &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2),                       &
     &      sph_bc_B%kr_in, sph_bc_B%kr_out, band_bp_evo%mat)
        call set_unit_mat_4_poisson                                     &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2),                       &
     &      sph_bc_B%kr_in, sph_bc_B%kr_out, band_bt_evo%mat)
      else
        coef_dbt = coef_imp_b * coef_d_magne * dt
        call set_unit_mat_4_time_evo                                    &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2), band_bp_evo%mat)
        call set_unit_mat_4_time_evo                                    &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2), band_bt_evo%mat)
      end if
!
      call add_vector_poisson_mat_sph                                   &
     &   (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2), sph_rj%ar_1d_rj,        &
     &    sph_bc_B%kr_in, sph_bc_B%kr_out, coef_dbt, band_bp_evo%mat)
      call add_vector_poisson_mat_sph                                   &
     &   (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2), sph_rj%ar_1d_rj,        &
     &    sph_bc_B%kr_in, sph_bc_B%kr_out, coef_dbt, band_bt_evo%mat)
!
!
      if(sph_bc_B%iflag_icb .eq. iflag_sph_fill_center) then
        call add_vector_poisson_mat_center                              &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2), sph_bc_B%r_ICB,       &
     &      fdm2_fix_fld_ctr1, coef_dbt, band_bp_evo%mat)
        call add_vector_poisson_mat_center                              &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2), sph_bc_B%r_ICB,       &
     &      fdm2_fix_fld_ctr1, coef_dbt, band_bt_evo%mat)
      else if(sph_bc_B%iflag_icb .eq. iflag_radial_magne) then
        call add_fix_flux_icb_poisson_mat                               &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2),                       &
     &      sph_bc_B%kr_in, sph_bc_B%r_ICB, sph_bc_B%fdm2_fix_dr_ICB,   &
     &      coef_dbt, band_bp_evo%mat)
        call set_fix_fld_icb_poisson_mat                                &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2),                       &
     &      sph_bc_B%kr_in, band_bt_evo%mat)
      else
        call set_ins_magne_icb_rmat_sph                                 &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2),                       &
     &      sph_bc_B%kr_in, sph_bc_B%r_ICB, sph_bc_B%fdm2_fix_dr_ICB,   &
     &      coef_dbt, band_bp_evo%mat)
        call set_fix_fld_icb_poisson_mat                                &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2),                       &
     &      sph_bc_B%kr_in, band_bt_evo%mat)
      end if
!
      if(sph_bc_B%iflag_cmb .eq. iflag_radial_magne) then
        call add_fix_flux_cmb_poisson_mat                               &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2),                       &
     &      sph_bc_B%kr_out, sph_bc_B%r_CMB, sph_bc_B%fdm2_fix_dr_CMB,  &
     &      coef_dbt, band_bp_evo%mat)
      else
        call set_ins_magne_cmb_rmat_sph                                 &
     &     (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2),                       &
     &      sph_bc_B%kr_out, sph_bc_B%r_CMB, sph_bc_B%fdm2_fix_dr_CMB,  &
     &      coef_dbt, band_bp_evo%mat)
      end if
      call set_fix_fld_cmb_poisson_mat                                  &
     &   (sph_rj%nidx_rj(1), sph_rj%nidx_rj(2),                         &
     &    sph_bc_B%kr_out, band_bt_evo%mat)
!
!
      call ludcmp_3band_mul_t                                           &
     &   (np_smp, sph_rj%istack_rj_j_smp, band_bp_evo)
      call ludcmp_3band_mul_t                                           &
     &   (np_smp, sph_rj%istack_rj_j_smp, band_bt_evo)
!
      if(i_debug .eq. iflag_full_msg) then
        write(50+my_rank,'(a)') 'evolution matrix for poloidal magne'
        call check_radial_band_mat(my_rank, sph_rj, band_bp_evo)
!
        write(50+my_rank,'(a)') 'evolution matrix for toroidal magne'
        call check_radial_band_mat(my_rank, sph_rj, band_bt_evo)
      end if
!
      end subroutine const_radial_mat_4_magne_sph
!
! -----------------------------------------------------------------------
!
      end module const_r_mat_4_vector_sph
