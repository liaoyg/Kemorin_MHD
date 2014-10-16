!>@file   copy_MHD_4_sph_trans.f90
!!@brief  module copy_MHD_4_sph_trans
!!
!!@author H. Matsui
!!@date    programmed by H.Matsui in Oct., 2009
!
!>@brief Copy spectrum data and field data to spherical transform buffer
!!       for dynamo simulation
!!
!!@verbatim
!!  routines for backward transform
!!      subroutine copy_mhd_vec_fld_from_trans
!!      subroutine copy_mhd_scl_fld_from_trans
!!
!!  routines for forward transform
!!      subroutine copy_mhd_scalar_spec_from_trans
!!
!!      subroutine copy_mhd_vec_fld_to_trans
!!@endverbatim
!
      module copy_MHD_4_sph_trans
!
      use m_precision
!
      use m_sph_phys_address
      use m_addresses_trans_sph_MHD
!
      implicit  none
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine copy_mhd_vec_fld_from_trans
!
      use copy_sph_field_4_sph_trans
!
!
!$omp parallel
      call copy_vec_fld_from_trans(irtp%i_velo, b_trns%i_velo)
      call copy_vec_fld_from_trans(irtp%i_vort, b_trns%i_vort)
      call copy_vec_fld_from_trans(irtp%i_magne, b_trns%i_magne)
      call copy_vec_fld_from_trans(irtp%i_current, b_trns%i_current)
!$omp end parallel
!
      end subroutine copy_mhd_vec_fld_from_trans
!
!-----------------------------------------------------------------------
!
      subroutine copy_mhd_scl_fld_from_trans
!
      use copy_sph_field_4_sph_trans
!
!
!$omp parallel
      call copy_scalar_fld_from_trans(irtp%i_temp, b_trns%i_temp)
      call copy_scalar_fld_from_trans(irtp%i_light, b_trns%i_light)
!$omp end parallel
!
      end subroutine copy_mhd_scl_fld_from_trans
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine copy_mhd_vec_fld_to_trans
!
      use copy_sph_field_4_sph_trans
!
!
!$omp parallel
!   advection flag
      call copy_vec_fld_to_trans(irtp%i_m_advect, f_trns%i_m_advect)
!   Coriolis flag
      call copy_vec_fld_to_trans(irtp%i_coriolis, f_trns%i_coriolis)
!   Lorentz flag
      call copy_vec_fld_to_trans(irtp%i_lorentz, f_trns%i_lorentz)
!
!   induction flag
      call copy_vec_fld_to_trans(irtp%i_vp_induct, f_trns%i_vp_induct)
!   divergence of heat flux flag
      call copy_vec_fld_to_trans(irtp%i_h_flux, f_trns%i_h_flux)
!
!   divergence of composition flux flag
      call copy_vec_fld_to_trans(irtp%i_c_flux, f_trns%i_c_flux)
!$omp end parallel
!
      end  subroutine copy_mhd_vec_fld_to_trans
!
!-----------------------------------------------------------------------
!
      end module copy_MHD_4_sph_trans
