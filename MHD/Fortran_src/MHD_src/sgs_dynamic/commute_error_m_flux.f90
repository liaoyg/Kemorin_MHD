!
!     module commute_error_m_flux
!
!     Written by H. Matsui
!
!       subroutine cal_commute_error_4_m_flux(i_filter)
!       subroutine cal_commute_error_4_filter_mf(i_filter)
!
!         i_filter: ID for filter function
!         i_sgs: field ID for obtained difference term
!         i_flux: field ID for SGS term
!         i_vect: field ID for origianl vector field
!
      module commute_error_m_flux
!
      use m_precision
!
      implicit none
!
      private :: cal_commute_error_4_mf
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine cal_commute_error_4_mf(i_filter, i_sgs, i_flux,        &
     &          i_vect)
!
      use m_control_parameter
      use m_geometry_data
      use m_geometry_data_MHD
      use m_group_data
      use m_phys_constants
      use m_node_phys_data
      use m_finite_element_matrix
      use m_int_vol_data
!
      use cal_ff_smp_to_ffs
      use cal_for_ffs
      use int_vol_commute_1st
      use int_surf_div_fluxes_sgs
!
      integer(kind = kint), intent(in) :: i_flux, i_vect
      integer(kind = kint), intent(in) :: i_sgs, i_filter
!
!
      call reset_ff_smps
!
      call int_vol_commute_div_m_flux(iele_fl_smp_stack,                &
     &    intg_point_t_evo, i_filter, i_flux, i_vect)
!
      call int_surf_commute_sgs_m_flux(sf_grp1, intg_point_t_evo,       &
     &    i_filter, i_flux, i_vect)
!
      call set_ff_nl_smp_2_ff(node1, rhs_tbl1, n_vector)
      call cal_ff_2_vector(node1%numnod, node1%istack_nod_smp,          &
     &    ff_nl, ml_fl, nod_fld1%ntot_phys, i_sgs, nod_fld1%d_fld)
!
      end subroutine cal_commute_error_4_mf
!
!-----------------------------------------------------------------------
!
      subroutine cal_commute_error_4_m_flux(i_filter)
!
      use m_node_phys_address
!
      integer(kind = kint), intent(in) :: i_filter
!
       call cal_commute_error_4_mf(i_filter, iphys%i_sgs_grad,          &
     &     iphys%i_SGS_m_flux, iphys%i_velo)
!
      end subroutine cal_commute_error_4_m_flux
!
!-----------------------------------------------------------------------
!
      subroutine cal_commute_error_4_filter_mf(i_filter)
!
      use m_node_phys_address
!
      integer(kind = kint), intent(in) :: i_filter
!
       call cal_commute_error_4_mf(i_filter, iphys%i_sgs_grad_f,        &
     &     iphys%i_sgs_grad_f, iphys%i_filter_velo)
!
       end subroutine cal_commute_error_4_filter_mf
!
!-----------------------------------------------------------------------
!
      end module commute_error_m_flux
