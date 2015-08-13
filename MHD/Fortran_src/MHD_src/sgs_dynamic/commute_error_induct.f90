!commute_error_induct.g90
!     module commute_error_induct
!
!     Written by H. Matsui
!
!      subroutine cal_commute_error_4_mf(i_filter, i_sgs, i_flux,       &
!     &          i_vect)
!         i_filter: ID for filter function
!         i_sgs: field ID for obtained difference term
!         i_flux: field ID for SGS term
!         i_vect: field ID for origianl vector field
!         i_scalar: field ID for origianl scalar field
!
!       subroutine cal_commute_error_4_m_flux(i_filter)
!       subroutine cal_commute_error_4_filter_idct(i_filter)
!
      module commute_error_induct
!
      use m_precision
!
      implicit none
!
      private :: cal_commute_error_4_idct
!
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine cal_commute_error_4_idct(i_filter, i_sgs,              &
     &          i_flux, i_v, i_b)
!
      use m_control_parameter
      use m_geometry_data
      use m_geometry_data_MHD
      use m_group_data
      use m_phys_constants
      use m_node_phys_address
      use m_node_phys_data
      use m_finite_element_matrix
      use m_int_vol_data
!
      use cal_ff_smp_to_ffs
      use cal_for_ffs
      use int_vol_commute_1st
      use int_surf_div_induct_tsr_sgs
!
       integer(kind = kint), intent(in) :: i_flux, i_v, i_b
       integer(kind = kint), intent(in) :: i_sgs, i_filter
!
!
      call reset_ff_smps
!
      call int_vol_commute_induct_t(iele_cd_smp_stack,                  &
     &    intg_point_t_evo, i_filter, i_flux, i_v, i_b)
!
      call int_surf_commute_induct_t(sf_grp1, intg_point_t_evo, i_flux, &
     &    i_filter, i_v, i_b)
!
      call set_ff_nl_smp_2_ff(n_vector)
      call cal_ff_2_vector(node1%numnod, node1%istack_nod_smp,          &
     &    d_nod(1,i_sgs), ff_nl, ml_fl)
!
      end subroutine cal_commute_error_4_idct
!
!-----------------------------------------------------------------------
!
      subroutine cal_commute_error_4_induct(i_filter)
!
      use m_node_phys_address
!
      integer(kind = kint), intent(in) :: i_filter
!
       call cal_commute_error_4_idct(i_filter, iphys%i_sgs_grad,        &
     &     iphys%i_SGS_induct_t, iphys%i_velo, iphys%i_magne)
!
      end subroutine cal_commute_error_4_induct
!
!-----------------------------------------------------------------------
!
      subroutine cal_commute_error_4_filter_idct(i_filter)
!
      use m_node_phys_address
!
      integer(kind = kint), intent(in) :: i_filter
!
       call cal_commute_error_4_idct(i_filter, iphys%i_sgs_grad_f,      &
     &     iphys%i_sgs_grad_f, iphys%i_filter_velo,                     &
     &     iphys%i_filter_magne)
!
      end subroutine cal_commute_error_4_filter_idct
!
!-----------------------------------------------------------------------
!
      end module commute_error_induct
