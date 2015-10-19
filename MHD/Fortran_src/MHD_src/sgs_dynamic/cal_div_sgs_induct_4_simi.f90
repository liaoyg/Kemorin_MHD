!
!     module cal_div_sgs_induct_4_simi
!
!     Written by H. Matsui
!
!      subroutine cal_div_sgs_induct_simi
!      subroutine cal_div_sgs_filter_idct_simi
!
      module cal_div_sgs_induct_4_simi
!
      use m_precision
!
      implicit none
!
      private :: cal_div_sgs_idct_simi
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine cal_div_sgs_induct_simi
!
      use m_node_phys_address
!
       call cal_div_sgs_idct_simi(iphys%i_sgs_grad,                     &
     &     iphys%i_SGS_induct_t, iphys%i_velo, iphys%i_magne)
!
      end subroutine cal_div_sgs_induct_simi
!
!-----------------------------------------------------------------------
!
      subroutine cal_div_sgs_filter_idct_simi
!
      use m_node_phys_address
!
       call cal_div_sgs_idct_simi(iphys%i_sgs_simi, iphys%i_sgs_grad_f, &
     &     iphys%i_filter_velo, iphys%i_filter_magne)
!
      end subroutine cal_div_sgs_filter_idct_simi
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine cal_div_sgs_idct_simi(i_sgs, i_flux, i_v, i_b)
!
      use m_control_parameter
      use m_nod_comm_table
      use m_geometry_data
      use m_phys_constants
      use m_finite_element_matrix
      use m_int_vol_data
      use m_sorted_node
      use m_node_phys_data
      use m_element_phys_data
!
      use cal_ff_smp_to_ffs
      use cal_for_ffs
      use nod_phys_send_recv
      use int_div_sgs_idct_simi
!
      integer(kind = kint), intent(in) :: i_flux, i_v, i_b
      integer(kind = kint), intent(in) :: i_sgs
!
!
      call reset_ff_smps(node1%max_nod_smp, f1_l, f1_nl)
!
      if ( iflag_mag_supg .gt. id_turn_OFF) then
        call int_div_sgs_idct_simi_upw(i_flux, i_v, i_b,                &
     &      fld_ele1%ntot_phys, iphys_ele%i_velo, fld_ele1%d_fld)
      else
        call int_div_sgs_idct_simi_pg(i_flux, i_v, i_b)
      end if
!
      call set_ff_nl_smp_2_ff(node1, rhs_tbl1, n_vector)
      call cal_ff_2_vector(node1%numnod, node1%istack_nod_smp,          &
     &    f1_nl%ff, mhd_fem1_wk%mlump_fl%ml, nod_fld1%ntot_phys,        &
     &    i_sgs, nod_fld1%d_fld)
!
! ----------   communications
!
      call vector_send_recv(i_sgs, node1, nod_comm, nod_fld1)
!
      end subroutine cal_div_sgs_idct_simi
!
!-----------------------------------------------------------------------
!
      end module cal_div_sgs_induct_4_simi
