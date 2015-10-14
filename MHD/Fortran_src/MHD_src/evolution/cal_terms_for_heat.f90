!
!     module cal_terms_for_heat
!
!     Written by H. Matsui on June, 2005
!
!      subroutine cal_terms_4_heat(i_field)
!      subroutine cal_viscous_diffusion
!
      module cal_terms_for_heat
!
      use m_precision
!
      use m_control_parameter
      use m_phys_constants
      use m_finite_element_matrix
!
      use cal_multi_pass
      use cal_ff_smp_to_ffs
      use cal_for_ffs
      use nod_phys_send_recv
      use set_boundary_scalars
      use check_finite_element_mat
      use int_surf_temp
!
      implicit none
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine cal_terms_4_heat(i_field)
!
      use m_geometry_data
      use m_node_phys_address
      use m_node_phys_data
      use m_int_vol_data
!
      use int_vol_temp_monitor
!
      integer (kind = kint), intent(in) :: i_field
!
!
      call reset_ff_smps
!
      if (iflag_temp_supg .gt. id_turn_OFF) then
       call int_vol_ene_monitor_upw(i_field)
      else
       call int_vol_ene_monitor(i_field)
      end if
!
      call int_surf_temp_monitor(i_field)
!
      call cal_t_evo_4_scalar_fl(iflag_temp_supg)
!
      call set_boundary_ene_4_rhs
!
!       call check_ff_nl(n_scalar)
!
      call cal_ff_2_scalar(node1%numnod, node1%istack_nod_smp,          &
     &    ff_nl, mhd_fem1_wk%ml_fl,                                     &
     &    nod_fld1%ntot_phys, i_field, nod_fld1%d_fld)
!
!   communication
!
      call scalar_send_recv                                             &
     &   (nod_fld1%ntot_phys, i_field, nod_fld1%d_fld)
!
      end subroutine cal_terms_4_heat
!
!-----------------------------------------------------------------------
!
      subroutine cal_thermal_diffusion
!
      use m_geometry_data
      use m_node_phys_address
      use m_node_phys_data
      use m_int_vol_data
!
      use int_vol_diffusion_ele
!
!
      call reset_ff_smps
!
      call int_vol_t_diffuse_ele_monitor
!
      call int_surf_temp_monitor(iphys%i_t_diffuse)
!
      call set_ff_nl_smp_2_ff(node1, rhs_tbl1, n_scalar)
!
      call set_boundary_ene_4_rhs
!
      call cal_ff_2_scalar(node1%numnod, node1%istack_nod_smp,          &
     &    f1_l%ff, mhd_fem1_wk%ml_fl,                                   &
     &    nod_fld1%ntot_phys, iphys%i_t_diffuse, nod_fld1%d_fld)
!
!   communication
!
      call scalar_send_recv                                             &
     &   (nod_fld1%ntot_phys, iphys%i_t_diffuse, nod_fld1%d_fld)
!
      end subroutine cal_thermal_diffusion
!
!-----------------------------------------------------------------------
!
      end module cal_terms_for_heat
