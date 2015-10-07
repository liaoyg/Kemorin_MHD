!
!     module cal_rotation
!
!     Written by H. Matsui
!
!      subroutine cal_rotation_whole(iflag_4_supg, i_res, i_vector)
!      subroutine cal_rotation_in_fluid(iflag_4_supg, i_res, i_vector)
!      subroutine cal_rotation_in_conduct(iflag_4_supg,                 &
!     &          i_res, i_vector)
!
!      private :: choose_int_vol_rotations
!
      module cal_rotation
!
      use m_precision
!
      use m_phys_constants
      use m_finite_element_matrix
      use m_int_vol_data
!
      use cal_ff_smp_to_ffs
      use cal_for_ffs
      use nod_phys_send_recv
!
      implicit none
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine cal_rotation_whole(iflag_4_supg, i_res, i_vector)
!
      use m_geometry_data
      use m_node_phys_data
      use m_sorted_node
!
      integer(kind = kint), intent(in) :: iflag_4_supg
      integer(kind = kint), intent(in) :: i_vector, i_res
!
!
      call reset_ff_smps
!
      call choose_int_vol_rotations(iflag_4_supg,                       &
     &     ele1%istack_ele_smp, i_vector)
!
      call set_ff_nl_smp_2_ff(node1, rhs_tbl1, n_vector)
      call cal_ff_2_vector(node1%numnod, node1%istack_nod_smp,          &
     &    ff_nl, ml, nod_fld1%ntot_phys, i_res, nod_fld1%d_fld)
!
! ----------   communications
!
      call vector_send_recv(nod_fld1%ntot_phys, i_res, nod_fld1%d_fld)
!
      end subroutine cal_rotation_whole
!
!-----------------------------------------------------------------------
!
      subroutine cal_rotation_in_fluid(iflag_4_supg, i_res, i_vector)
!
      use m_geometry_data
      use m_geometry_data_MHD
      use m_node_phys_data
      use m_sorted_node
!
      integer(kind = kint), intent(in) :: iflag_4_supg
      integer(kind = kint), intent(in) :: i_vector, i_res
!
!
      call reset_ff_smps
!
      call choose_int_vol_rotations(iflag_4_supg,                       &
     &    iele_fl_smp_stack, i_vector)
!
      call set_ff_nl_smp_2_ff(node1, rhs_tbl1, n_vector)
      call cal_ff_2_vector(node1%numnod, node1%istack_nod_smp,          &
     &    ff_nl, ml_fl, nod_fld1%ntot_phys, i_res, nod_fld1%d_fld)
!
! ----------   communications
!
      call vector_send_recv(nod_fld1%ntot_phys, i_res, nod_fld1%d_fld)
!
      end subroutine cal_rotation_in_fluid
!
!-----------------------------------------------------------------------
!------- subroutine cal_rotation_in_conduct ---------------------
!
      subroutine cal_rotation_in_conduct(iflag_4_supg, i_res, i_vector)
!
      use m_geometry_data_MHD
      use m_node_phys_data
      use m_sorted_node
!
      integer(kind = kint), intent(in) :: iflag_4_supg
      integer(kind = kint), intent(in) :: i_vector, i_res
!
!
       call reset_ff_smps
!
       call choose_int_vol_rotations(iflag_4_supg,                      &
     &     iele_cd_smp_stack, i_vector)
!
       call cal_ff_smp_2_vector(node1, rhs_tbl1,                        &
     &     ff_nl_smp, ml_cd, nod_fld1%ntot_phys, i_res, nod_fld1%d_fld)
!
! ----------   communications
!
      call vector_send_recv(nod_fld1%ntot_phys, i_res, nod_fld1%d_fld)
!
      end subroutine cal_rotation_in_conduct
!
!-----------------------------------------------------------------------
!
      subroutine choose_int_vol_rotations(iflag_4_supg,                 &
     &          iele_fsmp_stack, i_vector)
!
      use m_control_parameter
      use m_element_phys_data
!
      use int_vol_vect_diff_1st
      use int_vol_vect_diff_upw_1st
!
       integer(kind = kint), intent(in) :: iflag_4_supg, i_vector
       integer(kind = kint), intent(in) :: iele_fsmp_stack(0:np_smp)
!
       if ( iflag_4_supg .eq. id_magnetic_SUPG) then
        call int_vol_rot_upw_1st                                        &
     &     (iele_fsmp_stack, intg_point_t_evo, i_vector,                &
     &      fld_ele1%ntot_phys, iphys_ele%i_magne, fld_ele1%d_fld)
       else if ( iflag_4_supg .eq. id_turn_ON) then
        call int_vol_rot_upw_1st                                        &
     &     (iele_fsmp_stack, intg_point_t_evo, i_vector,                &
     &      fld_ele1%ntot_phys, iphys_ele%i_velo, fld_ele1%d_fld)
       else
        call int_vol_rot_1st(iele_fsmp_stack, intg_point_t_evo,         &
     &      i_vector)
       end if
!
      end subroutine choose_int_vol_rotations
!
!-----------------------------------------------------------------------
!
      end module cal_rotation
