!
!     module int_sgs_induction
!
!        programmed by H.Matsui on July, 2005
!        modified by H.Matsui on AUg., 2007
!
!      subroutine int_vol_sgs_induction
!
      module int_sgs_induction
!
      use m_precision
!
      implicit none
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine int_vol_sgs_induction
!
      use m_control_parameter
      use m_geometry_parameter
      use m_machine_parameter
      use m_geometry_data_MHD
      use m_geometry_data
      use m_phys_constants
      use m_node_phys_address
      use m_node_phys_data
      use m_finite_element_matrix
!
      use int_vol_vect_diff_1st
      use cal_ff_smp_to_ffs
      use nod_phys_send_recv
!
!
      ff_nl_smp = 0.0d0
!
      call int_vol_rot_1st(iele_cd_smp_stack, intg_point_t_evo,         &
          iphys%i_SGS_vp_induct)
!
!      call cal_multi_pass_4_vector_ff
!      call cal_ff_2_vector(node1%numnod, inod_smp_stack,               &
!    &     d_nod(1,iphys%i_magne), ff, ml_cd)
       call cal_ff_smp_2_vector(d_nod(1,iphys%i_SGS_induction),         &
     &    ff_nl_smp, ml_cd)
!
       call vector_send_recv(iphys%i_SGS_induction)
!
      end subroutine int_vol_sgs_induction
!
!-----------------------------------------------------------------------
!
      end module int_sgs_induction
