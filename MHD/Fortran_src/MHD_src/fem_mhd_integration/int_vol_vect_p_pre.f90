!
!     module int_vol_vect_p_pre
!
!     numerical integration for finite elememt equations of induction
!
!        programmed by H.Matsui and H.Okuda
!                              on July 2000 (ver 1.1)
!        modified by H. Matsui on Oct., 2005
!        modified by H. Matsui on Aug., 2007
!
!      subroutine int_vol_vect_p_pre_ele
!      subroutine int_vol_vect_p_pre_ele_upm
!
      module int_vol_vect_p_pre
!
      use m_precision
!
      use m_machine_parameter
      use m_control_parameter
      use m_geometry_data
      use m_phys_constants
      use m_geometry_data_MHD
      use m_node_phys_address
      use m_element_phys_address
      use m_element_phys_data
!
      use m_physical_property
!
      implicit none
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine int_vol_vect_p_pre_ele
!
      use m_finite_element_matrix
      use m_int_vol_data
!
      use cal_add_smp
      use nodal_fld_cst_to_ele_1st
      use cal_skv_to_ff_smp_1st
      use fem_skv_nonlinear_1st
!
      integer(kind=kint) :: k2
!
!
      call reset_sk6(n_vector)
!
!   include external magnetic field
!$omp parallel
      call add_const_to_vector_smp(np_smp, ele1%numele, iele_smp_stack, &
     &      d_ele(1,iphys_ele%i_magne), ex_magne, vect_e)
!$omp end parallel
!
! -------- loop for shape function for the phsical values
      do k2=1, nnod_4_ele
        call vector_cst_phys_2_each_ele(k2, iphys%i_velo,               &
     &      coef_induct, velo_1)
!
        call fem_skv_rot_inertia_1st(iele_cd_smp_stack,                 &
     &      intg_point_t_evo, k2, velo_1, vect_e, sk6)
      end do
!
      call sub3_skv_to_ff_v_smp_1st(ff_nl_smp, sk6)
!
      end subroutine int_vol_vect_p_pre_ele
!
!-----------------------------------------------------------------------
!
      subroutine int_vol_vect_p_pre_ele_upm
!
      use m_finite_element_matrix
      use m_int_vol_data
!
      use cal_add_smp
      use nodal_fld_cst_to_ele_1st
      use cal_skv_to_ff_smp_1st
      use fem_skv_nonlinear_upw_1st
!
      integer(kind = kint) :: k2
!
!
      call reset_sk6(n_vector)
!
!$omp parallel
      call add_const_to_vector_smp(np_smp, ele1%numele, iele_smp_stack, &
     &    d_ele(1,iphys_ele%i_magne), ex_magne, vect_e)
!$omp end parallel
!
! -------- loop for shape function for the phsical values
      do k2 = 1, nnod_4_ele
        call vector_cst_phys_2_each_ele(k2, iphys%i_velo,               &
     &      coef_induct, velo_1)
!
        call fem_skv_rot_inertia_upw_1st(iele_cd_smp_stack,             &
     &      intg_point_t_evo, k2, velo_1, vect_e,                       &
     &      d_ele(1,iphys_ele%i_magne), sk6)
      end do
!
      call sub3_skv_to_ff_v_smp_1st(ff_nl_smp, sk6)
!
      end subroutine int_vol_vect_p_pre_ele_upm
!
!-----------------------------------------------------------------------
!
      end module int_vol_vect_p_pre
