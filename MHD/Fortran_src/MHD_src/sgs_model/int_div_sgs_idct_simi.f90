!
!      module int_div_sgs_idct_simi
!
!     Written by H. Matsui on July, 2005
!     modified by H. Matsui on July, 2007
!
!      subroutine int_div_sgs_idct_simi_pg(i_flux, i_v, i_b)
!      subroutine int_div_sgs_idct_simi_upw(i_flux, i_v, i_b)
!
      module int_div_sgs_idct_simi
!
      use m_precision
!
      use m_machine_parameter
      use m_control_parameter
      use m_geometry_data
      use m_geometry_data_MHD
      use m_phys_constants
      use m_finite_element_matrix
      use m_int_vol_data
!
      implicit none
!
! ----------------------------------------------------------------------
!
       contains
!
! ----------------------------------------------------------------------
!
      subroutine int_div_sgs_idct_simi_pg(i_flux, i_v, i_b)
!
      use sgs_terms_2_each_ele
      use fem_skv_vector_diff_1st
      use cal_skv_to_ff_smp_1st
!
      integer(kind = kint), intent(in) :: i_flux, i_v, i_b
!
      integer(kind=kint) :: k2
!
!
      call reset_sk6(n_vector)
      do k2 = 1, ele1%nnod_4_ele
        call SGS_induct_2_each_element                                  &
     &     (ele1%numele, ele1%nnod_4_ele, ele1%ie, ele1%istack_ele_smp, &
     &      k2, i_b, i_v, i_flux, vect_e)
        call fem_skv_div_asym_tsr(iele_cd_smp_stack,                    &
     &      intg_point_t_evo, k2, vect_e, sk6)
      end do
!
      call add3_skv_to_ff_v_smp_1st(ff_nl_smp, sk6)
!
      end subroutine int_div_sgs_idct_simi_pg
!
! ----------------------------------------------------------------------
!
      subroutine int_div_sgs_idct_simi_upw(i_flux, i_v, i_b)
!
      use m_element_phys_address
      use m_element_phys_data
!
      use sgs_terms_2_each_ele
      use fem_skv_vect_diff_upw_1st
      use cal_skv_to_ff_smp_1st
!
      integer(kind = kint), intent(in) :: i_flux, i_v, i_b
!
      integer(kind=kint) :: k2
!
!
      call reset_sk6(n_vector)
      do k2 = 1, ele1%nnod_4_ele
        call SGS_induct_2_each_element                                  &
     &     (ele1%numele, ele1%nnod_4_ele, ele1%ie, ele1%istack_ele_smp, &
     &      k2, i_b, i_v, i_flux, vect_e)
        call fem_skv_div_as_tsr_upw(iele_cd_smp_stack,                  &
     &      intg_point_t_evo, k2, d_ele(1,iphys_ele%i_velo),            &
     &      vect_e, sk6)
      end do
!
      call add3_skv_to_ff_v_smp_1st(ff_nl_smp, sk6)
!
      end subroutine int_div_sgs_idct_simi_upw
!
! ----------------------------------------------------------------------
!
      end module int_div_sgs_idct_simi
