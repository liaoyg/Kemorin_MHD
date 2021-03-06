!int_vol_vect_diff_upw.f90
!      module int_vol_vect_diff_upw
!
!     programmed by H.Matsui on July 2005
!     Modified by H. Matsui on Oct., 2006
!
!!      subroutine int_vol_gradient_upw                                 &
!!     &         (node, ele, g_FEM, jac_3d, rhs_tbl, nod_fld,           &
!!     &          iele_fsmp_stack, num_int, dt, i_field,                &
!!     &          ncomp_ele, iv_up, d_ele, fem_wk, f_nl)
!!      subroutine int_vol_divergence_upw                               &
!!     &         (node, ele, g_FEM, jac_3d, rhs_tbl, nod_fld,           &
!!     &          iele_fsmp_stack, num_int, dt, i_field,                &
!!     &          ncomp_ele, iv_up, d_ele, fem_wk, f_nl)
!!      subroutine int_vol_rotation_upw                                 &
!!     &         (node, ele, g_FEM, jac_3d, rhs_tbl, nod_fld,           &
!!     &          iele_fsmp_stack, num_int, dt, i_field,                &
!!     &          ncomp_ele, iv_up, d_ele, fem_wk, f_nl)
!!
!!      subroutine int_vol_div_tsr_upw                                  &
!!     &         (node, ele, g_FEM, jac_3d, rhs_tbl, nod_fld,           &
!!     &          iele_fsmp_stack, num_int, dt, i_field,                &
!!     &          ncomp_ele, iv_up, d_ele, fem_wk, f_nl)
!!      subroutine int_vol_div_as_tsr_upw                               &
!!     &         (node, ele, g_FEM, jac_3d, rhs_tbl,  nod_fld,          &
!!     &          iele_fsmp_stack, num_int, dt, i_field,                &
!!     &          ncomp_ele, iv_up, d_ele, fem_wk, f_nl)
!!        type(node_data), intent(in) :: node
!!        type(element_data), intent(in) :: ele
!!        type(FEM_gauss_int_coefs), intent(in) :: g_FEM
!!        type(jacobians_3d), intent(in) :: jac_3d
!!        type(phys_data),    intent(in) :: nod_fld
!!        type(tables_4_FEM_assembles), intent(in) :: rhs_tbl
!!        type(work_finite_element_mat), intent(inout) :: fem_wk
!!        type(finite_ele_mat_node), intent(inout) :: f_nl
!
      module int_vol_vect_diff_upw
!
      use m_precision
      use m_phys_constants
!
      use t_geometry_data
      use t_phys_data
      use t_fem_gauss_int_coefs
      use t_jacobians
      use t_table_FEM_const
      use t_finite_element_mat
!
      use nodal_fld_2_each_element
      use cal_skv_to_ff_smp
      use fem_skv_vect_diff_upw_type
!
      implicit none
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine int_vol_gradient_upw                                   &
     &         (node, ele, g_FEM, jac_3d, rhs_tbl, nod_fld,             &
     &          iele_fsmp_stack, num_int, dt, i_field,                  &
     &          ncomp_ele, iv_up, d_ele, fem_wk, f_nl)
!
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
      type(FEM_gauss_int_coefs), intent(in) :: g_FEM
      type(jacobians_3d), intent(in) :: jac_3d
      type(phys_data),    intent(in) :: nod_fld
      type(tables_4_FEM_assembles), intent(in) :: rhs_tbl
!
      integer(kind=kint), intent(in) :: num_int
      integer(kind=kint), intent(in) :: i_field
      integer(kind=kint), intent(in) :: iele_fsmp_stack(0:np_smp)
!
      integer(kind = kint), intent(in) :: ncomp_ele, iv_up
      real(kind = kreal), intent(in) :: d_ele(ele%numele,ncomp_ele)
      real(kind = kreal), intent(in) :: dt
!
      type(work_finite_element_mat), intent(inout) :: fem_wk
      type(finite_ele_mat_node), intent(inout) :: f_nl
!
      integer(kind=kint) :: k2
!
!
      call reset_sk6(n_vector, ele, fem_wk%sk6)
!
! -------- loop for shape function for the field values
!
      do k2 = 1, ele%nnod_4_ele
        call scalar_phys_2_each_element(node, ele, nod_fld,             &
     &      k2, i_field, fem_wk%scalar_1)
        call fem_skv_gradient_upw                                       &
     &     (iele_fsmp_stack, num_int, k2, dt, d_ele(1,iv_up),           &
     &      ele, g_FEM, jac_3d, fem_wk%scalar_1, fem_wk%sk6)
      end do
!
      call add3_skv_to_ff_v_smp                                         &
     &   (node, ele, rhs_tbl, fem_wk%sk6, f_nl%ff_smp)
!
      end subroutine int_vol_gradient_upw
!
!-----------------------------------------------------------------------
!
      subroutine int_vol_divergence_upw                                 &
     &         (node, ele, g_FEM, jac_3d, rhs_tbl, nod_fld,             &
     &          iele_fsmp_stack, num_int, dt, i_field,                  &
     &          ncomp_ele, iv_up, d_ele, fem_wk, f_nl)
!
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
      type(FEM_gauss_int_coefs), intent(in) :: g_FEM
      type(jacobians_3d), intent(in) :: jac_3d
      type(phys_data),    intent(in) :: nod_fld
      type(tables_4_FEM_assembles), intent(in) :: rhs_tbl
!
      integer(kind=kint), intent(in) :: num_int
      integer(kind=kint), intent(in) :: i_field
      integer(kind=kint), intent(in) :: iele_fsmp_stack(0:np_smp)
!
      integer(kind = kint), intent(in) :: ncomp_ele, iv_up
      real(kind = kreal), intent(in) :: d_ele(ele%numele,ncomp_ele)
      real(kind = kreal), intent(in) :: dt
!
      type(work_finite_element_mat), intent(inout) :: fem_wk
      type(finite_ele_mat_node), intent(inout) :: f_nl
!
      integer(kind=kint) :: k2
!
!
      call reset_sk6(n_scalar, ele, fem_wk%sk6)
!
! -------- loop for shape function for the field values
!
      do k2 = 1, ele%nnod_4_ele
        call vector_phys_2_each_element(node, ele, nod_fld,             &
     &      k2, i_field, fem_wk%vector_1)
        call fem_skv_divergence_upw                                     &
     &     (iele_fsmp_stack, num_int, k2, dt, d_ele(1,iv_up),           &
     &      ele, g_FEM, jac_3d, fem_wk%vector_1, fem_wk%sk6)
      end do
!
      call add1_skv_to_ff_v_smp                                         &
     &   (node, ele, rhs_tbl, fem_wk%sk6, f_nl%ff_smp)
!
      end subroutine int_vol_divergence_upw
!
!-----------------------------------------------------------------------
!
      subroutine int_vol_rotation_upw                                   &
     &         (node, ele, g_FEM, jac_3d, rhs_tbl, nod_fld,             &
     &          iele_fsmp_stack, num_int, dt, i_field,                  &
     &          ncomp_ele, iv_up, d_ele, fem_wk, f_nl)
!
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
      type(FEM_gauss_int_coefs), intent(in) :: g_FEM
      type(jacobians_3d), intent(in) :: jac_3d
      type(phys_data),    intent(in) :: nod_fld
      type(tables_4_FEM_assembles), intent(in) :: rhs_tbl
!
      integer(kind=kint), intent(in) :: num_int
      integer(kind=kint), intent(in) :: i_field
      integer(kind=kint), intent(in) :: iele_fsmp_stack(0:np_smp)
!
      integer(kind = kint), intent(in) :: ncomp_ele, iv_up
      real(kind = kreal), intent(in) :: d_ele(ele%numele,ncomp_ele)
      real(kind = kreal), intent(in) :: dt
!
      type(work_finite_element_mat), intent(inout) :: fem_wk
      type(finite_ele_mat_node), intent(inout) :: f_nl
!
      integer(kind=kint) :: k2
!
!
      call reset_sk6(n_vector, ele, fem_wk%sk6)
!
! -------- loop for shape function for the field values
!
      do k2 = 1, ele%nnod_4_ele
        call vector_phys_2_each_element(node, ele, nod_fld,             &
     &      k2, i_field, fem_wk%vector_1)
        call fem_skv_rotation_upw                                       &
     &     (iele_fsmp_stack, num_int, k2, dt, d_ele(1,iv_up),           &
     &      ele, g_FEM, jac_3d, fem_wk%vector_1, fem_wk%sk6)
      end do
!
      call add3_skv_to_ff_v_smp                                         &
     &   (node, ele, rhs_tbl, fem_wk%sk6, f_nl%ff_smp)
!
      end subroutine int_vol_rotation_upw
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine int_vol_div_tsr_upw                                    &
     &         (node, ele, g_FEM, jac_3d, rhs_tbl, nod_fld,             &
     &          iele_fsmp_stack, num_int, dt, i_field,                  &
     &          ncomp_ele, iv_up, d_ele, fem_wk, f_nl)
!
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
      type(FEM_gauss_int_coefs), intent(in) :: g_FEM
      type(jacobians_3d), intent(in) :: jac_3d
      type(phys_data),    intent(in) :: nod_fld
      type(tables_4_FEM_assembles), intent(in) :: rhs_tbl
!
      integer(kind=kint), intent(in) :: num_int
      integer(kind=kint), intent(in) :: i_field
      integer(kind=kint), intent(in) :: iele_fsmp_stack(0:np_smp)
!
      integer(kind = kint), intent(in) :: ncomp_ele, iv_up
      real(kind = kreal), intent(in) :: d_ele(ele%numele,ncomp_ele)
      real(kind = kreal), intent(in) :: dt
!
      type(work_finite_element_mat), intent(inout) :: fem_wk
      type(finite_ele_mat_node), intent(inout) :: f_nl
!
      integer(kind=kint) :: k2
!
!
      call reset_sk6(n_vector, ele, fem_wk%sk6)
!
! -------- loop for shape function for the field values
!
      do k2 = 1, ele%nnod_4_ele
        call tensor_phys_2_each_element(node, ele, nod_fld,             &
     &      k2, i_field, fem_wk%tensor_1)
        call fem_skv_div_tsr_upw                                        &
     &     (iele_fsmp_stack, num_int, k2, dt, d_ele(1,iv_up),           &
     &      ele, g_FEM, jac_3d, fem_wk%tensor_1, fem_wk%sk6)
      end do
!
      call add3_skv_to_ff_v_smp                                         &
     &   (node, ele, rhs_tbl, fem_wk%sk6, f_nl%ff_smp)
!
      end subroutine int_vol_div_tsr_upw
!
!-----------------------------------------------------------------------
!
      subroutine int_vol_div_as_tsr_upw                                 &
     &         (node, ele, g_FEM, jac_3d, rhs_tbl,  nod_fld,            &
     &          iele_fsmp_stack, num_int, dt, i_field,                  &
     &          ncomp_ele, iv_up, d_ele, fem_wk, f_nl)
!
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
      type(FEM_gauss_int_coefs), intent(in) :: g_FEM
      type(jacobians_3d), intent(in) :: jac_3d
      type(phys_data),    intent(in) :: nod_fld
      type(tables_4_FEM_assembles), intent(in) :: rhs_tbl
!
      integer(kind=kint), intent(in) :: num_int
      integer(kind=kint), intent(in) :: i_field
      integer(kind=kint), intent(in) :: iele_fsmp_stack(0:np_smp)
!
      integer(kind = kint), intent(in) :: ncomp_ele, iv_up
      real(kind = kreal), intent(in) :: d_ele(ele%numele,ncomp_ele)
      real(kind = kreal), intent(in) :: dt
!
      type(work_finite_element_mat), intent(inout) :: fem_wk
      type(finite_ele_mat_node), intent(inout) :: f_nl
!
      integer(kind=kint) :: k2
!
!
      call reset_sk6(n_vector, ele, fem_wk%sk6)
!
! -------- loop for shape function for the field values
!
      do k2 = 1, ele%nnod_4_ele
        call vector_phys_2_each_element(node, ele, nod_fld,             &
     &      k2, i_field, fem_wk%vector_1)
        call fem_skv_div_as_tsr_upw                                     &
     &     (iele_fsmp_stack, num_int, k2, dt, d_ele(1,iv_up),           &
     &      ele, g_FEM, jac_3d, fem_wk%vector_1, fem_wk%sk6)
      end do
!
      call add3_skv_to_ff_v_smp                                         &
     &   (node, ele, rhs_tbl, fem_wk%sk6, f_nl%ff_smp)
!
      end subroutine int_vol_div_as_tsr_upw
!
!-----------------------------------------------------------------------
!
      end module int_vol_vect_diff_upw
