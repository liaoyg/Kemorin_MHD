!fem_skv_vect_diff_upw_type.f90
!      module fem_skv_vect_diff_upw_type
!
!     programmed by H.Matsui on May 2012
!
!!      subroutine fem_skv_gradient_upw(iele_fsmp_stack, n_int,         &
!!     &          k2, dt, vxe, ele, g_FEM, jac_3d, scalar_1, sk_v)
!!      subroutine fem_skv_divergence_upw(iele_fsmp_stack, n_int,       &
!!     &          k2, dt, vxe, ele, g_FEM, jac_3d, vector_1, sk_v)
!!      subroutine fem_skv_rotation_upw(iele_fsmp_stack, n_int,         &
!!     &          k2, dt, vxe, ele, g_FEM, jac_3d, vector_1, sk_v)
!!
!!      subroutine fem_skv_div_tsr_upw(iele_fsmp_stack,                 &
!!     &          n_int, k2, dt, vxe, ele, g_FEM, jac_3d, flux_1, sk_v)
!!      subroutine fem_skv_div_as_tsr_upw(iele_fsmp_stack,              &
!!     &          n_int, k2, dt, vxe, ele, g_FEM, jac_3d, flux_1, sk_v)
!!
!!      subroutine fem_skv_grp_gradient_upw                             &
!!     &         (iele_fsmp_stack, nele_grp, iele_grp, n_int, k2, dt,   &
!!     &          vxe, ele, g_FEM, jac_3d, scalar_1, sk_v)
!!      subroutine fem_skv_grp_divergence_upw                           &
!!     &         (iele_fsmp_stack, nele_grp, iele_grp, n_int, k2, dt,   &
!!     &          vxe, ele, g_FEM, jac_3d, vector_1, sk_v)
!!      subroutine fem_skv_grp_rotation_upw                             &
!!     &         (iele_fsmp_stack, nele_grp, iele_grp, n_int, k2,       &
!!     &          dt, vxe, ele, g_FEM, jac_3d, vector_1, sk_v)
!!
!!      subroutine fem_skv_grp_div_tsr_upw                              &
!!     &         (iele_fsmp_stack, nele_grp, iele_grp, n_int, k2,       &
!!     &         dt, vxe, ele, g_FEM, jac_3d, flux_1, sk_v)
!!      subroutine fem_skv_grp_div_as_tsr_upw                           &
!!     &         (iele_fsmp_stack, nele_grp, iele_grp, n_int, k2,       &
!!     &          dt, vxe, ele, g_FEM, jac_3d, flux_1, sk_v)
!!
!!      subroutine fem_skv_linear_gradient_upw                          &
!!     &         (iele_fsmp_stack, n_int, k2, dt, vxe,                  &
!!     &          ele, g_FEM, jac_3d, jac_3d_l, scalar_1, sk_v)
!!        type(element_data), intent(in) :: ele
!!        type(FEM_gauss_int_coefs), intent(in) :: g_FEM
!!        type(jacobians_3d), intent(in) :: jac_3d
!!        type(jacobians_3d), intent(in) :: jac_3d_l
!!        type(work_finite_element_mat), intent(inout) :: fem_wk
!!        integer(kind=kint), intent(in) :: n_int, k2
!!        integer(kind=kint), intent(in) :: iele_fsmp_stack(0:np_smp)
!!        real (kind=kreal), intent(in) :: vxe(ele%numele,3)
!
!
      module fem_skv_vect_diff_upw_type
!
      use m_precision
      use m_machine_parameter
      use m_geometry_constants
!
      use t_geometry_data
      use t_fem_gauss_int_coefs
      use t_jacobians
!
      implicit none
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine fem_skv_gradient_upw(iele_fsmp_stack, n_int,           &
     &          k2, dt, vxe, ele, g_FEM, jac_3d, scalar_1, sk_v)
!
      use fem_skv_grad_upw
!
      type(element_data), intent(in) :: ele
      type(FEM_gauss_int_coefs), intent(in) :: g_FEM
      type(jacobians_3d), intent(in) :: jac_3d
!
      integer(kind=kint), intent(in) :: n_int, k2
      integer(kind=kint), intent(in) :: iele_fsmp_stack(0:np_smp)
      real(kind=kreal), intent(in) :: dt
      real(kind=kreal), intent(in) :: vxe(ele%numele,3)
      real(kind=kreal), intent(in) :: scalar_1(ele%numele)
!
      real (kind=kreal), intent(inout)                                  &
     &            :: sk_v(ele%numele,n_sym_tensor,ele%nnod_4_ele)
!
!
      call fem_skv_all_grad_upw                                         &
     &   (ele%numele, ele%nnod_4_ele, ele%nnod_4_ele,                   &
     &    np_smp, iele_fsmp_stack, g_FEM%max_int_point,                 &
     &    g_FEM%maxtot_int_3d, g_FEM%int_start3, g_FEM%owe3d,           &
     &    n_int, k2, dt, jac_3d%ntot_int, jac_3d%xjac,                  &
     &    jac_3d%an, jac_3d%dnx, jac_3d%dnx, vxe, scalar_1, sk_v)
!
      end subroutine fem_skv_gradient_upw
!
!-----------------------------------------------------------------------
!
      subroutine fem_skv_divergence_upw(iele_fsmp_stack, n_int,         &
     &          k2, dt, vxe, ele, g_FEM, jac_3d, vector_1, sk_v)
!
      use fem_skv_div_upw
!
      type(element_data), intent(in) :: ele
      type(FEM_gauss_int_coefs), intent(in) :: g_FEM
      type(jacobians_3d), intent(in) :: jac_3d
!
      integer(kind=kint), intent(in) :: n_int, k2
      integer(kind=kint), intent(in) :: iele_fsmp_stack(0:np_smp)
      real(kind=kreal), intent(in) :: dt
      real(kind=kreal), intent(in) :: vxe(ele%numele,3)
      real(kind=kreal), intent(in) :: vector_1(ele%numele,3)
!
      real (kind=kreal), intent(inout)                                  &
     &            :: sk_v(ele%numele,n_sym_tensor,ele%nnod_4_ele)
!
!
      call fem_skv_all_div_upw                                          &
     &   (ele%numele, ele%nnod_4_ele, ele%nnod_4_ele,                   &
     &    np_smp, iele_fsmp_stack, g_FEM%max_int_point,                 &
     &    g_FEM%maxtot_int_3d, g_FEM%int_start3, g_FEM%owe3d,           &
     &    n_int, k2, dt, jac_3d%ntot_int, jac_3d%xjac,                  &
     &    jac_3d%an, jac_3d%dnx, jac_3d%dnx, vxe, vector_1, sk_v)
!
      end subroutine fem_skv_divergence_upw
!
!-----------------------------------------------------------------------
!
      subroutine fem_skv_rotation_upw(iele_fsmp_stack, n_int,           &
     &          k2, dt, vxe, ele, g_FEM, jac_3d, vector_1, sk_v)
!
      use fem_skv_rot_upw
!
      type(element_data), intent(in) :: ele
      type(FEM_gauss_int_coefs), intent(in) :: g_FEM
      type(jacobians_3d), intent(in) :: jac_3d
!
      integer(kind=kint), intent(in) :: n_int, k2
      integer(kind=kint), intent(in) :: iele_fsmp_stack(0:np_smp)
      real(kind=kreal), intent(in) :: dt
      real(kind=kreal), intent(in) :: vxe(ele%numele,3)
      real(kind=kreal), intent(in) :: vector_1(ele%numele,3)
!
      real (kind=kreal), intent(inout)                                  &
     &            :: sk_v(ele%numele,n_sym_tensor,ele%nnod_4_ele)
!
!
      call fem_all_skv_rot_upw                                          &
     &   (ele%numele, ele%nnod_4_ele, ele%nnod_4_ele,                   &
     &    np_smp, iele_fsmp_stack, g_FEM%max_int_point,                 &
     &    g_FEM%maxtot_int_3d, g_FEM%int_start3, g_FEM%owe3d,           &
     &    n_int, k2, dt, jac_3d%ntot_int, jac_3d%xjac,                  &
     &    jac_3d%an, jac_3d%dnx, jac_3d%dnx, vxe, vector_1, sk_v)
!
      end subroutine fem_skv_rotation_upw
!
!-----------------------------------------------------------------------
!
      subroutine fem_skv_div_tsr_upw(iele_fsmp_stack,                   &
     &          n_int, k2, dt, vxe, ele, g_FEM, jac_3d, flux_1, sk_v)
!
      use fem_skv_div_flux_upw
!
      type(element_data), intent(in) :: ele
      type(FEM_gauss_int_coefs), intent(in) :: g_FEM
      type(jacobians_3d), intent(in) :: jac_3d
!
      integer(kind=kint), intent(in) :: n_int, k2
      integer(kind=kint), intent(in) :: iele_fsmp_stack(0:np_smp)
      real(kind=kreal), intent(in) :: dt
      real(kind=kreal), intent(in) :: vxe(ele%numele,3)
      real(kind=kreal), intent(in) :: flux_1(ele%numele,n_sym_tensor)
!
      real (kind=kreal), intent(inout)                                  &
     &            :: sk_v(ele%numele,n_sym_tensor,ele%nnod_4_ele)
!
!
      call fem_skv_all_div_flux_upw                                     &
     &   (ele%numele, ele%nnod_4_ele, ele%nnod_4_ele,                   &
     &    np_smp, iele_fsmp_stack, g_FEM%max_int_point,                 &
     &    g_FEM%maxtot_int_3d, g_FEM%int_start3, g_FEM%owe3d,           &
     &    n_int, k2, dt, jac_3d%ntot_int, jac_3d%xjac,                  &
     &    jac_3d%an, jac_3d%dnx, jac_3d%dnx, vxe, flux_1, sk_v)
!
      end subroutine fem_skv_div_tsr_upw
!
!-----------------------------------------------------------------------
!
      subroutine fem_skv_div_as_tsr_upw(iele_fsmp_stack,                &
     &          n_int, k2, dt, vxe, ele, g_FEM, jac_3d, flux_1, sk_v)
!
      use fem_skv_div_asym_t_upw
!
      type(element_data), intent(in) :: ele
      type(FEM_gauss_int_coefs), intent(in) :: g_FEM
      type(jacobians_3d), intent(in) :: jac_3d
!
      integer(kind=kint), intent(in) :: n_int, k2
      integer(kind=kint), intent(in) :: iele_fsmp_stack(0:np_smp)
      real(kind=kreal), intent(in) :: dt
      real(kind=kreal), intent(in) :: vxe(ele%numele,3)
      real(kind=kreal), intent(in) :: flux_1(ele%numele,3)
!
      real (kind=kreal), intent(inout)                                  &
     &            :: sk_v(ele%numele,n_sym_tensor,ele%nnod_4_ele)
!
!
      call fem_skv_all_div_asym_t_upw                                   &
     &   (ele%numele, ele%nnod_4_ele, ele%nnod_4_ele,                   &
     &    np_smp, iele_fsmp_stack, g_FEM%max_int_point,                 &
     &    g_FEM%maxtot_int_3d, g_FEM%int_start3, g_FEM%owe3d,           &
     &    n_int, k2, dt, jac_3d%ntot_int, jac_3d%xjac,                  &
     &    jac_3d%an, jac_3d%dnx, jac_3d%dnx, vxe, flux_1, sk_v)
!
      end subroutine fem_skv_div_as_tsr_upw
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine fem_skv_grp_gradient_upw                               &
     &         (iele_fsmp_stack, nele_grp, iele_grp, n_int, k2, dt,     &
     &          vxe, ele, g_FEM, jac_3d, scalar_1, sk_v)
!
      use fem_skv_grad_upw
!
      type(element_data), intent(in) :: ele
      type(FEM_gauss_int_coefs), intent(in) :: g_FEM
      type(jacobians_3d), intent(in) :: jac_3d
!
      integer(kind=kint), intent(in) :: n_int, k2
      integer(kind=kint), intent(in) :: iele_fsmp_stack(0:np_smp)
      integer(kind=kint), intent(in) :: nele_grp
      integer(kind=kint), intent(in) :: iele_grp(nele_grp)
      real(kind=kreal), intent(in) :: dt
      real(kind=kreal), intent(in) :: vxe(ele%numele,3)
      real(kind=kreal),   intent(in) :: scalar_1(ele%numele)
!
      real (kind=kreal), intent(inout)                                  &
     &            :: sk_v(ele%numele,n_sym_tensor,ele%nnod_4_ele)
!
!
      call fem_skv_grp_grad_upw                                         &
     &   (ele%numele, ele%nnod_4_ele, ele%nnod_4_ele,                   &
     &    np_smp, iele_fsmp_stack, nele_grp, iele_grp,                  &
     &    g_FEM%max_int_point, g_FEM%maxtot_int_3d, g_FEM%int_start3,   &
     &    g_FEM%owe3d, n_int, k2, dt, jac_3d%ntot_int, jac_3d%xjac,     &
     &    jac_3d%an, jac_3d%dnx, jac_3d%dnx,                            &
     &    vxe, scalar_1, sk_v)
!
      end subroutine fem_skv_grp_gradient_upw
!
!-----------------------------------------------------------------------
!
      subroutine fem_skv_grp_divergence_upw                             &
     &         (iele_fsmp_stack, nele_grp, iele_grp, n_int, k2, dt,     &
     &          vxe, ele, g_FEM, jac_3d, vector_1, sk_v)
!
      use fem_skv_div_upw
!
      type(element_data), intent(in) :: ele
      type(FEM_gauss_int_coefs), intent(in) :: g_FEM
      type(jacobians_3d), intent(in) :: jac_3d
!
      integer(kind=kint), intent(in) :: n_int, k2
      integer(kind=kint), intent(in) :: iele_fsmp_stack(0:np_smp)
      integer(kind=kint), intent(in) :: nele_grp
      integer(kind=kint), intent(in) :: iele_grp(nele_grp)
      real(kind=kreal), intent(in) :: dt
      real(kind=kreal), intent(in) :: vxe(ele%numele,3)
      real(kind=kreal), intent(in) :: vector_1(ele%numele,3)
!
      real (kind=kreal), intent(inout)                                  &
     &            :: sk_v(ele%numele,n_sym_tensor,ele%nnod_4_ele)
!
!
      call fem_skv_grp_div_upw                                          &
     &   (ele%numele, ele%nnod_4_ele, ele%nnod_4_ele,                   &
     &    np_smp, iele_fsmp_stack, nele_grp, iele_grp,                  &
     &    g_FEM%max_int_point, g_FEM%maxtot_int_3d, g_FEM%int_start3,   &
     &    g_FEM%owe3d, n_int, k2, dt, jac_3d%ntot_int, jac_3d%xjac,     &
     &    jac_3d%an, jac_3d%dnx, jac_3d%dnx, vxe, vector_1, sk_v)
!
      end subroutine fem_skv_grp_divergence_upw
!
!-----------------------------------------------------------------------
!
      subroutine fem_skv_grp_rotation_upw                               &
     &         (iele_fsmp_stack, nele_grp, iele_grp, n_int, k2,         &
     &          dt, vxe, ele, g_FEM, jac_3d, vector_1, sk_v)
!
      use fem_skv_rot_upw
!
      type(element_data), intent(in) :: ele
      type(FEM_gauss_int_coefs), intent(in) :: g_FEM
      type(jacobians_3d), intent(in) :: jac_3d
!
      integer(kind=kint), intent(in) :: n_int, k2
      integer(kind=kint), intent(in) :: iele_fsmp_stack(0:np_smp)
      integer(kind=kint), intent(in) :: nele_grp
      integer(kind=kint), intent(in) :: iele_grp(nele_grp)
      real(kind=kreal), intent(in) :: dt
      real(kind=kreal), intent(in) :: vxe(ele%numele,3)
      real(kind=kreal), intent(in) :: vector_1(ele%numele,3)
!
      real (kind=kreal), intent(inout)                                  &
     &            :: sk_v(ele%numele,n_sym_tensor,ele%nnod_4_ele)
!
!
      call fem_skv_grp_rot_upw                                          &
     &   (ele%numele, ele%nnod_4_ele, ele%nnod_4_ele,                   &
     &    np_smp, iele_fsmp_stack, nele_grp, iele_grp,                  &
     &    g_FEM%max_int_point, g_FEM%maxtot_int_3d, g_FEM%int_start3,   &
     &    g_FEM%owe3d, n_int, k2, dt, jac_3d%ntot_int, jac_3d%xjac,     &
     &    jac_3d%an, jac_3d%dnx, jac_3d%dnx, vxe, vector_1, sk_v)
!
      end subroutine fem_skv_grp_rotation_upw
!
!-----------------------------------------------------------------------
!
      subroutine fem_skv_grp_div_tsr_upw                                &
     &         (iele_fsmp_stack, nele_grp, iele_grp, n_int, k2,         &
     &         dt, vxe, ele, g_FEM, jac_3d, flux_1, sk_v)
!
      use fem_skv_div_flux_upw
!
      type(element_data), intent(in) :: ele
      type(FEM_gauss_int_coefs), intent(in) :: g_FEM
      type(jacobians_3d), intent(in) :: jac_3d
!
      integer(kind=kint), intent(in) :: n_int, k2
      integer(kind=kint), intent(in) :: iele_fsmp_stack(0:np_smp)
      integer(kind=kint), intent(in) :: nele_grp
      integer(kind=kint), intent(in) :: iele_grp(nele_grp)
      real(kind=kreal), intent(in) :: dt
      real(kind=kreal), intent(in) :: vxe(ele%numele,3)
      real(kind=kreal), intent(in) :: flux_1(ele%numele,n_sym_tensor)
!
      real (kind=kreal), intent(inout)                                  &
     &            :: sk_v(ele%numele,n_sym_tensor,ele%nnod_4_ele)
!
!
      call fem_skv_grp_div_flux_upw                                     &
     &   (ele%numele, ele%nnod_4_ele, ele%nnod_4_ele,                   &
     &    np_smp, iele_fsmp_stack, nele_grp, iele_grp,                  &
     &    g_FEM%max_int_point, g_FEM%maxtot_int_3d, g_FEM%int_start3,   &
     &    g_FEM%owe3d, n_int, k2, dt, jac_3d%ntot_int, jac_3d%xjac,     &
     &    jac_3d%an, jac_3d%dnx, jac_3d%dnx, vxe, flux_1, sk_v)
!
      end subroutine fem_skv_grp_div_tsr_upw
!
!-----------------------------------------------------------------------
!
      subroutine fem_skv_grp_div_as_tsr_upw                             &
     &         (iele_fsmp_stack, nele_grp, iele_grp, n_int, k2,         &
     &          dt, vxe, ele, g_FEM, jac_3d, flux_1, sk_v)
!
      use fem_skv_div_asym_t_upw
!
      type(element_data), intent(in) :: ele
      type(FEM_gauss_int_coefs), intent(in) :: g_FEM
      type(jacobians_3d), intent(in) :: jac_3d
!
      integer(kind=kint), intent(in) :: n_int, k2
      integer(kind=kint), intent(in) :: iele_fsmp_stack(0:np_smp)
      integer(kind=kint), intent(in) :: nele_grp
      integer(kind=kint), intent(in) :: iele_grp(nele_grp)
      real(kind=kreal), intent(in) :: dt
      real(kind=kreal), intent(in) :: vxe(ele%numele,3)
      real(kind=kreal), intent(in) :: flux_1(ele%numele,3)
!
      real (kind=kreal), intent(inout)                                  &
     &            :: sk_v(ele%numele,n_sym_tensor,ele%nnod_4_ele)
!
!
      call fem_skv_grp_div_asym_t_upw                                   &
     &   (ele%numele, ele%nnod_4_ele, ele%nnod_4_ele,                   &
     &    np_smp, iele_fsmp_stack, nele_grp, iele_grp,                  &
     &    g_FEM%max_int_point, g_FEM%maxtot_int_3d, g_FEM%int_start3,   &
     &    g_FEM%owe3d, n_int, k2, dt, jac_3d%ntot_int, jac_3d%xjac,     &
     &    jac_3d%an, jac_3d%dnx, jac_3d%dnx, vxe, flux_1, sk_v)
!
      end subroutine fem_skv_grp_div_as_tsr_upw
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine fem_skv_linear_gradient_upw                            &
     &         (iele_fsmp_stack, n_int, k2, dt, vxe,                    &
     &          ele, g_FEM, jac_3d, jac_3d_l, scalar_1, sk_v)
!
      use fem_skv_grad_upw
!
      type(element_data), intent(in) :: ele
      type(FEM_gauss_int_coefs), intent(in) :: g_FEM
      type(jacobians_3d), intent(in) :: jac_3d
      type(jacobians_3d), intent(in) :: jac_3d_l
!
      integer(kind=kint), intent(in) :: n_int, k2
      integer(kind=kint), intent(in) :: iele_fsmp_stack(0:np_smp)
      real(kind=kreal), intent(in) :: dt
      real(kind=kreal), intent(in) :: vxe(ele%numele,3)
      real(kind=kreal), intent(in) :: scalar_1(ele%numele)
!
      real (kind=kreal), intent(inout)                                  &
     &           :: sk_v(ele%numele,n_sym_tensor,ele%nnod_4_ele)
!
!
      call fem_skv_all_grad_upw                                         &
     &   (ele%numele, ele%nnod_4_ele, num_t_linear,                     &
     &    np_smp, iele_fsmp_stack, g_FEM%max_int_point,                 &
     &    g_FEM%maxtot_int_3d, g_FEM%int_start3, g_FEM%owe3d,           &
     &    n_int, k2, dt, jac_3d%ntot_int, jac_3d%xjac,                  &
     &    jac_3d%an, jac_3d%dnx, jac_3d_l%dnx, vxe, scalar_1, sk_v)
!
      end subroutine fem_skv_linear_gradient_upw
!
!-----------------------------------------------------------------------
!
      end module fem_skv_vect_diff_upw_type
