!fem_skv_diffusion_sgs_1st.f90
!      module fem_skv_diffusion_sgs_1st
!
!     programmed by H.Matsui on July 2005
!     Modified by H. Matsui on Oct., 2006
!
!      subroutine fem_skv_scalar_diffuse_1st(iele_fsmp_stack,           &
!     &          n_int, k2, ak_d, scalar_1, sk_v)
!      subroutine fem_skv_vector_diffuse_sgs_1st(iele_fsmp_stack,       &
!     &          n_int, k2, i_filter, ak_diff, ak_d, vect_1, sk_v)
!
!      subroutine fem_skv_poisson_sgs_1st(iele_fsmp_stack,              &
!     &          n_int, k2, i_filter, ak_diff, sk_v)
!      subroutine fem_skv_poisson_linear_sgs_1st(iele_fsmp_stack,       &
!     &          n_int, k2, i_filter, ak_diff, sk_v)
!
!
      module fem_skv_diffusion_sgs_1st
!
      use m_precision
!
      use m_machine_parameter
      use m_geometry_constants
      use m_geometry_data
      use m_fem_gauss_int_coefs
      use m_jacobians
      use m_filter_elength
!
      implicit none
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine fem_skv_scalar_diffuse_sgs_1st(iele_fsmp_stack,        &
     &          n_int, k2, i_filter, ak_diff, ak_d, scalar_1, sk_v)
!
      use fem_skv_diffusion_sgs
!
      integer(kind=kint), intent(in) :: n_int, k2, i_filter
      integer(kind=kint), intent(in) :: iele_fsmp_stack(0:np_smp)
!
      real (kind=kreal), intent(in) :: ak_d(ele1%numele)
      real (kind=kreal), intent(in) :: ak_diff(ele1%numele)
      real(kind=kreal),   intent(in) :: scalar_1(ele1%numele)
!
      real (kind=kreal), intent(inout)                                  &
     &            :: sk_v(ele1%numele,n_sym_tensor,ele1%nnod_4_ele)
!
!
      call fem_skv_scalar_diffuse_sgs                                   &
     & (ele1%numele, ele1%nnod_4_ele, ele1%nnod_4_ele,                  &
     &  np_smp, iele_fsmp_stack, n_int, k2,                             &
     &  jac1_3d_q%ntot_int, jac1_3d_q%xjac, dwx, dwx,                   &
     &  FEM1_elen%filter_conf%xmom_1d_org(i_filter,2),                  &
     &  FEM1_elen%nele_filter_mom,                                      &
     &  FEM1_elen%elen_ele%diff2%df_x2, FEM1_elen%elen_ele%diff2%df_y2, &
     &  FEM1_elen%elen_ele%diff2%df_z2, FEM1_elen%elen_ele%diff2%df_xy, &
     &  FEM1_elen%elen_ele%diff2%df_yz, FEM1_elen%elen_ele%diff2%df_zx, &
     &  ak_diff, ak_d, scalar_1, sk_v)
!
      end subroutine fem_skv_scalar_diffuse_sgs_1st
!
!-----------------------------------------------------------------------
!
      subroutine fem_skv_vector_diffuse_sgs_1st(iele_fsmp_stack,        &
     &          n_int, k2, i_filter, ak_diff, ak_d, vect_1, sk_v)
!
      use fem_skv_diffusion_sgs
!
      integer(kind=kint), intent(in) :: n_int, k2, i_filter
      integer(kind=kint), intent(in) :: iele_fsmp_stack(0:np_smp)
!
      real (kind=kreal), intent(in) :: ak_d(ele1%numele)
      real (kind=kreal), intent(in) :: ak_diff(ele1%numele)
      real(kind=kreal), intent(in) :: vect_1(ele1%numele,3)
!
      real (kind=kreal), intent(inout)                                  &
     &            :: sk_v(ele1%numele,n_sym_tensor,ele1%nnod_4_ele)
!
!
      call fem_skv_vector_diffuse_sgs                                   &
     & (ele1%numele, ele1%nnod_4_ele, ele1%nnod_4_ele,                  &
     &  np_smp, iele_fsmp_stack, n_int, k2,                             &
     &  jac1_3d_q%ntot_int, jac1_3d_q%xjac, dwx, dwx,                   &
     &  FEM1_elen%filter_conf%xmom_1d_org(i_filter,2),                  &
     &  FEM1_elen%nele_filter_mom,                                      &
     &  FEM1_elen%elen_ele%diff2%df_x2, FEM1_elen%elen_ele%diff2%df_y2, &
     &  FEM1_elen%elen_ele%diff2%df_z2, FEM1_elen%elen_ele%diff2%df_xy, &
     &  FEM1_elen%elen_ele%diff2%df_yz, FEM1_elen%elen_ele%diff2%df_zx, &
     &  ak_diff, ak_d, vect_1, sk_v)
!
      end subroutine fem_skv_vector_diffuse_sgs_1st
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine fem_skv_poisson_sgs_1st(iele_fsmp_stack,               &
     &          n_int, k2, i_filter, ak_diff, sk_v)
!
      use fem_skv_poisson_sgs
!
      integer(kind=kint), intent(in) :: n_int, k2, i_filter
      integer(kind=kint), intent(in) :: iele_fsmp_stack(0:np_smp)
!
      real (kind=kreal), intent(in) :: ak_diff(ele1%numele)
!
      real (kind=kreal), intent(inout)                                  &
     &            :: sk_v(ele1%numele,n_sym_tensor,ele1%nnod_4_ele)
!
!
      call fem_skv_poisson_sgs_pg                                       &
     & (ele1%numele, ele1%nnod_4_ele, ele1%nnod_4_ele,                  &
     &  np_smp, iele_fsmp_stack, n_int, k2,                             &
     &  jac1_3d_q%ntot_int, jac1_3d_q%xjac, dwx, dwx,                   &
     &  FEM1_elen%filter_conf%xmom_1d_org(i_filter,2),                  &
     &  FEM1_elen%nele_filter_mom,                                      &
     &  FEM1_elen%elen_ele%diff2%df_x2, FEM1_elen%elen_ele%diff2%df_y2, &
     &  FEM1_elen%elen_ele%diff2%df_z2, FEM1_elen%elen_ele%diff2%df_xy, &
     &  FEM1_elen%elen_ele%diff2%df_yz, FEM1_elen%elen_ele%diff2%df_zx, &
     &  ak_diff, sk_v)
!
      end subroutine fem_skv_poisson_sgs_1st
!
!-----------------------------------------------------------------------
!
      subroutine fem_skv_poisson_linear_sgs_1st(iele_fsmp_stack,        &
     &          n_int, k2, i_filter, ak_diff, sk_v)
!
      use fem_skv_poisson_sgs
!
      integer(kind=kint), intent(in) :: n_int, k2, i_filter
      integer(kind=kint), intent(in) :: iele_fsmp_stack(0:np_smp)
!
      real (kind=kreal), intent(in) :: ak_diff(ele1%numele)
!
      real (kind=kreal), intent(inout)                                  &
     &            :: sk_v(ele1%numele,n_sym_tensor,ele1%nnod_4_ele)
!
!
      call fem_skv_poisson_sgs_pg                                       &
     & (ele1%numele, num_t_linear, num_t_linear,                        &
     &  np_smp, iele_fsmp_stack, n_int, k2, jac1_3d_l%ntot_int,         &
     &  jac1_3d_l%xjac, jac1_3d_l%dnx, jac1_3d_l%dnx,                   &
     &  FEM1_elen%filter_conf%xmom_1d_org(i_filter,2),                  &
     &  FEM1_elen%nele_filter_mom,                                      &
     &  FEM1_elen%elen_ele%diff2%df_x2, FEM1_elen%elen_ele%diff2%df_y2, &
     &  FEM1_elen%elen_ele%diff2%df_z2, FEM1_elen%elen_ele%diff2%df_xy, &
     &  FEM1_elen%elen_ele%diff2%df_yz, FEM1_elen%elen_ele%diff2%df_zx, &
     &  ak_diff, sk_v)
!
      end subroutine fem_skv_poisson_linear_sgs_1st
!
!-----------------------------------------------------------------------
!
      end module fem_skv_diffusion_sgs_1st
