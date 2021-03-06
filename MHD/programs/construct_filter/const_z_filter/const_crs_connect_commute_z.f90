!
!      module const_crs_connect_commute_z
!
!     Written by Hiroaki Matsui
!
!      subroutine set_connect_2_n_filter(node)
!      subroutine set_spatial_difference(numele, n_int, g_FEM, jac_1d)
!      subroutine set_crs_connect_commute_z(node, tbl_crs)
!
      module const_crs_connect_commute_z
!
      use m_precision
!
      use m_constants
      use m_commute_filter_z
!
      use t_crs_connect
!
      implicit none
!
      private :: set_num_off_diag_z_commute
      private :: set_stack_crs_z_commute, set_item_crs_z_commute
!
! ----------------------------------------------------------------------
!
      contains
!
!   --------------------------------------------------------------------
!
      subroutine set_connect_2_n_filter(node)
!
      use t_geometry_data
      use m_commute_filter_z
      use m_neibor_data_z
!
      type(node_data), intent(in) :: node
      integer(kind = kint) :: inod
      integer(kind = kint_gl) :: i
!
!
      do inod = 1, node%numnod
        i = node%inod_global(inod)
        ncomp_st(inod) = max(1, 1+nneib_nod(i,1) - (ncomp_mat-1)/2 )
        ncomp_st(inod) = min(ncomp_st(inod)+ncomp_mat-1, nfilter2_3)    &
                  - ncomp_mat + 1
      end do
!
      end subroutine set_connect_2_n_filter
!
!   --------------------------------------------------------------------
!   --------------------------------------------------------------------
!
      subroutine set_spatial_difference(numele, n_int, g_FEM, jac_1d)
!
      use m_commute_filter_z
      use m_int_edge_data
!
      use t_fem_gauss_int_coefs
      use t_jacobian_1d
!
      integer (kind = kint), intent(in) :: numele, n_int
      type(FEM_gauss_int_coefs), intent(in) :: g_FEM
      type(jacobians_1d), intent(in) :: jac_1d
      integer (kind = kint) :: iele, k, ix
!
!
       do k = 1, n_int
         ix = k + g_FEM%int_start1(n_int)
         do iele = 1, numele
           dz(iele) = dz(iele)                                          &
     &               + jac_1d%xeg_edge(iele,ix,3) * g_FEM%owe(ix)
         end do
       end do
       do iele = 1, numele
         dz(iele) = half * dz(iele)
       end do
!
      end subroutine set_spatial_difference
!
!   --------------------------------------------------------------------
!   --------------------------------------------------------------------
!
      subroutine set_crs_connect_commute_z(node, tbl_crs)
!
      use t_geometry_data
!
      type(node_data), intent(inout) :: node
      type(CRS_matrix_connect), intent(inout) :: tbl_crs
!
!
      call set_num_off_diag_z_commute(node%internal_node, tbl_crs)
!
      call alloc_crs_stack(node%numnod, tbl_crs)
      call alloc_crs_connect(tbl_crs)
!
      call set_stack_crs_z_commute(node%internal_node, tbl_crs)
!
      call set_item_crs_z_commute(node%internal_node, tbl_crs)
!
      end subroutine set_crs_connect_commute_z
!
! ----------------------------------------------------------------------
!
      subroutine set_num_off_diag_z_commute(internal_node, tbl_crs)
!
      integer(kind = kint), intent(in) :: internal_node
      type(CRS_matrix_connect), intent(inout) :: tbl_crs
!
!
      tbl_crs%ntot_l = internal_node - 1
      tbl_crs%ntot_u = internal_node - 1
!
      end subroutine set_num_off_diag_z_commute
!
! ----------------------------------------------------------------------
!
      subroutine set_stack_crs_z_commute(internal_node, tbl_crs)
!
      integer(kind = kint), intent(in) :: internal_node
      type(CRS_matrix_connect), intent(inout) :: tbl_crs
      integer(kind = kint) :: i
!
!    set lower component
!
      do i = 1, internal_node
        tbl_crs%istack_l(i) = i-1
      end do
!
!   set upper component
!
      do i = 1, internal_node-1
        tbl_crs%istack_u(i) = i
      end do
      tbl_crs%istack_u(internal_node) = internal_node-1
!
      end subroutine set_stack_crs_z_commute
!
! ----------------------------------------------------------------------
!
      subroutine set_item_crs_z_commute(internal_node, tbl_crs)
!
      integer(kind = kint), intent(in) :: internal_node
      type(CRS_matrix_connect), intent(inout) :: tbl_crs
      integer(kind = kint) :: i
!
!    set lower component
!
      do i = 1, internal_node-1
        tbl_crs%item_l(i) = i
      end do
!
!   set upper component
!
      do i = 1, internal_node-1
        tbl_crs%item_u(i) = i+1
      end do
!
      end subroutine set_item_crs_z_commute
!
! ----------------------------------------------------------------------
!
      end module const_crs_connect_commute_z
