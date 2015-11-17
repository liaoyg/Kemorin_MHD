!
!      module set_ele_4_nodal_bc
!
!      Written by H. Matsui nad H. Okuda
!      Modified by H. Matsui on Oct., 2005
!
!!      subroutine set_ele_4_vector_nodal_bc                            &
!!     &         (node, ele, ibc, ibc2, nmax_idx_ibc,                   &
!!     &          num_idx_ibc, ele_bc_id, nod_bc_id, nmax_idx_ibc2,     &
!!     &          num_idx_ibc2, ele_bc2_id, nod_bc2_id,                 &
!!     &          ibc_end, ibc_shape, ibc_stack, ibc_stack_smp)
!!      subroutine set_ele_4_vector_nodal_bc_fl                         &
!!     &         (node, ele, ibc, ibc2, nmax_idx_ibc,                   &
!!     &          num_idx_ibc, ele_bc_id, nod_bc_id, nmax_idx_ibc2,     &
!!     &          num_idx_ibc2, ele_bc2_id, nod_bc2_id,                 &
!!     &          ibc_end, ibc_shape, ibc_stack, ibc_stack_smp)
!!
!!      subroutine set_ele_4_scalar_nodal_bc                            &
!!     &         (node, ele, ibc, ibc2, num_idx_ibc,                    &
!!     &          ele_bc_id, nod_bc_id, num_idx_ibc2, ele_bc2_id,       &
!!     &          nod_bc2_id, ibc_end, ibc_shape, ibc_stack,            &
!!     &          ibc_stack_smp)
!!      subroutine set_ele_4_scalar_nodal_bc_fl                         &
!!     &         (node, ele, ibc, ibc2, num_idx_ibc,                    &
!!     &          ele_bc_id, nod_bc_id, num_idx_ibc2, ele_bc2_id,       &
!!     &          nod_bc2_id, ibc_end, ibc_shape, ibc_stack,            &
!!     &          ibc_stack_smp)
!!        type(node_data), intent(in) :: node
!!        type(element_data), intent(in) :: ele
!
      module set_ele_4_nodal_bc
!
      use m_precision
      use m_machine_parameter
      use t_geometry_data
!
      use set_bc_element
      use ordering_ele_4_fix_bd
!
      implicit none
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine set_ele_4_vector_nodal_bc                              &
     &         (node, ele, ibc, ibc2, nmax_idx_ibc,                     &
     &          num_idx_ibc, ele_bc_id, nod_bc_id, nmax_idx_ibc2,       &
     &          num_idx_ibc2, ele_bc2_id, nod_bc2_id,                   &
     &          ibc_end, ibc_shape, ibc_stack, ibc_stack_smp)
!
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
!
      integer (kind=kint), intent(in) :: ibc(node%numnod,3)
      integer (kind=kint), intent(in) :: ibc2(node%numnod,3)
!
      integer (kind=kint), intent(in) :: nmax_idx_ibc
      integer (kind=kint), intent(in) :: num_idx_ibc(3)
      integer (kind=kint), intent(inout) :: ele_bc_id(nmax_idx_ibc,3)
      integer (kind=kint), intent(inout) :: nod_bc_id(nmax_idx_ibc,3)
!
      integer (kind=kint), intent(in) :: nmax_idx_ibc2
      integer (kind=kint), intent(in) :: num_idx_ibc2(3)
      integer (kind=kint), intent(inout) :: ele_bc2_id(nmax_idx_ibc2,3)
      integer (kind=kint), intent(inout) :: nod_bc2_id(nmax_idx_ibc2,3)
!
!
      integer (kind=kint), intent(inout) :: ibc_end(3)
      integer (kind=kint), intent(inout) :: ibc_shape(ele%nnod_4_ele,3)
      integer (kind=kint), intent(inout)                                &
     &        :: ibc_stack(0:ele%nnod_4_ele,3)
      integer (kind=kint), intent(inout)                                &
     &        :: ibc_stack_smp(0:ele%nnod_4_ele*np_smp,3)
!
      integer(kind = kint) :: nd
!
!
      do nd = 1, 3
        call set_bc_element_whole(node, ele,                            &
     &      num_idx_ibc(nd), ibc(1,nd), ele_bc_id(1,nd),                &
     &      nod_bc_id(1,nd),  ele%nnod_4_ele)
        call set_bc_element_whole(node, ele,                            &
     &      num_idx_ibc2(nd), ibc2(1,nd), ele_bc2_id(1,nd),             &
     &      nod_bc2_id(1,nd), ele%nnod_4_ele)
!
        call reordering_ele_4_fix_bd(np_smp, nmax_idx_ibc,              &
     &      num_idx_ibc(nd), ele_bc_id(1,nd), nod_bc_id(1,nd),          &
     &      ibc_end(nd), ibc_shape(1,nd), ibc_stack(0,nd),              &
     &      ibc_stack_smp(0,nd), ele%nnod_4_ele)
      end do
!
      end subroutine set_ele_4_vector_nodal_bc
!
!  ---------------------------------------------------------------------
!
      subroutine set_ele_4_vector_nodal_bc_fl                           &
     &         (node, ele, ibc, ibc2, nmax_idx_ibc,                     &
     &          num_idx_ibc, ele_bc_id, nod_bc_id, nmax_idx_ibc2,       &
     &          num_idx_ibc2, ele_bc2_id, nod_bc2_id,                   &
     &          ibc_end, ibc_shape, ibc_stack, ibc_stack_smp)
!
      use m_geometry_data_MHD
!
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
!
      integer (kind=kint), intent(in) :: ibc(node%numnod,3)
      integer (kind=kint), intent(in) :: ibc2(node%numnod,3)
!
      integer (kind=kint), intent(in) :: nmax_idx_ibc
      integer (kind=kint), intent(in) :: num_idx_ibc(3)
      integer (kind=kint), intent(inout) :: ele_bc_id(nmax_idx_ibc,3)
      integer (kind=kint), intent(inout) :: nod_bc_id(nmax_idx_ibc,3)
!
      integer (kind=kint), intent(in) :: nmax_idx_ibc2
      integer (kind=kint), intent(in) :: num_idx_ibc2(3)
      integer (kind=kint), intent(inout) :: ele_bc2_id(nmax_idx_ibc2,3)
      integer (kind=kint), intent(inout) :: nod_bc2_id(nmax_idx_ibc2,3)
!
!
      integer (kind=kint), intent(inout) :: ibc_end(3)
      integer (kind=kint), intent(inout) :: ibc_shape(ele%nnod_4_ele,3)
      integer (kind=kint), intent(inout)                                &
     &         :: ibc_stack(0:ele%nnod_4_ele,3)
      integer (kind=kint), intent(inout)                                &
     &        :: ibc_stack_smp(0:ele%nnod_4_ele*np_smp,3)
!
      integer(kind = kint) :: nd
!
!
      do nd = 1, 3
        call set_bc_element_layer                                       &
     &     (node, ele, iele_fl_start, iele_fl_end,                      &
     &      num_idx_ibc(nd), ibc(1,nd),                                 &
     &      ele_bc_id(1,nd),  nod_bc_id(1,nd),  ele%nnod_4_ele)
        call set_bc_element_layer                                       &
     &     (node, ele, iele_fl_start, iele_fl_end,                      &
     &      num_idx_ibc2(nd), ibc2(1,nd),                               &
     &      ele_bc2_id(1,nd), nod_bc2_id(1,nd), ele%nnod_4_ele)
!
        call reordering_ele_4_fix_bd(np_smp, nmax_idx_ibc,              &
     &      num_idx_ibc(nd), ele_bc_id(1,nd), nod_bc_id(1,nd),          &
     &      ibc_end(nd), ibc_shape(1,nd), ibc_stack(0,nd),              &
     &      ibc_stack_smp(0,nd), ele%nnod_4_ele)
      end do
!
      end subroutine set_ele_4_vector_nodal_bc_fl
!
!  ---------------------------------------------------------------------
!
      subroutine set_ele_4_scalar_nodal_bc                              &
     &         (node, ele, ibc, ibc2, num_idx_ibc,                      &
     &          ele_bc_id, nod_bc_id, num_idx_ibc2, ele_bc2_id,         &
     &          nod_bc2_id, ibc_end, ibc_shape, ibc_stack,              &
     &          ibc_stack_smp)
!
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
!
      integer (kind=kint), intent(in) :: ibc(node%numnod)
      integer (kind=kint), intent(in) :: ibc2(node%numnod)
!
      integer (kind=kint), intent(in) :: num_idx_ibc
      integer (kind=kint), intent(inout) :: ele_bc_id(num_idx_ibc)
      integer (kind=kint), intent(inout) :: nod_bc_id(num_idx_ibc)
!
      integer (kind=kint), intent(in) :: num_idx_ibc2
      integer (kind=kint), intent(inout) :: ele_bc2_id(num_idx_ibc2)
      integer (kind=kint), intent(inout) :: nod_bc2_id(num_idx_ibc2)
!
!
      integer (kind=kint), intent(inout) :: ibc_end
      integer (kind=kint), intent(inout) :: ibc_shape(ele%nnod_4_ele)
      integer (kind=kint), intent(inout) :: ibc_stack(0:ele%nnod_4_ele)
      integer (kind=kint), intent(inout)                                &
     &        :: ibc_stack_smp(0:ele%nnod_4_ele*np_smp)
!
!
        call set_bc_element_whole(node, ele,                            &
     &      num_idx_ibc, ibc, ele_bc_id, nod_bc_id, ele%nnod_4_ele)
        call set_bc_element_whole(node, ele,                            &
     &      num_idx_ibc2, ibc2, ele_bc2_id, nod_bc2_id, ele%nnod_4_ele)
!
        call reordering_ele_4_fix_bd(np_smp, num_idx_ibc,               &
     &      num_idx_ibc, ele_bc_id, nod_bc_id, ibc_end, ibc_shape,      &
     &      ibc_stack, ibc_stack_smp, ele%nnod_4_ele)
!
      end subroutine set_ele_4_scalar_nodal_bc
!
!  ---------------------------------------------------------------------
!
      subroutine set_ele_4_scalar_nodal_bc_fl                           &
     &         (node, ele, ibc, ibc2, num_idx_ibc,                      &
     &          ele_bc_id, nod_bc_id, num_idx_ibc2, ele_bc2_id,         &
     &          nod_bc2_id, ibc_end, ibc_shape, ibc_stack,              &
     &          ibc_stack_smp)
!
      use m_geometry_data_MHD
!
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
!
      integer(kind=kint), intent(in) :: ibc(node%numnod)
      integer(kind=kint), intent(in) :: ibc2(node%numnod)
!
      integer(kind=kint), intent(in) :: num_idx_ibc
      integer(kind=kint), intent(inout) :: ele_bc_id(num_idx_ibc)
      integer(kind=kint), intent(inout) :: nod_bc_id(num_idx_ibc)
!
      integer(kind=kint), intent(in) :: num_idx_ibc2
      integer(kind=kint), intent(inout) :: ele_bc2_id(num_idx_ibc2)
      integer(kind=kint), intent(inout) :: nod_bc2_id(num_idx_ibc2)
!
      integer(kind=kint), intent(inout) :: ibc_end
      integer(kind=kint), intent(inout) :: ibc_shape(ele%nnod_4_ele)
      integer(kind=kint), intent(inout) :: ibc_stack(0:ele%nnod_4_ele)
      integer(kind=kint), intent(inout)                                &
     &        :: ibc_stack_smp(0:ele%nnod_4_ele*np_smp)
!
!
      call set_bc_element_layer                                         &
     &   (node, ele, iele_fl_start, iele_fl_end,                        &
     &    num_idx_ibc, ibc, ele_bc_id, nod_bc_id, ele%nnod_4_ele)
      call set_bc_element_layer                                         &
     &   (node, ele, iele_fl_start, iele_fl_end,                        &
     &    num_idx_ibc2, ibc2,  ele_bc2_id, nod_bc2_id, ele%nnod_4_ele)
!
      call reordering_ele_4_fix_bd(np_smp, num_idx_ibc, num_idx_ibc,    &
     &      ele_bc_id, nod_bc_id, ibc_end, ibc_shape, ibc_stack,        &
     &      ibc_stack_smp, ele%nnod_4_ele)
!
      end subroutine set_ele_4_scalar_nodal_bc_fl
!
!  ---------------------------------------------------------------------
!
      end module set_ele_4_nodal_bc
