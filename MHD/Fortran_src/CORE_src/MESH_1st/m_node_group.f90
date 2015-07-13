!m_node_group.f90
!
!     module m_node_group
!
!> @brief node group data
!
!     Written by H. Matsui
!
!      subroutine allocate_boundary_data
!      subroutine clear_boundary_data
!      subroutine allocate_boundary_param_smp
!
!      subroutine deallocate_boundary_data
!      subroutine deallocate_boundary_param_smp
!
!      subroutine check_bc_4_sheard_para(my_rank)
!
      module m_node_group
!
      use m_precision
      use t_group_data
!
      implicit  none
!
!>  Structure for node and node group
      type(group_data), save :: nod_grp1
!
!nod_grp1%grp_name
!
!      integer (kind=kint) :: num_bc
!<      number of node group
!      integer (kind=kint) :: num_nod_bc
!<      total number of nodes for node group
!
!      integer (kind=kint), allocatable, target :: bc_istack(:)
!<      end address of each node group
      integer (kind=kint), allocatable, target :: bc_item(:)
!<      local node ID for node group
!
!      character (len=kchara), allocatable, target :: bc_name(:)
!<      node group name
!
!      integer( kind=kint )  ::  num_bc_smp
!<      number of node group for SMP process
!      integer( kind=kint ), allocatable :: ibc_smp_stack(:)
!<      end address of each node group for SMP process
!
!      integer( kind=kint )  ::  max_bc_4_smp
!<      maximum number of node group for SMP process
!
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine allocate_boundary_data
!
       allocate(nod_grp1%istack_grp(0:nod_grp1%num_grp))
       allocate(nod_grp1%grp_name(nod_grp1%num_grp))
       allocate(bc_item(nod_grp1%num_item))
!
      call clear_boundary_data
!
      end subroutine allocate_boundary_data
!
! ----------------------------------------------------------------------
!
      subroutine clear_boundary_data
!
       nod_grp1%istack_grp=0
       bc_item=0
!
      end subroutine clear_boundary_data
!
! ----------------------------------------------------------------------
!
      subroutine deallocate_boundary_data
!
       deallocate(nod_grp1%istack_grp)
       deallocate(nod_grp1%grp_name)
       deallocate(bc_item)
!
      end subroutine deallocate_boundary_data
!
! ----------------------------------------------------------------------
!
       subroutine allocate_boundary_param_smp
!
       allocate( nod_grp1%istack_grp_smp(0:nod_grp1%num_grp_smp))
       nod_grp1%istack_grp_smp = 0
!
       end subroutine allocate_boundary_param_smp
!
!-----------------------------------------------------------------------
!
       subroutine deallocate_boundary_param_smp
!
       deallocate(nod_grp1%istack_grp_smp)
!
       end subroutine deallocate_boundary_param_smp
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine check_bc_4_sheard_para(my_rank)
!
      integer(kind = kint), intent(in) :: my_rank
!
       write(*,*) 'PE: ', my_rank, 'num_bc ', nod_grp1%num_grp
       write(*,*) 'PE: ', my_rank, 'num_bc_smp ', nod_grp1%num_grp_smp
       write(*,*) 'PE: ', my_rank,                                      &
     &            'ibc_smp_stack ', nod_grp1%istack_grp_smp
!
      end subroutine check_bc_4_sheard_para
!
!-----------------------------------------------------------------------
!
      end module m_node_group
