!>@file   t_geometry_data.f90
!!@brief  module t_geometry_data
!!
!!@author  H. Matsui
!!@date Programmed in 2008
!
!>@brief structure of geometry data for FEM mesh
!!  including node and element position, connectivities
!!
!!@verbatim
!!      subroutine alloc_numnod_stack(nprocs, node)
!!      subroutine alloc_numele_stack(nprocs, ele)
!!      subroutine allocate_node_geometry_type(node)
!!      subroutine alloc_node_geometry_base(node)
!!      subroutine allocate_sph_node_geometry(node)
!!
!!      subroutine allocate_ele_connect_type(ele)
!!      subroutine alloc_element_types(ele)
!!      subroutine alloc_ele_connectivity(ele)
!!      subroutine alloc_overlaped_ele(ele)
!!      subroutine alloc_ele_geometry(ele)
!!      subroutine allocate_node_param_smp_type(node)
!!      subroutine allocate_ele_param_smp_type(ele)
!!        type(element_data), intent(inout) :: ele
!!
!!      subroutine dealloc_numnod_stack(node)
!!      subroutine dealloc_numele_stack(ele)
!!      subroutine deallocate_node_geometry_type(node)
!!      subroutine dealloc_node_geometry_base(node)
!!      subroutine deallocate_sph_node_geometry(node)
!!
!!      subroutine deallocate_ele_connect_type(ele)
!!      subroutine dealloc_overlaped_ele(ele)
!!      subroutine deallocate_ele_geometry_type(ele)
!!      subroutine deallocate_node_param_smp_type(node)
!!      subroutine deallocate_ele_param_smp_type(ele)
!!        type(element_data), intent(inout) :: ele
!!
!!      subroutine check_nod_size_smp_type(node, my_rank)
!!      subroutine check_ele_size_smp_type(ele, my_rank)
!!@endverbatim
!
      module t_noise_node_data
!
      use m_precision
!
      implicit  none
!
!
!>  structure for noise node data (position)
      type noise_node
!>        node value (current level noise value)
        real( kind=kreal )  ::  n_value
!>        dimesion size of sub node, eg(2, total sub node 2*2*2 = 8)
        integer( kind=kint )  ::  node_dim
!>        node hierarchical level
        integer( kind=kint )  ::  node_level
!>        subnode list for next level noise node
        type( noise_node ), dimension(:), pointer ::  sub_node
!
      end type noise_node
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine alloc_noise_node(n_node, dim, level)
!
      type(noise_node), intent(inout) :: n_node
      integer(kind=kint), intent(in) :: level, dim
      integer(kind=kint) :: i, size
      real(kind=kreal) :: rand_v

      rand_v = rand(1)
      if(rand_v .gt. 0.5) then
        n_node%n_value = 1.0
      else
        n_node%n_value = 0.0
      end if
      n_node%node_dim = dim
      n_node%node_level = level
      if(level .gt. 0) then
        size = dim*dim*dim
        allocate(n_node%sub_node(size))
        do i = 1, size
          call alloc_noise_node(n_node%sub_node(i), dim, level-1)
        end do
      end if

      end subroutine alloc_noise_node
!
!  ---------------------------------------------------------------------
!
      subroutine dealloc_noise_node(n_node)
!
      type(noise_node), intent(inout) :: n_node
      integer(kind=kint) :: i, size, dim

      if(n_node%node_level .gt. 0) then
        dim = n_node%node_dim
        size = dim*dim*dim
        do i = 1, size
          call dealloc_noise_node(n_node%sub_node(i))
        end do
        deallocate(n_node%sub_node)
      end if

      end subroutine dealloc_noise_node
!
!-----------------------------------------------------------------------
!
      end module t_noise_node_data
