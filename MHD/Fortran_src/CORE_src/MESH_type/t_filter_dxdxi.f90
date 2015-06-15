!>@file  t_filter_dxdxi.f90
!!       module t_filter_dxdxi
!!
!!@author H. Matsui
!!@date   Programmed on March, 2009
!!@n      Modified by H. Matsui on Feb., 2012
!
!> @brief Strucure for dxdxi for SGS model
!!
!!@verbatim
!!      subroutine alloc_jacobians_on_node_type(nnod_filter_mom, dxi_nod)
!!        integer (kind = kint), intent(in) :: nnod_filter_mom
!!        type(dxdxi_direction_type), intent(inout) :: dxi_nod
!!      subroutine alloc_jacobians_ele_type(nele_filter_mom, dxi_ele)
!!        integer (kind = kint), intent(in) :: nele_filter_mom
!!        type(dxdxi_direction_type), intent(inout) :: dxi_ele
!!      subroutine dealloc_jacobians_on_node_type(dxi_nod)
!!        type(dxdxi_direction_type), intent(inout) :: dxi_nod
!!      subroutine dealloc_jacobians_ele_type(dxi_ele)
!!        type(dxdxi_direction_type), intent(inout) :: dxi_ele
!!
!!          1st difference of elengh_nod
!!              (node ID, direction of diffrence)
!!      dxdxi_nod(:)... dxi_nod%dx%df_dxi(:)
!!      dxdei_nod(:)... dxi_nod%dx%df_dei(:)
!!      dxdzi_nod(:)... dxi_nod%dx%df_dzi(:)
!!      dydxi_nod(:)... dxi_nod%dy%df_dxi(:)
!!      dydei_nod(:)... dxi_nod%dy%df_dei(:)
!!      dydzi_nod(:)... dxi_nod%dy%df_dzi(:)
!!      dzdxi_nod(:)... dxi_nod%dz%df_dxi(:)
!!      dzdei_nod(:)... dxi_nod%dz%df_dei(:)
!!      dzdzi_nod(:)... dxi_nod%dz%df_dzi(:)
!!
!!          1st difference of elengh_nod
!!              (node ID, direction of diffrence)
!!      dxdxi_ele(:)... dxi_ele%dx%df_dxi(:)
!!      dxdei_ele(:)... dxi_ele%dx%df_dei(:)
!!      dxdzi_ele(:)... dxi_ele%dx%df_dzi(:)
!!      dydxi_ele(:)... dxi_ele%dy%df_dxi(:)
!!      dydei_ele(:)... dxi_ele%dy%df_dei(:)
!!      dydzi_ele(:)... dxi_ele%dy%df_dzi(:)
!!      dzdxi_ele(:)... dxi_ele%dz%df_dxi(:)
!!      dzdei_ele(:)... dxi_ele%dz%df_dei(:)
!!      dzdzi_ele(:)... dxi_ele%dz%df_dzi(:)
!!@endverbatim
!
module t_filter_dxdxi
!
      use m_precision
!
      implicit none
!
!
!>        Structure for dxi/dxi
      type dxdxi_diff_type
!>        delivative in \xi-diretion
        real(kind=kreal), pointer :: df_dxi(:)
!>        delivative in \eta-diretion
        real(kind=kreal), pointer :: df_dei(:)
!>        delivative in \zeta-diretion
        real(kind=kreal), pointer :: df_dzi(:)
      end  type dxdxi_diff_type
!
!>        Structure for dx/dxi
      type dxdxi_direction_type
!>        delivative for x (dx/dxi)
        type(dxdxi_diff_type) :: dx
!>        delivative for y (dy/dxi)
        type(dxdxi_diff_type) :: dy
!>        delivative for z (dz/dxi)
        type(dxdxi_diff_type) :: dz
      end type dxdxi_direction_type
!
!>        Structure for dx/dxi at nodes and elements
      type gradient_model_data_type
!>        Structure for dx/dxi at nodes
        type(dxdxi_direction_type) :: dxi_nod
!>        Structure for dx/dxi in elements
        type(dxdxi_direction_type) :: dxi_ele
      end type gradient_model_data_type
!
      private :: alloc_dxi_diff_type, dealloc_dxi_diff_type
!      private :: alloc_dxdxi_diff_type, dealloc_dxdxi_diff_type
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine alloc_dxi_diff_type(num_fmom, dx)
!
      integer (kind = kint), intent(in) :: num_fmom
      type(dxdxi_diff_type), intent(inout) :: dx
!
!
      allocate(dx%df_dxi(num_fmom))
      allocate(dx%df_dei(num_fmom))
      allocate(dx%df_dzi(num_fmom))
!
      if(num_fmom .gt. 0) then
        dx%df_dxi = 0.0d0
        dx%df_dei = 0.0d0
        dx%df_dzi = 0.0d0
      end if
!
      end subroutine alloc_dxi_diff_type
!
!  ---------------------------------------------------------------------
!
      subroutine dealloc_dxi_diff_type(dx)
!
      type(dxdxi_diff_type), intent(inout) :: dx
!
!
      deallocate(dx%df_dxi, dx%df_dei, dx%df_dzi)
!
      end subroutine dealloc_dxi_diff_type
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine alloc_dxdxi_diff_type(num_fmom, dxi)
!
      integer (kind = kint), intent(in) :: num_fmom
      type(dxdxi_direction_type), intent(inout) :: dxi
!
!
      call alloc_dxi_diff_type(num_fmom, dxi%dx)
      call alloc_dxi_diff_type(num_fmom, dxi%dy)
      call alloc_dxi_diff_type(num_fmom, dxi%dz)
!
      end subroutine alloc_dxdxi_diff_type
!
!  ---------------------------------------------------------------------
!
      subroutine dealloc_dxdxi_diff_type(dxi)
!
      type(dxdxi_direction_type), intent(inout) :: dxi
!
!
      call dealloc_dxi_diff_type(dxi%dx)
      call dealloc_dxi_diff_type(dxi%dy)
      call dealloc_dxi_diff_type(dxi%dz)
!
      end subroutine dealloc_dxdxi_diff_type
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine alloc_jacobians_on_node_type(nnod_filter_mom, dxi_nod)
!
      integer (kind = kint), intent(in) :: nnod_filter_mom
      type(dxdxi_direction_type), intent(inout) :: dxi_nod
!
!
      call alloc_dxdxi_diff_type(nnod_filter_mom, dxi_nod)
!
      end subroutine alloc_jacobians_on_node_type
!
!  ---------------------------------------------------------------------
!
      subroutine alloc_jacobians_ele_type(nele_filter_mom, dxi_ele)
!
      integer (kind = kint), intent(in) :: nele_filter_mom
      type(dxdxi_direction_type), intent(inout) :: dxi_ele
!
!
      call alloc_dxdxi_diff_type(nele_filter_mom, dxi_ele)
!
      end subroutine alloc_jacobians_ele_type
!
!  ---------------------------------------------------------------------
!
      subroutine dealloc_jacobians_on_node_type(dxi_nod)
!
      type(dxdxi_direction_type), intent(inout) :: dxi_nod
!
!
      call dealloc_dxdxi_diff_type(dxi_nod)
!
      end subroutine dealloc_jacobians_on_node_type
!
!  ---------------------------------------------------------------------
!
      subroutine dealloc_jacobians_ele_type(dxi_ele)
!
      type(dxdxi_direction_type), intent(inout) :: dxi_ele
!
!
      call dealloc_dxdxi_diff_type(dxi_ele)
!
      end subroutine dealloc_jacobians_ele_type
!
!  ---------------------------------------------------------------------
!
      end module t_filter_dxdxi
