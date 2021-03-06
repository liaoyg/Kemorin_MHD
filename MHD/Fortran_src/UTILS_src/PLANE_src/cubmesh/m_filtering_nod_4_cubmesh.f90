!
!      module m_filtering_nod_4_cubmesh
!
!     Written by H. Matsui
!
!       subroutine allocate_work_4_filter_nod
!       subroutine reset_work_4_filter_nod
!       subroutine allocate_filters_nod(nf_type)
!       subroutine deallocate_work_4_filter_nod
!       subroutine deallocate_filters_nod
!
      module m_filtering_nod_4_cubmesh
!
      use m_precision
!
      implicit none
!
!
      integer (kind = kint) :: inod_f_total
!
      integer (kind = kint), dimension (:), allocatable :: inod_f
      integer (kind = kint), dimension (:), allocatable :: inod_f_stack
      integer (kind = kint), dimension (:), allocatable :: inod_f_item
      integer (kind = kint), dimension (:,:), allocatable               &
     &      :: inod_f_dist
      real (kind = kreal), dimension (:,:), allocatable :: filter_c
!
!
      integer (kind = kint), dimension (:,:,:), allocatable             &
     &      :: nnod_neib_x
      integer (kind = kint), dimension (:,:,:,:), allocatable           &
     &      :: inod_f_item_x
      integer (kind = kint), dimension (:,:,:,:), allocatable           &
     &      :: inod_f_dist_x
!
      integer (kind = kint), dimension (:,:,:), allocatable             &
     &      :: nnod_neib_y
      integer (kind = kint), dimension (:,:,:,:), allocatable           &
     &      :: inod_f_item_y
      integer (kind = kint), dimension (:,:,:,:), allocatable           &
     &      :: inod_f_dist_y
!
      integer (kind = kint), dimension (:,:,:), allocatable             &
     &      :: nnod_neib_z
      integer (kind = kint), dimension (:,:,:,:), allocatable           &
     &      :: inod_f_item_z
      integer (kind = kint), dimension (:,:,:,:), allocatable           &
     &      :: inod_f_dist_z
!
      integer (kind = kint), dimension (:,:,:), allocatable             &
     &      :: nnod_neib_xy
      integer (kind = kint), dimension (:,:,:,:,:), allocatable         &
     &      :: inod_f_item_xy
      integer (kind = kint), dimension (:,:,:,:,:), allocatable         &
     &      :: inod_f_dist_xy
!
      integer (kind = kint), dimension (:,:,:), allocatable             &
     &      :: nnod_neib
      integer (kind = kint), dimension (:,:,:,:,:), allocatable         &
     &      :: inod_f_item_3d
      integer (kind = kint), dimension (:,:,:,:,:), allocatable         &
     &      :: inod_f_dist_3d
!
!
      integer (kind = kint) :: iflag_z_filter
      real(kind = kreal)    :: eps_filter
!
      real(kind = kreal), dimension(:,:), allocatable :: f_mom_1d
      real(kind = kreal), dimension(:,:), allocatable :: df_mom_1d
!
      real(kind = kreal), dimension(:,:,:,:,:), allocatable             &
     &      :: filter_c_x, filter_c_y, filter_c_z
      real(kind = kreal), dimension(:,:,:,:,:), allocatable             &
     &      :: filter_c_xy, filter_c_3d
!
      integer (kind = kint), dimension (:,:), allocatable               &
     &      :: inod_f_dist_new
!
!  ----------------------------------------------------------------------
!
      contains
!
!  ----------------------------------------------------------------------
!
       subroutine allocate_work_4_filter_nod
!
       use m_size_of_cube
       use m_comm_data_cube_kemo
!
!
       allocate( nnod_neib_x (numnod_x,numnod_y,numnod_z)    )
       allocate( nnod_neib_y (numnod_x,numnod_y,numnod_z)    )
       allocate( nnod_neib_z (numnod_x,numnod_y,numnod_z)    )
       allocate( nnod_neib_xy(numnod_x,numnod_y,numnod_z)    )
       allocate( nnod_neib   (numnod_x,numnod_y,numnod_z)    )
!
       allocate( inod_f_item_x(ndep_1,numnod_x,numnod_y,numnod_z)    )
       allocate( inod_f_dist_x(ndep_1,numnod_x,numnod_y,numnod_z)    )
       allocate( inod_f_item_y(ndep_1,numnod_x,numnod_y,numnod_z)    )
       allocate( inod_f_dist_y(ndep_1,numnod_x,numnod_y,numnod_z)    )
       allocate( inod_f_item_z(ndep_1,numnod_x,numnod_y,numnod_z)    )
       allocate( inod_f_dist_z(ndep_1,numnod_x,numnod_y,numnod_z)    )
!
       allocate( inod_f_item_xy(ndep_2,numnod_x,numnod_y,numnod_z,2) )
       allocate( inod_f_dist_xy(ndep_2,numnod_x,numnod_y,numnod_z,2) )
!
       allocate( inod_f_item_3d(ndep_3,numnod_x,numnod_y,numnod_z,3) )
       allocate( inod_f_dist_3d(ndep_3,numnod_x,numnod_y,numnod_z,3) )
!
       call reset_work_4_filter_nod
!
       end subroutine allocate_work_4_filter_nod
!
!  ----------------------------------------------------------------------
!
       subroutine reset_work_4_filter_nod
!
       nnod_neib =    0
       nnod_neib_x =  0
       nnod_neib_y =  0
       nnod_neib_z =  0
       nnod_neib_xy = 0
!
       inod_f_item_x =  0
       inod_f_dist_x =  0
       inod_f_item_y =  0
       inod_f_dist_y =  0
       inod_f_item_z =  0
       inod_f_dist_z =  0
       inod_f_item_xy = 0
       inod_f_dist_xy = 0
       inod_f_item_3d = 0
       inod_f_dist_3d = 0
!
       end subroutine reset_work_4_filter_nod
!
!  ----------------------------------------------------------------------
!
       subroutine allocate_filters_nod(nf_type)
!
       use m_size_of_cube
       use m_comm_data_cube_kemo
!
       integer(kind = kint), intent(in) :: nf_type
!
       allocate(filter_c_x(ndep_1,numnod_x,numnod_y,numnod_z,nf_type))
       allocate(filter_c_y(ndep_1,numnod_x,numnod_y,numnod_z,nf_type))
       allocate(filter_c_z(ndep_1,numnod_x,numnod_y,numnod_z,nf_type))
       allocate(filter_c_xy(ndep_2,numnod_x,numnod_y,numnod_z,nf_type))
       allocate(filter_c_3d(ndep_3,numnod_x,numnod_y,numnod_z,nf_type))
!
       filter_c_x = 0.0d0
       filter_c_y = 0.0d0
       filter_c_z = 0.0d0
       filter_c_xy = 0.0d0
       filter_c_3d = 0.0d0
!
       end subroutine allocate_filters_nod
!
!  ----------------------------------------------------------------------
!
       subroutine deallocate_work_4_filter_nod
!
       deallocate( nnod_neib_x     )
       deallocate( nnod_neib_y     )
       deallocate( nnod_neib_z     )
       deallocate( nnod_neib_xy    )
       deallocate( nnod_neib       )
!
       deallocate( inod_f_item_x  )
       deallocate( inod_f_dist_x  )
       deallocate( inod_f_item_y  )
       deallocate( inod_f_dist_y  )
       deallocate( inod_f_item_z  )
       deallocate( inod_f_dist_z  )
!
       deallocate( inod_f_item_xy )
       deallocate( inod_f_dist_xy )
!
       deallocate( inod_f_item_3d )
       deallocate( inod_f_dist_3d )
!
       end subroutine deallocate_work_4_filter_nod
!
!  ----------------------------------------------------------------------
!
       subroutine deallocate_filters_nod
!
       deallocate (filter_c_x)
       deallocate (filter_c_y)
       deallocate (filter_c_z)
       deallocate (filter_c_xy)
       deallocate (filter_c_3d)
!
       end subroutine deallocate_filters_nod
!
!  ----------------------------------------------------------------------
!
      end module m_filtering_nod_4_cubmesh
