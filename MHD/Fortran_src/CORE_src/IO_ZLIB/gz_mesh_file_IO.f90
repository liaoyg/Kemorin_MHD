!>@file   gz_mesh_file_IO.f90
!!@brief  module gz_mesh_file_IO
!!
!!@author H.Matsui
!!@date     Programmed by H.Matsui in Apr., 2006
!
!>@brief  Mesh file IO for gxipped format
!!
!!@verbatim
!!      subroutine gz_read_mesh(my_rank_IO, file_name, fem_IO, ierr)
!!        type(mesh_data), intent(inout) :: fem_IO
!!
!!      subroutine gz_read_mesh_geometry                                &
!!     &         (my_rank_IO, file_name, mesh_IO, ierr)
!!      subroutine gz_read_node_size                                    &
!!     &         (my_rank_IO, file_name, mesh_IO, ierr)
!!      subroutine gz_read_geometry_size                                &
!!     &         (my_rank_IO, file_name, mesh_IO, ierr)
!!        type(mesh_geometry), intent(inout) :: mesh_IO
!!
!!      subroutine gz_write_mesh_file(my_rank_IO, file_name, fem_IO)
!!        type(mesh_data), intent(inout) :: fem_IO
!!@endverbatim
!!
      module gz_mesh_file_IO
!
      use m_precision
      use m_machine_parameter
!
      use t_mesh_data
      use skip_gz_comment
      use gz_mesh_data_IO
!
      implicit none
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine gz_read_mesh(my_rank_IO, file_name, fem_IO, ierr)
!
      integer(kind = kint), intent(in) :: my_rank_IO
      character(len=kchara), intent(in) :: file_name
!
      type(mesh_data), intent(inout) :: fem_IO
      integer(kind = kint), intent(inout) :: ierr
!
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &   'Read gzipped mesh file: ', trim(file_name)
!
      call open_rd_gzfile_f(file_name)
!
      call gz_read_geometry_data(my_rank_IO, fem_IO%mesh, ierr)
      call gz_read_mesh_groups(fem_IO%group)
!
      call close_gzfile_f
!
      end subroutine gz_read_mesh
!
!  ---------------------------------------------------------------------
!
      subroutine gz_read_mesh_geometry                                  &
     &         (my_rank_IO, file_name, mesh_IO, ierr)
!
      integer(kind = kint), intent(in) :: my_rank_IO
      character(len=kchara), intent(in) :: file_name
!
      type(mesh_geometry), intent(inout) :: mesh_IO
      integer(kind = kint), intent(inout) :: ierr
!
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &   'Read gzipped mesh file: ', trim(file_name)
!
      call open_rd_gzfile_f(file_name)
      call gz_read_geometry_data(my_rank_IO, mesh_IO, ierr)
      call close_gzfile_f
!
      end subroutine gz_read_mesh_geometry
!
!  ---------------------------------------------------------------------
!
      subroutine gz_read_node_size                                      &
     &         (my_rank_IO, file_name, mesh_IO, ierr)
!
      integer(kind = kint), intent(in) :: my_rank_IO
      character(len=kchara), intent(in) :: file_name
!
      type(mesh_geometry), intent(inout) :: mesh_IO
      integer(kind = kint), intent(inout) :: ierr
!
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &    'Read gzipped mesh file: ', trim(file_name)
!
      call open_rd_gzfile_f(file_name)
!
      call gz_read_num_node(my_rank_IO, mesh_IO, ierr)
      call close_gzfile_f
!
      end subroutine gz_read_node_size
!
!------------------------------------------------------------------
!
      subroutine gz_read_geometry_size                                  &
     &         (my_rank_IO, file_name, mesh_IO, ierr)
!
      integer(kind = kint), intent(in) :: my_rank_IO
      character(len=kchara), intent(in) :: file_name
!
      type(mesh_geometry), intent(inout) :: mesh_IO
      integer(kind = kint), intent(inout) :: ierr
!
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &   'Read gzipped mesh file: ', trim(file_name)
!
      call open_rd_gzfile_f(file_name)
!
      call gz_read_num_node_ele(my_rank_IO, mesh_IO, ierr)
      call close_gzfile_f
!
      end subroutine gz_read_geometry_size
!
!------------------------------------------------------------------
!------------------------------------------------------------------
!
      subroutine gz_write_mesh_file(my_rank_IO, file_name, fem_IO)
!
      integer(kind = kint), intent(in) :: my_rank_IO
      character(len=kchara), intent(in) :: file_name
!
      type(mesh_data), intent(inout) :: fem_IO
!
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &   'Write gzipped mesh file: ', trim(file_name)
!
      call open_wt_gzfile_f(file_name)
!
      call gz_write_geometry_data(my_rank_IO, fem_IO%mesh)
!
      call gz_write_mesh_groups(fem_IO%group)
!
      call close_gzfile_f
!
      end subroutine gz_write_mesh_file
!
!  ---------------------------------------------------------------------
!
      end module gz_mesh_file_IO
