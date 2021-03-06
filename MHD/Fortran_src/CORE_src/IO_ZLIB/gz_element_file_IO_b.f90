!>@file  gz_element_file_IO_b.f90
!!      module gz_element_file_IO_b
!!
!!@author  H. Matsui
!!@date Programmed in Aug., 2006
!
!>@brief File IO for element communication table
!!
!!@verbatim
!!      subroutine gz_input_element_file_b                              &
!!     &         (my_rank_IO, file_name, ele_mesh_IO, ierr)
!!      subroutine gz_input_surface_file_b                              &
!!     &         (my_rank_IO, file_name, surf_mesh_IO, ierr)
!!      subroutine gz_input_edge_file_b                                 &
!!     &         (my_rank_IO, file_name, edge_mesh_IO, ierr)
!!
!!      subroutine gz_output_element_file_b                             &
!!     &         (my_rank_IO, ele_mesh_IO)
!!      subroutine gz_output_surface_file_b                             &
!!     &         (my_rank_IO, file_name, surf_mesh_IO)
!!      subroutine gz_output_edge_file_b                                &
!!     &         (my_rank_IO, file_name, edge_mesh_IO)
!!        type(surf_edge_IO_file), intent(inout) :: ele_mesh_IO
!!        type(surf_edge_IO_file), intent(inout) :: surf_mesh_IO
!!        type(surf_edge_IO_file), intent(inout) :: edge_mesh_IO
!!@endverbatim
!!
!!@param my_rank_IO  MPI rank
!
      module gz_element_file_IO_b
!
      use m_precision
      use m_machine_parameter
!
      use m_file_format_switch
      use t_read_mesh_data
      use set_mesh_file_names
      use skip_gz_comment
!
      implicit none
!
!------------------------------------------------------------------
!
       contains
!
!------------------------------------------------------------------
!
      subroutine gz_input_element_file_b                                &
     &         (my_rank_IO, file_name, ele_mesh_IO, ierr)
!
      use gz_element_data_IO_b
!
      character(len=kchara), intent(in) :: file_name
      integer(kind = kint), intent(in) :: my_rank_IO
      type(surf_edge_IO_file), intent(inout) :: ele_mesh_IO
      integer(kind = kint), intent(inout) :: ierr
!
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &  'Read gzipped binary element comm file: ', trim(file_name)
!
      call open_rd_gzfile_f(file_name)
      call gz_read_element_comm_table_b                                 &
     &   (my_rank_IO, ele_mesh_IO%comm, ierr)
!      call gz_read_element_geometry_b                                  &
!     &   (ele_mesh_IO%node, ele_mesh_IO%sfed)
      call close_gzfile_f
!
      end subroutine gz_input_element_file_b
!
!------------------------------------------------------------------
!
      subroutine gz_input_surface_file_b                                &
     &         (my_rank_IO, file_name, surf_mesh_IO, ierr)
!
      use gz_surface_data_IO_b
!
      character(len=kchara), intent(in) :: file_name
      integer(kind = kint), intent(in) :: my_rank_IO
      type(surf_edge_IO_file), intent(inout) :: surf_mesh_IO
      integer(kind = kint), intent(inout) :: ierr
!
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &  'Read gzipped binary surface mesh file: ', trim(file_name)
!
      call open_rd_gzfile_f(file_name)
      call gz_read_surface_connection_b(my_rank_IO, surf_mesh_IO%comm,  &
     &   surf_mesh_IO%ele, surf_mesh_IO%sfed, ierr)
!      call gz_read_surface_geometry_b                                  &
!     &   (surf_mesh_IO%node, surf_mesh_IO%sfed)
      call close_gzfile_f
!
      end subroutine gz_input_surface_file_b
!
!------------------------------------------------------------------
!
      subroutine gz_input_edge_file_b                                   &
     &         (my_rank_IO, file_name, edge_mesh_IO, ierr)
!
      use gz_edge_data_IO_b
!
      character(len=kchara), intent(in) :: file_name
      integer(kind = kint), intent(in) :: my_rank_IO
      type(surf_edge_IO_file), intent(inout) :: edge_mesh_IO
      integer(kind = kint), intent(inout) :: ierr
!
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &  'Read gzipped binary edge mesh file: ', trim(file_name)
!
      call open_rd_gzfile_f(file_name)
      call gz_read_edge_connection_b(my_rank_IO, edge_mesh_IO%comm,     &
     &    edge_mesh_IO%ele, edge_mesh_IO%sfed, ierr)
!      call gz_read_edge_geometry_b                                     &
!     &   (edge_mesh_IO%node, edge_mesh_IO%sfed)
      call close_gzfile_f
!
      end subroutine gz_input_edge_file_b
!
!------------------------------------------------------------------
!------------------------------------------------------------------
!
      subroutine gz_output_element_file_b                               &
     &         (my_rank_IO, file_name, ele_mesh_IO)
!
      use gz_element_data_IO_b
!
      character(len=kchara), intent(in) :: file_name
      integer(kind = kint), intent(in) :: my_rank_IO
      type(surf_edge_IO_file), intent(inout) :: ele_mesh_IO
!
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &  'Write gzipped binary element comm file: ', trim(file_name)
!
      call open_wt_gzfile_f(file_name)
      call gz_write_element_comm_table_b(my_rank_IO, ele_mesh_IO%comm)
!      call gz_write_element_geometry_b                                 &
!     &   (ele_mesh_IO%node, ele_mesh_IO%sfed)
      call close_gzfile_f
!
      end subroutine gz_output_element_file_b
!
!------------------------------------------------------------------
!
      subroutine gz_output_surface_file_b                               &
     &         (my_rank_IO, file_name, surf_mesh_IO)
!
      use gz_surface_data_IO_b
!
      character(len=kchara), intent(in) :: file_name
      integer(kind = kint), intent(in) :: my_rank_IO
      type(surf_edge_IO_file), intent(inout) :: surf_mesh_IO
!
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &  'Write gzipped binary surface mesh file: ', trim(file_name)
!
      call open_wt_gzfile_f(file_name)
      call gz_write_surface_connection_b(my_rank_IO, surf_mesh_IO%comm, &
     &   surf_mesh_IO%ele, surf_mesh_IO%sfed)
!      call gz_write_surface_geometry_b                                 &
!     &   (surf_mesh_IO%node, surf_mesh_IO%sfed)
      call close_gzfile_f
!
      end subroutine gz_output_surface_file_b
!
!------------------------------------------------------------------
!
      subroutine gz_output_edge_file_b                                  &
     &         (my_rank_IO, file_name, edge_mesh_IO)
!
      use gz_edge_data_IO_b
!
      character(len=kchara), intent(in) :: file_name
      integer(kind = kint), intent(in) :: my_rank_IO
      type(surf_edge_IO_file), intent(inout) :: edge_mesh_IO
!
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &  'Write gzipped binary edge mesh file: ', trim(file_name)
!
      call open_wt_gzfile_f(file_name)
      call gz_write_edge_connection_b(my_rank_IO, edge_mesh_IO%comm,    &
     &   edge_mesh_IO%ele, edge_mesh_IO%sfed)
!      call gz_write_edge_geometry_b                                    &
!     &   (edge_mesh_IO%node, edge_mesh_IO%sfed)
      call close_gzfile_f
!
      end subroutine gz_output_edge_file_b
!
!------------------------------------------------------------------
!
      end module gz_element_file_IO_b
