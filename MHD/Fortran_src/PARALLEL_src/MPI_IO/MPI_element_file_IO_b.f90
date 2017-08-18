!>@file  MPI_element_file_IO_b.f90
!!      module MPI_element_file_IO_b
!!
!!@author  H. Matsui
!!@date Programmed in Aug., 2006
!
!>@brief File IO for element communication table
!!
!!@verbatim
!!      subroutine mpi_input_element_file_b                             &
!!     &         (nprocs_in, my_rank_IO, file_prefix, ele_mesh_IO, ierr)
!!      subroutine mpi_input_surface_file_b                             &
!!     &         (nprocs_in, my_rank_IO, file_prefix, surf_mesh_IO, ierr)
!!      subroutine mpi_input_edge_geometries_b                          &
!!     &         (nprocs_in, my_rank_IO, file_prefix, edge_mesh_IO, ierr)
!!
!!      subroutine mpi_output_element_file_b                            &
!!     &         (nprocs_in, my_rank_IO, ele_mesh_IO)
!!      subroutine mpi_output_surface_file_b                            &
!!     &         (nprocs_in, my_rank_IO, file_prefix, surf_mesh_IO)
!!      subroutine mpi_output_edge_geometries_b                         &
!!     &         (nprocs_in, my_rank_IO, file_prefix, edge_mesh_IO)
!!        type(surf_edge_IO_file), intent(inout) :: ele_mesh_IO
!!        type(surf_edge_IO_file), intent(inout) :: surf_mesh_IO
!!        type(surf_edge_IO_file), intent(inout) :: edge_mesh_IO
!!@endverbatim
!!
!!@param my_rank_IO  MPI rank
!
      module MPI_element_file_IO_b
!
      use m_precision
      use m_machine_parameter
!
      use m_file_format_switch
      use t_calypso_mpi_IO_param
      use t_read_mesh_data
      use set_mesh_file_names
!
      implicit none
!
      type(calypso_MPI_IO_params), save, private :: IO_param
      character(len=kchara), private :: file_name
!
!------------------------------------------------------------------
!
       contains
!
!------------------------------------------------------------------
!
      subroutine mpi_input_element_file_b                               &
     &         (nprocs_in, my_rank_IO, file_prefix, ele_mesh_IO, ierr)
!
      use MPI_element_data_IO_b
!
      character(len=kchara), intent(in) :: file_prefix
      integer(kind = kint), intent(in) :: nprocs_in, my_rank_IO
      type(surf_edge_IO_file), intent(inout) :: ele_mesh_IO
      integer(kind = kint), intent(inout) :: ierr
!
!
      call set_ele_comm_file_name(file_prefix, iflag_single+id_binary_file_fmt,            &
     &    my_rank_IO, file_name)
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &  'Read binary element comm file: ', trim(file_name)
!
      call open_read_mpi_file                                           &
     &   (file_name, nprocs_in, my_rank_IO, IO_param)
!
      call mpi_read_element_comm_table_b(IO_param, ele_mesh_IO%comm)
      call mpi_read_element_geometry_b(IO_param,                        &
     &    ele_mesh_IO%node, ele_mesh_IO%sfed)
      call close_mpi_file(IO_param)
!
      end subroutine mpi_input_element_file_b
!
!------------------------------------------------------------------
!
      subroutine mpi_input_surface_file_b                               &
     &         (nprocs_in, my_rank_IO, file_prefix, surf_mesh_IO, ierr)
!
      use MPI_surface_data_IO_b
!
      character(len=kchara), intent(in) :: file_prefix
      integer(kind = kint), intent(in) :: nprocs_in, my_rank_IO
      type(surf_edge_IO_file), intent(inout) :: surf_mesh_IO
      integer(kind = kint), intent(inout) :: ierr
!
!
      call set_surf_mesh_file_name(file_prefix, iflag_single+id_binary_file_fmt,           &
     &    my_rank_IO, file_name)
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &  'Read binary surface mesh file: ', trim(file_name)
!
      call open_read_mpi_file                                           &
     &   (file_name, nprocs_in, my_rank_IO, IO_param)
!
      call mpi_read_surface_connection_b(IO_param, surf_mesh_IO%comm,   &
     &   surf_mesh_IO%ele, surf_mesh_IO%sfed)
      call mpi_read_surface_geometry_b(IO_param,                        &
     &    surf_mesh_IO%node, surf_mesh_IO%sfed)
      call close_mpi_file(IO_param)
!
      end subroutine mpi_input_surface_file_b
!
!------------------------------------------------------------------
!
      subroutine mpi_input_edge_geometries_b                            &
     &         (nprocs_in, my_rank_IO, file_prefix, edge_mesh_IO, ierr)
!
      use MPI_edge_data_IO_b
!
      character(len=kchara), intent(in) :: file_prefix
      integer(kind = kint), intent(in) :: nprocs_in, my_rank_IO
      type(surf_edge_IO_file), intent(inout) :: edge_mesh_IO
      integer(kind = kint), intent(inout) :: ierr
!
!
      call set_edge_mesh_file_name(file_prefix, iflag_single+id_binary_file_fmt,           &
     &    my_rank_IO, file_name)
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &  'Read binary edge mesh file: ', trim(file_name)
!
      call open_read_mpi_file                                           &
     &   (file_name, nprocs_in, my_rank_IO, IO_param)
!
      call mpi_read_edge_connection_b(IO_param, edge_mesh_IO%comm,      &
     &    edge_mesh_IO%ele, edge_mesh_IO%sfed)
      call mpi_read_edge_geometry_b(IO_param,                           &
     &    edge_mesh_IO%node, edge_mesh_IO%sfed)
      call close_mpi_file(IO_param)
!
      end subroutine mpi_input_edge_geometries_b
!
!------------------------------------------------------------------
!------------------------------------------------------------------
!
      subroutine mpi_output_element_file_b                              &
     &         (nprocs_in, my_rank_IO, file_prefix, ele_mesh_IO)
!
      use MPI_element_data_IO_b
!
      character(len=kchara), intent(in) :: file_prefix
      integer(kind = kint), intent(in) :: nprocs_in, my_rank_IO
      type(surf_edge_IO_file), intent(inout) :: ele_mesh_IO
!
!
      call set_ele_comm_file_name(file_prefix, iflag_single+id_binary_file_fmt,            &
     &    my_rank_IO, file_name)
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &  'Write binary element comm file: ', trim(file_name)
!
      call open_write_mpi_file                                          &
     &   (file_name, nprocs_in, my_rank_IO, IO_param)
!
      call mpi_write_element_comm_table_b(IO_param, ele_mesh_IO%comm)
      call mpi_write_element_geometry_b(IO_param,                       &
     &    ele_mesh_IO%node, ele_mesh_IO%sfed)
      call close_mpi_file(IO_param)
!
      end subroutine mpi_output_element_file_b
!
!------------------------------------------------------------------
!
      subroutine mpi_output_surface_file_b                              &
     &         (nprocs_in, my_rank_IO, file_prefix, surf_mesh_IO)
!
      use MPI_surface_data_IO_b
!
      character(len=kchara), intent(in) :: file_prefix
      integer(kind = kint), intent(in) :: nprocs_in, my_rank_IO
      type(surf_edge_IO_file), intent(inout) :: surf_mesh_IO
!
!
      call set_surf_mesh_file_name(file_prefix, iflag_single+id_binary_file_fmt,           &
     &    my_rank_IO, file_name)
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &  'Write binary surface mesh file: ', trim(file_name)
!
      call open_write_mpi_file                                          &
     &   (file_name, nprocs_in, my_rank_IO, IO_param)
!
      call mpi_write_surface_connection_b(IO_param, surf_mesh_IO%comm,  &
     &   surf_mesh_IO%ele, surf_mesh_IO%sfed)
      call mpi_write_surface_geometry_b(IO_param,                       &
     &    surf_mesh_IO%node, surf_mesh_IO%sfed)
      call close_mpi_file(IO_param)
!
      end subroutine mpi_output_surface_file_b
!
!------------------------------------------------------------------
!
      subroutine mpi_output_edge_geometries_b                           &
     &         (nprocs_in, my_rank_IO, file_prefix, edge_mesh_IO)
!
      use MPI_edge_data_IO_b
!
      character(len=kchara), intent(in) :: file_prefix
      integer(kind = kint), intent(in) :: nprocs_in, my_rank_IO
      type(surf_edge_IO_file), intent(inout) :: edge_mesh_IO
!
!
      call set_edge_mesh_file_name(file_prefix, iflag_single+id_binary_file_fmt,           &
     &    my_rank_IO, file_name)
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &  'Write binary edge mesh file: ', trim(file_name)
!
      call open_write_mpi_file                                          &
     &   (file_name, nprocs_in, my_rank_IO, IO_param)
!
      call mpi_write_edge_connection_b(IO_param, edge_mesh_IO%comm,     &
     &   edge_mesh_IO%ele, edge_mesh_IO%sfed)
      call mpi_write_edge_geometry_b(IO_param,                          &
     &   edge_mesh_IO%node, edge_mesh_IO%sfed)
      call close_mpi_file(IO_param)
!
      end subroutine mpi_output_edge_geometries_b
!
!------------------------------------------------------------------
!
      end module MPI_element_file_IO_b