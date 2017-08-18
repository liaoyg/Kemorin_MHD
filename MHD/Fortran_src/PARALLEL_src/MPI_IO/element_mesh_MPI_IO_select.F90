!>@file   element_mesh_MPI_IO_select.f90
!!@brief  module element_mesh_MPI_IO_select
!!
!!@author H.Matsui
!!@date     Programmed by H.Matsui in Apr., 2006
!
!>@brief  Choose mesh file to read
!!
!!@verbatim
!!      subroutine sel_mpi_read_ele_mesh(mesh_file, fem_IO)
!!      subroutine sel_mpi_read_surf_mesh(mesh_file, fem_IO)
!!      subroutine sel_mpi_read_edge_mesh(mesh_file, fem_IO)
!!        type(field_IO_params), intent(in) ::  mesh_file
!!        type(mesh_data), intent(inout) :: fem_IO
!!
!!      subroutine sel_mpi_write_ele_mesh_file(mesh_file, fem_IO)
!!      subroutine sel_mpi_write_surf_mesh_file(mesh_file, fem_IO)
!!      subroutine sel_mpi_write_edge_mesh_file(mesh_file, fem_IO)
!!        type(field_IO_params), intent(in) ::  mesh_file
!!        type(mesh_data), intent(inout) :: fem_IO
!!@endverbatim
!
      module element_mesh_MPI_IO_select
!
      use m_precision
      use calypso_mpi
!
      use m_file_format_switch
      use t_file_IO_parameter
      use t_mesh_data
!
      use mesh_IO_select
!
      use MPI_mesh_file_IO
      use MPI_mesh_file_IO_b
      use set_mesh_file_names
!
#ifdef ZLIB_IO
      use gz_MPI_mesh_file_IO
      use gz_MPI_mesh_file_IO_b
#endif
!
      implicit none
!
      character(len=kchara), private :: file_name
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine sel_mpi_read_ele_mesh(mesh_file, fem_IO)
!
      type(field_IO_params), intent(in) ::  mesh_file
      type(mesh_data), intent(inout) :: fem_IO
!
      integer(kind = kint) :: ierr = 0
!
!
      call set_mesh_file_name_by_param(mesh_file, my_rank, file_name)
!
      if(mesh_file%iflag_format                                         &
     &     .eq. iflag_single+id_binary_file_fmt) then
        call mpi_input_element_file_b                                   &
     &     (nprocs, my_rank, file_name, fem_IO)
      else if(mesh_file%iflag_format .eq. iflag_single) then
        call mpi_input_element_file                                     &
     &     (nprocs, my_rank, file_name, fem_IO)
!
#ifdef ZLIB_IO
      else if(mesh_file%iflag_format                                    &
     &        .eq. iflag_single+id_gzip_bin_file_fmt) then
        call gz_mpi_input_element_file_b                                &
     &     (nprocs, my_rank, file_name, fem_IO)
      else if(mesh_file%iflag_format                                    &
     &        .eq. iflag_single+id_gzip_txt_file_fmt) then
        call gz_mpi_input_element_file                                  &
     &     (nprocs, my_rank, file_name, fem_IO)
#endif
!
      else
        call sel_read_ele_mesh(mesh_file, my_rank, fem_IO, ierr)
      end if 
!
      if(ierr .gt. 0) then
        call calypso_mpi_abort(ierr, 'Mesh data is wrong!!')
      end if
!
      end subroutine sel_mpi_read_ele_mesh
!
!------------------------------------------------------------------
!
      subroutine sel_mpi_read_surf_mesh(mesh_file, fem_IO)
!
      type(field_IO_params), intent(in) ::  mesh_file
      type(mesh_data), intent(inout) :: fem_IO
!
      integer(kind = kint) :: ierr = 0
!
!
      call set_mesh_file_name_by_param(mesh_file, my_rank, file_name)
!
      if(mesh_file%iflag_format                                         &
     &     .eq. iflag_single+id_binary_file_fmt) then
        call mpi_input_surface_file_b                                   &
     &     (nprocs, my_rank, file_name, fem_IO)
      else if(mesh_file%iflag_format .eq. iflag_single) then
        call mpi_input_surface_file                                     &
     &     (nprocs, my_rank, file_name, fem_IO)
!
#ifdef ZLIB_IO
      else if(mesh_file%iflag_format                                    &
     &        .eq. iflag_single+id_gzip_bin_file_fmt) then
        call gz_mpi_input_surface_file_b                                &
     &     (nprocs, my_rank, file_name, fem_IO)
      else if(mesh_file%iflag_format                                    &
     &        .eq. iflag_single+id_gzip_txt_file_fmt) then
        call gz_mpi_input_surface_file                                  &
     &     (nprocs, my_rank, file_name, fem_IO)
#endif
!
      else
        call sel_read_surf_mesh(mesh_file, my_rank, fem_IO, ierr)
      end if 
!
      if(ierr .gt. 0) then
        call calypso_mpi_abort(ierr, 'Mesh data is wrong!!')
      end if
!
      end subroutine sel_mpi_read_surf_mesh
!
!------------------------------------------------------------------
!
      subroutine sel_mpi_read_edge_mesh(mesh_file, fem_IO)
!
      type(field_IO_params), intent(in) ::  mesh_file
      type(mesh_data), intent(inout) :: fem_IO
!
      integer(kind = kint) :: ierr = 0
!
!
      call set_mesh_file_name_by_param(mesh_file, my_rank, file_name)
!
      if(mesh_file%iflag_format                                         &
     &     .eq. iflag_single+id_binary_file_fmt) then
        call mpi_input_edge_file_b                                      &
     &     (nprocs, my_rank, file_name, fem_IO)
      else if(mesh_file%iflag_format .eq. iflag_single) then
        call mpi_input_edge_file                                        &
     &     (nprocs, my_rank, file_name, fem_IO)
!
#ifdef ZLIB_IO
      else if(mesh_file%iflag_format                                    &
     &        .eq. iflag_single+id_gzip_bin_file_fmt) then
        call gz_mpi_input_edge_file_b                                   &
     &     (nprocs, my_rank, file_name, fem_IO)
      else if(mesh_file%iflag_format                                    &
     &        .eq. iflag_single+id_gzip_txt_file_fmt) then
        call gz_mpi_input_edge_file                                     &
     &     (nprocs, my_rank, file_name, fem_IO)
#endif
!
      else
        call sel_read_edge_mesh(mesh_file, my_rank, fem_IO, ierr)
      end if 
!
      if(ierr .gt. 0) then
        call calypso_mpi_abort(ierr, 'Mesh data is wrong!!')
      end if
!
      end subroutine sel_mpi_read_edge_mesh
!
!------------------------------------------------------------------
!------------------------------------------------------------------
!
      subroutine sel_mpi_write_ele_mesh_file(mesh_file, fem_IO)
!
      type(field_IO_params), intent(in) ::  mesh_file
      type(mesh_data), intent(inout) :: fem_IO
!
!
      call set_mesh_file_name_by_param(mesh_file, my_rank, file_name)
!
      if(mesh_file%iflag_format                                         &
     &     .eq. iflag_single+id_binary_file_fmt) then
        call mpi_output_element_file_b                                  &
     &     (nprocs, my_rank, file_name, fem_IO)
      else if(mesh_file%iflag_format .eq. iflag_single) then
        call mpi_output_element_file                                    &
     &     (nprocs, my_rank, file_name, fem_IO)
!
#ifdef ZLIB_IO
      else if(mesh_file%iflag_format                                    &
     &        .eq. iflag_single+id_gzip_bin_file_fmt) then
        call gz_mpi_output_element_file_b                               &
     &     (nprocs, my_rank, file_name, fem_IO)
      else if(mesh_file%iflag_format                                    &
     &        .eq. iflag_single+id_gzip_txt_file_fmt) then
        call gz_mpi_output_element_file                                 &
     &     (nprocs, my_rank, file_name, fem_IO)
#endif
!
      else
        call sel_write_ele_mesh_file(mesh_file, my_rank, fem_IO)
      end if
!
      end subroutine sel_mpi_write_ele_mesh_file
!
!  ---------------------------------------------------------------------
!
      subroutine sel_mpi_write_surf_mesh_file(mesh_file, fem_IO)
!
      type(field_IO_params), intent(in) ::  mesh_file
      type(mesh_data), intent(inout) :: fem_IO
!
!
      call set_mesh_file_name_by_param(mesh_file, my_rank, file_name)
!
      if(mesh_file%iflag_format                                         &
     &     .eq. iflag_single+id_binary_file_fmt) then
        call mpi_output_surface_file_b                                  &
     &     (nprocs, my_rank, file_name, fem_IO)
      else if(mesh_file%iflag_format .eq. iflag_single) then
        call mpi_output_surface_file                                    &
     &     (nprocs, my_rank, file_name, fem_IO)
!
#ifdef ZLIB_IO
      else if(mesh_file%iflag_format                                    &
     &        .eq. iflag_single+id_gzip_bin_file_fmt) then
        call gz_mpi_output_surface_file_b                               &
     &     (nprocs, my_rank, file_name, fem_IO)
      else if(mesh_file%iflag_format                                    &
     &        .eq. iflag_single+id_gzip_txt_file_fmt) then
        call gz_mpi_output_surface_file                                 &
     &     (nprocs, my_rank, file_name, fem_IO)
#endif
!
      else
        call sel_write_surf_mesh_file(mesh_file, my_rank, fem_IO)
      end if
!
      end subroutine sel_mpi_write_surf_mesh_file
!
!  ---------------------------------------------------------------------
!
      subroutine sel_mpi_write_edge_mesh_file(mesh_file, fem_IO)
!
      type(field_IO_params), intent(in) ::  mesh_file
      type(mesh_data), intent(inout) :: fem_IO
!
!
      call set_mesh_file_name_by_param(mesh_file, my_rank, file_name)
!
      if(mesh_file%iflag_format                                         &
     &     .eq. iflag_single+id_binary_file_fmt) then
        call mpi_output_edge_file_b                                     &
     &     (nprocs, my_rank, file_name, fem_IO)
      else if(mesh_file%iflag_format .eq. iflag_single) then
        call mpi_output_edge_file                                       &
     &     (nprocs, my_rank, file_name, fem_IO)
!
#ifdef ZLIB_IO
      else if(mesh_file%iflag_format                                    &
     &        .eq. iflag_single+id_gzip_bin_file_fmt) then
        call gz_mpi_output_edge_file_b                                  &
     &     (nprocs, my_rank, file_name, fem_IO)
      else if(mesh_file%iflag_format                                    &
     &        .eq. iflag_single+id_gzip_txt_file_fmt) then
        call gz_mpi_output_edge_file                                    &
     &     (nprocs, my_rank, file_name, fem_IO)
#endif
!
      else
        call sel_write_edge_mesh_file(mesh_file, my_rank, fem_IO)
      end if
!
      end subroutine sel_mpi_write_edge_mesh_file
!
!  ---------------------------------------------------------------------
!
      end module element_mesh_MPI_IO_select