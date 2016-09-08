!>@file  gz_MPI_groups_IO_b.f90
!!       module gz_MPI_groups_IO_b
!!
!!@author H. Matsui
!!@date   Programmed in July, 2007
!
!> @brief Binary output routines for group data
!!
!!@verbatim
!!      subroutine gz_mpi_read_group_data_b                             &
!!     &         (id_file, nprocs_in, id_rank, ioff_gl, group_IO)
!!      subroutine gz_mpi_read_surf_grp_data_b                          &
!!     &         (id_file, nprocs_in, id_rank, ioff_gl, surf_grp_IO)
!!
!!      subroutine gz_mpi_write_grp_data_b(id_file, ioff_gl, group_IO)
!!      subroutine gz_mpi_write_surf_grp_data_b                         &
!!     &         (id_file, ioff_gl, surf_grp_IO)
!!        type(group_data), intent(inout) :: group_IO
!!        type(surface_group_data), intent(inout) :: surf_grp_IO
!!@endverbatim
!
      module gz_MPI_groups_IO_b
!
      use m_precision
      use m_machine_parameter
!
      use t_group_data
!
      use gz_MPI_binary_data_IO
      use gz_MPI_binary_head_IO
      use gz_MPI_binary_datum_IO
!
      implicit none
!
!------------------------------------------------------------------
!
      contains
!
!------------------------------------------------------------------
!
      subroutine gz_mpi_read_group_data_b                               &
     &         (id_file, nprocs_in, id_rank, ioff_gl, group_IO)
!
      integer, intent(in) ::  id_file
      integer(kind=kint), intent(in) :: id_rank, nprocs_in
      integer(kind = kint_gl), intent(inout) :: ioff_gl
!
      type(group_data), intent(inout) :: group_IO
!
!
      call gz_mpi_read_one_inthead_b                                    &
     &   (id_file, ioff_gl, group_IO%num_grp)
      call allocate_grp_type_num(group_IO)
!
      if (group_IO%num_grp .gt. 0) then
        call gz_mpi_read_mul_charahead_b                                &
     &     (id_file, ioff_gl, group_IO%num_grp, group_IO%grp_name)
        call gz_mpi_read_integer_stack_b                                &
     &     (id_file, nprocs_in, id_rank, ioff_gl,                       &
     &      group_IO%num_grp, group_IO%istack_grp, group_IO%num_item)
!
        call allocate_grp_type_item(group_IO)
!
        call gz_mpi_read_int_vector_b                                   &
     &     (id_file, nprocs_in, id_rank, ioff_gl,                       &
     &      group_IO%num_item, group_IO%item_grp)
      else
        group_IO%num_item = 0
        call allocate_grp_type_item(group_IO)
      end if
!
      end subroutine gz_mpi_read_group_data_b
!
!------------------------------------------------------------------
!
      subroutine gz_mpi_read_surf_grp_data_b                            &
     &         (id_file, nprocs_in, id_rank, ioff_gl, surf_grp_IO)
!
      integer, intent(in) ::  id_file
      integer(kind=kint), intent(in) :: id_rank, nprocs_in
      integer(kind = kint_gl), intent(inout) :: ioff_gl
!
      type(surface_group_data), intent(inout) :: surf_grp_IO
!
      integer(kind = kint) :: nitem
!
!
      call gz_mpi_read_one_inthead_b                                    &
     &   (id_file, ioff_gl, surf_grp_IO%num_grp)
      call allocate_sf_grp_type_num(surf_grp_IO)
!
      if (surf_grp_IO%num_grp .gt. 0) then
        call gz_mpi_read_mul_charahead_b(id_file, ioff_gl,              &
     &      surf_grp_IO%num_grp, surf_grp_IO%grp_name)
        call gz_mpi_read_integer_stack_b                                &
     &     (id_file, nprocs_in, id_rank, ioff_gl,                       &
     &      surf_grp_IO%num_grp, surf_grp_IO%istack_grp,                &
     &      surf_grp_IO%num_item)
!
        call allocate_sf_grp_type_item(surf_grp_IO)
!
        nitem = 2 * surf_grp_IO%num_item
        call gz_mpi_read_int_vector_b                                   &
     &     (id_file, nprocs_in, id_rank, ioff_gl, nitem,                &
     &      surf_grp_IO%item_sf_grp)
      else
        call allocate_sf_grp_type_item(surf_grp_IO)
      end if
!
      end subroutine gz_mpi_read_surf_grp_data_b
!
!------------------------------------------------------------------
!------------------------------------------------------------------
!
      subroutine gz_mpi_write_grp_data_b(id_file, ioff_gl, group_IO)
!
      integer, intent(in) ::  id_file
      integer(kind = kint_gl), intent(inout) :: ioff_gl
!
      type(group_data), intent(inout) :: group_IO
!
!
      call gz_mpi_write_one_inthead_b                                   &
     &   (id_file, ioff_gl, group_IO%num_grp)
!
      if (group_IO%num_grp .gt. 0) then
        call gz_mpi_write_mul_charahead_b                               &
     &     (id_file, ioff_gl, group_IO%num_grp, group_IO%grp_name)
        call gz_mpi_write_integer_stack_b                               &
     &     (id_file, ioff_gl, group_IO%num_grp, group_IO%istack_grp)
        call gz_mpi_write_int_vector_b                                  &
     &     (id_file, ioff_gl, group_IO%num_item, group_IO%item_grp)
      end if
!
      call deallocate_grp_type(group_IO)
!
      end subroutine gz_mpi_write_grp_data_b
!
!------------------------------------------------------------------
!
      subroutine gz_mpi_write_surf_grp_data_b                           &
     &         (id_file, ioff_gl, surf_grp_IO)
!
      integer, intent(in) ::  id_file
      integer(kind = kint_gl), intent(inout) :: ioff_gl
!
      type(surface_group_data), intent(inout) :: surf_grp_IO
!
      integer(kind = kint) :: nitem
!
!
      call gz_mpi_write_one_inthead_b                                   &
     &   (id_file, ioff_gl, surf_grp_IO%num_grp)
!
      if (surf_grp_IO%num_grp .gt. 0) then
        call gz_mpi_write_mul_charahead_b(id_file, ioff_gl,             &
     &     surf_grp_IO%num_grp, surf_grp_IO%grp_name)
        call gz_mpi_write_integer_stack_b(id_file, ioff_gl,             &
     &      surf_grp_IO%num_grp, surf_grp_IO%istack_grp)
!
        nitem = 2 * surf_grp_IO%num_item
        call gz_mpi_write_int_vector_b                                  &
     &     (id_file, ioff_gl, nitem, surf_grp_IO%item_sf_grp)
      end if
!
      call deallocate_sf_grp_type(surf_grp_IO)
!
      end subroutine gz_mpi_write_surf_grp_data_b
!
!------------------------------------------------------------------
!
      end module gz_MPI_groups_IO_b
