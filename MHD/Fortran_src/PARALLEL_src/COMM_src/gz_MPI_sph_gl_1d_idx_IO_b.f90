!>@file   gz_MPI_sph_gl_1d_idx_IO_b.f90
!!@brief  module gz_MPI_sph_gl_1d_idx_IO_b
!!
!!@author H.Matsui
!!@date      Programmed in Aug., 2016
!
!>@brief  Mesh file IO for gxipped format
!!
!!@verbatim
!!      subroutine gz_mpi_read_rtp_gl_1d_table_b                        &
!!     &         (id_file, nprocs_in, id_rank, ioff_gl)
!!      subroutine gz_mpi_read_rj_gl_1d_table_b                         &
!!     &         (id_file, nprocs_in, id_rank, ioff_gl)
!!
!!      subroutine gz_mpi_write_rtp_gl_1d_table_b(id_file, ioff_gl)
!!      subroutine gz_mpi_write_rj_gl_1d_table_b(id_file, ioff_gl)
!!@endverbatim
!
      module gz_MPI_sph_gl_1d_idx_IO_b
!
      use m_precision
      use m_constants
!
      use m_node_id_spherical_IO
      use gz_MPI_binary_data_IO
!
      implicit none
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------!
      subroutine gz_mpi_read_rtp_gl_1d_table_b                          &
     &         (id_file, nprocs_in, id_rank, ioff_gl)
!
      integer, intent(in) ::  id_file
      integer(kind = kint_gl), intent(inout) :: ioff_gl
      integer(kind=kint), intent(in) :: id_rank, nprocs_in
!
      integer(kind = kint) :: nvect
!
!
      sph_IO1%numdir_sph = 3
      sph_IO1%ncomp_table_1d(1) = 1
      sph_IO1%ncomp_table_1d(2) = 1
      sph_IO1%ncomp_table_1d(3) = 2
!
      call gz_mpi_read_int_vector_b                                     &
     &   (id_file, nprocs_in, id_rank, ioff_gl,                         &
     &    sph_IO1%numdir_sph, sph_IO1%nidx_sph)
      call gz_mpi_read_int_vector_b                                     &
     &   (id_file, nprocs_in, id_rank, ioff_gl,                         &
     &    sph_IO1%numdir_sph, sph_IO1%ist_sph)
      call gz_mpi_read_int_vector_b                                     &
     &   (id_file, nprocs_in, id_rank, ioff_gl,                         &
     &    sph_IO1%numdir_sph, sph_IO1%ied_sph)
!
      call alloc_idx_sph_1d1_IO(sph_IO1)
      call alloc_idx_sph_1d2_IO(sph_IO1)
      call alloc_idx_sph_1d3_IO(sph_IO1)
!
      call gz_mpi_read_int_vector_b                                     &
     &   (id_file, nprocs_in, id_rank, ioff_gl,                         &
     &    sph_IO1%nidx_sph(1), sph_IO1%idx_gl_1)
      call gz_mpi_read_1d_vector_b                                      &
     &   (id_file, nprocs_in, id_rank, ioff_gl,                         &
     &    sph_IO1%nidx_sph(1), sph_IO1%r_gl_1)
!
      nvect = sph_IO1%nidx_sph(2) * sph_IO1%ncomp_table_1d(2)
      call gz_mpi_read_int_vector_b                                     &
     &   (id_file, nprocs_in, id_rank, ioff_gl, nvect, sph_IO1%idx_gl_2)
!
      nvect = sph_IO1%nidx_sph(3) * sph_IO1%ncomp_table_1d(3)
      call gz_mpi_read_int_vector_b                                     &
     &   (id_file, nprocs_in, id_rank, ioff_gl, nvect, sph_IO1%idx_gl_3)
!
      end subroutine gz_mpi_read_rtp_gl_1d_table_b
!
! -----------------------------------------------------------------------
!
      subroutine gz_mpi_read_rj_gl_1d_table_b                           &
     &         (id_file, nprocs_in, id_rank, ioff_gl)
!
      integer, intent(in) ::  id_file
      integer(kind = kint_gl), intent(inout) :: ioff_gl
      integer(kind=kint), intent(in) :: id_rank, nprocs_in
!
      integer(kind = kint) :: nvect
!
!
      sph_IO1%numdir_sph = 2
      sph_IO1%ncomp_table_1d(1) = 1
      sph_IO1%ncomp_table_1d(2) = 3
!
      call gz_mpi_read_int_vector_b                                     &
     &   (id_file, nprocs_in, id_rank, ioff_gl,                         &
     &    sph_IO1%numdir_sph, sph_IO1%nidx_sph)
      call gz_mpi_read_int_vector_b                                     &
     &   (id_file, nprocs_in, id_rank, ioff_gl,                         &
     &    sph_IO1%numdir_sph, sph_IO1%ist_sph)
      call gz_mpi_read_int_vector_b                                     &
     &   (id_file, nprocs_in, id_rank, ioff_gl,                         &
     &    sph_IO1%numdir_sph, sph_IO1%ied_sph)
!
      call alloc_idx_sph_1d1_IO(sph_IO1)
      call alloc_idx_sph_1d2_IO(sph_IO1)
!
      call gz_mpi_read_int_vector_b                                     &
     &   (id_file, nprocs_in, id_rank, ioff_gl,                         &
     &    sph_IO1%nidx_sph(1), sph_IO1%idx_gl_1)
      call gz_mpi_read_1d_vector_b                                      &
     &   (id_file, nprocs_in, id_rank, ioff_gl,                         &
     &    sph_IO1%nidx_sph(1), sph_IO1%r_gl_1)
!
      nvect = sph_IO1%nidx_sph(2) * sph_IO1%ncomp_table_1d(2)
      call gz_mpi_read_int_vector_b                                     &
     &   (id_file, nprocs_in, id_rank, ioff_gl, nvect, sph_IO1%idx_gl_2)
!
      end subroutine gz_mpi_read_rj_gl_1d_table_b
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine gz_mpi_write_rtp_gl_1d_table_b(id_file, ioff_gl)
!
      integer, intent(in) ::  id_file
      integer(kind = kint_gl), intent(inout) :: ioff_gl
!
      integer(kind = kint) :: nvect
!
!
      call gz_mpi_write_int_vector_b                                    &
     &   (id_file, ioff_gl, sph_IO1%numdir_sph, sph_IO1%nidx_sph)
      call gz_mpi_write_int_vector_b                                    &
     &   (id_file, ioff_gl, sph_IO1%numdir_sph, sph_IO1%ist_sph)
      call gz_mpi_write_int_vector_b                                    &
     &   (id_file, ioff_gl, sph_IO1%numdir_sph, sph_IO1%ied_sph)
!
      call gz_mpi_write_int_vector_b                                    &
     &   (id_file, ioff_gl, sph_IO1%nidx_sph(1), sph_IO1%idx_gl_1)
      call gz_mpi_write_1d_vector_b                                     &
     &   (id_file, ioff_gl, sph_IO1%nidx_sph(1), sph_IO1%r_gl_1)
!
      nvect = sph_IO1%nidx_sph(2) * sph_IO1%ncomp_table_1d(2)
      call gz_mpi_write_int_vector_b                                    &
     &   (id_file, ioff_gl, nvect, sph_IO1%idx_gl_2)
!
      nvect = sph_IO1%nidx_sph(3) * sph_IO1%ncomp_table_1d(3)
      call gz_mpi_write_int_vector_b                                    &
     &   (id_file, ioff_gl, nvect, sph_IO1%idx_gl_3)
!
      call dealloc_idx_sph_1d1_IO(sph_IO1)
      call dealloc_idx_sph_1d2_IO(sph_IO1)
      call dealloc_idx_sph_1d3_IO(sph_IO1)
!
      end subroutine gz_mpi_write_rtp_gl_1d_table_b
!
! ----------------------------------------------------------------------
!
      subroutine gz_mpi_write_rj_gl_1d_table_b(id_file, ioff_gl)
!
      integer, intent(in) ::  id_file
      integer(kind = kint_gl), intent(inout) :: ioff_gl
!
      integer(kind = kint) :: nvect
!
!
      call gz_mpi_write_int_vector_b                                    &
     &   (id_file, ioff_gl, sph_IO1%numdir_sph, sph_IO1%nidx_sph)
      call gz_mpi_write_int_vector_b                                    &
     &   (id_file, ioff_gl, sph_IO1%numdir_sph, sph_IO1%ist_sph)
      call gz_mpi_write_int_vector_b                                    &
     &   (id_file, ioff_gl, sph_IO1%numdir_sph, sph_IO1%ied_sph)
!
      call gz_mpi_write_int_vector_b                                    &
     &   (id_file, ioff_gl, sph_IO1%nidx_sph(1), sph_IO1%idx_gl_1)
      call gz_mpi_write_1d_vector_b                                     &
     &   (id_file, ioff_gl, sph_IO1%nidx_sph(1), sph_IO1%r_gl_1)
!
      nvect = sph_IO1%nidx_sph(2) * sph_IO1%ncomp_table_1d(2)
      call gz_mpi_write_int_vector_b                                    &
     &   (id_file, ioff_gl, nvect, sph_IO1%idx_gl_2)
!
      call dealloc_idx_sph_1d1_IO(sph_IO1)
      call dealloc_idx_sph_1d2_IO(sph_IO1)
!
      end subroutine gz_mpi_write_rj_gl_1d_table_b
!
! ----------------------------------------------------------------------
!
      end module gz_MPI_sph_gl_1d_idx_IO_b
