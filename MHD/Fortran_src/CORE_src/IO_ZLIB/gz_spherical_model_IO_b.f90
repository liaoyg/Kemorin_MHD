!>@file  gz_spherical_model_IO_b.f90
!!       module gz_spherical_model_IO_b
!!
!!@author H. Matsui
!!@date        programmed by H.Matsui in July, 2007
!
!> @brief  Data IO routines for spectrum data
!!
!!@verbatim
!!      subroutine gz_read_rank_4_sph_b
!!      subroutine gz_read_gl_resolution_sph_b
!!      subroutine gz_read_gl_nodes_sph_b
!!
!!      subroutine gz_write_rank_4_sph_b
!!      subroutine gz_write_gl_resolution_sph_b
!!      subroutine gz_write_gl_nodes_sph_b
!!@endverbatim
!
      module gz_spherical_model_IO_b
!
      use m_precision
!
      use m_node_id_spherical_IO
      use gz_binary_IO
!
      implicit none
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine gz_read_rank_4_sph_b
!
!
      call gz_read_fld_mul_inthead_b(ndir_sph_IO, sph_rank_IO)
!
      end subroutine gz_read_rank_4_sph_b
!
! -----------------------------------------------------------------------
!
      subroutine gz_read_gl_resolution_sph_b
!
!
      call gz_read_fld_mul_inthead_b(ndir_sph_IO, nidx_gl_sph_IO)
      call gz_read_fld_inthead_b(ltr_gl_IO)
!
      end subroutine gz_read_gl_resolution_sph_b
!
! -----------------------------------------------------------------------
!
      subroutine gz_read_gl_nodes_sph_b
!
      integer(kind = kint) :: nvect
!
!
      call gz_read_fld_inthead_b(nnod_sph_IO)
      call allocate_nod_id_sph_IO
!
      call gz_read_fld_mul_i8head_b(nnod_sph_IO, inod_gl_sph_IO)
      nvect = nnod_sph_IO * ndir_sph_IO
      call gz_read_fld_mul_inthead_b(nvect, idx_gl_sph_IO)
!
      end subroutine gz_read_gl_nodes_sph_b
!
! -----------------------------------------------------------------------
!
      subroutine gz_write_rank_4_sph_b
!
!
      call gz_write_fld_mul_inthead_b(ndir_sph_IO, sph_rank_IO)
!
      end subroutine gz_write_rank_4_sph_b
!
! -----------------------------------------------------------------------
!
      subroutine gz_write_gl_resolution_sph_b
!
!
      call gz_write_fld_mul_inthead_b(ndir_sph_IO, nidx_gl_sph_IO)
      call gz_write_fld_inthead_b(ltr_gl_IO)
!
      end subroutine gz_write_gl_resolution_sph_b
!
! -----------------------------------------------------------------------
!
      subroutine gz_write_gl_nodes_sph_b
!
      integer(kind = kint) ::  nvect
!
!
      call gz_write_fld_inthead_b(nnod_sph_IO)
      call gz_write_fld_mul_i8head_b(nnod_sph_IO, inod_gl_sph_IO)
      nvect = nnod_sph_IO * ndir_sph_IO
      call gz_write_fld_mul_inthead_b(nvect, idx_gl_sph_IO)
!
      call deallocate_nod_id_sph_IO
!
      end subroutine gz_write_gl_nodes_sph_b
!
! -----------------------------------------------------------------------
!
      end module gz_spherical_model_IO_b
