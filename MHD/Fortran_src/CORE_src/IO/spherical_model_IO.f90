!>@file  spherical_model_IO.f90
!!       module spherical_model_IO
!!
!!@author H. Matsui
!!@date        programmed by H.Matsui in July, 2007
!
!> @brief  Data IO routines for spectrum data
!!
!!@verbatim
!!      subroutine read_rank_4_sph(id_file)
!!      subroutine read_gl_resolution_sph(id_file)
!!      subroutine read_gl_nodes_sph(id_file)
!!
!!      subroutine write_rank_4_sph(id_file)
!!      subroutine write_gl_resolution_sph(id_file)
!!      subroutine write_gl_nodes_sph(id_file)
!!@endverbatim
!
      module spherical_model_IO
!
      use m_precision
!
      use m_node_id_spherical_IO
!
      implicit none
!
      character(len=255) :: character_4_read = ''
      private :: character_4_read
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine read_rank_4_sph(id_file)
!
      use skip_comment_f
!
      integer(kind = kint), intent(in) :: id_file
!
      call skip_comment(character_4_read,id_file)
      read(character_4_read,*) sph_IO1%sph_rank(1:ndir_sph_IO)
!
      end subroutine read_rank_4_sph
!
! -----------------------------------------------------------------------
!
      subroutine read_gl_resolution_sph(id_file)
!
      use skip_comment_f
!
      integer(kind = kint), intent(in) :: id_file
!
      call skip_comment(character_4_read,id_file)
      read(character_4_read,*) nidx_gl_sph_IO(1:ndir_sph_IO)
      call skip_comment(character_4_read,id_file)
      read(character_4_read,*) ltr_gl_IO
!
      end subroutine read_gl_resolution_sph
!
! -----------------------------------------------------------------------
!
      subroutine read_gl_nodes_sph(id_file)
!
      use skip_comment_f
!
      integer(kind = kint), intent(in) :: id_file
      integer(kind = kint) :: i
!
!
      call skip_comment(character_4_read,id_file)
      read(character_4_read,*) nnod_sph_IO
!
      call allocate_nod_id_sph_IO
!
      call skip_comment(character_4_read,id_file)
      read(character_4_read,*)                                          &
     &           inod_gl_sph_IO(1), idx_gl_sph_IO(1,1:ndir_sph_IO)
      do i = 2, nnod_sph_IO
        read(id_file,*)                                                 &
     &           inod_gl_sph_IO(i), idx_gl_sph_IO(i,1:ndir_sph_IO)
      end do
!
      end subroutine read_gl_nodes_sph
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine write_rank_4_sph(id_file)
!
      use m_sph_modes_grid_labels
!
      integer(kind = kint), intent(in) :: id_file
!
!
      write(id_file,'(a)', advance='NO') hd_segment()
      write(id_file,'(10i16)') sph_IO1%sph_rank(1:ndir_sph_IO)
!
      end subroutine write_rank_4_sph
!
! -----------------------------------------------------------------------
!
      subroutine write_gl_resolution_sph(id_file)
!
      use m_sph_modes_grid_labels
!
      integer(kind = kint), intent(in) :: id_file
!
!
      write(id_file,'(a)', advance='NO') hd_trunc()
      write(id_file,'(3i16)') nidx_gl_sph_IO(1:ndir_sph_IO)
      write(id_file,'(i16)') ltr_gl_IO
!
      end subroutine write_gl_resolution_sph
!
! -----------------------------------------------------------------------
!
      subroutine write_gl_nodes_sph(id_file)
!
      integer(kind = kint), intent(in) :: id_file
      integer(kind = kint) :: i
!
!
      write(id_file,'(i16)') nnod_sph_IO
      do i = 1, nnod_sph_IO
        write(id_file,'(20i16)')                                        &
     &              inod_gl_sph_IO(i), idx_gl_sph_IO(i,1:ndir_sph_IO)
      end do
!
      call deallocate_nod_id_sph_IO
!
      end subroutine write_gl_nodes_sph
!
! -----------------------------------------------------------------------
!
      end module spherical_model_IO
