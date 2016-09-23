!>@file   MPI_sph_gl_1d_idx_IO.f90
!!@brief  module MPI_sph_gl_1d_idx_IO
!!
!!@author H.Matsui
!!@date      Programmed in Aug., 2016
!
!>@brief  Mesh file IO for gxipped format
!!
!!@verbatim
!!      subroutine mpi_read_rtp_gl_1d_table(IO_param, sph_IO)
!!      subroutine mpi_read_rj_gl_1d_table(IO_param, sph_IO)
!!
!!      subroutine mpi_write_rtp_gl_1d_table(IO_param, sph_IO)
!!      subroutine mpi_write_rj_gl_1d_table(IO_param, sph_IO)
!!        type(calypso_MPI_IO_params), intent(inout) :: IO_param
!!        type(sph_IO_data), intent(inout) :: sph_IO
!!@endverbatim
!
      module MPI_sph_gl_1d_idx_IO
!
      use m_precision
      use m_constants
!
      use t_node_id_spherical_IO
      use t_calypso_mpi_IO_param
      use m_sph_modes_grid_labels
      use MPI_ascii_data_IO
      use MPI_domain_data_IO
      use MPI_ascii_data_IO
!
      implicit none
!
      integer(kind = kint_gl), allocatable :: idx_gl_tmp(:)
      private :: idx_gl_tmp
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------!
      subroutine mpi_read_rtp_gl_1d_table(IO_param, sph_IO)
!
      type(calypso_MPI_IO_params), intent(inout) :: IO_param
      type(sph_IO_data), intent(inout) :: sph_IO
!
!
      sph_IO%numdir_sph = 3
      sph_IO%ncomp_table_1d(1) = 1
      sph_IO%ncomp_table_1d(2) = 1
      sph_IO%ncomp_table_1d(3) = 2
!
      call alloc_num_idx_sph_IO(sph_IO)
!
      call mpi_skip_read(IO_param, len(hd_rgrid()))
      call mpi_read_num_of_data(IO_param, sph_IO%ist_sph(1))
      call mpi_read_num_of_data(IO_param, sph_IO%ied_sph(1))
      call mpi_read_num_of_data(IO_param, sph_IO%nidx_sph(1))
!
      call alloc_idx_sph_1d1_IO(sph_IO)
!
      allocate(idx_gl_tmp(sph_IO%nidx_sph(1)))
      call mpi_read_node_position(IO_param,                             &
     &   sph_IO%nidx_sph(1), sph_IO%ncomp_table_1d(1),                  &
     &   idx_gl_tmp, sph_IO%r_gl_1)
      sph_IO%idx_gl_1(1:sph_IO%nidx_sph(1))                             &
     &       = int(idx_gl_tmp(1:sph_IO%nidx_sph(1)))
      deallocate(idx_gl_tmp)
!
!
      call mpi_skip_read(IO_param, len(hd_tgrid()))
      call mpi_read_num_of_data(IO_param, sph_IO%ist_sph(2))
      call mpi_read_num_of_data(IO_param, sph_IO%ied_sph(2))
      call mpi_read_num_of_data(IO_param, sph_IO%nidx_sph(2))
!
      call alloc_idx_sph_1d2_IO(sph_IO)
!
      call mpi_read_1d_gl_address(IO_param,                             &
     &   sph_IO%nidx_sph(2), sph_IO%ncomp_table_1d(2), sph_IO%idx_gl_2)
!
!
      call mpi_skip_read(IO_param, len(hd_pgrid()))
      call mpi_read_num_of_data(IO_param, sph_IO%ist_sph(3))
      call mpi_read_num_of_data(IO_param, sph_IO%ied_sph(3))
      call mpi_read_num_of_data(IO_param, sph_IO%nidx_sph(3))
!
      call alloc_idx_sph_1d3_IO(sph_IO)
!
      call mpi_read_1d_gl_address(IO_param,                             &
     &   sph_IO%nidx_sph(3), sph_IO%ncomp_table_1d(3), sph_IO%idx_gl_3)
!
      end subroutine mpi_read_rtp_gl_1d_table
!
! -----------------------------------------------------------------------
!
      subroutine mpi_read_rj_gl_1d_table(IO_param, sph_IO)
!
      type(calypso_MPI_IO_params), intent(inout) :: IO_param
      type(sph_IO_data), intent(inout) :: sph_IO
!
!
      sph_IO%numdir_sph = 2
      sph_IO%ncomp_table_1d(1) = 1
      sph_IO%ncomp_table_1d(2) = 3
!
      call alloc_num_idx_sph_IO(sph_IO)
!
      call mpi_skip_read(IO_param, len(hd_rgrid()))
      call mpi_read_num_of_data(IO_param, sph_IO%ist_sph(1))
      call mpi_read_num_of_data(IO_param, sph_IO%ied_sph(1))
      call mpi_read_num_of_data(IO_param, sph_IO%nidx_sph(1))
!
      call alloc_idx_sph_1d1_IO(sph_IO)
!
      allocate(idx_gl_tmp(sph_IO%nidx_sph(1)))
      call mpi_read_node_position(IO_param,                             &
     &   sph_IO%nidx_sph(1), sph_IO%ncomp_table_1d(1),                  &
     &   idx_gl_tmp, sph_IO%r_gl_1)
      sph_IO%idx_gl_1(1:sph_IO%nidx_sph(1))                             &
     &       = int(idx_gl_tmp(1:sph_IO%nidx_sph(1)))
      deallocate(idx_gl_tmp)
!
!
      call mpi_skip_read(IO_param, len(hd_jmode()))
      call mpi_read_num_of_data(IO_param, sph_IO%ist_sph(2))
      call mpi_read_num_of_data(IO_param, sph_IO%ied_sph(2))
      call mpi_read_num_of_data(IO_param, sph_IO%nidx_sph(2))
!
      call alloc_idx_sph_1d2_IO(sph_IO)
!
      call mpi_read_1d_gl_address(IO_param,                             &
     &   sph_IO%nidx_sph(2), sph_IO%ncomp_table_1d(2), sph_IO%idx_gl_2)
!
!
      end subroutine mpi_read_rj_gl_1d_table
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine mpi_write_rtp_gl_1d_table(IO_param, sph_IO)
!
      type(calypso_MPI_IO_params), intent(inout) :: IO_param
      type(sph_IO_data), intent(inout) :: sph_IO
!
!
      call mpi_write_charahead                                          &
     &   (IO_param, len(hd_rgrid()), hd_rgrid())
!
      call mpi_write_num_of_data(IO_param, sph_IO%ist_sph(1))
      call mpi_write_num_of_data(IO_param, sph_IO%ied_sph(1))
!
      allocate(idx_gl_tmp(sph_IO%nidx_sph(1)))
      idx_gl_tmp(1:sph_IO%nidx_sph(1))                                  &
     &       =  sph_IO%idx_gl_1(1:sph_IO%nidx_sph(1))
      call mpi_write_node_position(IO_param,                            &
     &   sph_IO%nidx_sph(1), ione, idx_gl_tmp, sph_IO%r_gl_1)
      deallocate(idx_gl_tmp)
!
!
      call mpi_write_charahead                                          &
     &   (IO_param, len(hd_tgrid()), hd_tgrid())
!
      call mpi_write_num_of_data(IO_param, sph_IO%ist_sph(2))
      call mpi_write_num_of_data(IO_param, sph_IO%ied_sph(2))
!
      call mpi_write_1d_gl_address(IO_param,                            &
     &    sph_IO%nidx_sph(2), sph_IO%ncomp_table_1d(2),                 &
     &    sph_IO%idx_gl_2)
!
!
      call mpi_write_charahead                                          &
     &   (IO_param, len(hd_pgrid()), hd_pgrid())
!
      call mpi_write_num_of_data(IO_param, sph_IO%ist_sph(3))
      call mpi_write_num_of_data(IO_param, sph_IO%ied_sph(3))
!
      call mpi_write_1d_gl_address(IO_param,                            &
     &    sph_IO%nidx_sph(3), sph_IO%ncomp_table_1d(3),                 &
     &    sph_IO%idx_gl_3)
!
!
      call dealloc_num_idx_sph_IO(sph_IO)
      call dealloc_idx_sph_1d1_IO(sph_IO)
      call dealloc_idx_sph_1d2_IO(sph_IO)
      call dealloc_idx_sph_1d3_IO(sph_IO)
!
      end subroutine mpi_write_rtp_gl_1d_table
!
! ----------------------------------------------------------------------
!
      subroutine mpi_write_rj_gl_1d_table(IO_param, sph_IO)
!
      type(calypso_MPI_IO_params), intent(inout) :: IO_param
      type(sph_IO_data), intent(inout) :: sph_IO
!
!
      call mpi_write_charahead                                          &
     &   (IO_param, len(hd_rgrid()), hd_rgrid())
!
      call mpi_write_num_of_data(IO_param, sph_IO%ist_sph(1))
      call mpi_write_num_of_data(IO_param, sph_IO%ied_sph(1))
!
      allocate(idx_gl_tmp(sph_IO%nidx_sph(1)))
      idx_gl_tmp(1:sph_IO%nidx_sph(1))                                  &
     &       =  sph_IO%idx_gl_1(1:sph_IO%nidx_sph(1))
      call mpi_write_node_position(IO_param,                            &
     &    sph_IO%nidx_sph(1), sph_IO%ncomp_table_1d(1),                 &
     &    idx_gl_tmp, sph_IO%r_gl_1)
      deallocate(idx_gl_tmp)
!
!
      call mpi_write_charahead                                          &
     &   (IO_param, len(hd_jmode()), hd_jmode())
!
      call mpi_write_num_of_data(IO_param, sph_IO%ist_sph(2))
      call mpi_write_num_of_data(IO_param, sph_IO%ied_sph(2))
!
      call mpi_write_1d_gl_address(IO_param,                            &
     &    sph_IO%nidx_sph(2), sph_IO%ncomp_table_1d(2),                 &
     &    sph_IO%idx_gl_2)
!
      call dealloc_num_idx_sph_IO(sph_IO)
      call dealloc_idx_sph_1d1_IO(sph_IO)
      call dealloc_idx_sph_1d2_IO(sph_IO)
!
      end subroutine mpi_write_rj_gl_1d_table
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine mpi_read_node_position                                 &
     &         (IO_param, nnod, numdir, id_global, xx)
!
      type(calypso_MPI_IO_params), intent(inout) :: IO_param
      integer(kind=kint), intent(in) :: nnod, numdir
      integer(kind=kint_gl), intent(inout) :: id_global(nnod)
      real(kind=kreal), intent(inout) :: xx(nnod, numdir)
!
      real(kind = kreal) :: xx_tmp(numdir)
!
      integer(kind = kint) :: i, led, n_item, ilength
      integer(kind = MPI_OFFSET_KIND) :: ioffset
!
!
      call mpi_skip_read                                                &
     &   (IO_param, len_multi_int_textline(IO_param%nprocs_in))
!
      IO_param%istack_merged(0) = 0
      do i = 1, IO_param%nprocs_in
        n_item = int(IO_param%istack_merged(i))
        if(n_item .le. 0) then
          led = ione
        else if(n_item .gt. 0) then
          led = len_int8_and_vector_textline(numdir) * n_item
        end if
        IO_param%istack_merged(i) = IO_param%istack_merged(i-1)         &
     &                             + led
      end do
!
      if(nnod .eq. 0) then
        led = ione
      else
        ioffset = IO_param%ioff_gl                                      &
     &           + IO_param%istack_merged(IO_param%id_rank)
!
        ilength = len_int8_and_vector_textline(numdir)
        do i = 1, nnod
          call read_int8_and_vector_textline                            &
     &       (calypso_mpi_seek_read_chara(IO_param%id_file,             &
     &                                    ioffset, ilength),            &
     &        id_global(i), numdir, xx_tmp)
          xx(i,1:numdir) = xx_tmp(1:numdir)
        end do
      end if
      IO_param%ioff_gl = IO_param%ioff_gl                               &
     &         + IO_param%istack_merged(IO_param%nprocs_in)
!
      end subroutine mpi_read_node_position
!
! -----------------------------------------------------------------------
!
      subroutine mpi_read_1d_gl_address                                 &
     &         (IO_param, nnod, numdir, idx)
!
      type(calypso_MPI_IO_params), intent(inout) :: IO_param
      integer(kind=kint), intent(in) :: nnod, numdir
      integer(kind=kint), intent(inout) :: idx(nnod, numdir)
!
      integer(kind = kint) :: idx_tmp(numdir)

      integer(kind = kint) :: i, led, ilength, n_item
      integer(kind = MPI_OFFSET_KIND) :: ioffset
!
!
      call mpi_skip_read                                                &
     &   (IO_param, len_multi_int_textline(IO_param%nprocs_in))
!
      IO_param%istack_merged(0) = 0
      do i = 1, IO_param%nprocs_in
        n_item = int(IO_param%istack_merged(i))
        if(n_item .le. 0) then
          led = ione
        else if(n_item .gt. 0) then
          led = len_multi_int_textline(numdir) * n_item
        end if
        IO_param%istack_merged(i) = IO_param%istack_merged(i-1)         &
     &                             + led
      end do
!
      if(nnod .eq. 0) then
        led = ione
      else
        ioffset = IO_param%ioff_gl                                      &
     &           + IO_param%istack_merged(IO_param%id_rank)
!
        ilength = len_multi_int_textline(numdir)
        do i = 1, nnod
          call read_multi_int_textline                                  &
     &       (calypso_mpi_seek_read_chara(IO_param%id_file,             &
     &                                    ioffset, ilength),            &
     &        numdir, idx_tmp)
          idx(i,1:numdir) = idx_tmp(1:numdir)
        end do
      end if
      IO_param%ioff_gl = IO_param%ioff_gl                               &
     &         + IO_param%istack_merged(IO_param%nprocs_in)
!
      end subroutine mpi_read_1d_gl_address
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine mpi_write_node_position                                &
     &         (IO_param, nnod, numdir, id_global, xx)
!
      type(calypso_MPI_IO_params), intent(inout) :: IO_param
      integer(kind=kint), intent(in) :: nnod, numdir
      integer(kind=kint_gl), intent(in) :: id_global(nnod)
      real(kind=kreal), intent(in) :: xx(nnod, numdir)
!
      integer(kind = kint) :: i, led, ilength
      real(kind = kreal) :: xx_tmp(numdir)
      integer(kind = MPI_OFFSET_KIND) :: ioffset
!
!
      call mpi_write_num_of_data(IO_param, nnod)
!
      ilength = len_int8_and_vector_textline(numdir)
      led = nnod * len_int8_and_vector_textline(numdir)
      call mpi_write_stack_over_domain(IO_param, led)
!
      ioffset = IO_param%ioff_gl                                        &
     &         + IO_param%istack_merged(IO_param%id_rank)
      IO_param%ioff_gl = IO_param%ioff_gl                               &
     &         + IO_param%istack_merged(IO_param%nprocs_in)
!
      if(IO_param%id_rank .ge. IO_param%nprocs_in) return
      if(nnod .le. 0) then
        call calypso_mpi_seek_write_chara                               &
     &     (IO_param%id_file, ioffset, ione, char(10))
      else
        do i = 1, nnod
          xx_tmp(1:numdir) = xx(i,1:numdir)
          call calypso_mpi_seek_write_chara                             &
     &       (IO_param%id_file, ioffset, ilength,                       &
     &        int8_and_vector_textline(id_global(i), numdir, xx_tmp))
        end do
      end if
      call calypso_mpi_barrier
!
      end subroutine mpi_write_node_position
!
! -----------------------------------------------------------------------
!
      subroutine mpi_write_1d_gl_address                                &
     &         (IO_param, nnod, numdir, idx)
!
      type(calypso_MPI_IO_params), intent(inout) :: IO_param
      integer(kind=kint), intent(in) :: nnod, numdir
      integer(kind=kint), intent(in) :: idx(nnod, numdir)
!
      integer(kind = kint) :: i, led, ilength
      integer(kind = kint) :: idx_tmp(numdir)
      integer(kind = MPI_OFFSET_KIND) :: ioffset
!
!
      call mpi_write_num_of_data(IO_param, nnod)
!
      ilength = len_multi_int_textline(numdir)
      if(nnod .le. 0) then
        led = ione
      else if(nnod .gt. 0) then
        led = nnod * len_multi_int_textline(numdir)
      end if
!
      call mpi_write_stack_over_domain(IO_param, led)
!
      ioffset = IO_param%ioff_gl                                        &
     &         + IO_param%istack_merged(IO_param%id_rank)
      IO_param%ioff_gl = IO_param%ioff_gl                               &
     &         + IO_param%istack_merged(IO_param%nprocs_in)
!
      if(nnod .le. 0) then
        call calypso_mpi_seek_write_chara                               &
     &     (IO_param%id_file, ioffset, ione, char(10))
      else
        do i = 1, nnod
          idx_tmp(1:numdir) = idx(i,1:numdir)
          call calypso_mpi_seek_write_chara                             &
     &       (IO_param%id_file, ioffset, ilength,                       &
     &        multi_int_textline(numdir, idx_tmp))
        end do
      end if
      call calypso_mpi_barrier
!
      end subroutine mpi_write_1d_gl_address
!
! -----------------------------------------------------------------------
!
      end module MPI_sph_gl_1d_idx_IO
