!>@file  gz_field_file_IO_b.f90
!!       module gz_field_file_IO_b
!!
!!@author H. Matsui
!!@date   Programmed in May, 2015
!
!> @brief Output merged binary field file using MPI-IO
!!
!!@verbatim
!!      subroutine gz_write_step_fld_file_b(gzip_name, my_rank, fld_IO)
!!
!!      subroutine gz_read_step_field_file_b(gzip_name, my_rank, fld_IO)
!!      subroutine gz_rd_alloc_st_fld_file_b(gzip_name, my_rank, fld_IO)
!!
!!      subroutine gz_rd_alloc_st_fld_head_b(gzip_name, my_rank, fld_IO)
!!@endverbatim
!
      module gz_field_file_IO_b
!
      use m_precision
      use m_constants
      use m_machine_parameter
!
      use m_time_data_IO
      use t_field_data_IO
!
      implicit none
!
      integer(kind = kint), private :: iflag_endian_swap
!
      private :: gz_read_step_data_mpi_b, gz_read_field_data_mpi_b
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine gz_write_step_fld_file_b(gzip_name, my_rank, fld_IO)
!
      character(len=kchara), intent(in) :: gzip_name
      integer(kind = kint), intent(in) :: my_rank
!
      type(field_IO), intent(in) :: fld_IO
!
!
      if(my_rank.eq.0 .or. i_debug .gt. 0) write(*,*)                   &
     &   'Write gzipped binary data file: ', trim(gzip_name)
!
      call open_wt_gzfile(gzip_name)
!
      call gz_write_field_data_b(my_rank,                               &
     &    fld_IO%nnod_IO, fld_IO%num_field_IO, fld_IO%ntot_comp_IO,     &
     &    fld_IO%num_comp_IO, fld_IO%fld_name, fld_IO%d_IO)
!
      call close_gzfile()
!
      end subroutine gz_write_step_fld_file_b
!
! -----------------------------------------------------------------------
!
      subroutine gz_read_step_field_file_b(gzip_name, my_rank, fld_IO)
!
      character(len=kchara), intent(in) :: gzip_name
      integer(kind = kint), intent(in) :: my_rank
!
      type(field_IO), intent(inout) :: fld_IO
!
      integer(kind = kint_gl) :: istack_merged(1)
!
!
      if(my_rank.eq.0 .or. i_debug .gt. 0) write(*,*)                   &
     &   'Read gzipped binary data file: ', trim(gzip_name)
!
      call open_rd_gzfile(gzip_name)
!
      call gz_read_step_data_mpi_b(my_rank)
      call gz_read_fld_mul_i8head_b(ione, istack_merged)
      call gz_read_fld_inthead_b(fld_IO%num_field_IO)
!
      call gz_read_fld_mul_inthead_b                                    &
     &   (fld_IO%num_field_IO, fld_IO%num_comp_IO)
!
      call gz_read_field_data_mpi_b                                     &
     &   (fld_IO%nnod_IO, fld_IO%num_field_IO, fld_IO%ntot_comp_IO,     &
     &    fld_IO%fld_name, fld_IO%d_IO)
!
      call close_gzfile()
!
      end subroutine gz_read_step_field_file_b
!
! -----------------------------------------------------------------------
!
      subroutine gz_rd_alloc_st_fld_file_b(gzip_name, my_rank, fld_IO)
!
      character(len=kchara), intent(in) :: gzip_name
      integer(kind = kint), intent(in) :: my_rank
!
      type(field_IO), intent(inout) :: fld_IO
!
      integer(kind = kint_gl) :: istack_merged(1)
!
!
      if(my_rank.eq.0 .or. i_debug .gt. 0) write(*,*)                   &
     &   'Read gzipped binary data file: ', trim(gzip_name)
!
      call open_rd_gzfile(gzip_name)
!
      call gz_read_step_data_mpi_b(my_rank)
      call gz_read_fld_mul_i8head_b(ione, istack_merged)
      call gz_read_fld_inthead_b(fld_IO%num_field_IO)
!
      call alloc_phys_name_IO(fld_IO)
      call gz_read_fld_mul_inthead_b                                    &
     &   (fld_IO%num_field_IO, fld_IO%num_comp_IO)
!
      fld_IO%nnod_IO = int(istack_merged(1))
      call cal_istack_phys_comp_IO(fld_IO)
      call alloc_phys_data_IO(fld_IO)
      call gz_read_field_data_mpi_b                                     &
     &   (fld_IO%nnod_IO, fld_IO%num_field_IO, fld_IO%ntot_comp_IO,     &
     &    fld_IO%fld_name, fld_IO%d_IO)
!
      call close_gzfile()
!
      end subroutine gz_rd_alloc_st_fld_file_b
!
! -----------------------------------------------------------------------
!
      subroutine gz_rd_alloc_st_fld_head_b(gzip_name, my_rank, fld_IO)
!
      character(len=kchara), intent(in) :: gzip_name
      integer(kind = kint), intent(in) :: my_rank
!
      type(field_IO), intent(inout) :: fld_IO
!
      integer(kind = kint_gl) :: istack_merged(1)
!
!
      if(my_rank.eq.0 .or. i_debug .gt. 0) write(*,*)                   &
     &   'Read gzipped binary data file: ', trim(gzip_name)
!
      call open_rd_gzfile(gzip_name)
!
      call gz_read_step_data_mpi_b(my_rank)
      call gz_read_fld_mul_i8head_b(ione, istack_merged)
      call gz_read_fld_inthead_b(fld_IO%num_field_IO)
!
      call alloc_phys_name_IO(fld_IO)
      call gz_read_fld_mul_inthead_b                                    &
     &   (fld_IO%num_field_IO, fld_IO%num_comp_IO)
!
      call close_gzfile()
!
      call cal_istack_phys_comp_IO(fld_IO)
!
      end subroutine gz_rd_alloc_st_fld_head_b
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine gz_write_field_data_b(id_rank, nnod, num_field,        &
     &          ntot_comp, ncomp_field, field_name, d_nod)
!
      use m_phys_constants
!
      integer(kind=kint), intent(in) :: id_rank, nnod
      integer(kind=kint), intent(in) :: num_field, ntot_comp
      integer(kind=kint), intent(in) :: ncomp_field(num_field)
      character(len=kchara), intent(in) :: field_name(num_field)
      real(kind = kreal), intent(in) :: d_nod(nnod,ntot_comp)
!
      integer(kind = kint_gl) :: istack_merged(1)
!
!
      call write_step_data_mpi_b(id_rank)
!
      istack_merged(1) = nnod
      call gz_write_fld_mul_i8head_b(ione, istack_merged)
      call gz_write_fld_inthead_b(num_field)
      call gz_write_fld_mul_inthead_b(num_field, ncomp_field)
!
      call gz_write_fld_mul_charhead_b(num_field, field_name)
      call gz_write_fld_realarray2_b(nnod, ntot_comp, d_nod)
!
      end subroutine gz_write_field_data_b
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine write_step_data_mpi_b(id_rank)
!
      integer(kind=kint), intent(in) :: id_rank
!
!
      call gz_write_fld_inthead_b(i_UNIX)
      call gz_write_fld_inthead_b(id_rank)
      call gz_write_fld_inthead_b(i_time_step_IO)
!
      call gz_write_fld_realhead_b(time_IO)
      call gz_write_fld_realhead_b(delta_t_IO)
!
      end subroutine write_step_data_mpi_b
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine gz_read_step_data_mpi_b(my_rank)
!
      integer(kind = kint), intent(in) :: my_rank
      integer(kind = kint) :: id_rank, ie_read
!
!
      call gz_read_fld_inthead_b(ie_read)
      call gz_read_fld_inthead_b(id_rank)
      call gz_read_fld_inthead_b(i_time_step_IO)
      call gz_read_fld_realhead_b(time_IO)
      call gz_read_fld_realhead_b(delta_t_IO)
!
      iflag_endian_swap = 0
      if(i_UNIX .ne. ie_read) then
        if(my_rank .eq. 0) write(*,*)                                   &
     &                            'binary data have opposite endian!'
        iflag_endian_swap = 1
      end if
!
      end subroutine gz_read_step_data_mpi_b
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine gz_read_field_data_mpi_b(nnod, num_field, ncomp,       &
     &          field_name, vect)
!
      integer(kind=kint), intent(in) :: nnod
      integer(kind=kint), intent(in) :: num_field, ncomp
      character(len=kchara), intent(inout) :: field_name(num_field)
      real(kind = kreal), intent(inout) :: vect(nnod,ncomp)
!
!
      call gz_read_fld_mul_charhead_b(num_field, field_name)
      call gz_read_fld_realarray2_b(nnod, ncomp, vect)
!
      end subroutine gz_read_field_data_mpi_b
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine gz_write_fld_inthead_b(int_dat)
!
      integer(kind = kint), intent(in) :: int_dat
!
      integer(kind = kint) :: ierr
!
!
      call gzwrite_f(kint, kint,    ierr)
      call gzwrite_f(kint, int_dat, ierr)
      call gzwrite_f(kint, kint,    ierr)
!
      end subroutine gz_write_fld_inthead_b
!
! -----------------------------------------------------------------------
!
      subroutine gz_write_fld_realhead_b(real_dat)
!
      real(kind = kreal), intent(in) :: real_dat
!
      integer(kind = kint) :: ierr
!
!
      call gzwrite_f(kint,  kreal,    ierr)
      call gzwrite_f(kreal, real_dat, ierr)
      call gzwrite_f(kint,  kreal,    ierr)
!
      end subroutine gz_write_fld_realhead_b
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine gz_write_fld_mul_i8head_b(num, int_gl_dat)
!
      integer(kind = kint), intent(in) :: num
      integer(kind = kint_gl), intent(in) :: int_gl_dat(num)
!
      integer(kind = kint) :: ierr, ilength
!
!
      ilength = num *  kint_gl
      call gzwrite_f(kint,    ilength,    ierr)
      call gzwrite_f(ilength, int_gl_dat, ierr)
      call gzwrite_f(kint,    ilength,    ierr)
!
      end subroutine gz_write_fld_mul_i8head_b
!
! -----------------------------------------------------------------------
!
      subroutine gz_write_fld_mul_inthead_b(num, int_dat)
!
      integer(kind = kint), intent(in) :: num
      integer(kind = kint), intent(in) :: int_dat(num)
!
      integer(kind = kint) :: ierr, ilength
!
!
      ilength = num *  kint
      call gzwrite_f(kint,    ilength, ierr)
      call gzwrite_f(ilength, int_dat, ierr)
      call gzwrite_f(kint,    ilength, ierr)
!
      end subroutine gz_write_fld_mul_inthead_b
!
! -----------------------------------------------------------------------
!
      subroutine gz_write_fld_mul_charhead_b(num, chara_dat)
!
      integer(kind = kint), intent(in) :: num
      character(len=kchara), intent(in) :: chara_dat(num)
!
      integer(kind = kint) :: ierr, ilength
!
!
      ilength = num *  kchara
      call gzwrite_f(kint,    ilength,   ierr)
      call gzwrite_f(ilength, chara_dat, ierr)
      call gzwrite_f(kint,    ilength,   ierr)
!
      end subroutine gz_write_fld_mul_charhead_b
!
! -----------------------------------------------------------------------
!
      subroutine gz_write_fld_realarray2_b(n1, n2, real_dat)
!
      integer(kind = kint), intent(in) :: n1, n2
      real(kind = kreal), intent(in) :: real_dat(n1,n2)
!
      integer(kind = kint) :: ierr, ilength
!
!
      ilength = n1 * n2 * kreal
      call gzwrite_f(kint,    ilength,  ierr)
      call gzwrite_f(ilength, real_dat, ierr)
      call gzwrite_f(kint, ilength,     ierr)
!
      end subroutine gz_write_fld_realarray2_b
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine gz_read_fld_inthead_b(int_dat)
!
      integer(kind = kint), intent(inout) :: int_dat
!
      integer(kind = kint) :: ilength, ierr
!
!
      call gzread_f(kint, ilength, ierr)
      call gzread_f(kint, int_dat, ierr)
      call gzread_f(kint, ilength, ierr)
!
      end subroutine gz_read_fld_inthead_b
!
! -----------------------------------------------------------------------
!
      subroutine gz_read_fld_realhead_b(real_dat)
!
      real(kind = kreal), intent(inout) :: real_dat
!
      integer(kind = kint) :: ilength, ierr
!
!
      call gzread_f(kint,  ilength,  ierr)
      call gzread_f(kreal, real_dat, ierr)
      call gzread_f(kint,  ilength,  ierr)
!
      end subroutine gz_read_fld_realhead_b
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine gz_read_fld_mul_i8head_b(num, int_gl_dat)
!
      integer(kind = kint), intent(in) :: num
      integer(kind = kint_gl), intent(inout) :: int_gl_dat(num)
!
      integer(kind = kint) :: ilength, ierr
!
!
      call gzread_f(kint,    ilength,    ierr)
      ilength = num * kint_gl
      call gzread_f(ilength, int_gl_dat, ierr)
      call gzread_f(kint,    ilength,    ierr)
!
      end subroutine gz_read_fld_mul_i8head_b
!
! -----------------------------------------------------------------------
!
      subroutine gz_read_fld_mul_inthead_b(num, int_dat)
!
      integer(kind = kint), intent(in) :: num
      integer(kind = kint), intent(inout) :: int_dat(num)
!
      integer(kind = kint) :: ilength, ierr
!
!
      call gzread_f(kint,    ilength, ierr)
      ilength = num * kint
      call gzread_f(ilength, int_dat, ierr)
      call gzread_f(kint,    ilength, ierr)
!
      end subroutine gz_read_fld_mul_inthead_b
!
! -----------------------------------------------------------------------
!
      subroutine gz_read_fld_mul_charhead_b(num, chara_dat)
!
      integer(kind = kint), intent(in) :: num
      character(len=kchara), intent(inout) :: chara_dat(num)
!
      integer(kind = kint) :: ilength, ierr
!
!
      call gzread_f(kint,    ilength,   ierr)
      ilength = num * kchara
      call gzread_f(ilength, chara_dat, ierr)
      call gzread_f(kint,    ilength,   ierr)
!
      end subroutine gz_read_fld_mul_charhead_b
!
! -----------------------------------------------------------------------
!
      subroutine gz_read_fld_realarray2_b(n1, n2, real_dat)
!
      integer(kind = kint), intent(in) :: n1, n2
      real(kind = kreal), intent(inout) :: real_dat(n1,n2)
!
      integer(kind = kint) :: ilength, ierr
!
!
      call gzread_f(kint,    ilength,  ierr)
      ilength =  n1 * n2 * kreal
      call gzread_f(ilength, real_dat, ierr)
      call gzread_f(kint,    ilength,  ierr)
!
      end subroutine gz_read_fld_realarray2_b
!
! -----------------------------------------------------------------------
!
      end module gz_field_file_IO_b