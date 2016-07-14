!>@file   itp_table_file_IO.f90
!!@brief  module itp_table_file_IO
!!
!!@author H. Matsui
!!@date Programmed in Sep. 2006
!
!>@brief ASCII Interpolation table file IO
!!
!!@verbatim
!!      subroutine write_itp_table_file_a(file_name, my_rank, IO_itp_org)
!!      subroutine read_itp_table_file_a                                &
!!     &         (file_name, my_rank, IO_itp_org, ierr)
!!
!!      subroutine write_itp_coefs_dest_file_a(file_name, my_rank)
!!      subroutine read_itp_coefs_dest_file_a(file_name, my_rank, ierr)
!!      subroutine read_itp_table_dest_file_a(file_name, my_rank, ierr)
!!      subroutine read_itp_domain_dest_file_a(file_name, my_rank, ierr)
!!@endverbatim
!
      module itp_table_file_IO
!
      use m_precision
      use m_error_IDs
!
      use m_interpolate_table_dest_IO
      use t_interpolate_tbl_org
!
      use itp_table_data_IO
!
      implicit none
!
      integer(kind = kint), parameter :: id_tbl_file = 19
      private :: id_tbl_file
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine write_itp_table_file_a(file_name, my_rank, IO_itp_org)
!
      character(len=kchara), intent(in) :: file_name
      integer(kind= kint), intent(in) :: my_rank
      type(interpolate_table_org), intent(inout) :: IO_itp_org
!
      open (id_tbl_file, file = file_name, form = 'formatted')
      call write_interpolate_table_dest(id_tbl_file, my_rank)
!
      call write_interpolate_table_org                                  &
     &   (id_tbl_file, my_rank, IO_itp_org)
      call write_interpolate_coefs_org(id_tbl_file, IO_itp_org)
!
      close(id_tbl_file)
!
      if (IO_itp_org%num_dest_domain .gt. 0) then
        call dealloc_itp_table_org(IO_itp_org)
        call dealloc_itp_num_org(IO_itp_org)
      end if
!
      if (IO_itp_dest%num_org_domain .gt. 0) then
        call deallocate_itp_nod_dst_IO
        call deallocate_itp_num_dst_IO
      end if
!
      end subroutine write_itp_table_file_a
!
!-----------------------------------------------------------------------
!
      subroutine read_itp_table_file_a                                  &
     &         (file_name, my_rank, IO_itp_org, ierr)
!
      character(len=kchara), intent(in) :: file_name
      integer(kind = kint), intent(in) :: my_rank
!
      type(interpolate_table_org), intent(inout) :: IO_itp_org
      integer(kind = kint), intent(inout) :: ierr
!
      integer(kind = kint) :: n_rank_file
!
!
      open (id_tbl_file, file = file_name, form = 'formatted')
!        write(*,*) 'read_interpolate_domain_dest', trim(file_name)
      call read_interpolate_domain_dest(id_tbl_file, n_rank_file)
!        write(*,*) 'read_interpolate_table_dest'
      call read_interpolate_table_dest(id_tbl_file)
!
!        write(*,*) 'read_interpolate_domain_org'
      call read_interpolate_domain_org                                  &
     &   (id_tbl_file, n_rank_file, IO_itp_org)
!        write(*,*) 'read_interpolate_table_org'
      call read_interpolate_table_org(id_tbl_file, IO_itp_org)
!        write(*,*) 'read_interpolate_coefs_org'
      call read_interpolate_coefs_org(id_tbl_file, IO_itp_org)
!
      close(id_tbl_file)
!
      ierr = 0
      if (n_rank_file .ne. my_rank) ierr = n_rank_file
!
      end subroutine read_itp_table_file_a
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine write_itp_coefs_dest_file_a(file_name, my_rank)
!
      character(len=kchara), intent(in) :: file_name
      integer(kind= kint), intent(in) :: my_rank
!
!
      open (id_tbl_file, file = file_name, form = 'formatted')
      call write_interpolate_table_dest(id_tbl_file, my_rank)
      call write_interpolate_coefs_dest(id_tbl_file)
      close(id_tbl_file)
!
      if (IO_itp_dest%num_org_domain .gt. 0) then
        call deallocate_itp_coefs_dst_IO
        call deallocate_itp_nod_dst_IO
        call deallocate_itp_num_dst_IO
      end if
!
      end subroutine write_itp_coefs_dest_file_a
!
!-----------------------------------------------------------------------
!
      subroutine read_itp_coefs_dest_file_a(file_name, my_rank, ierr)
!
      character(len=kchara), intent(in) :: file_name
      integer(kind= kint), intent(in) :: my_rank
      integer(kind = kint), intent(inout) :: ierr
!
      integer(kind = kint) :: n_rank_file
!
!
      open (id_tbl_file, file = file_name, form = 'formatted')
      call read_interpolate_domain_dest(id_tbl_file, n_rank_file)
      call read_interpolate_table_dest(id_tbl_file)
      call read_interpolate_coefs_dest(id_tbl_file)
      close(id_tbl_file)
!
      ierr = 0
      if (n_rank_file .ne. my_rank) ierr = ierr_file
!
      end subroutine read_itp_coefs_dest_file_a
!
!-----------------------------------------------------------------------
!
      subroutine read_itp_table_dest_file_a(file_name, my_rank, ierr)
!
      character(len=kchara), intent(in) :: file_name
      integer(kind= kint), intent(in) :: my_rank
      integer(kind = kint), intent(inout) :: ierr
!
      integer(kind = kint) :: n_rank_file
!
!
      open (id_tbl_file, file = file_name, form = 'formatted')
      call read_interpolate_domain_dest(id_tbl_file, n_rank_file)
      call read_interpolate_table_dest(id_tbl_file)
      close(id_tbl_file)
!
      ierr = 0
      if (n_rank_file .ne. my_rank) ierr = ierr_file
!
      end subroutine read_itp_table_dest_file_a
!
!-----------------------------------------------------------------------
!
      subroutine read_itp_domain_dest_file_a(file_name, my_rank, ierr)
!
      character(len=kchara), intent(in) :: file_name
      integer(kind = kint), intent(in) :: my_rank
      integer(kind = kint), intent(inout) :: ierr
!
      integer(kind = kint) :: n_rank_file
!
!
      open (id_tbl_file, file = file_name, form = 'formatted')
      call read_interpolate_domain_dest(id_tbl_file, n_rank_file)
      close(id_tbl_file)
!
      ierr = 0
      if (n_rank_file .ne. my_rank) ierr = ierr_file
!
      end subroutine read_itp_domain_dest_file_a
!
!-----------------------------------------------------------------------
!
      end module itp_table_file_IO
