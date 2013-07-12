!>@file  gz_udt_type_file_IO.f90
!!       module gz_udt_type_file_IO
!!
!! @author H. Matsui
!! @date   Programmed in July, 2006
!
!> @brief gzipped UCD format data IO
!!
!!@verbatim
!!      subroutine write_gz_ucd_type_file(my_rank, istep, ucd)
!!      subroutine write_gz_udt_type_file(my_rank, istep, ucd)
!!      subroutine write_gz_grd_type_file(my_rank, ucd)
!!
!!      subroutine read_udt_type_file_gz(my_rank, istep, ucd)
!!      subroutine read_alloc_udt_type_head_gz(my_rank, istep, ucd)
!!      subroutine read_alloc_udt_type_file_gz(my_rank, istep, ucd)
!!        type(ucd_data), intent(inout) :: ucd
!!@endverbatim
!!
!!@param my_rank  process ID
!!@param istep    step number for output
!!@param ucd      Structure for FEM field data IO
!
      module gz_udt_type_file_IO
!
      use m_precision
      use m_constants
      use m_machine_parameter
      use m_field_file_format
!
      use t_ucd_data
!
      use set_ucd_file_names
      use skip_gz_comment
      use gz_ucd_data_IO
!
      implicit none
!
      private :: write_gz_udt_type_fields, write_gz_ucd_type_mesh
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine write_gz_ucd_type_file(my_rank, istep, ucd)
!
      integer(kind=kint), intent(in) :: my_rank, istep
      type(ucd_data), intent(in) :: ucd
!
      character(len=kchara) :: gzip_name
!
!
      call set_parallel_ucd_file_name(ucd%header_name, iflag_ucd_gz,    &
     &    my_rank, istep, gzip_name)
!
      if(i_debug.gt.0 .or. my_rank.eq.0) write(*,*)                     &
     &      'Write gzipped ucd file: ', trim(gzip_name)
      call open_wt_gzfile(gzip_name)
!
      call write_gz_ucd_type_mesh(ucd)
      call write_gz_udt_type_fields(ucd)
      call close_gzfile
!
      end subroutine write_gz_ucd_type_file
!
!-----------------------------------------------------------------------
!
      subroutine write_gz_udt_type_file(my_rank, istep, ucd)
!
      integer(kind=kint), intent(in) :: my_rank, istep
      type(ucd_data), intent(in) :: ucd
!
      character(len=kchara) :: gzip_name
!
!
      call set_parallel_ucd_file_name(ucd%header_name, iflag_udt_gz,    &
     &    my_rank, istep, gzip_name)
!
      if(i_debug.gt.0 .or. my_rank.eq.0) write(*,*)                     &
     &      'Write gzipped ucd file: ', trim(gzip_name)
      call open_wt_gzfile(gzip_name)
!
      call write_gz_udt_type_fields(ucd)
      call close_gzfile
!
      end subroutine write_gz_udt_type_file
!
!-----------------------------------------------------------------------
!
      subroutine write_gz_grd_type_file(my_rank, ucd)
!
      integer(kind=kint), intent(in) :: my_rank
      type(ucd_data), intent(in) :: ucd
!
      character(len=kchara) :: gzip_name
!
!
      call set_parallel_grd_file_name(ucd%header_name, iflag_udt_gz,    &
     &    my_rank, gzip_name)
!
      if(i_debug.gt.0 .or. my_rank.eq.0) write(*,*)                     &
     &      'Write gzipped ucd grid file: ', trim(gzip_name)
      call open_wt_gzfile(gzip_name)
!
      call write_gz_ucd_type_mesh(ucd)
      call close_gzfile
!
      end subroutine write_gz_grd_type_file
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine read_udt_type_file_gz(my_rank, istep, ucd)
!
      integer(kind=kint), intent(in) :: my_rank, istep
      type(ucd_data), intent(inout) :: ucd
!
      character(len=kchara) :: gzip_name
!
!
      call set_parallel_ucd_file_name(ucd%header_name, iflag_udt_gz,    &
     &    my_rank, istep, gzip_name)
!
      if(i_debug.gt.0 .or. my_rank.eq.0) write(*,*)                     &
     &     'Read gzipped data file: ', trim(gzip_name)
!
      call open_rd_gzfile(gzip_name)
!
      call read_gz_udt_field_header(ucd%num_field, ucd%num_comp,        &
     &    ucd%phys_name)
!
      call cal_istack_ucd_comp_type(ucd)
!
      call read_gz_single_udt_data(ucd%nnod, ucd%ntot_comp, ucd%d_ucd)
!
      call close_gzfile
!
      end subroutine read_udt_type_file_gz
!
!-----------------------------------------------------------------------
!
      subroutine read_alloc_udt_type_head_gz(my_rank, istep, ucd)
!
      integer(kind=kint), intent(in) :: my_rank, istep
      type(ucd_data), intent(inout) :: ucd
!
      character(len=kchara) :: gzip_name
!
!
      call set_parallel_ucd_file_name(ucd%header_name, iflag_udt_gz,    &
     &    my_rank, istep, gzip_name)
!
      if(i_debug.gt.0 .or. my_rank.eq.0) write(*,*)                     &
     &     'Read gzipped data file: ', trim(gzip_name)
!
      call open_rd_gzfile(gzip_name)
!
      call read_gz_udt_field_num(ucd%num_field)
      call alloc_ucd_phys_name_t(ucd)
!
      call read_gz_udt_field_name(ucd%num_field, ucd%num_comp,          &
     &    ucd%phys_name)
!
      call close_gzfile
!
      call cal_istack_ucd_comp_type(ucd)
      call alloc_ucd_phys_data_t(ucd)
!
      end subroutine read_alloc_udt_type_head_gz
!
!-----------------------------------------------------------------------
!
      subroutine read_alloc_udt_type_file_gz(my_rank, istep, ucd)
!
      integer(kind=kint), intent(in) :: my_rank, istep
      type(ucd_data), intent(inout) :: ucd
!
      character(len=kchara) :: gzip_name
!
!
      call set_parallel_ucd_file_name(ucd%header_name, iflag_udt_gz,    &
     &    my_rank, istep, gzip_name)
!
      if(i_debug.gt.0 .or. my_rank.eq.0) write(*,*)                     &
     &     'Write gzipped data file: ', trim(gzip_name)
      call open_rd_gzfile(gzip_name)
!
      call read_gz_udt_field_num(ucd%num_field)
      call alloc_ucd_phys_name_t(ucd)
!
      call read_gz_udt_field_name(ucd%num_field, ucd%num_comp,          &
     &    ucd%phys_name)
!
      call cal_istack_ucd_comp_type(ucd)
      call alloc_ucd_phys_data_t(ucd)
!
      call read_gz_single_udt_data(ucd%nnod, ucd%ntot_comp, ucd%d_ucd)
      call close_gzfile
!
      end subroutine read_alloc_udt_type_file_gz
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine write_gz_udt_type_fields(ucd)
!
      type(ucd_data), intent(in) :: ucd
!
!
      if(ucd%num_field .gt. 0) then
        call write_gz_udt_field_header(ucd%num_field,                   &
     &      ucd%num_comp, ucd%phys_name)
        call write_gz_ucd_field_data(ucd%nnod, ucd%ntot_comp,           &
     &      ucd%nnod, ucd%inod_global, ucd%d_ucd)
      end if
!
      end subroutine write_gz_udt_type_fields
!
!-----------------------------------------------------------------------
!
      subroutine write_gz_ucd_type_mesh(ucd)
!
      type(ucd_data), intent(in) :: ucd
!
!
      call write_gz_udt_mesh_header(ucd%nnod, ucd%nele, ucd%ntot_comp)
!
      call write_gz_ucd_field_data(ucd%nnod, ithree, ucd%nnod,          &
     &    ucd%inod_global, ucd%xx)
      call write_gz_ucd_mesh_connect(ucd%nele, ucd%nnod_4_ele,          &
     &    ucd%nele, ucd%iele_global, ucd%ie)
!
      end subroutine write_gz_ucd_type_mesh
!
! -----------------------------------------------------------------------
!
      end module gz_udt_type_file_IO
