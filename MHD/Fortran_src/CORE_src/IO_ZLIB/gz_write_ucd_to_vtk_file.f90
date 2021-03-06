!>@file   gz_write_ucd_to_vtk_file.f90
!!@brief  module gz_write_ucd_to_vtk_file
!!
!!@author H. Matsui
!!@date    programmed by H.Matsui on July, 2006
!!@n       Modified by H.Matsui on March, 2013
!
!>@brief Output FEM field data to distributed VTK file (gzipped)
!!
!!@verbatim
!!      subroutine write_gz_parallel_vtk_file(my_rank, nprocs,          &
!!     &          istep, file_prefix)
!!
!!      subroutine write_ucd_data_2_gz_vtk(my_rank, gzip_name, ucd)
!!      subroutine write_ucd_data_2_gz_vtk_phys(my_rank, gzip_name, ucd)
!!      subroutine write_ucd_data_2_gz_vtk_grid(my_rank, gzip_name, ucd)
!!@endverbatim
!!
!!@param my_rank    subdomain ID
!!@param gzip_name  file name
!
      module gz_write_ucd_to_vtk_file
!
      use m_precision
      use m_constants
      use m_machine_parameter
      use m_field_file_format
!
      use t_ucd_data
      use set_ucd_file_names
      use skip_gz_comment
!
      implicit none
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine write_gz_parallel_vtk_file(my_rank, nprocs,            &
     &          istep, file_prefix)
!
      use set_parallel_file_name
      use skip_gz_comment
!
      character(len=kchara), intent(in) :: file_prefix
      integer(kind=kint), intent(in) :: my_rank, nprocs, istep
!
      character(len=kchara) :: gzip_name, file_name, fname_tmp
      character(len=kchara) :: fname_nodir
      integer(kind = kint) :: ip
!
!
      if(my_rank .gt. 0) return
!
      call delete_directory_name(file_prefix, fname_nodir)
      call add_int_suffix(istep, file_prefix, fname_tmp)
      call add_pvtk_extension(fname_tmp, file_name)
      call add_gzip_extension(file_name, gzip_name)
!
      write(*,*) 'Write gzipped parallel VTK file: ', trim(gzip_name)
      call open_wt_gzfile_f(gzip_name)
!
      write(textbuf,'(a,a1)') '<File version="pvtk-1.0"', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,a1)')                                           &
     &     '       dataType="vtkUnstructuredGrid"', char(0)
      call gz_write_textbuf_w_lf
      write(textbuf,'(a,i6,a,a1)')                                      &
     &     '       numberOfPieces="', nprocs, '" >', char(0)
      call gz_write_textbuf_w_lf
      do ip = 0, nprocs-1
        call set_parallel_ucd_file_name(fname_nodir, iflag_vtk,         &
     &      ip, istep, file_name)
        write(textbuf,'(3a,a1)') '   <Piece fileName="',                &
     &                       trim(file_name), '" />', char(0)
        call gz_write_textbuf_w_lf
      end do
      write(textbuf,'(a,a1)') '</File>', char(0)
      call gz_write_textbuf_w_lf
!
      call close_gzfile_f
!
      end subroutine write_gz_parallel_vtk_file
!
!  ---------------------------------------------------------------------
!
      subroutine write_ucd_data_2_gz_vtk(my_rank, gzip_name, ucd)
!
      use gz_vtk_file_IO
!
      character(len=kchara), intent(in) :: gzip_name
      integer(kind = kint), intent(in) ::  my_rank
      type(ucd_data), intent(in) :: ucd
!
!
      if(my_rank.le.0) write(*,*)                                       &
     &    'Write gzipped VTK data: ', trim(gzip_name)
!
      call write_gz_vtk_file(gzip_name,                                 &
     &    ucd%nnod, ucd%nele, ucd%nnod_4_ele, ucd%xx, ucd%ie,           &
     &    ucd%num_field, ucd%ntot_comp, ucd%num_comp, ucd%phys_name,    &
     &    ucd%d_ucd)
!
      end subroutine write_ucd_data_2_gz_vtk
!
! -----------------------------------------------------------------------
!
      subroutine write_ucd_data_2_gz_vtk_phys(my_rank, gzip_name, ucd)
!
      use gz_vtk_file_IO
!
      character(len=kchara), intent(in) :: gzip_name
      integer(kind = kint), intent(in) ::  my_rank
      type(ucd_data), intent(in) :: ucd
!
!
      if(my_rank.le.0) write(*,*)                                       &
     &    'Write gzipped VTK field: ', trim(gzip_name)
!
      call write_gz_vtk_phys(gzip_name, ucd%nnod, ucd%num_field,        &
     &    ucd%ntot_comp, ucd%num_comp, ucd%phys_name, ucd%d_ucd)
!
      end subroutine write_ucd_data_2_gz_vtk_phys
!
! -----------------------------------------------------------------------
!
      subroutine write_ucd_data_2_gz_vtk_grid(my_rank, gzip_name, ucd)
!
      use gz_vtk_file_IO
!
      character(len=kchara), intent(in) :: gzip_name
      integer(kind = kint), intent(in) ::  my_rank
      type(ucd_data), intent(in) :: ucd
!
!
      if(my_rank.le.0) write(*,*)                                       &
     &    'Write gzipped VTK grid: ', trim(gzip_name)
!
      call write_gz_vtk_grid(gzip_name,                                 &
     &    ucd%nnod, ucd%nele, ucd%nnod_4_ele, ucd%xx, ucd%ie)
!
      end subroutine write_ucd_data_2_gz_vtk_grid
!
! -----------------------------------------------------------------------
!
      end module gz_write_ucd_to_vtk_file
