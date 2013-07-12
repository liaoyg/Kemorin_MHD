!>@file  gz_vtk_file_IO.f90
!!       module gz_vtk_file_IO
!!
!!@author H. Matsui
!!@date   Programmed by H. Matsui in Feb., 2007
!
!> @brief Output routine for gzipped VTK data segments
!!
!!@verbatim
!!      subroutine write_gz_vtk_file(gzip_name,                         &
!!     &          nnod, nele, nnod_ele, xx, ie, num_field,  ntot_comp,  &
!!     &          ncomp_field, field_name, d_nod)
!!      subroutine write_gz_vtk_phys(gzip_name, nnod, num_field,        &
!!     &          ntot_comp, ncomp_field, field_name, d_nod)
!!      subroutine write_gz_vtk_grid(gzip_name, nnod, nele, nnod_ele,   &
!!     &          xx, ie)
!!@endverbatim
!
      module gz_vtk_file_IO
!
      use m_precision
      use m_constants
!
      use gz_vtk_data_IO
!
      implicit none
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine write_gz_vtk_file(gzip_name,                           &
     &          nnod, nele, nnod_ele, xx, ie, num_field,  ntot_comp,    &
     &          ncomp_field, field_name, d_nod)
!
      use set_parallel_file_name
!
      character(len=kchara), intent(in) :: gzip_name
!
      integer(kind = kint), intent(in) :: nnod, nele
      integer(kind = kint), intent(in) :: nnod_ele
      integer(kind = kint), intent(in) :: ie(nele,nnod_ele)
      real(kind = kreal), intent(in) :: xx(nnod,3)
!
      integer (kind=kint), intent(in) :: num_field, ntot_comp
      integer(kind=kint ), intent(in) :: ncomp_field(num_field)
      character(len=kchara), intent(in) :: field_name(num_field)
      real(kind = kreal), intent(in) :: d_nod(nnod,ntot_comp)
!
!
      call open_wt_gzfile(gzip_name)
      call write_gz_vtk_mesh(nnod, nele, nnod_ele, xx, ie)
      call write_gz_vtk_data(nnod, num_field, ntot_comp, ncomp_field,   &
     &    field_name, d_nod)
!
      call close_gzfile
!
      end subroutine write_gz_vtk_file
!
! -----------------------------------------------------------------------
!
      subroutine write_gz_vtk_phys(gzip_name, nnod, num_field,          &
     &          ntot_comp, ncomp_field, field_name, d_nod)
!
      use set_parallel_file_name
!
      character(len=kchara), intent(in) :: gzip_name
!
      integer (kind=kint), intent(in) :: nnod
      integer (kind=kint), intent(in) :: num_field, ntot_comp
      integer(kind=kint ), intent(in) :: ncomp_field(num_field)
      character(len=kchara), intent(in) :: field_name(num_field)
      real(kind = kreal), intent(in) :: d_nod(nnod,ntot_comp)
!
!
      call open_wt_gzfile(gzip_name)
      call write_gz_vtk_data(nnod, num_field, ntot_comp, ncomp_field,   &
     &    field_name, d_nod)
      call close_gzfile
!
      end subroutine write_gz_vtk_phys
!
! -----------------------------------------------------------------------
!
      subroutine write_gz_vtk_grid(gzip_name, nnod, nele, nnod_ele,     &
     &          xx, ie)
!
      use set_parallel_file_name
!
      character(len=kchara), intent(in) :: gzip_name
!
      integer(kind = kint), intent(in) :: nnod, nele
      integer(kind = kint), intent(in) :: nnod_ele
      integer(kind = kint), intent(in) :: ie(nele,nnod_ele)
      real(kind = kreal), intent(in) :: xx(nnod,3)
!
!
      call open_wt_gzfile(gzip_name)
      call write_gz_vtk_mesh(nnod, nele, nnod_ele, xx, ie)
      call close_gzfile
!
      end subroutine write_gz_vtk_grid
!
! -----------------------------------------------------------------------
!
      end module gz_vtk_file_IO
