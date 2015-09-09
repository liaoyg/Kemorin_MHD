!>@file  ucd_file_MPI_IO.f90
!!       module ucd_file_MPI_IO
!!
!!@author H. Matsui
!!@date   Programmed in Ma, 2015
!
!> @brief Output merged VTK file usgin MPI-IO
!!
!!@verbatim
!!      subroutine write_ucd_file_mpi(file_name, ucd, m_ucd)
!!      subroutine write_ucd_phys_mpi(file_name, ucd, m_ucd)
!!      subroutine write_ucd_grid_mpi(file_name, ucd, m_ucd)
!!@endverbatim
!
      module ucd_file_MPI_IO
!
      use m_precision
      use m_constants
!
      use calypso_mpi
      use m_calypso_mpi_IO
      use t_ucd_data
      use vtk_data_to_buffer
      use ucd_data_to_buffer
!
      implicit none
!
      private :: write_ucd_data_mpi
      private :: write_ucd_node_mpi, write_ucd_connect_mpi
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine write_ucd_file_mpi(file_name, ucd, m_ucd)
!
      character(len=kchara), intent(in) :: file_name
!
      type(ucd_data), intent(in) :: ucd
      type(merged_ucd_data), intent(in) :: m_ucd
!
      integer :: id_vtk
      integer(kind = kint_gl) :: ioff_gl
!
!
      call calypso_mpi_write_file_open(file_name, nprocs, id_vtk)
!
      ioff_gl = 0
      call write_ucd_node_mpi(id_vtk, ioff_gl,                          &
     &    ucd%nnod, ucd%ntot_comp, ucd%xx,                              &
     &    m_ucd%istack_merged_intnod, m_ucd%istack_merged_ele)
      call write_ucd_connect_mpi(id_vtk, ioff_gl,                       &
     &    ucd%nele, ucd%nnod_4_ele, ucd%ie, m_ucd%istack_merged_ele)
!
      call write_ucd_data_mpi(id_vtk, ioff_gl,                          &
     &    ucd%nnod, ucd%num_field, ucd%ntot_comp, ucd%num_comp,         &
     &    ucd%phys_name, ucd%d_ucd, m_ucd%istack_merged_intnod)
!
      call calypso_close_mpi_file(id_vtk)
!
      end subroutine write_ucd_file_mpi
!
! -----------------------------------------------------------------------
!
      subroutine write_ucd_phys_mpi(file_name, ucd, m_ucd)
!
      character(len=kchara), intent(in) :: file_name
!
      type(ucd_data), intent(in) :: ucd
      type(merged_ucd_data), intent(in) :: m_ucd
!
      integer :: id_vtk
      integer(kind = kint_gl) :: ioff_gl
!
!
      call calypso_mpi_write_file_open(file_name, nprocs, id_vtk)
!
      ioff_gl = 0
      call write_ucd_data_mpi(id_vtk, ioff_gl,                          &
     &    ucd%nnod, ucd%num_field, ucd%ntot_comp, ucd%num_comp,         &
     &    ucd%phys_name, ucd%d_ucd, m_ucd%istack_merged_intnod)
!
      call calypso_close_mpi_file(id_vtk)
!
      end subroutine write_ucd_phys_mpi
!
! -----------------------------------------------------------------------
!
      subroutine write_ucd_grid_mpi(file_name, ucd, m_ucd)
!
      character(len=kchara), intent(in) :: file_name
!
      type(ucd_data), intent(in) :: ucd
      type(merged_ucd_data), intent(in) :: m_ucd
!
      integer :: id_vtk
      integer(kind = kint_gl) :: ioff_gl
!
!
      call calypso_mpi_write_file_open(file_name, nprocs, id_vtk)
!
      ioff_gl = 0
      call write_ucd_node_mpi(id_vtk, ioff_gl,                          &
     &    ucd%nnod, ucd%ntot_comp, ucd%xx,                              &
     &    m_ucd%istack_merged_intnod, m_ucd%istack_merged_ele)
      call write_ucd_connect_mpi(id_vtk, ioff_gl,                       &
     &    ucd%nele, ucd%nnod_4_ele, ucd%ie, m_ucd%istack_merged_ele)
!
      call calypso_close_mpi_file(id_vtk)
!
      end subroutine write_ucd_grid_mpi
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine write_ucd_data_mpi(id_vtk, ioff_gl,                    &
     &          nnod, num_field, ntot_comp, ncomp_field,                &
     &          field_name, d_nod, istack_merged_intnod)
!
      use m_phys_constants
!
      integer(kind = kint_gl), intent(inout) :: ioff_gl
      integer(kind = kint_gl), intent(in)                               &
     &         :: istack_merged_intnod(0:nprocs)
      integer(kind=kint_gl), intent(in) :: nnod
      integer(kind=kint), intent(in) :: num_field, ntot_comp
      integer(kind=kint), intent(in) :: ncomp_field(num_field)
      character(len=kchara), intent(in) :: field_name(num_field)
      real(kind = kreal), intent(in) :: d_nod(nnod,ntot_comp)
!
      integer, intent(in) ::  id_vtk
!
      integer(kind = kint) :: j
      integer(kind = kint_gl) :: inod, num, nt_nod, inod_gl
      integer(kind = MPI_OFFSET_KIND) :: ioffset
      real(kind = kreal)  :: dat_1(ntot_comp)
!
      integer(kind = kint) :: ilen_n
      character(len=16+ntot_comp*23+1), allocatable, target             &
     &                                :: textbuf_n(:)
!
      ilen_n = 16+ntot_comp*23+1
      nt_nod = istack_merged_intnod(nprocs)
      num = istack_merged_intnod(my_rank+1)                             &
     &     - istack_merged_intnod(my_rank)
!
      call calypso_mpi_seek_write_head_c(id_vtk, ioff_gl,               &
     &   ucd_num_comps(num_field, ncomp_field))
!
      do j = 1, num_field
        call calypso_mpi_seek_write_head_c(id_vtk, ioff_gl,             &
     &      ucd_field_name(field_name(j)))
      end do
!
      ioffset = int(ioff_gl + ilen_n * istack_merged_intnod(my_rank))
      ioff_gl = ioff_gl + ilen_n * nt_nod
!
      if(num .le. 0) return
!
      allocate(textbuf_n(num))
      do inod = 1, num
        inod_gl =    inod + istack_merged_intnod(my_rank)
        dat_1(1:ntot_comp) = d_nod(inod,1:ntot_comp)
        textbuf_n(inod) = ucd_each_field(inod_gl, ntot_comp, dat_1)
      end do
      call calypso_mpi_seek_wrt_mul_chara(id_vtk, ioffset, ilen_n,      &
     &    num, textbuf_n)
      deallocate(textbuf_n)
!
      end subroutine write_ucd_data_mpi
!
! -----------------------------------------------------------------------
!
      subroutine write_ucd_node_mpi(id_vtk, ioff_gl, nnod,              &
     &          ntot_comp, xx, istack_merged_intnod, istack_merged_ele)
!
      use m_phys_constants
!
      integer(kind = kint_gl), intent(inout) :: ioff_gl
      integer(kind = kint_gl), intent(in)                               &
     &         :: istack_merged_intnod(0:nprocs)
      integer(kind = kint_gl), intent(in)                               &
     &         :: istack_merged_ele(0:nprocs)
      integer(kind = kint), intent(in) :: ntot_comp
      integer(kind = kint_gl), intent(in) :: nnod
      real(kind = kreal), intent(in) :: xx(nnod,3)
!
      integer, intent(in) ::  id_vtk
!
      integer(kind = kint_gl) :: inod_gl
      integer(kind = kint_gl) :: inod, nt_nod, nt_ele, num
      real(kind = kreal)  :: dat_1(n_vector)
!
      integer(kind = MPI_OFFSET_KIND) :: ioffset
!
      integer(kind = kint), parameter :: ilen_n = 16+n_vector*23+1
      character(len=ilen_n), allocatable, target  :: textbuf_n(:)
!
!
      nt_nod = istack_merged_intnod(nprocs)
      nt_ele = istack_merged_ele(nprocs)
      num = istack_merged_intnod(my_rank+1)                             &
     &     - istack_merged_intnod(my_rank)
!
      call calypso_mpi_seek_write_head_c(id_vtk, ioff_gl,               &
     &   ucd_connect_head(nt_nod, nt_ele, ntot_comp))

      ioffset = int(ioff_gl + ilen_n * istack_merged_intnod(my_rank))
      ioff_gl = ioff_gl + ilen_n * nt_nod
!
      if(num .le. 0) return
!
      allocate(textbuf_n(num))
      do inod = 1, num
        inod_gl =    inod + istack_merged_intnod(my_rank)
        dat_1(1:n_vector) = xx(inod,1:n_vector)
        textbuf_n(inod) = ucd_each_field(inod_gl, n_vector, dat_1)
      end do
      call calypso_mpi_seek_wrt_mul_chara(id_vtk, ioffset, ilen_n,      &
     &    num, textbuf_n)
      deallocate(textbuf_n)
!
      end subroutine write_ucd_node_mpi
!
! -----------------------------------------------------------------------
!
      subroutine write_ucd_connect_mpi(id_vtk, ioff_gl,                 &
     &          nele, nnod_ele, ie, istack_merged_ele)
!
      use m_phys_constants
!
      integer(kind = kint_gl), intent(inout) :: ioff_gl
      integer(kind = kint_gl), intent(in)                               &
     &         :: istack_merged_ele(0:nprocs)
      integer(kind = kint), intent(in) :: nnod_ele
      integer(kind = kint_gl), intent(in) :: nele
      integer(kind = kint_gl), intent(in) :: ie(nele,nnod_ele)
!
      integer, intent(in) ::  id_vtk
!
      integer(kind = kint_gl) :: ie0(nnod_ele), iele_gl
      integer(kind = kint_gl) :: iele, nt_ele
!
      integer(kind = MPI_OFFSET_KIND) :: ioffset
!
      integer(kind = kint) :: ilen_e
      character(len=16+3+6+16*nnod_ele+1), allocatable, target          &
     &                                :: textbuf_e(:)
!
!
      nt_ele = istack_merged_ele(nprocs)
!
      ilen_e = 16+3+6+16*nnod_ele+1
      ioffset = int(ioff_gl + ilen_e * istack_merged_ele(my_rank))
      ioff_gl = ioff_gl + ilen_e * nt_ele
!
      if(nele .le. 0) return
!
      allocate(textbuf_e(nele))
      do iele = 1, nele
        iele_gl = iele + istack_merged_ele(my_rank)
        ie0(1:nnod_ele) = ie(iele,1:nnod_ele)
        textbuf_e(iele) = ucd_each_connect(iele_gl, nnod_ele,ie0)
      end do
      call calypso_mpi_seek_wrt_mul_chara(id_vtk, ioffset, ilen_e,      &
     &    nele, textbuf_e)
      deallocate(textbuf_e)
!
      end subroutine write_ucd_connect_mpi
!
! -----------------------------------------------------------------------
!
      end module ucd_file_MPI_IO
