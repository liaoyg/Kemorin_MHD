!>@file  udt_file_IO.f90
!!       module udt_file_IO
!!
!! @author H. Matsui
!! @date   Programmed in July, 2006
!
!> @brief UCD format data IO
!!
!!@verbatim
!!      subroutine write_ucd_file(my_rank, istep)
!!      subroutine write_udt_file(my_rank, istep)
!!      subroutine write_grd_file(my_rank)
!!
!!      subroutine read_udt_file(my_rank, istep)
!!      subroutine read_and_alloc_udt_params(my_rank, istep)
!!      subroutine read_and_alloc_udt_file(my_rank, istep)
!!@endverbatim
!!
!!@param my_rank  process ID
!!@param istep    step number for output
!
      module udt_file_IO
!
      use m_precision
      use m_constants
      use m_machine_parameter
!
      use m_field_file_format
      use m_ucd_data
      use set_ucd_file_names
!
      implicit none
!
!>      file ID for UCD file
      integer(kind = kint), parameter, private :: id_ucd_file = 16
!
      private :: write_udt_fields, write_ucd_mesh
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine write_ucd_file(my_rank, istep)
!
      use set_parallel_file_name
!
      integer(kind=kint), intent(in) :: my_rank, istep
      character(len=kchara) :: file_name
!
!
      call set_parallel_ucd_file_name(ucd_header_name, iflag_ucd,       &
     &    my_rank, istep, file_name)
!
      if(my_rank.eq.0 .or. i_debug .gt. 0) write(*,*)                   &
     &     'Write ascii UCD file: ', trim(file_name)
!
      open(id_ucd_file,file=file_name, form='formatted')
      call write_ucd_mesh
      call write_udt_fields
      close(id_ucd_file)
!
      end subroutine write_ucd_file
!
!-----------------------------------------------------------------------
!
      subroutine write_udt_file(my_rank, istep)
!
      use set_parallel_file_name
!
      integer(kind=kint), intent(in) :: my_rank, istep
      character(len=kchara) :: file_name
!
!
      call set_parallel_ucd_file_name(ucd_header_name, iflag_udt,       &
     &    my_rank, istep, file_name)
!
      if(my_rank.eq.0 .or. i_debug .gt. 0) write(*,*)                   &
     &     'Write ascii UCD field: ', trim(file_name)
!
      open(id_ucd_file,file=file_name, form='formatted')
      call write_udt_fields
      close(id_ucd_file)
!
      end subroutine write_udt_file
!
!-----------------------------------------------------------------------
!
      subroutine write_grd_file(my_rank)
!
      use set_parallel_file_name
!
      integer(kind=kint), intent(in) :: my_rank
      character(len=kchara) :: file_name
!
!
      call set_parallel_grd_file_name(ucd_header_name, iflag_udt,       &
     &    my_rank, file_name)
!
      if(my_rank.eq.0 .or. i_debug .gt. 0) write(*,*)                   &
     &     'Write ascii UCD mesh: ', trim(file_name)
!
      open (id_ucd_file, file=file_name, status='replace')
      call write_ucd_mesh
      close(id_ucd_file)
!
      end subroutine write_grd_file
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine write_udt_fields
!
      use udt_data_IO
!
!
!   =====================
!   [1] data on nodes
!   =====================
!    (A) header for fields
!     *  number of fields, number of components, ...
!   ===================================================================
!    (B) nodal data name
!   =====================
      if(num_field_ucd .gt. 0) then
        call write_udt_field_header(id_ucd_file, num_field_ucd,         &
     &      num_comp_ucd, phys_name_ucd)
!
!    (C) fields
!     *  global node ID, results
!   ===================================================
        call write_ucd_field_data(id_ucd_file, nnod_ucd,                &
     &      ntot_comp_ucd, nnod_ucd, inod_gl_ucd, d_nod_ucd)
      end if
!
      end subroutine write_udt_fields
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine read_udt_file(my_rank, istep)
!
      use udt_data_IO
!
      integer(kind=kint), intent(in) :: my_rank, istep
      character(len=kchara) :: file_name
!
!
      call set_parallel_ucd_file_name(ucd_header_name, iflag_udt,       &
     &    my_rank, istep, file_name)
!
      open (id_ucd_file, file=file_name, status='old')
!
      call read_udt_field_header(id_ucd_file, num_field_ucd,            &
     &    num_comp_ucd, phys_name_ucd)
!
      call cal_istack_ucd_component
!
      call read_ucd_field_data(id_ucd_file, nnod_ucd,                   &
     &    ntot_comp_ucd, d_nod_ucd)
!
      close(id_ucd_file)
!
      end subroutine read_udt_file
!
!-----------------------------------------------------------------------
!
      subroutine read_and_alloc_udt_params(my_rank, istep)
!
      use udt_data_IO
!
      integer(kind=kint), intent(in) :: my_rank, istep
      character(len=kchara) :: file_name
!
!
      call set_parallel_ucd_file_name(ucd_header_name, iflag_udt,       &
     &    my_rank, istep, file_name)
!
      open (id_ucd_file, file=file_name, status='old')
      read(id_ucd_file,'(i10)') num_field_ucd
      close(id_ucd_file)
!
      call allocate_ucd_phys_name
!
!
      open (id_ucd_file, file=file_name, status='old')
      call read_udt_field_header(id_ucd_file, num_field_ucd,            &
     &    num_comp_ucd, phys_name_ucd)
      close(id_ucd_file)
!
      call cal_istack_ucd_component
      call allocate_ucd_phys_data
!
      end subroutine read_and_alloc_udt_params
!
! -----------------------------------------------------------------------
!
      subroutine read_and_alloc_udt_file(my_rank, istep)
!
      use udt_data_IO
!
      integer(kind=kint), intent(in) :: my_rank, istep
      character(len=kchara) :: file_name
!
!
      call set_parallel_ucd_file_name(ucd_header_name, iflag_udt,       &
     &    my_rank, istep, file_name)
!
      open (id_ucd_file, file=file_name, status='old')
      read(id_ucd_file,'(i10)') num_field_ucd
      close(id_ucd_file)
!
      call allocate_ucd_phys_name
!
!
      open (id_ucd_file, file=file_name, status='old')
      call read_udt_field_header(id_ucd_file, num_field_ucd,            &
     &    num_comp_ucd, phys_name_ucd)
!
      call cal_istack_ucd_component
      call allocate_ucd_phys_data
!
      call read_ucd_field_data(id_ucd_file, nnod_ucd,                   &
     &    ntot_comp_ucd, d_nod_ucd)
      close(id_ucd_file)
!
      end subroutine read_and_alloc_udt_file
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine write_ucd_mesh
!
      use udt_data_IO
!
!
!   =====================
!   [1] node data
!   =====================
!     * number of data
!     * num. of node, elements, nodal field components, 
!       elemental fields components, and num of model
!   ====================================================
      call write_udt_mesh_header(id_ucd_file, nnod_ucd,                 &
     &    nele_ucd, ntot_comp_ucd)
!
!   =====================
!   [2] position of node
!   =====================
!     global node ID, position
!   ==================================
      call write_ucd_field_data(id_ucd_file, nnod_ucd, ithree,          &
     &    nnod_ucd, inod_gl_ucd, xx_ucd)
!
!   =====================
!   [3] element data
!   =====================
!     * global element ID, node connection (by global ID)
!   ===========================================================
      call write_ucd_mesh_connect(id_ucd_file, nele_ucd,                &
     &    nnod_4_ele_ucd, nele_ucd, iele_gl_ucd, ie_ucd)
!
      end subroutine write_ucd_mesh
!
! -----------------------------------------------------------------------
!
      end module udt_file_IO
