!set_med_patch_4_element_grp.f90
!      module set_med_patch_4_element_grp
!
!      Written by H. Matsui on March, 2012
!
!      subroutine set_med_patch_ele_grp(file_head)
!
!      data format
!       number of element groups
!         number of node and element for each group
!         node id and node position
!         ....
!         element id and element connection
!         ....
!     (node and element id starts for each group)
!
      module set_med_patch_4_element_grp
!
      use m_precision
      use m_constants
!
      implicit none
!
      integer(kind = kint), parameter :: id_file = 20
!
      integer(kind = kint) :: npatch_grp, npatch_g
      integer(kind = kint), allocatable :: npatch_l(:)
      integer(kind = kint), allocatable :: istack_npatch_l(:)
      real(kind = kreal), allocatable :: xyz_med(:,:)
!
      integer(kind = kint), allocatable :: mark_elegrp(:)
!
      integer, allocatable :: sta1(:,:)
      integer, allocatable :: sta2(:,:)
      integer, allocatable :: req1(:  )
      integer, allocatable :: req2(:  )
!
      private :: id_file, mark_elegrp
      private :: npatch_grp, npatch_g
      private :: npatch_l, istack_npatch_l, xyz_med
      private :: sta1, sta2, req1, req2
!
      private :: count_med_grp_patch, set_med_grp_patch
      private :: init_ele_grp_med_patch, clean_ele_grp_med_patch
      private :: collect_ele_grp_patch
      private :: set_triangle_at_meridian
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine set_med_patch_ele_grp(file_head)
!
      use calypso_mpi
      use m_geometry_constants
      use m_geometry_parameter
      use m_geometry_data
      use m_element_group
      use m_psf_case_table
      use set_psf_case_table
      use patch_4_psf
!
      use t_psf_geometry_list
!
      character(len=kchara), intent(in) :: file_head
!
      integer(kind = kint) :: ip, igrp, inum
      type(sect_search_list) :: ele_search
!
!
      call init_ele_grp_med_patch                                       &
     &   (file_head, ele_grp1%num_grp, ele_grp1%num_item)
!
      call set_sectioning_case_table
!
      do igrp = 1, ele_grp1%num_grp
        call alloc_num_psf_search_list(np_smp, ele_search)
!
        do ip = 1, np_smp
          ele_search%istack_search_smp(ip)                              &
     &            = imat_smp_stack(ip+(igrp-1)*np_smp)                  &
     &             - imat_smp_stack((igrp-1)*np_smp)
        end do
        ele_search%num_search = ele_grp1%istack_grp(igrp)               &
     &                         - ele_grp1%istack_grp(igrp-1)
!
        call alloc_psf_search_list(ele_search)
        do inum = 1, ele_search%num_search
          ele_search%id_search(inum)                                    &
     &            = ele_grp1%item_grp(inum+ele_grp1%istack_grp(igrp-1))
        end do
        call dealloc_psf_search_list(ele_search)
        call dealloc_num_psf_search_list(ele_search)
!
        call set_psf_type_id(numnod, numele, nnod_4_ele, ie,            &
     &      ele_search, mark_elegrp, xx(1,2))
!
        call count_med_grp_patch(igrp, npatch_grp)
        allocate(xyz_med(3,3*npatch_grp+3))
!
        call set_med_grp_patch(igrp, npatch_grp, xyz_med)
!
!   collect patch data
!
        call collect_ele_grp_patch(ele_grp1%grp_name(igrp))
!
        deallocate(xyz_med)
      end do
!
      call clean_ele_grp_med_patch
!
      end subroutine set_med_patch_ele_grp
!
!  ---------------------------------------------------------------------
!
      subroutine count_med_grp_patch(igrp, npatch)
!
      use m_geometry_constants
      use m_geometry_parameter
      use m_geometry_data
      use m_element_group
      use m_psf_case_table
!
      integer(kind = kint), intent(in) :: igrp
      integer(kind = kint), intent(inout) :: npatch
!
      integer(kind = kint) :: ie_1ele(8), iflag
      real(kind = kreal) :: xx_1ele(3,8)
!
      integer(kind = kint) :: ist, ied, inum, iele, mark
      integer(kind = kint) :: n, n_patch
      real(kind = kreal) :: x_tmp(3,3)
!
!
      npatch = 0
      ist = ele_grp1%istack_grp(igrp-1) + 1
      ied = ele_grp1%istack_grp(igrp)
      do inum = ist, ied
        iele = ele_grp1%item_grp(inum)
        ie_1ele(1:8) = ie(iele,1:8)
        xx_1ele(1,1:8) = xx(ie_1ele(1:8),1)
        xx_1ele(2,1:8) = xx(ie_1ele(1:8),2)
        xx_1ele(3,1:8) = xx(ie_1ele(1:8),3)
!
        mark = mark_elegrp(inum) * int(e_multi(iele))
        n_patch = psf_case_tbl(mark)%npatch
        if(n_patch .gt. 0) then
          do n = 1, n_patch
            call set_triangle_at_meridian(xx_1ele, mark, n,             &
     &            iflag, x_tmp)
            npatch = npatch+iflag
          end do
        end if
      end do
!
      end subroutine count_med_grp_patch
!
!  ---------------------------------------------------------------------
!
      subroutine set_med_grp_patch(igrp, npatch, x_patch)
!
      use m_geometry_constants
      use m_geometry_parameter
      use m_geometry_data
      use m_element_group
      use m_psf_case_table
!
      integer(kind = kint), intent(in) :: igrp
      integer(kind = kint), intent(in) :: npatch
      real(kind = kreal), intent(inout) :: x_patch(3,3*npatch)
!
      integer(kind = kint) :: ie_1ele(8), iflag
      real(kind = kreal) :: xx_1ele(3,8)
!
      integer(kind = kint) :: ist, ied, inum, iele, mark, icou
      integer(kind = kint) :: n, n_patch
!
!
      icou = 0
      ist = ele_grp1%istack_grp(igrp-1) + 1
      ied = ele_grp1%istack_grp(igrp)
      do inum = ist, ied
        iele = ele_grp1%item_grp(inum)
        ie_1ele(1:8) =   ie(iele,1:8)
        xx_1ele(1,1:8) = xx(ie_1ele(1:8),1)
        xx_1ele(2,1:8) = xx(ie_1ele(1:8),2)
        xx_1ele(3,1:8) = xx(ie_1ele(1:8),3)
!
        mark = mark_elegrp(inum) * int(e_multi(iele))
        n_patch = psf_case_tbl(mark)%npatch
        if(n_patch .gt. 0) then
          do n = 1, n_patch
            call set_triangle_at_meridian(xx_1ele, mark, n,             &
     &            iflag, x_patch(1,icou+1) )
            icou = icou + 3*iflag
          end do
        end if
      end do
!
      end subroutine set_med_grp_patch
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine set_triangle_at_meridian(xx_1ele, mark_psf, i_1ele,    &
     &          iflag, x_patch)
!
      use m_geometry_constants
      use m_psf_case_table
!
      integer(kind = kint), intent(in) :: mark_psf, i_1ele
      real(kind = kreal), intent(in) :: xx_1ele(3,8)
!
      integer(kind = kint), intent(inout) :: iflag
      real(kind = kreal), intent(inout) :: x_patch(3,3)
!
      integer(kind = kint) :: ie_p(3), il(3,2)
      integer(kind = kint) :: i1, i2, k1
      real(kind = kreal) :: dy1, dy2
!
!
      iflag = 0
!
      ie_p(1:3) = psf_case_tbl(mark_psf)%iedge(i_1ele,1:3)
      il(1:3,1) = node_on_edge_l(1,ie_p(1:3))
      il(1:3,2) = node_on_edge_l(2,ie_p(1:3))
!
      do k1 = 1, 3
        i1 = il(k1,1)
        i2 = il(k1,2)
        dy1 = xx_1ele(2,i1)
        dy2 = xx_1ele(2,i2)
        x_patch(1:3,k1) = (dy2*xx_1ele(1:3,i1) - dy1*xx_1ele(1:3,i2))   &
     &                     / (dy1-dy2)
      end do
!
      if     (x_patch(1,1).gt.zero .or. x_patch(1,2).gt.zero            &
     &   .or. x_patch(1,3).gt.zero) iflag = 1
!
      end subroutine set_triangle_at_meridian
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine init_ele_grp_med_patch(file_head, num_mat, num_mat_bc)
!
      use calypso_mpi
      use set_parallel_file_name
!
      character(len=kchara), intent(in) :: file_head
      integer(kind = kint), intent(in) :: num_mat, num_mat_bc
!
      character(len=kchara) :: file_name
!
!
      allocate(mark_elegrp(num_mat_bc))
!
      allocate (sta1(MPI_STATUS_SIZE,ione))
      allocate (req1(ione))
      allocate (sta2(MPI_STATUS_SIZE,nprocs))
      allocate (req2(nprocs))
!
      if(my_rank .eq. 0) then
        allocate(npatch_l(nprocs))
        allocate(istack_npatch_l(0:nprocs))
        istack_npatch_l = 0

        call add_dat_extension(file_head, file_name)
        write(*,*) 'ascii mesh file: ', trim(file_name)
        open (id_file, file = file_name, form = 'formatted')
!
        write(id_file,'(a)') '#'
        write(id_file,'(a)') '# number of element groups'
        write(id_file,'(a)') '#'
        write(id_file,'(i16)') num_mat
      end if
!
      end subroutine init_ele_grp_med_patch
!
!  ---------------------------------------------------------------------
!
      subroutine clean_ele_grp_med_patch
!
      use calypso_mpi
!
!
      deallocate(mark_elegrp)
      deallocate (sta1, req1, sta2, req2)
!
      if(my_rank .eq. 0) then
        close(id_file)
        deallocate(istack_npatch_l, npatch_l)
      end if
!
      end subroutine clean_ele_grp_med_patch
!
!  ---------------------------------------------------------------------
!
      subroutine collect_ele_grp_patch(grp_name)
!
      use calypso_mpi
!
      character(len=kchara), intent(in) :: grp_name
      integer(kind = kint) :: np_recv, nneib_recv
      integer(kind = kint) :: ip, num, inum, ist, k1
!
      real(kind = kreal), allocatable :: xyz_med_g(:,:)
!
!
      call MPI_Gather(npatch_grp, ione, CALYPSO_INTEGER, npatch_l,      &
     &    ione, CALYPSO_INTEGER, izero, CALYPSO_COMM, ierr_MPI)
!
      np_recv = 0
      if(my_rank .eq. 0) then
        do ip = 1, nprocs
          istack_npatch_l(ip) = istack_npatch_l(ip-1) + npatch_l(ip)
          if(np_recv .gt. 0) np_recv = np_recv + 1
        end do
        npatch_g = istack_npatch_l(nprocs)
!
        allocate(xyz_med_g(3,3*npatch_g))
      end if
!
      nneib_recv = 0
      call MPI_Isend(xyz_med(1,1), 9*npatch_grp, CALYPSO_REAL,          &
     &      izero, 0, CALYPSO_COMM, req1, ierr_MPI)
!
      if(my_rank .eq. 0) then
        nneib_recv = nprocs
        do ip = 1, nprocs
          ist = 3*istack_npatch_l(ip-1) + 1
          num = 9*npatch_l(ip)
          call MPI_Irecv(xyz_med_g(1,ist), num, CALYPSO_REAL,           &
     &        (ip-1), 0, CALYPSO_COMM, req2(ip), ierr_MPI)
        end do
      end if
      call MPI_WAITALL(nneib_recv, req2, sta2, ierr_MPI)
      call MPI_WAITALL(ione, req1, sta1, ierr_MPI)
!
      if(my_rank .eq. 0) then
          write(id_file,'(a)') '#'
          write(id_file,'(a)') '# element group name, node lists'
          write(id_file,'(a)') '#'
          write(id_file,'(a)') trim(grp_name)
          write(id_file,'(2i16)') (3*npatch_g), npatch_g
!
          do inum = 1, 3*npatch_g
            write(id_file,'(i16,1p3e23.13e3)')                          &
     &                               inum, xyz_med_g(1:3,inum)
          end do
          write(id_file,'(a)') '# element lists'
          do inum = 1, npatch_g
            write(id_file,'(4i16)') inum, ((3*inum+k1-3), k1=1,3)
          end do
!
          deallocate(xyz_med_g)
      end if
!
      end subroutine collect_ele_grp_patch
!
!  ---------------------------------------------------------------------
!
      end module set_med_patch_4_element_grp
