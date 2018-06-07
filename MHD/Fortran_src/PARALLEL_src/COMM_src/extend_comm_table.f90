!> @file  extend_comm_table.f90
!!      module extend_comm_table
!!
!! @author  H. Matsui
!! @date Programmed in Jan., 2009
!
!> @brief Mark node and element to extend export table
!!
!!@verbatim
!!@endverbatim
!
      module extend_comm_table
!
      use m_precision
      use m_constants
      use m_phys_constants
      use m_machine_parameter
!
      use t_comm_table
!
      implicit none
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine extend_node_comm_table(nod_comm, org_node, neib_nod,   &
     &          new_comm, new_node)
!
      use t_geometry_data
      use t_next_node_ele_4_node
      use m_merged_ucd_data
      use solver_SR_type
      use extend_comm_table_SR
      use mark_export_nod_ele_extend
!
      type(communication_table), intent(in) :: nod_comm
      type(node_data), intent(in) :: org_node
      type(next_nod_id_4_nod), intent(in) :: neib_nod
!
      type(communication_table), intent(inout) :: new_comm
      type(node_data), intent(inout) :: new_node
!
!>      Structure of double numbering
      type(parallel_double_numbering) :: dbl_id1
      type(parallel_double_numbering) :: dbl_id2
!
      integer(kind = kint), allocatable :: num_send_added(:)
      integer(kind = kint), allocatable :: istack_send_added(:)
!
      integer(kind = kint) :: ntot_send_added
      integer(kind = kint), allocatable :: inod_send_added(:)
      integer(kind = kint), allocatable :: irank_send_added(:)
      integer(kind = kint), allocatable :: iflag_send_2_del(:)
      integer(kind = kint_gl), allocatable :: inod_gl_send_added(:)
      real(kind = kreal), allocatable :: xx_send_added(:,:)
!
      integer(kind = kint), allocatable :: num_recv_added(:)
      integer(kind = kint), allocatable :: istack_recv_added(:)
!
      integer(kind = kint), allocatable :: ntot_recv_added
      integer(kind = kint), allocatable :: inod_recv_added(:)
      integer(kind = kint), allocatable :: irank_recv_added(:)
      integer(kind = kint), allocatable :: inod_new_id(:)
      integer(kind = kint_gl), allocatable :: inod_gl_recv_added(:)
      real(kind = kreal), allocatable :: xx_recv_added(:,:)
!
      integer(kind = kint), allocatable :: iflag_recv(:)
      integer(kind = kint), allocatable :: iflag_send(:)
      integer(kind = kint), allocatable :: iflag_node(:)
!
      integer(kind = kint), allocatable :: inod_import_new(:)
      integer(kind = kint), allocatable :: irank_import_new(:)
      integer(kind = kint), allocatable :: inod_export_new(:)
      integer(kind = kint), allocatable :: irank_export_new(:)
!
      integer(kind = kint), allocatable :: inod_lc_check(:)
      integer(kind = kint), allocatable :: irank_lc_check(:)
!
      integer(kind = kint) :: inum, inod, i, ist, ied, icou, ip, num
      integer(kind = kint) :: jnum, jnod, j, jst, jed
!
!
      call alloc_double_numbering(org_node%numnod, dbl_id1)
      call set_para_double_numbering                                    &
     &   (org_node%internal_node, nod_comm, dbl_id1)
      call calypso_mpi_barrier
!
      allocate(iflag_node(org_node%numnod))
      iflag_node(1:org_node%numnod) = 0
!
      allocate(num_send_added(nod_comm%num_neib))
      allocate(istack_send_added(0:nod_comm%num_neib))
      if(nod_comm%num_neib .gt. 0) num_send_added = 0
      istack_send_added = 0
!
      do i = 1, nod_comm%num_neib
        call mark_next_node_of_export(i, nod_comm%num_neib,             &
     &      nod_comm%istack_import, nod_comm%item_import,               &
     &      nod_comm%istack_export, nod_comm%item_export,               &
     &      org_node%numnod, neib_nod%ntot, neib_nod%istack_next,       &
     &      neib_nod%inod_next, iflag_node)
!
        do inod = 1, org_node%numnod
          num_send_added(i) = num_send_added(i) + iflag_node(inod)
        end do
        istack_send_added(i) = istack_send_added(i-1)                   &
     &                          + num_send_added(i)
      end do
!
      ntot_send_added = istack_send_added(nod_comm%num_neib)
!
      allocate(inod_send_added(ntot_send_added))
      allocate(xx_send_added(ntot_send_added,3))
      allocate(irank_send_added(ntot_send_added))
      allocate(iflag_send_2_del(ntot_send_added))
      allocate(inod_gl_send_added(ntot_send_added))
!
!$omp parallel workshare
      xx_send_added(1:ntot_send_added,1) = 0.0d0
      xx_send_added(1:ntot_send_added,2) = 0.0d0
      xx_send_added(1:ntot_send_added,3) = 0.0d0
      inod_send_added(1:ntot_send_added) =     0
      irank_send_added(1:ntot_send_added) =    0
      inod_gl_send_added(1:ntot_send_added) =  0
      iflag_send_2_del(1:ntot_send_added) =    0
!$omp end parallel workshare
!
      do i = 1, nod_comm%num_neib
        call mark_next_node_of_export(i, nod_comm%num_neib,             &
     &      nod_comm%istack_import, nod_comm%item_import,               &
     &      nod_comm%istack_export, nod_comm%item_export,               &
     &      org_node%numnod, neib_nod%ntot, neib_nod%istack_next,       &
     &      neib_nod%inod_next, iflag_node)
!
        icou = istack_send_added(i-1)
        do inod = 1, org_node%numnod
          if(iflag_node(inod) .gt. 0) then
            icou = icou + 1
            inod_send_added(icou) =    dbl_id1%inod_local(inod)
            irank_send_added(icou) =   dbl_id1%irank_home(inod)
            inod_gl_send_added(icou) = org_node%inod_global(inod)
            xx_send_added(icou,1) =    org_node%xx(inod,1)
            xx_send_added(icou,2) =    org_node%xx(inod,2)
            xx_send_added(icou,3) =    org_node%xx(inod,3)
          end if
        end do
      end do
!
!
      allocate(num_recv_added(nod_comm%num_neib))
      allocate(istack_recv_added(0:nod_comm%num_neib))
      if(nod_comm%num_neib .gt. 0) num_recv_added = 0
      istack_recv_added = 0
!
      call SOLVER_SEND_RECV_num_type                                    &
     &   (nod_comm, num_send_added, num_recv_added)
!
      do i = 1, nod_comm%num_neib
        istack_recv_added(i) = istack_recv_added(i-1)                   &
     &                          + num_recv_added(i)
      end do
      ntot_recv_added = istack_recv_added(nod_comm%num_neib)
!
      write(*,*) 'istack_send_added', my_rank, istack_send_added
      write(*,*) 'ntot_send_added', my_rank, ntot_send_added
      write(*,*) 'istack_recv_added', my_rank, istack_recv_added
      write(*,*) 'ntot_recv_added', my_rank, ntot_recv_added
!
      allocate(inod_recv_added(ntot_recv_added))
      allocate(irank_recv_added(ntot_recv_added))
      allocate(inod_new_id(ntot_recv_added))
      allocate(inod_gl_recv_added(ntot_recv_added))
      allocate(xx_recv_added(ntot_recv_added,3))
!
!$omp parallel workshare
      xx_recv_added(1:ntot_recv_added,1) = 0.0d0
      xx_recv_added(1:ntot_recv_added,2) = 0.0d0
      xx_recv_added(1:ntot_recv_added,3) = 0.0d0
      inod_recv_added(1:ntot_recv_added) =     0
      irank_recv_added(1:ntot_recv_added) =    0
      inod_gl_recv_added(1:ntot_recv_added) =  0
      inod_new_id(1:ntot_recv_added) =    0
!$omp end parallel workshare
!
      call added_geometry_send_recv(nod_comm,                           &
     &    istack_send_added, ntot_send_added, xx_send_added,            &
     &    istack_recv_added, ntot_recv_added, xx_recv_added)
      call added_global_id_send_recv(nod_comm,                          &
     &    istack_send_added, ntot_send_added, inod_gl_send_added,       &
     &    istack_recv_added, ntot_recv_added, inod_gl_recv_added)
      call added_nod_id_send_recv(nod_comm,                             &
     &    istack_send_added, ntot_send_added, inod_send_added,          &
     &    istack_recv_added, ntot_recv_added, inod_recv_added)
      call added_nod_id_send_recv(nod_comm,                             &
     &    istack_send_added, ntot_send_added, irank_send_added,         &
     &    istack_recv_added, ntot_recv_added, irank_recv_added)
!
!
      do i = 1, nod_comm%num_neib
        ist = istack_recv_added(i-1) + 1
        ied = istack_recv_added(i)
!
        do inum = ist, ied
          if(inod_new_id(inum) .lt. 0) cycle
          if(irank_recv_added(inum) .eq. nod_comm%id_neib(i)) cycle
!
          do j = 1, nod_comm%num_neib
            if(i .eq. j) cycle
!
            jst = istack_recv_added(j-1) + 1
            jed = istack_recv_added(j)
            do jnum = jst, jed
              if(irank_recv_added(jnum) .eq. irank_recv_added(inum)     &
     &         .and. inod_recv_added(jnum).eq.inod_recv_added(inum)     &
     &         .and. inod_new_id(jnum).eq.0) then
                inod_new_id(inum) = -1
                exit
              end if
            end do
            if(inod_new_id(inum) .lt. 0) exit
          end do
        end do
      end do
!
      do i = 1, nod_comm%num_neib
        ist = istack_recv_added(i-1) + 1
        ied = istack_recv_added(i)
!
        do inum = ist, ied
          if(inod_new_id(inum) .lt. 0) cycle
          if(irank_recv_added(inum) .eq. nod_comm%id_neib(i)) cycle
!
          do jnum = 1, nod_comm%istack_import(nod_comm%num_neib)
            jnod = nod_comm%item_import(jnum)
            if(dbl_id1%irank_home(jnod) .eq. irank_recv_added(inum)     &
     &       .and. dbl_id1%inod_local(jnod).eq.inod_recv_added(inum)    &
     &       ) then
              inod_new_id(inum) = -1
              exit
            end if
          end do
        end do
      end do
!
!
!      do i = 1, nod_comm%num_neib
!        ist = nod_comm%istack_import(i-1) + 1
!        ied = nod_comm%istack_import(i)
!        write(120+my_rank,*) 'import', nod_comm%id_neib(i), ist, ied
!
!        do inum = ist, ied
!          inod = nod_comm%item_import(inum)
!          write(120+my_rank,*) inum, inod,        &
!     &        dbl_id1%irank_home(inod), dbl_id1%inod_local(inod), '  '
!        end do
!      end do
!      do i = 1, nod_comm%num_neib
!        ist = istack_recv_added(i-1) + 1
!        ied = istack_recv_added(i)
!        write(120+my_rank,*) 'istack_recv_added',   &
!     &                        nod_comm%id_neib(i), ist, ied
!
!        do inum = ist, ied
!          write(120+my_rank,*) inum, irank_recv_added(inum),   &
!     &           inod_recv_added(inum), inod_new_id(inum) 
!        end do
!      end do
!
!
      call added_nod_id_send_recv(nod_comm,                             &
     &    istack_recv_added, ntot_recv_added, inod_new_id,              &
     &    istack_send_added, ntot_send_added, iflag_send_2_del)
!
!      do inum = 1, ntot_recv_added
!        if(inod_new_id(inum) .lt. 0) write(*,*) 'recv delete', my_rank,&
!     &      inum, irank_recv_added(inum), inod_recv_added(inum)
!      end do
!      do inum = 1, ntot_send_added
!        if(iflag_send_2_del(inum) .lt. 0) write(*,*) 'send delete',    &
!     &    my_rank, inum, irank_send_added(inum), inod_send_added(inum)
!      end do
!
!
      new_node%internal_node = org_node%internal_node
      new_node%numnod = org_node%numnod
      do inum = 1, ntot_recv_added
        if(inod_new_id(inum) .eq. 0) then
          new_node%numnod = new_node%numnod + 1
        end if
      end do
!
      call alloc_node_geometry_base(new_node)
      call alloc_double_numbering(new_node%numnod, dbl_id2)
!
!$omp parallel do
      do inod = 1, org_node%numnod
        new_node%inod_global(inod) = org_node%inod_global(inod)
        new_node%xx(inod,1) = org_node%xx(inod,1)
        new_node%xx(inod,2) = org_node%xx(inod,2)
        new_node%xx(inod,3) = org_node%xx(inod,3)
        dbl_id2%inod_local(inod) = dbl_id1%inod_local(inod)
        dbl_id2%irank_home(inod) = dbl_id1%irank_home(inod)
      end do
!$omp end parallel do
!
      icou = org_node%numnod
      do inum = 1, ntot_recv_added
        if(inod_new_id(inum) .eq. 0) then
          icou = icou + 1
          inod_new_id(inum) = icou
          new_node%inod_global(icou) = inod_gl_recv_added(inum)
          new_node%xx(icou,1) = xx_recv_added(inum,1)
          new_node%xx(icou,2) = xx_recv_added(inum,2)
          new_node%xx(icou,3) = xx_recv_added(inum,3)
          dbl_id2%inod_local(icou) = inod_recv_added(inum)
          dbl_id2%irank_home(icou) = irank_recv_added(inum)
        end if
      end do
!
!      write(100+my_rank,*) new_node%numnod, &
!     &             new_node%internal_node, org_node%numnod
!      do inod = 1, new_node%numnod
!        write(100+my_rank,*) inod, dbl_id2%inod_local(inod),  &
!     &         dbl_id2%irank_home(inod), new_node%inod_global(inod) 
!      end do
!
      allocate(iflag_recv(0:nprocs-1))
      allocate(iflag_send(0:nprocs-1))
      iflag_recv(0:nprocs-1) = 0
      iflag_send(0:nprocs-1) = 0
!
      do inum = 1, ntot_recv_added
        ip = irank_recv_added(inum)
        iflag_recv(ip) = 1
      end do
!
      do i = 1, nod_comm%num_neib
        ip = nod_comm%id_neib(i)
        iflag_recv(ip) = -1
      end do
!
      do ip = 0, nprocs-1
        call MPI_Scatter(iflag_recv(0), ione, CALYPSO_INTEGER,          &
     &                   iflag_send(ip), ione, CALYPSO_INTEGER,         &
     &                   ip, CALYPSO_COMM, ierr_MPI)
      end do
!

      new_comm%num_neib = 0
      do ip = 0, nprocs-1
        if(iflag_recv(ip).ne.0 .or. iflag_send(ip).ne.0) then
          new_comm%num_neib = new_comm%num_neib + 1
        end if
      end do
!
!
      call allocate_type_comm_tbl_num(new_comm)
!
      new_comm%id_neib(1:nod_comm%num_neib)                             &
     &              = nod_comm%id_neib(1:nod_comm%num_neib)
      icou = nod_comm%num_neib
      do i = 0, nprocs-1
        ip = mod(i+my_rank,nprocs)
        if(iflag_recv(ip).gt.0 .or. iflag_send(ip).gt.0) then
          icou = icou + 1
          new_comm%id_neib(i) = ip
        end if
      end do
!
      do i = 1, nod_comm%num_neib
        new_comm%num_import(i)                                          &
     &       = nod_comm%istack_import(i) - nod_comm%istack_import(i-1)
      end do
      do i = nod_comm%num_neib+1, new_comm%num_neib
        new_comm%num_import(i) = 0
      end do
      do i = 1, new_comm%num_neib
        ip = new_comm%id_neib(i)
        do inum = 1, ntot_recv_added
          if(irank_recv_added(inum).eq.ip                               &
     &         .and. inod_new_id(inum).gt. 0) then
            new_comm%num_import(i) = new_comm%num_import(i) + 1
          end if
        end do
      end do
!
      call SOLVER_SEND_RECV_num_type                                    &
     &   (new_comm, new_comm%num_import, new_comm%num_export)
!
      do i = 1, new_comm%num_neib
        new_comm%istack_import(i) = new_comm%istack_import(i-1)         &
     &                             + new_comm%num_import(i)
        new_comm%istack_export(i) = new_comm%istack_export(i-1)         &
     &                             + new_comm%num_export(i)
      end do
      new_comm%ntot_import = new_comm%istack_import(new_comm%num_neib)
      new_comm%ntot_export = new_comm%istack_export(new_comm%num_neib)
!
      call allocate_type_comm_tbl_item(new_comm)
!
      do i = 1, nod_comm%num_neib
        ist = new_comm%istack_import(i-1)
        jst = nod_comm%istack_import(i-1)
        num = nod_comm%istack_import(i) - nod_comm%istack_import(i-1)
        do inum = 1, num
          new_comm%item_import(ist+inum)                                &
     &       = nod_comm%item_import(jst+inum)
        end do
!
        ist = new_comm%istack_export(i-1)
        jst = nod_comm%istack_export(i-1)
        num = nod_comm%istack_export(i) - nod_comm%istack_export(i-1)
        do inum = 1, num
          new_comm%item_export(ist+inum)                                &
     &       = nod_comm%item_export(jst+inum)
        end do
      end do
!
      do i = 1, new_comm%num_neib
        ip = new_comm%id_neib(i)
        icou = new_comm%istack_import(i-1)                              &
     &        + nod_comm%istack_import(i) - nod_comm%istack_import(i-1)
        do inum = 1, ntot_recv_added
          if(irank_recv_added(inum).eq.ip                               &
     &         .and. inod_new_id(inum).gt.0) then
            icou = icou + 1
            new_comm%item_import(icou) = inod_new_id(inum)
          end if
        end do
      end do
!
      allocate(inod_import_new(new_comm%ntot_import))
      allocate(irank_import_new(new_comm%ntot_import))
      inod_import_new = 0
      irank_import_new = -1
!
      do inum = 1, new_comm%ntot_import
        inod =  new_comm%item_import(inum)
        inod_import_new(inum) =  dbl_id2%inod_local(inod)
        irank_import_new(inum) = dbl_id2%irank_home(inod)
      end do
!
      allocate(inod_export_new(new_comm%ntot_export))
      allocate(irank_export_new(new_comm%ntot_export))
      inod_export_new = 0
      irank_export_new = -1
!
      call added_nod_id_send_recv(new_comm,                             &
     &  new_comm%istack_import, new_comm%ntot_import, inod_import_new,  &
     &  new_comm%istack_export, new_comm%ntot_export, inod_export_new)
      call added_nod_id_send_recv(new_comm,                             &
     &  new_comm%istack_import, new_comm%ntot_import, irank_import_new, &
     &  new_comm%istack_export, new_comm%ntot_export, irank_export_new)
!
      do i = 1, new_comm%num_neib
        ip = new_comm%id_neib(i)
        ist = new_comm%istack_export(i-1) + 1
        ied = new_comm%istack_export(i)
        do inum = ist, ied
          if(irank_export_new(inum) .eq. my_rank) then
            new_comm%item_export(inum) = inod_export_new(inum)
          end if
        end do
      end do
!
!
      allocate(inod_lc_check(new_node%numnod))
      allocate(irank_lc_check(new_node%numnod))
      inod_lc_check =   0
      irank_lc_check = -1
!
      do inod = 1, new_node%internal_node
        inod_lc_check(inod) = inod
        irank_lc_check(inod) = my_rank
      end do
!
      call SOLVER_SEND_RECV_int_type                                    &
     &   (new_node%numnod, new_comm, inod_lc_check)
      call SOLVER_SEND_RECV_int_type                                    &
     &   (new_node%numnod, new_comm, irank_lc_check)
!
!      write(50+my_rank,*) 'error list'
!      do inod = new_node%internal_node+1, new_node%numnod
!        if(dbl_id2%irank_home(inod) .ne. irank_lc_check(inod)          &
!     &    .and. dbl_id2%inod_local(inod) .ne. inod_lc_check(inod)) then
!          write(50+my_rank,*) inod, my_rank,                           &
!     &     dbl_id2%irank_home(inod), irank_lc_check(inod),             &
!     &     dbl_id2%inod_local(inod), inod_lc_check(inod)
!        end if
!      end do
!      close(50+my_rank)
      call calypso_mpi_barrier
!
!      write(*,*) 'iflag_recv', my_rank, iflag_recv
!      write(*,*) 'iflag_send', my_rank, iflag_send
!      write(*,*) 'num_neib', my_rank,                                  &
!     &       nod_comm%num_neib, new_comm%num_neib
!
      end subroutine extend_node_comm_table
!
!  ---------------------------------------------------------------------
!
      subroutine extend_ele_comm_table                                  &
     &         (nod_comm, ele_comm, org_node, ele, neib_ele, neib_nod,  &
     &          new_comm, new_node, new_ele_comm, new_ele)
!
      use t_geometry_data
      use t_next_node_ele_4_node
      use m_merged_ucd_data
      use solver_SR_type
      use extend_comm_table_SR
      use mark_export_nod_ele_extend
      use const_mesh_information
      use const_element_comm_tables
!
      type(communication_table), intent(in) :: nod_comm
      type(communication_table), intent(in) :: ele_comm
      type(node_data), intent(in) :: org_node
      type(element_data), intent(in) :: ele
      type(element_around_node), intent(in) :: neib_ele
      type(next_nod_id_4_nod), intent(in) :: neib_nod
!
      type(communication_table), intent(in) :: new_comm
!
      type(node_data), intent(inout) :: new_node
      type(communication_table), intent(inout) :: new_ele_comm
      type(element_data), intent(inout) :: new_ele
!
!>      Structure of double numbering
      type(parallel_double_numbering) :: dbl_id1
      type(parallel_double_numbering) :: dbl_id2
      type(parallel_double_numbering) :: dbl_ele
!
      type(belonged_table) :: new_blng_tbl
!
      integer(kind = kint), allocatable :: num_send_added(:)
      integer(kind = kint), allocatable :: istack_send_added(:)
!
      integer(kind = kint) :: ntot_send_added
      integer(kind = kint), allocatable :: iele_lc_added(:)
      integer(kind = kint), allocatable :: iele_send_added(:)
      integer(kind = kint), allocatable :: irank_send_added(:)
      integer(kind = kint), allocatable :: iflag_send_2_del(:)
      integer(kind = kint_gl), allocatable :: iele_gl_send_added(:)
      integer(kind = kint), allocatable :: ie_send_added(:,:)
      integer(kind = kint), allocatable :: ip_send_added(:,:)
!
      integer(kind = kint), allocatable :: num_recv_added(:)
      integer(kind = kint), allocatable :: istack_recv_added(:)
!
      integer(kind = kint), allocatable :: ntot_recv_added
      integer(kind = kint), allocatable :: iele_recv_added(:)
      integer(kind = kint), allocatable :: irank_recv_added(:)
      integer(kind = kint), allocatable :: iflag_recv_2_del(:)
      integer(kind = kint_gl), allocatable :: iele_gl_recv_added(:)
      integer(kind = kint), allocatable :: ie_recv_added(:,:)
      integer(kind = kint), allocatable :: ip_recv_added(:,:)
!
      integer(kind = kint), allocatable :: iflag_node(:)
      integer(kind = kint), allocatable :: iflag_ele(:)
!
      integer(kind = kint), allocatable :: istack_ele_ip(:)
      integer(kind = kint), allocatable :: iele_by_ip(:)
!
      integer(kind = kint), allocatable :: iflag_sleeve1(:)
      integer(kind = kint), allocatable :: iflag_sleeve2(:)
!
      integer(kind = kint) :: inum, inod, i, ist, ied, icou
      integer(kind = kint) :: jnum, jnod, j, jst, jed, jele
!
      integer(kind = kint) :: iele, k1, i1, iflag, ip
!
!
      call alloc_double_numbering(new_node%numnod, dbl_id1)
      call set_para_double_numbering                                    &
     &   (new_node%internal_node, nod_comm, dbl_id1)
      call alloc_double_numbering(new_node%numnod, dbl_id2)
      call set_para_double_numbering                                    &
     &   (new_node%internal_node, new_comm, dbl_id2)
      call alloc_double_numbering(ele%numele, dbl_ele)
      call set_para_ele_double_numbering                                &
     &   (new_node%internal_node, ele_comm, ele, dbl_ele)
      call calypso_mpi_barrier
!
      allocate(iflag_node(org_node%numnod))
      allocate(iflag_ele(ele%numele))
      iflag_node(1:org_node%numnod) = 0
      iflag_ele(1:ele%numele) = 0
!
      allocate(num_send_added(nod_comm%num_neib))
      allocate(istack_send_added(0:nod_comm%num_neib))
      if(nod_comm%num_neib .gt. 0) num_send_added = 0
      istack_send_added = 0
!
      allocate(iflag_sleeve1(nod_comm%num_neib))
      allocate(iflag_sleeve2(new_comm%num_neib))
      if(nod_comm%num_neib .gt. 0) num_send_added = 0
      istack_send_added = 0
!
      do i = 1, nod_comm%num_neib
        call mark_used_ele_of_export(i, nod_comm%num_neib,              &
     &      new_comm%istack_export, new_comm%item_export,               &
     &      nod_comm%istack_export, nod_comm%item_export,               &
     &      org_node%numnod, neib_ele%ntot, neib_ele%istack_4_node,     &
     &      neib_ele%iele_4_node, ele%numele, ele%nnod_4_ele, ele%ie,   &
     &      iflag_node, iflag_ele)
!
        do inod = 1, ele%numele
          num_send_added(i) = num_send_added(i) + iflag_ele(inod)
        end do
        istack_send_added(i) = istack_send_added(i-1)                   &
     &                          + num_send_added(i)
      end do
      ntot_send_added = istack_send_added(nod_comm%num_neib)
!
      write(*,*) 'istack_send_added ele', istack_send_added
!
      allocate(iele_lc_added(ntot_send_added))
      allocate(iele_send_added(ntot_send_added))
      allocate(ie_send_added(ntot_send_added,ele%nnod_4_ele))
      allocate(ip_send_added(ntot_send_added,ele%nnod_4_ele))
      allocate(irank_send_added(ntot_send_added))
      allocate(iflag_send_2_del(ntot_send_added))
      allocate(iele_gl_send_added(ntot_send_added))
!
!$omp parallel workshare
      ie_send_added(1:ntot_send_added,1:ele%nnod_4_ele) =    0
      ip_send_added(1:ntot_send_added,1:ele%nnod_4_ele) =    0
!$omp end parallel workshare
!$omp parallel workshare
      iele_lc_added(1:ntot_send_added) =       0
      iele_send_added(1:ntot_send_added) =     0
      irank_send_added(1:ntot_send_added) =    0
      iele_gl_send_added(1:ntot_send_added) =  0
      iflag_send_2_del(1:ntot_send_added) =    0
!$omp end parallel workshare
!
      do i = 1, nod_comm%num_neib
        call mark_used_ele_of_export(i, nod_comm%num_neib,              &
     &      new_comm%istack_export, new_comm%item_export,               &
     &      nod_comm%istack_export, nod_comm%item_export,               &
     &      org_node%numnod, neib_ele%ntot, neib_ele%istack_4_node,     &
     &      neib_ele%iele_4_node, ele%numele, ele%nnod_4_ele, ele%ie,   &
     &      iflag_node, iflag_ele)
!
        icou = istack_send_added(i-1)
        do iele = 1, ele%numele
          if(iflag_ele(iele) .gt. 0) then
            icou = icou + 1
            iele_lc_added(icou) = iele
            iele_send_added(icou) =    dbl_ele%inod_local(iele)
            irank_send_added(icou) =   dbl_ele%irank_home(iele)
            iele_gl_send_added(icou) = ele%iele_global(iele)
            do k1 = 1, ele%nnod_4_ele
              inod = ele%ie(iele,k1)
              ie_send_added(icou,k1) =    dbl_id1%inod_local(inod)
              ip_send_added(icou,k1) =    dbl_id1%irank_home(inod)
            end do
          end if
        end do
      end do
!
      do i = 1, nod_comm%num_neib
        ist = istack_send_added(i-1) + 1
        ied = istack_send_added(i)
        write(50+my_rank,*) 'istack_send_added',  &
     &                      i, nod_comm%id_neib(i), ist, ied
        do inum = ist, ied
          if(irank_send_added(inum) .eq. nod_comm%id_neib(i)) then
              write(50+my_rank,*) inum, iele_lc_added(inum),       &
     &         ie_send_added(inum,1:ele%nnod_4_ele)
              write(50+my_rank,*) inum, irank_send_added(inum),    &
     &         ip_send_added(inum,1:ele%nnod_4_ele)
          end if
        end do
      end do
!
      allocate(num_recv_added(nod_comm% num_neib))
      allocate(istack_recv_added(0:nod_comm%num_neib))
      if(nod_comm%num_neib .gt. 0) num_recv_added = 0
      istack_recv_added = 0
!
      call SOLVER_SEND_RECV_num_type                                    &
     &   (nod_comm, num_send_added, num_recv_added)
!
      do i = 1, nod_comm%num_neib
        istack_recv_added(i) = istack_recv_added(i-1)                   &
     &                          + num_recv_added(i)
      end do
      ntot_recv_added = istack_recv_added(nod_comm%num_neib)
!
      write(*,*) 'istack_send_added', my_rank, istack_send_added
      write(*,*) 'ntot_send_added', my_rank, ntot_send_added
      write(*,*) 'istack_recv_added', my_rank, istack_recv_added
      write(*,*) 'ntot_recv_added', my_rank, ntot_recv_added
!
      allocate(iele_recv_added(ntot_recv_added))
      allocate(irank_recv_added(ntot_recv_added))
      allocate(iflag_recv_2_del(ntot_recv_added))
      allocate(iele_gl_recv_added(ntot_recv_added))
      allocate(ie_recv_added(ntot_recv_added,ele%nnod_4_ele))
      allocate(ip_recv_added(ntot_recv_added,ele%nnod_4_ele))
!
!$omp parallel workshare
      iele_recv_added(1:ntot_recv_added) =     0
      irank_recv_added(1:ntot_recv_added) =    0
      iele_gl_recv_added(1:ntot_recv_added) =  0
      iflag_recv_2_del(1:ntot_recv_added) =    0
!$omp end parallel workshare
!$omp parallel workshare
      ie_recv_added(1:ntot_recv_added,1:ele%nnod_4_ele) =    0
      ip_recv_added(1:ntot_recv_added,1:ele%nnod_4_ele) =    0
!$omp end parallel workshare
!
      call added_global_id_send_recv(nod_comm,                          &
     &    istack_send_added, ntot_send_added, iele_gl_send_added,       &
     &    istack_recv_added, ntot_recv_added, iele_gl_recv_added)
      call added_nod_id_send_recv(nod_comm,                             &
     &    istack_send_added, ntot_send_added, iele_send_added,          &
     &    istack_recv_added, ntot_recv_added, iele_recv_added)
      call added_nod_id_send_recv(nod_comm,                             &
     &    istack_send_added, ntot_send_added, irank_send_added,         &
     &    istack_recv_added, ntot_recv_added, irank_recv_added)
      do k1 = 1, ele%nnod_4_ele
        call added_nod_id_send_recv(nod_comm,                           &
     &      istack_send_added, ntot_send_added, ie_send_added(1,k1),    &
     &      istack_recv_added, ntot_recv_added, ie_recv_added(1,k1))
        call added_nod_id_send_recv(nod_comm,                           &
     &      istack_send_added, ntot_send_added, ip_send_added(1,k1),    &
     &      istack_recv_added, ntot_recv_added, ip_recv_added(1,k1))
      end do
!
!      do i = 1, nod_comm%num_neib
!        ist = istack_recv_added(i-1) + 1
!        ied = istack_recv_added(i)
!        write(50+my_rank,*) 'istack_recv_added',  &
!     &                      i, nod_comm%id_neib(i), ist, ied
!        do inum = ist, ied
!          if(irank_send_added(inum) .eq. nod_comm%id_neib(i)) then
!              write(50+my_rank,*) inum,                            &
!     &         ie_recv_added(inum,1:ele%nnod_4_ele)
!              write(50+my_rank,*) inum,                            &
!     &         ip_recv_added(inum,1:ele%nnod_4_ele)
!          end if
!        end do
!      end do
!
      allocate(istack_ele_ip(0:nprocs))
      allocate(iele_by_ip(ele%numele))
      istack_ele_ip = 0
      iele_by_ip =    0
!
      icou = 0
      do ip = 1, nprocs
        do iele = 1, ele%numele
          if(dbl_ele%irank_home(iele) .eq. (ip-1)) then
            icou = icou + 1
            iele_by_ip(icou) = iele
          end if
        end do
        istack_ele_ip(ip) = icou
      end do
!
      do i = 1, nod_comm%num_neib
        ist = istack_recv_added(i-1) + 1
        ied = istack_recv_added(i)
        do iele = ist, ied
          ip = irank_recv_added(iele) + 1
          jst = istack_ele_ip(ip-1) + 1
          jed = istack_ele_ip(ip)
          do jnum = jst, jed
            jele = iele_by_ip(jnum)
            if(iele_recv_added(iele) .eq. dbl_ele%inod_local(jele))     &
     &       then
              iflag_recv_2_del(iele) = 1
              exit
            end if
          end do
!
        end do
      end do
!
      do i = 1, nod_comm%num_neib
        ist = istack_recv_added(i-1) + 1
        ied = istack_recv_added(i)
        do iele = ist, ied
          if(iflag_recv_2_del(iele) .gt. 0) cycle
!
          do jele = 1, ntot_recv_added
            if(iele .eq. jele) cycle
            if(irank_recv_added(jele) .eq. irank_recv_added(iele)       &
     &       .and. iele_recv_added(jele) .eq. iele_recv_added(iele)     &
     &       .and. iflag_recv_2_del(jele) .eq. 0) then
              iflag_recv_2_del(iele) = 1
              exit
            end if
          end do
        end do
      end do
!
      call added_nod_id_send_recv(nod_comm,                             &
     &    istack_recv_added, ntot_recv_added, iflag_recv_2_del,         &
     &    istack_send_added, ntot_send_added, iflag_send_2_del)
!
      new_ele%numele =     ele%numele
      new_ele%nnod_4_ele = ele%nnod_4_ele
      new_ele%internal_ele = ele%internal_ele
      do inum = 1, ntot_recv_added
        if(iflag_recv_2_del(inum) .eq. 0) then
          new_ele%numele = new_ele%numele + 1
        end if
      end do
!
      call allocate_ele_connect_type(new_ele)
!
!$omp parallel workshare
      new_ele%elmtyp(1:ele%numele) = ele%elmtyp(1:ele%numele)
      new_ele%nodelm(1:ele%numele) = ele%nodelm(1:ele%numele)
      new_ele%iele_global(1:ele%numele) = ele%iele_global(1:ele%numele)
!$omp end parallel workshare
!$omp parallel workshare
      new_ele%elmtyp(ele%numele+1:new_ele%numele) = ele%elmtyp(1)
      new_ele%nodelm(ele%numele+1:new_ele%numele) = ele%nodelm(1)
!$omp end parallel workshare
!
!$omp parallel
      do k1 = 1, ele%nnod_4_ele
!$omp workshare
        new_ele%ie(1:ele%numele,k1) = ele%ie(1:ele%numele,k1)
!$omp end workshare nowait
      end do
!$omp end parallel
!
      icou = ele%numele
      do i = 1, nod_comm%num_neib
        ist = istack_recv_added(i-1) + 1
        ied = istack_recv_added(i)
        do inum = ist, ied
          if(iflag_recv_2_del(inum) .gt. 0)  cycle
!
          icou = icou + 1
          new_ele%iele_global(icou) = iele_gl_recv_added(inum)
          do k1 = 1, ele%nnod_4_ele
            do jnod = new_node%numnod, 1, -1
              if(ip_recv_added(inum,k1) .eq. dbl_id2%irank_home(jnod)   &
     &          .and. ie_recv_added(inum,k1)                            &
     &              .eq. dbl_id2%inod_local(jnod)) then
                new_ele%ie(icou,k1) = jnod
                exit
              end if
            end do
          end do
        end do
      end do
!
!      write(50+my_rank,*) 'new_ele%iele', new_ele%numele
!      do i = 1, new_ele%numele
!         write(50+my_rank,*) new_ele%iele_global(i), new_ele%ie(i,1:8)
!      end do
!
      call calypso_mpi_barrier
!
      call allocate_sph_node_geometry(new_node)
      call set_nod_and_ele_infos(new_node, new_ele)
      call const_ele_comm_tbl                                           &
     &   (new_node, new_ele, new_comm, new_blng_tbl, new_ele_comm)
!
      end subroutine extend_ele_comm_table
!
!  ---------------------------------------------------------------------
!
      end module extend_comm_table

