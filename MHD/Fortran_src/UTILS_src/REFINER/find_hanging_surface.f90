!find_hanging_surface.f90
!      module find_hanging_surface
!
      module find_hanging_surface
!
      use m_precision
!
      implicit none
!
      integer(kind = kint), allocatable :: iflag_hang_sf(:)
      integer(kind = kint), allocatable :: iflag_hang_ed(:)
!
      integer(kind = kint) :: nnod_hang_4
      integer(kind = kint) :: nnod_hang_2
      integer(kind = kint), allocatable :: inod_hang_4(:,:)
      integer(kind = kint), allocatable :: inod_hang_2(:,:)
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine allocate_iflag_hangings
!
      use m_geometry_parameter
!
      allocate( iflag_hang_sf(numsurf) )
      allocate( iflag_hang_ed(numedge) )
      iflag_hang_sf = 0
      iflag_hang_ed = 0
!
      end subroutine allocate_iflag_hangings
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine deallocate_iflag_hangings
!
!
      deallocate( iflag_hang_sf, iflag_hang_ed)
!
      end subroutine deallocate_iflag_hangings
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine check_hanging_surface
!
      use m_geometry_constants
      use m_geometry_parameter
      use m_geometry_data
      use m_refine_flag_parameters
      use m_refined_element_data
!
      integer(kind = kint) :: iele, isurf, iedge, iflag1, iflag2
      integer(kind = kint) :: iele1,  iele2, k1, k2
!
!
      call allocate_iflag_hangings
!
      do iele = 1, numele
        do k1 = 1, nsurf_4_ele
          isurf = abs( isf_4_ele(iele,k1) )
          if(iflag_refine_sf_lcl(k1,iele) .eq. iflag_dbl_sf             &
     &      .and.  iflag_refine_surf(isurf) .eq. iflag_nothing_sf) then
            iflag_refine_surf(isurf) = iflag_refine_sf_lcl(k1,iele)
          end if
        end do
      end do
!
      do iele = 1, numele
        do k1 = 1, nedge_4_ele
          iedge = abs(iedge_4_ele(iele,k1))
          if(iflag_refine_ed_lcl(k1,iele) .eq. iflag_dbl_ed             &
     &      .and.  iflag_refine_edge(iedge) .eq. iflag_nothing_ed) then
            iflag_refine_edge(iedge) = iflag_dbl_ed
          end if
        end do
      end do
!
!
      do isurf = 1, numsurf
        if(iele_4_surf(isurf,2,1) .gt. 0) then
          iele1 = iele_4_surf(isurf,1,1)
          k1 =    iele_4_surf(isurf,1,2)
          iele2 = iele_4_surf(isurf,2,1)
          k2 =    iele_4_surf(isurf,2,2)
          iflag1 = iflag_refine_sf_lcl(k1,iele1)
          iflag2 = iflag_refine_sf_lcl(k2,iele2)
          if(     (iflag1.eq.iflag_dbl_sf                               &
     &       .and. iflag2.eq.iflag_nothing_sf)                          &
     &       .or. (iflag1.eq.iflag_nothing_sf                           &
     &       .and. iflag2.eq.iflag_dbl_sf) ) then
              write(*,*) 'full hanginig surface', iele1, k1, isurf
              iflag_hang_sf(isurf) = 1
          end if
        end if
      end do
!
      do iele = 1, numele
        do k1 = 1, nedge_4_ele
          iedge = abs(iedge_4_ele(iele,k1))
          if(iflag_refine_edge(iedge) .eq. iflag_dbl_sf                 &
     &       .and. iflag_refine_ed_lcl(k1,iele).eq.iflag_nothing_ed)    &
     &       then
              write(*,*) 'hanginig edge', iedge, iele, k1
              iflag_hang_ed(iedge) = 1
          end if
        end do
      end do
!
      end subroutine check_hanging_surface
!
! -----------------------------------------------------------------------
!
      subroutine set_hanging_nodes
!
      use m_geometry_constants
      use m_geometry_parameter
      use m_geometry_data
      use m_refine_flag_parameters
      use m_refined_element_data
      use m_refined_node_id
!
      integer(kind = kint) :: isurf, iedge
      integer(kind = kint) :: icou, ist, k1
!
!
      nnod_hang_4 = 0
      do isurf = 1, numsurf
        nnod_hang_4 = nnod_hang_4 + iflag_hang_sf(isurf)
      end do
      nnod_hang_2 = 0
      do iedge = 1, numedge
        nnod_hang_2 = nnod_hang_2 + iflag_hang_ed(iedge)
      end do
!
      allocate(inod_hang_4(5,nnod_hang_4))
      allocate(inod_hang_2(3,nnod_hang_2))
!
      icou = 0
      do isurf = 1, numsurf
        if(iflag_hang_sf(isurf) .eq. 1) then
          icou = icou + 1
          ist = istack_nod_refine_surf(isurf-1)
          inod_hang_4(1,icou) = inod_refine_surf(ist+1)
          do k1 = 1, 4
            inod_hang_4(k1+1,icou) = ie_surf(isurf,k1)
          end do
        end if
      end do
!
!
      icou = 0
      do iedge = 1, numedge
        if(iflag_hang_ed(iedge) .eq. 1) then
          icou = icou + 1
          ist = istack_nod_refine_edge(iedge-1)
          inod_hang_2(1,icou) = inod_refine_edge(ist+1)
          do k1 = 1, 2
            inod_hang_2(k1+1,icou) = ie_edge(iedge,k1)
          end do
        end if
      end do
!
!
      write(*,*) 'inod_hang_4(1:5,icou)', nnod_hang_4
      do icou = 1, nnod_hang_4
        write(*,*) inod_hang_4(1:5,icou)
      end do
      write(*,*) 'inod_hang_2(1:3,icou)', nnod_hang_2
      do icou = 1, nnod_hang_2
        write(*,*) inod_hang_2(1:3,icou)
      end do
!
      end subroutine set_hanging_nodes
!
! -----------------------------------------------------------------------
!
      subroutine add_hanging_node_group_num
!
      use m_node_group
      use m_2nd_group_data
!
!
      if(nnod_hang_4 .gt. 0) num_bc_2nd = num_bc_2nd + 2
      if(nnod_hang_2 .gt. 0) num_bc_2nd = num_bc_2nd + 2
!
      end subroutine add_hanging_node_group_num
!
! -----------------------------------------------------------------------
!
      subroutine add_hanging_node_group_name
!
      use m_node_group
      use m_2nd_group_data
!
      integer(kind = kint) :: icou
!
      icou = num_bc
      if(nnod_hang_4 .gt. 0) then
        icou = icou+1
        bc_name_2nd(icou) = 'HANGING_NODE_SURF'
        bc_istack_2nd(icou) = bc_istack_2nd(icou-1) + nnod_hang_4
        icou = icou+1
        bc_name_2nd(icou) = 'HANGING_SOURCE_SURF'
        bc_istack_2nd(icou) = bc_istack_2nd(icou-1) + 4 * nnod_hang_4
!
        num_nod_bc_2nd = bc_istack_2nd(icou)
      end if
!
      if(nnod_hang_2 .gt. 0) then
        icou = icou+1
        bc_name_2nd(icou) = 'HANGING_NODE_EDGE'
        bc_istack_2nd(icou) = bc_istack_2nd(icou-1) + nnod_hang_2
        icou = icou+1
        bc_name_2nd(icou) = 'HANGING_SOURCE_EDGE'
        bc_istack_2nd(icou) = bc_istack_2nd(icou-1) + 2 * nnod_hang_2
!
        num_nod_bc_2nd = bc_istack_2nd(icou)
      end if
!
      end subroutine add_hanging_node_group_name
!
! -----------------------------------------------------------------------
!
      subroutine add_hanging_node_group_item
!
      use m_node_group
      use m_2nd_group_data
!
      integer(kind = kint) :: icou, inum, ist
!
      icou = num_bc
      if(nnod_hang_4 .gt. 0) then
        ist =  bc_istack_2nd(icou)
        do inum = 1, nnod_hang_4
          bc_item_2nd(ist+inum) = inod_hang_4(1,inum)
        end do
        icou = icou + 1
!
        ist =  bc_istack_2nd(icou)
        do inum = 1, nnod_hang_4
          bc_item_2nd(ist+4*inum-3) = inod_hang_4(2,inum)
          bc_item_2nd(ist+4*inum-2) = inod_hang_4(3,inum)
          bc_item_2nd(ist+4*inum-1) = inod_hang_4(4,inum)
          bc_item_2nd(ist+4*inum  ) = inod_hang_4(5,inum)
        end do
        icou = icou + 1
      end if
!
      if(nnod_hang_2 .gt. 0) then
        ist =  bc_istack_2nd(icou)
        do inum = 1, nnod_hang_2
          bc_item_2nd(ist+inum) = inod_hang_2(1,inum)
        end do
        icou = icou + 1
!
        ist =  bc_istack_2nd(icou)
        do inum = 1, nnod_hang_2
          bc_item_2nd(ist+2*inum-1) = inod_hang_2(2,inum)
          bc_item_2nd(ist+2*inum  ) = inod_hang_2(3,inum)
        end do
        icou = icou + 1
      end if
!
      end subroutine add_hanging_node_group_item
!
! -----------------------------------------------------------------------
!
      end module find_hanging_surface
