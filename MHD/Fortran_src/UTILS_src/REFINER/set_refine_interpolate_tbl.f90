!set_refine_interpolate_tbl.f90
!      module set_refine_interpolate_tbl
!
!     Written by H. Matsui on Apr., 2010
!
!!      subroutine set_itp_course_to_fine_origin                        &
!!     &         (ele, surf, edge, itp_org)
!      subroutine set_itp_course_to_fine_dest(nnod_2, itp_dest)
!
!      subroutine set_itp_fine_to_course_origin(nnod_4_ele)
!      subroutine set_itp_fine_to_course_dest(itp_dest)
!
      module set_refine_interpolate_tbl
!
      use m_precision
      use m_machine_parameter
!
      use m_constants
      use m_geometry_constants
      use m_refined_node_id
!
      use t_interpolate_tbl_org
      use t_interpolate_tbl_dest
!
      use copy_local_position_2_ele
!
      implicit none
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine set_itp_course_to_fine_origin                          &
     &         (ele, surf, edge, itp_org)
!
      use t_geometry_data
      use t_surface_data
      use t_edge_data
!
      type(element_data), intent(in) :: ele
      type(surface_data), intent(in) :: surf
      type(edge_data), intent(in) :: edge
!
      type(interpolate_table_org), intent(inout) :: itp_org
!
      integer(kind = kint) :: iele, isurf, iedge, inod, isig
      integer(kind = kint) :: isf_ele, ied_ele, ist, ied, inum, k1
      real(kind = kreal) :: xi_surf(2), xi_ele(3)
!
!
      itp_org%iflag_self_itp_send = ione
      call set_num_dest_domain(ione, itp_org)
      call alloc_itp_num_org(np_smp, itp_org)
!
      itp_org%id_dest_domain(1) =             izero
      itp_org%istack_itp_type_org(0) = izero
      itp_org%istack_itp_type_org(1) =  ntot_nod_refine_nod
      itp_org%istack_itp_type_org(2) = itp_org%istack_itp_type_org(1)   &
     &                               + ntot_nod_refine_edge
      itp_org%istack_itp_type_org(3) = itp_org%istack_itp_type_org(2)   &
     &                               + ntot_nod_refine_surf
      itp_org%istack_itp_type_org(4) = itp_org%istack_itp_type_org(3)   &
     &                               + ntot_nod_refine_ele
!
      itp_org%istack_nod_tbl_org(0) = izero
      itp_org%istack_nod_tbl_org(1) = itp_org%istack_itp_type_org(4)
      itp_org%ntot_table_org =        itp_org%istack_nod_tbl_org(1)
!
      call alloc_itp_table_org(itp_org)
      itp_org%itype_inter_org(1:itp_org%ntot_table_org) = -1
!
      do iele = 1, ele%numele
!
        do k1 = 1, ele%nnod_4_ele
          inod = ele%ie(iele,k1)
!
          if(itp_org%itype_inter_org(inod) .eq. -1) then
            call copy_node_local_posi_2_element(k1, xi_ele)
!
            itp_org%inod_gl_dest_4_org(inod) = inod
            itp_org%iele_org_4_org(inod) =     iele
            itp_org%itype_inter_org(inod) =    k1
            itp_org%coef_inter_org(inod,1:3) = xi_ele(1:3)
          end if
        end do
!
        do ied_ele = 1, nedge_4_ele
          iedge = abs(edge%iedge_4_ele(iele,ied_ele))
          isig = edge%iedge_4_ele(iele,ied_ele) / iedge
!
            ist = istack_nod_refine_edge(iedge-1) + 1
            ied = istack_nod_refine_edge(iedge)
            do inum = ist, ied
              inod = inum + ntot_nod_refine_nod
!
              if(itp_org%itype_inter_org(inod) .eq. -1) then
                call copy_edge_local_posi_2_element(ied_ele, isig,      &
     &              xi_refine_edge(inum,1), xi_ele)
!
                itp_org%inod_gl_dest_4_org(inod) = inod
                itp_org%iele_org_4_org(inod) =     iele
                itp_org%itype_inter_org(inod) =    100 + ied_ele
                itp_org%coef_inter_org(inod,1:3) = xi_ele(1:3)
              end if
            end do
        end do
!
        do isf_ele = 1, nsurf_4_ele
          isurf = abs(surf%isf_4_ele(iele,isf_ele))
          if(surf%isf_4_ele(iele,isf_ele) .gt. 0) then
            ist = istack_nod_refine_surf(isurf-1) + 1
            ied = istack_nod_refine_surf(isurf)
            do inum = ist, ied
              inod = inum + ntot_nod_refine_nod + ntot_nod_refine_edge
!
              if(itp_org%itype_inter_org(inod) .eq. -1) then
                xi_surf(1) = xi_refine_surf(inum,1)
                xi_surf(2) = xi_refine_surf(inum,2)
                call copy_surf_local_posi_2_element(isf_ele,            &
     &              surf%isf_rot_ele(iele,isf_ele), xi_surf, xi_ele)
!
                itp_org%inod_gl_dest_4_org(inod) = inod
                itp_org%iele_org_4_org(inod) =     iele
                itp_org%itype_inter_org(inod) =    200 + isf_ele
                itp_org%coef_inter_org(inod,1:3) = xi_ele(1:3)
              end if
            end do
          end if
        end do
!
        ist = istack_nod_refine_ele(iele-1) + 1
        ied = istack_nod_refine_ele(iele)
        do inum = ist, ied
          inod = inum + ntot_nod_refine_nod + ntot_nod_refine_edge      &
     &                + ntot_nod_refine_surf
!
          itp_org%inod_gl_dest_4_org(inod) = inod
          itp_org%iele_org_4_org(inod) =     iele
          itp_org%itype_inter_org(inod) =    izero
          itp_org%coef_inter_org(inod,1:3) = xi_refine_ele(inum,1:3)
        end do
!
      end do
!
      end subroutine set_itp_course_to_fine_origin
!
! ----------------------------------------------------------------------
!
      subroutine set_itp_course_to_fine_dest(nnod_2, itp_dest)
!
      integer(kind = kint), intent(in) :: nnod_2
      type(interpolate_table_dest), intent(inout) :: itp_dest
!
      integer(kind = kint) :: inod
!
!
      itp_dest%iflag_self_itp_recv = ione
!
      call set_num_org_domain(ione, itp_dest)
      call alloc_itp_num_dest(itp_dest)
!
      itp_dest%id_org_domain(1) =       izero
      itp_dest%istack_nod_tbl_dest(0) = izero
      itp_dest%istack_nod_tbl_dest(1) = nnod_2
      itp_dest%ntot_table_dest =        nnod_2
!
      call alloc_itp_table_dest(itp_dest)
!
      do inod = 1, nnod_2
        itp_dest%inod_dest_4_dest(inod) = inod
      end do
!
      end subroutine set_itp_course_to_fine_dest
!
! ----------------------------------------------------------------------
!
      subroutine set_itp_fine_to_course_origin(nnod_4_ele, itp_org)
!
      use m_refined_element_data
!
      integer(kind = kint), intent(in) :: nnod_4_ele
      type(interpolate_table_org), intent(inout) :: itp_org
!
      integer(kind = kint) :: iele, k1, inod
      real(kind = kreal) :: xi_ele(3)
!
!
      itp_org%iflag_self_itp_send = ione
      call set_num_dest_domain(ione, itp_org)
      call alloc_itp_num_org(np_smp, itp_org)
!
      itp_org%id_dest_domain(1) = izero
      itp_org%istack_itp_type_org(0) = izero
      itp_org%istack_itp_type_org(1) = ntot_nod_refine_nod
      itp_org%istack_itp_type_org(2) = ntot_nod_refine_nod
      itp_org%istack_itp_type_org(3) = ntot_nod_refine_nod
      itp_org%istack_itp_type_org(4) = ntot_nod_refine_nod
      itp_org%istack_nod_tbl_org(0) = izero
      itp_org%istack_nod_tbl_org(1) = ntot_nod_refine_nod
      itp_org%ntot_table_org =        itp_org%istack_nod_tbl_org(1)
!
      call alloc_itp_table_org(itp_org)
      itp_org%itype_inter_org(1:itp_org%ntot_table_org) = -1
!
      do iele = 1, ntot_ele_refined
!
        do k1 = 1, nnod_4_ele
          inod = ie_refined(iele,k1)
!
          if(inod.le.ntot_nod_refine_nod) then
            if(itp_org%itype_inter_org(inod).eq.-1) then
              call copy_node_local_posi_2_element(k1, xi_ele)
!
              itp_org%inod_gl_dest_4_org(inod) = inod
              itp_org%iele_org_4_org(inod) =     iele
              itp_org%itype_inter_org(inod) =     k1
              itp_org%coef_inter_org(inod,1:3) = xi_ele(1:3)
            end if
          end if
        end do
      end do
!
      end subroutine set_itp_fine_to_course_origin
!
! ----------------------------------------------------------------------
!
      subroutine set_itp_fine_to_course_dest(itp_dest)
!
      type(interpolate_table_dest), intent(inout) :: itp_dest
      integer(kind = kint) :: inod
!
!
      itp_dest%iflag_self_itp_recv = ione
!
      call set_num_org_domain(ione, itp_dest)
      call alloc_itp_num_dest(itp_dest)
!
      itp_dest%id_org_domain(1) =       izero
      itp_dest%istack_nod_tbl_dest(0) = izero
      itp_dest%istack_nod_tbl_dest(1) = ntot_nod_refine_nod
      itp_dest%ntot_table_dest =        ntot_nod_refine_nod
!
      call alloc_itp_table_dest(itp_dest)
!
      do inod = 1, ntot_nod_refine_nod
        itp_dest%inod_dest_4_dest(inod) = inod
      end do
!
      end subroutine set_itp_fine_to_course_dest
!
! ----------------------------------------------------------------------
!
      end module set_refine_interpolate_tbl
