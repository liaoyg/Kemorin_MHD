!set_merged_refine_itp.f90
!      module set_merged_refine_itp
!
!      subroutine set_merged_itp_course_to_fine                         &
!     &          (nnod_4_ele, nnod_2, nnod_4_ele_2, xx_2)
!      subroutine set_merged_itp_fine_to_course(nnod_4_ele)
!
      module set_merged_refine_itp
!
      use m_precision
!
      use m_constants
      use m_machine_parameter
      use m_work_merge_refine_itp
!
      implicit none
!
      private :: copy_merge_itp_table_refine
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine set_merged_itp_course_to_fine                          &
     &          (nnod_4_ele, nnod_2, nnod_4_ele_2, xx_2)
!
      use t_interpolate_tbl_org
      use set_refine_interpolate_tbl
      use merge_refine_itp_table
!
      integer(kind = kint), intent(in) :: nnod_4_ele
      integer(kind = kint), intent(in) :: nnod_2, nnod_4_ele_2
      real(kind = kreal), intent(in) :: xx_2(nnod_2,3)
!
!
      write(*,*) 'set_merged_refine_data_org'
      call set_merged_refine_data_org                                   &
     &   (nnod_4_ele, nnod_2, nnod_4_ele_2, xx_2)
      write(*,*) 'set_itp_course_to_fine_dest'
      call set_itp_course_to_fine_dest(nnod_2)
!
      write(*,*) 'copy_merge_itp_table_refine'
      call copy_merge_itp_table_refine                                  &
     &   (nnod_4_ele, nnod_2, nnod_4_ele_2, xx_2)
!
      call dealloc_type_itp_table_org(c2f_mgd%tbl_org)
      call deallocate_itp_num_org(c2f_mgd%tbl_org)
!
      end subroutine set_merged_itp_course_to_fine
!
! ----------------------------------------------------------------------
!
      subroutine copy_merge_itp_table_refine                            &
     &         (nnod_4_ele, nnod_2, nnod_4_ele_2, xx_2)
!
      use m_interpolate_table_orgin
      use merge_refine_itp_table
!
      integer(kind = kint), intent(in) :: nnod_4_ele
      integer(kind = kint), intent(in) :: nnod_2, nnod_4_ele_2
      real(kind = kreal), intent(in) :: xx_2(nnod_2,3)
!
!
      itp1_org%iflag_self_itp_send = c2f_mgd%tbl_org%iflag_self_itp_send
      call set_num_dest_domain                                          &
     &   (c2f_mgd%tbl_org%num_dest_domain, itp1_org)
      call alloc_type_itp_num_org(np_smp, itp1_org)
!
      ntot_table_org = c2f_mgd%tbl_org%ntot_table_org
!
      call allocate_itp_table_org
!
!
      itp1_org%id_dest_domain =     c2f_mgd%tbl_org%id_dest_domain
      itp1_org%istack_nod_tbl_org = c2f_mgd%tbl_org%istack_nod_tbl_org
!
      itp1_org%istack_itp_type_org(0) = izero
!
      call set_merged_refine_data_org                                   &
     &   (nnod_4_ele, nnod_2, nnod_4_ele_2, xx_2)
!
      end subroutine copy_merge_itp_table_refine
!
! ----------------------------------------------------------------------
!
      subroutine set_merged_itp_fine_to_course(nnod_4_ele)
!
      use m_refined_element_data
      use m_interpolate_table_orgin
      use copy_local_position_2_ele
!
      integer(kind = kint), intent(in) :: nnod_4_ele
!
      integer(kind = kint) :: iele, k1, inod
      real(kind = kreal) :: xi_ele(3)
!
!
      itp1_org%iflag_self_itp_send = ione
      call set_num_dest_domain(ione, itp1_org)
      call alloc_type_itp_num_org(np_smp, itp1_org)
!
      itp1_org%id_dest_domain(1) = izero
      itp1_org%istack_itp_type_org(0) = izero
      itp1_org%istack_itp_type_org(1) = node_org_refine%numnod
      itp1_org%istack_itp_type_org(2) = node_org_refine%numnod
      itp1_org%istack_itp_type_org(3) = node_org_refine%numnod
      itp1_org%istack_itp_type_org(4) = node_org_refine%numnod
      itp1_org%istack_nod_tbl_org(0) = izero
      itp1_org%istack_nod_tbl_org(1) = node_org_refine%numnod
      ntot_table_org =        node_org_refine%numnod
!
      call allocate_itp_table_org
      itype_inter_org(1:ntot_table_org) = -1
!
      do iele = 1, ntot_ele_refined
!
        do k1 = 1, nnod_4_ele
          inod = ie_refined(iele,k1)
!
          if(inod.le.node_org_refine%numnod) then
            if(itype_inter_org(inod).eq.-1) then
              call copy_node_local_posi_2_element(k1, xi_ele)
!
              inod_gl_dest_4_org(inod) = inod
              iele_org_4_org(inod) =     iele
              itype_inter_org(inod) =     k1
              coef_inter_org(inod,1:3) = xi_ele(1:3)
            end if
          end if
        end do
      end do
!
      end subroutine set_merged_itp_fine_to_course
!
! ----------------------------------------------------------------------
!
      end module set_merged_refine_itp
