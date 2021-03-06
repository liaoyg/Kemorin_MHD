!
!      module set_element_list_4_filter
!
!     Written by H. Matsui on Oct., 2006
!
!      subroutine s_set_element_list_4_filter(ele, ele_grp)
!
      module set_element_list_4_filter
!
      use m_precision
      use m_constants
!
      use m_machine_parameter
      use m_element_list_4_filter
!
      use t_geometry_data
!
      implicit none
!
      integer(kind = kint), allocatable :: imark_ele_filter(:)
!
      private :: imark_ele_filter
!
      private :: allocate_mark_list_4_filter
      private :: deallocate_mark_list_4_filter
      private :: mark_ele_list_4_filter
      private :: count_ele_list_4_filter, set_ele_list_4_filter
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine s_set_element_list_4_filter(ele, ele_grp)
!
      use t_group_data
      use cal_minmax_and_stacks
!
      type(element_data), intent(in) :: ele
      type(group_data), intent(in) :: ele_grp
!
!
      call allocate_mark_list_4_filter(ele)
!
      call mark_ele_list_4_filter(ele, ele_grp)
!
      call count_ele_list_4_filter(ele)
!
      call allocate_ele_list_4_filter
!
      call set_ele_list_4_filter(ele)
!
      call deallocate_mark_list_4_filter
!
!
      call allocate_ele_smp_stk_filter(np_smp)
!
      call count_number_4_smp( np_smp, ione, nele_4_filter,             &
     &    iele_filter_smp_stack, maxele_filter_4_smp)

!
      end subroutine s_set_element_list_4_filter
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine allocate_mark_list_4_filter(ele)
!
      type(element_data), intent(in) :: ele
!
!
      allocate(imark_ele_filter(ele%numele) )
      imark_ele_filter(1:ele%numele) = 0
!
      end subroutine allocate_mark_list_4_filter
!
! ----------------------------------------------------------------------
!
      subroutine deallocate_mark_list_4_filter
!
      deallocate(imark_ele_filter)
!
      end subroutine deallocate_mark_list_4_filter
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine mark_ele_list_4_filter(ele, ele_grp)
!
      use t_group_data
      use m_ctl_params_4_gen_filter
      use skip_comment_f
!
      type(element_data), intent(in) :: ele
      type(group_data), intent(in) :: ele_grp
      integer(kind = kint) :: i, igrp, inum, ist, ied, iele
!
!
      if (cmp_no_case(filter_area_name(1), 'all')) then
        id_filter_area_grp(1) = -1
        imark_ele_filter(1:ele%numele) = 1
      else
!
        do igrp = 1, num_filtering_grp
          do i = 1, ele_grp%num_grp
            if(filter_area_name(igrp) .eq. ele_grp%grp_name(i)) then
              id_filter_area_grp(igrp) = i
              ist = ele_grp%istack_grp(i-1) + 1
              ied = ele_grp%istack_grp(i)
              do inum = ist, ied
                iele = ele_grp%item_grp(inum)
                imark_ele_filter(iele) = 1
              end do
            end if
          end do
        end do
!
      end if
!
      end subroutine mark_ele_list_4_filter
!
! ----------------------------------------------------------------------
!
      subroutine count_ele_list_4_filter(ele)
!
      type(element_data), intent(in) :: ele
      integer(kind = kint) :: iele
!
!
      nele_4_filter = 0
      do iele = 1, ele%numele
        nele_4_filter = nele_4_filter + imark_ele_filter(iele)
      end do
!
      end subroutine count_ele_list_4_filter
!
! ----------------------------------------------------------------------
!
      subroutine set_ele_list_4_filter(ele)
!
      type(element_data), intent(in) :: ele
      integer(kind = kint) :: inum, iele
!
!
      inum = 0
      do iele = 1, ele%numele
        if (imark_ele_filter(iele) .eq. 1) then
          inum = inum + 1
          iele_4_filter(inum) = iele
        end if
      end do
!
      end subroutine set_ele_list_4_filter
!
! ----------------------------------------------------------------------
!
      end module set_element_list_4_filter
