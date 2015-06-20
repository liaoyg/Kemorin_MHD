!
!      module correct_wrong_filters
!
      module correct_wrong_filters
!
!     Written by H. Matsui on Nov., 2008
!
      use m_precision
!
      use m_constants
      use m_ctl_params_4_gen_filter
      use m_geometry_parameter
      use m_filter_coefs
!
      use expand_filter_area_4_1node
      use copy_moments_2_matrix
      use cal_filter_func_each_node
      use cal_simple_filter_each_node
      use cal_3d_filter_4_each_node
      use cal_filter_moments_again
      use write_filters_4_each_node
!
!
      implicit none
!
!      subroutine s_correct_wrong_filters(id_filter_coef)
!      subroutine correct_wrong_fluid_filters(id_filter_coef)
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine s_correct_wrong_filters(id_filter_coef, dxidxs)
!
      use t_filter_dxdxi
!
      integer(kind = kint), intent(in) :: id_filter_coef
      type(dxidx_data_type), intent(inout) :: dxidxs
      integer(kind = kint) :: inod, ierr2, ierr
!
!
      if (inod_end_filter .eq. -1) then
        inod_end_filter = internal_node
      end if
!
      call init_4_cal_fileters
!
      write(70+my_rank,*) ' Best condition for filter'
!
      do inod = inod_start_filter, inod_end_filter
        call read_each_filter_stack_coef(id_filter_coef)
!
        call cal_rms_filter_coefs(min_rms_weight, ierr2)
!
        if (min_rms_weight .gt. max_rms_weight_limit) then
          iflag_make_whole_filter(inod) = 1
          iflag_make_fluid_filter(inod) = 1
        end if
!
        if (ierr2 .eq. -1) then
          iflag_make_whole_filter(inod) = 1
          iflag_make_fluid_filter(inod) = 1
        end if
!
        if( iflag_make_whole_filter(inod) .eq. 0) then
          call write_each_filter_stack_coef(inod)
!
          if(iflag_tgt_filter_type .ge. -4                              &
     &      .and. iflag_tgt_filter_type.le. -2) then
            call s_cal_filter_moments_again(inod)
          end if
        else
!
          if (iflag_tgt_filter_type .eq. -1) then
            call copy_filter_coefs_to_tmp
            call const_filter_func_nod_by_nod(inod, ierr)
          else if(iflag_tgt_filter_type .ge. -4                         &
     &      .and. iflag_tgt_filter_type.le. -2) then
            call set_simple_filter_nod_by_nod(inod, dxidxs%dx_nod)
          end if
!
        end if
!
      end do
!
      end subroutine s_correct_wrong_filters
!
! -----------------------------------------------------------------------
!
      subroutine correct_wrong_fluid_filters(id_filter_coef, dxidxs)
!
      use t_filter_dxdxi
!
      integer(kind = kint), intent(in) :: id_filter_coef
      type(dxidx_data_type), intent(inout) :: dxidxs
!
      integer(kind = kint) :: inod, ierr2, ierr
!
!
!
      call init_4_cal_fluid_fileters
!
      write(70+my_rank,*) ' Best condition for fluid filter'
!
      do inod = inod_start_filter, inod_end_filter
!
        call read_each_filter_stack_coef(id_filter_coef)
!
        if ( nnod_near_1nod_weight .gt. 0) then
          call cal_rms_filter_coefs(min_rms_weight, ierr2)
!
          if (min_rms_weight .gt. max_rms_weight_limit) then
            iflag_make_fluid_filter(inod) = 1
          end if
!
        end if
!
!     no correction
!
        if( iflag_make_fluid_filter(inod) .eq. 0) then
!
          if (nnod_near_1nod_weight .eq. 0) then
            call write_each_no_filter_coef(inod)
          else if (nnod_near_1nod_weight .lt. 0) then
            nnod_near_1nod_weight = -nnod_near_1nod_weight
            call write_each_same_filter_coef(inod)
          else
            call write_each_filter_stack_coef(inod)
          end if
!
!       correct fluid filter
!
        else
!
          if (iflag_tgt_filter_type .eq. -1) then
            call copy_filter_coefs_to_tmp
            call const_fluid_filter_nod_by_nod(inod, ierr)
          else if(iflag_tgt_filter_type .ge. -4                         &
     &      .and. iflag_tgt_filter_type.le. -2) then
            call set_simple_fl_filter_nod_by_nod(inod, dxidxs%dx_nod)
          end if
!
        end if
!
      end do
!
      end subroutine correct_wrong_fluid_filters
!
! -----------------------------------------------------------------------
!
      end module correct_wrong_filters
