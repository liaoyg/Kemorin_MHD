!>@file   set_control_4_filtering.f90
!!@brief  module set_control_4_filtering
!!
!!@author H. Matsui
!!@date Programmed by H. Matsui in July, 2012
!
!> @brief set filtering parameters for SGS model
!!        from control data
!!
!!@verbatim
!!      subroutine s_set_control_4_filtering                            &
!!     &         (SGS_filter_name_ctl, ffile_ctl, s3df_ctl)
!!        type(read_character_item), intent(in) :: SGS_filter_name_ctl
!!        type(filter_file_control), intent(in) :: ffile_ctl
!!        type(SGS_3d_filter_control), intent(inout) :: s3df_ctl
!!@endverbatim
!
      module set_control_4_filtering
!
      use m_precision
!
      implicit  none
!
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine s_set_control_4_filtering                              &
     &         (SGS_filter_name_ctl, ffile_ctl, s3df_ctl)
!
      use calypso_mpi
      use m_constants
      use m_error_IDs
      use m_machine_parameter
      use m_file_format_switch
      use m_phys_labels
      use m_control_parameter
      use m_filter_file_names
      use t_ctl_data_SGS_filter
      use t_ctl_data_filter_files
      use t_read_control_arrays
      use sgs_ini_model_coefs_IO
      use set_control_ele_layering
      use skip_comment_f
!
      type(read_character_item), intent(in) :: SGS_filter_name_ctl
      type(filter_file_control), intent(in) :: ffile_ctl
      type(SGS_3d_filter_control), intent(inout) :: s3df_ctl
!
      integer(kind = kint) :: i
      character(len=kchara) :: tmpchara
!
!
      if (SGS_filter_name_ctl%iflag .gt. 0) then
        tmpchara = SGS_filter_name_ctl%charavalue
        if   (cmp_no_case(tmpchara, '3D')                               &
     &   .or. cmp_no_case(tmpchara, '3-D')                              &
     &   .or. cmp_no_case(tmpchara, '3-Dimensional')                    &
     &  )  iflag_SGS_filter = id_SGS_3D_FILTERING
!
        if   (cmp_no_case(tmpchara, 'Line')                             &
     &  )  iflag_SGS_filter = id_SGS_LINE_FILTERING
!
        if   (cmp_no_case(tmpchara, 'Plane')                            &
     &  )  iflag_SGS_filter = id_SGS_PLANE_FILTERING
!
        if   (cmp_no_case(tmpchara, '3D_easy')                          &
     &  )  iflag_SGS_filter = id_SGS_3D_EZ_FILTERING
!
        if   (cmp_no_case(tmpchara, '3D_smp')                           &
     &   .or. cmp_no_case(tmpchara, '3-D_smp')                          &
     &   .or. cmp_no_case(tmpchara, '3-Dimensional_smp')                &
     &       )  iflag_SGS_filter = id_SGS_3D_SMP_FILTERING
!
        if   (cmp_no_case(tmpchara, '3D_easy_sm')                       &
     &  )  iflag_SGS_filter = id_SGS_3D_EZ_SMP_FILTERING
!
        if (iflag_debug .gt. 0)  write(*,*)                             &
     &       'iflag_SGS_filter',     iflag_SGS_filter
      end if
!
      if (iflag_SGS_model.eq.id_SGS_similarity                          &
     &     .or. iflag_dynamic_SGS .ne. id_SGS_DYNAMIC_OFF) then
        if (iflag_SGS_filter .eq. id_SGS_NO_FILTERING) then
          e_message = 'Set filtering type for dynamic model'
          call calypso_MPI_abort(ierr_SGS, e_message)
        end if
      end if
!
!
!
      if     (iflag_SGS_filter .eq. id_SGS_3D_FILTERING                 &
     &   .or. iflag_SGS_filter .eq. id_SGS_3D_EZ_FILTERING              &
     &   .or. iflag_SGS_filter .eq. id_SGS_3D_SMP_FILTERING             &
     &   .or. iflag_SGS_filter .eq. id_SGS_3D_EZ_SMP_FILTERING ) then
        if (s3df_ctl%whole_filter_grp_ctl%icou .gt. 0) then
          num_whole_filter_grp =   s3df_ctl%whole_filter_grp_ctl%num
        else
          num_whole_filter_grp =   1
        end if
        num_whole_w_filter_grp = num_whole_filter_grp
!
        call allocate_whole_filter_groups
!
        if (s3df_ctl%whole_filter_grp_ctl%icou .gt. 0) then
          do i = 1, num_whole_filter_grp
            whole_filter_grp(i)                                         &
     &         = s3df_ctl%whole_filter_grp_ctl%c_tbl(i)
          end do
        else
          whole_filter_grp(1) =   'all'
          whole_w_filter_grp(1) = 'all'
        end if
        whole_w_filter_grp(1:num_whole_filter_grp)                      &
     &         = whole_filter_grp(1:num_whole_filter_grp)
!
        call dealloc_control_array_chara(s3df_ctl%whole_filter_grp_ctl)
!
        if (s3df_ctl%fluid_filter_grp_ctl%icou .gt. 0) then
          num_fluid_filter_grp = s3df_ctl%fluid_filter_grp_ctl%num
        else
          num_fluid_filter_grp = 1
        end if
        num_fluid_w_filter_grp = num_fluid_filter_grp
!
        call allocate_fluid_filter_groups
!
        if (s3df_ctl%fluid_filter_grp_ctl%icou .gt. 0) then
          do i = 1, num_fluid_filter_grp
            fluid_filter_grp(i)                                         &
     &         = s3df_ctl%fluid_filter_grp_ctl%c_tbl(i)
          end do
        else
          fluid_filter_grp(1) = 'all'
        end if
        fluid_w_filter_grp(1:num_fluid_filter_grp)                      &
     &         = fluid_filter_grp(1:num_fluid_filter_grp)
!
        call dealloc_control_array_chara(s3df_ctl%fluid_filter_grp_ctl)
!
        if (evo_temp%iflag_scheme .gt. id_no_evolution) then
          iflag_heat_filtering = 0
          if (cmp_no_case(s3df_ctl%heat_filter_ctl%charavalue,          &
     &        'Whole_filtering')) iflag_heat_filtering = 0
          if (cmp_no_case(s3df_ctl%heat_filter_ctl%charavalue,          &
     &        'Fluid_filtering')) iflag_heat_filtering = 1
        end if
!
        if ( evo_velo%iflag_scheme .gt. id_no_evolution) then
          iflag_momentum_filtering = 0
          if (cmp_no_case(s3df_ctl%momentum_filter_ctl%charavalue,      &
     &        'Whole_filtering')) iflag_momentum_filtering = 0
          if (cmp_no_case(s3df_ctl%momentum_filter_ctl%charavalue,      &
     &        'Fluid_filtering')) iflag_momentum_filtering = 1
        end if
!
        if (evo_magne%iflag_scheme .gt. id_no_evolution                 &
     &      .or. evo_vect_p%iflag_scheme .gt. id_no_evolution) then
          iflag_induction_filtering = 0
          if (cmp_no_case(s3df_ctl%induction_filter_ctl%charavalue,     &
     &        'Whole_filtering')) iflag_induction_filtering = 0
          if (cmp_no_case(s3df_ctl%induction_filter_ctl%charavalue,     &
     &        'Fluid_filtering')) iflag_induction_filtering = 1
        end if
!
        if (iflag_debug.eq.1) then
          write(*,*) 'num_whole_filter_grp', num_whole_filter_grp
          write(*,*) 'whole_filter_grp'
          do i = 1, num_whole_filter_grp
            write(*,*) i, trim(whole_filter_grp(i))
          end do
          write(*,*) 'num_fluid_filter_grp', num_fluid_filter_grp
          write(*,*) 'fluid_filter_grp'
          do i = 1, num_fluid_filter_grp
            write(*,*) i, trim(fluid_filter_grp(i))
          end do
          write(*,*) 'iflag_heat_filtering',                            &
     &              iflag_heat_filtering
          write(*,*) 'iflag_momentum_filtering',                        &
     &              iflag_momentum_filtering
          write(*,*) 'iflag_induction_filtering',                       &
     &              iflag_induction_filtering
        end if
!
!     set filter file header
!
        if (ffile_ctl%filter_head_ctl%iflag .eq. 1) then
          filter_3d_head = ffile_ctl%filter_head_ctl%charavalue
        end if
!
        if (ffile_ctl%filter_wide_head_ctl%iflag .eq. 1) then
          filter_wide_head = ffile_ctl%filter_wide_head_ctl%charavalue
        end if
!
        call choose_file_format                                         &
     &     (ffile_ctl%filter_3d_format, ifmt_3d_filter)
        call choose_file_format                                         &
     &     (ffile_ctl%filter_wide_format, ifmt_wide_filter)
!
        if (iflag_debug .gt. 0)  then
          write(*,*) 'filter_3d_head: ',     trim(filter_3d_head)
          write(*,*) 'filter_3d_format: ',   ifmt_3d_filter
          write(*,*) 'filter_wide_head: ',   trim(filter_wide_head)
          write(*,*) 'filter_wide_format: ', ifmt_wide_filter
        end if
      end if
!
      if (iflag_SGS_filter .eq. id_SGS_LINE_FILTERING) then
        if (ffile_ctl%filter_head_ctl%iflag .eq. 1) then
          filter_line_head = ffile_ctl%filter_head_ctl%charavalue
        else
          filter_line_head = filter_l_def_hd
        end if
      end if
!
!
      end subroutine s_set_control_4_filtering
!
! -----------------------------------------------------------------------
!
      end module set_control_4_filtering
