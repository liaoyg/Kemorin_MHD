!m_control_param_4_refiner.f90
!      module m_control_param_4_refiner
!
!      Written by Kemorin on Oct., 2007
!
!      subroutine allocate_refine_param
!      subroutine deallocate_refine_param
!      subroutine deallocate_refine_param_chara
!
!      subroutine set_control_4_refiner
!
      module m_control_param_4_refiner
!
      use m_precision
      use t_file_IO_parameter
!
      implicit    none
!
!
      integer (kind = kint) :: iflag_interpolate_type
!
      integer(kind = kint) :: num_refine_type = 0
      character (len = kchara), allocatable :: refined_ele_grp(:)
      character (len = kchara), allocatable :: refined_ele_type(:)
      integer (kind=kint), allocatable :: id_refined_ele_grp(:)
      integer (kind=kint), allocatable :: iflag_refine_type(:)
      integer (kind=kint) :: iflag_redefine_tri
!
      integer(kind = kint) :: iflag_small_tri_refine = 1
!
      type(field_IO_params), save ::  original_mesh_file
      type(field_IO_params), save ::  refined_mesh_file
!
      character(len=kchara) :: course_2_fine_head = 'course_2_fine'
      character(len=kchara) :: fine_2_course_head = 'fine_2_course'
!
      character(len = kchara) :: refine_info_head = 'refine_info'
      character(len = kchara) :: old_refine_info_head = 'refine_info.0'
      integer(kind = kint) :: iflag_read_old_refine_file = 0
!
      private :: allocate_refine_param
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine allocate_refine_param
!
      allocate(refined_ele_grp(num_refine_type) )
      allocate(refined_ele_type(num_refine_type) )
      allocate(id_refined_ele_grp(num_refine_type) )
      allocate(iflag_refine_type(num_refine_type) )
!
      id_refined_ele_grp = 0
      iflag_refine_type = 0
!
      end subroutine allocate_refine_param
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine deallocate_refine_param
!
      deallocate(id_refined_ele_grp)
      deallocate(iflag_refine_type)
!
      end subroutine deallocate_refine_param
!
! -----------------------------------------------------------------------
!
      subroutine deallocate_refine_param_chara
!
      deallocate(refined_ele_grp)
      deallocate(refined_ele_type)
!
      end subroutine deallocate_refine_param_chara
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine set_control_4_refiner
!
      use m_control_data_4_refine
      use m_default_file_prefix
      use skip_comment_f
      use set_control_platform_data
!
      character(len = kchara) :: tmpchara
!
!
      call set_control_mesh_def(source_plt, original_mesh_file)
      call set_control_mesh_file_def                                    &
     &   (def_new_mesh_head, refined_plt, refined_mesh_file)
!
!
      if (coarse_2_fine_head_ctl%iflag .gt. 0) then
        course_2_fine_head = coarse_2_fine_head_ctl%charavalue
      end if
!
      if (fine_2_course_head_ctl%iflag .gt. 0) then
        fine_2_course_head = fine_2_course_head_ctl%charavalue
      end if
!
      if (refine_info_head_ctl%iflag .gt. 0) then
        refine_info_head = refine_info_head_ctl%charavalue
      end if
!
      iflag_read_old_refine_file = old_refine_info_head_ctl%iflag
      if (iflag_read_old_refine_file .gt. 0) then
        old_refine_info_head = old_refine_info_head_ctl%charavalue
      end if
!
!
      iflag_interpolate_type = 0
      if (interpolate_type_ctl%iflag .gt. 0) then
        tmpchara = interpolate_type_ctl%charavalue
        if (   cmp_no_case(tmpchara, 'project_sphere')                  &
     &    .or. cmp_no_case(tmpchara, 'project_to_sphere')) then
          iflag_interpolate_type = 2
!        else if (cmp_no_case(tmpchara, 'rtp')                          &
!     &      .or. cmp_no_case(tmpchara, 'spherical')) then
!          iflag_interpolate_type = 1
        end if
      end if
!
!
      if (refined_ele_grp_ctl%icou .gt. 0) then
        num_refine_type = refined_ele_grp_ctl%num
      else if (refine_i_ele_grp_ctl%icou .gt. 0) then
        num_refine_type = refine_i_ele_grp_ctl%num
      else
        write(*,*) 'set refine type and area'
        stop
      end if
!
      if (num_refine_type .gt. 0) then
        call allocate_refine_param
!
        if (refined_ele_grp_ctl%icou .gt. 0) then
!
          iflag_redefine_tri = 1
          refined_ele_grp(1:num_refine_type)                            &
     &         = refined_ele_grp_ctl%c1_tbl(1:num_refine_type)
          refined_ele_type(1:num_refine_type)                           &
     &         = refined_ele_grp_ctl%c2_tbl(1:num_refine_type)
          call alloc_control_array_c2(refined_ele_grp_ctl)
!
        else if (refine_i_ele_grp_ctl%icou .gt. 0) then
!
          iflag_redefine_tri = 0
          refined_ele_grp(1:num_refine_type)                            &
     &         = refine_i_ele_grp_ctl%c_tbl(1:num_refine_type)
          iflag_refine_type(1:num_refine_type)                          &
     &         = refine_i_ele_grp_ctl%ivec(1:num_refine_type)
          call dealloc_control_array_c_i(refine_i_ele_grp_ctl)
        end if
!
      else
        write(*,*) 'set refine type and area'
        stop
      end if
!
      end subroutine set_control_4_refiner
!
! -----------------------------------------------------------------------
!
      end module m_control_param_4_refiner
