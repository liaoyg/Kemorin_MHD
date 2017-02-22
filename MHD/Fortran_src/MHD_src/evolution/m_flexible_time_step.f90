!m_flexible_time_step.f90
!     module m_flexible_time_step
!
!      Written by H. Matsui on Nov., 2009
!
!!      subroutine set_new_time_and_step(cd_prop, iphys, nod_fld)
!!      subroutine s_check_flexible_time_step(node, ele, fluid, cd_prop,&
!!     &          iphys, nod_fld, jac_3d_q, jac_3d_l, fem_wk)
!!        type(conductive_property), intent(in) :: cd_prop
!!        type(node_data), intent(in) :: node
!!        type(element_data), intent(in) :: ele
!!        type(field_geometry_data), intent(in) :: fluid
!!        type(phys_address), intent(in) :: iphys
!!        type(phys_data), intent(in) :: nod_fld
!!        type(jacobians_3d), intent(in) :: jac_3d_q, jac_3d_l
!!        type(work_finite_element_mat), intent(inout) :: fem_wk
!
      module m_flexible_time_step
!
      use m_precision
      use calypso_mpi
!
      use m_constants
      use m_machine_parameter
      use m_t_step_parameter
      use m_t_int_parameter
      use m_control_parameter
!
      use t_flex_delta_t_data
!
      implicit  none
!
      integer(kind=kint), parameter :: dt_check_max_code =   15
      integer(kind=kint), parameter :: dt_check_min_code =   17
!
      character(len=kchara), parameter                                  &
     &      :: dt_check_max_name = 'maximum_dt_chack.dat'
      character(len=kchara), parameter                                  &
     &      :: dt_check_min_name = 'minimum_dt_chack.dat'
!
      type(flexible_steppind_data), save :: flex_data
!
!
      private :: dt_check_max_code, dt_check_min_code
      private :: dt_check_max_name, dt_check_min_name
      private :: shrink_delta_t, extend_delta_t
      private :: open_flex_step_monitor
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine set_new_time_and_step(cd_prop, iphys, nod_fld)
!
      use t_material_property
      use t_phys_data
      use t_phys_address
!
      use copy_field_data_4_dt_check
!
      type(conductive_property), intent(in) :: cd_prop
      type(phys_address), intent(in) :: iphys
      type(phys_data), intent(inout) :: nod_fld
!
!
      time = time + dt
      i_step_MHD = i_step_MHD + 1
!
      if (iflag_flexible_step .eq. iflag_fixed_step) then
        istep_max_dt = i_step_MHD
      else
        istep_flex_to_max = istep_flex_to_max + 1
        istep_flex_to_max = mod(istep_flex_to_max,i_interval_flex_2_max)
        if(istep_flex_to_max .eq. 0) istep_max_dt = istep_max_dt + 1
      end if
!
      if (iflag_debug.eq.1) write(*,*) 's_copy_field_data_for_dt_check'
      call s_copy_field_data_for_dt_check(cd_prop, iphys, nod_fld)
!
      end subroutine set_new_time_and_step
!
! -----------------------------------------------------------------------
!
      subroutine s_check_flexible_time_step(node, ele, fluid, cd_prop,  &
     &          iphys, nod_fld, jac_3d_q, jac_3d_l, fem_wk, flex_data)
!
      use t_geometry_data_MHD
      use t_geometry_data
      use t_phys_data
      use t_phys_address
      use t_jacobian_3d
      use t_finite_element_mat
      use t_flex_delta_t_data
      use t_physical_property
!
      use check_deltat_by_prev_rms
!
      type(node_data), intent(in) :: node
      type(element_data), intent(in) :: ele
      type(field_geometry_data), intent(in) :: fluid
      type(conductive_property), intent(in) :: cd_prop
      type(phys_address), intent(in) :: iphys
      type(phys_data), intent(in) :: nod_fld
      type(jacobians_3d), intent(in) :: jac_3d_q, jac_3d_l
      type(work_finite_element_mat), intent(inout) :: fem_wk
      type(flexible_steppind_data), intent(inout) :: flex_data
!
!
      if( mod(istep_flex_to_max,itwo) .eq. izero) then
!        call s_check_deltat_by_previous                                &
!     &     (node, cd_prop1, iphys, nod_fld, flex_data)
        call s_check_deltat_by_prev_rms(node, ele, fluid, cd_prop,      &
     &      iphys, nod_fld, jac_3d_q, jac_3d_l, fem_wk, flex_data)
!
        if(flex_data%d_ratio_allmax .gt. min_eps_to_expand_dt) then
          call shrink_delta_t
          iflag_flex_step_changed = 1
          return
        else
          iflag_flex_step_changed = 0
        end if
!
        if(istep_flex_to_max .eq. 0) then
          if(flex_data%d_ratio_allmax .lt. max_eps_to_shrink_dt) then
            call extend_delta_t
            iflag_flex_step_changed = 1
          end if
!
          if(my_rank .eq. izero) then
            call open_flex_step_monitor
            call write_rms_delta_t_check(dt_check_max_code,             &
     &        i_step_MHD, time, flex_data)
            close(dt_check_max_code)
          end if
        end if
      end if
!
      end subroutine s_check_flexible_time_step
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine shrink_delta_t
!
!
      if(my_rank .eq. izero) then
        write(*,*) 'Shrink Delta t from ', dt, ' at ', i_step_MHD,      &
     &            'd_ratio_max', flex_data%d_ratio_allmax
      end if
      if(iflag_debug .gt. izero) then
        write(*,*) 'Old temporal step is ', istep_flex_to_max,          &
     &             'to ', istep_flex_to_max
      end if
!
      if     (dt_fact .eq. one) then
        dt_fact =   five
        idt_digit = idt_digit - 1
        i_interval_flex_2_max = i_interval_flex_2_max * 2
        istep_flex_to_max = istep_flex_to_max * 2
      else if(dt_fact .eq. two) then
        dt_fact = one
        i_interval_flex_2_max = (i_interval_flex_2_max * 5) / 2
        istep_flex_to_max = (istep_flex_to_max * 5) / 2
      else if(dt_fact .eq. five) then
        dt_fact = two
        i_interval_flex_2_max = i_interval_flex_2_max * 2
        istep_flex_to_max = istep_flex_to_max * 2
      end if
!
      dt = dt_fact * ten**(idt_digit)
      ddt = one / dt
!
      if(my_rank .eq. izero) then
        write(*,*) 'New Delta t is ', dt
      end if
      if(iflag_debug .gt. izero) then
        write(*,*) 'New temporal step is ', istep_flex_to_max,          &
     &             'to ', istep_flex_to_max
      end if
!
      end subroutine shrink_delta_t
!
! -----------------------------------------------------------------------
!
      subroutine extend_delta_t
!
!
      if(my_rank .eq. izero) then
        write(*,*) 'Extend Delta t from ', dt, ' at ', i_step_MHD,      &
     &            'd_ratio_max', flex_data%d_ratio_allmax
      end if
      if(iflag_debug .gt. izero) then
        write(*,*) 'Old temporal step is ', istep_flex_to_max,          &
     &             'to ', istep_flex_to_max
      end if
!
      if     (dt_fact .eq. one) then
        dt_fact =   two
        i_interval_flex_2_max = i_interval_flex_2_max / 2
        istep_flex_to_max =     istep_flex_to_max / 2
      else if(dt_fact .eq. two) then
        dt_fact = five
        i_interval_flex_2_max = (i_interval_flex_2_max * 2) / 5
        istep_flex_to_max =     (istep_flex_to_max * 2) / 5
      else if(dt_fact .eq. five) then
        dt_fact = one
        idt_digit = idt_digit + 1
        i_interval_flex_2_max = i_interval_flex_2_max / 2
        istep_flex_to_max =     istep_flex_to_max / 2
      end if
!
      dt = dt_fact * ten**(idt_digit)
      ddt = one / dt
!
      if(my_rank .eq. izero) then
        write(*,*) 'New Delta t is ', dt
      end if
      if(iflag_debug .gt. izero) then
        write(*,*) 'New temporal step is ', istep_flex_to_max,          &
     &             'to ', istep_flex_to_max
      end if
!
      end subroutine extend_delta_t
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine open_flex_step_monitor
!
!
      open(dt_check_max_code, file=dt_check_max_name,                   &
     &    form='formatted', status='old', position='append', err = 99)
      return
!
  99  continue
      open(dt_check_max_code, file = dt_check_max_name)
      call write_delta_t_check_head(dt_check_max_code, flex_data)
!
      end subroutine open_flex_step_monitor
!
! -----------------------------------------------------------------------
!
      end module m_flexible_time_step
