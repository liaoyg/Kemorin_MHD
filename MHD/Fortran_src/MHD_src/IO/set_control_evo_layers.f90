!>@file   set_control_evo_layers.f90
!!@brief  module set_control_evo_layers
!!
!!@author H. Matsui
!!@date Programmed in 2002
!
!> @brief Set parameters for simulation areas from control data
!!
!!@verbatim
!!      subroutine s_set_control_evo_layers                             &
!!     &        (earea_ctl, fl_prop, cd_prop, ht_prop, cp_prop, FEM_prm)
!!        type(fluid_property), intent(in) :: fl_prop
!!        type(scalar_property), intent(in) :: ht_prop, cp_prop
!!        type(conductive_property), intent(inout)  :: cd_prop
!!        type(mhd_evo_area_control), intent(inout) :: earea_ctl
!!        type(FEM_MHD_paremeters), intent(inout) :: FEM_prm
!!@endverbatim
!
!
      module set_control_evo_layers
!
      use m_precision
!
      use m_machine_parameter
      use t_physical_property
      use t_ctl_data_mhd_evolution
      use t_FEM_control_parameter
!
      implicit  none
!
      private :: set_fluid_layer_egrp_name
      private :: set_conduct_layer_egrp_name
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine s_set_control_evo_layers                               &
     &        (earea_ctl, fl_prop, cd_prop, ht_prop, cp_prop, FEM_prm)
!
      type(fluid_property), intent(in) :: fl_prop
      type(conductive_property), intent(inout)  :: cd_prop
      type(scalar_property), intent(in) :: ht_prop, cp_prop
      type(mhd_evo_area_control), intent(inout) :: earea_ctl
      type(FEM_MHD_paremeters), intent(inout) :: FEM_prm
!
!
      if       (fl_prop%iflag_scheme .eq. id_no_evolution               &
     &    .and. ht_prop%iflag_scheme .eq. id_no_evolution               &
     &    .and. cp_prop%iflag_scheme .eq. id_no_evolution) then
          call alloc_area_group_name(ione, FEM_prm%fluid_group)
          FEM_prm%fluid_group%group_name = 'none'
!
          call set_conduct_layer_egrp_name(earea_ctl, FEM_prm)
!
      else
        call set_fluid_layer_egrp_name(earea_ctl, FEM_prm)
!
        if     (cd_prop%iflag_Bevo_scheme .eq. id_no_evolution          &
     &    .and. cd_prop%iflag_Aevo_scheme .eq. id_no_evolution) then
          call alloc_area_group_name(ione, FEM_prm%condutive_group)
          FEM_prm%condutive_group%group_name = 'none'
!
        else
          call set_conduct_layer_egrp_name(earea_ctl, FEM_prm)
        end if
      end if
!
      call alloc_area_group_name(izero, FEM_prm%inner_core_group)
!
      end subroutine s_set_control_evo_layers
!
! -----------------------------------------------------------------------
!
      subroutine set_fluid_layer_egrp_name(earea_ctl, FEM_prm)
!
      type(mhd_evo_area_control), intent(inout) :: earea_ctl
      type(FEM_MHD_paremeters), intent(inout) :: FEM_prm
!
!
      if (earea_ctl%evo_fluid_group_ctl%icou .eq. 0) then
        call alloc_area_group_name(ione, FEM_prm%fluid_group)
        FEM_prm%fluid_group%group_name = 'all'
      else
        call alloc_area_group_name(earea_ctl%evo_fluid_group_ctl%num,   &
     &      FEM_prm%fluid_group)
        FEM_prm%fluid_group%group_name                                  &
     &        =  earea_ctl%evo_fluid_group_ctl%c_tbl
        call dealloc_ele_fl_grp_ctl(earea_ctl)
      end if
!
      end subroutine set_fluid_layer_egrp_name
!
! -----------------------------------------------------------------------
!
      subroutine set_conduct_layer_egrp_name(earea_ctl, FEM_prm)
!
      type(mhd_evo_area_control), intent(inout) :: earea_ctl
      type(FEM_MHD_paremeters), intent(inout) :: FEM_prm
!
!
      if (earea_ctl%evo_conduct_group_ctl%icou .eq. 0) then
        call alloc_area_group_name(ione, FEM_prm%condutive_group)
        FEM_prm%condutive_group%group_name =  'all'
      else
        call alloc_area_group_name(earea_ctl%evo_conduct_group_ctl%num, &
     &      FEM_prm%condutive_group)
        FEM_prm%condutive_group%group_name                              &
     &            =  earea_ctl%evo_conduct_group_ctl%c_tbl
        call dealloc_ele_cd_grp_ctl(earea_ctl)
      end if
!
      end subroutine set_conduct_layer_egrp_name
!
! -----------------------------------------------------------------------
!
      end module set_control_evo_layers
