!>@file   t_SGS_control_parameter.f90
!!@brief  module t_SGS_control_parameter
!!
!!@author H. Matsui and H. Okuda
!!@date Programmed in 2000
!!@n modified in Feb., 2009
!
!> @brief control flags for MHD dynamo model
!!
!!@verbatim
!!      subroutine alloc_icore_ele_grp_name
!!      subroutine alloc_whole_filter_groups
!!      subroutine alloc_fluid_filter_groups
!!
!!      subroutine dealloc_icore_ele_grp_name
!!      subroutine dealloc_whole_filter_groups
!!      subroutine dealloc_fluid_filter_groups
!!@endverbatim
!
      module   m_control_parameter
!
      use m_precision
      use t_time_stepping_parameter
!
      implicit  none
!
!
!>      Turn OFF flag
      integer (kind=kint), parameter :: id_turn_OFF = 0
!>      Turn ON flag
      integer (kind=kint), parameter :: id_turn_ON =  1
!
!
      integer (kind=kint), parameter :: id_SGS_none =       0
      integer (kind=kint), parameter :: id_SGS_NL_grad =    1
      integer (kind=kint), parameter :: id_SGS_similarity = 2
      integer (kind=kint), parameter :: id_SGS_diffusion =  3
!
      integer (kind=kint), parameter :: id_SGS_DYNAMIC_OFF =   0
      integer (kind=kint), parameter :: id_SGS_DYNAMIC_ON =    1
!
      type SGS_model_control_params
        integer (kind=kint) :: iflag_SGS_model = id_SGS_none
        integer (kind=kint) :: iflag_dynamic_SGS = id_SGS_DYNAMIC_OFF
!
        integer (kind=kint) :: iflag_SGS_heat =      id_SGS_none
        integer (kind=kint) :: iflag_SGS_inertia =   id_SGS_none
        integer (kind=kint) :: iflag_SGS_lorentz =   id_SGS_none
        integer (kind=kint) :: iflag_SGS_induction = id_SGS_none
        integer (kind=kint) :: iflag_SGS_comp_flux = id_SGS_none
        integer (kind=kint) :: iflag_SGS_gravity =   id_SGS_none
!
        real (kind = kreal) :: SGS_hf_factor =      1.0d0
        real (kind = kreal) :: SGS_mf_factor =      1.0d0
        real (kind = kreal) :: SGS_mawell_factor =  1.0d0
        real (kind = kreal) :: SGS_uxb_factor =     1.0d0
        real (kind = kreal) :: SGS_cf_factor =      1.0d0
!
        integer (kind=kint) :: min_step_dynamic =  1
        integer (kind=kint) :: max_step_dynamic =  1
        real (kind = kreal) :: delta_to_shrink_dynamic = 1.0d5
        real (kind = kreal) :: delta_to_extend_dynamic = 1.0d-5
!
        integer (kind=kint) :: iflag_SGS_parterbuation = 0
!
        integer (kind=kint) :: itype_SGS_model_coef =  0
        integer (kind=kint) :: icoord_SGS_model_coef = 0
!
        integer (kind=kint) :: itype_SGS_h_flux_coef =   0
        integer (kind=kint) :: itype_SGS_c_flux_coef =   0
        integer (kind=kint) :: itype_SGS_m_flux_coef =   0
        integer (kind=kint) :: itype_SGS_maxwell_coef =  0
        integer (kind=kint) :: itype_SGS_uxb_coef =      0
!
        integer (kind=kint) :: iset_SGS_nagetive_clip = 0
        integer (kind=kint) :: iset_SGS_coef_marging =  0
        real (kind = kreal) :: SGS_clipping_limit = 0.0d0
      end type SGS_model_control_params
!
!>      ID not to apply commutation error correction
      integer (kind=kint), parameter :: id_SGS_commute_OFF = 0
!>      ID to apply commutation error correction
      integer (kind=kint), parameter :: id_SGS_commute_ON =  1
!
      type commutation_control_params
!>      commutation error correction flag for system
        integer (kind=kint) :: iflag_commute_correction                 &
     &                      = id_SGS_commute_OFF
!>      commutation error correction flag for linear terms
        integer (kind=kint) :: iflag_commute_linear                     &
     &                      = id_SGS_commute_OFF
!>      commutation error correction flag for nonlinear terms
        integer (kind=kint) :: iflag_commute_nonlinar                   &
     &                      = id_SGS_commute_OFF
!
!>      commutation error correction flag for temperature
        integer (kind=kint) :: iflag_commute_temp                       &
     &                      = id_SGS_commute_OFF
!>      commutation error correction flag for velocity
        integer (kind=kint) :: iflag_commute_velo                       &
     &                      = id_SGS_commute_OFF
!>      commutation error correction flag for magnetic field
        integer (kind=kint) :: iflag_commute_magne                      &
     &                      = id_SGS_commute_OFF
!>      commutation error correction flag for composition variation
        integer (kind=kint) :: iflag_commute_composit                   &
     &                      = id_SGS_commute_OFF
!
!>      commutation error correction flag for heat flux
        integer (kind=kint) :: iflag_commute_heat                       &
     &                      = id_SGS_commute_OFF
!>      commutation error correction flag for momentum flux
        integer (kind=kint) :: iflag_commute_inertia                    &
     &                      = id_SGS_commute_OFF
!>      commutation error correction flag for heat flux
        integer (kind=kint) :: iflag_commute_lorentz                    &
     &                      = id_SGS_commute_OFF
!>      commutation error correction flag for magnetic induction
        integer (kind=kint) :: iflag_commute_induction                  &
     &                      = id_SGS_commute_OFF
!>      commutation error correction flag for composition flux
        integer (kind=kint) :: iflag_commute_c_flux                     &
     &                      = id_SGS_commute_OFF
      end type commutation_control_params
!
      integer (kind=kint), parameter :: id_SGS_NO_FILTERING =         0
      integer (kind=kint), parameter :: id_SGS_3D_FILTERING =         1
      integer (kind=kint), parameter :: id_SGS_3D_EZ_FILTERING =     11
      integer (kind=kint), parameter :: id_SGS_3D_SMP_FILTERING =    21
      integer (kind=kint), parameter :: id_SGS_3D_EZ_SMP_FILTERING = 31
!
      integer (kind=kint), parameter :: id_SGS_LINE_FILTERING =       2
      integer (kind=kint), parameter :: id_SGS_PLANE_FILTERING =      3
      integer (kind=kint), parameter :: id_SGS_IDEAL_SPH_LOWPASS =    4
!
!>      filter ID for @f$ s\Delta @f$  filter
      integer (kind=kint), parameter :: ifilter_2delta = 1
!>      filter ID for @f$ 4\Delta @f$  filter
      integer (kind=kint), parameter :: ifilter_4delta = 2
!>      filter ID to obtain SGS terms
!
!
        integer (kind=kint) :: ifilter_final = ifilter_2delta
!
        integer (kind=kint) :: iflag_SGS_filter = id_SGS_3D_FILTERING
        integer (kind=kint) :: iset_DIFF_model_coefs =  0
!
!
      type SGS_filter_area_params
        integer (kind=kint) :: num_f_group = 0
        integer (kind=kint), allocatable :: id_f_group(:)
        character (len=kchara), allocatable :: gourp_name(:)
      end type SGS_filter_area_params
!
      type SGS_filtering_params
        type(SGS_filter_area_params) :: whole
        type(SGS_filter_area_params) :: fluid
!
        type(SGS_filter_area_params) :: whole_wide
        type(SGS_filter_area_params) :: fluid_wide
!
        integer (kind=kint) :: iflag_heat_filtering =        0
        integer (kind=kint) :: iflag_composition_filtering = 0
        integer (kind=kint) :: iflag_momentum_filtering =    0
        integer (kind=kint) :: iflag_induction_filtering =   0
      end type SGS_filtering_params
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine alloc_whole_filter_groups
!
!
      call alloc_filter_group_param(num_whole_filter_grp, whole)
      call alloc_filter_group_param(num_whole_w_filter_grp, fluid_wide)
!
      end subroutine alloc_whole_filter_groups
!
!  ---------------------------------------------------------------------
!
      subroutine alloc_fluid_filter_groups
!
!
      call alloc_filter_group_param(num_fluid_filter_grp, fluid)
      call alloc_filter_group_param(num_fluid_w_filter_grp, fluid_wide)
!
      end subroutine alloc_fluid_filter_groups
!
!  ---------------------------------------------------------------------
!
      subroutine alloc_filter_group_param(num_grp, f_area)
!
      integer(kind = kint), intent(in) :: num_grp
      type(SGS_filter_area_params), intent(inout) :: f_area
!
!
      f_area%num_f_group = num_grp
      allocate(f_area%gourp_name(f_area%num_f_group))
      allocate(f_area%id_f_group(f_area%num_f_group))
      if(f_area%num_f_group .gt. 0) f_area%id_f_group = 0
!
      end subroutine alloc_icore_ele_grp_name
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine dealloc_whole_filter_groups(filter_param)
!
      type(SGS_filtering_params), intent(inout) :: filter_param
!
!
      call dealloc_filter_group_param(num_whole_filter_grp, whole)
      call dealloc_filter_group_param                                   &
     &   (num_whole_w_filter_grp, fluid_wide)
!
      end subroutine dealloc_whole_filter_groups
!
!  ---------------------------------------------------------------------
!
      subroutine dealloc_fluid_filter_groups(filter_param)
!
      type(SGS_filtering_params), intent(inout) :: filter_param
!
!
      call dealloc_filter_group_param(num_fluid_filter_grp, fluid)
      call dealloc_filter_group_param                                   &
     &   (num_fluid_w_filter_grp, fluid_wide)
!
      end subroutine dealloc_fluid_filter_groups
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine dealloc_filter_group_param(f_area)
!
      type(SGS_filter_area_params), intent(inout) :: f_area
!
!
      deallocate(f_area%gourp_name, f_area%id_f_group)
!
      end subroutine dealloc_icore_ele_grp_name
!
!  ---------------------------------------------------------------------
!
      end module t_SGS_control_parameter
