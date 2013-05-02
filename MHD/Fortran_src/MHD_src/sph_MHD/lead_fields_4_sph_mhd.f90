!>@file   lead_fields_4_sph_mhd.f90
!!@brief  module lead_fields_4_sph_mhd
!!
!!@author H. Matsui
!!@date Programmed in Aug, 2007
!
!>@brief  Evaluate pressure and energy fluxes for snapshots
!!
!!@verbatim
!!      subroutine s_lead_fields_4_sph_mhd
!!      subroutine pressure_4_sph_mhd
!!      subroutine enegy_fluxes_4_sph_mhd
!!@endverbatim
!
      module lead_fields_4_sph_mhd
!
      use m_precision
      use m_machine_parameter
!
!
      implicit none
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine s_lead_fields_4_sph_mhd
!
      use m_t_step_parameter
      use output_viz_file_control
!
      integer (kind =kint) :: iflag
!
!
      call set_lead_physical_values_flag(iflag)
!
      if ( (iflag*mod(istep_max_dt,i_step_output_rst)) .eq.0 ) then
         call pressure_4_sph_mhd
      end if
!
      if(iflag .eq. 0) call enegy_fluxes_4_sph_mhd
!
      end subroutine s_lead_fields_4_sph_mhd
!
! ----------------------------------------------------------------------
!
      subroutine pressure_4_sph_mhd
!
      use cal_sol_sph_fluid_crank
!
      use cal_sph_field_by_rotation
      use const_radial_forces_on_bc
      use cal_div_of_forces
!
!
      if (iflag_debug.eq.1) write(*,*) 'cal_div_of_forces_sph_2'
      call cal_div_of_forces_sph_2
!
      call s_const_radial_forces_on_bc
      call s_cal_div_of_forces
!
      if (iflag_debug.eq.1) write(*,*) 'cal_sol_pressure_by_div_v'
      call cal_sol_pressure_by_div_v
!
      end subroutine pressure_4_sph_mhd
!
! ----------------------------------------------------------------------
!
      subroutine enegy_fluxes_4_sph_mhd
!
      use sph_transforms_4_MHD
      use cal_energy_flux_rtp
!
!
      if (iflag_debug.eq.1) write(*,*) 'sph_back_trans_snapshot_MHD'
      call sph_back_trans_snapshot_MHD
!
!
!
      if (iflag_debug.eq.1) write(*,*) 's_cal_energy_flux_rtp'
      call s_cal_energy_flux_rtp
      if (iflag_debug.eq.1) write(*,*)                                  &
     &                          'sph_forward_trans_snapshot_MHD'
      call sph_forward_trans_snapshot_MHD
!
      end subroutine enegy_fluxes_4_sph_mhd
!
! ----------------------------------------------------------------------
!
      end module lead_fields_4_sph_mhd
