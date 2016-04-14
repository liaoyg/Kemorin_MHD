!>@file   cal_sol_sph_MHD_crank.f90
!!@brief  module cal_sol_sph_MHD_crank
!!
!!@author H. Matsui
!!@date   Programmed  H. Matsui in Oct., 2009
!
!>@brief  Update fields for MHD dynamo model
!!
!!@verbatim
!!      subroutine s_cal_sol_sph_MHD_crank(rj_fld)
!!      subroutine set_sph_field_to_start(rj_fld)
!!        type(phys_data), intent(inout) :: rj_fld
!!
!!      subroutine check_ws_spectr((rj_fld)
!!        type(phys_data), intent(in) :: rj_fld
!!@endverbatim
!
      module cal_sol_sph_MHD_crank
!
      use m_precision
!
      use calypso_mpi
      use m_machine_parameter
      use m_control_parameter
      use m_spheric_parameter
      use m_spheric_param_smp
      use m_radial_matrices_sph
      use m_sph_phys_address
      use m_physical_property
      use const_sph_radial_grad
      use const_sph_rotation
      use const_sph_diffusion
!
      use t_phys_data
!
      implicit none
!
      private :: update_after_vorticity_sph
      private :: update_after_heat_sph, update_after_composit_sph
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine s_cal_sol_sph_MHD_crank(rj_fld)
!
      use m_boundary_params_sph_MHD
      use cal_rot_buoyancies_sph_MHD
      use cal_sol_sph_fluid_crank
      use const_sph_radial_grad
!
      type(phys_data), intent(inout) :: rj_fld
!
!      integer(kind = kint) :: j, k, inod
!
!*-----  time evolution   -------------
!*
!      call check_ws_spectr(rj_fld)
!
      if(iflag_t_evo_4_velo .gt. id_no_evolution) then
!         Input:    ipol%i_vort, itor%i_vort
!         Solution: ipol%i_velo, itor%i_velo, idpdr%i_velo
        if (iflag_debug .gt. 0)                                         &
     &       write(*,*) 'cal_sol_velo_by_vort_sph_crank'
        call cal_sol_velo_by_vort_sph_crank                             &
     &     (rj_fld%ntot_phys, rj_fld%d_fld)
        call const_grad_vp_and_vorticity                                &
     &   (ipol%i_velo, ipol%i_vort, rj_fld)
      end if
!
!  Input: ipol%i_temp,  Solution: ipol%i_temp
      if(iflag_debug.gt.0) write(*,*) 'cal_sol_temperature_sph_crank'
      if(iflag_t_evo_4_temp .gt. id_no_evolution) then
        call cal_sol_temperature_sph_crank                              &
     &     (rj_fld%ntot_phys, rj_fld%d_fld)
      end if
!
!  Input: ipol%i_light,  Solution: ipol%i_light
      if(iflag_debug.gt.0) write(*,*) 'cal_sol_composition_sph_crank'
      if(iflag_t_evo_4_composit .gt. id_no_evolution) then
        call cal_sol_composition_sph_crank                              &
     &     (rj_fld%ntot_phys, rj_fld%d_fld)
      end if
!
!  Input: ipol%i_magne, itor%i_magne
!  Solution: ipol%i_magne, itor%i_magne, idpdr%i_magne
      if(iflag_debug.gt.0) write(*,*) 'cal_sol_magne_sph_crank'
      if(iflag_t_evo_4_magne .gt. id_no_evolution) then
        call cal_sol_magne_sph_crank(rj_fld%ntot_phys, rj_fld%d_fld)
        call const_grad_bp_and_current                                  &
     &     (sph_bc_B, ipol%i_magne, ipol%i_current, rj_fld)
      end if
!
!*  ---- update after evolution ------------------
!      call check_vs_spectr(rj_fld)
!
      if(iflag_t_evo_4_velo .gt. id_no_evolution) then
        call update_after_vorticity_sph(rj_fld)
        call cal_rot_radial_self_gravity(sph_bc_U, rj_fld)
      end if
!
      if(iflag_t_evo_4_temp .gt.     id_no_evolution) then
        call update_after_heat_sph(rj_fld)
      end if
      if(iflag_t_evo_4_composit .gt. id_no_evolution) then
        call update_after_composit_sph(rj_fld)
      end if
      if(iflag_t_evo_4_magne .gt.    id_no_evolution) then
        call update_after_magne_sph(rj_fld)
      end if
!
      end subroutine s_cal_sol_sph_MHD_crank
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine set_sph_field_to_start(rj_fld)
!
      use m_boundary_params_sph_MHD
      use const_sph_radial_grad
      use cal_rot_buoyancies_sph_MHD
!
      type(phys_data), intent(inout) :: rj_fld
!
!
      if(ipol%i_velo*ipol%i_vort .gt. 0) then
        call const_grad_vp_and_vorticity                                &
     &     (ipol%i_velo, ipol%i_vort, rj_fld)
      end if
!
      if(iflag_t_evo_4_velo .gt. id_no_evolution) then
        if(iflag_debug.gt.0) write(*,*) 'update_after_vorticity_sph'
        call update_after_vorticity_sph(rj_fld)
        if(iflag_debug.gt.0) write(*,*) 'cal_rot_radial_self_gravity'
        call cal_rot_radial_self_gravity(sph_bc_U, rj_fld)
      end if
!
      if(iflag_debug.gt.0) write(*,*) 'update_after_heat_sph'
      call update_after_heat_sph(rj_fld)
      if(iflag_debug.gt.0) write(*,*) 'update_after_composit_sph'
      call update_after_composit_sph(rj_fld)
!
      if(ipol%i_magne*ipol%i_current .gt. 0) then
        call const_grad_bp_and_current                                  &
     &     (sph_bc_B, ipol%i_magne, ipol%i_current, rj_fld)
      end if
!
      call update_after_magne_sph(rj_fld)
!
      end subroutine set_sph_field_to_start
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine update_after_vorticity_sph(rj_fld)
!
      use m_physical_property
      use m_boundary_params_sph_MHD
      use cal_inner_core_rotation
!
      type(phys_data), intent(inout) :: rj_fld
!
!
      if(sph_bc_U%iflag_icb .eq. iflag_rotatable_ic) then
        call set_inner_core_rotation(sph_bc_U%kr_in, rj_fld)
      end if
!
!       Input: ipol%i_vort, itor%i_vort
!       Solution: ipol%i_v_diffuse, itor%i_v_diffuse, idpdr%i_v_diffuse
      if(ipol%i_v_diffuse .gt. 0) then
        if(iflag_debug.gt.0) write(*,*) 'const_sph_viscous_by_vort2'
        call const_sph_viscous_by_vort2(sph_bc_U, coef_d_velo,          &
     &     ipol%i_velo, ipol%i_vort, ipol%i_v_diffuse, rj_fld)
      end if
!
!       Input:    ipol%i_vort, itor%i_vort
!       Solution: ipol%i_w_diffuse, itor%i_w_diffuse, idpdr%i_w_diffuse
      if(ipol%i_w_diffuse .gt. 0) then
        if(iflag_debug.gt.0) write(*,*)'const_sph_vorticirty_diffusion'
        call const_sph_vorticirty_diffusion(sph_bc_U, coef_d_velo,      &
     &      ipol%i_vort, ipol%i_w_diffuse, rj_fld)
      end if
!
      end subroutine update_after_vorticity_sph
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine update_after_magne_sph(rj_fld)
!
      use m_physical_property
      use m_boundary_params_sph_MHD
!
      type(phys_data), intent(inout) :: rj_fld
!
!
!       Input:    ipol%i_current, itor%i_current
!       Solution: ipol%i_b_diffuse, itor%i_b_diffuse, idpdr%i_b_diffuse
      if(ipol%i_b_diffuse .gt. 0) then
        if(iflag_debug .gt. 0) write(*,*) 'const_sph_mag_diffuse_by_j'
        call const_sph_mag_diffuse_by_j(sph_bc_B, coef_d_magne,         &
     &      ipol%i_magne, ipol%i_current, ipol%i_b_diffuse, rj_fld)
      end if
!
      end subroutine update_after_magne_sph
!
! -----------------------------------------------------------------------
!
      subroutine update_after_heat_sph(rj_fld)
!
      use m_boundary_params_sph_MHD
      use m_physical_property
!
      type(phys_data), intent(inout) :: rj_fld
!
!         Input: ipol%i_temp,  Solution: ipol%i_grad_t
      if(iflag_debug .gt. 0)  write(*,*)                                &
     &           'const_radial_grad_temp', ipol%i_grad_t
      if(ipol%i_grad_t .gt. 0) then
        call const_radial_grad_scalar                                   &
     &     (sph_bc_T, ipol%i_temp, ipol%i_grad_t, rj_fld)
      end if
!
!         Input: ipol%i_temp,  Solution: ipol%i_t_diffuse
      if(ipol%i_t_diffuse .gt. 0) then
        if(iflag_debug .gt. 0)  write(*,*)                              &
     &           'const_sph_scalar_diffusion', ipol%i_t_diffuse
        call const_sph_scalar_diffusion(sph_bc_T, coef_d_temp,          &
     &      ipol%i_temp, ipol%i_t_diffuse, rj_fld)
      end if
!
      end subroutine update_after_heat_sph
!
! -----------------------------------------------------------------------
!
      subroutine update_after_composit_sph(rj_fld)
!
      use m_boundary_params_sph_MHD
      use m_physical_property
!
      type(phys_data), intent(inout) :: rj_fld
!
!
!         Input: ipol%i_light,  Solution: ipol%i_grad_composit
      if(ipol%i_grad_composit .gt. 0) then
        call const_radial_grad_scalar(sph_bc_C,                         &
     &     ipol%i_light, ipol%i_grad_composit, rj_fld)
      end if
!
!         Input: ipol%i_light,  Solution: ipol%i_c_diffuse
      if(ipol%i_c_diffuse .gt. 0) then
        if(iflag_debug .gt. 0)  write(*,*)                              &
     &           'const_sph_scalar_diffusion', ipol%i_c_diffuse
        call const_sph_scalar_diffusion(sph_bc_C, coef_d_light,         &
     &      ipol%i_light, ipol%i_c_diffuse, rj_fld)
      end if
!
      end subroutine update_after_composit_sph
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine check_vs_spectr(rj_fld)
!
      type(phys_data), intent(in) :: rj_fld
!
      integer(kind = kint) :: j, k, inod
!
      write(150+my_rank,*) 'j, k, s_velo, ds_velo, t_velo, w_diffuse'
      do j = 1, nidx_rj(2)
         do k = 1, nidx_rj(1)
          inod = j + (k-1) * nidx_rj(2)
          write(150+my_rank,'(2i16,1p20E25.15e3)') j, k,                &
     &        rj_fld%d_fld(inod,ipol%i_velo),                           &
     &        rj_fld%d_fld(inod,idpdr%i_velo),                          &
     &        rj_fld%d_fld(inod,itor%i_velo)
        end do
      end do
!
      end subroutine check_vs_spectr
!
! -----------------------------------------------------------------------
!
      subroutine check_ws_spectr(rj_fld)
!
      type(phys_data), intent(in) :: rj_fld
!
      integer(kind = kint) :: j, k, inod
!
      write(150+my_rank,*) 'j, k, s_vort, ds_vort, t_vort'
      do j = 1, nidx_rj(2)
         do k = 1, nidx_rj(1)
          inod = j + (k-1) * nidx_rj(2)
          write(150+my_rank,'(2i16,1p20E25.15e3)') j, k,                &
     &        rj_fld%d_fld(inod,ipol%i_vort),                           &
     &        rj_fld%d_fld(inod,idpdr%i_vort),                          &
     &        rj_fld%d_fld(inod,itor%i_vort)
        end do
      end do
!
      end subroutine check_ws_spectr
!
! -----------------------------------------------------------------------
!
      end module cal_sol_sph_MHD_crank
