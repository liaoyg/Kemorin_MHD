!const_sph_diffusion.f90
!      module const_sph_diffusion
!
!      modified by H. Matsui on Oct., 2009
!
!
!      subroutine const_sph_viscous_diffusion
!        Input:    ipol%i_velo, itor%i_velo
!        Solution: ipol%i_v_diffuse, itor%i_v_diffuse, idpdr%i_v_diffuse
!      subroutine const_sph_vorticirty_diffusion
!        Input:    ipol%i_vort, itor%i_vort
!        Solution: ipol%i_w_diffuse, itor%i_w_diffuse, idpdr%i_w_diffuse
!      subroutine const_sph_magnetic_diffusion
!        Input:    ipol%i_magne, itor%i_magne
!        Solution: ipol%i_b_diffuse, itor%i_b_diffuse, idpdr%i_b_diffuse
!      subroutine const_sph_thermal_diffusion
!      subroutine const_sph_dscalar_diffusion
!
      module const_sph_diffusion
!
      use m_precision
!
      use m_constants
      use m_spheric_parameter
      use m_sph_spectr_data
      use m_sph_phys_address
      use m_control_params_sph_MHD
      use cal_sph_exp_diffusion
      use cal_sph_exp_1st_diff
!
      implicit none
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine const_sph_viscous_diffusion
!
      use set_sph_exp_rigid_ICB
      use set_sph_exp_rigid_CMB
      use set_sph_exp_free_ICB
      use set_sph_exp_free_CMB
      use cal_sph_exp_fixed_scalar
!
      integer(kind = kint) :: kr_st, kr_ed
!
!
      kr_st = nlayer_ICB+1
      kr_ed = nlayer_CMB-1
      call cal_sph_nod_vect_diffuse2(kr_st, kr_ed, ipol%i_velo,         &
     &    ipol%i_v_diffuse)
      call cal_sph_nod_vect_dr_2(kr_st, kr_ed,                          &
     &    d_rj(1,ipol%i_v_diffuse), d_rj(1,idpdr%i_v_diffuse) )
!
      if(iflag_icb_velocity .eq. iflag_free_slip) then
        call cal_sph_nod_cmb_rigid_diffuse2(ipol%i_velo,                &
     &      ipol%i_v_diffuse)
      else
        call cal_sph_nod_icb_rigid_diffuse2(ipol%i_velo,                &
     &      ipol%i_v_diffuse)
      end if
!
      call cal_dsdr_sph_icb_nobc_2(ipol%i_v_diffuse, idpdr%i_v_diffuse)
!
      if(iflag_cmb_velocity .eq. iflag_free_slip) then
        call cal_sph_nod_cmb_free_diffuse2(ipol%i_velo,                 &
     &      ipol%i_v_diffuse)
      else
        call cal_sph_nod_cmb_rigid_diffuse2(ipol%i_velo,                &
     &      ipol%i_v_diffuse)
      end if
!
      call cal_dsdr_sph_cmb_nobc_2(ipol%i_v_diffuse, idpdr%i_v_diffuse)
!
      end subroutine const_sph_viscous_diffusion
!
! -----------------------------------------------------------------------
!
      subroutine const_sph_vorticirty_diffusion
!
      use set_sph_exp_rigid_ICB
      use set_sph_exp_rigid_CMB
      use set_sph_exp_free_ICB
      use set_sph_exp_free_CMB
      use cal_sph_exp_fixed_scalar
      use cal_inner_core_rotation
!
      integer(kind = kint) :: kr_st, kr_ed
!
!
      kr_st = nlayer_ICB+1
      kr_ed = nlayer_CMB-1
      call cal_sph_nod_vect_diffuse2(kr_st, kr_ed, ipol%i_vort,         &
     &    ipol%i_w_diffuse)
      call cal_sph_nod_vect_dr_2(kr_st, kr_ed,                          &
     &    d_rj(1,ipol%i_w_diffuse), d_rj(1,idpdr%i_w_diffuse) )
!
!
      if(iflag_icb_velocity .eq. iflag_free_slip) then
        call cal_sph_nod_icb_free_w_diffuse2(ipol%i_vort,               &
     &      ipol%i_w_diffuse)
      else
        call cal_sph_nod_icb_rgd_w_diffuse2(ipol%i_vort,                &
     &      ipol%i_w_diffuse)
      end if
!
      if(iflag_icb_velocity .eq. iflag_rotatable_ic) then
        call cal_icore_viscous_drag_explicit
      end if
!
      call cal_dsdr_sph_icb_nobc_2(ipol%i_w_diffuse, idpdr%i_w_diffuse)
!
      if(iflag_cmb_velocity .eq. iflag_free_slip) then
        call cal_sph_nod_cmb_free_w_diffuse2(ipol%i_vort,               &
     &      ipol%i_w_diffuse)
      else
        call cal_sph_nod_cmb_rgd_w_diffuse2(ipol%i_vort,                &
     &      ipol%i_w_diffuse)
      end if
!
      call cal_dsdr_sph_cmb_nobc_2(ipol%i_w_diffuse, idpdr%i_w_diffuse)
!
      end subroutine const_sph_vorticirty_diffusion
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine const_sph_magnetic_diffusion
!
      use cal_sph_exp_nod_icb_ins
      use cal_sph_exp_nod_cmb_ins
      use cal_sph_exp_nod_icb_qvac
      use cal_sph_exp_nod_cmb_qvac
      use set_sph_exp_nod_center
      use cal_sph_exp_fixed_scalar
!
      integer(kind = kint) :: kr_st, kr_ed
!
!
      if(iflag_icb_magne .eq. iflag_sph_fill_center) then
        kr_st = itwo
        call cal_sph_nod_center_diffuse2(ipol%i_magne,                  &
     &      ipol%i_b_diffuse)
        call cal_dsdr_sph_center_2(ipol%i_b_diffuse)
      else if(iflag_icb_magne .eq. iflag_pseudo_vacuum) then
        kr_st = nlayer_ICB+1
        call cal_sph_nod_icb_qvc_diffuse2(ipol%i_magne,                 &
     &      ipol%i_b_diffuse)
        call cal_dsdr_sph_icb_nobc_2(ipol%i_b_diffuse,                  &
     &      idpdr%i_b_diffuse)
      else
        kr_st = nlayer_ICB+1
        call cal_sph_nod_icb_ins_diffuse2(ipol%i_magne,                 &
     &      ipol%i_b_diffuse)
        call cal_dsdr_sph_icb_nobc_2(ipol%i_b_diffuse,                  &
     &      idpdr%i_b_diffuse)
      end if
!
!
      kr_ed = nlayer_CMB-1
      call cal_sph_nod_vect_diffuse2(kr_st, kr_ed, ipol%i_magne,        &
     &    ipol%i_b_diffuse)
      call cal_sph_nod_vect_dr_2(kr_st, kr_ed,                          &
     &    d_rj(1,ipol%i_b_diffuse), d_rj(1,idpdr%i_b_diffuse) )
!
      if(iflag_cmb_magne .eq. iflag_pseudo_vacuum) then
        call cal_sph_nod_cmb_qvc_diffuse2(ipol%i_magne,                 &
     &      ipol%i_b_diffuse)
      else
        call cal_sph_nod_cmb_ins_diffuse2(ipol%i_magne,                 &
     &      ipol%i_b_diffuse)
      end if
      call cal_dsdr_sph_cmb_nobc_2(ipol%i_b_diffuse, idpdr%i_b_diffuse)
!
      end subroutine const_sph_magnetic_diffusion
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine const_sph_thermal_diffusion
!
      use m_machine_parameter
      use m_control_params_sph_MHD
      use cal_sph_exp_fixed_scalar
      use cal_sph_exp_fixed_flux
!
      integer(kind = kint) :: kr_st, kr_ed
!
!
      kr_st = nlayer_ICB+1
      kr_ed = nlayer_CMB-1
      call cal_sph_nod_scalar_diffuse2(kr_st, kr_ed,                    &
     &    ipol%i_temp, ipol%i_t_diffuse)
!
      if (iflag_icb_temp .eq. iflag_fixed_flux) then
        call cal_sph_icb_fix_flux_diffuse2(nidx_rj(2), h_flux_ICB_bc,   &
     &      ipol%i_temp, ipol%i_t_diffuse)
      else
        call cal_sph_icb_fix_scalar_diffuse2(nidx_rj(2), temp_ICB_bc,   &
     &     ipol%i_temp, ipol%i_t_diffuse)
      end if
!
      if (iflag_cmb_temp .eq. iflag_fixed_flux) then
        call cal_sph_cmb_fix_flux_diffuse2(nidx_rj(2), h_flux_CMB_bc,   &
     &      ipol%i_temp, ipol%i_t_diffuse)
      else
        call cal_sph_cmb_fix_scalar_diffuse2(nidx_rj(2), temp_CMB_bc,   &
     &     ipol%i_temp, ipol%i_t_diffuse)
      end if
!
      end subroutine const_sph_thermal_diffusion
!
! -----------------------------------------------------------------------
!
      subroutine const_sph_dscalar_diffusion
!
      use m_control_params_sph_MHD
      use cal_sph_exp_fixed_scalar
      use cal_sph_exp_fixed_flux
!
      integer(kind = kint) :: kr_st, kr_ed
!
!
      kr_st = nlayer_ICB+1
      kr_ed = nlayer_CMB-1
      call cal_sph_nod_scalar_diffuse2(kr_st, kr_ed, ipol%i_light,      &
     &    ipol%i_c_diffuse)
!
      if (iflag_icb_composition .eq. iflag_fixed_flux) then
        call cal_sph_icb_fix_flux_diffuse2(nidx_rj(2), c_flux_ICB_bc,   &
     &      ipol%i_light, ipol%i_c_diffuse)
      else
        call cal_sph_icb_fix_scalar_diffuse2(nidx_rj(2),                &
     &      composition_ICB_bc, ipol%i_light, ipol%i_c_diffuse)
      end if
!
      if (iflag_cmb_composition .eq. iflag_fixed_flux) then
        call cal_sph_cmb_fix_flux_diffuse2(nidx_rj(2), c_flux_CMB_bc,   &
     &      ipol%i_light, ipol%i_c_diffuse)
      else
        call cal_sph_cmb_fix_scalar_diffuse2(nidx_rj(2),                &
     &      composition_CMB_bc, ipol%i_light, ipol%i_c_diffuse)
      end if
!
      end subroutine const_sph_dscalar_diffusion
!
! -----------------------------------------------------------------------
!
      end module const_sph_diffusion
