!>@file   set_address_sph_trans_MHD.f90
!!@brief  module set_address_sph_trans_MHD
!!
!!@author H. Matsui
!!@date Programmed in Jan., 2010
!
!>@brief Field addresses for spherical harmonics transform
!!       in MHD dynamo simulation
!!
!!@verbatim
!!      subroutine set_addresses_trans_sph_MHD(ipol, trns_MHD,          &
!!     &          ncomp_sph_trans, nvector_sph_trans, nscalar_sph_trans)
!!        type(phys_address), intent(in) :: ipol
!!        type(address_4_sph_trans), intent(inout) :: trns_MHD
!!      subroutine check_address_trans_sph_MHD                          &
!!     &         (ipol, idpdr, itor, trns_MHD, ncomp_sph_trans)
!!        type(phys_address), intent(in) :: ipol, idpdr, itor
!!        type(address_4_sph_trans), intent(in) :: trns_MHD
!!@endverbatim
!
      module set_address_sph_trans_MHD
!
      use m_precision
!
      use t_phys_address
      use t_addresses_sph_transform
!
      implicit none
!
      private :: b_trans_address_vector_MHD
      private :: b_trans_address_scalar_MHD
      private :: f_trans_address_vector_MHD
      private :: f_trans_address_scalar_MHD
      private :: check_add_trans_sph_MHD
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine set_addresses_trans_sph_MHD(ipol, trns_MHD,            &
     &          ncomp_sph_trans, nvector_sph_trans, nscalar_sph_trans)
!
      use m_control_parameter
!
      type(phys_address), intent(in) :: ipol
      type(address_4_sph_trans), intent(inout) :: trns_MHD
      integer(kind = kint), intent(inout) :: ncomp_sph_trans
      integer(kind = kint), intent(inout) :: nvector_sph_trans
      integer(kind = kint), intent(inout) :: nscalar_sph_trans
!
      integer(kind = kint) :: nscltsr_rtp_2_rj, nscltsr_rj_2_rtp
!
      call b_trans_address_vector_MHD                                   &
     &   (ipol, trns_MHD%nvector_rj_2_rtp, trns_MHD%b_trns)
      call b_trans_address_scalar_MHD                                   &
     &   (ipol, trns_MHD%nvector_rj_2_rtp, trns_MHD%nscalar_rj_2_rtp,   &
     &    trns_MHD%b_trns)
      trns_MHD%ntensor_rj_2_rtp = 0
!
      call f_trans_address_vector_MHD                                   &
     &   (ipol, trns_MHD%nvector_rtp_2_rj, trns_MHD%f_trns)
      call f_trans_address_scalar_MHD                                   &
     &   (trns_MHD%nvector_rtp_2_rj, trns_MHD%nscalar_rtp_2_rj,         &
     &    trns_MHD%f_trns)
      trns_MHD%ntensor_rtp_2_rj = 0
!
      nscltsr_rtp_2_rj                                                  &
     &      = trns_MHD%nscalar_rj_2_rtp + 6*trns_MHD%ntensor_rj_2_rtp
      trns_MHD%ncomp_rj_2_rtp                                           &
     &      = 3*trns_MHD%nvector_rj_2_rtp + nscltsr_rtp_2_rj
!
      nscltsr_rj_2_rtp                                                  &
     &      = trns_MHD%nscalar_rtp_2_rj + 6*trns_MHD%ntensor_rtp_2_rj
      trns_MHD%ncomp_rtp_2_rj                                           &
     &      = 3*trns_MHD%nvector_rtp_2_rj + nscltsr_rj_2_rtp
!
      ncomp_sph_trans                                                   &
     &      = max(trns_MHD%ncomp_rj_2_rtp, trns_MHD%ncomp_rtp_2_rj)
      nvector_sph_trans                                                 &
     &      = max(trns_MHD%nvector_rj_2_rtp, trns_MHD%nvector_rtp_2_rj)
      nscalar_sph_trans = max(nscltsr_rtp_2_rj, nscltsr_rj_2_rtp)
!
      end subroutine set_addresses_trans_sph_MHD
!
!-----------------------------------------------------------------------
!
      subroutine check_address_trans_sph_MHD                            &
     &         (ipol, idpdr, itor, trns_MHD, ncomp_sph_trans)
!
      type(phys_address), intent(in) :: ipol, idpdr, itor
      type(address_4_sph_trans), intent(in) :: trns_MHD
      integer(kind = kint), intent(in) :: ncomp_sph_trans
!
!
      call check_add_trans_sph_MHD                                      &
     &   (ipol, idpdr, itor, trns_MHD%b_trns, trns_MHD%f_trns,          &
     &    trns_MHD%ncomp_rj_2_rtp, trns_MHD%nvector_rj_2_rtp,           &
     &    trns_MHD%nscalar_rj_2_rtp, trns_MHD%ncomp_rtp_2_rj,           &
     &    trns_MHD%nvector_rtp_2_rj, trns_MHD%nscalar_rtp_2_rj,         &
     &    ncomp_sph_trans)
!
      end subroutine check_address_trans_sph_MHD
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine b_trans_address_vector_MHD                             &
     &         (ipol, nvector_rj_2_rtp, b_trns)
!
      use m_control_parameter
!
      type(phys_address), intent(in) :: ipol
      integer(kind = kint), intent(inout) :: nvector_rj_2_rtp
      type(phys_address), intent(inout) :: b_trns
!
      nvector_rj_2_rtp = 0
!   velocity flag
      if(       iflag_t_evo_4_velo .gt.     id_no_evolution             &
     &     .or. iflag_t_evo_4_magne .gt.    id_no_evolution             &
     &     .or. iflag_t_evo_4_temp .gt.     id_no_evolution             &
     &     .or. iflag_t_evo_4_composit .gt. id_no_evolution) then
        nvector_rj_2_rtp = nvector_rj_2_rtp + 1
        b_trns%i_velo = 3*nvector_rj_2_rtp - 2
      end if
!   vorticity flag
      if(       iflag_t_evo_4_velo .gt. id_no_evolution) then
        nvector_rj_2_rtp = nvector_rj_2_rtp + 1
        b_trns%i_vort = 3*nvector_rj_2_rtp - 2
      end if
!   magnetic field flag
      if(       iflag_t_evo_4_magne .gt. id_no_evolution                &
     &     .or. iflag_4_lorentz .gt.     id_turn_OFF) then
        nvector_rj_2_rtp = nvector_rj_2_rtp + 1
        b_trns%i_magne = 3*nvector_rj_2_rtp - 2
      end if
!   current density flag
      if(iflag_4_lorentz .gt. id_turn_OFF) then
        nvector_rj_2_rtp = nvector_rj_2_rtp + 1
        b_trns%i_current = 3*nvector_rj_2_rtp - 2
      end if
!
!
!   filtered velocity
      if(ipol%i_filter_velo .gt. izero) then
        nvector_rj_2_rtp = nvector_rj_2_rtp + 1
        b_trns%i_filter_velo =    3*nvector_rj_2_rtp - 2
      end if
!   filtered vorticity
      if(ipol%i_filter_vort .gt. izero) then
        nvector_rj_2_rtp = nvector_rj_2_rtp + 1
        b_trns%i_filter_vort =    3*nvector_rj_2_rtp - 2
      end if
!   filtered magnetic field
      if(ipol%i_filter_magne .gt. izero) then
        nvector_rj_2_rtp = nvector_rj_2_rtp + 1
        b_trns%i_filter_magne =   3*nvector_rj_2_rtp - 2
      end if
!   filtered current density
      if(ipol%i_filter_current .gt. izero) then
        nvector_rj_2_rtp = nvector_rj_2_rtp + 1
        b_trns%i_filter_current = 3*nvector_rj_2_rtp - 2
      end if
!
!   dual filtered velocity
      if(ipol%i_wide_fil_velo .gt. izero) then
        nvector_rj_2_rtp = nvector_rj_2_rtp + 1
        b_trns%i_wide_fil_velo =    3*nvector_rj_2_rtp - 2
      end if
!   dual filtered vorticity
      if(ipol%i_wide_fil_vort .gt. izero) then
        nvector_rj_2_rtp = nvector_rj_2_rtp + 1
        b_trns%i_wide_fil_vort =    3*nvector_rj_2_rtp - 2
      end if
!   dual filtered magnetic field
      if(ipol%i_wide_fil_magne .gt. izero) then
        nvector_rj_2_rtp = nvector_rj_2_rtp + 1
        b_trns%i_wide_fil_magne =   3*nvector_rj_2_rtp - 2
      end if
!   dual filtered current density
      if(ipol%i_wide_fil_current .gt. izero) then
        nvector_rj_2_rtp = nvector_rj_2_rtp + 1
        b_trns%i_wide_fil_current = 3*nvector_rj_2_rtp - 2
      end if
!
      end subroutine b_trans_address_vector_MHD
!
!-----------------------------------------------------------------------
!
      subroutine b_trans_address_scalar_MHD                             &
     &         (ipol, nvector_rj_2_rtp, nscalar_rj_2_rtp, b_trns)
!
      use m_control_parameter
!
      type(phys_address), intent(in) :: ipol
      integer(kind = kint), intent(in) :: nvector_rj_2_rtp
      integer(kind = kint), intent(inout) :: nscalar_rj_2_rtp
      type(phys_address), intent(inout) :: b_trns
!
!
      nscalar_rj_2_rtp = 0
!   temperature flag
      if(iflag_t_evo_4_temp .gt. id_no_evolution) then
        nscalar_rj_2_rtp = nscalar_rj_2_rtp + 1
        b_trns%i_temp = nscalar_rj_2_rtp + 3*nvector_rj_2_rtp
      end if
!   composition flag
      if(iflag_t_evo_4_composit .gt. id_no_evolution) then
        nscalar_rj_2_rtp = nscalar_rj_2_rtp + 1
        b_trns%i_light = nscalar_rj_2_rtp + 3*nvector_rj_2_rtp
      end if
!
!   filtered temperature
      if(ipol%i_filter_temp .gt. izero) then
        nscalar_rj_2_rtp = nscalar_rj_2_rtp + 1
        b_trns%i_filter_temp = nscalar_rj_2_rtp + 3*nvector_rj_2_rtp
      end if
!   filtered composition
      if(ipol%i_filter_comp .gt. izero) then
        nscalar_rj_2_rtp = nscalar_rj_2_rtp + 1
        b_trns%i_filter_comp = nscalar_rj_2_rtp + 3*nvector_rj_2_rtp
      end if
!
!   dual filtered temperature
      if(ipol%i_wide_fil_temp .gt. izero) then
        nscalar_rj_2_rtp = nscalar_rj_2_rtp + 1
        b_trns%i_wide_fil_temp = nscalar_rj_2_rtp + 3*nvector_rj_2_rtp
      end if
!   dual filtered composition
      if(ipol%i_wide_fil_comp .gt. izero) then
        nscalar_rj_2_rtp = nscalar_rj_2_rtp + 1
        b_trns%i_wide_fil_comp = nscalar_rj_2_rtp + 3*nvector_rj_2_rtp
      end if
!
      end subroutine b_trans_address_scalar_MHD
!
!-----------------------------------------------------------------------
!
      subroutine f_trans_address_vector_MHD                             &
     &         (ipol, nvector_rtp_2_rj, f_trns)
!
      use m_control_parameter
!
      type(phys_address), intent(in) :: ipol
      type(phys_address), intent(inout) :: f_trns
      integer(kind = kint), intent(inout) :: nvector_rtp_2_rj
!
!
      nvector_rtp_2_rj = 0
!   advection flag
      if(iflag_t_evo_4_velo .gt. id_no_evolution) then
        nvector_rtp_2_rj = nvector_rtp_2_rj + 1
        f_trns%i_m_advect = 3*nvector_rtp_2_rj - 2
!   Coriolis flag
        if(iflag_4_coriolis .gt. id_turn_OFF) then
          nvector_rtp_2_rj = nvector_rtp_2_rj + 1
          f_trns%i_coriolis = 3*nvector_rtp_2_rj - 2
        end if
        if(iflag_4_coriolis .gt. id_turn_OFF) then
          nvector_rtp_2_rj =      nvector_rtp_2_rj + 1
          f_trns%i_rot_Coriolis = 3*nvector_rtp_2_rj - 2
        end if
!   Lorentz flag
        if(iflag_4_lorentz .gt. id_turn_OFF) then
          nvector_rtp_2_rj = nvector_rtp_2_rj + 1
          f_trns%i_lorentz = 3*nvector_rtp_2_rj - 2
        end if
      end if
!
!   induction flag
      if(iflag_t_evo_4_magne .gt. id_no_evolution) then
        nvector_rtp_2_rj = nvector_rtp_2_rj + 1
        f_trns%i_vp_induct =  3*nvector_rtp_2_rj - 2
      end if
!
!   heat flux flag
      if(iflag_t_evo_4_temp .gt. id_no_evolution) then
        nvector_rtp_2_rj = nvector_rtp_2_rj + 1
        f_trns%i_h_flux = 3*nvector_rtp_2_rj - 2
      end if
!
!   composition flux flag
      if(iflag_t_evo_4_composit .gt. id_no_evolution) then
        nvector_rtp_2_rj = nvector_rtp_2_rj + 1
        f_trns%i_c_flux = 3*nvector_rtp_2_rj - 2
      end if
!
!
!   filtered advection flag
      if(ipol%i_SGS_inertia .gt. izero) then
        nvector_rtp_2_rj = nvector_rtp_2_rj + 1
        f_trns%i_SGS_inertia = 3*nvector_rtp_2_rj - 2
      end if
!
!   filtered Lorentz force flag
      if(ipol%i_SGS_Lorentz .gt. izero) then
        nvector_rtp_2_rj = nvector_rtp_2_rj + 1
        f_trns%i_SGS_Lorentz = 3*nvector_rtp_2_rj - 2
      end if
!
!   filtered induction flag
      if(ipol%i_SGS_vp_induct .gt. izero) then
        nvector_rtp_2_rj = nvector_rtp_2_rj + 1
        f_trns%i_SGS_vp_induct = 3*nvector_rtp_2_rj - 2
      end if
!
!   filtered heat flux flag
      if(ipol%i_SGS_h_flux .gt. izero) then
        nvector_rtp_2_rj = nvector_rtp_2_rj + 1
        f_trns%i_SGS_h_flux = 3*nvector_rtp_2_rj - 2
      end if
!
!   filtered composition flux flag
      if(ipol%i_SGS_c_flux .gt. izero) then
        nvector_rtp_2_rj = nvector_rtp_2_rj + 1
        f_trns%i_SGS_c_flux = 3*nvector_rtp_2_rj - 2
      end if
!
      end subroutine f_trans_address_vector_MHD
!
!-----------------------------------------------------------------------
!
      subroutine f_trans_address_scalar_MHD(nvector_rtp_2_rj,           &
     &          nscalar_rtp_2_rj, f_trns)
!
      use m_control_parameter
!
      integer(kind = kint), intent(in) :: nvector_rtp_2_rj
      integer(kind = kint), intent(inout) :: nscalar_rtp_2_rj
      type(phys_address), intent(inout) :: f_trns
!
!
      nscalar_rtp_2_rj = 0
!   divergence of Coriolis flux flag
      if(iflag_4_coriolis .gt. id_turn_OFF) then
        nscalar_rtp_2_rj = nscalar_rtp_2_rj + 1
        f_trns%i_div_Coriolis = nscalar_rtp_2_rj + 3*nvector_rtp_2_rj
      end if
!
      end subroutine f_trans_address_scalar_MHD
!
!-----------------------------------------------------------------------
!
      subroutine check_add_trans_sph_MHD                                &
     &         (ipol, idpdr, itor, b_trns, f_trns,                      &
     &          ncomp_rj_2_rtp, nvector_rj_2_rtp, nscalar_rj_2_rtp,     &
     &          ncomp_rtp_2_rj, nvector_rtp_2_rj, nscalar_rtp_2_rj,     &
     &          ncomp_sph_trans)
!
      type(phys_address), intent(in) :: ipol, idpdr, itor
      type(phys_address), intent(in) :: b_trns, f_trns
      integer(kind = kint), intent(in) :: ncomp_rj_2_rtp
      integer(kind = kint), intent(in) :: nvector_rj_2_rtp
      integer(kind = kint), intent(in) :: nscalar_rj_2_rtp
      integer(kind = kint), intent(in) :: ncomp_rtp_2_rj
      integer(kind = kint), intent(in) :: nvector_rtp_2_rj
      integer(kind = kint), intent(in) :: nscalar_rtp_2_rj
      integer(kind = kint), intent(in) :: ncomp_sph_trans
!
!
      write(*,*) 'ncomp_sph_trans ', ncomp_sph_trans
      write(*,*) 'ncomp_rj_2_rtp  ', ncomp_rj_2_rtp
      write(*,*) 'ncomp_rtp_2_rj  ', ncomp_rtp_2_rj
!
      write(*,*) 'nvector_rj_2_rtp  ', nvector_rj_2_rtp
      if(b_trns%i_velo .gt. 0) write(*,*)                               &
     &        'b_trns%i_velo  ', b_trns%i_velo,                         &
     &        ipol%i_velo, itor%i_velo, idpdr%i_velo
      if(b_trns%i_vort .gt. 0) write(*,*)                               &
     &        'b_trns%i_vort  ', b_trns%i_vort,                         &
     &        ipol%i_vort, itor%i_vort, idpdr%i_vort
      if(b_trns%i_magne .gt. 0) write(*,*)                              &
     &        'b_trns%i_magne ', b_trns%i_magne,                        &
     &        ipol%i_magne, itor%i_magne, idpdr%i_magne
      if(b_trns%i_current .gt. 0) write(*,*)                            &
     &        'b_trns%i_current ', b_trns%i_current, ipol%i_current,    &
     &        itor%i_current, idpdr%i_current
!
      if(b_trns%i_filter_velo .gt. 0) write(*,*)                        &
     &        'b_trns%i_filter_velo  ', b_trns%i_filter_velo,           &
     &        ipol%i_filter_velo, itor%i_filter_velo,                   &
     &        idpdr%i_filter_velo
      if(b_trns%i_filter_vort .gt. 0) write(*,*)                        &
     &        'b_trns%i_filter_vort  ', b_trns%i_filter_vort,           &
     &        ipol%i_filter_vort, itor%i_filter_vort,                   &
     &        idpdr%i_filter_vort
      if(b_trns%i_filter_magne .gt. 0) write(*,*)                       &
     &        'b_trns%i_filter_magne ', b_trns%i_filter_magne,          &
     &        ipol%i_filter_magne, itor%i_filter_magne,                 &
     &       idpdr%i_filter_magne
      if(b_trns%i_filter_current .gt. 0) write(*,*)                     &
     &        'b_trns%i_filter_current ', b_trns%i_filter_current,      &
     &        ipol%i_filter_current, itor%i_filter_current,             &
     &        idpdr%i_filter_current
!
      if(b_trns%i_wide_fil_velo .gt. 0) write(*,*)                      &
     &        'b_trns%i_wide_fil_velo  ', b_trns%i_wide_fil_velo,       &
     &        ipol%i_wide_fil_velo, itor%i_wide_fil_velo,               &
     &        idpdr%i_wide_fil_velo
      if(b_trns%i_wide_fil_vort .gt. 0) write(*,*)                      &
     &        'b_trns%i_wide_fil_vort  ', b_trns%i_wide_fil_vort,       &
     &        ipol%i_wide_fil_vort, itor%i_wide_fil_vort,               &
     &        idpdr%i_wide_fil_vort
      if(b_trns%i_wide_fil_magne .gt. 0) write(*,*)                     &
     &        'b_trns%i_wide_fil_magne ', b_trns%i_wide_fil_magne,      &
     &        ipol%i_wide_fil_magne, itor%i_wide_fil_magne,             &
     &        idpdr%i_wide_fil_magne
      if(b_trns%i_wide_fil_current .gt. 0) write(*,*)                   &
     &        'b_trns%i_wide_fil_current ', b_trns%i_wide_fil_current,  &
      &       ipol%i_wide_fil_current, itor%i_wide_fil_current,         &
     &        idpdr%i_wide_fil_current
      write(*,*)
!
      write(*,*) 'nscalar_rj_2_rtp  ', nscalar_rj_2_rtp
      if(b_trns%i_temp .gt. 0) write(*,*) 'b_trns%i_temp   ',           &
     &         b_trns%i_temp, ipol%i_temp
      if(b_trns%i_light .gt. 0) write(*,*) 'b_trns%i_light  ',          &
     &         b_trns%i_light, ipol%i_light
!
      if(b_trns%i_filter_temp .gt. 0) write(*,*)                        &
     &         'b_trns%i_filter_temp   ',                               &
     &         b_trns%i_filter_temp, ipol%i_filter_temp
      if(b_trns%i_filter_comp .gt. 0) write(*,*)                        &
     &         'b_trns%i_filter_comp  ',                                &
     &         b_trns%i_filter_comp, ipol%i_filter_comp
!
      if(b_trns%i_wide_fil_temp .gt. 0) write(*,*)                      &
     &         'b_trns%i_wide_fil_temp   ',                             &
     &         b_trns%i_wide_fil_temp, ipol%i_wide_fil_temp
      if(b_trns%i_wide_fil_comp .gt. 0) write(*,*)                      &
     &         'b_trns%i_wide_fil_comp  ',                              &
     &         b_trns%i_wide_fil_comp, ipol%i_wide_fil_comp
      write(*,*)
!
      write(*,*) 'nvector_rtp_2_rj  ', nvector_rtp_2_rj
      if(f_trns%i_m_advect .gt. 0) write(*,*) 'f_trns%i_m_advect ',     &
     &        f_trns%i_m_advect, ipol%i_m_advect,                       &
     &        itor%i_m_advect, idpdr%i_m_advect
      if(f_trns%i_coriolis .gt. 0) write(*,*) 'f_trns%i_coriolis  ',    &
     &        f_trns%i_coriolis, ipol%i_coriolis,                       &
     &        itor%i_coriolis, idpdr%i_coriolis
      if(f_trns%i_rot_Coriolis .gt. 0) write(*,*)                       &
     &       'f_trns%i_rot_Coriolis  ', f_trns%i_rot_Coriolis,          &
     &        ipol%i_rot_Coriolis, itor%i_rot_Coriolis,                 &
     &        idpdr%i_rot_Coriolis
      if(f_trns%i_lorentz .gt. 0) write(*,*) 'f_trns%i_lorentz  ',      &
     &        f_trns%i_lorentz, ipol%i_lorentz,                         &
     &        itor%i_lorentz, idpdr%i_lorentz
      if(f_trns%i_vp_induct .gt. 0) write(*,*) 'f_trns%i_vp_induct ',   &
     &        f_trns%i_vp_induct, ipol%i_vp_induct,                     &
     &        itor%i_vp_induct, idpdr%i_vp_induct
      if(f_trns%i_h_flux .gt. 0) write(*,*) 'f_trns%i_h_flux',          &
     &        f_trns%i_h_flux, ipol%i_h_flux,                           &
     &        itor%i_h_flux, idpdr%i_h_flux
      if(f_trns%i_c_flux .gt. 0) write(*,*) 'f_trns%i_c_flux',          &
     &        f_trns%i_c_flux, ipol%i_c_flux,                           &
     &        itor%i_c_flux, idpdr%i_c_flux
!
      if(f_trns%i_SGS_inertia .gt. 0) write(*,*)                        &
     &        'f_trns%i_SGS_inertia ', f_trns%i_SGS_inertia,            &
     &        ipol%i_SGS_inertia, itor%i_SGS_inertia,                   &
     &        idpdr%i_SGS_inertia
      if(f_trns%i_SGS_Lorentz .gt. 0) write(*,*)                        &
     &        'f_trns%i_SGS_Lorentz  ', f_trns%i_SGS_Lorentz,           &
     &        ipol%i_SGS_Lorentz, itor%i_SGS_Lorentz,                   &
     &        idpdr%i_SGS_Lorentz
      if(f_trns%i_SGS_vp_induct .gt. 0) write(*,*)                      &
     &        'f_trns%i_SGS_vp_induct ', f_trns%i_SGS_vp_induct,        &
     &        ipol%i_SGS_vp_induct, itor%i_SGS_vp_induct,               &
     &        idpdr%i_SGS_vp_induct
      if(f_trns%i_SGS_h_flux .gt. 0) write(*,*) 'f_trns%i_SGS_h_flux',  &
     &        f_trns%i_SGS_h_flux, ipol%i_SGS_h_flux,                   &
     &        itor%i_SGS_h_flux,idpdr%i_SGS_h_flux
      if(f_trns%i_SGS_c_flux .gt. 0) write(*,*) 'f_trns%i_SGS_c_flux',  &
     &        f_trns%i_SGS_c_flux, ipol%i_SGS_c_flux,                   &
     &        itor%i_SGS_c_flux, idpdr%i_SGS_c_flux
!
!
      write(*,*) 'nscalar_rtp_2_rj  ', nscalar_rtp_2_rj
      if(f_trns%i_div_Coriolis .gt. 0) write(*,*)                       &
     &       'f_trns%i_div_Coriolis  ', f_trns%i_div_Coriolis,          &
     &        ipol%i_div_Coriolis
      write(*,*)
!
      end subroutine check_add_trans_sph_MHD
!
!-----------------------------------------------------------------------
!
      end module set_address_sph_trans_MHD
