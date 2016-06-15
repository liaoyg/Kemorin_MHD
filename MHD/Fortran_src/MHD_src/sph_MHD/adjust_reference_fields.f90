!>@file   adjust_reference_fields.f90
!!@brief  module adjust_reference_fields
!!
!!@author H. Matsui
!!@date    programmed by H.Matsui in Oct., 2015
!
!>@brief Set boundary conditions for MHD dynamo simulation
!!
!!@verbatim
!!      subroutine adjust_press_by_average_on_CMB                       &
!!     &         (kr_in, kr_out, sph_rj, ipol, rj_fld)
!!      subroutine sync_temp_by_per_temp_sph                            &
!!     &         (reftemp_rj, sph_rj, ipol, idpdr, rj_fld)
!!        d_rj(inod,ipol%i_temp):        T => \Theta = T - T0
!!        d_rj(inod,ipol%i_par_temp):    \Theta = T - T0
!!        d_rj(inod,ipol%i_grad_t):      T => d \Theta / dr
!!        d_rj(inod,ipol%i_grad_part_t): d \Theta / dr
!!      subroutine trans_per_temp_to_temp_sph                           &
!!     &         (reftemp_rj, sph_rj, ipol, idpdr, rj_fld)
!!        d_rj(inod,ipol%i_temp):        \Theta = T - T0 => T
!!        d_rj(inod,ipol%i_par_temp):    \Theta = T - T0
!!        d_rj(inod,ipol%i_grad_t):      d \Theta / dr   => dT / dr
!!        d_rj(inod,ipol%i_grad_part_t): d \Theta / dr
!!
!!        type(sph_shell_parameters), intent(in) :: sph_params
!!        type(sph_rj_grid), intent(in) ::  sph_rj
!!        type(phys_address), intent(in) :: ipol, idpdr
!!@endverbatim
!
      module adjust_reference_fields
!
      use m_precision
      use m_machine_parameter
!
      use t_spheric_parameter
      use t_spheric_rj_data
      use t_phys_data
      use t_phys_address
!
      implicit  none
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine adjust_press_by_average_on_CMB                         &
     &         (kr_in, kr_out, sph_rj, ipol, rj_fld)
!
      use set_reference_sph_mhd
!
      integer(kind = kint), intent(in) :: kr_in, kr_out
      type(sph_rj_grid), intent(in) ::  sph_rj
      type(phys_address), intent(in) :: ipol
!
      type(phys_data), intent(inout) :: rj_fld
!
!
      call adjust_by_ave_pressure_on_CMB(kr_in, kr_out,                 &
     &    sph_rj%idx_rj_degree_zero, sph_rj%nidx_rj, ipol%i_press,      &
     &    rj_fld%n_point, rj_fld%ntot_phys, rj_fld%d_fld)
!
      end subroutine adjust_press_by_average_on_CMB
!
! -----------------------------------------------------------------------
!
      subroutine sync_temp_by_per_temp_sph                              &
     &         (reftemp_rj, sph_rj, ipol, idpdr, rj_fld)
!
      use set_reference_sph_mhd
!
      type(sph_rj_grid), intent(in) ::  sph_rj
      type(phys_address), intent(in) :: ipol, idpdr
      real(kind=kreal), intent(in)                                      &
     &               :: reftemp_rj(sph_rj%nidx_rj(1),0:1)
!
      type(phys_data), intent(inout) :: rj_fld
!
!
      call chenge_temp_to_per_temp_sph(sph_rj%idx_rj_degree_zero,       &
     &    sph_rj%nidx_rj, sph_rj%radius_1d_rj_r, reftemp_rj,            &
     &    ipol, idpdr, rj_fld%n_point, rj_fld%ntot_phys, rj_fld%d_fld)
!
      end subroutine sync_temp_by_per_temp_sph
!
! -----------------------------------------------------------------------
!
      subroutine trans_per_temp_to_temp_sph                             &
     &         (reftemp_rj, sph_rj, ipol, idpdr, rj_fld)
!
      use set_reference_sph_mhd
!
      type(sph_rj_grid), intent(in) ::  sph_rj
      type(phys_address), intent(in) :: ipol, idpdr
      real(kind=kreal), intent(in)                                      &
     &        :: reftemp_rj(sph_rj%nidx_rj(1),0:1)
!
      type(phys_data), intent(inout) :: rj_fld
!
!
      call transfer_per_temp_to_temp_sph(sph_rj%idx_rj_degree_zero,     &
     &    sph_rj%nidx_rj, sph_rj%radius_1d_rj_r, reftemp_rj,            &
     &    ipol, idpdr, rj_fld%n_point, rj_fld%ntot_phys, rj_fld%d_fld)
!
      end subroutine trans_per_temp_to_temp_sph
!
! -----------------------------------------------------------------------
!
      end module adjust_reference_fields
