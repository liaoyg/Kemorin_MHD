!>@file   m_neutral_pt_by_pick_sph.f90
!!@brief      module m_neutral_pt_by_pick_sph
!!
!!@author H. Matsui and H. Okuda
!!@date Programmed in  Dec., 2012
!
!> @brief choose spectr data to output
!!
!!@verbatim
!!      subroutine alloc_neutral_point(num_layer)
!!      subroutine dealloc_neutral_point
!!      subroutine find_field_address(picked)
!!      subroutine set_radius_for_fdm(picked, sph_params, sph_rj, r_2nd)
!!        type(sph_shell_parameters), intent(inout) :: sph_params
!!        type(sph_rj_grid), intent(inout) ::  sph_rj
!!        type(fdm_matrices), intent(inout) :: r_2nd
!!      subroutine set_radial_grad_scalars(istep, time,                 &
!!     &          nri, radius_1d_rj_r, d1nod_mat_fdm_2, buo_ratio,      &
!!     &          picked)
!!@endverbatim
!
      module m_neutral_pt_by_pick_sph
!
      use m_precision
      use m_constants
      use m_phys_labels
!
      use t_pickup_sph_spectr_data
!
      use set_radius_func_noequi
!
      implicit  none
!
      integer(kind = kint), parameter :: id_neutral_pt = 98
      integer(kind = kint), parameter :: id_ave_den = 99
      character(len=kchara), parameter                                  &
     &               :: fname_neutral_pt = "neutral_point.dat"
      character(len=kchara), parameter                                  &
     &               :: fname_ave_den = "ave_density.dat"
!
      real(kind = kreal), allocatable :: temp00(:)
      real(kind = kreal), allocatable :: comp00(:)
      real(kind = kreal), allocatable :: grad_temp00(:)
      real(kind = kreal), allocatable :: grad_comp00(:)
      real(kind = kreal), allocatable :: freq(:)
      real(kind = kreal), allocatable :: freq2(:)
!
      integer(kind = kint) :: icomp_temp, icomp_light, ipick_l0m0
!
      private :: id_neutral_pt, fname_neutral_pt
      private :: id_ave_den, fname_ave_den
      private :: temp00, comp00, grad_temp00, grad_comp00, freq, freq2
      private :: icomp_temp, icomp_light, ipick_l0m0
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine alloc_neutral_point(num_layer)
!
      integer(kind = kint), intent(in) :: num_layer
!
      allocate(temp00(num_layer))
      allocate(comp00(num_layer))
      allocate(grad_temp00(num_layer))
      allocate(grad_comp00(num_layer))
      allocate(freq(num_layer))
      allocate(freq2(num_layer))
      temp00 = 0.0d0
      comp00 = 0.0d0
      grad_temp00 = 0.0d0
      grad_comp00 = 0.0d0
      freq =  0.0d0
      freq2 = 0.0d0
!
      open(id_neutral_pt, file=fname_neutral_pt)
      open(id_ave_den, file=fname_ave_den)
!
      write(id_neutral_pt,'(a)') 'time_step  time  neutral_radius'
      write(id_ave_den,'(3a)') 'time_step  time  radius_ID  radius  ',  &
     &                      '  temperature  composition  heat_flux',    &
     &                      '   comp_flux  density_flux  freq'
!
      end subroutine alloc_neutral_point
!
! ----------------------------------------------------------------------
!
      subroutine dealloc_neutral_point
!
      deallocate(temp00, grad_temp00)
      deallocate(comp00, grad_comp00, freq, freq2)
!
      close(id_neutral_pt)
      close(id_ave_den)
!
      end subroutine dealloc_neutral_point
!
! ----------------------------------------------------------------------
!
      subroutine find_field_address(picked)
!
      type(picked_spectrum_data), intent(in) :: picked
!
      integer(kind = kint) :: i
!
      do i = 1, picked%ntot_comp_rj
        if(picked%spectr_name(i) .eq. fhd_temp)  icomp_temp =  i
        if(picked%spectr_name(i) .eq. fhd_light) icomp_light = i
      end do
!
      end subroutine find_field_address
!
! ----------------------------------------------------------------------
!
      subroutine set_radius_for_fdm(picked, sph_params, sph_rj, r_2nd)
!
      use t_spheric_parameter
      use t_fdm_coefs
      use const_fdm_coefs
!
      type(picked_spectrum_data), intent(in) :: picked
!
      type(sph_shell_parameters), intent(inout) :: sph_params
      type(sph_rj_grid), intent(inout) ::  sph_rj
      type(fdm_matrices), intent(inout) :: r_2nd
!
      integer(kind = kint) :: i
!
!
      sph_rj%nidx_rj(1) = picked%num_layer
      sph_rj%nidx_rj(2) = 1
      call alloc_type_sph_1d_index_rj(sph_rj)
!
      do i = 1, picked%num_layer
        sph_rj%radius_1d_rj_r(i) = picked%radius_gl(i)
      end do
      do i = 1, picked%num_sph_mode
        if(picked%idx_gl(i,1) .eq. 0) ipick_l0m0 = i
      end do
!
      call allocate_dr_rj_noequi(sph_rj%nidx_rj(1))
      call set_dr_for_nonequi(sph_params%nlayer_CMB,                    &
     &   sph_rj%nidx_rj(1), sph_rj%radius_1d_rj_r)
      call const_2nd_fdm_matrices(sph_params, sph_rj, r_2nd)
!
      write(*,*) 'icomp_temp, icomp_light', icomp_temp, icomp_light
      write(*,*) 'ipick_l0m0', ipick_l0m0
!
      end subroutine set_radius_for_fdm
!
! ----------------------------------------------------------------------
!
      subroutine set_radial_grad_scalars(istep, time,                   &
     &          nri, radius_1d_rj_r, d1nod_mat_fdm_2, buo_ratio,        &
     &          picked)
!
      integer(kind = kint), intent(in) :: istep
      real(kind = kreal), intent(in) :: time, buo_ratio
!
      integer(kind = kint), intent(in) :: nri
      real(kind = kreal), intent(in) :: radius_1d_rj_r(nri)
      real(kind = kreal), intent(in) :: d1nod_mat_fdm_2(nri,-1:1)
      type(picked_spectrum_data), intent(in) :: picked
!
      integer(kind = kint) :: k, ipick
!
      real(kind = kreal) :: r_neut
!
      do k = 1, picked%num_layer
        ipick = k + (ipick_l0m0-1) * picked%num_layer
        temp00(k) = picked%d_rj_gl(icomp_temp, ipick)
        comp00(k) = picked%d_rj_gl(icomp_light,ipick)
      end do
!
      do k = 2, picked%num_layer - 1
        grad_temp00(k) =  d1nod_mat_fdm_2(k,-1) * temp00(k-1)           &
     &                  + d1nod_mat_fdm_2(k, 0) * temp00(k  )           &
     &                  + d1nod_mat_fdm_2(k, 1) * temp00(k+1)
        grad_comp00(k) =  d1nod_mat_fdm_2(k,-1) * comp00(k-1)           &
     &                  + d1nod_mat_fdm_2(k, 0) * comp00(k  )           &
     &                  + d1nod_mat_fdm_2(k, 1) * comp00(k+1)
!
        freq2(k) = buo_ratio * grad_comp00(k) + grad_temp00(k)
        if(freq2(k) .gt. 0.0d0) freq(k) = sqrt(freq2(k))
        freq2(k) = freq2(k) * radius_1d_rj_r(k  )**2
      end do
!
      do k = picked%num_layer - 2, 2, - 1
        if( freq2(k).lt.0.0d0 .and. freq2(k+1).ge.0.0d0) then
          r_neut = (radius_1d_rj_r(k  )*abs(freq2(k+1))                 &
     &            + radius_1d_rj_r(k+1)*abs(freq2(k)  ) )               &
     &             / (abs(freq2(k+1) - freq2(k)))
          write(id_neutral_pt,'(i15,1p2E25.15e3)') istep, time, r_neut
        end if
      end do
!
      do k = 1, picked%num_layer
            write(id_ave_den,'(i15,1pE25.15e3,i15,1p7E25.15e3)')        &
     &           istep, time, k, radius_1d_rj_r(k),                     &
     &           temp00(k), comp00(k), grad_temp00(k), grad_comp00(k),  &
     &           freq2(k), freq(k)
      end do
!
       end subroutine set_radial_grad_scalars
!
! ----------------------------------------------------------------------
!
      end module m_neutral_pt_by_pick_sph
