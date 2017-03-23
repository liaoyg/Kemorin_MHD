!> @file  cal_diff_adv_src_explicit.f90
!!      module cal_diff_adv_src_explicit
!!
!! @author  H. Matsui
!! @date Programmed in Oct. 2009
!
!> @brief Evaluate time evolution explicitly
!!
!!@verbatim
!!      subroutine scalar_diff_advect_adams(kr_st, kr_ed, jmax,         &
!!     &          ipol_diffuse, ipol_advect, ipol_scalar, ipol_pre,     &
!!     &          dt, coef_exp, n_point, ntot_phys_rj, d_rj)
!!      subroutine scalar_diff_advect_euler(kr_st, kr_ed, jmax,         &
!!     &          ipol_diffuse, ipol_advect, ipol_scalar, dt, coef_exp, &
!!     &          n_point, ntot_phys_rj, d_rj)
!!      subroutine set_ini_adams_scalar(kr_st, kr_ed, jmax,             &
!!     &          ipol_advect, ipol_pre, n_point, ntot_phys_rj, d_rj)
!!
!!      subroutine scalar_diff_adv_src_adams                            &
!!     &         (kr_st, kr_ed, jmax, inod_rj_center, ipol_diffuse,     &
!!     &          ipol_advect, ipol_source, ipol_scalar, ipol_pre, dt,  &
!!     &          coef_exp, coef_src, n_point, ntot_phys_rj, d_rj)
!!      subroutine scalar_diff_adv_src_euler                            &
!!     &         (kr_st, kr_ed, jmax, inod_rj_center,                   &
!!     &          ipol_diffuse, ipol_advect, ipol_source, ipol_scalar,  &
!!     &          dt, coef_exp, coef_src, n_point, ntot_phys_rj, d_rj)
!!      subroutine set_ini_adams_scalar_w_src                           &
!!     &         (kr_st, kr_ed, jmax, inod_rj_center,                   &
!!     &          ipol_advect, ipol_source, ipol_pre, coef_src,         &
!!     &          n_point, ntot_phys_rj, d_rj)
!!      subroutine scalar_stable_diff_src                               &
!!     &         (kr_st, kr_ed, jmax, inod_rj_center,                   &
!!     &          ipol_source, ipol_scalar, coef_src,                   &
!!     &          n_point, ntot_phys_rj, d_rj)
!!      subroutine scalar_stable_diffusion(kr_st, kr_ed, jmax,          &
!!     &          ipol_scalar, n_point, ntot_phys_rj, d_rj)
!!@endverbatim
!!
!!@param kr_st         Radial address for inner boundary
!!@param kr_ed         Radial address for outer boundary
!!@param ipol_diffuse  address for diffusion term
!!@param ipol_advect   address for advection term
!!@param ipol_source   address for source term
!!@param ipol_scalar   address for scalar field to update
!!@param ipol_pre      address for storeing previous evolution
!!@param coef_exp      coeefient for expilict evolution for diffusion
!!@param coef_src      coefficient for source term
!
      module cal_diff_adv_src_explicit
!
      use m_precision
      use m_constants
      use m_t_step_parameter
!
      implicit  none
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine scalar_diff_advect_adams(kr_st, kr_ed, jmax,           &
     &          ipol_diffuse, ipol_advect, ipol_scalar, ipol_pre,       &
     &          dt, coef_exp, n_point, ntot_phys_rj, d_rj)
!
      integer(kind = kint), intent(in) :: kr_st, kr_ed, jmax
      integer(kind = kint), intent(in) :: ipol_diffuse, ipol_advect
      integer(kind = kint), intent(in) :: ipol_scalar, ipol_pre
      integer(kind = kint), intent(in) :: n_point, ntot_phys_rj
      real(kind = kreal), intent(in) :: coef_exp
      real(kind = kreal), intent(in) :: dt
!
      real (kind=kreal), intent(inout) :: d_rj(n_point,ntot_phys_rj)
!
      integer(kind = kint) :: inod, ist, ied
!
!
      ist = (kr_st-1)*jmax + 1
      ied = kr_ed * jmax
!$omp do private (inod)
      do inod = ist, ied
        d_rj(inod,ipol_scalar) = d_rj(inod,ipol_scalar)                 &
     &          + dt * (coef_exp * d_rj(inod,ipol_diffuse)              &
     &              - adam_0 * d_rj(inod,ipol_advect)                   &
     &              + adam_1 * d_rj(inod,ipol_pre) )
!
         d_rj(inod,ipol_pre) = - d_rj(inod,ipol_advect)
       end do
!$omp end do
!
      end subroutine scalar_diff_advect_adams
!
! ----------------------------------------------------------------------
!
      subroutine scalar_diff_advect_euler(kr_st, kr_ed, jmax,           &
     &          ipol_diffuse, ipol_advect, ipol_scalar, dt, coef_exp,   &
     &          n_point, ntot_phys_rj, d_rj)
!
      integer(kind = kint), intent(in) :: kr_st, kr_ed, jmax
      integer(kind = kint), intent(in) :: ipol_diffuse, ipol_advect
      integer(kind = kint), intent(in) :: ipol_scalar
      integer(kind = kint), intent(in) :: n_point, ntot_phys_rj
      real(kind = kreal), intent(in) :: coef_exp
      real(kind = kreal), intent(in) :: dt
!
      real (kind=kreal), intent(inout) :: d_rj(n_point,ntot_phys_rj)
!
      integer(kind = kint) :: inod, ist, ied
!
!
      ist = (kr_st-1) * jmax + 1
      ied = kr_ed * jmax
!$omp do private (inod)
      do inod = ist, ied
        d_rj(inod,ipol_scalar) = d_rj(inod,ipol_scalar)                 &
     &         + dt * (coef_exp*d_rj(inod,ipol_diffuse)                 &
     &             - d_rj(inod,ipol_advect) )
       end do
!$omp end do
!
      end subroutine scalar_diff_advect_euler
!
! ----------------------------------------------------------------------
!
      subroutine set_ini_adams_scalar(kr_st, kr_ed, jmax,               &
     &          ipol_advect, ipol_pre, n_point, ntot_phys_rj, d_rj)
!
      integer(kind = kint), intent(in) :: kr_st, kr_ed, jmax
      integer(kind = kint), intent(in) :: ipol_advect, ipol_pre
      integer(kind = kint), intent(in) :: n_point, ntot_phys_rj
!
      real (kind=kreal), intent(inout) :: d_rj(n_point,ntot_phys_rj)
!
      integer(kind = kint) :: inod, ist, ied
!
!
      ist = (kr_st-1) * jmax + 1
      ied = kr_ed * jmax
!$omp do private (inod)
      do inod = ist, ied
         d_rj(inod,ipol_pre) =  -d_rj(inod,ipol_advect)
      end do
!$omp end do
!
      end subroutine set_ini_adams_scalar
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine scalar_diff_adv_src_adams                              &
     &         (kr_st, kr_ed, jmax, inod_rj_center, ipol_diffuse,       &
     &          ipol_advect, ipol_source, ipol_scalar, ipol_pre, dt,    &
     &          coef_exp, coef_src, n_point, ntot_phys_rj, d_rj)
!
      integer(kind = kint), intent(in) :: inod_rj_center
      integer(kind = kint), intent(in) :: kr_st, kr_ed, jmax
      integer(kind = kint), intent(in) :: ipol_diffuse, ipol_advect
      integer(kind = kint), intent(in) :: ipol_source
      integer(kind = kint), intent(in) :: ipol_scalar, ipol_pre
      integer(kind = kint), intent(in) :: n_point, ntot_phys_rj
      real(kind = kreal), intent(in) :: coef_exp, coef_src
      real(kind = kreal), intent(in) :: dt
!
      real (kind=kreal), intent(inout) :: d_rj(n_point,ntot_phys_rj)
!
      integer(kind = kint) :: inod, ist, ied
!
!
      ist = (kr_st-1) * jmax + 1
      ied = kr_ed * jmax
!$omp do private (inod)
      do inod = ist, ied
        d_rj(inod,ipol_scalar) = d_rj(inod,ipol_scalar)                 &
     &          + dt * (coef_exp * d_rj(inod,ipol_diffuse)              &
     &              + adam_0 * ( -d_rj(inod,ipol_advect)                &
     &                 + coef_src * d_rj(inod,ipol_source) )            &
     &              + adam_1 * d_rj(inod,ipol_pre) )
!
         d_rj(inod,ipol_pre) = - d_rj(inod,ipol_advect)                 &
     &               + coef_src * d_rj(inod,ipol_source)
       end do
!$omp end do
!
      if(inod_rj_center .eq. 0) return
        inod = inod_rj_center
        d_rj(inod,ipol_scalar) = d_rj(inod,ipol_scalar)                 &
     &          + dt * (coef_exp * d_rj(inod,ipol_diffuse)              &
     &              + adam_0 * ( -d_rj(inod,ipol_advect)                &
     &                 + coef_src * d_rj(inod,ipol_source) )            &
     &              + adam_1 * d_rj(inod,ipol_pre) )
!
         d_rj(inod,ipol_pre) = - d_rj(inod,ipol_advect)                 &
     &               + coef_src * d_rj(inod,ipol_source)
!
      end subroutine scalar_diff_adv_src_adams
!
! ----------------------------------------------------------------------
!
      subroutine scalar_diff_adv_src_euler                              &
     &         (kr_st, kr_ed, jmax, inod_rj_center,                     &
     &          ipol_diffuse, ipol_advect, ipol_source, ipol_scalar,    &
     &          dt, coef_exp, coef_src, n_point, ntot_phys_rj, d_rj)
!
      integer(kind = kint), intent(in) :: inod_rj_center
      integer(kind = kint), intent(in) :: kr_st, kr_ed, jmax
      integer(kind = kint), intent(in) :: ipol_diffuse, ipol_advect
      integer(kind = kint), intent(in) :: ipol_source
      integer(kind = kint), intent(in) :: ipol_scalar
      integer(kind = kint), intent(in) :: n_point, ntot_phys_rj
      real(kind = kreal), intent(in) :: coef_exp, coef_src
      real(kind = kreal), intent(in) :: dt
!
      real (kind=kreal), intent(inout) :: d_rj(n_point,ntot_phys_rj)
!
      integer(kind = kint) :: inod, ist, ied
!
!
      ist = (kr_st-1) * jmax + 1
      ied = kr_ed * jmax
!$omp do private (inod)
      do inod = ist, ied
        d_rj(inod,ipol_scalar) = d_rj(inod,ipol_scalar)                 &
     &         + dt * (coef_exp*d_rj(inod,ipol_diffuse)                 &
     &             - d_rj(inod,ipol_advect)                             &
     &              + coef_src * d_rj(inod,ipol_source) )
       end do
!$omp end do
!
      if(inod_rj_center .eq. 0) return
        inod = inod_rj_center
        d_rj(inod,ipol_scalar) = d_rj(inod,ipol_scalar)                 &
     &         + dt * (coef_exp*d_rj(inod,ipol_diffuse)                 &
     &             - d_rj(inod,ipol_advect)                             &
     &              + coef_src * d_rj(inod,ipol_source) )
!
      end subroutine scalar_diff_adv_src_euler
!
! ----------------------------------------------------------------------
!
      subroutine set_ini_adams_scalar_w_src                             &
     &         (kr_st, kr_ed, jmax, inod_rj_center,                     &
     &          ipol_advect, ipol_source, ipol_pre, coef_src,           &
     &          n_point, ntot_phys_rj, d_rj)
!
      integer(kind = kint), intent(in) :: inod_rj_center
      integer(kind = kint), intent(in) :: kr_st, kr_ed, jmax
      integer(kind = kint), intent(in) :: ipol_advect, ipol_source
      integer(kind = kint), intent(in) :: ipol_pre
      integer(kind = kint), intent(in) :: n_point, ntot_phys_rj
      real(kind = kreal), intent(in) :: coef_src
!
      real (kind=kreal), intent(inout) :: d_rj(n_point,ntot_phys_rj)
!
      integer(kind = kint) :: inod, ist, ied
!
!
      ist = (kr_st-1) * jmax + 1
      ied = kr_ed * jmax
!$omp do private (inod)
      do inod = ist, ied
         d_rj(inod,ipol_pre) =  -d_rj(inod,ipol_advect)                 &
     &                        + coef_src * d_rj(inod,ipol_source)
      end do
!$omp end do
!
      if(inod_rj_center .eq. 0) return
      d_rj(inod_rj_center,ipol_pre) = -d_rj(inod_rj_center,ipol_advect) &
     &                    + coef_src * d_rj(inod_rj_center,ipol_source)
!
      end subroutine set_ini_adams_scalar_w_src
!
! ----------------------------------------------------------------------
!
      subroutine scalar_stable_diff_src                                 &
     &         (kr_st, kr_ed, jmax, inod_rj_center,                     &
     &          ipol_source, ipol_scalar, coef_src,                     &
     &          n_point, ntot_phys_rj, d_rj)
!
      integer(kind = kint), intent(in) :: inod_rj_center
      integer(kind = kint), intent(in) :: kr_st, kr_ed, jmax
      integer(kind = kint), intent(in) :: ipol_source, ipol_scalar
      integer(kind = kint), intent(in) :: n_point, ntot_phys_rj
      real(kind = kreal), intent(in) :: coef_src
!
      real (kind=kreal), intent(inout) :: d_rj(n_point,ntot_phys_rj)
!
      integer(kind = kint) :: inod, ist, ied
!
!
      ist = (kr_st-1) * jmax + 1
      ied = kr_ed * jmax
!$omp do private (inod)
      do inod = ist, ied
        d_rj(inod,ipol_scalar) = coef_src * d_rj(inod,ipol_source)
      end do
!$omp end do
!
      if(inod_rj_center .eq. 0) return
      d_rj(inod_rj_center,ipol_scalar)                                  &
     &     = coef_src * d_rj(inod_rj_center,ipol_source)
!
      end subroutine scalar_stable_diff_src
!
! ----------------------------------------------------------------------
!
      subroutine scalar_stable_diffusion                                &
     &         (kr_st, kr_ed, jmax, inod_rj_center,                     &
     &          ipol_scalar, n_point, ntot_phys_rj, d_rj)
!
      integer(kind = kint), intent(in) :: inod_rj_center
      integer(kind = kint), intent(in) :: kr_st, kr_ed, jmax
      integer(kind = kint), intent(in) :: ipol_scalar
      integer(kind = kint), intent(in) :: n_point, ntot_phys_rj
!
      real (kind=kreal), intent(inout) :: d_rj(n_point,ntot_phys_rj)
!
      integer(kind = kint) :: inod, ist, ied
!
!
      ist = (kr_st-1) * jmax + 1
      ied = kr_ed * jmax
!$omp do private (inod)
      do inod = ist, ied
        d_rj(inod,ipol_scalar) = zero
      end do
!$omp end do
!
      if(inod_rj_center .gt. 0) d_rj(inod_rj_center,ipol_scalar) = zero
!
      end subroutine scalar_stable_diffusion
!
! ----------------------------------------------------------------------
!
      end module cal_diff_adv_src_explicit
