!merge_dynamic_coefs.f90
!      module merge_dynamic_coefs
!
!      written by H. Matsui on Aug., 2007
!
!!      subroutine merge_coefs_4_dynamic(iflag_Csim_marging,            &
!!     &         numdir, n_layer, cor, sgs_les, dnum, c_comps, c_fields)
!!      subroutine cal_Csim_buo_by_Reynolds_ratio(n_layer, irms_buo,    &
!!     &          num_fld, num_comp, i_fld, i_comp,                     &
!!     &          sgs_les, sgs_c_coef, sgs_f_coef)
!!      subroutine single_Csim_buo_by_mf_ratio(n_layer, irms_buo,       &
!!     &          num_fld, num_comp, i_fld, i_comp,                     &
!!     &          sgs_les, sgs_c_coef, sgs_f_coef)
!
      module merge_dynamic_coefs
!
      use m_precision
      use m_constants
!
      implicit none
!
      private :: sum_by_direction_4_dynamic, ave_by_direction_4_dynamic
      private :: ave_by_correlate_4_dynamic
      private :: cal_each_components_m_coefs
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine merge_coefs_4_dynamic(iflag_Csim_marging,              &
     &         numdir, n_layer, cor, sgs_les, dnum, c_comps, c_fields)
!
      use calypso_mpi
      use t_SGS_control_parameter
!
      integer (kind = kint), intent(in) :: iflag_Csim_marging
      integer (kind = kint), intent(in) :: numdir, n_layer
      real(kind = kreal), intent(in) :: cor(n_layer,numdir)
!
      real(kind = kreal), intent(inout) :: sgs_les(n_layer,18)
      real(kind = kreal), intent(inout) :: c_fields(n_layer)
      real(kind = kreal), intent(inout) :: c_comps(n_layer, numdir)
      real(kind = kreal), intent(inout) :: dnum(n_layer)
!
!
      call cal_each_components_m_coefs                                  &
     &   (numdir, n_layer, sgs_les, c_comps)
!
      if(iflag_Csim_marging .eq. id_SGS_DIR_AVERAGE) then
        call ave_by_direction_4_dynamic(numdir, n_layer,                &
     &      sgs_les, c_comps, c_fields, dnum)
      else if(iflag_Csim_marging .eq. id_SGS_DIR_CORRELATE) then
        call ave_by_correlate_4_dynamic(numdir, n_layer,                &
     &      sgs_les, cor, c_comps, c_fields, dnum)
      else
        call sum_by_direction_4_dynamic                                 &
     &     (numdir, n_layer, sgs_les, c_fields)
      end if
!
      end subroutine merge_coefs_4_dynamic
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine cal_each_components_m_coefs                            &
     &         (numdir, n_layer, sgs_les, c_comps)
!
      integer (kind = kint), intent(in) :: numdir, n_layer
      real(kind = kreal), intent(in) :: sgs_les(n_layer,18)
!
      real(kind = kreal), intent(inout) :: c_comps(n_layer, numdir)
      integer (kind = kint) :: nd, igrp
!
!
!$omp parallel do private(igrp,nd)
      do nd = 1, numdir
        do igrp = 1, n_layer
          if (sgs_les(igrp,nd+9) .eq. 0.0d0) then
            c_comps(igrp,nd) = 0.0d0
          else
            c_comps(igrp,nd) = sgs_les(igrp,nd) / sgs_les(igrp,nd+9)
          end if
        end do
      end do
!$omp end parallel do
!
      end subroutine cal_each_components_m_coefs
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine sum_by_direction_4_dynamic                             &
     &         (numdir, n_layer, sgs_les, c_fields)
!
      integer (kind = kint), intent(in) :: numdir, n_layer
!
      real(kind = kreal), intent(inout) :: sgs_les(n_layer,18)
      real(kind = kreal), intent(inout) :: c_fields(n_layer)
!
      integer (kind = kint) :: nd, igrp
!
!
!$omp parallel
      do nd = 2, numdir
!$omp do
        do igrp = 1, n_layer
          sgs_les(igrp,1 ) = sgs_les(igrp,1 ) + sgs_les(igrp,nd  )
          sgs_les(igrp,10) = sgs_les(igrp,10) + sgs_les(igrp,9+nd)
        end do
!$omp end do nowait
      end do
!
!$omp do
      do igrp = 1, n_layer
        if (sgs_les(igrp,10) .eq. 0.0d0) then
          c_fields(igrp) = 0.0d0
        else
          c_fields(igrp) = sgs_les(igrp,1) / sgs_les(igrp,10)
        end if
      end do
!$omp end do
!$omp end parallel
!
      end subroutine sum_by_direction_4_dynamic
!
!  ---------------------------------------------------------------------
!
      subroutine ave_by_direction_4_dynamic(numdir, n_layer,            &
     &          sgs_les, c_comps, c_fields, dnum)
!
      integer (kind = kint), intent(in) :: numdir, n_layer
      real(kind = kreal), intent(in) :: sgs_les(n_layer,18)
      real(kind = kreal), intent(in) :: c_comps(n_layer, numdir)
!
      real(kind = kreal), intent(inout) :: c_fields(n_layer)
      real(kind = kreal), intent(inout) :: dnum(n_layer)
!
      integer (kind = kint) :: nd, igrp
!
!
      dnum(1:n_layer) = 0.0d0
      do nd = 1, numdir
        do igrp = 1, n_layer
          if(sgs_les(igrp,nd+9) .ne. zero) then
            dnum(igrp) = dnum(igrp) + one
          end if
        end do
      end do
      do igrp = 1, n_layer
        if (dnum(igrp) .eq. 0.0d0) dnum(igrp) = 1.0d0
      end do
!
!
      if (numdir .eq. 1) then
!
!$omp parallel do
        do igrp = 1, n_layer
          c_fields(igrp) = c_comps(igrp,1)
        end do
!$omp end parallel do
!
      else if (numdir .eq. 3) then
!
!$omp parallel do
        do igrp = 1, n_layer
          c_fields(igrp) = ( c_comps(igrp,1) + c_comps(igrp,2)          &
     &                   + c_comps(igrp,3) ) / dnum(igrp)
        end do
!$omp end parallel do
!
      else if (numdir .eq. 6) then
!
!$omp parallel do
        do igrp = 1, n_layer
          c_fields(igrp) = (  c_comps(igrp,1) +     c_comps(igrp,2)     &
     &                  +     c_comps(igrp,3) +     c_comps(igrp,4)     &
     &                  +     c_comps(igrp,5) +     c_comps(igrp,6) )   &
     &                  / dnum(igrp)
        end do
!$omp end parallel do
!
      else if (numdir .eq. 9) then
!
!$omp parallel do
        do igrp = 1, n_layer
          c_fields(igrp) = ( c_comps(igrp,1) + c_comps(igrp,2)          &
     &                    +  c_comps(igrp,3) + c_comps(igrp,4)          &
     &                    +  c_comps(igrp,5) + c_comps(igrp,6)          &
     &                    +  c_comps(igrp,7) + c_comps(igrp,8)          &
     &                    +  c_comps(igrp,9) ) / dnum(igrp)
        end do
!$omp end parallel do
!
      end if
!
      end subroutine ave_by_direction_4_dynamic
!
!  ---------------------------------------------------------------------
!
      subroutine ave_by_correlate_4_dynamic(numdir, n_layer,            &
     &          sgs_les, cor, c_comps, c_fields, dnum)
!
      integer (kind = kint), intent(in) :: numdir, n_layer
      real(kind = kreal), intent(in) :: cor(n_layer,numdir)
      real(kind = kreal), intent(in) :: sgs_les(n_layer,18)
!
      real(kind = kreal), intent(in) :: c_comps(n_layer, numdir)
      real(kind = kreal), intent(inout) :: c_fields(n_layer)
      real(kind = kreal), intent(inout) :: dnum(n_layer)
!
      integer (kind = kint) :: nd, igrp
!
!
      dnum(1:n_layer) = 0.0d0
      do nd = 1, numdir
        do igrp = 1, n_layer
          if(sgs_les(igrp,nd+9) .ne. 0.0d0) then
            dnum(igrp) = dnum(igrp) + cor(igrp,nd)
          end if
        end do
      end do
      do igrp = 1, n_layer
        if (dnum(igrp) .eq. 0.0d0) dnum(igrp) = 1.0d0
      end do
!
!
      if (numdir .eq. 1) then
!$omp parallel do
        do igrp = 1, n_layer
          c_fields(igrp) = c_comps(igrp,1)
        end do
!$omp end parallel do
      else if (numdir .eq. 3) then
!$omp parallel do
        do igrp = 1, n_layer
          c_fields(igrp) = ( c_comps(igrp,1) * cor(igrp,1)              &
     &                     + c_comps(igrp,2) * cor(igrp,2)              &
     &                     + c_comps(igrp,3) * cor(igrp,3) )            &
     &                      / dnum(igrp)
        end do
!$omp end parallel do
      else if (numdir .eq. 6) then
!$omp parallel do
        do igrp = 1, n_layer
          c_fields(igrp) = (    c_comps(igrp,1) * cor(igrp,1)           &
     &                    + two*c_comps(igrp,2) * cor(igrp,2)           &
     &                    + two*c_comps(igrp,3) * cor(igrp,3)           &
     &                    +     c_comps(igrp,4) * cor(igrp,4)           &
     &                    + two*c_comps(igrp,5) * cor(igrp,5)           &
     &                    +     c_comps(igrp,6) * cor(igrp,6) )         &
     &                     / dnum(igrp)
        end do
!$omp end parallel do
      else if (numdir .eq. 9) then
!$omp parallel do
        do igrp = 1, n_layer
          c_fields(igrp) = ( c_comps(igrp,1) * cor(igrp,1)              &
     &                  +  c_comps(igrp,2) * cor(igrp,2)                &
     &                  +  c_comps(igrp,3) * cor(igrp,3)                &
     &                  +  c_comps(igrp,4) * cor(igrp,4)                &
     &                  +  c_comps(igrp,5) * cor(igrp,5)                &
     &                  +  c_comps(igrp,6) * cor(igrp,6)                &
     &                  +  c_comps(igrp,7) * cor(igrp,7)                &
     &                  +  c_comps(igrp,8) * cor(igrp,8)                &
     &                  +  c_comps(igrp,9) * cor(igrp,9) ) / dnum(igrp)
        end do
!$omp end parallel do
      end if
!
      end subroutine ave_by_correlate_4_dynamic
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine cal_Csim_buo_by_Reynolds_ratio(n_layer, irms_buo,      &
     &          num_fld, num_comp, i_fld, i_comp,                       &
     &          sgs_les, sgs_c_coef, sgs_f_coef)
!
      integer(kind = kint), intent(in) :: irms_buo
      integer (kind = kint), intent(in) :: n_layer
      integer (kind = kint), intent(in) :: num_fld,  i_fld
      integer (kind = kint), intent(in) :: num_comp, i_comp
      real(kind = kreal), intent(in) :: sgs_les(n_layer,18)
!
      real(kind = kreal), intent(inout) :: sgs_c_coef(n_layer,num_comp)
      real(kind = kreal), intent(inout) :: sgs_f_coef(n_layer,num_fld)
!
      integer(kind = kint) :: igrp
!
!
!$omp parallel do
      do igrp = 1, n_layer
        sgs_f_coef(igrp,i_fld) = sgs_les(igrp,irms_buo) / sgs_les(igrp,4)
        sgs_f_coef(igrp,i_fld) = sgs_f_coef(igrp,i_fld) ** (two/three)
!
        sgs_c_coef(igrp,i_comp  ) = sgs_f_coef(igrp,i_fld)
        sgs_c_coef(igrp,i_comp+1) = sgs_f_coef(igrp,i_fld)
        sgs_c_coef(igrp,i_comp+2) = sgs_f_coef(igrp,i_fld)
        sgs_c_coef(igrp,i_comp+3) = sgs_f_coef(igrp,i_fld)
        sgs_c_coef(igrp,i_comp+4) = sgs_f_coef(igrp,i_fld)
        sgs_c_coef(igrp,i_comp+5) = sgs_f_coef(igrp,i_fld)
      end do
!$omp end parallel do
!
      end subroutine cal_Csim_buo_by_Reynolds_ratio
!
! ----------------------------------------------------------------------
!
      subroutine single_Csim_buo_by_mf_ratio(n_layer, irms_buo,         &
     &          num_fld, num_comp, i_fld, i_comp,                       &
     &          sgs_les, sgs_c_coef, sgs_f_coef)
!
      integer(kind = kint), intent(in) :: irms_buo
      integer (kind = kint), intent(in) :: n_layer
      integer (kind = kint), intent(in) :: num_fld,  i_fld
      integer (kind = kint), intent(in) :: num_comp, i_comp
      real(kind = kreal), intent(in) :: sgs_les(n_layer,18)
!
      real(kind = kreal), intent(inout) :: sgs_c_coef(n_layer,num_comp)
      real(kind = kreal), intent(inout) :: sgs_f_coef(n_layer,num_fld)
!
      integer(kind = kint) :: igrp
      real(kind = kreal) :: sgs_sum, rey_sum
!
!
      sgs_sum = sgs_les(1,irms_buo)
      rey_sum = sgs_les(1,4)
      do igrp = 2, n_layer
        sgs_sum = sgs_sum + sgs_les(igrp,irms_buo)
        rey_sum = rey_sum + sgs_les(igrp,4)
      end do
      sgs_sum = sgs_sum / rey_sum
!
!$omp parallel do
      do igrp = 1, n_layer
        sgs_f_coef(igrp,i_fld) = sgs_sum
!
        sgs_c_coef(igrp,i_comp  ) = sgs_sum
        sgs_c_coef(igrp,i_comp+1) = sgs_sum
        sgs_c_coef(igrp,i_comp+2) = sgs_sum
        sgs_c_coef(igrp,i_comp+3) = sgs_sum
        sgs_c_coef(igrp,i_comp+4) = sgs_sum
        sgs_c_coef(igrp,i_comp+5) = sgs_sum
      end do
!$omp end parallel do
!
      end subroutine single_Csim_buo_by_mf_ratio
!
! ----------------------------------------------------------------------
!
      end module merge_dynamic_coefs
