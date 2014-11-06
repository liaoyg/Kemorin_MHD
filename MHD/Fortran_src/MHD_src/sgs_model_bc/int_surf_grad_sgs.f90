!
!      module int_surf_grad_sgs
!
!      Written by H. Matsui on Sep., 2005
!
!      subroutine int_surf_gradient_sgs(n_int, ngrp_sf, id_grp_sf,      &
!     &         i_filter, iak_diff, i_scalar)
!      subroutine int_surf_grad_commute_sgs(n_int, ngrp_sf, id_grp_sf,  &
!     &         i_filter, i_scalar)
!
      module int_surf_grad_sgs
!
      use m_precision
!
      use m_constants
      use m_phys_constants
      use m_geometry_parameter
      use m_surface_group
      use m_finite_element_matrix
!
      implicit none
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine int_surf_gradient_sgs(n_int, ngrp_sf, id_grp_sf,       &
     &         i_filter, iak_diff, i_scalar)
!
      use m_int_surface_data
      use m_SGS_model_coefs
      use m_int_surface_data
!
      use delta_phys_2_each_surface
      use fem_surf_skv_sgs_commute_1
      use cal_skv_to_ff_smp_1st
!
!
      integer(kind = kint), intent(in) :: n_int, ngrp_sf
      integer(kind = kint), intent(in) :: id_grp_sf(ngrp_sf)
      integer(kind = kint), intent(in) :: i_scalar, iak_diff, i_filter
!
      integer(kind=kint) :: k2, i, igrp
!
!
!  ---------  set number of integral points
!
      if (ngrp_sf.eq.0) return
      call reset_sk6(n_vector)
!
      do i = 1, ngrp_sf
        igrp = id_grp_sf(i)
        if ((surf_istack(igrp) - surf_istack(igrp-1)) .gt. 0) then
!
          do k2=1, nnod_4_surf
            call dlt_scl_phys_2_each_surface(igrp, k2, i_scalar,        &
     &                vect_sf(1,1) )
            call fem_sf_skv_sgs_grad_p1(igrp, k2, n_int, i_filter,      &
     &               dxe_sf, scalar_sf, ak_diff(1,iak_diff), one, sk6)
          end do
!
        end if
      end do
!
      call add3_skv_to_ff_v_smp_1st(ff_nl_smp, sk6)
!
      end subroutine int_surf_gradient_sgs
!
!-----------------------------------------------------------------------
!
      subroutine int_surf_grad_commute_sgs(n_int, ngrp_sf, id_grp_sf,   &
     &         i_filter, i_scalar)
!
      use m_int_surface_data
      use delta_phys_2_each_surface
      use fem_surf_skv_sgs_commute_1
      use cal_skv_to_ff_smp_1st
!
!
       integer(kind = kint), intent(in) :: n_int, ngrp_sf
       integer(kind = kint), intent(in) :: id_grp_sf(ngrp_sf)
       integer(kind = kint), intent(in) :: i_scalar, i_filter
!
       integer(kind=kint) :: k2, i, igrp
!
!
!  ---------  set number of integral points
!
      if (ngrp_sf.eq.0) return
      call reset_sk6(n_vector)
!
      do i = 1, ngrp_sf
        igrp = id_grp_sf(i)
        if ((surf_istack(igrp) - surf_istack(igrp-1)) .gt. 0) then
!
          do k2=1, nnod_4_surf
            call dlt_scl_phys_2_each_surface(igrp, k2, i_scalar,        &
     &              scalar_sf )
            call fem_sf_skv_grad_commute_p1(igrp, k2, n_int, i_filter,  &
     &          dxe_sf, scalar_sf, sk6)
          end do
!
        end if
      end do
!
      call add3_skv_to_ff_v_smp_1st(ff_nl_smp, sk6)
!
      end subroutine int_surf_grad_commute_sgs
!
!-----------------------------------------------------------------------
!
      end module int_surf_grad_sgs
