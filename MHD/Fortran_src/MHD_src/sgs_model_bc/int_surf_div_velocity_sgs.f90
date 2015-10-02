!int_surf_div_velocity_sgs.f90
!      module int_surf_div_velocity_sgs
!
!      Written by H. Matsui on Sep., 2005
!
!      subroutine int_surf_sgs_div_velo
!      subroutine int_surf_sgs_div_vect_p
!      subroutine int_surf_sgs_div_magne
!
      module int_surf_div_velocity_sgs
!
      use m_precision
!
      use m_control_parameter
!
      implicit none
!
      private :: int_surf_sgs_div_velo_ele
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine int_surf_sgs_div_velo
!
      use m_group_data
      use m_node_phys_address
      use m_SGS_address
      use m_surf_data_torque
!
      if (iflag_commute_velo .ne. id_SGS_commute_ON) return
      call int_surf_sgs_div_velo_ele(sf_grp1, intg_point_poisson,       &
     &      nmax_sf_sgs_velo, ngrp_sf_sgs_velo, id_grp_sf_sgs_velo,     &
     &      ifilter_final, iak_diff_v, iphys%i_velo)
!
      end subroutine int_surf_sgs_div_velo
!
!-----------------------------------------------------------------------
!
      subroutine int_surf_sgs_div_vect_p
!
      use m_group_data
      use m_node_phys_address
      use m_SGS_address
      use m_surf_data_vector_p
!
      if (iflag_commute_magne .ne. id_SGS_commute_ON) return
      call int_surf_sgs_div_velo_ele(sf_grp1, intg_point_poisson,       &
     &       nmax_sf_sgs_vect_p, ngrp_sf_sgs_vect_p,                    &
     &       id_grp_sf_sgs_vect_p, ifilter_final, iak_diff_b,           &
     &       iphys%i_vecp)
!
      end subroutine int_surf_sgs_div_vect_p
!
!-----------------------------------------------------------------------
!
      subroutine int_surf_sgs_div_magne
!
      use m_group_data
      use m_node_phys_address
      use m_SGS_address
      use m_surf_data_magne
!
      if (iflag_commute_magne .ne. id_SGS_commute_ON) return
      call int_surf_sgs_div_velo_ele(sf_grp1, intg_point_poisson,       &
     &       nmax_sf_sgs_magne, ngrp_sf_sgs_magne,                      &
     &       id_grp_sf_sgs_magne, ifilter_final, iak_diff_b,            &
     &       iphys%i_magne)
!
      end subroutine int_surf_sgs_div_magne
!
!-----------------------------------------------------------------------
!
      subroutine int_surf_sgs_div_velo_ele                              &
     &          (sf_grp, n_int, nmax_grp_sf, ngrp_sf,                   &
     &           id_grp_sf, i_filter, iak_diff, i_vect)
!
      use m_phys_constants
      use m_geometry_data
      use m_SGS_model_coefs
      use m_SGS_address
      use m_finite_element_matrix
      use m_int_surface_data
      use m_jacobians
      use m_jacobian_sf_grp
      use m_filter_elength
!
      use t_group_data
!
      use delta_phys_2_each_surface
      use fem_surf_skv_sgs_commute_t
      use cal_skv_to_ff_smp_1st
!
      type(surface_group_data), intent(in) :: sf_grp
      integer(kind = kint), intent(in) :: n_int, nmax_grp_sf
      integer(kind = kint), intent(in) :: ngrp_sf(3)
      integer(kind = kint), intent(in) :: id_grp_sf(nmax_grp_sf,3)
      integer(kind = kint), intent(in) :: i_vect, iak_diff, i_filter
!
       integer(kind=kint) :: k2, nd, i, igrp, i_comp, num
!
!
      if(nmax_grp_sf .eq. 0) return
      call reset_sk6(n_scalar)
!
! --------- set vector at each node in an element
!
      do nd = 1, n_vector
        i_comp = i_vect + nd - 1
        do i = 1, ngrp_sf(nd)
          igrp = id_grp_sf(i,nd)
          num = sf_grp%istack_grp(igrp) - sf_grp%istack_grp(igrp-1)
          if(num .gt. 0) then
!
            do k2 = 1, surf1%nnod_4_surf
              call dlt_scl_phys_2_each_surface(sf_grp, igrp, k2,        &
     &            i_comp, scalar_sf)
              call fem_sf_grp_skv_sgs_div_lin_p(ele1, surf1, sf_grp,    &
     &            jac1_sf_grp_2d_q, jac1_sf_grp_2d_l, FEM1_elen,        &
     &            igrp, k2, nd, n_int, i_filter, dxe_sf, scalar_sf,     &
     &            ak_diff(1,iak_diff), fem1_wk%sk6)
            end do
!
          end if
        end do
      end do
!
      call add1_skv_to_ff_v_smp_1st(ff_smp, fem1_wk%sk6)
!
      end subroutine int_surf_sgs_div_velo_ele
!
!-----------------------------------------------------------------------
!
      end module int_surf_div_velocity_sgs
