!
!      module cal_magnetic_pre
!
!        programmed by H.Matsui and H.Okuda
!                                    on July 2000 (ver 1.1)
!        modieied by H. Matsui on Sep., 2005
!
!      subroutine cal_magnetic_field_pre
!
      module cal_magnetic_pre
!
      use m_precision
!
      use m_machine_parameter
      use calypso_mpi
      use m_control_parameter
      use m_t_int_parameter
!
      implicit none
!
      private :: cal_magne_pre_euler, cal_magne_pre_adams
      private :: cal_magne_pre_crank, cal_magne_pre_consist_crank
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine cal_magnetic_field_pre
!
      use m_nod_comm_table
      use m_node_phys_address
      use m_node_phys_data
      use m_element_phys_data
      use m_finite_element_matrix
!
      use set_velocity_boundary
      use nod_phys_send_recv
      use cal_sgs_fluxes
      use set_magne_boundary
      use int_vol_diffusion_ele
      use int_vol_magne_pre
      use int_surf_magne_pre
!
!      use check_surface_groups
!
!
!      call check_surface_param_smp('cal_magnetic_field_pre start',     &
!     &    my_rank, sf_grp1, sf_grp_nod1)
      if ( iflag_SGS_induction .ne. id_SGS_none) then
        call cal_sgs_magne_induction
      end if
!
      call reset_ff_smps(node1%max_nod_smp, f1_l, f1_nl)
!
! lead diffusion term
!
      call int_vol_magne_diffuse_ele
!
! lead induction terms
!
      if (iflag_debug .eq. 0 ) write(*,*) 'coefs_4_time_evolution'
      if (iflag_mag_supg .gt. id_turn_OFF) then
       call int_vol_magne_pre_ele_upm                                   &
     &    (fld_ele1%ntot_phys, fld_ele1%d_fld, iphys_ele)
      else
       call int_vol_magne_pre_ele                                       &
     &    (fld_ele1%ntot_phys, fld_ele1%d_fld, iphys_ele)
      end if
!
!
      call int_surf_magne_pre_ele
!
      if (iflag_t_evo_4_magne .eq. id_explicit_euler) then
       call cal_magne_pre_euler
!
      else if (iflag_t_evo_4_magne .eq. id_explicit_adams2) then
       call cal_magne_pre_adams
!
      else if (iflag_t_evo_4_magne .eq. id_Crank_nicolson) then
       call cal_magne_pre_crank
!
      else if (iflag_t_evo_4_magne .eq. id_Crank_nicolson_cmass) then 
       call cal_magne_pre_consist_crank
      end if
!
      call set_boundary_vect                                            &
     &   (nod_bc1_b, bc_b_id_apt, iphys%i_magne, nod_fld1)
!
      call vector_send_recv(iphys%i_magne, node1, nod_comm, nod_fld1)
!
      end subroutine cal_magnetic_field_pre
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!  --------  subroutine cal_magne_pre_euler  -------
!
       subroutine cal_magne_pre_euler
!
       use cal_sol_vector_explicit
       use cal_multi_pass
!
!
       call cal_t_evo_4_vector_cd(iflag_mag_supg)
       call cal_sol_magne_pre_euler
!
       end subroutine cal_magne_pre_euler
!
! ----------------------------------------------------------------------
!  --------  subroutine cal_magne_pre_adams  -------
!
       subroutine cal_magne_pre_adams
!
       use cal_sol_vector_explicit
       use cal_multi_pass
!
!
       call cal_t_evo_4_vector_cd(iflag_mag_supg)
       call cal_sol_magne_pre_adams
!
       end subroutine cal_magne_pre_adams
!
! ----------------------------------------------------------------------
!
       subroutine cal_magne_pre_crank
!
      use m_phys_constants
      use m_t_step_parameter
      use m_finite_element_matrix
      use cal_multi_pass
      use cal_sol_vector_pre_crank
      use set_magne_boundary
      use int_sk_4_fixed_boundary
      use cal_solver_MHD
      use set_velocity_boundary
!
!
       if (coef_imp_b.gt.0.0d0) then
         call int_sk_4_fixed_magne
!         if (iflag_initial_step.eq.1) coef_imp_b = 1.0d0 / coef_imp_b
       end if
!
       call cal_t_evo_4_vector_cd(iflag_mag_supg)
!
      if (iflag_debug .eq. 0 ) write(*,*) 'bc_4_magne_rhs'
      call delete_vector_ffs_on_bc(node1, nod_bc1_b, f1_l, f1_nl)
!
      call cal_sol_magne_pre_linear
!
      if (iflag_debug .eq. 0 ) write(*,*) 'time_evolution'
      call cal_sol_magne_pre_crank
!
       end subroutine cal_magne_pre_crank
!
! ----------------------------------------------------------------------
!
       subroutine cal_magne_pre_consist_crank
!
      use m_t_step_parameter
      use m_phys_constants
      use m_finite_element_matrix
      use cal_sol_vector_pre_crank
      use set_magne_boundary
      use int_sk_4_fixed_boundary
      use cal_ff_smp_to_ffs
      use int_vol_initial_MHD
      use cal_solver_MHD
      use set_velocity_boundary
!
!
      if (coef_imp_b.gt.0.0d0) then
        call int_sk_4_fixed_magne
!         if (iflag_initial_step.eq.1) coef_imp_b = 1.0d0 / coef_imp_b
      end if
!
      call int_vol_initial_magne
      call set_ff_nl_smp_2_ff(node1, rhs_tbl1, n_vector)
!
      if (iflag_debug.eq.1) write(*,*) 'bc_4_magne_rhs'
      call delete_vector_ffs_on_bc(node1, nod_bc1_b, f1_l, f1_nl)
!
      call cal_sol_vect_p_pre_consist
!
      if (iflag_debug.eq.1)                                             &
     &        write(*,*) 'time_evolution for magnetic field'
      call cal_sol_magne_pre_crank
!
       end subroutine cal_magne_pre_consist_crank
!
! ----------------------------------------------------------------------
!
      end module cal_magnetic_pre
