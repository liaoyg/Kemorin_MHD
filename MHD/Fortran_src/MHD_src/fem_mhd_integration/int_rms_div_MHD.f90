!
!     module int_rms_div_MHD
!
!     Written by H. Matsui on June, 2005
!
!      subroutine int_rms_div_v_monitor(iloop, rsig)
!      subroutine int_rms_div_b_monitor(iloop, rsig)
!      subroutine int_rms_div_a_monitor(iloop, rsig)
!      subroutine int_rms_div_v
!      subroutine int_rms_div_b
!      subroutine int_rms_div_a
!      subroutine int_rms_div_filter_v
!      subroutine int_rms_div_filter_b
!      subroutine int_rms_div_filter_a
!
      module int_rms_div_MHD
!
      use m_precision
!
      use calypso_mpi
      use m_control_parameter
      use m_geometry_data
      use m_bulk_values
!
      implicit none
!
      real (kind=kreal) :: rms_div_a_sig,  rms_div_a_sig0
      real (kind=kreal) :: rms_div_b_sig,  rms_div_b_sig0
      real (kind=kreal) :: rms_div_v_sig,  rms_div_v_sig0
!
      private :: int_rms_divergence
!
! ----------------------------------------------------------------------
!
       contains
!
! ----------------------------------------------------------------------
!
      subroutine int_rms_div_v_monitor(iloop, rsig)
!
      use m_node_phys_data
      use m_geometry_data_MHD
!
      integer(kind = kint), intent(in) :: iloop
      real(kind = kreal), intent(inout) :: rsig
!
!
      call int_rms_divergence                                           &
     &   (fluid1%istack_ele_fld_smp, ir_divv, iphys%i_velo)
!
      call MPI_allREDUCE (rms_local(ir_divv) , rms_div_v_sig, 1,        &
     &    CALYPSO_REAL, MPI_SUM, CALYPSO_COMM, ierr_MPI)
!
      rms_div_v_sig = sqrt(rms_div_v_sig / fluid1%volume)
!
      if (rms_div_v_sig .ne. 0.0d0 .and. iloop .ge.0) then
        rsig = ( rms_div_v_sig0-rms_div_v_sig ) / rms_div_v_sig
      end if
      rms_div_v_sig0 = rms_div_v_sig
!
      if (my_rank.eq.0)                                                 &
     &  write(12,*) iloop, ' : RMS(div v) = ', rms_div_v_sig0
!
      end subroutine int_rms_div_v_monitor
!
! ----------------------------------------------------------------------
!
      subroutine int_rms_div_b_monitor(iloop, rsig)
!
      use m_node_phys_data
      use m_geometry_data
!
      integer(kind = kint), intent(in) :: iloop
      real(kind = kreal), intent(inout) :: rsig
!
!
      call int_rms_divergence                                           &
     &   (ele1%istack_ele_smp, ir_divb, iphys%i_magne)
!
      call MPI_allREDUCE (rms_local(ir_divb) , rms_div_b_sig, 1,        &
     &    CALYPSO_REAL, MPI_SUM, CALYPSO_COMM, ierr_MPI)
!
      rms_div_b_sig = sqrt(rms_div_b_sig / ele1%volume)
!
      if (rms_div_b_sig .ne. 0.0d0 .and. iloop .ge.0) then
        rsig = ( rms_div_b_sig0-rms_div_b_sig ) / rms_div_b_sig
      end if
      rms_div_b_sig0 = rms_div_b_sig
!
      if (my_rank.eq.0)                                                 &
     &  write(12,*) iloop, ' : RMS(div B) = ', rms_div_b_sig0
!
!
      end subroutine int_rms_div_b_monitor
!
! ----------------------------------------------------------------------
!
      subroutine int_rms_div_a_monitor(iloop, rsig)
!
      use m_node_phys_data
      use m_geometry_data
!
      integer(kind = kint), intent(in) :: iloop
      real(kind = kreal), intent(inout) :: rsig
!
!
      call int_rms_divergence                                           &
     &   (ele1%istack_ele_smp, ir_diva, iphys%i_vecp)
!
      call MPI_allREDUCE ( rms_local(ir_diva) , rms_div_a_sig, 1,       &
     &    CALYPSO_REAL, MPI_SUM, CALYPSO_COMM, ierr_MPI)
!
      rms_div_a_sig = sqrt(rms_div_a_sig / ele1%volume)
!
      if (rms_div_a_sig .ne. 0.0d0 .and. iloop .ge.0) then
        rsig = ( rms_div_a_sig0-rms_div_a_sig ) / rms_div_a_sig
      end if
      rms_div_a_sig0 = rms_div_a_sig
!
      if (my_rank.eq.0)                                                 &
     &  write(12,*) iloop, ' : RMS(div A) = ', rms_div_a_sig0
!
!
      end subroutine int_rms_div_a_monitor
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine int_rms_div_v
!
      use m_geometry_data_MHD
      use m_node_phys_data
!
      call int_rms_divergence                                           &
     &   (fluid1%istack_ele_fld_smp, ir_divv, iphys%i_velo)
!
      end subroutine int_rms_div_v
!
! ----------------------------------------------------------------------
!
      subroutine int_rms_div_b
!
      use m_geometry_data
      use m_node_phys_data
!
!
      call int_rms_divergence                                           &
     &   (ele1%istack_ele_smp, ir_divb, iphys%i_magne)
!
      end subroutine int_rms_div_b
!
! ----------------------------------------------------------------------
!
      subroutine int_rms_div_a
!
      use m_geometry_data
      use m_node_phys_data
!
!
      call int_rms_divergence                                           &
     &   (ele1%istack_ele_smp, ir_diva, iphys%i_vecp)
!
      end subroutine int_rms_div_a
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine int_rms_div_filter_v
!
      use m_geometry_data_MHD
      use m_node_phys_data
!
      call int_rms_divergence(fluid1%istack_ele_fld_smp, ir_divv_f,     &
     &    iphys%i_filter_velo)
!
      end subroutine int_rms_div_filter_v
!
! ----------------------------------------------------------------------
!
      subroutine int_rms_div_filter_b
!
      use m_geometry_data
      use m_node_phys_data
!
!
      call int_rms_divergence(ele1%istack_ele_smp, ir_divb_f,           &
     &    iphys%i_filter_magne)
!
      end subroutine int_rms_div_filter_b
!
! ----------------------------------------------------------------------
!
      subroutine int_rms_div_filter_a
!
      use m_geometry_data
      use m_node_phys_data
!
!
      call int_rms_divergence(ele1%istack_ele_smp, ir_diva_f,           &
     &    iphys%i_filter_vecp)
!
      end subroutine int_rms_div_filter_a
!
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
      subroutine int_rms_divergence(iele_fsmp_stack, i_res, i_field)
!
      use m_machine_parameter
      use m_geometry_data
      use m_node_phys_data
      use m_finite_element_matrix
      use m_int_vol_data
!
      use fem_div_4_norm
      use sum_normalized_div
      use nodal_fld_2_each_element
!
       integer(kind = kint), intent(in)    :: iele_fsmp_stack(0:np_smp)
       integer(kind = kint), intent(in)    :: i_field
       integer(kind = kint), intent(inout) :: i_res
!
       integer(kind = kint) :: num_int, k2
!
!
      num_int = 1
!$omp workshare
      fem1_wk%scalar_1(1:ele1%numele) =  0.0d0
!$omp end workshare
!
! -------- loop for shape function for phsical values
!
      do k2=1, ele1%nnod_4_ele
!
! ---------  set vector at each node in an element
       call vector_phys_2_each_element(node1, ele1, nod_fld1,           &
     &     k2, i_field, fem1_wk%vector_1)
       call fem_rms_flux_pg(iele_fsmp_stack, num_int, k2,               &
     &     fem1_wk%vector_1, fem1_wk%scalar_1)
      end do
!
! --------- caliculate total divergence of velocity
!
      call sum_norm_of_div(ele1%numele, np_smp, iele_fsmp_stack,        &
     &    ele1%interior_ele, fem1_wk%scalar_1, rms_local(i_res) )
!
      end subroutine int_rms_divergence
!
! ----------------------------------------------------------------------
!
      end module int_rms_div_MHD
