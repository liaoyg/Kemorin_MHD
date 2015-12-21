!
!      module cal_rms_potentials
!
!     Written by H. Matsui on March, 2006
!
!      subroutine cal_rms_pressure_4_loop(iloop, rsig)
!      subroutine cal_rms_scalar_potential(iloop, rsig)
!
      module cal_rms_potentials
!
      use m_precision
!
      use calypso_mpi
      use m_constants
      use m_machine_parameter
      use m_geometry_data
      use m_geometry_data_MHD
      use m_bulk_values
      use int_all_energy
!
      implicit none
!
      real(kind = kreal) :: ave_pr, rms_pr, ave_pr0, rms_pr0
      real(kind = kreal) :: ave_mp, rms_mp, ave_mp0, rms_mp0
!
      private :: ave_pr, rms_pr, ave_pr0, rms_pr0
      private :: ave_mp, rms_mp, ave_mp0, rms_mp0
      private :: ione, zero
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine cal_rms_pressure_4_loop(iloop, rsig)
!
      use m_node_phys_data
!
      integer(kind = kint), intent(in) :: iloop
      real(kind = kreal), intent(inout) :: rsig
!
      integer(kind = kint) :: num_int
!
!
      num_int = ione
!
      rms_local(i_rms%i_press) = zero
      bulk_local(j_ave%i_press) = zero
      call int_all_4_scalar (iele_fl_smp_stack, num_int,                &
     &         i_rms%i_press, j_ave%i_press, iphys%i_press)
!
      call MPI_allREDUCE ( bulk_local(j_ave%i_press) , ave_pr, ione,    &
     &    CALYPSO_REAL, MPI_SUM, CALYPSO_COMM, ierr_MPI)
      call MPI_allREDUCE ( rms_local(i_rms%i_press) , rms_pr, ione,     &
     &    CALYPSO_REAL, MPI_SUM, CALYPSO_COMM, ierr_MPI)
!
      if (iloop .eq. 0) then
        ave_pr0 = ave_pr
        rms_pr0 = rms_pr
      end if
      rsig = (rms_pr - rms_pr0) / rms_pr0
!
      if (iflag_debug.eq.1)                                             &
     &         write(12,*) 'average and RMS of presssur correction: ',  &
     &         iloop, ave_pr, rms_pr
!
      end subroutine cal_rms_pressure_4_loop
!
! ----------------------------------------------------------------------
!
      subroutine cal_rms_scalar_potential(iloop, rsig)
!
      use m_node_phys_data
!
      integer(kind = kint), intent(in) :: iloop
      real(kind = kreal), intent(inout) :: rsig
!
      integer(kind = kint) :: num_int
!
!
      num_int = ione
!
      rms_local(i_rms%i_mag_p) = zero
      bulk_local(j_ave%i_mag_p) = zero
      call int_all_4_scalar (ele1%istack_ele_smp, num_int,              &
     &    i_rms%i_mag_p, j_ave%i_mag_p, iphys%i_mag_p)
!
      call MPI_allREDUCE ( bulk_local(j_ave%i_mag_p) , ave_mp, ione,    &
     &    CALYPSO_REAL, MPI_SUM, CALYPSO_COMM, ierr_MPI)
      call MPI_allREDUCE ( rms_local(i_rms%i_mag_p) , rms_mp, ione,     &
     &    CALYPSO_REAL, MPI_SUM, CALYPSO_COMM, ierr_MPI)
!
      if (iloop .eq. 0) then
        ave_mp0 = ave_mp
        rms_mp0 = rms_mp
      end if
      rsig = (rms_mp - rms_mp0) / rms_mp0
!
      if (iflag_debug.eq.1)                                             &
     &         write(12,*) 'average and RMS of potential correction: ', &
     &         iloop, ave_mp, rms_mp
!
      end subroutine cal_rms_scalar_potential
!
! ----------------------------------------------------------------------
!
      end module cal_rms_potentials
