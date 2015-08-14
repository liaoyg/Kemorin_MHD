!cal_average_mag_potential.f90
!      module cal_average_mag_potential
!
!        programmed by H.Matsui and H.Okuda
!                                    on July 2000 (ver 1.1)
!      Modified by H. Matsui on Aug, 2007
!
!      subroutine s_cal_average_mag_potential                           &
!     &          (numele, nnod_4_ele, ie, interior_ele)
!
      module cal_average_mag_potential
!
      use m_precision
!
      implicit none
!
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine s_cal_average_mag_potential                            &
     &          (numele, nnod_4_ele, ie, interior_ele)
!
      use calypso_mpi
      use m_control_parameter
      use m_geometry_constants
      use m_machine_parameter
      use m_geometry_data_MHD
      use m_finite_element_matrix
      use m_fem_gauss_int_coefs
      use m_jacobians
      use m_int_vol_data
      use m_node_phys_address
      use m_node_phys_data
      use m_bulk_values
!
      integer(kind=kint), intent(in) :: numele, nnod_4_ele
      integer(kind=kint), intent(in) :: ie(numele,nnod_4_ele)
      integer(kind=kint), intent(in) :: interior_ele(numele)
!
      real (kind=kreal) :: bulk_e_smp(np_smp)
!
      integer (kind=kint):: inod1, inod2, n_int, ist, ied, k1
      integer (kind = kint) :: iproc, ii, ix, inod, iele, inum
!
!
      if ( numele_in_core .gt. 0 ) then
!
!  ---------  set number of integral points
!
        n_int = intg_point_t_evo
!
! ---------  initialize
!
        ave_mp_core_local = 0.0d0
        bulk_e_smp = 0.0d0
!
!$omp parallel do private(k1,iele,ist,ied,inod1,inod2) 
        do iproc = 1, np_smp
          ist = iele_in_core_smp_stack(iproc-1)+1
          ied = iele_in_core_smp_stack(iproc)
!
          do k1 = 1, num_t_linear
            do ii= 1, n_int * n_int * n_int 
              ix = int_start3(n_int) + ii
!
!cdir nodep
              do inum = ist, ied
                iele = iele_in_core(inum)
                inod = ie(iele,k1)
!
                bulk_e_smp(iproc) = bulk_e_smp(iproc)                   &
     &                      + dble(interior_ele(iele))                  &
     &                       * d_nod(inod,iphys%i_mag_p)                &
     &                       * an(k1,ix)*xjac(iele,ix)*owe3d(ix)
              end do
            end do
          end do
        end do
!$omp end parallel do
!
!cdir noconcur
        do iproc = 1, np_smp
          ave_mp_core_local = ave_mp_core_local + bulk_e_smp(iproc)
        end do
!
        call MPI_allREDUCE (ave_mp_core_local, ave_mp_core, 1,          &
     &       CALYPSO_REAL, MPI_SUM, CALYPSO_COMM, ierr_MPI)
!
        ave_mp_core = ave_mp_core * vol_i_core
        if (my_rank.eq.0) write(84,*) ' ave_mp: ', ave_mp_core
      end if
!
      end subroutine s_cal_average_mag_potential
!
! ----------------------------------------------------------------------
!
      end module cal_average_mag_potential