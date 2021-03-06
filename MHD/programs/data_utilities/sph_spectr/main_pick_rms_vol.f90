!main_pick_rms_vol.f90
!
!     program  pick_rms_vol
!
!-----------------------------------------------------------------------
!    main routine for GeoFEM/Tiger version       on mar. 2000 (ver 1.0)
!
      program pick_rms_vol
!
      use m_precision
!
      use calypso_mpi
      use analyzer_pick_rms_vol
!
      implicit none
!
!
      call calypso_MPI_init
!
      call initialyze_pick_rms_vol
      call analyze_pick_rms_vol

      call calypso_MPI_finalize
!
      write(*,*) '***** program finished *****'
      stop
!
      end program pick_rms_vol
