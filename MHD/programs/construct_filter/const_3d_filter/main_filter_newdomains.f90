!main_filter_newdomains.f90
!
!     program  filter_newdomains

!-----------------------------------------------------------------------
      program filter_newdomains
!
      use m_precision
!
      use calypso_mpi
      use analyzer_filter_newdomains
!
      implicit none
!

      call calypso_MPI_init
!
      call filter_to_newdomain_init

      call filter_to_newdomain_analyze

      call calypso_MPI_finalize
!
      write(*,*) '***** program finished *****'
      stop
!
      end program filter_newdomains
