!
!     program  product_udt
!
!      Main routine for taking products of field data
!      Programmed by H. Matsui on 2009
!
!-----------------------------------------------------------------------
!
      program product_udt
!
      use m_precision
!
      use calypso_mpi
      use analyzer_udt_product

      implicit none
!

      call calypso_MPI_init
!
      call initialize_udt_product
      call analyze_udt_product
!
      call  calypso_MPI_finalize
!
      write(*,*) '***** program finished *****'
      stop
!
      end program product_udt
