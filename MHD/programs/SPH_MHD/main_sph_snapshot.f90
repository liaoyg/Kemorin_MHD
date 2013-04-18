!>@file   main_sph_snapshot.f90
!!@brief  program kemorin_sph_snapshot
!!
!!@author H. Matsui
!!@date Programmed by H. Okuda in 2000
!!@n    Modified by H. Matsui in May, 2003 (ver 2.0)
!!@n    Connect to vizs  by H. Matsui in July 2006 (ver 2.0)
!
!>@brief  Main program to evaluate snapshots from spectr data
!
      program kemorin_sph_snapshot
!
      use m_precision
!
      use m_parallel_var_dof
      use analyzer_sph_snap
!
      implicit none
!
!
      call parallel_cal_init
!
      call initialize_sph_snap
!
      call evolution_sph_snap
!
      call parallel_cal_fin
!
      write(*,*) '***** program finished *****'
      stop
!
      end program kemorin_sph_snapshot
