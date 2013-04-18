!>@file   main_sph_MHD_noviz.f90
!!@brief  program kemorin_sph_MHD_noviz
!!
!!@author H. Matsui
!!@date Programmed by H. Okuda in 2000
!!@n    Modified by H. Matsui in May, 2003 (ver 2.0)
!!@n    Connect to vizs  by H. Matsui in July 2006 (ver 2.0)
!
!>@brief  Main program for MHD dynamo simulation
!!        without visualization routines
!
     program kemorin_sph_MHD_noviz
!
      use m_precision
!
      use m_parallel_var_dof
      use analyzer_noviz_sph_MHD
!
      implicit none
!
!
      call parallel_cal_init
!
      call initialize_noviz_sph_MHD
!
      call evolution_noviz_sph_MHD
!
      call  parallel_cal_fin
!
      stop
      end program kemorin_sph_MHD_noviz
