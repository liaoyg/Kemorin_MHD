!>@file   MPI_sph_modes_file_IO.f90
!!@brief  module MPI_sph_modes_file_IO
!!
!!@author H. Matsui
!!@date Programmed in July, 2007
!
!>@brief merged ASCII spectr data IO routines
!!
!!@verbatim
!!      subroutine mpi_read_geom_rtp_file                               &
!!     &         (file_name, nprocs_in, my_rank_IO, sph_file)
!!      subroutine mpi_read_spectr_rj_file                              &
!!     &         (file_name, nprocs_in, my_rank_IO, sph_file)
!!      subroutine mpi_read_geom_rtm_file                               &
!!     &         (file_name, nprocs_in, my_rank_IO, sph_file)
!!      subroutine mpi_read_modes_rlm_file                              &
!!     &         (file_name, nprocs_in, my_rank_IO, sph_file)
!!        type(sph_file_data_type), intent(inout) :: sph_file
!!
!!      subroutine mpi_write_geom_rtp_file                              &
!!     &         (file_name, nprocs_in, my_rank_IO, sph_file)
!!      subroutine mpi_write_spectr_rj_file                             &
!!     &         (file_name, nprocs_in, my_rank_IO, sph_file)
!!      subroutine mpi_write_geom_rtm_file                              &
!!     &         (file_name, nprocs_in, my_rank_IO, sph_file)
!!      subroutine mpi_write_modes_rlm_file                             &
!!     &         (file_name, nprocs_in, my_rank_IO, sph_file)
!!        type(sph_file_data_type), intent(inout) :: sph_file
!!@endverbatim
!!
!!@param my_rank_IO    Process ID
!!@param file_name  file name for IO (.gz is appended in this module)
!
      module MPI_sph_modes_file_IO
!
      use m_precision
      use m_machine_parameter
!
      use t_spheric_data_IO
      use MPI_sph_modes_data_IO
      use MPI_ascii_data_IO
!
      implicit none
!
      type(calypso_MPI_IO_params), save, private :: IO_param
!
!------------------------------------------------------------------
!
      contains
!
!------------------------------------------------------------------
!
      subroutine mpi_read_geom_rtp_file                                 &
     &         (file_name, nprocs_in, my_rank_IO, sph_file)
!
      character(len=kchara), intent(in) :: file_name
      integer(kind = kint), intent(in) :: nprocs_in, my_rank_IO
      type(sph_file_data_type), intent(inout) :: sph_file
!
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &     'Read merged ascii grid file: ', trim(file_name)
      call open_read_mpi_file                                           &
     &   (file_name, nprocs_in, my_rank_IO, IO_param)
!
      call mpi_read_geom_rtp_data(IO_param,                             &
     &    sph_file%comm_IO, sph_file%sph_IO, sph_file%sph_grp_IO)
!
      call close_mpi_file(IO_param)
!
      end subroutine mpi_read_geom_rtp_file
!
!------------------------------------------------------------------
!
      subroutine mpi_read_spectr_rj_file                                &
     &         (file_name, nprocs_in, my_rank_IO, sph_file)
!
      character(len=kchara), intent(in) :: file_name
      integer(kind = kint), intent(in) :: nprocs_in, my_rank_IO
      type(sph_file_data_type), intent(inout) :: sph_file
!
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &     'Read merged ascii spectr modes file: ', trim(file_name)
      call open_read_mpi_file                                           &
     &   (file_name, nprocs_in, my_rank_IO, IO_param)
!
      call mpi_read_spectr_rj_data(IO_param,                            &
     &    sph_file%comm_IO, sph_file%sph_IO, sph_file%sph_grp_IO)
!
      call close_mpi_file(IO_param)
!
      end subroutine mpi_read_spectr_rj_file
!
!------------------------------------------------------------------
!
      subroutine mpi_read_geom_rtm_file                                 &
     &         (file_name, nprocs_in, my_rank_IO, sph_file)
!
      character(len=kchara), intent(in) :: file_name
      integer(kind = kint), intent(in) :: nprocs_in, my_rank_IO
      type(sph_file_data_type), intent(inout) :: sph_file
!
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &     'Read merged ascii grid file: ', trim(file_name)
      call open_read_mpi_file                                           &
     &   (file_name, nprocs_in, my_rank_IO, IO_param)
!
      call mpi_read_geom_rtm_data(IO_param,                             &
     &    sph_file%comm_IO, sph_file%sph_IO)
!
      call close_mpi_file(IO_param)
!
      end subroutine mpi_read_geom_rtm_file
!
!------------------------------------------------------------------
!
      subroutine mpi_read_modes_rlm_file                                &
     &         (file_name, nprocs_in, my_rank_IO, sph_file)
!
      character(len=kchara), intent(in) :: file_name
      integer(kind = kint), intent(in) :: nprocs_in, my_rank_IO
      type(sph_file_data_type), intent(inout) :: sph_file
!
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &     'Read merged ascii spectr modes file: ', trim(file_name)
      call open_read_mpi_file                                           &
     &   (file_name, nprocs_in, my_rank_IO, IO_param)
!
      call mpi_read_modes_rlm_data(IO_param,                            &
     &    sph_file%comm_IO, sph_file%sph_IO)
!
      call close_mpi_file(IO_param)
!
      end subroutine mpi_read_modes_rlm_file
!
!------------------------------------------------------------------
!------------------------------------------------------------------
!
      subroutine mpi_write_geom_rtp_file                                &
     &         (file_name, nprocs_in, my_rank_IO, sph_file)
!
      character(len=kchara), intent(in) :: file_name
      integer(kind = kint), intent(in) :: nprocs_in, my_rank_IO
      type(sph_file_data_type), intent(inout) :: sph_file
!
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &     'Write merged ascii grid file: ', trim(file_name)
      call open_write_mpi_file                                          &
     &   (file_name, nprocs_in, my_rank_IO, IO_param)
!
      call mpi_write_geom_rtp_data(IO_param,                            &
     &    sph_file%comm_IO, sph_file%sph_IO, sph_file%sph_grp_IO)
!
      call close_mpi_file(IO_param)
!
      end subroutine mpi_write_geom_rtp_file
!
!------------------------------------------------------------------
!
      subroutine mpi_write_spectr_rj_file                               &
     &         (file_name, nprocs_in, my_rank_IO, sph_file)
!
      character(len=kchara), intent(in) :: file_name
      integer(kind = kint), intent(in) :: nprocs_in, my_rank_IO
      type(sph_file_data_type), intent(inout) :: sph_file
!
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &     'Write merged ascii spectr modes file: ', trim(file_name)
      call open_write_mpi_file                                          &
     &   (file_name, nprocs_in, my_rank_IO, IO_param)
!
      call mpi_write_spectr_rj_data(IO_param,                           &
     &    sph_file%comm_IO, sph_file%sph_IO, sph_file%sph_grp_IO)
!
      call close_mpi_file(IO_param)
!
      end subroutine mpi_write_spectr_rj_file
!
!------------------------------------------------------------------
!
      subroutine mpi_write_geom_rtm_file                                &
     &         (file_name, nprocs_in, my_rank_IO, sph_file)
!
      character(len=kchara), intent(in) :: file_name
      integer(kind = kint), intent(in) :: nprocs_in, my_rank_IO
      type(sph_file_data_type), intent(inout) :: sph_file
!
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &     'Write merged ascii grid file: ', trim(file_name)
      call open_write_mpi_file                                          &
     &   (file_name, nprocs_in, my_rank_IO, IO_param)
!
      call mpi_write_geom_rtm_data(IO_param,                            &
     &    sph_file%comm_IO, sph_file%sph_IO)
!
      call close_mpi_file(IO_param)
!
      end subroutine mpi_write_geom_rtm_file
!
!------------------------------------------------------------------
!
      subroutine mpi_write_modes_rlm_file                               &
     &         (file_name, nprocs_in, my_rank_IO, sph_file)
!
      character(len=kchara), intent(in) :: file_name
      integer(kind = kint), intent(in) :: nprocs_in, my_rank_IO
      type(sph_file_data_type), intent(inout) :: sph_file
!
!
      if(my_rank_IO.eq.0 .or. i_debug .gt. 0) write(*,*)                &
     &     'Write merged ascii spectr modes file: ', trim(file_name)
      call open_write_mpi_file                                          &
     &   (file_name, nprocs_in, my_rank_IO, IO_param)
!
      call mpi_write_modes_rlm_data(IO_param,                          &
!
     &   sph_file%comm_IO, sph_file%sph_IO)
      call close_mpi_file(IO_param)
!
      end subroutine mpi_write_modes_rlm_file
!
!------------------------------------------------------------------
!
      end module MPI_sph_modes_file_IO
