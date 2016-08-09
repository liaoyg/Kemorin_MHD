!>@file   gauss_coefs_monitor_IO.f90
!!@brief  module gauss_coefs_monitor_IO
!!
!!@author H. Matsui
!!@date Programmed in Dec., 2012
!
!>@brief  Data arrays for gauss coefficients output
!!
!!@verbatim
!!      subroutine write_gauss_coefs_4_monitor                          &
!!     &         (my_rank, i_step, time, gauss_coefs_file_head, gauss)
!!      subroutine open_gauss_coefs_read_monitor                        &
!!     &         (id_pick, gauss_coefs_file_head, gauss)
!!      subroutine read_gauss_coefs_4_monitor(id_pick, i_step, time,    &
!!     &          gauss, ierr)
!!@endverbatim
!!
!!@n @param  my_rank   Process ID
!!@n @param  i_step    time step
!!@n @param  time      time
!!@n @param  id_pick   file ID
!!@n @param  ierr      Error flag (0:success, 1:error)
!
      module gauss_coefs_monitor_IO
!
      use m_precision
      use m_constants
      use t_pickup_sph_spectr_data
!
      implicit  none
!
!>      File ID for Gauss coefficients IO
      integer(kind = kint), parameter :: id_gauss_coef = 23
!
      private :: open_gauss_coefs_4_monitor
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine open_gauss_coefs_4_monitor                             &
     &         (gauss_coefs_file_head, gauss)
!
      use set_parallel_file_name
      use write_field_labels
!
      character(len = kchara), intent(in) :: gauss_coefs_file_head
      type(picked_spectrum_data), intent(in) :: gauss
      character(len = kchara) :: gauss_coefs_file_name
!
!
      call add_dat_extension(gauss_coefs_file_head,                     &
     &    gauss_coefs_file_name)
      open(id_gauss_coef, file = gauss_coefs_file_name,                 &
     &    form='formatted', status='old', position='append', err = 99)
      return
!
   99 continue
      open(id_gauss_coef, file = gauss_coefs_file_name,                 &
     &    form='formatted', status='replace')
!
!
      write(id_gauss_coef,'(a)')    'num_spectr, reference_radius'
      write(id_gauss_coef,'(i16,1pe25.15e3)')                           &
     &     gauss%num_sph_mode, gauss%radius_gl(1)
!
      write(id_gauss_coef,'(a)',advance='NO')    't_step    time    '
!
      call write_multi_labels(id_gauss_coef, gauss%num_sph_mode,        &
     &    gauss%gauss_mode_name)
      write(id_gauss_coef,'(a)') ''
!
      end subroutine open_gauss_coefs_4_monitor
!
! -----------------------------------------------------------------------
!
      subroutine write_gauss_coefs_4_monitor                            &
     &         (my_rank, i_step, time, gauss_coefs_file_head, gauss)
!
      character(len = kchara), intent(in) :: gauss_coefs_file_head
      integer(kind = kint), intent(in) :: my_rank
      integer(kind = kint), intent(in) :: i_step
      real(kind = kreal), intent(in) :: time
      type(picked_spectrum_data), intent(in) :: gauss
!
      integer(kind = kint) :: inum
!
!
      if(gauss%num_sph_mode .eq. izero) return
      if(my_rank .gt. izero) return
!
      call open_gauss_coefs_4_monitor(gauss_coefs_file_head, gauss)
!
      write(id_gauss_coef,'(i16,1pe23.14e3)', advance='NO')             &
     &       i_step, time
      do inum = 1, gauss%num_sph_mode
        write(id_gauss_coef,'(1pe23.14e3)', advance='NO')               &
     &       gauss%d_rj_gl(1,inum)
      end do
      write(id_gauss_coef,'(a)') ''
!
      close(id_gauss_coef)
!
      end subroutine write_gauss_coefs_4_monitor
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine open_gauss_coefs_read_monitor                          &
     &         (id_pick, gauss_coefs_file_head, gauss)
!
      use m_phys_labels
      use skip_comment_f
      use set_parallel_file_name
!
      character(len = kchara), intent(in) :: gauss_coefs_file_head
      integer(kind = kint), intent(in) :: id_pick
      type(picked_spectrum_data), intent(inout) :: gauss
!
      integer(kind = kint) :: i
      character(len = kchara) :: gauss_coefs_file_name
      character(len=255) :: tmpchara
!
!
      call add_dat_extension(gauss_coefs_file_head,                     &
     &    gauss_coefs_file_name)
      open(id_pick, file = gauss_coefs_file_name)
!
      call skip_comment(tmpchara,id_pick)
      read(id_pick,*) gauss%num_sph_mode
!
      gauss%num_layer =    1
      call alloc_num_pick_layer(gauss)
      gauss%radius_gl(1) = 2.91
!
      gauss%num_field_rj = 1
      gauss%ntot_comp_rj = 1
      call alloc_pick_sph_monitor(gauss)
!
      gauss%istack_comp_rj(1) = 1
      gauss%ifield_monitor_rj(1) = 1
      gauss%spectr_name(1) = fhd_magne
!
      call alloc_gauss_coef_monitor(gauss)
      read(id_pick,*) (tmpchara,i=1,2),                                 &
     &                 gauss%gauss_mode_name(1:gauss%num_sph_mode)
!
      end subroutine open_gauss_coefs_read_monitor
!
! -----------------------------------------------------------------------
!
      subroutine read_gauss_coefs_4_monitor(id_pick, i_step, time,      &
     &          gauss, ierr)
!
      integer(kind = kint), intent(in) :: id_pick
      integer(kind = kint), intent(inout) :: i_step, ierr
!
      real(kind = kreal), intent(inout) :: time
      type(picked_spectrum_data), intent(inout) :: gauss
!
!
      ierr = 0
      read(id_pick,*,err=99,end=99) i_step, time,                       &
     &       gauss%d_rj_gl(1,1:gauss%num_sph_mode)
      write(*,*) 'i_step', i_step, time
!
      return
!
   99 continue
      ierr = 1
      return
!
      end subroutine read_gauss_coefs_4_monitor
!
! -----------------------------------------------------------------------
!
      end module gauss_coefs_monitor_IO