!>@file   m_time_data_IO.f90
!!@brief  module m_time_data_IO
!!
!!@author H. Matsui
!!@date Programmed in Oct., 2007
!
!>@brief  time and time step data for data IO
!!
!!@verbatim
!!      integer(kind = kint) function len_step_data_buf()
!!      function step_data_buffer(my_rank)
!!      subroutine read_step_data_buffer(textbuf, id_rank)
!!
!!      subroutine write_step_data(id_file, my_rank)
!!      subroutine read_step_data(id_file)
!!
!!      subroutine write_step_data_b(id_file, my_rank)
!!      subroutine read_step_data_b(id_file)
!!@endverbatim
!!
!!@n @param  my_rank   Process ID
!!@n @param  id_file   file ID for data IO
!
      module m_time_data_IO
!
      use m_precision
!
      implicit none
!
!>      Time step
      integer(kind = kint) :: i_time_step_IO
!>      Time                  @f$ t @f$
      real(kind = kreal) :: time_IO
!>      Length of time step   @f$ \Delta t @f$
      real(kind = kreal) :: delta_t_IO
!
!
      character(len=12), parameter :: TIME_HD1 = '!  domain ID'
      character(len=19), parameter :: TIME_HD2 = '!  time step number'
      character(len=16), parameter :: TIME_HD3 = '!  time, Delta t'
!
      integer(kind = kint), parameter :: l_hd = 12 + 19 + 16 + 3
      integer(kind = kint), parameter :: l_dt = 2*16 + 2*25 + 3
      integer(kind = kint), parameter :: n_buffer = l_hd+l_dt
!
      private :: TIME_HD1, TIME_HD2, TIME_HD3
      private :: l_hd, l_dt, n_buffer
!
! -------------------------------------------------------------------
!
      contains
!
! -------------------------------------------------------------------
!
      integer(kind = kint) function len_step_data_buf()
!
      len_step_data_buf = n_buffer
!
      end function len_step_data_buf
!
! -------------------------------------------------------------------
!
      function step_data_buffer(my_rank)
!
      integer(kind = kint), intent(in) :: my_rank
!
      character(len=n_buffer) :: step_data_buffer
!
      character(len=16) :: buf_pe, buf_step
      character(len=2*25) :: buf_time
!
!
      write(buf_pe,'(i16)')      my_rank
      write(buf_step,'(i16)')    i_time_step_IO
      write(buf_time,'(1p2E25.15e3)') time_IO, delta_t_IO
!
      step_data_buffer =   TIME_HD1 // char(10)                         &
     &                  // buf_pe   // char(10)                         &
     &                  // TIME_HD2 // char(10)                         &
     &                  // buf_step // char(10)                         &
     &                  // TIME_HD3 // char(10)                         &
     &                  // buf_time // char(10)
!
      end function step_data_buffer
!
! -------------------------------------------------------------------
!
      subroutine read_step_data_buffer(textbuf, id_rank)
!
      character(len=n_buffer), intent(in) :: textbuf
      integer(kind = kint), intent(inout) :: id_rank
!
      character(len=kchara) :: tmpchara(6)
!
!
      read(textbuf,'(a13,a17,a20,a17,a17,a51)') tmpchara(1:6)
      read(tmpchara(2),*) id_rank
      read(tmpchara(4),*) i_time_step_IO
      read(tmpchara(6),*) time_IO, delta_t_IO
!
      end subroutine read_step_data_buffer
!
! -------------------------------------------------------------------
! -------------------------------------------------------------------
!
      subroutine write_step_data(id_file, my_rank)
!
      integer(kind = kint), intent(in) :: id_file, my_rank
!
!
      write(id_file,'(a)'   )   TIME_HD1
      write(id_file,'(i16)') my_rank
      write(id_file,'(a)'   )   TIME_HD2
      write(id_file,'(i16)') i_time_step_IO
      write(id_file,'(a)'   )   TIME_HD3
      write(id_file,'(1p20E25.15e3)') time_IO, delta_t_IO
!
      end subroutine write_step_data
!
! -------------------------------------------------------------------
!
      subroutine read_step_data(id_file)
!
      use skip_comment_f
!
      integer(kind = kint), intent(in) :: id_file
!
      character(len=255) :: character_4_read
      integer(kind = kint) :: itmp
!
!
      call skip_comment(character_4_read,id_file)
      read(character_4_read,*) itmp
      call skip_comment(character_4_read,id_file)
      read(character_4_read,*) i_time_step_IO
      call skip_comment(character_4_read,id_file)
      read(character_4_read,*,err=99, end=99)  time_IO, delta_t_IO
!
      go to 10
  99    write(*,*) 'no delta t data... continue'
        delta_t_IO = 0.0d0
  10  continue
!
      end subroutine read_step_data
!
! -------------------------------------------------------------------
! -------------------------------------------------------------------
!
      subroutine write_step_data_b(id_file, my_rank)
!
      integer(kind = kint), intent(in) :: id_file, my_rank
!
!
      write(id_file) my_rank
      write(id_file) i_time_step_IO
      write(id_file)  time_IO, delta_t_IO
!
      end subroutine write_step_data_b
!
! -------------------------------------------------------------------
!
      subroutine read_step_data_b(id_file)
!
      integer(kind = kint), intent(in) :: id_file
!
      integer(kind = kint) :: itmp
!
!
      read(id_file) itmp
      read(id_file) i_time_step_IO
      read(id_file) time_IO, delta_t_IO
!
      end subroutine read_step_data_b
!
! -------------------------------------------------------------------
!
      end module m_time_data_IO
