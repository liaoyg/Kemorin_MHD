!>@file   t_ctl_data_filter_files.f90
!!@brief  module t_ctl_data_filter_files
!!
!!@author H. Matsui
!!@date Programmed in July, 2006
!
!>@brief  Structure for reading parameters for filtering files
!!
!!@verbatim
!!      subroutine read_filter_fnames_control                           &
!!     &         (hd_block, iflag, ffile_ctl)
!!
!!  ---------------------------------------------------------------------
!!
!!      begin filter_files_def
!!        filter_file_prefix           'filter/filter_node'
!!        filter_elength_prefix        'filter/filter_elength'
!!        filter_moment_prefix         'filter/filter_moms'
!!        filter_coefs_prefix          'filter/filter_coef'
!!        wider_filter_prefix          'filter/filter_coef_2'
!!
!!        filter_elen_format        'ascii'
!!        filter_3d_format          'binary'
!!        filter_wide_format        'gzip'
!!
!!        model_coef_rst_prefix       'model_coefs_ini'
!!        commutel_coef_rst_prefix    'commute_coefs_ini'
!!        model_coef_rst_format      'merged_gz'
!!        commute_coef_rst_format    'merged_gz'
!!      end filter_files_def
!!
!!  ---------------------------------------------------------------------
!!@endverbatim
!
      module t_ctl_data_filter_files
!
      use m_precision
      use t_control_elements
!
      implicit  none
!
!
!>      Structure for filtering files
      type filter_file_control
!>        Structure for filter file for nodes
        type(read_character_item) :: filter_head_ctl
!>        Structure for filter coefficients file for nodes
        type(read_character_item) :: filter_coef_head_ctl
!>        Structure for filter size file for nodes
        type(read_character_item) :: filter_elen_head_ctl
!>        Structure for filter moments file for nodes
        type(read_character_item) :: filter_moms_head_ctl
!
!>        Structure for wider filter file for nodes
        type(read_character_item) :: filter_wide_head_ctl
!
!>        Structure for model coefficients file for nodes
        type(read_character_item) :: model_coef_ini_head_ctl
!>        Structure for commutation coefficients file for nodes
        type(read_character_item) :: commute_coef_ini_head_ctl
!
!>        Structure for file format of element length
        type(read_character_item) :: filter_elen_format
!>        Structure for file format of 3D filter file
        type(read_character_item) :: filter_3d_format
!>        Structure for file format of wider filter file
        type(read_character_item) :: filter_wide_format
!
!>        Structure for file format of model coefficient
        type(read_character_item) :: model_coef_rst_format
!>        Structure for file format of commutation coefficient
        type(read_character_item) :: commute_coef_rst_format
      end type filter_file_control
!
!     flags for filter file headers
!
      character(len=kchara), parameter                                  &
     &         :: hd_filter_head_ctl =       'filter_file_prefix'
      character(len=kchara), parameter                                  &
     &         :: hd_filter_elen_head_ctl =  'filter_elength_prefix'
      character(len=kchara), parameter                                  &
     &         :: hd_filter_moms_head_ctl =  'filter_moment_prefix'
      character(len=kchara), parameter                                  &
     &         :: hd_filter_coef_head_ctl =  'filter_coefs_prefix'
      character(len=kchara), parameter                                  &
     &         :: hd_filter_wide_head =      'wider_filter_prefix'
      character(len=kchara), parameter                                  &
     &         :: hd_model_coef_ini_head =   'model_coef_rst_prefix'
      character(len=kchara), parameter                                  &
     &         :: hd_commute_coef_ini_head = 'commutel_coef_rst_prefix'
!
      character(len=kchara), parameter                                  &
     &         :: hd_filter_elen_fmt = 'filter_elen_format'
      character(len=kchara), parameter                                  &
     &         :: hd_filter_3d_fmt =   'filter_3d_format'
      character(len=kchara), parameter                                  &
     &         :: hd_filter_wide_fmt = 'filter_wide_format'
      character(len=kchara), parameter                                  &
     &         :: hd_model_coef_rst_format = 'model_coef_rst_format'
      character(len=kchara), parameter                                  &
     &         :: hd_commute_c_rst_format = 'commute_coef_rst_format'
!
      private :: hd_filter_head_ctl, hd_filter_elen_head_ctl
      private :: hd_filter_moms_head_ctl, hd_filter_coef_head_ctl
      private :: hd_filter_wide_head, hd_commute_c_rst_format
      private :: hd_model_coef_ini_head, hd_commute_coef_ini_head
      private :: hd_filter_elen_fmt, hd_filter_3d_fmt
      private :: hd_filter_wide_fmt, hd_model_coef_rst_format
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine read_filter_fnames_control                             &
     &         (hd_block, iflag, ffile_ctl)
!
      use m_machine_parameter
      use skip_comment_f
!
      character(len=kchara), intent(in) :: hd_block
!
      integer(kind = kint), intent(inout) :: iflag
      type(filter_file_control), intent(inout) :: ffile_ctl
!
!
      if(right_begin_flag(hd_block) .eq. 0) return
      if (iflag .gt. 0) return
      do
        call load_ctl_label_and_line
!
        call find_control_end_flag(hd_block, iflag)
        if(iflag .gt. 0) exit
!
!
        call read_chara_ctl_type                                        &
     &     (hd_filter_head_ctl, ffile_ctl%filter_head_ctl)
        call read_chara_ctl_type                                        &
     &     (hd_filter_coef_head_ctl, ffile_ctl%filter_coef_head_ctl)
        call read_chara_ctl_type                                        &
     &     (hd_filter_elen_head_ctl, ffile_ctl%filter_elen_head_ctl)
        call read_chara_ctl_type                                        &
     &     (hd_filter_moms_head_ctl, ffile_ctl%filter_moms_head_ctl)
        call read_chara_ctl_type                                        &
     &     (hd_filter_wide_head, ffile_ctl%filter_wide_head_ctl)
        call read_chara_ctl_type                                        &
     &     (hd_model_coef_ini_head, ffile_ctl%model_coef_ini_head_ctl)
        call read_chara_ctl_type(hd_commute_coef_ini_head,              &
     &      ffile_ctl%commute_coef_ini_head_ctl)
!
        call read_chara_ctl_type                                        &
     &     (hd_filter_elen_fmt, ffile_ctl%filter_elen_format)
        call read_chara_ctl_type                                        &
     &     (hd_filter_3d_fmt, ffile_ctl%filter_3d_format)
        call read_chara_ctl_type                                        &
     &     (hd_filter_wide_fmt, ffile_ctl%filter_wide_format)
        call read_chara_ctl_type                                        &
     &     (hd_model_coef_rst_format, ffile_ctl%model_coef_rst_format)
        call read_chara_ctl_type(hd_commute_c_rst_format,               &
     &      ffile_ctl%commute_coef_rst_format)
      end do
!
      end subroutine read_filter_fnames_control
!
!  ---------------------------------------------------------------------
!
      end module t_ctl_data_filter_files
