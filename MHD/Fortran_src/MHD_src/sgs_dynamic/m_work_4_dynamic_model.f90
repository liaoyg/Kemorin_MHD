!
!      module m_work_4_dynamic_model
!
!  Written by Kemorin
!
      module m_work_4_dynamic_model
!
      use m_precision
      use t_ele_info_4_dynamic
      use t_work_4_dynamic_model
      use t_work_layer_correlate
!
      implicit none
!
!
        type(dynamic_model_data), save :: wk_sgs1
!
        type(dynamic_model_data), save :: wk_diff1
!
        type(dynamis_least_suare_data), save :: wk_lsq1
!
        type(dynamis_correlation_data), save :: wk_cor1
!
! ----------------------------------------------------------------------
!
!      contains
!
! ----------------------------------------------------------------------
!
      end module m_work_4_dynamic_model