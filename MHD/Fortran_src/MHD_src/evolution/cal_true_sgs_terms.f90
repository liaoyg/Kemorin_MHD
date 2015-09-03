!
!      module cal_true_sgs_terms
!
!      Written by H. Matsui on Oct., 2005
!
!      subroutine cal_true_sgs_terms_pre
!      subroutine cal_true_sgs_terms_post
!
      module cal_true_sgs_terms
!
      use m_precision
!
      use m_machine_parameter
      use m_control_parameter
!
      use cal_fluxes
      use copy_nodal_fields
!
      use cal_filtering_scalars
      use cal_filtering_vectors
      use subtract_nodal_fields
!
      implicit none
!
      private :: cal_div_sgs_h_flux_true_pre
      private :: cal_div_sgs_m_flux_true_pre
      private :: cal_div_sgs_maxwell_true_pre
      private :: cal_div_sgs_induct_true_pre
      private :: cal_div_sgs_h_flux_true_post
      private :: cal_div_sgs_m_flux_true_post
      private :: cal_div_sgs_maxwell_true_post
      private :: cal_div_sgs_induct_true_post
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine cal_true_sgs_terms_pre
!
      use calypso_mpi
      use m_node_phys_data
!
      integer(kind = kint) :: i
!
       do i = 1, num_nod_phys
         if ( phys_nod_name(i).eq.'SGS_div_h_flux_true') then
           if(iflag_debug.gt.0) write(*,*)                              &
     &                         'lead  ', trim(phys_nod_name(i) )
           call cal_div_sgs_h_flux_true_pre
         else if ( phys_nod_name(i).eq.'SGS_div_m_flux_true') then
           if(iflag_debug.gt.0) write(*,*)                              &
     &                         'lead  ', trim(phys_nod_name(i) )
           call cal_div_sgs_m_flux_true_pre
         else if ( phys_nod_name(i).eq.'SGS_Lorentz_true') then
           if(iflag_debug.gt.0) write(*,*)                              &
     &                         'lead  ', trim(phys_nod_name(i) )
           call cal_div_sgs_maxwell_true_pre
         else if ( phys_nod_name(i).eq.'SGS_mag_induct_true') then
           if(iflag_debug.gt.0) write(*,*)                              &
     &                         'lead  ', trim(phys_nod_name(i) )
           call cal_div_sgs_induct_true_pre
         end if
       end do
!
      end subroutine cal_true_sgs_terms_pre
!
!-----------------------------------------------------------------------
!
      subroutine cal_true_sgs_terms_post
!
      use m_node_phys_data
!
      integer(kind = kint) :: i
!
       do i = 1, num_nod_phys
         if ( phys_nod_name(i).eq.'SGS_div_h_flux_true') then
           if(iflag_debug.gt.0) write(*,*)                              &
     &                         'lead  ', trim(phys_nod_name(i) )
           call cal_div_sgs_h_flux_true_post
         else if ( phys_nod_name(i).eq.'SGS_div_m_flux_true') then
           if(iflag_debug.gt.0) write(*,*)                              &
     &                         'lead  ', trim(phys_nod_name(i) )
           call cal_div_sgs_m_flux_true_post
         else if ( phys_nod_name(i).eq.'SGS_Lorentz_true') then
           if(iflag_debug.gt.0) write(*,*)                              &
     &                         'lead  ', trim(phys_nod_name(i) )
           call cal_div_sgs_maxwell_true_post
         else if ( phys_nod_name(i).eq.'SGS_mag_induct_true') then
           if(iflag_debug.gt.0) write(*,*)                              &
     &                         'lead  ', trim(phys_nod_name(i) )
           call cal_div_sgs_induct_true_post
         end if
       end do
!
      end subroutine cal_true_sgs_terms_post
!
!-----------------------------------------------------------------------
!
      subroutine cal_div_sgs_h_flux_true_pre
!
      use m_node_phys_address
      use cal_terms_for_heat
!
!
      call cal_flux_vector(iphys%i_h_flux, iphys%i_filter_velo,         &
     &    iphys%i_filter_temp)
      call cal_terms_4_heat(iphys%i_h_flux_div)
      call copy_scalar_component                                        &
     &   (iphys%i_h_flux_div, iphys%i_SGS_div_hf_true)
!
      end subroutine cal_div_sgs_h_flux_true_pre
!
!-----------------------------------------------------------------------
!
      subroutine cal_div_sgs_m_flux_true_pre
!
      use m_node_phys_address
      use cal_momentum_terms
!
!
      call cal_flux_tensor(iphys%i_m_flux, iphys%i_filter_velo,         &
     &    iphys%i_filter_velo)
      call cal_terms_4_momentum(iphys%i_m_flux_div)
      call copy_vector_component                                        &
     &   (iphys%i_m_flux_div, iphys%i_SGS_div_mf_true)
!
      end subroutine cal_div_sgs_m_flux_true_pre
!
!-----------------------------------------------------------------------
!
      subroutine cal_div_sgs_maxwell_true_pre
!
      use m_node_phys_address
      use cal_momentum_terms
!
!
      call cal_maxwell_tensor(iphys%i_maxwell, iphys%i_filter_magne)
      call cal_terms_4_momentum(iphys%i_maxwell_div)
      call copy_vector_component                                        &
     &   (iphys%i_maxwell_div, iphys%i_SGS_Lor_true)
!
      end subroutine cal_div_sgs_maxwell_true_pre
!
!-----------------------------------------------------------------------
!
      subroutine cal_div_sgs_induct_true_pre
!
      use m_node_phys_address
      use cal_magnetic_terms
!
!
      call cal_induction_tensor(iphys%i_induct_t,                       &
     &    iphys%i_filter_magne, iphys%i_filter_velo)
      call cal_terms_4_magnetic(iphys%i_induct_div)
      call copy_vector_component                                        &
     &   (iphys%i_induct_div, iphys%i_SGS_idct_true)
!
      end subroutine cal_div_sgs_induct_true_pre
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine cal_div_sgs_h_flux_true_post
!
      use m_node_phys_address
!
!
      call copy_scalar_component                                        &
     &   (iphys%i_SGS_div_hf_true, iphys%i_sgs_simi)
      call cal_filtered_scalar(iphys%i_SGS_div_hf_true,                 &
     &    iphys%i_h_flux_div)
      call subtract_2_nod_scalars                                       &
     &   (iphys%i_SGS_div_hf_true, iphys%i_sgs_simi,                    &
     &    iphys%i_SGS_div_hf_true)
!
      end subroutine cal_div_sgs_h_flux_true_post
!
!-----------------------------------------------------------------------
!
      subroutine cal_div_sgs_m_flux_true_post
!
      use m_node_phys_address
!
!
      call copy_vector_component                                        &
     &   (iphys%i_SGS_div_mf_true, iphys%i_sgs_simi)
      call cal_filtered_vector(iphys%i_SGS_div_mf_true,                 &
     &    iphys%i_m_flux_div)
      call subtract_2_nod_vectors                                       &
     &   (iphys%i_SGS_div_mf_true, iphys%i_sgs_simi,                    &
     &    iphys%i_SGS_div_mf_true)
!
      end subroutine cal_div_sgs_m_flux_true_post
!
!-----------------------------------------------------------------------
!
      subroutine cal_div_sgs_maxwell_true_post
!
      use m_node_phys_address
!
!
      call copy_vector_component                                        &
     &   (iphys%i_SGS_Lor_true, iphys%i_sgs_simi)
      call cal_filtered_vector(iphys%i_SGS_Lor_true,                    &
     &    iphys%i_maxwell_div)
      call subtract_2_nod_vectors                                       &
     &   (iphys%i_SGS_Lor_true, iphys%i_sgs_simi,                       &
     &    iphys%i_SGS_Lor_true)
!
      end subroutine cal_div_sgs_maxwell_true_post
!
!-----------------------------------------------------------------------
!
      subroutine cal_div_sgs_induct_true_post
!
      use m_node_phys_address
!
!
      call copy_vector_component                                        &
     &   (iphys%i_SGS_idct_true, iphys%i_sgs_simi)
      call cal_filtered_vector(iphys%i_SGS_idct_true,                   &
     &    iphys%i_induct_div)
      call subtract_2_nod_vectors                                       &
     &   (iphys%i_SGS_idct_true, iphys%i_sgs_simi,                      &
     &    iphys%i_SGS_idct_true)
!
      end subroutine cal_div_sgs_induct_true_post
!
!-----------------------------------------------------------------------
!
      end module cal_true_sgs_terms
