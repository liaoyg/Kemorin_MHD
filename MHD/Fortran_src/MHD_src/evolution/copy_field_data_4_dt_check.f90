!copy_field_data_4_dt_check.f90
!     module copy_field_data_4_dt_check
!
      module copy_field_data_4_dt_check
!
!      Written by H. Matsui on Nov., 2009
!
      use m_constants
!
!      subroutine s_copy_field_data_for_dt_check
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine s_copy_field_data_for_dt_check
!
      use m_control_parameter
      use m_geometry_data
      use m_node_phys_data
      use copy_nodal_fields
!
!
      if((iphys%i_chk_mom_2*iphys%i_chk_mom) .gt. izero) then
        call copy_vector_component(node1, nod_fld1,                     &
     &      iphys%i_chk_mom, iphys%i_chk_mom_2)
      end if
!
      if((iphys%i_chk_press_2*iphys%i_chk_press) .gt. izero) then
        call copy_scalar_component(node1, nod_fld1,                     &
     &      iphys%i_chk_press, iphys%i_chk_press_2)
      end if
!
!
      if((iphys%i_chk_uxb_2*iphys%i_chk_uxb) .gt. izero) then
        call copy_vector_component(node1, nod_fld1,                     &
     &      iphys%i_chk_uxb, iphys%i_chk_uxb_2)
      end if
!
      if((iphys%i_chk_potential_2*iphys%i_chk_potential) .gt. izero)    &
     &      then
        call copy_scalar_component(node1, nod_fld1,                     &
     &      iphys%i_chk_potential, iphys%i_chk_potential_2)
      end if
!
!
      if((iphys%i_chk_heat_2*iphys%i_chk_heat) .gt. izero) then
        call copy_scalar_component(node1, nod_fld1,                     &
     &      iphys%i_chk_heat, iphys%i_chk_heat_2)
      end if
!
      if((iphys%i_chk_composit_2*iphys%i_chk_composit) .gt. izero) then
        call copy_scalar_component(node1, nod_fld1,                     &
     &      iphys%i_chk_composit, iphys%i_chk_composit_2)
      end if
!
!
!
      if( (iphys%i_chk_mom*iphys%i_velo) .gt. izero) then
        call copy_vector_component(node1, nod_fld1,                     &
     &      iphys%i_velo, iphys%i_chk_mom)
      end if
!
      if( (iphys%i_chk_press*iphys%i_press) .gt. izero) then
        call copy_scalar_component(node1, nod_fld1,                     &
     &      iphys%i_press, iphys%i_chk_press)
      end if
!
!
      if(iflag_t_evo_4_vect_p .gt. id_no_evolution) then
        if( (iphys%i_chk_uxb*iphys%i_vecp) .gt. izero) then
          call copy_vector_component(node1, nod_fld1,                   &
     &        iphys%i_vecp, iphys%i_chk_uxb)
        end if
      else
        if( (iphys%i_chk_uxb*iphys%i_magne) .gt. izero) then
          call copy_vector_component(node1, nod_fld1,                   &
     &        iphys%i_magne, iphys%i_chk_uxb)
        end if
      end if
!
      if( (iphys%i_chk_potential*iphys%i_mag_p) .gt. izero) then
        call copy_scalar_component(node1, nod_fld1,                     &
     &      iphys%i_mag_p, iphys%i_chk_potential)
      end if
!
!
      if( (iphys%i_chk_heat*iphys%i_temp) .gt. izero) then
        call copy_scalar_component(node1, nod_fld1,                     &
     &      iphys%i_temp, iphys%i_chk_heat)
      end if
!
      if( (iphys%i_chk_composit*iphys%i_light) .gt. izero) then
        call copy_scalar_component(node1, nod_fld1,                     &
     &      iphys%i_light, iphys%i_chk_composit)
      end if
!
      end subroutine s_copy_field_data_for_dt_check
!
! ----------------------------------------------------------------------
!
      end module copy_field_data_4_dt_check
