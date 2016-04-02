!
!     module count_sgs_components
!
!      Written by H. Matsui on 2004
!      Modified by H. Matsui on July, 2007
!
!     subroutine s_count_sgs_components(numnod, numele, layer_tbl)
!        type(layering_tbl), intent(in) :: layer_tbl
!
      module count_sgs_components
!
      use m_precision
      use m_machine_parameter
!
      implicit none
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine s_count_sgs_components(numnod, numele, layer_tbl)
!
      use calypso_mpi
      use m_phys_labels
      use m_control_parameter
      use m_ele_info_4_dynamical
      use m_SGS_model_coefs
      use m_SGS_address
!
      use t_layering_ele_list
!
      integer(kind = kint), intent(in) :: numnod, numele
      type(layering_tbl), intent(in) :: layer_tbl
!
      integer(kind = kint) :: i, j, id, jd
!
!
!    count coefficients for SGS terms
!
      sgs_coefs%num_field = 0
      sgs_coefs%ntot_comp = 0
      diff_coefs%num_field = 0
      num_diff_coefs = 0
      if (iflag_t_evo_4_temp .gt. id_no_evolution) then
        if (iflag_SGS_heat .ne. id_SGS_none) then
          sgs_coefs%num_field = sgs_coefs%num_field + 1
          sgs_coefs%ntot_comp = sgs_coefs%ntot_comp + 3
          if (iflag_commute_heat .eq. id_SGS_commute_ON) then
            diff_coefs%num_field = diff_coefs%num_field + 1
            num_diff_coefs = num_diff_coefs + 3
          end if
        end if
      end if
!
      if (iflag_t_evo_4_velo .gt. id_no_evolution) then
        if (iflag_SGS_inertia .ne. id_SGS_none) then
          sgs_coefs%num_field = sgs_coefs%num_field + 1
          sgs_coefs%ntot_comp = sgs_coefs%ntot_comp + 6
          if (iflag_commute_inertia .eq. id_SGS_commute_ON) then
            diff_coefs%num_field = diff_coefs%num_field + 1
            num_diff_coefs = num_diff_coefs + 9
          end if
        end if
!
        if (iflag_SGS_lorentz .ne. id_SGS_none) then
          sgs_coefs%num_field = sgs_coefs%num_field + 1
          sgs_coefs%ntot_comp = sgs_coefs%ntot_comp + 6
          if (iflag_commute_lorentz .eq. id_SGS_commute_ON) then
            diff_coefs%num_field = diff_coefs%num_field + 1
            num_diff_coefs = num_diff_coefs + 9
          end if
        end if
!
        if (iflag_SGS_gravity .ne. id_SGS_none) then
          if(iflag_4_gravity .gt. id_turn_OFF) then
            sgs_coefs%num_field = sgs_coefs%num_field + 1
            sgs_coefs%ntot_comp = sgs_coefs%ntot_comp + 6
          end if
          if(iflag_4_composit_buo .gt. id_turn_OFF) then
            sgs_coefs%num_field = sgs_coefs%num_field + 1
            sgs_coefs%ntot_comp = sgs_coefs%ntot_comp + 6
          end if
        end if
      end if
!
      if (iflag_t_evo_4_vect_p .gt. id_no_evolution) then
        if (iflag_SGS_induction .ne. id_SGS_none) then
          sgs_coefs%num_field = sgs_coefs%num_field + 1
          sgs_coefs%ntot_comp = sgs_coefs%ntot_comp + 3
        end if
      end if
      if (iflag_t_evo_4_magne .gt. id_no_evolution) then
        if (iflag_SGS_induction .ne. id_SGS_none) then
          sgs_coefs%num_field = sgs_coefs%num_field + 1
          sgs_coefs%ntot_comp = sgs_coefs%ntot_comp + 3
          if(iflag_commute_induction .eq. id_SGS_commute_ON) then
            diff_coefs%num_field = diff_coefs%num_field + 1
            num_diff_coefs = num_diff_coefs + 9
          end if
        end if
      end if
!
      if (iflag_t_evo_4_composit .gt. id_no_evolution) then
        if (iflag_SGS_comp_flux .ne. id_SGS_none) then
          sgs_coefs%num_field = sgs_coefs%num_field + 1
          sgs_coefs%ntot_comp = sgs_coefs%ntot_comp + 3
          if (iflag_commute_c_flux .eq. id_SGS_commute_ON) then
            diff_coefs%num_field = diff_coefs%num_field + 1
            num_diff_coefs = num_diff_coefs + 3
          end if
        end if
      end if
!
      if (iflag_t_evo_4_temp .gt. id_no_evolution) then
        if(iflag_SGS_model .ne. id_SGS_none                             &
     &      .and. iflag_commute_temp .eq. id_SGS_commute_ON) then
          diff_coefs%num_field = diff_coefs%num_field + 1
          num_diff_coefs = num_diff_coefs + 3
        end if
      end if
!
      if (iflag_t_evo_4_composit .gt. id_no_evolution) then
        if(iflag_SGS_model .ne. id_SGS_none                             &
     &      .and. iflag_commute_composit .eq. id_SGS_commute_ON) then
          diff_coefs%num_field = diff_coefs%num_field + 1
          num_diff_coefs = num_diff_coefs + 3
        end if
      end if
!
      if (iflag_t_evo_4_velo .gt. id_no_evolution) then
        if(iflag_SGS_model .ne. id_SGS_none                             &
     &        .and. iflag_commute_velo .eq. id_SGS_commute_ON) then
          diff_coefs%num_field = diff_coefs%num_field + 1
          num_diff_coefs = num_diff_coefs + 9
        end if
      end if
!
      if (iflag_t_evo_4_vect_p .gt. id_no_evolution) then
        if(iflag_SGS_model .ne. id_SGS_none                             &
     &       .and. iflag_commute_magne .eq. id_SGS_commute_ON) then
          diff_coefs%num_field = diff_coefs%num_field + 1
          num_diff_coefs = num_diff_coefs + 9
        end if
      end if
!
      if (iflag_t_evo_4_magne .gt. id_no_evolution) then
        if(iflag_SGS_model .gt. id_SGS_none                             &
     &      .and. iflag_commute_magne .eq. id_SGS_commute_ON) then
          diff_coefs%num_field = diff_coefs%num_field + 1
          num_diff_coefs = num_diff_coefs + 9
        end if
      end if
!
!   set index for model coefficients
!
      call allocate_sgs_coefs_layer(layer_tbl%e_grp%num_grp,            &
     &    sgs_coefs%num_field, sgs_coefs%ntot_comp)
!
      call alloc_MHD_num_coefs(sgs_coefs)
      call alloc_MHD_coefs(numele, sgs_coefs)
!
      call alloc_MHD_num_coefs(diff_coefs)
      call alloc_MHD_coefs(numele, diff_coefs)
!
      if (iflag_dynamic_SGS .ne. id_SGS_DYNAMIC_OFF                     &
     &      .or. iflag_SGS_model.eq.id_SGS_similarity)  then
        call copy_MHD_num_coefs(sgs_coefs, sgs_coefs_nod)
        call alloc_MHD_coefs(numnod, sgs_coefs_nod)
      end if
!
       i = 1
       j = 1
       id = 1
       jd = 1
       if (iflag_t_evo_4_temp .gt. id_no_evolution) then
         if (iflag_SGS_heat .ne. id_SGS_none) then
           icomp_sgs_hf = i
           iak_sgs_hf =   j
           name_ak_sgs(j) = fhd_SGS_h_flux
           sgs_coefs%num_comps(j) = 3
           i = i + sgs_coefs%num_comps(j)
           j = j + 1
           if (iflag_commute_heat .eq. id_SGS_commute_ON) then
             icomp_diff_hf = id
             iak_diff_hf = jd
             name_ak_diff(jd) = fhd_SGS_h_flux
             diff_coefs%num_comps(jd) = 3
             id = id + diff_coefs%num_comps(jd)
             jd = jd + 1
           end if
         end if
       end if
!
       if (iflag_t_evo_4_velo .gt. id_no_evolution) then
         if (iflag_SGS_inertia .ne. id_SGS_none) then
           icomp_sgs_mf = i
           iak_sgs_mf =   j
           name_ak_sgs(j) = fhd_SGS_m_flux
           sgs_coefs%num_comps(j) = 6
           i = i + sgs_coefs%num_comps(j)
           j = j + 1
           if (iflag_commute_inertia .eq. id_SGS_commute_ON) then
             icomp_diff_mf = id
             iak_diff_mf = jd
             name_ak_diff(jd) = fhd_SGS_m_flux
             diff_coefs%num_comps(jd) = 9
             id = id + diff_coefs%num_comps(jd)
             jd = jd + 1
           end if
         end if
!
         if (iflag_SGS_lorentz .ne. id_SGS_none) then
           icomp_sgs_lor = i
           iak_sgs_lor =   j
           name_ak_sgs(j) = fhd_SGS_maxwell_t
           sgs_coefs%num_comps(j) = 6
           i = i + sgs_coefs%num_comps(j)
           j = j + 1
           if (iflag_commute_lorentz .eq. id_SGS_commute_ON) then
             icomp_diff_lor = id
             iak_diff_lor = jd
             name_ak_diff(jd) = fhd_SGS_Lorentz
             diff_coefs%num_comps(jd) = 9
             id = id + diff_coefs%num_comps(jd)
             jd = jd + 1
           end if
         end if
!
        if (iflag_SGS_gravity .ne. id_SGS_none) then
          if(iflag_4_gravity .gt. 0) then
            icomp_sgs_tbuo = i
            iak_sgs_tbuo =   j
            name_ak_sgs(j) = fhd_SGS_buoyancy
            sgs_coefs%num_comps(j) = 6
            i = i + sgs_coefs%num_comps(j)
            j = j + 1
          end if
          if(iflag_4_composit_buo .gt. id_turn_OFF) then
            icomp_sgs_cbuo = i
            iak_sgs_cbuo =   j
            name_ak_sgs(j) = fhd_SGS_comp_buo
            sgs_coefs%num_comps(j) = 6
            i = i + sgs_coefs%num_comps(j)
            j = j + 1
          end if
        end if
       end if
!
       if (iflag_t_evo_4_vect_p .gt. id_no_evolution) then
         if (iflag_SGS_induction .ne. id_SGS_none) then
           icomp_sgs_uxb = i
           iak_sgs_uxb =   j
           name_ak_sgs(j) = fhd_SGS_induction
           sgs_coefs%num_comps(j) = 3
           i = i + sgs_coefs%num_comps(j)
           j = j + 1
         end if
       end if
       if (iflag_t_evo_4_magne .gt. id_no_evolution) then
         if (iflag_SGS_induction .ne. id_SGS_none) then
           icomp_sgs_uxb = i
           iak_sgs_uxb =   j
           name_ak_sgs(j) = fhd_SGS_induction
           sgs_coefs%num_comps(j) = 3
           i = i + sgs_coefs%num_comps(j)
           j = j + 1
           if (iflag_commute_induction .eq. id_SGS_commute_ON) then
             icomp_diff_uxb = id
             iak_diff_uxb = jd
             name_ak_diff(jd) = fhd_SGS_induction
             diff_coefs%num_comps(jd) = 9
             id = id + diff_coefs%num_comps(jd)
             jd = jd + 1
           end if
         end if
       end if
!
       if (iflag_t_evo_4_composit .gt. id_no_evolution) then
         if (iflag_SGS_comp_flux .ne. id_SGS_none) then
           icomp_sgs_cf = i
           iak_sgs_cf =   j
           name_ak_sgs(j) = fhd_SGS_c_flux
           sgs_coefs%num_comps(j) = 3
           i = i + sgs_coefs%num_comps(j)
           j = j + 1
           if (iflag_commute_c_flux .eq. id_SGS_commute_ON) then
             icomp_diff_cf = id
             iak_diff_cf = jd
             name_ak_diff(jd) = fhd_SGS_c_flux
             diff_coefs%num_comps(jd) = 3
             id = id + diff_coefs%num_comps(jd)
             jd = jd + 1
           end if
         end if
       end if
!
       sgs_coefs%istack_comps(0) = 0
       do i = 1, sgs_coefs%num_field
         sgs_coefs%istack_comps(i) = sgs_coefs%istack_comps(i-1)        &
     &                              + sgs_coefs%num_comps(i)
       end do
!
!
      if (iflag_t_evo_4_temp .gt. id_no_evolution) then
        if(iflag_SGS_model .ne. id_SGS_none                             &
     &      .and. iflag_commute_temp .eq. id_SGS_commute_ON) then
            icomp_diff_t = id
            iak_diff_t =   jd
            name_ak_diff(jd) = fhd_temp
            diff_coefs%num_comps(jd) = 3
            id = id + diff_coefs%num_comps(jd)
            jd = jd + 1
        end if
      end if
!
      if (iflag_t_evo_4_composit .gt. id_no_evolution) then
        if(iflag_SGS_model .ne. id_SGS_none                             &
     &        .and. iflag_commute_composit .eq. id_SGS_commute_ON) then
            icomp_diff_c = id
            iak_diff_c =   jd
            name_ak_diff(jd) = fhd_light
            diff_coefs%num_comps(jd) = 3
            id = id + diff_coefs%num_comps(jd)
            jd = jd + 1
        end if
      end if
!
      if (iflag_t_evo_4_velo .gt. id_no_evolution) then
        if(iflag_SGS_model .ne. id_SGS_none                             &
     &      .and. iflag_commute_velo .eq. id_SGS_commute_ON) then
            icomp_diff_v = id
            iak_diff_v =   jd
            name_ak_diff(jd) = fhd_velo
            diff_coefs%num_comps(jd) = 9
            id = id + diff_coefs%num_comps(jd)
            jd = jd + 1
        end if
      end if
!
      if (iflag_t_evo_4_vect_p .gt. id_no_evolution) then
        if(iflag_SGS_model .ne. id_SGS_none                             &
     &      .and. iflag_commute_magne .eq. id_SGS_commute_ON) then
            icomp_diff_b = id
            iak_diff_b =   jd
            name_ak_diff(jd) = fhd_magne
            diff_coefs%num_comps(jd) = 9
            id = id + diff_coefs%num_comps(jd)
            jd = jd + 1
        end if
      end if
!
      if (iflag_t_evo_4_magne .gt. id_no_evolution) then
        if(iflag_SGS_model .ne. id_SGS_none                             &
     &      .and. iflag_commute_magne .eq. id_SGS_commute_ON) then
            icomp_diff_b = id
            iak_diff_b =   jd
            name_ak_diff(jd) = fhd_magne
            diff_coefs%num_comps(jd) = 9
            id = id + diff_coefs%num_comps(jd)
            jd = jd + 1
        end if
      end if
!
       diff_coefs%istack_comps(0) = 0
       do i = 1, diff_coefs%num_field
         diff_coefs%istack_comps(i) = diff_coefs%istack_comps(i-1)      &
     &                               + diff_coefs%num_comps(i)
       end do
!
!
      if(iflag_debug .gt. 0) then
        write(*,*) 'num_sgs_kinds', sgs_coefs%num_field
        write(*,*) 'num_sgs_coefs', sgs_coefs%ntot_comp
!
        if(iak_sgs_hf .gt. 0) then
          write(*,*) 'iak_sgs_hf', iak_sgs_hf, icomp_sgs_hf,            &
     &              sgs_coefs%num_comps(iak_sgs_hf),                    &
     &              trim(name_ak_sgs(iak_sgs_hf))
        end if
        if(iak_sgs_mf .gt. 0) then
          write(*,*) 'iak_sgs_mf', iak_sgs_mf, icomp_sgs_mf,            &
     &              sgs_coefs%num_comps(iak_sgs_mf),                    &
     &              trim(name_ak_sgs(iak_sgs_mf))
        end if
        if(iak_sgs_lor .gt. 0) then
          write(*,*) 'iak_sgs_lor', iak_sgs_lor, icomp_sgs_lor,         &
     &              sgs_coefs%num_comps(iak_sgs_lor),                   &
     &              trim(name_ak_sgs(iak_sgs_lor))
        end if
        if(iak_sgs_tbuo .gt. 0) then
          write(*,*) 'iak_sgs_tbuo', iak_sgs_tbuo, icomp_sgs_tbuo,      &
     &              sgs_coefs%num_comps(iak_sgs_tbuo),                  &
     &              trim(name_ak_sgs(iak_sgs_tbuo))
        end if
        if(iak_sgs_cbuo .gt. 0) then
          write(*,*) 'iak_sgs_cbuo', iak_sgs_cbuo, icomp_sgs_cbuo,      &
     &               sgs_coefs%num_comps(iak_sgs_cbuo),                 &
     &               trim(name_ak_sgs(iak_sgs_cbuo))
        end if
        if(iak_sgs_uxb .gt. 0) then
          write(*,*) 'iak_sgs_uxb', iak_sgs_uxb, icomp_sgs_uxb,         &
     &               sgs_coefs%num_comps(iak_sgs_uxb),                  &
     &               trim(name_ak_sgs(iak_sgs_uxb))
        end if
        if(iak_sgs_cf .gt. 0) then
          write(*,*) 'iak_sgs_cf', iak_sgs_cf, icomp_sgs_cf,            &
     &               sgs_coefs%num_comps(iak_sgs_cf),                   &
     &               trim(name_ak_sgs(iak_sgs_cf))
        end if
        diff_coefs%ntot_comp = diff_coefs%num_field
!
        write(*,*) 'num_diff_coefs', num_diff_coefs
        write(*,*) 'diff_coefs%num_field', diff_coefs%num_field
!
        if(iak_diff_hf .gt. 0) then
          write(*,*) 'iak_diff_hf', iak_diff_hf, icomp_diff_hf,         &
     &             diff_coefs%num_comps(iak_diff_hf),                   &
     &             trim(name_ak_diff(iak_diff_hf))
        end if
        if(iak_diff_mf .gt. 0) then
          write(*,*) 'iak_diff_mf', iak_diff_mf, icomp_diff_mf,         &
     &             diff_coefs%num_comps(iak_diff_mf),                   &
     &             trim(name_ak_diff(iak_diff_mf))
        end if
        if(iak_diff_lor .gt. 0) then
          write(*,*) 'iak_diff_lor', iak_diff_lor, icomp_diff_lor,      &
     &             diff_coefs%num_comps(iak_diff_lor),                  &
     &             trim(name_ak_diff(iak_diff_lor))
        end if
        if(iak_diff_uxb .gt. 0) then
          write(*,*) 'iak_diff_uxb', iak_diff_uxb, icomp_diff_uxb,      &
     &             diff_coefs%num_comps(iak_diff_uxb),                  &
     &             trim(name_ak_diff(iak_diff_uxb))
        end if
!
        if(iak_diff_t .gt. 0) then
          write(*,*) 'iak_diff_t', iak_diff_t, icomp_diff_t,            &
     &             diff_coefs%num_comps(iak_diff_t),                    &
     &             trim(name_ak_diff(iak_diff_t))
        end if
        if(iak_diff_v .gt. 0) then
          write(*,*) 'iak_diff_v', iak_diff_v, icomp_diff_v,            &
     &             diff_coefs%num_comps(iak_diff_v),                    &
     &             trim(name_ak_diff(iak_diff_v))
        end if
        if(iak_diff_b .gt. 0) then
          write(*,*) 'iak_diff_b', iak_diff_b, icomp_diff_b,            &
     &             diff_coefs%num_comps(iak_diff_b),                    &
     &             trim(name_ak_diff(iak_diff_b))
        end if
        if(iak_diff_c .gt. 0) then
          write(*,*) 'iak_diff_c', iak_diff_c, icomp_diff_c,            &
     &             diff_coefs%num_comps(iak_diff_c),                    &
     &             trim(name_ak_diff(iak_diff_c))
        end if
      end if
!
      end subroutine s_count_sgs_components
!
! ----------------------------------------------------------------------
!
      end module count_sgs_components
      