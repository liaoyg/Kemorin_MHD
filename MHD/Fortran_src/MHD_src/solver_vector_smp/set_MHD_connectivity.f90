!
!      module set_MHD_connectivity
!
!      Written by H. Matsui on Jan., 2006
!      modified by H. Matsui on Aug., 2007
!
!      subroutine set_connectivity_whole
!      subroutine set_connectivity_fluid
!      subroutine set_connectivity_conduct
!      subroutine set_connectivity_insulate
!
!      subroutine set_connectivity_linear
!      subroutine set_connectivity_linear_fl
!      subroutine set_connectivity_linear_cd
!      subroutine set_connectivity_linear_ins
!
      module set_MHD_connectivity
!
      use m_precision
      use m_constants
!
      use m_machine_parameter
      use m_geometry_constants
      use m_geometry_parameter
      use m_crs_connect
!
      implicit none
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine set_connectivity_whole
!
      use t_crs_connect
      use t_solver_djds
!
      use m_element_id_4_node
      use m_next_node_id_4_node
      use m_solver_djds_MHD
!
      use set_crs_connect_type
      use set_geometry_to_types
      use reordering_djds_smp_type
      use set_djds_smp_ordering_type
!
!C +-------------------------------+
!C | set connectivity in CRS array |
!C +-------------------------------+
!C===
      call s_set_crs_connect_type(np_smp, numnod, inod_smp_stack,       &
     &          ntot_next_nod_4_node, inod_next_stack_4_node,           &
     &          inod_next_4_node, MHD_CRS)
!
!C +-----------------+
!C | DJDS reordering |
!C +-----------------+
!C===
!C
      call s_reordering_djds_smp_type(np_smp, numnod, internal_node,    &
     &    inter_smp_stack, MHD_CRS, DJDS_entire)
!C
!      write(*,*) 'STACKmc', size(DJDS_entire%STACKmc)
!      write(*,*) 'NLmaxHYP', size(DJDS_entire%NLmaxHYP),               &
!     &          DJDS_entire%NHYP
!      write(*,*) 'NUmaxHYP', size(DJDS_entire%NUmaxHYP),               &
!     &          DJDS_entire%NHYP
!      write(*,*) 'OLDtoNEW', size(DJDS_entire%OLDtoNEW),               &
!     &          DJDS_entire%NP
!      write(*,*) 'OLDtoNEW_DJDS_L', size(DJDS_entire%OLDtoNEW_DJDS_L)
!      write(*,*) 'OLDtoNEW_DJDS_U', size(DJDS_entire%OLDtoNEW_DJDS_U)
!      write(*,*) 'indexDJDS_L', size(DJDS_entire%indexDJDS_L),         &
!     &          DJDS_entire%PEsmpTOT, DJDS_entire%NLmax,NHYP
!      write(*,*) 'indexDJDS_U', size(DJDS_entire%indexDJDS_U),         &
!     &          DJDS_entire%PEsmpTOT, DJDS_entire%NUmax,NHYP
!      write(*,*) 'itemDJDS_L', size(DJDS_entire%itemDJDS_L),           &
!     &          DJDS_entire%itotal_l
!      write(*,*) 'itemDJDS_U', size(DJDS_entire%itemDJDS_U),           &
!     &          DJDS_entire%itotal_u
!      write(*,*) 'PEon', size(DJDS_entire%PEon)
!      write(*,*) 'COLORon', size(DJDS_entire%COLORon)
!
!C +--------------------------------------+
!C | set new communication table 4 solver |
!C +--------------------------------------+
!C===
!C
      call set_nod_comm_tbl_2_type(DJDS_comm_etr)
      call set_new_comm_table_type(numnod, DJDS_comm_etr, DJDS_entire)
!
      call dealloc_type_crs_connect(MHD_CRS)
!
      end subroutine set_connectivity_whole
!
!-----------------------------------------------------------------------
!
      subroutine set_connectivity_fluid
!
      use m_geometry_data_MHD
      use m_solver_djds_MHD
      use reordering_djds_smp_type
      use set_djds_smp_ordering_type
      use m_element_id_4_node
      use m_next_node_id_4_node
      use set_element_id_4_node
      use set_crs_connect_type
!
!
      call set_layerd_ele_id_4_node(nnod_4_ele,                         &
     &    iele_fl_start, iele_fl_end)
      call const_next_nod_id_4_node
!
      call s_set_crs_connect_type(np_smp, numnod, inod_smp_stack,       &
     &    ntot_next_nod_4_node, inod_next_stack_4_node,                 &
     &    inod_next_4_node, MHD_CRS)
!
      call s_reordering_djds_smp_type(np_smp, numnod, internal_node,    &
     &    inter_smp_stack, MHD_CRS, DJDS_fluid)
      call set_new_comm_table_type(numnod, DJDS_comm_fl, DJDS_fluid)
!
      call dealloc_type_crs_connect(MHD_CRS)
      call deallocate_iele_belonged
      call deallocate_inod_next_node
!
      end subroutine set_connectivity_fluid
!
!-----------------------------------------------------------------------
!
      subroutine set_connectivity_conduct
!
      use m_geometry_data_MHD
      use m_solver_djds_MHD
      use m_element_id_4_node
      use m_next_node_id_4_node
      use set_element_id_4_node
      use reordering_djds_smp_type
      use set_djds_smp_ordering_type
      use set_crs_connect_type
!
!
      call set_layerd_ele_id_4_node(nnod_4_ele,                         &
     &    iele_cd_start, iele_cd_end)
      call const_next_nod_id_4_node
!
      call s_set_crs_connect_type(np_smp, numnod, inod_smp_stack,       &
     &    ntot_next_nod_4_node, inod_next_stack_4_node,                 &
     &    inod_next_4_node, MHD_CRS)
!
      call s_reordering_djds_smp_type(np_smp, numnod, internal_node,    &
     &    inter_smp_stack, MHD_CRS, DJDS_conduct)
      call set_new_comm_table_type(numnod, DJDS_comm_etr, DJDS_conduct)
!
      call dealloc_type_crs_connect(MHD_CRS)
      call deallocate_iele_belonged
      call deallocate_inod_next_node
!
      end subroutine set_connectivity_conduct
!
!-----------------------------------------------------------------------
!
      subroutine set_connectivity_insulate
!
      use m_geometry_data_MHD
      use m_solver_djds_MHD
      use m_element_id_4_node
      use m_next_node_id_4_node
      use set_element_id_4_node
      use reordering_djds_smp_type
      use set_djds_smp_ordering_type
      use set_crs_connect_type
!
!
      call set_layerd_ele_id_4_node(nnod_4_ele,                         &
     &    iele_ins_start, iele_ins_end)
      call const_next_nod_id_4_node
!
      call s_set_crs_connect_type(np_smp, numnod, inod_smp_stack,       &
     &    ntot_next_nod_4_node, inod_next_stack_4_node,                 &
     &    inod_next_4_node, MHD_CRS)
!
!
      call s_reordering_djds_smp_type(np_smp, numnod, internal_node,    &
     &    inter_smp_stack, MHD_CRS, DJDS_insulator)
      call set_new_comm_table_type(numnod,                              &
     &    DJDS_comm_etr, DJDS_insulator)
!
      call dealloc_type_crs_connect(MHD_CRS)
      call deallocate_iele_belonged
      call deallocate_inod_next_node
!
      end subroutine set_connectivity_insulate
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine set_connectivity_linear
!
      use t_crs_connect
      use t_solver_djds
!
      use m_nod_comm_table
      use m_element_id_4_node
      use m_next_node_id_4_node
      use m_solver_djds_MHD
      use set_element_id_4_node
      use set_crs_connect_type
      use reordering_djds_smp_type
      use set_djds_smp_ordering_type
!
!
      if ( nnod_4_ele .ne. num_t_linear) then
        call set_layerd_ele_id_4_node(num_t_linear, ione, numele)
        call const_next_nod_id_4_node
!
        call s_set_crs_connect_type(np_smp, numnod, inod_smp_stack,     &
     &      ntot_next_nod_4_node, inod_next_stack_4_node,               &
     &      inod_next_4_node, MHD_CRS)
!
        call s_reordering_djds_smp_type(np_smp, numnod, internal_node,  &
     &      inter_smp_stack, MHD_CRS, DJDS_linear)
        call set_new_comm_table_type(numnod,                            &
     &      DJDS_comm_etr, DJDS_linear)
!
        call dealloc_type_crs_connect(MHD_CRS)
        call deallocate_iele_belonged
        call deallocate_inod_next_node
      else
        call link_djds_connect_structs(DJDS_entire, DJDS_linear)
      end if
!
      end subroutine set_connectivity_linear
!
!-----------------------------------------------------------------------
!
      subroutine set_connectivity_linear_fl
!
      use m_geometry_data_MHD
      use m_solver_djds_MHD
      use m_element_id_4_node
      use m_next_node_id_4_node
      use set_element_id_4_node
      use set_crs_connect_type
      use reordering_djds_smp_type
      use set_djds_smp_ordering_type
!
!
      if ( nnod_4_ele .ne. num_t_linear) then
        call set_layerd_ele_id_4_node(num_t_linear,                     &
     &      iele_fl_start, iele_fl_end)
        call const_next_nod_id_4_node
!
        call s_set_crs_connect_type(np_smp, numnod, inod_smp_stack,     &
     &      ntot_next_nod_4_node, inod_next_stack_4_node,               &
     &      inod_next_4_node, MHD_CRS)
!
        call s_reordering_djds_smp_type(np_smp, numnod, internal_node,  &
     &      inter_smp_stack, MHD_CRS, DJDS_fl_l)
        call set_new_comm_table_type(numnod, DJDS_comm_fl, DJDS_fl_l)
!
        call dealloc_type_crs_connect(MHD_CRS)
        call deallocate_iele_belonged
        call deallocate_inod_next_node
      else
        call link_djds_connect_structs(DJDS_fluid, DJDS_fl_l)
      end if
!
!
      end subroutine set_connectivity_linear_fl
!
!-----------------------------------------------------------------------
!
      subroutine set_connectivity_linear_cd

      use m_geometry_data_MHD
      use m_solver_djds_MHD
      use m_element_id_4_node
      use m_next_node_id_4_node
      use set_element_id_4_node
!
      use reordering_djds_smp_type
      use set_djds_smp_ordering_type
      use set_crs_connect_type


      if ( nnod_4_ele .ne. num_t_linear) then
        call set_layerd_ele_id_4_node(num_t_linear,                     &
     &      iele_cd_start, iele_cd_end)
        call const_next_nod_id_4_node
!
        call s_set_crs_connect_type(np_smp, numnod, inod_smp_stack,     &
     &      ntot_next_nod_4_node, inod_next_stack_4_node,               &
     &      inod_next_4_node, MHD_CRS)
!
        call s_reordering_djds_smp_type(np_smp, numnod, internal_node,  &
     &      inter_smp_stack, MHD_CRS, DJDS_cd_l)
        call set_new_comm_table_type(numnod,                            &
     &      DJDS_comm_etr, DJDS_cd_l)
!
        call dealloc_type_crs_connect(MHD_CRS)
        call deallocate_iele_belonged
        call deallocate_inod_next_node
      else
        call link_djds_connect_structs(DJDS_conduct, DJDS_cd_l)
      end if

      end subroutine set_connectivity_linear_cd
!
!-----------------------------------------------------------------------
!
      subroutine set_connectivity_linear_ins

      use m_geometry_data_MHD
      use m_solver_djds_MHD
      use m_element_id_4_node
      use m_next_node_id_4_node
      use set_element_id_4_node
!
      use reordering_djds_smp_type
      use set_djds_smp_ordering_type
      use set_crs_connect_type


      if ( nnod_4_ele .ne. num_t_linear) then
        call set_layerd_ele_id_4_node(num_t_linear,                     &
     &      iele_ins_start, iele_ins_end)
        call const_next_nod_id_4_node
!
        call s_set_crs_connect_type(np_smp, numnod, inod_smp_stack,     &
     &      ntot_next_nod_4_node, inod_next_stack_4_node,               &
     &      inod_next_4_node, MHD_CRS)
!
        call s_reordering_djds_smp_type(np_smp, numnod, internal_node,  &
     &      inter_smp_stack, MHD_CRS, DJDS_ins_l)
        call set_new_comm_table_type(numnod,                            &
     &      DJDS_comm_etr, DJDS_ins_l)
!
        call dealloc_type_crs_connect(MHD_CRS)
        call deallocate_iele_belonged
        call deallocate_inod_next_node
      else
        call link_djds_connect_structs(DJDS_insulator, DJDS_ins_l)
      end if
!
      end subroutine set_connectivity_linear_ins
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      end module set_MHD_connectivity
