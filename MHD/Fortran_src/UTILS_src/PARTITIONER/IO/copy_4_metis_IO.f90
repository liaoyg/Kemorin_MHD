!
!     module copy_4_metis_IO
!
!     Written by H. Matsui on Aug., 2007
!
!      subroutine copy_graph_4_metis_IO(numnod, internal_node)
!
      module copy_4_metis_IO
!
      use m_precision
!
      implicit  none
!
!   --------------------------------------------------------------------
!
      contains
!
!   --------------------------------------------------------------------
!
      subroutine copy_graph_4_metis_IO(numnod, internal_node)
!
      use m_geometry_graph
      use m_metis_IO
!
      integer(kind = kint), intent(in) :: numnod, internal_node
!
!
      nnod_metis_IO =  internal_node
      ntot_metis_IO =  ntot_graph_inter_nod
      nedge_metis_IO = ntot_graph_inter_nod / 2
!
      call allocate_metis_grp_stack_IO
      call allocate_metis_graph_IO
!
      istack_metis_IO(0:internal_node)                                  &
     &      = istack_graph_inter_nod(0:internal_node)
      igraph_metis_IO(1:ntot_graph_inter_nod)                           &
     &      = igraph_inter_nod(1:ntot_graph_inter_nod)
!
      call deallocate_internod_graph
!
      end subroutine copy_graph_4_metis_IO
!
!   --------------------------------------------------------------------
!
      end module copy_4_metis_IO
