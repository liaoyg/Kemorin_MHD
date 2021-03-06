!
!     module count_export_peri
!
!     Written by H. Matsui
!     modified by H. Matsui on Aug., 2007
!
!      subroutine count_export_peri_linear(ipe, jpe, inod)
!      subroutine count_export_peri_quad(ipe, jpe, kpe, inod)
!
      module count_export_peri
!
      use m_precision
!
      use m_size_4_plane
      use m_size_of_cube
      use m_comm_data_cube_kemo
      use m_neighb_range_cube
      use m_sleeve_cube
      use set_comm_nod_4_cube
!
      implicit none
!
! ----------------------------------------------------------------------
!
      contains
!
! ----------------------------------------------------------------------
!
      subroutine count_export_peri_linear(ipe, jpe, inod)
!
      integer (kind = kint) :: ipe, jpe
      integer (kind = kint) :: inod
!
      integer (kind = kint) :: inp, jnp, knp
!
!
!  outdside (x>xmax)
!
!                                     .... count nodes 
!
            if (ipe .eq. ndx) then
             inp = 1
             do knp=knp_st,knp_end
              do jnp=jnp_st,jnp_end
!
               call set_boundary_size(inp, jnp, knp)
               is = nxi+1
               ie = nxi+ndepth

               neibpetot = neibpetot  + 1
               call count_node_id(inod)

               stack_export(neibpetot) = inod

              enddo
             enddo
            endif

!
!    ---   outside wall (x<xmin)
!
!                                     .... count nodes 
            if (ipe .eq. 1) then
             do knp=knp_st,knp_end
              do jnp=jnp_st,jnp_end

               call set_boundary_size(inp, jnp, knp)
               is = ndepth+1
               ie = 2*ndepth

               neibpetot = neibpetot  + 1
               call count_node_id(inod)

               stack_export(neibpetot) = inod

              enddo
             enddo
            endif
!
!  outdside (y<ymax)
!
!                                     .... count nodes 
            if ( jpe .eq. ndy ) then
             jnp = 1
             do knp=knp_st,knp_end
              do inp=inp_st,inp_end

               call set_boundary_size(inp, jnp, knp)
               js = nyi+1
               je = nyi+ndepth

               neibpetot = neibpetot  + 1
               call count_node_id(inod)

               stack_export(neibpetot) = inod

              enddo
             enddo
            endif

!  outdside (y<ymin)
!
!                                     .... count nodes 
            if ( jpe .eq. 1 ) then
             do knp=knp_st,knp_end
              do inp=inp_st,inp_end

               call set_boundary_size(inp, jnp, knp)
               js = ndepth+1
               je = 2*ndepth

               neibpetot = neibpetot  + 1
               call count_node_id(inod)

               stack_export(neibpetot) = inod

              enddo
             enddo
            endif
!
!  outdside (x>xmax, y>ymax)
!
            if ( ipe .eq. ndx  .and. jpe .eq. ndy ) then
             inp = 1
             jnp = 1
              do knp=knp_st,knp_end

               call set_boundary_size(inp, jnp, knp)
               is = nxi+1
               ie = nxi+ndepth
               js = nyi+1
               je = nyi+ndepth

               neibpetot = neibpetot  + 1
               call count_node_id(inod)

               stack_export(neibpetot) = inod

              enddo
            endif

!
!  outdside (x<xmin, y>ymax)
!
            if ( ipe .eq. 1  .and. jpe .eq. ndy ) then
             jnp = 1
              do knp=knp_st,knp_end

               call set_boundary_size(inp, jnp, knp)
               is = ndepth+1
               ie = 2*ndepth
               js = nyi+1
               je = nyi+ndepth

               neibpetot = neibpetot  + 1
               call count_node_id(inod)

               stack_export(neibpetot) = inod

              enddo
            endif
!
!  outdside (x<xmin, y<ymin)
!
            if ( ipe .eq. 1  .and. jpe .eq. 1 ) then
              do knp=knp_st,knp_end

               call set_boundary_size(inp, jnp, knp)
               is = ndepth+1
               ie = 2*ndepth
               js = ndepth+1
               je = 2*ndepth

               neibpetot = neibpetot  + 1
               call count_node_id(inod)

               stack_export(neibpetot) = inod

              enddo
            endif


!  outdside (x>xmax, y<ymin)
!
            if ( ipe .eq. ndx  .and. jpe .eq. 1 ) then
             inp = 1
              do knp=knp_st,knp_end

               call set_boundary_size(inp, jnp, knp)
               is = nxi+1
               ie = nxi+ndepth
               js = ndepth+1
               je = 2*ndepth

               neibpetot = neibpetot  + 1
               call count_node_id(inod)

               stack_export(neibpetot) = inod

              enddo
            endif

          end subroutine count_export_peri_linear
!
! ----------------------------------------------------------------------
!
      subroutine count_export_peri_quad(ipe, jpe, kpe, inod)
!
      use set_comm_edge_4_cube
!
      integer (kind = kint) :: ipe, jpe, kpe
      integer (kind = kint) :: inod
!
      integer (kind = kint) :: nd
      integer (kind = kint) :: inp, jnp, knp
!
!
!  outdside (x>xmax)
!                                     .... count nodes 
!
            if (ipe .eq. ndx) then
             inp = 1
             do knp=knp_st,knp_end
              do jnp=jnp_st,jnp_end
!
               call set_boundary_size(inp, jnp, knp)
               is = nxi+1
               ie = nxi+ndepth

               neibpetot = neibpetot  + 1
               call count_node_id(inod)

               nd = 1
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 2
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 3
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)

               stack_export(neibpetot) = inod

              enddo
             enddo
            endif
!
!    ---   outside wall (x<xmin)
!                                     .... count nodes 
            if (ipe .eq. 1) then
             do knp=knp_st,knp_end
              do jnp=jnp_st,jnp_end

               call set_boundary_size(inp, jnp, knp)
               is = ndepth+1

               neibpetot = neibpetot  + 1
               ie = 2*ndepth
               call count_node_id(inod)

               nd = 1
               ie = 2*ndepth - 1
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 2
               ie = 2*ndepth
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 3
               ie = 2*ndepth
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)

               stack_export(neibpetot) = inod

              enddo
             enddo
            endif
!
!  outdside (y>ymax)
!                                     .... count nodes 
            if ( jpe .eq. ndy ) then
             jnp = 1
             do knp=knp_st,knp_end
              do inp=inp_st,inp_end

               call set_boundary_size(inp, jnp, knp)
               js = nyi+1
               je = nyi+ndepth

               neibpetot = neibpetot  + 1
               call count_node_id(inod)

               nd = 1
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 2
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 3
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)

               stack_export(neibpetot) = inod

                enddo
              enddo
            endif
!
!  outdside (y<ymin)
!                                     .... count nodes 
            if ( jpe .eq. 1 ) then
             do knp=knp_st,knp_end
              do inp=inp_st,inp_end

               call set_boundary_size(inp, jnp, knp)
               js = ndepth+1

               neibpetot = neibpetot  + 1
               je = 2*ndepth
               call count_node_id(inod)

               nd = 1
               je = 2*ndepth
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 2
               je = 2*ndepth - 1
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 3
               je = 2*ndepth
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)

               stack_export(neibpetot) = inod

                enddo
              enddo
            endif
!
!  outdside (x>xmax, y>ymax)
!
            if ( ipe .eq. ndx  .and. jpe .eq. ndy ) then
             inp = 1
             jnp = 1
              do knp=knp_st,knp_end

               call set_boundary_size(inp, jnp, knp)
               is = nxi+1
               ie = nxi+ndepth
               js = nyi+1
               je = nyi+ndepth

               neibpetot = neibpetot  + 1
               call count_node_id(inod)

               nd = 1
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 2
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 3
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)

               stack_export(neibpetot) = inod

              enddo
            endif
!
!  outdside (x<xmin, y>ymax)
!
            if ( ipe .eq. 1  .and. jpe .eq. ndy ) then
             jnp = 1
              do knp=knp_st,knp_end

               call set_boundary_size(inp, jnp, knp)
               is = ndepth+1
               js = nyi+1
               je = nyi+ndepth

               neibpetot = neibpetot  + 1
               ie = 2*ndepth
               call count_node_id(inod)

               nd = 1
               ie = 2*ndepth - 1
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 2
               ie = 2*ndepth
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 3
               ie = 2*ndepth
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)
!
               stack_export(neibpetot) = inod

              enddo
            endif
!
!  outdside (x<xmin, y<ymin)
!
            if ( ipe .eq. 1  .and. jpe .eq. 1 ) then
              do knp=knp_st,knp_end

               call set_boundary_size(inp, jnp, knp)
               is = ndepth+1
               js = ndepth+1

               neibpetot = neibpetot  + 1
               ie = 2*ndepth
               je = 2*ndepth
               call count_node_id(inod)

               nd = 1
               ie = 2*ndepth - 1
               je = 2*ndepth
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 2
               ie = 2*ndepth
               je = 2*ndepth - 1
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 3
               ie = 2*ndepth
               je = 2*ndepth
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)

               stack_export(neibpetot) = inod

              enddo
            endif
!
!  outdside (x>xmax, y<ymin)
!
            if ( ipe .eq. ndx  .and. jpe .eq. 1 ) then
             inp = 1
              do knp=knp_st,knp_end

               call set_boundary_size(inp, jnp, knp)
               is = nxi+1
               ie = nxi+ndepth
               js = ndepth+1

               neibpetot = neibpetot  + 1
               je = 2*ndepth
               call count_node_id(inod)

               nd = 1
               je = 2*ndepth
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 2
               je = 2*ndepth - 1
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 3
               je = 2*ndepth
               call count_ex_edge(kpe, inp, jnp, knp, inod, nd)

               stack_export(neibpetot) = inod

              enddo
            endif

          end subroutine count_export_peri_quad
!
! ----------------------------------------------------------------------
!
      end module count_export_peri
