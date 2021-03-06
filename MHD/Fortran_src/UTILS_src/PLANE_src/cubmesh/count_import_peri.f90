!count_import_peri.f90
!     module count_import_peri
!
!     Written by H. Matsui
!     modified by H. Matsui on Aug., 2007
!
!      subroutine count_import_peri_linear(ipe, jpe, kpe, inod)
!      subroutine count_import_peri_quad(ipe, jpe, kpe, inod)
!
      module count_import_peri
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
      subroutine count_import_peri_linear(ipe, jpe, kpe, inod)
!
      integer (kind = kint) :: ipe, jpe, kpe
      integer (kind = kint) :: inod
!
      integer (kind = kint) :: inp, jnp, knp
!
!    ---   outside wall (x<xmin)
!                                     .... count nodes 
            if (ipe .eq. 1) then
             do knp=knp_st,knp_end
              do jnp=jnp_st,jnp_end

               call set_sleeve_size(inp, jnp, knp)
               is = 1
               ie = ndepth

               neibpetot = neibpetot  + 1
               call count_node_id(inod)

               stack_import(neibpetot) = inod

              enddo
             enddo
            endif
!
!  outdside (x>xmax)
!                                     .... count nodes 
            if (ipe .eq. ndx) then
             do knp=knp_st,knp_end
              do jnp=jnp_st,jnp_end
!
               call set_sleeve_size(inp, jnp, knp)
                is = nxi+ndepth+1
                ie = nxi+2*ndepth

               neibpetot = neibpetot  + 1
               call count_node_id(inod)

               stack_import(neibpetot) = inod

              enddo
             enddo
            endif
!
!  outdside (y<ymin)
!                                     .... count nodes 
            if ( jpe .eq. 1 ) then
             do knp=knp_st,knp_end
              do inp=inp_st,inp_end

               call set_sleeve_size(inp, jnp, knp)
               js = 1
               je = ndepth

               neibpetot = neibpetot  + 1
               call count_node_id(inod)

               stack_import(neibpetot) = inod

              enddo
             enddo
            endif
!
!  outdside (y<ymax)
!                                     .... count nodes 
            if ( jpe .eq. ndy ) then
             do knp=knp_st,knp_end
              do inp=inp_st,inp_end

               call set_sleeve_size(inp, jnp, knp)
               js = nyi+ndepth+1
               je = nyi+2*ndepth

               neibpetot = neibpetot  + 1
               call count_node_id(inod)

               stack_import(neibpetot) = inod

              enddo
             enddo
            endif
!
!  outdside (x<xmin, y<ymin)
!
            if ( ipe .eq. 1  .and. jpe .eq. 1 ) then
              do knp=knp_st,knp_end

               call set_sleeve_size(inp, jnp, knp)
               is = 1
               ie = ndepth
               js = 1
               je = ndepth

               neibpetot = neibpetot  + 1
               call count_node_id(inod)

               stack_import(neibpetot) = inod

              enddo
            endif
!
!  outdside (x>xmax, y<ymin)
!
            if ( ipe .eq. ndx  .and. jpe .eq. 1 ) then
              do knp=knp_st,knp_end

               call set_sleeve_size(inp, jnp, knp)
               is = nxi+ndepth+1
               ie = nxi+2*ndepth
               js = 1
               je = ndepth

               neibpetot = neibpetot  + 1
               call count_node_id(inod)

               stack_import(neibpetot) = inod

              enddo
            endif
!
!  outdside (x>xmax, y<ymax)
!
            if ( ipe .eq. ndx  .and. jpe .eq. ndy ) then
              do knp=knp_st,knp_end

               call set_sleeve_size(inp, jnp, knp)
               is = nxi+ndepth+1
               ie = nxi+2*ndepth
               js = nyi+ndepth+1
               je = nyi+2*ndepth

               neibpetot = neibpetot  + 1
               call count_node_id(inod)

               stack_import(neibpetot) = inod

              enddo
            endif
!
!  outdside (x>xmin, y<ymax)
!
            if ( ipe .eq. 1  .and. jpe .eq. ndy ) then
              do knp=knp_st,knp_end

               call set_sleeve_size(inp, jnp, knp)
               is = 1
               ie = ndepth
               js = nyi+ndepth+1
               je = nyi+2*ndepth

               neibpetot = neibpetot  + 1
               call count_node_id(inod)

               stack_import(neibpetot) = inod

              enddo
            endif
!
      end subroutine  count_import_peri_linear
!
! ----------------------------------------------------------------------
!
      subroutine count_import_peri_quad(ipe, jpe, kpe, inod)
!
      use set_comm_edge_4_cube
!
      integer (kind = kint) :: ipe, jpe, kpe
      integer (kind = kint) :: inod
!
      integer (kind = kint) :: inp, jnp, knp
      integer (kind = kint) :: nd
!
!    ---   outside wall (x<xmin)
!                                     .... count nodes 
            if (ipe .eq. 1) then
             inp = -1
             do knp=knp_st,knp_end
              do jnp=jnp_st,jnp_end

               call set_sleeve_size(inp, jnp, knp)
               is = 1
               ie = ndepth

               neibpetot = neibpetot  + 1
               call count_node_id(inod)

               nd = 1
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 2
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 3
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               stack_import(neibpetot) = inod

              enddo
             enddo
            endif
!
!  outdside (x>xmax)
!                                     .... count nodes 
!
            if (ipe .eq. ndx) then
             do knp=knp_st,knp_end
              do jnp=jnp_st,jnp_end
!
               call set_sleeve_size(inp, jnp, knp)
               is = nxi+ndepth+1

               neibpetot = neibpetot  + 1
               ie = nxi+2*ndepth
               call count_node_id(inod)

               nd = 1
               ie = nxi+2*ndepth - 1
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 2
               ie = nxi+2*ndepth
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 3
               ie = nxi+2*ndepth
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               stack_import(neibpetot) = inod

              enddo
             enddo
            endif
!
!  outdside (y<ymin)
!                                     .... count nodes 
            if ( jpe .eq. 1 ) then
             jnp = -1
             do knp=knp_st,knp_end
              do inp=inp_st,inp_end

               call set_sleeve_size(inp, jnp, knp)
               js = 1
               je = ndepth

               neibpetot = neibpetot  + 1
               call count_node_id(inod)

               nd = 1
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 2
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 3
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               stack_import(neibpetot) = inod

              enddo
             enddo
            endif
!
!  outdside (y>ymax)
!                                     .... count nodes 
            if ( jpe .eq. ndy ) then
             do knp=knp_st,knp_end
              do inp=inp_st,inp_end

               call set_sleeve_size(inp, jnp, knp)
               js = nyi+ndepth+1

               neibpetot = neibpetot  + 1
               je = nyi+2*ndepth
               call count_node_id(inod)

               nd = 1
               je = nyi+2*ndepth
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 2
               je = nyi+2*ndepth - 1
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 3
               je = nyi+2*ndepth
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               stack_import(neibpetot) = inod

              enddo
             enddo
            endif
!
!  outdside (x<xmin, y<ymin)
!
            if ( ipe .eq. 1  .and. jpe .eq. 1 ) then
             inp = -1
             jnp = -1
              do knp=knp_st,knp_end

               call set_sleeve_size(inp, jnp, knp)
               is = 1
               ie = ndepth
               js = 1
               je = ndepth

               neibpetot = neibpetot  + 1
               call count_node_id(inod)

               nd = 1
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 2
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 3
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               stack_import(neibpetot) = inod

              enddo
            endif
!
!  outdside (x>xmax, y<ymin)
!
            if ( ipe .eq. ndx  .and. jpe .eq. 1 ) then
             jnp = -1
              do knp=knp_st,knp_end

               call set_sleeve_size(inp, jnp, knp)
               js = 1
               je = ndepth
               is = nxi+ndepth+1

               neibpetot = neibpetot  + 1
               ie = nxi+2*ndepth
               call count_node_id(inod)

               nd = 1
               ie = nxi+2*ndepth - 1
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 2
               ie = nxi+2*ndepth
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 3
               ie = nxi+2*ndepth
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               stack_import(neibpetot) = inod

              enddo
            endif
!
!  outdside (x>xmax, y<ymax)
!
            if ( ipe .eq. ndx  .and. jpe .eq. ndy ) then
              do knp=knp_st,knp_end

               call set_sleeve_size(inp, jnp, knp)
               is = nxi+ndepth+1
               js = nyi+ndepth+1

               neibpetot = neibpetot  + 1
               ie = nxi+2*ndepth
               je = nyi+2*ndepth
               call count_node_id(inod)

               nd = 1
               ie = nxi+2*ndepth - 1
               je = nyi+2*ndepth
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 2
               ie = nxi+2*ndepth
               je = nyi+2*ndepth - 1
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 3
               ie = nxi+2*ndepth
               je = nyi+2*ndepth
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               stack_import(neibpetot) = inod

              enddo
            endif
!
!  outdside (x>xmin, y<ymax)
!
            if ( ipe .eq. 1  .and. jpe .eq. ndy ) then
             inp = -1
              do knp=knp_st,knp_end

               call set_sleeve_size(inp, jnp, knp)
               is = 1
               ie = ndepth
               js = nyi+ndepth+1

               neibpetot = neibpetot  + 1
               je = nyi+2*ndepth
               call count_node_id(inod)

               nd = 1
               je = nyi+2*ndepth
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 2
               je = nyi+2*ndepth - 1
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               nd = 3
               je = nyi+2*ndepth
               call count_im_edge(kpe, inp, jnp, knp, inod, nd)

               stack_import(neibpetot) = inod

              enddo
            endif
!
          end subroutine count_import_peri_quad
!
! ----------------------------------------------------------------------
!
      end module count_import_peri
