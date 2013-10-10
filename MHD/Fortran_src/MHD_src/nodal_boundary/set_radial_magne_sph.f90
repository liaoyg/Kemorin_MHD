!set_radial_magne_sph.f90
!     module set_radial_magne_sph
!
!        programmed by H.Matsui on July 2000 (ver 1.1)
!        modified by H.Matsui on Aug., 2007
!
!      subroutine set_r_magne_sph(l_f, i, j)
!
      module set_radial_magne_sph
!
      use m_precision
!
      implicit none
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine set_r_magne_sph(l_f, i, j)
!
      use m_node_group
      use m_bc_data_list
      use m_bc_data_magne
      use m_geometry_data
      use m_geometry_parameter
      use m_node_phys_address
      use m_node_phys_data
      use m_schmidt_polynomial
!
      integer(kind = kint) :: i, j, k, l_f(3)
!
      integer(kind = kint) :: inod, nd, i_comp
      integer(kind = kint) :: jj,ll,mm
      real ( kind = kreal) :: bmag, dph
!
!
      call allocate_schmidt_polynomial
!
      jj = int( aint( magne_nod%bc_magnitude(j)) )
      ll = int( aint( sqrt(magne_nod%bc_magnitude(j)) ))
      mm = jj - ll
!
      do k=1, bc_istack(i)-bc_istack(i-1)
!
       inod = bc_item(k+bc_istack(i-1))
!
       do nd = 1, 3
         l_f(nd) = l_f(nd) + 1
         ibc_b_id(l_f(nd),nd) = inod
       end do
!
       dth = colatitude(inod) 
       dph = longitude(inod)
       call dschmidt
!
       if (mm.ge.0) then
         bmag = p(mm,ll) * cos( dph*dble(mm) )
       else
         bmag = p(mm,ll) * sin( dph*dble(mm) )
       end if
!
       do nd = 1, 3
         i_comp = iphys%i_magne + nd - 1
         ibc_magne(inod,nd) = 1
         ibc2_magne(inod,nd) = 1
         bc_b_id_apt(l_f(nd),nd)= bmag * xx(inod,1) *a_radius(inod)
         d_nod(inod,i_comp) = bc_b_id_apt(l_f(nd),nd)
       end do
!
      end do
!
      call deallocate_schmidt_polynomial
!
      end subroutine set_r_magne_sph
!
!-----------------------------------------------------------------------
!
      end module set_radial_magne_sph
