!set_node_on_edge_quad_psf.f90
!      module set_node_on_edge_quad_psf
!
!      Written by H. Matsui on Apr., 2006
!
!!      subroutine set_node_on_edge_4_quad_psf(numnod, numedge,         &
!!     &          nnod_4_edge, ie_edge, xx, const_psf, nnod_on_edge,    &
!!     &          np_smp, istack_n_on_e_smp, iedge_4_nod, coef_on_edge)
!
!
      module set_node_on_edge_quad_psf
!
      use m_precision
      use m_constants
!
      implicit none
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine set_node_on_edge_4_quad_psf(numnod, numedge,           &
     &          nnod_4_edge, ie_edge, xx, const_psf, nnod_on_edge,      &
     &          np_smp, istack_n_on_e_smp, iedge_4_nod, coef_on_edge)
!
      integer(kind = kint), intent(in) :: numnod, numedge, nnod_4_edge
      integer(kind = kint), intent(in) :: ie_edge(numedge,nnod_4_edge)
      real(kind = kreal), intent(in) :: xx(numnod,3)
!
      integer(kind = kint), intent(in) :: nnod_on_edge, np_smp
      integer(kind = kint), intent(in) :: istack_n_on_e_smp(0:np_smp)
      integer(kind = kint), intent(in) :: iedge_4_nod(nnod_on_edge)
!
      real(kind = kreal), intent(in) :: const_psf(10)
!
      real(kind= kreal), intent(inout) :: coef_on_edge(nnod_on_edge,2)
!
      integer(kind = kint) :: ip, ist, ied, inum
      integer(kind = kint) :: inod1, inod2, iedge, ierr
      real(kind = kreal) :: x_mid, y_mid, z_mid
      real(kind = kreal) :: x_dif, y_dif, z_dif
      real(kind = kreal) :: c_xi2, c_xi1, c_xi0, diag, Xi, X1, X2
!
!
!$omp parallel do                                                       &
!$omp& private(xi,x1,x2,diag,c_xi0,c_xi1,c_xi2,x_dif,y_dif,z_dif,       &
!$omp&         x_mid,y_mid,z_mid,inod1,inod2,iedge,ierr,inum,ist,ied)
      do ip = 1, np_smp
        ist = istack_n_on_e_smp(ip-1) + 1
        ied = istack_n_on_e_smp(ip)
        do inum = ist, ied
          iedge = iedge_4_nod(inum)
          inod1 = ie_edge(iedge,1)
          inod2 = ie_edge(iedge,2)
          x_mid = half*(xx(inod2,1) + xx(inod1,1))
          y_mid = half*(xx(inod2,2) + xx(inod1,2))
          z_mid = half*(xx(inod2,3) + xx(inod1,3))
          x_dif = half*(xx(inod2,1) - xx(inod1,1))
          y_dif = half*(xx(inod2,2) - xx(inod1,2))
          z_dif = half*(xx(inod2,3) - xx(inod1,3))
!
          c_xi2 =  const_psf( 1) * x_dif**2                             &
     &         + const_psf( 2) * y_dif**2                               &
     &         + const_psf( 3) * z_dif**2                               &
     &         + const_psf( 4) * x_dif*y_dif                            &
     &         + const_psf( 5) * y_dif*z_dif                            &
     &         + const_psf( 6) * z_dif*x_dif
!
          c_xi1 =  const_psf( 1) * two*x_dif*x_mid                      &
     &         + const_psf( 2) * two*y_dif*y_mid                        &
     &         + const_psf( 3) * two*z_dif*z_mid                        &
     &         + const_psf( 4) * (x_dif*y_mid + x_mid*y_dif)            &
     &         + const_psf( 5) * (y_dif*z_mid + y_mid*z_dif)            &
     &         + const_psf( 6) * (z_dif*x_mid + z_mid*x_dif)            &
     &         + const_psf( 7) * x_dif                                  &
     &         + const_psf( 8) * y_dif                                  &
     &         + const_psf( 9) * z_dif
!
          c_xi0 =  const_psf( 1) * x_mid**2                             &
     &         + const_psf( 2) * y_mid**2                               &
     &         + const_psf( 3) * z_mid**2                               &
     &         + const_psf( 4) * x_mid*y_mid                            &
     &         + const_psf( 5) * y_mid*z_mid                            &
     &         + const_psf( 6) * z_mid*x_mid                            &
     &         + const_psf( 7) * y_mid                                  &
     &         + const_psf( 8) * y_mid                                  &
     &         + const_psf( 9) * z_mid                                  &
     &         + const_psf(10)
!
          if(c_xi2 .ne. zero) then
            diag = c_xi1*c_xi1 - four*c_xi2*c_xi0
!
            if(diag .ge. zero) then
              X1 = half * (-c_xi1 + sqrt(diag) ) / c_xi2
              X2 = half * (-c_xi1 - sqrt(diag) ) / c_xi2
            else
              ierr = 1
              write(*,*) 'Detarminant is  wrong in solution'
              cycle
            end if
!
            if(X1 .ge. -one .and. X1 .le. one) then
              xi = X1
            else if(X2 .ge. -one .and. X2 .le. one) then
              xi = X2
            else
              ierr = 1
              write(*,*) 'PSF solution is  wrong'
              cycle
            end if
          else if(c_xi1 .ne. zero) then
            xi = -c_xi0 / c_xi1
          else
            write(*,*) 'PSF solution is  wrong'
            cycle
          end if
!
          coef_on_edge(inum,1) = half*(one - xi)
          coef_on_edge(inum,2) = half*(one + xi)
        end do
      end do
!$omp end parallel do
!
      end subroutine set_node_on_edge_4_quad_psf
!
!  ---------------------------------------------------------------------
!
      end module set_node_on_edge_quad_psf
