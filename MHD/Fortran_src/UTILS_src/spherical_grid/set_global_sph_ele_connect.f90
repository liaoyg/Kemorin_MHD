!set_global_sph_ele_connect.f90
!      module set_global_sph_ele_connect
!
!     Written by H. Matsui on July, 2007
!
!      subroutine set_sph_ele_connect_no_pole
!
      module set_global_sph_ele_connect
!
      use m_precision
!
      use m_geometry_parameter
      use m_geometry_data
      use m_spheric_parameter
!
      implicit none
!
      private :: set_sph_ele_connect
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine set_sph_ele_connect_no_pole
!
      use m_geometry_constants
!
      integer(kind = kint) :: icou
!
      numele =  (nidx_global_rtp(1) - 1)                                &
     &        * (nidx_global_rtp(2) - 1)                                &
     &        *  nidx_global_rtp(3)
!
      nnod_4_ele =  num_t_linear
      nnod_4_surf = num_linear_sf
      nnod_4_edge = num_linear_edge
!
      call allocate_element_connection
!
      elmtyp(1:numele) = 331
      nodelm(1:numele) = nnod_4_ele
!
      icou = 0
      call set_sph_ele_connect(icou)
!
      end subroutine set_sph_ele_connect_no_pole
!
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
!
      subroutine set_sph_ele_connect(icou)
!
      integer(kind = kint), intent(inout) :: icou
      integer(kind = kint) :: kr, kt, kp
!
      do kp = 1, nidx_global_rtp(3) - 1
        do kt = 1, nidx_global_rtp(2) - 1
          do kr = 1, nidx_global_rtp(1) - 1
            icou = icou + 1
            iele_global(icou) = icou
            ie(icou,1) = (kr  ) + (kt-1)*nidx_global_rtp(1)             &
     &                  + (kp-1)*nidx_global_rtp(1)*nidx_global_rtp(2)
            ie(icou,2) = (kr  ) + (kt-1)*nidx_global_rtp(1)             &
     &                  + (kp  )*nidx_global_rtp(1)*nidx_global_rtp(2)
            ie(icou,3) = (kr  ) + (kt  )*nidx_global_rtp(1)             &
     &                  + (kp  )*nidx_global_rtp(1)*nidx_global_rtp(2)
            ie(icou,4) = (kr  ) + (kt  )*nidx_global_rtp(1)             &
     &                  + (kp-1)*nidx_global_rtp(1)*nidx_global_rtp(2)
            ie(icou,5) = (kr+1) + (kt-1)*nidx_global_rtp(1)             &
     &                  + (kp-1)*nidx_global_rtp(1)*nidx_global_rtp(2)
            ie(icou,6) = (kr+1) + (kt-1)*nidx_global_rtp(1)             &
     &                  + (kp  )*nidx_global_rtp(1)*nidx_global_rtp(2)
            ie(icou,7) = (kr+1) + (kt  )*nidx_global_rtp(1)             &
     &                  + (kp  )*nidx_global_rtp(1)*nidx_global_rtp(2)
            ie(icou,8) = (kr+1) + (kt  )*nidx_global_rtp(1)             &
     &                  + (kp-1)*nidx_global_rtp(1)*nidx_global_rtp(2)
          end do
        end do
      end do
!
      kp = nidx_global_rtp(3)
        do kt = 1, nidx_global_rtp(2) - 1
          do kr = 1, nidx_global_rtp(1) - 1
            icou = icou + 1
            iele_global(icou) = icou
            ie(icou,1) = (kr  ) + (kt-1)*nidx_global_rtp(1)             &
     &                  + (kp-1)*nidx_global_rtp(1)*nidx_global_rtp(2)
            ie(icou,2) = (kr  ) + (kt-1)*nidx_global_rtp(1)
            ie(icou,3) = (kr  ) + (kt  )*nidx_global_rtp(1)
            ie(icou,4) = (kr  ) + (kt  )*nidx_global_rtp(1)             &
     &                  + (kp-1)*nidx_global_rtp(1)*nidx_global_rtp(2)
            ie(icou,5) = (kr+1) + (kt-1)*nidx_global_rtp(1)             &
     &                  + (kp-1)*nidx_global_rtp(1)*nidx_global_rtp(2)
            ie(icou,6) = (kr+1) + (kt-1)*nidx_global_rtp(1)
            ie(icou,7) = (kr+1) + (kt  )*nidx_global_rtp(1)
            ie(icou,8) = (kr+1) + (kt  )*nidx_global_rtp(1)             &
     &                  + (kp-1)*nidx_global_rtp(1)*nidx_global_rtp(2)
          end do
        end do
!
      end subroutine set_sph_ele_connect
!
! -----------------------------------------------------------------------
!
      end module set_global_sph_ele_connect
