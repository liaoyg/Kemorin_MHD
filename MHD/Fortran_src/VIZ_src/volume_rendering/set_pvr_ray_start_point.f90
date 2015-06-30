!set_pvr_ray_start_point.f90
!      module set_pvr_ray_start_point
!
!        programmed by H.Matsui on Aug., 2011
!
!!      subroutine count_each_pvr_ray_start (numnod, numele, numsurf,   &
!!     &          nnod_4_surf, ie_surf, isf_4_ele,                      &
!!     &          x_nod_screen, npixel_x, npixel_y,                     &
!!     &          pixel_point_x, pixel_point_y,num_pvr_surf,            &
!!     &          item_pvr_surf_domain, screen_norm_pvr_domain,         &
!!     &          isurf_xrng_pvr_domain, jsurf_yrng_pvr_domain, ray_vec,&
!!     &          num_pvr_ray, istack_pvr_ray_sf)
!!      subroutine set_each_pvr_ray_start                               &
!!     &         (numnod, numele, numsurf, nnod_4_surf,                 &
!!     &          xx, ie_surf, isf_4_ele, x_nod_screen,                 &
!!     &          npixel_x, npixel_y, pixel_point_x, pixel_point_y,     &
!!     &          num_pvr_surf, item_pvr_surf_domain,                   &
!!     &          screen_norm_pvr_domain, isurf_xrng_pvr_domain,        &
!!     &          jsurf_yrng_pvr_domain, viewpoint_vec, ray_vec,        &
!!     &          istack_pvr_ray_sf, num_pvr_ray, id_pixel_start,       &
!!     &          icount_pvr_trace, isf_pvr_ray_start, xi_pvr_start,    &
!!     &          xx_pvr_start, xx_pvr_ray_start, pvr_ray_dir)
!
      module set_pvr_ray_start_point
!
      use m_precision
!
      use calypso_mpi
      use m_constants
      use m_geometry_constants
!
      implicit  none
!
      private :: cal_coefs_on_surf
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine count_temporal_pvr_ray_start                           &
     &         (num_pvr_surf, screen_norm_pvr_domain,                   &
     &          isurf_xrng_pvr_domain, jsurf_yrng_pvr_domain, ray_vec,  &
     &          ntot_tmp_pvr_ray_sf, istack_tmp_pvr_ray_st)
!
      integer(kind = kint), intent(in) :: num_pvr_surf
      integer(kind = kint), intent(in)                                  &
     &                    :: isurf_xrng_pvr_domain(2,num_pvr_surf)
      integer(kind = kint), intent(in)                                  &
     &                    :: jsurf_yrng_pvr_domain(2,num_pvr_surf)
      real(kind = kreal), intent(in)                                    &
     &                    :: screen_norm_pvr_domain(3,num_pvr_surf)
!
      real(kind = kreal), intent(in) :: ray_vec(3)
!
      integer(kind = kint), intent(inout) :: ntot_tmp_pvr_ray_sf
      integer(kind = kint), intent(inout)                               &
     &                   :: istack_tmp_pvr_ray_st(0:num_pvr_surf)
!
      integer(kind = kint) :: inum
!
!
!$omp parallel do private(inum)
      do inum = 1, num_pvr_surf
        if((screen_norm_pvr_domain(3,inum)*ray_vec(3)) .gt. zero) then
          istack_tmp_pvr_ray_st(inum) = (isurf_xrng_pvr_domain(2,inum)  &
     &                            - isurf_xrng_pvr_domain(1,inum)+1)    &
     &                           * (jsurf_yrng_pvr_domain(2,inum)       &
     &                            - jsurf_yrng_pvr_domain(1,inum)+1)
        else
          istack_tmp_pvr_ray_st(inum) = 0
        end if
      end do
!$omp end parallel do
!
      istack_tmp_pvr_ray_st(0) = 0
      do inum = 1, num_pvr_surf
        istack_tmp_pvr_ray_st(inum) = istack_tmp_pvr_ray_st(inum-1)     &
     &                               + istack_tmp_pvr_ray_st(inum)
      end do
      ntot_tmp_pvr_ray_sf = istack_tmp_pvr_ray_st(num_pvr_surf)
!
      end subroutine count_temporal_pvr_ray_start
!
!  ---------------------------------------------------------------------
!
      subroutine count_each_pvr_ray_start (numnod, numele, numsurf,     &
     &          nnod_4_surf, ie_surf, isf_4_ele,                        &
     &          x_nod_screen, npixel_x, npixel_y,                       &
     &          pixel_point_x, pixel_point_y,num_pvr_surf,              &
     &          item_pvr_surf_domain, screen_norm_pvr_domain,           &
     &          isurf_xrng_pvr_domain, jsurf_yrng_pvr_domain, ray_vec,  &
     &          num_pvr_ray, istack_pvr_ray_sf, ntot_tmp_pvr_ray,       &
     &          istack_tmp_pvr_ray_st, iflag_start_tmp,                 &
     &          xi_pvr_start_tmp)
!
      integer(kind = kint), intent(in) :: numnod, numele, numsurf
      integer(kind = kint), intent(in) :: nnod_4_surf
      integer(kind = kint), intent(in) :: ie_surf(numsurf,nnod_4_surf)
      integer(kind = kint), intent(in) :: isf_4_ele(numele,nsurf_4_ele)
!
      real(kind = kreal), intent(in) :: x_nod_screen(numnod,4)
!
      integer(kind = kint), intent(in) :: npixel_x, npixel_y
      real(kind = kreal), intent(in) :: pixel_point_x(npixel_x)
      real(kind = kreal), intent(in) :: pixel_point_y(npixel_y)
!
      integer(kind = kint), intent(in) :: num_pvr_surf
      integer(kind = kint), intent(in)                                  &
     &                    :: item_pvr_surf_domain(2,num_pvr_surf)
      real(kind = kreal), intent(in)                                    &
     &                    :: screen_norm_pvr_domain(3,num_pvr_surf)
      integer(kind = kint), intent(in)                                  &
     &                    :: isurf_xrng_pvr_domain(2,num_pvr_surf)
      integer(kind = kint), intent(in)                                  &
     &                    :: jsurf_yrng_pvr_domain(2,num_pvr_surf)
!
      real(kind = kreal), intent(in) :: ray_vec(3)
!
      integer(kind = kint), intent(in) :: ntot_tmp_pvr_ray
      integer(kind = kint), intent(in)                                  &
     &                   :: istack_tmp_pvr_ray_st(0:num_pvr_surf)
!
      integer(kind = kint), intent(inout) :: num_pvr_ray
      integer(kind = kint), intent(inout)                               &
     &                   :: iflag_start_tmp(ntot_tmp_pvr_ray)
      real(kind = kreal), intent(inout)                                 &
     &                   :: xi_pvr_start_tmp(2,ntot_tmp_pvr_ray)
      integer(kind = kint), intent(inout)                               &
     &                   :: istack_pvr_ray_sf(0:num_pvr_surf)
!
      integer(kind = kint) :: inum, iele, k1, isurf, iflag, itmp, icou
      integer(kind = kint) :: ist_pix, ied_pix, jst_pix, jed_pix
      integer(kind = kint) :: ipix, jpix, i1, i2, i3, i4
      real(kind = kreal) :: x_pix(2), xi(2)
      real(kind = kreal) :: x_surf(2,4)
!
!
!$omp parallel do private(inum,iele,k1,isurf,iflag,x_surf,i1,i2,i3,i4,  &
!$omp&        icou,ist_pix,ied_pix,jst_pix,jed_pix,ipix,jpix,x_pix,xi)
      do inum = 1, num_pvr_surf
        istack_pvr_ray_sf(inum) = 0
        iele = item_pvr_surf_domain(1,inum)
        k1 =   item_pvr_surf_domain(2,inum)
        isurf = abs(isf_4_ele(iele,k1))
        icou = istack_tmp_pvr_ray_st(inum-1)
!
        if((screen_norm_pvr_domain(3,inum)*ray_vec(3)) .gt. zero) then
          i1 = ie_surf(isurf,1)
          i2 = ie_surf(isurf,2)
          i3 = ie_surf(isurf,3)
          i4 = ie_surf(isurf,4)
          x_surf(1:2,1) = x_nod_screen(i1,1:2)
          x_surf(1:2,2) = x_nod_screen(i2,1:2)
          x_surf(1:2,3) = x_nod_screen(i3,1:2)
          x_surf(1:2,4) = x_nod_screen(i4,1:2)
!
          ist_pix = isurf_xrng_pvr_domain(1,inum)
          ied_pix = isurf_xrng_pvr_domain(2,inum)
          jst_pix = jsurf_yrng_pvr_domain(1,inum)
          jed_pix = jsurf_yrng_pvr_domain(2,inum)
          do ipix = ist_pix, ied_pix
            do jpix = jst_pix, jed_pix
              icou = icou + 1
              x_pix(1) = pixel_point_x(ipix)
              x_pix(2) = pixel_point_y(jpix)
              call cal_coefs_on_surf(x_surf, x_pix,                     &
     &            iflag_start_tmp(icou), xi_pvr_start_tmp(1,icou))
!              write(100+my_rank,*) inum, ipix, jpix, iflag, xi
!
              if(iflag_start_tmp(icou) .gt. 0) then
                istack_pvr_ray_sf(inum) = istack_pvr_ray_sf(inum) + 1
              end if
            end do
          end do
        end if
      end do
!$omp end parallel do
!
      istack_pvr_ray_sf(0) = 0
      do inum = 1, num_pvr_surf
        istack_pvr_ray_sf(inum) = istack_pvr_ray_sf(inum-1)             &
     &                           + istack_pvr_ray_sf(inum)
      end do
      num_pvr_ray = istack_pvr_ray_sf(num_pvr_surf)
!
      end subroutine count_each_pvr_ray_start
!
!  ---------------------------------------------------------------------
!
      subroutine set_each_pvr_ray_start                                 &
     &         (numnod, numele, numsurf, nnod_4_surf,                   &
     &          xx, ie_surf, isf_4_ele, x_nod_screen,                   &
     &          npixel_x, npixel_y, pixel_point_x, pixel_point_y,       &
     &          num_pvr_surf, item_pvr_surf_domain,                     &
     &          screen_norm_pvr_domain, isurf_xrng_pvr_domain,          &
     &          jsurf_yrng_pvr_domain, viewpoint_vec, ray_vec,          &
     &          ntot_tmp_pvr_ray, istack_tmp_pvr_ray_st,                &
     &          iflag_start_tmp, xi_pvr_start_tmp,                      &
     &          istack_pvr_ray_sf, num_pvr_ray, id_pixel_start,         &
     &          icount_pvr_trace, isf_pvr_ray_start, xi_pvr_start,      &
     &          xx_pvr_start, xx_pvr_ray_start, pvr_ray_dir)
!
      use cal_field_on_surf_viz
!
      integer(kind = kint), intent(in) :: numnod, numele, numsurf
      integer(kind = kint), intent(in) :: nnod_4_surf
      integer(kind = kint), intent(in) :: ie_surf(numsurf,nnod_4_surf)
      integer(kind = kint), intent(in) :: isf_4_ele(numele,nsurf_4_ele)
      real(kind = kreal), intent(in)  :: xx(numnod,3)
!
      real(kind = kreal), intent(in) :: x_nod_screen(numnod,4)
!
      integer(kind = kint), intent(in) :: npixel_x, npixel_y
      real(kind = kreal), intent(in) :: pixel_point_x(npixel_x)
      real(kind = kreal), intent(in) :: pixel_point_y(npixel_y)
!
      integer(kind = kint), intent(in) :: num_pvr_surf
      integer(kind = kint), intent(in)                                  &
     &                    :: item_pvr_surf_domain(2,num_pvr_surf)
      integer(kind = kint), intent(in)                                  &
     &                    :: isurf_xrng_pvr_domain(2,num_pvr_surf)
      integer(kind = kint), intent(in)                                  &
     &                    :: jsurf_yrng_pvr_domain(2,num_pvr_surf)
      real(kind = kreal), intent(in)                                    &
     &                    :: screen_norm_pvr_domain(3,num_pvr_surf)
!
      real(kind = kreal), intent(in)  :: viewpoint_vec(3)
      real(kind = kreal), intent(in) :: ray_vec(3)
!
      integer(kind = kint), intent(in) :: ntot_tmp_pvr_ray
      integer(kind = kint), intent(in)                                  &
     &                   :: istack_tmp_pvr_ray_st(0:num_pvr_surf)
      integer(kind = kint), intent(in)                                  &
     &                   :: iflag_start_tmp(ntot_tmp_pvr_ray)
      real(kind = kreal), intent(in)                                    &
     &                   :: xi_pvr_start_tmp(2,ntot_tmp_pvr_ray)
      integer(kind = kint), intent(in)                                  &
     &                    :: istack_pvr_ray_sf(0:num_pvr_surf)
!
      integer(kind = kint), intent(in) ::  num_pvr_ray
      integer(kind = kint), intent(inout)                               &
     &                    :: id_pixel_start(num_pvr_ray)
      integer(kind = kint), intent(inout)                               &
     &                    :: icount_pvr_trace(num_pvr_ray)
      integer(kind = kint), intent(inout)                               &
     &                    :: isf_pvr_ray_start(3,num_pvr_ray)
      real(kind = kreal), intent(inout) :: xi_pvr_start(2,num_pvr_ray)
      real(kind = kreal), intent(inout) :: xx_pvr_start(3,num_pvr_ray)
      real(kind = kreal), intent(inout)                                 &
     &                   :: xx_pvr_ray_start(3,num_pvr_ray)
      real(kind = kreal), intent(inout) :: pvr_ray_dir(3,num_pvr_ray)
!
      integer(kind = kint) :: inum, icou, jcou, iele, k1, isurf
      integer(kind = kint) :: ist_pix, ied_pix, jst_pix, jed_pix
      integer(kind = kint) :: ipix, jpix
!
!
!      write(*,*) 'set_each_pvr_ray_start loop '
!$omp parallel do private(inum,icou,jcou,iele,k1,isurf,                 &
!$omp&                    ipix,jpix,ist_pix,ied_pix,jst_pix,jed_pix)
      do inum = 1, num_pvr_surf
        icou = istack_tmp_pvr_ray_st(inum-1)
        jcou = istack_pvr_ray_sf(inum-1)
        iele = item_pvr_surf_domain(1,inum)
        k1 =   item_pvr_surf_domain(2,inum)
        isurf = abs(isf_4_ele(iele,k1))
!
        if( (screen_norm_pvr_domain(3,inum)*ray_vec(3)) .gt. zero) then
          ist_pix = isurf_xrng_pvr_domain(1,inum)
          ied_pix = isurf_xrng_pvr_domain(2,inum)
          jst_pix = jsurf_yrng_pvr_domain(1,inum)
          jed_pix = jsurf_yrng_pvr_domain(2,inum)
          do ipix = ist_pix, ied_pix
            do jpix = jst_pix, jed_pix
              icou = icou + 1
              if(iflag_start_tmp(icou) .gt. 0) then
                jcou = jcou + 1
                if(jcou .gt. num_pvr_ray) write(*,*) 'aho', my_rank,    &
     &                         jcou, num_pvr_ray, inum, num_pvr_surf
!
                icount_pvr_trace(jcou) = 0
                id_pixel_start(jcou) = ipix + (jpix-1)*npixel_x
                isf_pvr_ray_start(1,jcou) = iele
                isf_pvr_ray_start(2,jcou) = k1
                isf_pvr_ray_start(3,jcou) = ie_surf(isurf,1)
                xi_pvr_start(1:2,jcou) =   xi_pvr_start_tmp(1:2,icou)
                xx_pvr_ray_start(1,jcou) = pixel_point_x(ipix)
                xx_pvr_ray_start(2,jcou) = pixel_point_y(jpix)
!
                call cal_field_on_surf_scalar(numnod, numsurf,          &
     &              nnod_4_surf, ie_surf, isurf, xi_pvr_start(1,jcou),  &
     &              x_nod_screen(1,3), xx_pvr_ray_start(3,jcou) )
                call cal_field_on_surf_vector(numnod, numsurf,          &
     &              nnod_4_surf, ie_surf, isurf, xi_pvr_start(1,jcou),  &
     &              xx(1,1), xx_pvr_start(1,jcou) )
!
                pvr_ray_dir(1,jcou) = viewpoint_vec(1)                  &
     &                                 - xx_pvr_start(1,jcou)
                pvr_ray_dir(2,jcou) = viewpoint_vec(2)                  &
     &                                 - xx_pvr_start(2,jcou)
                pvr_ray_dir(3,jcou) = viewpoint_vec(3)                  &
     &                                 - xx_pvr_start(3,jcou)
              end if
            end do
          end do
        end if
!
      end do
!$omp end parallel do
!       write(*,*) 'set_each_pvr_ray_start end '
!
      end subroutine set_each_pvr_ray_start
!
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!
      subroutine cal_coefs_on_surf(x_surf, x_tgt, iflag_inside, xi)
!
      real(kind = kreal), intent(in) ::    x_tgt(2)
      real(kind = kreal), intent(in) ::    x_surf(2,4)
      integer(kind = kint), intent(inout) :: iflag_inside
      real(kind = kreal), intent(inout) :: xi(2)
!
      real(kind = kreal) :: xt1(2), a(2,2), coef(2,2)
      real(kind = kreal) :: c1, c3, aj
!
!
      xt1(1:2) = x_tgt(1:2) - x_surf(1:2,1)
      a(1:2,1) = x_surf(1:2,2) - x_surf(1:2,1)
      a(1:2,2) = x_surf(1:2,4) - x_surf(1:2,1)
      aj = one / (a(1,1)*a(2,2) - a(2,1)*a(1,2))
!
      coef(1,1) = ( a(2,2)*xt1(1) - a(1,2)*xt1(2) ) * aj
      coef(2,1) = (-a(2,1)*xt1(1) + a(1,1)*xt1(2) ) * aj
!
      xt1(1:2) = x_tgt(1:2) - x_surf(1:2,3)
      a(1:2,1) = x_surf(1:2,2) - x_surf(1:2,3)
      a(1:2,2) = x_surf(1:2,4) - x_surf(1:2,3)
      aj = one / (a(1,1)*a(2,2) - a(2,1)*a(1,2))
!
      coef(1,2) = ( a(2,2)*xt1(1) - a(1,2)*xt1(2) ) * aj
      coef(2,2) = (-a(2,1)*xt1(1) + a(1,1)*xt1(2) ) * aj
!
      c1 = one - (coef(1,1) + coef(2,1))
      c3 = one - (coef(1,2) + coef(2,2))
!
      if(coef(1,1).ge.zero .and. coef(2,1).ge.zero                      &
     &      .and. c1 .ge. zero) then
        iflag_inside = 1
        xi(1) = -one + two*coef(1,1)
        xi(2) = -one + two*coef(2,1)
!
      else if(coef(1,2).ge.zero .and. coef(2,2).ge.zero                 &
     &      .and. c3 .ge. zero) then
        iflag_inside = 1
        xi(1) = one - two*coef(2,2)
        xi(2) = one - two*coef(1,2)
      else
        xi(1) = -two
        xi(2) = -two
        iflag_inside = 0
      end if
!
      end subroutine cal_coefs_on_surf
!
!-----------------------------------------------------------------------
!
      end module set_pvr_ray_start_point
