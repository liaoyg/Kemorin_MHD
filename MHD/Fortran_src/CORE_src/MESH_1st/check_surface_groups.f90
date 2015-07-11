!check_surface_groups.f90
!      module check_surface_groups
!
!      Written by H> Matsui on Aug., 2006
!
!      subroutine check_center_of_surface_grp(id_check)
!      subroutine check_center_of_surface_grp_sph(id_check)
!      subroutine check_surface_param_smp(txt, id_check)
!
!       subroutine check_norm_surface_grp(id_check)
!
      module check_surface_groups
!
      use m_precision
!
      implicit  none
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine check_center_of_surface_grp(id_check)
!
      use m_surface_group
      use m_surface_group_geometry
!
      integer(kind = kint), intent(in) :: id_check
!
      integer(kind = kint) :: i_grp, ist, ied, inum
!
      write(id_check,*) ' inum, center of surface'
      do i_grp = 1, sf_grp1%num_grp
        ist = sf_grp1%istack_grp(i_grp-1) + 1
        ied = sf_grp1%istack_grp(i_grp)
        do inum = ist, ied
          write(id_check,'(i16,1p3e23.12)') inum,                       &
     &              sf_grp_v1%x_sf_grp(inum,1:3)
        end do
      end do
!
      end subroutine check_center_of_surface_grp
!
!-----------------------------------------------------------------------
!
      subroutine check_center_of_surface_grp_sph(id_check)
!
      use m_surface_group
      use m_surface_group_geometry
!
      integer(kind = kint), intent(in) :: id_check
!
      integer(kind = kint) :: i_grp, ist, ied, inum
!
      write(id_check,*) ' inum, center of surface'
      write(id_check,*) '          (r, theta, phi, cyl_r)'
      do i_grp = 1, sf_grp1%num_grp
        ist = sf_grp1%istack_grp(i_grp-1) + 1
        ied = sf_grp1%istack_grp(i_grp)
        do inum = ist, ied
          write(id_check,'(i16,1p4e23.12)') inum,                       &
     &        sf_grp_v1%r_sf_grp(inum), sf_grp_v1%theta_sf_grp(inum),   &
     &        sf_grp_v1%phi_sf_grp(inum), sf_grp_v1%s_sf_grp(inum)
        end do
      end do
!
      end subroutine check_center_of_surface_grp_sph
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
       subroutine check_surface_param_smp(txt, id_check)
!
       use m_machine_parameter
       use m_surface_group
       use m_surface_group_connect
!
       integer(kind = kint), intent(in) :: id_check
       character(*), intent(in) :: txt
       integer(kind = kint) :: i, isurf
!
       write(50+id_check,*) txt
       write(50+id_check,*) 'surf_istack'
       write(50+id_check,*) sf_grp1%istack_grp(0:sf_grp1%num_grp)
       write(50+id_check,*) 'inod_stack_sf_grp'
       write(50+id_check,*)                                             &
     &     sf_grp_nod1%inod_stack_sf_grp(0:sf_grp1%num_grp)
!
       write(50+id_check,*) 'isurf_grp_smp_stack'
       do isurf = 1, sf_grp1%num_grp
         write(50+id_check,*) isurf,                                    &
     &    (isurf_grp_smp_stack(i),i=(np_smp*(isurf-1)+1),np_smp*isurf)
       end do
      write(50+id_check,*) 'isurf_nod_smp_stack'
      call check_surf_nod_4_sheard_para                                 &
     &   (id_check, sf_grp1%num_grp, sf_grp_nod1)
!
       end subroutine check_surface_param_smp
!
!-----------------------------------------------------------------------
!
       subroutine check_norm_surface_grp(id_check)
!
       use m_geometry_data
       use m_surface_group
       use m_surface_group_geometry
!
       integer(kind = kint), intent(in) :: id_check
       integer(kind = kint) :: isurf, iele
!
!
       do isurf = 1, sf_grp1%num_item
         iele = sf_grp1%item_sf_grp(1,isurf)
         write(id_check,'(2i16, 1p6E25.15e3)') isurf, iele,             &
     &           sf_grp_v1%vnorm_sf_grp(isurf,1:3),                     &
     &           sf_grp_v1%area_sf_grp(isurf),                          &
     &           volume_ele(iele)
       end do
!
       end subroutine check_norm_surface_grp
!
!-----------------------------------------------------------------------
!
      end module check_surface_groups
