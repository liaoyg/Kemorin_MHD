!>@file   share_spectr_index_data.f90
!!@brief  module share_spectr_index_data
!!
!!@author H. Matsui
!!@date Programmed in Jan., 2014
!
!>@brief routines for parallel spectr data assemble
!!
!!@verbatim
!!      subroutine share_sph_rj_data(ip_org, sph_mesh)
!!        type(sph_mesh_data), intent(inout) :: sph_mesh
!!      subroutine share_r_interpolation_tbl                            &
!!     &         (new_sph_mesh, r_itp, nlayer_ICB_org, nlayer_CMB_org,  &
!!     &          nlayer_ICB_new, nlayer_CMB_new)
!!        type(sph_mesh_data), intent(in) :: new_sph_mesh
!!        type(sph_radial_itp_data), intent(inout) :: r_itp
!!@endverbatim
!!
!
      module share_spectr_index_data
!
      use m_precision
      use m_constants
      use calypso_mpi
!
      use t_SPH_mesh_field_data
      use r_interpolate_marged_sph
!
      implicit none
!
!
! -----------------------------------------------------------------------
!
      contains
!
! -----------------------------------------------------------------------
!
      subroutine share_sph_rj_data(ip_org, sph_mesh)
!
      integer(kind = kint), intent(in) :: ip_org
      type(sph_mesh_data), intent(inout) :: sph_mesh
!
      integer(kind = kint) :: irank_org
!
!
      irank_org = mod(ip_org - 1,nprocs)
!      write(*,*) 'MPI_Bcast irank_sph_rj', ip_org
      call MPI_Bcast(sph_mesh%sph%sph_rj%irank_sph_rj,                  &
     &    itwo, CALYPSO_INTEGER, irank_org, CALYPSO_COMM, ierr_MPI)
!      write(*,*) 'MPI_Bcast nidx_global_rj', ip_org
      call MPI_Bcast(sph_mesh%sph%sph_rj%nidx_global_rj,                &
     &    itwo, CALYPSO_INTEGER, irank_org, CALYPSO_COMM, ierr_MPI)
!      write(*,*) 'MPI_Bcast nnod_rj', ip_org
      call MPI_Bcast(sph_mesh%sph%sph_rj%nnod_rj,                       &
     &    ione, CALYPSO_INTEGER, irank_org, CALYPSO_COMM, ierr_MPI)
!      write(*,*) 'MPI_Bcast nidx_rj', ip_org
      call MPI_Bcast(sph_mesh%sph%sph_rj%nidx_rj,                       &
     &    itwo, CALYPSO_INTEGER, irank_org, CALYPSO_COMM, ierr_MPI)
!      write(*,*) 'MPI_Bcast ist_rj', ip_org
      call MPI_Bcast(sph_mesh%sph%sph_rj%ist_rj,                        &
     &    itwo, CALYPSO_INTEGER, irank_org, CALYPSO_COMM, ierr_MPI)
!      write(*,*) 'MPI_Bcast ied_rj', ip_org
      call MPI_Bcast(sph_mesh%sph%sph_rj%ied_rj,                        &
     &    itwo, CALYPSO_INTEGER, irank_org, CALYPSO_COMM, ierr_MPI)
!
      if(mod(ip_org-1,nprocs) .ne. my_rank) then
        call alloc_type_spheric_param_rj(sph_mesh%sph%sph_rj)
        call alloc_type_sph_1d_index_rj(sph_mesh%sph%sph_rj)
      end if
!
!      write(*,*) 'MPI_Bcast idx_global_rj', ip_org
      call MPI_Bcast(sph_mesh%sph%sph_rj%idx_global_rj,                 &
     &    (2*sph_mesh%sph%sph_rj%nnod_rj),                              &
     &    CALYPSO_INTEGER, irank_org, CALYPSO_COMM, ierr_MPI)
!
!      write(*,*) 'MPI_Bcast radius_1d_rj_r', ip_org
      call MPI_Bcast(sph_mesh%sph%sph_rj%radius_1d_rj_r,                &
     &    sph_mesh%sph%sph_rj%nidx_rj(1),                               &
     &    CALYPSO_REAL, irank_org, CALYPSO_COMM, ierr_MPI)
!      write(*,*) 'MPI_Bcast idx_gl_1d_rj_r', ip_org
      call MPI_Bcast(sph_mesh%sph%sph_rj%idx_gl_1d_rj_r,                &
     &    sph_mesh%sph%sph_rj%nidx_rj(1),                               &
     &    CALYPSO_INTEGER, irank_org, CALYPSO_COMM, ierr_MPI)
!      write(*,*) 'MPI_Bcast idx_gl_1d_rj_j', ip_org
      call MPI_Bcast(sph_mesh%sph%sph_rj%idx_gl_1d_rj_j,                &
     &    (3*sph_mesh%sph%sph_rj%nidx_rj(2)),                           &
     &    CALYPSO_INTEGER, irank_org, CALYPSO_COMM, ierr_MPI)
!
      end subroutine share_sph_rj_data
!
! -----------------------------------------------------------------------
!
      subroutine share_r_interpolation_tbl                              &
     &         (new_sph_mesh, r_itp, nlayer_ICB_org, nlayer_CMB_org,    &
     &          nlayer_ICB_new, nlayer_CMB_new)
!
      type(sph_mesh_data), intent(in) :: new_sph_mesh
!
      integer(kind = kint), intent(inout) :: nlayer_ICB_org
      integer(kind = kint), intent(inout) :: nlayer_CMB_org
      integer(kind = kint), intent(inout) :: nlayer_ICB_new
      integer(kind = kint), intent(inout) :: nlayer_CMB_new
      type(sph_radial_itp_data), intent(inout) :: r_itp
!
!
      call MPI_Bcast(nlayer_ICB_org, ione, CALYPSO_INTEGER, izero,      &
     &    CALYPSO_COMM, ierr_MPI)
      call MPI_Bcast(nlayer_CMB_org, ione, CALYPSO_INTEGER, izero,      &
     &    CALYPSO_COMM, ierr_MPI)
      call MPI_Bcast(nlayer_ICB_new, ione, CALYPSO_INTEGER, izero,      &
     &    CALYPSO_COMM, ierr_MPI)
      call MPI_Bcast(nlayer_CMB_new, ione, CALYPSO_INTEGER, izero,      &
     &    CALYPSO_COMM, ierr_MPI)
!
      if(my_rank .eq. 0) then
        write(*,*) 'nlayer_ICB_org: ', nlayer_ICB_org, nlayer_CMB_org
        write(*,*) 'nlayer_ICB_new: ', nlayer_ICB_new, nlayer_CMB_new
      end if
!
      call MPI_Bcast(r_itp%iflag_same_rgrid, ione, CALYPSO_INTEGER,     &
     &    izero, CALYPSO_COMM, ierr_MPI)
      call MPI_Bcast(new_sph_mesh%sph%sph_rj%nidx_rj(1),                &
     &    ione, CALYPSO_INTEGER, izero, CALYPSO_COMM, ierr_MPI)
      if(my_rank .eq. 0) write(*,*) 'iflag_same_rgrid: ',               &
     &            r_itp%iflag_same_rgrid,                               &
     &            new_sph_mesh%sph%sph_rj%nidx_rj(1)
!
      if(r_itp%iflag_same_rgrid .eq. 0) then
        if(my_rank .ne. 0)  call allocate_radial_itp_tbl                &
     &             (new_sph_mesh%sph%sph_rj%nidx_rj(1), r_itp)
!
        call MPI_Bcast(r_itp%nri_old2new, ione, CALYPSO_INTEGER,        &
     &      izero, CALYPSO_COMM, ierr_MPI)
        call MPI_Bcast(r_itp%kr_inner_domain, ione, CALYPSO_INTEGER,    &
     &      izero, CALYPSO_COMM, ierr_MPI)
        call MPI_Bcast(r_itp%kr_outer_domain, ione, CALYPSO_INTEGER,    &
     &      izero, CALYPSO_COMM, ierr_MPI)
        call MPI_Bcast(r_itp%k_old2new_in, r_itp%nri_old2new,           &
     &      CALYPSO_INTEGER, izero, CALYPSO_COMM, ierr_MPI)
        call MPI_Bcast(r_itp%k_old2new_out, r_itp%nri_old2new,          &
     &      CALYPSO_INTEGER, izero, CALYPSO_COMM, ierr_MPI)
        call MPI_Bcast(r_itp%coef_old2new_in, r_itp%nri_old2new,        &
     &      CALYPSO_REAL, izero, CALYPSO_COMM, ierr_MPI)
      end if
!
      end subroutine share_r_interpolation_tbl
!
! -----------------------------------------------------------------------
!
      end module share_spectr_index_data