!>@file   t_sph_field_FFTW.F90
!!@brief  module t_sph_field_FFTW
!!
!!@author H. Matsui
!!@date Programmed in Apr., 2013
!!@n    Modified on Oct., 2014
!
!>@brief  Fourier transform using FFTW Ver.3
!!
!!@verbatim
!! ------------------------------------------------------------------
!!      subroutine init_sph_field_FFTW(nidx_rtp, irt_rtp_smp_stack,     &
!!     &          FFTW_t)
!!      subroutine finalize_sph_field_FFTW(FFTW_t)
!!      subroutine verify_sph_field_FFTW(nnod_rtp, nidx_rtp,            &
!!     &          irt_rtp_smp_stack, FFTW_t)
!!
!!   wrapper subroutine for initierize FFT by FFTW
!! ------------------------------------------------------------------
!!
!!      subroutine sph_field_fwd_FFTW_to_send                           &
!!     &         (nnod_rtp, nidx_rtp, irt_rtp_smp_stack, ncomp,         &
!!     &          n_WS, irev_sr_rtp, X_rtp, WS, FFTW_t)
!! ------------------------------------------------------------------
!!
!! wrapper subroutine for forward Fourier transform by FFTW3
!!
!!   a_{k} = \frac{2}{Nfft} \sum_{j=0}^{Nfft-1} x_{j} \cos (\frac{2\pi j k}{Nfft})
!!   b_{k} = \frac{2}{Nfft} \sum_{j=0}^{Nfft-1} x_{j} \cos (\frac{2\pi j k}{Nfft})
!!
!!   a_{0} = \frac{1}{Nfft} \sum_{j=0}^{Nfft-1} x_{j}
!!    K = Nfft/2....
!!   a_{k} = \frac{1}{Nfft} \sum_{j=0}^{Nfft-1} x_{j} \cos (\frac{2\pi j k}{Nfft})
!!
!! ------------------------------------------------------------------
!!
!!      subroutine sph_field_back_FFTW_from_recv                        &
!!     &         (nnod_rtp, nidx_rtp, irt_rtp_smp_stack, ncomp,         &
!!     &          n_WR, irev_sr_rtp, WR, X_rtp, FFTW_t)
!! ------------------------------------------------------------------
!!
!! wrapper subroutine for backward Fourier transform by FFTW3
!!
!!   x_{k} = a_{0} + (-1)^{j} a_{Nfft/2} + sum_{k=1}^{Nfft/2-1}
!!          (a_{k} \cos(2\pijk/Nfft) + b_{k} \sin(2\pijk/Nfft))
!!
!! ------------------------------------------------------------------
!!
!!       i = 1:     a_{0}
!!       i = 2:     a_{Nfft/2}
!!       i = 3:     a_{1}
!!       i = 4:     b_{1}
!!       ...
!!       i = 2*k+1: a_{k}
!!       i = 2*k+2: b_{k}
!!       ...
!!       i = Nfft-1:   a_{Nfft/2-1}
!!       i = Nfft:     b_{Nfft/2-1}
!!
!! ------------------------------------------------------------------
!!@endverbatim
!!
!!@n @param Nstacksmp(0:np_smp)   End number for each SMP process
!!@n @param Ncomp           Number of components for Fourier transforms
!!@n @param Nfft        Data length for eadh FFT
!!@n @param X(Ncomp, Nfft)  Data for Fourier transform
!
      module t_sph_field_FFTW
!
      use m_precision
      use m_constants
      use m_machine_parameter
      use t_sph_single_FFTW
!
      implicit none
!
      real(kind = kreal) :: elapsed_fftw(3) = (/0.0,0.0,0.0/)
!
      private :: alloc_fld_FFTW_plan
!
! ------------------------------------------------------------------
!
      contains
!
! ------------------------------------------------------------------
!
      subroutine init_sph_field_FFTW(nidx_rtp, irt_rtp_smp_stack,       &
     &          FFTW_t)
!
      integer(kind = kint), intent(in) :: nidx_rtp(3)
      integer(kind = kint), intent(in) :: irt_rtp_smp_stack(0:np_smp)
      type(work_for_sgl_FFTW), intent(inout) :: FFTW_t
!
      integer(kind = kint) :: ip, ist
      integer(kind = 4) :: Nfft4, howmany, idist_r, idist_c
!
      integer, parameter :: IONE_4 = 1
      integer, parameter :: inembed = 0
      integer, parameter :: istride = 1
!
!
!
      Nfft4 = int(nidx_rtp(3))
      call alloc_fld_FFTW_plan                                          &
     &   (irt_rtp_smp_stack(np_smp), nidx_rtp(3), FFTW_t)
!
      do ip = 1, np_smp
        ist = irt_rtp_smp_stack(ip-1) + 1
        howmany = int(irt_rtp_smp_stack(ip  )                           &
     &           - irt_rtp_smp_stack(ip-1))
        idist_r = int(nidx_rtp(3))
        idist_c = int(nidx_rtp(3)/2+1)
!
#ifdef FFTW3_C
        call kemo_fftw_plan_many_dft_r2c                                &
     &     (FFTW_t%plan_fwd(ip), IONE_4, Nfft4, howmany,                &
     &      FFTW_t%X(1,ist), inembed, istride, idist_r,                 &
     &      FFTW_t%C(1,ist), inembed, istride, idist_c, FFTW_ESTIMATE)
        call kemo_fftw_plan_many_dft_c2r                                &
     &     (FFTW_t%plan_bwd(ip), IONE_4, Nfft4, howmany,                &
     &      FFTW_t%C(1,ist), inembed, istride, idist_c,                 &
     &      FFTW_t%X(1,ist), inembed, istride, idist_r, FFTW_ESTIMATE)
#else
        call dfftw_plan_many_dft_r2c                                    &
     &     (FFTW_t%plan_fwd(ip), IONE_4, Nfft4, howmany,                &
     &      FFTW_t%X(1,ist), inembed, istride, idist_r,                 &
     &      FFTW_t%C(1,ist), inembed, istride, idist_c, FFTW_ESTIMATE)
        call dfftw_plan_many_dft_c2r                                    &
     &     (FFTW_t%plan_bwd(ip), IONE_4, Nfft4, howmany,                &
     &      FFTW_t%C(1,ist), inembed, istride, idist_c,                 &
     &      FFTW_t%X(1,ist), inembed, istride, idist_r, FFTW_ESTIMATE)
#endif
      end do
      FFTW_t%aNfft = one / dble(nidx_rtp(3))
!
      end subroutine init_sph_field_FFTW
!
! ------------------------------------------------------------------
!
      subroutine finalize_sph_field_FFTW(FFTW_t)
!
      type(work_for_sgl_FFTW), intent(inout) :: FFTW_t
!
      integer(kind = kint) :: j
!
!
#ifdef FFTW3_C
      do j = 1, np_smp
        call kemo_fftw_destroy_plan(FFTW_t%plan_fwd(j))
        call kemo_fftw_destroy_plan(FFTW_t%plan_bwd(j))
        call kemo_fftw_cleanup
      end do
#else
      do j = 1, np_smp
        call dfftw_destroy_plan(FFTW_t%plan_fwd(j))
        call dfftw_destroy_plan(FFTW_t%plan_bwd(j))
        call dfftw_cleanup
      end do
#endif
!
      call dealloc_FFTW_plan(FFTW_t)
!
      end subroutine finalize_sph_field_FFTW
!
! ------------------------------------------------------------------
!
      subroutine verify_sph_field_FFTW(nnod_rtp, nidx_rtp,              &
     &          irt_rtp_smp_stack, FFTW_t)
!
      integer(kind = kint), intent(in) :: nnod_rtp
      integer(kind = kint), intent(in) :: nidx_rtp(3)
      integer(kind = kint), intent(in) :: irt_rtp_smp_stack(0:np_smp)
      type(work_for_sgl_FFTW), intent(inout) :: FFTW_t
!
!
      if(allocated(FFTW_t%X) .eqv. .false.) then
        call init_sph_field_FFTW(nidx_rtp, irt_rtp_smp_stack, FFTW_t)
        return
      end if
!
      if(size(FFTW_t%X) .ne. nnod_rtp) then
        call finalize_sph_field_FFTW(FFTW_t)
        call init_sph_field_FFTW(nidx_rtp, irt_rtp_smp_stack, FFTW_t)
      end if
!
      end subroutine verify_sph_field_FFTW
!
! ------------------------------------------------------------------
! ------------------------------------------------------------------
!
      subroutine sph_field_fwd_FFTW_to_send                             &
     &         (nnod_rtp, nidx_rtp, irt_rtp_smp_stack, ncomp,           &
     &          n_WS, irev_sr_rtp, X_rtp, WS, FFTW_t)
!
      integer(kind = kint), intent(in) :: nnod_rtp
      integer(kind = kint), intent(in) :: nidx_rtp(3)
      integer(kind = kint), intent(in) :: irt_rtp_smp_stack(0:np_smp)
!
      integer(kind = kint), intent(in) :: ncomp
      real(kind = kreal), intent(in)                                    &
     &         :: X_rtp(irt_rtp_smp_stack(np_smp),nidx_rtp(3),ncomp)
!
      integer(kind = kint), intent(in) :: n_WS
      integer(kind = kint), intent(in) :: irev_sr_rtp(nnod_rtp)
      real (kind=kreal), intent(inout):: WS(n_WS)
      type(work_for_sgl_FFTW), intent(inout) :: FFTW_t
!
      integer(kind = kint) ::  m, j, ip, ist, ied, nd
      integer(kind = kint) :: ic_rtp, is_rtp, ic_send, is_send
!      real :: dummy(3), rtmp(3)
!
!
!      call cpu_time(dummy(1))
!$omp parallel do schedule(static)                                      &
!$omp&         private(nd,m,j,ip,ist,ied,ic_rtp,is_rtp,ic_send,is_send)
      do ip = 1, np_smp
        ist = irt_rtp_smp_stack(ip-1) + 1
        ied = irt_rtp_smp_stack(ip)
        do nd = 1, ncomp
          do j = ist, ied
            FFTW_t%X(1:nidx_rtp(3),j) = X_rtp(j,1:nidx_rtp(3),nd)
          end do
!        call cpu_time(rtmp(1))
!
!        call cpu_time(dummy(2))
#ifdef FFTW3_C
          call kemo_fftw_execute_dft_r2c(FFTW_t%plan_fwd(ip),           &
     &        FFTW_t%X(1,ist), FFTW_t%C(1,ist))
#else
          call dfftw_execute_dft_r2c(FFTW_t%plan_fwd(ip),               &
     &        FFTW_t%X(1,ist), FFTW_t%C(1,ist))
#endif
!
!      call cpu_time(rtmp(2))
!   normalization
!      call cpu_time(dummy(3))
          do j = ist, ied
            ic_send = nd + (irev_sr_rtp(j) - 1) * ncomp
            WS(ic_send) = FFTW_t%aNfft * real(FFTW_t%C(1,j))
            do m = 2, (nidx_rtp(3)+1)/2
              ic_rtp = j + (2*m-2) * irt_rtp_smp_stack(np_smp)
              is_rtp = j + (2*m-1) * irt_rtp_smp_stack(np_smp)
              ic_send = nd + (irev_sr_rtp(ic_rtp) - 1) * ncomp
              is_send = nd + (irev_sr_rtp(is_rtp) - 1) * ncomp
              WS(ic_send) = two * FFTW_t%aNfft * real(FFTW_t%C(m,j))
              WS(is_send) = two * FFTW_t%aNfft * real(FFTW_t%C(m,j)*iu)
            end do 
            m = (nidx_rtp(3)+1)/2 + 1
            ic_rtp = j + irt_rtp_smp_stack(np_smp)
            ic_send = nd + (irev_sr_rtp(ic_rtp) - 1) * ncomp
            WS(ic_send) = two * FFTW_t%aNfft * real(FFTW_t%C(m,j))
          end do
!      call cpu_time(rtmp(3))
        end do
      end do
!$omp end parallel do
!      elapsed_fftw(1:3) = elapsed_fftw(1:3) + rtmp(1:3) - dummy(1:3)
!
      end subroutine sph_field_fwd_FFTW_to_send
!
! ------------------------------------------------------------------
!
      subroutine sph_field_back_FFTW_from_recv                          &
     &         (nnod_rtp, nidx_rtp, irt_rtp_smp_stack, ncomp,           &
     &          n_WR, irev_sr_rtp, WR, X_rtp, FFTW_t)
!
      integer(kind = kint), intent(in) :: nnod_rtp
      integer(kind = kint), intent(in) :: nidx_rtp(3)
      integer(kind = kint), intent(in) :: irt_rtp_smp_stack(0:np_smp)
!
!
      integer(kind = kint), intent(in) :: ncomp
      integer(kind = kint), intent(in) :: n_WR
      integer(kind = kint), intent(in) :: irev_sr_rtp(nnod_rtp)
      real (kind=kreal), intent(inout):: WR(n_WR)
!
      real(kind = kreal), intent(inout)                                 &
     &          :: X_rtp(irt_rtp_smp_stack(np_smp),nidx_rtp(3),ncomp)
      type(work_for_sgl_FFTW), intent(inout) :: FFTW_t
!
      integer(kind = kint) :: m, j, ip, ist, ied, nd
      integer(kind = kint) :: ic_rtp, is_rtp, ic_recv, is_recv
!      real :: dummy(3), rtmp(3)
!
!
!      call cpu_time(dummy(3))
!$omp parallel do schedule(static)                                      &
!$omp&         private(nd,m,j,ip,ist,ied,ic_rtp,is_rtp,ic_recv,is_recv)
      do ip = 1, np_smp
        ist = irt_rtp_smp_stack(ip-1) + 1
        ied = irt_rtp_smp_stack(ip)
        do nd = 1, ncomp
!   normalization
          do j = ist, ied
            ic_recv = nd + (irev_sr_rtp(j) - 1) * ncomp
            FFTW_t%C(1,j) = cmplx(WR(ic_recv), zero, kind(0d0))
            do m = 2, (nidx_rtp(3)+1)/2
              ic_rtp = j + (2*m-2) * irt_rtp_smp_stack(np_smp)
              is_rtp = j + (2*m-1) * irt_rtp_smp_stack(np_smp)
              ic_recv = nd + (irev_sr_rtp(ic_rtp) - 1) * ncomp
              is_recv = nd + (irev_sr_rtp(is_rtp) - 1) * ncomp
              FFTW_t%C(m,j)                                             &
     &            = half * cmplx(WR(ic_recv), -WR(is_recv),kind(0d0))
            end do
            m = (nidx_rtp(3)+1)/2 + 1
            ic_rtp = j + irt_rtp_smp_stack(np_smp)
            ic_recv = nd + (irev_sr_rtp(ic_rtp) - 1) * ncomp
            FFTW_t%C(m,j) = half * cmplx(WR(ic_recv), zero, kind(0d0))
          end do
!        call cpu_time(rtmp(3))
!
!        call cpu_time(dummy(2))
#ifdef FFTW3_C
          call kemo_fftw_execute_dft_c2r(FFTW_t%plan_bwd(ip),           &
     &        FFTW_t%C(1,ist), FFTW_t%X(1,ist))
#else
          call dfftw_execute_dft_c2r(FFTW_t%plan_bwd(ip),               &
     &        FFTW_t%C(1,ist), FFTW_t%X(1,ist))
#endif
!        call cpu_time(rtmp(2))
!
!        call cpu_time(dummy(1))
          do j = ist, ied
            X_rtp(j,1:nidx_rtp(3),nd) = FFTW_t%X(1:nidx_rtp(3),j)
          end do
        end do
!      call cpu_time(rtmp(1))
      end do
!$omp end parallel do
!
!      elapsed_fftw(1:3) = elapsed_fftw(1:3) + rtmp(1:3) - dummy(1:3)
!
      end subroutine sph_field_back_FFTW_from_recv
!
! ------------------------------------------------------------------
! ------------------------------------------------------------------
!
      subroutine alloc_fld_FFTW_plan(nnod_rt, Nfft, FFTW_t)
!
      integer(kind = kint), intent(in) :: nnod_rt, Nfft
      type(work_for_sgl_FFTW), intent(inout) :: FFTW_t
!
!
      allocate(FFTW_t%plan_bwd(np_smp))
      allocate(FFTW_t%plan_fwd(np_smp))
!
      allocate(FFTW_t%X(Nfft,nnod_rt))
      allocate(FFTW_t%C(Nfft/2+1,nnod_rt))
      FFTW_t%X = 0.0d0
      FFTW_t%C = 0.0d0
!
      end subroutine alloc_fld_FFTW_plan
!
! ------------------------------------------------------------------
!
      end module t_sph_field_FFTW
