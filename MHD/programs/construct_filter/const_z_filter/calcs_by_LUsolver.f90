!
!      module calcs_by_LUsolver
!
!      Written by H. Matsui
!
!      subroutine solve_z_commute_LU
!
!      subroutine solve_delta_z_LU
!      subroutine solve_delta_dz_LU
!      subroutine solve_delta_d2z_LU
!
      module calcs_by_LUsolver
!
      use m_precision
!
      use m_matrix_4_LU
      use m_ludcmp
!
      implicit none
!
      private :: solve_delta_z_etc_LU
!
!  ---------------------------------------------------------------------
!
      contains
!
!  ---------------------------------------------------------------------
!
      subroutine solve_z_commute_LU
!
      use m_geometry_parameter
      use m_crs_matrix
      use m_commute_filter_z
!
      integer(kind = kint) :: inod, i, j
!
       ncomp_lu = ncomp_mat
       call allocate_matrix_4_LU
!
       do inod = 1, numnod
!
         do i = 1, ncomp_mat
           do j = 1, ncomp_mat
             a_nod(j,i) = D_crs(j,i,inod)
           end do
           b_nod(i) = b_crs( ncomp_lu*(inod-1)+i )
         end do
!
!c decompose A = LU
         call ludcmp(a_nod,ncomp_lu,ncomp_lu,indx,d_nod)
!c solve Ax=LUx=b
         call lubksb(a_nod,ncomp_lu,ncomp_lu,indx,b_nod)
!
         do i = 1, ncomp_mat
           X_crs( ncomp_mat*(inod-1)+i ) = b_nod(i)
           d_nod = d_nod*a_nod(i,i)
         end do
!
         write(*,*) 'det A', inod, d_nod
!
       end do
!
      end subroutine solve_z_commute_LU
!
!  ---------------------------------------------------------------------
!
      subroutine solve_delta_z_LU
!
      use m_int_edge_vart_width
!
      call solve_delta_z_etc_LU(delta_z)
!
      end subroutine solve_delta_z_LU
!
!  ---------------------------------------------------------------------
!
      subroutine solve_delta_dz_LU
!
      use m_int_edge_vart_width
!
      call solve_delta_z_etc_LU(delta_dz)
!
      end subroutine solve_delta_dz_LU
!
!  ---------------------------------------------------------------------
!
      subroutine solve_delta_d2z_LU
!
      use m_int_edge_vart_width
!
      call solve_delta_z_etc_LU(d2_dz)
!
      end subroutine solve_delta_d2z_LU
!
!  ---------------------------------------------------------------------
!
      subroutine solve_delta_z_etc_LU(X_lu)
!
      use m_geometry_parameter
      use m_int_edge_data
      use m_int_edge_vart_width
!
      real(kind = kreal), intent(inout) :: X_lu(numnod)
!
      integer(kind = kint) :: inod, i, j
!
         ncomp_lu = numnod
!
         call allocate_matrix_4_LU
!
         do i = 1, numnod
           do j = 1, numnod
             a_nod(j,i) = mk_c(i,j)
           end do
           b_nod(i) = rhs_dz(i)
         end do
!
!c decompose A = LU
         call ludcmp(a_nod,ncomp_lu,ncomp_lu,indx,d_nod)
!c solve Ax=LUx=b
         call lubksb(a_nod,ncomp_lu,ncomp_lu,indx,b_nod)
!
         do i = 1, numnod
           X_lu(i) = b_nod(i)
           d_nod = d_nod*a_nod(i,i)
         end do
!
         write(*,*) 'det A', d_nod
!
         call deallocate_matrix_4_LU
!
      end subroutine solve_delta_z_etc_LU
!
!  ---------------------------------------------------------------------
!
      end module calcs_by_LUsolver