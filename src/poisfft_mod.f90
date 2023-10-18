!> @file poisfft_mod.f90
!--------------------------------------------------------------------------------------------------!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the terms of the GNU General
! Public License as published by the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
! implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
! Public License for more details.
!
! You should have received a copy of the GNU General Public License along with PALM. If not, see
! <http://www.gnu.org/licenses/>.
!
! Copyright 1997-2021 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
!
! Description:
! ------------
!> Solves the Poisson equation with a 2D spectral method
!>        d^2 p / dx^2 + d^2 p / dy^2 + d^2 p / dz^2 = s
!>
!> Input:
!> real   ar   contains (nnz,nny,nnx) elements of the velocity divergence, starting from (1,nys,nxl)
!>
!> Output:
!> real   ar   contains the solution for perturbation pressure p
!--------------------------------------------------------------------------------------------------!
 MODULE poisfft_mod

#if defined( __parallel )
    USE MPI
#endif

    USE control_parameters,                                                                        &
        ONLY:  bc_lr_cyc,                                                                          &
               bc_ns_cyc,                                                                          &
               message_string,                                                                     &
               temperton_fft_vec

    USE fft_xy,                                                                                    &
        ONLY:  fft_init,                                                                           &
               fft_y,                                                                              &
               fft_y_1d,                                                                           &
               fft_y_m,                                                                            &
               fft_x,                                                                              &
               fft_x_1d,                                                                           &
               fft_x_m

    USE indices,                                                                                   &
        ONLY:  nnx,                                                                                &
               nny,                                                                                &
               nx,                                                                                 &
               nxl,                                                                                &
               nxr,                                                                                &
               ny,                                                                                 &
               nys,                                                                                &
               nyn,                                                                                &
               nz

    USE pegrid,                                                                                    &
        ONLY:  non_uniform_subdomain

    USE transpose_mod,                                                                             &
        ONLY:  nnx_x_max,                                                                          &
               nnx_y_max,                                                                          &
               nnz_x_max,                                                                          &
               nnz_z_max,                                                                          &
               nxl_y,                                                                              &
               nxl_z,                                                                              &
               nxr_x_max,                                                                          &
               nxr_y_max,                                                                          &
               nxr_y,                                                                              &
               nxr_z,                                                                              &
               nx_y_max,                                                                           &
               nys_x,                                                                              &
               nys_z,                                                                              &
               nyn_x,                                                                              &
               nyn_x_max,                                                                          &
               nyn_z,                                                                              &
               nyn_z_max,                                                                          &
               ny_z_max,                                                                           &
               nzb_x,                                                                              &
               nzb_y,                                                                              &
               nzt_x,                                                                              &
               nzt_x_max,                                                                          &
               nzt_y_max,                                                                          &
               nzt_y,                                                                              &
               nz_x_max,                                                                           &
               resort_for_xy, transpose_xy,                                                        &
               resort_for_xz, transpose_xz,                                                        &
               resort_for_yx, transpose_yx,                                                        &
               resort_for_yz, transpose_yz,                                                        &
               resort_for_zx, transpose_zx,                                                        &
               resort_for_zy, transpose_zy

    USE tridia_solver,                                                                             &
        ONLY:  tridia_1dd,                                                                         &
               tridia_init,                                                                        &
               tridia_substi


    IMPLICIT NONE

    LOGICAL, SAVE ::  poisfft_initialized = .FALSE.  !<

    PRIVATE

    PUBLIC  poisfft, poisfft_init

    INTERFACE poisfft
       MODULE PROCEDURE poisfft
    END INTERFACE poisfft

    INTERFACE poisfft_init
       MODULE PROCEDURE poisfft_init
    END INTERFACE poisfft_init


 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Setup coefficients for FFT and the tridiagonal solver
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE poisfft_init

    IMPLICIT NONE


    CALL fft_init

    CALL tridia_init

    poisfft_initialized = .TRUE.

 END SUBROUTINE poisfft_init



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Two-dimensional Fourier Transformation in x- and y-direction.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE poisfft( ar )

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE indices,                                                                                   &
        ONLY:                                                                                      &
               nxl,                                                                                &
               nxr,                                                                                &
               nyn,                                                                                &
               nys,                                                                                &
               nz

    USE kinds

    USE pegrid

    IMPLICIT NONE

    REAL(wp), DIMENSION(1:nz,nys:nyn,nxl:nxr) ::  ar      !<

#define __acc_fft_device ( defined( _OPENACC ) && ( defined ( __cuda_fft ) ) )
#if __acc_fft_device
    !$ACC DECLARE CREATE(ar_inv)
#endif

    REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  xy_in    !<
    REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  xy_out   !<
    REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  yz_in    !<
    REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  yz_out   !<
    REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  zx_in    !<
    REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  zx_out   !<


    IF ( .NOT. poisfft_initialized )  CALL poisfft_init

    CALL cpu_log( log_point_s(3), 'poisfft', 'start' )

#if !__acc_fft_device
    !$ACC UPDATE HOST(ar)
#endif

#ifndef _OPENACC
!
!-- Two-dimensional Fourier Transformation in x- and y-direction.
    IF ( npey == 1  .AND.  npex > 1 )  THEN

!
!--    1d-domain-decomposition along x:
!--    FFT along y and transposition y --> x
       CALL ffty_tr_yx( ar, ar )

!
!--    FFT along x, solving the tridiagonal system and backward FFT
       CALL fftx_tri_fftx( ar )

!
!--    Transposition x --> y and backward FFT along y
       CALL tr_xy_ffty( ar, ar )

    ELSEIF ( npex == 1 .AND. npey > 1 )  THEN

!
!--    1d-domain-decomposition along y:
!--    FFT along x and transposition x --> y
       CALL fftx_tr_xy( ar, ar )

!
!--    FFT along y, solving the tridiagonal system and backward FFT
       CALL ffty_tri_ffty( ar )

!
!--    Transposition y --> x and backward FFT along x
       CALL tr_yx_fftx( ar, ar )

    ELSE
#endif

       ALLOCATE( xy_in(nys_x:nyn_x_max,nzb_x:nzt_x,0:nx_y_max) )
       ALLOCATE( zx_in(nys:nyn,nxl:nxr_x_max,1:nz_x_max) )
       ALLOCATE( zx_out(0:nx,nys_x:nyn_x,nzb_x:nzt_x) )
!
!--    2d-domain-decomposition or no decomposition (1 PE run)
!--    Transposition z --> x
       CALL cpu_log( log_point_s(5), 'transpo forward', 'start' )
       CALL resort_for_zx( ar, zx_in )
!
!--    In case of temperton_fft_vec, zx_out is bypassed by f_vec_x.
       CALL transpose_zx( zx_in, zx_out )
       CALL cpu_log( log_point_s(5), 'transpo forward', 'pause' )

       CALL cpu_log( log_point_s(4), 'fft_x', 'start' )
       IF ( temperton_fft_vec  .AND. .NOT.  non_uniform_subdomain)  THEN
!
!--       Vector version outputs a transformed array ar_inv that does not require resorting
!--       (which is done for ar further below)
          CALL fft_x( zx_out, 'forward',  ar_inv=xy_in)
       ELSE
          CALL fft_x( zx_out, 'forward')
       ENDIF
       CALL cpu_log( log_point_s(4), 'fft_x', 'pause' )

!
!--    Transposition x --> y
       ALLOCATE( xy_out(0:ny,nxl_y:nxr_y,nzb_y:nzt_y) )
       ALLOCATE( yz_in(nxl_y:nxr_y,nzb_y:nzt_y_max,0:ny_z_max) )

       CALL cpu_log( log_point_s(5), 'transpo forward', 'continue' )
       IF ( .NOT. temperton_fft_vec  .OR.  non_uniform_subdomain )  THEN
          CALL resort_for_xy( zx_out, xy_in )
       ENDIF
       CALL transpose_xy( xy_in, xy_out )
       CALL cpu_log( log_point_s(5), 'transpo forward', 'pause' )

       CALL cpu_log( log_point_s(7), 'fft_y', 'start' )
       IF ( temperton_fft_vec .AND. .NOT.  non_uniform_subdomain)  THEN
!
!--       Input array ar_inv from fft_x can be directly used here.
!--       The output (also in array ar_inv) does not require resorting below.

!--       This is currently only programmed for uniform sundomains.
!--       TODO: Please check performance on NEC to decide if this branch should also be
!--             implemented for nonuniform subdomains.
!--       This branch saves one resort call, the vector version of Temperton-fft is active in
!--       both cases of this IF condition.
          CALL fft_y( xy_out, 'forward', ar_inv = yz_in, nxl_y_bound = nxl_y, nxr_y_bound = nxr_y, &
                      nxl_y_l = nxl_y, nxr_y_l = nxr_y )
       ELSE
          CALL fft_y( xy_out, 'forward', ar_tr = xy_out, nxl_y_bound = nxl_y, nxr_y_bound = nxr_y, &
                      nxl_y_l = nxl_y, nxr_y_l = nxr_y )
       ENDIF
       CALL cpu_log( log_point_s(7), 'fft_y', 'pause' )

!
!--    Transposition y --> z
       ALLOCATE( yz_out(nxl_z:nxr_z,nys_z:nyn_z,1:nz) )

       CALL cpu_log( log_point_s(5), 'transpo forward', 'continue' )
       IF ( .NOT. temperton_fft_vec  .OR.  non_uniform_subdomain )  THEN
          CALL resort_for_yz( xy_out, yz_in )
       ENDIF
       CALL transpose_yz( yz_in, yz_out )
       CALL cpu_log( log_point_s(5), 'transpo forward', 'stop' )

!
!--    Solve the tridiagonal equation system along z
       CALL cpu_log( log_point_s(6), 'tridia', 'start' )
       CALL tridia_substi( yz_out )
       CALL cpu_log( log_point_s(6), 'tridia', 'stop' )
!
!
!--    Inverse Fourier Transformation
!--    Transposition z --> y
       CALL cpu_log( log_point_s(8), 'transpo invers', 'start' )
       CALL transpose_zy( yz_out, yz_in )
!
!--    The fft_y below (vector branch) can directly process ar_inv (i.e. does not require a
!--    resorting)
       IF ( .NOT. temperton_fft_vec  .OR. non_uniform_subdomain )  THEN
          CALL resort_for_zy( yz_in, xy_out )
       ENDIF
       CALL cpu_log( log_point_s(8), 'transpo invers', 'pause' )

       DEALLOCATE( yz_out )

       CALL cpu_log( log_point_s(7), 'fft_y', 'continue' )
       IF ( temperton_fft_vec  .AND.  .NOT. non_uniform_subdomain )  THEN
!
!--       Output array ar_inv can be used as input to the below fft_x routine without resorting
          CALL fft_y( xy_out, 'backward', ar_inv = yz_in, nxl_y_bound = nxl_y, nxr_y_bound = nxr_y,&
                      nxl_y_l = nxl_y, nxr_y_l = nxr_y )
       ELSE
          CALL fft_y( xy_out, 'backward', ar_tr = xy_out, nxl_y_bound = nxl_y, nxr_y_bound = nxr_y,&
                      nxl_y_l = nxl_y, nxr_y_l = nxr_y )
       ENDIF

       DEALLOCATE( yz_in )

       CALL cpu_log( log_point_s(7), 'fft_y', 'stop' )

!
!--    Transposition y --> x
       CALL cpu_log( log_point_s(8), 'transpo invers', 'continue' )
       CALL transpose_yx( xy_out, xy_in )
       IF ( .NOT. temperton_fft_vec  .OR.  non_uniform_subdomain )  THEN
          CALL resort_for_yx( xy_in, zx_out )
       ENDIF
       CALL cpu_log( log_point_s(8), 'transpo invers', 'pause' )

       DEALLOCATE( xy_out )

       CALL cpu_log( log_point_s(4), 'fft_x', 'continue' )
       IF ( temperton_fft_vec  .AND.  .NOT. non_uniform_subdomain )  THEN
          CALL fft_x( zx_out, 'backward',  ar_inv=xy_in )
       ELSE
          CALL fft_x( zx_out, 'backward' )
       ENDIF

       DEALLOCATE( xy_in )

       CALL cpu_log( log_point_s(4), 'fft_x', 'stop' )

!
!--    Transposition x --> z
       CALL cpu_log( log_point_s(8), 'transpo invers', 'continue' )
       CALL transpose_xz( zx_out, zx_in )
       CALL resort_for_xz( zx_in, ar )
       CALL cpu_log( log_point_s(8), 'transpo invers', 'stop' )

       DEALLOCATE( zx_in )
       DEALLOCATE( zx_out )

#ifndef _OPENACC
    ENDIF
#endif

#if !__acc_fft_device
    !$ACC UPDATE DEVICE(ar)
#endif

    CALL cpu_log( log_point_s(3), 'poisfft', 'stop' )

 END SUBROUTINE poisfft


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Fourier-transformation along y with subsequent transposition y --> x for a 1d-decomposition
!> along x.
!>
!> @attention The performance of this routine is much faster on the NEC-SX6, if the first index of
!>            work_ffty_vec is odd. Otherwise memory bank conflicts may occur (especially if the
!>            index is a multiple of 128). That's why work_ffty_vec is dimensioned as 0:ny+1.
!>            Of course, this will not work if users are using an odd number of gridpoints along y.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE ffty_tr_yx( f_in, f_out )

    USE control_parameters,                                                                        &
        ONLY:  loop_optimization

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE kinds

    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  i       !<
    INTEGER(iwp) ::  iend    !<
    INTEGER(iwp) ::  iouter  !<
    INTEGER(iwp) ::  ir      !<
    INTEGER(iwp) ::  j       !<
    INTEGER(iwp) ::  k       !<

    INTEGER(iwp), PARAMETER ::  stridex = 4  !<

    REAL(wp), DIMENSION(1:nz,0:ny,nxl:nxr)         ::  f_in   !<
    REAL(wp), DIMENSION(nnx,1:nz,nys_x:nyn_x,npex) ::  f_out  !<
    REAL(wp), DIMENSION(nxl:nxr,1:nz,0:ny)         ::  work   !<

    REAL(wp), DIMENSION(:,:), ALLOCATABLE   ::  work_ffty      !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  work_ffty_vec  !<

!
!-- Carry out the FFT along y, where all data are present due to the 1d-decomposition along x.
!-- Resort the data in a way that x becomes the first index.
    CALL cpu_log( log_point_s(7), 'fft_y_1d', 'start' )

    IF ( loop_optimization == 'vector' )  THEN

       ALLOCATE( work_ffty_vec(0:ny+1,1:nz,nxl:nxr) )
!
!--    Code optimized for vector processors
       !$OMP PARALLEL PRIVATE ( i, j, k )
       !$OMP DO
       DO  i = nxl, nxr

          DO  j = 0, ny
             DO  k = 1, nz
                work_ffty_vec(j,k,i) = f_in(k,j,i)
             ENDDO
          ENDDO

          CALL fft_y_m( work_ffty_vec(:,:,i), ny+1, 'forward' )

       ENDDO

       !$OMP DO
       DO  k = 1, nz
          DO  j = 0, ny
             DO  i = nxl, nxr
                work(i,k,j) = work_ffty_vec(j,k,i)
             ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL

       DEALLOCATE( work_ffty_vec )

    ELSE
!
!--    Cache optimized code.
       ALLOCATE( work_ffty(0:ny,stridex) )
!
!--    The i-(x-)direction is split into a strided outer loop and an inner loop for better cache
!--    performance
       !$OMP PARALLEL PRIVATE (i,iend,iouter,ir,j,k,work_ffty)
       !$OMP DO
       DO  iouter = nxl, nxr, stridex

          iend = MIN( iouter+stridex-1, nxr )  ! Upper bound for inner i loop

          DO  k = 1, nz

             DO  i = iouter, iend

                ir = i-iouter+1  ! Counter within a stride
                DO  j = 0, ny
                   work_ffty(j,ir) = f_in(k,j,i)
                ENDDO
!
!--             FFT along y
                CALL fft_y_1d( work_ffty(:,ir), 'forward' )

             ENDDO

!
!--          Resort
             DO  j = 0, ny
                DO  i = iouter, iend
                   work(i,k,j) = work_ffty(j,i-iouter+1)
                ENDDO
             ENDDO

          ENDDO

       ENDDO
       !$OMP END PARALLEL

       DEALLOCATE( work_ffty )

    ENDIF

    CALL cpu_log( log_point_s(7), 'fft_y_1d', 'pause' )

!
!-- Transpose array
#if defined( __parallel )
    CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'start' )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLTOALL( work(nxl,1,0), sendrecvcount_xy, MPI_REAL, f_out(1,1,nys_x,1),              &
                       sendrecvcount_xy, MPI_REAL, comm1dx, ierr )
    CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'stop' )
#else
!
!-- Next line required to avoid compile error about unused dummy argument in serial mode
    i = SIZE( f_out )
#endif

 END SUBROUTINE ffty_tr_yx


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Transposition x --> y with a subsequent backward Fourier transformation for a 1d-decomposition
!> along x
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE tr_xy_ffty( f_in, f_out )

    USE control_parameters,                                                                        &
        ONLY:  loop_optimization

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE kinds

    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  i       !<
    INTEGER(iwp) ::  iend    !<
    INTEGER(iwp) ::  iouter  !<
    INTEGER(iwp) ::  ir      !<
    INTEGER(iwp) ::  j       !<
    INTEGER(iwp) ::  k       !<

    INTEGER(iwp), PARAMETER ::  stridex = 4  !<

    REAL(wp), DIMENSION(nnx,1:nz,nys_x:nyn_x,npex) ::  f_in   !<
    REAL(wp), DIMENSION(1:nz,0:ny,nxl:nxr)         ::  f_out  !<
    REAL(wp), DIMENSION(nxl:nxr,1:nz,0:ny)         ::  work   !<

    REAL(wp), DIMENSION(:,:),   ALLOCATABLE ::  work_ffty      !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  work_ffty_vec  !<

!
!-- Transpose array
#if defined( __parallel )
    CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'start' )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLTOALL( f_in(1,1,nys_x,1), sendrecvcount_xy, MPI_REAL, work(nxl,1,0),               &
                       sendrecvcount_xy, MPI_REAL, comm1dx, ierr )
    CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'stop' )
#else
!
!-- Next line required to avoid compile error about unused dummy argument in serial mode
    i = SIZE( f_in )
#endif

!
!-- Resort the data in a way that y becomes the first index and carry out the backward fft along y.
    CALL cpu_log( log_point_s(7), 'fft_y_1d', 'continue' )

    IF ( loop_optimization == 'vector' )  THEN

       ALLOCATE( work_ffty_vec(0:ny+1,1:nz,nxl:nxr) )
!
!--    Code optimized for vector processors
       !$OMP PARALLEL PRIVATE ( i, j, k )
       !$OMP DO
       DO  k = 1, nz
          DO  j = 0, ny
             DO  i = nxl, nxr
                work_ffty_vec(j,k,i) = work(i,k,j)
             ENDDO
          ENDDO
       ENDDO

       !$OMP DO
       DO  i = nxl, nxr

          CALL fft_y_m( work_ffty_vec(:,:,i), ny+1, 'backward' )

          DO  j = 0, ny
             DO  k = 1, nz
                f_out(k,j,i) = work_ffty_vec(j,k,i)
             ENDDO
          ENDDO

       ENDDO
       !$OMP END PARALLEL

       DEALLOCATE( work_ffty_vec )

    ELSE
!
!--    Cache optimized code.
       ALLOCATE( work_ffty(0:ny,stridex) )
!
!--    The i-(x-)direction is split into a strided outer loop and an inner loop for better cache
!--    performance
       !$OMP PARALLEL PRIVATE ( i, iend, iouter, ir, j, k, work_ffty )
       !$OMP DO
       DO  iouter = nxl, nxr, stridex

          iend = MIN( iouter+stridex-1, nxr )  ! Upper bound for inner i loop

          DO  k = 1, nz
!
!--          Resort
             DO  j = 0, ny
                DO  i = iouter, iend
                   work_ffty(j,i-iouter+1) = work(i,k,j)
                ENDDO
             ENDDO

             DO  i = iouter, iend

!
!--             FFT along y
                ir = i-iouter+1  ! Counter within a stride
                CALL fft_y_1d( work_ffty(:,ir), 'backward' )

                DO  j = 0, ny
                   f_out(k,j,i) = work_ffty(j,ir)
                ENDDO
             ENDDO

          ENDDO

       ENDDO
       !$OMP END PARALLEL

       DEALLOCATE( work_ffty )

    ENDIF

    CALL cpu_log( log_point_s(7), 'fft_y_1d', 'stop' )

 END SUBROUTINE tr_xy_ffty


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> FFT along x, solution of the tridiagonal system and backward FFT for a 1d-decomposition along x
!>
!> @warning This subroutine may still not work for hybrid parallelization with OpenMP (for possible
!>          necessary changes see the original routine poisfft_hybrid, developed by Klaus Ketelsen,
!>          May 2002)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE fftx_tri_fftx( ar )

    USE control_parameters,                                                                        &
        ONLY:  loop_optimization

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE grid_variables,                                                                            &
        ONLY:  ddx2,                                                                               &
               ddy2

    USE kinds

    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  i                   !<
    INTEGER(iwp) ::  j                   !<
    INTEGER(iwp) ::  k                   !<
    INTEGER(iwp) ::  m                   !<
    INTEGER(iwp) ::  n                   !<
!$  INTEGER(iwp) ::  omp_get_thread_num  !<
    INTEGER(iwp) ::  tn                  !<

    REAL(wp), DIMENSION(0:nx)                      ::  work_fftx  !<
    REAL(wp), DIMENSION(0:nx,1:nz)                 ::  work_trix  !<
    REAL(wp), DIMENSION(nnx,1:nz,nys_x:nyn_x,npex) ::  ar         !<
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE      ::  tri        !<


    CALL cpu_log( log_point_s(33), 'fft_x_1d + tridia', 'start' )

    ALLOCATE( tri(5,0:nx,0:nz-1,0:threads_per_task-1) )

    tn = 0              ! Default thread number in case of one thread
    !$OMP  PARALLEL DO PRIVATE ( i, j, k, m, n, tn, work_fftx, work_trix )
    DO  j = nys_x, nyn_x

!$     tn = omp_get_thread_num()

       IF ( loop_optimization == 'vector' )  THEN
!
!--       Code optimized for vector processors
          DO  k = 1, nz

             m = 0
             DO  n = 1, npex
                DO  i = 1, nnx
                   work_trix(m,k) = ar(i,k,j,n)
                   m = m + 1
                ENDDO
             ENDDO

          ENDDO

          CALL fft_x_m( work_trix, 'forward' )

       ELSE
!
!--       Cache optimized code
          DO  k = 1, nz

             m = 0
             DO  n = 1, npex
                DO  i = 1, nnx
                   work_fftx(m) = ar(i,k,j,n)
                   m = m + 1
                ENDDO
             ENDDO

             CALL fft_x_1d( work_fftx, 'forward' )

             DO  i = 0, nx
                work_trix(i,k) = work_fftx(i)
             ENDDO

          ENDDO

       ENDIF

!
!--    Solve the linear equation system
       CALL tridia_1dd( ddx2, ddy2, nx, ny, j, work_trix, tri(:,:,:,tn) )

       IF ( loop_optimization == 'vector' )  THEN
!
!--       Code optimized for vector processors
          CALL fft_x_m( work_trix, 'backward' )

          DO  k = 1, nz

             m = 0
             DO  n = 1, npex
                DO  i = 1, nnx
                   ar(i,k,j,n) = work_trix(m,k)
                   m = m + 1
                ENDDO
             ENDDO

          ENDDO

       ELSE
!
!--       Cache optimized code
          DO  k = 1, nz

             DO  i = 0, nx
                work_fftx(i) = work_trix(i,k)
             ENDDO

             CALL fft_x_1d( work_fftx, 'backward' )

             m = 0
             DO  n = 1, npex
                DO  i = 1, nnx
                   ar(i,k,j,n) = work_fftx(m)
                   m = m + 1
                ENDDO
             ENDDO

          ENDDO

       ENDIF

    ENDDO

    DEALLOCATE( tri )

    CALL cpu_log( log_point_s(33), 'fft_x_1d + tridia', 'stop' )

 END SUBROUTINE fftx_tri_fftx


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Fourier-transformation along x with subsequent transposition x --> y for a 1d-decomposition
!> along y.
!>
!> @attention NEC-branch of this routine may significantly profit from further optimizations. So
!>            far, performance is much worse than for routine ffty_tr_yx (more than three times
!>            slower).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE fftx_tr_xy( f_in, f_out )


    USE control_parameters,                                                                        &
        ONLY:  loop_optimization

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE kinds

    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<

    REAL(wp), DIMENSION(0:nx,1:nz,nys:nyn)         ::  work_fftx  !<
    REAL(wp), DIMENSION(1:nz,nys:nyn,0:nx)         ::  f_in       !<
    REAL(wp), DIMENSION(nny,1:nz,nxl_y:nxr_y,npey) ::  f_out      !<
    REAL(wp), DIMENSION(nys:nyn,1:nz,0:nx)         ::  work       !<

!
!--  Carry out the FFT along x, where all data are present due to the 1d-decomposition along y.
!--  Resort the data in a way that y becomes the first index.
     CALL cpu_log( log_point_s(4), 'fft_x_1d', 'start' )

     IF ( loop_optimization == 'vector' )  THEN
!
!--    Code for vector processors
       !$OMP  PARALLEL PRIVATE ( i, j, k )
       !$OMP  DO
       DO  i = 0, nx

          DO  j = nys, nyn
             DO  k = 1, nz
                work_fftx(i,k,j) = f_in(k,j,i)
             ENDDO
          ENDDO

       ENDDO

       !$OMP  DO
       DO  j = nys, nyn

          CALL fft_x_m( work_fftx(:,:,j), 'forward' )

          DO  k = 1, nz
             DO  i = 0, nx
                work(j,k,i) = work_fftx(i,k,j)
             ENDDO
          ENDDO

       ENDDO
       !$OMP  END PARALLEL

    ELSE

!
!--    Cache optimized code (there might still be a potential for better optimization).
       !$OMP  PARALLEL PRIVATE (i,j,k)
       !$OMP  DO
       DO  i = 0, nx

          DO  j = nys, nyn
             DO  k = 1, nz
                work_fftx(i,k,j) = f_in(k,j,i)
             ENDDO
          ENDDO

       ENDDO

       !$OMP  DO
       DO  j = nys, nyn
          DO  k = 1, nz

             CALL fft_x_1d( work_fftx(0:nx,k,j), 'forward' )

             DO  i = 0, nx
                work(j,k,i) = work_fftx(i,k,j)
             ENDDO
          ENDDO

       ENDDO
       !$OMP  END PARALLEL

    ENDIF
    CALL cpu_log( log_point_s(4), 'fft_x_1d', 'pause' )

!
!-- Transpose array
#if defined( __parallel )
    CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'start' )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLTOALL( work(nys,1,0), sendrecvcount_xy, MPI_REAL, f_out(1,1,nxl_y,1),              &
                       sendrecvcount_xy, MPI_REAL, comm1dy, ierr )
    CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'stop' )
#else
!
!-- Next line required to avoid compile error about unused dummy argument in serial mode
    i = SIZE( f_out )
#endif

 END SUBROUTINE fftx_tr_xy


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Transposition y --> x with a subsequent backward Fourier transformation for a 1d-decomposition
!> along x.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE tr_yx_fftx( f_in, f_out )


    USE control_parameters,                                                                        &
        ONLY:  loop_optimization

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE kinds

    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<

    REAL(wp), DIMENSION(nny,1:nz,nxl_y:nxr_y,npey) ::  f_in       !<
    REAL(wp), DIMENSION(1:nz,nys:nyn,0:nx)         ::  f_out      !<
    REAL(wp), DIMENSION(nys:nyn,1:nz,0:nx)         ::  work       !<
    REAL(wp), DIMENSION(0:nx,1:nz,nys:nyn)         ::  work_fftx  !<


!
!-- Transpose array
#if defined( __parallel )
    CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'start' )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLTOALL( f_in(1,1,nxl_y,1), sendrecvcount_xy, MPI_REAL, work(nys,1,0),               &
                       sendrecvcount_xy, MPI_REAL, comm1dy, ierr )
    CALL cpu_log( log_point_s(32), 'mpi_alltoall', 'stop' )
#else
!
!-- Next line required to avoid compile error about unused dummy argument in serial mode
    i = SIZE( f_in )
#endif

!
!-- Carry out the FFT along x, where all data are present due to the 1d-decomposition along y.
!-- Resort the data in a way that y becomes the first index.
    CALL cpu_log( log_point_s(4), 'fft_x_1d', 'continue' )

    IF ( loop_optimization == 'vector' )  THEN
!
!--    Code optimized for vector processors
!$OMP  PARALLEL PRIVATE ( i, j, k )
!$OMP  DO
       DO  j = nys, nyn

          DO  k = 1, nz
             DO  i = 0, nx
                work_fftx(i,k,j) = work(j,k,i)
             ENDDO
          ENDDO

          CALL fft_x_m( work_fftx(:,:,j), 'backward' )

       ENDDO

!$OMP  DO
       DO  i = 0, nx
          DO  j = nys, nyn
             DO  k = 1, nz
                f_out(k,j,i) = work_fftx(i,k,j)
             ENDDO
          ENDDO
       ENDDO
!$OMP  END PARALLEL

    ELSE

!
!--    Cache optimized code (there might be still a potential for better optimization).
!$OMP  PARALLEL PRIVATE (i,j,k)
!$OMP  DO
       DO  j = nys, nyn
          DO  k = 1, nz

             DO  i = 0, nx
                work_fftx(i,k,j) = work(j,k,i)
             ENDDO

             CALL fft_x_1d( work_fftx(0:nx,k,j), 'backward' )

          ENDDO
       ENDDO

!$OMP  DO
       DO  i = 0, nx
          DO  j = nys, nyn
             DO  k = 1, nz
                f_out(k,j,i) = work_fftx(i,k,j)
             ENDDO
          ENDDO
       ENDDO
!$OMP  END PARALLEL

    ENDIF
    CALL cpu_log( log_point_s(4), 'fft_x_1d', 'stop' )

 END SUBROUTINE tr_yx_fftx


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> FFT along y, solution of the tridiagonal system and backward FFT for a 1d-decomposition along y.
!>
!> @warning This subroutine may still not work for hybrid parallelization with OpenMP (for possible
!>          necessary changes see the original routine poisfft_hybrid, developed by Klaus Ketelsen,
!>          May 2002)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE ffty_tri_ffty( ar )


    USE control_parameters,                                                                        &
        ONLY:  loop_optimization

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE grid_variables,                                                                            &
        ONLY:  ddx2,                                                                               &
               ddy2

    USE kinds

    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  i                   !<
    INTEGER(iwp) ::  j                   !<
    INTEGER(iwp) ::  k                   !<
    INTEGER(iwp) ::  m                   !<
    INTEGER(iwp) ::  n                   !<
!$  INTEGER(iwp) ::  omp_get_thread_num  !<
    INTEGER(iwp) ::  tn                  !<

    REAL(wp), DIMENSION(0:ny)                      ::  work_ffty  !<
    REAL(wp), DIMENSION(0:ny,1:nz)                 ::  work_triy  !<
    REAL(wp), DIMENSION(nny,1:nz,nxl_y:nxr_y,npey) ::  ar         !<
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE      ::  tri        !<


    CALL cpu_log( log_point_s(39), 'fft_y_1d + tridia', 'start' )

    ALLOCATE( tri(5,0:ny,0:nz-1,0:threads_per_task-1) )

    tn = 0           ! Default thread number in case of one thread
!$OMP PARALLEL DO PRIVATE ( i, j, k, m, n, tn, work_ffty, work_triy )
    DO  i = nxl_y, nxr_y

!$     tn = omp_get_thread_num()

       IF ( loop_optimization == 'vector' )  THEN
!
!--       Code optimized for vector processors
          DO  k = 1, nz

             m = 0
             DO  n = 1, npey
                DO  j = 1, nny
                   work_triy(m,k) = ar(j,k,i,n)
                   m = m + 1
                ENDDO
             ENDDO

          ENDDO

          CALL fft_y_m( work_triy, ny, 'forward' )

       ELSE
!
!--       Cache optimized code
          DO  k = 1, nz

             m = 0
             DO  n = 1, npey
                DO  j = 1, nny
                   work_ffty(m) = ar(j,k,i,n)
                   m = m + 1
                ENDDO
             ENDDO

             CALL fft_y_1d( work_ffty, 'forward' )

             DO  j = 0, ny
                work_triy(j,k) = work_ffty(j)
             ENDDO

          ENDDO

       ENDIF

!
!--    Solve the linear equation system
       CALL tridia_1dd( ddy2, ddx2, ny, nx, i, work_triy, tri(:,:,:,tn) )

       IF ( loop_optimization == 'vector' )  THEN
!
!--       Code optimized for vector processors
          CALL fft_y_m( work_triy, ny, 'backward' )

          DO  k = 1, nz

             m = 0
             DO  n = 1, npey
                DO  j = 1, nny
                   ar(j,k,i,n) = work_triy(m,k)
                   m = m + 1
                ENDDO
             ENDDO

          ENDDO

       ELSE
!
!--       Cache optimized code
          DO  k = 1, nz

             DO  j = 0, ny
                work_ffty(j) = work_triy(j,k)
             ENDDO

             CALL fft_y_1d( work_ffty, 'backward' )

             m = 0
             DO  n = 1, npey
                DO  j = 1, nny
                   ar(j,k,i,n) = work_ffty(m)
                   m = m + 1
                ENDDO
             ENDDO

          ENDDO

       ENDIF

    ENDDO

    DEALLOCATE( tri )

    CALL cpu_log( log_point_s(39), 'fft_y_1d + tridia', 'stop' )

 END SUBROUTINE ffty_tri_ffty

 END MODULE poisfft_mod
