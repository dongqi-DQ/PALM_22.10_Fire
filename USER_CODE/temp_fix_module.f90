!> @file user_module.f90
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
!> Declaration of user-defined variables. This module may only be used in the user-defined routines
!> (contained in user_interface.f90).
!--------------------------------------------------------------------------------------------------!
 MODULE temp_fix

    USE arrays_3d

    USE control_parameters

    USE cpulog

    USE indices

    USE kinds

    USE pegrid

    USE statistics

    USE surface_mod

    IMPLICIT NONE

    INTEGER(iwp) ::  dots_num_palm      !<
    INTEGER(iwp) ::  dots_num_user = 0  !<
    INTEGER(iwp) ::  user_idummy        !<

    LOGICAL ::  temp_fix_module_enabled = .TRUE.  !<

    REAL(wp) ::  user_rdummy  !<

    SAVE

    PRIVATE

!
!- Public functions
    PUBLIC                                                                                         &
       temp_fix_actions, temp_fix_parin 

!
!- Public parameters, constants and initial values
   PUBLIC                                                                                          &
      temp_fix_module_enabled

    INTERFACE temp_fix_actions
       MODULE PROCEDURE temp_fix_actions
    END INTERFACE temp_fix_actions

    INTERFACE temp_fix_parin
       MODULE PROCEDURE temp_fix_parin
    END INTERFACE temp_fix_parin   


 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE temp_fix_actions( location )


    CHARACTER(LEN=*) ::  location  !<

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<
    INTEGER(iwp) ::  n
    REAL(wp)     ::  pt_sum, pt_save
    REAL(wp) ::  fire_start_time

    CALL cpu_log( log_point(25), 'temp_fix_actions', 'start' )

!
!-- Here the user-defined actions follow. No calls for single grid points are allowed at locations
!-- before and after the timestep, since these calls are not within an i,j-loop
    SELECT CASE ( location )

       CASE ( 'before_timestep' )
!
!--       Enter actions to be done before every timestep here


       CASE ( 'before_prognostic_equations' )
!
!--       Enter actions to be done before all prognostic equations here

       CASE ( 'after_integration' )
!




fire_start_time = 7200.0_wp
        IF (coupling_char .EQ. '_N04') THEN !only correct in nest 
         IF (simulated_time .ge.fire_start_time)  THEN
          DO  i = nxlg, nxrg
              DO  j = nysg, nyn
              ! don't take top layers              
                 DO  k = nzb, nzt
! undershoot (less than IC) but not topography                 
                    IF ( (pt(k,j,i) <= t_min) .AND. (pt(k,j,i)>= -6000.0_wp) ) THEN
                            n = 0
                            pt_sum = 0.0
!                            pt_save=pt(k,j,i)
! pick cells without fire only ( <350) and with no undershoots                         
                            IF ((pt(k,j,i+1) <t_tran) .AND.(pt(k,j,i+1) >t_min))  THEN
                                    n = n+1
                                    pt_sum = pt_sum + pt(k,j,i+1)
                            ENDIF
                            IF ((pt(k,j,i-1) <t_tran) .AND.(pt(k,j,i-1) >t_min))  THEN
                                    n = n+1
                                    pt_sum = pt_sum + pt(k,j,i-1)
                            ENDIF
                            IF ((pt(k,j+1,i) <t_tran) .AND.(pt(k,j+1,i) >t_min))  THEN
                                    n = n+1
                                    pt_sum = pt_sum + pt(k,j+1,i)
                            ENDIF
                            IF ((pt(k,j-1,i) <t_tran).AND.(pt(k,j-1,i) >t_min))  THEN
                                    n = n+1
                                    pt_sum = pt_sum + pt(k,j-1,i)
                            ENDIF
                            ! if it's above fire
                            IF (pt(k-1,j,i) >=t_tran) THEN
                                    n = n+1
                                    pt_sum = pt_sum + pt(k-1,j,i) 
                            ENDIF

                            IF ((pt(k-1,j,i) >=t_tran) .AND. (pt(k+1,j,i) >t_min)) THEN
                                    n = n+1
                                    pt_sum = pt_sum + pt(k+1,j,i) 
                            ENDIF
                            ! if it's at the fire front
                            !! if both i-1 and i+1 are greater than Tfire_min
                            !!! check if they are overshoots
                            IF ((pt(k,j,i+1) >t_tran) .AND.(pt(k,j,i-1) >t_tran))  THEN
                                    IF (pt(k,j,i+1) < t_max) THEN
                                      n = n+1
                                      pt_sum = pt_sum + pt(k,j,i+1)
                                    ENDIF
                                    IF (pt(k,j,i-1) < t_max) THEN
                                      n = n+1
                                      pt_sum = pt_sum + pt(k,j,i-1)
                                    ENDIF

                            ENDIF
                            !! Do the same for j-1 and j+1
                            IF ((pt(k,j+1,i) >t_tran) .AND.(pt(k,j-1,i) >t_tran))  THEN
                                    IF (pt(k,j+1,i) < t_max) THEN
                                      n = n+1
                                      pt_sum = pt_sum + pt(k,j+1,i)
                                    ENDIF
                                    IF (pt(k,j-1,i) < t_max) THEN
                                      n = n+1
                                      pt_sum = pt_sum + pt(k,j-1,i)
                                    ENDIF

                            ENDIF

                            IF (n>0) THEN
                                    pt(k,j,i) =  pt_sum / n
                            ELSE
                                    pt(k,j,i) = 300
                            ENDIF


                            IF ( debug_output )  THEN
                                   CALL debug_message( 'user position 1', 'start' )
                                   WRITE( 9, * ) 'undershoot at ijk:', i,j,k,'old ',pt_save, &
                                   '4 values next to it: ',pt(k,j,i+1),pt(k,j,i-1),pt(k,j+1,i), &
                                   pt(k,j-1,i),  'new ',pt(k,j,i) ,'sum and n: ',pt_sum,n
                                   FLUSH ( 9)
                            ENDIF
                    ENDIF !end undershoot correction
! overshoots (more than 900)
                    IF  (pt(k,j,i) > t_max) THEN
                            n = 0
                            pt_sum = 0.0
!                             pt_save=pt(k,j,i)
! pick cells with fire only ( >500) and with no overshoots (<900)
                            IF ((pt(k,j,i+1) >t_tran).AND.(pt(k,j,i+1) .le.t_max))  THEN
                                    n = n+1
                                    pt_sum = pt_sum + pt(k,j,i+1)
                            ENDIF
                            IF ((pt(k,j,i-1) >t_tran).AND.(pt(k,j,i-1) .le.t_max)) THEN
                                    n = n+1
                                    pt_sum = pt_sum + pt(k,j,i-1)
                            ENDIF
                            IF ((pt(k,j+1,i) >t_tran).AND.(pt(k,j+1,i) .le.t_max)) THEN
                                    n = n+1
                                    pt_sum = pt_sum + pt(k,j+1,i)
                            ENDIF
                            IF ((pt(k,j-1,i) >t_tran).AND.(pt(k,j-1,i) .le.t_max)) THEN
                                    n = n+1
                                    pt_sum = pt_sum + pt(k,j-1,i)
                            ENDIF
                            IF (n>0) THEN
                                    pt(k,j,i) =  pt_sum / n
                            ELSE
                                    pt(k,j,i) = 900
                            ENDIF

                            IF ( debug_output )  THEN
                                   CALL debug_message( 'user position 1', 'start' )
                                   WRITE( 9, * ) 'overshoot at ijk:', i,j,k,'old ',pt_save, &
                                   '4 values next to it: ',pt(k,j,i+1),pt(k,j,i-1),pt(k,j+1,i), &
                                   pt(k,j-1,i),  'new ',pt(k,j,i) ,'sum and n: ',pt_sum,n, t_max
                                   FLUSH ( 9)
                            ENDIF


                    ENDIF ! end overshoot correction


                 ENDDO
              ENDDO
           ENDDO
         ENDIF !time condition 
        ENDIF  ! nest condition


       CASE ( 'after_timestep' )
!


       CASE ( 'u-tendency' )

       CASE ( 'v-tendency' )


       CASE ( 'w-tendency' )


       CASE ( 'pt-tendency' )


       CASE ( 'sa-tendency' )


       CASE ( 'e-tendency' )


       CASE ( 'q-tendency' )


       CASE ( 's-tendency' )


       CASE DEFAULT
          CONTINUE

    END SELECT

    CALL cpu_log( log_point(25), 'temp_fix_actions', 'stop' )

 END SUBROUTINE temp_fix_actions



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &user_parameters for user module
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE temp_fix_parin

    CHARACTER (LEN=80) ::  line  !< string containing the last line read from namelist file

    INTEGER(iwp) ::  i          !<
    INTEGER(iwp) ::  io_status  !< status after reading the namelist file
    INTEGER(iwp) ::  j          !<

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file

    NAMELIST /temp_fix_parameters/ switch_off_module, t_min, t_max, t_tran

!
!-- Next statement is to avoid compiler warnings about unused variables. Please remove in case
!-- that you are using them.
    IF ( dots_num_palm == 0  .OR.  dots_num_user == 0  .OR.  user_idummy == 0  .OR.                &
         user_rdummy == 0.0_wp )  CONTINUE

!
!-- Set revision number of this default interface version. It will be checked within the main
!-- program (palm). Please change the revision number in case that the current revision does not
!-- match with previous revisions (e.g. if routines have been added/deleted or if parameter lists
!-- in subroutines have been changed).
!    user_interface_current_revision = 'r4495'

!
!-- Position the namelist-file at the beginning (it has already been opened in parin), and try to
!-- read (find) a namelist named "user_parameters".
    REWIND ( 11 )
    READ( 11, temp_fix_parameters, IOSTAT=io_status )

!
!-- Actions depending on the READ status
    IF ( io_status == 0 )  THEN
!
!--    User namelist found and correctly read. Set default module switch to true. This activates
!--    calls of the user-interface subroutines.
       IF ( .NOT. switch_off_module )  temp_fix_module_enabled = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    User namelist was found, but contained errors. Print an error message containing the line
!--    that caused the problem
       BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'temp_fix', line )

    ENDIF

!
!-- Determine the number of user-defined profiles and append them to the standard data output
!-- (data_output_pr)
!    IF ( temp_fix_module_enabled )  THEN
!       IF ( data_output_pr_user(1) /= ' ' )  THEN
!          i = 1
!          DO WHILE ( data_output_pr(i) /= ' '  .AND.  i <= 200 )
!             i = i + 1
!          ENDDO
!          j = 1
!          DO WHILE ( data_output_pr_user(j) /= ' '  .AND.  j <= 300 )
!             data_output_pr(i) = data_output_pr_user(j)
!             max_pr_user_tmp   = max_pr_user_tmp + 1
!             i = i + 1
!             j = j + 1
!          ENDDO
!       ENDIF
!    ENDIF


 END SUBROUTINE temp_fix_parin


 END MODULE temp_fix
