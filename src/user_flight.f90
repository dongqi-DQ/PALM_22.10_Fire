!> @file user_flight.f90
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
!> Calculation of user-defined output quantity for flight measurements after each timestep.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_flight( var, id )

    USE control_parameters

    USE grid_variables

    USE indices

    USE kinds

    USE user

    USE arrays_3d

    IMPLICIT NONE

!    INTEGER(iwp) ::  i  !< index along x
!    INTEGER(iwp) ::  j  !< index along y
!    INTEGER(iwp) ::  k  !< index along z
    INTEGER(iwp) ::  id  !< variable identifyer, according to the settings in user_init_flight

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  var  !< treated variable

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( id == 0  .OR.  var(nzb,nysg,nxlg) == 0.0_wp )  CONTINUE

!
!-- Here, the respective variable is calculated. There is no call of exchange_horiz necessary.
!-- The variable identifyer (id) must be set according to the settings in user_init_flight.
!-- Please note, so far, variable must be located at the center of a grid box.
!     var = 0.0_wp

!     SELECT CASE ( id )
!
!        CASE ( 1 )
!           DO  i = nxl-1, nxr+1
!              DO  j = nys-1, nyn+1
!                 DO  k = nzb, nzt
!                    var(k,j,i) = ABS( u(k,j,i )
!                 ENDDO
!              ENDDO
!           ENDDO
!
!        CASE ( 2 )
!           DO  i = nxl-1, nxr+1
!              DO  j = nys-1, nyn+1
!                 DO  k = nzb, nzt
!                    var(k,j,i) = ABS( v(k,j,i) )
!                 ENDDO
!              ENDDO
!           ENDDO
!
!     END SELECT


 END SUBROUTINE user_flight
