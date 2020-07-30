!=================================================================================================================================
! This file is NOT part of the original FLEXI project. Instead, it is an addition made by the LES group at
! Instituto Tecnológico de Aeronáutica, DCTA/IAE, Brazil, in order to implement wall modeling for LES computations.
!
! For more information on the original project, see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! ORIGINAL DISCLAIMER:
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================

#include "flexi.h"

!==================================================================================================================================
!> Module to implement wall modeling for LES computations, hence dealing with boundary condition operations
!==================================================================================================================================
MODULE MOD_WMLES
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE

INTEGER, PARAMETER :: WMLES_SCHUMANN = 1
INTEGER, PARAMETER :: WMLES_WERNERWANGLE = 2
INTEGER, PARAMETER :: WMLES_EQTBLE = 3


PUBLIC :: DefineParametersWMLES
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersWMLES()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Wall-Modeled LES")
CALL prms%CreateIntFromStringOption('WMLES', "Wall model to be used on walls defined with approximate boundary conditions")
CALL addStrListEntry('WMLES', 'Schumann', WMLES_SCHUMANN)
CALL addStrListEntry('WMLES', 'WernerWangle', WMLES_WERNERWANGLE)
CALL addStrListEntry('WMLES', 'EquilibriumTBLE', WMLES_EQTBLE)

CALL prms%CreateRealOption('H_WM', "Distance from the wall at which LES and Wall Model "// &
                                    "exchange instantaneous flow information", "0.2").

END SUBROUTINE DefineParametersWMLES

!==================================================================================================================================
!> Read and initialize parameters of WMLES computation
!==================================================================================================================================
SUBROUTINE InitWMLES()
! MODULES
USE MOD_Globals
USE MOD_WMLES_Vars
USE MOD_Equation_Vars
IMPLICIT NONE

!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

!==================================================================================================================================


END SUBROUTINE InitWMLES


!==================================================================================================================================
END MODULE MOD_WMLES