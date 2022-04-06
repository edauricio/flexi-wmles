!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz 
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "flexi.h"

!===================================================================================================================================
!> This program prepares the "mesh" needed for the wall model. It needs to connect each boundary gauss point that is on a
!> wall modelled surface to a point in a user specified wall-normal distance. The solution at this point is taken as an input
!> for the wall model. So during runtime, we need to communicate the solution at this point to the respective boundary point.
!===================================================================================================================================
PROGRAM PrepareWMmesh
! MODULES
USE MOD_Globals
USE MOD_Commandline_Arguments
USE MOD_PrepareWMmesh_Vars
USE MOD_ReadInTools
USE MOD_Mesh,                    ONLY: DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_MPI,                     ONLY: DefineParametersMPI,InitMPI
USE MOD_Interpolation,           ONLY: DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_IO_HDF5,                 ONLY: DefineParametersIO_HDF5,InitIOHDF5
#if USE_MPI
USE MOD_MPI,                     ONLY: InitMPIvars,FinalizeMPI
#endif
USE MOD_Connect,                 ONLY: InitConnect,Connect,ProjectToGaussPoints,FinalizeConnect
USE MOD_Connect_Output,          ONLY: ConnectHDF5Output,ConnectVisu
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL InitMPI() ! NO PARALLELIZATION, ONLY FOR COMPILING WITH MPI FLAGS ON SOME MACHINES OR USING MPI-DEPENDANT HDF5
 IF (nProcessors.GT.1) CALL abort(__STAMP__, &
      'This tool is designed only for single execution!')

CALL ParseCommandlineArguments()

! Define parameters needed
CALL DefineParametersInterpolation()
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
CALL DefineParametersMesh()

CALL prms%SetSection("prepareWMmesh")
CALL prms%CreateLogicalOption(      "doVisu"           , "Visualize the generated connection" ,".TRUE.")
CALL prms%CreateLogicalOption(      "cartesianMode"    , "Assume the elements are aligned with the wall-normal direction" ,".TRUE.")
CALL prms%CreateLogicalOption(      "useGaussPoints"   , "Use the nearest gauss point as the interface location, except for first"// &
                                                          " or last point GL point" ,".FALSE.")
CALL prms%CreateIntFromStringOption("boundaryLayerRegime", "Regime of boundary layer. Laminar/Transition/Turbulent")
CALL addStrListEntry("boundaryLayerRegime", "Laminar", BLREGIME_LAMINAR)
CALL addStrListEntry("boundaryLayerRegime", "Transition", BLREGIME_TRANSITION)
CALL addStrListEntry("boundaryLayerRegime", "Turbulent", BLREGIME_TURBULENT)
CALL prms%CreateRealOption("xTransitionPoint", "Streamwise location of boundary layer laminar-turbulent transition (for Transition BL Regime)")
CALL prms%CreateIntFromStringOption("interfaceShape"   , "Shape of the interface", multiple=.TRUE.)
CALL addStrListEntry(               'interfaceShape'   , 'constant', INTERFACESHAPE_CONSTANT)
CALL addStrListEntry(               'interfaceShape'   , 'linearx',  INTERFACESHAPE_LINEARX)
CALL addStrListEntry(               'interfaceShape'   , 'naca64418',INTERFACESHAPE_NACA64418)
CALL addStrListEntry(               'interfaceShape'   , 'blasius',INTERFACESHAPE_BLASIUS)
CALL addStrListEntry(               'interfaceShape'   , 'naca0012',INTERFACESHAPE_NACA0012)
CALL addStrListEntry(               'interfaceShape'   , 'blaslin',INTERFACESHAPE_BLASLIN)
CALL prms%CreateRealOption(         "interfaceDistance", "Distance of the interface from the wall (interface shape constant)", multiple=.TRUE.)
CALL prms%CreateRealOption(         "xStart"           , "First x position of linear function (interface shape linearX)", multiple=.TRUE.)
CALL prms%CreateRealOption(         "xEnd"             , "Last x position of linear function (interface shape linearX)", multiple=.TRUE.)
CALL prms%CreateRealOption(         "distanceStart"    , "Interface distance at first x position of linear function "//&
                                                         "(interface shape linearX)", multiple=.TRUE.)
CALL prms%CreateRealOption(         "distanceEnd"      , "Interface distance at last x position of linear function "//&
                                                         "(interface shape linearX)", multiple=.TRUE.)
CALL prms%CreateRealOption(         "Reynolds"         , "Unit-length Reynolds number to be considered for Blasius BL interface shape", multiple=.TRUE.)
CALL prms%CreateIntOption(          "NSuper",           "Supersampling polynomial degree","5")

! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF

! Parse parameters
CALL prms%read_options(Args(1))


CALL InitIOHDF5()


SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') &
" __          ____  __    _____ ____  _   _ _   _ ______ _____ _______    "
SWRITE(UNIT_stdOut,'(A)') &
" \\ \\        / /  \\/  |  / ____/ __ \\| \\ | | \\ | |  ____/ ____|__   __|   "
SWRITE(UNIT_stdOut,'(A)') &
"  \\ \\  /\\  / /| \\  / | | |   | |  | |  \\| |  \\| | |__ | |       | |      "
SWRITE(UNIT_stdOut,'(A)') &
"   \\ \\/  \\/ / | |\\/| | | |   | |  | | . ` | . ` |  __|| |       | |      "
SWRITE(UNIT_stdOut,'(A)') &
"    \\  /\\  /  | |  | | | |___| |__| | |\\  | |\\  | |___| |____   | |      "
SWRITE(UNIT_stdOut,'(A)') &
"     \\/  \\/   |_|  |_|  \\_____\\____/|_| \\_|_| \\_|______\\_____|  |_|      "
                                                                      
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(132("="))')

! Initialization
CALL InitInterpolation()
CALL InitMesh(meshMode=2)
#if USE_MPI
CALL InitMPIvars()
#endif
CALL  InitConnect()

! Call routine to build the actual connection
CALL Connect(cartesianMode)
IF (useGaussPoints) CALL ProjectToGaussPoints()
! Call routine to write the connection information to a HDF5 file
CALL ConnectHDF5Output()
! Call routine to visualize the result if needed
IF (doVisu) CALL ConnectVisu()

! Finalize
CALL FinalizeConnect()
CALL FinalizeInterpolation()
CALL FinalizeMesh()
CALL FinalizeParameters()
CALL FinalizeCommandlineArguments()
#if USE_MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) STOP 'MPI finalize error'
#endif

WRITE(UNIT_stdOut,'(132("="))')
WRITE(UNIT_stdOut,'(A)') ' PREPARE WM MESH TOOL FINISHED! '
WRITE(UNIT_stdOut,'(132("="))')

END PROGRAM PrepareWMmesh
