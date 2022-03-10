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

!===================================================================================================================================
!> Contains global variables for the prepare WM mesh tool
!===================================================================================================================================
MODULE MOD_PrepareWMmesh_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! User defined parameters
LOGICAL             :: doVisu                            !< Visualize connection or not
LOGICAL             :: CartesianMode                     !< In cartesian mode, the elements are assumed to be aligned with with the 
                                                         !< wall-normal direction. This means we only interpolate the solution along
                                                         !< the wall-normal line.
LOGICAL             :: UseGaussPoints                    !< Use a gauss point as the interface location (except the first or last GL
                                                         !< point)
INTEGER             :: boundaryRegime                    !< Boundary layer regime (parameters below)
INTEGER             :: NSuper                            !< Polynomial degree for supersampling

INTEGER,PARAMETER   :: BLREGIME_LAMINAR = 1              !< Full Laminar BL
INTEGER,PARAMETER   :: BLREGIME_TRANSITION = 2           !< External BL: laminar near leading edge, transition to turbulent
INTEGER,PARAMETER   :: BLREGIME_TURBULENT = 3            !< Full Turbulent BL

INTEGER,PARAMETER   :: INTERFACESHAPE_CONSTANT    = 1    !< Interface in constant wall distance
INTEGER,PARAMETER   :: INTERFACESHAPE_LINEARX     = 2    !< Interface with linear variation in x
INTEGER,PARAMETER   :: INTERFACESHAPE_NACA64418   = 3    !< Interface for the NACA64418 interface
INTEGER,PARAMETER   :: INTERFACESHAPE_BLASIUS = 4        !< Interface considering Blasius laminar BL

! Parameters of interface shapes
REAL,ALLOCATABLE    :: interfaceShapeInfo(:,:)         !< Stores inferface shape info
                                                       !< (INFO,1:nShape), where INFO refers to shape, h_wm distance, xStart ,etc
                                                       !< Depending on INTERFACESHAPE, some of INFO might be dummy variables/numbers
REAL                :: xTransition                     !< Location of transition point for Transition BL Regime
! INTERFACESHAPE_CONSTANT 
REAL                :: interfaceDistance                 !< Constant distance of the interface away from the wall

! Connection information
REAL,ALLOCATABLE    :: wallconnect(:,:,:,:)              !< Stores the actual connection information
                                                         !< (INFO,0:PP_N+1,0:PP_N+1,1:nModelledBCSides)
INTEGER,ALLOCATABLE :: mapBCSideToModelledSide(:)        !< Takes side index and return index of modelled BC side
INTEGER             :: nModelledBCSides                  !< Number of modelled BC sides
INTEGER,ALLOCATABLE :: sideShapeInfo(:,:)                !< For each modelled side, contains info about the boundary
                                                         !< layer regime for this side, and whether it is a boundary
                                                         !< of the shape or not (needed for laminar BL model)
                                                         !< (2, nModelledBCSides): 1: Interface shape number; 2: Is Shape Boundary?

! Data needed for point search in mesh
REAL,ALLOCATABLE    :: XCL_Ngeo(:,:,:,:,:)               !< Mesh coordinates on GL points (Ngeo), later needed to find
                                                         !< parametric coordinates.
REAL,ALLOCATABLE    :: DCL_NGeo(:,:)                     !< Derivative matrix for GL mesh points
REAL,ALLOCATABLE    :: Xi_CLNgeo(:)                      !< Reference coordinates for CL points on Ngeo
REAL,ALLOCATABLE    :: wBary_CLNgeo(:)                   !< Barycentric weights for CL points on Ngeo
REAL,ALLOCATABLE    :: Xi_NSuper(:)                      !< Reference coordinates for supersampling points
REAL,ALLOCATABLE    :: Vdm_NGeo_NSuper(:,:)              !< Vandermonde matrix from geometry to supersampling points
REAL,ALLOCATABLE    :: Xi_N(:)                           !< Reference coordinates on calculation mesh

! Variable names of wallconnect array
CHARACTER(LEN=255),DIMENSION(17),PARAMETER  :: StrVarNamesWallconnect = (/ CHARACTER(LEN=255) :: &
                                                       'ElemSend',&
                                                       'ElemRcv',&
                                                       'XI',&
                                                       'ETA',&
                                                       'ZETA',&
                                                       'TangVec1_x',&
                                                       'TangVec1_y',&
                                                       'TangVec1_z',&
                                                       'TangVec2_x',&
                                                       'TangVec2_y',&
                                                       'TangVec2_z',&
                                                       'dir',&
                                                       'p',&
                                                       'q',&
                                                       'l',&
                                                       'BLRegime',&
                                                       'ShapeBoundary'/)

END MODULE MOD_PrepareWMmesh_Vars

