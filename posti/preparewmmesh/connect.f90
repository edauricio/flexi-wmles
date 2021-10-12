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
!> Module containing the main procedures needed to build the connection between the WM interface and the boundary points.
!===================================================================================================================================
MODULE MOD_Connect
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES

!INTEGER,PARAMETER :: INTERFACE_L = N_INTERFACE_PARAMS + 1
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitConnect
  MODULE PROCEDURE InitConnect
END INTERFACE

INTERFACE FinalizeConnect
  MODULE PROCEDURE FinalizeConnect
END INTERFACE

INTERFACE Connect
  MODULE PROCEDURE Connect
END INTERFACE

INTERFACE ProjectToGaussPoints
  MODULE PROCEDURE ProjectToGaussPoints
END INTERFACE

PUBLIC:: InitConnect,FinalizeConnect
PUBLIC:: Connect,ProjectToGaussPoints

CONTAINS

!===================================================================================================================================
!> Initialization of the connect routines. Get user defined parameters and allocate the necessary arrays.
!> Also build all the nodes, weights, Vdm's and so one needed by the point search later.
!===================================================================================================================================
SUBROUTINE InitConnect()
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_PreProc
USE MOD_PrepareWMmesh_Vars
USE MOD_ReadInTools
USE MOD_Mesh_Vars,            ONLY: nBCSides,BoundaryType,BC,NGeo,useCurveds,MeshFile,nElems,offsetElem
USE MOD_Interpolation_Vars,   ONLY: NodeTypeCL,NodeTypeVISU,NodeType
USE MOD_Interpolation,        ONLY: GetVandermonde,GetDerivativeMatrix,GetNodesAndWeights
USE MOD_IO_HDF5,              ONLY: OpenDataFile,CloseDataFile
USE MOD_HDF5_Input,           ONLY: ReadArray
USE MOD_ChangeBasis,          ONLY: ChangeBasis3D
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
! Space-separated list of input and output types. Use: (int|real|logical|...)_(in|out|inout)_dim(n)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iSide,locType,iElem
REAL,ALLOCATABLE  :: NodeCoordsTmp(:,:,:,:,:),NodeCoords(:,:,:,:,:)
REAL,ALLOCATABLE  :: Vdm_EQNGeo_CLNGeo(:,:)
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PREPARE WM MESH TOOL...'

! Get user defined parameters
doVisu         = GETLOGICAL("doVisu")
interfaceShape = GETINTFROMSTR("interfaceShape")
CartesianMode  = GETLOGICAL("cartesianMode")
IF (cartesianMode) UseGaussPoints  = GETLOGICAL("UseGaussPoints")

NSuper         = GETINT("NSuper")

SELECT CASE(interfaceShape)
CASE(INTERFACESHAPE_CONSTANT)
  interfaceDistance = GETREAL("interfaceDistance")
CASE(INTERFACESHAPE_LINEARX)
  xStart        = GETREAL("xStart")               ! First x position of linear function
  xEnd          = GETREAL("xEnd")                 ! Last x position of linear function
  distanceStart = GETREAL("distanceStart")        ! Interface distance at first x position of linear function
  distanceEnd   = GETREAL("distanceEnd")          ! Interface distance at last x position of linear function
CASE(INTERFACESHAPE_NACA64418)
  ! Nothing to init here
CASE DEFAULT
  CALL Abort(__STAMP__,"Selected interface shape not implemented!")
END SELECT

! Build mapping from all BC sides to modelled BC sides
ALLOCATE(mapBCSideToModelledSide(nBCSides))
mapBCSideToModelledSide = 0
nModelledBCSides = 0
DO iSide=1,nBCSides
  locType =BoundaryType(BC(iSide),BC_TYPE)
  IF (locType.EQ.5) THEN
    nModelledBCSides = nModelledBCSides + 1
    mapBCSideToModelledSide(iSide) = nModelledBCSides
  END IF
END DO

SWRITE(*,*) 'Number of wall modelled sides: ',nModelledBCSides

! Allocate the array that is used to store the wallconnection information
ALLOCATE(wallconnect(15,0:PP_N,0:PP_N,nModelledBCSides))
wallconnect = 0.

!==================== Prepare data for point search ========================!

! Store the mesh coordinates for later
CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

IF(useCurveds)THEN
  ALLOCATE(NodeCoords(3,0:NGeo,0:NGeo,0:NGeo,nElems))
  CALL ReadArray('NodeCoords',2,(/3,nElems*(NGeo+1)**3/),offsetElem*(NGeo+1)**3,2,RealArray=NodeCoords)
ELSE
  ALLOCATE(NodeCoords(   3,0:1,   0:1,   0:1,   nElems))
  ALLOCATE(NodeCoordsTmp(3,0:NGeo,0:NGeo,0:NGeo,nElems))
  CALL ReadArray('NodeCoords',2,(/3,nElems*(NGeo+1)**3/),offsetElem*(NGeo+1)**3,2,RealArray=NodeCoordsTmp)
  NodeCoords(:,0,0,0,:)=NodeCoordsTmp(:,0,   0,   0,   :)
  NodeCoords(:,1,0,0,:)=NodeCoordsTmp(:,NGeo,0,   0,   :)
  NodeCoords(:,0,1,0,:)=NodeCoordsTmp(:,0,   NGeo,0,   :)
  NodeCoords(:,1,1,0,:)=NodeCoordsTmp(:,NGeo,NGeo,0,   :)
  NodeCoords(:,0,0,1,:)=NodeCoordsTmp(:,0,   0,   NGeo,:)
  NodeCoords(:,1,0,1,:)=NodeCoordsTmp(:,NGeo,0,   NGeo,:)
  NodeCoords(:,0,1,1,:)=NodeCoordsTmp(:,0,   NGeo,NGeo,:)
  NodeCoords(:,1,1,1,:)=NodeCoordsTmp(:,NGeo,NGeo,NGeo,:)
  DEALLOCATE(NodeCoordsTmp)
  NGeo=1
ENDIF

CALL CloseDataFile()

! Vandermonde from EQUI NGeo to CLNGeo
ALLOCATE(Vdm_EQNGeo_CLNGeo(0:NGeo,0:NGeo))
CALL GetVandermonde(NGeo,NodeTypeVISU,NGeo,NodeTypeCL,Vdm_EQNGeo_CLNGeo,modal=.FALSE.)

! Get node coordinates on CL points
ALLOCATE(XCL_NGeo(3,0:NGeo,0:NGeo,0:NGeo,nElems))
DO iElem=1,nElems
  CALL ChangeBasis3D(3,NGeo,NGeo,Vdm_EQNGeo_CLNGeo,NodeCoords(:,:,:,:,iElem),XCL_NGeo(:,:,:,:,iElem))
END DO ! iElem
SDEALLOCATE(Vdm_EQNGeo_CLNGeo)

! Derivative matrix for dXCL_NGeo
ALLOCATE(DCL_NGeo(0:NGeo,0:NGeo))
CALL GetDerivativeMatrix(NGeo,NodeTypeCL,DCL_NGeo)

! Reference nodes and weights for the geometry on CL points
ALLOCATE(Xi_CLNGeo(   0:NGeo))
ALLOCATE(wBary_CLNGeo(0:NGeo))
CALL GetNodesAndWeights(NGeo,NodeTypeCL,Xi_CLNGeo,wIPBary=wBary_CLNGeo)

! Reference nodes on supersampling points
ALLOCATE(Xi_NSuper(0:NSuper))
CALL GetNodesAndWeights(NSuper,NodeTypeCL,Xi_NSuper)

! Vandermonde from geometry to supersampling points
ALLOCATE(Vdm_NGeo_NSuper(0:NSuper,0:NGeo))
CALL GetVandermonde(NGeo,NodeTypeCL,NSuper,NodeTypeCL,Vdm_NGeo_NSuper,modal=.FALSE.)

! Reference nodes on calculation mesh
ALLOCATE(Xi_N(0:PP_N))
CALL GetNodesAndWeights(PP_N,NodeType,Xi_N)

SWRITE(UNIT_stdOut,'(A)')' INIT PREPARE WM MESH TOOL DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')


END SUBROUTINE InitConnect

!===================================================================================================================================
!> Build the connection between the wall modelled boundary points and the respective interface point.
!> To do this we:
!>   * Loop over all wall modelled boundary points
!>   * Store the tangential vectors for this point as well as the element it belongs to (= the element the data needs to
!>      be send to)
!>   * Calculate the height of the wall model layer at this point, depending on the choosen shape of the layer
!>   * Calculate the physical position of the interface, which is the position of the boundary point plus the height of the wall
!>      model layer in the wall-normal direction
!>   * Find the element and the reference coordinates of the calculated interface using a point search algorithm
!>   * If we are in cartesian mode, the interpolation of the values at the interface will be done in only one reference direction.
!>     This direction and the other two interfaces will be calculated. Otherwise, a full 3D interpolation must be performed.
!===================================================================================================================================
SUBROUTINE Connect(cartesian)
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_PreProc
USE MOD_PrepareWMmesh_Vars
USE MOD_Mesh_Vars,          ONLY: nBCSides,SideToElem,TangVec1,TangVec2,NormVec,Face_xGP,NGeo,nElems,Elem_xGP,ElemToSide
USE MOD_Mesh_Vars,          ONLY: BoundaryName,BC
USE MOD_Stringtools,        ONLY: STRICMP
USE MOD_TestCase_Vars,      ONLY: testcase
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
LOGICAL,INTENT(IN) :: cartesian
! Space-separated list of input and output types. Use: (int|real|logical|...)_(in|out|inout)_dim(n)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iSide,iSideModelled,ElemID,p,q,ElemID_send,i,refDirection,SideID,ElemID2
INTEGER           :: nextLocSide
REAL              :: InterfaceCoordinates(3),wm_l
REAL              :: ParamCoords(3)
REAL              :: biggestScalProd,scalProd(3),xi_vec(3),eta_vec(3),zeta_vec(3)
REAL              :: x
REAL              :: percent
INTEGER           :: nSidesDone
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)')' Starting to connect the boundary points to the interface...'
nSidesDone = 0

! Loop over all BC sides
DO iSide = 1,nBCSides
  ! Map to index of the modelled sides
  iSideModelled = mapBCSideToModelledSide(iSide)
  ! Check if this side is modelled
  IF (iSideModelled.GT.0) THEN
    nSidesDone = nSidesDone + 1
    percent = REAL(nSidesDone)/REAL(nModelledBCSides)*100.
    WRITE(UNIT_stdOut,'(F7.2,A16,A1)',ADVANCE='NO') percent, '% of sides done ', ACHAR(13)

    ElemID = SideToElem(S2E_ELEM_ID,iSide)
    DO q=0,PP_N; DO p=0,PP_N
      ! ID of element that is on the recieving side
      wallconnect(INTERFACE_ELEMRECV,p,q,iSideModelled) = ElemID
      ! Both tangential vectors
      wallconnect(INTERFACE_TANGVEC1,p,q,iSideModelled) = TangVec1(:,p,q,0,iSide)
      wallconnect(INTERFACE_TANGVEC2,p,q,iSideModelled) = TangVec2(:,p,q,0,iSide)
      ! Distance to the interface
      SELECT CASE(interfaceShape)
      CASE(INTERFACESHAPE_CONSTANT)
        wm_l = interfaceDistance
      CASE(INTERFACESHAPE_LINEARX)
        x = Face_xGP(1,p,q,0,iSide)
        wm_l = distanceStart + (x-xStart)/(xEnd-xStart)*(distanceEnd-distanceStart)
      CASE(INTERFACESHAPE_NACA64418)
        x = Face_xGP(1,p,q,0,iSide)
        ! Check if we are on the lower or upper wall by comparing the BC name
        ! On both sides, linear variation of the interface height such that h \approx 0.15 \delta_99
        IF (STRICMP(BoundaryName(BC(iSide)),'BC_wall_U')) THEN
          ! Upper wall
          wm_l = 9.45811E-3 * x - 3.45811E-3
        ELSE IF (STRICMP(BoundaryName(BC(iSide)),'BC_wall_L')) THEN
          ! Lower wall
          wm_l = 4.39865E-3 * x - 1.09865E-3
        ELSE
          CALL Abort(__STAMP__,"Wrong BC names for NACA64418 test case!")
        END IF
      CASE DEFAULT
        CALL Abort(__STAMP__,"Selected interface shape not implemented!")
      END SELECT
      wallconnect(INTERFACE_L,p,q,iSideModelled) = wm_l
      ! Coordinates of the interface = coordinate of boundary points plus distance in wall-normal direction (NormVec is outwards)
      InterfaceCoordinates = Face_xGP(:,p,q,0,iSide) - wm_l*NormVec(:,p,q,0,iSide)

#if (PP_NodeType==2)
      ! For the channel testcase with GL points, to avoid a "jump" of the interface on the double-valued cell borders, we
      ! determine the element that has the interface "manually" and only search in this single element
      IF (testcase.EQ."channel") THEN
        ! Check if we are on the upper or lower half of the channel
        IF (Elem_xGP(2,0,INT(PP_N/2),0,ElemID).GT.0.) THEN
          ! Upper half, move towards ETA_MINUS sides
          nextLocSide = ETA_MINUS
        ELSE
          ! Lower half, move towards ETA_PLUS sides
          nextLocSide = ETA_PLUS
        END IF

        ! Loop until we reach the element containing the interface
        ElemID2 = ElemID
        DO WHILE ((1.-ABS(Elem_xGP(2,0,0,0,ElemID2))-wm_l)*(1.-ABS(Elem_xGP(2,0,PP_N,0,ElemID2))-wm_l).GT.0)
          SideID = ElemToSide(E2S_SIDE_ID,nextLocSide,ElemID2)
          IF (SideToElem(S2E_ELEM_ID,SideID).EQ.ElemID2) THEN
            ! The current element is master to this side, move to the slave element
            ElemID2 = SideToElem(S2E_NB_ELEM_ID,SideID)
          ELSE
            ! The current element is slave to this side, move to the master element
            ElemID2 = SideToElem(S2E_ELEM_ID,SideID)
          END IF
        END DO
        ! Search for the interface in the mesh
        CALL SearchForParametricCoordinates(XCL_NGeo(:,:,:,:,ElemID2:ElemID2),NGeo,1,DCL_NGeo,Xi_CLNGeo,wBary_CLNGeo,NSuper,Vdm_NGeo_NSuper,&
                Xi_NSuper,InterfaceCoordinates,ParamCoords,ElemID_send)
        ! ID of element that has to send
        wallconnect(INTERFACE_ELEMSEND,p,q,iSideModelled) = ElemID2
      ELSE
#endif
        ! Search for the interface in the mesh
        CALL SearchForParametricCoordinates(XCL_NGeo,NGeo,nElems,DCL_NGeo,Xi_CLNGeo,wBary_CLNGeo,NSuper,Vdm_NGeo_NSuper,&
                Xi_NSuper,InterfaceCoordinates,ParamCoords,ElemID_send,ElemID)
        ! ID of element that has to send
        wallconnect(INTERFACE_ELEMSEND,p,q,iSideModelled) = ElemID_send
#if (PP_NodeType==2)
      END IF
#endif

      ! Reference coordinates
      wallconnect(INTERFACE_XI  ,p,q,iSideModelled) = ParamCoords(1)
      wallconnect(INTERFACE_ETA ,p,q,iSideModelled) = ParamCoords(2)
      wallconnect(INTERFACE_ZETA,p,q,iSideModelled) = ParamCoords(3)

      IF (cartesian) THEN

        ! For the cartesian mesh, the interpolation will be done on a single line. We now need to find out which reference
        ! direction we need to do the interpolation along. To do this, we compare unit vectors of the reference directions in
        ! physical space with the normal vector of the boundary point we are looking at. The direction which the greatest alignment
        ! along the normal vector will be the direction to interpolate along.

        ! Get reference unit vectors
        xi_vec   = XCL_NGeo(:,NGeo,0,0,ElemID_send) - XCL_NGeo(:,0,0,0,ElemID_send)
        eta_vec  = XCL_NGeo(:,0,NGeo,0,ElemID_send) - XCL_NGeo(:,0,0,0,ElemID_send)
        zeta_vec = XCL_NGeo(:,0,0,NGeo,ElemID_send) - XCL_NGeo(:,0,0,0,ElemID_send)

        ! Normalize reference unit vectors
        xi_vec   = xi_vec/  SQRT(SUM(xi_vec(:)*  xi_vec(:)))
        eta_vec  = eta_vec/ SQRT(SUM(eta_vec(:)* eta_vec(:)))
        zeta_vec = zeta_vec/SQRT(SUM(zeta_vec(:)*zeta_vec(:)))

        ! Scalar product between reference unit vectors and normal vector of boundary point = alignment between the two
        scalProd(1) = ABS(SUM(xi_vec(:)  *NormVec(:,p,q,0,iSide)))
        scalProd(2) = ABS(SUM(eta_vec(:) *NormVec(:,p,q,0,iSide)))
        scalProd(3) = ABS(SUM(zeta_vec(:)*NormVec(:,p,q,0,iSide)))

        biggestScalProd = 0.
        DO i=1,3
          IF (scalProd(i).GT.biggestScalProd) THEN
            refDirection = i
            biggestScalProd = scalProd(i)
          END IF
        END DO ! i=0,3

        ! Save the reference direction
        wallconnect(INTERFACE_DIR,p,q,iSideModelled) = refDirection

        ! Now get the other two volume indizes
        ! For this, we compare the reference coordinate against the positions of the gauss points and take the index of the closest
        ! one
        SELECT CASE(refDirection)
        CASE(1)
          ! Interpolation is done along xi direction, get eta and zeta indizes
          wallconnect(INTERFACE_P,p,q,iSideModelled) = MINLOC(ABS(Xi_N-ParamCoords(2)),1) - 1
          wallconnect(INTERFACE_Q,p,q,iSideModelled) = MINLOC(ABS(Xi_N-ParamCoords(3)),1) - 1
        CASE(2)
          ! Interpolation is done along eta direction, get xi and zeta indizes
          wallconnect(INTERFACE_P,p,q,iSideModelled) = MINLOC(ABS(Xi_N-ParamCoords(1)),1) - 1
          wallconnect(INTERFACE_Q,p,q,iSideModelled) = MINLOC(ABS(Xi_N-ParamCoords(3)),1) - 1
        CASE(3)
          ! Interpolation is done along zeta direction, get xi and eta indizes
          wallconnect(INTERFACE_P,p,q,iSideModelled) = MINLOC(ABS(Xi_N-ParamCoords(1)),1) - 1
          wallconnect(INTERFACE_Q,p,q,iSideModelled) = MINLOC(ABS(Xi_N-ParamCoords(2)),1) - 1
        END SELECT

        ! In the calculation, to get the solution at the interface, we must do a interpolation of (e.g. for refDirection = XI)
        ! U(:,INTERFACE_P,INTERFACE_Q,INTERFACE_ELEMSEND) at the position INTERFACE_XI.

      END IF ! cartesian

      ! For non-cartesian mode, the interpolation will be done in full 3D, so U(:,:,:,INTERFACE_ELEMSEND) at
      ! (INTERFACE_XI,INTERFACE_ETA,INTERFACE_ZETA)

    END DO; END DO

  END IF ! Modelled side
END DO ! iSideModelled = 1,nModelledBCSides

SWRITE(UNIT_stdOut,'(A)')' DONE!'

END SUBROUTINE Connect

!===================================================================================================================================
!> After a cartesian mesh connect has been done, this routine will take the result and project the interface onto the nearest gauss
!> point. For GL points, the points on the interface are excluded.
!===================================================================================================================================
SUBROUTINE ProjectToGaussPoints()
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_PreProc
USE MOD_PrepareWMmesh_Vars
USE MOD_Mesh_Vars,          ONLY: nBCSides,NormVec,Face_xGP,Elem_xGP
USE MOD_Interpolation,      ONLY: GetNodesAndWeights
USE MOD_Interpolation_Vars, ONLY: NodeType
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iSide,iSideModelled,p,q,i,refDirection
REAL              :: xi(0:PP_N),best,dist,newInterface,refCoordinate
INTEGER           :: refPoint,NormVecDir
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)')' Starting to project the interface onto gauss points...'

! Get 1D solution nodes
CALL GetNodesAndWeights(PP_N,NodeType,xi)

! Loop over all BC sides
DO iSide = 1,nBCSides
  ! Map to index of the modelled sides
  iSideModelled = mapBCSideToModelledSide(iSide)
  ! Check if this side is modelled
  IF (iSideModelled.GT.0) THEN
    DO q=0,PP_N; DO p=0,PP_N
      ! Read cartesian direction of interpolation and the corresponding reference coordinate
      refDirection = NINT(wallconnect(INTERFACE_DIR,p,q,iSideModelled))
      SELECT CASE(refDirection)
      CASE(1)
        refCoordinate = wallconnect(INTERFACE_XI  ,p,q,iSideModelled)
      CASE(2)
        refCoordinate = wallconnect(INTERFACE_ETA ,p,q,iSideModelled)
      CASE(3)
        refCoordinate = wallconnect(INTERFACE_ZETA,p,q,iSideModelled)
      END SELECT
      ! Compare the reference coordinate with the 1D nodes and find the closest one
      best = HUGE(1.)
      DO i = 0, PP_N
        dist = ABS(refCoordinate-xi(i))
        IF (dist.LT.best) THEN
          refPoint = i
          best = dist
        END IF
      END DO ! i = 0, PP_N
#if (PP_NodeType==2)
      ! For GL, exclude the points on the boundary
      refPoint = MAX(1,refPoint)
      refPoint = MIN(refPoint,PP_N-1)
#endif
      ! Overwrite the reference coordinate and calculate the new wall distance
      NormVecDir = MAXLOC(ABS(NormVec(:,p,q,0,iSide)),1)
      SELECT CASE(refDirection)
      CASE(1)
        wallconnect(INTERFACE_XI  ,p,q,iSideModelled) = xi(refPoint)
        newInterface = Elem_xGP(NormVecDir,refPoint,0,0,NINT(wallconnect(INTERFACE_ELEMSEND,p,q,iSideModelled)))
      CASE(2)
        wallconnect(INTERFACE_ETA ,p,q,iSideModelled) = xi(refPoint)
        newInterface = Elem_xGP(NormVecDir,0,refPoint,0,NINT(wallconnect(INTERFACE_ELEMSEND,p,q,iSideModelled)))
      CASE(3)
        wallconnect(INTERFACE_ZETA,p,q,iSideModelled) = xi(refPoint)
        newInterface = Elem_xGP(NormVecDir,0,0,refPoint,NINT(wallconnect(INTERFACE_ELEMSEND,p,q,iSideModelled)))
      END SELECT
      wallconnect(INTERFACE_L,p,q,iSideModelled) = (Face_xGP(NormVecDir,p,q,0,iSideModelled) - newInterface) * &
                                             SIGN(1.,NormVec(NormVecDir,p,q,0,iSide))
    END DO; END DO
  END IF ! Modelled side
END DO ! iSideModelled = 1,nModelledBCSides

SWRITE(UNIT_stdOut,'(A)')' DONE!'

END SUBROUTINE ProjectToGaussPoints

!=================================================================================================================================
!> This routine takes a mesh (coordinates in CL nodes) and the physical coordinates of a single point and searches for the
!> element that contains this single point. It will then return the ID of this element and the coordinates in reference space of
!> this point.
!=================================================================================================================================
SUBROUTINE SearchForParametricCoordinates(NodeCoords,NGeo,nElems,DCL_NGeo,Xi_CLNGeo,wBary_CLNGeo,NSuper,Vdm_NGeo_NSuper,&
                                          Xi_NSuper,Coordinates,ParametricCoords,ElemID,firstGuess)
! MODULES
USE MOD_Globals
USE MOD_ChangeBasis,     ONLY: ChangeBasis3D
USE MOD_Basis,           ONLY: LagrangeInterpolationPolys
USE MOD_Mathtools,       ONLY: INVERSE
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT)  :: NodeCoords(3,0:NGeo,0:NGeo,0:NGeo,nElems) !< Mesh coordinates
INTEGER,INTENT(IN)  :: NGeo                                      !< Polynomial degree of mesh representation
INTEGER,INTENT(IN)  :: nElems                                    !< Number of elements in mesh
REAL,INTENT(IN)     :: DCL_NGeo(0:NGeo,0:NGeo)                   !< Derivative matrix on NGeo
REAL,INTENT(IN)     :: Xi_CLNGeo(0:NGeo)                         !< Reference coordinates on geometry
REAL,INTENT(IN)     :: wBary_CLNGeo(0:NGeo)                      !< Barycentric weights on geometry
INTEGER,INTENT(IN)  :: NSuper                                    !< Polynomial degree for supersampling
REAL,INTENT(IN)     :: Vdm_NGeo_NSuper(0:NSuper,0:NGeo)          !< Vandermonde to get mesh points on NSuper
REAL,INTENT(IN)     :: Xi_NSuper(0:NSuper)                       !< Reference coordinates of supersampled values
REAL,INTENT(IN)     :: Coordinates(3)                            !< Physical coordinates to search for
INTEGER,INTENT(IN),OPTIONAL :: firstGuess                        !< Guess of the correct element id
!---------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: ParametricCoords(3)                       !< Coordinates of the points in reference space
INTEGER,INTENT(OUT) :: ElemID                                    !< ID of the element that contains the point
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,i,j,k,iter,ll
REAL                :: center(3),radius,dist,smallesDist
REAL                :: NodeCoordsSuperSamp(3,0:NSuper,0:NSuper,0:NSuper),CoordinateSuperSamp(3)
REAL                :: Xi(3),Lag(1:3,0:NGeo),Jac(1:3,1:3),sJac(1:3,1:3),F(3),eps_F,sdetJac
LOGICAL             :: found
REAL                :: dXCL_NGeo(  3,3,0:NGeo,0:NGeo,0:NGeo)          ! jacobi matrix on CL NGeo
!=================================================================================================================================
found = .FALSE.
! First test of a guessed element
IF (PRESENT(firstGuess)) THEN
  CALL SearchAlgorithm(firstGuess)
END IF
IF (.NOT.found) THEN
  ! Loop over all elements
  DO iElem=1,nElems
    CALL SearchAlgorithm(iElem)
    IF (Found) EXIT
  END DO ! iElem
END IF

IF (.NOT.found) THEN
  WRITE(*,*) 'Physical coordinates to find: ',Coordinates
  CALL Abort(__STAMP__,&
                      "Could not find interface location!")
END IF

CONTAINS
  !===================================================================================================================================
  !> The actual search algorithm, based on a coarse search and Newton's method
  !===================================================================================================================================
  SUBROUTINE SearchAlgorithm(iElem)
  !----------------------------------------------------------------------------------------------------------------------------------!
  IMPLICIT NONE
  ! INPUT / OUTPUT VARIABLES
  INTEGER,INTENT(IN) :: iElem
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
  !===================================================================================================================================
  ! Get the center and radius of the current element. This tells us if we need to consider the element in our search.
  CALL getCentroidAndRadius(NodeCoords(:,:,:,:,iElem),NGeo,center,radius)

  ! Check if the search point is in the radius of the current element - if not, cycle to next element
  IF ((ABS(Coordinates(1)-center(1)).GT.radius*1.05).OR.&
      (ABS(Coordinates(2)-center(2)).GT.radius*1.05).OR.&
      (ABS(Coordinates(3)-center(3)).GT.radius*1.05)) RETURN

  !Jacobi Matrix of d/dxi_dd(X_nn): dXCL_NGeo(dd,nn,i,j,k))
  dXCL_NGeo=0.
  DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
    ! Matrix-vector multiplication
    DO ll=0,NGeo
      dXCL_NGeo(:,1,i,j,k)=dXCL_NGeo(:,1,i,j,k) + DCL_NGeo(i,ll)*NodeCoords(:,ll,j,k,iElem)
      dXCL_NGeo(:,2,i,j,k)=dXCL_NGeo(:,2,i,j,k) + DCL_NGeo(j,ll)*NodeCoords(:,i,ll,k,iElem)
      dXCL_NGeo(:,3,i,j,k)=dXCL_NGeo(:,3,i,j,k) + DCL_NGeo(k,ll)*NodeCoords(:,i,j,ll,iElem)
    END DO !l=0,N
  END DO; END DO; END DO !i,j,k=0,NGeo

  ! Get a initial guess for the Newton algorithm by finding the closest supersampled grid point
  CALL ChangeBasis3D(3,NGeo,NSuper,Vdm_NGeo_NSuper,NodeCoords(:,:,:,:,iElem),NodeCoordsSuperSamp(:,:,:,:))
  smallesDist = HUGE(1.)
  DO k=0,NSuper; DO j=0,NSuper; DO i=0,NSuper
    CoordinateSuperSamp = NodeCoordsSuperSamp(:,i,j,k)
    ! Calculate distance from current supersampling point to search coordinate
    dist = SUM((Coordinates(:)-CoordinateSuperSamp(:))*(Coordinates(:)-CoordinateSuperSamp(:)))
    ! Save this point if it has the smallest distance
    IF (dist.LT.smallesDist) THEN
      smallesDist = Dist
      Xi(:)=(/Xi_NSuper(i),Xi_NSuper(j),Xi_NSuper(k)/)
    END IF
  END DO; END DO; END DO! i,j,k=0,PP_N

  ! Evaluate the values of the lagrange polynomials at the inital best position
  CALL LagrangeInterpolationPolys(Xi(1),NGeo,Xi_CLNGeo,wBary_CLNGeo,Lag(1,:))
  CALL LagrangeInterpolationPolys(Xi(2),NGeo,Xi_CLNGeo,wBary_CLNGeo,Lag(2,:))
  CALL LagrangeInterpolationPolys(Xi(3),NGeo,Xi_CLNGeo,wBary_CLNGeo,Lag(3,:))

  ! Now start a Newton iteration to find the root of the function F(xi) = x(xi) - Coordinate
  ! For this, do a first evaluation of the function F
  F=-Coordinates
  DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
    F=F+NodeCoords(:,i,j,k,iElem)*Lag(1,i)*Lag(2,j)*Lag(3,k)
  END DO; END DO; END DO

  ! Start the loop until convergence is reached
  eps_F=1.E-8*SUM(F*F) ! relative error to initial guess
  iter=0
  DO WHILE ((SUM(F*F).GT.eps_F).AND.(iter.LT.100))
    iter=iter+1
    ! Compute F Jacobian dx/dXi
    Jac=0.
    DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
      Jac=Jac+dXCL_NGeo(:,:,i,j,k)*Lag(1,i)*Lag(2,j)*Lag(3,k)
    END DO; END DO; END DO

    ! Compute inverse of Jacobian
    sdetJac=getDet(Jac)
    IF(sdetJac.NE.0.) THEN
     sdetJac=1./sdetJac
    ELSE
      ! Seems like newton has not converged
      ! allow Newton to fail without aborting, may happen when far outside of reference space [-1,1]
      WRITE(UNIT_stdOut,*)' Newton has not converged after ',iter, ' iterations! skipping...'
      WRITE(UNIT_stdOut,*)' Xi,Eta,Zeta = ', Xi
      !CALL abort(__STAMP__, &
          !'Newton method is singular!')
      EXIT
    ENDIF
    sJac=INVERSE(Jac)

    ! Iterate Xi using Newton step
    Xi = Xi - MATMUL(sJac,F)
    ! if Newton gets outside reference space range [-1,1], exit.
    ! But allow for some oscillation in the first couple of iterations, as we may discard the correct point/element!!
    IF((iter.GT.4).AND.(ANY(ABS(Xi).GT.1.2))) EXIT

    ! Compute function value
    CALL LagrangeInterpolationPolys(Xi(1),NGeo,Xi_CLNGeo,wBary_CLNGeo,Lag(1,:))
    CALL LagrangeInterpolationPolys(Xi(2),NGeo,Xi_CLNGeo,wBary_CLNGeo,Lag(2,:))
    CALL LagrangeInterpolationPolys(Xi(3),NGeo,Xi_CLNGeo,wBary_CLNGeo,Lag(3,:))
    ! F(xi) = x(xi) - xInter
    F=-Coordinates
    DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
      F=F+NodeCoords(:,i,j,k,iElem)*Lag(1,i)*Lag(2,j)*Lag(3,k)
    END DO; END DO; END DO
  END DO !newton

  ! check if Newton has converged and found a point inside the current element
  IF((SUM(F*F).LT.eps_F).OR.(SUM(F*F).LT.1E-10)) THEN
    IF(MAXVAL(ABS(Xi)).LE.1.0000001) THEN
      ! If so, return this point
      ParametricCoords = Xi
      ElemID = iElem
      found = .TRUE.
      RETURN
    END IF
  END IF
  END SUBROUTINE


END SUBROUTINE SearchForParametricCoordinates

!=================================================================================================================================
!> Computes the centroid and the radius of a single element
!=================================================================================================================================
SUBROUTINE getCentroidAndRadius(Coordinates,NGeo,center,radius)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: NGeo                                 !< Polynomial degree of mesh representation
REAL,INTENT(IN)    :: Coordinates(3,0:NGeo,0:NGeo,0:NGeo)  !< Mesh coordinates of a single element
!---------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: center(3)                            !< Coordinate of the center of the element
REAL,INTENT(OUT)   :: radius                               !< Biggest radius of the element
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,nNodes
LOGICAL            :: onSide(3)
!=================================================================================================================================
! Compute centroid of Coordinate
nNodes= (NGeo+1)**3
center(1) = SUM(Coordinates(1,:,:,:))/nNodes
center(2) = SUM(Coordinates(2,:,:,:))/nNodes
center(3) = SUM(Coordinates(3,:,:,:))/nNodes

! Compute max distance from bary to surface nodes
radius=0.
DO k=0,NGeo
  onSide(3)=((k.EQ.0).OR.(k.EQ.NGeo))
  DO j=0,NGeo
    onSide(2)=((j.EQ.0).OR.(j.EQ.NGeo))
    DO i=0,NGeo
      onSide(1)=((i.EQ.0).OR.(i.EQ.NGeo))
      IF(.NOT.ANY(onSide)) CYCLE
      radius=MAX(radius,NORM2(Coordinates(:,i,j,k)-center))
    END DO
  END DO
END DO

END SUBROUTINE getCentroidAndRadius

!=================================================================================================================================
!> Compute determinant of 3x3 matrix
!=================================================================================================================================
PURE FUNCTION getDet(Mat)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: Mat(3,3)
!---------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL             :: getDet
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!=================================================================================================================================
getDet=   ( Mat(1,1) * Mat(2,2) - Mat(1,2) * Mat(2,1) ) * Mat(3,3) &
        + ( Mat(1,2) * Mat(2,3) - Mat(1,3) * Mat(2,2) ) * Mat(3,1) &
        + ( Mat(1,3) * Mat(2,1) - Mat(1,1) * Mat(2,3) ) * Mat(3,2)
END FUNCTION getDet


!===================================================================================================================================
!> Finalization of connect routines, deallocate etc.
!===================================================================================================================================
SUBROUTINE FinalizeConnect()
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_PrepareWMmesh_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
! Space-separated list of input and output types. Use: (int|real|logical|...)_(in|out|inout)_dim(n)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(wallconnect)
SDEALLOCATE(mapBCSideToModelledSide)
SDEALLOCATE(DCL_NGeo)
SDEALLOCATE(Xi_CLNGeo)
SDEALLOCATE(wBary_CLNGeo)
SDEALLOCATE(XCL_NGeo)
SDEALLOCATE(Xi_NSuper)
SDEALLOCATE(Vdm_NGeo_NSuper)

END SUBROUTINE FinalizeConnect

END MODULE MOD_Connect
