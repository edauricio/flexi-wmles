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
#if USE_MPI
USE MPI 
#endif

IMPLICIT NONE
PRIVATE

PUBLIC :: DefineParametersWMLES, InitWMLES, FinalizeWMLES
PUBLIC :: ComputeWallStress, EvalDiffFlux3D_WMLES
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersWMLES()
! MODULES
USE MOD_Globals
USE MOD_WMLES_Vars
USE MOD_ReadInTools             ,ONLY: prms, addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Wall-Modeled LES")
CALL prms%CreateIntFromStringOption('WMLES', "Wall model to be used on walls defined with approximate boundary conditions")
CALL addStrListEntry('WMLES', 'Schumann', WMLES_SCHUMANN)
CALL addStrListEntry('WMLES', 'WernerWangle', WMLES_WERNERWANGLE)
CALL addStrListEntry('WMLES', 'EquilibriumTBLE', WMLES_EQTBLE)

CALL prms%CreateRealOption('h_wm', "Distance from the wall at which LES and Wall Model "// &
                                    "exchange instantaneous flow information", "0.2")
CALL prms%CreateRealOption('delta', "Estimated Boundary Layer Thickness")

END SUBROUTINE DefineParametersWMLES

!==================================================================================================================================
!> Read and initialize parameters of WMLES computation
!==================================================================================================================================
SUBROUTINE InitWMLES()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_HDF5_Input
USE MOD_WMLES_Vars
USE MOD_Interpolation_Vars      ,ONLY: InterpolationInitIsDone
USE MOD_Mesh_Vars              ! ,ONLY: nBCSides, BC, BoundaryType, MeshInitIsDone, ElemToSide, SideToElem, SideToGlobalSide
USE MOD_ReadInTools             ,ONLY: GETINTFROMSTR, GETREAL
USE MOD_StringTools             ,ONLY: STRICMP
USE MOD_Testcase_Vars           ,ONLY: Testcase
#if USE_MPI
USE MOD_MPI_Vars                !,ONLY:
#endif

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=10)               :: Channel="channel"
REAL                            :: WMLES_Tol, RealTmp
INTEGER                         :: SideID, iSide, LocSide, MineCnt, nSideIDs, offsetSideID
INTEGER                         :: Loc_hwmElemID, Glob_hwmElemID, OppSideID
REAL                            :: hwmElem_NodeCoords(3,0:NGeo,0:NGeo,0:NGeo), OppSideNodeCoords(3,4)
INTEGER                         :: TotalNWMLESSides, GlobalOppSideID, InnerOppSideID
REAL                            :: OppSideEdge(3,2), OppSideNormal(3), hwmDistVec(3), TolVec(3)
INTEGER                         :: i,j,k,p,q,r

INTEGER, ALLOCATABLE            :: WMLESToBCSide_tmp(:), TauW_Proc_tmp(:,:,:,:)
REAL, ALLOCATABLE               :: OthersPointInfo(:,:,:) ! indices-- 1-3: hwm_Coords, 4-6: TangVec1, 7: Glob_hwmElemID
INTEGER                         :: WallStressCount_local(0:nProcessors-1), WallStressCount(0:nProcessors-1)
INTEGER                         :: FirstElemInd, LastElemInd, FirstSideInd, LastSideInd
REAL, ALLOCATABLE               :: PointInfo(:,:)
REAL                            :: h_wm_Coords(3), DistanceVect(3), Distance
LOGICAL                         :: FoundhwmElem, FoundhwmPoint
INTEGER                         :: nPoints_MINE_tmp(2), nPoints_YOUR_tmp(2) ! indices -- 1: my rank, 2: how many points u calculate for me
INTEGER                         :: nPoints_MINE_tmp2(0:nProcessors-1), Proc_RecvTauW_tmp(nProcessors), Proc_SendTauW_tmp(nProcessors)

LOGICAL, ALLOCATABLE            :: TauW_MINE_IsFace(:,:), TauW_MINE_IsInterior(:,:), TauW_MINE_IsInterpolation(:,:)

#if USE_MPI
INTEGER, ALLOCATABLE            :: OthersSideInfo(:,:), OthersElemInfo(:)
INTEGER                         :: iStat(MPI_STATUS_SIZE), CalcInfoRequests(0:nProcessors-1), PointInfoRequests(0:nProcessors-1)
INTEGER                         :: hwmRank
#endif
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.(.NOT.MeshInitIsDone).OR.WMLESInitDone) THEN
  CALL CollectiveStop(__STAMP__,&
    'Wall-Modeled LES not ready to be called or already called.')
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT Wall-Modeled LES...'

WallModel = GETINTFROMSTR('WMLES')

! If Schumann's model is selected, check if we are in a channel flow
IF ((WallModel .EQ. WMLES_SCHUMANN) .AND. (.NOT.STRICMP(Testcase, Channel))) THEN
    CALL CollectiveStop(__STAMP__,&
        "Schumann's wall model can only be applied to the Channel flow testcase.")
END IF


h_wm = GETREAL('h_wm')
IF (H_WM .LE. 0) &
    CALL CollectiveStop(__STAMP__,&
        'h_wm distance from the wall must be positive.')

delta = GETREAL('delta')

abs_h_wm = h_wm*delta

! =-=-=-=--=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! We now set up the h_wm locations, determine in which element they lie,
! and who is the MPI proc. responsible for the element. 
! This is done by reading the hdf5 mesh file to avoid complications by 
! communications AND DEADLOCKS!! when determining such information.
!
! Also, we set up the information for communication between the proc. responsible
! for the calculation (owns h_wm element) and the proc. responsible for BC
! imposition (owns WMLESSide)
!
! Finally, we check whether interpolation is needed or if the h_wm point may be approximated
! by an interior/face point in the element. If interpolation is needed, all the info
! for this interpolation is also set up (i.e. mapping matrices, basis functions, etc)
! =-=-=-=--=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

ALLOCATE(BCSideToWMLES(nBCSides))
ALLOCATE(WMLESToBCSide_tmp(nBCSides))
ALLOCATE(TauW_Proc_tmp(2,0:PP_N, 0:PP_NZ, nBCSides))
ALLOCATE(OthersPointInfo(7,nBCSides*(PP_N+1)*(PP_NZ+1),0:nProcessors-1))

nWMLESSides = 0
BCSideToWMLES = 0
WMLESToBCSide_tmp = 0
TauW_Proc_tmp = -1
OthersPointInfo(7,:,:) = -1
WallStressCount_local = 0
WallStressCount = 0
WMLES_Tol = 1.E6

DO iSide=1,nBCSides
    IF (BoundaryType(BC(iSide),BC_TYPE) .EQ. 5) THEN ! WMLES side

        nWMLESSides = nWMLESSides + 1
        BCSideToWMLES(iSide) = nWMLESSides
        WMLESToBCSide_tmp(nWMLESSides) = iSide

        ! Since not every MPI proc gets into this IF, we cannot open the file 
        ! collectively, unfortunately. Hence, single = TRUE so that every proc.
        ! within this IF opens the file on its own.
        ! If this becomes too heavy a burden, we may pre-check which procs have
        ! WMLES Sides and create a communicator for them to open the file collectively.
        CALL OpenDataFile(MeshFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)

        DO q=0, PP_NZ; DO p=0, PP_N
            ! Calculate h_wm vector coords
            h_wm_Coords = -1.*abs_h_wm*NormVec(:,p,q,0,iSide) + Face_xGP(:,p,q,0,iSide)

            ! Set vars for first iteration on finding h_wm element
            Loc_hwmElemID = SideToElem(S2E_ELEM_ID, iSide) ! Wall elem is always master of BC side
            Glob_hwmElemID = Loc_hwmElemID+offsetElem ! +offset to put it in mesh-global numbers
            GlobalOppSideID = SideToGlobalSide(iSide)

            
            
            FoundhwmElem = .FALSE.

            DO WHILE (.NOT. FoundhwmElem)
                !> Given global element and global "lower side", get the "upper side"
                
                ! Check whose element is this, so that we read correct info from file
                hwmRank = ELEMIPROC(Glob_hwmElemID)
                SDEALLOCATE(OthersElemInfo)
                ! FirstElemInd = offsetElemMPI(hwmRank)+1
                ! LastElemInd = offsetElemMPI(hwmRank+1)
                ! ALLOCATE(OthersElemInfo(6,FirstElemInd:LastElemInd))
                ALLOCATE(OthersElemInfo(6))

                IF (.NOT. (hwmRank.EQ.myRank)) THEN
                    ! Read the correct portion of mesh file
                    ! CALL ReadArray('ElemInfo',2,(/6,(LastElemInd-FirstElemInd+1)/),offsetElemMPI(hwmRank),2,IntArray=OthersElemInfo)
                    CALL ReadArray('ElemInfo',2,(/6,1/),Glob_hwmElemID-1,2,IntArray=OthersElemInfo)
                ELSE
                    ! Use the ElemInfo already in memory from mesh read-in
                    ! OthersElemInfo = ElemInfo
                    OthersElemInfo = ElemInfo(:,Glob_hwmElemID)
                END IF

                SDEALLOCATE(OthersSideInfo)
                ! offsetSideID=OthersElemInfo(3,FirstElemInd) ! hdf5 array starts at 0-> -1
                ! nSideIDs    =OthersElemInfo(4,LastElemInd)-OthersElemInfo(3,FirstElemInd)
                offsetSideID=OthersElemInfo(3) ! hdf5 array starts at 0-> -1
                nSideIDs    =OthersElemInfo(4)-OthersElemInfo(3)
                    
                FirstSideInd=offsetSideID+1
                LastSideInd =offsetSideID+nSideIDs
                ALLOCATE(OthersSideInfo(5,FirstSideInd:LastSideInd))

                IF (.NOT. (hwmRank.EQ.myRank)) THEN
                    ! Read the correct portion of mesh file
                    CALL ReadArray('SideInfo',2,(/5,(nSideIDs)/),offsetSideID,2,IntArray=OthersSideInfo)
                ELSE
                    ! Use the SideInfo already in memory from mesh read-in
                    ! OthersSideInfo = SideInfo
                    OthersSideInfo = SideInfo(:,FirstSideInd:LastSideInd)
                END IF

                ! Scan the list of sides looking for a side with the Ind = GlobalOppSideID
                ! If found, check if this side is inner to the GlobalhwmElement. (deprecated)
                ! If so, we found our guy. Otherwise, we found the inner side of a neighbor. Keep looking. (deprecated)
                DO i=FirstSideInd, LastSideInd
                    IF (ABS(OthersSideInfo(2,i)) .EQ. GlobalOppSideID) THEN
                        ! Ok, we found a candidate.
                        ! IF ((i.GT.OthersElemInfo(3,Glob_hwmElemID)) .AND. (i.LE.OthersElemInfo(4,Glob_hwmElemID))) THEN
                            ! Yes, this is our guy.
                            InnerOppSideID = i
                            EXIT
                        ! END IF
                    END IF
                END DO

                ! NodeCoords array was deallocated after mesh read-in, so we read it again regardless
                CALL ReadArray('NodeCoords',2,(/3,(NGeo+1)**3/),(Glob_hwmElemID-1)*(NGeo+1)**3,2,RealArray=hwmElem_NodeCoords)
                ! Now get opposite side of WMLES Side or of the previous side being analyzed and its coords
                ! LocSide = InnerOppSideID - OthersElemInfo(3,Glob_hwmElemID)
                LocSide = InnerOppSideID - OthersElemInfo(3)
                SELECT CASE(LocSide)
                CASE (ETA_MINUS)
                    ! InnerOppSideID = OthersElemInfo(3,Glob_hwmElemID) + ETA_PLUS
                    InnerOppSideID = OthersElemInfo(3) + ETA_PLUS
                    OppSideNodeCoords(:,1) = hwmElem_NodeCoords(:,0,NGeo,0)
                    OppSideNodeCoords(:,2) = hwmElem_NodeCoords(:,0,NGeo,NGeo)
                    OppSideNodeCoords(:,3) = hwmElem_NodeCoords(:,NGeo,NGeo,0)
                    OppSideNodeCoords(:,4) = hwmElem_NodeCoords(:,NGeo,NGeo,NGeo)
                    ! Compute the vector to check node approximation tolerance
                    TolVec = hwmElem_NodeCoords(:,0,NGeo,0) - hwmElem_NodeCoords(:,0,0,0)
                CASE (ETA_PLUS)
                    ! InnerOppSideID = OthersElemInfo(3,Glob_hwmElemID) + ETA_MINUS
                    InnerOppSideID = OthersElemInfo(3) + ETA_MINUS
                    OppSideNodeCoords(:,1) = hwmElem_NodeCoords(:,0,0,0)
                    OppSideNodeCoords(:,2) = hwmElem_NodeCoords(:,NGeo,0,0)
                    OppSideNodeCoords(:,3) = hwmElem_NodeCoords(:,0,0,NGeo)
                    OppSideNodeCoords(:,4) = hwmElem_NodeCoords(:,NGeo,0,NGeo)
                    ! Compute the vector to check node approximation tolerance
                    TolVec = hwmElem_NodeCoords(:,0,NGeo,0) - hwmElem_NodeCoords(:,0,0,0)
                CASE (XI_MINUS)
                    ! InnerOppSideID = OthersElemInfo(3,Glob_hwmElemID) + XI_PLUS
                    InnerOppSideID = OthersElemInfo(3) + XI_PLUS
                    OppSideNodeCoords(:,1) = hwmElem_NodeCoords(:,NGeo,0,0)
                    OppSideNodeCoords(:,2) = hwmElem_NodeCoords(:,NGeo,NGeo,0)
                    OppSideNodeCoords(:,3) = hwmElem_NodeCoords(:,NGeo,0,NGeo)
                    OppSideNodeCoords(:,4) = hwmElem_NodeCoords(:,NGeo,NGeo,NGeo)
                    ! Compute the vector to check node approximation tolerance
                    TolVec = hwmElem_NodeCoords(:,NGeo,0,0) - hwmElem_NodeCoords(:,0,0,0)
                CASE (XI_PLUS)
                    ! InnerOppSideID = OthersElemInfo(3,Glob_hwmElemID) + XI_MINUS
                    InnerOppSideID = OthersElemInfo(3) + XI_MINUS
                    OppSideNodeCoords(:,1) = hwmElem_NodeCoords(:,0,0,0)
                    OppSideNodeCoords(:,2) = hwmElem_NodeCoords(:,0,0,NGeo)
                    OppSideNodeCoords(:,3) = hwmElem_NodeCoords(:,0,NGeo,0)
                    OppSideNodeCoords(:,4) = hwmElem_NodeCoords(:,0,NGeo,NGeo)
                    ! Compute the vector to check node approximation tolerance
                    TolVec = hwmElem_NodeCoords(:,NGeo,0,0) - hwmElem_NodeCoords(:,0,0,0)
                CASE (ZETA_MINUS)
                    ! InnerOppSideID = OthersElemInfo(3,Glob_hwmElemID) + ZETA_PLUS
                    InnerOppSideID = OthersElemInfo(3) + ZETA_PLUS
                    OppSideNodeCoords(:,1) = hwmElem_NodeCoords(:,0,0,NGeo)
                    OppSideNodeCoords(:,2) = hwmElem_NodeCoords(:,NGeo,0,NGeo)
                    OppSideNodeCoords(:,3) = hwmElem_NodeCoords(:,0,NGeo,NGeo)
                    OppSideNodeCoords(:,4) = hwmElem_NodeCoords(:,NGeo,NGeo,NGeo)
                    ! Compute the vector to check node approximation tolerance
                    TolVec = hwmElem_NodeCoords(:,0,0,NGeo) - hwmElem_NodeCoords(:,0,0,0)
                CASE (ZETA_PLUS)
                    ! InnerOppSideID = OthersElemInfo(3,Glob_hwmElemID) + ZETA_MINUS
                    InnerOppSideID = OthersElemInfo(3) + ZETA_MINUS
                    OppSideNodeCoords(:,1) = hwmElem_NodeCoords(:,0,0,0)
                    OppSideNodeCoords(:,2) = hwmElem_NodeCoords(:,0,NGeo,0)
                    OppSideNodeCoords(:,3) = hwmElem_NodeCoords(:,NGeo,0,0)
                    OppSideNodeCoords(:,4) = hwmElem_NodeCoords(:,NGeo,NGeo,0)
                    ! Compute the vector to check node approximation tolerance
                    TolVec = hwmElem_NodeCoords(:,0,0,NGeo) - hwmElem_NodeCoords(:,0,0,0)
                END SELECT

                GlobalOppSideID = ABS(OthersSideInfo(2,InnerOppSideID))

                ! Set node approximation tolerance (sphere radius) to 10% of the distance
                ! between interpolation nodes based on equidistant distribution
                RealTmp = 0
                DO k=1,3
                    RealTmp = RealTmp + TolVec(k)**2
                END DO
                RealTmp = 0.1*SQRT(RealTmp)/DBLE(PP_N)
                IF (WMLES_Tol .GT. RealTmp) WMLES_Tol = RealTmp
                
                ! Check if h_wmVec lies above or below this opposite side

                !> Calculate the normal vector of the opposite side (assume 4 coplanar vertices)
                !> -- Since the 4 vertices are assumed coplanar, take arbitrarily 3 and compute normal
                ! In this case, it is not strictly arbitrary 3 points. Which ones we take is set
                ! within the case selection above. These guarantee an outward normal
                OppSideEdge(:,1) = OppSideNodeCoords(:,2) - OppSideNodeCoords(:,1)
                OppSideEdge(:,2) = OppSideNodeCoords(:,3) - OppSideNodeCoords(:,1)

                OppSideNormal(1) = OppSideEdge(2,1)*OppSideEdge(3,2) - OppSideEdge(3,1)*OppSideEdge(2,2)
                OppSideNormal(2) = OppSideEdge(3,1)*OppSideEdge(1,2) - OppSideEdge(1,1)*OppSideEdge(3,2)
                OppSideNormal(3) = OppSideEdge(1,1)*OppSideEdge(2,2) - OppSideEdge(2,1)*OppSideEdge(1,2)

                hwmDistVec(:) = h_wm_Coords(:) - OppSideNodeCoords(:,1)

                RealTmp = DOT_PRODUCT(OppSideNormal, hwmDistVec)

                IF (ABS(RealTmp) .LE. WMLES_Tol) THEN ! Lies on face (angle between vecs ~ 90 degrees)
                    ! Mark as "next element" so info is retrieved from the approximation further away from the wall
                    Glob_hwmElemID = ABS(OthersSideInfo(3,InnerOppSideID))
                    FoundhwmElem = .TRUE.
                ELSE
                    IF (RealTmp .LE. 0) THEN ! Lies within this element
                        ! Mark as within this element.
                        ! Then set the loop var to false so that we can go to the next step, of finding
                        ! which MPI proc is the owner of this element.
                        FoundhwmElem = .TRUE.

                    ELSE ! Is not within this element
                        ! Mark as probably being within the next element -- needs check after looping
                        Glob_hwmElemID = ABS(OthersSideInfo(3,InnerOppSideID))
                    END IF
                END IF

                IF (Glob_hwmElemID .EQ. 0) THEN
                    ! We reached the opposite boundary!! Something is wrong...
                    CALL CollectiveStop(__STAMP__, 'Error defining h_wm element -- Possibly reached the opposite boundary!')
                END IF

            END DO ! do while

            ! Found h_wm element for this p,q face point
            LOGWRITE(*,'(2(I4,2X),(I15,2X))') p, q, Glob_hwmElemID, ELEMIPROC(Glob_hwmElemID)

            hwmRank = ELEMIPROC(Glob_hwmElemID)
            WallStressCount_local(hwmRank) = WallStressCount_local(hwmRank) + 1
            TauW_Proc_tmp(1,p,q,nWMLESSides) = hwmRank
            TauW_Proc_tmp(2,p,q,nWMLESSides) = WallStressCount_local(hwmRank)
            OthersPointInfo(1:3,WallStressCount_local(hwmRank),hwmRank) = h_wm_Coords

            ! TODO:
            ! The normal and tangent vectors are oriented in outward-facing coordinates
            ! with the tangent vector in z pointing to z+. Hence, for a lower wall, the 
            ! tangent vector in the x direction is poiting to x-, so we need the "-" sign.
            ! However, for upper walls, for example, the tangent vector in x direction will 
            ! be pointing to x+, so we do not need the "-" sign.
            ! Thus, a "check" must be implemented here so that we pass the tang. vec. pointing
            ! to the correct direction, that is, streamwise.

            ! One solution would be to take the dot product of this vector and the initial
            ! velocity so that we know if it is streamwise or not, and change sign accordingly.
            OthersPointInfo(4:6,WallStressCount_local(hwmRank),hwmRank) = -TangVec2(:,p,q,0,iSide)
            OthersPointInfo(7,WallStressCount_local(hwmRank),hwmRank) = Glob_hwmElemID

        END DO; END DO ! p, q

        CALL CloseDataFile()
    END IF ! 
END DO
SDEALLOCATE(OthersSideInfo)
SDEALLOCATE(OthersElemInfo)

CALL MPI_REDUCE(nWMLESSides, TotalNWMLESSides, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_FLEXI, iError)
IF (MPIRoot .AND. (TotalNWMLESSides .EQ. 0)) THEN
    CALL PrintWarning('A Wall Model is set up, but no BoundaryConditions'// &
                                ' of WMLES Type has been found on the Mesh/Parameter File!')
END IF 

#if USE_MPI
!> Reduce and broadcast the minimun node approximation tolerance found in the previous steps
IF (MPIRoot) THEN
    CALL MPI_REDUCE(MPI_IN_PLACE, WMLES_Tol, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_FLEXI, iError)
ELSE
    CALL MPI_REDUCE(WMLES_Tol, WMLES_Tol, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_FLEXI, iError)
END IF
CALL MPI_BCAST(WMLES_Tol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FLEXI, iError)

!> Sum up all tau_w counts and broadcast this info so that
!> every MPI proc knows it (and not only those who actually computed it)
CALL MPI_REDUCE(WallStressCount_local, WallStressCount, nProcessors, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_FLEXI, iError)
CALL MPI_BCAST(WallStressCount, nProcessors, MPI_INTEGER, 0, MPI_COMM_FLEXI, iError)
#endif

!> Mount data and send it to other procs (non-blocking)
CalcInfoRequests = MPI_REQUEST_NULL
PointInfoRequests = MPI_REQUEST_NULL
nProcs_RecvTauW = 0
IF (nWMLESSides .NE. 0) THEN
    ALLOCATE(Proc_RecvTauW_Inv(0:nProcessors-1))
    Proc_RecvTauW_Inv = -1
    Proc_RecvTauW_Inv(myRank) = 0
    DO i=0,nProcessors-1
        IF (i .EQ. myRank) CYCLE
        IF (WallStressCount_local(i) .GT. 0) THEN
            nProcs_RecvTauW = nProcs_RecvTauW + 1
            Proc_RecvTauW_tmp(nProcs_RecvTauW) = i
            Proc_RecvTauW_Inv(i) = nProcs_RecvTauW
            ! First, send a message for each proc. that need to calculate tau_w for you
            ! containing your rank and how many points they will calculate for you (tag: 0)
            nPoints_YOUR_tmp(1) = myRank
            nPoints_YOUR_tmp(2) = WallStressCount_local(i)
            CALL MPI_ISEND(nPoints_YOUR_tmp,2,MPI_INTEGER,i,0,MPI_COMM_FLEXI,CalcInfoRequests(i),iError)
            !CALL MPI_Request_Free(sendRequest, iError)
            ! Then, send the message itself (tag: 1)
            CALL MPI_ISEND(OthersPointInfo(:,:,i),7*WallStressCount_local(i),MPI_DOUBLE_PRECISION,i,1,MPI_COMM_FLEXI,PointInfoRequests(i),iError)
            !CALL MPI_Request_Free(sendRequest, iError)
        END IF
    END DO
END IF

!> Check if we need to receive any data from other procs
nProcs_SendTauW = 0
nPoints_MINE_tmp2 = 0
IF (WallStressCount(myRank) .NE. 0) THEN
    IF (WallStressCount_local(myRank) .NE. WallStressCount(myRank)) THEN ! I calculate info for other procs
        MineCnt = WallStressCount_local(myRank)
        DO WHILE (MineCnt .LT. WallStressCount(myRank))
            CALL MPI_RECV(nPoints_MINE_tmp,2,MPI_INTEGER,MPI_ANY_SOURCE,0,MPI_COMM_FLEXI,iStat,iError)
            nPoints_MINE_tmp2(nPoints_MINE_tmp(1)) = nPoints_MINE_tmp(2)
            MineCnt = MineCnt + nPoints_MINE_tmp(2)
            nProcs_SendTauW = nProcs_SendTauW + 1
            Proc_SendTauW_tmp(nProcs_SendTauW) = nPoints_MINE_tmp(1)
        END DO

    END IF
END IF

!> Free some space: allocate and populate permanent information, delete temporaries
ALLOCATE(WMLESToBCSide(nWMLESSides))
ALLOCATE(TauW_Proc(2,0:PP_N, 0:PP_NZ, nWMLESSides))
ALLOCATE(Proc_RecvTauW(nProcs_RecvTauW))
ALLOCATE(Proc_SendTauW(nProcs_SendTauW))
ALLOCATE(nTauW_MINE(0:nProcs_SendTauW)) ! index 0 are my LOCAL points
ALLOCATE(nTauW_YOURS(nProcs_RecvTauW))

DO i=1,nWMLESSides
    WMLESToBCSide(i) = WMLESToBCSide_tmp(i)
    DO q=0,PP_NZ; DO p=0,PP_N
        TauW_Proc(:,p,q,i) = TauW_Proc_tmp(:,p,q,i)
    END DO; END DO ! p,q
END DO
DO i=1,nProcs_RecvTauW
    Proc_RecvTauW(i) = Proc_RecvTauW_tmp(i)
    nTauW_YOURS(i) = WallStressCount_local(Proc_RecvTauW(i))
END DO
DO i=1,nProcs_SendTauW
    Proc_SendTauW(i) = Proc_SendTauW_tmp(i)
    nTauW_MINE(i) = nPoints_MINE_tmp2(Proc_SendTauW(i))
END DO
nTauW_MINE(0) = WallStressCount_local(myRank)
DEALLOCATE(WMLESToBCSide_tmp)
DEALLOCATE(TauW_Proc_tmp)

!> Those procs who need to receive data, do it now.
!> Once received, check if the points may be approximated by any interpolation node,
!> or if we need to interpolate the solution
ALLOCATE(PointInfo(7,MAXVAL(nTauW_MINE)))
ALLOCATE(TauW_MINE_FacePoint(3,MAXVAL(nTauW_MINE),0:nProcs_SendTauW))
ALLOCATE(TauW_MINE_InteriorPoint(4,MAXVAL(nTauW_MINE),0:nProcs_SendTauW))
ALLOCATE(TauW_MINE_TangVec(3,MAXVAL(nTauW_MINE),0:nProcs_SendTauW))

IF (Logging) THEN ! Logging info only
ALLOCATE(TauW_MINE_IsFace(MAXVAL(nTauW_MINE),0:nProcs_SendTauW))
ALLOCATE(TauW_MINE_IsInterior(MAXVAL(nTauW_MINE),0:nProcs_SendTauW))
ALLOCATE(TauW_MINE_IsInterpolation(MAXVAL(nTauW_MINE),0:nProcs_SendTauW))
TauW_MINE_IsFace = .FALSE.
TauW_MINE_IsInterior = .FALSE.
TauW_MINE_IsInterpolation = .FALSE.
END IF

ALLOCATE(nTauW_MINE_FacePoint(0:nProcs_SendTauW))
ALLOCATE(nTauW_MINE_InteriorPoint(0:nProcs_SendTauW))
ALLOCATE(nTauW_MINE_Interpolate(0:nProcs_SendTauW))
ALLOCATE(FaceToLocalPoint(MAXVAL(nTauW_MINE)))
ALLOCATE(InteriorToLocalPoint(MAXVAL(nTauW_MINE)))
ALLOCATE(InterpToLocalPoint(MAXVAL(nTauW_MINE)))

nTauW_MINE_FacePoint = 0
nTauW_MINE_InteriorPoint = 0
nTauW_MINE_Interpolate = 0
FaceToLocalPoint = 0
InteriorToLocalPoint = 0
InterpToLocalPoint = 0



DO i=0,nProcs_SendTauW
    IF (i.NE.0) THEN 
        CALL MPI_RECV(PointInfo,7*nTauW_MINE(i),MPI_DOUBLE_PRECISION,Proc_SendTauW(i),1,MPI_COMM_FLEXI,iStat,iError)
    ELSE
        IF (nTauW_MINE(0).NE.0) PointInfo = OthersPointInfo(:,1:WallStressCount_local(myRank),myRank)
    END IF

    ! Store wall-tangent vector
    IF (nTauW_MINE(i).NE.0) TauW_MINE_TangVec(1:3,:,i) = PointInfo(4:6,:)
    
    DO j=1,nTauW_MINE(i)

        ! For each point to calculate tau_w for this proc, check whether it may be
        ! approximated by a point inside the element or if interpolation is needed.

        ! Global to local element ID
        Loc_hwmElemID = INT(PointInfo(7,j))-offsetElem

        FoundhwmPoint = .FALSE.
        ! Start from faces
        DO LocSide=1,6
            IF (FoundhwmPoint) EXIT

            SideID = ElemToSide(E2S_SIDE_ID, LocSide, Loc_hwmElemID)
            CALL GetOppositeSide(SideID, Loc_hwmElemID, OppSideID)
            DO q=0,PP_NZ
                IF (FoundhwmPoint) EXIT
                DO p=0,PP_N
                    IF (FoundhwmPoint) EXIT
                    DistanceVect = PointInfo(1:3,j) - Face_xGP(:,p,q,0,SideID)
                    Distance = 0
                    DO k=1,3
                        Distance = Distance + DistanceVect(k)**2
                    END DO
                    Distance = SQRT(Distance)
                    IF (Distance .LE. WMLES_Tol) THEN ! May be approximated by this point, in this face
                        IF (Logging) TauW_MINE_IsFace(j,i) = .TRUE.
                        nTauW_MINE_FacePoint(i) = nTauW_MINE_FacePoint(i) + 1
                        FaceToLocalPoint(nTauW_MINE_FacePoint(i)) = j
                        TauW_MINE_FacePoint(:,nTauW_MINE_FacePoint(i),i) = (/p,q,SideID/)
                        FoundhwmPoint = .TRUE.
                    END IF
                END DO ! p
            END DO ! q
        END DO ! LocSide

        ! If the point could not be approximated by a face point, 
        ! then check within the element
        DO r=0,PP_NZ
            IF (FoundhwmPoint) EXIT
            DO q=0,PP_N
                IF (FoundhwmPoint) EXIT
                DO p=0,PP_N
                    IF (FoundhwmPoint) EXIT
                    DistanceVect(:) = PointInfo(1:3,j) - Elem_xGP(:,p,q,r,Loc_hwmElemID)
                    Distance = 0
                    DO k=1,3
                        Distance = Distance + DistanceVect(k)**2
                    END DO
                    Distance = SQRT(Distance)
                    ! Tolerance has been set in the previous loop, we stick to it
                    IF (Distance .LE. WMLES_Tol) THEN ! May be approximated by this point
                        IF (Logging) TauW_MINE_IsInterior(j,i) = .TRUE.
                        nTauW_MINE_InteriorPoint(i) = nTauW_MINE_InteriorPoint(i) + 1
                        InteriorToLocalPoint(nTauW_MINE_InteriorPoint(i)) = j
                        TauW_MINE_InteriorPoint(:,nTauW_MINE_InteriorPoint(i),i) = (/p,q,r,Loc_hwmElemID/)
                        FoundhwmPoint = .TRUE.
                    END IF
                END DO ! p
            END DO ! q
        END DO ! r

        ! If the point may not be approximated neither by a face nor an interior point,
        ! then we must interpolate.
        IF (.NOT.FoundhwmPoint) THEN ! Interpolate
            WRITE(*,*) myRank, 'INTERPOLATION POINT', LocSide, Distance, WMLES_Tol, Loc_hwmElemID, SideID, OppSideID
            IF (Logging) TauW_MINE_IsInterpolation(j,i) = .TRUE.
            nTauW_MINE_Interpolate(i) = nTauW_MINE_Interpolate(i) + 1
            InterpToLocalPoint(nTauW_MINE_Interpolate(i)) = j
            ! TODO
            ! Set up the interpolation info
            ! - Use mappings to map h_wm coords to local coords and save this info

        END IF
    END DO ! j -- each point that we must calculate tau_w, for this MPI process

    ! Sanity check
    IF ((nTauW_MINE_FacePoint(i)+nTauW_MINE_Interpolate(i)+nTauW_MINE_InteriorPoint(i)) &
                .NE. nTauW_MINE(i)) CALL Abort(__STAMP__,"Dangling points?")

END DO ! i -- message from each process having a BC side which needs our calculation of tau_w

!> Set up variables for MPI and actual TauW values 
! IF (nWMLESSides .NE. 0) ALLOCATE(WMLES_TauW(2,0:PP_N,0:PP_NZ,nWMLESSides))
! IF (nProcs_SendTauW .NE. 0) ALLOCATE(TauW_MINE(2,MAXVAL(nTauW_MINE),0:nProcs_SendTauW))
! IF (nProcs_RecvTauW .NE. 0) THEN
!     ALLOCATE(TauW_YOURS(2,MAXVAL(WallStressCount_local),nProcs_RecvTauW))
!     ALLOCATE(WMLES_Requests(nProcs_RecvTauW))
! END IF
ALLOCATE(WMLES_TauW(2,0:PP_N,0:PP_NZ,nWMLESSides))
ALLOCATE(TauW_MINE(2,MAXVAL(nTauW_MINE),0:nProcs_SendTauW))
ALLOCATE(TauW_YOURS(2,MAXVAL(WallStressCount_local),0:nProcs_RecvTauW))
ALLOCATE(WMLES_RecvRequests(nProcs_RecvTauW))
ALLOCATE(WMLES_SendRequests(nProcs_SendTauW))

!> =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!> Dump Logging/Debugging Information
!> =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

LOGWRITE(*,*) "====================== START WMLES Info ===================="
LOGWRITE(*,'(A35,ES16.6)') 'Distance h_wm from the wall:', abs_h_wm
LOGWRITE(*,'(A35,I16)') 'Local WMLES BC Sides:', nWMLESSides
LOGWRITE(*,'(A35,I16)') 'Local WMLES BC Poins:', nWMLESSides*(PP_N+1)*(PP_NZ+1)
LOGWRITE(*,'(A35,I6)') 'Local Tau_w Count:', WallStressCount_local(myRank)
LOGWRITE(*,'(A35,I6)') 'Total Tau_w Count:', WallStressCount(myRank)
LOGWRITE(*,*) '---------------------------------------'
LOGWRITE(*,*) "======= Tau_w Points which are my responbility ======"
LOGWRITE(*,'(2(A15),A15,A20,A25)') 'Origin MPI Rank', 'Local Point', 'Face Point', 'Interior Point', 'Interpolation Point'
DO i=0,nProcs_SendTauW
    MineCnt=0
    DO k=i-1,0,-1
        MineCnt=MineCnt+nTauW_MINE(k)
    END DO
    DO j=1,nTauW_MINE(i)
        IF (i.EQ.0) THEN
            LOGWRITE(*,'(I10,I15,L15,L20,L25)') myRank, j+MineCnt, TauW_MINE_IsFace(j,i), TauW_MINE_IsInterior(j,i), TauW_MINE_IsInterpolation(j,i)
        ELSE
            LOGWRITE(*,'(I10,I15,L15,L20,L25)') Proc_SendTauW(i), j+MineCnt, TauW_MINE_IsFace(j,i), TauW_MINE_IsInterior(j,i), TauW_MINE_IsInterpolation(j,i)
        END IF
    END DO
END DO
LOGWRITE(*,*) '---------------------------------------'
LOGWRITE(*,*) "====================== END WMLES Info ===================="
IF (Logging) FLUSH(UNIT_logOut)

! Before allowing an MPI proc. to leave this subroutine (and thus possibly clean sending buffers/variables)
! we make sure that all non-blocking Send operations are completed (i.e. the send buffer may be modified)
CALL MPI_Waitall(nProcessors,CalcInfoRequests,MPI_STATUSES_IGNORE,iError)
CALL MPI_Waitall(nProcessors,PointInfoRequests,MPI_STATUSES_IGNORE,iError)

SDEALLOCATE(TauW_MINE_IsFace)
SDEALLOCATE(TauW_MINE_IsInterior)
SDEALLOCATE(TauW_MINE_IsInterpolation)
SDEALLOCATE(OthersPointInfo)
SDEALLOCATE(PointInfo)

!!!>> IMPORTANT INFO:
! During the calculation (ComputeWallStress), each proc will run a loop from
! i = 1,nProcs_SendTauW and within that from j = 1,nTauW_MINE(i) so that it calculates
! the points for MPI proc. 'i' one by one, and put it in an ordered array of 
! TauW_MINE_Send(1:2,nTauW_MINE(i)). Then, communicate this result to the proc. that imposes the BC,
! i.e., proc. 'i'.
! Then, proc 'i' receives this info and maps it to the actual WMLES array TauW(1:2,p,q,nWMLESSide)
! using the info on TauW_Proc array, which contains info on p, q, and WMLES Side, and the MPI proc.
! responsible for calculating wall stress for this point.
WMLESInitDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT Wall-Modeled LES DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitWMLES


!==================================================================================================================================
!> Get the ID of the side opposite of that indicated by SideID, in the same element
!==================================================================================================================================
SUBROUTINE GetOppositeSide(SideID, ElemID, OppositeSide)
! MODULES
USE MOD_Mesh_Vars,              ONLY: SideToElem, ElemToSide
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, INTENT(IN)                 :: ElemID, SideID
INTEGER, INTENT(OUT)                :: OppositeSide
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER                             :: LocSide, tmpElemID
!==================================================================================================================================

tmpElemID = SideToElem(S2E_ELEM_ID, SideID)
IF (tmpElemID .EQ. ElemID) THEN ! ElemID is master of this side
    LocSide = SideToElem(S2E_LOC_SIDE_ID, SideID)
ELSE ! ElemID is slave of this side
    LocSide = SideToElem(S2E_NB_LOC_SIDE_ID, SideID)
END IF

SELECT CASE(LocSide)
CASE (ETA_MINUS)
    OppositeSide = ElemToSide(E2S_SIDE_ID, ETA_PLUS, ElemID)
CASE (ETA_PLUS)
    OppositeSide = ElemToSide(E2S_SIDE_ID, ETA_MINUS, ElemID)
CASE (XI_MINUS)
    OppositeSide = ElemToSide(E2S_SIDE_ID, XI_PLUS, ElemID)
CASE (XI_PLUS)
    OppositeSide = ElemToSide(E2S_SIDE_ID, XI_MINUS, ElemID)
CASE (ZETA_MINUS)
    OppositeSide = ElemToSide(E2S_SIDE_ID, ZETA_PLUS, ElemID)
CASE (ZETA_PLUS)
    OppositeSide = ElemToSide(E2S_SIDE_ID, ZETA_MINUS, ElemID)
END SELECT

END SUBROUTINE

!==================================================================================================================================
!> Get the ID of element where the exchange location point, h_wm (p,q), lies.
!==================================================================================================================================
! SUBROUTINE FindHwmElement(p, q, SideID, hwmElemID)
! ! MODULES
! USE MOD_Globals
! USE MOD_Mesh_Vars,              !ONLY: NormVec, Face_xGP, SideToElem, nElems
! USE MOD_WMLES_Vars
! #if USE_MPI
! USE MOD_MPI_Vars
! #endif
! IMPLICIT NONE
! !----------------------------------------------------------------------------------------------------------------------------------
! ! LOCAL VARIABLES
! INTEGER, INTENT(IN)                 :: p, q
! INTEGER, INTENT(INOUT)              :: hwmElemID, SideID!, OppSideID
! !----------------------------------------------------------------------------------------------------------------------------------

!==================================================================================================================================


! END SUBROUTINE FindHwmElement

!==================================================================================================================================
!> Compute the wall stress, tau_w, at each point in a WMLES BC surface.
!==================================================================================================================================
SUBROUTINE ComputeWallStress()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_WMLES_Vars
USE MOD_Mesh_Vars
USE MOD_DG_Vars                     
USE MOD_EOS_Vars                    ,ONLY: mu0
IMPLICIT NONE

!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                             :: iSide, sProc, FPInd, IPInd, IntPInd
INTEGER                             :: p,q,r,SideID,ElemID
REAL                             :: u_tau, u_mean, utang
REAL, PARAMETER                  :: vKarman=0.41, B=5.5
!==================================================================================================================================
IF (.NOT.WMLESInitDone) THEN
    CALL CollectiveStop(__STAMP__,&
    'Cannot compute Wall Stress -- Wall-Modeled LES INIT is not done.')
END IF

! Start non-blockingly receiving WMLES info (for those who have anything to receive)
CALL StartReceiveTauWMPI()

! Compute tau_w for points of our responsibility
SELECT CASE(WallModel)

  CASE(WMLES_SCHUMANN)
    
    ! Schumann model
    ! Computation is done using instantaneous information, projected onto the wall-tangent direction
    
    ! The channel flow testcase assumes <tau_w> = 1.0

    DO sProc=0,nProcs_SendTauW
        ! Calculate tau_w for each h_wm that is approximated as a face node
        DO FPInd=1,nTauW_MINE_FacePoint(sProc)
            p = TauW_MINE_FacePoint(1,FPInd,sProc)
            q = TauW_MINE_FacePoint(2,FPInd,sProc)
            SideID = TauW_MINE_FacePoint(3,FPInd,sProc)

            ! Face_xGP is only populated for master sides, so that only master sides have approximation nodes (check InitWMLES above)
            ! hence, use of UPrim_master is guaranteed to be correct here
            utang = DOT_PRODUCT(UPrim_master(2:4,p,q,SideID),TauW_MINE_TangVec(:,FaceToLocalPoint(FPInd),sProc))
            ! utang = UPrim_master(2,p,q,SideID)

            u_tau = SQRT(1.0/UPrim_master(1,p,q,SideID)) 
            u_mean = u_tau*( (1./vKarman) * LOG((abs_h_wm*u_tau)/mu0) + B )
            TauW_MINE(1,FaceToLocalPoint(FPInd),sProc) = utang/u_mean * 1.0 ! <tau_w> = 1.0
            TauW_MINE(2,FaceToLocalPoint(FPInd),sProc) = (2.0*mu0)*(UPrim_master(3,p,q,SideID)/abs_h_wm)
        END DO

        ! Calculate tau_w for each h_wm that is approximated as an interior node
        DO IPInd=1,nTauW_MINE_InteriorPoint(sProc)
            p = TauW_MINE_InteriorPoint(1,IPInd,sProc)
            q = TauW_MINE_InteriorPoint(2,IPInd,sProc)
            r = TauW_MINE_InteriorPoint(3,IPInd,sProc)
            ElemID = TauW_MINE_InteriorPoint(4,IPInd,sProc)

            utang = DOT_PRODUCT(U(2:4,p,q,r,ElemID),TauW_MINE_TangVec(:,InteriorToLocalPoint(IPInd),sProc))
            ! utang = U(2,p,q,r,ElemID)

            u_tau = SQRT(1.0/U(1,p,q,r,ElemID)) 
            u_mean = u_tau*( (1./vKarman) * LOG((abs_h_wm*u_tau)/mu0) + B )
            TauW_MINE(1,InteriorToLocalPoint(IPInd),sProc) = utang/u_mean * 1.0 ! <tau_w> = 1.0
            TauW_MINE(2,InteriorToLocalPoint(IPInd),sProc) = (2.0*mu0)*(U(3,p,q,r,ElemID)/abs_h_wm)
        END DO

        ! Calculate tau_w for each h_wm that must be interpolated
        DO IntPInd=1,nTauW_MINE_Interpolate(sProc)

        END DO
    END DO


    ! Send computed tau_w for processess responsible for the BC imposition on such points


    ! Slave top sides first
    ! LOGWRITE(*,*) '================= Wall Stress Calc ==================='
    ! DO iSide=1,nSlaveSides
    !   LOGWRITE(*,*) '------- Slave Side', iSide, '--------'
    !   LOGWRITE(*,*) 'SideID', SlaveToTSide(iSide)
    !   LOGWRITE(*,*) 'Wall SideID', WMLES_SideInv(SlaveToWMLESSide(iSide))
    !   LOGWRITE(*,*) 'WMLES SideID', SlaveToWMLESSide(iSide)
    !   LOGWRITE(*,'(2(A4,2X),9(A15,2X))') 'p', 'q', 'offdist', 'Prim(1)', 'u_tau', 'mu0', 'u_mean', 'Prim(2)', 'Prim(3)', 'tau_xy', 'tau_yz'
    !   DO q=0,PP_NZ; DO p=0,PP_N
    !     ! This is constant, we can calculate it in pre-processing stage (one more vector, 1:nWMLESSides)
    !     offdist = ABS(Face_xGP(2,p,q,0,SlaveToTSide(iSide)) - Face_xGP(2,p,q,0, WMLES_SideInv(SlaveToWMLESSide(iSide)) ))
    !     u_tau = SQRT(1.0/UPrim_slave(1,p,q, SlaveToTSide(iSide)) )
    !     u_mean = u_tau*( (1./vKarman) * LOG((offdist*u_tau)/mu0) + B )
    !     WMLES_Tauw(1,p,q, SlaveToWMLESSide(iSide) ) = (UPrim_slave(2,p,q, SlaveToTSide(iSide))/u_mean) * 1.0
    !     WMLES_Tauw(2,p,q, SlaveToWMLESSide(iSide) ) = (2.0*mu0)*(UPrim_slave(3,p,q, SlaveToTSide(iSide))/offdist)

    !     LOGWRITE(*,'(2(I4,2X),9(E15.8,2X))') p, q, offdist, UPrim_slave(1,p,q, SlaveToTSide(iSide)), u_tau, mu0, u_mean, UPrim_slave(2,p,q, SlaveToTSide(iSide)), UPrim_slave(3,p,q, SlaveToTSide(iSide)), WMLES_Tauw(1,p,q,SlaveToWMLESSide(iSide)), WMLES_Tauw(2,p,q,SlaveToWMLESSide(iSide))

    !   END DO; END DO ! p,q
    ! END DO
    ! ! Master top sides next
    ! DO iSide=1,nMasterSides
    !   LOGWRITE(*,*) '------- Master Side', iSide, '--------'
    !   LOGWRITE(*,*) 'SideID', MasterToTSide(iSide)
    !   LOGWRITE(*,*) 'Wall SideID', WMLES_SideInv(MasterToWMLESSide(iSide))
    !   LOGWRITE(*,*) 'WMLES SideID', MasterToWMLESSide(iSide)
    !   LOGWRITE(*,'(2(A4,2X),9(A15,2X))') 'p', 'q', 'offdist', 'Prim(1)', 'u_tau', 'mu0', 'u_mean', 'Prim(2)', 'Prim(3)', 'tau_xy', 'tau_yz'
    !   DO q=0,PP_NZ; DO p=0,PP_N
    !     offdist = ABS(Face_xGP(2,p,q,0,MasterToTSide(iSide)) - Face_xGP(2,p,q,0, WMLES_SideInv(MasterToWMLESSide(iSide)) ))
    !     u_tau = SQRT(1.0/UPrim_master(1,p,q, MasterToTSide(iSide)) )
    !     u_mean = u_tau*( (1./vKarman) * LOG((offdist*u_tau)/mu0) + B )
    !     WMLES_Tauw(1,p,q, MasterToWMLESSide(iSide) ) = (UPrim_master(2,p,q, MasterToTSide(iSide))/u_mean) * 1.0
    !     WMLES_Tauw(2,p,q, MasterToWMLESSide(iSide) ) = (2.0*mu0)*(UPrim_master(3,p,q, MasterToTSide(iSide))/offdist)

    !     LOGWRITE(*,'(2(I4,2X),9(E15.8,2X))') p, q, offdist, UPrim_master(1,p,q, MasterToTSide(iSide)), u_tau, mu0, u_mean, UPrim_master(2,p,q, MasterToTSide(iSide)), UPrim_master(3,p,q, MasterToTSide(iSide)), WMLES_Tauw(1,p,q,MasterToWMLESSide(iSide)), WMLES_Tauw(2,p,q,MasterToWMLESSide(iSide))

    !   END DO; END DO ! p,q
    ! END DO
    ! LOGWRITE(*,*) '================= End of Wall Stress Calc ==================='
    ! IF (Logging) FLUSH(UNIT_logOut)

  CASE DEFAULT
    CALL abort(__STAMP__,&
         'Unknown definition of Wall Model.')

END SELECT

! Start non-blockingly sending WMLES info (for those who have anything to send)
CALL StartSendTauWMPI()

! Finish receiving WMLES info (now this is a blocking op.)
CALL FinishReceiveTauWMPI()


! Set up WMLES_TauW vector with all necessary, computed values of wall shear stress
DO iSide=1,nWMLESSides
    DO q=0,PP_NZ; DO p=0,PP_N
        WMLES_TauW(:,p,q,iSide) = TauW_YOURS(:,TauW_Proc(2,p,q,iSide),Proc_RecvTauW_Inv(TauW_Proc(1,p,q,iSide)))
    END DO; END DO ! p, q
END DO


END SUBROUTINE ComputeWallStress


!==================================================================================================================================
!> Evaluate the viscous fluxes through the Boundary side of WMLES.
!> This is done simply setting the appropriate components from the previously calculated wall stress, WMLES_Tauw
!==================================================================================================================================
SUBROUTINE EvalDiffFlux3D_WMLES(Nloc, UPrim_master, Tauw, f, g, h)
! MODULES
USE MOD_Globals
IMPLICIT NONE

!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, INTENT(IN)                    :: Nloc
REAL, INTENT(IN)                       :: Tauw(1:2, 0:Nloc, 0:ZDIM(Nloc))
REAL, INTENT(IN)                       :: UPrim_master( PP_nVarPrim,0:Nloc,0:ZDIM(Nloc))
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: f,g,h
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER                                 :: p,q
REAL                                    :: tau_xy, tau_yz
!==================================================================================================================================

LOGWRITE(*,*) '-------- DIFFUSIVE FLUX (VISCOUS) ------------'
LOGWRITE(*,'(2(A4,2X),15(A15,2X))') 'p', 'q', 'f(1)', 'f(2)', 'f(3)', 'f(4)', 'f(5)',&
            'g(1)', 'g(2)', 'g(3)', 'g(4)', 'g(5)', 'h(1)', 'h(2)', 'h(3)', 'h(4)', 'h(5)'

DO q=0,ZDIM(Nloc); DO p=0,Nloc

  tau_xy = Tauw(1,p,q)
  tau_yz = Tauw(2,p,q)

  ! Assuming adiabatic wall, i.e., dT_dx = 0.
  f(1,p,q) = 0.
  f(2,p,q) = 0. ! tau_xx = 0
  f(3,p,q) = -tau_xy
  f(4,p,q) = 0. ! tau_xz = 0
  f(5,p,q) = -(UPrim_master(3,p,q)*tau_xy) ! -(u*tau_xx + v*tau_xy + w*tau_xz - kdTdx)

  g(1,p,q) = 0.
  g(2,p,q) = -tau_xy ! tau_yx = tau_xy
  g(3,p,q) = 0. ! tau_yy = 0
  g(4,p,q) = -tau_yz ! tau_yz = tau_zy
  g(5,p,q) = -(UPrim_master(2,p,q)*tau_xy + UPrim_master(4,p,q)*tau_yz) ! -(u*tau_yx + v*tau_yy + w*tau_yz - kdTdx)

  h(1,p,q) = 0.
  h(2,p,q) = 0. ! tau_zx = 0
  h(3,p,q) = -tau_yz ! tau_zy = tau_yz
  h(4,p,q) = 0. ! tau_zz = 0
  h(5,p,q) = -(UPrim_master(3,p,q)*tau_yz) ! -(u*tau_zx + v*tau_zy + w*tau_zz - kdTdx)


  LOGWRITE(*,'(2(I4,2X),15(E15.8,2X))') p, q, f(1,p,q), f(2,p,q), f(3,p,q), f(4,p,q), f(5,p,q),&
     g(1,p,q), g(2,p,q), g(3,p,q), g(4,p,q), g(5,p,q), h(1,p,q), h(2,p,q), h(3,p,q), h(4,p,q), h(5,p,q)

END DO; END DO ! p,q

LOGWRITE(*,'(X)')

END SUBROUTINE EvalDiffFlux3D_WMLES

!==================================================================================================================================
!> Subroutine that controls the receive operations for the Tau_W to be exchanged between processors.
!==================================================================================================================================
SUBROUTINE StartReceiveTauWMPI()
! MODULES
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_WMLES_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
! Every VAR we'll be using is declared in MOD_WMLES_Vars, loaded here.
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iProc
!==================================================================================================================================
DO iProc=1,nProcs_RecvTauW
    CALL MPI_IRECV(TauW_YOURS(:,:,iProc),2*nTauW_YOURS(iProc), MPI_DOUBLE_PRECISION, &
                    Proc_RecvTauW(iProc), 0, MPI_COMM_FLEXI, WMLES_RecvRequests(iProc),iError)
END DO

END SUBROUTINE StartReceiveTauWMPI


!==================================================================================================================================
!> Subroutine that controls the send operations for the Tau_W to be exchanged between processors.
!==================================================================================================================================
SUBROUTINE StartSendTauWMPI()
! MODULES
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_WMLES_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
! Every VAR we'll be using is declared in MOD_WMLES_Vars, loaded here.
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iProc
!==================================================================================================================================
DO iProc=1,nProcs_SendTauW
    CALL MPI_ISEND(TauW_MINE(:,:,iProc),2*nTauW_MINE(iProc), MPI_DOUBLE_PRECISION, &
                    Proc_SendTauW(iProc), 0, MPI_COMM_FLEXI, WMLES_SendRequests(iProc),iError)
END DO

END SUBROUTINE StartSendTauWMPI

!==================================================================================================================================
!> Subroutine that completes the receive operations for the Tau_W to be exchanged between processors.
!==================================================================================================================================
SUBROUTINE FinishReceiveTauWMPI()
! MODULES
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_WMLES_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
! Every VAR we'll be using is declared in MOD_WMLES_Vars, loaded here.
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                             :: i
!==================================================================================================================================
CALL MPI_Waitall(nProcs_RecvTauW,WMLES_RecvRequests,MPI_STATUSES_IGNORE,iError)

! In addition to guaranteeing the receive operations so that TauW_YOURS is populated,
! we check if there is any tau_w point where the MPI rank that calculates is the same
! as the rank responsible for BC imposition (i.e. nTauW_MINE(0) .NE. 0), and hence 
! no comm. is needed, since the info is local to this MPI rank.

IF (nTauW_MINE(0).NE.0) THEN
    DO i=1,nTauW_MINE(0)
        TauW_YOURS(:,i,0) = TauW_MINE(:,i,0) ! Populate TauW_YOURS(:,:,0)
    END DO
END IF

END SUBROUTINE FinishReceiveTauWMPI


#if USE_MPI
!==================================================================================================================================
!> Find the id of a processor on which an element with a given ElemID lies, based on the MPI element offsets defined earlier.
!> Use a bisection algorithm for faster search.
!==================================================================================================================================
FUNCTION ELEMIPROC(ElemID)
! MODULES
USE MOD_Globals,   ONLY:nProcessors
USE MOD_MPI_vars,  ONLY:offsetElemMPI
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER, INTENT(IN)                :: ElemID     !< (IN)  NodeID to search for
INTEGER                            :: ELEMIPROC  !< (OUT) processor id
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: i,maxSteps,low,up,mid
!==================================================================================================================================
ELEMIPROC=0
maxSteps=INT(LOG(REAL(nProcessors))*1.4426950408889634556)+1    !1/LOG(2.)=1.4426950408889634556
low=0
up=nProcessors-1
IF((ElemID.GT.offsetElemMPI(low)).AND.(ElemID.LE.offsetElemMPI(low+1)))THEN
  ELEMIPROC=low
ELSEIF((ElemID.GT.offsetElemMPI(up)).AND.(ElemID.LE.offsetElemMPI(up+1)))THEN
  ELEMIPROC=up
ELSE
  !bisection
  DO i=1,maxSteps
    mid=(up-low)/2+low
    IF((ElemID.GT.offsetElemMPI(mid)).AND.(ElemID.LE.offsetElemMPI(mid+1)))THEN
      ELEMIPROC=mid                     !index found!
      EXIT
    ELSEIF(ElemID .GT. offsetElemMPI(mid+1))THEN ! seek in upper half
      low=mid+1
    ELSE
      up=mid
    END IF
  END DO
END IF
END FUNCTION ELEMIPROC
#endif /*USE_MPI*/

!==================================================================================================================================
!> Finalize parameters of WMLES computation
!==================================================================================================================================
SUBROUTINE FinalizeWMLES()
! MODULES
USE MOD_WMLES_Vars
IMPLICIT NONE

!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

!==================================================================================================================================
SDEALLOCATE(WMLES_Tauw)


END SUBROUTINE FinalizeWMLES



!==================================================================================================================================
END MODULE MOD_WMLES