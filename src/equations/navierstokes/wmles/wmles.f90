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
CALL prms%CreateIntFromStringOption('WallModel', "Wall model to be used on walls defined with approximate boundary conditions", 'LogLaw')
CALL addStrListEntry('WallModel', 'Schumann', WMLES_SCHUMANN)
CALL addStrListEntry('WallModel', 'LogLaw', WMLES_LOGLAW)
CALL addStrListEntry('WallModel', 'WernerWangle', WMLES_WERNERWANGLE)
CALL addStrListEntry('WallModel', 'Reichardt', WMLES_REICHARDT)
CALL addStrListEntry('WallModel', 'Spalding', WMLES_SPALDING)
CALL addStrListEntry('WallModel', 'EquilibriumTBLE', WMLES_EQTBLE)
CALL addStrListEntry('WallModel', 'Couette', WMLES_COUETTE)

CALL prms%CreateRealOption('h_wm', "Distance from the wall at which LES and Wall Model "// &
                                    "exchange instantaneous flow information", "0.2")
CALL prms%CreateRealOption('delta', "Estimated Boundary Layer Thickness")
CALL prms%CreateRealOption('vKarman', "von Karman constant", "0.41")
CALL prms%CreateRealOption('B', "Log-law intercept coefficient", "5.2")

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
USE MOD_Mesh_Vars              ! ,ONLY: nBCSides, BC, BoundaryType, MeshInitIsDone, ElemToSide, SideToElem, SideToGlobalSide
USE MOD_Interpolation_Vars      ,ONLY: InterpolationInitIsDone,xGP,wBary
USE MOD_Basis                   ,ONLY: LagrangeInterpolationPolys
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
!IF((.NOT.InterpolationInitIsDone).OR.(.NOT.MeshInitIsDone).OR.WMLESInitDone) THEN
!  CALL CollectiveStop(__STAMP__,&
!    'Wall-Modeled LES not ready to be called or already called.')
!END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT Wall-Modeled LES...'

WallModel = GETINTFROMSTR('WallModel')

! If Schumann's model is selected, check if we are in a channel flow
IF ((WallModel .EQ. WMLES_SCHUMANN) .AND. (.NOT.STRICMP(Testcase, Channel))) THEN
    CALL CollectiveStop(__STAMP__,&
        "Schumann's wall model can only be applied to the Channel flow testcase.")
END IF

vKarman = GETREAL('vKarman')
B = GETREAL('B')

h_wm = GETREAL('h_wm')
IF (H_WM .LE. 0) &
    CALL CollectiveStop(__STAMP__,&
        'h_wm distance from the wall must be positive.')

delta = GETREAL('delta')

abs_h_wm = h_wm*delta
NSuper = 3*PP_N

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
            ! Inward normal vector (hence the minus sign)
            OthersPointInfo(4:6,WallStressCount_local(hwmRank),hwmRank) = -NormVec(:,p,q,0,iSide)
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
ALLOCATE(TauW_MINE_Interpolate(4,MAXVAL(nTauW_MINE),0:nProcs_SendTauW))
ALLOCATE(TauW_MINE_NormVec(3,MAXVAL(nTauW_MINE),0:nProcs_SendTauW))

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
        IF (nTauW_MINE(0).NE.0) PointInfo(:,1:WallStressCount_local(myRank)) = OthersPointInfo(:,1:WallStressCount_local(myRank),myRank)
    END IF

    ! Store wall-normal vector
    IF (nTauW_MINE(i).NE.0) TauW_MINE_NormVec(1:3,1:nTauW_MINE(i),i) = PointInfo(4:6,1:nTauW_MINE(i))
    
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
            IF (Logging) TauW_MINE_IsInterpolation(j,i) = .TRUE.
            nTauW_MINE_Interpolate(i) = nTauW_MINE_Interpolate(i) + 1
            InterpToLocalPoint(nTauW_MINE_Interpolate(i)) = j
            !> Map h_wm coords from physical to local (standard) coords and save this info
            CALL PhysToStdCoords(PointInfo(1:3,j),Loc_hwmElemID,TauW_MINE_Interpolate(1:3,nTauW_MINE_Interpolate(i),i))
            TauW_MINE_Interpolate(4,nTauW_MINE_Interpolate(i),i) = Loc_hwmElemID
        END IF
    END DO ! j -- each point that we must calculate tau_w, for this MPI process i

    ! Sanity check
    IF ((nTauW_MINE_FacePoint(i)+nTauW_MINE_Interpolate(i)+nTauW_MINE_InteriorPoint(i)) &
                .NE. nTauW_MINE(i)) CALL Abort(__STAMP__,"Dangling points?")

END DO ! i -- message from each process having a BC side which needs our calculation of tau_w

! Check if there are interpolation points to calculate in this MPI proc., so that we set up
! and store the Lagrange polynomials calculate at each point
IF (MAXVAL(nTauW_MINE_Interpolate).GT.0) THEN 
    ! Set up the Lagrangian interpolating polynomials for each interp. node.
    ALLOCATE(Lag_xi(0:PP_N,MAXVAL(nTauW_MINE_Interpolate),0:nProcs_SendTauW))
    ALLOCATE(Lag_eta(0:PP_N,MAXVAL(nTauW_MINE_Interpolate),0:nProcs_SendTauW))
    ALLOCATE(Lag_zeta(0:PP_N,MAXVAL(nTauW_MINE_Interpolate),0:nProcs_SendTauW))
    DO i=0,nProcs_SendTauW
        DO j=1,nTauW_MINE_Interpolate(i)
            CALL LagrangeInterpolationPolys(TauW_MINE_Interpolate(1,j,i),PP_N,xGP,wBary,Lag_xi(:,j,i))
            CALL LagrangeInterpolationPolys(TauW_MINE_Interpolate(2,j,i),PP_N,xGP,wBary,Lag_eta(:,j,i))
            CALL LagrangeInterpolationPolys(TauW_MINE_Interpolate(3,j,i),PP_N,xGP,wBary,Lag_zeta(:,j,i))
        END DO
    END DO
END IF

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
CALL MPI_Waitall(nProcessors,CalcInfoRequests(0:nProcessors-1),MPI_STATUSES_IGNORE,iError)
CALL MPI_Waitall(nProcessors,PointInfoRequests(0:nProcessors-1),MPI_STATUSES_IGNORE,iError)

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
INTEGER                             :: p,q,r,i,SideID,ElemID
REAL                                :: u_tau, u_mean, tau_w_mag, utang, VelMag
REAL                                :: vel_inst(3), rho_inst, tangvec(3), tau_w_vec(3)
!==================================================================================================================================
IF (.NOT.WMLESInitDone) THEN
    CALL CollectiveStop(__STAMP__,&
    'Cannot compute Wall Stress -- Wall-Modeled LES INIT is not done.')
END IF

! Start non-blockingly receiving WMLES info (for those who have anything to receive)
CALL StartReceiveTauWMPI()

! TODO ?
! For algebraic models (except for Schumann), the outcome is really u_tau instead of tau_w directly
! Hence, to compute tau_w, we take tau_w = u_tau^2 * rho. However, this rho should be 
! the density AT THE WALL, and we are currently assuming the density at the exchange location, h_wm,
! which is conceptually wrong for a compressible flow.
! Thus, maybe we should send u_tau instead of tau_w to the MPI proc. responsible for BC imposition,
! and then this MPI proc. computes tau_w.
! (See Bocquet, Sagaut & Jouhaud).

! Compute tau_w for points of our responsibility
SELECT CASE(WallModel)

  CASE(WMLES_SCHUMANN)
    
    ! Schumann model
    ! Computation is done using instantaneous information, projected onto the wall-tangent direction
    ! Mean wall shear stress is an INPUT to this model ! Therefore, only works for periodic channel flows.
    ! The channel flow testcase assumes <tau_w> = 1.0
    ! u_tau is used in the log-law eq. to extract <u> (u_mean), which is then used to calculate tau_w

    DO sProc=0,nProcs_SendTauW
        ! Calculate tau_w for each h_wm that is approximated as a face node
        DO FPInd=1,nTauW_MINE_FacePoint(sProc)
            p = TauW_MINE_FacePoint(1,FPInd,sProc)
            q = TauW_MINE_FacePoint(2,FPInd,sProc)
            SideID = TauW_MINE_FacePoint(3,FPInd,sProc)

            ! Face_xGP is only populated for master sides, so that only master sides have approximation nodes (check InitWMLES above)
            ! hence, use of UPrim_master is guaranteed to be correct here
            !utang = DOT_PRODUCT(UPrim_master(2:4,p,q,SideID),TauW_MINE_NormVec(:,FaceToLocalPoint(FPInd),sProc))
            utang = UPrim_master(2,p,q,SideID)

            u_tau = SQRT(1.0/UPrim_master(1,p,q,SideID)) 
            u_mean = u_tau*( (1./vKarman) * LOG((abs_h_wm*u_tau)/mu0) + B )
            TauW_MINE(1,FaceToLocalPoint(FPInd),sProc) = utang/u_mean * 1.0 ! <tau_w> = 1.0
            TauW_MINE(2,FaceToLocalPoint(FPInd),sProc) = (2.0*mu0)*(UPrim_master(4,p,q,SideID)/abs_h_wm)
        END DO

        ! Calculate tau_w for each h_wm that is approximated as an interior node
        DO IPInd=1,nTauW_MINE_InteriorPoint(sProc)
            p = TauW_MINE_InteriorPoint(1,IPInd,sProc)
            q = TauW_MINE_InteriorPoint(2,IPInd,sProc)
            r = TauW_MINE_InteriorPoint(3,IPInd,sProc)
            ElemID = TauW_MINE_InteriorPoint(4,IPInd,sProc)

            !utang = DOT_PRODUCT(UPrim(2:4,p,q,r,ElemID),TauW_MINE_NormVec(:,InteriorToLocalPoint(IPInd),sProc))
            utang = UPrim(2,p,q,r,ElemID)

            u_tau = SQRT(1.0/UPrim(1,p,q,r,ElemID)) 
            u_mean = u_tau*( (1./vKarman) * LOG((abs_h_wm*u_tau)/mu0) + B )
            TauW_MINE(1,InteriorToLocalPoint(IPInd),sProc) = utang/u_mean * 1.0 ! <tau_w> = 1.0
            TauW_MINE(2,InteriorToLocalPoint(IPInd),sProc) = (2.0*mu0)*(UPrim(4,p,q,r,ElemID)/abs_h_wm)
        END DO

        ! Calculate tau_w for each h_wm that must be interpolated
        DO IntPInd=1,nTauW_MINE_Interpolate(sProc)
            ElemID = INT(TauW_MINE_Interpolate(4,IntPInd,sProc))
            vel_inst(1) = InterpolateHwm(ElemID,Lag_xi(:,IntPInd,sProc),Lag_eta(:,IntPInd,sProc),Lag_zeta(:,IntPInd,sProc),2,prim=.TRUE.)
            vel_inst(2) = InterpolateHwm(ElemID,Lag_xi(:,IntPInd,sProc),Lag_eta(:,IntPInd,sProc),Lag_zeta(:,IntPInd,sProc),3,prim=.TRUE.)
            vel_inst(3) = InterpolateHwm(ElemID,Lag_xi(:,IntPInd,sProc),Lag_eta(:,IntPInd,sProc),Lag_zeta(:,IntPInd,sProc),4,prim=.TRUE.)
            rho_inst = InterpolateHwm(ElemID,Lag_xi(:,IntPInd,sProc),Lag_eta(:,IntPInd,sProc),Lag_zeta(:,IntPInd,sProc),1)

            !utang = DOT_PRODUCT(vel_inst(:),TauW_MINE_NormVec(:,InterpToLocalPoint(IntPInd),sProc))
            utang = vel_inst(1)

            u_tau = SQRT(1.0/rho_inst) 
            u_mean = u_tau*( (1./vKarman) * LOG((abs_h_wm*u_tau)/mu0) + B )
            TauW_MINE(1,InterpToLocalPoint(IntPInd),sProc) = utang/u_mean * 1.0 ! <tau_w> = 1.0
            TauW_MINE(2,InterpToLocalPoint(IntPInd),sProc) = (2.0*mu0)*(vel_inst(3)/abs_h_wm)
        END DO
    END DO

  CASE (WMLES_LOGLAW)

    ! In this case, it is assumed that the log-law is an instantaneous relation
    ! Then, u (projected onto the wall-tangent streamwise direction) is used in the log-law
    ! and u_tau is computed (using Newton iterations).
    ! Then, u_tau is used to compute wall shear stresses.

    DO sProc=0,nProcs_SendTauW
        ! Calculate tau_w for each h_wm that is approximated as a face node
        DO FPInd=1,nTauW_MINE_FacePoint(sProc)
            p = TauW_MINE_FacePoint(1,FPInd,sProc)
            q = TauW_MINE_FacePoint(2,FPInd,sProc)
            SideID = TauW_MINE_FacePoint(3,FPInd,sProc)

            ! Face_xGP is only populated for master sides, so that only master sides have approximation nodes (check InitWMLES above)
            ! hence, use of UPrim_master is guaranteed to be correct here
            tangvec = UPrim_master(2:4,p,q,SideID) - DOT_PRODUCT(UPrim_master(2:4,p,q,SideID),TauW_MINE_NormVec(:,FaceToLocalPoint(FPInd),sProc))*TauW_MINE_NormVec(:,FaceToLocalPoint(FPInd),sProc)
            VelMag = 0.
            DO i=1,3
                VelMag = VelMag + tangvec(i)**2
            END DO
            VelMag = SQRT(VelMag)
            tangvec = tangvec/VelMag
            !<-=-=-=-=- DEBUG; REMOVE LATER IF IT DOES NOT WORK =-=-=-=-=>!
            ! Force tangvec to be in the x dir
            tangvec = (/1.,0.,0./)
            !<-=-=-=-=- END OF DEBUG; REMOVE LATER IF IT DOES NOT WORK =-=-=-=-=>!

            utang = DOT_PRODUCT(UPrim_master(2:4,p,q,SideID),tangvec)

            u_tau = NewtonLogLaw(utang,(mu0/UPrim_master(1,p,q,SideID)))
            tau_w_mag = UPrim_master(1,p,q,SideID)*(u_tau**2) ! CHECK TODO ABOVE
            tau_w_vec = tau_w_mag*tangvec

            !<-=-=-=-=- DEBUG; REMOVE LATER IF IT DOES NOT WORK =-=-=-=-=>!
            ! Create vector aligned with wall-normal velocity with tau_xy
            tau_w_vec = (/0.,tau_w_mag,0./)
            ! Project it onto the normal direction so that the correct signal is imposed
            TauW_MINE(1,FaceToLocalPoint(FPInd),sProc) = -1.*DOT_PRODUCT(tau_w_vec(1:3),-TauW_MINE_NormVec(1:3,FaceToLocalPoint(FPInd),sProc))
            !<-=-=-=-=- END OF DEBUG; REMOVE LATER IF IT DOES NOT WORK =-=-=-=-=>!

            ! TODO
            ! Associated with the TODO above: should we pass u_tau or tau_w_vec to the MPI proc.
            ! responsible for BC imposition, so that it can calculate TauW with the knowledge of 
            ! the face-tangent vectors?
            ! Right now, as a work-around, we assume that TauW_MINE(1) (i.e. tau_xy) = tau_w_mag
            ! and TauW_MINE(2) = 0, which will be OK for the channel testcase, but definitely NOT OK
            ! for any other flow that is not completely aligned with the x direction.
            !TauW_MINE(1,FaceToLocalPoint(FPInd),sProc) = tau_w_mag
            TauW_MINE(2,FaceToLocalPoint(FPInd),sProc) = 0.
        END DO

        ! Calculate tau_w for each h_wm that is approximated as an interior node
        DO IPInd=1,nTauW_MINE_InteriorPoint(sProc)
            p = TauW_MINE_InteriorPoint(1,IPInd,sProc)
            q = TauW_MINE_InteriorPoint(2,IPInd,sProc)
            r = TauW_MINE_InteriorPoint(3,IPInd,sProc)
            ElemID = TauW_MINE_InteriorPoint(4,IPInd,sProc)

            tangvec = UPrim(2:4,p,q,r,ElemID) - DOT_PRODUCT(UPrim(2:4,p,q,r,ElemID),TauW_MINE_NormVec(:,InteriorToLocalPoint(IPInd),sProc))*TauW_MINE_NormVec(:,InteriorToLocalPoint(IPInd),sProc)
            VelMag = 0.
            DO i=1,3
                VelMag = VelMag + tangvec(i)**2
            END DO
            VelMag = SQRT(VelMag)
            tangvec = tangvec/VelMag
            !<-=-=-=-=- DEBUG; REMOVE LATER IF IT DOES NOT WORK =-=-=-=-=>!
            ! Force tangvec to be in the x dir
            tangvec = (/1.,0.,0./)
            !<-=-=-=-=- END OF DEBUG; REMOVE LATER IF IT DOES NOT WORK =-=-=-=-=>!

            utang = DOT_PRODUCT(UPrim(2:4,p,q,r,ElemID),tangvec)

            u_tau = NewtonLogLaw(utang,(mu0/UPrim(1,p,q,r,ElemID)))
            tau_w_mag = UPrim(1,p,q,r,ElemID)*(u_tau**2) ! CHECK TODO ABOVE
            tau_w_vec = tau_w_mag*tangvec

            !<-=-=-=-=- DEBUG; REMOVE LATER IF IT DOES NOT WORK =-=-=-=-=>!
            ! Create vector aligned with wall-normal velocity with tau_xy
            tau_w_vec = (/0.,tau_w_mag,0./)
            ! Project it onto the normal direction so that the correct signal is imposed
            TauW_MINE(1,InteriorToLocalPoint(IPInd),sProc) = -1.*DOT_PRODUCT(tau_w_vec(1:3),-TauW_MINE_NormVec(1:3,InteriorToLocalPoint(IPInd),sProc))
            !<-=-=-=-=- END OF DEBUG; REMOVE LATER IF IT DOES NOT WORK =-=-=-=-=>!

            ! TODO
            ! Associated with the TODO above: should we pass u_tau or tau_w_vec to the MPI proc.
            ! responsible for BC imposition, so that it can calculate TauW with the knowledge of 
            ! the face-tangent vectors?
            ! Right now, as a work-around, we assume that TauW_MINE(1) (i.e. tau_xy) = tau_w_mag
            ! and TauW_MINE(2) = 0, which will be OK for the channel testcase, but definitely NOT OK
            ! for any other flow that is not completely aligned with the x direction.
            !TauW_MINE(1,InteriorToLocalPoint(IPInd),sProc) = tau_w_mag
            TauW_MINE(2,InteriorToLocalPoint(IPInd),sProc) = 0.
        END DO

        ! Calculate tau_w for each h_wm that must be interpolated
        DO IntPInd=1,nTauW_MINE_Interpolate(sProc)
            ElemID = INT(TauW_MINE_Interpolate(4,IntPInd,sProc))
            vel_inst(1) = InterpolateHwm(ElemID,Lag_xi(:,IntPInd,sProc),Lag_eta(:,IntPInd,sProc),Lag_zeta(:,IntPInd,sProc),2,prim=.TRUE.)
            vel_inst(2) = InterpolateHwm(ElemID,Lag_xi(:,IntPInd,sProc),Lag_eta(:,IntPInd,sProc),Lag_zeta(:,IntPInd,sProc),3,prim=.TRUE.)
            vel_inst(3) = InterpolateHwm(ElemID,Lag_xi(:,IntPInd,sProc),Lag_eta(:,IntPInd,sProc),Lag_zeta(:,IntPInd,sProc),4,prim=.TRUE.)
            rho_inst = InterpolateHwm(ElemID,Lag_xi(:,IntPInd,sProc),Lag_eta(:,IntPInd,sProc),Lag_zeta(:,IntPInd,sProc),1)

            tangvec = vel_inst(:) - DOT_PRODUCT(vel_inst(:),TauW_MINE_NormVec(:,InterpToLocalPoint(IntPInd),sProc))*TauW_MINE_NormVec(:,InterpToLocalPoint(IntPInd),sProc)
            VelMag = 0.
            DO i=1,3
                VelMag = VelMag + tangvec(i)**2
            END DO
            VelMag = SQRT(VelMag)
            tangvec = tangvec/VelMag

            !<-=-=-=-=- DEBUG; REMOVE LATER IF IT DOES NOT WORK =-=-=-=-=>!
            ! Force tangvec to be in the x dir
            tangvec = (/1.,0.,0./)
            !<-=-=-=-=- END OF DEBUG; REMOVE LATER IF IT DOES NOT WORK =-=-=-=-=>!

            utang = DOT_PRODUCT(vel_inst,tangvec)

            u_tau = NewtonLogLaw(utang,(mu0/rho_inst))
            tau_w_mag = rho_inst*(u_tau**2) ! CHECK TODO ABOVE
            tau_w_vec = tau_w_mag*tangvec

            !<-=-=-=-=- DEBUG; REMOVE LATER IF IT DOES NOT WORK =-=-=-=-=>!
            ! Create vector aligned with wall-normal velocity with tau_xy
            tau_w_vec = (/0.,tau_w_mag,0./)
            ! Project it onto the normal direction so that the correct signal is imposed
            TauW_MINE(1,InterpToLocalPoint(IntPInd),sProc) = -1.*DOT_PRODUCT(tau_w_vec(1:3),-TauW_MINE_NormVec(1:3,InterpToLocalPoint(IntPInd),sProc))
            !<-=-=-=-=- END OF DEBUG; REMOVE LATER IF IT DOES NOT WORK =-=-=-=-=>!

            ! TODO
            ! Associated with the TODO above: should we pass u_tau or tau_w_vec to the MPI proc.
            ! responsible for BC imposition, so that it can calculate TauW with the knowledge of 
            ! the face-tangent vectors?
            ! Right now, as a work-around, we assume that TauW_MINE(1) (i.e. tau_xy) = tau_w_mag
            ! and TauW_MINE(2) = 0, which will be OK for the channel testcase, but definitely NOT OK
            ! for any other flow that is not completely aligned with the x direction.
            !TauW_MINE(1,InterpToLocalPoint(IntPInd),sProc) = tau_w_mag
            TauW_MINE(2,InterpToLocalPoint(IntPInd),sProc) = 0.
        END DO
    END DO

  CASE(WMLES_COUETTE)

    ! Test for a laminar flow with known velocity profile, so that we can debug the implementation

    DO sProc=0,nProcs_SendTauW
        ! Calculate tau_w for each h_wm that is approximated as a face node
        DO FPInd=1,nTauW_MINE_FacePoint(sProc)
            p = TauW_MINE_FacePoint(1,FPInd,sProc)
            q = TauW_MINE_FacePoint(2,FPInd,sProc)
            SideID = TauW_MINE_FacePoint(3,FPInd,sProc)

            ! Hard coded value (analytic solution)
            ! TauW_MINE(1,FaceToLocalPoint(FPInd),sProc) = 0.5
            ! IF (Face_xGP(2,p,q,0,SideID) .GE. 0) TauW_MINE(1,FaceToLocalPoint(FPInd),sProc) = -0.5

            ! Calculated as if a model was employed
            utang = UPrim_master(2,p,q,SideID)

            u_tau = NewtonCouette(utang,mu0) ! dpdx, obviously.

            ! Create vector aligned with wall-normal direction and magnitude of tau_xy
            tau_w_vec = (/0.,-0.5*u_tau,0./)
            ! Project it onto the normal direction so that the correct signal is imposed
            TauW_MINE(1,FaceToLocalPoint(FPInd),sProc) = -1.*DOT_PRODUCT(tau_w_vec(1:3),-TauW_MINE_NormVec(1:3,FaceToLocalPoint(FPInd),sProc))

            ! TauW_MINE(1,FaceToLocalPoint(FPInd),sProc) = -0.5*u_tau
            ! IF(Face_xGP(2,p,q,0,SideID) .GE. 0) TauW_MINE(1,FaceToLocalPoint(FPInd),sProc) = 0.5*u_tau
            TauW_MINE(2,FaceToLocalPoint(FPInd),sProc) = 0.

        END DO

        ! Calculate tau_w for each h_wm that is approximated as an interior node
        DO IPInd=1,nTauW_MINE_InteriorPoint(sProc)
            p = TauW_MINE_InteriorPoint(1,IPInd,sProc)
            q = TauW_MINE_InteriorPoint(2,IPInd,sProc)
            r = TauW_MINE_InteriorPoint(3,IPInd,sProc)
            ElemID = TauW_MINE_InteriorPoint(4,IPInd,sProc)

            ! Hard coded value (analytic solution)
            ! TauW_MINE(1,InteriorToLocalPoint(IPInd),sProc) = 0.5
            ! IF (Elem_xGP(2,p,q,r,ElemID) .GE. 0) TauW_MINE(1,InteriorToLocalPoint(IPInd),sProc) = -0.5

            ! Calculated as if a model was employed
            utang = UPrim(2,p,q,r,ElemID)

            u_tau = NewtonCouette(utang,mu0) ! dpdx, obviously.

            ! Create vector aligned with wall-normal direction and magnitude of tau_xy
            tau_w_vec = (/0.,-0.5*u_tau,0./)
            ! Project it onto the normal direction so that the correct signal is imposed
            TauW_MINE(1,InteriorToLocalPoint(IPInd),sProc) = -1.*DOT_PRODUCT(tau_w_vec(1:3),-TauW_MINE_NormVec(1:3,InteriorToLocalPoint(IPInd),sProc))

            ! TauW_MINE(1,InteriorToLocalPoint(IPInd),sProc) = -0.5*u_tau
            ! IF(Elem_xGP(2,p,q,r,ElemID) .GE. 0) TauW_MINE(1,InteriorToLocalPoint(IPInd),sProc) = 0.5*u_tau
            TauW_MINE(2,InteriorToLocalPoint(IPInd),sProc) = 0.

        END DO

        ! Calculate tau_w for each h_wm that must be interpolated
        DO IntPInd=1,nTauW_MINE_Interpolate(sProc)
            ElemID = INT(TauW_MINE_Interpolate(4,IntPInd,sProc))
            vel_inst(1) = InterpolateHwm(ElemID,Lag_xi(:,IntPInd,sProc),Lag_eta(:,IntPInd,sProc),Lag_zeta(:,IntPInd,sProc),2,prim=.TRUE.)

            ! Hard coded value (analytic solution)
            ! TauW_MINE(1,InterpToLocalPoint(IntPInd),sProc) = 0.5
            ! IF (Elem_xGP(2,1,1,1,ElemID) .GE. 0) TauW_MINE(1,InterpToLocalPoint(IntPInd),sProc) = -0.5

            u_tau = NewtonCouette(vel_inst(1),mu0) ! dpdx, obviously.

            ! Create vector aligned with wall-normal direction and magnitude of tau_xy
            tau_w_vec = (/0.,-0.5*u_tau,0./)
            ! Project it onto the normal direction so that the correct signal is imposed
            TauW_MINE(1,InterpToLocalPoint(IntPInd),sProc) = -1.*DOT_PRODUCT(tau_w_vec(1:3),-TauW_MINE_NormVec(1:3,InterpToLocalPoint(IntPInd),sProc))

            ! TauW_MINE(1,InterpToLocalPoint(IntPInd),sProc) = -0.5*u_tau
            ! IF(Elem_xGP(2,p,q,r,ElemID) .GE. 0) TauW_MINE(1,InterpToLocalPoint(IntPInd),sProc) = 0.5*u_tau
            TauW_MINE(2,InterpToLocalPoint(IntPInd),sProc) = 0.
        END DO
    END DO


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
  f(5,p,q) = 0. ! -(UPrim_master(3,p,q)*tau_xy) ! -(u*tau_xx + v*tau_xy + w*tau_xz - kdTdx)

  g(1,p,q) = 0.
  g(2,p,q) = -tau_xy ! tau_yx = tau_xy
  g(3,p,q) = 0. ! tau_yy = 0
  g(4,p,q) = -tau_yz ! tau_yz = tau_zy
  g(5,p,q) = 0. ! -(UPrim_master(2,p,q)*tau_xy + UPrim_master(4,p,q)*tau_yz) ! -(u*tau_yx + v*tau_yy + w*tau_yz - kdTdx)

  h(1,p,q) = 0.
  h(2,p,q) = 0. ! tau_zx = 0
  h(3,p,q) = -tau_yz ! tau_zy = tau_yz
  h(4,p,q) = 0. ! tau_zz = 0
  h(5,p,q) = 0. !-(UPrim_master(3,p,q)*tau_yz) ! -(u*tau_zx + v*tau_zy + w*tau_zz - kdTdx)


  LOGWRITE(*,'(2(I4,2X),15(E15.8,2X))') p, q, f(1,p,q), f(2,p,q), f(3,p,q), f(4,p,q), f(5,p,q),&
     g(1,p,q), g(2,p,q), g(3,p,q), g(4,p,q), g(5,p,q), h(1,p,q), h(2,p,q), h(3,p,q), h(4,p,q), h(5,p,q)

END DO; END DO ! p,q

LOGWRITE(*,'(X)')

END SUBROUTINE EvalDiffFlux3D_WMLES

!===================================================================================================================================
!> Computes parametric coordinates of exchange location point (h_wm), given the physical coordinates, iteratively.
!> This is a copy/paste version of GetParametricCoordinates from posti/recordpoints, modified to use only the parts we need.
!===================================================================================================================================
SUBROUTINE PhysToStdCoords(hwmPhys,iElem,hwmStd)
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Mathtools,         ONLY: INVERSE
USE MOD_Mesh_Vars,         ONLY: NGeo,Elem_xGP,SideToElem,nBCSides,Face_xGP,NormVec
USE MOD_Basis,             ONLY: LagrangeInterpolationPolys,ChebyGaussLobNodesAndWeights,BarycentricWeights
USE MOD_Basis,             ONLY: PolynomialDerivativeMatrix
USE MOD_Interpolation,     ONLY: GetVandermonde,GetNodesAndWeights,GetDerivativeMatrix
USE MOD_Interpolation_Vars
USE MOD_ChangeBasis,       ONLY: ChangeBasis3D,ChangeBasis2D
USE MOD_WMLES_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL, INTENT(IN)            :: hwmPhys(3)
INTEGER, INTENT(IN)         :: iElem
REAL, INTENT(OUT)           :: hwmStd(3)
!----------------------------------------------------------------------------------------------------------------------------------
REAL                  :: DCL_NGeo(0:Ngeo,0:Ngeo)
REAL                  :: XCL_Ngeo(3,0:Ngeo,0:Ngeo,0:Ngeo)          !< mapping X(xi) P\in Ngeo
REAL                  :: hmax2
REAL                  :: X_NSuper(1:3,0:NSuper,0:NSuper,0:NSuper)
REAL                  :: Winner_Dist2,Dist2
REAL                  :: dXCL_NGeo(1:3,1:3,0:NGeo,0:NGeo,0:NGeo)
REAL                  :: Xi_CLNGeo(0:NGeo),wBary_CLNGeo(0:NGeo)
REAL                  :: Xi_NSuper(0:NSuper)
REAL                  :: xC(3)
REAL                  :: Vdm_N_CLNGeo(0:NGeo,0:PP_N)
REAL                  :: Vdm_CLNGeo_EquiNSuper(0:NSuper,0:NGeo)
REAL                  :: Lag(1:3,0:NGeo)
REAL                  :: F(1:3),eps_F,Xi(1:3),Jac(1:3,1:3),sJac(1:3,1:3)
INTEGER               :: i,j,k,l,nNodes
INTEGER               :: NewtonIter
CHARACTER(LEN=255)    :: NodeType_Super
LOGICAL               :: changeBasisDone,calcJacobianDone
LOGICAL               :: hwmFound
REAL, PARAMETER       :: maxTol = 1.+1.E-3 ! This souldn't be too restricting
! boundary projection
INTEGER               :: SideID,locSideID,iWinner,jWinner
REAL                  :: dist2RP
REAL                  :: Vdm_N_EquiNSuper(0:NSuper,0:PP_N),xBC_NSuper(3,0:NSuper,0:NSuper)
REAL                  :: NormVec_NSuper(3,0:NSuper,0:NSuper)
!newton method
REAL                  :: wBary_NSuper(0:NSuper), D_NSuper(0:NSuper,0:NSuper), Lag_NSuper(1:2,0:NSuper)
REAL                  :: dxBC_NSuper(3,2,0:NSuper,0:NSuper)
REAL                  :: Gmat(2,0:NSuper,0:NSuper),dGmat(2,2,0:NSuper,0:NSuper)
REAL                  :: G(2),Xi2(2),Jac2(2,2),sJac2(2,2),xWinner(3),NormVecWinner(3)
!===================================================================================================================================
! Prepare CL basis evaluation
CALL GetNodesAndWeights( NGeo, NodeTypeCL, Xi_CLNGeo,wIPBary=wBary_CLNGeo)
CALL GetDerivativeMatrix(Ngeo, NodeTypeCL, DCL_Ngeo)

! Set NSuper to three times the current polynomial degree
DO i=0,NSuper
  Xi_NSuper(i) = 2./REAL(NSuper) * REAL(i) - 1.
END DO

NodeType_Super='VISU'
CALL GetVanderMonde(PP_N,NodeType  ,NGeo  ,NodeTypeCL    ,Vdm_N_CLNGeo)
CALL GetVanderMonde(NGeo,NodeTypeCL,NSuper,NodeType_Super,Vdm_CLNGeo_EquiNSuper)
! Define maximum mesh size (square) hmax2
nNodes=(NGeo+1)**3

hwmFound=.FALSE.
hwmStd = HUGE(9.) ! First std coords "estimate"

! We have the Elem_xGP on N, compute XCL_NGeo for minimal changes
CALL ChangeBasis3D(3,PP_N,NGeo,Vdm_N_CLNGeo,Elem_xGP(:,:,:,:,iElem),XCL_NGeo)

changeBasisDone=.FALSE.
calcJacobianDone=.FALSE.


! Find initial guess for Newton, get supersampled element for first recordpoint
IF(.NOT.changeBasisDone) &
CALL ChangeBasis3D(3,NGeo,NSuper,Vdm_CLNGeo_EquiNSuper,XCL_NGeo,X_NSuper)
changeBasisDone=.TRUE.

Winner_Dist2=HUGE(1.)
DO i=0,NSuper; DO j=0,NSuper; DO k=0,NSuper
  Dist2=SUM((hwmPhys-X_NSuper(:,i,j,k))*(hwmPhys-X_NSuper(:,i,j,k)))
  IF (Dist2.LT.Winner_Dist2) THEN
    Winner_Dist2=Dist2
    Xi=(/Xi_NSuper(i),Xi_NSuper(j),Xi_NSuper(k)/)
END IF
END DO; END DO; END DO

! Compute Jacobian of Winner element Mapping for each CL point
IF(.NOT.calcJacobianDone)THEN
  dXCL_NGeo=0.
  DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
    ! Matrix-vector multiplication
    DO l=0,NGeo
      dXCL_NGeo(:,1,i,j,k)=dXCL_NGeo(:,1,i,j,k) + DCL_NGeo(i,l)*XCL_NGeo(:,l,j,k)
      dXCL_NGeo(:,2,i,j,k)=dXCL_NGeo(:,2,i,j,k) + DCL_NGeo(j,l)*XCL_NGeo(:,i,l,k)
      dXCL_NGeo(:,3,i,j,k)=dXCL_NGeo(:,3,i,j,k) + DCL_NGeo(k,l)*XCL_NGeo(:,i,j,l)
  END DO !l=0,NGeo
END DO; END DO; END DO
calcJacobianDone=.TRUE.
END IF

CALL LagrangeInterpolationPolys(Xi(1),NGeo,Xi_CLNGeo,wBary_CLNGeo,Lag(1,:))
CALL LagrangeInterpolationPolys(Xi(2),NGeo,Xi_CLNGeo,wBary_CLNGeo,Lag(2,:))
CALL LagrangeInterpolationPolys(Xi(3),NGeo,Xi_CLNGeo,wBary_CLNGeo,Lag(3,:))
! F(xi) = x(xi) - hwmPhys
F=-hwmPhys
DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
  F=F+XCL_NGeo(:,i,j,k)*Lag(1,i)*Lag(2,j)*Lag(3,k)
END DO; END DO; END DO

eps_F=1.E-10
NewtonIter=0
DO WHILE ((SUM(F*F).GT.eps_F).AND.(NewtonIter.LT.50))
  NewtonIter=NewtonIter+1
  ! Compute F Jacobian dx/dXi
  Jac=0.
  DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
    Jac=Jac+dXCL_NGeo(:,:,i,j,k)*Lag(1,i)*Lag(2,j)*Lag(3,k)
  END DO; END DO; END DO

  ! Compute inverse of Jacobian
  sJac=INVERSE(Jac)

  ! Iterate Xi using Newton step
  Xi = Xi - MATMUL(sJac,F)
  !  if Newton gets outside reference space range [-1,1], exit.
  ! But allow for some oscillation in the first couple of iterations, as we may discard the correct point/element!!
  IF((NewtonIter.GT.4).AND.(ANY(ABS(Xi).GT.1.2))) EXIT

  ! Compute function value
  CALL LagrangeInterpolationPolys(Xi(1),NGeo,Xi_CLNGeo,wBary_CLNGeo,Lag(1,:))
  CALL LagrangeInterpolationPolys(Xi(2),NGeo,Xi_CLNGeo,wBary_CLNGeo,Lag(2,:))
  CALL LagrangeInterpolationPolys(Xi(3),NGeo,Xi_CLNGeo,wBary_CLNGeo,Lag(3,:))
  ! F(xi) = x(xi) - hwmPhys
  F=-hwmPhys
  DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
      F=F+XCL_NGeo(:,i,j,k)*Lag(1,i)*Lag(2,j)*Lag(3,k)
  END DO; END DO; END DO
END DO !newton

! check if result is better then previous result
IF(MAXVAL(ABS(Xi)).LT.MAXVAL(ABS(hwmStd))) THEN
  IF(MAXVAL(ABS(Xi)).LE.1.) hwmFound = .TRUE. ! if point is inside element, stop searching
  hwmStd=Xi
  F=0.
  DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
    F=F+XCL_NGeo(:,i,j,k)*Lag(1,i)*Lag(2,j)*Lag(3,k)
  END DO; END DO; END DO
END IF



IF(.NOT.hwmFound)THEN
    ! mark as valid if greater than max tolerance
    IF(MAXVAL(ABS(hwmStd)).LT.maxTol) hwmFound=.TRUE.
END IF

! Remaining points: If the point is close to a boundary, project the point on the boundary
IF(.NOT.hwmFound) THEN
  WRITE(*,*) myRank, '!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!'
  LOGWRITE(UNIT_logOut,'(A,I4,A)') 'h_wm for',iElem,' have not been found inside the mesh. THIS SHOULD NOT HAPPEN.'
  LOGWRITE(UNIT_logOut,'(A)')' Attempting to project them on the closest boundary...'


  dist2RP=HUGE(1.)
  IF(NSuper.LT.Ngeo*2) THEN
    LOGWRITE(UNIT_logOut,*)'Warning: NSuper<2Ngeo, derivative may be wrong.'
END IF
! Prepare equidistant basis
CALL BarycentricWeights(NSuper,Xi_NSuper,wBary_NSuper)
CALL GetVandermonde(PP_N,NodeType,NSuper,NodeType_Super,Vdm_N_EquiNSuper)
CALL PolynomialDerivativeMatrix(NSuper,Xi_NSuper,D_NSuper)
nNodes=(NSuper+1)**2
DO SideID=1,nBCSides
    calcJacobianDone=.FALSE.
    ! Supersampling of the side to equidistant grid for search
    CALL ChangeBasis2D(3,PP_N,NSuper,Vdm_N_EquiNSuper,Face_xGP(:,:,:,0,SideID),xBC_NSuper)
    CALL ChangeBasis2D(3,PP_N,NSuper,Vdm_N_EquiNSuper,NormVec( :,:,:,0,SideID),NormVec_NSuper)
    ! get the BCSide centroid
    xC(1) = SUM(xBC_NSuper(1,:,:))/nNodes
    xC(2) = SUM(xBC_NSuper(2,:,:))/nNodes
    xC(3) = SUM(xBC_NSuper(3,:,:))/nNodes
    ! calculate the max. distance within the side to the centroid
    hmax2=0.
    DO j=0,NSuper
      DO i=0,NSuper
        hmax2=MAX(hmax2,SUM((xBC_NSuper(:,i,j)-xC)*(xBC_NSuper(:,i,j)-xC)))
    END DO
END DO
hmax2=hmax2*1.20 ! 20% tolerance

! Check if the point is close to the boundary side
! Coarse search of possible elems
IF(SUM((xC-hwmPhys)*((xC-hwmPhys))).GT.hmax2) CYCLE
! Get closest point on the supersampled side as starting value
locSideID= SideToElem(S2E_LOC_SIDE_ID,SideID)
Winner_Dist2=HUGE(1.)
DO i=0,NSuper; DO j=0,NSuper
    Dist2=SUM((hwmPhys-xBC_NSuper(:,i,j))*(hwmPhys-xBC_NSuper(:,i,j)))
    IF (Dist2.LT.Winner_Dist2) THEN
      Winner_Dist2=Dist2
      iWinner=i
      jWinner=j
  END IF
END DO; END DO
Xi2=(/Xi_NSuper(iWinner),Xi_NSuper(jWinner)/)
xWinner=xBC_NSuper(:,iWinner,jWinner)
NormVecWinner=NormVec_NSuper(:,iWinner,jWinner)
F=hwmPhys-xWinner

! Newton to find the minimum distance
! Calculate the surface jacobian
IF(.NOT.calcJacobianDone)THEN
    dxBC_NSuper=0.
    DO j=0,NSuper
      DO i=0,NSuper
          ! Matrix-vector multiplication
          DO l=0,NSuper
              dxBC_NSuper(:,1,i,j)=dxBC_NSuper(:,1,i,j) + D_NSuper(i,l)*xBC_NSuper(:,l,j)
              dxBC_NSuper(:,2,i,j)=dxBC_NSuper(:,2,i,j) + D_NSuper(j,l)*xBC_NSuper(:,i,l)
          END DO !l=0,NSuper
      END DO !i=0,NSuper
  END DO !j=0,NSuper
  calcJacobianDone=.TRUE.
END IF

! for Newton we first need the function Gmat(:,i,j) and its gradient in parameter space
! G= d/dXi((xBC-hwmPhys)²)=0, degree of G 2NGeo
Gmat=0.
DO j=0,NSuper
    DO i=0,NSuper
      Gmat(:,i,j)=Gmat(:,i,j)+dxBC_nSuper(1,:,i,j)*2*(xBC_NSuper(1,i,j)-hwmPhys(1)) &
      +dxBC_nSuper(2,:,i,j)*2*(xBC_NSuper(2,i,j)-hwmPhys(2)) &
      +dxBC_nSuper(3,:,i,j)*2*(xBC_NSuper(3,i,j)-hwmPhys(3))
  END DO! i=0,NSuper
END DO! j=0,NSuper

dGmat=0.
DO j=0,NSuper
    DO i=0,NSuper
      ! Matrix-vector multiplication
      DO l=0,NSuper
        dGmat(:,1,i,j)=dGmat(:,1,i,j) + D_NSuper(i,l)*Gmat(:,l,j)
        dGmat(:,2,i,j)=dGmat(:,2,i,j) + D_NSuper(j,l)*Gmat(:,i,l)
    END DO !l=0,NSuper
END DO! i=0,NSuper
END DO! j=0,NSuper
! get initial value of the functional G
CALL LagrangeInterpolationPolys(Xi2(1),NSuper,Xi_NSuper,wBary_NSuper,Lag_NSuper(1,:))
CALL LagrangeInterpolationPolys(Xi2(2),NSuper,Xi_NSuper,wBary_NSuper,Lag_NSuper(2,:))
G=0.
DO j=0,NSuper
    DO i=0,NSuper
      G=G+Gmat(:,i,j)*Lag_NSuper(1,i)*Lag_NSuper(2,j)
  END DO! i=0,NSuper
END DO! j=0,NSuper
eps_F=1.E-10
NewtonIter=0
DO WHILE ((SUM(G*G).GT.eps_F).AND.(NewtonIter.LT.50))
    NewtonIter=NewtonIter+1
    ! Compute G Jacobian dG/dXi

    Jac2=0.
    DO j=0,NSuper
      DO i=0,NSuper
        Jac2=Jac2 + dGmat(:,:,i,j)*Lag_NSuper(1,i)*Lag_NSuper(2,j)
    END DO !l=0,NSuper
END DO !i=0,NSuper

! Compute inverse of Jacobian
sJac2=INVERSE(Jac2)

! Iterate Xi using Newton step
Xi2 = Xi2 - MATMUL(sJac2,G)
! if Newton gets outside reference space range [-1,1], exit.
! But allow for some oscillation in the first couple of iterations, as we may discard the correct point/element!!
IF((NewtonIter.GT.4).AND.(ANY(ABS(Xi2).GT.1.2))) EXIT

! Compute function value
CALL LagrangeInterpolationPolys(Xi2(1),NSuper,Xi_NSuper,wBary_NSuper,Lag_NSuper(1,:))
CALL LagrangeInterpolationPolys(Xi2(2),NSuper,Xi_NSuper,wBary_NSuper,Lag_NSuper(2,:))
G=0.
DO j=0,NSuper
   DO i=0,NSuper
     G=G+Gmat(:,i,j)*Lag_NSuper(1,i)*Lag_NSuper(2,j)
 END DO! i=0,NSuper
END DO! j=0,NSuper
END DO !newton

! use Newton result if minimum is within parameter range, else see if supersampled
! initial guess is better than previous result
IF(MAXVAL(ABS(Xi2)).LE.1.) THEN ! use newton result
    ! calculate new distance
    xWinner=0.
    NormVecWinner=0.
    DO j=0,NSuper
      DO i=0,NSuper
        xWinner=xWinner+xBC_NSuper(:,i,j)*Lag_NSuper(1,i)*Lag_NSuper(2,j)
        NormVecWinner=NormVecWinner+NormVec_NSuper(:,i,j)*Lag_NSuper(1,i)*Lag_NSuper(2,j)
    END DO! i=0,NSuper
END DO! j=0,NSuper
Winner_Dist2=SUM((xWinner-hwmPhys)*(xWinner-hwmPhys))
END IF

NormVecWinner=NormVecWinner/NORM2(NormVecWinner)
F=(hwmPhys-xWinner)/NORM2(hwmPhys-xWinner)
IF(Winner_Dist2.LE.dist2RP)THEN
    SELECT CASE(locSideID)
    CASE(XI_MINUS)
      ! switch to right hand system
      hwmStd=(/-1.,Xi2(2),Xi2(1)/)
  CASE(ETA_MINUS)
      hwmStd=(/Xi2(1),-1.,Xi2(2)/)
  CASE(ZETA_MINUS)
      ! switch to right hand system
      hwmStd=(/Xi2(2),Xi2(1),-1./)
  CASE(XI_PLUS)
      hwmStd=(/1.,Xi2(1),Xi2(2)/)
  CASE(ETA_PLUS)
      ! switch to right hand system
      hwmStd=(/-Xi2(1),1.,Xi2(2)/)
  CASE(ZETA_PLUS)
      hwmStd=(/Xi2(1),Xi2(2),1./)
  END SELECT
  ! keep the RP in the not found list, find the side with the best angle and minimum distance
  dist2RP=Winner_Dist2
END IF
END DO! SideID=1,nBCSides
LOGWRITE(UNIT_logOut,'(A)')' done.'
LOGWRITE(UNIT_logOut,'(A,F15.8)')'  Max. distance: ',SQRT(dist2RP)
END IF!(.NOT.hwmFound)

IF(.NOT.hwmFound)THEN
    ! Only mark as invalid if greater then max tolerance
    IF(MAXVAL(ABS(hwmStd)).GT.maxTol)THEN
      ! h_wm has not been found
      WRITE(*,*) 'Exchange location with Coordinates ',hwmPhys, ' is a troublemaker!'
      CALL abort(__STAMP__, &
         'Newton has reached 50 Iter, Point not found')
  END IF
END IF

END SUBROUTINE PhysToStdCoords

!==================================================================================================================================
!> Evaluate solution at current time t at recordpoint positions and fill output buffer
!==================================================================================================================================
FUNCTION InterpolateHwm(iElem,L_xi,L_eta,L_zeta,comp,prim)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_DG_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER, INTENT(IN)            :: iElem                                      !< Mesh element containing the solution to be interpolated
INTEGER,INTENT(IN)             :: comp                                       !< Component of the Cons/Prim vector to be interpolated
REAL,INTENT(IN)                :: L_xi(0:PP_N),L_eta(0:PP_N),L_zeta(0:PP_NZ) !< Lagrange interpolating polynomials
LOGICAL,INTENT(IN),OPTIONAL    :: prim                                       !< .TRUE. so that comp refers to the components of Primitive var. vector
REAL                           :: InterpolateHwm
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: i,j,k
LOGICAL                        :: primitive
REAL                           :: tmpSum
REAL                           :: L_eta_zeta
REAL                           :: Var(0:PP_N,0:PP_N,0:PP_NZ)
!----------------------------------------------------------------------------------------------------------------------------------
primitive=.FALSE.

IF (PRESENT(prim)) primitive=prim

IF (primitive) THEN
    Var(:,:,:) = UPrim(comp,:,:,:,iElem)
ELSE
    Var(:,:,:) = U(comp,:,:,:,iElem)
END IF

tmpSum=0.
DO k=0,PP_NZ; DO j=0,PP_N
  L_eta_zeta=L_eta(j)*L_zeta(k)
  DO i=0,PP_N
    tmpSum=tmpSum + Var(i,j,k)*L_xi(i)*L_eta_zeta
  END DO !i
END DO; END DO !k

InterpolateHwm = tmpSum

END FUNCTION InterpolateHwm


!==================================================================================================================================
!> Evaluate u_tau from a log-law given instantaneous velocity, using Newton method
!==================================================================================================================================
FUNCTION NewtonLogLaw(velx,nu)
! MODULES
USE MOD_WMLES_Vars                  
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL, INTENT(IN)            :: velx         ! Tangential velocity to feed wall model (log-law function)
REAL, INTENT(IN)            :: nu           ! kinematic viscosity at h_wm
REAL                        :: NewtonLogLaw ! Computed friction velocity u_tau
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iter
REAL                        :: f,fprime
!----------------------------------------------------------------------------------------------------------------------------------

NewtonLogLaw = 1.0 ! Initial guess, usually not so bad for u_tau
iter = 0

SELECT CASE(WallModel)
    CASE(WMLES_LOGLAW)
        DO WHILE(iter.LT.10) ! Maximum number of iterations. Usually, u_tau is found with 1~4 iters
            iter = iter+1
            f = NewtonLogLaw*( (1./vKarman)*LOG(abs_h_wm*NewtonLogLaw/nu) + B ) - velx
            fprime = (1./vKarman) * (LOG(abs_h_wm*NewtonLogLaw/nu) + 1.) + B
            ! IF(ABS(fprime).LE.1.E-5) EXIT ! fprime ~ 0 -- INCOMING OVERFLOW BY DIVISION
            NewtonLogLaw = NewtonLogLaw - f/fprime
            IF (ABS(f/fprime).LE.1.E-5) EXIT ! 1.E-5 = stop criterion (tolerance)
        END DO

    CASE(WMLES_REICHARDT)

END SELECT


END FUNCTION NewtonLogLaw

!==================================================================================================================================
!> Evaluate u_tau from a log-law given instantaneous velocity, using Newton method
!==================================================================================================================================
FUNCTION NewtonCouette(velx,nu)
! MODULES
USE MOD_WMLES_Vars                  
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL, INTENT(IN)            :: velx         ! Tangential velocity to feed wall model (log-law function)
REAL, INTENT(IN)            :: nu           ! kinematic viscosity at h_wm
REAL                        :: NewtonCouette ! Computed friction velocity u_tau
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iter
REAL                        :: f,fprime
!----------------------------------------------------------------------------------------------------------------------------------

NewtonCouette = 0.5 ! Initial guess, usually not so bad for u_tau
iter = 0

DO WHILE(iter.LT.10) ! Maximum number of iterations. Usually, u_tau is found with 1~4 iters
    iter = iter+1
    f = velx + (1./(2.*nu))*NewtonCouette*abs_h_wm*(1.-abs_h_wm)
    fprime = (1./(2.*nu))*abs_h_wm*(1.-abs_h_wm)
    ! IF(ABS(fprime).LE.1.E-5) EXIT ! fprime ~ 0 -- INCOMING OVERFLOW BY DIVISION
    NewtonCouette = NewtonCouette - f/fprime
    IF (ABS(f/fprime).LE.1.E-5) EXIT ! 1.E-5 = stop criterion (tolerance)
END DO

END FUNCTION NewtonCouette

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
SDEALLOCATE(BCSideToWMLES)
SDEALLOCATE(WMLESToBCSide)
SDEALLOCATE(TauW_Proc)
SDEALLOCATE(Proc_RecvTauW)
SDEALLOCATE(Proc_RecvTauW_Inv)
SDEALLOCATE(Proc_SendTauW)
SDEALLOCATE(nTauW_MINE)
SDEALLOCATE(nTauW_YOURS)
SDEALLOCATE(TauW_MINE)
SDEALLOCATE(TauW_YOURS)
SDEALLOCATE(TauW_MINE_FacePoint)
SDEALLOCATE(TauW_MINE_InteriorPoint)
SDEALLOCATE(TauW_MINE_Interpolate)
SDEALLOCATE(TauW_MINE_NormVec)
SDEALLOCATE(nTauW_MINE_FacePoint)
SDEALLOCATE(nTauW_MINE_InteriorPoint)
SDEALLOCATE(nTauW_MINE_Interpolate)
SDEALLOCATE(FaceToLocalPoint)
SDEALLOCATE(InteriorToLocalPoint)
SDEALLOCATE(InterpToLocalPoint)
SDEALLOCATE(WMLES_RecvRequests)
SDEALLOCATE(WMLES_SendRequests)


END SUBROUTINE FinalizeWMLES



!==================================================================================================================================
END MODULE MOD_WMLES