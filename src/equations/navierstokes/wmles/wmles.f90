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
CALL prms%CreateStringOption( 'WMConnectionFile', "(relative) path to the file containing the wall-modeling connection information.")
CALL prms%CreateIntFromStringOption('WallModel', "Wall model to be used on walls defined with approximate boundary conditions", 'LogLaw')
CALL addStrListEntry('WallModel', 'Schumann', WMLES_SCHUMANN)
CALL addStrListEntry('WallModel', 'LogLaw', WMLES_LOGLAW)
CALL addStrListEntry('WallModel', 'WernerWangle', WMLES_WERNERWANGLE)
CALL addStrListEntry('WallModel', 'Reichardt', WMLES_REICHARDT)
CALL addStrListEntry('WallModel', 'Spalding', WMLES_SPALDING)
CALL addStrListEntry('WallModel', 'EquilibriumTBLE', WMLES_EQTBLE)
CALL addStrListEntry('WallModel', 'Couette', WMLES_COUETTE)

CALL prms%CreateRealOption('vKarman', "von Karman constant for the log-law model", "0.41")
CALL prms%CreateRealOption('B', "intercept coefficient for the log-law model", "5.2")

END SUBROUTINE DefineParametersWMLES

!==================================================================================================================================
!> This subroutine reads and initialize parameters of WMLES computation, and also populates all the necessary
!> mappings for the solution interpolation to the h_wm point and the communications between processors with
!> BC side and h_wm point.
!==================================================================================================================================
SUBROUTINE InitWMLES()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_HDF5_Input
USE MOD_WMLES_Vars
USE MOD_Mesh_Vars              ! ,ONLY: nBCSides, BC, BoundaryType, MeshInitIsDone, ElemToSide, SideToElem, SideToGlobalSide
USE MOD_Interpolation_Vars      ,ONLY: InterpolationInitIsDone,xGP,wBary,NodeType
USE MOD_Basis                   ,ONLY: LagrangeInterpolationPolys
USE MOD_ReadInTools             ,ONLY: GETINTFROMSTR, GETREAL, GETSTR
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
INTEGER                         :: SideID, iSide, iGlobalSide, LocSide, MineCnt, nSideIDs, offsetSideID
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

! NEW ONES
INTEGER                         :: iExt, nWMLESSides_global, N_WMConnection, nBCSides_global
INTEGER                         :: iElem, iProc, iWMLESSide
CHARACTER(LEN=255)              :: WMConnectionFileProposal, NodeType_WMConnection
INTEGER, ALLOCATABLE            :: BCSideToWMLES_global(:)
REAL, ALLOCATABLE               :: WMConnection(:,:,:,:), HWMInterpInfo_tmp(:,:)

#if USE_MPI
INTEGER, ALLOCATABLE            :: nBCSidesPerProc(:)
REAL, ALLOCATABLE               :: HWMSendInfo_tmp(:,:,:), HWMLocalInfo_tmp(:,:), HWMRecvInfo_tmp(:,:,:)
INTEGER, ALLOCATABLE            :: nHWMSendPoints_tmp(:), nHWMRecvPoints_tmp(:)
INTEGER, ALLOCATABLE            :: SendToInterpPoint_tmp(:,:), LocalToInterpPoint_tmp(:)
LOGICAL, ALLOCATABLE            :: WMLESRecvFromProc(:), WMLESSendToProc(:)
#endif
!==================================================================================================================================
! IF((.NOT.InterpolationInitIsDone).OR.(.NOT.MeshInitIsDone).OR.WMLESInitDone) THEN
!  CALL CollectiveStop(__STAMP__,&
!    'Wall-Modeled LES not ready to be called or already called.')
! END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT Wall-Modeled LES...'


!> Read the parameters that are independent from the WM connection

WallModel = GETINTFROMSTR('WallModel')

SELECT CASE(WallModel)
    CASE(WMLES_LAMINAR_INTEGRAL)
        nHWMPropSend = 4 ! velocity vector (3 components) and dp/dx

    CASE DEFAULT
        nHWMPropSend = 3 ! velocity vector only (3 components)
END SELECT

! If Schumann's model is selected, check if we are in a channel flow
IF ((WallModel .EQ. WMLES_SCHUMANN) .AND. (.NOT.STRICMP(Testcase, Channel))) THEN
    CALL CollectiveStop(__STAMP__,&
        "Schumann's wall model can only be applied to the Channel flow testcase.")
END IF

! Log-law model parameters
vKarman = GETREAL('vKarman')
B = GETREAL('B')

! =-=-=-=--=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! We now set up the necessary information for the wall-stress computation and 
! message communication between the relevant processors for the procedure.
! In other words, we do the following:
! (0) Bear in mind that all info in the WM connection file is in terms of GLOBAL elements and BC side numbers)
! 1) Get the total number of BC sides (global), so that each MPI proc. may scan the whole array in the HDF5 file
! 2) Check wheter a Modelled Side is within this proc. (using the mapping in the HDF5 file)
! 3) For all modeled sides within this proc, we store wall-tangent vector info, as well as the
!    processors (and count) that need to send us information for the wall-stress calculation
! 4) Looping over all modeled sides in the HDF5 file, check whether there is any element in this proc.
!    responsible for sending information to a proc. responsible for modeled BC imposition. In such a case,
!    store the relevant information about the processor to receive your information, as well as the h_wm
!    point (in standard coordinates, \xi, \eta, \zeta) where the solution within your element must be interpolated.
! =-=-=-=--=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! Read WM connection file

iExt=INDEX(MeshFile,'.',BACK = .TRUE.) ! Position of file extension
WMConnectionFileProposal = MeshFile(:iExt-6)
WMConnectionFileProposal = TRIM(WMConnectionFileProposal) // '_WM.h5'

WMConnectionFile = GETSTR('WMConnectionFile', WMConnectionFileProposal)
IF (.NOT. FILEEXISTS(WMConnectionFile)) THEN
    CALL Abort(__STAMP__, 'WM Connection File not found! Check if the parameter WmConnnectionFile' &
                            // 'is correctly set in the parameter file, or if the filename is the default one (_WM.h5)')
END IF

CALL OpenDataFile(TRIM(WMConnectionFile),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadAttribute(File_ID,'nModelledBCSides',1,IntScalar=nWMLESSides_global)

! Safety checks if interface preparation information is compatible with current simulation
CALL ReadAttribute(File_ID,'N',1,IntScalar=N_WMConnection)
IF (N_WMConnection.NE.PP_N) &
    CALL Abort(__STAMP__,'Polynomial degree of calculation not equal to WM connection!')
CALL ReadAttribute(File_ID,'NodeType',1,StrScalar=NodeType_WMConnection)
IF (.NOT.STRICMP(NodeType,TRIM(NodeType_WMConnection))) &
        CALL Abort(__STAMP__,'Node type of calculation not equal to WM connection!')

! Get global number of BC Sides
#if USE_MPI
CALL MPI_ALLREDUCE(nBCSides,nBCSides_global,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,iError)
#else
nBCSides_global =  nBCSides
#endif

! Read global mapping from BCside to WMSide
ALLOCATE(BCSideToWMLES_global(nBCSides_global))
CALL ReadArray('MappingWM',1,(/nBCSides_global/),0,1,IntArray=BCSideToWMLES_global)

! Read connection information
ALLOCATE(WMConnection(N_INTERFACE_PARAMS,0:PP_N,0:PP_N,nWMLESSides_global))
CALL ReadArray('WM',4,&
               (/N_INTERFACE_PARAMS,PP_N+1,PP_N+1,nWMLESSides_global/),&
               0,4,RealArray=WMConnection)

CALL CloseDataFile()

! We now count the number of WM sides in this proc. (i.e., local info)
ALLOCATE(BCSideToWMLES(nBCSides))
ALLOCATE(WMLESToBCSide_tmp(nBCSides))

BCSideToWMLES = 0
nWMLESSides = 0
WMLESToBCSide_tmp = 0
DO iSide=1,nBCSides
    IF(BoundaryType(BC(iSide),BC_TYPE).EQ.5) THEN ! WMLES BC Side
        nWMLESSides = nWMLESSides + 1
        BCSideToWMLES(iSide) = nWMLESSides
        WMLESToBCSide_tmp(nWMLESSides) = iSide
    END IF
END DO

! Allocate the final array, with the right number of elements, and free the excess memory
ALLOCATE(WMLESToBCSide(nWMLESSides))
DO iSide=1,nWMLESSides
    WMLESToBCSide(iSide) = WMLESToBCSide_tmp(iSide)
END DO
DEALLOCATE(WMLESToBCSide_tmp)

! Sanity check
#if USE_MPI
CALL MPI_REDUCE(nWMLESSides,TotalNWMLESSides,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iError)
IF (myRank.EQ.0) THEN
    IF (TotalNWMLESSides .NE. nWMLESSides_global) CALL Abort(__STAMP__, &
        'Total number of WMLES Sides does not match the total modeled sides in the connection file!')
END IF
#endif /* USE_MPI */

! We now loop over all modeled sides from the WM Connection file, and check the sides that correspond
! to any of the WMLES side in this proc., so that we can store information about receiving communication
! In order to do so, we use a trick available due to the way the mesh is partitioned.
! From the partitioning method, we know that each consecutive MPI process holds consecutive
! BC sides (because they are sorted in increasing order in the mesh process). 
! Hence, we gather information about the range of BC sides for each of the current MPI procs., and store
! the offset values (in global terms) where the BC sides in this proc. starts.
! Then, we check the global side number of the modeled side (WM Conn file) against these ranges (given by
! the offset) to see if this proc. owns the side locally. If so, we store receiving information

#if USE_MPI
ALLOCATE(nBCSidesPerProc(0:nProcessors-1))
nBCSidesPerProc = 0
nBCSidesPerProc(myRank) = nBCSides
CALL MPI_ALLREDUCE(MPI_IN_PLACE,nBCSidesPerProc,nProcessors,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,iError)
offsetBCSides = SUM(nBCSidesPerProc(0:myRank-1))
DEALLOCATE(nBCSidesPerProc)
#else
! No MPI means no offset for the sides
offsetBCSides = 0
#endif /* USE_MPI */

ALLOCATE(HWMSendInfo_tmp(nHWMPropSend+4, nWMLESSides_global*(PP_N+1)*(PP_N+1), 0:nProcessors-1))
ALLOCATE(HWMRecvInfo_tmp(4, nWMLESSides_global*(PP_N+1)*(PP_N+1), 0:nProcessors-1))
ALLOCATE(HWMInterpInfo_tmp(4, nWMLESSides_global*(PP_N+1)*(PP_N+1)))
ALLOCATE(nHWMSendPoints_tmp(0:nProcessors-1))
ALLOCATE(nHWMRecvPoints_tmp(0:nProcessors-1))
ALLOCATE(HWMLocalInfo_tmp(4, nWMLESSides_global*(PP_N+1)*(PP_N+1)))
ALLOCATE(LocalToInterpPoint_tmp(nWMLESSides_global*(PP_N+1)*(PP_N+1)))
ALLOCATE(SendToInterpPoint_tmp(nWMLESSides_global*(PP_N+1)*(PP_N+1), 0:nProcessors-1))
ALLOCATE(WMLESRecvFromProc(0:nProcessors-1))
ALLOCATE(WMLESSendToProc(0:nProcessors-1))

! Initializing the allocated memory values

WMLESRecvFromProc = .FALSE. ! Flag to check whether we receive anything from this proc.
WMLESSendToProc = .FALSE. ! Flag to check whether we send anything to this proc.
HWMSendInfo_tmp = 0.
HWMRecvInfo_tmp = 0
HWMInterpInfo_tmp = 0.
nHWMSendPoints_tmp = 0
nHWMRecvPoints_tmp = 0
HWMLocalInfo_tmp = 0
LocalToInterpPoint_tmp = 0
SendToInterpPoint_tmp = 0
nHWMLocalPoints = 0
nHWMInterpPoints = 0
nWMLESRecvProcs = 0
nWMLESSendProcs = 0


DO iGlobalSide=1, nWMLESSides_global
    IF ((iGlobalSide.GT.offsetBCSides) .AND. (iGlobalSide.LE.(offsetBCSides+nBCSides))) THEN ! Local BC side!
        ! Sanity check
        IF (BCSideToWMLES(iGlobalSide-offsetBCSides) .EQ. 0) CALL Abort(__STAMP__, &
            'Modeled local BC found in WM Connection file, but it has not been marked as a WMLES Side using the current config./mesh!')

        iWMLESSide = BCSideToWMLES(iGlobalSide-offsetBCSides)

        ! We now check if we need to receive h_wm point info from another proc., or if h_wm is within a local element
        ! This check is done for each GL point within this BC side
        DO p=0,PP_N; DO q=0,PP_N

            iElem = INT(WMConnection(INTERFACE_ELEMSEND,p,q,iGlobalSide)) ! Global ElemID
            IF ((iElem.GT.offsetElem) .AND. (iElem.LE.(offsetElem+nElems))) THEN
                ! Element where h_wm lies is local to this proc
                ! Hence, no communication will be needed. 
                ! We count this h_wm point as "interpolation needed", and store standard coordinate 
                ! information for later interpolation
                ! We also store the BC global side and the p,q, point that this interpolation corresponds to

                nHWMLocalPoints = nHWMLocalPoints + 1
                nHWMInterpPoints = nHWMInterpPoints + 1

                LocalToInterpPoint_tmp(nHWMLocalPoints) = nHWMInterpPoints

                HWMLocalInfo_tmp(1:4, nHWMLocalPoints) = (/WMConnection(INTERFACE_L,p,q,iGlobalSide), REAL(p), REAL(q), REAL(iWMLESSide)/)
                HWMInterpInfo_tmp(:,nHWMInterpPoints) = (/REAL(iElem-offsetElem), & ! All elements in the array construction must be of the same type.
                                                        WMConnection(INTERFACE_XI, p, q, iGlobalSide), &
                                                        WMConnection(INTERFACE_ETA,p, q, iGlobalSide), &
                                                        WMConnection(INTERFACE_ZETA,p,q, iGlobalSide)/)

            ELSE
                ! Element where h_wm lies is in another proc.
                ! Hence, we set up the communication info to receive h_wm later

                ! Get the proc. rank responsible for the element containing h_wm
                iProc = ELEMIPROC(iElem)

                nHWMRecvPoints_tmp(iProc) = nHWMRecvPoints_tmp(iProc) + 1 ! Number of flow props to receive from this proc.
                IF (.NOT.WMLESRecvFromProc(iProc)) THEN ! Check if we've come across this receiving rank yet
                    nWMLESRecvProcs = nWMLESRecvProcs + 1 ! Sum up the number of procs to receive info from
                    !WMLESRecvProc_tmp(nWMLESRecvProcs) = iProc
                    WMLESRecvFromProc(iProc) = .TRUE.
                END IF

                HWMRecvInfo_tmp(1:4, nHWMRecvPoints_tmp(iProc), iProc) = (/WMConnection(INTERFACE_L, p, q, iGlobalSide), REAL(p), REAL(q), REAL(iWMLESSide)/)

            END IF ! IF ELSE element containing h_wm is local to this proc.
        END DO; END DO ! p, q

    ELSE 
        ! This global WMLES Side is not within this proc.
        ! However, this proc. may need to SEND h_wm info to the rank responsible for the BC imposition on this side

        DO p=0,PP_N; DO q=0,PP_N

            iElem = INT(WMConnection(INTERFACE_ELEMSEND,p,q,iGlobalSide)) ! Global ElemID
            IF ((iElem.GT.offsetElem) .AND. (iElem.LE.(offsetElem+nElems))) THEN
                ! This proc. is not responsible for WMLES BC imposition, but it contains
                ! the info in h_wm point.
                ! Hence, we need to SEND this info to the proc. responsible for BC imposition

                iProc = ELEMIPROC(INT(WMConnection(INTERFACE_ELEMRECV,p,q,iGlobalSide)))

                nHWMSendPoints_tmp(iProc) = nHWMSendPoints_tmp(iProc) + 1 ! Number of flow props to send to the responsible rank
                IF (.NOT.WMLESSendToProc(iProc)) THEN
                    nWMLESSendProcs = nWMLESSendProcs + 1 ! Sum up the number of procs to send the info
                    !WMLESSendProc_tmp(nWMLESSendProcs) = iProc
                    WMLESSendToProc(iProc) = .TRUE.
                END IF

                ! We now store the information needed for sending flow props at h_wm location
                ! Besides sending the flow properties needed for the model, we also send
                ! the corresponding BC Global Side and the p,q point where the receiving rank
                ! should use these flow properties 
                DO i=1,nHWMPropSend
                    HWMSendInfo_tmp(i,nHWMSendPoints_tmp(iProc),iProc) = -10.
                END DO

                ! IMPORTANT NOTE:
                ! There are two ways to design the send-receiving communication, in terms of the information
                ! contained in each message.
                ! The first one is to make the info about which side/point corresponds to the flow property
                ! being sent as part of the message itself. This would be the general (and least bugged)
                ! way to do it.
                ! The second way is to create a mapping from the receiving info (flow properties) to the 
                ! corresponding modeled side/point. Depending on the way the parallel framework is designed,
                ! this may be rather tricky. However, in this case, it is not. Since all MPI procs. are reading info 
                ! from the WM Connection file in the same order (i.e., from 1:N global side, and for each side, from 
                ! p=0, q=0,...,Np, p=1, q=0,...,Np etc), then, we are guaranteed that whenever an MPI proc. identify 
                ! a side/point that needs to receive information, this side/point will be identified by the 
                ! counterpart MPI proc (the sending proc.) at the same iteration step through the WM connection info. 
                ! Hence, we are sure that, e.g., the first receiving info from MPI proc. #X corresponds to the first 
                ! point we identified as being owned by MPI proc. #X, in the receiving counterpart.
                ! Anyway, the first approach is being used here since it prevents bugs from this "automatic" and
                ! fortunate mapping (due to the way the parallel domain decomposition is done and how we read the info)
                
                HWMSendInfo_tmp(nHWMPropSend+1:nHWMPropSend+4,nHWMSendPoints_tmp(iProc),iProc) = (/WMConnection(INTERFACE_L,p,q,iGlobalSide), REAL(p), REAL(q), REAL(iGlobalSide)/)

                nHWMInterpPoints = nHWMInterpPoints + 1
                HWMInterpInfo_tmp(:,nHWMInterpPoints) = (/REAL(iElem-offsetElem), &
                                                        WMConnection(INTERFACE_XI, p, q, iGlobalSide), &
                                                        WMConnection(INTERFACE_ETA,p, q, iGlobalSide), &
                                                        WMConnection(INTERFACE_ZETA,p,q, iGlobalSide)/)

                SendToInterpPoint_tmp(nHWMSendPoints_tmp(iProc), iProc) = nHWMInterpPoints

            END IF ! IF element containing h_wm is local to this proc.

        END DO; END DO ! p, q
    END IF ! IF ELSE Local BC Side
END DO

! In order to reduce memory usage, we re-map the SendToInterp vector, ordering Send Points from the
! lowest to the highest MPI rank, sequentially. Therefore, we also need to re-map HWMInterpInfo.
ALLOCATE(nHWMSendPoints(0:nWMLESSendProcs))
ALLOCATE(SendToInterpPoint(SUM(nHWMSendPoints_tmp)))
ALLOCATE(WMLESSendProc(nWMLESSendProcs))
ALLOCATE(HWMSendInfo(nHWMPropSend+4, SUM(nHWMSendPoints_tmp)))
ALLOCATE(WMLESSendRange(2,nWMLESSendProcs))

nHWMSendPoints(0) = 0
nWMLESSendProcs = 0 ! We will recount this
DO iProc=0, nProcessors-1
    IF (WMLESSendToProc(iProc)) THEN
        nWMLESSendProcs = nWMLESSendProcs + 1
        WMLESSendProc(nWMLESSendProcs) = iProc
        nHWMSendPoints(nWMLESSendProcs) = nHWMSendPoints_tmp(iProc)

        WMLESSendRange(1,nWMLESSendProcs) = SUM(nHWMSendPoints(:nWMLESSendProcs-1))+1
        WMLESSendRange(2,nWMLESSendProcs) = SUM(nHWMSendPoints(:nWMLESSendProcs))

        HWMSendInfo(:, WMLESSendRange(1,nWMLESSendProcs):WMLESSendRange(2,nWMLESSendProcs)) = HWMSendInfo_tmp(:, 1:nHWMSendPoints_tmp(iProc), iProc)

        SendToInterpPoint(WMLESSendRange(1,nWMLESSendProcs):WMLESSendRange(2,nWMLESSendProcs)) = SendToInterpPoint_tmp(1:nHWMSendPoints_tmp(iProc), iProc)
    END IF
END DO

! The same logic above is done for the receiving buffer, just so that we store the receiving ranges
! according to each MPI rank to receive from
ALLOCATE(WMLESRecvProc(nWMLESRecvProcs))
ALLOCATE(nHWMRecvPoints(0:nWMLESRecvProcs))
ALLOCATE(HWMRecvInfo(nHWMPropSend+4, SUM(nHWMRecvPoints_tmp)))
ALLOCATE(WMLESRecvRange(2,nWMLESRecvProcs))

nWMLESRecvProcs = 0 ! We will recount this
nHWMRecvPoints(0) = 0
DO iProc=0, nProcessors-1
    IF (WMLESRecvFromProc(iProc)) THEN
        nWMLESRecvProcs = nWMLESRecvProcs + 1
        WMLESRecvProc(nWMLESRecvProcs) = iProc
        nHWMRecvPoints(nWMLESRecvProcs) = nHWMRecvPoints_tmp(iProc)

        WMLESRecvRange(1,nWMLESRecvProcs) = SUM(nHWMRecvPoints(:nWMLESRecvProcs-1))+1
        WMLESRecvRange(2,nWMLESRecvProcs) = SUM(nHWMRecvPoints(:nWMLESRecvProcs))

        HWMRecvInfo(nHWMPropSend+1:nHWMPropSend+4, WMLESRecvRange(1,nWMLESRecvProcs):WMLESRecvRange(2,nWMLESRecvProcs)) = HWMRecvInfo_tmp(1:4, 1:nHWMRecvPoints_tmp(iProc), iProc)
    END IF
END DO

! Memory bookeeping (allocate the minimum necessary memory, and free the excess from the temp variables)
ALLOCATE(HWMInterpInfo(4, nHWMInterpPoints))
ALLOCATE(LocalToInterpPoint(nHWMLocalPoints))
ALLOCATE(HWMLocalInfo(4, nHWMLocalPoints))

DO i=1,nHWMInterpPoints
    HWMInterpInfo(1:4, i) = HWMInterpInfo_tmp(1:4, i)
END DO

DO i=1,nHWMLocalPoints
    LocalToInterpPoint(i) = LocalToInterpPoint_tmp(i)
    HWMLocalInfo(1:4, i) = HWMLocalInfo_tmp(1:4, i)
END DO


!> =-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=--=-
!> Logging Information (mainly for debugging purposes)
!> =-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=--=-

LOGWRITE(*,'(20("=-"))') 
LOGWRITE(*,*) "WMLES Connection Information"
LOGWRITE(*,'(20("=-"))')
LOGWRITE(*,'(A55, I6)') 'Number of WMLESSides: ', nWMLESSides
LOGWRITE(*,'(A55, I6)') 'Number of Procs to Send Info: ', nWMLESSendProcs
LOGWRITE(*,'(A55, I6)') 'Number of Procs to Receive Info: ', nWMLESRecvProcs
LOGWRITE(*,'(A55, I6)') 'Number of h_wm Sending Points: ', SUM(nHWMSendPoints)
LOGWRITE(*,'(A55, I6)') 'Number of h_wm Receiving Points: ', SUM(nHWMRecvPoints)
LOGWRITE(*,'(A55, I6)') 'Number of h_wm Local Points: ', nHWMLocalPoints
LOGWRITE(*,'(A55, I6)') 'Number of h_wm Interpolation Points: ', nHWMInterpPoints
LOGWRITE(*,'(20("--"))')
LOGWRITE(*,'(A30)') 'Local Information/Mapping'
LOGWRITE(*,'(A10,1X, 2(A2,1X), 3(A10,1X), A5,1X, 3(A10,1X))') "localPoint", "p", "q", "WMLESBCSide", "globalSide", "interpPoint", "iElem", "xi", "eta", "zeta"
DO i=1,nHWMLocalPoints
    LOGWRITE(*,'(I10,1X, 2(I2,1X), 3(I10,1X), I5,1X, 3(E10.4,1X))') &
    i, INT(HWMLocalInfo(2, i)), INT(HWMLocalInfo(3, i)), INT(HWMLocalInfo(4, i)), WMLESToBCSide(INT(HWMLocalInfo(4, i)))+offsetBCSides, &
    LocalToInterpPoint(i), INT(HWMInterpInfo(1, LocalToInterpPoint(i))), HWMInterpInfo(2, LocalToInterpPoint(i)), HWMInterpInfo(3, LocalToInterpPoint(i)), HWMInterpInfo(4, LocalToInterpPoint(i))
END DO
LOGWRITE(*,'(20("--"))')
LOGWRITE(*,'(A30)') 'Sending Information/Mapping'
LOGWRITE(*,'(2(A5,2X), A10,2X, 2(A2,1X), 2(A10,1X), A5,1X, 3(A10,1X))') "iProc", "nProc", "sendPoint", "p", "q", "globalSide", "interpPoint", "iElem", "xi", "eta", "zeta"
DO i=1,nWMLESSendProcs
    DO j=WMLESSendRange(1,i), WMLESSendRange(2,i)
    LOGWRITE(*,'(2(I5,2X), I10,2X, 2(I2,1X), 2(I10,1X), I5,1X, 3(E10.4,1X))') &
            i, WMLESSendProc(i), j, INT(HWMSendInfo(nHWMPropSend+2, j)), INT(HWMSendInfo(nHWMPropSend+3, j)), INT(HWMSendInfo(nHWMPropSend+4, j)), &
            SendToInterpPoint(j), INT(HWMInterpInfo(1, SendToInterpPoint(j))), HWMInterpInfo(2, SendToInterpPoint(j)), HWMInterpInfo(3, SendToInterpPoint(j)), HWMInterpInfo(4, SendToInterpPoint(j))
    END DO
END DO
LOGWRITE(*,'(A30)') 'Receiving Information/Mapping'
LOGWRITE(*,'(2(A5,2X), A10,2X, 2(A2,1X), A10,1X)') "iProc", "nProc", "sendPoint", "p", "q", "globalSide"
DO i=1,nWMLESRecvProcs
    DO j=WMLESRecvRange(1,i), WMLESRecvRange(2,i)
    LOGWRITE(*,'(2(I5,2X), I10,2X, 2(I2,1X), I10,1X)') &
            i, WMLESRecvProc(i), j, INT(HWMRecvInfo(nHWMPropSend+2, j)), INT(HWMRecvInfo(nHWMPropSend+3, j)), INT(HWMRecvInfo(nHWMPropSend+4, j))+offsetBCSides
    END DO
END DO
IF (LOGGING) CALL FLUSH(UNIT_logOut)

SDEALLOCATE(HWMSendInfo_tmp)
SDEALLOCATE(HWMRecvInfo_tmp)
SDEALLOCATE(HWMInterpInfo_tmp)
SDEALLOCATE(nHWMSendPoints_tmp)
SDEALLOCATE(nHWMRecvPoints_tmp)
SDEALLOCATE(HWMLocalInfo_tmp)
SDEALLOCATE(LocalToInterpPoint_tmp)
SDEALLOCATE(SendToInterpPoint_tmp)
SDEALLOCATE(WMLESRecvFromProc)
SDEALLOCATE(WMLESSendToProc)

ALLOCATE(WMLES_TauW(2, 0:PP_N, 0:PP_NZ, nWMLESSides))
ALLOCATE(HWMInfo(nHWMPropSend+1, 0:PP_N, 0:PP_NZ, nWMLESSides))
ALLOCATE(WMLES_RecvRequests(nWMLESRecvProcs))
ALLOCATE(WMLES_SendRequests(nWMLESSendProcs))

! Sanity check:
! The number of WMLES points (WMLESSides * (PP_N+1) * (PP_NZ+1)) in each MPI proc. must equal the
! sum of the Local plus the Receiving h_wm points
IF ((nWMLESSides*(PP_N+1)*(PP_NZ+1)) .NE. (nHWMLocalPoints+SUM(nHWMRecvPoints))) &
    CALL Abort(__STAMP__, "The number of wall-modeled points for this proc. is not equal" &
                // "to the number of local + receiving h_wm points! Inconsistency!", IntInfo=myRank)

! Just a checkpoint to make sure no MPI proc. leaves the subroutine before any other, so that
! shared info (buffers) are not erased (well, in this new implementation this is 
! not strictly necessary, but let us keep it this way.)
CALL MPI_BARRIER(MPI_COMM_FLEXI, iError)

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
!> Compute the wall stress, tau_w, at each point in a WMLES BC surface.
!==================================================================================================================================
SUBROUTINE ComputeWallStress()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Basis
USE MOD_WMLES_Vars
USE MOD_Mesh_Vars
USE MOD_Interpolation_Vars
USE MOD_DG_Vars                     
USE MOD_EOS_Vars                    ,ONLY: mu0
IMPLICIT NONE

!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                             :: sProc, FPInd, IPInd, IntPInd
INTEGER                             :: p,q,r,ElemID
REAL                                :: u_tau, u_mean, tau_w_mag, utang, VelMag
REAL                                :: vel_inst(3), rho_inst, tangvec(3), tau_w_vec(3)
INTEGER                             :: SurfLowUp ! Debug

INTEGER                             :: i, k, SideID
INTEGER                             :: iSendPoint, iLocalPoint
REAL                                :: xi, eta, zeta
REAL                                :: L_xi(0:PP_N), L_eta(0:PP_N), L_zeta(0:PP_N)
!==================================================================================================================================
IF (.NOT.WMLESInitDone) THEN
    CALL CollectiveStop(__STAMP__,&
    'Cannot compute Wall Stress -- Wall-Modeled LES INIT is not done.')
END IF

! TODO ?
! For algebraic models (except for Schumann), the outcome is really u_tau instead of tau_w directly
! Hence, to compute tau_w, we take tau_w = u_tau^2 * rho. However, this rho should be 
! the density AT THE WALL, and we are currently assuming the density at the exchange location, h_wm,
! which is conceptually wrong for a compressible flow.
! Thus, maybe we should send u_tau instead of tau_w to the MPI proc. responsible for BC imposition,
! and then this MPI proc. computes tau_w.
! (See Bocquet, Sagaut & Jouhaud).

! Start non-blockingly receiving h_wm info (for those who have anything to receive)
CALL StartReceiveHwmMPIData()

! First, we interpolate all h_wm information in the elements that this proc. owns, and that
! needs to be sent to another proc. (for BC imposition and calculation)
DO iSendPoint=1, WMLESSendRange(2,nWMLESSendProcs) ! loop from 1 to the last send point
    xi = HWMInterpInfo(2, SendToInterpPoint(iSendPoint))
    eta = HWMInterpInfo(3, SendToInterpPoint(iSendPoint))
    zeta = HWMInterpInfo(4, SendToInterpPoint(iSendPoint))

    ! Evaluate all the 1D Lagrange polynomials at each specific standard coordinate
    CALL LagrangeInterpolationPolys(xi,  PP_N, xGP, wBary, L_xi)
    CALL LagrangeInterpolationPolys(eta, PP_N, xGP, wBary, L_eta)
    CALL LagrangeInterpolationPolys(zeta,PP_N, xGP, wBary, L_zeta)

    ! Interpolate the flow properties at h_wm point
    DO k=1,3 ! Loop for velocity components
        HWMSendInfo(k, iSendPoint) = InterpolateHwm(INT(HWMInterpInfo(1, SendToInterpPoint(iSendPoint))), & ! iElem
                                        L_xi, L_eta, L_zeta, k+1, prim=.TRUE.)
    END DO
    ! Adjust for Laminar
    !IF (nHWMPropSend.GT.3) HWMSendInfo(4, i) = InterpolateHwm(INT(HWMInterpInfo(1, SendToInterpPoint(i))), & ! iElem
    !                                    L_xi, L_eta, L_zeta, 5, prim=.TRUE.)
END DO

! Send the information to the appropriate MPI procs.
CALL StartSendHwmMPIData()

! We now interpolate every h_wm point local to this proc. and arrange this info directly
! into the HWMInfo array
DO iLocalPoint=1, nHWMLocalPoints
    xi = HWMInterpInfo(2, LocalToInterpPoint(iLocalPoint))
    eta = HWMInterpInfo(3, LocalToInterpPoint(iLocalPoint))
    zeta = HWMInterpInfo(4, LocalToInterpPoint(iLocalPoint))

    ! Evaluate all the 1D Lagrange polynomials at each specific standard coordinate
    CALL LagrangeInterpolationPolys(xi, PP_N, xGP, wBary, L_xi)
    CALL LagrangeInterpolationPolys(eta, PP_N, xGP, wBary, L_eta)
    CALL LagrangeInterpolationPolys(zeta, PP_N, xGP, wBary, L_zeta)

    ! Interpolate the flow properties at h_wm point
    HWMInfo(1, INT(HWMLocalInfo(2, iLocalPoint)), & ! p
                INT(HWMLocalInfo(3, iLocalPoint)), & ! q
                INT(HWMLocalInfo(4, iLocalPoint))) & ! iWMLESSide
            = HWMLocalInfo(1, iLocalPoint) ! h_wm value

    DO k=2,4 ! Loop for velocity components
        HWMInfo(k, INT(HWMLocalInfo(2, iLocalPoint)), & ! p
                INT(HWMLocalInfo(3, iLocalPoint)), & ! q
                INT(HWMLocalInfo(4, iLocalPoint))) & ! iWMLESSide
            = InterpolateHwm(INT(HWMInterpInfo(1, LocalToInterpPoint(iLocalPoint))), & ! iElem
                            L_xi, L_eta, L_zeta, k, prim=.TRUE.)
    END DO
    ! HWMInfo(2, INT(HWMLocalInfo(2, iLocalPoint)), & ! p
    !             INT(HWMLocalInfo(3, iLocalPoint)), & ! q
    !             INT(HWMLocalInfo(4, iLocalPoint))) & ! iWMLESSide
    !         = InterpolateHwm(INT(HWMInterpInfo(1, LocalToInterpPoint(iLocalPoint))), & ! iElem
    !                         L_xi, L_eta, L_zeta, 2, prim=.TRUE.)
    ! HWMInfo(3, INT(HWMLocalInfo(2, iLocalPoint)), & ! p
    !             INT(HWMLocalInfo(3, iLocalPoint)), & ! q
    !             INT(HWMLocalInfo(4, iLocalPoint))) & ! iWMLESSide
    !         = InterpolateHwm(INT(HWMInterpInfo(1, LocalToInterpPoint(iLocalPoint))), & ! iElem
    !                         L_xi, L_eta, L_zeta, 3, prim=.TRUE.)
    ! HWMInfo(4, INT(HWMLocalInfo(2, iLocalPoint)), & ! p
    !             INT(HWMLocalInfo(3, iLocalPoint)), & ! q
    !             INT(HWMLocalInfo(4, iLocalPoint))) & ! iWMLESSide
    !         = InterpolateHwm(INT(HWMInterpInfo(1, LocalToInterpPoint(iLocalPoint))), & ! iElem
    !                         L_xi, L_eta, L_zeta, 4, prim=.TRUE.)
END DO

! Finish receiving h_wm info from other procs.
CALL FinishExchangeHWMData()

! Finally, with all h_wm info properly arranged into HWMInfo array, we can
! calculate the wall stress according to the chosen model

SELECT CASE(WallModel)

  CASE (WMLES_LOGLAW)

    LOGWRITE(*,'(30("=-"))')
    LOGWRITE(*,'(A60)') 'WMLES_TauW Calculation Logging/Debug Information'
    LOGWRITE(*,'(30("=-"))')
    LOGWRITE(*,'(2(A9,1X), A7,1X, 2(A2,1X), 3(A10,1X), 3(A10,1X), 2(A10,1X))') 'globalSide', 'WMLESSide', 'BC Side', 'p', 'q', 'HWMInfo_u', 'HWMInfo_v', 'HWMInfo_w', 'utang', 'h_wm', 'u_tau', 'TauW_x', 'TauW_y'
    DO SideID=1,nWMLESSides
        DO p=0,PP_N; DO q=0,PP_N

            ! Project the velocity vector onto the normal vector and subtract this projection (in the
            ! wall-normal direction) from the velocity. The result is the velocity aligned with the
            ! tangential direction
            tangvec = HWMInfo(2:4,p,q,SideID) - DOT_PRODUCT(HWMInfo(2:4,p,q,SideID),NormVec(1:3,p,q,0,WMLESToBCSide(SideID)))*NormVec(1:3,p,q,0,WMLESToBCSide(SideID))

            VelMag = 0.
                DO i=1,3
                    VelMag = VelMag + tangvec(i)**2
                END DO
            VelMag = SQRT(VelMag)
            tangvec = tangvec/VelMag ! Unit tangential vector

            !<-=-=-=-=- DEBUG; REMOVE LATER IF IT DOES NOT WORK =-=-=-=-=>!
            ! Force tangvec to be in the x dir
            tangvec = (/1.,0.,0./)
            !<-=-=-=-=- END OF DEBUG; REMOVE LATER IF IT DOES NOT WORK =-=-=-=-=>!
            utang = DOT_PRODUCT(HWMInfo(2:4,p,q,SideID),tangvec)

            u_tau = NewtonLogLaw(utang, (mu0/UPrim_master(1,p,q,WMLESToBCSide(SideID))), HWMInfo(1,p,q,SideID)) ! rho_wall is used here (through UPrim_master). Read TODO above
            tau_w_mag = UPrim_master(1,p,q,WMLESToBCSide(SideID))*(u_tau**2) ! CHECK TODO ABOVE
            tau_w_vec = tau_w_mag*tangvec

            !<-=-=-=-=- DEBUG; REMOVE LATER IF IT DOES NOT WORK =-=-=-=-=>!
            ! Create vector aligned with wall-normal velocity with tau_xy
            tau_w_vec = (/0.,tau_w_mag,0./)
            ! Project it onto the normal direction so that the correct signal is imposed
            WMLES_TauW(1,p,q,SideID) = -1.*DOT_PRODUCT(tau_w_vec(1:3),NormVec(1:3,p,q,0,WMLESToBCSide(SideID)))
            !<-=-=-=-=- END OF DEBUG; REMOVE LATER IF IT DOES NOT WORK =-=-=-=-=>!
            WMLES_TauW(2,p,q,SideID) = 0.

            LOGWRITE(*,'(2(I9,1X), I7,1X, 2(I2,1X), 3(E10.4,1X), 3(E10.4,1X), 2(E10.4,1X))') &
                    WMLESToBCSide(SideID)+offsetBCSides, SideID, WMLESToBCSide(SideID), p, q, HWMInfo(2,p,q,SideID), HWMInfo(3,p,q,SideID), HWMInfo(4,p,q,SideID), &
                    utang, HWMInfo(1,p,q,SideID), u_tau, WMLES_TauW(1,p,q,SideID), WMLES_TauW(2,p,q,SideID)
        END DO; END DO ! p,q
    END DO ! iSide

  CASE DEFAULT
    CALL Abort(__STAMP__,&
         'Unknown definition of Wall Model.')

END SELECT


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
REAL,INTENT(IN)                :: L_xi(0:PP_N),L_eta(0:PP_N),L_zeta(0:PP_N)  !< Lagrange interpolating polynomials
LOGICAL,INTENT(IN),OPTIONAL    :: prim                                       !< .TRUE. so that comp refers to the components of Primitive var. vector
REAL                           :: InterpolateHwm
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: i,j,k
LOGICAL                        :: primitive
REAL                           :: tmpSum
REAL                           :: L_eta_zeta
REAL                           :: Var(0:PP_N,0:PP_N,0:PP_N)
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
FUNCTION NewtonLogLaw(velx,nu,abs_h_wm)
! MODULES
USE MOD_WMLES_Vars
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL, INTENT(IN)            :: velx         ! Tangential velocity to feed wall model (log-law function)
REAL, INTENT(IN)            :: nu           ! kinematic viscosity at h_wm
REAL, INTENT(IN)            :: abs_h_wm     ! Wall-model height (h_wm)
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
        DO WHILE(iter.LT.10) ! Maximum number of iterations. Usually, u_tau is found with 1~4 iters
            iter = iter+1
            f = NewtonLogLaw*( (1./vKarman)*LOG(1. + vKarman*abs_h_wm*NewtonLogLaw/nu) + (B - (1./vKarman)*LOG(vKarman))*(1. - EXP(-abs_h_wm*NewtonLogLaw/(11.*nu)) - (abs_h_wm*NewtonLogLaw/(11.*nu))*EXP(-abs_h_wm*NewtonLogLaw/(3.*nu)) ) ) - velx
            fprime = (1./vKarman)*LOG(1. + vKarman*abs_h_wm*NewtonLogLaw/nu) + (B - (1./vKarman)*LOG(vKarman))*(1. - EXP(-abs_h_wm*NewtonLogLaw/(11.*nu)) - (abs_h_wm*NewtonLogLaw/(11.*nu))*EXP(-abs_h_wm*NewtonLogLaw/(3.*nu)) ) &
                + abs_h_wm*NewtonLogLaw/(nu + vKarman*abs_h_wm*NewtonLogLaw) &
                + NewtonLogLaw*(B - (1./vKarman)*LOG(vKarman))*((abs_h_wm/(11.*nu))*(EXP(-abs_h_wm*NewtonLogLaw/(11.*nu)) - EXP(-abs_h_wm*NewtonLogLaw/(3.*nu)) + (abs_h_wm*NewtonLogLaw/(11.*nu))*EXP(-abs_h_wm*NewtonLogLaw/(3.*nu)) ) )
            ! IF(ABS(fprime).LE.1.E-5) EXIT ! fprime ~ 0 -- INCOMING OVERFLOW BY DIVISION
            NewtonLogLaw = NewtonLogLaw - f/fprime            
            IF (ABS(f/fprime).LE.1.E-5) EXIT ! 1.E-5 = stop criterion (tolerance)
        END DO

END SELECT


END FUNCTION NewtonLogLaw


!==================================================================================================================================
!> Evaluate u_tau from a log-law given instantaneous velocity, using Newton method
!==================================================================================================================================
FUNCTION NewtonCouette(velx,nu,abs_h_wm)
! MODULES
USE MOD_WMLES_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL, INTENT(IN)            :: velx         ! Tangential velocity to feed wall model (log-law function)
REAL, INTENT(IN)            :: nu           ! kinematic viscosity at h_wm
REAL, INTENT(IN)            :: abs_h_wm     ! Wall-model height (h_wm)
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
SUBROUTINE StartReceiveHwmMPIData()
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
INTEGER                     :: startIdx, endIdx, bufferSize
!==================================================================================================================================
DO iProc=1,nWMLESRecvProcs
    startIdx = WMLESRecvRange(1, iProc)
    endIdx = WMLESRecvRange(2, iProc)
    bufferSize = (nHWMPropSend+4)*nHWMRecvPoints(iProc)
    CALL MPI_IRECV(HWMRecvInfo(:, startIdx:endIdx),bufferSize, MPI_DOUBLE_PRECISION, &
                    WMLESRecvProc(iProc), 0, MPI_COMM_FLEXI, WMLES_RecvRequests(iProc),iError)
END DO

END SUBROUTINE StartReceiveHwmMPIData


!==================================================================================================================================
!> Subroutine that controls the send operations for the Tau_W to be exchanged between processors.
!==================================================================================================================================
SUBROUTINE StartSendHwmMPIData()
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
INTEGER                     :: startIdx, endIdx, bufferSize
!==================================================================================================================================
DO iProc=1,nWMLESSendProcs
    startIdx = WMLESSendRange(1, iProc)
    endIdx = WMLESSendRange(2, iProc)
    bufferSize = (nHWMPropSend+4)*nHWMSendPoints(iProc)
    CALL MPI_ISEND(HWMSendInfo(:,startIdx:endIdx),bufferSize, MPI_DOUBLE_PRECISION, &
                    WMLESSendProc(iProc), 0, MPI_COMM_FLEXI, WMLES_SendRequests(iProc),iError)
END DO

END SUBROUTINE StartSendHwmMPIData

!==================================================================================================================================
!> Subroutine that completes the receive operations for the Tau_W to be exchanged between processors.
!==================================================================================================================================
SUBROUTINE FinishExchangeHWMData()
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
INTEGER                             :: i, k
!==================================================================================================================================
CALL MPI_Waitall(nWMLESRecvProcs,WMLES_RecvRequests,MPI_STATUSES_IGNORE,iError)

! In addition to guaranteeing that the info has been received and HWMRecvInfo is populated,
! we also manage the entries in HWMInfo, using the mapping received in the message itself.

DO i=1,WMLESRecvRange(2,nWMLESRecvProcs)
    HWMInfo(1, INT(HWMRecvInfo(nHWMPropSend+2, i)), & ! p
                INT(HWMRecvInfo(nHWMPropSend+3, i)), & ! q
                INT(HWMRecvInfo(nHWMPropSend+4, i))-offsetBCSides) & ! iWMLESSide
            = HWMRecvInfo(nHWMPropSend+1, i) ! h_wm value

    DO k=1,nHWMPropSend
        HWMInfo(k+1, INT(HWMRecvInfo(nHWMPropSend+2, i)), & ! p
                INT(HWMRecvInfo(nHWMPropSend+3, i)), & ! q
                INT(HWMRecvInfo(nHWMPropSend+4, i))-offsetBCSides) & ! iWMLESSide
            = HWMRecvInfo(k, i) ! flow property corresponding to index k
    END DO
END DO

END SUBROUTINE FinishExchangeHWMData


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
SDEALLOCATE(nHWMSendPoints)
SDEALLOCATE(SendToInterpPoint)
SDEALLOCATE(WMLESSendProc)
SDEALLOCATE(HWMSendInfo)
SDEALLOCATE(WMLESSendRange)
SDEALLOCATE(WMLESRecvProc)
SDEALLOCATE(nHWMRecvPoints)
SDEALLOCATE(HWMRecvInfo)
SDEALLOCATE(WMLESRecvRange)
SDEALLOCATE(HWMInterpInfo)
SDEALLOCATE(LocalToInterpPoint)
SDEALLOCATE(HWMLocalInfo)
SDEALLOCATE(HWMInfo)
SDEALLOCATE(WMLES_RecvRequests)
SDEALLOCATE(WMLES_SendRequests)


END SUBROUTINE FinalizeWMLES



!==================================================================================================================================
END MODULE MOD_WMLES
