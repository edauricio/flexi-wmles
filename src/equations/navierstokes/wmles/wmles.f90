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
PUBLIC :: ComputeWallStress, EvalDiffFlux3D_WMLES, FalknerSkan
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
CALL addStrListEntry('WallModel', 'FalknerSkan', WMLES_FALKNER_SKAN)
CALL prms%CreateRealOption('Reynolds', "Reynolds number of the simulation -- important for Laminar wall-modeling")
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
REAL                            :: WMLES_Tol, RealTmp, SSig, Re
INTEGER                         :: SideID, iSide, GlobalWMLESSide, LocSide, MineCnt, nSideIDs, offsetSideID, iGlobalBCSide, NMax
INTEGER                         :: Loc_hwmElemID, Glob_hwmElemID, OppSideID
REAL                            :: hwmElem_NodeCoords(3,0:NGeo,0:NGeo,0:NGeo), OppSideNodeCoords(3,4)
INTEGER                         :: TotalNWMLESSides, GlobalOppSideID, InnerOppSideID
REAL                            :: OppSideEdge(3,2), OppSideNormal(3), hwmDistVec(3), TolVec(3)
INTEGER                         :: i,j,k,p,q,r

INTEGER, ALLOCATABLE            :: WMLESToBCSide_tmp(:), TauW_Proc_tmp(:,:,:,:)
REAL, ALLOCATABLE               :: OthersPointInfo(:,:,:) ! indices-- 1-3: hwm_Coords, 4-6: TangVec1, 7: Glob_hwmElemID
INTEGER                         :: WallStressCount_local(0:nProcessors-1), WallStressCount(0:nProcessors-1)
INTEGER                         :: FirstElemInd, LastElemInd, FirstSideInd, LastSideInd
REAL, ALLOCATABLE               :: xi(:), fps(:), FSBeta_tmp(:,:,:), FSAlpha_tmp(:,:,:)
REAL                            :: h_wm_Coords(3), DistanceVect(3), Distance, inner_prod, beta_l, etainf, ddfddn
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
INTEGER, ALLOCATABLE            :: BCSideToWMLES_global(:), WMLESShapeInfo_global(:,:), WMLESToTurbulentSide_tmp(:), WMLESToLaminarSide_tmp(:)
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
! (0) Bear in mind that all info in the WM connection file is in terms of GLOBAL elements and BC side numbers
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

! Read boundary shape information
ALLOCATE(WMLESShapeInfo_global(2,nWMLESSides_global))
CALL ReadArray('sideShapeInfo',2,(/2,nWMLESSides_global/),0,1,IntArray=WMLESShapeInfo_global)

! Read connection information
ALLOCATE(WMConnection(N_INTERFACE_PARAMS,0:PP_N,0:PP_N,nWMLESSides_global))
CALL ReadArray('WM',4,&
               (/N_INTERFACE_PARAMS,PP_N+1,PP_N+1,nWMLESSides_global/),&
               0,4,RealArray=WMConnection)

CALL CloseDataFile()

! Consistency check between given WM Connection file and chosen WMLES Model
! (for complete models, file must have 2 shapes; for single models, file must have only 1 shape)
SELECT CASE(WallModel)
    CASE(WMLES_FSLL)
        IF (MAXVAL(WMLESShapeInfo_global(1,:)).NE.2) &
            CALL Abort(__STAMP__, 'Complete modeling strategy (laminar+turbulent regions) was selected, but a WM Connection File with a single WMLES shape was given.')

    CASE DEFAULT
        IF (MAXVAL(WMLESShapeInfo_global(1,:)).EQ.2) &
            CALL Abort(__STAMP__, 'WM Connection File with two shapes for a complete model (laminar+turbulent) was given, but a single modeling strategy was selected.')
END SELECT

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
ALLOCATE(WMLESShapeInfo(2,nWMLESSides))

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
WMLESShapeInfo = 0


DO iGlobalBCSide=1,nBCSides_global
    GlobalWMLESSide = BCSideToWMLES_global(iGlobalBCSide)
    IF (GlobalWMLESSide.EQ.0) CYCLE ! Not a modelled side
    IF (GlobalWMLESSide.GT.nWMLESSides_global) CALL Abort(__STAMP__, "Something's wrong with the global WMLES Side number information.")

    IF ((iGlobalBCSide.GT.offsetBCSides) .AND. (iGlobalBCSide.LE.(offsetBCSides+nBCSides))) THEN ! Local BC side!
        ! Sanity check
        IF (BCSideToWMLES(iGlobalBCSide-offsetBCSides) .EQ. 0) CALL Abort(__STAMP__, &
            'Modeled local BC found in WM Connection file, but it has not been marked as a WMLES Side using the current config/mesh!')

        iWMLESSide = BCSideToWMLES(iGlobalBCSide-offsetBCSides)
        WMLESShapeInfo(:,iWMLESSide) = WMLESShapeInfo_global(:,GlobalWMLESSide)

        ! We now check if we need to receive h_wm point info from another proc., or if h_wm is within a local element
        ! This check is done for each GL point within this BC side
        DO p=0,PP_N; DO q=0,PP_N

            iElem = INT(WMConnection(INTERFACE_ELEMSEND,p,q,GlobalWMLESSide)) ! Global ElemID
            IF ((iElem.GT.offsetElem) .AND. (iElem.LE.(offsetElem+nElems))) THEN
                ! Element where h_wm lies is local to this proc
                ! Hence, no communication will be needed. 
                ! We count this h_wm point as "interpolation needed", and store standard coordinate 
                ! information for later interpolation
                ! We also store the BC global side and the p,q, point that this interpolation corresponds to

                nHWMLocalPoints = nHWMLocalPoints + 1
                nHWMInterpPoints = nHWMInterpPoints + 1

                LocalToInterpPoint_tmp(nHWMLocalPoints) = nHWMInterpPoints

                HWMLocalInfo_tmp(1:4, nHWMLocalPoints) = (/WMConnection(INTERFACE_L,p,q,GlobalWMLESSide), REAL(p), REAL(q), REAL(iWMLESSide)/)
                HWMInterpInfo_tmp(:,nHWMInterpPoints) = (/REAL(iElem-offsetElem), & ! All elements in the array construction must be of the same type.
                                                        WMConnection(INTERFACE_XI, p, q, GlobalWMLESSide), &
                                                        WMConnection(INTERFACE_ETA,p, q, GlobalWMLESSide), &
                                                        WMConnection(INTERFACE_ZETA,p,q, GlobalWMLESSide)/)

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

                HWMRecvInfo_tmp(1:4, nHWMRecvPoints_tmp(iProc), iProc) = (/WMConnection(INTERFACE_L, p, q, GlobalWMLESSide), REAL(p), REAL(q), REAL(iWMLESSide)/)

            END IF ! IF ELSE element containing h_wm is local to this proc.
        END DO; END DO ! p, q

    ELSE 
        ! This global WMLES Side is not within this proc.
        ! However, this proc. may need to SEND h_wm info to the rank responsible for the BC imposition on this side

        DO p=0,PP_N; DO q=0,PP_N

            iElem = INT(WMConnection(INTERFACE_ELEMSEND,p,q,GlobalWMLESSide)) ! Global ElemID
            IF ((iElem.GT.offsetElem) .AND. (iElem.LE.(offsetElem+nElems))) THEN
                ! This proc. is not responsible for WMLES BC imposition, but it contains
                ! the info in h_wm point.
                ! Hence, we need to SEND this info to the proc. responsible for BC imposition

                iProc = ELEMIPROC(INT(WMConnection(INTERFACE_ELEMRECV,p,q,GlobalWMLESSide)))

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
                
                HWMSendInfo_tmp(nHWMPropSend+1:nHWMPropSend+4,nHWMSendPoints_tmp(iProc),iProc) = (/WMConnection(INTERFACE_L,p,q,GlobalWMLESSide), REAL(p), REAL(q), REAL(iGlobalBCSide)/)

                nHWMInterpPoints = nHWMInterpPoints + 1
                HWMInterpInfo_tmp(:,nHWMInterpPoints) = (/REAL(iElem-offsetElem), &
                                                        WMConnection(INTERFACE_XI, p, q, GlobalWMLESSide), &
                                                        WMConnection(INTERFACE_ETA,p, q, GlobalWMLESSide), &
                                                        WMConnection(INTERFACE_ZETA,p,q, GlobalWMLESSide)/)

                SendToInterpPoint_tmp(nHWMSendPoints_tmp(iProc), iProc) = nHWMInterpPoints

            END IF ! IF element containing h_wm is local to this proc.

        END DO; END DO ! p, q
    END IF ! IF ELSE Local BC Side
END DO

SDEALLOCATE(BCSideToWMLES_global)
SDEALLOCATE(WMLESShapeInfo_global)
SDEALLOCATE(WMConnection)

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
SDEALLOCATE(HWMSendInfo_tmp)
SDEALLOCATE(nHWMSendPoints_tmp)
SDEALLOCATE(SendToInterpPoint_tmp)
SDEALLOCATE(WMLESSendToProc)

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
SDEALLOCATE(HWMRecvInfo_tmp)
SDEALLOCATE(nHWMRecvPoints_tmp)
SDEALLOCATE(WMLESRecvFromProc)

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
SDEALLOCATE(HWMInterpInfo_tmp)
SDEALLOCATE(HWMLocalInfo_tmp)
SDEALLOCATE(LocalToInterpPoint_tmp)

! We now calculate the local, point-wise beta's defining the Falkner-Skan equation
! profile to be solved for each point in the boundary layer pertaining to the laminar
! region being modeled by the FS model, if any
! In addition, since the final solution of the FS equation only depends on this beta,
! we cache the solution points (along the similarity variable -- eta --  direction)
! for later use in the model
SELECT CASE(WallModel)
    CASE(WMLES_FSLL,WMLES_FALKNER_SKAN)
        Re = GETREAL("Reynolds")
        ALLOCATE(FSBeta_tmp(0:PP_N,0:PP_N,nWMLESSides))!,FSDelta_tmp(0:PP_N,0:PP_N,nWMLESSides))
        ALLOCATE(FSAlpha_tmp(0:PP_N,0:PP_N,nWMLESSides))!,FSEtaInf_tmp(0:PP_N,0:PP_N,nWMLESSides))
        ALLOCATE(WMLESToLaminarSide(nWMLESSides),WMLESToTurbulentSide(nWMLESSides))
        NMax = 0
        nWMLaminarSides = 0
        nWMTurbulentSides = 0
        WMLESToLaminarSide = 0
        WMLESToTurbulentSide = 0
        ! Allocate generous temporary buffer for FS solution at each laminar boundary point
        !ALLOCATE(FSEta_tmp(1000,0:PP_N,0:PP_N,nWMLESSides),FSPrime_tmp(1000,0:PP_N,0:PP_N,nWMLESSides))
        !FSEta_tmp = 0
        !FSPrime_tmp = 0
       ! OPEN(UNIT=211,  &
       ! FILE='FSDebug_Pre_Beta9.csv',      &
       ! STATUS='UNKNOWN',  &
       ! ACTION='WRITE', &
       ! POSITION='APPEND')
       !  WRITE(211,*) "wmside,bcside,surf,x,p,q,beta,delta_eta,alpha,etainf"
        DO iSide=1,nWMLESSides
            IF (WMLESShapeInfo(1,iSide).EQ.1) THEN ! Laminar region
                nWMLaminarSides = nWMLaminarSides + 1
                WMLESToLaminarSide(iSide) = nWMLaminarSides
                DO p=0,PP_N; DO q=0,PP_N
                    inner_prod = DOT_PRODUCT(NormVec(1:3,p,q,0,WMLESToBCSide(iSide)), (/0., 1., 0./))
                    ! TODO: Adjust for angle of attack (freestream flow angle)
                    IF (inner_prod.LT.0) THEN ! upper surface
                        beta_l = SIGN(1.0,-TangVec2(2,p,q,0,WMLESToBCSide(iSide)))*ACOS(DOT_PRODUCT(-TangVec2(1:3,p,q,0,WMLESToBCSide(iSide)), (/1., 0. ,0./)))
                    ELSE
                        beta_l = SIGN(1.0,-TangVec2(2,p,q,0,WMLESToBCSide(iSide)))*ACOS(DOT_PRODUCT(TangVec2(1:3,p,q,0,WMLESToBCSide(iSide)), (/1., 0. ,0./)))
                    END IF
                    ! Three betas:
                    ! 1) Only geometric
                    ! 2) Geometric + taking into account x-position
                    ! 3) Geometric + taking x-position + taking beta itself into account
                    ! 4) Geometric + taking into account x-position (differently than 2)
                    ! The logic is that, since there is a boundary layer growing in the x direction,
                    ! we must have a "smoothing" of the wedge beta, because the local flow will feel
                    ! it differently in terms of acceleration and wall shear stress due to this boundary layer
                    ! And also, if we have a positive beta, the boundary layer is growing slower, whereas
                    ! for negative betas, it is the contrary.
                    
                    !1) 
                    !FSBeta_tmp(p,q,nWMLaminarSides) = beta_l*(2./PI)
                    !2)
                    !FSBeta_tmp(p,q,nWMLaminarSides) = beta_l*(2./PI)*EXP(-Face_xGP(1,p,q,0,WMLESToBCSide(iSide)))
                    !3)
                    !FSBeta_tmp(p,q,nWMLaminarSides) = 0.39*beta_l*(2./PI)*EXP(0.3*(beta_l*(2./PI)) - Face_xGP(1,p,q,0,WMLESToBCSide(iSide)))
                    !4)
                    !FSBeta_tmp(p,q,nWMLaminarSides) = beta_l*(2./PI)*EXP(Face_xGP(1,p,q,0,WMLESToBCSide(iSide))-1.)
                    !5)
                    !FSBeta_tmp(p,q,nWMLaminarSides) = beta_l*(2./PI)*EXP(Face_xGP(1,p,q,0,WMLESToBCSide(iSide)) - 1. - 0.3*(beta_l*(2./PI)))
                    !6) Based on "stretched" Sigmoid function -- https://math.stackexchange.com/questions/1214167/how-to-remodel-sigmoid-function-so-as-to-move-stretch-enlarge-it
                    ! We use the function below to decrease Beta based on the Reynolds number of the flow: 
                    ! Z(t) = 1 / (1 + EXP(-(1/Re_c) * LOG(b) * Re)), where Re_c is the Reynolds where we want the inflection point in the function,
                    ! and b is a constant which also determines such inflection point. We are currently using Sigmoid b = 2 + SQRT(3)
                    ! SSig = 1. / (1. + EXP(-(1./50000.) * LOG(2.+SQRT(3.)) *  Re))
                    ! FSBeta_tmp(p,q,nWMLaminarSides) = beta_l*(2./PI) - ABS(beta_l*(2./PI)) * (1. - SSig)
                    !7) This is (6) slightly modified to take into account x-distance (growth of boundary layer)
                    ! SSig = 1. / (1. + EXP(-(1./80000.) * LOG(2. + SQRT(3.) - (5./2.)*Face_xGP(1,p,q,0,WMLESToBCSide(iSide))) *  Re))
                    ! FSBeta_tmp(p,q,nWMLaminarSides) = beta_l*(2./PI) - ABS(beta_l*(2./PI)) * (1. - SSig)
                    !8) This is (7) slightly modified
                    ! SSig = (1. - Face_xGP(1,p,q,0,WMLESToBCSide(iSide))*EXP(-(1./50000.) * LOG(2. + SQRT(3.)) *  Re)) / (1. + EXP(-(1./50000.) * LOG(2. + SQRT(3.)) *  Re))
                    ! FSBeta_tmp(p,q,nWMLaminarSides) = beta_l*(2./PI) - ABS(beta_l*(2./PI)) * (1. - SSig)
                    !9) This is (8) slightly modified -- numerator now is in terms of SQRT(x) instead of linear in x
                    SSig = (1. - SQRT(Face_xGP(1,p,q,0,WMLESToBCSide(iSide)))*EXP(-(1./50000.) * LOG(2. + SQRT(3.)) *  Re)) / (1. + EXP(-(1./50000.) * LOG(2. + SQRT(3.)) *  Re))
                    FSBeta_tmp(p,q,nWMLaminarSides) = beta_l*(2./PI) - ABS(beta_l*(2./PI)) * (1. - SSig)
                    IF (FSBeta_tmp(p,q,nWMLaminarSides) .LE. -0.195) FSBeta_tmp(p,q,nWMLaminarSides) = -0.19
                    
                    ! Solve FS once for each point and cache the solution (needed later)                    
                    CALL FalknerSkan(1.0, FSBeta_tmp(p,q,nWMLaminarSides), etainf, ddfddn, xi, fps)
                    FSAlpha_tmp(p,q,nWMLaminarSides) = ddfddn
                    ! FSEtaInf_tmp(p,q,nWMLaminarSides) = etainf
                    IF (SIZE(xi,1).GT.NMax) NMax = SIZE(xi,1)
                    IF (NMax.GT.1000) CALL Abort(__STAMP__,'Unexpected size for "XI" vector of FalknerSkan solution')
                    !FSEta_tmp(1:SIZE(xi,1),p,q,nWMLaminarSides) = xi(:)
                    !FSPrime_tmp(1:SIZE(fps,1),p,q,nWMLaminarSides) = fps(:)
                    

                    ! Look for eta_delta (i.e., eta such that fprime ~ 0.99)
                    ! FSDelta_tmp(p,q,nWMLaminarSides) = etainf
                    ! DO i=1,SIZE(fps,1)
                    !     IF (fps(i)>=0.99) THEN
                    !         FSDelta_tmp(p,q,nWMLaminarSides) = xi(i)*etainf
                    !         EXIT
                    !     END IF
                    ! END DO
                    ! IF (BC(WMLESToBCSide(iSide)).EQ.1) THEN
                    !     WRITE(211,*) iSide&
                    !                 ,",",WMLESToBCSide(iSide)&
                    !                 ,",upper,"&
                    !                 ,Face_xGP(1,p,q,0,WMLESToBCSide(iSide)) &
                    !                 ,",",p &
                    !                 ,",",q &
                    !                 ,",",FSBeta_tmp(p,q,iSide) &
                    !                 ,",",FSDelta_tmp(p,q,iSide) &
                    !                 ,",",FSAlpha_tmp(p,q,iSide) &
                    !                 ,",",FSEtaInf_tmp(p,q,iSide)
                    ! ELSE IF (BC(WMLESToBCSide(iSide)).EQ.2) THEN
                    !     WRITE(211,*) iSide&
                    !                 ,",",WMLESToBCSide(iSide)&
                    !                 ,",lower,"&
                    !                 ,Face_xGP(1,p,q,0,WMLESToBCSide(iSide)) &
                    !                 ,",",p &
                    !                 ,",",q &
                    !                 ,",",FSBeta_tmp(p,q,iSide) &
                    !                 ,",",FSDelta_tmp(p,q,iSide) &
                    !                 ,",",FSAlpha_tmp(p,q,iSide) &
                    !                 ,",",FSEtaInf_tmp(p,q,iSide)
                    ! ELSE 
                    !     WRITE(211,*) "none"
                    ! END IF
                END DO; END DO
            ELSE
                nWMTurbulentSides = nWMTurbulentSides + 1
                WMLESToTurbulentSide(iSide) = nWMTurbulentSides
            END IF
        END DO

        ! Allocate permanent memory and copy data from buffer; free buffer next
        ALLOCATE(FSBeta(0:PP_N,0:PP_N,nWMLaminarSides))!,FSDelta(0:PP_N,0:PP_N,nWMLaminarSides))
        ALLOCATE(FSWallShear(0:PP_N,0:PP_N,nWMLaminarSides))!,FSEtaInf(0:PP_N,0:PP_N,nWMLaminarSides),)
        !ALLOCATE(FSEta(NMax,0:PP_N,0:PP_N,nWMLaminarSides),FSFPrime(NMax,0:PP_N,0:PP_N,nWMLaminarSides))
        

        DO iSide=1,nWMLaminarSides
            DO p=0,PP_N; DO q=0,PP_N
                FSBeta(p,q,iSide) = FSBeta_tmp(p,q,iSide)
                ! FSDelta(p,q,iSide) = FSDelta_tmp(p,q,iSide)
                FSWallShear(p,q,iSide) = FSAlpha_tmp(p,q,iSide)
                ! FSEtaInf(p,q,iSide) = FSEtaInf_tmp(p,q,iSide)
                ! DO i=1,NMax
                !     FSEta(i,p,q,iSide) = FSEta_tmp(i,p,q,iSide)
                !     FSFPrime(i,p,q,iSide) = FSPrime_tmp(i,p,q,iSide)
                ! END DO
            END DO; END DO
        END DO
        DEALLOCATE(FSBeta_tmp)
        !DEALLOCATE(FSDelta_tmp)
        !DEALLOCATE(FSEta_tmp)
        !DEALLOCATE(FSPrime_tmp)
        DEALLOCATE(FSAlpha_tmp)
        !DEALLOCATE(FSEtaInf_tmp)

END SELECT
SDEALLOCATE(WMLESShapeInfo)

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
USE MOD_WMLES_Utils                 ,ONLY: GetParams
USE MOD_Mesh_Vars
USE MOD_Interpolation_Vars
USE MOD_TimeDisc_Vars               ,ONLY: t
USE MOD_Lifting_Vars                ,ONLY: gradUy_master
USE MOD_DG_Vars                     
USE MOD_EOS_Vars                    ,ONLY: mu0
IMPLICIT NONE

!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                             :: sProc, FPInd, IPInd, IntPInd, adj_uinv
INTEGER                             :: p,q,r,ElemID
REAL                                :: u_tau, u_mean, tau_w_mag, utang, VelMag, fp, eta_wm, inner_prod, eps
REAL                                :: vel_inst(3), rho_inst, tangvec(3), tau_w_vec(3), eta_root
INTEGER                             :: SurfLowUp ! Debug

INTEGER                             :: i, k, SideID, locSide
INTEGER                             :: iSendPoint, iLocalPoint
REAL                                :: xi, eta, zeta
REAL                                :: L_xi(0:PP_N), L_eta(0:PP_N), L_zeta(0:PP_N)
LOGICAL                             :: isLE, fopen
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
IF (nWMLESSendProcs.NE.0) THEN
    ! IF above avoid errors when indexing WMLESSendRange -- even though when compiled in RELEASE mode everything works fine....
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
END IF

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
END DO

! Finish receiving h_wm info from other procs.
CALL FinishExchangeHWMData()

! Finally, with all h_wm info properly arranged into HWMInfo array, we can
! calculate the wall stress according to the chosen model

SELECT CASE(WallModel)

  CASE (WMLES_LOGLAW,WMLES_REICHARDT)

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

            utang = DOT_PRODUCT(HWMInfo(2:4,p,q,SideID),tangvec)

            u_tau = NewtonLogLaw(utang, (mu0/UPrim_master(1,p,q,WMLESToBCSide(SideID))), HWMInfo(1,p,q,SideID)) ! rho_wall is used here (through UPrim_master). Read TODO above
            tau_w_mag = UPrim_master(1,p,q,WMLESToBCSide(SideID))*(u_tau**2) ! CHECK TODO ABOVE
            tau_w_vec = tau_w_mag*tangvec
            ! We now project tau_w_vec onto local coordinates, since WMLES_TauW is used in a local context
            SELECT CASE(SideToElem(S2E_LOC_SIDE_ID, WMLESToBCSide(SideID)))
                CASE (ZETA_PLUS) ! Bottom wall
                    WMLES_TauW(1,p,q,SideID) = DOT_PRODUCT(tau_w_vec(1:3),TangVec1(1:3,p,q,0,WMLESToBCSide(SideID)))
                    WMLES_TauW(2,p,q,SideID) = DOT_PRODUCT(tau_w_vec(1:3),TangVec2(1:3,p,q,0,WMLESToBCSide(SideID)))
                CASE (XI_PLUS,ETA_MINUS) ! Bottom wall
                    WMLES_TauW(1,p,q,SideID) = -1.*DOT_PRODUCT(tau_w_vec(1:3),TangVec2(1:3,p,q,0,WMLESToBCSide(SideID)))
                    WMLES_TauW(2,p,q,SideID) = DOT_PRODUCT(tau_w_vec(1:3),TangVec1(1:3,p,q,0,WMLESToBCSide(SideID)))
                CASE (ZETA_MINUS) ! Top wall
                    WMLES_TauW(1,p,q,SideID) = -1.*DOT_PRODUCT(tau_w_vec(1:3),TangVec1(1:3,p,q,0,WMLESToBCSide(SideID)))
                    WMLES_TauW(2,p,q,SideID) = DOT_PRODUCT(tau_w_vec(1:3),TangVec2(1:3,p,q,0,WMLESToBCSide(SideID)))
                CASE (XI_MINUS,ETA_PLUS) ! Top wall
                    WMLES_TauW(1,p,q,SideID) = -1.*DOT_PRODUCT(tau_w_vec(1:3),TangVec2(1:3,p,q,0,WMLESToBCSide(SideID)))
                    WMLES_TauW(2,p,q,SideID) = -1.*DOT_PRODUCT(tau_w_vec(1:3),TangVec1(1:3,p,q,0,WMLESToBCSide(SideID)))
            END SELECT

            ! LOGWRITE(*,'(2(I9,1X), I7,1X, 2(I2,1X), 3(E10.4,1X), 3(E10.4,1X), 2(E10.4,1X))') &
            !         WMLESToBCSide(SideID)+offsetBCSides, SideID, WMLESToBCSide(SideID), p, q, HWMInfo(2,p,q,SideID), HWMInfo(3,p,q,SideID), HWMInfo(4,p,q,SideID), &
            !         utang, HWMInfo(1,p,q,SideID), u_tau, WMLES_TauW(1,p,q,SideID), WMLES_TauW(2,p,q,SideID)
            ! LOGWRITE(*,*) 'loc, NormVec', SideToElem(S2E_LOC_SIDE_ID,WMLESToBCSide(SideID)), NormVec(:,p, q, FV_ENABLED, SideID)
            ! LOGWRITE(*,*) 'loc, TangVec1', SideToElem(S2E_LOC_SIDE_ID,WMLESToBCSide(SideID)), TangVec1(:,p, q, FV_ENABLED, WMLESToBCSide(SideID))
            ! LOGWRITE(*,*) 'loc, TangVec2', SideToElem(S2E_LOC_SIDE_ID,WMLESToBCSide(SideID)), TangVec2(:,p, q, FV_ENABLED, WMLESToBCSide(SideID))
            ! LOGWRITE(*,*) 'loc, Face_xGP', SideToElem(S2E_LOC_SIDE_ID,WMLESToBCSide(SideID)), Face_xGP(:,p, q, FV_ENABLED, WMLESToBCSide(SideID))
        END DO; END DO ! p,q
    END DO ! iSide

  CASE (WMLES_FALKNER_SKAN)
    ! IF (ALMOSTEQUALABSOLUTE(t,1E-2,1E-6)) THEN
    ! INQUIRE(unit=211, opened=fopen)
    ! OPEN(UNIT=211,  &
    !    FILE='FSDebug.csv',      &
    !    STATUS='UNKNOWN',  &
    !    ACTION='WRITE', &
    !    POSITION='APPEND')
    !  IF(.NOT.fopen)   WRITE(211,*) "wmside,bcside,surf,x,p,q,beta,delta_eta,alpha,t,eta_wm,eta_root,mu0,utang,adj_uinv,velx,vely,dudy,hwvx,hwvy,tangx,tangy,innerprod,wmxy,wmxz"
    ! END IF
    DO SideID=1,nWMLESSides
        DO p=0,PP_N; DO q=0,PP_N
            isLE = .FALSE.
            IF (ABS(Face_xGP(1,p,q,0,WMLESToBCSide(SideID))).LE.1E-6) isLE = .TRUE.

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

            utang = DOT_PRODUCT(HWMInfo(2:4,p,q,SideID),tangvec)
            ! Should we just get the x-component of the velocity in this case?
            ! Since tau_xy is not being computed directly from this projection... 
            ! I mean, the computed tau_xy is indeed aligned with the x-direction according to the model...
            ! adj_uinv = 0
            !! IF (.NOT.isLE) THEN ! If leading edge, then about any eta_wm will suffice, and utang is taken as above
                ! eps = 0.
                ! IF (isLE) eps = 1E-6
                ! DO WHILE(adj_uinv.LE.100)
                !     eta_root = SQRT((utang**3) / ((2.-FSBeta(p,q,WMLESToLaminarSide(SideID))) * (mu0/UPrim_master(1,p,q,WMLESToBCSide(SideID))) * (Face_xGP(1,p,q,0,WMLESToBCSide(SideID)) + eps)))
                !     !eta_wm = HWMInfo(1,p,q,SideID)*eta_root
                !     eta_wm = HWMInfo(1,p,q,SideID)*eta_root/utang
                !     IF (eta_wm.GE.FSDelta(p,q,WMLESToLaminarSide(SideID)) .OR. ALMOSTEQUALABSOLUTE(eta_wm, FSDelta(p,q,WMLESToLaminarSide(SideID)), 1E-3)) EXIT
                !     fp = GetNewFP(eta_wm/FSEtaInf(p,q,WMLESToLaminarSide(SideID)), FSEta(:,p,q,WMLESToLaminarSide(SideID)), FSFPrime(:,p,q,WMLESToLaminarSide(SideID)))
                !     utang = utang/fp
                !     adj_uinv = adj_uinv + 1
                ! END DO
            !! END IF
            ! STRATEGY #2: No adjustment; only guarantee that we do not have a division-by-zero at the leading edge
            eps=0.
            IF (isLE) eps = 1E-6
            eta_root = SQRT((utang**3) / ((2.-FSBeta(p,q,WMLESToLaminarSide(SideID))) * (mu0/UPrim_master(1,p,q,WMLESToBCSide(SideID))) * (Face_xGP(1,p,q,0,WMLESToBCSide(SideID)) + eps)))
            IF (isLE) eta_root = 0.
            tau_w_vec = mu0*FSWallShear(p,q,WMLESToLaminarSide(SideID))*eta_root*tangvec

            ! We now project tau_w_vec onto local coordinates, since WMLES_TauW is used in a local context
            !IF (isLE) THEN ! leading edge (x ~ 0)
                ! Check tangents for upper/lower surface
            !ELSE
                ! sign of inner_prod takes into account upper/lower surfaces
                inner_prod = DOT_PRODUCT(NormVec(1:3,p,q,0,WMLESToBCSide(SideID)), (/0., 1., 0./))
                WMLES_TauW(1,p,q,SideID) = -1.*DOT_PRODUCT(tau_w_vec(1:3),TangVec2(1:3,p,q,0,WMLESToBCSide(SideID)))
                WMLES_TauW(2,p,q,SideID) = SIGN(1.,-inner_prod)*DOT_PRODUCT(tau_w_vec(1:3),TangVec1(1:3,p,q,0,WMLESToBCSide(SideID)))
                ! TEST BELOW!! (See last commit)
                ! WMLES_TauW(1,p,q,SideID) = DOT_PRODUCT(tau_w_vec(1:3),TangVec2(1:3,p,q,0,WMLESToBCSide(SideID)))
                ! WMLES_TauW(2,p,q,SideID) = 0. ! TEST!!!
            ! END IF
            ! IF (ALMOSTEQUALABSOLUTE(t,1E-2,1E-6)) THEN
            ! WRITE(211,*) SideID &
            !         ,",",WMLESToBCSide(SideID) &
            !         ,",",BC(WMLESToBCSide(SideID))  &
            !         ,",",Face_xGP(1,p,q,0,WMLESToBCSide(SIdeID))  &
            !         ,",",p  &
            !         ,",",q  &
            !         ,",",FSBeta(p,q,WMLESToLaminarSide(SideID))  &
            !         ,",",FSDelta(p,q,WMLESToLaminarSide(SideID))  &
            !         ,",",FSWallShear(p,q,WMLESToLaminarSide(SideID))  &
            !         ,",",t  &
            !         ,",",0.  & !eta_wm
            !         ,",",eta_root  &
            !         ,",",mu0  &
            !         ,",",utang  &
            !         ,",",0. & !adj_uinv
            !         ,",",UPrim_master(2,p,q,WMLESToBCSide(SideID)) &
            !         ,",",UPrim_master(3,p,q,WMLESToBCSide(SideID)) &
            !         ,",",gradUy_master(2,p,q,WMLESToBCSide(SideID)) &
            !         ,",",HWMInfo(2,p,q,SideID) &
            !         ,",",HWMInfo(3,p,q,SideID) &
            !         ,",",tangvec(1) &
            !         ,",",tangvec(2) &
            !         ,",",inner_prod &
            !         ,",",WMLES_TauW(1,p,q,SideID)&
            !         ,",",WMLES_TauW(2,p,q,SideID)
            ! END IF
        END DO; END DO
    END DO


  CASE DEFAULT
    CALL Abort(__STAMP__,&
         'Wall Model not yet implemented!')

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
!> Evaluate solution at the point corresponding to the wall model height h_wm within element iElem
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
!> Interpolate Falkner-Skan solution for f'
!==================================================================================================================================
FUNCTION GetNewFP(xi, steps, solut)
! MODULES

IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL, INTENT(IN)               :: xi                    ! Point in terms of similarity variable where f' should be interpolate to
REAL, INTENT(IN)               :: steps(:)              ! Vector containing the steps from Runge-Kutta method
REAL, INTENT(IN)               :: solut(:)              ! Solution f' for each point of the similarity variable in steps(:)
REAL                           :: GetNewFP
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: i, sz
REAL                           :: mx, fprime
!----------------------------------------------------------------------------------------------------------------------------------

! Naive, linear search will be used
! Obviously, this can get much better, although vector steps(:) tends to be small-sized
sz = SIZE(steps,1)
DO i=1,sz
    IF (steps(i).GE.xi) EXIT
END DO

IF (i.EQ.1) THEN
    fprime = solut(1)
ELSE
    mx = (xi - steps(i-1)) * (solut(i) - solut(i-1)) / (steps(i) - steps(i-1))
    fprime = mx + solut(i-1)
END IF

GetNewFP = fprime

END FUNCTION GetNewFP


!==================================================================================================================================
!> Evaluate u_tau from a log-law given instantaneous velocity, using Newton method
!==================================================================================================================================
FUNCTION NewtonLogLaw(velx,nu,abs_h_wm)
! MODULES
USE MOD_WMLES_Vars
USE MOD_Testcase_Vars       ,ONLY: dpdx
USE MOD_TimeDisc_Vars       ,ONLY: t
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
            LOGWRITE(*,'(3(E15.8,2X))') NewtonLogLaw, f, fprime
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

IF (iter.EQ.10) THEN 
    LOGWRITE(*,*) "NEWTON METHOD FAILED TO CONVERGE!"
END IF
IF (Logging) FLUSH(UNIT_logOut)


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
!> Solves the Falkner-Skan Equation using an iterative procedure
!> See Zhang, J., Chen, B. - An Iterative Method for Solving the Falkner-Skan Equation
!==================================================================================================================================
SUBROUTINE FalknerSkan(beta0, beta1, etainf, ddfddn, xis, dfdn, tol1, tol2)
! MODULES
USE MOD_WMLES_Vars
USE MOD_WMLES_Utils
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL, INTENT(IN)                            :: beta0        ! Parameter that defines the actual equation to be solved
REAL, INTENT(IN)                            :: beta1        ! Parameter that defines the actual equation to be solved
REAL, INTENT(OUT)                           :: etainf       ! Value of the free boundary, where fprime = 1, i.e. eta -> inf
REAL, INTENT(OUT)                           :: ddfddn       ! Value of the second derivative of f that matches the boundary conditions at eta -> 1
REAL, INTENT(OUT), OPTIONAL, ALLOCATABLE    :: xis(:)       ! Final time-steps xis
REAL, INTENT(OUT), OPTIONAL, ALLOCATABLE    :: dfdn(:)      ! Final solution for all "time-steps" xis
REAL, INTENT(IN), OPTIONAL                  :: tol1         ! tolerance for the convergence of |p|
REAL, INTENT(IN), OPTIONAL                  :: tol2         ! tolerance for the convergence of |q|
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
PROCEDURE(RHS), POINTER      :: func
REAL                         :: f, u, v, tol_p, tol_q, d_alpha, d23_alpha, d23_eta, l_s, r_s
REAL, ALLOCATABLE            :: xi_25(:), sol_25(:,:), xi_27(:), sol_27(:,:)
INTEGER                      :: N16, N25, N27, k, n, i
LOGICAL                      :: k_conv, n_conv, out_solution, out_tsteps
REAL, DIMENSION(3)           :: y0
REAL, DIMENSION(2)           :: rhs_23
!----------------------------------------------------------------------------------------------------------------------------------
! Set the parameters of the equation ("globally")
beta_0 = beta0
beta = beta1

! Set tolerances for convergence
tol_p = 1E-10
tol_q = 1E-6
out_tsteps = .FALSE.
out_solution = .FALSE.
IF (PRESENT(tol1)) tol_p = tol1
IF (PRESENT(tol2)) tol_q = tol2
IF (PRESENT(xis)) out_tsteps = .TRUE.
IF (PRESENT(dfdn)) out_solution = .TRUE.

! First guess of the global solution
alpha = 1.0
IF (beta .LT. 10) THEN
    eta_inf = 1.0
ELSE
    eta_inf = 0.1
END IF
k=0
k_conv = .FALSE.
DO WHILE (.NOT.k_conv)
    n_conv = .FALSE.
    n=0
    DO WHILE (.NOT.n_conv)
        func => RHS16
        y0(1) = 0.; y0(2) = 0.; y0(3) = alpha
        CALL SOLVE_ODE(func, 0., y0, 1., sol_16, xi_16, max_h=1./100.)
        func => RHS25
        y0(1) = 0.; y0(2) = 0.; y0(3) = 1.
        CALL SOLVE_ODE(func, 0., y0, 1., sol_25, xi_25, max_h=1./100., preStep=.TRUE.)
        N16 = SIZE(xi_16,1)
        N25 = SIZE(xi_25,1)
        d_alpha = - (sol_16(2,N16) - 1.) / sol_25(2,N25)
        alpha = alpha + d_alpha
        n = n+1
        IF (ABS(sol_16(2,N16)-1) .LE. tol_p) n_conv = .TRUE.
    END DO

    func => RHS27
    y0(1) = 0.; y0(2) = 0.; y0(3) = 0.
    CALL SOLVE_ODE(func, 0., y0, 1., sol_27, xi_27, max_h=1./100., preStep=.TRUE.)
    N27 = SIZE(xi_27,1)

    rhs_23(1) =  1. - sol_16(2,N16)
    rhs_23(2) = -1. * sol_16(3,N16)
    l_s = sol_25(2,N25) - (sol_25(3,N25) * sol_27(2,N27) / sol_27(3,N27))
    r_s = (-1. * rhs_23(2) * sol_27(2,N27) / sol_27(3,N27)) + rhs_23(1)
    d23_alpha = r_s / l_s

    l_s = sol_27(2,N27)
    r_s = rhs_23(1) - sol_25(2,N25)*d23_alpha
    d23_eta = r_s / l_s

    eta_inf = eta_inf + d23_eta
    alpha = alpha + d23_alpha
    k = k+1
    IF (ABS(sol_16(2,N16)-1).LE.tol_p .AND. ABS(sol_16(3,N16)).LE.tol_q) k_conv=.TRUE.
END DO

ddfddn = alpha
etainf = eta_inf

IF (out_tsteps) THEN
    SDEALLOCATE(xis)
    ALLOCATE(xis(N16))
    xis(:) = xi_16(:)
END IF
IF (out_solution) THEN
    SDEALLOCATE(dfdn)
    ALLOCATE(dfdn(N16))
    dfdn(:) = sol_16(2,:)
END IF


END SUBROUTINE FalknerSkan

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

IF (nWMLESRecvProcs.NE.0) THEN
    ! IF above avoid errors when indexing WMLESRecvRange -- even though when compiled in RELEASE mode everything works fine....
    DO i=1,WMLESRecvRange(2,nWMLESRecvProcs)
        HWMInfo(1, INT(HWMRecvInfo(nHWMPropSend+2, i)), & ! p
                    INT(HWMRecvInfo(nHWMPropSend+3, i)), & ! q
                    BCSideToWMLES(INT(HWMRecvInfo(nHWMPropSend+4, i))-offsetBCSides)) & ! iWMLESSide
                = HWMRecvInfo(nHWMPropSend+1, i) ! h_wm value

        DO k=1,nHWMPropSend
            HWMInfo(k+1, INT(HWMRecvInfo(nHWMPropSend+2, i)), & ! p
                    INT(HWMRecvInfo(nHWMPropSend+3, i)), & ! q
                    BCSideToWMLES(INT(HWMRecvInfo(nHWMPropSend+4, i))-offsetBCSides)) & ! iWMLESSide
                = HWMRecvInfo(k, i) ! flow property corresponding to index k
        END DO
    END DO
END IF

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

SDEALLOCATE(WMLESToLaminarSide)
SDEALLOCATE(WMLESToTurbulentSide)
SDEALLOCATE(FSBeta)
!SDEALLOCATE(FSDelta)
!SDEALLOCATE(FSEtaInf)
SDEALLOCATE(FSWallShear)
!SDEALLOCATE(FSEta)
!SDEALLOCATE(FSFPrime)

END SUBROUTINE FinalizeWMLES



!==================================================================================================================================
END MODULE MOD_WMLES
