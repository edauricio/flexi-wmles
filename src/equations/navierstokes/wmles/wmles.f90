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

CALL prms%CreateRealOption('vKarman', "von Karman constant", "0.41")
CALL prms%CreateRealOption('B', "Log-law intercept coefficient", "5.2")
CALL prms%CreateLogicalOption('UseSemiLocal', 'Set true to compute Wall Stress using information from the element above from the wall-adjacent element', '.TRUE.')
CALL prms%CreateIntOption('WMLES_NFilter','Number of high order modes to cut-off/attenuate while filtering the LES solution (< PP_N!)', '0')

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
USE MOD_Interpolation_Vars      ,ONLY: Vdm_Leg,sVdm_Leg
USE MOD_ReadInTools             ,ONLY: GETINTFROMSTR, GETREAL, GETLOGICAL, GETINT
USE MOD_StringTools             ,ONLY: STRICMP
USE MOD_Testcase_Vars           ,ONLY: Testcase
#if USE_MPI
USE MOD_MPI
USE MOD_MPI_Vars                !,ONLY:
#endif

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=10)               :: Channel="channel"
REAL                            :: hwmDistVec(3), Distance
INTEGER                         :: i,j,p,q,flip
INTEGER, ALLOCATABLE            :: WMLESToBCSide_tmp(:), WMLESFlip_tmp(:)
INTEGER                         :: iSide, WallElemID, OppSideID, OppLocSide
INTEGER,ALLOCATABLE             :: MasterToWMLESSide_tmp(:), SlaveToWMLESSide_tmp(:)
INTEGER,ALLOCATABLE             :: MasterToOppSide_tmp(:), SlaveToOppSide_tmp(:)
LOGICAL,ALLOCATABLE             :: IsMaster(:), IsSlave(:)
!==================================================================================================================================
!IF((.NOT.InterpolationInitIsDone).OR.(.NOT.MeshInitIsDone).OR.WMLESInitDone) THEN
!  CALL CollectiveStop(__STAMP__,&
!    'Wall-Modeled LES not ready to be called or already called.')
!END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT Wall-Modeled LES...'

WallModel = GETINTFROMSTR('WallModel')
UseSemiLocal = GETLOGICAL('UseSemiLocal')

! If Schumann's model is selected, check if we are in a channel flow
IF ((WallModel .EQ. WMLES_SCHUMANN) .AND. (.NOT.STRICMP(Testcase, Channel))) THEN
    CALL CollectiveStop(__STAMP__,&
        "Schumann's wall model can only be applied to the Channel flow testcase.")
END IF

vKarman = GETREAL('vKarman')
B = GETREAL('B')

WMLES_Filter = GETINT('WMLES_NFilter', '0')


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
ALLOCATE(MasterToWMLESSide_tmp(nBCSides),SlaveToWMLESSide_tmp(nBCSides))
ALLOCATE(MasterToOppSide_tmp(nBCSides),SlaveToOppSide_tmp(nBCSides))
ALLOCATE(WMLESFlip_tmp(nBCSides))

nWMLESSides = 0
nMasterWMLESSide = 0
nSlaveWMLESSide = 0
BCSideToWMLES = 0
WMLESToBCSide_tmp = 0
MasterToWMLESSide_tmp = 0
MasterToOppSide_tmp = 0
SlaveToWMLESSide_tmp = 0
SlaveToOppSide_tmp = 0
WMLESFlip_tmp = 0

IF (UseSemiLocal) THEN
    ! Use LES information from the element above the wall-adjacent element
    DO iSide=1,nBCSides
        IF (BoundaryType(BC(iSide),BC_TYPE) .EQ. 5) THEN ! WMLES side

            nWMLESSides = nWMLESSides + 1
            BCSideToWMLES(iSide) = nWMLESSides
            WMLESToBCSide_tmp(nWMLESSides) = iSide

            WallElemID = SideToElem(S2E_ELEM_ID, iSide)
            CALL GetOppositeSide(iSide, WallElemID, OppSideID)

            IF (SideToElem(S2E_ELEM_ID, OppSideID) .EQ. WallElemID) THEN 
                ! WallElem is master of opposite side
                ! So we use the info from slave, i.e., the element just above WallElem
                nSlaveWMLESSide = nSlaveWMLESSide + 1
                SlaveToWMLESSide_tmp(nSlaveWMLESSide) = nWMLESSides
                SlaveToOppSide_tmp(nSlaveWMLESSide) = OppSideID
            ELSE
                ! Otherwise, if WallElem is slave, we use the master info, i.e.,
                ! coming from the element just above WallElem
                nMasterWMLESSide = nMasterWMLESSide + 1
                MasterToWMLESSide_tmp(nMasterWMLESSide) = nWMLESSides
                MasterToOppSide_tmp(nMasterWMLESSide) = OppSideID
            END IF
        END IF
    END DO

ELSE
    ! In this case, we use information from the wall-adjacent element ONLY
    DO iSide=1,nBCSides
        IF (BoundaryType(BC(iSide),BC_TYPE) .EQ. 5) THEN ! WMLES side

            nWMLESSides = nWMLESSides + 1
            BCSideToWMLES(iSide) = nWMLESSides
            WMLESToBCSide_tmp(nWMLESSides) = iSide

            WallElemID = SideToElem(S2E_ELEM_ID, iSide)
            CALL GetOppositeSide(iSide, WallElemID, OppSideID, OppLocSide)

            IF (SideToElem(S2E_ELEM_ID, OppSideID) .EQ. WallElemID) THEN 
                ! WallElem is master of opposite side
                ! So we use the info from master variables
                nMasterWMLESSide = nMasterWMLESSide + 1
                MasterToWMLESSide_tmp(nMasterWMLESSide) = nWMLESSides
                MasterToOppSide_tmp(nMasterWMLESSide) = OppSideID

                ! Also, since we are taking information from the element-local opposite side, 
                ! we need to "flip" coordinates from the BC to the OppSide
                SELECT CASE(OppLocSide)
                    CASE (XI_MINUS,XI_PLUS,ZETA_MINUS,ZETA_PLUS)
                        WMLESFlip_tmp(nWMLESSides) = 1
                    CASE (ETA_MINUS,ETA_PLUS)
                        WMLESFlip_tmp(nWMLESSides) = 2
                END SELECT
            ELSE
                ! Otherwise, if WallElem is slave, we use the slave info
                nSlaveWMLESSide = nSlaveWMLESSide + 1
                SlaveToWMLESSide_tmp(nSlaveWMLESSide) = nWMLESSides
                SlaveToOppSide_tmp(nSlaveWMLESSide) = OppSideID

                ! We also need to "flip" coordinates for slave sides,
                ! but this info is directly available in SideToElem array
                WMLESFlip_tmp(nWMLESSides) = SideToElem(S2E_FLIP, OppSideID)
            END IF
        END IF
    END DO

END IF ! UseSemiLocal

!> Allocate permanent space and free temporary ones
ALLOCATE(MasterToWMLESSide(nMasterWMLESSide))
ALLOCATE(MasterToOppSide(nMasterWMLESSide))
ALLOCATE(SlaveToWMLESSide(nSlaveWMLESSide))
ALLOCATE(SlaveToOppSide(nSlaveWMLESSide))
ALLOCATE(WMLESToBCSide(nWMLESSides))
ALLOCATE(WMLESFlip(nWMLESSides))
IF (Logging) THEN
    ALLOCATE(IsMaster(nWMLESSides),IsSlave(nWMLESSides))
    IsMaster = .FALSE.
    IsSlave = .FALSE.
END IF

DO i=1,nMasterWMLESSide
    MasterToWMLESSide(i) = MasterToWMLESSide_tmp(i)
    MasterToOppSide(i) = MasterToOppSide_tmp(i)
    IF (Logging) IsMaster(MasterToWMLESSide(i)) = .TRUE.
END DO
DO i=1,nSlaveWMLESSide
    SlaveToWMLESSide(i) = SlaveToWMLESSide_tmp(i)
    SlaveToOppSide(i) = SlaveToOppSide_tmp(i)
    IF (Logging) IsSlave(SlaveToWMLESSide(i)) = .TRUE.
END DO
DO i=1,nWMLESSides
    WMLESToBCSide(i) = WMLESToBCSide_tmp(i)
    WMLESFlip(i) = WMLESFlip_tmp(i)
END DO
SDEALLOCATE(WMLESToBCSide_tmp)
SDEALLOCATE(MasterToWMLESSide_tmp)
SDEALLOCATE(SlaveToWMLESSide_tmp)
SDEALLOCATE(MasterToOppSide_tmp)
SDEALLOCATE(SlaveToOppSide_tmp)
SDEALLOCATE(WMLESFlip_tmp)

! Opposite side's Slave MPI proc responsible for BC imposition needs Face_xGP data 
! to calculate the wall model height, h_wm
CALL StartSendMPIData(Face_xGP, 3*(PP_N+1)*(PP_NZ+1), 1, nSides, MPIRequest_U(:,SEND), SendID=1)
CALL StartReceiveMPIData(Face_xGP, 3*(PP_N+1)*(PP_NZ+1), 1, nSides, MPIRequest_U(:,RECV), SendID=1)
CALL FinishExchangeMPIData(2*nNbProcs, MPIRequest_U)

ALLOCATE(h_wm(0:PP_N,0:PP_NZ,nWMLESSides))

IF (nWMLESSides.GT.0) THEN
    IF (UseSemiLocal) THEN
        DO i=1,nMasterWMLESSide
            ! No flip needed, we're taking info from the element above the wall-adjacent element
            ! which has the same coordinate orientation than the BC side
            DO q=0,PP_NZ; DO p=0,PP_N
                hwmDistVec(1:3) = Face_xGP(1:3,p,q,0,MasterToOppSide(i)) - Face_xGP(1:3,p,q,0,WMLESToBCSide(MasterToWMLESSide(i)))
                Distance = 0
                DO j=1,3
                    Distance = Distance + hwmDistVec(j)**2
                END DO
                h_wm(p,q,MasterToWMLESSide(i)) = SQRT(Distance)
            END DO; END DO ! p, q                
        END DO

        DO i=1,nSlaveWMLESSide
            ! Here we need flip, because since the element above is slave of face,
            ! it means wall-adjacent element is master of face, and Face_xGP info
            ! comes from master
            OppLocSide = SideToElem(S2E_LOC_SIDE_ID, SlaveToOppSide(i))
            IF (OppLocSide.EQ.-1) CALL Abort(__STAMP__, "OPPPSS, -1")
            SELECT CASE(OppLocSide)
                CASE (XI_MINUS,XI_PLUS,ZETA_MINUS,ZETA_PLUS)
                    flip = 1
                CASE (ETA_MINUS,ETA_PLUS)
                    flip = 2
            END SELECT
            DO q=0,PP_NZ; DO p=0,PP_N
                hwmDistVec(1:3) = Face_xGP(1:3,FS2M(1,p,q,flip),FS2M(2,p,q,flip),0,SlaveToOppSide(i)) - Face_xGP(1:3,p,q,0,WMLESToBCSide(SlaveToWMLESSide(i)))
                Distance = 0
                DO j=1,3
                    Distance = Distance + hwmDistVec(j)**2
                END DO
                h_wm(p,q,SlaveToWMLESSide(i)) = SQRT(Distance)
            END DO; END DO ! p, q            
        END DO

    ELSE
        DO i=1,nMasterWMLESSide
            ! Here a flip is needed, because opposite side is Master of interface, and Face_xGP
            ! is a master information
            OppLocSide = SideToElem(S2E_LOC_SIDE_ID, MasterToOppSide(i))
            IF (OppLocSide.EQ.-1) CALL Abort(__STAMP__, "OPPPSS, -1")
            SELECT CASE(OppLocSide)
                CASE (XI_MINUS,XI_PLUS,ZETA_MINUS,ZETA_PLUS)
                    flip = 1
                CASE (ETA_MINUS,ETA_PLUS)
                    flip = 2
            END SELECT
            DO q=0,PP_NZ; DO p=0,PP_N
                hwmDistVec(1:3) = Face_xGP(1:3,FS2M(1,p,q,flip),FS2M(2,p,q,flip),0,MasterToOppSide(i)) - Face_xGP(1:3,p,q,0,WMLESToBCSide(MasterToWMLESSide(i)))
                Distance = 0
                DO j=1,3
                    Distance = Distance + hwmDistVec(j)**2
                END DO
                h_wm(p,q,MasterToWMLESSide(i)) = SQRT(Distance)
            END DO; END DO ! p, q
        END DO

        DO i=1,nSlaveWMLESSide
            ! No flip is needed here, since wall-adjacent is slave of face and
            ! Face_xGP info comes from the master (element above), which has
            ! the same orientation already
            DO q=0,PP_NZ; DO p=0,PP_N
                hwmDistVec(1:3) = Face_xGP(1:3,p,q,0,SlaveToOppSide(i)) - Face_xGP(1:3,p,q,0,WMLESToBCSide(SlaveToWMLESSide(i)))
                Distance = 0
                DO j=1,3
                    Distance = Distance + hwmDistVec(j)**2
                END DO
                h_wm(p,q,SlaveToWMLESSide(i)) = SQRT(Distance)
            END DO; END DO ! p, q
        END DO
    END IF
END IF

IF (nWMLESSides.GT.0 .AND. WMLES_Filter.GT.0) THEN
    IF (WMLES_Filter .GT. PP_N) CALL CollectiveStop(__STAMP__, "WMLES_Filter must be <= N")
    ALLOCATE(WMLES_FilterMat(0:PP_N, 0:PP_N))
    WMLES_FilterMat = 0.
    ! Using a simple cut-off filter
    DO i=0,(PP_N-WMLES_Filter)
        WMLES_FilterMat(i,i) = 1.
    END DO
    WMLES_FilterMat = MATMUL(MATMUL(Vdm_Leg,WMLES_FilterMat),sVdm_Leg)
END IF

ALLOCATE(WMLES_TauW(2,0:PP_N,0:PP_NZ,nWMLESSides))
WMLES_TauW = 0.

!> Display debug information in Logfile
IF (nWMLESSides.GT.0) THEN
    LOGWRITE(*,'(/,A72)') '========================================================================'
    LOGWRITE(*,'(20X,A32,20X)') 'WMLES Pre-Computation Stage Info'
    LOGWRITE(*,'(A72)') '========================================================================'
    LOGWRITE(*,'(10X,A9)') "Summary:"
    LOGWRITE(*,'(20X,A18,X,I6)') "nWMLESSides:", nWMLESSides
    LOGWRITE(*,'(20X,A18,X,I6)') "nMasterWMLESSide:", nMasterWMLESSide
    LOGWRITE(*,'(20X,A18,X,I6)') "nSlaveWMLESSide:", nSlaveWMLESSide
    LOGWRITE(*,'(A72)') '------------------------------------------------------------------------'
    LOGWRITE(*,'(2(A4,2X),A10,2X,2(A6,2X),A9,2X,A15)') "p","q","nWMLESSide","Master","Slave","WMLESFlip", "h_wm Distance"
    DO i=1,nWMLESSides
        DO q=0,PP_NZ; DO p=0,PP_N
            LOGWRITE(*,'(2(I4,2X),I10,2X,2(L6,2X),I9,2XE15.8)') p, q, i, IsMaster(i), IsSlave(i), WMLESFlip(i), h_wm(p,q,i)
        END DO; END DO
    END DO
    LOGWRITE(*,'(A18,X,A35,X,A18,/)') '==================','END OF WMLES PRE-COMPUT. STAGE INFO','=================='
    IF (Logging) CALL FLUSH(UNIT_logOut)
END IF

WMLESInitDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT Wall-Modeled LES DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitWMLES


!==================================================================================================================================
!> Get the ID of the side opposite of that indicated by SideID, in the same element
!==================================================================================================================================
SUBROUTINE GetOppositeSide(SideID, ElemID, OppositeSide, SideFlip)
! MODULES
USE MOD_Mesh_Vars,              ONLY: SideToElem, ElemToSide
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, INTENT(IN)                 :: ElemID, SideID
INTEGER, INTENT(OUT)                :: OppositeSide
INTEGER, INTENT(OUT), OPTIONAL      :: SideFlip
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
    IF (PRESENT(SideFlip)) SideFlip = ETA_PLUS
CASE (ETA_PLUS)
    OppositeSide = ElemToSide(E2S_SIDE_ID, ETA_MINUS, ElemID)
    IF (PRESENT(SideFlip)) SideFlip = ETA_MINUS
CASE (XI_MINUS)
    OppositeSide = ElemToSide(E2S_SIDE_ID, XI_PLUS, ElemID)
    IF (PRESENT(SideFlip)) SideFlip = XI_PLUS
CASE (XI_PLUS)
    OppositeSide = ElemToSide(E2S_SIDE_ID, XI_MINUS, ElemID)
    IF (PRESENT(SideFlip)) SideFlip = XI_MINUS
CASE (ZETA_MINUS)
    OppositeSide = ElemToSide(E2S_SIDE_ID, ZETA_PLUS, ElemID)
    IF (PRESENT(SideFlip)) SideFlip = ZETA_PLUS
CASE (ZETA_PLUS)
    OppositeSide = ElemToSide(E2S_SIDE_ID, ZETA_MINUS, ElemID)
    IF (PRESENT(SideFlip)) SideFlip = ZETA_MINUS
END SELECT

END SUBROUTINE


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
USE MOD_MPI
USE MOD_MPI_Vars
USE MOD_ChangeBasisByDim            ,ONLY: ChangeBasisSurf
IMPLICIT NONE

!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                             :: p,q,i,j,SideID,flip
REAL                                :: u_tau, tau_w_mag, utang, VelMag
REAL                                :: tangvec(3), tau_w_vec(3), FaceData(PP_nVarPrim,0:PP_N,0:PP_NZ)
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

! Compute tau_w for points of our responsibility
SELECT CASE(WallModel)

  CASE(WMLES_SCHUMANN)
    
    ! Schumann model
    ! Computation is done using instantaneous information, projected onto the wall-tangent direction
    ! Mean wall shear stress is an INPUT to this model ! Therefore, only works for periodic channel flows.
    ! The channel flow testcase assumes <tau_w> = 1.0
    ! u_tau is used in the log-law eq. to extract <u> (u_mean), which is then used to calculate tau_w

    ! DO sProc=0,nProcs_SendTauW
    !     ! Calculate tau_w for each h_wm that is approximated as a face node
    !     DO FPInd=1,nTauW_MINE_FacePoint(sProc)
    !         p = TauW_MINE_FacePoint(1,FPInd,sProc)
    !         q = TauW_MINE_FacePoint(2,FPInd,sProc)
    !         SideID = TauW_MINE_FacePoint(3,FPInd,sProc)

    !         ! Face_xGP is only populated for master sides, so that only master sides have approximation nodes (check InitWMLES above)
    !         ! hence, use of UPrim_master is guaranteed to be correct here
    !         !utang = DOT_PRODUCT(UPrim_master(2:4,p,q,SideID),TauW_MINE_NormVec(:,FaceToLocalPoint(FPInd,sProc),sProc))
    !         utang = UPrim_master(2,p,q,SideID)

    !         u_tau = SQRT(1.0/UPrim_master(1,p,q,SideID)) 
    !         u_mean = u_tau*( (1./vKarman) * LOG((abs_h_wm*u_tau)/mu0) + B )
    !         TauW_MINE(1,FaceToLocalPoint(FPInd,sProc),sProc) = utang/u_mean * 1.0 ! <tau_w> = 1.0
    !         TauW_MINE(2,FaceToLocalPoint(FPInd,sProc),sProc) = (2.0*mu0)*(UPrim_master(4,p,q,SideID)/abs_h_wm)
    !     END DO

    !     ! Calculate tau_w for each h_wm that is approximated as an interior node
    !     DO IPInd=1,nTauW_MINE_InteriorPoint(sProc)
    !         p = TauW_MINE_InteriorPoint(1,IPInd,sProc)
    !         q = TauW_MINE_InteriorPoint(2,IPInd,sProc)
    !         r = TauW_MINE_InteriorPoint(3,IPInd,sProc)
    !         ElemID = TauW_MINE_InteriorPoint(4,IPInd,sProc)

    !         !utang = DOT_PRODUCT(UPrim(2:4,p,q,r,ElemID),TauW_MINE_NormVec(:,InteriorToLocalPoint(IPInd),sProc))
    !         utang = UPrim(2,p,q,r,ElemID)

    !         u_tau = SQRT(1.0/UPrim(1,p,q,r,ElemID)) 
    !         u_mean = u_tau*( (1./vKarman) * LOG((abs_h_wm*u_tau)/mu0) + B )
    !         TauW_MINE(1,InteriorToLocalPoint(IPInd,sProc),sProc) = utang/u_mean * 1.0 ! <tau_w> = 1.0
    !         TauW_MINE(2,InteriorToLocalPoint(IPInd,sProc),sProc) = (2.0*mu0)*(UPrim(4,p,q,r,ElemID)/abs_h_wm)
    !     END DO

    !     ! Calculate tau_w for each h_wm that must be interpolated
    !     DO IntPInd=1,nTauW_MINE_Interpolate(sProc)
    !         ElemID = INT(TauW_MINE_Interpolate(4,IntPInd,sProc))
    !         !vel_inst(1) = InterpolateHwm(ElemID,Lag_xi(:,IntPInd,sProc),Lag_eta(:,IntPInd,sProc),Lag_zeta(:,IntPInd,sProc),2,prim=.TRUE.)
    !         !vel_inst(2) = InterpolateHwm(ElemID,Lag_xi(:,IntPInd,sProc),Lag_eta(:,IntPInd,sProc),Lag_zeta(:,IntPInd,sProc),3,prim=.TRUE.)
    !         !vel_inst(3) = InterpolateHwm(ElemID,Lag_xi(:,IntPInd,sProc),Lag_eta(:,IntPInd,sProc),Lag_zeta(:,IntPInd,sProc),4,prim=.TRUE.)
    !         !rho_inst = InterpolateHwm(ElemID,Lag_xi(:,IntPInd,sProc),Lag_eta(:,IntPInd,sProc),Lag_zeta(:,IntPInd,sProc),1)

    !         !utang = DOT_PRODUCT(vel_inst(:),TauW_MINE_NormVec(:,InterpToLocalPoint(IntPInd),sProc))
    !         utang = vel_inst(1)

    !         u_tau = SQRT(1.0/rho_inst) 
    !         u_mean = u_tau*( (1./vKarman) * LOG((abs_h_wm*u_tau)/mu0) + B )
    !         TauW_MINE(1,InterpToLocalPoint(IntPInd,sProc),sProc) = utang/u_mean * 1.0 ! <tau_w> = 1.0
    !         TauW_MINE(2,InterpToLocalPoint(IntPInd,sProc),sProc) = (2.0*mu0)*(vel_inst(3)/abs_h_wm)
    !     END DO
    ! END DO

  CASE (WMLES_LOGLAW)

    ! In this case, it is assumed that the log-law is an instantaneous relation
    ! Then, u (projected onto the wall-tangent streamwise direction) is used in the log-law
    ! and u_tau is computed (using Newton iterations).
    ! Then, u_tau is used to compute wall shear stresses.
    DO i=1,nMasterWMLESSide
        SideID = MasterToOppSide(i)

        ! Whenever we are using information from the local side, we need to flip the reference
        ! so that the indices p,q, coincide between the BC and the upper face.
        ! This is not needed when we're taking non-local information, i.e., from the
        ! face owned by the element just above the wall element
        flip = WMLESFlip(MasterToWMLESSide(i))

        ! If filtering is ON, we should filter the face solution now and store it in LES_Sol array
        IF (WMLES_Filter .GT. 0) THEN
            !CALL WMLESFilter(UPrim_master(:,:,:,SideID), FaceData, WMLES_FilterMat)
            CALL ChangeBasisSurf(PP_nVarPrim,PP_N,PP_N,WMLES_FilterMat,UPrim_master(:,:,:,SideID),FaceData)
        ELSE
            FaceData(:,:,:) = UPrim_master(:,:,:,SideID)
        END IF

        DO q=0,PP_NZ; DO p=0,PP_N
            tangvec = FaceData(2:4,FS2M(1,p,q,flip),FS2M(2,p,q,flip)) - DOT_PRODUCT(FaceData(2:4,FS2M(1,p,q,flip),FS2M(2,p,q,flip)),NormVec(1:3,p,q,0,WMLESToBCSide(MasterToWMLESSide(i))))*NormVec(1:3,p,q,0,WMLESToBCSide(MasterToWMLESSide(i)))
            VelMag = 0.
            DO j=1,3
                VelMag = VelMag + tangvec(j)**2
            END DO
            VelMag = SQRT(VelMag)
            tangvec = tangvec/VelMag

            utang = DOT_PRODUCT(FaceData(2:4,FS2M(1,p,q,flip),FS2M(2,p,q,flip)),tangvec)

            u_tau = NewtonLogLaw(utang,(mu0/FaceData(1,FS2M(1,p,q,flip),FS2M(2,p,q,flip))),h_wm(p,q,MasterToWMLESSide(i)))
            tau_w_mag = FaceData(1,FS2M(1,p,q,flip),FS2M(2,p,q,flip))*(u_tau**2) ! CHECK TODO ABOVE
            tau_w_vec = (/0.,tau_w_mag,0./)

            WMLES_Tauw(1,p,q,MasterToWMLESSide(i)) = -1.*DOT_PRODUCT(tau_w_vec(1:3),NormVec(1:3,p,q,0,WMLESToBCSide(MasterToWMLESSide(i))))
            WMLES_Tauw(2,p,q,MasterToWMLESSide(i)) = 0.
        END DO; END DO ! p, q
    END DO

    DO i=1,nSlaveWMLESSide
        SideID = SlaveToOppSide(i)

        ! Whenever we are using information from the local side, we need to flip the reference
        ! so that the indices p,q, coincide between the BC and the upper face.
        ! This is not needed when we're taking non-local information, i.e., from the
        ! face owned by the element just above the wall element
        flip = WMLESFlip(SlaveToWMLESSide(i))

        IF (WMLES_Filter .GT. 0) THEN
            !CALL WMLESFilter(UPrim_slave, FaceData, WMLES_FilterMat)
            CALL ChangeBasisSurf(PP_nVarPrim,PP_N,PP_N,WMLES_FilterMat,UPrim_slave(:,:,:,SideID),FaceData)
        ELSE
            FaceData(:,:,:) = UPrim_slave(:,:,:,SideID)
        END IF

        DO q=0,PP_NZ; DO p=0,PP_N
            tangvec = FaceData(2:4,FS2M(1,p,q,flip),FS2M(2,p,q,flip)) - DOT_PRODUCT(FaceData(2:4,FS2M(1,p,q,flip),FS2M(2,p,q,flip)),NormVec(1:3,p,q,0,WMLESToBCSide(SlaveToWMLESSide(i))))*NormVec(1:3,p,q,0,WMLESToBCSide(SlaveToWMLESSide(i)))
            VelMag = 0.
            DO j=1,3
                VelMag = VelMag + tangvec(j)**2
            END DO
            VelMag = SQRT(VelMag)
            tangvec = tangvec/VelMag

            utang = DOT_PRODUCT(FaceData(2:4,FS2M(1,p,q,flip),FS2M(2,p,q,flip)),tangvec)

            u_tau = NewtonLogLaw(utang,(mu0/FaceData(1,FS2M(1,p,q,flip),FS2M(2,p,q,flip))),h_wm(p,q,SlaveToWMLESSide(i)))
            tau_w_mag = FaceData(1,FS2M(1,p,q,flip),FS2M(2,p,q,flip))*(u_tau**2) ! CHECK TODO ABOVE
            tau_w_vec = (/0.,tau_w_mag,0./)

            WMLES_Tauw(1,p,q,SlaveToWMLESSide(i)) = -1.*DOT_PRODUCT(tau_w_vec(1:3),NormVec(1:3,p,q,0,WMLESToBCSide(SlaveToWMLESSide(i))))
            WMLES_Tauw(2,p,q,SlaveToWMLESSide(i)) = 0.

        END DO; END DO ! p, q
    END DO

    IF(Logging) CALL FLUSH(UNIT_logOut)

  CASE(WMLES_COUETTE)

    ! Test for a laminar flow with known velocity profile, so that we can debug the implementation

    DO i=1,nMasterWMLESSide
        SideID = MasterToOppSide(i)

        DO q=0,PP_NZ; DO p=0,PP_N
            ! Calculated as if a model was employed
            utang = UPrim_master(2,p,q,SideID)

            u_tau = NewtonCouette(utang,mu0,h_wm(p,q,MasterToWMLESSide(i))) ! dpdx, obviously.

            ! Create vector aligned with wall-normal direction and magnitude of tau_xy
            tau_w_vec = (/0.,-0.5*u_tau,0./)
            ! Project it onto the normal direction so that the correct signal is imposed
            WMLES_Tauw(1,p,q,MasterToWMLESSide(i)) = -1.*DOT_PRODUCT(tau_w_vec(1:3),NormVec(1:3,p,q,0,WMLESToBCSide(MasterToWMLESSide(i))))
            WMLES_Tauw(2,p,q,MasterToWMLESSide(i)) = 0.

        END DO; END DO ! p, q
    END DO

    DO i=1,nSlaveWMLESSide
        SideID = SlaveToOppSide(i)

        DO q=0,PP_NZ; DO p=0,PP_N
            ! Calculated as if a model was employed
            utang = UPrim_slave(2,p,q,SideID)

            u_tau = NewtonCouette(utang,mu0,h_wm(p,q,SlaveToWMLESSide(i))) ! dpdx, obviously.

            ! Create vector aligned with wall-normal direction and magnitude of tau_xy
            tau_w_vec = (/0.,-0.5*u_tau,0./)
            ! Project it onto the normal direction so that the correct signal is imposed
            WMLES_Tauw(1,p,q,SlaveToWMLESSide(i)) = -1.*DOT_PRODUCT(tau_w_vec(1:3),NormVec(1:3,p,q,0,WMLESToBCSide(SlaveToWMLESSide(i))))
            WMLES_Tauw(2,p,q,SlaveToWMLESSide(i)) = 0.

        END DO; END DO ! p, q
    END DO


  CASE DEFAULT
    CALL abort(__STAMP__,&
         'Unknown definition of Wall Model.')

END SELECT

! Start non-blockingly sending WMLES info (for those who have anything to send)
!CALL StartSendTauWMPI()

! Finish receiving WMLES info (now this is a blocking op.)
!CALL FinishReceiveTauWMPI()


! Set up WMLES_TauW vector with all necessary, computed values of wall shear stress
! DO iSide=1,nWMLESSides
!     DO q=0,PP_NZ; DO p=0,PP_N
!         WMLES_TauW(:,p,q,iSide) = TauW_YOURS(:,TauW_Proc(2,p,q,iSide),Proc_RecvTauW_Inv(TauW_Proc(1,p,q,iSide)))
!     END DO; END DO ! p, q
! END DO

!### DEBUG -- REMOVE
! IF (Logging) THEN
!     LOGWRITE(*,*) "=-=-=-=-=-=-=-=-=-=- CHECK POINTS EQUIVALENCE =-=-=-=-=-=-=-=-=-=-"
!     LOGWRITE(*,'((I4,2X),3(I15,2X))') "N#", "X", "Y", "Z"
!     DO iSide=1,nWMLESSides
!         DO q=0,PP_NZ; DO p=0,PP_N
!             LOGWRITE(*,'((I4,2X),3(E15.8,2X))') 
!         END DO; END DO ! p, q
!     END DO
! END IF


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
!> Evaluate u_tau from a log-law given instantaneous velocity, using Newton method
!==================================================================================================================================
FUNCTION NewtonLogLaw(velx,nu,hwm)
! MODULES
USE MOD_WMLES_Vars                  
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL, INTENT(IN)            :: velx         ! Tangential velocity to feed wall model (log-law function)
REAL, INTENT(IN)            :: nu           ! kinematic viscosity at h_wm
REAL, INTENT(IN)            :: hwm          ! Wall model height
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
            f = NewtonLogLaw*( (1./vKarman)*LOG(hwm*NewtonLogLaw/nu) + B ) - velx
            fprime = (1./vKarman) * (LOG(hwm*NewtonLogLaw/nu) + 1.) + B
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
FUNCTION NewtonCouette(velx,nu,hwm)
! MODULES
USE MOD_WMLES_Vars                  
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL, INTENT(IN)            :: velx         ! Tangential velocity to feed wall model (log-law function)
REAL, INTENT(IN)            :: nu           ! kinematic viscosity at h_wm
REAL, INTENT(IN)            :: hwm          ! Wall model height
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
    f = velx + (1./(2.*nu))*NewtonCouette*hwm*(1.-hwm)
    fprime = (1./(2.*nu))*hwm*(1.-hwm)
    ! IF(ABS(fprime).LE.1.E-5) EXIT ! fprime ~ 0 -- INCOMING OVERFLOW BY DIVISION
    NewtonCouette = NewtonCouette - f/fprime
    IF (ABS(f/fprime).LE.1.E-5) EXIT ! 1.E-5 = stop criterion (tolerance)
END DO
END FUNCTION NewtonCouette

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
SDEALLOCATE(h_wm)
SDEALLOCATE(WMLESFlip)
SDEALLOCATE(MasterToWMLESSide)
SDEALLOCATE(SlaveToWMLESSide)
SDEALLOCATE(MasterToOppSide)
SDEALLOCATE(SlaveToOppSide)
SDEALLOCATE(WMLES_FilterMat)

END SUBROUTINE FinalizeWMLES



!==================================================================================================================================
END MODULE MOD_WMLES
