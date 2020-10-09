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
REAL                            :: WMLES_Tol
INTEGER                         :: SideID, WSideID, iSide
INTEGER                         :: hwmElemID, OppSideID, OppSideFlip, WLocSide
INTEGER                         :: TotalNWMLESSides
REAL                            :: ElemHeight, DiffVect(1:3), DistanceFromPoint
INTEGER                         :: hOnOppFace, nextElem, firstElem
INTEGER                         :: i,j,k,p,q,r
INTEGER, ALLOCATABLE            :: tmpMasterToTSide(:), tmpMasterToWMLES(:), tmpSlaveToTSide(:), tmpSlaveToWMLES(:)
INTEGER                         :: nTauW_MINE, LocSide
REAL, ALLOCATABLE               :: tmpTauW_NormVec_MINE(:,:), tmpTauW_FacexGP_MINE(:,:)
INTEGER, ALLOCATABLE            :: tmpTauW_Element_MINE(:)
REAL                            :: h_wm_Coords(3), DistanceVect(3), InwardNorm(3), Distance
REAL                            :: abs_h_wm
LOGICAL                         :: Found

#if USE_MPI
REAL, ALLOCATABLE               :: dataToSend(:,:), dataToRecv(:)
INTEGER, ALLOCATABLE            :: commCntSend(:), commCntRcvd(:), commRequests(:), maxCommRecv(:)
INTEGER                         :: iErr, iStat(MPI_STATUS_SIZE), sendCnt
INTEGER                         :: iNbProc
LOGICAL                         :: WaitComm(nNbProcs), FoundOrNomB
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

! Count how many WMLES BC Sides we have defined,
! and populate the mappings:
! - Mesh SideID to WMLES SideID (For use in GetBoundaryFlux.f90 in conjunction with WMLES_Tauw)
! - WMLES SideID to TopSideID master/slave (side where info is gathered) (For use in ComputeWallStress)

! REMEMBER: It does not matter wheter the information (i.e. adjacent element) is on this or
! another MPI process. By the time we call ComputeWallStress, all the face information needed
! has already been communicated (step 3 of DG_WeakTimeDerivative).

! However, IF a mapping between side -> neighbor process rank is needed for further communication
! (e.g. when solving the TBL equations on the same mesh with p-refinement), it may be
! accomplished by checking SideToElem(S2E_(NB_)ELEM_ID, OppSideID) == -1 (i.e., the needed elem is not on this proc)

ALLOCATE(WMLES_Side(1:nBCSides))

ALLOCATE(tmpMasterToTSide(1:nBCSides))
ALLOCATE(tmpSlaveToTSide(1:nBCSides))
ALLOCATE(tmpMasterToWMLES(1:nBCSides))
ALLOCATE(tmpSlaveToWMLES(1:nBCSides))
WMLES_Side = 0
tmpMasterToTSide = 0
tmpSlaveToTSide = 0
tmpMasterToWMLES = 0
tmpSlaveToWMLES = 0

nWMLESSides = 0
nSlaveSides = 0
nMasterSides = 0
!---
ALLOCATE(tmpTauW_NormVec_MINE(1:3,(PP_N+1)*(PP_NZ+1)*nSides))
ALLOCATE(tmpTauW_FacexGP_MINE(1:3,(PP_N+1)*(PP_NZ+1)*nSides))
ALLOCATE(tmpTauW_Element_MINE((PP_N+1)*(PP_NZ+1)*nSides))

tmpTauW_NormVec_MINE = 0
tmpTauW_FacexGP_MINE = 0
tmpTauW_Element_MINE = 0

nTauW_MINE = 0
!---
! Create and set some temp MPI variables
#if USE_MPI
ALLOCATE(dataToSend(10,(PP_N+1)*(PP_NZ+1)*nMPISides))
ALLOCATE(dataToRecv(10))
ALLOCATE(commCntSend(nNbProcs))
ALLOCATE(commCntRcvd(nNbProcs))
ALLOCATE(commRequests(nNbProcs))
ALLOCATE(maxCommRecv(nNbProcs))

commCntSend = 0
commCntRcvd = 0
maxCommRecv = 0
commRequests = MPI_REQUEST_NULL
sendCnt = 0
#endif

LOGWRITE(*,*) "====================== WMLES Exchange Location Info ===================="
DO WSideID=1,nBCSides
    IF (BoundaryType(BC(WSideID),BC_TYPE) .EQ. 5) THEN ! WMLES side
        nWMLESSides = nWMLESSides + 1
        WMLES_Side(WSideID) = nWMLESSides


        LOGWRITE(*,*) "====== Side ", SideID, " ======="
        LOGWRITE(*,'(2(A4,2X),(A15,2X))') 'p', 'q', 'Within Element'
        DO q=0, PP_NZ; DO p=0, PP_N
            SideID = WSideID
            CALL FindHwmElement(p, q, SideID, hwmElemID)
#if USE_MPI
            ! Check if h_wm for this p,q is within the element of another MPI proc
            IF (((hwmElemID - nElems) .GT. 0) .AND. ((hwmElemID - nElems) .LE. nNbProcs)) THEN ! hwmElem is in another proc
                sendCnt = sendCnt + 1
                ! Mount data to be sent
                dataToSend(1,sendCnt) = myRank ! Send my rank so the other MPI proc knowns whom to return to
                dataToSend(2,sendCnt) = SideToGlobalSide(SideID) ! Global SideID of "lower side" at adj. element, so that the receiving MPI proc.
                                                                            ! can get his elemID, in order to get the oppositeSide and thus its info (face_xgp etc)
                                                                            ! Unfortunately, this info is all tha THIS proc. have.
                dataToSend(3:4,sendCnt) = (/p,q/)
                dataToSend(5:7,sendCnt) = Face_xGP(:,p,q,0,WSideID) ! Position vector of p,q at wall side
                dataToSend(8:10,sendCnt) = NormVec(:,p,q,0,WSideID) ! Normal vector of p,q at wall side

                CALL MPI_ISEND(dataToSend(:,sendCnt),10, MPI_DOUBLE_PRECISION, NbProc(hwmElemID-nElems), &
                    0, MPI_COMM_FLEXI, commRequests(hwmElemID-nElems), iErr)
                commCntSend(hwmElemID-nElems) = commCntSend(hwmElemID-nElems)+1
            ELSE ! Element is within my partition
#endif

            ! Mark this point's calculation (p,q,WSideID) as "my responsibility"
            ! Also, store the NormVec and Face_xGP positions so that we can, once again,
            ! see where in the physical domain h_wm lies, to check if we can approximate it
            ! as any point inside this element, or if we need to interpolate.
            nTauW_MINE = nTauW_MINE+1
            tmpTauW_NormVec_MINE(:,nTauW_MINE) = NormVec(:,p,q,0,WSideID)
            tmpTauW_FacexGP_MINE(:,nTauW_MINE) = Face_xGP(:,p,q,0,WSideID)
            tmpTauW_Element_MINE(nTauW_MINE) = hwmElemID

#if USE_MPI
            END IF
#endif  
        END DO; END DO ! p, q
    END IF
END DO
LOGWRITE(*,*) "====================== WMLES Exchange Location Info END ===================="
IF (Logging) FLUSH(UNIT_logOut)


#if USE_MPI
WaitComm = .TRUE.

! WMLES imposition procs send the "ending comm" message to each neighboring MPI proc
! and will only wait for messages different from that on the WHILE loop below (that is, different data/count/tag)
IF (nWMLESSides .NE. 0) THEN
    WaitComm = .FALSE.
    DO i=1,nNbProcs
        sendCnt = sendCnt + 1 ! Just to make sure we do not modify any other buffer sent before
        dataToSend(1,sendCnt) = myRank
        dataToSend(2,sendCnt) = -(commCntSend(i)+1)
        CALL MPI_ISEND(dataToSend(:,sendCnt), 10, MPI_DOUBLE_PRECISION, NbProc(i), 0, MPI_COMM_FLEXI, commRequests(i), iErr)
    END DO
END IF

DO WHILE (ANY(WaitComm))
    CALL MPI_RECV(dataToRecv(:), 10, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 0, MPI_COMM_FLEXI, iStat, iErr)
    ! Map rank of sender to my local neighboring proc list
    DO i=1,nNbProcs
        IF (INT(dataToRecv(1)) .EQ. NbProc(i)) THEN
            iNbProc=i
            EXIT
        END IF
    END DO
    commCntRcvd(iNbProc) = commCntRcvd(iNbProc) + 1

    IF (INT(dataToRecv(2)) .LT. 0) THEN ! IF to check the ending comm WITH THIS PROC
        maxCommRecv(iNbProc) = ABS(INT(dataToRecv(2))) ! Number of maximum comms to expect from this proc

        IF (maxCommRecv(iNbProc) .EQ. commCntRcvd(iNbProc)) THEN ! Reached max cons: do not wait anymore
            WaitComm(iNbProc) = .FALSE.
        END IF

        CYCLE
    END IF

    IF (maxCommRecv(iNbProc) .EQ. commCntRcvd(iNbProc)) THEN ! Reached max cons: do not wait anymore
        WaitComm(iNbProc) = .FALSE.
    END IF

    ! Find my SideID which maps to the received globalSideID
    DO i=1,nSides
        IF (SideToGlobalSide(i) .EQ. INT(dataToRecv(2))) THEN
            SideID = i
            EXIT
        END IF
    END DO

    hwmElemID = SideToElem(S2E_ELEM_ID, SideID)
    IF (hwmElemID .LT. 0) THEN ! This partition is actually slave for this side, not master
        hwmElemID = SideToElem(S2E_NB_ELEM_ID, SideID)
    END IF

    ! Search for h_wm element
    CALL FindHwmElementMPI(INT(dataToRecv(3)), INT(dataToRecv(4)), dataToRecv(5:7), dataToRecv(8:10), SideID, hwmElemID)

    IF (((hwmElemID - nElems) .GT. 0) .AND. ((hwmElemID - nElems) .LE. nNbProcs)) THEN ! hwmElem is in another proc
        sendCnt = sendCnt + 1
        ! Mount data to be sent
        ! dataToSend(1,...), i.e., the rank of the sender, is not changed, because we want
        ! the mpi proc with the element to communicate his ownership directly to the 
        ! mpi proc responsible for boundary imposition
        dataToSend(1,sendCnt) = dataToRecv(1)
        dataToSend(2,sendCnt) = SideToGlobalSide(SideID) ! Global SideID of "lower side" at adj. element, so that the receiving MPI proc.
                                                        ! can get his elemID, in order to get the oppositeSide and thus its info (face_xgp etc)
                                                        ! Unfortunately, this info is all tha THIS proc. have.
        dataToSend(3:4,sendCnt) = dataToRecv(3:4)
        dataToSend(5:7,sendCnt) = dataToRecv(5:7) ! Position vector of p,q at wall side
        dataToSend(8:10,sendCnt) = dataToRecv(8:10) ! Normal vector of p,q at wall side

        CALL MPI_ISEND(dataToSend(:,sendCnt),10, MPI_DOUBLE_PRECISION, NbProc(hwmElemID-nElems), &
                0, MPI_COMM_FLEXI, commRequests(hwmElemID-nElems), iErr)
        commCntSend(hwmElemID-nElems) = commCntSend(hwmElemID-nElems)+1
    ELSE ! Element is within my partition

        ! Mark this point's calculation (p,q,WSideID) as "my responsibility"
        ! Also, store the NormVec and Face_xGP positions so that we can, once again,
        ! see where in the physical domain h_wm lies, to check if we can approximate it
        ! as any point inside this element, or if we need to interpolate.
        nTauW_MINE = nTauW_MINE+1
        tmpTauW_NormVec_MINE(:,nTauW_MINE) = dataToRecv(8:10)
        tmpTauW_FacexGP_MINE(:,nTauW_MINE) = dataToRecv(5:7)
        tmpTauW_Element_MINE(nTauW_MINE) = hwmElemID

    END IF
END DO

! By now, every MPI proc finished sending any info needed about finding the h_wm element.
! Therefore, 

#endif /* USE_MPI */


!-----------------------------------------------------
! Transform temp variables into permanent ones
!-----------------------------------------------------
ALLOCATE(TauW_NormVec_MINE(3,nTauW_MINE))
ALLOCATE(TauW_FacexGP_MINE(3,nTauW_MINE))
ALLOCATE(TauW_Element_MINE(nTauW_MINE))

DO i=1,nTauW_MINE
    TauW_Element_MINE(i) = tmpTauW_Element_MINE(i)
    TauW_NormVec_MINE(:,i) = tmpTauW_NormVec_MINE(:,i)
    TauW_FacexGP_MINE(:,i) = tmpTauW_FacexGP_MINE(:,i)
END DO

DEALLOCATE(tmpTauW_Element_MINE)
DEALLOCATE(tmpTauW_NormVec_MINE)
DEALLOCATE(tmpTauW_FacexGP_MINE)


! From this point on, every partition that will calculate any of the WMLES_Tauw already knowns
! their responsibility.
! What is needed, now, is to check whether we can approximate the h_wm point as any of OUR points

ALLOCATE(TauW_CalcInfo_MINE(5,nTauW_MINE))
TauW_CalcInfo_MINE = 0

! Check if any points that are my responsibility can be approximated by a point inside my element (or faces of this element)
! Otherwise, there is a need for interpolation.
DO i=1,nTauW_MINE
    ! Calc x,y,z coordinates of h_wm within my element
    h_wm_Coords = -1.*(h_wm*delta)*TauW_NormVec_MINE(:,i) + TauW_FacexGP_MINE(:,i)

    WMLES_Tol = 0.15 ! temporary. Find a way to calculate this tolerance better (inscribed sphere? minimum distance between sides?)
    Found = .FALSE.

    ! Start searching on faces
    DO LocSide=1,6
        DO q=0,PP_NZ
            DO p=0,PP_N
                DistanceVect = h_wm_Coords - Face_xGP(:,p,q,0, ElemToSide(E2S_SIDE_ID, LocSide, TauW_Element_MINE(i)))
                Distance = 0
                DO j=1,3
                    Distance = Distance + DistanceVect(j)**2
                END DO
                Distance = SQRT(Distance)

                IF (Distance .LE. WMLES_Tol) THEN
                    ! Approximate this h_wm point to the Face_xGP point. No interpolation needed.
                    TauW_CalcInfo_MINE(1,i) = 1
                    TauW_CalcInfo_MINE(2:4,i) = (/p,q,0/)
                    TauW_CalcInfo_MINE(5,i) = ElemToSide(E2S_SIDE_ID, LocSide, TauW_Element_MINE(i))
                    Found = .TRUE.
                    EXIT
                END IF
            END DO
            IF (Found) EXIT
        END DO
        IF (Found) EXIT
    END DO

    ! Nothing found on faces. Check inside elements
    DO r=0,PP_NZ
        IF (Found) EXIT
        DO q=0,PP_N
            IF (Found) EXIT
            DO p=0,PP_N
                IF (Found) EXIT

                DistanceVect = h_wm_Coords - Elem_xGP(:,p,q,r,TauW_Element_MINE(i))
                Distance = 0
                DO j=1,3
                    Distance = Distance + DistanceVect(j)**2
                END DO
                Distance = SQRT(Distance)

                IF (Distance .LE. WMLES_Tol) THEN
                    ! Approximate this h_wm point to the Elem_xGP point. No interpolation needed.
                    TauW_CalcInfo_MINE(1,i) = 2
                    TauW_CalcInfo_MINE(2:4,i) = (/p,q,r/)
                    TauW_CalcInfo_MINE(5,i) = TauW_Element_MINE(i)
                    Found = .TRUE.
                END IF
            END DO
        END DO
    END DO

    ! If Found is still .FALSE., there is a need for interpolation.
    ! Hence, set up the interpolation vars for this point.

    ! Steps:
    ! 1. Calculate Lagrange polynomials, in each direction, for the point in question
    ! (Should this be on standard or physical domain? .... Check it out)
    ! 2. Set up this matrix of l_xi,l_eta,l_zeta.
    ! 3. Map somehow this new matrix to the current point (possibily in an index with TauW_MINE (i) index)
END DO


CALL Abort(__STAMP__, ".")

IF (nMasterSides+nSlaveSides .NE. nWMLESSides) CALL Abort(__STAMP__,&
                                        'Number of WMLES sides in mappings do not match!')

ALLOCATE(SlaveToTSide(1:nSlaveSides))
ALLOCATE(MasterToTSide(1:nMasterSides))
ALLOCATE(SlaveToWMLESSide(1:nSlaveSides))
ALLOCATE(MasterToWMLESSide(1:nMasterSides))

DO iSide=1,nMasterSides
    IF ((tmpMasterToTSide(iSide) .EQ. 0) .OR. (tmpMasterToWMLES(iSide) .EQ. 0)) THEN
        CALL Abort(__STAMP__, "Hooo ho ho Master")
    END IF
  MasterToTSide(iSide) = tmpMasterToTSide(iSide)
  MasterToWMLESSide(iSide) = tmpMasterToWMLES(iSide)
END DO
DO iSide=1,nSlaveSides
    IF ((tmpSlaveToTSide(iSide) .EQ. 0) .OR. (tmpSlaveToWMLES(iSide) .EQ. 0)) THEN
        CALL Abort(__STAMP__, "Hooo ho ho Slave")
    END IF
  SlaveToTSide(iSide) = tmpSlaveToTSide(iSide)
  SlaveToWMLESSide(iSide) = tmpSlaveToWMLES(iSide)
END DO
DEALLOCATE(tmpMasterToTSide)
DEALLOCATE(tmpSlaveToTSide)
DEALLOCATE(tmpMasterToWMLES)
DEALLOCATE(tmpSlaveToWMLES)

ALLOCATE(WMLES_Tauw(1:2, 0:PP_N, 0:PP_NZ, 1:nWMLESSides))
WMLES_Tauw = 0.

! Inverse mapping, from WMLES_Side ID -> SideID of BC
! Note that neither map has any particular ordering of elements
ALLOCATE(WMLES_SideInv(1:nWMLESSides))
WMLES_SideInv = 0
DO SideID=1,nBCSides
    IF (WMLES_Side(SideID) .GT. 0) THEN
        WMLES_SideInv(WMLES_Side(SideID)) = SideID
    END IF
END DO

! Warn the user if a BoundaryType has not been set in parameter file according to WMLES
#if USE_MPI
CALL MPI_REDUCE(nWMLESSides, TotalNWMLESSides, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_FLEXI, iError)
#else
TotalNWMLESSides = nWMLESSides
#endif
IF (MPIRoot .AND. (TotalNWMLESSides .EQ. 0)) THEN
    CALL PrintWarning('A Wall Model is set up, but no BoundaryConditions'// &
                                ' of WMLES Type has been found on the Mesh/Paramter File!')
END IF 

! Debugging / Logging info

LOGWRITE(UNIT_logOut,*) '========= WMLES INFO =========='
LOGWRITE(*,*) 'nWMLESSides', nWMLESSides
IF (nWMLESSides .NE. 0) THEN
  LOGWRITE(*,*) 'nMasterSides', nMasterSides
  LOGWRITE(*,*) 'nSlaveSides', nSlaveSides
  DO iSide=1,nMasterSides
    LOGWRITE(*,*) '=========== MasterSide', iSide, '============='
    LOGWRITE(*,*) 'SideID', MasterToTSide(iSide)
    LOGWRITE(*,*) 'Global SideID', SideToGlobalSide(MasterToTSide(iSide))
    LOGWRITE(*,*) 'Wall SideID', WMLES_SideInv(MasterToWMLESSide(iSide))
    LOGWRITE(*,*) 'Wall WMLESSideID', MasterToWMLESSide(iSide)
    LOGWRITE(*,*) 'Wall ElemID', SideToElem(S2E_NB_ELEM_ID, MasterToTSide(iSide))
    LOGWRITE(*,*) 'Adj. ElemID', SideToElem(S2E_ELEM_ID, MasterToTSide(iSide))
  END DO
  DO iSide=1,nSlaveSides
    LOGWRITE(*,*) '=========== SlaveSide', iSide, '============='
    LOGWRITE(*,*) 'SideID', SlaveToTSide(iSide)
    LOGWRITE(*,*) 'Global SideID', SideToGlobalSide(SlaveToTSide(iSide))
    LOGWRITE(*,*) 'Wall SideID', WMLES_SideInv(SlaveToWMLESSide(iSide))
    LOGWRITE(*,*) 'Wall WMLESSideID', SlaveToWMLESSide(iSide)
    LOGWRITE(*,*) 'Wall ElemID', SideToElem(S2E_ELEM_ID, SlaveToTSide(iSide))
    LOGWRITE(*,*) 'Adj. ElemID', SideToElem(S2E_NB_ELEM_ID, SlaveToTSide(iSide))
  END DO
END IF
LOGWRITE(*,*) '========= END OF WMLES INFO =========='
IF (Logging) FLUSH(UNIT_logOut)


! ALLOCATE(nbElemID(nWMLESSides))
! nLocalNbElem=nWMLESSides
! ! Count how many slave and master top sides there are (for looping purposes)
! ! and populate the mapping between WMLESSideID to OppSideID
! DO iWSide=1,nWMLESSides
!   WSideID = WMLES_SideInv(iWSide)
!   WElemID = SideToElem(S2E_ELEM_ID, WSideID)
!   TSideID = ElemToSide(E2S_SIDE_ID, ETA_PLUS, WElemID)
!   OppSideFlip = ElemToSide(E2S_FLIP, ETA_PLUS, WElemID)

!   ! Get ID of neighbor element of the wall adjacent element
!   ! If TSideFlip == 0, WElem is the root of the face, so S2E_NB_ELEM_ID gives the adjacent element
!   IF (TSideFlip .EQ. 0) THEN ! WElem is "master" on this side
!     nbElemID(iWSide) = SideToElem(S2E_NB_ELEM_ID, TSideID)
!   ELSE ! WElem is "slave" on this side
!     nbElemID(iWSide) = SideToElem(S2E_ELEM_ID, TSideID)
!   END IF

! #if USE_MPI
!   ! Check if neighbor (adjacent) element is in another MPI processor (i.e. nbElemID = -1)
!   ! If so, find out who is the master/slave of the side so the info is properly gathered on U_master/slave arrays
!   IF (nbElemID(iWSide) .LT. 0) THEN
!     nLocalNbElem = nLocalNbElem - 1 ! Decrease local (relative to processor) neighbor element count
!     IF (TSideID .LT. firstMPISide_YOUR) WRITE(*,*) 'Side is MINE'
!     IF (TSideID .GE. firstMPISide_YOUR) WRITE(*,*) 'Side is YOUR'
!     DO iNbProc=1,nNbProcs
!         IF ((TSideID .GT. offsetMPISides_MINE(iNbProc-1)) .AND. (TSideID .LE. offsetMPISides_MINE(iNbProc))) THEN ! Element adjacent to this side is from NbProc(iNbProc)
!           WRITE(*,*) 'myRank, TopSideID, nSides, offseti-1, offset; is MINE Side, Adjacent element is from'
!           WRITE(*,*) myRank, TSideID, nSides, offsetMPISides_MINE(iNbProc-1), offsetMPISides_MINE(iNbProc), NbProc(iNbProc)
!         ELSE IF ((TSideID .GT. offsetMPISides_YOUR(iNbProc-1)) .AND. (TSideID .LE. offsetMPISides_YOUR(iNbProc))) THEN ! Element adjacent to this side is from NbProc(iNbProc)
!           WRITE(*,*) 'myRank, TopSideID, nSides, offseti-1, offset; is YOUR Side, Adjacent element is from'
!           WRITE(*,*) myRank, TSideID, nSides, offsetMPISides_MINE(iNbProc-1), offsetMPISides_MINE(iNbProc), NbProc(iNbProc)
!         END IF
!     END DO
!   END IF
! #endif /* USE_MPI */
!END DO

WMLESInitDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT Wall-Modeled LES DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitWMLES


!==================================================================================================================================
!> Get the ID of the side opposite of that indicated by SideID, in the same element
!==================================================================================================================================
SUBROUTINE GetOppositeSide(SideID, ElemID)
! MODULES
USE MOD_Mesh_Vars,              ONLY: SideToElem, ElemToSide
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, INTENT(IN)                 :: ElemID
INTEGER, INTENT(INOUT)              :: SideID
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
    SideID = ElemToSide(E2S_SIDE_ID, ETA_PLUS, ElemID)
CASE (ETA_PLUS)
    SideID = ElemToSide(E2S_SIDE_ID, ETA_MINUS, ElemID)
CASE (XI_MINUS)
    SideID = ElemToSide(E2S_SIDE_ID, XI_PLUS, ElemID)
CASE (XI_PLUS)
    SideID = ElemToSide(E2S_SIDE_ID, XI_MINUS, ElemID)
CASE (ZETA_MINUS)
    SideID = ElemToSide(E2S_SIDE_ID, ZETA_PLUS, ElemID)
CASE (ZETA_PLUS)
    SideID = ElemToSide(E2S_SIDE_ID, ZETA_MINUS, ElemID)
END SELECT

END SUBROUTINE

!==================================================================================================================================
!> Get the ID of element where the exchange location point, h_wm (p,q), lies.
!==================================================================================================================================
SUBROUTINE FindHwmElement(p, q, SideID, hwmElemID)
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,              ONLY: NormVec, Face_xGP, SideToElem, ElemToSide, nElems
USE MOD_WMLES_Vars
#if USE_MPI
USE MOD_MPI_Vars
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, INTENT(IN)                 :: p, q
INTEGER, INTENT(INOUT)              :: hwmElemID, SideID!, OppSideID
!----------------------------------------------------------------------------------------------------------------------------------
REAL                                :: Height, abs_h_wm, Tol
REAL                                :: DistanceVect(1:3), InwardNorm(1:3)
LOGICAL                             :: FoundOrNomB, Master
INTEGER                             :: tmpElem, WallSideID, cnt
INTEGER                             :: iNbProc
!==================================================================================================================================
FoundOrNomB = .FALSE. ! Logical to check whether h_wm has been Found or if it is None of my Bussiness (i.e. not in my partition)

abs_h_wm = h_wm*delta ! get absolute value of off-wall distance
WallSideID = SideID

! Get wall adjacent element from SideID
hwmElemID = SideToElem(S2E_ELEM_ID, WallSideID) ! ELEM_ID: Wall elements are always masters of BC sides

! Calc inward normal vector with magnitude h_wm, at point p,q. Use it to check in which
! element does h_wm lies.
InwardNorm = -1.*abs_h_wm*NormVec(1:3,p,q,0,WallSideID)

! This element check is naive: for highly curved walls with low polynomial degrees, there
! might exist a h_wm with y-component < d_y, but lieing within the next element, not
! the current one. This might happen for either concave or convex walls
cnt = 0
DO WHILE (FoundOrNomB .EQV. .FALSE.)
    cnt = cnt + 1
    CALL GetOppositeSide(SideID, hwmElemID)
    DistanceVect = Face_xGP(1:3,p,q,0,SideID) - Face_xGP(1:3,p,q,0,WallSideID)
    tmpElem = hwmElemID
    Tol = 1E-3 ! Tol should be set to 10% of distance between nodes.
    ! In order to do that, calc the element "height" (mag(distanceVect)/nElems)/PP_N
    ! where nElems is the number of elements we've advanced in this loop.

    IF (ABS(DistanceVect(2) - InwardNorm(2)) .GT. Tol) THEN ! IF TRUE: NOT THE SAME POINT, so check which element
        IF ((DistanceVect(2) - InwardNorm(2)) .LT. 0) THEN ! Negative distance: h_wm "higher" than d: next elem
            hwmElemID = SideToElem(S2E_ELEM_ID, SideID) ! Get adjacent element interfacing at SideID
            ! If hwmElemID is not equal to tmpElem (or -1), it means that THIS element is slave, not master, of OppSide
            ! since ELEM_ID was used.
            Master = .FALSE.
            IF (hwmElemID .EQ. tmpElem) THEN ! Element is master of OppSide: get the nb element then
                Master = .TRUE.
                hwmElemID = SideToElem(S2E_NB_ELEM_ID, SideID)
            END IF
#if USE_MPI
            IF (hwmElemID .EQ. -1) THEN ! Element in another MPI proc.
                ! Thus, it is none of my business -- get out of the loop
                FoundOrNomB = .TRUE.

                ! Identify which neighboring processor owns the element.
                IF (Master .EQV. .TRUE.) THEN
                    DO iNbProc=1,nNbProcs
                        IF ((SideID .GT. OffsetMPISides_MINE(iNbProc-1)) .AND. (SideID .LE. OffsetMPISides_MINE(iNbProc))) THEN
                            hwmElemID = nElems + iNbProc ! The information on which nbProc holds the next element goes in hwmElemID
                            EXIT
                        END IF
                    END DO
                ELSE
                    DO iNbProc=1,nNbProcs
                        IF ((SideID .GT. OffsetMPISides_YOUR(iNbProc-1)) .AND. (SideID .LE. OffsetMPISides_YOUR(iNbProc))) THEN
                            hwmElemID = nElems + iNbProc ! The information on which nbProc holds the next element goes in hwmElemID
                            EXIT
                        END IF
                    END DO
                END IF
            END IF
#endif
        ELSE ! Positive: h_wm "smaller" than d. Consider being within this element.
            FoundOrNomB = .TRUE.
        END IF
    ELSE ! h_wm_y and d_y are very close: h_wm lies within this element.
        ! TODO: Since it lies on the face, say that it is within the upper element, so we get
        ! info from the upper polynomial approximation, which is generally better (mesh reqs.)
        FoundOrNomB = .TRUE.
    END IF
END DO

END SUBROUTINE FindHwmElement

#if USE_MPI
!==================================================================================================================================
!> Get the ID of element where the exchange location point, h_wm (p,q), lies.
!==================================================================================================================================
SUBROUTINE FindHwmElementMPI(p, q, NormVecWall, FacexGPWall, SideID, hwmElemID)
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,              ONLY: NormVec, Face_xGP, SideToElem, ElemToSide, nElems
USE MOD_WMLES_Vars
#if USE_MPI
USE MOD_MPI_Vars
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, INTENT(IN)                 :: p, q
REAL, INTENT(IN)                    :: NormVecWall(3), FacexGPWall(3)
INTEGER, INTENT(INOUT)              :: hwmElemID, SideID!, OppSideID
!----------------------------------------------------------------------------------------------------------------------------------
REAL                                :: Height, abs_h_wm, Tol
REAL                                :: DistanceVect(1:3), InwardNorm(1:3)
LOGICAL                             :: FoundOrNomB, Master
INTEGER                             :: tmpElem, WallSideID, cnt
INTEGER                             :: iNbProc
!==================================================================================================================================
FoundOrNomB = .FALSE.
InwardNorm = -1.*(h_wm*delta)*NormVecWall(1:3)

DO WHILE (FoundOrNomB .EQV. .FALSE.)
    CALL GetOppositeSide(SideID, hwmElemID)
    DistanceVect = Face_xGP(:,p,q,0,SideID) - FacexGPWall(:)
    tmpElem = hwmElemID
    Tol = 1e-3
    IF (ABS(DistanceVect(2) - InwardNorm(2)) .GT. Tol) THEN ! It can't be approximate as the face point
        IF ((DistanceVect(2) - InwardNorm(2)) .LT. 0) THEN ! h_wm above this element
            hwmElemID = SideToElem(S2E_ELEM_ID, SideID)
            Master = .FALSE.
            IF (hwmElemID .EQ. tmpElem) THEN
                Master = .TRUE.
                hwmElemID = SideToElem(S2E_NB_ELEM_ID, SideID)
            END IF

            IF (hwmElemID .EQ. -1) THEN ! Next element in another MPI proc.
                FoundOrNomB = .TRUE. ! This is none of my business then.

                ! Identify which neighboring processor owns the element.
                IF (Master .EQV. .TRUE.) THEN
                    DO iNbProc=1,nNbProcs
                        IF ((SideID .GT. OffsetMPISides_MINE(iNbProc-1)) .AND. (SideID .LE. OffsetMPISides_MINE(iNbProc))) THEN
                            hwmElemID = nElems + iNbProc ! The information on which nbProc holds the next element goes in hwmElemID
                            EXIT
                        END IF
                    END DO
                ELSE
                    DO iNbProc=1,nNbProcs
                        IF ((SideID .GT. OffsetMPISides_YOUR(iNbProc-1)) .AND. (SideID .LE. OffsetMPISides_YOUR(iNbProc))) THEN
                            hwmElemID = nElems + iNbProc ! The information on which nbProc holds the next element goes in hwmElemID
                            EXIT
                        END IF
                    END DO
                END IF
            END IF

        ELSE ! h_wm is within this element
            FoundOrNomB = .TRUE.
        END IF

    ELSE ! This h_wm can be approximated as the face point
        ! TODO (in both functions)
        ! Mark this as the element above, so the info on the face of the next element is used
        ! In order to do so, a check whether it is in another MPI proc or MINE should be done
        FoundOrNomB = .TRUE.
    END IF
END DO

END SUBROUTINE FindHwmElementMPI
#endif /* USE_MPI */

!==================================================================================================================================
!> Compute the wall stress, tau_w, at each point in a WMLES BC surface.
!==================================================================================================================================
SUBROUTINE ComputeWallStress()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_WMLES_Vars
USE MOD_Mesh_Vars
USE MOD_DG_Vars                     ,ONLY:UPrim_master, UPrim_slave
USE MOD_EOS_Vars                    ,ONLY: mu0
IMPLICIT NONE

!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                             :: iSide
INTEGER                             :: p,q
REAL                             :: u_tau, u_mean, offdist
REAL, PARAMETER                  :: vKarman=0.41, B=5.5
!==================================================================================================================================
IF (.NOT.WMLESInitDone) THEN
    CALL CollectiveStop(__STAMP__,&
    'Cannot compute Wall Stress -- Wall-Modeled LES INIT is not done.')
END IF
IF (nWMLESSides .EQ. 0) RETURN ! No WMLES Boundary Sides defined on this proc. No business here.

SELECT CASE(WallModel)

  CASE(WMLES_SCHUMANN)
    
    ! Schumann model
    ! Computation is done using instantaneous information, retrieved from the bottom boundary
    ! of the top neighbor element of the element adjacent to the wall, as in Frere et al.
    
    ! The channel flow testcase assumes <tau_w> = 1.0

    ! Slave top sides first
    LOGWRITE(*,*) '================= Wall Stress Calc ==================='
    DO iSide=1,nSlaveSides
      LOGWRITE(*,*) '------- Slave Side', iSide, '--------'
      LOGWRITE(*,*) 'SideID', SlaveToTSide(iSide)
      LOGWRITE(*,*) 'Wall SideID', WMLES_SideInv(SlaveToWMLESSide(iSide))
      LOGWRITE(*,*) 'WMLES SideID', SlaveToWMLESSide(iSide)
      LOGWRITE(*,'(2(A4,2X),9(A15,2X))') 'p', 'q', 'offdist', 'Prim(1)', 'u_tau', 'mu0', 'u_mean', 'Prim(2)', 'Prim(3)', 'tau_xy', 'tau_yz'
      DO q=0,PP_NZ; DO p=0,PP_N
        ! This is constant, we can calculate it in pre-processing stage (one more vector, 1:nWMLESSides)
        offdist = ABS(Face_xGP(2,p,q,0,SlaveToTSide(iSide)) - Face_xGP(2,p,q,0, WMLES_SideInv(SlaveToWMLESSide(iSide)) ))
        u_tau = SQRT(1.0/UPrim_slave(1,p,q, SlaveToTSide(iSide)) )
        u_mean = u_tau*( (1./vKarman) * LOG((offdist*u_tau)/mu0) + B )
        WMLES_Tauw(1,p,q, SlaveToWMLESSide(iSide) ) = (UPrim_slave(2,p,q, SlaveToTSide(iSide))/u_mean) * 1.0
        WMLES_Tauw(2,p,q, SlaveToWMLESSide(iSide) ) = (2.0*mu0)*(UPrim_slave(3,p,q, SlaveToTSide(iSide))/offdist)

        LOGWRITE(*,'(2(I4,2X),9(E15.8,2X))') p, q, offdist, UPrim_slave(1,p,q, SlaveToTSide(iSide)), u_tau, mu0, u_mean, UPrim_slave(2,p,q, SlaveToTSide(iSide)), UPrim_slave(3,p,q, SlaveToTSide(iSide)), WMLES_Tauw(1,p,q,SlaveToWMLESSide(iSide)), WMLES_Tauw(2,p,q,SlaveToWMLESSide(iSide))

      END DO; END DO ! p,q
    END DO
    ! Master top sides next
    DO iSide=1,nMasterSides
      LOGWRITE(*,*) '------- Master Side', iSide, '--------'
      LOGWRITE(*,*) 'SideID', MasterToTSide(iSide)
      LOGWRITE(*,*) 'Wall SideID', WMLES_SideInv(MasterToWMLESSide(iSide))
      LOGWRITE(*,*) 'WMLES SideID', MasterToWMLESSide(iSide)
      LOGWRITE(*,'(2(A4,2X),9(A15,2X))') 'p', 'q', 'offdist', 'Prim(1)', 'u_tau', 'mu0', 'u_mean', 'Prim(2)', 'Prim(3)', 'tau_xy', 'tau_yz'
      DO q=0,PP_NZ; DO p=0,PP_N
        offdist = ABS(Face_xGP(2,p,q,0,MasterToTSide(iSide)) - Face_xGP(2,p,q,0, WMLES_SideInv(MasterToWMLESSide(iSide)) ))
        u_tau = SQRT(1.0/UPrim_master(1,p,q, MasterToTSide(iSide)) )
        u_mean = u_tau*( (1./vKarman) * LOG((offdist*u_tau)/mu0) + B )
        WMLES_Tauw(1,p,q, MasterToWMLESSide(iSide) ) = (UPrim_master(2,p,q, MasterToTSide(iSide))/u_mean) * 1.0
        WMLES_Tauw(2,p,q, MasterToWMLESSide(iSide) ) = (2.0*mu0)*(UPrim_master(3,p,q, MasterToTSide(iSide))/offdist)

        LOGWRITE(*,'(2(I4,2X),9(E15.8,2X))') p, q, offdist, UPrim_master(1,p,q, MasterToTSide(iSide)), u_tau, mu0, u_mean, UPrim_master(2,p,q, MasterToTSide(iSide)), UPrim_master(3,p,q, MasterToTSide(iSide)), WMLES_Tauw(1,p,q,MasterToWMLESSide(iSide)), WMLES_Tauw(2,p,q,MasterToWMLESSide(iSide))

      END DO; END DO ! p,q
    END DO
    LOGWRITE(*,*) '================= End of Wall Stress Calc ==================='
    IF (Logging) FLUSH(UNIT_logOut)

  CASE DEFAULT
    CALL abort(__STAMP__,&
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
SDEALLOCATE(WMLES_Side)
SDEALLOCATE(WMLES_SideInv)

END SUBROUTINE FinalizeWMLES

!==================================================================================================================================
END MODULE MOD_WMLES