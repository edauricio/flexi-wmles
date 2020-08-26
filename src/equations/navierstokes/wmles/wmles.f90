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
INTEGER                         :: SideID, iSide
INTEGER                         :: WElemID, TSideID, TSideFlip
INTEGER                         :: TotalNWMLESSides
! temp
INTEGER, ALLOCATABLE            :: tmpMasterToTSide(:), tmpMasterToWMLES(:), tmpSlaveToTSide(:), tmpSlaveToWMLES(:)
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

! Count how many WMLES BC Sides we have defined,
! and populate the mappings:
! - Mesh SideID to WMLES SideID (For use in GetBoundaryFlux.f90 in conjunction with WMLES_Tauw)
! - WMLES SideID to TopSideID master/slave (side where info is gathered) (For use in ComputeWallStress)

! REMEMBER: It does not matter wheter the information (i.e. adjacent element) is on this or
! another MPI process. By the time we call ComputeWallStress, all the face information needed
! has already been communicated (step 3 of DG_WeakTimeDerivative).

! However, IF a mapping between side -> neighbor process rank is needed for further communication
! (e.g. when solving the TBL equations on the same mesh with p-refinement), it may be
! accomplished by checking SideToElem(S2E_(NB_)ELEM_ID, TSideID) == -1 (i.e., the needed elem is not on this proc)

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

DO SideID=1,nBCSides
    IF (BoundaryType(BC(SideID),BC_TYPE) .EQ. 5) THEN ! WMLES side
        nWMLESSides = nWMLESSides + 1
        WMLES_Side(SideID) = nWMLESSides

        WElemID = SideToElem(S2E_ELEM_ID, SideID)
        IF (ElemToSide(E2S_SIDE_ID, ETA_MINUS, WElemID) .EQ. SideID) THEN ! lower wall
          TSideID = ElemToSide(E2S_SIDE_ID, ETA_PLUS, WElemID)
          TSideFlip = ElemToSide(E2S_FLIP, ETA_PLUS, WElemID)
        ELSE IF (ElemToSide(E2S_SIDE_ID, ETA_PLUS, WElemID) .EQ. SideID) THEN ! upper wall
          TSideID = ElemToSide(E2S_SIDE_ID, ETA_MINUS, WElemID)
          TSideFlip = ElemToSide(E2S_FLIP, ETA_MINUS, WElemID)
        ELSE
          CALL Abort(__STAMP__, 'Vertical walls?')
        END IF

        IF (TSideFlip .EQ. 0) THEN ! WElem is master of Top side. Hence, we should use info from UPrim_slave (adjacent element)
          nSlaveSides = nSlaveSides + 1
          tmpSlaveToTSide(nSlaveSides) = TSideID
          tmpSlaveToWMLES(nSlaveSides) = nWMLESSides
        ELSE
          nMasterSides = nMasterSides + 1
          tmpMasterToTSide(nMasterSides) = TSideID
          tmpMasterToWMLES(nMasterSides) = nWMLESSides
        END IF
    END IF
END DO
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
! ! and populate the mapping between WMLESSideID to TSideID
! DO iWSide=1,nWMLESSides
!   WSideID = WMLES_SideInv(iWSide)
!   WElemID = SideToElem(S2E_ELEM_ID, WSideID)
!   TSideID = ElemToSide(E2S_SIDE_ID, ETA_PLUS, WElemID)
!   TSideFlip = ElemToSide(E2S_FLIP, ETA_PLUS, WElemID)

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
!USE MOD_WMLES_Vars
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

END DO; END DO ! p,q

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