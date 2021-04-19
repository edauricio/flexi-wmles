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

MODULE MOD_WMLES_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE


INTEGER, PARAMETER :: WMLES_SCHUMANN = 1
INTEGER, PARAMETER :: WMLES_LOGLAW = 2
INTEGER, PARAMETER :: WMLES_WERNERWANGLE = 3
INTEGER, PARAMETER :: WMLES_REICHARDT = 4
INTEGER, PARAMETER :: WMLES_SPALDING = 5
INTEGER, PARAMETER :: WMLES_EQTBLE = 6
INTEGER, PARAMETER :: WMLES_COUETTE = 7

!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER                     :: WallModel        ! Integer corresponding to the WallModel
REAL                        :: delta            ! Approximate, or characteristic boundary layer thickness
REAL                        :: vKarman          ! von Karman constant
REAL                        :: B                ! Intercept coefficient for log-law-based models
REAL,ALLOCATABLE            :: h_wm(:,:,:)
INTEGER,ALLOCATABLE         :: WMLESFlip(:)
REAL,ALLOCATABLE            :: WMLES_Tauw(:,:,:,:) ! Wall stress tensor.
                                                   ! First index: 1 or 2, where 1 is tau_xy and 2 is tau_yz
                                                   ! Second and third indices: indices "i,j" of the BC face
                                                   ! Fourth index: WMLES Side

INTEGER                     :: nWMLESSides, nMasterWMLESSide, nSlaveWMLESSide ! Number of WMLES BC Sides                                                   
INTEGER,ALLOCATABLE         :: BCSideToWMLES(:) ! Mapping between WMLES Side and Mesh BCSide.
                                           ! Usage: BCSideToWMLES(SideID), SideID \in [1:nBCSides]
                                           ! OUTPUT: [1:nWMLESSides]
INTEGER,ALLOCATABLE         :: WMLESToBCSide(:) ! Inverse of BCSideToWMLES mapping, that is,
                                                ! get SideID of BC from WMLESSideID
INTEGER,ALLOCATABLE         :: MasterToWMLESSide(:), SlaveToWMLESSide(:) ! Mapping between master/slave opposite sides and the boundary
                                                                         ! sides in terms of nWMLESSide number.
INTEGER,ALLOCATABLE         :: MasterToOppSide(:), SlaveToOppSide(:) ! Mapping between master/slave sides, in terms of nMaster/SlaveWMLESSides,
                                                                     ! and the OppSide in terms of "SideID"
LOGICAL                     :: WMLESInitDone = .FALSE.
LOGICAL                     :: UseSemiLocal
INTEGER                     :: WMLES_Filter
REAL, ALLOCATABLE           :: WMLES_FilterMat(:,:)

END MODULE MOD_WMLES_Vars