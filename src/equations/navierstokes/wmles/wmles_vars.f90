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
INTEGER, PARAMETER :: WMLES_WERNERWANGLE = 2
INTEGER, PARAMETER :: WMLES_EQTBLE = 3

!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER                     :: WallModel ! Integer corresponding to the WallModel
REAL                        :: h_wm
REAL,ALLOCATABLE            :: WMLES_Tauw(:,:,:,:) ! Wall stress tensor.
                                                   ! First index: 1 or 2, where 1 is tau_xy and 2 is tau_yz
                                                   ! Second and third indices: indices "i,j" of the BC face
                                                   ! Fourth index: SideID
REAL,ALLOCATABLE            :: Side_OffWallDist(:,:,:) ! Off-wall distance h_wm for each point of each WMLES Side
                                                       ! First and second indices: "i,j" of face
                                                       ! Third index: WMLESSideID
INTEGER                     :: nWMLESSides ! Number of WMLES BC Sides                                                   
INTEGER,ALLOCATABLE         :: WMLES_Side(:) ! Mapping between WMLES BC Side and Mesh BCSide.
                                           ! Usage: WMLES_Side(SideID), SideID \in [1:nBCSides]
                                           ! OUTPUT: [1:nWMLESSides]
INTEGER,ALLOCATABLE         :: WMLES_SideInv(:) ! Inverse of WMLES_Side mapping, that is,
                                                ! get SideID of BC from WMLESSideID
INTEGER,ALLOCATABLE         :: SlaveToTSide(:), MasterToTSide(:) ! Mapping between WMLES wall side to the top side of element
                                                    ! arranged in master/slave so that we known which to use for the info on the adjacent element,
                                                    ! that is, UPrim_master/slave.
INTEGER,ALLOCATABLE         :: SlaveToWMLESSide(:), MasterToWMLESSide(:)
INTEGER                     :: nSlaveSides, nMasterSides ! number of slave and master (to this MPI proc) sides where info is to be exchanged
#if USE_MPI
INTEGER                     :: nLocalNbElem ! number of local neighbor elements
INTEGER                     :: nMPINbElem ! number of neighbor elements in other MPI partitions
#endif                                                    
LOGICAL                     :: WMLESInitDone = .FALSE.

!=================================================================================================================================

END MODULE MOD_WMLES_Vars