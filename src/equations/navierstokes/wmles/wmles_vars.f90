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
INTEGER, PARAMETER :: WMLES_FALKNER_SKAN = 8
INTEGER, PARAMETER :: WMLES_LAMINAR_INTEGRAL = 9
INTEGER, PARAMETER :: WMLES_FSLL = 10 ! Falker-Skan + Loglaw (complete regime model)

REAL, PARAMETER    :: PI = 3.141592653589

!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255)          :: WMConnectionFile ! Name of file containing the WM mapping/connection information
INTEGER                     :: WallModel        ! Integer corresponding to the WallModel
INTEGER                     :: offsetBCSides    ! offset for BC sides of this proc., related to the global BC side number
INTEGER                     :: nHWMInterpPoints, nHWMLocalPoints ! NuPmber of points where h_wm should be interpolated,
                                                ! and number of local h_wm points, i.e., where both h_wm and BC info is in
                                                ! the same MPI proc.
REAL                        :: vKarman          ! von Karman constant
REAL                        :: B                ! Intercept coefficient for log-law-based models
INTEGER                     :: NSuper           ! Parameter to find h_wm in standard coordinates, when interpolation is needed.
REAL,ALLOCATABLE            :: WMLES_Tauw(:,:,:,:) ! Wall stress tensor.
                                                   ! First index: 1 or 2, where 1 is tau_xy and 2 is tau_yz
                                                   ! Second and third indices: indices "i,j" of the BC face
                                                   ! Fourth index: WMLES Side
REAL,ALLOCATABLE            :: FSBeta(:,:,:)        ! Beta of Falkner-Skan equations for each (p,q,nWMLESSide) of the laminar region
!REAL,ALLOCATABLE            :: FSDelta(:,:,:)       ! Eta value (similarity variable) corresponding to the boundary layer edge, according
                                                    ! to the Falkner-Skan solution, for each point (p,q,nWMLESSide)
!REAL,ALLOCATABLE            :: FSEtaInf(:,:,:)    ! "Free boundary" eta_inf for each point (p,q,nWMLESSide) according to Falkner-Skan solution
REAL,ALLOCATABLE            :: FSWallShear(:,:,:)   ! Solution of f'' at eta=0 ( i.e., "FS Wall Shear") from Falkner Skan equation
                                                    ! for each point (p,q,nWMLESSide)
!REAL,ALLOCATABLE            :: FSEta(:,:,:,:)         ! Cached eta's for each step when advancing the Falkner-Skan solution in eta-direcion
                                                      ! using the RK scheme, for each point (p,q,nWMLESSide)

!REAL,ALLOCATABLE            :: FSFPrime(:,:,:,:)    ! Cached solution (f') of the Falkner-Skan equation at each "time-step" (steps along the
                                                    ! similarity variable direction), for each point (p,q,nWMLESSide)

INTEGER                     :: nWMLESSides ! Number of WMLES BC Sides                                                   
INTEGER                     :: nWMLaminarSides ! Number of Laminar WM BC Sides                                                   
INTEGER                     :: nWMTurbulentSides ! Number of Turbulent WM BC Sides                                                   
INTEGER,ALLOCATABLE         :: BCSideToWMLES(:) ! Mapping between WMLES BC Side and Mesh BCSide.
                                           ! Usage: BCSideToWMLES(SideID), SideID \in [1:nBCSides]
                                           ! OUTPUT: [1:nWMLESSides]
INTEGER,ALLOCATABLE         :: WMLESToBCSide(:) ! Inverse of BCSideToWMLES mapping, that is,
                                                ! get SideID of BC from WMLESSideID
INTEGER,ALLOCATABLE         :: WMLESToLaminarSide(:)      ! Mapping from WMLES Side ID to Laminar Side number
INTEGER,ALLOCATABLE         :: WMLESToTurbulentSide(:)    ! Mapping from WMLES Side ID to Turbulent Side number
INTEGER,ALLOCATABLE         :: WMLESShapeInfo(:,:)  ! Information about the number of the shape (when two shapes are present, the first is laminar and the second turbulent)
                                                    ! Also, second position for first index indicates whether this side is a boundary between shapes

LOGICAL                     :: WMLESInitDone = .FALSE.

! Variables for the solution of the Falkner-Skan equation
REAL                        :: beta_0, beta            ! Parameters defining the actual equation to be solved
REAL                        :: eta_inf, alpha          ! Final solutions
REAL, ALLOCATABLE           :: xi_16(:)                ! Vector of time points when solving the RHS16 system
REAL, ALLOCATABLE           :: sol_16(:,:)             ! Solution vector of the RHS16 system for each time point xi_16
REAL                        :: fs_f, fs_u, fs_v        ! Solution variables of the first ODE system solved, used subsequently in the other two systems
                                                       ! These variables actually correspond to f, f' and f''

#if USE_MPI
!----------------------------------------------------------------------------------------------------------------------------------
! MPI VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER                     :: nWMLESRecvProcs, nWMLESSendProcs ! Number of sending/receiving procs for each MPI rank (local info)
INTEGER                     :: nHWMPropSend ! Number of flow properties in h_wm to be sent to the element responsible for WMLES Side
INTEGER, ALLOCATABLE        :: nHWMSendPoints(:), nHWMRecvPoints(:) ! Arrays containing the number of points to send/receive for each
                                                                    ! sending/receiving proc.
INTEGER, ALLOCATABLE        :: SendToInterpPoint(:), LocalToInterpPoint(:) ! Arrays containing the number mapping between the send
                                                                    ! or the local point to the interpolation point of h_wm, in the local
                                                                    ! proc.
INTEGER, ALLOCATABLE        :: WMLESSendProc(:), WMLESRecvProc(:)   ! Mapping between the send/receive list local to proc (i.e., 1:nWMLESSendProcs
                                                                    ! and 1:nWMLESRecvProc) and the global MPI rank of that proc
REAL, ALLOCATABLE           :: HWMSendInfo(:,:), HWMRecvInfo(:,:) ! Arrays containing the sending/receiving information. For the first one,
                                                                  ! info about the flow property is sent, and the second contains the 
                                                                  ! corresponding BC Side / p,q point at which this received flow property
                                                                  ! should be use in the modelling procedure
REAL, ALLOCATABLE           :: HWMInterpInfo(:,:) ! Array containing the informataion for interpolation of flow properties at h_wm
                                                  ! point, i.e., local element ID and the corresponding standard coordinates
REAL, ALLOCATABLE        :: HWMLocalInfo(:,:) ! Array containing the information about BC Side / p,q point of the correspongin
                                                 ! h_wm point that is interpolated locally.
REAL, ALLOCATABLE           :: HWMInfo(:,:,:,:) ! Array containing the necessary information from h_wm in order to feed the wall model
INTEGER, ALLOCATABLE        :: WMLESSendRange(:,:), WMLESRecvRange(:,:) ! Interval ranges (first index: left (1)/right (2)) for each
                                                  ! sending/receiving MPI proc. (local list, i.e. nWMLESSendProcs/nWMLESRecvProcs)
INTEGER, ALLOCATABLE        :: WMLES_RecvRequests(:), WMLES_SendRequests(:) ! Requests for the non-blocking send/receive operations











INTEGER, ALLOCATABLE        :: nTauW_MINE(:) ! Number of points to calculate tau_w, for each MPI proc. (nProcs_SendTauW)
                                             ! this partition needs to send the results later
INTEGER, ALLOCATABLE        :: nTauW_YOURS(:) !


REAL, ALLOCATABLE           :: TauW_MINE(:,:,:) ! Array to store wall shear stress computed by this MPI proc. that must be sent
                                                ! to the MPI proc. responsible for BC imposition.
                                                ! First index: 1 or 2, where 1 is tau_xy and 2 is tau_zy
                                                ! Second index: Local points of my responsibility, with respect to each MPI proc that we must communicate later
                                                ! Third index: MPI proc responsible for BC imposition -- the one we should send the results later

REAL, ALLOCATABLE           :: TauW_YOURS(:,:,:) ! Array to store wall shear stress received from MPI procs. responsible for tau_w calculation
                                                 ! First index: 1 or 2, where 1 is tau_xy and 2 is tau_zy
                                                 ! Second index: Local points, associated to a WMLES Side, and with respect to another MPI proc.
                                                 ! Third index: MPI proc to receive tau_w from.

REAL, ALLOCATABLE           :: TauW_MINE_NormVec(:,:,:) ! Wall-normal vectors (1:3,:,:) for each point that this MPI proc.
                                                        ! is responsible for computing tau_w (:,j,:), and each MPI proc. to send
                                                        ! the data (:,:,i)

INTEGER, ALLOCATABLE        :: TauW_Proc(:,:,:,:) ! Mapping between the node (p,q) of BC imposition within a WMLES Side
                                                  ! and the local point/MPI proc. responsible for calculating Tau_W for this node/side

INTEGER                     :: nProcs_RecvTauW, nProcs_SendTauW ! Number of MPI procs. to send/receive tau_w info

INTEGER, ALLOCATABLE        :: Proc_RecvTauW(:), Proc_SendTauW(:) ! Mapping between local to global index of MPI procs to send/receive TauW
                                                                  ! [1:nProcs_(Send/Recv)TauW]
INTEGER, ALLOCATABLE        :: Proc_RecvTauW_Inv(:) ! Inverse mapping of the above, that is, from global MPI rank to local receiving rank.

INTEGER, ALLOCATABLE        :: nTauW_MINE_FacePoint(:), nTauW_MINE_InteriorPoint(:), nTauW_MINE_Interpolate(:) ! Number of face/interior points that may
                                                                                                      ! be used to approximate h_wm, or number of
                                                                                                      ! h_wm points that must be interpolated, for each MPI proc.
                                                                                                      ! that this partition is responsible for calculating/comm.
                                                                                                      ! tau_w values

INTEGER, ALLOCATABLE        :: FaceToLocalPoint(:,:), InteriorToLocalPoint(:,:), InterpToLocalPoint(:,:) ! Mappings between face/interior/interpolation point and
                                                                                                   ! the local points for each MPI proc. that THIS proc. should 
                                                                                                   ! calculate tau_w for

INTEGER, ALLOCATABLE        :: TauW_MINE_FacePoint(:,:,:) ! If h_wm may be approximated as a face point, then this array stores
                                                              ! mapping from the local tau_w point with respect to each MPI proc. to its p,q and SideID
                                                              ! First index: 1-3 (/p,q,SideID/)
                                                              ! Second index: Local Tau_W calc point (with respect to each MPI proc. responsible for imposition)
                                                              ! Third index: MPI proc. to send info

INTEGER, ALLOCATABLE        :: TauW_MINE_InteriorPoint(:,:,:) ! If h_wm may be approximated as an interior element point, then this array stores
                                                              ! mapping from the local tau_w point with respect to each MPI proc. to its p,q,r and ElemID
                                                              ! First index: 1-4 (/p,q,r,ElemID/)
                                                              ! Second index: Local Tau_W calc point (with respect to each MPI proc. responsible for imposition)
                                                              ! Third index: MPI proc. to send info
REAL, ALLOCATABLE           :: TauW_MINE_Interpolate(:,:,:) ! h_wm coords in standard region for solution interpolation
                                                            ! First index: 1-3 (/xi,eta,zeta/) [-1,1]^3
                                                            ! Second index: Local Tau_W calc point (with respect to each MPI proc. responsible for imposition)
                                                            ! Third index: MPI proc. to send info
REAL, ALLOCATABLE           :: Lag_xi(:,:,:), Lag_eta(:,:,:), Lag_zeta(:,:,:) ! Pre-calculated Lagrangian polynomials to interpolate each h_wm point that cannot
                                                                              ! be approximated as a face or interior node.
                                                                              ! First index: 0-PP_N for each \ell_i, i \in [0,PP_N]
                                                                              ! Second index: Local Tau_W calc point (with respect to each MPI proc. responsible for imposition)
                                                                              ! Third index: MPI proc. to send info



!=================================================================================================================================
#endif                                                    

END MODULE MOD_WMLES_Vars