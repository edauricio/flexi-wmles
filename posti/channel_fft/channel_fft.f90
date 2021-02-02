!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "flexi.h"

!===================================================================================================================================
!> This tool will read in state files (single or multiple at once) from channel flow simulations and will perform
!> a FFT (using external libraries) to generate spectra as well as mean profiles for fluctuations and the mean velocity.
!> If several state files are given, an average over all of them is calculated.
!> The mesh needs to be ijk sorted.
!> To perform the FFT, the mesh and the solution will be interpolated to an equidistant FFT grid.
!===================================================================================================================================
PROGRAM channel_fft
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Commandline_Arguments
USE MOD_ReadInTools
USE MOD_StringTools,             ONLY: STRICMP,GetFileExtension
USE MOD_Mesh,                    ONLY: DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_Mesh_Readin,             ONLY: ReadIJKSorting
USE MOD_Mesh_Vars,               ONLY: nElems,OffsetElem,nBCSides,BC,BoundaryName,NormVec
USE MOD_MPI,                     ONLY: DefineParametersMPI,InitMPI
USE MOD_Interpolation,           ONLY: DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_IO_HDF5,                 ONLY: DefineParametersIO_HDF5,InitIOHDF5,File_ID
#if USE_MPI
USE MOD_MPI,                     ONLY: InitMPIvars,FinalizeMPI
#endif
USE MOD_HDF5_Input,              ONLY: OpenDataFile,CloseDataFile,GetDataProps,ReadAttribute,ReadArray
USE MOD_Interpolation_Vars,      ONLY: NodeType
USE MOD_DG_Vars,                 ONLY: U, UPrim_master
USE MOD_Lifting_Vars,            ONLY: gradUx_master, gradUy_master, gradUz_master
USE MOD_EOS_Vars,                ONLY: mu0
USE MOD_WMLES_Vars
USE MOD_FFT,                     ONLY: InitFFT,PerformFFT,FFTOutput,FinalizeFFT,PrimStateAtFFTCoords,ReadState
USE MOD_FFT_Vars,                ONLY: ProjectName,Time,Normalize,prmfile,WallBCName,Re_tauReal,u_tau
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iArg,iSide,i,k
INTEGER                            :: nElems_State,nVar_State,N_State
CHARACTER(LEN=255)                 :: MeshFile_state,NodeType_State
INTEGER                            :: nWallSides, MeshModeNorm=1, OppSide
INTEGER, ALLOCATABLE               :: WallToBCSide(:),WallToBCSide_tmp(:)
REAL                               :: tauSurf(3,3),GradV(3,3),DivV,rho_w,tau_w,sideCnt
REAL                               :: tauSurf2
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
IF (nProcessors.GT.1) CALL CollectiveStop(__STAMP__, &
     'This tool is designed only for single execution!')

CALL ParseCommandlineArguments()

! Define parameters needed
CALL DefineParametersInterpolation()
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
CALL DefineParametersMesh()

CALL prms%SetSection("channelFFT")
CALL prms%CreateIntOption( "OutputFormat",  "Choose the main format for output. 0: Tecplot, 2: HDF5")
CALL prms%CreateIntOption( "NCalc",  "Polynomial degree to perform DFFT on.")
CALL prms%CreateRealOption("Re_tau", "Reynolds number based on friction velocity and channel half height.")
CALL prms%CreateLogicalOption("Normalize", "Normalize the mean and fluctuations output with the friction velocity.", ".FALSE.")
CALL prms%CreateStringOption("WallBCName", "BC Name of the Wall to Compute Shear Stress")

! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF
! check if parameter file is given
IF ((nArgs.LT.1).OR.(.NOT.(STRICMP(GetFileExtension(Args(1)),'ini')))) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: channel_fft [prm-file] statefile [statefiles]')
END IF
! Parse parameters
CALL prms%read_options(Args(1))

! Set meshMode according to Normalize option. 
! If normalization is needed, mesh metrics must be computed
Normalize = GETLOGICAL('Normalize')
IF (Normalize) MeshModeNorm=2
IF (Normalize) WallBCName = GETSTR('WallBCName')

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(A)') &
"   _____ _    _          _   _ _   _ ______ _        ______ ______ _______ "
SWRITE(UNIT_stdOut,'(A)') &
"  / ____| |  | |   /\   | \ | | \ | |  ____| |      |  ____|  ____|__   __|"
SWRITE(UNIT_stdOut,'(A)') &
" | |    | |__| |  /  \  |  \| |  \| | |__  | |      | |__  | |__     | |   "
SWRITE(UNIT_stdOut,'(A)') &
" | |    |  __  | / /\ \ | . ` | . ` |  __| | |      |  __| |  __|    | |   "
SWRITE(UNIT_stdOut,'(A)') &
" | |____| |  | |/ ____ \| |\  | |\  | |____| |____  | |    | |       | |   "
SWRITE(UNIT_stdOut,'(A)') &
"  \_____|_|  |_/_/    \_\_| \_|_| \_|______|______| |_|    |_|       |_|   "
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(132("="))')


! Initialization
CALL InitIOHDF5()

! Open the first statefile to read necessary attributes
CALL OpenDataFile(Args(2),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadAttribute(File_ID,'MeshFile',    1,StrScalar=MeshFile_state)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)
CALL ReadAttribute(File_ID,'Time',        1,RealScalar=Time)
CALL GetDataProps(nVar_State,N_State,nElems_State,NodeType_State)
CALL CloseDataFile()
IF (.NOT.STRICMP(NodeType_State,NodeType)) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Node type of statefile not equal to node type FLEXI is compiled with')
END IF
postiMode=.TRUE.
CALL InitInterpolation(N_State)
CALL InitMesh(meshMode=MeshModeNorm,MeshFile_IN=MeshFile_state)
CALL ReadIJKSorting()
#if USE_MPI
CALL InitMPIvars()
#endif

! Now that the mesh has been read, we can populate our wall information
IF (Normalize) THEN
! Populate mapping between mesh wall sides and DG sides
  nWallSides=0
  ALLOCATE(WallToBCSide_tmp(nBCSides))
  DO iSide=1,nBCSides
    IF (STRICMP(BoundaryName(BC(iSide)),WallBCName)) THEN
      nWallSides = nWallSides + 1
      WallToBCSide_tmp(nWallSides) = iSide
    END IF
  END DO
  ALLOCATE(WallToBCSide(nWallSides))
  WallToBCSide(1:nWallSides) = WallToBCSide_tmp(1:nWallSides)
  DEALLOCATE(WallToBCSide_tmp)
END IF

CALL InitFFT()

! Loop over all statefiles
DO iArg=2,nArgs
  SWRITE(UNIT_stdOut,'(132("="))')
  SWRITE(UNIT_stdOut,'(A,I5,A,I5,A)') ' PROCESSING FILE ',iArg-1,' of ',nArgs-1,' FILES.'
  SWRITE(UNIT_stdOut,'(A,A,A)') ' ( "',TRIM(Args(iArg)),'" )'
  SWRITE(UNIT_stdOut,'(132("="))')

  ! Get Solution
  CALL OpenDataFile(Args(iArg),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  CALL ReadArray('DG_Solution',5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,nElems/),OffsetElem,5,RealArray=U)
  CALL CloseDataFile()

  ! Interpolate solution to FFT grid
  ! FFT grid is globally equidistant
  CALL PrimStateAtFFTCoords()

  ! Perform actual FFT
  CALL PerformFFT()
END DO

IF (Normalize) THEN
  WRITE(*,*) 'NORMALIZE '
  tauSurf = 0.
  rho_w = 0.
  tau_w = 0.
  tauSurf2 = 0.
  sideCnt = 0.

  ! Loop over all statefiles
  DO iArg=2,nArgs
    ! Read the state and compute gradients.
    ! We use a copy-and-paste version of this procedure from posti_visu
    ! The downside is all the colored output from each function; we will deal with this latter. LRP first
    CALL ReadState(prmfile,Args(iArg),MeshFile_state)
    
    ! Calculate shear stress tensor for each point on each wall side and sum them all up
    DO iSide=1,nWMLESSides
      IF (STRICMP(BoundaryName(BC(WMLESToBCSide(iSide))),WallBCName)) THEN
        IF (iArg.EQ.2) sideCnt = sideCnt + 1
        DO i=0,PP_N; DO k=0,PP_NZ
          tau_w = tau_w + WMLES_TauW(1,i,k,nWMLESSides)
          rho_w = rho_w + UPrim_master(1,i,k,WMLESToBCSide(iSide))
        END DO; END DO
      END IF
    END DO
    ! DO iSide=1,nWallSides
    !   WRITE(*,*) WMLES_TauW(1:2,1,1,1)
    !   CALL Abort(__STAMP__, 'Test done')
    !   DO i=0,PP_N
    !   ! DO i=1,N_FFT(1) !<-- used if we interpolate our solution to the FFT grid
    !     DO k=0,PP_NZ
    !     ! DO k=1,N_FFT(3) !<-- used if we interpolate our solution to the FFT grid
    !       GradV(1:3,1) = gradUx_master(2:4,i,k,WallToBCSide(iSide))
    !       GradV(1:3,2) = gradUy_master(2:4,i,k,WallToBCSide(iSide))
    !       GradV(1:3,3) = gradUz_master(2:4,i,k,WallToBCSide(iSide))
    !       DivV = GradV(1,1) + GradV(2,2) + GradV(3,3)
    !       tauSurf2 =  (mu0*(gradUy_master(2,i,k,WallToBCSide(iSide)) + gradUx_master(3,i,k,WallToBCSide(iSide))))
    !       tauSurf = tauSurf + mu0*(GradV + TRANSPOSE(GradV))
    !       tauSurf(1,1) = tauSurf(1,1) - (2./3.)*mu0*DivV
    !       tauSurf(2,2) = tauSurf(2,2) - (2./3.)*mu0*DivV
    !       tauSurf(3,3) = tauSurf(3,3) - (2./3.)*mu0*DivV
          
    !       !tau_w = tau_w  -1.*DOT_PRODUCT(tauSurf(1,1:3),NormVec(1:3,i,k,0,WallToBCSide(iSide)))
    !       tau_w = tau_w + tauSurf2!tauSurf(1,2)
    !       rho_w = rho_w + UPrim_master(1,i,k,WallToBCSide(iSide))
    !     END DO
    !   END DO
    ! END DO

  END DO
  ! Average tau_w
  !tauSurf = tauSurf / ((nArgs-1)*PP_N*PP_NZ*nWallSides)
  !tau_w = tau_w / ((nArgs-1)*(PP_N+1)*(PP_NZ+1)*nWMLESSides)
  !rho_w = rho_w / ((nArgs-1)*(PP_N+1)*(PP_NZ+1)*nWMLESSides)
  tau_w = ABS(tau_w) / ((nArgs-1)*(PP_N+1)*(PP_NZ+1)*sideCnt)
  rho_w = rho_w / ((nArgs-1)*(PP_N+1)*(PP_NZ+1)*sideCnt)
  !u_tau = SQRT(tauSurf/rho_w)
  u_tau = SQRT(tau_w/rho_w)
  Re_tauReal = 1. / ((mu0/rho_w)/u_tau)
  WRITE(*,*) "tau_w, rho_w, u_tau, Re_tau: ", tau_w, rho_w, u_tau, Re_tauReal

END IF

! Do output of results
CALL FFTOutput()

! Finalize
CALL FinalizeFFT()
CALL FinalizeMesh()
#if USE_MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL abort(__STAMP__,'MPI finalize error',iError)
#endif
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') ' Channel_FFT TOOL FINISHED! '
SWRITE(UNIT_stdOut,'(132("="))')

END PROGRAM channel_fft

