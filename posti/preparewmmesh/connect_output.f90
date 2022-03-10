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
!> Module containing the routines to write the already build connection information to a HDF5 file which can later be read
!> by FLEXI to set up the wall model calculation.
!===================================================================================================================================
MODULE MOD_Connect_Output
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE ConnectHDF5Output
  MODULE PROCEDURE ConnectHDF5Output
END INTERFACE

INTERFACE ConnectVisu
  MODULE PROCEDURE ConnectVisu
END INTERFACE

PUBLIC:: ConnectHDF5Output,ConnectVisu

CONTAINS

!===================================================================================================================================
!> Write the array wallconnect which contains all the necessary information to set up the wall model to a HDF5 file.
!> Some additional information like the number of wall modelled points, the polynomial degree and the node type are written as
!> well.
!===================================================================================================================================
SUBROUTINE ConnectHDF5Output()
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_PreProc
USE MOD_PrepareWMmesh_Vars
USE MOD_IO_HDF5
USE MOD_HDF5_Output,         ONLY: WriteAttribute,GatheredWriteArray
USE MOD_Mesh_Vars,           ONLY: MeshFile,nBCSides
USE MOD_Interpolation_Vars,  ONLY: NodeType
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
! Space-separated list of input and output types. Use: (int|real|logical|...)_(in|out|inout)_dim(n)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)        :: FileName
INTEGER                   :: iExt
INTEGER(HID_T)            :: DSet_ID,FileSpace,HDF5DataType
INTEGER(HSIZE_T)          :: Dimsf(4)
INTEGER(HSIZE_T)          :: Dimsf_2(1)
INTEGER(HSIZE_T)          :: Dimsf_3(2)
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)') ' Writing connection information to HDF5 file...'

! Create name of the connection file from the name of the mesh file
iExt=INDEX(MeshFile,'.',BACK = .TRUE.) ! Position of file extension
FileName = MeshFile(:iExt-6)
FileName = TRIM(FileName) // '_WM.h5'

!============ Create the HDF5 file and preallocate storage space ===========!

! Wallconnect arry

! Create file
CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)

! Preallocate the data space for the dataset.
Dimsf=(/15,PP_N+1,PP_N+1,nModelledBCSides/)

CALL H5SCREATE_SIMPLE_F(4, Dimsf, FileSpace, iError)
! Create the dataset with default properties.
HDF5DataType=H5T_NATIVE_DOUBLE
CALL H5DCREATE_F(File_ID,'WM', HDF5DataType, FileSpace, DSet_ID, iError)
! Close the filespace and the dataset
CALL H5DCLOSE_F(Dset_id, iError)
CALL H5SCLOSE_F(FileSpace, iError)

! Write dataset properties "N","MeshFile","NodeType" and number of modelled sides
CALL WriteAttribute(File_ID,'N',1,IntScalar=PP_N)
CALL WriteAttribute(File_ID,'MeshFile',1,StrScalar=(/MeshFile/))
CALL WriteAttribute(File_ID,'NodeType',1,StrScalar=(/NodeType/))
CALL WriteAttribute(File_ID,'nModelledBCSides',1,IntScalar=nModelledBCSides)

CALL CloseDataFile()

! Mapping

! Create file
CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)

! Preallocate the data space for the dataset.
Dimsf_2=(/nBCSides/)

CALL H5SCREATE_SIMPLE_F(1, Dimsf_2, FileSpace, iError)
! Create the dataset with default properties.
HDF5DataType=H5T_NATIVE_INTEGER
CALL H5DCREATE_F(File_ID,'MappingWM', HDF5DataType, FileSpace, DSet_ID, iError)
! Close the filespace and the dataset
CALL H5DCLOSE_F(Dset_id, iError)
CALL H5SCLOSE_F(FileSpace, iError)

CALL CloseDataFile()

! Shape Boundary Info

! Create file
CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)

! Preallocate the data space for the dataset.
Dimsf_3=(/2,nModelledBCSides/)

CALL H5SCREATE_SIMPLE_F(2, Dimsf_3, FileSpace, iError)
! Create the dataset with default properties.
HDF5DataType=H5T_NATIVE_INTEGER
CALL H5DCREATE_F(File_ID,'sideShapeInfo', HDF5DataType, FileSpace, DSet_ID, iError)
! Close the filespace and the dataset
CALL H5DCLOSE_F(Dset_id, iError)
CALL H5SCLOSE_F(FileSpace, iError)

CALL CloseDataFile()

!============================ Actual data output ===========================!

CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName='WM', rank=4,&
                        nValGlobal=(/15,PP_N+1,PP_N+1,nModelledBCSides/),&
                        nVal=      (/15,PP_N+1,PP_N+1,nModelledBCSides/),&
                        offset=    (/0,      0,     0,               0/),&
                        collective=.TRUE.,RealArray=wallconnect)


CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName='MappingWM', rank=1,&
                        nValGlobal=(/nBCSides/),&
                        nVal=      (/nBCSides/),&
                        offset=    (/       0/),&
                        collective=.TRUE.,IntArray=mapBCSideToModelledSide)

CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName='sideShapeInfo',rank=2,&
                        nValGlobal=(/2,nModelledBCSides/),&
                        nVal=      (/2,nModelledBCSides/),&
                        offset=    (/0,               0/),&
                        collective=.TRUE.,IntArray=sideShapeInfo)
SWRITE(UNIT_stdOut,'(A)')' DONE!'


END SUBROUTINE ConnectHDF5Output


!===================================================================================================================================
!> Visualize the result of the wall connection process. Since we need to visualize surface data here, this is restricted
!> in some ways. To visualize the data, the surface data will be stored in the cell next to boundary point (constant along the
!> wall-normal direction. Thus only one wall-modelled boundary side per cell can be visualized!
!===================================================================================================================================
SUBROUTINE ConnectVisu()
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_PreProc
USE MOD_PrepareWMmesh_Vars
USE MOD_Mesh_Vars,             ONLY: SideToElem,nBCSides,MeshFile,Elem_xGP,nElems,S2V2
USE MOD_VTK,                   ONLY: WriteDataToVTK
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,POINTER         :: Wallconnect_Volume(:,:,:,:,:)
REAL,POINTER         :: Coord_p(:,:,:,:,:)
INTEGER              :: ElemID,locSideID,iBC,iModelledSide,iExt,i,p,q
CHARACTER(LEN=255)   :: FileName
INTEGER              :: pq(2)
!===================================================================================================================================
! Initialize the wallconnect array with zeros
ALLOCATE(Wallconnect_Volume(17,0:PP_N,0:PP_N,0:PP_N,nElems))
Wallconnect_Volume = 0.

! Loop over all BC sides
DO iBC=1,nBCSides
  iModelledSide = mapBCSideToModelledSide(iBC)
  IF (iModelledSide.GT.0) THEN
    ! Write the surface solution to the next cell, extruded along the wall-normal direction
    ! We need to rotate the array from side to volume coordinate system!
    ElemID    = SideToElem(S2E_ELEM_ID,iBC)
    locSideID = SideToElem(S2E_LOC_SIDE_ID,iBC)
    SELECT CASE(locSideID)
    CASE(XI_MINUS,XI_PLUS)
      DO q=0,PP_N; DO p=0,PP_N
        pq(1) = S2V2(1,p,q,0,locSideID)
        pq(2) = S2V2(2,p,q,0,locSideID)
        DO i=1,15
          Wallconnect_Volume(i,:,pq(1),pq(2),ElemID) = wallconnect(i,p,q,iModelledSide)
        END DO ! i=1,15
        Wallconnect_Volume(16,:,pq(1),pq(2),ElemID) = sideShapeInfo(INTERFACE_SHAPE, iModelledSide)
        Wallconnect_Volume(17,:,pq(1),pq(2),ElemID) = sideShapeInfo(SHAPE_BOUNDARY, iModelledSide)
      END DO; END DO ! p,q=0,PP_N
    CASE(ETA_MINUS,ETA_PLUS)
      DO q=0,PP_N; DO p=0,PP_N
        pq(1) = S2V2(1,p,q,0,locSideID)
        pq(2) = S2V2(2,p,q,0,locSideID)
        DO i=1,15
          Wallconnect_Volume(i,pq(1),:,pq(2),ElemID) = wallconnect(i,p,q,iModelledSide)
        END DO ! i=1,15
        Wallconnect_Volume(16,pq(1),:,pq(2),ElemID) = sideShapeInfo(INTERFACE_SHAPE, iModelledSide)
        Wallconnect_Volume(17,pq(1),:,pq(2),ElemID) = sideShapeInfo(SHAPE_BOUNDARY, iModelledSide)
      END DO; END DO ! p,q=0,PP_N
    CASE(ZETA_MINUS,ZETA_PLUS)
      DO q=0,PP_N; DO p=0,PP_N
        pq(1) = S2V2(1,p,q,0,locSideID)
        pq(2) = S2V2(2,p,q,0,locSideID)
        DO i=1,15
          Wallconnect_Volume(i,pq(1),:,pq(2),ElemID) = wallconnect(i,p,q,iModelledSide)
        END DO ! i=1,15
        Wallconnect_Volume(16,pq(1),:,pq(2),ElemID) = sideShapeInfo(INTERFACE_SHAPE, iModelledSide)
        Wallconnect_Volume(17,pq(1),:,pq(2),ElemID) = sideShapeInfo(SHAPE_BOUNDARY, iModelledSide)
      END DO; END DO ! p,q=0,PP_N
    END SELECT
  END IF
END DO

! Create name of the connection file from the name of the mesh file
iExt=INDEX(MeshFile,'.',BACK = .TRUE.) ! Position of file extension
FileName = MeshFile(:iExt-6)
FileName = TRIM(FileName) // '_WM.vtu'

! Write the wallconnect information to vtu
Coord_p => Elem_xGP
CALL WriteDataToVTK(17,PP_N,nElems,StrVarNamesWallconnect,Coord_p,Wallconnect_Volume,TRIM(FileName),dim=3)
DEALLOCATE(Wallconnect_Volume)

END SUBROUTINE ConnectVisu

END MODULE MOD_Connect_Output
