!> \file
!> \author Chris Bradley
!> \brief This module handles all generated mesh routines.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
!>
!> Contributor(s): Chris Bradley, Jack Lee
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> This module handles all generated mesh routines.
MODULE GeneratedMeshRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE ComputationRoutines
  USE ComputationAccessRoutines
  USE Constants
  USE CoordinateSystemAccessRoutines
  USE DecompositionRoutines
  USE DecompositionAccessRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE GeneratedMeshAccessRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE Maths
  USE MeshRoutines
  USE MeshAccessRoutines
  USE NodeRoutines
  USE Strings
  USE RegionAccessRoutines
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters


  !Module types

  !Interfaces

  !>Starts the process of creating a generated mesh
  INTERFACE GeneratedMesh_CreateStart
    MODULE PROCEDURE GeneratedMesh_CreateStartInterface
    MODULE PROCEDURE GeneratedMesh_CreateStartRegion
  END INTERFACE GeneratedMesh_CreateStart

  !>Initialises the generated meshes for a region or interface.
  INTERFACE GeneratedMeshes_Initialise
    MODULE PROCEDURE GeneratedMeshes_InitialiseInterface
    MODULE PROCEDURE GeneratedMeshes_InitialiseRegion
  END INTERFACE GeneratedMeshes_Initialise

  PUBLIC GeneratedMesh_BaseVectorsSet

  PUBLIC GeneratedMesh_BasisSet

  PUBLIC GeneratedMesh_CreateStart,GeneratedMesh_CreateFinish

  PUBLIC GeneratedMesh_Destroy

  PUBLIC GeneratedMesh_ExtentSet

  PUBLIC GeneratedMesh_NumberOfElementsSet

  PUBLIC GeneratedMesh_OriginSet

  PUBLIC GeneratedMesh_TypeSet

  PUBLIC GeneratedMesh_GeometricParametersCalculate

  PUBLIC GeneratedMesh_SurfaceGet
  
  PUBLIC GeneratedMeshes_Initialise,GeneratedMeshes_Finalise

CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets/changes the basis of a generated mesh. \see OpenCMISS::Iron::cmfe_GeneratedMesh_BasisSet
  SUBROUTINE GeneratedMesh_BasisSet(generatedMesh,bases,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to set the basis of
    TYPE(BasisPtrType) :: bases(:) !<An array of pointers to the basis to generate the mesh with
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coordinateDimension,basisIdx,firstNumberOfXi,firstMeshBasisType,numberOfBases,numberOfXi,meshBasisType
    TYPE(BasisType), POINTER :: basis
    TYPE(BasisPtrType), ALLOCATABLE :: newBases(:)
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ellipsoidMesh
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_BasisSet",err,error,*999)

    CALL GeneratedMesh_AssertNotFinished(generatedMesh,err,error,*999)
    
    NULLIFY(coordinateSystem)
    CALL GeneratedMesh_CoordinateSystemGet(generatedMesh,coordinateSystem,err,error,*999)
    CALL CoordinateSystem_DimensionGet(coordinateSystem,coordinateDimension,err,error,*999)
    NULLIFY(basis)
    CALL GeneratedMeshRegular_BasisGet(regularMesh,1,basis,err,error,*999)
    numberOfBases=SIZE(bases)
    CALL Basis_NumberOfXiGet(basis,firstNumberOfXi,err,error,*999)
    CALL Basis_TypeGet(basis,firstMeshBasisType,err,error,*999)
    DO basisIdx=2,numberOfBases
      NULLIFY(basis)
      CALL GeneratedMeshRegular_BasisGet(regularMesh,basisIdx,basis,err,error,*999)
      CALL Basis_NumberOfXiGet(basis,numberOfXi,err,error,*999)
      CALL Basis_TypeGet(basis,meshBasisType,err,error,*999)
      IF(numberOfXi>coordinateDimension) THEN
        localError="The basis number of xi dimensions of "// &
          & TRIM(NumberToVString(numberOfXi,"*",err,error))// &
          & " is invalid. The number of xi dimensions must be <= the number of coordinate dimensions of "// &
          & TRIM(NumberToVString(coordinateDimension,"*",err,error))
        CALL FlagError(localError,err,error,*999)
      ENDIF      
      IF(numberOfXi/=firstNumberOfXi) CALL FlagError("All bases must have the same number of xi.",err,error,*999)
      IF(meshBasisType/=firstMeshBasisType) &
        & CALL FlagError("Using different basis types is not supported for generated meshes.",err,error,*999)
    ENDDO !basisIdx
    ALLOCATE(newBases(numberOfBases),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate bases.",err,error,*999)
    DO basisIdx=1,numberOfBases
      newBases(basisIdx)%ptr=>bases(basisIdx)%ptr
    ENDDO !basisIdx    
    SELECT CASE(generatedMesh%generatedType)
    CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
      NULLIFY(regularMesh)
      CALL GeneratedMesh_RegularMeshGet(generatedMesh,regularMesh,err,error,*999)
      IF(ALLOCATED(regularMesh%baseVectors)) &
        & CALL FlagError("Can not reset the basis if base vectors have been specified.",err,error,*999)
      CALL MOVE_ALLOC(newBases,regularMesh%bases)
    CASE(GENERATED_MESH_POLAR_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
      NULLIFY(cylinderMesh)
      CALL GeneratedMesh_CylinderMeshGet(generatedMesh,cylinderMesh,err,error,*999)
      CALL MOVE_ALLOC(newBases,cylinderMesh%bases)
    CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
      NULLIFY(ellipsoidMesh)
      CALL GeneratedMesh_EllipsoidMeshGet(generatedMesh,ellipsoidMesh,err,error,*999)
      CALL MOVE_ALLOC(newBases,ellipsoidMesh%bases)
    CASE DEFAULT
      localError="The generated mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("GeneratedMesh_BasisSet")
    RETURN
999 IF(ALLOCATED(newBases)) DEALLOCATE(newBases)
    ERRORSEXITS("GeneratedMesh_BasisSet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_BasisSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the base vectors of a generated mesh. \see OpenCMISS::Iron::cmfe_GeneratedMesh_BaseVectorsSet
  SUBROUTINE GeneratedMesh_BaseVectorsSet(generatedMesh,baseVectors,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to set the base vectors fo
    REAL(DP), INTENT(IN) :: baseVectors(:,:) !<baseVectors(coordinateIdx,xiIdx). The base vectors for the generated mesh to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coordinateDimension,numberOfXi
    REAL(DP), ALLOCATABLE :: newBaseVectors(:,:)
    TYPE(BasisType), POINTER :: basis
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_BaseVectorsSet",err,error,*999)

    CALL GeneratedMesh_AssertNotFinished(generatedMesh,err,error,*999)
    
    NULLIFY(coordinateSystem)
    CALL GeneratedMesh_CoordinateSystemGet(generatedMesh,coordinateSystem,err,error,*999)
    CALL CoordinateSystem_DimensionGet(coordinateSystem,coordinateDimension,err,error,*999)
    IF(SIZE(baseVectors,1)/=coordinateDimension) THEN
      localError="The size of the first dimension of base vectors of "// &
        & TRIM(NumberToVString(SIZE(baseVectors,1),"*",err,error))// &
        & " is invalid. The first dimension size must match the coordinate system dimension of "// &
        & TRIM(NumberToVString(coordinateDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(generatedMesh%generatedType)
    CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
      NULLIFY(regularMesh)
      CALL GeneratedMesh_RegularMeshGet(generatedMesh,regularMesh,err,error,*999)
      CALL GeneratedMeshRegular_BasisGet(regularMesh,1,basis,err,error,*999)
      CALL Basis_NumberOfXiGet(basis,numberOfXi,err,error,*999)
      IF(SIZE(baseVectors,2)==numberOfXi) THEN
        localError="The size of the second dimension of base vectors of "// &
          & TRIM(NumberToVString(SIZE(baseVectors,2),"*",err,error))// &
          & " is invalid. The second dimension size must match the number of mesh dimensions of "// &
          & TRIM(NumberToVString(numberOfXi,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      ALLOCATE(newBaseVectors(SIZE(baseVectors,1),SIZE(baseVectors,2)),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate base vectors.",err,error,*999)      
      newBaseVectors=baseVectors
      CALL MOVE_ALLOC(newBaseVectors,regularMesh%baseVectors)
    CASE(GENERATED_MESH_POLAR_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The generated mesh mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("GeneratedMesh_BaseVectorsSet")
    RETURN
999 ERRORSEXITS("GeneratedMesh_BaseVectorsSet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_BaseVectorsSet

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a generated mesh. \see OpenCMISS::Iron::cmfe_GeneratedMesh_CreateFinish
  SUBROUTINE GeneratedMesh_CreateFinish(generatedMesh,meshUserNumber,mesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to finish the creation of
    INTEGER(INTG), INTENT(IN) :: meshUserNumber !<The mesh's user number
    TYPE(MeshType), POINTER :: mesh !<On exit, a pointer to the generated mesh. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_CreateFinish",err,error,*998)

    IF(ASSOCIATED(mesh)) CALL FlagError("Mesh is already associated.",err,error,*998)
    CALL GeneratedMesh_AssertNotFinished(generatedMesh,err,error,*999)
    
    SELECT CASE(generatedMesh%generatedType)
    CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
      CALL GeneratedMesh_RegularCreateFinish(generatedMesh,meshUserNumber,err,error,*999)
    CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
      CALL GeneratedMesh_CylinderCreateFinish(generatedMesh,meshUserNumber,err,error,*999)
    CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
      CALL GeneratedMesh_EllipsoidCreateFinish(generatedMesh,meshUserNumber,err,error,*999)
    CASE(GENERATED_MESH_POLAR_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The generated mesh mesh type of "// &
        & TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    !Return the pointers
    mesh=>generatedMesh%mesh
    mesh%generatedMesh=>generatedMesh
    generatedMesh%generatedMeshFinished=.TRUE.

    EXITS("GeneratedMesh_CreateFinish")
    RETURN
999 NULLIFY(mesh)
998 ERRORSEXITS("GeneratedMesh_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_CreateFinish

  !
  !================================================================================================================================
  !

  !>Starts the creation of a generic generated mesh.
  SUBROUTINE GeneratedMesh_CreateStartGeneric(generatedMeshes,userNumber,generatedMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshesType), POINTER :: generatedMeshes !<A pointer to the generated meshes
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the generated mesh to create
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<On exit, a pointer to the created generated mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,generatedMeshIdx
    TYPE(GeneratedMeshType), POINTER :: newGeneratedMesh
    TYPE(GeneratedMeshPtrType), ALLOCATABLE :: newGeneratedMeshes(:)
    TYPE(VARYING_STRING) :: dummyError

    NULLIFY(newGeneratedMesh)

    ENTERS("GeneratedMesh_CreateStartGeneric",err,error,*997)

    IF(ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is already associated.",err,error,*997)
    IF(ASSOCIATED(generatedMeshes)) CALL FlagError("Generated meshes is not associated.",err,error,*997)
    
    !Initialise generated mesh
    CALL GeneratedMesh_Initialise(newGeneratedMesh,err,error,*999)
    !Set default generated mesh values
    newGeneratedMesh%userNumber=userNumber
    newGeneratedMesh%globalNumber=generatedMeshes%numberOfGeneratedMeshes+1
    newGeneratedMesh%generatedMeshes=>generatedMeshes
    !Add new generated mesh into list of generated meshes
    ALLOCATE(newGeneratedMeshes(generatedMeshes%numberOfGeneratedMeshes+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new generated meshes.",err,error,*999)
    DO generatedMeshIdx=1,generatedMeshes%numberOfGeneratedMeshes
      newGeneratedMeshes(generatedMeshIdx)%ptr=>generatedMeshes%generatedMeshes(generatedMeshIdx)%ptr
    ENDDO !generatedMeshIdx
    newGeneratedMeshes(generatedMeshes%numberOfGeneratedMeshes+1)%ptr=>newGeneratedMesh
    CALL MOVE_ALLOC(newGeneratedMeshes,generatedMeshes%generatedMeshes)
    generatedMeshes%numberOfGeneratedMeshes=generatedMeshes%numberOfGeneratedMeshes+1
    !Return the pointer
    generatedMesh=>newGeneratedMesh
 
    EXITS("GeneratedMesh_CreateStartGeneric")
    RETURN
999 CALL GeneratedMesh_Finalise(newGeneratedMesh,dummyErr,dummyError,*998)
998 IF(ALLOCATED(newGeneratedMeshes)) DEALLOCATE(newGeneratedMeshes)
    NULLIFY(generatedMesh)
997 ERRORSEXITS("GeneratedMesh_CreateStartGeneric",err,error)
    RETURN 1

  END SUBROUTINE GeneratedMesh_CreateStartGeneric

  !
  !================================================================================================================================
  !

  !>Starts the creation of a generated mesh. \see OpenCMISS::Iron::cmfe_GeneratedMesh_CreateFinish
  SUBROUTINE GeneratedMesh_CreateStartInterface(userNumber,interface,generatedMesh,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the generated mesh to create
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to create the generated mesh on
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<On exit, a pointer to the created generated mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_CreateStartInterface",err,error,*998)

    IF(ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*999)
    
    CALL GeneratedMesh_UserNumberFind(userNumber,interface,generatedMesh,err,error,*999)
    IF(ASSOCIATED(generatedMesh)) THEN
      localError="The specified user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " has already been used for a generated mesh on interface number "// &
        & TRIM(NumberToVString(interface%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    CALL GeneratedMesh_CreateStartGeneric(interface%generatedMeshes,userNumber,generatedMesh,err,error,*999)
    generatedMesh%interface=>interface
 
    EXITS("GeneratedMesh_CreateStartInterface")
    RETURN
999 NULLIFY(generatedMesh)
998 ERRORSEXITS("GeneratedMesh_CreateStartInterface",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_CreateStartInterface

  !
  !================================================================================================================================
  !

  !>Starts the creation of a generated mesh. \see OpenCMISS::Iron::cmfe_GeneratedMesh_CreateFinish
  SUBROUTINE GeneratedMesh_CreateStartRegion(userNumber,region,generatedMesh,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the generated mesh to create
    TYPE(RegionType), POINTER :: region !<A pointer to the region to create the generated mesh on
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<On exit, a pointer to the created generated mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_CreateStartRegion",err,error,*998)

    IF(ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is already associated.",err,error,*998)
    IF(ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)
    
    NULLIFY(generatedMesh)
    CALL GeneratedMesh_UserNumberFind(userNumber,region,generatedMesh,err,error,*999)
    IF(ASSOCIATED(generatedMesh)) THEN
      localError="The specified user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " has already been used for a generated mesh on region number "// &
        &  TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    CALL GeneratedMesh_CreateStartGeneric(region%generatedMeshes,userNumber,generatedMesh,err,error,*999)
    generatedMesh%region=>region

    EXITS("GeneratedMesh_CreateStartRegion")
    RETURN
999 NULLIFY(generatedMesh)
998 ERRORSEXITS("GeneratedMesh_CreateStartRegion",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_CreateStartRegion

  !
  !================================================================================================================================
  !

  !>Destroys a generated mesh. \see OpenCMISS::Iron::cmfe_GeneratedMesh_CreateDestroy
  SUBROUTINE GeneratedMesh_Destroy(generatedMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: generatedMeshIdx,generatedMeshPosition
    TYPE(GeneratedMeshPtrType), ALLOCATABLE :: newGeneratedMeshes(:)
    TYPE(GeneratedMeshesType), POINTER :: generatedMeshes

    ENTERS("GeneratedMesh_Destroy",err,error,*998)

    IF(ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is not associated.",err,error,*998)

    NULLIFY(generatedMeshes)
    CALL GeneratedMesh_GeneratedMeshesGet(generatedMesh,generatedMeshes,err,error,*999)    
    IF(.NOT.ALLOCATED(generatedMeshes%generatedMeshes))  &
      & CALL FlagError("Generated mesh generated meshes is not allocated.",err,error,*998)
    
    generatedMeshPosition=generatedMesh%globalNumber
    CALL GeneratedMesh_Finalise(generatedMesh,err,error,*999)
    !Remove the generated mesh from the list of generated meshes
    IF(generatedMeshes%numberOfGeneratedMeshes>1) THEN
      ALLOCATE(newGeneratedMeshes(generatedMeshes%numberOfGeneratedMeshes-1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new generated meshes.",err,error,*999)
      DO generatedMeshIdx=1,generatedMeshes%numberOfGeneratedMeshes
        IF(generatedMeshIdx<generatedMeshPosition) THEN
          newGeneratedMeshes(generatedMeshIdx)%ptr=>generatedMeshes%generatedMeshes(generatedMeshIdx)%ptr
        ELSE IF(generatedMeshIdx>generatedMeshPosition) THEN
          generatedMeshes%generatedMeshes(generatedMeshIdx)%ptr%globalNumber=generatedMeshes% &
            & generatedMeshes(generatedMeshIdx)%ptr%globalNumber-1
          newGeneratedMeshes(generatedMeshIdx-1)%ptr=>generatedMeshes%generatedMeshes(generatedMeshIdx)%ptr
        ENDIF
      ENDDO !generatedMeshIdx
      CALL MOVE_ALLOC(newGeneratedMeshes,generatedMeshes%generatedMeshes)
      generatedMeshes%numberOfGeneratedMeshes=generatedMeshes%numberOfGeneratedMeshes-1
    ELSE
      DEALLOCATE(generatedMeshes%generatedMeshes)
      generatedMeshes%numberOfGeneratedMeshes=0
    ENDIF
 
    EXITS("GeneratedMesh_Destroy")
    RETURN
999 IF(ALLOCATED(newGeneratedMeshes)) DEALLOCATE(newGeneratedMeshes)
998 ERRORSEXITS("GeneratedMesh_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_Destroy

  !
  !================================================================================================================================
  !

  !>Sets/changes the extent of a generated mesh. \see OpenCMISS::Iron::cmfe_GeneratedMesh_ExtentSet
  SUBROUTINE GeneratedMesh_ExtentSet(generatedMesh,extent,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to set the type of
    REAL(DP), INTENT(IN) :: extent(:) !<The extent of the generated mesh (MAXIMUM for regular type, CYLINDER for cylinder type)
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coordinateDimension
    REAL(DP) :: norm
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ellipsoidMesh
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_ExtentSet",err,error,*999)

    CALL GeneratedMesh_AssertNotFinished(generatedMesh,err,error,*999)
    NULLIFY(coordinateSystem)
    CALL GeneratedMesh_CoordinateSystemGet(generatedMesh,coordinateSystem,err,error,*999)
    CALL CoordinateSystem_DimensionGet(coordinateSystem,coordinateDimension,err,error,*999)
    SELECT CASE(generatedMesh%generatedType)
    CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
      IF(SIZE(extent,1)/=coordinateDimension) THEN
        localError="The extent size of "//TRIM(NumberToVString(SIZE(extent,1),"*",err,error))// &
          & " is invalid. The extent size must match the coordinate system dimension of "// &
          & TRIM(NumberToVString(coordinateDimension,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      CALL L2Norm(extent,norm,err,error,*999)
      IF(norm<=ZERO_TOLERANCE) CALL FlagError("The norm of the mesh extent is zero.",err,error,*999)
      NULLIFY(regularMesh)
      CALL GeneratedMesh_RegularMeshGet(generatedMesh,regularMesh,err,error,*999)
      IF(ALLOCATED(regularMesh%maximumExtent)) DEALLOCATE(regularMesh%maximumExtent)
      ALLOCATE(regularMesh%maximumExtent(coordinateDimension),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate maximum extent.",err,error,*999)
      regularMesh%maximumExtent=extent
    CASE(GENERATED_MESH_POLAR_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
      IF(SIZE(extent,1)/=coordinateDimension) THEN
        localError="The extent size of "//TRIM(NumberToVString(SIZE(extent,1),"*",err,error))// &
          & " is invalid. The extent size must match the coordinate system dimension of "// &
          & TRIM(NumberToVString(coordinateDimension,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(cylinderMesh)
      CALL GeneratedMesh_CylinderMeshGet(generatedMesh,cylinderMesh,err,error,*999)
      IF(ALLOCATED(cylinderMesh%cylinderExtent)) DEALLOCATE(cylinderMesh%cylinderExtent)
      ALLOCATE(cylinderMesh%cylinderExtent(SIZE(extent)),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate maximum extent.",err,error,*999)
      cylinderMesh%cylinderExtent=extent
    CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
      IF((SIZE(extent,1)-1)/=coordinateDimension) THEN
        localError="The extent size of "//TRIM(NumberToVString(SIZE(extent,1),"*",err,error))// &
          & " is invalid. The extent size must be equal one plus the coordinate system dimension of "// &
          & TRIM(NumberToVString(coordinateDimension,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(ellipsoidMesh)
      CALL GeneratedMesh_EllipsoidMeshGet(generatedMesh,ellipsoidMesh,err,error,*999)
      IF(ALLOCATED(ellipsoidMesh%ellipsoidExtent)) DEALLOCATE(ellipsoidMesh%ellipsoidExtent)
      ALLOCATE(ellipsoidMesh%ellipsoidExtent(SIZE(extent)),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate maximum extent.",err,error,*999)
      ellipsoidMesh%ellipsoidExtent=extent
    CASE DEFAULT
      localError="The generated mesh mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("GeneratedMesh_ExtentSet")
    RETURN
999 ERRORSEXITS("GeneratedMesh_ExtentSet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_ExtentSet

  !
  !================================================================================================================================
  !

  !>Get one of the surfaces of a generated mesh. \see OpenCMISS::Iron::cmfe_GeneratedMesh_SurfaceGet
  SUBROUTINE GeneratedMesh_SurfaceGet(generatedMesh,meshComponent,surfaceType,surfaceNodes,normalXi,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to get the type of
    INTEGER(INTG), INTENT(IN) :: meshComponent !<The mesh component to get the surface for.
    INTEGER(INTG), INTENT(IN) :: surfaceType !<The surface you are interested in
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: surfaceNodes (:) !<The nodes on the specified surface
    INTEGER(INTG), INTENT(OUT) :: normalXi !<The normal outward pointing xi direction
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ellipsoidMesh
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_SurfaceGet",err,error,*999)

    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is not associated",err,error,*999)
    SELECT CASE(generatedMesh%generatedType)
    CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
      NULLIFY(regularMesh)
      CALL GeneratedMesh_RegularMeshGet(generatedMesh,regularMesh,err,error,*999)
      CALL GeneratedMesh_RegularSurfaceGet(regularMesh,meshComponent,surfaceType,surfaceNodes,normalXi,err,error,*999)
    CASE(GENERATED_MESH_POLAR_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
      NULLIFY(cylinderMesh)
      CALL GeneratedMesh_CylinderMeshGet(generatedMesh,cylinderMesh,err,error,*999)
      CALL GeneratedMesh_CylinderSurfaceGet(cylinderMesh,meshComponent,surfaceType,surfaceNodes,normalXi,err,error,*999)
    CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
      NULLIFY(ellipsoidMesh)
      CALL GeneratedMesh_EllipsoidMeshGet(generatedMesh,ellipsoidMesh,err,error,*999)
      CALL GeneratedMesh_EllipsoidSurfaceGet(ellipsoidMesh,meshComponent,surfaceType,surfaceNodes,normalXi,err,error,*999)
    CASE DEFAULT
      localError="The generated mesh mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("GeneratedMesh_SurfaceGet")
    RETURN
999 ERRORSEXITS("GeneratedMesh_SurfaceGet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_SurfaceGet

  !
  !================================================================================================================================
  !

  !>Finalises a generated mesh and dellocates all memory.
  SUBROUTINE GeneratedMesh_Finalise(generatedMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("GeneratedMesh_Finalise",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) THEN
      CALL GeneratedMesh_RegularFinalise(generatedMesh%regularMesh,err,error,*999)
      CALL GeneratedMesh_CylinderFinalise(generatedMesh%cylinderMesh,err,error,*999)
      CALL GeneratedMesh_EllipsoidFinalise(generatedMesh%ellipsoidMesh,err,error,*999)
      DEALLOCATE(generatedMesh)
    ENDIF

    EXITS("GeneratedMesh_Finalise")
    RETURN
999 ERRORSEXITS("GeneratedMesh_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises a generated mesh.
  SUBROUTINE GeneratedMesh_Initialise(generatedMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("GeneratedMesh_Initialise",err,error,*999)

    IF(ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is already associated.",err,error,*999)
    
    ALLOCATE(generatedMesh,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate generated mesh.",err,error,*999)
    generatedMesh%userNumber=0
    generatedMesh%globalNumber=0
    generatedMesh%generatedMeshFinished=.FALSE.
    NULLIFY(generatedMesh%region)
    NULLIFY(generatedMesh%interface)
    generatedMesh%generatedType=0
    NULLIFY(generatedMesh%regularMesh)
    NULLIFY(generatedMesh%cylinderMesh)
    NULLIFY(generatedMesh%ellipsoidMesh)
    NULLIFY(generatedMesh%MESH)
    !Default to a regular mesh.
    CALL GeneratedMesh_RegularInitialise(generatedMesh,err,error,*999)

    EXITS("GeneratedMesh_Initialise")
    RETURN
999 ERRORSEXITS("GeneratedMesh_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_Initialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of elements in a generated mesh. \see OpenCMISS::Iron::cmfe_GeneratedMesh_NumberOfElementsSet
  SUBROUTINE GeneratedMesh_NumberOfElementsSet(generatedMesh,numberOfElementsXi,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to set the type of
    INTEGER(INTG), INTENT(IN) :: numberOfElementsXi(:) !<numberOfElementsXi(ni). The number of elements in the ni'th xi direction (or, r,theta,z for cylinder) to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfXi
    TYPE(BasisType), POINTER :: basis
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ellipsoidMesh
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_NumberOfElementsSet",err,error,*999)

    CALL GeneratedMesh_AssertNotFinished(generatedMesh,err,error,*999)

    SELECT CASE(generatedMesh%generatedType)
    CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
      NULLIFY(regularMesh)
      CALL GeneratedMesh_RegularMeshGet(generatedMesh,regularMesh,err,error,*999)
      NULLIFY(basis)
      CALL GeneratedMeshRegular_BasisGet(regularMesh,1,basis,err,error,*999)
      CALL Basis_NumberOfXiGet(basis,numberOfXi,err,error,*999)
      IF(SIZE(numberOfElementsXi,1)/=numberOfXi) THEN
        localError="The number of elements xi size of "// &
          & TRIM(NumberToVString(SIZE(numberOfElementsXi,1),"*",err,error))// &
          & " is invalid. The number of elements xi size must match the basis number of xi dimensions of "// &
          & TRIM(NumberToVString(numberOfXi,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(ANY(numberOfElementsXi<=0)) CALL FlagError("Must have 1 or more elements in all directions.",err,error,*999)
      IF(ALLOCATED(regularMesh%numberOfElementsXi)) DEALLOCATE(regularMesh%numberOfElementsXi)
      ALLOCATE(regularMesh%numberOfElementsXi(numberOfXi),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate number of elements xi.",err,error,*999)
      regularMesh%numberOfElementsXi(1:numberOfXi)=numberOfElementsXi(1:numberOfXi)
    CASE(GENERATED_MESH_POLAR_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
      NULLIFY(cylinderMesh)
      CALL GeneratedMesh_CylinderMeshGet(generatedMesh,cylinderMesh,err,error,*999)
      IF(ALLOCATED(cylinderMesh%numberOfElementsXi)) DEALLOCATE(cylinderMesh%numberOfElementsXi)
      ALLOCATE(cylinderMesh%numberOfElementsXi(SIZE(numberOfElementsXi,1)),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate number of elements xi.",err,error,*999)
      cylinderMesh%numberOfElementsXi(1:SIZE(numberOfElementsXi,1))=numberOfElementsXi(1:SIZE(numberOfElementsXi,1))
    CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
      NULLIFY(ellipsoidMesh)
      CALL GeneratedMesh_EllipsoidMeshGet(generatedMesh,ellipsoidMesh,err,error,*999)
      IF(ALLOCATED(ellipsoidMesh%numberOfElementsXi)) DEALLOCATE(ellipsoidMesh%numberOfElementsXi)
      ALLOCATE(ellipsoidMesh%numberOfElementsXi(SIZE(numberOfElementsXi,1)),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate number of elements xi.",err,error,*999)
      ellipsoidMesh%numberOfElementsXi(1:SIZE(numberOfElementsXi,1))=numberOfElementsXi(1:SIZE(numberOfElementsXi,1))
    CASE DEFAULT
      localError="The generated mesh mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("GeneratedMesh_NumberOfElementsSet")
    RETURN
999 ERRORSEXITS("GeneratedMesh_NumberOfElementsSet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_NumberOfElementsSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the origin of a generated mesh. \see OpenCMISS::Iron::cmfe_GeneratedMesh_OriginSet
  SUBROUTINE GeneratedMesh_OriginSet(generatedMesh,origin,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to set the type of
    REAL(DP), INTENT(IN) :: origin(:) !<origin(nj). The nj'th coordinate origin for the generated mesh
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coordinateDimension
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ellipsoidMesh
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_OriginSet",err,error,*999)

    CALL GeneratedMesh_AssertNotFinished(generatedMesh,err,error,*999)

    NULLIFY(coordinateSystem)
    CALL GeneratedMesh_CoordinateSystemGet(generatedMesh,coordinateSystem,err,error,*999)
    CALL CoordinateSystem_DimensionGet(coordinateSystem,coordinateDimension,err,error,*999)
    IF(SIZE(origin,1)<coordinateDimension) THEN
      localError="The specified origin size of "//TRIM(NumberToVString(SIZE(origin,1),"*",err,error))// &
        & " is invalid. The extent size must be >= coordinate system dimension of "// &
        & TRIM(NumberToVString(coordinateDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(generatedMesh%generatedType)
    CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
      NULLIFY(regularMesh)
      CALL GeneratedMesh_RegularMeshGet(generatedMesh,regularMesh,err,error,*999)
      IF(ALLOCATED(regularMesh%origin)) DEALLOCATE(regularMesh%origin)
      ALLOCATE(regularMesh%origin(coordinateDimension),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate origin.",err,error,*999)             
      regularMesh%origin(1:coordinateDimension)=origin(1:coordinateDimension)
    CASE(GENERATED_MESH_POLAR_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
      NULLIFY(cylinderMesh)
      CALL GeneratedMesh_CylinderMeshGet(generatedMesh,cylinderMesh,err,error,*999)
      IF(ALLOCATED(cylinderMesh%origin)) DEALLOCATE(cylinderMesh%origin)
      ALLOCATE(cylinderMesh%origin(coordinateDimension),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate origin.",err,error,*999)             
      cylinderMesh%origin(1:coordinateDimension)=origin(1:coordinateDimension)
    CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
      NULLIFY(ellipsoidMesh)
      CALL GeneratedMesh_EllipsoidMeshGet(generatedMesh,ellipsoidMesh,err,error,*999)
      IF(ALLOCATED(ellipsoidMesh%origin)) DEALLOCATE(ellipsoidMesh%origin)
      ALLOCATE(ellipsoidMesh%origin(coordinateDimension),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate origin.",err,error,*999)             
      ellipsoidMesh%origin(1:coordinateDimension)=origin(1:coordinateDimension)
    CASE DEFAULT
      localError="The generated mesh mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("GeneratedMesh_OriginSet")
    RETURN
999 ERRORSEXITS("GeneratedMesh_OriginSet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_OriginSet

  !
  !================================================================================================================================
  !
  
  !>Start to create the regular generated mesh type
  SUBROUTINE GeneratedMesh_RegularCreateFinish(generatedMesh,meshUserNumber,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(IN) :: meshUserNumber !<The user number for the mesh to generate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coordinateIdx,coordinateType,count,elementFactor,gridElementNumber,gridNumberOfElements,elementNumber, &
      & elementIdx1,elementIdx2,elementIdx3,localNodeIdx,localNodeIdx1,localNodeIdx2,localNodeIdx3,nodeNumber, &
      & numberOfElementsXi(3),totalNumberOfNodesXi(3),totalNumberOfNodes,numberOfCornerNodes, &
      & totalNumberOfElements,xiIdx,numberOfBases,basisIdx,basisNumberOfNodes(10)
    INTEGER(INTG), ALLOCATABLE :: elementNodes(:),elementNodesUserNumbers(:)
    TYPE(BasisType), POINTER :: basis
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh
    TYPE(InterfaceType), POINTER :: interface
    TYPE(MeshElementsType), POINTER :: meshElements
    TYPE(NodesType), POINTER :: nodes
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_RegularCreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated Mesh is not associated.",err,error,*999)
    NULLIFY(coordinateSystem)
    CALL GeneratedMesh_CoordinateSystemGet(generatedMesh,coordinateSystem,err,error,*999)
    CALL CoordinateSystem_TypeGet(coordinateSystem,coordinateType,err,error,*999)
    NULLIFY(regularMesh)
    CALL GeneratedMesh_RegularMeshGet(generatedMesh,regularMesh,err,error,*999)
    region=>generatedMesh%region
    INTERFACE=>generatedMesh%INTERFACE
      
    SELECT CASE(coordinateType)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      NULLIFY(basis)
      CALL GeneratedMeshRegular_BasisGet(regularMesh,1,basis,err,error,*999)
      SELECT CASE(basis%TYPE)
      CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE,BASIS_SIMPLEX_TYPE)
        IF(.NOT.ALL(basis%collapsedXi==BASIS_NOT_COLLAPSED))  &
          & CALL FlagError("Degenerate (collapsed) basis not implemented.",err,error,*999)
        !Determine the coordinate system and create the regular mesh for that system
        regularMesh%coordinateDimension=coordinateSystem%numberOfDimensions
        regularMesh%meshDimension=basis%numberOfXi
        IF(.NOT.ALLOCATED(regularMesh%origin)) THEN
          ALLOCATE(regularMesh%origin(regularMesh%coordinateDimension),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate origin.",err,error,*999)
          regularMesh%origin=0.0_DP
        ENDIF
        IF(.NOT.ALLOCATED(regularMesh%maximumExtent)) THEN
          ALLOCATE(regularMesh%maximumExtent(regularMesh%coordinateDimension),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate maximum extent.",err,error,*999)
          regularMesh%maximumExtent=1.0_DP
        ENDIF
        IF(.NOT.ALLOCATED(regularMesh%numberOfElementsXi)) THEN
          ALLOCATE(regularMesh%numberOfElementsXi(regularMesh%meshDimension),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate number of elements xi.",err,error,*999)
          regularMesh%numberOfElementsXi=1
        ENDIF
        IF(ALLOCATED(regularMesh%baseVectors)) THEN
!!TODO: Check base vectors
        ELSE
          !Calculate base vectors
          ALLOCATE(regularMesh%baseVectors(regularMesh%coordinateDimension,regularMesh%meshDimension),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate number of elements xi.",err,error,*999)
          regularMesh%baseVectors=0.0_DP
          IF(regularMesh%meshDimension==1) THEN
            !The base vector is just the extent vector
            regularMesh%baseVectors(:,1)=regularMesh%maximumExtent
          ELSE
            IF(regularMesh%meshDimension<regularMesh%coordinateDimension) THEN
              !Find the first number of mesh dimensions for which the extent is non-zero.
              count=0
              coordinateIdx=1
              DO xiIdx=1,regularMesh%meshDimension
                DO WHILE(ABS(regularMesh%maximumExtent(coordinateIdx))<=ZERO_TOLERANCE)
                  coordinateIdx=coordinateIdx+1
                ENDDO
                regularMesh%baseVectors(coordinateIdx,xiIdx)=regularMesh%maximumExtent(coordinateIdx)
                coordinateIdx=coordinateIdx+1
                count=count+1
              ENDDO !xiIdx
              IF(count/=regularMesh%meshDimension)  &
                & CALL FlagError("Invalid mesh extent. There number of non-zero components is < the mesh dimension.", &
                & err,error,*999)
            ELSE IF(regularMesh%meshDimension==regularMesh%coordinateDimension) THEN
              !The default base vectors are aligned with the coordinate vectors
              DO coordinateIdx=1,regularMesh%coordinateDimension
                regularMesh%baseVectors(coordinateIdx,coordinateIdx)=regularMesh%maximumExtent(coordinateIdx)
              ENDDO !coordinateIdx
            ELSE
              CALL FlagError("The mesh dimension is greater than the coordinate dimension.",err,error,*999)
            ENDIF
          ENDIF
        ENDIF
        !Calculate the sizes of a regular grid of elements with the appropriate number of basis nodes in each dimension of
        !the grid element
        totalNumberOfNodes=1
        gridNumberOfElements=1
        numberOfElementsXi=1
        numberOfBases=SIZE(regularMesh%bases)
        DO xiIdx=1,regularMesh%meshDimension
          !Set total number of nodes to corner nodes only
          totalNumberOfNodes=totalNumberOfNodes*(regularMesh%numberOfElementsXi(xiIdx)+1)
          numberOfElementsXi(xiIdx)=regularMesh%numberOfElementsXi(xiIdx)
          gridNumberOfElements=gridNumberOfElements*regularMesh%numberOfElementsXi(xiIdx)
        ENDDO !xiIdx
        numberOfCornerNodes=totalNumberOfNodes
        !Add extra nodes for each basis
        !Will end up with some duplicate nodes if bases have the same interpolation in one direction
        basisNumberOfNodes=0
        DO basisIdx=1,numberOfBases
          NULLIFY(basis)
          CALL GeneratedMeshRegular_BasisGet(regularMesh,basisIdx,basis,err,error,*999)
          basisNumberOfNodes(basisIdx)=1
          DO xiIdx=1,regularMesh%meshDimension
            basisNumberOfNodes(basisIdx)=basisNumberOfNodes(basisIdx)*((basis%numberOfNodesXiC(xiIdx)-1)* &
              & regularMesh%numberOfElementsXi(xiIdx)+1)
          ENDDO !xiIdx
          basisNumberOfNodes(basisIdx)=totalNumberOfNodes+basisNumberOfNodes(basisIdx)-numberOfCornerNodes
          ! basisNumberOfNodes=1
          ! DO xiIdx=1,regularMesh%meshDimension
          !   basisNumberOfNodes=basisNumberOfNodes*((basis%numberOfNodesXiC(xiIdx)-1)* &
          !       & regularMesh%numberOfElementsXi(xiIdx)+1)
          ! ENDDO !xiIdx
          ! totalNumberOfNodes=totalNumberOfNodes+basisNumberOfNodes-numberOfCornerNodes
        ENDDO
        totalNumberOfNodes=MAXVAL(basisNumberOfNodes)
        !Compute the element factor i.e., the number of sub elements each grid element will be split into.
        IF(basis%type==BASIS_LAGRANGE_HERMITE_TP_TYPE) THEN
          elementFactor=1
        ELSE
          SELECT CASE(regularMesh%meshDimension)
          CASE(1)
            elementFactor=1
          CASE(2)
            elementFactor=2
          CASE(3)
            elementFactor=6
          CASE DEFAULT
            localError="The mesh dimension dimension of "//TRIM(NumberToVString(regularMesh%meshDimension,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDIF
        totalNumberOfElements=elementFactor*gridNumberOfElements
        !Create the default node set
        NULLIFY(nodes)
        IF(ASSOCIATED(region)) THEN
          CALL Nodes_CreateStart(region,totalNumberOfNodes,nodes,err,error,*999)
        ELSE
          CALL Nodes_CreateStart(INTERFACE,totalNumberOfNodes,nodes,err,error,*999)
        ENDIF
        !Finish the nodes creation
        CALL Nodes_CreateFinish(nodes,err,error,*999)
        !Create the mesh
        IF(ASSOCIATED(region)) THEN
          CALL Mesh_CreateStart(meshUserNumber,region,regularMesh%meshDimension,generatedMesh%mesh,err,error,*999)
        ELSE
          CALL Mesh_CreateStart(meshUserNumber,INTERFACE,regularMesh%meshDimension,generatedMesh%mesh,err,error,*999)
        ENDIF
        !Set the number of mesh components
        CALL Mesh_NumberOfComponentsSet(generatedMesh%mesh,numberOfBases,err,error,*999)
        !Create the elements
        CALL Mesh_NumberOfElementsSet(generatedMesh%mesh,totalNumberOfElements,err,error,*999)
        DO basisIdx=1,numberOfBases
          NULLIFY(basis)
          CALL GeneratedMeshRegular_BasisGet(regularMesh,basisIdx,basis,err,error,*999)
          !Get number of nodes in each xi direction for this basis
          DO xiIdx=1,regularMesh%meshDimension
            totalNumberOfNodesXi(xiIdx)=(basis%numberOfNodesXiC(xiIdx)-1)*regularMesh%numberOfElementsXi(xiIdx)+1
          ENDDO !xiIdx
          NULLIFY(meshElements)
          CALL MeshElements_CreateStart(generatedMesh%mesh,basisIdx,basis,meshElements,err,error,*999)
          !Set the elements for the regular mesh
          IF (ALLOCATED(elementNodes)) DEALLOCATE(elementNodes)
          ALLOCATE(elementNodes(basis%numberOfNodes),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate element nodes.",err,error,*999)
          IF (ALLOCATED(elementNodesUserNumbers)) DEALLOCATE(elementNodesUserNumbers)
          ALLOCATE(elementNodesUserNumbers(basis%numberOfNodes),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate element nodes user numbers.",err,error,*999)
          !Step in the xi(3) direction
          DO elementIdx3=1,numberOfElementsXi(3)+1
            DO elementIdx2=1,numberOfElementsXi(2)+1
              DO elementIdx1=1,numberOfElementsXi(1)+1
                IF(basis%numberOfXi<3.OR.elementIdx3<=numberOfElementsXi(3)) THEN
                  IF(basis%numberOfXi<2.OR.elementIdx2<=numberOfElementsXi(2)) THEN
                    IF(elementIdx1<=numberOfElementsXi(1)) THEN
                      gridElementNumber=elementIdx1
                      nodeNumber=1+(elementIdx1-1)*(basis%numberOfNodesXiC(1)-1)
                      IF(basis%numberOfXi>1) THEN
                        gridElementNumber=gridElementNumber+(elementIdx2-1)*numberOfElementsXi(1)
                        nodeNumber=nodeNumber+(elementIdx2-1)*totalNumberOfNodesXi(1)*(basis%numberOfNodesXiC(2)-1)
                        IF(basis%numberOfXi>2) THEN
                          gridElementNumber=gridElementNumber+(elementIdx3-1)*numberOfElementsXi(1)*numberOfElementsXi(2)
                          nodeNumber=nodeNumber+(elementIdx3-1)*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)* &
                            & (basis%numberOfNodesXiC(3)-1)
                        ENDIF
                      ENDIF
                      IF(basis%TYPE==basis_LAGRANGE_HERMITE_TP_TYPE) THEN
                        !Lagrange Hermite TP elements
                        elementNumber=gridElementNumber
                        localNodeIdx=0
                        DO localNodeIdx1=1,basis%numberOfNodesXiC(1)
                          localNodeIdx=localNodeIdx+1
                          elementNodes(localNodeIdx)=nodeNumber+(localNodeIdx1-1)
                        ENDDO !localNodeIdx1
                        IF(basis%numberOfXi>1) THEN
                          DO localNodeIdx2=2,basis%numberOfNodesXiC(2)
                            DO localNodeIdx1=1,basis%numberOfNodesXiC(1)
                              localNodeIdx=localNodeIdx+1
                              elementNodes(localNodeIdx)=nodeNumber+(localNodeIdx1-1)+(localNodeIdx2-1)*totalNumberOfNodesXi(1)
                            ENDDO !localNodeIdx1
                          ENDDO !localNodeIdx2
                          IF(basis%numberOfXi>2) THEN
                            DO localNodeIdx3=2,basis%numberOfNodesXiC(3)
                              DO localNodeIdx2=1,basis%numberOfNodesXiC(2)
                                DO localNodeIdx1=1,basis%numberOfNodesXiC(1)
                                  localNodeIdx=localNodeIdx+1
                                  elementNodes(localNodeIdx)=nodeNumber+(localNodeIdx1-1)+ &
                                    & (localNodeIdx2-1)*totalNumberOfNodesXi(1)+ &
                                    & (localNodeIdx3-1)*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                                ENDDO !localNodeIdx1
                              ENDDO !localNodeIdx2
                            ENDDO !localNodeIdx3
                          ENDIF
                        ENDIF
                        CALL GeneratedMesh_RegularComponentNodesToUserNumbers(regularMesh%generatedMesh, &
                          & basisIdx,elementNodes,elementNodesUserNumbers,err,error,*999)
                        CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                      ELSE
                        !Simplex elements
                        SELECT CASE(basis%numberOfXi)
                        CASE(1)
                          !Line element
                          elementNumber=gridElementNumber
                          localNodeIdx=0
                          DO localNodeIdx1=1,basis%numberOfNodesXiC(1)
                            localNodeIdx=localNodeIdx+1
                            elementNodes(localNodeIdx)=nodeNumber+(localNodeIdx1-1)
                          ENDDO !localNodeIdx1
                          CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                            & elementNodesUserNumbers,err,error,*999)
                          CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                        CASE(2)
                          !Triangular element
                          !Break the grid square element into 2 triangles. The 2 triangles are
                          !Element 1: vertices {(0,0);(1,0);(1,1)}
                          !Element 2: vertices {(0,0);(1,1);(0,1)}
                          SELECT CASE(basis%interpolationOrder(1))
                          CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
                            !First sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+1
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+1
                            elementNodes(3)=nodeNumber+1+totalNumberOfNodesXi(1)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                            !Second sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+2
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+1+totalNumberOfNodesXi(1)
                            elementNodes(3)=nodeNumber+totalNumberOfNodesXi(1)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                          CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
                            !First sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+1
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+2
                            elementNodes(3)=nodeNumber+2+2*totalNumberOfNodesXi(1)
                            elementNodes(4)=nodeNumber+1
                            elementNodes(5)=nodeNumber+2+totalNumberOfNodesXi(1)
                            elementNodes(6)=nodeNumber+1+totalNumberOfNodesXi(1)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                            !Second sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+2
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+2+2*totalNumberOfNodesXi(1)
                            elementNodes(3)=nodeNumber+2*totalNumberOfNodesXi(1)
                            elementNodes(4)=nodeNumber+1+totalNumberOfNodesXi(1)
                            elementNodes(5)=nodeNumber+1+2*totalNumberOfNodesXi(1)
                            elementNodes(6)=nodeNumber+totalNumberOfNodesXi(1)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                          CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
                            !First sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+1
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+3
                            elementNodes(3)=nodeNumber+3+3*totalNumberOfNodesXi(1)
                            elementNodes(4)=nodeNumber+1
                            elementNodes(5)=nodeNumber+2
                            elementNodes(6)=nodeNumber+3+totalNumberOfNodesXi(1)
                            elementNodes(7)=nodeNumber+3+2*totalNumberOfNodesXi(1)
                            elementNodes(8)=nodeNumber+2+2*totalNumberOfNodesXi(1)
                            elementNodes(9)=nodeNumber+1+totalNumberOfNodesXi(1)
                            elementNodes(10)=nodeNumber+2+totalNumberOfNodesXi(1)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                            !Second sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+2
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+3+3*totalNumberOfNodesXi(1)
                            elementNodes(3)=nodeNumber+3*totalNumberOfNodesXi(1)
                            elementNodes(4)=nodeNumber+2+2*totalNumberOfNodesXi(1)
                            elementNodes(5)=nodeNumber+1+totalNumberOfNodesXi(1)
                            elementNodes(6)=nodeNumber+2+3*totalNumberOfNodesXi(1)
                            elementNodes(7)=nodeNumber+1+3*totalNumberOfNodesXi(1)
                            elementNodes(8)=nodeNumber+totalNumberOfNodesXi(1)
                            elementNodes(9)=nodeNumber+2*totalNumberOfNodesXi(1)
                            elementNodes(10)=nodeNumber+1+2*totalNumberOfNodesXi(1)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                          CASE DEFAULT
                            localError="The simplex basis interpolation order of "// &
                              & TRIM(NumberToVString(basis%interpolationOrder(1),"*",err,error))// &
                              & " is invalid."
                            CALL FlagError(localError,err,error,*999)
                          END SELECT
                        CASE(3)
                          !Tetrahedra element
                          !Break the grid cube element into 6 tetrahedra (so that we have a break down the main diagonal of the
                          !cube in order to allow for the middle node in quadratics to be included). The 6 tetrahedra are
                          !Element 1: vertices {(0,0,0);(1,0,0);(1,1,0);(1,1,1)}
                          !Element 2: vertices {(0,0,0);(1,1,0);(0,1,0);(1,1,1)}
                          !Element 3: vertices {(0,0,0);(1,0,1);(1,0,0);(1,1,1)}
                          !Element 4: vertices {(0,0,0);(0,0,1);(1,0,1);(1,1,1)}
                          !Element 5: vertices {(0,0,0);(0,1,0);(0,1,1);(1,1,1)}
                          !Element 6: vertices {(0,0,0);(0,1,1);(0,0,1);(1,1,1)}
                          SELECT CASE(basis%interpolationOrder(1))
                          CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
                            !First sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+1
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+1
                            elementNodes(3)=nodeNumber+1+totalNumberOfNodesXi(1)
                            elementNodes(4)=nodeNumber+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                            !Second sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+2
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+1+totalNumberOfNodesXi(1)
                            elementNodes(3)=nodeNumber+totalNumberOfNodesXi(1)
                            elementNodes(4)=nodeNumber+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                            !Third sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+3
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+1+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(3)=nodeNumber+1
                            elementNodes(4)=nodeNumber+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                            !Fourth sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+4
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(3)=nodeNumber+1+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(4)=nodeNumber+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                            !Fifth sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+5
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+totalNumberOfNodesXi(1)
                            elementNodes(3)=nodeNumber+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(4)=nodeNumber+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                            !Sixth sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+6
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(3)=nodeNumber+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(4)=nodeNumber+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                          CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
                            !First sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+1
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+2
                            elementNodes(3)=nodeNumber+2+2*totalNumberOfNodesXi(1)
                            elementNodes(4)=nodeNumber+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(5)=nodeNumber+1
                            elementNodes(6)=nodeNumber+1+totalNumberOfNodesXi(1)
                            elementNodes(7)=nodeNumber+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(8)=nodeNumber+2+totalNumberOfNodesXi(1)
                            elementNodes(9)=nodeNumber+2+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(10)=nodeNumber+2+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                            !Second sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+2
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+2+2*totalNumberOfNodesXi(1)
                            elementNodes(3)=nodeNumber+2*totalNumberOfNodesXi(1)
                            elementNodes(4)=nodeNumber+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(5)=nodeNumber+1+totalNumberOfNodesXi(1)
                            elementNodes(6)=nodeNumber+totalNumberOfNodesXi(1)
                            elementNodes(7)=nodeNumber+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(8)=nodeNumber+1+2*totalNumberOfNodesXi(1)
                            elementNodes(9)=nodeNumber+1+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(10)=nodeNumber+2+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                            !Third sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+3
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+2+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(3)=nodeNumber+2
                            elementNodes(4)=nodeNumber+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(5)=nodeNumber+1+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(6)=nodeNumber+1
                            elementNodes(7)=nodeNumber+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(8)=nodeNumber+2+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(9)=nodeNumber+2+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(10)=nodeNumber+2+totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                            !Fourth sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+4
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(3)=nodeNumber+2+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(4)=nodeNumber+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(5)=nodeNumber+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(6)=nodeNumber+1+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(7)=nodeNumber+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(8)=nodeNumber+1+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(9)=nodeNumber+2+totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(10)=nodeNumber+1+totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                            !Fifth sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+5
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+2*totalNumberOfNodesXi(1)
                            elementNodes(3)=nodeNumber+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(4)=nodeNumber+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(5)=nodeNumber+totalNumberOfNodesXi(1)
                            elementNodes(6)=nodeNumber+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(7)=nodeNumber+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(8)=nodeNumber+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(9)=nodeNumber+1+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(10)=nodeNumber+1+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                            !Sixth sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+6
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(3)=nodeNumber+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(4)=nodeNumber+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(5)=nodeNumber+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(6)=nodeNumber+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(7)=nodeNumber+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(8)=nodeNumber+totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(9)=nodeNumber+1+totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(10)=nodeNumber+1+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                          CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
                            !First sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+1
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+3
                            elementNodes(3)=nodeNumber+3+3*totalNumberOfNodesXi(1)
                            elementNodes(4)=nodeNumber+3+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(5)=nodeNumber+1
                            elementNodes(6)=nodeNumber+2
                            elementNodes(7)=nodeNumber+1+totalNumberOfNodesXi(1)
                            elementNodes(8)=nodeNumber+2+2*totalNumberOfNodesXi(1)
                            elementNodes(9)=nodeNumber+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(10)=nodeNumber+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            elementNodes(11)=nodeNumber+3+totalNumberOfNodesXi(1)
                            elementNodes(12)=nodeNumber+3+2*totalNumberOfNodesXi(1)
                            elementNodes(13)=nodeNumber+3+3*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(14)=nodeNumber+3+3*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            elementNodes(15)=nodeNumber+3+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(16)=nodeNumber+3+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            elementNodes(17)=nodeNumber+2+totalNumberOfNodesXi(1)
                            elementNodes(18)=nodeNumber+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(19)=nodeNumber+2+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(20)=nodeNumber+3+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                            !Second sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+2
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+3+3*totalNumberOfNodesXi(1)
                            elementNodes(3)=nodeNumber+3*totalNumberOfNodesXi(1)
                            elementNodes(4)=nodeNumber+3+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(5)=nodeNumber+1+totalNumberOfNodesXi(1)
                            elementNodes(6)=nodeNumber+2+2*totalNumberOfNodesXi(1)
                            elementNodes(7)=nodeNumber+totalNumberOfNodesXi(1)
                            elementNodes(8)=nodeNumber+2*totalNumberOfNodesXi(1)
                            elementNodes(9)=nodeNumber+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(10)=nodeNumber+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            elementNodes(11)=nodeNumber+2+3*totalNumberOfNodesXi(1)
                            elementNodes(12)=nodeNumber+1+3*totalNumberOfNodesXi(1)
                            elementNodes(13)=nodeNumber+1+3*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(14)=nodeNumber+2+3*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            elementNodes(15)=nodeNumber+3+3*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(16)=nodeNumber+3+3*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            elementNodes(17)=nodeNumber+1+2*totalNumberOfNodesXi(1)
                            elementNodes(18)=nodeNumber+2+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(19)=nodeNumber+1+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(20)=nodeNumber+2+3*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                            !Third sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+3
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+3+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(3)=nodeNumber+3
                            elementNodes(4)=nodeNumber+3+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(5)=nodeNumber+1+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(6)=nodeNumber+2+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(7)=nodeNumber+1
                            elementNodes(8)=nodeNumber+2
                            elementNodes(9)=nodeNumber+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(10)=nodeNumber+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            elementNodes(11)=nodeNumber+3+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(12)=nodeNumber+3+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(13)=nodeNumber+3+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(14)=nodeNumber+3+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            elementNodes(15)=nodeNumber+3+totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(16)=nodeNumber+3+2*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            elementNodes(17)=nodeNumber+2+2*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            elementNodes(18)=nodeNumber+2+totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(19)=nodeNumber+2+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(20)=nodeNumber+3+totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                            !Fourth sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+4
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(3)=nodeNumber+3+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(4)=nodeNumber+3+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(5)=nodeNumber+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(6)=nodeNumber+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(7)=nodeNumber+1+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(8)=nodeNumber+2+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(9)=nodeNumber+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(10)=nodeNumber+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            elementNodes(11)=nodeNumber+1+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(12)=nodeNumber+2+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(13)=nodeNumber+3+totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(14)=nodeNumber+3+2*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            elementNodes(15)=nodeNumber+1+totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(16)=nodeNumber+2+2*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            elementNodes(17)=nodeNumber+1+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(18)=nodeNumber+1+totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(19)=nodeNumber+2+totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(20)=nodeNumber+2+totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                            !Fifth sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+5
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+3*totalNumberOfNodesXi(1)
                            elementNodes(3)=nodeNumber+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(4)=nodeNumber+3+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(5)=nodeNumber+totalNumberOfNodesXi(1)
                            elementNodes(6)=nodeNumber+2*totalNumberOfNodesXi(1)
                            elementNodes(7)=nodeNumber+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(8)=nodeNumber+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(9)=nodeNumber+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(10)=nodeNumber+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            elementNodes(11)=nodeNumber+3*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(12)=nodeNumber+3*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(13)=nodeNumber+1+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            elementNodes(14)=nodeNumber+2+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            elementNodes(15)=nodeNumber+1+3*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(16)=nodeNumber+2+3*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            elementNodes(17)=nodeNumber+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(18)=nodeNumber+1+2*totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(19)=nodeNumber+1+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            elementNodes(20)=nodeNumber+1+3*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                            !Sixth sub-element
                            elementNumber=(gridElementNumber-1)*elementFactor+6
                            elementNodes(1)=nodeNumber
                            elementNodes(2)=nodeNumber+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(3)=nodeNumber+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(4)=nodeNumber+3+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(5)=nodeNumber+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(6)=nodeNumber+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(7)=nodeNumber+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(8)=nodeNumber+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(9)=nodeNumber+1+totalNumberOfNodesXi(1)+totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(10)=nodeNumber+2+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            elementNodes(11)=nodeNumber+2*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(12)=nodeNumber+totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(13)=nodeNumber+1+totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(14)=nodeNumber+2+2*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            elementNodes(15)=nodeNumber+1+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            elementNodes(16)=nodeNumber+2+3*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            elementNodes(17)=nodeNumber+totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(18)=nodeNumber+1+2*totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            elementNodes(19)=nodeNumber+1+totalNumberOfNodesXi(1)+2*totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
                            elementNodes(20)=nodeNumber+1+2*totalNumberOfNodesXi(1)+3*totalNumberOfNodesXi(1)* &
                              & totalNumberOfNodesXi(2)
                            CALL GeneratedMeh_ComponentNodesToUserNumbers(regularMesh%generatedMesh,basisIdx,elementNodes, &
                              & elementNodesUserNumbers,err,error,*999)
                            CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
                          CASE DEFAULT
                            localError="The simplex basis interpolation order of "// &
                              & TRIM(NumberToVString(basis%interpolationOrder(1),"*",err,error))// &
                              & " is invalid."
                            CALL FlagError(localError,err,error,*999)
                          END SELECT
                        CASE DEFAULT
                          localError="The simplex number of xi directions of "// &
                            & TRIM(NumberToVString(basis%numberOfXi,"*",err,error))// &
                            & " is invalid."
                          CALL FlagError(localError,err,error,*999)
                        END SELECT
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO !elementIdx1
            ENDDO !elementIdx2
          ENDDO !elementIdx3
          CALL MeshElements_CreateFinish(meshElements,err,error,*999)
        ENDDO !basisIdx
        !Finish the mesh
        CALL Mesh_CreateFinish(generatedMesh%mesh,err,error,*999)
      CASE DEFAULT
        CALL FlagError("Basis type is either invalid or not implemented.",err,error,*999)
      END SELECT
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The coordinate system type of "//TRIM(NumberToVString(coordinateSystem%type,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    IF(ALLOCATED(elementNodes)) DEALLOCATE(elementNodes)
    IF(ALLOCATED(elementNodesUserNumbers)) DEALLOCATE(elementNodesUserNumbers)

    EXITS("GeneratedMesh_RegularCreateFinish")
    RETURN
    ! TODO invalidate other associations
999 IF(ALLOCATED(elementNodes)) DEALLOCATE(elementNodes)
    IF(ALLOCATED(elementNodesUserNumbers)) DEALLOCATE(elementNodesUserNumbers)
    ERRORSEXITS("GeneratedMesh_RegularCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_RegularCreateFinish

  !
  !================================================================================================================================
  !

  !>Start to create the regular generated mesh type
  SUBROUTINE GeneratedMesh_EllipsoidCreateFinish(generatedMesh,meshUserNumber,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(IN) :: meshUserNumber !<The user number for the mesh to generate.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: basisNumberOfNodes,cornerNumberOfNodes,elementIdx1,elementIdx2,elementIdx3,elementNumber, &
      & from1,from2,from3,localNodeIdx,localNodeIdx1,localNodeIdx2,localNodeIdx3,numberOfDimensions,totalNumberOfElements, &
      & totalNumberOfNodes,basesIdx
    INTEGER(INTG), ALLOCATABLE :: apexElementNodes(:),apexElementNodesUserNumbers(:),cornerNodes(:,:,:),elementIndices(:,:,:), &
      & nodeIndices(:,:,:),numberOfElementsXi(:),wallElementNodes(:),wallElementNodesUserNumbers(:)
    REAL(DP) :: delta(3),deltai(3)
    TYPE(BasisType), POINTER :: basis1,basis2
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ellipsoidMesh
    TYPE(MeshElementsType), POINTER :: meshElements
    TYPE(NodesType), POINTER :: nodes
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_EllipsoidCreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated Mesh is not associated.",err,error,*999)

    NULLIFY(ellipsoidMesh)
    CALL GeneratedMesh_EllipsoidMeshGet(generatedMesh,ellipsoidMesh,err,error,*999)
    NULLIFY(region)
    CALL GeneratedMesh_RegionGet(generatedMesh,region,err,error,*999)
    NULLIFY(coordinateSystem)
    CALL Region_CoordinateSystemGet(region,coordinateSystem,err,error,*999)
    SELECT CASE(coordinateSystem%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      !Determine the coordinate system and create the regular mesh for that system
      ellipsoidMesh%meshDimension=coordinateSystem%numberOfDimensions
      numberOfDimensions=ellipsoidMesh%meshDimension
      IF(numberOfDimensions/=3) CALL FlagError("Ellipsoid mesh requires a 3 dimensional coordinate system.",err,error,*999)! hard-coded for 3D only
      IF(.NOT.ALLOCATED(ellipsoidMesh%origin)) THEN
        ALLOCATE(ellipsoidMesh%origin(numberOfDimensions),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate origin.",err,error,*999)
        ellipsoidMesh%origin=0.0_DP
      ENDIF
      IF(SIZE(ellipsoidMesh%origin)/=ellipsoidMesh%meshDimension) &
        & CALL FlagError("The number of dimensions of the given regular mesh does not match the size of the origin.",err,error,*999)
      IF(SIZE(ellipsoidMesh%ellipsoidExtent)/=4) THEN
        localError="For an ellipsoid mesh the following measures need to be given: longAxis, shortAxis, wallThickness "// &
          & "and cutoffAngle."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(.NOT.ALLOCATED(ellipsoidMesh%bases)) CALL FlagError("Bases is not allocated.",err,error,*999)
      IF(MOD(SIZE(ellipsoidMesh%bases),2)/=0) &
        & CALL FlagError("An ellipsoid mesh requires a collapsed basis for each basis, so there must be n*2 bases.",err,error,*999)
      ALLOCATE(numberOfElementsXi(SIZE(ellipsoidMesh%numberOfElementsXi)),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate number of elements xi.",err,error,*999)
      numberOfElementsXi(1:SIZE(ellipsoidMesh%numberOfElementsXi))= &
        & ellipsoidMesh%numberOfElementsXi(1:SIZE(ellipsoidMesh%numberOfElementsXi))
      !Calculate total number of nodes from all bases and start mesh
      cornerNumberOfNodes=numberOfElementsXi(1)*(numberOfElementsXi(2)+1)*(numberOfElementsXi(3)+1)- &
        & (numberOfElementsXi(1)-1)*(numberOfElementsXi(3)+1)
      totalNumberOfNodes=cornerNumberOfNodes
      DO basesIdx=1,SIZE(ellipsoidMesh%bases),2
        NULLIFY(basis1)
        CALL GeneratedMeshEllipsoid_BasisGet(ellipsoidMesh,basesIdx,basis1,err,error,*999)
        basisNumberOfNodes=numberOfElementsXi(1)*(basis1%numberOfNodesXiC(1)-1)* &
          & (numberOfElementsXi(2)*(basis1%numberOfNodesXiC(2)-1)+1)* &
          & (numberOfElementsXi(3)*(basis1%numberOfNodesXiC(3)-1)+1)- &
          & (numberOfElementsXi(1)*(basis1%numberOfNodesXiC(1)-1)-1)* &
          & (numberOfElementsXi(3)*(basis1%numberOfNodesXiC(3)-1)+1)
        totalNumberOfNodes=totalNumberOfNodes+basisNumberOfNodes-cornerNumberOfNodes
      ENDDO !basesIdx
      NULLIFY(nodes)
      CALL Nodes_CreateStart(region,totalNumberOfNodes,nodes,err,error,*999)
      !Finish the nodes creation
      CALL Nodes_CreateFinish(nodes,err,error,*999)
      !Create the mesh
      CALL Mesh_CreateStart(meshUserNumber,generatedMesh%region,SIZE(numberOfElementsXi,1),generatedMesh%mesh,err,error,*999)
      !Create the elements
      CALL Mesh_NumberOfComponentsSet(generatedMesh%mesh,SIZE(ellipsoidMesh%bases)/2,err,error,*999)
      DO basesIdx=1,SIZE(ellipsoidMesh%bases),2
        NULLIFY(basis1)
        CALL GeneratedMeshEllipsoid_BasisGet(ellipsoidMesh,basesIdx,basis1,err,error,*999)
        NULLIFY(basis2)
        CALL GeneratedMeshEllipsoid_BasisGet(ellipsoidMesh,basesIdx+1,basis2,err,error,*999)
        IF((basis1%numberOfCollapsedXi/=0).OR.(basis2%numberOfCollapsedXi<=0)) THEN
          !test for collapsed nodes and force non collapsed to wall elements and collapsed to apex elements
          CALL FlagError("For each basis, one non collapsed version (basis1) and one collapsed version (basis2) is needed.", &
            & err,error,*999)
        ENDIF
        SELECT CASE(basis1%TYPE)
          !should also test for basis2
        CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
          IF(basis1%numberOfXi/=SIZE(numberOfElementsXi,1).OR.basis2%numberOfXi/=SIZE(numberOfElementsXi,1)) THEN
            localError="The number of xi directions of the given basis does not match the size of the number of elements "// &
              & "for the mesh."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          IF(.NOT.ALL(numberOfElementsXi>0)) CALL FlagError("Must have 1 or more elements in all directions.",err,error,*999)
          IF(numberOfElementsXi(1)<3) CALL FlagError("Need >2 elements around the circumferential direction.",err,error,*999)
          !IF(.NOT.ALL(basis%collapsedXi==basis_NOT_COLLAPSED))  &
          !    & CALL FlagError("Degenerate (collapsed) basis not implemented.",err,error,*999)
          IF(ALLOCATED(nodeIndices)) DEALLOCATE(nodeIndices)
          IF(ALLOCATED(elementIndices)) DEALLOCATE(elementIndices)
          IF(ALLOCATED(cornerNodes)) DEALLOCATE(cornerNodes)
          CALL GeneratedMesh_EllipsoidBuildNodeIndices(numberOfElementsXi,basis1% &
            & numberOfNodesXiC, ellipsoidMesh%ellipsoidExtent, totalNumberOfNodes, &
            & totalNumberOfElements, nodeIndices,cornerNodes,elementIndices,delta,deltai,err,error,*999)
          IF(basesIdx==1) CALL Mesh_NumberOfElementsSet(generatedMesh%mesh,totalNumberOfElements,err,error,*999)
          !Create the default node set
          !TODO we finish create after the nodes are initialised?
          NULLIFY(meshElements)
          CALL MeshElements_CreateStart(generatedMesh%mesh,basesIdx/2+1,basis1,meshElements,err,error,*999)
          !Set the elements for the ellipsoid mesh
          IF(ALLOCATED(wallElementNodes)) DEALLOCATE(wallElementNodes)
          IF(ALLOCATED(apexElementNodes)) DEALLOCATE(apexElementNodes)
          IF(ALLOCATED(wallElementNodesUserNumbers)) DEALLOCATE(wallElementNodesUserNumbers)
          IF(ALLOCATED(apexElementNodesUserNumbers)) DEALLOCATE(apexElementNodesUserNumbers)
          ALLOCATE(wallElementNodes(basis1%numberOfNodes),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate wall element nodes.",err,error,*999)
          ALLOCATE(apexElementNodes(basis2%numberOfNodes),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate apex element nodes.",err,error,*999)
          ALLOCATE(wallElementNodesUserNumbers(basis1%numberOfNodes),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate wall element nodes.",err,error,*999)
          ALLOCATE(apexElementNodesUserNumbers(basis2%numberOfNodes),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate apex element nodes.",err,error,*999)
          ! calculate element topology (nodes per each element)
          ! the idea is to translate given (r,theta,z) to nodeIndices equivalents, which include interior nodes
          elementNumber=0
          localNodeIdx=0
          !fromJ=global J direction counting number of first node in element in J direction
          DO elementIdx3=1,numberOfElementsXi(3)
            from3=NINT(delta(3)*(elementIdx3-1)/deltai(3)+1)
            elementIdx2=1
            from2=NINT(delta(2)*(elementIdx2-1)/deltai(2)+1)
            !apex elements
            DO elementIdx1=1,numberOfElementsXi(1)
              from1=NINT(delta(1)*(elementIdx1-1)/deltai(1)+1)
              localNodeIdx=0
              ! number of nodes in an element is dependent on basis used
              DO localNodeIdx3=from3,from3+basis2%numberOfNodesXiC(3)-1
                localNodeIdx2=1
                localNodeIdx1=1
                !central axis nodes
                localNodeIdx=localNodeIdx+1
                apexElementNodes(localNodeIdx)=nodeIndices(localNodeIdx1,localNodeIdx2,localNodeIdx3)
                DO localNodeIdx2=from2+1,from2+basis2%numberOfNodesXiC(2)-1
                  DO localNodeIdx1=from1,from1+basis2%numberOfNodesXiC(1)-1
                    localNodeIdx=localNodeIdx+1
                    ! circumferential loop-around
                    IF(localNodeIdx1>SIZE(nodeIndices,1)) THEN
                      apexElementNodes(localNodeIdx)=nodeIndices(1,localNodeIdx2,localNodeIdx3)
                    ELSE
                      apexElementNodes(localNodeIdx)=nodeIndices(localNodeIdx1,localNodeIdx2,localNodeIdx3)
                    ENDIF
                  ENDDO !localNodeIdx1
                ENDDO !localNodeIdx2
              ENDDO !localNodeIdx3
              elementNumber=elementNumber+1
              CALL MeshElements_ElementBasisSet(meshElements,elementNumber,basis2,err,error,*999)
              CALL GeneratedMeh_ComponentNodesToUserNumbers(ellipsoidMesh%generatedMesh,basesIdx,apexElementNodes, &
                & apexElementNodesUserNumbers,err,error,*999)
              CALL MeshElements_ElementNodesSet(meshElements,elementNumber,apexElementNodesUserNumbers,err,error,*999)
            ENDDO ! elementIdx1
            !wall elements
            DO elementIdx2=2,numberOfElementsXi(2)
              from2=NINT(delta(2)*(elementIdx2-1)/deltai(2)+1)
              DO elementIdx1=1,numberOfElementsXi(1)
                from1=NINT(delta(1)*(elementIdx1-1)/deltai(1)+1)
                localNodeIdx=0
                ! number of nodes in an element is dependent on basis used
                DO localNodeIdx3=from3,from3+basis1%numberOfNodesXiC(3)-1
                  DO localNodeIdx2=from2,from2+basis1%numberOfNodesXiC(2)-1
                    DO localNodeIdx1=from1,from1+basis1%numberOfNodesXiC(1)-1
                      localNodeIdx=localNodeIdx+1
                      ! circumferential loop-around
                      IF(localNodeIdx1>SIZE(nodeIndices,1)) THEN
                        wallElementNodes(localNodeIdx)=nodeIndices(1,localNodeIdx2,localNodeIdx3)
                      ELSE
                        wallElementNodes(localNodeIdx)=nodeIndices(localNodeIdx1,localNodeIdx2,localNodeIdx3)
                      ENDIF
                    ENDDO ! localNodeIdx1
                  ENDDO ! localNodeIdx2
                ENDDO ! localNodeIdx3
                elementNumber=elementNumber+1
                CALL GeneratedMeh_ComponentNodesToUserNumbers(ellipsoidMesh%generatedMesh,basesIdx,wallElementNodes, &
                  & wallElementNodesUserNumbers,err,error,*999)
                CALL MeshElements_ElementNodesSet(meshElements,elementNumber,wallElementNodesUserNumbers,err,error,*999)
              ENDDO ! elementIdx1
            ENDDO ! elementIdx2
          ENDDO ! elementIdx3
          CALL MeshElements_CreateFinish(meshElements,err,error,*999)
        CASE(BASIS_SIMPLEX_TYPE)
          CALL FlagError("Ellipsoid meshes with simplex basis types is not implemented.",err,error,*999)
        CASE DEFAULT
          CALL FlagError("Basis type is either invalid or not implemented.",err,error,*999)
        END SELECT
      ENDDO
      !Finish the mesh
      CALL Mesh_CreateFinish(generatedMesh%mesh,err,error,*999)
    CASE DEFAULT
      CALL FlagError("Coordinate type is either invalid or not implemented.",err,error,*999)
    END SELECT
    
    EXITS("GeneratedMesh_EllipsoidCreateFinish")
    RETURN
    ! TODO invalidate other associations
999 IF(ALLOCATED(nodeIndices)) DEALLOCATE(nodeIndices)
    IF(ALLOCATED(elementIndices)) DEALLOCATE(elementIndices)
    IF(ALLOCATED(cornerNodes)) DEALLOCATE(cornerNodes)
    IF(ALLOCATED(numberOfElementsXi)) DEALLOCATE(numberOfElementsXi)
    IF(ALLOCATED(wallElementNodes)) DEALLOCATE(wallElementNodes)
    IF(ALLOCATED(apexElementNodes)) DEALLOCATE(apexElementNodes)
    ERRORSEXITS("GeneratedMesh_EllipsoidCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_EllipsoidCreateFinish
  !
  !================================================================================================================================
  !

  !>Start to create the regular generated mesh type
  SUBROUTINE GeneratedMesh_CylinderCreateFinish(generatedMesh,meshUserNumber,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(IN) :: meshUserNumber !<The user number for the mesh to generate.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: basisIdx,basisNumberOfNodes,cornerNumberOfNodes, elementIdx1,elementIdx2,elementIdx3,elementNumber, &
      & from1,from2,from3,localNodeIdx,localNodeIdx1,localNodeIdx2,localNodeIdx3,numberOfDimensions,totalNumberOfNodes, &
      & totalNumberOfElements
    INTEGER(INTG), ALLOCATABLE :: elementNodes(:),elementNodesUserNumbers(:),nodeIndices(:,:,:),elementIndices(:,:,:)
    INTEGER(INTG), ALLOCATABLE :: numberOfElementsXi(:)
    REAL(DP) :: delta(3),deltai(3)
    TYPE(BasisType), POINTER :: basis
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh
    TYPE(NodesType), POINTER :: nodes
    TYPE(MeshElementsType), POINTER :: meshElements
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_CylinderCreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated Mesh is not associated.",err,error,*999)

    NULLIFY(cylinderMesh)
    CALL GeneratedMesh_CylinderMeshGet(generatedMesh,cylinderMesh,err,error,*999)
    NULLIFY(region)
    CALL GeneratedMesh_RegionGet(generatedMesh,region,err,error,*999)
    NULLIFY(coordinateSystem)
    CALL Region_CoordinateSystemGet(region,coordinateSystem,err,error,*999)
    !TODO is regular type only for COORDINATE_RECTANGULAR_CARTESIAN_TYPE?
    !If that, should we use IF rather than select?
    SELECT CASE(region%coordinateSystem%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      !Determine the coordinate system and create the regular mesh for that system
      cylinderMesh%meshDimension=coordinateSystem%numberOfDimensions
      numberOfDimensions=cylinderMesh%meshDimension
      IF(numberOfDimensions/=3) CALL FlagError("Cylinder mesh requires a 3 dimensional coordinate system.",err,error,*999) ! hard-coded for 3D only
      IF(.NOT.ALLOCATED(cylinderMesh%origin)) THEN
        ALLOCATE(cylinderMesh%origin(numberOfDimensions),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate origin.",err,error,*999)
        cylinderMesh%origin=0.0_DP
      ENDIF
      IF(SIZE(cylinderMesh%origin)/=cylinderMesh%meshDimension) &
        & CALL FlagError("The number of dimensions of the given cylinder  mesh does not match the size of the origin.", &
        & err,error,*999)
      IF(SIZE(cylinderMesh%cylinderExtent)/=cylinderMesh%meshDimension) &
        & CALL FlagError("The number of dimensions of the given cylinder mesh does not match the size of the cylinder extent.", &
        & err,error,*999)
      IF(.NOT.ALLOCATED(cylinderMesh%bases)) CALL FlagError("Bases are not allocated.",err,error,*999)
      ALLOCATE(numberOfElementsXi(SIZE(cylinderMesh%numberOfElementsXi)),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate number of elements xi.",err,error,*999)
      numberOfElementsXi(1:SIZE(cylinderMesh%numberOfElementsXi))= &
        & cylinderMesh%numberOfElementsXi(1:SIZE(cylinderMesh%numberOfElementsXi))
      CALL Mesh_CreateStart(meshUserNumber,generatedMesh%region,SIZE(numberOfElementsXi,1),generatedMesh%mesh,err,error,*999)
      CALL Mesh_NumberOfComponentsSet(generatedMesh%mesh,SIZE(cylinderMesh%bases),err,error,*999)
      !Calculate number of nodes
      cornerNumberOfNodes=(numberOfElementsXi(3)+1)*numberOfElementsXi(2)*(numberOfElementsXi(1)+1)
      totalNumberOfNodes=cornerNumberOfNodes
      DO basisIdx=1,SIZE(cylinderMesh%bases)
        NULLIFY(basis)
        CALL GeneratedMeshCylinder_BasisGet(cylinderMesh,basisIdx,basis,err,error,*999)
        basisNumberOfNodes=((basis%numberOfNodesXiC(3)-1)*numberOfElementsXi(3)+1)* &
          & ((basis%numberOfNodesXiC(2)-1)*numberOfElementsXi(2))* &
          & ((basis%numberOfNodesXiC(1)-1)*numberOfElementsXi(1)+1)
        totalNumberOfNodes=totalNumberOfNodes+basisNumberOfNodes-cornerNumberOfNodes
      ENDDO !basisIdx
      NULLIFY(nodes)
      CALL Nodes_CreateStart(region,totalNumberOfNodes,nodes,err,error,*999)
      !Finish the nodes creation
      CALL Nodes_CreateFinish(nodes,err,error,*999)
      !Set the total number of elements
      totalNumberOfElements=numberOfElementsXi(1)*numberOfElementsXi(2)*numberOfElementsXi(3)
      CALL Mesh_NumberOfElementsSet(generatedMesh%mesh,totalNumberOfElements,err,error,*999)
      DO basisIdx=1,SIZE(cylinderMesh%bases)
        NULLIFY(basis)
        CALL GeneratedMeshCylinder_BasisGet(cylinderMesh,basisIdx,basis,err,error,*999)
        SELECT CASE(basis%TYPE)
        CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
          IF(basis%numberOfXi/=SIZE(numberOfElementsXi,1)) THEN
            localError="The number of xi directions of the given basis does not match the size of the number of elements "// &
              & "for the mesh."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          IF(.NOT.ALL(numberOfElementsXi>0)) CALL FlagError("Must have 1 or more elements in all directions.",err,error,*999)
          IF(numberOfElementsXi(2)<3) CALL FlagError("Need >2 elements around the circumferential direction.",err,error,*999)
          IF(.NOT.ALL(basis%collapsedXi==BASIS_NOT_COLLAPSED))  &
            & CALL FlagError("Degenerate (collapsed) basis not implemented.",err,error,*999)
          !Calculate nodes and element sizes
          IF(ALLOCATED(nodeIndices)) DEALLOCATE(nodeIndices)
          IF(ALLOCATED(elementIndices)) DEALLOCATE(elementIndices)
          CALL GeneratedMesh_CylinderBuildNodeIndices(numberOfElementsXi,basis%numberOfNodesXiC, &
            & cylinderMesh%cylinderExtent, totalNumberOfNodes,totalNumberOfElements, &
            & nodeIndices,elementIndices,delta,deltai,err,error,*999)
          !Set the elements for the cylinder mesh
          IF(ALLOCATED(elementNodes)) DEALLOCATE(elementNodes)
          IF(ALLOCATED(elementNodesUserNumbers)) DEALLOCATE(elementNodesUserNumbers)
          ALLOCATE(elementNodesUserNumbers(basis%numberOfNodes),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate element nodes user numbers.",err,error,*999)
          ALLOCATE(elementNodes(basis%numberOfNodes),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate element nodes.",err,error,*999)
          !Create the elements
          NULLIFY(meshElements)
          CALL MeshElements_CreateStart(generatedMesh%mesh,basisIdx,basis,meshElements,err,error,*999)
          ! calculate element topology (nodes per each element)
          ! the idea is to translate given (r,theta,z) to nodeIndices equivalents, which include interior nodes
          elementNumber=0
          DO elementIdx3=1,numberOfElementsXi(3)
            from3=NINT(delta(3)*(elementIdx3-1)/deltai(3)+1)
            DO elementIdx2=1,numberOfElementsXi(2)
              from2=NINT(delta(2)*(elementIdx2-1)/deltai(2)+1)
              DO elementIdx1=1,numberOfElementsXi(1)
                from1=NINT(delta(1)*(elementIdx1-1)/deltai(1)+1)
                localNodeIdx=0
                ! number of nodes in an element is dependent on basis used
                DO localNodeIdx3=from3,from3+basis%numberOfNodesXiC(3)-1
                  DO localNodeIdx2=from2,from2+basis%numberOfNodesXiC(2)-1
                    DO localNodeIdx1=from1,from1+basis%numberOfNodesXiC(1)-1
                      localNodeIdx=localNodeIdx+1
                      ! compensate for circumferential loop-around
                      IF(localNodeIdx2>SIZE(nodeIndices,2)) THEN
                        ! DEBUG: little check here
                        IF(localNodeIdx2>SIZE(nodeIndices,2)+1) CALL FlagError("nodeIndices needs debugging.",err,error,*999)
                        elementNodes(localNodeIdx)=nodeIndices(localNodeIdx1,1,localNodeIdx3)
                      ELSE
                        elementNodes(localNodeIdx)=nodeIndices(localNodeIdx1,localNodeIdx2,localNodeIdx3)
                      ENDIF
                    ENDDO ! localNodeIdx1
                  ENDDO ! localNodeIdx2
                ENDDO ! localNodeIdx3
                elementNumber=elementNumber+1
                CALL GeneratedMeh_ComponentNodesToUserNumbers(cylinderMesh%generatedMesh,basisIdx,elementNodes, &
                  & elementNodesUserNumbers,err,error,*999)
                CALL MeshElements_ElementNodesSet(meshElements,elementNumber,elementNodesUserNumbers,err,error,*999)
              ENDDO ! elementIdx1
            ENDDO ! elementIdx2
          ENDDO ! elementIdx3
          CALL MeshElements_CreateFinish(meshElements,err,error,*999)
        CASE(BASIS_SIMPLEX_TYPE)
          CALL FlagError("Cylinder meshes with simplex basis types is not implemented.",err,error,*999)
        CASE DEFAULT
          CALL FlagError("Basis type is either invalid or not implemented.",err,error,*999)
        END SELECT
      ENDDO !basisIdx
      !Finish the mesh
      CALL Mesh_CreateFinish(generatedMesh%mesh,err,error,*999)
    CASE DEFAULT
      CALL FlagError("Coordinate type is either invalid or not implemented.",err,error,*999)
    END SELECT

    EXITS("GeneratedMesh_CylinderCreateFinish")
    RETURN
    ! TODO invalidate other associations
999 IF(ALLOCATED(nodeIndices)) DEALLOCATE(nodeIndices)
    IF(ALLOCATED(elementIndices)) DEALLOCATE(elementIndices)
    IF(ALLOCATED(numberOfElementsXi)) DEALLOCATE(numberOfElementsXi)
    IF(ALLOCATED(elementNodes)) DEALLOCATE(elementNodes)
    ERRORSEXITS("GeneratedMesh_CylinderCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_CylinderCreateFinish

  !
  !================================================================================================================================
  !

  !>Finalise the cylinder mesh type
  SUBROUTINE GeneratedMesh_CylinderFinalise(cylinderMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("GeneratedMesh_CylinderFinalise",err,error,*999)

    IF(ASSOCIATED(cylinderMesh)) THEN
      IF(ALLOCATED(cylinderMesh%origin)) DEALLOCATE(cylinderMesh%origin)
      IF(ALLOCATED(cylinderMesh%cylinderExtent)) DEALLOCATE(cylinderMesh%cylinderExtent)
      IF(ALLOCATED(cylinderMesh%numberOfElementsXi)) DEALLOCATE(cylinderMesh%numberOfElementsXi)
      IF(ALLOCATED(cylinderMesh%bases)) DEALLOCATE(cylinderMesh%bases)
      DEALLOCATE(cylinderMesh)
    ENDIF

    EXITS("GeneratedMesh_CylinderFinalise")
    RETURN
    ! TODO invalidate other associations
999 ERRORSEXITS("GeneratedMesh_CylinderFinalise",err,error)
    RETURN 1
  END SUBROUTINE GeneratedMesh_CylinderFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the cylinder generated mesh type
  SUBROUTINE GeneratedMesh_CylinderInitialise(generatedMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("GeneratedMesh_CylinderInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is not associated.",err,error,*998)
    IF(ASSOCIATED(generatedMesh%cylinderMesh)) &
      & CALL FlagError("Cylinder mesh is already associated for this generated mesh.",err,error,*998)
      
    ALLOCATE(generatedMesh%cylinderMesh,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate cylinder generated mesh.",err,error,*999)
    generatedMesh%cylinderMesh%generatedMesh=>generatedMesh
    generatedMesh%generatedType=GENERATED_MESH_CYLINDER_MESH_TYPE

    EXITS("GeneratedMesh_CylinderInitialise")
    RETURN
999 CALL GeneratedMesh_CylinderFinalise(generatedMesh%cylinderMesh,dummyErr,dummyError,*998)
998 ERRORSEXITS("GeneratedMesh_CylinderInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_CylinderInitialise

  !
  !================================================================================================================================
  !

  !>Finalise the regular mesh type
  SUBROUTINE GeneratedMesh_RegularFinalise(regularMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("GeneratedMesh_RegularFinalise",err,error,*999)

    IF(ASSOCIATED(regularMesh)) THEN
      IF(ALLOCATED(regularMesh%origin)) DEALLOCATE(regularMesh%origin)
      IF(ALLOCATED(regularMesh%maximumExtent)) DEALLOCATE(regularMesh%maximumExtent)
      IF(ALLOCATED(regularMesh%numberOfElementsXi)) DEALLOCATE(regularMesh%numberOfElementsXi)
      IF(ALLOCATED(regularMesh%baseVectors)) DEALLOCATE(regularMesh%baseVectors)
      IF(ALLOCATED(regularMesh%bases)) DEALLOCATE(regularMesh%bases)
      DEALLOCATE(regularMesh)
    ENDIF

    EXITS("GeneratedMesh_RegularFinalise")
    RETURN
999 ERRORSEXITS("GeneratedMesh_RegularFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_RegularFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the regular generated mesh type
  SUBROUTINE GeneratedMesh_RegularInitialise(generatedMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("GeneratedMesh_RegularInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is not associated.",err,error,*998)
    IF(ASSOCIATED(generatedMesh%regularMesh)) &
      & CALL FlagError("Regular mesh is already associated for this generated mesh.",err,error,*998)
      
    ALLOCATE(generatedMesh%regularMesh,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate regular generated mesh.",err,error,*999)
    generatedMesh%regularMesh%generatedMesh=>generatedMesh
    generatedMesh%generatedType=GENERATED_MESH_REGULAR_MESH_TYPE

    EXITS("GeneratedMesh_RegularInitialise")
    RETURN
999 CALL GeneratedMesh_RegularFinalise(generatedMesh%regularMesh,dummyErr,dummyError,*998)
998 ERRORSEXITS("GeneratedMesh_RegularInitialise",err,error)
    RETURN 1
  END SUBROUTINE GeneratedMesh_RegularInitialise

  !
  !================================================================================================================================
  !

  !>Finalise ellipsoid mesh type
  SUBROUTINE GeneratedMesh_EllipsoidFinalise(ellipsoidMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ellipsoidMesh !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("GeneratedMesh_EllipsoidFinalise",err,error,*999)

    IF(ASSOCIATED(ellipsoidMesh)) THEN
      IF(ALLOCATED(ellipsoidMesh%origin)) DEALLOCATE(ellipsoidMesh%origin)
      IF(ALLOCATED(ellipsoidMesh%ellipsoidExtent)) DEALLOCATE(ellipsoidMesh%ellipsoidExtent)
      IF(ALLOCATED(ellipsoidMesh%numberOfElementsXi)) DEALLOCATE(ellipsoidMesh%numberOfElementsXi)
      IF(ALLOCATED(ellipsoidMesh%bases)) DEALLOCATE(ellipsoidMesh%bases)
      DEALLOCATE(ellipsoidMesh)
    ENDIF

    EXITS("GeneratedMesh_EllipsoidFinalise")
    RETURN
    ! TODO invalidate other associations
999 ERRORSEXITS("GeneratedMesh_EllipsoidFinalise",err,error)
    RETURN 1
  END SUBROUTINE GeneratedMesh_EllipsoidFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the ellipsoid generated mesh type
  SUBROUTINE GeneratedMesh_EllipsoidInitialise(generatedMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("GeneratedMesh_EllipsoidInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is not associated.",err,error,*998)
    IF(ASSOCIATED(generatedMesh%ellipsoidMesh)) &
      & CALL FlagError("Ellipsoid mesh is already associated for this generated mesh.",err,error,*998)
    
    ALLOCATE(generatedMesh%ellipsoidMesh,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate ellipsoid generated mesh.",err,error,*999)
    generatedMesh%ellipsoidMesh%generatedMesh=>generatedMesh
    generatedMesh%generatedType=GENERATED_MESH_ELLIPSOID_MESH_TYPE

    EXITS("GeneratedMesh_EllipsoidInitialise")
    RETURN
999 CALL GeneratedMesh_EllipsoidFinalise(generatedMesh%ellipsoidMesh,dummyErr,dummyError,*998)
998 ERRORSEXITS("GeneratedMesh_EllipsoidInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_EllipsoidInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of a generated mesh. \see OpenCMISS::Iron::cmfe_GeneratedMesh_TypeSet
  SUBROUTINE GeneratedMesh_TypeSet(generatedMesh,generatedType,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to set the type of
    INTEGER(INTG), INTENT(IN) :: generatedType !<The type of mesh to generate \see GeneratedMeshRoutines_GeneratedMeshTypes,GeneratedMeshRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: oldGeneratedType
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_TypeSet",err,error,*999)

    CALL GeneratedMesh_AssertNotFinished(generatedMesh,err,error,*999)
    
    oldGeneratedType=generatedMesh%generatedType
    IF(oldGeneratedType/=generatedType) THEN
      !Initialise the new generated mesh type
      SELECT CASE(generatedType)
      CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
        CALL GeneratedMesh_RegularInitialise(generatedMesh,err,error,*999)
      CASE(GENERATED_MESH_POLAR_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
        CALL GeneratedMesh_CylinderInitialise(generatedMesh,err,error,*999)
      CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
        CALL GeneratedMesh_EllipsoidInitialise(generatedMesh,err,error,*999)
      CASE DEFAULT
        localError="The specified generated mesh mesh type of "//TRIM(NumberToVString(generatedType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Finalise the new generated mesh type
      SELECT CASE(oldGeneratedType)
      CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
        CALL GeneratedMesh_RegularFinalise(generatedMesh%regularMesh,err,error,*999)
      CASE(GENERATED_MESH_POLAR_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
        CALL GeneratedMesh_CylinderFinalise(generatedMesh%cylinderMesh,err,error,*999)
      CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
        CALL GeneratedMesh_EllipsoidFinalise(generatedMesh%ellipsoidMesh,err,error,*999)
      CASE DEFAULT
        localError="The generated mesh mesh type of "//TRIM(NumberToVString(oldGeneratedType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF

    EXITS("GeneratedMesh_TypeSet")
    RETURN
999 ERRORSEXITS("GeneratedMesh_TypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_TypeSet

  !
  !================================================================================================================================
  !

  !>Finalises all generated meshes and deallocates all memory.
  SUBROUTINE GeneratedMeshes_Finalise(generatedMeshes,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshesType), POINTER :: generatedMeshes !<A pointer to the generated meshes to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh

    ENTERS("GeneratedMeshes_Finalise",err,error,*999)

    IF(ASSOCIATED(generatedMeshes)) THEN
      DO WHILE(generatedMeshes%numberOfGeneratedMeshes>0)
        generatedMesh=>generatedMeshes%generatedMeshes(1)%ptr
        CALL GeneratedMesh_Destroy(generatedMesh,err,error,*999)
      ENDDO !generatedMeshIdx
      DEALLOCATE(generatedMeshes)
    ENDIF

    EXITS("GeneratedMeshes_Finalise")
    RETURN
999 ERRORSEXITS("GeneratedMeshes_Finalise",err,error)
    RETURN 1

  END SUBROUTINE GeneratedMeshes_Finalise

  !
  !================================================================================================================================
  !

  !>Intialises all generated meshes.
  SUBROUTINE GeneratedMeshes_InitialiseGeneric(generatedMeshes,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshesType), POINTER :: generatedMeshes !<A pointer to the generated meshes to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("GeneratedMeshes_InitialiseGeneric",err,error,*998)

    IF(ASSOCIATED(generatedMeshes)) CALL FlagError("Generated meshes is already associated.",err,error,*998)
    
    ALLOCATE(generatedMeshes,STAT=err)
    IF(err/=0) CALL FlagError("Generated meshes is not associated.",err,error,*999)
    generatedMeshes%numberOfGeneratedMeshes=0
    NULLIFY(generatedMeshes%region)
    NULLIFY(generatedMeshes%interface)

    EXITS("GeneratedMeshes_InitialiseGeneric")
    RETURN
999 CALL GeneratedMeshes_Finalise(generatedMeshes,dummyErr,dummyError,*998)
998 ERRORSEXITS("GeneratedMeshes_InitialiseGeneric",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMeshes_InitialiseGeneric

  !
  !================================================================================================================================
  !

  !>Intialises the generated meshes for an interface.
  SUBROUTINE GeneratedMeshes_InitialiseInterface(interface,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to initialise the generated meshes for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("GeneratedMeshes_InitialiseInterface",err,error,*999)

    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*999)
    IF(ASSOCIATED(INTERFACE%generatedMeshes)) CALL FlagError("Interface generated meshes is already associated.",err,error,*999)
     
    CALL GeneratedMeshes_InitialiseGeneric(INTERFACE%generatedMeshes,err,error,*999)
    INTERFACE%generatedMeshes%INTERFACE=>interface

    EXITS("GeneratedMeshes_InitialiseInterface")
    RETURN
999 ERRORSEXITS("GeneratedMeshes_InitialiseInterface",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMeshes_InitialiseInterface

  !
  !================================================================================================================================
  !

  !>Intialises the generated meshes for a region.
  SUBROUTINE GeneratedMeshes_InitialiseRegion(region,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("GeneratedMeshes_InitialiseRegion",err,error,*999)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)
    IF(ASSOCIATED(region%generatedMeshes)) CALL FlagError("Region generated meshes is already associated.",err,error,*999)
    
    CALL GeneratedMeshes_InitialiseGeneric(region%generatedMeshes,err,error,*999)
    region%generatedMeshes%region=>region

    EXITS("GeneratedMeshes_InitialiseRegion")
    RETURN
999 ERRORSEXITS("GeneratedMeshes_InitialiseRegion",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMeshes_InitialiseRegion

  !
  !================================================================================================================================
  !

  !>Updates the geometric field parameters from the initial nodal positions of the mesh. \see OpenCMISS::Iron::cmfe_GeneratedMesh_GeometricParametersCalculate
  SUBROUTINE GeneratedMesh_GeometricParametersCalculate(generatedMesh,field,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<The mesh which is generated by the generated mesh
    TYPE(FieldType), POINTER :: field !<A pointer to the field to update the geometric parameters for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_GeometricParametersCalculate",err,error,*999)

    CALL GeneratedMesh_AssertIsFinished(generatedMesh,err,error,*999)
    CALL Field_AssertIsFinished(field,err,error,*999)
    
    SELECT CASE(generatedMesh%generatedType)
    CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
      CALL GeneratedMesh_RegularGeometricParametersCalculate(generatedMesh%regularMesh,field,err,error,*999)
    CASE(GENERATED_MESH_POLAR_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
      CALL GeneratedMesh_CylinderGeometricParametersCalculate(generatedMesh%cylinderMesh,field,err,error,*999)
    CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
      CALL GeneratedMesh_EllipsoidGeometricParametersCalculate(generatedMesh%ellipsoidMesh,field,err,error,*999)
    CASE DEFAULT
      localError="The generated mesh mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("GeneratedMesh_GeometricParametersCalculate")
    RETURN
999 ERRORSEXITS("GeneratedMesh_GeometricParametersCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_GeometricParametersCalculate

  !
  !================================================================================================================================
  !

  !>Updates the geometric field parameters from the initial nodal positions of the regular mesh.
  SUBROUTINE GeneratedMesh_RegularGeometricParametersCalculate(regularMesh,field,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh !<A pointer to the regular mesh object
    TYPE(FieldType), POINTER :: field !<A pointer to the field to update the geometric parameters for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: componentIdx,componentNode,derivativeIdx,interpolationType,meshComponent,nodeIdx,nodePositionIdx(3), &
      & nodeUserNumber,totalNumberOfNodesXi(3),xiIdx
    REAL(DP) :: deltaCoordinate(3,3),myOrigin(3),value,norm
    REAL(DP) :: derivativeValues(MAXIMUM_GLOBAL_DERIV_NUMBER)
    LOGICAL :: nodeExists,ghostNode
    TYPE(BasisType), POINTER :: basis
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_RegularGeometricParametersCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(regularMesh)) CALL FlagError("Regular mesh is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(field%type/=FIELD_GEOMETRIC_TYPE) THEN
      localError="Field number "//TRIM(NumberToVString(FIELD%userNumber,"*",err,error))//" is not a geometric field."
      CALL FlagError(localError,err,error,*999)
    ENDIF
   
    NULLIFY(coordinateSystem)
    CALL Field_CoordinateSystemGet(field,coordinateSystem,err,error,*999)
    IF(coordinateSystem%type/=COORDINATE_RECTANGULAR_CARTESIAN_TYPE) &
      & CALL FlagError("Non rectangular Cartesian coordinate systems are not implemented.",err,error,*999)

    myOrigin=0.0_DP
    myOrigin(1:regularMesh%coordinateDimension)=regularMesh%origin(1:regularMesh%coordinateDimension)
    deltaCoordinate=0.0_DP
    totalNumberOfNodesXi=1
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
    DO componentIdx=1,fieldVariable%numberOfComponents
      CALL FieldVariable_ComponentInterpolationGet(fieldVariable,componentIdx,interpolationType,err,error,*999)
      IF(interpolationType==FIELD_NODE_BASED_INTERPOLATION) THEN
        localError="Component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
          & " of field number "//TRIM(NumberToVString(field%userNumber,"*",err,error))// &
          & " does not have node based interpolation."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      CALL FieldVariable_ComponentMeshComponentGet(fieldVariable,componentIdx,meshComponent,err,error,*999)
      NULLIFY(basis)
      CALL GeneratedMeshRegular_BasisGet(regularMesh,meshComponent,basis,err,error,*999)
      DO xiIdx=1,regularMesh%meshDimension
        totalNumberOfNodesXi(xiIdx)=(basis%numberOfNodesXiC(xiIdx)-2)*regularMesh%numberOfElementsXi(xiIdx)+ &
          & regularMesh%numberOfElementsXi(xiIdx)+1
      ENDDO !xiIdx
      DO xiIdx=1,regularMesh%meshDimension
        deltaCoordinate(1:regularMesh%coordinateDimension,xiIdx)= &
          & regularMesh%baseVectors(1:regularMesh%coordinateDimension,xiIdx)/REAL(totalNumberOfNodesXi(xiIdx)-1,DP)
      ENDDO !xiIdx
      SELECT CASE(field%scalings%scalingType)
      CASE(FIELD_NO_SCALING,FIELD_UNIT_SCALING)
        derivativeValues=0.0_DP
        IF(regularMesh%numberOfElementsXi(1)>0) &
          & derivativeValues(GLOBAL_DERIV_S1)=regularMesh%baseVectors(componentIdx,1)/regularMesh%numberOfElementsXi(1)
        IF(regularMesh%meshDimension>1) THEN
          IF(regularMesh%numberOfElementsXi(2)>0) &
            & derivativeValues(GLOBAL_DERIV_S2)=regularMesh%baseVectors(componentIdx,2)/regularMesh%numberOfElementsXi(2)
          IF(regularMesh%meshDimension>2) THEN
            IF(regularMesh%numberOfElementsXi(3)>0) &
              & derivativeValues(GLOBAL_DERIV_S3)=regularMesh%baseVectors(componentIdx,3)/regularMesh%numberOfElementsXi(3)
          ENDIF
        ENDIF
      CASE DEFAULT
        !Arc length or arithmetic mean scaling
        derivativeValues=0.0_DP
        IF(regularMesh%numberOfElementsXi(1)>0) THEN
          CALL L2Norm(regularMesh%baseVectors(:,1),norm,err,error,*999)
          derivativeValues(GLOBAL_DERIV_S1)=regularMesh%baseVectors(componentIdx,1)/norm
        END IF
        IF(regularMesh%meshDimension>1) THEN
          IF(regularMesh%numberOfElementsXi(2)>0) THEN
            CALL L2Norm(regularMesh%baseVectors(:,2),norm,err,error,*999)
            derivativeValues(GLOBAL_DERIV_S2)=regularMesh%baseVectors(componentIdx,2)/norm
          ENDIF
          IF(regularMesh%meshDimension>2) THEN
            IF(regularMesh%numberOfElementsXi(3)>0) THEN
              CALL L2Norm(regularMesh%baseVectors(:,3),norm,err,error,*999)
              derivativeValues(GLOBAL_DERIV_S3)=regularMesh%baseVectors(componentIdx,3)/norm
            ENDIF
          ENDIF
        ENDIF
      END SELECT
      !Update geometric parameters in this computation domain only
      NULLIFY(domain)
      CALL FieldVariable_ComponentDomainGet(fieldVariable,componentIdx,domain,err,error,*999)
      NULLIFY(domainTopology)
      CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
      NULLIFY(domainNodes)
      CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
      DO componentNode=1,totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)*totalNumberOfNodesXi(3)
        !Regular meshes with Lagrange/Hermite elements use different node numberings to other mesh types
        IF(basis%type==BASIS_LAGRANGE_HERMITE_TP_TYPE) THEN
          CALL GeneratedMesh_RegularComponentNodeToUserNumber(regularMesh%generatedMesh,meshComponent,componentNode, &
            & nodeUserNumber,err,error,*999)
        ELSE
          nodeUserNumber=GeneratedMesh_ComponentNodeToUserNumber(regularMesh%generatedMesh,meshComponent,componentNode,err,error)
        ENDIF
        CALL DomainNodes_NodeCheckExists(domainNodes,nodeUserNumber,nodeExists,nodeIdx,ghostNode,err,error,*999)
        IF(nodeExists.AND..NOT.ghostNode) THEN
          nodePositionIdx(3)=(componentNode-1)/(totalNumberOfNodesXi(2)*totalNumberOfNodesXi(1))+1
          nodePositionIdx(2)=MOD(componentNode-1,totalNumberOfNodesXi(2)*totalNumberOfNodesXi(1))/totalNumberOfNodesXi(1)+1
          nodePositionIdx(1)=MOD(MOD(componentNode-1,totalNumberOfNodesXi(2)*totalNumberOfNodesXi(1)),totalNumberOfNodesXi(1))+1
          value=0.0_DP
          DO xiIdx=1,regularMesh%meshDimension
            value=value+REAL(nodePositionIdx(xiIdx)-1,DP)*deltaCoordinate(componentIdx,xiIdx)
          ENDDO !xiIdx
          value=myOrigin(componentIdx)+value
          CALL FieldVariable_ParameterSetUpdateNode(fieldVariable,FIELD_VALUES_SET_TYPE,1,1,nodeUserNumber,componentIdx, &
            & value,err,error,*999)
          !Set derivatives
          DO derivativeIdx=2,domainNodes%nodes(nodeIdx)%numberOfDerivatives
            CALL FieldVariable_ParameterSetUpdateNode(fieldVariable,FIELD_VALUES_SET_TYPE,1,derivativeIdx,nodeUserNumber, &
              & componentIdx,derivativeValues(derivativeIdx),err,error,*999)
          ENDDO !derivativeIdx
        ENDIF !node_exists
      ENDDO !nodeIdx
    ENDDO !componentIdx
!!TODO: do boundary nodes first then start the update to overlap computation and computation.
    CALL FieldVariable_ParameterSetUpdateStart(fieldVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
    CALL FieldVariable_ParameterSetUpdateFinish(fieldVariable,FIELD_VALUES_SET_TYPE,err,error,*999)

    EXITS("GeneratedMesh_RegularGeometricParametersCalculate")
    RETURN
999 ERRORS("GeneratedMesh_RegularGeometricParametersCalculate",err,error)
    EXITS("GeneratedMesh_RegularGeometricParametersCalculate")
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_RegularGeometricParametersCalculate

  !
  !================================================================================================================================
  !

  !>Updates the geometric field parameters from the initial nodal positions of the mesh. Derivatives are averaged via straight line approximation, except for circumferential component
  SUBROUTINE GeneratedMesh_CylinderGeometricParametersCalculate(cylinderMesh,field,err,error,*)

    ! Argument variables
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh !<A pointer to the cylinder mesh object
    TYPE(FieldType), POINTER :: field !<A pointer to the field to update the geometric parameters for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    ! Local variables
    INTEGER(INTG) :: numberOfElementsXi(3),numberOfNodesXic(3)
    INTEGER(INTG) :: totalNumberOfNodesXi(3),interpolationTypes(3)
    INTEGER(INTG) :: componentIdx,xiIdx
    INTEGER(INTG) :: nodeNumber,globalNodeNumber,componentNodeNumber,localDOF,derivativeIdx
    INTEGER(INTG) :: numberOfPlanerNodes,scalingType,meshComponent
    INTEGER(INTG), ALLOCATABLE :: nodeIndices(:,:,:),elementIndices(:,:,:)
    INTEGER(INTG) :: nodeIdx(3) ! holds r,theta,z indices
    REAL(DP) :: delta(3),deltai(3),polarCoordinates(3),rcCoordinates(3)
    REAL(DP) :: cylinderExtent(3),derivative
    TYPE(BasisType), POINTER :: basis
    TYPE(DomainType),POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_CylinderGeometricParametersCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(cylinderMesh)) CALL FlagError("Cylinder mesh is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(field%TYPE/=FIELD_GEOMETRIC_TYPE) THEN
      localError="Field number "//TRIM(NumberToVString(field%userNumber,"*",err,error))//" is not a geometric field."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
    IF(fieldVariable%numberOfComponents/=3) CALL FlagError("Geometric field must be three dimensional.",err,error,*999)    
    CALL Field_ScalingTypeGet(field,scalingType,err,error,*999)
    IF(scalingType/=FIELD_UNIT_SCALING) &
      & CALL WriteString(GENERAL_OUTPUT_TYPE,"  Note: If the cylinder looks wonky, set field scaling to&
      & unit scaling type.",err,error,*999)
    DO componentIdx=1,3
      CALL FieldVariable_ComponentInterpolationGet(fieldVariable,componentIdx,interpolationTypes(componentIdx),err,error,*999)
    ENDDO !componentIdx
    IF(.NOT.ALL(interpolationTypes==FIELD_NODE_BASED_INTERPOLATION)) &
      & CALL FlagError("All field variable components must have node-based interpolation.",err,error,*999)
    
    DO componentIdx=1,fieldVariable%numberOfComponents
      CALL FieldVariable_ComponentMeshComponentGet(fieldVariable,componentIdx,meshComponent,err,error,*999)
      NULLIFY(basis)
      CALL GeneratedMeshCylinder_BasisGet(cylinderMesh,meshComponent,basis,err,error,*999)
      numberOfElementsXi=cylinderMesh%numberOfElementsXi
      numberOfNodesXic=basis%numberOfNodesXiC
      DO xiIdx=1,3
        totalNumberOfNodesXi(xiIdx)=(numberOfNodesXic(xiIdx)-1)*numberOfElementsXi(xiIdx)+1
      ENDDO !xiIdx
      totalNumberOfNodesXi(2)=totalNumberOfNodesXi(2)-1 ! theta loops around so slightly different
      numberOfPlanerNodes=totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
      NULLIFY(domain)
      CALL FieldVariable_ComponentDomainGet(fieldVariable,componentIdx,domain,err,error,*999)
      NULLIFY(domainTopology)
      CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
      NULLIFY(domainNodes)
      CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
      ! calculate deltai now
      cylinderExtent=cylinderMesh%cylinderExtent
      delta(1)=(cylinderExtent(2)-cylinderExtent(1))/numberOfElementsXi(1)
      delta(2)=TWOPI/numberOfElementsXi(2)
      delta(3)=cylinderExtent(3)/numberOfElementsXi(3)
      DO xiIdx=1,3
        deltai(xiIdx)=delta(xiIdx)/(numberOfNodesXic(xiIdx)-1)
      ENDDO !xiIdx
      DO nodeNumber=1,domainNodes%numberOfNodes
        globalNodeNumber=domainNodes%nodes(nodeNumber)%globalNumber
        componentNodeNumber=GeneratedMesh_UserNumberToComponentNode(cylinderMesh%generatedMesh,meshComponent,globalNodeNumber, &
          & err,error)
        ! calculate nodeIdx which will be used to calculate (r,theta,z) then (x,y,z)
        componentNodeNumber=componentNodeNumber-1 ! let's go 0-based index for a bit
        nodeIdx(3)=componentNodeNumber/numberOfPlanerNodes
        nodeIdx(2)=(componentNodeNumber-(nodeIdx(3))*numberOfPlanerNodes)/totalNumberOfNodesXi(1)
        nodeIdx(1)=MOD(componentNodeNumber-(nodeIdx(3))*numberOfPlanerNodes,totalNumberOfNodesXi(1))
        DO xiIdx=1,3
          polarCoordinates(xiIdx)=nodeIdx(xiIdx)*deltai(xiIdx)
        ENDDO !xiIdx
        polarCoordinates(1)=nodeIdx(1)*deltai(1)+cylinderExtent(1) ! add the inner radius
        rcCoordinates(1)=polarCoordinates(1)*COS(polarCoordinates(2))
        rcCoordinates(2)=polarCoordinates(1)*SIN(polarCoordinates(2))
        rcCoordinates(3)=polarCoordinates(3)
        rcCoordinates=rcCoordinates+cylinderMesh%origin
        !Default to version 1 of each node derivative
        CALL FieldVariable_LocalNodeDOFGet(fieldVariable,1,1,nodeNumber,componentIdx,localDOF,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,localDOF,rcCoordinates(componentIdx), &
          & err,error,*999)
        ! Do derivatives: if there are derivatives, we can assume it's cubic hermite
        !   given that quadratic hermites are only used for collapsed hex elements,
        !   but NB mixed bases have to be handled (e.g. CH-CH-linear combinations)
        IF(domainNodes%nodes(nodeNumber)%numberOfDerivatives>1) THEN
          ! Since I decided how xi 1,2,3 line up with the cylinder polar coordinates,
          ! we know a priori that only some of the derivatives are nonzero (analytically).
          ! NOTE: if hermite type used, should assign FIELD_UNIT_SCALING type for this to work
          DO derivativeIdx=2,domainNodes%nodes(nodeNumber)%numberOfDerivatives
            SELECT CASE(domainNodes%nodes(nodeNumber)%derivatives(derivativeIdx)%globalDerivativeIndex)
            CASE(GLOBAL_DERIV_S1)
              SELECT CASE(componentIdx)
              CASE(1)
                derivative=COS(polarCoordinates(2))*delta(1)
              CASE(2)
                derivative=SIN(polarCoordinates(2))*delta(1)
              CASE DEFAULT
                derivative=0.0_DP
              END SELECT
            CASE(GLOBAL_DERIV_S2)
              SELECT CASE(componentIdx)
              CASE(1)
                derivative=-polarCoordinates(1)*SIN(polarCoordinates(2))*delta(2)
              CASE(2)
                derivative=polarCoordinates(1)*COS(polarCoordinates(2))*delta(2)
              CASE DEFAULT
                derivative=0.0_DP
              END SELECT
            CASE(GLOBAL_DERIV_S3)
              IF(componentIdx==3) THEN
                derivative=delta(3)
              ELSE
                derivative=0.0_DP
              ENDIF
            CASE(GLOBAL_DERIV_S1_S2)
              SELECT CASE(componentIdx)
              CASE(1)
                derivative=-SIN(polarCoordinates(2))*delta(1)*delta(2)
              CASE(2)
                derivative=COS(polarCoordinates(2))*delta(1)*delta(2)
              CASE DEFAULT
                derivative=0.0_DP
              END SELECT
            CASE DEFAULT  ! all other non-xy-planar cross derivatives
              derivative=0.0_DP
            END SELECT
            ! assign derivative
            !Default to version 1 of each node derivative
            CALL FieldVariable_LocalNodeDOFGet(fieldVariable,1,derivativeIdx,nodeNumber,componentIdx,localDOF,err,error,*999)
            CALL FieldVariable_ParameterSetUpdateLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,localDOF,derivative,err,error,*999)
          ENDDO !derivativeIdx
        ENDIF !derivatives
      ENDDO !nodeNumber
    ENDDO !componentIdx
    CALL FieldVariable_ParameterSetUpdateStart(fieldVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
    CALL FieldVariable_ParameterSetUpdateFinish(fieldVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
    
    ! all done
    IF(ALLOCATED(nodeIndices)) DEALLOCATE(nodeIndices)
    IF(ALLOCATED(elementIndices)) DEALLOCATE(elementIndices)

    EXITS("GeneratedMesh_CylinderGeometricParametersCalculate")
    RETURN
999 IF(ALLOCATED(nodeIndices)) DEALLOCATE(nodeIndices)
    IF(ALLOCATED(elementIndices)) DEALLOCATE(elementIndices)
    ERRORS("GeneratedMesh_CylinderGeometricParametersCalculate",err,error)
    EXITS("GeneratedMesh_CylinderGeometricParametersCalculate")
    RETURN 1

  END SUBROUTINE GeneratedMesh_CylinderGeometricParametersCalculate

  !
  !================================================================================================================================
  !

  !>Updates the geometric field parameters from the initial nodal positions of the mesh.
  !>Derivatives are averaged via straight line approximation, except for circumferential component
  SUBROUTINE GeneratedMesh_EllipsoidGeometricParametersCalculate(ellipsoidMesh,field,err,error,*)

    ! Argument variables
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ellipsoidMesh !<A pointer to the ellipsoid mesh object
    TYPE(FieldType), POINTER :: field !<A pointer to the field to update the geometric parameters for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    ! Local variables
    INTEGER(INTG) :: basisIdx, componentIdx,domainNumber,globalNodeNumber,interpolationTypes(3),i,j,k,localNode,meshComponent, &
      & meshComponent2,myComputationNode,nodeNumber,numberOfElementsXi(3),numberOfNodesXic(3),scalingType, &
      & totalNumberOfNodesXi(3),xiIdx
    INTEGER(INTG), ALLOCATABLE :: elementIndices(:,:,:),nodeIndices(:,:,:)
    REAL(DP) :: alpha,delta(3),deltai(3),ellipsoidExtent(4),nu,phi,rcCoordinates(3),t,xi,x,y,z
    TYPE(BasisType), POINTER :: basis
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType),POINTER :: domain
    TYPE(DomainMappingType), POINTER :: nodesMapping
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError
    TYPE(WorkGroupType), POINTER :: workGroup

    ENTERS("GeneratedMesh_EllipsoidGeometricParametersCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(ellipsoidMesh)) CALL FlagError("Ellipsoid mesh is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(field%TYPE/=FIELD_GEOMETRIC_TYPE) THEN
      localError="Field number "//TRIM(NumberToVString(field%userNumber,"*",err,error))//" is not a geometric field."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    CALL Field_ScalingTypeGet(field,scalingType,err,error,*999)
    IF(scalingType/=FIELD_UNIT_SCALING) &
      & CALL WriteString(GENERAL_OUTPUT_TYPE,"  Note: If the ellipsoid looks wonky, set field scaling to &
      & unit scaling type.",err,error,*999)
     
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(field,decomposition,err,error,*999)
    NULLIFY(workGroup)
    CALL Decomposition_WorkGroupGet(decomposition,workGroup,err,error,*999)    
    CALL WorkGroup_GroupNodeNumberGet(workGroup,myComputationNode,err,error,*999)
    
    ! assign to the field
    nodeNumber=0
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
    IF(fieldVariable%numberOfComponents==3) CALL FlagError("Geometric field must be three dimensional.",err,error,*999)
    CALL FieldVariable_ComponentMeshComponentGet(fieldVariable,1,meshComponent,err,error,*999)
    DO componentIdx=2,3
      CALL FieldVariable_ComponentMeshComponentGet(fieldVariable,componentIdx,meshComponent2,err,error,*999)
      IF(meshComponent2/=meshComponent) &
        & CALL FlagError("Multiple mesh components for geometric components is not implemented.",err,error,*999)
    ENDDO !componentIdx
    DO componentIdx=1,3
      CALL FieldVariable_ComponentInterpolationGet(fieldVariable,componentIdx,interpolationTypes(componentIdx),err,error,*999)
    ENDDO !componentIdx
    IF(.NOT.ALL(interpolationTypes==FIELD_NODE_BASED_INTERPOLATION)) &
      & CALL FlagError("All field variable components must have node-based interpolation.",err,error,*999)
    
    basisIdx=meshComponent*2-1
    NULLIFY(basis)
    CALL GeneratedMeshEllipsoid_BasisGet(ellipsoidMesh,basisIdx,basis,err,error,*999)
    !< Ellipsoid_extent= inner long axis, inner short axis, wall thickness, top angle (from 0)
    ! calculate the total number of nodes in each xi direction
    numberOfElementsXi=ellipsoidMesh%numberOfElementsXi
    numberOfNodesXic=basis%numberOfNodesXiC
    DO xiIdx=1,3
      totalNumberOfNodesXi(xiIdx)=(numberOfNodesXic(xiIdx)-1)*numberOfElementsXi(xiIdx)+1
    ENDDO !xiIdx
    totalNumberOfNodesXi(1)=totalNumberOfNodesXi(1)-1 ! theta loops around so slightly different
    ! calculate deltai now
    ellipsoidExtent=ellipsoidMesh%ellipsoidExtent
    delta(1)=TWOPI/numberOfElementsXi(1)
    delta(2)=(PI-ellipsoidExtent(4))/numberOfElementsXi(2)
    delta(3)=ellipsoidExtent(3)/numberOfElementsXi(3)
    DO xiIdx=1,3
      deltai(xiIdx)=delta(xiIdx)/(numberOfNodesXic(xiIdx)-1)
    ENDDO !xiIdx
    ! numberOfPlanerNodes=totalNumberOfNodesXi(1)*totalNumberOfNodesXi(2)
    NULLIFY(domain)
    CALL FieldVariable_ComponentDomainGet(fieldVariable,1,domain,err,error,*999)
    NULLIFY(domainTopology)
    CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(domainMappings)
    CALL Domain_DomainMappingsGet(domain,domainMappings,err,error,*999)
    NULLIFY(domainNodes)
    CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
    NULLIFY(nodesMapping)
    CALL DomainMappings_NodesMappingGet(domainMappings,nodesMapping,err,error,*999)
    IF (ellipsoidExtent(1)>ellipsoidExtent(2)) THEN
      !Prolate spheroid
      k=1
      !inner surface
      alpha=SQRT((ellipsoidExtent(1))**2-(ellipsoidExtent(2))**2)
      !xi=log(ellipsoidExtent(1)/alpha+sqrt((ellipsoidExtent(1)/alpha)**2+1))
      xi=ACOSH(ellipsoidExtent(1)/alpha)
      
      j=1
      !apex node
      nodeNumber=1
      globalNodeNumber=GeneratedMesh_ComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,basisIdx,nodeNumber,err,error)
      CALL Decomposition_NodeDomainGet(decomposition,globalNodeNumber,meshComponent,domainNumber,err,error,*999)
      IF(domainNumber==myComputationNode) THEN
        rcCoordinates(1)=0.0_DP
        rcCoordinates(2)=0.0_DP
        rcCoordinates(3)=-ellipsoidExtent(1)
        DO componentIdx=1,fieldVariable%numberOfComponents
          !Default to version 1 of each node derivative
          CALL FieldVariable_ParameterSetUpdateNode(fieldVariable,FIELD_VALUES_SET_TYPE,1,1,globalNodeNumber,componentIdx, &
            & rcCoordinates(componentIdx),err,error,*999)
          localNode=nodesMapping%globalToLocalMap(globalNodeNumber)%localNumber(1)
          IF(domainNodes%nodes(localNode)%numberOfDerivatives>1) &
            & CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
        ENDDO !componentIdx
      ENDIF
      
      DO j=2,totalNumberOfNodesXi(2)
        !longitudinal loop
        nu=PI-deltai(2)*(j-1)
        DO i=1,totalNumberOfNodesXi(1)
          !circumferential loop
          phi=deltai(1)*(i-1)
          rcCoordinates(1)=alpha*(SINH(xi)*SIN(nu)*COS(phi))
          rcCoordinates(2)=alpha*(SINH(xi)*SIN(nu)*SIN(phi))
          rcCoordinates(3)=alpha*(COSH(xi)*COS(nu))
          nodeNumber=nodeNumber+1
          globalNodeNumber=GeneratedMesh_ComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,basisIdx,nodeNumber,err,error)
          CALL Decomposition_NodeDomainGet(decomposition,globalNodeNumber,meshComponent,domainNumber,err,error,*999)
          IF(domainNumber==myComputationNode) THEN
            DO componentIdx=1,fieldVariable%numberOfComponents
              CALL FieldVariable_ParameterSetUpdateNode(fieldvariable,FIELD_VALUES_SET_TYPE,1,1,globalNodeNumber, &
                & componentIdx,rcCoordinates(componentIdx),err,error,*999)
              localNode=nodesMapping%globalToLocalMap(globalNodeNumber)%localNumber(1)
              IF(domainNodes%nodes(localNode)%numberOfDerivatives>1) &
                & CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
            ENDDO !componentIdx
          ENDIF
        ENDDO !i
      ENDDO !j

      DO k=2,totalNumberOfNodesXi(3)
        !transmural loop
        j=1
        !apex nodes
        rcCoordinates(1)=0.0_DP
        rcCoordinates(2)=0.0_DP
        rcCoordinates(3)=-ellipsoidExtent(1)-(k-1)*(deltai(3))
        nodeNumber=nodeNumber+1
        globalNodeNumber=GeneratedMesh_ComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,basisIdx,nodeNumber,err,error)
        CALL Decomposition_NodeDomainGet(decomposition,globalNodeNumber,meshComponent,domainNumber,err,error,*999)
        IF(domainNumber==myComputationNode) THEN
          DO componentIdx=1,fieldVariable%numberOfComponents
            CALL FieldVariable_ParameterSetUpdateNode(fieldVariable,FIELD_VALUES_SET_TYPE,1,1,globalNodeNumber, &
              & componentIdx,rcCoordinates(componentIdx),err,error,*999)
            localNode=nodesMapping%globalToLocalMap(globalNodeNumber)%localNumber(1)
            IF(domainNodes%nodes(localNode)%numberOfDerivatives>1) &
              & CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
          ENDDO !componentIdx
        ENDIF
        
        DO j=2,totalNumberOfNodesXi(2)
          !longitudinal loop
          nu=PI-deltai(2)*(j-1)
          DO i=1,totalNumberOfNodesXi(1)
            !circumferential loop
            phi=deltai(1)*(i-1)
            x=alpha*(SINH(xi)*SIN(nu)*COS(phi))
            y=alpha*(SINH(xi)*SIN(nu)*SIN(phi))
            z=alpha*(COSH(xi)*COS(nu))
            !Normal vector from inner surface with length deltai(3)(k-1)
            ! Finney&Thomas: Calculus, second edition, Addison-Wesley Publishing Company, 1994, page 847
            !X=x(1+2t/a^2) Y=y(1+2t/a^2) Z=z(1+2t/c^2)
            t=(deltai(3)*(k-1))/SQRT((4*x**2/(ellipsoidExtent(2))**4)+ &
              & (4*y**2/(ellipsoidExtent(2))**4)+(4*z**2/(ellipsoidExtent(1))**4))
            rcCoordinates(1)=x*(1+2*t/(ellipsoidExtent(2))**2)
            rcCoordinates(2)=y*(1+2*t/(ellipsoidExtent(2))**2)
            rcCoordinates(3)=z*(1+2*t/(ellipsoidExtent(1))**2)
            nodeNumber=nodeNumber+1
            globalNodeNumber=GeneratedMesh_ComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,basisIdx,nodeNumber,err,error)
            CALL Decomposition_NodeDomainGet(decomposition,globalNodeNumber,meshComponent,domainNumber,err,error,*999)
            IF(domainNumber==myComputationNode) THEN
              DO componentIdx=1,fieldVariable%numberOfComponents
                CALL FieldVariable_ParameterSetUpdateNode(fieldVariable,FIELD_VALUES_SET_TYPE,1,1,globalNodeNumber, &
                  & componentIdx,rcCoordinates(componentIdx),err,error,*999)
                localNode=nodesMapping%globalToLocalMap(globalNodeNumber)%localNumber(1)
                IF(domainNodes%nodes(localNode)%numberOfDerivatives>1) &
                  & CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
               ENDDO !componentIdx
            ENDIF
          ENDDO !i
        ENDDO !j
      ENDDO !k
    ELSE IF(ABS(ellipsoidExtent(1)-ellipsoidExtent(2))<ZERO_TOLERANCE) THEN
      !Sphere
      nodeNumber=0
      DO k=1,totalNumberOfNodesXi(3)
        !transmural loop
        alpha=ellipsoidExtent(1)+(k-1)*(deltai(3))
        j=1
        !apex nodes
        rcCoordinates(1)=0.0_DP
        rcCoordinates(2)=0.0_DP
        rcCoordinates(3)=-alpha
        nodeNumber=nodeNumber+1
        globalNodeNumber=GeneratedMesh_ComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,basisIdx,nodeNumber,err,error)
        CALL Decomposition_NodeDomainGet(decomposition,globalNodeNumber,meshComponent,domainNumber,err,error,*999)
        IF(domainNumber==myComputationNode) THEN
          DO componentIdx=1,fieldVariable%numberOfComponents
            CALL FieldVariable_ParameterSetUpdateNode(fieldvariable,FIELD_VALUES_SET_TYPE,1,1,globalNodeNumber, &
              & componentIdx,rcCoordinates(componentIdx),err,error,*999)
            localNode=nodesMapping%globalToLocalMap(globalNodeNumber)%localNumber(1)
            IF(domainNodes%nodes(localNode)%numberOfDerivatives>1) &
              & CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
          ENDDO !componentIdx
        ENDIF
        
        DO j=2,totalNumberOfNodesXi(2)
          !longitudinal loop
          nu=PI-deltai(2)*(j-1)
          DO i=1,totalNumberOfNodesXi(1)
            !circumferential loop
            phi=deltai(1)*(i-1)
            rcCoordinates(1)=alpha*SIN(nu)*COS(phi)
            rcCoordinates(2)=alpha*SIN(nu)*SIN(phi)
            rcCoordinates(3)=alpha*COS(nu)
            nodeNumber=nodeNumber+1
            globalNodeNumber=GeneratedMesh_ComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,basisIdx,nodeNumber,err,error)
            CALL Decomposition_NodeDomainGet(decomposition,globalNodeNumber,meshComponent,domainNumber,err,error,*999)
            IF(domainNumber==myComputationNode) THEN
              DO componentIdx=1,fieldVariable%numberOfComponents
                CALL FieldVariable_ParameterSetUpdateNode(fieldVariable,FIELD_VALUES_SET_TYPE,1,1,globalNodeNumber, &
                  & componentIdx,rcCoordinates(componentIdx),err,error,*999)
                localNode=nodesMapping%globalToLocalMap(globalNodeNumber)%localNumber(1)
                IF(domainNodes%nodes(localNode)%numberOfDerivatives>1) &
                  & CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
              ENDDO !componentIdx
            ENDIF
          ENDDO !i
        ENDDO !j
      ENDDO !k

    ELSE IF(ellipsoidExtent(1)<ellipsoidExtent(2)) THEN
      !Oblate spheroid
      k=1
      !inner surface
      alpha=SQRT((ellipsoidExtent(2))**2-(ellipsoidExtent(1))**2)
      !xi=log(ellipsoidExtent(1)/alpha+sqrt((ellipsoidExtent(1)/alpha)**2+1))
      xi=ACOSH(ellipsoidExtent(2)/alpha)
      
      j=1
      !apex node
      nodeNumber=1
      globalNodeNumber=GeneratedMesh_ComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,basisIdx,nodeNumber,err,error)
      CALL Decomposition_NodeDomainGet(decomposition,globalNodeNumber,meshComponent,domainNumber,err,error,*999)
      IF(domainNumber==myComputationNode) THEN
        rcCoordinates(1)=0.0_DP
        rcCoordinates(2)=0.0_DP
        rcCoordinates(3)=-ellipsoidExtent(1)
        DO componentIdx=1,fieldVariable%numberOfComponents
          CALL FieldVariable_ParameterSetUpdateNode(fieldVariable,FIELD_VALUES_SET_TYPE,1,1,globalNodeNumber, &
            & componentIdx,rcCoordinates(componentIdx),err,error,*999)
          localNode=nodesMapping%globalToLocalMap(globalNodeNumber)%localNumber(1)
          IF(domainNodes%nodes(localNode)%numberOfDerivatives>1) &
            & CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
        ENDDO !componentIdx
      ENDIF
      
      DO j=2,totalNumberOfNodesXi(2)
        !longitudinal loop
        nu=-PI/2+deltai(2)*(j-1)
        DO i=1,totalNumberOfNodesXi(1)
          !circumferential loop
          phi=deltai(1)*(i-1)
          rcCoordinates(1)=alpha*(COSH(xi)*COS(nu)*COS(phi))
          rcCoordinates(2)=alpha*(COSH(xi)*COS(nu)*SIN(phi))
          rcCoordinates(3)=alpha*(SINH(xi)*SIN(nu))
          nodeNumber=nodeNumber+1
          globalNodeNumber=GeneratedMesh_ComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,basisIdx,nodeNumber,err,error)
          CALL Decomposition_NodeDomainGet(decomposition,globalNodeNumber,meshComponent,domainNumber,err,error,*999)
          IF(domainNumber==myComputationNode) THEN
            DO componentIdx=1,fieldVariable%numberOfComponents
              CALL FieldVariable_ParameterSetUpdateNode(fieldVariable,FIELD_VALUES_SET_TYPE,1,1,globalNodeNumber, &
                & componentIdx,rcCoordinates(componentIdx),err,error,*999)
              localNode=nodesMapping%globalToLocalMap(globalNodeNumber)%localNumber(1)
              IF(domainNodes%nodes(localNode)%numberOfDerivatives>1) &
                & CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
            ENDDO !componentIdx
          ENDIF
        ENDDO !i
      ENDDO !j
      
      DO k=2,totalNumberOfNodesXi(3)
        !transmural loop
        j=1
        !apex nodes
        rcCoordinates(1)=0.0_DP
        rcCoordinates(2)=0.0_DP
        rcCoordinates(3)=-ellipsoidExtent(1)-(k-1)*(deltai(3))
        nodeNumber=nodeNumber+1
        globalNodeNumber=GeneratedMesh_ComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,basisIdx,nodeNumber,err,error)
        CALL Decomposition_NodeDomainGet(decomposition,globalNodeNumber,meshComponent,domainNumber,err,error,*999)
        IF(domainNumber==myComputationNode) THEN
          DO componentIdx=1,fieldVariable%numberOfComponents
            CALL FieldVariable_ParameterSetUpdateNode(fieldVariable,FIELD_VALUES_SET_TYPE,1,1,globalNodeNumber, &
              & componentIdx,rcCoordinates(componentIdx),err,error,*999)
            localNode=nodesMapping%globalToLocalMap(globalNodeNumber)%localNumber(1)
            IF(domainNodes%nodes(localNode)%numberOfDerivatives>1) &
              & CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
          ENDDO !componentIdx
        ENDIF
        
        DO j=2,totalNumberOfNodesXi(2)
          !longitudinal loop
          nu=-PI/2+deltai(2)*(j-1)
          DO i=1,totalNumberOfNodesXi(1)
            !circumferential loop
            phi=deltai(1)*(i-1)
            x=alpha*(COSH(xi)*COS(nu)*COS(phi))
            y=alpha*(COSH(xi)*COS(nu)*SIN(phi))
            z=alpha*(SINH(xi)*SIN(nu))
            !Normal vector from inner surface with length deltai(3)(k-1)
            ! Finney&Thomas: Calculus, second edition, Addison-Wesley Publishing Company, 1994, page 847
            !X=x(1+2t/a^2) Y=y(1+2t/a^2) Z=z(1+2t/c^2)
            t=(deltai(3)*(k-1))/SQRT((4*x**2/(ellipsoidExtent(2))**4)+ &
              & (4*y**2/(ellipsoidExtent(2))**4)+(4*z**2/(ellipsoidExtent(1))**4))
            rcCoordinates(1)=x*(1+2*t/(ellipsoidExtent(2))**2)
            rcCoordinates(2)=y*(1+2*t/(ellipsoidExtent(2))**2)
            rcCoordinates(3)=z*(1+2*t/(ellipsoidExtent(1))**2)
            nodeNumber=nodeNumber+1
            globalNodeNumber=GeneratedMesh_ComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,basisIdx,nodeNumber,err,error)
            CALL Decomposition_NodeDomainGet(decomposition,globalNodeNumber,meshComponent,domainNumber,err,error,*999)
            IF(domainNumber==myComputationNode) THEN
              DO componentIdx=1,fieldVariable%numberOfComponents
                CALL FieldVariable_ParameterSetUpdateNode(fieldVariable,FIELD_VALUES_SET_TYPE,1,1,globalNodeNumber, &
                  & componentIdx,rcCoordinates(componentIdx),err,error,*999)
                localNode=nodesMapping%globalToLocalMap(globalNodeNumber)%localNumber(1)
                IF(domainNodes%nodes(localNode)%numberOfDerivatives>1) &
                  & CALL FlagError("Not generalized to hermittean elements.",err,error,*999)
              ENDDO !componentIdx
            ENDIF
          ENDDO !i
        ENDDO !j
      ENDDO !k
    ELSE
      CALL FlagError("Not valid long axis - short axis relation",err,error,*999)
    ENDIF
    CALL FieldVariable_ParameterSetUpdateStart(fieldVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
    CALL FieldVariable_ParameterSetUpdateFinish(fieldVariable,FIELD_VALUES_SET_TYPE,err,error,*999)

    ! all done
    IF(ALLOCATED(nodeIndices)) DEALLOCATE(nodeIndices)
    IF(ALLOCATED(elementIndices)) DEALLOCATE(elementIndices)

    EXITS("GeneratedMesh_EllipsoidGeometricParametersCalculate")
    RETURN
999 IF(ALLOCATED(nodeIndices)) DEALLOCATE(nodeIndices)
    IF(ALLOCATED(elementIndices)) DEALLOCATE(elementIndices)
    ERRORS("GeneratedMesh_EllipsoidGeometricParametersCalculate",err,error)
    EXITS("GeneratedMesh_EllipsoidGeometricParametersCalculate")
    RETURN 1

  END SUBROUTINE GeneratedMesh_EllipsoidGeometricParametersCalculate

  !
  !================================================================================================================================
  !

  !>Provides an easy way to grab surfaces for boundary condition assignment
  SUBROUTINE GeneratedMesh_RegularSurfaceGet(regularMesh,meshComponent,surfaceType,surfaceNodes,normalXi,err,error,*)

    ! Argument variables
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh !<A pointer to the regular mesh object
    INTEGER(INTG), INTENT(IN) :: meshComponent !<The mesh component to get the surface for.
    INTEGER(INTG), INTENT(IN) :: surfaceType !<A constant identifying the type of surface to get \see GeneratedMeshRoutines_GeneratedMeshRegularSurfaces,GeneratedMeshRoutines
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: surfaceNodes(:) !<On exit, contains the list of nodes belonging to the surface
    INTEGER(INTG), INTENT(OUT) :: normalXi !<On exit, contains the xi direction of the outward pointing normal of the surface
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    ! Local variables
    INTEGER(INTG) :: nodeCounter,i,j,k
    INTEGER(INTG) :: numberOfElementsXi(3) !Specified number of elements in each xi direction
    INTEGER(INTG) :: numberOfNodesXi(3) ! Number of nodes per element in each xi direction (basis property)
    INTEGER(INTG) :: numberOfDimensions,totalNumberOfNodes,totalNumberOfElements,nodeUserNumber
    INTEGER(INTG),ALLOCATABLE :: nodeIndices(:,:,:),elementIndices(:,:,:)
    REAL(DP) :: delta(3),deltaI(3)
    TYPE(BasisType), POINTER :: basis
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("GeneratedMesh_RegularSurfaceGet",err,error,*999)

    IF(ALLOCATED(surfaceNodes)) CALL FlagError("Surface nodes array is already allocated.",err,error,*999)
    IF(.NOT.ASSOCIATED(regularMesh)) CALL FlagError("Regular mesh is not associated.",err,error,*999)
    
    IF(.NOT.ALLOCATED(regularMesh%numberOfElementsXi)) &
      & CALL FlagError("Regular mesh object does not have number of elements property specified.",err,error,*999)

    numberOfDimensions=SIZE(regularMesh%numberOfElementsXi)
    IF(numberOfDimensions==2) THEN
      numberOfElementsXi(1:2)=regularMesh%numberOfElementsXi(1:2)
      numberOfElementsXi(3)=1
    ELSE IF (numberOfDimensions==1) THEN
      numberOfElementsXi(1)=regularMesh%numberOfElementsXi(1)
      numberOfElementsXi(2)=1
      numberOfElementsXi(3)=1
    ELSE
      numberOfElementsXi=regularMesh%numberOfElementsXi
    ENDIF
    NULLIFY(basis)
    CALL GeneratedMeshRegular_BasisGet(regularMesh,meshComponent,basis,err,error,*999)
    !Node that simplex bases have an extra area coordinate so size of number_of_nodes_xic=numberOfDimensions+1
    numberOfNodesXi(1:numberOfDimensions)=basis%numberOfNodesXiC(1:numberOfDimensions)
    numberOfNodesXi(numberOfDimensions+1:3)=1

    ! build indices first (some of these are dummy arguments)
    CALL GeneratedMesh_RegularBuildNodeIndices(numberOfElementsXi,numberOfNodesXi,regularMesh%maximumExtent, &
      & totalNumberOfNodes,totalNumberOfElements,nodeIndices,elementIndices,delta,deltaI,err,error,*999)
    nodeCounter=0
    SELECT CASE(surfaceType)
    CASE(GENERATED_MESH_REGULAR_LEFT_SURFACE)
      ALLOCATE(surfaceNodes((SIZE(nodeIndices,2))*(SIZE(nodeIndices,3))),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate surface nodes array.",err,error,*999)
      DO k=1,SIZE(nodeIndices,3)
        DO j=1,SIZE(nodeIndices,2)
          nodeCounter=nodeCounter+1
          surfaceNodes(nodeCounter)=nodeIndices(1,j,k)
        ENDDO !j
      ENDDO !k
      normalXi=-1
    CASE(GENERATED_MESH_REGULAR_RIGHT_SURFACE)
      ALLOCATE(surfaceNodes((SIZE(nodeIndices,2))*(SIZE(nodeIndices,3))),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate surface nodes array.",err,error,*999)
      DO k=1,SIZE(nodeIndices,3)
        DO j=1,SIZE(nodeIndices,2)
          nodeCounter=nodeCounter+1
          surfaceNodes(nodeCounter)=nodeIndices(SIZE(nodeIndices,1),j,k)
        ENDDO !j
      ENDDO !k
      normalXi=1
    CASE(GENERATED_MESH_REGULAR_TOP_SURFACE)
      ALLOCATE(surfaceNodes((SIZE(nodeIndices,1))*(SIZE(nodeIndices,2))),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate surface nodes array.",err,error,*999)
      DO j=1,SIZE(nodeIndices,2)
        DO i=1,SIZE(nodeIndices,1)
          nodeCounter=nodeCounter+1
          surfaceNodes(nodeCounter)=nodeIndices(i,j,SIZE(nodeIndices,3))
        ENDDO !i
      ENDDO !j
      normalXi=3
    CASE(GENERATED_MESH_REGULAR_BOTTOM_SURFACE)
      ALLOCATE(surfaceNodes((SIZE(nodeIndices,1))*(SIZE(nodeIndices,2))),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate surface nodes array.",err,error,*999)
      DO j=1,SIZE(nodeIndices,2)
        DO i=1,SIZE(nodeIndices,1)
          nodeCounter=nodeCounter+1
          surfaceNodes(nodeCounter)=nodeIndices(i,j,1)
        ENDDO !i
      ENDDO !j
      normalXi=-3
    CASE(GENERATED_MESH_REGULAR_FRONT_SURFACE)
      ALLOCATE(surfaceNodes((SIZE(nodeIndices,1))*(SIZE(nodeIndices,3))),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate surface nodes array.",err,error,*999)
      DO j=1,SIZE(nodeIndices,3)
        DO i=1,SIZE(nodeIndices,1)
          nodeCounter=nodeCounter+1
          surfaceNodes(nodeCounter)=nodeIndices(i,1,j)
        ENDDO !i
      ENDDO !j
      normalXi=-2
    CASE(GENERATED_MESH_REGULAR_BACK_SURFACE)
      ALLOCATE(surfaceNodes((SIZE(nodeIndices,1))*(SIZE(nodeIndices,3))),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate NODES array.",err,error,*999)
      DO j=1,SIZE(nodeIndices,3)
        DO i=1,SIZE(nodeIndices,1)
          nodeCounter=nodeCounter+1
          surfaceNodes(nodeCounter)=nodeIndices(i,SIZE(nodeIndices,2),j)
        ENDDO !i
      ENDDO !j
      normalXi=2
    CASE DEFAULT
      localError="The specified surface type of "//TRIM(NumberToVString(surfaceType,"*",err,error))// &
        & " is invalid for a regular mesh."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !Now convert the component node numbering to user numbers if a mesh has multiple components
    DO nodeCounter=1,SIZE(surfaceNodes,1)
      SELECT CASE(basis%TYPE)
      CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
        CALL GeneratedMesh_RegularComponentNodeToUserNumber(regularMesh%generatedMesh,meshComponent, &
          & surfaceNodes(nodeCounter),nodeUserNumber,err,error,*999)
        surfaceNodes(nodeCounter)=nodeUserNumber
      CASE(BASIS_SIMPLEX_TYPE)
        surfaceNodes(nodeCounter)=GeneratedMesh_ComponentNodeToUserNumber(regularMesh%generatedMesh,meshComponent, &
          & surfaceNodes(nodeCounter),err,error)
        IF(err/=0) GOTO 999
      CASE DEFAULT
        CALL FlagError("The basis type of "//TRIM(NumberToVString(basis%TYPE,"*",err,error))// &
          & " is not implemented when getting a regular mesh surface.",err,error,*999)
      END SELECT
    ENDDO !nodeCounter

    EXITS("GeneratedMesh_RegularSurfaceGet")
    RETURN
999 ERRORSEXITS("GeneratedMesh_RegularSurfaceGet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_RegularSurfaceGet

  !
  !================================================================================================================================
  !

  !>Provides an easy way to grab surfaces for boundary condition assignment
  SUBROUTINE GeneratedMesh_CylinderSurfaceGet(cylinderMesh,meshComponent,surfaceType,surfaceNodes,normalXi,err,error,*)

    ! Argument variables
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh !<A pointer to the cylinder mesh object
    INTEGER(INTG), INTENT(IN) :: meshComponent !<The mesh component to get the surface for.
    INTEGER(INTG), INTENT(IN) :: surfaceType !<A constant identifying the type of surface to get \see GeneratedMeshRoutines_GeneratedMeshCylinderSurfaces,GeneratedMeshRoutines
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: surfaceNodes(:) !<On exit, contains the list of nodes belonging to the surface
    INTEGER(INTG), INTENT(OUT) :: normalXi !<On exit, contains the xi direction of the outward pointing normal of the surface
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    ! Local variables
    INTEGER(INTG) :: nodeCounter,i, j, k
    INTEGER(INTG) :: numberOfElementsXi(3) !Specified number of elements in each xi direction
    INTEGER(INTG) :: numberOfNodesXi(3) ! Number of nodes per element in each xi direction (basis property)
    INTEGER(INTG) :: totalNumberOfNodes,totalNumberOfElements
    INTEGER(INTG), ALLOCATABLE :: nodeIndices(:,:,:),elementIndices(:,:,:)
    REAL(DP) :: delta(3),deltai(3)
    TYPE(BasisType), POINTER :: basis
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_CylinderSurfaceGet",err,error,*999)

    IF(ALLOCATED(surfaceNodes)) CALL FlagError("Surface nodes array is already allocated.",err,error,*999)
    IF(.NOT.ASSOCIATED(cylinderMesh)) CALL FlagError("Cylinder mesh is not associated.",err,error,*999)
    
    IF(.NOT.ALLOCATED(cylinderMesh%numberOfElementsXi))  &
      & CALL FlagError("Cylinder mesh object does not have number of elements property specified.",err,error,*999)
    
    numberOfElementsXi=cylinderMesh%numberOfElementsXi
    NULLIFY(basis)
    CALL GeneratedMeshCylinder_BasisGet(cylinderMesh,meshComponent,basis,err,error,*999)
    numberOfNodesXi=basis%numberOfNodesXiC
    ! build indices first (some of these are dummy arguments)
    CALL GeneratedMesh_CylinderBuildNodeIndices(numberOfElementsXi,numberOfNodesXi,cylinderMesh%cylinderExtent, &
      & totalNumberOfNodes,totalNumberOfElements,nodeIndices,elementIndices,delta,deltai,err,error,*999)
    nodeCounter=0
    SELECT CASE(surfaceType)
    CASE(GENERATED_MESH_CYLINDER_INNER_SURFACE)
      ALLOCATE(surfaceNodes((SIZE(nodeIndices,2))*(SIZE(nodeIndices,3))),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate surface nodes array.",err,error,*999)
      DO k=1,SIZE(nodeIndices,3)
        DO j=1,SIZE(nodeIndices,2)
          nodeCounter=nodeCounter+1
          surfaceNodes(nodeCounter)=GeneratedMesh_ComponentNodeToUserNumber(cylinderMesh%generatedMesh,meshComponent, &
            & nodeIndices(1,j,k),err,error)
        ENDDO !j
      ENDDO !k
      normalXi=-1
    CASE(GENERATED_MESH_CYLINDER_OUTER_SURFACE)
      ALLOCATE(surfaceNodes((SIZE(nodeIndices,2))*(SIZE(nodeIndices,3))),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate surface nodes array.",err,error,*999)
      DO k=1,SIZE(nodeIndices,3)
        DO j=1,SIZE(nodeIndices,2)
          nodeCounter=nodeCounter+1
          surfaceNodes(nodeCounter)=GeneratedMesh_ComponentNodeToUserNumber(cylinderMesh%generatedMesh,meshComponent, &
            & nodeIndices(SIZE(nodeIndices,1),j,k),err,error)
        ENDDO !j
      ENDDO !k
      normalXi=1
    CASE(GENERATED_MESH_CYLINDER_TOP_SURFACE)
      ALLOCATE(surfaceNodes((SIZE(nodeIndices,1))*(SIZE(nodeIndices,2))),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate surface nodes array.",err,error,*999)
      DO j=1,SIZE(nodeIndices,2)
        DO i=1,SIZE(nodeIndices,1)
          nodeCounter=nodeCounter+1
          surfaceNodes(nodeCounter)=GeneratedMesh_ComponentNodeToUserNumber(cylinderMesh%generatedMesh,meshComponent, &
            & nodeIndices(i,j,SIZE(nodeIndices,3)),err,error)
        ENDDO !i
      ENDDO !j
      normalXi=3
    CASE(GENERATED_MESH_CYLINDER_BOTTOM_SURFACE)
      ALLOCATE(surfaceNodes((SIZE(nodeIndices,1))*(SIZE(nodeIndices,2))),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate surface nodes array.",err,error,*999)
      DO j=1,SIZE(nodeIndices,2)
        DO i=1,SIZE(nodeIndices,1)
          nodeCounter=nodeCounter+1
          surfaceNodes(nodeCounter)=GeneratedMesh_ComponentNodeToUserNumber(cylinderMesh%generatedMesh,meshComponent, &
            & nodeIndices(i,j,1),err,error)
        ENDDO !i
      ENDDO !j
      normalXi=-3
    CASE DEFAULT
      localError="The specified surface type of "//TRIM(NumberToVString(surfaceType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("GeneratedMesh_CylinderSurfaceGet")
    RETURN
999 ERRORSEXITS("GeneratedMesh_CylinderSurfaceGet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_CylinderSurfaceGet
  !
  !================================================================================================================================
  !

  !>Provides an easy way to grab surfaces for boundary condition assignment
  SUBROUTINE GeneratedMesh_EllipsoidSurfaceGet(ellipsoidMesh,meshComponent,surfaceType,surfaceNodes,normalXi,err,error,*)

    ! Argument variables
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ellipsoidMesh !<A pointer to the ellipsoid mesh object
    INTEGER(INTG), INTENT(IN) :: meshComponent !<The mesh component to get the surface for.
    INTEGER(INTG), INTENT(IN) :: surfaceType !<A constant identifying the type of surface to get \see GeneratedMeshRoutines_GeneratedMeshEllipsoidSurfaces,GeneratedMeshRoutines
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: surfaceNodes(:) !<On exit, contains the list of nodes belonging to the surface
    INTEGER(INTG), INTENT(OUT) :: normalXi !<On exit, contains the xi direction of the outward pointing normal of the surface
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    ! Local variables
    TYPE(BasisType), POINTER :: basis
    INTEGER(INTG),ALLOCATABLE :: nodeIndices(:,:,:),elementIndices(:,:,:),cornerNodes(:,:,:)
    INTEGER(INTG) :: numberOfElementsXi(3) !Specified number of elements in each xi direction
    INTEGER(INTG) :: numberOfNodesXi(3) ! Number of nodes per element in each xi direction (basis property)
    INTEGER(INTG) :: totalNumberOfNodes,totalNumberOfElements
    REAL(DP) :: delta(3),deltai(3)
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: nodeCounter,i, j, k

    ENTERS("GeneratedMesh_EllipsoidSurfaceGet",err,error,*999)

    IF(ALLOCATED(surfaceNodes)) CALL FlagError("Surface nodes array is already allocated.",err,error,*999)
    IF(.NOT.ASSOCIATED(ellipsoidMesh)) CALL FlagError("Ellipsoid mesh is not associated.",err,error,*999)
    
    IF(.NOT.ALLOCATED(ellipsoidMesh%numberOfElementsXi)) &
      & CALL FlagError("Ellipsoid mesh object does not have number of elements property specified.",err,error,*999)

    numberOfElementsXi=ellipsoidMesh%numberOfElementsXi
    NULLIFY(basis)
    CALL GeneratedMeshEllipsoid_BasisGet(ellipsoidMesh,meshComponent,basis,err,error,*999)
    numberOfNodesXi=basis%numberOfNodesXiC
    ! build indices first (some of these are dummy arguments)
    CALL GeneratedMesh_EllipsoidBuildNodeIndices(numberOfElementsXi,numberOfNodesXi, &
      & ellipsoidMesh%ellipsoidExtent,totalNumberOfNodes,totalNumberOfElements,nodeIndices, &
      & cornerNodes,elementIndices,delta,deltai,err,error,*999)
    nodeCounter=0

    SELECT CASE(surfaceType)
    CASE(GENERATED_MESH_ELLIPSOID_INNER_SURFACE)
      ALLOCATE(surfaceNodes((SIZE(nodeIndices,1))*(SIZE(nodeIndices,2)-1)+1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate surface nodes array.",err,error,*999)
      j=1
      i=1
      nodeCounter=nodeCounter+1
      surfaceNodes(nodeCounter)=GeneratedMesh_ComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,meshComponent, &
        nodeIndices(i,j,1),err,error)
      DO j=2,SIZE(nodeIndices,2)
        DO i=1, SIZE(nodeIndices,1)
          nodeCounter=nodeCounter+1
          IF (nodeIndices(i,j,1)/=0) THEN
            surfaceNodes(nodeCounter)=GeneratedMesh_ComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,meshComponent, &
              & nodeIndices(i,j,1),err,error)
          ELSE
            nodeCounter=nodeCounter-1
          ENDIF
        ENDDO !i
      ENDDO !j
      normalXi=-3      
    CASE(GENERATED_MESH_ELLIPSOID_OUTER_SURFACE)
      ALLOCATE(surfaceNodes((SIZE(nodeIndices,1))*(SIZE(nodeIndices,2)-1)+1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate surface nodes array.",err,error,*999)
      j=1
      i=1
      nodeCounter=nodeCounter+1
      surfaceNodes(nodeCounter)=GeneratedMesh_ComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,meshComponent, &
        & nodeIndices(i,j,SIZE(nodeIndices,3)),err,error)
      DO j=2,SIZE(nodeIndices,2)
        DO i=1, SIZE(nodeIndices,1)
          nodeCounter=nodeCounter+1
          IF (nodeIndices(i,j,SIZE(nodeIndices,3))/=0) THEN
            surfaceNodes(nodeCounter)=GeneratedMesh_ComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,meshComponent, &
              & nodeIndices(i,j,SIZE(nodeIndices,3)),err,error)
          ELSE
            nodeCounter=nodeCounter-1
          ENDIF
        ENDDO !i
      ENDDO !j
      normalXi=3      
    CASE(GENERATED_MESH_ELLIPSOID_TOP_SURFACE)
      ALLOCATE(surfaceNodes((SIZE(nodeIndices,1))*(SIZE(nodeIndices,3))),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate surface nodes array.",err,error,*999)
      DO k=1,SIZE(nodeIndices,3)
        DO i=1, SIZE(nodeIndices,1)
          nodeCounter=nodeCounter+1
          IF (nodeIndices(i,SIZE(nodeIndices,2),k)/=0) THEN
            surfaceNodes(nodeCounter)=GeneratedMesh_ComponentNodeToUserNumber(ellipsoidMesh%generatedMesh,meshComponent, &
              & nodeIndices(i,SIZE(nodeIndices,2),k),err,error)
          ELSE
            nodeCounter=nodeCounter-1
          ENDIF
        ENDDO !i
      ENDDO !k
      normalXi=2
    CASE DEFAULT
      localError="The specified surface type of "//TRIM(NumberToVString(surfaceType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("GeneratedMesh_EllipsoidSurfaceGet")
    RETURN
999 ERRORSEXITS("GeneratedMesh_EllipsoidSurfaceGet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_EllipsoidSurfaceGet

  !
  !================================================================================================================================
  !

  !>Calculates the mesh topology information for a given regular mesh (Not to be called by user)
  SUBROUTINE GeneratedMesh_RegularBuildNodeIndices(numberOfElementsXi,numberOfNodesXic,maximumExtent, &
    & totalNumberOfNodes,totalNumberOfElements,nodeIndices,elementIndices,delta,deltai,err,error,*)

    !Argument variables
    INTEGER(INTG),INTENT(IN) :: numberOfElementsXi(3) !<Specified number of elements in each xi direction
    INTEGER(INTG),INTENT(IN) :: numberOfNodesXic(3) !<Number of nodes per element in each xi direction for this basis
    REAL(DP),INTENT(IN) :: maximumExtent(3)         !<width, length and height of regular mesh
    INTEGER(INTG),INTENT(OUT) :: totalNumberOfNodes    !<On exit, contains total number of nodes in regular mesh component
    INTEGER(INTG),INTENT(OUT) :: totalNumberOfElements !<On exit, contains total number of elements in regular mesh
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: nodeIndices(:,:,:)  !<Mapping array to find a node number for a given (x,y,z)
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: elementIndices(:,:,:)  !<Mapping array to find an element number for a given (x,y,z)
    REAL(DP),INTENT(OUT) :: delta(3)  !<Step sizes in each of (x,y,z) for elements
    REAL(DP),INTENT(OUT) :: deltai(3) !<Step sizes in each of (x,y,z) for node (identical to delta if 2 nodes per xi direction)
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: xiIdx,elementIdx1,elementIdx2,elementIdx3,localNodeIdx1,localNodeIdx2,localNodeIdx3,localNodeNumber, &
      & localElementNumber
    INTEGER(INTG) :: totalNumberOfNodesXi(3) !<Total number of nodes per element in each xi direction for this basis

    ENTERS("GeneratedMesh_RegularBuildNodeIndices",err,error,*999)

    IF(ALLOCATED(nodeIndices)) CALL FlagError("Node indices array is already allocated.",err,error,*999)
    IF(ALLOCATED(elementIndices)) CALL FlagError("Element indices array is already allocated.",err,error,*999)
    
    !Calculate delta and deltai
    delta(1)=maximumExtent(1)/numberOfElementsXi(1)
    delta(2)=maximumExtent(2)/numberOfElementsXi(2)
    delta(3)=maximumExtent(3)/numberOfElementsXi(3)
    DO xiIdx=1,3
      IF(numberOfNodesXic(xiIdx)>1) deltai(xiIdx)=delta(xiIdx)/(numberOfNodesXic(xiIdx)-1)
    ENDDO !xiIdx

    !Calculate total elements and nodes
    DO xiIdx=1,3
      totalNumberOfNodesXi(xiIdx)=(numberOfNodesXic(xiIdx)-1)*numberOfElementsXi(xiIdx)+1
    ENDDO !xiIdx
    totalNumberOfElements=PRODUCT(numberOfElementsXi)

    !Calculate nodeIndices first
    ALLOCATE(nodeIndices(totalNumberOfNodesXi(1),totalNumberOfNodesXi(2),totalNumberOfNodesXi(3)),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate nodeIndices array.",err,error,*999)
    localNodeNumber=0
    DO localNodeIdx3=1,totalNumberOfNodesXi(3)
      DO localNodeIdx2=1,totalNumberOfNodesXi(2)
        DO localNodeIdx1=1,totalNumberOfNodesXi(1)
          localNodeNumber=localNodeNumber+1
          nodeIndices(localNodeIdx1,localNodeIdx2,localNodeIdx3)=localNodeNumber
        ENDDO !localNodeIdx1
      ENDDO !localNodeIdx2
    ENDDO !localNodeIdx3
    totalNumberOfNodes=localNodeNumber

    !Now do elementIndices
    ALLOCATE(elementIndices(numberOfElementsXi(1),numberOfElementsXi(2),numberOfElementsXi(3)),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate elementIndices array.",err,error,*999)
    localElementNumber=0
    DO elementIdx3=1,numberOfElementsXi(3)
      DO elementIdx2=1,numberOfElementsXi(2)
        DO elementIdx1=1,numberOfElementsXi(1)
          localElementNumber=localElementNumber+1
          elementIndices(elementIdx1,elementIdx2,elementIdx3)=localElementNumber
        ENDDO !elementIdx1
      ENDDO !elementIdx2
    ENDDO !elementIdx3
    totalNumberOfElements=localElementNumber

    EXITS("GeneratedMesh_RegularBuildNodeIndices")
    RETURN
999 ERRORSEXITS("GeneratedMesh_RegularBuildNodeIndices",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_RegularBuildNodeIndices

  !
  !================================================================================================================================
  !

  !>Calculates the mesh topology information for a given cylinder (Not to be called by user)
  SUBROUTINE GeneratedMesh_CylinderBuildNodeIndices(numberOfElementsXi,numberOfNodesXic,cylinderExtent, &
      & totalNumberOfNodes,totalNumberOfElements,nodeIndices,elementIndices,delta,deltai,err,error,*)

    ! Argument variables
    INTEGER(INTG),INTENT(IN) :: numberOfElementsXi(3) !<Specified number of elements in each xi direction
    INTEGER(INTG),INTENT(IN) :: numberOfNodesXic(3) !<Number of nodes per element in each xi direction (basis property)
    REAL(DP),INTENT(IN) :: cylinderExtent(3)         !<inner & outer radii and height of cylinder
    INTEGER(INTG),INTENT(OUT) :: totalNumberOfNodes    !<On exit, contains total number of nodes in cylinder mesh
    INTEGER(INTG),INTENT(OUT) :: totalNumberOfElements !<On exit, contains total number of elements in cylinder mesh
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: nodeIndices(:,:,:)  !<Mapping array to find a node number for a given (r,theta,z)
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: elementIndices(:,:,:)  !<Mapping array to find an element number for a given (r,theta,z)
    REAL(DP),INTENT(OUT) :: delta(3)  !<Step sizes in each of (r,theta,z) for elements
    REAL(DP),INTENT(OUT) :: deltai(3) !<Step sizes in each of (r,theta,z) for node (identical to delta if 2 nodes per xi direction)
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING) :: error !<The error string

    ! Local variables
    INTEGER(INTG) :: xiIdx,elementIdx1,elementIdx2,elementIdx3,localNodeIdx1,localNodeIdx2,localNodeIdx3,localNodeNumber, &
      & localElementNumber
    INTEGER(INTG) :: totalNumberOfNodesXi(3) ! total number of nodes in each xi direction

    ENTERS("GeneratedMesh_CylinderBuildNodeIndices",err,error,*999)

    IF(ALLOCATED(nodeIndices)) CALL FlagError("Node indices array is already allocated.",err,error,*999)
    IF(ALLOCATED(elementIndices)) CALL FlagError("Element indices array is already allocated.",err,error,*999)
    
    !Calculate delta and deltai
    delta(1)=(cylinderExtent(2)-cylinderExtent(1))/numberOfElementsXi(1)
    delta(2)=TWOPI/numberOfElementsXi(2)
    delta(3)=cylinderExtent(3)/numberOfElementsXi(3)
    DO xiIdx=1,3
      deltai(xiIdx)=delta(xiIdx)/(numberOfNodesXic(xiIdx)-1)
    ENDDO !xiIdx

    !Calculate total elements and nodes
    DO xiIdx=1,3
      totalNumberOfNodesXi(xiIdx)=(numberOfNodesXic(xiIdx)-1)*numberOfElementsXi(xiIdx)+1
    ENDDO !xiIdx
    totalNumberOfNodesXi(2)=totalNumberOfNodesXi(2)-1 ! theta loops around so slightly different
    !totalNumberOfElements=PRODUCT(numberOfElementsXi)
    
    !Calculate nodeIndices first
    ALLOCATE(nodeIndices(totalNumberOfNodesXi(1),totalNumberOfNodesXi(2),totalNumberOfNodesXi(3)),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate nodeIndices array.",err,error,*999)
    localNodeNumber=0
    DO localNodeIdx3=1,totalNumberOfNodesXi(3)
      DO localNodeIdx2=1,totalNumberOfNodesXi(2)
        DO localNodeIdx1=1,totalNumberOfNodesXi(1)
          localNodeNumber=localNodeNumber+1
          nodeIndices(localNodeIdx1,localNodeIdx2,localNodeIdx3)=localNodeNumber
        ENDDO ! localNodeIdx1
      ENDDO ! localNodeIdx2
    ENDDO ! localNodeIdx3
    totalNumberOfNodes=localNodeNumber
    
    !Now do elementIndices
    ALLOCATE(elementIndices(numberOfElementsXi(1),numberOfElementsXi(2),numberOfElementsXi(3)),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate elementIndices array.",err,error,*999)
    localElementNumber=0
    DO elementIdx3=1,numberOfElementsXi(3)
      DO elementIdx2=1,numberOfElementsXi(2)
        DO elementIdx1=1,numberOfElementsXi(1)
          localElementNumber=localElementNumber+1
          elementIndices(elementIdx1,elementIdx2,elementIdx3)=localElementNumber
        ENDDO !elementIdx1
      ENDDO !elementIdx2
    ENDDO !elementIdx3
    totalNumberOfElements=localElementNumber

    EXITS("GeneratedMesh_CylinderBuildNodeIndices")
    RETURN
999 ERRORSEXITS("GeneratedMesh_CylinderBuildNodeIndices",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_CylinderBuildNodeIndices

  !
  !================================================================================================================================
  !

  !>Calculate the mesh topology information for a given ellipsoid (Not to be called by user)
  SUBROUTINE GeneratedMesh_EllipsoidBuildNodeIndices(numberOfElementsXi,numberOfNodesXi,ellipsoidExtent, &
    & totalNumberOfNodes,totalNumberOfElements,nodeIndices,cornerNodes,elementIndices,delta,deltai,err,error,*)

    ! Argument variables
    INTEGER(INTG),INTENT(IN) :: numberOfElementsXi(3) !<Specified number of elements in each xi direction
    INTEGER(INTG),INTENT(IN) :: numberOfNodesXi(3) !<Number of nodes per element in each xi direction (basis property)
    REAL(DP),INTENT(IN) :: ellipsoidExtent(4)         !< long axis, short axis, wall thickness, top angle
    INTEGER(INTG),INTENT(OUT) :: totalNumberOfNodes    !<On exit, contains total number of nodes in ellipsoid mesh
    INTEGER(INTG),INTENT(OUT) :: totalNumberOfElements !<On exit, contains total number of elements in ellipsoid mesh
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: nodeIndices(:,:,:)  !<Mapping array to find a node number for a given (r,theta,z)
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: cornerNodes(:,:,:) ! Returns the array of corner nodes numbered
    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: elementIndices(:,:,:)  !<Mapping array to find an element number for a given (r,theta,z)
    REAL(DP),INTENT(OUT) :: delta(3)  !<Step sizes in each of (r,theta,z) for elements
    REAL(DP),INTENT(OUT) :: deltai(3) !<Step sizes in each of (r,theta,z) for node (identical to delta if 2 nodes per xi direction)
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING) :: error !<The error string
    ! Local variables
    INTEGER(INTG) :: xiIdx,elementIdx1,elementIdx2,elementIdx3,localNodeIdx1,localNodeIdx2,localNodeIdx3, &
      & globalNodeIdx1,globalNodeIdx2,globalNodeIdx3,localNodeNumber,localElementNumber
    INTEGER(INTG) :: totalNumberOfNodesXi(3) ! total number of nodes in each xi direction

    ENTERS("GeneratedMesh_EllipsoidBuildNodeIndices",err,error,*999)

    IF(ALLOCATED(nodeIndices)) CALL FlagError("Node indices array is already allocated.",err,error,*999)
    IF(ALLOCATED(cornerNodes)) CALL FlagError("Corner nodes array is already allocated.",err,error,*999)
    IF(ALLOCATED(elementIndices)) CALL FlagError("Element indices array is already allocated.",err,error,*999)
    
    !Calculate delta and deltai
    delta(1)=TWOPI/numberOfElementsXi(1)
    delta(2)=(PI-ellipsoidExtent(4))/numberOfElementsXi(2)
    delta(3)=ellipsoidExtent(3)/numberOfElementsXi(3)
    DO xiIdx=1,3
      deltai(xiIdx)=delta(xiIdx)/(numberOfNodesXi(xiIdx)-1)
    ENDDO !xiIdx
    
    !Calculate total elements and nodes
    DO xiIdx=1,3
      totalNumberOfNodesXi(xiIdx)=(numberOfNodesXi(xiIdx)-1)*numberOfElementsXi(xiIdx)+1
    ENDDO !xiIdx
    totalNumberOfNodesXi(1)=totalNumberOfNodesXi(1)-1 ! circumferential loops around so slightly different
    totalNumberOfElements=PRODUCT(numberOfElementsXi)

    !Calculate nodeIndices first
    ALLOCATE(nodeIndices(totalNumberOfNodesXi(1),totalNumberOfNodesXi(2),totalNumberOfNodesXi(3)),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate node indices array.",err,error,*999)
    ALLOCATE(cornerNodes(numberOfElementsXi(1),numberOfElementsXi(2)+1,numberOfElementsXi(3)+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate corner nodes array.",err,error,*999)
    
    !localNodeNumber: node number inside element in certain direction
    !localElementNumber: element number in certain direction
    !globalNodeIdx: global node number in certain direction
    !localNodeNumber: Node counter
    !Due to one more corner node than elements in transmural direction, first shell is taken separatly
    localNodeNumber=0
    elementIdx3=1
    localNodeIdx3=1
    !Due to one more corner node than elements in longitudinal direction, apex elements are taken separatly
    elementIdx2=1
    localNodeIdx2=1
    elementIdx1=1
    localNodeIdx1=1
    !Apex nodes
    localNodeNumber=localNodeNumber+1
    globalNodeIdx3=1
    globalNodeIdx2=1
    globalNodeIdx1=1
    nodeIndices(globalNodeIdx1,globalNodeIdx2,globalNodeIdx3)=localNodeNumber
    cornerNodes(elementIdx1,elementIdx2,elementIdx3)=localNodeNumber
    DO elementIdx2=1,numberOfElementsXi(2)
      DO localNodeIdx2=2,(numberOfNodesXi(2))
        globalNodeIdx2=globalNodeIdx2+1
        globalNodeIdx1=0
        DO elementIdx1=1,numberOfElementsXi(1)
          DO localNodeIdx1=1,(numberOfNodesXi(1)-1)
            globalNodeIdx1=globalNodeIdx1+1
            localNodeNumber=localNodeNumber+1
            nodeIndices(globalNodeIdx1,globalNodeIdx2,globalNodeIdx3)=localNodeNumber
            IF ((localNodeIdx1==1).AND.(localNodeIdx2==numberOfNodesXi(2))) THEN
              cornerNodes(elementIdx1,elementIdx2+1,elementIdx3)=localNodeNumber
            ENDIF
          ENDDO !localNodeIdx1
        ENDDO !elementIdx1
      ENDDO !localNodeIdx2
    ENDDO !elementIdx2
    DO elementIdx3=1,numberOfElementsXi(3)
      DO localNodeIdx3=2,numberOfNodesXi(3)
        elementIdx2=1
        localNodeIdx2=1
        elementIdx1=1
        localNodeIdx1=1
        !apex nodes
        localNodeNumber=localNodeNumber+1
        globalNodeIdx3=globalNodeIdx3+1
        globalNodeIdx2=1
        globalNodeIdx1=1
        nodeIndices(globalNodeIdx1,globalNodeIdx2,globalNodeIdx3)=localNodeNumber
        IF (localNodeIdx3==numberOfNodesXi(3)) cornerNodes(elementIdx1,elementIdx2,elementIdx3+1)=localNodeNumber
        DO elementIdx2=1,numberOfElementsXi(2)
          DO localNodeIdx2=2,(numberOfNodesXi(2))
            globalNodeIdx2=globalNodeIdx2+1
            globalNodeIdx1=0
            DO elementIdx1=1,numberOfElementsXi(1)
              DO localNodeIdx1=1,(numberOfNodesXi(1)-1)
                globalNodeIdx1=globalNodeIdx1+1
                localNodeNumber=localNodeNumber+1
                nodeIndices(globalNodeIdx1,globalNodeIdx2,globalNodeIdx3)=localNodeNumber
                IF ((localNodeIdx1==1).AND.(localNodeIdx3==numberOfNodesXi(3)).AND.(localNodeIdx2==numberOfNodesXi(2))) THEN
                  cornerNodes(elementIdx1,elementIdx2+1,elementIdx3+1)=localNodeNumber
                ENDIF
              ENDDO !localNodeIdx1
            ENDDO !elementIdx1
          ENDDO !localNodeIdx2
        ENDDO !elementIdx2
      ENDDO !localNodeIdx3
    ENDDO !elementIdx3
    totalNumberOfNodes=localNodeNumber

    !Now do elementIndices
    ALLOCATE(elementIndices(numberOfElementsXi(1),numberOfElementsXi(2),numberOfElementsXi(3)),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate elementIndices array.",err,error,*999)
    localElementNumber=0
    DO elementIdx3=1,numberOfElementsXi(3)
      DO elementIdx2=1,numberOfElementsXi(2)
        DO elementIdx1=1,numberOfElementsXi(1)
          localElementNumber=localElementNumber+1
          elementIndices(elementIdx1,elementIdx2,elementIdx3)=localElementNumber
        ENDDO !elementIdx1
      ENDDO !elementIdx2
    ENDDO !elementIdx3
    totalNumberOfElements=localElementNumber
    
    EXITS("GeneratedMesh_EllipsoidBuildNodeIndices")
    RETURN
999 ERRORSEXITS("GeneratedMesh_EllipsoidBuildNodeIndices",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_EllipsoidBuildNodeIndices

  !
  !================================================================================================================================
  !

  !>Calculates the user node numbers for an array of nodes numbered using one basis
  SUBROUTINE GeneratedMeh_ComponentNodesToUserNumbers(generatedMesh,basisIndex,nodeComponentNumbers, &
    & nodeUserNumbers,err,error,*)
    
    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh object
    INTEGER(INTG),INTENT(IN) :: basisIndex !<The number of the basis being used
    INTEGER(INTG),INTENT(IN) :: nodeComponentNumbers(:) !<The node numbers for this component basis
    INTEGER(INTG),INTENT(INOUT) :: nodeUserNumbers(:) !<On return, the corresponding user numbers
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: nodeIdx

    ENTERS("GeneratedMeh_ComponentNodesToUserNumbers",err,error,*999)

    IF(SIZE(nodeUserNumbers)/=SIZE(nodeComponentNumbers)) &
      & CALL FlagError("Node component numbers and node user numbers arrays have different sizes.",err,error,*999)
    
    DO nodeIdx=1,SIZE(nodeComponentNumbers,1)
      nodeUserNumbers(nodeIdx)=GeneratedMesh_ComponentNodeToUserNumber(generatedMesh,basisIndex,nodeComponentNumbers(nodeIdx), &
        & err,error)
    ENDDO !nodeIdx

    EXITS("GeneratedMeh_ComponentNodesToUserNumbers")
    RETURN
999 ERRORSEXITS("GeneratedMeh_ComponentNodesToUserNumbers",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMeh_ComponentNodesToUserNumbers

  !
  !================================================================================================================================
  !

  !>Calculates the user node number for a node numbered using one basis
  FUNCTION GeneratedMesh_ComponentNodeToUserNumber(generatedMesh,basisIndex,nodeComponentNumber,err,error)
    
    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh  !<A pointer to the generated mesh object
    INTEGER(INTG),INTENT(IN) :: basisIndex !<The number of the basis being used
    INTEGER(INTG),INTENT(IN) :: nodeComponentNumber!<The node number for this component basis
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING) :: error !<The error string
    !Function variable
    INTEGER(INTG) :: GeneratedMesh_ComponentNodeToUserNumber !<On return, the corresponding user node number
    !Local variables
    INTEGER(INTG) :: numberOfBases,numberOfDimensions,basisIdx,xiIdx,remainder,remainder2,tempTerm,numberOfCornerNodes, &
      & nodeOffset,basisNumberOfNodes
    INTEGER(INTG) :: position(3),position2(3),cornerNodeFactor(3),basisNodeFactor(3),basisElementFactor(3), &
      & numberOfPreviousCorners,step
    INTEGER(INTG), POINTER :: numberOfElementsXi(:)
    LOGICAL :: cornerNode,finishedCount
    TYPE(BasisType), POINTER :: basis
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ellipsoidMesh
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_ComponentNodeToUserNumber",err,error,*999)

    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is not associated",err,error,*999)
      
    numberOfCornerNodes=1
    remainder=nodeComponentNumber-1 !use zero based numbering
    remainder2=nodeComponentNumber-1
    GeneratedMesh_ComponentNodeToUserNumber=0
    position=0
    position2=0

    SELECT CASE(generatedMesh%generatedType)
    CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
      NULLIFY(regularMesh)
      CALL GeneratedMesh_RegularMeshGet(generatedMesh,regularMesh,err,error,*999)
      numberOfBases=SIZE(regularMesh%bases)
      numberOfDimensions=regularMesh%meshDimension
      numberOfElementsXi=>regularMesh%numberOfElementsXi
    CASE(GENERATED_MESH_POLAR_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
      NULLIFY(cylinderMesh)
      CALL GeneratedMesh_CylinderMeshGet(generatedMesh,cylinderMesh,err,error,*999)
      numberOfBases=SIZE(cylinderMesh%bases)
      numberOfDimensions=cylinderMesh%meshDimension
      numberOfElementsXi=>cylinderMesh%numberOfElementsXi
    CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
      NULLIFY(ellipsoidMesh)
      CALL GeneratedMesh_EllipsoidMeshGet(generatedMesh,ellipsoidMesh,err,error,*999)
      numberOfBases=SIZE(ellipsoidMesh%bases)
      numberOfDimensions=ellipsoidMesh%meshDimension
      numberOfElementsXi=>ellipsoidMesh%numberOfElementsXi
    CASE DEFAULT
      localError="The generated mesh generated type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    IF(basisIndex<1.OR.basisIndex>numberOfBases) THEN
      localError="Mesh component must be >=1 and <= "//(NumberToVString(numberOfBases,"*",err,error))// &
        & " but it is "//(NumberToVString(basisIndex,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    IF(numberOfBases==1) THEN
      !If is the only basis, don't do anything
      GeneratedMesh_ComponentNodeToUserNumber=nodeComponentNumber
    ELSE
      tempTerm=1
      numberOfCornerNodes=1
      DO xiIdx=1,numberOfDimensions
        numberOfCornerNodes=numberOfCornerNodes*(numberOfElementsXi(xiIdx)+1)
        cornerNodeFactor(xiIdx)=1
        IF(xiIdx>1) THEN
          tempTerm=tempTerm*(numberOfElementsXi(xiIdx-1)+1)
          cornerNodeFactor(xiIdx)=cornerNodeFactor(xiIdx)*tempTerm
        ENDIF
      ENDDO !xiIdx
      !Adjust for other mesh types
      IF(generatedMesh%generatedType==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
        cornerNodeFactor(3)=cornerNodeFactor(3)-numberOfElementsXi(1)-1
        numberOfCornerNodes=numberOfCornerNodes-(numberOfElementsXi(1)+1)*(numberOfElementsXi(3)+1)
      ELSE IF(generatedMesh%generatedType==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
        cornerNodeFactor(3)=cornerNodeFactor(3)-numberOfElementsXi(1)-numberOfElementsXi(2)
        cornerNodeFactor(2)=cornerNodeFactor(2)-1
        numberOfCornerNodes=numberOfCornerNodes-(numberOfElementsXi(2)+1)*(numberOfElementsXi(3)+1)- &
          & (numberOfElementsXi(1)-1)*(numberOfElementsXi(3)+1)
      ENDIF
      nodeOffset=numberOfCornerNodes
      IF(generatedMesh%generatedType==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
        !Every second mesh component is the collapsed node version
        step=2
      ELSE
        step=1
      ENDIF
      DO basisIdx=1,basisIndex-1,step
        NULLIFY(basis)
        IF(generatedMesh%generatedType==GENERATED_MESH_REGULAR_MESH_TYPE) THEN
          CALL GeneratedMeshRegular_BasisGet(regularMesh,basisIdx,basis,err,error,*999)
        ELSE IF(generatedMesh%generatedType==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
          CALL GeneratedMeshCylinder_BasisGet(cylinderMesh,basisIdx,basis,err,error,*999)
        ELSE IF(generatedMesh%generatedType==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
          CALL GeneratedMeshEllipsoid_BasisGet(ellipsoidMesh,basisIdx,basis,err,error,*999)
        ELSE
          CALL FlagError("Could not get a basis.",err,error,*999)
        ENDIF
        basisNumberOfNodes=1
        DO xiIdx=1,numberOfDimensions
          basisNumberOfNodes=basisNumberOfNodes*(numberOfElementsXi(xiIdx)*(basis%numberOfNodesXiC(xiIdx)-1)+1)
        ENDDO !xiIdx
        !Adjust for other mesh types
        IF(generatedMesh%generatedType==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
          basisNumberOfNodes=basisNumberOfNodes-(numberOfElementsXi(1)+1)*(basis%numberOfNodesXiC(1)-1)* &
            & (numberOfElementsXi(3)+1)*(basis%numberOfNodesXiC(3)-1)
        ELSE IF(generatedMesh%generatedType==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
          basisNumberOfNodes=basisNumberOfNodes-(numberOfElementsXi(2)*(basis%numberOfNodesXiC(2)-1)+1)* &
            & (numberOfElementsXi(3)*(basis%numberOfNodesXiC(3)-1)+1)- &
            & (numberOfElementsXi(1)*(basis%numberOfNodesXiC(1)-1)-1)* &
            & (numberOfElementsXi(3)*(basis%numberOfNodesXiC(3)-1)+1)
        ENDIF
        nodeOffset=nodeOffset+basisNumberOfNodes-numberOfCornerNodes
      ENDDO !basisIdx
      NULLIFY(basis)
      IF(generatedMesh%generatedType==GENERATED_MESH_REGULAR_MESH_TYPE) THEN
        CALL GeneratedMeshRegular_BasisGet(regularMesh,basisIndex,basis,err,error,*999)
      ELSE IF(generatedMesh%generatedType==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
        CALL GeneratedMeshCylinder_BasisGet(cylinderMesh,basisIndex,basis,err,error,*999)
      ELSE IF(generatedMesh%generatedType==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
        CALL GeneratedMeshEllipsoid_BasisGet(ellipsoidMesh,basisIndex,basis,err,error,*999)
      ELSE
        CALL FlagError("Could not get a basis.",err,error,*999)
      ENDIF
      tempTerm=1
      DO xiIdx=1,numberOfDimensions
        basisNodeFactor(xiIdx)=1
        basisElementFactor(xiIdx)=basis%numberOfNodesXiC(xiIdx)-1
        IF(xiIdx>1) THEN
          tempTerm=tempTerm*((basis%numberOfNodesXiC(xiIdx-1)-1)*numberOfElementsXi(xiIdx-1)+1)
          basisNodeFactor(xiIdx)=basisNodeFactor(xiIdx)*tempTerm
          basisElementFactor(xiIdx)=basisElementFactor(xiIdx)*tempTerm
        ENDIF
      ENDDO !xiIdx
      !Adjust for other mesh types
      IF(generatedMesh%generatedType==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
        !subtract nodes along line where y wraps around
        basisNodeFactor(3)=basisNodeFactor(3)-numberOfElementsXi(1)*(basis%numberOfNodesXiC(1)-1)-1
        basisElementFactor(3)=basisElementFactor(3)-(numberOfElementsXi(1)* &
          & (basis%numberOfNodesXiC(1)-1)+1)*(basis%numberOfNodesXiC(3)-1)
      ELSE IF(generatedMesh%generatedType==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
        !subtract missing nodes at apex
        basisNodeFactor(3)=basisNodeFactor(3)-numberOfElementsXi(1)*(basis%numberOfNodesXiC(1)-1)+1
        basisElementFactor(3)=basisElementFactor(3)-(numberOfElementsXi(1)* &
          & (basis%numberOfNodesXiC(1)-1)+1)*(basis%numberOfNodesXiC(3)-1)
        !subtract nodes along line where x wraps around
        basisNodeFactor(3)=basisNodeFactor(3)-numberOfElementsXi(2)*(basis%numberOfNodesXiC(2)-1)-1
        basisElementFactor(3)=basisElementFactor(3)-(numberOfElementsXi(2)*(basis%numberOfNodesXiC(2)-1)-1)* &
          & (basis%numberOfNodesXiC(3)-1)
        basisNodeFactor(2)=basisNodeFactor(2)-1
        basisElementFactor(2)=basisElementFactor(2)-(basis%numberOfNodesXiC(2)-1)
      ENDIF
      !Work out if we have a corner node, otherwise add node numbers used by corners and
      !previous basis interpolations and subtract number of corner nodes used before the
      !given component node number to get the user number
      cornerNode=.TRUE.
      IF(numberOfDimensions>2) THEN
        position(3)=remainder/basisNodeFactor(3)
        position2(3)=remainder2/basisElementFactor(3)
        remainder=MOD(remainder,basisNodeFactor(3))
        remainder2=MOD(remainder2,basisElementFactor(3))
        IF(MOD(position(3),basis%numberOfNodesXiC(3)-1)/=0) cornerNode=.FALSE.
      ENDIF
      IF(numberOfDimensions>1) THEN
        IF(generatedMesh%generatedType==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
          !Need to account for missing nodes at apex
          IF(remainder>0) THEN
            remainder=remainder+numberOfElementsXi(1)*(basis%numberOfNodesXiC(1)-1)-1
            remainder2=remainder2+numberOfElementsXi(1)*(basis%numberOfNodesXiC(1)-1)-1
          ENDIF
        ENDIF
        position(2)=remainder/basisNodeFactor(2)
        position2(2)=remainder2/basisElementFactor(2)
        remainder=MOD(remainder,basisNodeFactor(2))
        remainder2=MOD(remainder2,basisElementFactor(2))
        IF(MOD(position(2),basis%numberOfNodesXiC(2)-1)/=0) cornerNode=.FALSE.
      ENDIF
      position(1)=remainder/basisNodeFactor(1)
      position2(1)=remainder2/basisElementFactor(1)
      IF(MOD(position(1),basis%numberOfNodesXiC(1)-1)/=0) cornerNode=.FALSE.
      IF(cornerNode) THEN
        GeneratedMesh_ComponentNodeToUserNumber=position2(1)*cornerNodeFactor(1)+position2(2)*cornerNodeFactor(2)+ &
          & position2(3)*cornerNodeFactor(3)
        IF(generatedMesh%generatedType==GENERATED_MESH_ELLIPSOID_MESH_TYPE.AND.position2(2)/=0) THEN
          !Subtract off non-existent nodes at apex
          GeneratedMesh_ComponentNodeToUserNumber=GeneratedMesh_ComponentNodeToUserNumber-(numberOfElementsXi(1)-1)
        ENDIF
        GeneratedMesh_ComponentNodeToUserNumber=GeneratedMesh_ComponentNodeToUserNumber+1
      ELSE
        !subtract previous corner nodes from node offset
        numberOfPreviousCorners=0
        finishedCount=.FALSE.
        IF(numberOfDimensions>2) THEN
          IF(MOD(position(3),basis%numberOfNodesXiC(3)-1)/=0) THEN
            numberOfPreviousCorners=numberOfPreviousCorners+cornerNodeFactor(3)*(position2(3)+1)
            finishedCount=.TRUE.
          ELSE
            numberOfPreviousCorners=numberOfPreviousCorners+cornerNodeFactor(3)*position2(3)
          ENDIF
        ENDIF
        IF((numberOfDimensions>1).AND.(.NOT.finishedCount)) THEN
          IF(MOD(position(2),basis%numberOfNodesXiC(2)-1)/=0) THEN
            numberOfPreviousCorners=numberOfPreviousCorners+cornerNodeFactor(2)*(position2(2)+1)
            finishedCount=.TRUE.
          ELSE
            numberOfPreviousCorners=numberOfPreviousCorners+cornerNodeFactor(2)*position2(2)
          ENDIF
          IF(generatedMesh%generatedType==GENERATED_MESH_ELLIPSOID_MESH_TYPE) THEN
            numberOfPreviousCorners=numberOfPreviousCorners-(numberOfElementsXi(1)-1)
          ENDIF
        ENDIF
        IF(.NOT.finishedCount) THEN
          IF(MOD(position(1),basis%numberOfNodesXiC(1)-1)/=0) THEN
            numberOfPreviousCorners=numberOfPreviousCorners+cornerNodeFactor(1)*(position2(1)+1)
          ELSE
            numberOfPreviousCorners=numberOfPreviousCorners+cornerNodeFactor(1)*position2(1)
          ENDIF
        ENDIF
        nodeOffset=nodeOffset-numberOfPreviousCorners
        GeneratedMesh_ComponentNodeToUserNumber=nodeOffset+nodeComponentNumber
      ENDIF
    ENDIF

    EXITS("GeneratedMesh_ComponentNodeToUserNumber")
    RETURN
999 ERRORSEXITS("GeneratedMesh_ComponentNodeToUserNumber",err,error)
    RETURN
    
  END FUNCTION GeneratedMesh_ComponentNodeToUserNumber

  !
  !================================================================================================================================
  !
  !>Calculates the user node numbers for an array of nodes numbered using one basis for regular mesh type.
  !>1. For the current mesh component/basis, search previous basis to see if the current basis has occurred.
  !>2(1). If occurred, reuse user node number (i.e. same mesh topology)--> finish.
  !>2(2). If not occurred (i.e. different mesh topology), reuse corner nodes
  !> 3. Search previous basis to see if current interpolation scheme in xi1/2/3 direction has occurred in the same xi direction if previous basis.
  !>4(1). If occurred in xi1/2/3 direction, reuse node user numbers on corresponding edges/faces. e.g. linear-quadratic scheme v.s. biquadratic scheme, then node user numbers on edges alone xi2 direction can be reused.
  !>4(2). If never occurred (i.e. completely different basis. e.g. biquadratic v.s. bicubic), do nothing.
  !>5. Search previous basis to find the largest node user number, any new node user number will increment based on the current largest.
  !>6. Give node user numbers to nodes that have never appeared in previous basis.--> finish.

  SUBROUTINE GeneratedMesh_RegularComponentNodesToUserNumbers(generatedMesh,basisIndex,nodeComponentNumbers,nodeUserNumbers, &
    & err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh object
    INTEGER(INTG),INTENT(IN) :: basisIndex !<The number of the basis being used
    INTEGER(INTG),INTENT(IN) :: nodeComponentNumbers(:) !<The node numbers for this component basis
    INTEGER(INTG),INTENT(INOUT) :: nodeUserNumbers(:) !<On return, the corresponding user numbers
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: numberOfBases,numberOfDimensions,nodeOffsetLastBasis,lastElementNumber,nodeOffsetElement,offsetUnit
    INTEGER(INTG) :: nodeOffsetXi2Accum,nodeOffsetXi2,nodeOffset,nodeOffsetXi3Accum,elementNumber
    INTEGER(INTG) :: nodeIndexCurrent,nodeIndexFirst,nodeIndexPrevious
    INTEGER(INTG) :: nodeIdx,localNodeIdx1,localNodeIdx2,localNodeIdx3,xiIdx,basisIdx
    INTEGER(INTG) :: elementIndex(3),sameBasis(3),numberOfNodesXic(3),numberOfElementsXi(3),reminderTemp
    INTEGER(INTG) :: numberOfNodesTemp,nodeIndexTemp,nodeCount,indexCount,zeroCountXi1(16)
    INTEGER(INTG) :: zeroCountXi12(4),edgeNode(16),totalZeroNode,nodeOffsetElementXi12
    INTEGER(INTG) :: numberOfNodesLayer
    LOGICAL :: basisAppeared
    TYPE(BasisType), POINTER :: basis,basis2,basisFirstComponent,basisPrevious
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshElementsType), POINTER :: meshElements
    TYPE(MeshTopologyType), POINTER :: meshTopology
    

    ENTERS("GeneratedMesh_RegularComponentNodesToUserNumbers",err,error,*999)

    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is not associated.",err,error,*999)

    NULLIFY(regularMesh)
    CALL GeneratedMesh_RegularMeshGet(generatedMesh,regularMesh,err,error,*999)
    
    IF(SIZE(nodeUserNumbers)/=SIZE(nodeComponentNumbers)) &
      & CALL FlagError("Node component numbers and node user numbers arrays have different sizes.",err,error,*999)

    nodeUserNumbers=0
    numberOfBases=SIZE(regularMesh%bases)
    numberOfDimensions=regularMesh%meshDimension
    numberOfElementsXi=1
    DO xiIdx=1,numberOfDimensions
      numberOfElementsXi(xiIdx)=regularMesh%numberOfElementsXi(xiIdx)
    ENDDO !xiIdx
    !Number of nodes in each xi direction
    NULLIFY(basis)
    CALL GeneratedMeshRegular_BasisGet(regularMesh,basisIndex,basis,err,error,*999)
    numberOfNodesXic=1
    DO xiIdx=1,numberOfDimensions
      numberOfNodesXic(xiIdx)=basis%numberOfNodesXiC(xiIdx)
    ENDDO !xiIdx
    !Calculate current element indices and number
    reminderTemp=0;
    elementIndex=1;
    SELECT CASE(numberOfDimensions)
    CASE(1)
      !Calculate xi1 element index
      elementIndex(1)=(nodeComponentNumbers(1)-1)/(numberOfNodesXic(1)-1)+1
      !Calculate element number
      elementNumber=elementIndex(1)
    CASE(2)
      !Calculate xi2 element index
      numberOfNodesLayer=((numberOfNodesXic(1)-1)*numberOfElementsXi(1)+1)*(numberOfNodesXic(2)-1)
      elementIndex(2)=nodeComponentNumbers(1)/numberOfNodesLayer+1
      reminderTemp=MOD(nodeComponentNumbers(1),numberOfNodesLayer)
      !Calculate xi1 element index
      elementIndex(1)=(reminderTemp-1)/(numberOfNodesXic(1)-1)+1
      !Calculate element number
      elementNumber=(elementIndex(2)-1)*numberOfElementsXi(1)+elementIndex(1)
    CASE(3)
      !Calculate xi3 element index
      numberOfNodesLayer=((numberOfNodesXic(1)-1)*numberOfElementsXi(1)+1)*((numberOfNodesXic(2)-1)* &
        & numberOfElementsXi(2)+1)*(numberOfNodesXic(3)-1)
      elementIndex(3)=nodeComponentNumbers(1)/numberOfNodesLayer+1
      reminderTemp=MOD(nodeComponentNumbers(1),numberOfNodesLayer)
      !Calculate xi2 element index
      numberOfNodesLayer=((numberOfNodesXic(1)-1)*numberOfElementsXi(1)+1)*(numberOfNodesXic(2)-1)
      elementIndex(2)=reminderTemp/numberOfNodesLayer+1
      reminderTemp=MOD(reminderTemp,numberOfNodesLayer)
      !Calculate xi1 element index
      elementIndex(1)=(reminderTemp-1)/(numberOfNodesXic(1)-1)+1
      !Calculate element number
      elementNumber=(elementIndex(3)-1)*numberOfElementsXi(1)*numberOfElementsXi(2)+ &
        & (elementIndex(2)-1)*numberOfElementsXi(1)+elementIndex(1)
    END SELECT
    

    !If not the first basis, check if previous basis have same interpolation order in each xi direction
    !sameBasis(3) is initialised to have zeros in all entries. If an interpolation scheme has been
    !found to have appeared in previous basis, then record the basis number in the corresponding
    !xi direction. e.g. First basis: bi-quadratic, Second basis: quadratic-cubic, then sameBasis(3)
    !for the second basis will be [1,0,0]
    sameBasis=0
    DO xiIdx=1,numberOfDimensions
      DO basisIdx=1,basisIndex-1
        NULLIFY(basis2)
        CALL GeneratedMeshRegular_BasisGet(regularMesh,basisIdx,basis2,err,error,*999)
        IF(basis%numberOfNodesXiC(xiIdx)==basis2%numberOfNodesXiC(xiIdx)) sameBasis(xiIdx)=basisIdx
      ENDDO !basisIdx
    ENDDO !xiIdx
    !Check if the interpolation scheme has appeared in previous basis
    basisAppeared=.FALSE.
    IF(sameBasis(1)/=0) THEN
      SELECT CASE(numberOfDimensions)
      CASE(1)
        basisAppeared=.TRUE.
      CASE(2)
        IF(sameBasis(1)==sameBasis(2)) basisAppeared=.TRUE.
      CASE(3)
        IF(sameBasis(1)==sameBasis(2).AND.sameBasis(1)==sameBasis(3)) basisAppeared=.TRUE.
      END SELECT
    ENDIF
    NULLIFY(mesh)
    CALL GeneratedMesh_MeshGet(generatedMesh,mesh,err,error,*999)
    NULLIFY(meshTopology)
    CALL Mesh_MeshTopologyGet(mesh,sameBasis(1),meshTopology,err,error,*999)
    NULLIFY(meshElements)
    CALL MeshTopology_MeshElementsget(meshTopology,meshElements,err,error,*999)
    IF(basisIndex==1) THEN
      !If this is the first basis, don't do anything
      DO nodeIdx=1,SIZE(nodeComponentNumbers)
        nodeUserNumbers(nodeIdx)=nodeComponentNumbers(nodeIdx)
      ENDDO !nodeIdx
    ELSE IF(basisAppeared) THEN
      !If the basis has appeared before, reuse node user numbers
      DO nodeIdx=1,SIZE(nodeComponentNumbers)
        nodeUserNumbers(nodeIdx)=meshElements%elements(elementNumber)%userElementNodes(nodeIdx)
      ENDDO !nodeIdx
    ELSE
      !If the basis has never appeared exactly in previous basis

      !Find corner node user number from the first basis
      NULLIFY(basisFirstComponent)
      CALL GeneratedMeshRegular_BasisGet(regularMesh,1,basisFirstComponent,err,error,*999)
      NULLIFY(meshTopology)
      CALL Mesh_MeshTopologyGet(mesh,1,meshTopology,err,error,*999)
      NULLIFY(meshElements)
      CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
      lastElementNumber=meshElements%numberOfElements !The mesh has the same topology regardless of mesh components
       DO localNodeIdx3=1,2
        DO localNodeIdx2=1,2
          DO localNodeIdx1=1,2
            nodeIndexCurrent=localNodeIdx1
            nodeIndexFirst=localNodeIdx1
            IF(localNodeIdx1==2) THEN
              nodeIndexCurrent=numberOfNodesXic(1)
              nodeIndexFirst=basisFirstComponent%numberOfNodesXiC(1)
            ENDIF
            IF(numberOfDimensions>1 .AND. localNodeIdx2==2) THEN
              nodeIndexCurrent=nodeIndexCurrent+(numberOfNodesXic(2)-1)*numberOfNodesXic(1)
              nodeIndexFirst=nodeIndexFirst+(basisFirstComponent%numberOfNodesXiC(2)-1)* &
                & basisFirstComponent%numberOfNodesXiC(1)
            ENDIF
            IF(numberOfDimensions>2 .AND. localNodeIdx3==2) THEN
              nodeIndexCurrent=nodeIndexCurrent+numberOfNodesXic(1)* &
                & numberOfNodesXic(2)*(numberOfNodesXic(3)-1)
              nodeIndexFirst=nodeIndexFirst+basisFirstComponent%numberOfNodesXiC(1)* &
                & basisFirstComponent%numberOfNodesXiC(2)*(basisFirstComponent%numberOfNodesXiC(3)-1)
            ENDIF
            nodeUserNumbers(nodeIndexCurrent)=meshElements%elements(elementNumber)%globalElementNodes(nodeIndexFirst)
          ENDDO !localNodeIdx1
        ENDDO !localNodeIdx2
      ENDDO !localNodeIdx3
      
      !Find edge node user number from previous basis
      IF(sameBasis(1)/=0.AND.numberOfDimensions>1) THEN !Do not consider 1D since it's a complete new basis
        NULLIFY(basisPrevious)
        CALL GeneratedMeshRegular_BasisGet(regularMesh,sameBasis(1),basisPrevious,err,error,*999)
        NULLIFY(meshTopology)
        CALL Mesh_MeshTopologyGet(mesh,sameBasis(1),meshTopology,err,error,*999)
        NULLIFY(meshElements)
        CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
        DO localNodeIdx3=1,2
          DO localNodeIdx2=1,2
            DO localNodeIdx1=2,numberOfNodesXic(1)-1
              nodeIndexCurrent=localNodeIdx1
              nodeIndexPrevious=localNodeIdx1
              IF(localNodeIdx2==2) THEN
                nodeIndexCurrent=nodeIndexCurrent+(numberOfNodesXic(2)-1)*numberOfNodesXic(1)
                nodeIndexPrevious=nodeIndexPrevious+(basisPrevious%numberOfNodesXiC(2)-1)*basisPrevious%numberOfNodesXiC(1)
              ENDIF
              IF(numberOfDimensions>2 .AND. localNodeIdx3==2) THEN
                nodeIndexCurrent=nodeIndexCurrent+numberOfNodesXic(1)*numberOfNodesXic(2)* &
                  & (numberOfNodesXic(3)-1)
                nodeIndexPrevious=nodeIndexPrevious+basisPrevious%numberOfNodesXiC(1)*basisPrevious% &
                  & numberOfNodesXiC(2)*(basisPrevious%numberOfNodesXiC(3)-1)
              ENDIF
              nodeUserNumbers(nodeIndexCurrent)=meshElements%elements(elementNumber)%globalElementNodes(nodeIndexPrevious)
            ENDDO !localNodeIdx1
          ENDDO !localNodeIdx2
        ENDDO !localNodeIdx3
      ENDIF
      IF(sameBasis(2)/=0) THEN
        NULLIFY(basisPrevious)
        CALL GeneratedMeshRegular_BasisGet(regularMesh,sameBasis(2),basisPrevious,err,error,*999)
        NULLIFY(meshTopology)
        CALL Mesh_MeshTopologyGet(mesh,sameBasis(2),meshTopology,err,error,*999)
        NULLIFY(meshElements)
        CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
        DO localNodeIdx3=1,2
          DO localNodeIdx2=2,numberOfNodesXic(2)-1
            DO localNodeIdx1=1,2
              IF(localNodeIdx1==1) THEN
                nodeIndexCurrent=localNodeIdx1+(localNodeIdx2-1)*numberOfNodesXic(1)
                nodeIndexPrevious=localNodeIdx1+(localNodeIdx2-1)*basisPrevious%numberOfNodesXiC(1)
              ELSE
                nodeIndexCurrent=localNodeIdx2*numberOfNodesXic(1)
                nodeIndexPrevious=localNodeIdx2*basisPrevious%numberOfNodesXiC(1)
              ENDIF
              IF(numberOfDimensions>2.AND.localNodeIdx3==2) THEN
                nodeIndexCurrent=nodeIndexCurrent+numberOfNodesXic(1)*numberOfNodesXic(2)* &
                  & (numberOfNodesXic(3)-1)
                nodeIndexPrevious=nodeIndexPrevious+basisPrevious%numberOfNodesXiC(1)*basisPrevious% &
                  & numberOfNodesXiC(2)*(basisPrevious%numberOfNodesXiC(3)-1)
              ENDIF
              nodeUserNumbers(nodeIndexCurrent)=meshElements%elements(elementNumber)%globalElementNodes(nodeIndexPrevious)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF(sameBasis(3)/=0) THEN !Must be 3D
        NULLIFY(basisPrevious)
        CALL GeneratedMeshRegular_BasisGet(regularMesh,sameBasis(3),basisPrevious,err,error,*999)
        NULLIFY(meshTopology)
        CALL Mesh_MeshTopologyGet(mesh,sameBasis(3),meshTopology,err,error,*999)
        NULLIFY(meshElements)
        CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
        nodeIndexCurrent=0
        nodeIndexPrevious=0
        DO localNodeIdx3=2,numberOfNodesXic(3)-1
          DO localNodeIdx2=1,2
            IF(localNodeIdx2==2) THEN
              nodeIndexCurrent=(numberOfNodesXic(2)-1)*numberOfNodesXic(1)+numberOfNodesXic(1)* &
                & numberOfNodesXic(2)*(numberOfNodesXic(3)-1)
              nodeIndexPrevious=(basisPrevious%numberOfNodesXiC(1)-1)*basisPrevious%numberOfNodesXiC(1)+ &
                & basisPrevious%numberOfNodesXiC(1)*basisPrevious%numberOfNodesXiC(2)* &
                & (basisPrevious%numberOfNodesXiC(3)-1)
            ENDIF
            DO localNodeIdx1=1,2
              IF(localNodeIdx1==1) THEN
                nodeIndexCurrent=1+nodeIndexCurrent
                nodeIndexPrevious=1+nodeIndexPrevious
              ELSE
                nodeIndexCurrent=numberOfNodesXic(1)+nodeIndexCurrent
                nodeIndexPrevious=basisPrevious%numberOfNodesXiC(1)+nodeIndexPrevious
              ENDIF
              nodeUserNumbers(nodeIndexCurrent)=meshElements%elements(elementNumber)%globalElementNodes(nodeIndexPrevious)
            ENDDO !localNodeIdx1
          ENDDO !localNodeIdx2
        ENDDO !localNodeIdx3
      ENDIF
      !The following code would only be executed if 3D (automatically satisfied, don't need to check,
      !since there must be at least 1 direction that has different interpolation scheme, if two direction
      ! has the same interpolation that has appeared before, then interpolation for the last direction
      ! must be different) and has same basis in 2 xi direction
      !i.e. find user node numbers for face nodes
      IF(sameBasis(1)==sameBasis(2).AND.sameBasis(1)/=0) THEN
        NULLIFY(basisPrevious)
        CALL GeneratedMeshRegular_BasisGet(regularMesh,sameBasis(1),basisPrevious,err,error,*999)
        NULLIFY(meshTopology)
        CALL Mesh_MeshTopologyGet(mesh,sameBasis(1),meshTopology,err,error,*999)
        NULLIFY(meshElements)
        CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
        DO localNodeIdx3=1,2
          DO localNodeIdx2=2,numberOfNodesXic(2)-1
            DO localNodeIdx1=2,numberOfNodesXic(1)-1
              nodeIndexCurrent=localNodeIdx1+(localNodeIdx2-1)*numberOfNodesXic(1)
              nodeIndexPrevious=localNodeIdx1+(localNodeIdx2-1)*basisPrevious%numberOfNodesXiC(1)
              IF(localNodeIdx3==2) THEN
                nodeIndexCurrent=nodeIndexCurrent+numberOfNodesXic(1)*numberOfNodesXic(2)* &
                  & (numberOfNodesXic(3)-1)
                nodeIndexPrevious=nodeIndexPrevious+basisPrevious%numberOfNodesXiC(1)*basisPrevious% &
                  & numberOfNodesXiC(2)*(basisPrevious%numberOfNodesXiC(3)-1)
              ENDIF
              nodeUserNumbers(nodeIndexCurrent)=meshElements%elements(elementNumber)%globalElementNodes(nodeIndexPrevious)
            ENDDO !localNodeIdx1
          ENDDO !localNodeIdx2
        ENDDO !localNodeIdx3
      ELSE IF(sameBasis(1)==sameBasis(3).AND.sameBasis(1)/=0) THEN
        NULLIFY(basisPrevious)
        CALL GeneratedMeshRegular_BasisGet(regularMesh,sameBasis(1),basisPrevious,err,error,*999)
        NULLIFY(meshTopology)
        CALL Mesh_MeshTopologyGet(mesh,sameBasis(1),meshTopology,err,error,*999)
        NULLIFY(meshElements)
        CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
        nodeIndexCurrent=0
        nodeIndexPrevious=0
        DO localNodeIdx3=2,numberOfNodesXic(3)-1
          DO localNodeIdx2=1,2
            IF(localNodeIdx2==2) THEN
              nodeIndexCurrent=(numberOfNodesXic(2)-1)*numberOfNodesXic(1)+numberOfNodesXic(1)* &
                & numberOfNodesXic(2)*(localNodeIdx3-1)
              nodeIndexPrevious=(basisPrevious%numberOfNodesXiC(2)-1)*basisPrevious%numberOfNodesXiC(1)+ &
                & basisPrevious%numberOfNodesXiC(1)*basisPrevious%numberOfNodesXiC(2)*(localNodeIdx3-1)
            ENDIF
            DO localNodeIdx1=2,numberOfNodesXic(1)-1
              nodeIndexCurrent=localNodeIdx1+nodeIndexCurrent
              nodeIndexPrevious=localNodeIdx1+nodeIndexPrevious
              nodeUserNumbers(nodeIndexCurrent)=meshElements%elements(elementNumber)%globalElementNodes(nodeIndexPrevious)
            ENDDO !localNodeIdx1
          ENDDO !localNodeIdx2
        ENDDO !localNodeIdx3
      ELSE IF(sameBasis(2)==sameBasis(3).AND.sameBasis(2)/=0) THEN
        NULLIFY(basisPrevious)
        CALL GeneratedMeshRegular_BasisGet(regularMesh,sameBasis(2),basisPrevious,err,error,*999)
        NULLIFY(meshTopology)
        CALL Mesh_MeshTopologyGet(mesh,sameBasis(2),meshTopology,err,error,*999)
        NULLIFY(meshElements)
        CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
        DO localNodeIdx3=2,numberOfNodesXic(3)-1
          DO localNodeIdx2=2,numberOfNodesXic(2)-1
            DO localNodeIdx1=1,2
              IF(localNodeIdx1==1) THEN
                nodeIndexCurrent=1+(localNodeIdx2-1)*numberOfNodesXic(1)+numberOfNodesXic(1)* &
                  & numberOfNodesXic(2)*(localNodeIdx3-1)
                nodeIndexPrevious=1+(localNodeIdx2-1)*basisPrevious%numberOfNodesXiC(1)+basisPrevious%numberOfNodesXiC(1)* &
                  & basisPrevious%numberOfNodesXiC(2)*(localNodeIdx3-1)
              ELSE
                nodeIndexCurrent=localNodeIdx2*numberOfNodesXic(1)+numberOfNodesXic(1)* &
                  & numberOfNodesXic(2)*(localNodeIdx3-1)
                nodeIndexPrevious=localNodeIdx2*basisPrevious%numberOfNodesXiC(1)+basisPrevious%numberOfNodesXiC(1)* &
                  & basisPrevious%numberOfNodesXiC(2)*(localNodeIdx3-1)
              ENDIF
              nodeUserNumbers(nodeIndexCurrent)=meshElements%elements(elementNumber)%globalElementNodes(nodeIndexPrevious)
            ENDDO !localNodeIdx1
          ENDDO !localNodeIdx2
        ENDDO !localNodeIdx3
      ENDIF
      
      !Find the largest node user number in the previous basis
      nodeOffsetLastBasis=0
      DO basisIdx=1,basisIndex-1
        NULLIFY(meshTopology)
        CALL Mesh_MeshTopologyGet(mesh,basisIdx,meshTopology,err,error,*999)
        NULLIFY(meshElements)
        CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
        numberOfNodesTemp=SIZE(meshElements%elements(lastElementNumber)%globalElementNodes,1)
        DO nodeIndexTemp=1,numberOfNodesTemp
          IF(meshElements%elements(lastElementNumber)%globalElementNodes(nodeIndexTemp)>nodeOffsetLastBasis) &
            & nodeOffsetLastBasis=meshElements%elements(lastElementNumber)%globalElementNodes(nodeIndexTemp)
        ENDDO !nodeIndexTemp
      ENDDO !basisIdx
      
      !Calculate number of zeros nodes in different dimensions
      indexCount=1
      zeroCountXi1=0
      zeroCountXi12=0
      totalZeroNode=0
      edgeNode=0
      DO localNodeIdx3=1,numberOfNodesXic(3)
        DO localNodeIdx2=1,numberOfNodesXic(2)
          nodeCount=0
          DO localNodeIdx1=1,numberOfNodesXic(1)
            nodeIdx=(localNodeIdx3-1)*numberOfNodesXic(1)*numberOfNodesXic(2)+(localNodeIdx2-1)*numberOfNodesXic(1)+localNodeIdx1
            IF(nodeUserNumbers(nodeIdx)==0) THEN
              nodeCount=nodeCount+1
              totalZeroNode=totalZeroNode+1 !Total number of zeros in an element
            ENDIF
          ENDDO !localNodeIdx1
          zeroCountXi1(indexCount)=nodeCount !Total number of zero summed up across xi1 direction.
          IF(nodeCount==numberOfNodesXic(1)) edgeNode(indexCount)=1 !Shared edge node (with zero value) in xi1 direction (1 number for each node in xi2 direction)
          zeroCountXi12(localNodeIdx3)=zeroCountXi12(localNodeIdx3)+zeroCountXi1(indexCount) !Total number of zero summed on xi1-xi2 faces
          indexCount=indexCount+1
        ENDDO !localNodeIdx2
      ENDDO !localNodeIdx3
      
      !Calculate how many zero nodes has occurred in previous elements
      nodeOffsetElement=0
      IF(numberOfDimensions==2.AND.elementIndex(2)/=1) THEN !Zero nodes occurred in the previous rows of elements
        offsetUnit=totalZeroNode-zeroCountXi1(1)-SUM(edgeNode(1:numberOfNodesXic(2)))+edgeNode(indexCount)
        !This is number of zero nodes in the elements before the current row of elements
        nodeOffsetElement=(elementIndex(2)-1)*numberOfElementsXi(1)*offsetUnit+(elementIndex(2)-1)* &
          & SUM(edgeNode(2:numberOfNodesXic(2)-1))
      ELSE IF(numberOfDimensions==3.AND.elementIndex(3)/=1) THEN !Zero nodes occurred in the previous layer of elements
        nodeOffsetXi3Accum=0
        DO localNodeIdx3=1,numberOfNodesXic(3)-1
          offsetUnit=zeroCountXi12(localNodeIdx3)-zeroCountXi1((localNodeIdx3-1)*numberOfNodesXic(2)+1)- &
            & SUM(edgeNode((localNodeIdx3-1)*numberOfNodesXic(2)+1:localNodeIdx3*numberOfNodesXic(2)))+ &
            & edgeNode((localNodeIdx3-1)*numberOfNodesXic(2)+1)
          nodeOffsetXi3Accum=nodeOffsetXi3Accum+offsetUnit*numberOfElementsXi(1)*numberOfElementsXi(2)+ &
            & (numberOfElementsXi(1)-1)*(zeroCountXi1((localNodeIdx3-1)*numberOfNodesXic(2)+1)- &
            & edgeNode((localNodeIdx3-1)*numberOfNodesXic(2)+1))+zeroCountXi1((localNodeIdx3-1)*numberOfNodesXic(2)+1)+ &
            & SUM(edgeNode((localNodeIdx3-1)*numberOfNodesXic(2)+2:localNodeIdx3*numberOfNodesXic(2)))* &
            & numberOfElementsXi(2)
        ENDDO !localNodeIdx3
        nodeOffsetElement=(elementIndex(3)-1)*nodeOffsetXi3Accum
      ENDIF
      
      !Compute other nodes which haven't appeared in previous basis
      indexCount=1
      nodeOffsetElementXi12=0
      nodeOffsetXi2=0 !Number of zero nodes in the current row
      nodeOffsetXi3Accum=0 !Number of zero nodes in the layers in xi3 direction (localNodeIdx3)
      DO localNodeIdx3=1,numberOfNodesXic(3)
        nodeOffsetXi2Accum=0 !Number of zero nodes in the previous rows
        offsetUnit=zeroCountXi12(localNodeIdx3)-zeroCountXi1((localNodeIdx3-1)*numberOfNodesXic(2)+1)- &
          & SUM(edgeNode((localNodeIdx3-1)*numberOfNodesXic(2)+1:localNodeIdx3*numberOfNodesXic(2)))+ &
          & edgeNode((localNodeIdx3-1)*numberOfNodesXic(2)+1)
        IF(elementIndex(2)/=1 .AND. numberOfDimensions==3) THEN
          nodeOffsetElementXi12=offsetUnit*(elementIndex(2)-1)*numberOfElementsXi(1)+ &
            & (elementIndex(2)-1)*SUM(edgeNode((localNodeIdx3-1)*numberOfNodesXic(2)+2:localNodeIdx3*numberOfNodesXic(2)))
        ENDIF
        DO localNodeIdx2=1,numberOfNodesXic(2)
          nodeOffsetXi2=(zeroCountXi1(indexCount)-edgeNode(indexCount))*(elementIndex(1)-1)
          nodeOffset=nodeOffsetLastBasis+nodeOffsetElement+nodeOffsetXi3Accum+ &
            & nodeOffsetElementXi12+nodeOffsetXi2Accum+nodeOffsetXi2
          DO localNodeIdx1=1,numberOfNodesXic(1)
            !Local node index in the current element
            nodeIdx=(localNodeIdx3-1)*numberOfNodesXic(1)*numberOfNodesXic(2)+(localNodeIdx2-1)* &
              & numberOfNodesXic(1)+localNodeIdx1
            IF(nodeUserNumbers(nodeIdx)==0) THEN
              !This is for 2D case
              nodeOffset=nodeOffset+1
              nodeUserNumbers(nodeIdx)=nodeOffset
            ENDIF
          ENDDO !localNodeIdx1
          nodeOffsetXi2Accum=nodeOffsetXi2Accum+(zeroCountXi1(indexCount)-edgeNode(indexCount))* &
            & numberOfElementsXi(1)+edgeNode(indexCount)
          indexCount=indexCount+1
        ENDDO !localNodeIdx2
        IF(numberOfDimensions==3) THEN
          nodeOffsetXi3Accum=nodeOffsetXi3Accum+offsetUnit*numberOfElementsXi(1)*numberOfElementsXi(2)+ &
            & (numberOfElementsXi(1)-1)*(zeroCountXi1((localNodeIdx3-1)*numberOfNodesXic(2)+1)- &
            & edgeNode((localNodeIdx3-1)*numberOfNodesXic(2)+1))+zeroCountXi1((localNodeIdx3-1)*numberOfNodesXic(2)+1)+ &
            & SUM(edgeNode((localNodeIdx3-1)*numberOfNodesXic(2)+2:localNodeIdx3*numberOfNodesXic(2)))* &
            & numberOfElementsXi(2)
        ENDIF
      ENDDO !localNodeIdx3
    ENDIF
    
    EXITS("GeneratedMesh_RegularComponentNodesToUserNumbers")
    RETURN
999 ERRORS("GeneratedMesh_RegularComponentNodesToUserNumbers",err,error)
    EXITS("GeneratedMesh_RegularComponentNodesToUserNumbers")
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_RegularComponentNodesToUserNumbers

  !
  !================================================================================================================================
  !

  !>Retrieve the user node number for a component number in a regular generated mesh
  !>This routine only works for Lagrange/Hermite elements
  SUBROUTINE GeneratedMesh_RegularComponentNodeToUserNumber(generatedMesh,basisIndex,nodeComponentNumber,nodeUserNumber, &
    & err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh  !<A pointer to the generated mesh object
    INTEGER(INTG),INTENT(IN) :: basisIndex  !<The number of the basis being used
    INTEGER(INTG),INTENT(IN) :: nodeComponentNumber  !<The node numbers for this component basis
    INTEGER(INTG),INTENT(OUT) :: nodeUserNumber  !<On return, the corresponding user numbers
    INTEGER(INTG) :: err  !<The error code
    TYPE(VARYING_STRING) :: error  !<The error string
    !Local variables
    INTEGER(INTG) :: numberOfBases,numberOfDimensions,elementNumber,localNodeNumber,numberOfNodesLayer,xiIdx
    INTEGER(INTG) :: elementIndex(3),nodeIdx(3),numberOfNodesXic(3),numberOfElementsXi(3),reminderTemp
    TYPE(BasisType), POINTER :: basis
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshElementsType), POINTER :: meshElements
    TYPE(MeshTopologyType), POINTER :: meshTopology

    ENTERS("GeneratedMesh_RegularComponentNodeToUserNumber",err,error,*999)

    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is not associated",err,error,*999)

    NULLIFY(regularMesh)
    CALL GeneratedMesh_RegularMeshGet(generatedMesh,regularMesh,err,error,*999)
    numberOfBases=SIZE(regularMesh%bases)
    numberOfDimensions=regularMesh%meshDimension
    numberOfElementsXi=1
    DO xiIdx=1,numberOfDimensions
      numberOfElementsXi(xiIdx)=regularMesh%numberOfElementsXi(xiIdx)
    ENDDO !xiIdx
    !Number of nodes in each xi direction
    NULLIFY(basis)
    CALL GeneratedMeshRegular_BasisGet(regularMesh,basisIndex,basis,err,error,*999)
    numberOfNodesXic=1
    DO xiIdx=1,numberOfDimensions
      numberOfNodesXic(xiIdx)=basis%numberOfNodesXiC(xiIdx)
    ENDDO !xiIdx
 
    !Calculate current element/node indices/number
    reminderTemp=0;
    elementIndex=1;
    nodeIdx=1;
    SELECT CASE(numberOfDimensions)
    CASE(1)
      !Calculate xi1 element index
      elementIndex(1)=(nodeComponentNumber-1)/(numberOfNodesXic(1)-1)+1
      nodeIdx(1)=MOD(nodeComponentNumber-1,numberOfNodesXic(1)-1)+1
      !If it's the last node in the line
      IF (elementIndex(1)>numberOfElementsXi(1)) THEN
        elementIndex(1)=elementIndex(1)-1
        nodeIdx(1)=numberOfNodesXic(1)
      ENDIF
      !Calculate element number
      elementNumber=elementIndex(1)
      localNodeNumber=nodeIdx(1)
    CASE(2)
      !Calculate xi2 element index
      numberOfNodesLayer=((numberOfNodesXic(1)-1)*numberOfElementsXi(1)+1)*(numberOfNodesXic(2)-1)
      elementIndex(2)=(nodeComponentNumber-1)/numberOfNodesLayer+1
      reminderTemp=MOD(nodeComponentNumber-1,numberOfNodesLayer)
      numberOfNodesLayer=((numberOfNodesXic(1)-1)*numberOfElementsXi(1)+1)
      nodeIdx(2)=reminderTemp/numberOfNodesLayer+1
      !If it's the last line of nodes in the line
      IF (elementIndex(2)>numberOfElementsXi(2)) THEN
        elementIndex(2)=elementIndex(2)-1
        nodeIdx(2)=numberOfNodesXic(2)
      ENDIF
      !Calculate xi1 element index
      reminderTemp=MOD(reminderTemp,numberOfNodesLayer)
      elementIndex(1)=reminderTemp/(numberOfNodesXic(1)-1)+1
      nodeIdx(1)=MOD(reminderTemp,numberOfNodesXic(1)-1)+1
      !If it's the last node in the line
      IF (elementIndex(1)>numberOfElementsXi(1)) THEN
        elementIndex(1)=elementIndex(1)-1
        nodeIdx(1)=numberOfNodesXic(1)
      ENDIF
      !Calculate element number
      elementNumber=(elementIndex(2)-1)*numberOfElementsXi(1)+elementIndex(1)
      localNodeNumber=(nodeIdx(2)-1)*numberOfNodesXic(1)+nodeIdx(1)
    CASE(3)
      !Calculate xi3 element index
      numberOfNodesLayer=((numberOfNodesXic(1)-1)*numberOfElementsXi(1)+1)*((numberOfNodesXic(2)-1)* &
        & numberOfElementsXi(2)+1)*(numberOfNodesXic(3)-1) !Multiple planes of nodes
      elementIndex(3)=(nodeComponentNumber-1)/numberOfNodesLayer+1
      reminderTemp=MOD(nodeComponentNumber-1,numberOfNodesLayer) !Multiple planes of nodes
      numberOfNodesLayer=((numberOfNodesXic(1)-1)*numberOfElementsXi(1)+1)*((numberOfNodesXic(2)-1)* &
        & numberOfElementsXi(2)+1) !One plane of nodes
      nodeIdx(3)=reminderTemp/numberOfNodesLayer+1
      IF (elementIndex(3)>numberOfElementsXi(3)) THEN
        elementIndex(3)=elementIndex(3)-1
        nodeIdx(3)=numberOfNodesXic(3)
      ENDIF
      reminderTemp=MOD(reminderTemp,numberOfNodesLayer) !One plane of nodes
      !Calculate xi2 element index
      numberOfNodesLayer=((numberOfNodesXic(1)-1)*numberOfElementsXi(1)+1)*(numberOfNodesXic(2)-1) !Multiple lines of nodes
      elementIndex(2)=reminderTemp/numberOfNodesLayer+1
      reminderTemp=MOD(reminderTemp,numberOfNodesLayer) !Multiple lines of nodes
      numberOfNodesLayer=(numberOfNodesXic(1)-1)*numberOfElementsXi(1)+1 !One line of nodes
      nodeIdx(2)=reminderTemp/numberOfNodesLayer+1
      reminderTemp=MOD(reminderTemp,numberOfNodesLayer) !One line of nodes
      IF (elementIndex(2)>numberOfElementsXi(2)) THEN
        elementIndex(2)=elementIndex(2)-1
        nodeIdx(2)=numberOfNodesXic(2)
      ENDIF
      !Calculate xi1 element index
      elementIndex(1)=reminderTemp/(numberOfNodesXic(1)-1)+1
      nodeIdx(1)=MOD(reminderTemp,numberOfNodesXic(1)-1)+1
      IF (elementIndex(1)>numberOfElementsXi(1)) THEN
        elementIndex(1)=elementIndex(1)-1
        nodeIdx(1)=numberOfNodesXic(1)
      ENDIF
      !Calculate element number
      elementNumber=(elementIndex(3)-1)*numberOfElementsXi(1)*numberOfElementsXi(2)+ &
        & (elementIndex(2)-1)*numberOfElementsXi(1)+elementIndex(1)
      localNodeNumber=(nodeIdx(3)-1)*numberOfNodesXic(1)*numberOfNodesXic(2)+(nodeIdx(2)-1)*numberOfNodesXic(1)+nodeIdx(1)
    END SELECT
    !Retrieve node user number
    NULLIFY(mesh)
    CALL GeneratedMesh_MeshGet(generatedMesh,mesh,err,error,*999)
    NULLIFY(meshTopology)
    CALL Mesh_MeshTopologyGet(mesh,basisIndex,meshTopology,err,error,*999)
    NULLIFY(meshElements)
    CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
    nodeUserNumber=meshElements%elements(elementNumber)%userElementNodes(localNodeNumber)

    EXITS("GeneratedMesh_RegularComponentNodeToUserNumber")
    RETURN
999 ERRORS("GeneratedMesh_RegularComponentNodeToUserNumber",err,error)
    EXITS("GeneratedMesh_RegularComponentNodeToUserNumber")
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_RegularComponentNodeToUserNumber

  !
  !================================================================================================================================
  !

  !>Calculates the user node number for a node numbered using one basis.
  !>This is currently only used for cylinder meshes, other mesh types don't require this.
  FUNCTION GeneratedMesh_UserNumberToComponentNode(generatedMesh,basisIndex,nodeUserNumber,err,error)
    TYPE(GeneratedMeshType), POINTER :: generatedMesh        !<A pointer to the generated mesh object
    INTEGER(INTG),INTENT(IN) :: basisIndex                     !<The number of the basis being used
    INTEGER(INTG),INTENT(IN) :: nodeUserNumber                !<The corresponding user node number
    INTEGER(INTG) :: err          !<The error code
    TYPE(VARYING_STRING) :: error !<The error string
    !function variable
    INTEGER(INTG) :: GeneratedMesh_UserNumberToComponentNode !<On return, the node number for this component basis
    !local variables
    INTEGER(INTG) :: numberOfBases,numberOfDimensions,basisIdx,xiIdx,remainder,tempTerm,numberOfCornerNodes,nodeOffset
    INTEGER(INTG) :: position(3),cornerNodeFactor(3),basisElementFactor(3),numberOfPreviousCorners,basisNumberOfNodes
    INTEGER(INTG), POINTER :: numberOfElementsXi(:)
    LOGICAL :: finishedCount,offEdge
    TYPE(BasisType), POINTER :: basis,basis2
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh
     TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_UserNumberToComponentNode",err,error,*999)

    NULLIFY(basis)

    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is not associated",err,error,*999)

    numberOfCornerNodes=1
    remainder=nodeUserNumber-1 !use zero based numbering
    position=0
    
    !Only cylinder mesh type uses this now, although it was previously used by regular
    !meshes so some things relate to that.
    SELECT CASE(generatedMesh%generatedType)
    CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_POLAR_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
      NULLIFY(cylinderMesh)
      CALL GeneratedMesh_CylinderMeshGet(generatedMesh,cylinderMesh,err,error,*999)
      numberOfBases=SIZE(cylinderMesh%bases)
      numberOfDimensions=cylinderMesh%meshDimension
      numberOfElementsXi=cylinderMesh%numberOfElementsXi
    CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The generated mesh generated type of "// &
        & TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    IF(basisIndex<1.OR.basisIndex>numberOfBases) THEN
      localError="Mesh component must be less than or equal to "//(NumberToVString(numberOfBases,"*",err,error))// &
        & " but it is "//(NumberToVString(basisIndex,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    IF(numberOfBases==1) THEN
      !If is the only basis, don't do anything
      GeneratedMesh_UserNumberToComponentNode=nodeUserNumber
    ELSE
      tempTerm=1
      numberOfCornerNodes=1
      DO xiIdx=1,numberOfDimensions
        numberOfCornerNodes=numberOfCornerNodes*(numberOfElementsXi(xiIdx)+1)
        cornerNodeFactor(xiIdx)=1
        IF(xiIdx>1) THEN
          tempTerm=tempTerm*(numberOfElementsXi(xiIdx-1)+1)
          cornerNodeFactor(xiIdx)=cornerNodeFactor(xiIdx)*tempTerm
        ENDIF
      ENDDO !xiIdx
      !Adjust for other mesh types
      IF(generatedMesh%generatedType==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
        cornerNodeFactor(3)=cornerNodeFactor(3)-numberOfElementsXi(1)-1
        numberOfCornerNodes=numberOfCornerNodes-(numberOfElementsXi(1)+1)*(numberOfElementsXi(3)+1)
      ENDIF
      nodeOffset=numberOfCornerNodes
      DO basisIdx=1,basisIndex-1
        NULLIFY(basis2)
        CALL GeneratedMeshCylinder_BasisGet(cylinderMesh,basisIdx,basis2,err,error,*999)
        basisNumberOfNodes=1
        DO xiIdx=1,numberOfDimensions
          basisNumberOfNodes=basisNumberOfNodes*(numberOfElementsXi(xiIdx)*(basis2%numberOfNodesXiC(xiIdx)-1)+1)
        ENDDO !xiIdx
        !Adjust for other mesh types
        IF(generatedMesh%generatedType==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
          basisNumberOfNodes=basisNumberOfNodes-(numberOfElementsXi(1)+1)*(basis2%numberOfNodesXiC(1)-1)* &
            & (numberOfElementsXi(3)+1)*(basis2%numberOfNodesXiC(3)-1)
        ENDIF
        nodeOffset=nodeOffset+basisNumberOfNodes-numberOfCornerNodes
      ENDDO !basisIdx
      NULLIFY(basis)
      CALL GeneratedMeshCylinder_BasisGet(cylinderMesh,basisIndex,basis,err,error,*999)
      tempTerm=1
      DO xiIdx=1,numberOfDimensions
        basisElementFactor(xiIdx)=basis%numberOfNodesXiC(xiIdx)-1
        IF(xiIdx>1) THEN
          tempTerm=tempTerm*((basis%numberOfNodesXiC(xiIdx-1)-1)*numberOfElementsXi(xiIdx-1)+1)
          basisElementFactor(xiIdx)=basisElementFactor(xiIdx)*tempTerm
        ENDIF
      ENDDO !xiIdx
      !Adjust for other mesh types
      IF(generatedMesh%generatedType==GENERATED_MESH_CYLINDER_MESH_TYPE) THEN
        !subtract nodes along line where y wraps around
        basisElementFactor(3)=basisElementFactor(3)-(numberOfElementsXi(1)* &
          & (basis%numberOfNodesXiC(1)-1)+1)*(basis%numberOfNodesXiC(3)-1)
      ENDIF
      IF(nodeUserNumber<=numberOfCornerNodes) THEN
        !we have a node on a corner
        IF(numberOfDimensions>2) THEN
          position(3)=remainder/cornerNodeFactor(3)
          remainder=MOD(remainder,cornerNodeFactor(3))
        ENDIF
        IF(numberOfDimensions>1) THEN
          position(2)=remainder/cornerNodeFactor(2)
          remainder=MOD(remainder,cornerNodeFactor(2))
        ENDIF
        position(1)=remainder/cornerNodeFactor(1)
        GeneratedMesh_UserNumberToComponentNode=position(1)*basisElementFactor(1)+position(2)*basisElementFactor(2)+ &
          & position(3)*basisElementFactor(3)
        GeneratedMesh_UserNumberToComponentNode=GeneratedMesh_UserNumberToComponentNode+1
      ELSE IF(nodeUserNumber>nodeOffset) THEN
        remainder=remainder-nodeOffset
        DO xiIdx=1,numberOfDimensions
          basisElementFactor(xiIdx)=basisElementFactor(xiIdx)-cornerNodeFactor(xiIdx)
        ENDDO
        numberOfPreviousCorners=0
        finishedCount=.FALSE.
        offEdge=.FALSE.
        IF(numberOfDimensions>2) THEN
          IF(generatedMesh%generatedType==GENERATED_MESH_CYLINDER_MESH_TYPE.AND. &
            & (MOD(remainder,basisElementFactor(3)) > basisElementFactor(2)*numberOfElementsXi(2)-1)) THEN
            offEdge=.TRUE.
          ELSE IF(generatedMesh%generatedType==GENERATED_MESH_REGULAR_MESH_TYPE.AND. &
            & MOD(remainder,basisElementFactor(3)) > (basisElementFactor(2)*numberOfElementsXi(2)+ &
            & basisElementFactor(1)*numberOfElementsXi(1)-1)) THEN
            offEdge=.TRUE.
          ENDIF
          IF(offEdge) THEN
            numberOfPreviousCorners=numberOfPreviousCorners+cornerNodeFactor(3)*(1+remainder/basisElementFactor(3))
            remainder=MOD(remainder,basisElementFactor(3))
            finishedCount=.TRUE.
          ELSE
            numberOfPreviousCorners=numberOfPreviousCorners+cornerNodeFactor(3)*(remainder/basisElementFactor(3))
            remainder=MOD(remainder,basisElementFactor(3))
          ENDIF
        ENDIF
        IF((numberOfDimensions>1) .AND. (finishedCount.NEQV..TRUE.)) THEN
          IF(MOD(remainder,basisElementFactor(2)) > &
            & basisElementFactor(1)*numberOfElementsXi(1)-1) THEN
            numberOfPreviousCorners=numberOfPreviousCorners+cornerNodeFactor(2)*(1+remainder/basisElementFactor(2))
            remainder=MOD(remainder,basisElementFactor(2))
            finishedCount=.TRUE.
          ELSE
            numberOfPreviousCorners=numberOfPreviousCorners+cornerNodeFactor(2)*(remainder/basisElementFactor(2))
            remainder=MOD(remainder,basisElementFactor(2))
          ENDIF
        ENDIF
        IF(finishedCount.NEQV..TRUE.) THEN
          numberOfPreviousCorners=numberOfPreviousCorners+cornerNodeFactor(1)*(remainder/basisElementFactor(1))+1
        ENDIF
        nodeOffset=nodeOffset-numberOfPreviousCorners
        GeneratedMesh_UserNumberToComponentNode=nodeUserNumber-nodeOffset
      ELSE
        CALL FlagError("Invalid node number specified.",err,error,*999)
      ENDIF
    ENDIF
 
    EXITS("GeneratedMesh_UserNumberToComponentNode")
    RETURN
999 ERRORSEXITS("GeneratedMesh_UserNumberToComponentNode",err,error)
    RETURN
    
  END FUNCTION GeneratedMesh_UserNumberToComponentNode

  !
  !================================================================================================================================
  !

END MODULE GeneratedMeshRoutines

