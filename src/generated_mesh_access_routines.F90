!> \file
!> \author Chris Bradley
!> \brief This module contains all generated mesh access method routines.
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
!> Contributor(s):
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

!> This module contains all generated mesh access method routines.
MODULE GeneratedMeshAccessRoutines
  
  USE BaseRoutines
  USE Kinds
  USE ISO_VARYING_STRING
  USE MeshAccessRoutines
  USE Strings
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup GeneratedMeshRoutines_GeneratedMeshTypes GeneratedMeshRoutines::GeneratedMeshTypes
  !> \brief Generated mesh types.
  !> \see GeneratedMeshRoutines,OPENCMISS_GeneratedMeshTypes
  !>@{
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_MESH_TYPE=1 !<A regular generated mesh. \see GeneratedMeshRoutines_GeneratedMeshTypes,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_POLAR_MESH_TYPE=2 !<A polar generated mesh. \see GeneratedMeshRoutines_GeneratedMeshTypes,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_FRACTAL_TREE_MESH_TYPE=3 !<A fractal tree generated mesh. \see GeneratedMeshRoutines_GeneratedMeshTypes,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_CYLINDER_MESH_TYPE=4 !<A cylinder generated mesh. \see GeneratedMeshRoutines_GeneratedMeshTypes,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_ELLIPSOID_MESH_TYPE=5 !<An ellipsoid generated mesh. \see GeneratedMeshRoutines_GeneratedMeshTypes,GeneratedMeshRoutines
  !>@}

  !> \addtogroup GeneratedMeshRoutines_GeneratedMeshCylinderSurfaces GeneratedMeshRoutines::GeneratedMeshCylinderSurfaces
  !> \brief Generated mesh cylinder type surface types.
  !>@{
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_CYLINDER_INNER_SURFACE=1  !<Inner surface of the cylinder. \see GeneratedMeshRoutines_GeneratedMeshCylinderSurfaces,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_CYLINDER_OUTER_SURFACE=2  !<Outer surface of the cylinder. \see GeneratedMeshRoutines_GeneratedMeshCylinderSurfaces,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_CYLINDER_TOP_SURFACE=3    !<Top surface of the cylinder. \see GeneratedMeshRoutines_GeneratedMeshCylinderSurfaces,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_CYLINDER_BOTTOM_SURFACE=4 !<Bottom surface of the cylinder. \see GeneratedMeshRoutines_GeneratedMeshCylinderSurfaces,GeneratedMeshRoutines
  !>@}

  !> \addtogroup GeneratedMeshRoutines_GeneratedMeshEllipsoidSurfaces GeneratedMeshRoutines::GeneratedMeshEllipsoidSurfaces
  !> \brief Generated mesh ellipsoid type surface types.
  !>@{
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_ELLIPSOID_INNER_SURFACE=5  !<Inner surface of the ellipsoid. \see GeneratedMeshRoutines_GeneratedMeshEllipsoidSurfaces,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_ELLIPSOID_OUTER_SURFACE=6  !<Outer surface of the ellipsoid. \see GeneratedMeshRoutines_GeneratedMeshEllipsoidSurfaces,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_ELLIPSOID_TOP_SURFACE=7    !<Top surface of the ellipsoid. \see GeneratedMeshRoutines_GeneratedMeshEllipsoidSurfaces,GeneratedMeshRoutines
  !>@}

  !> \addtogroup GeneratedMeshRoutines_GeneratedMeshRegularSurfaces GeneratedMeshRoutines::GeneratedMeshRegularSurfaces
  !> \brief Generated mesh regular type surface types.
  !>@{
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_LEFT_SURFACE=8    !<Left surface of the regular mesh. \see GeneratedMeshRoutines_GeneratedMeshRegularSurfaces,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_RIGHT_SURFACE=9   !<Right surface of the regular mesh. \see GeneratedMeshRoutines_GeneratedMeshRegularSurfaces,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_TOP_SURFACE=10    !<Top surface of the regular mesh. \see GeneratedMeshRoutines_GeneratedMeshRegularSurfaces,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_BOTTOM_SURFACE=11 !<Bottom surface of the regular mesh. \see GeneratedMeshRoutines_GeneratedMeshRegularSurfaces,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_FRONT_SURFACE=12  !<Front surface of the regular mesh. \see GeneratedMeshRoutines_GeneratedMeshRegularSurfaces,GeneratedMeshRoutines
  INTEGER(INTG), PARAMETER :: GENERATED_MESH_REGULAR_BACK_SURFACE=13   !<Back surface of the regular mesh. \see GeneratedMeshRoutines_GeneratedMeshRegularSurfaces,GeneratedMeshRoutines
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

  INTERFACE GeneratedMesh_UserNumberFind
    MODULE PROCEDURE GeneratedMesh_UserNumberFindInterface
    MODULE PROCEDURE GeneratedMesh_UserNumberFindRegion
  END INTERFACE GeneratedMesh_UserNumberFind

  PUBLIC GENERATED_MESH_REGULAR_MESH_TYPE,GENERATED_MESH_POLAR_MESH_TYPE,GENERATED_MESH_FRACTAL_TREE_MESH_TYPE, &
    & GENERATED_MESH_CYLINDER_MESH_TYPE,GENERATED_MESH_ELLIPSOID_MESH_TYPE
  
  PUBLIC GENERATED_MESH_CYLINDER_INNER_SURFACE,GENERATED_MESH_CYLINDER_OUTER_SURFACE,GENERATED_MESH_CYLINDER_TOP_SURFACE, &
    & GENERATED_MESH_CYLINDER_BOTTOM_SURFACE
  
  PUBLIC GENERATED_MESH_ELLIPSOID_INNER_SURFACE,GENERATED_MESH_ELLIPSOID_OUTER_SURFACE,GENERATED_MESH_ELLIPSOID_TOP_SURFACE
  
  PUBLIC GENERATED_MESH_REGULAR_LEFT_SURFACE,GENERATED_MESH_REGULAR_RIGHT_SURFACE,GENERATED_MESH_REGULAR_TOP_SURFACE, &
    & GENERATED_MESH_REGULAR_BOTTOM_SURFACE,GENERATED_MESH_REGULAR_FRONT_SURFACE,GENERATED_MESH_REGULAR_BACK_SURFACE
  
  PUBLIC GeneratedMesh_AssertIsFinished,GeneratedMesh_AssertNotFinished

  PUBLIC GeneratedMesh_BasisGet

  PUBLIC GeneratedMesh_CoordinateSystemGet

  PUBLIC GeneratedMesh_CylinderMeshGet

  PUBLIC GeneratedMesh_EllipsoidMeshGet

  PUBLIC GeneratedMesh_ExtentGet

  PUBLIC GeneratedMesh_GeneratedMeshesGet

  PUBLIC GeneratedMesh_MeshGet

  PUBLIC GeneratedMesh_NumberOfElementsGet

  PUBLIC GeneratedMesh_OriginGet

  PUBLIC GeneratedMesh_RegionGet
  
  PUBLIC GeneratedMesh_RegularMeshGet

  PUBLIC GeneratedMesh_TypeGet

  PUBLIC GeneratedMesh_UserNumberFind

  PUBLIC GeneratedMesh_UserNumberFindGeneric

  PUBLIC GeneratedMesh_UserNumberGet

  PUBLIC GeneratedMeshCylinder_BasisGet
  
  PUBLIC GeneratedMeshEllipsoid_BasisGet
  
  PUBLIC GeneratedMeshRegular_BasisGet

CONTAINS

  !
  !=================================================================================================================================
  !

  !>Assert that a generated mesh has been finished
  SUBROUTINE GeneratedMesh_AssertIsFinished(generatedMesh,err,error,*)

    !Argument Variables
    TYPE(GeneratedMeshType), POINTER, INTENT(INOUT) :: generatedMesh !<The generatedMesh to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(InterfaceType), POINTER :: interface
    TYPE(RegionType), POINTER :: parentRegion,region
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("GeneratedMesh_AssertIsFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is not associated.",err,error,*999)

    IF(.NOT.generatedMesh%generatedMeshFinished) THEN
      region=>generatedMesh%region
      IF(ASSOCIATED(region)) THEN
        localError="Generated mesh number "//TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))// &
          & " on region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))// &
          & " has not been finished."        
      ELSE
        interface=>generatedMesh%interface
        IF(ASSOCIATED(interface)) THEN
          parentRegion=>interface%parentRegion
          IF(ASSOCIATED(parentRegion)) THEN
            localError="Generated mesh number "//TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))// &
              & " on interface number "//TRIM(NumberToVString(interface%userNumber,"*",err,error))// &
              & " in parent region number "//TRIM(NumberToVString(parentRegion%userNumber,"*",err,error))// &
              & " has not been finished."        
          ELSE
            localError="Generated mesh number "//TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))// &
              & " on interface number "//TRIM(NumberToVString(interface%userNumber,"*",err,error))// &
              & " has not been finished."        
          ENDIF
        ELSE
          localError="Generated mesh number "//TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))// &
            & " has not been finished."
        ENDIF
      ENDIF
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("GeneratedMesh_AssertIsFinished")
    RETURN
999 ERRORSEXITS("GeneratedMesh_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a generated mesh has not been finished
  SUBROUTINE GeneratedMesh_AssertNotFinished(generatedMesh,err,error,*)

    !Argument Variables
    TYPE(GeneratedMeshType), POINTER, INTENT(INOUT) :: generatedMesh !<The generated mesh to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(InterfaceType), POINTER :: interface
    TYPE(RegionType), POINTER :: parentRegion,region
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("GeneratedMesh_AssertNotFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is not associated.",err,error,*999)

    IF(generatedMesh%generatedMeshFinished) THEN
      region=>generatedMesh%region
      IF(ASSOCIATED(region)) THEN
        localError="Generated mesh number "//TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))// &
          & " on region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))// &
          & " has been finished."        
      ELSE
        interface=>generatedMesh%interface
        IF(ASSOCIATED(interface)) THEN
          parentRegion=>interface%parentRegion
          IF(ASSOCIATED(parentRegion)) THEN
            localError="Generated mesh number "//TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))// &
              & " on interface number "//TRIM(NumberToVString(interface%userNumber,"*",err,error))// &
              & " in parent region number "//TRIM(NumberToVString(parentRegion%userNumber,"*",err,error))// &
              & " has been finished."        
          ELSE
            localError="Generated mesh number "//TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))// &
              & " on interface number "//TRIM(NumberToVString(interface%userNumber,"*",err,error))// &
              & " has been finished."        
          ENDIF
        ELSE
          localError="Generated mesh number "//TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))// &
            & " has been finished."
        ENDIF
      ENDIF
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("GeneratedMesh_AssertNotFinished")
    RETURN
999 ERRORSEXITS("GeneratedMesh_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Gets the basis of a generated mesh. \see OpenCMISS::Iron::cmfe_GeneratedMesh_BasisGet
  SUBROUTINE GeneratedMesh_BasisGet(generatedMesh,bases,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to get the bases of
    TYPE(BasisPtrType), POINTER :: bases(:) !<On return, the bases of mesh to generate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: basisIdx,numberOfBases
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ellipsoidMesh
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_BasisGet",err,error,*999)

    IF(ASSOCIATED(bases)) CALL FlagError("Bases is already associated.",err,error,*999)
    CALL GeneratedMesh_AssertIsFinished(generatedMesh,err,error,*999)

    SELECT CASE(generatedMesh%generatedType)
    CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
      NULLIFY(regularMesh)
      CALL GeneratedMesh_RegularMeshGet(generatedMesh,regularMesh,err,error,*999)
      IF(.NOT.ALLOCATED(regularMesh%bases)) CALL FlagError("Generated mesh bases are not allocated.",err,error,*999)
      numberOfBases=SIZE(regularMesh%bases)
      ALLOCATE(bases(numberOfBases),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate bases.",err,error,*999)
      DO basisIdx=1,numberOfBases
        bases(basisIdx)%ptr=>regularMesh%bases(basisIdx)%ptr
      ENDDO !basisIdx
    CASE(GENERATED_MESH_POLAR_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
      NULLIFY(cylinderMesh)
      CALL GeneratedMesh_CylinderMeshGet(generatedMesh,cylinderMesh,err,error,*999)
      IF(.NOT.ALLOCATED(cylinderMesh%bases)) CALL FlagError("Generated mesh bases are not allocated.",err,error,*999)
      numberOfBases=SIZE(cylinderMesh%bases)
      ALLOCATE(bases(numberOfBases),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate bases.",err,error,*999)
      DO basisIdx=1,numberOfBases
        bases(basisIdx)%ptr=>cylinderMesh%bases(basisIdx)%ptr
      ENDDO !basisIdx
    CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
      NULLIFY(ellipsoidMesh)
      CALL GeneratedMesh_EllipsoidMeshGet(generatedMesh,ellipsoidMesh,err,error,*999)
      IF(.NOT.ALLOCATED(ellipsoidMesh%bases)) CALL FlagError("Generated mesh bases are not allocated.",err,error,*999)
      numberOfBases=SIZE(generatedMesh%ellipsoidMesh%bases)
      ALLOCATE(bases(numberOfBases),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate bases.",err,error,*999)
      DO basisIdx=1,numberOfBases
        bases(basisIdx)%ptr=>ellipsoidMesh%bases(basisIdx)%ptr
      ENDDO !basisIdx
    CASE DEFAULT
      localError="The generated mesh generated type of "// &
        & TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
 
    EXITS("GeneratedMesh_BasisGet")
    RETURN
999 ERRORSEXITS("GeneratedMesh_BasisGet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_BasisGet

  !
  !================================================================================================================================
  !

  !>Returns the coordinate system for a generated mesh accounting for regions and interfaces
  SUBROUTINE GeneratedMesh_CoordinateSystemGet(generatedMesh,coordinateSystem,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to get the coordinate system for
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<On return, the generated meshes coordinate system. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(InterfaceType), POINTER :: interface
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_CoordinateSystemGet",err,error,*998)

    IF(ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is not associated.",err,error,*999)
    
    region=>generatedMesh%region
    IF(ASSOCIATED(region)) THEN
      coordinateSystem=>region%coordinateSystem
      IF(.NOT.ASSOCIATED(coordinateSystem)) THEN
        localError="The coordinate system is not associated for generated mesh number "// &
          & TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))//" of region number "// &
          & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      interface=>generatedMesh%interface
      IF(ASSOCIATED(interface)) THEN
        coordinateSystem=>interface%coordinateSystem
        IF(.NOT.ASSOCIATED(coordinateSystem)) THEN
          localError="The coordinate system is not associated for generated mesh number "// &
            & TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))//" of interface number "// &
            & TRIM(NumberToVString(interface%userNumber,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        localError="The interface is not associated for generated mesh number "// &
          & TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF

    EXITS("GeneratedMesh_CoordinateSystemGet")
    RETURN
999 NULLIFY(coordinateSystem)
998 ERRORSEXITS("GeneratedMesh_CoordinateSystemGet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_CoordinateSystemGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the cylinder mesh for a given generated mesh.
  SUBROUTINE GeneratedMesh_CylinderMeshGet(generatedMesh,cylinderMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to get the cylinder mesh for
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh !<On exit, a pointer to the cylinder mesh for the generated mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("GeneratedMesh_CylinderMeshGet",err,error,*998)

    IF(ASSOCIATED(cylinderMesh)) CALL FlagError("Cylinder mesh is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is not associated.",err,error,*998)
    
    cylinderMesh=>generatedMesh%cylinderMesh
    IF(.NOT.ASSOCIATED(cylinderMesh)) THEN
      localError="A cylinder mesh is not associated on generated mesh number "// &
        & TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("GeneratedMesh_CylinderMeshGet")
    RETURN
999 NULLIFY(cylinderMesh)
998 ERRORSEXITS("GeneratedMesh_CylinderMeshGet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_CylinderMeshGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the ellipsoid mesh for a given generated mesh.
  SUBROUTINE GeneratedMesh_EllipsoidMeshGet(generatedMesh,ellipsoidMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to get the ellipsoid mesh for
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ellipsoidMesh !<On exit, a pointer to the ellipsoid mesh for the generated mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("GeneratedMesh_EllipsoidMeshGet",err,error,*998)

    IF(ASSOCIATED(ellipsoidMesh)) CALL FlagError("Ellipsoid mesh is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is not associated.",err,error,*998)
    
    ellipsoidMesh=>generatedMesh%ellipsoidMesh
    IF(.NOT.ASSOCIATED(ellipsoidMesh)) THEN
      localError="An ellipsoid mesh is not associated on generated mesh number "// &
        & TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("GeneratedMesh_EllipsoidMeshGet")
    RETURN
999 NULLIFY(ellipsoidMesh)
998 ERRORSEXITS("GeneratedMesh_EllipsoidMeshGet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_EllipsoidMeshGet

  !
  !================================================================================================================================
  !

  !>Gets the extent of a generated mesh. \see OpenCMISS::Iron::cmfe_GeneratedMesh_ExtentGet
  SUBROUTINE GeneratedMesh_ExtentGet(generatedMesh,extent,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to get the type of
    REAL(DP), INTENT(OUT) :: extent(:) !<On return, maximum extent per axis, or inner & outer radii and length of cylinder
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ellipsoidMesh
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_ExtentGet",err,error,*999)

    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is not associated",err,error,*999)
      
    SELECT CASE(generatedMesh%generatedType)
    CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
      NULLIFY(regularMesh)
      CALL GeneratedMesh_RegularMeshGet(generatedMesh,regularMesh,err,error,*999)      
      IF(SIZE(extent,1)<SIZE(regularMesh%maximumExtent,1)) THEN
        localError="The size of extent is too small. The supplied size is "// &
          & TRIM(NumberToVString(SIZE(extent,1),"*",err,error))//" and it needs to be >= "// &
          & TRIM(NumberToVString(SIZE(regularMesh%maximumExtent,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      extent(1:SIZE(regularMesh%maximumExtent,1))=regularMesh%maximumExtent(1:SIZE(regularMesh%maximumExtent,1))
    CASE(GENERATED_MESH_POLAR_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
      NULLIFY(cylinderMesh)
      CALL GeneratedMesh_CylinderMeshGet(generatedMesh,cylinderMesh,err,error,*999)
      IF(SIZE(extent,1)<SIZE(cylinderMesh%cylinderExtent,1)) THEN
        localError="The size of extent is too small. The supplied size is "// &
          & TRIM(NumberToVString(SIZE(extent,1),"*",err,error))//" and it needs to be >= "// &
          & TRIM(NumberToVString(SIZE(cylinderMesh%cylinderExtent,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      extent(1:SIZE(cylinderMesh%cylinderExtent,1))=cylinderMesh%cylinderExtent(1:SIZE(cylinderMesh%cylinderExtent,1))
    CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
      NULLIFY(ellipsoidMesh)
      CALL GeneratedMesh_EllipsoidMeshGet(generatedMesh,ellipsoidMesh,err,error,*999)
      IF(SIZE(extent,1)<SIZE(ellipsoidMesh%ellipsoidExtent,1)) THEN
        localError="The size of extent is too small. The supplied size is "// &
          & TRIM(NumberToVString(SIZE(extent,1),"*",err,error))//" and it needs to be >= "// &
          & TRIM(NumberToVString(SIZE(ellipsoidMesh%ellipsoidExtent,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      extent(1:SIZE(ellipsoidMesh%ellipsoidExtent,1))=ellipsoidMesh%ellipsoidExtent(1:SIZE(ellipsoidMesh%ellipsoidExtent,1))
    CASE DEFAULT
      localError="The generated mesh mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("GeneratedMesh_ExtentGet")
    RETURN
999 ERRORSEXITS("GeneratedMesh_ExtentGet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_ExtentGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the generated meshes for a given generated mesh.
  SUBROUTINE GeneratedMesh_GeneratedMeshesGet(generatedMesh,generatedMeshes,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to get the generated meshes for
    TYPE(GeneratedMeshesType), POINTER :: generatedMeshes !<On exit, a pointer to the generated meshes for the generated mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("GeneratedMesh_GeneratedMeshesGet",err,error,*998)

    IF(ASSOCIATED(generatedMeshes)) CALL FlagError("Generated meshes is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is not associated.",err,error,*998)
    
    generatedMeshes=>generatedMesh%generatedMeshes
    IF(.NOT.ASSOCIATED(generatedMeshes)) THEN
      localError="The generated meshes is not associated for generated mesh number "// &
        & TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("GeneratedMesh_GeneratedMeshesMeshesGet")
    RETURN
999 NULLIFY(generatedMeshes)
998 ERRORSEXITS("GeneratedMesh_GeneratedMesheshGet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_GeneratedMeshesGet

  !
  !================================================================================================================================
  !

  !>Returns the mesh for a generated mesh.
  SUBROUTINE GeneratedMesh_MeshGet(generatedMesh,mesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to get the mesh for
    TYPE(MeshType), POINTER :: mesh !<On return, the generated meshes mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_MeshGet",err,error,*998)

    IF(ASSOCIATED(mesh)) CALL FlagError("Mesh is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is not associated.",err,error,*999)
    
    mesh=>generatedMesh%mesh
    IF(ASSOCIATED(mesh)) THEN
      localError="The mesh is not associated for generated mesh number "// &
        & TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("GeneratedMesh_MeshGet")
    RETURN
999 NULLIFY(mesh)
998 ERRORSEXITS("GeneratedMesh_MeshGet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_MeshGet

  !
  !================================================================================================================================
  !

  !>Gets the number of elements in a generated mesh.  \see OpenCMISS::Iron::cmfe_GeneratedMesh_NumberOfElementsGet
  SUBROUTINE GeneratedMesh_NumberOfElementsGet(generatedMesh,numberOfElements,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to set the type of
    INTEGER(INTG), INTENT(OUT) :: numberOfElements(:) !<On return, number of elements per axis
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ellipsoidMesh
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_NumberOfElementsGet",err,error,*999)

    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is not associated.",err,error,*999)
    
    SELECT CASE(generatedMesh%generatedType)
    CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
      NULLIFY(regularMesh)
      CALL GeneratedMesh_RegularMeshGet(generatedMesh,regularMesh,err,error,*999)
      IF(SIZE(numberOfElements,1)<SIZE(regularMesh%numberOfElementsXi,1)) THEN
        localError="The size of the specified number of elements array is too small. The supplied size is "// &
          & TRIM(NumberToVString(SIZE(numberOfElements,1),"*",err,error))//" and it needs to be >= "// &
          & TRIM(NumberToVString(SIZE(regularMesh%numberOfElementsXi,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      numberOfElements(1:SIZE(regularMesh%numberOfElementsXi,1))= &
        & regularMesh%numberOfElementsXi(1:SIZE(regularMesh%numberOfElementsXi,1))
    CASE(GENERATED_MESH_POLAR_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
      NULLIFY(cylinderMesh)
      CALL GeneratedMesh_CylinderMeshGet(generatedMesh,cylinderMesh,err,error,*999)
      IF(SIZE(numberOfElements,1)<SIZE(cylinderMesh%numberOfElementsXi,1)) THEN
        localError="The size of the specified number of elements array is too small. The supplied size is "// &
          & TRIM(NumberToVString(SIZE(numberOfElements,1),"*",err,error))//" and it needs to be >= "// &
          & TRIM(NumberToVString(SIZE(cylinderMesh%numberOfElementsXi,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      numberOfElements(1:SIZE(cylinderMesh%numberOfElementsXi,1))= &
        & cylinderMesh%numberOfElementsXi(1:SIZE(cylinderMesh%numberOfElementsXi,1))
    CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
      NULLIFY(ellipsoidMesh)
      CALL GeneratedMesh_EllipsoidMeshGet(generatedMesh,ellipsoidMesh,err,error,*999)
      IF(SIZE(numberOfElements,1)<SIZE(generatedMesh%ellipsoidMesh%numberOfElementsXi,1)) THEN
        localError="The size of the specified number of elements array is too small. The supplied size is "// &
          & TRIM(NumberToVString(SIZE(numberOfElements,1),"*",err,error))//" and it needs to be >= "// &
          & TRIM(NumberToVString(SIZE(ellipsoidMesh%numberOfElementsXi,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      numberOfElements(1:SIZE(ellipsoidMesh%numberOfElementsXi,1))= &
        & ellipsoidMesh%numberOfElementsXi(1:SIZE(ellipsoidMesh%numberOfElementsXi,1))
    CASE DEFAULT
      localError="The generated mesh mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("GeneratedMesh_NumberOfElementsGet")
    RETURN
999 ERRORSEXITS("GeneratedMesh_NumberOfElementsGet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_NumberOfElementsGet

  !
  !================================================================================================================================
  !

  !>Get the origin of a generated mesh. \see OpenCMISS::Iron::cmfe_GeneratedMesh_OriginGet
  SUBROUTINE GeneratedMesh_OriginGet(generatedMesh,origin,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to get the type of
    REAL(DP), INTENT(OUT) :: origin(:) !<On return, the origin coordinate for each axis
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ellipsoidMesh
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_OriginGet",err,error,*999)

    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is not associated",err,error,*999)
    
    SELECT CASE(generatedMesh%generatedType)
    CASE(GENERATED_MESH_REGULAR_MESH_TYPE)
      NULLIFY(regularMesh)
      CALL GeneratedMesh_RegularMeshGet(generatedMesh,regularMesh,err,error,*999)
      IF(SIZE(origin,1)<SIZE(regularMesh%origin,1)) THEN
        localError="The size of origin is too small. The supplied size is "// &
          & TRIM(NumberToVString(SIZE(origin,1),"*",err,error))//" and it needs to be >= "// &
          & TRIM(NumberToVString(SIZE(regularMesh%origin,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      origin(1:SIZE(regularMesh%origin,1))=regularMesh%origin(1:SIZE(regularMesh%origin,1))
    CASE(GENERATED_MESH_POLAR_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_FRACTAL_TREE_MESH_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(GENERATED_MESH_CYLINDER_MESH_TYPE)
      NULLIFY(cylinderMesh)
      CALL GeneratedMesh_CylinderMeshGet(generatedMesh,cylinderMesh,err,error,*999)
      IF(SIZE(origin,1) < SIZE(cylinderMesh%origin,1)) THEN
        localError="The size of origin is too small. The supplied size is "// &
          & TRIM(NumberToVString(SIZE(origin,1),"*",err,error))//" and it needs to be >= "// &
          & TRIM(NumberToVString(SIZE(cylinderMesh%origin,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      origin(1:SIZE(cylinderMesh%origin,1))=cylinderMesh%origin(1:SIZE(cylinderMesh%origin,1))
    CASE(GENERATED_MESH_ELLIPSOID_MESH_TYPE)
      NULLIFY(ellipsoidMesh)
      CALL GeneratedMesh_EllipsoidMeshGet(generatedMesh,ellipsoidMesh,err,error,*999)
      IF(SIZE(origin,1)<SIZE(ellipsoidMesh%origin,1)) THEN
        localError="The size of origin is too small. The supplied size is "// &
          & TRIM(NumberToVString(SIZE(origin,1),"*",err,error))//" and it needs to be >= "// &
          & TRIM(NumberToVString(SIZE(ellipsoidMesh%origin,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      origin(1:SIZE(ellipsoidMesh%origin,1))=ellipsoidMesh%origin(1:SIZE(ellipsoidMesh%origin,1))
    CASE DEFAULT
      localError="The generated mesh mesh type of "//TRIM(NumberToVString(generatedMesh%generatedType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("GeneratedMesh_OriginGet")
    RETURN
999 ERRORSEXITS("GeneratedMesh_OriginGet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_OriginGet

  !
  !================================================================================================================================
  !

  !>Returns the region for a generated mesh accounting for regions and interfaces
  SUBROUTINE GeneratedMesh_RegionGet(generatedMesh,region,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to get the region for
    TYPE(RegionType), POINTER :: region !<On return, the generated meshes region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(InterfaceType), POINTER :: interface
    TYPE(RegionType), POINTER :: parentRegion
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_RegionGet",err,error,*998)

    IF(ASSOCIATED(region)) CALL FlagError("Region is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is not associated.",err,error,*999)
    
    region=>generatedMesh%region
    IF(.NOT.ASSOCIATED(region)) THEN
      interface=>generatedMesh%interface
      IF(ASSOCIATED(interface)) THEN
        parentRegion=>interface%parentRegion
        IF(ASSOCIATED(parentRegion)) THEN
          region=>parentRegion
        ELSE
          localError="The parent region not associated for generated mesh number "// &
            & TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))//" of interface number "// &
            & TRIM(NumberToVString(interface%userNumber,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        localError="The region or interface is not associated for generated mesh number "// &
          & TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF

    EXITS("GeneratedMesh_RegionGet")
    RETURN
999 NULLIFY(region)
998 ERRORSEXITS("GeneratedMesh_RegionGet",err,error)    
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_RegionGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the regular mesh for a given generated mesh.
  SUBROUTINE GeneratedMesh_RegularMeshGet(generatedMesh,regularMesh,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to get the cylinder mesh for
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh !<On exit, a pointer to the regular mesh for the generated mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("GeneratedMesh_RegularMeshGet",err,error,*998)

    IF(ASSOCIATED(regularMesh)) CALL FlagError("Regular mesh is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is not associated.",err,error,*998)
    
    regularMesh=>generatedMesh%regularMesh
    IF(.NOT.ASSOCIATED(regularMesh)) THEN
      localError="A regular mesh is not associated on generated mesh number "// &
        & TRIM(NumberToVString(generatedMesh%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("GeneratedMesh_RegularMeshGet")
    RETURN
999 NULLIFY(regularMesh)
998 ERRORSEXITS("GeneratedMesh_RegularMeshGet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_RegularMeshGet

  !
  !================================================================================================================================
  !

  !>Gets the type of a generated mesh. \see OpenCMISS::Iron::cmfe_GeneratedMesh_TypeGet
  SUBROUTINE GeneratedMesh_TypeGet(generatedMesh,type,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh to set the type of
    INTEGER(INTG), INTENT(OUT) :: TYPE !<On return, the type of mesh to generate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("GeneratedMesh_TypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is not associated.",err,error,*999)
    
    type=generatedMesh%generatedType

    EXITS("GeneratedMesh_TypeGet")
    RETURN
999 ERRORSEXITS("GeneratedMesh_TypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_TypeGet

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the generated mesh identified by user number in the given list of generated meshes. If no generated mesh with that number exits generated mesh is left nullified.
  SUBROUTINE GeneratedMesh_UserNumberFindGeneric(userNumber,generatedMeshes,generatedMesh,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the generated mesh to find
    TYPE(GeneratedMeshesType), POINTER :: generatedMeshes !<The list of generated meshes containing the generated mesh.
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<On return, a pointer to the generated mesh of the specified user number. In no generated mesh with the specified user number exists the pointer is returned NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: generatedMeshIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("GeneratedMesh_UserNumberFindGeneric",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(generatedMeshes)) CALL FlagError("Generated meshes is not associated",err,error,*999)
    IF(ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is already associated.",err,error,*999)
    
    !Get the generated mesh from the user number
    NULLIFY(generatedMesh)
    IF(ALLOCATED(generatedMeshes%generatedMeshes)) THEN
      DO generatedMeshIdx=1,generatedMeshes%numberOfGeneratedMeshes
        IF(ASSOCIATED(generatedMeshes%generatedMeshes(generatedMeshIdx)%ptr)) THEN
          IF(generatedMeshes%generatedMeshes(generatedMeshIdx)%ptr%userNumber==userNumber) THEN
            generatedMesh=>generatedMeshes%generatedMeshes(generatedMeshIdx)%ptr
            EXIT
          ENDIF
        ELSE
          localError="The generated mesh pointer in generated meshes is not associated for generated mesh index "// &
            & TRIM(NumberToVString(generatedMeshIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)          
        ENDIF
      ENDDO !generatedMeshIdx      
    ENDIF
    
    EXITS("GeneratedMesh_UserNumberFindGeneric")
    RETURN
999 ERRORSEXITS("GeneratedMesh_UserNumberFindGeneric",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_UserNumberFindGeneric

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the generated mesh identified by user number in the given interface. If no generated mesh with that number exists generated mesh is left nullified.
  SUBROUTINE GeneratedMesh_UserNumberFindInterface(userNumber,interface,generatedMesh,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the generated mesh to find
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface containing the generated mesh
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<On return, a pointer to the generated mesh of the specified user number. In no generated mesh with the specified user number exists the pointer is returned NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
  
    ENTERS("GeneratedMesh_UserNumberFindInterface",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated",err,error,*999)
    
    CALL GeneratedMesh_UserNumberFindGeneric(userNumber,interface%generatedMeshes,generatedMesh,err,error,*999)
     
    EXITS("GeneratedMesh_UserNumberFindInterface")
    RETURN
999 ERRORSEXITS("GeneratedMesh_UserNumberFindInterface",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_UserNumberFindInterface

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the generated mesh identified by user number in the given. If no generated mesh with that number exists generated mesh is left nullified.
  SUBROUTINE GeneratedMesh_UserNumberFindRegion(userNumber,region,generatedMesh,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the generated mesh to find
    TYPE(RegionType), POINTER :: region !<The region containing the generated  mesh
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<On return, a pointer to the generated mesh of the specified user number. In no generated mesh with the specified user number exists the pointer is returned NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("GeneratedMesh_UserNumberFindRegion",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated",err,error,*999)
    
    CALL GeneratedMesh_UserNumberFindGeneric(userNumber,region%generatedMeshes,generatedMesh,err,error,*999)
     
    EXITS("GeneratedMesh_UserNumberFindRegion")
    RETURN
999 ERRORSEXITS("GeneratedMesh_UserNumberFindRegion",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_UserNumberFindRegion

  !
  !================================================================================================================================
  !

  !>Returns the user number for a generatedMesh.
  SUBROUTINE GeneratedMesh_UserNumberGet(generatedMesh,userNumber,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generatedMesh to get the user number for
    INTEGER(INTG), INTENT(OUT) :: userNumber !<On exit, the user number of the generatedMesh.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("GeneratedMesh_UserNumberGet",err,error,*999)

    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("GeneratedMesh is not associated.",err,error,*999)

    userNumber=generatedMesh%userNumber
  
    EXITS("GeneratedMesh_UserNumberGet")
    RETURN
999 ERRORSEXITS("GeneratedMesh_UserNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_UserNumberGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the basis for a cylinder mesh for a given basis index.
  SUBROUTINE GeneratedMeshCylinder_BasisGet(cylinderMesh,basisIndex,basis,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh !<A pointer to the cylinder mesh to get the basis for.
    INTEGER(INTG), INTENT(IN) :: basisIndex !<The index of the basis in the cylinder mesh to get.
    TYPE(BasisType), POINTER :: basis !<On exit, the cylinder mesh basis. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("GeneratedMeshCylinder_BasisGet",err,error,*998)

    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cylinderMesh)) CALL FlagError("Cylinder mesh is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(cylinderMesh%bases)) CALL FlagError("Cylinder mesh bases are not allocated.",err,error,*999)
    IF(basisIndex<1.OR.basisIndex>SIZE(cylinderMesh%bases,1)) THEN
      localError="The specified basis index of "//TRIM(NumberToVString(basisIndex,"*",err,error))// &
        & " is invalid. The basis index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(cylinderMesh%bases,1),"*",err,error))
      IF(ASSOCIATED(cylinderMesh%generatedMesh)) &
        & localError=localError//" for the cylinder mesh of generated mesh number "// &
        & TRIM(NumberToVString(cylinderMesh%generatedMesh%userNumber,"*",err,error))
      localError=localError//"."
    ENDIF
    
    basis=>cylinderMesh%bases(basisIndex)%ptr
    IF(.NOT.ASSOCIATED(basis)) THEN
      localError="The basis corresponding to basis index "//TRIM(NumberToVString(basisIndex,"*",err,error))
      IF(ASSOCIATED(cylinderMesh%generatedMesh)) &
        & localError=localError//" for the cylinder mesh of generated mesh number "// &
        & TRIM(NumberToVString(cylinderMesh%generatedMesh%userNumber,"*",err,error))
      localError=localError//" is not associated."
    ENDIF
    
    EXITS("GeneratedMeshCylinder_BasisGet")
    RETURN
999 NULLIFY(basis)
998 ERRORSEXITS("GeneratedMeshCylinder_BasisGet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMeshCylinder_BasisGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the basis for a ellipsoid mesh for a given basis index.
  SUBROUTINE GeneratedMeshEllipsoid_BasisGet(ellipsoidMesh,basisIndex,basis,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ellipsoidMesh !<A pointer to the ellipsoid mesh to get the basis for.
    INTEGER(INTG), INTENT(IN) :: basisIndex !<The index of the basis in the ellipsoid mesh to get.
    TYPE(BasisType), POINTER :: basis !<On exit, the regular mesh basis. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("GeneratedMeshEllipsoid_BasisGet",err,error,*998)

    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(ellipsoidMesh)) CALL FlagError("Ellipsoid mesh is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(ellipsoidMesh%bases)) CALL FlagError("Ellipsoid mesh bases are not allocated.",err,error,*999)
    IF(basisIndex<1.OR.basisIndex>SIZE(ellipsoidMesh%bases,1)) THEN
      localError="The specified basis index of "//TRIM(NumberToVString(basisIndex,"*",err,error))// &
        & " is invalid. The basis index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(ellipsoidMesh%bases,1),"*",err,error))
      IF(ASSOCIATED(ellipsoidMesh%generatedMesh)) &
        & localError=localError//" for the ellipsoid mesh of generated mesh number "// &
        & TRIM(NumberToVString(ellipsoidMesh%generatedMesh%userNumber,"*",err,error))
      localError=localError//"."
    ENDIF
    
    basis=>ellipsoidMesh%bases(basisIndex)%ptr
    IF(.NOT.ASSOCIATED(basis)) THEN
      localError="The basis corresponding to basis index "//TRIM(NumberToVString(basisIndex,"*",err,error))
      IF(ASSOCIATED(ellipsoidMesh%generatedMesh)) &
        & localError=localError//" for the ellipsoid mesh of generated mesh number "// &
        & TRIM(NumberToVString(ellipsoidMesh%generatedMesh%userNumber,"*",err,error))
      localError=localError//" is not associated."
    ENDIF
    
    EXITS("GeneratedMeshEllipsoid_BasisGet")
    RETURN
999 NULLIFY(basis)
998 ERRORSEXITS("GeneratedMeshEllipsoid_BasisGet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMeshEllipsoid_BasisGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the basis for a regular mesh for a given basis index.
  SUBROUTINE GeneratedMeshRegular_BasisGet(regularMesh,basisIndex,basis,err,error,*)

    !Argument variables
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh !<A pointer to the regular mesh to get the basis for.
    INTEGER(INTG), INTENT(IN) :: basisIndex !<The index of the basis in the regular mesh to get.
    TYPE(BasisType), POINTER :: basis !<On exit, the regular mesh basis. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("GeneratedMeshRegular_BasisGet",err,error,*998)

    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(regularMesh)) CALL FlagError("Regular mesh is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(regularMesh%bases)) CALL FlagError("Regular mesh bases are not allocated.",err,error,*999)
    IF(basisIndex<1.OR.basisIndex>SIZE(regularMesh%bases,1)) THEN
      localError="The specified basis index of "//TRIM(NumberToVString(basisIndex,"*",err,error))// &
        & " is invalid. The basis index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(regularMesh%bases,1),"*",err,error))
      IF(ASSOCIATED(regularMesh%generatedMesh)) &
        & localError=localError//" for the regular mesh of generated mesh number "// &
        & TRIM(NumberToVString(regularMesh%generatedMesh%userNumber,"*",err,error))
      localError=localError//"."
    ENDIF
    
    basis=>regularMesh%bases(basisIndex)%ptr
    IF(.NOT.ASSOCIATED(basis)) THEN
      localError="The basis corresponding to basis index "//TRIM(NumberToVString(basisIndex,"*",err,error))
      IF(ASSOCIATED(regularMesh%generatedMesh)) &
        & localError=localError//" for the regular mesh of generated mesh number "// &
        & TRIM(NumberToVString(regularMesh%generatedMesh%userNumber,"*",err,error))
      localError=localError//" is not associated."
    ENDIF
    
    EXITS("GeneratedMeshRegular_BasisGet")
    RETURN
999 NULLIFY(basis)
998 ERRORSEXITS("GeneratedMeshRegular_BasisGet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMeshRegular_BasisGet

  !
  !================================================================================================================================
  !

END MODULE GeneratedMeshAccessRoutines
