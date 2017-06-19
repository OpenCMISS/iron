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
  USE MeshAccessRoutines
  USE Strings
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  INTERFACE GeneratedMesh_UserNumberFind
    MODULE PROCEDURE GeneratedMesh_UserNumberFindInterface
    MODULE PROCEDURE GeneratedMesh_UserNumberFindRegion
  END INTERFACE GeneratedMesh_UserNumberFind

  INTERFACE GENERATED_MESH_USER_NUMBER_FIND
    MODULE PROCEDURE GeneratedMesh_UserNumberFindInterface
    MODULE PROCEDURE GeneratedMesh_UserNumberFindRegion
  END INTERFACE GENERATED_MESH_USER_NUMBER_FIND

  PUBLIC GeneratedMesh_UserNumberFind

  PUBLIC GeneratedMesh_UserNumberFindGeneric

  PUBLIC GENERATED_MESH_USER_NUMBER_FIND

  PUBLIC GeneratedMesh_UserNumberGet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the generated mesh identified by user number in the given list of generated meshes. If no generated mesh with that number exits generated mesh is left nullified.
  SUBROUTINE GeneratedMesh_UserNumberFindGeneric(userNumber,generatedMeshes,generatedMesh,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the generated mesh to find
    TYPE(GeneratedMeshesType), POINTER :: generatedMeshes !<The list of generated meshes containing the generated mesh.
    TYPE(GENERATED_MESH_TYPE), POINTER :: generatedMesh !<On return, a pointer to the generated mesh of the specified user number. In no generated mesh with the specified user number exists the pointer is returned NULL. Must not be associated on entry.
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
    IF(ASSOCIATED(generatedMeshes%generatedMeshes)) THEN
      DO generatedMeshIdx=1,generatedMeshes%numberOfGeneratedMeshes
        IF(ASSOCIATED(generatedMeshes%generatedMeshes(generatedMeshIdx)%ptr)) THEN
          IF(generatedMeshes%generatedMeshes(generatedMeshIdx)%ptr%USER_NUMBER==userNumber) THEN
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
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface containing the generated mesh
    TYPE(GENERATED_MESH_TYPE), POINTER :: generatedMesh !<On return, a pointer to the generated mesh of the specified user number. In no generated mesh with the specified user number exists the pointer is returned NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
  
    ENTERS("GeneratedMesh_UserNumberFindInterface",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(INTERFACE)) CALL FlagError("Interface is not associated",err,error,*999)
    
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
    TYPE(REGION_TYPE), POINTER :: region !<The region containing the generated  mesh
    TYPE(GENERATED_MESH_TYPE), POINTER :: generatedMesh !<On return, a pointer to the generated mesh of the specified user number. In no generated mesh with the specified user number exists the pointer is returned NULL. Must not be associated on entry.
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
    TYPE(GENERATED_MESH_TYPE), POINTER :: generatedMesh !<A pointer to the generatedMesh to get the user number for
    INTEGER(INTG), INTENT(OUT) :: userNumber !<On exit, the user number of the generatedMesh.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("GeneratedMesh_UserNumberGet",err,error,*999)

    IF(.NOT.ASSOCIATED(generatedMesh)) CALL FlagError("GeneratedMesh is not associated.",err,error,*999)

    userNumber=generatedMesh%USER_NUMBER
  
    EXITS("GeneratedMesh_UserNumberGet")
    RETURN
999 ERRORSEXITS("GeneratedMesh_UserNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE GeneratedMesh_UserNumberGet

  !
  !================================================================================================================================
  !

END MODULE GeneratedMeshAccessRoutines
