!> \file
!> \author Chris Bradley
!> \brief This module contains all mesh access method routines.
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

!> This module contains all mesh access method routines.
MODULE MeshAccessRoutines
  
  USE BaseRoutines
  USE Kinds
  USE Strings
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  INTERFACE DECOMPOSITION_USER_NUMBER_FIND
    MODULE PROCEDURE Decomposition_UserNumberFind
  END INTERFACE DECOMPOSITION_USER_NUMBER_FIND

  INTERFACE MESH_TOPOLOGY_ELEMENTS_GET
    MODULE PROCEDURE Mesh_MeshElementsGet
  END INTERFACE MESH_TOPOLOGY_ELEMENTS_GET

  INTERFACE MeshTopology_NodesGet
    MODULE PROCEDURE Mesh_MeshNodesGet
  END INTERFACE MeshTopology_NodesGet

  INTERFACE Mesh_UserNumberFind
    MODULE PROCEDURE Mesh_UserNumberFindInterface
    MODULE PROCEDURE Mesh_UserNumberFindRegion
  END INTERFACE Mesh_UserNumberFind

   INTERFACE MESH_USER_NUMBER_FIND
    MODULE PROCEDURE Mesh_UserNumberFindInterface
    MODULE PROCEDURE Mesh_UserNumberFindRegion
  END INTERFACE MESH_USER_NUMBER_FIND

  PUBLIC Decomposition_CoordinateSystemGet

  PUBLIC Decomposition_DomainGet

  PUBLIC Decomposition_TopologyGet
  
  PUBLIC Decomposition_UserNumberFind

  PUBLIC DECOMPOSITION_USER_NUMBER_FIND

  PUBLIC DecompositionTopology_DataPointsGet

  PUBLIC DecompositionTopology_ElementsGet

  PUBLIC DecompositionTopology_FacesGet
  
  PUBLIC DecompositionTopology_LinesGet

  PUBLIC Domain_MappingsGet

  PUBLIC Domain_TopologyGet

  PUBLIC DomainElements_BasisGet

  PUBLIC DomainFaces_FaceGet

  PUBLIC DomainLines_LineGet
  
  PUBLIC DomainMappings_ElementsGet

  PUBLIC DomainMappings_NodesGet

  PUBLIC DomainTopology_ElementsGet

  PUBLIC DomainTopology_FacesGet

  PUBLIC DomainTopology_LinesGet
  
  PUBLIC DomainTopology_NodesGet
  
  PUBLIC Mesh_DecompositionGet

  PUBLIC Mesh_MeshElementsGet

  PUBLIC Mesh_MeshNodesGet

  PUBLIC Mesh_RegionGet

  PUBLIC Mesh_UserNumberFind

  PUBLIC Mesh_UserNumberFindGeneric

  PUBLIC MESH_USER_NUMBER_FIND

  PUBLIC Mesh_UserNumberGet

  PUBLIC MESH_TOPOLOGY_ELEMENTS_GET

  PUBLIC MeshTopology_NodesGet


CONTAINS

  !
  !================================================================================================================================
  !

  !>Returns the coordinate system for a decomposition accounting for regions and interfaces. 
  SUBROUTINE Decomposition_CoordinateSystemGet(decomposition,coordinateSystem,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition !<The decomposition to get the coordinate system for.
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem !<On exit, a pointer to the coordinate system of decomposition. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Decomposition_CoordinateSystemGet",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
    IF(ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is already associated.",err,error,*999)

    IF(ASSOCIATED(decomposition%region)) THEN
      coordinateSystem=>decomposition%region%COORDINATE_SYSTEM
    ELSE IF(ASSOCIATED(decomposition%INTERFACE)) THEN
      coordinateSystem=>decomposition%interface%COORDINATE_SYSTEM
    ELSE
      CALL FlagError("Decomposition is not associated with a region or an interface.",err,error,*999)
    ENDIF
 
    !Check coordinate system is associated.
    IF(.NOT.ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is not associated.",err,error,*999)
    
    EXITS("Decomposition_CoordinateSystemGet")
    RETURN
999 ERRORSEXITS("Decomposition_CoordinateSystemGet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_CoordinateSystemGet

  !
  !================================================================================================================================
  !

  !>Gets a domain from a decomposition and a mesh component. If mesh component is 0 then the mesh component used for the decomposition is used. 
  SUBROUTINE Decomposition_DomainGet(decomposition,meshComponentNumber,domain,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition !<The decomposition to get the domain for.
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component number of the domain to get. If mesh component number is 0 the mesh component used to construct the decomposition is used.
    TYPE(DOMAIN_TYPE), POINTER :: domain !<On exit, a pointer to the domain of decomposition of the specified mesh component number. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshComponent
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Decomposition_DomainGet",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
    IF(meshComponentNumber<0.OR.meshComponentNumber>decomposition%numberOfComponents) THEN
      localError="The specified mesh component number of "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
        & " is invalid. The mesh component number must be >= 0 and <= "// &
        & TRIM(NumberToVString(decomposition%numberOfComponents,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(ASSOCIATED(domain)) CALL FlagError("Domain is already associated.",err,error,*999)

    !Get the domain
    IF(meshComponentNumber==0) THEN
      meshComponent=decomposition%MESH_COMPONENT_NUMBER
    ELSE
      meshComponent=meshComponentNumber
    ENDIF
    IF(.NOT.ASSOCIATED(decomposition%domain)) CALL FlagError("Decomposition domain is not associated.",err,error,*999)
    domain=>decomposition%domain(meshComponent)%ptr
    
    !Check domain is associated.
    IF(.NOT.ASSOCIATED(domain)) THEN
      localError="Domain is not associated for mesh component "//TRIM(NumberToVString(meshComponent,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Decomposition_DomainGet")
    RETURN
999 ERRORSEXITS("Decomposition_DomainGet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_DomainGet

  !
  !================================================================================================================================
  !

  !>Gets a decomposition topology from a decomposition.
  SUBROUTINE Decomposition_TopologyGet(decomposition,decompositionTopology,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition !<The decomposition to get the domain for.
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology !<On exit, a pointer to the decomposition topology. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Decomposition_TopologyGet",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
    IF(ASSOCIATED(decompositionTopology)) CALL FlagError("Decomposition topology is already associated.",err,error,*999)
 
    !Get the decomposition topology
    decompositionTopology=>decomposition%topology
    !Check decompositionTopology is associated.
    IF(.NOT.ASSOCIATED(decompositionTopology)) THEN
      localError="Decomposition topology is not associated for decomposition "// &
        & TRIM(NumberToVString(decomposition%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Decomposition_TopologyGet")
    RETURN
999 ERRORSEXITS("Decomposition_TopologyGet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_TopologyGet

  !
  !================================================================================================================================
  !

  !>Gets elements from a decomposition topology.
  SUBROUTINE DecompositionTopology_ElementsGet(decompositionTopology,decompositionElements,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology !<A pointer to the decomposition topology to get the elements for.
     TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: decompositionElements !<On exit, a pointer to the decomposition topology elements. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DecompositionTopology_ElementsGet",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(decompositionTopology)) CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    IF(ASSOCIATED(decompositionElements)) CALL FlagError("Decomposition elements is already associated.",err,error,*999)

    !Get the decomposition elements
    decompositionElements=>decompositionTopology%elements

    !Check decomposition elements is associated.
    IF(.NOT.ASSOCIATED(decompositionElements)) CALL FlagError("Decomposition topology elements is not associated.",err,error,*999)
    
    EXITS("DecompositionTopology_ElementsGet")
    RETURN
999 ERRORSEXITS("DecompositionTopology_ElementsGet",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_ElementsGet

  !
  !================================================================================================================================
  !

  !>Gets data points from a decomposition topology.
  SUBROUTINE DecompositionTopology_DataPointsGet(decompositionTopology,decompositionDataPoints,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology !<A pointer to the decomposition topology to get the data points for.
    TYPE(DecompositionDataPointsType), POINTER :: decompositionDataPoints !<On exit, a pointer to the decomposition topology data points. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DecompositionTopology_DataPointsGet",err,error,*998)

    !Check input arguments
    IF(ASSOCIATED(decompositionDataPoints)) CALL FlagError("Decomposition data points is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(decompositionTopology)) CALL FlagError("Decomposition topology is not associated.",err,error,*999)

    !Get the decomposition data points
    decompositionDataPoints=>decompositionTopology%dataPoints
    !Check decomposition data points is associated.
    IF(.NOT.ASSOCIATED(decompositionDataPoints)) &
      & CALL FlagError("Decomposition topology data points is not associated.",err,error,*999)
    
    EXITS("DecompositionTopology_DataPointsGet")
    RETURN
999 NULLIFY(decompositionDataPoints)
998 ERRORSEXITS("DecompositionTopology_DataPointsGet",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_DataPointsGet

  !
  !================================================================================================================================
  !

  !>Gets faces from a decomposition topology.
  SUBROUTINE DecompositionTopology_FacesGet(decompositionTopology,decompositionFaces,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology !<A pointer to the decomposition topology to get the faces for.
    TYPE(DECOMPOSITION_FACES_TYPE), POINTER :: decompositionFaces !<On exit, a pointer to the decomposition topology faces. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DecompositionTopology_FacesGet",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(decompositionTopology)) CALL FlagError("Decomposition is not associated.",err,error,*999)
    IF(ASSOCIATED(decompositionFaces)) CALL FlagError("Decomposition faces is already associated.",err,error,*999)

    !Get the decomposition faces
    decompositionFaces=>decompositionTopology%faces

    !Check decomposition faces is associated.
    IF(.NOT.ASSOCIATED(decompositionFaces)) CALL FlagError("Decomposition topology faces is not associated.",err,error,*999)
    
    EXITS("DecompositionTopology_FacesGet")
    RETURN
999 ERRORSEXITS("DecompositionTopology_FacesGet",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_FacesGet

  !
  !================================================================================================================================
  !

  !>Gets lines from a decomposition topology.
  SUBROUTINE DecompositionTopology_LinesGet(decompositionTopology,decompositionLines,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology !<A pointer to the decomposition topology to get the lines for.
    TYPE(DECOMPOSITION_LINES_TYPE), POINTER :: decompositionLines !<On exit, a pointer to the decomposition topology lines. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DecompositionTopology_LinesGet",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(decompositionTopology)) CALL FlagError("Decomposition is not associated.",err,error,*999)
    IF(ASSOCIATED(decompositionLines)) CALL FlagError("Decomposition lines is already associated.",err,error,*999)
 
    !Get the decomposition lines
    decompositionLines=>decompositionTopology%lines

    !Check decomposition lines is associated.
    IF(.NOT.ASSOCIATED(decompositionLines)) CALL FlagError("Decomposition topology lines is not associated.",err,error,*999)
    
    EXITS("DecompositionTopology_LinesGet")
    RETURN
999 ERRORSEXITS("DecompositionTopology_LinesGet",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_LinesGet
  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the decomposition identified by user number in the given mesh. If no decomposition with that user number exists decomposition is left nullified.
  SUBROUTINE Decomposition_UserNumberFind(userNumber,mesh,decomposition,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the decomposition to find
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh containing the decomposition to find
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition !<On return a pointer to the decomposition with the specified user number. If no decomposition with that user number exists then decomposition is returned NULL.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: decompositionIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("Decomposition_UserNumberFind",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)
    IF(ASSOCIATED(decomposition)) CALL FlagError("Decomposition is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(mesh%decompositions)) THEN
      localError="The decompositions on mesh number "//TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))// &
        & " are not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
     
    !Get the decomposition from the user number
    NULLIFY(decomposition)
    IF(ALLOCATED(mesh%decompositions%decompositions)) THEN
      DO decompositionIdx=1,mesh%decompositions%NUMBER_OF_DECOMPOSITIONS
        IF(ASSOCIATED(mesh%decompositions%decompositions(decompositionIdx)%ptr)) THEN
          IF(mesh%decompositions%decompositions(decompositionIdx)%ptr%USER_NUMBER==userNumber) THEN
            decomposition=>mesh%decompositions%decompositions(decompositionIdx)%ptr
            EXIT
          ENDIF
        ELSE
          localError="The decomposition pointer in mesh decompositions is not associated for decomposition index "// &
            & TRIM(NumberToVString(decompositionIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !decompositionIdx
    ENDIF
    
    EXITS("Decomposition_UserNumberFind")
    RETURN
999 ERRORSEXITS("Decomposition_UserNumberFind",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_UserNumberFind

  !
  !================================================================================================================================
  !

  !>Get the basis for an element in the domain elements identified by its local number
  SUBROUTINE DomainElements_BasisGet(domainElements,localElementNumber,basis,err,error,*)

    !Argument variables
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: domainElements !<A pointer to the domain elements to get the element basis for
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The element local number to get the basis for
    TYPE(BASIS_TYPE), POINTER, INTENT(OUT) :: basis !<On return, a pointer to the basis for the element. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("DomainElements_BasisGet",err,error,*999)

    IF(.NOT.ASSOCIATED(domainElements)) CALL FlagError("Domain elements is not associated.",err,error,*999)
    IF(localElementNumber<=0.OR.localElementNumber>domainElements%TOTAL_NUMBER_OF_ELEMENTS) THEN
      localError="The local element number of "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " is invalid. The local number must be >= 1 and <= "// &
        & TRIM(NumberToVString(domainElements%TOTAL_NUMBER_OF_ELEMENTS,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(domainElements%elements)) CALL FlagError("Domain elements elements is not associated.",err,error,*999)
      
    basis=>domainElements%elements(localElementNumber)%basis
    IF(.NOT.ASSOCIATED(basis)) THEN
      localError="The basis for local element number "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("DomainElements_BasisGet")
    RETURN
999 ERRORSEXITS("DomainElements_BasisGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainElements_BasisGet

  !
  !================================================================================================================================
  !

  !>Get the face in the domain faces identified by its local number
  SUBROUTINE DomainFaces_FaceGet(domainFaces,faceNumber,face,err,error,*)

    !Argument variables
    TYPE(DOMAIN_FACES_TYPE), POINTER :: domainFaces !<A pointer to the domain faces to get the face for
    INTEGER(INTG), INTENT(IN) :: faceNumber !<The face number to get
    TYPE(DOMAIN_FACE_TYPE), POINTER, INTENT(OUT) :: face !<On return, a pointer to the specified face. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("DomainFaces_FaceGet",err,error,*999)

    IF(ASSOCIATED(face)) CALL FlagError("Face is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainFaces)) CALL FlagError("Domain faces is not associated.",err,error,*999)
    IF(faceNumber<=0.OR.faceNumber>domainFaces%NUMBER_OF_FACES) THEN
      localError="The face number of "//TRIM(NumberToVString(faceNumber,"*",err,error))// &
        & " is invalid. The face number must be >= 1 and <= "// &
        & TRIM(NumberToVString(domainFaces%NUMBER_OF_FACES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainFaces%faces)) CALL FlagError("Domain faces faces is not associated.",err,error,*999)
      
    face=>domainFaces%faces(faceNumber)
    IF(.NOT.ASSOCIATED(face)) THEN
      localError="Face number "//TRIM(NumberToVString(faceNumber,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("DomainFaces_FaceGet")
    RETURN
999 NULLIFY(face)
998 ERRORSEXITS("DomainFaces_FaceGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainFaces_FaceGet    
  
  !
  !================================================================================================================================
  !

  !>Get the line in the domain faces identified by its local number
  SUBROUTINE DomainLines_LineGet(domainLines,lineNumber,line,err,error,*)

    !Argument variables
    TYPE(DOMAIN_LINES_TYPE), POINTER :: domainLines !<A pointer to the domain lines to get the line for
    INTEGER(INTG), INTENT(IN) :: lineNumber !<The line number to get
    TYPE(DOMAIN_LINE_TYPE), POINTER, INTENT(OUT) :: line !<On return, a pointer to the specified line. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("DomainLines_LineGet",err,error,*999)

    IF(ASSOCIATED(line)) CALL FlagError("Line is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainLines)) CALL FlagError("Domain lines is not associated.",err,error,*999)
    IF(lineNumber<=0.OR.lineNumber>domainLines%NUMBER_OF_LINES) THEN
      localError="The line number of "//TRIM(NumberToVString(lineNumber,"*",err,error))// &
        & " is invalid. The line number must be >= 1 and <= "// &
        & TRIM(NumberToVString(domainLines%NUMBER_OF_LINES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainLines%lines)) CALL FlagError("Domain lines lines is not associated.",err,error,*999)
      
    line=>domainLines%lines(lineNumber)
    IF(.NOT.ASSOCIATED(line)) THEN
      localError="Line number "//TRIM(NumberToVString(lineNumber,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("DomainLines_LineGet")
    RETURN
999 NULLIFY(line)
998 ERRORSEXITS("DomainLines_LineGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainLines_LineGet
  
  !
  !================================================================================================================================
  !

  !>Gets a domain mappings from a domain.
  SUBROUTINE Domain_MappingsGet(domain,domainMappings,err,error,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: domain !<The domain to get the mappings for.
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: domainMappings !<On exit, a pointer to the domain mappings. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Domain_MappingsGet",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(domain)) THEN
      CALL FlagError("Domain is not associated.",err,error,*999)
    ENDIF
    IF(ASSOCIATED(domainMappings)) THEN
      CALL FlagError("Domain mappings is already associated.",err,error,*999)
    ENDIF

    !Get the domain mappings
    domainMappings=>domain%mappings

    !Check domain mappings is associated.
    IF(.NOT.ASSOCIATED(domainMappings)) THEN
      CALL FlagError("Domain mappings is not associated.",err,error,*999)
    ENDIF
    
    EXITS("Domain_MappingsGet")
    RETURN
999 ERRORSEXITS("Domain_MappingsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Domain_MappingsGet

  !
  !================================================================================================================================
  !

  !>Gets a domain topology from a domain.
  SUBROUTINE Domain_TopologyGet(domain,domainTopology,err,error,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: domain !<The domain to get the topologya for.
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology !<On exit, a pointer to the domain topology. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Domain_TopologyGet",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(domain)) CALL FlagError("Domain is not associated.",err,error,*999)
    IF(ASSOCIATED(domainTopology)) CALL FlagError("Domain topology is already associated.",err,error,*999)

    !Get the domain topology
    domainTopology=>domain%topology

    !Check domain topology is associated.
    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain topology is not associated.",err,error,*999)
    
    EXITS("Domain_TopologyGet")
    RETURN
999 ERRORSEXITS("Domain_TopologyGet",err,error)
    RETURN 1
    
  END SUBROUTINE Domain_TopologyGet

  !
  !================================================================================================================================
  !

  !>Gets elements from a domain mappings.
  SUBROUTINE DomainMappings_ElementsGet(domainMappings,domainElements,err,error,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: domainMappings !<A pointer to the domain mappings to get the elements for.
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: domainElements !<On exit, a pointer to the domain mapping elements. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DomainMappings_ElementsGet",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(domainMappings)) CALL FlagError("Domain mappings is not associated.",err,error,*999)
    IF(ASSOCIATED(domainElements)) CALL FlagError("Domain mapping elements is already associated.",err,error,*999)

    !Get the domain elements
    domainElements=>domainMappings%elements
    !Check domain elements is associated.
    IF(.NOT.ASSOCIATED(domainElements)) CALL FlagError("Domain mappings elements is not associated.",err,error,*999)
    
    EXITS("DomainMappings_ElementsGet")
    RETURN
999 ERRORSEXITS("DomainMappings_ElementsGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMappings_ElementsGet

  !
  !================================================================================================================================
  !

  !>Gets nodes from a domain mappings.
  SUBROUTINE DomainMappings_NodesGet(domainMappings,domainNodes,err,error,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: domainMappings !<A pointer to the domain mappings to get the nodes for.
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: domainNodes !<On exit, a pointer to the domain mapping nodes. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DomainMappings_NodesGet",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(domainMappings)) CALL FlagError("Domain mappings is not associated.",err,error,*999)
    IF(ASSOCIATED(domainNodes)) CALL FlagError("Domain mapping nodes is already associated.",err,error,*999)

    !Get the domain nodes
    domainNodes=>domainMappings%nodes
    !Check domain nodes is associated.
    IF(.NOT.ASSOCIATED(domainNodes)) CALL FlagError("Domain mappings nodes is not associated.",err,error,*999)
    
    EXITS("DomainMappings_NodesGet")
    RETURN
999 ERRORSEXITS("DomainMappings_NodesGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMappings_NodesGet

  !
  !================================================================================================================================
  !

  !>Gets elements from a domain topology.
  SUBROUTINE DomainTopology_ElementsGet(domainTopology,domainElements,err,error,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology !<A pointer to the domain topology to get the elements for.
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: domainElements !<On exit, a pointer to the domain topology elements. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DomainTopology_ElementsGet",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain is not associated.",err,error,*999)
    IF(ASSOCIATED(domainElements)) CALL FlagError("Domain elements is already associated.",err,error,*999)

    !Get the domain elements
    domainElements=>domainTopology%elements
    !Check domain elements is associated.
    IF(.NOT.ASSOCIATED(domainElements)) CALL FlagError("Domain topology elements is not associated.",err,error,*999)
    
    EXITS("DomainTopology_ElementsGet")
    RETURN
999 ERRORSEXITS("DomainTopology_ElementsGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainTopology_ElementsGet

  !
  !================================================================================================================================
  !

  !>Gets lines from a domain topology.
  SUBROUTINE DomainTopology_LinesGet(domainTopology,domainLines,err,error,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology !<A pointer to the domain topology to get the lines for.
    TYPE(DOMAIN_LINES_TYPE), POINTER :: domainLines !<On exit, a pointer to the domain topology lines. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DomainTopology_LinesGet",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain is not associated.",err,error,*999)
    IF(ASSOCIATED(domainLines)) CALL FlagError("Domain lines is already associated.",err,error,*999)

    !Get the domain lines
    domainLines=>domainTopology%lines
    !Check domain lines is associated.
    IF(.NOT.ASSOCIATED(domainLines)) CALL FlagError("Domain topology lines is not associated.",err,error,*999)
    
    EXITS("DomainTopology_LinesGet")
    RETURN
999 ERRORSEXITS("DomainTopology_LinesGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainTopology_LinesGet

  !
  !================================================================================================================================
  !

  !>Gets nodes from a domain topology.
  SUBROUTINE DomainTopology_NodesGet(domainTopology,domainNodes,err,error,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology !<A pointer to the domain topology to get the nodes for.
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes !<On exit, a pointer to the domain topology nodes. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DomainTopology_NodesGet",err,error,*999)

    !Check input arguments
    IF(ASSOCIATED(domainNodes)) CALL FlagError("Domain lines is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain is not associated.",err,error,*999)

    !Get the domain nodes
    domainNodes=>domainTopology%nodes
    !Check domain nodes is associated.
    IF(.NOT.ASSOCIATED(domainNodes)) CALL FlagError("Domain topology nodes is not associated.",err,error,*999)
    
    EXITS("DomainTopology_NodesGet")
    RETURN
999 NULLIFY(domainNodes)
998 ERRORSEXITS("DomainTopology_NodesGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainTopology_NodesGet

  !
  !================================================================================================================================
  !

  !>Get the basis for an element in the domain identified by its local number
  SUBROUTINE DomainTopology_LocalElementBasisGet(domainTopology,localElementNumber,basis,err,error,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology !<A pointer to the domain topology to get the element basis for
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The element local number to get the basis for
    TYPE(BASIS_TYPE), POINTER, INTENT(OUT) :: basis !<On return, a pointer to the basis for the element.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DomainTopology_LocalElementBasisGet",err,error,*999)

    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain topology is not associated.",err,error,*999)
    IF(localElementNumber<=0.OR.localElementNumber>domainTopology%elements%TOTAL_NUMBER_OF_ELEMENTS) THEN
      localError="The local element number of "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " is invalid. The local number must be >= 1 and <= "// &
        & TRIM(NumberToVString(domainTopology%elements%TOTAL_NUMBER_OF_ELEMENTS,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(domainTopology%elements)) CALL FlagError("Domain topology elements is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(domainTopology%elements%elements)) &
      & CALL FlagError("Domain topology elements elements is not associated.",err,error,*999)
       
    basis=>domainTopology%elements%elements(localElementNumber)%basis
    IF(.NOT.ASSOCIATED(basis)) THEN
      localError="The basis for local element number "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("DomainTopology_LocalElementBasisGet")
    RETURN
999 ERRORSEXITS("DomainTopology_LocalElementBasisGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainTopology_LocalElementBasisGet

  !
  !================================================================================================================================
  !

  !>Gets faces from a domain topology.
  SUBROUTINE DomainTopology_FacesGet(domainTopology,domainFaces,err,error,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology !<A pointer to the domain topology to get the faces for.
    TYPE(DOMAIN_FACES_TYPE), POINTER :: domainFaces !<On exit, a pointer to the domain topology faces. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DomainTopology_FacesGet",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain is not associated.",err,error,*999)
    IF(ASSOCIATED(domainFaces)) CALL FlagError("Domain faces is already associated.",err,error,*999)

    !Get the domain faces
    domainFaces=>domainTopology%faces
    !Check domain faces is associated.
    IF(.NOT.ASSOCIATED(domainFaces)) CALL FlagError("Domain topology faces is not associated.",err,error,*999)
    
    EXITS("DomainTopology_FacesGet")
    RETURN
999 ERRORSEXITS("DomainTopology_FacesGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainTopology_FacesGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the decomposition for a given user number in a mesh. \see OPENCMISS::Iron::cmfe_Mesh_DecompositionGet
  SUBROUTINE Mesh_DecompositionGet(mesh,userNumber,decomposition,err,error,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh to get the decomposition for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the decomposition to get.
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition !<On exit, a pointer to the decomposition for the mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Mesh_DecompositionGet",err,error,*998)

    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*998)
    IF(.NOT.mesh%MESH_FINISHED) CALL FlagError("Mesh has not been finished.",err,error,*998)
    IF(ASSOCIATED(decomposition)) CALL FlagError("Decomposition is already associated.",err,error,*998)
    
    NULLIFY(decomposition)
    CALL Decomposition_UserNumberFind(userNumber,mesh,decomposition,err,error,*999)
    IF(.NOT.ASSOCIATED(decomposition)) THEN
      localError="A decomposition with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " do not exist on mesh number "//TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Mesh_DecompositionGet")
    RETURN
999 NULLIFY(decomposition)
998 ERRORSEXITS("Mesh_DecompositionGet",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_DecompositionGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the mesh elements for a given mesh component.
  SUBROUTINE Mesh_MeshElementsGet(mesh,meshComponentNumber,meshElements,err,error,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh to get the elements for
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component number to get the elements for
    TYPE(MeshElementsType), POINTER :: meshElements !<On return, a pointer to the mesh elements
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Mesh_MeshElementsGet",err,error,*998)
    
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated",err,error,*998)
    IF(ASSOCIATED(meshElements)) CALL FlagError("Elements is already associated.",err,error,*998)
    IF(meshComponentNumber<=0.OR.meshComponentNumber>mesh%NUMBER_OF_COMPONENTS) THEN
      localError="The specified mesh component number of "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
        & " is invalid for mesh number "//TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))// &
        & ". The component number must be between 1 and "//TRIM(NumberToVString(mesh%NUMBER_OF_COMPONENTS,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ASSOCIATED(mesh%topology(meshComponentNumber)%ptr)) THEN
      localError="The mesh topology is not associated for mesh component number "// &
        & TRIM(NumberToVString(meshComponentNumber,"*",err,error))//" of mesh number "// &
        & TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    meshElements=>mesh%topology(meshComponentNumber)%ptr%elements
    IF(.NOT.ASSOCIATED(meshElements)) THEN
      localError="The mesh topology elements is not associated for mesh component number "// &
        & TRIM(NumberToVString(meshComponentNumber,"*",err,error))//" of mesh number "// &
        & TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("Mesh_MeshElementsGet")
    RETURN
999 NULLIFY(meshElements)
998 ERRORSEXITS("Mesh_MeshElementsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_MeshElementsGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the mesh nodes for a given mesh component.
  SUBROUTINE Mesh_MeshNodesGet(mesh,meshComponentNumber,meshNodes,err,error,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh to get the nodes for
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component number to get the nodes for
    TYPE(MeshNodesType), POINTER :: meshNodes !<On return, a pointer to the mesh nodes. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Mesh_MeshNodesGet",err,error,*998)
    
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated",err,error,*998)
    IF(ASSOCIATED(meshNodes)) CALL FlagError("Nodes is already associated.",err,error,*998)    
    IF(meshComponentNumber<=0.OR.meshComponentNumber>mesh%NUMBER_OF_COMPONENTS) THEN
      localError="The specified mesh component number of "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
        & " is invalid for mesh number "//TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))// &
        & ". The component number must be between 1 and "//TRIM(NumberToVString(mesh%NUMBER_OF_COMPONENTS,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ASSOCIATED(mesh%topology(meshComponentNumber)%ptr)) THEN
      localError="The mesh topology is not associated for mesh component number "// &
        & TRIM(NumberToVString(meshComponentNumber,"*",err,error))//" of mesh number "// &
        & TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    meshNodes=>mesh%topology(meshComponentNumber)%ptr%nodes
    IF(.NOT.ASSOCIATED(meshNodes)) THEN
      localError="The mesh topology nodes is not associated for mesh component number "// &
        & TRIM(NumberToVString(meshComponentNumber,"*",err,error))//" of mesh number "// &
        & TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("Mesh_MeshNodesGet")
    RETURN
999 NULLIFY(meshNodes)
998 ERRORSEXITS("Mesh_MeshNodesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_MeshNodesGet

  !
  !================================================================================================================================
  !

  !>Returns the region for a mesh accounting for regions and interfaces
  SUBROUTINE Mesh_RegionGet(mesh,region,err,error,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh to get the region for
    TYPE(REGION_TYPE), POINTER :: region !<On return, the meshes region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(INTERFACE_TYPE), POINTER :: interface
    TYPE(REGION_TYPE), POINTER :: parentRegion
    TYPE(VARYING_STRING) :: localError

    ENTERS("Mesh_RegionGet",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)
    IF(ASSOCIATED(region)) CALL FlagError("Region is already associated.",err,error,*999)
    
    NULLIFY(region)
    NULLIFY(interface)
    region=>mesh%region
    IF(.NOT.ASSOCIATED(region)) THEN
      INTERFACE=>mesh%INTERFACE
      IF(ASSOCIATED(INTERFACE)) THEN
        parentRegion=>interface%PARENT_REGION
        IF(ASSOCIATED(parentRegion)) THEN
          region=>parentRegion
        ELSE
          localError="The parent region is not associated for mesh number "// &
            & TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))//" of interface number "// &
            & TRIM(NumberToVString(interface%USER_NUMBER,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        localError="The region or interface is not associated for mesh number "// &
          & TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF

    EXITS("Mesh_RegionGet")
    RETURN
999 ERRORSEXITS("Mesh_RegionGet",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_RegionGet

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the mesh identified by user number in the given list of meshes. If no mesh with that number exits mesh is left nullified.
  SUBROUTINE Mesh_UserNumberFindGeneric(userNumber,meshes,mesh,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the mesh to find
    TYPE(MESHES_TYPE), POINTER :: meshes !<The list of meshes containing the mesh.
    TYPE(MESH_TYPE), POINTER :: mesh !<On return, a pointer to the mesh of the specified user number. In no mesh with the specified user number exists the pointer is returned NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("Mesh_UserNumberFindGeneric",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(meshes)) CALL FlagError("Meshes is not associated",err,error,*999)
    IF(ASSOCIATED(mesh)) CALL FlagError("Mesh is already associated.",err,error,*999)
    
    !Get the mesh from the user number
    NULLIFY(mesh)
    IF(ASSOCIATED(meshes%meshes)) THEN
      DO meshIdx=1,meshes%NUMBER_OF_MESHES
        IF(ASSOCIATED(meshes%meshes(meshIdx)%ptr)) THEN
          IF(meshes%meshes(meshIdx)%ptr%USER_NUMBER==userNumber) THEN
            mesh=>meshes%meshes(meshIdx)%ptr
            EXIT
          ENDIF
        ELSE
          localError="The mesh pointer in meshes is not associated for mesh index "// &
            & TRIM(NumberToVString(meshIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)          
        ENDIF
      ENDDO !meshIdx      
    ENDIF
    
    EXITS("Mesh_UserNumberFindGeneric")
    RETURN
999 ERRORSEXITS("Mesh_UserNumberFindGeneric",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_UserNumberFindGeneric

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the mesh identified by user number in the given interface. If no mesh with that number exits mesh is left nullified.
  SUBROUTINE Mesh_UserNumberFindInterface(userNumber,interface,mesh,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the mesh to find
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface containing the mesh
    TYPE(MESH_TYPE), POINTER :: mesh !<On return, a pointer to the mesh of the specified user number. In no mesh with the specified user number exists the pointer is returned NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
  
    ENTERS("Mesh_UserNumberFindInterface",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated",err,error,*999)
    
    CALL Mesh_UserNumberFindGeneric(userNumber,interface%meshes,mesh,err,error,*999)
     
    EXITS("Mesh_UserNumberFindInterface")
    RETURN
999 ERRORSEXITS("Mesh_UserNumberFindInterface",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_UserNumberFindInterface

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the mesh identified by user number in the given. If no mesh with that number exits mesh is left nullified.
  SUBROUTINE Mesh_UserNumberFindRegion(userNumber,region,mesh,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the mesh to find
    TYPE(REGION_TYPE), POINTER :: REGION !<The region containing the mesh
    TYPE(MESH_TYPE), POINTER :: MESH !<On return, a pointer to the mesh of the specified user number. In no mesh with the specified user number exists the pointer is returned NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("Mesh_UserNumberFindRegion",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated",err,error,*999)
    
    CALL Mesh_UserNumberFindGeneric(userNumber,region%meshes,mesh,err,error,*999)
     
    EXITS("Mesh_UserNumberFindRegion")
    RETURN
999 ERRORSEXITS("Mesh_UserNumberFindRegion",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_UserNumberFindRegion

  !
  !================================================================================================================================
  !

  !>Returns the user number for a mesh.
  SUBROUTINE Mesh_UserNumberGet(mesh,userNumber,err,error,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh to get the user number for
    INTEGER(INTG), INTENT(OUT) :: userNumber !<On exit, the user number of the mesh.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Mesh_UserNumberGet",err,error,*999)

    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)

    userNumber=mesh%USER_NUMBER
  
    EXITS("Mesh_UserNumberGet")
    RETURN
999 ERRORSEXITS("Mesh_UserNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_UserNumberGet

  !
  !================================================================================================================================
  !
  
END MODULE MeshAccessRoutines
