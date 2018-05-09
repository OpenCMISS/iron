!> \file
!> \author Chris Bradley
!> \brief This module contains all decomposition access method routines.
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
!> Contributor(s): Chris Bradley
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

!> This module contains all decomposition access method routines.
MODULE DecompositionAccessRoutines
  
  USE BaseRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE Strings
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters
   
  !> \addtogroup DecompositionRoutines_DecompositionTypes DecompositionRoutines::DecompositionTypes
  !> \brief The Decomposition types parameters
  !> \see DecompositionRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: DECOMPOSITION_ALL_TYPE=1 !<The decomposition contains all elements. \see DecompositionRoutines_DecompositionTypes,DecompositionRoutines
  INTEGER(INTG), PARAMETER :: DECOMPOSITION_CALCULATED_TYPE=2 !<The element decomposition is calculated by graph partitioning. \see DecompositionRoutines_DecompositionTypes,DecompositionRoutines
  INTEGER(INTG), PARAMETER :: DECOMPOSITION_USER_DEFINED_TYPE=3 !<The user will set the element decomposition. \see DecompositionRoutines_DecompositionTypes,DecompositionRoutines
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

  INTERFACE DECOMPOSITION_USER_NUMBER_FIND
    MODULE PROCEDURE Decomposition_UserNumberFind
  END INTERFACE DECOMPOSITION_USER_NUMBER_FIND

  PUBLIC DECOMPOSITION_ALL_TYPE,DECOMPOSITION_CALCULATED_TYPE,DECOMPOSITION_USER_DEFINED_TYPE

  PUBLIC Decomposer_AssertIsFinished,Decomposer_AssertNotFinished

  PUBLIC Decomposer_RegionGet

  PUBLIC Decomposer_UserNumberFind
  
  PUBLIC Decomposition_AssertIsFinished,Decomposition_AssertNotFinished

  PUBLIC Decomposition_CoordinateSystemGet

  PUBLIC Decomposition_DomainGet

  PUBLIC Decomposition_MeshGet

  PUBLIC Decomposition_RegionGet
  
  PUBLIC Decomposition_TopologyGet
  
  PUBLIC Decomposition_UserNumberFind

  PUBLIC DECOMPOSITION_USER_NUMBER_FIND

  PUBLIC Decomposition_WorkGroupGet

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
  
  PUBLIC DomainTopology_LocalElementBasisGet

  PUBLIC DomainTopology_NodesGet
  
  
CONTAINS

  !
  !=================================================================================================================================
  !

  !>Assert that a decomposer has been finished
  SUBROUTINE Decomposer_AssertIsFinished(decomposer,err,error,*)

    !Argument Variables
    TYPE(DecomposerType), POINTER, INTENT(INOUT) :: decomposer !<The decomposer to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Decomposer_AssertIsFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(decomposer)) CALL FlagError("Decomposer is not associated.",err,error,*999)

    IF(.NOT.decomposer%decomposerFinished) THEN
      localError="Decomposer number "//TRIM(NumberToVString(decomposer%userNumber,"*",err,error))// &
        & " has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Decomposer_AssertIsFinished")
    RETURN
999 ERRORSEXITS("Decomposer_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a decomposer has not been finished
  SUBROUTINE Decomposer_AssertNotFinished(decomposer,err,error,*)

    !Argument Variables
    TYPE(DecomposerType), POINTER, INTENT(INOUT) :: decomposer !<The decomposer to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Decomposer_AssertNotFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(decomposer)) CALL FlagError("Decomposer is not associated.",err,error,*999)

    IF(decomposer%decomposerFinished) THEN
      localError="Decomposer number "//TRIM(NumberToVString(decomposer%userNumber,"*",err,error))// &
        & " has already been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Decomposer_AssertNotFinished")
    RETURN
999 ERRORSEXITS("Decomposer_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Returns the region for a decomposer
  SUBROUTINE Decomposer_RegionGet(decomposer,region,err,error,*)

    !Argument variables
    TYPE(DecomposerType), POINTER :: decomposer !<A pointer to the decomposer to get the region for
    TYPE(RegionType), POINTER :: region !<On return, the decomposer region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Decomposer_RegionGet",err,error,*999)

    !Check input arguments
    IF(ASSOCIATED(region)) CALL FlagError("Region is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(decomposer)) CALL FlagError("Decomposer is not associated.",err,error,*999)

    region=>decomposer%region
    IF(.NOT.ASSOCIATED(region)) THEN
      localError="The region associated with decomposer number "// &
        & TRIM(NumberToVString(decomposer%userNumber,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
 
    EXITS("Decomposer_RegionGet")
    RETURN
999 ERRORSEXITS("Decomposer_RegionGet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_RegionGet

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the decomposer identified by user number in the given region. If no decomposer with that user number exists the decomposer is left nullified.
  SUBROUTINE Decomposer_UserNumberFind(userNumber,region,decomposer,err,error,*)

    !Argument variables 
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the equation set to find.
    TYPE(RegionType), POINTER :: region !<The region to find the equations set in
    TYPE(DecomposerType), POINTER :: decomposer !<On return, a pointer to the decomposer if an equations set with the specified user number exists in the given region. If no decompser with the specified number exists a NULL pointer is returned. The pointer must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposerIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("Decomposer_UserNumberFind",err,error,*999)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)
    IF(ASSOCIATED(decomposer)) CALL FlagError("Equations set is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(region%decomposers)) THEN
      localError="The decomposers on region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    NULLIFY(decomposer)
    IF(ALLOCATED(region%decomposers%decomposers)) THEN
      DO decomposerIdx=1,region%decomposers%numberOfDecomposers
        IF(ASSOCIATED(region%decomposers%decomposers(decomposerIdx)%ptr)) THEN
          IF(region%decomposers%decomposers(decomposerIdx)%ptr%userNumber==userNumber) THEN
            decomposer=>region%decomposers%decomposers(decomposerIdx)%ptr
            EXIT
          ENDIF
        ELSE
          localError="The decomposer pointer in region decomposers is not associated for decomposer index "// &
            & TRIM(NumberToVString(decomposerIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !decomposerIdx
    ENDIF
    
    EXITS("Decomposer_UserNumberFind")
    RETURN
999 ERRORSEXITS("Decomposer_UserNumberFind",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_UserNumberFind

  !
  !=================================================================================================================================
  !

  !>Assert that a decomposition has been finished
  SUBROUTINE Decomposition_AssertIsFinished(decomposition,err,error,*)

    !Argument Variables
    TYPE(DecompositionType), POINTER, INTENT(INOUT) :: decomposition !<The decomposition to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Decomposition_AssertIsFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)

    IF(.NOT.decomposition%decompositionFinished) THEN
      localError="Decomposition number "//TRIM(NumberToVString(decomposition%userNumber,"*",err,error))// &
        & " has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Decomposition_AssertIsFinished")
    RETURN
999 ERRORSEXITS("Decomposition_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a decomposition has not been finished
  SUBROUTINE Decomposition_AssertNotFinished(decomposition,err,error,*)

    !Argument Variables
    TYPE(DecompositionType), POINTER, INTENT(INOUT) :: decomposition !<The decomposition to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Decomposition_AssertNotFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)

    IF(decomposition%decompositionFinished) THEN
      localError="Decomposition number "//TRIM(NumberToVString(decomposition%userNumber,"*",err,error))// &
        & " has already been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Decomposition_AssertNotFinished")
    RETURN
999 ERRORSEXITS("Decomposition_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Returns the coordinate system for a decomposition accounting for regions and interfaces. 
  SUBROUTINE Decomposition_CoordinateSystemGet(decomposition,coordinateSystem,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<The decomposition to get the coordinate system for.
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<On exit, a pointer to the coordinate system of decomposition. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Decomposition_CoordinateSystemGet",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
    IF(ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is already associated.",err,error,*999)

    IF(ASSOCIATED(decomposition%region)) THEN
      coordinateSystem=>decomposition%region%coordinateSystem
    ELSE IF(ASSOCIATED(decomposition%INTERFACE)) THEN
      coordinateSystem=>decomposition%interface%coordinateSystem
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
    TYPE(DecompositionType), POINTER :: decomposition !<The decomposition to get the domain for.
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component number of the domain to get. If mesh component number is 0 the mesh component used to construct the decomposition is used.
    TYPE(DomainType), POINTER :: domain !<On exit, a pointer to the domain of decomposition of the specified mesh component number. Must not be associated on entry.
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
      meshComponent=decomposition%meshComponentNumber
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
  
  !>Returns a pointer to the mesh for a decomposition. 
  SUBROUTINE Decomposition_MeshGet(decomposition,mesh,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to get the mesh for.
    TYPE(MeshType), POINTER :: mesh !<On exit, A pointer to the mesh for the decomposition. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Decomposition_MeshGet",err,error,*998)

    IF(ASSOCIATED(mesh)) CALL FlagError("Mesh is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
      
    mesh=>decomposition%mesh
    IF(.NOT.ASSOCIATED(mesh)) THEN
      localError="The mesh associated with decomposition number "// &
        & TRIM(NumberToVString(decomposition%userNumber,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("Decomposition_MeshGet")
    RETURN
999 NULLIFY(mesh)
998 ERRORSEXITS("Decomposition_MeshGet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_MeshGet

  !
  !================================================================================================================================
  !

  !>Returns the region for a decomposition
  SUBROUTINE Decomposition_RegionGet(decomposition,region,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to get the region for
    TYPE(RegionType), POINTER :: region !<On return, the decomposition region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Decomposition_RegionGet",err,error,*999)

    !Check input arguments
    IF(ASSOCIATED(region)) CALL FlagError("Region is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)

    IF(ASSOCIATED(decomposition%region)) THEN
      region=>decomposition%region
    ELSE IF(ASSOCIATED(decomposition%interface)) THEN
      region=>decomposition%interface%parentRegion
    ELSE
      localError="Decomposition number "//TRIM(NumberToVString(decomposition%userNumber,"*",err,error))// &
        & " is not associated with a region or an interface."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ASSOCIATED(region)) THEN
      localError="The region associated with decomposition number "// &
        & TRIM(NumberToVString(decomposition%userNumber,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
 
    EXITS("Decomposition_RegionGet")
    RETURN
999 ERRORSEXITS("Decomposition_RegionGet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_RegionGet

  !
  !================================================================================================================================
  !

  !>Gets a decomposition topology from a decomposition.
  SUBROUTINE Decomposition_TopologyGet(decomposition,decompositionTopology,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<The decomposition to get the domain for.
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<On exit, a pointer to the decomposition topology. Must not be associated on entry.
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
        & TRIM(NumberToVString(decomposition%userNumber,"*",err,error))//"."
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

  !>Finds and returns a pointer to the decomposition identified by user number in the given mesh. If no decomposition with that user number exists decomposition is left nullified.
  SUBROUTINE Decomposition_UserNumberFind(userNumber,mesh,decomposition,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the decomposition to find
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh containing the decomposition to find
    TYPE(DecompositionType), POINTER :: decomposition !<On return a pointer to the decomposition with the specified user number. If no decomposition with that user number exists then decomposition is returned NULL.
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
      localError="The decompositions on mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))// &
        & " are not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
     
    !Get the decomposition from the user number
    NULLIFY(decomposition)
    IF(ALLOCATED(mesh%decompositions%decompositions)) THEN
      DO decompositionIdx=1,mesh%decompositions%numberOfDecompositions
        IF(ASSOCIATED(mesh%decompositions%decompositions(decompositionIdx)%ptr)) THEN
          IF(mesh%decompositions%decompositions(decompositionIdx)%ptr%userNumber==userNumber) THEN
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
  
  !>Returns a pointer to the work group for a decomposition. 
  SUBROUTINE Decomposition_WorkGroupGet(decomposition,workGroup,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to get the work group for.
    TYPE(WorkGroupType), POINTER :: workGroup !<On exit, A pointer to the work group for the decomposition. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Decomposition_WorkGroupGet",err,error,*998)

    IF(ASSOCIATED(workGroup)) CALL FlagError("Work group is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
      
    workGroup=>decomposition%workGroup
    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("The decomposition work group is not associated.",err,error,*999)
       
    EXITS("Decomposition_WorkGroupGet")
    RETURN
999 NULLIFY(workGroup)
998 ERRORSEXITS("Decomposition_WorkGroupGet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_WorkGroupGet

  !
  !================================================================================================================================
  !

  !>Gets elements from a decomposition topology.
  SUBROUTINE DecompositionTopology_ElementsGet(decompositionTopology,decompositionElements,err,error,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to get the elements for.
     TYPE(DecompositionElementsType), POINTER :: decompositionElements !<On exit, a pointer to the decomposition topology elements. Must not be associated on entry.
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
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to get the data points for.
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
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to get the faces for.
    TYPE(DecompositionFacesType), POINTER :: decompositionFaces !<On exit, a pointer to the decomposition topology faces. Must not be associated on entry.
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
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to get the lines for.
    TYPE(DecompositionLinesType), POINTER :: decompositionLines !<On exit, a pointer to the decomposition topology lines. Must not be associated on entry.
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

  !>Get the basis for an element in the domain elements identified by its local number
  SUBROUTINE DomainElements_BasisGet(domainElements,localElementNumber,basis,err,error,*)

    !Argument variables
    TYPE(DomainElementsType), POINTER :: domainElements !<A pointer to the domain elements to get the element basis for
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The element local number to get the basis for
    TYPE(BasisType), POINTER, INTENT(OUT) :: basis !<On return, a pointer to the basis for the element. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("DomainElements_BasisGet",err,error,*999)

    IF(.NOT.ASSOCIATED(domainElements)) CALL FlagError("Domain elements is not associated.",err,error,*999)
    IF(localElementNumber<=0.OR.localElementNumber>domainElements%totalNumberOfElements) THEN
      localError="The local element number of "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " is invalid. The local number must be >= 1 and <= "// &
        & TRIM(NumberToVString(domainElements%totalNumberOfElements,"*",err,error))//"."
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
    TYPE(DomainFacesType), POINTER :: domainFaces !<A pointer to the domain faces to get the face for
    INTEGER(INTG), INTENT(IN) :: faceNumber !<The face number to get
    TYPE(DomainFaceType), POINTER, INTENT(OUT) :: face !<On return, a pointer to the specified face. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("DomainFaces_FaceGet",err,error,*999)

    IF(ASSOCIATED(face)) CALL FlagError("Face is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainFaces)) CALL FlagError("Domain faces is not associated.",err,error,*999)
    IF(faceNumber<=0.OR.faceNumber>domainFaces%numberOfFaces) THEN
      localError="The face number of "//TRIM(NumberToVString(faceNumber,"*",err,error))// &
        & " is invalid. The face number must be >= 1 and <= "// &
        & TRIM(NumberToVString(domainFaces%numberOfFaces,"*",err,error))//"."
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
    TYPE(DomainLinesType), POINTER :: domainLines !<A pointer to the domain lines to get the line for
    INTEGER(INTG), INTENT(IN) :: lineNumber !<The line number to get
    TYPE(DomainLineType), POINTER, INTENT(OUT) :: line !<On return, a pointer to the specified line. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("DomainLines_LineGet",err,error,*999)

    IF(ASSOCIATED(line)) CALL FlagError("Line is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainLines)) CALL FlagError("Domain lines is not associated.",err,error,*999)
    IF(lineNumber<=0.OR.lineNumber>domainLines%numberOfLines) THEN
      localError="The line number of "//TRIM(NumberToVString(lineNumber,"*",err,error))// &
        & " is invalid. The line number must be >= 1 and <= "// &
        & TRIM(NumberToVString(domainLines%numberOfLines,"*",err,error))//"."
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
    TYPE(DomainType), POINTER :: domain !<The domain to get the mappings for.
    TYPE(DomainMappingsType), POINTER :: domainMappings !<On exit, a pointer to the domain mappings. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Domain_MappingsGet",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(domain)) CALL FlagError("Domain is not associated.",err,error,*999)
    IF(ASSOCIATED(domainMappings)) CALL FlagError("Domain mappings is already associated.",err,error,*999)

    !Get the domain mappings
    domainMappings=>domain%mappings

    !Check domain mappings is associated.
    IF(.NOT.ASSOCIATED(domainMappings)) CALL FlagError("Domain mappings is not associated.",err,error,*999)
    
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
    TYPE(DomainType), POINTER :: domain !<The domain to get the topologya for.
    TYPE(DomainTopologyType), POINTER :: domainTopology !<On exit, a pointer to the domain topology. Must not be associated on entry.
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
    TYPE(DomainMappingsType), POINTER :: domainMappings !<A pointer to the domain mappings to get the elements for.
    TYPE(DomainMappingType), POINTER :: domainElements !<On exit, a pointer to the domain mapping elements. Must not be associated on entry.
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
    TYPE(DomainMappingsType), POINTER :: domainMappings !<A pointer to the domain mappings to get the nodes for.
    TYPE(DomainMappingType), POINTER :: domainNodes !<On exit, a pointer to the domain mapping nodes. Must not be associated on entry.
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
    TYPE(DomainTopologyType), POINTER :: domainTopology !<A pointer to the domain topology to get the elements for.
    TYPE(DomainElementsType), POINTER :: domainElements !<On exit, a pointer to the domain topology elements. Must not be associated on entry.
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
    TYPE(DomainTopologyType), POINTER :: domainTopology !<A pointer to the domain topology to get the lines for.
    TYPE(DomainLinesType), POINTER :: domainLines !<On exit, a pointer to the domain topology lines. Must not be associated on entry.
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
    TYPE(DomainTopologyType), POINTER :: domainTopology !<A pointer to the domain topology to get the nodes for.
    TYPE(DomainNodesType), POINTER :: domainNodes !<On exit, a pointer to the domain topology nodes. Must not be associated on entry.
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
    TYPE(DomainTopologyType), POINTER :: domainTopology !<A pointer to the domain topology to get the element basis for
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The element local number to get the basis for
    TYPE(BasisType), POINTER, INTENT(OUT) :: basis !<On return, a pointer to the basis for the element.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DomainTopology_LocalElementBasisGet",err,error,*999)

    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain topology is not associated.",err,error,*999)
    IF(localElementNumber<=0.OR.localElementNumber>domainTopology%elements%totalNumberOfElements) THEN
      localError="The local element number of "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " is invalid. The local number must be >= 1 and <= "// &
        & TRIM(NumberToVString(domainTopology%elements%totalNumberOfElements,"*",err,error))//"."
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
    TYPE(DomainTopologyType), POINTER :: domainTopology !<A pointer to the domain topology to get the faces for.
    TYPE(DomainFacesType), POINTER :: domainFaces !<On exit, a pointer to the domain topology faces. Must not be associated on entry.
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

END MODULE DecompositionAccessRoutines
