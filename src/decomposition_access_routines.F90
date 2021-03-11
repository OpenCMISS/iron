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
  USE Trees
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
  
  !> \addtogroup DecompositionRoutines_DecomposerOutputTypes DecompositionRoutines::DecomposerOutputTypes
  !> \brief The Decomposer output type parameters
  !> \see DecompositionRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: DECOMPOSER_NO_OUTPUT=0 !<No decomposer output. \see DecompositionRoutines_DecomposerOutputTypes,DecompositionRoutines
  INTEGER(INTG), PARAMETER :: DECOMPOSER_TIMING_OUTPUT=1 !<Timing decomposer output. \see DecompositionRoutines_DecomposerOutputTypes,DecompositionRoutines
  INTEGER(INTG), PARAMETER :: DECOMPOSER_ALL_OUTPUT=2 !<All decomposer output. \see DecompositionRoutines_DecomposerOutputTypes,DecompositionRoutines
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

  PUBLIC DECOMPOSITION_ALL_TYPE,DECOMPOSITION_CALCULATED_TYPE,DECOMPOSITION_USER_DEFINED_TYPE

  PUBLIC DECOMPOSER_NO_OUTPUT,DECOMPOSER_TIMING_OUTPUT,DECOMPOSER_ALL_OUTPUT

  PUBLIC Decomposer_AssertIsFinished,Decomposer_AssertNotFinished

  PUBLIC Decomposer_DecomposerGraphGet

  PUBLIC Decomposer_DecompositionGet

  PUBLIC Decomposer_RegionGet

  PUBLIC Decomposer_UserNumberFind

  PUBLIC DecomposerGraph_DecomposerGet

  PUBLIC DecomposerGraph_RootNodeGet

  PUBLIC DecomposerGraphLink_DecompositionGet

  PUBLIC DecomposerGraphLink_LinkedNodeGet

  PUBLIC DecomposerGraphNode_DecomposerGraphGet

  PUBLIC DecomposerGraphNode_DecompositionGet

  PUBLIC DecomposerGraphNode_DecomposerGraphLinkGet
  
  PUBLIC Decomposition_AssertIsDecomposed,Decomposition_AssertNotDecomposed

  PUBLIC Decomposition_AssertIsFinished,Decomposition_AssertNotFinished

  PUBLIC Decomposition_AssertCalculateFaces,Decomposition_AssertNotCalculateFaces

  PUBLIC Decomposition_AssertCalculateLines,Decomposition_AssertNotCalculateLines

  PUBLIC Decomposition_CalculateFacesGet

  PUBLIC Decomposition_CalculateLinesGet

  PUBLIC Decomposition_CoordinateSystemGet

  PUBLIC Decomposition_DecompositionsGet

  PUBLIC Decomposition_DomainGet

  PUBLIC Decomposition_InterfaceGet

  PUBLIC Decomposition_IsInterfaceDecomposition

  PUBLIC Decomposition_IsRegionDecomposition

  PUBLIC Decomposition_MeshGet

  PUBLIC Decomposition_RegionGet
  
  PUBLIC Decomposition_DecompositionTopologyGet
  
  PUBLIC Decomposition_UserNumberFind

  PUBLIC Decomposition_WorkGroupGet

  PUBLIC DecompositionDataPoints_DataPointCheckExists

  PUBLIC DecompositionDataPoints_ElementDataGlobalNumberGet
  
  PUBLIC DecompositionDataPoints_ElementDataLocalNumberGet

  PUBLIC DecompositionDataPoints_ElementDataNumbersGet

  PUBLIC DecompositionDataPoints_ElementDataUserNumberGet

  PUBLIC DecompositionDataPoints_ElementNumberOfDataPointsGet

  PUBLIC DecompositionDataPoints_LocalDataPointNumberGet

  PUBLIC DecompositionElements_ElementAdjacentNumberGet

  PUBLIC DecompositionElements_ElementBoundaryElementGet

  PUBLIC DecompositionElements_ElementCheckExists

  PUBLIC DecompositionElements_ElementFaceNumberGet

  PUBLIC DecompositionElements_ElementGet

  PUBLIC DecompositionElements_ElementLineNumberGet

  PUBLIC DecompositionElements_ElementNumberAdjacentGet

  PUBLIC DecompositionElements_GlobalElementNumberGet

  PUBLIC DecompositionElements_LocalElementNumberGet

  PUBLIC DecompositionElements_NumberOfElementsGet

  PUBLIC DecompositionElements_TotalNumberOfElementsGet

  PUBLIC DecompositionElements_UserElementNumberGet

  PUBLIC DecompositionFaces_FaceGet

  PUBLIC DecompositionFaces_FaceBoundaryFaceGet

  PUBLIC DecompositionFaces_FaceXiNormalDirectionGet

  PUBLIC DecompositionFaces_NumberOfFacesGet

  PUBLIC DecompositionLines_LineGet

  PUBLIC DecompositionLines_LineBoundaryLineGet

  PUBLIC DecompositionLines_LineXiDirectionGet

  PUBLIC DecompositionLines_NumberOfLinesGet

  PUBLIC DecompositionTopology_DecompositionGet

  PUBLIC DecompositionTopology_DecompositionDataPointsGet

  PUBLIC DecompositionTopology_DecompositionElementsGet

  PUBLIC DecompositionTopology_DecompositionFacesGet
  
  PUBLIC DecompositionTopology_DecompositionLinesGet

  PUBLIC Domain_DecompositionGet

  PUBLIC Domain_DomainMappingsGet

  PUBLIC Domain_DomainTopologyGet

  PUBLIC Domain_MeshComponentNumberGet

  PUBLIC Domain_NumberOfDimensionsGet

  PUBLIC Domain_RegionGet

  PUBLIC DomainElement_BasisGet

  PUBLIC DomainElement_ElementDerivativeGet
  
  PUBLIC DomainElement_ElementNodeGet

  PUBLIC DomainElement_ElementVersionGet
  
  PUBLIC DomainElements_ElementBasisGet

  PUBLIC DomainElements_DomainElementGet

  PUBLIC DomainElements_ElementDerivativeGet

  PUBLIC DomainElements_ElementNodeGet

  PUBLIC DomainElements_ElementVersionGet

  PUBLIC DomainElements_MaxElementParametersGet

  PUBLIC DomainElements_NumberOfElementsGet
  
  PUBLIC DomainElements_TotalNumberOfElementsGet
  
  PUBLIC DomainFaces_DerivativeGlobalIndexGet

  PUBLIC DomainFaces_DerivativeVersionNumberGet

  PUBLIC DomainFaces_FaceBasisGet

  PUBLIC DomainFaces_FaceBoundaryFaceGet

  PUBLIC DomainFaces_FaceGet

  PUBLIC DomainFaces_FaceNodeNumberGet

  PUBLIC DomainLines_DerivativeGlobalIndexGet

  PUBLIC DomainLines_DerivativeVersionNumberGet

  PUBLIC DomainLines_LineBasisGet

  PUBLIC DomainLines_LineBoundaryLineGet

  PUBLIC DomainLines_LineGet

  PUBLIC DomainLines_LineNodeNumberGet

  PUBLIC DomainMappings_DofsMappingGet
  
  PUBLIC DomainMappings_DomainGet
  
  PUBLIC DomainMappings_ElementsMappingGet

  PUBLIC DomainMappings_NodesMappingGet

  PUBLIC DomainNode_BoundaryNodeGet

  PUBLIC DomainNode_NodeDerivativeGet

  PUBLIC DomainNode_GlobalNodeNumberGet

  PUBLIC DomainNode_LocalNodeNumberGet

  PUBLIC DomainNode_MeshNodeNumberGet

  PUBLIC DomainNode_NodeFaceGet

  PUBLIC DomainNode_NodeLineGet

  PUBLIC DomainNode_NumberOfDerivativesGet

  PUBLIC DomainNode_NumberOfNodeFacesGet

  PUBLIC DomainNode_NumberOfNodeLinesGet

  PUBLIC DomainNode_NumberOfSurroundingElementsGet

  PUBLIC DomainNode_SurroundingElementGet

  PUBLIC DomainNode_UserNodeNumberGet

  PUBLIC DomainNodeDerivative_DerivativeGlobalIndexGet

  PUBLIC DomainNodeDerivative_NumberOfVersionsGet

  PUBLIC DomainNodeDerivative_PartialDerivativeIndexGet
  
  PUBLIC DomainNodeDerivative_VersionNumberGet
  
  PUBLIC DomainNodes_DerivativeGlobalIndexGet

  PUBLIC DomainNodes_DerivativeNumberOfVersionsGet

  PUBLIC DomainNodes_DerivativePartialIndexGet

  PUBLIC DomainNodes_DerivativeVersionNumberGet

  PUBLIC DomainNodes_GlobalNodeNumberGet

  PUBLIC DomainNodes_LocalNodeNumberGet

  PUBLIC DomainNodes_NodeBoundaryNodeGet

  PUBLIC DomainNodes_NodeCheckExists

  PUBLIC DomainNodes_NodeFaceNumberGet

  PUBLIC DomainNodes_NodeLineNumberGet

  PUBLIC DomainNodes_NodeGet

  PUBLIC DomainNodes_NodeNumberOfDerivativesGet

  PUBLIC DomainNodes_NodeNumberOfFacesGet

  PUBLIC DomainNodes_NodeNumberOfLinesGet

  PUBLIC DomainNodes_NodeNumberOfSurroundingElementsGet
  
  PUBLIC DomainNodes_NodeSurroundingElementGet
  
  PUBLIC DomainNodes_NumberOfNodesGet
  
  PUBLIC DomainNodes_TotalNumberOfNodesGet

  PUBLIC DomainNodes_UserNodeNumberGet
  
  PUBLIC DomainTopology_DomainGet

  PUBLIC DomainTopology_DomainDofsGet

  PUBLIC DomainTopology_DomainElementsGet

  PUBLIC DomainTopology_DomainFacesGet

  PUBLIC DomainTopology_DomainLinesGet
  
  PUBLIC DomainTopology_LocalElementBasisGet

  PUBLIC DomainTopology_DomainNodesGet
  
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

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(decomposer)) CALL FlagError("Decomposer is not associated.",err,error,*999)
#endif    

    IF(.NOT.decomposer%decomposerFinished) THEN
      localError="Decomposer number "//TRIM(NumberToVString(decomposer%userNumber,"*",err,error))
      IF(ASSOCIATED(decomposer%region)) localError=localError// &
        & " of region number "//TRIM(NumberToVString(decomposer%region%userNumber,"*",err,error))
      localError=localError//" has not been finished."
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

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(decomposer)) CALL FlagError("Decomposer is not associated.",err,error,*999)
#endif    

    IF(decomposer%decomposerFinished) THEN
      localError="Decomposer number "//TRIM(NumberToVString(decomposer%userNumber,"*",err,error))
      IF(ASSOCIATED(decomposer%region)) localError=localError// &
        & " of region number "//TRIM(NumberToVString(decomposer%region%userNumber,"*",err,error))
      localError=localError//" has already been finished."
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

  !>Returns the decomposition for a decomposer index
  SUBROUTINE Decomposer_DecompositionGet(decomposer,decompositionIndex,decomposition,err,error,*)

    !Argument variables
    TYPE(DecomposerType), POINTER :: decomposer !<A pointer to the decomposer to get the decomposition for
    INTEGER(INTG), INTENT(IN) :: decompositionIndex !<The index of the decomposition to get.
    TYPE(DecompositionType), POINTER :: decomposition !<On return, the specified decomposition. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("Decomposer_DecompositionGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(decomposition)) CALL FlagError("Decomposition is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(decomposer)) CALL FlagError("Decomposer is not associated.",err,error,*999)
    IF(decompositionIndex<1.OR.decompositionIndex>decomposer%numberOfDecompositions) THEN
      localError="The specified decompoisition index of "//TRIM(NumberToVString(decompositionIndex,"*",err,error))// &
        & " is invalid. The decomposition index should be >= 1 and <= "// &
        & TRIM(NumberToVString(decomposer%numberOfDecompositions,"*",err,error))// &
        & " of decomposer number "//TRIM(NumberToVString(decomposer%userNumber,"*",err,error))
      IF(ASSOCIATED(decomposer%region)) localError=localError// &
        & " of region number "//TRIM(NumberToVString(decomposer%region%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(decomposer%decompositions)) THEN
      localError="Decompositions is not allocated for decomposer number "// &
        & TRIM(NumberToVString(decomposer%userNumber,"*",err,error))
      IF(ASSOCIATED(decomposer%region)) localError=localError// &
        & " of region number "//TRIM(NumberToVString(decomposer%region%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    decomposition=>decomposer%decompositions(decompositionIndex)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(decomposition)) THEN
      localError="The decomposition is not associated for decomposition index "// &
        & TRIM(NumberToVString(decompositionIndex,"*",err,error))// &
        & " of decomposer number "//TRIM(NumberToVString(decomposer%userNumber,"*",err,error))
      IF(ASSOCIATED(decomposer%region)) localError=localError// &
        & " of region number "//TRIM(NumberToVString(decomposer%region%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
 
    EXITS("Decomposer_DecompositionGet")
    RETURN
999 NULLIFY(decomposition)
998 ERRORSEXITS("Decomposer_DecompositionGet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_DecompositionGet

  !
  !================================================================================================================================
  !

  !>Returns the decomposer graph for a decomposer 
  SUBROUTINE Decomposer_DecomposerGraphGet(decomposer,decomposerGraph,err,error,*)

    !Argument variables
    TYPE(DecomposerType), POINTER :: decomposer !<A pointer to the decomposer to get the decomposition for
    TYPE(DecomposerGraphType), POINTER :: decomposerGraph !<On return, the decomposer graph for the decomposer. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("Decomposer_DecomposerGraphGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(decomposerGraph)) CALL FlagError("Decomposer graph is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(decomposer)) CALL FlagError("Decomposer is not associated.",err,error,*999)
#endif    

    decomposerGraph=>decomposer%decomposerGraph

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(decomposerGraph)) THEN
      localError="The decomposer graph is not associated for decomposer number "// &
        & TRIM(NumberToVString(decomposer%userNumber,"*",err,error))
      IF(ASSOCIATED(decomposer%region)) localError=localError// &
        & " of region number "//TRIM(NumberToVString(decomposer%region%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
 
    EXITS("Decomposer_DecomposerGraphGet")
    RETURN
999 NULLIFY(decomposerGraph)
998 ERRORSEXITS("Decomposer_DecomposerGraphGet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_DecomposerGraphGet

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
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("Decomposer_RegionGet",err,error,*998 )

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(region)) CALL FlagError("Region is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(decomposer)) CALL FlagError("Decomposer is not associated.",err,error,*999)
#endif    

    region=>decomposer%region

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(region)) THEN
      localError="The region associated with decomposer number "// &
        & TRIM(NumberToVString(decomposer%userNumber,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
 
    EXITS("Decomposer_RegionGet")
    RETURN
999 NULLIFY(region)
998 ERRORSEXITS("Decomposer_RegionGet",err,error)
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

    ENTERS("Decomposer_UserNumberFind",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(decomposer)) CALL FlagError("Decomposer is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(region%decomposers)) THEN
      localError="The decomposers on region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
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
999 NULLIFY(decomposer)
998 ERRORSEXITS("Decomposer_UserNumberFind",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_UserNumberFind

  !
  !================================================================================================================================
  !

  !>Returns the decomposer for a decomposer graph.
  SUBROUTINE DecomposerGraph_DecomposerGet(decomposerGraph,decomposer,err,error,*)

    !Argument variables
    TYPE(DecomposerGraphType), POINTER :: decomposerGraph !<A pointer to the decomposer graph to get the decomposer for
    TYPE(DecomposerType), POINTER :: decomposer !<On return, the decomposer for the specified decomposer graph. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DecomposerGraph_DecomposerGet",err,error,*998)

#ifdef WITH_PRECHECKS
    !Check input arguments
    IF(ASSOCIATED(decomposer)) CALL FlagError("Decomposer is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(decomposerGraph)) CALL FlagError("Decomposer graph is not associated.",err,error,*999)
#endif    

    decomposer=>decomposerGraph%decomposer

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(decomposer)) CALL FlagError("The decomposer is not associated for the decomposer graph.",err,error,*999)
#endif    
 
    EXITS("DecomposerGraph_DecomposerGet")
    RETURN
999 NULLIFY(decomposer)
998 ERRORSEXITS("DecomposerGraph_DecomposereGet",err,error)
    RETURN 1
    
  END SUBROUTINE DecomposerGraph_DecomposerGet

  !
  !================================================================================================================================
  !

  !>Returns a root graph node for a decomposer root node index.
  SUBROUTINE DecomposerGraph_RootNodeGet(decomposerGraph,rootNodeIndex,rootNode,err,error,*)

    !Argument variables
    TYPE(DecomposerGraphType), POINTER :: decomposerGraph !<A pointer to the decomposer graph to get the root node for
    INTEGER(INTG), INTENT(IN) :: rootNodeIndex !<The index of the root graph node to get.
    TYPE(DecomposerGraphNodeType), POINTER :: rootNode !<On return, the specified root graph node. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("DecomposerGraph_RootNodeGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(rootNode)) CALL FlagError("Root node is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(decomposerGraph)) CALL FlagError("Decomposer graph is not associated.",err,error,*999)
    IF(rootNodeIndex<1.OR.rootNodeIndex>decomposerGraph%numberOfGraphRoots) THEN
      localError="The specified root node index of "//TRIM(NumberToVString(rootNodeIndex,"*",err,error))// &
        & " is invalid. The root node index should be >= 1 and <= "// &
        & TRIM(NumberToVString(decomposerGraph%numberOfGraphRoots,"*",err,error))// &
        & " for the decomposer graph"
      IF(ASSOCIATED(decomposerGraph%decomposer)) THEN
        localError=localError// " of decomposer number "// &
          & TRIM(NumberToVString(decomposerGraph%decomposer%userNumber,"*",err,error))
        IF(ASSOCIATED(decomposerGraph%decomposer%region)) localError=localError// &
          & " of region number "//TRIM(NumberToVString(decomposerGraph%decomposer%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(decomposerGraph%graphRoots)) THEN
      localError="Graph roots is not allocated for the decomposer graph"
      IF(ASSOCIATED(decomposerGraph%decomposer)) THEN
        localError=localError// " of decomposer number "// &
          & TRIM(NumberToVString(decomposerGraph%decomposer%userNumber,"*",err,error))
        IF(ASSOCIATED(decomposerGraph%decomposer%region)) localError=localError// &
          & " of region number "//TRIM(NumberToVString(decomposerGraph%decomposer%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    rootNode=>decomposerGraph%graphRoots(rootNodeIndex)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(rootNode)) THEN
      localError="The root node is not associated for root node index "// &
        & TRIM(NumberToVString(rootNodeIndex,"*",err,error))// &
        & " of the decomposer graph"
      IF(ASSOCIATED(decomposerGraph%decomposer)) THEN
        localError=localError// " of decomposer number "// &
          & TRIM(NumberToVString(decomposerGraph%decomposer%userNumber,"*",err,error))
        IF(ASSOCIATED(decomposerGraph%decomposer%region)) localError=localError// &
          & " of region number "//TRIM(NumberToVString(decomposerGraph%decomposer%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
 
    EXITS("DecomposerGraph_RootNodeGet")
    RETURN
999 NULLIFY(rootNode)
998 ERRORSEXITS("DecomposerGraph_RootNodeGet",err,error)
    RETURN 1
    
  END SUBROUTINE DecomposerGraph_RootNodeGet

  !
  !================================================================================================================================
  !

  !>Returns the decomposition for a decomposer graph link
  SUBROUTINE DecomposerGraphLink_DecompositionGet(graphLink,decomposition,err,error,*)

    !Argument variables
    TYPE(DecomposerGraphLinkType), POINTER :: graphLink !<A pointer to the decomposer graph link to get the decomposition for
    TYPE(DecompositionType), POINTER :: decomposition !<On return, the decomposer graph link decomposition. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("DecomposerGraphLink_DecompositionGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments    
    IF(ASSOCIATED(decomposition)) CALL FlagError("Decomposition is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(graphLink)) CALL FlagError("Graph link is not associated.",err,error,*999)
#endif    

    decomposition=>graphLink%decomposition

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(decomposition)) THEN
      localError="The decomposition for the graph link is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    EXITS("DecomposerGraphLink_DecompositionGet")
    RETURN
999 NULLIFY(decomposition)
998 ERRORSEXITS("DecomposerGraphLink_DecompositionGet",err,error)
    RETURN 1
    
  END SUBROUTINE DecomposerGraphLink_DecompositionGet

  !
  !================================================================================================================================
  !

  !>Returns the linked node for a decomposer graph link
  SUBROUTINE DecomposerGraphLink_LinkedNodeGet(graphLink,linkedNode,err,error,*)

    !Argument variables
    TYPE(DecomposerGraphLinkType), POINTER :: graphLink !<A pointer to the decomposer graph link to get the linked node for
    TYPE(DecomposerGraphNodeType), POINTER :: linkedNode !<On return, the decomposer graph link linked node. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("DecomposerGraphLink_LinkedNodeGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(linkedNode)) CALL FlagError("Linked node is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(graphLink)) CALL FlagError("Graph link is not associated.",err,error,*999)
#endif    

    linkedNode=>graphLink%linkedGraphNode

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(linkedNode)) THEN
      localError="The linked node for the graph link is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
 
    EXITS("DecomposerGraphLink_LinkedNodeGet")
    RETURN
999 NULLIFY(linkedNode)
998 ERRORSEXITS("DecomposerGraphLink_LinkedNodeGet",err,error)
    RETURN 1
    
  END SUBROUTINE DecomposerGraphLink_LinkedNodeGet

  !
  !================================================================================================================================
  !

  !>Returns the decomposer graph for a decomposer graph node
  SUBROUTINE DecomposerGraphNode_DecomposerGraphGet(graphNode,decomposerGraph,err,error,*)

    !Argument variables
    TYPE(DecomposerGraphNodeType), POINTER :: graphNode !<A pointer to the decomposer graph node to get the decomposer graph for
    TYPE(DecomposerGraphType), POINTER :: decomposerGraph !<On return, the decomposer graph for the decomposer graph node. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DecomposerGraphNode_DecomposerGraphGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(decomposerGraph)) CALL FlagError("Decomposer graph is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(graphNode)) CALL FlagError("Graph node is not associated.",err,error,*999)
#endif    

    decomposerGraph=>graphNode%decomposerGraph

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(decomposerGraph)) THEN
      localError="The decomposer graph for the graph node is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
 
    EXITS("DecomposerGraphLink_DecomposerGraphGet")
    RETURN
999 NULLIFY(decomposerGraph)
998 ERRORSEXITS("DecomposerGraphLink_DecomposerGraphGet",err,error)
    RETURN 1
    
  END SUBROUTINE DecomposerGraphNode_DecomposerGraphGet

  !
  !================================================================================================================================
  !

  !>Returns the decomposition for a decomposer graph node
  SUBROUTINE DecomposerGraphNode_DecompositionGet(graphNode,decomposition,err,error,*)

    !Argument variables
    TYPE(DecomposerGraphNodeType), POINTER :: graphNode !<A pointer to the decomposer graph node to get the decomposition for
    TYPE(DecompositionType), POINTER :: decomposition !<On return, the decomposer graph node decomposition. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DecomposerGraphNode_DecompositionGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(decomposition)) CALL FlagError("Decomposition is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(graphNode)) CALL FlagError("Graph node is not associated.",err,error,*999)
#endif    

    decomposition=>graphNode%decomposition

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(decomposition)) THEN
      localError="The decomposition for the graph node is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
 
    EXITS("DecomposerGraphNode_DecompositionGet")
    RETURN
999 NULLIFY(decomposition)
998 ERRORSEXITS("DecomposerGraphNode_DecompositionGet",err,error)
    RETURN 1
    
  END SUBROUTINE DecomposerGraphNode_DecompositionGet

  !
  !================================================================================================================================
  !

  !>Returns a decomposer graph link for a graph node index.
  SUBROUTINE DecomposerGraphNode_DecomposerGraphLinkGet(graphNode,graphLinkIndex,graphLink,err,error,*)

    !Argument variables
    TYPE(DecomposerGraphNodeType), POINTER :: graphNode !<A pointer to the decomposer graph node to get the graph linkfor
    INTEGER(INTG), INTENT(IN) :: graphLinkIndex !<The index of the graph link to get.
    TYPE(DecomposerGraphLinkType), POINTER :: graphLink !<On return, the specified graph link. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DecomposerGraphNode_GraphLinkGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(graphLink)) CALL FlagError("Graph link is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(graphNode)) CALL FlagError("Graph node is not associated.",err,error,*999)
    IF(graphLinkIndex<1.OR.graphLinkIndex>graphNode%numberOfGraphLinks) THEN
      localError="The specified graph link index of "//TRIM(NumberToVString(graphLinkIndex,"*",err,error))// &
        & " is invalid. The graph link index should be >= 1 and <= "// &
        & TRIM(NumberToVString(graphNode%numberOfGraphLinks,"*",err,error))// &
        & " for the decomposer graph node"
      IF(ASSOCIATED(graphNode%decomposerGraph)) THEN
        localError=localError//" of the decomposer graph"
        IF(ASSOCIATED(graphNode%decomposerGraph%decomposer)) THEN
          localError=localError// " of decomposer number "// &
            & TRIM(NumberToVString(graphNode%decomposerGraph%decomposer%userNumber,"*",err,error))
          IF(ASSOCIATED(graphNode%decomposerGraph%decomposer%region)) localError=localError// &
            & " of region number "// &
            & TRIM(NumberToVString(graphNode%decomposerGraph%decomposer%region%userNumber,"*",err,error))
        ENDIF
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(graphNode%graphLinks)) THEN
      localError="Graph links is not allocated for the decomposer graph node"
      IF(ASSOCIATED(graphNode%decomposerGraph)) THEN
        localError=localError//" of the decomposer graph"
        IF(ASSOCIATED(graphNode%decomposerGraph%decomposer)) THEN
          localError=localError// " of decomposer number "// &
            & TRIM(NumberToVString(graphNode%decomposerGraph%decomposer%userNumber,"*",err,error))
          IF(ASSOCIATED(graphNode%decomposerGraph%decomposer%region)) localError=localError// &
            & " of region number "// &
            & TRIM(NumberToVString(graphNode%decomposerGraph%decomposer%region%userNumber,"*",err,error))
        ENDIF
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    graphLink=>graphNode%graphLinks(graphLinkIndex)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(graphLink)) THEN
      localError="The graph link is not associated for graph link index "// &
        & TRIM(NumberToVString(graphLinkIndex,"*",err,error))// &
        & " of the decomposer graph node"
      IF(ASSOCIATED(graphNode%decomposerGraph)) THEN
        localError=localError//" of the decomposer graph"
        IF(ASSOCIATED(graphNode%decomposerGraph%decomposer)) THEN
          localError=localError// " of decomposer number "// &
            & TRIM(NumberToVString(graphNode%decomposerGraph%decomposer%userNumber,"*",err,error))
          IF(ASSOCIATED(graphNode%decomposerGraph%decomposer%region)) localError=localError// &
            & " of region number "// &
            & TRIM(NumberToVString(graphNode%decomposerGraph%decomposer%region%userNumber,"*",err,error))
        ENDIF
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
 
    EXITS("DecomposerGraphNode_DecomposerGraphLinkGet")
    RETURN
999 NULLIFY(graphLink)
998 ERRORSEXITS("DecomposerGraph_DecomposerGraphLinkGet",err,error)
    RETURN 1
    
  END SUBROUTINE DecomposerGraphNode_DecomposerGraphLinkGet

  !
  !=================================================================================================================================
  !

  !>Assert that a decomposition has been decomposed
  SUBROUTINE Decomposition_AssertIsDecomposed(decomposition,err,error,*)

    !Argument Variables
    TYPE(DecompositionType), POINTER, INTENT(INOUT) :: decomposition !<The decomposition to assert the decomposed status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Decomposition_AssertIsDecomposed",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
#endif    

    IF(.NOT.ASSOCIATED(decomposition%decomposer)) THEN
      localError="Decomposition number "//TRIM(NumberToVString(decomposition%userNumber,"*",err,error))// &
        & " has not been decomposed."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Decomposition_AssertIsDecomposed")
    RETURN
999 ERRORSEXITS("Decomposition_AssertIsDecomposed",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_AssertIsDecomposed

  !
  !=================================================================================================================================
  !

  !>Assert that a decomposition has not been decomposed
  SUBROUTINE Decomposition_AssertNotDecomposed(decomposition,err,error,*)

    !Argument Variables
    TYPE(DecompositionType), POINTER, INTENT(INOUT) :: decomposition !<The decomposition to assert the decomposition status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Decomposition_AssertNotDecomposed",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
#endif    

    IF(ASSOCIATED(decomposition%decomposer)) THEN
      localError="Decomposition number "//TRIM(NumberToVString(decomposition%userNumber,"*",err,error))// &
        & " has already been decomposed."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Decomposition_AssertNotDecomposed")
    RETURN
999 ERRORSEXITS("Decomposition_AssertNotDecomposed",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_AssertNotDecomposed

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

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
#endif    

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

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
#endif    

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
  !=================================================================================================================================
  !

  !>Assert that a decomposition has the calculate faces flag set
  SUBROUTINE Decomposition_AssertCalculateFaces(decomposition,err,error,*)

    !Argument Variables
    TYPE(DecompositionType), POINTER, INTENT(INOUT) :: decomposition !<The decomposition to assert the calculate faces status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Decomposition_AssertCalculateFaces",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
#endif    

    IF(.NOT.decomposition%calculateFaces) THEN
      localError="Decomposition number "//TRIM(NumberToVString(decomposition%userNumber,"*",err,error))// &
        & " does not have the calculate faces flag set."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Decomposition_AssertCalculateFaces")
    RETURN
999 ERRORSEXITS("Decomposition_AssertCalculateFaces",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_AssertCalculateFaces

  !
  !=================================================================================================================================
  !

  !>Assert that a decomposition does not have the calculate faces flag set
  SUBROUTINE Decomposition_AssertNotCalculateFaces(decomposition,err,error,*)

    !Argument Variables
    TYPE(DecompositionType), POINTER, INTENT(INOUT) :: decomposition !<The decomposition to assert the calculate faces status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Decomposition_AssertNotCalculateFaces",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
#endif    

    IF(decomposition%calculateFaces) THEN
      localError="Decomposition number "//TRIM(NumberToVString(decomposition%userNumber,"*",err,error))// &
        & " has the calculate faces flag set."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Decomposition_AssertNotCalculateFaces")
    RETURN
999 ERRORSEXITS("Decomposition_AssertNotCalculateFaces",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_AssertNotCalculateFaces

  !
  !=================================================================================================================================
  !

  !>Assert that a decomposition has the calculate lines flag set
  SUBROUTINE Decomposition_AssertCalculateLines(decomposition,err,error,*)

    !Argument Variables
    TYPE(DecompositionType), POINTER, INTENT(INOUT) :: decomposition !<The decomposition to assert the calculate lines status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Decomposition_AssertCalculateLines",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
#endif    

    IF(.NOT.decomposition%calculateLines) THEN
      localError="Decomposition number "//TRIM(NumberToVString(decomposition%userNumber,"*",err,error))// &
        & " does not have the calculate lines flag set."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Decomposition_AssertCalculateLines")
    RETURN
999 ERRORSEXITS("Decomposition_AssertCalculateLines",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_AssertCalculateLines

  !
  !=================================================================================================================================
  !

  !>Assert that a decomposition does not have the calculate lines flag set
  SUBROUTINE Decomposition_AssertNotCalculateLines(decomposition,err,error,*)

    !Argument Variables
    TYPE(DecompositionType), POINTER, INTENT(INOUT) :: decomposition !<The decomposition to assert the calculate lines status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Decomposition_AssertNotCalculateLines",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
#endif    

    IF(decomposition%calculateLines) THEN
      localError="Decomposition number "//TRIM(NumberToVString(decomposition%userNumber,"*",err,error))// &
        & " has the calculate lines flag set."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Decomposition_AssertNotCalculateLines")
    RETURN
999 ERRORSEXITS("Decomposition_AssertNotCalculateLines",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_AssertNotCalculateLines

  !
  !=================================================================================================================================
  !

  !>Return the calculate faces flag for a decomposition
  SUBROUTINE Decomposition_CalculateFacesGet(decomposition,calculateFaces,err,error,*)

    !Argument Variables
    TYPE(DecompositionType), POINTER, INTENT(INOUT) :: decomposition !<The decomposition to get the calculate faces flag for
    LOGICAL, INTENT(OUT) :: calculateFaces !<On return, the calculate faces flag for the decompsotion
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Decomposition_CalculateFacesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
#endif

    calculateFaces=decomposition%calculateFaces

    EXITS("Decomposition_CalculateFacesGet")
    RETURN
999 ERRORSEXITS("Decomposition_CalculateFacesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_CalculateFacesGet

  !
  !=================================================================================================================================
  !

  !>Return the calculate lines flag for a decomposition
  SUBROUTINE Decomposition_CalculateLinesGet(decomposition,calculateLines,err,error,*)

    !Argument Variables
    TYPE(DecompositionType), POINTER, INTENT(INOUT) :: decomposition !<The decomposition to get the calculate lines flag for
    LOGICAL, INTENT(OUT) :: calculateLines !<On return, the calculate lines flag for the decompsotion
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Decomposition_CalculateLinesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
#endif

    calculateLines=decomposition%calculateLines

    EXITS("Decomposition_CalculateLinesGet")
    RETURN
999 ERRORSEXITS("Decomposition_CalculateLinesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_CalculateLinesGet

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

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
    IF(ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is already associated.",err,error,*999)
#endif
    
    IF(ASSOCIATED(decomposition%region)) THEN
      coordinateSystem=>decomposition%region%coordinateSystem
    ELSE IF(ASSOCIATED(decomposition%INTERFACE)) THEN
      coordinateSystem=>decomposition%interface%coordinateSystem
    ELSE
      CALL FlagError("Decomposition is not associated with a region or an interface.",err,error,*999)
    ENDIF

#ifdef WITH_POSTCHECKS    
    !Check coordinate system is associated.
    IF(.NOT.ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is not associated.",err,error,*999)
#endif    
    
    EXITS("Decomposition_CoordinateSystemGet")
    RETURN
999 ERRORSEXITS("Decomposition_CoordinateSystemGet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_CoordinateSystemGet

  !
  !================================================================================================================================
  !

  !>Gets the decompositions from a decomposition. 
  SUBROUTINE Decomposition_DecompositionsGet(decomposition,decompositions,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<The decomposition to get the decompositions for.
    TYPE(DecompositionsType), POINTER :: decompositions !<On exit, a pointer to the decompositions for the decomposition. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Decomposition_DecompositionsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(decompositions)) CALL FlagError("Decompositions is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
#endif    

    !Get the decompositions
    decompositions=>decomposition%decompositions

#ifdef WITH_POSTCHECKS    
    !Check decompositions is associated.
    IF(.NOT.ASSOCIATED(decompositions)) CALL FlagError("Decomposition decompositions is not associated.",err,error,*999)
#endif    
    
    EXITS("Decomposition_DecompositionsGet")
    RETURN
999 NULLIFY(decompositions)
998 ERRORSEXITS("Decomposition_DecompositionsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_DecompositionsGet

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
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Decomposition_DomainGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(domain)) CALL FlagError("Domain is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
    IF(meshComponentNumber<0.OR.meshComponentNumber>decomposition%numberOfComponents) THEN
      localError="The specified mesh component number of "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
        & " is invalid. The mesh component number must be >= 0 and <= "// &
        & TRIM(NumberToVString(decomposition%numberOfComponents,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(decomposition%domain)) CALL FlagError("Decomposition domain is not allocated.",err,error,*999)
#endif    

    !Get the domain
    IF(meshComponentNumber==0) THEN
      meshComponent=decomposition%meshComponentNumber
    ELSE
      meshComponent=meshComponentNumber
    ENDIF
    domain=>decomposition%domain(meshComponent)%ptr

#ifdef WITH_POSTCHECKS    
    !Check domain is associated.
    IF(.NOT.ASSOCIATED(domain)) THEN
      localError="Domain is not associated for mesh component "//TRIM(NumberToVString(meshComponent,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("Decomposition_DomainGet")
    RETURN
999 ERRORSEXITS("Decomposition_DomainGet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_DomainGet

  !
  !================================================================================================================================
  !

  !>Returns the interface for a decomposition
  SUBROUTINE Decomposition_InterfaceGet(decomposition,interface,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to get the interface for
    TYPE(InterfaceType), POINTER :: interface !<On return, the decomposition interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("Decomposition_InterfaceGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(interface)) CALL FlagError("Interface is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
#endif    

    INTERFACE=>decomposition%INTERFACE

#ifdef WITH_POSTCHECKS      
    IF(.NOT.ASSOCIATED(INTERFACE)) THEN
      localError="The interface for decomposition number "// &
        & TRIM(NumberToVString(decomposition%userNumber,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
 
    EXITS("Decomposition_InterfaceGet")
    RETURN
999 NULLIFY(interface)
998 ERRORSEXITS("Decomposition_InterfaceGet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_InterfaceGet

  !
  !================================================================================================================================
  !

  !>Determines if the given decomposition is an interface decomposition or not. 
  SUBROUTINE Decomposition_IsInterfaceDecomposition(decomposition,interfaceDecomposition,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to determine if it is an interface decomposition or not.
    LOGICAL :: interfaceDecomposition !<On exit, .TRUE. if the given decomposition is in an interface region, .FALSE. if not. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Decomposition_IsInterfaceDecomposition",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
#endif    

    interfaceDecomposition = ASSOCIATED(decomposition%interface)
    
    EXITS("Decomposition_IsInterfaceDecomposition")
    RETURN
999 ERRORSEXITS("Decomposition_IsInterfaceDecomposition",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_IsInterfaceDecomposition

  !
  !================================================================================================================================
  !

  !>Determines if the given decomposition is a region decomposition or not. 
  SUBROUTINE Decomposition_IsRegionDecomposition(decomposition,regionDecomposition,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to determine if it is an region decomposition or not.
    LOGICAL :: regionDecomposition !<On exit, .TRUE. if the given decomposition is in a region, .FALSE. if not. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Decomposition_IsRegionDecomposition",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
#endif    

    regionDecomposition = ASSOCIATED(decomposition%region)
    
    EXITS("Decomposition_IsRegionDecomposition")
    RETURN
999 ERRORSEXITS("Decomposition_IsRegionDecomposition",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_IsRegionDecomposition

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
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Decomposition_MeshGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(mesh)) CALL FlagError("Mesh is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
#endif    
      
    mesh=>decomposition%mesh

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(mesh)) THEN
      localError="The mesh associated with decomposition number "// &
        & TRIM(NumberToVString(decomposition%userNumber,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
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

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(region)) CALL FlagError("Region is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
#endif    

    IF(ASSOCIATED(decomposition%region)) THEN
      region=>decomposition%region
    ELSE IF(ASSOCIATED(decomposition%interface)) THEN
      region=>decomposition%interface%parentRegion
    ELSE
      localError="Decomposition number "//TRIM(NumberToVString(decomposition%userNumber,"*",err,error))// &
        & " is not associated with a region or an interface."
      CALL FlagError(localError,err,error,*999)
    ENDIF

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(region)) THEN
      localError="The region associated with decomposition number "// &
        & TRIM(NumberToVString(decomposition%userNumber,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    EXITS("Decomposition_RegionGet")
    RETURN
999 ERRORSEXITS("Decomposition_RegionGet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_RegionGet

  !
  !================================================================================================================================
  !

  !>Gets a decomposition topology from a decomposition.
  SUBROUTINE Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<The decomposition to get the domain for.
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<On exit, a pointer to the decomposition topology. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif
    
    ENTERS("Decomposition_DecompositionTopologyGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
    IF(ASSOCIATED(decompositionTopology)) CALL FlagError("Decomposition topology is already associated.",err,error,*999)
#endif    
 
    !Get the decomposition topology
    decompositionTopology=>decomposition%topology

#ifdef WITH_POSTCHECKS    
    !Check decompositionTopology is associated.
    IF(.NOT.ASSOCIATED(decompositionTopology)) THEN
      localError="Decomposition topology is not associated for decomposition "// &
        & TRIM(NumberToVString(decomposition%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    EXITS("Decomposition_DecompositionTopologyGet")
    RETURN
999 ERRORSEXITS("Decomposition_DecompositionTopologyGet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_DecompositionTopologyGet

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

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)
    IF(ASSOCIATED(decomposition)) CALL FlagError("Decomposition is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(mesh%decompositions)) THEN
      localError="The decompositions on mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))// &
        & " are not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    !Get the decomposition from the user number
    NULLIFY(decomposition)
    IF(ALLOCATED(mesh%decompositions%decompositions)) THEN
      DO decompositionIdx=1,mesh%decompositions%numberOfDecompositions
        IF(.NOT.ASSOCIATED(mesh%decompositions%decompositions(decompositionIdx)%ptr)) THEN
          localError="The decomposition pointer in mesh decompositions is not associated for decomposition index "// &
            & TRIM(NumberToVString(decompositionIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(mesh%decompositions%decompositions(decompositionIdx)%ptr%userNumber==userNumber) THEN
          decomposition=>mesh%decompositions%decompositions(decompositionIdx)%ptr
          EXIT
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

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(workGroup)) CALL FlagError("Work group is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
#endif    
      
    workGroup=>decomposition%workGroup

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("The decomposition work group is not associated.",err,error,*999)
#endif    
       
    EXITS("Decomposition_WorkGroupGet")
    RETURN
999 NULLIFY(workGroup)
998 ERRORSEXITS("Decomposition_WorkGroupGet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_WorkGroupGet

  !
  !================================================================================================================================
  !

  !>Checks that a user data point number exists. 
  SUBROUTINE DecompositionDataPoints_DataPointCheckExists(decompositionDataPoints,userDataPointNumber,userDataPointExists, &
    & localDataPointNumber,ghostDataPoint,err,error,*)

    !Argument variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionDataPoints !<A pointer to the decomposition data points to check the data point exists on
    INTEGER(INTG), INTENT(IN) :: userDataPointNumber !<The user data point number to check if it exists
    LOGICAL, INTENT(OUT) :: userDataPointExists !<On exit, is .TRUE. if the data point user number exists in the decomposition topology, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: localDataPointNumber !<On exit, if the data point exists the local number corresponding to the user data point number. If the data point does not exist then local number will be 0.
    LOGICAL, INTENT(OUT) :: ghostDataPoint !<On exit, is .TRUE. if the local data point (if it exists) is a ghost data point, .FALSE. if not.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(TreeNodeType), POINTER :: treeNode

    ENTERS("DecompositionDataPoints_DataPointCheckExists",ERR,error,*999)

    userDataPointExists=.FALSE.
    localDataPointNumber=0
    ghostDataPoint=.FALSE.
#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(decompositionDataPoints)) CALL FlagError("Decomposition data points is not associated.",err,error,*999)
#endif    
    
    NULLIFY(treeNode)
    CALL Tree_Search(decompositionDataPoints%dataPointsTree,userDataPointNumber,treeNode,err,error,*999)
    IF(ASSOCIATED(treeNode)) THEN
      CALL Tree_NodeValueGet(decompositionDataPoints%dataPointsTree,treeNode,localDataPointNumber,err,error,*999)
      userDataPointExists=.TRUE.
      ghostDataPoint=localDataPointNumber>decompositionDataPoints%numberOfDataPoints
    ENDIF
 
    EXITS("DecompositionDataPoints_DataPointCheckExists")
    RETURN
999 ERRORS("DecompositionDataPoints_DataPointCheckExists",err,error)
    EXITS("DecompositionDataPoints_DataPointCheckExists")
    RETURN 1

  END SUBROUTINE DecompositionDataPoints_DataPointCheckExists

  !
  !================================================================================================================================
  !

  !>Gets the global data point number corresponding to an element data point index for decomposition data points. 
  SUBROUTINE DecompositionDataPoints_ElementDataGlobalNumberGet(decompositionDataPoints,elementDataPointIdx,localElementNumber, &
    globalDataPointNumber,err,error,*)

    !Argument variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionDataPoints !<A pointer to the decomposition data points to get the global data point for
    INTEGER(INTG), INTENT(IN) :: elementDataPointIdx !<The element data point index to get the global data point number for
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The local element data point number containing the element data point to get the global data point number for
    INTEGER(INTG), INTENT(OUT) :: globalDataPointNumber !<On exit the global data point number corresponding to the element data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif

    ENTERS("DecompositionDataPoints_ElementDataGlobalNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(decompositionDataPoints)) CALL FlagError("Decomposition data points is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(decompositionDataPoints%elementDataPoints)) &
      & CALL FlagError("The element data points array is not allocated for the decomposition data points.",err,error,*999)
    IF(localElementNumber<1.OR.localElementNumber>SIZE(decompositionDataPoints%elementDataPoints,1)) THEN
      localError="The specified local element number of "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " is invalid for the decomposition data points. The local element number should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(decompositionDataPoints%elementDataPoints,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(elementDataPointIdx<1.OR.elementDataPointIdx>decompositionDataPoints%elementDataPoints(localElementNumber)% &
      & numberOfProjectedData) THEN
      localError="The specified element data point index of "//TRIM(NumberToVString(elementDataPointIdx,"*",err,error))// &
        & " is invalid for local element number "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " for the decomposition data points. The element data point index should be >= 1 and <= "// &
        & TRIM(NumberToVString(decompositionDataPoints%elementDataPoints(localElementNumber)%numberOfProjectedData, &
        & "*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(ALLOCATED(decompositionDataPoints%elementDataPoints(localElementNumber)%dataIndices)) THEN
      localError="The data indices array is not allocated for local element number "// &
        & TRIM(NumberToVString(localElementNumber,"*",err,error))//" of the decomposition data points."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    globalDataPointNumber=decompositionDataPoints%elementDataPoints(localElementNumber)%dataIndices(elementDataPointIdx)% &
      & globalNumber
    
    EXITS("DecompositionDataPoints_ElementDataGlobalNumberGet")
    RETURN
999 ERRORS("DecompositionDataPoints_ElementDataGlobalNumberGet",err,error)
    EXITS("DecompositionDataPoints_ElementDataGlobalNumberGet")
    RETURN 1

  END SUBROUTINE DecompositionDataPoints_ElementDataGlobalNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the local data point number corresponding to an element data point index for decomposition data points. 
  SUBROUTINE DecompositionDataPoints_ElementDataLocalNumberGet(decompositionDataPoints,elementDataPointIdx,localElementNumber, &
    localDataPointNumber,err,error,*)

    !Argument variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionDataPoints !<A pointer to the decomposition data points to get the local data point for
    INTEGER(INTG), INTENT(IN) :: elementDataPointIdx !<The element data point index to get the local data point number for
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The local element data point number containing the element data point to get the local data point number for
    INTEGER(INTG), INTENT(OUT) :: localDataPointNumber !<On exit the local data point number corresponding to the element data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif

    ENTERS("DecompositionDataPoints_ElementDataLocalNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(decompositionDataPoints)) CALL FlagError("Decomposition data points is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(decompositionDataPoints%elementDataPoints)) &
      & CALL FlagError("The element data points array is not allocated for the decomposition data points.",err,error,*999)
    IF(localElementNumber<1.OR.localElementNumber>SIZE(decompositionDataPoints%elementDataPoints,1)) THEN
      localError="The specified local element number of "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " is invalid for the decomposition data points. The local element number should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(decompositionDataPoints%elementDataPoints,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(elementDataPointIdx<1.OR.elementDataPointIdx>decompositionDataPoints%elementDataPoints(localElementNumber)% &
      & numberOfProjectedData) THEN
      localError="The specified element data point index of "//TRIM(NumberToVString(elementDataPointIdx,"*",err,error))// &
        & " is invalid for local element number "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " for the decomposition data points. The element data point index should be >= 1 and <= "// &
        & TRIM(NumberToVString(decompositionDataPoints%elementDataPoints(localElementNumber)%numberOfProjectedData, &
        & "*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(ALLOCATED(decompositionDataPoints%elementDataPoints(localElementNumber)%dataIndices)) THEN
      localError="The data indices array is not allocated for local element number "// &
        & TRIM(NumberToVString(localElementNumber,"*",err,error))//" of the decomposition data points."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    localDataPointNumber=decompositionDataPoints%elementDataPoints(localElementNumber)%dataIndices(elementDataPointIdx)%localNumber
    
    EXITS("DecompositionDataPoints_ElementDataLocalNumberGet")
    RETURN
999 ERRORS("DecompositionDataPoints_ElementDataLocalNumberGet",err,error)
    EXITS("DecompositionDataPoints_ElementDataLocalNumberGet")
    RETURN 1

  END SUBROUTINE DecompositionDataPoints_ElementDataLocalNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the data point numbers (local, global and user) corresponding to an element data point index for decomposition data points. 
  SUBROUTINE DecompositionDataPoints_ElementDataNumbersGet(decompositionDataPoints,elementDataPointIdx,localElementNumber, &
    localDataPointNumber,globalDataPointNumber,userDataPointNumber,err,error,*)

    !Argument variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionDataPoints !<A pointer to the decomposition data points to get the data point numbers for
    INTEGER(INTG), INTENT(IN) :: elementDataPointIdx !<The element data point index to get the data point numbers for
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The local element data point number containing the element data point to get the data point numbers for
    INTEGER(INTG), INTENT(OUT) :: localDataPointNumber !<On exit the local data point number corresponding to the element data point
    INTEGER(INTG), INTENT(OUT) :: globalDataPointNumber !<On exit the global data point number corresponding to the element data point
    INTEGER(INTG), INTENT(OUT) :: userDataPointNumber !<On exit the user data point number corresponding to the element data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif

    ENTERS("DecompositionDataPoints_ElementDataNumbersGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(decompositionDataPoints)) CALL FlagError("Decomposition data points is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(decompositionDataPoints%elementDataPoints)) &
      & CALL FlagError("The element data points array is not allocated for the decomposition data points.",err,error,*999)
    IF(localElementNumber<1.OR.localElementNumber>SIZE(decompositionDataPoints%elementDataPoints,1)) THEN
      localError="The specified local element number of "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " is invalid for the decomposition data points. The local element number should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(decompositionDataPoints%elementDataPoints,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(elementDataPointIdx<1.OR.elementDataPointIdx>decompositionDataPoints%elementDataPoints(localElementNumber)% &
      & numberOfProjectedData) THEN
      localError="The specified element data point index of "//TRIM(NumberToVString(elementDataPointIdx,"*",err,error))// &
        & " is invalid for local element number "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " for the decomposition data points. The element data point index should be >= 1 and <= "// &
        & TRIM(NumberToVString(decompositionDataPoints%elementDataPoints(localElementNumber)%numberOfProjectedData, &
        & "*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(ALLOCATED(decompositionDataPoints%elementDataPoints(localElementNumber)%dataIndices)) THEN
      localError="The data indices array is not allocated for local element number "// &
        & TRIM(NumberToVString(localElementNumber,"*",err,error))//" of the decomposition data points."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    localDataPointNumber=decompositionDataPoints%elementDataPoints(localElementNumber)% &
      & dataIndices(elementDataPointIdx)%localNumber
    globalDataPointNumber=decompositionDataPoints%elementDataPoints(localElementNumber)% &
      & dataIndices(elementDataPointIdx)%globalNumber
    userDataPointNumber=decompositionDataPoints%elementDataPoints(localElementNumber)% &
      & dataIndices(elementDataPointIdx)%userNumber
    
    EXITS("DecompositionDataPoints_ElementDataNumbersGet")
    RETURN
999 ERRORS("DecompositionDataPoints_ElementDataNumbersGet",err,error)
    EXITS("DecompositionDataPoints_ElementDataNumbersGet")
    RETURN 1

  END SUBROUTINE DecompositionDataPoints_ElementDataNumbersGet

  !
  !================================================================================================================================
  !

  !>Gets the user data point number corresponding to an element data point index for decomposition data points. 
  SUBROUTINE DecompositionDataPoints_ElementDataUserNumberGet(decompositionDataPoints,elementDataPointIdx,localElementNumber, &
    userDataPointNumber,err,error,*)

    !Argument variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionDataPoints !<A pointer to the decomposition data points to get the user data point for
    INTEGER(INTG), INTENT(IN) :: elementDataPointIdx !<The element data point index to get the user data point number for
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The local element data point number containing the element data point to get the user data point number for
    INTEGER(INTG), INTENT(OUT) :: userDataPointNumber !<On exit the user data point number corresponding to the element data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif

    ENTERS("DecompositionDataPoints_ElementDataUserNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(decompositionDataPoints)) CALL FlagError("Decomposition data points is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(decompositionDataPoints%elementDataPoints)) &
      & CALL FlagError("The element data points array is not allocated for the decomposition data points.",err,error,*999)
    IF(localElementNumber<1.OR.localElementNumber>SIZE(decompositionDataPoints%elementDataPoints,1)) THEN
      localError="The specified local element number of "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " is invalid for the decomposition data points. The local element number should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(decompositionDataPoints%elementDataPoints,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(elementDataPointIdx<1.OR.elementDataPointIdx>decompositionDataPoints%elementDataPoints(localElementNumber)% &
      & numberOfProjectedData) THEN
      localError="The specified element data point index of "//TRIM(NumberToVString(elementDataPointIdx,"*",err,error))// &
        & " is invalid for local element number "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " for the decomposition data points. The element data point index should be >= 1 and <= "// &
        & TRIM(NumberToVString(decompositionDataPoints%elementDataPoints(localElementNumber)%numberOfProjectedData, &
        & "*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(ALLOCATED(decompositionDataPoints%elementDataPoints(localElementNumber)%dataIndices)) THEN
      localError="The data indices array is not allocated for local element number "// &
        & TRIM(NumberToVString(localElementNumber,"*",err,error))//" of the decomposition data points."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    userDataPointNumber=decompositionDataPoints%elementDataPoints(localElementNumber)%dataIndices(elementDataPointIdx)%userNumber
    
    EXITS("DecompositionDataPoints_ElementDataUserNumberGet")
    RETURN
999 ERRORS("DecompositionDataPoints_ElementDataUserNumberGet",err,error)
    EXITS("DecompositionDataPoints_ElementDataUserNumberGet")
    RETURN 1

  END SUBROUTINE DecompositionDataPoints_ElementDataUserNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the number of data points in a local element number for decomposition data points. 
  SUBROUTINE DecompositionDataPoints_ElementNumberOfDataPointsGet(decompositionDataPoints,localElementNumber, &
    numberOfElementDataPoints,err,error,*)

    !Argument variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionDataPoints !<A pointer to the decomposition data points to get the number of data points for
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The local element data point number to get the number of element data points for
    INTEGER(INTG), INTENT(OUT) :: numberOfElementDataPoints !<On exit, the number of data points in the specified element in the decomposition data points.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif

    ENTERS("DecompositionDataPoints_ElementNumberOfDataPointsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(decompositionDataPoints)) CALL FlagError("Decomposition data points is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(decompositionDataPoints%elementDataPoints)) &
      & CALL FlagError("The element data points array is not allocated for the decomposition data points.",err,error,*999)
    IF(localElementNumber<1.OR.localElementNumber>SIZE(decompositionDataPoints%elementDataPoints,1)) THEN
      localError="The specified local element number of "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " is invalid for the decomposition data points. The local element number should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(decompositionDataPoints%elementDataPoints,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    numberOfElementDataPoints=decompositionDataPoints%elementDataPoints(localElementNumber)%numberOfProjectedData
    
    EXITS("DecompositionDataPoints_ElementNumberOfDataPointsGet")
    RETURN
999 ERRORS("DecompositionDataPoints_ElementNumberOfDataPointsGet",err,error)
    EXITS("DecompositionDataPoints_ElementNumberOfDataPointsGet")
    RETURN 1

  END SUBROUTINE DecompositionDataPoints_ElementNumberOfDataPointsGet

  !
  !================================================================================================================================
  !

  !>Gets a local data point number that corresponds to a user data point number from a decomposition. An error will be raised if the user data point number does not exist.
  SUBROUTINE DecompositionDataPoints_LocalDataPointNumberGet(decompositionDataPoints,userDataPointNumber,localDataPointNumber, &
    & ghostDataPoint,err,error,*)

    !Argument variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionDataPoints !<A pointer to the decomposition data points to get the data point on
    INTEGER(INTG), INTENT(IN) :: userDataPointNumber !<The user data point number to get
    INTEGER(INTG), INTENT(OUT) :: localDataPointNumber !<On exit, the local number corresponding to the user data point number.
    LOGICAL, INTENT(OUT) :: ghostDataPoint !<On exit, is .TRUE. if the local data point is a ghost data point, .FALSE. if not.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    LOGICAL :: dataPointExists
#ifdef WITH_POSTCHECKS    
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(VARYING_STRING) :: localError
#endif
    
    ENTERS("DecompositionDataPoints_LocalDataPointNumberGet",err,error,*999)

    CALL DecompositionDataPoints_DataPointCheckExists(decompositionDataPoints,userDataPointNumber,dataPointExists, &
      & localDataPointNumber,ghostDataPoint,err,error,*999)

#ifdef WITH_POSTCHECKS    
    IF(.NOT.dataPointExists) THEN
      decompositionTopology=>decompositionDataPoints%decompositionTopology
      IF(ASSOCIATED(decompositionTopology)) THEN
        decomposition=>decompositionTopology%decomposition
        IF(ASSOCIATED(decomposition)) THEN
          localError="The user data point number "//TRIM(NumberToVString(userDataPointNumber,"*",err,error))// &
          & " does not exist in decomposition number "//TRIM(NumberToVString(decomposition%userNumber,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ELSE
          localError="The user element number "//TRIM(NumberToVString(userDataPointNumber,"*",err,error))//" does not exist."
        ENDIF
      ELSE
        localError="The user data point number "//TRIM(NumberToVString(userDataPointNumber,"*",err,error))//" does not exist."
      ENDIF
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("DecompositionDataPoints_LocalDataPointNumberGet")
    RETURN
999 ERRORS("DecompositionDataPoints_LocalDataPointNumberGet",err,error)
    EXITS("DecompositionDataPoints_LocalDataPointNumberGet")
    RETURN 1

  END SUBROUTINE DecompositionDataPoints_LocalDataPointNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the adjacent element number corresponding to an xi coordinate index for a local element number in decomposition elements. 
  SUBROUTINE DecompositionElements_ElementAdjacentNumberGet(decompositionElements,adjacentElementIdx,xiCoordinateIdx, &
    & localElementNumber,adjacentElementNumber,err,error,*)

    !Argument variables
    TYPE(DecompositionElementsType), POINTER :: decompositionElements !<A pointer to the decomposition elements to get the adjacent element number for
    INTEGER(INTG), INTENT(IN) :: adjacentElementIdx !<The adjacent element index to get the adjacent element number for.
    INTEGER(INTG), INTENT(IN) :: xiCoordinateIdx !<The xi coordinate index to get the adjacent element number for.
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The local element number to get the adjacent element number for.
    INTEGER(INTG), INTENT(OUT) :: adjacentElementNumber !<On exit, the adjacent element number for the specified xi coordinate index in the decomposition elements.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    INTEGER(INTG) :: numberOfXiCoordinates
    TYPE(DecompositionElementType), POINTER :: decompositionElement
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DecompositionElements_ElementAdjacentNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(decompositionElement)
    CALL DecompositionElements_ElementGet(decompositionElements,localElementNumber,decompositionElement,err,error,*999)
    IF(ALLOCATED(decompositionElement%adjacentElements)) THEN
      localError="Adjacent elements is not allocated for local element number "// &
        & TRIM(NumberToVString(localElementNumber,"*",err,error))//" of the decomposition elements."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    numberOfXiCoordinates=(SIZE(decompositionElement%adjacentElements,1)-1)/2
    IF(xiCoordinateIdx<-numberOfXiCoordinates.OR.xiCoordinateIdx>numberOfXiCoordinates) THEN
      localError="The specified xi coordinate index of "//TRIM(NumberToVString(xiCoordinateIdx,"*",err,error))// &
        & " is invalid for element number "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " of the decomposition elements. The element face number should be >= "// &
        & TRIM(NumberToVString(-numberOfXiCoordinates,"*",err,error))//" and <= "// &
        & TRIM(NumberToVString(numberOfXiCoordinates,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(adjacentElementIdx<0.AND. &
      & adjacentElementIdx>decompositionElement%adjacentElements(xiCoordinateIdx)%numberOfAdjacentElements) THEN
      localError="The specified adjacent element index of "//TRIM(NumberToVString(adjacentElementIdx,"*",err,error))// &
        & " is invalid for xi coordinate index "//TRIM(NumberToVString(xiCoordinateIdx,"*",err,error))// &
        & " of local element number "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & ". The adjacent element index should be >= 0 and <= "// &
        & TRIM(NumberToVString(decompositionElement%adjacentElements(xiCoordinateIdx)%numberOfAdjacentElements,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(decompositionElement%adjacentElements(xiCoordinateIdx)%adjacentElements)) THEN
      localError="The adjacent elements array is not allocated for xi coordinate index "// &
        & TRIM(NumberToVString(xiCoordinateIdx,"*",err,error))// &
        & " of local element number "//TRIM(NumberToVString(localElementNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    adjacentElementNumber=decompositionElements%elements(localElementNumber)%adjacentElements(xiCoordinateIdx)% &
      & adjacentElements(adjacentElementIdx)

    EXITS("DecompositionElements_ElementAdjacentNumberGet")
    RETURN
999 ERRORS("DecompositionElements_ElementAdjacentNumberGet",err,error)
    EXITS("DecompositionElements_ElementAdjacentNumberGet")
    RETURN 1

  END SUBROUTINE DecompositionElements_ElementAdjacentNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the boundary node status for a local element number from a decomposition. 
  SUBROUTINE DecompositionElements_ElementBoundaryElementGet(decompositionElements,localElementNumber,boundaryElement,err,error,*)

    !Argument variables
    TYPE(DecompositionElementsType), POINTER :: decompositionElements !<A pointer to the decomposition elements to get the boundary element status for
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The local element number to get the boundary element status for
    LOGICAL, INTENT(OUT) :: boundaryElement !<On exit, the boundary element flag for the local element number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DecompositionElementType), POINTER :: decompositionElement
#endif    

    ENTERS("DecompositionElements_ElementBoundaryElementGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(decompositionElement)
    CALL DecompositionElements_ElementGet(decompositionElements,localElementNumber,decompositionElement,err,error,*999)
#endif

    boundaryElement=decompositionElements%elements(localElementNumber)%boundaryElement

    EXITS("DecompositionElements_ElementBoundaryElementGet")
    RETURN
999 ERRORS("DecompositionElements_ElementBoundaryElementGet",err,error)
    EXITS("DecompositionElements_ElementBoundaryElementGet")
    RETURN 1
    
  END SUBROUTINE DecompositionElements_ElementBoundaryElementGet
  
  !
  !================================================================================================================================
  !

  !>Checks that a user element number exists in a decomposition. 
  SUBROUTINE DecompositionElements_ElementCheckExists(decompositionElements,userElementNumber,elementExists,localElementNumber, &
    & ghostElement,err,error,*)

    !Argument variables
    TYPE(DecompositionElementsType), POINTER :: decompositionElements !<A pointer to the decomposition elements to check the element exists on
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to check if it exists
    LOGICAL, INTENT(OUT) :: elementExists !<On exit, is .TRUE. if the element user number exists in the decomposition topology, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: localElementNumber !<On exit, if the element exists the local number corresponding to the user element number. If the element does not exist then local number will be 0.
    LOGICAL, INTENT(OUT) :: ghostElement !<On exit, is .TRUE. if the local element (if it exists) is a ghost element, .FALSE. if not.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(TreeNodeType), POINTER :: treeNode

    ENTERS("DecompositionElements_ElementCheckExists",err,error,*999)

    elementExists=.FALSE.
    localElementNumber=0
    ghostElement=.FALSE.
#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(decompositionElements)) CALL FlagError("Decomposition elements is not associated.",err,error,*999)
#endif    

    NULLIFY(treeNode)
    CALL Tree_Search(decompositionElements%elementsTree,userElementNumber,treeNode,err,error,*999)
    IF(ASSOCIATED(treeNode)) THEN
      CALL Tree_NodeValueGet(decompositionElements%elementsTree,treeNode,localElementNumber,err,error,*999)
      elementExists=.TRUE.
      ghostElement=localElementNumber>decompositionElements%numberOfElements
    ENDIF

    EXITS("DecompositionElements_ElementCheckExists")
    RETURN
999 ERRORSEXITS("DecompositionElements_ElementCheckExists",err,error)
    RETURN 1

  END SUBROUTINE DecompositionElements_ElementCheckExists

  !
  !================================================================================================================================
  !

  !>Gets an element for a local element number in decomposition elements 
  SUBROUTINE DecompositionElements_ElementGet(decompositionElements,localElementNumber,decompositionElement,err,error,*)

    !Argument variables
    TYPE(DecompositionElementsType), POINTER :: decompositionElements !<A pointer to the decomposition elements to get the decomposition element for
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The local element number to get the decomposition element for
    TYPE(DecompositionElementType), POINTER, INTENT(OUT) :: decompositionElement !<On exit, a pointer to the specified decomposition element. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DecompositionElements_ElementGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(decompositionElement)) CALL FlagError("Decomposition element is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(decompositionElements)) CALL FlagError("Decomposition elements is not associated.",err,error,*999)
    IF(localElementNumber<1.OR.localElementNumber>decompositionElements%totalNumberOfElements) THEN
      localError="The specified local element number of "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " is invalid. The local element number should be >= 1 and <= "// &
        & TRIM(NumberToVString(decompositionElements%totalNumberOfElements,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(decompositionElements%elements)) &
      & CALL FlagError("The elements array is not allocated for the decomposition elements.",err,error,*999)
#endif

    decompositionElement=>decompositionElements%elements(localElementNumber)

#ifdef WITH_POSTCHECKS
    IF(.NOT.ASSOCIATED(decompositionElement)) THEN
      localError="Decomposition element is not associated for local element number "// &
        & TRIM(NumberToVString(localElementNumber,"*",err,error))//" of the decomposition elements."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("DecompositionElements_ElementGet")
    RETURN
999 NULLIFY(decompositionElement)
998 ERRORSEXITS("DecompositionElements_ElementGet",err,error)
    RETURN 1

  END SUBROUTINE DecompositionElements_ElementGet

  !
  !================================================================================================================================
  !

  !>Gets the local face number corresponding to an element face number for a local element number in decomposition elements. 
  SUBROUTINE DecompositionElements_ElementFaceNumberGet(decompositionElements,elementFaceNumber,localElementNumber, &
    & localFaceNumber,err,error,*)

    !Argument variables
    TYPE(DecompositionElementsType), POINTER :: decompositionElements !<A pointer to the decomposition elements to get the face number for
    INTEGER(INTG), INTENT(IN) :: elementFaceNumber !<The element face number to get the face number for.
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The local element number to get the face number for.
    INTEGER(INTG), INTENT(OUT) :: localFaceNumber !<On exit, the local face number for the specified element face number in the decomposition elements.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DecompositionElementType), POINTER :: decompositionElement
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DecompositionElements_ElementFaceNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(decompositionElement)
    CALL DecompositionElements_ElementGet(decompositionElements,localElementNumber,decompositionElement,err,error,*999)
    IF(ALLOCATED(decompositionElement%elementFaces)) THEN
      localError="Element faces is not allocated for local element number "// &
        & TRIM(NumberToVString(localElementNumber,"*",err,error))//" of the decomposition elements."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(elementFaceNumber<0.OR.elementFaceNumber>SIZE(decompositionElement%elementFaces,1)) THEN
      localError="The specified element face number of "//TRIM(NumberToVString(elementFaceNumber,"*",err,error))// &
        & " is invalid for element number "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " of the decomposition elements. The element face number should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(decompositionElement%elementFaces,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    localFaceNumber=decompositionElements%elements(localElementNumber)%elementFaces(elementFaceNumber)

    EXITS("DecompositionElements_ElementFaceNumberGet")
    RETURN
999 ERRORS("DecompositionElements_ElementFaceNumberGet",err,error)
    EXITS("DecompositionElements_ElementFaceNumberGet")
    RETURN 1

  END SUBROUTINE DecompositionElements_ElementFaceNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the global element number for a local element number in decomposition elements. 
  SUBROUTINE DecompositionElements_GlobalElementNumberGet(decompositionElements,localElementNumber,globalElementNumber,err,error,*)

    !Argument variables
    TYPE(DecompositionElementsType), POINTER :: decompositionElements !<A pointer to the decomposition elements to get the global element number for
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The local element number to get the global element number for.
    INTEGER(INTG), INTENT(OUT) :: globalElementNumber !<On exit, the global element number for the specified local element number in the decomposition elements.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DecompositionElementType), POINTER :: decompositionElement
#endif    

    ENTERS("DecompositionElements_GlobalElementNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(decompositionElement)
    CALL DecompositionElements_ElementGet(decompositionElements,localElementNumber,decompositionElement,err,error,*999)
#endif

    globalElementNumber=decompositionElements%elements(localElementNumber)%globalNumber

    EXITS("DecompositionElements_GlobalElementNumberGet")
    RETURN
999 ERRORS("DecompositionElements_GlobalElementNumberGet",err,error)
    EXITS("DecompositionElements_GlobalElementNumberGet")
    RETURN 1

  END SUBROUTINE DecompositionElements_GlobalElementNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the local line number corresponding to an element line number for a local element number in decomposition elements. 
  SUBROUTINE DecompositionElements_ElementLineNumberGet(decompositionElements,elementLineNumber,localElementNumber, &
    & localLineNumber,err,error,*)

    !Argument variables
    TYPE(DecompositionElementsType), POINTER :: decompositionElements !<A pointer to the decomposition elements to get the line number for
    INTEGER(INTG), INTENT(IN) :: elementLineNumber !<The element line number to get the line number for.
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The local element number to get the line number for.
    INTEGER(INTG), INTENT(OUT) :: localLineNumber !<On exit, the local line number for the specified element line number in the decomposition elements.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DecompositionElementType), POINTER :: decompositionElement
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DecompositionElements_ElementLineNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(decompositionElement)
    CALL DecompositionElements_ElementGet(decompositionElements,localElementNumber,decompositionElement,err,error,*999)
    IF(ALLOCATED(decompositionElement%elementLines)) THEN
      localError="Element lines is not allocated for local element number "// &
        & TRIM(NumberToVString(localElementNumber,"*",err,error))//" of the decomposition elements."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(elementLineNumber<0.OR.elementLineNumber>SIZE(decompositionElement%elementLines,1)) THEN
      localError="The specified element line number of "//TRIM(NumberToVString(elementLineNumber,"*",err,error))// &
        & " is invalid for element number "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " of the decomposition elements. The element line number should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(decompositionElement%elementLines,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    localLineNumber=decompositionElements%elements(localElementNumber)%elementLines(elementLineNumber)

    EXITS("DecompositionElements_ElementLineNumberGet")
    RETURN
999 ERRORS("DecompositionElements_ElementLineNumberGet",err,error)
    EXITS("DecompositionElements_ElementLineNumberGet")
    RETURN 1

  END SUBROUTINE DecompositionElements_ElementLineNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the number of adjacent elements corresponding to an xi coordinate index for a local element number in decomposition elements. 
  SUBROUTINE DecompositionElements_ElementNumberAdjacentGet(decompositionElements,xiCoordinateIdx,localElementNumber, &
    & numberOfAdjacentElements,err,error,*)

    !Argument variables
    TYPE(DecompositionElementsType), POINTER :: decompositionElements !<A pointer to the decomposition elements to get the number of adjacent elements for
    INTEGER(INTG), INTENT(IN) :: xiCoordinateIdx !<The xi coordinate index to get the number of adjacent elements for.
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The local element number to get the number of adjacent elements for.
    INTEGER(INTG), INTENT(OUT) :: numberOfAdjacentElements !<On exit, the number of adjacent elements for the specified xi coordinate index in the decomposition elements.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    INTEGER(INTG) :: numberOfXiCoordinates
    TYPE(DecompositionElementType), POINTER :: decompositionElement
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DecompositionElements_ElementNumberAdjacentGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(decompositionElement)
    CALL DecompositionElements_ElementGet(decompositionElements,localElementNumber,decompositionElement,err,error,*999)
    IF(ALLOCATED(decompositionElement%adjacentElements)) THEN
      localError="Adjacent elements is not allocated for local element number "// &
        & TRIM(NumberToVString(localElementNumber,"*",err,error))//" of the decomposition elements."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    numberOfXiCoordinates=(SIZE(decompositionElement%adjacentElements,1)-1)/2
    IF(xiCoordinateIdx<-numberOfXiCoordinates.OR.xiCoordinateIdx>numberOfXiCoordinates) THEN
      localError="The specified xi coordinate index of "//TRIM(NumberToVString(xiCoordinateIdx,"*",err,error))// &
        & " is invalid for element number "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " of the decomposition elements. The xi coordinate index should be >= "// &
        & TRIM(NumberToVString(-numberOfXiCoordinates,"*",err,error))//" and <= "// &
        & TRIM(NumberToVString(numberOfXiCoordinates,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    numberOfAdjacentElements=decompositionElements%elements(localElementNumber)%adjacentElements(xiCoordinateIdx)% &
      & numberOfAdjacentElements

    EXITS("DecompositionElements_ElementNumberAdjacentGet")
    RETURN
999 ERRORS("DecompositionElements_ElementNumberAdjacentGet",err,error)
    EXITS("DecompositionElements_ElementNumberAdjacentGet")
    RETURN 1

  END SUBROUTINE DecompositionElements_ElementNumberAdjacentGet

  !
  !================================================================================================================================
  !

  !>Gets a local element number that corresponds to a user element number from a decomposition. An error will be raised if the user element number does not exist.
  SUBROUTINE DecompositionElements_LocalElementNumberGet(decompositionElements,userElementNumber,localElementNumber,ghostElement, &
    & err,error,*)

    !Argument variables
    TYPE(DecompositionElementsType), POINTER :: decompositionElements !<A pointer to the decomposition elements to get the element on
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to get
    INTEGER(INTG), INTENT(OUT) :: localElementNumber !<On exit, the local number corresponding to the user element number.
    LOGICAL, INTENT(OUT) :: ghostElement !<On exit, is .TRUE. if the local element is a ghost element, .FALSE. if not.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    LOGICAL :: elementExists
#ifdef WITH_POSTCHECKS    
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DecompositionElements_LocalElementNumberGet",err,error,*999)

    CALL DecompositionElements_ElementCheckExists(decompositionElements,userElementNumber,elementExists,localElementNumber, &
      & ghostElement,err,error,*999)

#ifdef WITH_POSTCHECKS    
    IF(.NOT.elementExists) THEN
      decompositionTopology=>decompositionElements%decompositionTopology
      IF(ASSOCIATED(decompositionTopology)) THEN
        decomposition=>decompositionTopology%decomposition
        IF(ASSOCIATED(decomposition)) THEN
          localError="The user element number "//TRIM(NumberToVString(userElementNumber,"*",err,error))// &
            & " does not exist in decomposition number "//TRIM(NumberToVString(decomposition%userNumber,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ELSE
          localError="The user element number "//TRIM(NumberToVString(userElementNumber,"*",err,error))//" does not exist."
        ENDIF
      ELSE
        localError="The user element number "//TRIM(NumberToVString(userElementNumber,"*",err,error))//" does not exist."
      ENDIF
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("DecompositionElements_LocalElementNumberGet")
    RETURN
999 ERRORSEXITS("DecompositionElements_LocalElementNumberGet",err,error)
    RETURN 1

  END SUBROUTINE DecompositionElements_LocalElementNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the number of elements from a decomposition elements. 
  SUBROUTINE DecompositionElements_NumberOfElementsGet(decompositionElements,numberOfElements,err,error,*)

    !Argument variables
    TYPE(DecompositionElementsType), POINTER :: decompositionElements !<A pointer to the decomposition elements to get the number of elements for
    INTEGER(INTG), INTENT(OUT) :: numberOfElements !<On exit, the number of elements for the decomposition elements.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DecompositionElements_NumberOfElementsGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(decompositionElements)) CALL FlagError("Decomposition elements is not associated.",err,error,*999)
#endif

    numberOfElements=decompositionElements%numberOfElements

    EXITS("DecompositionElements_NumberOfElementsGet")
    RETURN
999 ERRORS("DecompositionElements_NumberOfElementsGet",err,error)
    EXITS("DecompositionElements_NumberOfElementsGet")
    RETURN 1

  END SUBROUTINE DecompositionElements_NumberOfElementsGet

  !
  !================================================================================================================================
  !

  !>Gets the total number of elements from a decomposition elements. 
  SUBROUTINE DecompositionElements_TotalNumberOfElementsGet(decompositionElements,totalNumberOfElements,err,error,*)

    !Argument variables
    TYPE(DecompositionElementsType), POINTER :: decompositionElements !<A pointer to the decomposition elements to get the number of elements for
    INTEGER(INTG), INTENT(OUT) :: totalNumberOfElements !<On exit, the total number of elements for the decomposition elements.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DecompositionElements_TotalNumberOfElementsGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(decompositionElements)) CALL FlagError("Decomposition elements is not associated.",err,error,*999)
#endif

    totalNumberOfElements=decompositionElements%totalNumberOfElements

    EXITS("DecompositionElements_TotalNumberOfElementsGet")
    RETURN
999 ERRORS("DecompositionElements_TotalNumberOfElementsGet",err,error)
    EXITS("DecompositionElements_TotalNumberOfElementsGet")
    RETURN 1

  END SUBROUTINE DecompositionElements_TotalNumberOfElementsGet

  !
  !================================================================================================================================
  !

  !>Gets the user element number for a local element number in decomposition elements. 
  SUBROUTINE DecompositionElements_UserElementNumberGet(decompositionElements,localElementNumber,userElementNumber,err,error,*)

    !Argument variables
    TYPE(DecompositionElementsType), POINTER :: decompositionElements !<A pointer to the decomposition elements to get the user element number for
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The local element number to get the user element number for.
    INTEGER(INTG), INTENT(OUT) :: userElementNumber !<On exit, the user element number for the specified local element number in the decomposition elements.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DecompositionElementType), POINTER :: decompositionElement
#endif    

    ENTERS("DecompositionElements_UserElementNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(decompositionElement)
    CALL DecompositionElements_ElementGet(decompositionElements,localElementNumber,decompositionElement,err,error,*999)
#endif

    userElementNumber=decompositionElements%elements(localElementNumber)%userNumber

    EXITS("DecompositionElements_UserElementNumberGet")
    RETURN
999 ERRORS("DecompositionElements_UserElementNumberGet",err,error)
    EXITS("DecompositionElements_UserElementNumberGet")
    RETURN 1

  END SUBROUTINE DecompositionElements_UserElementNumberGet

  !
  !================================================================================================================================
  !

  !>Gets a face for a local face number in decomposition faces 
  SUBROUTINE DecompositionFaces_FaceGet(decompositionFaces,localFaceNumber,decompositionFace,err,error,*)

    !Argument variables
    TYPE(DecompositionFacesType), POINTER :: decompositionFaces !<A pointer to the decomposition faces to get the decomposition face for
    INTEGER(INTG), INTENT(IN) :: localFaceNumber !<The local face number to get the decomposition face for
    TYPE(DecompositionFaceType), POINTER, INTENT(OUT) :: decompositionFace !<On exit, a pointer to the specified decomposition face. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DecompositionFaces_FaceGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(decompositionFace)) CALL FlagError("Decomposition face is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(decompositionFaces)) CALL FlagError("Decomposition faces is not associated.",err,error,*999)
    IF(localFaceNumber<1.OR.localFaceNumber>decompositionFaces%numberOfFaces) THEN
      localError="The specified local face number of "//TRIM(NumberToVString(localFaceNumber,"*",err,error))// &
        & " is invalid. The local face number should be >= 1 and <= "// &
        & TRIM(NumberToVString(decompositionFaces%numberOfFaces,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(decompositionFaces%faces)) &
      & CALL FlagError("The faces array is not allocated for the decomposition faces.",err,error,*999)
#endif

    decompositionFace=>decompositionFaces%faces(localFaceNumber)

#ifdef WITH_POSTCHECKS
    IF(.NOT.ASSOCIATED(decompositionFace)) THEN
      localError="Decomposition face is not associated for local face number "// &
        & TRIM(NumberToVString(localFaceNumber,"*",err,error))//" of the decomposition faces."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("DecompositionFaces_FaceGet")
    RETURN
999 NULLIFY(decompositionFace)
998 ERRORSEXITS("DecompositionFaces_FaceGet",err,error)
    RETURN 1

  END SUBROUTINE DecompositionFaces_FaceGet

  !
  !================================================================================================================================
  !

  !>Gets the boundary face status for a local face number from a decomposition. 
  SUBROUTINE DecompositionFaces_FaceBoundaryFaceGet(decompositionFaces,localFaceNumber,boundaryFace,err,error,*)

    !Argument variables
    TYPE(DecompositionFacesType), POINTER :: decompositionFaces !<A pointer to the decomposition faces to get the boundary face status for
    INTEGER(INTG), INTENT(IN) :: localFaceNumber !<The local face number to get the boundary face status for
    LOGICAL, INTENT(OUT) :: boundaryFace !<On exit, the boundary face flag for the local face number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DecompositionFaceType), POINTER :: decompositionFace
#endif    

    ENTERS("DecompositionFaces_FaceBoundaryFaceGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(decompositionFace)
    CALL DecompositionFaces_FaceGet(decompositionFaces,localFaceNumber,decompositionFace,err,error,*999)
#endif

    boundaryFace=decompositionFaces%faces(localFaceNumber)%boundaryFace

    EXITS("DecompositionFaces_FaceBoundaryFaceGet")
    RETURN
999 ERRORS("DecompositionFaces_FaceBoundaryFaceGet",err,error)
    EXITS("DecompositionFaces_FaceBoundaryFaceGet")
    RETURN 1
    
  END SUBROUTINE DecompositionFaces_FaceBoundaryFaceGet
  
  !
  !================================================================================================================================
  !

  !>Gets the xi normal direction for a local face number from a decomposition. 
  SUBROUTINE DecompositionFaces_FaceXiNormalDirectionGet(decompositionFaces,localFaceNumber,faceXiNormalDirection,err,error,*)

    !Argument variables
    TYPE(DecompositionFacesType), POINTER :: decompositionFaces !<A pointer to the decomposition faces to get the xi normal direction for
    INTEGER(INTG), INTENT(IN) :: localFaceNumber !<The local face number to get the xi normal direction for
    INTEGER(INTG), INTENT(OUT) :: faceXiNormalDirection !<On exit, the face xi normal direction for the local face number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DecompositionFaceType), POINTER :: decompositionFace
#endif    

    ENTERS("DecompositionFaces_FaceXiNormalDirectionGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(decompositionFace)
    CALL DecompositionFaces_FaceGet(decompositionFaces,localFaceNumber,decompositionFace,err,error,*999)
#endif

    faceXiNormalDirection=decompositionFaces%faces(localFaceNumber)%xiNormalDirection

    EXITS("DecompositionFaces_FaceXiNormalDirectionGet")
    RETURN
999 ERRORS("DecompositionFaces_FaceXiNormalDirectionGet",err,error)
    EXITS("DecompositionFaces_FaceXiNormalDirectionGet")
    RETURN 1
    
  END SUBROUTINE DecompositionFaces_FaceXiNormalDirectionGet
  
  !
  !================================================================================================================================
  !

  !>Gets the number of faces from a decomposition. 
  SUBROUTINE DecompositionFaces_NumberOfFacesGet(decompositionFaces,numberOfFaces,err,error,*)

    !Argument variables
    TYPE(DecompositionFacesType), POINTER :: decompositionFaces !<A pointer to the decomposition faces to get the number of faces for
    INTEGER(INTG), INTENT(OUT) :: numberOfFaces !<On exit, the number of faces in the decomposition faces.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DecompositionFaces_NumberOfFacesGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(decompositionFaces)) CALL FlagError("Decomposition faces is not associated.",err,error,*999)
#endif

    numberOfFaces=decompositionFaces%numberOfFaces

    EXITS("DecompositionFaces_NumberOfFacesGet")
    RETURN
999 ERRORS("DecompositionFaces_NumberOfFacesGet",err,error)
    EXITS("DecompositionFaces_NumberOfFacesGet")
    RETURN 1
    
  END SUBROUTINE DecompositionFaces_NumberOfFacesGet
  
  !
  !================================================================================================================================
  !

  !>Gets a line for a local line number in decomposition lines 
  SUBROUTINE DecompositionLines_LineGet(decompositionLines,localLineNumber,decompositionLine,err,error,*)

    !Argument variables
    TYPE(DecompositionLinesType), POINTER :: decompositionLines !<A pointer to the decomposition lines to get the decomposition line for
    INTEGER(INTG), INTENT(IN) :: localLineNumber !<The local line number to get the decomposition line for
    TYPE(DecompositionLineType), POINTER, INTENT(OUT) :: decompositionLine !<On exit, a pointer to the specified decomposition line. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DecompositionLines_LineGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(decompositionLine)) CALL FlagError("Decomposition line is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(decompositionLines)) CALL FlagError("Decomposition lines is not associated.",err,error,*999)
    IF(localLineNumber<1.OR.localLineNumber>decompositionLines%numberOfLines) THEN
      localError="The specified local line number of "//TRIM(NumberToVString(localLineNumber,"*",err,error))// &
        & " is invalid. The local line number should be >= 1 and <= "// &
        & TRIM(NumberToVString(decompositionLines%numberOfLines,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(decompositionLines%lines)) &
      & CALL FlagError("The lines array is not allocated for the decomposition lines.",err,error,*999)
#endif

    decompositionLine=>decompositionLines%lines(localLineNumber)

#ifdef WITH_POSTCHECKS
    IF(.NOT.ASSOCIATED(decompositionLine)) THEN
      localError="Decomposition line is not associated for local line number "// &
        & TRIM(NumberToVString(localLineNumber,"*",err,error))//" of the decomposition lines."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("DecompositionLines_LineGet")
    RETURN
999 NULLIFY(decompositionLine)
998 ERRORSEXITS("DecompositionLines_LineGet",err,error)
    RETURN 1

  END SUBROUTINE DecompositionLines_LineGet

  !
  !================================================================================================================================
  !

  !>Gets the boundary line status for a local line number from a decomposition. 
  SUBROUTINE DecompositionLines_LineBoundaryLineGet(decompositionLines,localLineNumber,boundaryLine,err,error,*)

    !Argument variables
    TYPE(DecompositionLinesType), POINTER :: decompositionLines !<A pointer to the decomposition lines to get the boundary line status for
    INTEGER(INTG), INTENT(IN) :: localLineNumber !<The local line number to get the boundary line status for
    LOGICAL, INTENT(OUT) :: boundaryLine !<On exit, the boundary line flag for the local line number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DecompositionLineType), POINTER :: decompositionLine
#endif    

    ENTERS("DecompositionLines_LineBoundaryLineGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(decompositionLine)
    CALL DecompositionLines_LineGet(decompositionLines,localLineNumber,decompositionLine,err,error,*999)
#endif

    boundaryLine=decompositionLines%lines(localLineNumber)%boundaryLine

    EXITS("DecompositionLines_LineBoundaryLineGet")
    RETURN
999 ERRORS("DecompositionLines_LineBoundaryLineGet",err,error)
    EXITS("DecompositionLines_LineBoundaryLineGet")
    RETURN 1
    
  END SUBROUTINE DecompositionLines_LineBoundaryLineGet
  
  !
  !================================================================================================================================
  !

  !>Gets the xi direction for a local line number from a decomposition. 
  SUBROUTINE DecompositionLines_LineXiDirectionGet(decompositionLines,localLineNumber,lineXiDirection,err,error,*)

    !Argument variables
    TYPE(DecompositionLinesType), POINTER :: decompositionLines !<A pointer to the decomposition lines to get the xi direction for
    INTEGER(INTG), INTENT(IN) :: localLineNumber !<The local line number to get the xi direction for
    INTEGER(INTG), INTENT(OUT) :: lineXiDirection !<On exit, the line xi direction for the local line number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DecompositionLineType), POINTER :: decompositionLine
#endif    

    ENTERS("DecompositionLines_LineXiDirectionGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(decompositionLine)
    CALL DecompositionLines_LineGet(decompositionLines,localLineNumber,decompositionLine,err,error,*999)
#endif

    lineXiDirection=decompositionLines%lines(localLineNumber)%xiDirection

    EXITS("DecompositionLines_LineXiDirectionGet")
    RETURN
999 ERRORS("DecompositionLines_LineXiDirectionGet",err,error)
    EXITS("DecompositionLines_LineXiDirectionGet")
    RETURN 1
    
  END SUBROUTINE DecompositionLines_LineXiDirectionGet
  
  !
  !================================================================================================================================
  !

  !>Gets the number of lines from a decomposition. 
  SUBROUTINE DecompositionLines_NumberOfLinesGet(decompositionLines,numberOfLines,err,error,*)

    !Argument variables
    TYPE(DecompositionLinesType), POINTER :: decompositionLines !<A pointer to the decomposition lines to get the number of lines for
    INTEGER(INTG), INTENT(OUT) :: numberOfLines !<On exit, the number of lines in the decomposition lines.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DecompositionLines_NumberOfLinesGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(decompositionLines)) CALL FlagError("Decomposition lines is not associated.",err,error,*999)
#endif

    numberOfLines=decompositionLines%numberOfLines

    EXITS("DecompositionLines_NumberOfLinesGet")
    RETURN
999 ERRORS("DecompositionLines_NumberOfLinesGet",err,error)
    EXITS("DecompositionLines_NumberOfLinesGet")
    RETURN 1
    
  END SUBROUTINE DecompositionLines_NumberOfLinesGet
  
  !
  !================================================================================================================================
  !

  !>Gets elements from a decomposition topology.
  SUBROUTINE DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to get the elements for.
     TYPE(DecompositionElementsType), POINTER :: decompositionElements !<On exit, a pointer to the decomposition topology elements. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DecompositionTopology_DecompositionElementsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(decompositionElements)) CALL FlagError("Decomposition elements is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(decompositionTopology)) CALL FlagError("Decomposition topology is not associated.",err,error,*998)
#endif
    
    !Get the decomposition elements
    decompositionElements=>decompositionTopology%elements

#ifdef WITH_POSTCHECKS    
    !Check decomposition elements is associated.
    IF(.NOT.ASSOCIATED(decompositionElements)) CALL FlagError("Decomposition topology elements is not associated.",err,error,*999)
#endif    
    
    EXITS("DecompositionTopology_DecompositionElementsGet")
    RETURN
999 NULLIFY(decompositionElements)
998 ERRORS("DecompositionTopology_DecompositionElementsGet",err,error)
    EXITS("DecompositionTopology_DecompositionElementsGet")
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_DecompositionElementsGet

  !
  !================================================================================================================================
  !

  !>Gets data points from a decomposition topology.
  SUBROUTINE DecompositionTopology_DecompositionDataPointsGet(decompositionTopology,decompositionDataPoints,err,error,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to get the data points for.
    TYPE(DecompositionDataPointsType), POINTER :: decompositionDataPoints !<On exit, a pointer to the decomposition topology data points. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DecompositionTopology_DecompositionDataPointsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(decompositionDataPoints)) CALL FlagError("Decomposition data points is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(decompositionTopology)) CALL FlagError("Decomposition topology is not associated.",err,error,*999)
#endif    

    !Get the decomposition data points
    decompositionDataPoints=>decompositionTopology%dataPoints

#ifdef WITH_POSTCHECKS    
    !Check decomposition data points is associated.
    IF(.NOT.ASSOCIATED(decompositionDataPoints)) &
      & CALL FlagError("Decomposition topology data points is not associated.",err,error,*999)
#endif    
    
    EXITS("DecompositionTopology_DecompositionDataPointsGet")
    RETURN
999 NULLIFY(decompositionDataPoints)
998 ERRORS("DecompositionTopology_DecompositionDataPointsGet",err,error)
    EXITS("DecompositionTopology_DecompositionDataPointsGet")
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_DecompositionDataPointsGet

  !
  !================================================================================================================================
  !

  !>Gets the decomposition from a decomposition topology.
  SUBROUTINE DecompositionTopology_DecompositionGet(decompositionTopology,decomposition,err,error,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to get the decomposition for.
    TYPE(DecompositionType), POINTER :: decomposition !<On exit, a pointer to the decomposition topology decomposition. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DecompositionTopology_DecompositionGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(decomposition)) CALL FlagError("Decomposition is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(decompositionTopology)) CALL FlagError("Decomposition topology is not associated.",err,error,*999)
#endif    

    !Get the decomposition data points
    decomposition=>decompositionTopology%decomposition

#ifdef WITH_POSTCHECKS    
    !Check decomposition is associated.
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition topology decomposition is not associated.",err,error,*999)
#endif    
    
    EXITS("DecompositionTopology_DecompositionGet")
    RETURN
999 NULLIFY(decomposition)
998 ERRORS("DecompositionTopology_DecompositionGet",err,error)
    EXITS("DecompositionTopology_DecompositionGet")
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_DecompositionGet

  !
  !================================================================================================================================
  !

  !>Gets faces from a decomposition topology.
  SUBROUTINE DecompositionTopology_DecompositionFacesGet(decompositionTopology,decompositionFaces,err,error,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to get the faces for.
    TYPE(DecompositionFacesType), POINTER :: decompositionFaces !<On exit, a pointer to the decomposition topology faces. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DecompositionTopology_DecompositionFacesGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(decompositionFaces)) CALL FlagError("Decomposition faces is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(decompositionTopology)) CALL FlagError("Decomposition is not associated.",err,error,*999)
#endif    

    !Get the decomposition faces
    decompositionFaces=>decompositionTopology%faces

#ifdef WITH_POSTCHECKS    
    !Check decomposition faces is associated.
    IF(.NOT.ASSOCIATED(decompositionFaces)) CALL FlagError("Decomposition topology faces is not associated.",err,error,*999)
#endif    
    
    EXITS("DecompositionTopology_DecompositionFacesGet")
    RETURN
999 NULLIFY(decompositionFaces)
998 ERRORSEXITS("DecompositionTopology_DecompositionFacesGet",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_DecompositionFacesGet

  !
  !================================================================================================================================
  !

  !>Gets lines from a decomposition topology.
  SUBROUTINE DecompositionTopology_DecompositionLinesGet(decompositionTopology,decompositionLines,err,error,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to get the lines for.
    TYPE(DecompositionLinesType), POINTER :: decompositionLines !<On exit, a pointer to the decomposition topology lines. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DecompositionTopology_DecompositionLinesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(.NOT.ASSOCIATED(decompositionTopology)) CALL FlagError("Decomposition is not associated.",err,error,*999)
    IF(ASSOCIATED(decompositionLines)) CALL FlagError("Decomposition lines is already associated.",err,error,*999)
#endif    
 
    !Get the decomposition lines
    decompositionLines=>decompositionTopology%lines

#ifdef WITH_POSTCHECKS    
    !Check decomposition lines is associated.
    IF(.NOT.ASSOCIATED(decompositionLines)) CALL FlagError("Decomposition topology lines is not associated.",err,error,*999)
#endif    
    
    EXITS("DecompositionTopology_DecompositionLinesGet")
    RETURN
999 ERRORSEXITS("DecompositionTopology_DecompositionLinesGet",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_DecompositionLinesGet
  
  !  
  !================================================================================================================================
  !

  !>Get the basis for a domain element.
  SUBROUTINE DomainElement_BasisGet(domainElement,basis,err,error,*)

    !Argument variables
    TYPE(DomainElementType), POINTER :: domainElement !<A pointer to the domain element to get the basis for
    TYPE(BasisType), POINTER :: basis !<On return, a pointer to the basis for the domain element. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainElement_BasisGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainElement)) CALL FlagError("Domain element is not associated.",err,error,*999)
#endif    
      
    basis=>domainElement%basis

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated for the domain element.",err,error,*999)
#endif    
    
    EXITS("DomainElement_BasisGet")
    RETURN
999 NULLIFY(basis)
998 ERRORSEXITS("DomainElement_BasisGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainElement_BasisGet

  !  
  !================================================================================================================================
  !

  !>Get the global derivative index for a local node and derivative index in a domain element.
  SUBROUTINE DomainElement_ElementDerivativeGet(domainElement,localDerivativeIdx,localNodeIdx,globalDerivativeIdx,err,error,*)

    !Argument variables
    TYPE(DomainElementType), POINTER :: domainElement !<A pointer to the domain element to get the global derivative index for
    INTEGER(INTG), INTENT(IN) :: localDerivativeIdx !<The local derivative index to get the element global derivative index for
    INTEGER(INTG), INTENT(IN) :: localNodeIdx !<The local node index to get the global derivative index for
    INTEGER(INTG), INTENT(OUT) :: globalDerivativeIdx !<On return, the global derivative index for the local node and derivative index of the domain element.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainElement_ElementDerivativeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainElement)) CALL FlagError("Domain element is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(domainElement%elementDerivatives)) &
      & CALL FlagError("The element derivatives array is not allocated for the domain element.",err,error,*999)
    IF(localDerivativeIdx<1.OR.localDerivativeIdx>SIZE(domainElement%elementDerivatives,1)) THEN
      localError="The specified local derivative index of "//TRIM(NumberToVString(localDerivativeIdx,"*",err,error))// &
        & " is invalid. The local derivative index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(domainElement%elementDerivatives,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localNodeIdx<1.OR.localNodeIdx>SIZE(domainElement%elementDerivatives,2)) THEN
      localError="The specified local node index of "//TRIM(NumberToVString(localNodeIdx,"*",err,error))// &
        & " is invalid. The local node index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(domainElement%elementDerivatives,2),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    globalDerivativeIdx=domainElement%elementDerivatives(localDerivativeIdx,localNodeIdx)
    
    EXITS("DomainElement_ElementDerivativeGet")
    RETURN
999 ERRORSEXITS("DomainElement_ElementDerivativeGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainElement_ElementDerivativeGet

  !  
  !================================================================================================================================
  !

  !>Get the local element node number for a local node in a domain element.
  SUBROUTINE DomainElement_ElementNodeGet(domainElement,localNodeIdx,localNodeNumber,err,error,*)

    !Argument variables
    TYPE(DomainElementType), POINTER :: domainElement !<A pointer to the domain element to get the local element node for
    INTEGER(INTG), INTENT(IN) :: localNodeIdx !<The local node index to get the local element node for
    INTEGER(INTG), INTENT(OUT) :: localNodeNumber !<On return, the local node number for the local node index of the domain element.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainElement_ElementNodeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainElement)) CALL FlagError("Domain element is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(domainElement%elementNodes)) &
      & CALL FlagError("The element nodes array is not allocated for the domain element.",err,error,*999)
    IF(localNodeIdx<1.OR.localNodeIdx>SIZE(domainElement%elementNodes,1)) THEN
      localError="The specified local node index of "//TRIM(NumberToVString(localNodeIdx,"*",err,error))// &
        & " is invalid. The local node index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(domainElement%elementNodes,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    localNodeNumber=domainElement%elementNodes(localNodeIdx)
    
    EXITS("DomainElement_ElementNodeGet")
    RETURN
999 ERRORSEXITS("DomainElement_ElementNodeGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainElement_ElementNodeGet

  !  
  !================================================================================================================================
  !

  !>Get the element version number for a local node and derivative index in a domain element.
  SUBROUTINE DomainElement_ElementVersionGet(domainElement,localDerivativeIdx,localNodeIdx,versionNumber,err,error,*)

    !Argument variables
    TYPE(DomainElementType), POINTER :: domainElement !<A pointer to the domain element to get the element version number for
    INTEGER(INTG), INTENT(IN) :: localDerivativeIdx !<The local derivative index to get the element version number for
    INTEGER(INTG), INTENT(IN) :: localNodeIdx !<The local node index to get the element version number for
    INTEGER(INTG), INTENT(OUT) :: versionNumber !<On return, the version number for the local node and derivative index of the domain element.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainElement_ElementVersionGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainElement)) CALL FlagError("Domain element is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(domainElement%elementVersions)) &
      & CALL FlagError("The element versions array is not allocated for the domain element.",err,error,*999)
    IF(localDerivativeIdx<1.OR.localDerivativeIdx>SIZE(domainElement%elementVersions,1)) THEN
      localError="The specified local derivative index of "//TRIM(NumberToVString(localDerivativeIdx,"*",err,error))// &
        & " is invalid. The local derivative index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(domainElement%elementVersions,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localNodeIdx<1.OR.localNodeIdx>SIZE(domainElement%elementVersions,2)) THEN
      localError="The specified local node index of "//TRIM(NumberToVString(localNodeIdx,"*",err,error))// &
        & " is invalid. The local node index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(domainElement%elementVersions,2),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    versionNumber=domainElement%elementVersions(localDerivativeIdx,localNodeIdx)
    
    EXITS("DomainElement_ElementVersionGet")
    RETURN
999 ERRORSEXITS("DomainElement_ElementVersionGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainElement_ElementVersionGet
  
  !  
  !================================================================================================================================
  !

  !>Get the basis for an element in the domain elements identified by its local number
  SUBROUTINE DomainElements_ElementBasisGet(domainElements,localElementNumber,basis,err,error,*)

    !Argument variables
    TYPE(DomainElementsType), POINTER :: domainElements !<A pointer to the domain elements to get the element basis for
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The element local number to get the basis for
    TYPE(BasisType), POINTER, INTENT(OUT) :: basis !<On return, a pointer to the basis for the element. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DomainElementType), POINTER :: domainElement
#endif
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainElements_ElementBasisGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated.",err,error,*998)
    NULLIFY(domainElement)
    CALL DomainElements_DomainElementGet(domainElements,localElementNumber,domainElement,err,error,*999)
#endif    
      
    basis=>domainElements%elements(localElementNumber)%basis

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(basis)) THEN
      localError="The basis for local element number "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("DomainElements_ElementBasisGet")
    RETURN
999 NULLIFY(basis)
998 ERRORSEXITS("DomainElements_ElementBasisGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainElements_ElementBasisGet

  !  
  !================================================================================================================================
  !

  !>Get a domain element for a local element number in the domain elements
  SUBROUTINE DomainElements_DomainElementGet(domainElements,localElementNumber,domainElement,err,error,*)

    !Argument variables
    TYPE(DomainElementsType), POINTER :: domainElements !<A pointer to the domain elements to get the domain element for
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The element local number to get the domain element for
    TYPE(DomainElementType), POINTER, INTENT(OUT) :: domainElement !<On return, a pointer to the domain element for the element. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainElements_DomainElementGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(domainElement)) CALL FlagError("Domain element is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainElements)) CALL FlagError("Domain elements is not associated.",err,error,*999)
    IF(localElementNumber<=0.OR.localElementNumber>domainElements%totalNumberOfElements) THEN
      localError="The local element number of "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " is invalid. The local element number must be >= 1 and <= "// &
        & TRIM(NumberToVString(domainElements%totalNumberOfElements,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainElements%elements)) CALL FlagError("Domain elements elements is not allocated.",err,error,*999)
#endif    
      
    domainElement=>domainElements%elements(localElementNumber)

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(domainElement)) THEN
      localError="The domain element is not associated for local element number "// &
        & TRIM(NumberToVString(localElementNumber,"*",err,error))//" of the domain elements."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("DomainElements_DomainElementGet")
    RETURN
999 NULLIFY(domainElement)
998 ERRORSEXITS("DomainElements_DomainElementGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainElements_DomainElementGet

  !  
  !================================================================================================================================
  !

  !>Get the element global derivative index for the local node and derivative index of an element in the domain.
  SUBROUTINE DomainElements_ElementDerivativeGet(domainElements,localDerivativeIdx,localNodeIdx,localElementNumber, &
    & globalDerivativeIndex,err,error,*)

    !Argument variables
    TYPE(DomainElementsType), POINTER :: domainElements !<A pointer to the domain elements to get the element global derivative index for
    INTEGER(INTG), INTENT(IN) :: localDerivativeIdx !<The local derivative index to get the element global derivative index for.
    INTEGER(INTG), INTENT(IN) :: localNodeIdx !<The local node index to get the element global derivative index for.
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The element local number to get the elment global derivative index for
    INTEGER(INTG), INTENT(OUT) :: globalDerivativeIndex !<On return, the global derivative index for the local node and derivative index for the element.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DomainElementType), POINTER :: domainElement
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainElements_ElementDerivativeGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(domainElement)
    CALL DomainElements_DomainElementGet(domainElements,localElementNumber,domainElement,err,error,*999)
    IF(.NOT.ALLOCATED(domainElement%elementDerivatives)) THEN
      localError="The element derivatives array is not allocated for element number "// &
        & TRIM(NumberToVString(localElementNumber,"*",err,error))//" of the domain elements."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localDerivativeIdx<1.OR.localDerivativeIdx>SIZE(domainElement%elementDerivatives,1)) THEN
      localError="The specified local derivative index of "//TRIM(NumberToVString(localDerivativeIdx,"*",err,error))// &
        & " is invalid for element number "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & ". The local node index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(domainElement%elementDerivatives,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localNodeIdx<1.OR.localNodeIdx>SIZE(domainElement%elementDerivatives,2)) THEN
      localError="The specified local node index of "//TRIM(NumberToVString(localNodeIdx,"*",err,error))// &
        & " is invalid for element number "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & ". The local node index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(domainElement%elementDerivatives,2),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    globalDerivativeIndex=domainElements%elements(localElementNumber)%elementDerivatives(localDerivativeIdx,localNodeIdx)
    
    EXITS("DomainElements_ElementDerivativeGet")
    RETURN
999 ERRORSEXITS("DomainElements_ElementDerivativeGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainElements_ElementDerivativeGet
  
  !  
  !================================================================================================================================
  !

  !>Get the local node number for the local node index of an element in the domain.
  SUBROUTINE DomainElements_ElementNodeGet(domainElements,localNodeIdx,localElementNumber,localNodeNumber,err,error,*)

    !Argument variables
    TYPE(DomainElementsType), POINTER :: domainElements !<A pointer to the domain elements to get the element basis for
    INTEGER(INTG), INTENT(IN) :: localNodeIdx !<The local node index to get the element local node for.
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The element local number to get the local node for
    INTEGER(INTG), INTENT(OUT) :: localNodeNumber !<On return, the local node number for the local node index for the element.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DomainElementType), POINTER :: domainElement
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainElements_ElementNodeGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(domainElement)
    CALL DomainElements_DomainElementGet(domainElements,localElementNumber,domainElement,err,error,*999)
    IF(.NOT.ALLOCATED(domainElement%elementNodes)) THEN
      localError="The element nodes array is not allocated for element number "// &
        & TRIM(NumberToVString(localElementNumber,"*",err,error))//" of the domain elements."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localNodeIdx<1.OR.localNodeIdx>SIZE(domainElement%elementNodes,1)) THEN
      localError="The specified local node index of "//TRIM(NumberToVString(localNodeIdx,"*",err,error))// &
        & " is invalid for element number "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & ". The local node index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(domainElement%elementNodes,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    localNodeNumber=domainElements%elements(localElementNumber)%elementNodes(localNodeIdx)
    
    EXITS("DomainElements_ElementNodeGet")
    RETURN
999 ERRORSEXITS("DomainElements_ElementNodeGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainElements_ElementNodeGet
  
  !  
  !================================================================================================================================
  !

  !>Get the element version number for the local node and derivative index of an element in the domain.
  SUBROUTINE DomainElements_ElementVersionGet(domainElements,localDerivativeIdx,localNodeIdx,localElementNumber, &
    & versionNumber,err,error,*)

    !Argument variables
    TYPE(DomainElementsType), POINTER :: domainElements !<A pointer to the domain elements to get the element version number for
    INTEGER(INTG), INTENT(IN) :: localDerivativeIdx !<The local derivative index to get the element version number for.
    INTEGER(INTG), INTENT(IN) :: localNodeIdx !<The local node index to get the element version number for.
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The element local number to get the elment version number for
    INTEGER(INTG), INTENT(OUT) :: versionNumber !<On return, the version number for the local node and derivative index for the element.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DomainElementType), POINTER :: domainElement
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainElements_ElementVersionGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(domainElement)
    CALL DomainElements_DomainElementGet(domainElements,localElementNumber,domainElement,err,error,*999)
    IF(.NOT.ALLOCATED(domainElement%elementVersions)) THEN
      localError="The element versions array is not allocated for element number "// &
        & TRIM(NumberToVString(localElementNumber,"*",err,error))//" of the domain elements."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localDerivativeIdx<1.OR.localDerivativeIdx>SIZE(domainElement%elementVersions,1)) THEN
      localError="The specified local derivative index of "//TRIM(NumberToVString(localDerivativeIdx,"*",err,error))// &
        & " is invalid for element number "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & ". The local node index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(domainElement%elementVersions,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localNodeIdx<1.OR.localNodeIdx>SIZE(domainElement%elementVersions,2)) THEN
      localError="The specified local node index of "//TRIM(NumberToVString(localNodeIdx,"*",err,error))// &
        & " is invalid for element number "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & ". The local node index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(domainElement%elementVersions,2),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    versionNumber=domainElements%elements(localElementNumber)%elementVersions(localDerivativeIdx,localNodeIdx)
    
    EXITS("DomainElements_ElementVersionGet")
    RETURN
999 ERRORSEXITS("DomainElements_ElementVersionGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainElements_ElementVersionGet
  
  !
  !================================================================================================================================
  !

  !>Gets the maximum number of element parameter for element in a domain. 
  SUBROUTINE DomainElements_MaxElementParametersGet(domainElements,maxElementParameters,err,error,*)

    !Argument variables
    TYPE(DomainElementsType), POINTER :: domainElements !<A pointer to the domain elements to get the maximum number of element parameters for
    INTEGER(INTG), INTENT(OUT) :: maxElementParameters !<On exit, the maximum number of element parameters for the domain elements.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainElements_MaxElementParametersGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(domainElements)) CALL FlagError("Domain elements is not associated.",err,error,*999)
#endif

    maxElementParameters=domainElements%maximumNumberOfElementParameters

    EXITS("DomainElements_MaxElementParametersGet")
    RETURN
999 ERRORSEXITS("DomainElements_MaxElementParametersGet",err,error)
    RETURN 1

  END SUBROUTINE DomainElements_MaxElementParametersGet

  !
  !================================================================================================================================
  !

  !>Gets the number of elements from a domain. 
  SUBROUTINE DomainElements_NumberOfElementsGet(domainElements,numberOfElements,err,error,*)

    !Argument variables
    TYPE(DomainElementsType), POINTER :: domainElements !<A pointer to the domain elements to get the number of elements for
    INTEGER(INTG), INTENT(OUT) :: numberOfElements !<On exit, the number of elements for the domain elements.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainElements_NumberOfElementsGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(domainElements)) CALL FlagError("Domain elements is not associated.",err,error,*999)
#endif

    numberOfElements=domainElements%numberOfElements

    EXITS("DomainElements_NumberOfElementssGet")
    RETURN
999 ERRORSEXITS("DomainElements_NumberOfElementsGet",err,error)
    RETURN 1

  END SUBROUTINE DomainElements_NumberOfElementsGet

  !
  !================================================================================================================================
  !

  !>Gets the total number of elements from a domain. 
  SUBROUTINE DomainElements_TotalNumberOfElementsGet(domainElements,totalNumberOfElements,err,error,*)

    !Argument variables
    TYPE(DomainElementsType), POINTER :: domainElements !<A pointer to the domain elements to get the total number of elements for
    INTEGER(INTG), INTENT(OUT) :: totalNumberOfElements !<On exit, the total number of elements for the domain elements.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainElements_TotalNumberOfElementsGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(domainElements)) CALL FlagError("Domain elements is not associated.",err,error,*999)
#endif

    totalNumberOfElements=domainElements%totalNumberOfElements

    EXITS("DomainElements_TotalNumberOfElementsGet")
    RETURN
999 ERRORSEXITS("DomainElements_TotalNumberOfElementsGet",err,error)
    RETURN 1

  END SUBROUTINE DomainElements_TotalNumberOfElementsGet

  !
  !================================================================================================================================
  !

  !>Gets the global derivative index for a derivative index of a local node index in the domain faces identified by its local number
  SUBROUTINE DomainFaces_DerivativeGlobalIndexGet(domainFaces,derivativeIdx,localNodeIdx,faceNumber,globalDerivativeIndex, &
    & err,error,*)

    !Argument variables
    TYPE(DomainFacesType), POINTER :: domainFaces !<A pointer to the domain faces to get the global derivative index for
    INTEGER(INTG), INTENT(IN) :: derivativeIdx !<The local derivative index to get the global derivative index for.
    INTEGER(INTG), INTENT(IN) :: localNodeIdx !<The local node index to get the global derivative index.
    INTEGER(INTG), INTENT(IN) :: faceNumber !<The face number to get the node number for.
    INTEGER(INTG), INTENT(OUT) :: globalDerivativeIndex !<On return, the global derivative index for the derivative index of the local node index in the specified face.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainFaces_DerivativeGlobalIndexGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainFaces)) CALL FlagError("Domain faces is not associated.",err,error,*999)
    IF(faceNumber<1.OR.faceNumber>domainFaces%numberOfFaces) THEN
      localError="The face number of "//TRIM(NumberToVString(faceNumber,"*",err,error))// &
        & " is invalid. The face number must be >= 1 and <= "// &
        & TRIM(NumberToVString(domainFaces%numberOfFaces,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainFaces%faces)) CALL FlagError("Domain faces faxes is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(domainFaces%faces(faceNumber)%derivativesInFace)) THEN
      localError="The derivatives in face is not allocated for face number "//TRIM(NumberToVString(faceNumber,"*",err,error))// &
        & " in the domain faces."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(derivativeIdx<1.OR.derivativeIdx>SIZE(domainFaces%faces(faceNumber)%derivativesInFace,2)) THEN
      localError="The specified derivative index of "//TRIM(NumberToVString(derivativeIdx,"*",err,error))// &
        & " is invalid for face number "//TRIM(NumberToVString(faceNumber,"*",err,error))// &
        & " in the domain faces. The derivative index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(domainFaces%faces(faceNumber)%derivativesInFace,2),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localNodeIdx<1.OR.localNodeIdx>SIZE(domainFaces%faces(faceNumber)%derivativesInFace,3)) THEN
      localError="The specified local node index of "//TRIM(NumberToVString(localNodeIdx,"*",err,error))// &
        & " is invalid for face number "//TRIM(NumberToVString(faceNumber,"*",err,error))// &
        & " in the domain faces. The local node index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(domainFaces%faces(faceNumber)%derivativesInFace,3),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    globalDerivativeIndex=domainFaces%faces(faceNumber)%derivativesInFace(1,derivativeIdx,localNodeIdx)
    
    EXITS("DomainFaces_DerivativeGlobalIndexGet")
    RETURN
999 ERRORSEXITS("DomainFaces_DerivativeGlobalIndexGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainFaces_DerivativeGlobalIndexGet
  
  !
  !================================================================================================================================
  !

  !>Gets the version number for a derivative index of a local node index in the domain faces identified by its local number
  SUBROUTINE DomainFaces_DerivativeVersionNumberGet(domainFaces,derivativeIdx,localNodeIdx,faceNumber,versionNumber, &
    & err,error,*)

    !Argument variables
    TYPE(DomainFacesType), POINTER :: domainFaces !<A pointer to the domain faces to get the version number for
    INTEGER(INTG), INTENT(IN) :: derivativeIdx !<The local derivative index to get the version number for.
    INTEGER(INTG), INTENT(IN) :: localNodeIdx !<The local node index to get the version number for.
    INTEGER(INTG), INTENT(IN) :: faceNumber !<The face number to get the version number for.
    INTEGER(INTG), INTENT(OUT) :: versionNumber !<On return, the version number for the derivative index of the local node index in the specified face.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainFaces_DerivativeVersionNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainFaces)) CALL FlagError("Domain faces is not associated.",err,error,*999)
    IF(faceNumber<1.OR.faceNumber>domainFaces%numberOfFaces) THEN
      localError="The face number of "//TRIM(NumberToVString(faceNumber,"*",err,error))// &
        & " is invalid. The face number must be >= 1 and <= "// &
        & TRIM(NumberToVString(domainFaces%numberOfFaces,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainFaces%faces)) CALL FlagError("Domain faces faces is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(domainFaces%faces(faceNumber)%derivativesInFace)) THEN
      localError="The derivatives in face is not allocated for face number "//TRIM(NumberToVString(faceNumber,"*",err,error))// &
        & " in the domain faces."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(derivativeIdx<1.OR.derivativeIdx>SIZE(domainFaces%faces(faceNumber)%derivativesInFace,2)) THEN
      localError="The specified derivative index of "//TRIM(NumberToVString(derivativeIdx,"*",err,error))// &
        & " is invalid for face number "//TRIM(NumberToVString(faceNumber,"*",err,error))// &
        & " in the domain faces. The derivative index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(domainFaces%faces(faceNumber)%derivativesInFace,2),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localNodeIdx<1.OR.localNodeIdx>SIZE(domainFaces%faces(faceNumber)%derivativesInFace,3)) THEN
      localError="The specified local node index of "//TRIM(NumberToVString(localNodeIdx,"*",err,error))// &
        & " is invalid for face number "//TRIM(NumberToVString(faceNumber,"*",err,error))// &
        & " in the domain faces. The local node index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(domainFaces%faces(faceNumber)%derivativesInFace,3),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    versionNumber=domainFaces%faces(faceNumber)%derivativesInFace(2,derivativeIdx,localNodeIdx)
    
    EXITS("DomainFaces_DerivativeVersionNumberGet")
    RETURN
999 ERRORSEXITS("DomainFaces_DerivativeVersionNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainFaces_DerivativeVersionNumberGet
  
  !
  !================================================================================================================================
  !

  !>Get the basis for a face in the domain faces identified by its local number
  SUBROUTINE DomainFaces_FaceBasisGet(domainFaces,localFaceNumber,basis,err,error,*)

    !Argument variables
    TYPE(DomainFacesType), POINTER :: domainFaces !<A pointer to the domain faces to get the face basis for
    INTEGER(INTG), INTENT(IN) :: localFaceNumber !<The face local number to get the basis for
    TYPE(BasisType), POINTER, INTENT(OUT) :: basis !<On return, a pointer to the basis for the face. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainFaces_FaceBasisGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainFaces)) CALL FlagError("Domain faces is not associated.",err,error,*999)    
    IF(localFaceNumber<=0.OR.localFaceNumber>domainFaces%numberOfFaces) THEN
      localError="The local face number of "//TRIM(NumberToVString(localFaceNumber,"*",err,error))// &
        & " is invalid. The local number must be >= 1 and <= "// &
        & TRIM(NumberToVString(domainFaces%numberOfFaces,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainFaces%faces)) CALL FlagError("Domain faces faces is not allocated.",err,error,*999)
#endif    
      
    basis=>domainFaces%faces(localFaceNumber)%basis

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(basis)) THEN
      localError="The basis for local face number "//TRIM(NumberToVString(localFaceNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("DomainFaces_FaceBasisGet")
    RETURN
999 NULLIFY(basis)
998 ERRORSEXITS("DomainFaces_FaceBasisGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainFaces_FaceBasisGet

  !  
  !================================================================================================================================
  !

  !>Get the boundary face status for a face in the domain faces identified by its local number
  SUBROUTINE DomainFaces_FaceBoundaryFaceGet(domainFaces,localFaceNumber,boundaryFace,err,error,*)

    !Argument variables
    TYPE(DomainFacesType), POINTER :: domainFaces !<A pointer to the domain faces to get the face boundary status for
    INTEGER(INTG), INTENT(IN) :: localFaceNumber !<The face local number to get the boundary status for
    LOGICAL, INTENT(OUT) :: boundaryFace !<On return, the boundary status for the face.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainFaces_FaceBoundarFaceGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainFaces)) CALL FlagError("Domain faces is not associated.",err,error,*999)    
    IF(localFaceNumber<=0.OR.localFaceNumber>domainFaces%numberOfFaces) THEN
      localError="The local face number of "//TRIM(NumberToVString(localFaceNumber,"*",err,error))// &
        & " is invalid. The local number must be >= 1 and <= "// &
        & TRIM(NumberToVString(domainFaces%numberOfFaces,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainFaces%faces)) CALL FlagError("Domain faces faces is not allocated.",err,error,*999)
#endif    
      
    boundaryFace=domainFaces%faces(localFaceNumber)%boundaryFace

    EXITS("DomainFaces_FaceBoundaryFaceGet")
    RETURN
999 ERRORSEXITS("DomainFaces_FaceBoundaryFaceGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainFaces_FaceBoundaryFaceGet

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
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainFaces_FaceGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(face)) CALL FlagError("Face is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainFaces)) CALL FlagError("Domain faces is not associated.",err,error,*999)
    IF(faceNumber<=0.OR.faceNumber>domainFaces%numberOfFaces) THEN
      localError="The face number of "//TRIM(NumberToVString(faceNumber,"*",err,error))// &
        & " is invalid. The face number must be >= 1 and <= "// &
        & TRIM(NumberToVString(domainFaces%numberOfFaces,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainFaces%faces)) CALL FlagError("Domain faces faces is not associated.",err,error,*999)
#endif    
      
    face=>domainFaces%faces(faceNumber)

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(face)) THEN
      localError="Face number "//TRIM(NumberToVString(faceNumber,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("DomainFaces_FaceGet")
    RETURN
999 NULLIFY(face)
998 ERRORSEXITS("DomainFaces_FaceGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainFaces_FaceGet    
  
  !
  !================================================================================================================================
  !

  !>Gets the node number for a local node index in the domain faces identified by its local number
  SUBROUTINE DomainFaces_FaceNodeNumberGet(domainFaces,localNodeIdx,faceNumber,nodeNumber,err,error,*)

    !Argument variables
    TYPE(DomainFacesType), POINTER :: domainFaces !<A pointer to the domain faces to get the node number for
    INTEGER(INTG), INTENT(IN) :: localNodeIdx !<The local node index to get the face node number.
    INTEGER(INTG), INTENT(IN) :: faceNumber !<The face number to get the node number for.
    INTEGER(INTG), INTENT(OUT) :: nodeNumber !<On return, the node number for the local node index in the specified face.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainFaces_FaceNodeNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainFaces)) CALL FlagError("Domain faces is not associated.",err,error,*999)
    IF(faceNumber<1.OR.faceNumber>domainFaces%numberOfFaces) THEN
      localError="The face number of "//TRIM(NumberToVString(faceNumber,"*",err,error))// &
        & " is invalid. The face number must be >= 1 and <= "// &
        & TRIM(NumberToVString(domainFaces%numberOfFaces,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainFaces%faces)) CALL FlagError("Domain faces faces is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(domainFaces%faces(faceNumber)%nodesInFace)) THEN
      localError="The nodes in face is not allocated for face number "//TRIM(NumberToVString(faceNumber,"*",err,error))// &
        & " in the domain faces."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localNodeIdx<1.OR.localNodeIdx>SIZE(domainFaces%faces(faceNumber)%nodesInFace,1)) THEN
      localError="The specified local node index of "//TRIM(NumberToVString(localNodeIdx,"*",err,error))// &
        & " is invalid for face number "//TRIM(NumberToVString(faceNumber,"*",err,error))// &
        & " in the domain faces. The local node index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(domainFaces%faces(faceNumber)%nodesInFace,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    nodeNumber=domainFaces%faces(faceNumber)%nodesInFace(localNodeIdx)
    
    EXITS("DomainFaces_FaceNodeNumberGet")
    RETURN
999 ERRORSEXITS("DomainFaces_FaceNodeNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainFaces_FaceNodeNumberGet
  
  !
  !================================================================================================================================
  !

  !>Gets the global derivative index for a derivative index of a local node index in the domain lines identified by its local number
  SUBROUTINE DomainLines_DerivativeGlobalIndexGet(domainLines,derivativeIdx,localNodeIdx,lineNumber,globalDerivativeIndex, &
    & err,error,*)

    !Argument variables
    TYPE(DomainLinesType), POINTER :: domainLines !<A pointer to the domain lines to get the global derivative index for
    INTEGER(INTG), INTENT(IN) :: derivativeIdx !<The local derivative index to get the global derivative index for.
    INTEGER(INTG), INTENT(IN) :: localNodeIdx !<The local node index to get the global derivative index.
    INTEGER(INTG), INTENT(IN) :: lineNumber !<The line number to get the node number for.
    INTEGER(INTG), INTENT(OUT) :: globalDerivativeIndex !<On return, the global derivative index for the derivative index of the local node index in the specified line.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainLines_DerivativeGlobalIndexGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainLines)) CALL FlagError("Domain lines is not associated.",err,error,*999)
    IF(lineNumber<1.OR.lineNumber>domainLines%numberOfLines) THEN
      localError="The line number of "//TRIM(NumberToVString(lineNumber,"*",err,error))// &
        & " is invalid. The line number must be >= 1 and <= "// &
        & TRIM(NumberToVString(domainLines%numberOfLines,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainLines%lines)) CALL FlagError("Domain lines lines is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(domainLines%lines(lineNumber)%derivativesInLine)) THEN
      localError="The derivatives in line is not allocated for line number "//TRIM(NumberToVString(lineNumber,"*",err,error))// &
        & " in the domain lines."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(derivativeIdx<1.OR.derivativeIdx>SIZE(domainLines%lines(lineNumber)%derivativesInLine,2)) THEN
      localError="The specified derivative index of "//TRIM(NumberToVString(derivativeIdx,"*",err,error))// &
        & " is invalid for line number "//TRIM(NumberToVString(lineNumber,"*",err,error))// &
        & " in the domain lines. The derivative index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(domainLines%lines(lineNumber)%derivativesInLine,2),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localNodeIdx<1.OR.localNodeIdx>SIZE(domainLines%lines(lineNumber)%derivativesInLine,3)) THEN
      localError="The specified local node index of "//TRIM(NumberToVString(localNodeIdx,"*",err,error))// &
        & " is invalid for line number "//TRIM(NumberToVString(lineNumber,"*",err,error))// &
        & " in the domain lines. The local node index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(domainLines%lines(lineNumber)%derivativesInLine,3),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    globalDerivativeIndex=domainLines%lines(lineNumber)%derivativesInLine(1,derivativeIdx,localNodeIdx)
    
    EXITS("DomainLines_DerivativeGlobalIndexGet")
    RETURN
999 ERRORSEXITS("DomainLines_DerivativeGlobalIndexGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainLines_DerivativeGlobalIndexGet
  
  !
  !================================================================================================================================
  !

  !>Gets the version number for a derivative index of a local node index in the domain lines identified by its local number
  SUBROUTINE DomainLines_DerivativeVersionNumberGet(domainLines,derivativeIdx,localNodeIdx,lineNumber,versionNumber, &
    & err,error,*)

    !Argument variables
    TYPE(DomainLinesType), POINTER :: domainLines !<A pointer to the domain lines to get the version number for
    INTEGER(INTG), INTENT(IN) :: derivativeIdx !<The local derivative index to get the version number for.
    INTEGER(INTG), INTENT(IN) :: localNodeIdx !<The local node index to get the version number for.
    INTEGER(INTG), INTENT(IN) :: lineNumber !<The line number to get the version number for.
    INTEGER(INTG), INTENT(OUT) :: versionNumber !<On return, the version number for the derivative index of the local node index in the specified line.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainLines_DerivativeVersionNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainLines)) CALL FlagError("Domain lines is not associated.",err,error,*999)
    IF(lineNumber<1.OR.lineNumber>domainLines%numberOfLines) THEN
      localError="The line number of "//TRIM(NumberToVString(lineNumber,"*",err,error))// &
        & " is invalid. The line number must be >= 1 and <= "// &
        & TRIM(NumberToVString(domainLines%numberOfLines,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainLines%lines)) CALL FlagError("Domain lines lines is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(domainLines%lines(lineNumber)%derivativesInLine)) THEN
      localError="The derivatives in line is not allocated for line number "//TRIM(NumberToVString(lineNumber,"*",err,error))// &
        & " in the domain lines."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(derivativeIdx<1.OR.derivativeIdx>SIZE(domainLines%lines(lineNumber)%derivativesInLine,2)) THEN
      localError="The specified derivative index of "//TRIM(NumberToVString(derivativeIdx,"*",err,error))// &
        & " is invalid for line number "//TRIM(NumberToVString(lineNumber,"*",err,error))// &
        & " in the domain lines. The derivative index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(domainLines%lines(lineNumber)%derivativesInLine,2),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localNodeIdx<1.OR.localNodeIdx>SIZE(domainLines%lines(lineNumber)%derivativesInLine,3)) THEN
      localError="The specified local node index of "//TRIM(NumberToVString(localNodeIdx,"*",err,error))// &
        & " is invalid for line number "//TRIM(NumberToVString(lineNumber,"*",err,error))// &
        & " in the domain lines. The local node index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(domainLines%lines(lineNumber)%derivativesInLine,3),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    versionNumber=domainLines%lines(lineNumber)%derivativesInLine(2,derivativeIdx,localNodeIdx)
    
    EXITS("DomainLines_DerivativeVersionNumberGet")
    RETURN
999 ERRORSEXITS("DomainLines_DerivativeVersionNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainLines_DerivativeVersionNumberGet
  
  !  
  !================================================================================================================================
  !

  !>Get the basis for a line in the domain lines identified by its local number
  SUBROUTINE DomainLines_LineBasisGet(domainLines,localLineNumber,basis,err,error,*)

    !Argument variables
    TYPE(DomainLinesType), POINTER :: domainLines !<A pointer to the domain lines to get the line basis for
    INTEGER(INTG), INTENT(IN) :: localLineNumber !<The line local number to get the basis for
    TYPE(BasisType), POINTER, INTENT(OUT) :: basis !<On return, a pointer to the basis for the line. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainLiness_LineBasisGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainLines)) CALL FlagError("Domain lines is not associated.",err,error,*999)    
    IF(localLineNumber<=0.OR.localLineNumber>domainLines%numberOfLines) THEN
      localError="The local line number of "//TRIM(NumberToVString(localLineNumber,"*",err,error))// &
        & " is invalid. The local number must be >= 1 and <= "// &
        & TRIM(NumberToVString(domainLines%numberOfLines,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainLines%lines)) CALL FlagError("Domain lines lines is not allocated.",err,error,*999)
#endif    
      
    basis=>domainLines%lines(localLineNumber)%basis

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(basis)) THEN
      localError="The basis for local line number "//TRIM(NumberToVString(localLineNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("DomainLines_LineBasisGet")
    RETURN
999 NULLIFY(basis)
998 ERRORSEXITS("DomainLines_LineBasisGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainLines_LineBasisGet

  !  
  !================================================================================================================================
  !

  !>Get the boundary line status for a line in the domain lines identified by its local number
  SUBROUTINE DomainLines_LineBoundaryLineGet(domainLines,localLineNumber,boundaryLine,err,error,*)

    !Argument variables
    TYPE(DomainLinesType), POINTER :: domainLines !<A pointer to the domain lines to get the line boundary status for
    INTEGER(INTG), INTENT(IN) :: localLineNumber !<The line local number to get the boundary status for
    LOGICAL, INTENT(OUT) :: boundaryLine !<On return, the boundary status for the line.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainLiness_LineBoundaryLineGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainLines)) CALL FlagError("Domain lines is not associated.",err,error,*999)    
    IF(localLineNumber<=0.OR.localLineNumber>domainLines%numberOfLines) THEN
      localError="The local line number of "//TRIM(NumberToVString(localLineNumber,"*",err,error))// &
        & " is invalid. The local number must be >= 1 and <= "// &
        & TRIM(NumberToVString(domainLines%numberOfLines,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainLines%lines)) CALL FlagError("Domain lines lines is not allocated.",err,error,*999)
#endif    
      
    boundaryLine=domainLines%lines(localLineNumber)%boundaryLine

    EXITS("DomainLines_LineBoundaryLineGet")
    RETURN
999 ERRORSEXITS("DomainLines_LineBoundaryLineGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainLines_LineBoundaryLineGet

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
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainLines_LineGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(line)) CALL FlagError("Line is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainLines)) CALL FlagError("Domain lines is not associated.",err,error,*999)
    IF(lineNumber<=0.OR.lineNumber>domainLines%numberOfLines) THEN
      localError="The line number of "//TRIM(NumberToVString(lineNumber,"*",err,error))// &
        & " is invalid. The line number must be >= 1 and <= "// &
        & TRIM(NumberToVString(domainLines%numberOfLines,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainLines%lines)) CALL FlagError("Domain lines lines is not associated.",err,error,*999)
#endif    
      
    line=>domainLines%lines(lineNumber)

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(line)) THEN
      localError="Line number "//TRIM(NumberToVString(lineNumber,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("DomainLines_LineGet")
    RETURN
999 NULLIFY(line)
998 ERRORSEXITS("DomainLines_LineGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainLines_LineGet
  
  !
  !================================================================================================================================
  !

  !>Gets the node number for a local node index in the domain lines identified by its local number
  SUBROUTINE DomainLines_LineNodeNumberGet(domainLines,localNodeIdx,lineNumber,nodeNumber,err,error,*)

    !Argument variables
    TYPE(DomainLinesType), POINTER :: domainLines !<A pointer to the domain lines to get the node number for
    INTEGER(INTG), INTENT(IN) :: localNodeIdx !<The local node index to get the line node number.
    INTEGER(INTG), INTENT(IN) :: lineNumber !<The line number to get the node number for.
    INTEGER(INTG), INTENT(OUT) :: nodeNumber !<On return, the node number for the local node index in the specified line.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainLines_LineNodeNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainLines)) CALL FlagError("Domain lines is not associated.",err,error,*999)
    IF(lineNumber<1.OR.lineNumber>domainLines%numberOfLines) THEN
      localError="The line number of "//TRIM(NumberToVString(lineNumber,"*",err,error))// &
        & " is invalid. The line number must be >= 1 and <= "// &
        & TRIM(NumberToVString(domainLines%numberOfLines,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainLines%lines)) CALL FlagError("Domain lines lines is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(domainLines%lines(lineNumber)%nodesInLine)) THEN
      localError="The nodes in line is not allocated for line number "//TRIM(NumberToVString(lineNumber,"*",err,error))// &
        & " in the domain lines."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localNodeIdx<1.OR.localNodeIdx>SIZE(domainLines%lines(lineNumber)%nodesInLine,1)) THEN
      localError="The specified local node index of "//TRIM(NumberToVString(localNodeIdx,"*",err,error))// &
        & " is invalid for line number "//TRIM(NumberToVString(lineNumber,"*",err,error))// &
        & " in the domain lines. The local node index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(domainLines%lines(lineNumber)%nodesInLine,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    nodeNumber=domainLines%lines(lineNumber)%nodesInLine(localNodeIdx)
    
    EXITS("DomainLines_LineNodeNumberGet")
    RETURN
999 ERRORSEXITS("DomainLines_LineNodeNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainLines_LineNodeNumberGet
  
  !
  !================================================================================================================================
  !

  !>Gets a domain decomposition for a domain.
  SUBROUTINE Domain_DecompositionGet(domain,decomposition,err,error,*)

    !Argument variables
    TYPE(DomainType), POINTER :: domain !<The domain to get the decomposition for.
    TYPE(DecompositionType), POINTER :: decomposition !<On exit, a pointer to the domain decomposition. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Domain_DecompositionGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(decomposition)) CALL FlagError("Decomposition is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domain)) CALL FlagError("Domain is not associated.",err,error,*999)
#endif    

    !Get the domain decomposition
    decomposition=>domain%decomposition

#ifdef WITH_POSTCHECKS    
    !Check decomposition is associated.
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
#endif    
    
    EXITS("Domain_DecompositionGet")
    RETURN
999 NULLIFY(decomposition)
998 ERRORSEXITS("Domain_DecompositionGet",err,error)
    RETURN 1
    
  END SUBROUTINE Domain_DecompositionGet

  !
  !================================================================================================================================
  !

  !>Gets a domain mappings from a domain.
  SUBROUTINE Domain_DomainMappingsGet(domain,domainMappings,err,error,*)

    !Argument variables
    TYPE(DomainType), POINTER :: domain !<The domain to get the mappings for.
    TYPE(DomainMappingsType), POINTER :: domainMappings !<On exit, a pointer to the domain mappings. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Domain_DomainMappingsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(.NOT.ASSOCIATED(domain)) CALL FlagError("Domain is not associated.",err,error,*999)
    IF(ASSOCIATED(domainMappings)) CALL FlagError("Domain mappings is already associated.",err,error,*999)
#endif    

    !Get the domain mappings
    domainMappings=>domain%mappings

#ifdef WITH_POSTCHECKS    
    !Check domain mappings is associated.
    IF(.NOT.ASSOCIATED(domainMappings)) CALL FlagError("Domain mappings is not associated.",err,error,*999)
#endif    
    
    EXITS("Domain_DomainMappingsGet")
    RETURN
999 ERRORSEXITS("Domain_DomainMappingsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Domain_DomainMappingsGet

  !
  !================================================================================================================================
  !

  !>Gets a domain topology from a domain.
  SUBROUTINE Domain_DomainTopologyGet(domain,domainTopology,err,error,*)

    !Argument variables
    TYPE(DomainType), POINTER :: domain !<The domain to get the topologya for.
    TYPE(DomainTopologyType), POINTER :: domainTopology !<On exit, a pointer to the domain topology. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Domain_DomainTopologyGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(.NOT.ASSOCIATED(domain)) CALL FlagError("Domain is not associated.",err,error,*999)
    IF(ASSOCIATED(domainTopology)) CALL FlagError("Domain topology is already associated.",err,error,*999)
#endif    

    !Get the domain topology
    domainTopology=>domain%topology

#ifdef WITH_POSTCHECKS    
    !Check domain topology is associated.
    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain topology is not associated.",err,error,*999)
#endif    
    
    EXITS("Domain_DomainTopologyGet")
    RETURN
999 ERRORSEXITS("Domain_DomainTopologyGet",err,error)
    RETURN 1
    
  END SUBROUTINE Domain_DomainTopologyGet

  !
  !================================================================================================================================
  !

  !>Gets a mesh component number from a domain.
  SUBROUTINE Domain_MeshComponentNumberGet(domain,meshComponentNumber,err,error,*)

    !Argument variables
    TYPE(DomainType), POINTER :: domain !<The domain to get the mesh component number for.
    INTEGER(INTG), INTENT(OUT) :: meshComponentNumber !<On exit, the mesh component number for the domain.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Domain_MeshComponentNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(.NOT.ASSOCIATED(domain)) CALL FlagError("Domain is not associated.",err,error,*999)
#endif    

    !Get the domain mesh component number
    meshComponentNumber=domain%meshComponentNumber

    EXITS("Domain_MeshComponentNumberGet")
    RETURN
999 ERRORSEXITS("Domain_MeshComponentNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE Domain_MeshComponentNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the number of dimensions from a domain.
  SUBROUTINE Domain_NumberOfDimensionsGet(domain,numberOfDimensions,err,error,*)

    !Argument variables
    TYPE(DomainType), POINTER :: domain !<The domain to get the number of dimensions for.
    INTEGER(INTG), INTENT(OUT) :: numberOfDimensions !<On exit, the number of dimensions for the domain.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Domain_NumberOfDimensionsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(.NOT.ASSOCIATED(domain)) CALL FlagError("Domain is not associated.",err,error,*999)
#endif    

    !Get the domain number of dimensions
    numberOfDimensions=domain%numberOfDimensions

    EXITS("Domain_NumberOfDimensionsGet")
    RETURN
999 ERRORSEXITS("Domain_NumberOfDimensionsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Domain_NumberOfDimensionsGet

  !
  !================================================================================================================================
  !

  !>Gets the region from a domain.
  SUBROUTINE Domain_RegionGet(domain,region,err,error,*)

    !Argument variables
    TYPE(DomainType), POINTER :: domain !<The domain to get the region for.
    TYPE(RegionType), POINTER :: region !<On exit, a pointer to the domain region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Domain_RegionGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(region)) CALL FlagError("Region is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domain)) CALL FlagError("Domain is not associated.",err,error,*999)
#endif    

    !Get the domain region
    region=>domain%region

#ifdef WITH_POSTCHECKS    
    !Check domain region is associated.
    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Domain region is not associated.",err,error,*999)
#endif    
    
    EXITS("Domain_RegionGet")
    RETURN
999 NULLIFY(region)
998 ERRORSEXITS("Domain_RegionGet",err,error)
    RETURN 1
    
  END SUBROUTINE Domain_RegionGet

  !
  !================================================================================================================================
  !

  !>Gets dofs from a domain mappings.
  SUBROUTINE DomainMappings_DOFsMappingGet(domainMappings,domainDofs,err,error,*)

    !Argument variables
    TYPE(DomainMappingsType), POINTER :: domainMappings !<A pointer to the domain mappings to get the dofs for.
    TYPE(DomainMappingType), POINTER :: domainDofs !<On exit, a pointer to the domain mapping dofs. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DomainMappings_DOFsMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(domainDofs)) CALL FlagError("Domain mapping dofs is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainMappings)) CALL FlagError("Domain mappings is not associated.",err,error,*999)
#endif    

    !Get the domain dofs
    domainDofs=>domainMappings%dofs

#ifdef WITH_POSTCHECKS    
    !Check domain dofs is associated.
    IF(.NOT.ASSOCIATED(domainDofs)) CALL FlagError("Domain mappings dofs is not associated.",err,error,*999)
#endif    
    
    EXITS("DomainMappings_DOFsMappingGet")
    RETURN
999 NULLIFY(domainDofs)
998 ERRORSEXITS("DomainMappings_DOFsMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMappings_DOFsMappingGet

  !
  !================================================================================================================================
  !

  !>Gets the domain from a domain mappings.
  SUBROUTINE DomainMappings_DomainGet(domainMappings,domain,err,error,*)

    !Argument variables
    TYPE(DomainMappingsType), POINTER :: domainMappings !<A pointer to the domain mappings to get the domain for.
    TYPE(DomainType), POINTER :: domain !<On exit, a pointer to the domain mapping domain. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DomainMappings_DomainGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(domain)) CALL FlagError("Domain is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainMappings)) CALL FlagError("Domain mappings is not associated.",err,error,*999)
#endif    

    !Get the domain
    domain=>domainMappings%domain

#ifdef WITH_POSTCHECKS    
    !Check domain is associated.
    IF(.NOT.ASSOCIATED(domain)) CALL FlagError("Domain mappings domain is not associated.",err,error,*999)
#endif    
    
    EXITS("DomainMappings_DomainGet")
    RETURN
999 NULLIFY(domain)
998 ERRORSEXITS("DomainMappings_DomainGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMappings_DomainGet

  !
  !================================================================================================================================
  !

  !>Gets elements from a domain mappings.
  SUBROUTINE DomainMappings_ElementsMappingGet(domainMappings,domainElements,err,error,*)

    !Argument variables
    TYPE(DomainMappingsType), POINTER :: domainMappings !<A pointer to the domain mappings to get the elements for.
    TYPE(DomainMappingType), POINTER :: domainElements !<On exit, a pointer to the domain mapping elements. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DomainMappings_ElementsMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(domainElements)) CALL FlagError("Domain mapping elements is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainMappings)) CALL FlagError("Domain mappings is not associated.",err,error,*999)
#endif    

    !Get the domain elements
    domainElements=>domainMappings%elements

#ifdef WITH_POSTCHECKS    
    !Check domain elements is associated.
    IF(.NOT.ASSOCIATED(domainElements)) CALL FlagError("Domain mappings elements is not associated.",err,error,*999)
#endif    
    
    EXITS("DomainMappings_ElementsMappingGet")
    RETURN
999 NULLIFY(domainElements)
998 ERRORSEXITS("DomainMappings_ElementsMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMappings_ElementsMappingGet

  !
  !================================================================================================================================
  !

  !>Gets nodes from a domain mappings.
  SUBROUTINE DomainMappings_NodesMappingGet(domainMappings,domainNodes,err,error,*)

    !Argument variables
    TYPE(DomainMappingsType), POINTER :: domainMappings !<A pointer to the domain mappings to get the nodes for.
    TYPE(DomainMappingType), POINTER :: domainNodes !<On exit, a pointer to the domain mapping nodes. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DomainMappings_NodesMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(domainNodes)) CALL FlagError("Domain mapping nodes is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainMappings)) CALL FlagError("Domain mappings is not associated.",err,error,*999)
#endif    

    !Get the domain nodes
    domainNodes=>domainMappings%nodes

#ifdef WITH_POSTCHECKS    
    !Check domain nodes is associated.
    IF(.NOT.ASSOCIATED(domainNodes)) CALL FlagError("Domain mappings nodes is not associated.",err,error,*999)
#endif    
    
    EXITS("DomainMappings_NodesMappingGet")
    RETURN
999 NULLIFY(domainNodes)
998 ERRORSEXITS("DomainMappings_NodesMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMappings_NodesMappingGet

  !
  !================================================================================================================================
  !

  !>Gets the boundary node flag of a domain node. 
  SUBROUTINE DomainNode_BoundaryNodeGet(domainNode,boundaryNode,err,error,*)

    !Argument variables
    TYPE(DomainNodeType), POINTER :: domainNode !<A pointer to the domain node to get the boundary node flag for
    LOGICAL, INTENT(OUT) :: boundaryNode !<On exit, the boundary node flag of the domain node.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainNode_BoundaryNodeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainNode)) CALL FlagError("Domain node is not associated.",err,error,*999)
#endif

    boundaryNode=domainNode%boundaryNode

    EXITS("DomainNode_BoundaryNodeGet")
    RETURN
999 ERRORSEXITS("DomainNode_BoundaryNodeGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNode_BoundaryNodeGet

  !
  !================================================================================================================================
  !

  !>Gets a node derivative for a domain node. 
  SUBROUTINE DomainNode_NodeDerivativeGet(domainNode,derivativeIdx,domainNodeDerivative,err,error,*)

    !Argument variables
    TYPE(DomainNodeType), POINTER :: domainNode !<A pointer to the domain node to get the node derivative for
    INTEGER(INTG), INTENT(IN) :: derivativeIdx !<The derivative index to get the node deriviative for
    TYPE(DomainNodeDerivativeType), POINTER :: domainNodeDerivative !<On exit, a pointer the specified node derivative for the domain node. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainNode_NodeDerivativeGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(domainNodeDerivative)) CALL FlagError("Domain node derivative is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainNode)) CALL FlagError("Domain node is not associated.",err,error,*999)
    IF(derivativeIdx<1.OR.derivativeIdx>domainNode%numberOfDerivatives) THEN
      localError="The specified derivative index of "//TRIM(NumberToVString(derivativeIdx,"*",err,error))// &
        & " is invalid. The derivative index should be >= 1 and <= "// &
        & TRIM(NumberToVString(domainNode%numberOfDerivatives,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    IF(.NOT.ALLOCATED(domainNode%derivatives)) &
      & CALL FlagError("Derivatives is not allocated for the domain node.",err,error,*999)
#endif

    domainNodeDerivative=>domainNode%derivatives(derivativeIdx)

#ifdef WITH_POSTCHECKS
    IF(.NOT.ASSOCIATED(domainNodeDerivative)) THEN
      localError="The domain node derivative is not associated for derivative index "// &
        & TRIM(NumberToVString(derivativeIdx,"*",err,error))//" of the domain node."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("DomainNode_NodeDerivativeGet")
    RETURN
999 NULLIFY(domainNodeDerivative)
998 ERRORSEXITS("DomainNode_NodeDerivativeGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNode_NodeDerivativeGet

  !
  !================================================================================================================================
  !

  !>Gets a global node number of a domain node. 
  SUBROUTINE DomainNode_GlobalNodeNumberGet(domainNode,globalNodeNumber,err,error,*)

    !Argument variables
    TYPE(DomainNodeType), POINTER :: domainNode !<A pointer to the domain node to get the global node number for
    INTEGER(INTG), INTENT(OUT) :: globalNodeNumber !<On exit, the global node number of the domain node.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainNode_GlobalNodeNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainNode)) CALL FlagError("Domain node is not associated.",err,error,*999)
#endif

    globalNodeNumber=domainNode%globalNumber

    EXITS("DomainNode_GlobalNodeNumberGet")
    RETURN
999 ERRORSEXITS("DomainNode_GlobalNodeNumberGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNode_GlobalNodeNumberGet

  !
  !================================================================================================================================
  !

  !>Gets a local node number of a domain node. 
  SUBROUTINE DomainNode_LocalNodeNumberGet(domainNode,localNodeNumber,err,error,*)

    !Argument variables
    TYPE(DomainNodeType), POINTER :: domainNode !<A pointer to the domain node to get the local node number for
    INTEGER(INTG), INTENT(OUT) :: localNodeNumber !<On exit, the local node number of the domain node.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainNode_LocalNodeNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainNode)) CALL FlagError("Domain node is not associated.",err,error,*999)
#endif

    localNodeNumber=domainNode%localNumber

    EXITS("DomainNode_LocalNodeNumberGet")
    RETURN
999 ERRORSEXITS("DomainNode_LocalNodeNumberGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNode_LocalNodeNumberGet

  !
  !================================================================================================================================
  !

  !>Gets a mesh node number of a domain node. 
  SUBROUTINE DomainNode_MeshNodeNumberGet(domainNode,meshNodeNumber,err,error,*)

    !Argument variables
    TYPE(DomainNodeType), POINTER :: domainNode !<A pointer to the domain node to get the mesh node number for
    INTEGER(INTG), INTENT(OUT) :: meshNodeNumber !<On exit, the mesh node number of the domain node.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Mesh Variables

    ENTERS("DomainNode_MeshNodeNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainNode)) CALL FlagError("Domain node is not associated.",err,error,*999)
#endif

    meshNodeNumber=domainNode%meshNumber

    EXITS("DomainNode_MeshNodeNumberGet")
    RETURN
999 ERRORSEXITS("DomainNode_MeshNodeNumberGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNode_MeshNodeNumberGet

  !
  !================================================================================================================================
  !

  !>Gets a node face number for a domain node. 
  SUBROUTINE DomainNode_NodeFaceGet(domainNode,nodeFaceIdx,nodeFaceNumber,err,error,*)

    !Argument variables
    TYPE(DomainNodeType), POINTER :: domainNode !<A pointer to the domain node to get the node face number for
    INTEGER(INTG), INTENT(IN) :: nodeFaceIdx !<The node face index to get the node face number for
    INTEGER(INTG), INTENT(OUT) :: nodeFaceNumber !<On exit, the specified node face number for the domain node.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainNode_NodeFaceGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainNode)) CALL FlagError("Domain node is not associated.",err,error,*999)
    IF(nodeFaceIdx<1.OR.nodeFaceIdx>domainNode%numberOfNodeFaces) THEN
      localError="The specified node face index of "//TRIM(NumberToVString(nodeFaceIdx,"*",err,error))// &
        & " is invalid. The node face index should be >= 1 and <= "// &
        & TRIM(NumberToVString(domainNode%numberOfNodeFaces,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    IF(.NOT.ALLOCATED(domainNode%nodeFaces)) &
      & CALL FlagError("Node faces is not allocated for the domain node.",err,error,*999)
#endif

    nodeFaceNumber=domainNode%nodeFaces(nodeFaceIdx)

    EXITS("DomainNode_NodeFaceGet")
    RETURN
999 ERRORSEXITS("DomainNode_NodeFaceGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNode_NodeFaceGet

  !
  !================================================================================================================================
  !

  !>Gets a node line number for a domain node. 
  SUBROUTINE DomainNode_NodeLineGet(domainNode,nodeLineIdx,nodeLineNumber,err,error,*)

    !Argument variables
    TYPE(DomainNodeType), POINTER :: domainNode !<A pointer to the domain node to get the node line number for
    INTEGER(INTG), INTENT(IN) :: nodeLineIdx !<The node line index to get the node line number for
    INTEGER(INTG), INTENT(OUT) :: nodeLineNumber !<On exit, the specified node line number for the domain node.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainNode_NodeLineGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainNode)) CALL FlagError("Domain node is not associated.",err,error,*999)
    IF(nodeLineIdx<1.OR.nodeLineIdx>domainNode%numberOfNodeLines) THEN
      localError="The specified node line index of "//TRIM(NumberToVString(nodeLineIdx,"*",err,error))// &
        & " is invalid. The node line index should be >= 1 and <= "// &
        & TRIM(NumberToVString(domainNode%numberOfNodeLines,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    IF(.NOT.ALLOCATED(domainNode%nodeLines)) &
      & CALL FlagError("Node lines is not allocated for the domain node.",err,error,*999)
#endif

    nodeLineNumber=domainNode%nodeLines(nodeLineIdx)

    EXITS("DomainNode_NodeLineGet")
    RETURN
999 ERRORSEXITS("DomainNode_NodeLineGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNode_NodeLineGet

  !
  !================================================================================================================================
  !

  !>Gets the number of derivatives for a domain node. 
  SUBROUTINE DomainNode_NumberOfDerivativesGet(domainNode,numberOfDerivatives,err,error,*)

    !Argument variables
    TYPE(DomainNodeType), POINTER :: domainNode !<A pointer to the domain node to get the number of derivatives for
    INTEGER(INTG), INTENT(OUT) :: numberOfDerivatives !<On exit, the number of derivatives for the domain node.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Mesh Variables

    ENTERS("DomainNode_NumberOfDerivativesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainNode)) CALL FlagError("Domain node is not associated.",err,error,*999)
#endif

    numberOfDerivatives=domainNode%numberOfDerivatives

    EXITS("DomainNode_NumberOfDerivativesGet")
    RETURN
999 ERRORSEXITS("DomainNode_NumberOfDerivativesGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNode_NumberOfDerivativesGet

  !
  !================================================================================================================================
  !

  !>Gets the number of node faces for a domain node. 
  SUBROUTINE DomainNode_NumberOfNodeFacesGet(domainNode,numberOfNodeFaces,err,error,*)

    !Argument variables
    TYPE(DomainNodeType), POINTER :: domainNode !<A pointer to the domain node to get the number of node faces for
    INTEGER(INTG), INTENT(OUT) :: numberOfNodeFaces !<On exit, the number of node faces for the domain node.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Mesh Variables

    ENTERS("DomainNode_NumberOfNodeFacesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainNode)) CALL FlagError("Domain node is not associated.",err,error,*999)
#endif

    numberOfNodeFaces=domainNode%numberOfNodeFaces

    EXITS("DomainNode_NumberOfNodeFacesGet")
    RETURN
999 ERRORSEXITS("DomainNode_NumberOfNodeFacesGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNode_NumberOfNodeFacesGet

  !
  !================================================================================================================================
  !

  !>Gets the number of node lines for a domain node. 
  SUBROUTINE DomainNode_NumberOfNodeLinesGet(domainNode,numberOfNodeLines,err,error,*)

    !Argument variables
    TYPE(DomainNodeType), POINTER :: domainNode !<A pointer to the domain node to get the number of node lines for
    INTEGER(INTG), INTENT(OUT) :: numberOfNodeLines !<On exit, the number of node lines for the domain node.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Mesh Variables

    ENTERS("DomainNode_NumberOfNodeLinesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainNode)) CALL FlagError("Domain node is not associated.",err,error,*999)
#endif

    numberOfNodeLines=domainNode%numberOfNodeLines

    EXITS("DomainNode_NumberOfNodeLinesGet")
    RETURN
999 ERRORSEXITS("DomainNode_NumberOfNodeLinesGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNode_NumberOfNodeLinesGet

  !
  !================================================================================================================================
  !

  !>Gets the number of surrounding elements for a domain node. 
  SUBROUTINE DomainNode_NumberOfSurroundingElementsGet(domainNode,numberOfSurroundingElements,err,error,*)

    !Argument variables
    TYPE(DomainNodeType), POINTER :: domainNode !<A pointer to the domain node to get the number of surrounding elements for
    INTEGER(INTG), INTENT(OUT) :: numberOfSurroundingElements !<On exit, the number of surrounding for the domain node.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Mesh Variables

    ENTERS("DomainNode_NumberOfSurroundingElementsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainNode)) CALL FlagError("Domain node is not associated.",err,error,*999)
#endif

    numberOfSurroundingElements=domainNode%numberOfSurroundingElements

    EXITS("DomainNode_NumberOfSurroundingElementsGet")
    RETURN
999 ERRORSEXITS("DomainNode_NumberOfSurroundingElementsGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNode_NumberOfSurroundingElementsGet

  !
  !================================================================================================================================
  !

  !>Gets a surrounding element number for a domain node. 
  SUBROUTINE DomainNode_SurroundingElementGet(domainNode,surroundingElementIdx,surroundingElementNumber,err,error,*)

    !Argument variables
    TYPE(DomainNodeType), POINTER :: domainNode !<A pointer to the domain node to get the surrounding element number for
    INTEGER(INTG), INTENT(IN) :: surroundingElementIdx !<The surrounding element index to get the surrounding element number for
    INTEGER(INTG), INTENT(OUT) :: surroundingElementNumber !<On exit, the specified surrounding element number for the domain node.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainNode_SurroundingElementGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainNode)) CALL FlagError("Domain node is not associated.",err,error,*999)
    IF(surroundingElementIdx<1.OR.surroundingElementIdx>domainNode%numberOfSurroundingElements) THEN
      localError="The specified surrounding element index of "//TRIM(NumberToVString(surroundingElementIdx,"*",err,error))// &
        & " is invalid. The surrounding element index should be >= 1 and <= "// &
        & TRIM(NumberToVString(domainNode%numberOfSurroundingElements,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    IF(.NOT.ALLOCATED(domainNode%surroundingElements)) &
      & CALL FlagError("Surrounding elements is not allocated for the domain node.",err,error,*999)
#endif

    surroundingElementNumber=domainNode%surroundingElements(surroundingElementIdx)

    EXITS("DomainNode_SurroundingElementGet")
    RETURN
999 ERRORSEXITS("DomainNode_SurroundingElementGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNode_SurroundingElementGet

  !
  !================================================================================================================================
  !

  !>Gets a user node number of a domain node. 
  SUBROUTINE DomainNode_UserNodeNumberGet(domainNode,userNodeNumber,err,error,*)

    !Argument variables
    TYPE(DomainNodeType), POINTER :: domainNode !<A pointer to the domain node to get the user node number for
    INTEGER(INTG), INTENT(OUT) :: userNodeNumber !<On exit, the user node number of the domain node.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainNode_UserNodeNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainNode)) CALL FlagError("Domain node is not associated.",err,error,*999)
#endif

    userNodeNumber=domainNode%userNumber

    EXITS("DomainNode_UserNodeNumberGet")
    RETURN
999 ERRORSEXITS("DomainNode_UserNodeNumberGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNode_UserNodeNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the derivative number for a domain node derivative. 
  SUBROUTINE DomainNodeDerivative_DerivativeGlobalIndexGet(domainNodeDerivative,derivativeNumber,err,error,*)

    !Argument variables
    TYPE(DomainNodeDerivativeType), POINTER :: domainNodeDerivative !<A pointer to the domain node derivative to get the derivative number for
    INTEGER(INTG), INTENT(OUT) :: derivativeNumber !<On exit, the derivative number for the domain node derivative.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Mesh Variables

    ENTERS("DomainNodeDerivative_DerivativeGlobalIndexGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainNodeDerivative)) CALL FlagError("Domain node derivative is not associated.",err,error,*999)
#endif

    derivativeNumber=domainNodeDerivative%globalDerivativeIndex

    EXITS("DomainNodeDerivatives_DerivativeNumberGet")
    RETURN
999 ERRORSEXITS("DomainNodeDerivative_DerivativeGlobalIndexGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNodeDerivative_DerivativeGlobalIndexGet

  !
  !================================================================================================================================
  !

  !>Gets the number of versions for a domain node derivative. 
  SUBROUTINE DomainNodeDerivative_NumberOfVersionsGet(domainNodeDerivative,numberOfVersions,err,error,*)

    !Argument variables
    TYPE(DomainNodeDerivativeType), POINTER :: domainNodeDerivative !<A pointer to the domain node derivative to get the number of versions for
    INTEGER(INTG), INTENT(OUT) :: numberOfVersions !<On exit, the number of versions for the domain node derivative.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Mesh Variables

    ENTERS("DomainNodeDerivative_NumberOfVersionsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainNodeDerivative)) CALL FlagError("Domain node derivative is not associated.",err,error,*999)
#endif

    numberOfVersions=domainNodeDerivative%numberOfVersions

    EXITS("DomainNodeDerivatives_NumberOfVersionsGet")
    RETURN
999 ERRORSEXITS("DomainNodeDerivative_NumberOfVersionsGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNodeDerivative_NumberOfVersionsGet

  !
  !================================================================================================================================
  !

  !>Gets the partial derivative index for a domain node derivative. 
  SUBROUTINE DomainNodeDerivative_PartialDerivativeIndexGet(domainNodeDerivative,partialDerivativeIndex,err,error,*)

    !Argument variables
    TYPE(DomainNodeDerivativeType), POINTER :: domainNodeDerivative !<A pointer to the domain node derivative to get the partial derivative index for
    INTEGER(INTG), INTENT(OUT) :: partialDerivativeIndex !<On exit, the partial derivative index for the domain node derivative.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Mesh Variables

    ENTERS("DomainNodeDerivative_ParitalDerivativeIndexGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainNodeDerivative)) CALL FlagError("Domain node derivative is not associated.",err,error,*999)
#endif

    partialDerivativeIndex=domainNodeDerivative%partialDerivativeIndex

    EXITS("DomainNodeDerivatives_PartialDerivativeIndexGet")
    RETURN
999 ERRORSEXITS("DomainNodeDerivative_PartialDerivativeIndexGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNodeDerivative_PartialDerivativeIndexGet

  !
  !================================================================================================================================
  !

  !>Gets a version number for a domain node derivative. 
  SUBROUTINE DomainNodeDerivative_VersionNumberGet(domainNodeDerivative,versionIdx,versionNumber,err,error,*)

    !Argument variables
    TYPE(DomainNodeDerivativeType), POINTER :: domainNodeDerivative !<A pointer to the domain node derivative to get the version number for
    INTEGER(INTG), INTENT(IN) :: versionIdx !<The surrounding element index to get the surrounding element number for
    INTEGER(INTG), INTENT(OUT) :: versionNumber !<On exit, the specified version number for the domain node derivative.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainNodeDerivative_VersionNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainNodeDerivative)) CALL FlagError("Domain node derivative is not associated.",err,error,*999)
    IF(versionIdx<1.OR.versionIdx>domainNodeDerivative%numberOfversions) THEN
      localError="The specified version index of "//TRIM(NumberToVString(versionIdx,"*",err,error))// &
        & " is invalid. The version index should be >= 1 and <= "// &
        & TRIM(NumberToVString(domainNodeDerivative%numberOfVersions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    IF(.NOT.ALLOCATED(domainNodeDerivative%userVersionNumbers)) &
      & CALL FlagError("User version numbers is not allocated for the domain node derivative.",err,error,*999)
#endif

    versionNumber=domainNodeDerivative%userVersionNumbers(versionIdx)

    EXITS("DomainNodeDerivative_VersionNumberGet")
    RETURN
999 ERRORSEXITS("DomainNodeDerivative_VersionNumberGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNodeDerivative_VersionNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the derivative number for a derivative index of local node number from a domain. 
  SUBROUTINE DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,localNodeNumber,derivativeNumber,err,error,*)

    !Argument variables
    TYPE(DomainNodesType), POINTER :: domainNodes !<A pointer to the domain nodes to get the derivative number for
    INTEGER(INTG), INTENT(IN) :: derivativeIdx !<The derivative index to get the derivative number for
    INTEGER(INTG), INTENT(IN) :: localNodeNumber !<The local node number to get the derivative number for
    INTEGER(INTG), INTENT(OUT) :: derivativeNumber !<On exit, the derivative number for the derivative index of the local node number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DomainNodeType), POINTER :: domainNode
    TYPE(DomainNodeDerivativeType), POINTER :: domainNodeDerivative
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainNodes_DerivativeGlobalIndexGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(domainNode)
    CALL DomainNodes_NodeGet(domainNodes,localNodeNumber,domainNode,err,error,*999)
    NULLIFY(domainNodeDerivative)
    CALL DomainNode_NodeDerivativeGet(domainNode,derivativeIdx,domainNodeDerivative,err,error,*999)
#endif

    derivativeNumber=domainNodes%nodes(localNodeNumber)%derivatives(derivativeIdx)%globalDerivativeIndex

    EXITS("DomainNodes_DerivativeGlobalIndexGet")
    RETURN
999 ERRORSEXITS("DomainNodes_DerivativeGlobalIndexGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNodes_DerivativeGlobalIndexGet

  !
  !================================================================================================================================
  !

  !>Gets the number of versions for a derivative of local node number from a domain. 
  SUBROUTINE DomainNodes_DerivativeNumberOfVersionsGet(domainNodes,derivativeIdx,localNodeNumber,numberOfVersions,err,error,*)

    !Argument variables
    TYPE(DomainNodesType), POINTER :: domainNodes !<A pointer to the domain nodes to get the number of versions for
    INTEGER(INTG), INTENT(IN) :: derivativeIdx !<The derivative index to get the number of versions for
    INTEGER(INTG), INTENT(IN) :: localNodeNumber !<The local node number to get the derivative number for
    INTEGER(INTG), INTENT(OUT) :: numberOfVersions !<On exit, the number of versions for the derivative index of the local node number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DomainNodeType), POINTER :: domainNode
    TYPE(DomainNodeDerivativeType), POINTER :: domainNodeDerivative
#endif    

    ENTERS("DomainNodes_DerivativeNumberOfVersionsGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(domainNode)
    CALL DomainNodes_NodeGet(domainNodes,localNodeNumber,domainNode,err,error,*999)
    NULLIFY(domainNodeDerivative)
    CALL DomainNode_NodeDerivativeGet(domainNode,derivativeIdx,domainNodeDerivative,err,error,*999)
#endif

    numberOfVersions=domainNodes%nodes(localNodeNumber)%derivatives(derivativeIdx)%numberOfVersions

    EXITS("DomainNodes_DerivativeNumberOfVersionsGet")
    RETURN
999 ERRORSEXITS("DomainNodes_DerivativeNumberOfVersionsGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNodes_DerivativeNumberOfVersionsGet

  !
  !================================================================================================================================
  !

  !>Gets the partial derivative index for a derivative index of local node number from a domain. 
  SUBROUTINE DomainNodes_DerivativePartialIndexGet(domainNodes,derivativeIdx,localNodeNumber,partialDerivativeIndex,err,error,*)

    !Argument variables
    TYPE(DomainNodesType), POINTER :: domainNodes !<A pointer to the domain nodes to get the partial derivative index for
    INTEGER(INTG), INTENT(IN) :: derivativeIdx !<The derivative index to get the partial derivative index for
    INTEGER(INTG), INTENT(IN) :: localNodeNumber !<The local node number to get the partial derivative index for
    INTEGER(INTG), INTENT(OUT) :: partialDerivativeIndex !<On exit, the partial derivative index for the derivative index of the local node number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DomainNodeType), POINTER :: domainNode
    TYPE(DomainNodeDerivativeType), POINTER :: domainNodeDerivative
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainNodes_DerivativePartialIndexGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(domainNode)
    CALL DomainNodes_NodeGet(domainNodes,localNodeNumber,domainNode,err,error,*999)
    NULLIFY(domainNodeDerivative)
    CALL DomainNode_NodeDerivativeGet(domainNode,derivativeIdx,domainNodeDerivative,err,error,*999)
#endif

    partialDerivativeIndex=domainNodes%nodes(localNodeNumber)%derivatives(derivativeIdx)%partialDerivativeIndex

    EXITS("DomainNodes_DerivativePartialIndexGet")
    RETURN
999 ERRORSEXITS("DomainNodes_DerivativePartialIndexGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNodes_DerivativePartialIndexGet

  !
  !================================================================================================================================
  !

  !>Gets the version number for a derivative of local node number from a domain. 
  SUBROUTINE DomainNodes_DerivativeVersionNumberGet(domainNodes,versionIdx,derivativeIdx,localNodeNumber,versionNumber,err,error,*)

    !Argument variables
    TYPE(DomainNodesType), POINTER :: domainNodes !<A pointer to the domain nodes to get the version number for
    INTEGER(INTG), INTENT(IN) :: versionIdx !<The version index to get the version number for
    INTEGER(INTG), INTENT(IN) :: derivativeIdx !<The derivative index to get the version number for
    INTEGER(INTG), INTENT(IN) :: localNodeNumber !<The local node number to get the version number for
    INTEGER(INTG), INTENT(OUT) :: versionNumber !<On exit, the version number for the derivative version index of the local node number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DomainNodeType), POINTER :: domainNode
    TYPE(DomainNodeDerivativeType), POINTER :: domainNodeDerivative
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainNodes_DerivativeVersionNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(domainNode)
    CALL DomainNodes_NodeGet(domainNodes,localNodeNumber,domainNode,err,error,*999)
    NULLIFY(domainNodeDerivative)
    CALL DomainNode_NodeDerivativeGet(domainNode,derivativeIdx,domainNodeDerivative,err,error,*999)
    IF(versionIdx<1.OR.versionIdx>domainNodeDerivative%numberOfVersions) THEN
      localError="The specified version index of "//TRIM(NumberToVString(versionIdx,"*",err,error))// &
        & " is invalid. The version index for derivative index "//TRIM(NumberToVString(derivativeIdx,"*",err,error))// &
        & " of local node number "//TRIM(NumberToVString(localNodeNumber,"*",err,error))//" should be >= 1 and <= "// &
        & TRIM(NumberToVString(domainNodeDerivative%numberOfVersions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainNodeDerivative%userVersionNumbers)) THEN
      localError="The user version numbers array is not allocated for derivative index "// &
        & TRIM(NumberToVString(derivativeIdx,"*",err,error))//" of local node number "// &
        & TRIM(NumberToVString(localNodeNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    versionNumber=domainNodes%nodes(localNodeNumber)%derivatives(derivativeIdx)%userVersionNumbers(versionIdx)

    EXITS("DomainNodes_DerivativeVersionNumberGet")
    RETURN
999 ERRORSEXITS("DomainNodes_DerivativeVersionNumberGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNodes_DerivativeVersionNumberGet

  !
  !================================================================================================================================
  !

  !>Gets a global node number that corresponds to a local node number from a domain. 
  SUBROUTINE DomainNodes_GlobalNodeNumberGet(domainNodes,localNodeNumber,globalNodeNumber,err,error,*)

    !Argument variables
    TYPE(DomainNodesType), POINTER :: domainNodes !<A pointer to the domain nodes to get the global node number for
    INTEGER(INTG), INTENT(IN) :: localNodeNumber !<The local node number to get the global node number for
    INTEGER(INTG), INTENT(OUT) :: globalNodeNumber !<On exit, the global node number corresponding to the local node number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DomainNodeType), POINTER :: domainNode
#endif    

    ENTERS("DomainNodes_GlobalNodeNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(domainNode)
    CALL DomainNodes_NodeGet(domainNodes,localNodeNumber,domainNode,err,error,*999)
#endif

    globalNodeNumber=domainNodes%nodes(localNodeNumber)%globalNumber

    EXITS("DomainNodes_GlobalNodeNumberGet")
    RETURN
999 ERRORSEXITS("DomainNodes_GlobalNodeNumberGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNodes_GlobalNodeNumberGet

  !
  !================================================================================================================================
  !

  !>Gets a local node number that corresponds to a user node number from a domain. An error will be raised if the user node number does not exist.
  SUBROUTINE DomainNodes_LocalNodeNumberGet(domainNodes,userNodeNumber,localNodeNumber,ghostNode,err,error,*)

    !Argument variables
    TYPE(DomainNodesType), POINTER :: domainNodes !<A pointer to the domain nodes to get the node on
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number to get
    INTEGER(INTG), INTENT(OUT) :: localNodeNumber !<On exit, the local number corresponding to the user node number.
    LOGICAL, INTENT(OUT) :: ghostNode !<On exit, is .TRUE. if the local node is a ghost node, .FALSE. if not.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    LOGICAL :: nodeExists
#ifdef WITH_POSTCHECKS    
    INTEGER(INTG) :: meshComponentNumber
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainNodes_LocalNodeNumberGet",err,error,*999)

    CALL DomainNodes_NodeCheckExists(domainNodes,userNodeNumber,nodeExists,localNodeNumber,ghostNode,err,error,*999)
    
#ifdef WITH_POSTCHECKS    
    IF(.NOT.nodeExists) THEN
      domainTopology=>domainNodes%domainTopology
      IF(ASSOCIATED(domainTopology)) THEN
        domain=>domainTopology%domain
        IF(ASSOCIATED(domain)) THEN
          meshComponentNumber=domain%meshComponentNumber
          decomposition=>domain%decomposition
          IF(ASSOCIATED(decomposition)) THEN
            localError="The user node number "//TRIM(NumberToVString(userNodeNumber,"*",err,error))// &
              & " does not exist in the domain from mesh component number "// &
              & TRIM(NumberToVString(meshComponentNumber,"*",err,error))//" of decomposition number "// &
              & TRIM(NumberToVString(decomposition%userNumber,"*",err,error))//"."
          ELSE
            localError="The user node number "//TRIM(NumberToVString(userNodeNumber,"*",err,error))// &
              & " does not exist in the domain from mesh component number "// &
              & TRIM(NumberToVString(meshComponentNumber,"*",err,error))//"."
          ENDIF
        ELSE
          localError="The user node number "//TRIM(NumberToVString(userNodeNumber,"*",err,error))//" does not exist."
        ENDIF
      ELSE
        localError="The user node number "//TRIM(NumberToVString(userNodeNumber,"*",err,error))//" does not exist."        
      ENDIF
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("DomainNodes_LocalNodeNumberGet")
    RETURN
999 ERRORSEXITS("DomainNodes_LocalNodeNumberGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNodes_LocalNodeNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the boundary node stauts for a local node number from a domain. 
  SUBROUTINE DomainNodes_NodeBoundaryNodeGet(domainNodes,localNodeNumber,boundaryNode,err,error,*)

    !Argument variables
    TYPE(DomainNodesType), POINTER :: domainNodes !<A pointer to the domain nodes to get the boundary node status for
    INTEGER(INTG), INTENT(IN) :: localNodeNumber !<The local node number to get the boundary node status for
    LOGICAL, INTENT(OUT) :: boundaryNode !<On exit, the boundary node flag for the local node number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DomainNodeType), POINTER :: domainNode
#endif    

    ENTERS("DomainNodes_NodeBoundaryNodeGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(domainNode)
    CALL DomainNodes_NodeGet(domainNodes,localNodeNumber,domainNode,err,error,*999)
#endif

    boundaryNode=domainNodes%nodes(localNodeNumber)%boundaryNode

    EXITS("DomainNodes_NodeBoundaryNodeGet")
    RETURN
999 ERRORSEXITS("DomainNodes_NodeBoundaryNodeGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainNodes_NodeBoundaryNodeGet
  
  !
  !================================================================================================================================
  !

  !>Checks that a user node number exists in domain nodes. 
  SUBROUTINE DomainNodes_NodeCheckExists(domainNodes,userNodeNumber,nodeExists,localNodeNumber,ghostNode,err,error,*)

    !Argument variables
    TYPE(DomainNodesType), POINTER :: domainNodes !<A pointer to the domain nodes to check the node exists on
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number to check if it exists
    LOGICAL, INTENT(OUT) :: nodeExists !<On exit, is .TRUE. if the node user number exists in the domain nodes topolgoy (even if it is a ghost node), .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: localNodeNumber !<On exit, if the node exists the local number corresponding to the user node number. If the node does not exist then global number will be 0.
    LOGICAL, INTENT(OUT) :: ghostNode !<On exit, is .TRUE. if the local node (if it exists) is a ghost node, .FALSE. if not. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(TreeNodeType), POINTER :: treeNode
    
    ENTERS("DomainNodes_NodeCheckExists",err,error,*999)

    nodeExists=.FALSE.
    localNodeNumber=0
    ghostNode=.FALSE.
    NULLIFY(treeNode)
    CALL Tree_Search(domainNodes%nodesTree,userNodeNumber,treeNode,err,error,*999)
    IF(ASSOCIATED(treeNode)) THEN
      CALL Tree_NodeValueGet(domainNodes%nodesTree,treeNode,localNodeNumber,err,error,*999)
      nodeExists=.TRUE.
      ghostNode=localNodeNumber>domainNodes%numberOfNodes
    ENDIF
    
    EXITS("DomainNodes_NodeCheckExists")
    RETURN
999 ERRORSEXITS("DomainNodes_NodeCheckExists",err,error)
    RETURN 1
    
  END SUBROUTINE DomainNodes_NodeCheckExists
  
  !
  !================================================================================================================================
  !

  !>Gets a face number for a local node number from a domain. 
  SUBROUTINE DomainNodes_NodeFaceNumberGet(domainNodes,faceIdx,localNodeNumber,localFaceNumber,err,error,*)

    !Argument variables
    TYPE(DomainNodesType), POINTER :: domainNodes !<A pointer to the domain nodes to get a face numberfor
    INTEGER(INTG), INTENT(IN) :: faceIdx !<The index of the face containing the node to get
    INTEGER(INTG), INTENT(IN) :: localNodeNumber !<The local node number containing the face to get
    INTEGER(INTG), INTENT(OUT) :: localFaceNumber !<On exit, the face number of the face containing the local node number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DomainNodeType), POINTER :: domainNode
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainNodes_NodeFaceNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(domainNode)
    CALL DomainNodes_NodeGet(domainNodes,localNodeNumber,domainNode,err,error,*999)
    IF(faceIdx<0.OR.faceIdx>domainNode%numberOfNodeFaces) THEN
      localError="The specified face index of "//TRIM(NumberToVString(faceIdx,"*",err,error))// &
        & " is invalid for local node number "//TRIM(NumberToVString(localNodeNumber,"*",err,error))// &
        & ". The face index should be >= 1 and <= "// &
        & TRIM(NumberToVString(domainNode%numberOfNodeFaces,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainNode%nodeFaces)) THEN
      localError="The node faces array is not allocated for local node number "// &
        & TRIM(NumberToVString(localNodeNumber,"*",err,error))//" of the domain nodes."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    localFaceNumber=domainNodes%nodes(localNodeNumber)%nodeFaces(faceIdx)

    EXITS("DomainNodes_NodeFaceNumberGet")
    RETURN
999 ERRORSEXITS("DomainNodes_NodeFaceNumberGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNodes_NodeFaceNumberGet

  !
  !================================================================================================================================
  !

  !>Gets a line number for a local node number from a domain. 
  SUBROUTINE DomainNodes_NodeLineNumberGet(domainNodes,lineIdx,localNodeNumber,localLineNumber,err,error,*)

    !Argument variables
    TYPE(DomainNodesType), POINTER :: domainNodes !<A pointer to the domain nodes to get a line numberfor
    INTEGER(INTG), INTENT(IN) :: lineIdx !<The index of the line containing the node to get
    INTEGER(INTG), INTENT(IN) :: localNodeNumber !<The local node number containing the line to get
    INTEGER(INTG), INTENT(OUT) :: localLineNumber !<On exit, the line number of the line containing the local node number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DomainNodeType), POINTER :: domainNode
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainNodes_NodeLineNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(domainNode)
    CALL DomainNodes_NodeGet(domainNodes,localNodeNumber,domainNode,err,error,*999)
    IF(lineIdx<0.OR.lineIdx>domainNode%numberOfNodeLines) THEN
      localError="The specified line index of "//TRIM(NumberToVString(lineIdx,"*",err,error))// &
        & " is invalid for local node number "//TRIM(NumberToVString(localNodeNumber,"*",err,error))// &
        & ". The line index should be >= 1 and <= "// &
        & TRIM(NumberToVString(domainNode%numberOfNodeLines,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainNode%nodeLines)) THEN
      localError="The node lines array is not allocated for local node number "// &
        & TRIM(NumberToVString(localNodeNumber,"*",err,error))//" of the domain nodes."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    localLineNumber=domainNodes%nodes(localNodeNumber)%nodeLines(lineIdx)

    EXITS("DomainNodes_NodeLineNumberGet")
    RETURN
999 ERRORSEXITS("DomainNodes_NodeLineNumberGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNodes_NodeLineNumberGet

  !
  !================================================================================================================================
  !

  !>Gets a node for a local node number in domain nodes 
  SUBROUTINE DomainNodes_NodeGet(domainNodes,localNodeNumber,domainNode,err,error,*)

    !Argument variables
    TYPE(DomainNodesType), POINTER :: domainNodes !<A pointer to the domain nodes to get the domain node for
    INTEGER(INTG), INTENT(IN) :: localNodeNumber !<The local node number to get the domain node for
    TYPE(DomainNodeType), POINTER, INTENT(OUT) :: domainNode !<On exit, a pointer to the specified domain node. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainNodes_NodeGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(domainNode)) CALL FlagError("Domain node is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainNodes)) CALL FlagError("Domain nodes is not associated.",err,error,*999)
    IF(localNodeNumber<1.OR.localNodeNumber>domainNodes%totalNumberOfNodes) THEN
      localError="The specified local node number of "//TRIM(NumberToVString(localNodeNumber,"*",err,error))// &
        & " is invalid. The local node number should be >= 1 and <= "// &
        & TRIM(NumberToVString(domainNodes%totalNumberOfNodes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainNodes%nodes)) CALL FlagError("The nodes array is not allocated for the domain nodes.",err,error,*999)
#endif

    domainNode=>domainNodes%nodes(localNodeNumber)

#ifdef WITH_POSTCHECKS
    IF(.NOT.ASSOCIATED(domainNode)) THEN
      localError="Domain node is not associated for local node number "//TRIM(NumberToVString(localNodeNumber,"*",err,error))// &
        & " of the domain nodes."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("DomainNodes_NodeGet")
    RETURN
999 NULLIFY(domainNode)
998 ERRORSEXITS("DomainNodes_NodeGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNodes_NodeGet

  !
  !================================================================================================================================
  !

  !>Gets the number of derivatives for a local node number from a domain. 
  SUBROUTINE DomainNodes_NodeNumberOfDerivativesGet(domainNodes,localNodeNumber,numberOfDerivatives,err,error,*)

    !Argument variables
    TYPE(DomainNodesType), POINTER :: domainNodes !<A pointer to the domain nodes to get the number of derivatives for
    INTEGER(INTG), INTENT(IN) :: localNodeNumber !<The local node number to get the number of derivatives for
    INTEGER(INTG), INTENT(OUT) :: numberOfDerivatives !<On exit, the number of derivatives for the local node number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DomainNodeType), POINTER :: domainNode
#endif    

    ENTERS("DomainNodes_NodeNumberOfDerivativesGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(domainNode)
    CALL DomainNodes_NodeGet(domainNodes,localNodeNumber,domainNode,err,error,*999)
#endif

    numberOfDerivatives=domainNodes%nodes(localNodeNumber)%numberOfDerivatives

    EXITS("DomainNodes_NodeNumberOfDerivativesGet")
    RETURN
999 ERRORSEXITS("DomainNodes_NodeNumberOfDerivativesGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNodes_NodeNumberOfDerivativesGet

  !
  !================================================================================================================================
  !

  !>Gets the number of faces for a local node number from a domain. 
  SUBROUTINE DomainNodes_NodeNumberOfFacesGet(domainNodes,localNodeNumber,numberOfFaces,err,error,*)

    !Argument variables
    TYPE(DomainNodesType), POINTER :: domainNodes !<A pointer to the domain nodes to get the number of node faces for
    INTEGER(INTG), INTENT(IN) :: localNodeNumber !<The local node number to get the number of faces for
    INTEGER(INTG), INTENT(OUT) :: numberOfFaces !<On exit, the number of faces for the local node number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DomainNodeType), POINTER :: domainNode
#endif    

    ENTERS("DomainNodes_NodeNumberOfFacesGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(domainNode)
    CALL DomainNodes_NodeGet(domainNodes,localNodeNumber,domainNode,err,error,*999)
#endif

    numberOfFaces=domainNodes%nodes(localNodeNumber)%numberOfNodeFaces

    EXITS("DomainNodes_NodeNumberOfFacesGet")
    RETURN
999 ERRORSEXITS("DomainNodes_NodeNumberOfFacesGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNodes_NodeNumberOfFacesGet

  !
  !================================================================================================================================
  !

  !>Gets the number of lines for a local node number from a domain. 
  SUBROUTINE DomainNodes_NodeNumberOfLinesGet(domainNodes,localNodeNumber,numberOfLines,err,error,*)

    !Argument variables
    TYPE(DomainNodesType), POINTER :: domainNodes !<A pointer to the domain nodes to get the number of node lines for
    INTEGER(INTG), INTENT(IN) :: localNodeNumber !<The local node number to get the number of lines for
    INTEGER(INTG), INTENT(OUT) :: numberOfLines !<On exit, the number of lines for the local node number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DomainNodeType), POINTER :: domainNode
#endif    

    ENTERS("DomainNodes_NodeNumberOfLinesGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(domainNode)
    CALL DomainNodes_NodeGet(domainNodes,localNodeNumber,domainNode,err,error,*999)
#endif

    numberOfLines=domainNodes%nodes(localNodeNumber)%numberOfNodeLines

    EXITS("DomainNodes_NodeNumberOfLinesGet")
    RETURN
999 ERRORSEXITS("DomainNodes_NodeNumberOfLinesGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNodes_NodeNumberOfLinesGet

  !
  !================================================================================================================================
  !

  !>Gets the number of surrounding elements for a local node number from a domain. 
  SUBROUTINE DomainNodes_NodeNumberOfSurroundingElementsGet(domainNodes,localNodeNumber,numberOfSurroundingElements,err,error,*)

    !Argument variables
    TYPE(DomainNodesType), POINTER :: domainNodes !<A pointer to the domain nodes to get the number of surrounding elements for
    INTEGER(INTG), INTENT(IN) :: localNodeNumber !<The local node number to get the number of surrounding elements for
    INTEGER(INTG), INTENT(OUT) :: numberOfSurroundingElements !<On exit, the number of surrounding elements for the local node number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DomainNodeType), POINTER :: domainNode
#endif    

    ENTERS("DomainNodes_NodeNumberOfSurroundingElementsGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(domainNode)
    CALL DomainNodes_NodeGet(domainNodes,localNodeNumber,domainNode,err,error,*999)
#endif

    numberOfSurroundingElements=domainNodes%nodes(localNodeNumber)%numberOfSurroundingElements

    EXITS("DomainNodes_NodeNumberOfSurroundingElementsGet")
    RETURN
999 ERRORS("DomainNodes_NodeNumberOfSurroundingElementsGet",err,error)
    EXITS("DomainNodes_NodeNumberOfSurroundingElementsGet")
    RETURN 1

  END SUBROUTINE DomainNodes_NodeNumberOfSurroundingElementsGet

  !
  !================================================================================================================================
  !

  !>Gets a surrounding elements for a local node number from a domain. 
  SUBROUTINE DomainNodes_NodeSurroundingElementGet(domainNodes,surroundingElementIdx,localNodeNumber,localElementNumber,err,error,*)

    !Argument variables
    TYPE(DomainNodesType), POINTER :: domainNodes !<A pointer to the domain nodes to get a surroundig element for
    INTEGER(INTG), INTENT(IN) :: surroundingElementIdx !<The index of the surrounding element to get
    INTEGER(INTG), INTENT(IN) :: localNodeNumber !<The local node number to get the surrounding element for
    INTEGER(INTG), INTENT(OUT) :: localElementNumber !<On exit, the surrounding element number of the surrounding element index and local node number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DomainNodeType), POINTER :: domainNode
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainNodes_NodeSurroundingElementGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(domainNode)
    CALL DomainNodes_NodeGet(domainNodes,localNodeNumber,domainNode,err,error,*999)
    IF(surroundingElementIdx<0.OR.surroundingElementIdx>domainNode%numberOfSurroundingElements) THEN
      localError="The specified surrounding element index of "//TRIM(NumberToVString(surroundingElementIdx,"*",err,error))// &
        & " is invalid for local node number "//TRIM(NumberToVString(localNodeNumber,"*",err,error))// &
        & ". The surrounding element index should be >= 1 and <= "// &
        & TRIM(NumberToVString(domainNode%numberOfSurroundingElements,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainNode%surroundingElements)) THEN
      localError="The surrounding elements array is not allocated for local node number "// &
        & TRIM(NumberToVString(localNodeNumber,"*",err,error))//" of the domain nodes."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    localElementNumber=domainNodes%nodes(localNodeNumber)%surroundingElements(surroundingElementIdx)

    EXITS("DomainNodes_NodeSurroundingElementGet")
    RETURN
999 ERRORSEXITS("DomainNodes_NodeSurroundingElementGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNodes_NodeSurroundingElementGet

  !
  !================================================================================================================================
  !

  !>Gets the number of nodes from a domain. 
  SUBROUTINE DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*)

    !Argument variables
    TYPE(DomainNodesType), POINTER :: domainNodes !<A pointer to the domain nodes to get the number of nodes for
    INTEGER(INTG), INTENT(OUT) :: numberOfNodes !<On exit, the number of nodes for the domain nodes.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainNodes_NumberOfNodesGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(domainNodes)) CALL FlagError("Domain nodes is not associated.",err,error,*999)
#endif

    numberOfNodes=domainNodes%numberOfNodes

    EXITS("DomainNodes_NumberOfNodesGet")
    RETURN
999 ERRORSEXITS("DomainNodes_NumberOfNodesGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNodes_NumberOfNodesGet

  !
  !================================================================================================================================
  !

  !>Gets the total number of nodes from a domain. 
  SUBROUTINE DomainNodes_TotalNumberOfNodesGet(domainNodes,totalNumberOfNodes,err,error,*)

    !Argument variables
    TYPE(DomainNodesType), POINTER :: domainNodes !<A pointer to the domain nodes to get the total number of nodes for
    INTEGER(INTG), INTENT(OUT) :: totalNumberOfNodes !<On exit, the total number of nodes for the domain nodes.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainNodes_TotalNumberOfNodesGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(domainNodes)) CALL FlagError("Domain nodes is not associated.",err,error,*999)
#endif

    totalNumberOfNodes=domainNodes%totalNumberOfNodes

    EXITS("DomainNodes_TotalNumberOfNodesGet")
    RETURN
999 ERRORSEXITS("DomainNodes_TotalNumberOfNodesGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNodes_TotalNumberOfNodesGet

  !
  !================================================================================================================================
  !

  !>Gets a user node number that corresponds to a local node number from a domain. 
  SUBROUTINE DomainNodes_UserNodeNumberGet(domainNodes,localNodeNumber,userNodeNumber,err,error,*)

    !Argument variables
    TYPE(DomainNodesType), POINTER :: domainNodes !<A pointer to the domain nodes to get the user node number for
    INTEGER(INTG), INTENT(IN) :: localNodeNumber !<The local node number to get the user node number for
    INTEGER(INTG), INTENT(OUT) :: userNodeNumber !<On exit, the user node number corresponding to the local node number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DomainNodeType), POINTER :: domainNode
#endif    

    ENTERS("DomainNodes_UserNodeNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS
    NULLIFY(domainNode)
    CALL DomainNodes_NodeGet(domainNodes,localNodeNumber,domainNode,err,error,*999)
#endif

    userNodeNumber=domainNodes%nodes(localNodeNumber)%userNumber

    EXITS("DomainNodes_UserNodeNumberGet")
    RETURN
999 ERRORSEXITS("DomainNodes_UserNodeNumberGet",err,error)
    RETURN 1

  END SUBROUTINE DomainNodes_UserNodeNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the domain from a domain topology.
  SUBROUTINE DomainTopology_DomainGet(domainTopology,domain,err,error,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: domainTopology !<A pointer to the domain topology to get the domain for.
    TYPE(DomainType), POINTER :: domain !<On exit, a pointer to the domain topology domain. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DomainTopology_DomainGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(domain)) CALL FlagError("Domain is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain topology is not associated.",err,error,*999)
#endif    

    !Get the domain
    domain=>domainTopology%domain

#ifdef WITH_POSTCHECKS    
    !Check domain is associated.
    IF(.NOT.ASSOCIATED(domain)) CALL FlagError("Domain topology domain is not associated.",err,error,*999)
#endif    
    
    EXITS("DomainTopology_DomainGet")
    RETURN
999 NULLIFY(domain)
998 ERRORSEXITS("DomainTopology_DomainGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainTopology_DomainGet

  !
  !================================================================================================================================
  !

  !>Gets dofs from a domain topology.
  SUBROUTINE DomainTopology_DomainDOFsGet(domainTopology,domainDofs,err,error,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: domainTopology !<A pointer to the domain topology to get the dofs for.
    TYPE(DomainDofsType), POINTER :: domainDofs !<On exit, a pointer to the domain topology dofs. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DomainTopology_DomainDOFsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(domainDofs)) CALL FlagError("Domain dofs is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain topology is not associated.",err,error,*999)
#endif    

    !Get the domain dofs
    domainDofs=>domainTopology%dofs

#ifdef WITH_POSTCHECKS    
    !Check domain dofs is associated.
    IF(.NOT.ASSOCIATED(domainDofs)) CALL FlagError("Domain topology dofs is not associated.",err,error,*999)
#endif    
    
    EXITS("DomainTopology_DomainDOFsGet")
    RETURN
999 NULLIFY(domainDofs)
998 ERRORSEXITS("DomainTopology_DomainDOFsGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainTopology_DomainDOFsGet

  !
  !================================================================================================================================
  !

  !>Gets elements from a domain topology.
  SUBROUTINE DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: domainTopology !<A pointer to the domain topology to get the elements for.
    TYPE(DomainElementsType), POINTER :: domainElements !<On exit, a pointer to the domain topology elements. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DomainTopology_DomainElementsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(domainElements)) CALL FlagError("Domain elements is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain topology is not associated.",err,error,*999)
#endif    

    !Get the domain elements
    domainElements=>domainTopology%elements

#ifdef WITH_POSTCHECKS    
    !Check domain elements is associated.
    IF(.NOT.ASSOCIATED(domainElements)) CALL FlagError("Domain topology elements is not associated.",err,error,*999)
#endif    
    
    EXITS("DomainTopology_DomainElementsGet")
    RETURN
999 ERRORSEXITS("DomainTopology_DomainElementsGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainTopology_DomainElementsGet

  !
  !================================================================================================================================
  !

  !>Gets lines from a domain topology.
  SUBROUTINE DomainTopology_DomainLinesGet(domainTopology,domainLines,err,error,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: domainTopology !<A pointer to the domain topology to get the lines for.
    TYPE(DomainLinesType), POINTER :: domainLines !<On exit, a pointer to the domain topology lines. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DomainTopology_DomainLinesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain is not associated.",err,error,*999)
    IF(ASSOCIATED(domainLines)) CALL FlagError("Domain lines is already associated.",err,error,*999)
#endif    

    !Get the domain lines
    domainLines=>domainTopology%lines

#ifdef WITH_POSTCHECKS    
    !Check domain lines is associated.
    IF(.NOT.ASSOCIATED(domainLines)) CALL FlagError("Domain topology lines is not associated.",err,error,*999)
#endif    
    
    EXITS("DomainTopology_DomainLinesGet")
    RETURN
999 ERRORSEXITS("DomainTopology_DomainLinesGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainTopology_DomainLinesGet

  !
  !================================================================================================================================
  !

  !>Gets nodes from a domain topology.
  SUBROUTINE DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: domainTopology !<A pointer to the domain topology to get the nodes for.
    TYPE(DomainNodesType), POINTER :: domainNodes !<On exit, a pointer to the domain topology nodes. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DomainTopology_DomainNodesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(domainNodes)) CALL FlagError("Domain lines is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain is not associated.",err,error,*999)
#endif    

    !Get the domain nodes
    domainNodes=>domainTopology%nodes

#ifdef WITH_POSTCHECKS    
    !Check domain nodes is associated.
    IF(.NOT.ASSOCIATED(domainNodes)) CALL FlagError("Domain topology nodes is not associated.",err,error,*999)
#endif    
    
    EXITS("DomainTopology_DomainNodesGet")
    RETURN
999 NULLIFY(domainNodes)
998 ERRORSEXITS("DomainTopology_DomainNodesGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainTopology_DomainNodesGet

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
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("DomainTopology_LocalElementBasisGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain topology is not associated.",err,error,*999)
    IF(localElementNumber<=0.OR.localElementNumber>domainTopology%elements%totalNumberOfElements) THEN
      localError="The local element number of "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " is invalid. The local number must be >= 1 and <= "// &
        & TRIM(NumberToVString(domainTopology%elements%totalNumberOfElements,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(domainTopology%elements)) CALL FlagError("Domain topology elements is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(domainTopology%elements%elements)) &
      & CALL FlagError("Domain topology elements elements is not associated.",err,error,*999)
#endif    
       
    basis=>domainTopology%elements%elements(localElementNumber)%basis

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(basis)) THEN
      localError="The basis for local element number "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("DomainTopology_LocalElementBasisGet")
    RETURN
999 ERRORSEXITS("DomainTopology_LocalElementBasisGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainTopology_LocalElementBasisGet

  !
  !================================================================================================================================
  !

  !>Gets faces from a domain topology.
  SUBROUTINE DomainTopology_DomainFacesGet(domainTopology,domainFaces,err,error,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: domainTopology !<A pointer to the domain topology to get the faces for.
    TYPE(DomainFacesType), POINTER :: domainFaces !<On exit, a pointer to the domain topology faces. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("DomainTopology_DomainFacesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain is not associated.",err,error,*999)
    IF(ASSOCIATED(domainFaces)) CALL FlagError("Domain faces is already associated.",err,error,*999)
#endif    

    !Get the domain faces
    domainFaces=>domainTopology%faces

#ifdef WITH_POSTCHECKS    
    !Check domain faces is associated.
    IF(.NOT.ASSOCIATED(domainFaces)) CALL FlagError("Domain topology faces is not associated.",err,error,*999)
#endif    
    
    EXITS("DomainTopology_DomainFacesGet")
    RETURN
999 ERRORSEXITS("DomainTopology_DomainFacesGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainTopology_DomainFacesGet

  !
  !================================================================================================================================
  !

END MODULE DecompositionAccessRoutines
