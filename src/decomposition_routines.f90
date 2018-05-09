!> \file
!> \author Chris Bradley
!> \brief This module handles all decomposition routines.
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

!>This module handles all decomposition routines.
MODULE DecompositionRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE CmissMPI
  USE CMISS_ParMETIS
  USE ComputationRoutines
  USE ComputationAccessRoutines
  USE Constants
  USE ContextAccessRoutines
  USE DecompositionAccessRoutines
  USE DOMAIN_MAPPINGS
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Kinds
  USE Lists
  USE MeshAccessRoutines
#ifndef NOMPIMOD
  USE MPI
#endif
  USE RegionAccessRoutines
  USE Strings
  USE Trees
  USE Types

#include "macros.h"  

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

   PRIVATE

  !Module parameters
  
  !Module types

  !Module variables

  !Interfaces
  
  INTERFACE DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS
    MODULE PROCEDURE DecompositionTopology_ElementCheckExists
  END INTERFACE DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS
  
  PUBLIC Decomposer_CreateStart,Decomposer_CreateFinish

  PUBLIC Decomposer_DecompositionAdd

  PUBLIC Decomposer_Destroy

  PUBLIC Decomposers_Initialise,Decomposers_Finalise

  PUBLIC DECOMPOSITIONS_INITIALISE,DECOMPOSITIONS_FINALISE

  PUBLIC DECOMPOSITION_CREATE_START,DECOMPOSITION_CREATE_FINISH

  PUBLIC DECOMPOSITION_DESTROY

  PUBLIC DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE

  PUBLIC DECOMPOSITION_ELEMENT_DOMAIN_GET,DECOMPOSITION_ELEMENT_DOMAIN_SET

  PUBLIC DECOMPOSITION_MESH_COMPONENT_NUMBER_GET,DECOMPOSITION_MESH_COMPONENT_NUMBER_SET
  
  PUBLIC DECOMPOSITION_NUMBER_OF_DOMAINS_GET,DECOMPOSITION_NUMBER_OF_DOMAINS_SET

  PUBLIC Decomposition_WorkGroupSet

  PUBLIC DecompositionTopology_DataPointCheckExists
  
  PUBLIC DecompositionTopology_DataProjectionCalculate

  PUBLIC DecompositionTopology_ElementCheckExists,DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS

  PUBLIC DecompositionTopology_ElementDataPointLocalNumberGet

  PUBLIC DecompositionTopology_ElementDataPointUserNumberGet

  PUBLIC DecompositionTopology_ElementGet

  PUBLIC DecompositionTopology_NumberOfElementDataPointsGet
  
  PUBLIC DECOMPOSITION_TYPE_GET,DECOMPOSITION_TYPE_SET
  
  PUBLIC DECOMPOSITION_NODE_DOMAIN_GET

  PUBLIC DECOMPOSITION_CALCULATE_LINES_SET,DECOMPOSITION_CALCULATE_FACES_SET

  PUBLIC DomainTopology_ElementBasisGet

  PUBLIC DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finish the creation of decomposer.
  SUBROUTINE Decomposer_CreateFinish(decomposer,err,error,*)

    !Argument variables
    TYPE(DecomposerType), POINTER :: decomposer !<A pointer to the decomposer to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Decomposer_CreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(decomposer)) CALL FlagError("Decomposer is not associated.",err,error,*999)
    IF(decomposer%decomposerFinished) CALL FlagError("Decomposer have already been finished.",err,error,*999)
    
    !Set the finished flag
    decomposer%decomposerFinished=.TRUE.
       
    EXITS("Decomposer_CreateFinish")
    RETURN
999 ERRORSEXITS("Decomposer_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_CreateFinish

  !
  !================================================================================================================================
  !

  !>Start the creation of decomposer for a region.
  SUBROUTINE Decomposer_CreateStart(decomposerUserNumber,region,workGroup,decomposer,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: decomposerUserNumber !<The user number of the decomposer to create
    TYPE(RegionType), POINTER :: region !<A pointer to the region to create decomposer for
    TYPE(WorkGroupType), POINTER :: workGroup !<A pointer to the work group for the created decomposer
    TYPE(DecomposerType), POINTER :: decomposer !<On exit, a pointer to the created decomposer. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposerIdx,dummyErr
    TYPE(DecomposerType), POINTER :: newDecomposer
    TYPE(DecomposerPtrType), ALLOCATABLE :: newDecomposers(:)
    TYPE(DecomposersType), POINTER :: decomposers
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("Decomposer_CreateStart",err,error,*998)

    IF(ASSOCIATED(decomposer)) CALL FlagError("Decomposer is already associated.",err,error,*998)
    CALL Region_AssertIsFinished(region,err,error,*998)
    CALL WorkGroup_AssertIsFinished(workGroup,err,error,*998)
    NULLIFY(decomposers)
    CALL Region_DecomposersGet(region,decomposers,err,error,*998)

    !Check if the user number has already been used.
    NULLIFY(newDecomposer)
    CALL Decomposer_UserNumberFind(decomposerUserNumber,region,newDecomposer,err,error,*998)
    IF(ASSOCIATED(newDecomposer)) THEN
      localError="Decomposer user number "//TRIM(NumberToVString(decomposerUserNumber,"*",err,error))// &
        & " has already been created on region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*998)
    ENDIF

    !Initialise the decomposer
    CALL Decomposer_Initialise(newDecomposer,err,error,*999)
    !Set default decomposer values
    newDecomposer%userNumber=decomposerUserNumber
    newDecomposer%decomposers=>decomposers
    newDecomposer%region=>region
    newDecomposer%workGroup=>workGroup
    !Add the decomposer to the list of decomposers for the region
    ALLOCATE(newDecomposers(decomposers%numberOfDecomposers+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new decomposers.",err,error,*999)
    DO decomposerIdx=1,decomposers%numberOfDecomposers
      newDecomposers(decomposerIdx)%ptr=>decomposers%decomposers(decomposerIdx)%ptr
    ENDDO !decomposerIdx
    newDecomposers(decomposers%numberOfDecomposers+1)%ptr=>newDecomposer
    CALL MOVE_ALLOC(newDecomposers,decomposers%decomposers)
    decomposers%numberOfDecomposers=decomposers%numberOfDecomposers+1
    !Return the pointer
    decomposer=>newDecomposer
       
    EXITS("Decomposer_CreateStart")
    RETURN
999 IF(ALLOCATED(newDecomposers)) DEALLOCATE(newDecomposers)
    CALL Decomposer_Finalise(newDecomposer,dummyErr,dummyError,*998)
998 ERRORSEXITS("Decomposer_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_CreateStart

  !
  !================================================================================================================================
  !

  !>Add a decomposition to a decomposer.
  SUBROUTINE Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err,error,*)

    !Argument variables
    TYPE(DecomposerType), POINTER :: decomposer !<A pointer to the decomposer to add decompositions to
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to add to the decomposer
    INTEGER(INTG), INTENT(OUT) :: decompositionIndex !<On return, the index of the decomposition that has been added to the decomposer
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: decompositionIdx
    LOGICAL :: isSubRegion
    TYPE(DecompositionPtrType), ALLOCATABLE :: newDecompositions(:)
    TYPE(RegionType), POINTER :: decomposerRegion,decompositionRegion
    TYPE(VARYING_STRING) :: localError

    ENTERS("Decomposer_DecompositionAdd",err,error,*999)

    CALL Decomposer_AssertNotFinished(decomposer,err,error,*999)
    CALL Decomposition_AssertIsFinished(decomposition,err,error,*999)

    NULLIFY(decomposerRegion)
    CALL Decomposer_RegionGet(decomposer,decomposerRegion,err,error,*999)
    NULLIFY(decompositionRegion)
    CALL Decomposition_RegionGet(decomposition,decompositionRegion,err,error,*999)
    CALL Region_IsSubRegion(decomposerRegion,decompositionRegion,isSubRegion,err,error,*999)
    IF(.NOT.isSubRegion) THEN
      localError="Decomposition number "//TRIM(NumberToVString(decomposition%userNumber,"*",err,error))// &
        & " is in region number "//TRIM(NumberToVString(decompositionRegion%userNumber,"*",err,error))// &
        & " which is not a sub-region of region number "//TRIM(NumberToVString(decomposerRegion%userNumber,"*",err,error))// &
        & " which has decomposer number "//TRIM(NumberToVString(decomposer%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(ASSOCIATED(decomposition%decomposer)) THEN
      localError="Decomposition number "//TRIM(NumberToVString(decomposition%userNumber,"*",err,error))// &
        & " already has an decomposer associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    ALLOCATE(newDecompositions(decomposer%numberOfDecompositions+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new decomposition.",err,error,*999)
    DO decompositionIdx=1,decomposer%numberOfDecompositions
      newDecompositions(decompositionIdx)%ptr=>decomposer%decompositions(decompositionIdx)%ptr
    ENDDO !decompositionIdx
    newDecompositions(decomposer%numberOfDecompositions+1)%ptr=>decomposition
    CALL MOVE_ALLOC(newDecompositions,decomposer%decompositions)
    decomposer%numberOfDecompositions=decomposer%numberOfDecompositions+1
    decompositionIndex=decomposer%numberOfDecompositions
       
    EXITS("Decomposer_DecompositionAdd")
    RETURN
999 ERRORSEXITS("Decomposer_DecompositionAdd",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_DecompositionAdd

  !
  !================================================================================================================================
  !

  !>Destroys a decomposer
  SUBROUTINE Decomposer_Destroy(decomposer,err,error,*)

    !Argument variables
    TYPE(DecomposerType), POINTER :: decomposer !<A pointer to the decomposer to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Decomposer_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(decomposer)) CALL FlagError("Decomposer is not associated.",err,error,*999)
    
    CALL Decomposer_Finalise(decomposer,err,error,*999)
       
    EXITS("Decomposer_Destroy")
    RETURN
999 ERRORSEXITS("Decomposer_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_Destroy

  !
  !================================================================================================================================
  !

  !>Finalise the decomposer and deallocate all memory.
  SUBROUTINE Decomposer_Finalise(decomposer,err,error,*)

    !Argument variables
    TYPE(DecomposerType), POINTER :: decomposer !<A pointer to the decomposer to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: decompositionIdx

    ENTERS("Decomposer_Finalise",err,error,*999)

    IF(ASSOCIATED(decomposer)) THEN
      IF(ALLOCATED(decomposer%decompositions)) THEN
        IF(decomposer%decomposerFinished) THEN
          DO decompositionIdx=1,SIZE(decomposer%decompositions,1)
            IF(ASSOCIATED(decomposer%decompositions(decompositionIdx)%ptr)) THEN
              NULLIFY(decomposer%decompositions(decompositionIdx)%ptr%decomposer)
            ENDIF
          ENDDO !decompositionIdx
        ENDIF
        DEALLOCATE(decomposer%decompositions)
      ENDIF
      DEALLOCATE(decomposer)
    ENDIF
       
    EXITS("Decomposer_Finalise")
    RETURN
999 ERRORSEXITS("Decomposer_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the decomposer for a region.
  SUBROUTINE Decomposer_Initialise(decomposer,err,error,*)

    !Argument variables
    TYPE(DecomposerType), POINTER :: decomposer !<A pointer to the decomposer to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("Decomposer_Initialise",err,error,*998)

    IF(ASSOCIATED(decomposer)) CALL FlagError("Decomposer is already associated.",err,error,*998)

    ALLOCATE(decomposer,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate decomposer.",err,error,*999)
    decomposer%userNumber=0
    NULLIFY(decomposer%decomposers)
    NULLIFY(decomposer%region)
    decomposer%decomposerFinished=.FALSE.
    decomposer%numberOfDecompositions=0
    NULLIFY(decomposer%workGroup)
        
    EXITS("Decomposer_Initialise")
    RETURN
999 CALL Decomposer_Finalise(decomposer,dummyErr,dummyError,*998)
998 ERRORSEXITS("Decomposer_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_Initialise
  
  !
  !================================================================================================================================
  !

  !>Finalise the decomposers and deallocate all memory.
  SUBROUTINE Decomposers_Finalise(decomposers,err,error,*)

    !Argument variables
    TYPE(DecomposersType), POINTER :: decomposers !<A pointer to the decomposers to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposerIdx

    ENTERS("Decomposers_Finalise",err,error,*999)

    IF(ASSOCIATED(decomposers)) THEN
      DO decomposerIdx=1,SIZE(decomposers%decomposers,1)
        CALL Decomposer_Finalise(decomposers%decomposers(decomposeridx)%ptr,err,error,*999)
      ENDDO !decomposerIdx
      DEALLOCATE(decomposers)
    ENDIF
       
    EXITS("Decomposers_Finalise")
    RETURN
999 ERRORSEXITS("Decomposers_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposers_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the decomposers for a region.
  SUBROUTINE Decomposers_Initialise(region,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to initialise the decomposers for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("Decomposers_Initialise",err,error,*998)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*998)
    IF(ASSOCIATED(region%decomposers)) CALL FlagError("Decomposers is already associated for this region.",err,error,*998)

    ALLOCATE(region%decomposers,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate decomposers.",err,error,*999)
    region%decomposers%region=>region
    region%decomposers%numberOfDecomposers=0
       
    EXITS("Decomposers_Initialise")
    RETURN
999 CALL Decomposers_Finalise(region%decomposers,dummyErr,dummyError,*998)
998 ERRORSEXITS("Decomposers_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposers_Initialise
  
  !
  !================================================================================================================================
  !

  !>Finialises a decomposition adjacent element information and deallocates all memory
  SUBROUTINE DECOMPOSITION_ADJACENT_ELEMENT_FINALISE(DECOMPOSITION_ADJACENT_ELEMENT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(DecompositionAdjacentElementType) :: DECOMPOSITION_ADJACENT_ELEMENT !<The decomposition adjacent element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_ADJACENT_ELEMENT_FINALISE",ERR,ERROR,*999)

    DECOMPOSITION_ADJACENT_ELEMENT%numberOfAdjacentElements=0
    IF(ALLOCATED(DECOMPOSITION_ADJACENT_ELEMENT%adjacentElements)) DEALLOCATE(DECOMPOSITION_ADJACENT_ELEMENT%adjacentElements)
       
    EXITS("DECOMPOSITION_ADJACENT_ELEMENT_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_ADJACENT_ELEMENT_FINALISE",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE DECOMPOSITION_ADJACENT_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !
  !>Initalises a decomposition adjacent element information.
  SUBROUTINE DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE(DECOMPOSITION_ADJACENT_ELEMENT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(DecompositionAdjacentElementType) :: DECOMPOSITION_ADJACENT_ELEMENT !<The decomposition adjacent element to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE",ERR,ERROR,*999)

    DECOMPOSITION_ADJACENT_ELEMENT%numberOfAdjacentElements=0
       
    EXITS("DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a domain decomposition on a given mesh. \see OpenCMISS::Iron::cmfe_Decomposition_CreateFinish
  SUBROUTINE DECOMPOSITION_CREATE_FINISH(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: DECOMPOSITION !<A pointer to the decomposition to finish creating
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposition_no
    TYPE(MeshType), POINTER :: MESH

    ENTERS("DECOMPOSITION_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      !Calculate which elements belong to which domain
      CALL DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE(DECOMPOSITION,ERR,ERROR,*999)
      !Initialise the topology information for this decomposition
      CALL DECOMPOSITION_TOPOLOGY_INITIALISE(DECOMPOSITION,ERR,ERROR,*999)
      !Initialise the domain for this computation node
      CALL DOMAIN_INITIALISE(DECOMPOSITION,ERR,ERROR,*999)
      !Calculate the decomposition topology
      CALL DECOMPOSITION_TOPOLOGY_CALCULATE(DECOMPOSITION,ERR,ERROR,*999)
      DECOMPOSITION%decompositionFinished=.TRUE.
      ! 
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      MESH=>DECOMPOSITION%MESH
      IF(ASSOCIATED(MESH)) THEN
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Mesh = ",MESH%userNumber,ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of decompositions = ", &
          & MESH%DECOMPOSITIONS%numberOfDecompositions,ERR,ERROR,*999)
        DO decomposition_no=1,MESH%DECOMPOSITIONS%numberOfDecompositions
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Decomposition number = ",decomposition_no,ERR,ERROR,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global number        = ", &
            & MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_no)%PTR%globalNumber,ERR,ERROR,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    User number          = ", &
            & MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_no)%PTR%userNumber,ERR,ERROR,*999)
        ENDDO !decomposition_no
      ELSE
        CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
      ENDIF
    ENDIF
    
    EXITS("DECOMPOSITION_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_CREATE_FINISH",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE DECOMPOSITION_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the creation of a domain decomposition for a given mesh. \see OpenCMISS::Iron::cmfe_Decomposition_CreateStart
  SUBROUTINE DECOMPOSITION_CREATE_START(USER_NUMBER,MESH,DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the decomposition 
    TYPE(MeshType), POINTER :: MESH !<A pointer to the mesh to decompose
    TYPE(DecompositionType), POINTER :: DECOMPOSITION !<On return a pointer to the created decomposition. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: decompositionIdx,dummyErr
    TYPE(ComputationEnvironmentType), POINTER :: computationEnvironment
    TYPE(ContextType), POINTER :: context
    TYPE(DecompositionType), POINTER :: newDecomposition
    TYPE(DecompositionPtrType), ALLOCATABLE :: newDecompositions(:)
    TYPE(RegionType), POINTER :: region    
    TYPE(VARYING_STRING) :: dummyError,LOCAL_ERROR
    TYPE(WorkGroupType), POINTER :: worldWorkGroup

    NULLIFY(newDecomposition)

    ENTERS("DECOMPOSITION_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(mesh)) THEN
      IF(mesh%meshFinished) THEN
        IF(ASSOCIATED(mesh%topology)) THEN
          IF(ASSOCIATED(mesh%decompositions)) THEN
            IF(ASSOCIATED(decomposition)) THEN
              CALL FlagError("Decomposition is already associated.",err,error,*999)
            ELSE
              NULLIFY(decomposition)
              CALL Decomposition_UserNumberFind(USER_NUMBER,MESH,decomposition,ERR,ERROR,*999)
              IF(ASSOCIATED(DECOMPOSITION)) THEN
                LOCAL_ERROR="Decomposition number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
                  & " has already been created on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%userNumber,"*",ERR,ERROR))//"."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              ELSE
                NULLIFY(newDecomposition)
                CALL Decomposition_Initialise(newDecomposition,err,error,*999)
                !Set default decomposition properties              
                newDecomposition%globalNumber=MESH%DECOMPOSITIONS%numberOfDecompositions+1
                newDecomposition%userNumber=USER_NUMBER
                newDecomposition%decompositions=>MESH%DECOMPOSITIONS
                newDecomposition%mesh=>mesh
                newDecomposition%region=>mesh%region
                newDecomposition%interface=>mesh%interface
                newDecomposition%numberOfDimensions=mesh%numberOfDimensions
                newDecomposition%numberOfComponents=mesh%NUMBER_OF_COMPONENTS
                !By default, the process of decompostion was done on the first mesh components. But the decomposition is the
                !same for all mesh components, since the decomposition is element-based.
                newDecomposition%meshComponentNumber=1
                !Default decomposition is all the mesh with one domain.
                newDecomposition%domainDecompositionType=DECOMPOSITION_ALL_TYPE
                NULLIFY(region)
                CALL Mesh_RegionGet(mesh,region,err,error,*999)
                NULLIFY(context)
                CALL Region_ContextGet(region,context,err,error,*999)
                NULLIFY(computationEnvironment)
                CALL Context_ComputationEnvironmentGet(context,computationEnvironment,err,error,*999)
                NULLIFY(worldWorkGroup)
                CALL ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err,error,*999)
                newDecomposition%workGroup=>worldWorkGroup
                newDecomposition%numberOfDomains=1
                newDecomposition%numberOfElements=mesh%numberOfElements
                ALLOCATE(newDecomposition%elementDomain(MESH%numberOfElements),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate new decomposition element domain.",ERR,ERROR,*999)
                newDecomposition%elementDomain=0          
                !\todo change this to use move alloc.
                !Add new decomposition into list of decompositions on the mesh
                ALLOCATE(newDecompositions(mesh%decompositions%numberOfDecompositions+1),STAT=err)
                IF(ERR/=0) CALL FlagError("Could not allocate new decompositions.",ERR,ERROR,*999)
                DO decompositionIdx=1,MESH%DECOMPOSITIONS%numberOfDecompositions
                  newDecompositions(decompositionIdx)%ptr=>mesh%decompositions%decompositions(decompositionIdx)%ptr
                ENDDO !decompositionIdx
                newDecompositions(mesh%decompositions%numberOfDecompositions+1)%ptr=>newDecomposition
                CALL MOVE_ALLOC(newDecompositions,mesh%decompositions%decompositions)
                mesh%decompositions%numberOfDecompositions=mesh%decompositions%numberOfDecompositions+1        
                decomposition=>newDecomposition
              ENDIF
            ENDIF
          ELSE
            LOCAL_ERROR="The decompositions on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%userNumber,"*",ERR,ERROR))// &
              & " are not associated."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Mesh topology is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_CREATE_START")
    RETURN
999 IF(ASSOCIATED(newDecomposition)) CALL Decomposition_Finalise(newDecomposition,dummyErr,dummyError,*998)
998 IF(ALLOCATED(newDecompositions)) DEALLOCATE(newDecompositions)
    ERRORSEXITS("DECOMPOSITION_CREATE_START",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE DECOMPOSITION_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroys a domain decomposition identified by a user number and deallocates all memory. \see OpenCMISS::Iron::cmfe_Decomposition_Destroy
  SUBROUTINE DECOMPOSITION_DESTROY_NUMBER(USER_NUMBER,MESH,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the decomposition to destroy.
    TYPE(MeshType), POINTER :: MESH !<A pointer to the mesh containing the decomposition.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposition_idx,decomposition_position
    LOGICAL :: FOUND    
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionPtrType), ALLOCATABLE :: newDecompositions(:)

    ENTERS("DECOMPOSITION_DESTROY_NUMBER",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN

        !Find the decomposition identified by the user number
        FOUND=.FALSE.
        decomposition_position=0
        DO WHILE(decomposition_position<MESH%DECOMPOSITIONS%numberOfDecompositions.AND..NOT.FOUND)
          decomposition_position=decomposition_position+1
          IF(MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_position)%PTR%userNumber==USER_NUMBER) FOUND=.TRUE.
        ENDDO
        
        IF(FOUND) THEN
          
          decomposition=>MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_position)%PTR

          !Finalise the decomposition
          CALL Decomposition_Finalise(decomposition,err,error,*999)

          !Remove the decomposition from the list of decompositions
          IF(MESH%DECOMPOSITIONS%numberOfDecompositions>1) THEN
            ALLOCATE(newDecompositions(MESH%DECOMPOSITIONS%numberOfDecompositions-1),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate new decompositions.",ERR,ERROR,*999)
            DO decomposition_idx=1,MESH%DECOMPOSITIONS%numberOfDecompositions
              IF(decomposition_idx<decomposition_position) THEN
                newDecompositions(decomposition_idx)%ptr=>mesh%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR
              ELSE IF(decomposition_idx>decomposition_position) THEN
                MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR%globalNumber= &
                  & MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR%globalNumber-1
                newDecompositions(decomposition_idx-1)%PTR=>MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR
              ENDIF
            ENDDO !decomposition_idx
            CALL MOVE_ALLOC(newDecompositions,mesh%decompositions%decompositions)
            MESH%DECOMPOSITIONS%numberOfDecompositions=MESH%DECOMPOSITIONS%numberOfDecompositions-1
          ELSE
            DEALLOCATE(MESH%DECOMPOSITIONS%DECOMPOSITIONS)
            MESH%DECOMPOSITIONS%numberOfDecompositions=0
          ENDIF
          
        ELSE
          LOCAL_ERROR="Decomposition number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has not been created on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%userNumber,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The decompositions on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%userNumber,"*",ERR,ERROR))// &
          & " are not associated."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_DESTROY_NUMBER")
    RETURN
999 IF(ALLOCATED(newDecompositions)) DEALLOCATE(newDecompositions)
    ERRORSEXITS("DECOMPOSITION_DESTROY_NUMBER",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE DECOMPOSITION_DESTROY_NUMBER

  !
  !================================================================================================================================
  !

  !>Destroys a domain decomposition identified by an object and deallocates all memory. \see OpenCMISS::Iron::cmfe_Decomposition_Destroy
  SUBROUTINE DECOMPOSITION_DESTROY(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: DECOMPOSITON !<The decomposition to destroy.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposition_idx,decomposition_position
    LOGICAL :: FOUND    
    TYPE(MeshType), POINTER :: MESH !
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionPtrType), ALLOCATABLE :: newDecompositions(:)
    TYPE(DecompositionsType), POINTER :: decompositions

    ENTERS("DECOMPOSITION_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      DECOMPOSITIONS=>DECOMPOSITION%DECOMPOSITIONS
      IF(ASSOCIATED(DECOMPOSITIONS)) THEN
        mesh=>decomposition%mesh
        IF(ASSOCIATED(mesh)) THEN
          !Find the decomposition identified by the user number
          FOUND=.FALSE.
          decomposition_position=0
          DO WHILE(decomposition_position<MESH%DECOMPOSITIONS%numberOfDecompositions.AND..NOT.FOUND)
            decomposition_position=decomposition_position+1
            IF(MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_position)%PTR%userNumber==DECOMPOSITION%userNumber) FOUND=.TRUE.
          ENDDO

          !Finalise the decomposition
          CALL Decomposition_Finalise(decomposition,err,error,*999)
          
          !Remove the decomposition from the list of decompositions
          IF(DECOMPOSITIONS%numberOfDecompositions>1) THEN
            ALLOCATE(newDecompositions(DECOMPOSITIONS%numberOfDecompositions-1),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate new decompositions.",ERR,ERROR,*999)
            DO decomposition_idx=1,DECOMPOSITIONS%numberOfDecompositions
              IF(decomposition_idx<decomposition_position) THEN
                newDecompositions(decomposition_idx)%PTR=>DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR
              ELSE IF(decomposition_idx>decomposition_position) THEN
                DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR%globalNumber= &
                  & DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR%globalNumber-1
                newDecompositions(decomposition_idx-1)%PTR=>DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR
              ENDIF
            ENDDO !decomposition_idx
            CALL MOVE_ALLOC(newDecompositions,DECOMPOSITIONS%DECOMPOSITIONS)
            DECOMPOSITIONS%numberOfDecompositions=DECOMPOSITIONS%numberOfDecompositions-1
          ELSE
            DEALLOCATE(DECOMPOSITIONS%DECOMPOSITIONS)
            DECOMPOSITIONS%numberOfDecompositions=0
          ENDIF
        ELSE
          CALL FlagError("Decomposition mesh is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition decompositions is not associated.",ERR,ERROR,*999)
       ENDIF
    ELSE
      CALL FlagError("Decompositions is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_DESTROY")
    RETURN
999 IF(ALLOCATED(newDecompositions)) DEALLOCATE(newDecompositions)
    ERRORSEXITS("DECOMPOSITION_DESTROY",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_DESTROY

  !
  !================================================================================================================================ !

  !>Calculates the element domains for a decomposition of a mesh. \see OpenCMISS::Iron::cmfe_Decomposition_ElementDomainCalculate
  SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: DECOMPOSITION !<A pointer to the decomposition to calculate the element domains for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: number_elem_indicies,elem_index,elem_count,ne,nn,myComputationNodeNumber,numberOfComputationNodes, &
      & no_computation_node,ELEMENT_START,ELEMENT_STOP,MY_ELEMENT_START,MY_ELEMENT_STOP,numberOfElements, &
      & MY_NUMBER_OF_ELEMENTS,MPI_IERROR,MAX_NUMBER_ELEMENTS_PER_NODE,component_idx,minNumberXi
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_COUNT(:),ELEMENT_PTR(:),ELEMENT_INDICIES(:),ELEMENT_DISTANCE(:),DISPLACEMENTS(:), &
      & RECEIVE_COUNTS(:)
    INTEGER(INTG) :: ELEMENT_WEIGHT(1),WEIGHT_FLAG,NUMBER_FLAG,NUMBER_OF_CONSTRAINTS, &
      & NUMBER_OF_COMMON_NODES,PARMETIS_OPTIONS(0:2),randomSeedsSize
    INTEGER(INTG) :: groupCommunicator
    INTEGER(INTG), ALLOCATABLE :: randomSeeds(:)
    REAL(DP) :: UBVEC(1)
    REAL(DP), ALLOCATABLE :: TPWGTS(:)
    REAL(DP) :: NUMBER_ELEMENTS_PER_NODE
    TYPE(ContextType), POINTER :: context    
    TYPE(BasisType), POINTER :: BASIS
    TYPE(MeshType), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
        MESH=>DECOMPOSITION%MESH
        IF(ASSOCIATED(MESH%TOPOLOGY)) THEN

          component_idx=decomposition%meshComponentNumber

          CALL WorkGroup_GroupCommunicatorGet(decomposition%workGroup,groupCommunicator,err,error,*999)
          CALL WorkGroup_NumberOfGroupNodesGet(decomposition%workGroup,numberOfComputationNodes,err,error,*999)
          CALL WorkGroup_GroupNodeNumberGet(decomposition%workGroup,myComputationNodeNumber,err,error,*999)
          
          SELECT CASE(DECOMPOSITION%domainDecompositionType)
          CASE(DECOMPOSITION_ALL_TYPE)
            !Do nothing. Decomposition checked below.
           CASE(DECOMPOSITION_CALCULATED_TYPE)
            !Calculate the general decomposition

            IF(decomposition%numberOfDomains==1) THEN
              DECOMPOSITION%elementDomain=0
            ELSE              
              NUMBER_ELEMENTS_PER_NODE=REAL(MESH%numberOfElements,DP)/REAL(numberOfComputationNodes,DP)
              ELEMENT_START=1
              ELEMENT_STOP=0
              MAX_NUMBER_ELEMENTS_PER_NODE=-1
              ALLOCATE(RECEIVE_COUNTS(0:numberOfComputationNodes-1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate recieve counts.",ERR,ERROR,*999)
              ALLOCATE(DISPLACEMENTS(0:numberOfComputationNodes-1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate displacements.",ERR,ERROR,*999)
              ALLOCATE(ELEMENT_DISTANCE(0:numberOfComputationNodes),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate element distance.",ERR,ERROR,*999)
              ELEMENT_DISTANCE(0)=0
              DO no_computation_node=0,numberOfComputationNodes-1
                ELEMENT_START=ELEMENT_STOP+1
                IF(no_computation_node==numberOfComputationNodes-1) THEN
                  ELEMENT_STOP=MESH%numberOfElements
                ELSE
                  ELEMENT_STOP=ELEMENT_START+NINT(NUMBER_ELEMENTS_PER_NODE,INTG)-1
                ENDIF
                IF((numberOfComputationNodes-1-no_computation_node)>(MESH%numberOfElements-ELEMENT_STOP)) &
                  & ELEMENT_STOP=MESH%numberOfElements-(numberOfComputationNodes-1-no_computation_node)
                IF(ELEMENT_START>MESH%numberOfElements) ELEMENT_START=MESH%numberOfElements
                IF(ELEMENT_STOP>MESH%numberOfElements) ELEMENT_STOP=MESH%numberOfElements
                DISPLACEMENTS(no_computation_node)=ELEMENT_START-1
                ELEMENT_DISTANCE(no_computation_node+1)=ELEMENT_STOP !C numbering
                numberOfElements=ELEMENT_STOP-ELEMENT_START+1
                RECEIVE_COUNTS(no_computation_node)=numberOfElements
                IF(numberOfElements>MAX_NUMBER_ELEMENTS_PER_NODE) MAX_NUMBER_ELEMENTS_PER_NODE=numberOfElements
                IF(no_computation_node==myComputationNodeNumber) THEN
                  MY_ELEMENT_START=ELEMENT_START
                  MY_ELEMENT_STOP=ELEMENT_STOP
                  MY_NUMBER_OF_ELEMENTS=ELEMENT_STOP-ELEMENT_START+1
                  number_elem_indicies=0
                  DO ne=MY_ELEMENT_START,MY_ELEMENT_STOP
                    BASIS=>MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)%BASIS
                    number_elem_indicies=number_elem_indicies+BASIS%numberOfNodes
                  ENDDO !ne
                ENDIF
              ENDDO !no_computation_node
              
              ALLOCATE(ELEMENT_PTR(0:MY_NUMBER_OF_ELEMENTS),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate element pointer list.",ERR,ERROR,*999)
              ALLOCATE(ELEMENT_INDICIES(0:number_elem_indicies-1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate element indicies list.",ERR,ERROR,*999)
              ALLOCATE(TPWGTS(1:decomposition%numberOfDomains),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate tpwgts.",ERR,ERROR,*999)
              elem_index=0
              elem_count=0
              ELEMENT_PTR(0)=0
              minNumberXi=99999
              DO ne=MY_ELEMENT_START,MY_ELEMENT_STOP
                elem_count=elem_count+1
                BASIS=>MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)%BASIS
                IF(BASIS%numberOfXi<minNumberXi) minNumberXi=BASIS%numberOfXi
                DO nn=1,BASIS%numberOfNodes
                  ELEMENT_INDICIES(elem_index)=MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)% &
                    & meshElementNodes(nn)-1 !C numbering
                  elem_index=elem_index+1
                ENDDO !nn
                ELEMENT_PTR(elem_count)=elem_index !C numbering
              ENDDO !ne
              
              !Set up ParMETIS variables
              NULLIFY(context)
              CALL WorkGroup_ContextGet(decomposition%workGroup,context,err,error,*999)
              CALL Context_RandomSeedsSizeGet(context,randomSeedsSize,err,error,*999)
              ALLOCATE(randomSeeds(randomSeedsSize),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate random seeds.",err,error,*999)
              CALL Context_RandomSeedsGet(context,randomSeeds,err,error,*999)              
              WEIGHT_FLAG=0 !No weights
              ELEMENT_WEIGHT(1)=1 !Isn't used due to weight flag
              NUMBER_FLAG=0 !C Numbering as there is a bug with Fortran numbering              
              NUMBER_OF_CONSTRAINTS=1
              IF(minNumberXi==1) THEN
                NUMBER_OF_COMMON_NODES=1
              ELSE
                NUMBER_OF_COMMON_NODES=2
              ENDIF
              !ParMETIS now has doule precision for these
             !TPWGTS=1.0_SP/REAL(decomposition%numberOfDomains,SP)
              !UBVEC=1.05_SP
              TPWGTS=1.0_DP/REAL(decomposition%numberOfDomains,DP)
              UBVEC=1.05_DP
              PARMETIS_OPTIONS(0)=1 !If zero, defaults are used, otherwise next two values are used
              PARMETIS_OPTIONS(1)=7 !Level of information to output
              PARMETIS_OPTIONS(2)=randomSeeds(1) !Seed for random number generator
              IF(ALLOCATED(randomSeeds)) DEALLOCATE(randomSeeds)
              !Call ParMETIS to calculate the partitioning of the mesh graph.
              CALL PARMETIS_PARTMESHKWAY(ELEMENT_DISTANCE,ELEMENT_PTR,ELEMENT_INDICIES,ELEMENT_WEIGHT,WEIGHT_FLAG,NUMBER_FLAG, &
                & NUMBER_OF_CONSTRAINTS,NUMBER_OF_COMMON_NODES,decomposition%numberOfDomains,TPWGTS,UBVEC,PARMETIS_OPTIONS, &
                & DECOMPOSITION%numberOfEdgesCut,DECOMPOSITION%elementDomain(DISPLACEMENTS(myComputationNodeNumber)+1:), &
                & groupCommunicator,ERR,ERROR,*999)
              
              !Transfer all the element domain information to the other computation nodes so that each rank has all the info
              IF(numberOfComputationNodes>1) THEN
                !This should work on a single processor but doesn't for mpich2 under windows. Maybe a bug? Avoid for now.
                CALL MPI_ALLGATHERV(MPI_IN_PLACE,MAX_NUMBER_ELEMENTS_PER_NODE,MPI_INTEGER,DECOMPOSITION%elementDomain, &
                  & RECEIVE_COUNTS,DISPLACEMENTS,MPI_INTEGER,groupCommunicator,MPI_IERROR)
                CALL MPI_ERROR_CHECK("MPI_ALLGATHERV",MPI_IERROR,ERR,ERROR,*999)
              ENDIF
              
              DEALLOCATE(DISPLACEMENTS)
              DEALLOCATE(RECEIVE_COUNTS)
              DEALLOCATE(ELEMENT_DISTANCE)
              DEALLOCATE(ELEMENT_PTR)
              DEALLOCATE(ELEMENT_INDICIES)
              DEALLOCATE(TPWGTS)

            ENDIF
            
          CASE(DECOMPOSITION_USER_DEFINED_TYPE)
            !Do nothing. Decomposition checked below.          
          CASE DEFAULT
            CALL FlagError("Invalid domain decomposition type.",ERR,ERROR,*999)            
          END SELECT

          !Check decomposition and check that each domain has an element in it.
          ALLOCATE(ELEMENT_COUNT(0:numberOfComputationNodes-1),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate element count.",ERR,ERROR,*999)
          ELEMENT_COUNT=0
          DO elem_index=1,MESH%numberOfElements
            no_computation_node=DECOMPOSITION%elementDomain(elem_index)
            IF(no_computation_node>=0.AND.no_computation_node<numberOfComputationNodes) THEN
              ELEMENT_COUNT(no_computation_node)=ELEMENT_COUNT(no_computation_node)+1
            ELSE
              LOCAL_ERROR="The computation node number of "//TRIM(NUMBER_TO_VSTRING(no_computation_node,"*",ERR,ERROR))// &
                & " for element number "//TRIM(NUMBER_TO_VSTRING(elem_index,"*",ERR,ERROR))// &
                & " is invalid. The computation node number must be between 0 and "// &
                & TRIM(NUMBER_TO_VSTRING(numberOfComputationNodes-1,"*",ERR,ERROR))//"."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDDO !elem_index
          DO no_computation_node=0,numberOfComputationNodes-1
            IF(ELEMENT_COUNT(no_computation_node)==0) THEN
              LOCAL_ERROR="Invalid decomposition. There are no elements in computation node "// &
                & TRIM(NUMBER_TO_VSTRING(no_computation_node,"*",ERR,ERROR))//"."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDDO !no_computation_node
          DEALLOCATE(ELEMENT_COUNT)
          
        ELSE
          CALL FlagError("Decomposition mesh topology is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Decomposition for mesh number ",DECOMPOSITION%MESH%userNumber, &
        & ERR,ERROR,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of domains = ", decomposition%numberOfDomains, &
        & ERR,ERROR,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Element domains:",ERR,ERROR,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Decomposition type = ", DECOMPOSITION%domainDecompositionType, &
        & ERR,ERROR,*999)
      IF(DECOMPOSITION%domainDecompositionType==DECOMPOSITION_CALCULATED_TYPE) THEN
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of edges cut = ",DECOMPOSITION%numberOfEdgesCut, &
          & ERR,ERROR,*999)
      ENDIF
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of elements = ",DECOMPOSITION%numberOfElements, &
        & ERR,ERROR,*999)
      DO ne=1,DECOMPOSITION%numberOfElements
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Element = ",ne,ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Domain = ",DECOMPOSITION%elementDomain(ne), &
          & ERR,ERROR,*999)
      ENDDO !ne
    ENDIF
    
    EXITS("DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE")
    RETURN
999 IF(ALLOCATED(RECEIVE_COUNTS)) DEALLOCATE(RECEIVE_COUNTS)
    IF(ALLOCATED(DISPLACEMENTS)) DEALLOCATE(DISPLACEMENTS)
    IF(ALLOCATED(ELEMENT_DISTANCE)) DEALLOCATE(ELEMENT_DISTANCE)
    IF(ALLOCATED(ELEMENT_PTR)) DEALLOCATE(ELEMENT_PTR)
    IF(ALLOCATED(ELEMENT_INDICIES)) DEALLOCATE(ELEMENT_INDICIES)
    IF(ALLOCATED(TPWGTS)) DEALLOCATE(TPWGTS)
    IF(ALLOCATED(randomSeeds)) DEALLOCATE(randomSeeds)
    ERRORSEXITS("DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Gets the domain for a given element in a decomposition of a mesh. \todo should be able to specify lists of elements. \see OpenCMISS::Iron::cmfe_Decomposition_ElementDomainGet
  SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_GET(DECOMPOSITION,USER_ELEMENT_NUMBER,DOMAIN_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the element domain for
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NUMBER !<The user element number to set the domain for.
    INTEGER(INTG), INTENT(OUT) :: DOMAIN_NUMBER !<On return, the domain of the global element.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables`
    TYPE(MeshType), POINTER :: MESH
    TYPE(MeshComponentTopologyType), POINTER :: MESH_TOPOLOGY
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: GLOBAL_ELEMENT_NUMBER
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    TYPE(MeshElementsType), POINTER :: MESH_ELEMENTS


    ENTERS("DECOMPOSITION_ELEMENT_DOMAIN_GET",ERR,ERROR,*999)

    GLOBAL_ELEMENT_NUMBER=0
    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%decompositionFinished) THEN
        MESH=>DECOMPOSITION%MESH
        IF(ASSOCIATED(MESH)) THEN
          MESH_TOPOLOGY=>MESH%TOPOLOGY(decomposition%meshComponentNumber)%PTR
          IF(ASSOCIATED(MESH_TOPOLOGY)) THEN
            MESH_ELEMENTS=>MESH_TOPOLOGY%ELEMENTS
            IF(ASSOCIATED(MESH_ELEMENTS)) THEN
              NULLIFY(TREE_NODE)
              CALL TREE_SEARCH(MESH_ELEMENTS%elementsTree,USER_ELEMENT_NUMBER,TREE_NODE,ERR,ERROR,*999)
              IF(ASSOCIATED(TREE_NODE)) THEN
                CALL TREE_NODE_VALUE_GET(MESH_ELEMENTS%elementsTree,TREE_NODE,GLOBAL_ELEMENT_NUMBER,ERR,ERROR,*999)
                IF(GLOBAL_ELEMENT_NUMBER>0.AND.GLOBAL_ELEMENT_NUMBER<=MESH_TOPOLOGY%ELEMENTS%numberOfElements) THEN
                  DOMAIN_NUMBER=DECOMPOSITION%elementDomain(GLOBAL_ELEMENT_NUMBER)
                ELSE
                  LOCAL_ERROR="Global element number found "//TRIM(NUMBER_TO_VSTRING(GLOBAL_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                    & " is invalid. The limits are 1 to "// &
                    & TRIM(NUMBER_TO_VSTRING(MESH_TOPOLOGY%ELEMENTS%numberOfElements,"*",ERR,ERROR))//"."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="Decomposition mesh element corresponding to user element number "// &
                  & TRIM(NumberToVString(USER_ELEMENT_NUMBER,"*",err,error))//" is not found."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Decomposition mesh elements are not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Decomposition mesh topology is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_ELEMENT_DOMAIN_GET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_ELEMENT_DOMAIN_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_GET
  
  !
  !================================================================================================================================
  !

  !>Sets the domain for a given element in a decomposition of a mesh. \todo move to user number, should be able to specify lists of elements. \see OpenCMISS::Iron::cmfe_Decomposition_ElementDomainSet 
  SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_SET(DECOMPOSITION,GLOBAL_ELEMENT_NUMBER,DOMAIN_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the element domain for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_ELEMENT_NUMBER !<The global element number to set the domain for.
    INTEGER(INTG), INTENT(IN) :: DOMAIN_NUMBER !<The domain of the global element.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfComputationNodes
    TYPE(MeshType), POINTER :: MESH
    TYPE(MeshComponentTopologyType), POINTER :: MESH_TOPOLOGY
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_ELEMENT_DOMAIN_SET",ERR,ERROR,*999)

!!TODO: interface should specify user element number ???
    
    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%decompositionFinished) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        MESH=>DECOMPOSITION%MESH
        IF(ASSOCIATED(MESH)) THEN
          MESH_TOPOLOGY=>MESH%TOPOLOGY(decomposition%meshComponentNumber)%PTR
          IF(ASSOCIATED(MESH_TOPOLOGY)) THEN
            IF(GLOBAL_ELEMENT_NUMBER>0.AND.GLOBAL_ELEMENT_NUMBER<=MESH_TOPOLOGY%ELEMENTS%numberOfElements) THEN
              CALL WorkGroup_NumberOfGroupNodesGet(decomposition%workGroup,numberOfComputationNodes,err,error,*999)
              IF(DOMAIN_NUMBER>=0.AND.DOMAIN_NUMBER<numberOfComputationNodes) THEN
                DECOMPOSITION%elementDomain(GLOBAL_ELEMENT_NUMBER)=DOMAIN_NUMBER
              ELSE
                LOCAL_ERROR="Domain number "//TRIM(NUMBER_TO_VSTRING(DOMAIN_NUMBER,"*",ERR,ERROR))// &
                  & " is invalid. The limits are 0 to "//TRIM(NUMBER_TO_VSTRING(numberOfComputationNodes,"*",ERR,ERROR))//"."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                & " is invalid. The limits are 1 to "// &
                & TRIM(NUMBER_TO_VSTRING(MESH_TOPOLOGY%ELEMENTS%numberOfElements,"*",ERR,ERROR))//"."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Decomposition mesh topology is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_ELEMENT_DOMAIN_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_ELEMENT_DOMAIN_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_SET
  
  !
  !================================================================================================================================
  !

  !>Finalises a decomposition and deallocates all memory.
  SUBROUTINE Decomposition_Finalise(decomposition,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Decomposition_Finalise",err,error,*999)

    
    IF(ASSOCIATED(decomposition)) THEN
      !Destroy all the decomposition components          
      IF(ALLOCATED(decomposition%elementDomain)) DEALLOCATE(decomposition%elementDomain)
      CALL DECOMPOSITION_TOPOLOGY_FINALISE(decomposition,err,error,*999)
      CALL DOMAIN_FINALISE(decomposition,err,error,*999)          
      DEALLOCATE(decomposition)
    ENDIF
    
    EXITS("Decomposition_Finalise")
    RETURN
999 ERRORSEXITS("Decomposition_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_Finalise
  
  !
  !================================================================================================================================
  !

  !>Inintialises a decomposition and allocates all memory.
  SUBROUTINE Decomposition_Initialise(decomposition,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to initialise. Must not be allocated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Decomposition_Initialise",err,error,*998)

    IF(ASSOCIATED(decomposition)) THEN
      CALL FlagError("Decomposition is already associated.",err,error,*998)
    ELSE
      ALLOCATE(decomposition,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate decomposition.",err,error,*999)
      decomposition%globalNumber=0
      decomposition%userNumber=0
      decomposition%decompositionFinished=.FALSE.
      NULLIFY(decomposition%decompositions)
      NULLIFY(decomposition%decomposer)
      NULLIFY(decomposition%mesh)
      NULLIFY(decomposition%region)
      NULLIFY(decomposition%INTERFACE)
      decomposition%numberOfDimensions=0
      decomposition%numberOfComponents=0
      decomposition%domainDecompositionType=DECOMPOSITION_ALL_TYPE
      NULLIFY(decomposition%workGroup)
      decomposition%numberOfDomains=0
      decomposition%numberOfEdgesCut=0
      decomposition%numberOfElements=0
      NULLIFY(decomposition%topology)
      NULLIFY(decomposition%domain)
      decomposition%calculateLines=.TRUE. 
      decomposition%calculateFaces=.TRUE.
    ENDIF
    
    EXITS("Decomposition_Initialise")
    RETURN
999 CALL Decomposition_Finalise(decomposition,dummyErr,dummyError,*998)
998 ERRORSEXITS("Decomposition_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_Initialise
  
  !
  !================================================================================================================================
  !

  !!MERGE: ditto
  
  !>Gets the mesh component number which will be used for the decomposition of a mesh. \see OpenCMISS::Iron::cmfe_DecompositionMeshComponentGet
  SUBROUTINE DECOMPOSITION_MESH_COMPONENT_NUMBER_GET(DECOMPOSITION,MESH_COMPONENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: DECOMPOSITION !<A pointer to the decomposition to get the mesh component for
    INTEGER(INTG), INTENT(OUT) :: MESH_COMPONENT_NUMBER !<On return, the mesh component number to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_MESH_COMPONENT_NUMBER_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN     
      IF(DECOMPOSITION%decompositionFinished) THEN
        IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
          MESH_COMPONENT_NUMBER=decomposition%meshComponentNumber
        ELSE
          CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_MESH_COMPONENT_NUMBER_GET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_MESH_COMPONENT_NUMBER_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_MESH_COMPONENT_NUMBER_GET
  
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the mesh component number which will be used for the decomposition of a mesh. \see OpenCMISS::Iron::cmfe_DecompositionMeshComponentSet
  SUBROUTINE DECOMPOSITION_MESH_COMPONENT_NUMBER_SET(DECOMPOSITION,MESH_COMPONENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the mesh component for
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component number to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_MESH_COMPONENT_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN     
      IF(DECOMPOSITION%decompositionFinished) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
          IF(MESH_COMPONENT_NUMBER>0.AND.MESH_COMPONENT_NUMBER<=DECOMPOSITION%numberOfComponents) THEN
            decomposition%meshComponentNumber=MESH_COMPONENT_NUMBER
          ELSE
            LOCAL_ERROR="The specified mesh component number of "//TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
              & "is invalid. The component number must be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(DECOMPOSITION%numberOfComponents,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_MESH_COMPONENT_NUMBER_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_MESH_COMPONENT_NUMBER_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_MESH_COMPONENT_NUMBER_SET
  
  !
  !================================================================================================================================
  !

  !>Gets the domain for a given node in a decomposition of a mesh. \todo should be able to specify lists of elements. \see OpenCMISS::Iron::cmfe_Decomposition_NodeDomainGet
  SUBROUTINE DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,USER_NODE_NUMBER,MESH_COMPONENT_NUMBER,DOMAIN_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the element domain for
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<The global element number to set the domain for.
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component number to set the domain for.
    INTEGER(INTG), INTENT(OUT) :: DOMAIN_NUMBER !<On return, the domain of the global element.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables`
    TYPE(MeshType), POINTER :: MESH
    TYPE(MeshComponentTopologyType), POINTER :: MESH_TOPOLOGY
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: GLOBAL_NODE_NUMBER
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    TYPE(MeshNodesType), POINTER :: MESH_NODES
    TYPE(DomainType), POINTER :: MESH_DOMAIN

    ENTERS("DECOMPOSITION_NODE_DOMAIN_GET",ERR,ERROR,*999)

!!TODO: interface should specify user element number ???
    GLOBAL_NODE_NUMBER=0
    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%decompositionFinished) THEN
        MESH=>DECOMPOSITION%MESH
        IF(ASSOCIATED(MESH)) THEN
          MESH_TOPOLOGY=>MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR
          IF(ASSOCIATED(MESH_TOPOLOGY)) THEN
            MESH_NODES=>MESH_TOPOLOGY%nodes
            IF(ASSOCIATED(MESH_NODES)) THEN
              NULLIFY(TREE_NODE)
              CALL TREE_SEARCH(MESH_NODES%nodesTree,USER_NODE_NUMBER,TREE_NODE,ERR,ERROR,*999)
              IF(ASSOCIATED(TREE_NODE)) THEN
                CALL TREE_NODE_VALUE_GET(MESH_NODES%nodesTree,TREE_NODE,GLOBAL_NODE_NUMBER,ERR,ERROR,*999)
                IF(GLOBAL_NODE_NUMBER>0.AND.GLOBAL_NODE_NUMBER<=MESH_TOPOLOGY%NODES%numberOfNodes) THEN
                  IF(MESH_COMPONENT_NUMBER>0.AND.MESH_COMPONENT_NUMBER<=MESH%NUMBER_OF_COMPONENTS) THEN
                    MESH_DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR
                    IF(ASSOCIATED(MESH_DOMAIN)) THEN
                      DOMAIN_NUMBER=MESH_DOMAIN%nodeDomain(GLOBAL_NODE_NUMBER)
                    ELSE
                      CALL FlagError("Decomposition domain is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="Mesh Component number "//TRIM(NumberToVString(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid. The limits are 1 to "// &
                      & TRIM(NumberToVString(MESH%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))//"."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="Global element number found "//TRIM(NumberToVString(GLOBAL_NODE_NUMBER,"*",ERR,ERROR))// &
                    & " is invalid. The limits are 1 to "// &
                    & TRIM(NumberToVString(MESH_TOPOLOGY%NODES%numberOfNodes,"*",ERR,ERROR))//"."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="Decomposition mesh node corresponding to user node number "// &
                  & TRIM(NumberToVString(USER_NODE_NUMBER,"*",err,error))//" is not found."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Decomposition mesh nodes are not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Decomposition mesh topology is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("DECOMPOSITION_NODE_DOMAIN_GET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_NODE_DOMAIN_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_NODE_DOMAIN_GET

  !
  !================================================================================================================================
  !

  !!MERGE: ditto
  
  !>Gets the number of domains for a decomposition. \see OpenCMISS::Iron::cmfe_DecompositionNumberOfDomainsGet
  SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_GET(DECOMPOSITION,numberOfDomains,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: DECOMPOSITION !<A pointer to the decomposition to get the number of domains for.
    INTEGER(INTG), INTENT(OUT) :: numberOfDomains !<On return, the number of domains to get.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_NUMBER_OF_DOMAINS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%decompositionFinished) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        numberOfDomains=decomposition%numberOfDomains
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_GET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_GET
  

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of domains for a decomposition. \see OpenCMISS::Iron::cmfe_DecompositionNumberOfDomainsSet
  SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_SET(DECOMPOSITION,numberOfDomains,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the number of domains for.
    INTEGER(INTG), INTENT(IN) :: numberOfDomains !<The number of domains to set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfGroupComputationNodes
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_NUMBER_OF_DOMAINS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%decompositionFinished) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(DECOMPOSITION%domainDecompositionType)
        CASE(DECOMPOSITION_ALL_TYPE)
          IF(numberOfDomains==1) THEN
            decomposition%numberOfDomains=1
          ELSE
            CALL FlagError("Can only have one domain for all decomposition type.",ERR,ERROR,*999)
          ENDIF
        CASE(DECOMPOSITION_CALCULATED_TYPE,DECOMPOSITION_USER_DEFINED_TYPE)
          IF(numberOfDomains>=1) THEN
            !wolfye???<=?
            IF(numberOfDomains<=DECOMPOSITION%numberOfElements) THEN
              !Get the number of computation nodes
              CALL WorkGroup_NumberOfGroupNodesGet(decomposition%workGroup,numberOfGroupComputationNodes, &
                & err,error,*999)
              !!TODO: relax this later
              !IF(numberOfDomains==numberOfGroupComputationNodes) THEN
                decomposition%numberOfDomains=numberOfDomains             
              !ELSE
              !  LOCAL_ERROR="The number of domains ("//TRIM(NUMBER_TO_VSTRING(numberOfDomainsS,"*",ERR,ERROR))// &
              !    & ") is not equal to the number of computation nodes ("// &
              !    & TRIM(NUMBER_TO_VSTRING(numberOfGroupComputationNodes,"*",ERR,ERROR))//")"
              !  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              !ENDIF
            ELSE
              LOCAL_ERROR="The number of domains ("//TRIM(NUMBER_TO_VSTRING(numberOfDomains,"*",ERR,ERROR))// &
                & ") must be <= the number of global elements ("// &
                & TRIM(NUMBER_TO_VSTRING(DECOMPOSITION%numberOfElements,"*",ERR,ERROR))//") in the mesh."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
             CALL FlagError("Number of domains must be >= 1.",ERR,ERROR,*999)
           ENDIF
         CASE DEFAULT
          LOCAL_ERROR="Decomposition type "//TRIM(NUMBER_TO_VSTRING(DECOMPOSITION%domainDecompositionType,"*",ERR,ERROR))// &
            & " is not valid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_SET

  !
  !================================================================================================================================
  !

  !>Sets the workgroup to use for a decomposition on a given mesh. \see OpenCMISS::Iron::cmfe_Decomposition_WorkGroupSet
  SUBROUTINE Decomposition_WorkGroupSet(decomposition,workGroup,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to set the work group for
    TYPE(WorkGroupType), POINTER :: workGroup !<A pointer to the work group to use for the decomposition
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposition_no
    TYPE(MeshType), POINTER :: MESH

    ENTERS("Decomposition_WorkGroupSet",err,error,*999)

    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
    IF(decomposition%decompositionFinished) CALL FlagError("Decomposition has already been finished.",err,error,*999)
    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("Work group is not associated.",err,error,*999)
    IF(.NOT.workGroup%workGroupFinished) CALL FlagError("Work group has not been finished.",err,error,*999)

    decomposition%workGroup=>workGroup
       
    EXITS("Decomposition_WorkGroupSet")
    RETURN
999 ERRORSEXITS("Decomposition_WorkGroupSet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_WorkGroupSet

  !
  !================================================================================================================================
  !

  !>Calculates the topology for a decomposition.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_CALCULATE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: DECOMPOSITION !<A pointer to the decomposition to finish creating.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: meshComponentNumber

    ENTERS("DECOMPOSITION_TOPOLOGY_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION%TOPOLOGY)) THEN
      !Calculate the elements topology
      CALL DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      !Calculate the line topology
      IF(DECOMPOSITION%calculateLines)THEN
        CALL DECOMPOSITION_TOPOLOGY_LINES_CALCULATE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
      !Calculate the face topology
      IF(DECOMPOSITION%calculateFaces) THEN
        CALL DECOMPOSITION_TOPOLOGY_FACES_CALCULATE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
      meshComponentNumber=decomposition%meshComponentNumber
      IF(ALLOCATED(DECOMPOSITION%MESH%TOPOLOGY(meshComponentNumber)%PTR%dataPoints%dataPoints)) THEN
          CALL DecompositionTopology_DataPointsCalculate(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        ENDIF   
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_CALCULATE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_CALCULATE
  

  !
  !================================================================================================================================
  !

  !>Calculates the decomposition element topology.
  SUBROUTINE DecompositionTopology_DataPointsCalculate(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: localElement,globalElement,dataPointIdx,localData,meshComponentNumber,groupCommunicator
    INTEGER(INTG) :: INSERT_STATUS,MPI_IERROR,numberOfComputationNodes,myGroupComputationNodeNumber,NUMBER_OF_GHOST_DATA, &
      & NUMBER_OF_LOCAL_DATA
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    TYPE(DomainMappingType), POINTER :: elementsMapping
    TYPE(MeshDataPointsType), POINTER :: meshData

    ENTERS("DecompositionTopology_DataPointsCalculate",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      decompositionData=>TOPOLOGY%dataPoints
      IF(ASSOCIATED(decompositionData)) THEN
        decomposition=>decompositionData%DECOMPOSITION
        IF(ASSOCIATED(decomposition)) THEN
          CALL WorkGroup_GroupCommunicatorGet(decomposition%workGroup,groupCommunicator,err,error,*999)
          decompositionElements=>TOPOLOGY%ELEMENTS
          IF(ASSOCIATED(decompositionElements)) THEN
            elementsMapping=>decomposition%DOMAIN(decomposition%meshComponentNumber)%PTR%MAPPINGS%ELEMENTS
            IF(ASSOCIATED(elementsMapping)) THEN
              meshComponentNumber=decomposition%meshComponentNumber
              meshData=>decomposition%MESH%TOPOLOGY(meshComponentNumber)%PTR%dataPoints
              IF(ASSOCIATED(meshData)) THEN
                CALL WorkGroup_NumberOfGroupNodesGet(decomposition%workGroup,numberOfComputationNodes,err,error,*999)
                CALL WorkGroup_GroupNodeNumberGet(decomposition%workGroup,myGroupComputationNodeNumber,err,error,*999)
                ALLOCATE(decompositionData%numberOfDomainLocal(0:numberOfComputationNodes-1),STAT=ERR)
                ALLOCATE(decompositionData%numberOfDomainGhost(0:numberOfComputationNodes-1),STAT=ERR)
                ALLOCATE(decompositionData%numberOfElementDataPoints(decompositionElements%numberOfGlobalElements),STAT=ERR)
                ALLOCATE(decompositionData%elementDataPoint(decompositionElements%totalNumberOfElements),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate decomposition element data points.",ERR,ERROR,*999)
                CALL TREE_CREATE_START(decompositionData%dataPointsTree,ERR,ERROR,*999)
                CALL TREE_INSERT_TYPE_SET(decompositionData%dataPointsTree,TREE_NO_DUPLICATES_ALLOWED,ERR,ERROR,*999)
                CALL TREE_CREATE_FINISH(decompositionData%dataPointsTree,ERR,ERROR,*999)
                decompositionData%numberOfGlobalDataPoints=meshData%totalNumberOfProjectedData 
                DO globalElement=1,decompositionElements%numberOfGlobalElements
                  decompositionData%numberOfElementDataPoints(globalElement)= &
                    & meshData%elementDataPoint(globalElement)%numberOfProjectedData
                ENDDO !globalElement
                localData=0;
                DO localElement=1,decompositionElements%totalNumberOfElements
                  globalElement=decompositionElements%ELEMENTS(localElement)%globalNumber
                  decompositionData%elementDataPoint(localElement)%numberOfProjectedData= &
                    & meshData%elementDataPoint(globalElement)%numberOfProjectedData
                  decompositionData%elementDataPoint(localElement)%globalElementNumber=globalElement
                  IF(localElement<elementsMapping%ghostStart) THEN
                    decompositionData%numberOfDataPoints=decompositionData%numberOfDataPoints+ &
                      & decompositionData%elementDataPoint(localElement)%numberOfProjectedData
                  ENDIF               
                  decompositionData%totalNumberOfDataPoints=decompositionData%totalNumberOfDataPoints+ &
                    & decompositionData%elementDataPoint(localElement)%numberOfProjectedData
                  ALLOCATE(decompositionData%elementDataPoint(localElement)%dataIndices(decompositionData% &
                    & elementDataPoint(localElement)%numberOfProjectedData),STAT=ERR)
                  DO dataPointIdx=1,decompositionData%elementDataPoint(localElement)%numberOfProjectedData
                    decompositionData%elementDataPoint(localElement)%dataIndices(dataPointIdx)%userNumber= &
                      & meshData%elementDataPoint(globalElement)%dataIndices(dataPointIdx)%userNumber
                    decompositionData%elementDataPoint(localElement)%dataIndices(dataPointIdx)%globalNumber= &
                      & meshData%elementDataPoint(globalElement)%dataIndices(dataPointIdx)%globalNumber
                    localData=localData+1
                    decompositionData%elementDataPoint(localElement)%dataIndices(dataPointIdx)%localNumber=localData
                    CALL TREE_ITEM_INSERT(decompositionData%dataPointsTree,decompositionData% &
                      & elementDataPoint(localElement)%dataIndices(dataPointIdx)%userNumber,localData, &
                      & INSERT_STATUS,ERR,ERROR,*999)
                  ENDDO !dataPointIdx
                ENDDO !localElement   
                !Calculate number of ghost data points on the current computation domain
                NUMBER_OF_LOCAL_DATA=decompositionData%numberOfDataPoints
                NUMBER_OF_GHOST_DATA=decompositionData%totalNumberOfDataPoints-decompositionData%numberOfDataPoints
                !Gather number of local data points on all computation nodes
                CALL MPI_ALLGATHER(NUMBER_OF_LOCAL_DATA,1,MPI_INTEGER,decompositionData% &
                  & numberOfDomainLocal,1,MPI_INTEGER,groupCommunicator,MPI_IERROR)
                CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)
                !Gather number of ghost data points on all computation nodes
                CALL MPI_ALLGATHER(NUMBER_OF_GHOST_DATA,1,MPI_INTEGER,decompositionData% &
                  & numberOfDomainGhost,1,MPI_INTEGER,groupCommunicator,MPI_IERROR)
                CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)
              ELSE
                CALL FlagError("Mesh data points topology is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Element mapping  is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Decomposition elements topology is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition data points topology is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DecompositionTopology_DataPointsCalculate")
    RETURN
999 ERRORSEXITS("DecompositionTopology_DataPointsCalculate",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_DataPointsCalculate

  !
  !================================================================================================================================
  !

  !>Calculates the decomposition element topology for a data projection (for data projections on fields).
  SUBROUTINE DecompositionTopology_DataProjectionCalculate(decompositionTopology,err,error,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("DecompositionTopology_DataProjectionCalculate",err,error,*999)

    IF(ASSOCIATED(decompositionTopology)) THEN
      CALL DecompositionTopology_DataPointsInitialise(decompositionTopology,err,error,*999)
      CALL DecompositionTopology_DataPointsCalculate(decompositionTopology,err,error,*999)
    ELSE
      CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DecompositionTopology_DataProjectionCalculate")
    RETURN
999 ERRORS("DecompositionTopology_DataProjectionCalculate",err,error)
    EXITS("DecompositionTopology_DataProjectionCalculate")
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_DataProjectionCalculate

  !
  !================================================================================================================================
  !

  !>Gets the local data point number for data points projected on an element
  SUBROUTINE DecompositionTopology_ElementDataPointLocalNumberGet(decompositionTopology,elementNumber,dataPointIndex, &
       & dataPointLocalNumber,err,error,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to get the data point for
    INTEGER(INTG), INTENT(IN) :: dataPointIndex !<The data point index to get the number for
    INTEGER(INTG), INTENT(OUT) :: dataPointLocalNumber !<The data point local number to return
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    INTEGER(INTG) :: numberOfDataPoints
    TYPE(VARYING_STRING) :: localError

    ENTERS("DecompositionTopology_ElementDataPointLocalNumberGet",err,error,*999)

    IF(ASSOCIATED(decompositionTopology)) THEN
      decompositionData=>decompositionTopology%dataPoints
      IF(ASSOCIATED(decompositionData)) THEN
        numberOfDataPoints = decompositionData%elementDataPoint(elementNumber)%numberOfProjectedData
        IF(dataPointIndex > 0 .AND. dataPointIndex <= numberOfDataPoints) THEN
          dataPointLocalNumber = decompositionData%elementDataPoint(elementNumber)%dataIndices(dataPointIndex)%localNumber
        ELSE
          localError="Element data point index "//TRIM(NUMBER_TO_VSTRING(dataPointIndex,"*",ERR,ERROR))// &
           & " out of range for element "//TRIM(NUMBER_TO_VSTRING(elementNumber,"*",ERR,ERROR))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition topology data points are not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DecompositionTopology_ElementDataPointLocalNumberGet")
    RETURN
999 ERRORS("DecompositionTopology_ElementDataPointLocalNumberGet",err,error)
    EXITS("DecompositionTopology_ElementDataPointLocalNumberGet")
    RETURN 1
  END SUBROUTINE DecompositionTopology_ElementDataPointLocalNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the user (global) data point number for data points projected on an element
  SUBROUTINE DecompositionTopology_ElementDataPointUserNumberGet(decompositionTopology,userElementNumber,dataPointIndex, &
       & dataPointUserNumber,err,error,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The element number to get the data point for
    INTEGER(INTG), INTENT(IN) :: dataPointIndex !<The data point index to get the number for
    INTEGER(INTG), INTENT(OUT) :: dataPointUserNumber !<The data point user number to return
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    INTEGER(INTG) :: numberOfDataPoints,decompositionLocalElementNumber
    LOGICAL :: ghostElement,userElementExists
    TYPE(VARYING_STRING) :: localError

    ENTERS("DecompositionTopology_ElementDataPointUserNumberGet",err,error,*999)

    IF(ASSOCIATED(decompositionTopology)) THEN
      decompositionData=>decompositionTopology%dataPoints
      IF(ASSOCIATED(decompositionData)) THEN
        CALL DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS(decompositionTopology,userElementNumber, &
          & userElementExists,decompositionLocalElementNumber,ghostElement,err,error,*999)
        IF(userElementExists) THEN
          IF(ghostElement) THEN
            localError="Cannot update by data point for user element "// &
              & TRIM(NUMBER_TO_VSTRING(userElementNumber,"*",err,error))//" as it is a ghost element."
            CALL FlagError(localError,err,error,*999)
          ELSE
            numberOfDataPoints = decompositionData%elementDataPoint(decompositionLocalElementNumber)%numberOfProjectedData
            IF(dataPointIndex > 0 .AND. dataPointIndex <= numberOfDataPoints) THEN
              dataPointUserNumber = decompositionData%elementDataPoint(decompositionLocalElementNumber)% &
                & dataIndices(dataPointIndex)%userNumber
            ELSE
              localError="Element data point index "//TRIM(NUMBER_TO_VSTRING(dataPointIndex,"*",ERR,ERROR))// &
               & " out of range for element "//TRIM(NUMBER_TO_VSTRING(userElementNumber,"*",ERR,ERROR))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDIF
        ELSE
          localError="The specified user element number of "// &
            & TRIM(NUMBER_TO_VSTRING(userElementNumber,"*",err,error))// &
            & " does not exist."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition topology data points are not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DecompositionTopology_ElementDataPointUserNumberGet")
    RETURN
999 ERRORS("DecompositionTopology_ElementDataPointUserNumberGet",err,error)
    EXITS("DecompositionTopology_ElementDataPointUserNumberGet")
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_ElementDataPointUserNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the number of data points projected on an element
  SUBROUTINE DecompositionTopology_NumberOfElementDataPointsGet(decompositionTopology,userElementNumber, &
       & numberOfDataPoints,err,error,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The element number to get the data point for
    INTEGER(INTG), INTENT(OUT) :: numberOfDataPoints !<The data point local number to return
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    INTEGER(INTG) :: decompositionLocalElementNumber
    LOGICAL :: ghostElement,userElementExists
    TYPE(VARYING_STRING) :: localError

    ENTERS("DecompositionTopology_NumberOfElementDataPointsGet",err,error,*999)

    IF(ASSOCIATED(decompositionTopology)) THEN
      decompositionData=>decompositionTopology%dataPoints
      IF(ASSOCIATED(decompositionData)) THEN
        CALL DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS(decompositionTopology,userElementNumber, &
          & userElementExists,decompositionLocalElementNumber,ghostElement,err,error,*999)
        IF(userElementExists) THEN
          IF(ghostElement) THEN
            localError="Cannot update by data point for user element "// &
              & TRIM(NUMBER_TO_VSTRING(userElementNumber,"*",err,error))//" as it is a ghost element."
            CALL FlagError(localError,err,error,*999)
          ELSE
            numberOfDataPoints = decompositionData%elementDataPoint(decompositionLocalElementNumber)%numberOfProjectedData
          ENDIF
        ELSE
          localError="The specified user element number of "// &
            & TRIM(NUMBER_TO_VSTRING(userElementNumber,"*",err,error))// &
            & " does not exist."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition topology data points are not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DecompositionTopology_NumberOfElementDataPointsGet")
    RETURN
999 ERRORS("DecompositionTopology_NumberOfElementDataPointsGet",err,error)
    EXITS("DecompositionTopology_NumberOfElementDataPointsGet")
    RETURN 1
  END SUBROUTINE DecompositionTopology_NumberOfElementDataPointsGet
  
  !
  !================================================================================================================================
  !

  !>Checks that a user element number exists in a decomposition. 
  SUBROUTINE DecompositionTopology_DataPointCheckExists(decompositionTopology,userDataPointNumber,userDataPointExists, &
        & decompositionLocalDataPointNumber,ghostDataPoint,err,error,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to check the data point exists on
    INTEGER(INTG), INTENT(IN) :: userDataPointNumber !<The user data point number to check if it exists
    LOGICAL, INTENT(OUT) :: userDataPointExists !<On exit, is .TRUE. if the data point user number exists in the decomposition topology, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: decompositionLocalDataPointNumber !<On exit, if the data point exists the local number corresponding to the user data point number. If the data point does not exist then local number will be 0.
    LOGICAL, INTENT(OUT) :: ghostDataPoint !<On exit, is .TRUE. if the local data point (if it exists) is a ghost data point, .FALSE. if not.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
    
    ENTERS("DecompositionTopology_DataPointCheckExists",ERR,error,*999)

    userDataPointExists=.FALSE.
    decompositionLocalDataPointNumber=0
    ghostDataPoint=.FALSE.
    IF(ASSOCIATED(decompositionTopology)) THEN
      decompositionData=>decompositionTopology%dataPoints
      IF(ASSOCIATED(decompositionData)) THEN
        NULLIFY(treeNode)
        CALL TREE_SEARCH(decompositionData%dataPointsTree,userDataPointNumber,treeNode,err,error,*999)
        IF(ASSOCIATED(treeNode)) THEN
          CALL TREE_NODE_VALUE_GET(decompositionData%dataPointsTree,treeNode,decompositionLocalDataPointNumber,err,error,*999)
          userDataPointExists=.TRUE.
          ghostDataPoint=decompositionLocalDataPointNumber>decompositionData%numberOfDataPoints
        ENDIF
      ELSE
        CALL FlagError("Decomposition data point topology is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DecompositionTopology_DataPointCheckExists")
    RETURN
999 ERRORS("DecompositionTopology_DataPointCheckExists",err,error)
    EXITS("DecompositionTopology_DataPointCheckExists")
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_DataPointCheckExists
  
  !
  !================================================================================================================================
  !

  !>Checks that a user element number exists in a decomposition. 
  SUBROUTINE DecompositionTopology_ElementCheckExists(decompositionTopology,userElementNumber,elementExists,localElementNumber, &
    & ghostElement,err,error,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to check the element exists on
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to check if it exists
    LOGICAL, INTENT(OUT) :: elementExists!<On exit, is .TRUE. if the element user number exists in the decomposition topology, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: localElementNumber !<On exit, if the element exists the local number corresponding to the user element number. If the element does not exist then local number will be 0.
    LOGICAL, INTENT(OUT) :: ghostElement !<On exit, is .TRUE. if the local element (if it exists) is a ghost element, .FALSE. if not.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
    
    ENTERS("DecompositionTopology_ElementCheckExists",err,error,*999)

    elementExists=.FALSE.
    localElementNumber=0
    ghostElement=.FALSE.
    IF(.NOT.ASSOCIATED(decompositionTopology)) CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_ElementsGet(decompositionTopology,decompositionElements,err,error,*999)
    NULLIFY(treeNode)
    CALL Tree_Search(decompositionElements%elementsTree,userElementNumber,treeNode,err,error,*999)
    IF(ASSOCIATED(treeNode)) THEN
      CALL Tree_NodeValueGet(decompositionElements%elementsTree,treeNode,localElementNumber,err,error,*999)
      elementExists=.TRUE.
      ghostElement=localElementNumber>decompositionElements%numberOfElements
    ENDIF
    
    EXITS("DecompositionTopology_ElementCheckExists")
    RETURN
999 ERRORSEXITS("DecompositionTopology_ElementCheckExists",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_ElementCheckExists

  !
  !================================================================================================================================
  !

!!TODO: Should this be decompositionElements???
  !>Gets a local element number that corresponds to a user element number from a decomposition. An error will be raised if the user element number does not exist.
  SUBROUTINE DecompositionTopology_ElementGet(decompositionTopology,userElementNumber,localElementNumber,ghostElement,err,error,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to get the element on
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to get
    INTEGER(INTG), INTENT(OUT) :: localElementNumber !<On exit, the local number corresponding to the user element number.
    LOGICAL, INTENT(OUT) :: ghostElement !<On exit, is .TRUE. if the local element is a ghost element, .FALSE. if not.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    LOGICAL :: elementExists
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DecompositionTopology_ElementGet",err,error,*999)

    CALL DecompositionTopology_ElementCheckExists(decompositionTopology,userElementNumber,elementExists,localElementNumber, &
      & ghostElement,err,error,*999)
    IF(.NOT.elementExists) THEN
      decomposition=>decompositionTopology%decomposition
      IF(ASSOCIATED(decomposition)) THEN
        localError="The user element number "//TRIM(NumberToVString(userElementNumber,"*",err,error))// &
          & " does not exist in decomposition number "//TRIM(NumberToVString(decomposition%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ELSE
        localError="The user element number "//TRIM(NumberToVString(userElementNumber,"*",err,error))//" does not exist."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF
   
    EXITS("DecompositionTopology_ElementGet")
    RETURN
999 ERRORSEXITS("DecompositionTopology_ElementGet",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_ElementGet

  !
  !================================================================================================================================
  !

  !>Finalises the given decomposition topology element.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE(ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionElementType) :: ELEMENT !<The decomposition element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nic !\todo add comment

    ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(ELEMENT%adjacentElements)) THEN
      DO nic=LBOUND(ELEMENT%adjacentElements,1),UBOUND(ELEMENT%adjacentElements,1)
        CALL DECOMPOSITION_ADJACENT_ELEMENT_FINALISE(ELEMENT%adjacentElements(nic),ERR,ERROR,*999)
      ENDDO !nic
      DEALLOCATE(ELEMENT%adjacentElements)
    ENDIF
    IF(ALLOCATED(ELEMENT%elementLines)) DEALLOCATE(ELEMENT%elementLines)
    IF(ALLOCATED(ELEMENT%elementFaces)) DEALLOCATE(ELEMENT%elementFaces)
 
    EXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given decomposition topology element.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE(ELEMENT,ERR,ERROR,*)
 
    !Argument variables
    TYPE(DecompositionElementType) :: ELEMENT !<The decomposition element to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR,*999)

    ELEMENT%userNumber=0
    ELEMENT%localNumber=0
    ELEMENT%globalNumber=0
    ELEMENT%boundaryElement=.FALSE.
  
    EXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Calculates the element numbers adjacent to an element in a decomposition topology.
  SUBROUTINE DecompositionTopology_ElementAdjacentElementCalculate(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the adjacent elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: j,ne,ne1,nep1,ni,nic,nn,nn1,nn2,nn3,np,np1,DUMMY_ERR,FACE_XI(2),FACE_XIC(3),nodePositionIndex(4)
    INTEGER(INTG) :: xi_direction,direction_index,xi_dir_check,xi_dir_search,NUMBER_NODE_MATCHES
    INTEGER(INTG) :: candidate_idx,face_node_idx,node_idx,surrounding_el_idx,candidate_el,idx
    INTEGER(INTG) :: NUMBER_SURROUNDING,NUMBER_OF_NODES_XIC(4),numberSurroundingElements
    INTEGER(INTG), ALLOCATABLE :: NODE_MATCHES(:),adjacentElements(:), surroundingElements(:)
    LOGICAL :: XI_COLLAPSED,FACE_COLLAPSED(-3:3),SUBSET
    TYPE(LIST_TYPE), POINTER :: NODE_MATCH_LIST, surroundingElementsList
    TYPE(LIST_PTR_TYPE) :: ADJACENT_ELEMENTS_LIST(-4:4)
    TYPE(BasisType), POINTER :: BASIS
    TYPE(DecompositionType), POINTER :: DECOMPOSITION
    TYPE(DecompositionElementsType), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(DomainType), POINTER :: DOMAIN
    TYPE(DomainElementsType), POINTER :: DOMAIN_ELEMENTS
    TYPE(DomainNodesType), POINTER :: DOMAIN_NODES
    TYPE(DomainTopologyType), POINTER :: DOMAIN_TOPOLOGY
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    NULLIFY(NODE_MATCH_LIST)
    DO nic=-4,4
      NULLIFY(ADJACENT_ELEMENTS_LIST(nic)%PTR)
    ENDDO !nic
    
    ENTERS("DecompositionTopology_ElementAdjacentElementCalculate",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
      IF(ASSOCIATED(DECOMPOSITION)) THEN
        DECOMPOSITION_ELEMENTS=>TOPOLOGY%ELEMENTS
        IF(ASSOCIATED(DECOMPOSITION_ELEMENTS)) THEN
          DOMAIN=>DECOMPOSITION%DOMAIN(decomposition%meshComponentNumber)%PTR
          IF(ASSOCIATED(DOMAIN)) THEN
            DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
            IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
              DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
              IF(ASSOCIATED(DOMAIN_NODES)) THEN
                DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                  !Loop over the elements in the decomposition
                  DO ne=1,DECOMPOSITION_ELEMENTS%totalNumberOfElements
                    BASIS=>DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS
                    !Create a list for every xi direction (plus and minus)
                    DO nic=-BASIS%numberOfXiCoordinates,BASIS%numberOfXiCoordinates
                      NULLIFY(ADJACENT_ELEMENTS_LIST(nic)%PTR)
                      CALL LIST_CREATE_START(ADJACENT_ELEMENTS_LIST(nic)%PTR,ERR,ERROR,*999)
                      CALL LIST_DATA_TYPE_SET(ADJACENT_ELEMENTS_LIST(nic)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                      CALL LIST_INITIAL_SIZE_SET(ADJACENT_ELEMENTS_LIST(nic)%PTR,5,ERR,ERROR,*999)
                      CALL LIST_CREATE_FINISH(ADJACENT_ELEMENTS_LIST(nic)%PTR,ERR,ERROR,*999)
                    ENDDO !nic
                    NUMBER_OF_NODES_XIC=1
                    NUMBER_OF_NODES_XIC(1:BASIS%numberOfXiCoordinates)= &
                      & BASIS%numberOfNodesXiC(1:BASIS%numberOfXiCoordinates)
                    !Place the current element in the surrounding list
                    CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(0)%PTR,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%localNumber, &
                      & ERR,ERROR,*999)
                    SELECT CASE(BASIS%TYPE)
                    CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
!!TODO: Calculate this and set it as part of the basis type
                      !Determine the collapsed "faces" if any
                      nodePositionIndex=1
                      !Loop over the face normals of the element
                      DO ni=1,BASIS%numberOfXi
                        !Determine the face xi directions that lie in this xi direction
                        FACE_XI(1)=OTHER_XI_DIRECTIONS3(ni,2,1)
                        FACE_XI(2)=OTHER_XI_DIRECTIONS3(ni,3,1)
                        !Reset the node_position_index in this xi direction
                        nodePositionIndex(ni)=1
                        !Loop over the two faces with this normal
                        DO direction_index=-1,1,2
                          xi_direction=direction_index*ni
                          FACE_COLLAPSED(xi_direction)=.FALSE.
                          DO j=1,2
                            xi_dir_check=FACE_XI(j)
                            IF(xi_dir_check<=BASIS%numberOfXi) THEN
                              xi_dir_search=FACE_XI(3-j)
                              nodePositionIndex(xi_dir_search)=1
                              XI_COLLAPSED=.TRUE.
                              DO WHILE(nodePositionIndex(xi_dir_search)<=NUMBER_OF_NODES_XIC(xi_dir_search).AND.XI_COLLAPSED)
                                !Get the first local node along the xi check direction
                                nodePositionIndex(xi_dir_check)=1
                                nn1=BASIS%nodePositionIndexInv(nodePositionIndex(1),nodePositionIndex(2), &
                                  & nodePositionIndex(3),1)
                                !Get the second local node along the xi check direction
                                nodePositionIndex(xi_dir_check)=2
                                nn2=BASIS%nodePositionIndexInv(nodePositionIndex(1),nodePositionIndex(2), &
                                  & nodePositionIndex(3),1)
                                IF(nn1/=0.AND.nn2/=0) THEN
                                  IF(DOMAIN_ELEMENTS%ELEMENTS(ne)%elementNodes(nn1)/= &
                                    & DOMAIN_ELEMENTS%ELEMENTS(ne)%elementNodes(nn2)) XI_COLLAPSED=.FALSE.
                                ENDIF
                                nodePositionIndex(xi_dir_search)=nodePositionIndex(xi_dir_search)+1
                              ENDDO !xi_dir_search
                              IF(XI_COLLAPSED) FACE_COLLAPSED(xi_direction)=.TRUE.
                            ENDIF
                          ENDDO !j
                          nodePositionIndex(ni)=NUMBER_OF_NODES_XIC(ni)
                        ENDDO !direction_index
                      ENDDO !ni
                      !Loop over the xi directions and calculate the surrounding elements
                      DO ni=1,BASIS%numberOfXi
                        !Determine the xi directions that lie in this xi direction
                        FACE_XI(1)=OTHER_XI_DIRECTIONS3(ni,2,1)
                        FACE_XI(2)=OTHER_XI_DIRECTIONS3(ni,3,1)
                        !Loop over the two faces
                        DO direction_index=-1,1,2
                          xi_direction=direction_index*ni                  
                          !Find nodes in the element on the appropriate face/line/point
                          NULLIFY(NODE_MATCH_LIST)
                          CALL LIST_CREATE_START(NODE_MATCH_LIST,ERR,ERROR,*999)
                          CALL LIST_DATA_TYPE_SET(NODE_MATCH_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                          CALL LIST_INITIAL_SIZE_SET(NODE_MATCH_LIST,16,ERR,ERROR,*999)
                          CALL LIST_CREATE_FINISH(NODE_MATCH_LIST,ERR,ERROR,*999)
                          IF(direction_index==-1) THEN
                            nodePositionIndex(ni)=1
                          ELSE
                            nodePositionIndex(ni)=NUMBER_OF_NODES_XIC(ni)
                          ENDIF
                          !If the face is collapsed then don't look in this xi direction. The exception is if the opposite face is
                          !also collapsed. This may indicate that we have a funny element in non-rc coordinates that goes around the
                          !central axis back to itself
                          IF(FACE_COLLAPSED(xi_direction).AND..NOT.FACE_COLLAPSED(-xi_direction)) THEN
                            !Do nothing - the match lists are already empty
                          ELSE
                            !Find the nodes to match and add them to the node match list
                            SELECT CASE(BASIS%numberOfXi)
                            CASE(1)
                              nn=BASIS%nodePositionIndexInv(nodePositionIndex(1),1,1,1)
                              IF(nn/=0) THEN
                                np=DOMAIN_ELEMENTS%ELEMENTS(ne)%elementNodes(nn)
                                CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                              ENDIF
                            CASE(2)
                              DO nn1=1,NUMBER_OF_NODES_XIC(FACE_XI(1)),NUMBER_OF_NODES_XIC(FACE_XI(1))-1
                                nodePositionIndex(FACE_XI(1))=nn1
                                nn=BASIS%nodePositionIndexInv(nodePositionIndex(1),nodePositionIndex(2),1,1)
                                IF(nn/=0) THEN
                                  np=DOMAIN_ELEMENTS%ELEMENTS(ne)%elementNodes(nn)
                                  CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                                ENDIF
                              ENDDO !nn1
                            CASE(3)
                              DO nn1=1,NUMBER_OF_NODES_XIC(FACE_XI(1)),NUMBER_OF_NODES_XIC(FACE_XI(1))-1
                                nodePositionIndex(FACE_XI(1))=nn1
                                DO nn2=1,NUMBER_OF_NODES_XIC(FACE_XI(2)),NUMBER_OF_NODES_XIC(FACE_XI(2))-1
                                  nodePositionIndex(FACE_XI(2))=nn2
                                  nn=BASIS%nodePositionIndexInv(nodePositionIndex(1),nodePositionIndex(2), &
                                    & nodePositionIndex(3),1)
                                  IF(nn/=0) THEN
                                    np=DOMAIN_ELEMENTS%ELEMENTS(ne)%elementNodes(nn)
                                    CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                                  ENDIF
                                ENDDO !nn2
                              ENDDO !nn1
                            CASE DEFAULT
                              LOCAL_ERROR="The number of xi directions in the basis of "// &
                                & TRIM(NUMBER_TO_VSTRING(BASIS%numberOfXi,"*",ERR,ERROR))//" is invalid."
                              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                            END SELECT
                          ENDIF
                          CALL LIST_REMOVE_DUPLICATES(NODE_MATCH_LIST,ERR,ERROR,*999)
                          CALL LIST_DETACH_AND_DESTROY(NODE_MATCH_LIST,NUMBER_NODE_MATCHES,NODE_MATCHES,ERR,ERROR,*999)
                          NUMBER_SURROUNDING=0
                          IF(NUMBER_NODE_MATCHES>0) THEN
                            !NODE_MATCHES now contain the list of corner nodes in the current face with normal_xi of ni.
                            !Look at the surrounding elements of each of these nodes, if there is a repeated element that
                            !is not the current element ne, it's an adjacent element.
                            candidate_idx=0
                            NULLIFY(surroundingElementsList)
                            CALL LIST_CREATE_START(surroundingElementsList,ERR,ERROR,*999)
                            CALL LIST_DATA_TYPE_SET(surroundingElementsList,LIST_INTG_TYPE,ERR,ERROR,*999)
                            CALL LIST_INITIAL_SIZE_SET(surroundingElementsList,2,ERR,ERROR,*999)
                            CALL LIST_CREATE_FINISH(surroundingElementsList,ERR,ERROR,*999)
                            DO face_node_idx=1,NUMBER_NODE_MATCHES
                              !Dump all the surrounding elements into an array, see if any are repeated
                              node_idx=NODE_MATCHES(face_node_idx)
                              DO surrounding_el_idx=1,DOMAIN_NODES%NODES(node_idx)%numberOfSurroundingElements
                                candidate_el=DOMAIN_NODES%NODES(node_idx)%surroundingElements(surrounding_el_idx)
                                IF(candidate_el/=ne) THEN
                                  candidate_idx=candidate_idx+1
                                  CALL LIST_ITEM_ADD(surroundingElementsList,candidate_el,ERR,ERROR,*999)
                                ENDIF
                              ENDDO
                            ENDDO !face_node_idx
                            CALL LIST_DETACH_AND_DESTROY(surroundingElementsList,numberSurroundingElements,surroundingElements, &
                              & ERR,ERROR,*999)
                            DO idx=1,candidate_idx
                              ne1=surroundingElements(idx)
                              IF(COUNT(surroundingElements(1:numberSurroundingElements)==ne1)>=BASIS%numberOfXi) THEN
                                !Found it, just exit
                                CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(xi_direction)%PTR,ne1,ERR,ERROR,*999)
                                NUMBER_SURROUNDING=NUMBER_SURROUNDING+1
                                EXIT
                              ENDIF
                            ENDDO
                          ENDIF
                          IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
                          IF(ALLOCATED(surroundingElements)) DEALLOCATE(surroundingElements)
                        ENDDO !direction_index
                      ENDDO !ni
                    CASE(BASIS_SIMPLEX_TYPE)
                      !Loop over the xi coordinates and calculate the surrounding elements
                      DO nic=1,BASIS%numberOfXiCoordinates
                        !Find the other coordinates of the face/line/point
                        FACE_XIC(1)=OTHER_XI_DIRECTIONS4(nic,1)
                        FACE_XIC(2)=OTHER_XI_DIRECTIONS4(nic,2)
                        FACE_XIC(3)=OTHER_XI_DIRECTIONS4(nic,3)
                        !Find nodes in the element on the appropriate face/line/point
                        NULLIFY(NODE_MATCH_LIST)
                        CALL LIST_CREATE_START(NODE_MATCH_LIST,ERR,ERROR,*999)
                        CALL LIST_DATA_TYPE_SET(NODE_MATCH_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                        CALL LIST_INITIAL_SIZE_SET(NODE_MATCH_LIST,16,ERR,ERROR,*999)
                        CALL LIST_CREATE_FINISH(NODE_MATCH_LIST,ERR,ERROR,*999)
                        nodePositionIndex(nic)=1 !Furtherest away from node with the nic'th coordinate
                        !Find the nodes to match and add them to the node match list
                        DO nn1=1,NUMBER_OF_NODES_XIC(FACE_XIC(1))
                          nodePositionIndex(FACE_XIC(1))=nn1
                          DO nn2=1,NUMBER_OF_NODES_XIC(FACE_XIC(2))
                            nodePositionIndex(FACE_XIC(2))=nn2
                            DO nn3=1,NUMBER_OF_NODES_XIC(FACE_XIC(3))
                              nodePositionIndex(FACE_XIC(3))=nn3
                              nn=BASIS%nodePositionIndexInv(nodePositionIndex(1),nodePositionIndex(2), &
                                & nodePositionIndex(3),nodePositionIndex(4))
                              IF(nn/=0) THEN
                                np=DOMAIN_ELEMENTS%ELEMENTS(ne)%elementNodes(nn)
                                CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                              ENDIF
                            ENDDO !nn3
                          ENDDO !nn2
                        ENDDO !nn1
                        CALL LIST_REMOVE_DUPLICATES(NODE_MATCH_LIST,ERR,ERROR,*999)
                        CALL LIST_DETACH_AND_DESTROY(NODE_MATCH_LIST,NUMBER_NODE_MATCHES,NODE_MATCHES,ERR,ERROR,*999)
                        IF(NUMBER_NODE_MATCHES>0) THEN
                          !Find list of elements surrounding those nodes
                          DO node_idx=1,NUMBER_NODE_MATCHES
                            np1=NODE_MATCHES(node_idx)
                            DO nep1=1,DOMAIN_NODES%NODES(np1)%numberOfSurroundingElements
                              ne1=DOMAIN_NODES%NODES(np1)%surroundingElements(nep1)
                              IF(ne1/=ne) THEN !Don't want the current element
                                ! grab the nodes list for current and this surrouding elements
                                ! current face : NODE_MATCHES
                                ! candidate elem : TOPOLOGY%ELEMENTS%ELEMENTS(ne1)%meshElementNodes 
                                ! if all of current face belongs to the candidate element, we will have found the neighbour
                                CALL LIST_SUBSET_OF(NODE_MATCHES(1:NUMBER_NODE_MATCHES),DOMAIN_ELEMENTS%ELEMENTS(ne1)% &
                                  & elementNodes,SUBSET,ERR,ERROR,*999)
                                IF(SUBSET) THEN
                                  CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(nic)%PTR,ne1,ERR,ERROR,*999)
                                ENDIF
                              ENDIF
                            ENDDO !nep1
                          ENDDO !node_idx
                        ENDIF
                        IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
                      ENDDO !nic
                    CASE(BASIS_SERENDIPITY_TYPE)
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    CASE(BASIS_AUXILLIARY_TYPE)
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    CASE(BASIS_B_SPLINE_TP_TYPE)
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    CASE(BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE)
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    CASE(BASIS_EXTENDED_LAGRANGE_TP_TYPE)
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The basis type of "//TRIM(NUMBER_TO_VSTRING(BASIS%TYPE,"*",ERR,ERROR))// &
                        & " is invalid."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)                     
                    END SELECT
                    !Set the surrounding elements for this element
                    ALLOCATE(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%adjacentElements(-BASIS%numberOfXiCoordinates: &
                      BASIS%numberOfXiCoordinates),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate adjacent elements.",ERR,ERROR,*999)
                    DO nic=-BASIS%numberOfXiCoordinates,BASIS%numberOfXiCoordinates
                      CALL DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%adjacentElements(nic), &
                        & ERR,ERROR,*999)
                      CALL LIST_DETACH_AND_DESTROY(ADJACENT_ELEMENTS_LIST(nic)%PTR,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)% &
                        & adjacentElements(nic)%numberOfAdjacentElements,adjacentElements,ERR,ERROR,*999)
                      ALLOCATE(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%adjacentElements(nic)%adjacentElements( &
                        DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%adjacentElements(nic)%numberOfAdjacentElements),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate element adjacent elements.",ERR,ERROR,*999)
                      DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%adjacentElements(nic)%adjacentElements(1: &
                        & DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%adjacentElements(nic)%numberOfAdjacentElements)= &
                        adjacentElements(1:DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%adjacentElements(nic)%numberOfAdjacentElements)
                      IF(ALLOCATED(adjacentElements)) DEALLOCATE(adjacentElements)
                    ENDDO !nic
                  ENDDO !ne           
                ELSE
                  CALL FlagError("Domain topology elements is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Domain topology nodes is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Topology decomposition domain topology is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Topology decomposition domain is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Topology elements is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Topology decomposition is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Topology is not allocated.",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Total number of elements = ",DECOMPOSITION_ELEMENTS% &
        & totalNumberOfElements,ERR,ERROR,*999)
      DO ne=1,DECOMPOSITION_ELEMENTS%totalNumberOfElements
        BASIS=>DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Local element number : ",ne,ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of xi coordinates = ",BASIS%numberOfXiCoordinates, &
          & ERR,ERROR,*999)
        DO nic=-BASIS%numberOfXiCoordinates,BASIS%numberOfXiCoordinates
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Xi coordinate : ",nic,ERR,ERROR,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of adjacent elements = ", &
            & DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%adjacentElements(nic)%numberOfAdjacentElements,ERR,ERROR,*999)
          IF(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%adjacentElements(nic)%numberOfAdjacentElements>0) THEN
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)% &
              & adjacentElements(nic)%numberOfAdjacentElements,8,8,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%adjacentElements(nic)% &
              & adjacentElements,'("        Adjacent elements :",8(X,I6))','(30x,8(X,I6))',ERR,ERROR,*999)
          ENDIF
        ENDDO !nic
      ENDDO !ne
    ENDIF

    EXITS("DecompositionTopology_ElementAdjacentElementCalculate")
    RETURN
999 IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
    IF(ALLOCATED(surroundingElements)) DEALLOCATE(surroundingElements)
    IF(ALLOCATED(adjacentElements)) DEALLOCATE(adjacentElements)
    IF(ASSOCIATED(NODE_MATCH_LIST)) CALL LIST_DESTROY(NODE_MATCH_LIST,DUMMY_ERR,DUMMY_ERROR,*998)
998 IF(ASSOCIATED(surroundingElementsList)) CALL LIST_DESTROY(surroundingElementsList,DUMMY_ERR,DUMMY_ERROR,*997)
997 DO nic=-4,4
      IF(ASSOCIATED(ADJACENT_ELEMENTS_LIST(nic)%PTR)) CALL LIST_DESTROY(ADJACENT_ELEMENTS_LIST(nic)%PTR,DUMMY_ERR,DUMMY_ERROR,*996)
    ENDDO !ni
996 ERRORS("DecompositionTopology_ElementAdjacentElementCalculate",ERR,ERROR)
    EXITS("DecompositionTopology_ElementAdjacentElementCalculate")
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_ElementAdjacentElementCalculate
  
  !
  !================================================================================================================================
  !

  !>Calculates the decomposition element topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: global_element,INSERT_STATUS,local_element
    TYPE(DecompositionType), POINTER :: DECOMPOSITION
    TYPE(DecompositionElementsType), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(DomainType), POINTER :: DOMAIN
    TYPE(DomainElementsType), POINTER :: DOMAIN_ELEMENTS
    TYPE(DomainMappingType), POINTER :: DOMAIN_ELEMENTS_MAPPING
    TYPE(DomainMappingsType), POINTER :: DOMAIN_MAPPINGS
    TYPE(DomainTopologyType), POINTER :: DOMAIN_TOPOLOGY
    TYPE(MeshType), POINTER :: MESH
    TYPE(MeshElementsType), POINTER :: MESH_ELEMENTS
    TYPE(MeshComponentTopologyType), POINTER :: MESH_TOPOLOGY

    ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      DECOMPOSITION_ELEMENTS=>TOPOLOGY%ELEMENTS
      IF(ASSOCIATED(DECOMPOSITION_ELEMENTS)) THEN
        DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
        IF(ASSOCIATED(DECOMPOSITION)) THEN
          DOMAIN=>DECOMPOSITION%DOMAIN(decomposition%meshComponentNumber)%PTR
          IF(ASSOCIATED(DOMAIN)) THEN
            DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
            IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
              DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
              IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
                IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
                  DOMAIN_ELEMENTS_MAPPING=>DOMAIN_MAPPINGS%ELEMENTS
                  IF(ASSOCIATED(DOMAIN_ELEMENTS_MAPPING)) THEN
                    MESH=>DECOMPOSITION%MESH
                    IF(ASSOCIATED(MESH)) THEN
                      MESH_TOPOLOGY=>MESH%TOPOLOGY(decomposition%meshComponentNumber)%PTR
                      IF(ASSOCIATED(MESH_TOPOLOGY)) THEN
                        MESH_ELEMENTS=>MESH_TOPOLOGY%ELEMENTS
                        IF(ASSOCIATED(MESH_ELEMENTS)) THEN
                          !Allocate the element topology arrays
                          ALLOCATE(DECOMPOSITION_ELEMENTS%ELEMENTS(DOMAIN_ELEMENTS%totalNumberOfElements),STAT=ERR)
                          IF(ERR/=0) CALL FlagError("Could not allocate decomposition elements elements.",ERR,ERROR,*999)
                          DECOMPOSITION_ELEMENTS%numberOfElements=DOMAIN_ELEMENTS%numberOfElements
                          DECOMPOSITION_ELEMENTS%totalNumberOfElements=DOMAIN_ELEMENTS%totalNumberOfElements
                          DECOMPOSITION_ELEMENTS%numberOfGlobalElements=DOMAIN_ELEMENTS%numberOfGlobalElements
                          CALL TREE_CREATE_START(DECOMPOSITION_ELEMENTS%elementsTree,ERR,ERROR,*999)
                          CALL TREE_INSERT_TYPE_SET(DECOMPOSITION_ELEMENTS%elementsTree,TREE_NO_DUPLICATES_ALLOWED,ERR,ERROR,*999)
                          CALL TREE_CREATE_FINISH(DECOMPOSITION_ELEMENTS%elementsTree,ERR,ERROR,*999)
                          DO local_element=1,DECOMPOSITION_ELEMENTS%totalNumberOfElements
                            CALL DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE(DECOMPOSITION_ELEMENTS%ELEMENTS(local_element), &
                              & ERR,ERROR,*999)
                            global_element=DOMAIN_ELEMENTS_MAPPING%localToGlobalMap(local_element)
                            DECOMPOSITION_ELEMENTS%ELEMENTS(local_element)%userNumber=MESH_ELEMENTS%ELEMENTS(global_element)% &
                              & userNumber
                            DECOMPOSITION_ELEMENTS%ELEMENTS(local_element)%localNumber=local_element
                            CALL TREE_ITEM_INSERT(DECOMPOSITION_ELEMENTS%elementsTree,DECOMPOSITION_ELEMENTS% &
                              & ELEMENTS(local_element)%userNumber,local_element,INSERT_STATUS,ERR,ERROR,*999)
                            DECOMPOSITION_ELEMENTS%ELEMENTS(local_element)%globalNumber=global_element
                            DECOMPOSITION_ELEMENTS%ELEMENTS(local_element)%boundaryElement=MESH_ELEMENTS% &
                              & ELEMENTS(global_element)%boundaryElement
                          ENDDO !local_element
                          !Calculate the elements surrounding the elements in the decomposition topology
                          CALL DecompositionTopology_ElementAdjacentElementCalculate(TOPOLOGY,ERR,ERROR,*999)
                        ELSE
                          CALL FlagError("Mesh elements is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Mesh topology is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Domain mappings elements is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Domain mappings is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Domain topology elements is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Topology decomposition domain topology is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Topology decomposition domain is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Topology decomposition is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Topology elements is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Finalises the elements in the given decomposition topology. \todo Pass in the decomposition elements pointer.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to finalise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne

    ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        DO ne=1,TOPOLOGY%ELEMENTS%totalNumberOfElements
          CALL DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE(TOPOLOGY%ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
        ENDDO !ne
        IF(ASSOCIATED(TOPOLOGY%ELEMENTS%ELEMENTS)) DEALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS)
        IF(ASSOCIATED(TOPOLOGY%ELEMENTS%elementsTree)) CALL TREE_DESTROY(TOPOLOGY%ELEMENTS%elementsTree,ERR,ERROR,*999)
        DEALLOCATE(TOPOLOGY%ELEMENTS)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the element data structures for a decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        CALL FlagError("Decomposition already has topology elements associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%ELEMENTS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology elements.",ERR,ERROR,*999)
        TOPOLOGY%ELEMENTS%numberOfElements=0
        TOPOLOGY%ELEMENTS%totalNumberOfElements=0
        TOPOLOGY%ELEMENTS%numberOfGlobalElements=0
        TOPOLOGY%ELEMENTS%DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS)
        NULLIFY(TOPOLOGY%ELEMENTS%elementsTree)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the topology in the given decomposition. \todo Pass in a pointer to the decomposition topology
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FINALISE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: DECOMPOSITION !<A pointer to the decomposition to finalise the topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      CALL DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      IF(DECOMPOSITION%calculateLines) THEN
        CALL DECOMPOSITION_TOPOLOGY_LINES_FINALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
      IF(DECOMPOSITION%calculateFaces) THEN
        CALL DECOMPOSITION_TOPOLOGY_FACES_FINALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
      DEALLOCATE(DECOMPOSITION%TOPOLOGY)
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DECOMPOSITION_TOPOLOGY_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the topology for a given decomposition.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_INITIALISE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: DECOMPOSITION !<A pointer to the decomposition to initialise the topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: meshComponentNumber

    ENTERS("DECOMPOSITION_TOPOLOGY_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%TOPOLOGY)) THEN
        CALL FlagError("Decomposition already has topology associated.",ERR,ERROR,*999)
      ELSE
        !Allocate decomposition topology
        ALLOCATE(DECOMPOSITION%TOPOLOGY,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Decomposition topology could not be allocated.",ERR,ERROR,*999)
        DECOMPOSITION%TOPOLOGY%DECOMPOSITION=>DECOMPOSITION
        NULLIFY(DECOMPOSITION%TOPOLOGY%ELEMENTS)
        NULLIFY(DECOMPOSITION%TOPOLOGY%LINES)
        NULLIFY(DECOMPOSITION%TOPOLOGY%FACES)
        NULLIFY(DECOMPOSITION%TOPOLOGY%dataPoints)
        !Initialise the topology components
        CALL DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        IF(DECOMPOSITION%calculateLines) THEN !Default is currently true
          CALL DECOMPOSITION_TOPOLOGY_LINES_INITIALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        ENDIF
        IF(DECOMPOSITION%calculateFaces) THEN !Default is currently false
          CALL DECOMPOSITION_TOPOLOGY_FACES_INITIALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        ENDIF      
        meshComponentNumber=decomposition%meshComponentNumber
        IF(ALLOCATED(DECOMPOSITION%MESH%TOPOLOGY(meshComponentNumber)%PTR%dataPoints%dataPoints)) THEN
          CALL DecompositionTopology_DataPointsInitialise(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        ENDIF   
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises a line in the given decomposition topology and deallocates all memory.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINE_FINALISE(LINE,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionLineType) :: LINE !<The decomposition line to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_LINE_FINALISE",ERR,ERROR,*999)

    LINE%NUMBER=0
    line%xiDirection=0
    LINE%numberOfSurroundingElements=0
    IF(ALLOCATED(LINE%surroundingElements)) DEALLOCATE(LINE%surroundingElements)
    IF(ALLOCATED(LINE%elementLines)) DEALLOCATE(LINE%elementLines)
    LINE%adjacentLines=0
 
    EXITS("DECOMPOSITION_TOPOLOGY_LINE_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_LINE_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the line data structure for a decomposition topology line.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINE_INITIALISE(LINE,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionLineType) :: LINE !<The decomposition line to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_LINE_INITIALISE",ERR,ERROR,*999)

    LINE%NUMBER=0
    line%xiDirection=0
    LINE%numberOfSurroundingElements=0
    LINE%adjacentLines=0
    LINE%boundaryLine=.FALSE.

    EXITS("DECOMPOSITION_TOPOLOGY_LINE_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_LINE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINE_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Calculates the lines in the given decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,element_idx,surrounding_element_idx,basis_local_line_idx, &
      & surrounding_element_basis_local_line_idx,element_local_node_idx,basis_local_line_node_idx,derivative_idx,version_idx, &
      & local_line_idx,surrounding_element_local_line_idx,node_idx,local_node_idx,elem_idx,line_end_node_idx,basis_node_idx, &
      & nodesInLine(4),numberOfLines,MAX_NUMBER_OF_LINES,NEW_MAX_NUMBER_OF_LINES,LINE_NUMBER,COUNT
    INTEGER(INTG), ALLOCATABLE :: NODES_NUMBER_OF_LINES(:)
    INTEGER(INTG), POINTER :: TEMP_LINES(:,:),NEW_TEMP_LINES(:,:)
    REAL(DP) :: APPROX_DIMENSION
    LOGICAL :: FOUND
    TYPE(BasisType), POINTER :: BASIS,BASIS2
    TYPE(DecompositionType), POINTER :: DECOMPOSITION
    TYPE(DecompositionElementType), POINTER :: DECOMPOSITION_ELEMENT
    TYPE(DecompositionElementsType), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(DecompositionLineType), POINTER :: DECOMPOSITION_LINE,DECOMPOSITION_LINE2
    TYPE(DecompositionLinesType), POINTER :: DECOMPOSITION_LINES
    TYPE(DomainType), POINTER :: DOMAIN
    TYPE(DomainElementType), POINTER :: DOMAIN_ELEMENT
    TYPE(DomainElementsType), POINTER :: DOMAIN_ELEMENTS
    TYPE(DomainLineType), POINTER :: DOMAIN_LINE,DOMAIN_LINE2
    TYPE(DomainLinesType), POINTER :: DOMAIN_LINES
    TYPE(DomainNodeType), POINTER :: DOMAIN_NODE
    TYPE(DomainNodesType), POINTER :: DOMAIN_NODES
    TYPE(DomainTopologyType), POINTER :: DOMAIN_TOPOLOGY
    TYPE(MeshType), POINTER :: MESH

    NULLIFY(TEMP_LINES)
    NULLIFY(NEW_TEMP_LINES)
    
    ENTERS("DECOMPOSITION_TOPOLOGY_LINES_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      DECOMPOSITION_LINES=>TOPOLOGY%LINES
      IF(ASSOCIATED(DECOMPOSITION_LINES)) THEN
        DECOMPOSITION_ELEMENTS=>TOPOLOGY%ELEMENTS
        IF(ASSOCIATED(DECOMPOSITION_ELEMENTS)) THEN
          DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
          IF(ASSOCIATED(DECOMPOSITION)) THEN
            !Process the mesh component number (component number the decomposition was calculated from) first to establish line
            !topology then process the other mesh components.
            DOMAIN=>DECOMPOSITION%DOMAIN(decomposition%meshComponentNumber)%PTR
            IF(ASSOCIATED(DOMAIN)) THEN
              DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
              IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
                IF(ASSOCIATED(DOMAIN_NODES)) THEN
                  DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                  IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                    !Guestimate the number of lines
                    SELECT CASE(DOMAIN%numberOfDimensions)
                    CASE(1)
                      MAX_NUMBER_OF_LINES=DOMAIN_ELEMENTS%totalNumberOfElements
                    CASE(2)
                      APPROX_DIMENSION=SQRT(REAL(DOMAIN_ELEMENTS%totalNumberOfElements,DP))
                      !This should give the maximum and will over estimate the number of lines for a "square mesh" by approx 33%
                      MAX_NUMBER_OF_LINES=NINT(3.0_DP*APPROX_DIMENSION*(APPROX_DIMENSION+1),INTG)
                    CASE(3)
                      !This should give the maximum and will over estimate the number of lines for a "cube mesh" by approx 73%
                      APPROX_DIMENSION=REAL(DOMAIN_ELEMENTS%totalNumberOfElements,DP)**(1.0_DP/3.0_DP)
                      MAX_NUMBER_OF_LINES=NINT(11.0_DP*APPROX_DIMENSION*APPROX_DIMENSION*(APPROX_DIMENSION+1),INTG)
                    CASE DEFAULT
                      CALL FlagError("Invalid number of dimensions for a topology domain.",ERR,ERROR,*999)
                    END SELECT
                    DOMAIN_LINES=>DOMAIN_TOPOLOGY%LINES
                    IF(ASSOCIATED(DOMAIN_LINES)) THEN
                      ALLOCATE(TEMP_LINES(4,MAX_NUMBER_OF_LINES),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate temporary lines array.",ERR,ERROR,*999)
                      ALLOCATE(NODES_NUMBER_OF_LINES(DOMAIN_NODES%totalNumberOfNodes),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate nodes number of lines array.",ERR,ERROR,*999)
                      NODES_NUMBER_OF_LINES=0
                      numberOfLines=0
                      TEMP_LINES=0
                      !Loop over the elements in the topology
                      DO element_idx=1,DOMAIN_ELEMENTS%totalNumberOfElements
                        DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                        DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(element_idx)
                        BASIS=>DOMAIN_ELEMENT%BASIS
                        ALLOCATE(DECOMPOSITION_ELEMENT%elementLines(BASIS%numberOfLocalLines),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate element element lines.",ERR,ERROR,*999)
                        !Loop over the local lines of the element
                        DO basis_local_line_idx=1,BASIS%numberOfLocalLines
                          !Calculate the topology node numbers that make up the line
                          nodesInLine=0
                          DO basis_local_line_node_idx=1,BASIS%numberOfNodesInLocalLine(basis_local_line_idx)
                            nodesInLine(basis_local_line_node_idx)=DOMAIN_ELEMENT%elementNodes( &
                              & BASIS%nodeNumbersInLocalLine(basis_local_line_node_idx,basis_local_line_idx))
                          ENDDO !basis_local_line_node_idx
                          !Try and find a previously created line that matches in the adjacent elements
                          FOUND=.FALSE.
                          node_idx=nodesInLine(1)
                          DO elem_idx=1,DOMAIN_NODES%NODES(node_idx)%numberOfSurroundingElements
                            surrounding_element_idx=DOMAIN_NODES%NODES(node_idx)%surroundingElements(elem_idx)
                            IF(surrounding_element_idx/=element_idx) THEN
                              IF(ALLOCATED(DECOMPOSITION_ELEMENTS%ELEMENTS(surrounding_element_idx)%elementLines)) THEN
                                BASIS2=>DOMAIN_ELEMENTS%ELEMENTS(surrounding_element_idx)%BASIS
                                DO surrounding_element_basis_local_line_idx=1,BASIS2%numberOfLocalLines
                                  local_line_idx=DECOMPOSITION_ELEMENTS%ELEMENTS(surrounding_element_idx)% &
                                    & elementLines(surrounding_element_basis_local_line_idx)
                                  IF(ALL(nodesInLine(1:BASIS%numberOfNodesInLocalLine(basis_local_line_idx))== &
                                    & TEMP_LINES(1:BASIS%numberOfNodesInLocalLine(basis_local_line_idx),local_line_idx))) THEN
                                    FOUND=.TRUE.
                                    EXIT
                                  ENDIF
                                ENDDO !surrounding_element_basis_local_line_idx
                                IF(FOUND) EXIT
                              ENDIF
                            ENDIF
                          ENDDO !elem_idx
                          IF(FOUND) THEN
                            !Line has already been created
                            DECOMPOSITION_ELEMENT%elementLines(basis_local_line_idx)=local_line_idx
                          ELSE
                            !Line has not been created
                            IF(numberOfLines==MAX_NUMBER_OF_LINES) THEN
                              !We are at maximum. Reallocate the LINES array to be 20% bigger and try again.
                              NEW_MAX_NUMBER_OF_LINES=NINT(1.20_DP*REAL(MAX_NUMBER_OF_LINES,DP),INTG)
                              ALLOCATE(NEW_TEMP_LINES(4,NEW_MAX_NUMBER_OF_LINES),STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate new number of lines.",ERR,ERROR,*999)
                              NEW_TEMP_LINES(:,1:numberOfLines)=TEMP_LINES(:,1:numberOfLines)
                              NEW_TEMP_LINES(:,numberOfLines+1:NEW_MAX_NUMBER_OF_LINES)=0
                              DEALLOCATE(TEMP_LINES)
                              TEMP_LINES=>NEW_TEMP_LINES
                              NULLIFY(NEW_TEMP_LINES)
                              MAX_NUMBER_OF_LINES=NEW_MAX_NUMBER_OF_LINES
                            ENDIF
                            numberOfLines=numberOfLines+1
                            TEMP_LINES(:,numberOfLines)=nodesInLine
                            DECOMPOSITION_ELEMENT%elementLines(basis_local_line_idx)=numberOfLines
                            DO basis_local_line_node_idx=1,SIZE(nodesInLine,1)
                              IF(nodesInLine(basis_local_line_node_idx)/=0) &
                                & NODES_NUMBER_OF_LINES(nodesInLine(basis_local_line_node_idx))= &
                                & NODES_NUMBER_OF_LINES(nodesInLine(basis_local_line_node_idx))+1
                            ENDDO !basis_local_line_node_idx
                          ENDIF
                        ENDDO !basis_local_line_idx
                      ENDDO !element_idx
                      !Allocate the line arrays and set them from the LINES and nodeLines arrays
                      DO node_idx=1,DOMAIN_NODES%totalNumberOfNodes
                        ALLOCATE(DOMAIN_NODES%NODES(node_idx)%nodeLines(NODES_NUMBER_OF_LINES(node_idx)),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate node lines array.",ERR,ERROR,*999)
                        DOMAIN_NODES%NODES(node_idx)%numberOfNodeLines=0
                      ENDDO !node_idx
                      DEALLOCATE(NODES_NUMBER_OF_LINES)
                      ALLOCATE(DECOMPOSITION_LINES%LINES(numberOfLines),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate decomposition topology lines.",ERR,ERROR,*999)
                      DECOMPOSITION_LINES%numberOfLines=numberOfLines
                      ALLOCATE(DOMAIN_LINES%LINES(numberOfLines),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate domain topology lines.",ERR,ERROR,*999)
                      DOMAIN_LINES%numberOfLines=numberOfLines
                      DO local_line_idx=1,DOMAIN_LINES%numberOfLines
                        CALL DECOMPOSITION_TOPOLOGY_LINE_INITIALISE(DECOMPOSITION_LINES%LINES(local_line_idx),ERR,ERROR,*999)
                        CALL DOMAIN_TOPOLOGY_LINE_INITIALISE(DOMAIN_LINES%LINES(local_line_idx),ERR,ERROR,*999)
                        DO basis_local_line_node_idx=1,SIZE(TEMP_LINES,1)
                          IF(TEMP_LINES(basis_local_line_node_idx,local_line_idx)/=0) THEN
                            node_idx=TEMP_LINES(basis_local_line_node_idx,local_line_idx)
                            DOMAIN_NODES%NODES(node_idx)%numberOfNodeLines=DOMAIN_NODES%NODES(node_idx)%numberOfNodeLines+1
                            DOMAIN_NODES%NODES(node_idx)%nodeLines(DOMAIN_NODES%NODES(node_idx)%numberOfNodeLines)= &
                              & local_line_idx
                          ENDIF
                        ENDDO !basis_local_line_node_idx  
                      ENDDO !local_line_idx
                      DO element_idx=1,DECOMPOSITION_ELEMENTS%totalNumberOfElements
                        DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(element_idx)
                        DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                        BASIS=>DOMAIN_ELEMENT%BASIS
                        DO basis_local_line_idx=1,BASIS%numberOfLocalLines
                          LINE_NUMBER=DECOMPOSITION_ELEMENT%elementLines(basis_local_line_idx)
                          DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(LINE_NUMBER)
                          DOMAIN_LINE=>DOMAIN_LINES%LINES(LINE_NUMBER)
                          DECOMPOSITION_LINE%numberOfSurroundingElements=DECOMPOSITION_LINE%numberOfSurroundingElements+1
                          IF(.NOT.ASSOCIATED(DOMAIN_LINE%BASIS)) THEN
                            DECOMPOSITION_LINE%NUMBER=LINE_NUMBER
                            DOMAIN_LINE%NUMBER=LINE_NUMBER
                            DOMAIN_LINE%elementNumber=element_idx !Needs checking
                            DECOMPOSITION_line%xiDirection=BASIS%localLineXiDirection(basis_local_line_idx)
                            IF(ALLOCATED(BASIS%localLineBasis)) THEN
                              DOMAIN_LINE%BASIS=>BASIS%localLineBasis(basis_local_line_idx)%PTR
                            ELSE
                              !Basis is only 1D
                              DOMAIN_LINE%BASIS=>BASIS
                            ENDIF
                            ALLOCATE(DOMAIN_LINE%nodesInLine(BASIS%numberOfNodesInLocalLine(basis_local_line_idx)),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate line nodes in line.",ERR,ERROR,*999)
                            ALLOCATE(DOMAIN_LINE%derivativesInLine(2,DOMAIN_LINE%BASIS%maximumNumberOfDerivatives, &
                              & BASIS%numberOfNodesInLocalLine(basis_local_line_idx)),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate line derivatives in line.",ERR,ERROR,*999)
                            DOMAIN_LINE%nodesInLine(1:BASIS%numberOfNodesInLocalLine(basis_local_line_idx))= &
                              & TEMP_LINES(1:BASIS%numberOfNodesInLocalLine(basis_local_line_idx),LINE_NUMBER)
                            DO basis_local_line_node_idx=1,BASIS%numberOfNodesInLocalLine(basis_local_line_idx)
                              !Set derivative number of u (NO_GLOBAL_DERIV) for the domain line
                              DOMAIN_LINE%derivativesInLine(1,1,basis_local_line_node_idx)=NO_GLOBAL_DERIV
                              !Set version number of u (NO_GLOBAL_DERIV) for the domain line
                              version_idx=DOMAIN_ELEMENT%elementVersions(1,BASIS%nodeNumbersInLocalLine( &
                                & basis_local_line_node_idx,basis_local_line_idx))
                              DOMAIN_LINE%derivativesInLine(2,1,basis_local_line_node_idx)=version_idx
                              IF(DOMAIN_LINE%BASIS%maximumNumberOfDerivatives>1) THEN
                                derivative_idx=DOMAIN_ELEMENT%elementDerivatives( &
                                  & BASIS%derivativeNumbersInLocalLine(basis_local_line_node_idx,basis_local_line_idx), &
                                  & BASIS%nodeNumbersInLocalLine(basis_local_line_node_idx,basis_local_line_idx))
                                DOMAIN_LINE%derivativesInLine(1,2,basis_local_line_node_idx)=derivative_idx
                                version_idx=DOMAIN_ELEMENT%elementVersions(BASIS%derivativeNumbersInLocalLine( &
                                  & basis_local_line_node_idx,basis_local_line_idx),BASIS%nodeNumbersInLocalLine( &
                                  & basis_local_line_node_idx,basis_local_line_idx))
                                DOMAIN_LINE%derivativesInLine(2,2,basis_local_line_node_idx)=version_idx
                              ENDIF
                            ENDDO !basis_local_line_node_idx
                          ENDIF
                        ENDDO !basis_local_line_idx
                      ENDDO !element_idx
                      DEALLOCATE(TEMP_LINES)
                      !Calculate adjacent lines and the surrounding elements for each line
                      DO local_line_idx=1,DECOMPOSITION_LINES%numberOfLines
                        DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(local_line_idx)
                        DOMAIN_LINE=>DOMAIN_LINES%LINES(local_line_idx)
                        BASIS=>DOMAIN_LINE%BASIS
                        IF(DECOMPOSITION_LINE%numberOfSurroundingElements==1) THEN
                          DECOMPOSITION_LINE%boundaryLine=.TRUE.
                          DOMAIN_LINE%boundaryLine=.TRUE.
                        ENDIF
                        !Allocate the elements surrounding the line
                        ALLOCATE(DECOMPOSITION_LINE%surroundingElements(DECOMPOSITION_LINE%numberOfSurroundingElements), &
                          & STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate line surrounding elements.",ERR,ERROR,*999)
                        ALLOCATE(DECOMPOSITION_LINE%elementLines(DECOMPOSITION_LINE%numberOfSurroundingElements), &
                          & STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate line element lines.",ERR,ERROR,*999)
                        DECOMPOSITION_LINE%numberOfSurroundingElements=0
                        DECOMPOSITION_LINE%adjacentLines=0
                        !Loop over the nodes at each end of the line
                        DO line_end_node_idx=0,1
                          FOUND=.FALSE.
                          node_idx=DOMAIN_LINE%nodesInLine(line_end_node_idx*(BASIS%numberOfNodes-1)+1)
                          !Loop over the elements surrounding the node.
                          DO elem_idx=1,DOMAIN_NODES%NODES(node_idx)%numberOfSurroundingElements
                            element_idx=DOMAIN_NODES%NODES(node_idx)%surroundingElements(elem_idx)
                            DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(element_idx)
                            DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                            !Loop over the local lines of the element
                            DO basis_local_line_idx=1,DOMAIN_ELEMENT%BASIS%numberOfLocalLines
                              surrounding_element_local_line_idx=DECOMPOSITION_ELEMENT%elementLines(basis_local_line_idx)
                              IF(surrounding_element_local_line_idx/=local_line_idx) THEN
                                DECOMPOSITION_LINE2=>DECOMPOSITION_LINES%LINES(surrounding_element_local_line_idx)
                                DOMAIN_LINE2=>DOMAIN_LINES%LINES(surrounding_element_local_line_idx)
                                IF(DECOMPOSITION_LINE2%xiDirection==DECOMPOSITION_line%xiDirection) THEN
                                  !Lines run in the same direction.
                                  BASIS2=>DOMAIN_LINE2%BASIS
                                  IF(line_end_node_idx==0) THEN
                                    local_node_idx=DOMAIN_LINE2%nodesInLine(BASIS2%numberOfNodes)
                                  ELSE
                                    local_node_idx=DOMAIN_LINE2%nodesInLine(1)
                                  ENDIF
                                  IF(local_node_idx==node_idx) THEN
                                    !The node at the 'other' end of this line matches the node at the current end of the line.
                                    !Check it is not a coexistant line running the other way
                                    IF(BASIS2%interpolationOrder(1)==BASIS%interpolationOrder(1)) THEN
                                      COUNT=0
                                      DO basis_node_idx=1,BASIS%numberOfNodes
                                        IF(DOMAIN_LINE2%nodesInLine(basis_node_idx)== &
                                          & DOMAIN_LINE%nodesInLine(BASIS2%numberOfNodes-basis_node_idx+1)) &
                                          & COUNT=COUNT+1
                                      ENDDO !basis_node_idx
                                      IF(COUNT<BASIS%numberOfNodes) THEN
                                        FOUND=.TRUE.
                                        EXIT
                                      ENDIF
                                    ELSE
                                      FOUND=.TRUE.
                                      EXIT
                                    ENDIF
                                  ENDIF
                                ENDIF
                              ENDIF
                            ENDDO !basis_local_line_idx
                            IF(FOUND) EXIT
                          ENDDO !element_idx
                          IF(FOUND) DECOMPOSITION_LINE%adjacentLines(line_end_node_idx)=surrounding_element_local_line_idx
                        ENDDO !line_end_node_idx
                      ENDDO !local_line_idx
                      !Set the surrounding elements
                      DO element_idx=1,DECOMPOSITION_ELEMENTS%totalNumberOfElements
                        DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(element_idx)
                        DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                        BASIS=>DOMAIN_ELEMENT%BASIS
                        DO basis_local_line_idx=1,BASIS%numberOfLocalLines
                          LINE_NUMBER=DECOMPOSITION_ELEMENT%elementLines(basis_local_line_idx)
                          DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(LINE_NUMBER)
                          DECOMPOSITION_LINE%numberOfSurroundingElements=DECOMPOSITION_LINE%numberOfSurroundingElements+1
                          DECOMPOSITION_LINE%surroundingElements(DECOMPOSITION_LINE%numberOfSurroundingElements)=element_idx
                          DECOMPOSITION_LINE%elementLines(DECOMPOSITION_LINE%numberOfSurroundingElements)=basis_local_line_idx
                        ENDDO !basis_local_line_idx
                      ENDDO !element_idx
                    ELSE
                      CALL FlagError("Domain topology lines is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Domain topology elements is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Domain topology nodes is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Topology decomposition domain topology is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Topology decomposition domain is not associated.",ERR,ERROR,*999)
            ENDIF
            !Now loop over the other mesh components in the decomposition and calculate the domain lines
            MESH=>DECOMPOSITION%MESH
            IF(ASSOCIATED(MESH)) THEN
              DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
                IF(component_idx/=decomposition%meshComponentNumber) THEN
                  DOMAIN=>DECOMPOSITION%DOMAIN(component_idx)%PTR
                  IF(ASSOCIATED(DOMAIN)) THEN
                    DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
                    IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                      DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
                      IF(ASSOCIATED(DOMAIN_NODES)) THEN
                        DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                        IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                          DOMAIN_LINES=>DOMAIN_TOPOLOGY%LINES                      
                          IF(ASSOCIATED(DOMAIN_LINES)) THEN
                            ALLOCATE(DOMAIN_LINES%LINES(DECOMPOSITION_LINES%numberOfLines),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate domain lines lines.",ERR,ERROR,*999)
                            DOMAIN_LINES%numberOfLines=DECOMPOSITION_LINES%numberOfLines
                            ALLOCATE(NODES_NUMBER_OF_LINES(DOMAIN_NODES%totalNumberOfNodes),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate nodes number of lines array.",ERR,ERROR,*999)
                            NODES_NUMBER_OF_LINES=0
                            !Loop over the lines in the topology
                            DO local_line_idx=1,DECOMPOSITION_LINES%numberOfLines
                              DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(local_line_idx)
                              DOMAIN_LINE=>DOMAIN_LINES%LINES(local_line_idx)
                              IF(DECOMPOSITION_LINE%numberOfSurroundingElements>0) THEN
                                element_idx=DECOMPOSITION_LINE%surroundingElements(1)
                                basis_local_line_idx=DECOMPOSITION_LINE%elementLines(1)
                                CALL DOMAIN_TOPOLOGY_LINE_INITIALISE(DOMAIN_LINES%LINES(local_line_idx),ERR,ERROR,*999)
                                DOMAIN_LINE%NUMBER=local_line_idx
                                DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                                BASIS=>DOMAIN_ELEMENT%BASIS
                                DOMAIN_LINE%elementNumber=DOMAIN_ELEMENT%NUMBER
                                IF(ALLOCATED(BASIS%localLineBasis)) THEN
                                  DOMAIN_LINE%BASIS=>BASIS%localLineBasis(basis_local_line_idx)%PTR
                                ELSE
                                  !Basis is only 1D
                                  DOMAIN_LINE%BASIS=>BASIS
                                ENDIF
                                ALLOCATE(DOMAIN_LINE%nodesInLine(BASIS%numberOfNodesInLocalLine(basis_local_line_idx)), &
                                  & STAT=ERR)
                                IF(ERR/=0) CALL FlagError("Could not allocate nodes in line.",ERR,ERROR,*999)
                                ALLOCATE(DOMAIN_LINE%derivativesInLine(2,DOMAIN_LINE%BASIS%maximumNumberOfDerivatives, &
                                  & BASIS%numberOfNodesInLocalLine(basis_local_line_idx)),STAT=ERR)
                                IF(ERR/=0) CALL FlagError("Could not allocate derivatives in line.",ERR,ERROR,*999)
                                DO basis_local_line_node_idx=1,BASIS%numberOfNodesInLocalLine(basis_local_line_idx)
                                  element_local_node_idx=BASIS%nodeNumbersInLocalLine(basis_local_line_node_idx, &
                                    & basis_local_line_idx)
                                  node_idx=DOMAIN_ELEMENT%elementNodes(element_local_node_idx)
                                  DOMAIN_LINE%nodesInLine(basis_local_line_node_idx)=node_idx
                                  !Set derivative number of u (NO_GLOBAL_DERIV) for the domain line
                                  DOMAIN_LINE%derivativesInLine(1,1,basis_local_line_node_idx)=NO_GLOBAL_DERIV
                                  !Set version number of u (NO_GLOBAL_DERIV) for the domain line
                                  version_idx=DOMAIN_ELEMENT%elementVersions(1,BASIS%nodeNumbersInLocalLine( &
                                    & basis_local_line_node_idx,basis_local_line_idx))
                                  DOMAIN_LINE%derivativesInLine(2,1,basis_local_line_node_idx)=version_idx
                                  IF(DOMAIN_LINE%BASIS%maximumNumberOfDerivatives>1) THEN
                                    derivative_idx=DOMAIN_ELEMENT%elementDerivatives(BASIS%derivativeNumbersInLocalLine( &
                                      & basis_local_line_node_idx,basis_local_line_idx),element_local_node_idx)
                                    DOMAIN_LINE%derivativesInLine(1,2,basis_local_line_node_idx)=derivative_idx
                                    version_idx=DOMAIN_ELEMENT%elementVersions(BASIS%derivativeNumbersInLocalLine( &
                                      & basis_local_line_node_idx,basis_local_line_idx),BASIS%nodeNumbersInLocalLine( &
                                      & basis_local_line_node_idx,basis_local_line_idx))
                                    DOMAIN_LINE%derivativesInLine(2,2,basis_local_line_node_idx)=version_idx
                                  ENDIF
                                  NODES_NUMBER_OF_LINES(node_idx)=NODES_NUMBER_OF_LINES(node_idx)+1
                                ENDDO !basis_local_line_node_idx
                              ELSE
                                CALL FlagError("Line is not surrounded by any elements?",ERR,ERROR,*999)
                              ENDIF
                            ENDDO !local_line_idx
                            DO node_idx=1,DOMAIN_NODES%totalNumberOfNodes
                              ALLOCATE(DOMAIN_NODES%NODES(node_idx)%nodeLines(NODES_NUMBER_OF_LINES(node_idx)),STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate node lines.",ERR,ERROR,*999)
                              DOMAIN_NODES%NODES(node_idx)%numberOfNodeLines=0
                            ENDDO !node_idx
                            DEALLOCATE(NODES_NUMBER_OF_LINES)
                            DO local_line_idx=1,DOMAIN_LINES%numberOfLines
                              DOMAIN_LINE=>DOMAIN_LINES%LINES(local_line_idx)
                              BASIS=>DOMAIN_LINE%BASIS
                              DO basis_local_line_node_idx=1,BASIS%numberOfNodes
                                node_idx=DOMAIN_LINE%nodesInLine(basis_local_line_node_idx)
                                DOMAIN_NODE=>DOMAIN_NODES%NODES(node_idx)
                                DOMAIN_NODE%numberOfNodeLines=DOMAIN_NODE%numberOfNodeLines+1
                                DOMAIN_NODE%nodeLines(DOMAIN_NODE%numberOfNodeLines)=local_line_idx
                              ENDDO !basis_local_line_node_idx
                            ENDDO !local_line_idx
                          ELSE
                            CALL FlagError("Domain lines is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Domain elements is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Domain nodes is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Domain topology is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Decomposition mesh is not associated",ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ENDDO !component_idx
            ELSE
              CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
            ENDIF                        
          ELSE
            CALL FlagError("Topology decomposition is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Topology decomposition elements is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Topology lines is not associated.",ERR,ERROR,*999)

      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Decomposition topology lines:",ERR,ERROR,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of mesh components = ",MESH%NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of lines = ",DECOMPOSITION_LINES%numberOfLines,ERR,ERROR,*999)
      DO local_line_idx=1,DECOMPOSITION_LINES%numberOfLines
        DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(local_line_idx)
        DOMAIN_LINE=>DOMAIN_LINES%LINES(local_line_idx)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Line number = ",DECOMPOSITION_LINE%NUMBER,ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Xi direction = ",DECOMPOSITION_line%xiDirection,ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of surrounding elements = ", &
          & DECOMPOSITION_LINE%numberOfSurroundingElements,ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_LINE%numberOfSurroundingElements,4,4, &
          & DECOMPOSITION_LINE%surroundingElements,'("      Surrounding elements :",4(X,I8))','(28X,4(X,I8))',ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_LINE%numberOfSurroundingElements,4,4, &
          & DECOMPOSITION_LINE%elementLines,'("      Element lines        :",4(X,I8))','(28X,4(X,I8))',ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,2,2,DECOMPOSITION_LINE%adjacentLines, &
          & '("      Adjacent lines       :",2(X,I8))','(28X,2(X,I8))',ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Boundary line = ",DECOMPOSITION_LINE%boundaryLine,ERR,ERROR,*999)
        DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Mesh component : ",component_idx,ERR,ERROR,*999)
          DOMAIN=>DECOMPOSITION%DOMAIN(component_idx)%PTR
          DOMAIN_LINE=>DOMAIN%TOPOLOGY%LINES%LINES(local_line_idx)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Basis user number = ",DOMAIN_LINE%BASIS%userNumber, &
            & ERR,ERROR,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Basis family number = ",DOMAIN_LINE%BASIS%familyNumber, &
            & ERR,ERROR,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation type = ",DOMAIN_LINE%BASIS% &
            & interpolationType(1),ERR,ERROR,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation order = ",DOMAIN_LINE%BASIS% &
            & interpolationOrder(1),ERR,ERROR,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of nodes in lines = ",DOMAIN_LINE%BASIS%numberOfNodes, &
            & ERR,ERROR,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_LINE%BASIS%numberOfNodes,4,4,DOMAIN_LINE%nodesInLine, &
            & '("        Nodes in line        :",4(X,I8))','(30X,4(X,I8))',ERR,ERROR,*999)
          DO basis_local_line_node_idx=1,DOMAIN_LINE%BASIS%numberOfNodes
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Node : ",basis_local_line_node_idx,ERR,ERROR,*999)
            !/TODO::Loop over local_derivative index so this output makes more sense !<derivativesInLine(i,local_derivative_idx,local_node_idx)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & DOMAIN_LINE%BASIS%numberOfDerivatives(basis_local_line_node_idx),4,4, &
              & DOMAIN_LINE%derivativesInLine(1,:,basis_local_line_node_idx),'("            Derivatives in line  :",4(X,I8))', &
              & '(34X,4(X,I8))',ERR,ERROR,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & DOMAIN_LINE%BASIS%numberOfDerivatives(basis_local_line_node_idx),4,4, &
              & DOMAIN_LINE%derivativesInLine(2,:,basis_local_line_node_idx), &
              & '("            Derivatives Versions in line  :",4(X,I8))','(34X,4(X,I8))',ERR,ERROR,*999)
          ENDDO !basis_local_line_node_idx
        ENDDO !component_idx
      ENDDO !local_line_idx
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_LINES_CALCULATE")
    RETURN
999 IF(ASSOCIATED(TEMP_LINES)) DEALLOCATE(TEMP_LINES)
    IF(ASSOCIATED(NEW_TEMP_LINES)) DEALLOCATE(NEW_TEMP_LINES)
    IF(ALLOCATED(NODES_NUMBER_OF_LINES)) DEALLOCATE(NODES_NUMBER_OF_LINES)
    ERRORSEXITS("DECOMPOSITION_TOPOLOGY_LINES_CALCULATE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalises the lines in the given decomposition topology. \todo Pass in the topology lines
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nl
    
    ENTERS("DECOMPOSITION_TOPOLOGY_LINES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%LINES)) THEN
        DO nl=1,TOPOLOGY%LINES%numberOfLines
          CALL DECOMPOSITION_TOPOLOGY_LINE_FINALISE(TOPOLOGY%LINES%LINES(nl),ERR,ERROR,*999)
        ENDDO !nl
        IF(ALLOCATED(TOPOLOGY%LINES%LINES)) DEALLOCATE(TOPOLOGY%LINES%LINES)
        DEALLOCATE(TOPOLOGY%LINES)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DECOMPOSITION_TOPOLOGY_LINES_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_LINES_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the line data structures for a decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to initialise the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_LINES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%LINES)) THEN
        CALL FlagError("Decomposition already has topology lines associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%LINES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology lines.",ERR,ERROR,*999)
        TOPOLOGY%LINES%numberOfLines=0
        TOPOLOGY%LINES%DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_LINES_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_LINES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the line data structures for a decomposition topology.
  SUBROUTINE DecompositionTopology_DataPointsInitialise(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to initialise the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DecompositionTopology_DataPointsInitialise",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%dataPoints)) THEN
        CALL FlagError("Decomposition already has topology data points associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%dataPoints,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology data points.",ERR,ERROR,*999)
        TOPOLOGY%dataPoints%numberOfDataPoints=0
        TOPOLOGY%dataPoints%totalNumberOfDataPoints=0
        TOPOLOGY%dataPoints%numberOfGlobalDataPoints=0
        NULLIFY(TOPOLOGY%dataPoints%dataPointsTree)
        TOPOLOGY%dataPoints%DECOMPOSITION=>TOPOLOGY%DECOMPOSITION      
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DecompositionTopology_DataPointsInitialise")
    RETURN
999 ERRORSEXITS("DecompositionTopology_DataPointsInitialise",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DecompositionTopology_DataPointsInitialise
  
  !
  !================================================================================================================================
  !
 
  !>Finalises a face in the given decomposition topology and deallocates all memory.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACE_FINALISE(FACE,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionFaceType) :: FACE !<The decomposition face to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_FACE_FINALISE",ERR,ERROR,*999)

    FACE%NUMBER=0
    face%xiNormalDirection=0
    FACE%numberOfSurroundingElements=0
    IF(ALLOCATED(FACE%surroundingElements)) DEALLOCATE(FACE%surroundingElements)
    IF(ALLOCATED(FACE%elementFaces)) DEALLOCATE(FACE%elementFaces)
!    FACE%ADJACENT_FACES=0
 
    EXITS("DECOMPOSITION_TOPOLOGY_FACE_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FACE_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the face data structure for a decomposition topology face.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACE_INITIALISE(FACE,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionFaceType) :: FACE !<The decomposition face to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_FACE_INITIALISE",ERR,ERROR,*999)

    FACE%NUMBER=0
    face%xiNormalDirection=0
    FACE%numberOfSurroundingElements=0
!    FACE%ADJACENT_FACES=0
    FACE%boundaryFace=.FALSE.
    
    EXITS("DECOMPOSITION_TOPOLOGY_FACE_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FACE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACE_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Calculates the faces in the given decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the faces for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,ne,surrounding_element_idx,basis_local_face_idx,surrounding_element_basis_local_face_idx, &
      & element_local_node_idx,basis_local_face_node_idx,basis_local_face_derivative_idx,derivative_idx,version_idx,face_idx, &
      & node_idx,elem_idx,nodesInFace(16),numberOfFaces,MAX_NUMBER_OF_FACES,NEW_MAX_NUMBER_OF_FACES,FACE_NUMBER
    INTEGER(INTG), ALLOCATABLE :: NODES_NUMBER_OF_FACES(:)
    INTEGER(INTG), POINTER :: TEMP_FACES(:,:),NEW_TEMP_FACES(:,:)
    LOGICAL :: FOUND
    TYPE(BasisType), POINTER :: BASIS,BASIS2
    TYPE(DecompositionType), POINTER :: DECOMPOSITION
    TYPE(DecompositionElementType), POINTER :: DECOMPOSITION_ELEMENT
    TYPE(DecompositionElementsType), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(DecompositionFaceType), POINTER :: DECOMPOSITION_FACE!,DECOMPOSITION_FACE2
    TYPE(DecompositionFacesType), POINTER :: DECOMPOSITION_FACES
    TYPE(DomainType), POINTER :: DOMAIN
    TYPE(DomainElementType), POINTER :: DOMAIN_ELEMENT
    TYPE(DomainElementsType), POINTER :: DOMAIN_ELEMENTS
    TYPE(DomainFaceType), POINTER :: DOMAIN_FACE!,DOMAIN_FACE2
    TYPE(DomainFacesType), POINTER :: DOMAIN_FACES
    TYPE(DomainNodeType), POINTER :: DOMAIN_NODE
    TYPE(DomainNodesType), POINTER :: DOMAIN_NODES
    TYPE(DomainTopologyType), POINTER :: DOMAIN_TOPOLOGY
    TYPE(MeshType), POINTER :: MESH

    NULLIFY(TEMP_FACES)
    NULLIFY(NEW_TEMP_FACES)
    
    ENTERS("DECOMPOSITION_TOPOLOGY_FACES_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      DECOMPOSITION_FACES=>TOPOLOGY%FACES
      IF(ASSOCIATED(DECOMPOSITION_FACES)) THEN
        DECOMPOSITION_ELEMENTS=>TOPOLOGY%ELEMENTS
        IF(ASSOCIATED(DECOMPOSITION_ELEMENTS)) THEN
          DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
          IF(ASSOCIATED(DECOMPOSITION)) THEN
            !Process the mesh component number (component number the decomposition was calculated from) first to establish face
            !topology then process the other mesh components.
            DOMAIN=>DECOMPOSITION%DOMAIN(decomposition%meshComponentNumber)%PTR
            IF(ASSOCIATED(DOMAIN)) THEN
              DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
              IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
                IF(ASSOCIATED(DOMAIN_NODES)) THEN
                  DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                  IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                    !Estimate the number of faces
                    SELECT CASE(DOMAIN%numberOfDimensions)
                    CASE(1)
                      ! Faces not calculated in 1D 
                    CASE(2)
                      ! Faces not calculated in 2D
                    CASE(3)
                      !This should give the maximum and will over estimate the number of faces for a "cube mesh" by approx 33%
                      MAX_NUMBER_OF_FACES= &
                        & NINT(((REAL(DOMAIN_ELEMENTS%totalNumberOfElements,DP)*5.0_DP)+1.0_DP)*(4.0_DP/3.0_DP),INTG)

                      DOMAIN_FACES=>DOMAIN_TOPOLOGY%FACES
                      IF(ASSOCIATED(DOMAIN_FACES)) THEN
                        ALLOCATE(TEMP_FACES(16,MAX_NUMBER_OF_FACES),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate temporary faces array",ERR,ERROR,*999)
                        ALLOCATE(NODES_NUMBER_OF_FACES(DOMAIN_NODES%totalNumberOfNodes),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate nodes number of faces array",ERR,ERROR,*999)
                        NODES_NUMBER_OF_FACES=0
                        numberOfFaces=0
                        TEMP_FACES=0
                        !Loop over the elements in the topology and fill temp_faces with node numbers for each element
                        DO ne=1,DOMAIN_ELEMENTS%totalNumberOfElements
                          DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                          DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
                          BASIS=>DOMAIN_ELEMENT%BASIS
                          ALLOCATE(DECOMPOSITION_ELEMENT%elementFaces(BASIS%numberOfLocalFaces),STAT=ERR)
                          IF(ERR/=0) CALL FlagError("Could not allocate element faces of element",ERR,ERROR,*999)
                          !Loop over the local faces of the element
                          DO basis_local_face_idx=1,BASIS%numberOfLocalFaces
                            !Calculate the topology node numbers that make up the face
                            nodesInFace=0
                            !Check whether face has already been read out
                            DO basis_local_face_node_idx=1,BASIS%numberOfNodesInLocalFace(basis_local_face_idx)
                              !Read out node numbers of local face from ELEMENT_NODES
                              nodesInFace(basis_local_face_node_idx)=DOMAIN_ELEMENT%elementNodes( &
                                & BASIS%nodeNumbersInLocalFace(basis_local_face_node_idx,basis_local_face_idx))
                            ENDDO !basis_local_face_node_idx
                            !Try and find a previously created face that matches in the adjacent elements
                            FOUND=.FALSE.
                            node_idx=nodesInFace(1)
                            DO elem_idx=1,DOMAIN_NODES%NODES(node_idx)%numberOfSurroundingElements
                              surrounding_element_idx=DOMAIN_NODES%NODES(node_idx)%surroundingElements(elem_idx)
                              IF(surrounding_element_idx/=ne) THEN
                                IF(ALLOCATED(DECOMPOSITION_ELEMENTS%ELEMENTS(surrounding_element_idx)%elementFaces)) THEN
                                  BASIS2=>DOMAIN_ELEMENTS%ELEMENTS(surrounding_element_idx)%BASIS
                                  DO surrounding_element_basis_local_face_idx=1,BASIS2%numberOfLocalFaces
                                    face_idx=DECOMPOSITION_ELEMENTS%ELEMENTS(surrounding_element_idx)%elementFaces( &
                                      & surrounding_element_basis_local_face_idx)
                                    IF(ALL(nodesInFace(1:BASIS%numberOfNodesInLocalFace(basis_local_face_idx))== &
                                      & TEMP_FACES(1:BASIS%numberOfNodesInLocalFace(basis_local_face_idx),face_idx))) THEN
                                      FOUND=.TRUE.
                                      EXIT
                                    ENDIF
                                  ENDDO !surrounding_element_basis_local_face_idx
                                  IF(FOUND) EXIT
                                ENDIF
                              ENDIF
                            ENDDO !elem_idx
                            IF(FOUND) THEN
                              !Face has already been created
                              DECOMPOSITION_ELEMENT%elementFaces(basis_local_face_idx)=face_idx
                            ELSE
                              !Face has not been created
                              IF(numberOfFaces==MAX_NUMBER_OF_FACES) THEN
                                !We are at maximum. Reallocate the FACES array to be 20% bigger and try again.
                                NEW_MAX_NUMBER_OF_FACES=NINT(1.20_DP*REAL(MAX_NUMBER_OF_FACES,DP),INTG)
                                !\todo: Change 16 to a variable and above for nodesInFace
                                ALLOCATE(NEW_TEMP_FACES(16,NEW_MAX_NUMBER_OF_FACES),STAT=ERR)
                                IF(ERR/=0) CALL FlagError("Could not allocate new number of faces",ERR,ERROR,*999)
                                NEW_TEMP_FACES(:,1:numberOfFaces)=TEMP_FACES(:,1:numberOfFaces)
                                NEW_TEMP_FACES(:,numberOfFaces+1:NEW_MAX_NUMBER_OF_FACES)=0
                                DEALLOCATE(TEMP_FACES)
                                TEMP_FACES=>NEW_TEMP_FACES
                                NULLIFY(NEW_TEMP_FACES)
                                MAX_NUMBER_OF_FACES=NEW_MAX_NUMBER_OF_FACES
                              ENDIF
                              numberOfFaces=numberOfFaces+1
                              TEMP_FACES(:,numberOfFaces)=nodesInFace(:)
                              DECOMPOSITION_ELEMENT%elementFaces(basis_local_face_idx)=numberOfFaces
                              DO basis_local_face_node_idx=1,SIZE(nodesInFace,1)
                                IF(nodesInFace(basis_local_face_node_idx)/=0) &
                                  & NODES_NUMBER_OF_FACES(nodesInFace(basis_local_face_node_idx))= &
                                  & NODES_NUMBER_OF_FACES(nodesInFace(basis_local_face_node_idx))+1
                              ENDDO !basis_local_face_node_idx
                            ENDIF
                          ENDDO !basis_local_face_idx
                        ENDDO !ne

                        !Allocate the face arrays and set them from the FACES and nodeFaces arrays
                        DO node_idx=1,DOMAIN_NODES%totalNumberOfNodes
                          ALLOCATE(DOMAIN_NODES%NODES(node_idx)%nodeFaces(NODES_NUMBER_OF_FACES(node_idx)),STAT=ERR)
                          IF(ERR/=0) CALL FlagError("Could not allocate node faces array",ERR,ERROR,*999)
                          DOMAIN_NODES%NODES(node_idx)%numberOfNodeFaces=0
                        ENDDO !node_idx
                        DEALLOCATE(NODES_NUMBER_OF_FACES)
                        ALLOCATE(DECOMPOSITION_FACES%FACES(numberOfFaces),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate decomposition topology faces",ERR,ERROR,*999)
                        DECOMPOSITION_FACES%numberOfFaces=numberOfFaces
                        ALLOCATE(DOMAIN_FACES%FACES(numberOfFaces),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate domain topology faces",ERR,ERROR,*999)
                        DOMAIN_FACES%numberOfFaces=numberOfFaces
                        DO face_idx=1,DOMAIN_FACES%numberOfFaces
                          CALL DECOMPOSITION_TOPOLOGY_FACE_INITIALISE(DECOMPOSITION_FACES%FACES(face_idx),ERR,ERROR,*999)
                          CALL DOMAIN_TOPOLOGY_FACE_INITIALISE(DOMAIN_FACES%FACES(face_idx),ERR,ERROR,*999)
                          DO basis_local_face_node_idx=1,SIZE(TEMP_FACES,1)
                            IF(TEMP_FACES(basis_local_face_node_idx,face_idx)/=0) THEN
                              node_idx=TEMP_FACES(basis_local_face_node_idx,face_idx)
                              DOMAIN_NODES%NODES(node_idx)%numberOfNodeFaces=DOMAIN_NODES%NODES(node_idx)%numberOfNodeFaces+1
                              DOMAIN_NODES%NODES(node_idx)%nodeFaces(DOMAIN_NODES%NODES(node_idx)%numberOfNodeFaces)=face_idx
                            ENDIF
                          ENDDO !basis_local_face_node_idx  
                        ENDDO !face_idx

                        !Set nodes in face and derivatives of nodes in face for domain faces
                        DO ne=1,DECOMPOSITION_ELEMENTS%totalNumberOfElements
                          DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
                          DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                          BASIS=>DOMAIN_ELEMENT%BASIS
                          !Loop over local faces of element
                          DO basis_local_face_idx=1,BASIS%numberOfLocalFaces
                            FACE_NUMBER=DECOMPOSITION_ELEMENT%elementFaces(basis_local_face_idx)
                            DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(FACE_NUMBER)
                            DOMAIN_FACE=>DOMAIN_FACES%FACES(FACE_NUMBER)
                            DECOMPOSITION_FACE%numberOfSurroundingElements=DECOMPOSITION_FACE%numberOfSurroundingElements+1
                            IF(.NOT.ASSOCIATED(DOMAIN_FACE%BASIS)) THEN
                              DECOMPOSITION_FACE%NUMBER=FACE_NUMBER
                              DOMAIN_FACE%NUMBER=FACE_NUMBER
                              DOMAIN_FACE%elementNumber=ne !! Needs checking
!                              DECOMPOSITION_FACE%elementNumber=DECOMPOSITION_ELEMENT%NUMBER
!                              DOMAIN_FACE%elementNumber=DOMAIN_ELEMENT%NUMBER
                              DECOMPOSITION_face%xiNormalDirection=BASIS%localFaceXiNormal(basis_local_face_idx)
                              IF(ALLOCATED(BASIS%localFaceBasis)) THEN
                                DOMAIN_FACE%BASIS=>BASIS%localFaceBasis(basis_local_face_idx)%PTR
                              ELSE
                                !Basis is only 2D
                                DOMAIN_FACE%BASIS=>BASIS
                              ENDIF
                              ALLOCATE(DOMAIN_FACE%nodesInFace(BASIS%numberOfNodesInLocalFace(basis_local_face_idx)), &
                                & STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate face nodes in face",ERR,ERROR,*999)
                              ALLOCATE(DOMAIN_FACE%derivativesInFace(2,DOMAIN_FACE%BASIS%maximumNumberOfDerivatives, &
                                & BASIS%numberOfNodesInLocalFace(basis_local_face_idx)),STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate face derivatives in face",ERR,ERROR,*999)
                              DOMAIN_FACE%derivativesInFace=0
                              !Set nodes in face based upon face number
                              DOMAIN_FACE%nodesInFace(1:BASIS%numberOfNodesInLocalFace(basis_local_face_idx))= &
                                & TEMP_FACES(1:BASIS%numberOfNodesInLocalFace(basis_local_face_idx),FACE_NUMBER)
                              !Set derivatives of nodes in domain face from derivatives of nodes in element
                              DO basis_local_face_node_idx=1,BASIS%numberOfNodesInLocalFace(basis_local_face_idx)
                                element_local_node_idx=BASIS%nodeNumbersInLocalFace(basis_local_face_node_idx, &
                                  & basis_local_face_idx)
                                !Set derivative number of u (NO_GLOBAL_DERIV) for the domain face
                                DOMAIN_FACE%derivativesInFace(1,1,basis_local_face_node_idx)=NO_GLOBAL_DERIV
                                !Set version number of u (NO_GLOBAL_DERIV) for the domain face
                                version_idx=DOMAIN_ELEMENT%elementVersions(1,BASIS%nodeNumbersInLocalFace( &
                                  & basis_local_face_node_idx,basis_local_face_idx))
                                DOMAIN_FACE%derivativesInFace(2,1,basis_local_face_node_idx)=version_idx
                                IF(DOMAIN_FACE%BASIS%maximumNumberOfDerivatives>1) THEN
                                  DO basis_local_face_derivative_idx=2,DOMAIN_FACE%BASIS%maximumNumberOfDerivatives
                                    derivative_idx=DOMAIN_ELEMENT%elementDerivatives(BASIS%derivativeNumbersInLocalFace( &
                                      & basis_local_face_derivative_idx,basis_local_face_node_idx,basis_local_face_idx), &
                                      & element_local_node_idx)
                                    DOMAIN_FACE%derivativesInFace(1,basis_local_face_derivative_idx, &
                                      & basis_local_face_node_idx)=derivative_idx
                                    version_idx=DOMAIN_ELEMENT%elementVersions(BASIS%derivativeNumbersInLocalFace( &
                                      & basis_local_face_derivative_idx,basis_local_face_node_idx,basis_local_face_idx), &
                                      & element_local_node_idx)
                                    DOMAIN_FACE%derivativesInFace(2,basis_local_face_derivative_idx, &
                                      & basis_local_face_node_idx)=version_idx
                                  ENDDO !basis_local_face_derivative_idx
                                ENDIF
                              ENDDO !basis_local_face_node_idx
                            ENDIF
                          ENDDO !basis_local_face_idx
                        ENDDO !ne

                        DEALLOCATE(TEMP_FACES)
                        !\todo Note: Adjacency will be left out of faces calculation for the time being
                        !Calculate adjacent faces and the surrounding elements for each face
                        DO face_idx=1,DECOMPOSITION_FACES%numberOfFaces
                          DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(face_idx)
                          DOMAIN_FACE=>DOMAIN_FACES%FACES(face_idx)
                          BASIS=>DOMAIN_FACE%BASIS
                          IF(DECOMPOSITION_FACE%numberOfSurroundingElements==1) THEN
                            DECOMPOSITION_FACE%boundaryFace=.TRUE.
                            DOMAIN_FACE%boundaryFace=.TRUE.
                          ENDIF
                          !Allocate the elements surrounding the face
                          ALLOCATE(DECOMPOSITION_FACE%surroundingElements(DECOMPOSITION_FACE%numberOfSurroundingElements), &
                            & STAT=ERR)
                          IF(ERR/=0) CALL FlagError("Could not allocate face surrounding elements",ERR,ERROR,*999)

                          ALLOCATE(DECOMPOSITION_FACE%elementFaces(DECOMPOSITION_FACE%numberOfSurroundingElements), &
                            & STAT=ERR)
                          IF(ERR/=0) CALL FlagError("Could not allocate face element faces",ERR,ERROR,*999)
!                          DECOMPOSITION_FACE%numberOfSurroundingElements=0
!                          DECOMPOSITION_FACE%ADJACENT_FACES=0

                           !Loop over the nodes at each end of the face
!                          DO node_idx1=0,1
!                           DO node_idx2=0,1
!                            FOUND=.FALSE.
!                            node_idx=DOMAIN_FACE%nodesInFace((node_idx2*BASIS%numberOfNodes_IN_XI_DIRECTION*(BASIS%numberOfFaces-1))&
!                                                                             &+(node_idx1*(BASIS%numberOfNodes_IN_XI_DIRECTION-1))+1)
                             !Loop over the elements surrounding the node.
!                            DO elem_idx=1,DOMAIN_NODES%NODES(node_idx)%numberOfSurroundingElements
!                              ne=DOMAIN_NODES%NODES(node_idx)%surroundingElements(elem_idx)
!                              DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
!                              DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                               !Loop over the local faces of the element
!                              DO basis_local_face_idx=1,DOMAIN_ELEMENT%BASIS%numberOfLocalFaces
!                                nf2=DECOMPOSITION_ELEMENT%elementFaces(basis_local_face_idx)
!                                IF(nf2/=face_idx) THEN
!                                  DECOMPOSITION_FACE2=>DECOMPOSITION_FACES%FACES(nf2)
!                                  DOMAIN_FACE2=>DOMAIN_FACES%FACES(nf2)
                                   !Check whether XI of face have same direction
!                                  IF ((OTHER_XI_DIRECTIONS3(BASIS%localFaceXiNormal(basis_local_face_idx),2,1)==&
!                                     &OTHER_XI_DIRECTIONS3(BASIS2%localFaceXiNormal(basis_local_face_idx),2,1)).OR.&
!                                     &(OTHER_XI_DIRECTIONS3(BASIS%localFaceXiNormal(basis_local_face_idx),3,1)==&
!                                     &OTHER_XI_DIRECTIONS3(BASIS2%localFaceXiNormal(basis_local_face_idx),3,1))) THEN
                                     !Loop over nodes in face of surrounding element
!                                    BASIS2=>DOMAIN_FACE2%BASIS
!                                    IF(BASIS2%interpolationOrder(1)==BASIS%interpolationOrder(1)) THEN
!                                      NODE_COUNT=0
!                                      DO node_idx3=1,BASIS%numberOfNodes_IN_XI_DIRECTION
!                                        DO node_idx4=1,BASIS%numberOfNodes_IN_XI_DIRECTION
!                                          np2=DOMAIN_FACE2%nodesInFace((node_idx4*(BASIS2%numberOfFaces-1))&
!                                                                      &+(node_idx3*(BASIS2%numberOfNodes_IN_XI_DIRECTION-1))+1)
!                                          IF(np2==node_idx) NODE_COUNT=NODE_COUNT+1
!                                        ENDDO !node_idx4
!                                      ENDDO !node_idx3
!                                      IF(NODE_COUNT<BASIS%numberOfNodes) THEN
!                                        FOUND=.TRUE.
!                                        EXIT
!                                      ENDIF
!                                    ENDIF
!                                  ENDIF
!                                ENDIF
!                              ENDDO !basis_local_face_idx
!                                IF(FOUND) EXIT
!                            ENDDO !elem_idx
!                            IF(FOUND) DECOMPOSITION_FACE%ADJACENT_FACES(node_idx2)=nf2
!                           ENDDO !node_idx2
!                           IF(FOUND) DECOMPOSITION_FACE%ADJACENT_FACES(node_idx1)=nf2
!                          ENDDO !node_idx1
                        ENDDO !face_idx

                        !Set the surrounding elements
                        DO ne=1,DECOMPOSITION_ELEMENTS%totalNumberOfElements
                          DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
                          DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                          BASIS=>DOMAIN_ELEMENT%BASIS
                          DO basis_local_face_idx=1,BASIS%numberOfLocalFaces
                            FACE_NUMBER=DECOMPOSITION_ELEMENT%elementFaces(basis_local_face_idx)
                            DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(FACE_NUMBER)
                            DO face_idx=1,DECOMPOSITION_FACE%numberOfSurroundingElements
                              DECOMPOSITION_FACE%surroundingElements(face_idx)=ne
                              DECOMPOSITION_FACE%elementFaces(face_idx)=basis_local_face_idx
                            ENDDO
                          ENDDO !basis_local_face_idx
                        ENDDO !ne
                      ELSE
                        CALL FlagError("Domain topology faces is not associated",ERR,ERROR,*999)
                      ENDIF
                    CASE DEFAULT
                      CALL FlagError("Invalid number of dimensions for a topology domain",ERR,ERROR,*999)
                    END SELECT
                 ELSE
                    CALL FlagError("Domain topology elements is not associated",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Domain topology nodes is not associated",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Topology decomposition domain topology is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Topology decomposition domain is not associated",ERR,ERROR,*999)
            ENDIF
            !Now loop over the other mesh components in the decomposition and calculate the domain faces
            MESH=>DECOMPOSITION%MESH
            IF(ASSOCIATED(MESH)) THEN
              DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
                IF(component_idx/=decomposition%meshComponentNumber) THEN
                  DOMAIN=>DECOMPOSITION%DOMAIN(component_idx)%PTR
                  IF(ASSOCIATED(DOMAIN)) THEN
                    DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
                    IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                      DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
                      IF(ASSOCIATED(DOMAIN_NODES)) THEN
                        DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                        IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                          DOMAIN_FACES=>DOMAIN_TOPOLOGY%FACES
                          IF(ASSOCIATED(DOMAIN_FACES)) THEN
                            ALLOCATE(DOMAIN_FACES%FACES(DECOMPOSITION_FACES%numberOfFaces),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate domain faces faces",ERR,ERROR,*999)
                            DOMAIN_FACES%numberOfFaces=DECOMPOSITION_FACES%numberOfFaces
                            ALLOCATE(NODES_NUMBER_OF_FACES(DOMAIN_NODES%totalNumberOfNodes),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate nodes number of faces array",ERR,ERROR,*999)
                            NODES_NUMBER_OF_FACES=0
                            !Loop over the faces in the topology
                            DO face_idx=1,DECOMPOSITION_FACES%numberOfFaces
                              DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(face_idx)
                              DOMAIN_FACE=>DOMAIN_FACES%FACES(face_idx)
                              IF(DECOMPOSITION_FACE%numberOfSurroundingElements>0) THEN
                                ne=DECOMPOSITION_FACE%surroundingElements(1)
                                basis_local_face_idx=DECOMPOSITION_FACE%elementFaces(1)
                                CALL DOMAIN_TOPOLOGY_FACE_INITIALISE(DOMAIN_FACES%FACES(face_idx),ERR,ERROR,*999)
                                DOMAIN_FACE%NUMBER=face_idx
                                DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                                BASIS=>DOMAIN_ELEMENT%BASIS
                                IF(ALLOCATED(BASIS%localFaceBasis)) THEN
                                  DOMAIN_FACE%BASIS=>BASIS%localFaceBasis(basis_local_face_idx)%PTR
                                ELSE
                                  !Basis is only 2D
                                  DOMAIN_FACE%BASIS=>BASIS
                                ENDIF
                                ALLOCATE(DOMAIN_FACE%nodesInFace(BASIS%numberOfNodesInLocalFace(basis_local_face_idx)), &
                                  & STAT=ERR)
                                IF(ERR/=0) CALL FlagError("Could not allocate nodes in face",ERR,ERROR,*999)
                                ALLOCATE(DOMAIN_FACE%derivativesInFace(2,DOMAIN_FACE%BASIS%maximumNumberOfDerivatives, &
                                  & BASIS%numberOfNodesInLocalFace(basis_local_face_idx)),STAT=ERR)
                                IF(ERR/=0) CALL FlagError("Could not allocate derivatives in face",ERR,ERROR,*999)
                                !Set derivatives of nodes in domain face from derivatives of nodes in element
                                DO basis_local_face_node_idx=1,BASIS%numberOfNodesInLocalFace(basis_local_face_idx)
                                  element_local_node_idx=BASIS%nodeNumbersInLocalFace(basis_local_face_node_idx, &
                                    & basis_local_face_idx)
                                  node_idx=DOMAIN_ELEMENT%elementNodes(element_local_node_idx)
                                  DOMAIN_FACE%nodesInFace(basis_local_face_node_idx)=node_idx
                                  !Set derivative number of u (NO_GLOBAL_DERIV) for the domain face
                                  DOMAIN_FACE%derivativesInFace(1,1,basis_local_face_node_idx)=NO_GLOBAL_DERIV
                                  !Set version number of u (NO_GLOBAL_DERIV) for the domain face
                                  version_idx=DOMAIN_ELEMENT%elementVersions(1,BASIS%nodeNumbersInLocalFace( &
                                    & basis_local_face_node_idx,basis_local_face_idx))
                                  DOMAIN_FACE%derivativesInFace(2,1,basis_local_face_node_idx)=version_idx
                                  IF(DOMAIN_FACE%BASIS%maximumNumberOfDerivatives>1) THEN
                                    DO basis_local_face_derivative_idx=2,DOMAIN_FACE%BASIS%maximumNumberOfDerivatives
                                      derivative_idx=DOMAIN_ELEMENT%elementDerivatives(BASIS%derivativeNumbersInLocalFace( &
                                        & basis_local_face_derivative_idx,basis_local_face_node_idx,basis_local_face_idx), &
                                        & element_local_node_idx)
                                      DOMAIN_FACE%derivativesInFace(1,basis_local_face_derivative_idx, &
                                        & basis_local_face_node_idx)=derivative_idx
                                      version_idx=DOMAIN_ELEMENT%elementVersions(BASIS%derivativeNumbersInLocalFace( &
                                        & basis_local_face_derivative_idx,basis_local_face_node_idx,basis_local_face_idx), &
                                        & element_local_node_idx)
                                      DOMAIN_FACE%derivativesInFace(2,basis_local_face_derivative_idx, &
                                        & basis_local_face_node_idx)=version_idx
                                    ENDDO !basis_local_face_derivative_idx
                                  ENDIF
                                  NODES_NUMBER_OF_FACES(node_idx)=NODES_NUMBER_OF_FACES(node_idx)+1
                                ENDDO !basis_local_face_node_idx
                              ELSE
                                CALL FlagError("Face is not surrounded by any elements?",ERR,ERROR,*999)
                              ENDIF
                            ENDDO !face_idx
                            DO node_idx=1,DOMAIN_NODES%totalNumberOfNodes
                              ALLOCATE(DOMAIN_NODES%NODES(node_idx)%nodeFaces(NODES_NUMBER_OF_FACES(node_idx)),STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate node faces",ERR,ERROR,*999)
                              DOMAIN_NODES%NODES(node_idx)%numberOfNodeFaces=0
                            ENDDO !node_idx
                            DEALLOCATE(NODES_NUMBER_OF_FACES)
                            DO face_idx=1,DOMAIN_FACES%numberOfFaces
                              DOMAIN_FACE=>DOMAIN_FACES%FACES(face_idx)
                              BASIS=>DOMAIN_FACE%BASIS
                              DO basis_local_face_node_idx=1,BASIS%numberOfNodes
                                node_idx=DOMAIN_FACE%nodesInFace(basis_local_face_node_idx)
                                DOMAIN_NODE=>DOMAIN_NODES%NODES(node_idx)
                                DOMAIN_NODE%numberOfNodeFaces=DOMAIN_NODE%numberOfNodeFaces+1
                                !Set the face numbers a node is on
                                DOMAIN_NODE%nodeFaces(DOMAIN_NODE%numberOfNodeFaces)=face_idx
                              ENDDO !basis_local_face_node_idx
                            ENDDO !face_idx
                          ELSE
                            CALL FlagError("Domain faces is not associated",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Domain elements is not associated",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Domain nodes is not associated",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Domain topology is not associated",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Decomposition mesh is not associated",ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ENDDO !component_idx
            ELSE
              CALL FlagError("Decomposition mesh is not associated",ERR,ERROR,*999)
            ENDIF                        
          ELSE
            CALL FlagError("Topology decomposition is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Topology decomposition elements is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Topology faces is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Decomposition topology faces:",ERR,ERROR,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of mesh components = ",MESH%NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of faces = ",DECOMPOSITION_FACES%numberOfFaces,ERR,ERROR,*999)
      DO face_idx=1,DECOMPOSITION_FACES%numberOfFaces
        DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(face_idx)
        DOMAIN_FACE=>DOMAIN_FACES%FACES(face_idx)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Face number = ",DECOMPOSITION_FACE%NUMBER,ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Xi direction (Normal to Face) = &
                                                                         &",DECOMPOSITION_face%xiNormalDirection,ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of surrounding elements = ", &
          & DECOMPOSITION_FACE%numberOfSurroundingElements,ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_FACE%numberOfSurroundingElements,4,4, &
          & DECOMPOSITION_FACE%surroundingElements,'("      Surrounding elements :",4(X,I8))','(28X,4(X,I8))',ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_FACE%numberOfSurroundingElements,4,4, &
          & DECOMPOSITION_FACE%elementFaces,'("      Element faces        :",4(X,I8))','(28X,4(X,I8))',ERR,ERROR,*999)
!        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,2,2,DECOMPOSITION_FACE%ADJACENT_FACES, &
!          & '("      Adjacent faces       :",2(X,I8))','(28X,2(X,I8))',ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Boundary face = ",DECOMPOSITION_FACE%boundaryFace,ERR,ERROR,*999)
        DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Mesh component : ",component_idx,ERR,ERROR,*999)
          DOMAIN=>DECOMPOSITION%DOMAIN(component_idx)%PTR
          DOMAIN_FACE=>DOMAIN%TOPOLOGY%FACES%FACES(face_idx)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Basis user number = ",DOMAIN_FACE%BASIS%userNumber, &
            & ERR,ERROR,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Basis family number = ",DOMAIN_FACE%BASIS%familyNumber, &
            & ERR,ERROR,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation type = ",DOMAIN_FACE%BASIS% &
            & interpolationType(1),ERR,ERROR,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation order = ",DOMAIN_FACE%BASIS% &
            & interpolationOrder(1),ERR,ERROR,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of nodes in faces = ",DOMAIN_FACE%BASIS%numberOfNodes, &
            & ERR,ERROR,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_FACE%BASIS%numberOfNodes,4,4,DOMAIN_FACE%nodesInFace, &
            & '("        Nodes in face        :",4(X,I8))','(30X,4(X,I8))',ERR,ERROR,*999)
          DO basis_local_face_node_idx=1,DOMAIN_FACE%BASIS%numberOfNodes
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Node : ",basis_local_face_node_idx,ERR,ERROR,*999)
            !/TODO::Loop over local_derivative index so this output makes more sense !<derivativesInLine(i,local_derivative_idx,local_node_idx)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & DOMAIN_FACE%BASIS%numberOfDerivatives(basis_local_face_node_idx),4,4,DOMAIN_FACE% &
              & derivativesInFace(1,:,basis_local_face_node_idx),'("            Derivatives in face  :",4(X,I8))', &
              & '(34X,4(X,I8))',ERR,ERROR,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & DOMAIN_FACE%BASIS%numberOfDerivatives(basis_local_face_node_idx),4,4,DOMAIN_FACE% &
              & derivativesInFace(2,:,basis_local_face_node_idx),'("            Derivatives Versions in face  :",4(X,I8))', &
              & '(34X,4(X,I8))',ERR,ERROR,*999)
          ENDDO !basis_local_face_node_idx
        ENDDO !component_idx
      ENDDO !face_idx
    ENDIF

    EXITS("DECOMPOSITION_TOPOLOGY_FACES_CALCULATE")
    RETURN
999 IF(ASSOCIATED(TEMP_FACES)) DEALLOCATE(TEMP_FACES)
    IF(ASSOCIATED(NEW_TEMP_FACES)) DEALLOCATE(NEW_TEMP_FACES)
    IF(ALLOCATED(NODES_NUMBER_OF_FACES)) DEALLOCATE(NODES_NUMBER_OF_FACES)
    ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FACES_CALCULATE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalises the faces in the given decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nf
    
    ENTERS("DECOMPOSITION_TOPOLOGY_FACES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%FACES)) THEN
        DO nf=1,TOPOLOGY%FACES%numberOfFaces
          CALL DECOMPOSITION_TOPOLOGY_FACE_FINALISE(TOPOLOGY%FACES%FACES(nf),ERR,ERROR,*999)
        ENDDO !nf
        IF(ALLOCATED(TOPOLOGY%FACES%FACES)) DEALLOCATE(TOPOLOGY%FACES%FACES)
        DEALLOCATE(TOPOLOGY%FACES)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DECOMPOSITION_TOPOLOGY_FACES_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FACES_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the face data structures for a decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to initialise the faces for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_FACES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%FACES)) THEN
        CALL FlagError("Decomposition already has topology faces associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%FACES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology faces",ERR,ERROR,*999)
        TOPOLOGY%FACES%numberOfFaces=0
        TOPOLOGY%FACES%DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_FACES_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FACES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_INITIALISE

  !
  !================================================================================================================================
  !
  
  !>Gets the decomposition type for a decomposition. \see OpenCMISS::Iron::cmfe_DecompositionTypeGet
  SUBROUTINE DECOMPOSITION_TYPE_GET(DECOMPOSITION,TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: DECOMPOSITION !<A pointer to the decomposition to get the type for
    INTEGER(INTG), INTENT(OUT) :: TYPE !<On return, the decomposition type for the specified decomposition \see MESH_ROUTINES_DecompositionTypes,MESH_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%decompositionFinished) THEN
        TYPE=DECOMPOSITION%domainDecompositionType
      ELSE
        CALL FlagError("Decomposition has not finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TYPE_GET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TYPE_GET",ERR,ERROR)
    RETURN
  END SUBROUTINE DECOMPOSITION_TYPE_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the decomposition type for a decomposition.  \see OpenCMISS::Iron::cmfe_DecompositionTypeSet
  SUBROUTINE DECOMPOSITION_TYPE_SET(DECOMPOSITION,TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the type for
    INTEGER(INTG), INTENT(IN) :: TYPE !<The decomposition type to set \see MESH_ROUTINES_DecompositionTypes,MESH_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%decompositionFinished) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(TYPE)
        CASE(DECOMPOSITION_ALL_TYPE)
          !heye: three types for decomposition--decompostion_all_type means no decomposition
          DECOMPOSITION%domainDecompositionType=DECOMPOSITION_ALL_TYPE
        CASE(DECOMPOSITION_CALCULATED_TYPE)
          DECOMPOSITION%domainDecompositionType=DECOMPOSITION_CALCULATED_TYPE
        CASE(DECOMPOSITION_USER_DEFINED_TYPE)
          DECOMPOSITION%domainDecompositionType=DECOMPOSITION_USER_DEFINED_TYPE
        CASE DEFAULT
          LOCAL_ERROR="Decomposition type "//TRIM(NUMBER_TO_VSTRING(TYPE,"*",ERR,ERROR))//" is not valid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TYPE_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TYPE_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes whether lines should be calculated in the the decomposition. \see OpenCMISS::Iron::cmfe_DecompositionCalculateLinesSet
  SUBROUTINE DECOMPOSITION_CALCULATE_LINES_SET(DECOMPOSITION,CALCULATE_LINES_FLAG,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: DECOMPOSITION !<A pointer to the decomposition
    LOGICAL, INTENT(IN) :: CALCULATE_LINES_FLAG !<The boolean flag to determine whether the lines should be calculated or not
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("DECOMPOSITION_CALCULATE_LINES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%decompositionFinished) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        DECOMPOSITION%calculateLines=CALCULATE_LINES_FLAG
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_CALCULATE_LINES_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_CALCULATE_LINES_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_CALCULATE_LINES_SET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes whether faces should be calculated in the the decomposition. \see OpenCMISS::Iron::cmfe_DecompositionCalculateFacesSet
  SUBROUTINE DECOMPOSITION_CALCULATE_FACES_SET(DECOMPOSITION,CALCULATE_FACES_FLAG,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: DECOMPOSITION !<A pointer to the decomposition
    LOGICAL, INTENT(IN) :: CALCULATE_FACES_FLAG !<The boolean flag to determine whether the faces should be calculated or not
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("DECOMPOSITION_CALCULATE_FACES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%decompositionFinished) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        DECOMPOSITION%calculateFaces=CALCULATE_FACES_FLAG
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_CALCULATE_FACES_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_CALCULATE_FACES_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_CALCULATE_FACES_SET
  
  !
  !================================================================================================================================
  !

  !>Finalises the domain decompositions for a given mesh. \todo Pass in a pointer to the decompositions.
  SUBROUTINE DECOMPOSITIONS_FINALISE(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshType), POINTER :: MESH !<A pointer to the mesh to finalise the decomposition for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITIONS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN
        DO WHILE(MESH%DECOMPOSITIONS%numberOfDecompositions>0)
          CALL DECOMPOSITION_DESTROY(MESH%DECOMPOSITIONS%DECOMPOSITIONS(1)%PTR,ERR,ERROR,*999)
        ENDDO !no_decomposition
       DEALLOCATE(MESH%DECOMPOSITIONS)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITIONS_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITIONS_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITIONS_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the domain decompositions for a given mesh.
  SUBROUTINE DECOMPOSITIONS_INITIALISE(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshType), POINTER :: MESH !<A pointer to the mesh to initialise the decompositions for

    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITIONS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN
        CALL FlagError("Mesh already has decompositions associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(MESH%DECOMPOSITIONS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Mesh decompositions could not be allocated.",ERR,ERROR,*999)
        MESH%DECOMPOSITIONS%numberOfDecompositions=0
        MESH%DECOMPOSITIONS%MESH=>MESH
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITIONS_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITIONS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITIONS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the domain for a given decomposition and deallocates all memory. \todo Pass in a pointer to the domain
  SUBROUTINE DOMAIN_FINALISE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: DECOMPOSITION !<A pointer to the decomposition to finalise the domain for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx
    
    ENTERS("DOMAIN_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
        IF(ASSOCIATED(DECOMPOSITION%DOMAIN)) THEN
          DO component_idx=1,SIZE(DECOMPOSITION%DOMAIN,1)
            IF(ALLOCATED(DECOMPOSITION%DOMAIN(component_idx)%PTR%nodeDomain))  &
              & DEALLOCATE(DECOMPOSITION%DOMAIN(component_idx)%PTR%nodeDomain)
            CALL DOMAIN_MAPPINGS_FINALISE(DECOMPOSITION%DOMAIN(component_idx)%PTR,ERR,ERROR,*999)        
            CALL DOMAIN_TOPOLOGY_FINALISE(DECOMPOSITION%DOMAIN(component_idx)%PTR,ERR,ERROR,*999)
            DEALLOCATE(DECOMPOSITION%DOMAIN(component_idx)%PTR)
          ENDDO !component_idx
          DEALLOCATE(DECOMPOSITION%DOMAIN)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the domain for a given decomposition.
  SUBROUTINE DOMAIN_INITIALISE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: DECOMPOSITION !<A pointer to the decomposition to initialise the domain for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx

    ENTERS("DOMAIN_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
        IF(ASSOCIATED(DECOMPOSITION%DOMAIN)) THEN
          CALL FlagError("Decomposition already has a domain associated.",ERR,ERROR,*999)
        ELSE
          ALLOCATE(DECOMPOSITION%DOMAIN(DECOMPOSITION%numberOfComponents),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Decomposition domain could not be allocated.",ERR,ERROR,*999)
          DO component_idx=1,DECOMPOSITION%numberOfComponents!Mesh component
            ALLOCATE(DECOMPOSITION%DOMAIN(component_idx)%PTR,STAT=ERR)
            IF(ERR/=0) CALL FlagError("Decomposition domain component could not be allocated.",ERR,ERROR,*999)
            DECOMPOSITION%DOMAIN(component_idx)%PTR%DECOMPOSITION=>DECOMPOSITION
            DECOMPOSITION%DOMAIN(component_idx)%PTR%MESH=>DECOMPOSITION%MESH
            DECOMPOSITION%DOMAIN(component_idx)%PTR%meshComponentNumber=component_idx
            DECOMPOSITION%DOMAIN(component_idx)%PTR%REGION=>DECOMPOSITION%region
            DECOMPOSITION%DOMAIN(component_idx)%PTR%numberOfDimensions=DECOMPOSITION%numberOfDimensions
            !DECOMPOSITION%DOMAIN(component_idx)%PTR%numberOfElements=0
            !DECOMPOSITION%DOMAIN(component_idx)%PTR%numberOfFaces=0
            !DECOMPOSITION%DOMAIN(component_idx)%PTR%numberOfLines=0
            !DECOMPOSITION%DOMAIN(component_idx)%PTR%numberOfNodes=0
            !DECOMPOSITION%DOMAIN(component_idx)%PTR%NUMBER_OF_MESH_DOFS=0
            NULLIFY(DECOMPOSITION%DOMAIN(component_idx)%PTR%MAPPINGS)
            NULLIFY(DECOMPOSITION%DOMAIN(component_idx)%PTR%TOPOLOGY)
            CALL DOMAIN_MAPPINGS_INITIALISE(DECOMPOSITION%DOMAIN(component_idx)%PTR,ERR,ERROR,*999)
            CALL DOMAIN_TOPOLOGY_INITIALISE(DECOMPOSITION%DOMAIN(component_idx)%PTR,ERR,ERROR,*999)
          ENDDO !component_idx
        ENDIF
      ELSE
        CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_INITIALISE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE DOMAIN_INITIALISE
  
 !
  !================================================================================================================================
  !

  !>Finalises the dofs mapping in the given domain mappings. \todo Pass in the domain mappings dofs
  SUBROUTINE DOMAIN_MAPPINGS_DOFS_FINALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainMappingsType), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mappings to finalise the dofs for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("DOMAIN_MAPPINGS_DOFS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%DOFS)) THEN
        CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(DOMAIN_MAPPINGS%DOFS,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_DOFS_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_DOFS_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_DOFS_FINALISE

  !
  !================================================================================================================================
  !

  !>Intialises the dofs mapping in the given domain mapping.
  SUBROUTINE DOMAIN_MAPPINGS_DOFS_INITIALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainMappingsType), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mappings to initialise the dofs for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("DOMAIN_MAPPINGS_DOFS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%DOFS)) THEN
        CALL FlagError("Domain dofs mappings are already associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(DOMAIN_MAPPINGS%DOFS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate domain mappings dofs.",ERR,ERROR,*999)
        CALL DomainMappings_MappingInitialise(DOMAIN_MAPPINGS%DOMAIN%DECOMPOSITION%workGroup, &
          & DOMAIN_MAPPINGS%DOFS,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_DOFS_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_DOFS_INITIALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_DOFS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Calculates the local/global element mappings for a domain decomposition.
  SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_CALCULATE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainType), POINTER :: DOMAIN !<A pointer to the domain to calculate the element mappings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,no_adjacent_element,adjacent_element,domain_no,domain_idx,ne,nn,np,numberOfDomains, &
      & numberOfAdjacentElements,component_idx
    INTEGER(INTG), ALLOCATABLE :: adjacentElements(:),DOMAINS(:),LOCAL_ELEMENT_NUMBERS(:)
    TYPE(LIST_TYPE), POINTER :: adjacentDomainsList
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: ADJACENT_ELEMENTS_LIST(:)
    TYPE(BasisType), POINTER :: BASIS
    TYPE(MeshType), POINTER :: MESH
    TYPE(DecompositionType), POINTER :: DECOMPOSITION
    TYPE(DomainMappingType), POINTER :: ELEMENTS_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS("DOMAIN_MAPPINGS_ELEMENTS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%MAPPINGS)) THEN
        IF(ASSOCIATED(DOMAIN%MAPPINGS%ELEMENTS)) THEN
          ELEMENTS_MAPPING=>DOMAIN%MAPPINGS%ELEMENTS
          IF(ASSOCIATED(DOMAIN%DECOMPOSITION)) THEN
            DECOMPOSITION=>DOMAIN%DECOMPOSITION
            IF(ASSOCIATED(DOMAIN%MESH)) THEN
              MESH=>DOMAIN%MESH
              component_idx=DOMAIN%meshComponentNumber
                            
              !Calculate the local and global numbers and set up the mappings                           
              ALLOCATE(ELEMENTS_MAPPING%globalToLocalMap(MESH%numberOfElements),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate element mapping global to local map.",ERR,ERROR,*999)
              ELEMENTS_MAPPING%numberOfGlobal=MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%numberOfElements
              !Loop over the global elements and calculate local numbers
              ALLOCATE(LOCAL_ELEMENT_NUMBERS(0:decomposition%numberOfDomains-1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate local element numbers.",ERR,ERROR,*999)
              LOCAL_ELEMENT_NUMBERS=0
              ALLOCATE(ADJACENT_ELEMENTS_LIST(0:decomposition%numberOfDomains-1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate adjacent elements list.",ERR,ERROR,*999)
              DO domain_idx=0,decomposition%numberOfDomains-1
                NULLIFY(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR)
                CALL LIST_CREATE_START(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_INITIAL_SIZE_SET(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,MAX(INT(MESH%numberOfElements/2),1), &
                  & ERR,ERROR,*999)
                CALL LIST_CREATE_FINISH(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,ERR,ERROR,*999)
              ENDDO !domain_idx
            
              DO ne=1,MESH%numberOfElements
                !Calculate the local numbers
                domain_no=DECOMPOSITION%elementDomain(ne)
                LOCAL_ELEMENT_NUMBERS(domain_no)=LOCAL_ELEMENT_NUMBERS(domain_no)+1
                !Calculate the adjacent elements to the computation domains and the adjacent domain numbers themselves
                BASIS=>MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)%BASIS
                NULLIFY(adjacentDomainsList)
                CALL LIST_CREATE_START(adjacentDomainsList,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(adjacentDomainsList,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_INITIAL_SIZE_SET(adjacentDomainsList,decomposition%numberOfDomains,ERR,ERROR,*999)
                CALL LIST_CREATE_FINISH(adjacentDomainsList,ERR,ERROR,*999)
                CALL LIST_ITEM_ADD(adjacentDomainsList,domain_no,ERR,ERROR,*999)
                DO nn=1,BASIS%numberOfNodes
                  np=MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)%meshElementNodes(nn)
                  DO no_adjacent_element=1,MESH%TOPOLOGY(component_idx)%PTR%NODES%NODES(np)%numberOfSurroundingElements
                    adjacent_element=MESH%TOPOLOGY(component_idx)%PTR%NODES%NODES(np)%surroundingElements(no_adjacent_element)
                    IF(DECOMPOSITION%elementDomain(adjacent_element)/=domain_no) THEN
                      CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(domain_no)%PTR,adjacent_element,ERR,ERROR,*999)
                      CALL LIST_ITEM_ADD(adjacentDomainsList,DECOMPOSITION%elementDomain(adjacent_element),ERR,ERROR,*999)
                    ENDIF
                  ENDDO !no_adjacent_element
                ENDDO !nn
                CALL LIST_REMOVE_DUPLICATES(adjacentDomainsList,ERR,ERROR,*999)
                CALL LIST_DETACH_AND_DESTROY(adjacentDomainsList,numberOfDomains,DOMAINS,ERR,ERROR,*999)
                DEALLOCATE(DOMAINS)
                CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(ELEMENTS_MAPPING%globalToLocalMap(ne),ERR,ERROR,*999)
                ALLOCATE(ELEMENTS_MAPPING%globalToLocalMap(ne)%localNumber(numberOfDomains),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate element global to local map local number.",ERR,ERROR,*999)
                ALLOCATE(ELEMENTS_MAPPING%globalToLocalMap(ne)%domainNumber(numberOfDomains),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate element global to local map domain number.",ERR,ERROR,*999)
                ALLOCATE(ELEMENTS_MAPPING%globalToLocalMap(ne)%localType(numberOfDomains),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate element global to local map local type.",ERR,ERROR,*999)
                ELEMENTS_MAPPING%globalToLocalMap(ne)%numberOfDomains=1
                ELEMENTS_MAPPING%globalToLocalMap(ne)%localNumber(1)=LOCAL_ELEMENT_NUMBERS(domain_no)
                ELEMENTS_MAPPING%globalToLocalMap(ne)%domainNumber(1)=DECOMPOSITION%elementDomain(ne) 
                IF(numberOfDomains==1) THEN
                  !Element is an internal element
                  ELEMENTS_MAPPING%globalToLocalMap(ne)%localType(1)=DOMAIN_LOCAL_INTERNAL
                ELSE
                  !Element is on the boundary of computation domains
                  ELEMENTS_MAPPING%globalToLocalMap(ne)%localType(1)=DOMAIN_LOCAL_BOUNDARY
                ENDIF
              ENDDO !ne
              
              !Compute ghost element mappings
              DO domain_idx=0,decomposition%numberOfDomains-1
                CALL LIST_REMOVE_DUPLICATES(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                CALL LIST_DETACH_AND_DESTROY(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,numberOfAdjacentElements, &
                  & adjacentElements,ERR,ERROR,*999)
                DO no_adjacent_element=1,numberOfAdjacentElements
                  adjacent_element=adjacentElements(no_adjacent_element)
                  LOCAL_ELEMENT_NUMBERS(domain_idx)=LOCAL_ELEMENT_NUMBERS(domain_idx)+1
                  ELEMENTS_MAPPING%globalToLocalMap(adjacent_element)%numberOfDomains= &
                    & ELEMENTS_MAPPING%globalToLocalMap(adjacent_element)%numberOfDomains+1
                  ELEMENTS_MAPPING%globalToLocalMap(adjacent_element)%localNumber( &
                    & ELEMENTS_MAPPING%globalToLocalMap(adjacent_element)%numberOfDomains)=LOCAL_ELEMENT_NUMBERS(domain_idx)
                  ELEMENTS_MAPPING%globalToLocalMap(adjacent_element)%domainNumber( &
                    & ELEMENTS_MAPPING%globalToLocalMap(adjacent_element)%numberOfDomains)=domain_idx
                  ELEMENTS_MAPPING%globalToLocalMap(adjacent_element)%localType( &
                    & ELEMENTS_MAPPING%globalToLocalMap(adjacent_element)%numberOfDomains)= &
                    & DOMAIN_LOCAL_GHOST
                ENDDO !no_adjacent_element
                IF(ALLOCATED(adjacentElements)) DEALLOCATE(adjacentElements)
              ENDDO !domain_idx
              
              DEALLOCATE(ADJACENT_ELEMENTS_LIST)
              DEALLOCATE(LOCAL_ELEMENT_NUMBERS)
              
              !Calculate element local to global maps from global to local map
              CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(ELEMENTS_MAPPING,ERR,ERROR,*999)
                            
            ELSE
              CALL FlagError("Domain mesh is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Domain decomposition is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Domain mappings elements is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Domain mappings is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain is not associated.",ERR,ERROR,*998)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Element mappings :",ERR,ERROR,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Global to local map :",ERR,ERROR,*999)
      DO ne=1,MESH%numberOfElements
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global element = ",ne,ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of domains  = ", &
          & ELEMENTS_MAPPING%globalToLocalMap(ne)%numberOfDomains,ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%globalToLocalMap(ne)% &
          & numberOfDomains,8,8,ELEMENTS_MAPPING%globalToLocalMap(ne)%localNumber, &
          & '("      Local number :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%globalToLocalMap(ne)% &
          & numberOfDomains,8,8,ELEMENTS_MAPPING%globalToLocalMap(ne)%domainNumber, &
          & '("      Domain number:",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%globalToLocalMap(ne)% &
          & numberOfDomains,8,8,ELEMENTS_MAPPING%globalToLocalMap(ne)%localType, &
          & '("      Local type   :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
      ENDDO !ne
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map :",ERR,ERROR,*999)
      DO ne=1,ELEMENTS_MAPPING%totalNumberOfLocal
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Local element = ",ne,ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Global element = ", &
          & ELEMENTS_MAPPING%localToGlobalMap(ne),ERR,ERROR,*999)
      ENDDO !ne
      IF(DIAGNOSTICS2) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Internal elements :",ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of internal elements = ", &
          & ELEMENTS_MAPPING%numberOfInternal,ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%numberOfInternal,8,8, &
          & ELEMENTS_MAPPING%domainList(ELEMENTS_MAPPING%internalStart:ELEMENTS_MAPPING%internalFinish), &
          & '("    Internal elements:",8(X,I7))','(22X,8(X,I7))',ERR,ERROR,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary elements :",ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary elements = ", &
          & ELEMENTS_MAPPING%numberOfBoundary,ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%numberOfBoundary,8,8, &
          & ELEMENTS_MAPPING%domainList(ELEMENTS_MAPPING%boundaryStart:ELEMENTS_MAPPING%boundaryFinish), &
          & '("    Boundary elements:",8(X,I7))','(22X,8(X,I7))',ERR,ERROR,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Ghost elements :",ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost elements = ", &
          & ELEMENTs_MAPPING%numberOfGhost,ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%numberOfGhost,8,8, &
          & ELEMENTS_MAPPING%domainList(ELEMENTS_MAPPING%ghostStart:ELEMENTS_MAPPING%ghostFinish), &
          & '("    Ghost elements   :",8(X,I7))','(22X,8(X,I7))',ERR,ERROR,*999)
      ENDIF
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains :",ERR,ERROR,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
        & ELEMENTS_MAPPING%numberOfAdjacentDomains,ERR,ERROR,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%numberOfDomains+1,8,8, &
        & ELEMENTS_MAPPING%adjacentDomainsPtr,'("    Adjacent domains ptr  :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%adjacentDomainsPtr( &
        & ELEMENTS_MAPPING%numberOfDomains)-1,8,8,ELEMENTS_MAPPING%adjacentDomainsList, &
        '("    Adjacent domains list :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
      DO domain_idx=1,ELEMENTS_MAPPING%numberOfAdjacentDomains
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",domain_idx,ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
          & ELEMENTS_MAPPING%adjacentDomains(domain_idx)%domainNumber,ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
          & ELEMENTS_MAPPING%adjacentDomains(domain_idx)%numberOfSendGhosts,ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%adjacentDomains(domain_idx)% &
          & numberOfSendGhosts,6,6,ELEMENTS_MAPPING%adjacentDomains(domain_idx)%localGhostSendIndices, &
          & '("      Local send ghost indicies       :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)      
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of recieve ghosts = ", &
          & ELEMENTS_MAPPING%adjacentDomains(domain_idx)%numberOfReceiveGhosts,ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%adjacentDomains(domain_idx)% &
          & numberOfReceiveGhosts,6,6,ELEMENTS_MAPPING%adjacentDomains(domain_idx)%localGhostReceiveIndices, &
          & '("      Local receive ghost indicies    :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)              
      ENDDO !domain_idx
    ENDIF
    
    EXITS("DOMAIN_MAPPINGS_ELEMENTS_CALCULATE")
    RETURN
999 IF(ALLOCATED(DOMAINS)) DEALLOCATE(DOMAINS)
    IF(ALLOCATED(adjacentElements)) DEALLOCATE(adjacentElements)    
    IF(ASSOCIATED(DOMAIN%MAPPINGS%ELEMENTS)) CALL DOMAIN_MAPPINGS_ELEMENTS_FINALISE(DOMAIN%MAPPINGS,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("DOMAIN_MAPPINGS_ELEMENTS_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Finalises the mappings in the given domain. \todo pass in the domain mappings
  SUBROUTINE DOMAIN_MAPPINGS_FINALISE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainType), POINTER :: DOMAIN !<A pointer to the domain to finalise the mappings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_MAPPINGS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      CALL DOMAIN_MAPPINGS_ELEMENTS_FINALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
      CALL DOMAIN_MAPPINGS_NODES_FINALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
      CALL DOMAIN_MAPPINGS_DOFS_FINALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
      DEALLOCATE(DOMAIN%MAPPINGS)
    ELSE
      CALL FlagError("Domain is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_FINALISE

  !
  !================================================================================================================================
  !

  !>Finalises the element mapping in the given domain mapping. \todo pass in the domain mappings elements
  SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_FINALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainMappingsType), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mappings to finalise the elements for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("DOMAIN_MAPPINGS_ELEMENTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%ELEMENTS)) THEN
        CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(DOMAIN_MAPPINGS%ELEMENTS,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_ELEMENTS_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_ELEMENTS_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the element mapping in the given domain mapping.
  SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_INITIALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainMappingsType), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mappings to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("DOMAIN_MAPPINGS_ELEMENTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%ELEMENTS)) THEN
        CALL FlagError("Domain elements mappings are already associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(DOMAIN_MAPPINGS%ELEMENTS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate domain mappings elements.",ERR,ERROR,*999)
        CALL DomainMappings_MappingInitialise(DOMAIN_MAPPINGS%DOMAIN%DECOMPOSITION%workGroup, &
          & DOMAIN_MAPPINGS%ELEMENTS,ERR,ERROR,*999)
       ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_ELEMENTS_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_ELEMENTS_INITIALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the mappings for a domain decomposition. \todo finalise on error.
  SUBROUTINE DOMAIN_MAPPINGS_INITIALISE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainType), POINTER :: DOMAIN !<A pointer to the domain to initialise the mappings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_MAPPINGS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%MAPPINGS)) THEN
        CALL FlagError("Domain already has mappings associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(DOMAIN%MAPPINGS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate domain mappings.",ERR,ERROR,*999)
        DOMAIN%MAPPINGS%DOMAIN=>DOMAIN
        NULLIFY(DOMAIN%MAPPINGS%ELEMENTS)
        NULLIFY(DOMAIN%MAPPINGS%NODES)
        NULLIFY(DOMAIN%MAPPINGS%DOFS)
        !Calculate the node and element mappings
        CALL DOMAIN_MAPPINGS_ELEMENTS_INITIALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_NODES_INITIALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_DOFS_INITIALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_ELEMENTS_CALCULATE(DOMAIN,ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE(DOMAIN,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_MAPPINGS_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Calculates the local/global node and dof mappings for a domain decomposition.
  SUBROUTINE DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainType), POINTER :: DOMAIN !<A pointer to the domain to calculate the node dofs for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,no_adjacent_element,no_computation_node,no_ghost_node,adjacent_element,ghost_node, &
      & NUMBER_OF_NODES_PER_DOMAIN,domain_idx,domain_idx2,domain_no,node_idx,derivative_idx,version_idx,ny,numberOfDomains, &
      & MAX_NUMBER_DOMAINS,NUMBER_OF_GHOST_NODES,myGroupComputationNodeNumber,numberOfGroupComputationNodes,component_idx
    INTEGER(INTG), ALLOCATABLE :: LOCAL_NODE_NUMBERS(:),LOCAL_DOF_NUMBERS(:),NODE_COUNT(:),NUMBER_INTERNAL_NODES(:), &
      & NUMBER_BOUNDARY_NODES(:)
    INTEGER(INTG), ALLOCATABLE :: DOMAINS(:),ALL_DOMAINS(:),GHOST_NODES(:)
    LOGICAL :: BOUNDARY_DOMAIN
    TYPE(LIST_TYPE), POINTER :: adjacentDomainsList,ALL_ADJACENT_DOMAINS_LIST
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: GHOST_NODES_LIST(:)
    TYPE(MeshType), POINTER :: MESH
    TYPE(MeshComponentTopologyType), POINTER :: MESH_TOPOLOGY
    TYPE(DecompositionType), POINTER :: DECOMPOSITION
    TYPE(DomainMappingType), POINTER :: ELEMENTS_MAPPING
    TYPE(DomainMappingType), POINTER :: NODES_MAPPING
    TYPE(DomainMappingType), POINTER :: DOFS_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    ENTERS("DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%MAPPINGS)) THEN
        IF(ASSOCIATED(DOMAIN%MAPPINGS%NODES)) THEN
          NODES_MAPPING=>DOMAIN%MAPPINGS%NODES
          IF(ASSOCIATED(DOMAIN%MAPPINGS%DOFS)) THEN
            DOFS_MAPPING=>DOMAIN%MAPPINGS%DOFS
            IF(ASSOCIATED(DOMAIN%MAPPINGS%ELEMENTS)) THEN
              ELEMENTS_MAPPING=>DOMAIN%MAPPINGS%ELEMENTS
              IF(ASSOCIATED(DOMAIN%DECOMPOSITION)) THEN
                DECOMPOSITION=>DOMAIN%DECOMPOSITION
                IF(ASSOCIATED(DOMAIN%MESH)) THEN
                  MESH=>DOMAIN%MESH
                  component_idx=DOMAIN%meshComponentNumber
                  MESH_TOPOLOGY=>MESH%TOPOLOGY(component_idx)%PTR

                  CALL WorkGroup_NumberOfGroupNodesGet(NODES_MAPPING%workGroup,numberOfGroupComputationNodes, &
                    & err,error,*999)
                  CALL WorkGroup_GroupNodeNumberGet(NODES_MAPPING%workGroup,myGroupComputationNodeNumber, &
                    & err,error,*999)
                   
                  !Calculate the local and global numbers and set up the mappings
                  ALLOCATE(NODES_MAPPING%globalToLocalMap(MESH_TOPOLOGY%NODES%numberOfNodes),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate node mapping global to local map.",ERR,ERROR,*999)
                  NODES_MAPPING%numberOfGlobal=MESH_TOPOLOGY%NODES%numberOfNodes
                  ALLOCATE(LOCAL_NODE_NUMBERS(0:decomposition%numberOfDomains-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate local node numbers.",ERR,ERROR,*999)
                  LOCAL_NODE_NUMBERS=0
                  ALLOCATE(DOFS_MAPPING%globalToLocalMap(MESH_TOPOLOGY%dofs%numberOfDofs),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate dofs mapping global to local map.",ERR,ERROR,*999)
                  DOFS_MAPPING%numberOfGlobal=MESH_TOPOLOGY%DOFS%numberOfDofs
                  ALLOCATE(LOCAL_DOF_NUMBERS(0:decomposition%numberOfDomains-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate local dof numbers.",ERR,ERROR,*999)
                  LOCAL_DOF_NUMBERS=0
                  ALLOCATE(GHOST_NODES_LIST(0:decomposition%numberOfDomains-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate ghost nodes list.",ERR,ERROR,*999)
                  DO domain_idx=0,decomposition%numberOfDomains-1
                    NULLIFY(GHOST_NODES_LIST(domain_idx)%PTR)
                    CALL LIST_CREATE_START(GHOST_NODES_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(GHOST_NODES_LIST(domain_idx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_INITIAL_SIZE_SET(GHOST_NODES_LIST(domain_idx)%PTR,INT(MESH_TOPOLOGY%NODES%numberOfNodes/2), &
                      & ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(GHOST_NODES_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                  ENDDO !domain_idx
                  ALLOCATE(NUMBER_INTERNAL_NODES(0:decomposition%numberOfDomains-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate number of internal nodes.",ERR,ERROR,*999)
                  NUMBER_INTERNAL_NODES=0
                  ALLOCATE(NUMBER_BOUNDARY_NODES(0:decomposition%numberOfDomains-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate number of boundary nodes.",ERR,ERROR,*999)
                  NUMBER_BOUNDARY_NODES=0

                  !For the first pass just determine the internal and boundary nodes
                  DO node_idx=1,MESH_TOPOLOGY%NODES%numberOfNodes
                    NULLIFY(adjacentDomainsList)
                    CALL LIST_CREATE_START(adjacentDomainsList,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(adjacentDomainsList,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_INITIAL_SIZE_SET(adjacentDomainsList,decomposition%numberOfDomains,ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(adjacentDomainsList,ERR,ERROR,*999)
                    NULLIFY(ALL_ADJACENT_DOMAINS_LIST)
                    CALL LIST_CREATE_START(ALL_ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(ALL_ADJACENT_DOMAINS_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_INITIAL_SIZE_SET(ALL_ADJACENT_DOMAINS_LIST,decomposition%numberOfDomains,ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(ALL_ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    DO no_adjacent_element=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%numberOfSurroundingElements
                      adjacent_element=MESH_TOPOLOGY%NODES%NODES(node_idx)%surroundingElements(no_adjacent_element)
                      domain_no=DECOMPOSITION%elementDomain(adjacent_element)
                      CALL LIST_ITEM_ADD(adjacentDomainsList,domain_no,ERR,ERROR,*999)
                      DO domain_idx=1,ELEMENTS_MAPPING%globalToLocalMap(adjacent_element)%numberOfDomains
                        CALL LIST_ITEM_ADD(ALL_ADJACENT_DOMAINS_LIST,ELEMENTS_MAPPING%globalToLocalMap(adjacent_element)% &
                          & domainNumber(domain_idx),ERR,ERROR,*999)
                      ENDDO !domain_idx
                    ENDDO !no_adjacent_element
                    CALL LIST_REMOVE_DUPLICATES(adjacentDomainsList,ERR,ERROR,*999)
                    CALL LIST_DETACH_AND_DESTROY(adjacentDomainsList,numberOfDomains,DOMAINS,ERR,ERROR,*999)
                    CALL LIST_REMOVE_DUPLICATES(ALL_ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    CALL LIST_DETACH_AND_DESTROY(ALL_ADJACENT_DOMAINS_LIST,MAX_NUMBER_DOMAINS,ALL_DOMAINS,ERR,ERROR,*999)
                    IF(numberOfDomains/=MAX_NUMBER_DOMAINS) THEN !Ghost node
                      DO domain_idx=1,MAX_NUMBER_DOMAINS
                        domain_no=ALL_DOMAINS(domain_idx)
                        BOUNDARY_DOMAIN=.FALSE.
                        DO domain_idx2=1,numberOfDomains
                          IF(domain_no==DOMAINS(domain_idx2)) THEN
                            BOUNDARY_DOMAIN=.TRUE.
                            EXIT
                          ENDIF
                        ENDDO !domain_idx2
                        IF(.NOT.BOUNDARY_DOMAIN) CALL LIST_ITEM_ADD(GHOST_NODES_LIST(domain_no)%PTR,node_idx,ERR,ERROR,*999)
                      ENDDO !domain_idx
                    ENDIF
                    ALLOCATE(NODES_MAPPING%globalToLocalMap(node_idx)%localNumber(MAX_NUMBER_DOMAINS),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate node global to local map local number.",ERR,ERROR,*999)
                    ALLOCATE(NODES_MAPPING%globalToLocalMap(node_idx)%domainNumber(MAX_NUMBER_DOMAINS),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate node global to local map domain number.",ERR,ERROR,*999)
                    ALLOCATE(NODES_MAPPING%globalToLocalMap(node_idx)%localType(MAX_NUMBER_DOMAINS),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate node global to local map local type.",ERR,ERROR,*999)
                    DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%numberOfDerivatives
                      DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%numberOfVersions
                        ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                        ALLOCATE(DOFS_MAPPING%globalToLocalMap(ny)%localNumber(MAX_NUMBER_DOMAINS),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate dof global to local map local number.",ERR,ERROR,*999)
                        ALLOCATE(DOFS_MAPPING%globalToLocalMap(ny)%domainNumber(MAX_NUMBER_DOMAINS),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate dof global to local map domain number.",ERR,ERROR,*999)
                        ALLOCATE(DOFS_MAPPING%globalToLocalMap(ny)%localType(MAX_NUMBER_DOMAINS),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate dof global to local map local type.",ERR,ERROR,*999)
                      ENDDO !version_idx
                    ENDDO !derivative_idx
                    IF(numberOfDomains==1) THEN
                      !Node is an internal node
                      domain_no=DOMAINS(1)
                      NUMBER_INTERNAL_NODES(domain_no)=NUMBER_INTERNAL_NODES(domain_no)+1
                      !LOCAL_NODE_NUMBERS(domain_no)=LOCAL_NODE_NUMBERS(domain_no)+1
                      NODES_MAPPING%globalToLocalMap(node_idx)%numberOfDomains=1
                      !NODES_MAPPING%globalToLocalMap(node_idx)%localNumber(1)=LOCAL_NODE_NUMBERS(DOMAINS(1))
                      NODES_MAPPING%globalToLocalMap(node_idx)%localNumber(1)=-1
                      NODES_MAPPING%globalToLocalMap(node_idx)%domainNumber(1)=DOMAINS(1) 
                      NODES_MAPPING%globalToLocalMap(node_idx)%localType(1)=DOMAIN_LOCAL_INTERNAL
                      DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%numberOfDerivatives
                        DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%numberOfVersions
                          ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                          DOFS_MAPPING%globalToLocalMap(ny)%numberOfDomains=1
                          DOFS_MAPPING%globalToLocalMap(ny)%localNumber(1)=-1
                          DOFS_MAPPING%globalToLocalMap(ny)%domainNumber(1)=domain_no
                          DOFS_MAPPING%globalToLocalMap(ny)%localType(1)=DOMAIN_LOCAL_INTERNAL
                        ENDDO !version_idx
                      ENDDO !derivative_idx
                    ELSE
                      !Node is on the boundary of computation domains
                      NODES_MAPPING%globalToLocalMap(node_idx)%numberOfDomains=numberOfDomains
                      DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%numberOfDerivatives
                        DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%numberOfVersions
                          ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                          DOFS_MAPPING%globalToLocalMap(ny)%numberOfDomains=numberOfDomains
                        ENDDO !version_idx
                      ENDDO !derivative_idx
                      DO domain_idx=1,numberOfDomains
                        domain_no=DOMAINS(domain_idx)
                        !LOCAL_NODE_NUMBERS(domain_no)=LOCAL_NODE_NUMBERS(domain_no)+1
                        !NODES_MAPPING%globalToLocalMap(node_idx)%localNumber(domain_idx)=LOCAL_NODE_NUMBERS(domain_no)
                        NODES_MAPPING%globalToLocalMap(node_idx)%localNumber(domain_idx)=-1
                        NODES_MAPPING%globalToLocalMap(node_idx)%domainNumber(domain_idx)=domain_no
                        NODES_MAPPING%globalToLocalMap(node_idx)%localType(domain_idx)=DOMAIN_LOCAL_BOUNDARY
                        DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%numberOfDerivatives
                          DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%numberOfVersions
                            ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                            DOFS_MAPPING%globalToLocalMap(ny)%localNumber(domain_idx)=-1
                            DOFS_MAPPING%globalToLocalMap(ny)%domainNumber(domain_idx)=domain_no
                            DOFS_MAPPING%globalToLocalMap(ny)%localType(domain_idx)=DOMAIN_LOCAL_BOUNDARY
                          ENDDO !version_idx
                        ENDDO !derivative_idx
                      ENDDO !domain_idx
                    ENDIF
                    DEALLOCATE(DOMAINS) 
                    DEALLOCATE(ALL_DOMAINS)
                  ENDDO !node_idx

                  !For the second pass assign boundary nodes to one domain on the boundary and set local node numbers.
                  NUMBER_OF_NODES_PER_DOMAIN=FLOOR(REAL(MESH_TOPOLOGY%NODES%numberOfNodes,DP)/ &
                    & REAL(decomposition%numberOfDomains,DP))
                  ALLOCATE(DOMAIN%nodeDomain(MESH_TOPOLOGY%NODES%numberOfNodes),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate node domain",ERR,ERROR,*999)
                  DOMAIN%nodeDomain=-1
                  DO node_idx=1,MESH_TOPOLOGY%NODES%numberOfNodes
                    IF(NODES_MAPPING%globalToLocalMap(node_idx)%numberOfDomains==1) THEN !Internal node
                      domain_no=NODES_MAPPING%globalToLocalMap(node_idx)%domainNumber(1)
                      DOMAIN%nodeDomain(node_idx)=domain_no
                      LOCAL_NODE_NUMBERS(domain_no)=LOCAL_NODE_NUMBERS(domain_no)+1
                      NODES_MAPPING%globalToLocalMap(node_idx)%localNumber(1)=LOCAL_NODE_NUMBERS(domain_no)
                      DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%numberOfDerivatives
                        DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%numberOfVersions
                          ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                          LOCAL_DOF_NUMBERS(domain_no)=LOCAL_DOF_NUMBERS(domain_no)+1
                          DOFS_MAPPING%globalToLocalMap(ny)%localNumber(1)=LOCAL_DOF_NUMBERS(domain_no)
                        ENDDO !version_idx
                      ENDDO !derivative_idx
                    ELSE !Boundary node
                      numberOfDomains=NODES_MAPPING%globalToLocalMap(node_idx)%numberOfDomains
                      DO domain_idx=1,numberOfDomains
                        domain_no=NODES_MAPPING%globalToLocalMap(node_idx)%domainNumber(domain_idx)
                        IF(DOMAIN%nodeDomain(node_idx)<0) THEN
                          IF((NUMBER_INTERNAL_NODES(domain_no)+NUMBER_BOUNDARY_NODES(domain_no)<NUMBER_OF_NODES_PER_DOMAIN).OR. &
                            & (domain_idx==NODES_MAPPING%globalToLocalMap(node_idx)%numberOfDomains)) THEN
                            !Allocate the node to this domain
                            DOMAIN%nodeDomain(node_idx)=domain_no
                            NUMBER_BOUNDARY_NODES(domain_no)=NUMBER_BOUNDARY_NODES(domain_no)+1
                            LOCAL_NODE_NUMBERS(domain_no)=LOCAL_NODE_NUMBERS(domain_no)+1
                            !Reset the boundary information to be in the first domain index. The remaining domain indicies will
                            !be overwritten when the ghost nodes are calculated below. 
                            NODES_MAPPING%globalToLocalMap(node_idx)%numberOfDomains=1
                            NODES_MAPPING%globalToLocalMap(node_idx)%localNumber(1)=LOCAL_NODE_NUMBERS(domain_no)
                            NODES_MAPPING%globalToLocalMap(node_idx)%domainNumber(1)=domain_no
                            NODES_MAPPING%globalToLocalMap(node_idx)%localType(1)=DOMAIN_LOCAL_BOUNDARY
                            DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%numberOfDerivatives
                              DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%numberOfVersions
                                ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                                LOCAL_DOF_NUMBERS(domain_no)=LOCAL_DOF_NUMBERS(domain_no)+1
                                DOFS_MAPPING%globalToLocalMap(ny)%numberOfDomains=1
                                DOFS_MAPPING%globalToLocalMap(ny)%localNumber(1)=LOCAL_DOF_NUMBERS(domain_no)
                                DOFS_MAPPING%globalToLocalMap(ny)%domainNumber(1)=domain_no
                                DOFS_MAPPING%globalToLocalMap(ny)%localType(1)=DOMAIN_LOCAL_BOUNDARY
                              ENDDO !version_idx
                            ENDDO !derivative_idx
                          ELSE
                            !The node as already been assigned to a domain so it must be a ghost node in this domain
                            CALL LIST_ITEM_ADD(GHOST_NODES_LIST(domain_no)%PTR,node_idx,ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          !The node as already been assigned to a domain so it must be a ghost node in this domain
                          CALL LIST_ITEM_ADD(GHOST_NODES_LIST(domain_no)%PTR,node_idx,ERR,ERROR,*999)
                        ENDIF
                      ENDDO !domain_idx
                    ENDIF
                  ENDDO !node_idx
                  DEALLOCATE(NUMBER_INTERNAL_NODES)
                  
                  !Calculate ghost node and dof mappings
                  DO domain_idx=0,decomposition%numberOfDomains-1
                    CALL LIST_REMOVE_DUPLICATES(GHOST_NODES_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                    CALL LIST_DETACH_AND_DESTROY(GHOST_NODES_LIST(domain_idx)%PTR,NUMBER_OF_GHOST_NODES,GHOST_NODES,ERR,ERROR,*999)
                    DO no_ghost_node=1,NUMBER_OF_GHOST_NODES
                      ghost_node=GHOST_NODES(no_ghost_node)
                      LOCAL_NODE_NUMBERS(domain_idx)=LOCAL_NODE_NUMBERS(domain_idx)+1
                      NODES_MAPPING%globalToLocalMap(ghost_node)%numberOfDomains= &
                        & NODES_MAPPING%globalToLocalMap(ghost_node)%numberOfDomains+1
                      NODES_MAPPING%globalToLocalMap(ghost_node)%localNumber( &
                        & NODES_MAPPING%globalToLocalMap(ghost_node)%numberOfDomains)= &
                        & LOCAL_NODE_NUMBERS(domain_idx)
                      NODES_MAPPING%globalToLocalMap(ghost_node)%domainNumber( &
                        & NODES_MAPPING%globalToLocalMap(ghost_node)%numberOfDomains)=domain_idx
                      NODES_MAPPING%globalToLocalMap(ghost_node)%localType( &
                        & NODES_MAPPING%globalToLocalMap(ghost_node)%numberOfDomains)= &
                        & DOMAIN_LOCAL_GHOST
                      DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(ghost_node)%numberOfDerivatives
                        DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(ghost_node)%DERIVATIVES(derivative_idx)%numberOfVersions
                          ny=MESH_TOPOLOGY%NODES%NODES(ghost_node)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                          LOCAL_DOF_NUMBERS(domain_idx)=LOCAL_DOF_NUMBERS(domain_idx)+1
                          DOFS_MAPPING%globalToLocalMap(ny)%numberOfDomains= &
                            & DOFS_MAPPING%globalToLocalMap(ny)%numberOfDomains+1
                          DOFS_MAPPING%globalToLocalMap(ny)%localNumber( &
                            & DOFS_MAPPING%globalToLocalMap(ny)%numberOfDomains)= &
                            & LOCAL_DOF_NUMBERS(domain_idx)
                          DOFS_MAPPING%globalToLocalMap(ny)%domainNumber( &
                            & DOFS_MAPPING%globalToLocalMap(ny)%numberOfDomains)=domain_idx
                          DOFS_MAPPING%globalToLocalMap(ny)%localType( &
                            & DOFS_MAPPING%globalToLocalMap(ny)%numberOfDomains)= &
                            & DOMAIN_LOCAL_GHOST
                        ENDDO !version_idx
                      ENDDO !derivative_idx
                    ENDDO !no_ghost_node
                    DEALLOCATE(GHOST_NODES)
                  ENDDO !domain_idx
                  
                  !Check decomposition and check that each domain has a node in it.
                  ALLOCATE(NODE_COUNT(0:numberOfGroupComputationNodes-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate node count.",ERR,ERROR,*999)
                  NODE_COUNT=0
                  DO node_idx=1,MESH_TOPOLOGY%NODES%numberOfNodes
                    no_computation_node=DOMAIN%nodeDomain(node_idx)
                    IF(no_computation_node>=0.AND.no_computation_node<numberOfGroupComputationNodes) THEN
                      NODE_COUNT(no_computation_node)=NODE_COUNT(no_computation_node)+1
                    ELSE
                      LOCAL_ERROR="The computation node number of "// &
                        & TRIM(NUMBER_TO_VSTRING(no_computation_node,"*",ERR,ERROR))// &
                        & " for node number "//TRIM(NUMBER_TO_VSTRING(node_idx,"*",ERR,ERROR))// &
                        & " is invalid. The computation node number must be between 0 and "// &
                        & TRIM(NUMBER_TO_VSTRING(numberOfGroupComputationNodes-1,"*",ERR,ERROR))//"."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDDO !node_idx
                  DO no_computation_node=0,numberOfGroupComputationNodes-1
                    IF(NODE_COUNT(no_computation_node)==0) THEN
                      LOCAL_ERROR="Invalid decomposition. There are no nodes in computation node "// &
                        & TRIM(NUMBER_TO_VSTRING(no_computation_node,"*",ERR,ERROR))//"."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDDO !no_computation_node
                  DEALLOCATE(NODE_COUNT)
          
                  DEALLOCATE(GHOST_NODES_LIST)
                  DEALLOCATE(LOCAL_NODE_NUMBERS)
                  
                  !Calculate node and dof local to global maps from global to local map
                  CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(NODES_MAPPING,ERR,ERROR,*999)
                  CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(DOFS_MAPPING,ERR,ERROR,*999)
                  
                ELSE
                  CALL FlagError("Domain mesh is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Domain decomposition is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Domain mappings elements is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Domain mappings dofs is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Domain mappings nodes is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Domain mappings is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain is not associated.",ERR,ERROR,*998)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Node decomposition :",ERR,ERROR,*999)
      DO node_idx=1,MESH_TOPOLOGY%NODES%numberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Node = ",node_idx,ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Domain = ",DOMAIN%nodeDomain(node_idx),ERR,ERROR,*999)
      ENDDO !node_idx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Node mappings :",ERR,ERROR,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Global to local map :",ERR,ERROR,*999)
      DO node_idx=1,MESH_TOPOLOGY%NODES%numberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global node = ",node_idx,ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of domains  = ", &
          & NODES_MAPPING%globalToLocalMap(node_idx)%numberOfDomains,ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%globalToLocalMap(node_idx)% &
          & numberOfDomains,8,8,NODES_MAPPING%globalToLocalMap(node_idx)%localNumber, &
          & '("      Local number :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%globalToLocalMap(node_idx)% &
            & numberOfDomains,8,8,NODES_MAPPING%globalToLocalMap(node_idx)%domainNumber, &
            & '("      Domain number:",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%globalToLocalMap(node_idx)% &
          & numberOfDomains,8,8,NODES_MAPPING%globalToLocalMap(node_idx)%localType, &
          & '("      Local type   :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
      ENDDO !node_idx     
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map :",ERR,ERROR,*999)
      DO node_idx=1,NODES_MAPPING%totalNumberOfLocal
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Local node = ",node_idx,ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Global node = ", &
          & NODES_MAPPING%localToGlobalMap(node_idx),ERR,ERROR,*999)
      ENDDO !node_idx
      IF(DIAGNOSTICS2) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Internal nodes :",ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of internal nodes = ", &
          & NODES_MAPPING%numberOfInternal,ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%numberOfInternal,8,8, &
          & NODES_MAPPING%domainList(NODES_MAPPING%internalStart:NODES_MAPPING%internalFinish), &
          & '("    Internal nodes:",8(X,I7))','(19X,8(X,I7))',ERR,ERROR,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary nodes :",ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary nodes = ", &
          & NODES_MAPPING%numberOfBoundary,ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%numberOfBoundary,8,8, &
          & NODES_MAPPING%domainList(NODES_MAPPING%boundaryStart:NODES_MAPPING%boundaryFinish), &
          & '("    Boundary nodes:",8(X,I7))','(19X,8(X,I7))',ERR,ERROR,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Ghost nodes :",ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost nodes = ", &
          & NODES_MAPPING%numberOfGhost,ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%numberOfGhost,8,8, &
          & NODES_MAPPING%domainList(NODES_MAPPING%ghostStart:NODES_MAPPING%ghostFinish), &
          & '("    Ghost nodes   :",8(X,I7))','(19X,8(X,I7))',ERR,ERROR,*999)
      ENDIF
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains :",ERR,ERROR,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
        & NODES_MAPPING%numberOfAdjacentDomains,ERR,ERROR,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%numberOfDomains+1,8,8, &
        & NODES_MAPPING%adjacentDomainsPtr,'("    Adjacent domains ptr  :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%adjacentDomainsPtr( &
        & NODES_MAPPING%numberOfDomains)-1,8,8,NODES_MAPPING%adjacentDomainsList, &
        '("    Adjacent domains list :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
      DO domain_idx=1,NODES_MAPPING%numberOfAdjacentDomains
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",domain_idx,ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
          & NODES_MAPPING%adjacentDomains(domain_idx)%domainNumber,ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
          & NODES_MAPPING%adjacentDomains(domain_idx)%numberOfSendGhosts,ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%adjacentDomains(domain_idx)% &
          & numberOfSendGhosts,6,6,NODES_MAPPING%adjacentDomains(domain_idx)%localGhostSendIndices, &
          & '("      Local send ghost indicies       :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)      
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of recieve ghosts = ", &
          & NODES_MAPPING%adjacentDomains(domain_idx)%numberOfReceiveGhosts,ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%adjacentDomains(domain_idx)% &
          & numberOfReceiveGhosts,6,6,NODES_MAPPING%adjacentDomains(domain_idx)%localGhostReceiveIndices, &
          & '("      Local receive ghost indicies    :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)              
      ENDDO !domain_idx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Dofs mappings :",ERR,ERROR,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Global to local map :",ERR,ERROR,*999)
      DO ny=1,MESH_TOPOLOGY%DOFS%numberOfDofs
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global dof = ",ny,ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of domains  = ", &
          & DOFS_MAPPING%globalToLocalMap(ny)%numberOfDomains,ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%globalToLocalMap(ny)% &
          & numberOfDomains,8,8,DOFS_MAPPING%globalToLocalMap(ny)%localNumber, &
          & '("      Local number :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%globalToLocalMap(ny)% &
            & numberOfDomains,8,8,DOFS_MAPPING%globalToLocalMap(ny)%domainNumber, &
            & '("      Domain number:",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%globalToLocalMap(ny)% &
          & numberOfDomains,8,8,DOFS_MAPPING%globalToLocalMap(ny)%localType, &
          & '("      Local type   :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
      ENDDO !ny   
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map :",ERR,ERROR,*999)
      DO ny=1,DOFS_MAPPING%totalNumberOfLocal
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Local dof = ",ny,ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Global dof = ", &
          & DOFS_MAPPING%localToGlobalMap(ny),ERR,ERROR,*999)
      ENDDO !node_idx
      IF(DIAGNOSTICS2) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Internal dofs :",ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of internal dofs = ", &
          & DOFS_MAPPING%numberOfInternal,ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%numberOfInternal,8,8, &
          & DOFS_MAPPING%domainList(DOFS_MAPPING%internalStart:DOFS_MAPPING%internalFinish), &
          & '("    Internal dofs:",8(X,I7))','(18X,8(X,I7))',ERR,ERROR,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary dofs :",ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary dofs = ", &
          & DOFS_MAPPING%numberOfBoundary,ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%numberOfBoundary,8,8, &
          & DOFS_MAPPING%domainList(DOFS_MAPPING%boundaryStart:DOFS_MAPPING%boundaryFinish), &
          & '("    Boundary dofs:",8(X,I7))','(18X,8(X,I7))',ERR,ERROR,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Ghost dofs :",ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost dofs = ", &
          & DOFS_MAPPING%numberOfGhost,ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%numberOfGhost,8,8, &
          & DOFS_MAPPING%domainList(DOFS_MAPPING%ghostStart:DOFS_MAPPING%ghostFinish), &
          & '("    Ghost dofs   :",8(X,I7))','(18X,8(X,I7))',ERR,ERROR,*999)
      ENDIF
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains :",ERR,ERROR,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
        & DOFS_MAPPING%numberOfAdjacentDomains,ERR,ERROR,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%numberOfDomains+1,8,8, &
        & DOFS_MAPPING%adjacentDomainsPtr,'("    Adjacent domains ptr  :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%adjacentDomainsPtr( &
        & DOFS_MAPPING%numberOfDomains)-1,8,8,DOFS_MAPPING%adjacentDomainsList, &
        '("    Adjacent domains list :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
      DO domain_idx=1,DOFS_MAPPING%numberOfAdjacentDomains
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",domain_idx,ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
          & DOFS_MAPPING%adjacentDomains(domain_idx)%domainNumber,ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
          & DOFS_MAPPING%adjacentDomains(domain_idx)%numberOfSendGhosts,ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%adjacentDomains(domain_idx)% &
          & numberOfSendGhosts,6,6,DOFS_MAPPING%adjacentDomains(domain_idx)%localGhostSendIndices, &
          & '("      Local send ghost indicies       :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)      
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of recieve ghosts = ", &
          & DOFS_MAPPING%adjacentDomains(domain_idx)%numberOfReceiveGhosts,ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%adjacentDomains(domain_idx)% &
          & numberOfReceiveGhosts,6,6,DOFS_MAPPING%adjacentDomains(domain_idx)%localGhostReceiveIndices, &
          & '("      Local receive ghost indicies    :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)              
      ENDDO !domain_idx
    ENDIF
    
    EXITS("DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE")
    RETURN
999 IF(ALLOCATED(DOMAINS)) DEALLOCATE(DOMAINS)
    IF(ALLOCATED(ALL_DOMAINS)) DEALLOCATE(ALL_DOMAINS)
    IF(ALLOCATED(GHOST_NODES)) DEALLOCATE(GHOST_NODES)
    IF(ALLOCATED(NUMBER_INTERNAL_NODES)) DEALLOCATE(NUMBER_INTERNAL_NODES)
    IF(ALLOCATED(NUMBER_BOUNDARY_NODES)) DEALLOCATE(NUMBER_BOUNDARY_NODES)
    IF(ASSOCIATED(DOMAIN%MAPPINGS%NODES)) CALL DOMAIN_MAPPINGS_NODES_FINALISE(DOMAIN%MAPPINGS,DUMMY_ERR,DUMMY_ERROR,*998)
998 IF(ASSOCIATED(DOMAIN%MAPPINGS%DOFS)) CALL DOMAIN_MAPPINGS_DOFS_FINALISE(DOMAIN%MAPPINGS,DUMMY_ERR,DUMMY_ERROR,*997)
997 ERRORSEXITS("DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Finalises the node mapping in the given domain mappings. \todo pass in the nodes mapping
  SUBROUTINE DOMAIN_MAPPINGS_NODES_FINALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainMappingsType), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mapping to finalise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("DOMAIN_MAPPINGS_NODES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%NODES)) THEN
        CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(DOMAIN_MAPPINGS%NODES,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_NODES_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_NODES_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_NODES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the node mapping in the given domain mapping. \todo finalise on error
  SUBROUTINE DOMAIN_MAPPINGS_NODES_INITIALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainMappingsType), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mappings to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("DOMAIN_MAPPINGS_NODES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%NODES)) THEN
        CALL FlagError("Domain nodes mappings are already associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(DOMAIN_MAPPINGS%NODES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate domain mappings nodes.",ERR,ERROR,*999)
        CALL DomainMappings_MappingInitialise(DOMAIN_MAPPINGS%DOMAIN%DECOMPOSITION%workGroup, &
          & DOMAIN_MAPPINGS%NODES,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_NODES_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_NODES_INITIALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_NODES_INITIALISE

  !
  !================================================================================================================================
  !

  !>Calculates the domain topology.
  SUBROUTINE DOMAIN_TOPOLOGY_CALCULATE(TOPOLOGY,ERR,ERROR,*)

   !Argument variables
    TYPE(DomainTopologyType), POINTER :: TOPOLOGY !<A pointer to the domain topology to calculate.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne,np
    TYPE(BasisType), POINTER :: BASIS

    ENTERS("DOMAIN_TOPOLOGY_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      !Find maximum number of element parameters for all elements in the domain toplogy
      TOPOLOGY%ELEMENTS%maximumNumberOfElementParameters=-1
      DO ne=1,TOPOLOGY%ELEMENTS%totalNumberOfElements
        BASIS=>TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS
        IF(ASSOCIATED(BASIS)) THEN
          IF(BASIS%numberOfElementParameters>TOPOLOGY%ELEMENTS%maximumNumberOfElementParameters) &
            & TOPOLOGY%ELEMENTS%maximumNumberOfElementParameters=BASIS%numberOfElementParameters
        ELSE
          CALL FlagError("Basis is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDDO !ne
      !Find maximum number of derivatives for all nodes in the domain toplogy
      TOPOLOGY%NODES%maximumNumberOfDerivatives=-1
      DO np=1,TOPOLOGY%NODES%totalNumberOfNodes
        IF(TOPOLOGY%NODES%NODES(np)%numberOfDerivatives>TOPOLOGY%NODES%maximumNumberOfDerivatives) &
            & TOPOLOGY%NODES%maximumNumberOfDerivatives=TOPOLOGY%NODES%NODES(np)%numberOfDerivatives
      ENDDO !np
      !Calculate the elements surrounding the nodes in the domain topology
      CALL DomainTopology_NodesSurroundingElementsCalculate(TOPOLOGY,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_CALCULATE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Initialises the local domain topology from the mesh topology.
  SUBROUTINE DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainType), POINTER :: DOMAIN !<A pointer to the domain to initialise the domain topology from the mesh topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: local_element,global_element,local_node,global_node,version_idx,derivative_idx,node_idx,dof_idx, &
      & component_idx
    INTEGER(INTG) :: ne,nn,nkk,INSERT_STATUS
    LOGICAL :: FOUND
    TYPE(BasisType), POINTER :: BASIS
    TYPE(MeshType), POINTER :: MESH
    TYPE(DomainElementsType), POINTER :: DOMAIN_ELEMENTS
    TYPE(MeshElementsType), POINTER :: MESH_ELEMENTS
    TYPE(DomainNodesType), POINTER :: DOMAIN_NODES
    TYPE(MeshNodesType), POINTER :: MESH_NODES
    TYPE(DomainDOFsType), POINTER :: DOMAIN_DOFS

    ENTERS("DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH",ERR,ERROR,*999)
    
    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
        IF(ASSOCIATED(DOMAIN%MAPPINGS)) THEN
          IF(ASSOCIATED(DOMAIN%MESH)) THEN
            MESH=>DOMAIN%MESH
            component_idx=DOMAIN%meshComponentNumber
            IF(ASSOCIATED(MESH%TOPOLOGY(component_idx)%PTR)) THEN
              MESH_ELEMENTS=>MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS
              DOMAIN_ELEMENTS=>DOMAIN%TOPOLOGY%ELEMENTS
              MESH_NODES=>MESH%TOPOLOGY(component_idx)%PTR%NODES
              DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
              DOMAIN_DOFS=>DOMAIN%TOPOLOGY%DOFS
              ALLOCATE(DOMAIN_ELEMENTS%ELEMENTS(DOMAIN%MAPPINGS%ELEMENTS%totalNumberOfLocal),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate domain elements elements.",ERR,ERROR,*999)
              DOMAIN_ELEMENTS%numberOfElements=DOMAIN%MAPPINGS%ELEMENTS%numberOfLocal
              DOMAIN_ELEMENTS%totalNumberOfElements=DOMAIN%MAPPINGS%ELEMENTS%totalNumberOfLocal
              DOMAIN_ELEMENTS%numberOfGlobalElements=DOMAIN%MAPPINGS%ELEMENTS%numberOfGlobal
              ALLOCATE(DOMAIN_NODES%NODES(DOMAIN%MAPPINGS%NODES%totalNumberOfLocal),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate domain nodes nodes.",ERR,ERROR,*999)
              DOMAIN_NODES%numberOfNodes=DOMAIN%MAPPINGS%NODES%numberOfLocal
              DOMAIN_NODES%totalNumberOfNodes=DOMAIN%MAPPINGS%NODES%totalNumberOfLocal
              DOMAIN_NODES%numberOfGlobalNodes=DOMAIN%MAPPINGS%NODES%numberOfGlobal
              ALLOCATE(DOMAIN_DOFS%dofIndex(3,DOMAIN%MAPPINGS%DOFS%totalNumberOfLocal),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate domain dofs dof index.",ERR,ERROR,*999)
              DOMAIN_DOFS%numberOfDofs=DOMAIN%MAPPINGS%DOFS%numberOfLocal
              DOMAIN_DOFS%totalNumberOfDofs=DOMAIN%MAPPINGS%DOFS%totalNumberOfLocal
              DOMAIN_DOFS%numberOfGlobalDofs=DOMAIN%MAPPINGS%DOFS%numberOfGlobal
              !Loop over the domain nodes and calculate the parameters from the mesh nodes
              CALL TREE_CREATE_START(DOMAIN_NODES%nodesTree,ERR,ERROR,*999)
              CALL TREE_INSERT_TYPE_SET(DOMAIN_NODES%nodesTree,TREE_NO_DUPLICATES_ALLOWED,ERR,ERROR,*999)
              CALL TREE_CREATE_FINISH(DOMAIN_NODES%nodesTree,ERR,ERROR,*999)
              dof_idx=0
              DO local_node=1,DOMAIN_NODES%totalNumberOfNodes
                CALL DOMAIN_TOPOLOGY_NODE_INITIALISE(DOMAIN_NODES%NODES(local_node),ERR,ERROR,*999)
                global_node=DOMAIN%MAPPINGS%NODES%localToGlobalMap(local_node)
                DOMAIN_NODES%NODES(local_node)%localNumber=local_node
                DOMAIN_NODES%NODES(local_node)%meshNumber=global_node
                DOMAIN_NODES%NODES(local_node)%globalNumber=MESH_NODES%NODES(global_node)%globalNumber
                DOMAIN_NODES%NODES(local_node)%userNumber=MESH_NODES%NODES(global_node)%userNumber
                CALL TREE_ITEM_INSERT(DOMAIN_NODES%nodesTree,DOMAIN_NODES%NODES(local_node)%userNumber,local_node, &
                  & INSERT_STATUS,ERR,ERROR,*999)
                DOMAIN_NODES%NODES(local_node)%numberOfSurroundingElements=0
                NULLIFY(DOMAIN_NODES%NODES(local_node)%surroundingElements)
                DOMAIN_NODES%NODES(local_node)%numberOfDerivatives=MESH_NODES%NODES(global_node)%numberOfDerivatives
                ALLOCATE(DOMAIN_NODES%NODES(local_node)%DERIVATIVES(MESH_NODES%NODES(global_node)%numberOfDerivatives),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate domain node derivatives.",ERR,ERROR,*999)
                DO derivative_idx=1,DOMAIN_NODES%NODES(local_node)%numberOfDerivatives
                  CALL DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE( &
                    & DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx),ERR,ERROR,*999)
                  DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%globalDerivativeIndex= & 
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%globalDerivativeIndex
                  DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%partialDerivativeIndex= &
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%partialDerivativeIndex
                  DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%numberOfVersions= & 
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%numberOfVersions
                  ALLOCATE(DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%userVersionNumbers( &
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%numberOfVersions),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate node derivative version numbers.",ERR,ERROR,*999)
                  DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%userVersionNumbers(1: &
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%numberOfVersions)= &
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%userVersionNumbers(1: &
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%numberOfVersions)
                  ALLOCATE(DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%dofIndex( &
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%numberOfVersions),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate node dervative versions dof index.",ERR,ERROR,*999)
                  DO version_idx=1,DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%numberOfVersions
                    dof_idx=dof_idx+1
                    DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)=dof_idx
                    DOMAIN_DOFS%dofIndex(1,dof_idx)=version_idx
                    DOMAIN_DOFS%dofIndex(2,dof_idx)=derivative_idx
                    DOMAIN_DOFS%dofIndex(3,dof_idx)=local_node
                  ENDDO !version_idx
                ENDDO !derivative_idx
                DOMAIN_NODES%NODES(local_node)%boundaryNode=MESH_NODES%NODES(global_node)%boundaryNode
              ENDDO !local_node
              !Loop over the domain elements and renumber from the mesh elements
              DO local_element=1,DOMAIN_ELEMENTS%totalNumberOfElements
                CALL DOMAIN_TOPOLOGY_ELEMENT_INITIALISE(DOMAIN_ELEMENTS%ELEMENTS(local_element),ERR,ERROR,*999)
                global_element=DOMAIN%MAPPINGS%ELEMENTS%localToGlobalMap(local_element)               
                BASIS=>MESH_ELEMENTS%ELEMENTS(global_element)%BASIS
                DOMAIN_ELEMENTS%ELEMENTS(local_element)%NUMBER=local_element
                DOMAIN_ELEMENTS%ELEMENTS(local_element)%BASIS=>BASIS
                ALLOCATE(DOMAIN_ELEMENTS%ELEMENTS(local_element)%elementNodes(BASIS%numberOfNodes),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate domain elements element nodes.",ERR,ERROR,*999)
                ALLOCATE(DOMAIN_ELEMENTS%ELEMENTS(local_element)%elementDerivatives(BASIS%maximumNumberOfDerivatives, &
                  & BASIS%numberOfNodes),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate domain elements element derivatives.",ERR,ERROR,*999)
                ALLOCATE(DOMAIN_ELEMENTS%ELEMENTS(local_element)%elementVersions(BASIS%maximumNumberOfDerivatives, &
                  & BASIS%numberOfNodes),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate domain elements element versions.",ERR,ERROR,*999)
                DO nn=1,BASIS%numberOfNodes
                  global_node=MESH_ELEMENTS%ELEMENTS(global_element)%meshElementNodes(nn)
                  local_node=DOMAIN%MAPPINGS%NODES%globalToLocalMap(global_node)%localNumber(1)
                  DOMAIN_ELEMENTS%ELEMENTS(local_element)%elementNodes(nn)=local_node
                  DO derivative_idx=1,BASIS%numberOfDerivatives(nn)
                    !Find equivalent node derivative by matching partial derivative index
                    !/todo Take a look at this - is it needed any more?
                    FOUND=.FALSE.
                    DO nkk=1,DOMAIN_NODES%NODES(local_node)%numberOfDerivatives
                      IF(DOMAIN_NODES%NODES(local_node)%DERIVATIVES(nkk)%partialDerivativeIndex == &
                        & BASIS%partialDerivativeIndex(derivative_idx,nn)) THEN
                        FOUND=.TRUE.
                        EXIT
                      ENDIF
                    ENDDO !nkk
                    IF(FOUND) THEN
                      DOMAIN_ELEMENTS%ELEMENTS(local_element)%elementDerivatives(derivative_idx,nn)=nkk
                      DOMAIN_ELEMENTS%ELEMENTS(local_element)%elementVersions(derivative_idx,nn) = & 
                        & MESH_ELEMENTS%ELEMENTS(global_element)%userElementNodeVersions(derivative_idx,nn)
                    ELSE
                      CALL FlagError("Could not find equivalent node derivative",ERR,ERROR,*999)
                    ENDIF
                  ENDDO !derivative_idx
                ENDDO !nn
              ENDDO !local_element                       
            ELSE
              CALL FlagError("Mesh topology is not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Mesh is not associated",ERR,ERROR,*999)

          ENDIF
        ELSE
          CALL FlagError("Domain mapping is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Domain topology is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Initialised domain topology :",ERR,ERROR,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of domain nodes = ",DOMAIN_NODES%totalNumberOfNodes, &
        & ERR,ERROR,*999)
      DO node_idx=1,DOMAIN_NODES%totalNumberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Node number = ",DOMAIN_NODES%NODES(node_idx)%localNumber, &
          & ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Node mesh number = ",DOMAIN_NODES%NODES(node_idx)%meshNumber, &
          & ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Node global number = ",DOMAIN_NODES%NODES(node_idx)%globalNumber, &
          & ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Node user number = ",DOMAIN_NODES%NODES(node_idx)%userNumber, &
          & ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of derivatives = ", &
          & DOMAIN_NODES%NODES(node_idx)%numberOfDerivatives,ERR,ERROR,*999)
        DO derivative_idx=1,DOMAIN_NODES%NODES(local_node)%numberOfDerivatives
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Node local derivative number = ",derivative_idx,ERR,ERROR,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Global derivative index = ", &
            & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%globalDerivativeIndex,ERR,ERROR,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Partial derivative index = ", &
            & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%partialDerivativeIndex,ERR,ERROR,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
            & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%numberOfVersions,4,4, &
            & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%dofIndex, &
            & '("        Degree-of-freedom index(version_idx)  :",4(X,I9))','(36X,4(X,I9))',ERR,ERROR,*999)
        ENDDO !derivative_idx
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Boundary node = ", &
          & DOMAIN_NODES%NODES(node_idx)%boundaryNode,ERR,ERROR,*999)
     ENDDO !node_idx
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Total Number of domain dofs = ",DOMAIN_DOFS%totalNumberOfDofs, &
        & ERR,ERROR,*999)
      DO dof_idx=1,DOMAIN_DOFS%totalNumberOfDofs
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Dof number = ",dof_idx,ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,3,3, &
          & DOMAIN_DOFS%dofIndex(:,dof_idx),'("    Degree-of-freedom index :",3(X,I9))','(29X,3(X,I9))', &
          & ERR,ERROR,*999)
      ENDDO !dof_idx
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of domain elements = ", &
        & DOMAIN_ELEMENTS%totalNumberOfElements,ERR,ERROR,*999)
      DO ne=1,DOMAIN_ELEMENTS%totalNumberOfElements
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Element number = ",DOMAIN_ELEMENTS%ELEMENTS(ne)%NUMBER,ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Basis user number = ", &
          & DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%userNumber,ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of local nodes = ", &
          & DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%numberOfNodes,ERR,ERROR,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%numberOfNodes,8,8, &
          & DOMAIN_ELEMENTS%ELEMENTS(ne)%elementNodes,'("    Element nodes(nn) :",8(X,I9))','(23X,8(X,I9))', &
          & ERR,ERROR,*999)
        DO nn=1,DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%numberOfNodes
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Local node number : ",nn,ERR,ERROR,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%numberOfDerivatives(nn),8,8, &
            & DOMAIN_ELEMENTS%ELEMENTS(ne)%elementDerivatives(:,nn), &
            & '("        Element derivatives :",8(X,I2))','(29X,8(X,I2))',ERR,ERROR,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%numberOfDerivatives(nn),8,8, &
            & DOMAIN_ELEMENTS%ELEMENTS(ne)%elementVersions(:,nn), &
            & '("        Element versions    :",8(X,I2))','(29X,8(X,I2))',ERR,ERROR,*999)
        ENDDO !nn
      ENDDO !ne
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH
  
  !
  !================================================================================================================================
  !

  !>Finalises the dofs in the given domain topology. \todo pass in the dofs topolgy
  SUBROUTINE DOMAIN_TOPOLOGY_DOFS_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: TOPOLOGY !<A pointer to the domain topology to finalise the dofs for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_DOFS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%DOFS)) THEN
        IF(ALLOCATED(TOPOLOGY%DOFS%dofIndex)) DEALLOCATE(TOPOLOGY%DOFS%dofIndex)
        DEALLOCATE(TOPOLOGY%DOFS)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_TOPOLOGY_DOFS_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_DOFS_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_DOFS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the dofs data structures for a domain topology. \todo finalise on exit
  SUBROUTINE DOMAIN_TOPOLOGY_DOFS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the dofs for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_DOFS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%DOFS)) THEN
        CALL FlagError("Decomposition already has topology dofs associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%DOFS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology dofs",ERR,ERROR,*999)
        TOPOLOGY%DOFS%numberOfDofs=0
        TOPOLOGY%DOFS%totalNumberOfDofs=0
        TOPOLOGY%DOFS%numberOfGlobalDofs=0
        TOPOLOGY%DOFS%DOMAIN=>TOPOLOGY%DOMAIN
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_DOFS_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_DOFS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_DOFS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the given domain topology element.
  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENT_FINALISE(ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainElementType) :: ELEMENT !<The domain element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(ELEMENT%elementNodes)) DEALLOCATE(ELEMENT%elementNodes)
    IF(ALLOCATED(ELEMENT%elementDerivatives)) DEALLOCATE(ELEMENT%elementDerivatives)
    IF(ALLOCATED(ELEMENT%elementVersions)) DEALLOCATE(ELEMENT%elementVersions)
 
    EXITS("DOMAIN_TOPOLOGY_ELEMENT_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given domain topology element.
  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENT_INITIALISE(ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainElementType) :: ELEMENT !<The domain element to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR,*999)

    ELEMENT%NUMBER=0
    NULLIFY(ELEMENT%BASIS)
  
    EXITS("DOMAIN_TOPOLOGY_ELEMENT_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the elements in the given domain topology. \todo pass in the domain elements
  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENTS_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: TOPOLOGY !<A pointer to the domain topology to finalise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne

    ENTERS("DOMAIN_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        DO ne=1,TOPOLOGY%ELEMENTS%totalNumberOfElements
          CALL DOMAIN_TOPOLOGY_ELEMENT_FINALISE(TOPOLOGY%ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
        ENDDO !ne
        IF(ASSOCIATED(TOPOLOGY%ELEMENTS%ELEMENTS)) DEALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS)
        DEALLOCATE(TOPOLOGY%ELEMENTS)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_TOPOLOGY_ELEMENTS_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_ELEMENTS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the element data structures for a domain topology. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        CALL FlagError("Domain already has topology elements associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%ELEMENTS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology elements",ERR,ERROR,*999)
        TOPOLOGY%ELEMENTS%numberOfElements=0
        TOPOLOGY%ELEMENTS%totalNumberOfElements=0
        TOPOLOGY%ELEMENTS%numberOfGlobalElements=0
        TOPOLOGY%ELEMENTS%DOMAIN=>TOPOLOGY%DOMAIN
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS)
        TOPOLOGY%ELEMENTS%maximumNumberOfElementParameters=0
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE
  
  
  !
  !================================================================================================================================
  !

  !>Finalises the topology in the given domain. \todo pass in domain topology
  SUBROUTINE DOMAIN_TOPOLOGY_FINALISE(DOMAIN,ERR,ERROR,*)

   !Argument variables
    TYPE(DomainType), POINTER :: DOMAIN !<A pointer to the domain to finalise the topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      CALL DOMAIN_TOPOLOGY_NODES_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      CALL DOMAIN_TOPOLOGY_DOFS_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      CALL DOMAIN_TOPOLOGY_ELEMENTS_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      CALL DOMAIN_TOPOLOGY_LINES_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      CALL DOMAIN_TOPOLOGY_FACES_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      DEALLOCATE(DOMAIN%TOPOLOGY)
    ELSE
      CALL FlagError("Domain is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_TOPOLOGY_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the topology for a given domain. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_INITIALISE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainType), POINTER :: DOMAIN !A pointer to the domain to initialise the topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
        CALL FlagError("Domain already has topology associated",ERR,ERROR,*999)
      ELSE
        !Allocate domain topology
        ALLOCATE(DOMAIN%TOPOLOGY,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Domain topology could not be allocated",ERR,ERROR,*999)
        DOMAIN%TOPOLOGY%DOMAIN=>DOMAIN
        NULLIFY(DOMAIN%TOPOLOGY%ELEMENTS)
        NULLIFY(DOMAIN%TOPOLOGY%NODES)
        NULLIFY(DOMAIN%TOPOLOGY%DOFS)
        NULLIFY(DOMAIN%TOPOLOGY%LINES)
        NULLIFY(DOMAIN%TOPOLOGY%FACES)
        !Initialise the topology components
        CALL DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        CALL DOMAIN_TOPOLOGY_NODES_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        CALL DOMAIN_TOPOLOGY_DOFS_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        CALL DOMAIN_TOPOLOGY_LINES_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        CALL DOMAIN_TOPOLOGY_FACES_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        !Initialise the domain topology from the domain mappings and the mesh it came from
        CALL DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH(DOMAIN,ERR,ERROR,*999)
        !Calculate the topological information.
        CALL DOMAIN_TOPOLOGY_CALCULATE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Get the basis for an element in the domain identified by its user number
  SUBROUTINE DomainTopology_ElementBasisGet(domainTopology,userElementNumber,basis,err,error,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: domainTopology !<A pointer to the domain topology to get the element basis for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The element user number to get the basis for
    TYPE(BasisType), POINTER, INTENT(OUT) :: basis !<On return, a pointer to the basis for the element.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    LOGICAL :: userElementExists,ghostElement
    INTEGER(INTG) :: localElementNumber
    TYPE(VARYING_STRING) :: localError

    ENTERS("DomainTopology_ElementBasisGet",err,error,*999)

    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain topology is not associated.",err,error,*999)
    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(domainTopology%elements)) CALL FlagError("Domain topology elements is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(domainTopology%elements%elements)) &
      & CALL FlagError("Domain topology elements elements is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(domainTopology%domain)) CALL FlagError("Domain topology domain is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(domainTopology%domain%decomposition)) &
      & CALL FlagError("Domain topology domain decomposition is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(domainTopology%domain%decomposition%topology)) &
      & CALL FlagError("Domain topology domain decomposition topology is not associated.",err,error,*999)
    
    CALL DecompositionTopology_ElementCheckExists(domainTopology%domain%decomposition%topology,userElementNumber, &
      & userElementExists,localElementNumber,ghostElement,err,error,*999)
    IF(.NOT.userElementExists) THEN
      CALL FlagError("The specified user element number of "// &
        & TRIM(NumberToVstring(userElementNumber,"*",err,error))// &
        & " does not exist in the domain decomposition.",err,error,*999)
    END IF
    
    basis=>domainTopology%elements%elements(localElementNumber)%basis
    IF(.NOT.ASSOCIATED(basis)) THEN
      localError="The basis for user element number "//TRIM(NumberToVString(userElementNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("DomainTopology_ElementBasisGet")
    RETURN
999 ERRORSEXITS("DomainTopology_ElementBasisGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainTopology_ElementBasisGet

  !
  !================================================================================================================================
  !

  !>Finalises a line in the given domain topology and deallocates all memory.
  SUBROUTINE DOMAIN_TOPOLOGY_LINE_FINALISE(LINE,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainLineType) :: LINE !<The domain line to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_LINE_FINALISE",ERR,ERROR,*999)

    LINE%NUMBER=0
    NULLIFY(LINE%BASIS)
    IF(ALLOCATED(LINE%nodesInLine)) DEALLOCATE(LINE%nodesInLine)
    IF(ALLOCATED(LINE%derivativesInLine)) DEALLOCATE(LINE%derivativesInLine)
 
    EXITS("DOMAIN_TOPOLOGY_LINE_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_LINE_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_LINE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the line data structure for a domain topology line.
  SUBROUTINE DOMAIN_TOPOLOGY_LINE_INITIALISE(LINE,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainLineType) :: LINE !<The domain line to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_LINE_INITIALISE",ERR,ERROR,*999)

    LINE%NUMBER=0
    NULLIFY(LINE%BASIS)
    LINE%boundaryLine=.FALSE.
    
    EXITS("DOMAIN_TOPOLOGY_LINE_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_LINE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_LINE_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the lines in the given domain topology. \todo pass in domain lines
  SUBROUTINE DOMAIN_TOPOLOGY_LINES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: TOPOLOGY !<A pointer to the domain topology to finalise the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nl
    
    ENTERS("DOMAIN_TOPOLOGY_LINES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%LINES)) THEN
        DO nl=1,TOPOLOGY%LINES%numberOfLines
          CALL DOMAIN_TOPOLOGY_LINE_FINALISE(TOPOLOGY%LINES%LINES(nl),ERR,ERROR,*999)
        ENDDO !nl
        IF(ALLOCATED(TOPOLOGY%LINES%LINES)) DEALLOCATE(TOPOLOGY%LINES%LINES)
        DEALLOCATE(TOPOLOGY%LINES)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_TOPOLOGY_LINES_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_LINES_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_LINES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the line data structures for a domain topology. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_LINES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_LINES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%LINES)) THEN
        CALL FlagError("Decomposition already has topology lines associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%LINES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology lines",ERR,ERROR,*999)
        TOPOLOGY%LINES%numberOfLines=0
        TOPOLOGY%LINES%DOMAIN=>TOPOLOGY%DOMAIN
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_LINES_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_LINES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_LINES_INITIALISE
  
  !
  !================================================================================================================================
  !
  !>Finalises a face in the given domain topology and deallocates all memory.
  SUBROUTINE DOMAIN_TOPOLOGY_FACE_FINALISE(FACE,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainFaceType) :: FACE !<The domain face to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_FACE_FINALISE",ERR,ERROR,*999)

    FACE%NUMBER=0
    NULLIFY(FACE%BASIS)
    IF(ALLOCATED(FACE%nodesInFace)) DEALLOCATE(FACE%nodesInFace)
    IF(ALLOCATED(FACE%derivativesInFace)) DEALLOCATE(FACE%derivativesInFace)
 
    EXITS("DOMAIN_TOPOLOGY_FACE_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_FACE_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_FACE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the face data structure for a domain topology face.
  SUBROUTINE DOMAIN_TOPOLOGY_FACE_INITIALISE(FACE,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainFaceType) :: FACE !<The domain face to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_FACE_INITIALISE",ERR,ERROR,*999)

    FACE%NUMBER=0
    NULLIFY(FACE%BASIS)
    FACE%boundaryFace=.FALSE.
    
    EXITS("DOMAIN_TOPOLOGY_FACE_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_FACE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_FACE_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the faces in the given domain topology. \todo pass in domain faces
  SUBROUTINE DOMAIN_TOPOLOGY_FACES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: TOPOLOGY !<A pointer to the domain topology to finalise the faces for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nf
    
    ENTERS("DOMAIN_TOPOLOGY_FACES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%FACES)) THEN
        DO nf=1,TOPOLOGY%FACES%numberOfFaces
          CALL DOMAIN_TOPOLOGY_FACE_FINALISE(TOPOLOGY%FACES%FACES(nf),ERR,ERROR,*999)
        ENDDO !nf
        IF(ALLOCATED(TOPOLOGY%FACES%FACES)) DEALLOCATE(TOPOLOGY%FACES%FACES)
        DEALLOCATE(TOPOLOGY%FACES)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_TOPOLOGY_FACES_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_FACES_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_FACES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the face data structures for a domain topology. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_FACES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the faces for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_FACES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%FACES)) THEN
        CALL FlagError("Decomposition already has topology faces associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%FACES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology faces",ERR,ERROR,*999)
        TOPOLOGY%FACES%numberOfFaces=0
        TOPOLOGY%FACES%DOMAIN=>TOPOLOGY%DOMAIN
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_FACES_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_FACES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_FACES_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Checks that a user node number exists in a domain nodes topolgoy. 
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS(DOMAIN_TOPOLOGY,USER_NODE_NUMBER,NODE_EXISTS,DOMAIN_LOCAL_NODE_NUMBER, &
    & GHOST_NODE,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: DOMAIN_TOPOLOGY !<A pointer to the domain topology to check the node exists on
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<The user node number to check if it exists
    LOGICAL, INTENT(OUT) :: NODE_EXISTS !<On exit, is .TRUE. if the node user number exists in the domain nodes topolgoy (even if it is a ghost node), .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: DOMAIN_LOCAL_NODE_NUMBER !<On exit, if the node exists the local number corresponding to the user node number. If the node does not exist then global number will be 0.
    LOGICAL, INTENT(OUT) :: GHOST_NODE !<On exit, is .TRUE. if the local node (if it exists) is a ghost node, .FALSE. if not. 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DomainNodesType), POINTER :: DOMAIN_NODES
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    
    ENTERS("DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS",ERR,ERROR,*999)

    NODE_EXISTS=.FALSE.
    DOMAIN_LOCAL_NODE_NUMBER=0
    GHOST_NODE=.FALSE.
    IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
      DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
      IF(ASSOCIATED(DOMAIN_NODES)) THEN
        NULLIFY(TREE_NODE)
        CALL TREE_SEARCH(DOMAIN_NODES%nodesTree,USER_NODE_NUMBER,TREE_NODE,ERR,ERROR,*999)
        IF(ASSOCIATED(TREE_NODE)) THEN
          CALL TREE_NODE_VALUE_GET(DOMAIN_NODES%nodesTree,TREE_NODE,DOMAIN_LOCAL_NODE_NUMBER,ERR,ERROR,*999)
          NODE_EXISTS=.TRUE.
          GHOST_NODE=DOMAIN_LOCAL_NODE_NUMBER>DOMAIN_NODES%numberOfNodes
        ENDIF
      ELSE
        CALL FlagError("Domain topology nodes is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS
  
  !
  !================================================================================================================================
  !

  !>Finalises the given domain topology node derivative and deallocates all memory.
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE(NODE_DERIVATIVE,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainNodeDerivativeType) :: NODE_DERIVATIVE !<The domain node to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(NODE_DERIVATIVE%userVersionNumbers)) DEALLOCATE(NODE_DERIVATIVE%userVersionNumbers)
    IF(ALLOCATED(NODE_DERIVATIVE%dofIndex)) DEALLOCATE(NODE_DERIVATIVE%dofIndex)

    EXITS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given mesh topology node.
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE(NODE_DERIVATIVE,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainNodeDerivativeType) :: NODE_DERIVATIVE !<The domain node to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE",ERR,ERROR,*999)

    NODE_DERIVATIVE%numberOfVersions=0
    NODE_DERIVATIVE%globalDerivativeIndex=0
    NODE_DERIVATIVE%partialDerivativeIndex=0

    EXITS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the given domain topology node and deallocates all memory.
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_FINALISE(NODE,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainNodeType) :: NODE !<The domain node to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: derivative_idx

    ENTERS("DOMAIN_TOPOLOGY_NODE_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(NODE%DERIVATIVES)) THEN
      DO derivative_idx=1,NODE%numberOfDerivatives
        CALL DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE(NODE%DERIVATIVES(derivative_idx),ERR,ERROR,*999)
      ENDDO !derivative_idx
      DEALLOCATE(NODE%DERIVATIVES)
    ENDIF
    IF(ASSOCIATED(NODE%surroundingElements)) DEALLOCATE(NODE%surroundingElements)
    IF(ALLOCATED(NODE%nodeLines)) DEALLOCATE(NODE%nodeLines)
 
    EXITS("DOMAIN_TOPOLOGY_NODE_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODE_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given domain topology node.
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_INITIALISE(NODE,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainNodeType) :: NODE !<The domain node to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_NODE_INITIALISE",ERR,ERROR,*999)

    NODE%localNumber=0
    NODE%meshNumber=0
    NODE%globalNumber=0
    NODE%userNumber=0
    NODE%numberOfSurroundingElements=0
    NODE%numberOfNodeLines=0
    NODE%boundaryNode=.FALSE.
    
    EXITS("DOMAIN_TOPOLOGY_NODE_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODE_INITIALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the nodees in the given domain topology. \todo pass in domain nodes
  SUBROUTINE DOMAIN_TOPOLOGY_NODES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: np

    ENTERS("DOMAIN_TOPOLOGY_NODES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
        DO np=1,TOPOLOGY%NODES%totalNumberOfNodes
          CALL DOMAIN_TOPOLOGY_NODE_FINALISE(TOPOLOGY%NODES%NODES(np),ERR,ERROR,*999)
        ENDDO !np
        IF(ASSOCIATED(TOPOLOGY%NODES%NODES)) DEALLOCATE(TOPOLOGY%NODES%NODES)
        IF(ASSOCIATED(TOPOLOGY%NODES%nodesTree)) CALL TREE_DESTROY(TOPOLOGY%NODES%nodesTree,ERR,ERROR,*999)
        DEALLOCATE(TOPOLOGY%NODES)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_TOPOLOGY_NODES_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODES_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_NODES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the nodes data structures for a domain topology. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_NODES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_NODES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
        CALL FlagError("Decomposition already has topology nodes associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%NODES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology nodes",ERR,ERROR,*999)
        TOPOLOGY%NODES%numberOfNodes=0
        TOPOLOGY%NODES%totalNumberOfNodes=0
        TOPOLOGY%NODES%numberOfGlobalNodes=0
        TOPOLOGY%NODES%maximumNumberOfDerivatives=0
        TOPOLOGY%NODES%DOMAIN=>TOPOLOGY%DOMAIN
        NULLIFY(TOPOLOGY%NODES%NODES)
        NULLIFY(TOPOLOGY%NODES%nodesTree)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_NODES_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_NODES_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Calculates the element numbers surrounding a node for a domain.
  SUBROUTINE DomainTopology_NodesSurroundingElementsCalculate(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: TOPOLOGY !<A pointer to the domain topology to calculate the elements surrounding the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_no,insert_position,ne,nn,np,surrounding_elem_no
    INTEGER(INTG), POINTER :: NEW_SURROUNDING_ELEMENTS(:)
    LOGICAL :: FOUND_ELEMENT
    TYPE(BasisType), POINTER :: BASIS

    NULLIFY(NEW_SURROUNDING_ELEMENTS)

    ENTERS("DomainTopology_NodesSurroundingElementsCalculate",ERR,ERROR,*999)
    
    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
          IF(ASSOCIATED(TOPOLOGY%NODES%NODES)) THEN
            DO np=1,TOPOLOGY%NODES%totalNumberOfNodes
              TOPOLOGY%NODES%NODES(np)%numberOfSurroundingElements=0
              IF(ASSOCIATED(TOPOLOGY%NODES%NODES(np)%surroundingElements)) &
                & DEALLOCATE(TOPOLOGY%NODES%NODES(np)%surroundingElements)
            ENDDO !np
            DO ne=1,TOPOLOGY%ELEMENTS%totalNumberOfElements
              BASIS=>TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS
              DO nn=1,BASIS%numberOfNodes
                np=TOPOLOGY%ELEMENTS%ELEMENTS(ne)%elementNodes(nn)
                FOUND_ELEMENT=.FALSE.
                element_no=1
                insert_position=1
                DO WHILE(element_no<=TOPOLOGY%NODES%NODES(np)%numberOfSurroundingElements.AND..NOT.FOUND_ELEMENT)
                  surrounding_elem_no=TOPOLOGY%NODES%NODES(np)%surroundingElements(element_no)
                  IF(surrounding_elem_no==ne) THEN
                    FOUND_ELEMENT=.TRUE.
                  ENDIF
                  element_no=element_no+1
                  IF(ne>=surrounding_elem_no) THEN
                    insert_position=element_no
                  ENDIF
                ENDDO
                IF(.NOT.FOUND_ELEMENT) THEN
                  !Insert element into surrounding elements
                  ALLOCATE(NEW_SURROUNDING_ELEMENTS(TOPOLOGY%NODES%NODES(np)%numberOfSurroundingElements+1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate new surrounding elements",ERR,ERROR,*999)
                  IF(ASSOCIATED(TOPOLOGY%NODES%NODES(np)%surroundingElements)) THEN
                    NEW_SURROUNDING_ELEMENTS(1:insert_position-1)=TOPOLOGY%NODES%NODES(np)% &
                      & surroundingElements(1:insert_position-1)
                    NEW_SURROUNDING_ELEMENTS(insert_position)=ne
                    NEW_SURROUNDING_ELEMENTS(insert_position+1:TOPOLOGY%NODES%NODES(np)% &
                      & numberOfSurroundingElements+1)=TOPOLOGY%NODES%NODES(np)%surroundingElements(insert_position: &
                      & TOPOLOGY%NODES%NODES(np)%numberOfSurroundingElements)
                    DEALLOCATE(TOPOLOGY%NODES%NODES(np)%surroundingElements)
                  ELSE
                    NEW_SURROUNDING_ELEMENTS(1)=ne
                  ENDIF
                  TOPOLOGY%NODES%NODES(np)%surroundingElements=>NEW_SURROUNDING_ELEMENTS
                  TOPOLOGY%NODES%NODES(np)%numberOfSurroundingElements= &
                    & TOPOLOGY%NODES%NODES(np)%numberOfSurroundingElements+1
                ENDIF
              ENDDO !nn
            ENDDO !ne
          ELSE
            CALL FlagError("Domain topology nodes nodes are not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Domain topology nodes are not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Domain topology elements is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain topology is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("DomainTopology_NodesSurroundingElementsCalculate")
    RETURN
999 IF(ASSOCIATED(NEW_SURROUNDING_ELEMENTS)) DEALLOCATE(NEW_SURROUNDING_ELEMENTS)
    ERRORS("DomainTopology_NodesSurroundingElementsCalculate",ERR,ERROR)
    EXITS("DomainTopology_NodesSurroundingElementsCalculate")
    RETURN 1   
  END SUBROUTINE DomainTopology_NodesSurroundingElementsCalculate
  
  !
  !================================================================================================================================
  !

  
END MODULE DecompositionRoutines
