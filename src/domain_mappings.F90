!> \file
!> \author Chris Bradley
!> \brief This module handles all domain mappings routines.
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

!> This module handles all domain mappings routines.
MODULE DomainMappings

  USE BaseRoutines
  USE ComputationRoutines
  USE ComputationAccessRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE Lists
  USE Strings
  USE Types

#include "macros.h"  


  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup DomainMappings_DomainType DomainMappings::DomainType
  !> \see DomainMappings
  !>@{
  INTEGER(INTG), PARAMETER :: DOMAIN_LOCAL_INTERNAL=1 !<The domain item is internal to the domain \see DomainMappings_DomainType,DomainMappings
  INTEGER(INTG), PARAMETER :: DOMAIN_LOCAL_BOUNDARY=2 !<The domain item is on the boundary of the domain \see DomainMappings_DomainType,DomainMappings
  INTEGER(INTG), PARAMETER :: DOMAIN_LOCAL_GHOST=3 !<The domain item is ghosted from another domain \see DomainMappings_DomainType,DomainMappings
  !>@}
  
  !Module types

  !Module variables

  !Interfaces
  
  PUBLIC DOMAIN_LOCAL_INTERNAL,DOMAIN_LOCAL_BOUNDARY,DOMAIN_LOCAL_GHOST
  
  PUBLIC DomainGlobalMapping_Initialise

  PUBLIC DomainMapping_BoundaryFinishGet

  PUBLIC DomainMapping_BoundaryStartGet
  
  PUBLIC DomainMapping_DomainNumberFromGlobalGet
  
  PUBLIC DomainMapping_Finalise,DomainMapping_Initialise

  PUBLIC DomainMapping_GhostFinishGet

  PUBLIC DomainMapping_GhostStartGet

  PUBLIC DomainMapping_InternalFinishGet

  PUBLIC DomainMapping_InternalStartGet

  PUBLIC DomainMapping_LocalToGlobalGet
  
  PUBLIC DomainMapping_LocalFromGlobalCalculate

  PUBLIC DomainMapping_LocalNumberFromGlobalGet

  PUBLIC DomainMapping_LocalTypeFromGlobalGet

  PUBLIC DomainMapping_NumberGet

  PUBLIC DomainMapping_NumberOfBoundaryGet

  PUBLIC DomainMapping_NumberOfDomainsFromGlobalGet

  PUBLIC DomainMapping_NumberOfInternalGet

  PUBLIC DomainMapping_NumberOfGhostGet

  PUBLIC DomainMapping_NumberOfGlobalGet

  PUBLIC DomainMapping_NumberOfLocalGet

  PUBLIC DomainMapping_TotalNumberOfLocalGet

  PUBLIC DomainMapping_WorkGroupGet,DomainMapping_WorkGroupSet

  !PUBLIC DomainMappings_MappingCalculate

CONTAINS
  
  !
  !================================================================================================================================
  !

  !>Finalises the adjacent domain and deallocates all memory for a domain mapping.
  SUBROUTINE DomainAdjacentDomain_Finalise(adjacentDomain,err,error,*)

    !Argument variables
    TYPE(DomainAdjacentDomainType) :: adjacentDomain !<The adjacent domain to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainAdjacentDomain_Finalise",err,error,*999)

    IF(ALLOCATED(adjacentDomain%localGhostSendIndices)) DEALLOCATE(adjacentDomain%localGhostSendIndices)
    IF(ALLOCATED(adjacentDomain%localGhostReceiveIndices)) DEALLOCATE(adjacentDomain%localGhostReceiveIndices)
    adjacentDomain%numberOfSendGhosts=0
    adjacentDomain%numberOfReceiveGhosts=0
    adjacentDomain%domainNumber=0
    
    EXITS("DomainAdjacentDomain_Finalise")
    RETURN
999 ERRORSEXITS("DomainAdjacentDomain_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE DomainAdjacentDomain_Finalise
  
  !
  !================================================================================================================================
  !

  !>Initialise the adjacent domain for a domain mapping.
  SUBROUTINE DomainAdjacentDomain_Initialise(adjacentDomain,err,error,*)

    !Argument variables
    TYPE(DomainAdjacentDomainType) :: adjacentDomain !<The adjacent domain to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainAdjacentDomain_Initialise",err,error,*999)

    adjacentDomain%numberOfSendGhosts=0
    adjacentDomain%numberOfReceiveGhosts=0
    adjacentDomain%domainNumber=0
    
    EXITS("DomainAdjacentDomain_Initialise")
    RETURN
999 ERRORSEXITS("DomainAdjacentDomain_Initialise",err,error)
    RETURN 1
  END SUBROUTINE DomainAdjacentDomain_Initialise
  
  !
  !================================================================================================================================
  !

  !> Finalises the global mapping in the given domain mappings.
  SUBROUTINE DomainGlobalMapping_Finalise(globalMapping,err,error,*)

    !Argument variables
    TYPE(DomainGlobalMappingType) :: globalMapping !<The domain global mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<Th error string
    !Local Variables

    ENTERS("DomainGlobalMapping_Finalise",err,error,*999)

    IF(ALLOCATED(globalMapping%localNumber)) DEALLOCATE(globalMapping%localNumber)
    IF(ALLOCATED(globalMapping%domainNumber)) DEALLOCATE(globalMapping%domainNumber)
    IF(ALLOCATED(globalMapping%localType)) DEALLOCATE(globalMapping%localType)
 
    EXITS("DomainGlobalMapping_Finalise")
    RETURN
999 ERRORSEXITS("DomainGlobalMapping_Finalise",err,error)
    RETURN 1
   
  END SUBROUTINE DomainGlobalMapping_Finalise

  !
  !================================================================================================================================
  !

  !>Gets the boundary finish for a domain mapping.
  SUBROUTINE DomainMapping_BoundaryFinishGet(domainMapping,boundaryFinish,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainMapping !<A pointer to the domain mapping to get the boundary finish for
    INTEGER(INTG), INTENT(OUT) :: boundaryFinish !<On exit, the boundary finish number for the domain mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainMapping_BoundaryFinishGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*999)
#endif    

    boundaryFinish=domainMapping%boundaryFinish
    
    EXITS("DomainMapping_BoundaryFinishGet")
    RETURN
999 ERRORSEXITS("DomainMapping_BoundaryFinishGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMapping_BoundaryFinishGet
  
  !
  !================================================================================================================================
  !

  !>Gets the boundary start for a domain mapping.
  SUBROUTINE DomainMapping_BoundaryStartGet(domainMapping,boundaryStart,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainMapping !<A pointer to the domain mapping to get the boundary start for
    INTEGER(INTG), INTENT(OUT) :: boundaryStart !<On exit, the boundary start number for the domain mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainMapping_BoundaryStartGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*999)
#endif    

    boundaryStart=domainMapping%boundaryStart
    
    EXITS("DomainMapping_BoundaryStartGet")
    RETURN
999 ERRORSEXITS("DomainMapping_BoundaryStartGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMapping_BoundaryStartGet
  
  !
  !================================================================================================================================
  !

  !>Gets the domain number from a global number in a domain in a domain mapping.
  SUBROUTINE DomainMapping_DomainNumberFromGlobalGet(domainMapping,globalNumber,domainIdx,domainNumber,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainMapping !<A pointer to the domain mapping to get the domain number for
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number to get the domain number for
    INTEGER(INTG), INTENT(IN) :: domainIdx !<The domain index to get the domain number for
    INTEGER(INTG), INTENT(OUT) :: domainNumber !<On exit, the domain number for the global number and domain index the domain mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainMapping_DomainNumberFromGlobalGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(domainMapping%globalToLocalMap)) &
      & CALL FlagError("The global to local map is not allocated in the domain mapping.",err,error,*999)
    IF(globalNumber<1.OR.globalNumber>domainMapping%numberOfGlobal) THEN
      localError="The specified global number of "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
        & " is invalid for the domain mapping. The global number should be >= 1 and <= "// &
        & TRIM(NumberToVString(domainMapping%numberOfGlobal,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(domainIdx<1.OR.domainIdx>domainMapping%globalToLocalMap(globalNumber)%numberOfDomains) THEN
      localError="The specified domain index of "//TRIM(NumberToVString(domainIdx,"*",err,error))// &
        & " is invalid for global number "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
        & " in the domain mapping. The domain index should be >= 1 and <= "// &
        & TRIM(NumberToVString(domainMapping%globalToLocalMap(globalNumber)%numberOfDomains,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    domainNumber=domainMapping%globalToLocalMap(globalNumber)%domainNumber(domainIdx)
    
    EXITS("DomainMapping_DomainNumberFromGlobalGet")
    RETURN
999 ERRORSEXITS("DomainMapping_DomainNumberFromGlobalGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMapping_DomainNumberFromGlobalGet
  
  !
  !================================================================================================================================
  !

  !>Initialises the global mapping in the given domain mappings.
  SUBROUTINE DomainGlobalMapping_Initialise(globalMapping,err,error,*)

    !Argument variables
    TYPE(DomainGlobalMappingType) :: globalMapping !<The domain global mapping to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainGlobalMapping_Initialise",err,error,*999)

    globalMapping%numberOfDomains=0
    
    EXITS("DomainGlobalMapping_Initialise")
    RETURN
999 ERRORSEXITS("DomainGlobalMapping_Initialise",err,error)
    RETURN 1
   
  END SUBROUTINE DomainGlobalMapping_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises a domain mapping and deallocates all memory.
  SUBROUTINE DomainMapping_Finalise(domainMapping,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainMapping !<A pointer to the domain mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: idx

    ENTERS("DomainMapping_Finalise",err,error,*999)

    IF(ASSOCIATED(domainMapping)) THEN
      IF(ALLOCATED(domainMapping%numberOfDomainLocal)) DEALLOCATE(domainMapping%numberOfDomainLocal)
      IF(ALLOCATED(domainMapping%numberOfDomainGhost)) DEALLOCATE(domainMapping%numberOfDomainGhost)
      IF(ALLOCATED(domainMapping%domainList)) DEALLOCATE(domainMapping%domainList)
      IF(ALLOCATED(domainMapping%localToGlobalMap)) DEALLOCATE(domainMapping%localToGlobalMap)
      IF(ALLOCATED(domainMapping%globalToLocalMap)) THEN
        DO idx=1,SIZE(domainMapping%globalToLocalMap,1)
          CALL DomainGlobalMapping_Finalise(domainMapping%globalToLocalMap(idx),err,error,*999)
        ENDDO !idx
        DEALLOCATE(domainMapping%globalToLocalMap)
      ENDIF
      IF(ALLOCATED(domainMapping%adjacentDomainsPtr)) DEALLOCATE(domainMapping%adjacentDomainsPtr)
      IF(ALLOCATED(domainMapping%adjacentDomainsList)) DEALLOCATE(domainMapping%adjacentDomainsList)
      IF(ALLOCATED(domainMapping%adjacentDomains)) THEN
        DO idx=1,SIZE(domainMapping%adjacentDomains,1)
          CALL DomainAdjacentDomain_Finalise(domainMapping%adjacentDomains(idx),err,error,*999)
        ENDDO !idx
        DEALLOCATE(domainMapping%adjacentDomains)
      ENDIF
      DEALLOCATE(domainMapping)
    ENDIF
    
    EXITS("DomainMapping_Finalise")
    RETURN
999 ERRORSEXITS("DomainMapping_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMapping_Finalise
  
  !
  !================================================================================================================================
  !

  !>Initialises the domain mapping.
  SUBROUTINE DomainMapping_Initialise(domainMapping,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainMapping !<A pointer to the domain mapping to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
  
    ENTERS("DomainMapping_Initialise",err,error,*998)

    IF(ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is already associated.",err,error,*998)

    ALLOCATE(domainMapping,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate domain mapping.",err,error,*999)

    NULLIFY(domainMapping%workGroup)
    domainMapping%totalNumberOfLocal=0
    domainMapping%numberOfLocal=0
    domainMapping%numberOfGlobal=0
    domainMapping%numberOfInternal=0
    domainMapping%numberOfBoundary=0
    domainMapping%numberOfGhost=0
    domainMapping%internalStart=0
    domainMapping%internalFinish=0
    domainMapping%boundaryStart=0
    domainMapping%boundaryFinish=0
    domainMapping%ghostStart=0
    domainMapping%ghostFinish=0
    domainMapping%numberOfAdjacentDomains=0
    
    EXITS("DomainMapping_Initialise")
    RETURN
999 CALL DomainMapping_Finalise(domainMapping,dummyErr,dummyError,*999)
998 ERRORSEXITS("DomainMapping_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMapping_Initialise

  !
  !================================================================================================================================
  !

  !>Gets the ghost finish for a domain mapping.
  SUBROUTINE DomainMapping_GhostFinishGet(domainMapping,ghostFinish,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainMapping !<A pointer to the domain mapping to get the ghost finish for
    INTEGER(INTG), INTENT(OUT) :: ghostFinish !<On exit, the ghost finish number for the domain mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainMapping_GhostFinishGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*999)
#endif    

    ghostFinish=domainMapping%ghostFinish
    
    EXITS("DomainMapping_GhostFinishGet")
    RETURN
999 ERRORSEXITS("DomainMapping_GhostFinishGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMapping_GhostFinishGet
  
  !
  !================================================================================================================================
  !

  !>Gets the ghost start for a domain mapping.
  SUBROUTINE DomainMapping_GhostStartGet(domainMapping,ghostStart,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainMapping !<A pointer to the domain mapping to get the ghost start for
    INTEGER(INTG), INTENT(OUT) :: ghostStart !<On exit, the ghost start number for the domain mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainMapping_GhostStartGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*999)
#endif    

    ghostStart=domainMapping%ghostStart
    
    EXITS("DomainMapping_GhostStartGet")
    RETURN
999 ERRORSEXITS("DomainMapping_GhostStartGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMapping_GhostStartGet
  
  !
  !================================================================================================================================
  !

  !>Gets the internal finish for a domain mapping.
  SUBROUTINE DomainMapping_InternalFinishGet(domainMapping,internalFinish,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainMapping !<A pointer to the domain mapping to get the internal finish for
    INTEGER(INTG), INTENT(OUT) :: internalFinish !<On exit, the internal finish number for the domain mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainMapping_InternalFinishGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*999)
#endif    

    internalFinish=domainMapping%internalFinish
    
    EXITS("DomainMapping_InternalFinishGet")
    RETURN
999 ERRORSEXITS("DomainMapping_InternalFinishGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMapping_InternalFinishGet
  
  !
  !================================================================================================================================
  !

  !>Gets the internal start for a domain mapping.
  SUBROUTINE DomainMapping_InternalStartGet(domainMapping,internalStart,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainMapping !<A pointer to the domain mapping to get the internal start for
    INTEGER(INTG), INTENT(OUT) :: internalStart !<On exit, the internal start number for the domain mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainMapping_InternalStartGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*999)
#endif    

    internalStart=domainMapping%internalStart
    
    EXITS("DomainMapping_InternalStartGet")
    RETURN
999 ERRORSEXITS("DomainMapping_InternalStartGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMapping_InternalStartGet
  
  !
  !================================================================================================================================
  !

  !>Gets the global number from a local number for a domain mapping.
  SUBROUTINE DomainMapping_LocalToGlobalGet(domainMapping,localNumber,globalNumber,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainMapping !<A pointer to the domain mapping to get the global number for
    INTEGER(INTG), INTENT(IN) :: localNumber !<The local number to get the global number for
    INTEGER(INTG), INTENT(OUT) :: globalNumber !<On exit, the global number corresponding to the local number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainMapping_LocalToGlobalGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*999)
    IF(localNumber<1.OR.localNumber>domainMapping%totalNumberOfLocal) THEN
      localError="The specified local number of "//TRIM(NumberToVString(localNumber,"*",err,error))// &
        & " is invalid for the domain mapping. The local number should be >= 1 and <= "// &
        & TRIM(NumberToVString(domainMapping%totalNumberOfLocal,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainMapping%localToGlobalMap)) &
      & CALL FlagError("The local to global map array is not allocated for the domain mapping.",err,error,*999)
#endif    

    globalNumber=domainMapping%localToGlobalMap(localNumber)
    
    EXITS("DomainMapping_LocalToGlobalGet")
    RETURN
999 ERRORSEXITS("DomainMapping_LocalToGlobalGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMapping_LocalToGlobalGet
  
  !
  !================================================================================================================================
  !

  !>Calculates the domain mappings local map from a domain mappings global map.
  SUBROUTINE DomainMapping_LocalFromGlobalCalculate(domainMapping,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainMapping !<The domain mapping to calculate the local mappings
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: domainIdx,domainIdx2,domainNumber,domainNumber2,globalNumber,idx,localNumber,localNumber2,numberOfInternal, &
      & numberOfBoundary,numberOfGhost,myGroupComputationNodeNumber,myDomainIndex,temp,numberOfAdjacentDomains, &
      & numberOfGroupComputationNodes,receiveFromDomain,dummyErr,numberOfGhostReceive,numberOfGhostSend,localType,count, &
      & totalNumberOfAdjacentDomains
    INTEGER(INTG), ALLOCATABLE :: adjacentDomainMap(:),adjacentDomains(:,:),sendList(:),receiveList(:)
    LOGICAL :: ownedByAll,sendGlobal
    TYPE(ListPtrType), ALLOCATABLE :: ghostSendLists(:),ghostReceiveLists(:)
    TYPE(VARYING_STRING) :: localError,dummyError
    TYPE(WorkGroupType), POINTER :: workGroup
    
    ENTERS("DomainMapping_LocalFromGlobalCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*999)

    NULLIFY(workGroup)
    CALL DomainMapping_WorkGroupGet(domainMapping,workGroup,err,error,*999)
    CALL WorkGroup_GroupNodeNumberGet(workGroup,myGroupComputationNodeNumber,err,error,*999)
    CALL WorkGroup_NumberOfGroupNodesGet(workGroup,numberOfGroupComputationNodes,err,error,*999)

    !Calculate local to global maps from global to local map
    ALLOCATE(domainMapping%numberOfDomainLocal(0:numberOfGroupComputationNodes-1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate number of domain local.",err,error,*999)
    domainMapping%numberOfDomainLocal=0

    ALLOCATE(domainMapping%numberOfDomainGhost(0:numberOfGroupComputationNodes-1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate number of domain ghost.",err,error,*999)

    domainMapping%numberOfDomainGhost=0
    numberOfInternal=0   !counters for my computational node
    numberOfBoundary=0
    numberOfGhost=0

    ALLOCATE(adjacentDomains(0:numberOfGroupComputationNodes-1,0:numberOfGroupComputationNodes-1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate adjacent domains.",err,error,*999)
    adjacentDomains=0

    !Loop over global elements
    DO globalNumber=1,domainMapping%numberOfGlobal
      !-------- If necessary, reset global domain index so that my computational node is in the first index position ------------
      !find out current domain index of my computational node
      myDomainIndex=1
      DO domainIdx=2,domainMapping%globalToLocalMap(globalNumber)%numberOfDomains
        domainNumber=domainMapping%globalToLocalMap(globalNumber)%domainNumber(domainIdx)
        IF(domainNumber==myGroupComputationNodeNumber) THEN
          myDomainIndex=domainIdx
          EXIT
        ENDIF
      ENDDO !domainIdx

      !if necessary, swap the data structure at the current domain index for my computational node with domain index 1
      IF(myDomainIndex/=1) THEN
        !Swap domain index in the global to local map, 1 <-> myDomainIndex
        temp=domainMapping%globalToLocalMap(globalNumber)%localNumber(1)
        domainMapping%globalToLocalMap(globalNumber)%localNumber(1) = &
          & domainMapping%globalToLocalMap(globalNumber)%localNumber(myDomainIndex)
        domainMapping%globalToLocalMap(globalNumber)%localNumber(myDomainIndex) = temp

        temp=domainMapping%globalToLocalMap(globalNumber)%domainNumber(1)
        domainMapping%globalToLocalMap(globalNumber)%domainNumber(1) = &
          & domainMapping%globalToLocalMap(globalNumber)%domainNumber(myDomainIndex)
        domainMapping%globalToLocalMap(globalNumber)%domainNumber(myDomainIndex) = temp

        temp=domainMapping%globalToLocalMap(globalNumber)%localType(1)
        domainMapping%globalToLocalMap(globalNumber)%localType(1) = &
          & domainMapping%globalToLocalMap(globalNumber)%localType(myDomainIndex)
        domainMapping%globalToLocalMap(globalNumber)%localType(myDomainIndex) = temp
      ENDIF

      !Set the global adjacent_domains array to 1 for domains that share the current element
      DO domainIdx=1,domainMapping%globalToLocalMap(globalNumber)%numberOfDomains
        domainNumber=domainMapping%globalToLocalMap(globalNumber)%domainNumber(domainIdx)
        DO domainIdx2=1,domainMapping%globalToLocalMap(globalNumber)%numberOfDomains
          domainNumber2=domainMapping%globalToLocalMap(globalNumber)%domainNumber(domainIdx2)
          adjacentDomains(domainNumber,domainNumber2)=1
        ENDDO !domainIdx2
      ENDDO !domainIdx

      !Loop over domains where current global is
      DO domainIdx=1,domainMapping%globalToLocalMap(globalNumber)%numberOfDomains
        domainNumber=domainMapping%globalToLocalMap(globalNumber)%domainNumber(domainIdx)
        localType=domainMapping%globalToLocalMap(globalNumber)%localType(domainIdx)
        !increment counter of total ghost/local per domain
        IF(localType==DOMAIN_LOCAL_GHOST) THEN
          domainMapping%numberOfDomainGhost(domainNumber)=domainMapping%numberOfDomainGhost(domainNumber)+1
        ELSE
          domainMapping%numberOfDomainLocal(domainNumber)=domainMapping%numberOfDomainLocal(domainNumber)+1
        ENDIF

        !Increment counter of internal, boundary and ghost elements on my domain
        IF(domainNumber==myGroupComputationNodeNumber) THEN
          SELECT CASE(localType)
          CASE(DOMAIN_LOCAL_INTERNAL)
            numberOfInternal=numberOfInternal+1
          CASE(DOMAIN_LOCAL_BOUNDARY)
            numberOfBoundary=numberOfBoundary+1
          CASE(DOMAIN_LOCAL_GHOST)
            numberOfGhost=numberOfGhost+1
          CASE DEFAULT
            localError="The domain local type of "//TRIM(NumberToVString(domainMapping%globalToLocalMap(globalNumber)% &
              & localType(domainIdx),"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDIF
      ENDDO !domainIdx
    ENDDO !globalNumber

!!TODO: move adjacent domains calculation back to where the global to local array is set up????
    !count number of adjacent domains for my domain and in total
    numberOfAdjacentDomains=0
    totalNumberOfAdjacentDomains=0
    DO domainNumber=0,numberOfGroupComputationNodes-1
      DO domainNumber2=0,numberOfGroupComputationNodes-1
        IF(domainNumber/=domainNumber2) THEN 
          IF(adjacentDomains(domainNumber,domainNumber2)>0) THEN
            totalNumberOfAdjacentDomains=totalNumberOfAdjacentDomains+1
            IF(domainNumber==myGroupComputationNodeNumber) numberOfAdjacentDomains=numberOfAdjacentDomains+1
          ENDIF
        ENDIF
      ENDDO !domainNumber2
    ENDDO !domainNumber

    ALLOCATE(domainMapping%adjacentDomainsPtr(0:numberOfGroupComputationNodes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate adjacent domains ptr.",err,error,*999)
    ALLOCATE(domainMapping%adjacentDomainsList(totalNumberOfAdjacentDomains),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate adjacent domains list.",err,error,*999)

    !Store the adjacent domains for a domain in adjacentDomainsList starting at index count and store the starting index count in
    !adjacentDomainsPtr the adjacent domains for a domain are then adjacentDomainsList(adjacentDomainsPtr(domain) : adjacentDomainsPtr(domain+1))
    count=1
    DO domainNumber=0,numberOfGroupComputationNodes-1
      domainMapping%adjacentDomainsPtr(domainNumber)=count
      DO domainNumber2=0,numberOfGroupComputationNodes-1
        IF(domainNumber/=domainNumber2) THEN
          IF(adjacentDomains(domainNumber,domainNumber2)>0) THEN
            domainMapping%adjacentDomainsList(count)=domainNumber2
            count=count+1
          ENDIF
        ENDIF
      ENDDO !domainNumber2
    ENDDO !domainNumber

    domainMapping%adjacentDomainsPtr(numberOfGroupComputationNodes)=count
    DEALLOCATE(adjacentDomains)

    !Compute domain list for my computational node
    ALLOCATE(domainMapping%domainList(numberOfInternal+numberOfBoundary+numberOfGhost),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate domain map domain list.",err,error,*999)
    ALLOCATE(domainMapping%localToGlobalMap(numberOfInternal+numberOfBoundary+numberOfGhost),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate domain map local to global list.",err,error,*999)

    !Set constants
    domainMapping%totalNumberOfLocal=numberOfInternal+numberOfBoundary+numberOfGhost
    domainMapping%numberOfLocal=numberOfInternal+numberOfBoundary
    domainMapping%numberOfInternal=numberOfInternal
    domainMapping%numberOfBoundary=numberOfBoundary
    domainMapping%numberOfGhost=numberOfGhost
    domainMapping%internalStart=1
    domainMapping%internalFinish=numberOfInternal
    domainMapping%boundaryStart=numberOfInternal+1
    domainMapping%boundaryFinish=numberOfInternal+numberOfBoundary
    domainMapping%ghostStart=numberOfInternal+numberOfBoundary+1
    domainMapping%ghostFinish=numberOfInternal+numberOfBoundary+numberOfGhost

    !Adjacent_domains maps a domain index (index between 1 and domainMapping%numberOfAdjacentDomains) to the domain number
    ALLOCATE(domainMapping%adjacentDomains(numberOfAdjacentDomains),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate adjacent domains.",err,error,*999)
    domainMapping%numberOfAdjacentDomains=numberOfAdjacentDomains

    !Adjacent_domain_map maps a domain number back to its index (index between 1 and domainMapping%numberOfAdjacentDomains)
    ALLOCATE(adjacentDomainMap(0:numberOfGroupComputationNodes-1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate adjacent domain map.",err,error,*999)

    !GhostSendLists contains a list for each adjacent domain with elements in my domain (by local element numbers) that will be
    !sent to that foreign domain
    ALLOCATE(ghostSendLists(domainMapping%numberOfAdjacentDomains),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate ghost send list.",err,error,*999)

    !GhostRecieveLists contains a list for each adjacent domain and contains the local numbers of ghost elements that can be
    !received from that foreign domain
    ALLOCATE(ghostReceiveLists(domainMapping%numberOfAdjacentDomains),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate ghost recieve list.",err,error,*999)

    !Set adjacent domains data structures and initialize ghost send and receive lists
    !Loop over adjacent domains of my computational node
    DO domainIdx=1,domainMapping%numberOfAdjacentDomains
      ! set variables to 0
      CALL DomainAdjacentDomain_Initialise(domainMapping%adjacentDomains(domainIdx),err,error,*999)

      !Get number of current adjacent domain
      domainNumber=domainMapping%adjacentDomainsList(domainMapping%adjacentDomainsPtr(myGroupComputationNodeNumber)+domainIdx-1)
      !Set number in adjacent_domains and adjacent_domain_map
      domainMapping%adjacentDomains(domainIdx)%domainNumber=domainNumber
      adjacentDomainMap(domainNumber)=domainIdx

      !Initialize send and receive lists for ghosts
      NULLIFY(ghostSendLists(domainIdx)%ptr)
      CALL List_CreateStart(ghostSendLists(domainIdx)%ptr,err,error,*999)
      CALL List_DataTypeSet(ghostSendLists(domainIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
      CALL List_InitialSizeSet(ghostSendLists(domainIdx)%ptr,MAX(domainMapping%numberOfGhost,1),err,error,*999)
      CALL List_CreateFinish(ghostSendLists(domainIdx)%ptr,err,error,*999)

      NULLIFY(ghostReceiveLists(domainIdx)%ptr)
      CALL List_CreateStart(ghostReceiveLists(domainIdx)%ptr,err,error,*999)
      CALL List_DataTypeSet(ghostReceiveLists(domainIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
      CALL List_InitialSizeSet(ghostReceiveLists(domainIdx)%ptr,MAX(domainMapping%numberOfGhost,1),err,error,*999)
      CALL List_CreateFinish(ghostReceiveLists(domainIdx)%ptr,err,error,*999)
    ENDDO !domainIdx

    numberOfInternal=0
    numberOfBoundary=0
    numberOfGhost=0

    !Loop over global numbers
    DO globalNumber=1,domainMapping%numberOfGlobal
      sendGlobal=.FALSE.

      !Set receiveFromDomain and sendGlobal
      !If global number is on multiple domains
      IF(domainMapping%globalToLocalMap(globalNumber)%numberOfDomains>1) THEN

        !Check if we have a special case where the global number is owned by all domains e.g., as in a constant field
        IF(domainMapping%globalToLocalMap(globalNumber)%numberOfDomains==numberOfGroupComputationNodes) THEN
          ownedByAll=.TRUE.

          !Check if the element is an internal element on all domains
          DO domainIdx=1,domainMapping%globalToLocalMap(globalNumber)%numberOfDomains
            localType=domainMapping%globalToLocalMap(globalNumber)%localType(domainIdx)
            ownedByAll=ownedByAll.AND.localType==DOMAIN_LOCAL_INTERNAL
          ENDDO !domainIdx
        ELSE
          ownedByAll=.FALSE.
        ENDIF

        !If the element is not owned by all domains
        IF(.NOT.ownedByAll) THEN
          receiveFromDomain=-1

          !Loop over the domains where the current number is present
          DO domainIdx=1,domainMapping%globalToLocalMap(globalNumber)%numberOfDomains
            domainNumber=domainMapping%globalToLocalMap(globalNumber)%domainNumber(domainIdx)
            localType=domainMapping%globalToLocalMap(globalNumber)%localType(domainIdx)

            IF(localType/=DOMAIN_LOCAL_GHOST) THEN
              IF(domainNumber==myGroupComputationNodeNumber) sendGlobal=.TRUE.
              IF(receiveFromDomain==-1) THEN
                receiveFromDomain=domainNumber
              ELSE
                localError="Invalid domain mapping. Global number "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
                  & " is owned by domain number "//TRIM(NumberToVString(receiveFromDomain,"*",err,error))// &
                  & " as well as domain number "//TRIM(NumberToVString(domainNumber,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ENDIF
          ENDDO !domainIdx

          IF(receiveFromDomain==-1) THEN
            localError="Invalid domain mapping. Global number "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
              & " is not owned by any domain."
            CALL FlagError(localError,err,error,*999)          
          ENDIF
        ENDIF
      ENDIF

      !Loop over domains of current global number
      DO domainIdx=1,domainMapping%globalToLocalMap(globalNumber)%numberOfDomains
        domainNumber=domainMapping%globalToLocalMap(globalNumber)%domainNumber(domainIdx)
        localNumber=domainMapping%globalToLocalMap(globalNumber)%localNumber(domainIdx)
        localType=domainMapping%globalToLocalMap(globalNumber)%localType(domainIdx)

        IF(domainNumber==myGroupComputationNodeNumber) THEN
          !Set local number
          domainMapping%localToGlobalMap(localNumber)=globalNumber

          !Set entry of current element in domain_list
          SELECT CASE(domainMapping%globalToLocalMap(globalNumber)%localType(domainIdx))
          CASE(DOMAIN_LOCAL_INTERNAL)
            numberOfInternal=numberOfInternal+1
            domainMapping%domainList(numberOfInternal)=localNumber
          CASE(DOMAIN_LOCAL_BOUNDARY)
            numberOfBoundary=numberOfBoundary+1
            domainMapping%domainList(domainMapping%internalFinish+numberOfBoundary)=localNumber
          CASE(DOMAIN_LOCAL_GHOST)
            numberOfGhost=numberOfGhost+1
            domainMapping%domainList(domainMapping%boundaryFinish+numberOfGhost)=localNumber

            !Add local number of ghost element to receive list of domain from which to receive
            CALL List_ItemAdd(ghostReceiveLists(adjacentDomainMap(receiveFromDomain))%ptr,localNumber,err,error,*999)              
          CASE DEFAULT
            localError="The domain local type of "//TRIM(NumberToVString(domainMapping%globalToLocalMap( &
              & globalNumber)%localType(domainIdx),"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT

        ELSE IF(sendGlobal.AND.localType==DOMAIN_LOCAL_GHOST) THEN
          localNumber2=domainMapping%globalToLocalMap(globalNumber)%localNumber(1) !The local number for this node
          CALL List_ItemAdd(ghostSendLists(adjacentDomainMap(domainNumber))%ptr,localNumber2,err,error,*999)            
        ENDIF
      ENDDO !domainIdx
    ENDDO !globalNumber

    !Loop over adjacent domains
    DO domainIdx=1,domainMapping%numberOfAdjacentDomains

      ! transfer the ghost_send_list for adjacent domain to localGhostSendIndices
      CALL List_RemoveDuplicates(ghostSendLists(domainIdx)%ptr,err,error,*999)
      CALL List_DetachAndDestroy(ghostSendLists(domainIdx)%ptr,numberOfGhostSend,sendList,err,error,*999)

      ALLOCATE(domainMapping%adjacentDomains(domainIdx)%localGhostSendIndices(numberOfGhostSend),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate local ghost send inidices.",err,error,*999)

      domainMapping%adjacentDomains(domainIdx)%localGhostSendIndices(1:numberOfGhostSend)=sendList(1:numberOfGhostSend)
      domainMapping%adjacentDomains(domainIdx)%numberOfSendGhosts=numberOfGhostSend
      
      DEALLOCATE(sendList)

      !Transfer the ghost_receive_lists for the current adjacent domain to localGhostReceiveIndices
      CALL List_RemoveDuplicates(ghostReceiveLists(domainIdx)%ptr,err,error,*999)
      CALL List_DetachAndDestroy(ghostReceiveLists(domainIdx)%ptr,numberOfGhostReceive,receiveList,err,error,*999)

      ALLOCATE(domainMapping%adjacentDomains(domainIdx)%localGhostReceiveIndices(numberOfGhostReceive),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate local ghost receive inidices.",err,error,*999)

      domainMapping%adjacentDomains(domainIdx)%localGhostReceiveIndices(1:numberOfGhostReceive)=receiveList(1:numberOfGhostReceive)
      domainMapping%adjacentDomains(domainIdx)%numberOfReceiveGhosts=numberOfGhostReceive
      
      DEALLOCATE(receiveList)

    ENDDO !domainIdx

    DEALLOCATE(adjacentDomainMap)
    DEALLOCATE(ghostSendLists)
    DEALLOCATE(ghostReceiveLists)

    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Domain mappings:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of domains  = ",numberOfGroupComputationNodes,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of global = ",domainMapping%numberOfGlobal,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of local = ",domainMapping%numberOfLocal,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of local = ",domainMapping%totalNumberOfLocal,err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Domain numbers:",err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfGroupComputationNodes,8,8,domainMapping% &
        & numberOfDomainLocal,'("    Number of domain local :",8(X,I10))','(26X,8(X,I10))',err,error,*999)      
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfGroupComputationNodes,8,8,domainMapping% &
        & numberOfDomainGhost,'("    Number of domain ghost :",8(X,I10))','(26X,8(X,I10))',err,error,*999)      
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Domain list:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of internal = ",domainMapping%numberOfInternal,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary = ",domainMapping%numberOfBoundary,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost = ",domainMapping%numberOfGhost,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Internal start = ",domainMapping%internalStart,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Internal finish = ",domainMapping%internalFinish,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Boundary start = ",domainMapping%boundaryStart,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Boundary finish = ",domainMapping%boundaryFinish,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Ghost start = ",domainMapping%ghostStart,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Ghost finish = ",domainMapping%ghostFinish,err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,domainMapping%internalStart,1,domainMapping%internalFinish,8,8, &
        & domainMapping%domainList,'("    Internal list :",8(X,I10))','(19X,8(X,I10))',err,error,*999)      
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,domainMapping%boundaryStart,1,domainMapping%boundaryFinish,8,8, &
        & domainMapping%domainList,'("    Boundary list :",8(X,I10))','(19X,8(X,I10))',err,error,*999)      
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,domainMapping%ghostStart,1,domainMapping%ghostFinish,8,8, &
        & domainMapping%domainList,'("    Ghost list    :",8(X,I10))','(19X,8(X,I10))',err,error,*999)      
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map:",err,error,*999)
      DO idx=1,domainMapping%totalNumberOfLocal
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Local index : ",idx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Global index = ",domainMapping%localToGlobalMap(idx),err,error,*999)
      ENDDO !idx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Global to local map:",err,error,*999)
      DO idx=1,domainMapping%numberOfGlobal
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global idx : ",idx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of domains  = ", &
          & domainMapping%globalToLocalMap(idx)%numberOfDomains,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,domainMapping%globalToLocalMap(idx)% &
          & numberOfDomains,8,8,domainMapping%globalToLocalMap(idx)%localNumber, &
          & '("      Local number  :",8(X,I10))','(21X,8(X,I10))',err,error,*999)      
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,domainMapping%globalToLocalMap(idx)% &
          & numberOfDomains,8,8,domainMapping%globalToLocalMap(idx)%domainNumber, &
          & '("      Domain number :",8(X,I10))','(21X,8(X,I10))',err,error,*999)      
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,domainMapping%globalToLocalMap(idx)% &
          & numberOfDomains,8,8,domainMapping%globalToLocalMap(IDX)%localType, &
          & '("      Local type    :",8(X,I10))','(21X,8(X,I10))',err,error,*999)      
      ENDDO !ne
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
        & domainMapping%numberOfAdjacentDomains,err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfGroupComputationNodes+1,8,8, &
        & domainMapping%adjacentDomainsPtr,'("    Adjacent domains ptr  :",8(X,I5))','(27X,8(X,I5))',err,error,*999)
      IF(domainMapping%numberOfAdjacentDomains>0) THEN
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,domainMapping%adjacentDomainsPtr( &
          & numberOfGroupComputationNodes)-1,8,8,domainMapping%adjacentDomainsList, &
          '("    Adjacent domains list :",8(X,I5))','(27X,8(X,I5))',err,error,*999)
        DO domainIdx=1,domainMapping%numberOfAdjacentDomains
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",domainIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
            & domainMapping%adjacentDomains(domainIdx)%domainNumber,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
            & domainMapping%adjacentDomains(domainIdx)%numberOfSendGhosts,err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,domainMapping%adjacentDomains(domainIdx)% &
            & numberOfSendGhosts,8,8,domainMapping%adjacentDomains(domainIdx)%localGhostSendIndices, &
            & '("      Local send ghost indices        :",8(X,I10))','(39X,8(X,I10))',err,error,*999)      
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of receive ghosts = ", &
            & domainMapping%adjacentDomains(domainIdx)%numberOfReceiveGhosts,err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,domainMapping%adjacentDomains(domainIdx)% &
            & numberOfReceiveGhosts,8,8,domainMapping%adjacentDomains(domainIdx)%localGhostReceiveIndices, &
            & '("      Local receive ghost indices     :",8(X,I10))','(39X,8(X,I10))',err,error,*999)              
        ENDDO !domainIdx
      ENDIF
    ENDIF
    
    EXITS("DomainMapping_LocalFromGlobalCalculate")
    RETURN
999 IF(ALLOCATED(sendList)) DEALLOCATE(sendList)
    IF(ALLOCATED(receiveList)) DEALLOCATE(receiveList)
    IF(ALLOCATED(adjacentDomainMap)) DEALLOCATE(adjacentDomainMap)
    IF(ALLOCATED(adjacentDomains)) DEALLOCATE(adjacentDomains)
    IF(ALLOCATED(ghostSendLists)) THEN
      DO domainIdx=1,SIZE(ghostSendLists)
        IF(ASSOCIATED(ghostSendLists(domainIdx)%ptr)) &
          & CALL LIST_DESTROY(ghostSendLists(domainIdx)%ptr,dummyErr,dummyError,*998)
      ENDDO !domainIdx
998   DEALLOCATE(ghostSendLists)
    ENDIF
    IF(ALLOCATED(ghostReceiveLists)) THEN
      DO domainIdx=1,SIZE(ghostReceiveLists)
        IF(ASSOCIATED(ghostReceiveLists(domainIdx)%ptr)) &
          & CALL LIST_DESTROY(ghostReceiveLists(domainIdx)%ptr,dummyErr,dummyError,*997)
      ENDDO !domainIdx
997   DEALLOCATE(ghostReceiveLists)
    ENDIF
    ERRORSEXITS("DomainMapping_LocalFromGlobalCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMapping_LocalFromGlobalCalculate
  
  !
  !================================================================================================================================
  !

  !>Gets the local number from a global number in a domain in a domain mapping.
  SUBROUTINE DomainMapping_LocalNumberFromGlobalGet(domainMapping,globalNumber,domainIdx,localNumber,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainMapping !<A pointer to the domain mapping to get the local number for
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number to get the local number for
    INTEGER(INTG), INTENT(IN) :: domainIdx !<The domain index to get the local number for
    INTEGER(INTG), INTENT(OUT) :: localNumber !<On exit, the local number for the global number and domain index the domain mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainMapping_LocalNumberFromGlobalGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(domainMapping%globalToLocalMap)) &
      & CALL FlagError("The global to local map is not allocated in the domain mapping.",err,error,*999)
    IF(globalNumber<1.OR.globalNumber>domainMapping%numberOfGlobal) THEN
      localError="The specified global number of "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
        & " is invalid for the domain mapping. The global number should be >= 1 and <= "// &
        & TRIM(NumberToVString(domainMapping%numberOfGlobal,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(domainIdx<1.OR.domainIdx>domainMapping%globalToLocalMap(globalNumber)%numberOfDomains) THEN
      localError="The specified domain index of "//TRIM(NumberToVString(domainIdx,"*",err,error))// &
        & " is invalid for global number "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
        & " in the domain mapping. The domain index should be >= 1 and <= "// &
        & TRIM(NumberToVString(domainMapping%globalToLocalMap(globalNumber)%numberOfDomains,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    localNumber=domainMapping%globalToLocalMap(globalNumber)%localNumber(domainIdx)
    
    EXITS("DomainMapping_LocalNumberFromGlobalGet")
    RETURN
999 ERRORSEXITS("DomainMapping_LocalNumberFromGlobalGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMapping_LocalNumberFromGlobalGet
  
  !
  !================================================================================================================================
  !

  !>Gets the local number type from a global number in a domain in a domain mapping.
  SUBROUTINE DomainMapping_LocalTypeFromGlobalGet(domainMapping,globalNumber,domainIdx,localType,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainMapping !<A pointer to the domain mapping to get the local type for
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number to get the local type for
    INTEGER(INTG), INTENT(IN) :: domainIdx !<The domain index to get the local type for
    INTEGER(INTG), INTENT(OUT) :: localType !<On exit, the local type for the global number and domain index the domain mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainMapping_LocalTypeFromGlobalGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(domainMapping%globalToLocalMap)) &
      & CALL FlagError("The global to local map is not allocated in the domain mapping.",err,error,*999)
    IF(globalNumber<1.OR.globalNumber>domainMapping%numberOfGlobal) THEN
      localError="The specified global number of "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
        & " is invalid for the domain mapping. The global number should be >= 1 and <= "// &
        & TRIM(NumberToVString(domainMapping%numberOfGlobal,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(domainIdx<1.OR.domainIdx>domainMapping%globalToLocalMap(globalNumber)%numberOfDomains) THEN
      localError="The specified domain index of "//TRIM(NumberToVString(domainIdx,"*",err,error))// &
        & " is invalid for global number "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
        & " in the domain mapping. The domain index should be >= 1 and <= "// &
        & TRIM(NumberToVString(domainMapping%globalToLocalMap(globalNumber)%numberOfDomains,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    localType=domainMapping%globalToLocalMap(globalNumber)%localType(domainIdx)
    
    EXITS("DomainMapping_LocalTypeFromGlobalGet")
    RETURN
999 ERRORSEXITS("DomainMapping_LocalTypeFromGlobalGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMapping_LocalTypeFromGlobalGet
  
  !
  !================================================================================================================================
  !

  !>Gets the number in a domain mapping.
  SUBROUTINE DomainMapping_NumberGet(domainMapping,mappingIdx,mappingNumber,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainMapping !<A pointer to the domain mapping to get the number for
    INTEGER(INTG), INTENT(IN) :: mappingIdx !<The index in the mapping to get the mapping number for
    INTEGER(INTG), INTENT(OUT) :: mappingNumber !<On exit, the mapping number for the mapping index for the domain mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainMapping_NumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*999)
    IF(mappingIdx<1.OR.mappingIdx>domainMapping%totalNumberOfLocal) THEN
      localError="The specified mapping index of "//TRIM(NumberToVString(mappingIdx,"*",err,error))// &
        & " is invalid. The mapping index should be >= 1 and <= "// &
        & TRIM(NumberToVString(domainMapping%totalNumberOfLocal,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainMapping%domainList)) &
      & CALL FlagError("The domain list is not allocated for the domain mapping.",err,error,*999)
#endif    

    mappingNumber=domainMapping%domainList(mappingIdx)
    
    EXITS("DomainMapping_NumberGet")
    RETURN
999 ERRORSEXITS("DomainMapping_NumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMapping_NumberGet
  
  !
  !================================================================================================================================
  !

  !>Gets the number of boundary for a domain mapping.
  SUBROUTINE DomainMapping_NumberOfBoundaryGet(domainMapping,numberOfBoundary,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainMapping !<A pointer to the domain mapping to get the number of boundary for
    INTEGER(INTG), INTENT(OUT) :: numberOfBoundary !<On exit, the number of boundary for the domain mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainMapping_NumberOfBoundaryGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*999)
#endif    

    numberOfBoundary=domainMapping%numberOfBoundary
    
    EXITS("DomainMapping_NumberOfBoundaryGet")
    RETURN
999 ERRORSEXITS("DomainMapping_NumberOfBoundaryGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMapping_NumberOfBoundaryGet
  
  !
  !================================================================================================================================
  !

  !>Gets the number of domains that a global number in a domain mapping is in.
  SUBROUTINE DomainMapping_NumberOfDomainsFromGlobalGet(domainMapping,globalNumber,numberOfDomains,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainMapping !<A pointer to the domain mapping to get the number of domains for
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number to get the number of domains for
    INTEGER(INTG), INTENT(OUT) :: numberOfDomains !<On exit, the number of domains for the global number the domain mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("DomainMapping_NumberOfDomainsFromGlobalGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(domainMapping%globalToLocalMap)) &
      & CALL FlagError("The global to local map is not allocated in the domain mapping.",err,error,*999)
    IF(globalNumber<1.OR.globalNumber>domainMapping%numberOfGlobal) THEN
      localError="The specified global number of "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
        & " is invalid for the domain mapping. The global number should be >= 1 and <= "// &
        & TRIM(NumberToVString(domainMapping%numberOfGlobal,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    numberOfDomains=domainMapping%globalToLocalMap(globalNumber)%numberOfDomains
    
    EXITS("DomainMapping_NumberOfDomainsFromGlobalGet")
    RETURN
999 ERRORSEXITS("DomainMapping_NumberOfDomainsFromGlobalGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMapping_NumberOfDomainsFromGlobalGet
  
  !
  !================================================================================================================================
  !

  !>Gets the number of internal for a domain mapping.
  SUBROUTINE DomainMapping_NumberOfInternalGet(domainMapping,numberOfInternal,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainMapping !<A pointer to the domain mapping to get the number of internal for
    INTEGER(INTG), INTENT(OUT) :: numberOfInternal !<On exit, the number of internal for the domain mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainMapping_NumberOfInternalGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*999)
#endif    

    numberOfInternal=domainMapping%numberOfInternal
    
    EXITS("DomainMapping_NumberOfInternalGet")
    RETURN
999 ERRORSEXITS("DomainMapping_NumberOfInternalGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMapping_NumberOfInternalGet
  
  !
  !================================================================================================================================
  !

  !>Gets the number of ghost for a domain mapping.
  SUBROUTINE DomainMapping_NumberOfGhostGet(domainMapping,numberOfGhost,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainMapping !<A pointer to the domain mapping to get the number of ghost for
    INTEGER(INTG), INTENT(OUT) :: numberOfGhost !<On exit, the number of ghost for the domain mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainMapping_NumberOfGhostGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*999)
#endif    

    numberOfGhost=domainMapping%numberOfGhost
    
    EXITS("DomainMapping_NumberOfGhostGet")
    RETURN
999 ERRORSEXITS("DomainMapping_NumberOfGhostGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMapping_NumberOfGhostGet
  
  !
  !================================================================================================================================
  !

  !>Gets the number of global for a domain mapping.
  SUBROUTINE DomainMapping_NumberOfGlobalGet(domainMapping,numberOfGlobal,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainMapping !<A pointer to the domain mapping to get the number of global for
    INTEGER(INTG), INTENT(OUT) :: numberOfGlobal !<On exit, the number of global for the domain mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainMapping_NumberOfGlobalGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*999)
#endif    

    numberOfGlobal=domainMapping%numberOfGlobal
    
    EXITS("DomainMapping_NumberOfGlobalGet")
    RETURN
999 ERRORSEXITS("DomainMapping_NumberOfGlobalGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMapping_NumberOfGlobalGet
  
  !
  !================================================================================================================================
  !

  !>Gets the number of local for a domain mapping.
  SUBROUTINE DomainMapping_NumberOfLocalGet(domainMapping,numberOfLocal,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainMapping !<A pointer to the domain mapping to get the number of local for
    INTEGER(INTG), INTENT(OUT) :: numberOfLocal !<On exit, the number of local for the domain mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainMapping_NumberOfLocalGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*999)
#endif    

    numberOfLocal=domainMapping%numberOfLocal
    
    EXITS("DomainMapping_NumberOfLocalGet")
    RETURN
999 ERRORSEXITS("DomainMapping_NumberOfLocalGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMapping_NumberOfLocalGet
  
  !
  !================================================================================================================================
  !

  !>Gets the total number of local for a domain mapping.
  SUBROUTINE DomainMapping_TotalNumberOfLocalGet(domainMapping,totalNumberOfLocal,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainMapping !<A pointer to the domain mapping to get the total number of local for
    INTEGER(INTG), INTENT(OUT) :: totalNumberOfLocal !<On exit, the total number of local for the domain mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainMapping_TotalNumberOfLocalGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*999)
#endif    

    totalNumberOfLocal=domainMapping%totalNumberOfLocal
    
    EXITS("DomainMapping_TotalNumberOfLocalGet")
    RETURN
999 ERRORSEXITS("DomainMapping_TotalNumberOfLocalGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMapping_TotalNumberOfLocalGet
  
  !
  !================================================================================================================================
  !

  !>Gets the work group for a domain mapping.
  SUBROUTINE DomainMapping_WorkGroupGet(domainMapping,workGroup,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainMapping !<A pointer to the domain mapping to get the work group for
    TYPE(WorkGroupType), POINTER :: workGroup !<On return, a pointer to the work group for the domain mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainMapping_WorkGroupGet",err,error,*998)

    IF(ASSOCIATED(workGroup)) CALL FlagError("Work group is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*999)

    workGroup=>domainMapping%workGroup
    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("The work group for the domain mapping is not associated.",err,error,*999)
    
    EXITS("DomainMapping_WorkGroupGet")
    RETURN
999 NULLIFY(workGroup)
998 ERRORSEXITS("DomainMapping_WorkGroupGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMapping_WorkGroupGet
  
  !
  !================================================================================================================================
  !

  !>Sets the work group for a domain mapping.
  SUBROUTINE DomainMapping_WorkGroupSet(domainMapping,workGroup,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainMapping !<A pointer to the domain mapping to set the work group for
    TYPE(WorkGroupType), POINTER :: workGroup !<A pointer to the work group to set for the domain mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainMapping_WorkGroupSet",err,error,*999)

    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("Work group is not associated.",err,error,*999)

    domainMapping%workGroup=>workGroup
    
    EXITS("DomainMapping_WorkGroupSet")
    RETURN
999 ERRORSEXITS("DomainMapping_WorkGroupSet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMapping_WorkGroupSet
  
  !
  !================================================================================================================================
  !
 
END MODULE DomainMappings
