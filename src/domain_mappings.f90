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
MODULE DOMAIN_MAPPINGS

  USE BaseRoutines
  USE ComputationRoutines
  USE ComputationAccessRoutines
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE LISTS
  USE STRINGS
  USE TYPES

#include "macros.h"  


  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup DOMAIN_MAPPINGS_DomainType DOMAIN_MAPPINGS::DomainType
  !> \see DOMAIN_MAPPINGS
  !>@{
  INTEGER(INTG), PARAMETER :: DOMAIN_LOCAL_INTERNAL=1 !<The domain item is internal to the domain \see DOMAIN_MAPPINGS_DomainType,DOMAIN_MAPPINGS
  INTEGER(INTG), PARAMETER :: DOMAIN_LOCAL_BOUNDARY=2 !<The domain item is on the boundary of the domain \see DOMAIN_MAPPINGS_DomainType,DOMAIN_MAPPINGS
  INTEGER(INTG), PARAMETER :: DOMAIN_LOCAL_GHOST=3 !<The domain item is ghosted from another domain \see DOMAIN_MAPPINGS_DomainType,DOMAIN_MAPPINGS
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

  PUBLIC DOMAIN_LOCAL_INTERNAL,DOMAIN_LOCAL_BOUNDARY,DOMAIN_LOCAL_GHOST
  
  PUBLIC DOMAIN_MAPPINGS_MAPPING_FINALISE,DomainMappings_MappingInitialise,DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE, &
    & DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET,DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE

CONTAINS
  
  !
  !================================================================================================================================
  !

  !>Finalises the adjacent domain and deallocates all memory for a domain mapping.
  SUBROUTINE DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE(ADJACENT_DOMAIN,err,error,*)

    !Argument variables
    TYPE(DomainAdjacentDomainType) :: ADJACENT_DOMAIN !<The adjacent domain to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE",err,error,*999)

    IF(ALLOCATED(ADJACENT_DOMAIN%localGhostSendIndices)) DEALLOCATE(ADJACENT_DOMAIN%localGhostSendIndices)
    IF(ALLOCATED(ADJACENT_DOMAIN%localGhostReceiveIndices)) DEALLOCATE(ADJACENT_DOMAIN%localGhostReceiveIndices)
    ADJACENT_DOMAIN%numberOfSendGhosts=0
    ADJACENT_DOMAIN%numberOfReceiveGhosts=0
    ADJACENT_DOMAIN%domainNumber=0
    
    EXITS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE",err,error)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialise the adjacent domain for a domain mapping.
  SUBROUTINE DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE(ADJACENT_DOMAIN,err,error,*)

    !Argument variables
    TYPE(DomainAdjacentDomainType) :: ADJACENT_DOMAIN !<The adjacent domain to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE",err,error,*999)

    ADJACENT_DOMAIN%numberOfSendGhosts=0
    ADJACENT_DOMAIN%numberOfReceiveGhosts=0
    ADJACENT_DOMAIN%domainNumber=0
    
    EXITS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE",err,error)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Returns the local number, if it exists on the rank, for the specifed global number
  SUBROUTINE DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(DOMAIN_MAPPING,GLOBAL_NUMBER,LOCAL_EXISTS,localNumber,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: DOMAIN_MAPPING !<A pointer to the domain mapping to get the local number from
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to get the local number for
    LOGICAL, INTENT(OUT) :: LOCAL_EXISTS !<On exit, is .TRUE. if the specifed global number exists on the local rank, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: localNumber !<On exit, the local number corresponding to the global number if it exists. If it doesn't exist then 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: myGroupComputationNodeNumber
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET",err,error,*999)

    LOCAL_EXISTS=.FALSE.
    localNumber=0
    IF(ASSOCIATED(DOMAIN_MAPPING)) THEN
      IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=DOMAIN_MAPPING%numberOfGlobal) THEN
        CALL WorkGroup_GroupNodeNumberGet(DOMAIN_MAPPING%workGroup,myGroupComputationNodeNumber,err,error,*999)
        IF(DOMAIN_MAPPING%globalToLocalMap(GLOBAL_NUMBER)%domainNumber(1)==myGroupComputationNodeNumber) THEN
          localNumber=DOMAIN_MAPPING%globalToLocalMap(GLOBAL_NUMBER)%localNumber(1)
          LOCAL_EXISTS=.TRUE.
        ENDIF
      ELSE
        LOCAL_ERROR="The specified global number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",err,error))// &
          & " is invalid. The number must be between 1 and "// &
          & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING%numberOfGlobal,"*",err,error))//"."
        CALL FlagError(LOCAL_ERROR,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET",err,error)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET
  
  !
  !================================================================================================================================
  !

  !>Calculates the domain mappings local map from a domain mappings global map.
  SUBROUTINE DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(DOMAIN_MAPPING,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: DOMAIN_MAPPING !<The domain mapping to calculate the local mappings
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: domain_idx,domain_idx2,domain_no,domain_no2,global_number,idx,local_number,local_number2,NUMBER_INTERNAL, &
      & NUMBER_BOUNDARY,NUMBER_GHOST,myGroupComputationNodeNumber,MY_DOMAIN_INDEX,TEMP,numberOfAdjacentDomains, &
      & RECEIVE_FROM_DOMAIN,DUMMY_ERR,NUMBER_OF_GHOST_RECEIVE,NUMBER_OF_GHOST_SEND,local_type,COUNT, &
      & TOTAL_numberOfAdjacentDomains
    INTEGER(INTG), ALLOCATABLE :: ADJACENT_DOMAIN_MAP(:),adjacentDomains(:,:),SEND_LIST(:),RECEIVE_LIST(:)
    LOGICAL :: OWNED_BY_ALL,SEND_GLOBAL
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: GHOST_SEND_LISTS(:),GHOST_RECEIVE_LISTS(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR,DUMMY_ERROR
    
    ENTERS("DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE",err,error,*999)

    IF(ASSOCIATED(DOMAIN_MAPPING)) THEN
      CALL WorkGroup_GroupNodeNumberGet(DOMAIN_MAPPING%workGroup,myGroupComputationNodeNumber,err,error,*999)
      IF(ERR/=0) GOTO 999        
      
      !Calculate local to global maps from global to local map
      ALLOCATE(DOMAIN_MAPPING%numberOfDomainLocal(0:DOMAIN_MAPPING%numberOfDomains-1),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate number of domain local.",err,error,*999)
      DOMAIN_MAPPING%numberOfDomainLocal=0
      
      ALLOCATE(DOMAIN_MAPPING%numberOfDomainGhost(0:DOMAIN_MAPPING%numberOfDomains-1),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate number of domain ghost.",err,error,*999)
      
      DOMAIN_MAPPING%numberOfDomainGhost=0
      NUMBER_INTERNAL=0   ! counters for my computational node
      NUMBER_BOUNDARY=0
      NUMBER_GHOST=0
      
      ALLOCATE(adjacentDomains(0:DOMAIN_MAPPING%numberOfDomains-1,0:DOMAIN_MAPPING%numberOfDomains-1),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate adjacent domains.",err,error,*999)
      adjacentDomains=0
      
      ! loop over global elements
      DO global_number=1,DOMAIN_MAPPING%numberOfGlobal
        !-------- If necessary, reset global domain index so that my computational node is in the first index position ------------
        !find out current domain index of my computational node
        MY_DOMAIN_INDEX=1
        DO domain_idx=2,DOMAIN_MAPPING%globalToLocalMap(global_number)%numberOfDomains
          domain_no=DOMAIN_MAPPING%globalToLocalMap(global_number)%domainNumber(domain_idx)
          IF(domain_no==myGroupComputationNodeNumber) THEN
            MY_DOMAIN_INDEX=domain_idx
            EXIT
          ENDIF
        ENDDO !domain_idx
        
        !if necessary, swap the data structure at the current domain index for my computational node with domain index 1
        IF(MY_DOMAIN_INDEX/=1) THEN
          !Swap domain index in the global to local map, 1 <-> MY_DOMAIN_INDEX
          TEMP=DOMAIN_MAPPING%globalToLocalMap(global_number)%localNumber(1)
          DOMAIN_MAPPING%globalToLocalMap(global_number)%localNumber(1) = &
            & DOMAIN_MAPPING%globalToLocalMap(global_number)%localNumber(MY_DOMAIN_INDEX)
          DOMAIN_MAPPING%globalToLocalMap(global_number)%localNumber(MY_DOMAIN_INDEX) = TEMP
          
          TEMP=DOMAIN_MAPPING%globalToLocalMap(global_number)%domainNumber(1)
          DOMAIN_MAPPING%globalToLocalMap(global_number)%domainNumber(1) = &
            & DOMAIN_MAPPING%globalToLocalMap(global_number)%domainNumber(MY_DOMAIN_INDEX)
          DOMAIN_MAPPING%globalToLocalMap(global_number)%domainNumber(MY_DOMAIN_INDEX) = TEMP
          
          TEMP=DOMAIN_MAPPING%globalToLocalMap(global_number)%localType(1)
          DOMAIN_MAPPING%globalToLocalMap(global_number)%localType(1) = &
            & DOMAIN_MAPPING%globalToLocalMap(global_number)%localType(MY_DOMAIN_INDEX)
          DOMAIN_MAPPING%globalToLocalMap(global_number)%localType(MY_DOMAIN_INDEX) = TEMP
        ENDIF
        
        !set the global adjacent_domains array to 1 for domains that share the current element
        DO domain_idx=1,DOMAIN_MAPPING%globalToLocalMap(global_number)%numberOfDomains
          domain_no=DOMAIN_MAPPING%globalToLocalMap(global_number)%domainNumber(domain_idx)
          DO domain_idx2=1,DOMAIN_MAPPING%globalToLocalMap(global_number)%numberOfDomains
            domain_no2=DOMAIN_MAPPING%globalToLocalMap(global_number)%domainNumber(domain_idx2)
            adjacentDomains(domain_no,domain_no2)=1
          ENDDO !domain_idx2
        ENDDO !domain_idx
        
        !loop over domains where current element is
        DO domain_idx=1,DOMAIN_MAPPING%globalToLocalMap(global_number)%numberOfDomains
          domain_no=DOMAIN_MAPPING%globalToLocalMap(global_number)%domainNumber(domain_idx)
          local_type=DOMAIN_MAPPING%globalToLocalMap(global_number)%localType(domain_idx)
          !increment counter of total ghost/local elements per domain
          IF(local_type==DOMAIN_LOCAL_GHOST) THEN
            DOMAIN_MAPPING%numberOfDomainGhost(domain_no)=DOMAIN_MAPPING%numberOfDomainGhost(domain_no)+1
          ELSE
            DOMAIN_MAPPING%numberOfDomainLocal(domain_no)=DOMAIN_MAPPING%numberOfDomainLocal(domain_no)+1
          ENDIF
          
          !increment counter of internal, boundary and ghost elements on my domain
          IF(domain_no==myGroupComputationNodeNumber) THEN
            SELECT CASE(local_type)
            CASE(DOMAIN_LOCAL_INTERNAL)
              NUMBER_INTERNAL=NUMBER_INTERNAL+1
            CASE(DOMAIN_LOCAL_BOUNDARY)
              NUMBER_BOUNDARY=NUMBER_BOUNDARY+1
            CASE(DOMAIN_LOCAL_GHOST)
              NUMBER_GHOST=NUMBER_GHOST+1
            CASE DEFAULT
              LOCAL_ERROR="The domain local type of "//TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING%globalToLocalMap( &
                & global_number)%localType(domain_idx),"*",err,error))//" is invalid."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          ENDIF
        ENDDO !domain_idx
      ENDDO !global_number

      !!TODO: move adjacent domains calculation back to where the global to local array is set up????
      !count number of adjacent domains for my domain and in total
      numberOfAdjacentDomains=0
      TOTAL_numberOfAdjacentDomains=0
      DO domain_no=0,DOMAIN_MAPPING%numberOfDomains-1
        DO domain_no2=0,DOMAIN_MAPPING%numberOfDomains-1
          IF(domain_no/=domain_no2) THEN 
            IF(adjacentDomains(domain_no,domain_no2)>0) THEN
              TOTAL_numberOfAdjacentDomains=TOTAL_numberOfAdjacentDomains+1
              IF(domain_no==myGroupComputationNodeNumber) numberOfAdjacentDomains=numberOfAdjacentDomains+1
            ENDIF
          ENDIF
        ENDDO !domain_no2
      ENDDO !domain_no
      
      ALLOCATE(DOMAIN_MAPPING%adjacentDomainsPtr(0:DOMAIN_MAPPING%numberOfDomains),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate adjacent domains ptr.",err,error,*999)
      ALLOCATE(DOMAIN_MAPPING%adjacentDomainsList(TOTAL_numberOfAdjacentDomains),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate adjacent domains list.",err,error,*999)
      
      ! store the adjacent domains for a domain in adjacentDomainsList starting at index COUNT and store the starting index COUNT in adjacentDomainsPtr
      ! the adjacent domains for a domain are then adjacentDomainsList(adjacentDomainsPtr(domain) : adjacentDomainsPtr(domain+1))
      COUNT=1
      DO domain_no=0,DOMAIN_MAPPING%numberOfDomains-1
        DOMAIN_MAPPING%adjacentDomainsPtr(domain_no)=COUNT
        DO domain_no2=0,DOMAIN_MAPPING%numberOfDomains-1
          IF(domain_no/=domain_no2) THEN
            IF(adjacentDomains(domain_no,domain_no2)>0) THEN
              DOMAIN_MAPPING%adjacentDomainsList(COUNT)=domain_no2
              COUNT=COUNT+1
            ENDIF
          ENDIF
        ENDDO !domain_no2
      ENDDO !domain_no
      
      DOMAIN_MAPPING%adjacentDomainsPtr(DOMAIN_MAPPING%numberOfDomains)=COUNT
      DEALLOCATE(adjacentDomains)
      
      !compute domain list for my computational node
      ALLOCATE(DOMAIN_MAPPING%domainList(NUMBER_INTERNAL+NUMBER_BOUNDARY+NUMBER_GHOST),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate domain map domain list.",err,error,*999)
      ALLOCATE(DOMAIN_MAPPING%localToGlobalMap(NUMBER_INTERNAL+NUMBER_BOUNDARY+NUMBER_GHOST),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate domain map local to global list.",err,error,*999)
      
      !set constants
      DOMAIN_MAPPING%totalNumberOfLocal=NUMBER_INTERNAL+NUMBER_BOUNDARY+NUMBER_GHOST
      DOMAIN_MAPPING%numberOfLocal=NUMBER_INTERNAL+NUMBER_BOUNDARY
      DOMAIN_MAPPING%numberOfInternal=NUMBER_INTERNAL
      DOMAIN_MAPPING%numberOfBoundary=NUMBER_BOUNDARY
      DOMAIN_MAPPING%numberOfGhost=NUMBER_GHOST
      DOMAIN_MAPPING%internalStart=1
      DOMAIN_MAPPING%internalFinish=NUMBER_INTERNAL
      DOMAIN_MAPPING%boundaryStart=NUMBER_INTERNAL+1
      DOMAIN_MAPPING%boundaryFinish=NUMBER_INTERNAL+NUMBER_BOUNDARY
      DOMAIN_MAPPING%ghostStart=NUMBER_INTERNAL+NUMBER_BOUNDARY+1
      DOMAIN_MAPPING%ghostFinish=NUMBER_INTERNAL+NUMBER_BOUNDARY+NUMBER_GHOST
      
      ! adjacent_domains maps a domain index (index between 1 and DOMAIN_MAPPING%numberOfAdjacentDomains) to the domain number
      ALLOCATE(DOMAIN_MAPPING%adjacentDomains(numberOfAdjacentDomains),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate adjacent domains.",err,error,*999)
      DOMAIN_MAPPING%numberOfAdjacentDomains=numberOfAdjacentDomains
      
      ! adjacent_domain_map maps a domain number back to its index (index between 1 and DOMAIN_MAPPING%numberOfAdjacentDomains)
      ALLOCATE(ADJACENT_DOMAIN_MAP(0:DOMAIN_MAPPING%numberOfDomains-1),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate adjacent domain map.",err,error,*999)
      
      ! ghost_send_lists contains a list for each adjacent domain with elements in my domain (by local element numbers) that will be sent to that foreign domain
      ALLOCATE(GHOST_SEND_LISTS(DOMAIN_MAPPING%numberOfAdjacentDomains),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate ghost send list.",err,error,*999)
      
      ! ghost_receive_lists contains a list for each adjacent domain and contains the local numbers of ghost elements that can be received from that foreign domain
      ALLOCATE(GHOST_RECEIVE_LISTS(DOMAIN_MAPPING%numberOfAdjacentDomains),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate ghost recieve list.",err,error,*999)
      
      ! set adjacent domains data structures and initialize ghost send and receive lists
      ! loop over adjacent domains of my computational node
      DO domain_idx=1,DOMAIN_MAPPING%numberOfAdjacentDomains
        ! set variables to 0
        CALL DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE(DOMAIN_MAPPING%adjacentDomains(domain_idx),err,error,*999)
        
        ! get number of current adjacent domain
        domain_no= &
          & DOMAIN_MAPPING%adjacentDomainsList(DOMAIN_MAPPING%adjacentDomainsPtr(myGroupComputationNodeNumber)+domain_idx-1)
        ! set number in adjacent_domains and adjacent_domain_map
        DOMAIN_MAPPING%adjacentDomains(domain_idx)%domainNumber=domain_no
        ADJACENT_DOMAIN_MAP(domain_no)=domain_idx
        
        ! initialize send and receive lists for ghosts
        NULLIFY(GHOST_SEND_LISTS(domain_idx)%PTR)
        CALL LIST_CREATE_START(GHOST_SEND_LISTS(domain_idx)%PTR,err,error,*999)
        CALL LIST_DATA_TYPE_SET(GHOST_SEND_LISTS(domain_idx)%PTR,LIST_INTG_TYPE,err,error,*999)
        CALL LIST_INITIAL_SIZE_SET(GHOST_SEND_LISTS(domain_idx)%PTR,MAX(DOMAIN_MAPPING%numberOfGhost,1),err,error,*999)
        CALL LIST_CREATE_FINISH(GHOST_SEND_LISTS(domain_idx)%PTR,err,error,*999)
        
        NULLIFY(GHOST_RECEIVE_LISTS(domain_idx)%PTR)
        CALL LIST_CREATE_START(GHOST_RECEIVE_LISTS(domain_idx)%PTR,err,error,*999)
        CALL LIST_DATA_TYPE_SET(GHOST_RECEIVE_LISTS(domain_idx)%PTR,LIST_INTG_TYPE,err,error,*999)
        CALL LIST_INITIAL_SIZE_SET(GHOST_RECEIVE_LISTS(domain_idx)%PTR,MAX(DOMAIN_MAPPING%numberOfGhost,1),err,error,*999)
        CALL LIST_CREATE_FINISH(GHOST_RECEIVE_LISTS(domain_idx)%PTR,err,error,*999)
      ENDDO !domain_idx
      
      NUMBER_INTERNAL=0
      NUMBER_BOUNDARY=0
      NUMBER_GHOST=0
      
      ! loop over global elements
      DO global_number=1,DOMAIN_MAPPING%numberOfGlobal
        SEND_GLOBAL=.FALSE.
        
        ! set RECEIVE_FROM_DOMAIN and SEND_GLOBAL
        ! if element is on multiple domains
        IF(DOMAIN_MAPPING%globalToLocalMap(global_number)%numberOfDomains>1) THEN
          
          !Check if we have a special case where the global number is owned by all domains e.g., as in a constant field
          IF(DOMAIN_MAPPING%globalToLocalMap(global_number)%numberOfDomains==DOMAIN_MAPPING%numberOfDomains) THEN
            OWNED_BY_ALL=.TRUE.
          
            ! check if the element is an internal element on all domains
            DO domain_idx=1,DOMAIN_MAPPING%globalToLocalMap(global_number)%numberOfDomains
              local_type=DOMAIN_MAPPING%globalToLocalMap(global_number)%localType(domain_idx)
              OWNED_BY_ALL=OWNED_BY_ALL.AND.local_type==DOMAIN_LOCAL_INTERNAL
            ENDDO !domain_idx
          ELSE
            OWNED_BY_ALL=.FALSE.
          ENDIF
          
          ! if the element is not owned by all domains
          IF(.NOT.OWNED_BY_ALL) THEN
            RECEIVE_FROM_DOMAIN=-1
            
            ! loop over the domains where the current element is present
            DO domain_idx=1,DOMAIN_MAPPING%globalToLocalMap(global_number)%numberOfDomains
              domain_no=DOMAIN_MAPPING%globalToLocalMap(global_number)%domainNumber(domain_idx)
              local_type=DOMAIN_MAPPING%globalToLocalMap(global_number)%localType(domain_idx)
              
              IF(local_type/=DOMAIN_LOCAL_GHOST) THEN
                IF(domain_no==myGroupComputationNodeNumber) SEND_GLOBAL=.TRUE.
                IF(RECEIVE_FROM_DOMAIN==-1) THEN
                  RECEIVE_FROM_DOMAIN=domain_no
                ELSE
                  LOCAL_ERROR="Invalid domain mapping. Global number "//TRIM(NUMBER_TO_VSTRING(global_number,"*",err,error))// &
                    & " is owned by domain number "//TRIM(NUMBER_TO_VSTRING(RECEIVE_FROM_DOMAIN,"*",err,error))// &
                    & " as well as domain number "//TRIM(NUMBER_TO_VSTRING(domain_no,"*",err,error))//"."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                ENDIF
              ENDIF
            ENDDO !domain_idx
            
            IF(RECEIVE_FROM_DOMAIN==-1) THEN
              LOCAL_ERROR="Invalid domain mapping. Global number "//TRIM(NUMBER_TO_VSTRING(global_number,"*",err,error))// &
                & " is not owned by any domain."
              CALL FlagError(LOCAL_ERROR,err,error,*999)          
            ENDIF
          ENDIF
        ENDIF
        
        ! loop over domains of current element
        DO domain_idx=1,DOMAIN_MAPPING%globalToLocalMap(global_number)%numberOfDomains
          domain_no=DOMAIN_MAPPING%globalToLocalMap(global_number)%domainNumber(domain_idx)
          local_number=DOMAIN_MAPPING%globalToLocalMap(global_number)%localNumber(domain_idx)
          local_type=DOMAIN_MAPPING%globalToLocalMap(global_number)%localType(domain_idx)
          
          IF(domain_no==myGroupComputationNodeNumber) THEN
            ! set local number
            DOMAIN_MAPPING%localToGlobalMap(local_number)=global_number
            
            ! set entry of current element in domain_list
            SELECT CASE(DOMAIN_MAPPING%globalToLocalMap(global_number)%localType(domain_idx))
            CASE(DOMAIN_LOCAL_INTERNAL)
              NUMBER_INTERNAL=NUMBER_INTERNAL+1
              DOMAIN_MAPPING%domainList(NUMBER_INTERNAL)=local_number
            CASE(DOMAIN_LOCAL_BOUNDARY)
              NUMBER_BOUNDARY=NUMBER_BOUNDARY+1
              DOMAIN_MAPPING%domainList(DOMAIN_MAPPING%internalFinish+NUMBER_BOUNDARY)=local_number
            CASE(DOMAIN_LOCAL_GHOST)
              NUMBER_GHOST=NUMBER_GHOST+1
              DOMAIN_MAPPING%domainList(DOMAIN_MAPPING%boundaryFinish+NUMBER_GHOST)=local_number
              
              ! add local number of ghost element to receive list of domain from which to receive
              CALL LIST_ITEM_ADD(GHOST_RECEIVE_LISTS(ADJACENT_DOMAIN_MAP(RECEIVE_FROM_DOMAIN))%PTR,local_number,err,error,*999)              
            CASE DEFAULT
              LOCAL_ERROR="The domain local type of "//TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING%globalToLocalMap( &
                & global_number)%localType(domain_idx),"*",err,error))//" is invalid."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
            
          ELSE IF(SEND_GLOBAL.AND.local_type==DOMAIN_LOCAL_GHOST) THEN
            local_number2=DOMAIN_MAPPING%globalToLocalMap(global_number)%localNumber(1) !The local number for this node
            CALL LIST_ITEM_ADD(GHOST_SEND_LISTS(ADJACENT_DOMAIN_MAP(domain_no))%PTR,local_number2,err,error,*999)            
          ENDIF
        ENDDO !domain_idx
      ENDDO !global_number
      
      ! loop over adjacent domains
      DO domain_idx=1,DOMAIN_MAPPING%numberOfAdjacentDomains
        
        ! transfer the ghost_send_list for adjacent domain to localGhostSendIndices
        CALL LIST_REMOVE_DUPLICATES(GHOST_SEND_LISTS(domain_idx)%PTR,err,error,*999)
        CALL LIST_DETACH_AND_DESTROY(GHOST_SEND_LISTS(domain_idx)%PTR,NUMBER_OF_GHOST_SEND,SEND_LIST,err,error,*999)
        
        ALLOCATE(DOMAIN_MAPPING%adjacentDomains(domain_idx)%localGhostSendIndices(NUMBER_OF_GHOST_SEND),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate local ghost send inidices.",err,error,*999)
        
        DOMAIN_MAPPING%adjacentDomains(domain_idx)%localGhostSendIndices(1:NUMBER_OF_GHOST_SEND)= &
          & SEND_LIST(1:NUMBER_OF_GHOST_SEND)
          
        DOMAIN_MAPPING%adjacentDomains(domain_idx)%numberOfSendGhosts=NUMBER_OF_GHOST_SEND
        DEALLOCATE(SEND_LIST)
        
        ! transfer the ghost_receive_lists for the current adjacent domain to localGhostReceiveIndices
        CALL LIST_REMOVE_DUPLICATES(GHOST_RECEIVE_LISTS(domain_idx)%PTR,err,error,*999)
        CALL LIST_DETACH_AND_DESTROY(GHOST_RECEIVE_LISTS(domain_idx)%PTR,NUMBER_OF_GHOST_RECEIVE,RECEIVE_LIST,err,error,*999)
        
        ALLOCATE(DOMAIN_MAPPING%adjacentDomains(domain_idx)%localGhostReceiveIndices(NUMBER_OF_GHOST_RECEIVE),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate local ghost receive inidices.",err,error,*999)
        
        DOMAIN_MAPPING%adjacentDomains(domain_idx)%localGhostReceiveIndices(1:NUMBER_OF_GHOST_RECEIVE)= &
          & RECEIVE_LIST(1:NUMBER_OF_GHOST_RECEIVE)
        
        DOMAIN_MAPPING%adjacentDomains(domain_idx)%numberOfReceiveGhosts=NUMBER_OF_GHOST_RECEIVE
        DEALLOCATE(RECEIVE_LIST)
        
      ENDDO !domain_idx

      DEALLOCATE(ADJACENT_DOMAIN_MAP)
      DEALLOCATE(GHOST_SEND_LISTS)
      DEALLOCATE(GHOST_RECEIVE_LISTS)
      
      IF(DIAGNOSTICS1) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Domain mappings:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of domains  = ",DOMAIN_MAPPING%numberOfDomains,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of global = ",DOMAIN_MAPPING%numberOfGlobal,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of local = ",DOMAIN_MAPPING%numberOfLocal,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of local = ",DOMAIN_MAPPING%totalNumberOfLocal, &
          & err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Domain numbers:",err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%numberOfDomains,8,8,DOMAIN_MAPPING% &
          & numberOfDomainLocal,'("    Number of domain local :",8(X,I10))','(26X,8(X,I10))',err,error,*999)      
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%numberOfDomains,8,8,DOMAIN_MAPPING% &
          & numberOfDomainGhost,'("    Number of domain ghost :",8(X,I10))','(26X,8(X,I10))',err,error,*999)      
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Domain list:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of internal = ",DOMAIN_MAPPING%numberOfInternal,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary = ",DOMAIN_MAPPING%numberOfBoundary,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost = ",DOMAIN_MAPPING%numberOfGhost,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Internal start = ",DOMAIN_MAPPING%internalStart,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Internal finish = ",DOMAIN_MAPPING%internalFinish,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Boundary start = ",DOMAIN_MAPPING%boundaryStart,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Boundary finish = ",DOMAIN_MAPPING%boundaryFinish,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Ghost start = ",DOMAIN_MAPPING%ghostStart,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Ghost finish = ",DOMAIN_MAPPING%ghostFinish,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,DOMAIN_MAPPING%internalStart,1,DOMAIN_MAPPING%internalFinish,8,8, &
          & DOMAIN_MAPPING%domainList,'("    Internal list :",8(X,I10))','(19X,8(X,I10))',err,error,*999)      
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,DOMAIN_MAPPING%boundaryStart,1,DOMAIN_MAPPING%boundaryFinish,8,8, &
          & DOMAIN_MAPPING%domainList,'("    Boundary list :",8(X,I10))','(19X,8(X,I10))',err,error,*999)      
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,DOMAIN_MAPPING%ghostStart,1,DOMAIN_MAPPING%ghostFinish,8,8, &
          & DOMAIN_MAPPING%domainList,'("    Ghost list    :",8(X,I10))','(19X,8(X,I10))',err,error,*999)      
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map:",err,error,*999)
        DO idx=1,DOMAIN_MAPPING%totalNumberOfLocal
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Local index : ",idx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Global index = ",DOMAIN_MAPPING%localToGlobalMap(idx), &
            & err,error,*999)
        ENDDO !idx
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Global to local map:",err,error,*999)
        DO idx=1,DOMAIN_MAPPING%numberOfGlobal
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global idx : ",idx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of domains  = ", &
            & DOMAIN_MAPPING%globalToLocalMap(idx)%numberOfDomains,err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%globalToLocalMap(idx)% &
            & numberOfDomains,8,8,DOMAIN_MAPPING%globalToLocalMap(idx)%localNumber, &
            & '("      Local number  :",8(X,I10))','(21X,8(X,I10))',err,error,*999)      
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%globalToLocalMap(idx)% &
            & numberOfDomains,8,8,DOMAIN_MAPPING%globalToLocalMap(idx)%domainNumber, &
            & '("      Domain number :",8(X,I10))','(21X,8(X,I10))',err,error,*999)      
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%globalToLocalMap(idx)% &
            & numberOfDomains,8,8,DOMAIN_MAPPING%globalToLocalMap(IDX)%localType, &
            & '("      Local type    :",8(X,I10))','(21X,8(X,I10))',err,error,*999)      
        ENDDO !ne
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
          & DOMAIN_MAPPING%numberOfAdjacentDomains,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%numberOfDomains+1,8,8, &
          & DOMAIN_MAPPING%adjacentDomainsPtr,'("    Adjacent domains ptr  :",8(X,I5))','(27X,8(X,I5))',err,error,*999)
        IF(DOMAIN_MAPPING%numberOfAdjacentDomains>0) THEN
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%adjacentDomainsPtr( &
            & DOMAIN_MAPPING%numberOfDomains)-1,8,8,DOMAIN_MAPPING%adjacentDomainsList, &
            '("    Adjacent domains list :",8(X,I5))','(27X,8(X,I5))',err,error,*999)
          DO domain_idx=1,DOMAIN_MAPPING%numberOfAdjacentDomains
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",domain_idx,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
              & DOMAIN_MAPPING%adjacentDomains(domain_idx)%domainNumber,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
              & DOMAIN_MAPPING%adjacentDomains(domain_idx)%numberOfSendGhosts,err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%adjacentDomains(domain_idx)% &
              & numberOfSendGhosts,8,8,DOMAIN_MAPPING%adjacentDomains(domain_idx)%localGhostSendIndices, &
              & '("      Local send ghost indices        :",8(X,I10))','(39X,8(X,I10))',err,error,*999)      
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of receive ghosts = ", &
              & DOMAIN_MAPPING%adjacentDomains(domain_idx)%numberOfReceiveGhosts,err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%adjacentDomains(domain_idx)% &
              & numberOfReceiveGhosts,8,8,DOMAIN_MAPPING%adjacentDomains(domain_idx)%localGhostReceiveIndices, &
              & '("      Local receive ghost indices     :",8(X,I10))','(39X,8(X,I10))',err,error,*999)              
          ENDDO !domain_idx
        ENDIF
      ENDIF
      
    ELSE
      CALL FlagError("Domain mapping is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE")
    RETURN
999 IF(ALLOCATED(SEND_LIST)) DEALLOCATE(SEND_LIST)
    IF(ALLOCATED(RECEIVE_LIST)) DEALLOCATE(RECEIVE_LIST)
    IF(ALLOCATED(ADJACENT_DOMAIN_MAP)) DEALLOCATE(ADJACENT_DOMAIN_MAP)
    IF(ALLOCATED(adjacentDomains)) DEALLOCATE(adjacentDomains)
    IF(ALLOCATED(GHOST_SEND_LISTS)) THEN
      DO domain_idx=1,SIZE(GHOST_SEND_LISTS)
        IF(ASSOCIATED(GHOST_SEND_LISTS(domain_idx)%PTR)) &
          & CALL LIST_DESTROY(GHOST_SEND_LISTS(domain_idx)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
      ENDDO ! domain_idx
998   DEALLOCATE(GHOST_SEND_LISTS)
    ENDIF
    IF(ALLOCATED(GHOST_RECEIVE_LISTS)) THEN
      DO domain_idx=1,SIZE(GHOST_RECEIVE_LISTS)
        IF(ASSOCIATED(GHOST_RECEIVE_LISTS(domain_idx)%PTR)) &
          & CALL LIST_DESTROY(GHOST_RECEIVE_LISTS(domain_idx)%PTR,DUMMY_ERR,DUMMY_ERROR,*997)
      ENDDO ! domain_idx
997   DEALLOCATE(GHOST_RECEIVE_LISTS)
    ENDIF
    ERRORSEXITS("DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE",err,error)
    RETURN 1
    
  END SUBROUTINE DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Finalises the mapping for a domain mappings mapping and deallocates all memory.
  SUBROUTINE DOMAIN_MAPPINGS_MAPPING_FINALISE(DOMAIN_MAPPING,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: DOMAIN_MAPPING !<A pointer to the domain mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: idx

    ENTERS("DOMAIN_MAPPINGS_MAPPING_FINALISE",err,error,*999)

    IF(ASSOCIATED(DOMAIN_MAPPING)) THEN
      IF(ALLOCATED(DOMAIN_MAPPING%numberOfDomainLocal)) DEALLOCATE(DOMAIN_MAPPING%numberOfDomainLocal)
      IF(ALLOCATED(DOMAIN_MAPPING%numberOfDomainGhost)) DEALLOCATE(DOMAIN_MAPPING%numberOfDomainGhost)
      IF(ALLOCATED(DOMAIN_MAPPING%domainList)) DEALLOCATE(DOMAIN_MAPPING%domainList)
      IF(ALLOCATED(DOMAIN_MAPPING%localToGlobalMap)) DEALLOCATE(DOMAIN_MAPPING%localToGlobalMap)
      IF(ALLOCATED(DOMAIN_MAPPING%globalToLocalMap)) THEN
        DO idx=1,SIZE(DOMAIN_MAPPING%globalToLocalMap,1)
          CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_FINALISE(DOMAIN_MAPPING%globalToLocalMap(idx),err,error,*999)
        ENDDO !idx
        DEALLOCATE(DOMAIN_MAPPING%globalToLocalMap)
      ENDIF
      IF(ALLOCATED(DOMAIN_MAPPING%adjacentDomainsPtr)) DEALLOCATE(DOMAIN_MAPPING%adjacentDomainsPtr)
      IF(ALLOCATED(DOMAIN_MAPPING%adjacentDomainsList)) DEALLOCATE(DOMAIN_MAPPING%adjacentDomainsList)
      IF(ALLOCATED(DOMAIN_MAPPING%adjacentDomains)) THEN
        DO idx=1,SIZE(DOMAIN_MAPPING%adjacentDomains,1)
          CALL DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE(DOMAIN_MAPPING%adjacentDomains(idx),err,error,*999)
        ENDDO !idx
        DEALLOCATE(DOMAIN_MAPPING%adjacentDomains)
      ENDIF
      DEALLOCATE(DOMAIN_MAPPING)
    ENDIF
    
    EXITS("DOMAIN_MAPPINGS_MAPPING_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_MAPPING_FINALISE",err,error)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_MAPPING_FINALISE
  
  !
  !================================================================================================================================
  !

  !> Finalises the global mapping in the given domain mappings.
  SUBROUTINE DOMAIN_MAPPINGS_MAPPING_GLOBAL_FINALISE(MAPPING_GLOBAL_MAP,err,error,*)

    !Argument variables
    TYPE(DomainGlobalMappingType) :: MAPPING_GLOBAL_MAP !<The domain global mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<Th error string
    !Local Variables

    ENTERS("DOMAIN_MAPPINGS_MAPPING_GLOBAL_FINALISE",err,error,*999)

    IF(ALLOCATED(MAPPING_GLOBAL_MAP%localNumber)) DEALLOCATE(MAPPING_GLOBAL_MAP%localNumber)
    IF(ALLOCATED(MAPPING_GLOBAL_MAP%domainNumber)) DEALLOCATE(MAPPING_GLOBAL_MAP%domainNumber)
    IF(ALLOCATED(MAPPING_GLOBAL_MAP%localType)) DEALLOCATE(MAPPING_GLOBAL_MAP%localType)
 
    EXITS("DOMAIN_MAPPINGS_MAPPING_GLOBAL_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_MAPPING_GLOBAL_FINALISE",err,error)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_MAPPING_GLOBAL_FINALISE

  !
  !================================================================================================================================
  !

  !>Finalises the global mapping in the given domain mappings.
  SUBROUTINE DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(MAPPING_GLOBAL_MAP,err,error,*)

    !Argument variables
    TYPE(DomainGlobalMappingType) :: MAPPING_GLOBAL_MAP !<The domain global mapping to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE",err,error,*999)

    MAPPING_GLOBAL_MAP%numberOfDomains=0
    
    EXITS("DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE",err,error)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE

  !
  !================================================================================================================================
  !

  !>Initialises the mapping for a domain mappings mapping.
  SUBROUTINE DomainMappings_MappingInitialise(workGroup,domainMapping,err,error,*)

    !Argument variables
    TYPE(WorkGroupType), POINTER :: workGroup !<A pointer to the work group to initialise the mappings for
    TYPE(DomainMappingType), POINTER :: domainMapping !<A pointer to the domain mapping to initialise the mappings for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfDomains
 
    ENTERS("DomainMappings_MappingInitialise",err,error,*999)

    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("WorkGroup is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*999)

    CALL WorkGroup_NumberOfGroupNodesGet(workGroup,numberOfDomains,err,error,*999)

    domainMapping%workGroup=>workGroup
    domainMapping%totalNumberOfLocal=0
    domainMapping%numberOfLocal=0
    domainMapping%numberOfGlobal=0
    domainMapping%numberOfDomains=numberOfDomains
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
    
    EXITS("DomainMappings_MappingInitialise")
    RETURN
999 ERRORSEXITS("DomainMappings_MappingInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMappings_MappingInitialise

  !
  !================================================================================================================================
  !
  
END MODULE DOMAIN_MAPPINGS
