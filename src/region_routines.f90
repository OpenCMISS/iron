!> \file
!> \author Chris Bradley
!> \brief This module contains all region routines.
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

!> This module contains all region routines.
MODULE REGION_ROUTINES

  USE BaseRoutines
  USE ContextAccessRoutines
  USE COORDINATE_ROUTINES
  USE CoordinateSystemAccessRoutines
  USE CMISS_CELLML
  USE DataPointRoutines
  USE EQUATIONS_SET_ROUTINES
  USE FIELD_ROUTINES
  USE GENERATED_MESH_ROUTINES
  USE INPUT_OUTPUT
  USE INTERFACE_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE MESH_ROUTINES
  USE NodeRoutines
  USE RegionAccessRoutines
  USE STRINGS
  USE TYPES

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  INTERFACE REGION_LABEL_GET
    MODULE PROCEDURE REGION_LABEL_GET_C
    MODULE PROCEDURE REGION_LABEL_GET_VS
  END INTERFACE REGION_LABEL_GET
  
  INTERFACE REGION_LABEL_SET
    MODULE PROCEDURE REGION_LABEL_SET_C
    MODULE PROCEDURE REGION_LABEL_SET_VS
  END INTERFACE REGION_LABEL_SET

  PUBLIC REGION_COORDINATE_SYSTEM_SET

  PUBLIC REGION_CREATE_START,REGION_CREATE_FINISH

  PUBLIC REGION_DESTROY

  PUBLIC REGION_INITIALISE,REGION_FINALISE

  PUBLIC REGION_LABEL_GET,REGION_LABEL_SET
  
  PUBLIC Regions_Initialise,Regions_Finalise

CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets the coordinate system of region.  \see OPENCMISS::Iron::cmfe_RegionCoordinateSystemSet
  SUBROUTINE REGION_COORDINATE_SYSTEM_SET(REGION,COORDINATE_SYSTEM,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: REGION !<A pointer to the region to set the coordinate system for
    TYPE(CoordinateSystemType), POINTER :: COORDINATE_SYSTEM !<The coordinate system to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("REGION_COORDINATE_SYSTEM_SET",err,error,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(REGION%regionFinished) THEN
        CALL FlagError("Region has been finished.",err,error,*999)
      ELSE
        IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
          IF(COORDINATE_SYSTEM%coordinateSystemFinished) THEN
            REGION%coordinateSystem=>COORDINATE_SYSTEM
          ELSE
            CALL FlagError("Coordinate system has not been finished.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Coordinate system is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",err,error,*999)
    ENDIF
    
    EXITS("REGION_COORDINATE_SYSTEM_SET")
    RETURN
999 ERRORSEXITS("REGION_COORDINATE_SYSTEM_SET",err,error)
    RETURN 1
  END SUBROUTINE REGION_COORDINATE_SYSTEM_SET
  
  !
  !================================================================================================================================
  !

  !>Finishes the creation of a region. \see OPENCMISS::Iron::cmfe_RegionCreateFinish
  SUBROUTINE REGION_CREATE_FINISH(REGION,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: REGION !<A pointer to the region to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
     
    ENTERS("REGION_CREATE_FINISH",err,error,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(REGION%regionFinished) THEN
        CALL FlagError("Region has already been finished.",err,error,*999)
      ELSE
        REGION%regionFinished=.TRUE.
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",err,error,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Region : ",REGION%userNumber,err,error,*999)      
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Label = ",REGION%label,err,error,*999)
      IF(ASSOCIATED(REGION%parentRegion)) THEN
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Parent region user number = ",REGION%parentRegion%userNumber, &
          & err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Parent region label = ",REGION%parentRegion%label, &
          & err,error,*999)        
      ENDIF
    ENDIF
    
    EXITS("REGION_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("REGION_CREATE_FINISH",err,error)
    RETURN 1
  END SUBROUTINE REGION_CREATE_FINISH

  !
  !================================================================================================================================
  !
  
  !>Starts the creation a new region number USER_NUMBER as a sub region to the given PARENT_REGION, initialises all
  !>variables and inherits the PARENT_REGIONS coordinate system. \see OPENCMISS::Iron::cmfe_RegionCreateFinish
  !>Default values set for the REGION's attributes are:
  !>- COORDINATE_SYSTEM: parent coordinate system. See \ref COORDINATE_SYSTEM_TYPE
  !>- DATA_POINTS: null
  !>- NODES: null
  !>- MESHES: 0 mesh
  !>- FIELDS: 0 field
  !>- EQUATIONS_SETS: 0 equation set
  !>- PARENT_REGION: global region
  !>- numberOfSubRegions: 0
  !>- subRegions: 0 region
  SUBROUTINE REGION_CREATE_START(USER_NUMBER,PARENT_REGION,REGION,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the region to create
    TYPE(RegionType), POINTER :: PARENT_REGION !<A pointer to the parent region
    TYPE(RegionType), POINTER :: REGION !<On exit, a pointer to the created region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,region_idx
    TYPE(RegionType), POINTER :: NEW_REGION
    TYPE(RegionPtrType), POINTER :: NEW_subRegions(:)
    TYPE(RegionsType), POINTER :: regions
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR,LOCAL_STRING

    NULLIFY(NEW_REGION)
    NULLIFY(NEW_subRegions)
    
    ENTERS("REGION_CREATE_START",err,error,*997)

    NULLIFY(regions)
    CALL Region_RegionsGet(PARENT_REGION,regions,err,error,*999)
    CALL REGION_USER_NUMBER_FIND(regions,USER_NUMBER,NEW_REGION,err,error,*997)
    IF(ASSOCIATED(NEW_REGION)) THEN
      LOCAL_ERROR="Region number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",err,error))// &
        & " has already been created."
      CALL FlagError(LOCAL_ERROR,err,error,*997)
    ELSE
      IF(ASSOCIATED(REGION)) THEN
        CALL FlagError("Region is already associated.",err,error,*997)
      ELSE
        NULLIFY(REGION)
        IF(ASSOCIATED(PARENT_REGION)) THEN
          IF(PARENT_REGION%regionFinished) THEN
            IF(ASSOCIATED(PARENT_REGION%coordinateSystem)) THEN
              !Initialise the region
              CALL REGION_INITIALISE(REGION,err,error,*999)
              !Set the user number
              REGION%userNumber=USER_NUMBER
              region%regions=>PARENT_REGION%regions
              !CPB 21/02/07 The vstring operation crashes the AIX compiler so put a CHAR() etc. around it.
              !REGION%label="Region "//NUMBER_TO_VSTRING(USER_NUMBER,"*",err,error)
              LOCAL_STRING="Region "//NUMBER_TO_VSTRING(USER_NUMBER,"*",err,error)
              REGION%label=CHAR(LOCAL_STRING)
              IF(err/=0) GOTO 999
              REGION%coordinateSystem=>PARENT_REGION%coordinateSystem
              !Adjust the parent region to include this new daughter
              ALLOCATE(NEW_subRegions(PARENT_REGION%numberOfSubRegions+1),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate new sub-regions.",err,error,*999)
              DO region_idx=1,PARENT_REGION%numberOfSubRegions
                NEW_subRegions(region_idx)%PTR=>PARENT_REGION%subRegions(region_idx)%PTR
              ENDDO !region_no
              PARENT_REGION%numberOfSubRegions=PARENT_REGION%numberOfSubRegions+1
              NEW_subRegions(PARENT_REGION%numberOfSubRegions)%PTR=>REGION
              IF(ASSOCIATED(PARENT_REGION%subRegions)) DEALLOCATE(PARENT_REGION%subRegions)
              PARENT_REGION%subRegions=>NEW_subRegions
              !Set the new regions parent region to the parent region
              REGION%parentRegion=>PARENT_REGION
            ELSE
              CALL FlagError("Parent region does not have an associated coordinate system.",err,error,*997)
            ENDIF
          ELSE
            CALL FlagError("Parent region has not been finished.",err,error,*997)
          ENDIF
        ELSE
          CALL FlagError("Parent region is not associated.",err,error,*997)
        ENDIF
      ENDIF
    ENDIF
    
    EXITS("REGION_CREATE_START")
    RETURN
999 CALL REGION_FINALISE(REGION,DUMMY_ERR,DUMMY_ERROR,*998)
998 IF(ASSOCIATED(NEW_subRegions)) DEALLOCATE(NEW_subRegions)
997 ERRORSEXITS("REGION_CREATE_START",err,error)
    RETURN 1
  END SUBROUTINE REGION_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroys a region given by USER_NUMBER and all sub-regions under it. \todo create destroy by pointer method. \see OpenCMISS::Iron::cmfe_Region_Destroy
  RECURSIVE SUBROUTINE REGION_DESTROY_NUMBER(regions,USER_NUMBER,err,error,*)

    !Argument variables
    TYPE(RegionsType), POINTER :: regions !<A pointer to the regions
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the region to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: count,nr
    TYPE(RegionType), POINTER :: REGION
    TYPE(RegionPtrType), POINTER :: NEW_subRegions(:)

    ENTERS("REGION_DESTROY_NUMBER",err,error,*999)

    NULLIFY(REGION)
    CALL REGION_USER_NUMBER_FIND(regions,USER_NUMBER,REGION,err,error,*999)
    IF(ASSOCIATED(REGION)) THEN

!!NOTE: We have to find a pointer to the region to destroy within this routine rather than passing in a pointer to a
!!DESTROY_REGION_PTR type routine because we need to change REGION%subRegions of the PARENT region and this would violate section
!!12.4.1.6 of the Fortran standard if the dummy REGION pointer argument was associated with the subRegions(x)%PTR actual
!!argument.
      
      IF(REGION%numberOfSubRegions==0) THEN
        !No more daughter sub regions so delete this instance
        IF(ASSOCIATED(REGION%parentRegion)) THEN
          NULLIFY(NEW_subRegions)
          IF(REGION%parentRegion%numberOfSubRegions>1) THEN
            !If the parent region has more than one sub regions then remove this instance from its sub-regions list 
            ALLOCATE(NEW_subRegions(REGION%parentRegion%numberOfSubRegions-1),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate new sub-regions.",err,error,*999)
            count=0
            DO nr=1,REGION%parentRegion%numberOfSubRegions
              IF(REGION%parentRegion%subRegions(nr)%PTR%userNumber/=REGION%userNumber) THEN
                count=count+1
                NEW_subRegions(count)%PTR=>REGION%parentRegion%subRegions(nr)%PTR
              ENDIF
            ENDDO !nr
          ENDIF
          REGION%parentRegion%numberOfSubRegions=REGION%parentRegion%numberOfSubRegions-1
          IF(ASSOCIATED(REGION%parentRegion%subRegions)) DEALLOCATE(REGION%parentRegion%subRegions)
          REGION%parentRegion%subRegions=>NEW_subRegions
          !Finalise the region
          CALL REGION_FINALISE(REGION,err,error,*999)
        ELSE
          CALL FlagError("Parent region is not associated.",err,error,*999)
        ENDIF
      ELSE
        !Recursively delete sub regions first
        DO WHILE(REGION%numberOfSubRegions>0)
          CALL REGION_DESTROY_NUMBER(regions,REGION%subRegions(1)%PTR%userNumber,err,error,*999)
        ENDDO
        !Now delete this instance
        CALL REGION_DESTROY_NUMBER(regions,REGION%userNumber,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Region number does not exist.",err,error,*999)
    ENDIF
    
    EXITS("REGION_DESTROY_NUMBER")
    RETURN
999 ERRORSEXITS("REGION_DESTROY_NUMBER",err,error)
    RETURN 1
  END SUBROUTINE REGION_DESTROY_NUMBER

  !
  !================================================================================================================================
  !

  !>Destroys a region identified by a pointer and all sub-regions under it. \see OPENCMISS::Iron::cmfe_RegionDestroy
  SUBROUTINE REGION_DESTROY(REGION,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: REGION !<A pointer to the region to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: USER_NUMBER
    TYPE(RegionsType), POINTER :: regions
 
    ENTERS("REGION_DESTROY",err,error,*999)
    
    IF(ASSOCIATED(REGION)) THEN
      NULLIFY(regions)
      CALL Region_RegionsGet(region,regions,err,error,*999)
      USER_NUMBER=REGION%userNumber
      CALL REGION_DESTROY_NUMBER(regions,USER_NUMBER,err,error,*999)
    ELSE
      CALL FlagError("Region is not associated.",err,error,*999)
    ENDIF
    
    EXITS("REGION_DESTROY")
    RETURN
999 ERRORSEXITS("REGION_DESTROY",err,error)
    RETURN 1
  END SUBROUTINE REGION_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises a region and deallocates all memory
  SUBROUTINE REGION_FINALISE(REGION,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: REGION !<A pointer to the region to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("REGION_FINALISE",err,error,*999)
    
    IF(ASSOCIATED(REGION)) THEN
      REGION%label=""
      CALL CELLML_ENVIRONMENTS_FINALISE(REGION%CELLML_ENVIRONMENTS,err,error,*999)
      CALL EQUATIONS_SETS_FINALISE(REGION,err,error,*999)
      CALL FIELDS_FINALISE(REGION%FIELDS,err,error,*999)
      CALL MESHES_FINALISE(REGION%MESHES,err,error,*999)
      CALL DataPointSets_Finalise(region%dataPointSets,err,error,*999)
      IF(ASSOCIATED(REGION%NODES)) CALL NODES_DESTROY(REGION%NODES,err,error,*999)
      IF(ASSOCIATED(REGION%subRegions)) DEALLOCATE(REGION%subRegions)
      IF(ASSOCIATED(REGION%INTERFACES)) CALL INTERFACES_FINALISE(REGION%INTERFACES,err,error,*999)
      IF(ASSOCIATED(REGION%generatedMeshes)) CALL GENERATED_MESHES_FINALISE(REGION%generatedMeshes,err,error,*999)
      DEALLOCATE(REGION)
    ENDIF
    
    EXITS("REGION_FINALISE")
    RETURN
999 ERRORSEXITS("REGION_FINALISE",err,error)
    RETURN 1
  END SUBROUTINE REGION_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a region.
  SUBROUTINE REGION_INITIALISE(REGION,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: REGION !<A pointer to the region to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
   
    ENTERS("REGION_INITIALISE",err,error,*998)

    IF(ASSOCIATED(REGION)) THEN
      CALL FlagError("Region is already associated.",err,error,*998)
    ELSE
      ALLOCATE(REGION,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate region.",err,error,*999)
      REGION%userNumber=0
      NULLIFY(region%regions)
      REGION%regionFinished=.FALSE.
      REGION%label=""
      NULLIFY(REGION%coordinateSystem)
      NULLIFY(region%dataPointSets)
      NULLIFY(REGION%NODES)
      NULLIFY(REGION%MESHES)
      NULLIFY(REGION%generatedMeshes)
      NULLIFY(REGION%FIELDS)
      NULLIFY(REGION%EQUATIONS_SETS)
      NULLIFY(REGION%CELLML_ENVIRONMENTS)
      NULLIFY(REGION%parentRegion)
      REGION%numberOfSubRegions=0
      NULLIFY(REGION%subRegions)
      NULLIFY(REGION%INTERFACES)
      CALL DataPointSets_Initialise(region,err,error,*999)
      CALL MESHES_INITIALISE(REGION,err,error,*999)
      CALL GENERATED_MESHES_INITIALISE(REGION,err,error,*999)
      CALL FIELDS_INITIALISE(REGION,err,error,*999)
      CALL EQUATIONS_SETS_INITIALISE(REGION,err,error,*999)
      CALL CELLML_ENVIRONMENTS_INITIALISE(REGION,err,error,*999)
      CALL INTERFACES_INITIALISE(REGION,err,error,*999)
    ENDIF
    
    EXITS("REGION_INITIALISE")
    RETURN
999 CALL REGION_FINALISE(REGION,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("REGION_INITIALISE",err,error)
    RETURN 1
    
  END SUBROUTINE REGION_INITIALISE

  !
  !================================================================================================================================
  !

  !>Returns the label of a region. \see OPENCMISS::Iron::cmfe_RegionLabelGet
  SUBROUTINE REGION_LABEL_GET_C(REGION,LABEL,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: REGION !<A pointer to the region to get the label for
    CHARACTER(LEN=*), INTENT(OUT) :: LABEL !<On return the region label.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: C_LENGTH,VS_LENGTH

    ENTERS("REGION_LABEL_GET_C",err,error,*999)

    IF(ASSOCIATED(REGION)) THEN
      C_LENGTH=LEN(LABEL)
      VS_LENGTH=LEN_TRIM(REGION%label)
      IF(C_LENGTH>VS_LENGTH) THEN
        LABEL=CHAR(REGION%label,VS_LENGTH)
      ELSE
        LABEL=CHAR(REGION%label,C_LENGTH)
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",err,error,*999)
    ENDIF
    
    EXITS("REGION_LABEL_GET_C")
    RETURN
999 ERRORSEXITS("REGION_LABEL_GET_C",err,error)
    RETURN 1
    
  END SUBROUTINE REGION_LABEL_GET_C

   !
  !================================================================================================================================
  !

  !>Returns the label of a region. \see OPENCMISS::Iron::cmfe_RegionLabelGet
  SUBROUTINE REGION_LABEL_GET_VS(REGION,LABEL,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: REGION !<A pointer to the region to get the label for
    TYPE(VARYING_STRING), INTENT(OUT) :: LABEL !<On return the region label.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("REGION_LABEL_GET_VS",err,error,*999)

    IF(ASSOCIATED(REGION)) THEN
      !CPB 20/2/07 The following line crashes the AIX compiler unless it has a VAR_STR(CHAR()) around it
      LABEL=VAR_STR(CHAR(REGION%label))
    ELSE
      CALL FlagError("Region is not associated.",err,error,*999)
    ENDIF
    
    EXITS("REGION_LABEL_GET_VS")
    RETURN
999 ERRORSEXITS("REGION_LABEL_GET_VS",err,error)
    RETURN 1
    
  END SUBROUTINE REGION_LABEL_GET_VS

  !
  !================================================================================================================================
  !

  !>Sets the label of a region. \see OPENCMISS::Iron::cmfe_RegionLabelSet
  SUBROUTINE REGION_LABEL_SET_C(REGION,LABEL,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: REGION !<A pointer to the region to set the label for 
    CHARACTER(LEN=*), INTENT(IN) :: LABEL !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("REGION_LABEL_SET_C",err,error,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(REGION%regionFinished) THEN
        CALL FlagError("Region has been finished.",err,error,*999)
      ELSE
        REGION%label=LABEL
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",err,error,*999)
    ENDIF
    
    EXITS("REGION_LABEL_SET_C")
    RETURN
999 ERRORSEXITS("REGION_LABEL_SET_C",err,error)
    RETURN 1
  END SUBROUTINE REGION_LABEL_SET_C

  !
  !================================================================================================================================
  !

  !>Sets the label of a region. \see OPENCMISS::Iron::cmfe_RegionLabelSet
  SUBROUTINE REGION_LABEL_SET_VS(REGION,LABEL,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: REGION !<A pointer to the region to set the label for 
    TYPE(VARYING_STRING), INTENT(IN) :: LABEL !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("REGION_LABEL_SET_VS",err,error,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(REGION%regionFinished) THEN
        CALL FlagError("Region has been finished.",err,error,*999)
      ELSE
        REGION%label=LABEL
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",err,error,*999)
    ENDIF
    
    EXITS("REGION_LABEL_SET_VS")
    RETURN
999 ERRORSEXITS("REGION_LABEL_SET_VS",err,error)
    RETURN 1
  END SUBROUTINE REGION_LABEL_SET_VS

  !
  !================================================================================================================================
  !

  !>Finalises the regions and destroys any current regions.
  SUBROUTINE Regions_Finalise(regions,err,error,*)

    !Argument variables
    TYPE(RegionsType), POINTER :: regions !<A pointer to the regions to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Regions_Finalise",err,error,*999)

    IF(ASSOCIATED(regions)) THEN
      IF(ASSOCIATED(regions%worldRegion)) THEN
        !Destroy any global region daughter regions first
        DO WHILE(regions%worldRegion%numberOfSubRegions>0)
          CALL REGION_DESTROY_NUMBER(regions,regions%worldRegion%subRegions(1)%PTR%userNumber,err,error,*999)
        ENDDO !region
        !Destroy global region and deallocated any memory allocated in the global region
        CALL REGION_FINALISE(regions%worldRegion,err,error,*999)
      ENDIF
      DEALLOCATE(regions)
    ENDIF
   
    EXITS("Regions_Finalise")
    RETURN
999 ERRORSEXITS("Regions_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Regions_Finalise
  
  !
  !================================================================================================================================
  !

  !>Initialises the regions and creates the global world region for a context.
  SUBROUTINE Regions_Initialise(context,err,error,*)    

    !Argument variables
    TYPE(ContextType), POINTER :: context !<A pointer to the context to initialise the regions for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(CoordinateSystemType), POINTER :: worldCoordinateSystem
    TYPE(CoordinateSystemsType), POINTER :: coordinateSystems
    TYPE(VARYING_STRING) :: dummyError

    NULLIFY(worldCoordinateSystem)
    
    ENTERS("Regions_Initialise",err,error,*998)

    IF(.NOT.ASSOCIATED(context)) CALL FlagError("Context is not associated.",err,error,*998)
    IF(ASSOCIATED(context%regions)) CALL FlagError("Context regions is already associated.",err,error,*998)

    NULLIFY(coordinateSystems)
    CALL Context_CoordinateSystemsGet(context,coordinateSystems,err,error,*998)
    CALL CoordinateSystem_UserNumberFind(coordinateSystems,0,worldCoordinateSystem,err,error,*998)
    IF(.NOT.ASSOCIATED(worldCoordinateSystem)) CALL FlagError("World coordinate system has not been created.",err,error,*998)

    !Allocate context regions
    ALLOCATE(context%regions,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate regions.",err,error,*999)
    !Initialise
    context%regions%context=>context
    NULLIFY(context%regions%worldRegion)
    !Create world region 
    CALL Region_Initialise(context%regions%worldRegion,err,error,*999)
    context%regions%worldRegion%userNumber=0
    context%regions%worldRegion%regions=>context%regions
    context%regions%worldRegion%label="World Region"
    context%regions%worldRegion%coordinateSystem=>worldCoordinateSystem
    context%regions%worldRegion%regionFinished=.TRUE.
   
    EXITS("Regions_Initialise")
    RETURN
999 CALL Regions_Finalise(context%regions,dummyErr,dummyError,*998)
998 ERRORSEXITS("Regions_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Regions_Initialise

  !
  !================================================================================================================================
  !
  
END MODULE REGION_ROUTINES
