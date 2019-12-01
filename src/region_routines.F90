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
MODULE RegionRoutines

  USE BaseRoutines
  USE ContextAccessRoutines
  USE CoordinateSystemRoutines
  USE CoordinateSystemAccessRoutines
  USE CMISS_CELLML
  USE DataPointRoutines
  USE DecompositionRoutines
  USE EquationsSetRoutines
  USE FieldRoutines
  USE GENERATED_MESH_ROUTINES
  USE InputOutput
  USE InterfaceRoutines
  USE ISO_VARYING_STRING
  USE KINDS
  USE MeshRoutines
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

  INTERFACE Region_LabelSet
    MODULE PROCEDURE Region_LabelSetC
    MODULE PROCEDURE Region_LabelSetVS
  END INTERFACE Region_LabelSet

  PUBLIC Region_CoordinateSystemSet

  PUBLIC Region_CreateStart,Region_CreateFinish

  PUBLIC Region_Destroy

  PUBLIC Region_Initialise,Region_Finalise

  PUBLIC Region_LabelSet
  
  PUBLIC Regions_Initialise,Regions_Finalise

CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets the coordinate system of region.  \see OPENCMISS::Iron::cmfe_Region_CoordinateSystemSet
  SUBROUTINE Region_CoordinateSystemSet(region,coordinateSystem,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to set the coordinate system for
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<The coordinate system to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Region_CoordinateSystemSet",err,error,*999)

    CALL Region_AssertNotFinished(region,err,error,*999)
    CALL CoordinateSystem_AssertIsFinished(coordinateSystem,err,error,*999)
    
    region%coordinateSystem=>coordinateSystem
     
    EXITS("Region_CoordinateSystemSet")
    RETURN
999 ERRORSEXITS("Region_CoordinateSystemSet",err,error)
    RETURN 1
    
  END SUBROUTINE Region_CoordinateSystemSet
  
  !
  !================================================================================================================================
  !

  !>Finishes the creation of a region. \see OPENCMISS::Iron::cmfe_Region_CreateFinish
  SUBROUTINE Region_CreateFinish(region,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
     
    ENTERS("Region_CreateFinish",err,error,*999)

    CALL Region_AssertNotFinished(region,err,error,*999)
    
    region%regionFinished=.TRUE.
    
    IF(diagnostics1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Region : ",region%userNumber,err,error,*999)      
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Label = ",region%label,err,error,*999)
      IF(ASSOCIATED(region%parentRegion)) THEN
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Parent region user number = ",region%parentRegion%userNumber, &
          & err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Parent region label = ",region%parentRegion%label, &
          & err,error,*999)        
      ENDIF
    ENDIF
    
    EXITS("Region_CreateFinish")
    RETURN
999 ERRORSEXITS("Region_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE Region_CreateFinish

  !
  !================================================================================================================================
  !
  
  !>Starts the creation a new region number userNumber as a sub region to the given parentRegion, initialises all
  !>variables and inherits the parentRegionS coordinate system. \see OPENCMISS::Iron::cmfe_Region_CreateFinish
  !>Default values set for the region's attributes are:
  !>- coordinateSystem: parent coordinate system. See \ref COORDINATE_SYSTEM_TYPE
  !>- DATA_POINTS: null
  !>- NODES: null
  !>- MESHES: 0 mesh
  !>- FIELDS: 0 field
  !>- equationsSets: 0 equation set
  !>- parentRegion: global region
  !>- numberOfSubRegions: 0
  !>- subRegions: 0 region
  SUBROUTINE Region_CreateStart(userNumber,parentRegion,region,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the region to create
    TYPE(RegionType), POINTER :: parentRegion !<A pointer to the parent region
    TYPE(RegionType), POINTER :: region !<On exit, a pointer to the created region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,regionIdx
    TYPE(RegionType), POINTER :: newRegion
    TYPE(RegionPtrType), POINTER :: newSubRegions(:)
    TYPE(RegionsType), POINTER :: regions
    TYPE(VARYING_STRING) :: dummyError,localError,localString

    NULLIFY(newRegion)
    NULLIFY(newSubRegions)
    
    ENTERS("Region_CreateStart",err,error,*997)

    IF(ASSOCIATED(region)) CALL FlagError("Region is already associated.",err,error,*997)
    IF(.NOT.ASSOCIATED(parentRegion)) CALL FlagError("Parent region is not associated.",err,error,*997)
    IF(.NOT.parentRegion%regionFinished)  CALL FlagError("Parent region has not been finished.",err,error,*997)
    IF(.NOT.ASSOCIATED(parentRegion%coordinateSystem)) &
      & CALL FlagError("Parent region does not have an associated coordinate system.",err,error,*997)      
    NULLIFY(regions)
    CALL Region_RegionsGet(parentRegion,regions,err,error,*999)
    CALL Region_UserNumberFind(regions,userNumber,newRegion,err,error,*997)
    IF(ASSOCIATED(newRegion)) THEN
      localError="Region number "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " has already been created."
      CALL FlagError(localError,err,error,*997)
    ENDIF
    
    NULLIFY(region)
    !Initialise the region
    CALL Region_Initialise(region,err,error,*999)
    !Set the user number
    region%userNumber=userNumber
    region%regions=>parentRegion%regions
    !CPB 21/02/07 The vstring operation crashes the AIX compiler so put a CHAR() etc. around it.
    !region%label="Region "//NumberToVString(userNumber,"*",err,error)
    localString="Region "//NumberToVString(userNumber,"*",err,error)
    region%label=CHAR(localString)
    IF(err/=0) GOTO 999
    region%coordinateSystem=>parentRegion%coordinateSystem
    !Adjust the parent region to include this new daughter
    ALLOCATE(newSubRegions(parentRegion%numberOfSubRegions+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new sub-regions.",err,error,*999)
    DO regionIdx=1,parentRegion%numberOfSubRegions
      newSubRegions(regionIdx)%PTR=>parentRegion%subRegions(regionIdx)%PTR
    ENDDO !region_no
    parentRegion%numberOfSubRegions=parentRegion%numberOfSubRegions+1
    newSubRegions(parentRegion%numberOfSubRegions)%PTR=>region
    IF(ASSOCIATED(parentRegion%subRegions)) DEALLOCATE(parentRegion%subRegions)
    parentRegion%subRegions=>newSubRegions
    !Set the new regions parent region to the parent region
    region%parentRegion=>parentRegion
   
    EXITS("Region_CreateStart")
    RETURN
999 CALL Region_Finalise(region,dummyErr,dummyError,*998)
998 IF(ASSOCIATED(newSubRegions)) DEALLOCATE(newSubRegions)
997 ERRORSEXITS("Region_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE Region_CreateStart

  !
  !================================================================================================================================
  !

  !>Destroys a region given by userNumber and all sub-regions under it. \todo create destroy by pointer method. \see OpenCMISS::Iron::cmfe_Region_Destroy
  RECURSIVE SUBROUTINE Region_DestroyNumber(regions,userNumber,err,error,*)

    !Argument variables
    TYPE(RegionsType), POINTER :: regions !<A pointer to the regions
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the region to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: count,regionIdx
    TYPE(RegionType), POINTER :: region
    TYPE(RegionPtrType), POINTER :: newSubRegions(:)
    TYPE(VARYING_STRING) :: localError

    ENTERS("Region_DestroyNumber",err,error,*999)

    IF(.NOT.ASSOCIATED(regions)) CALL FlagError("Regions is not associated.",err,error,*999)
    NULLIFY(region)
    CALL Region_UserNumberFind(regions,userNumber,region,err,error,*999)
    IF(.NOT.ASSOCIATED(region)) THEN
      localError="Region number "//TRIM(NumberToVString(userNumber,"*",err,error))//" does not exist."
      CALL FlagError(localError,err,error,*999)
    ENDIF

!!NOTE: We have to find a pointer to the region to destroy within this routine rather than passing in a pointer to a
!!DESTROY_REGION_PTR type routine because we need to change region%subRegions of the PARENT region and this would violate section
!!12.4.1.6 of the Fortran standard if the dummy region pointer argument was associated with the subRegions(x)%PTR actual
!!argument.
      
    IF(region%numberOfSubRegions==0) THEN
      !No more daughter sub regions so delete this instance
      IF(.NOT.ASSOCIATED(region%parentRegion)) CALL FlagError("Parent region is not associated.",err,error,*999)
      NULLIFY(newSubRegions)
      IF(region%parentRegion%numberOfSubRegions>1) THEN
        !If the parent region has more than one sub regions then remove this instance from its sub-regions list 
        ALLOCATE(newSubRegions(region%parentRegion%numberOfSubRegions-1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new sub-regions.",err,error,*999)
        count=0
        DO regionIdx=1,region%parentRegion%numberOfSubRegions
          IF(region%parentRegion%subRegions(regionIdx)%PTR%userNumber/=region%userNumber) THEN
            count=count+1
            newSubRegions(count)%PTR=>region%parentRegion%subRegions(regionIdx)%PTR
          ENDIF
        ENDDO !regionIdx
      ENDIF
      region%parentRegion%numberOfSubRegions=region%parentRegion%numberOfSubRegions-1
      IF(ASSOCIATED(region%parentRegion%subRegions)) DEALLOCATE(region%parentRegion%subRegions)
      region%parentRegion%subRegions=>newSubRegions
      !Finalise the region
      CALL Region_Finalise(region,err,error,*999)
   ELSE
      !Recursively delete sub regions first
      DO WHILE(region%numberOfSubRegions>0)
        CALL Region_DestroyNumber(regions,region%subRegions(1)%PTR%userNumber,err,error,*999)
      ENDDO
      !Now delete this instance
      CALL Region_DestroyNumber(regions,region%userNumber,err,error,*999)
    ENDIF
    
    EXITS("Region_DestroyNumber")
    RETURN
999 ERRORSEXITS("Region_DestroyNumber",err,error)
    RETURN 1
    
  END SUBROUTINE Region_DestroyNumber

  !
  !================================================================================================================================
  !

  !>Destroys a region identified by a pointer and all sub-regions under it. \see OPENCMISS::Iron::cmfe_Region_Destroy
  SUBROUTINE Region_Destroy(region,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: userNumber
    TYPE(RegionsType), POINTER :: regions
 
    ENTERS("Region_Destroy",err,error,*999)
    
    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)
    
    NULLIFY(regions)
    CALL Region_RegionsGet(region,regions,err,error,*999)
    userNumber=region%userNumber
    CALL Region_DestroyNumber(regions,userNumber,err,error,*999)
    
    EXITS("Region_Destroy")
    RETURN
999 ERRORSEXITS("Region_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE Region_Destroy

  !
  !================================================================================================================================
  !

  !>Finalises a region and deallocates all memory
  SUBROUTINE Region_Finalise(region,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("Region_Finalise",err,error,*999)
    
    IF(ASSOCIATED(region)) THEN
      region%label=""
      CALL CELLML_ENVIRONMENTS_FINALISE(region%cellMLEnvironments,err,error,*999)
      CALL EquationsSets_Finalise(region,err,error,*999)
      CALL Fields_Finalise(region%fields,err,error,*999)
      CALL Decomposers_Finalise(region%decomposers,err,error,*999)
      IF(ASSOCIATED(region%generatedMeshes)) CALL GENERATED_MESHES_FINALISE(region%generatedMeshes,err,error,*999)
      CALL Meshes_Finalise(region%meshes,err,error,*999)
      CALL DataPointSets_Finalise(region%dataPointSets,err,error,*999)
      IF(ASSOCIATED(region%nodes)) CALL Nodes_Destroy(region%nodes,err,error,*999)
      IF(ASSOCIATED(region%subRegions)) DEALLOCATE(region%subRegions)
      IF(ASSOCIATED(region%interfaces)) CALL Interfaces_Finalise(region%interfaces,err,error,*999)
      DEALLOCATE(region)
    ENDIF
    
    EXITS("Region_Finalise")
    RETURN
999 ERRORSEXITS("Region_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Region_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises a region.
  SUBROUTINE Region_Initialise(region,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
   
    ENTERS("Region_Initialise",err,error,*998)

    IF(ASSOCIATED(region)) CALL FlagError("Region is already associated.",err,error,*998)
    
    ALLOCATE(region,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate region.",err,error,*999)
    region%userNumber=0
    NULLIFY(region%regions)
    region%regionFinished=.FALSE.
    region%label=""
    NULLIFY(region%coordinateSystem)
    NULLIFY(region%dataPointSets)
    NULLIFY(region%nodes)
    NULLIFY(region%meshes)
    NULLIFY(region%generatedMeshes)
    NULLIFY(region%decomposers)
    NULLIFY(region%fields)
    NULLIFY(region%equationsSets)
    NULLIFY(region%cellMLEnvironments)
    NULLIFY(region%parentRegion)
    region%numberOfSubRegions=0
    NULLIFY(region%subRegions)
    NULLIFY(region%interfaces)
    CALL DataPointSets_Initialise(region,err,error,*999)
    CALL Meshes_Initialise(region,err,error,*999)
    CALL GENERATED_MESHES_INITIALISE(region,err,error,*999)
    CALL Decomposers_Initialise(region,err,error,*999)
    CALL Fields_Initialise(region,err,error,*999)
    CALL EquationsSets_Initialise(region,err,error,*999)
    CALL CELLML_ENVIRONMENTS_INITIALISE(region,err,error,*999)
    CALL Interfaces_Initialise(region,err,error,*999)
    
    EXITS("Region_Initialise")
    RETURN
999 CALL Region_Finalise(region,dummyErr,dummyError,*998)
998 ERRORSEXITS("Region_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Region_Initialise

  !
  !================================================================================================================================
  !

  !>Sets the label of a region. \see OPENCMISS::Iron::cmfe_Region_LabelSet
  SUBROUTINE Region_LabelSetC(region,label,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to set the label for 
    CHARACTER(LEN=*), INTENT(IN) :: label !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Region_LabelSetC",err,error,*999)

    CALL Region_AssertNotFinished(region,err,error,*999)
    
    region%label=label
        
    EXITS("Region_LabelSetC")
    RETURN
999 ERRORSEXITS("Region_LabelSetC",err,error)
    RETURN 1
    
  END SUBROUTINE Region_LabelSetC

  !
  !================================================================================================================================
  !

  !>Sets the label of a region. \see OPENCMISS::Iron::cmfe_Region_LabelSet
  SUBROUTINE Region_LabelSetVS(region,label,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to set the label for 
    TYPE(VARYING_STRING), INTENT(IN) :: label !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Region_LabelSetVS",err,error,*999)

    CALL Region_AssertNotFinished(region,err,error,*999)
    
    region%label=label
    
    EXITS("Region_LabelSetVS")
    RETURN
999 ERRORSEXITS("Region_LabelSetVS",err,error)
    RETURN 1
    
  END SUBROUTINE Region_LabelSetVS

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
          CALL Region_DestroyNumber(regions,regions%worldRegion%subRegions(1)%ptr%userNumber,err,error,*999)
        ENDDO !region
        !Destroy global region and deallocated any memory allocated in the global region
        CALL Region_Finalise(regions%worldRegion,err,error,*999)
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
    
    ENTERS("Regions_Initialise",err,error,*998)

    IF(.NOT.ASSOCIATED(context)) CALL FlagError("Context is not associated.",err,error,*998)
    IF(ASSOCIATED(context%regions)) CALL FlagError("Context regions is already associated.",err,error,*998)

    NULLIFY(coordinateSystems)
    CALL Context_CoordinateSystemsGet(context,coordinateSystems,err,error,*998)
    NULLIFY(worldCoordinateSystem)
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
  
END MODULE RegionRoutines
