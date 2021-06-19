!> \file
!> \author Chris Bradley
!> \brief This module contains all interface conditions routines.
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

!>This module contains all interface conditions routines.
MODULE InterfaceConditionRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE CoordinateSystemAccessRoutines
  USE DomainMappings
  USE EquationsSetAccessRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE InputOutput
  USE InterfaceAccessRoutines
  USE InterfaceConditionAccessRoutines
  USE InterfaceEquationsRoutines
  USE InterfaceEquationsAccessRoutines
  USE InterfaceMappingRoutines
  USE InterfaceMatricesRoutines
  USE InterfaceMatricesAccessRoutines
  USE InterfaceOperatorsRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE MatrixVector
  USE ProfilingRoutines
  USE RegionAccessRoutines
  USE Strings
  USE Timer
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  !Module types

  !Module variables

  !Interfaces

  PUBLIC InterfaceCondition_Assemble

  PUBLIC InterfaceCondition_CreateFinish,InterfaceCondition_CreateStart

  PUBLIC InterfaceCondition_DependentVariableAdd

  PUBLIC InterfaceCondition_Destroy

  PUBLIC InterfaceCondition_EquationsCreateFinish,InterfaceCondition_EquationsCreateStart

  PUBLIC InterfaceCondition_EquationsDestroy
  
  PUBLIC InterfaceCondition_IntegrationTypeSet

  PUBLIC InterfaceCondition_LagrangeFieldCreateFinish,InterfaceCondition_LagrangeFieldCreateStart

  PUBLIC InterfaceCondition_MethodSet

  PUBLIC InterfaceCondition_OperatorSet

  PUBLIC InterfaceCondition_OutputTypeSet
  
  PUBLIC InterfaceCondition_PenaltyFieldCreateFinish,InterfaceCondition_PenaltyFieldCreateStart

  PUBLIC InterfaceConditions_Finalise,InterfaceConditions_Initialise

CONTAINS

  !
  !================================================================================================================================
  !

  !>Assembles the equations for an interface condition.
  SUBROUTINE InterfaceCondition_Assemble(interfaceCondition,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to assemble the equations for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceCondition_Assemble",err,error,*999)
    
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
    
    IF(interfaceCondition%outputType>=INTERFACE_CONDITION_PROGRESS_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Interface condition assemble: ",interfaceCondition%label,err,error,*999)
    ENDIF
    NULLIFY(interfaceEquations)
    CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
    CALL InterfaceEquations_AssertIsFinished(interfaceEquations,err,error,*999)
    SELECT CASE(interfaceCondition%method)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
      CALL InterfaceCondition_AssembleFEM(interfaceCondition,err,error,*999)
    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The interface condition method of "//TRIM(NumberToVString(interfaceCondition%method,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("InterfaceCondition_Assemble")
    RETURN
999 ERRORSEXITS("InterfaceCondition_Assemble",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_Assemble

  !
  !================================================================================================================================
  !
  
  !>Assembles the interface matricesand rhs for using the finite element method.
  SUBROUTINE InterfaceCondition_AssembleFEM(interfaceCondition,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to assemble the equations for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: boundaryStart,elementIdx,elementNumber,ghostFinish,internalFinish,internalStart,numberOfTimes,outputType
    REAL(SP) :: elementUserElapsed,elementSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1),userTime4(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1),systemTime4(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainMappingType), POINTER :: elementsMapping
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices
    TYPE(FieldType), POINTER :: lagrangeField
    
    ENTERS("InterfaceCondition_AssembleFEM",err,error,*999)

    NULLIFY(interfaceEquations)
    CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
    CALL InterfaceEquations_OutputTypeGet(interfaceEquations,outputType,err,error,*999)
    NULLIFY(interfaceMatrices)
    CALL InterfaceEquations_InterfaceMatricesGet(interfaceEquations,interfaceMatrices,err,error,*999)
    IF(outputType>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
    !Initialise the matrices and rhs vector
    CALL InterfaceMatrices_ValueInitialise(interfaceMatrices,0.0_DP,err,error,*999)
    !Assemble the elements
    !Allocate the element matrices 
    CALL InterfaceMatrices_ElementInitialise(interfaceMatrices,err,error,*999)
    NULLIFY(lagrangeField)
    CALL InterfaceCondition_LagrangeFieldGet(interfaceCondition,lagrangeField,err,error,*999)
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(lagrangeField,decomposition,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainMappings)
    CALL Domain_DomainMappingsGet(domain,domainMappings,err,error,*999)
    NULLIFY(elementsMapping)
    CALL DomainMappings_ElementsMappingGet(domainMappings,elementsMapping,err,error,*999)
    CALL DomainMapping_InternalStartGet(elementsMapping,internalStart,err,error,*999)
    CALL DomainMapping_InternalFinishGet(elementsMapping,internalFinish,err,error,*999)
    CALL DomainMapping_BoundaryStartGet(elementsMapping,boundaryStart,err,error,*999)
    CALL DomainMapping_GhostFinishGet(elementsMapping,ghostFinish,err,error,*999)
    
    !Output timing information if required
    IF(outputType>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime2,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime2,err,error,*999)
      userElapsed=userTime2(1)-userTime1(1)
      systemElapsed=systemTime2(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Setup and initialisation",userElapsed,systemElapsed,err,error,*999)
      elementUserElapsed=0.0_SP
      elementSystemElapsed=0.0_SP
    ENDIF
    
    numberOfTimes=0
    !Loop over the internal elements
    DO elementIdx=internalStart,internalFinish
      CALL DomainMapping_NumberGet(elementsMapping,elementIdx,elementNumber,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL InterfaceMatrices_ElementCalculate(interfaceMatrices,elementNumber,err,error,*999)
      CALL InterfaceCondition_FiniteElementCalculate(interfaceCondition,elementNumber,err,error,*999)
      CALL InterfaceMatrices_ElementAdd(interfaceMatrices,err,error,*999)
    ENDDO !elementIdx

    !Output timing information if required
    IF(outputType>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime3,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime3,err,error,*999)
      userElapsed=userTime3(1)-userTime2(1)
      systemElapsed=systemTime3(1)-systemTime2(1)
      elementUserElapsed=userElapsed
      elementSystemElapsed=systemElapsed
      CALL Profiling_TimingsOutput(1,"Internal elements equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    
    !Loop over the boundary and ghost elements
    DO elementIdx=boundaryStart,ghostFinish
      CALL DomainMapping_NumberGet(elementsMapping,elementIdx,elementNumber,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL InterfaceMatrices_ElementCalculate(interfaceMatrices,elementNumber,err,error,*999)
      CALL InterfaceCondition_FiniteElementCalculate(interfaceCondition,elementNumber,err,error,*999)
      CALL InterfaceMatrices_ElementAdd(interfaceMatrices,err,error,*999)
    ENDDO !elementIdx
    
    !Output timing information if required
    IF(outputType>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      elementUserElapsed=elementUserElapsed+userElapsed
      elementSystemElapsed=elementSystemElapsed+systemElapsed
      CALL Profiling_TimingsOutput(1,"Boundary+ghost elements equations assembly",userElapsed,systemElapsed,err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average element equations assembly", &
        & elementUserElapsed/numberOfTimes,elementSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    
    !Finalise the element matrices
    CALL InterfaceMatrices_ElementFinalise(interfaceMatrices,err,error,*999)
    
    !Output timing information if required
    IF(outputType>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations matrices and vector if required
    IF(outputType>=INTERFACE_EQUATIONS_MATRIX_OUTPUT) THEN
      CALL InterfaceMatrices_Output(GENERAL_OUTPUT_TYPE,interfaceMatrices,err,error,*999)
    ENDIF
       
    EXITS("InterfaceCondition_AssembleFEM")
    RETURN
999 ERRORSEXITS("InterfaceCondition_AssembleFEM",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_AssembleFEM

  !
  !==================================================================================================================================
  !

  !>Finishes the process of creating an interface condition. \see OpenCMISS::Iron::cmfe_InterfaceCondition_CreateFinish
  SUBROUTINE InterfaceCondition_CreateFinish(interfaceCondition,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to finish creating
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceConditionMethod,interfaceOperator,meshIdx,meshIdxCount,numberOfComponents,numberOfComponents2, &
      & numberOfCoupledMeshes,numberOfDependentVariables,variableIdx,variableMeshIdx
    INTEGER(INTG), POINTER :: newVariableMeshIndices(:)
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(FieldVariablePtrType), POINTER :: newFieldVariables(:)
    TYPE(InterfaceType), POINTER :: INTERFACE
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("InterfaceCondition_CreateFinish",err,error,*999)

    CALL InterfaceCondition_AssertNotFinished(interfaceCondition,err,error,*999)

    NULLIFY(INTERFACE)
    CALL InterfaceCondition_InterfaceGet(interfaceCondition,INTERFACE,err,error,*999)
    
    !Test various inputs have been set up.
    CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
    SELECT CASE(interfaceConditionMethod)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
      NULLIFY(interfaceDependent)
      CALL InterfaceCondition_InterfaceDependentGet(interfaceCondition,interfaceDependent,err,error,*999)
      !Check the dependent field variables have been set.
      CALL InterfaceDependent_NumberOfDependentVariablesGet(interfaceDependent,numberOfDependentVariables,err,error,*999)
      IF(numberOfDependentVariables<2) THEN
        localError="The number of added dependent variables of "// &
          & TRIM(NumberToVString(numberOfDependentVariables,"*",err,error))//" is invalid. The number must be >= 2."
        CALL FlagError(localError,err,error,*999)
      ENDIF

      !\todo check if interface mesh connectivity basis has same number of gauss points as interface geometric field
 
      !Note There is no need to check that the dependent variables have the same number of components.
      !The user will need to set a fixed BC on the interface dof relating to the field components 
      !not present in each of the coupled bodies, eliminating this dof from the solver matrices
      CALL InterfaceCondition_OperatorGet(interfaceCondition,interfaceOperator,err,error,*999)
      SELECT CASE(interfaceOperator)
      CASE(INTERFACE_CONDITION_FIELD_CONTINUITY_OPERATOR,INTERFACE_CONDITION_FLS_CONTACT_OPERATOR, &
        & INTERFACE_CONDITION_FLS_CONTACT_REPROJECT_OPERATOR,INTERFACE_CONDITION_SOLID_FLUID_OPERATOR)
        !Check that the dependent variables have the same number of components
        NULLIFY(fieldVariable)
        CALL InterfaceDependent_DependentVariableGet(interfaceDependent,1,fieldVariable,err,error,*999)
        CALL FieldVariable_NumberOfComponentsGet(fieldVariable,numberOfComponents,err,error,*999)
        DO variableIdx=2,numberOfDependentVariables
          NULLIFY(fieldVariable)
          CALL InterfaceDependent_DependentVariableGet(interfaceDependent,variableIdx,fieldVariable,err,error,*999)
          CALL FieldVariable_NumberOfComponentsGet(fieldVariable,numberOfComponents2,err,error,*999)
          IF(numberOfComponents/=numberOfComponents2) THEN
            localError="The number of components of "//TRIM(NumberToVString(numberOfComponents2,"*",err,error))// &
              & " for interface dependent variable "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
              & " does not match the number of components of "//TRIM(NumberToVString(numberOfComponents,"*",err,error))// &
              & " for the first interface dependent variable."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDDO !variableIdx 
      CASE(INTERFACE_CONDITION_FIELD_NORMAL_CONTINUITY_OPERATOR)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(INTERFACE_CONDITION_SOLID_FLUID_NORMAL_OPERATOR)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The interface condition operator of "// &
          & TRIM(NumberToVString(interfaceOperator,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT

      !Reorder the dependent variables based on mesh index order
      ALLOCATE(newFieldVariables(numberOfDependentVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new field variables.",err,error,*999)
      ALLOCATE(newVariableMeshIndices(numberOfDependentVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new variable mesh indices.",err,error,*999)
      newVariableMeshIndices=0
      meshIdxCount=0
      CALL Interface_NumberOfCoupledMeshesGet(INTERFACE,numberOfCoupledMeshes,err,error,*999)
      DO meshIdx=1,numberOfCoupledMeshes
        DO variableIdx=1,numberOfDependentVariables
          CALL InterfaceDependent_VariableMeshIndexGet(interfaceDependent,variableIdx,variableMeshIdx,err,error,*999)
          IF(variableMeshIdx==meshIdx) THEN
            meshIdxCount=meshIdxCount+1
            NULLIFY(fieldVariable)
            CALL InterfaceDependent_DependentVariableGet(interfaceDependent,variableIdx,fieldVariable,err,error,*999)
            newFieldVariables(meshIdxCount)%ptr=>fieldVariable
            newVariableMeshIndices(meshIdxCount)=meshIdx
          ENDIF
        ENDDO !variableIdx
      ENDDO !meshIdx
      IF(meshIdxCount/=numberOfDependentVariables) &
        & CALL FlagError("Invalid dependent variable mesh index setup.",err,error,*999)
      IF(ASSOCIATED(interfaceDependent%fieldVariables)) DEALLOCATE(interfaceDependent%fieldVariables)
      IF(ASSOCIATED(interfaceDependent%variableMeshIndices)) DEALLOCATE(interfaceDependent%variableMeshIndices)
      interfaceDependent%fieldVariables=>newFieldVariables
      interfaceDependent%variableMeshIndices=>newVariableMeshIndices
    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The interface condition method of "//TRIM(NumberToVString(interfaceCondition%method,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !Finish the interface condition creation
    interfaceCondition%interfaceConditionFinished=.TRUE.
       
    EXITS("InterfaceCondition_CreateFinish")
    RETURN
999 IF(ASSOCIATED(newFieldVariables)) DEALLOCATE(newFieldVariables)
    IF(ASSOCIATED(newVariableMeshIndices)) DEALLOCATE(newVariableMeshIndices)
    ERRORSEXITS("InterfaceCondition_CreateFinish",err,error)    
    RETURN 1
   
  END SUBROUTINE InterfaceCondition_CreateFinish
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating an interface condition on an interface. \see OpenCMISS::Iron::cmfe_InterfaceCondition_CreateStart
  SUBROUTINE InterfaceCondition_CreateStart(userNumber,interface,geometricField,interfaceCondition,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the interface condition
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to create the interface condition on
    TYPE(FieldType), POINTER :: geometricField !<A pointer to the geometric field for the interface condition.
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<On return, a pointer to the interface condition. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,interfaceConditionIdx,interfaceUserNumber
    TYPE(InterfaceType), POINTER :: geometricInterface
    TYPE(InterfaceConditionType), POINTER :: newInterfaceCondition
    TYPE(InterfaceConditionPtrType), POINTER :: newInterfaceConditions(:)
    TYPE(RegionType), POINTER :: geometricRegion,geometricInterfaceParentRegion,interfaceParentRegion
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("InterfaceCondition_CreateStart",err,error,*997)

    IF(ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is already associated.",err,error,*997)
    IF(.NOT.ASSOCIATED(INTERFACE)) CALL FlagError("Interface is not associated.",err,error,*997)
    CALL Interface_UserNumberGet(INTERFACE,interfaceUserNumber,err,error,*997)
    CALL InterfaceCondition_UserNumberFind(userNumber,INTERFACE,newInterfaceCondition,err,error,*997)
    IF(ASSOCIATED(newInterfaceCondition)) THEN
      localError="Interface condition user number "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " has already been created on interface number "//TRIM(NumberToVString(interfaceUserNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*997)
    ENDIF
    
    CALL Field_AssertIsFinished(geometricField,err,error,*997)
    !Check the geometric field is defined on the interface
    CALL Field_InterfaceCheck(geometricField,INTERFACE,err,error,*997)
    !OK, so crete new interface condition
    NULLIFY(newInterfaceCondition)
    !Initialise the new interface condition
    CALL InterfaceCondition_Initialise(newInterfaceCondition,err,error,*998)
    !Set default interface condition values
    newInterfaceCondition%userNumber=userNumber
    newInterfaceCondition%globalNumber=INTERFACE%interfaceConditions%numberOfInterfaceConditions+1
    newInterfaceCondition%interfaceConditions=>interface%interfaceConditions
    newInterfaceCondition%label="Interface Condition "//TRIM(NumberToVString(userNumber,"*",err,error))
    newInterfaceCondition%INTERFACE=>interface
    !Default attributes
    newInterfaceCondition%geometry%geometricField=>geometricField
    newInterfaceCondition%method=INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD
    newInterfaceCondition%operator=INTERFACE_CONDITION_FIELD_CONTINUITY_OPERATOR
    IF(ASSOCIATED(INTERFACE%pointsConnectivity)) THEN
      newInterfaceCondition%integrationType=INTERFACE_CONDITION_DATA_POINTS_INTEGRATION
    ELSE
      newInterfaceCondition%integrationType=INTERFACE_CONDITION_GAUSS_INTEGRATION
    ENDIF
    CALL InterfaceCondition_DependentInitialise(newInterfaceCondition,err,error,*998)
    !Add new interface condition into list of interface conditions in the interface
    ALLOCATE(newInterfaceConditions(INTERFACE%interfaceConditions%numberOfInterfaceConditions+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new interface conditions.",err,error,*999)
    DO interfaceConditionIdx=1,INTERFACE%interfaceConditions%numberOfInterfaceConditions
      newInterfaceConditions(interfaceConditionIdx)%ptr=> &
        & interface%interfaceConditions%interfaceConditions(interfaceConditionIdx)%ptr
    ENDDO !interfaceConditionIdx
    newInterfaceConditions(INTERFACE%interfaceConditions%numberOfInterfaceConditions+1)%ptr=>newInterfaceCondition
    IF(ASSOCIATED(INTERFACE%interfaceConditions%interfaceConditions)) DEALLOCATE(INTERFACE%interfaceConditions%interfaceConditions)
    INTERFACE%interfaceConditions%interfaceConditions=>newInterfaceConditions
    INTERFACE%interfaceConditions%numberOfInterfaceConditions=INTERFACE%interfaceConditions%numberOfInterfaceConditions+1
    !Return the pointer
    interfaceCondition=>newInterfaceCondition
    
    EXITS("InterfaceCondition_CreateStart")
    RETURN
999 IF(ASSOCIATED(newInterfaceConditions)) DEALLOCATE(newInterfaceConditions)
998 IF(ASSOCIATED(newInterfaceCondition)) CALL InterfaceCondition_Finalise(newInterfaceCondition,dummyErr,dummyError,*997)
997 ERRORSEXITS("InterfaceCondition_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_CreateStart
  
  !
  !================================================================================================================================
  !

  !>Finalise the interface condition dependent field information and deallocate all memory.
  SUBROUTINE InterfaceCondition_DependentFinalise(interfaceDependent,err,error,*)

    !Argument variables
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent !<A pointer to the interface condition dependent field information to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfaceCondition_DependentFinalise",err,error,*999)

    IF(ASSOCIATED(interfaceDependent)) THEN
      IF(ASSOCIATED(interfaceDependent%equationsSets)) DEALLOCATE(interfaceDependent%equationsSets)
      IF(ASSOCIATED(interfaceDependent%fieldVariables)) DEALLOCATE(interfaceDependent%fieldVariables)
      IF(ASSOCIATED(interfaceDependent%variableMeshIndices)) DEALLOCATE(interfaceDependent%variableMeshIndices)
      DEALLOCATE(interfaceDependent)
    ENDIF
       
    EXITS("InterfaceCondition_DependentFinalise")
    RETURN
999 ERRORSEXITS("InterfaceCondition_DependentFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_DependentFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an interface condition dependent field information.
  SUBROUTINE InterfaceCondition_DependentInitialise(interfaceCondition,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<The pointer to the interface condition to initialise to initialise the dependent field information for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("InterfaceCondition_DependentInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*998)
    IF(ASSOCIATED(interfaceCondition%dependent)) &
      & CALL FlagError("Interface condition dependent is already associated.",err,error,*998)
    
    ALLOCATE(interfaceCondition%dependent,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface condition dependent.",err,error,*999)
    interfaceCondition%dependent%interfaceCondition=>interfaceCondition
    interfaceCondition%dependent%numberOfDependentVariables=0
    NULLIFY(interfaceCondition%dependent%equationsSets)
    NULLIFY(interfaceCondition%dependent%fieldVariables)
    NULLIFY(interfaceCondition%dependent%variableMeshIndices)
       
    EXITS("InterfaceCondition_DependentInitialise")
    RETURN
999 CALL InterfaceCondition_DependentFinalise(interfaceCondition%dependent,dummyErr,dummyError,*998)
998 ERRORSEXITS("InterfaceCondition_DependentInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_DependentInitialise

  !
  !================================================================================================================================
  !

  !>Adds an equations set to an interface condition. \see OpenCMISS::Iron::cmfe_InterfaceCondition_DependentVariableAdd
  SUBROUTINE InterfaceCondition_DependentVariableAdd(interfaceCondition,meshIndex,equationsSet,variableType,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to add the dependent variable to
    INTEGER(INTG), INTENT(IN) :: meshIndex !<The mesh index in the interface conditions interface that the dependent variable corresponds to
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set containing the dependent field to add the variable from.
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type of the dependent field to add \see FieldRoutines_VariableTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esUserNumber,fieldUserNumber,interfaceUserNumber,numberOfCoupledMeshes,variableIdx
    INTEGER(INTG), POINTER :: newVariableMeshIndices(:)
    LOGICAL :: foundMeshIndex
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(EquationsSetPtrType), POINTER :: newEquationsSets(:)
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldVariableType), POINTER :: dependentVariable,fieldVariable,interfaceVariable
    TYPE(FieldVariablePtrType), POINTER :: newFieldVariables(:)
    TYPE(InterfaceType), POINTER :: INTERFACE
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent
    TYPE(MeshType), POINTER :: decompositionMesh,dependentMesh,interfaceMesh
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceCondition_DependentVariableAdd",err,error,*999)

    NULLIFY(interfaceDependent)
    CALL InterfaceCondition_InterfaceDependentGet(interfaceCondition,interfaceDependent,err,error,*999)
    CALL EquationsSet_UserNumberGet(equationsSet,esUserNumber,err,error,*999)
    NULLIFY(INTERFACE)
    CALL InterfaceCondition_InterfaceGet(interfaceCondition,INTERFACE,err,error,*999)
    CALL Interface_UserNumberGet(INTERFACE,interfaceUserNumber,err,error,*999)
    CALL Interface_NumberOfCoupledMeshesGet(INTERFACE,numberOfCoupledMeshes,err,error,*999)
    IF(meshIndex<1.OR.meshIndex>numberOfCoupledMeshes) THEN
      localError="The specified mesh index of "//TRIM(NumberToVString(meshIndex,"*",err,error))// &
        & " is invalid. The mesh index must be >= 1 and <= "//TRIM(NumberToVString(numberOfCoupledMeshes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    CALL Field_UserNumberGet(dependentField,fieldUserNumber,err,error,*999)
    NULLIFY(dependentVariable)
    CALL Field_VariableGet(dependentField,variableType,dependentVariable,err,error,*999)
    
    !Check that the field variable hasn't already been added.
    variableIdx=1
    NULLIFY(interfaceVariable)
    DO WHILE(variableIdx<=interfaceDependent%numberOfDependentVariables.AND..NOT.ASSOCIATED(interfaceVariable))
      IF(ASSOCIATED(fieldVariable,interfaceDependent%fieldVariables(variableIdx)%ptr)) THEN
        interfaceVariable=>interfaceDependent%fieldVariables(variableIdx)%ptr
      ELSE
        variableIdx=variableIdx+1
      ENDIF
    ENDDO
    IF(ASSOCIATED(interfaceVariable)) THEN
      !Check if we are dealing with the same mesh index.
      IF(meshIndex/=interfaceDependent%variableMeshIndices(variableIdx)) THEN
        localError="Dependent variable type "//TRIM(NumberToVString(variableType,"*",err,error))// &
          & " of field number "//TRIM(NumberToVString(fieldUserNumber,"*",err,error))// &
          & " has already been added to the interface condition "// &
          & TRIM(NumberToVString(interfaceCondition%userNumber,"*",err,error))//" at position index "// &
          & TRIM(NumberToVString(variableIdx,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      !Check the dependent variable and the mesh index match.
      NULLIFY(interfaceMesh)
      CALL Interface_CoupledMeshGet(INTERFACE,meshIndex,interfaceMesh,err,error,*999)
      NULLIFY(decomposition)
      CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
      NULLIFY(decompositionMesh)
      CALL Decomposition_MeshGet(decomposition,decompositionMesh,err,error,*999)
      IF(.NOT.ASSOCIATED(interfaceMesh,decompositionMesh)) THEN
        localError="The mesh for dependent field number "//TRIM(NumberToVString(fieldUserNumber,"*",err,error))// &
          & " does not match the mesh for interface number "//TRIM(NumberToVString(interfaceUserNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !The meshes match. Check if the dependent variable has already been added for the mesh index.
      foundMeshIndex=.FALSE.
      DO variableIdx=1,interfaceDependent%numberOfDependentVariables
        IF(interfaceDependent%variableMeshIndices(variableIdx)==meshIndex) THEN
          foundMeshIndex=.TRUE.
          EXIT
        ENDIF
      ENDDO !variableIdx
      IF(foundMeshIndex) THEN
        !The mesh index has already been added to replace the dependent variable with the specified variable
        interfaceDependent%fieldVariables(variableIdx)%ptr=>dependentVariable
      ELSE
        !The mesh index has not been found so add a new dependent variable.
        ALLOCATE(newEquationsSets(interfaceDependent%numberOfDependentVariables+1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new equations sets.",err,error,*999)
        ALLOCATE(newFieldVariables(interfaceDependent%numberOfDependentVariables+1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new field variables.",err,error,*999)
        ALLOCATE(newVariableMeshIndices(interfaceDependent%numberOfDependentVariables+1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new variable mesh indices.",err,error,*999)
        DO variableIdx=1,interfaceDependent%numberOfDependentVariables
          newEquationsSets(variableIdx)%ptr=>interfaceDependent%equationsSets(variableIdx)%ptr
          newFieldVariables(variableIdx)%ptr=>interfaceDependent%fieldVariables(variableIdx)%ptr
          newVariableMeshIndices(variableIdx)=interfaceDependent%variableMeshIndices(variableIdx)
        ENDDO !variableIdx
        newEquationsSets(interfaceDependent%numberOfDependentVariables+1)%ptr=>equationsSet
        newFieldVariables(interfaceDependent%numberOfDependentVariables+1)%ptr=>dependentVariable
        newVariableMeshIndices(interfaceDependent%numberOfDependentVariables+1)=meshIndex
        IF(ASSOCIATED(interfaceDependent%equationsSets)) DEALLOCATE(interfaceDependent%equationsSets)
        IF(ASSOCIATED(interfaceDependent%fieldVariables)) DEALLOCATE(interfaceDependent%fieldVariables)
        IF(ASSOCIATED(interfaceDependent%variableMeshIndices)) DEALLOCATE(interfaceDependent%variableMeshIndices)
        interfaceDependent%equationsSets=>newEquationsSets
        interfaceDependent%fieldVariables=>newFieldVariables
        interfaceDependent%variableMeshIndices=>newVariableMeshIndices
        interfaceDependent%numberOfDependentVariables=interfaceDependent%numberOfDependentVariables+1
      ENDIF
    ENDIF
    
    EXITS("InterfaceCondition_DependentVariableAdd")
    RETURN
999 ERRORSEXITS("InterfaceCondition_DependentVariableAdd",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_DependentVariableAdd
  
  !
  !================================================================================================================================
  !

  !>Destroys an interface condition. \see OpenCMISS::Iron::cmfe_InterfaceCondition_Destroy
  SUBROUTINE InterfaceCondition_Destroy(interfaceCondition,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceConditionIdx,interfaceConditionPosition
    TYPE(InterfaceConditionPtrType), POINTER :: newInterfaceConditions(:)
    TYPE(InterfaceConditionsType), POINTER :: interfaceConditions

    NULLIFY(newInterfaceConditions)

    ENTERS("InterfaceCondition_Destroy",err,error,*999)

    IF(ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*998)

    interfaceConditionPosition=interfaceCondition%globalNumber

    !Destroy all the interface condition components
    CALL InterfaceCondition_Finalise(interfaceCondition,err,error,*998)
        
    !Remove the interface condition from the list of interface conditions
    IF(interfaceConditions%numberOfInterfaceConditions>1) THEN
      ALLOCATE(newInterfaceConditions(interfaceConditionS%numberOfInterfaceConditions-1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new interface conditions.",err,error,*999)
      DO interfaceConditionIdx=1,interfaceConditions%numberOfInterfaceConditions
        IF(interfaceConditionIdx<interfaceConditionPosition) THEN
          newInterfaceConditions(interfaceConditionIdx)%ptr=>interfaceConditions%interfaceConditions(interfaceConditionIdx)%ptr
        ELSE IF(interfaceConditionIdx>interfaceConditionPosition) THEN
          interfaceConditions%interfaceConditions(interfaceConditionIdx)%ptr%globalNumber= &
            & interfaceConditions%interfaceConditions(interfaceConditionIdx)%ptr%globalNumber-1
          newInterfaceConditions(interfaceConditionIdx-1)%ptr=>interfaceConditions%interfaceConditions(interfaceConditionIdx)%ptr
        ENDIF
      ENDDO !interfaceConditionIdx
      IF(ASSOCIATED(interfaceConditions%interfaceConditions)) DEALLOCATE(interfaceConditions%interfaceConditions)
      interfaceConditions%interfaceConditions=>newInterfaceConditions
      interfaceConditions%numberOfInterfaceConditions=interfaceConditions%numberOfInterfaceConditions-1
    ELSE
      DEALLOCATE(interfaceConditions%interfaceConditions)
      interfaceConditions%numberOfInterfaceConditions=0
    ENDIF

    EXITS("InterfaceCondition_Destroy")
    RETURN
999 IF(ASSOCIATED(newInterfaceConditions)) DEALLOCATE(newInterfaceConditions)
998 ERRORSEXITS("InterfaceCondition_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_Destroy
  
  !
  !================================================================================================================================
  !

  !>Finish the creation of interface equations for the interface condition. \see OpenCMISS::Iron::cmfe_InterfaceCondition_EquationsCreateFinish
  SUBROUTINE InterfaceCondition_EquationsCreateFinish(interfaceCondition,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to finish the creation of the interface equations for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceConditionMethod,matrixIdx,numberOfDependentVariables,numberOfInterfaceMatrices,sparsityType
    INTEGER(INTG), ALLOCATABLE :: storageTypes(:),structureTypes(:)
    REAL(DP) :: matrixCoefficient
    LOGICAL :: hasTranspose
    LOGICAL, ALLOCATABLE :: matricesTranspose(:)
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("InterfaceCondition_EquationsCreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface conditions is not associated.",err,error,*999)

    CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
    SELECT CASE(interfaceConditionMethod)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
      !Finish the interface equations creation
      NULLIFY(interfaceEquations)
      CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
      CALL InterfaceEquations_AssertNotFinished(interfaceEquations,err,error,*999)
      CALL InterfaceEquations_CreateFinish(interfaceEquations,err,error,*999)
      NULLIFY(interfaceDependent)
      CALL InterfaceCondition_InterfaceDependentGet(interfaceCondition,interfaceDependent,err,error,*999)
      CALL InterfaceDependent_NumberOfDependentVariablesGet(interfaceDependent,numberOfDependentVariables,err,error,*999)
      !Create the interface mapping.
      NULLIFY(interfaceMapping)
      CALL InterfaceMapping_CreateStart(interfaceEquations,interfaceMapping,err,error,*999)
      CALL InterfaceMapping_LagrangeVariableTypeSet(interfaceMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
      SELECT CASE(interfaceConditionMethod)
      CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
        numberOfInterfaceMatrices=numberOfDependentVariables
      CASE(INTERFACE_CONDITION_PENALTY_METHOD)
        numberOfInterfaceMatrices=numberOfDependentVariables+1
      ENDSELECT
      CALL InterfaceMapping_NumberOfMatricesSet(interfaceMapping,numberOfInterfaceMatrices,err,error,*999)
      ALLOCATE(matricesTranspose(numberOfInterfaceMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate matrices transpose.",err,error,*999)
      matricesTranspose=.TRUE.
      SELECT CASE(interfaceConditionMethod)
      CASE(INTERFACE_CONDITION_PENALTY_METHOD)
        !Set the last interface matrix to have no transpose
        matricesTranspose(numberOfInterfaceMatrices)=.FALSE.
      END SELECT
      CALL InterfaceMapping_MatricesTransposeSet(interfaceMapping,matricesTranspose,err,error,*999)
      IF(ALLOCATED(matricesTranspose)) DEALLOCATE(matricesTranspose)
      CALL InterfaceMapping_RHSVariableTypeSet(interfaceMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
      
      CALL InterfaceMapping_CreateFinish(interfaceMapping,err,error,*999)
      !Create the interface matrices
      NULLIFY(interfaceMatrices)
      CALL InterfaceMatrices_CreateStart(interfaceEquations,interfaceMatrices,err,error,*999)
      CALL InterfaceMatrices_NumberOfInterfaceMatricesGet(interfaceMatrices,numberOfInterfaceMatrices,err,error,*999)
      ALLOCATE(storageTypes(numberOfInterfaceMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate storage type.",err,error,*999)
      CALL InterfaceEquations_SparsityTypeGet(interfaceEquations,sparsityType,err,error,*999)
      SELECT CASE(sparsityType)
      CASE(INTERFACE_MATRICES_FULL_MATRICES) 
        storageTypes=MATRIX_BLOCK_STORAGE_TYPE
        CALL InterfaceMatrices_StorageTypeSet(interfaceMatrices,storageTypes,err,error,*999)
      CASE(INTERFACE_MATRICES_SPARSE_MATRICES) 
        ALLOCATE(structureTypes(numberOfInterfaceMatrices),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate structure type.",err,error,*999)
        storageTypes=MATRIX_COMPRESSED_ROW_STORAGE_TYPE
        structureTypes=INTERFACE_MATRIX_FEM_STRUCTURE
        CALL InterfaceMatrices_StorageTypeSet(interfaceMatrices,storageTypes,err,error,*999)
        CALL InterfaceMatrices_StructureTypeSet(interfaceMatrices,structureTypes,err,error,*999)
        IF(ALLOCATED(structureTypes)) DEALLOCATE(structureTypes)
      CASE DEFAULT
        localError="The interface equations sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      IF(ALLOCATED(storageTypes)) DEALLOCATE(storageTypes)
      DO matrixIdx=1,numberOfInterfaceMatrices
        NULLIFY(interfaceMatrix)
        CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,matrixIdx,interfaceMatrix,err,error,*999)
        IF((matrixIdx==1).OR. &
          & (interfaceConditionMethod==INTERFACE_CONDITION_PENALTY_METHOD.AND.matrixIdx==numberOfInterfaceMatrices)) THEN
          matrixCoefficient=1.0_DP
        ELSE
          matrixCoefficient=-1.0_DP
        ENDIF
        CALL InterfaceMatrix_MatrixCoefficientSet(interfaceMatrix,matrixCoefficient,err,error,*999)
        CALL InterfaceMatrix_HasTransposeGet(interfaceMatrix,hasTranspose,err,error,*999)
        IF(hasTranspose) CALL InterfaceMatrix_TransposeMatrixCoefficientSet(interfaceMatrix,matrixCoefficient,err,error,*999)
      ENDDO
      CALL InterfaceMatrices_CreateFinish(interfaceMatrices,err,error,*999)
    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The interface condition method of "//TRIM(NumberToVString(interfaceConditionMethod,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    EXITS("InterfaceCondition_EquationsCreateFinish")
    RETURN
999 IF(ALLOCATED(matricesTranspose)) DEALLOCATE(matricesTranspose)
    IF(ALLOCATED(storageTypes)) DEALLOCATE(storageTypes)
    IF(ALLOCATED(structureTypes)) DEALLOCATE(structureTypes)
    ERRORSEXITS("InterfaceCondition_EquationsCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_EquationsCreateFinish

  !
  !================================================================================================================================
  !

  !>Start the creation of interface equations for the interface condition. \see OpenCMISS::Iron::cmfe_InterfaceCondition_EquationsCreateStart
  !>Default values set for the interfaceEquations's attributes are:
  !>- OUTPUT_TYPE: 0 (INTERFACE_EQUATIONS_NO_OUTPUT)
  !>- SPARSITY_TYPE: 1 (INTERFACE_EQUATIONS_SPARSE_MATRICES)
  SUBROUTINE InterfaceCondition_EquationsCreateStart(interfaceCondition,interfaceEquations,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to create the interface equations for
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<On exit, a pointer to the created interface equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceConditionMethod,variableIdx
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent
    TYPE(InterfaceEquationsInterpolationType), POINTER :: interfaceEquationsInterpolation
    TYPE(InterfaceLagrangeType), POINTER :: interfaceLagrange
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceCondition_EquationsCreateStart",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
    IF(ASSOCIATED(interfaceEquations)) CALL FlagError("Interface equations is already associated.",err,error,*999)
      
    NULLIFY(interfaceEquations)
    CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
    SELECT CASE(interfaceConditionMethod)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
      NULLIFY(interfaceLagrange)
      CALL InterfaceCondition_InterfaceLagrangeGet(interfaceCondition,interfaceLagrange,err,error,*999)
      CALL InterfaceLagrange_AssertIsFinished(interfaceLagrange,err,error,*999)
      NULLIFY(interfaceDependent)
      CALL InterfaceCondition_InterfaceDependentGet(interfaceCondition,interfaceDependent,err,error,*999)
      !Initialise the setup
      CALL InterfaceEquations_CreateStart(interfaceCondition,interfaceEquations,err,error,*999)
      !Set the number of interpolation sets
      NULLIFY(interfaceEquationsInterpolation)
      CALL InterfaceEquations_EquationsInterpolationGet(interfaceEquations,interfaceEquationsInterpolation,err,error,*999)
      CALL InterfaceEquations_InterpolationSetsNumberSet(interfaceEquationsInterpolation%interfaceInterpolation,1,1,1, &
        & err,error,*999)
      DO variableIdx=1,interfaceDependent%numberOfDependentVariables
        CALL InterfaceEquations_InterpolationSetsNumberSet(interfaceEquationsInterpolation%variableInterpolation(variableIdx), &
          & 1,1,0,err,error,*999)
      ENDDO !variableIdx
      !Return the pointer
      interfaceEquations=>interfaceCondition%interfaceEquations
    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The interface condition method of "//TRIM(NumberToVString(interfaceConditionMethod,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("InterfaceCondition_EquationsCreateStart")
    RETURN
999 ERRORSEXITS("InterfaceCondition_EquationsCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_EquationsCreateStart

  !
  !================================================================================================================================
  !

  !>Destroy the interface equations for an interface condition. \see OpenCMISS::Iron::cmfe_InterfaceCondition_EquationsDestroy
  SUBROUTINE InterfaceCondition_EquationsDestroy(interfaceCondition,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface conditions to destroy the interface equations for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfaceCondition_EquationsDestroy",err,error,*999)

    IF(ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
    
    CALL InterfaceEquations_Destroy(interfaceCondition%interfaceEquations,err,error,*999)
        
    EXITS("InterfaceCondition_EquationsDestroy")
    RETURN
999 ERRORSEXITS("InterfaceCondition_EquationsDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_EquationsDestroy

  !
  !================================================================================================================================
  !

  !>Finalise the interface condition and deallocate all memory.
  SUBROUTINE InterfaceCondition_Finalise(interfaceCondition,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfaceCondition_Finalise",err,error,*999)

    IF(ASSOCIATED(interfaceCondition)) THEN
      CALL InterfaceCondition_GeometryFinalise(interfaceCondition%geometry,err,error,*999)
      CALL InterfaceCondition_LagrangeFinalise(interfaceCondition%lagrange,err,error,*999)
      CALL InterfaceCondition_PenaltyFinalise(interfaceCondition%penalty,err,error,*999)
      CALL InterfaceCondition_DependentFinalise(interfaceCondition%dependent,err,error,*999)
      IF(ASSOCIATED(interfaceCondition%interfaceEquations)) &
        & CALL InterfaceEquations_Destroy(interfaceCondition%interfaceEquations,err,error,*999)
      DEALLOCATE(interfaceCondition)
    ENDIF
       
    EXITS("InterfaceCondition_Finalise")
    RETURN
999 ERRORSEXITS("InterfaceCondition_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_Finalise

  !
  !================================================================================================================================
  !

  !>Sets/changes the interface condition integration type 
  SUBROUTINE InterfaceCondition_IntegrationTypeSet(interfaceCondition,interfaceConditionIntegrationType,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to set the operator for
    INTEGER(INTG), INTENT(IN) :: interfaceConditionIntegrationType !<The interface condition integration type to set. \see interfaceConditions_IntegrationType,interfaceConditions 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceCondition_IntegrationTypeSet",err,error,*999)

    CALL InterfaceCondition_AssertNotFinished(interfaceCondition,err,error,*999)
    
    SELECT CASE(interfaceConditionIntegrationType)
    CASE(INTERFACE_CONDITION_GAUSS_INTEGRATION)
      interfaceCondition%integrationType=INTERFACE_CONDITION_GAUSS_INTEGRATION
    CASE(INTERFACE_CONDITION_DATA_POINTS_INTEGRATION)
      interfaceCondition%integrationType=INTERFACE_CONDITION_DATA_POINTS_INTEGRATION
    CASE DEFAULT
      localError="The specified interface condition integration type of "// &
        & TRIM(NumberToVString(interfaceConditionIntegrationType,"*",err,error))//" is not valid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("InterfaceCondition_IntegrationTypeSet")
    RETURN
999 ERRORSEXITS("InterfaceCondition_IntegrationTypeSet",err,ERROR)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_IntegrationTypeSet

  !
  !================================================================================================================================
  !

  !>Finalise the interface condition geometry information and deallocate all memory.
  SUBROUTINE InterfaceCondition_GeometryFinalise(interfaceGeometry,err,error,*)

    !Argument variables
    TYPE(InterfaceGeometryType) :: interfaceGeometry !<The interface condition geometry information to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfaceCondition_GeometryFinalise",err,error,*999)

    NULLIFY(interfaceGeometry%interfaceCondition)
    NULLIFY(interfaceGeometry%geometricField)
       
    EXITS("InterfaceCondition_GeometryFinalise")
    RETURN
999 ERRORSEXITS("InterfaceCondition_GeometryFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_GeometryFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an interface condition geometry information.
  SUBROUTINE InterfaceCondition_GeometryInitialise(interfaceCondition,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<The pointer to the interface condition to initialise to initialise the geometry information for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("InterfaceCondition_GeometryInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*998)
    
    interfaceCondition%geometry%interfaceCondition=>interfaceCondition
    NULLIFY(interfaceCondition%geometry%geometricField)
      
    EXITS("InterfaceCondition_GeometryInitialise")
    RETURN
999 CALL InterfaceCondition_GeometryFinalise(interfaceCondition%geometry,dummyErr,dummyError,*998)
998 ERRORSEXITS("InterfaceCondition_GeometryInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_GeometryInitialise

  !
  !================================================================================================================================
  !

  !>Initialises an interface condition.
  SUBROUTINE InterfaceCondition_Initialise(interfaceCondition,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<The pointer to the interface condition to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("InterfaceCondition_Initialise",err,error,*998)

    IF(ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is already associated.",err,error,*998)
    
    ALLOCATE(interfaceCondition,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface condition.",err,error,*999)
    interfaceCondition%userNumber=0
    interfaceCondition%globalNumber=0
    interfaceCondition%interfaceConditionFinished=.FALSE.
    NULLIFY(interfaceCondition%interfaceConditions)
    interfaceCondition%label=""
    NULLIFY(interfaceCondition%INTERFACE)
    interfaceCondition%outputType=INTERFACE_CONDITION_NO_OUTPUT
    interfaceCondition%method=0
    interfaceCondition%OPERATOR=0
    NULLIFY(interfaceCondition%lagrange)
    NULLIFY(interfaceCondition%penalty)
    NULLIFY(interfaceCondition%dependent)
    NULLIFY(interfaceCondition%interfaceEquations)
    CALL InterfaceCondition_GeometryInitialise(interfaceCondition,err,error,*999)
    NULLIFY(interfaceCondition%boundaryConditions)
       
    EXITS("InterfaceCondition_Initialise")
    RETURN
999 CALL InterfaceCondition_Finalise(interfaceCondition,dummyErr,dummyError,*998)
998 ERRORSEXITS("InterfaceCondition_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_Initialise

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating an interface condition's Lagrange multiplier field \see OpenCMISS::Iron::cmfe_InterfaceCondition_LagrangeFieldCreateFinish
  SUBROUTINE InterfaceCondition_LagrangeFieldCreateFinish(interfaceCondition,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to finish creating the Lagrange field for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: LagrangeFieldUVariableNumberOfComponents,LagrangeFieldDelUDelNVariableNumberOfComponents
    TYPE(InterfaceLagrangeType), POINTER :: interfaceLagrange
    
    ENTERS("InterfaceCondition_LagrangeFieldCreateFinish",err,error,*999)

    NULLIFY(interfaceLagrange)
    CALL InterfaceCondition_InterfaceLagrangeGet(interfaceCondition,interfaceLagrange,err,error,*999)
    CALL InterfaceLagrange_AssertNotFinished(interfaceLagrange,err,error,*999)
    
    !Finish the Lagrange field creation
    IF(interfaceLagrange%lagrangeFieldAutoCreated) CALL Field_CreateFinish(interfaceLagrange%lagrangeField,err,error,*999)
    interfaceLagrange%lagrangeFinished=.TRUE.
       
    EXITS("InterfaceCondition_LagrangeFieldCreateFinish")
    RETURN
999 ERRORS("InterfaceCondition_LagrangeFieldCreateFinish",err,error)
    EXITS("InterfaceCondition_LagrangeFieldCreateFinish")
    RETURN 1
   
  END SUBROUTINE InterfaceCondition_LagrangeFieldCreateFinish
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating the Lagrange multiplyer field for interface condition. \see OpenCMISS::Iron::cmfe_InterfaceCondition_LagrangeFieldCreateStart
  SUBROUTINE InterfaceCondition_LagrangeFieldCreateStart(interfaceCondition,lagrangeFieldUserNumber,lagrangeField,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to create the Lagrange field on
    INTEGER(INTG), INTENT(IN) :: lagrangeFieldUserNumber !<The user specified Lagrange field number
    TYPE(FieldType), POINTER :: lagrangeField !<If associated on entry, a pointer to the user created Lagrange field which has the same user number as the specified Lagrange field user number. If not associated on entry, on exit, a pointer to the created Lagrange field for the interface condition.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,dependentVariableNumber,fieldUserNumber,fieldRegionUserNumber,geometricScalingType, &
      & interfaceUserNumber,interfaceRegionUserNumber,interpolationType,numberOfDimensions
    TYPE(CoordinateSystemType), POINTER :: interfaceCoordinateSystem
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(FieldType), POINTER :: field,geometricField
    TYPE(InterfaceType), POINTER :: interface
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent
    TYPE(RegionType), POINTER :: interfaceRegion,lagrangeFieldRegion
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceCondition_LagrangeFieldCreateStart",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
    IF(ASSOCIATED(interfaceCondition%lagrange)) CALL FlagError("Interface condition Lagrange is already associated.",err,error,*999)
    
    NULLIFY(interfaceDependent)
    CALL InterfaceCondition_InterfaceDependentGet(interfaceCondition,interfaceDependent,err,error,*999)
    NULLIFY(interface)
    CALL InterfaceCondition_InterfaceGet(interfaceCondition,INTERFACE,err,error,*999)
    CALL Interface_UserNumberGet(INTERFACE,interfaceUserNumber,err,error,*999)
    NULLIFY(interfaceRegion)
    CALL Interface_ParentRegionGet(INTERFACE,interfaceRegion,err,error,*999)
    NULLIFY(interfaceCoordinateSystem)
    CALL Interface_CoordinateSystemGet(INTERFACE,interfaceCoordinateSystem,err,error,*999)
    CALL CoordinateSystem_DimensionGet(interfaceCoordinateSystem,numberOfDimensions,err,error,*999)
    CALL InterfaceCondition_GeometricFieldGet(interfaceCondition,geometricField,err,error,*999)
    IF(ASSOCIATED(lagrangeField)) THEN
      !Check the Lagrange field has been finished
      CALL Field_AssertIsFinished(lagrangeField,err,error,*999)
      !Check the user numbers match
      CALL Field_UserNumberGet(lagrangeField,fieldUserNumber,err,error,*999)
      IF(lagrangeFieldUserNumber/=fieldUserNumber) THEN
        localError="The specified Lagrange field user number of "//TRIM(NumberToVString(lagrangeFieldUserNumber,"*",err,error))// &
          & " does not match the user number of the specified Lagrange field of "// &
          & TRIM(NumberToVString(fieldUserNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(lagrangeFieldRegion)
      CALL Field_RegionGet(lagrangeField,lagrangeFieldRegion,err,error,*999)
      !Check the field is defined on the same region as the interface
      CALL Region_UserNumberGet(lagrangeFieldRegion,fieldRegionUserNumber,err,error,*999)
      CALL Region_UserNumberGet(interfaceRegion,interfaceRegionUserNumber,err,error,*999)
      IF(fieldRegionUserNumber/=interfaceRegionUserNumber) THEN
        localError="Invalid region setup. The specified Lagrange field has been created on interface number "// &
          & TRIM(NumberToVString(interfaceUserNumber,"*",err,error))//" in parent region number "// &
          & TRIM(NumberToVString(fieldRegionUserNumber,"*",err,error))// &
          & " and the specified interface has been created in parent region number "// &
          & TRIM(NumberToVString(interfaceRegionUserNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      !Check the user number has not already been used for a field in this region.
      NULLIFY(field)
      CALL Field_UserNumberFind(lagrangeFieldUserNumber,interface,field,err,error,*999)
      IF(ASSOCIATED(field)) THEN
        localError="The specified Lagrange field user number of "// &
          & TRIM(NumberToVString(lagrangeFieldUserNumber,"*",err,error))// &
          & " has already been used to create a field on interface number "// &
          & TRIM(NumberToVString(interfaceUserNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF
    CALL InterfaceCondition_LagrangeInitialise(interfaceCondition,err,error,*999)
    IF(.NOT.ASSOCIATED(lagrangeField)) THEN
      !Create the Lagrange field
      interfaceCondition%lagrange%lagrangeFieldAutoCreated=.TRUE.
      CALL Field_CreateStart(lagrangeFieldUserNumber,interfaceCondition%INTERFACE,interfaceCondition%lagrange%lagrangeField, &
        & err,error,*999)
      CALL Field_LabelSet(interfaceCondition%lagrange%lagrangeField,"Lagrange Multipliers Field",err,error,*999)
      CALL Field_TypeSetAndLock(interfaceCondition%lagrange%lagrangeField,FIELD_GENERAL_TYPE,err,error,*999)
      CALL Field_DependentTypeSetAndLock(interfaceCondition%lagrange%lagrangeField,FIELD_DEPENDENT_TYPE,err,error,*999)
      NULLIFY(geometricDecomposition)
      CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
      CALL Field_DecompositionSetAndLock(interfaceCondition%lagrange%lagrangeField,geometricDecomposition,err,error,*999)
      CALL Field_GeometricFieldSetAndLock(interfaceCondition%lagrange%lagrangeField,geometricField,err,error,*999)
      CALL Field_NumberOfVariablesSetAndLock(interfaceCondition%lagrange%lagrangeField,2,err,error,*999)
      CALL Field_VariableTypesSetAndLock(interfaceCondition%lagrange%lagrangeField,[FIELD_U_VARIABLE_TYPE, &
        & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
      CALL Field_VariableLabelSet(interfaceCondition%lagrange%lagrangeField,FIELD_U_VARIABLE_TYPE,"Lambda",err,error,*999)
      CALL Field_VariableLabelSet(interfaceCondition%lagrange%lagrangeField,FIELD_DELUDELN_VARIABLE_TYPE,"Lambda RHS", &
        & err,error,*999)
      CALL Field_DimensionSetAndLock(interfaceCondition%lagrange%lagrangeField,FIELD_U_VARIABLE_TYPE, &
        & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
      CALL Field_DimensionSetAndLock(interfaceCondition%lagrange%lagrangeField,FIELD_DELUDELN_VARIABLE_TYPE, &
        & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
      CALL Field_DataTypeSetAndLock(interfaceCondition%lagrange%lagrangeField,FIELD_U_VARIABLE_TYPE, &
        & FIELD_DP_TYPE,err,error,*999)
      CALL Field_DataTypeSetAndLock(interfaceCondition%lagrange%lagrangeField,FIELD_DELUDELN_VARIABLE_TYPE, &
        & FIELD_DP_TYPE,err,error,*999)
      IF (interfaceCondition%OPERATOR==INTERFACE_CONDITION_SOLID_FLUID_OPERATOR) THEN
        ! Remove pressure component from number of coupled components
        ! INTERFACE_CONDITION_SOLID_FLUID_OPERATOR might not be used as it is equivalent to
        ! INTERFACE_CONDITION_FIELD_CONTINUITY_OPERATOR if set up correctly
        interfaceCondition%lagrange%numberOfComponents=numberOfDimensions
      ELSE
        !Note that only components present in both the coupled meshes interface dependent fields can be coupled
        !Default the number of component to be the minimum number of components across all the coupled dependent 
        !variables.
        !\todo Check ordering of variable components which are coupled and uncoupled are handled correctly to ensure that
        !coupled variable components don't have to always come before the uncoupled variable components
        interfaceCondition%lagrange%numberOfComponents=0
        DO dependentVariableNumber=1,interfaceDependent%numberOfDependentVariables
          IF(interfaceDependent%fieldVariables(dependentVariableNumber)%ptr%numberOfComponents< &
            & interfaceCondition%lagrange%numberOfComponents) THEN
            interfaceCondition%lagrange%numberOfComponents= &
              & interfaceDependent%fieldVariables(dependentVariableNumber)%ptr%numberOfComponents
          ELSE IF(interfaceCondition%lagrange%numberOfComponents==0) THEN
            interfaceCondition%lagrange%numberOfComponents= &
              & interfaceDependent%fieldVariables(dependentVariableNumber)%ptr%numberOfComponents
          ENDIF
        ENDDO !dependentVariableNumber
      ENDIF
      CALL Field_NumberOfComponentsSet(interfaceCondition%lagrange%lagrangeField,FIELD_U_VARIABLE_TYPE, &
        & interfaceCondition%lagrange%numberOfComponents,err,error,*999)
      CALL Field_NumberOfComponentsSet(interfaceCondition%lagrange%lagrangeField,FIELD_DELUDELN_VARIABLE_TYPE, &
        & interfaceCondition%lagrange%numberOfComponents,err,error,*999)
      DO componentIdx=1,interfaceCondition%lagrange%numberOfComponents
        CALL Field_ComponentInterpolationGet(interfaceDependent%fieldVariables(1)%ptr%field,FIELD_U_VARIABLE_TYPE, &
          & componentIdx,interpolationType,err,error,*999)
        CALL Field_ComponentInterpolationSet(interfaceCondition%lagrange%lagrangeField, &
          & FIELD_U_VARIABLE_TYPE,componentIdx,interpolationType,err,error,*999)
        CALL Field_ComponentInterpolationSet(interfaceCondition%lagrange%lagrangeField, &
          & FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,interpolationType,err,error,*999)
      ENDDO !componentIdx
      CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
      CALL Field_ScalingTypeSet(interfaceCondition%lagrange%lagrangeField,geometricScalingType,err,error,*999)
    ELSE
      !Check the Lagrange field
      CALL FlagError("Not implemented.",err,error,*999)
    ENDIF
    !Set pointers
    IF(interfaceCondition%lagrange%lagrangeFieldAutoCreated) THEN
      lagrangeField=>interfaceCondition%lagrange%lagrangeField
    ELSE
      interfaceCondition%lagrange%lagrangeField=>lagrangeField
    ENDIF
    
    EXITS("InterfaceCondition_LagrangeFieldCreateStart")
    RETURN
999 ERRORSEXITS("InterfaceCondition_LagrangeFieldCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_LagrangeFieldCreateStart
  
  !
  !================================================================================================================================
  !

  !>Finalise the interface condition Lagrange information and deallocate all memory.
  SUBROUTINE InterfaceCondition_LagrangeFinalise(interfaceLagrange,err,error,*)

    !Argument variables
    TYPE(InterfaceLagrangeType), POINTER :: interfaceLagrange !<A pointer to the interface condition Lagrange information to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfaceCondition_LagrangeFinalise",err,error,*999)

    IF(ASSOCIATED(interfaceLagrange)) THEN
      DEALLOCATE(interfaceLagrange)
    ENDIF
       
    EXITS("InterfaceCondition_LagrangeFinalise")
    RETURN
999 ERRORSEXITS("InterfaceCondition_LagrangeFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_LagrangeFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an interface condition Lagrange information.
  SUBROUTINE InterfaceCondition_LagrangeInitialise(interfaceCondition,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<The pointer to the interface condition to initialise to initialise the Lagrange information for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("InterfaceCondition_LagrangeInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*998)
    IF(ASSOCIATED(interfaceCondition%lagrange))  &
      & CALL FlagError("Interface condition Lagrange is already associated.",err,error,*998)
      
    ALLOCATE(interfaceCondition%lagrange,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface condition Lagrange.",err,error,*999)
    interfaceCondition%lagrange%interfaceCondition=>interfaceCondition
    interfaceCondition%lagrange%lagrangeFinished=.FALSE.
    interfaceCondition%lagrange%lagrangeFieldAutoCreated=.FALSE.
    NULLIFY(interfaceCondition%lagrange%lagrangeField)
    interfaceCondition%lagrange%numberOfComponents=0
       
    EXITS("InterfaceCondition_LagrangeInitialise")
    RETURN
999 CALL InterfaceCondition_LagrangeFinalise(interfaceCondition%lagrange,dummyErr,dummyError,*998)
998 ERRORSEXITS("InterfaceCondition_LagrangeInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_LagrangeInitialise

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating an interface condition's penalty field'. \see OpenCMISS::Iron::cmfe_InterfaceCondition_PenaltyConditionCreateFinish
  SUBROUTINE InterfaceCondition_PenaltyFieldCreateFinish(interfaceCondition,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to finish creating the penalty field for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(InterfacePenaltyType), POINTER :: interfacePenalty
    
    ENTERS("InterfaceCondition_PenaltyFieldCreateFinish",err,error,*999)

    CALL InterfaceCondition_InterfacePenaltyGet(interfaceCondition,interfacePenalty,err,error,*999)
    CALL InterfacePenalty_AssertNotFinished(interfacePenalty,err,error,*999)
    
    !Finish the penalty field creation
    IF(interfacePenalty%penaltyFieldAutoCreated) CALL Field_CreateFinish(interfacePenalty%penaltyField,err,error,*999)
    interfacePenalty%penaltyFinished=.TRUE.
       
    EXITS("InterfaceCondition_PenaltyFieldCreateFinish")
    RETURN
999 ERRORSEXITS("InterfaceCondition_PenaltyFieldCreateFinish",err,error)
    RETURN 1
   
  END SUBROUTINE InterfaceCondition_PenaltyFieldCreateFinish

  !
  !================================================================================================================================
  !

  !>Starts the process of creating the penalty field for interface condition. \see OpenCMISS::Iron::cmfe_InterfaceCondition_PenaltyFieldCreateStart
  SUBROUTINE InterfaceCondition_PenaltyFieldCreateStart(interfaceCondition,penaltyFieldUserNumber,penaltyField,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to create the penalty field on
    INTEGER(INTG), INTENT(IN) :: penaltyFieldUserNumber !<The user specified penalty field number
    TYPE(FieldType), POINTER :: penaltyField !<If associated on entry, a pointer to the user created penalty field which has the same user number as the specified penalty field user number. If not associated on entry, on exit, a pointer to the created penalty field for the interface condition.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,fieldUserNumber,fieldRegionUserNumber,geometricScalingType,interfaceUserNumber, &
      & numberOfComponents,parentRegionUserNumber
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(FieldType), POINTER :: field,geometricField
    TYPE(FieldVariableType), POINTER :: dependentVariable
    TYPE(InterfaceType), POINTER :: INTERFACE
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent
    TYPE(InterfacePenaltyType), POINTER :: interfacePenalty
    TYPE(RegionType), POINTER :: interfaceRegion,parentRegion,penaltyFieldRegion
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceCondition_PenaltyFieldCreateStart",err,error,*999)

    NULLIFY(interfacePenalty)
    CALL InterfaceCondition_InterfacePenaltyGet(interfaceCondition,interfacePenalty,err,error,*999)
    NULLIFY(interfaceDependent)
    CALL InterfaceCondition_InterfaceDependentGet(interfaceCondition,interfaceDependent,err,error,*999)
    NULLIFY(INTERFACE)
    CALL InterfaceCondition_InterfaceGet(interfaceCondition,INTERFACE,err,error,*999)
    CALL Interface_UserNumberGet(INTERFACE,interfaceUserNumber,err,error,*999)
    NULLIFY(interfaceRegion)
    CALL Interface_ParentRegionGet(INTERFACE,parentRegion,err,error,*999)
    CALL Region_UserNumberGet(parentRegion,parentRegionUserNumber,err,error,*999)
    NULLIFY(geometricField)
    CALL InterfaceCondition_GeometricFieldGet(interfaceCondition,geometricField,err,error,*999)
    IF(ASSOCIATED(penaltyField)) THEN
      !Check the penalty field has been finished
      CALL Field_AssertIsFinished(penaltyField,err,error,*999)
      !Check the user numbers match
      CALL Field_UserNumberGet(penaltyField,fieldUserNumber,err,error,*999)
      IF(penaltyFieldUserNumber/=fieldUserNumber) THEN
        localError="The specified penalty field user number of "//TRIM(NumberToVString(penaltyFieldUserNumber,"*",err,error))// &
          & " does not match the user number of the specified penalty field of "// &
          & TRIM(NumberToVString(fieldUserNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Check that the penalty field is defined on the same interface as the interface condition
      CALL Field_InterfaceCheck(penaltyField,INTERFACE,err,error,*999)
      !Check the user number has not already been used for a field in this region.
      NULLIFY(field)
      CALL Field_UserNumberFind(penaltyFieldUserNumber,INTERFACE,field,err,error,*999)
      IF(ASSOCIATED(field)) THEN
        localError="The specified penalty field user number of "//TRIM(NumberToVString(penaltyFieldUserNumber,"*",err,error))// &
          & " has already been used to create a field on interface number "// &
          & TRIM(NumberToVString(interfaceUserNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF
    CALL InterfaceCondition_PenaltyInitialise(interfaceCondition,err,error,*999)
    IF(.NOT.ASSOCIATED(penaltyField)) THEN
      !Create the penalty field
      interfacePenalty%penaltyFieldAutoCreated=.TRUE.
      CALL Field_CreateStart(penaltyFieldUserNumber,interface,interfacePenalty%penaltyField,err,error,*999)
      CALL Field_LabelSet(interfacePenalty%penaltyField,"Penalty Field",err,error,*999)
      CALL Field_TypeSetAndLock(interfacePenalty%penaltyField,FIELD_GENERAL_TYPE,err,error,*999)
      CALL Field_DependentTypeSetAndLock(interfacePenalty%penaltyField,FIELD_DEPENDENT_TYPE,  err,error,*999)
      NULLIFY(geometricDecomposition)
      CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
      CALL Field_DecompositionSetAndLock(interfacePenalty%penaltyField,geometricDecomposition,err,error,*999)
      CALL Field_GeometricFieldSetAndLock(interfacePenalty%penaltyField,geometricField,err,error,*999)
      CALL Field_NumberOfVariablesSetAndLock(interfacePenalty%penaltyField,1,err,error,*999)
      CALL Field_VariableTypesSetAndLock(interfacePenalty%penaltyField,FIELD_U_VARIABLE_TYPE,err,error,*999)
      CALL Field_VariableLabelSet(interfacePenalty%penaltyField,FIELD_U_VARIABLE_TYPE,"Alpha",err,error,*999)
      CALL Field_DimensionSetAndLock(interfacePenalty%penaltyField,FIELD_U_VARIABLE_TYPE, &
        & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
      CALL Field_DataTypeSetAndLock(interfacePenalty%penaltyField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
      IF(interfaceCondition%OPERATOR==INTERFACE_CONDITION_FLS_CONTACT_OPERATOR .OR. &
        & interfaceCondition%OPERATOR==INTERFACE_CONDITION_FLS_CONTACT_REPROJECT_OPERATOR) THEN
        !Default 1 component for the contact lagrange variable in a frictionless contact problem
        CALL Field_NumberOfComponentsSet(interfacePenalty%penaltyField,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
        CALL Field_ComponentInterpolationSet(interfacePenalty%penaltyField,FIELD_U_VARIABLE_TYPE,1, &
          & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
      ELSE
        !Default the number of component to the first variable of the interface dependent field's number of components,
        NULLIFY(dependentVariable)
        CALL InterfaceDependent_DependentVariableGet(interfaceDependent,1,dependentVariable,err,error,*999)
        CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfComponents,err,error,*999)
        CALL Field_NumberOfComponentsSet(interfacePenalty%penaltyField,FIELD_U_VARIABLE_TYPE,numberOfComponents, &
          & err,error,*999)
        DO componentIdx=1,numberOfComponents
          CALL Field_ComponentInterpolationSetAndLock(interfacePenalty%penaltyField,FIELD_U_VARIABLE_TYPE,componentIdx, &
            & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
        ENDDO !componentIdx
      ENDIF
      CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
      CALL Field_ScalingTypeSet(interfacePenalty%penaltyField,geometricScalingType,err,error,*999)
    ELSE
      !Check the penalty field
      CALL FlagError("Not implemented.",err,error,*999)
    ENDIF
    !Set pointers
    IF(interfacePenalty%penaltyFieldAutoCreated) THEN
      penaltyField=>interfacePenalty%penaltyField
    ELSE
      interfacePenalty%penaltyField=>penaltyField
    ENDIF
    
    EXITS("InterfaceCondition_PenaltyFieldCreateStart")
    RETURN
999 ERRORSEXITS("InterfaceCondition_PenaltyFieldCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_PenaltyFieldCreateStart
  
  !
  !================================================================================================================================
  !

  !>Finalise the interface condition penalty information and deallocate all memory.
  SUBROUTINE InterfaceCondition_PenaltyFinalise(interfacePenalty,err,error,*)

    !Argument variables
    TYPE(InterfacePenaltyType), POINTER :: interfacePenalty !<A pointer to the interface condition penalty information to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfaceCondition_PenaltyFinalise",err,error,*999)

    IF(ASSOCIATED(interfacePenalty)) THEN
      DEALLOCATE(interfacePenalty)
    ENDIF
       
    EXITS("InterfaceCondition_PenaltyFinalise")
    RETURN
999 ERRORSEXITS("InterfaceCondition_PenaltyFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_PenaltyFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an interface condition penalty information.
  SUBROUTINE InterfaceCondition_PenaltyInitialise(interfaceCondition,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<The pointer to the interface condition to initialise to initialise the penalty information for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("InterfaceCondition_PenaltyInitialise",err,error,*998)

    IF(ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*998)
    IF(ASSOCIATED(interfaceCondition%penalty)) CALL FlagError("Interface condition penalty is already associated.",err,error,*998)
      
    ALLOCATE(interfaceCondition%penalty,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface condition penalty.",err,error,*999)
    interfaceCondition%penalty%interfaceCondition=>interfaceCondition
    interfaceCondition%penalty%penaltyFinished=.FALSE.
    interfaceCondition%penalty%penaltyFieldAutoCreated=.FALSE.
    NULLIFY(interfaceCondition%penalty%penaltyField)
       
    EXITS("InterfaceCondition_PenaltyInitialise")
    RETURN
999 CALL InterfaceCondition_PenaltyFinalise(interfaceCondition%penalty,dummyErr,dummyError,*998)
998 ERRORSEXITS("InterfaceCondition_PenaltyInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_PenaltyInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the interface condition method \see OpenCMISS::Iron::cmfe_InterfaceCondition_MethodSet
  SUBROUTINE InterfaceCondition_MethodSet(interfaceCondition,interfaceConditionMethod,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to set the method for
    INTEGER(INTG), INTENT(IN) :: interfaceConditionMethod !<The interface condition method to set. \see interfaceConditions_Methods,interfaceConditions
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceCondition_MethodSet",err,error,*999)

    CALL InterfaceCondition_AssertNotFinished(interfaceCondition,err,error,*999)
    
    SELECT CASE(interfaceConditionMethod)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      interfaceCondition%method=INTERFACE_CONDITION_POINT_TO_POINT_METHOD
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
      interfaceCondition%method=INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD
    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      interfaceCondition%method=INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD
    CASE(INTERFACE_CONDITION_PENALTY_METHOD)
      interfaceCondition%method=INTERFACE_CONDITION_PENALTY_METHOD
    CASE DEFAULT
      localError="The specified interface condition method of "//TRIM(NumberToVString(interfaceConditionMethod,"*",err,error))// &
        & " is not valid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("InterfaceCondition_MethodSet")
    RETURN
999 ERRORSEXITS("InterfaceCondition_MethodSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_MethodSet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the interface condition operator \see OpenCMISS::Iron::cmfe_InterfaceCondition_OperatorSet
  SUBROUTINE InterfaceCondition_OperatorSet(interfaceCondition,interfaceConditionOperator,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to set the operator for
    INTEGER(INTG), INTENT(IN) :: interfaceConditionOperator !<The interface condition operator to set. \see interfaceConditions_Operators,interfaceConditions 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceCondition_OperatorSet",err,error,*999)

    CALL InterfaceCondition_AssertNotFinished(interfaceCondition,err,error,*999)
    
    SELECT CASE(interfaceConditionOperator)
    CASE(INTERFACE_CONDITION_FIELD_CONTINUITY_OPERATOR)
      interfaceCondition%OPERATOR=INTERFACE_CONDITION_FIELD_CONTINUITY_OPERATOR
    CASE(INTERFACE_CONDITION_FIELD_NORMAL_CONTINUITY_OPERATOR)
      interfaceCondition%OPERATOR=INTERFACE_CONDITION_FIELD_NORMAL_CONTINUITY_OPERATOR
    CASE(INTERFACE_CONDITION_FLS_CONTACT_OPERATOR)
      interfaceCondition%OPERATOR=INTERFACE_CONDITION_FLS_CONTACT_OPERATOR
    CASE(INTERFACE_CONDITION_FLS_CONTACT_REPROJECT_OPERATOR)
      interfaceCondition%OPERATOR=INTERFACE_CONDITION_FLS_CONTACT_REPROJECT_OPERATOR
    CASE(INTERFACE_CONDITION_SOLID_FLUID_OPERATOR)
      interfaceCondition%OPERATOR=INTERFACE_CONDITION_SOLID_FLUID_OPERATOR
    CASE(INTERFACE_CONDITION_SOLID_FLUID_NORMAL_OPERATOR)
      interfaceCondition%OPERATOR=INTERFACE_CONDITION_SOLID_FLUID_NORMAL_OPERATOR
    CASE DEFAULT
      localError="The specified interface condition operator of "// &
        & TRIM(NumberToVString(interfaceConditionOperator,"*",err,error))//" is not valid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("InterfaceCondition_OperatorSet")
    RETURN
999 ERRORSEXITS("InterfaceCondition_OperatorSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_OperatorSet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the output type for an interface condition.
  SUBROUTINE InterfaceCondition_OutputTypeSet(interfaceCondition,outputType,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to set the output type for
    INTEGER(INTG), INTENT(IN) :: outputType !<The output type to set \see InterfaceCondition_OutputTypes,InterfaceCondition
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceCondition_OutputTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)

    SELECT CASE(outputType)
    CASE(INTERFACE_CONDITION_NO_OUTPUT)
      interfaceCondition%outputType=INTERFACE_CONDITION_NO_OUTPUT
    CASE(INTERFACE_CONDITION_PROGRESS_OUTPUT)
      interfaceCondition%outputType=INTERFACE_CONDITION_PROGRESS_OUTPUT
    CASE DEFAULT
      localError="The specified output type of "//TRIM(NumberToVString(outputType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("InterfaceCondition_OutputTypeSet")
    RETURN
999 ERRORSEXITS("InterfaceCondition_OutputTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_OutputTypeSet

  !
  !================================================================================================================================
  !

  !>Evaluates the residual for an interface condition.
  SUBROUTINE InterfaceCondition_ResidualEvaluate(interfaceCondition,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceConditionMethod
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceCondition_ResidualEvaluate",err,error,*999)

    NULLIFY(interfaceEquations)
    CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
    CALL InterfaceEquations_AssertIsFinished(interfaceEquations,err,error,*999)

    CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
    SELECT CASE(interfaceConditionMethod)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
      CALL InterfaceCondition_ResidualEvaluateFEM(interfaceCondition,err,error,*999)
    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The interface condition method of "//TRIM(NumberToVString(interfaceConditionMethod,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("InterfaceCondition_ResidualEvaluate")
    RETURN
999 ERRORSEXITS("InterfaceCondition_ResidualEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_ResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the residual for an interface condition using the finite element method
  SUBROUTINE InterfaceCondition_ResidualEvaluateFEM(interfaceCondition,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: boundaryStart,elementIdx,elementNumber,ghostFinish,internalFinish,internalStart,numberOfTimes,outputType
    REAL(SP) :: elementUserElapsed,elementSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1),userTime4(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1),systemTime4(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainMappingType), POINTER :: elementsMapping
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(FieldType), POINTER :: lagrangeField
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceLagrangeType), POINTER :: interfaceLagrange
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices
 
    ENTERS("InterfaceCondition_ResidualEvaluateFEM",err,error,*999)

    NULLIFY(interfaceLagrange)
    CALL InterfaceCondition_InterfaceLagrangeGet(interfaceCondition,interfaceLagrange,err,error,*999)
    NULLIFY(lagrangeField)
    CALL InterfaceLagrange_LagrangeFieldGet(interfaceLagrange,lagrangeField,err,error,*999)
    NULLIFY(interfaceEquations)
    CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
    CALL InterfaceEquations_OutputTypeGet(interfaceEquations,outputType,err,error,*999)
    NULLIFY(interfaceMatrices)
    CALL InterfaceEquations_InterfaceMatricesGet(interfaceEquations,interfaceMatrices,err,error,*999)
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(lagrangeField,decomposition,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainMappings)
    CALL Domain_DomainMappingsGet(domain,domainMappings,err,error,*999)
    NULLIFY(elementsMapping)
    CALL DomainMappings_ElementsMappingGet(domainMappings,elementsMapping,err,error,*999)
    CALL DomainMapping_InternalStartGet(elementsMapping,internalStart,err,error,*999)
    CALL DomainMapping_InternalFinishGet(elementsMapping,internalFinish,err,error,*999)
    CALL DomainMapping_BoundaryStartGet(elementsMapping,boundaryStart,err,error,*999)
    CALL DomainMapping_GhostFinishGet(elementsMapping,ghostFinish,err,error,*999)
    
    IF(outputType>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
!!Do we need to transfer parameter sets???
    !Initialise the matrices and rhs vector
    CALL InterfaceMatrices_ValueInitialise(interfaceMatrices,0.0_DP,err,error,*999)
    !Assemble the elements
    !Allocate the element matrices 
    CALL InterfaceMatrices_ElementInitialise(interfaceMatrices,err,error,*999)
    
    !Output timing information if required
    IF(interfaceEquations%outputType>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime2,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime2,err,error,*999)
      userElapsed=userTime2(1)-userTime1(1)
      systemElapsed=systemTime2(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Setup and initialisation",userElapsed,systemElapsed,err,error,*999)
      elementUserElapsed=0.0_SP
      elementSystemElapsed=0.0_SP
    ENDIF
    numberOfTimes=0
    
    !Loop over the internal elements
    DO elementIdx=internalStart,internalFinish
      CALL DomainMapping_NumberGet(elementsMapping,elementIdx,elementNumber,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL InterfaceMatrices_ElementCalculate(interfaceMatrices,elementNumber,err,error,*999)
      CALL InterfaceCondition_FiniteElementCalculate(interfaceCondition,elementNumber,err,error,*999)
      CALL InterfaceMatrices_ElementAdd(interfaceMatrices,err,error,*999)
    ENDDO !elementIdx
    
    !Output timing information if required
    IF(outputType>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime3,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime3,err,error,*999)
      userElapsed=userTime3(1)-userTime2(1)
      systemElapsed=systemTime3(1)-systemTime2(1)
      elementUserElapsed=userElapsed
      elementSystemElapsed=systemElapsed
      CALL Profiling_TimingsOutput(1,"Internal elements equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    
    !Loop over the boundary and ghost elements
    DO elementIdx=boundaryStart,ghostFinish
      CALL DomainMapping_NumberGet(elementsMapping,elementIdx,elementNumber,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL InterfaceMatrices_ElementCalculate(interfaceMatrices,elementNumber,err,error,*999)
      CALL InterfaceCondition_FiniteElementCalculate(interfaceCondition,elementNumber,err,error,*999)
      CALL InterfaceMatrices_ElementAdd(interfaceMatrices,err,error,*999)
    ENDDO !elementIdx
    
    !Output timing information if required
    IF(interfaceEquations%outputType>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      elementUserElapsed=elementUserElapsed+userElapsed
      elementSystemElapsed=elementSystemElapsed+systemElapsed
      CALL Profiling_TimingsOutput(1,"Boundary+ghost elements equations assembly",userElapsed,systemElapsed, &
        & err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average element equations assembly", &
        & elementUserElapsed/numberOfTimes,elementSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    
    !Finalise the element matrices
    CALL InterfaceMatrices_ElementFinalise(interfaceMatrices,err,error,*999)
    
    !Output timing information if required
    IF(interfaceEquations%outputType>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations matrices and RHS vector if required
    IF(interfaceEquations%outputType>=INTERFACE_EQUATIONS_MATRIX_OUTPUT) THEN
      CALL InterfaceMatrices_Output(GENERAL_OUTPUT_TYPE,interfaceMatrices,err,error,*999)
    ENDIF
       
    EXITS("InterfaceCondition_ResidualEvaluateFEM")
    RETURN
999 ERRORSEXITS("InterfaceCondition_ResidualEvaluateFEM",err,ERROR)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_ResidualEvaluateFEM
  
  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries for the given element number for a finite element interface equations.
  SUBROUTINE InterfaceCondition_FiniteElementCalculate(interfaceCondition,interfaceElementNumber,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition
    INTEGER(INTG), INTENT(IN) :: interfaceElementNumber !<The element number to calcualte
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceConditionOperator,interfaceMatrixIdx,numberOfInterfaceMatrices,outputType
    LOGICAL :: updateMatrix
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("InterfaceCondition_FiniteElementCalculate",err,error,*999)

    NULLIFY(interfaceEquations)
    CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
    CALL InterfaceEquations_OutputTypeGet(interfaceEquations,outputType,err,error,*999)
    
    CALL InterfaceCondition_OperatorGet(interfaceCondition,interfaceConditionOperator,err,error,*999)    
    SELECT CASE(interfaceConditionOperator)
    CASE(INTERFACE_CONDITION_FIELD_CONTINUITY_OPERATOR)
      CALL FieldContinuity_FiniteElementCalculate(interfaceCondition,interfaceElementNumber,err,error,*999)
    CASE(INTERFACE_CONDITION_FIELD_NORMAL_CONTINUITY_OPERATOR)
      CALL FlagError("Not implemented!",err,error,*999)
    CASE(INTERFACE_CONDITION_FLS_CONTACT_OPERATOR,INTERFACE_CONDITION_FLS_CONTACT_REPROJECT_OPERATOR)
      CALL FrictionlessContact_FiniteElementCalculate(interfaceCondition,interfaceElementNumber,err,error,*999)
    CASE(INTERFACE_CONDITION_SOLID_FLUID_OPERATOR)
      CALL SolidFluidOperator_FiniteElementCalculate(interfaceCondition,interfaceElementNumber,err,error,*999)
    CASE(INTERFACE_CONDITION_SOLID_FLUID_NORMAL_OPERATOR)
      CALL FlagError("Not implemented!",err,error,*999)
    CASE DEFAULT
      localError="The interface condition operator of "//TRIM(NumberToVString(interfaceConditionOperator,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    IF(outputType>=INTERFACE_EQUATIONS_ELEMENT_MATRIX_OUTPUT) THEN
      NULLIFY(interfaceMatrices)
      CALL InterfaceEquations_InterfaceMatricesGet(interfaceEquations,interfaceMatrices,err,error,*999)
      CALL InterfaceMatrices_NumberOfInterfaceMatricesGet(interfaceMatrices,numberOfInterfaceMatrices,err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"Finite element interface matrices:",err,error,*999)          
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Element number = ",interfaceElementNumber,err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of interface matrices = ",numberOfInterfaceMatrices,err,error,*999)
      DO interfaceMatrixIdx=1,numberOfInterfaceMatrices
        NULLIFY(interfaceMatrix)
        CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,interfaceMatrixIdx,interfaceMatrix,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Interface matrix : ",interfaceMatrixIdx,err,error,*999)
        CALL InterfaceMatrix_UpdateMatrixGet(interfaceMatrix,updateMatrix,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update matrix = ",updateMatrix,err,error,*999)
        IF(updateMatrix) CALL InterfaceMatrix_ElementMatrixOutput(GENERAL_OUTPUT_TYPE,interfaceMatrix,err,error,*999)
      ENDDO !interfaceMatrixIdx
    ENDIF
       
    EXITS("InterfaceCondition_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("InterfaceCondition_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Finalises an interface conditions and deallocates all memory.
  SUBROUTINE InterfaceConditions_Finalise(interfaceConditions,err,error,*) 

    !Argument variables
    TYPE(InterfaceConditionsType), POINTER :: interfaceConditions !<A pointer to the interface conditions to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    
    ENTERS("InterfaceConditions_Finalise",err,error,*999)
    
    IF(ASSOCIATED(interfaceConditions)) THEN
      DO WHILE(interfaceConditions%numberOfInterfaceConditions>0)
        interfaceCondition=>interfaceConditions%interfaceConditions(1)%ptr
        CALL InterfaceCondition_Destroy(interfaceCondition,err,error,*999)
      ENDDO
      IF(ASSOCIATED(interfaceConditions%interfaceConditions)) DEALLOCATE(interfaceConditions%interfaceConditions)
      DEALLOCATE(interfaceConditions)
    ENDIF
    
    EXITS("InterfaceConditions_Finalise")
    RETURN
999 ERRORSEXITS("InterfaceConditions_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceConditions_Finalise

  !
  !================================================================================================================================
  !
  
  !>Initialises an interface conditions for an interface.
  SUBROUTINE InterfaceConditions_Initialise(interface,err,error,*) 

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to initialise the conditions for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,userNumber
    TYPE(VARYING_STRING) :: dummyError,localError
     
    ENTERS("InterfaceConditions_Initialise",err,error,*998)

    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*998)
    IF(ASSOCIATED(interface%interfaceConditions)) THEN
      CALL Interface_UserNumberGet(interface,userNumber,err,error,*998)
      localError="Interface conditions is already associated for interface number "// &
        & TRIM(NumberToVString(userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    
    ALLOCATE(interface%interfaceConditions,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface interface conditions.",err,error,*999)
    interface%interfaceConditions%interface=>interface
    interface%interfaceConditions%numberOfInterfaceConditions=0
    NULLIFY(interface%interfaceConditions%interfaceConditions)
    
    EXITS("InterfaceConditions_Initialise")
    RETURN
999 CALL InterfaceConditions_Finalise(interface%interfaceConditions,dummyErr,dummyError,*998)
998 ERRORSEXITS("InterfaceConditions_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceConditions_Initialise

  !
  !================================================================================================================================
  !

END MODULE InterfaceConditionRoutines
