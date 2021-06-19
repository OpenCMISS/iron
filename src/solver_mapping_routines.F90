!> \file
!> \author Chris Bradley
!> \brief This module handles all solver mapping routines.
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

!> This module handles all solver mapping routines.
MODULE SolverMappingRoutines

  USE BaseRoutines
  USE BoundaryConditionsRoutines
  USE BoundaryConditionAccessRoutines
  USE ComputationRoutines
  USE ComputationAccessRoutines
  USE DistributedMatrixVector
  USE DomainMappings
  USE EquationsAccessRoutines
  USE EquationsMappingAccessRoutines
  USE EquationsSetAccessRoutines
  USE FieldAccessRoutines
  USE InputOutput
  USE InterfaceAccessRoutines
  USE InterfaceConditionAccessRoutines
  USE InterfaceEquationsAccessRoutines
  USE InterfaceMappingAccessRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE Lists
  USE MatrixVector
  USE ProblemAccessRoutines
  USE RegionAccessRoutines
  USE SolverAccessRoutines
  USE SolverMappingAccessRoutines
  USE Strings
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters
  !Module types

  !Module variables

  !Interfaces

  PUBLIC SolverMapping_CreateFinish,SolverMapping_CreateStart

  PUBLIC SolverMapping_Destroy
  
  PUBLIC SolverMapping_EquationsSetAdd

  PUBLIC SolverMapping_InterfaceConditionAdd

  PUBLIC SolverMapping_EquatsVarsToSolverMatrixSet

  PUBLIC SolverMapping_NumberOfSolverMatricesSet
  
  CONTAINS

  !
  !=================================================================================================================================
  !
  
  !>Calculates the solver mappings
  SUBROUTINE SolverMapping_Calculate(solverMapping,err,error,*)
    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to calcualte
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,columnListItem(5),columnRank,dofIdx,dofType,equationType, &
      & equationsColumn,equationsIdx,equationsIdx2,equationMatrixIdx,equationsMatrixNumber,equationsRowNumber, &
      & equationsSetIdx,equationsSetUserNumber,equationsVariableListItem(3),globalColumn,globalDOF,globalDOFIdx,globalDOFsOffset, &
      & globalRow,globalRowIdx,interfaceColumn,interfaceColumnNumber,interfaceConditionIdx,interfaceConditionIdx2, &
      & interfaceConditionIndex,interfaceEquationsListItem(2),interfaceMatrixIdx,interfaceMatrixNumber,interfaceRow, &
      & interfaceRowNumber,interfaceConditionUserNumber,interfaceUserNumber,jacobianColumn,jacobianMatrixNumber, &
      & larangeVariableType,localColumn,localDOF,localDOFsOffset,localRow,matricesType,matrixNumber,matrixType,matrixTypeIdx, &
      & matrixVariableIdx,myrank,numberOfColumns,numberOfConstraints,numberOfDOFs,numberOfDependentVariables, &
      & numberOfDynamicMatrices,numberOfEquations,numberOfEquationsColumns,numberOfEquationDOFs,numberOfEquationsSets, &
      & numberOfEquationsSetRows,numberOfEquationsVariables,numberOfGlobalDOFs,numberOfInterfaceConditions, &
      & numberOfInterfaceColumns,numberOfInterfaceMatrices,numberOfInterfaceRows,numberOfInterfaceVariables,numberOfJacobians, &
      & numberOfJacobianMatrices,numberOfGlobalSolverDOFs,numberOfGlobalSolverRows,numberOfLinearMatrices,numberOfLocalSolverDOFs, &
      & numberOfLocalSolverRows,numberOfSolverMatrices,numberOfRankColumns,numberOfRankRows,numberOfResiduals,numberOfRows, &
      & numberOfVariables,parentRegionUserNumber,rank,rankIdx,regionUserNumber,residualIdx,rowIdx,rowListItem(4),rowRank, &
      & solverGlobalDOF,solverMatrixIdx,solverMatrixNumber,solverVariableIdx,totalNumberOfColumns, &
      & totalNumberOfDOFs,totalNumberOfLocalSolverDOFs,totalNumberOfRows, &
      & variableIdx,variableIdx2,variableIndex,variableListItem(3),variablePositionIdx,variableType,numberOfGroupComputationNodes, &
      & numberRowEquationsRows,numberColEquationsCols,rowEquationsRowIdx,colEquationsColIdx, &
      & globalDofCouplingNumber,equationsRow,eqnLocalDof,numberOfEquationsRHSVariables,rhsVariableType
    INTEGER(INTG) :: tempOffset, tempSolverVariableIdx
    INTEGER(INTG), ALLOCATABLE :: equationsSetVariables(:,:),equationsVariables(:,:),interfaceEquationsList(:,:), &
      & interfaceVariables(:,:),rankGlobalRowsList(:,:),rankGlobalColumnsList(:,:),solverLocalDOF(:), &
      & equationsRHSVariables(:,:)
    INTEGER(INTG), ALLOCATABLE :: numberOfVariableGlobalSolverDOFs(:),numberOfVariableLocalSolverDOFs(:), &
      & totalNumberOfVariableLocalSolverDOFs(:),subMatrixInformation(:,:,:),subMatrixList(:,:,:),variableTypes(:)
    REAL(DP) :: additiveConstant,couplingCoefficient,variableCoefficient
    LOGICAL :: found,haveDynamic,haveLinear,haveNonlinear,includeColumn,includeRow,constrainedDOF
    LOGICAL, ALLOCATABLE :: variableProcessed(:),variableRankProcessed(:,:)
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(BoundaryConditionsCoupledDofsType), POINTER :: colEquationsCols,dofCoupling,rowEquationsRows
    TYPE(BoundaryConditionsCoupledDofsType), TARGET :: dummyDofCoupling
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable
    TYPE(DomainMappingType), POINTER :: columnDomainMapping,columnDOFsMapping,rowDomainMapping,rowDOFsMapping
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap
    TYPE(EquationsMatrixToSolverMatricesMapType), POINTER :: equationsMatrixToSolverMatricesMap
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: equationsMatrixToSolverMatrixMap
    TYPE(EquationsMatrixToVarMapType), POINTER :: equationsMatrixToVariableMap
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap
    TYPE(FieldType), POINTER :: dependentField,lagrangeField
    TYPE(FieldVariableType), POINTER :: dependentVariable,fieldVariable,lagrangeVariable,variable,rhsVariable
    TYPE(IntegerIntgPtrType), POINTER :: dofMap(:)
    TYPE(InterfaceType), POINTER :: interface
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMap
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap
    TYPE(InterfaceMatrixToSolverMatricesMapType), POINTER :: interfaceMatrixToSolverMatricesMap
    TYPE(InterfaceMatrixToSolverMatrixMapType), POINTER :: interfaceMatrixToSolverMatrixMap
    TYPE(JacobianMatrixToSolverMatrixMapType), POINTER :: jacobianMatrixToSolverMatrixMap
    TYPE(JacobianMatrixToVarMapType), POINTER :: jacobianMatrixToVariableMap
    TYPE(ListType), POINTER :: equationsSetVariableList
    TYPE(ListPtrType), ALLOCATABLE :: interfaceEquationsLists(:),rankGlobalRowsLists(:,:), &
      & rankGlobalColumnLists(:,:,:,:),variablesList(:)
    TYPE(MatrixRowColCouplingType), POINTER :: equationsRowToSolverRowsMap(:),interfaceColToSolverColsMap(:), &
      & interfaceColToSolverRowsMap(:),interfaceRowToSolverColsMap(:)
    TYPE(RegionType), POINTER :: parentRegion,region
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverColToDynamicEquationsMapType), POINTER :: solverColToDynamicEquationsMap
    TYPE(SolverColToLinearEquationsMapType), POINTER :: solverColToLinearEquationsMap
    TYPE(SolverColToNonlinearEquationsMapType), POINTER :: solverColToNonlinearEquationsMap
    TYPE(SolverDOFToVariableDOFsMapType), POINTER :: solverDOFToVariableDOFsMap
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(SolverMappingDOFCouplingsType) :: columnCouplings,rowCouplings
    TYPE(SolverMappingVariableType), POINTER :: solverMappingVariable
    TYPE(SolverMappingVariablesType), POINTER :: rhsVariablesList,solverVariablesList
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap
    TYPE(SolverMatrixToEquationsSetMapType), POINTER :: solverMatrixToEquationsSetMap
    TYPE(SolverRowToEquationsMapType), POINTER :: solverRowToEquationsMap
    TYPE(VariableDOFToSolverDOFsMapType), POINTER :: dependentVarDOFToSolverDOFsMap,lagrangeVarDOFToSolverDOFsMap, &
      & varDOFToSolverDOFsMap
    TYPE(VARYING_STRING) :: localError
    TYPE(WorkGroupType), POINTER :: workGroup

    ENTERS("SolverMapping_Calculate",err,error,*999)

    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)

    NULLIFY(createValuesCache)
    CALL SolverMapping_CreateValuesCacheGet(solverMapping,createValuesCache,err,error,*999)
    NULLIFY(solverEquations)
    CALL SolverMapping_SolverEquationsGet(solverMapping,solverEquations,err,error,*999)
    NULLIFY(boundaryConditions)
    CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
    NULLIFY(solver)
    CALL SolverEquations_SolverGet(solverEquations,solver,err,error,*999)
    NULLIFY(workGroup)
    CALL Solver_WorkGroupGet(solver,workGroup,err,error,*999)
          
    !
    !--- Equations set <-> interface conditions  ---
    !
    ! 1. Calculate the list interface conditions that influence an equations set and vice versa.
    !
          
    !Allocate equations set to solver map
    ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(solverMapping%numberOfEquationsSets),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver mapping equations set to solver map.",err,error,*999)      
    !Allocate interface condition to solver map
    ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(solverMapping%numberOfInterfaceConditions),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver mapping interface condition to solver map.",err,error,*999)
    !
    ! Allocate and initialise
    !
    ALLOCATE(interfaceEquationsLists(solverMapping%numberOfInterfaceConditions),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set list.",err,error,*999)
    DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
      NULLIFY(interfaceCondition)
      CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
      NULLIFY(interfaceEquations)
      CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
      NULLIFY(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr)
      CALL SolverMappingICToSMSMap_Initialise(solverMapping% &
        & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr,err,error,*999)      
      solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr%interfaceConditionIndex=interfaceConditionIdx
      solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr%solverMapping=>solverMapping
      solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr%interfaceEquations=>interfaceEquations
      NULLIFY(interfaceEquationsLists(interfaceConditionIdx)%ptr)
      CALL List_CreateStart(interfaceEquationsLists(interfaceConditionIdx)%PTR,err,error,*999)
      CALL List_DataTypeSet(interfaceEquationsLists(interfaceConditionIdx)%PTR,LIST_INTG_TYPE,err,error,*999)
      CALL List_DataDimensionSet(interfaceEquationsLists(interfaceConditionIdx)%PTR,2,err,error,*999)
      CALL List_CreateFinish(interfaceEquationsLists(interfaceConditionIdx)%PTR,err,error,*999)
    ENDDO !interfaceConditionIdx
    !
    ! Loop over equations sets
    !
    DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      NULLIFY(equations)
      CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
      NULLIFY(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr)
      CALL SolverMappingESToSMSMap_Initialise(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr,err,error,*999)
      solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%equationsSetIndex=equationsSetIdx
      solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%solverMapping=>solverMapping
      solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%equations=>equations
      !Set up list of interface conditions affecting this equations set
      CALL List_DetachAndDestroy(createValuesCache%interfaceIndices(equationsSetIdx)%ptr,numberOfInterfaceConditions, &
        & interfaceEquationsList,err,error,*999)
      ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
        & interfaceConditionToEquationsSetMaps(numberOfInterfaceConditions),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations to solver maps interface.",err,error,*999)
      solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%numberOfInterfaceConditions=numberOfInterfaceConditions
      DO interfaceConditionIdx=1,numberOfInterfaceConditions
        CALL SolverMappingSMICToESMap_Initialise(solverMapping%equationsSetToSolverMatricesMaps( &
          & equationsSetIdx)%ptr%interfaceConditionToEquationsSetMaps(interfaceConditionIdx),err,error,*999)
        interfaceConditionIdx2=interfaceEquationsList(1,interfaceConditionIdx)
        interfaceMatrixIdx=interfaceEquationsList(2,interfaceConditionIdx)
        NULLIFY(interfaceCondition)
        CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx2,interfaceCondition,err,error,*999)
        solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%interfaceConditionToEquationsSetMaps( &
          & interfaceConditionIdx)%interfaceConditionIndex=interfaceConditionIdx2
        solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%interfaceConditionToEquationsSetMaps( &
          & interfaceConditionIdx)%interfaceCondition=>interfaceCondition
        solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%interfaceConditionToEquationsSetMaps( &
          & interfaceConditionIdx)%interfaceMatrixNumber=interfaceMatrixIdx
        interfaceEquationsListItem(1)=equationsSetIdx
        interfaceEquationsListItem(2)=interfaceMatrixIdx
        CALL List_ItemAdd(interfaceEquationsLists(interfaceConditionIdx2)%ptr,interfaceEquationsListItem,err,error,*999)
      ENDDO !interfaceConditionIdx
      IF(ALLOCATED(interfaceEquationsList)) DEALLOCATE(interfaceEquationsList)
    ENDDO !equationsSetIdx
    DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
      CALL List_DetachAndDestroy(interfaceEquationsLists(interfaceConditionIdx)%ptr,numberOfEquationsSets, &
        & interfaceEquationsList,err,error,*999)
      ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
        & equationsSetToInterfaceConditionMaps(numberOfEquationsSets),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate interface to solver maps equations.",err,error,*999)
      solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr%numberOfEquationsSets=numberOfEquationsSets
      DO equationsIdx=1,numberOfEquationsSets
        CALL SolverMappingSMESToICMap_Initialise(solverMapping%interfaceConditionToSolverMatricesMaps( &
          & interfaceConditionIdx)%ptr%equationsSetToInterfaceConditionMaps(equationsIdx),err,error,*999)
        equationsSetIdx=interfaceEquationsList(1,equationsIdx)
        interfaceMatrixIdx=interfaceEquationsList(2,equationsIdx)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
          & equationsSetToInterfaceConditionMaps(equationsIdx)%equationsSetIndex=equationsSetIdx
        solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
          & equationsSetToInterfaceConditionMaps(equationsIdx)%equationsSet=>equationsSet
        solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
          & equationsSetToInterfaceConditionMaps(equationsIdx)%interfaceMatrixIndex=interfaceMatrixIdx
      ENDDO !equationsIdx
      IF(ALLOCATED(interfaceEquationsList)) DEALLOCATE(interfaceEquationsList)
    ENDDO !interfaceConditionIdx
    !
    !--- Row mappings ---
    !
    ! 2. Determine the number of rows in the solver matrix. Do this the by setting up a list of rows for each rank.
    !    We can then later arrange the rows in rank order by looping over the ranks in the list and then the rows
    !    for each rank.
    !
    !Calculate the row mappings.
    CALL WorkGroup_NumberOfGroupNodesGet(workGroup,numberOfGroupComputationNodes,err,error,*999)
    CALL WorkGroup_GroupNodeNumberGet(workGroup,myrank,err,error,*999)
    numberOfGlobalSolverRows=0
    numberOfLocalSolverRows=0
    !Add in the rows from any equations sets that have been added to the solver equations
    !Presort the row numbers by rank.
    !
    !Allocate and initialise the rank lists.
    ALLOCATE(rankGlobalRowsLists(solverMapping%numberOfEquationsSets+solverMapping%numberOfInterfaceConditions, &
      & 0:numberOfGroupComputationNodes-1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate rank global rows lists.",err,error,*999)
    CALL SolverMappingDOFCouplings_Initialise(rowCouplings,err,error,*999)
    DO rank=0,numberOfGroupComputationNodes-1
      equationsIdx=0
      DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
        equationsIdx=equationsIdx+1
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(equations)
        CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        NULLIFY(vectorMapping)
        CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
        NULLIFY(lhsMapping)
        CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
        NULLIFY(rankGlobalRowsLists(equationsIdx,rank)%ptr)
        CALL List_CreateStart(rankGlobalRowsLists(equationsIdx,rank)%ptr,err,error,*999)
        CALL List_DataTypeSet(rankGlobalRowsLists(equationsIdx,rank)%ptr,LIST_INTG_TYPE,err,error,*999)
        CALL List_InitialSizeSet(rankGlobalRowsLists(equationsIdx,rank)%ptr, &
          & INT(lhsMapping%numberOfGlobalRows/numberOfGroupComputationNodes,INTG),err,error,*999)
        CALL List_DataDimensionSet(rankGlobalRowsLists(equationsIdx,rank)%ptr,4,err,error,*999)
        CALL List_KeyDimensionSet(rankGlobalRowsLists(equationsIdx,rank)%ptr,1,err,error,*999)
        CALL List_CreateFinish(rankGlobalRowsLists(equationsIdx,rank)%ptr,err,error,*999)  
      ENDDO !equationsSetIdx
      DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
        equationsIdx=equationsIdx+1
        NULLIFY(interfaceCondition)
        CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
        SELECT CASE(interfaceCondition%method)
        CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
          NULLIFY(interfaceEquations)
          CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
          NULLIFY(interfaceMapping)
          CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
          NULLIFY(rankGlobalRowsLists(equationsIdx,rank)%ptr)
          CALL List_CreateStart(rankGlobalRowsLists(equationsIdx,rank)%ptr,err,error,*999)
          CALL List_DataTypeSet(rankGlobalRowsLists(equationsIdx,rank)%ptr,LIST_INTG_TYPE,err,error,*999)
          CALL List_InitialSizeSet(rankGlobalRowsLists(equationsIdx,rank)%ptr, &
            & INT(interfaceMapping%numberOfGlobalColumns/numberOfGroupComputationNodes,INTG),err,error,*999)
          CALL List_DataDimensionSet(rankGlobalRowsLists(equationsIdx,rank)%ptr,4,err,error,*999)
          CALL List_KeyDimensionSet(rankGlobalRowsLists(equationsIdx,rank)%ptr,1,err,error,*999)
          CALL List_CreateFinish(rankGlobalRowsLists(equationsIdx,rank)%ptr,err,error,*999)
        CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The interface condition method of "//TRIM(NumberToVString(interfaceCondition%method,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !interfaceConditionIdx
    ENDDO !rank
    !Calculate the number of local and global rows. Do this by looking at the boundary conditions for field variables
    !involved in the row. If all the variables are set as a fixed boundary condition then do not include the row. If
    !any variable is not fixed then include the row.
    equationsIdx=0
    DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
      equationsIdx=equationsIdx+1
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      NULLIFY(equations)
      CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
      NULLIFY(vectorEquations)
      CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
      NULLIFY(vectorMapping)
      CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
      NULLIFY(lhsMapping)
      CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
      NULLIFY(rowDOFsMapping)
      CALL EquationsMappingLHS_RowDOFsMappingGet(lhsMapping,rowDOFsMapping,err,error,*999)
      NULLIFY(dynamicMapping)
      CALL EquationsMappingVector_DynamicMappingExists(vectorMapping,dynamicMapping,err,error,*999)
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
      NULLIFY(nonlinearMapping)
      CALL EquationsMappingVector_NonlinearMappingExists(vectorMapping,nonlinearMapping,err,error,*999)
      NULLIFY(rhsMapping)
      CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
      NULLIFY(sourceMapping)
      CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*999)
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      !Loop over the global rows for this equations set
      DO globalRow=1,lhsMapping%numberOfGlobalRows
        !Find the rank that owns this global row
        rowRank=-1
        DO rankIdx=1,rowDOFsMapping%globalToLocalMap(globalRow)%numberOfDomains
          IF(rowDOFsMapping%globalToLocalMap(globalRow)%localType(rankIdx)/=DOMAIN_LOCAL_GHOST) THEN
            rowRank=rowDOFsMapping%globalToLocalMap(globalRow)%domainNumber(rankIdx)
            localRow=rowDOFsMapping%globalToLocalMap(globalRow)%localNumber(rankIdx)
            EXIT
          ENDIF
        ENDDO !rankIdx
        IF(rowRank<0) THEN
          localError="Global row "//TRIM(NumberToVString(globalRow,"*",err,error))//" is not owned by a domain."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        includeRow=.TRUE.
        constrainedDOF=.FALSE.
        globalDofCouplingNumber=0
        IF(ASSOCIATED(dynamicMapping)) THEN
          NULLIFY(dependentVariable)
          CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,dependentVariable,err,error,*999)
          NULLIFY(boundaryConditionsVariable)
          CALL BoundaryConditions_VariableGet(boundaryConditions,dependentVariable,boundaryConditionsVariable,err,error,*999)
          !This is wrong as we only have the mappings for the local rank not the global ranks.
          !For now assume 1-1 mapping between rows and dofs.
          globalDOF=globalRow
          includeRow=includeRow.AND.(boundaryConditionsVariable%DOFTypes(globalDOF)==BOUNDARY_CONDITION_DOF_FREE)
          constrainedDOF=constrainedDOF.OR.(boundaryConditionsVariable%DOFTypes(globalDOF)==BOUNDARY_CONDITION_DOF_CONSTRAINED)
          NULLIFY(dofConstraints)
          CALL BoundaryConditionsVariable_DOFConstraintsExists(boundaryConditionsVariable,dofConstraints,err,error,*999)
          IF(ASSOCIATED(dofConstraints)) THEN
            CALL BoundaryConditionsVariableDOFConstraints_NumberOfConstraintsGet(dofConstraints,numberOfConstraints,err,error,*999)
            IF(numberOfConstraints>0) THEN
              NULLIFY(dofCoupling)
              CALL BoundaryConditionsVariableDOFConstraints_DOFCouplingExists(dofConstraints,globalDOF,dofCoupling,err,error,*999)
              IF(ASSOCIATED(dofCoupling)) THEN
                !This equations row is the owner of a solver row that is mapped to
                !multiple other equations rows, add it to the list of global row
                !couplings and remember the index into the global list for this solver row
                CALL SolverDofCouplings_AddCoupling(rowCouplings,dofCoupling,globalDofCouplingNumber,err,error,*999)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
        IF(ASSOCIATED(nonlinearMapping)) THEN
          !Look at the boundary conditions for nonlinear variables for this row
          DO residualIdx=1,nonlinearMapping%numberOfResiduals
            NULLIFY(residualMapping)
            CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
            CALL EquationsMappingResidual_NumberOfResidualVariablesGet(residualMapping,numberOfVariables,err,error,*999)
            DO variableIdx=1,numberOfVariables
              NULLIFY(dependentVariable)
              CALL EquationsMappingResidual_VariableGet(residualMapping,variableIdx,dependentVariable,err,error,*999)
              NULLIFY(boundaryConditionsVariable)
              CALL BoundaryConditions_VariableGet(boundaryConditions,dependentVariable,boundaryConditionsVariable,err,error,*999)
              globalDOF=globalRow
              includeRow=includeRow.AND.(boundaryConditionsVariable%DOFTypes(globalDOF)==BOUNDARY_CONDITION_DOF_FREE)
              constrainedDOF=constrainedDOF.OR.(boundaryConditionsVariable%DOFTypes(globalDOF)==BOUNDARY_CONDITION_DOF_CONSTRAINED)
              NULLIFY(dofConstraints)
              CALL BoundaryConditionsVariable_DOFConstraintsExists(boundaryConditionsVariable,dofConstraints,err,error,*999)
              IF(ASSOCIATED(dofConstraints)) THEN
                CALL BoundaryConditionsVariableDOFConstraints_NumberOfConstraintsGet(dofConstraints,numberOfConstraints, &
                  & err,error,*999)
                IF(numberOfConstraints>0) THEN
                  NULLIFY(dofCoupling)
                  CALL BoundaryConditionsVariableDOFConstraints_DOFCouplingExists(dofConstraints,globalDOF,dofCoupling, &
                    & err,error,*999)
                  IF(ASSOCIATED(dofCoupling)) THEN
                    CALL SolverDofCouplings_AddCoupling(rowCouplings,dofCoupling,globalDofCouplingNumber,err,error,*999)
                  ENDIF
                ENDIF
              ENDIF
            ENDDO !variableIdx
          ENDDO !residualIdx
        ENDIF
        IF(ASSOCIATED(linearMapping)) THEN
          !Loop over the variables in the equations set. Don't include the row in the solver matrices if
          !all the variable dofs associated with this equations row are fixed.
          DO variableIdx=1,linearMapping%numberOfLinearVariables
            NULLIFY(dependentVariable)
            CALL EquationsMappingLinear_LinearVariableGet(linearMapping,variableIdx,dependentVariable,err,error,*999)
            NULLIFY(boundaryConditionsVariable)
            CALL BoundaryConditions_VariableGet(boundaryConditions,dependentVariable,boundaryConditionsVariable,err,error,*999)
            !\todo This is wrong as we only have the mappings for the local rank not the global ranks. See below
            !\todo For now assume 1-1 mapping between rows and dofs.
            globalDOF=globalRow
            includeRow=includeRow.AND.(boundaryConditionsVariable%DOFTypes(globalDOF)==BOUNDARY_CONDITION_DOF_FREE)
            constrainedDOF=constrainedDOF.OR.(boundaryConditionsVariable%DOFTypes(globalDOF)==BOUNDARY_CONDITION_DOF_CONSTRAINED)
            NULLIFY(dofConstraints)
            CALL BoundaryConditionsVariable_DOFConstraintsExists(boundaryConditionsVariable,dofConstraints,err,error,*999)
            IF(ASSOCIATED(dofConstraints)) THEN
              CALL BoundaryConditionsVariableDOFConstraints_NumberOfConstraintsGet(dofConstraints,numberOfConstraints, &
                & err,error,*999)
              IF(numberOfConstraints>0) THEN
                NULLIFY(dofCoupling)
                CALL BoundaryConditionsVariableDOFConstraints_DOFCouplingExists(dofConstraints,globalDOF,dofCoupling, &
                  & err,error,*999)
                IF(ASSOCIATED(dofCoupling)) THEN
                  CALL SolverDofCouplings_AddCoupling(rowCouplings,dofCoupling,globalDofCouplingNumber,err,error,*999)
                ENDIF
              ENDIF
            ENDIF
          ENDDO !variableIdx
        ENDIF
        rowListItem(1)=globalRow
        rowListItem(2)=localRow
        rowListItem(4)=globalDofCouplingNumber
        IF(includeRow) THEN
          rowListItem(3)=1
          numberOfGlobalSolverRows=numberOfGlobalSolverRows+1
          !Don't need to worry about ghosted rows.
          IF(rowRank==myrank) numberOfLocalSolverRows=numberOfLocalSolverRows+1 !1-1 mapping
        ELSE IF(constrainedDOF) THEN
          rowListItem(3)=2
        ELSE
          rowListItem(3)=0
        ENDIF !include row
        CALL List_ItemAdd(rankGlobalRowsLists(equationsIdx,rowRank)%ptr,rowListItem,err,error,*999)
      ENDDO !globalRow
    ENDDO !equations set idx
    !Now add in rows from any interface matrices
    DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
      equationsIdx=equationsIdx+1
      NULLIFY(interfaceCondition)
      CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
      SELECT CASE(interfaceCondition%method)
      CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
        NULLIFY(interfaceEquations)
        CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
        NULLIFY(interfaceMapping)
        CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfacemapping,err,error,*999)
        NULLIFY(columnDOFsMapping)
        CALL InterfaceMapping_ColumnDOFsMappingGet(interfaceMapping,columnDOFsMapping,err,error,*999)
        NULLIFY(lagrangeField)
        CALL InterfaceCondition_LagrangeFieldGet(interfaceCondition,lagrangeField,err,error,*999)
        NULLIFY(boundaryConditions)
        CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
        !\todo Lagrange variable type set to the first variable type for now
        variableType=1
        NULLIFY(lagrangeVariable)
        CALL Field_VariableGet(lagrangeField,variableType,lagrangeVariable,err,error,*999)
        NULLIFY(boundaryConditionsVariable)
        CALL BoundaryConditions_VariableGet(boundaryConditions,lagrangeVariable,boundaryConditionsVariable,err,error,*999)
        !Loop over the global columns for this interface equation
        DO globalColumn=1,interfaceMapping%numberOfGlobalColumns
          !Find the rank that owns this global column
          columnRank=-1
          DO rankIdx=1,columnDOFsMapping%globalToLocalMap(globalColumn)%numberOfDomains
            IF(columnDOFsMapping%globalToLocalMap(globalColumn)%localType(rankIdx)/=DOMAIN_LOCAL_GHOST) THEN
              columnRank=columnDOFsMapping%globalToLocalMap(globalColumn)%domainNumber(rankIdx)
              localColumn=columnDOFsMapping%globalToLocalMap(globalColumn)%localNumber(rankIdx)
              EXIT
            ENDIF
          ENDDO !rankIdx
          IF(columnRank<0) THEN
            localError="Global column "//TRIM(NumberToVString(globalColumn,"*",err,error))//" is not owned by a domain."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          includeColumn=.TRUE.
          !\todo This is wrong as we only have the mappings for the local rank not the global ranks. See above
          !\todo For now assume 1-1 mapping between rows and dofs.
          globalDOF=globalColumn
          includeColumn=includeColumn.AND.(boundaryConditionsVariable%DOFTypes(globalDOF)==BOUNDARY_CONDITION_DOF_FREE)
          rowListItem(1)=globalColumn
          rowListItem(2)=localColumn
          IF(includeColumn) THEN
            rowListItem(3)=1
            numberOfGlobalSolverRows=numberOfGlobalSolverRows+1
            !Don't need to worry about ghosted rows.
            IF(columnRank==myrank) numberOfLocalSolverRows=numberOfLocalSolverRows+1 !1-1 mapping
          ELSE
            rowListItem(3)=0
          ENDIF !include column
          CALL List_ItemAdd(rankGlobalRowsLists(equationsIdx,columnRank)%ptr,rowListItem,err,error,*999)
        ENDDO !globalColumn
      CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The interface condition method of "// &
          & TRIM(NumberToVString(interfaceCondition%method,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !interfaceConditionIdx

    !Sanity check.
    IF(numberOfLocalSolverRows==0) CALL FlagError("Invalid problem setup. The number of local solver rows is zero.",err,error,*999)
    IF(numberOfGlobalSolverRows==0) &
      & CALL FlagError("Invalid problem setup. The number of global solver rows is zero.",err,error,*999)

    !
    ! 3. We now know how many local and global rows are in the solver matrix. Loop over the rows in rank order and calculate
    !    the row mappings.
    !
    ! 3a Allocate and initialise the data structures
    !
    
    !Allocate memory for the rows mapping
    !Allocate the solver rows to equations set maps
    ALLOCATE(solverMapping%solverRowToEquationsMaps(numberOfLocalSolverRows),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver mapping solver row to equation rows map.",err,error,*999)
    !Set the number of rows
    solverMapping%numberOfRows=numberOfLocalSolverRows
    solverMapping%numberOfGlobalRows=numberOfGlobalSolverRows
    !Allocate the solver rows domain mapping
    CALL DomainMapping_Initialise(solverMapping%rowDOFsMapping,err,error,*999)
    CALL DomainMapping_WorkGroupSet(solverMapping%rowDOFsMapping,workGroup,err,error,*999)
    NULLIFY(rowDomainMapping)
    CALL SolverMapping_RowDOFsMappingGet(solverMapping,rowDomainMapping,err,error,*999)
    ALLOCATE(rowDomainMapping%globalToLocalMap(numberOfGlobalSolverRows),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate row dofs mapping global to local map.",err,error,*999)
    rowDomainMapping%numberOfGlobal=numberOfGlobalSolverRows

    !Initialise the equations sets to solver matrices maps
    DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      NULLIFY(equations)
      CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
      NULLIFY(vectorEquations)
      CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
      NULLIFY(vectorMapping)
      CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
      NULLIFY(lhsMapping)
      CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
            
      !Allocate the equations matrices to solver matrix maps
      ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
        & equationsMatricesToSolverMatrixMaps(solverMapping%numberOfSolverMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations set to solver matrices equations matrices to solver matrix maps.", &
        & err,error,*999)
!!TODO: This should just be the actual number of solver matrices that the equations set matrices are mapped to rather than the total number.
      solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%numberOfSolverMatrices= &
        & solverMapping%numberOfSolverMatrices
      DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
        NULLIFY(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
          & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr)
        CALL SolverMappingEMSToSMMap_Initialise(solverMapping% &
          & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
          & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr,err,error,*999)
      ENDDO !solverMatrixIdx
      
      !Allocate the equations row to solver rows maps
      ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
        & equationsRowToSolverRowsMap(lhsMapping%totalNumberOfRows),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations set set to solver matrices map equations row to solver rows maps.", &
        & err,error,*999)
      DO equationsRowNumber=1,lhsMapping%totalNumberOfRows
        !Initialise
        CALL MatrixRowColCoupling_Initialise(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
          & equationsRowToSolverRowsMap(equationsRowNumber),err,error,*999)
      ENDDO !equationsRowIdx
             
    ENDDO !equationsSetIdx
          
    !Initialise the interface condition to solver matrices maps
    DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
      NULLIFY(interfaceCondition)
      CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
      NULLIFY(interfaceEquations)
      CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
      NULLIFY(interfaceMapping)
      CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
      
      !Allocate the interface matrices to solver matrix maps
      ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
        & interfaceMatricesToSolverMatrixMaps(solverMapping%numberOfSolverMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate interface condition to solver matrices map interface matrices"// &
        & " to solver matrix maps.",err,error,*999)
      solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr%numberOfSolverMatrices= &
        & solverMapping%numberOfSolverMatrices
      DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
        NULLIFY(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
          & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr)
        CALL SolverMappingIMSToSMMap_Initialise(solverMapping% &
          & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
          & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr,err,error,*999)
      ENDDO !solverMatrixIdx
            
      !Allocate the interface matrix to solver matrices maps
      ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
        & interfaceMatrixToSolverMatricesMaps(interfaceMapping%numberOfInterfaceMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate interface condition to solver matrices map interface matrix"// &
        & " to solver matrices maps.",err,error,*999)
      solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr%numberOfInterfaceMatrices= &
        & interfaceMapping%numberOfInterfaceMatrices
      DO interfaceMatrixIdx=1,interfaceMapping%numberOfInterfaceMatrices
        NULLIFY(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
          & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr)
        CALL SolverMappingIMToSMSMap_Initialise(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
          & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr,err,error,*999)
                      
        !Allocate the interface row to solver row maps
        ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
          & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr%&
          & interfaceRowToSolverRowsMap(interfaceMapping%interfaceMatrixToVarMaps(interfaceMatrixIdx)%ptr%totalNumberOfRows), &
          & STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate interface condition to solver matrices map"// &
          & " interface matrix to solver matrices maps interface row to solver rows map.",err,error,*999)
        DO interfaceRowNumber=1,interfaceMapping%interfaceMatrixToVarMaps(interfaceMatrixIdx)%ptr%totalNumberOfRows
          CALL MatrixRowColCoupling_Initialise(solverMapping% &
            & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr% &
            & interfaceRowToSolverRowsMap(interfaceRowNumber),err,error,*999)
        ENDDO !interfaceRowNumber
        
      ENDDO !interfaceMatrixIdx

      !Allocate the interface column to solver row maps
      ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
        & interfaceColToSolverRowsMap(interfaceMapping%totalNumberOfColumns),STAT=err)
      IF(err/=0)  &
        & CALL FlagError("Could not allocate interface condition to solver matrices map interface column to solver rows map.", &
        & err,error,*999)
      DO interfaceColumnNumber=1,interfaceMapping%totalNumberOfColumns
        !Initialise
        CALL MatrixRowColCoupling_Initialise(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
          & interfaceColToSolverRowsMap(interfaceColumnNumber),err,error,*999)
      ENDDO !interfaceColumnNumber
      
    ENDDO !interfaceConditionIdx
    
    !
    ! 3b Now calculate the mappings for each solver row <-> equations row & interface row/column in rank order.
    !
    
    numberOfGlobalSolverRows=0
    !Make a "dof coupling" for rows that aren't coupled
    ALLOCATE(dummyDofCoupling%globalDofs(1),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate dummy DOF coupling DOFs.",err,error,*999)
    ALLOCATE(dummyDofCoupling%localDofs(1),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate dummy DOF coupling local Dofs.",err,error,*999)
    ALLOCATE(dummyDofCoupling%coefficients(1),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate dummy DOF coupling values.",err,error,*999)
    dummyDofCoupling%numberOfDofs=1
    !Loop over the ranks to ensure that the lowest ranks have the lowest numbered solver variables
    DO rank=0,numberOfGroupComputationNodes-1
      numberOfLocalSolverRows=0

      !Calculate the solver row <-> equations row & interface row mappings.
      equationsIdx=0
      DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
        equationsIdx=equationsIdx+1
        
        !Get rows list
        CALL List_Sort(rankGlobalRowsLists(equationsIdx,rank)%ptr,err,error,*999)
        CALL List_DetachAndDestroy(rankGlobalRowsLists(equationsIdx,rank)%ptr,numberOfRankRows,rankGlobalRowsList,err,error,*999)

        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(equations)
        CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        NULLIFY(vectorMapping)
        CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
        NULLIFY(dynamicMapping)
        CALL EquationsMappingVector_DynamicMappingExists(vectorMapping,dynamicMapping,err,error,*999)
        NULLIFY(linearMapping)
        CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
        NULLIFY(nonlinearMapping)
        CALL EquationsMappingVector_NonlinearMappingExists(vectorMapping,nonlinearMapping,err,error,*999)

        !Loop over the global rows for this rank.
        DO globalRowIdx=1,numberOfRankRows
          globalRow=rankGlobalRowsList(1,globalRowIdx)
          localRow=rankGlobalRowsList(2,globalRowIdx)
          includeRow=(rankGlobalRowsList(3,globalRowIdx)==1)
          constrainedDOF=(rankGlobalRowsList(3,globalRowIdx)==2)
          globalDofCouplingNumber=rankGlobalRowsList(4,globalRowIdx)
          IF(globalDofCouplingNumber>0) THEN
            rowEquationsRows=>rowCouplings%dofCouplings(globalDofCouplingNumber)%ptr
            IF(.NOT.ASSOCIATED(rowEquationsRows)) THEN
              localError="DOF coupling is not associated for global DOF coupling number "// &
                & TRIM(NumberToVstring(globalDofCouplingNumber,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            numberRowEquationsRows=rowEquationsRows%numberOfDofs
          ELSE
            numberRowEquationsRows=1
            dummyDofCoupling%globalDofs(1)=globalRow
            dummyDofCoupling%localDofs(1)=localRow
            dummyDofCoupling%coefficients(1)=1.0_DP
            rowEquationsRows=>dummyDofCoupling
          ENDIF

          IF(includeRow) THEN
            numberOfGlobalSolverRows=numberOfGlobalSolverRows+1
            numberOfLocalSolverRows=numberOfLocalSolverRows+1
            !Set up the row domain mappings.
            !There are no ghosted rows for the solver matrices so there is only one domain for the global to local map.
            !Initialise
            CALL DomainGlobalMapping_Initialise(rowDomainMapping%globalToLocalMap(numberOfGlobalSolverRows),err,error,*999)
            !Allocate the global to local map arrays
            ALLOCATE(rowDomainMapping%globalToLocalMap(numberOfGlobalSolverRows)%localNumber(1),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate row global to local map local number.",err,error,*999)
            ALLOCATE(rowDomainMapping%globalToLocalMap(numberOfGlobalSolverRows)%domainNumber(1),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate row global to local map domain number.",err,error,*999)
            ALLOCATE(rowDomainMapping%globalToLocalMap(numberOfGlobalSolverRows)%localType(1),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate row global to local map local type.",err,error,*999)
            !Set the global to local mappings
            rowDomainMapping%globalToLocalMap(numberOfGlobalSolverRows)%numberOfDomains=1
            rowDomainMapping%globalToLocalMap(numberOfGlobalSolverRows)%localNumber(1)=numberOfLocalSolverRows
            rowDomainMapping%globalToLocalMap(numberOfGlobalSolverRows)%domainNumber(1)=rank
            rowDomainMapping%globalToLocalMap(numberOfGlobalSolverRows)%localType(1)=DOMAIN_LOCAL_INTERNAL
            IF(rank==myrank) THEN
              !If this is my rank then set up the solver->equations and equations->solver row mappings
                    
              !Set up the solver row -> equations row mappings. Will need to look
              !At the interface conditions for this equations set later.
              !Initialise
              NULLIFY(solverMapping%solverRowToEquationsMaps(numberOfLocalSolverRows)%ptr)
              CALL SolverMappingSRToEQSMap_Initialise(solverMapping%solverRowToEquationsMaps(numberOfLocalSolverRows)%ptr, &
                & err,error,*999)
              
              ALLOCATE(solverMapping%solverRowToEquationsMaps(numberOfLocalSolverRows)%ptr% &
                & equationsIndex(numberRowEquationsRows),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate solver row to equations rows equations index.",err,error,*999)
              ALLOCATE(solverMapping%solverRowToEquationsMaps(numberOfLocalSolverRows)%ptr% &
                & rowColNumber(numberRowEquationsRows),STAT=err)
              IF(err/=0) &
                & CALL FlagError("Could not allocate solver row to equations rows row/col number.",err,error,*999)
              ALLOCATE(solverMapping%solverRowToEquationsMaps(numberOfLocalSolverRows)%ptr% &
                & couplingCoefficients(numberRowEquationsRows),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate solver row to equations rows coupling coefficients.", &
                & err,error,*999)
              !Set the mappings for the first equations DOF, the rest will be set up later using the DOF constraints
              solverMapping%solverRowToEquationsMaps(numberOfLocalSolverRows)%ptr% &
                & numberOfEquationsSetRows=numberRowEquationsRows
              DO rowEquationsRowIdx=1,numberRowEquationsRows
                solverMapping%solverRowToEquationsMaps(numberOfLocalSolverRows)%ptr% &
                  & equationsIndex(rowEquationsRowIdx)=equationsSetIdx
                solverMapping%solverRowToEquationsMaps(numberOfLocalSolverRows)%ptr% &
                  & rowColNumber(rowEquationsRowIdx)=rowEquationsRows%localDofs(rowEquationsRowIdx)
                solverMapping%solverRowToEquationsMaps(numberOfLocalSolverRows)%ptr% &
                  & couplingCoefficients(rowEquationsRowIdx)=rowEquationsRows%coefficients(rowEquationsRowIdx)
              ENDDO !rowEquationsRowIdx
              !Set up the equations row -> solver row mappings
              DO rowEquationsRowIdx=1,numberRowEquationsRows
                equationsRow=rowEquationsRows%localDofs(rowEquationsRowIdx)
                !Allocate the equations row to solver row mappings arrays
                IF(ALLOCATED(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & equationsRowToSolverRowsMap(equationsRow)%rowCols)) THEN
                  CALL FlagError("Equations row to solver row map is already allocated, "// &
                    & "mapping an equations row to multiple solver rows is not yet implemented.",err,error,*999)
                ENDIF
                ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & equationsRowToSolverRowsMap(equationsRow)%rowCols(1),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate equations row to solver rows maps solver rows.",err,error,*999)
                ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & equationsRowToSolverRowsMap(equationsRow)%couplingCoefficients(1),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate equations row to solver rows maps coupling coefficients.", &
                  & err,error,*999)
                
                !Set the mappings
                solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & equationsRowToSolverRowsMap(equationsRow)%numberOfRowCols=1
                solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%&
                  & equationsRowToSolverRowsMap(equationsRow)%rowCols(1)=numberOfLocalSolverRows
                solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & equationsRowToSolverRowsMap(equationsRow)%couplingCoefficients(1)= &
                  & rowEquationsRows%coefficients(rowEquationsRowIdx)
                !Now set up any interface condition rows to solver rows that affect this equations set.
                DO interfaceConditionIdx=1,solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & numberOfInterfaceConditions
                  interfaceConditionIdx2=solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                    & interfaceConditionToEquationsSetMaps(interfaceConditionIdx)%interfaceConditionIndex
                  interfaceMatrixIdx=solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                    & interfaceConditionToEquationsSetMaps(interfaceConditionIdx)%interfaceMatrixNumber
                  ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx2)%ptr% &
                    & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr%&
                    & interfaceRowToSolverRowsMap(equationsRow)%rowCols(1),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate interface row to solver rows map solver row columns.", &
                    & err,error,*999)
                  ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx2)%ptr% &
                    & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr% &
                    & interfaceRowToSolverRowsMap(equationsRow)%couplingCoefficients(1),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate interface row to solver rows map coupling coefficients.", &
                    & err,error,*999)
                  !Set the mappings
                  solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx2)%ptr% &
                    & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr% &
                    & interfaceRowToSolverRowsMap(equationsRow)%numberOfRowCols=1
                  solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx2)%ptr% &
                    & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr% &
                    & interfaceRowToSolverRowsMap(equationsRow)%rowCols(1)=numberOfLocalSolverRows
                  solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx2)%ptr% &
                    & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr% &
                    & interfaceRowToSolverRowsMap(equationsRow)%couplingCoefficients(1)= &
                    & rowEquationsRows%coefficients(rowEquationsRowIdx)
                ENDDO !interfaceConditionIdx
              ENDDO !rowEquationsRowIdx
                    
            ENDIF !rank==my rank
          ELSE IF(constrainedDOF) THEN
            !Do nothing, the row mapping for this equations row will be set up with
            !the first equations row mapped to the solver row
          ELSE
            !Note that in this case these mappings are set to zero since these equation rows don't appear in the solver matrices
            IF(rank==myrank) THEN
              !Set the mappings
              CALL MatrixRowColCoupling_Initialise(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsRowToSolverRowsMap(localRow),err,error,*999)
              !Now set up any interface condition rows to solver rows that affect this equations set.
              DO interfaceConditionIdx=1,solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & numberOfInterfaceConditions
                interfaceConditionIdx2=solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & interfaceConditionToEquationsSetMaps(interfaceConditionIdx)%interfaceConditionIndex
                interfaceMatrixIdx=solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & interfaceConditionToEquationsSetMaps(interfaceConditionIdx)%interfaceMatrixNumber
                !Set the mappings
                CALL MatrixRowColCoupling_Initialise(solverMapping% &
                  & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx2)%ptr% &
                  & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr% &
                  & interfaceRowToSolverRowsMap(localRow),err,error,*999)
              ENDDO !interfaceConditionIdx
            ENDIF !rank==my rank
          ENDIF !include row
        ENDDO !globalRowIdx
        IF(ALLOCATED(rankGlobalRowsList)) DEALLOCATE(rankGlobalRowsList)
      ENDDO !equationsSetIdx

      !Calculate the solver row <-> interface column/row mappings
      DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
        equationsIdx=equationsIdx+1

        !Get rows list
        CALL List_Sort(rankGlobalRowsLists(equationsIdx,rank)%ptr,err,error,*999)
        CALL List_DetachAndDestroy(rankGlobalRowsLists(equationsIdx,rank)%ptr,numberOfRankRows,rankGlobalRowsList,err,error,*999)

        NULLIFY(interfaceCondition)
        CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
        NULLIFY(interfaceEquations)
        CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
        NULLIFY(interfaceMapping)
        CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
        
        !Loop over the global rows for this rank.
        DO globalRowIdx=1,numberOfRankRows
          globalColumn=rankGlobalRowsList(1,globalRowIdx)
          localColumn=rankGlobalRowsList(2,globalRowIdx)
          includeColumn=(rankGlobalRowsList(3,globalRowIdx)==1)
          constrainedDOF=(rankGlobalRowsList(3,globalRowIdx)==2)
          IF(includeColumn) THEN
            numberOfGlobalSolverRows=numberOfGlobalSolverRows+1
            numberOfLocalSolverRows=numberOfLocalSolverRows+1
            !Set up the row domain mappings.
            !There are no ghosted rows for the solver matrices so there is only one domain for the global to local map.
            !Initialise
            CALL DomainGlobalMapping_Initialise(rowDomainMapping%globalToLocalMap(numberOfGlobalSolverRows),err,error,*999)
            !Allocate the global to local map arrays
            ALLOCATE(rowDomainMapping%globalToLocalMap(numberOfGlobalSolverRows)%localNumber(1),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate row global to local map local number.",err,error,*999)
            ALLOCATE(rowDomainMapping%globalToLocalMap(numberOfGlobalSolverRows)%domainNumber(1),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate row global to local map domain number.",err,error,*999)
            ALLOCATE(rowDomainMapping%globalToLocalMap(numberOfGlobalSolverRows)%localType(1),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate row global to local map local type.",err,error,*999)
            !Set the global to local mappings
            rowDomainMapping%globalToLocalMap(numberOfGlobalSolverRows)%numberOfDomains=1
            rowDomainMapping%globalToLocalMap(numberOfGlobalSolverRows)%localNumber(1)=numberOfLocalSolverRows
            rowDomainMapping%globalToLocalMap(numberOfGlobalSolverRows)%domainNumber(1)=rank
            rowDomainMapping%globalToLocalMap(numberOfGlobalSolverRows)%localType(1)=DOMAIN_LOCAL_INTERNAL
            !If this is my rank then set up the solver->equations and equations->solver row mappings
            IF(rank==myrank) THEN
              !Set the interface column/row -> solver row mappings
              !Note that for populating solverMapping%solverRowToEquationsMaps(numberOfLocalSolverRows)%ptr%rowColNumber(i)
              !If the row are equations set rows this is the i'th row number that the solver row is mapped to.
              
              !If the rows are interface rows (which is the case here) then this is the i'th column number that the solver
              !row is mapped to.
              
              !Initialise
              NULLIFY(solverMapping%solverRowToEquationsMaps(numberOfLocalSolverRows)%ptr)
              CALL SolverMappingSRToEQSMap_Initialise(solverMapping%solverRowToEquationsMaps(numberOfLocalSolverRows)%ptr, &
                & err,error,*999)
              ALLOCATE(solverMapping%solverRowToEquationsMaps(numberOfLocalSolverRows)%ptr%rowColNumber(1),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate solver row to equations rows row/col number.",err,error,*999)
              ALLOCATE(solverMapping%solverRowToEquationsMaps(numberOfLocalSolverRows)%ptr%couplingCoefficients(1),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate solver row to equations rows coupling coefficients.",err,error,*999)
              !Set up the interface column -> solver row mappings
              !/todo the solverRowToEquationsMaps may need to be renamed for clarity
              solverMapping%solverRowToEquationsMaps(numberOfLocalSolverRows)%ptr%interfaceConditionIndex=interfaceConditionIdx
              solverMapping%solverRowToEquationsMaps(numberOfLocalSolverRows)%ptr%rowColNumber(1)=localColumn
              solverMapping%solverRowToEquationsMaps(numberOfLocalSolverRows)%ptr%couplingCoefficients(1)=1.0_DP
              !Set the mappings
              ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                & interfaceColToSolverRowsMap(localColumn)%rowCols(1),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate interface column to solver rows map solver row columns.", &
                & err,error,*999)
              ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                & interfaceColToSolverRowsMap(localColumn)%couplingCoefficients(1),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate interface column to solver rows map coupling coefficients.", &
                & err,error,*999)
              solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                & interfaceColToSolverRowsMap(localColumn)%numberOfRowCols=1
              solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                & interfaceColToSolverRowsMap(localColumn)%rowCols(1)=numberOfLocalSolverRows
              solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                & interfaceColToSolverRowsMap(localColumn)%couplingCoefficients(1)=1.0_DP
              SELECT CASE(interfaceCondition%method)
              CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                !Set up the solver row <-> interface row mappings
                !Penalty matrix is the last interface matrix
                interfaceMatrixIdx=interfaceMapping%numberOfInterfaceMatrices
                ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                  & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr% &
                  & interfaceRowToSolverRowsMap(equationsRow)%rowCols(1),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate interface row to solver rows map solver row columns.",err,error,*999)
                ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                  & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr% &
                  & interfaceRowToSolverRowsMap(equationsRow)%couplingCoefficients(1),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate interface row to solver rows map coupling coefficients.", &
                  & err,error,*999)
                !Set the mappings
                solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                  & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr% &
                  & interfaceRowToSolverRowsMap(localColumn)%numberOfRowCols=1
                solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                  & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr% &
                  & interfaceRowToSolverRowsMap(localColumn)%rowCols(1)=numberOfLocalSolverRows
                solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                  & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr% &
                  & interfaceRowToSolverRowsMap(localColumn)%couplingCoefficients(1)=1.0_DP
              ENDSELECT
            ENDIF !rank==my rank
          ELSE IF(constrainedDOF) THEN
            CALL FlagError("Constrained DOFs have not been implemented for Lagrange variables.",err,error,*999)
          ELSE
            !Set the interface column/row -> solver row mappings
            IF(rank==myrank) THEN
              !Set up the solver row <-> interface column mappings
              CALL MatrixRowColCoupling_Initialise(solverMapping% &
                & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                & interfaceColToSolverRowsMap(localColumn),err,error,*999)
              SELECT CASE(interfaceCondition%method)
              CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                !Set up the solver row <-> interface row mappings
                !Penalty matrix is the last interface matrix
                interfaceMatrixIdx=interfaceMapping%numberOfInterfaceMatrices
                !Set the mappings
                CALL MatrixRowColCoupling_Initialise(solverMapping% &
                  & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                  & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr% &
                  & interfaceRowToSolverRowsMap(localColumn),err,error,*999)
              ENDSELECT
            ENDIF
          ENDIF
        ENDDO !globalRowIdx
        IF(ALLOCATED(rankGlobalRowsList)) DEALLOCATE(rankGlobalRowsList)
      ENDDO !interfaceConditionIdx
    ENDDO !rank
    CALL SolverMappingDOFCouplings_Finalise(rowCouplings,err,error,*999)

    CALL DomainMapping_LocalFromGlobalCalculate(rowDomainMapping,err,error,*999)
    !
    !--- Column mappings ---
    !
    ! 4. Calculate the number of local and global columns in the solver matrix. Do this by calculating the list of columns
    !    for each rank so that we can determine the column numbers in rank order.
    !
    !Allocate solver matrix to equations mapping array
    ALLOCATE(solverMapping%solverMatricesToEquationsMaps(solverMapping%numberOfSolverMatrices),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver mapping solver matrices to equations maps.",err,error,*999)
    CALL SolverMappingDOFCouplings_Initialise(columnCouplings,err,error,*999)

    !Calculate the column mappings for each solver matrix
    DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices

      !Initialise the solver matrix to equations map
      NULLIFY(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr)
      CALL SolverMappingSMToEquationsMap_Initialise(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr,err,error,*999)
      !Initialise the variables list
      NULLIFY(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList)
      CALL SolverMappingVariables_Initialise(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList, &
        & err,error,*999)
      !
      ! 4a Calculate the list of field variables involved in the columns of the solver matrix
      !
      !Compute the order of variables for the solver matrices
      CALL List_DetachAndDestroy(createValuesCache%equationsVariableList(solverMatrixIdx)%ptr,numberOfEquationsVariables, &
        & equationsVariables,err,error,*999)
      CALL List_DetachAndDestroy(createValuesCache%interfaceVariableList(solverMatrixIdx)%ptr,numberOfInterfaceVariables, &
        & interfaceVariables,err,error,*999)
      ALLOCATE(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList%variables( &
        & numberOfEquationsVariables+numberOfInterfaceVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate variables list variables.",err,error,*999)
      solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList%numberOfVariables= &
        & numberOfEquationsVariables+numberOfInterfaceVariables
      ALLOCATE(variablesList(numberOfEquationsVariables+numberOfInterfaceVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate variables list.",err,error,*999)
      ALLOCATE(variableProcessed(numberOfEquationsVariables+numberOfInterfaceVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate variable processed.",err,error,*999)
      variableProcessed=.FALSE.
      solverVariableIdx=0
      DO variableIdx=1,numberOfEquationsVariables
        solverVariableIdx=solverVariableIdx+1
        NULLIFY(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList%variables(solverVariableIdx)%ptr)
        CALL SolverMappingVariable_Initialise(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
          & variablesList%variables(solverVariableIdx)%ptr,err,error,*999)
        equationsSetIdx=equationsVariables(1,variableIdx)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        variableType=equationsVariables(2,variableIdx)
        NULLIFY(variable)
        CALL Field_VariableGet(dependentField,variableType,variable,err,error,*999)
        solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList%variables(solverVariableIdx)%ptr% &
          & variable=>variable
        solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList%variables(solverVariableIdx)%ptr% &
          & variableType=variableType
        NULLIFY(variablesList(solverVariableIdx)%ptr)
        CALL List_CreateStart(variablesList(solverVariableIdx)%ptr,err,error,*999)
        CALL List_DataTypeSet(variablesList(solverVariableIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
        CALL List_DataDimensionSet(variablesList(solverVariableIdx)%ptr,3,err,error,*999)
        CALL List_KeyDimensionSet(variablesList(solverVariableIdx)%ptr,1,err,error,*999)
        CALL List_CreateFinish(variablesList(solverVariableIdx)%ptr,err,error,*999)
      ENDDO !variableIdx
      IF(ALLOCATED(equationsVariables)) DEALLOCATE(equationsVariables)
      DO variableIdx=1,numberOfInterfaceVariables
        solverVariableIdx=solverVariableIdx+1
        NULLIFY(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList%variables(solverVariableIdx)%ptr)
        CALL SolverMappingVariable_Initialise(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList% &
          & variables(solverVariableIdx)%ptr,err,error,*999)
        interfaceConditionIdx=interfaceVariables(1,variableIdx)
        NULLIFY(interfaceCondition)
        CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
        NULLIFY(lagrangeField)
        CALL InterfaceCondition_LagrangeFieldGet(interfaceCondition,lagrangeField,err,error,*999)
        variableType=interfaceVariables(2,variableIdx)
        NULLIFY(variable)
        CALL Field_VariableGet(lagrangeField,variableType,variable,err,error,*999)
        solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList%variables(solverVariableIdx)%ptr% &
          & variable=>variable
        solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList%variables(solverVariableIdx)%ptr% &
          & variableType=variableType
        NULLIFY(variablesList(solverVariableIdx)%ptr)
        CALL List_CreateStart(variablesList(solverVariableIdx)%ptr,err,error,*999)
        CALL List_DataTypeSet(variablesList(solverVariableIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
        CALL List_DataDimensionSet(variablesList(solverVariableIdx)%ptr,3,err,error,*999)
        CALL List_KeyDimensionSet(variablesList(solverVariableIdx)%ptr,1,err,error,*999)
        CALL List_CreateFinish(variablesList(solverVariableIdx)%ptr,err,error,*999)
      ENDDO !variableIdx
      IF(ALLOCATED(interfaceVariables)) DEALLOCATE(interfaceVariables)
      !Handle the RHS variables
      CALL List_DetachAndDestroy(createValuesCache%equationsRHSVariableList,numberOfEquationsRHSVariables, &
        & equationsRHSVariables,err,error,*999)
      NULLIFY(solverMapping%rhsVariablesList)
      CALL SolverMappingVariables_Initialise(solverMapping%rhsVariablesList,err,error,*999)
      solverMapping%rhsVariablesList%numberOfVariables=numberOfEquationsRHSVariables
      ALLOCATE(solverMapping%rhsVariablesList%variables(numberOfEquationsRHSVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate RHS variables list.",err,error,*999)
      DO variableIdx=1,numberOfEquationsRHSVariables
        NULLIFY(solverMapping%rhsVariablesList%variables(variableIdx)%ptr)
        CALL SolverMappingVariable_Initialise(solverMapping%rhsVariablesList%variables(variableIdx)%ptr,err,error,*999)
        equationsSetIdx=equationsRHSVariables(1,variableIdx)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        rhsVariableType=equationsRHSVariables(2,variableIdx)
        NULLIFY(rhsVariable)
        CALL Field_VariableGet(dependentField,rhsVariableType,rhsVariable,err,error,*999)
        solverMapping%rhsVariablesList%variables(variableIdx)%ptr%variable=>rhsVariable
        solverMapping%rhsVariablesList%variables(variableIdx)%ptr%variableType=rhsVariableType
!!TODO: Redo this properly with left hand side mapping. For now just add this equations set.
        solverMapping%rhsVariablesList%variables(variableIdx)%ptr%numberOfEquations=0
      ENDDO !variableIdx
      IF(ALLOCATED(equationsRHSVariables)) DEALLOCATE(equationsRHSVariables)
      !Set up solver matrix to equations mapping array
      solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%solverMatrixNumber=solverMatrixIdx
      solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%solverMapping=>solverMapping
      !Allocate the solver matrix to equations set maps array
      ALLOCATE(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
        & solverMatrixToEquationsSetMaps(solverMapping%numberOfEquationsSets),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate solver matrix to equations map solver matrix to equation set maps.", &
        & err,error,*999)
      solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%numberOfEquationsSets=solverMapping%numberOfEquationsSets
      !Allocate the solver matrix to interface conditions maps array
      ALLOCATE(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
        & solverMatrixToInterfaceConditionMaps(solverMapping%numberOfInterfaceConditions),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate solver matrix to equations map solver matrix to interface condition maps.", &
        & err,error,*999)
      solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%numberOfInterfaceConditions= &
        & solverMapping%numberOfInterfaceConditions
      !Presort the column numbers by rank.
      !rankGlobalColumnLists(dofType, equationsIdx, variableIdx, rankIdx)
      !dofType is 1 for domain local DOFs and 2 for ghost DOFs
      ALLOCATE(rankGlobalColumnLists(2,solverMapping%numberOfEquationsSets+solverMapping%numberOfInterfaceConditions, &
        & solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList%numberOfVariables, &
        & 0:numberOfGroupComputationNodes-1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate rank global columns lists.",err,error,*999)
      DO rank=0,numberOfGroupComputationNodes-1
        DO solverVariableIdx=1,solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList%numberOfVariables
          DO equationsIdx=1,solverMapping%numberOfEquationsSets+solverMapping%numberOfInterfaceConditions
            DO dofType=1,2
              NULLIFY(rankGlobalColumnLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr)
              CALL List_CreateStart(rankGlobalColumnLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr,err,error,*999)
              CALL List_DataTypeSet(rankGlobalColumnLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr,LIST_INTG_TYPE, &
                & err,error,*999)
              !Set approximate size for the number of columns per variable.
              CALL List_InitialSizeSet(rankGlobalColumnLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr, &
                & solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList%variables(solverVariableIdx)% &
                & ptr%variable%totalNumberOfDofs,err,error,*999)
              CALL List_DataDimensionSet(rankGlobalColumnLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr,5,err,error,*999)
              CALL List_KeyDimensionSet(rankGlobalColumnLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr,1,err,error,*999)
              CALL List_CreateFinish(rankGlobalColumnLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr,err,error,*999)
            ENDDO !dofType
          ENDDO !equationsIdx
        ENDDO !solverVariableIdx
      ENDDO !rank

      !Allocate sub-matrix information
      !subMatrixInformation(1,equationsIdx,variableIdx) = The equations type, see SolverMapping_EquationsTypes
      !subMatrixInformation(2,equationsIdx,variableIdx) = The equations set or interface condition index
      !subMatrixInformation(3,equationsIdx,variableIdx) = The interface matrix index, or 0 for an equations set matrix
      !equationsIdx goes from 1 to the number of equations sets + interface conditions
      !variableIdx goes from 1 to the number of variables mapped to this solver matrix
      ALLOCATE(subMatrixInformation(3,solverMapping%numberOfEquationsSets+solverMapping%numberOfInterfaceConditions, &
        & solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList%numberOfVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate sub matrix information.",err,error,*999)
      subMatrixInformation=0
      !Allocate sub-matrix list information
      ALLOCATE(subMatrixList(0:3,solverMapping%numberOfEquationsSets+solverMapping%numberOfInterfaceConditions, &
        & solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList%numberOfVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate sub matrix list.",err,error,*999)
      subMatrixList=0

      !
      ! 4b Calculate the number of columns
      !
      
      !Calculate the number of solver dofs
      ALLOCATE(numberOfVariableGlobalSolverDOFs(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
        & variablesList%numberOfVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate number of global solver dofs.",err,error,*999)
      ALLOCATE(numberOfVariableLocalSolverDOFs(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
        & variablesList%numberOfVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate number of local solver dofs.",err,error,*999)
      ALLOCATE(totalNumberOfVariableLocalSolverDOFs(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
        & variablesList%numberOfVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate total number of local solver dofs.",err,error,*999)

      numberOfVariableGlobalSolverDOFs=0
      numberOfVariableLocalSolverDOFs=0
      totalNumberOfVariableLocalSolverDOFs=0
      
      equationsIdx=0
      !Loop over the equations sets
      DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
        equationsIdx=equationsIdx+1
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(equations)
        CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        NULLIFY(vectorMapping)
        CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
        NULLIFY(dynamicMapping)
        CALL EquationsMappingVector_DynamicMappingExists(vectorMapping,dynamicMapping,err,error,*999)
        NULLIFY(linearMapping)
        CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
        NULLIFY(nonlinearMapping)
        CALL EquationsMappingVector_NonlinearMappingExists(vectorMapping,nonlinearMapping,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        NULLIFY(equationsSetVariableList)
        CALL List_CreateStart(equationsSetVariableList,err,error,*999)
        CALL List_DataTypeSet(equationsSetVariableList,LIST_INTG_TYPE,err,error,*999)
        CALL List_DataDimensionSet(equationsSetVariableList,3,err,error,*999)
        CALL List_KeyDimensionSet(equationsSetVariableList,1,err,error,*999)
        CALL List_CreateFinish(equationsSetVariableList,err,error,*999)
        IF(ASSOCIATED(dynamicMapping)) THEN
          equationsVariableListItem(1)=createValuesCache%dynamicVariableType(equationsSetIdx)
          equationsVariableListItem(2)=SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX
          equationsVariableListItem(3)=0
          CALL List_ItemAdd(equationsSetVariableList,equationsVariableListItem,err,error,*999)
        ENDIF
        IF(ASSOCIATED(linearMapping)) THEN
          DO variableIdx=1,solverMapping%createValuesCache%matrixVariableTypes(0,equationsSetIdx,solverMatrixIdx)
            equationsVariableListItem(1)=createValuesCache%matrixVariableTypes(variableIdx,equationsSetIdx,solverMatrixIdx)
            equationsVariableListItem(2)=SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX
            equationsVariableListItem(3)=variableIdx
            CALL List_ItemAdd(equationsSetVariableList,equationsVariableListItem,err,error,*999)
          ENDDO !variableIdx
        ENDIF
        IF(ASSOCIATED(nonlinearMapping)) THEN
          CALL EquationsMappingNonlinear_NumberOfResidualsGet(nonlinearMapping,numberOfResiduals,err,error,*999)
          DO residualIdx=1,numberOfResiduals
            DO variableIdx=1,createValuesCache%residualVariableTypes(0,residualIdx,equationsSetIdx)
              equationsVariableListItem(1)=createValuesCache%residualVariableTypes(variableIdx,residualIdx,equationsSetIdx)
              equationsVariableListItem(2)=SOLVER_MAPPING_EQUATIONS_NONLINEAR_MATRIX
              equationsVariableListItem(3)=variableIdx
              CALL List_ItemAdd(equationsSetVariableList,equationsVariableListItem,err,error,*999)
            ENDDO !variableIdx
          ENDDO !residualIdx
        ENDIF
        CALL List_RemoveDuplicates(equationsSetVariableList,err,error,*999)
        CALL List_DetachAndDestroy(equationsSetVariableList,numberOfVariables,equationsSetVariables,err,error,*999)
        !Initialise equations set to solver matrices equations matrices to solver matrix map
        NULLIFY(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
          & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr)
        CALL SolverMappingEMSToSMMap_Initialise(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
          & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr,err,error,*999)
        solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
          & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%solverMatrixNumber=solverMatrixIdx
        !Allocate the equations set to solver map variables arrays
        ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
          & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%variableTypes(numberOfVariables),STAT=err)        
        IF(err/=0) CALL FlagError("Could not allocate equations matrices to solver matrix maps variable types.",err,error,*999)
        ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
          & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%variables(numberOfVariables),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate equations matrices to solver matrix maps variables.",err,error,*999)
        ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
          & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%varDOFToSolverDOFsMaps(numberOfVariables),STAT=err)
        IF(err/=0)  &
          & CALL FlagError("Could not allocate equations matrices to solver matrix maps variable to solver col maps.", &
          & err,error,*999)
        !Setup
        solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
          & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfVariables=numberOfVariables
        numberOfDynamicMatrices=0
        numberOfLinearMatrices=0
        !Loop over the variables in this equations set.
        DO variableIdx=1,numberOfVariables
          variableType=equationsSetVariables(1,variableIdx)
          matricesType=equationsSetVariables(2,variableIdx)
          matrixVariableIdx=equationsSetVariables(3,variableIdx)
          NULLIFY(dependentVariable)
          CALL Field_VariableGet(dependentField,variableType,dependentVariable,err,error,*999)
          !Find the variable in the list of solver variables
          found=.FALSE.
          DO variablePositionIdx=1,solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList%numberOfVariables
            IF(ASSOCIATED(dependentVariable,solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList% &
              & variables(variablePositionIdx)%ptr%variable)) THEN
              found=.TRUE.
              EXIT
            ENDIF
          ENDDO !variablePositionIdx
          IF(.NOT.found) THEN
            localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
              & " in equations set index "//TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
              & " does not exist in the list of solver variables."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          !Add the equations set variable to the list of equations involving the solver variable
          variableListItem(1)=equationsSetIdx
          variableListItem(2)=variableType
          variableListItem(3)=SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET
          CALL List_ItemAdd(variablesList(variablePositionIdx)%ptr,variableListItem,err,error,*999)
          NULLIFY(columnDOFsMapping)
          CALL FieldVariable_DomainMappingGet(dependentVariable,columnDOFsMapping,err,error,*999)
          NULLIFY(boundaryConditionsVariable)
          CALL BoundaryConditions_VariableGet(boundaryConditions,dependentVariable,boundaryConditionsVariable,err,error,*999)
          solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%variableTypes(variableIdx)=variableType
          solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%variables(variableIdx)%ptr=>dependentVariable
          !Allocate the variable DOF to solver DOFs maps arrays
          NULLIFY(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%varDOFToSolverDOFsMaps(variableIdx)%ptr)
          CALL SolverMappingVDOFToSDOFsMap_Initialise(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%varDOFToSolverDOFsMaps(variableIdx)%ptr,err,error,*999)
          ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%varDOFToSolverDOFsMaps(variableIdx)%ptr% &
            & dofNumbers(dependentVariable%totalNumberOfDofs),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate variable DOF to solver DOFs maps DOF numbers.",err,error,*999)
          ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%varDOFToSolverDOFsMaps(variableIdx)%ptr% &
            & couplingCoefficients(dependentVariable%totalNumberOfDofs),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate variable DOF to solver DOFs maps coupling coefficients.",err,error,*999)
          ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%varDOFToSolverDOFsMaps(variableIdx)%ptr% &
            & additiveConstants(dependentVariable%totalNumberOfDofs),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate variable DOF to solver DOFs maps additive constants.",err,error,*999)
          !Setup
          !Set the sub-matrix information
          subMatrixInformation(1,equationsIdx,variablePositionIdx)=SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET
          subMatrixInformation(2,equationsIdx,variablePositionIdx)=equationsSetIdx
          !Set the sub-matrix lists
          IF(ASSOCIATED(dynamicMapping)) THEN
            numberOfDynamicMatrices=numberOfDynamicMatrices+ &
              & dynamicMapping%varToEquationsMatricesMap%numberOfEquationsMatrices
            IF(dynamicMapping%varToEquationsMatricesMap%numberOfEquationsMatrices>0) THEN
              subMatrixList(0,equationsIdx,variablePositionIdx)=subMatrixList(0,equationsIdx,variablePositionIdx)+1
              subMatrixList(subMatrixList(0,equationsIdx,variablePositionIdx),equationsIdx,variablePositionIdx)= &
                & SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX
            ENDIF
          ENDIF
          IF(ASSOCIATED(linearMapping)) THEN
            variableIdx2=linearMapping%linearVariableTypesMap(variableType)
            numberOfLinearMatrices=numberOfLinearMatrices+ &
              & linearMapping%varToEquationsMatricesMaps(variableIdx2)%ptr%numberOfEquationsMatrices
            IF(linearMapping%varToEquationsMatricesMaps(variableIdx2)%ptr%numberOfEquationsMatrices>0) THEN
              subMatrixList(0,equationsIdx,variablePositionIdx)=subMatrixList(0,equationsIdx,variablePositionIdx)+1
              subMatrixList(subMatrixList(0,equationsIdx,variablePositionIdx),equationsIdx,variablePositionIdx)= &
                & SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX
            ENDIF
          ENDIF
          IF(ASSOCIATED(nonlinearMapping)) THEN
            DO residualIdx=1,nonlinearMapping%numberOfResiduals
              NULLIFY(residualMapping)
              CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
              DO variableIdx2=1,residualMapping%numberOfVariables
                IF(residualMapping%variableTypes(variableIdx2)==variableType) THEN
                  subMatrixList(0,equationsIdx,variablePositionIdx)=subMatrixList(0,equationsIdx,variablePositionIdx)+1
                  subMatrixList(subMatrixList(0,equationsIdx,variablePositionIdx),equationsIdx,variablePositionIdx)= &
                    & SOLVER_MAPPING_EQUATIONS_NONLINEAR_MATRIX
                ENDIF
              ENDDO !variableIdx2
            ENDDO !residualIdx
          ENDIF
          !Loop over the global dofs for this variable.
          DO globalDOF=1,dependentVariable%numberOfGlobalDofs
            DO rankIdx=1,columnDOFsMapping%globalToLocalMap(globalDOF)%numberOfDomains
              localDOF=columnDOFsMapping%globalToLocalMap(globalDOF)%localNumber(rankIdx)
              dofType=columnDOFsMapping%globalToLocalMap(globalDOF)%localType(rankIdx)
              columnRank=columnDOFsMapping%globalToLocalMap(globalDOF)%domainNumber(rankIdx)
              includeColumn=(boundaryConditionsVariable%DOFTypes(globalDOF)==BOUNDARY_CONDITION_DOF_FREE)
              constrainedDOF=(boundaryConditionsVariable%DOFTypes(globalDOF)==BOUNDARY_CONDITION_DOF_CONSTRAINED)
              globalDofCouplingNumber=0
              NULLIFY(dofConstraints)
              CALL BoundaryConditionsVariable_DOFConstraintsExists(boundaryConditionsVariable,dofConstraints,err,error,*999)
              IF(ASSOCIATED(boundaryConditionsVariable%dofConstraints)) THEN
                IF(dofConstraints%numberOfConstraints>0) THEN
                  NULLIFY(dofCoupling)
                  CALL BoundaryConditionsVariableDOFConstraints_DOFCouplingExists(dofConstraints,globalDOF,dofCoupling, &
                    & err,error,*999)
                  IF(ASSOCIATED(dofCoupling)) &
                    & CALL SolverDofCouplings_AddCoupling(columnCouplings,dofCoupling,globalDofCouplingNumber,err,error,*999)
                ENDIF
              ENDIF
              columnListItem(1)=globalDOF
              columnListItem(2)=localDOF
              columnListItem(5)=globalDofCouplingNumber
              IF(dofType/=DOMAIN_LOCAL_GHOST) THEN
                !DOF is not a ghost dof
                IF(includeColumn) THEN
                  columnListItem(3)=1
                  IF(.NOT.variableProcessed(variablePositionIdx)) THEN
                    numberOfVariableGlobalSolverDOFs(variablePositionIdx)= &
                      & numberOfVariableGlobalSolverDOFs(variablePositionIdx)+1
                    IF(columnRank==myrank) THEN
                      numberOfVariableLocalSolverDOFs(variablePositionIdx)= &
                        & numberOfVariableLocalSolverDOFs(variablePositionIdx)+1
                      totalNumberOfVariableLocalSolverDOFs(variablePositionIdx)= &
                        & totalNumberOfVariableLocalSolverDOFs(variablePositionIdx)+1
                    ENDIF
                  ENDIF
                ELSE IF(constrainedDOF) THEN
                  columnListItem(3)=2
                ELSE
                  columnListItem(3)=0
                ENDIF
                columnListItem(4)=variableIdx
                CALL List_ItemAdd(rankGlobalColumnLists(1,equationsIdx,variablePositionIdx,columnRank)%ptr,columnListItem, &
                  & err,error,*999)
              ELSE
                !DOF is a ghost dof
                IF(includeColumn) THEN
                  columnListItem(3)=1
                  IF(.NOT.variableProcessed(variablePositionIdx)) THEN
                    IF(columnRank==myrank) totalNumberOfVariableLocalSolverDOFs(variablePositionIdx)= &
                      & totalNumberOfVariableLocalSolverDOFs(variablePositionIdx)+1
                  ENDIF
                ELSE IF(constrainedDOF) THEN
                  columnListItem(3)=2
                ELSE
                  columnListItem(3)=0
                ENDIF
                columnListItem(4)=variableIdx
                CALL List_ItemAdd(rankGlobalColumnLists(2,equationsIdx,variablePositionIdx,columnRank)%ptr, &
                  & columnListItem,err,error,*999)
              ENDIF
            ENDDO !rankIdx
          ENDDO !globalDOF
          variableProcessed(variablePositionIdx)=.TRUE.
        ENDDO !variableIdx
        IF(ALLOCATED(equationsSetVariables)) DEALLOCATE(equationsSetVariables)
      ENDDO !equationsSetIdx
      !Calculate the number of columns for the interface conditions
      DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
        equationsIdx=equationsIdx+1
        NULLIFY(interfaceCondition)
        CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
        SELECT CASE(interfaceCondition%method)
        CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
          NULLIFY(interfaceEquations)
          CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
          NULLIFY(interfaceMapping)
          CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
          NULLIFY(interfaceDependent)
          CALL InterfaceCondition_InterfaceDependentGet(interfaceCondition,interfaceDependent,err,error,*999)
          !Initialise interface matrices to solver matrix map
          NULLIFY(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr)
          CALL SolverMappingIMSToSMMap_Initialise(solverMapping% &
            & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr,err,error,*999)
          solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%solverMatrixNumber=solverMatrixIdx
          !Allocate the interface to solver map variables arrays
          ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
            & dependentVariableTypes(interfaceMapping%numberOfInterfaceMatrices),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate interface matrices to solver matrix maps dependent variable types.", &
            & err,error,*999)
          ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
            & dependentVariables(interfaceMapping%numberOfInterfaceMatrices),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate interface matrices to solver matrix maps dependent variables.", &
            & err,error,*999)
          ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
            & dependentVarDOFToSolverDOFsMaps(interfaceMapping%numberOfInterfaceMatrices),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate interface matrices to solver matrix maps dependent variables "// &
            & "to solver col maps.",err,error,*999)          
          !First add in the Lagrange to solver variables
          NULLIFY(lagrangeField)
          CALL InterfaceCondition_LagrangeFieldGet(interfaceCondition,lagrangeField,err,error,*999)
          !\todo Lagrange variable type set to the first variable type for now
          variableType=1
          NULLIFY(lagrangeVariable)
          CALL Field_VariableGet(lagrangeField,variableType,lagrangeVariable,err,error,*999)
          !Find the variable in the list of solver variables
          found=.FALSE.
          DO variablePositionIdx=1,solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList%numberOfVariables
            IF(ASSOCIATED(lagrangeVariable,solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList% &
              & variables(variablePositionIdx)%ptr%variable)) THEN
              found=.TRUE.
              EXIT
            ENDIF
          ENDDO !variablePositionIdx
          IF(.NOT.found) THEN
            localError="The Lagrange variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
              & " in interface conditions index "//TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))// &
              & " does not exist in the list of solver variables."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          !Add the interface condition variable to the list of equations involving the solver variable
          variableListItem(1)=interfaceConditionIdx
          variableListItem(2)=variableType
          variableListItem(3)=SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION
          CALL List_ItemAdd(variablesList(variablePositionIdx)%ptr,variableListItem,err,error,*999)
          NULLIFY(columnDOFsMapping)
          CALL FieldVariable_DomainMappingGet(lagrangeVariable,columnDOFsMapping,err,error,*999)
          NULLIFY(boundaryConditionsVariable)
          CALL BoundaryConditions_VariableGet(boundaryConditions,lagrangeVariable,boundaryConditionsVariable,err,error,*999)
          solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%lagrangeVariableType=variableType
          solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%lagrangeVariable=>lagrangeVariable
          !Allocate the variable to solver col maps arrays
          NULLIFY(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%lagrangeVarDOFToSolverDOFsMap)
          CALL SolverMappingVDOFToSDOFsMap_Initialise(solverMapping%interfaceConditionToSolverMatricesMaps( &
            & interfaceConditionIdx)%ptr%interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
            & lagrangeVarDOFToSolverDOFsMap,err,error,*999)
          ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
            & lagrangeVarDOFToSolverDOFsMap%dofNumbers(lagrangeVariable%totalNumberOfDofs),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate Lagrange variable DOF to solver DOFs maps DOF numbers.", &
            & err,error,*999)
          ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
            & lagrangeVarDOFToSolverDOFsMap%couplingCoefficients(lagrangeVariable%totalNumberOfDofs),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate Lagrange variable DOF to solver DOFs maps coupling coefficients.", &
            & err,error,*999)
          ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
            & lagrangeVarDOFToSolverDOFsMap%additiveConstants(lagrangeVariable%totalNumberOfDofs),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate Lagrange variable DOF to solver columnDOFs maps additive constants.", &
            & err,error,*999)
          DO equationsIdx2=1,solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & numberOfEquationsSets
            equationsSetIdx=solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & equationsSetToInterfaceConditionMaps(equationsIdx2)%equationsSetIndex
            interfaceMatrixIdx=solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & equationsSetToInterfaceConditionMaps(equationsIdx2)%interfaceMatrixIndex
            !Set the sub-matrix information
            subMatrixInformation(1,equationsSetIdx,variablePositionIdx)=SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION
            subMatrixInformation(2,equationsSetIdx,variablePositionIdx)=interfaceConditionIdx
            subMatrixInformation(3,equationsSetIdx,variablePositionIdx)=interfaceMatrixIdx
            !Loop over the global dofs for this variable.
            DO globalDOF=1,lagrangeVariable%numberOfGlobalDofs
              DO rankIdx=1,columnDOFsMapping%globalToLocalMap(globalDOF)%numberOfDomains
                localDOF=columnDOFsMapping%globalToLocalMap(globalDOF)%localNumber(rankIdx)
                dofType=columnDOFsMapping%globalToLocalMap(globalDOF)%localType(rankIdx)
                columnRank=columnDOFsMapping%globalToLocalMap(globalDOF)%domainNumber(rankIdx)
                includeColumn=(boundaryConditionsVariable%DOFTypes(globalDOF)==BOUNDARY_CONDITION_DOF_FREE)
                columnListItem(1)=globalDOF
                columnListItem(2)=localDOF
                IF(dofType/=DOMAIN_LOCAL_GHOST) THEN
                  !DOF is not a ghost dof
                  IF(includeColumn) THEN
                    columnListItem(3)=1
                    IF(.NOT.variableProcessed(variablePositionIdx)) THEN
                      numberOfVariableGlobalSolverDOFs(variablePositionIdx)= &
                        & numberOfVariableGlobalSolverDOFs(variablePositionIdx)+1
                      IF(columnRank==myrank) THEN
                        numberOfVariableLocalSolverDOFs(variablePositionIdx)= &
                          & numberOfVariableLocalSolverDOFs(variablePositionIdx)+1
                        totalNumberOfVariableLocalSolverDOFs(variablePositionIdx)= &
                          & totalNumberOfVariableLocalSolverDOFs(variablePositionIdx)+1
                      ENDIF
                    ENDIF
                  ELSE
                    columnListItem(3)=0
                  ENDIF
                  columnListItem(4)=variableIdx
                  CALL List_ItemAdd(rankGlobalColumnLists(1,equationsSetIdx,variablePositionIdx,columnRank)%ptr, &
                    & columnListItem,err,error,*999)
                ELSE
                  !DOF is a ghost dof
                  IF(includeColumn) THEN
                    columnListItem(3)=1
                    IF(.NOT.variableProcessed(variablePositionIdx)) THEN
                      IF(columnRank==myrank) totalNumberOfVariableLocalSolverDOFs(variablePositionIdx)= &
                        & totalNumberOfVariableLocalSolverDOFs(variablePositionIdx)+1
                    ENDIF
                  ELSE
                    columnListItem(3)=0
                  ENDIF
                  columnListItem(4)=variableIdx
                  CALL List_ItemAdd(rankGlobalColumnLists(2,equationsSetIdx,variablePositionIdx,columnRank)%ptr, &
                    & columnListItem,err,error,*999)
                ENDIF
              ENDDO !rankIdx
            ENDDO !globalDOF
            variableProcessed(variablePositionIdx)=.TRUE.
          ENDDO !equationsIdx2
          !Now add in the Dependent variables
          solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfDependentVariables= &
            & interfaceMapping%numberOfInterfaceMatrices
          DO interfaceMatrixIdx=1,interfaceMapping%numberOfInterfaceMatrices
            NULLIFY(dependentVariable)
            CALL InterfaceMapping_MatrixVariableGet(interfaceMapping,interfaceMatrixIdx,dependentVariable,err,error,*999)
            variableType=dependentVariable%variableType
            !Find the variable in the list of solver variables
            found=.FALSE.
            DO variablePositionIdx=1,solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
              & variablesList%numberOfVariables
              IF(ASSOCIATED(dependentVariable,solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
                & variablesList%variables(variablePositionIdx)%ptr%variable)) THEN
                found=.TRUE.
                EXIT
              ENDIF
            ENDDO !variablePositionIdx
            IF(.NOT.found) THEN
              localError="The dependent variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
                & " in interface conditions index "//TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))// &
                & " does not exist in the list of solver variables."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            NULLIFY(equationsSet)
            CALL InterfaceMapping_MatrixEquationsSetGet(interfaceMapping,interfaceMatrixIdx,equationsSet,err,error,*999)
            NULLIFY(columnDOFsMapping)
            CALL FieldVariable_DomainMappingGet(dependentVariable,columnDOFsMapping,err,error,*999)
            NULLIFY(boundaryConditionsVariable)
            CALL BoundaryConditions_VariableGet(boundaryConditions,dependentVariable,boundaryConditionsVariable,err,error,*999)
            !Setup
            solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & dependentVariableTypes(interfaceMatrixIdx)=variableType
            solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & dependentVariables(interfaceMatrixIdx)%ptr=>dependentVariable
            NULLIFY(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & dependentVarDOFToSolverDOFsMaps(interfaceMatrixIdx)%ptr)
            CALL SolverMappingVDOFToSDOFsMap_Initialise(solverMapping%interfaceConditionToSolverMatricesMaps( &
              & interfaceConditionIdx)%ptr%interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & dependentVarDOFToSolverDOFsMaps(interfaceMatrixIdx)%ptr,err,error,*999)
            ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & dependentVarDOFToSolverDOFsMaps(interfaceMatrixIdx)%ptr%dofNumbers(dependentVariable%totalNumberOfDofs), &
              & STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate dependent variable DOF to solver DOFs maps DOF numbers.", &
              & err,error,*999)
            ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & dependentVarDOFToSolverDOFsMaps(interfaceMatrixIdx)%ptr% &
              & couplingCoefficients(dependentVariable%totalNumberOfDofs),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate dependent variable DOF to solver DOFs maps coupling coefficients.", &
              & err,error,*999)
            ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & dependentVarDOFToSolverDOFsMaps(interfaceMatrixIdx)%ptr% &
              & additiveConstants(dependentVariable%totalNumberOfDofs),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate dependent variable DOF to solver DOFs maps additive constants.", &
              & err,error,*999)
            !Set the sub-matrix information
            subMatrixInformation(1,equationsIdx,variablePositionIdx)=SOLVER_MAPPING_EQUATIONS_INTERFACE_TRANSPOSE
            subMatrixInformation(2,equationsIdx,variablePositionIdx)=interfaceConditionIdx
            subMatrixInformation(3,equationsIdx,variablePositionIdx)=interfaceMatrixIdx
            !Loop over the global dofs for this variable.
            DO globalDOF=1,dependentVariable%numberOfGlobalDofs
              DO rankIdx=1,columnDOFsMapping%globalToLocalMap(globalDOF)%numberOfDomains
                localDOF=columnDOFsMapping%globalToLocalMap(globalDOF)%localNumber(rankIdx)
                dofType=columnDOFsMapping%globalToLocalMap(globalDOF)%localType(rankIdx)
                columnRank=columnDOFsMapping%globalToLocalMap(globalDOF)%domainNumber(rankIdx)
                includeColumn=(boundaryConditionsVariable%DOFTypes(globalDOF)==BOUNDARY_CONDITION_DOF_FREE)
                constrainedDOF=(boundaryConditionsVariable%DOFTypes(globalDOF)==BOUNDARY_CONDITION_DOF_CONSTRAINED)
                globalDofCouplingNumber=0
                NULLIFY(dofConstraints)
                CALL BoundaryConditionsVariable_DOFConstraintsExists(boundaryConditionsVariable,dofConstraints,err,error,*999)
                IF(ASSOCIATED(dofConstraints)) THEN
                  IF(dofConstraints%numberOfConstraints>0) THEN
                    NULLIFY(dofCoupling)
                    CALL BoundaryConditionsVariableDOFConstraints_DOFCouplingExists(dofConstraints,globalDOF,dofCoupling, &
                      & err,error,*999)
                    IF(ASSOCIATED(dofCoupling)) &
                      & CALL SolverDofCouplings_AddCoupling(columnCouplings,dofCoupling,globalDofCouplingNumber,err,error,*999)
                  ENDIF
                ENDIF
                columnListItem(1)=globalDOF
                columnListItem(2)=localDOF
                columnListItem(5)=globalDofCouplingNumber
                IF(dofType/=DOMAIN_LOCAL_GHOST) THEN
                  !DOF is not a ghost dof
                  IF(includeColumn) THEN
                    columnListItem(3)=1
                    IF(.NOT.variableProcessed(variablePositionIdx)) THEN
                      numberOfVariableGlobalSolverDOFs(variablePositionIdx)= &
                        & numberOfVariableGlobalSolverDOFs(variablePositionIdx)+1
                      IF(columnRank==myrank) THEN
                        numberOfVariableLocalSolverDOFs(variablePositionIdx)= &
                          & numberOfVariableLocalSolverDOFs(variablePositionIdx)+1
                        totalNumberOfVariableLocalSolverDOFs(variablePositionIdx)= &
                          & totalNumberOfVariableLocalSolverDOFs(variablePositionIdx)+1
                      ENDIF
                    ENDIF
                  ELSE IF(constrainedDOF) THEN
                    columnListItem(3)=2
                  ELSE
                    columnListItem(3)=0
                  ENDIF
                  columnListItem(4)=variableIdx
                  CALL List_ItemAdd(rankGlobalColumnLists(1,equationsIdx,variablePositionIdx, &
                    & columnRank)%ptr,columnListItem,err,error,*999)
                ELSE
                  !DOF is a ghost dof
                  IF(includeColumn) THEN
                    columnListItem(3)=1
                    IF(.NOT.variableProcessed(variablePositionIdx)) THEN
                      IF(columnRank==myrank) totalNumberOfVariableLocalSolverDOFs(variablePositionIdx)= &
                        & totalNumberOfVariableLocalSolverDOFs(variablePositionIdx)+1
                    ENDIF
                  ELSE IF(constrainedDOF) THEN
                    columnListItem(3)=2
                  ELSE
                    columnListItem(3)=0
                  ENDIF
                  columnListItem(4)=variableIdx
                  CALL List_ItemAdd(rankGlobalColumnLists(2,equationsIdx,variablePositionIdx,columnRank)%ptr, &
                    & columnListItem,err,error,*999)
                ENDIF
              ENDDO !rankIdx
            ENDDO !globalDOF
            variableProcessed(variablePositionIdx)=.TRUE.
          ENDDO !matrix_idx
        CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The interface condition method of "//TRIM(NumberToVString(interfaceCondition%method,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !interfaceConditionIdx

      IF(ALLOCATED(variableProcessed)) DEALLOCATE(variableProcessed)

      numberOfLocalSolverDOFs=0
      totalNumberOfLocalSolverDOFs=0
      numberOfGlobalSolverDOFs=0
      DO solverVariableIdx=1,solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList%numberOfVariables
        numberOfLocalSolverDOFs=numberOfLocalSolverDOFs+numberOfVariableLocalSolverDOFs(solverVariableIdx)
        totalNumberOfLocalSolverDOFs=totalNumberOfLocalSolverDOFs+totalNumberOfVariableLocalSolverDOFs(solverVariableIdx)
        numberOfGlobalSolverDOFs=numberOfGlobalSolverDOFs+numberOfVariableGlobalSolverDOFs(solverVariableIdx)
      ENDDO !solverVariableIdx

      !Sanity check
      IF(numberOfLocalSolverDOFs==0) THEN
        localError="Invalid problem setup. The number of local solver DOFs for solver matrix "// &
          & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))//" is zero."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(numberOfGlobalSolverDOFs==0) THEN
        localError="Invalid problem setup. The number of global solver DOFs for solver matrix "// &
          & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))//" is zero."
        CALL FlagError(localError,err,error,*999)
      ENDIF

      !
      ! 4c Set up, allocate and initialise column mappings
      !
      
      !Allocate memory for this solver matrix
      !Allocate solver columns to equations sets maps
      ALLOCATE(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
        & solverDOFToVariableDOFsMaps(totalNumberOfLocalSolverDOFs),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps.",err,error,*999)
      !Set the number of columns
      solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%numberOfColumns=numberOfGlobalSolverDOFs
      !Set the number of variables
      solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%numberOfDofs=numberOfLocalSolverDOFs
      solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%totalNumberOfDofs=totalNumberOfLocalSolverDOFs
      solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%numberOfGlobalDofs=numberOfGlobalSolverDOFs
      !Allocate the columns domain mapping
      CALL DomainMapping_Initialise(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%columnDOFsMapping, &
        & err,error,*999)
      columnDomainMapping=>solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%columnDOFsMapping
      CALL DomainMapping_WorkGroupSet(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%columnDOFsMapping, &
        & workGroup,err,error,*999)
      ALLOCATE(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%columnDOFsMapping% &
        & globalToLocalMap(numberOfGlobalSolverDOFs),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate column dofs mapping global to local.",err,error,*999)
      solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%columnDOFsMapping%numberOfGlobal=numberOfGlobalSolverDOFs
      ALLOCATE(variableRankProcessed(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList% &
        & numberOfVariables,0:numberOfGroupComputationNodes-1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate variable rank processed.",err,error,*999)
      variableRankProcessed=.FALSE.
      !Calculate the column mappings
      numberOfColumns=numberOfGlobalSolverDOFs

      !Initialise
      DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(equations)
        CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        NULLIFY(vectorMapping)
        CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
        NULLIFY(dynamicMapping)
        CALL EquationsMappingVector_DynamicMappingExists(vectorMapping,dynamicMapping,err,error,*999)
        NULLIFY(linearMapping)
        CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
        NULLIFY(nonlinearMapping)
        CALL EquationsMappingVector_NonlinearMappingExists(vectorMapping,nonlinearMapping,err,error,*999)
!!TODO: we could have additinal linear matrices in the equations matrices that have the same field variable as the dynamic variable. Thus it is not just dynamic vs linear but dynamic and linear. 
        IF(ASSOCIATED(dynamicMapping)) THEN
          !Allocate the equations set to solver matrices maps equations matrix to solver matrices maps
          ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & equationsMatrixToSolverMatricesMaps(dynamicMapping%numberOfDynamicMatrices),STAT=err)
          solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%numberOfEquationsMatrices= &
            & dynamicMapping%numberOfDynamicMatrices
          IF(err/=0) &
            & CALL FlagError("Could not allocate equations set to solver matrix map equations matrix to solver matrices maps.", &
            & err,error,*999)
          DO equationMatrixIdx=1,dynamicMapping%numberOfDynamicMatrices
            NULLIFY(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr)
            CALL SolverMappingEMToSMSMap_Initialise(solverMapping% &
              & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr,err,error,*999)
          ENDDO !equationMatrixIdx
          IF(ASSOCIATED(nonlinearMapping)) THEN
            numberOfJacobians=0
            DO residualIdx=1,nonlinearMapping%numberOfResiduals
              NULLIFY(residualMapping)
              CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
              numberOfJacobians=numberOfJacobians+residualMapping%numberOfJacobianMatrices
            ENDDO !residualIdx            
            !Allocate the equations set to solver matrices maps Jacobian matrix to solver matrix maps
            ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%  &
              & jacobianMatrixToSolverMatrixMaps(numberOfJacobians),STAT=err)
            IF(err/=0) &
              & CALL FlagError("Could not allocate equations set to solver matrices map Jacobian matrix to solver matrix maps.", &
              & err,error,*999)
            solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%numberOfJacobianMatrices=numberOfJacobians
            numberOfJacobians=0
            DO residualIdx=1,nonlinearMapping%numberOfResiduals
              NULLIFY(residualMapping)
              CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
              DO equationMatrixIdx=1,residualMapping%numberOfJacobianMatrices
                numberOfJacobians=numberOfJacobians+1
                NULLIFY(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & jacobianMatrixToSolverMatrixMaps(numberOfJacobians)%ptr)
              ENDDO !equationsMatrixIdx
            ENDDO !residualIdx
          ENDIF
        ELSE
          IF(ASSOCIATED(linearMapping)) THEN
            !Allocate the equations set to solver matrices maps equations matrix to solver matrices maps
            ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatrixToSolverMatricesMaps(linearMapping%numberOfLinearMatrices),STAT=err)
            IF(err/=0) &
              & CALL FlagError("Could not allocate equations set to solver matrices map equations matrix to solver matrix maps.", &
              & err,error,*999)
            solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%numberOfEquationsMatrices= &
              & linearMapping%numberOfLinearMatrices
            DO equationMatrixIdx=1,linearMapping%numberOfLinearMatrices
              NULLIFY(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr)
              CALL SolverMappingEMToSMSMap_Initialise(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr,err,error,*999)
            ENDDO !equationMatrixIdx
          ENDIF
          IF(ASSOCIATED(nonlinearMapping)) THEN
            numberOfJacobians=0
            DO residualIdx=1,nonlinearMapping%numberOfResiduals
              NULLIFY(residualMapping)
              CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
              numberOfJacobians=numberOfJacobians+residualMapping%numberOfJacobianMatrices
            ENDDO !residualIdx
            !Allocate the equations set to solver maps for Jacobian matrix (jm) indexing
            ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & jacobianMatrixToSolverMatrixMaps(numberOfJacobians),STAT=err)
            IF(err/=0) &
              & CALL FlagError("Could not allocate equations set to solver matrices map equations matrix to solver matrix maps.", &
              & err,error,*999)
            solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%numberOfJacobianMatrices=numberOfJacobians
            numberOfJacobians=0
            DO residualIdx=1,nonlinearMapping%numberOfResiduals
              NULLIFY(residualMapping)
              CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
              DO equationMatrixIdx=1,residualMapping%numberOfJacobianMatrices
                numberOfJacobians=numberOfJacobians+1
                NULLIFY(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & jacobianMatrixToSolverMatrixMaps(numberOfJacobians)%ptr)
              ENDDO !equationsMatrixIdx
            ENDDO !residualIdx
          ENDIF
        ENDIF
        
        !Initialise solver columns to equations set map
        NULLIFY(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
          & solverMatrixToEquationsSetMaps(equationsSetIdx)%ptr)
        CALL SolverMappingSMToESMap_Initialise(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
          & solverMatrixToEquationsSetMaps(equationsSetIdx)%ptr,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        IF(ASSOCIATED(dynamicMapping)) THEN
          numberOfVariables=1
        ELSE
          IF(ASSOCIATED(nonlinearMapping)) THEN
            !!TODO: Allow for multiple residuals
            numberOfVariables=createValuesCache%residualVariableTypes(0,1,equationsSetIdx)
          ELSE
            numberOfVariables=createValuesCache%matrixVariableTypes(0,equationsSetIdx,solverMatrixIdx)
          ENDIF
        ENDIF

!!TODO: These are not needed or calculated at the moment. Comment out for now
        !Allocate the solver columns to equations set map arrays
        !IF(ASSOCIATED(dynamicMapping)) THEN
        !  solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%solverMatrixToEquationsSetMaps(equationsSetIdx)%ptr% &
        !    & haveDynamic=.TRUE.
        !  ALLOCATE(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
        !    & solverMatrixToEquationsSetMaps(equationsSetIdx)%ptr%solverColToDynamicEquationsMaps(numberOfColumns),STAT=err)
        !  IF(err/=0) CALL FlagError("Could not allocate solver columns to dynamic equations map.",err,error,*999)
        !ENDIF
        !IF(ASSOCIATED(linearMapping)) THEN
        !  solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%solverMatrixToEquationsSetMaps(equationsSetIdx)%ptr% &
        !    & haveLinear=.TRUE.
        !  ALLOCATE(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
        !    & solverMatrixToEquationsSetMaps(equationsSetIdx)%ptr%solverColToLinearEquationsMaps(numberOfColumns),STAT=err)
        !  IF(err/=0) CALL FlagError("Could not allocate solver columns to linear equations map.",err,error,*999)
        !ENDIF
        !IF(ASSOCIATED(nonlinearMapping)) THEN
        !  solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%solverMatrixToEquationsSetMaps(equationsSetIdx)%ptr% &
        !    & haveNonlinear=.TRUE.
        !  ALLOCATE(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
        !    & solverMatrixToEquationsSetMaps(equationsSetIdx)%ptr%solverColToNonlinearEquationsMaps(numberOfColumns),STAT=err)
        !  IF(err/=0) CALL FlagError("Could not allocate solver columns to nonlinear equations map.",err,error,*999)
        !ENDIF
        !Set the solver column to equations set map
        solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%solverMatrixToEquationsSetMaps(equationsSetIdx)%ptr% &
          & equations=>equations
              
        !Allocate the equations matrices to solver matrix maps equations to solver maps
        IF(ASSOCIATED(dynamicMapping)) THEN
          ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%&
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
            & dynamicMatrixToSolverMatrixMaps(dynamicMapping%numberOfDynamicMatrices),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate equations matrices to solver matrix maps dynamic equations "// &
            & "to solver matrix maps.",err,error,*999)
          !Set up dynamic arrays
          solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
            & numberOfDynamicMatrices=dynamicMapping%numberOfDynamicMatrices
          DO equationMatrixIdx=1,dynamicMapping%numberOfDynamicMatrices
            NULLIFY(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%&
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & dynamicMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr)
            CALL SolverMappingEMToSMMap_Initialise(solverMapping% &
              & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & dynamicMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr,err,error,*999)
            solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & dynamicMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr%equationsMatrixType= &
              & SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX
            solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr%equationsMatrixNumber=equationMatrixIdx
            solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr%numberOfSolverMatrices= &
              & solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr%numberOfSolverMatrices+1
          ENDDO !equationMatrixIdx       
          !Set up nonlinear arrays
          IF(ASSOCIATED(nonlinearMapping)) THEN
            numberOfJacobians=0
            DO residualIdx=1,nonlinearMapping%numberOfResiduals
              NULLIFY(residualMapping)
              CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
              numberOfJacobians=numberOfJacobians+residualMapping%numberOfJacobianMatrices
            ENDDO !residualIdx
            solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfJacobianMatrices=numberOfJacobians
            ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & jacobianMatrixToSolverMatrixMaps(numberOfJacobians),STAT=err)            
            IF(err/=0) CALL FlagError("Could not allocate Jacobian matrix to solver matrix maps.",err,error,*999)
            numberOfJacobians=0
            DO residualIdx=1,nonlinearMapping%numberOfResiduals
              NULLIFY(residualMapping)
              CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
              DO equationMatrixIdx=1,residualMapping%numberOfJacobianMatrices
                numberOfJacobians=numberOfJacobians+1
                NULLIFY(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                  & jacobianMatrixToSolverMatrixMaps(numberOfJacobians)%ptr)
                CALL SolverMappingJMToSMMap_Initialise(solverMapping% &
                  & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                  & jacobianMatrixToSolverMatrixMaps(numberOfJacobians)%ptr,err,error,*999)
                solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & jacobianMatrixToSolverMatrixMaps(numberOfJacobians)%ptr=>solverMapping% &
                  & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                  & jacobianMatrixToSolverMatrixMaps(numberOfJacobians)%ptr
              ENDDO !equationsMatrixIdx
            ENDDO !residualIdx
          ENDIF
        ELSE
          IF(ASSOCIATED(linearMapping)) THEN
            ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & linearMatrixToSolverMatrixMaps(numberOfLinearMatrices),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate equations matrices to solver matrix maps linear equations"// &
              & " matrix to solver matrix maps.",err,error,*999)
            !Set up linear arrays
            solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfLinearMatrices=numberOfLinearMatrices
            DO equationMatrixIdx=1,numberOfLinearMatrices
              NULLIFY(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & linearMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr)
              CALL SolverMappingEMToSMMap_Initialise(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & linearMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr,err,error,*999)
              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%equationsMatricesToSolverMatrixMaps( &
                & solverMatrixIdx)%ptr%linearMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                & equationsMatrixType=SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX
              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr%equationsMatrixNumber=equationMatrixIdx
              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr%numberOfSolverMatrices= &
                & solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr%numberOfSolverMatrices+1
            ENDDO !equationMatrixIdx
          ELSE
            ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%linearMatrixToSolverMatrixMaps(0),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate equations matrices to solver matrix maps linear equations "// &
              & "matrix to solver matrix maps.",err,error,*999)
            !Set up linear arrays`
            solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfLinearMatrices=0
          ENDIF
          !Set up nonlinear arrays
          IF(ASSOCIATED(nonlinearMapping)) THEN
            numberOfJacobians=0
            DO residualIdx=1,nonlinearMapping%numberOfResiduals
              NULLIFY(residualMapping)
              CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
              numberOfJacobians=numberOfJacobians+residualMapping%numberOfJacobianMatrices
            ENDDO !residualIdx
            solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfJacobianMatrices=numberOfJacobians
            ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & jacobianMatrixToSolverMatrixMaps(numberOfJacobians),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate Jacobian matrix to solver matrix maps.",err,error,*999)
            numberOfJacobians=0
            DO residualIdx=1,nonlinearMapping%numberOfResiduals
              NULLIFY(residualMapping)
              CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
              DO equationMatrixIdx=1,residualMapping%numberOfJacobianMatrices
                numberOfJacobians=numberOfJacobians+1
                NULLIFY(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                  & jacobianMatrixToSolverMatrixMaps(numberOfJacobians)%ptr)
                CALL SolverMappingJMToSMMap_Initialise(solverMapping% &
                  & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                  & jacobianMatrixToSolverMatrixMaps(numberOfJacobians)%ptr,err,error,*999)
                solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & jacobianMatrixToSolverMatrixMaps(numberOfJacobians)%ptr=> &
                  & solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                  & jacobianMatrixToSolverMatrixMaps(numberOfJacobians)%ptr
              ENDDO !equationsMatrixIdx
            ENDDO !residualIdx
          ENDIF
        ENDIF
        DO variableIdx=1,numberOfVariables
          IF(ASSOCIATED(dynamicMapping)) THEN
            variableType=createValuesCache%dynamicVariableType(equationsSetIdx)
            IF(variableType==dynamicMapping%dynamicVariableType) THEN
              numberOfDynamicMatrices=dynamicMapping%numberOfDynamicMatrices
            ELSE
              numberOfDynamicMatrices=0
            ENDIF
          ELSE
            IF(ASSOCIATED(nonlinearMapping)) THEN
              !!TODO: Allow for multiple residuals
              variableType=createValuesCache%residualVariableTypes(variableIdx,1,equationsSetIdx)
            ELSE
              variableType=createValuesCache%matrixVariableTypes(variableIdx,equationsSetIdx,solverMatrixIdx)
              variableIndex=linearMapping%linearVariableTypesMap(variableType)
              IF(variableIndex==0) THEN
                numberOfLinearMatrices=0
              ELSE
                numberOfLinearMatrices=linearMapping%varToEquationsMatricesMaps(variableIndex)%ptr%numberOfEquationsMatrices
              ENDIF
            ENDIF
          ENDIF

          NULLIFY(dependentVariable)
          CALL Field_VariableGet(dependentField,variableType,dependentVariable,err,error,*999)
          NULLIFY(columnDOFsMapping)
          CALL FieldVariable_DomainMappingGet(dependentVariable,columnDOFsMapping,err,error,*999)
          IF(ASSOCIATED(dynamicMapping)) THEN
            !Allocate dynamic equations to solver matrix maps equations column to solver columns maps
            DO equationMatrixIdx=1,numberOfDynamicMatrices
              matrixNumber=dynamicMapping%varToEquationsMatricesMap%equationsMatrixNumbers(equationMatrixIdx)
              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & dynamicMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr%solverMatrixNumber=solverMatrixIdx
              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & dynamicMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr%equationsMatrixNumber=matrixNumber
              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & dynamicMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                & equationsMatrix=>dynamicMapping%equationsMatrixToVarMaps(matrixNumber)%ptr%equationsMatrix
              numberOfEquationsColumns=dynamicMapping%equationsMatrixToVarMaps(matrixNumber)%ptr%numberOfColumns
              ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & dynamicMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                & equationsColToSolverColsMap(numberOfEquationsColumns),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate dynamic equations equations matrix to solver matrix "// &
                & "column to solver columns map.",err,error,*999)
            ENDDO !equationMatrixIdx
          ELSE
            !Allocate linear equations to solver matrix maps equations column to solver columns maps
            IF(ASSOCIATED(linearMapping)) THEN
              variableIndex=linearMapping%linearVariableTypesMap(variableType)
              DO equationMatrixIdx=1,numberOfLinearMatrices
                matrixNumber=linearMapping%varToEquationsMatricesMaps(variableIndex)%ptr%equationsMatrixNumbers(equationMatrixIdx)
                solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                  & linearMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr%solverMatrixNumber=solverMatrixIdx
                solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                  & linearMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr%equationsMatrixNumber=matrixNumber
                solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                  & linearMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr%equationsMatrix=> &
                  & linearMapping%equationsMatrixToVarMaps(matrixNumber)%ptr%equationsMatrix
                numberOfEquationsColumns=linearMapping%equationsMatrixToVarMaps(matrixNumber)%ptr%numberOfColumns                
                ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                  & linearMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                  & equationsColToSolverColsMap(numberOfEquationsColumns),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate linear equations matrix to solver matrix maps equations "// &
                  & "column to solver columns map.",err,error,*999)
              ENDDO !equationMatrixIdx
            ENDIF
          ENDIF
        ENDDO !variableIdx
        IF(ASSOCIATED(nonlinearMapping)) THEN
          numberOfJacobians=0
          DO residualIdx=1,nonlinearMapping%numberOfResiduals
            NULLIFY(residualMapping)
            CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
            DO equationMatrixIdx=1,residualMapping%numberOfJacobianMatrices
              numberOfJacobians=numberOfJacobians+1
              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & jacobianMatrixToSolverMatrixMaps(numberOfJacobians)%ptr%solverMatrixNumber=solverMatrixIdx
              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & jacobianMatrixToSolverMatrixMaps(numberOfJacobians)%ptr%jacobianMatrixNumber=equationMatrixIdx
              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & jacobianMatrixToSolverMatrixMaps(numberOfJacobians)%ptr%jacobianMatrix=> &
                & residualMapping%jacobianMatrixToVarMaps(equationMatrixIdx)%ptr%jacobian
              numberOfEquationsColumns=residualMapping%jacobianMatrixToVarMaps(equationMatrixIdx)%ptr%numberOfColumns
              ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & jacobianMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                & jacobianColToSolverColsMap(numberOfEquationsColumns),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate Jacobian matrix to solver matrix map Jacobian column "// &
                & "to solver columns map.",err,error,*999)
            ENDDO !equationsMatrixIdx
          ENDDO !residualIdx
        ENDIF
      ENDDO !equationsSetIdx
      
      DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
        NULLIFY(interfaceCondition)
        CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
        NULLIFY(interfaceEquations)
        CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
        NULLIFY(interfaceMapping)
        CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
        SELECT CASE(interfaceCondition%method)
        CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
                
          !Initialise solver columns to interface condition map
          NULLIFY(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
            & solverMatrixToInterfaceConditionMaps(interfaceConditionIdx)%ptr)
          CALL SolverMappingSMToICMap_Initialise(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
            & solverMatrixToInterfaceConditionMaps(interfaceConditionIdx)%ptr,err,error,*999)

!!TODO: Not used at the moment. Comment for now.
          !Allocate the solver columns to equations set map arrays
          !ALLOCATE(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
          !  & solverMatrixToInterfaceConditionMaps(interfaceConditionIdx)%ptr% &
          !  & solverColToInterfaceEquationsMaps(numberOfColumns),STAT=err)
          !IF(err/=0) CALL FlagError("Could not allocate solver columns to interface equations map.",err,error,*999)
          
          !Allocate the interface to solver matrix maps sm interface to solver maps
          ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
            & interfaceMatrixToSolverMatrixMaps(interfaceMapping%numberOfInterfaceMatrices),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate interface matrices to solver matrix maps interface matrix "// &
            & "to solver matrix maps.",err,error,*999)
                
          !Set up interface arrays
          solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfInterfaceMatrices= &
            & interfaceMapping%numberOfInterfaceMatrices
          DO interfaceMatrixIdx=1,interfaceMapping%numberOfInterfaceMatrices
            NULLIFY(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & interfaceMatrixToSolverMatrixMaps(interfaceMatrixIdx)%ptr)
            CALL SolverMappingIMToSMMap_Initialise(solverMapping% &
              & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & interfaceMatrixToSolverMatrixMaps(interfaceMatrixIdx)%ptr,err,error,*999)
            solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr%interfaceMatrixNumber=interfaceMatrixIdx
            solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr%numberOfSolverMatrices=1
            NULLIFY(dependentVariable)
            CALL InterfaceMapping_MatrixVariableGet(interfaceMapping,interfaceMatrixIdx,dependentVariable,err,error,*999)
            solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & dependentVariables(interfaceMatrixIdx)%ptr=>dependentVariable
            solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & interfaceMatrixToSolverMatrixMaps(interfaceMatrixIdx)%ptr%solverMatrixNumber=solverMatrixIdx
            solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & interfaceMatrixToSolverMatrixMaps(interfaceMatrixIdx)%ptr%interfaceMatrixNumber=interfaceMatrixIdx
            solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & interfaceMatrixToSolverMatrixMaps(interfaceMatrixIdx)%ptr%interfaceMatrix=> &
              & interfaceMapping%interfaceMatrixToVarMaps(interfaceMatrixIdx)%ptr%interfaceMatrix
            numberOfInterfaceRows=interfaceMapping%interfaceMatrixToVarMaps(interfaceMatrixIdx)%ptr%numberOfRows
            ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & interfaceMatrixToSolverMatrixMaps(interfaceMatrixIdx)%ptr% &
              & interfaceRowToSolverColsMap(numberOfInterfaceRows),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate interface matrix to solver matrix maps interface row "// &
              & "to solver columns map.",err,error,*999)
          ENDDO !interfaceMatrixIdx
          numberOfInterfaceColumns=interfaceMapping%numberOfColumns
          ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
            & interfaceColToSolverColsMap(numberOfInterfaceColumns),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate interface matrices to solver matrix maps interface column "// &
            & "to solver columns map.",err,error,*999)
        CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The interface condition method of "//TRIM(NumberToVString(interfaceCondition%method,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !interfaceConditionIdx
      
      !Loop over the ranks to ensure that the lowest ranks have the lowest numbered solver variables

      !Allocate dof map to record column reordering
      ALLOCATE(dofMap(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList%numberOfVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate dof map.",err,error,*999)
      DO solverVariableIdx=1,solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList%numberOfVariables
        ALLOCATE(dofMap(solverVariableIdx)%ptr(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList% &
          & variables(solverVariableIdx)%ptr%variable%numberOfGlobalDofs),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate dof map global dof map.",err,error,*999)
        dofMap(solverVariableIdx)%ptr=0
      ENDDO !solverVariableIdx

      ALLOCATE(solverLocalDOF(0:numberOfGroupComputationNodes-1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate solver local dof array.",err,error,*999)
      
      !
      ! 4d Now calculate the solver mappings for each column in rank order
      !
      
      numberOfGlobalSolverDOFs=0
      solverGlobalDOF=0
      solverLocalDOF=0
      DO dofType=1,2
        DO rank=0,numberOfGroupComputationNodes-1
          
          DO solverVariableIdx=1,solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList%numberOfVariables
            
            IF (solverMapping%numberOfInterfaceConditions>0) THEN
              ! Ensure that the dofOffset is calculated as a sum of the number of dofs in the diagonal entries of the solver
              ! matrix (ie the sum of the number of solver dofs in each equation set).
              ! Note that this may not work when running problems in parallel, however, note that interfaces can not currently
              ! be used in parallel either, and the following code only executes if there are interface conditions present.
              tempOffset = 0
              DO tempSolverVariableIdx=1,solverVariableIdx
                DO globalDOF=1,SIZE(dofMap(tempSolverVariableIdx)%ptr)
                  IF (dofMap(tempSolverVariableIdx)%ptr(globalDOF)>0) THEN
                    tempOffset=tempOffset+1
                  ENDIF
                ENDDO !globalDOF
              ENDDO !tempSolverVariableIdx
              globalDOFsOffset=tempOffset
              localDOFsOffset=tempOffset
            ELSE
              globalDOFsOffset=solverGlobalDOF
              localDOFsOffset=solverLocalDOF(rank)
            ENDIF
            
            variableType=solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList% &
              & variables(solverVariableIdx)%ptr%variableType
                  
            DO equationsIdx=1,solverMapping%numberOfEquationsSets+solverMapping%numberOfInterfaceConditions
                    
              !Get columns list
              CALL List_Sort(rankGlobalColumnLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr,err,error,*999)
              CALL List_DetachAndDestroy(rankGlobalColumnLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr, &
                & numberOfRankColumns,rankGlobalColumnsList,err,error,*999)
              
              IF(numberOfRankColumns>0) THEN
                
                solverGlobalDOF=globalDOFsOffset
                solverLocalDOF(rank)=localDOFsOffset
                
                equationType=subMatrixInformation(1,equationsIdx,solverVariableIdx)
                SELECT CASE(equationType)
                CASE(SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET)
                  
                  equationsSetIdx=subMatrixInformation(2,equationsIdx,solverVariableIdx)
                  
                  NULLIFY(equationsSet)
                  CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
                  NULLIFY(equations)
                  CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
                  NULLIFY(vectorEquations)
                  CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                  NULLIFY(vectorMapping)
                  CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
                  NULLIFY(dynamicMapping)
                  CALL EquationsMappingVector_DynamicMappingExists(vectorMapping,dynamicMapping,err,error,*999)
                  NULLIFY(linearMapping)
                  CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
                  NULLIFY(nonlinearMapping)
                  CALL EquationsMappingVector_NonlinearMappingExists(vectorMapping,nonlinearMapping,err,error,*999)
                  NULLIFY(dependentField)
                  CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
                 
                  numberOfDynamicMatrices=0
                  numberOfLinearMatrices=0
                  IF(ASSOCIATED(dynamicMapping)) THEN
                    IF(variableType==dynamicMapping%dynamicVariableType) THEN
                      numberOfDynamicMatrices=dynamicMapping%numberOfDynamicMatrices
                    ELSE
                      numberOfDynamicMatrices=0
                    ENDIF
                  ENDIF
                  IF(ASSOCIATED(linearMapping)) THEN
                    variableIndex=linearMapping%linearVariableTypesMap(variableType)
                    IF(variableIndex==0) THEN
                      numberOfLinearMatrices=0
                    ELSE
                      numberOfLinearMatrices=linearMapping%varToEquationsMatricesMaps(variableIndex)%ptr%numberOfEquationsMatrices
                    ENDIF
                  ENDIF
                  
                  !Loop over the variables
                        
                  dependentVariable=>solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList% &
                    & variables(solverVariableIdx)%ptr%variable
                  NULLIFY(columnDOFsMapping)
                  CALL FieldVariable_DomainMappingGet(dependentVariable,columnDOFsMapping,err,error,*999)
                        
                  DO globalDOFIdx=1,numberOfRankColumns
                    globalDOF=rankGlobalColumnsList(1,globalDOFIdx)
                    localDOF=rankGlobalColumnsList(2,globalDOFIdx)
                    !dofType=rankGlobalColumnsList(3,globalDOFIdx)
                    includeColumn=(rankGlobalColumnsList(3,globalDOFIdx)==1)
                    constrainedDOF=(rankGlobalColumnsList(3,globalDOFIdx)==2)
                    variableIdx=rankGlobalColumnsList(4,globalDOFIdx)
                    globalDofCouplingNumber=rankGlobalColumnsList(5,globalDOFIdx)
                    IF(globalDofCouplingNumber>0) THEN
                      NULLIFY(colEquationsCols)
                      !Fix this!!!!
                      colEquationsCols=>columnCouplings%dofCouplings(globalDofCouplingNumber)%ptr
                      !CALL BoundaryConditionsVariableDOFConstraints_DOFCouplingGet(columnCouplings,globalDOFCouplingNumber, &
                      !  & colEquationsCols,err,error,*999)
                      numberColEquationsCols=colEquationsCols%numberOfDofs
                    ELSE
                      numberColEquationsCols=1
                      dummyDofCoupling%globalDofs(1)=globalDOF
                      dummyDofCoupling%localDofs(1)=localDOF
                      dummyDofCoupling%coefficients(1)=1.0_DP
                      colEquationsCols=>dummyDofCoupling
                    ENDIF

                    IF(includeColumn) THEN
                      !DOF is not fixed so map the variable/equation dof to a new solver dof
                      
                      IF(dofType==2) THEN
                        solverGlobalDOF=dofMap(solverVariableIdx)%ptr(globalDOF)
                      ELSE
                        solverGlobalDOF=solverGlobalDOF+1
                        dofMap(solverVariableIdx)%ptr(globalDOF)=solverGlobalDOF
                      ENDIF
                            
                      solverLocalDOF(rank)=solverLocalDOF(rank)+1
                      
                      IF(rank==myrank) THEN
                        
                        IF(.NOT.variableRankProcessed(solverVariableIdx,rank)) THEN
                          
                          !Set up the column domain mappings.
                          CALL DomainGlobalMapping_Initialise(columnDomainMapping%globalToLocalMap(solverGlobalDOF),err,error,*999)
                          !There are no ghosted cols for the solver matrices so there is only 1 domain for the global to
                          !local map.
                          !Allocate the global to local map arrays
                          ALLOCATE(columnDomainMapping%globalToLocalMap(solverGlobalDOF)%localNumber(1),STAT=err)
                          IF(err/=0) CALL FlagError("Could not allocate column domain global to local map local number.", &
                            & err,error,*999)
                          ALLOCATE(columnDomainMapping%globalToLocalMap(solverGlobalDOF)%domainNumber(1),STAT=err)
                          IF(err/=0) CALL FlagError("Could not allocate column domain global to local map domain number.", &
                            & err,error,*999)
                          ALLOCATE(columnDomainMapping%globalToLocalMap(solverGlobalDOF)%localType(1),STAT=err)
                          IF(err/=0) CALL FlagError("Could not allocate column domain global to local map domain number.", &
                            & err,error,*999)
                          !Set up the global to local mappings
                          columnDomainMapping%globalToLocalMap(solverGlobalDOF)%numberOfDomains=1
                          columnDomainMapping%globalToLocalMap(solverGlobalDOF)%localNumber(1)=solverLocalDOF(rank)
                          columnDomainMapping%globalToLocalMap(solverGlobalDOF)%domainNumber(1)=rank
                          columnDomainMapping%globalToLocalMap(solverGlobalDOF)%localType(1)=DOMAIN_LOCAL_INTERNAL

                          !Set up the solver DOF -> variable DOFs map
                          !Initialise
                          NULLIFY(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
                            & solverDOFToVariableDOFsMaps(solverLocalDOF(rank))%ptr)
                          CALL SolverMappingSDOFToVDOFsMap_Initialise(solverMapping%solverMatricesToEquationsMaps( &
                            & solverMatrixIdx)%ptr%solverDOFToVariableDOFsMaps(solverLocalDOF(rank))%ptr,err,error,*999)
                          !Allocate the solver DOF to variable DOFs arrays
                          ALLOCATE(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
                            & solverDOFToVariableDOFsMaps(solverLocalDOF(rank))%ptr% &
                            & equationTypes(numberColEquationsCols),STAT=err)
                          IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps equations types.", &
                            & err,error,*999)
                          ALLOCATE(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
                            & solverDOFToVariableDOFsMaps(solverLocalDOF(rank))%ptr% &
                            & equationIndices(numberColEquationsCols),STAT=err)
                          IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps equations indices.", &
                            & err,error,*999)
                          ALLOCATE(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
                            & solverDOFToVariableDOFsMaps(solverLocalDOF(rank))%ptr%variable(numberColEquationsCols),STAT=err)
                          IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps variable type.", &
                            & err,error,*999)
                          ALLOCATE(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
                            & solverDOFToVariableDOFsMaps(solverLocalDOF(rank))%ptr% &
                            & variableDOF(numberColEquationsCols),STAT=err)
                          IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps variable dof.", &
                            & err,error,*999)
                          ALLOCATE(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
                            & solverDOFToVariableDOFsMaps(solverLocalDOF(rank))%ptr% &
                            & variableCoefficient(numberColEquationsCols),STAT=err)
                          IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps variable coefficient.", &
                            & err,error,*999)
                          ALLOCATE(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
                            & solverDOFToVariableDOFsMaps(solverLocalDOF(rank))%ptr% &
                            & additiveConstant(numberColEquationsCols),STAT=err)
                          IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps additive constant.", &
                            & err,error,*999)
                          solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%solverDOFToVariableDOFsMaps( &
                            & solverLocalDOF(rank))%ptr%numberOfEquationDOFs=numberColEquationsCols
                          !Set the solver -> equations mappings
                          DO colEquationsColIdx=1,numberColEquationsCols
                            solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
                              & solverDOFToVariableDOFsMaps(solverLocalDOF(rank))%ptr%equationTypes(colEquationsColIdx)= &
                              & SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET
                            solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
                              & solverDOFToVariableDOFsMaps(solverLocalDOF(rank))%ptr%equationIndices(colEquationsColIdx)= &
                              & equationsSetIdx
                            solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
                              & solverDOFToVariableDOFsMaps(solverLocalDOF(rank))%ptr%variable(colEquationsColIdx)%ptr=> &
                              & dependentVariable
                            solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
                              & solverDOFToVariableDOFsMaps(solverLocalDOF(rank))%ptr%variableDOF(colEquationsColIdx)= &
                              & colEquationsCols%localDofs(colEquationsColIdx)
                            solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
                              & solverDOFToVariableDOFsMaps(solverLocalDOF(rank))%ptr%variableCoefficient(colEquationsColIdx)= &
                              & colEquationsCols%coefficients(colEquationsColIdx)
                            solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
                              & solverDOFToVariableDOFsMaps(solverLocalDOF(rank))%ptr%additiveConstant(colEquationsColIdx)=0.0_DP
                          ENDDO !colEquationsColIdx
                        ENDIF
                        !Set up the equations variables -> solver columns mapping
                        DO colEquationsColIdx=1,numberColEquationsCols
                          eqnLocalDof=colEquationsCols%localDofs(colEquationsColIdx)
                          couplingCoefficient=colEquationsCols%coefficients(colEquationsColIdx)
                          solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                            & varDOFToSolverDOFsMaps(variableIdx)%ptr%dofNumbers(eqnLocalDof)=solverGlobalDOF
                          solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                            & varDOFToSolverDOFsMaps(variableIdx)%ptr%couplingCoefficients(eqnLocalDof)=couplingCoefficient
                          solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                            & varDOFToSolverDOFsMaps(variableIdx)%ptr%additiveConstants(eqnLocalDof)=0.0_DP
                          !Set up the equations columns -> solver columns mapping
                          DO matrixTypeIdx=1,subMatrixList(0,equationsIdx,solverVariableIdx)
                            matrixType=subMatrixList(matrixTypeIdx,equationsIdx,solverVariableIdx)
                            SELECT CASE(matrixType)
                            CASE(SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX)
                              !Dynamic matrix
                              DO equationMatrixIdx=1,numberOfDynamicMatrices
                                matrixNumber=dynamicMapping%varToEquationsMatricesMap%equationsMatrixNumbers(equationMatrixIdx)
                                equationsColumn=dynamicMapping%varToEquationsMatricesMap%dofToColumnsMaps(equationMatrixIdx)% &
                                  & columnDOF(eqnLocalDof)
                                !Allocate the equation to solver map column items.
                                ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & dynamicMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                  & equationsColToSolverColsMap(equationsColumn)%rowCols(1),STAT=err)
                                IF(err/=0) CALL FlagError("Could not allocate dynamic equations matrix to solver matrix "//&
                                  & "equations column to solver columns map solver colums rowcols.",err,error,*999)
                                ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & dynamicMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                  & equationsColToSolverColsMap(equationsColumn)%couplingCoefficients(1),STAT=err)
                                IF(err/=0) CALL FlagError("Could not allocate dynamic equations matrix to solver matrix "//&
                                  & "equations column to solver columns map solver colums coupling coefficients.",err,error,*999)
                                solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & dynamicMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                  & equationsColToSolverColsMap(equationsColumn)%numberOfRowCols=1
                                solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & dynamicMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                  & equationsColToSolverColsMap(equationsColumn)%rowCols(1)=solverGlobalDOF
                                solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & dynamicMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                  & equationsColToSolverColsMap(equationsColumn)%couplingCoefficients(1)=couplingCoefficient
                              ENDDO !equationMatrixIdx
                            CASE(SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX)
                              variableIndex=linearMapping%linearVariableTypesMap(variableType)
                              DO equationMatrixIdx=1,numberOfLinearMatrices
                                matrixNumber=linearMapping%varToEquationsMatricesMaps(variableIndex)%ptr% &
                                  & equationsMatrixNumbers(equationMatrixIdx)
                                equationsColumn=linearMapping%varToEquationsMatricesMaps(variableIndex)%ptr% &
                                  & dofToColumnsMaps(equationMatrixIdx)%columnDOF(eqnLocalDof)
                                !Allocate the equation to solver map column items.
                                !No coupling yet so the mapping is 1-1
                                ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & linearMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                  & equationsColToSolverColsMap(equationsColumn)%rowCols(1),STAT=err)
                                IF(err/=0) CALL FlagError("Could not allocate linear equations matrix to solver matrix "//&
                                  & "equations column to solver columns map solver colums rowcols.",err,error,*999)
                                ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & linearMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                  & equationsColToSolverColsMap(equationsColumn)%couplingCoefficients(1),STAT=err)
                                IF(err/=0) CALL FlagError("Could not allocate linear equations matrix to solver matrix "//&
                                  & "equations column to solver columns map solver colums coupling coefficients.",err,error,*999)
                                solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & linearMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                  & equationsColToSolverColsMap(equationsColumn)%numberOfRowCols=1
                                solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & linearMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                  & equationsColToSolverColsMap(equationsColumn)%rowCols(1)=solverGlobalDOF
                                solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & linearMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                  & equationsColToSolverColsMap(equationsColumn)%couplingCoefficients(1)=couplingCoefficient
                              ENDDO !equationMatrixIdx
                            CASE(SOLVER_MAPPING_EQUATIONS_NONLINEAR_MATRIX)
                              residualLoop1: DO residualIdx=1,nonlinearMapping%numberOfResiduals
                                NULLIFY(residualMapping)
                                CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping, &
                                  & err,error,*999)
                                residualVariableLoop1: DO equationMatrixIdx=1,residualMapping%numberOfVariables
                                  IF(residualMapping%varToJacobianMatrixMaps(equationMatrixIdx)%ptr%variableType==variableType) &
                                    & EXIT residualLoop1
                                ENDDO residualVariableLoop1 !equationsMatrixIdx
                              ENDDO residualLoop1 !residualIdx
                              jacobianColumn=residualMapping%varToJacobianMatrixMaps(equationMatrixIdx)%ptr% &
                                & dofToColumnsMap(eqnLocalDof)
                              !Allocate the Jacobian to solver map column items.
                              !No coupling yet so the mapping is 1-1
                              ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                & jacobianMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                & jacobianColToSolverColsMap(jacobianColumn)%rowCols(1),STAT=err)
                              IF(err/=0) CALL FlagError("Could not allocate Jacobian equations matrix to solver matrix "//&
                                  & "Jacobian column to solver columns map solver colums rowcols.",err,error,*999)
                              ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                & jacobianMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                & jacobianColToSolverColsMap(jacobianColumn)%couplingCoefficients(1),STAT=err)
                              IF(err/=0) CALL FlagError("Could not allocate Jacobian equations matrix to solver matrix "//&
                                & "Jacobian column to solver columns map solver colums coupling coefficients.",err,error,*999)
                              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                & jacobianMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                & jacobianColToSolverColsMap(jacobianColumn)%numberOfRowCols=1
                              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                & jacobianMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                & jacobianColToSolverColsMap(jacobianColumn)%rowCols(1)=solverGlobalDOF
                              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                & jacobianMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                & jacobianColToSolverColsMap(jacobianColumn)%couplingCoefficients(1)=couplingCoefficient
                            CASE DEFAULT
                              localError="The equations matrix type of "//TRIM(NumberToVString(matrixType,"*",err,error))// &
                                & " is invalid."
                              CALL FlagError(localError,err,error,*999)
                            END SELECT
                          ENDDO !matrixTypeIdx
                        END DO !colEquationsColIdx
                      ELSE !rank /= myrank
                              
                        IF(.NOT.variableRankProcessed(solverVariableIdx,rank)) THEN
                                
                          !Set up the column domain mappings.
                          CALL DomainGlobalMapping_Initialise(columnDomainMapping%globalToLocalMap(solverGlobalDOF),err,error,*999)
                          !There are no ghosted cols for the solver matrices so there is only 1 domain for the global to
                          !local map.
                          !Allocate the global to local map arrays
                          ALLOCATE(columnDomainMapping%globalToLocalMap(solverGlobalDOF)%localNumber(1),STAT=err)
                          IF(err/=0) CALL FlagError("Could not allocate column domain global to local map local number.", &
                            & err,error,*999)
                          ALLOCATE(columnDomainMapping%globalToLocalMap(solverGlobalDOF)%domainNumber(1),STAT=err)
                          IF(err/=0) CALL FlagError("Could not allocate column domain global to local map domain number.", &
                            & err,error,*999)
                          ALLOCATE(columnDomainMapping%globalToLocalMap(solverGlobalDOF)%localType(1),STAT=err)
                          IF(err/=0) CALL FlagError("Could not allocate column domain global to local map domain number.", &
                            & err,error,*999)
                          !Set up the global to local mappings
                          columnDomainMapping%globalToLocalMap(solverGlobalDOF)%numberOfDomains=1
                          columnDomainMapping%globalToLocalMap(solverGlobalDOF)%localNumber(1)=solverLocalDOF(rank)
                          columnDomainMapping%globalToLocalMap(solverGlobalDOF)%domainNumber(1)=rank
                          columnDomainMapping%globalToLocalMap(solverGlobalDOF)%localType(1)=DOMAIN_LOCAL_INTERNAL
                          
                        ENDIF
                              
                      ENDIF !rank == myrank
                    ELSE IF(constrainedDOF) THEN
                      !Do nothing, this is set up above
                    ELSE
                      IF(rank==myrank) THEN
                        !Set up the equations variables -> solver columns mapping
                        solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%equationsMatricesToSolverMatrixMaps( &
                          & solverMatrixIdx)%ptr%varDOFToSolverDOFsMaps(variableIdx)%ptr%dofNumbers(localDOF)=0
                        solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%equationsMatricesToSolverMatrixMaps( &
                          & solverMatrixIdx)%ptr%varDOFToSolverDOFsMaps(variableIdx)%ptr%couplingCoefficients(localDOF)=0.0_DP
                        solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%equationsMatricesToSolverMatrixMaps( &
                          & solverMatrixIdx)%ptr%varDOFToSolverDOFsMaps(variableIdx)%ptr%additiveConstants(localDOF)=0.0_DP
                        DO matrixTypeIdx=1,subMatrixList(0,equationsIdx,solverVariableIdx)
                          matrixType=subMatrixList(matrixTypeIdx,equationsIdx,solverVariableIdx)
                          SELECT CASE(matrixType)
                          CASE(SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX)
                            !Set up the equations columns -> solver columns mapping
                            DO equationMatrixIdx=1,numberOfDynamicMatrices
                              equationsColumn=dynamicMapping%varToEquationsMatricesMap% &
                                & dofToColumnsMaps(equationMatrixIdx)%columnDOF(localDOF)
                              CALL MatrixRowColCoupling_Initialise(solverMapping% &
                                & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                & dynamicMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                & equationsColToSolverColsMap(equationsColumn),err,error,*999)
                            ENDDO !equationMatrixIdx
                          CASE(SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX)
                            !Set up the equations columns -> solver columns mapping
                            variableIndex=linearMapping%linearVariableTypesMap(variableType)
                            DO equationMatrixIdx=1,numberOfLinearMatrices
                              equationsColumn=linearMapping%varToEquationsMatricesMaps(variableIndex)%ptr% &
                                & dofToColumnsMaps(equationMatrixIdx)%columnDOF(localDOF)
                              CALL MatrixRowColCoupling_Initialise(solverMapping% &
                                & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                & linearMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                & equationsColToSolverColsMap(equationsColumn),err,error,*999)
                            ENDDO !equationMatrixIdx
                          CASE(SOLVER_MAPPING_EQUATIONS_NONLINEAR_MATRIX)
                            residualLoop2: DO residualIdx=1,nonlinearMapping%numberOfResiduals
                              NULLIFY(residualMapping)
                              CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping, &
                                & err,error,*999)
                              residualVariableLoop2: DO equationMatrixIdx=1,residualMapping%numberOfVariables
                                IF(residualMapping%varToJacobianMatrixMaps(equationMatrixIdx)%ptr%variableType==variableType) &
                                  & EXIT residualLoop2
                              ENDDO residualVariableLoop2 !equationsMatrixIdx
                            ENDDO residualLoop2 !residualIdx
                            jacobianColumn=residualMapping%varToJacobianMatrixMaps(equationMatrixIdx)%ptr%dofToColumnsMap(localDOF)
                            CALL MatrixRowColCoupling_Initialise(solverMapping% &
                              & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                              & jacobianMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                              & jacobianColToSolverColsMap(jacobianColumn),err,error,*999)
                          CASE DEFAULT
                            localError="The equations matrix type of "// &
                              & TRIM(NumberToVString(matrixType,"*",err,error))//" is invalid."
                            CALL FlagError(localError,err,error,*999)
                          END SELECT
                        ENDDO !matrixTypeIdx
                      ENDIF !rank == myrank
                    ENDIF !include_column
                  ENDDO !globalDOF
                  
                CASE(SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION,SOLVER_MAPPING_EQUATIONS_INTERFACE_TRANSPOSE)
                        
                  !Now handle the interface condition rows and columns
                  
                  interfaceConditionIdx=subMatrixInformation(2,equationsIdx,solverVariableIdx)
                  interfaceMatrixIdx=subMatrixInformation(3,equationsIdx,solverVariableIdx)

                  NULLIFY(interfaceCondition)
                  CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
                        
                  SELECT CASE(interfaceCondition%method)
                  CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
                    NULLIFY(interfaceEquations)
                    CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
                    NULLIFY(interfaceMapping)
                    CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
                    
                    !Loop over the variables
                    ! This is not only a Lagrange variable (it could be a equationset variable) - rename for clarity.
                    
                    lagrangeVariable=>solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList% &
                      & variables(solverVariableIdx)%ptr%variable
                         
                    DO globalDOFIdx=1,numberOfRankColumns
                      globalDOF=rankGlobalColumnsList(1,globalDOFIdx)
                      localDOF=rankGlobalColumnsList(2,globalDOFIdx)
                      !dofType=rankGlobalColumnsList(3,globalDOFIdx)
                      includeColumn=(rankGlobalColumnsList(3,globalDOFIdx)==1)
                      constrainedDOF=(rankGlobalColumnsList(3,globalDOFIdx)==2)
                      globalDofCouplingNumber=rankGlobalColumnsList(5,globalDOFIdx)
                      
                      IF(globalDofCouplingNumber>0) THEN
                        NULLIFY(colEquationsCols)
                        !Fix this!!!!
                        colEquationsCols=>columnCouplings%dofCouplings(globalDofCouplingNumber)%ptr
                        !CALL BoundaryConditionsVariableDOFConstraints_DOFCouplingGet(columnCouplings,globalDOFCouplingNumber, &
                        !  & colEquationsCols,err,error,*999)
                        numberColEquationsCols=colEquationsCols%numberOfDofs
                      ELSE
                        numberColEquationsCols=1
                        dummyDofCoupling%globalDofs(1)=globalDOF
                        dummyDofCoupling%localDofs(1)=localDOF
                        dummyDofCoupling%coefficients(1)=1.0_DP
                        colEquationsCols=>dummyDofCoupling
                      END IF

                      IF(includeColumn) THEN
                        !DOF is not fixed so map the variable/equation dof to a new solver dof
                              
                        IF(dofType==2) THEN
                          !Ghosted, reuse global dof
                          solverGlobalDOF=dofMap(solverVariableIdx)%ptr(globalDOF)
                        ELSE
                          solverGlobalDOF=solverGlobalDOF+1
                          dofMap(solverVariableIdx)%ptr(globalDOF)=solverGlobalDOF
                        ENDIF
                              
                        solverLocalDOF(rank)=solverLocalDOF(rank)+1
                              
                        IF(rank==myrank) THEN
                              
                          IF(.NOT.variableRankProcessed(solverVariableIdx,rank)) THEN
                                
                            !Set up the column domain mappings.
                            CALL DomainGlobalMapping_Initialise(columnDomainMapping%globalToLocalMap(solverGlobalDOF), &
                              & err,error,*999)
                            !There are no ghosted cols for the solver matrices so there is only 1 domain for the global to
                            !local map.
                            !Allocate the global to local map arrays
                            ALLOCATE(columnDomainMapping%globalToLocalMap(solverGlobalDOF)%localNumber(1),STAT=err)
                            IF(err/=0) CALL FlagError("Could not allocate column domain global to local map local number.", &
                              & err,error,*999)
                            ALLOCATE(columnDomainMapping%globalToLocalMap(solverGlobalDOF)%domainNumber(1),STAT=err)
                            IF(err/=0) CALL FlagError("Could not allocate column domain global to local map domain number.", &
                              & err,error,*999)
                            ALLOCATE(columnDomainMapping%globalToLocalMap(solverGlobalDOF)%localType(1),STAT=err)
                            IF(err/=0) CALL FlagError("Could not allocate column domain global to local map domain number.", &
                              & err,error,*999)
                            !Set up the global to local mappings
                            columnDomainMapping%globalToLocalMap(solverGlobalDOF)%numberOfDomains=1
                            columnDomainMapping%globalToLocalMap(solverGlobalDOF)%localNumber(1)=solverLocalDOF(rank)
                            columnDomainMapping%globalToLocalMap(solverGlobalDOF)%domainNumber(1)=rank
                            columnDomainMapping%globalToLocalMap(solverGlobalDOF)%localType(1)=DOMAIN_LOCAL_INTERNAL
                                  
                            !Set up the solver column -> equations column mappings.
                            !Set up the solver dofs -> variable dofs map
                            !Initialise
                            NULLIFY(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
                              & solverDOFToVariableDOFsMaps(solverLocalDOF(rank))%ptr)
                            CALL SolverMappingSDOFToVDOFsMap_Initialise(solverMapping% &
                              & solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
                              & solverDOFToVariableDOFsMaps(solverLocalDOF(rank))%ptr,err,error,*999)
                            !Allocate the solver dofs to variable dofs arrays
                            ALLOCATE(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
                              & solverDOFToVariableDOFsMaps(solverLocalDOF(rank))%ptr% &
                              & equationTypes(numberColEquationsCols),STAT=err)
                            IF(err/=0) CALL FlagError("Could not allocate solver DOF to variable maps equations types.", &
                              & err,error,*999)
                            ALLOCATE(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
                              & solverDOFToVariableDOFsMaps(solverLocalDOF(rank))%ptr% &
                              & equationIndices(numberColEquationsCols),STAT=err)
                            IF(err/=0) CALL FlagError("Could not allocate solver DOF to variable maps equations indices.", &
                              & err,error,*999)
                            ALLOCATE(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
                              & solverDOFToVariableDOFsMaps(solverLocalDOF(rank))%ptr%variable(numberColEquationsCols),STAT=err)
                            IF(err/=0) CALL FlagError("Could not allocate solver DOF to variable maps variable type.", &
                              & err,error,*999)
                            ALLOCATE(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
                              & solverDOFToVariableDOFsMaps(solverLocalDOF(rank))%ptr%variableDOF(numberColEquationsCols),STAT=err)
                            IF(err/=0) CALL FlagError("Could not allocate solver DOF to variable maps variable DOF.", &
                              & err,error,*999)
                            ALLOCATE(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
                              & solverDOFToVariableDOFsMaps(solverLocalDOF(rank))%ptr%variableCoefficient(numberColEquationsCols), &
                              & STAT=err)
                            IF(err/=0) CALL FlagError("Could not allocate solver DOF to variable maps variable coefficient.", &
                              & err,error,*999)
                            ALLOCATE(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr% &
                              & solverDOFToVariableDOFsMaps(solverLocalDOF(rank))%ptr% &
                              & additiveConstant(numberColEquationsCols),STAT=err)
                            IF(err/=0) CALL FlagError("Could not allocate solver DOF to variable maps additive constant.", &
                              & err,error,*999)
                            !Setup
                            solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%solverDOFToVariableDOFsMaps( &
                              & solverLocalDOF(rank))%ptr%numberOfEquationDOFs=numberColEquationsCols
                            DO colEquationsColIdx=1,numberColEquationsCols
                              eqnLocalDof=colEquationsCols%localDofs(colEquationsColIdx)
                              couplingCoefficient=colEquationsCols%coefficients(colEquationsColIdx)
                              solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%solverDOFToVariableDOFsMaps(&
                                & solverLocalDOF(rank))%ptr%equationTypes(colEquationsColIdx)= &
                                & SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION
                              solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%solverDOFToVariableDOFsMaps(&
                                & solverLocalDOF(rank))%ptr%equationIndices(colEquationsColIdx)=interfaceConditionIdx
                              solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%solverDOFToVariableDOFsMaps(&
                                & solverLocalDOF(rank))%ptr%variable(colEquationsColIdx)%ptr=>lagrangeVariable
                              solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%solverDOFToVariableDOFsMaps(&
                                & solverLocalDOF(rank))%ptr%variableDOF(colEquationsColIdx)=eqnLocalDof
                              solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%solverDOFToVariableDOFsMaps(&
                                & solverLocalDOF(rank))%ptr%variableCoefficient(colEquationsColIdx)=couplingCoefficient
                              solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%solverDOFToVariableDOFsMaps(&
                                & solverLocalDOF(rank))%ptr%additiveConstant(1)=0.0_DP
                            ENDDO !colEquationsColIdx
                          ENDIF
                          DO colEquationsColIdx=1,numberColEquationsCols
                            eqnLocalDof=colEquationsCols%localDofs(colEquationsColIdx)
                            couplingCoefficient=colEquationsCols%coefficients(colEquationsColIdx)
                            IF(equationType==SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION) THEN
                              IF(.NOT.variableRankProcessed(solverVariableIdx,rank)) THEN
                                !Set up the equations variables -> solver columns mapping
                                !No coupling yet so the mapping is 1-1
                                solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                                  & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%&
                                  & lagrangeVarDOFToSolverDOFsMap%dofNumbers(eqnLocalDof)=solverGlobalDOF
                                solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                                  & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & lagrangeVarDOFToSolverDOFsMap%couplingCoefficients(eqnLocalDof)=couplingCoefficient
                                solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                                  & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & lagrangeVarDOFToSolverDOFsMap%additiveConstants(eqnLocalDof)=0.0_DP
                                !Set up the equations columns -> solver columns mapping
                                interfaceColumn=interfaceMapping%lagrangeDOFToColumnMap(localDOF)
                                !Allocate the equation to solver map column items.
                                ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                                  & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & interfaceColToSolverColsMap(interfaceColumn)%rowCols(1),STAT=err)
                                IF(err/=0) CALL FlagError("Could not allocate interface matrices to solver matrix maps "// &
                                  & "interface column to solver columns map rowcols.",err,error,*999)
                                ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                                  & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & interfaceColToSolverColsMap(interfaceColumn)%couplingCoefficients(1),STAT=err)
                                IF(err/=0) CALL FlagError("Could not allocate interface matrices to solver matrix maps "// &
                                  & "interface column to solver columns map coupling coefficients.",err,error,*999)
                                solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                                  & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & interfaceColToSolverColsMap(interfaceColumn)%numberOfRowCols=1
                                solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                                  & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & interfaceColToSolverColsMap(interfaceColumn)%rowCols(1)=solverGlobalDOF
                                solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                                  & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & interfaceColToSolverColsMap(interfaceColumn)%couplingCoefficients(1)=couplingCoefficient
                              ENDIF
                            ELSE
                              !Set up the equations variables -> solver columns mapping
                              solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                                & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                & dependentVarDOFToSolverDOFsMaps(interfaceMatrixIdx)%ptr%dofNumbers(eqnLocalDof)= &
                                & solverGlobalDOF
                              solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                                & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                & dependentVarDOFToSolverDOFsMaps(interfaceMatrixIdx)%ptr%couplingCoefficients(eqnLocalDof)= &
                                & couplingCoefficient
                              solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                                & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                & dependentVarDOFToSolverDOFsMaps(interfaceMatrixIdx)%ptr%additiveConstants(eqnLocalDof)=0.0_DP
                              !Set up the equations columns -> solver columns mapping
                              interfaceRow=interfaceMapping%interfaceMatrixToVarMaps(interfaceMatrixIdx)%ptr% &
                                & variableDOFToRowMap(eqnLocalDof)
                              !Allocate the equation to solver map column items.
                              ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                                & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                & interfaceMatrixToSolverMatrixMaps(interfaceMatrixIdx)%ptr% &
                                & interfaceRowToSolverColsMap(interfaceRow)%rowCols(1),STAT=err)
                              IF(err/=0) CALL FlagError("Could not allocate interface matrix to solver matrix maps "// &
                                & "interface row to solver columns map rowcols.",err,error,*999)
                              ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                                & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                & interfaceMatrixToSolverMatrixMaps(interfaceMatrixIdx)%ptr% &
                                & interfaceRowToSolverColsMap(interfaceRow)%couplingCoefficients(1),STAT=err)
                              IF(err/=0) CALL FlagError("Could not allocate interface matrix to solver matrix maps "// &
                                & "interface row to solver columns map coupling coefficients.",err,error,*999)
                              solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                                & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                & interfaceMatrixToSolverMatrixMaps(interfaceMatrixIdx)%ptr% &
                                & interfaceRowToSolverColsMap(interfaceRow)%numberOfRowCols=1
                              solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                                & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                & interfaceMatrixToSolverMatrixMaps(interfaceMatrixIdx)%ptr% &
                                & interfaceRowToSolverColsMap(interfaceRow)%rowCols(1)=solverGlobalDOF
                              solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                                & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                & interfaceMatrixToSolverMatrixMaps(interfaceMatrixIdx)%ptr% &
                                & interfaceRowToSolverColsMap(interfaceRow)%couplingCoefficients(1)= &
                                & interfaceMapping%interfaceMatrixToVarMaps(interfaceMatrixIdx)%ptr%matrixCoefficient
                            ENDIF
                          ENDDO !colEquationsColIdx
                        ELSE !rank /= myrank
                          IF(.NOT.variableRankProcessed(solverVariableIdx,rank)) THEN
                              
                            !Set up the column domain mappings.
                            CALL DomainGlobalMapping_Initialise(columnDomainMapping%globalToLocalMap(solverGlobalDOF), &
                              & err,error,*999)
                            !There are no ghosted cols for the solver matrices so there is only 1 domain for the global to
                            !local map.
                            !Allocate the global to local map arrays
                            ALLOCATE(columnDomainMapping%globalToLocalMap(solverGlobalDOF)%localNumber(1),STAT=err)
                            IF(err/=0) CALL FlagError("Could not allocate column domain global to local map local number.", &
                              & err,error,*999)
                            ALLOCATE(columnDomainMapping%globalToLocalMap(solverGlobalDOF)%domainNumber(1),STAT=err)
                            IF(err/=0) CALL FlagError("Could not allocate column domain global to local map domain number.", &
                              & err,error,*999)
                            ALLOCATE(columnDomainMapping%globalToLocalMap(solverGlobalDOF)%localType(1),STAT=err)
                            IF(err/=0) CALL FlagError("Could not allocate column domain global to local map domain number.", &
                              & err,error,*999)
                            !Set up the global to local mappings
                            columnDomainMapping%globalToLocalMap(solverGlobalDOF)%numberOfDomains=1
                            columnDomainMapping%globalToLocalMap(solverGlobalDOF)%localNumber(1)=solverLocalDOF(rank)
                            columnDomainMapping%globalToLocalMap(solverGlobalDOF)%domainNumber(1)=rank
                            columnDomainMapping%globalToLocalMap(solverGlobalDOF)%localType(1)=DOMAIN_LOCAL_INTERNAL
                            
                          ENDIF
                        ENDIF !rank == myrank
                      ELSE IF(constrainedDOF) THEN
                        !Do nothing, this is set up above
                      ELSE
                        IF(rank==myrank) THEN
                          IF(equationType==SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION) THEN
                            !Set up the equations variables -> solver columns mapping
                            solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                              & lagrangeVarDOFToSolverDOFsMap%dofNumbers(localDOF)=0
                            solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                              & lagrangeVarDOFToSolverDOFsMap%couplingCoefficients(localDOF)=0.0_DP
                            solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                              & lagrangeVarDOFToSolverDOFsMap%additiveConstants(localDOF)=0.0_DP
                            interfaceColumn=interfaceMapping%lagrangeDOFToColumnMap(localDOF)
                            CALL MatrixRowColCoupling_Initialise(solverMapping% &
                              & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                              & interfaceColToSolverColsMap(interfaceColumn),err,error,*999)
                          ELSE
                            !Set up the equations variables -> solver columns mapping
                            solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                              & dependentVarDOFToSolverDOFsMaps(interfaceMatrixIdx)%ptr%dofNumbers(localDOF)=0
                            solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                              & dependentVarDOFToSolverDOFsMaps(interfaceMatrixIdx)%ptr%couplingCoefficients(localDOF)=0.0_DP
                            solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                              & dependentVarDOFToSolverDOFsMaps(interfaceMatrixIdx)%ptr%additiveConstants(localDOF)=0.0_DP
                            !Set up the equations columns -> solver columns mapping
                            interfaceRow=interfaceMapping%interfaceMatrixToVarMaps(interfaceMatrixIdx)%ptr% &
                              & variableDOFToRowMap(localDOF)
                            CALL MatrixRowColCoupling_Initialise(solverMapping% &
                              & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr%&
                              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                              & interfaceMatrixToSolverMatrixMaps(interfaceMatrixIdx)%ptr% &
                              & interfaceRowToSolverColsMap(interfaceRow),err,error,*999)
                          ENDIF
                        ENDIF !rank==myrank
                      ENDIF !includeColumn
                    ENDDO !globalDOF
                  CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE DEFAULT
                    localError="The interface condition method of "// &
                      & TRIM(NumberToVString(interfaceCondition%method,"*",err,error))//" is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE DEFAULT
                  localError="The equation type of "//TRIM(NumberToVString(equationType,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
                variableRankProcessed(solverVariableIdx,rank)=.TRUE.
              ENDIF !Number of rank columns > 0
              IF(ALLOCATED(rankGlobalColumnsList)) DEALLOCATE(rankGlobalColumnsList)
            ENDDO !equationIdx
            
          ENDDO !solverVariableIdx
        ENDDO !rank
      ENDDO !dofType
            
      !Deallocate dof map
      DO solverVariableIdx=1,solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr%variablesList%numberOfVariables
        DEALLOCATE(dofMap(solverVariableIdx)%ptr)
      ENDDO !solverVariableIdx
      DEALLOCATE(dofMap)
      !Deallocate solver local dof
      IF(ALLOCATED(solverLocalDOF)) DEALLOCATE(solverLocalDOF)
            
      CALL DomainMapping_LocalFromGlobalCalculate(columnDomainMapping,err,error,*999)
            
      IF(ALLOCATED(subMatrixInformation)) DEALLOCATE(subMatrixInformation)
      IF(ALLOCATED(subMatrixList)) DEALLOCATE(subMatrixList)
      IF(ALLOCATED(variableRankProcessed)) DEALLOCATE(variableRankProcessed)
      IF(ALLOCATED(numberOfVariableGlobalSolverDOFs)) DEALLOCATE(numberOfVariableGlobalSolverDOFs)
      IF(ALLOCATED(numberOfVariableLocalSolverDOFs)) DEALLOCATE(numberOfVariableLocalSolverDOFs)
      IF(ALLOCATED(totalNumberOfVariableLocalSolverDOFs)) DEALLOCATE(totalNumberOfVariableLocalSolverDOFs)
      DO solverVariableIdx=1,numberOfEquationsVariables+numberOfInterfaceVariables
        CALL List_Destroy(variablesList(solverVariableIdx)%ptr,err,error,*999)
      ENDDO !solverVariableIdx
    ENDDO !solverMatrixIdx
    IF(ALLOCATED(dummyDofCoupling%globalDofs)) DEALLOCATE(dummyDofCoupling%globalDofs)
    IF(ALLOCATED(dummyDofCoupling%localDofs)) DEALLOCATE(dummyDofCoupling%localDofs)
    IF(ALLOCATED(dummyDofCoupling%coefficients)) DEALLOCATE(dummyDofCoupling%coefficients)
    CALL SolverMappingDOFCouplings_Finalise(columnCouplings,err,error,*999)

    !
    ! 5. Set up the column mappings such that the solver matrix and equations matrix orderings are the same.
    !
    
          
    DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      NULLIFY(equations)
      CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
      NULLIFY(vectorEquations)
      CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
      NULLIFY(vectorMapping)
      CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
      NULLIFY(dynamicMapping)
      CALL EquationsMappingVector_DynamicMappingExists(vectorMapping,dynamicMapping,err,error,*999)
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
      IF(ASSOCIATED(dynamicMapping)) THEN
        DO equationMatrixIdx=1,dynamicMapping%numberOfDynamicMatrices
          ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr% &
            & equationsMatrixToSolverMatrixMaps(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr%numberOfSolverMatrices),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate equations matrix to solver matrix maps.",err,error,*999)
          solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr%numberOfSolverMatrices=0
        ENDDO !equationMatrixIdx
        DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
          DO equationMatrixIdx=1,solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfDynamicMatrices
            IF(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & dynamicMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr%solverMatrixNumber==solverMatrixIdx) THEN
              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr%numberOfSolverMatrices= &
                & solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr%numberOfSolverMatrices+1
              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr% &
                & equationsMatrixToSolverMatrixMaps(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr%numberOfSolverMatrices)%ptr=> &
                & solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & dynamicMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr
            ENDIF
          ENDDO !equationMatrixIdx
        ENDDO !solverMatrixIdx
      ELSE IF(ASSOCIATED(linearMapping)) THEN
        DO equationMatrixIdx=1,linearMapping%numberOfLinearMatrices
          ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr% &
            & equationsMatrixToSolverMatrixMaps(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr%numberOfSolverMatrices),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate equations matrix to solver matrix maps.",err,error,*999)
          solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr%numberOfSolverMatrices=0
        ENDDO !equationMatrixIdx
        DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
          DO equationMatrixIdx=1,solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfLinearMatrices
            IF(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & linearMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr%solverMatrixNumber==solverMatrixIdx) THEN
              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr%numberOfSolverMatrices= &
                & solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr%numberOfSolverMatrices+1
              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr% &
                & equationsMatrixToSolverMatrixMaps(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr%numberOfSolverMatrices)%ptr=> &
                & solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & linearMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr
            ENDIF
          ENDDO !equationMatrixIdx
        ENDDO !solverMatrixIdx
      ENDIF
    ENDDO !equationsSetIdx
    DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
      NULLIFY(interfaceCondition)
      CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
      SELECT CASE(interfaceCondition%method)
      CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
        NULLIFY(interfaceEquations)
        CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
        NULLIFY(interfaceMapping)
        CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
        DO interfaceMatrixIdx=1,interfaceMapping%numberOfInterfaceMatrices
          ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr% &
            & interfaceMatrixToSolverMatrixMaps(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr%numberOfSolverMatrices),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate interface matrix to solver matrix maps.",err,error,*999)
          solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr%numberOfSolverMatrices=0
        ENDDO !interfaceMatrixIdx
        DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
          DO interfaceMatrixIdx=1,solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfInterfaceMatrices
            IF(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & interfaceMatrixToSolverMatrixMaps(interfaceMatrixIdx)%ptr%solverMatrixNumber==solverMatrixIdx) THEN
              solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr%numberOfSolverMatrices= &
                & solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr%&
                & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr%numberOfSolverMatrices+1              
              solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr% &
                & interfaceMatrixToSolverMatrixMaps(solverMapping% &
                & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr%numberOfSolverMatrices)%ptr=> &
                & solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & interfaceMatrixToSolverMatrixMaps(interfaceMatrixIdx)%ptr
            ENDIF
          ENDDO !interfaceMatrixIdx
        ENDDO !solverMatrixIdx
      CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The interface condition method of "// &
          & TRIM(NumberToVString(interfaceCondition%method,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
                        
    ENDDO !interfaceConditionIdx

    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Solver mappings:",err,error,*999)
      CALL SolverMapping_NumberOfSolverMatricesGet(solverMapping,numberOfSolverMatrices,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of solver matrices = ",numberOfSolverMatrices,err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Equation sets:",err,error,*999)
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of equations sets = ",numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Equations set index : ",equationsSetIdx,err,error,*999)
        NULLIFY(region)
        CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
        CALL Region_UserNumberGet(region,regionUserNumber,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Region user number        = ",regionUserNumber,err,error,*999)
        CALL EquationsSet_UserNumberGet(equationsSet,equationsSetUserNumber,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Equations set user number = ",equationsSetUserNumber,err,error,*999)
      ENDDO !equationsSetIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Interface conditions:",err,error,*999)
      CALL SolverMapping_NumberOfInterfaceConditionsGet(solverMapping,numberOfInterfaceConditions,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of interface conditions = ",numberOfInterfaceConditions, &
        & err,error,*999)
      DO interfaceConditionIdx=1,numberOfInterfaceConditions
        NULLIFY(interfaceCondition)
        CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Interface condition index : ",interfaceConditionIdx, &
          & err,error,*999)
        NULLIFY(INTERFACE)
        CALL InterfaceCondition_InterfaceGet(interfaceCondition,INTERFACE,err,error,*999)
        CALL Interface_UserNumberGet(INTERFACE,interfaceUserNumber,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Interface user number           = ",interfaceUserNumber, &
          & err,error,*999)
        NULLIFY(parentRegion)
        CALL Interface_ParentRegionGet(INTERFACE,parentRegion,err,error,*999)
        CALL Region_UserNumberGet(parentRegion,parentRegionUserNumber,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Parent region user number       = ",parentRegionUserNumber, &
          & err,error,*999)
        CALL InterfaceCondition_UserNumberGet(interfaceCondition,interfaceConditionUserNumber,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Interface condition user number = ",interfaceConditionUserNumber, &
          & err,error,*999)
      ENDDO !interfaceConditionIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Equations variables list:",err,error,*999)
      DO solverMatrixIdx=1,numberOfSolverMatrices
        NULLIFY(solverMatrixToEquationsMap)
        CALL SolverMapping_SolverMatrixToEquationsMapGet(solverMapping,solverMatrixIdx,solverMatrixToEquationsMap,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Solver matrix : ",solverMatrixIdx,err,error,*999)
        NULLIFY(solverVariablesList)
        CALL SolverMappingSMToEQSMap_VariablesListGet(solverMatrixToEquationsMap,solverVariablesList,err,error,*999)
        CALL SolverMappingVariables_NumberOfVariablesGet(solverVariablesList,numberOfVariables,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of variables = ",numberOfVariables,err,error,*999)
        DO variableIdx=1,numberOfVariables
          NULLIFY(solverMappingVariable)
          CALL SolverMappingVariables_VariableGet(solverVariablesList,variableIdx,solverMappingVariable,err,error,*999)
          NULLIFY(fieldVariable)
          CALL SolverMappingVariable_FieldVariableGet(solverMappingVariable,fieldVariable,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Variable : ",variableIdx,err,error,*999)
          CALL FieldVariable_VariableTypeGet(fieldVariable,variableType,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Variable type = ",variableType,err,error,*999)
          CALL SolverMappingVariable_NumberOfEquationsGet(solverMappingVariable,numberOfEquations,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of equations = ",numberOfEquations,err,error,*999)          
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfEquations,5,5,solverVariablesList%variables(variableIdx)%ptr% &
            & equationTypes,'("        Equation types   :",5(X,I13))','(26X,5(X,I13))',err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfEquations,5,5,solverVariablesList%variables(variableIdx)%ptr% &
            & equationIndices,'("        Equation indices :",5(X,I13))','(26X,5(X,I13))',err,error,*999)
        ENDDO !variableIdx
      ENDDO !solverMatrixIdx
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    RHS vector : ",solverMatrixIdx,err,error,*999)
      NULLIFY(rhsVariablesList)
      CALL SolverMapping_RHSVariablesListGet(solverMapping,rhsVariablesList,err,error,*999)
      CALL SolverMappingVariables_NumberOfVariablesGet(rhsVariablesList,numberOfVariables,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of variables = ",numberOfVariables,err,error,*999)
      DO variableIdx=1,numberOfVariables
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Variable : ",variableIdx,err,error,*999)
        NULLIFY(solverMappingVariable)
        CALL SolverMappingVariables_VariableGet(rhsVariablesList,variableIdx,solverMappingVariable,err,error,*999)
        NULLIFY(fieldVariable)
        CALL SolverMappingVariable_FieldVariableGet(solverMappingVariable,fieldVariable,err,error,*999)
        CALL FieldVariable_VariableTypeGet(fieldVariable,variableType,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Variable type = ",variableType,err,error,*999)
        CALL SolverMappingVariable_NumberOfEquationsGet(solverMappingVariable,numberOfEquations,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of equations = ",numberOfEquations,err,error,*999)        
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfEquations,5,5,rhsVariablesList%variables(variableIdx)%ptr% &
          & equationTypes,'("        Equation types   :",5(X,I13))','(26X,5(X,I13))',err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfEquations,5,5,rhsVariablesList%variables(variableIdx)%ptr% &
          & equationIndices,'("        Equation indices :",5(X,I13))','(26X,5(X,I13))',err,error,*999)
      ENDDO !variableIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Solver row to equations rows mappings:",err,error,*999)      
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of rows = ",solverMapping%numberOfRows,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of global rows = ",solverMapping%numberOfGlobalRows, &
        & err,error,*999)
      IF(diagnostics2) THEN
        CALL SolverMapping_NumberOfRowsGet(solverMapping,numberOfRows,err,error,*999)
        DO rowIdx=1,numberOfRows
          NULLIFY(solverRowToEquationsMap)
          CALL SolverMapping_SolverRowToEquationsMapGet(solverMapping,rowIdx,solverRowToEquationsMap,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Solver row : ",rowIdx,err,error,*999)
          CALL SolverMappingSRowToEQSMap_NumberOfEquationsSetRowsGet(solverRowToEquationsMap,numberOfEquationsSetRows, &
            & err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of equations set rows mapped to = ", &
            & numberOfEquationsSetRows,err,error,*999)
          CALL SolverMappingSRowToEQSMap_InterfaceConditionIndexGet(solverRowToEquationsMap,interfaceConditionIndex, &
            & err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Interface condition index          = ",interfaceConditionIndex, &
            & err,error,*999)
          IF(interfaceConditionIndex==0) THEN
            !Row is an equations set row            
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfEquationsSetRows,5,5, &
              & solverRowToEquationsMap%equationsIndex,'("      Equations indices      :",5(X,I13))','(30X,5(X,I13))', &
              & err,error,*999) 
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfEquationsSetRows,5,5, &
              & solverRowToEquationsMap%rowColNumber,'("      Equations row numbers  :",5(X,I13))','(30X,5(X,I13))', &
              & err,error,*999) 
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfEquationsSetRows,5,5, &
              & solverRowToEquationsMap%couplingCoefficients,'("      Coupling coefficients  :",5(X,E13.6))', &
              & '(30X,5(X,E13.6))',err,error,*999)
          ELSE
            !Row is an interface condition row
!!TODO: format better
            CALL WriteStringFmtValue(DIAGNOSTIC_OUTPUT_TYPE,"        Interface col numbers : ", &
              & solverRowToEquationsMap%rowColNumber(1),"(I13)",err,error,*999) 
            CALL WriteStringFmtValue(DIAGNOSTIC_OUTPUT_TYPE,"        Coupling coefficients : ", &
              & solverRowToEquationsMap%couplingCoefficients(1),"(E13.6)",err,error,*999)
          ENDIF
        ENDDO !rowIdx
      ENDIF
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Solver column to equations column mappings:",err,error,*999)      
      DO solverMatrixIdx=1,numberOfSolverMatrices
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Solver matrix : ",solverMatrixIdx,err,error,*999)
        NULLIFY(solverMatrixToEquationsMap)
        CALL SolverMapping_SolverMatrixToEquationsMapGet(solverMapping,solverMatrixIdx,solverMatrixToEquationsMap,err,error,*999)
        CALL SolverMappingSMToEQSMap_NumberOfColumnsGet(solverMatrixToEquationsMap,numberOfColumns,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of columns = ",numberOfColumns,err,error,*999)
        IF(diagnostics2) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      Solver column to equations set columns mappings:",err,error,*999)
          DO equationsSetIdx=1,numberOfEquationsSets
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Equations set index : ",equationsSetIdx,err,error,*999)
            NULLIFY(solverMatrixToEquationsSetMap)
            CALL SolverMappingSMToEQSMap_SolverMatrixToEquationsSetMapGet(solverMatrixToEquationsMap,equationsSetIdx, &
              & solverMatrixToEquationsSetMap,err,error,*999)
            CALL SolverMappingSMToESMap_HaveEquationsGet(solverMatrixToEquationsSetMap,haveDynamic,haveLinear,haveNonlinear, &
              & err,error,*999)
            DO columnIdx=1,numberOfColumns           
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE, "          Solver column : ",columnIdx,err,error,*999)
              IF(haveDynamic) THEN
                NULLIFY(solverColToDynamicEquationsMap)
                CALL SolverMappingSMToESMap_SolverColToDynamicEquationsMapGet(solverMatrixToEquationsSetMap,columnIdx, &
                  & solverColToDynamicEquationsMap,err,error,*999)
                CALL SolverMappingSColToDEQSMap_NumberOfDynamicMatricesGet(solverColToDynamicEquationsMap, &
                  & numberOfDynamicMatrices,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of dynamic equations matrices mapped to = ", &
                  & numberOfDynamicMatrices,err,error,*999)
                IF(numberOfDynamicMatrices>0) THEN
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfDynamicMatrices,5,5, &
                    & solverColToDynamicEquationsMap%equationsMatrixNumbers, &
                    & '("            Equation matrices numbers :",5(X,I13))','(39X,5(X,I13))',err,error,*999)                   
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfDynamicMatrices,5,5, &
                    & solverColToDynamicEquationsMap%equationsColumnNumbers, &
                    & '("            Equation column numbers   :",5(X,I13))','(39X,5(X,I13))',err,error,*999)
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfDynamicMatrices,5,5, &
                    & solverColToDynamicEquationsMap%couplingCoefficients, &
                    & '("            Coupling coefficients     :",5(X,E13.6))','(39X,5(X,E13.6))',err,error,*999)
                ENDIF
              ELSE
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of dynamic equations matrices mapped to = ", &
                  & 0_INTG,err,error,*999)
              ENDIF
              IF(haveLinear) THEN
                NULLIFY(solverColToLinearEquationsMap)
                CALL SolverMappingSMToESMap_SolverColToLinearEquationsMapGet(solverMatrixToEquationsSetMap,columnIdx, &
                  & solverColToLinearEquationsMap,err,error,*999)
                CALL SolverMappingSColToLEQSMap_NumberOfLinearMatricesGet(solverColToLinearEquationsMap, &
                  & numberOfLinearMatrices,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE, "            Number of linear equations matrices mapped to  = ", &
                  & numberOfLinearMatrices,err,error,*999)
                IF(numberOfLinearMatrices>0) THEN
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfLinearMatrices,5,5, &
                    & solverColToLinearEquationsMap%equationsMatrixNumbers, &
                    & '("            Equation matrices numbers :",5(X,I13))','(36X,5(X,I13))',err,error,*999)
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfLinearMatrices,5,5, &
                    & solverColToLinearEquationsMap%equationsColumnNumbers, &
                    & '("            Equation column numbers   :",5(X,I13))','(36X,5(X,I13))',err,error,*999)
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfLinearMatrices,5,5, &
                    & solverColToLinearEquationsMap%couplingCoefficients, &
                    & '("            Coupling coefficients     :",5(X,E13.6))','(36X,5(X,E13.6))',err,error,*999)
                ENDIF
              ELSE
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of linear equations matrices mapped to  = ", &
                  & 0_INTG,err,error,*999)
              ENDIF
              IF(haveNonlinear) THEN
                NULLIFY(solverColToNonlinearEquationsMap)
                CALL SolverMappingSMToESMap_SolverColToNonlinearEquationsMapGet(solverMatrixToEquationsSetMap,columnIdx, &
                  & solverColToNonlinearEquationsMap,err,error,*999)
                CALL SolverMappingSColToNLEQSMap_NumberOfJacobianMatricesGet(solverColToNonlinearEquationsMap, &
                  & numberOfJacobianMatrices,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE, "            Number of Jacobian equations matrices mapped to  = ", &
                  & numberOfJacobianMatrices,err,error,*999)
                IF(numberOfJacobianMatrices>0) THEN
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfJacobianMatrices,5,5, &
                    & solverColToNonlinearEquationsMap%jacobianMatrixNumbers, &
                    & '("            Jacobian matrices numbers :",5(X,I13))','(36X,5(X,I13))',err,error,*999)
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfJacobianMatrices,5,5, &
                    & solverColToNonlinearEquationsMap%jacobianColumnNumbers, &
                    & '("            Jacobian column numbers   :",5(X,I13))','(36X,5(X,I13))',err,error,*999)
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfJacobianMatrices,5,5, &
                    & solverColToNonlinearEquationsMap%couplingCoefficients, &
                    & '("            Coupling coefficients     :",5(X,E13.6))','(36X,5(X,E13.6))',err,error,*999)
                ENDIF
              ELSE
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of Jacobian equations matrices mapped to  = ", &
                  & 0_INTG,err,error,*999)
              ENDIF
            ENDDO !columnIdx
          ENDDO !equationsSetIdx
        ENDIF
      ENDDO !solverMatrixIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Solver DOF to field DOFs mappings:",err,error,*999)
      DO solverMatrixIdx=1,numberOfSolverMatrices
        NULLIFY(solverMatrixToEquationsMap)
        CALL SolverMapping_SolverMatrixToEquationsMapGet(solverMapping,solverMatrixIdx,solverMatrixToEquationsMap,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Solver matrix : ",solverMatrixIdx,err,error,*999)
        CALL SolverMappingSMToEQSMap_NumberOfDOFsGet(solverMatrixToEquationsMap,numberOfDOFs,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of DOFs = ",numberOfDofs,err,error,*999)
        CALL SolverMappingSMToEQSMap_TotalNumberOfDOFsGet(solverMatrixToEquationsMap,totalNumberOfDOFs,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Total number of DOFs = ",totalNumberOfDofs,err,error,*999)
        CALL SolverMappingSMToEQSMap_NumberOfGlobalDOFsGet(solverMatrixToEquationsMap,numberOfGlobalDOFs,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of global DOFs = ",numberOfGlobalDofs,err,error,*999)
        ALLOCATE(variableTypes(solverMapping%numberOfEquationsSets+solverMapping%numberOfInterfaceConditions),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate variable types.",err,error,*999)
        DO dofIdx=1,totalNumberOfDofs
          NULLIFY(solverDOFToVariableDOFsMap)
          CALL SolverMappingSMToEQSMap_SolverDOFToVariableDOFsMapGet(solverMatrixToEquationsMap,dofIdx, &
            & solverDOFToVariableDOFsMap,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Solver local DOF : ",dofIdx,err,error,*999)
          CALL SolverMappingSDOFToVDOFsMap_NumberOfEquationDOFsGet(solverDOFToVariableDOFsMap,numberOfEquationDOFs, &
            & err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of equations DOFs mapped to     = ",numberOfEquationDOFs, &
            & err,error,*999)
          IF(numberOfEquationDOFs>0) THEN
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfEquationDOFs,5,5, &
              & solverDOFToVariableDOFsMap%equationTypes,'("        Equations types       :",5(X,I13))','(31X,5(X,I13))', &
              & err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfEquationDOFs,5,5, &
              & solverDOFToVariableDOFsMap%equationIndices,'("        Equations indices     :",5(X,I13))','(31X,5(X,I13))', &
              & err,error,*999)
            DO variableIdx=1,numberOfEquationDOFs
              NULLIFY(fieldVariable)
              CALL SolverMappingSDOFToVDOFsMap_VariableGet(solverDOFToVariableDOFsMap,variableIdx,fieldVariable,err,error,*999)
              CALL FieldVariable_VariableTypeGet(fieldVariable,variableTypes(variableIdx),err,error,*999)
            ENDDO !variableIdx
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfEquationDOFs,5,5,variableTypes, &
              & '("        Variable types        :",5(X,I13))','(31X,5(X,I13))',err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfEquationDOFs,5,5, &
              & solverDOFToVariableDOFsMap%variableDOF,'("        Variable DOFs         :",5(X,I13))','(31X,5(X,I13))', &
              & err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfEquationDOFs,5,5, &
              & solverDOFToVariableDOFsMap%variableCoefficient,'("        Variable coefficients :",5(X,E13.6))', &
              & '(31X,5(X,E13.6))',err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfEquationDOFs,5,5, &
              & solverDOFToVariableDOFsMap%additiveConstant,'("        Additive constants    :",5(X,E13.6))', &
              & '(31X,5(X,E13.6))',err,error,*999)
          ENDIF
        ENDDO !dofIdx
        IF(ALLOCATED(variableTypes)) DEALLOCATE(variableTypes)
      ENDDO !solverMatrixIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Equation sets to solver mappings:",err,error,*999)
      
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(equations)
        CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        NULLIFY(vectorMapping)
        CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
        NULLIFY(lhsMapping)
        CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
        CALL EquationsMappingLHS_TotalNumberOfRowsGet(lhsMapping,totalNumberOfRows,err,error,*999)
        NULLIFY(dynamicMapping)
        CALL EquationsMappingVector_DynamicMappingExists(vectorMapping,dynamicMapping,err,error,*999)
        NULLIFY(linearMapping)
        CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
        NULLIFY(nonlinearMapping)
        CALL EquationsMappingVector_NonlinearMappingExists(vectorMapping,nonlinearMapping,err,error,*999)
        NULLIFY(rhsMapping)
        CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
        NULLIFY(sourceMapping)
        CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Equations set index : ",equationsSetIdx,err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      Equations sets rows to solver rows mappings:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE, "        Number of equations set rows = ",totalNumberOfRows,err,error,*999)
        NULLIFY(equationsSetToSolverMatricesMap)
        CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
          & err,error,*999)
        NULLIFY(equationsRowToSolverRowsMap)
        CALL SolverMappingESToSMSMap_EquationsRowToSolverRowsMapGet(equationsSetToSolverMatricesMap,equationsRowToSolverRowsMap, &
          & err,error,*999)
        DO rowIdx=1,totalNumberOfRows
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Equations set row : ",rowIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Number of solver rows mapped to   = ", &
            & equationsRowToSolverRowsMap(rowIdx)%numberOfRowCols,err,error,*999)
          IF(equationsRowToSolverRowsMap(rowIdx)%numberOfRowCols>0) THEN
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,equationsRowToSolverRowsMap(rowIdx)%numberOfRowCols,5,5, &
              & equationsRowToSolverRowsMap(rowIdx)%rowCols,'("          Solver row numbers    :",5(X,I13))','(33X,5(X,I13))', &
              & err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,equationsRowToSolverRowsMap(rowIdx)%numberOfRowCols,5,5, &
              & equationsRowToSolverRowsMap(rowIdx)%couplingCoefficients, &
              & '("          Coupling coefficients :",5(X,E13.6))','(33X,5(X,E13.6))',err,error,*999)
          ENDIF
        ENDDO !rowIdx
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      Solver matrix indexing:",err,error,*999)
        DO solverMatrixIdx=1,numberOfSolverMatrices
          NULLIFY(equationsMatricesToSolverMatrixMap)
          CALL SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet(equationsSetToSolverMatricesMap,solverMatrixIdx, &
            & equationsMatricesToSolverMatrixMap,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Solver matrix : ",solverMatrixIdx,err,error,*999)
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"          Interface conditions affecting:",err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of interface conditions = ", &
            & equationsSetToSolverMatricesMap%numberOfInterfaceConditions,err,error,*999)
          DO interfaceConditionIdx=1,equationsSetToSolverMatricesMap%numberOfInterfaceConditions
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Interface condition : ",interfaceConditionIdx,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Interface condition index = ", &
              & equationsSetToSolverMatricesMap%interfaceConditionToEquationsSetMaps(interfaceConditionIdx)% &
              & interfaceConditionIndex,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Interface matrix number   = ", &
              & equationsSetToSolverMatricesMap%interfaceConditionToEquationsSetMaps(interfaceConditionIdx)% &
              & interfaceMatrixNumber,err,error,*999)
          ENDDO !interfaceConditionIdx                    
          IF(ASSOCIATED(dynamicMapping)) THEN
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"          Dynamic equations matrix columns to solver matrix columns:", &
              & err,error,*999)
            CALL SolverMappingEMSToSMMap_NumberOfDynamicMatricesGet(equationsMatricesToSolverMatrixMap,numberOfDynamicMatrices, &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of dynamic equations matrices = ", &
              & numberOfDynamicMatrices,err,error,*999)
            DO equationMatrixIdx=1,numberOfDynamicMatrices
              NULLIFY(equationsMatrixToSolverMatrixMap)
              CALL SolverMappingEMSToSMMap_DynamicMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
                & equationMatrixIdx,equationsMatrixToSolverMatrixMap,err,error,*999)
              NULLIFY(equationsMatrixToVariableMap)
              CALL EquationsMappingDynamic_EquationsMatrixToVarMapGet(dynamicMapping,equationMatrixIdx, &
                & equationsMatrixToVariableMap,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Equations matrix index : ",equationMatrixIdx, &
                & err,error,*999)
              equationsMatrixNumber=equationsMatrixToSolverMatrixMap%equationsMatrixNumber
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Equations matrix number = ",equationsMatrixNumber, &
                & err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Solver matrix number    = ", &
                & equationsMatrixToSolverMatrixMap%solverMatrixNumber,err,error,*999)
              DO columnIdx=1,equationsMatrixToVariableMap%numberOfColumns
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Equations matrix column : ",columnIdx, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                  Number of solver columns mapped to = ", &
                  & equationsMatrixToSolverMatrixMap%equationsColToSolverColsMap(columnIdx)%numberOfRowCols,err,error,*999)
                IF(equationsMatrixToSolverMatrixMap%equationsColToSolverColsMap(columnIdx)%numberOfRowCols>0) THEN
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,equationsMatrixToSolverMatrixMap%equationsColToSolverColsMap( &
                    & columnIdx)%numberOfRowCols,5,5,equationsMatrixToSolverMatrixMap%equationsColToSolverColsMap(columnIdx)% &
                    & rowCols,'("                  Solver columns         :",5(X,I13))','(42X,5(X,I13))',err,error,*999)
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,equationsMatrixToSolverMatrixMap%equationsColToSolverColsMap( &
                    & columnIdx)%numberOfRowCols,5,5,equationsMatrixToSolverMatrixMap%equationsColToSolverColsMap(columnIdx)% &
                    & couplingCoefficients,'("                  Coupling coefficients  :",5(X,E13.6))','(42X,5(X,E13.6))', &
                    & err,error,*999)
                ENDIF
              ENDDO !columnIdx
            ENDDO !equationMatrixIdx
          ELSE
            IF(ASSOCIATED(linearMapping)) THEN
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"          Linear equations matrix columns to solver matrix columns:", &
                & err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of liner equations matrices = ", &
                & equationsMatricesToSolverMatrixMap%numberOfLinearMatrices,err,error,*999)
              DO equationMatrixIdx=1,equationsMatricesToSolverMatrixMap%numberOfLinearMatrices
                NULLIFY(equationsMatrixToVariableMap)
                CALL EquationsMappingLinear_EquationsMatrixToVarMapGet(linearMapping,equationMatrixIdx, &
                  & equationsMatrixToVariableMap,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Equations matrix index : ",equationMatrixIdx, &
                  & err,error,*999)
                equationsMatrixNumber=equationsMatrixToSolverMatrixMap%equationsMatrixNumber
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Equations matrix number = ",equationsMatrixNumber, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Solver matrix number    = ", &
                  & equationsMatrixToSolverMatrixMap%solverMatrixNumber,err,error,*999)
                DO columnIdx=1,equationsMatrixToVariableMap%numberOfColumns
                  CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Equations matrix column : ",columnIdx, &
                    & err,error,*999)
                  CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Number of solver columns mapped to = ", &
                    & equationsMatrixToSolverMatrixMap%equationsColToSolverColsMap(columnIdx)%numberOfRowCols,err,error,*999)
                  IF(equationsMatrixToSolverMatrixMap%equationsColToSolverColsMap(columnIdx)%numberOfRowCols>0) THEN
                    CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,equationsMatrixToSolverMatrixMap% &
                      & equationsColToSolverColsMap(columnIdx)%numberOfRowCols,5,5,equationsMatrixToSolverMatrixMap% &
                      & equationsColToSolverColsMap(columnIdx)%rowCols, &
                      & '("                Solver columns         :",5(X,I13))','(40X,5(X,I13))',err,error,*999)
                    CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,equationsMatrixToSolverMatrixMap% &
                      & equationsColToSolverColsMap(columnIdx)%numberOfRowCols,5,5,equationsMatrixToSolverMatrixMap% &
                      & equationsColToSolverColsMap(columnIdx)%couplingCoefficients, &
                      & '("                Coupling coefficients  :",5(X,E13.6))','(40X,5(X,E13.6))',err,error,*999)
                  ENDIF
                ENDDO !columnIdx
              ENDDO !equationMatrixIdx
            ENDIF
            IF(ASSOCIATED(nonlinearMapping)) THEN
              CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"          Jacobian equations matrix columns to solver matrix columns:", &
                & err,error,*999)
              DO residualIdx=1,nonlinearMapping%numberOfResiduals
                NULLIFY(residualMapping)
                CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Residual number        : ",residualIdx,err,error,*999)
                DO equationMatrixIdx=1,residualMapping%numberOfVariables
                  NULLIFY(jacobianMatrixToSolverMatrixMap)
                  CALL SolverMappingEMSToSMMap_JacobianMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
                    & equationMatrixIdx,jacobianMatrixToSolverMatrixMap,err,error,*999)
                  NULLIFY(jacobianMatrixToVariableMap)
                  CALL EquationsMappingResidual_JacobianMatrixToVarMapGet(residualMapping,equationMatrixIdx, &
                    & jacobianMatrixToVariableMap,err,error,*999)
                  CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Jacobian matrix number  = ", &
                    & jacobianMatrixToSolverMatrixMap%jacobianMatrixNumber,err,error,*999)
                  CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Solver matrix number    = ", &
                    & jacobianMatrixToSolverMatrixMap%solverMatrixNumber,err,error,*999)
                  DO columnIdx=1,jacobianMatrixToVariableMap%numberOfColumns
                    CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Equations matrix column : ",columnIdx,err,error,*999)
                    CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Number of solver columns mapped to = ", &
                      & jacobianMatrixToSolverMatrixMap%jacobianColToSolverColsMap(columnIdx)%numberOfRowCols,err,error,*999)
                    IF(jacobianMatrixToSolverMatrixMap%jacobianColToSolverColsMap(columnIdx)%numberOfRowCols>0) THEN
                      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,jacobianMatrixToSolverMatrixMap% &
                        & jacobianColToSolverColsMap(columnIdx)%numberOfRowCols,5,5,jacobianMatrixToSolverMatrixMap% &
                        & jacobianColToSolverColsMap(columnIdx)%rowCols,'("              Solver columns         :",5(X,I13))', &
                        & '(38X,5(X,I13))',err,error,*999)
                      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,jacobianMatrixToSolverMatrixMap% &
                        & jacobianColToSolverColsMap(columnIdx)%numberOfRowCols,5,5,jacobianMatrixToSolverMatrixMap% &
                        & jacobianColToSolverColsMap(columnIdx)%couplingCoefficients, &
                        & '("              Coupling coefficients  :",5(X,E13.6))','(38X,5(X,E13.6))',err,error,*999)
                    ENDIF
                  ENDDO !columnIdx
                ENDDO !equationMatrixIdx
              ENDDO !residualIdx
            ENDIF
          ENDIF
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"          Variable DOFs to solver matrix DOFs:",err,error,*999)
          CALL SolverMappingEMSToSMMap_NumberOfVariablesGet(equationsMatricesToSolverMatrixMap,numberOfVariables,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of variables = ",numberOfVariables,err,error,*999) 
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfVariables,5,5,equationsMatricesToSolverMatrixMap% &
            & variableTypes,'("            Variable types :",5(X,I13))','(28X,5(X,I13))',err,error,*999)
          DO variableIdx=1,numberOfVariables
            NULLIFY(dependentVariable)
            CALL SolverMappingEMSToSMMap_VariableGet(equationsMatricesToSolverMatrixMap,variableIdx,dependentVariable, &
              & err,error,*999)            
            CALL FieldVariable_NumberOfDOFsGet(dependentVariable,numberOfDOFs,err,error,*999)
            CALL FieldVariable_TotalNumberOfDOFsGet(dependentVariable,totalNumberOfDOFs,err,error,*999)
            NULLIFY(varDOFToSolverDOFsMap)
            CALL SolverMappingEMSToSMMap_VariableDOFToSolverDOFsMapGet(equationsMatricesToSolverMatrixMap,variableIdx, &
              & varDOFToSolverDOFsMap,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Variable index : ",variableIdx,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Number of variable DOFs = ",totalNumberOfDofs, &
              & err,error,*999)
            DO localDOF=1,totalNumberOfDofs
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Variable DOF : ",localDOF,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Solver DOF number    = ",varDOFToSolverDOFsMap% &
                & dofNumbers(localDOF),err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Coupling coefficient = ",varDOFToSolverDOFsMap% &
                & couplingCoefficients(localDOF),err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Additive constant    = ",varDOFToSolverDOFsMap% &
                & additiveConstants(localDOF),err,error,*999) 
            ENDDO !localDOF
          ENDDO !variableIdx
        ENDDO !solverMatrixIdx
        IF(ASSOCIATED(dynamicMapping)) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      Dynamic equations matrix indexing:",err,error,*999)
          CALL EquationsMappingDynamic_NumberOfDynamicMatricesGet(dynamicMapping,numberOfDynamicMatrices,err,error,*999)
          DO equationMatrixIdx=1,numberOfDynamicMatrices
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Equations matrix : ",equationMatrixIdx,err,error,*999)
            NULLIFY(equationsMatrixToSolverMatricesMap)
            CALL SolverMappingESToSMSMap_EquationsMatrixToSolverMatricesMapGet(equationsSetToSolverMatricesMap,equationMatrixIdx, &
              & equationsMatrixToSolverMatricesMap,err,error,*999)
            CALL SolverMappingEMToSMSMap_NumberOfSolverMatricesGet(equationsMatrixToSolverMatricesMap,numberOfSolverMatrices, &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Number of solver matrices = ",numberOfSolverMatrices, &
              & err,error,*999)
            DO solverMatrixIdx=1,numberOfSolverMatrices
              NULLIFY(equationsMatrixToSolverMatrixMap)
              CALL SolverMappingEMToSMSMap_EquationsMatrixToSolverMatrixMapGet(equationsMatrixToSolverMatricesMap, &
                & solverMatrixIdx,equationsMatrixToSolverMatrixMap,err,error,*999)
              CALL SolverMappingEMToSMMap_EquationsMatrixNumberGet(equationsMatrixToSolverMatrixMap,equationsMatrixNumber, &
                & err,error,*999)
              CALL SolverMappingEMToSMMap_SolverMatrixNumberGet(equationsMatrixToSolverMatrixMap,solverMatrixNumber, &
                & err,error,*999)              
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Solver matrix index : ",solverMatrixIdx,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Equations matrix number = ",equationsMatrixNumber, &
                & err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Solver matrix number    = ",solverMatrixNumber, &
                & err,error,*999)
            ENDDO !solverMatrixIdx
          ENDDO !equationMatrixIdx
        ELSE
          IF(ASSOCIATED(linearMapping)) THEN
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      Linear equations matrix indexing:",err,error,*999)
            CALL EquationsMappingLinear_NumberOfLinearMatricesGet(linearMapping,numberOfLinearMatrices,err,error,*999)
            DO equationMatrixIdx=1,numberOfLinearMatrices
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Equations matrix : ",equationMatrixIdx,err,error,*999)
              NULLIFY(equationsMatrixToSolverMatricesMap)
              CALL SolverMappingESToSMSMap_EquationsMatrixToSolverMatricesMapGet(equationsSetToSolverMatricesMap, &
                & equationMatrixIdx,equationsMatrixToSolverMatricesMap,err,error,*999)
              CALL SolverMappingEMToSMSMap_NumberOfSolverMatricesGet(equationsMatrixToSolverMatricesMap,numberOfSolverMatrices, &
                & err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Number of solver matrices = ",numberOfSolverMatrices, &
                & err,error,*999)
              DO solverMatrixIdx=1,numberOfSolverMatrices
                NULLIFY(equationsMatrixToSolverMatrixMap)
                CALL SolverMappingEMToSMSMap_EquationsMatrixToSolverMatrixMapGet(equationsMatrixToSolverMatricesMap, &
                  & solverMatrixIdx,equationsMatrixToSolverMatrixMap,err,error,*999)
                CALL SolverMappingEMToSMMap_EquationsMatrixNumberGet(equationsMatrixToSolverMatrixMap,equationsMatrixNumber, &
                  & err,error,*999)
                CALL SolverMappingEMToSMMap_SolverMatrixNumberGet(equationsMatrixToSolverMatrixMap,SolverMatrixNumber, &
                  & err,error,*999)              
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Solver matrix index : ",solverMatrixIdx,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Equations matrix number = ",equationsMatrixNumber, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Solver matrix number    = ",solverMatrixNumber, &
                  & err,error,*999)
              ENDDO !solverMatrixIdx
            ENDDO !equationMatrixIdx
          ENDIF
          IF(ASSOCIATED(nonlinearMapping)) THEN
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      Jacobian matrix indexing:",err,error,*999)
            CALL EquationsMappingNonlinear_NumberOfResidualsGet(nonlinearMapping,numberOfResiduals,err,error,*999)
            DO residualIdx=1,numberOfResiduals
              NULLIFY(residualMapping)
              CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
              CALL EquationsMappingResidual_NumberOfResidualVariablesGet(residualMapping,numberOfVariables,err,error,*999)
              DO equationMatrixIdx=1,numberOfVariables
                NULLIFY(jacobianMatrixToSolverMatrixMap)
                CALL SolverMappingESToSMSMap_JacobianMatrixToSolverMatrixMapGet(equationsSetToSolverMatricesMap, &
                  & equationMatrixIdx,jacobianMatrixToSolverMatrixMap,err,error,*999)
                CALL SolverMappingJMToSMMap_JacobianMatrixNumberGet(jacobianMatrixToSolverMatrixMap,jacobianMatrixNumber, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Jacobian matrix number  = ",jacobianMatrixNumber, &
                  & err,error,*999)
                CALL SolverMappingJMToSMMap_SolverMatrixNumberGet(jacobianMatrixToSolverMatrixMap,solverMatrixNumber, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Solver matrix number    = ",solverMatrixNumber, &
                  & err,error,*999)
              ENDDO !equationsMatrixIdx
            ENDDO !residualIdx
          ENDIF
        ENDIF
      ENDDO !equationsSetIdx
      IF(solverMapping%numberOfInterfaceConditions>0) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Interface conditions to solver mappings:",err,error,*999)
        DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
          NULLIFY(interfaceCondition)
          CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
          NULLIFY(interfaceEquations)
          CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
          NULLIFY(interfaceMapping)
          CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
          NULLIFY(interfaceConditionToSolverMatricesMap)
          CALL SolverMapping_InterfaceConditionToSolverMatricesMapGet(solverMapping,interfaceConditionIdx, &
            & interfaceConditionToSolverMatricesMap,err,error,*999)
          CALL SolverMappingICToSMSMap_NumberOfEquationsSetsGet(interfaceConditionToSolverMatricesMap,numberOfEquationsSets, &
            & err,error,*999)
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Interface to equations sets mapping:",err,error,*999)
          DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Solver matrix : ",solverMatrixIdx,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of equations sets = ",numberOfEquationsSets, &
              & err,error,*999)
            DO equationsSetIdx=1,numberOfEquationsSets           
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Equations set : ",equationsSetIdx,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Equations set index     = ", &
                & interfaceConditionToSolverMatricesMap%equationsSetToInterfaceConditionMaps(equationsSetIdx)%equationsSetIndex, &
                & err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Interface matrix number = ", &
                & interfaceConditionToSolverMatricesMap%equationsSetToInterfaceConditionMaps(equationsSetIdx)% &
                & interfaceMatrixIndex,err,error,*999)
            ENDDO !equationsSetIdx
          ENDDO !solverMatrixIdx
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Interface condition rows to solver rows mappings:",err,error,*999)
          CALL InterfaceMapping_NumberOfInterfaceMatricesGet(interfaceMapping,numberOfInterfaceMatrices,err,error,*999)
          DO interfaceMatrixIdx=1,interfaceMapping%numberOfInterfaceMatrices
            NULLIFY(interfaceMatrixToVarMap)
            CALL InterfaceMapping_InterfaceMatrixToVarMapGet(interfaceMapping,interfaceMatrixIdx,interfaceMatrixToVarMap, &
              & err,error,*999)
            NULLIFY(interfaceMatrixToSolverMatricesMap)
            CALL SolverMappingICToSMSMap_InterfaceMatrixToSolverMatricesMapGet(interfaceConditionToSolverMatricesMap, &
              & interfaceMatrixIdx,interfaceMatrixToSolverMatricesMap,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Interface matrix idx : ",interfaceMatrixIdx,err,error,*999)
            CALL InterfaceMappingIMToVMap_TotalNumberOfRowsGet(interfaceMatrixToVarMap,totalNumberOfRows,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of interface rows = ",totalNumberOfRows,err,error,*999)
            DO rowIdx=1,totalNumberOfRows
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Interface row : ",rowIdx,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Number of solver rows mapped to = ", &
                & interfaceMatrixToSolverMatricesMap%interfaceRowToSolverRowsMap(rowIdx)%numberOfRowCols,err,error,*999)
              IF(interfaceMatrixToSolverMatricesMap%interfaceRowToSolverRowsMap(rowIdx)%numberOfRowCols>0) THEN                
                CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,interfaceMatrixToSolverMatricesMap% &
                  & interfaceRowToSolverRowsMap(rowIdx)%numberOfRowCols,5,5,interfaceMatrixToSolverMatricesMap% &
                  & interfaceRowToSolverRowsMap(rowIdx)%rowCols, &
                  & '("          Solver row numbers    :",5(X,I13))','(33X,5(X,I13))',err,error,*999)                 
                CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,interfaceMatrixToSolverMatricesMap% &
                  & interfaceRowToSolverRowsMap(rowIdx)%numberOfRowCols,5,5,interfaceMatrixToSolverMatricesMap% &
                  & interfaceRowToSolverRowsMap(rowIdx)%couplingCoefficients, &
                  & '("          Coupling coefficients :",5(X,E13.6))','(33X,5(X,E13.6))',err,error,*999)
              ENDIF
            ENDDO !rowIdx
          ENDDO !interfaceMatrixIdx
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Interface condition column to solver rows mappings:", &
            & err,error,*999)
          CALL InterfaceMapping_TotalNumberOfColumnsGet(interfaceMapping,totalNumberOfColumns,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of interface condition columns = ",totalNumberOfColumns, &
            & err,error,*999)
          NULLIFY(interfaceColToSolverRowsMap)
          CALL SolverMappingICToSMSMap_InterfaceColToSolverRowsMapGet(interfaceConditionToSolverMatricesMap, &
            & interfaceColToSolverRowsMap,err,error,*999)
          DO columnIdx=1,totalNumberOfColumns
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Interface condition column : ",columnIdx,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of rows mapped to = ", &
              & interfaceColToSolverRowsMap(columnIdx)%numberOfRowCols,err,error,*999)
            IF(interfaceColToSolverRowsMap(columnIdx)%numberOfRowCols>0) THEN
              CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,interfaceColToSolverRowsMap(columnIdx)%numberOfRowCols,5,5, &
                & interfaceColToSolverRowsMap(columnIdx)%rowCols,'("              Solver rows           :",5(X,I13))', &
                & '(38X,5(X,I13))',err,error,*999)
              CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,interfaceColToSolverRowsMap(columnIdx)%numberOfRowCols,5,5, &
                & interfaceColToSolverRowsMap(columnIdx)%couplingCoefficients, &
                & '("              Coupling coefficients  :",5(X,E13.6))','(38X,5(X,E13.6))',err,error,*999)
            ENDIF
          ENDDO !columnIdx
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Solver matrix indexing:",err,error,*999)
          DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
            NULLIFY(interfaceMatricesToSolverMatrixMap)
            CALL SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet(interfaceConditionToSolverMatricesMap, &
              & solverMatrixIdx,interfaceMatricesToSolverMatrixMap,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Solver matrix : ",solverMatrixIdx,err,error,*999)        
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"        Interface equations matrix rows to solver matrix columns:", &
              & err,error,*999)
            CALL SolverMappingIMSToSMMap_NumberOfInterfaceMatricesGet(interfaceMatricesToSolverMatrixMap, &
              & numberOfInterfaceMatrices,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Number of interface matrices = ",numberOfInterfaceMatrices, &
              & err,error,*999)
            DO interfaceMatrixIdx=1,numberOfInterfaceMatrices
              NULLIFY(interfaceMatrixToSolverMatrixMap)
              CALL SolverMappingIMSToSMMap_InterfaceMatrixToSolverMatrixMapGet(interfaceMatricesToSolverMatrixMap, &
                & interfaceMatrixIdx,interfaceMatrixToSolverMatrixMap,err,error,*999)
              NULLIFY(interfaceRowToSolverColsMap)
              CALL SolverMappingIMToSMMap_InterfaceRowToSolverColsMapGet(interfaceMatrixToSolverMatrixMap, &
                & interfaceRowToSolverColsMap,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Interface matrix index : ",interfaceMatrixIdx, &
                & err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Interface matrix number = ", &
                & interfaceMatrixToSolverMatrixMap%interfaceMatrixNumber,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Solver matrix number    = ", &
                & interfaceMatrixToSolverMatrixMap%solverMatrixNumber,err,error,*999)
              NULLIFY(interfaceMatrixToVarMap)
              CALL InterfaceMapping_InterfaceMatrixToVarMapGet(interfaceMapping,interfaceMatrixIdx, &
                & interfaceMatrixToVarMap,err,error,*999)
              CALL InterfaceMappingIMToVMap_NumberOfRowsGet(interfaceMatrixToVarMap,numberOfRows,err,error,*999)
              DO rowIdx=1,numberOfRows
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Interface matrix row : ",rowIdx,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Number of solver columns mapped to = ", &
                  & interfaceRowToSolverColsMap(rowIdx)%numberOfRowCols,err,error,*999)
                IF(interfaceRowToSolverColsMap(rowIdx)%numberOfRowCols>0) THEN
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,interfaceRowToSolverColsMap(rowIdx)%numberOfRowCols,5,5, &
                    & interfaceRowToSolverColsMap(rowIdx)%rowCols,'("              Solver columns         :",5(X,I13))', &
                    & '(38X,5(X,I13))',err,error,*999)
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,interfaceRowToSolverColsMap(rowIdx)%numberOfRowCols,5,5, &
                    & interfaceRowToSolverColsMap(rowIdx)%couplingCoefficients, &
                    & '("              Coupling coefficients  :",5(X,E13.6))','(38X,5(X,E13.6))',err,error,*999)
                ENDIF
              ENDDO !rowIdx
            ENDDO !interfaceMatrixIdx
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"        Variable DOFs to solver matrix DOFs:",err,error,*999)
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"          Lagrange variables:",err,error,*999)
            NULLIFY(lagrangeVariable)
            CALL SolverMappingIMSToSMMap_LagrangeVariableGet(interfaceMatricesToSolverMatrixMap,lagrangeVariable,err,error,*999)
            CALL FieldVariable_VariableTypeGet(lagrangeVariable,variableType,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Lagrange variable type = ",variableType,err,error,*999)
            CALL FieldVariable_TotalNumberOfDOFsGet(lagrangeVariable,totalNumberOfDOFs,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of Lagrange variable DOFs = ",totalNumberOfDofs, &
              & err,error,*999)
            NULLIFY(lagrangeVarDOFToSolverDOFsMap)
            CALL SolverMappingIMSToSMMap_LagrangeVarDOFToSolverDOFsMapGet(interfaceMatricesToSolverMatrixMap, &
              & lagrangeVarDOFToSolverDOFsMap,err,error,*999)
            DO localDOF=1,totalNumberOfDofs
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Variable DOF : ",localDOF,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Solver DOF number    = ", &
                & lagrangeVarDOFToSolverDOFsMap%dofNumbers(localDOF),err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Coupling coefficient = ", &
                & lagrangeVarDOFToSolverDOFsMap%couplingCoefficients(localDOF),err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Additive constant    = ", &
                & lagrangeVarDOFToSolverDOFsMap%additiveConstants(localDOF),err,error,*999)              
            ENDDO !localDOF
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"          Dependent variables:",err,error,*999)
            CALL SolverMappingIMSToSMMap_NumberOfDependentVariablesGet(interfaceMatricesToSolverMatrixMap, &
              & numberOfDependentVariables,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of dependent variables = ", &
              & numberOfDependentVariables,err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfDependentVariables, &
              & 5,5,interfaceMatricesToSolverMatrixMap%dependentVariableTypes, &
              & '("            Dependent variable types :",5(X,I13))','(38X,5(X,I13))',err,error,*999) 
            DO variableIdx=1,numberOfDependentVariables
              NULLIFY(dependentVariable)
              CALL SolverMappingIMSToSMMap_DependentVariableGet(interfaceMatricesToSolverMatrixMap,variableIdx, &
                & dependentVariable,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Dependent variable index : ",variableIdx,err,error,*999)
              CALL FieldVariable_TotalNumberOfDOFsGet(dependentVariable,totalNumberOfDOFs,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Number of dependent variable DOFs = ", &
                & totalNumberOfDofs,err,error,*999)
              NULLIFY(dependentVarDOFToSolverDOFsMap)
              CALL SolverMappingIMSToSMMap_DependentVarDOFToSolverDOFsMapGet(interfaceMatricesToSolverMatrixMap,variableIdx, &
                & dependentVarDOFToSolverDOFsMap,err,error,*999)
              DO localDOF=1,totalNumberOfDofs
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Variable DOF : ",localDOF,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Solver DOF number    = ", &
                  & dependentVarDOFToSolverDOFsMap%dofNumbers(localDOF),err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Coupling coefficient = ", &
                  & dependentVarDOFToSolverDOFsMap%couplingCoefficients(localDOF),err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Additive constant    = ", &
                  & dependentVarDOFToSolverDOFsMap%additiveConstants(localDOF),err,error,*999)              
              ENDDO !localDOF
            ENDDO !variableIdx
          ENDDO !solverMatrixIdx        
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Interface equations matrix indexing:",err,error,*999)
          DO interfaceMatrixIdx=1,numberOfInterfaceMatrices
            NULLIFY(interfaceMatrixToSolverMatricesMap)
            CALL SolverMappingICToSMSMap_InterfaceMatrixToSolverMatricesMapGet(interfaceConditionToSolverMatricesMap, &
              & interfaceMatrixIdx,interfaceMatrixToSolverMatricesMap,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Interface matrix : ",interfaceMatrixIdx,err,error,*999)
            CALL SolverMappingIMToSMSMap_NumberOfSolverMatricesGet(interfaceMatrixToSolverMatricesMap,numberOfSolverMatrices, &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of solver matrices = ",numberOfSolverMatrices, &
              & err,error,*999)
            DO solverMatrixIdx=1,numberOfSolverMatrices
              NULLIFY(interfaceMatrixToSolverMatrixMap)
              CALL SolverMappingIMToSMSMap_InterfaceMatrixToSolverMatrixMapGet(interfaceMatrixToSolverMatricesMap, &
                & solverMatrixIdx,interfaceMatrixToSolverMatrixMap,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Solver matrix index : ",solverMatrixIdx,err,error,*999)
              CALL SolverMappingIMToSMMap_InterfaceMatrixNumberGet(interfaceMatrixToSolverMatrixMap,interfaceMatrixNumber, &
                & err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Interface matrix number = ",interfaceMatrixNumber, &
                & err,error,*999)
              CALL SolverMappingIMToSMMap_SolverMatrixNumberGet(interfaceMatrixToSolverMatrixMap,solverMatrixNumber, &
                & err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Solver matrix number    = ",solverMatrixNumber, &
                & err,error,*999)
            ENDDO !solverMatrixIdx
          ENDDO !equationMatrixIdx
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Interface column to solver rows mapping:",err,error,*999)
          CALL InterfaceMapping_NumberOfColumnsGet(interfaceMapping,numberOfInterfaceColumns,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of columns = ",numberOfColumns,err,error,*999)
          NULLIFY(interfaceColToSolverRowsMap)
          CALL SolverMappingICToSMSMap_InterfaceColToSolverRowsMapGet(interfaceConditionToSolverMatricesMap, &
            & interfaceColToSolverRowsMap,err,error,*999)
          DO columnIdx=1,numberOfColumns
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Column : ",columnIdx, err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of solver rows = ", &
              & interfaceColToSolverRowsMap(columnIdx)%numberOfRowCols,err,error,*999)
            IF(interfaceColToSolverRowsMap(columnIdx)%numberOfRowCols>0) THEN
              CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,interfaceColToSolverRowsMap(columnIdx)%numberOfRowCols,5,5, &
                & interfaceColToSolverRowsMap(columnIdx)%rowCols,'("        Solver rows           :",5(X,I13))', &
                & '(32X,5(X,I13))',err,error,*999)
              CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,interfaceColToSolverRowsMap(columnIdx)%numberOfRowCols,5,5, &
                & interfaceColToSolverRowsMap(columnIdx)%couplingCoefficients,'("        Coupling coefficients  :",5(X,E13.6))', &
                & '(32X,5(X,E13.6))',err,error,*999)
            ENDIF
          ENDDO !columnIdx
        ENDDO !interfaceConditionIdx
      ENDIF
    ENDIF
    
    EXITS("SolverMapping_Calculate")
    RETURN
999 IF(ALLOCATED(subMatrixInformation)) DEALLOCATE(subMatrixInformation)
    IF(ALLOCATED(subMatrixList)) DEALLOCATE(subMatrixList)
    IF(ALLOCATED(equationsRHSVariables)) DEALLOCATE(equationsRHSVariables)
    IF(ALLOCATED(variableRankProcessed)) DEALLOCATE(variableRankProcessed)
    IF(ALLOCATED(numberOfVariableGlobalSolverDOFs)) DEALLOCATE(numberOfVariableGlobalSolverDOFs)
    IF(ALLOCATED(numberOfVariableLocalSolverDOFs)) DEALLOCATE(numberOfVariableLocalSolverDOFs)
    IF(ALLOCATED(totalNumberOfVariableLocalSolverDOFs)) DEALLOCATE(totalNumberOfVariableLocalSolverDOFs)    
    IF(ALLOCATED(dummyDofCoupling%globalDofs)) DEALLOCATE(dummyDofCoupling%globalDofs)
    IF(ALLOCATED(dummyDofCoupling%localDofs)) DEALLOCATE(dummyDofCoupling%localDofs)
    IF(ALLOCATED(dummyDofCoupling%coefficients)) DEALLOCATE(dummyDofCoupling%coefficients)
    CALL SolverMappingDOFCouplings_Finalise(rowCouplings,err,error,*998)
998 CALL SolverMappingDOFCouplings_Finalise(columnCouplings,err,error,*997)
997 ERRORSEXITS("SolverMapping_Calculate",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_Calculate

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a solver mapping
  SUBROUTINE SolverMapping_CreateFinish(solverMapping,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverMapping_CreateFinish",err,error,*998)

    CALL SolverMapping_AssertNotFinished(solverMapping,err,error,*998)

    IF(ASSOCIATED(solverMapping%createValuesCache)) THEN
      CALL SolverMapping_Calculate(solverMapping,err,error,*999)
      CALL SolverMapping_CreateValuesCacheFinalise(solverMapping%createValuesCache,err,error,*999)
      solverMapping%solverMappingFinished=.TRUE.
    ELSE
      CALL FlagError("Solver mapping create values cache is not associated",err,error,*999)
    ENDIF
       
    EXITS("SolverMapping_CreateFinish")
    RETURN
999 CALL SolverMapping_Finalise(solverMapping,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverMapping_CreateFinish",err,error)
    RETURN 1
  END SUBROUTINE SolverMapping_CreateFinish

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a solver mapping for a problem solver
  SUBROUTINE SolverMapping_CreateStart(solverEquations,solverMapping,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to create the solver mapping on.
    TYPE(SolverMappingType), POINTER :: solverMapping !<On return, a pointer to the solver mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverMapping_CreateStart",err,error,*998)

    IF(ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is already associated.",err,error,*998)
    CALL SolverEquations_AssertNotFinished(solverEquations,err,error,*999)
    
    CALL SolverMapping_Initialise(solverEquations,err,error,*999)
    solverMapping=>solverEquations%solverMapping
      
    EXITS("SolverMapping_CreateStart")
    RETURN
999 NULLIFY(solverMapping)
998 ERRORSEXITS("SolverMapping_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_CreateStart

  !
  !================================================================================================================================
  !

  !>Finalises a solver mapping create values cache and deallocates all memory
  SUBROUTINE SolverMapping_CreateValuesCacheFinalise(createValuesCache,err,error,*)

    !Argument variables
    TYPE(SolverMappingCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the create values cache
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,solverMatrixIdx

    ENTERS("SolverMapping_CreateValuesCacheFinalise",err,error,*999)

    IF(ASSOCIATED(createValuesCache)) THEN
      IF(ASSOCIATED(createValuesCache%equationsVariableList)) THEN
        DO solverMatrixIdx=1,SIZE(createValuesCache%equationsVariableList,1)
          IF(ASSOCIATED(createValuesCache%equationsVariableList(solverMatrixIdx)%ptr)) & 
            & CALL List_Destroy(createValuesCache%equationsVariableList(solverMatrixIdx)%ptr,err,error,*999)
        ENDDO !solverMatrixIdx
        DEALLOCATE(createValuesCache%equationsVariableList)
      ENDIF
      IF(ASSOCIATED(createValuesCache%equationsRHSVariableList)) &
        & CALL List_Destroy(createValuesCache%equationsRHSVariableList,err,error,*999)
      IF(ALLOCATED(createValuesCache%dynamicVariableType)) DEALLOCATE(createValuesCache%dynamicVariableType)
      IF(ALLOCATED(createValuesCache%matrixVariableTypes)) DEALLOCATE(createValuesCache%matrixVariableTypes)
      IF(ALLOCATED(createValuesCache%residualVariableTypes)) DEALLOCATE(createValuesCache%residualVariableTypes)
      IF(ALLOCATED(createValuesCache%rhsVariableType)) DEALLOCATE(createValuesCache%rhsVariableType)
      IF(ALLOCATED(createValuesCache%sourceVariableTypes)) DEALLOCATE(createValuesCache%sourceVariableTypes)
      IF(ASSOCIATED(createValuesCache%interfaceVariableList)) THEN
        DO solverMatrixIdx=1,SIZE(createValuesCache%interfaceVariableList,1)
          IF(ASSOCIATED(createValuesCache%interfaceVariableList(solverMatrixIdx)%ptr)) & 
            & CALL List_Destroy(createValuesCache%interfaceVariableList(solverMatrixIdx)%ptr,err,error,*999)
        ENDDO !solverMatrixIdx
        DEALLOCATE(createValuesCache%interfaceVariableList)
      ENDIF
      IF(ASSOCIATED(createValuesCache%interfaceIndices)) THEN
        DO equationsSetIdx=1,SIZE(createValuesCache%interfaceIndices,1)
          IF(ASSOCIATED(createValuesCache%interfaceIndices(equationsSetIdx)%ptr)) &
            & CALL List_Destroy(createValuesCache%interfaceIndices(equationsSetIdx)%ptr, &
            & err,error,*999)
        ENDDO !equaitons_set_idx
        DEALLOCATE(createValuesCache%interfaceIndices)
      ENDIF
      DEALLOCATE(createValuesCache)
    ENDIF
       
    EXITS("SolverMapping_CreateValuesCacheFinalise")
    RETURN
999 ERRORSEXITS("SolverMapping_CreateValuesCacheFinalise",err,error)
    RETURN 1
  END SUBROUTINE SolverMapping_CreateValuesCacheFinalise

  !
  !================================================================================================================================
  !

  !>Initialises a solver mapping create values cache 
  SUBROUTINE SolverMapping_CreateValuesCacheInitialise(solverMapping,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the create values cache
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,equationsSetIdx,solverMatrixIdx
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverMapping_CreateValuesCacheInitialise",err,error,*998)
    
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*998)
    IF(ASSOCIATED(solverMapping%createValuesCache)) &
      & CALL FlagError("Solver mapping create values cache is already associated.",err,error,*998)
      
    ALLOCATE(solverMapping%createValuesCache,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver mapping create values cache.",err,error,*999)
    ALLOCATE(solverMapping%createValuesCache%equationsVariableList(solverMapping%numberOfSolverMatrices),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver mapping create values cache equations variable list.",err,error,*999)
    ALLOCATE(solverMapping%createValuesCache%dynamicVariableType(solverMapping%numberOfEquationsSets),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver mapping create values cache dynamic variable type.",err,error,*999)
    ALLOCATE(solverMapping%createValuesCache%matrixVariableTypes(0:FIELD_NUMBER_OF_VARIABLE_TYPES, &
      & solverMapping%numberOfEquationsSets,solverMapping%numberOfSolverMatrices),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver mapping create values cache matrix variable types.",err,error,*999)
    ALLOCATE(solverMapping%createValuesCache%residualVariableTypes(0:FIELD_NUMBER_OF_VARIABLE_TYPES,0, &
      & solverMapping%numberOfEquationsSets),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver mapping create values cache residual variable type.",err,error,*999)
    ALLOCATE(solverMapping%createValuesCache%rhsVariableType(solverMapping%numberOfEquationsSets),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver mapping create values cache RHS variable type.",err,error,*999)
    ALLOCATE(solverMapping%createValuesCache%sourceVariableTypes(0:FIELD_NUMBER_OF_VARIABLE_TYPES,0, &
      & solverMapping%numberOfEquationsSets),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver mapping create values cache source variable types.",err,error,*999)
    ALLOCATE(solverMapping%createValuesCache%interfaceVariableList(solverMapping%numberOfSolverMatrices),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver mapping create values cache interface variable list.",err,error,*999)
    ALLOCATE(solverMapping%createValuesCache%interfaceIndices(solverMapping%numberOfEquationsSets),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver mapping create values cache interface condition indices.",err,error,*999)
    DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
      NULLIFY(solverMapping%createValuesCache%equationsVariableList(solverMatrixIdx)%ptr)
      CALL List_CreateStart(solverMapping%createValuesCache%equationsVariableList(solverMatrixIdx)%ptr,err,error,*999)
      CALL List_DataTypeSet(solverMapping%createValuesCache%equationsVariableList(solverMatrixIdx)%ptr,LIST_INTG_TYPE, &
        & err,error,*999)
      CALL List_DataDimensionSet(solverMapping%createValuesCache%equationsVariableList(solverMatrixIdx)%ptr,2,err,error,*999)
      CALL List_CreateFinish(solverMapping%createValuesCache%equationsVariableList(solverMatrixIdx)%ptr,err,error,*999)
      NULLIFY(solverMapping%createValuesCache%interfaceVariableList(solverMatrixIdx)%ptr)
      CALL List_CreateStart(solverMapping%createValuesCache%interfaceVariableList(solverMatrixIdx)%ptr,err,error,*999)
      CALL List_DataTypeSet(solverMapping%createValuesCache%interfaceVariableList(solverMatrixIdx)%ptr,LIST_INTG_TYPE, &
        & err,error,*999)
      CALL List_DataDimensionSet(solverMapping%createValuesCache%interfaceVariableList(solverMatrixIdx)%ptr,2,err,error,*999)
      CALL List_CreateFinish(solverMapping%createValuesCache%interfaceVariableList(solverMatrixIdx)%ptr,err,error,*999)
    ENDDO !solver_idx
    NULLIFY(solverMapping%createValuesCache%equationsRHSVariableList)
    CALL List_CreateStart(solverMapping%createValuesCache%equationsRHSVariableList,err,error,*999)
    CALL List_DataTypeSet(solverMapping%createValuesCache%equationsRHSVariableList,LIST_INTG_TYPE,err,error,*999)
    CALL List_DataDimensionSet(solverMapping%createValuesCache%equationsRHSVariableList,2,err,error,*999)
    CALL List_CreateFinish(solverMapping%createValuesCache%equationsRHSVariableList,err,error,*999)
    DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
      NULLIFY(solverMapping%createValuesCache%interfaceIndices(equationsSetIdx)%ptr)
      CALL List_CreateStart(solverMapping%createValuesCache%interfaceIndices(equationsSetIdx)%ptr,err,error,*999)
      CALL List_DataTypeSet(solverMapping%createValuesCache%interfaceIndices(equationsSetIdx)%ptr,LIST_INTG_TYPE, &
        & err,error,*999)
      CALL List_DataDimensionSet(solverMapping%createValuesCache%interfaceIndices(equationsSetIdx)%ptr,2,err,error,*999)
      CALL List_KeyDimensionSet(solverMapping%createValuesCache%interfaceIndices(equationsSetIdx)%ptr,1,err,error,*999)
      CALL List_CreateFinish(solverMapping%createValuesCache%interfaceIndices(equationsSetIdx)%ptr,err,error,*999)            
    ENDDO !equationsSetIdx
    solverMapping%createValuesCache%dynamicVariableType=0
    solverMapping%createValuesCache%matrixVariableTypes=0
    solverMapping%createValuesCache%residualVariableTypes=0
    solverMapping%createValuesCache%rhsVariableType=0
    solverMapping%createValuesCache%sourceVariableTypes=0
       
    EXITS("SolverMapping_CreateValuesCacheInitialise")
    RETURN
999 CALL SolverMapping_CreateValuesCacheFinalise(solverMapping%createValuesCache,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverMapping_CreateValuesCacheInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_CreateValuesCacheInitialise

  !
  !================================================================================================================================
  !

  !>Adds a variable type from an equations set dependent field to the list of variables for a particular solver matrix of a solver mapping.
  SUBROUTINE SolverMapping_CreateValuesCacheEqnVarListAdd(solverMapping,solverMatrixIdx,equationsSetIdx,variableType,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer the solver mapping to add to the var list for.
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index of the variable list
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index of the variable to add
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx2,numberOfVariables,variableIdx,variableItem(2)
    LOGICAL :: variableFound
    TYPE(EquationsSetType), POINTER :: equationsSet,variableEquationsSet
    TYPE(FieldType), POINTER :: dependentField,variableDependentField
    TYPE(ListType),  POINTER :: equationsVariableList
    TYPE(RegionType), POINTER :: dependentRegion,variableRegion,variableDependentRegion
    TYPE(SolverMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverMapping_CreateValuesCacheEqnVarListAdd",err,error,*999)

    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)
    NULLIFY(createValuesCache)
    CALL SolverMapping_CreateValuesCacheGet(solverMapping,createValuesCache,err,error,*999)
    NULLIFY(equationsVariableList)
    CALL SolverMappingCVC_EquationsVariableListGet(createValuesCache,solverMatrixIdx,equationsVariableList, &
      & err,error,*999)
    NULLIFY(equationsSet)
    CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(dependentRegion)
    CALL Field_RegionGet(dependentField,dependentRegion,err,error,*999)    
    IF(variableType/=0) THEN
      variableFound=.FALSE.
      CALL List_NumberOfItemsGet(equationsVariableList,numberOfVariables,err,error,*999)
      DO variableIdx=1,numberOfVariables
        CALL List_ItemGet(equationsVariableList,variableIdx,variableItem,err,error,*999)
        equationsSetIdx2=variableItem(1)
        IF(variableType==variableItem(2)) THEN
          NULLIFY(variableEquationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx2,variableEquationsSet,err,error,*999)
          NULLIFY(variableDependentField)
          CALL EquationsSet_DependentFieldGet(variableEquationsSet,variableDependentField,err,error,*999)
          IF(ASSOCIATED(dependentField,variableDependentField)) THEN
            NULLIFY(variableDependentRegion)
            CALL Field_RegionGet(variableDependentField,variableDependentRegion,err,error,*999)
            IF(ASSOCIATED(dependentRegion,variableDependentRegion)) THEN
              variableFound=.TRUE.
              EXIT
            ENDIF
          ENDIF
        ENDIF
      ENDDO !variableIdx
      IF(.NOT.variableFound) THEN
        variableItem(1)=equationsSetIdx
        variableItem(2)=variableType
        CALL List_ItemAdd(equationsVariableList,variableItem,err,error,*999)
      ENDIF
    ENDIF
       
    EXITS("SolverMapping_CreateValuesCacheEqnVarListAdd")
    RETURN
999 ERRORS("SolverMapping_CreateValuesCacheEqnVarListAdd",err,error)    
    EXITS("SolverMapping_CreateValuesCacheEqnVarListAdd")    
    RETURN 1
   
  END SUBROUTINE SolverMapping_CreateValuesCacheEqnVarListAdd

  !
  !================================================================================================================================
  !

  !>Adds a variable type from an equations set dependent field to the list of RHS variables for a particular solver mapping.
  SUBROUTINE SolverMapping_CreateValuesCacheEqnRHSVarListAdd(solverMapping,equationsSetIdx,rhsVariableType,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer the solver mapping to add to the variable list for.
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index of the variable to add
    INTEGER(INTG), INTENT(IN) :: rhsVariableType !<The RHS variable type to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx2,numberOfVariables,variableIdx,variableItem(2)
    LOGICAL :: variableFound
    TYPE(EquationsSetType), POINTER :: equationsSet,variableEquationsSet
    TYPE(FieldType), POINTER :: dependentField,variableDependentField
    TYPE(FieldVariableType), POINTER :: rhsVariable,variableRhsVariable
    TYPE(ListType), POINTER :: rhsVariableList
    TYPE(RegionType), POINTER :: dependentRegion
    TYPE(SolverMappingCreateValuesCacheType), POINTER :: createValuesCache

    ENTERS("SolverMapping_CreateValuesCacheEqnRHSVarListAdd",err,error,*999)

    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)

    NULLIFY(createValuesCache)
    CALL SolverMapping_CreateValuesCacheGet(solverMapping,createValuesCache,err,error,*999)
    NULLIFY(rhsVariableList)
    CALL SolverMappingCVC_RHSVariableListGet(createValuesCache,rhsVariableList,err,error,*999)
    NULLIFY(equationsSet)
    CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(rhsVariable)
    CALL Field_VariableGet(dependentField,rhsVariableType,rhsVariable,err,error,*999)
    NULLIFY(dependentRegion)
    CALL Field_RegionGet(dependentField,dependentRegion,err,error,*999)
    variableFound=.FALSE.
    CALL List_NumberOfItemsGet(rhsVariableList,numberOfVariables,err,error,*999)
    DO variableIdx=1,numberOfVariables
      CALL List_ItemGet(rhsVariableList,variableIdx,variableItem,err,error,*999)
      equationsSetIdx2=variableItem(1)
      NULLIFY(variableEquationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx2,variableEquationsSet,err,error,*999)
      NULLIFY(variableDependentField)
      CALL EquationsSet_DependentFieldGet(variableEquationsSet,variableDependentField,err,error,*999)
      NULLIFY(variableRHSVariable)
      CALL Field_VariableGet(variableDependentField,variableItem(2),variableRHSVariable,err,error,*999)
      IF(ASSOCIATED(rhsVariable,variableRHSVariable)) THEN
        variableFound=.TRUE.
        EXIT
      ENDIF
    ENDDO !variableIdx
    IF(.NOT.variableFound) THEN
      variableItem(1)=equationsSetIdx
      variableItem(2)=rhsVariableType
      CALL List_ItemAdd(rhsVariableList,variableItem,err,error,*999)
    ENDIF
       
    EXITS("SolverMapping_CreateValuesCacheEqnRHSVarListAdd")
    RETURN
999 ERRORS("SolverMapping_CreateValuesCacheEqnRHSVarListAdd",err,error)    
    EXITS("SolverMapping_CreateValuesCacheEqnRHSVarListAdd")    
    RETURN 1
   
  END SUBROUTINE SolverMapping_CreateValuesCacheEqnRHSVarListAdd

  !
  !================================================================================================================================
  !

  !>Adds a variable type from an interface condition Lagrange field to the list of variables for a particular solver matrix of a solver mapping.
  SUBROUTINE SolverMapping_CreateValuesCacheInterfVarListAdd(solverMapping,solverMatrixIdx,interfaceConditionIdx, &
    & variableType,err,error,*)
    
    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer the solver mapping to add to the var list for.
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index of the variable list
    INTEGER(INTG), INTENT(IN) :: interfaceConditionIdx !<The interface condition index of the variable to add
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceConditionIdx2,numberOfVariables,variableIdx,variableItem(2)
    LOGICAL :: variableFound
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition,variableInterfaceCondition
    TYPE(FieldType), POINTER :: lagrangeField,variableLagrangeField
    TYPE(InterfaceType), POINTER :: lagrangeInterface,variableLagrangeInterface
    TYPE(ListType), POINTER :: interfaceVariableList
    TYPE(SolverMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverMapping_CreateValuesCacheInterfVarListAdd",err,error,*999)

    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)
    CALL SolverMapping_CreateValuesCacheGet(solverMapping,createValuesCache,err,error,*999)
    NULLIFY(interfaceVariableList)
    CALL SolverMappingCVC_InterfaceVariableListGet(createValuesCache,solverMatrixIdx,interfaceVariableList, &
      & err,error,*999)
    NULLIFY(interfaceCondition)
    CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
    SELECT CASE(interfaceCondition%method)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
      NULLIFY(lagrangeField)
      CALL InterfaceCondition_LagrangeFieldGet(interfaceCondition,lagrangeField,err,error,*999)
      NULLIFY(lagrangeInterface)
      CALL Field_InterfaceGet(lagrangeField,lagrangeInterface,err,error,*999)
      IF(variableType/=0) THEN
        variableFound=.FALSE.
        CALL List_NumberOfItemsGet(interfaceVariableList,numberOfVariables,err,error,*999)
        DO variableIdx=1,numberOfVariables
          CALL List_ItemGet(interfaceVariableList,variableIdx,variableItem,err,error,*999)
          IF(variableType==variableItem(2)) THEN
            interfaceConditionIdx2=variableItem(1)
            NULLIFY(variableInterfaceCondition)
            CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx2,variableInterfaceCondition,err,error,*999)
            SELECT CASE(variableInterfaceCondition%method)
            CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
              NULLIFY(variableLagrangeField)
              CALL InterfaceCondition_LagrangeFieldGet(variableInterfaceCondition,variableLagrangeField,err,error,*999)
              IF(ASSOCIATED(lagrangeField,variableLagrangeField)) THEN
                NULLIFY(variableLagrangeInterface)
                CALL Field_InterfaceGet(variableLagrangeField,variableLagrangeInterface,err,error,*999)
                IF(ASSOCIATED(lagrangeInterface,variableLagrangeInterface)) THEN
                  variableFound=.TRUE.
                  EXIT
                ENDIF
              ENDIF
            CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="The interface condition method of "// &
                & TRIM(NumberToVString(variableInterfaceCondition%METHOD,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDIF
        ENDDO !variableIdx
        IF(.NOT.variableFound) THEN
          variableItem(1)=interfaceConditionIdx
          variableItem(2)=variableType
          CALL List_ItemAdd(InterfaceVariableList,variableItem,err,error,*999)
        ENDIF
      ENDIF
    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The interface condition method of "// &
        & TRIM(NumberToVString(interfaceCondition%METHOD,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("SolverMapping_CreateValuesCacheInterfVarListAdd")
    RETURN
999 ERRORS("SolverMapping_CreateValuesCacheInterfVarListAdd",err,error)
    EXITS("SolverMapping_CreateValuesCacheInterfVarListAdd")
    RETURN 1
   
  END SUBROUTINE SolverMapping_CreateValuesCacheInterfVarListAdd

  !
  !================================================================================================================================
  !

  !>Destroy a solver mapping.
  SUBROUTINE SolverMapping_Destroy(solverMapping,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer the solver mapping to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverMapping_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated",err,error,*999)
    
    CALL SolverMapping_Finalise(solverMapping,err,error,*999)
        
    EXITS("SolverMapping_Destroy")
    RETURN
999 ERRORSEXITS("SolverMapping_Destroy",err,error)    
    RETURN 1
   
  END SUBROUTINE SolverMapping_Destroy

  !
  !================================================================================================================================
  !

  !>Sets/changes the mapping of global variables to a solver matrix for the solver mapping
  SUBROUTINE SolverMapping_EquatsVarsToSolverMatrixSet(solverMapping,solverMatrixNumber,equationsSetIdx,variableTypes,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping
    INTEGER(INTG), INTENT(IN) :: solverMatrixNumber !<The solver matrix number to set the equations variables for
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index in the solver mapping to specify the variable types for
    INTEGER(INTG), INTENT(IN) :: variableTypes(:) !<The variable types to map to the solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIndex,linearMappingIndex,linearVariableIndex,variableIdx
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(SolverMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverMapping_EquatsVarsToSolverMatrixSet",err,error,*999)

    CALL SolverMapping_AssertNotFinished(solverMapping,err,error,*999)
    IF(SIZE(variableTypes,1)<1.OR.SIZE(variableTypes,1)>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The size of the specified variable types array of "// &
        & TRIM(NumberToVString(SIZE(variableTypes,1),"*",err,error))// &
        & " is invalid. The size must be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(solverMatrixNumber<1.OR.solverMatrixNumber>solverMapping%numberOfSolverMatrices) THEN
      localError="The solver matrix number of "//TRIM(NumberToVString(solverMatrixNumber,"*",err,error))// &
        & " is invalid. The number must be >= 1 and <= "// &
        & TRIM(NumberToVString(solverMapping%numberOfSolverMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(createValuesCache)
    CALL SolverMapping_CreateValuesCacheGet(solverMapping,createValuesCache,err,error,*999)
    NULLIFY(equationsSet)
    CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(linearMapping)
    CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
    DO variableIdx=1,SIZE(variableTypes,1)
!!TODO: CHECK THAT THE VARIABLE TYPE IS NOT REPEATED
      CALL EquationsMappingLinear_LinearVariableIndexGet(linearMapping,variableTypes(variableIdx),linearMappingIndex, &
        & err,error,*999)
       IF(linearMappingIndex==0) THEN
         localError="The variable type of "//TRIM(NumberToVString(variableTypes(variableIdx),"*",err,error))// &
           & " at position "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
           & " in the array is invalid. That variable type is not mapped to any equations matrices."
       ENDIF
     ENDDO !variableIdx
     createValuesCache%matrixVariableTypes(0,equationsSetIndex,solverMatrixNumber)=SIZE(variableTypes,1)
     createValuesCache%matrixVariableTypes(1:SIZE(variableTypes,1),equationsSetIndex,solverMatrixNumber)=variableTypes
   
    EXITS("SolverMapping_EquatsVarsToSolverMatrixSet")
    RETURN
999 ERRORSEXITS("SolverMapping_EquatsVarsToSolverMatrixSet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_EquatsVarsToSolverMatrixSet
  
  !
  !================================================================================================================================
  !

  !>Adds an equations set to a solver mapping
  SUBROUTINE SolverMapping_EquationsSetAdd(solverMapping,equationsSet,equationsSetIndex,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer the solver mapping to add the equations set to
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to add
    INTEGER(INTG), INTENT(OUT) :: equationsSetIndex !<On exit, the index of the equations set in the solver mapping
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dynamicVariableType,equationsSetIdx,linearVariableType,matrixIdx,newNumberOfResiduals,newNumberOfSources, &
      & numberOfEquationsMatrices,numberOfLinearVariables,numberOfResiduals,numberOfResidualVariables,numberOfSources, &
      & numberOfSourceVariables,residualIdx,residualVariableType,rhsVariableType,solverMatrixIdx,sourceIdx,sourceVariableType, &
      & variableIdx,variableIndex,variableTypeIdx
    INTEGER(INTG), ALLOCATABLE :: newDynamicVariableTypes(:),newMatrixVariableTypes(:,:,:),newRHSVariableTypes(:), &
      & newResidualVariableTypes(:,:,:),newSourceVariableTypes(:,:,:)
    LOGICAL :: matrixDone
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsSetType), POINTER :: equationsSet2
    TYPE(EquationsSetPtrType), ALLOCATABLE :: newEquationsSets(:)
    TYPE(ListPtrType), POINTER :: newInterfaceIndices(:)
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: localError

    NULLIFY(newInterfaceIndices)
    
    ENTERS("SolverMapping_EquationsSetAdd",err,error,*999)

    equationsSetIndex=0

    CALL SolverMapping_AssertNotFinished(solverMapping,err,error,*999)
    CALL EquationsSet_AssertIsFinished(equationsSet,err,error,*999)
    NULLIFY(createValuesCache)
    CALL SolverMapping_CreateValuesCacheGet(solverMapping,createValuesCache,err,error,*999)
    NULLIFY(solverEquations)
    CALL SolverMapping_SolverEquationsGet(solverMapping,solverEquations,err,error,*999)
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    !Check that the equations set has not already been added
    DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
      NULLIFY(equationsSet2)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet2,err,error,*999)
      IF(ASSOCIATED(equationsSet,equationsSet2)) THEN
        localError="Equations set number "//TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
        IF(ASSOCIATED(equationsSet%region)) &
          & localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
        localError=localError//" has already been added at equations set index position "// &
          & TRIM(NumberToVString(equationsSetIdx,"*",err,error))//" for the solver mapping."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDDO !equationsSetIdx
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingExists(vectorMapping,nonlinearMapping,err,error,*999)
    newNumberOfResiduals=SIZE(solverMapping%createValuesCache%residualVariableTypes,2)
    IF(ASSOCIATED(nonlinearMapping)) THEN
      CALL EquationsMappingNonlinear_NumberOfResidualsGet(nonlinearMapping,numberOfResiduals,err,error,*999)
      IF(numberOfResiduals>newNumberOfResiduals) newNumberOfResiduals=numberOfResiduals
    ENDIF
    NULLIFY(sourcesMapping)
    CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*999)
    newNumberOfSources=SIZE(solverMapping%createValuesCache%sourceVariableTypes,2)
    IF(ASSOCIATED(sourcesMapping)) THEN
      CALL EquationsMappingSources_NumberOfSourcesGet(sourcesMapping,numberOfSources,err,error,*999)
      IF(numberOfSources>newNumberOfSources) newNumberOfSources=numberOfSources
    ENDIF
    ALLOCATE(newDynamicVariableTypes(solverMapping%numberOfEquationsSets+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new dynamic variable types.",err,error,*999)
    ALLOCATE(newMatrixVariableTypes(0:FIELD_NUMBER_OF_VARIABLE_TYPES,solverMapping%numberOfEquationsSets+1, &
      & solverMapping%numberOfSolverMatrices),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new matrix variable types.",err,error,*999)
    ALLOCATE(newResidualVariableTypes(0:FIELD_NUMBER_OF_VARIABLE_TYPES,newNumberOfResiduals, &
      & solverMapping%numberOfEquationsSets+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new residual variable types.",err,error,*999)
    ALLOCATE(newRHSVariableTypes(solverMapping%numberOfEquationsSets+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new RHS variable types.",err,error,*999)
    ALLOCATE(newSourceVariableTypes(0:FIELD_NUMBER_OF_VARIABLE_TYPES,newNumberOfSources, &
      & solverMapping%numberOfEquationsSets+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new source variable types.",err,error,*999)
    ALLOCATE(newInterfaceIndices(solverMapping%numberOfEquationsSets+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new interface indices.",err,error,*999)
    ALLOCATE(newEquationsSets(solverMapping%numberOfEquationsSets+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new equations sets.",err,error,*999)              
    IF(solverMapping%numberOfEquationsSets>0) THEN
      DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
        newInterfaceIndices(equationsSetIdx)%ptr=>createValuesCache%interfaceIndices(equationsSetIdx)%ptr
      ENDDO !equations_sets
      newDynamicVariableTypes(1:solverMapping%numberOfEquationsSets)=createValuesCache%dynamicVariableType
      newMatrixVariableTypes(:,1:solverMapping%numberOfEquationsSets,:)=createValuesCache%matrixVariableTypes
      newResidualVariableTypes(:,:,1:solverMapping%numberOfEquationsSets)=createValuesCache%residualVariableTypes
      newRHSVariableTypes(1:solverMapping%numberOfEquationsSets)=createValuesCache%rhsVariableType
      newSourceVariableTypes(:,:,1:solverMapping%numberOfEquationsSets)=createValuesCache%sourceVariableTypes
      DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
        newEquationsSets(equationsSetIdx)%ptr=>solverMapping%equationsSets(equationsSetIdx)%ptr
      ENDDO !equationsSetIdx              
      newDynamicVariableTypes(solverMapping%numberOfEquationsSets+1)=0
      newMatrixVariableTypes(:,solverMapping%numberOfEquationsSets+1,:)=0
      newResidualVariableTypes(:,:,solverMapping%numberOfEquationsSets+1)=0
      newRHSVariableTypes(solverMapping%numberOfEquationsSets+1)=0
      newSourceVariableTypes(:,:,solverMapping%numberOfEquationsSets+1)=0            
    ELSE IF(solverMapping%numberOfEquationsSets==0) THEN      
      newDynamicVariableTypes=0
      newMatrixVariableTypes=0
      newResidualVariableTypes=0
      newRHSVariableTypes=0
      newSourceVariableTypes=0      
    ELSE
      CALL FlagError("The number of equations sets is < 0.",err,error,*999)
    ENDIF
    IF(ASSOCIATED(createValuesCache%interfaceIndices)) DEALLOCATE(createValuesCache%interfaceIndices)              
    CALL MOVE_ALLOC(newDynamicVariableTypes,createValuesCache%dynamicVariableType)
    CALL MOVE_ALLOC(newMatrixVariableTypes,createValuesCache%matrixVariableTypes)
    CALL MOVE_ALLOC(newResidualVariableTypes,createValuesCache%residualVariableTypes)
    CALL MOVE_ALLOC(newRHSVariableTypes,createValuesCache%rhsVariableType)
    CALL MOVE_ALLOC(newSourceVariableTypes,createValuesCache%sourceVariableTypes)
    CALL MOVE_ALLOC(newEquationsSets,solverMapping%equationsSets)
    createValuesCache%interfaceIndices=>newInterfaceIndices
    NULLIFY(createValuesCache%interfaceIndices(solverMapping%numberOfEquationsSets+1)%ptr)
    CALL List_CreateStart(createValuesCache%interfaceIndices(solverMapping%numberOfEquationsSets+1)%ptr,err,error,*999)
    CALL List_DataTypeSet(createValuesCache%interfaceIndices(solverMapping%numberOfEquationsSets+1)%ptr,LIST_INTG_TYPE, &
      & err,error,*999)
    CALL List_DataDimensionSet(createValuesCache%interfaceIndices(solverMapping%numberOfEquationsSets+1)%ptr,2,err,error,*999)
    CALL List_KeyDimensionSet(createValuesCache%interfaceIndices(solverMapping%numberOfEquationsSets+1)%ptr,1,err,error,*999)
    CALL List_CreateFinish(createValuesCache%interfaceIndices(solverMapping%numberOfEquationsSets+1)%ptr,err,error,*999)
    numberOfResiduals=0
    SELECT CASE(equations%timeDependence)
    CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC)
      SELECT CASE(equations%linearity)
      CASE(EQUATIONS_LINEAR)
        NULLIFY(linearMapping)
        CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
        !Linear matrices to map.

!! TODO: work out what variables should be mapped and what variables go to the RHS???
        
        !Map the first matrix variable found in the equations set to the first solver matrix, the second
        !variable found to the second, etc.
        linearVariableType=1
        DO matrixIdx=1,solverMapping%numberOfSolverMatrices
          matrixDone=.FALSE.
          DO WHILE(linearVariableType<=FIELD_NUMBER_OF_VARIABLE_TYPES.AND..NOT.matrixDone)
            CALL EquationsMappingLinear_LinearVariableIndexGet(linearMapping,linearVariableType,variableIndex,err,error,*999)
            IF(variableIndex==0) THEN
              linearVariableType=linearVariableType+1
            ELSE
              CALL EquationsMappingLinear_VariableNumberOfMatricesGet(linearMapping,variableIndex,numberOfEquationsMatrices, &
                & err,error,*999)
              IF(numberOfEquationsMatrices>0) THEN
                createValuesCache%matrixVariableTypes(0,solverMapping%numberOfEquationsSets+1,matrixIdx)=1
                createValuesCache%matrixVariableTypes(linearVariableType,solverMapping%numberOfEquationsSets+1,matrixIdx)= &
                  & linearVariableType
                matrixDone=.TRUE.
              ELSE
                linearVariableType=linearVariableType+1
              ENDIF
            ENDIF
          ENDDO
          IF(.NOT.matrixDone) THEN
            !Error - could not find any more variables to map to this solver matrix
            localError="Could not find any unmapped variables for solver matrix "// &
              & TRIM(NumberToVString(matrixIdx,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDDO !matrixIdx
      CASE(EQUATIONS_NONLINEAR)
        NULLIFY(nonlinearMapping)
        CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
        CALL EquationsMappingNonlinear_NumberOfResidualsGet(nonlinearMapping,numberOfResiduals,err,error,*999)
        NULLIFY(linearMapping)
        CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
        DO residualIdx=1,numberOfResiduals
          NULLIFY(residualMapping)
          CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
          CALL EquationsMappingResidual_NumberOfResidualVariablesGet(residualMapping,numberOfResidualVariables,err,error,*999)
          !Set the number of residual variables for this equations set
          createValuesCache%residualVariableTypes(0,residualIdx,solverMapping%numberOfEquationsSets+1)=numberOfResidualVariables
          !Map the residual variables to the solver Jacobian
          DO variableIdx=1,numberOfResidualVariables
            CALL EquationsMappingResidual_VariableTypeGet(residualMapping,variableIdx,residualVariableType,err,error,*999)
            createValuesCache%residualVariableTypes(residualVariableType,residualIdx,solverMapping%numberOfEquationsSets+1)= &
              & residualVariableType
          ENDDO !variableIdx
          IF(ASSOCIATED(linearMapping)) THEN
            !If there are linear matrices operating on the residual variable then map them to the
            !solver matrix (Jacobian)
            IF(solverMapping%numberOfSolverMatrices/=1) THEN
              localError="Invalid number of solver matrices. For nonlinear solver equations there should "// &
                & "be 1 solver matrix and there are "// &
                & TRIM(NumberToVString(solverMapping%numberOfSolverMatrices,"*",err,error))// &
                & " solver matrices."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            DO variableIdx=1,numberOfResidualVariables
              CALL EquationsMappingResidual_VariableTypeGet(residualMapping,variableIdx,residualVariableType,err,error,*999)
              CALL EquationsMappingResidual_VariableIndexGet(residualMapping,residualVariableType,variableIndex,err,error,*999)
              IF(variableIndex>0) THEN
                !There is a linear matrix with the same variable as the residual variable so map it to the solver matrix
                CALL EquationsMappingLinear_VariableNumberOfMatricesGet(linearMapping,variableIndex,numberOfEquationsMatrices, &
                  & err,error,*999)
                IF(numberOfEquationsMatrices>0) THEN
                  createValuesCache%matrixVariableTypes(residualVariableType,solverMapping%numberOfEquationsSets+1,1)= &
                    & residualVariableType
                ENDIF
              ENDIF
            ENDDO !matrixIdx
          ENDIF
        ENDDO !residualIdx
      CASE DEFAULT
        localError="The equations linearity type of "//TRIM(NumberToVString(solverEquations%linearity,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
      IF(solverMapping%numberOfSolverMatrices/=1) THEN
        localError="Invalid number of solver matrices. For dynamic solver equations there should "// &
          & "be 1 solver matrix and there are "//TRIM(NumberToVString(solverMapping%numberOfSolverMatrices,"*",err,error))// &
          & " solver matrices."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(dynamicMapping)
      CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      CALL EquationsMappingDynamic_DynamicVariableTypeGet(dynamicMapping,dynamicVariableType,err,error,*999)
      createValuesCache%dynamicVariableType(solverMapping%numberOfEquationsSets+1)=dynamicVariableType
      SELECT CASE(equations%linearity)
      CASE(EQUATIONS_LINEAR)
        !Look at linear equations below
      CASE(EQUATIONS_NONLINEAR)
        NULLIFY(nonlinearMapping)
        CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
        CALL EquationsMappingNonlinear_NumberOfResidualsGet(nonlinearMapping,numberOfResiduals,err,error,*999)
        DO residualIdx=1,numberOfResiduals
          NULLIFY(residualMapping)
          CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
          CALL EquationsMappingResidual_NumberOfResidualVariablesGet(residualMapping,numberOfResidualVariables,err,error,*999)
          DO variableIdx=1,numberOfResidualVariables
            CALL EquationsMappingResidual_VariableIndexGet(residualMapping,dynamicVariableType,variableIndex,err,error,*999)
            IF(variableIndex>0) THEN
              createValuesCache%residualVariableTypes(0,residualIdx,solverMapping%numberOfEquationsSets+1)=1
              createValuesCache%residualVariableTypes(dynamicVariableType,residualIdx,solverMapping%numberOfEquationsSets+1)= &
                & dynamicVariableType
            ENDIF
          ENDDO !variableIdx
        ENDDO !residualIdx
      CASE DEFAULT
        localError="The equations linearity type of "//TRIM(NumberToVString(solverEquations%linearity,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
      IF(ASSOCIATED(linearMapping)) THEN
        CALL EquationsMappingLinear_NumberOfLinearVariablesGet(linearMapping,numberOfLinearVariables,err,error,*999)
        DO variableIdx=1,numberOfLinearVariables
          CALL EquationsMappingLinear_LinearVariableIndexGet(linearMapping,dynamicVariableType,variableIndex,err,error,*999)
          IF(variableIndex>0) THEN
            !There is a linear variable that is the same as the dynamic variable so map the variable type
            createValuesCache%matrixVariableTypes(0,solverMapping%numberOfEquationsSets+1,1)=1
            createValuesCache%matrixVariableTypes(dynamicVariableType,solverMapping%numberOfEquationsSets+1,1)= &
              & dynamicVariableType
          ENDIF
        ENDDO !variableIdx
      ENDIF
    CASE DEFAULT
      localError="The equations time dependence type of "// &
        & TRIM(NumberToVString(solverEquations%timeDependence,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    NULLIFY(rhsMapping)
    CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
    IF(ASSOCIATED(rhsMapping)) THEN
      CALL EquationsMappingRHS_RHSVariableTypeGet(rhsMapping,rhsVariableType,err,error,*999)
      createValuesCache%rhsVariableType(solverMapping%numberOfEquationsSets+1)=rhsVariableType
    ENDIF
    NULLIFY(sourcesMapping)
    CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*999)
    IF(ASSOCIATED(sourcesMapping)) THEN
      CALL EquationsMappingSources_NumberOfSourcesGet(sourcesMapping,numberOfSources,err,error,*999)
      DO sourceIdx=1,numberOfSources
        NULLIFY(sourceMapping)
        CALL EquationsMappingSources_SourceMappingGet(sourcesMapping,sourceIdx,sourceMapping,err,error,*999)
        CALL EquationsMappingSource_SourceVariableTypeGet(sourceMapping,sourceVariableType,err,error,*999)
        createValuesCache%sourceVariableTypes(sourceVariableType,sourceIdx,solverMapping%numberOfEquationsSets+1)= &
          & sourceVariableType
      ENDDO !sourceIdx
    ENDIF
    !Add the equations set to the list of equations sets
    solverMapping%equationsSets(solverMapping%numberOfEquationsSets+1)%ptr=>equationsSet
    solverMapping%numberOfEquationsSets=solverMapping%numberOfEquationsSets+1
    equationsSetIndex=solverMapping%numberOfEquationsSets
    
    !Add the variables to the list of variables
    dynamicVariableType=createValuesCache%dynamicVariableType(equationsSetIndex)
    CALL SolverMapping_CreateValuesCacheEqnVarListAdd(solverMapping,1,equationsSetIndex,dynamicVariableType,err,error,*999)
    DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
      DO variableTypeIdx=1,FIELD_NUMBER_OF_VARIABLE_TYPES
        linearVariableType=createValuesCache%matrixVariableTypes(variableTypeIdx,equationsSetIndex,solverMatrixIdx)
        CALL SolverMapping_CreateValuesCacheEqnVarListAdd(solverMapping,1,equationsSetIndex,linearVariableType,err,error,*999)
      ENDDO !matrixIdx
    ENDDO !solverMatrixIdx
    DO residualIdx=1,numberOfResiduals
      DO variableTypeIdx=1,FIELD_NUMBER_OF_VARIABLE_TYPES
        residualVariableType=createValuesCache%residualVariableTypes(variableTypeIdx,residualIdx,equationsSetIndex)
        CALL SolverMapping_CreateValuesCacheEqnVarListAdd(solverMapping,1,equationsSetIndex,residualVariableType,err,error,*999)
      ENDDO !variableTypeIdx
    ENDDO !residualIdx
    rhsVariableType=createValuesCache%rhsVariableType(equationsSetIndex)
    IF(rhsVariableType/=0) &
      & CALL SolverMapping_CreateValuesCacheEqnRHSVarListAdd(solverMapping,equationsSetIndex,rhsVariableType,err,error,*999)
        
    EXITS("SolverMapping_EquationsSetAdd")
    RETURN
999 IF(ALLOCATED(newMatrixVariableTypes)) DEALLOCATE(newMatrixVariableTypes)
    IF(ALLOCATED(newResidualVariableTypes)) DEALLOCATE(newResidualVariableTypes)
    IF(ALLOCATED(newRHSVariableTypes)) DEALLOCATE(newRHSVariableTypes)
    IF(ALLOCATED(newSourceVariableTypes)) DEALLOCATE(newSourceVariableTypes)
    IF(ALLOCATED(newEquationsSets)) DEALLOCATE(newEquationsSets)
    ERRORSEXITS("SolverMapping_EquationsSetAdd",err,error)    
    RETURN 1
   
  END SUBROUTINE SolverMapping_EquationsSetAdd

  !
  !================================================================================================================================
  !

  !>Finalises a equations set to solver matrices map and deallocates all memory.
  SUBROUTINE SolverMappingESToSMSMap_Finalise(equationsSetToSolverMatricesMap,err,error,*)

    !Argument variables
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap !<The equations set to solver matrices map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationMatrixIdx,interfaceConditionIdx,jacobianMatrixIdx,rowIdx,solverMatrixIdx
    
    ENTERS("SolverMappingESToSMSMap_Finalise",err,error,*999)

    IF(ASSOCIATED(equationsSetToSolverMatricesMap)) THEN
      IF(ALLOCATED(equationsSetToSolverMatricesMap%interfaceConditionToEquationsSetMaps)) THEN
        DO interfaceConditionIdx=1,SIZE(equationsSetToSolverMatricesMap%interfaceConditionToEquationsSetMaps,1)
          CALL SolverMappingSMICToESMap_Finalise(equationsSetToSolverMatricesMap% &
            & interfaceConditionToEquationsSetMaps(interfaceConditionIdx),err,error,*999)
        ENDDO !interfaceConditionIdx
        DEALLOCATE(equationsSetToSolverMatricesMap%interfaceConditionToEquationsSetMaps)
      ENDIF
      IF(ALLOCATED(equationsSetToSolverMatricesMap%equationsMatricesToSolverMatrixMaps)) THEN
        DO solverMatrixIdx=1,SIZE(equationsSetToSolverMatricesMap%equationsMatricesToSolverMatrixMaps,1)
          CALL SolverMappingEMSToSMMap_Finalise(equationsSetToSolverMatricesMap% &
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr,err,error,*999)
        ENDDO !solverMatrixIdx
        DEALLOCATE(equationsSetToSolverMatricesMap%equationsMatricesToSolverMatrixMaps)
      ENDIF
      IF(ALLOCATED(equationsSetToSolverMatricesMap%equationsMatrixToSolverMatricesMaps)) THEN
        DO equationMatrixIdx=1,SIZE(equationsSetToSolverMatricesMap%equationsMatrixToSolverMatricesMaps,1)
          CALL SolverMappingEMToSMSMap_Finalise(equationsSetToSolverMatricesMap% &
            & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr,err,error,*999)
        ENDDO !equationMatrixIdx
        DEALLOCATE(equationsSetToSolverMatricesMap%equationsMatrixToSolverMatricesMaps)
      ENDIF
      IF(ALLOCATED(equationsSetToSolverMatricesMap%jacobianMatrixToSolverMatrixMaps)) THEN
        DO jacobianMatrixIdx=1,SIZE(equationsSetToSolverMatricesMap%jacobianMatrixToSolverMatrixMaps,1)
          CALL SolverMappingJMToSMMap_Finalise(equationsSetToSolverMatricesMap% &
            & jacobianMatrixToSolverMatrixMaps(jacobianMatrixIdx)%ptr,err,error,*999)
        ENDDO !jacobianMatrixIdx
        DEALLOCATE(equationsSetToSolverMatricesMap%jacobianMatrixToSolverMatrixMaps)
      ENDIF
      IF(ASSOCIATED(equationsSetToSolverMatricesMap%equationsRowToSolverRowsMap)) THEN
        DO rowIdx=1,SIZE(equationsSetToSolverMatricesMap%equationsRowToSolverRowsMap,1)
          CALL MatrixRowColCoupling_Finalise(equationsSetToSolverMatricesMap%equationsRowToSolverRowsMap(rowIdx),err,error,*999)
        ENDDO !rowIdx
        DEALLOCATE(equationsSetToSolverMatricesMap%equationsRowToSolverRowsMap)
      ENDIF
      DEALLOCATE(equationsSetToSolverMatricesMap)
    ENDIF
        
    EXITS("SolverMappingESToSMSMap_Finalise")
    RETURN
999 ERRORS("SolverMappingESToSMSMap_Finalise",err,error)    
    EXITS("SolverMappingESToSMSMap_Finalise")    
    RETURN 1
   
  END SUBROUTINE SolverMappingESToSMSMap_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises a equations set to solver map.
  SUBROUTINE SolverMappingESToSMSMap_Initialise(equationsSetToSolverMatricesMap,err,error,*)

    !Argument variables
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap !<The equations set to solver map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverMappingESToSMSMap_Initialise",err,error,*998)

    IF(ASSOCIATED(equationsSetToSolverMatricesMap)) &
      & CALL FlagError("Equations set to solver matrices map is already associated.",err,error,*998)

    ALLOCATE(equationsSetToSolverMatricesMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set to solver matrices map.",err,error,*999)
    equationsSetToSolverMatricesMap%equationsSetIndex=0
    NULLIFY(equationsSetToSolverMatricesMap%solverMapping)
    NULLIFY(equationsSetToSolverMatricesMap%equations)
    equationsSetToSolverMatricesMap%numberOfInterfaceConditions=0
    equationsSetToSolverMatricesMap%numberOfSolverMatrices=0
    equationsSetToSolverMatricesMap%numberOfEquationsMatrices=0
    equationsSetToSolverMatricesMap%numberOfJacobianMatrices=0

    EXITS("SolverMappingESToSMSMap_Initialise")
    RETURN
999 CALL SolverMappingESToSMSMap_Finalise(equationsSetToSolverMatricesMap,dummyErr,dummyError,*998)
998 ERRORS("SolverMappingESToSMSMap_Initialise",err,error)    
    EXITS("SolverMappingESToSMSMap_Initialise")    
    RETURN 1
   
  END SUBROUTINE SolverMappingESToSMSMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises a equations matrix to solver matrix map and deallocates all memory.
  SUBROUTINE SolverMappingEMToSMMap_Finalise(equationsMatrixToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: equationsMatrixToSolverMatrixMap !<The equations set to solver map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx
    
    ENTERS("SolverMappingEMToSMMap_Finalise",err,error,*999)

    IF(ASSOCIATED(equationsMatrixToSolverMatrixMap)) THEN
      IF(ASSOCIATED(equationsMatrixToSolverMatrixMap%equationsColToSolverColsMap)) THEN
        DO columnIdx=1,SIZE(equationsMatrixToSolverMatrixMap%equationsColToSolverColsMap,1)
          CALL MatrixRowColCoupling_Finalise(equationsMatrixToSolverMatrixMap%equationsColToSolverColsMap(columnIdx),err,error,*999)
        ENDDO !columnIdx
        DEALLOCATE(equationsMatrixToSolverMatrixMap%equationsColToSolverColsMap)
      ENDIF
      DEALLOCATE(equationsMatrixToSolverMatrixMap)
    ENDIF
        
    EXITS("SolverMappingEMToSMMap_Finalise")
    RETURN
999 ERRORSEXITS("SolverMappingEMToSMMap_Finalise",err,error)    
    RETURN 1
   
  END SUBROUTINE SolverMappingEMToSMMap_Finalise

  !
  !================================================================================================================================
  !

  !>Allocates and initialises an equations matrix to solver matrix map
  SUBROUTINE SolverMappingEMToSMMap_Initialise(equationsMatrixToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: equationsMatrixToSolverMatrixMap !<The equations matrix to solver matrix map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverMappingEMToSMMap_Initialise",err,error,*998)

    IF(ASSOCIATED(equationsMatrixToSolverMatrixMap)) &
      & CALL FlagError("Equations matrix to solver matrix map is already associated.",err,error,*998)

    ALLOCATE(equationsMatrixToSolverMatrixMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations matrix to solver matrix map.",err,error,*999)
    equationsMatrixToSolverMatrixMap%equationsMatrixType=0
    equationsMatrixToSolverMatrixMap%equationsMatrixNumber=0
    equationsMatrixToSolverMatrixMap%solverMatrixNumber=0
    NULLIFY(equationsMatrixToSolverMatrixMap%equationsMatrix)
    NULLIFY(equationsMatrixToSolverMatrixMap%solverMatrix)
    NULLIFY(equationsMatrixToSolverMatrixMap%equationsColToSolverColsMap)
       
    EXITS("SolverMappingEMToSMMap_Initialise")
    RETURN
999 CALL SolverMappingEMToSMMap_Finalise(equationsMatrixToSolverMatrixMap,dummyErr,dummyError,*998)
998 ERRORS("SolverMappingEMToSMMap_Initialise",err,error)    
    EXITS("SolverMappingEMToSMMap_Initialise")    
    RETURN 1
   
  END SUBROUTINE SolverMappingEMToSMMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises a solver matrix interface condition to equations set map and deallocates all memory.
  SUBROUTINE SolverMappingSMICToESMap_Finalise(interfaceConditionToEquationsSetMap,err,error,*)

    !Argument variables
    TYPE(SMInterfaceConditionToEquationsSetMapType), INTENT(INOUT) :: interfaceConditionToEquationsSetMap !<The solver matrix interface condition to equations set map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverMappingSMICToESMap_Finalise",err,error,*999)

    interfaceConditionToEquationsSetMap%interfaceConditionIndex=0
    NULLIFY(interfaceConditionToEquationsSetMap%interfaceCondition)
    interfaceConditionToEquationsSetMap%interfaceMatrixNumber=0
        
    EXITS("SolverMappingSMICToESMap_Finalise")
    RETURN
999 ERRORS("SolverMappingSMICToESMap_Finalise",err,error)    
    EXITS("SolverMappingSMICToESMap_Finalise")    
    RETURN 1
   
  END SUBROUTINE SolverMappingSMICToESMap_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises a equations set to solver matrix interface map.
  SUBROUTINE SolverMappingSMICToESMap_Initialise(interfaceConditionToEquationsSetMap,err,error,*)

    !Argument variables
    TYPE(SMInterfaceConditionToEquationsSetMapType), INTENT(INOUT) :: interfaceConditionToEquationsSetMap !<The equations set to solver map interface to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverMappingSMICToESMap_Initialise",err,error,*999)

    interfaceConditionToEquationsSetMap%interfaceConditionIndex=0
    NULLIFY(interfaceConditionToEquationsSetMap%interfaceCondition)
    interfaceConditionToEquationsSetMap%interfaceMatrixNumber=0
        
    EXITS("SolverMappingSMICToESMap_Initialise")
    RETURN
999 ERRORS("SolverMappingSMICToESMap_Initialise",err,error)    
    EXITS("SolverMappingSMICToESMap_Initialise")    
    RETURN 1
   
  END SUBROUTINE SolverMappingSMICToESMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises a equations matrix to solver matrices map em and deallocates all memory.
  SUBROUTINE SolverMappingEMToSMSMap_Finalise(equationsMatrixToSolverMatricesMap,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixToSolverMatricesMapType), POINTER :: equationsMatrixToSolverMatricesMap !<The equations matrix to solver matrices map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsMatrixIdx
    
    ENTERS("SolverMappingEMToSMSMap_Finalise",err,error,*999)

    IF(ASSOCIATED(equationsMatrixToSolverMatricesMap)) THEN
      IF(ALLOCATED(equationsMatrixToSolverMatricesMap%equationsMatrixToSolverMatrixMaps)) THEN
        DO equationsMatrixIdx=1,SIZE(equationsMatrixToSolverMatricesMap%equationsMatrixToSolverMatrixMaps,1)
          CALL SolverMappingEMToSMMap_Finalise(equationsMatrixToSolverMatricesMap% &
            & equationsMatrixToSolverMatrixMaps(equationsMatrixIdx)%ptr,err,error,*999)        
        ENDDO !equationsMatrixIdx
        DEALLOCATE(equationsMatrixToSolverMatricesMap%equationsMatrixToSolverMatrixMaps)
      ENDIF
      DEALLOCATE(equationsMatrixToSolverMatricesMap)
    ENDIF
    
    EXITS("SolverMappingEMToSMSMap_Finalise")
    RETURN
999 ERRORSEXITS("SolverMappingEMToSMSMap_Finalise",err,error)    
    RETURN 1
   
  END SUBROUTINE SolverMappingEMToSMSMap_Finalise

  !
  !================================================================================================================================
  !

  !>Allocates and initialises an equations matrix to solver matrices map.
  SUBROUTINE SolverMappingEMToSMSMap_Initialise(equationsMatrixToSolverMatricesMap,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixToSolverMatricesMapType), POINTER :: equationsMatrixToSolverMatricesMap !<A pointer to the equations matrices to solver matrices map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverMappingEMToSMSMap_Initialise",err,error,*998)

    IF(ASSOCIATED(equationsMatrixToSolverMatricesMap)) &
      & CALL FlagError("Equations matrix to solver matrices map is already associated.",err,error,*998)

    ALLOCATE(equationsMatrixToSolverMatricesMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocated equations matrix to solver matrices map.",err,error,*999)
    equationsMatrixToSolverMatricesMap%equationsMatrixNumber=0
    equationsMatrixToSolverMatricesMap%numberOfSolverMatrices=0
        
    EXITS("SolverMappingEMToSMSMap_Initialise")
    RETURN
999 CALL SolverMappingEMToSMSMap_Finalise(equationsMatrixToSolverMatricesMap,dummyErr,dummyError,*998)
998 ERRORS("SolverMappingEMToSMSMap_Initialise",err,error)    
    EXITS("SolverMappingEMToSMSMap_Initialise")    
    RETURN 1
   
  END SUBROUTINE SolverMappingEMToSMSMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises a equations matrices to solver matrix map and deallocates all memory.
  SUBROUTINE SolverMappingEMSToSMMap_Finalise(equationsMatricesToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap !<A pointer to the equations matrices to solver matrix map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,variableIdx
    
    ENTERS("SolverMappingEMSToSMMap_Finalise",err,error,*999)

    IF(ASSOCIATED(equationsMatricesToSolverMatrixMap)) THEN
      IF(ALLOCATED(equationsMatricesToSolverMatrixMap%variableTypes)) DEALLOCATE(equationsMatricesToSolverMatrixMap%variableTypes)
      IF(ALLOCATED(equationsMatricesToSolverMatrixMap%variables)) DEALLOCATE(equationsMatricesToSolverMatrixMap%variables)
      IF(ALLOCATED(equationsMatricesToSolverMatrixMap%varDOFToSolverDOFsMaps)) THEN
        DO variableIdx=1,SIZE(equationsMatricesToSolverMatrixMap%varDOFToSolverDOFsMaps,1)
          CALL SolverMappingVDOFToSDOFsMap_Finalise(equationsMatricesToSolverMatrixMap%varDOFToSolverDOFsMaps(variableIdx)%ptr, &
            & err,error,*999)
        ENDDO !variableIdx
        DEALLOCATE(equationsMatricesToSolverMatrixMap%varDOFToSolverDOFsMaps)
      ENDIF
      IF(ALLOCATED(equationsMatricesToSolverMatrixMap%dynamicMatrixToSolverMatrixMaps)) THEN
        DO matrixIdx=1,SIZE(equationsMatricesToSolverMatrixMap%dynamicMatrixToSolverMatrixMaps,1)
          CALL SolverMappingEMToSMMap_Finalise(equationsMatricesToSolverMatrixMap% &
            & dynamicMatrixToSolverMatrixMaps(matrixIdx)%ptr,err,error,*999)        
        ENDDO !matrixIdx
        DEALLOCATE(equationsMatricesToSolverMatrixMap%dynamicMatrixToSolverMatrixMaps)
      ENDIF
      IF(ALLOCATED(equationsMatricesToSolverMatrixMap%linearMatrixToSolverMatrixMaps)) THEN
        DO matrixIdx=1,SIZE(equationsMatricesToSolverMatrixMap%linearMatrixToSolverMatrixMaps,1)
          CALL SolverMappingEMToSMMap_Finalise(equationsMatricesToSolverMatrixMap% &
            & linearMatrixToSolverMatrixMaps(matrixIdx)%ptr,err,error,*999)        
        ENDDO !matrixIdx
        DEALLOCATE(equationsMatricesToSolverMatrixMap%linearMatrixToSolverMatrixMaps)
      ENDIF
      IF(ALLOCATED(equationsMatricesToSolverMatrixMap%jacobianMatrixToSolverMatrixMaps)) THEN
        DO matrixIdx=1,SIZE(equationsMatricesToSolverMatrixMap%jacobianMatrixToSolverMatrixMaps,1)
          CALL SolverMappingJMToSMMap_Finalise(equationsMatricesToSolverMatrixMap% &
            & jacobianMatrixToSolverMatrixMaps(matrixIdx)%ptr,err,error,*999)
        ENDDO !matrixIdx
        DEALLOCATE(equationsMatricesToSolverMatrixMap%jacobianMatrixToSolverMatrixMaps)
      ENDIF
      DEALLOCATE(equationsMatricesToSolverMatrixMap)
    ENDIF

    EXITS("SolverMappingEMSToSMMap_Finalise")
    RETURN
999 ERRORSEXITS("SolverMappingEMSToSMMap_Finalise",err,error)    
    RETURN 1
   
  END SUBROUTINE SolverMappingEMSToSMMap_Finalise

  !
  !================================================================================================================================
  !

  !>Allocate and initialise an equations matrices to solver matrix map.
  SUBROUTINE SolverMappingEMSToSMMap_Initialise(equationsMatricesToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap !<A pointer to the equations matrices to solver matrix map to initialise. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverMappingEMSToSMMap_Initialise",err,error,*998)

    IF(ASSOCIATED(equationsMatricesToSolverMatrixMap)) &
      & CALL FlagError("Equations matrices to solver matrix map is already associated.",err,error,*998)

    ALLOCATE(equationsMatricesToSolverMatrixMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations matrices to solver matrix map.",err,error,*999)
    equationsMatricesToSolverMatrixMap%solverMatrixNumber=0
    equationsMatricesToSolverMatrixMap%numberOfVariables=0
    equationsMatricesToSolverMatrixMap%numberOfDynamicMatrices=0
    equationsMatricesToSolverMatrixMap%numberOfLinearMatrices=0
    equationsMatricesToSolverMatrixMap%numberOfJacobianMatrices=0

    EXITS("SolverMappingEMSToSMMap_Initialise")
    RETURN
999 CALL SolverMappingEMSToSMMap_Finalise(equationsMatricesToSolverMatrixMap,dummyErr,dummyError,*998)
998 ERRORS("SolverMappingEMSToSMMap_Initialise",err,error)    
    EXITS("SolverMappingEMSToSMMap_Initialise")    
    RETURN 1
   
  END SUBROUTINE SolverMappingEMSToSMMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver mapping and deallocates all memory.
  SUBROUTINE SolverMapping_Finalise(solverMapping,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,interfaceConditionIdx,rowIdx,solverMatrixIdx

    ENTERS("SolverMapping_Finalise",err,error,*999)

    IF(ASSOCIATED(solverMapping)) THEN
      IF(ALLOCATED(solverMapping%equationsSets)) DEALLOCATE(solverMapping%equationsSets)        
      IF(ALLOCATED(solverMapping%equationsSetToSolverMatricesMaps)) THEN
        DO equationsSetIdx=1,SIZE(solverMapping%equationsSetToSolverMatricesMaps,1)
          CALL SolverMappingESToSMSMap_Finalise(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr, &
            & err,error,*999)          
        ENDDO !equationsSetIdx
        DEALLOCATE(solverMapping%equationsSetToSolverMatricesMaps)
      ENDIF
      IF(ALLOCATED(solverMapping%interfaceConditions)) DEALLOCATE(solverMapping%interfaceConditions)
      IF(ALLOCATED(solverMapping%interfaceConditionToSolverMatricesMaps)) THEN
        DO interfaceConditionIdx=1,SIZE(solverMapping%interfaceConditionToSolverMatricesMaps,1)
          CALL SolverMappingICToSMSMap_Finalise(solverMapping% &
            & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr,err,error,*999)
        ENDDO !interfaceConditionIdx
        DEALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps)
      ENDIF
      IF(ALLOCATED(solverMapping%solverMatricesToEquationsMaps)) THEN
        DO solverMatrixIdx=1,SIZE(solverMapping%solverMatricesToEquationsMaps,1)
          CALL SolverMappingSMToEquationsMap_Finalise(solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr, &
            & err,error,*999)
        ENDDO !solverMatrixIdx
        DEALLOCATE(solverMapping%solverMatricesToEquationsMaps)
      ENDIF
      IF(ALLOCATED(solverMapping%solverRowToEquationsMaps)) THEN
        DO rowIdx=1,SIZE(solverMapping%solverRowToEquationsMaps,1)
          CALL SolverMappingSRToEQSMap_Finalise(solverMapping%solverRowToEquationsMaps(rowIdx)%ptr,err,error,*999)
        ENDDO !rowIdx
        DEALLOCATE(solverMapping%solverRowToEquationsMaps)
      ENDIF
      CALL DomainMapping_Finalise(solverMapping%rowDOFsMapping,err,error,*999)
      CALL SolverMappingVariables_Finalise(solverMapping%rhsVariablesList,err,error,*999)      
      CALL SolverMapping_CreateValuesCacheFinalise(solverMapping%createValuesCache,err,error,*999)
      DEALLOCATE(solverMapping)
    ENDIF
       
    EXITS("SolverMapping_Finalise")
    RETURN
999 ERRORSEXITS("SolverMapping_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the solver mapping and deallocates all memory.
  SUBROUTINE SolverMapping_Initialise(solverEquations,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to initialise the solver mapping on.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverMapping_Initialise",err,error,*998)

    IF(.NOT.ASSOCIATED(solverEquations))  CALL FlagError("Solver equations is not associated.",err,error,*998)
    IF(ASSOCIATED(solverEquations%solverMapping)) &
      & CALL FlagError("Solver equations solver mapping is already associated.",err,error,*998)
      
    ALLOCATE(solverEquations%solverMapping,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver equations solver mapping.",err,error,*999)
    solverEquations%solverMapping%solverEquations=>solverEquations
    solverEquations%solverMapping%solverMappingFinished=.FALSE.
    solverEquations%solverMapping%numberOfSolverMatrices=1
    solverEquations%solverMapping%numberOfRows=0
    solverEquations%solverMapping%numberOfGlobalRows=0
    solverEquations%solverMapping%numberOfEquationsSets=0
    solverEquations%solverMapping%numberOfInterfaceConditions=0
    NULLIFY(solverEquations%solverMapping%rhsVariablesList)
    NULLIFY(solverEquations%solverMapping%rowDOFsMapping)
    NULLIFY(solverEquations%solverMapping%createValuesCache)
    CALL SolverMapping_CreateValuesCacheInitialise(solverEquations%solverMapping,err,error,*999)
    
    EXITS("SolverMapping_Initialise")
    RETURN
999 CALL SolverMapping_Finalise(solverEquations%solverMapping,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverMapping_Initialise",err,error)
    RETURN 1
  END SUBROUTINE SolverMapping_Initialise

  !
  !================================================================================================================================
  !

  !>Adds an interface condition to a solver mapping
  SUBROUTINE SolverMapping_InterfaceConditionAdd(solverMapping,interfaceCondition,interfaceConditionIndex,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer the solver mapping to add the interface condition to
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to add
    INTEGER(INTG), INTENT(OUT) :: interfaceConditionIndex !<On exit, the index of the interface condition in the solver mapping
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dependentVariableType,equationMatrixIdx,equationsSetIdx,interfaceConditionIdx,interfaceMatrixIdx, &
      & listItem(2),numberOfEquationsSets,numberOfInterfaceMatrices,residualIdx
    LOGICAL :: equationsSetFound,variableFound
    TYPE(EquationsSetType), POINTER :: equationsSet,equationsSet2
    TYPE(FieldVariableType), POINTER :: dependentVariable
    TYPE(InterfaceConditionPtrType), ALLOCATABLE :: newInterfaceConditions(:)
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap
    TYPE(ListType), POINTER :: interfaceIndicesList
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverMapping_InterfaceConditionAdd",err,error,*999)

    interfaceConditionIndex=0

    CALL SolverMapping_AssertNotFinished(solverMapping,err,error,*999)
    CALL InterfaceCondition_AssertIsFinished(interfaceCondition,err,error,*999)    
    NULLIFY(solverEquations)
    CALL SolverMapping_SolverEquationsGet(solverMapping,solverEquations,err,error,*999)
    NULLIFY(interfaceEquations)
    CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
    NULLIFY(interfaceMapping)
    CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
    NULLIFY(createValuesCache)
    CALL SolverMapping_CreateValuesCacheGet(solverMapping,createValuesCache,err,error,*999)
    !Check that the interface variables are already part of an added equations set.
    SELECT CASE(interfaceCondition%method)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
      NULLIFY(interfaceDependent)
      CALL InterfaceCondition_InterfaceDependentGet(interfaceCondition,interfaceDependent,err,error,*999)
      SELECT CASE(interfaceCondition%method)
      CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
        numberOfInterfaceMatrices=interfaceMapping%numberOfInterfaceMatrices
      CASE(INTERFACE_CONDITION_PENALTY_METHOD)
        numberOfInterfaceMatrices=interfaceMapping%numberOfInterfaceMatrices-1
      CASE DEFAULT
        localError="The interface condition method of "//TRIM(NumberToVString(interfaceCondition%method,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      DO interfaceMatrixIdx=1,numberOfInterfaceMatrices
        NULLIFY(interfaceMatrixToVarMap)
        CALL InterfaceMapping_InterfaceMatrixToVarMapGet(interfaceMapping,interfaceMatrixIdx,interfaceMatrixToVarMap, &
          & err,error,*999)
        NULLIFY(equationsSet)
        CALL InterfaceMappingIMToVMap_EquationsSetGet(interfaceMatrixToVarMap,equationsSet,err,error,*999)
        equationsSetFound=.FALSE.
        CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
        DO equationsSetIdx=1,numberOfEquationsSets
          NULLIFY(equationsSet2)
          CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet2,err,error,*999)
          IF(ASSOCIATED(equationsSet,equationsSet2)) THEN
            equationsSetFound=.TRUE.
            EXIT
          ENDIF
        ENDDO !equationsSetIdx
        IF(.NOT.equationsSetFound) THEN
          localError="The equations set for the dependent variable associated with interface "// &
            & "matrix number "//TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))// &
            & " has not been added to the solver equations."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        !See if the variable is in the equations set.
        NULLIFY(dependentVariable)
        CALL InterfaceMappingIMToVMap_VariableGet(interfaceMatrixToVarMap,dependentVariable,err,error,*999)
        CALL FieldVariable_VariableTypeGet(dependentVariable,dependentVariableType,err,error,*999)
        variableFound=.FALSE.
        !Check dynamic variables
        IF(createValuesCache%dynamicVariableType(equationsSetIdx)==dependentVariableType) THEN
          variableFound=.TRUE.
        ELSE
          !Check linear matrices. Just check for solver matrix 1 and the moment
          DO equationMatrixIdx=1,createValuesCache%matrixVariableTypes(0,equationsSetIdx,1)
            IF(createValuesCache%matrixVariableTypes(equationMatrixIdx,equationsSetIdx,1)==dependentVariableType) THEN
              variableFound=.TRUE.
              EXIT
            ENDIF
          ENDDO !equations matrixIdx
          IF(.NOT.variableFound) THEN
            !Check residual variable type
            DO residualIdx=1,SIZE(createValuesCache%residualVariableTypes,2)
              DO equationMatrixIdx=1,createValuesCache%residualVariableTypes(0,residualIdx,equationsSetIdx)
                IF(createValuesCache%residualVariableTypes(equationMatrixIdx,residualIdx,equationsSetIdx)== &
                  & dependentVariableType) THEN
                  variableFound=.TRUE.
                ENDIF
              ENDDO !equationsMatrixIdx
            ENDDO !residualIdx
          ENDIF
        ENDIF
        IF(.NOT.variableFound) THEN
          localError="The dependent variable associated with interface matrix number "// &
            & TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))//" is not mapped to the solver equations."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        !Add in interface condition to equations set (just for solver matrix 1 at the moment)
        listItem(1)=solverMapping%numberOfInterfaceConditions+1
        listItem(2)=interfaceMatrixIdx
        NULLIFY(interfaceIndicesList)
        CALL SolverMappingCVC_InterfaceIndicesListGet(createValuesCache,equationsSetIdx,interfaceIndicesList,err,error,*999)
        CALL List_ItemAdd(interfaceIndicesList,listItem,err,error,*999)
      ENDDO !interfaceMatrixIdx
    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The interface condition method of "//TRIM(NumberToVString(interfaceCondition%method,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    ALLOCATE(newInterfaceConditions(solverMapping%numberOfInterfaceConditions+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new interface conditions.",err,error,*999)
    IF(solverMapping%numberOfInterfaceConditions>0) THEN
      DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
        newInterfaceConditions(interfaceConditionIdx)%ptr=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr
      ENDDO !interfaceConditionIdx
    ENDIF
    newInterfaceConditions(solverMapping%numberOfInterfaceConditions+1)%ptr=>interfaceCondition
    CALL MOVE_ALLOC(newInterfaceConditions,solverMapping%interfaceConditions)
    solverMapping%numberOfInterfaceConditions=solverMapping%numberOfInterfaceConditions+1
    interfaceConditionIndex=solverMapping%numberOfInterfaceConditions
    
!!TODO: SORT OUT LAGRANGE FIELD VARIABLE
    CALL SolverMapping_CreateValuesCacheInterfVarListAdd(solverMapping,1,interfaceConditionIndex,1,err,error,*999)
                
    IF(ALLOCATED(newInterfaceConditions)) DEALLOCATE(newInterfaceConditions)
        
    EXITS("SolverMapping_InterfaceConditionAdd")
    RETURN
999 IF(ALLOCATED(newInterfaceConditions)) DEALLOCATE(newInterfaceConditions)
    ERRORSEXITS("SolverMapping_InterfaceConditionAdd",err,error)    
    RETURN 1
   
  END SUBROUTINE SolverMapping_InterfaceConditionAdd

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of solver matrices for the solver mapping
  SUBROUTINE SolverMapping_NumberOfSolverMatricesSet(solverMapping,numberOfSolverMatrices,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping
    INTEGER(INTG), INTENT(IN) :: numberOfSolverMatrices !<The number of solver matrices for the solver.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,matrixIdx,maximumNumberOfEquationsMatrices
    INTEGER(INTG), ALLOCATABLE :: newMatrixVariableTypes(:,:,:)
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(SolverMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverMapping_NumberOfSolverMatricesSet",err,error,*999)

    CALL SolverMapping_AssertNotFinished(solverMapping,err,error,*999)
    NULLIFY(createValuesCache)
    CALL SolverMapping_CreateValuesCacheGet(solverMapping,createValuesCache,err,error,*999)
    maximumNumberOfEquationsMatrices=1
    DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      NULLIFY(equations)
      CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
      NULLIFY(vectorEquations)
      CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
      NULLIFY(vectorMapping)
      CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
      IF(ASSOCIATED(linearMapping)) THEN
        IF(linearMapping%numberOfLinearMatrices>maximumNumberOfEquationsMatrices) &
          & maximumNumberOfEquationsMatrices=linearMapping%numberOfLinearMatrices
      ENDIF
    ENDDO !equationsSetIdx
    !Check number of matrices to set is valid
    IF(numberOfSolverMatrices<1.OR.numberOfSolverMatrices>maximumNumberOfEquationsMatrices) THEN
      localError="The specified number of solver matrices of "//TRIM(NumberToVString(numberOfSolverMatrices,"*",err,error))// &
        & " is invalid. The number must be >= 1 and <= "//TRIM(NumberToVString(maximumNumberOfEquationsMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    !If we need to reallocate and reset all the create values cache arrays and change the number of matrices
    IF(numberOfSolverMatrices/=solverMapping%numberOfSolverMatrices) THEN
      ALLOCATE(newMatrixVariableTypes(0:FIELD_NUMBER_OF_VARIABLE_TYPES,solverMapping%numberOfEquationsSets, &
        & solverMapping%numberOfSolverMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new matrix variable types.",err,error,*999)
      newMatrixVariableTypes(0:FIELD_NUMBER_OF_VARIABLE_TYPES,1:solverMapping%numberOfEquationsSets, &
        & 1:solverMapping%numberOfSolverMatrices)=solverMapping%createValuesCache%matrixVariableTypes( &
        & 0:FIELD_NUMBER_OF_VARIABLE_TYPES,1:solverMapping%numberOfEquationsSets, &
        & 1:solverMapping%numberOfSolverMatrices)
      IF(numberOfSolverMatrices>solverMapping%numberOfSolverMatrices) THEN
        DO matrixIdx=solverMapping%numberOfSolverMatrices+1,numberOfSolverMatrices
          newMatrixVariableTypes(0,:,matrixIdx)=1
          newMatrixVariableTypes(1,:,matrixIdx)=FIELD_U_VARIABLE_TYPE
          newMatrixVariableTypes(2:FIELD_NUMBER_OF_VARIABLE_TYPES,:,matrixIdx)=0
        ENDDO !matrixIdx
      ENDIF
      CALL MOVE_ALLOC(newMatrixVariableTypes,createValuesCache%matrixVariableTypes)
      solverMapping%numberOfSolverMatrices=numberOfSolverMatrices
    ENDIF
   
    EXITS("SolverMapping_NumberOfSolverMatricesSet")
    RETURN
999 IF(ALLOCATED(newMatrixVariableTypes)) DEALLOCATE(newMatrixVariableTypes)
    ERRORSEXITS("SolverMapping_NumberOfSolverMatricesSet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_NumberOfSolverMatricesSet
  
  !
  !================================================================================================================================
  !

  !>Finalises an interface condition to solver matrices map and deallocates all memory.
  SUBROUTINE SolverMappingICToSMSMap_Finalise(interfaceConditionToSolverMatricesMap,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap !<A pointer to the interface condition to solver matrices map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,equationsSetIdx,interfaceMatrixIdx,solverMatrixIdx
    
    ENTERS("SolverMappingICToSMSMap_Finalise",err,error,*999)

    IF(ASSOCIATED(interfaceConditionToSolverMatricesMap)) THEN
      IF(ALLOCATED(interfaceConditionToSolverMatricesMap%equationsSetToInterfaceConditionMaps)) THEN
        DO equationsSetIdx=1,SIZE(interfaceConditionToSolverMatricesMap%equationsSetToInterfaceConditionMaps,1)
          CALL SolverMappingSMESToICMap_Finalise(interfaceConditionToSolverMatricesMap% &
            & equationsSetToInterfaceConditionMaps(equationsSetIdx),err,error,*999)
        ENDDO !equationsSetIdx
        DEALLOCATE(interfaceConditionToSolverMatricesMap%equationsSetToInterfaceConditionMaps)
      ENDIF
      IF(ALLOCATED(interfaceConditionToSolverMatricesMap%interfaceMatricesToSolverMatrixMaps)) THEN
        DO solverMatrixIdx=1,SIZE(interfaceConditionToSolverMatricesMap%interfaceMatricesToSolverMatrixMaps,1)
          CALL SolverMappingIMSToSMMap_Finalise(interfaceConditionToSolverMatricesMap% &
            & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr,err,error,*999)
        ENDDO !solverMatrixIdx
        DEALLOCATE(interfaceConditionToSolverMatricesMap%interfaceMatricesToSolverMatrixMaps)
      ENDIF
      IF(ALLOCATED(interfaceConditionToSolverMatricesMap%interfaceMatrixToSolverMatricesMaps)) THEN
        DO interfaceMatrixIdx=1,SIZE(interfaceConditionToSolverMatricesMap%interfaceMatrixToSolverMatricesMaps,1)
          CALL SolverMappingIMToSMSMap_Finalise(interfaceConditionToSolverMatricesMap% &
            & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr,err,error,*999)
        ENDDO !interfaceMatrixIdx
        DEALLOCATE(interfaceConditionToSolverMatricesMap%interfaceMatrixToSolverMatricesMaps)
      ENDIF
      IF(ASSOCIATED(interfaceConditionToSolverMatricesMap%interfaceColToSolverRowsMap)) THEN
        DO columnIdx=1,SIZE(interfaceConditionToSolverMatricesMap%interfaceColToSolverRowsMap,1)
          CALL MatrixRowColCoupling_Finalise(interfaceConditionToSolverMatricesMap% &
            & interfaceColToSolverRowsMap(columnIdx),err,error,*999)
        ENDDO !columnIdx
        DEALLOCATE(interfaceConditionToSolverMatricesMap%interfaceColToSolverRowsMap)
      ENDIF
      DEALLOCATE(interfaceConditionToSolverMatricesMap)
    ENDIF
        
    EXITS("SolverMappingICToSMSMap_Finalise")
    RETURN
999 ERRORS("SolverMappingICToSMSMap_Finalise",err,error)    
    EXITS("SolverMappingICToSMSMap_Finalise")    
    RETURN 1
   
  END SUBROUTINE SolverMappingICToSMSMap_Finalise

  !
  !================================================================================================================================
  !

  !>Allocates and initialises an interface condition to solver matrices map.
  SUBROUTINE SolverMappingICToSMSMap_Initialise(interfaceConditionToSolverMatricesMap,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap !<A pointer to the interface condition to solver matrices map to initialise. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverMappingICToSMSMap_Initialise",err,error,*998)

    IF(ASSOCIATED(interfaceConditionToSolverMatricesMap)) &
      & CALL FlagError("Interface condition to solver matrices map is already associated.",err,error,*998)

    ALLOCATE(interfaceConditionToSolverMatricesMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface condition to solver matrices map.",err,error,*999)
    interfaceConditionToSolverMatricesMap%interfaceConditionIndex=0
    NULLIFY(interfaceConditionToSolverMatricesMap%solverMapping)
    NULLIFY(interfaceConditionToSolverMatricesMap%interfaceEquations)
    interfaceConditionToSolverMatricesMap%numberOfEquationsSets=0
    interfaceConditionToSolverMatricesMap%numberOfSolverMatrices=0
    interfaceConditionToSolverMatricesMap%numberOfInterfaceMatrices=0
    NULLIFY(interfaceConditionToSolverMatricesMap%interfaceColToSolverRowsMap)
    
    EXITS("SolverMappingICToSMSMap_Initialise")
    RETURN
999 CALL SolverMappingICToSMSMap_Finalise(interfaceConditionToSolverMatricesMap,dummyErr,dummyError,*998)
998 ERRORS("SolverMappingICToSMSMap_Initialise",err,error)    
    EXITS("SolverMappingICToSMSMap_Initialise")    
    RETURN 1
   
  END SUBROUTINE SolverMappingICToSMSMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises an interface matrix to solver matrix map and deallocates all memory.
  SUBROUTINE SolverMappingIMToSMMap_Finalise(interfaceMatrixToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToSolverMatrixMapType), POINTER :: interfaceMatrixToSolverMatrixMap !<A pointer to the interface matrix to solver matrix map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: rowIdx
    
    ENTERS("SolverMappingIMToSMMap_Finalise",err,error,*999)

    IF(ASSOCIATED(interfaceMatrixToSolverMatrixMap)) THEN
      IF(ASSOCIATED(interfaceMatrixToSolverMatrixMap%interfaceRowToSolverColsMap)) THEN
        DO rowIdx=1,SIZE(interfaceMatrixToSolverMatrixMap%interfaceRowToSolverColsMap,1)
          CALL MatrixRowColCoupling_Finalise(interfaceMatrixToSolverMatrixMap%interfaceRowToSolverColsMap(rowIdx),err,error,*999)
        ENDDO !rowIdx
        DEALLOCATE(interfaceMatrixToSolverMatrixMap%interfaceRowToSolverColsMap)
      ENDIF
      DEALLOCATE(interfaceMatrixToSolverMatrixMap)
    ENDIF
        
    EXITS("SolverMappingIMToSMMap_Finalise")
    RETURN
999 ERRORSEXITS("SolverMappingIMToSMMap_Finalise",err,error)    
    RETURN 1
   
  END SUBROUTINE SolverMappingIMToSMMap_Finalise

  !
  !================================================================================================================================
  !

  !>Allocates and initialises an interface matrix to solver matrix map
  SUBROUTINE SolverMappingIMToSMMap_Initialise(interfaceMatrixToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToSolverMatrixMapType), POINTER :: interfaceMatrixToSolverMatrixMap !<A pointer to the interface matrix to solver matrix map to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("SolverMappingIMToSMMap_Initialise",err,error,*998)

    IF(ASSOCIATED(interfaceMatrixToSolverMatrixMap)) &
      & CALL FlagError("Interface matrix to solver matrix map is already associated.",err,error,*998)

    ALLOCATE(interfaceMatrixToSolverMatrixMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface matrix to solver matrix map.",err,error,*999)
    interfaceMatrixToSolverMatrixMap%interfaceMatrixNumber=0
    interfaceMatrixToSolverMatrixMap%solverMatrixNumber=0
    NULLIFY(interfaceMatrixToSolverMatrixMap%interfaceMatrix)
    NULLIFY(interfaceMatrixToSolverMatrixMap%solverMatrix)
    NULLIFY(interfaceMatrixToSolverMatrixMap%interfaceRowToSolverColsMap)
        
    EXITS("SolverMappingIMToSMMap_Initialise")
    RETURN
999 CALL SolverMappingIMToSmMap_Finalise(interfaceMatrixToSolverMatrixMap,dummyErr,dummyError,*998)
998 ERRORS("SolverMappingIMToSMMap_Initialise",err,error)    
    EXITS("SolverMappingIMToSMMap_Initialise")    
    RETURN 1
   
  END SUBROUTINE SolverMappingIMToSMMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises an solver matrix equations set to interface condition map and deallocates all memory.
  SUBROUTINE SolverMappingSMESToICMap_Finalise(equationsSetToInterfaceConditionMap,err,error,*)

    !Argument variables
    TYPE(SMEquationsSetToInterfaceConditionMapType) :: equationsSetToInterfaceConditionMap !<The solver matrix equations set to interface condition map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverMappingSMESToICMap_Finalise",err,error,*999)

    equationsSetToInterfaceConditionMap%equationsSetIndex=0
    NULLIFY(equationsSetToInterfaceConditionMap%equationsSet)
    equationsSetToInterfaceConditionMap%interfaceMatrixIndex=0
         
    EXITS("SolverMappingSMESToICMap_Finalise")
    RETURN
999 ERRORS("SolverMappingSMESToICMap_Finalise",err,error)    
    EXITS("SolverMappingSMESToICMap_Finalise")    
    RETURN 1
   
  END SUBROUTINE SolverMappingSMESToICMap_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises an solver matrix equations set to interface condition map.
  SUBROUTINE SolverMappingSMESToICMap_Initialise(equationsSetToInterfaceConditionMap,err,error,*)

    !Argument variables
    TYPE(SMEquationsSetToInterfaceConditionMapType) :: equationsSetToInterfaceConditionMap !<The solver matrix equations set to interface condition map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverMappingSMESToICMap_Initialise",err,error,*999)

    equationsSetToInterfaceConditionMap%equationsSetIndex=0
    NULLIFY(equationsSetToInterfaceConditionMap%equationsSet)
    equationsSetToInterfaceConditionMap%interfaceMatrixIndex=0
         
    EXITS("SolverMappingSMESToICMap_Initialise")
    RETURN
999 ERRORS("SolverMappingSMESToICMap_Initialise",err,error)    
    EXITS("SolverMappingSMESToICMap_Initialise")    
    RETURN 1
   
  END SUBROUTINE SolverMappingSMESToICMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises an interface matrix to solver matrices map im and deallocates all memory.
  SUBROUTINE SolverMappingIMToSMSMap_Finalise(interfaceMatrixToSolverMatricesMap,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToSolverMatricesMapType), POINTER :: interfaceMatrixToSolverMatricesMap !<A pointer to the interface matrix to solver matrices maps to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: rowIdx,solverMatrixIdx
    
    ENTERS("SolverMappingIMToSMSMap_Finalise",err,error,*999)

    IF(ASSOCIATED(interfaceMatrixToSolverMatricesMap)) THEN
      IF(ALLOCATED(interfaceMatrixToSolverMatricesMap%interfaceMatrixToSolverMatrixMaps)) THEN
        DO solverMatrixIdx=1,SIZE(interfaceMatrixToSolverMatricesMap%interfaceMatrixToSolverMatrixMaps,1)
          CALL SolverMappingIMToSMMap_Finalise(interfaceMatrixToSolverMatricesMap% &
            & interfaceMatrixToSolverMatrixMaps(solverMatrixIdx)%ptr,err,error,*999)        
        ENDDO !solverMatrixIdx
        DEALLOCATE(interfaceMatrixToSolverMatricesMap%interfaceMatrixToSolverMatrixMaps)
      ENDIF
      IF(ASSOCIATED(interfaceMatrixToSolverMatricesMap%interfaceRowToSolverRowsMap)) THEN
        DO rowIdx=1,SIZE(interfaceMatrixToSolverMatricesMap%interfaceRowToSolverRowsMap,1)
          CALL MatrixRowColCoupling_Finalise(interfaceMatrixToSolverMatricesMap%interfaceRowToSolverRowsMap(rowIdx),err,error,*999)
        ENDDO !rowIdx
        DEALLOCATE(interfaceMatrixToSolverMatricesMap%interfaceRowToSolverRowsMap)
      ENDIF
      DEALLOCATE(interfaceMatrixToSolverMatricesMap)
    ENDIF
    
    EXITS("SolverMappingIMToSMSMap_Finalise")
    RETURN
999 ERRORSEXITS("SolverMappingIMToSMSMap_Finalise",err,error)    
    RETURN 1
   
  END SUBROUTINE SolverMappingIMToSMSMap_Finalise

  !
  !================================================================================================================================
  !

  !>Allocates and initialises an interface matrix to solver matrices map.
  SUBROUTINE SolverMappingIMToSMSMap_Initialise(interfaceMatrixToSolverMatricesMap,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToSolverMatricesMapType), POINTER :: interfaceMatrixToSolverMatricesMap !<A pointer to the interface matrix to solver matrices maps to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverMappingIMToSMSMap_Initialise",err,error,*998)

    IF(ASSOCIATED(interfaceMatrixToSolverMatricesMap)) &
      & CALL FlagError("Interface matrix to solver matrices map is already associated.",err,error,*998)

    ALLOCATE(interfaceMatrixToSolverMatricesMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate the interface matrix to solver matrices map.",err,error,*999)
    interfaceMatrixToSolverMatricesMap%interfaceMatrixNumber=0
    interfaceMatrixToSolverMatricesMap%numberOfSolverMatrices=0
    NULLIFY(interfaceMatrixToSolverMatricesMap%interfaceRowToSolverRowsMap)
    
    EXITS("SolverMappingIMToSMSMap_Initialise")
    RETURN
999 CALL SolverMappingIMToSMSMap_Finalise(interfaceMatrixToSolverMatricesMap,dummyErr,dummyError,*998)
998 ERRORS("SolverMappingIMToSMSMap_Initialise",err,error)    
    EXITS("SolverMappingIMToSMSMap_Initialise")    
    RETURN 1
   
  END SUBROUTINE SolverMappingIMToSMSMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises an interface matrix to solver matrices map and deallocates all memory.
  SUBROUTINE SolverMappingIMSToSMMap_Finalise(interfaceMatricesToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMap !<A pointer to the interface matrices to solver matrix map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,interfaceMatrixIdx
    
    ENTERS("SolverMappingIMSToSMMap_Finalise",err,error,*999)

    IF(ASSOCIATED(interfaceMatricesToSolverMatrixMap)) THEN
       CALL SolverMappingVDOFToSDOFsMap_Finalise(interfaceMatricesToSolverMatrixMap%lagrangeVarDOFToSolverDOFsMap,err,error,*999)
       IF(ALLOCATED(interfaceMatricesToSolverMatrixMap%dependentVariableTypes)) &
         & DEALLOCATE(interfaceMatricesToSolverMatrixMap%dependentVariableTypes)
       IF(ALLOCATED(interfaceMatricesToSolverMatrixMap%dependentVariables)) &
         & DEALLOCATE(interfaceMatricesToSolverMatrixMap%dependentVariables)
       IF(ALLOCATED(interfaceMatricesToSolverMatrixMap%dependentVarDOFToSolverDOFsMaps)) THEN
        DO interfaceMatrixIdx=1,SIZE(interfaceMatricesToSolverMatrixMap%dependentVarDOFToSolverDOFsMaps,1)
          CALL SolverMappingVDOFToSDOFsMap_Finalise(interfaceMatricesToSolverMatrixMap% &
            & dependentVarDOFToSolverDOFsMaps(interfaceMatrixIdx)%ptr,err,error,*999)
        ENDDO !interfaceMatrixIdx
        DEALLOCATE(interfaceMatricesToSolverMatrixMap%dependentVarDOFToSolverDOFsMaps)
      ENDIF
      IF(ALLOCATED(interfaceMatricesToSolverMatrixMap%interfaceMatrixToSolverMatrixMaps)) THEN
        DO interfaceMatrixIdx=1,SIZE(interfaceMatricesToSolverMatrixMap%interfaceMatrixToSolverMatrixMaps,1)
          CALL SolverMappingIMToSMMap_Finalise(interfaceMatricesToSolverMatrixMap% &
            interfaceMatrixToSolverMatrixMaps(interfaceMatrixIdx)%ptr,err,error,*999)
        ENDDO !interfaceMatrixIdx
        DEALLOCATE(interfaceMatricesToSolverMatrixMap%interfaceMatrixToSolverMatrixMaps)
      ENDIF
      IF(ASSOCIATED(interfaceMatricesToSolverMatrixMap%interfaceColToSolverColsMap)) THEN
        DO columnIdx=1,SIZE(interfaceMatricesToSolverMatrixMap%interfaceColToSolverColsMap,1)
          CALL MatrixRowColCoupling_Finalise(interfaceMatricesToSolverMatrixMap% &
            & interfaceColToSolverColsMap(columnIdx),err,error,*999)
        ENDDO !columnIdx
        DEALLOCATE(interfaceMatricesToSolverMatrixMap%interfaceColToSolverColsMap)
      ENDIF
      DEALLOCATE(interfaceMatricesToSolverMatrixMap)
    ENDIF
   
    EXITS("SolverMappingIMSToSMMap_Finalise")
    RETURN
999 ERRORSEXITS("SolverMappingIMSToSMMap_Finalise",err,error)    
    RETURN 1
   
  END SUBROUTINE SolverMappingIMSToSMMap_Finalise

  !
  !================================================================================================================================
  !

  !>Allocates and initialises an interface matrices to solver matrix map.
  SUBROUTINE SolverMappingIMSToSMMap_Initialise(interfaceMatricesToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMap !<A pointer to the interface matrices to solver matrix map to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverMappingIMSToSMMap_Initialise",err,error,*998)

    IF(ASSOCIATED(interfaceMatricesToSolverMatrixMap)) &
      & CALL FlagError("Interface matrices to solver matrix map is already associated.",err,error,*998)

    ALLOCATE(interfaceMatricesToSolverMatrixMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface matrices to solver matrix map.",err,error,*999)
    interfaceMatricesToSolverMatrixMap%solverMatrixNumber=0
    interfaceMatricesToSolverMatrixMap%lagrangeVariableType=0
    NULLIFY(interfaceMatricesToSolverMatrixMap%lagrangeVariable)
    NULLIFY(interfaceMatricesToSolverMatrixMap%lagrangeVarDOFToSolverDOFsMap)
    interfaceMatricesToSolverMatrixMap%numberOfDependentVariables=0
    interfaceMatricesToSolverMatrixMap%numberOfInterfaceMatrices=0
    NULLIFY(interfaceMatricesToSolverMatrixMap%interfaceColToSolverColsMap)
        
    EXITS("SolverMappingIMSToSMMap_Initialise")
    RETURN
999 CALL SolverMappingIMSToSMMap_Finalise(interfaceMatricesToSolverMatrixMap,dummyErr,dummyError,*998)
998 ERRORS("SolverMappingIMSToSMMap_Initialise",err,error)    
    EXITS("SolverMappingIMSToSMMap_Initialise")    
    RETURN 1
   
  END SUBROUTINE SolverMappingIMSToSMMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises a Jacobian matrix to solver matrix map and deallocates all memory.
  SUBROUTINE SolverMappingJMToSMMap_Finalise(jacobianMatrixToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixToSolverMatrixMapType), POINTER :: jacobianMatrixToSolverMatrixMap !<A pointer to the Jacobian matrix to solver matrix map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx
    
    ENTERS("SolverMappingJMToSMMap_Finalise",err,error,*999)

    IF(ASSOCIATED(jacobianMatrixToSolverMatrixMap)) THEN
      IF(ASSOCIATED(jacobianMatrixToSolverMatrixMap%jacobianColToSolverColsMap)) THEN
        DO columnIdx=1,SIZE(jacobianMatrixToSolverMatrixMap%jacobianColToSolverColsMap,1)
          CALL MatrixRowColCoupling_Finalise(jacobianMatrixToSolverMatrixMap%jacobianColToSolverColsMap(columnIdx),err,error,*999)
        ENDDO !columnIdx
        DEALLOCATE(jacobianMatrixToSolverMatrixMap%jacobianColToSolverColsMap)
      ENDIF
      DEALLOCATE(jacobianMatrixToSolverMatrixMap)
    ENDIF
        
    EXITS("SolverMappingJMToSMMap_Finalise")
    RETURN
999 ERRORSEXITS("SolverMappingJMToSMMap_Finalise",err,error)    
    RETURN 1
   
  END SUBROUTINE SolverMappingJMToSMMap_Finalise

  !
  !================================================================================================================================
  !

  !>Allocates an initialises a Jacobian matrix to solver matrix map
  SUBROUTINE SolverMappingJMToSMMap_Initialise(jacobianMatrixToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixToSolverMatrixMapType), POINTER :: jacobianMatrixToSolverMatrixMap !<A pointer to the Jacobian matrix to solver matrix map to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverMappingJMToSMMap_Initialise",err,error,*998)

    IF(ASSOCIATED(jacobianMatrixToSolverMatrixMap)) &
      & CALL FlagError("Jacobian matrix to solver matrix map is already associated.",err,error,*998)

    ALLOCATE(jacobianMatrixToSolverMatrixMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate Jacobian matrix to solver matrix map.",err,error,*999)
    jacobianMatrixToSolverMatrixMap%jacobianMatrixNumber=0
    jacobianMatrixToSolverMatrixMap%solverMatrixNumber=0
    NULLIFY(jacobianMatrixToSolverMatrixMap%jacobianMatrix)
    NULLIFY(jacobianMatrixToSolverMatrixMap%solverMatrix)
    NULLIFY(jacobianMatrixToSolverMatrixMap%jacobianColToSolverColsMap)
         
    EXITS("SolverMappingJMToSMMap_Initialise")
    RETURN
999 CALL SolverMappingJMToSMMap_Finalise(jacobianMatrixToSolverMatrixMap,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverMappingJMToSMMap_Initialise",err,error)    
    RETURN 1
   
  END SUBROUTINE SolverMappingJMToSMMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver column to dynamic equations map and deallocates all memory.
  SUBROUTINE SolverMappingSColToDEQSMap_Finalise(solverColToDynamicEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverColToDynamicEquationsMapType), POINTER :: solverColToDynamicEquationsMap !<A pointer to the solver column to dynamic equations map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverMappingSColToDEQSMap_Finalise",err,error,*999)

    IF(ASSOCIATED(solverColToDynamicEquationsMap)) THEN
      solverColToDynamicEquationsMap%numberOfDynamicMatrices=0
      IF(ALLOCATED(solverColToDynamicEquationsMap%equationsMatrixNumbers)) &
        & DEALLOCATE(solverColToDynamicEquationsMap%equationsMatrixNumbers)
      IF(ALLOCATED(solverColToDynamicEquationsMap%equationsColumnNumbers)) &
        & DEALLOCATE(solverColToDynamicEquationsMap%equationsColumnNumbers)
      IF(ALLOCATED(solverColToDynamicEquationsMap%couplingCoefficients)) &
        & DEALLOCATE(solverColToDynamicEquationsMap%couplingCoefficients)
      DEALLOCATE(solverColToDynamicEquationsMap)
    ENDIF
       
    EXITS("SolverMappingSColToDEQSMap_Finalise")
    RETURN
999 ERRORS("SolverMappingSColToDEQSMap_Finalise",err,error)
    EXITS("SolverMappingSColToDEQSMap_Finalise")
    RETURN 1
    
  END SUBROUTINE SolverMappingSColToDEQSMap_Finalise

  !
  !================================================================================================================================
  !

  !>Allocates and initialises the solver column to dynamic equations mapping.
  SUBROUTINE SolverMappingSColToDEQSMap_Initialise(solverColToDynamicEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverColToDynamicEquationsMapType), POINTER :: solverColToDynamicEquationsMap !<A pointer to the solver column to dynamic equations map to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("SolverMappingSColToDEQSMap_Initialise",err,error,*998)

    IF(ASSOCIATED(solverColToDynamicEquationsMap)) &
      & CALL FlagError("The solver column to dynamic equations map is already associated.",err,error,*998)

    ALLOCATE(solverColToDynamicEquationsMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate the solver column to dyanmic equations map.",err,error,*999)
    solverColToDynamicEquationsMap%numberOfDynamicMatrices=0
     
    EXITS("SolverMappingSColToDEQSMap_Initialise")
    RETURN
999 CALL SolverMappingSColToDEQSMap_Finalise(solverColToDynamicEquationsMap,dummyErr,dummyError,*998)
998 ERRORS("SolverMappingSColToDEQSMap_Initialise",err,error)
    EXITS("SolverMappingSColToDEQSMap_Initialise")
    RETURN 1
    
  END SUBROUTINE SolverMappingSColToDEQSMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver column to linear equations map and deallocates all memory.
  SUBROUTINE SolverMappingSColToLEQSMap_Finalise(solverColToLinearEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverColToLinearEquationsMapType), POINTER :: solverColToLinearEquationsMap !<A pointer to the solver column to linear equations map to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverMappingSColToLEQSMap_Finalise",err,error,*999)

    IF(ASSOCIATED(solverColToLinearEquationsMap)) THEN
      IF(ALLOCATED(solverColToLinearEquationsMap%equationsMatrixNumbers)) &
        & DEALLOCATE(solverColToLinearEquationsMap%equationsMatrixNumbers)
      IF(ALLOCATED(solverColToLinearEquationsMap%equationsColumnNumbers)) &
        & DEALLOCATE(solverColToLinearEquationsMap%equationsColumnNumbers)
      IF(ALLOCATED(solverColToLinearEquationsMap%couplingCoefficients)) &
        & DEALLOCATE(solverColToLinearEquationsMap%couplingCoefficients)
      DEALLOCATE(solverColToLinearEquationsMap)
    ENDIF
    
    EXITS("SolverMappingSColToLEQSMap_Finalise")
    RETURN
999 ERRORS("SolverMappingSColToLEQSMap_Finalise",err,error)
    EXITS("SolverMappingSColToLEQSMap_Finalise")
    RETURN 1
    
  END SUBROUTINE SolverMappingSColToLEQSMap_Finalise

  !
  !================================================================================================================================
  !

  !>Allocates and initialises the solver column to linear equations mapping.
  SUBROUTINE SolverMappingSColToLEQSMap_Initialise(solverColToLinearEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverColToLinearEquationsMapType), POINTER :: solverColToLinearEquationsMap !<A pointer to the solver column to linear equations map to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverMappingSColToLEQSMap_Initialise",err,error,*998)

    IF(ASSOCIATED(solverColToLinearEquationsMap)) &
      & CALL FlagError("The solver column to linear equations map is already associated.",err,error,*998)

    ALLOCATE(solverColToLinearEquationsMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate the solver column to linear equations map.",err,error,*998)
    solverColToLinearEquationsMap%numberOfLinearMatrices=0
    
    EXITS("SolverMappingSCToSEquationsMap_Initialise")
    RETURN
999 CALL  SolverMappingSColToLEQSMap_Finalise(solverColToLinearEquationsMap,dummyErr,dummyError,*998)
998 ERRORS("SolverMappingSColToLEQSMap_Initialise",err,error)
    EXITS("SolverMappingSColToLEQSMap_Initialise")
    RETURN 1
    
  END SUBROUTINE SolverMappingSColToLEQSMap_Initialise

   !
  !================================================================================================================================
  !

  !>Finalises the solver column to nonlinear equations map and deallocates all memory.
  SUBROUTINE SolverMappingSColToNLEQSMap_Finalise(solverColToNonlinearEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverColToNonlinearEquationsMapType), POINTER :: solverColToNonlinearEquationsMap !<A pointer to the solver column to nonlinear equations map to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverMappingSColToNLEQSMap_Finalise",err,error,*999)

    IF(ASSOCIATED(solverColToNonlinearEquationsMap)) THEN
      IF(ALLOCATED(solverColToNonlinearEquationsMap%jacobianMatrixNumbers)) &
        & DEALLOCATE(solverColToNonlinearEquationsMap%jacobianMatrixNumbers)
      IF(ALLOCATED(solverColToNonlinearEquationsMap%jacobianColumnNumbers)) &
        & DEALLOCATE(solverColToNonlinearEquationsMap%jacobianColumnNumbers)
      IF(ALLOCATED(solverColToNonlinearEquationsMap%couplingCoefficients)) &
        & DEALLOCATE(solverColToNonlinearEquationsMap%couplingCoefficients)
      DEALLOCATE(solverColToNonlinearEquationsMap)
    ENDIF
    
    EXITS("SolverMappingSColToNLEQSMap_Finalise")
    RETURN
999 ERRORS("SolverMappingSColToNLEQSMap_Finalise",err,error)
    EXITS("SolverMappingSColToNLEQSMap_Finalise")
    RETURN 1
    
  END SUBROUTINE SolverMappingSColToNLEQSMap_Finalise

  !
  !================================================================================================================================
  !

  !>Allocates and initialises the solver column to nonlinear equations mapping.
  SUBROUTINE SolverMappingSColToNLEQSMap_Initialise(solverColToNonlinearEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverColToNonlinearEquationsMapType), POINTER :: solverColToNonlinearEquationsMap !<A pointer to the solver column to nonlinear equations map to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverMappingSColToNLEQSMap_Initialise",err,error,*998)

    IF(ASSOCIATED(solverColToNonlinearEquationsMap)) &
      & CALL FlagError("The solver column to nonlinear equations map is already associated.",err,error,*998)

    ALLOCATE(solverColToNonlinearEquationsMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate the solver column to nonlinear equations map.",err,error,*998)
    solverColToNonlinearEquationsMap%numberOfJacobianMatrices=0
    
    EXITS("SolverMappingSColToNLEQSMap_Initialise")
    RETURN
999 CALL  SolverMappingSColToNLEQSMap_Finalise(solverColToNonlinearEquationsMap,dummyErr,dummyError,*998)
998 ERRORS("SolverMappingSColToNLEQSMap_Initialise",err,error)
    EXITS("SolverMappingSColToNLEQSMap_Initialise")
    RETURN 1
    
  END SUBROUTINE SolverMappingSColToNLEQSMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver matrix to equations set map and deallocates all memory.
  SUBROUTINE SolverMappingSMToESMap_Finalise(solverMatrixToEquationsSetMap,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsSetMapType), POINTER :: solverMatrixToEquationsSetMap !<A pointer to the solver matrix to equations set map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx
    
    ENTERS("SolverMappingSMToESMap_Finalise",err,error,*999)

    IF(ASSOCIATED(solverMatrixToEquationsSetMap)) THEN
      IF(ALLOCATED(solverMatrixToEquationsSetMap%solverColToDynamicEquationsMaps)) THEN
        DO columnIdx=1,SIZE(solverMatrixToEquationsSetMap%solverColToDynamicEquationsMaps,1)
          CALL SolverMappingSColToDEQSMap_Finalise(solverMatrixToEquationsSetMap% &
            & solverColToDynamicEquationsMaps(columnIdx)%ptr,err,error,*999)
        ENDDO !columnIdx
        DEALLOCATE(solverMatrixToEquationsSetMap%solverColToDynamicEquationsMaps)
      ENDIF
      IF(ALLOCATED(solverMatrixToEquationsSetMap%solverColToLinearEquationsMaps)) THEN
        DO columnIdx=1,SIZE(solverMatrixToEquationsSetMap%solverColToLinearEquationsMaps,1)
          CALL SolverMappingSColToLEQSMap_Finalise(solverMatrixToEquationsSetMap% &
            & solverColToLinearEquationsMaps(columnIdx)%ptr,err,error,*999)
        ENDDO !columnIdx
        DEALLOCATE(solverMatrixToEquationsSetMap%solverColToLinearEquationsMaps)
      ENDIF
      IF(ALLOCATED(solverMatrixToEquationsSetMap%solverColToNonlinearEquationsMaps)) THEN
        DO columnIdx=1,SIZE(solverMatrixToEquationsSetMap%solverColToNonlinearEquationsMaps,1)
          CALL SolverMappingSColToNLEQSMap_Finalise(solverMatrixToEquationsSetMap% &
            & solverColToNonlinearEquationsMaps(columnIdx)%ptr,err,error,*999)
        ENDDO !columnIdx
        DEALLOCATE(solverMatrixToEquationsSetMap%solverColToNonlinearEquationsMaps)
      ENDIF
      DEALLOCATE(solverMatrixToEquationsSetMap)
    ENDIF
       
    EXITS("SolverMappingSMToESMap_Finalise")
    RETURN
999 ERRORSEXITS("SolverMappingSMToESMap_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToESMap_Finalise

  !
  !================================================================================================================================
  !

  !>Allocates and initialises the solver matrix to equations set mapping and deallocates all memory.
  SUBROUTINE SolverMappingSMToESMap_Initialise(solverMatrixToEquationsSetMap,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsSetMapType), POINTER :: solverMatrixToEquationsSetMap !<A pointer to the solver matrix to equations set map to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverMappingSMToESMap_Initialise",err,error,*998)

    IF(ASSOCIATED(solverMatrixToEquationsSetMap)) &
      & CALL FlagError("Solver matrix to equations set map is already associated.",err,error,*998)

    ALLOCATE(solverMatrixToEquationsSetMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver matrix to equations set map.",err,error,*999)    
    NULLIFY(solverMatrixToEquationsSetMap%equations)
    solverMatrixToEquationsSetMap%haveDynamic=.FALSE.
    solverMatrixToEquationsSetMap%haveLinear=.FALSE.
    solverMatrixToEquationsSetMap%haveNonlinear=.FALSE.
    
    EXITS("SolverMappingSMToESMap_Initialise")
    RETURN
999 CALL SolverMappingSMToESMap_Finalise(solverMatrixtoEquationsSetMap,dummyErr,dummyError,*998)    
998 ERRORS("SolverMappingSMToESMap_Initialise",err,error)
    EXITS("SolverMappingSMToESMap_Initialise")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToESMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver matrix to equations map and deallocates all memory.
  SUBROUTINE SolverMappingSMToEquationsMap_Finalise(solverMatrixToEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap !<The solver matrix to equations map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dofIdx,equationsSetIdx,interfaceConditionIdx
    
    ENTERS("SolverMappingSMToEquationsMap_Finalise",err,error,*999)

    IF(ASSOCIATED(solverMatrixToEquationsMap)) THEN
      CALL SolverMappingVariables_Finalise(solverMatrixToEquationsMap%variablesList,err,error,*999)
      IF(ALLOCATED(solverMatrixToEquationsMap%solverMatrixToEquationsSetMaps)) THEN
        DO equationsSetIdx=1,SIZE(solverMatrixToEquationsMap%solverMatrixToEquationsSetMaps,1)
          CALL SolverMappingSMToESMap_Finalise(solverMatrixToEquationsMap%solverMatrixToEquationsSetMaps(equationsSetIdx)%ptr, &
            & err,error,*999)
        ENDDO !equationsSetIdx
        DEALLOCATE(solverMatrixToEquationsMap%solverMatrixToEquationsSetMaps)
      ENDIF
      IF(ALLOCATED(solverMatrixToEquationsMap%solverMatrixToInterfaceConditionMaps)) THEN
        DO interfaceConditionIdx=1,SIZE(solverMatrixToEquationsMap%solverMatrixToInterfaceConditionMaps,1)
          CALL SolverMappingSMToICMap_Finalise(solverMatrixToEquationsMap% &
            & solverMatrixToInterfaceConditionMaps(interfaceConditionIdx)%ptr,err,error,*999)
        ENDDO !interfaceConditionIdx
        DEALLOCATE(solverMatrixToEquationsMap%solverMatrixToInterfaceConditionMaps)
      ENDIF
      IF(ALLOCATED(solverMatrixToEquationsMap%solverDOFToVariableDOFsMaps)) THEN
        DO dofIdx=1,SIZE(solverMatrixToEquationsMap%solverDOFToVariableDOFsMaps,1)
          CALL SolverMappingSDOFToVDOFsMap_Finalise(solverMatrixToEquationsMap%solverDOFToVariableDOFsMaps(dofIdx)%ptr, &
            & err,error,*999)
        ENDDO !dofIdx
        DEALLOCATE(solverMatrixToEquationsMap%solverDOFToVariableDOFsMaps)
      ENDIF
      CALL DomainMapping_Finalise(solverMatrixToEquationsMap%columnDOFSMapping,err,error,*999)
      DEALLOCATE(solverMatrixToEquationsMap)
    ENDIF
    
    EXITS("SolverMappingSMToEquationsMap_Finalise")
    RETURN
999 ERRORSEXITS("SolverMappingSMToEquationsMap_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToEquationsMap_Finalise

  !
  !================================================================================================================================
  !

  !>Allocates and initialises the solver matrix to equations mapping .
  SUBROUTINE SolverMappingSMToEquationsMap_Initialise(solverMatrixToEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap !<A pointer to the solver matrix to equations map to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverMappingSMToEquationsMap_Initialise",err,error,*998)

    IF(ASSOCIATED(solverMatrixToEquationsMap)) &
      & CALL FlagError("Solver matrix to equations map is already associated.",err,error,*998)

    ALLOCATE(solverMatrixToEquationsMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver matrix to equations map.",err,error,*999)
    NULLIFY(solverMatrixToEquationsMap%solverMapping)
    solverMatrixToEquationsMap%solverMatrixNumber=0
    NULLIFY(solverMatrixToEquationsMap%solverMatrix)
    NULLIFY(solverMatrixToEquationsMap%variablesList)
    solverMatrixToEquationsMap%numberOfColumns=0
    solverMatrixToEquationsMap%numberOfDOFs=0
    solverMatrixToEquationsMap%totalNumberOfDOFs=0
    solverMatrixToEquationsMap%numberOfGlobalDofs=0
    solverMatrixToEquationsMap%numberOfEquationsSets=0
    solverMatrixToEquationsMap%numberOfInterfaceConditions=0
    NULLIFY(solverMatrixToEquationsMap%columnDOFSMapping)
    
    EXITS("SolverMappingSMToEquationsMap_Initialise")
    RETURN
999 CALL SolverMappingSMToEquationsMap_Finalise(solverMatrixToEquationsMap,dummyErr,dummyError,*998)
998 ERRORS("SolverMappingSMToEquationsMap_Initialise",err,error)
    EXITS("SolverMappingSMToEquationsMap_Initialise")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToEquationsMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver matrix to interface condition map and deallocates all memory.
  SUBROUTINE SolverMappingSMToICMap_Finalise(solverMatrixToInterfaceConditionMap,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToInterfaceConditionMapType), POINTER :: solverMatrixToInterfaceConditionMap !<A pointer to the solver matrix to interface condition map to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx
    
    ENTERS("SolverMappingSMToICMap_Finalise",err,error,*999)

    IF(ASSOCIATED(solverMatrixToInterfaceConditionMap)) THEN
      IF(ALLOCATED(solverMatrixToInterfaceConditionMap%solverColToInterfaceEquationsMaps)) THEN
        DO columnIdx=1,SIZE(solverMatrixToInterfaceConditionMap%solverColToInterfaceEquationsMaps,1)
          CALL SolverMappingSColToIEQSMap_Finalise(solverMatrixToInterfaceConditionMap% &
            & solverColToInterfaceEquationsMaps(columnIdx)%ptr,err,error,*999)
        ENDDO !columnIdx
        DEALLOCATE(solverMatrixToInterfaceConditionMap%solverColToInterfaceEquationsMaps)
      ENDIF
      DEALLOCATE(solverMatrixToInterfaceConditionMap)
    ENDIF
     
    EXITS("SolverMappingSMToICMap_Finalise")
    RETURN
999 ERRORSEXITS("SolverMappingSMToICMap_Finalise",err,error)
    RETURN 1
  END SUBROUTINE SolverMappingSMToICMap_Finalise

  !
  !================================================================================================================================
  !

  !>Allocates and initialises the solver matrix to interface condition mapping.
  SUBROUTINE SolverMappingSMToICMap_Initialise(solverMatrixToInterfaceConditionMap,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToInterfaceConditionMapType), POINTER :: solverMatrixToInterfaceConditionMap !<A pointer to the solver matrix to interface condition map to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverMappingSMToICMap_Initialise",err,error,*998)

    IF(ASSOCIATED(solverMatrixToInterfaceConditionMap)) &
      & CALL FlagError("The solver matrix to interface condition map is already associated.",err,error,*998)

    ALLOCATE(solverMatrixToInterfaceConditionMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate the solver matrix to interface condition map.",err,error,*999)
    NULLIFY(solverMatrixToInterfaceConditionMap%interfaceEquations)
    
    EXITS("SolverMappingSMToICMap_Initialise")
    RETURN
999 CALL SolverMappingSMToICMap_Finalise(solverMatrixToInterfaceConditionMap,dummyErr,dummyError,*998)
998 ERRORS("SolverMappingSMToICMap_Initialise",err,error)
    EXITS("SolverMappingSMToICMap_Initialise")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToICMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver column to interface equations map and deallocates all memory.
  SUBROUTINE SolverMappingSColToIEQSMap_Finalise(solverColToInterfaceEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverColToInterfaceEquationsMapType), POINTER :: solverColToInterfaceEquationsMap !<A pointer to the solver column to interface equatiosn map to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverMappingSColToIEQSMap_Finalise",err,error,*999)

    IF(ASSOCIATED(solverColToInterfaceEquationsMap)) THEN
      IF(ALLOCATED(solverColToInterfaceEquationsMap%interfaceMatrixNumbers))  &
        & DEALLOCATE(solverColToInterfaceEquationsMap%interfaceMatrixNumbers)
      IF(ALLOCATED(solverColToInterfaceEquationsMap%interfaceColumnNumbers))  &
        & DEALLOCATE(solverColToInterfaceEquationsMap%interfaceColumnNumbers)
      IF(ALLOCATED(solverColToInterfaceEquationsMap%couplingCoefficients))  &
        & DEALLOCATE(solverColToInterfaceEquationsMap%couplingCoefficients)
      DEALLOCATE(solverColToInterfaceEquationsMap)
    ENDIF
      
    EXITS("SolverMappingSColToIEQSMap_Finalise")
    RETURN
999 ERRORS("SolverMappingSColToIEQSMap_Finalise",err,error)
    EXITS("SolverMappingSColToIEQSMap_Finalise")
    RETURN 1
    
  END SUBROUTINE SolverMappingSColToIEQSMap_Finalise

  !
  !================================================================================================================================
  !

  !>Allocates and initialises the solver column to interface equations mapping.
  SUBROUTINE SolverMappingSColToIEQSMap_Initialise(solverColToInterfaceEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverColToInterfaceEquationsMapType), POINTER :: solverColToInterfaceEquationsMap !<A pointer to the solver column to interface equations map to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverMappingSColToIEQSMap_Initialise",err,error,*998)

    IF(ASSOCIATED(solverColToInterfaceEquationsMap)) &
      & CALL FlagError("The solver column to interface equations map is already associated.",err,error,*998)

    ALLOCATE(solverColToInterfaceEquationsMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocated the solver column to interface equayions map.",err,error,*999)
    solverColToInterfaceEquationsMap%numberOfInterfaceMatrices=0
    
    EXITS("SolverMappingSColToIEQSMap_Initialise")
    RETURN
999 CALL SolverMappingSColToIEQSMap_Finalise(solverColToInterfaceEquationsMap,dummyErr,dummyError,*998)
998 ERRORS("SolverMappingSColToIEQSMap_Initialise",err,error)
    EXITS("SolverMappingSColToIEQSMap_Initialise")
    RETURN 1
    
  END SUBROUTINE SolverMappingSColToIEQSMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver DOF to variable DOFs mapping and deallocates all memory.
  SUBROUTINE SolverMappingSDOFToVDOFsMap_Finalise(solverDOFToVariableDOFsMap,err,error,*)

    !Argument variables
    TYPE(SolverDOFToVariableDOFsMapType), POINTER :: solverDOFToVariableDOFsMap !<A pointer to the solver DOF to variable DOFs map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverMappingSDOFToVDOFsMap_Finalise",err,error,*999)

    IF(ASSOCIATED(solverDOFToVariableDOFsMap)) THEN
      IF(ALLOCATED(solverDOFToVariableDOFsMap%equationTypes)) DEALLOCATE(solverDOFToVariableDOFsMap%equationTypes)
      IF(ALLOCATED(solverDOFToVariableDOFsMap%equationIndices)) DEALLOCATE(solverDOFToVariableDOFsMap%equationIndices)
      IF(ALLOCATED(solverDOFToVariableDOFsMap%variable)) DEALLOCATE(solverDOFToVariableDOFsMap%variable)
      IF(ALLOCATED(solverDOFToVariableDOFsMap%variableDOF)) DEALLOCATE(solverDOFToVariableDOFsMap%variableDOF)
      IF(ALLOCATED(solverDOFToVariableDOFsMap%variableCoefficient)) DEALLOCATE(solverDOFToVariableDOFsMap%variableCoefficient)
      IF(ALLOCATED(solverDOFToVariableDOFsMap%additiveConstant)) DEALLOCATE(solverDOFToVariableDOFsMap%additiveConstant)
      DEALLOCATE(solverDOFToVariableDOFsMap)
    ENDIF
    
    EXITS("SolverMappingSDOFToVDOFsMap_Finalise")
    RETURN
999 ERRORS("SolverMappingSDOFToVDOFsMap_Finalise",err,error)
    EXITS("SolverMappingSDOFToVDOFsMap_Finalise")
    RETURN 1
    
  END SUBROUTINE SolverMappingSDOFToVDOFsMap_Finalise

  !
  !================================================================================================================================
  !

  !>Allocates and initialises the solver DOF to variable DOFs mapping.
  SUBROUTINE SolverMappingSDOFToVDOFsMap_Initialise(solverDOFToVariableDOFsMap,err,error,*)

    !Argument variables
    TYPE(SolverDOFToVariableDOFsMapType), POINTER :: solverDOFToVariableDOFsMap !<A pointer to the solver DOF to variable DOFs map to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverMappingSDOFToVDOFsMap_Initialise",err,error,*998)

    IF(ASSOCIATED(solverDOFToVariableDOFsMap)) &
      & CALL FlagError("The solver DOF to variable DOFs map is already associated.",err,error,*998)

    ALLOCATE(solverDOFToVariableDOFsMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver DOF to variable DOFs map.",err,error,*999)
    solverDOFToVariableDOFsMap%numberOfEquationDOFs=0
    
    EXITS("SolverMappingSDOFToVDOFsMap_Initialise")
    RETURN
999 CALL SolverMappingSDOFToVDOFsMap_Finalise(solverDOFToVariableDOFsMap,dummyErr,dummyError,*998)
998 ERRORS("SolverMappingSDOFToVDOFsMap_Initialise",err,error)
    EXITS("SolverMappingSDOFToVDOFsMap_Initialise")
    RETURN 1
    
  END SUBROUTINE SolverMappingSDOFToVDOFsMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises a solver row to equations map and deallocates all memory.
  SUBROUTINE SolverMappingSRToEQSMap_Finalise(solverRowToEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverRowToEquationsMapType), POINTER :: solverRowToEquationsMap !<A pointer to the solver row to equations map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverMappingSRToEQSMap_Finalise",err,error,*999)

    IF(ASSOCIATED(solverRowToEquationsMap)) THEN
      IF(ALLOCATED(solverRowToEquationsMap%equationsIndex)) DEALLOCATE(solverRowToEquationsMap%equationsIndex)
      IF(ALLOCATED(solverRowToEquationsMap%rowColNumber)) DEALLOCATE(solverRowToEquationsMap%rowColNumber)
      IF(ALLOCATED(solverRowToEquationsMap%couplingCoefficients)) DEALLOCATE(solverRowToEquationsMap%couplingCoefficients)
      DEALLOCATE(solverRowToEquationsMap)
    ENDIF
    
    EXITS("SolverMappingSRToEQSMap_Finalise")
    RETURN
999 ERRORSEXITS("SolverMappingSRToEQSMap_Finalise",err,error)    
    RETURN 1
    
  END SUBROUTINE SolverMappingSRToEQSMap_Finalise

  !
  !================================================================================================================================
  !

  !>Allocates and initialises a solver row to equations map.
  SUBROUTINE SolverMappingSRToEQSMap_Initialise(solverRowToEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverRowToEquationsMapType), POINTER :: solverRowToEquationsMap !<A pointer to the solver row to equations map to initialise. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("SolverMappingSRToEQSMap_Initialise",err,error,*998)

    IF(ASSOCIATED(solverRowToEquationsMap)) CALL FlagError("Solver row to equations map is already associated.",err,error,*998)

    ALLOCATE(solverRowToEquationsMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate the solver row to equations map.",err,error,*999)
    solverRowToEquationsMap%numberOfEquationsSetRows=0
    solverRowToEquationsMap%interfaceConditionIndex=0
        
    EXITS("SolverMappingSRToEQSMap_Initialise")
    RETURN
999 CALL SolverMappingSRToEQSMap_Finalise(solverRowToEquationsMap,dummyErr,dummyError,*998)    
998 ERRORS("SolverMappingSRToEQSMap_Initialise",err,error)    
    EXITS("SolverMappingSRToEQSMap_Initialise")    
    RETURN 1
    
  END SUBROUTINE SolverMappingSRToEQSMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver mapping variable and deallocates all memory.
  SUBROUTINE SolverMappingVariable_Finalise(solverMappingVariable,err,error,*)

    !Argument variables
    TYPE(SolverMappingVariableType), POINTER :: solverMappingVariable !<The solver mapping variable to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverMappingVariable_Finalise",err,error,*999)

    IF(ASSOCIATED(solverMappingVariable)) THEN
      NULLIFY(solverMappingVariable%variable)
      solverMappingVariable%variableType=0
      solverMappingVariable%numberOfEquations=0
      IF(ALLOCATED(solverMappingVariable%equationTypes)) DEALLOCATE(solverMappingVariable%equationTypes)
      IF(ALLOCATED(solverMappingVariable%equationIndices)) DEALLOCATE(solverMappingVariable%equationIndices)
      DEALLOCATE(solverMappingVariable)
    ENDIF
       
    EXITS("SolverMappingVariable_Finalise")
    RETURN
999 ERRORSEXITS("SolverMappingVariable_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMappingVariable_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the solver mapping and deallocates all memory.
  SUBROUTINE SolverMappingVariable_Initialise(solverMappingVariable,err,error,*)

    !Argument variables
    TYPE(SolverMappingVariableType), POINTER :: solverMappingVariable !<The solver mapping variable to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("SolverMappingVariable_Initialise",err,error,*998)

    IF(ASSOCIATED(solverMappingVariable)) CALL FlagError("Solver mapping variable is already associated.",err,error,*998)

    ALLOCATE(solverMappingVariable,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver mapping variable.",err,error,*999)
    NULLIFY(solverMappingVariable%variable)
    solverMappingVariable%variableType=0
    solverMappingVariable%numberOfEquations=0
    
    EXITS("SolverMappingVariable_Initialise")
    RETURN
999 CALL SolverMappingVariable_Finalise(solverMappingVariable,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverMappingVariable_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMappingVariable_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver mapping variables and deallocates all memory.
  SUBROUTINE SolverMappingVariables_Finalise(solverMappingVariables,err,error,*)

    !Argument variables
    TYPE(SolverMappingVariablesType), POINTER :: solverMappingVariables !<The solver mapping variables to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: variableIdx

    ENTERS("SolverMappingVariables_Finalise",err,error,*999)

    IF(ASSOCIATED(solverMappingVariables)) THEN
      solverMappingVariables%numberOfVariables=0
      DO variableIdx=1,SIZE(solverMappingVariables%variables,1)
        CALL SolverMappingVariable_Finalise(solverMappingVariables%variables(variableIdx)%ptr,err,error,*999)
      ENDDO !variableIdx
      DEALLOCATE(solverMappingVariables)
    ENDIF
     
    EXITS("SolverMappingVariables_Finalise")
    RETURN
999 ERRORSEXITS("SolverMappingVariables_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMappingVariables_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the solver mapping variables and deallocates all memory.
  SUBROUTINE SolverMappingVariables_Initialise(solverMappingVariables,err,error,*)

    !Argument variables
    TYPE(SolverMappingVariablesType), POINTER :: solverMappingVariables !<The solver mapping variables to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("SolverMappingVariables_Initialise",err,error,*998)

    IF(ASSOCIATED(solverMappingVariables)) CALL FlagError("Solver mapping variables is already associated.",err,error,*998)

    ALLOCATE(solverMappingVariables,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver mapping variables.",err,error,*999)

    solverMappingVariables%numberOfVariables=0
    
    EXITS("SolverMappingVariables_Initialise")
    RETURN
999 CALL SolverMappingVariables_Finalise(solverMappingVariables,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverMappingVariables_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMappingVariables_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises a variable to solver column map and deallocates all memory.
  SUBROUTINE SolverMappingVDOFToSDOFsMap_Finalise(varDOFToSolverDOFsMap,err,error,*)

    !Argument variables
    TYPE(VariableDOFToSolverDOFsMapType), POINTER :: varDOFToSolverDOFsMap !<The variable to solver column map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverMappingVDOFToSDOFsMap_Finalise",err,error,*999)

    IF(ASSOCIATED(varDOFToSolverDOFsMap)) THEN
      IF(ALLOCATED(varDOFToSolverDOFsMap%dofNumbers)) DEALLOCATE(varDOFToSolverDOFsMap%dofNumbers)
      IF(ALLOCATED(varDOFToSolverDOFsMap%couplingCoefficients)) DEALLOCATE(varDOFToSolverDOFsMap%couplingCoefficients)
      IF(ALLOCATED(varDOFToSolverDOFsMap%additiveConstants)) DEALLOCATE(varDOFToSolverDOFsMap%additiveConstants)
      DEALLOCATE(varDOFToSolverDOFsMap)
    ENDIF
      
    EXITS("SolverMappingVDOFToSDOFsMap_Finalise")
    RETURN
999 ERRORS("SolverMappingVDOFToSDOFsMap_Finalise",err,error)    
    EXITS("SolverMappingVDOFToSDOFsMap_Finalise")    
    RETURN 1
   
  END SUBROUTINE SolverMappingVDOFToSDOFsMap_Finalise
  
  !
  !================================================================================================================================
  !

  !>Initialises a variable to solver column map.
  SUBROUTINE SolverMappingVDOFToSDOFsMap_Initialise(varDOFToSolverDOFsMap,err,error,*)

    !Argument variables
    TYPE(VariableDOFToSolverDOFsMapType), POINTER :: varDOFToSolverDOFsMap !<The variable to solver column map to initalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("SolverMappingVDOFToSDOFsMap_Initialise",err,error,*998)

    IF(ASSOCIATED(varDOFToSolverDOFsMap)) CALL FlagError("Variable DOF to solver DOFs map is already associated.",err,error,*998)

    ALLOCATE(varDOFToSolverDOFsMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate variable DOF to solver DOFs map.",err,error,*999)
    !Nothing to do here, all members are allocatable

    EXITS("SolverMappingVDOFToSDOFsMap_Initialise")
    RETURN
999 CALL SolverMappingVDOFToSDOFsMap_Finalise(varDOFToSolverDOFsMap,dummyErr,dummyError,*998)
998 ERRORS("SolverMappingVDOFToSDOFsMap_Initialise",err,error)    
    EXITS("SolverMappingVDOFToSDOFsMap_Initialise")    
    RETURN 1
   
  END SUBROUTINE SolverMappingVDOFToSDOFsMap_Initialise

  !
  !================================================================================================================================
  !

  !>Add a new DOF coupling to the list of global solver couplings
  SUBROUTINE SolverDofCouplings_AddCoupling(solverDofCouplings,coupling,couplingIndex,err,error,*)

    !Argument variables
    TYPE(SolverMappingDofCouplingsType), INTENT(INOUT) :: solverDofCouplings !<The solver row or column couplings
    TYPE(BoundaryConditionsCoupledDofsType), POINTER, INTENT(IN) :: coupling !<The new DOF boundary condition coupling to add to the list
    INTEGER(INTG), INTENT(OUT) :: couplingIndex !<On return, the index of the coupling in the solver dof couplings
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: newSize
    TYPE(BoundaryConditionsCoupledDofsPtrType), ALLOCATABLE :: newDofCouplings(:)
    !Local Variables

    ENTERS("SolverDofCouplings_AddCoupling",err,error,*998)

    IF(solverDofCouplings%numberOfCouplings+1>solverDofCouplings%capacity) THEN
      !Resize or perform initial allocation if necessary
      IF(solverDofCouplings%capacity==0) THEN
        newSize=32
      ELSE
        newSize=solverDofCouplings%capacity*2
      END IF
      ALLOCATE(newDofCouplings(newSize),stat=err)
      IF(err/=0) CALL FlagError("Could not allocate new DOF couplings array.",err,error,*999)
      IF(solverDofCouplings%capacity>0) THEN
        newDofCouplings(1:solverDofCouplings%numberOfCouplings)= &
          & solverDofCouplings%dofCouplings(1:solverDofCouplings%numberOfCouplings)
      END IF
      solverDofCouplings%capacity=newSize
      CALL MOVE_ALLOC(newDofCouplings,solverDofCouplings%dofCouplings)
    END IF

    solverDofCouplings%dofCouplings(solverDofCouplings%numberOfCouplings+1)%ptr=>coupling
    solverDofCouplings%numberOfCouplings=solverDofCouplings%numberOfCouplings+1
    couplingIndex=solverDofCouplings%numberOfCouplings

    EXITS("SolverDofCouplings_AddCoupling")
    RETURN
999 IF(ALLOCATED(newDofCouplings)) DEALLOCATE(newDofCouplings)
998 ERRORSEXITS("SolverDofCouplings_AddCoupling",err,error)
    RETURN 1

  END SUBROUTINE SolverDofCouplings_AddCoupling

  !
  !================================================================================================================================
  !

  !>Initialise the list of solver row or column couplings
  SUBROUTINE SolverMappingDOFCouplings_Finalise(solverMappingDOFCouplings,err,error,*)

    !Argument variables
    TYPE(SolverMappingDofCouplingsType), INTENT(INOUT) :: solverMappingDOFCouplings !<The solver row or column couplings to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverMappingDOFCouplings_Finalise",err,error,*999)

    IF(ALLOCATED(solverMappingDOFCouplings%dofCouplings)) THEN
      !Don't finalise individual DOF couplings as these are owned
      !by the boundary conditions.
      DEALLOCATE(solverMappingDOFCouplings%dofCouplings)
    END IF
    solverMappingDOFCouplings%numberOfCouplings=0
    solverMappingDOFCouplings%capacity=0

    EXITS("SolverMappingDOFCouplings_Finalise")
    RETURN
999 ERRORSEXITS("SolverMappingDOFCouplings_Finalise",err,error)
    RETURN 1

  END SUBROUTINE SolverMappingDOFCouplings_Finalise

  !
  !================================================================================================================================
  !

  !>Initialise the list of solver row or column couplings
  SUBROUTINE SolverMappingDOFCouplings_Initialise(solverMappingDOFCouplings,err,error,*)

    !Argument variables
    TYPE(SolverMappingDofCouplingsType), INTENT(INOUT) :: solverMappingDOFCouplings !<The solver row or column couplings to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverMappingDOFCouplings_Initialise",err,error,*999)

    solverMappingDOFCouplings%numberOfCouplings=0
    solverMappingDOFCouplings%capacity=0

    EXITS("SolverMappingDOFCouplings_Initialise")
    RETURN
999 ERRORSEXITS("SolverMappingDOFCouplings_Initialise",err,error)
    RETURN 1

  END SUBROUTINE SolverMappingDOFCouplings_Initialise
  
  !
  !================================================================================================================================
  !

END MODULE SolverMappingRoutines
