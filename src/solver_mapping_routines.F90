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
  USE ComputationRoutines
  USE ComputationAccessRoutines
  USE DistributedMatrixVector
  USE EquationsAccessRoutines
  USE DomainMappings
  USE EquationsMappingAccessRoutines
  USE EquationsSetAccessRoutines
  USE FieldAccessRoutines
  USE InputOutput
  USE InterfaceConditionAccessRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE Lists
  USE MatrixVector
  USE ProblemAccessRoutines
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

  PUBLIC SOLVER_MAPPING_INTERFACE_CONDITION_ADD

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
      & equationsColumn,equationsIdx,equationsIdx2,equationsMatrixNumber,equationMatrixIdx,equationsRowNumber, &
      & equationsSetIdx,equationsVariableListItem(3),globalColumn,globalDOF,globalDOFIdx,globalDOFsOffset, &
      & globalRow,globalRowIdx,interfaceColumn,interfaceColumnNumber,interfaceConditionIdx,interfaceConditionIdx2, &
      & interfaceEquationsListItem(2),interfaceMatrixIdx,interfaceRow,interfaceRowNumber,jacobianColumn, &
      & localColumn,localDOF,localDOFsOffset,localRow,matricesType,matrixNumber,matrixType,matrixTypeIdx, &
      & matrixVariableIdx,myrank,numberOfColumns,numberOfDynamicEquationsMatrices,numberOfEquationsColumns, &
      & numberOfEquationsSets,numberOfEquationsVariables,numberOfInterfaceConditions,numberOfInterfaceColumns, &
      & numberOfInterfaceRows,numberOfInterfaceVariables,numberOfGlobalSolverDOFs,numberOfGlobalSolverRows, &
      & numberOfLinearMatrices,numberOfLocalSolverDOFs,numberOfLocalSolverRows,numberOfRankColumns, &
      & numberOfRankRows,numberOfVariables,rank,rankIdx,rowIdx,rowListItem(4),rowRank,solverGlobalDOF, &
      & solverMatrixIdx,solverVariableIdx,totalNumberOfLocalSolverDOFs,variableIdx, &
      & variableListItem(3),variablePositionIdx,variableType,numberOfGroupComputationNodes, &
      & numberRowEquationsRows,numberColEquationsCols,rowEquationsRowIdx,colEquationsColIdx, &
      & globalDofCouplingNumber,equationsRow,eqnLocalDof,numberOfEquationsRHSVariables,rhsVariableType,equationsSetIdx
    INTEGER(INTG) :: tempOffset, tempSolverVariableIdx
    INTEGER(INTG), ALLOCATABLE :: equationsSetVariables(:,:),equationsVariables(:,:),interfaceEquationsList(:,:), &
      & interfaceVariables(:,:),rankGlobalRowsList(:,:),rankGlobalColumnsList(:,:),solverLocalDOF(:), &
      & equationsRHSVariables(:,:)
    INTEGER(INTG), ALLOCATABLE :: numberOfVariableGlobalSolverDOFs(:),numberOfVariableLocalSolverDOFs(:), &
      & totalNumberOfVariableLocalSolverDOFs(:),subMatrixInformation(:,:,:),subMatrixList(:,:,:),variableTypes(:)
    REAL(DP) :: couplingCoefficient
    LOGICAL :: found,includeColumn,includeRow,constrainedDOF
    LOGICAL, ALLOCATABLE :: variableProcessed(:),variableRankProcessed(:,:)
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable
    TYPE(DomainMappingType), POINTER :: columnDomainMapping,columnDOFsMapping,rowDomainMapping,rowDOFsMapping
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLinearType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residulMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: equationsMatrixToSolverMatrixMap
    TYPE(FieldType), POINTER :: dependentField,lagrangeField
    TYPE(FieldVariableType), POINTER :: dependentVariable,lagrangeVariable,VARIABLE,rhsVariable
    TYPE(IntegerIntgPtrType), POINTER :: dofMap(:)
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(InterfaceMatrixToSolverMatrixMapType), POINTER :: interfaceMatrixToSolverMatrixMap
    TYPE(JacobianMatrixToSolverMatrixMapType), POINTER :: jacobianMatrixToSolverMatrixMap
    TYPE(ListType), POINTER :: equationsSetVariableList
    TYPE(ListPtrType), ALLOCATABLE :: interfaceEquationsLists(:),rankGlobalRowsLists(:,:), &
      & rankGlobalColumnLists(:,:,:,:),variablesList(:)
    TYPE(SOLVER_TYPE), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints
    TYPE(BoundaryConditionsCoupledDofsType), POINTER :: rowEquationRows,colEquationCols
    TYPE(BoundaryConditionsCoupledDofsType), TARGET :: dummyDofCoupling
    TYPE(SolverMappingDofCouplingsType) :: rowCouplings
    TYPE(SolverMappingDofCouplingsType) :: columnCouplings
    TYPE(VARYING_STRING) :: localError
    TYPE(WorkGroupType), POINTER :: workGroup

    ENTERS("SolverMapping_Calculate",err,error,*999)

    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)

    NULLIFY(createValuesCache)
    CALL SolverMapping_CreateValuesCacheGet(solverMapping,createValuesCache,err,error,*999)
    NULLIFY(solverEquations)
    CALL SolverMapping_SolverEquationsGet(solverMapping,solveEquations,err,error,*999)
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
      CALL SolverMapping_InterfaceConditionToSolverMatricesMapInitialise(solverMapping% &
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
      CALL EqutionsSet_EquationsGet(equationsSet,equations,err,error,*999)
      NULLIFY(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr)
      CALL SolverMapping_EquationsSetToSolverMapInitialise(solverMapping% &
        & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr,err,error,*999)
      solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%equationsSetIndex=equationsSetIdx
      solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%solverMapping=>solverMapping
      solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%equations=>equations
      !Set up list of interface conditions affecting this equations set
      CALL List_DetachAndDestroy(createValuesCache%interfaceIndices(equationsSetIdx)%ptr,numberOfInterfaceConditions, &
        & interfaceEquationsList,err,error,*999)
      ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr &
        & interfaceConditionToEquationsSetMaps(numberOfInterfaceConditions),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations to solver maps interface.",err,error,*999)
      solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%numberOfInterfaceConditions=numberOfInterfaceConditions
      DO interfaceConditionIdx=1,numberOfInterfaceConditions
        CALL SolverMapping_EquationsToSolverInterfaceInitialise(solverMapping%equationsSetToSolverMatricesMaps( &
          & equationsSetIdx)%ptr%interfaceConditionToEquationsSetMaps(interfaceConditionIdx),err,error,*999)
        interfaceConditionIdx=interfaceEquationsList(1,interfaceConditionIdx)
        interfaceMatrixIdx=interfaceEquationsList(2,interfaceConditionIdx)
        NULLIFY(interfaceCondition)
        CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondtion,err,error,*999)
        solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%interfaceConditionToEquationsSetMaps( &
          & interfaceConditionIdx)%interfaceConditionIndex=interfaceConditionIdx
        solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%interfaceConditionToEquationsSetMaps( &
          & interfaceConditionIdx)%interfaceCondition=>interfaceCondition
        solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%interfaceConditionToEquationsSetMaps( &
          & interfaceConditionIdx)%interfaceMatrixNumber=interfaceMatrixIdx
        interfaceEquationsListItem(1)=equationsSetIdx
        interfaceEquationsListItem(2)=interfaceMatrixIdx
        CALL List_ItemAdd(interfaceEquationsLists(interfaceConditionIdx)%ptr,interfaceEquationsListItem,err,error,*999)
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
        CALL SolverMapping_InterfaceToSolverEquationsInitialise(solverMapping%interfaceConditionToSolverMatricesMaps( &
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
    NULLIFY(rowCouplings)
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
      CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcseMapping,err,error,*999)
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
            NULLIFY(dofCoupling)
            CALL BoundaryConditionsVariableDOFConstraints_DOFCouplingExists(dofConstraints,globalDOF,dofCoupling,err,error,*999)
            IF(ASSOCIATED(dofConstraint)) THEN
              !This equations row is the owner of a solver row that is mapped to
              !multiple other equations rows, add it to the list of global row
              !couplings and remember the index into the global list for this solver row
              CALL SolverDofCouplings_AddCoupling(rowCouplings,dofConstraint,globalDofCouplingNumber,err,error,*999)
            ENDIF
          ENDIF
        ENDIF
        IF(ASSOCIATED(nonlinearMapping)) THEN
          !Look at the boundary conditions for nonlinear variables for this row
          DO residualIdx=1,nonlinearMapping%numberOfResiduals
            NULLIFY(residualMapping)
            CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualMapping,err,error,*999)
            DO variableIdx=1,residualMapping%numberOfVariables
              NULLIFY(dependentVariable)
              CALL EquationsMappingNonlinear_ResidualVariableGet(nonlinearMapping,variableIdx,residualIdx,dependentVariable, &
                & err,error,*999)
              NULLIFY(boundaryConditionsVariable)
              CALL BoundaryConditions_VariableGet(boundaryConditions,dependentVariable,boundaryConditionsVariable,err,error,*999)
              globalDOF=globalRow
              includeRow=includeRow.AND.(boundaryConditionsVariable%DOFTypes(globalDOF)==BOUNDARY_CONDITION_DOF_FREE)
              constrainedDOF=constrainedDOF.OR.(boundaryConditionsVariable%DOFTypes(globalDOF)==BOUNDARY_CONDITION_DOF_CONSTRAINED)
              NULLIFY(dofConstraints)
              CALL BoundaryConditionsVariable_DOFConstraintsExists(boundaryConditionsVariable,dofConstraints,err,error,*999)
              IF(ASSOCIATED(dofConstraints)) THEN
                NULLIFY(dofCoupling)
                CALL BoundaryConditionsVariableDOFConstraints_DOFCouplingExists(dofConstraints,globalDOF,dofCoupling, &
                  & err,error,*999)
                IF(ASSOCIATED(dofConstraint)) THEN
                  CALL SolverDofCouplings_AddCoupling(rowCouplings,dofConstraint,globalDofCouplingNumber,err,error,*999)
                ENDIF
              ENDIF
            ENDDO !variableIdx
          ENDDO !residualIdx
        ENDIF
        IF(ASSOCIATED(linearMapping)) THEN
          !Loop over the variables in the equations set. Don't include the row in the solver matrices if
          !all the variable dofs associated with this equations row are fixed.
          DO variableIdx=1,linearMapping%numberOfVariables
            NULLIFY(dependentVariable)
            CALL EquationsMappingLinear_VariableGet(linearMapping,variableIdx,dependentVariable,err,error,*999)
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
              NULLIFY(dofCoupling)
              CALL BoundaryConditionsVariableDOFConstraints_DOFCouplingExists(dofConstraints,globalDOF,dofCoupling,err,error,*999)
              IF(ASSOCIATED(dofConstraint)) THEN
                CALL SolverDofCouplings_AddCoupling(rowCouplings,dofConstraint,globalDofCouplingNumber,err,error,*999)
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
        CALL Field_VariableGet(langrangeField,variableType,lagrangeVariable,err,error,*999)
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
    ALLOCATE(solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver mapping solver row to equation rows map.",err,error,*999)
    !Set the number of rows
    solverMapping%numberOfRows=numberOfLocalSolverRows
    solverMapping%numberOfGlobalRows=numberOfGlobalSolverRows
    !Allocate the solver rows domain mapping
    CALL DomainMapping_Initialise(solverMapping%rowDOFsMapping,err,error,*999)
    CALL DomainMapping_WorkGroupSet(solverMapping%rowDOFsMapping,workGroup,err,error,*999)
    NULLIFY(rowDomainMapping)
    CALL SolverMapping_RowDomainMappingGet(solverMapping,rowDomainMapping,err,error,*999)
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
        CALL SolverMapping_EquatsToSolMatMapsSMInitialise(solverMapping% &
          & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
          & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr,err,error,*999)
      ENDDO !solverMatrixIdx
      
      !Allocate the equations row to solver rows maps
      ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
        & equationsRowToSolverRowsMap(vectorMapping%totalNumberOfRows),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations set set to solver matrices map equations row to solver rows maps.", &
        & err,error,*999)
      DO equationsRowNumber=1,vectorMapping%totalNumberOfRows
        !Initialise
        CALL MatrixRowColCoupling_Initialise(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
          & equationsRowToSolverRowsMap(equationsRowNumber),err,error,*999)
      ENDDO !equationsRowIdx
             
    ENDDO !equationsSetIdx
          
    !Initialise the interface condition to solver matrices maps
    DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
      NULLIFY(interfaceCondition)
      CALL SolverMatrices_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
      NULLIFY(interfaceEquations)
      CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
      NULLIFY(interfaceMapping)
      CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,err,error,*999)
      
      !Allocate the interface matrices to solver matrix maps
      ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
        & interfaceMatricesToSolverMatrixMaps(solverMapping%numberOfSolverMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate interface condition to solver matrices map interface matrices"// &
        & " to solver matrix maps.",err,error,*999)
      DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
        NULLIFY(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
          & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr)
        CALL SolverMapping_InterfToSolMatMapsSMInitialise(solverMapping% &
          & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
          & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr,err,error,*999)
      ENDDO !solverMatrixIdx
            
      !Allocate the interface matrix to solver matrices maps
      ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
        & interfaceMatrixToSolverMatricesMaps(interfaceMapping%numberOfInterfaceMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate interface condition to solver matrices map interface matrix"// &
        & " to solver matrices maps.",err,error,*999)
      DO interfaceMatrixIdx=1,interfaceMapping%numberOfInterfaceMatrices
        NULLIFY(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
          & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr
        CALL SolverMapping_InterfToSolMatMapsIMInitialise(solverMapping% &
          & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
          & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr,err,error,*999)
                      
        !Allocate the interface row to solver row maps
        ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
          & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr%&
          & interfaceRowToSolverRowsMap(interfaceMapping%interfaceMatrixRowsToVarMaps(interfaceMatrixIdx)%totalNumberOfRows), &
          & STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate interface condition to solver matrices map"// &
          & " interface matrix to solver matrices maps interface row to solver rows map.",err,error,*999)
        DO interfaceRowNumber=1,interfaceMapping%interfaceMatrixRowsToVarMaps(interfaceMatrixIdx)%totalNumberOfRows
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
            NULLIFY(rowEquationsRows)
            CALL BoundaryConditionsVariablesDOFCouplings_DOFCouplingGet(rowCoupling,globalDOFCouplingNumber,rowEquationRows, &
              & err,error,*999)
            numberRowEquationsRows=rowEquationRows%numberOfDofs
          ELSE
            numberRowEquationsRows=1
            dummyDofCoupling%globalDofs(1)=globalRow
            dummyDofCoupling%localDofs(1)=localRow
            dummyDofCoupling%coefficients(1)=1.0_DP
            rowEquationRows=>dummyDofCoupling
          END IF

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
              CALL SolverMappingSRToERSMap_Initialise(solverMapping% &
                & solverRowToEquationsRowsMap(numberOfLocalSolverRows),err,error,*999)
              
              ALLOCATE(solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)% &
                & equationsIndex(numberRowEquationsRows),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate solver row to equations rows equations index.",err,error,*999)
              ALLOCATE(solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)% &
                & rowColNumber(numberRowEquationsRows),STAT=err)
              IF(err/=0) &
                & CALL FlagError("Could not allocate solver row to equations rows row/col number.",err,error,*999)
              ALLOCATE(solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)% &
                & couplingCoefficients(numberRowEquationsRows),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate solver row to equations rows coupling coefficients.", &
                & err,error,*999)
              !Set the mappings for the first equations DOF, the rest will be set up later using the DOF constraints
              solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)%numberOfEquationsSetRows=numberRowEquationsRows
              DO rowEquationsRowIdx=1,numberRowEquationsRows
                solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)% &
                  & equationsIndex(rowEquationsRowIdx)=equationsSetIdx
                solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)% &
                  & rowColNumber(rowEquationsRowIdx)=rowEquationRows%localDofs(rowEquationsRowIdx)
                solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)% &
                  & couplingCoefficients(rowEquationsRowIdx)=rowEquationRows%coefficients(rowEquationsRowIdx)
              ENDDO !rowEquationsRowIdx
              !Set up the equations row -> solver row mappings
              DO rowEquationsRowIdx=1,numberRowEquationsRows
                equationsRow=rowEquationRows%localDofs(rowEquationsRowIdx)
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
                  & rowEquationRows%coefficients(rowEquationsRowIdx)
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
                    & rowEquationRows%coefficients(rowEquationsRowIdx)
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
        CALL InterfaceConditio_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
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
              !Note that for populating solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)%rowColNumber(i)
              !If the row are equations set rows this is the i'th row number that the solver row is mapped to.
              
              !If the rows are interface rows (which is the case here) then this is the i'th column number that the solver
              !row is mapped to.
              
              !Initialise
              CALL SolverMappingSRToERSMap_Initialise(solverMapping% &
                & solverRowToEquationsRowsMap(numberOfLocalSolverRows),err,error,*999)
              ALLOCATE(solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)%rowColNumber(1),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate solver row to equations rows row/col number.",err,error,*999)
              ALLOCATE(solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)%couplingCoefficients(1),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate solver row to equations rows coupling coefficients.",err,error,*999)
              !Set up the interface column -> solver row mappings
              !/todo the solverRowToEquationsRowsMap may need to be renamed for clarity
              solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)%interfaceConditionIndex=interfaceConditionIdx
              solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)%rowColNumber(1)=localColumn
              solverMapping%solverRowToEquationsRowsMap(numberOfLocalSolverRows)%couplingCoefficients(1)=1.0_DP
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
    !Allocate solver column to equations sets mapping array
    ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMapping%numberOfSolverMatrices),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver mapping solver column to equations column maps.",err,error,*999)
    ALLOCATE(solverMapping%variablesList(solverMapping%numberOfSolverMatrices),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver mapping variables list.",err,error,*999)
    CALL SolverMappingDOFCouplings_Initialise(columnCouplings,err,error,*999)

    !Calculate the column mappings for each solver matrix
    DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices

      !Initialise the variables list
      CALL SolverMappingVariables_Initialise(solverMapping%variablesList(solverMatrixIdx),err,error,*999)
      !
      ! 4a Calculate the list of field variables involved in the columns of the solver matrix
      !
      !Compute the order of variables for the solver matrices
      CALL List_DetachAndDestroy(createValuesCache%equationsVariableList(solverMatrixIdx)%ptr,numberOfEquationsVariables, &
        & equationsVariables,err,error,*999)
      CALL List_DetachAndDestroy(createValuesCache%interfaceVariableList(solverMatrixIdx)%ptr,numberOfInterfaceVariables, &
        & interfaceVariables,err,error,*999)
      ALLOCATE(solverMapping%variablesList(solverMatrixIdx)%variables(numberOfEquationsVariables+numberOfInterfaceVariables), &
        & STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate variables list variables.",err,error,*999)
      solverMapping%variablesList(solverMatrixIdx)%numberOfVariables=numberOfEquationsVariables+numberOfInterfaceVariables
      ALLOCATE(variablesList(numberOfEquationsVariables+numberOfInterfaceVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate variables list.",err,error,*999)
      ALLOCATE(variableProcessed(numberOfEquationsVariables+numberOfInterfaceVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate variable processed.",err,error,*999)
      variableProcessed=.FALSE.
      solverVariableIdx=0
      DO variableIdx=1,numberOfEquationsVariables
        solverVariableIdx=solverVariableIdx+1
        CALL SolverMappingVariable_Initialise(solverMapping%variablesList(solverMatrixIdx)%variables(solverVariableIdx), &
          & err,error,*999)
        equationsSetIdx=equationsVariables(1,variableIdx)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        variableType=equationsVariables(2,variableIdx)
        NULLIFY(variable)
        CALL Field_VariableGet(dependentField,variableType,variable,err,error,*999)
        solverMapping%variablesList(solverMatrixIdx)%variables(solverVariableIdx)%variable=>variable
        solverMapping%variablesList(solverMatrixIdx)%variables(solverVariableIdx)%variableType=variableType
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
        CALL SolverMappingVariable_Initialise(solverMapping%variablesList(solverMatrixIdx)%variable(solverVariableIdx), &
          & err,error,*999)
        interfaceConditionIdx=interfaceVariables(1,variableIdx)
        NULLIFY(interfaceCondition)
        CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
        NULLIFY(lagrangeField)
        CALL InterfaceCondition_LagrangeFieldGet(interfaceCondition,lagrangeField,err,error,*999)
        variableType=interfaceVariables(2,variableIdx)
        NULLIFY(variable)
        CALL Field_VariableGet(lagrangeField,variableType,variable,err,error,*999)
        solverMapping%variablesList(solverMatrixIdx)%variables(solverVariableIdx)%variable=>variable
        solverMapping%variablesList(solverMatrixIdx)%variables(solverVariableIdx)%variableType=variableType
        NULLIFY(variablesList(solverVariableIdx)%ptr)
        CALL List_CreateStart(variablesList(solverVariableIdx)%ptr,err,error,*999)
        CALL List_DataTypeSet(variablesList(solverVariableIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
        CALL List_DataDimensionSet(variablesList(solverVariableIdx)%ptr,3,err,error,*999)
        CALL List_KeyDimensionSet(variablesList(solverVariableIdx)%ptr,1,err,error,*999)
        CALL List_CreateFinish(variablesList(solverVariableIdx)%ptr,err,error,*999)
      ENDDO !variableIdx
      IF(ALLOCATED(interfaceVariables)) DEALLOCATE(interfaceVariables)
      !Handle the RHS variables
      CALL List_DetachAndDestroy(createValuesCache%equationsRHSVariablesList,numberOfEquationsRHSVariables, &
        & equationsRHSVariables,err,error,*999)
      solverMapping%rhsVariablesList%numberOfVariables=numberOfEquationsRHSVariables
      ALLOCATE(solverMapping%rhsVariablesList%variables(numberOfEquationsRHSVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate RHS variables list.",err,error,*999)
      DO variableIdx=1,numberOfEquationsRHSVariables
        CALL SolverMappingVariable_Initialise(solverMapping%rhsVariablesList%variables(variableIdx),err,error,*999)
        equationsSetIdx=equationsRHSVariables(1,variableIdx)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        rhsVariableType=equationsRHSVariables(2,variableIdx)
        NULLIFY(rhsVariable)
        CALL Field_VariableGet(dependentField,rhsVariableType,rhsVariable,err,error,*999)
        solverMapping%rhsVariablesList%variables(variableIdx)%variable=>rhsVariable
        solverMapping%rhsVariablesList%variables(variableIdx)%variableType=rhsVariableType
!!TODO: Redo this properly with left hand side mapping. For now just add this equations set.
        solverMapping%rhsVariablesList%variables(variableIdx)%numberOfEquations=0
      ENDDO !variableIdx
      IF(ALLOCATED(equationsRHSVariables)) DEALLOCATE(equationsRHSVariables)
      !Initialise solver column to equations sets mapping array
      CALL SolverMappingSMToEquationsMap_Initialise(solverMapping%solverColToEquationsColsMap(solverMatrixIdx),err,error,*999)
      solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverMatrixNumber=solverMatrixIdx
      solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverMapping=>solverMapping
      !Allocate the solver col to equations set maps array
      ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
        & solverMatrixToEquationsSetMaps(solverMapping%numberOfEquationsSets),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate solver col to equations map solver col to equation set maps.", &
        & err,error,*999)
      !Allocate the solver col to interface maps array
      ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
        & solverColToInterfaceConditionMaps(solverMapping%numberOfInterfaceConditions),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate solver col to equations map solver col to interface maps.", &
        & err,error,*999)
      !Presort the column numbers by rank.
      !rankGlobalColumnLists(dofType, equationsIdx, variableIdx, rankIdx)
      !dofType is 1 for domain local DOFs and 2 for ghost DOFs
      ALLOCATE(rankGlobalColumnLists(2,solverMapping%numberOfEquationsSets+solverMapping%numberOfInterfaceConditions, &
        & solverMapping%variablesList(solverMatrixIdx)%numberOfVariables,0:numberOfGroupComputationNodes-1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate rank global columns lists.",err,error,*999)
      DO rank=0,numberOfGroupComputationNodes-1
        DO solverVariableIdx=1,solverMapping%variablesList(solverMatrixIdx)%numberOfVariables
          DO equationsIdx=1,solverMapping%numberOfEquationsSets+solverMapping%numberOfInterfaceConditions
            DO dofType=1,2
              NULLIFY(rankGlobalColumnLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr)
              CALL List_CreateStart(rankGlobalColumnLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr,err,error,*999)
              CALL List_DataTypeSet(rankGlobalColumnLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr,LIST_INTG_TYPE, &
                & err,error,*999)
              !Set approximate size for the number of columns per variable.
              CALL List_InitialSizeSet(rankGlobalColumnLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr, &
                & solverMapping%variablesList(solverMatrixIdx)%variables(solverVariableIdx)%variable%totalNumberOfDofs, &
                & err,error,*999)
              CALL List_DataDimensionSet(rankGlobalColumnLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr,5,err,error,*999)
              CALL List_KeyDimensionSet(rankGlobalColumnLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr,1,err,error,*999)
              CALL List_CreateFinish(rankGlobalColumnLists(dofType,equationsIdx,solverVariableIdx,rank)%ptr,err,error,*999)
            ENDDO !dofType
          ENDDO !equationsIdx
        ENDDO !solverVariableIdx
      ENDDO !rank

      !Allocate sub-matrix information
      !subMatrixInformation(1,equationsIdx,variableIdx) = The equations type, see SOLVER_MAPPING_EquationsTypes
      !subMatrixInformation(2,equationsIdx,variableIdx) = The equations set or interface condition index
      !subMatrixInformation(3,equationsIdx,variableIdx) = The interface matrix index, or 0 for an equations set matrix
      !equationsIdx goes from 1 to the number of equations sets + interface conditions
      !variableIdx goes from 1 to the number of variables mapped to this solver matrix
      ALLOCATE(subMatrixInformation(3,solverMapping%numberOfEquationsSets+solverMapping%numberOfInterfaceConditions, &
        & solverMapping%variablesList(solverMatrixIdx)%numberOfVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate sub matrix information.",err,error,*999)
      subMatrixInformation=0
      !Allocate sub-matrix list information
      ALLOCATE(subMatrixList(0:3,solverMapping%numberOfEquationsSets+solverMapping%numberOfInterfaceConditions, &
        & solverMapping%variablesList(solverMatrixIdx)%numberOfVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate sub matrix list.",err,error,*999)
      subMatrixList=0

      !
      ! 4b Calculate the number of columns
      !
      
      !Calculate the number of solver dofs
      ALLOCATE(numberOfVariableGlobalSolverDOFs(solverMapping%variablesList(solverMatrixIdx)%numberOfVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate number of global solver dofs.",err,error,*999)
      ALLOCATE(numberOfVariableLocalSolverDOFs(solverMapping%variablesList(solverMatrixIdx)%numberOfVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate number of local solver dofs.",err,error,*999)
      ALLOCATE(totalNumberOfVariableLocalSolverDOFs(solverMapping%variablesList(solverMatrixIdx)%numberOfVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate total number of local solver dofs.",err,error,*999)

      numberOfVariableGlobalSolverDOFs=0
      numberOfVariableLocalSolverDOFs=0
      totalNumberOfVariableLocalSolverDOFs=0
      
      equationsIdx=0
      !Loop over the equations sets
      DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
        equationsIdx=equationsIdx+1
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equtionsSet,err,error,*999)
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
          DO variableIdx=1,createValuesCache%residualVariableTypes(0,equationsSetIdx)
            equationsVariableListItem(1)=createValuesCache%residualVariableTypes(variableIdx,equationsSetIdx)
            equationsVariableListItem(2)=SOLVER_MAPPING_EQUATIONS_NONLINEAR_MATRIX
            equationsVariableListItem(3)=variableIdx
            CALL List_ItemAdd(equationsSetVariableList,equationsVariableListItem,err,error,*999)
          ENDDO !variableIdx
        ENDIF
        CALL List_RemoveDuplicates(equationsSetVariableList,err,error,*999)
        CALL List_DetachAndDestroy(equationsSetVariableList,numberOfVariables,equationsSetVariables,err,error,*999)
        !Initialise equations set to solver matrices equations matrices to solver matrix map
        NULLIFY(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%, &
          & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr)
        CALL SolverMapping_EquatsToSolMatMapsSMInitialise(solverMapping% &
          & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
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
          & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%variableToSolverColMaps(numberOfVariables),STAT=err)
        IF(err/=0)  &
          & CALL FlagError("Could not allocate equations matrices to solver matrix maps variable to solver col maps.", &
          & err,error,*999)
        !Setup
        solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
          & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfVariables=numberOfVariables
        numberOfDynamicEquationsMatrices=0
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
          DO variablePositionIdx=1,solverMapping%variablesList(solverMatrixIdx)%numberOfVariables
            IF(ASSOCIATED(dependentVariable,solverMapping%variablesList(solverMatrixIdx)% &
              & variables(variablePositionIdx)%variable)) THEN
              found=.TRUE.
              EXIT
            ENDIF
          ENDDO !variablePositionIdx
          IF(.NOT.found) THEN
            localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
              & " in equations set index "//TRIM(NumberToVString(equationsSetIdx,"*",err.error))// &
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
          !Allocate the variable to solver col maps arrays
          ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%variableToSolverColMaps(variableIdx)% &
            & columnNumbers(dependentVariable%totalNumberOfDofs),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate variables to solver column maps column numbers.",err,error,*999)
          ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%variableToSolverColMaps(variableIdx)% &
            & couplingCoefficients(dependentVariable%totalNumberOfDofs),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate variables to solver column maps coupling coefficients.",err,error,*999)
          ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%variableToSolverColMaps(variableIdx)% &
            & additiveConstants(dependentVariable%totalNumberOfDofs),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate variables to solver column maps additive constants.",err,error,*999)
          !Setup
          !Set the sub-matrix information
          subMatrixInformation(1,equationsIdx,variablePositionIdx)=SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET
          subMatrixInformation(2,equationsIdx,variablePositionIdx)=equationsSetIdx
          !Set the sub-matrix lists
          IF(ASSOCIATED(dynamicMapping)) THEN
            numberOfDynamicEquationsMatrices=numberOfDynamicEquationsMatrices+ &
              & dynamicMapping%varToEquationsMatricesMap%numberOfEquationsMatrices
            IF(dynamicMapping%varToEquationsMatricesMap%numberOfEquationsMatrices>0) THEN
              subMatrixList(0,equationsIdx,variablePositionIdx)=subMatrixList(0,equationsIdx,variablePositionIdx)+1
              subMatrixList(subMatrixList(0,equationsIdx,variablePositionIdx),equationsIdx,variablePositionIdx)= &
                & SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX
            ENDIF
          ENDIF
          IF(ASSOCIATED(linearMapping)) THEN
            variableIdx=linearMapping%linearVariableTypesMap(variableType)
            numberOfLinearMatrices=numberOfLinearMatrices+ &
              & linearMapping%varToEquationsMatricesMaps(variableIdx)%numberOfEquationsMatrices
            IF(linearMapping%varToEquationsMatricesMaps(variableIdx)%numberOfEquationsMatrices>0) THEN
              subMatrixList(0,equationsIdx,variablePositionIdx)=subMatrixList(0,equationsIdx,variablePositionIdx)+1
              subMatrixList(subMatrixList(0,equationsIdx,variablePositionIdx),equationsIdx,variablePositionIdx)= &
                & SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX
            ENDIF
          ENDIF
          IF(ASSOCIATED(nonlinearMapping)) THEN
            DO equationMatrixIdx=1,nonlinearMapping%numberOfResidualVariables
              IF(nonlinearMapping%varToJacobianMap(equationMatrixIdx)%variableType==variableType) THEN
                subMatrixList(0,equationsIdx,variablePositionIdx)=subMatrixList(0,equationsIdx,variablePositionIdx)+1
                subMatrixList(subMatrixList(0,equationsIdx,variablePositionIdx),equationsIdx,variablePositionIdx)= &
                  & SOLVER_MAPPING_EQUATIONS_NONLINEAR_MATRIX
              ENDIF
            ENDDO !equationsMatrixIdx
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
              CALL BoundaryConditionsVariable_DOFConstratintsExists(boundaryConditionVariable,dofConstraints,err,error,*999)
              IF(ASSOCIATED(boundaryConditionsVariable%dofConstraints)) THEN
                IF(dofConstraints%numberOfConstraints>0) THEN
                  NULLIFY(dofCoupling)
                  CALL BoundaryConditionsVariableDOFContstraints_DOFCouplingExists(dofConstraints,globalDOF,dofCoupling, &
                    & err,error,*999)
                  IF(ASSOCIATED(dofConstraint)) THEN
                    CALL SolverDofCouplings_AddCoupling(columnCouplings,dofConstraint,globalDofCouplingNumber,err,error,*999)
                  ENDIF
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
          NULLIFY(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &, &
            & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr)
          CALL SolverMapping_InterfToSolMatMapsSMInitialise(solverMapping% & &
            & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &, &
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
            & dependentVariableToSolverColMaps(interfaceMapping%numberOfInterfaceMatrices),STAT=err)
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
          DO variablePositionIdx=1,solverMapping%variablesList(solverMatrixIdx)%numberOfVariables
            IF(ASSOCIATED(lagrangeVariable,solverMapping%variablesList(solverMatrixIdx)% &
              & variables(variablePositionIdx)%variable)) THEN
              found=.TRUE.
              EXIT
            ENDIF
          ENDDO !variablePositionIdx
          IF(.NOT.found) THEN
            localError="The Lagrange variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
              & " in interface conditions index "//TRIM(NumberToVString(interfaceConditionIdx,"*",err.error))// &
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
          ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
            & lagrangeVariableToSolverColMap%columnNumbers(lagrangeVariable%totalNumberOfDofs),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate Lagrange variables to solver column maps column numbers.", &
            & err,error,*999)
          ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
            & lagrangeVariableToSolverColMap%couplingCoefficients(lagrangeVariable%totalNumberOfDofs),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate Lagrange variables to solver column maps coupling coefficients.", &
            & err,error,*999)
          ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
            & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
            & lagrangeVariableToSolverColMap%additiveConstants(lagrangeVariable%totalNumberOfDofs),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate Lagrange variables to solver column maps additive constants.", &
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
            CALL InterfceMapping_MatrixVariableGet(interfaceMapping,interfaceMatrixIdx,dependentVariable,err,error,*999)
            variableType=dependentVariable%variableType
            !Find the variable in the list of solver variables
            found=.FALSE.
            DO variablePositionIdx=1,solverMapping%variablesList(solverMatrixIdx)%numberOfVariables
              IF(ASSOCIATED(dependentVariable,solverMapping%variablesList(solverMatrixIdx)% &
                & variables(variablePositionIdx)%variable)) THEN
                found=.TRUE.
                EXIT
              ENDIF
            ENDDO !variablePositionIdx
            IF(.NOT.found) THEN
              localError="The dependent variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
                & " in interface conditions index "//TRIM(NumberToVString(interfaceConditionIdx,"*",err.error))// &
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
            ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & dependentVariableToSolverColMaps(interfaceMatrixIdx)%columnNumbers(dependentVariable%totalNumberOfDofs),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate dependent variables to solver column maps column numbers.", &
              & err,error,*999)
            ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & dependentVariableToSolverColMaps(interfaceMatrixIdx)% &
              & couplingCoefficients(dependentVariable%totalNumberOfDofs),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate dependent variables to solver column maps coupling coefficients.", &
              & err,error,*999)
            ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & dependentVariableToSolverColMaps(interfaceMatrixIdx)% &
              & additiveConstants(dependentVariable%totalNumberOfDofs),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate dependent variables to solver column maps additive constants.", &
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
                CALL BoundaryConditionsVariable_DOFConstratintsExists(boundaryConditionVariable,dofConstraints,err,error,*999)
                IF(ASSOCIATED(boundaryConditionsVariable%dofConstraints)) THEN
                  IF(dofConstraints%numberOfConstraints>0) THEN
                    NULLIFY(dofCoupling)
                    CALL BoundaryConditionsVariableDOFContstraints_DOFCouplingExists(dofConstraints,globalDOF,dofCoupling, &
                      & err,error,*999)
                    IF(ASSOCIATED(dofConstraint)) THEN
                      CALL SolverDofCouplings_AddCoupling(columnCouplings,dofConstraint,globalDofCouplingNumber,err,error,*999)
                    ENDIF
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
      DO solverVariableIdx=1,solverMapping%variablesList(solverMatrixIdx)%numberOfVariables
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
      ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
        & solverDOFToVariableDOFsMap(totalNumberOfLocalSolverDOFs),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps.",err,error,*999)
      !Set the number of columns
      solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%numberOfColumns=numberOfGlobalSolverDOFs
      !Set the number of variables
      solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%numberOfDofs=numberOfLocalSolverDOFs
      solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%totalNumberOfDofs=totalNumberOfLocalSolverDOFs
      solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%numberOfGlobalDofs=numberOfGlobalSolverDOFs
      !Allocate the columns domain mapping
      CALL DomainMapping_Initialise(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%columnDOFSMapping,err,error,*999)
      CALL DomainMapping_WorkGroupSet(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%columnDOFSMapping, &
        & workGroup,err,error,*999)
      ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%columnDomainMapping% &
        & globalToLocalMap(numberOfGlobalSolverDOFs),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate column dofs mapping global to local.",err,error,*999)
      solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%columnDomainMapping%numberOfGlobal=numberOfGlobalSolverDOFs
      ALLOCATE(variableRankProcessed(solverMapping%variablesList(solverMatrixIdx)%numberOfVariables, &
        & 0:numberOfGroupComputationNodes-1),STAT=err)
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
              & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr
            CALL SolverMapping_EquatsToSolMatMapsEMInitialise(solverMapping% &
              & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr,err,error,*999)
          ENDDO !equationMatrixIdx
          IF(ASSOCIATED(nonlinearMapping)) THEN
            !Allocate the equations set to solver matrices maps Jacobian matrix to solver matrix maps
            ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%  &
              & jacobianMatrixToSolverMatrixMaps(nonlinearMapping%numberOfResidualVariables),STAT=err)
            IF(err/=0) &
              & CALL FlagError("Could not allocate equations set to solver matrices map Jacobian matrix to solver matrix maps.", &
              & err,error,*999)
            solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%numberOfJacobianMatrices= &
              & nonlinearMapping%numberOfResidualVariables
            DO equationMatrixIdx=1,nonlinearMapping%numberOfResidualVariables
              NULLIFY(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & jacobianMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr)
            ENDDO !equationsMatrixIdx
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
              CALL SolverMapping_EquatsToSolMatMapsEMInitialise(solverMapping% &
                & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr,err,error,*999)
            ENDDO !equationMatrixIdx
          ENDIF
          IF(ASSOCIATED(nonlinearMapping)) THEN
            !Allocate the equations set to solver maps for Jacobian matrix (jm) indexing
            ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & jacobianMatrixToSolverMatrixMaps(nonlinearMapping%numberOfResidualVariables),STAT=err)
            IF(err/=0) &
              & CALL FlagError("Could not allocate equations set to solver matrices map equations matrix to solver matrix maps.", &
              & err,error,*999)
            solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%numberOfJacobianMatrices= &
              & nonlinearMapping%numberOfResidualVariables
            DO equationMatrixIdx=1,nonlinearMapping%numberOfResidualVariables
              NULLIFY(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & jacobianMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr)
            ENDDO !equationsMatrixIdx
          ENDIF
        ENDIF
        
        !Initialise solver columns to equations set map
        CALL SolverMappingSMToESMap_Initialise(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
          & solverMatrixToEquationsSetMaps(equationsSetIdx),err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        IF(ASSOCIATED(dynamicMapping)) THEN
          numberOfVariables=1
        ELSE
          IF(ASSOCIATED(nonlinearMapping)) THEN
            numberOfVariables=createValuesCache%residualVariableTypes(0,equationsSetIdx)
          ELSE
            numberOfVariables=createValuesCache%matrixVariableTypes(0,equationsSetIdx,solverMatrixIdx)
          ENDIF
        ENDIF
              
        !Allocate the solver columns to equations set map arrays
        IF(ASSOCIATED(dynamicMapping)) THEN
          solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverMatrixToEquationsSetMaps(equationsSetIdx)% &
            & haveDynamic=.TRUE.
          ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
            & solverMatrixToEquationsSetMaps(equationsSetIdx)%solverColToDynamicEquationsMaps(numberOfColumns),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate solver columns to dynamic equations map.",err,error,*999)
        ELSE
          solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverMatrixToEquationsSetMaps(equationsSetIdx)% &
            & haveStatic=.TRUE.
          ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
            & solverMatrixToEquationsSetMaps(equationsSetIdx)%solverColToStaticEquationsMaps(numberOfColumns),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate solver columns to static equations map.",err,error,*999)
        ENDIF
        !Set the solver column to equations set map
        solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
          & solverMatrixToEquationsSetMaps(equationsSetIdx)%equations=>equations
              
        !Allocate the equations matrices to solver matrix maps equations to solver maps
        IF(ASSOCIATED(dynamicMapping)) THEN
          ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%&
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
            & dynamicEquationsMatrixToSolverMatrixMaps(dynamicMapping%numberOfDynamicMatrices),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate equations matrices to solver matrix maps dynamic equations "// &
            & "to solver matrix maps.",err,error,*999)
          !Set up dynamic arrays
          solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
            & numberOfDynamicEquationsMatrices=dynamicMapping%numberOfDynamicMatrices
          DO equationMatrixIdx=1,dynamicMapping%numberOfDynamicMatrices
            NULLIFY(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%&
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)% &
              & dynamicEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr)
            ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & dynamicEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr,STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate dynamic equations matrix to solver matrix maps.",err,error,*999)
            CALL SolverMapping_EquationsToSolverMapsInitialise(solverMapping% &
              & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & dynamicEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr,err,error,*999)
            solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & dynamicEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr%equationMatrixType= &
              & SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX
            solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr%equationsMatrixNumber=equationMatrixIdx
            solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%numberOfSolverMatrices= &
              & solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr%numberOfSolverMatrices+1
          ENDDO !equationMatrixIdx       
          !Set up nonlinear arrays
          IF(ASSOCIATED(nonlinearMapping)) THEN
            solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfEquationsJacobians= &
              & nonlinearMapping%numberOfResidualVariables            
            ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & jacobianMatrixToSolverMatrixMaps(nonlinearMapping%numberOfResidualVariables),STAT=err)            
            IF(err/=0) CALL FlagError("Could not allocate Jacobian matrix to solver matrix maps.",err,error,*999)
            DO equationMatrixIdx=1,nonlinearMapping%numberOfResidualVariables
              ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & jacobianMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr,STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate Jacobian matrix to solver matrix map.",err,error,*999)
              NULLIFY(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%
                & jacobianMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr
              CALL SolverMapping_JacobianToSolverMapInitialise(solverMapping% &
                & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%
                & jacobianMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr,err,error,*999)
              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & jacobianMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr=>solverMapping% &
                & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr5 &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & jacobianMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr
            ENDDO !equationsMtrixIdx
          ENDIF
        ELSE
          IF(ASSOCIATED(linearMapping)) THEN
            ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & linearEquationsMatrixToSolverMatrixMaps(numberOfLinearMatrices),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate equations matrices to solver matrix maps linear equations"// &
              & " matrix to solver matrix maps.",err,error,*999)
            !Set up linear arrays
            solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfLinearEquationsMatrices=numberOfLinearMatrices
            DO equationMatrixIdx=1,numberOfLinearMatrices
              NULLIFY(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & linearEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr)
              ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & linearEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr,STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate linear equations matrix to solver matrix maps.",err,error,*999)
              CALL SolverMapping_EquationsToSolverMapsInitialise(solverMapping% &
                & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & linearEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr,err,error,*999)
              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%equationsMatricesToSolverMatrixMaps( &
                & solverMatrixIdx)%linearEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                & equationMatrixType=SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX
              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr%equationsMatrixNumber=equationMatrixIdx
              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr%numberOfSolverMatrices= &
                & solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr%numberOfSolverMatrices+1
            ENDDO !equationMatrixIdx
          ELSE
            ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%linearEquationsMatrixToSolverMatrixMaps(0),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate equations matrices to solver matrix maps linear equations "// &
              & "matrix to solver matrix maps.",err,error,*999)
            !Set up linear arrays`
            solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfLinearEquationsMatrices=0
          ENDIF
          !Set up nonlinear arrays
          IF(ASSOCIATED(nonlinearMapping)) THEN
            solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfEquationsJacobians= &
              & nonlinearMapping%numberOfResidualVariables
            ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & jacobianMatrixToSolverMatrixMaps(nonlinearMapping%numberOfResidualVariables),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate Jacobian matrix to solver matrix maps.",err,error,*999)
            DO equationMatrixIdx=1,nonlinearMapping%numberOfResidualVariables
              ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & jacobianMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr,STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate Jacobian matrix to solver matrix map.",err,error,*999)
              CALL SolverMapping_JacobianToSolverMapInitialise(solverMapping% &
                & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & jacobianMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr,err,error,*999)
              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & jacobianMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr=> &
                & solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & jacobianMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr
            ENDDO !equationsMatrixIdx
          ENDIF
        ENDIF
        DO variableIdx=1,numberOfVariables
          IF(ASSOCIATED(dynamicMapping)) THEN
            variableType=createValuesCache%dynamicVariableType(equationsSetIdx)
            IF(variableType==dynamicMapping%dynamicVariableType) THEN
              numberOfDynamicEquationsMatrices=dynamicMapping%numberOfDynamicsMatrices
            ELSE
              numberOfDynamicEquationsMatrices=0
            ENDIF
          ELSE
            IF(ASSOCIATED(nonlinearMapping)) THEN
              variableType=createValuesCache%residualVariableTypes(variableIdx,equationsSetIdx)
            ELSE
              variableType=reateValuesCache%matrixVariableTypes(variableIdx,equationsSetIdx,solverMatrixIdx)
              variableIndex=linearMapping%linearVariableTypeMaps(variableType)
              IF(variableIndex==0) THEN
                numberOfLinearMatrices=0
              ELSE
                numberOfLinearMatrices=linearMapping%varToEquationsMatricesMaps(variableIndex)%numberOfEquationsMatrices
              ENDIF
            ENDIF
          ENDIF

          NULLIFY(dependentVariable)
          CALL Field_VariableGet(dependentField,variableType,dependentVarible,err,error,*999)
          NULLIFY(colunDOFsMapping)
          CALL FieldVariable_DomainMappingGet(dependentVariable,columnDOFsMapping,err,error,*999)
          IF(ASSOCIATED(dynamicMapping)) THEN
            !Allocate dynamic equations to solver matrix maps equations column to solver columns maps
            DO equationMatrixIdx=1,numberOfDynamicEquationsMatrices
              matrixNumber=dynamicMapping%varToEquationsMatricesMap%equationsMatrixNumbers(equationMatrixIdx)
              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & dynamicEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr%solverMatrixNumber=solverMatrixIdx
              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & dynamicEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr%equationsMatrixNumber=matrixNumber
              solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & dynamicEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                & equationsMatrix=>dynamicMapping%equationsMatrixToVarMaps(matrixNumber)%equationsMatrix
              numberOfEquationsColumns=dynamicMapping%equationsMatrixToVarMaps(matrixNumber)%numberOfColumns
              ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & dynamicEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                & equationsColToSolverColsMap(numberOfEquationsColumns),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate dynamic equations equations matrix to solver matrix "// &
                & "column to solver columns map.",err,error,*999)
            ENDDO !equationMatrixIdx
          ELSE
            !Allocate linear equations to solver matrix maps equations column to solver columns maps
            IF(ASSOCIATED(linearMapping)) THEN
              variableIndex=linearMapping%linearVaraibelTypesMap(variableType)
              DO equationMatrixIdx=1,numberOfLinearMatrices
                matrixNumber=linearMapping%varToEquationsMatricesMaps(variableIndex)%equationsMatrixNumbers(equationMatrixIdx)
                solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                  & linearEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr%solverMatrixNumber=solverMatrixIdx
                solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr &
                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                  & linearEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr%equationsMatrixNumber=matrixNumber
                solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                  & linearEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr%equationsMatrix=> &
                  & linearMapping%equationsMatrixToVarMaps(matrixNumber)%equationsMatrix
                numberOfEquationsColumns=linearMapping%equationsMatrixToVarMaps(matrixNumber)%numberOfColumns                
                ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                  & linearEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                  & equationsColToSolverColsMap(numberOfEquationsColumns),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate linear equations matrix to solver matrix maps equations "// &
                  & "column to solver columns map.",err,error,*999)
              ENDDO !equationMatrixIdx
            ENDIF
          ENDIF
        ENDDO !variableIdx
        IF(ASSOCIATED(nonlinearMapping)) THEN
          DO equationMatrixIdx=1,nonlinearMapping%numberOfResidualVariables
            solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & jacobianMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr%solverMatrixNumber=solverMatrixIdx
            solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & jacobianMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr%jacobianMatrix=> &
              & nonlinearMapping%jacobianMatrixToVarMaps(equationMatrixIdx)%jacobian
            numberOfEquationsColumns=nonlinearMapping%jacobianMatrixToVarMaps(equationMatrixIdx)%numberOfColumns
            ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & jacobianMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
              & jacobianColToSolverColsMap(numberOfEquationsColumns),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate Jacobian matrix to solver matrix map Jacobian column "// &
              & "to solver columns map.",err,error,*999)
          ENDDO !equationsMatrixIdx
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
          CALL SolverMappingSMToICMap_Initialise(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
            & solverColToInterfaceConditionMaps(interfaceConditionIdx),err,error,*999)
          
          !Allocate the solver columns to equations set map arrays
          ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
            & solverColToInterfaceConditionMaps(interfaceConditionIdx)% &
            & solverColToInterfaceEquationsMaps(numberOfColumns),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate solver columns to interface equations map.",err,error,*999)
          
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
            ALLOCATE(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & interfaceMatrixToSolverMatrixMaps(interfaceMatrixIdx)%ptr,STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate interface matrices to solver matrix maps interface matrix "// &
              & "to solver matrix map.",err,error,*999)
            CALL SolverMapping_InterfaceToSolverMapsInitialise(solverMapping% &
              & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & interfaceMatrixToSolverMatrixMaps(interfaceMatrixIdx)%ptr,err,error,*999)
            solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr%interfaceMatrixNumber=interfaceMatrixIdx
            solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
              & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%numberOfSolverMatrices=1
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
              & interfaceMapping%interfaceMatrixRowsToVarMaps(interfaceMatrixIdx)%interfaceMatrix
            numberOfInterfaceRows=interfaceMapping%interfaceMatrixRowsToVarMaps(interfaceMatrixIdx)%numberOfRows
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
      ALLOCATE(dofMap(solverMapping%variablesList(solverMatrixIdx)%numberOfVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate dof map.",err,error,*999)
      DO solverVariableIdx=1,solverMapping%variablesList(solverMatrixIdx)%numberOfVariables
        ALLOCATE(dofMap(solverVariableIdx)%ptr(solverMapping%variablesList(solverMatrixIdx)% &
          & variables(solverVariableIdx)%variable%numberOfGlobalDofs),STAT=err)
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
          
          DO solverVariableIdx=1,solverMapping%variablesList(solverMatrixIdx)%numberOfVariables
            
            IF (solverMapping%numberOfInterfaceConditions>0) THEN
              ! Ensure that the dof_offset is calculated as a sum of the number of dofs in the diagonal entries of the solver
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
            
            variableType=solverMapping%variablesList(solverMatrixIdx)%variables(solverVariableIdx)%variableType
                  
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
                  CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equtionsSet,err,error,*999)
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
                 
                  numberOfDynamicEquationsMatrices=0
                  numberOfLinearMatrices=0
                  IF(ASSOCIATED(dynamicMapping)) THEN
                    IF(variableType==dynamicMapping%dynamicVariableType) THEN
                      numberOfDynamicEquationsMatrices=dynamicMapping%numberOfDynamicMatrices
                    ELSE
                      numberOfDynamicEquationsMatrices=0
                    ENDIF
                  ENDIF
                  IF(ASSOCIATED(linearMapping)) THEN
                    variableIndex=linearMapping%linearVariableTypesMap(variableType)
                    IF(variableIndex==0) THEN
                      numberOfLinearMatrices=0
                    ELSE
                      numberOfLinearMatrices=linearMapping% &
                        & varToEquationsMatricesMaps(variableIndex)%numberOfEquationsMatrices
                    ENDIF
                  ENDIF
                  
                  !Loop over the variables
                        
                  dependentVariable=>solverMapping%variablesList(solverMatrixIdx)%variables(solverVariableIdx)%variable
                  NULLIFY(columnDOFsMapping)
                  CALL FieldVariable_DomanMappingGet(dependentVariable,columnDOFsMapping,err,error,*999)
                        
                  DO globalDOFIdx=1,numberOfRankColumns
                    globalDOF=rankGlobalColumnsList(1,globalDOFIdx)
                    localDOF=rankGlobalColumnsList(2,globalDOFIdx)
                    !dofType=rankGlobalColumnsList(3,globalDOFIdx)
                    includeColumn=(rankGlobalColumnsList(3,globalDOFIdx)==1)
                    constrainedDOF=(rankGlobalColumnsList(3,globalDOFIdx)==2)
                    variableIdx=rankGlobalColumnsList(4,globalDOFIdx)
                    globalDofCouplingNumber=rankGlobalColumnsList(5,globalDOFIdx)
                    IF(globalDofCouplingNumber>0) THEN
                      NULLIFY(colEquationCols)
                      CALL BoundaryConditionVariableDOFContraints_DOFCouplingGet(columnCouplin,globalDOFCouplingNumber, &
                        & colEquationsCols,err,error,*999)
                      numberColEquationsCols=colEquationCols%numberOfDofs
                    ELSE
                      numberColEquationsCols=1
                      dummyDofCoupling%globalDofs(1)=globalDOF
                      dummyDofCoupling%localDofs(1)=localDOF
                      dummyDofCoupling%coefficients(1)=1.0_DP
                      colEquationCols=>dummyDofCoupling
                    END IF

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

                          !Set up the solver dofs -> variable dofs map
                          !Initialise
                          CALL SolverMappingSDOFToVDOFs_Initialise(solverMapping% &
                            & solverColToEquationsColsMap(solverMatrixIdx)% &
                            & solverDOFToVariableDOFsMap(solverLocalDOF(rank)),err,error,*999)
                          !Allocate the solver dofs to variable dofs arrays
                          ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                            & solverDOFToVariableDOFsMap(solverLocalDOF(rank))% &
                            & equationTypes(numberColEquationsCols),STAT=err)
                          IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps equations types.", &
                            & err,error,*999)
                          ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                            & solverDOFToVariableDOFsMap(solverLocalDOF(rank))% &
                            & equationIndices(numberColEquationsCols),STAT=err)
                          IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps equations indices.", &
                            & err,error,*999)
                          ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                            & solverDOFToVariableDOFsMap(solverLocalDOF(rank))%variable(numberColEquationsCols),STAT=err)
                          IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps variable type.", &
                            & err,error,*999)
                          ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                            & solverDOFToVariableDOFsMap(solverLocalDOF(rank))% &
                            & variableDOF(numberColEquationsCols),STAT=err)
                          IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps variable dof.", &
                            & err,error,*999)
                          ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                            & solverDOFToVariableDOFsMap(solverLocalDOF(rank))% &
                            & variableCoefficient(numberColEquationsCols),STAT=err)
                          IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps variable coefficient.", &
                            & err,error,*999)
                          ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                            & solverDOFToVariableDOFsMap(solverLocalDOF(rank))% &
                            & additiveConstant(numberColEquationsCols),STAT=err)
                          IF(err/=0) CALL FlagError("Could not allocate solver dof to variable maps additive constant.", &
                            & err,error,*999)
                          solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDOFToVariableDOFsMap( &
                            & solverLocalDOF(rank))%numberOfEquationDOFs=numberColEquationsCols
                          !Set the solver -> equations mappings
                          DO colEquationsColIdx=1,numberColEquationsCols
                            solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                              & solverDOFToVariableDOFsMap(solverLocalDOF(rank))%equationTypes(colEquationsColIdx)= &
                              & SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET
                            solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                              & solverDOFToVariableDOFsMap(solverLocalDOF(rank))%equationIndices(colEquationsColIdx)= &
                              & equationsSetIdx
                            solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                              & solverDOFToVariableDOFsMap(solverLocalDOF(rank))%variable(colEquationsColIdx)%ptr=> &
                              & dependentVariable
                            solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                              & solverDOFToVariableDOFsMap(solverLocalDOF(rank))%variableDOF(colEquationsColIdx)= &
                              & colEquationCols%localDofs(colEquationsColIdx)
                            solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                              & solverDOFToVariableDOFsMap(solverLocalDOF(rank))%variableCoefficient(colEquationsColIdx)= &
                              & colEquationCols%coefficients(colEquationsColIdx)
                            solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                              & solverDOFToVariableDOFsMap(solverLocalDOF(rank))%additiveConstant(colEquationsColIdx)=0.0_DP
                          ENDDO !colEquationsColIdx
                        ENDIF
                        !Set up the equations variables -> solver columns mapping
                        DO colEquationsColIdx=1,numberColEquationsCols
                          eqnLocalDof=colEquationCols%localDofs(colEquationsColIdx)
                          couplingCoefficient=colEquationCols%coefficients(colEquationsColIdx)
                          solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                            & variableToSolverColMaps(variableIdx)%columnNumbers(eqnLocalDof)=solverGlobalDOF
                          solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%pr% &
                            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                            & variableToSolverColMaps(variableIdx)%couplingCoefficients(eqnLocalDof)=couplingCoefficient
                          solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                            & variableToSolverColMaps(variableIdx)%additiveConstants(eqnLocalDof)=0.0_DP
                          !Set up the equations columns -> solver columns mapping
                          DO matrixTypeIdx=1,subMatrixList(0,equationsIdx,solverVariableIdx)
                            matrixType=subMatrixList(matrixTypeIdx,equationsIdx,solverVariableIdx)
                            SELECT CASE(matrixType)
                            CASE(SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX)
                              !Dynamic matrix
                              DO equationMatrixIdx=1,numberOfDynamicEquationsMatrices
                                matrixNumber=dynamicMapping%varToEquationsMatricesMap%equationsMatrixNumbers(equationMatrixIdx)
                                equationsColumn=dynamicMapping%varToEquationsMatricesMap%dofToColumnsMaps(equationMatrixIdx)% &
                                  & columnDOF(eqnLocalDof)
                                !Allocate the equation to solver map column items.
                                ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & dynamicEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                  & equationsColToSolverColsMap(equationsColumn)%rowCols(1),STAT=err)
                                IF(err/=0) CALL FlagError("Could not allocate dynamic equations matrix to solver matrix "//&
                                  & "equations column to solver columns map solver colums rowcols.",err,error,*999)
                                ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & dynamicEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                  & equationsColToSolverColsMap(equationsColumn)%couplingCoefficients(1),STAT=err)
                                IF(err/=0) CALL FlagError("Could not allocate dynamic equations matrix to solver matrix "//&
                                  & "equations column to solver columns map solver colums coupling coefficients.",err,error,*999)
                                solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & dynamicEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                  & equationsColToSolverColsMap(equationsColumn)%numberOfRowCols=1
                                solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & dynamicEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                  & equationsColToSolverColsMap(equationsColumn)%rowCols(1)=solverGlobalDOF
                                solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & dynamicEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                  & equationsColToSolverColsMap(equationsColumn)%couplingCoefficients(1)=couplingCoefficient
                              ENDDO !equationMatrixIdx
                            CASE(SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX)
                              variableIndex=linearMapping%linearVariableTypesMap(variableType)
                              DO equationMatrixIdx=1,numberOfLinearMatrices
                                matrixNumber=linearMapping%varToEquationsMatricesMaps(variableIndex)% &
                                  & equationsMatrixNumbers(equationMatrixIdx)
                                equationsColumn=linearMapping%varToEquationsMatricesMaps(variableIndex)% &
                                  & dofToColumnsMaps(equationMatrixIdx)%columnDOF(eqnLocalDof)
                                !Allocate the equation to solver map column items.
                                !No coupling yet so the mapping is 1-1
                                ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & linearEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                  & equationsColToSolverColsMap(equationsColumn)%rowCols(1),STAT=err)
                                IF(err/=0) CALL FlagError("Could not allocate linear equations matrix to solver matrix "//&
                                  & "equations column to solver columns map solver colums rowcols.",err,error,*999)
                                ALLOCATE(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & linearEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                  & equationsColToSolverColsMap(equationsColumn)%couplingCoefficients(1),STAT=err)
                                IF(err/=0) CALL FlagError("Could not allocate linear equations matrix to solver matrix "//&
                                  & "equations column to solver columns map solver colums coupling coefficients.",err,error,*999)
                                solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & linearEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                  & equationsColToSolverColsMap(equationsColumn)%numberOfRowCols=1
                                solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & linearEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                  & equationsColToSolverColsMap(equationsColumn)%rowCols(1)=solverGlobalDOF
                                solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & linearEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                  & equationsColToSolverColsMap(equationsColumn)%couplingCoefficients(1)=couplingCoefficient
                              ENDDO !equationMatrixIdx
                            CASE(SOLVER_MAPPING_EQUATIONS_NONLINEAR_MATRIX)
                              DO equationMatrixIdx=1,nonlinearMapping%numberOfResidualVariables
                                IF(nonlinearMapping%varToJacobianMap(equationMatrixIdx)%variableType==variableType) EXIT
                              ENDDO !equationsMatrixIdx
                              jacobianColumn=nonlinearMapping%varToJacobianMap(equationMatrixIdx)%dofToColumnsMap(eqnLocalDof)
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
                        solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%equationsMatricesToSolverMatrixMaps( &
                          & solverMatrixIdx)%variableToSolverColMaps(variableIdx)%columnNumbers(localDOF)=0
                        solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%equationsMatricesToSolverMatrixMaps( &
                          & solverMatrixIdx)%variableToSolverColMaps(variableIdx)%couplingCoefficients(localDOF)=0.0_DP
                        solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%equationsMatricesToSolverMatrixMaps( &
                          & solverMatrixIdx)%variableToSolverColMaps(variableIdx)%additiveConstants(localDOF)=0.0_DP
                        DO matrixTypeIdx=1,subMatrixList(0,equationsIdx,solverVariableIdx)
                          matrixType=subMatrixList(matrixTypeIdx,equationsIdx,solverVariableIdx)
                          SELECT CASE(matrixType)
                          CASE(SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX)
                            !Set up the equations columns -> solver columns mapping
                            DO equationMatrixIdx=1,numberOfDynamicEquationsMatrices
                              equationsColumn=dynamicMapping%varToEquationsMatricesMap% &
                                & dofToColumnsMaps(equationMatrixIdx)%columnDOF(localDOF)
                              CALL MatrixRowColCoupling_Initialise(solverMapping% &
                                & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                & dynamicEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                & equationsColToSolverColsMap(equationsColumn),err,error,*999)
                            ENDDO !equationMatrixIdx
                          CASE(SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX)
                            !Set up the equations columns -> solver columns mapping
                            variableIndex=linearMapping%linearVariableTypesMap(variableType)
                            DO equationMatrixIdx=1,numberOfLinearMatrices
                              equationsColumn=linearMapping%varToEquationsMatricesMaps(variableIndex)% &
                                & dofToColumnsMaps(equationMatrixIdx)%columnDOF(localDOF)
                              CALL MatrixRowColCoupling_Initialise(solverMapping% &
                                & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr &
                                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                & linearEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr% &
                                & equationsColToSolverColsMap(equationsColumn),err,error,*999)
                            ENDDO !equationMatrixIdx
                          CASE(SOLVER_MAPPING_EQUATIONS_NONLINEAR_MATRIX)
                            DO equationMatrixIdx=1,nonlinearMapping%numberOfResidualVariables
                              IF(nonlinearMapping%varToJacobianMap(equationMatrixIdx)%variableType==variableType) EXIT
                            ENDDO !equationsMatrixIdx
                            jacobianColumn=nonlinearMapping%varToJacobianMap(equationMatrixIdx)%dofToColumnsMap(localDOF)
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
                    
                    lagrangeVariable=>solverMapping%variablesList(solverMatrixIdx)%variables(solverVariableIdx)%variable
                         
                    DO globalDOFIdx=1,numberOfRankColumns
                      globalDOF=rankGlobalColumnsList(1,globalDOFIdx)
                      localDOF=rankGlobalColumnsList(2,globalDOFIdx)
                      !dofType=rankGlobalColumnsList(3,globalDOFIdx)
                      includeColumn=(rankGlobalColumnsList(3,globalDOFIdx)==1)
                      constrainedDOF=(rankGlobalColumnsList(3,globalDOFIdx)==2)
                      globalDofCouplingNumber=rankGlobalColumnsList(5,globalDOFIdx)
                      
                      IF(globalDofCouplingNumber>0) THEN
                        NULLIFY(colEquationsCols)
                        CALL BoundaryConditionVariableDOFConstraints_DOFCouplingGet(columnCouplings,globalDOFCouplingNumber, &
                          & colEquationsCols,err,error,*999)
                        numberColEquationsCols=colEquationCols%numberOfDofs
                      ELSE
                        numberColEquationsCols=1
                        dummyDofCoupling%globalDofs(1)=globalDOF
                        dummyDofCoupling%localDofs(1)=localDOF
                        dummyDofCoupling%coefficients(1)=1.0_DP
                        colEquationCols=>dummyDofCoupling
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
                            CALL SolverMappingSDOFToVDOFs_Initialise(solverMapping% &
                              & solverColToEquationsColsMap(solverMatrixIdx)%solverDOFToVariableDOFsMap(solverLocalDOF(rank)), &
                              & err,error,*999)
                            !Allocate the solver dofs to variable dofs arrays
                            ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                              & solverDOFToVariableDOFsMap(solverLocalDOF(rank))%equationTypes(numberColEquationsCols),STAT=err)
                            IF(err/=0) CALL FlagError("Could not allocate solver DOF to variable maps equations types.", &
                              & err,error,*999)
                            ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                              & solverDOFToVariableDOFsMap(solverLocalDOF(rank))%equationIndices(numberColEquationsCols),STAT=err)
                            IF(err/=0) CALL FlagError("Could not allocate solver DOF to variable maps equations indices.", &
                              & err,error,*999)
                            ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                              & solverDOFToVariableDOFsMap(solverLocalDOF(rank))%variable(numberColEquationsCols),STAT=err)
                            IF(err/=0) CALL FlagError("Could not allocate solver DOF to variable maps variable type.", &
                              & err,error,*999)
                            ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                              & solverDOFToVariableDOFsMap(solverLocalDOF(rank))%variableDOF(numberColEquationsCols),STAT=err)
                            IF(err/=0) CALL FlagError("Could not allocate solver DOF to variable maps variable DOF.", &
                              & err,error,*999)
                            ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                              & solverDOFToVariableDOFsMap(solverLocalDOF(rank))%variableCoefficient(numberColEquationsCols), &
                              & STAT=err)
                            IF(err/=0) CALL FlagError("Could not allocate solver DOF to variable maps variable coefficient.", &
                              & err,error,*999)
                            ALLOCATE(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                              & solverDOFToVariableDOFsMap(solverLocalDOF(rank))%additiveConstant(numberColEquationsCols),STAT=err)
                            IF(err/=0) CALL FlagError("Could not allocate solver DOF to variable maps additive constant.", &
                              & err,error,*999)
                            !Setup
                            solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDOFToVariableDOFsMap( &
                              & solverLocalDOF(rank))%numberOfEquationDOFs=numberColEquationsCols
                            DO colEquationsColIdx=1,numberColEquationsCols
                              eqnLocalDof=colEquationCols%localDofs(colEquationsColIdx)
                              couplingCoefficient=colEquationCols%coefficients(colEquationsColIdx)
                              solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDOFToVariableDOFsMap(&
                                & solverLocalDOF(rank))%equationTypes(colEquationsColIdx)= &
                                & SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION
                              solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDOFToVariableDOFsMap(&
                                & solverLocalDOF(rank))%equationIndices(colEquationsColIdx)=interfaceConditionIdx
                              solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDOFToVariableDOFsMap(&
                                & solverLocalDOF(rank))%variable(colEquationsColIdx)%ptr=>lagrangeVariable
                              solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDOFToVariableDOFsMap(&
                                & solverLocalDOF(rank))%variableDOF(colEquationsColIdx)=eqnLocalDof
                              solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDOFToVariableDOFsMap(&
                                & solverLocalDOF(rank))%variableCoefficient(colEquationsColIdx)=couplingCoefficient
                              solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDOFToVariableDOFsMap(&
                                & solverLocalDOF(rank))%additiveConstant(1)=0.0_DP
                            ENDDO !colEquationsColIdx
                          ENDIF
                          DO colEquationsColIdx=1,numberColEquationsCols
                            eqnLocalDof=colEquationCols%localDofs(colEquationsColIdx)
                            couplingCoefficient=colEquationCols%coefficients(colEquationsColIdx)
                            IF(equationType==SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION) THEN
                              IF(.NOT.variableRankProcessed(solverVariableIdx,rank)) THEN
                                !Set up the equations variables -> solver columns mapping
                                !No coupling yet so the mapping is 1-1
                                solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                                  & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%&
                                  & lagrangeVariableToSolverColMap%columnNumbers(eqnLocalDof)=solverGlobalDOF
                                solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                                  & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & lagrangeVariableToSolverColMap%couplingCoefficients(eqnLocalDof)=couplingCoefficient
                                solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                                  & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                  & lagrangeVariableToSolverColMap%additiveConstants(eqnLocalDof)=0.0_DP
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
                                & dependentVariableToSolverColMaps(interfaceMatrixIdx)%columnNumbers(eqnLocalDof)=solverGlobalDOF
                              solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                                & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                & dependentVariableToSolverColMaps(interfaceMatrixIdx)%couplingCoefficients(eqnLocalDof)= &
                                & couplingCoefficient
                              solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                                & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                                & dependentVariableToSolverColMaps(interfaceMatrixIdx)%additiveConstants(eqnLocalDof)=0.0_DP
                              !Set up the equations columns -> solver columns mapping
                              interfaceRow=interfaceMapping%interfaceMatrixRowsToVarMaps(interfaceMatrixIdx)% &
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
                                & interfaceMapping%interfaceMatrixRowsToVarMaps(interfaceMatrixIdx)%matrixCoefficient
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
                              & lagrangeVariableToSolverColMap%columnNumbers(localDOF)=0
                            solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                              & lagrangeVariableToSolverColMap%couplingCoefficients(localDOF)=0.0_DP
                            solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                              & lagrangeVariableToSolverColMap%additiveConstants(localDOF)=0.0_DP
                            interfaceColumn=interfaceMapping%lagrangeDOFToColumnMap(localDOF)
                            CALL MatrixRowColCoupling_Initialise(solverMapping% &
                              & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                              & interfaceColToSolverColsMap(interfaceColumn),err,error,*999)
                          ELSE
                            !Set up the equations variables -> solver columns mapping
                            solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                              & dependentVariableToSolverColMaps(interfaceMatrixIdx)%columnNumbers(localDOF)=0
                            solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                              & dependentVariableToSolverColMaps(interfaceMatrixIdx)%couplingCoefficients(localDOF)=0.0_DP
                            solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
                              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                              & dependentVariableToSolverColMaps(interfaceMatrixIdx)%additiveConstants(localDOF)=0.0_DP
                            !Set up the equations columns -> solver columns mapping
                            interfaceRow=interfaceMapping%interfaceMatrixRowsToVarMaps(interfaceMatrixIdx)% &
                              & variableDOFToRowMap(localDOF)
                            CALL MatrixRowColCoupling_Initialise(solverMapping% &
                              & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr%&
                              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                              & interfaceMatrixToSolverMatrixMaps(interfaceMatrixIdx)%ptr% 
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
      DO solverVariableIdx=1,solverMapping%variablesList(solverMatrixIdx)%numberOfVariables
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
      CALL SolverMapping_EquationsSetGet(equationsSet,equationsSetIdx,equationsSet,err,error,*999)
      NULLIFY(equations0
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
            & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr% & &
            & equationsMatrixToSolverMatrixMaps(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)% &
            & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%numberOfSolverMatrices),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate equations matrix to solver matrix maps.",err,error,*999)
          solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%ptr%numberOfSolverMatrices=0
        ENDDO !equationMatrixIdx
        DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
          DO equationMatrixIdx=1,solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfDynamicEquationsMatrices
            IF(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & dynamicEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr%solverMatrixNumber==solverMatrixIdx) THEN
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
                & dynamicEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr
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
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfLinearEquationsMatrices
            IF(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
              & linearEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr%solverMatrixNumber==solverMatrixIdx) THEN
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
                & linearEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr
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
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of solver matrices = ",solverMapping%numberOfSolverMatrices, &
        & err,error,*999)               
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Equation sets:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of equations sets = ",solverMapping%numberOfEquationsSets, &
        & err,error,*999)
      DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Equations set index : ",equationsSetIdx,err,error,*999)
        NULLIFY(region)
        CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Region user number        = ",region%userNumber,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Equations set user number = ",equationsSet%userNumber,err,error,*999)
      ENDDO !equationsSetIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Interface conditions:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of interface conditions = ",solverMapping% &
        & numberOfInterfaceConditions,err,error,*999)
      DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
        NULLIFY(interfaceCondition)
        CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Interface condition index : ",interfaceConditionIdx, &
          & err,error,*999)
        NULLIFY(INTERFACE)
        CALL InterfaceCondition_InterfaceGet(interfaceCondition,INTERFACE,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Interface user number           = ",interface%userNumber, &
          & err,error,*999)
        NULLIFY(parentRegion)
        CALL Interface_ParentRegionGet(INTERFACE,parentRegion,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Parent region user number       = ",parentRegion%userNumber, &
          & err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Interface condition user number = ",interfaceCondition%userNumber, &
          & err,error,*999)
      ENDDO !equationsSetIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Equations variables list:",err,error,*999)
      DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Solver matrix : ",solverMatrixIdx,err,error,*999)        
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of variables = ",solverMapping%variablesList(solverMatrixIdx)% &
          & numberOfVariables,err,error,*999)
        DO variableIdx=1,solverMapping%variablesList(solverMatrixIdx)%numberOfVariables
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Variable : ",variableIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Variable type = ",solverMapping%variablesList(solverMatrixIdx)% &
            & variables(variableIdx)%variableType,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of equations = ",solverMapping% &
            & variablesList(solverMatrixIdx)%variables(variableIdx)%numberOfEquations,err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%variablesList(solverMatrixIdx)%variables(variableIdx)% &
            & numberOfEquations,5,5,solverMapping%variablesList(solverMatrixIdx)%variables(variableIdx)%equationTypes, &
            & '("        Equation types   :",5(X,I13))','(26X,5(X,I13))',err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%variablesList(solverMatrixIdx)%variables(variableIdx)% &
            & numberOfEquations,5,5,solverMapping%variablesList(solverMatrixIdx)%variables(variableIdx)%equationTypes, &
            & '("        Equation indices :",5(X,I13))','(26X,5(X,I13))',err,error,*999)
        ENDDO !variableIdx
      ENDDO !solverMatrixIdx
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    RHS vector : ",solverMatrixIdx,err,error,*999)        
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of variables = ",solverMapping%rhsVariablesList% &
        & numberOfVariables,err,error,*999)
      DO variableIdx=1,solverMapping%rhsVariablesList%numberOfVariables
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Variable : ",variableIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Variable type = ",solverMapping%rhsVariablesList% &
          & variables(variableIdx)%variableType,err,error,*999)        
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of equations = ",solverMapping%rhsVariablesList% &
          & variables(variableIdx)%numberOfEquations,err,error,*999)        
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%rhsVariablesList%variables(variableIdx)% &
          & numberOfEquations,5,5,solverMapping%rhsVariablesList%variables(variableIdx)%equationTypes, &
          & '("        Equation types   :",5(X,I13))','(26X,5(X,I13))',err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%rhsVariablesList%variables(variableIdx)% &
          & numberOfEquations,5,5,solverMapping%rhsVariablesList%variables(variableIdx)%equationTypes, &
          & '("        Equation indices :",5(X,I13))','(26X,5(X,I13))',err,error,*999)
      ENDDO !variableIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Solver row to equations rows mappings:",err,error,*999)      
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of rows = ",solverMapping%numberOfRows,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of global rows = ",solverMapping%numberOfGlobalRows, &
        & err,error,*999)
      IF(diagnostics2) THEN
        DO rowIdx=1,solverMapping%numberOfRows
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Solver row : ",rowIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of equations set rows mapped to = ",solverMapping% &
            & solverRowToEquationsRowsMap(rowIdx)%numberOfEquationsSetRows,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Interface condition index          = ",solverMapping% &
            & solverRowToEquationsRowsMap(rowIdx)%interfaceConditionIndex,err,error,*999)
          IF(solverMapping%solverRowToEquationsRowsMap(rowIdx)%interfaceConditionIndex==0) THEN
            !Row is an equations set row
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverRowToEquationsRowsMap(rowIdx)% &
              & numberOfEquationsSetRows,5,5,solverMapping%solverRowToEquationsRowsMap(rowIdx)%equationsIndex, &
              & '("      Equations indices      :",5(X,I13))','(30X,5(X,I13))',err,error,*999) 
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverRowToEquationsRowsMap(rowIdx)% &
              & numberOfEquationsSetRows,5,5,solverMapping%solverRowToEquationsRowsMap(rowIdx)%rowColNumber, &
              & '("      Equations row numbers  :",5(X,I13))','(30X,5(X,I13))',err,error,*999) 
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverRowToEquationsRowsMap(rowIdx)% &
              & numberOfEquationsSetRows,5,5,solverMapping%solverRowToEquationsRowsMap(rowIdx)%couplingCoefficients, &
              & '("      Coupling coefficients  :",5(X,E13.6))','(30X,5(X,E13.6))',err,error,*999)
          ELSE
            !Row is an interface condition row
!!TODO: format better
            CALL WriteStringFmtValue(DIAGNOSTIC_OUTPUT_TYPE,"        Interface col numbers : ",solverMapping% &
              & solverRowToEquationsRowsMap(rowIdx)%rowColNumber(1),"(I13)",err,error,*999) 
            CALL WriteStringFmtValue(DIAGNOSTIC_OUTPUT_TYPE,"        Coupling coefficients : ",solverMapping% &
              & solverRowToEquationsRowsMap(rowIdx)%couplingCoefficients(1),"(E13.6)",err,error,*999)
          ENDIF
        ENDDO !rowIdx
      ENDIF
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Solver column to equations column mappings:",err,error,*999)      
      DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Solver matrix : ",solverMatrixIdx,err,error,*999)
        NULLIFY(solverColToEquationsColsMap)
        CALL SolverMapping_SolverColToEquationsColMapGet(solverMapping,solverMatrixIdx,solverColToEquationsColMap,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of columns = ",solverColToEquationsColsMap%numberOfColumns, &
          & err,error,*999)
        IF(diagnostics2) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      Solver column to equations set columns mappings:",err,error,*999)
          DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Equations set index : ",equationsSetIdx,err,error,*999)
            DO columnIdx=1,solverColToEquationsColsMap%numberOfColumns           
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE, "          Solver column : ",columnIdx,err,error,*999)
              NULLIFY(solverColToEquationsSetMap)
              CALL SolverMappingSCToECSMap_SolverColToEquationsSetMapGet(solverColToEquationsSetMap,equationsSetIdx, &
                & solverColToEquationsSetMap,err,error,*999)
              IF(solverColToEquationsSetMap%haveDynamic) THEN
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of dynamic equations matrices mapped to = ", &
                  & solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverMatrixToEquationsSetMaps(equationsSetIdx)% &
                  & solverColToDynamicEquationsMaps(columnIdx)%numberOfDynamicEquationsMatrices,err,error,*999)
                IF(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverMatrixToEquationsSetMaps(equationsSetIdx)% &
                  & solverColToDynamicEquationsMaps(columnIdx)%numberOfDynamicEquationsMatrices>0) THEN
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                    & solverMatrixToEquationsSetMaps(equationsSetIdx)%solverColToDynamicEquationsMaps(columnIdx)% &
                    & numberOfDynamicEquationsMatrices,5,5,solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                    & solverMatrixToEquationsSetMaps(equationsSetIdx)%solverColToDynamicEquationsMaps(columnIdx)% &
                    & equationsMatrixNumbers,'("            Equation matrices numbers :",5(X,I13))','(39X,5(X,I13))',err,error,*999)
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                    & solverMatrixToEquationsSetMaps(equationsSetIdx)%solverColToDynamicEquationsMaps(columnIdx)% &
                    & numberOfDynamicEquationsMatrices,5,5,solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                    & solverMatrixToEquationsSetMaps(equationsSetIdx)%solverColToDynamicEquationsMaps(columnIdx)% &
                    & equationsColumnNumbers,'("            Equation column numbers   :",5(X,I13))','(39X,5(X,I13))',err,error,*999)
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                    & solverMatrixToEquationsSetMaps(equationsSetIdx)%solverColToDynamicEquationsMaps(columnIdx)% &
                    & numberOfDynamicEquationsMatrices,5,5,solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
                    & solverMatrixToEquationsSetMaps(equationsSetIdx)%solverColToDynamicEquationsMaps(columnIdx)% &
                    & couplingCoefficients,'("            Coupling coefficients     :",5(X,E13.6))','(39X,5(X,E13.6))', &
                    & err,error,*999)
                ENDIF
              ELSE
                 CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of dynamic equations matrices mapped to = ", &
                  & 0_INTG,err,error,*999)
              ENDIF
!!TODO what about dynamic nonlinear mappings???
              IF(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverMatrixToEquationsSetMaps(equationsSetIdx)% &
                & haveStatic) THEN
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE, "            Number of linear equations matrices mapped to  = ", &
                  & solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverMatrixToEquationsSetMaps(equationsSetIdx)% &
                  & solverColToStaticEquationsMaps(columnIdx)%numberOfLinearEquationsMatrices,err,error,*999)
                IF(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverMatrixToEquationsSetMaps(equationsSetIdx)% &
                  & solverColToStaticEquationsMaps(columnIdx)%numberOfLinearEquationsMatrices>0) THEN
                  !CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap( &
                  !  & solverMatrixIdx)%solverMatrixToEquationsSetMaps(equationsSetIdx)%solverColToStaticEquationsMaps( &
                  !  & columnIdx)%numberOfLinearMatrices,5,5,solverMapping%solverColToEquationsColsMap( &
                  !  & solverMatrixIdx)%solverMatrixToEquationsSetMaps(equationsSetIdx)%solverColToStaticEquationsMaps( &
                  !  & columnIdx)%equationsMatrixNumbers,'("            Equation matrices numbers :",5(X,I13))', &
                  !  & '(36X,5(X,I13))',err,error,*999)
                  !CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap( &
                  !  & solverMatrixIdx)%solverMatrixToEquationsSetMaps(equationsSetIdx)%solverColToStaticEquationsMaps( &
                  !  & columnIdx)%numberOfLinearMatrices,5,5,solverMapping%solverColToEquationsColsMap( &
                  !  & solverMatrixIdx)%solverMatrixToEquationsSetMaps(equationsSetIdx)%solverColToStaticEquationsMaps( &
                  !  & columnIdx)%equationsColumnNumbers,'("            Equation column numbers   :",5(X,I13))', &
                  !  & '(36X,5(X,I13))',err,error,*999)
                  !CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap( &
                  !  & solverMatrixIdx)%solverMatrixToEquationsSetMaps(equationsSetIdx)%solverColToStaticEquationsMaps( &
                  !  & columnIdx)%numberOfLinearMatrices,5,5,solverMapping%solverColToEquationsColsMap( &
                  !  & solverMatrixIdx)%solverMatrixToEquationsSetMaps(equationsSetIdx)%solverColToStaticEquationsMaps( &
                  !  & columnIdx)%couplingCoefficients,'("            Coupling coefficients     :",5(X,E13.6))', &
                  !  & '(36X,5(X,E13.6))',err,error,*999)
                ENDIF
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Jacobian column number     = ", &
                  & solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverMatrixToEquationsSetMaps(equationsSetIdx)% &
                  & solverColToStaticEquationsMaps(columnIdx)%jacobianColumnNumber,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Jacobian coupling coeff    = ", &
                  & solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverMatrixToEquationsSetMaps(equationsSetIdx)% &
                  & solverColToStaticEquationsMaps(columnIdx)%jacobianCouplingCoefficient,err,error,*999)
              ELSE
                 CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of static equations matrices mapped to  = ", &
                  & 0_INTG,err,error,*999)
              ENDIF
            ENDDO !columnIdx
          ENDDO !equationsSetIdx
        ENDIF
      ENDDO !solverMatrixIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Solver DOF to field DOFs mappings:",err,error,*999)
      DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Solver matrix : ",solverMatrixIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of DOFs = ",solverMapping% &
          & solverColToEquationsColsMap(solverMatrixIdx)%numberOfDofs,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Total number of DOFs = ",solverMapping% &
          & solverColToEquationsColsMap(solverMatrixIdx)%totalNumberOfDofs,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of global DOFs = ",solverMapping% &
          & solverColToEquationsColsMap(solverMatrixIdx)%numberOfGlobalDofs,err,error,*999)
        ALLOCATE(variableTypes(solverMapping%numberOfEquationsSets+solverMapping%numberOfInterfaceConditions),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate variable types.",err,error,*999)
        DO dofIdx=1,solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%totalNumberOfDofs     
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Solver local DOF : ",dofIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of equations DOFs mapped to     = ",solverMapping% &
            & solverColToEquationsColsMap(solverMatrixIdx)%solverDOFToVariableDOFsMap(dofIdx)%numberOfEquationDOFs,err,error,*999)
          IF(solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDOFToVariableDOFsMap(dofIdx)% &
            & numberOfEquationDOFs>0) THEN
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
              & solverDOFToVariableDOFsMap(dofIdx)%numberOfEquationDOFs,5,5,solverMapping% &
              & solverColToEquationsColsMap(solverMatrixIdx)%solverDOFToVariableDOFsMap(dofIdx)%equationIndices, &
              & '("        Equations types       :",5(X,I13))','(31X,5(X,I13))',err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
              & solverDOFToVariableDOFsMap(dofIdx)%numberOfEquationDOFs,5,5,solverMapping% &
              & solverColToEquationsColsMap(solverMatrixIdx)%solverDOFToVariableDOFsMap(dofIdx)%equationIndices, &
              & '("        Equations indices     :",5(X,I13))','(31X,5(X,I13))',err,error,*999)
            DO variableIdx=1,solverMapping%solverColToEquationsColsMap(solverMatrixIdx)%solverDOFToVariableDOFsMap(dofIdx)% &
              & numberOfEquationDOFs
              variableTypes(variableIdx)=solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
              & solverDOFToVariableDOFsMap(dofIdx)%variable(variableIdx)%ptr%variableType
            ENDDO !variableIdx
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
              & solverDOFToVariableDOFsMap(dofIdx)%numberOfEquationDOFs,5,5,variableTypes, &
              & '("        Variable types        :",5(X,I13))','(31X,5(X,I13))',err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
              & solverDOFToVariableDOFsMap(dofIdx)%numberOfEquationDOFs,5,5,solverMapping% &
              & solverColToEquationsColsMap(solverMatrixIdx)%solverDOFToVariableDOFsMap(dofIdx)%variableDOF, &
              & '("        Variable DOFs         :",5(X,I13))','(31X,5(X,I13))',err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
              & solverDOFToVariableDOFsMap(dofIdx)%numberOfEquationDOFs,5,5,solverMapping% &
              & solverColToEquationsColsMap(solverMatrixIdx)%solverDOFToVariableDOFsMap(dofIdx)% &
              & variableCoefficient, & 
              & '("        Variable coefficients :",5(X,E13.6))','(31X,5(X,E13.6))',err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%solverColToEquationsColsMap(solverMatrixIdx)% &
              & solverDOFToVariableDOFsMap(dofIdx)%numberOfEquationDOFs,5,5,solverMapping% &
              & solverColToEquationsColsMap(solverMatrixIdx)%solverDOFToVariableDOFsMap(dofIdx)%additiveConstant, &
              & '("        Additive constants    :",5(X,E13.6))','(31X,5(X,E13.6))',err,error,*999)
          ENDIF
        ENDDO !dofIdx
        IF(ALLOCATED(variableTypes)) DEALLOCATE(variableTypes)
      ENDDO !solverMatrixIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Equation sets to solver mappings:",err,error,*999)



      
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
        NULLIFY(rhsMapping)
        CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
        NULLIFY(sourceMapping)
        CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Equations set index : ",equationsSetIdx,err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      Equations sets rows to solver rows mappings:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE, "        Number of equations set rows = ",vectorMapping% &
         & totalNumberOfRows,err,error,*999)
        DO rowIdx=1,vectorMapping%totalNumberOfRows
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Equations set row : ",rowIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Number of solver rows mapped to   = ",solverMapping% &
            & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%equationsRowToSolverRowsMap(rowIdx)%numberOfRowCols, &
            & err,error,*999)
          IF(solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%equationsRowToSolverRowsMap(rowIdx)% &
            & numberOfRowCols>0) THEN
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping% &
              & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%equationsRowToSolverRowsMap(rowIdx)% &
              & numberOfRowCols,5,5,solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsRowToSolverRowsMap(rowIdx)%rowCols,'("          Solver row numbers    :",5(X,I13))','(33X,5(X,I13))', &
              & err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping% &
              & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%equationsRowToSolverRowsMap(rowIdx)% &
              & numberOfRowCols,5,5,solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsRowToSolverRowsMap(rowIdx)%couplingCoefficients, &
              & '("          Coupling coefficients :",5(X,E13.6))','(33X,5(X,E13.6))',err,error,*999)
          ENDIF
        ENDDO !rowIdx
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      Solver matrix indexing:",err,error,*999)
        DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Solver matrix : ",solverMatrixIdx,err,error,*999)
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"          Interface conditions affecting:",err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of interface conditions = ",solverMapping% &
            & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%numberOfInterfaceConditions,err,error,*999)
          DO interfaceConditionIdx=1,solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
            & numberOfInterfaceConditions
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Interface condition : ",interfaceConditionIdx, &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Interface condition index = ",solverMapping% &
              & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%interfaceConditionToEquationsSetMaps( &
              & interfaceConditionIdx)%interfaceConditionIndex,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Interface matrix number   = ",solverMapping% &
              & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%interfaceConditionToEquationsSetMaps( &
              & interfaceConditionIdx)%interfaceMatrixNumber,err,error,*999)
          ENDDO !interfaceConditionIdx                    
          IF(ASSOCIATED(dynamicMapping)) THEN
           CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"          Dynamic equations matrix columns to solver matrix columns:", &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of dynamic equations matrices = ",solverMapping% &
              & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfDynamicEquationsMatrices,err,error,*999)
            DO equationMatrixIdx=1,solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfDynamicEquationsMatrices
              equationsMatrixToSolverMatrixMap=>solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                & dynamicEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Equations matrix index : ",equationMatrixIdx, &
                & err,error,*999)
              equationsMatrixNumber=equationsMatrixToSolverMatrixMap%equationsMatrixNumber
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Equations matrix number = ",equationsMatrixNumber, &
                & err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Solver matrix number    = ", &
                & equationsMatrixToSolverMatrixMap%solverMatrixNumber,err,error,*999)
              DO columnIdx=1,dynamicMapping%equationsMatrixToVarMaps(equationsMatrixNumber)%numberOfColumns
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
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of liner equations matrices = ",solverMapping% &
                & equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfLinearEquationsMatrices,err,error,*999)
              DO equationMatrixIdx=1,solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfLinearEquationsMatrices
                equationsMatrixToSolverMatrixMap=>solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
                  & linearEquationsMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Equations matrix index : ",equationMatrixIdx, &
                  & err,error,*999)
                equationsMatrixNumber=equationsMatrixToSolverMatrixMap%equationsMatrixNumber
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Equations matrix number = ",equationsMatrixNumber, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Solver matrix number    = ", &
                  & equationsMatrixToSolverMatrixMap%solverMatrixNumber,err,error,*999)
                DO columnIdx=1,linearMapping%equationsMatrixToVarMaps(equationsMatrixNumber)%numberOfColumns
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
              DO equationMatrixIdx=1,nonlinearMapping%numberOfResidualVariables
                jacobianMatrixToSolverMatrixMap=>solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)% &
                  & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%jacobianMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Solver matrix number    = ", &
                  & jacobianMatrixToSolverMatrixMap%solverMatrixNumber,err,error,*999)
                DO columnIdx=1,nonlinearMapping%jacobianMatrixToVarMaps(equationMatrixIdx)%numberOfColumns
                  CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Equations matrix column : ",columnIdx,err,error,*999)
                  CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Number of solver columns mapped to = ", &
                    & jacobianMatrixToSolverMatrixMap%jacobianColToSolverColsMap(columnIdx)%numberOfRowCols,err,error,*999)
                  IF(jacobianMatrixToSolverMatrixMap%jacobianColToSolverColsMap(columnIdx)%numberOfRowCols>0) THEN
                    CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,jacobianMatrixToSolverMatrixMap%jacobianColToSolverColsMap( &
                      & columnIdx)%numberOfRowCols,5,5,jacobianMatrixToSolverMatrixMap%jacobianColToSolverColsMap(columnIdx)% &
                      & rowCols,'("              Solver columns         :",5(X,I13))','(38X,5(X,I13))',err,error,*999)
                    CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,jacobianMatrixToSolverMatrixMap%jacobianColToSolverColsMap( &
                      & columnIdx)%numberOfRowCols,5,5,jacobianMatrixToSolverMatrixMap%jacobianColToSolverColsMap(columnIdx)% &
                      & couplingCoefficients,'("              Coupling coefficients  :",5(X,E13.6))','(38X,5(X,E13.6))', &
                      & err,error,*999)
                  ENDIF
                ENDDO !columnIdx
              ENDDO !equationMatrixIdx
            ENDIF
          ENDIF
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"          Variable DOFs to solver matrix DOFs:",err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of variables = ",solverMapping% &
            & equationsSetToSolverMatricesMaps(equationsSetIdx)%equationsMatricesToSolverMatrixMaps(solverMatrixIdx)% &
            & numberOfVariables,err,error,*999) 
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)% &
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%numberOfVariables,5,5,solverMapping% &
            & equationsSetToSolverMatricesMaps(equationsSetIdx)%equationsMatricesToSolverMatrixMaps(solverMatrixIdx)% &
            & variableTypes,'("            Variable types :",5(X,I13))','(28X,5(X,I13))',err,error,*999)
          DO variableIdx=1,solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)% &
            & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%numberOfVariables
            dependentVariable=>solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)% &
              & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%variables(variableIdx)%ptr
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Variable index : ",variableIdx,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Number of variable DOFs = ",dependentVariable% &
              & numberOfDofs,err,error,*999)
            DO localDOF=1,dependentVariable%totalNumberOfDofs
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Variable DOF : ",localDOF,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Solver column number = ",solverMapping% &
                & equationsSetToSolverMatricesMaps(equationsSetIdx)%equationsMatricesToSolverMatrixMaps(solverMatrixIdx)% &
                & variableToSolverColMaps(variableIdx)%columnNumbers(localDOF),err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Coupling coefficient = ",solverMapping% &
                & equationsSetToSolverMatricesMaps(equationsSetIdx)%equationsMatricesToSolverMatrixMaps(solverMatrixIdx)% &
                & variableToSolverColMaps(variableIdx)%couplingCoefficients(localDOF),err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Additive constant    = ",solverMapping% &
                & equationsSetToSolverMatricesMaps(equationsSetIdx)%equationsMatricesToSolverMatrixMaps(solverMatrixIdx)% &
                & variableToSolverColMaps(variableIdx)%additiveConstants(localDOF),err,error,*999)              
            ENDDO !localDOF
          ENDDO !variableIdx
        ENDDO !solverMatrixIdx
        IF(ASSOCIATED(dynamicMapping)) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      Dynamic equations matrix indexing:",err,error,*999)
          DO equationMatrixIdx=1,dynamicMapping%numberOfDynamicMatrices
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Equations matrix : ",equationMatrixIdx,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Number of solver matrices = ",solverMapping% &
              & equationsSetToSolverMatricesMaps(equationsSetIdx)%equationsMatrixToSolverMatricesMaps(equationMatrixIdx)% &
              & numberOfSolverMatrices,err,error,*999)
            DO solverMatrixIdx=1,solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)% &
              & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%numberOfSolverMatrices
              equationsMatrixToSolverMatrixMap=>solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)% &
                & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%equationsMatrixToSolverMatrixMaps( &
                & solverMatrixIdx)%ptr
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Solver matrix index : ",solverMatrixIdx,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Equations matrix number = ", &
                & equationsMatrixToSolverMatrixMap%equationsMatrixNumber,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Solver matrix number    = ", &
                & equationsMatrixToSolverMatrixMap%solverMatrixNumber,err,error,*999)            
            ENDDO !solverMatrixIdx
          ENDDO !equationMatrixIdx
        ELSE
          IF(ASSOCIATED(linearMapping)) THEN
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      Linear equations matrix indexing:",err,error,*999)
            DO equationMatrixIdx=1,linearMapping%numberOfLinearMatrices
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Equations matrix : ",equationMatrixIdx,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Number of solver matrices = ",solverMapping% &
                & equationsSetToSolverMatricesMaps(equationsSetIdx)%equationsMatrixToSolverMatricesMaps(equationMatrixIdx)% &
                & numberOfSolverMatrices,err,error,*999)
              DO solverMatrixIdx=1,solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)% &
                & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%numberOfSolverMatrices
                equationsMatrixToSolverMatrixMap=>solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)% &
                  & equationsMatrixToSolverMatricesMaps(equationMatrixIdx)%equationsMatrixToSolverMatrixMaps( &
                  & solverMatrixIdx)%ptr
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Solver matrix index : ",solverMatrixIdx,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Equations matrix number = ", &
                  & equationsMatrixToSolverMatrixMap%equationsMatrixNumber,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Solver matrix number    = ", &
                  & equationsMatrixToSolverMatrixMap%solverMatrixNumber,err,error,*999)            
              ENDDO !solverMatrixIdx
            ENDDO !equationMatrixIdx
          ENDIF
          IF(ASSOCIATED(nonlinearMapping)) THEN
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      Jacobian matrix indexing:",err,error,*999)
            DO equationMatrixIdx=1,nonlinearMapping%numberOfResidualVariables
              jacobianMatrixToSolverMatrixMap=>solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)% &
                & jacobianMatrixToSolverMatrixMaps(equationMatrixIdx)%ptr
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Solver matrix number    = ", &
                & jacobianMatrixToSolverMatrixMap%solverMatrixNumber,err,error,*999)
            ENDDO
          ENDIF
        ENDIF
      ENDDO !equationsSetIdx
      IF(solverMapping%numberOfInterfaceConditions>0) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Interface conditions to solver mappings:",err,error,*999)
        DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
          interfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr
          interfaceEquations=>interfaceCondition%interfaceEquations
          interfaceMapping=>interfaceEquations%interfaceMapping
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Interface to equations sets mapping:",err,error,*999)
          DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Solver matrix : ",solverMatrixIdx,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of equations sets = ",solverMapping% &
              & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%numberOfEquationsSets,err,error,*999)
            DO equationsSetIdx=1,solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)% &
              & numberOfEquationsSets
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Equations set : ",equationsSetIdx,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Equations set index     = ",solverMapping% &
                & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%equationsSetToInterfaceConditionMaps( &
                & equationsSetIdx)%equationsSetIndex,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Interface matrix number = ",solverMapping% &
                & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%equationsSetToInterfaceConditionMaps( &
                & equationsSetIdx)%interfaceMatrixIndex,err,error,*999)
            ENDDO !equationsSetIdx
          ENDDO !solverMatrixIdx
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Interface condition rows to solver rows mappings:",err,error,*999)
          DO interfaceMatrixIdx=1,interfaceMapping%numberOfInterfaceMatrices
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Interface matrix idx : ",interfaceMatrixIdx,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of interface rows = ",interfaceMapping% &
              & interfaceMatrixRowsToVarMaps(interfaceMatrixIdx)%totalNumberOfRows,err,error,*999)
            DO rowIdx=1,interfaceMapping%interfaceMatrixRowsToVarMaps(interfaceMatrixIdx)%totalNumberOfRows
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Interface row : ",rowIdx,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Number of solver rows mapped to = ",solverMapping% &
                & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%interfaceMatrixToSolverMatricesMaps( &
                & interfaceMatrixIdx)%interfaceRowToSolverRowsMap(rowIdx)%numberOfRowCols,err,error,*999)
              IF(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%interfaceMatrixToSolverMatricesMaps( &
                & interfaceMatrixIdx)%interfaceRowToSolverRowsMap(rowIdx)%numberOfRowCols>0) THEN                
                CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%interfaceConditionToSolverMatricesMaps( &
                  & interfaceConditionIdx)%interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)% &
                  & interfaceRowToSolverRowsMap(rowIdx)%numberOfRowCols,5,5,solverMapping% &
                  & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%interfaceMatrixToSolverMatricesMaps( &
                  & interfaceMatrixIdx)%interfaceRowToSolverRowsMap(rowIdx)%rowCols, &
                  & '("          Solver row numbers    :",5(X,I13))','(33X,5(X,I13))',err,error,*999)
                CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%interfaceConditionToSolverMatricesMaps( &
                  & interfaceConditionIdx)%interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)% &
                  & interfaceRowToSolverRowsMap(rowIdx)%numberOfRowCols,5,5,solverMapping% &
                  & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%interfaceMatrixToSolverMatricesMaps( &
                  & interfaceMatrixIdx)%interfaceRowToSolverRowsMap(rowIdx)%couplingCoefficients, &
                  & '("          Coupling coefficients :",5(X,E13.6))','(33X,5(X,E13.6))',err,error,*999)
              ENDIF
            ENDDO !rowIdx
          ENDDO !interfaceMatrixIdx
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Interface condition column to solver rows mappings:", &
            & err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of interface condition columns = ",interfaceMapping% &
            & totalNumberOfColumns,err,error,*999)
          DO columnIdx=1,interfaceMapping%totalNumberOfColumns
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Interface condition column : ",columnIdx,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of rows mapped to = ",solverMapping% &
              & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%interfaceColToSolverRowsMap(columnIdx)% &
              & numberOfRowCols,err,error,*999)
            IF(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)% &
              & interfaceColToSolverRowsMap(columnIdx)%numberOfRowCols>0) THEN
              CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping% &
                & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%interfaceColToSolverRowsMap( &
                & rowIdx)%numberOfRowCols,5,5,solverMapping% &
                & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%interfaceColToSolverRowsMap(rowIdx)% &
                & rowCols,'("              Solver rows           :",5(X,I13))','(38X,5(X,I13))',err,error,*999)
              CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping% &
                & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%interfaceColToSolverRowsMap( &
                & rowIdx)%numberOfRowCols,5,5,solverMapping% &
                & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%interfaceColToSolverRowsMap(rowIdx)% &
                & couplingCoefficients,'("              Coupling coefficients  :",5(X,E13.6))','(38X,5(X,E13.6))', &
                & err,error,*999)
            ENDIF
          ENDDO !columnIdx
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Solver matrix indexing:",err,error,*999)
          DO solverMatrixIdx=1,solverMapping%numberOfSolverMatrices
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Solver matrix : ",solverMatrixIdx,err,error,*999)        
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"        Interface equations matrix rows to solver matrix columns:", &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Number of interface matrices = ",solverMapping% &
              & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%jacobianMatrixToSolverMatrixMap%numberOfInterfaceMatrices, &
              & err,error,*999)
            DO interfaceMatrixIdx=1,solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%numberOfInterfaceMatrices
              interfaceMatrixToSolverMatrixMap=>solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)% &
                & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%interfaceMatrixToSolverMatrixMaps( &
                & interfaceMatrixIdx)%ptr
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Interface matrix index : ",interfaceMatrixIdx, &
                & err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Interface matrix number = ", &
                & interfaceMatrixToSolverMatrixMap%interfaceMatrixNumber,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Solver matrix number    = ", &
                & interfaceMatrixToSolverMatrixMap%solverMatrixNumber,err,error,*999)
              DO rowIdx=1,interfaceMapping%interfaceMatrixRowsToVarMaps(interfaceMatrixIdx)%numberOfRows
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Interface matrix row : ",rowIdx, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Number of solver columns mapped to = ", &
                  & interfaceMatrixToSolverMatrixMap%interfaceRowToSolverColsMap(rowIdx)%numberOfRowCols,err,error,*999)
                IF(interfaceMatrixToSolverMatrixMap%interfaceRowToSolverColsMap(rowIdx)%numberOfRowCols>0) THEN
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,interfaceMatrixToSolverMatrixMap%interfaceRowToSolverColsMap( &
                    & rowIdx)%numberOfRowCols,5,5,interfaceMatrixToSolverMatrixMap%interfaceRowToSolverColsMap(rowIdx)% &
                    & rowCols,'("              Solver columns         :",5(X,I13))','(38X,5(X,I13))',err,error,*999)
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,interfaceMatrixToSolverMatrixMap%interfaceRowToSolverColsMap( &
                    & rowIdx)%numberOfRowCols,5,5,interfaceMatrixToSolverMatrixMap%interfaceRowToSolverColsMap(rowIdx)% &
                    & couplingCoefficients,'("              Coupling coefficients  :",5(X,E13.6))','(38X,5(X,E13.6))', &
                    & err,error,*999)
                ENDIF
              ENDDO !rowIdx
            ENDDO !interfaceMatrixIdx
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"        Variable dofs to solver matrix dofs:",err,error,*999)
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"          Lagrange variables:",err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Lagrange variable type = ",solverMapping% &
              & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%lagrangeVariableType,err,error,*999)
            lagrangeVariable=>solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%lagrangeVariable
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of Lagrange variable dofs = ",lagrangeVariable% &
              & numberOfDofs,err,error,*999)
            DO localDOF=1,lagrangeVariable%totalNumberOfDofs
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Variable dof : ",localDOF,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Solver column number = ",solverMapping% &
                & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%interfaceMatricesToSolverMatrixMaps( &
                & solverMatrixIdx)%lagrangeVariableToSolverColMap%columnNumbers(localDOF),err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Coupling coefficient = ",solverMapping% &
                & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%interfaceMatricesToSolverMatrixMaps( &
                & solverMatrixIdx)%lagrangeVariableToSolverColMap%couplingCoefficients(localDOF),err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Additive constant    = ",solverMapping% &
                & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%interfaceMatricesToSolverMatrixMaps( &
                & solverMatrixIdx)%lagrangeVariableToSolverColMap%additiveConstants(localDOF),err,error,*999)              
            ENDDO !localDOF
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"          Dependent variables:",err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Number of dependent variables = ",solverMapping% &
              & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%numberOfDependentVariables,err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%interfaceConditionToSolverMatricesMaps( &
              & interfaceConditionIdx)%interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%numberOfDependentVariables, &
              & 5,5,solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%dependentVariableTypes, &
              & '("            Dependent variable types :",5(X,I13))','(38X,5(X,I13))',err,error,*999) 
            DO variableIdx=1,solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)% &
              & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%numberOfDependentVariables
              dependentVariable=>solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)% &
                & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%dependentVariables(variableIdx)%ptr
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Dependent variable index : ",variableIdx,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Number of dependent variable dofs = ", &
                & dependentVariable%numberOfDofs,err,error,*999)
              DO localDOF=1,dependentVariable%totalNumberOfDofs
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"              Variable dof : ",localDOF,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Solver column number = ",solverMapping% &
                  & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%interfaceMatricesToSolverMatrixMaps( &
                  & solverMatrixIdx)%dependentVariableToSolverColMaps(variableIdx)%columnNumbers(localDOF), &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Coupling coefficient = ",solverMapping% &
                  & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%interfaceMatricesToSolverMatrixMaps( &
                  & solverMatrixIdx)%dependentVariableToSolverColMaps(variableIdx)%couplingCoefficients(localDOF), &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"                Additive constant    = ",solverMapping% &
                  & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%interfaceMatricesToSolverMatrixMaps( &
                  & solverMatrixIdx)%dependentVariableToSolverColMaps(variableIdx)%additiveConstants(localDOF), &
                  & err,error,*999)              
              ENDDO !localDOF
            ENDDO !variableIdx
          ENDDO !solverMatrixIdx        
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Interface equations matrix indexing:",err,error,*999)
          DO interfaceMatrixIdx=1,interfaceMapping%numberOfInterfaceMatrices
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Interface matrix : ",interfaceMatrixIdx,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of solver matrices = ",solverMapping% &
              & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%interfaceMatrixToSolverMatricesMaps( &
              & interfaceMatrixIdx)%numberOfSolverMatrices,err,error,*999)
            DO solverMatrixIdx=1,solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)% &
              & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%numberOfSolverMatrices
              interfaceMatrixToSolverMatrixMap=>solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)% &
                & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%interfaceMatrixToSolverMatrixMaps( &
                & solverMatrixIdx)%ptr
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Solver matrix index : ",solverMatrixIdx,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Interface matrix number = ", &
                & interfaceMatrixToSolverMatrixMap%interfaceMatrixNumber,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Solver matrix number    = ", &
                & interfaceMatrixToSolverMatrixMap%solverMatrixNumber,err,error,*999)            
            ENDDO !solverMatrixIdx
          ENDDO !equationMatrixIdx
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Interface column to solver rows mapping:",err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of columns = ",interfaceMapping%numberOfColumns, &
            & err,error,*999)            
          DO columnIdx=1,interfaceMapping%numberOfColumns
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Column : ",columnIdx, err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of solver rows = ",solverMapping% &
              & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%interfaceColToSolverRowsMap(columnIdx)% &
              & numberOfRowCols,err,error,*999)
            IF(solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%interfaceColToSolverRowsMap( &
              & columnIdx)%numberOfRowCols>0) THEN
              CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%interfaceConditionToSolverMatricesMaps( &
                & interfaceConditionIdx)%interfaceColToSolverRowsMap(rowIdx)%numberOfRowCols,5,5,solverMapping% &
                & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%interfaceColToSolverRowsMap(rowIdx)% &
                & rowCols,'("        Solver rows           :",5(X,I13))','(32X,5(X,I13))',err,error,*999)
              CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,solverMapping%interfaceConditionToSolverMatricesMaps( &
                & interfaceConditionIdx)%interfaceColToSolverRowsMap(rowIdx)%numberOfRowCols,5,5,solverMapping% &
                & interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%interfaceColToSolverRowsMap(rowIdx)% &
                & couplingCoefficients,'("        Coupling coefficients  :",5(X,E13.6))','(32X,5(X,E13.6))', &
                & err,error,*999)
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
    CALL SolverEquation_AssertNotFinished(solverEquations,err,error,*999)
    
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
      IF(ASSOCIATED(createValuesCache%equationsRHSVariablesList)) CALL List_Destroy(createValuesCache% &
        & equationsRHSVariablesList,err,error,*999)
      IF(ALLOCATED(createValuesCache%dynamicVariableType)) DEALLOCATE(createValuesCache%dynamicVariableType)
      IF(ALLOCATED(createValuesCache%matrixVariableTypes)) DEALLOCATE(createValuesCache%matrixVariableTypes)
      IF(ALLOCATED(createValuesCache%residualVariableTypes)) DEALLOCATE(createValuesCache%residualVariableTypes)
      IF(ALLOCATED(createValuesCache%rhsVariableType)) DEALLOCATE(createValuesCache%rhsVariableType)
      IF(ALLOCATED(createValuesCache%sourceVariableType)) DEALLOCATE(createValuesCache%sourceVariableType)
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
    ALLOCATE(solverMapping%createValuesCache%residualVariableTypes(0:FIELD_NUMBER_OF_VARIABLE_TYPES, &
      & solverMapping%numberOfEquationsSets),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver mapping create values cache residual variable type.",err,error,*999)
    ALLOCATE(solverMapping%createValuesCache%rhsVariableType(solverMapping%numberOfEquationsSets),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver mapping create values cache RHS variable type.",err,error,*999)
    ALLOCATE(solverMapping%createValuesCache%sourceVariableTypes(0:FIELD_NUMBER_OF_VARIABLE_TYPES, &
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
    NULLIFY(solverMapping%createValuesCache%equationsRHSVariablesList)
    CALL List_CreateStart(solverMapping%createValuesCache%equationsRHSVariablesList,err,error,*999)
    CALL List_DataTypeSet(solverMapping%createValuesCache%equationsRHSVariablesList,LIST_INTG_TYPE,err,error,*999)
    CALL List_DataDimensionSet(solverMapping%createValuesCache%equationsRHSVariablesList,2,err,error,*999)
    CALL List_CreateFinish(solverMapping%createValuesCache%equationsRHSVariablesList,err,error,*999)
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
    TYPE(EquationsSetType), POINTER :: equationsSet,variablEquationsSet
    TYPE(FieldType), POINTER :: dependentField,variableDependentField
    TYPE(ListType),  POINTER :: equationsVariableList
    TYPE(RegionType), POINTER :: dependentRegion,variableRegion
    TYPE(SolverMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverMapping_CreateValuesCacheEqnVarListAdd",err,error,*999)

    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)
    NULLIFY(createValuesCache)
    CALL SolverMapping_CreateValuesCacheGet(solverMapping,createValuesCache,err,error,*999)
    NULLIFY(equationsVariableList)
    CALL SolverMappingCreateValuesCache_EquationsVariableListGet(createValuesCache,solverMatrixIdx,equationsVariableList, &
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
          CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx2,variableEquationsList,err,error,*999)
          NULLIFY(variableDependentField)
          CALL EquationsSet_DependentFieldGet(variableEquationsSet,variableDependentField,err,error,*999)
          IF(ASSOCIATED(dependentField,variableDependentField)) THEN
            NULLIFY(variableRegion)
            CALL Field_RegionGet(variableDependentField,variableRegion,err,error,*999)
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
    LOGICAL :: varieableFound
    TYPE(EquationsSetType), POINTER :: equationsSet,variableEquationsSet
    TYPE(FieldType), POINTER :: dependentField,variableDependentField
    TYPE(FieldVariableType), POINTER :: rhsVariable,variableRhsVariable
    TYPE(ListType), POINTER :: rhsVariableList
    TYPE(SolverMappingCreateValuesCacheType), POINTER :: createValuesCache

    ENTERS("SolverMapping_CreateValuesCacheEqnRHSVarListAdd",err,error,*999)

    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)

    NULLIFY(createValuesCache)
    CALL SolverMapping_CreateValuesCacheGet(solverMapping,createValuesCache,err,error,*999)
    NULLIFY(equationsVariableList)
    CALL SolverMappingCreateValuesCache_RHSVariableListGet(createValuesCache,rhsVariableList,err,error,*999)
    NULLIFY(equationsSet)
    CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(rhsVariable)
    CALL Field_VariableGet(dependentField,rhsVariableType,rhsVariable,err,error,*999)
    NULLIFY(dependentRegion)
    CALL Field_RegionGet(dependentField,dependentRegion)
    variableFound=.FALSE.
    CALL List_NumberOfItemsGet(rhsVariablesList,numberOfVariables,err,error,*999)
    DO variableIdx=1,numberOfVariables
      CALL List_ItemGet(rhsVariablesList,variableIdx,variableItem,err,error,*999)
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
      CALL List_ItemAdd(rhsVariablesList,variableItem,err,error,*999)
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
    CALL SolverMappingCreateValuesCache_InterfaceVariableListGet(createValuesCache,solverMatrixIdx,interfaceVariableList, &
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
                IF(ASSOCIATED(lagrangeInterface,variableLagrangeField)) THEN
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
    INTEGER(INTG) :: linearVariableIndex,variableIdx
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
    INTEGER(INTG) :: equationsSetIdx,matrixIdx,solverMatrixIdx,variableIdx,variableType
    INTEGER(INTG), ALLOCATABLE :: newDynamicVariableTypes(:),newMatrixVariableTypes(:,:,:),newRHSVariableTypes(:), &
      & newResidualVariableTypes(:,:),newSourceVariableTypes(:,:)
    LOGICAL :: matrixDone
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
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
    ALLOCATE(newDynamicVariableTypes(solverMapping%numberOfEquationsSets+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new dynamic variable types.",err,error,*999)
    ALLOCATE(newMatrixVariableTypes(0:FIELD_NUMBER_OF_VARIABLE_TYPES,solverMapping%numberOfEquationsSets+1, &
      & solverMapping%numberOfSolverMatrices),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new matrix variable types.",err,error,*999)
    ALLOCATE(newResidualVariableTypes(0:FIELD_NUMBER_OF_VARIABLE_TYPES,solverMapping%numberOfEquationsSets+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new residual variable types.",err,error,*999)
    ALLOCATE(newRHSVariableTypes(solverMapping%numberOfEquationsSets+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new RHS variable types.",err,error,*999)
    ALLOCATE(newSourceVariableTypes(0:FIELD_NUMBER_OF_VARIABLE_TYPES,solverMapping%numberOfEquationsSets+1),STAT=err)
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
      newResidualVariableTypes(:,1:solverMapping%numberOfEquationsSets)=createValuesCache%residualVariableTypes
      newRHSVariableTypes(1:solverMapping%numberOfEquationsSets)=createValuesCache%rhsVariableType
      newSourceVariableTypes(:,1:solverMapping%numberOfEquationsSets)=createValuesCache%sourceVariableTypes
      DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
        newEquationsSets(equationsSetIdx)%ptr=>solverMapping%equationsSets(equationsSetIdx)%ptr
      ENDDO !equationsSetIdx              
      newDynamicVariableTypes(solverMapping%numberOfEquationsSets+1)=0
      newMatrixVariableTypes(:,solverMapping%numberOfEquationsSets+1,:)=0
      newResidualVariableTypes(:,solverMapping%numberOfEquationsSets+1)=0
      newRHSVariableTypes(solverMapping%numberOfEquationsSets+1)=0
      newSourceVariableTypes(:,solverMapping%numberOfEquationsSets+1)=0            
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
        variableType=1
        DO matrixIdx=1,solverMapping%numberOfSolverMatrices
          matrixDone=.FALSE.
          DO WHILE(variableType<=FIELD_NUMBER_OF_VARIABLE_TYPES.AND..NOT.matrixDone)
            CALL EquationsMappingLinear_LinearVariableIndexGet(linearMapping,variableType,variableIndex,err,error,*999)
            IF(variableIndex==0) THEN
              variableType=variableType+1
            ELSE
              IF(linearMapping%varToEquationsMatricesMaps(variableIndex)%numberOfEquationsMatrices>0) THEN
                createValuesCache%matrixVariableTypes(0,solverMapping%numberOfEquationsSets+1,matrixIdx)=1
                createValuesCache%matrixVariableTypes(variableType,solverMapping%numberOfEquationsSets+1,matrixIdx)=variableType
                matrixDone=.TRUE.
              ELSE
                variableType=variableType+1
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
          createValuesCache%residualVariableTypes(0,solverMapping%numberOfEquationsSets+1)=numberOfResidualVariables
          !Map the residual variables to the solver Jacobian
          DO variableIdx=1,numberOfResidualVariables
            CALL EquationsMappingResidual_VariableTypeGet(residualMapping,variableIdx,residualVariableType,err,error,*999)
            createValuesCache%residualVariableTypes(residualVariableType,solverMapping%numberOfEquationsSets+1)= &
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
              CALL EquationsMappingLinear_VariableIndexGet(residualMapping,residuaVariableType,variableIndex,err,error,*999)
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
        CALL EquationsMappingLinear_NumberOfVariablesGet(linearMapping,numberOfLinearVariables,err,error,*999)
        DO variableIdx=1,numberOfLinearVariables
          CALL EquationsMappingLinear_VariableIndexGet(residualMapping,dynamicVariableType,variableIndex,err,error,*999)
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
      CALL EquationsMappingRHS_VariableTypeGet(rhsMapping,rhsVariableType,err,error,*999)
      createValuesCache%rhsVariableType(solverMapping%numberOfEquationsSets+1)=rhsVariableType
    ENDIF
    NULLIFY(sourcesMapping)
    CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*999)
    IF(ASSOCIATED(sourcesMapping)) THEN
      CALL EquationsMappingSources_NumberOfSourcesGet(sourceMapping,numberOfSources,err,error,*999)
      DO sourceIdx=1,numberOfSources
        NULLIFY(sourceMapping)
        CALL EquationsMappingSources_SourceMappingGet(sourcesMapping,sourceIdx,sourceMapping,err,error,*999)
        CALL EquationsMappingSource_NumberOfVariablesGet(sourceMapping,numberOfVariables,err,error,*999)
        createValuesCache%sourceVariableTypes(0,solverMapping%numberOfEquationsSets+1)=numberOfVariables
        DO variableIdx=1,numberOfVariables
          CALL EquationsMappingSource_VariableTypeGet(sourceMapping,variableIdx,sourceVariableType,err,error,*999)
          createValuesCache%sourceVariableTypes(sourceVariableType,solverMapping%numberOfEquationsSets+1)=sourceVariableType
        ENDDO !variableIdx
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
        variableType=createValuesCache%matrixVariableTypes(variableTypeIdx,equationsSetIndex,solverMatrixIdx)
        CALL SolverMapping_CreateValuesCacheEqnVarListAdd(solverMapping,1,equationsSetIndex,variableType,err,error,*999)
      ENDDO !matrixIdx
    ENDDO !solverMatrixIdx
    DO variableTypeIdx=1,FIELD_NUMBER_OF_VARIABLE_TYPES
      variableType=createValuesCache%residualVariableTypes(variableTypeIdx,equationsSetIndex)
      CALL SolverMapping_CreateValuesCacheEqnVarListAdd(solverMapping,1,equationsSetIndex,variableType,err,error,*999)
    ENDDO
    rhsVariableType=createValuesCache%rhsVariableType(equationsSetIndex)
    CALL SolverMapping_CreateValuesCacheEqnRHSVarListAdd(solverMapping,equationsSetIndex,rhsVariableType,err,error,*999)
        
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
          CALL SolverMappingSMICToES_inalise(equationsSetToSolverMatricesMap% &
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

    ENTERS("SolverMapping_EquationsSetToSolverMapInitialise",err,error,*998)

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

    EXITS("SolverMapping_EquationsSetToSolverMapInitialise")
    RETURN
999 CALL SolverMapping_EquationsSetToSolverMapFinalise(equationsSetToSolverMatricesMap,dummyErr,dummyError,*998)
998 ERRORS("SolverMapping_EquationsSetToSolverMapInitialise",err,error)    
    EXITS("SolverMapping_EquationsSetToSolverMapInitialise")    
    RETURN 1
   
  END SUBROUTINE SolverMapping_EquationsSetToSolverMapInitialise

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
      IF(ALLOCATED(equationsMatrixToSolverMatrixMap%equationsColToSolverColsMap)) THEN
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
    equationsMatrixToSolverMatrixMap%equationMatrixType=0
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
        DO equationsMatrixIdx=1,SIZE(equationsMatrixToSolverMatricesMaps%equationsMatrixToSolverMatrixMaps,1)
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

    ALLOCATE(equationsMatrixToSolverMatricesMaps,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocated equations matrix to solver matrices map.",err,error,*999)
    equationsMatrixToSolverMatricesMaps%equationsMatrixNumber=0
    equationsMatrixToSolverMatricesMaps%numberOfSolverMatrices=0
        
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

  !>Finalises a Jacobian matrix to solver matrix map and deallocates all memory.
  SUBROUTINE SolverMappingJMToSMMap_Finalise(jacobianMatrixToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixToSolverMatrixMapType), POINTER :: jacobianMatrixToSolverMatrixMap !<The Jacobian matrix to solver matrix map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx
 
    ENTERS("SolverMappingJMToSMMap_Finalise",err,error,*999)

    IF(ASSOCIATED(jacobianMatrixToSolverMatrixMap)) THEN
      IF(ASSOCIATED(jacobianMatrixToSolverMatrixMap%jacobianColToSolverColsMap)) THEN
        DO columnIdx=1,SIZE(jacobianMatrixToSolverMatrixMaps,1)
          CALL MatrixRowColCoupling_Finalise(jacobianMatrixToSolverMatrixMap%jacobianColToSolverColsMap(columnIdx), &
            & err,error,*999)
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
      IF(ALLOCATED(equationsMatricesToSolverMatrixMap%variables)) DEALLOCATE(equationsMatricesToSolverMatrixMaps%variables)
      IF(ALLOCATED(equationsMatricesToSolverMatrixMap%variableToSolverColMaps)) THEN
        DO variableIdx=1,SIZE(equationsMatricesToSolverMatrixMap%variableToSolverColMaps,1)
          CALL SolverMappingVToSCMap_Finalise(equationsMatricesToSolverMatrixMap%variableToSolverColMaps(variableIdx)%ptr, &
            & err,error,*999)
        ENDDO !variableIdx
        DEALLOCATE(equationsMatricesToSolverMatrixMaps%variableToSolverColMaps)
      ENDIF
      IF(ALLOCATED(equationsMatricesToSolverMatrixMap%dynamicEquationsMatrixToSolverMatrixMaps)) THEN
        DO matrixIdx=1,SIZE(equationsMatricesToSolverMatrixMaps%dynamicEquationsMatrixToSolverMatrixMaps,1)
          CALL SolverMappingEMToSMMap_Finalise(equationsMatricesToSolverMatrixMap% &
            & dynamicEquationsMatrixToSolverMatrixMaps(matrixIdx)%ptr,err,error,*999)        
        ENDDO !matrixIdx
        DEALLOCATE(equationsMatricesToSolverMatrixMap%dynamicEquationsMatrixToSolverMatrixMaps)
      ENDIF
      IF(ALLOCATED(equationsMatricesToSolverMatrixMap%linearEquationsMatrixToSolverMatrixMaps)) THEN
        DO matrixIdx=1,SIZE(equationsMatricesToSolverMatrixMaps%linearEquationsMatrixToSolverMatrixMaps,1)
          CALL SolverMappingEMToSMMap_Finalise(equationsMatricesToSolverMatrixMap% &
            & linearEquationsMatrixToSolverMatrixMaps(matrixIdx)%ptr,err,error,*999)        
        ENDDO !matrixIdx
        DEALLOCATE(equationsMatricesToSolverMatrixMap%linearEquationsMatrixToSolverMatrixMaps)
      ENDIF
      IF(ALLOCATED(equationsMatricesToSolverMatrixMap%jacobianMatrixToSolverMatrixMaps)) THEN
        DO matrixIdx=1,SIZE(equationsMatricesToSolverMatrixMap%jacobianMatrixToSolverMatrixMaps,1)
          CALL SolverMappingJMToSMMap_Finalise(equationsMatricesToSolverMatrixMaps% &
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
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMaps !<A pointer to the equations matrices to solver matrix map to initialise. Must not be associated on entry
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
    equationsMatricesToSolverMatrixMap%numberOfDynamicEquationsMatrices=0
    equationsMatricesToSolverMatrixMap%numberOfLinearEquationsMatrices=0
    equationsMatricesToSolverMatrixMap%numberOfEquationsJacobians=0

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
      IF(ALLOCATED(solverMapping%solverMatrixToEquationsMaps)) THEN
        DO solverMatrixIdx=1,SIZE(solverMapping%solverMatrixToEquationsMaps,1)
          CALL SolverMappingSMToEquationsMaps_Finalise(solverMapping%solverMatrixToEquationsMaps(solverMatrixIdx)%ptr, &
            & err,error,*999)
        ENDDO !solverMatrixIdx
        DEALLOCATE(solverMapping%solverMatrixToEquationsMaps)
      ENDIF
      CALL SolverMappingVariables_Finalise(solverMapping%rhsVariablesList,err,error,*999)      
      IF(ALLOCATED(solverMapping%solverColToEquationsColsMap)) THEN
        DO solverMatrixIdx=1,SIZE(solverMapping%solverColToEquationsColsMap,1)
          CALL SolverMappingSMToEquationsMap_Finalise(solverMapping%solverColToEquationsColsMap( &
            & solverMatrixIdx),err,error,*999)
        ENDDO !solverMatrixIdx
        DEALLOCATE(solverMapping%solverColToEquationsColsMap)
      ENDIF
      IF(ALLOCATED(solverMapping%solverRowToEquationsRowsMap)) THEN
        DO rowIdx=1,SIZE(solverMapping%solverRowToEquationsRowsMap,1)
          CALL SolverMappingSRToERSMap_Finalise(solverMapping%solverRowToEquationsRowsMap( &
            & rowIdx),err,error,*999)
        ENDDO !rowIdx
        DEALLOCATE(solverMapping%solverRowToEquationsRowsMap)
      ENDIF
      CALL DomainMapping_Finalise(solverMapping%rowDOFsMapping,err,error,*999)
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
    solverEquations%solverMapping%numberOfInterfaceConditions=07
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
    INTEGER(INTG) :: equationMatrixIdx,equationsSetIdx,interfaceConditionIdx,interfaceMatrixIdx,listItem(2)
    INTEGER(INTG) :: numberOfInterfaceMatrices
    LOGICAL :: equationsSetFound,variableFound
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldVariableType), POINTER :: dependentVariable
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(InterfaceConditionPtrType), ALLOCATABLE :: newInterfaceConditions(:)
    TYPE(SolverEquationsType), POINTER :: solverEquations
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
        NULLIFY(interfaceMatrixRowsToVarMap)
        CALL InterfaceMapping_InterfaceMatrixRowsToVarMapsGet(interfaceMapping,interfaceMatrixIdx,interfaceMatrixRowsToVarMap, &
          & err,error,*999)
        CALL InterfaceMappingIMRSToVMap_EquationsSetGet(interfaceMatrixRowsToVarMap,equationsSet,err,error,*999)
        equationsSetFound=.FALSE.
        DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
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
        CALL InterfaceMappingIMRSToVMap_VariableGet(interfaceMatrixRowsToVarMap,dependentVariable,err,error,*999)
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
            DO equationMatrixIdx=1,createValuesCache%residualVariableTypes(0,equationsSetIdx)
              IF(createValuesCache%residualVariableTypes(equationMatrixIdx,equationsSetIdx)==dependentVariableType) THEN
                variableFound=.TRUE.
              ENDIF
            ENDDO !equationsMatrixIdx
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
    IF(solverMapping%numberOfInterfaceConditions>0) THEN
      ALLOCATE(newInterfaceConditions(solverMapping%numberOfInterfaceConditions+1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new interface conditions.",err,error,*999)
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
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(SolverMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverMapping_NumberOfSolverMatricesSet",err,error,*999)

    CALL SolverMapping_AssertNotFinished(solverMapping,err,error,*999)
    NULLIFY(createValuesCache)
    CALL SolverMapping_CreateValuesCacheGet(solverMapping,createValuesCache,err,error,*999)
    maximumNumberOfEquationsMatrices=1
    DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationSet,err,error,*999)
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
        DO equationsSetIdx=1,SIZE(interfaceConditionToSolverMatricesMaps%equationsSetToInterfaceConditionMaps,1)
          CALL SolverMappingsSMESToICMap_Finalise(interfaceConditionToSolverMatricesMap% &
            & equationsSetToInterfaceConditionMaps(equationsSetIdx),err,error,*999)
        ENDDO !equationsSetIdx
        DEALLOCATE(interfaceConditionToSolverMatricesMap%equationsSetToInterfaceConditionMaps)
      ENDIF
      IF(ALLOCATED(interfaceConditionToSolverMatricesMap%interfaceMatricesToSolverMatrixMaps)) THEN
        DO solverMatrixIdx=1,SIZE(interfaceConditionToSolverMatricesMaps%interfaceMatricesToSolverMatrixMaps,1)
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
      IF(ASSOCIATED(interfaceConditionToSolverMatricesMaps%interfaceColToSolverRowsMap)) THEN
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
  SUBROUTINE SolverMappingICToSMSMap_Initialise(interfaceConditionToSolverMatricesMas,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMaps !<A pointer to the interface condition to solver matrices map to initialise. Must not be associated on entry
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
    interfaceConditionToSolverMatricesMaps%interfaceConditionIndex=0
    NULLIFY(interfaceConditionToSolverMatricesMaps%solverMapping)
    NULLIFY(interfaceConditionToSolverMatricesMaps%interfaceEquations)
    interfaceConditionToSolverMatricesMaps%numberOfEquationsSets=0
    interfaceConditionToSolverMatricesMaps%numberOfSolverMatricess=0
    interfaceConditionToSolverMatricesMaps%numberOfInterfaceMatricess=0
    NULLIFY(interfaceConditionToSolverMatricesMaps%interfaceColToSolverRowsMap)
    
    EXITS("SolverMappingICToSMSMap_Initialise")
    RETURN
999 CALL SolverMapppingICToSMSMap_Finalise(interfaceCOnditionToSolverMatricesMap,dummyErr,dummyError,*998)
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
      IF(ALLOCATED(interfaceMatrixToSolverMatrixMap%interfaceRowToSolverColsMap)) THEN
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
    NULLIFY(interfaceMatrixToSolverMatrixMapinterfaceRowToSolverColsMap)
        
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
    equationsSetToInterfaceConditionMapy%interfaceMatrixIndex=0
         
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
        DO solverMatrixIdx=1,SIZE(interfaceMatrixToSolverMatricesMaps%interfaceMatrixToSolverMatrixMaps,1)
          CALL SolverMappingIMToSMMap_Finalise(interfaceMatrixToSolverMatricesMap% &
            & interfaceMatrixToSolverMatrixMaps(solverMatrixIdx)%ptr,err,error,*999)        
        ENDDO !solverMatrixIdx
        DEALLOCATE(interfaceMatrixToSolverMatricesMap%interfaceMatrixToSolverMatrixMaps)
      ENDIF
      IF(ALLOCATED(interfaceMatrixToSolverMatricesMaps%interfaceRowToSolverRowsMap)) THEN
        DO rowIdx=1,SIZE(interfaceMatrixToSolverMatricesMap%interfaceRowToSolverRowsMap,1)
          CALL MatrixRowColCoupling_Finalise(interfaceMatrixToSolverMatricesMap%interfaceRowToSolverRowsMap(rowIdx),err,error,*999)
        ENDDO !rowIdx
        DEALLOCATE(interfaceMatrixToSolverMatricesMaps%interfaceRowToSolverRowsMap)
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
999 ERRORS("SolverMappingIMToSMSMap_Initialise",err,error)    
    EXITS("SolverMappingIMToSMSMap_Initialise")    
    RETURN 1
   
  END SUBROUTINE SolverMappingIMToSMSMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises an interface matrix to solver matrices map and deallocates all memory.
  SUBROUTINE SolverMappingIMSToSMMap_Finalise(interfaceMatricesToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesToSolverMatrixMapsType), POINTER :: interfaceMatricesToSolverMatrixMap !<A pointer to the interface matrices to solver matrix map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,interfaceMatrixIdx
    
    ENTERS("SolverMappingIMSToSMMap_Finalise",err,error,*999)

    IF(ASSOCIATED(interfaceMatricesToSolverMatrixMap)) THEN
       CALL SolverMappingVToSCMap_Finalise(interfaceMatricesToSolverMatrixMap%lagrangeVariableToSolverColMap,err,error,*999)
       IF(ALLOCATED(interfaceMatricesToSolverMatrixMap%dependentVariableTypes)) &
         & DEALLOCATE(interfaceMatricesToSolverMatrixMap%dependentVariableTypes)
       IF(ALLOCATED(interfaceMatricesToSolverMatrixMap%dependentVariables)) &
         & DEALLOCATE(interfaceMatricesToSolverMatrixMap%dependentVariables)
       IF(ALLOCATED(interfaceMatricesToSolverMatrixMap%dependentVariableToSolverColMaps)) THEN
        DO interfaceMatrixIdx=1,SIZE(interfaceMatricesToSolverMatrixMaps%dependentVariableToSolverColMaps,1)
          CALL SolverMappingVToSCMap_Finalise(interfaceMatricesToSolverMatrixMap% &
            & dependentVariableToSolverColMaps(interfaceMatrixIdx)%ptr,err,error,*999)
        ENDDO !interfaceMatrixIdx
        DEALLOCATE(interfaceMatricesToSolverMatrixMap%dependentVariableToSolverColMaps)
      ENDIF
      IF(ALLOCATED(interfaceMatricesToSolverMatrixMaps%interfaceMatrixToSolverMatrixMaps)) THEN
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
    TYPE(InterfaceMatricesToSolverMatrixMapsType), POINTER :: interfaceMatricesToSolverMatrixMap !<A pointer to the interface matrices to solver matrix map to initialise. Must not be associated on entry.
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
    NULLIFY(interfaceMatricesToSolverMatrixMap%lagrangeVariableToSolverColMap)
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
      IF(ALLOCATED(jacobianMatrixToSolverMatrixMap%jacobianColToSolverColsMap)) THEN
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
  SUBROUTINE SolverMappingSCToDEquationsMap_Finalise(solverColToDynamicEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverColToDynamicEquationsMapType) :: solverColToDynamicEquationsMap !<The solver column to dynamic equations map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverMappingSCToDEquationsMap_Finalise",err,error,*999)

    solverColToDynamicEquationsMap%numberOfDynamicEquationsMatrices=0
    IF(ALLOCATED(solverColToDynamicEquationsMap%equationsMatrixNumbers)) &
      & DEALLOCATE(solverColToDynamicEquationsMap%equationsMatrixNumbers)
    IF(ALLOCATED(solverColToDynamicEquationsMap%equationsColumnNumbers)) &
      & DEALLOCATE(solverColToDynamicEquationsMap%equationsColumnNumbers)
    IF(ALLOCATED(solverColToDynamicEquationsMap%couplingCoefficients)) &
      & DEALLOCATE(solverColToDynamicEquationsMap%couplingCoefficients)
       
    EXITS("SolverMappingSCToDEquationsMap_Finalise")
    RETURN
999 ERRORS("SolverMappingSCToDEquationsMap_Finalise",err,error)
    EXITS("SolverMappingSCToDEquationsMap_Finalise")
    RETURN 1
    
  END SUBROUTINE SolverMappingSCToDEquationsMap_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the solver column to dynamic equations mapping and deallocates all memory.
  SUBROUTINE SolverMappingsSCToDEquationsMap_Initialise(solverColToDynamicEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverColToDynamicEquationsMapType) :: solverColToDynamicEquationsMap !<The solver column to dynamic equations map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverMappingsSCToDEquationsMap_Initialise",err,error,*999)

    solverColToDynamicEquationsMap%numberOfDynamicEquationsMatrices=0
     
    EXITS("SolverMappingsSCToDEquationsMap_Initialise")
    RETURN
999 ERRORS("SolverMappingsSCToDEquationsMap_Initialise",err,error)
    EXITS("SolverMappingsSCToDEquationsMap_Initialise")
    RETURN 1
    
  END SUBROUTINE SolverMappingsSCToDEquationsMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver column to static equations map and deallocates all memory.
  SUBROUTINE SolverMappingSCToSEquationsMap_Finalise(solverColToSEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverColToStaticEquationsMapType) :: solverColToSEquationsMap !<The solver column to static equations map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverMappingSCToSEquationsMap_Finalise",err,error,*999)

    solverColToSEquationsMap%numberOfLinearEquationsMatrices=0
    IF(ALLOCATED(solverColToSEquationsMap%equationsMatrixNumbers)) &
      & DEALLOCATE(solverColToSEquationsMap%equationsMatrixNumbers)
    IF(ALLOCATED(solverColToSEquationsMap%equationsColumnNumbers)) &
      & DEALLOCATE(solverColToSEquationsMap%equationsColumnNumbers)
    IF(ALLOCATED(solverColToSEquationsMap%couplingCoefficients)) &
      & DEALLOCATE(solverColToSEquationsMap%couplingCoefficients)
    solverColToSEquationsMap%jacobianColumnNumber=0    
    solverColToSEquationsMap%jacobianCouplingCoefficient=0.0_DP
    
    EXITS("SolverMappingSCToSEquationsMap_Finalise")
    RETURN
999 ERRORS("SolverMappingSCToSEquationsMap_Finalise",err,error)
    EXITS("SolverMappingSCToSEquationsMap_Finalise")
    RETURN 1
    
  END SUBROUTINE SolverMappingSCToSEquationsMap_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the solver column to static equations mapping and deallocates all memory.
  SUBROUTINE SolverMappingSCToSEquationsMap_Initialise(solverColToSEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverColToStaticEquationsMapType) :: solverColToSEquationsMap !<The solver column to static equations map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverMappingSCToSEquationsMap_Initialise",err,error,*999)

    solverColToSEquationsMap%numberOfLinearEquationsMatrices=0
    solverColToSEquationsMap%jacobianColumnNumber=0    
    solverColToSEquationsMap%jacobianCouplingCoefficient=0.0_DP
    
    EXITS("SolverMappingSCToSEquationsMap_Initialise")
    RETURN
999 ERRORS("SolverMappingSCToSEquationsMap_Initialise",err,error)
    EXITS("SolverMappingSCToSEquationsMap_Initialise")
    RETURN 1
    
  END SUBROUTINE SolverMappingSCToSEquationsMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver column to equations set map and deallocates all memory.
  SUBROUTINE SolverMappingSMToESMap_Finalise(solverMatrixToEquationsSetMap,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsSetMapType), POINTER :: solverMatrixToEquationsSetMap !<The solver column to equations set map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx
    
    ENTERS("SolverMappingSMToESMap_Finalise",err,error,*999)

    IF(ASSOCIATED(solverMatrixToEquationsSetMap)) THEN
      IF(ALLOCATED(solverMatrixToEquationsSetMap%solverColToDynamicEquationsMaps)) THEN
        DO columnIdx=1,SIZE(solverMatrixToEquationsSetMap%solverColToDynamicEquationsMaps,1)
          CALL SolverMappingSCToDEquationsMap_Finalise(solverMatrixToEquationsSetMap% &
            & solverColToDynamicEquationsMaps(columnIdx),err,error,*999)
        ENDDO !columnIdx
        DEALLOCATE(solverMatrixToEquationsSetMap%solverColToDynamicEquationsMaps)
        solverMatrixToEquationsSetMap%haveDynamic=.FALSE.
      ENDIF
      IF(ALLOCATED(solverMatrixToEquationsSetMap%solverColToStaticEquationsMaps)) THEN
        DO columnIdx=1,SIZE(solverMatrixToEquationsSetMap%solverColToStaticEquationsMaps,1)
          CALL SolverMappingSCToSEquationsMap_Finalise(solverMatrixToEquationsSetMap% &
            & solverColToStaticEquationsMaps(columnIdx),err,error,*999)
        ENDDO !columnIdx
        DEALLOCATE(solverMatrixToEquationsSetMap%solverColToStaticEquationsMaps)
        solverMatrixToEquationsSetMap%haveStatic=.FALSE.
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

  !>Initialises the solver column to equations set mapping and deallocates all memory.
  SUBROUTINE SolverMappingSMToESMap_Initialise(solverMatrixToEquationsSetMap,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsSetMapType), POINTER :: solverMatrixToEquationsSetMap !<The solver column to equations set map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverMappingSMToESMap_Initialise",err,error,*998)

    IF(ASSOCIATED(solverMatrixToEquationsSetMap)) &
      & CALL FlagError("Solver matrix to equations set map is already associated.",err,error,*998)

    ALLOCATE(solverMatrixToEquationsSetMap,STA7T=err)
    IF(err/=0) CALL FlagError("Could not allocate solver matrix to equations set map.",err,error,*999)    
    NULLIFY(solverMatrixToEquationsSetMap%equations)
    solverMatrixToEquationsSetMap%haveDynamic=.FALSE.
    solverMatrixToEquationsSetMap%haveStatic=.FALSE.
    
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
    TYPE(SolverColToEquationsMapsType), POINTER :: solverMatrixToEquationsMap !<The solver column to equations sets map to finalise
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
        DO interfaceConditionIdx=1,SIZE(solverMatrixToEquationsMap%solverMatrixToInterfaceConditionIdxMaps,1)
          CALL SolverMappingSMToICMap_Finalise(solverMatrixToEquationsMap% &
            & solverMatrixToInterfaceConditionMaps(iterfaceConditionIdx)%ptr,err,error,*999)
        ENDDO !interfaceConditionIdx
        DEALLOCATE(solverMatrixToEquationsMap%solverMatrixToInterfaceConditionMaps)
      ENDIF
      IF(ALLOCATED(solverMatrixToEquationsMap%solverDOFToVariableDOFsMap)) THEN
        DO dofIdx=1,SIZE(solverMatrixToEquationsMap%solverDOFToVariableDOFsMap,1)
          CALLSolverMappingSDOFToVDOFsMap_Finalise(solverMatrixToEquationsMap%solverDOFToVariableDOFsMap(dofIdx), &
            & err,error,*999)
        ENDDO !dofIdx
        DEALLOCATE(solverMatrixToEquationsMap%solverDOFToVariableDOFsMap)
      ENDIF
      CALL SolverMappingSCToECSMap_Finalise(solverMatrixToEquationsMap%solverColToEquationsColsMap,err,error,*999)
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

  !>Initialises the solver column to equations mapping and deallocates all memory.
  SUBROUTINE SolverMappingSMToEquationsMap_Initialise(solverMatrixToEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverColToEquationsMapsType), POINTER :: solverMatrixToEquationsMap !<A pointer to the solver matrix to equations map to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverMappingSMToEquationsMap_Initialise",err,error,*999)

    IF(ASSOCIATED(solverMatrixToEquationsMap)) &
      & CALL FlagError("Solver matrix to equations map is already associated.",err,error,*999)

    ALLOCATE(solverMatrixToEquationsMap%solverMatrixToEquationsMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver matrix to equations map.",err,error,*999)
    NULLIFY(solverMatrixToEquationsMap%solverMapping)
    solverMatrixToEquationsMap%solverMatrixNumber=0
    NULLIFY(solverMatrixToEquationsMap%solverMatrix)
    NULLIFY(solverMatrixToEquationsMap%variablesList)
    solverMatrixToEquationsMap%numberOfColumns=0
    solverMatrixToEquationsMap%numberOfDofs=0
    solverMatrixToEquationsMap%totalNumberOfDofs=0
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
    TYPE(SolverMatrixToInterfaceConditionMapType), POINTER :: solverMatrixToInterfaceConditionMap !<A pointer to the solver matrix to interface condition map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx
    
    ENTERS("SolverMappingSMToICMap_Finalise",err,error,*999)

    IF(ASSOCIATED(solverMatrixToInterfaceConditionMap)) THEN
      IF(ALLOCATED(solverMatrixToInterfaceConditionMap%solverColToInterfaceEquationsMaps)) THEN
        DO columnIdx=1,SIZE(solverMatrixToInterfaceConditionMap%solverColToInterfaceEquationsMaps,1)
          CALL SolverMappingSCToIEquationsMap_Finalise(solverMatrixToInterfaceConditionMap% &
            & solverColToInterfaceEquationsMaps(columnIdx),err,error,*999)
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

  !>Initialises the solver column to interface mapping.
  SUBROUTINE SolverMappingSMToICMap_Initialise(solverMatrixToInterfaceConditionMap,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToInterfaceConditionMapType) :: solverMatrixToInterfaceConditionMap !<The solver column to interface map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverMappingSMToICMap_Initialise",err,error,*999)

    NULLIFY(solverMatrixToInterfaceConditionMap%interfaceEquations)
    
    EXITS("SolverMappingSMToICMap_Initialise")
    RETURN
999 ERRORS("SolverMappingSMToICMap_Initialise",err,error)
    EXITS("SolverMappingSMToICMap_Initialise")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToICMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver column to interface equations map and deallocates all memory.
  SUBROUTINE SolverMappingSCToIEquationsMap_Finalise(solverColToInterfaceEquaitonsMap,err,error,*)

    !Argument variables
    TYPE(SolverColToInterfaceEquationsMapType) :: solverColToInterfaceEquaitonsMap !<The solver column to interface equatiosn map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverMappingSCToIEquationsMap_Finalise",err,error,*999)

    solverColToInterfaceEquaitonsMap%numberOfInterfaceMatrices=0
    IF(ALLOCATED(solverColToInterfaceEquaitonsMap%interfaceMatrixNumbers))  &
      & DEALLOCATE(solverColToInterfaceEquaitonsMap%interfaceMatrixNumbers)
    IF(ALLOCATED(solverColToInterfaceEquaitonsMap%interfaceColumnNumbers))  &
      & DEALLOCATE(solverColToInterfaceEquaitonsMap%interfaceColumnNumbers)
    IF(ALLOCATED(solverColToInterfaceEquaitonsMap%couplingCoefficients))  &
      & DEALLOCATE(solverColToInterfaceEquaitonsMap%couplingCoefficients)
    
    EXITS("SolverMappingSCToIEquationsMap_Finalise")
    RETURN
999 ERRORS("SolverMappingSCToIEquationsMap_Finalise",err,error)
    EXITS("SolverMappingSCToIEquationsMap_Finalise")
    RETURN 1
    
  END SUBROUTINE SolverMappingSCToIEquationsMap_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the solver column to interface equations mapping.
  SUBROUTINE SolverMappingSCToIEquationsMap_Initialise(solverColToInterfaceEquaitonsMap,err,error,*)

    !Argument variables
    TYPE(SolverColToInterfaceEquationsMapType) :: solverColToInterfaceEquaitonsMap !<The solver column to interface equations map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverMappingSCToIEquationsMap_Initialise",err,error,*999)

    solverColToInterfaceEquaitonsMap%numberOfInterfaceMatrices=0
    
    EXITS("SolverMappingSCToIEquationsMap_Initialise")
    RETURN
999 ERRORS("SolverMappingSCToIEquationsMap_Initialise",err,error)
    EXITS("SolverMappingSCToIEquationsMap_Initialise")
    RETURN 1
    
  END SUBROUTINE SolverMappingSCToIEquationsMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver dof to variable mapping and deallocates all memory.
  SUBROUTINE SolverMappingSDOFToVDOFsMap_Finalise(solverDOFToVariableDOFsMap,err,error,*)

    !Argument variables
    TYPE(SolverDOFToVariableDOFsMapType) :: solverDOFToVariableDOFsMap !<The solver dof to variable map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverMappingSDOFToVDOFsMap_Finalise",err,error,*999)

    solverDOFToVariableDOFsMap%numberOfEquationDOFs=0
    IF(ALLOCATED(solverDOFToVariableDOFsMap%equationTypes)) DEALLOCATE(solverDOFToVariableDOFsMap%equationTypes)
    IF(ALLOCATED(solverDOFToVariableDOFsMap%equationIndices)) DEALLOCATE(solverDOFToVariableDOFsMap%equationIndices)
    IF(ALLOCATED(solverDOFToVariableDOFsMap%variable)) DEALLOCATE(solverDOFToVariableDOFsMap%variable)
    IF(ALLOCATED(solverDOFToVariableDOFsMap%variableDOF)) DEALLOCATE(solverDOFToVariableDOFsMap%variableDOF)
    IF(ALLOCATED(solverDOFToVariableDOFsMap%variableCoefficient)) DEALLOCATE(solverDOFToVariableDOFsMap%variableCoefficient)
    IF(ALLOCATED(solverDOFToVariableDOFsMap%additiveConstant)) DEALLOCATE(solverDOFToVariableDOFsMap%additiveConstant)
    
    EXITS("SolverMappingSDOFToVDOFsMap_Finalise")
    RETURN
999 ERRORS("SolverMappingSDOFToVDOFsMap_Finalise",err,error)
    EXITS("SolverMappingSDOFToVDOFsMap_Finalise")
    RETURN 1
    
  END SUBROUTINE SolverMappingSDOFToVDOFsMap_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the solver dof to variable mapping.
  SUBROUTINE SolverMappingSDOFToVDOFs_Initialise(solverDOFToVariableDOFsMap,err,error,*)

    !Argument variables
    TYPE(SolverDOFToVariableDOFsMapType) :: solverDOFToVariableDOFsMap !<The solver dof to variable map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverMappingSDOFToVDOFs_Initialise",err,error,*999)

    solverDOFToVariableDOFsMap%numberOfEquationDOFs=0
    
    EXITS("SolverMappingSDOFToVDOFs_Initialise")
    RETURN
999 ERRORS("SolverMappingSDOFToVDOFs_Initialise",err,error)
    EXITS("SolverMappingSDOFToVDOFs_Initialise")
    RETURN 1
    
  END SUBROUTINE SolverMappingSDOFToVDOFs_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises a solver row to equations map and deallocates all memory.
  SUBROUTINE SolverMappingSRToERSMap_Finalise(solverRowToEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverRowToEquationsMapsType) :: solverRowToEquationsMap !<The solver row to equations map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverMappingSRToERSMap_Finalise",err,error,*999)

    solverRowToEquationsMap%numberOfEquationsSetRows=0
    solverRowToEquationsMap%interfaceConditionIndex=0
    IF(ALLOCATED(solverRowToEquationsMap%equationsIndex)) DEALLOCATE(solverRowToEquationsMap%equationsIndex)
    IF(ALLOCATED(solverRowToEquationsMap%rowColNumber)) DEALLOCATE(solverRowToEquationsMap%rowColNumber)
    IF(ALLOCATED(solverRowToEquationsMap%couplingCoefficients)) DEALLOCATE(solverRowToEquationsMap%couplingCoefficients)
    
    EXITS("SolverMappingSRToERSMap_Finalise")
    RETURN
999 ERRORSEXITS("SolverMappingSRToERSMap_Finalise",err,error)    
    RETURN 1
    
  END SUBROUTINE SolverMappingSRToERSMap_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises a solver row to equations map.
  SUBROUTINE SolverMappingSRToERSMap_Initialise(solverRowToEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverRowToEquationsMapsType) :: solverRowToEquationsMap !<The solver row to equations map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverMappingSRToERSMap_Initialise",err,error,*999)

    solverRowToEquationsMap%numberOfEquationsSetRows=0
    solverRowToEquationsMap%interfaceConditionIndex=0
        
    EXITS("SolverMappingSRToERSMap_Initialise")
    RETURN
999 ERRORS("SolverMappingSRToERSMap_Initialise",err,error)    
    EXITS("SolverMappingSRToERSMap_Initialise")    
    RETURN 1
    
  END SUBROUTINE SolverMappingSRToERSMap_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver mapping variable and deallocates all memory.
  SUBROUTINE SolverMappingVariable_Finalise(solverMappingVariable,err,error,*)

    !Argument variables
    TYPE(SolverMappingVariableType) :: solverMapping_VARIABLE !<The solver mapping variable to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverMappingVariable_Finalise",err,error,*999)

    NULLIFY(solverMappingVariable%variable)
    solverMappingVariable%variableType=0
    solverMappingVariable%numberOfEquations=0
    IF(ALLOCATED(solverMappingVariable%equationTypes)) DEALLOCATE(solverMappingVariable%equationTypes)
    IF(ALLOCATED(solverMappingVariable%equationIndices)) DEALLOCATE(solverMappingVariable%equationIndices)
       
    EXITS("SolverMappingVariable_Finalise")
    RETURN
999 ERRORSEXITS("SolverMappingVariable_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMappingVariable_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the solver mapping and deallocates all memory.
  SUBROUTINE SolverMappingVariable_Initiali7se(solverMappingVariable,err,error,*)

    !Argument variables
    TYPE(SolverMappingVariableType) :: solverMappingVariable !<The solver mapping variable to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingVariable_Initialise",err,error,*999)

    NULLIFY(solverMappingVariable%variable)
    solverMappingVariable%variableType=0
    solverMappingVariable%numberOfEquations=0
    
    EXITS("SolverMappingVariable_Initialise")
    RETURN
999 ERRORSEXITS("SolverMappingVariable_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMappingVariable_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the solver mapping variables and deallocates all memory.
  SUBROUTINE SolverMappingVariables_Finalise(solverMappingVariables,err,error,*)

    !Argument variables
    TYPE(SolverMappingVariablesType) :: solverMappingVariables !<The solver mapping variables to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: variableIdx

    ENTERS("SolverMappingVariables_Finalise",err,error,*999)

    solverMappingVariables%numberOfVariables=0
    DO variableIdx=1,SIZE(solverMappingVariables%variables,1)
      CALL SolverMappingVariable_Finalise(solverMappingVariables%variables(variableIdx),err,error,*999)
    ENDDO !variableIdx
     
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
    TYPE(SolverMappingVariablesType) :: solverMappingVariables !<The solver mapping variables to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingVariables_Initialise",err,error,*999)

    solverMappingVariables%numberOfVariables=0
    
    EXITS("SolverMappingVariables_Initialise")
    RETURN
999 ERRORSEXITS("SolverMappingVariables_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMappingVariables_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises a variable to solver column map and deallocates all memory.
  SUBROUTINE SolverMappingVToSCMap_Finalise(variableToSolverColMap,err,error,*)

    !Argument variables
    TYPE(VariableToSolverColMapType) :: variableToSolverColMap !<The variable to solver column map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverMappingVToSCMap_Finalise",err,error,*999)

    IF(ALLOCATED(variableToSolverColMap%columnNumbers)) DEALLOCATE(variableToSolverColMap%columnNumbers)
    IF(ALLOCATED(variableToSolverColMap%couplingCoefficients)) DEALLOCATE(variableToSolverColMap%couplingCoefficients)
    IF(ALLOCATED(variableToSolverColMap%additiveConstants)) DEALLOCATE(variableToSolverColMap%additiveConstants)
        
    EXITS("SolverMappingVToSCMap_Finalise")
    RETURN
999 ERRORS("SolverMappingVToSCMap_Finalise",err,error)    
    EXITS("SolverMappingVToSCMap_Finalise")    
    RETURN 1
   
  END SUBROUTINE SolverMappingVToSCMap_Finalise
  
  !
  !================================================================================================================================
  !

  !>Initialises a variable to solver column map.
  SUBROUTINE SolverMappingVToSCMap_Initialise(variableToSolverColMap,err,error,*)

    !Argument variables
    TYPE(VariableToSolverColMapType) :: variableToSolverColMap !<The variable to solver column map to initalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverMappingVToSCMap_Initialise",err,error,*999)

    !Nothing to do here, all members are allocatable

    EXITS("SolverMappingVToSCMap_Initialise")
    RETURN
999 ERRORS("SolverMappingVToSCMap_Initialise",err,error)    
    EXITS("SolverMappingVToSCMap_Initialise")    
    RETURN 1
   
  END SUBROUTINE SolverMappingVToSCMap_Initialise

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

  !>Finalise the solver matrix to equations map and deallocate all memory
  SUBROUTINE SolverMappingSMToEquationsMap_Finalise(solverMatrixToEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap !<The solver matrix to equations map to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

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
        DO interfaceConditionIdx=1,SIZE(solverMatrixToEquationsMap%solverMatrixToInterfaceConditionToEquationsMaps,1)
          CALL SolverMappingSMToICMap_Finalise(solverMatrixToEquationsMap% &
            & solverMatrixToInterfaceConditionMaps(interfaceConditionIdx)%ptr,err,error,*999)
        ENDDO !interfaceConditionIdx
        DEALLOCATE(solverMatrixToEquationsMap%solverMatrixToEquationsSetMaps)
      ENDIF
      IF(ALLOCATED(solverMatrixToEquationsMap%solverDOFToVariableDOFsMap,1)) THEN
        DO dofIdx=1,SIZE(solverMatrixToEquationsMap%,solverDOFToVariableDOFsMap,1)
          CALL SolverMappingSDOFToVDOFsMap_Finalise(solverMatrixToEquationsMap%solverDOFToVariableDOFsMap1(dofIdx),err,error,*999)
        ENDDO !interfaceConditionIdx
        DEALLOCATE(solverMatrixToEquationsMap%solverDOFToVariableDOFsMap)
      ENDIF
      CALL SolverMappingSCToECSMap_Finalise(solverMatrixToEquationsMap%solverColToEquationsColsMap,err,error,*999)
      CALL DomainMapping_Finalsie(solverMatrixToEquationsMap%columnDOFSMapping,err,error,*999)
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

  !>Allocate and initialise the solver matrix to equations map
  SUBROUTINE SolverMappingSMToEquationsMap_Initialise(solverMatrixToEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap !<The solver matrix to equations map to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverMappingSMToEquationsMap_Initialise",err,error,*998)

    IF(ASSOCIATED(solverMatrixToEquationsMap)) &
      & CALL FlagError("Solver mapping to equations map is already associated.",err,error,*998)

    ALLOCATE(solverMatrixToEquationsMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocated solver matrix to equations map.",err,error,*999)
    NULLIFY(solverMatrixToEquationsMap%solverMapping)
    solverMatrixToEquationsMap%solverMatrixNumber=0
    NULLIFY(solverMatrixToEquationsMap%solverMatrix)
    NULLIFY(solverMatrixToEquationsMap%variablesList)
    solverMatrixToEquationsMap%numberOfColumns=0
    solverMatrixToEquationsMap%numberOfDOFs=0
    solverMatrixToEquationsMap%totalNumberOfDOFs=0
    solverMatrixToEquationsMap%numberOfGlobalDOFs=0
    solverMatrixToEquationsMap%numberOfEquationsSets=0
    solverMatrixToEquationsMap%numberOfInterfaceConditions=0
    NULLIFY(solverMatrixToEquationsMap%solverColToEquationsColsMap)
    NULLIFY(solverMatrixToEquationsMap%columnDOFSMapping)    

    EXITS("SolverMappingSMToEquationsMap_Initialise")
    RETURN
999 CALL SolverMappingSMToEquationsMap_Finalise(solverMatrixToEquationsMap,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverMappingSMToEquationsMap_Initialise",err,error)
    RETURN 1

  END SUBROUTINE SolverMappingSMToEquationsMap_Initialise
  
  !
  !================================================================================================================================
  !

END MODULE SolverMappingRoutines
