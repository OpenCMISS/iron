!> \file
!> \author Ting Yu
!> \brief This module set the boundary conditions for the given equation set
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
!> Contributor(s): Ting Yu, Chris Bradley
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

!>This module handles all boundary conditions routines.
MODULE BoundaryConditionsRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BoundaryConditionAccessRoutines
  USE CmissMPI
  USE ComputationRoutines
  USE ComputationAccessRoutines
  USE Constants
  USE CoordinateSystemRoutines
  USE CoordinateSystemAccessRoutines
  USE DecompositionAccessRoutines
  USE DistributedMatrixVector
  USE DistributedMatrixVectorAccessRoutines
  USE DomainMappings
  USE EquationsAccessRoutines
  USE EquationsMappingAccessRoutines
  USE EquationsMatricesAccessRoutines
  USE EquationsSetAccessRoutines
  USE InterfaceConditionsAccessRoutines
  USE InterfaceEquationsAccessRoutines
  USE InterfaceMappingAccessRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE MeshAccessRoutines
#ifndef NOMPIMOD
  USE MPI
#endif
  USE NodeRoutines
  USE SolverAccessRoutines
  USE Strings
  USE Timer
  USE Types
  USE Lists
  USE LINKEDLIST_ROUTINES

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

  !>Adds to the value of the specified local DOF and sets this as a boundary condition on the specified local DOF.
  INTERFACE BoundaryConditions_AddLocalDOF
    MODULE PROCEDURE BoundaryConditions_AddLocalDOF0
    MODULE PROCEDURE BoundaryConditions_AddLocalDOF1
  END INTERFACE BoundaryConditions_AddLocalDOF

  INTERFACE BoundaryConditions_AddElement
    MODULE PROCEDURE BoundaryConditions_AddElement
  END INTERFACE BoundaryConditions_AddElement

  INTERFACE BoundaryConditions_AddNode
    MODULE PROCEDURE BoundaryConditions_AddNode
  END INTERFACE BoundaryConditions_AddNode

  INTERFACE BoundaryConditions_SetConstant
    MODULE PROCEDURE BoundaryConditions_SetConstant
  END INTERFACE BoundaryConditions_SetConstant

  !>Sets a boundary condition on the specified local DOF.
  INTERFACE BoundaryConditions_SetLocalDOF
    MODULE PROCEDURE BoundaryConditions_SetLocalDOF0
    MODULE PROCEDURE BoundaryConditions_SetLocalDOF1
  END INTERFACE BoundaryConditions_SetLocalDOF

  INTERFACE BoundaryConditions_SetElement
    MODULE PROCEDURE BoundaryConditions_SetElement
  END INTERFACE BoundaryConditions_SetElement

  INTERFACE BoundaryConditions_SetNode
    MODULE PROCEDURE BoundaryConditions_SetNode
  END INTERFACE BoundaryConditions_SetNode

  PUBLIC BoundaryConditions_CreateFinish,BoundaryConditions_CreateStart

  PUBLIC BoundaryConditions_Destroy

  PUBLIC BoundaryConditions_AddConstant

  PUBLIC BoundaryConditions_AddLocalDOF

  PUBLIC BoundaryConditions_AddElement

  PUBLIC BoundaryConditions_AddNode

  PUBLIC BoundaryConditions_VariableExists

  PUBLIC BoundaryConditions_VariableGet

  PUBLIC BoundaryConditions_SetConstant

  PUBLIC BoundaryConditions_SetLocalDOF
  
  PUBLIC BoundaryConditions_SetElement

  PUBLIC BoundaryConditions_SetNode

  PUBLIC BoundaryConditions_NeumannIntegrate

  PUBLIC BoundaryConditions_NeumannSparsityTypeSet

  PUBLIC BoundaryConditions_ConstrainNodeDofsEqual

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finish the creation of boundary conditions.
  SUBROUTINE BoundaryConditions_CreateFinish(boundaryConditions,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: mpiIError,sendCount,storageType, numberOfNonZeros, numberOfRows,count
    INTEGER(INTG) :: variableIdx,dofIdx, equationsMatrixIdx, dirichletIdx, rowIdx, dummy, last, dirichletDOF
    INTEGER(INTG) :: colIdx,equationsSetIdx,parameterSetIdx
    INTEGER(INTG) :: pressureIdx,neumannIdx,numberOfGroupComputationNodes,myGroupComputationNodeNumber,groupCommunicator
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionVariable
    TYPE(DomainMappingType), POINTER :: domainMapping
    TYPE(FieldVariableType), POINTER :: fieldVariable,matrixVariable
    TYPE(BoundaryConditionsDirichletType), POINTER :: boundaryConditionsDirichlet
    TYPE(BoundaryConditionsPressureIncrementedType), POINTER :: boundaryConditionsPressureIncremented
    TYPE(VARYING_STRING) :: localError
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverType), POINTER :: solver
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(BoundaryConditionsSparsityIndicesType), POINTER :: sparsityIndices
    TYPE(ListType), POINTER :: sparseIndices
    TYPE(LinkedList),POINTER :: list(:)
    INTEGER(INTG), ALLOCATABLE:: columnArray(:)
    TYPE(WorkGroupType), POINTER :: workGroup

    ENTERS("BoundaryConditions_CreateFinish",err,error,*999)

    NULLIFY(boundaryConditionsPressureIncremented)

    CALL BoundaryConditions_AssertNotFinished(boundaryConditions,err,error,*999)
    NULLIFY(solverEquations)
    CALL BoundaryConditions_SolverEquationsGet(boundaryConditions,solverEquations,err,error,*999)
    NULLIFY(solver)
    CALL SolverEquations_SolverGet(solverEquations,solver,err,error,*999)
    IF(ALLOCATED(boundaryConditions%boundaryConditionsVariables)) THEN
      NULLIFY(workGroup)
      CALL Solver_WorkGroupGet(solver,workGroup,err,error,*999)
      CALL WorkGroup_GroupCommunicatorGet(workGroup,groupCommunicator,err,error,*999)
      CALL WorkGroup_NumberOfGroupNodesGet(workGroup,numberOfGroupComputationNodes,err,error,*999)
      CALL WorkGroup_GroupNodeNumberGet(workGroup,myGroupComputationNodeNumber,err,error,*999)
      IF(numberOfGroupComputationNodes>0) THEN
        !Transfer all the boundary conditions to all the computation nodes.
        !\todo Look at this.
        DO variableIdx=1,boundaryConditions%numberOfBoundaryConditionsVariables
          boundaryConditionVariable=>boundaryConditions%boundaryConditionsVariables(variableIdx)%PTR
          IF(ASSOCIATED(boundaryConditionVariable)) THEN
            fieldVariable=>boundaryConditionVariable%VARIABLE
            IF(ASSOCIATED(fieldVariable)) THEN
              domainMapping=>fieldVariable%domainMapping
              IF(ASSOCIATED(domainMapping)) THEN
                sendCount=domainMapping%numberOfGlobal
                IF(numberOfGroupComputationNodes>1) THEN
                  !\todo This operation is a little expensive as we are doing an unnecessary sum across all the ranks in order to combin
                  !\todo the data from each rank into all ranks. We will see how this goes for now.
                  CALL MPI_ALLREDUCE(MPI_IN_PLACE,boundaryConditionVariable%DOFTypes, &
                    & sendCount,MPI_INTEGER,MPI_SUM,groupCommunicator,mpiIError)
                  CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",mpiIError,err,error,*999)
                  CALL MPI_ALLREDUCE(MPI_IN_PLACE,boundaryConditionVariable%conditionTypes, &
                    & sendCount,MPI_INTEGER,MPI_SUM,groupCommunicator,mpiIError)
                  CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",mpiIError,err,error,*999)
                ENDIF !mpi_in_place bug workaround - only do this when num comp nodes > 1

              ELSE
                localError="Field variable domain mapping is not associated for variable type "// &
                  & TRIM(NumberToVString(variableIdx,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF

              IF(numberOfGroupComputationNodes>1) THEN

                ! Update the total number of boundary condition types by summing across all nodes
                CALL MPI_ALLREDUCE(MPI_IN_PLACE,boundaryConditionVariable%dofCounts, &
                  & MAX_BOUNDARY_CONDITION_NUMBER,MPI_INTEGER,MPI_SUM,groupCommunicator,mpiIError)
                CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",mpiIError,err,error,*999)
                CALL MPI_ALLREDUCE(MPI_IN_PLACE,boundaryConditionVariable%numberOfDirichletConditions, &
                  & 1,MPI_INTEGER,MPI_SUM,groupCommunicator,mpiIError)
                CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",mpiIError,err,error,*999)
              ENDIF !mpi_in_place bug workaround - only do this when num comp nodes > 1

              ! Check that the boundary conditions set are appropriate for equations sets
              CALL BoundaryConditions_CheckEquations(boundaryConditionVariable,err,error,*999)

              IF(numberOfGroupComputationNodes>1) THEN
                !Make sure the required parameter sets are created on all computational nodes and begin updating them
                CALL MPI_ALLREDUCE(MPI_IN_PLACE,boundaryConditionVariable%parameterSetRequired, &
                  & FIELD_NUMBER_OF_SET_TYPES,MPI_LOGICAL,MPI_LOR,groupCommunicator,mpiIError)
                CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",mpiIError,err,error,*999)
                DO parameterSetIdx=1,FIELD_NUMBER_OF_SET_TYPES
                  IF(boundaryConditionVariable%parameterSetRequired(parameterSetIdx)) THEN
                    CALL Field_ParameterSetEnsureCreated(fieldVariable%FIELD,fieldVariable%variableType, &
                      & parameterSetIdx,err,error,*999)
                    CALL Field_ParameterSetUpdateStart(fieldVariable%FIELD,fieldVariable%variableType, &
                      & parameterSetIdx,err,error,*999)
                  END IF
                END DO
              ENDIF !mpi_in_place bug workaround - only do this when num comp nodes > 1

              ! Set up pressure incremented condition, if it exists
              IF(boundaryConditionVariable%dofCounts(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)>0) THEN
                CALL BoundaryConditions_PressureIncrementedInitialise(boundaryConditionVariable,err,error,*999)
                boundaryConditionsPressureIncremented=>boundaryConditionVariable%pressureIncrementedBoundaryConditions
              END IF

              ! Set up Neumann condition information if there are any Neumann conditions
              IF(boundaryConditionVariable%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT)>0.OR. &
                & boundaryConditionVariable%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)>0) THEN
                CALL BoundaryConditions_NeumannInitialise(boundaryConditionVariable,err,error,*999)
              END IF

              ! Loop over all global DOFs, keeping track of the dof indices of specific BC types where required
              pressureIdx=1
              neumannIdx=1
              DO dofIdx=1,fieldVariable%numberOfGlobalDofs
                IF(boundaryConditionVariable%conditionTypes(dofIdx)== BOUNDARY_CONDITION_PRESSURE_INCREMENTED) THEN
                  boundaryConditionsPressureIncremented%pressureIncrementedDOFIndices(pressureIdx)=dofIdx
                  pressureIdx=pressureIdx+1
                ELSE IF(boundaryConditionVariable%conditionTypes(dofIdx)==BOUNDARY_CONDITION_NEUMANN_POINT.OR. &
                  & boundaryConditionVariable%conditionTypes(dofIdx)==BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) THEN
                  boundaryConditionVariable%neumannBoundaryConditions%setDofs(neumannIdx)=dofIdx
                  neumannIdx=neumannIdx+1
                END IF
              END DO

              ! Now that we know where Neumann point DOFs are, we can calculate matrix structure
              IF(boundaryConditionVariable%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT)>0.OR. &
                & boundaryConditionVariable%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)>0) THEN
                CALL BoundaryConditions_NeumannMatricesInitialise(boundaryConditionVariable,err,error,*999)
              END IF

              ! Check that there is at least one dirichlet condition
              IF(boundaryConditionVariable%numberOfDirichletConditions>0) THEN
                CALL BoundaryConditions_DirichletInitialise(boundaryConditionVariable,err,error,*999)
                boundaryConditionsDirichlet=>boundaryConditionVariable%dirichletBoundaryConditions
                IF(ASSOCIATED(boundaryConditionsDirichlet)) THEN
                  ! Find dirichlet conditions
                  dirichletIdx=1
                  DO dofIdx=1,fieldVariable%numberOfGlobalDofs
                    IF(boundaryConditionVariable%DOFTypes(dofIdx)==BOUNDARY_CONDITION_DOF_FIXED) THEN
                      boundaryConditionsDirichlet%dirichletDOFIndices(dirichletIdx)=dofIdx
                      dirichletIdx=dirichletIdx+1
                    ENDIF
                  ENDDO

                  !Store Dirichlet dof indices
                  solverEquations=>boundaryConditions%solverEquations
                  IF(ASSOCIATED(solverEquations)) THEN
                    IF(ASSOCIATED(solverEquations%solverMapping)) THEN
                      DO equationsSetIdx=1,solverEquations%solverMapping%numberOfEquationsSets
                        equationsSet=>solverEquations%solverMapping%equationsSets(equationsSetIdx)%PTR
                        IF(ASSOCIATED(equationsSet)) THEN
                          NULLIFY(equations)
                          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
                          NULLIFY(vectorEquations)
                          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                          NULLIFY(vectorMatrices)
                          CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
                          NULLIFY(vectorMapping)
                          CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
                          linearMapping=>vectorMapping%linearMapping
                          IF(ASSOCIATED(linearMapping)) THEN
                            NULLIFY(linearMatrices)
                            CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
                            !Iterate through equations matrices
                            DO equationsMatrixIdx=1,linearMatrices%numberOfLinearMatrices
                              NULLIFY(equationsMatrix)
                              CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,equationsMatrixIdx,equationsMatrix, &
                                & err,error,*999)
                              matrixVariable=>linearMapping%equationsMatrixToVarMaps(equationsMatrixIdx)%variable
                              IF(ASSOCIATED(matrixVariable,fieldVariable)) THEN
                                IF(boundaryConditionVariable%numberOfDirichletConditions>0) THEN
                                  CALL DistributedMatrix_TransposeRowsColumnsSet(equationsMatrix%MATRIX, &
                                    & boundaryConditionsDirichlet%dirichletDOFIndices(1: &
                                    & boundaryConditionVariable%numberOfDirichletConditions),err,error,*999)
                                ENDIF
                                CALL DistributedMatrix_StorageTypeGet(equationsMatrix%MATRIX,storageType,err,error,*999)
                                SELECT CASE(storageType)
                                CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
                                  !Do nothing
                                CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                                  !Do nothing
                                CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                                  CALL FlagError("Not implemented for column major storage.",err,error,*999)
                                CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                                  CALL FlagError("Not implemented for row major storage.",err,error,*999)
                                CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                                  !Get Sparsity pattern, number of non zeros, number of rows
                                  !CALL DistributedMatrix_StorageLocationsGet(equationsMatrix%MATRIX,ROW_INDICES, &
                                  !  & COLUMN_INDICES,err,error,*999)
                                  CALL DistributedMatrix_NumberOfNonZerosGet(equationsMatrix%MATRIX,numberOfNonZeros, &
                                    & err,error,*999)
                                  !Get the matrix stored as a linked list
                                  CALL DistributedMatrix_LinkListGet(equationsMatrix%MATRIX,list,err,error,*999)
                                  numberOfRows=vectorMatrices%totalNumberOfRows
                                  !Initialise sparsity indices arrays
                                  CALL BoundaryConditions_SparsityIndicesInitialise(boundaryConditionsDirichlet% &
                                    & linearSparsityIndices(equationsSetIdx,equationsMatrixIdx)%PTR, &
                                    & boundaryConditionVariable%numberOfDirichletConditions,err,error,*999)
                                  !Find dirichlet columns and store the non zero indices (with respect to the 1D storage array)
                                  NULLIFY(sparsityIndices)
                                  sparsityIndices=>boundaryConditionsDirichlet%linearSparsityIndices( &
                                    & equationsSetIdx,equationsMatrixIdx)%PTR
                                  IF(ASSOCIATED(sparsityIndices)) THEN
                                    !Setup list for storing dirichlet non zero indices
                                    NULLIFY(sparseIndices)
                                    CALL LIST_CREATE_START(sparseIndices,err,error,*999)
                                    CALL LIST_DATA_TYPE_SET(sparseIndices,LIST_INTG_TYPE,err,error,*999)
                                    CALL LIST_INITIAL_SIZE_SET(sparseIndices, &
                                      & boundaryConditionVariable%numberOfDirichletConditions*( &
                                      & numberOfNonZeros/numberOfRows),err,error,*999)
                                    CALL LIST_CREATE_FINISH(sparseIndices,err,error,*999)
                                    count=0
                                    sparsityIndices%sparseColumnIndices(1)=1
                                    last=1
                                    DO dirichletIdx=1,boundaryConditionVariable%numberOfDirichletConditions
                                      dirichletDOF=boundaryConditionsDirichlet%dirichletDOFIndices(dirichletIdx)
                                      CALL LinkedList_to_Array(list(dirichletDOF),columnArray,err,error,*999)
                                      DO rowIdx=1,SIZE(columnArray)
                                        CALL LIST_ITEM_ADD(sparseIndices,columnArray(rowIdx),err,error,*999)
                                        count=count+1
                                        last=rowIdx+1
                                      ENDDO
                                      sparsityIndices%sparseColumnIndices(dirichletIdx+1)=count+1
                                    ENDDO
                                    CALL LIST_DETACH_AND_DESTROY(sparseIndices,dummy,sparsityIndices%sparseRowIndices, &
                                      & err,error,*999)
                                    DO colIdx =1,numberOfRows
                                      CALL LINKEDLIST_DESTROY(list(colIdx),err,error,*999)
                                    ENDDO
                                  ELSE
                                    localError="Sparsity indices arrays are not associated for this equations matrix."
                                    CALL FlagError(localError,err,error,*999)
                                  ENDIF
                                CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                                  CALL FlagError("Not implemented for compressed column storage.",err,error,*999)
                                CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                                  CALL FlagError("Not implemented for row column storage.",err,error,*999)
                                CASE DEFAULT
                                  localError="The storage type of "//TRIM(NumberToVString(storageType,"*",err,error)) &
                                    //" is invalid."
                                  CALL FlagError(localError,err,error,*999)
                                END SELECT
                              ENDIF
                            ENDDO
                          ENDIF

                          dynamicMapping=>vectorMapping%dynamicMapping
                          IF(ASSOCIATED(dynamicMapping)) THEN
                            NULLIFY(dynamicMatrices)
                            CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
                            !Iterate through equations matrices
                            DO equationsMatrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
                              NULLIFY(equationsMatrix)
                              CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,equationsMatrixIdx,equationsMatrix, &
                                & err,error,*999)
                              matrixVariable=>dynamicMapping%equationsMatrixToVarMaps(equationsMatrixIdx)%variable
                              IF(ASSOCIATED(matrixVariable,fieldVariable)) THEN
                                IF(boundaryConditionVariable%numberOfDirichletConditions>0) THEN
                                  CALL DistributedMatrix_TransposeRowsColumnsSet(equationsMatrix%MATRIX, &
                                    & boundaryConditionsDirichlet%dirichletDOFIndices(1: &
                                    & boundaryConditionVariable%numberOfDirichletConditions),err,error,*999)
                                ENDIF
                                CALL DistributedMatrix_StorageTypeGet(equationsMatrix%MATRIX,storageType,err,error,*999)
                                SELECT CASE(storageType)
                                CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
                                  !Do nothing
                                CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                                  !Do nothing
                                CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                                  CALL FlagError("Not implemented for column major storage.",err,error,*999)
                                CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                                  CALL FlagError("Not implemented for row major storage.",err,error,*999)
                                CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                                  !Get Sparsity pattern, number of non zeros, number of rows
                                  !CALL DistributedMatrix_StorageLocationsGet(equationsMatrix%MATRIX,ROW_INDICES, &
                                  !  & COLUMN_INDICES,err,error,*999)
                                  CALL DistributedMatrix_NumberOfNonZerosGet(equationsMatrix%MATRIX,numberOfNonZeros, &
                                    & err,error,*999)
                                  !Sparse matrix in a list
                                  CALL DistributedMatrix_LinkListGet(equationsMatrix%MATRIX,list,err,error,*999)
                                  numberOfRows=vectorMatrices%totalNumberOfRows
                                  !Intialise sparsity indices arrays
                                  CALL BoundaryConditions_SparsityIndicesInitialise(boundaryConditionsDirichlet% &
                                    & dynamicSparsityIndices(equationsSetIdx,equationsMatrixIdx)%PTR, &
                                    & boundaryConditionVariable%numberOfDirichletConditions,err,error,*999)
                                  !Find dirichlet columns and store the non zero indices (with respect to the 1D storage array)
                                  NULLIFY(sparsityIndices)
                                  sparsityIndices=>boundaryConditionsDirichlet%dynamicSparsityIndices( &
                                    & equationsSetIdx,equationsMatrixIdx)%PTR
                                  IF(ASSOCIATED(sparsityIndices)) THEN
                                    ! Setup list for storing dirichlet non zero indices
                                    NULLIFY(sparseIndices)
                                    CALL LIST_CREATE_START(sparseIndices,err,error,*999)
                                    CALL LIST_DATA_TYPE_SET(sparseIndices,LIST_INTG_TYPE,err,error,*999)
                                    CALL LIST_INITIAL_SIZE_SET(sparseIndices, &
                                      & boundaryConditionVariable%numberOfDirichletConditions*( &
                                      & numberOfNonZeros/numberOfRows),err,error,*999)
                                    CALL LIST_CREATE_FINISH(sparseIndices,err,error,*999)
                                    count=0
                                    sparsityIndices%sparseColumnIndices(1)=1
                                    last=1
                                    DO dirichletIdx=1,boundaryConditionVariable%numberOfDirichletConditions
                                      !Dirichlet columns
                                      dirichletDOF=boundaryConditionsDirichlet%dirichletDOFIndices(dirichletIdx)
                                      CALL LinkedList_to_Array(list(dirichletDOF),columnArray,err,error,*999)
                                      !The row indices
                                      DO rowIdx=1,SIZE(columnArray)
                                        CALL LIST_ITEM_ADD(sparseIndices,columnArray(rowIdx),err,error,*999)
                                        count=count+1
                                        last=rowIdx+1
                                      ENDDO
                                      sparsityIndices%sparseColumnIndices(dirichletIdx+1)=count+1
                                    ENDDO
                                    CALL LIST_DETACH_AND_DESTROY(sparseIndices,dummy,sparsityIndices%sparseRowIndices, &
                                      & err,error,*999)
                                    DO colIdx =1,numberOfRows
                                      CALL LINKEDLIST_DESTROY(list(colIdx),err,error,*999)
                                    ENDDO
                                  ELSE
                                    localError="Sparsity indices arrays are not associated for this equations matrix."
                                    CALL FlagError(localError,err,error,*999)
                                  ENDIF
                                CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                                  CALL FlagError("Not implemented for compressed column storage.",err,error,*999)
                                CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                                  CALL FlagError("Not implemented for row column storage.",err,error,*999)
                                CASE DEFAULT
                                  localError="The storage type of "//TRIM(NumberToVString(storageType,"*",err,error)) &
                                    //" is invalid."
                                  CALL FlagError(localError,err,error,*999)
                                END SELECT
                              ENDIF
                            ENDDO
                          ENDIF
                        ELSE
                          localError="Equations Set is not associated for boundary conditions variable "// &
                            & TRIM(NumberToVString(variableIdx,"*",err,error))//"."
                          CALL FlagError(localError,err,error,*999)
                        ENDIF
                      ENDDO !equationsSetIdx
                      !\todo Update interface sparsity structure calculate first then update code below.
                      !                          !Loop over interface conditions. Note that only linear interface matrices implemented so far.
                      !                          DO interface_condition_idx=1,solverEquations%solverMapping%numberOfInterfaceConditions
                      !                            INTERFACE_CONDITION=>solverEquations%solverMapping%interfaceConditions(interface_condition_idx)%PTR
                      !                            IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
                      !                              INTERFACE_EQUATIONS=>INTERFACE_CONDITION%interfaceEquations
                      !                              IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
                      !                                INTERFACE_MATRICES=>INTERFACE_EQUATIONS%interfaceMatrices
                      !                                IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
                      !                                  !Iterate through equations matrices
                      !                                  DO interface_matrix_idx=1,INTERFACE_MATRICES%numberOfInterfaceMatrices
                      !                                    INTERFACE_MATRIX=>INTERFACE_MATRICES%MATRICES(interface_matrix_idx)%PTR
                      !                                    IF(ASSOCIATED(INTERFACE_MATRIX)) THEN
                      !                                      CALL DistributedMatrix_StorageTypeGet(INTERFACE_MATRIX%MATRIX,storageType,err,error,*999)
                      !                                      SELECT CASE(storageType)
                      !                                      CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
                      !                                        !Do nothing
                      !                                      CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                      !                                        !Do nothing
                      !                                      CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                      !                                        CALL FlagError("Not implemented for column major storage.",err,error,*999)
                      !                                      CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                      !                                        CALL FlagError("Not implemented for row major storage.",err,error,*999)
                      !                                      CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                      !                                        !Get Sparsity pattern, number of non zeros, number of rows
                      !                                        CALL DistributedMatrix_StorageLocationsGet(INTERFACE_MATRIX%MATRIX,ROW_INDICES, &
                      !                                          & COLUMN_INDICES,err,error,*999)
                      !                                        CALL DistributedMatrix_NumberOfNonZerosGet(INTERFACE_MATRIX%MATRIX,numberOfNonZeros, &
                      !                                          & err,error,*999)
                      !                                        !Get the matrix stored as a linked list
                      !                                        CALL DistributedMatrix_LinkListGet(INTERFACE_MATRIX%MATRIX,list,err,error,*999)
                      !                                        numberOfRows=vectorMatrices%totalNumberOfRows
                      !                                        !Initialise sparsity indices arrays
                      !                                        CALL BoundaryConditions_SparsityIndicesInitialise(boundaryConditionsDirichlet% &
                      !                                          & linearSparsityIndices(interface_condition_idx,interface_matrix_idx)%PTR, &
                      !                                          & boundaryConditionVariable%numberOfDirichletConditions,err,error,*999)
                      !                                        !Find dirichlet columns and store the non zero indices (with respect to the 1D storage array)
                      !                                        NULLIFY(sparsityIndices)
                      !                                        sparsityIndices=>boundaryConditionsDirichlet%linearSparsityIndices( &
                      !                                            & interface_condition_idx,interface_matrix_idx)%PTR
                      !                                        IF(ASSOCIATED(sparsityIndices)) THEN
                      !                                          !Setup list for storing dirichlet non zero indices
                      !                                          NULLIFY(sparseIndices)
                      !                                          CALL LIST_CREATE_START(sparseIndices,err,error,*999)
                      !                                          CALL LIST_DATA_TYPE_SET(sparseIndices,LIST_INTG_TYPE,err,error,*999)
                      !                                          CALL LIST_INITIAL_SIZE_SET(sparseIndices, &
                      !                                            & numberOfDirichletConditions*(numberOfNonZeros/numberOfRows),err,error,*999)
                      !                                          CALL LIST_CREATE_FINISH(sparseIndices,err,error,*999)
                      !                                          count=0
                      !                                          sparsityIndices%sparseColumnIndices(1)=1
                      !                                          last=1
                      !                                          DO dirichletIdx=1,boundaryConditionVariable%numberOfDirichletConditions
                      !                                            dirichletDOF=boundaryConditionsDirichlet%dirichletDOFIndices(dirichletIdx)
                      !                                            CALL LinkedList_to_Array(list(dirichletDOF),columnArray)
                      !                                              DO rowIdx=1,size(columnArray)
                      !                                                CALL LIST_ITEM_ADD(sparseIndices,columnArray(rowIdx),err,error,*999)
                      !                                                count=count+1
                      !                                                last=rowIdx+1
                      !                                              ENDDO
                      !                                            sparsityIndices%sparseColumnIndices(dirichletIdx+1)=count+1
                      !                                          ENDDO
                      !                                          CALL LIST_DETACH_AND_DESTROY(sparseIndices,dummy,sparsityIndices%sparseRowIndices, &
                      !                                            & err,error,*999)
                      !                                          DO colIdx =1,numberOfRows
                      !                                            CALL LINKEDLIST_DESTROY(list(colIdx))
                      !                                          ENDDO
                      !                                        ELSE
                      !                                          localError="Sparsity indices arrays are not associated for this interface matrix."
                      !                                          CALL FlagError(localError,err,error,*999)
                      !                                        ENDIF
                      !                                      CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                      !                                        CALL FlagError("Not implemented for compressed column storage.",err,error,*999)
                      !                                      CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                      !                                        CALL FlagError("Not implemented for row column storage.",err,error,*999)
                      !                                      CASE DEFAULT
                      !                                        localError="The storage type of "//TRIM(NumberToVString(storageType,"*",err,error)) &
                      !                                          //" is invalid."
                      !                                        CALL FlagError(localError,err,error,*999)
                      !                                      END SELECT
                      !                                    ELSE
                      !                                      CALL FlagError("The interface matrix is not associated.",err,error,*999)
                      !                                    ENDIF
                      !                                  ENDDO
                      !                                ELSE
                      !                                  localError="Interface matrices is not associated for these interface equations."
                      !                                  CALL FlagError(localError,err,error,*999)
                      !                                ENDIF
                      !                              ELSE
                      !                                localError="Interface equations is not associated for this interface condition."
                      !                                CALL FlagError(localError,err,error,*999)
                      !                              ENDIF
                      !                            ELSE
                      !                              localError="Interface condition is not associated for boundary conditions variable "// &
                      !                                & TRIM(NumberToVString(variableIdx,"*",err,error))//"."
                      !                              CALL FlagError(localError,err,error,*999)
                      !                            ENDIF
                      !                          ENDDO !interface_condition_idx
                    ELSE
                      localError="Solver equations solver mapping is not associated."
                      CALL FlagError(localError,err,error,*999)
                    ENDIF
                  ELSE
                    localError="Solver equations is not associated."
                    CALL FlagError(localError,err,error,*999)
                  ENDIF
                ELSE
                  localError="Dirichlet Boundary Conditions type is not associated for boundary condition variable type "// &
                    & TRIM(NumberToVString(variableIdx,"*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
              ENDIF
              ! Finish field update
              DO parameterSetIdx=1,FIELD_NUMBER_OF_SET_TYPES
                IF(boundaryConditionVariable%parameterSetRequired(parameterSetIdx)) THEN
                  CALL Field_ParameterSetUpdateFinish(fieldVariable%FIELD,fieldVariable%variableType, &
                    & parameterSetIdx,err,error,*999)
                END IF
              END DO

              !Finish creating the boundary conditions DOF constraints
              CALL BoundaryConditions_DofConstraintsCreateFinish(boundary_condition_variable,err,error,*999)
            ELSE
              localError="Field variable is not associated for variable index "// &
                & TRIM(NumberToVString(variableIdx,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Boundary conditions variable is not associated for variable index "// &
              & TRIM(NumberToVString(variableIdx,"*",err,error))//".",err,error,*999)
          ENDIF
        ENDDO ! variableIdx

      ENDIF
      !Set the finished flag
      boundaryConditions%boundaryConditionsFinished=.TRUE.
    ELSE
      CALL FlagError("Boundary conditions variables array is not allocated.",err,error,*999)
    ENDIF
    IF(DIAGNOSTICS1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Boundary conditions:",err,error,*999)
      DO variableIdx=1,boundaryConditions%numberOfBoundaryConditionsVariables
        boundaryConditionVariable=>boundaryConditions%boundaryConditionsVariables(variableIdx)%PTR
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Variable type = ",boundaryConditionVariable%variableType, &
            & err,error,*999)
        IF(ASSOCIATED(boundaryConditionVariable)) THEN
          fieldVariable=>boundaryConditionVariable%VARIABLE
          domainMapping=>fieldVariable%domainMapping
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of global dofs = ",domainMapping% &
            & numberOfGlobal,err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,domainMapping%numberOfGlobal,8,8, &
            & boundaryConditionVariable%conditionTypes,'("    Global BCs:",8(X,I8))','(15X,8(X,I8))', &
            & err,error,*999)
        ELSE
          CALL FlagError("Boundary condition variable is not associated",err,error,*999)
        ENDIF
      ENDDO !variableIdx
    ENDIF

    EXITS("BoundaryConditions_CreateFinish")
    RETURN
999 ERRORSEXITS("BoundaryConditions_CreateFinish",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditions_CreateFinish

  !
  !================================================================================================================================
  !

  !>Start the creation of boundary conditions for the equation set.
  SUBROUTINE BoundaryConditions_CreateStart(solverEquations,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to create boundary conditions for
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<On exit, a pointer to the created boundary conditions. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_CreateStart",err,error,*999)

    IF(ASSOCIATED(solverEquations)) THEN
      IF(ASSOCIATED(solverEquations%boundaryConditions)) THEN
        CALL FlagError("Boundary conditions are already associated for the solver equations.",err,error,*999)
      ELSE
        IF(ASSOCIATED(boundaryConditions)) THEN
          CALL FlagError("Boundary conditions is already associated.",err,error,*999)
        ELSE
          IF(ASSOCIATED(solverEquations%solverMapping)) THEN
            !Initialise the boundary conditions
            CALL BoundaryConditions_Initialise(solverEquations,err,error,*999)
          ELSE
            localError="Solver equations solver mapping is not associated."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          !Return the pointer
          boundaryConditions=>solverEquations%boundaryConditions
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Solver equations is not associated.",err,error,*999)
    ENDIF

    EXITS("BoundaryConditions_CreateStart")
    RETURN
999 ERRORSEXITS("BoundaryConditions_CreateStart",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditions_CreateStart

  !
  !================================================================================================================================
  !

  !>Destroys boundary conditions
  SUBROUTINE BoundaryConditions_Destroy(boundaryConditions,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("BoundaryConditions_Destroy",err,error,*999)

    IF(ASSOCIATED(boundaryConditions)) THEN
      CALL BoundaryConditions_Finalise(boundaryConditions,err,error,*999)
    ELSE
      CALL FlagError("Boundary conditions is not associated.",err,error,*999)
    ENDIF

    EXITS("BoundaryConditions_Destroy")
    RETURN
999 ERRORSEXITS("BoundaryConditions_Destroy",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditions_Destroy

  !
  !================================================================================================================================
  !

  !>Finalise the boundary conditions and deallocate all memory.
  SUBROUTINE BoundaryConditions_Finalise(boundaryConditions,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: variableIdx

    ENTERS("BoundaryConditions_Finalise",err,error,*999)

    IF(ASSOCIATED(boundaryConditions)) THEN
      IF(ALLOCATED(boundaryConditions%boundaryConditionsVariables)) THEN
        DO variableIdx=1,boundaryConditions%numberOfBoundaryConditionsVariables
          IF(ASSOCIATED(boundaryConditions%boundaryConditionsVariables(variableIdx)%PTR)) THEN
            CALL BoundaryCondition_VariableFinalise(boundaryConditions%boundaryConditionsVariables(variableIdx)%PTR, &
                & err,error,*999)
          ELSE
            CALL FlagError("Boundary conditions variable number "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
                  & " is not associated",err,error,*999)
          ENDIF
        ENDDO !variableIdx
        NULLIFY(boundaryConditions%solverEquations%SOLVER%solverEquations)
        !boundaryConditions%solverEquations%SOLVER_equationsFinished = .FALSE.
        !boundaryConditions%solverEquations%solverMapping%solverMappingFinished = .FALSE.
        DEALLOCATE(boundaryConditions%boundaryConditionsVariables)
      ENDIF
      DEALLOCATE(boundaryConditions)
    ENDIF

    EXITS("BoundaryConditions_Finalise")
    RETURN
999 ERRORSEXITS("BoundaryConditions_Finalise",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the boundary conditions for an equations set.
  SUBROUTINE BoundaryConditions_Initialise(solverEquations,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to initialise the boundary conditions for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,variableIdx,variableType,equationsSetIdx,interfaceConditionIdx
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(InterfaceMappingRHSType), POINTER :: interfaceRHSMapping
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("BoundaryConditions_Initialise",err,error,*998)

    IF(ASSOCIATED(solverEquations)) THEN
      IF(ASSOCIATED(solverEquations%boundaryConditions)) THEN
        CALL FlagError("Boundary conditions is already associated for these solver equations.",err,error,*998)
      ELSE
        IF(ASSOCIATED(solverEquations%solverMapping)) THEN
          ALLOCATE(solverEquations%boundaryConditions,STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate boundary conditions.",err,error,*999)
          solverEquations%boundaryConditions%boundaryConditionsFinished=.FALSE.
          solverEquations%boundaryConditions%numberOfBoundaryConditionsVariables=0
          solverEquations%boundaryConditions%solverEquations=>solverEquations
          solverEquations%boundaryConditions%neumannMatrixSparsity=BOUNDARY_CONDITION_SPARSE_MATRICES
          DO equationsSetIdx=1,solverEquations%solverMapping%numberOfEquationsSets
            equationsSet=>solverEquations%solverMapping%equationsSets(equationsSetIdx)%PTR
            IF(ASSOCIATED(equationsSet)) THEN
              equations=>equationsSet%equations
              IF(ASSOCIATED(equations)) THEN
                NULLIFY(vectorEquations)
                CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                IF(equations%equationsFinished) THEN
                  NULLIFY(vectorMapping)
                  CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
                  IF(vectorMapping%vectorMappingFinished) THEN
                    equationsSet%boundaryConditions=>solverEquations%boundaryConditions
                    SELECT CASE(equations%timeDependence)
                    CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC)
                      SELECT CASE(equations%linearity)
                      CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
                        linearMapping=>vectorMapping%linearMapping
                        IF(ASSOCIATED(linearMapping)) THEN
                          DO variableIdx=1,linearMapping%numberOfLinearMatrixVariables
                            variableType=linearMapping%linearMatrixVariableTypes(variableIdx)
                            IF(linearMapping%varToEquationsMatricesMaps(variableType)%numberOfEquationsMatrices>0) THEN
                              CALL BoundaryConditions_VariableInitialise(solverEquations%boundaryConditions, &
                                & linearMapping%varToEquationsMatricesMaps(variableType)%VARIABLE,err,error,*999)
                            ENDIF
                          ENDDO !variableIdx
                        ELSE
                          CALL FlagError("Equations mapping linear mapping is not associated.",err,error,*999)
                        ENDIF
                        rhsMapping=>vectorMapping%rhsMapping
                        IF(ASSOCIATED(rhsMapping)) THEN
                          CALL BoundaryConditions_VariableInitialise(solverEquations%boundaryConditions, &
                            & rhsMapping%rhsVariable,err,error,*999)
                        ENDIF
                      CASE(EQUATIONS_NONLINEAR)
                        nonlinearMapping=>vectorMapping%nonlinearMapping
                        IF(ASSOCIATED(nonlinearMapping)) THEN
                          DO variableIdx=1,nonlinearMapping%numberOfResidualVariables
                            CALL BoundaryConditions_VariableInitialise(solverEquations%boundaryConditions, &
                              & nonlinearMapping%residualVariables(variableIdx)%PTR,err,error,*999)
                          ENDDO
                        ELSE
                          CALL FlagError("Equations mapping nonlinear mapping is not associated.",err,error,*999)
                        ENDIF
                        rhsMapping=>vectorMapping%rhsMapping
                        IF(ASSOCIATED(rhsMapping)) THEN
                          CALL BoundaryConditions_VariableInitialise(solverEquations%boundaryConditions, &
                            & rhsMapping%rhsVariable,err,error,*999)
                        ELSE
                          CALL FlagError("Equations mapping RHS mapping is not associated.",err,error,*999)
                        ENDIF
                      CASE DEFAULT
                        localError="The equations linearity type of "//TRIM(NumberToVString(equations%linearity,"*", &
                          & err,error))//" is invalid."
                        CALL FlagError(localError,err,error,*999)
                      END SELECT
                    CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
                      SELECT CASE(equations%linearity)
                      CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
                        dynamicMapping=>vectorMapping%dynamicMapping
                        IF(ASSOCIATED(dynamicMapping)) THEN
                          CALL BoundaryConditions_VariableInitialise(solverEquations%boundaryConditions, &
                            & dynamicMapping%dynamicVariable,err,error,*999)
                        ELSE
                          CALL FlagError("Equations mapping dynamic mapping is not associated.",err,error,*999)
                        ENDIF
                        rhsMapping=>vectorMapping%rhsMapping
                        IF(ASSOCIATED(rhsMapping)) THEN
                          CALL BoundaryConditions_VariableInitialise(solverEquations%boundaryConditions, &
                            & rhsMapping%rhsVariable,err,error,*999)
                        ELSE
                          CALL FlagError("Equations mapping RHS mapping is not associated.",err,error,*999)
                        ENDIF
                      CASE(EQUATIONS_NONLINEAR)
                        dynamicMapping=>vectorMapping%dynamicMapping
                        IF(ASSOCIATED(dynamicMapping)) THEN
                          CALL BoundaryConditions_VariableInitialise(solverEquations%boundaryConditions, &
                            & dynamicMapping%dynamicVariable,err,error,*999)
                        ELSE
                          CALL FlagError("Equations mapping dynamic mapping is not associated.",err,error,*999)
                        ENDIF
                        rhsMapping=>vectorMapping%rhsMapping
                        IF(ASSOCIATED(rhsMapping)) THEN
                          CALL BoundaryConditions_VariableInitialise(solverEquations%boundaryConditions, &
                            & rhsMapping%rhsVariable,err,error,*999)
                        ELSE
                          CALL FlagError("Equations mapping RHS mapping is not associated.",err,error,*999)
                        ENDIF
                      CASE DEFAULT
                        localError="The equations linearity type of "//TRIM(NumberToVString(equations%linearity,"*", &
                          & err,error))//" is invalid."
                        CALL FlagError(localError,err,error,*999)
                      END SELECT
                    CASE DEFAULT
                      localError="The equations time dependence type of "// &
                        & TRIM(NumberToVString(equations%timeDependence,"*",err,error))//" is invalid."
                      CALL FlagError(localError,err,error,*999)
                    END SELECT
                  ELSE
                    CALL FlagError("Equations mapping has not been finished.",err,error,*998)
                  ENDIF
                ELSE
                  CALL FlagError("Equations has not been finished.",err,error,*998)
                ENDIF
              ELSE
                CALL FlagError("Equations set equations is not associated.",err,error,*998)
              ENDIF
            ELSE
              CALL FlagError("Equations set is not associated.",err,error,*998)
            ENDIF
          ENDDO !equationsSetIdx
          DO interfaceConditionIdx=1,solverEquations%solverMapping%numberOfInterfaceConditions
            interfaceCondition=>solverEquations%solverMapping%interfaceConditions(interfaceConditionIdx)%PTR
            IF(ASSOCIATED(interfaceCondition)) THEN
              interfaceEquations=>interfaceCondition%interfaceEquations
              CALL InterfaceEquations_AssertIsFinished(interfaceEquations,err,error,*999)
              interfaceMapping=>interfaceEquations%interfaceMapping
              CALL InterfaceMapping_AssertIsFinished(interfaceMapping,err,error,*999)
              interfaceCondition%boundaryConditions=>solverEquations%boundaryConditions
              !Only linear interface equations implemented at the moment
              SELECT CASE(interfaceEquations%timeDependence)
              CASE(INTERFACE_EQUATIONS_STATIC,INTERFACE_EQUATIONS_QUASISTATIC)
                SELECT CASE(interfaceEquations%linearity)
                CASE(INTERFACE_EQUATIONS_LINEAR)
                  interfaceMapping=>interfaceEquations%interfaceMapping
                  IF(ASSOCIATED(interfaceMapping)) THEN
                    variableType=interfaceMapping%lagrangeVariableType
                    IF(interfaceMapping%numberOfInterfaceMatrices>0) THEN
                      CALL BoundaryConditions_VariableInitialise(solverEquations%boundaryConditions, &
                        & interfaceMapping%lagrangeVariable,err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Interface mapping mapping is not associated.",err,error,*999)
                  ENDIF
                  interfaceRHSMapping=>interfaceMapping%rhsMapping
                  IF(ASSOCIATED(interfaceRHSMapping)) THEN
                    CALL BoundaryConditions_VariableInitialise(solverEquations%boundaryConditions, &
                      & interfaceRHSMapping%rhsVariable,err,error,*999)
                  ENDIF
                CASE DEFAULT
                  localError="The interface equations linearity type of "// &
                    & TRIM(NumberToVString(interfaceEquations%linearity,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE DEFAULT
                localError="The interface equations time dependence type of "// &
                  & TRIM(NumberToVString(interfaceEquations%timeDependence,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ELSE
              CALL FlagError("Interface condition not associated.",err,error,*998)
            ENDIF
          ENDDO !interfaceConditionIdx
        ELSE
          CALL FlagError("Solver equations solver mapping is not associated.",err,error,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Solver equations is not associated",err,error,*998)
    ENDIF

    EXITS("BoundaryConditions_Initialise")
    RETURN
999 CALL BoundaryConditions_Finalise(solverEquations%boundaryConditions,dummyErr,dummyError,*998)
998 ERRORSEXITS("BoundaryConditions_Initialise",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditions_Initialise

  !
  !================================================================================================================================
  !

  !>Adds to the value of the specified constant and sets this as a boundary condition on the specified constant. \see OPENCMISS::CMISSBoundaryConditionAddConstant
  SUBROUTINE BoundaryConditions_AddConstant(boundaryConditions,field,variableType,componentNumber,condition,bcValue,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldType), POINTER :: field !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValue !<The value of the boundary condition to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: localDOF,globalDOF
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_AddConstant",err,error,*999)

    NULLIFY(boundaryConditionsVariable)
    NULLIFY(fieldVariable)

    !Note: This routine is for constant interpolation
    CALL BoundaryConditions_AssertNotFinished(boundaryConditions,err,error,*999)
    IF(ASSOCIATED(field)) THEN
      CALL Field_ComponentDOFGetConstant(field,variableType,componentNumber,localDOF,globalDOF, &
        & err,error,*999)
      CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
      NULLIFY(boundaryConditionsVariable)
      CALL BoundaryConditions_VariableGet(boundaryConditions,fieldVariable,boundaryConditionsVariable, &
        & err,error,*999)
      CALL BoundaryConditions_CheckInterpolationType(condition,field,variableType,componentNumber,err,error,*999)
      CALL BoundaryConditions_AddLocalDOF(boundaryConditions,field,variableType, &
        & localDOF,condition,bcValue,err,error,*999)
    ELSE
      CALL FlagError("The dependent field is not associated.",err,error,*999)
    ENDIF

    EXITS("BoundaryConditions_AddConstant")
    RETURN
999 ERRORSEXITS("BoundaryConditions_AddConstant",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_AddConstant

 !
  !================================================================================================================================
  !

  !>Sets a boundary condition on the specified constant. \see OPENCMISS::CMISSBoundaryConditionsSetConstant
  SUBROUTINE BoundaryConditions_SetConstant(boundaryConditions,field,variableType,componentNumber,condition,bcValue,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldType), POINTER :: field !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValue !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: localDOF,globalDOF
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_SetConstant",err,error,*999)

    !Note: This routine is for constant interpolation
    CALL BoundaryConditions_AssertNotFinished(boundaryConditions,err,error,*999)
    IF(ASSOCIATED(field)) THEN
      CALL Field_ComponentDOFGetConstant(field,variableType,componentNumber,localDOF,globalDOF, &
        & err,error,*999)
      CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
      NULLIFY(boundaryConditionsVariable)
      CALL BoundaryConditions_VariableGet(boundaryConditions,fieldVariable,boundaryConditionsVariable,err,error,*999)
      CALL BoundaryConditions_CheckInterpolationType(condition,field,variableType,componentNumber,err,error,*999)
      CALL BoundaryConditions_SetLocalDOF(boundaryConditions,field,variableType, &
        & localDOF,condition,bcValue,err,error,*999)
    ELSE
      CALL FlagError("The dependent field is not associated.",err,error,*999)
    ENDIF
    
    EXITS("BoundaryConditions_SetConstant")
    RETURN
999 ERRORSEXITS("BoundaryConditions_SetConstant",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditions_SetConstant

  !
  !================================================================================================================================
  !

  !>Adds to the value of the specified DOF and sets this as a boundary condition on the specified DOF.
  SUBROUTINE BoundaryConditions_AddLocalDOF0(boundaryConditions,field,variableType,dofIndex,condition,bcValue,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldType), POINTER :: field !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: dofIndex !<The local dof index to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValue !<The value of the boundary condition to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("BoundaryConditions_AddLocalDOF0",err,error,*999)

    CALL BoundaryConditions_AddLocalDOF1(boundaryConditions,field,variableType,[dofIndex],[condition],[bcValue], &
        & err,error,*999)

    EXITS("BoundaryConditions_AddLocalDOF0")
    RETURN
999 ERRORSEXITS("BoundaryConditions_AddLocalDOF0",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_AddLocalDOF0

  !
  !================================================================================================================================
  !

  !>Adds to the value of the specified DOF and sets this as a boundary condition on the specified DOFs.
  SUBROUTINE BoundaryConditions_AddLocalDOF1(boundaryConditions,field,variableType,dofIndices,conditions,bcValues,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldType), POINTER :: field !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: dofIndices(:) !<dofIndices(:). The local dof index for the i'th dof to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: conditions(:) !<conditions(:). The boundary condition type to set for the i'th dof \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValues(:) !<bcValues(:). The value of the boundary condition for the i'th dof to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i,globalDOF,localDOF
    REAL(DP) :: initialValue
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable
    TYPE(DomainMappingType), POINTER :: domainMapping
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_AddLocalDOF1",err,error,*999)
    NULLIFY(dependent_variable)

    CALL BoundaryConditions_AssertNotFinished(boundaryConditions,err,error,*999)
    IF(ASSOCIATED(field)) THEN
      NULLIFY(fieldVariable)
      CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
      IF(ASSOCIATED(fieldVariable)) THEN
        domainMapping=>fieldVariable%domainMapping
        IF(ASSOCIATED(domainMapping)) THEN
          NULLIFY(boundaryConditionsVariable)
          CALL BoundaryConditions_VariableGet(boundaryConditions,fieldVariable,boundaryConditionsVariable, &
            & err,error,*999)
          IF(SIZE(dofIndices,1)==SIZE(conditions,1)) THEN
            IF(SIZE(dofIndices,1)==SIZE(bcValues,1)) THEN
              DO i=1,SIZE(dofIndices,1)
                localDOF=dofIndices(i)
                IF(localDOF>=1.AND.localDOF<=domainMapping%numberOfLocal) THEN
                  globalDOF=domainMapping%localToGlobalMap(localDOF)
                  ! Set boundary condition and dof type, and make sure parameter sets are created
                  CALL BoundaryConditions_SetConditionType(boundaryConditionsVariable,globalDOF,conditions(i), &
                    & err,error,*999)
                  ! Update field sets by adding boundary condition values
                  SELECT CASE(conditions(i))
                  CASE(BOUNDARY_CONDITION_FREE)
                    ! No field update
                  CASE(BOUNDARY_CONDITION_FIXED)
                    CALL Field_ParameterSetAddLocalDOF(field,variableType,FIELD_VALUES_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_FIXED_INLET)
                    CALL Field_ParameterSetAddLocalDOF(field,variableType,FIELD_VALUES_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_FIXED_OUTLET)
                    CALL Field_ParameterSetAddLocalDOF(field,variableType,FIELD_VALUES_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_FIXED_WALL)
                    CALL Field_ParameterSetAddLocalDOF(field,variableType,FIELD_VALUES_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_MOVED_WALL)
                    CALL Field_ParameterSetAddLocalDOF(field,variableType,FIELD_VALUES_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_FREE_WALL)
                    CALL Field_ParameterSetAddLocalDOF(field,variableType,FIELD_VALUES_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED)
                    CALL Field_ParameterSetAddLocalDOF(field,variableType,FIELD_VALUES_SET_TYPE,localDOF,bcValues(i), &
                      & err,error,*999)
                    CALL Field_ParameterSetAddLocalDOF(field,variableType,FIELD_BOUNDARY_conditions_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_FIXED_INCREMENTED)
                    ! For increment loops, we need to set the full BC parameter set value by
                    ! getting the current value from the values parameter set
                    CALL Field_ParameterSetGetLocalDOF(field,variableType,FIELD_VALUES_SET_TYPE, &
                      & localDOF,initialValue,err,error,*999)
                    CALL Field_ParameterSetUpdateLocalDOF(field,variableType,FIELD_BOUNDARY_conditions_SET_TYPE, &
                      & localDOF,initialValue+bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_PRESSURE)
                    CALL Field_ParameterSetAddLocalDOF(field,variableType,FIELD_PRESSURE_VALUES_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
                    ! For pressure incremented, adding to the values_set parameter value doesn't make sense,
                    ! so just increment the value in the pressure values parameter set
                    CALL Field_ParameterSetAddLocalDOF(field,variableType,FIELD_PRESSURE_VALUES_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE)
                    ! No field update
                  CASE(BOUNDARY_CONDITION_IMPERMEABLE_WALL)
                    CALL Field_ParameterSetAddLocalDOF(field,variableType,FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_NEUMANN_POINT)
                    ! Point value is stored in boundary conditions field set, and is then integrated to
                    ! get the RHS variable value
                    CALL Field_ParameterSetAddLocalDOF(field,variableType,FIELD_BOUNDARY_conditions_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
                    CALL Field_ParameterSetAddLocalDOF(field,variableType,FIELD_BOUNDARY_conditions_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_NEUMANN_INTEGRATED,BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY)
                    ! For integrated Neumann condition, integration is already done, so set the RHS
                    ! dof value directly
                    CALL Field_ParameterSetAddLocalDOF(field,variableType,FIELD_VALUES_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_FIXED_FITTED,BOUNDARY_CONDITION_FIXED_NONREFLECTING, &
                    &  BOUNDARY_CONDITION_FIXED_CELLML,BOUNDARY_CONDITION_FIXED_STREE)
                    CALL Field_ParameterSetAddLocalDOF(field,variableType,FIELD_VALUES_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE DEFAULT
                    localError="The specified boundary condition type for dof index "// &
                      & TRIM(NumberToVString(i,"*",err,error))//" of "// &
                      & TRIM(NumberToVString(conditions(i),"*",err,error))//" is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                ELSE
                  localError="The local dof of  "//&
                    & TRIM(NumberToVString(localDOF,"*",err,error))//" at dof index "// &
                    & TRIM(NumberToVString(i,"*",err,error))// &
                    & " is invalid. The dof should be between 1 and "// &
                    & TRIM(NumberToVString(domainMapping%numberOfLocal,"*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
              ENDDO !i
            ELSE
              localError="The size of the dof indices array ("// &
                & TRIM(NumberToVString(SIZE(dofIndices,1),"*",err,error))// &
                & ") does not match the size of the values array ("// &
                & TRIM(NumberToVString(SIZE(bcValues,1),"*",err,error))//")."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="The size of the dof indices array ("// &
              & TRIM(NumberToVString(SIZE(dofIndices,1),"*",err,error))// &
              & ") does not match the size of the fixed conditions array ("// &
              & TRIM(NumberToVString(SIZE(conditions,1),"*",err,error))//")."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("The dependent field variable domain mapping is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("The dependent field variable is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("The dependent field is not associated..",err,error,*999)
    ENDIF

    EXITS("BoundaryConditions_AddLocalDOF1")
    RETURN
999 ERRORSEXITS("BoundaryConditions_AddLocalDOF1",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_AddLocalDOF1

  !
  !================================================================================================================================
  !

  !>Sets a boundary condition on the specified DOF.
  SUBROUTINE BoundaryConditions_SetLocalDOF0(boundaryConditions,field,variableType,dofIndex,condition,bcValue,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldType), POINTER :: field !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: dofIndex !<The local dof index to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValue !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("BoundaryConditions_SetLocalDOF0",err,error,*999)

    CALL BoundaryConditions_SetLocalDOF1(boundaryConditions,field,variableType,[dofIndex],[condition],[bcValue], &
      & err,error,*999)

    EXITS("BoundaryConditions_SetLocalDOF0")
    RETURN
999 ERRORSEXITS("BoundaryConditions_SetLocalDOF0",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_SetLocalDOF0

  !
  !================================================================================================================================
  !

  !>Sets a boundary condition on the specified DOFs.
  SUBROUTINE BoundaryConditions_SetLocalDOF1(boundaryConditions,field,variableType,dofIndices,conditions,bcValues,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldType), POINTER :: field !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: dofIndices(:) !<dofIndices(:). The local dof index for the i'th dof to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: conditions(:) !<conditions(:). The boundary condition type to set for the i'th dof \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValues(:) !<bcValues(:). The value of the boundary condition for the i'th dof to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i,globalDOF,localDOF
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable
    TYPE(DomainMappingType), POINTER :: domainMapping
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_SetLocalDOF1",err,error,*999)

    CALL BoundaryConditions_AssertNotFinished(boundaryConditions,err,error,*999)
    
    IF(ASSOCIATED(field)) THEN
      NULLIFY(fieldVariable)
      CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
      IF(ASSOCIATED(fieldVariable)) THEN
        domainMapping=>fieldVariable%domainMapping
        IF(ASSOCIATED(domainMapping)) THEN
          NULLIFY(boundaryConditionsVariable)
          CALL BoundaryConditions_VariableGet(boundaryConditions,fieldVariable,boundaryConditionsVariable, &
            & err,error,*999)
          IF(SIZE(dofIndices,1)==SIZE(conditions,1)) THEN
            IF(SIZE(dofIndices,1)==SIZE(bcValues,1)) THEN
              DO i=1,SIZE(dofIndices,1)
                localDOF=dofIndices(i)
                IF(localDOF>=1.AND.localDOF<=domainMapping%numberOfLocal) THEN
                  globalDOF=domainMapping%localToGlobalMap(localDOF)
                  ! Set boundary condition and dof type
                  CALL BoundaryConditions_SetConditionType(boundaryConditionsVariable,globalDOF,conditions(i), &
                    & err,error,*999)
                  ! Update field sets with boundary condition value
                  
                  SELECT CASE(conditions(i))
                  CASE(BOUNDARY_CONDITION_FREE)
                    ! No field update
                  CASE(BOUNDARY_CONDITION_FIXED)
                    CALL Field_ParameterSetUpdateLocalDOF(field,variableType,FIELD_VALUES_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_FIXED_INLET)
                    CALL Field_ParameterSetUpdateLocalDOF(field,variableType,FIELD_VALUES_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_FIXED_OUTLET)
                    CALL Field_ParameterSetUpdateLocalDOF(field,variableType,FIELD_VALUES_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_FIXED_WALL)
                    CALL Field_ParameterSetUpdateLocalDOF(field,variableType,FIELD_VALUES_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_MOVED_WALL)
                    CALL Field_ParameterSetUpdateLocalDOF(field,variableType,FIELD_VALUES_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_FREE_WALL)
                    CALL Field_ParameterSetUpdateLocalDOF(field,variableType,FIELD_VALUES_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED)
                    CALL Field_ParameterSetUpdateLocalDOF(field,variableType,FIELD_VALUES_SET_TYPE,localDOF,bcValues(i), &
                      & err,error,*999)
                    CALL Field_ParameterSetUpdateLocalDOF(field,variableType,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_FIXED_INCREMENTED) !For load increment loops
                    CALL Field_ParameterSetUpdateLocalDOF(field,variableType,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_PRESSURE)
                    CALL Field_ParameterSetUpdateLocalDOF(field,variableType,FIELD_PRESSURE_VALUES_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_FIXED_PRESSURE)
                    CALL Field_ParameterSetUpdateLocalDOF(field,variableType,FIELD_VALUES_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
                    CALL Field_ParameterSetUpdateLocalDOF(field,variableType,FIELD_PRESSURE_VALUES_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE)
                    ! No field update
                  CASE(BOUNDARY_CONDITION_NEUMANN_POINT)
                    CALL Field_ParameterSetUpdateLocalDOF(field,variableType,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
                    CALL Field_ParameterSetUpdateLocalDOF(field,variableType,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_NEUMANN_INTEGRATED,BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY)
                    CALL Field_ParameterSetUpdateLocalDOF(field,variableType,FIELD_VALUES_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_IMPERMEABLE_WALL)
                    CALL Field_ParameterSetUpdateLocalDOF(field,variableType,FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE(BOUNDARY_CONDITION_FIXED_FITTED,BOUNDARY_CONDITION_FIXED_NONREFLECTING, &
                    & BOUNDARY_CONDITION_FIXED_CELLML,BOUNDARY_CONDITION_FIXED_STREE)
                    CALL Field_ParameterSetUpdateLocalDOF(field,variableType,FIELD_VALUES_SET_TYPE, &
                      & localDOF,bcValues(i),err,error,*999)
                  CASE DEFAULT
                    localError="The specified boundary condition type for dof index "// &
                      & TRIM(NumberToVString(i,"*",err,error))//" of "// &
                      & TRIM(NumberToVString(conditions(i),"*",err,error))//" is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                ELSE
                  localError="The local dof of  "//&
                    & TRIM(NumberToVString(localDOF,"*",err,error))//" at dof index "// &
                    & TRIM(NumberToVString(i,"*",err,error))// &
                    & " is invalid. The dof should be between 1 and "// &
                    & TRIM(NumberToVString(domainMapping%numberOfLocal,"*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
              ENDDO !i
            ELSE
              localError="The size of the dof indices array ("// &
                & TRIM(NumberToVString(SIZE(dofIndices,1),"*",err,error))// &
                & ") does not match the size of the values array ("// &
                & TRIM(NumberToVString(SIZE(bcValues,1),"*",err,error))//")."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="The size of the dof indices array ("// &
              & TRIM(NumberToVString(SIZE(dofIndices,1),"*",err,error))// &
              & ") does not match the size of the fixed conditions array ("// &
              & TRIM(NumberToVString(SIZE(conditions,1),"*",err,error))//")."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("The dependent field variable domain mapping is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("The dependent field variable is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("The dependent field is not associated.",err,error,*999)
    ENDIF

    EXITS("BoundaryConditions_SetLocalDOF1")
    RETURN
999 ERRORSEXITS("BoundaryConditions_SetLocalDOF1",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_SetLocalDOF1

  !
  !================================================================================================================================
  !

  !> Checks the boundary condition type and sets the boundary condition type and dof type for the boundary conditions.
  !> Makes sure any field parameter sets required are created, and sets the parameter set required array value.
  SUBROUTINE BoundaryConditions_SetConditionType(boundaryConditionsVariable,globalDof,condition,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to set the boundary condition for
    INTEGER(INTG), INTENT(IN) :: globalDof !<The globalDof to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: dofType, previousCondition, previousDof

    ENTERS("BoundaryConditions_SetConditionType",err,error,*999)

    ! We won't do much checking here as this is only used internally and everything has been checked for
    ! association already
    ! Don't need to make sure field values set type is available as this will always be there, but need
    ! to make sure any other parameter sets required are.
    SELECT CASE(condition)
    CASE(BOUNDARY_CONDITION_FREE)
      dofType=BOUNDARY_CONDITION_DOF_FREE
    CASE(BOUNDARY_CONDITION_FIXED)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_LINEAR_CONSTRAINT)
      dofType=BOUNDARY_CONDITION_DOF_CONSTRAINED
    CASE(BOUNDARY_CONDITION_FIXED_INLET)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_FIXED_OUTLET)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_FIXED_WALL)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_MOVED_WALL)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_FREE_WALL)
      dofType=BOUNDARY_CONDITION_DOF_FREE
    CASE(BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
      CALL Field_ParameterSetEnsureCreated(boundaryConditionsVariable%VARIABLE%FIELD,boundaryConditionsVariable%variableType, &
        & FIELD_BOUNDARY_CONDITIONS_SET_TYPE,err,error,*999)
      boundaryConditionsVariable%parameterSetRequired(FIELD_BOUNDARY_CONDITIONS_SET_TYPE)=.TRUE.
    CASE(BOUNDARY_CONDITION_FIXED_INCREMENTED) !For load increment loops
      dofType=BOUNDARY_CONDITION_DOF_FIXED
      CALL Field_ParameterSetEnsureCreated(boundaryConditionsVariable%VARIABLE%FIELD,boundaryConditionsVariable%variableType, &
        & FIELD_BOUNDARY_CONDITIONS_SET_TYPE,err,error,*999)
      boundaryConditionsVariable%parameterSetRequired(FIELD_BOUNDARY_CONDITIONS_SET_TYPE)=.TRUE.
    CASE(BOUNDARY_CONDITION_PRESSURE)
      ! Pressure boundary conditions leave the RHS dof as free, as the Neumann terms
      ! are calculated in finite elasticity routines when calculating the element residual
      dofType=BOUNDARY_CONDITION_DOF_FREE
      CALL Field_ParameterSetEnsureCreated(boundaryConditionsVariable%VARIABLE%FIELD,boundaryConditionsVariable%variableType, &
        & FIELD_PRESSURE_VALUES_SET_TYPE,err,error,*999)
      boundaryConditionsVariable%parameterSetRequired(FIELD_PRESSURE_VALUES_SET_TYPE)=.TRUE.
    CASE(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
      dofType=BOUNDARY_CONDITION_DOF_FREE
      CALL Field_ParameterSetEnsureCreated(boundaryConditionsVariable%VARIABLE%FIELD,boundaryConditionsVariable%variableType, &
        & FIELD_PRESSURE_VALUES_SET_TYPE,err,error,*999)
      boundaryConditionsVariable%parameterSetRequired(FIELD_PRESSURE_VALUES_SET_TYPE)=.TRUE.
      CALL Field_ParameterSetEnsureCreated(boundaryConditionsVariable%VARIABLE%FIELD,boundaryConditionsVariable%variableType, &
        & FIELD_PREVIOUS_PRESSURE_SET_TYPE,err,error,*999)
      boundaryConditionsVariable%parameterSetRequired(FIELD_PREVIOUS_PRESSURE_SET_TYPE)=.TRUE.
    CASE(BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_IMPERMEABLE_WALL)
      dofType=BOUNDARY_CONDITION_DOF_FREE
      CALL Field_ParameterSetEnsureCreated(boundaryConditionsVariable%VARIABLE%FIELD,boundaryConditionsVariable%variableType, &
        & FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE,err,error,*999)
      boundaryConditionsVariable%parameterSetRequired(FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE)=.TRUE.
    CASE(BOUNDARY_CONDITION_NEUMANN_POINT,BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
      CALL Field_ParameterSetEnsureCreated(boundaryConditionsVariable%VARIABLE%FIELD,boundaryConditionsVariable%variableType, &
        & FIELD_BOUNDARY_CONDITIONS_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetEnsureCreated(boundaryConditionsVariable%VARIABLE%FIELD,boundaryConditionsVariable%variableType, &
        & FIELD_INTEGRATED_NEUMANN_SET_TYPE,err,error,*999)
      boundaryConditionsVariable%parameterSetRequired(FIELD_BOUNDARY_CONDITIONS_SET_TYPE)=.TRUE.
      boundaryConditionsVariable%parameterSetRequired(FIELD_INTEGRATED_NEUMANN_SET_TYPE)=.TRUE.
    CASE(BOUNDARY_CONDITION_NEUMANN_INTEGRATED,BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_FIXED_FITTED,BOUNDARY_CONDITION_FIXED_NONREFLECTING,BOUNDARY_CONDITION_FIXED_CELLML, &
      & BOUNDARY_CONDITION_FIXED_STREE)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE DEFAULT
      CALL FlagError("The specified boundary condition type for dof number "// &
        & TRIM(NumberToVString(globalDof,"*",err,error))//" of "// &
        & TRIM(NumberToVString(condition,"*",err,error))//" is invalid.", &
        & err,error,*999)
    END SELECT

    !We have a valid boundary condition type
    !Update condition type counts
    previousCondition=boundaryConditionsVariable%conditionTypes(globalDof)
    IF(previousCondition/=condition) THEN
      ! dofCounts array doesn't include a count for BOUNDARY_CONDITION_FREE, which equals zero
      IF(previousCondition/=BOUNDARY_CONDITION_FREE) THEN
        boundaryConditionsVariable%dofCounts(previousCondition)= &
          & boundaryConditionsVariable%dofCounts(previousCondition)-1
      END IF
      IF(condition/=BOUNDARY_CONDITION_FREE) THEN
        boundaryConditionsVariable%dofCounts(condition)= &
          & boundaryConditionsVariable%dofCounts(condition)+1
      END IF
    END IF
    !Update Dirichlet DOF count
    previousDof=boundaryConditionsVariable%DOFTypes(globalDof)
    IF(dofType==BOUNDARY_CONDITION_DOF_FIXED.AND.previousDof/=BOUNDARY_CONDITION_DOF_FIXED) THEN
      boundaryConditionsVariable%numberOfDirichletConditions= &
        & boundaryConditionsVariable%numberOfDirichletConditions+1
    ELSE IF(dofType/=BOUNDARY_CONDITION_DOF_FIXED.AND.previousDof==BOUNDARY_CONDITION_DOF_FIXED) THEN
      boundaryConditionsVariable%numberOfDirichletConditions= &
        & boundaryConditionsVariable%numberOfDirichletConditions-1
    END IF

    !Set the boundary condition type and DOF type
    boundaryConditionsVariable%conditionTypes(globalDof)=condition
    boundaryConditionsVariable%DOFTypes(globalDof)=dofType
    IF(DIAGNOSTICS1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Boundary Condition Being Set",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"global dof = ", globalDof,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Variable Type = ", &
        & boundaryConditionsVariable%variableType,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"New Condition = ", &
        & condition,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"dof type = ", &
        & dofType,err,error,*999)
    ENDIF
    EXITS("BoundaryConditions_SetConditionType")
    RETURN
999 ERRORSEXITS("BoundaryConditions_SetConditionType",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_SetConditionType

  !
  !================================================================================================================================
  !

  !>Adds to the value of the specified constant and sets this as a boundary condition on the specified user element. \see OPENCMISS_CMISSBoundaryConditionsAddElement
  SUBROUTINE BoundaryConditions_AddElement(boundaryConditions,field,variableType,userElementNumber,componentNumber, &
    & condition,VALUE,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldType), POINTER :: field !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValue !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: localDOF,globalDOF
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_AddElement",err,error,*999)

    !Note: this routine is for element based interpolation
    CALL BoundaryConditions_AssertNotFinished(boundaryConditions,err,error,*999)
    IF(ASSOCIATED(field)) THEN
      CALL Field_ComponentDOFGetUserElement(field,variableType,userElementNumber,componentNumber, &
        & localDOF,globalDOF,err,error,*999)
      NULLIFY(fieldVariable)
      CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
      IF(ASSOCIATED(fieldVariable)) THEN
        NULLIFY(boundaryConditionsVariable)        
        CALL BoundaryConditions_VariableGet(boundaryConditions,fieldVariable,boundaryConditionsVariable,err,error,*999)
        CALL BoundaryConditions_CheckInterpolationType(condition,field,variableType,componentNumber,err,error,*999)
        CALL BoundaryConditions_AddLocalDOF(boundaryConditions,field,variableType, &
          & localDOF,condition,bcValue,err,error,*999)
      ELSE
        CALL FlagError("The dependent field variable is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("The dependent field is not associated.",err,error,*999)
    ENDIF

    EXITS("BoundaryConditions_AddElement")
    RETURN
999 ERRORSEXITS("BoundaryConditions_AddElement",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_AddElement

  !
  !================================================================================================================================
  !

  !> Checks that the specified boundary condition is appropriate for the field variable interpolation type
  SUBROUTINE BoundaryConditions_CheckInterpolationType(condition,field,variableType,componentNumber,err,error,*)

    ! Argument variables
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type being set
    TYPE(FieldType), POINTER :: field !<A pointer to the field the boundary condition is set on
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type the boundary condition is set on
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number the boundary condition is set on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    ! Local variables
    INTEGER(INTG) :: interpolationType
    LOGICAL :: validCondition

    ENTERS("BoundaryConditions_CheckInterpolationType",err,error,*999)

    CALL Field_ComponentInterpolationGet(field,variableType,componentNumber,interpolationType,err,error,*999)

    validCondition=.TRUE.
    SELECT CASE(condition)
    CASE(BOUNDARY_CONDITION_FREE, &
        & BOUNDARY_CONDITION_FIXED, &
        & BOUNDARY_CONDITION_FIXED_INCREMENTED)
      ! Valid for all interpolation types
    CASE(BOUNDARY_CONDITION_FIXED_INLET, &
        & BOUNDARY_CONDITION_FIXED_OUTLET)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
        validCondition=.FALSE.
      END IF
    CASE(BOUNDARY_CONDITION_FIXED_WALL, &
        & BOUNDARY_CONDITION_MOVED_WALL, &
        & BOUNDARY_CONDITION_FREE_WALL, &
        & BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
        validCondition=.FALSE.
      END IF
    CASE(BOUNDARY_CONDITION_PRESSURE, &
        & BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
        validCondition=.FALSE.
      END IF
    CASE(BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
        validCondition=.FALSE.
      END IF
    CASE(BOUNDARY_CONDITION_IMPERMEABLE_WALL)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
        validCondition=.FALSE.
      END IF
    CASE(BOUNDARY_CONDITION_NEUMANN_POINT, &
        & BOUNDARY_CONDITION_NEUMANN_INTEGRATED, &
        & BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY, &
        & BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
        validCondition=.FALSE.
      END IF
    CASE(BOUNDARY_CONDITION_FIXED_FITTED,BOUNDARY_CONDITION_FIXED_NONREFLECTING,BOUNDARY_CONDITION_FIXED_CELLML, &
      & BOUNDARY_CONDITION_FIXED_STREE)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
        validCondition=.FALSE.
      END IF
    CASE DEFAULT
      CALL FlagError("The specified boundary condition type of "// &
        & TRIM(NumberToVString(condition,"*",err,error))//" is invalid.", &
        & err,error,*999)
    END SELECT
    IF(.NOT.validCondition) THEN
      CALL FlagError("The specified boundary condition type of "// &
        & TRIM(NumberToVString(condition,"*",err,error))//" is not valid for the field component "// &
        & "interpolation type of "//TRIM(NumberToVString(interpolationType,"*",err,error))//".", &
        & err,error,*999)
    END IF

    EXITS("BoundaryConditions_CheckInterpolationType")
    RETURN
999 ERRORSEXITS("BoundaryConditions_CheckInterpolationType",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_CheckInterpolationType

  !
  !================================================================================================================================
  !

  !> Checks that the applied boundary conditions are supported by the equations sets in the solver equations
  SUBROUTINE BoundaryConditions_CheckEquations(boundaryConditionsVariable,err,error,*)

    ! Argument variables
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to check
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    type(varying_string), intent(out) :: error !<The error string
    ! Local variables
    INTEGER(INTG) :: boundaryConditionType,equationsSetIdx,specificationSize
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(EquationsSetType), POINTER :: equationsSet
    LOGICAL :: validEquationsSetFound

    ENTERS("BoundaryConditions_CheckEquations",err,error,*999)

    !Get and check pointers we need
    solverEquations=>boundaryConditionsVariable%boundaryConditions%solverEquations
    IF(.NOT.ASSOCIATED(solverEquations)) THEN
      CALL FlagError("Boundary conditions solver equations are not associated.",err,error,*999)
    END IF
    solverMapping=>solverEquations%solverMapping
    IF(.NOT.ASSOCIATED(solverMapping)) THEN
      CALL FlagError("Solver equations solver mapping is not associated.",err,error,*999)
    END IF

    DO boundaryConditionType=1,MAX_BOUNDARY_CONDITION_NUMBER
      !Check if any DOFs have been set to this BC type
      IF(boundaryConditionsVariable%dofCounts(boundaryConditionType)>0) THEN
        validEquationsSetFound=.FALSE.
        DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
          equationsSet=>solverMapping%equationsSets(equationsSetIdx)%PTR
          IF(.NOT.ASSOCIATED(equationsSet)) THEN
            CALL FlagError("Solver equations equations set is not associated.",err,error,*999)
          END IF
          IF(.NOT.ALLOCATED(equationsSet%specification)) THEN
            CALL FlagError("Equations set specification is not allocated.",err,error,*999)
          END IF
          specificationSize=SIZE(equationsSet%specification,1)

          SELECT CASE(boundaryConditionType)
          CASE(BOUNDARY_CONDITION_FREE)
            ! Valid for any equations set
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_FIXED)
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_LINEAR_CONSTRAINT)
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_FIXED_INLET, &
              & BOUNDARY_CONDITION_FIXED_OUTLET)
            IF(specificationSize>=2) THEN
              IF(equationsSet%specification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS.AND. &
                  & (equationsSet%specification(2)==EQUATIONS_SET_STOKES_EQUATION_TYPE.OR. &
                  & equationsSet%specification(2)==EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE.OR. &
                  & equationsSet%specification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)) THEN
                validEquationsSetFound=.TRUE.
              END IF
            END IF
          CASE(BOUNDARY_CONDITION_FIXED_WALL,BOUNDARY_CONDITION_MOVED_WALL, &
              & BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED,BOUNDARY_CONDITION_FREE_WALL)
            IF(specificationSize>=2) THEN
              IF(equationsSet%specification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS.AND. &
                  & (equationsSet%specification(2)==EQUATIONS_SET_STOKES_EQUATION_TYPE.OR. &
                  & equationsSet%specification(2)==EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE.OR. &
                  & equationsSet%specification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE.OR. &
                  & equationsSet%specification(2)==EQUATIONS_SET_DARCY_EQUATION_TYPE)) THEN
                validEquationsSetFound=.TRUE.
              ELSE IF(specificationSize==3) THEN
                IF(equationsSet%specification(1)==EQUATIONS_SET_CLASSICAL_FIELD_CLASS.AND. &
                    & equationsSet%specification(2)==EQUATIONS_SET_LAPLACE_EQUATION_TYPE.AND. &
                    & equationsSet%specification(3)==EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE) THEN
                  validEquationsSetFound=.TRUE.
                END IF
              END IF
            END IF
          CASE(BOUNDARY_CONDITION_FIXED_INCREMENTED)
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_PRESSURE, &
              & BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
            IF(specificationSize>=2) THEN
              IF(equationsSet%specification(1)==EQUATIONS_SET_ELASTICITY_CLASS.AND. &
                & equationsSet%specification(2)==EQUATIONS_SET_FINITE_ELASTICITY_TYPE) THEN
                validEquationsSetFound=.TRUE.
              ELSE IF(equationsSet%specification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS .AND. &
                & (equationsSet%specification(2)==EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE.OR. &
                & equationsSet%specification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)) THEN
                validEquationsSetFound=.TRUE.
              END IF
            ENDIF
          CASE(BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE)
            !Not actually used anywhere? So keep it as invalid, although maybe it should be removed?
            validEquationsSetFound=.FALSE.
          CASE(BOUNDARY_CONDITION_IMPERMEABLE_WALL)
            IF(specificationSize>=3) THEN
              IF(equationsSet%specification(1)==EQUATIONS_SET_ELASTICITY_CLASS.AND. &
                  & equationsSet%specification(2)==EQUATIONS_SET_FINITE_ELASTICITY_TYPE.AND. &
                  & (equationsSet%specification(3)==EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE.OR. &
                  & equationsSet%specification(3)==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE.OR. &
                  & equationsSet%specification(3)==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)) THEN
                validEquationsSetFound=.TRUE.
              END IF
            END IF
          CASE(BOUNDARY_CONDITION_NEUMANN_POINT,BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_NEUMANN_INTEGRATED,BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY)
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_FIXED_FITTED)
            IF(equationsSet%specification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS.AND. &
              & (equationsSet%specification(2)==EQUATIONS_SET_STOKES_EQUATION_TYPE.OR. &
              & equationsSet%specification(2)==EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE.OR. &
              & equationsSet%specification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)) THEN
              validEquationsSetFound=.TRUE.
            END IF
          CASE(BOUNDARY_CONDITION_FIXED_NONREFLECTING,BOUNDARY_CONDITION_FIXED_CELLML,BOUNDARY_CONDITION_FIXED_STREE)
            IF(equationsSet%specification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS.AND. &
                & (equationsSet%specification(2)==EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE.OR. &
                & equationsSet%specification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)) THEN
              validEquationsSetFound=.TRUE.
            END IF
          CASE DEFAULT
            CALL FlagError("The specified boundary condition type of "// &
              & TRIM(NumberToVString(boundaryConditionType,"*",err,error))// &
              & " is invalid.",err,error,*999)
          END SELECT
        END DO

        IF(.NOT.validEquationsSetFound) THEN
            CALL FlagError("The specified boundary condition type of "// &
              & TRIM(NumberToVString(boundaryConditionType,"*",err,error))// &
              & " is invalid for the equations sets in the solver equations.",err,error,*999)
        END IF
      END IF
    END DO

    EXITS("BoundaryConditions_CheckEquations")
    RETURN
999 ERRORSEXITS("BoundaryConditions_CheckEquations",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_CheckEquations

  !
  !================================================================================================================================
  !

  !>Sets a boundary condition on the specified user element. \see OPENCMISS_CMISSBoundaryConditionsSetElement
  SUBROUTINE BoundaryConditions_SetElement(boundaryConditions,field,variableType,userElementNumber,componentNumber, &
    & condition,bcValue,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldType), POINTER :: field !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValue !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: localDOF,globalDOF
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_SetElement",err,error,*999)

    !Note: this routine is for element based interpolation
    CALL BoundaryConditions_AssertNotFinished(boundaryConditions,err,error,*999)
    IF(ASSOCIATED(field)) THEN
      CALL Field_ComponentDOFGetUserElement(field,variableType,userElementNumber,componentNumber, &
        & localDOF,globalDOF,err,error,*999)
      NULLIFY(fieldVariable)
      CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
      IF(ASSOCIATED(fieldVariable)) THEN
        NULLIFY(boundaryConditionsVariable)
        CALL BoundaryConditions_VariableGet(boundaryConditions,fieldVariable,boundaryConditionsVariable,err,error,*999)
        CALL BoundaryConditions_CheckInterpolationType(condition,field,variableType,componentNumber,err,error,*999)
        CALL BoundaryConditions_SetLocalDOF(boundaryConditions,field,variableType, &
          & localDOF,condition,bcValue,err,error,*999)
      ELSE
        CALL FlagError("The dependent field variable is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("The dependent field is not associated.",err,error,*999)
    ENDIF
 
    EXITS("BoundaryConditions_SetElement")
    RETURN
999 ERRORSEXITS("BoundaryConditions_SetElement",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_SetElement

  !
  !================================================================================================================================
  !

  !>Adds to the value of the specified constant and sets this as a boundary condition on the specified user node. \see OPENCMISS_CMISSBoundaryConditionsAddNode
  SUBROUTINE BoundaryConditions_AddNode(boundaryConditions,field,variableType,versionNumber,derivativeNumber, &
    & userNodeNumber,componentNumber,condition,bcValue,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldType), POINTER :: field !<A pointer to the field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: versionNumber !<The derivative version to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: derivativeNumber !<The derivative to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValue !<The value of the boundary condition to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: localDOF,globalDOF
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_AddNode",err,error,*999)

    NULLIFY(fieldVariable)

    CALL BoundaryConditions_AssertNotFinished(boundaryConditions,err,error,*999)
    IF(ASSOCIATED(field)) THEN
      CALL Field_ComponentDOFGetUserNode(field,variableType,versionNumber,derivativeNumber, &
        & userNodeNumber,componentNumber,localDOF,globalDOF,err,error,*999)
      CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
      IF(ASSOCIATED(fieldVariable)) THEN
        NULLIFY(boundaryConditionsVariable)
        CALL BoundaryConditions_VariableGet(boundaryConditions,fieldVariable,boundaryConditionsVariable,err,error,*999)
        CALL BoundaryConditions_CheckInterpolationType(condition,field,variableType,componentNumber,err,error,*999)
        CALL BoundaryConditions_AddLocalDOF(boundaryConditions,field,variableType, &
          & localDOF,condition,bcValue,err,error,*999)
      ELSE
        CALL FlagError("The dependent field variable is not associated",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("The dependent field is not associated.",err,error,*999)
    ENDIF
    
    EXITS("BoundaryConditions_AddNode")
    RETURN
999 ERRORSEXITS("BoundaryConditions_AddNode",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_AddNode

  !
  !================================================================================================================================
  !

  !>Initialise the Neumann boundary conditions information
  SUBROUTINE BoundaryConditions_NeumannInitialise(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to initialise Neumann conditions for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(BoundaryConditionsNeumannType), POINTER :: boundaryConditionsNeumann
    INTEGER(INTG) :: numberOfValues,numberOfLocalDofs
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("BoundaryConditions_NeumannInitialise",err,error,*998)

    IF(ASSOCIATED(boundaryConditionsVariable)) THEN
      numberOfValues=boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT)+ &
        & boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
      ALLOCATE(boundaryConditionsVariable%neumannBoundaryConditions,stat=err)
      IF(err/=0) CALL FlagError("Could not allocate Neumann Boundary Conditions",err,error,*998)
      boundaryConditionsNeumann=>boundaryConditionsVariable%neumannBoundaryConditions
      IF(ASSOCIATED(boundaryConditionsNeumann)) THEN
        NULLIFY(boundaryConditionsNeumann%integrationMatrix)
        NULLIFY(boundaryConditionsNeumann%pointValues)
        NULLIFY(boundaryConditionsNeumann%pointDofMapping)

        numberOfLocalDofs=boundaryConditionsVariable%VARIABLE%numberOfDofs
        ALLOCATE(boundaryConditionsNeumann%setDofs(numberOfValues),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate Neumann set DOFs.",err,error,*999)
        boundaryConditionsNeumann%setDofs=0
      ELSE
        CALL FlagError("The boundary condition Neumann is not associated",err,error,*998)
      END IF
    ELSE
      CALL FlagError("Boundary conditions variable is not associated.",err,error,*998)
    END IF

    EXITS("BoundaryConditions_NeumannInitialise")
    RETURN
999 CALL BoundaryConditions_NeumannFinalise(boundaryConditionsVariable,dummyErr,dummyError,*998)
998 ERRORSEXITS("BoundaryConditions_NeumannInitialise",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_NeumannInitialise

  !
  !================================================================================================================================
  !

  !>Initialise the Neumann boundary conditions matrices and vectors.
  !>This must be done after we know which DOFs have Neumann point conditions so
  !>that we can work out the matrix sparsity pattern.
  SUBROUTINE BoundaryConditions_NeumannMatricesInitialise(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to initialise Neumann condition matrices for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: myGroupNodeNumber,numberOfGroupNodes
    INTEGER(INTG) :: numberOfPointDofs, numberNonZeros, numberRowEntries, neumannConditionNumber, localNeumannConditionIdx
    INTEGER(INTG) :: neumannIdx, globalDof, localDof, localDofNyy, domainIdx, numberOfDomains, domainNumber, componentNumber
    INTEGER(INTG) :: nodeIdx, derivIdx, nodeNumber, versionNumber, derivativeNumber, columnNodeNumber, lineIdx, faceIdx, columnDof
    INTEGER(INTG) :: dummyErr
    INTEGER(INTG), ALLOCATABLE :: rowIndices(:), columnIndices(:), localDofNumbers(:)
    REAL(DP) :: pointValue
    TYPE(BoundaryConditionsNeumannType), POINTER :: boundaryConditionsNeumann
    TYPE(FieldVariableType), POINTER :: rhsVariable
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainMappingType), POINTER :: rowMapping, pointDofMapping
    TYPE(DomainTopologyType), POINTER :: topology
    TYPE(DomainLineType), POINTER :: line
    TYPE(DomainFaceType), POINTER :: face
    TYPE(ListType), POINTER :: columnIndicesList, rowColumnIndicesList
    TYPE(VARYING_STRING) :: dummyError
    TYPE(WorkGroupType), POINTER :: workGroup

    ENTERS("BoundaryConditions_NeumannMatricesInitialise",err,error,*999)

    IF(ASSOCIATED(boundaryConditionsVariable)) THEN
      rhsVariable=>boundaryConditionsVariable%variable
      IF(.NOT.ASSOCIATED(rhsVariable)) &
        & CALL FlagError("RHS boundary conditions variable field variable is not associated.",err,error,*999)
      numberOfPointDofs=boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT) + &
        & boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
      boundaryConditionsNeumann=>boundaryConditionsVariable%neumannBoundaryConditions
      IF(ASSOCIATED(boundaryConditionsNeumann)) THEN
        ! For rows we can re-use the RHS variable row mapping
        NULLIFY(rowMapping)
        CALL FieldVariable_DomainMappingGet(rhsVariable,rowMapping,err,error,*999)
        NULLIFY(workGroup)
        CALL DomainMapping_WorkGroupGet(rowMapping,workGroup,err,error,*999)
        CALL WorkGroup_NumberOfGroupNodesGet(workGroup,numberOfGroupNodes,err,error,*999)

        ! Create a domain mapping for the Neumann point DOFs, required for the distributed matrix columns
        NULLIFY(pointDOFMapping)
        CALL DomainMapping_Initialise(pointDofMapping,err,error,*999)
        CALL DomainMapping_WorkGroupSet(pointDofMapping,workGroup,err,error,*999)
        boundaryConditionsNeumann%pointDofMapping=>pointDofMapping
        ! Calculate global to local mapping for Neumann DOFs
        pointDofMapping%numberOfGlobal=numberOfPointDofs
        ALLOCATE(pointDofMapping%globalToLocalMap(numberOfPointDofs),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate Neumann point DOF global to local mapping.",err,error,*999)
        ALLOCATE(localDofNumbers(0:numberOfGroupNodes-1),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate local Neumann DOF numbers.",err,error,*999)
        localDofNumbers=0

        IF(DIAGNOSTICS2) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Local numbering",err,error,*999)
        END IF
        DO neumannIdx=1,numberOfPointDofs
          globalDof=boundaryConditionsNeumann%setDofs(neumannIdx)
          ! Get domain information from the RHS variable domain mapping, but set new local numbers.
          numberOfDomains=rhsVariable%domainMapping%globalToLocalMap(globalDof)%numberOfDomains
          pointDofMapping%globalToLocalMap(neumannIdx)%numberOfDomains=numberOfDomains
          ALLOCATE(pointDofMapping%globalToLocalMap(neumannIdx)%localNumber(numberOfDomains),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate Neumann DOF global to local map local number.",err,error,*999)
          ALLOCATE(pointDofMapping%globalToLocalMap(neumannIdx)%domainNumber(numberOfDomains),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate Neumann DOF global to local map domain number.",err,error,*999)
          ALLOCATE(pointDofMapping%globalToLocalMap(neumannIdx)%localType(numberOfDomains),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate Neumann DOF global to local map local type.",err,error,*999)
          IF(DIAGNOSTICS2) THEN
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Neumann point DOF index = ",neumannIdx,err,error,*999)
          END IF
          DO domainIdx=1,numberOfDomains
            domainNumber=rhsVariable%domainMapping%globalToLocalMap(globalDof)%domainNumber(domainIdx)
            pointDofMapping%globalToLocalMap(neumannIdx)%domainNumber(domainIdx)=domainNumber
            pointDofMapping%globalToLocalMap(neumannIdx)%localType(domainIdx)= &
              & rhsVariable%domainMapping%globalToLocalMap(globalDof)%localType(domainIdx)
            IF(pointDofMapping%globalToLocalMap(neumannIdx)%localType(domainIdx)==DOMAIN_LOCAL_INTERNAL.OR. &
                & pointDofMapping%globalToLocalMap(neumannIdx)%localType(domainIdx)==DOMAIN_LOCAL_BOUNDARY) THEN
              localDofNumbers(domainNumber)=localDofNumbers(domainNumber)+1
              pointDofMapping%globalToLocalMap(neumannIdx)%localNumber(domainIdx)=localDofNumbers(domainNumber)
              IF(DIAGNOSTICS2) THEN
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global rhs var DOF = ",globalDof,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Domain number = ",domainNumber,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Local type = ", &
                  & pointDofMapping%globalToLocalMap(neumannIdx)%localType(domainIdx),err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Local number = ",localDofNumbers(domainNumber),err,error,*999)
              END IF
            ENDIF
          END DO
        END DO
        !Local DOFs must be numbered before ghost DOFs, so loop though again, this time numbering GHOST DOFs
        IF(DIAGNOSTICS2) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Ghost numbering",err,error,*999)
        END IF
        DO neumannIdx=1,numberOfPointDofs
          globalDof=boundaryConditionsNeumann%setDofs(neumannIdx)
          numberOfDomains=rhsVariable%domainMapping%globalToLocalMap(globalDof)%numberOfDomains
          IF(DIAGNOSTICS2) THEN
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Neumann point DOF index = ",neumannIdx,err,error,*999)
          END IF
          DO domainIdx=1,numberOfDomains
            IF(pointDofMapping%globalToLocalMap(neumannIdx)%localType(domainIdx)==DOMAIN_LOCAL_GHOST) THEN
              domainNumber=rhsVariable%domainMapping%globalToLocalMap(globalDof)%domainNumber(domainIdx)
              localDofNumbers(domainNumber)=localDofNumbers(domainNumber)+1
              pointDofMapping%globalToLocalMap(neumannIdx)%localNumber(domainIdx)=localDofNumbers(domainNumber)
              IF(DIAGNOSTICS2) THEN
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global rhs var DOF = ",globalDof,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Domain number = ",domainNumber,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Local number = ",localDofNumbers(domainNumber),err,error,*999)
              END IF
            ENDIF
          END DO
        END DO

        CALL DomainMapping_LocalFromGlobalCalculate(pointDofMapping,err,error,*999)

        CALL DistributedMatrix_CreateStart(rowMapping,pointDofMapping,boundaryConditionsNeumann%integrationMatrix,err,error,*999)
        SELECT CASE(boundaryConditionsVariable%boundaryConditions%neumannMatrixSparsity)
        CASE(BOUNDARY_CONDITION_SPARSE_MATRICES)
          ! Work out integration matrix sparsity structure
          ! For a single process, compressed column would be more memory efficient, but with
          ! multiple processes the number of Neumann point DOFs could be more than the number
          ! of local row DOFs, and multiplying a compressed row matrix by a vector is faster,
          ! so we will use compressed row storage
          ALLOCATE(rowIndices(rowMapping%totalNumberOfLocal+1),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate Neumann integration matrix column indices.",err,error,*999)
          ! We don't know the number of non zeros before hand, so use a list to keep track of column indices
          NULLIFY(columnIndicesList)
          CALL LIST_CREATE_START(columnIndicesList,err,error,*999)
          CALL LIST_DATA_TYPE_SET(columnIndicesList,LIST_INTG_TYPE,err,error,*999)
          CALL LIST_CREATE_FINISH(columnIndicesList,err,error,*999)
          ! Stores the column indices for the current row
          NULLIFY(rowColumnIndicesList)
          CALL LIST_CREATE_START(rowColumnIndicesList,err,error,*999)
          CALL LIST_DATA_TYPE_SET(rowColumnIndicesList,LIST_INTG_TYPE,err,error,*999)
          CALL LIST_MUTABLE_SET(rowColumnIndicesList,.TRUE.,err,error,*999)
          CALL LIST_CREATE_FINISH(rowColumnIndicesList,err,error,*999)
          rowIndices(1)=1

          DO localDof=1,rhsVariable%domainMapping%totalNumberOfLocal
            localDofNyy=rhsVariable%dofToParamMap%DOFType(2,localDof)
            componentNumber=rhsVariable%dofToParamMap%nodeDOF2ParamMap(4,localDofNyy)
            ! Get topology for finding faces/lines
            topology=>rhsVariable%COMPONENTS(componentNumber)%DOMAIN%TOPOLOGY
            IF(.NOT.ASSOCIATED(topology)) THEN
              CALL FlagError("Field component topology is not associated.",err,error,*999)
            END IF

            SELECT CASE(rhsVariable%COMPONENTS(componentNumber)%interpolationType)
            CASE(FIELD_NODE_BASED_INTERPOLATION)
              nodeNumber=rhsVariable%dofToParamMap%nodeDOF2ParamMap(3,localDofNyy)
              IF(.NOT.ALLOCATED(topology%NODES%NODES)) THEN
                CALL FlagError("Topology nodes are not allocated.",err,error,*999)
              END IF
              IF(topology%NODES%NODES(nodeNumber)%boundaryNode) THEN
                SELECT CASE(rhsVariable%COMPONENTS(componentNumber)%DOMAIN%numberOfDimensions)
                CASE(1)
                  ! Only one column used, as this is the same as setting an integrated
                  ! value so no other DOFs are affected
                  globalDof=rhsVariable%domainMapping%localToGlobalMap(localDof)
                  IF(boundaryConditionsVariable%conditionTypes(globalDof)==BOUNDARY_CONDITION_NEUMANN_POINT.OR. &
                      & boundaryConditionsVariable%conditionTypes(globalDof)==BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) THEN
                    ! Find the Neumann condition number
                    neumannConditionNumber=0
                    DO neumannIdx=1,numberOfPointDofs
                      IF(boundaryConditionsNeumann%setDofs(neumannIdx)==globalDof) THEN
                        neumannConditionNumber=neumannIdx
                      END IF
                    END DO
                    IF(neumannConditionNumber==0) THEN
                      CALL FlagError("Could not find matching Neuamann condition number for global DOF "// &
                        & TRIM(NumberToVString(globalDof,"*",err,error))//" with Neumann condition set.",err,error,*999)
                    ELSE
                      CALL LIST_ITEM_ADD(rowColumnIndicesList,neumannConditionNumber,err,error,*999)
                    END IF
                  END IF
                CASE(2)
                  ! Loop over all lines for this node and find any DOFs that have a Neumann point condition set
                  DO lineIdx=1,topology%NODES%NODES(nodeNumber)%numberOfNodeLines
                    IF(.NOT.ALLOCATED(topology%LINES%LINES)) THEN
                      CALL FlagError("Topology lines have not been calculated.",err,error,*999)
                    END IF
                    line=>topology%LINES%LINES(topology%NODES%NODES(nodeNumber)%nodeLines(lineIdx))
                    IF(.NOT.line%boundaryLine) CYCLE
                    DO nodeIdx=1,line%BASIS%numberOfNodes
                      columnNodeNumber=line%nodesInLine(nodeIdx)
                      DO derivIdx=1,line%BASIS%numberOfDerivatives(nodeIdx)
                        derivativeNumber=line%derivativesInLine(1,derivIdx,nodeIdx)
                        versionNumber=line%derivativesInLine(2,derivIdx,nodeIdx)
                        columnDof=rhsVariable%COMPONENTS(componentNumber)%paramToDOFMap%nodeParam2DOFMap% &
                          & NODES(columnNodeNumber)%DERIVATIVES(derivativeNumber)%VERSIONS(versionNumber)
                        globalDof=rhsVariable%domainMapping%localToGlobalMap(columnDof)
                        IF(boundaryConditionsVariable%conditionTypes(globalDof)==BOUNDARY_CONDITION_NEUMANN_POINT.OR. &
                            & boundaryConditionsVariable%conditionTypes(globalDof)== &
                            & BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) THEN
                          neumannConditionNumber=0
                          DO neumannIdx=1,numberOfPointDofs
                            IF(boundaryConditionsNeumann%setDofs(neumannIdx)==globalDof) THEN
                              neumannConditionNumber=neumannIdx
                            END IF
                          END DO
                          IF(neumannConditionNumber==0) THEN
                            CALL FlagError("Could not find matching Neuamann condition number for global DOF "// &
                              & TRIM(NumberToVString(globalDof,"*",err,error))//" with Neumann condition set.",err,error,*999)
                          ELSE
                            CALL LIST_ITEM_ADD(rowColumnIndicesList,neumannConditionNumber,err,error,*999)
                          END IF
                        END IF
                      END DO
                    END DO
                  END DO
                CASE(3)
                  ! Loop over all faces for this node and find any DOFs that have a Neumann point condition set
                  DO faceIdx=1,topology%NODES%NODES(nodeNumber)%numberOfNodeFaces
                    IF(.NOT.ALLOCATED(topology%faces%faces)) THEN
                      CALL FlagError("Topology faces have not been calculated.",err,error,*999)
                    END IF
                    face=>topology%FACES%FACES(topology%NODES%NODES(nodeNumber)%nodeFaces(faceIdx))
                    IF(.NOT.face%boundaryFace) CYCLE
                    DO nodeIdx=1,face%BASIS%numberOfNodes
                      columnNodeNumber=face%nodesInFace(nodeIdx)
                      DO derivIdx=1,face%BASIS%numberOfDerivatives(nodeIdx)
                        derivativeNumber=face%derivativesInFace(1,derivIdx,nodeIdx)
                        versionNumber=face%derivativesInFace(2,derivIdx,nodeIdx)
                        columnDof=rhsVariable%COMPONENTS(componentNumber)%paramToDOFMap%nodeParam2DOFMap% &
                          & NODES(columnNodeNumber)%DERIVATIVES(derivativeNumber)%VERSIONS(versionNumber)
                        globalDof=rhsVariable%domainMapping%localToGlobalMap(columnDof)
                        IF(boundaryConditionsVariable%conditionTypes(globalDof)==BOUNDARY_CONDITION_NEUMANN_POINT.OR. &
                            & boundaryConditionsVariable%conditionTypes(globalDof)== &
                            & BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) THEN
                          neumannConditionNumber=0
                          DO neumannIdx=1,numberOfPointDofs
                            IF(boundaryConditionsNeumann%setDofs(neumannIdx)==globalDof) THEN
                              neumannConditionNumber=neumannIdx
                            END IF
                          END DO
                          IF(neumannConditionNumber==0) THEN
                            CALL FlagError("Could not find matching Neuamann condition number for global DOF "// &
                              & TRIM(NumberToVString(globalDof,"*",err,error))//" with Neumann condition set.",err,error,*999)
                          ELSE
                            CALL LIST_ITEM_ADD(rowColumnIndicesList,neumannConditionNumber,err,error,*999)
                          END IF
                        END IF
                      END DO
                    END DO
                  END DO
                CASE DEFAULT
                  CALL FlagError("The dimension is invalid for point Neumann conditions",err,error,*999)
                END SELECT !number of dimensions
              END IF
            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(FIELD_CONSTANT_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              CALL FlagError("The interpolation type of "// &
                & TRIM(NumberToVString(rhsVariable%COMPONENTS(componentNumber) &
                & %interpolationType,"*",err,error))//" is invalid for component number "// &
                & TRIM(NumberToVString(componentNumber,"*",err,error))//".", &
                & err,error,*999)
            END SELECT

            !Sort and remove duplicates
            CALL LIST_REMOVE_DUPLICATES(rowColumnIndicesList,err,error,*999)
            !Now add all column DOFs in this row that use Neumann conditions to the overall column indices
            CALL List_AppendList(columnIndicesList,rowColumnIndicesList,err,error,*999)
            CALL LIST_NUMBER_OF_ITEMS_GET(rowColumnIndicesList,numberRowEntries,err,error,*999)
            rowIndices(localDof+1)=rowIndices(localDof)+numberRowEntries
            CALL List_ClearItems(rowColumnIndicesList,err,error,*999)
          END DO !local DOFs

          CALL LIST_DESTROY(rowColumnIndicesList,err,error,*999)
          CALL LIST_DETACH_AND_DESTROY(columnIndicesList,numberNonZeros,columnIndices,err,error,*999)
          IF(DIAGNOSTICS1) THEN
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Neumann integration matrix sparsity",err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number non-zeros = ", numberNonZeros,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number columns = ",numberOfPointDofs,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number rows = ", &
              & rhsVariable%domainMapping%totalNumberOfLocal,err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfPointDofs+1,6,6, &
              & rowIndices,'("  Row indices: ",6(X,I6))', '(6X,6(X,I6))',err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberNonZeros,6,6, &
              & columnIndices,'("  Column indices: ",6(X,I6))', '(6X,6(X,I6))',err,error,*999)
          END IF

          CALL DistributedMatrix_StorageTypeSet(boundaryConditionsNeumann%integrationMatrix, &
            & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
          CALL DistributedMatrix_NumberOfNonZerosSet(boundaryConditionsNeumann%integrationMatrix,numberNonZeros,err,error,*999)
          CALL DistributedMatrix_StorageLocationsSet(boundaryConditionsNeumann%integrationMatrix, &
            & rowIndices,columnIndices(1:numberNonZeros),err,error,*999)

          DEALLOCATE(localDofNumbers)
          DEALLOCATE(rowIndices)
          DEALLOCATE(columnIndices)
        CASE(BOUNDARY_CONDITION_FULL_MATRICES)
          CALL DistributedMatrix_StorageTypeSet(boundaryConditionsNeumann%integrationMatrix, &
            & DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,err,error,*999)
        CASE DEFAULT
          CALL FlagError("The Neumann matrix sparsity type of "// &
              & TRIM(NumberToVString(boundaryConditionsVariable%boundaryConditions%neumannMatrixSparsity,"*",err,error))// &
              & " is invalid.",err,error,*999)
        END SELECT
        CALL DistributedMatrix_TransposeTypeSet(boundaryConditionsNeumann%integrationMatrix, &
          & DISTRIBUTED_MATRIX_NO_TRANSPOSE_REQUIRED,err,error,*999)
        CALL DistributedMatrix_CreateFinish(boundaryConditionsNeumann%integrationMatrix,err,error,*999)

        !Set up vector of Neumann point values
        CALL DistributedVector_CreateStart(pointDofMapping,boundaryConditionsNeumann%pointValues,err,error,*999)
        CALL DistributedVector_CreateFinish(boundaryConditionsNeumann%pointValues,err,error,*999)
        NULLIFY(decomposition)
        CALL Field_DecompositionGet(rhsVariable%field,decomposition,err,error,*999)
        NULLIFY(workGroup)
        CALL Decomposition_WorkGroupGet(decomposition,workGroup,err,error,*999)
        CALL WorkGroup_GroupNodeNumberGet(workGroup,myGroupNodeNumber,err,error,*999)
        !Set point values vector from boundary conditions field parameter set
        DO neumannIdx=1,numberOfPointDofs
          globalDof=boundaryConditionsNeumann%setDofs(neumannIdx)
          IF(rhsVariable%domainMapping%globalToLocalMap(globalDof)%domainNumber(1)==myGroupNodeNumber) THEN
            localDof=rhsVariable%domainMapping%globalToLocalMap(globalDof)%localNumber(1)
            ! Set point DOF vector value
            localNeumannConditionIdx=boundaryConditionsNeumann%pointDofMapping%globalToLocalMap(neumannIdx)%localNumber(1)
            CALL Field_ParameterSetGetLocalDOF(rhsVariable%FIELD,rhsVariable%variableType, &
              & FIELD_BOUNDARY_CONDITIONS_SET_TYPE,localDof,pointValue,err,error,*999)
            CALL DistributedVector_ValuesSet(boundaryConditionsNeumann%pointValues, &
              & localNeumannConditionIdx,pointValue,err,error,*999)
          END IF
        END DO
        CALL DistributedVector_UpdateStart(boundaryConditionsNeumann%pointValues,err,error,*999)
        CALL DistributedVector_UpdateFinish(boundaryConditionsNeumann%pointValues,err,error,*999)

      ELSE
        CALL FlagError("The boundary condition Neumann is not associated",err,error,*998)
      END IF
    ELSE
      CALL FlagError("Boundary conditions variable is not associated.",err,error,*998)
    END IF

    EXITS("BoundaryConditions_NeumannMatricesInitialise")
    RETURN
999 IF(ALLOCATED(rowIndices)) THEN
      DEALLOCATE(rowIndices)
    END IF
    IF(ALLOCATED(columnIndices)) THEN
      DEALLOCATE(columnIndices)
    END IF
    IF(ALLOCATED(localDofNumbers)) THEN
      DEALLOCATE(localDofNumbers)
    END IF
    CALL BoundaryConditions_NeumannMatricesFinalise(boundaryConditionsVariable,dummyErr,dummyError,*998)
998 ERRORS("BoundaryConditions_NeumannMatricesInitialise",err,error)
    EXITS("BoundaryConditions_NeumannMatricesInitialise")
    RETURN 1

  END SUBROUTINE BoundaryConditions_NeumannMatricesInitialise

  !
  !================================================================================================================================
  !

  !Finalise the Neumann condition information for a boundary conditions variable
  SUBROUTINE BoundaryConditions_NeumannFinalise(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to finalise the Neumann conditions for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(BoundaryConditionsNeumannType), POINTER :: boundaryConditionsNeumann

    ENTERS("BoundaryConditions_NeumannFinalise",err,error,*999)

    IF(ASSOCIATED(boundaryConditionsVariable)) THEN
      boundaryConditionsNeumann=>boundaryConditionsVariable%neumannBoundaryConditions
      IF(ASSOCIATED(boundaryConditionsNeumann)) THEN
        IF(ALLOCATED(boundaryConditionsNeumann%setDofs)) &
          & DEALLOCATE(boundaryConditionsNeumann%setDofs)
        CALL BoundaryConditions_NeumannMatricesFinalise(boundaryConditionsVariable,err,error,*999)
        DEALLOCATE(boundaryConditionsNeumann)
        NULLIFY(boundaryConditionsVariable%neumannBoundaryConditions)
      END IF
    ELSE
      CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
    END IF

    EXITS("BoundaryConditions_NeumannFinalise")
    RETURN
999 ERRORSEXITS("BoundaryConditions_NeumannFinalise",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_NeumannFinalise

  !
  !================================================================================================================================
  !

  !Finalise the Neumann condition matrices for a boundary conditions variable
  SUBROUTINE BoundaryConditions_NeumannMatricesFinalise(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to finalise Neumann condition matrices for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(BoundaryConditionsNeumannType), POINTER :: boundaryConditionsNeumann

    ENTERS("BoundaryConditions_NeumannMatricesFinalise",err,error,*999)

    IF(ASSOCIATED(boundaryConditionsVariable)) THEN
      boundaryConditionsNeumann=>boundaryConditionsVariable%neumannBoundaryConditions
      IF(ASSOCIATED(boundaryConditionsNeumann)) THEN
        IF(ASSOCIATED(boundaryConditionsNeumann%integrationMatrix)) &
          & CALL DistributedMatrix_Destroy(boundaryConditionsNeumann%integrationMatrix,err,error,*999)
        IF(ASSOCIATED(boundaryConditionsNeumann%pointValues)) &
          & CALL DistributedVector_Destroy(boundaryConditionsNeumann%pointValues,err,error,*999)
        CALL DomainMapping_Finalise(boundaryConditionsNeumann%pointDofMapping,err,error,*999)
      END IF
    ELSE
      CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
    END IF

    EXITS("BoundaryConditions_NeumannMatricesFinalise")
    RETURN
999 ERRORSEXITS("BoundaryConditions_NeumannMatricesFinalise",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_NeumannMatricesFinalise

  !
  !================================================================================================================================
  !

  !>Calculates integrated Neumann condition values from point values for a boundary conditions variable and
  !>updates the FIELD_INTEGRATED_NEUMANN_SET_TYPE parameter set for the field variable.
  SUBROUTINE BoundaryConditions_NeumannIntegrate(rhsBoundaryConditions,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionVariableType), POINTER, INTENT(IN) :: rhsBoundaryConditions !<The boundary conditions for the right hand side field variable
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    !Local variables
    INTEGER(INTG) :: componentNumber,globalDof,localDof,neumannDofIdx,myGroupNodeNumber
    INTEGER(INTG) :: numberOfNeumann,neumannLocalDof,neumannDofNyy
    INTEGER(INTG) :: neumannGlobalDof,neumannNodeNumber,neumannLocalNodeNumber,neumannLocalDerivNumber
    INTEGER(INTG) :: faceIdx,lineIdx,nodeIdx,derivIdx,gaussIdx
    INTEGER(INTG) :: faceNumber,lineNumber
    INTEGER(INTG) :: ms,os,nodeNumber,derivativeNumber,versionNumber
    LOGICAL :: dependentGeometry
    REAL(DP) :: integratedValue,phim,phio
    TYPE(BoundaryConditionsNeumannType), POINTER :: neumannConditions
    TYPE(BasisType), POINTER :: basis
    TYPE(FieldType), POINTER :: geometricField
    TYPE(FieldVariableType), POINTER :: rhsVariable
    TYPE(FieldInterpolatedPointMetricsPtrType), POINTER :: interpolatedPointMetrics(:)
    TYPE(FieldInterpolatedPointPtrType), POINTER :: interpolatedPoints(:)
    TYPE(FieldInterpolationParametersPtrType), POINTER :: interpolationParameters(:), scalingParameters(:)
    TYPE(DistributedVectorType), POINTER :: integratedValues
    TYPE(DomainTopologyType), POINTER :: topology
    TYPE(DomainFacesType), POINTER :: faces
    TYPE(DomainLinesType), POINTER :: lines
    TYPE(DomainFaceType), POINTER :: face
    TYPE(DomainLineType), POINTER :: line
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(QuadratureSchemeType), POINTER :: quadratureScheme
    TYPE(WorkGroupType), POINTER :: workGroup

    ENTERS("BoundaryConditions_NeumannIntegrate",err,error,*999)

    NULLIFY(scalingParameters)
    NULLIFY(interpolationParameters)
    NULLIFY(interpolatedPoints)
    NULLIFY(interpolatedPointMetrics)
    NULLIFY(integratedValues)

    neumannConditions=>rhsBoundaryConditions%neumannBoundaryConditions
    !Check that Neumann conditions are associated, otherwise do nothing
    IF(ASSOCIATED(neumannConditions)) THEN
      rhsVariable=>rhsBoundaryConditions%VARIABLE
      IF(.NOT.ASSOCIATED(rhsVariable)) THEN
        CALL FlagError("Field variable for RHS boundary conditions is not associated.",err,error,*999)
      END IF

      NULLIFY(geometricField)
      CALL Field_GeometricGeneralFieldGet(rhsVariable%field,geometricField,dependentGeometry,err,error,*999)

      CALL DistributedMatrix_AllValuesSet(neumannConditions%integrationMatrix,0.0_DP,err,error,*999)

      numberOfNeumann=rhsBoundaryConditions%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT) + &
        & rhsBoundaryConditions%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)

      NULLIFY(decomposition)
      CALL Field_DecompositionGet(rhsVariable%field,decomposition,err,error,*999)
      NULLIFY(workGroup)
      CALL Decomposition_WorkGroupGet(decomposition,workGroup,err,error,*999)
      CALL WorkGroup_GroupNodeNumberGet(workGroup,myGroupNodeNumber,err,error,*999)

      ! Initialise field interpolation parameters for the geometric field, which are required for the
      ! face/line Jacobian and scale factors
      CALL Field_InterpolationParametersInitialise(geometricField,interpolationParameters,err,error,*999)
      CALL Field_InterpolationParametersInitialise(rhsVariable%field,scalingParameters,err,error,*999)
      CALL Field_InterpolatedPointsInitialise(interpolationParameters,interpolatedPoints,err,error,*999)
      CALL Field_InterpolatedPointsMetricsInitialise(interpolatedPoints,interpolatedPointMetrics,err,error,*999)

      ! Loop over all Neumann point DOFs, finding the boundary lines or faces they are on
      ! and integrating over them
      DO neumannDofIdx=1,numberOfNeumann
        neumannGlobalDof=neumannConditions%setDofs(neumannDofIdx)
        IF(rhsVariable%domainMapping%globalToLocalMap(neumannGlobalDof)%domainNumber(1)==myGroupNodeNumber) THEN
          neumannLocalDof=rhsVariable%domainMapping%globalToLocalMap(neumannGlobalDof)%localNumber(1)
          ! Get Neumann DOF component and topology for that component
          neumannDofNyy=rhsVariable%dofToParamMap%DOFType(2,neumannLocalDof)
          componentNumber=rhsVariable%dofToParamMap%nodeDOF2ParamMap(4,neumannDofNyy)
          topology=>rhsVariable%COMPONENTS(componentNumber)%DOMAIN%TOPOLOGY
          IF(.NOT.ASSOCIATED(topology)) THEN
            CALL FlagError("Field component topology is not associated.",err,error,*999)
          END IF
          decomposition=>rhsVariable%COMPONENTS(componentNumber)%DOMAIN%DECOMPOSITION
          IF(.NOT.ASSOCIATED(decomposition)) THEN
            CALL FlagError("Field component decomposition is not associated.",err,error,*999)
          END IF
          SELECT CASE(rhsVariable%COMPONENTS(componentNumber)%interpolationType)
          CASE(FIELD_NODE_BASED_INTERPOLATION)
            neumannNodeNumber=rhsVariable%dofToParamMap%nodeDOF2ParamMap(3,neumannDofNyy)
            SELECT CASE(rhsVariable%COMPONENTS(componentNumber)%domain%numberOfDimensions)
            CASE(1)
              CALL DistributedMatrix_ValuesSet(neumannConditions%integrationMatrix,neumannLocalDof,neumannDofIdx, &
                & 1.0_DP,err,error,*999)
            CASE(2)
              IF(.NOT.decomposition%calculateLines) THEN
                CALL FlagError("Decomposition does not have lines calculated.",err,error,*999)
              END IF
              lines=>topology%LINES
              IF(.NOT.ASSOCIATED(lines)) THEN
                CALL FlagError("Mesh topology lines is not associated.",err,error,*999)
              END IF
              linesLoop: DO lineIdx=1,topology%NODES%NODES(neumannNodeNumber)%numberOfNodeLines
                lineNumber=topology%NODES%NODES(neumannNodeNumber)%nodeLines(lineIdx)
                line=>topology%lines%lines(lineNumber)
                IF(.NOT.line%boundaryLine) &
                  CYCLE linesLoop
                basis=>line%basis
                IF(.NOT.ASSOCIATED(basis)) THEN
                  CALL FlagError("Line basis is not associated.",err,error,*999)
                END IF
                neumannLocalNodeNumber=0
                neumannLocalDerivNumber=0
                ! Check all nodes in line to find the local numbers for the Neumann DOF, and
                ! make sure we don't have an integrated_only condition set on the line
                DO nodeIdx=1,line%BASIS%numberOfNodes
                  nodeNumber=line%nodesInLine(nodeIdx)
                  DO derivIdx=1,line%BASIS%numberOfDerivatives(nodeIdx)
                    derivativeNumber=line%derivativesInLine(1,derivIdx,nodeIdx)
                    versionNumber=line%derivativesInLine(2,derivIdx,nodeIdx)
                    localDof=rhsVariable%COMPONENTS(componentNumber)%paramToDOFMap%nodeParam2DOFMap% &
                      & NODES(nodeNumber)%DERIVATIVES(derivativeNumber)%VERSIONS(versionNumber)
                    globalDof=rhsVariable%domainMapping%localToGlobalMap(localDof)
                    IF(globalDof==neumannGlobalDof) THEN
                      neumannLocalNodeNumber=nodeIdx
                      neumannLocalDerivNumber=derivIdx
                    ELSE IF(rhsBoundaryConditions%conditionTypes(globalDof)==BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY) THEN
                      CYCLE linesLoop
                    END IF
                  END DO
                END DO
                IF(neumannLocalNodeNumber==0) THEN
                  CALL FlagError("Could not find local Neumann node and derivative numbers in line.",err,error,*999)
                END IF

                ! Now perform actual integration
                quadratureScheme=>basis%quadrature%quadratureSchemeMap(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                IF(.NOT.ASSOCIATED(quadratureScheme)) THEN
                  CALL FlagError("Line basis default quadrature scheme is not associated.",err,error,*999)
                END IF
                CALL Field_InterpolationParametersLineGet(FIELD_VALUES_SET_TYPE,lineNumber, &
                  & interpolationParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                IF(rhsVariable%FIELD%SCALINGS%scalingType/=FIELD_NO_SCALING) THEN
                  CALL Field_InterpolationParametersScaleFactorsLineGet(lineNumber, &
                    & scalingParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                END IF

                DO nodeIdx=1,line%BASIS%numberOfNodes
                  nodeNumber=line%nodesInLine(nodeIdx)
                  DO derivIdx=1,line%BASIS%numberOfDerivatives(nodeIdx)
                    derivativeNumber=line%derivativesInLine(1,derivIdx,nodeIdx)
                    versionNumber=line%derivativesInLine(2,derivIdx,nodeIdx)
                    localDof=rhsVariable%COMPONENTS(componentNumber)%paramToDOFMap%nodeParam2DOFMap% &
                      & NODES(nodeNumber)%DERIVATIVES(derivativeNumber)%VERSIONS(versionNumber)

                    ms=basis%elementParameterIndex(derivIdx,nodeIdx)
                    os=basis%elementParameterIndex(neumannLocalDerivNumber,neumannLocalNodeNumber)

                    integratedValue=0.0_DP
                    ! Loop over line gauss points, adding gauss weighted terms to the integral
                    DO gaussIdx=1,quadratureScheme%numberOfGauss
                      CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
                        & interpolatedPoints(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                      CALL Field_InterpolatedPointMetricsCalculate(COORDINATE_JACOBIAN_LINE_TYPE, &
                        & interpolatedPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)

                      !Get basis function values at guass points
                      phim=quadratureScheme%gaussBasisFunctions(ms,NO_PART_DERIV,gaussIdx)
                      phio=quadratureScheme%gaussBasisFunctions(os,NO_PART_DERIV,gaussIdx)

                      !Add gauss point value to total line integral
                      integratedValue=integratedValue+phim*phio* &
                        & quadratureScheme%gaussWeights(gaussIdx)* &
                        & interpolatedPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%jacobian
                    END DO

                    ! Multiply by scale factors for dependent variable
                    IF(rhsVariable%FIELD%SCALINGS%scalingType/=FIELD_NO_SCALING) THEN
                      integratedValue=integratedValue* &
                        & scalingParameters(FIELD_U_VARIABLE_TYPE)%ptr%scaleFactors(ms,componentNumber)* &
                        & scalingParameters(FIELD_U_VARIABLE_TYPE)%ptr%scaleFactors(os,componentNumber)
                    END IF

                    ! Add integral term to N matrix
                    CALL DistributedMatrix_ValuesAdd(neumannConditions%integrationMatrix,localDof,neumannDofIdx, &
                      & integratedValue,err,error,*999)
                  END DO
                END DO
              END DO linesLoop
            CASE(3)
              IF(.NOT.decomposition%calculateFaces) THEN
                CALL FlagError("Decomposition does not have faces calculated.",err,error,*999)
              END IF
              faces=>topology%FACES
              IF(.NOT.ASSOCIATED(faces)) THEN
                CALL FlagError("Mesh topology faces is not associated.",err,error,*999)
              END IF
              facesLoop: DO faceIdx=1,topology%NODES%NODES(neumannNodeNumber)%numberOfNodeFaces
                faceNumber=topology%NODES%NODES(neumannNodeNumber)%nodeFaces(faceIdx)
                face=>topology%FACES%FACES(faceNumber)
                IF(.NOT.face%boundaryFace) &
                  CYCLE facesLoop
                basis=>face%BASIS
                IF(.NOT.ASSOCIATED(basis)) THEN
                  CALL FlagError("Line face is not associated.",err,error,*999)
                END IF
                neumannLocalNodeNumber=0
                neumannLocalDerivNumber=0
                ! Check all nodes in the face to find the local numbers for the Neumann DOF, and
                ! make sure we don't have an integrated_only condition set on the face
                DO nodeIdx=1,basis%numberOfNodes
                  nodeNumber=face%nodesInFace(nodeIdx)
                  DO derivIdx=1,basis%numberOfDerivatives(nodeIdx)
                    derivativeNumber=face%derivativesInFace(1,derivIdx,nodeIdx)
                    versionNumber=face%derivativesInFace(2,derivIdx,nodeIdx)
                    localDof=rhsVariable%COMPONENTS(componentNumber)%paramToDOFMap%nodeParam2DOFMap% &
                      & NODES(nodeNumber)%DERIVATIVES(derivativeNumber)%VERSIONS(versionNumber)
                    globalDof=rhsVariable%domainMapping%localToGlobalMap(localDof)
                    IF(globalDof==neumannGlobalDof) THEN
                      neumannLocalNodeNumber=nodeIdx
                      neumannLocalDerivNumber=derivIdx
                    ELSE IF(rhsBoundaryConditions%conditionTypes(globalDof)==BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY) THEN
                      CYCLE facesLoop
                    END IF
                  END DO
                END DO
                IF(neumannLocalNodeNumber==0) THEN
                  CALL FlagError("Could not find local Neumann node and derivative numbers in line.",err,error,*999)
                END IF

                ! Now perform actual integration
                quadratureScheme=>basis%quadrature%quadratureSchemeMap(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                IF(.NOT.ASSOCIATED(quadratureScheme)) THEN
                  CALL FlagError("Face basis default quadrature scheme is not associated.",err,error,*999)
                END IF
                CALL Field_InterpolationParametersFaceGet(FIELD_VALUES_SET_TYPE,faceNumber, &
                  & interpolationParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                IF(rhsVariable%FIELD%SCALINGS%scalingType/=FIELD_NO_SCALING) THEN
                  CALL Field_InterpolationParametersScaleFactorsFaceGet(faceNumber, &
                    & scalingParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                END IF

                DO nodeIdx=1,basis%numberOfNodes
                  nodeNumber=face%nodesInFace(nodeIdx)
                  DO derivIdx=1,basis%numberOfDerivatives(nodeIdx)
                    derivativeNumber=face%derivativesInFace(1,derivIdx,nodeIdx)
                    versionNumber=face%derivativesInFace(2,derivIdx,nodeIdx)
                    localDof=rhsVariable%COMPONENTS(componentNumber)%paramToDOFMap%nodeParam2DOFMap% &
                      & NODES(nodeNumber)%DERIVATIVES(derivativeNumber)%VERSIONS(versionNumber)

                    ms=basis%elementParameterIndex(derivIdx,nodeIdx)
                    os=basis%elementParameterIndex(neumannLocalDerivNumber,neumannLocalNodeNumber)

                    integratedValue=0.0_DP
                    ! Loop over line gauss points, adding gauss weighted terms to the integral
                    DO gaussIdx=1,quadratureScheme%numberOfGauss
                      CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
                        & interpolatedPoints(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                      CALL Field_InterpolatedPointMetricsCalculate(COORDINATE_JACOBIAN_AREA_TYPE, &
                        & interpolatedPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)

                      !Get basis function values at guass points
                      phim=quadratureScheme%gaussBasisFunctions(ms,NO_PART_DERIV,gaussIdx)
                      phio=quadratureScheme%gaussBasisFunctions(os,NO_PART_DERIV,gaussIdx)

                      !Add gauss point value to total line integral
                      integratedValue=integratedValue+phim*phio* &
                        & quadratureScheme%gaussWeights(gaussIdx)* &
                        & interpolatedPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%jacobian
                    END DO

                    ! Multiply by scale factors
                    IF(rhsVariable%FIELD%SCALINGS%scalingType/=FIELD_NO_SCALING) THEN
                      integratedValue=integratedValue* &
                        & scalingParameters(FIELD_U_VARIABLE_TYPE)%ptr%scaleFactors(ms,componentNumber)* &
                        & scalingParameters(FIELD_U_VARIABLE_TYPE)%ptr%scaleFactors(os,componentNumber)
                    END IF

                    ! Add integral term to N matrix
                    CALL DistributedMatrix_ValuesAdd(neumannConditions%integrationMatrix,localDof,neumannDofIdx, &
                      & integratedValue,err,error,*999)
                  END DO
                END DO
              END DO facesLoop
            CASE DEFAULT
              CALL FlagError("The dimension is invalid for point Neumann conditions",err,error,*999)
            END SELECT
          CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(FIELD_CONSTANT_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            CALL FlagError("The interpolation type of "// &
              & TRIM(NumberToVString(rhsVariable%COMPONENTS(componentNumber) &
              & %interpolationType,"*",err,error))//" is invalid for component number "// &
              & TRIM(NumberToVString(componentNumber,"*",err,error))//".", &
              & err,error,*999)
          END SELECT
        END IF
      END DO

      CALL DistributedMatrix_UpdateStart(neumannConditions%integrationMatrix,err,error,*999)
      CALL DistributedMatrix_UpdateFinish(neumannConditions%integrationMatrix,err,error,*999)

      CALL Field_ParameterSetVectorGet(rhsVariable%field,rhsVariable%variableType,FIELD_INTEGRATED_NEUMANN_SET_TYPE, &
        & integratedValues,err,error,*999)
      CALL DistributedVector_AllValuesSet(integratedValues,0.0_DP,err,error,*999)
      ! Perform matrix multiplication, f = N q, to calculate force vector from integration matrix and point values
      CALL DistributedMatrix_MatrixByVectorAdd(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,1.0_DP, &
        & neumannConditions%integrationMatrix,.FALSE.,neumannConditions%pointValues,integratedValues, &
        & err,error,*999)

      CALL Field_ParameterSetUpdateStart(rhsVariable%FIELD,rhsVariable%variableType,FIELD_INTEGRATED_NEUMANN_SET_TYPE, &
        & err,error,*999)
      IF(DIAGNOSTICS1) THEN
        IF(dependentGeometry) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Using dependent field geometry",err,error,*999)
        ELSE
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Using undeformed geometry",err,error,*999)
        END IF
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfNeumann,6,6,neumannConditions%setDofs, &
          & '("  setDofs:",6(X,I8))', '(10X,6(X,I8))',err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Neumann point values",err,error,*999)
        CALL DistributedVector_Output(DIAGNOSTIC_OUTPUT_TYPE,neumannConditions%pointValues,err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Neumann integration matrix",err,error,*999)
        CALL DistributedMatrix_Output(DIAGNOSTIC_OUTPUT_TYPE,neumannConditions%integrationMatrix,err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Integrated values",err,error,*999)
        CALL DistributedVector_Output(DIAGNOSTIC_OUTPUT_TYPE,integratedValues,err,error,*999)
      END IF
      CALL Field_ParameterSetUpdateFinish(rhsVariable%FIELD,rhsVariable%variableType,FIELD_INTEGRATED_NEUMANN_SET_TYPE, &
        & err,error,*999)

    END IF !Neumann conditions associated

    EXITS("BoundaryConditions_NeumannIntegrate")
    RETURN
999 ERRORSEXITS("BoundaryConditions_NeumannIntegrate",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_NeumannIntegrate

  !
  !================================================================================================================================
  !

  !>Sets/changes the sparsity type for the Neumann integration matrices
  SUBROUTINE BoundaryConditions_NeumannSparsityTypeSet(boundaryConditions,sparsityType,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: sparsityType !<The matrix sparsity type to be set \see SolverRoutines_SparsityTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions

    ENTERS("BoundaryConditions_NeumannSparsityTypeSet",err,error,*999)

    IF(ASSOCIATED(boundaryConditions)) THEN
      SELECT CASE(sparsityType)
      CASE(BOUNDARY_CONDITION_SPARSE_MATRICES)
        boundaryConditions%neumannMatrixSparsity=BOUNDARY_CONDITION_SPARSE_MATRICES
      CASE(BOUNDARY_CONDITION_FULL_MATRICES)
        boundaryConditions%neumannMatrixSparsity=BOUNDARY_CONDITION_FULL_MATRICES
      CASE DEFAULT
        CALL FlagError("The specified Neumann integration matrix sparsity type of "// &
          & TRIM(NumberToVString(sparsityType,"*",err,error))//" is invalid.",err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Boundary conditions are not associated.",err,error,*999)
    END IF

    EXITS("BoundaryConditions_NeumannSparsityTypeSet")
    RETURN
999 ERRORSEXITS("BoundaryConditions_NeumannSparsityTypeSet",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditions_NeumannSparsityTypeSet

  !
  !================================================================================================================================
  !

  !>Sets a boundary condition on the specified user node. \see OPENCMISS_CMISSBoundaryConditionsSetNode
  SUBROUTINE BoundaryConditions_SetNode(boundaryConditions,field,variableType,versionNumber,derivativeNumber, &
    & userNodeNumber,componentNumber,condition,bcValue,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldType), POINTER :: field !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: versionNumber !<The derivative version to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: derivativeNumber !<The derivative to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValue !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: localDOF,globalDOF
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_SetNode",err,error,*999)

    NULLIFY(fieldVariable)

    CALL BoundaryConditions_AssertNotFinished(boundaryConditions,err,error,*999)
    IF(ASSOCIATED(field)) THEN
      CALL Field_ComponentDOFGetUserNode(field,variableType,versionNumber,derivativeNumber, &
        & userNodeNumber,componentNumber,localDOF,globalDOF,err,error,*999)
      CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
      IF(ASSOCIATED(fieldVariable)) THEN
        NULLIFY(boundaryConditionsVariable)
        CALL BoundaryConditions_VariableGet(boundaryConditions,fieldVariable,boundaryConditionsVariable, &
          & err,error,*999)
        CALL BoundaryConditions_CheckInterpolationType(condition,field,variableType,componentNumber,err,error,*999)
        CALL BoundaryConditions_SetLocalDOF(boundaryConditions,field,variableType, &
          & localDOF,condition,bcValue,err,error,*999)
      ELSE
        CALL FlagError("The dependent field variable is not associated",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("The dependent field is not associated",err,error,*999)
    ENDIF
    
    EXITS("BoundaryConditions_SetNode")
    RETURN
999 ERRORSEXITS("BoundaryConditions_SetNode",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_SetNode

  !
  !================================================================================================================================
  !

  !>Constrain multiple equations dependent field DOFs to be a single solver DOF in the solver equations
  SUBROUTINE BoundaryConditions_ConstrainDofsEqual(boundaryConditions,fieldVariable,globalDofs,coefficient,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER, INTENT(IN) :: boundaryConditions !<The boundary conditions for the solver equations in which to constrain the DOF.
    TYPE(FieldVariableType), POINTER, INTENT(IN) :: fieldVariable !<A pointer to the field variable containing the DOFs.
    INTEGER(INTG), INTENT(IN) :: globalDofs(:) !<The global DOFs to be constrained to be equal.
    REAL(DP), INTENT(IN) :: coefficient !<The coefficient of constraint.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message.
    !Local variables
    INTEGER(INTG) :: numberOfDofs,dofIdx,dofIdx2

    ENTERS("BoundaryConditions_ConstrainDofsEqual",err,error,*999)

    numberOfDofs=SIZE(globalDofs,1)
    IF(numberOfDofs<2) THEN
      CALL FlagError("Cannot constrain zero or 1 DOF to be equal.",err,error,*999)
    END IF

    !Check for duplicate DOFs
    DO dofIdx=1,numberOfDofs
      DO dofIdx2=dofIdx+1,numberOfDofs
        IF(globalDofs(dofIdx)==globalDofs(dofIdx2)) THEN
          CALL FlagError("DOF number "//TRIM(NumberToVstring(globalDofs(dofIdx),"*",err,error))// &
            & " is duplicated in the DOFs constrained to be equal.",err,error,*999)
        END IF
      END DO
    END DO

    !Add new DOF constraints
    !We set all DOFs except the first to be equal to coefficient * the first DOF
    !The first DOF is left unconstrained
    DO dofIdx=2,numberOfDofs
      CALL BoundaryConditions_DofConstraintSet( &
        & boundaryConditions,fieldVariable,globalDofs(dofIdx),[globalDofs(1)],[coefficient],err,error,*999)
    END DO

    EXITS("BoundaryConditions_ConstrainDofsEqual")
    RETURN
999 ERRORSEXITS("BoundaryConditions_ConstrainDofsEqual",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_ConstrainDofsEqual

  !
  !================================================================================================================================
  !

  !>Constrain multiple nodal equations dependent field DOFs to be a single solver DOF in the solver equations
  SUBROUTINE BoundaryConditions_ConstrainNodeDofsEqual(boundaryConditions,field,variableType,versionNumber,derivativeNumber, &
    & component,nodes,coefficient,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER, INTENT(IN) :: boundaryConditions !<The solver equations boundary conditions to constrain the DOFs for.
    TYPE(FieldType), POINTER, INTENT(IN) :: field !<The equations dependent field containing the field DOFs to be constrained.
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type of the DOFs to be constrained. \see OPENCMISS_FieldVariableTypes
    INTEGER(INTG), INTENT(IN) :: versionNumber !<The derivative version number.
    INTEGER(INTG), INTENT(IN) :: derivativeNumber !<The derivative number.
    INTEGER(INTG), INTENT(IN) :: component !<The field component number of the DOFs to be constrained.
    INTEGER(INTG), INTENT(IN) :: nodes(:) !<The user numbers of the nodes to be constrained to be equal.
    REAL(DP), INTENT(IN) :: coefficient !<The coefficient of constraint, applied to all but the first node.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message.
    !Local variables
    TYPE(FieldVariableType), POINTER :: fieldVariable
    INTEGER(INTG) :: numberOfNodes, nodeIdx, dof
    INTEGER(INTG), ALLOCATABLE :: globalDofs(:)

    ENTERS("BoundaryConditions_ConstrainNodeDofsEqual",err,error,*998)

    NULLIFY(fieldVariable)

    IF(.NOT.ASSOCIATED(boundaryConditions)) THEN
      CALL FlagError("Boundary conditions are not associated.",err,error,*998)
    END IF

    numberOfNodes=SIZE(nodes,1)
    ALLOCATE(globalDofs(numberOfNodes),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate equal global DOFs array.",err,error,*998)
    !Get field DOFs for nodes
    DO nodeIdx=1,numberOfNodes
      CALL Field_ComponentDOFGetUserNode(field,variableType,versionNumber,derivativeNumber,nodes(nodeIdx), &
        & component,dof,globalDofs(nodeIdx),err,error,*999)
    END DO
    !Get the field variable and boundary conditions variable for the field
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)

    !Now set DOF constraint
    CALL BoundaryConditions_ConstrainDofsEqual(boundaryConditions,fieldVariable,globalDofs,coefficient,err,error,*999)

    DEALLOCATE(globalDofs)

    EXITS("BoundaryConditions_ConstrainNodeDofsEqual")
    RETURN
999 IF(ALLOCATED(globalDofs)) DEALLOCATE(globalDofs)
998 ERRORSEXITS("BoundaryConditions_ConstrainNodeDofsEqual",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_ConstrainNodeDofsEqual

  !
  !================================================================================================================================
  !

  !>Constrain a DOF to be a linear combination of other DOFs.
  SUBROUTINE BoundaryConditions_DofConstraintSet(boundaryConditions,fieldVariable,globalDof,dofs,coefficients,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER, INTENT(IN) :: boundaryConditions !<The boundary conditions for the solver equations in which to constrain the DOF.
    TYPE(FieldVariableType), POINTER, INTENT(IN) :: fieldVariable !<A pointer to the field variable containing the DOFs.
    INTEGER(INTG), INTENT(IN) :: globalDof !<The global DOF to set the constraint on.
    INTEGER(INTG), INTENT(IN) :: dofs(:) !<The global DOFs that this DOF is constrained to depend on.
    REAL(DP), INTENT(IN) :: coefficients(:) !<The coefficient values in the DOF constraint.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message.
    !Local variables
    INTEGER(INTG) :: numberOfDofs,dofIdx,dofIdx2
    TYPE(BoundaryConditionsDofConstraintPtrType), ALLOCATABLE :: newConstraints(:)
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints
    TYPE(BoundaryConditionsDofConstraintType), POINTER :: dofConstraint

    NULLIFY(dofConstraint)
    NULLIFY(dofConstraints)

    ENTERS("BoundaryConditions_DofConstraintSet",err,error,*998)

    !Check pointers for association
    CALL BoundaryConditions_AssertNotFinished(boundaryConditions,err,error,*999)
    IF(.NOT.ASSOCIATED(fieldVariable)) THEN
      CALL FlagError("Field variable is not associated.",err,error,*998)
    END IF
    NULLIFY(boundaryConditionsVariable)
    CALL boundary_conditions_variable_get(boundaryConditions,fieldVariable,boundaryConditionsVariable,err,error,*998)
    dofConstraints=>boundaryConditionsVariable%dofConstraints
    IF(.NOT.ASSOCIATED(dofConstraints)) THEN
      CALL FlagError("Boundary conditions DOF constraints are not associated.",err,error,*998)
    END IF

    numberOfDofs=SIZE(dofs,1)
    IF(numberOfDofs==0) THEN
      CALL FlagError("Empty DOFs list.",err,error,*998)
    ELSE IF(numberOfDofs/=SIZE(coefficients,1)) THEN
      CALL FlagError("Length of coefficients does not match length of DOFs array.",err,error,*998)
    ELSE IF(numberOfDofs>1) THEN
      CALL FlagError("Support for constraining an equations DOF to be depended on multiple "// &
        & "other DOFs is not yet implemented.",err,error,*998)
    END IF

    !Check for duplicate DOFs
    DO dofIdx=1,numberOfDofs
      DO dofIdx2=dofIdx+1,numberOfDofs
        IF(dofs(dofIdx)==dofs(dofIdx2)) THEN
          CALL FlagError("DOF number "//TRIM(NumberToVstring(dofs(dofIdx),"*",err,error))// &
            & " is duplicated in the DOF constraint.",err,error,*998)
        END IF
      END DO
    END DO

    !Check DOFs are free
    DO dofIdx=1,numberOfDofs
      IF(boundaryConditionsVariable%DOFTypes(dofs(dofIdx))/=BOUNDARY_CONDITION_DOF_FREE) THEN
        CALL FlagError("DOF number "//TRIM(NumberToVstring(dofs(dofIdx),"*",err,error))// &
          & " is not free in the boundary conditions.",err,error,*998)
      END IF
    END DO

    !Allocate new DOF constraints and copy over old constraints
    ALLOCATE(newConstraints(dofConstraints%numberOfConstraints+1),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate new DOF constraints array.",err,error,*998)
    IF(dofConstraints%numberOfConstraints>0) THEN
      newConstraints(1:dofConstraints%numberOfConstraints)= &
        & dofConstraints%constraints(1:dofConstraints%numberOfConstraints)
    END IF

    !Set the new DOF constraint
    ALLOCATE(dofConstraint,stat=err)
    IF(err/=0) CALL FlagError("Could not allocate new DOF constraint.",err,error,*999)
    ALLOCATE(dofConstraint%dofs(numberOfDofs),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate constraint DOFs array.",err,error,*999)
    ALLOCATE(dofConstraint%coefficients(numberOfDofs),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate constraint coefficients array.",err,error,*999)
    dofConstraint%globalDof=globalDof
    dofConstraint%numberOfDofs=numberOfDofs
    dofConstraint%dofs(1:numberOfDofs)=dofs(1:numberOfDofs)
    dofConstraint%coefficients(1:numberOfDofs)=coefficients(1:numberOfDofs)

    !Add new DOF constraint to new array
    newConstraints(dofConstraints%numberOfConstraints+1)%ptr=>dofConstraint
    !Replace old DOF constraints with new ones
    CALL MOVE_ALLOC(newConstraints,dofConstraints%constraints)
    dofConstraints%numberOfConstraints=dofConstraints%numberOfConstraints+1

    !Set the DOF type and BC type of the constrained DOF
    CALL BoundaryConditions_SetConditionType(boundaryConditionsVariable,globalDof,BOUNDARY_CONDITION_LINEAR_CONSTRAINT, &
      & err,error,*999)

    EXITS("BoundaryConditions_DofConstraintSet")
    RETURN
999 IF(ASSOCIATED(dofConstraint)) THEN
      IF(ALLOCATED(dofConstraint%dofs)) DEALLOCATE(dofConstraint%dofs)
      IF(ALLOCATED(dofConstraint%coefficients)) DEALLOCATE(dofConstraint%coefficients)
      DEALLOCATE(dofConstraint)
    END IF
    IF(ALLOCATED(newConstraints)) DEALLOCATE(newConstraints)
998 ERRORSEXITS("BoundaryConditions_DofConstraintSet",err,error)
    RETURN 1
  END SUBROUTINE BoundaryConditions_DofConstraintSet

  !
  !================================================================================================================================
  !

  !>Finish the creation of the dof constraints for a boundary conditions variable
  SUBROUTINE BoundaryConditions_DofConstraintsCreateFinish(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable !<A pointer to boundary conditions variable to finish the dof constraints for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: constraintIdx,dofIdx,thisDofDomain,otherDofDomain
    INTEGER(INTG) :: globalDof,globalDof2,localDof,localDof2
    INTEGER(INTG) :: numberOfCoupledDofs,numberOfGroupNodes
    INTEGER(INTG), ALLOCATABLE :: newCoupledGlobalDofs(:),newCoupledLocalDofs(:)
    REAL(DP), ALLOCATABLE :: newCoefficients(:)
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints
    TYPE(BoundaryConditionsDofConstraintType), POINTER :: dofConstraint
    TYPE(BoundaryConditionsCoupledDofsType), POINTER :: dofCoupling
    TYPE(DomainMappingType), POINTER :: variableDomainMapping
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(WorkGroupType), POINTER :: workGroup

    ENTERS("BoundaryConditions_DofConstraintsCreateFinish",err,error,*998)

    NULLIFY(dofCoupling)

    !We have a list of DOF constraints, which give the values for a field variable
    !DOF as a linear combination of other DOFs.
    !In order to be able to construct the solver matrices in the solver mapping routines,
    !we create a set of couplings, where a coupling is a set of field variable DOFs
    !mapped to a single solver row or column.

    IF(ASSOCIATED(boundaryConditionsVariable)) THEN
      fieldVariable=>boundaryConditionsVariable%variable
      IF(ASSOCIATED(fieldVariable)) THEN
        IF(ASSOCIATED(boundaryConditionsVariable%dofConstraints)) THEN
          dofConstraints=>boundaryConditionsVariable%dofConstraints
        ELSE
          CALL FlagError("Boundary conditions DOF constraints are not associated.",err,error,*998)
        END IF

        NULLIFY(variableDomainMapping)
        CALL FieldVariable_DomainMappingGet(fieldVariable,variableDomainMapping,err,error,*999)
        NULLIFY(workGroup)
        CALL DomainMapping_WorkGroupGet(variableDomainMapping,workGroup,err,error,*999)
        CALL WorkGroup_NumberOfGroupNodesGet(workGroup,numberOfGroupNodes,err,error,*999)

        !Allocate an array of pointers to DOF couplings
        IF(dofConstraints%numberOfConstraints>0) THEN
          ALLOCATE(dofConstraints%dofCouplings(fieldVariable%numberOfGlobalDofs),stat=err)
          IF(err/=0) CALL FlagError( &
            & "Could not allocate dof constraints dof couplings array.",err,error,*998)
          dofConstraints%numberOfDofs=fieldVariable%numberOfGlobalDofs
          DO dofIdx=1,fieldVariable%numberOfGlobalDofs
            NULLIFY(dofConstraints%dofCouplings(dofIdx)%ptr)
          END DO
        END IF

        !Loop over all constraints
        DO constraintIdx=1,dofConstraints%numberOfConstraints
          dofConstraint=>dofConstraints%constraints(constraintIdx)%ptr
          IF(.NOT.ASSOCIATED(dofConstraint)) THEN
            CALL FlagError("DOF constraint number "// &
              & TRIM(NumberToVstring(constraintIdx,"*",err,error))// &
              & " is not associated.",err,error,*999)
          END IF

          globalDof=dofConstraint%globalDof
          localDof=variableDomainMapping%globalToLocalMap(globalDof)%localNumber(1)
          thisDofDomain=variableDomainMapping%globalToLocalMap(globalDof)%domainNumber(1)

          !Check that the constrained DOFs are still set to be constrained, as
          !subsequently setting a boundary condition would change the DOF type but
          !not update the DOF constraints structure.
          IF(boundaryConditionsVariable%DOFTypes(globalDof)/=BOUNDARY_CONDITION_DOF_CONSTRAINED) THEN
            CALL FlagError("Global DOF number "//TRIM(NumberToVstring(globalDof,"*",err,error))// &
              & " is part of a linear constraint but the DOF type has been changed"// &
              & " by applying a boundary condition.",err,error,*999)
          END IF

          DO dofIdx=1,dofConstraint%numberOfDofs
            globalDof2=dofConstraint%dofs(dofIdx)
            localDof2=variableDomainMapping%globalToLocalMap(globalDof2)%localNumber(1)
            !Check a Dirichlet conditions hasn't also been set on this DOF
            IF(boundaryConditionsVariable%DOFTypes(globalDof2)/=BOUNDARY_CONDITION_DOF_FREE) THEN
              CALL FlagError("A Dirichlet boundary condition has been set on DOF number "// &
                & TRIM(NumberToVstring(globalDof2,"*",err,error))// &
                & " which is part of a linear constraint.",err,error,*999)
            END IF

            !Check we don't have DOF constraints that are split over domains
            !\todo Implement support for DOF constraints that are split over domains
            IF(numberOfGroupNodes>1) THEN
              otherDofDomain=variableDomainMapping%globalToLocalMap(globalDof2)%domainNumber(1)
              IF(thisDofDomain/=otherDofDomain) THEN
                CALL FlagError("An equal DOF constraint is split over multiple domains, "// &
                  & "support for this has not yet been implemented.",err,error,*999)
              END IF
            END IF

            !Add to the DOFs that are coupled with globalDof2
            !This might be quite inefficient if there are many dofs mapped to a single row/column
            !due to the reallocation at each step
            IF(ASSOCIATED(dofConstraints%dofCouplings(globalDof2)%ptr)) THEN
              numberOfCoupledDofs=dofConstraints%dofCouplings(globalDof2)%ptr%numberOfDofs
              ALLOCATE(newCoupledGlobalDofs(numberOfCoupledDofs+1),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate new DOF coupling global DOFs.",err,error,*999)
              ALLOCATE(newCoupledLocalDofs(numberOfCoupledDofs+1),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate new DOF coupling local DOFs.",err,error,*999)
              ALLOCATE(newCoefficients(numberOfCoupledDofs+1),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate new DOF coupling values.",err,error,*999)
              newCoupledGlobalDofs(1:numberOfCoupledDofs)=dofCoupling%globalDofs(1:numberOfCoupledDofs)
              newCoupledLocalDofs(1:numberOfCoupledDofs)=dofCoupling%localDofs(1:numberOfCoupledDofs)
              newCoefficients(1:numberOfCoupledDofs)=dofCoupling%coefficients(1:numberOfCoupledDofs)
            ELSE
              !Set up a a new dofCoupling and set globalDof2 as the first DOF
              ALLOCATE(dofConstraints%dofCouplings(globalDof2)%ptr,stat=err)
              IF(err/=0) CALL FlagError("Could not allocate new DOF coupling type.",err,error,*999)
              ALLOCATE(newCoupledGlobalDofs(2),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate new DOF coupling global DOFs.",err,error,*999)
              ALLOCATE(newCoupledLocalDofs(2),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate new DOF coupling local DOFs.",err,error,*999)
              ALLOCATE(newCoefficients(2),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate new DOF coupling values.",err,error,*999)
              newCoupledGlobalDofs(1)=globalDof2
              newCoupledLocalDofs(1)=localDof2
              newCoefficients(1)=1.0_DP
              numberOfCoupledDofs=1
            END IF
            dofCoupling=>dofConstraints%dofCouplings(globalDof2)%ptr
            newCoupledGlobalDofs(numberOfCoupledDofs+1)=globalDof
            newCoupledLocalDofs(numberOfCoupledDofs+1)=localDof
            newCoefficients(numberOfCoupledDofs+1)=dofConstraint%coefficients(dofIdx)
            CALL MOVE_ALLOC(newCoupledGlobalDofs,dofCoupling%globalDofs)
            CALL MOVE_ALLOC(newCoupledLocalDofs,dofCoupling%localDofs)
            CALL MOVE_ALLOC(newCoefficients,dofCoupling%coefficients)
            dofCoupling%numberOfDofs=numberOfCoupledDofs+1
          END DO !dofIdx
        END DO !constraintIdx
      ELSE
        CALL FlagError("Field variable is not associated for this boundary conditions variable",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
    END IF

    EXITS("BoundaryConditions_DofConstraintsCreateFinish")
    RETURN
999 IF(ALLOCATED(newCoupledGlobalDofs)) DEALLOCATE(newCoupledGlobalDofs)
    IF(ALLOCATED(newCoupledLocalDofs)) DEALLOCATE(newCoupledLocalDofs)
    IF(ALLOCATED(newCoefficients)) DEALLOCATE(newCoefficients)
    CALL BoundaryConditions_DofConstraintsFinalise(dofConstraints,err,error,*998)
998 ERRORS("BoundaryConditions_DofConstraintsCreateFinish",err,error)
    EXITS("BoundaryConditions_DofConstraintsCreateFinish")
    RETURN 1

  END SUBROUTINE BoundaryConditions_DofConstraintsCreateFinish

  !
  !================================================================================================================================
  !

  !>Finalise the DOF constraints structure
  SUBROUTINE BoundaryConditions_DofConstraintsFinalise(dofConstraints,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints !<A pointer to the dof constraints to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: constraintIdx,dofIdx

    ENTERS("BoundaryConditions_DofConstraintsFinalise",err,error,*999)

    IF(ASSOCIATED(dofConstraints)) THEN
      IF(ALLOCATED(dofConstraints%constraints)) THEN
        DO constraintIdx=1,dofConstraints%numberOfConstraints
          IF(ASSOCIATED(dofConstraints%constraints(constraintIdx)%ptr)) THEN
            IF(ALLOCATED(dofConstraints%constraints(constraintIdx)%ptr%dofs)) THEN
              DEALLOCATE(dofConstraints%constraints(constraintIdx)%ptr%dofs)
            END IF
            IF(ALLOCATED(dofConstraints%constraints(constraintIdx)%ptr%coefficients)) THEN
              DEALLOCATE(dofConstraints%constraints(constraintIdx)%ptr%coefficients)
            END IF
            DEALLOCATE(dofConstraints%constraints(constraintIdx)%ptr)
          END IF
        END DO
        DEALLOCATE(dofConstraints%constraints)
      END IF
      IF(ALLOCATED(dofConstraints%dofCouplings)) THEN
        DO dofIdx=1,dofConstraints%numberOfDofs
          IF(ASSOCIATED(dofConstraints%dofCouplings(dofIdx)%ptr)) THEN
            IF(ALLOCATED(dofConstraints%dofCouplings(dofIdx)%ptr%globalDofs)) THEN
              DEALLOCATE(dofConstraints%dofCouplings(dofIdx)%ptr%globalDofs)
            END IF
            IF(ALLOCATED(dofConstraints%dofCouplings(dofIdx)%ptr%localDofs)) THEN
              DEALLOCATE(dofConstraints%dofCouplings(dofIdx)%ptr%localDofs)
            END IF
            IF(ALLOCATED(dofConstraints%dofCouplings(dofIdx)%ptr%coefficients)) THEN
              DEALLOCATE(dofConstraints%dofCouplings(dofIdx)%ptr%coefficients)
            END IF
          END IF
        END DO
        DEALLOCATE(dofConstraints%dofCouplings)
      END IF
    ELSE
      CALL FlagError("dofConstraints pointer is not associated.",err,error,*999)
    END IF

    EXITS("BoundaryConditions_DofConstraintsFinalise")
    RETURN
999 ERRORSEXITS("BoundaryConditions_DofConstraintsFinalise",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditions_DofConstraintsFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the DOF constraints structure
  SUBROUTINE BoundaryConditions_DofConstraintsInitialise(dofConstraints,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints !<A pointer to the dof constraints to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("BoundaryConditions_DofConstraintsInitialise",err,error,*999)

    IF(ASSOCIATED(dofConstraints)) THEN
      dofConstraints%numberOfConstraints=0
      dofConstraints%numberOfDofs=0
    ELSE
      CALL FlagError("dofConstraints pointer is not associated.",err,error,*999)
    END IF

    EXITS("BoundaryConditions_DofConstraintsInitialise")
    RETURN
999 ERRORSEXITS("BoundaryConditions_DofConstraintsInitialise",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditions_DofConstraintsInitialise

  !
  !================================================================================================================================
  !

  !>Finalise the boundary conditions variable and deallocate all memory.
  SUBROUTINE BoundaryCondition_VariableFinalise(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(BoundaryConditionsDirichletType), POINTER :: boundaryConditionsDirichlet

    ENTERS("BoundaryCondition_VariableFinalise",err,error,*999)

    IF(ASSOCIATED(boundaryConditionsVariable)) THEN
      IF(ALLOCATED(boundaryConditionsVariable%conditionTypes))  &
        & DEALLOCATE(boundaryConditionsVariable%conditionTypes)
      IF(ALLOCATED(boundaryConditionsVariable%DOFTypes))  &
        & DEALLOCATE(boundaryConditionsVariable%DOFTypes)
      IF(ASSOCIATED(boundaryConditionsVariable%dirichletBoundaryConditions)) THEN
        boundaryConditionsDirichlet=>boundaryConditionsVariable%dirichletBoundaryConditions
        CALL BoundaryConditions_SparsityIndicesArrayFinalise(boundaryConditionsDirichlet% &
            & linearSparsityIndices,err,error,*999)
        CALL BoundaryConditions_SparsityIndicesArrayFinalise(boundaryConditionsDirichlet% &
            & dynamicSparsityIndices,err,error,*999)
        IF(ALLOCATED(boundaryConditionsDirichlet%dirichletDOFIndices)) THEN
          DEALLOCATE(boundaryConditionsDirichlet%dirichletDOFIndices)
        ENDIF
        DEALLOCATE(boundaryConditionsDirichlet)
      ENDIF
      CALL BoundaryConditions_NeumannFinalise(boundaryConditionsVariable,err,error,*999)
      IF(ASSOCIATED(boundaryConditionsVariable%pressureIncrementedBoundaryConditions)) &
        & DEALLOCATE(boundaryConditionsVariable%pressureIncrementedBoundaryConditions)
      IF(ASSOCIATED(boundary_conditions_variable%dofConstraints)) THEN
        CALL BoundaryConditions_DofConstraintsFinalise(boundary_conditions_variable%dofConstraints,err,error,*999)
        DEALLOCATE(boundary_conditions_variable%dofConstraints)
      END IF
      DEALLOCATE(boundaryConditionsVariable)
    ENDIF

    EXITS("BoundaryCondition_VariableFinalise")
    RETURN
999 ERRORSEXITS("BoundaryCondition_VariableFinalise",err,error)
    RETURN 1
  END SUBROUTINE BoundaryCondition_VariableFinalise

  !
  !================================================================================================================================
  !

  !>Finalise an array of sparcity indices and deallocate all memory.
  SUBROUTINE BoundaryConditions_SparsityIndicesArrayFinalise(sparsityIndicesArray,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsSparsityIndicesPtrType), ALLOCATABLE :: sparsityIndicesArray(:,:)
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx, equationsMatrixIdx
    TYPE(BoundaryConditionsSparsityIndicesType), POINTER :: sparsityIndices

    ENTERS("BoundaryConditions_SparsityIndicesArrayFinalise",err,error,*999)

    IF (ALLOCATED(sparsityIndicesArray)) THEN
      DO equationsSetIdx=1,SIZE(sparsityIndicesArray,1)
        DO equationsMatrixIdx=1,SIZE(sparsityIndicesArray,2)
          sparsityIndices=>sparsityIndicesArray(equationsSetIdx,equationsMatrixIdx)%PTR
          IF(ASSOCIATED(sparsityIndices)) THEN
            IF(ALLOCATED(sparsityIndices%sparseRowIndices)) THEN
              DEALLOCATE(sparsityIndices%sparseRowIndices)
            ENDIF
            IF(ALLOCATED(sparsityIndices%sparseColumnIndices)) THEN
              DEALLOCATE(sparsityIndices%sparseColumnIndices)
            ENDIF
            DEALLOCATE(sparsityIndices)
          ENDIF
        ENDDO
      ENDDO
      DEALLOCATE(sparsityIndices_ARRAY)
    ENDIF

    EXITS("BoundaryConditions_SparsityIndicesArrayFinalise")
    RETURN
999 ERRORS("BoundaryConditions_SparsityIndicesArrayFinalise",err,error)
    EXITS("BoundaryConditions_SparsityIndicesArrayFinalise")
    RETURN 1

  END SUBROUTINE BoundaryConditions_SparsityIndicesArrayFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the boundary conditions variable for a variable type if that variable has not already been initialised, otherwise do nothing.
  SUBROUTINE BoundaryConditions_VariableInitialise(boundaryConditions,fieldVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to initialise a variable type for.
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to initialise the boundary conditions variable for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,variableIdx
    TYPE(DomainMappingType), POINTER :: domainMapping
    TYPE(VARYING_STRING) :: dummyError
    TYPE(BoundaryConditionsVariablePtrType), ALLOCATABLE :: newBoundaryConditionsVariables(:)
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable

    ENTERS("BoundaryConditions_VariableInitialise",err,error,*998)

    IF(ASSOCIATED(boundaryConditions)) THEN
      IF(ASSOCIATED(fieldVariable)) THEN
        domainMapping=>fieldVariable%domainMapping
        IF(ASSOCIATED(domainMapping)) THEN
          !Check if boundary conditions variable has already been added, if so then we don't do anything as different equations
          !sets can have the same dependent field variables and will both want to add the variable
          NULLIFY(boundaryConditionsVariable)
          CALL BoundaryConditions_VariableExists(boundaryConditions,fieldVariable,boundaryConditionsVariable,err,error,*999)
          IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) THEN
            ALLOCATE(newBoundaryConditionsVariables(boundaryConditions%numberOfBoundaryConditionsVariables+1),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate new boundary conditions variables array.",err,error,*998)
            IF(ALLOCATED(boundaryConditions%boundaryConditionsVariables)) THEN
              DO variableIdx=1,boundaryConditions%numberOfBoundaryConditionsVariables
                newBoundaryConditionsVariables(variableIdx)%PTR=> &
                    & boundaryConditions%boundaryConditionsVariables(variableIdx)%PTR
              ENDDO
            ENDIF

            ALLOCATE(newBoundaryConditionsVariables(boundaryConditions%numberOfBoundaryConditionsVariables+1)%PTR,STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate boundary condition variable.",err,error,*998)
            boundaryConditionsVariable=>newBoundaryConditionsVariables( &
                & boundaryConditions%numberOfBoundaryConditionsVariables+1)%PTR
            boundaryConditionsVariable%boundaryConditions=>boundaryConditions
            boundaryConditionsVariable%variableType=fieldVariable%variableType
            boundaryConditionsVariable%VARIABLE=>fieldVariable
            ALLOCATE(boundaryConditionsVariable%conditionTypes(domainMapping%numberOfGlobal),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate global boundary condition types.",err,error,*999)
            ALLOCATE(boundaryConditionsVariable%DOFTypes(domainMapping%numberOfGlobal),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate global boundary condition dof types.",err,error,*999)
            boundaryConditionsVariable%conditionTypes=BOUNDARY_CONDITION_FREE
            boundaryConditionsVariable%DOFTypes=BOUNDARY_CONDITION_DOF_FREE
            ALLOCATE(boundaryConditionsVariable%dofCounts(MAX_BOUNDARY_CONDITION_NUMBER),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate boundary condition DOF counts array.",err,error,*999)
            boundaryConditionsVariable%dofCounts=0
            NULLIFY(boundaryConditionsVariable%dirichletBoundaryConditions)
            boundaryConditionsVariable%numberOfDirichletConditions=0
            NULLIFY(boundaryConditionsVariable%neumannBoundaryConditions)
            NULLIFY(boundaryConditionsVariable%pressureIncrementedBoundaryConditions)
            ALLOCATE(boundaryConditionsVariable%parameterSetRequired(FIELD_NUMBER_OF_SET_TYPES),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate boundary condition parameter set required array.",err,error,*999)
            boundaryConditionsVariable%parameterSetRequired=.FALSE.
            boundaryConditionsVariable%parameterSetRequired(FIELD_VALUES_SET_TYPE)=.TRUE.

            CALL MOVE_ALLOC(newBoundaryConditionsVariables,boundaryConditions%boundaryConditionsVariables)
            boundaryConditions%numberOfBoundaryConditionsVariables= &
                & boundaryConditions%numberOfBoundaryConditionsVariables+1

            ALLOCATE(boundary_conditions_variable%DofConstraints,stat=err)
            IF(err/=0) CALL FlagError("Could not allocate boundary conditions dof constraints.",err,error,*999)
            CALL BoundaryConditions_DofConstraintsInitialise(boundary_conditions_variable%DofConstraints,err,error,*999)

          END IF
        ELSE
          CALL FlagError("Field variable domain mapping is not associated.",err,error,*998)
        ENDIF
      ELSE
        CALL FlagError("Field variable is not associated.",err,error,*998)
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions is not associated.",err,error,*998)
    ENDIF

    EXITS("BoundaryConditions_VariableInitialise")
    RETURN
999 CALL BoundaryCondition_VariableFinalise(boundaryConditionsVariable,dummyErr,dummyError,*998)
    DEALLOCATE(newBoundaryConditionsVariables)
998 ERRORSEXITS("BoundaryConditions_VariableInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_VariableInitialise

  !
  !================================================================================================================================
  !

  !>Find the boundary conditions variable for a given field variable if it exists
  SUBROUTINE BoundaryConditions_VariableExists(boundaryConditions,fieldVariable,boundaryConditionsVariable,err,error,*999)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to see if the boundary conditions variable exists for
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to check the boundary conditions variable for.
    TYPE(BoundaryConditionVariableType), POINTER, INTENT(OUT) :: boundaryConditionsVariable !<On return, a pointer to the boundary conditions variable, or NULL if it wasn't found. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: variableIdx
    TYPE(FieldVariableType), POINTER :: variable
    LOGICAL :: variableFound
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_VariableExists",err,error,*999)

    IF(ASSOCIATED(boundaryConditionsVariable)) CALL FlagError("Boundary conditions variable is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    
    IF(ALLOCATED(boundaryCondtions%boundaryConditionsVariables)) THEN    
      variableFound=.FALSE.
      variableIdx=1
      DO WHILE(variableIdx<=boundaryConditions%numberOfBoundaryConditionsVariables.AND..NOT.variableFound)
        IF(.NOT.ASSOCIATED(boundaryConditions%boundaryConditionsVariables(variableIdx)%ptr)) THEN
          localError="Boundary conditions variable pointer is not associated for variable index "// &
            & TRIM(NumberToVString(variableIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        variable=>boundaryConditions%boundaryConditionsVariables(variableIdx)%ptr%variable
        IF(.NOT.ASSOCIATED(variable)) THEN
          localError="Boundary conditions variable pointer variable is not associated for variable index "// &
            & TRIM(NumberToVString(variableIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        NULLIFY(variableField)
        CALL FieldVariable_FieldGet(variable,variableField,err,error,*999)
        NULLIFY(fieldVariableField)
        CALL FieldVariable_FieldGet(fieldVariable,fieldVariableField,err,error,*999)
        IF(variable%variableType==fieldVariable%variableType.AND. &
          & variableField%userNumber==fieldVaribleField%userNumber) THEN
          IF(ASSOCIATED(variableField%region)) THEN
            IF(ASSOCIATED(fieldVariableField%region)) THEN
              IF(variableField%region%userNumber==fieldVariableField%region%userNumber) THEN
                variableFound=.TRUE.
                boundaryConditionsVariable=>boundaryConditions%boundaryConditionsVariables(variableIdx)%ptr
              ENDIF
            ENDIF
          ELSE IF(ASSOCIATED(variableField%INTERFACE)) THEN
            !!TODO: should also check parent regions
            IF(ASSOCIATED(fieldVariableField%INTERFACE)) THEN
              IF(variableField%INTERFACE%userNumber==fieldVariableField%INTERFACE%userNumber) THEN
                VARIABLE_FOUND=.TRUE.
                boundaryConditionsVariable=>boundaryConditions%boundaryConditionsVariables(variableIdx)%PTR
              ENDIF
            ENDIF
          ELSE
            localError="Field number "//TRIM(NumberToVString(variableField%userNumber,"*",err,error))// &
              & " does not have a region or interface associated.",err,error,*999)
          ENDIF
        ENDIF
        variableIdx=variableIdx+1
      ENDIF
    ENDDO


    EXITS("BoundaryConditions_VariableExists")
    RETURN
999 ERRORSEXITS("BoundaryConditions_VariableExists",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_VariableExists

  !
  !================================================================================================================================
  !

  !>Find the boundary conditions variable for a given field variable
  SUBROUTINE BoundaryConditions_VariableGet(boundaryConditions,fieldVariable,boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to get the boundary conditions variable for.
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the boundary conditions variable for.
    TYPE(BoundaryConditionVariableType), POINTER, INTENT(OUT) :: boundaryConditionsVariable !<On return, a pointer to the boundary conditions variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_VariableGet",err,error,*999)

    CALL BoundaryConditions_VariableExists(boundaryConditions,fieldVariable,boundaryConditionsVariable,err,error,*999)

    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) THEN
      localError="Boundary condition variable is not associated"
      IF(ASSOCIATED(fieldVariable)) THEN
        localError=localError//" for variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
        IF(ASSOCIATED(fieldVariable%field)) THEN
          localError=localError//" of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
          IF(ASSOCIATED(fieldVariable%field%region)) THEN
            localError=localError//" of region number "//TRIM(NumberToVString(fieldVariable%field%region%userNumber,"*",err,error))
          ELSE IF(ASSOCIATED(fieldVariable%field%interface)) THEN
            localError=localError//" of interface number "// &
              & TRIM(NumberToVString(fieldVariable%field%interface%userNumber,"*",err,error))
            IF(ASSOCIATED(fieldVariable%field%interface%parentRegion)) &
              & localError=localError//" of parent region number "// &
              & TRIM(NumberToVString(fieldVariable%field%interface%parentRegion%userNumber,"*",err,error))
          ENDIF
        ENDIF
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("BoundaryConditions_VariableGet")
    RETURN
999 ERRORSEXITS("BoundaryConditions_VariableGet",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_VariableGet

  !
  !================================================================================================================================
  !

  !>Initialise dirichlet boundary conditions for a boundary conditions.
  SUBROUTINE BoundaryConditions_DirichletInitialise(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to initialise a boundary conditions dirichlet type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfDirichletConditions,numberOfLinearMatrices,numberOfDynamicMatrices,matrixIdx, &
      & maxNumberOfLinearMatrices,maxNumberOfDynamicMatrices,equationsSetIdx
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(BoundaryConditionsDirichletType), POINTER :: boundaryConditionsDirichlet
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping

    ENTERS("BoundaryConditions_DirichletInitialise",err,error,*999)

    IF(ASSOCIATED(boundaryConditionsVariable)) THEN
      IF(ASSOCIATED(boundaryConditionsVariable%dirichletBoundaryConditions)) THEN
        CALL FlagError("Dirichlet boundary conditions are already associated for this boundary conditions variable." &
           & ,err,error,*999)
      ELSE
        ALLOCATE(boundaryConditionsVariable%dirichletBoundaryConditions,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate Dirichlet Boundary Conditions",err,error,*999)
        boundaryConditionsDirichlet=>boundaryConditionsVariable%dirichletBoundaryConditions
        numberOfDirichletConditions=boundaryConditionsVariable%numberOfDirichletConditions
        ALLOCATE(boundaryConditionsDirichlet%dirichletDOFIndices(numberOfDirichletConditions),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate Dirichlet DOF indices array",err,error,*999)

        solverEquations=>boundaryConditionsVariable%boundaryConditions%solverEquations
        IF(ASSOCIATED(solverEquations)) THEN
          maxNumberOfLinearMatrices=0
          maxNumberOfDynamicMatrices=0
          DO equationsSetIdx=1,solverEquations%solverMapping%numberOfEquationsSets
            equationsSet=>solverEquations%solverMapping%equationsSets(equationsSetIdx)%PTR
            IF(ASSOCIATED(equationsSet)) THEN
              NULLIFY(equations)
              CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
              NULLIFY(vectorEquations)
              CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
              NULLIFY(vectorMapping)
              CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
              linearMapping=>vectorMapping%linearMapping
              dynamicMapping=>vectorMapping%dynamicMapping
              IF(ASSOCIATED(linearMapping)) THEN
                numberOfLinearMatrices=linearMapping%numberOfLinearMatrices
                IF(numberOfLinearMatrices>maxNumberOfLinearMatrices) &
                  & maxNumberOfLinearMatrices=numberOfLinearMatrices
              ENDIF
              IF(ASSOCIATED(dynamicMapping)) THEN
                numberOfDynamicMatrices=dynamicMapping%numberOfDynamicMatrices
                IF(numberOfDynamicMatrices>maxNumberOfDynamicMatrices) &
                  & maxNumberOfDynamicMatrices=numberOfDynamicMatrices
              ENDIF
            ELSE
              CALL FlagError("Equations set is not associated.",err,error,*999)
            ENDIF
          ENDDO
          ALLOCATE(boundaryConditionsDirichlet%linearSparsityIndices(solverEquations%solverMapping%numberOfEquationsSets, &
                & maxNumberOfLinearMatrices),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate Dirichlet linear sparsity indices array",err,error,*999)
          ALLOCATE(boundaryConditionsDirichlet%dynamicSparsityIndices(solverEquations%solverMapping%numberOfEquationsSets,&
                & maxNumberOfDynamicMatrices),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate Dirichlet dynamic sparsity indices array",err,error,*999)
          DO equationsSetIdx=1,solverEquations%solverMapping%numberOfEquationsSets
            DO matrixIdx=1,maxNumberOfLinearMatrices
              NULLIFY(boundaryConditionsDirichlet%linearSparsityIndices(equationsSetIdx,matrixIdx)%PTR)
            ENDDO
            DO matrixIdx=1,maxNumberOfDynamicMatrices
              NULLIFY(boundaryConditionsDirichlet%dynamicSparsityIndices(equationsSetIdx,matrixIdx)%PTR)
            ENDDO
          ENDDO
        ELSE
          CALL FlagError("Solver equations is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
    ENDIF

    EXITS("BoundaryConditions_DirichletInitialise")
    RETURN
!!TODO \todo write boundaryConditionsDirichlet_FINALISE
999 ERRORSEXITS("BoundaryConditions_DirichletInitialise",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditions_DirichletInitialise

  !
  !================================================================================================================================
  !

  !>Initialise Sparsity Indices type
  SUBROUTINE BoundaryConditions_SparsityIndicesInitialise(sparsityIndices,numberOfDirichlet,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsSparsityIndicesType), POINTER :: sparsityIndices !<A pointer to the Sparsity Indices type tp initialise
    INTEGER(INTG), INTENT(IN) :: numberOfDirichlet !<The number of dirichlet conditions this sparsity indices type will hold
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("BoundaryConditions_SparsityIndicesInitialise",err,error,*999)

    IF(ASSOCIATED(sparsityIndices)) THEN
     CALL FlagError("Sparsity Indices are already associated.",err,error,*999)
    ELSE
      ALLOCATE(sparsityIndices,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate sparsity indicies.",err,error,*999)
      ALLOCATE(sparsityIndices%sparseColumnIndices(numberOfDirichlet+1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate sparsity column indices array",err,error,*999)
    ENDIF

    EXITS("BoundaryConditions_SparsityIndicesInitialise")
    RETURN
!!TODO \todo write boundaryConditions_SPARSITY_INDICES_FINALISE
999 ERRORS("BoundaryConditions_SparsityIndicesInitialise",err,error)
    EXITS("BoundaryConditions_SparsityIndicesInitialise")
    RETURN 1

  END SUBROUTINE BoundaryConditions_SparsityIndicesInitialise

  !
  !================================================================================================================================
  !

  !>Initialises the pressure incremented boundary condition.
  SUBROUTINE BoundaryConditions_PressureIncrementedInitialise(boundaryConditionsVariable,err,error,*)
    !Argument variables
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to initialise a boundary conditions dirichlet type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(BoundaryConditionsPressureIncrementedType), POINTER :: boundaryConditionsPressureIncremented
    INTEGER(INTG) :: numberOfPressureIncrementedConditions

    ENTERS("BoundaryConditions_PressureIncrementedInitialise",err,error,*999)

    IF(ASSOCIATED(boundaryConditionsVariable)) THEN
      IF(ASSOCIATED(boundaryConditionsVariable%pressureIncrementedBoundaryConditions)) THEN
        CALL FlagError("Pressure incremented boundary conditions are already associated for this boundary conditions variable." &
           & ,err,error,*999)
      ELSE
        ALLOCATE(boundaryConditionsVariable%pressureIncrementedBoundaryConditions,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate Pressure incremented Boundary Conditions",err,error,*999)
        boundaryConditionsPressureIncremented=>boundaryConditionsVariable%pressureIncrementedBoundaryConditions
        numberOfPressureIncrementedConditions=boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
        ALLOCATE(boundaryConditionsPressureIncremented%pressureIncrementedDOFIndices &
          & (numberOfPressureIncrementedConditions),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate Pressure incremented DOF indices array",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
    ENDIF

    EXITS("BoundaryConditions_PressureIncrementedInitialise")
    RETURN
!!TODO \todo write boundaryConditionsPressureIncremented_FINALISE
999 ERRORS("BoundaryConditions_PressureIncrementedInitialise",err,error)
    EXITS("BoundaryConditions_PressureIncrementedInitialise")
    RETURN 1

  END SUBROUTINE BoundaryConditions_PressureIncrementedInitialise

  !
  !================================================================================================================================
  !

END MODULE BoundaryConditionsRoutines
