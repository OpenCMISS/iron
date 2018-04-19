!> \file
!> \author Chris Bradley
!> \brief This module handles all equations matrix and rhs routines.
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
!> Contributor(s): David Ladd
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

!> This module handles all equations matrix and rhs routines.
MODULE EquationsMatricesRoutines

  USE BaseRoutines
  USE DistributedMatrixVector
  USE EquationsAccessRoutines
  USE EquationsMappingAccessRoutines
  USE EquationsMatricesAccessRoutines
  USE EquationsSetAccessRoutines
  USE EquationsSetConstants
  USE FIELD_ROUTINES
  USE ISO_VARYING_STRING
  USE Kinds
  USE Lists
  USE MatrixVector
  USE Strings
  USE Types
  USE LINKEDLIST_ROUTINES

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup EquationsMatricesRoutines_EquationsMatrixStructureTypes EquationsMatricesRoutines::EquationsMatrixStructureTypes
  !> \brief Equations matrices structure (sparsity) types
  !> \see EquationsMatricesRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRIX_NO_STRUCTURE=1 !<No matrix structure - all elements can contain a value. \see EquationsMatricesRoutines_EquationsMatrixStructureTypes,EquationsMatricesRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRIX_FEM_STRUCTURE=2 !<Finite element matrix structure. \see EquationsMatricesRoutines_EquationsMatrixStructureTypes,EquationsMatricesRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRIX_DIAGONAL_STRUCTURE=3 !<Diagonal matrix structure. \see EquationsMatricesRoutines_EquationsMatrixStructureTypes,EquationsMatricesRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRIX_NODAL_STRUCTURE=4 !<Nodal matrix structure. \see EquationsMatricesRoutines_EquationsMatrixStructureTypes,EquationsMatricesRoutines
  !>@}


  !> \addtogroup EquationsMatricesRoutines_LumpingTypes EquationsMatricesRoutines::LumpingTypes
  !> \brief Equations matrix lumping types
  !> \see EquationsMatricesRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRIX_UNLUMPED=1 !<The matrix is not lumped \see EquationsMatricesRoutines_LumpingTypes,EquationsMatricesRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRIX_LUMPED=2 !<The matrix is "mass" lumped \see EquationsMatricesRoutines_LumpingTypes,EquationsMatricesRoutines
 !>@}
 
  !> \addtogroup EquationsMatricesRoutines_EquationsMatricesSparsityTypes EquationsMatricesRoutines::EquationsMatricesSparsityTypes
  !> \brief Equations matrices sparsity types
  !> \see EquationsMatricesRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_SPARSE_MATRICES=1 !<Use sparse equations matrices \see EquationsMatricesRoutines_EquationsMatricesSparsityTypes,EquationsMatricesRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_FULL_MATRICES=2 !<Use fully populated equation matrices \see EquationsMatricesRoutines_EquationsMatricesSparsityTypes,EquationsMatricesRoutines
  !>@}

  !> \addtogroup EquationsMatricesRoutines_SelectMatricesTypes EquationsMatricesRoutines::SelectMatricesTypes
  !> \brief The types of selection available for the equations matrices
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_ALL=1 !<Select all the equations matrices and vectors \see EquationsMatricesRoutines::SelectMatricesTypes,EquationsMatricesRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_DYNAMIC_ONLY=2 !<Select only the dynamic equations matrices and vectors \see EquationsMatricesRoutines::SelectMatricesTypes,EquationsMatricesRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_LINEAR_ONLY=3 !<Select only the linear equations matrices and vectors \see EquationsMatricesRoutines::SelectMatricesTypes,EquationsMatricesRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_NONLINEAR_ONLY=4 !<Select only the nonlinear equations matrices and vectors \see EquationsMatricesRoutines::SelectMatricesTypes,EquationsMatricesRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_JACOBIAN_ONLY=5 !<Select only the Jacobian equations matrix \see EquationsMatricesRoutines::SelectMatricesTypes,EquationsMatricesRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_RESIDUAL_ONLY=6 !<Select only the residual equations vector \see EquationsMatricesRoutines::SelectMatricesTypes,EquationsMatricesRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_RHS_ONLY=7 !<Select only the RHS equations vector \see EquationsMatricesRoutines::SelectMatricesTypes,EquationsMatricesRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_SOURCE_ONLY=8 !<Select only the RHS equations vector \see EquationsMatricesRoutines::SelectMatricesTypes,EquationsMatricesRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_RHS_RESIDUAL_ONLY=9 !<Select only the RHS and residual equations vectors \see EquationsMatricesRoutines::SelectMatricesTypes,EquationsMatricesRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_RHS_SOURCE_ONLY=10 !<Assemble only the RHS and source equations vectors \see EquationsMatricesRoutines::SelectMatricesTypes,EquationsMatricesRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_RESIDUAL_SOURCE_ONLY=11 !<Assemble only the residual and source equations vectors\see EquationsMatricesRoutines::SelectMatricesTypes,EquationsMatricesRoutines 
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_VECTORS_ONLY=12 !<Assemble only the equations vectors \see EquationsMatricesRoutines::SelectMatricesTypes,EquationsMatricesRoutines
  !>@}

  !> \addtogroup EquationsMatricesRoutines_GradientCalculationTypes EquationsMatricesRoutines::GradientCalculationTypes
  !> \brief Gradient calculation types
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_GRADIENT_FINITE_DIFFERENCE_CALCULATED=1 !<Use finite differencing to calculate the gradient
  INTEGER(INTG), PARAMETER :: EQUATIONS_GRADIENT_ANALYTIC_CALCULATED=2 !<Use an analytic gradient evaluation
  !>@}
  !> \addtogroup EquationsMatricesRoutines_JacobianCalculationTypes EquationsMatricesRoutines::JacobianCalculationTypes
  !> \brief Jacobian calculation types
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_JACOBIAN_FINITE_DIFFERENCE_CALCULATED=1 !<Use finite differencing to calculate the Jacobian
  INTEGER(INTG), PARAMETER :: EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED=2 !<Use an analytic Jacobian evaluation
  !>@}
  !> \addtogroup EquationsMatricesRoutines_HessianCalculationTypes EquationsMatricesRoutines::HessianCalculationTypes
  !> \brief Hessian calculation types
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_HESSIAN_FINITE_DIFFERENCE_CALCULATED=1 !<Use finite differencing to calculate the Hessian
  INTEGER(INTG), PARAMETER :: EQUATIONS_HESSIAN_ANALYTIC_CALCULATED=2 !<Use an analytic Hessian evaluation
  !>@}
  !>@}

  !Module types

  !Module variables

  !Interfaces

  !>Sets the storage type (sparsity) of the nonlinear (Jacobian) equations matrices
  INTERFACE EquationsMatrices_NonlinearStorageTypeSet
    MODULE PROCEDURE EquationsMatrices_NonlinearStorageTypeSet0
    MODULE PROCEDURE EquationsMatrices_NonlinearStorageTypeSet1
  END INTERFACE EquationsMatrices_NonlinearStorageTypeSet

  !>Sets the structure (sparsity) of the nonlinear (Jacobian) equations matrices
  INTERFACE EquationsMatrices_NonlinearStructureTypeSet
    MODULE PROCEDURE EquationsMatrices_NonlinearStructureTypeSet0
    MODULE PROCEDURE EquationsMatrices_NonlinearStructureTypeSet1
  END INTERFACE EquationsMatrices_NonlinearStructureTypeSet

  PUBLIC EQUATIONS_MATRIX_NO_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE, &
    & EQUATIONS_MATRIX_NODAL_STRUCTURE

  PUBLIC EQUATIONS_MATRIX_UNLUMPED,EQUATIONS_MATRIX_LUMPED

  PUBLIC EQUATIONS_MATRICES_SPARSE_MATRICES,EQUATIONS_MATRICES_FULL_MATRICES

  PUBLIC EQUATIONS_MATRICES_ALL,EQUATIONS_MATRICES_LINEAR_ONLY,EQUATIONS_MATRICES_NONLINEAR_ONLY,EQUATIONS_MATRICES_JACOBIAN_ONLY, &
    & EQUATIONS_MATRICES_RESIDUAL_ONLY,EQUATIONS_MATRICES_RHS_ONLY,EQUATIONS_MATRICES_SOURCE_ONLY, &
    & EQUATIONS_MATRICES_RHS_RESIDUAL_ONLY,EQUATIONS_MATRICES_RHS_SOURCE_ONLY,EQUATIONS_MATRICES_RESIDUAL_SOURCE_ONLY, &
    & EQUATIONS_MATRICES_VECTORS_ONLY
  
  PUBLIC EQUATIONS_GRADIENT_FINITE_DIFFERENCE_CALCULATED,EQUATIONS_GRADIENT_ANALYTIC_CALCULATED

  PUBLIC EQUATIONS_JACOBIAN_FINITE_DIFFERENCE_CALCULATED,EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED

  PUBLIC EQUATIONS_HESSIAN_FINITE_DIFFERENCE_CALCULATED,EQUATIONS_HESSIAN_ANALYTIC_CALCULATED

  PUBLIC EquationsMatrices_DynamicLumpingTypeSet

  PUBLIC EquationsMatrices_DynamicStorageTypeSet

  PUBLIC EquationsMatrices_DynamicStructureTypeSet

  !!TODO check if the elements should be create/destroy rather than initialise/finalise
  PUBLIC EquationsMatrices_ElementInitialise,EquationsMatrices_ElementFinalise

  PUBLIC EquationsMatrices_ElementAdd

  PUBLIC EquationsMatrices_JacobianElementAdd

  PUBLIC EquationsMatrices_ElementCalculate

  PUBLIC EquationsMatrices_ElementMatrixFinalise,EquationsMatrices_ElementMatrixInitialise
  
  PUBLIC EquationsMatrices_ElementMatrixCalculate

  PUBLIC EquationsMatrices_ElementMatrixSetup

  PUBLIC EquationsMatrices_ElementVectorFinalise,EquationsMatrices_ElementVectorInitialise

  PUBLIC EquationsMatrices_ElementVectorCalculate

  PUBLIC EquationsMatrices_ElementVectorSetup

  PUBLIC EquationsMatrices_HessianOutput

  PUBLIC EquationsMatrices_JacobianNodeAdd

  PUBLIC EquationsMatrices_JacobianOutput

  PUBLIC EquationsMatrices_JacobianTypesSet

  PUBLIC EquationsMatrices_LinearStorageTypeSet

  PUBLIC EquationsMatrices_LinearStructureTypeSet

  PUBLIC EquationsMatrices_NodalInitialise,EquationsMatrices_NodalFinalise

  PUBLIC EquationsMatrices_NodeAdd

  PUBLIC EquationsMatrices_NodalCalculate

  PUBLIC EquationsMatrices_NodalMatrixFinalise,EquationsMatrices_NodalMatrixInitialise

  PUBLIC EquationsMatrices_NodalMatrixCalculate

  PUBLIC EquationsMatrices_NodalMatrixSetup

  PUBLIC EquationsMatrices_NodalVectorFinalise,EquationsMatrices_NodalVectorInitialise

  PUBLIC EquationsMatrices_NodalVectorCalculate

  PUBLIC EquationsMatrices_NodalVectorSetup

  PUBLIC EquationsMatrices_NonlinearStorageTypeSet

  PUBLIC EquationsMatrices_NonlinearStructureTypeSet

  PUBLIC EquationsMatrices_ScalarDestroy

  PUBLIC EquationsMatrices_VectorCreateFinish,EquationsMatrices_VectorCreateStart

  PUBLIC EquationsMatrices_VectorDestroy

  PUBLIC EquationsMatrices_VectorOutput

  PUBLIC EquationsMatrices_VectorValuesInitialise


CONTAINS

  !
  !================================================================================================================================
  !

  !>Finalise the equations Hessian and deallocate all memory
  SUBROUTINE EquationsMatrices_HessianFinalise(hessian,err,error,*)

    !Argument variables
    TYPE(EquationsHessianType), POINTER :: hessian !<A pointer to the equations Hessian to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("EquationsMatrices_HessianFinalise",err,error,*999)

    IF(.NOT.ASSOCIATED(hessian)) THEN
      IF(ASSOCIATED(hessian%hessian)) CALL DistributedMatrix_Destroy(hessian%hessian,err,error,*999)
      CALL EquationsMatrices_ElementMatrixFinalise(hessian%elementHessian,err,error,*999)
      DEALLOCATE(hessian)
    ENDIF
    
    EXITS("EquationsMatrices_HessianFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_HessianFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_HessianFinalise

  !
  !================================================================================================================================
  !

  !>Finalise the equations Jacobian and deallocate all memory
  SUBROUTINE EquationsMatrices_JacobianFinalise(equationsJacobian,err,error,*)

    !Argument variables
    TYPE(EquationsJacobianType), POINTER :: equationsJacobian !<A pointer to the equations Jacobian to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("EquationsMatrices_JacobianFinalise",err,error,*999)

    IF(ASSOCIATED(equationsJacobian)) THEN
      IF(ASSOCIATED(equationsJacobian%jacobian)) CALL DistributedMatrix_Destroy(equationsJacobian%jacobian,err,error,*999)
      CALL EquationsMatrices_ElementMatrixFinalise(equationsJacobian%elementJacobian,err,error,*999)
      CALL EquationsMatrices_NodalMatrixFinalise(equationsJacobian%nodalJacobian,err,error,*999)
      DEALLOCATE(equationsJacobian)
    ENDIF
    
    EXITS("EquationsMatrices_JacobianFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_JacobianFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_JacobianFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the equations Jacobian.
  SUBROUTINE EquationsMatrices_JacobianInitialise(nonlinearMatrices,matrixNumber,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices !<A pointer to the equations matrices nonlinear matrices to initialise the Jacobian for
    INTEGER(INTG), INTENT(IN) :: matrixNumber !<The index of the Jacobian matrix to initialise for the nonlinear matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("EquationsMatrices_JacobianInitialise",err,error,*998)
 
    IF(.NOT.ASSOCIATED(nonlinearMatrices)) CALL FlagError("Nonlinear matrices is not associated.",err,error,*998)
    IF(.NOT.ALLOCATED(nonlinearMatrices%jacobians)) CALL FlagError("Nonlinear matrices Jacobians is not allocated.",err,error,*999)
    NULLIFY(nonlinearMapping)
    CALL EquationsMatricesNonlinear_NonlinearMappingGet(nonlinearMatrices,nonlinearMapping,err,error,*999)
    IF(matrixNumber<1.OR.matrixNumber>nonlinearMatrices%numberOfJacobians) THEN
      localError="The matrix number of "//TRIM(NumberToVString(matrixNumber,"*",err,error))// &
        & " is invalid. The number must be >= 1 and <= "// &
        & TRIM(NumberToVString(nonlinearMatrices%numberOfJacobians,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    IF(ASSOCIATED(nonlinearMatrices%jacobians(matrixNumber)%ptr)) THEN
      localError="Nonlinear matrices Jacobian is already associated for matrix number "// &
        & TRIM(NumberToVString(matrixNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*998)
    ENDIF
     
    ALLOCATE(nonlinearMatrices%jacobians(matrixNumber)%ptr,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations Jacobian.",err,error,*999)
    nonlinearMatrices%jacobians(matrixNumber)%ptr%jacobianNumber=matrixNumber
    nonlinearMatrices%jacobians(matrixNumber)%ptr%nonlinearMatrices=>nonlinearMatrices
    nonlinearMatrices%jacobians(matrixNumber)%ptr%storageType=MATRIX_BLOCK_STORAGE_TYPE
    nonlinearMatrices%jacobians(matrixNumber)%ptr%structureType=EQUATIONS_MATRIX_NO_STRUCTURE
    nonlinearMatrices%jacobians(matrixNumber)%ptr%numberOfColumns= &
      & nonlinearMapping%jacobianToVarMap(matrixNumber)%numberOfColumns
    nonlinearMatrices%jacobians(matrixNumber)%ptr%updateJacobian=.TRUE.
    nonlinearMatrices%jacobians(matrixNumber)%ptr%firstAssembly=.TRUE.
    nonlinearMapping%jacobianToVarMap(matrixNumber)%jacobian=>nonlinearMatrices%jacobians(matrixNumber)%ptr
    NULLIFY(nonlinearMatrices%jacobians(matrixNumber)%ptr%jacobian)
    CALL EquationsMatrices_ElementMatrixInitialise(nonlinearMatrices%jacobians(matrixNumber)%ptr%elementJacobian,err,error,*999)
    CALL EquationsMatrices_NodalMatrixInitialise(nonlinearMatrices%jacobians(matrixNumber)%ptr%nodalJacobian,err,error,*999)
    nonlinearMatrices%jacobians(matrixNumber)%ptr%jacobianCalculationType=EQUATIONS_JACOBIAN_FINITE_DIFFERENCE_CALCULATED
    
    EXITS("EquationsMatrices_JacobianInitialise")
    RETURN
999 CALL EquationsMatrices_JacobianFinalise(nonlinearMatrices%jacobians(matrixNumber)%ptr,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_JacobianInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_JacobianInitialise

  !
  !================================================================================================================================
  !

  !>Finishes the creation of the vector equations matrices and RHS for the vector equations
  SUBROUTINE EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<The pointer to the vector equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string  
    !Local Variables
    INTEGER(INTG) :: dummyErr,matrixIdx,numberOfNonZeros
    INTEGER(INTG), POINTER :: rowIndices(:),columnIndices(:)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: rowDomainMap,columnDomainMap
    TYPE(EquationsJacobianType), POINTER :: jacobianMatrix
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(VARYING_STRING) :: dummyError,localError
    TYPE(LinkedList), POINTER :: list(:)
    
    NULLIFY(rowIndices)
    NULLIFY(columnIndices)

    ENTERS("EquationsMatrices_VectorCreateFinish",err,error,*998)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*998)
    IF(vectorMatrices%vectorMatricesFinished) CALL FlagError("Vector equations matrices have already been finished.",err,error,*998)
    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*999)    
    rowDomainMap=>vectorMapping%rowDOFSMapping
    IF(.NOT.ASSOCIATED(rowDomainMap)) CALL FlagError("Row domain map is not associated.",err,error,*999)
    
    dynamicMatrices=>vectorMatrices%dynamicMatrices
    IF(ASSOCIATED(dynamicMatrices)) THEN
      !Dynamic matrices
      NULLIFY(dynamicMapping)
      CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      !Now create the individual dynamic equations matrices
      DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
        NULLIFY(equationsMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,equationsMatrix,err,error,*999)
        columnDomainMap=>dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%columnDOFSMapping
        IF(.NOT.ASSOCIATED(columnDomainMap)) THEN
          localError="Column domain map for dynamic matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
            & " is not associated."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        !Create the distributed equations matrix
        CALL DistributedMatrix_CreateStart(rowDomainMap,columnDomainMap,vectorMatrices%dynamicMatrices%matrices(matrixIdx)% &
          & ptr%matrix,err,error,*999)
        CALL DistributedMatrix_DataTypeSet(equationsMatrix%matrix,MATRIX_VECTOR_DP_TYPE,err,error,*999)
        CALL DistributedMatrix_StorageTypeSet(equationsMatrix%matrix,equationsMatrix%storageType,err,error,*999)
        !Calculate and set the matrix structure/sparsity pattern
        IF(equationsMatrix%storageType/=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE.AND. &
          & equationsMatrix%storageType/=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE) THEN
          NULLIFY(list)
          CALL EquationsMatrix_StructureCalculate(equationsMatrix,numberOfNonZeros,rowIndices,columnIndices,list,err,error,*999)
          CALL DistributedMatrix_LinkListSet(equationsMatrix%matrix,list,err,error,*999)
          CALL DistributedMatrix_NumberOfNonZerosSet(equationsMatrix%matrix,numberOfNonZeros,err,error,*999)
          CALL DistributedMatrix_StorageLocationsSet(equationsMatrix%matrix,rowIndices,columnIndices,err,error,*999)
          IF(ASSOCIATED(rowIndices)) DEALLOCATE(rowIndices)
          IF(ASSOCIATED(columnIndices)) DEALLOCATE(columnIndices)
        ENDIF
        CALL DistributedMatrix_CreateFinish(equationsMatrix%matrix,err,error,*999)
      ENDDO !matrixIdx                
    ENDIF
    linearMatrices=>vectorMatrices%linearMatrices
    IF(ASSOCIATED(linearMatrices)) THEN
      !Linear matrices
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      !Now create the individual linear equations matrices
      DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
        NULLIFY(equationsMatrix)
        CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,equationsMatrix,err,error,*999)
        columnDomainMap=>linearMapping%equationsMatrixToVarMaps(matrixIdx)%columnDOFSMapping
        IF(.NOT.ASSOCIATED(columnDomainMap)) THEN
          localError="Column domain map for linear matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
            & " is not associated."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        !Create the distributed equations matrix
        CALL DistributedMatrix_CreateStart(rowDomainMap,columnDomainMap,vectorMatrices%linearMatrices%matrices(matrixIdx)% &
          & ptr%matrix,err,error,*999)
        CALL DistributedMatrix_DataTypeSet(equationsMatrix%matrix,MATRIX_VECTOR_DP_TYPE,err,error,*999)
        CALL DistributedMatrix_StorageTypeSet(equationsMatrix%matrix,equationsMatrix%storageType,err,error,*999)
        !Calculate and set the matrix structure/sparsity pattern
        IF(equationsMatrix%storageType/=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE.AND. &
          & equationsMatrix%storageType/=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE) THEN
          NULLIFY(list)
          CALL EquationsMatrix_StructureCalculate(equationsMatrix,numberOfNonZeros,rowIndices,columnIndices,list,err,error,*999)
          CALL DistributedMatrix_LinkListSet(equationsMatrix%matrix,LIST,err,error,*999)
          CALL DistributedMatrix_NumberOfNonZerosSet(equationsMatrix%matrix,numberOfNonZeros,err,error,*999)
          CALL DistributedMatrix_StorageLocationsSet(equationsMatrix%matrix,rowIndices,columnIndices,err,error,*999)
          IF(ASSOCIATED(rowIndices)) DEALLOCATE(rowIndices)
          IF(ASSOCIATED(columnIndices)) DEALLOCATE(columnIndices)
        ENDIF
        CALL DistributedMatrix_CreateFinish(equationsMatrix%matrix,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    nonlinearMatrices=>vectorMatrices%nonlinearMatrices
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      !Nonlinear matrices
      NULLIFY(nonlinearMapping)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
      !Set up the Jacobian matrices
      DO matrixIdx=1,nonlinearMatrices%numberOfJacobians
        NULLIFY(jacobianMatrix)
        CALL EquationsMatricesNonlinear_JacobianMatrixGet(nonlinearMatrices,matrixIdx,jacobianMatrix,err,error,*999)
        columnDomainMap=>nonlinearMapping%jacobianToVarMap(matrixIdx)%columnDOFSMapping
        IF(.NOT.ASSOCIATED(columnDomainMap)) THEN
          localError="Column domain map for Jacobian matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
            & " is not associated."
          CALL FlagError(localError,err,error,*999)
        ENDIF
!!TODO: Set the distributed matrix not to allocate the data if the Jacobian is not calculated.
        !Create the distributed Jacobian matrix
        CALL DistributedMatrix_CreateStart(rowDomainMap,columnDomainMap,jacobianMatrix%jacobian,err,error,*999)
        CALL DistributedMatrix_DataTypeSet(jacobianMatrix%jacobian,MATRIX_VECTOR_DP_TYPE,err,error,*999)
        CALL DistributedMatrix_StorageTypeSet(jacobianMatrix%jacobian,jacobianMatrix%storageType,err,error,*999)
        !Calculate and set the matrix structure/sparsity pattern
        IF(jacobianMatrix%storageType/=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE.AND. &
          & jacobianMatrix%storageType/=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE) THEN
          CALL JacobianMatrix_StructureCalculate(jacobianMatrix,numberOfNonZeros,rowIndices,columnIndices,err,error,*999)
          CALL DistributedMatrix_NumberOfNonZerosSet(jacobianMatrix%jacobian,numberOfNonZeros,err,error,*999)
          CALL DistributedMatrix_StorageLocationsSet(jacobianMatrix%jacobian,rowIndices,columnIndices,err,error,*999)
          IF(ASSOCIATED(rowIndices)) DEALLOCATE(rowIndices)
          IF(ASSOCIATED(columnIndices)) DEALLOCATE(columnIndices)
        ENDIF
        CALL DistributedMatrix_CreateFinish(jacobianMatrix%jacobian,err,error,*999)
      ENDDO !matrixIdx
      !Set up the residual vector                
      CALL DistributedVector_CreateStart(rowDomainMap,vectorMatrices%nonlinearMatrices%residual,err,error,*999)
      CALL DistributedVector_DataTypeSet(nonlinearMatrices%residual,MATRIX_VECTOR_DP_TYPE,err,error,*999)
      CALL DistributedVector_CreateFinish(nonlinearMatrices%residual,err,error,*999)
      !Initialise the residual vector to zero for time dependent problems so that the previous residual is set to zero
      CALL DistributedVector_AllValuesSet(nonlinearMatrices%residual,0.0_DP,err,error,*999)
    ENDIF
    rhsVector=>vectorMatrices%rhsVector
    IF(ASSOCIATED(rhsVector)) THEN
      !Set up the equations RHS vector          
      CALL DistributedVector_CreateStart(rowDomainMap,vectorMatrices%rhsVector%vector,err,error,*999)
      CALL DistributedVector_DataTypeSet(rhsVector%vector,MATRIX_VECTOR_DP_TYPE,err,error,*999)
      CALL DistributedVector_CreateFinish(rhsVector%vector,err,error,*999)
    ENDIF
    sourceVector=>vectorMatrices%sourceVector
    IF(ASSOCIATED(sourceVector)) THEN
      !Set up the equations source vector          
      CALL DistributedVector_CreateStart(rowDomainMap,vectorMatrices%sourceVector%vector,err,error,*999)
      CALL DistributedVector_DataTypeSet(sourceVector%vector,MATRIX_VECTOR_DP_TYPE,err,error,*999)
      CALL DistributedVector_CreateFinish(sourceVector%vector,err,error,*999)
    ENDIF
    !Finish up
    vectorMatrices%vectorMatricesFinished=.TRUE.
        
    EXITS("EquationsMatrices_VectorCreateFinish")
    RETURN
999 IF(ASSOCIATED(rowIndices)) DEALLOCATE(rowIndices)
    IF(ASSOCIATED(columnIndices)) DEALLOCATE(columnIndices)
    CALL EquationsMatrices_VectorFinalise(vectorMatrices,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_VectorCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_VectorCreateFinish

  !
  !================================================================================================================================
  !

  !>Starts the creation of the vector equations matrices and rhs for the vector equations
  SUBROUTINE EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsVectorType), POINTER :: vectorEquations !<The pointer to the vector equations to create the vector equations matrices for
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<On return, a pointer to the vector equations matrices being created.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string  
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(EquationsType), POINTER :: equations
    TYPE(VARYING_STRING) :: dummyError    

    ENTERS("EquationsMatrices_VectorCreateStart",err,error,*998)

    IF(ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is not associated.",err,error,*998)   
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*999)
    IF(.NOT.equations%equationsFinished) CALL FlagError("Vector equations equations has not been finished.",err,error,*999)

    NULLIFY(vectorMatrices)
    !Initialise the equations matrices
    CALL EquationsMatrices_VectorInitialise(vectorEquations,err,error,*999)
    vectorMatrices=>vectorEquations%vectorMatrices
    
    EXITS("EquationsMatrices_VectorCreateStart")
    RETURN
999 CALL EquationsMatrices_VectorFinalise(vectorEquations%vectorMatrices,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_VectorCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_VectorCreateStart

  !
  !================================================================================================================================
  !

  !>Destroy the vector equations matrices
  SUBROUTINE EquationsMatrices_VectorDestroy(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer the vector equations matrices to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMatrices_VectorDestroy",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated",err,error,*999)
    
    CALL EquationsMatrices_VectorFinalise(vectorMatrices,err,error,*999)
        
    EXITS("EquationsMatrices_VectorDestroy")
    RETURN
999 ERRORSEXITS("EquationsMatrices_VectorDestroy",err,error)    
    RETURN 1
   
  END SUBROUTINE EquationsMatrices_VectorDestroy

  !
  !================================================================================================================================
  !

  !>Calculate the positions in the equations matrices of the element matrix. Old CMISS name MELGE.
  SUBROUTINE EquationsMatrices_ElementMatrixCalculate(elementMatrix,updateMatrix,rowElementNumbers,columnElementNumbers, &
    & rowsFieldVariable,colsFieldVariable,err,error,*)

    !Argument variables
    TYPE(ElementMatrixType) :: elementMatrix !<The element matrix to calculate
    LOGICAL :: updateMatrix !<Is .TRUE. if the element matrix is to be updated, .FALSE. if not.
    INTEGER(INTG), INTENT(IN) :: rowElementNumbers(:) !<The row element number to calculate
    INTEGER(INTG), INTENT(IN) :: columnElementNumbers(:) !<The column element number to calculate
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: colsFieldVariable !<A pointer to the field variable associated with the columns
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,dataPointIdx,derivative,derivativeIdx,globalDOFIdx,localDOFIdx,localNodeIdx,node, &
      & version,localDataPointNumber,elementIdx,rowElementNumber,colElementNumber
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: elementsTopology
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMatrices_ElementMatrixCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(rowsFieldVariable)) CALL FlagError("Rows field variable is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(colsFieldVariable)) CALL FlagError("Columns field variable is not associated.",err,error,*999)
    
    elementMatrix%numberOfRows=0
    elementMatrix%numberOfColumns=0
    IF(updateMatrix) THEN
      IF(ASSOCIATED(rowsFieldVariable,colsFieldVariable)) THEN
        !Row and columns variable is the same.
        DO componentIdx=1,rowsFieldVariable%NUMBER_OF_COMPONENTS
          elementsTopology=>rowsFieldVariable%components(componentIdx)%domain%topology%elements
          DO elementIdx=1,SIZE(rowElementNumbers)
            rowElementNumber=rowElementNumbers(elementIdx)
            IF(rowElementNumber<1.OR.rowElementNumber>elementsTopology%TOTAL_NUMBER_OF_ELEMENTS) THEN
              localError="Element number "//TRIM(NumberToVString(rowElementNumber,"*",err,error))// &
                & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
                & " of rows field variable type "//TRIM(NumberToVString(rowsFieldVariable%VARIABLE_TYPE,"*",err,error))// &
                & ". The element number must be between 1 and "// &
                & TRIM(NumberToVString(elementsTopology%TOTAL_NUMBER_OF_ELEMENTS,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            SELECT CASE(rowsFieldVariable%components(componentIdx)%INTERPOLATION_TYPE)
            CASE(FIELD_CONSTANT_INTERPOLATION)
              localDOFIdx=rowsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
              globalDOFIdx=rowsFieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localDOFIdx)
              elementMatrix%numberOfRows=elementMatrix%numberOfRows+1
              elementMatrix%numberOfColumns=elementMatrix%numberOfColumns+1
              elementMatrix%rowDOFS(elementMatrix%numberOfRows)=localDOFIdx
              elementMatrix%columnDOFS(elementMatrix%numberOfColumns)=globalDOFIdx
            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
              localDOFIdx=rowsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP% &
                & elements(rowElementNumber)
              globalDOFIdx=rowsFieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localDOFIdx)
              elementMatrix%numberOfRows=elementMatrix%numberOfRows+1
              elementMatrix%numberOfColumns=elementMatrix%numberOfColumns+1
              elementMatrix%rowDOFS(elementMatrix%numberOfRows)=localDOFIdx
              elementMatrix%columnDOFS(elementMatrix%numberOfColumns)=globalDOFIdx
            CASE(FIELD_NODE_BASED_INTERPOLATION)
              basis=>elementsTopology%elements(rowElementNumber)%basis
              DO localNodeIdx=1,basis%NUMBER_OF_NODES
                node=elementsTopology%elements(rowElementNumber)%ELEMENT_NODES(localNodeIdx)
                DO derivativeIdx=1,basis%NUMBER_OF_DERIVATIVES(localNodeIdx)
                  derivative=elementsTopology%elements(rowElementNumber)%ELEMENT_DERIVATIVES(derivativeIdx,localNodeIdx)
                  version=elementsTopology%elements(rowElementNumber)%elementVersions(derivativeIdx,localNodeIdx)
                  localDOFIdx=rowsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node)% &
                    & derivatives(derivative)%versions(version)
                  globalDOFIdx=rowsFieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localDOFIdx)
                  elementMatrix%numberOfRows=elementMatrix%numberOfRows+1
                  elementMatrix%numberOfColumns=elementMatrix%numberOfColumns+1
                  elementMatrix%rowDOFS(elementMatrix%numberOfRows)=localDOFIdx
                  elementMatrix%columnDOFS(elementMatrix%numberOfColumns)=globalDOFIdx
                ENDDO !derivativeIdx
              ENDDO !localNodeIdx
            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
              decompositionData=>rowsFieldVariable%components(componentIdx)%domain%decomposition%topology%dataPoints
              DO dataPointIdx=1,decompositionData%elementDataPoint(rowElementNumber)%numberOfProjectedData
                localDataPointNumber=decompositionData%elementDataPoint(rowElementNumber)%dataIndices(dataPointIdx)%localNumber
                localDOFIdx=rowsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%DATA_POINT_PARAM2DOF_MAP% &
                  & DATA_POINTS(localDataPointNumber)
                globalDOFIdx=rowsFieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localDOFIdx)
                elementMatrix%numberOfRows=elementMatrix%numberOfRows+1
                elementMatrix%numberOfColumns=elementMatrix%numberOfColumns+1
                elementMatrix%rowDOFS(elementMatrix%numberOfRows)=localDOFIdx
                elementMatrix%columnDOFS(elementMatrix%numberOfColumns)=globalDOFIdx
              ENDDO !dataPointIdx
            CASE DEFAULT
              localError="The interpolation type of "// &
                & TRIM(NumberToVString(rowsFieldVariable%components(componentIdx)%INTERPOLATION_TYPE,"*",err,error))// &
                & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
                & " of rows field variable type "//TRIM(NumberToVString(rowsFieldVariable%VARIABLE_TYPE,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)          
            END SELECT
          ENDDO !elementIdx
        ENDDO !componentIdx
      ELSE
        !Row and column variables are different
        !Row mapping
        DO componentIdx=1,rowsFieldVariable%NUMBER_OF_COMPONENTS
          elementsTopology=>rowsFieldVariable%components(componentIdx)%domain%topology%elements
          DO elementIdx=1,SIZE(rowElementNumbers)
            rowElementNumber=rowElementNumbers(elementIdx)
            IF(rowElementNumber<1.OR.rowElementNumber>elementsTopology%TOTAL_NUMBER_OF_ELEMENTS) THEN
              localError="Row element number "//TRIM(NumberToVString(rowElementNumber,"*",err,error))// &
                & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
                & " of rows field variable type "//TRIM(NumberToVString(rowsFieldVariable%VARIABLE_TYPE,"*",err,error))// &
                & ". The element number must be between 1 and "// &
                & TRIM(NumberToVString(elementsTopology%TOTAL_NUMBER_OF_ELEMENTS,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            SELECT CASE(rowsFieldVariable%components(componentIdx)%INTERPOLATION_TYPE)
            CASE(FIELD_CONSTANT_INTERPOLATION)
              localDOFIdx=rowsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
              elementMatrix%numberOfRows=elementMatrix%numberOfRows+1
              elementMatrix%rowDOFS(elementMatrix%numberOfRows)=localDOFIdx
            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
              localDOFIdx=rowsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP% &
                & elements(rowElementNumber)
              elementMatrix%numberOfRows=elementMatrix%numberOfRows+1
              elementMatrix%rowDOFS(elementMatrix%numberOfRows)=localDOFIdx
            CASE(FIELD_NODE_BASED_INTERPOLATION)
              basis=>elementsTopology%elements(rowElementNumber)%basis
              DO localNodeIdx=1,basis%NUMBER_OF_NODES
                node=elementsTopology%elements(rowElementNumber)%ELEMENT_NODES(localNodeIdx)
                DO derivativeIdx=1,basis%NUMBER_OF_DERIVATIVES(localNodeIdx)
                  derivative=elementsTopology%elements(rowElementNumber)%ELEMENT_DERIVATIVES(derivativeIdx,localNodeIdx)
                  version=elementsTopology%elements(rowElementNumber)%elementVersions(derivativeIdx,localNodeIdx)
                  localDOFIdx=rowsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%nodes(node)% &
                    & derivatives(derivative)%versions(version)
                  elementMatrix%numberOfRows=elementMatrix%numberOfRows+1
                  elementMatrix%rowDOFS(elementMatrix%numberOfRows)=localDOFIdx
                ENDDO !derivativeIdx
              ENDDO !localNodeIdx
            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
              decompositionData=>rowsFieldVariable%components(componentIdx)%domain%decomposition%topology%dataPoints
              DO dataPointIdx=1,decompositionData%elementDataPoint(colElementNumber)%numberOfProjectedData
                localDataPointNumber=decompositionData%elementDataPoint(colElementNumber)%dataIndices(dataPointIdx)%localNumber
                localDOFIdx=rowsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%DATA_POINT_PARAM2DOF_MAP% &
                  & DATA_POINTS(localDataPointNumber)
                globalDOFIdx=rowsFieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localDOFIdx)
                elementMatrix%numberOfColumns=elementMatrix%numberOfColumns+1
                elementMatrix%columnDOFS(elementMatrix%numberOfColumns)=globalDOFIdx
              ENDDO !dataPointIdx
            CASE DEFAULT
              localError="The interpolation type of "// &
                & TRIM(NumberToVString(rowsFieldVariable%components(componentIdx)%INTERPOLATION_TYPE,"*",err,error))// &
                & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
                & " of rows field variable type "//TRIM(NumberToVString(rowsFieldVariable%VARIABLE_TYPE,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)          
            END SELECT
          ENDDO !elementIdx
        ENDDO !componentIdx
        !Column mapping
        DO componentIdx=1,colsFieldVariable%NUMBER_OF_COMPONENTS
          elementsTopology=>colsFieldVariable%components(componentIdx)%domain%topology%elements
          DO elementIdx=1,SIZE(columnElementNumbers)
            colElementNumber=columnElementNumbers(elementIdx)
            IF(colElementNumber<1.AND.colElementNumber>elementsTopology%TOTAL_NUMBER_OF_ELEMENTS) THEN
              localError="Column element number "//TRIM(NumberToVString(colElementNumber,"*",err,error))// &
                & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
                & " of column field variable type "//TRIM(NumberToVString(colsFieldVariable%VARIABLE_TYPE,"*",err,error))// &
                & ". The element number must be between 1 and "// &
                & TRIM(NumberToVString(elementsTopology%TOTAL_NUMBER_OF_ELEMENTS,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            SELECT CASE(colsFieldVariable%components(componentIdx)%INTERPOLATION_TYPE)
            CASE(FIELD_CONSTANT_INTERPOLATION)
              localDOFIdx=colsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
              globalDOFIdx=colsFieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localDOFIdx)
              elementMatrix%numberOfColumns=elementMatrix%numberOfColumns+1
              elementMatrix%columnDOFS(elementMatrix%numberOfColumns)=globalDOFIdx
            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
              localDOFIdx=colsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP% &
                & elements(colElementNumber)
              globalDOFIdx=colsFieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localDOFIdx)
              elementMatrix%numberOfColumns=elementMatrix%numberOfColumns+1
              elementMatrix%columnDOFS(elementMatrix%numberOfColumns)=globalDOFIdx
            CASE(FIELD_NODE_BASED_INTERPOLATION)
              basis=>elementsTopology%elements(colElementNumber)%basis
              DO localNodeIdx=1,basis%NUMBER_OF_NODES
                node=elementsTopology%elements(colElementNumber)%ELEMENT_NODES(localNodeIdx)
                DO derivativeIdx=1,basis%NUMBER_OF_DERIVATIVES(localNodeIdx)
                  derivative=elementsTopology%elements(colElementNumber)%ELEMENT_DERIVATIVES(derivativeIdx,localNodeIdx)
                  version=elementsTopology%elements(colElementNumber)%elementVersions(derivativeIdx,localNodeIdx)
                  localDOFIdx=colsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node)% &
                    & derivatives(derivative)%versions(version)
                  globalDOFIdx=colsFieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localDOFIdx)
                  elementMatrix%numberOfColumns=elementMatrix%numberOfColumns+1
                  elementMatrix%columnDOFS(elementMatrix%numberOfColumns)=globalDOFIdx
                ENDDO !derivativeIdx
              ENDDO !localNodeIdx
            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
              decompositionData=>colsFieldVariable%components(componentIdx)%domain%decomposition%topology%dataPoints
              DO dataPointIdx=1,decompositionData%elementDataPoint(colElementNumber)%numberOfProjectedData
                localDataPointNumber=decompositionData%elementDataPoint(colElementNumber)%dataIndices(dataPointIdx)%localNumber
                localDOFIdx=colsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%DATA_POINT_PARAM2DOF_MAP% &
                  & DATA_POINTS(localDataPointNumber)
                globalDOFIdx=colsFieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localDOFIdx)
                elementMatrix%numberOfColumns=elementMatrix%numberOfColumns+1
                elementMatrix%columnDOFS(elementMatrix%numberOfColumns)=globalDOFIdx
              ENDDO !dataPointIdx
            CASE DEFAULT
              localError="The interpolation type of "// &
                & TRIM(NumberToVString(colsFieldVariable%components(componentIdx)%INTERPOLATION_TYPE,"*",err,error))// &
                & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
                & " of column field variable type "//TRIM(NumberToVString(colsFieldVariable%VARIABLE_TYPE,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)          
            END SELECT
          ENDDO !elementIdx
        ENDDO !componentIdx
      ENDIF
      elementMatrix%matrix=0.0_DP
    ENDIF
    
    EXITS("EquationsMatrices_ElementMatrixCalculate")
    RETURN
999 ERRORSEXITS("EquationsMatrices_ElementMatrixCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_ElementMatrixCalculate

  !
  !================================================================================================================================
  !

  !>Finalise an element matrix and deallocate all memory
  SUBROUTINE EquationsMatrices_ElementMatrixFinalise(elementMatrix,err,error,*)

    !Argument variables
    TYPE(ElementMatrixType):: elementMatrix!<The element matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("EquationsMatrices_ElementMatrixFinalise",err,error,*999)

    elementMatrix%maxNumberOfRows=0
    elementMatrix%maxNumberOfColumns=0
    IF(ALLOCATED(elementMatrix%rowDOFS)) DEALLOCATE(elementMatrix%rowDOFS)
    IF(ALLOCATED(elementMatrix%columnDOFS)) DEALLOCATE(elementMatrix%columnDOFS)
    IF(ALLOCATED(elementMatrix%matrix)) DEALLOCATE(elementMatrix%matrix)
    
    EXITS("EquationsMatrices_ElementMatrixFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_ElementMatrixFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_ElementMatrixFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the element matrix.
  SUBROUTINE EquationsMatrices_ElementMatrixInitialise(elementMatrix,err,error,*)

    !Argument variables
    TYPE(ElementMatrixType) :: elementMatrix !The element matrix to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMatrices_ElementMatrixInitialise",err,error,*999)

    elementMatrix%equationsMatrixNumber=0
    elementMatrix%numberOfRows=0
    elementMatrix%numberOfColumns=0
    elementMatrix%maxNumberOfRows=0
    elementMatrix%maxNumberOfColumns=0
       
    EXITS("EquationsMatrices_ElementMatrixInitialise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_ElementMatrixInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_ElementMatrixInitialise

  !
  !================================================================================================================================
  !

  !>Sets up the element matrix for the row and column field variables.
  SUBROUTINE EquationsMatrices_ElementMatrixSetup(elementMatrix,rowsFieldVariable,columnsFieldVariable,rowsNumberOfElements, &
    & colsNumberOfElements,err,error,*)

    !Argument variables
    TYPE(ElementMatrixType) :: elementMatrix !<The element matrix to setup
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: columnsFieldVariable !<A pointer to the field variable associated with the columns
    INTEGER(INTG), INTENT(IN)  :: rowsNumberOfElements !<Number of elements in the row variables whose dofs are present in this element matrix
    INTEGER(INTG), INTENT(IN)  :: colsNumberOfElements !<Number of elements in the col variables whose dofs are present in this element matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr, componentIdx
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMatrices_ElementMatrixSetup",err,error,*998)

    IF(.NOT.ASSOCIATED(rowsFieldVariable)) CALL FlagError("Rows field variable is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(columnsFieldVariable)) CALL FlagError("Columns field variable is not associated.",err,error,*999)
    IF(ALLOCATED(elementMatrix%rowDOFS)) CALL FlagError("Element matrix row dofs already allocated.",err,error,*999)
    IF(ALLOCATED(elementMatrix%columnDOFS)) CALL FlagError("Element matrix column dofs already allocated.",err,error,*999)
    IF(ALLOCATED(elementMatrix%matrix)) CALL FlagError("Element matrix already allocated.",err,error,*999)
     
    elementMatrix%maxNumberOfRows=0
    DO componentIdx=1,rowsFieldVariable%NUMBER_OF_COMPONENTS
      elementMatrix%maxNumberOfRows=elementMatrix%maxNumberOfRows+ &
        & rowsFieldVariable%components(componentIdx)%maxNumberElementInterpolationParameters
    ENDDO !componentIdx
    elementMatrix%maxNumberOfRows=elementMatrix%maxNumberOfRows*rowsNumberOfElements
    elementMatrix%maxNumberOfColumns=0
    DO componentIdx=1,columnsFieldVariable%NUMBER_OF_COMPONENTS
      elementMatrix%maxNumberOfColumns=elementMatrix%maxNumberOfColumns+ &
        & columnsFieldVariable%components(componentIdx)%maxNumberElementInterpolationParameters
    ENDDO !componentIdx
    elementMatrix%maxNumberOfColumns=elementMatrix%maxNumberOfColumns*colsNumberOfElements
    ALLOCATE(elementMatrix%rowDOFS(elementMatrix%maxNumberOfRows),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate element matrix row dofs.",err,error,*999)
    ALLOCATE(elementMatrix%columnDOFS(elementMatrix%maxNumberOfColumns),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate element matrix column dofs.",err,error,*999)
    ALLOCATE(elementMatrix%matrix(elementMatrix%maxNumberOfRows,elementMatrix%maxNumberOfColumns),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate element matrix.",err,error,*999)
    
    EXITS("EquationsMatrices_ElementMatrixSetup")
    RETURN
999 CALL EquationsMatrices_ElementMatrixFinalise(elementMatrix,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_ElementMatrixSetup",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_ElementMatrixSetup

  !
  !================================================================================================================================
  !

  !>Calculate the positions in the equations rhs of the element rhs vector. Old CMISS name MELGE.
  SUBROUTINE EquationsMatrices_ElementVectorCalculate(elementVector,updateVector,elementNumber,rowsFieldVariable,err,error,*)

    !Argument variables
    TYPE(ElementVectorType) :: elementVector !<The element vector to calculate.
    LOGICAL :: updateVector !<Is .TRUE. if the element vector is to be updated, .FALSE. if not.
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,derivative,derivativeIdx,localDOFIdx,node,localNodeIdx,version,dataPointIdx,localDataPointNumber
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: elementsTopology
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsMatrices_ElementVectorCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(rowsFieldVariable)) CALL FlagError("Rows field variable is not associated.",err,error,*999)
    
    !Calculate the rows for the element vector
    elementVector%numberOfRows=0
    IF(updateVector) THEN
      DO componentIdx=1,rowsFieldVariable%NUMBER_OF_COMPONENTS
        elementsTopology=>rowsFieldVariable%components(componentIdx)%domain%topology%elements
        IF(elementNumber<1.OR.elementNumber>elementsTopology%TOTAL_NUMBER_OF_ELEMENTS) THEN
          localError="Element number "//TRIM(NumberToVString(elementNumber,"*",err,error))// &
            & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
            & " of rows field variable type "//TRIM(NumberToVString(rowsFieldVariable%VARIABLE_TYPE,"*",err,error))// &
            & ". The element number must be between 1 and "// &
            & TRIM(NumberToVString(elementsTopology%TOTAL_NUMBER_OF_ELEMENTS,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        SELECT CASE(rowsFieldVariable%components(componentIdx)%INTERPOLATION_TYPE)
        CASE(FIELD_CONSTANT_INTERPOLATION)
          localDOFIdx=rowsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
          elementVector%numberOfRows=elementVector%numberOfRows+1
          elementVector%rowDOFS(elementVector%numberOfRows)=localDOFIdx
        CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
          localDOFIdx=rowsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%elements(elementNumber)
          elementVector%numberOfRows=elementVector%numberOfRows+1
          elementVector%rowDOFS(elementVector%numberOfRows)=localDOFIdx
        CASE(FIELD_NODE_BASED_INTERPOLATION)
          basis=>elementsTopology%elements(elementNumber)%basis
          DO localNodeIdx=1,basis%NUMBER_OF_NODES
            node=elementsTopology%elements(elementNumber)%ELEMENT_NODES(localNodeIdx)
            DO derivativeIdx=1,basis%NUMBER_OF_DERIVATIVES(localNodeIdx)
              derivative=elementsTopology%elements(elementNumber)%ELEMENT_DERIVATIVES(derivativeIdx,localNodeIdx)
              version=elementsTopology%elements(elementNumber)%elementVersions(derivativeIdx,localNodeIdx)
              localDOFIdx=rowsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%nodes(node)% &
                & derivatives(derivative)%versions(version)
              elementVector%numberOfRows=elementVector%numberOfRows+1
              elementVector%rowDOFS(elementVector%numberOfRows)=localDOFIdx
            ENDDO !derivativeIdx
          ENDDO !localNodeIdx
        CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
          decompositionData=>rowsFieldVariable%components(componentIdx)%domain%decomposition%topology%dataPoints
          DO dataPointIdx=1,decompositionData%elementDataPoint(elementNumber)%numberOfProjectedData
            localDataPointNumber=decompositionData%elementDataPoint(elementNumber)% &
              & dataIndices(dataPointIdx)%localNumber
            localDOFIdx=rowsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%DATA_POINT_PARAM2DOF_MAP% &
              & DATA_POINTS(localDataPointNumber)
            elementVector%numberOfRows=elementVector%numberOfRows+1
            elementVector%rowDOFS(elementVector%numberOfRows)=localDOFIdx
          ENDDO !dataPointIdx
        CASE DEFAULT
          localError="The interpolation type of "// &
            & TRIM(NumberToVString(rowsFieldVariable%components(componentIdx)%INTERPOLATION_TYPE,"*",err,error))// &
            & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
            & " of rows field variable type "//TRIM(NumberToVString(rowsFieldVariable%VARIABLE_TYPE,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)          
        END SELECT
      ENDDO !componentIdx
      elementVector%vector=0.0_DP
    ENDIF
    
    EXITS("EquationsMatrices_ElementVectorCalculate")
    RETURN
999 ERRORSEXITS("EquationsMatrices_ElementVectorCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_ElementVectorCalculate

  !
  !================================================================================================================================
  !

  !>Finalise an element vector and deallocate all memory
  SUBROUTINE EquationsMatrices_ElementVectorFinalise(elementVector,err,error,*)

    !Argument variables
    TYPE(ElementVectorType) :: elementVector !<The element vector to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("EquationsMatrices_ElementVectorFinalise",err,error,*999)

    IF(ALLOCATED(elementVector%rowDOFS)) DEALLOCATE(elementVector%rowDOFS)
    IF(ALLOCATED(elementVector%vector)) DEALLOCATE(elementVector%vector)
    
    EXITS("EquationsMatrices_ElementVectorFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_ElementVectorFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_ElementVectorFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the element vector
  SUBROUTINE EquationsMatrices_ElementVectorInitialise(elementVector,err,error,*)

    !Argument variables
    TYPE(ElementVectorType) :: elementVector !The element vector to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMatrices_ElementVectorInitialise",err,error,*999)

    elementVector%numberOfRows=0
    elementVector%maxNumberOfRows=0
       
    EXITS("EquationsMatrices_ElementVectorInitialise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_ElementVectorInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_ElementVectorInitialise

  !
  !================================================================================================================================
  !

  !>Sets up the element vector for the row field variables.
  SUBROUTINE EquationsMatrices_ElementVectorSetup(elementVector,rowsFieldVariable,err,error,*)

    !Argument variables
    TYPE(ElementVectorType) :: elementVector !<The element vector to setup
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,componentIdx
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMatrices_ElementVectorSetup",err,error,*998)

    IF(.NOT.ASSOCIATED(rowsFieldVariable)) CALL FlagError("Rows field variable is not associated.",err,error,*999)
    IF(ALLOCATED(elementVector%rowDOFS)) CALL FlagError("Element vector row dofs is already allocated.",err,error,*999)
    IF(ALLOCATED(elementVector%vector)) CALL FlagError("Element vector vector already allocated.",err,error,*999)
   
    elementVector%maxNumberOfRows = 0
    DO componentIdx=1,rowsFieldVariable%NUMBER_OF_COMPONENTS
      elementVector%maxNumberOfRows=elementVector%maxNumberOfRows+ &
        & rowsFieldVariable%components(componentIdx)%maxNumberElementInterpolationParameters
    ENDDO !componentIdx
    ALLOCATE(elementVector%rowDOFS(elementVector%maxNumberOfRows),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate element vector row dofs.",err,error,*999)
    ALLOCATE(elementVector%vector(elementVector%maxNumberOfRows),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate element vector vector.",err,error,*999)
    
    EXITS("EquationsMatrices_ElementVectorSetup")
    RETURN
999 CALL EquationsMatrices_ElementVectorFinalise(elementVector,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_ElementVectorSetup",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_ElementVectorSetup

  !
  !================================================================================================================================
  !

  !>Adds the element matrices and rhs vector into the vector equations matrices and rhs vector.
  SUBROUTINE EquationsMatrices_ElementAdd(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,matrixIdx,rowIdx
    REAL(DP) :: sum
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("EquationsMatrices_ElementAdd()")
#endif

    ENTERS("EquationsMatrices_ElementAdd",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not allocated.",err,error,*999)
    
    dynamicMatrices=>vectorMatrices%dynamicMatrices
    IF(ASSOCIATED(dynamicMatrices)) THEN
      !Add the element matrices
      DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
        NULLIFY(equationsMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,equationsMatrix,err,error,*999)
        IF(equationsMatrix%updateMatrix) THEN
          !Handle lumped matrices
          IF(equationsMatrix%lumped) THEN
            DO rowIdx=1,equationsMatrix%elementMatrix%numberOfRows
              sum=0.0_DP
              DO columnIdx=1,equationsMatrix%elementMatrix%numberOfColumns
                sum=sum+equationsMatrix%elementMatrix%matrix(rowIdx,columnIdx)
                equationsMatrix%elementMatrix%matrix(rowIdx,columnIdx)=0.0_DP
              ENDDO !columnIdx
              equationsMatrix%elementMatrix%matrix(rowIdx,rowIdx)=sum
              !Add the element matrices into the distributed equations matrix
              CALL DistributedMatrix_ValuesAdd(equationsMatrix%matrix,equationsMatrix%elementMatrix%rowDOFS(rowIdx), &
                & equationsMatrix%elementMatrix%columnDOFS(rowIdx),equationsMatrix%elementMatrix%matrix(rowIdx,rowIdx), &
                & err,error,*999)
            ENDDO !rowIdx
          ELSE
            !Add the element matrice into the distributed equations matrix
            CALL DistributedMatrix_ValuesAdd(equationsMatrix%matrix,equationsMatrix%elementMatrix%rowDOFS(1: &
              & equationsMatrix%elementMatrix%numberOfRows),equationsMatrix%elementMatrix%columnDOFS(1: &
              & equationsMatrix%elementMatrix%numberOfColumns),equationsMatrix%elementMatrix%matrix(1: &
              & equationsMatrix%elementMatrix%numberOfRows,1:equationsMatrix%elementMatrix%numberOfColumns), &
              & err,error,*999)
          ENDIF
        ENDIF
      ENDDO !matrixIdx
    ENDIF
    linearMatrices=>vectorMatrices%linearMatrices
    IF(ASSOCIATED(linearMatrices)) THEN
      !Add the element matrices
      DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
        NULLIFY(equationsMatrix)
        CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,equationsMatrix,err,error,*999)
        IF(equationsMatrix%updateMatrix) THEN
          !Handle lumped matrices
          IF(equationsMatrix%lumped) THEN
            DO rowIdx=1,equationsMatrix%elementMatrix%numberOfRows
              sum=0.0_DP
              DO columnIdx=1,equationsMatrix%elementMatrix%numberOfColumns
                sum=sum+equationsMatrix%elementMatrix%matrix(rowIdx,columnIdx)
                equationsMatrix%elementMatrix%matrix(rowIdx,columnIdx)=0.0_DP
              ENDDO !columnIdx
              equationsMatrix%elementMatrix%matrix(rowIdx,rowIdx)=sum
              !Add the element matrice into the distributed equations matrix
              CALL DistributedMatrix_ValuesAdd(equationsMatrix%matrix,equationsMatrix%elementMatrix%rowDOFS(rowIdx), &
                & equationsMatrix%elementMatrix%columnDOFS(rowIdx),equationsMatrix%elementMatrix%matrix(rowIdx, &
                & rowIdx),err,error,*999)
            ENDDO !rowIdx
          ELSE
            !Add the element matrice into the distributed equations matrix
            CALL DistributedMatrix_ValuesAdd(equationsMatrix%matrix,equationsMatrix%elementMatrix%rowDOFS(1: &
              & equationsMatrix%elementMatrix%numberOfRows),equationsMatrix%elementMatrix%columnDOFS(1: &
              & equationsMatrix%elementMatrix%numberOfColumns),equationsMatrix%elementMatrix%matrix(1: &
              & equationsMatrix%elementMatrix%numberOfRows,1:equationsMatrix%elementMatrix%numberOfColumns), &
              & err,error,*999)
          ENDIF
        ENDIF
      ENDDO !matrixIdx
    ENDIF
    nonlinearMatrices=>vectorMatrices%nonlinearMatrices
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      IF(nonlinearMatrices%updateResidual) THEN
        !Add the residual element vector
        CALL DistributedVector_ValuesAdd(nonlinearMatrices%residual,nonlinearMatrices%elementResidual%rowDOFS(1: &
          & nonlinearMatrices%elementResidual%numberOfRows),nonlinearMatrices%elementResidual%vector(1:nonlinearMatrices% &
          & elementResidual%numberOfRows),err,error,*999)
      ENDIF
    ENDIF
    rhsVector=>vectorMatrices%rhsVector
    IF(ASSOCIATED(rhsVector)) THEN
      IF(rhsVector%updateVector) THEN
        !Add the rhs element vector
        CALL DistributedVector_ValuesAdd(rhsVector%vector,rhsVector%elementVector%rowDOFS(1: &
          & rhsVector%elementVector%numberOfRows),rhsVector%elementVector%vector(1:rhsVector% &
          & elementVector%numberOfRows),err,error,*999)
      ENDIF
    ENDIF
    sourceVector=>vectorMatrices%sourceVector
    IF(ASSOCIATED(sourceVector)) THEN
      IF(sourceVector%updateVector) THEN
        !Add the rhs element vector
        CALL DistributedVector_ValuesAdd(sourceVector%vector,sourceVector%elementVector%rowDOFS(1: &
          & sourceVector%elementVector%numberOfRows),sourceVector%elementVector%vector(1:sourceVector% &
          & elementVector%numberOfRows),err,error,*999)
      ENDIF
    ENDIF

#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("EquationsMatrices_ElementAdd()")
#endif
    
    EXITS("EquationsMatrices_ElementAdd")
    RETURN
999 ERRORSEXITS("EquationsMatrices_ElementAdd",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_ElementAdd

  !
  !================================================================================================================================
  !

  !>Calculate the positions in the vector equations matrices and rhs of the element matrices and rhs vector. Old CMISS name MELGE.
  SUBROUTINE EquationsMatrices_ElementCalculate(vectorMatrices,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate the mappings for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsJacobianType), POINTER :: jacobianMatrix
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable,colFieldVariable

#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("EquationsMatrices_ElementCalculate()")
#endif

    ENTERS("EquationsMatrices_ElementCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated",err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*999)

    dynamicMatrices=>vectorMatrices%dynamicMatrices
    IF(ASSOCIATED(dynamicMatrices)) THEN
      !Calculate the row and columns for the dynamic equations matrices
      NULLIFY(dynamicMapping)
      CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
        NULLIFY(equationsMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,equationsMatrix,err,error,*999)
        fieldVariable=>dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%variable
        CALL EquationsMatrices_ElementMatrixCalculate(equationsMatrix%elementMatrix,equationsMatrix%updateMatrix, &
          & [elementNumber],[elementNumber],fieldVariable,fieldVariable,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    linearMatrices=>vectorMatrices%linearMatrices
    IF(ASSOCIATED(linearMatrices)) THEN
      !Calculate the row and columns for the linear equations matrices
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
        NULLIFY(equationsMatrix)
        CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,equationsMatrix,err,error,*999)
        fieldVariable=>linearMapping%equationsMatrixToVarMaps(matrixIdx)%VARIABLE
        CALL EquationsMatrices_ElementMatrixCalculate(equationsMatrix%elementMatrix,equationsMatrix%updateMatrix, &
          & [elementNumber],[elementNumber],fieldVariable,fieldVariable,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    nonlinearMatrices=>vectorMatrices%nonlinearMatrices
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      !Calculate the rows and columns of the Jacobian
      NULLIFY(nonlinearMapping)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
      fieldVariable=>nonlinearMapping%jacobianToVarMap(1)%VARIABLE !Row field variable
      DO matrixIdx=1,nonlinearMatrices%numberOfJacobians
        NULLIFY(jacobianMatrix)
        CALL EquationsMatricesNonlinear_JacobianMatrixGet(nonlinearMatrices,matrixIdx,jacobianMatrix,err,error,*999)
        colFieldVariable=>nonlinearMapping%jacobianToVarMap(matrixIdx)%variable
        CALL EquationsMatrices_ElementMatrixCalculate(jacobianMatrix%elementJacobian,jacobianMatrix%updateJacobian, &
          & [elementNumber],[elementNumber],fieldVariable,colFieldVariable,err,error,*999)
      ENDDO !matrixIdx
      !Calculate the rows of the equations residual
      rhsMapping=>vectorMapping%rhsMapping
      IF(ASSOCIATED(rhsMapping)) THEN
        fieldVariable=>rhsMapping%rhsVariable
      ELSE
        fieldVariable=>nonlinearMapping%jacobianToVarMap(1)%variable
      ENDIF
      CALL EquationsMatrices_ElementVectorCalculate(nonlinearMatrices%elementResidual,nonlinearMatrices%updateResidual, &
        & elementNumber,fieldVariable,err,error,*999)
      nonlinearMatrices%elementResidualCalculated=0
    ENDIF
    rhsVector=>vectorMatrices%rhsVector
    IF(ASSOCIATED(rhsVector)) THEN
      NULLIFY(rhsMapping)
      CALL EquationsMappingVector_RHSMappingGet(vectorMapping,rhsMapping,err,error,*999)
      !Calculate the rows for the vector equations RHS
      fieldVariable=>rhsMapping%rhsVariable
      CALL EquationsMatrices_ElementVectorCalculate(rhsVector%elementVector,rhsVector%updateVector,elementNumber,fieldVariable, &
        & err,error,*999)
    ENDIF
    sourceVector=>vectorMatrices%sourceVector
    IF(ASSOCIATED(sourceVector)) THEN
      !Calculate the rows the equations source. The number of rows is not set by the source field so take the number of rows
      !from the RHS vector in the first instance.
      NULLIFY(rhsMapping)
      CALL EquationsMappingVector_RHSMappingGet(vectorMapping,rhsMapping,err,error,*999)
      fieldVariable=>rhsMapping%rhsVariable
      CALL EquationsMatrices_ElementVectorCalculate(sourceVector%elementVector,sourceVector%updateVector,elementNumber, &
        & fieldVariable,err,error,*999)
    ENDIF
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("EquationsMatrices_ElementCalculate()")
#endif
    
    EXITS("EquationsMatrices_ElementCalculate")
    RETURN
999 ERRORSEXITS("EquationsMatrices_ElementCalculate",err,error)
    RETURN 1
  END SUBROUTINE EquationsMatrices_ElementCalculate

  !
  !================================================================================================================================
  !

  !>Calculate the positions in the vector equations matrices and rhs of the nodal matrices and rhs vector. Old CMISS name MELGE.
  SUBROUTINE EquationsMatrices_NodalCalculate(vectorMatrices,nodeNumber,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(IN) :: nodeNumber !<The nodal number to calculate the mappings for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsJacobianType), POINTER :: jacobianMatrix
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable,colFieldVariable

#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("EquationsMatrices_NodalCalculate()")
#endif

    ENTERS("EquationsMatrices_NodalCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not allocated",err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*999)

    dynamicMatrices=>vectorMatrices%dynamicMatrices
    IF(ASSOCIATED(dynamicMatrices)) THEN
      !Calculate the row and columns for the dynamic equations matrices
      NULLIFY(dynamicMapping)
      CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
        NULLIFY(equationsMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,equationsMatrix,err,error,*999)
        fieldVariable=>dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%VARIABLE
        CALL EquationsMatrices_NodalMatrixCalculate(equationsMatrix%nodalMatrix,equationsMatrix%updateMatrix, &
          & nodeNumber,nodeNumber,fieldVariable,fieldVariable,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    linearMatrices=>vectorMatrices%linearMatrices
    IF(ASSOCIATED(linearMatrices)) THEN
      !Calculate the row and columns for the linear equations matrices
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
        NULLIFY(equationsMatrix)
        CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,equationsMatrix,err,error,*999)
        fieldVariable=>linearMapping%equationsMatrixToVarMaps(matrixIdx)%VARIABLE
        CALL EquationsMatrices_NodalMatrixCalculate(equationsMatrix%nodalMatrix,equationsMatrix%updateMatrix, &
          & nodeNumber,nodeNumber,fieldVariable,fieldVariable,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    nonlinearMatrices=>vectorMatrices%nonlinearMatrices
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      !Calculate the rows and columns of the Jacobian
      NULLIFY(nonlinearMapping)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
      nonlinearMapping=>vectorMapping%nonlinearMapping
      DO matrixIdx=1,nonlinearMatrices%numberOfJacobians
        NULLIFY(jacobianMatrix)
        CALL EquationsMatricesNonlinear_JacobianMatrixGet(nonlinearMatrices,matrixIdx,jacobianMatrix,err,error,*999)
        colFieldVariable=>nonlinearMapping%jacobianToVarMap(matrixIdx)%VARIABLE
        CALL EquationsMatrices_NodalMatrixCalculate(jacobianMatrix%nodalJacobian,jacobianMatrix%updateJacobian, &
          & nodeNumber,nodeNumber,fieldVariable,colFieldVariable,err,error,*999)
      ENDDO !matrixIdx
      !Calculate the rows of the equations residual
      rhsMapping=>vectorMapping%rhsMapping
      IF(ASSOCIATED(rhsMapping)) THEN
        fieldVariable=>rhsMapping%rhsVariable
      ELSE
        fieldVariable=>nonlinearMapping%jacobianToVarMap(1)%variable
      ENDIF
      CALL EquationsMatrices_NodalVectorCalculate(nonlinearMatrices%nodalResidual,nonlinearMatrices%updateResidual, &
        & nodeNumber,fieldVariable,err,error,*999)
      nonlinearMatrices%nodalResidualCalculated=0
    ENDIF
    rhsVector=>vectorMatrices%rhsVector
    IF(ASSOCIATED(rhsVector)) THEN
      !Calculate the rows for the equations RHS
      NULLIFY(rhsMapping)
      CALL EquationsMappingVector_RHSMappingGet(vectorMapping,rhsMapping,err,error,*999)
      fieldVariable=>rhsMapping%rhsVariable
      CALL EquationsMatrices_NodalVectorCalculate(rhsVector%nodalVector,rhsVector%updateVector,nodeNumber, &
        & fieldVariable,err,error,*999)
    ENDIF
    sourceVector=>vectorMatrices%sourceVector
    IF(ASSOCIATED(sourceVector)) THEN
      !Calculate the rows the equations source. The number of rows is not set by the source field so take the number of rows
      !from the RHS vector in the first instance.
      NULLIFY(rhsMapping)
      CALL EquationsMappingVector_RHSMappingGet(vectorMapping,rhsMapping,err,error,*999)
      fieldVariable=>rhsMapping%rhsVariable
      CALL EquationsMatrices_NodalVectorCalculate(sourceVector%nodalVector,sourceVector%updateVector,nodeNumber, &
        & fieldVariable,err,error,*999)
    ENDIF
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("EquationsMatrices_NodalCalculate()")
#endif
    
    EXITS("EquationsMatrices_NodalCalculate")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalCalculate",err,error)
    RETURN 1
  END SUBROUTINE EquationsMatrices_NodalCalculate

  !
  !================================================================================================================================
  !

  !>Calculate the positions in the equations matrices of the nodal matrix.
  SUBROUTINE EquationsMatrices_NodalMatrixCalculate(nodalMatrix,updateMatrix,rowNodeNumber,columnNodeNumber, &
    & rowsFieldVariable,colsFieldVariable,err,error,*)

    !Argument variables
    TYPE(NodalMatrixType) :: nodalMatrix !<The nodal matrix to calculate
    LOGICAL :: updateMatrix !<Is .TRUE. if the nodal matrix is to be updated, .FALSE. if not.
    INTEGER(INTG), INTENT(IN) :: rowNodeNumber !<The row nodal number to calculate
    INTEGER(INTG), INTENT(IN) :: columnNodeNumber !<The column nodal number to calculate
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: colsFieldVariable !<A pointer to the field variable associated with the columns
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx
    INTEGER(INTG) :: localRow,globalRow,localColumn,globalColumn
    INTEGER(INTG) :: numberOfDerivatives,numberOfVersions,versionIdx,derivativeIdx
    TYPE(DOMAIN_NODES_TYPE), POINTER :: nodesTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMatrices_NodalMatrixCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(rowsFieldVariable)) CALL FlagError("Rows field variable is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(colsFieldVariable)) CALL FlagError("Columns field variable is not associated.",err,error,*999)
    
    nodalMatrix%numberOfRows=0
    nodalMatrix%numberOfColumns=0
    IF(updateMatrix) THEN
      IF(ASSOCIATED(rowsFieldVariable,colsFieldVariable)) THEN
        !Row and columns variable is the same.
        DO componentIdx=1,rowsFieldVariable%NUMBER_OF_COMPONENTS
          nodesTopology=>rowsFieldVariable%components(componentIdx)%domain%topology%nodes
          IF(rowNodeNumber<1.OR.rowNodeNumber>nodesTopology%TOTAL_NUMBER_OF_NODES) THEN
            localError="Nodal number "//TRIM(NumberToVString(rowNodeNumber,"*",err,error))// &
              & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
              & " of rows field variable type "//TRIM(NumberToVString(rowsFieldVariable%VARIABLE_TYPE,"*",err,error))// &
              & ". The nodal number must be between 1 and "// &
              & TRIM(NumberToVString(nodesTopology%TOTAL_NUMBER_OF_NODES,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          SELECT CASE(rowsFieldVariable%components(componentIdx)%INTERPOLATION_TYPE)
          CASE(FIELD_CONSTANT_INTERPOLATION)
            localRow=rowsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
            globalRow=rowsFieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localRow)
            nodalMatrix%numberOfRows=nodalMatrix%numberOfRows+1
            nodalMatrix%numberOfColumns=nodalMatrix%numberOfColumns+1
            nodalMatrix%rowDofs(nodalMatrix%numberOfRows)=localRow
            nodalMatrix%columnDofs(nodalMatrix%numberOfColumns)=globalRow
          CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
            localRow=rowsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%elements(rowNodeNumber)
            globalRow=rowsFieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localRow)
            nodalMatrix%numberOfRows=nodalMatrix%numberOfRows+1
            nodalMatrix%numberOfColumns=nodalMatrix%numberOfColumns+1
            nodalMatrix%rowDofs(nodalMatrix%numberOfRows)=localRow
            nodalMatrix%columnDofs(nodalMatrix%numberOfColumns)=globalRow
          CASE(FIELD_NODE_BASED_INTERPOLATION)
            numberOfDerivatives=rowsFieldVariable%components(componentIdx)%domain%topology%nodes%nodes(rowNodeNumber)% &
              & NUMBER_OF_DERIVATIVES
            DO derivativeIdx=1,numberOfDerivatives
              numberOfVersions=rowsFieldVariable%components(componentIdx)%domain%topology%NODES%NODES(rowNodeNumber)% &
                & derivatives(derivativeIdx)%numberOfVersions
              DO versionIdx=1,numberOfVersions
                localRow=rowsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% & 
                  & nodes(rowNodeNumber)%derivatives(derivativeIdx)%versions(versionIdx)
                globalRow=rowsFieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localRow)
                nodalMatrix%numberOfRows=nodalMatrix%numberOfRows+1
                nodalMatrix%numberOfColumns=nodalMatrix%numberOfColumns+1
                nodalMatrix%rowDofs(nodalMatrix%numberOfRows)=localRow
                nodalMatrix%columnDofs(nodalMatrix%numberOfColumns)=globalRow
              ENDDO !versionIdx
            ENDDO !derivativeIdx
          CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The interpolation type of "// &
              & TRIM(NumberToVString(rowsFieldVariable%components(componentIdx)%INTERPOLATION_TYPE,"*",err,error))// &
              & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
              & " of rows field variable type "//TRIM(NumberToVString(rowsFieldVariable%VARIABLE_TYPE,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)          
          END SELECT
        ENDDO !componentIdx
      ELSE
        !Row and column variables are different
        !Row mapping
        DO componentIdx=1,rowsFieldVariable%NUMBER_OF_COMPONENTS
          nodesTopology=>rowsFieldVariable%components(componentIdx)%domain%topology%nodes
          IF(rowNodeNumber<1.OR.rowNodeNumber>nodesTopology%TOTAL_NUMBER_OF_NODES) THEN
            localError="Row nodal number "//TRIM(NumberToVString(rowNodeNumber,"*",err,error))// &
              & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
              & " of rows field variable type "//TRIM(NumberToVString(rowsFieldVariable%VARIABLE_TYPE,"*",err,error))// &
              & ". The nodal number must be between 1 and "// &
              & TRIM(NumberToVString(nodesTopology%TOTAL_NUMBER_OF_NODES,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          SELECT CASE(rowsFieldVariable%components(componentIdx)%INTERPOLATION_TYPE)
          CASE(FIELD_CONSTANT_INTERPOLATION)
            localRow=rowsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
            nodalMatrix%numberOfRows=nodalMatrix%numberOfRows+1
            nodalMatrix%rowDofs(nodalMatrix%numberOfRows)=localRow
          CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
            localRow=rowsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%elements(rowNodeNumber)
            nodalMatrix%numberOfRows=nodalMatrix%numberOfRows+1
            nodalMatrix%rowDofs(nodalMatrix%numberOfRows)=localRow
          CASE(FIELD_NODE_BASED_INTERPOLATION)
            numberOfDerivatives=rowsFieldVariable%components(componentIdx)%domain%topology%nodes%nodes(rowNodeNumber)% &
              & NUMBER_OF_DERIVATIVES
            DO derivativeIdx=1,numberOfDerivatives
              numberOfVersions=colsFieldVariable%components(componentIdx)%domain%topology%nodes%nodes(rowNodeNumber)% &
                & derivatives(derivativeIdx)%numberOfVersions
              DO versionIdx=1,numberOfVersions
                localRow=rowsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% & 
                  & nodes(rowNodeNumber)%derivatives(derivativeIdx)%versions(versionIdx)
                nodalMatrix%numberOfRows=nodalMatrix%numberOfRows+1
                nodalMatrix%rowDofs(nodalMatrix%numberOfRows)=localRow
              ENDDO !versionIdx
            ENDDO !derivativeIdx
          CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The interpolation type of "// &
              & TRIM(NumberToVString(rowsFieldVariable%components(componentIdx)%INTERPOLATION_TYPE,"*",err,error))// &
              & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
              & " of rows field variable type "//TRIM(NumberToVString(rowsFieldVariable%VARIABLE_TYPE,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)          
          END SELECT
        ENDDO !componentIdx
        !Column mapping
        DO componentIdx=1,colsFieldVariable%NUMBER_OF_COMPONENTS
          nodesTopology=>colsFieldVariable%components(componentIdx)%domain%topology%nodes
          IF(columnNodeNumber<1.OR.columnNodeNumber>nodesTopology%TOTAL_NUMBER_OF_NODES) THEN
            localError="Column nodal number "//TRIM(NumberToVString(columnNodeNumber,"*",err,error))// &
              & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
              & " of column field variable type "//TRIM(NumberToVString(colsFieldVariable%VARIABLE_TYPE,"*",err,error))// &
              & ". The nodal number must be between 1 and "// &
              & TRIM(NumberToVString(nodesTopology%TOTAL_NUMBER_OF_NODES,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          SELECT CASE(colsFieldVariable%components(componentIdx)%INTERPOLATION_TYPE)
          CASE(FIELD_CONSTANT_INTERPOLATION)
            localColumn=colsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
            globalColumn=colsFieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localColumn)
            nodalMatrix%numberOfColumns=nodalMatrix%numberOfColumns+1
            nodalMatrix%columnDofs(nodalMatrix%numberOfColumns)=globalColumn
          CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
            localColumn=colsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%elements(columnNodeNumber)
            globalColumn=colsFieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localColumn)
            nodalMatrix%numberOfColumns=nodalMatrix%numberOfColumns+1
            nodalMatrix%columnDofs(nodalMatrix%numberOfColumns)=globalColumn
          CASE(FIELD_NODE_BASED_INTERPOLATION)
            numberOfDerivatives=colsFieldVariable%components(componentIdx)%domain%topology%nodes%nodes(rowNodeNumber)% &
              & NUMBER_OF_DERIVATIVES
            DO derivativeIdx=1,numberOfDerivatives
              numberOfVersions=colsFieldVariable%components(componentIdx)%domain%topology%nodes%nodes(rowNodeNumber)% &
                & derivatives(derivativeIdx)%numberOfVersions
              DO versionIdx=1,numberOfVersions
                localRow=colsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% & 
                  & nodes(rowNodeNumber)%derivatives(derivativeIdx)%versions(versionIdx)
                nodalMatrix%numberOfColumns=nodalMatrix%numberOfColumns+1
                nodalMatrix%columnDofs(nodalMatrix%numberOfColumns)=localRow
              ENDDO !versionIdx
            ENDDO !derivativeIdx
          CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The interpolation type of "// &
              & TRIM(NumberToVString(colsFieldVariable%components(componentIdx)%INTERPOLATION_TYPE,"*",err,error))// &
              & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
              & " of column field variable type "//TRIM(NumberToVString(colsFieldVariable%VARIABLE_TYPE,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)          
          END SELECT
        ENDDO !componentIdx
      ENDIF
      nodalMatrix%matrix=0.0_DP
    ENDIF
    
    EXITS("EquationsMatrices_NodalMatrixCalculate")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalMatrixCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NodalMatrixCalculate
  
  !
  !================================================================================================================================
  !

  !>Calculate the positions in the vector equations rhs of the nodal rhs vector.
  SUBROUTINE EquationsMatrices_NodalVectorCalculate(nodalVector,updateVector,rowNodeNumber,rowsFieldVariable,err,error,*)

    !Argument variables
    TYPE(NodalVectorType) :: nodalVector !<The nodal vector to calculate.
    LOGICAL :: updateVector !<Is .TRUE. if the nodal vector is to be updated, .FALSE. if not.
    INTEGER(INTG), INTENT(IN) :: rowNodeNumber !<The nodal number to calculate
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,localRow
    INTEGER(INTG) :: numberOfDerivatives,numberOfVersions,versionIdx,derivativeIdx
    TYPE(DOMAIN_NODES_TYPE), POINTER :: nodesTopology
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsMatrices_NodalVectorCalculate",err,error,*999)
    
    IF(.NOT.ASSOCIATED(rowsFieldVariable)) CALL FlagError("Rows field variable is not associated.",err,error,*999)
    
    !Calculate the rows for the nodal vector
    nodalVector%numberOfRows=0
    IF(updateVector) THEN
      DO componentIdx=1,rowsFieldVariable%NUMBER_OF_COMPONENTS
        nodesTopology=>rowsFieldVariable%components(componentIdx)%domain%topology%nodes
        IF(rowNodeNumber<1.OR.rowNodeNumber>nodesTopology%TOTAL_NUMBER_OF_NODES) THEN
          localError="Node number "//TRIM(NumberToVString(rowNodeNumber,"*",err,error))// &
            & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
            & " of rows field variable type "//TRIM(NumberToVString(rowsFieldVariable%VARIABLE_TYPE,"*",err,error))// &
            & ". The nodal number must be between 1 and "// &
            & TRIM(NumberToVString(nodesTopology%TOTAL_NUMBER_OF_NODES,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        SELECT CASE(rowsFieldVariable%components(componentIdx)%INTERPOLATION_TYPE)
        CASE(FIELD_CONSTANT_INTERPOLATION)
          localRow=rowsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
          nodalVector%numberOfRows=nodalVector%numberOfRows+1
          nodalVector%rowDofs(nodalVector%numberOfRows)=localRow
        CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
          localRow=rowsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%elements(rowNodeNumber)
          nodalVector%numberOfRows=nodalVector%numberOfRows+1
          nodalVector%rowDofs(nodalVector%numberOfRows)=localRow
        CASE(FIELD_NODE_BASED_INTERPOLATION)
          numberOfDerivatives=rowsFieldVariable%components(componentIdx)%domain%topology%nodes%nodes(rowNodeNumber)% &
            & NUMBER_OF_DERIVATIVES
          DO derivativeIdx=1,numberOfDerivatives
            numberOfVersions=rowsFieldVariable%components(componentIdx)%domain%topology%nodes%nodes(rowNodeNumber)% &
              & derivatives(derivativeIdx)%numberOfVersions
            DO versionIdx=1,numberOfVersions
              localRow=rowsFieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% & 
                & nodes(rowNodeNumber)%derivatives(derivativeIdx)%versions(versionIdx)
              nodalVector%numberOfRows=nodalVector%numberOfRows+1
              nodalVector%rowDofs(nodalVector%numberOfRows)=localRow
            ENDDO !versionIdx
          ENDDO !derivativeIdx
        CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The interpolation type of "// &
            & TRIM(NumberToVString(rowsFieldVariable%components(componentIdx)%INTERPOLATION_TYPE,"*",err,error))// &
            & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
            & " of rows field variable type "//TRIM(NumberToVString(rowsFieldVariable%VARIABLE_TYPE,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)          
        END SELECT
      ENDDO !componentIdx
      nodalVector%vector=0.0_DP
    ENDIF
    
    EXITS("EquationsMatrices_NodalVectorCalculate")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalVectorCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NodalVectorCalculate

  !
  !================================================================================================================================
  !

  !>Adds the nodal matrices and rhs vector into the equations matrices and rhs vector.
  SUBROUTINE EquationsMatrices_NodeAdd(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,matrixIdx,rowIdx
    REAL(DP) :: sum
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("EquationsMatrices_NodeAdd()")
#endif

    ENTERS("EquationsMatrices_NodeAdd",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)
    
    dynamicMatrices=>vectorMatrices%dynamicMatrices
    IF(ASSOCIATED(dynamicMatrices)) THEN
      !Add the nodal matrices
      DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
        NULLIFY(equationsMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,equationsMatrix,err,error,*999)
        IF(equationsMatrix%updateMatrix) THEN
          !Handle lumped matrices
          IF(equationsMatrix%lumped) THEN
            DO rowIdx=1,equationsMatrix%nodalMatrix%numberOfRows
              sum=0.0_DP
              DO columnIdx=1,equationsMatrix%nodalMatrix%numberOfColumns
                sum=sum+equationsMatrix%nodalMatrix%matrix(rowIdx,columnIdx)
                equationsMatrix%nodalMatrix%matrix(rowIdx,columnIdx)=0.0_DP
              ENDDO !columnIdx
              equationsMatrix%nodalMatrix%matrix(rowIdx,rowIdx)=sum
              !Add the nodal matrice into the distributed equations matrix
              CALL DistributedMatrix_ValuesAdd(equationsMatrix%matrix,equationsMatrix%nodalMatrix%rowDofs(rowIdx), &
                & equationsMatrix%nodalMatrix%columnDofs(rowIdx),equationsMatrix%nodalMatrix%matrix(rowIdx,rowIdx), &
                & err,error,*999)
            ENDDO !rowIdx
          ELSE
            !Add the nodal matrice into the distributed equations matrix
            CALL DistributedMatrix_ValuesAdd(equationsMatrix%matrix,equationsMatrix%nodalMatrix%rowDofs(1: &
              & equationsMatrix%nodalMatrix%numberOfRows),equationsMatrix%nodalMatrix%columnDofs(1: &
              & equationsMatrix%nodalMatrix%numberOfColumns),equationsMatrix%nodalMatrix%matrix(1: &
              & equationsMatrix%nodalMatrix%numberOfRows,1:equationsMatrix%nodalMatrix%numberOfColumns), &
              & err,error,*999)
          ENDIF
        ENDIF
      ENDDO !matrixIdx
    ENDIF
    linearMatrices=>vectorMatrices%linearMatrices
    IF(ASSOCIATED(linearMatrices)) THEN
      !Add the nodal matrices
      DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
        NULLIFY(equationsMatrix)
        CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,equationsMatrix,err,error,*999)
        IF(equationsMatrix%updateMatrix) THEN
          !Handle lumped matrices
          IF(equationsMatrix%lumped) THEN
            DO rowIdx=1,equationsMatrix%nodalMatrix%numberOfRows
              sum=0.0_DP
              DO columnIdx=1,equationsMatrix%nodalMatrix%numberOfColumns
                sum=sum+equationsMatrix%nodalMatrix%matrix(rowIdx,columnIdx)
                equationsMatrix%nodalMatrix%matrix(rowIdx,columnIdx)=0.0_DP
              ENDDO !columnIdx
              equationsMatrix%nodalMatrix%matrix(rowIdx,rowIdx)=sum
              !Add the nodal matrice into the distributed equations matrix
              CALL DistributedMatrix_ValuesAdd(equationsMatrix%matrix,equationsMatrix%nodalMatrix%rowDofs(rowIdx), &
                & equationsMatrix%nodalMatrix%columnDofs(rowIdx),equationsMatrix%nodalMatrix%matrix(rowIdx,rowIdx), &
                & err,error,*999)
            ENDDO !rowIdx
          ELSE
            !Add the nodal matrice into the distributed equations matrix
            CALL DistributedMatrix_ValuesAdd(equationsMatrix%matrix,equationsMatrix%nodalMatrix%rowDofs(1: &
              & equationsMatrix%nodalMatrix%numberOfRows),equationsMatrix%nodalMatrix%columnDofs(1: &
              & equationsMatrix%nodalMatrix%numberOfColumns),equationsMatrix%nodalMatrix%matrix(1: &
              & equationsMatrix%nodalMatrix%numberOfRows,1:equationsMatrix%nodalMatrix%numberOfColumns), &
              & err,error,*999)
          ENDIF
        ENDIF
      ENDDO !matrixIdx
    ENDIF
    nonlinearMatrices=>vectorMatrices%nonlinearMatrices
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      IF(nonlinearMatrices%updateResidual) THEN
        !Add the residual nodal vector
        CALL DistributedVector_ValuesAdd(nonlinearMatrices%residual,nonlinearMatrices%nodalResidual%rowDofs(1: &
          & nonlinearMatrices%nodalResidual%numberOfRows),nonlinearMatrices%nodalResidual%vector(1:nonlinearMatrices% &
          & NodalResidual%numberOfRows),err,error,*999)
      ENDIF
    ENDIF
    rhsVector=>vectorMatrices%rhsVector
    IF(ASSOCIATED(rhsVector)) THEN
      IF(rhsVector%updateVector) THEN
        !Add the rhs nodal vector
        CALL DistributedVector_ValuesAdd(rhsVector%vector,rhsVector%nodalVector%rowDofs(1: &
          & rhsVector%nodalVector%numberOfRows),rhsVector%nodalVector%vector(1:rhsVector% &
          & NodalVector%numberOfRows),err,error,*999)
      ENDIF
    ENDIF
    sourceVector=>vectorMatrices%sourceVector
    IF(ASSOCIATED(sourceVector)) THEN
      IF(sourceVector%updateVector) THEN
        !Add the rhs nodal vector
        CALL DistributedVector_ValuesAdd(sourceVector%vector,sourceVector%nodalVector%rowDofs(1: &
          & sourceVector%nodalVector%numberOfRows),sourceVector%nodalVector%vector(1:sourceVector% &
          & NodalVector%numberOfRows),err,error,*999)
      ENDIF
    ENDIF

#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("EquationsMatrices_NodeAdd()")
#endif
    
    EXITS("EquationsMatrices_NodeAdd")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodeAdd",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NodeAdd

  !
  !================================================================================================================================
  !

  !>Initialise the nodal calculation information for the equations matrices
  SUBROUTINE EquationsMatrices_NodalInitialise(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !The equations matrices to initialise the nodal information for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsJacobianType), POINTER :: jacobianMatrix
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable,colFieldVariable
    
    ENTERS("EquationsMatrices_NodalInitialise",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)

    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*999)
    
    dynamicMatrices=>vectorMatrices%dynamicMatrices
    IF(ASSOCIATED(dynamicMatrices)) THEN
      !Initialise the dynamic nodal matrices
      NULLIFY(dynamicMapping)
      CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
        NULLIFY(equationsMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,equationsMatrix,err,error,*999)
        fieldVariable=>dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%variable
        CALL EquationsMatrices_NodalMatrixSetup(equationsMatrix%nodalMatrix,fieldVariable,fieldVariable,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    linearMatrices=>vectorMatrices%linearMatrices
    IF(ASSOCIATED(linearMatrices)) THEN
      !Initialise the linear nodal matrices
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
        NULLIFY(equationsMatrix)
        CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,equationsMatrix,err,error,*999)
        fieldVariable=>linearMapping%equationsMatrixToVarMaps(matrixIdx)%variable
        CALL EquationsMatrices_NodalMatrixSetup(equationsMatrix%nodalMatrix,fieldVariable,fieldVariable,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    nonlinearMatrices=>vectorMatrices%nonlinearMatrices
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      !Initialise the Jacobian nodal matrices
      NULLIFY(nonlinearMapping)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
      fieldVariable=>nonlinearMapping%jacobianToVarMap(1)%variable
      DO matrixIdx=1,nonlinearMatrices%numberOfJacobians
        NULLIFY(jacobianMatrix)
        CALL EquationsMatricesNonlinear_JacobianMatrixGet(nonlinearMatrices,matrixIdx,jacobianMatrix,err,error,*999)
        colFieldVariable=>nonlinearMapping%jacobianToVarMap(matrixIdx)%variable
        CALL EquationsMatrices_NodalMatrixSetup(jacobianMatrix%nodalJacobian,fieldVariable,colFieldVariable,err,error,*999)
      ENDDO !matrixIdx
      !Use RHS variable for residual vector, otherwise first nonlinear variable if no RHS
      rhsMapping=>vectorMapping%rhsMapping
      IF(ASSOCIATED(rhsMapping)) THEN
        fieldVariable=>rhsMapping%rhsVariable
      ELSE
        fieldVariable=>nonlinearMapping%jacobianToVarMap(1)%variable
      ENDIF
      CALL EquationsMatrices_NodalVectorSetup(nonlinearMatrices%nodalResidual,fieldVariable,err,error,*999)
      nonlinearMatrices%nodalResidualCalculated=0
    ENDIF
    rhsVector=>vectorMatrices%rhsVector
    IF(ASSOCIATED(rhsVector)) THEN
      !Initialise the RHS nodal vector
      NULLIFY(rhsMapping)
      CALL EquationsMappingVector_RHSMappingGet(vectorMapping,rhsMapping,err,error,*999)
      fieldVariable=>rhsMapping%rhsVariable
      CALL EquationsMatrices_NodalVectorSetup(rhsVector%nodalVector,fieldVariable,err,error,*999)
    ENDIF
    sourceVector=>vectorMatrices%sourceVector
    IF(ASSOCIATED(sourceVector)) THEN
      !Initialise the source nodal vector. Note that the number of rows in the source vector is taken, for now, from the RHS
      !vector
      IF(ASSOCIATED(rhsVector)) THEN
        !Initialise the RHS nodal vector
        rhsMapping=>vectorMapping%rhsMapping
        NULLIFY(rhsMapping)
        CALL EquationsMappingVector_RHSMappingGet(vectorMapping,rhsMapping,err,error,*999)
        fieldVariable=>rhsMapping%rhsVariable
        CALL EquationsMatrices_NodalVectorSetup(sourceVector%nodalVector,fieldVariable,err,error,*999)
      ELSE
        CALL FlagError("Not implemented.",err,error,*999)
      ENDIF
    ENDIF
    
    EXITS("EquationsMatrices_NodalInitialise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NodalInitialise

  !
  !================================================================================================================================
  !

  !>Sets up the nodal matrix for the row and column field variables.
  SUBROUTINE EquationsMatrices_NodalMatrixSetup(nodalMatrix,rowsFieldVariable,colsFieldVariable,err,error,*)

    !Argument variables
    TYPE(NodalMatrixType) :: nodalMatrix !<The nodal matrix to setup
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: colsFieldVariable !<A pointer to the field variable associated with the columns
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMatrices_NodalMatrixSetup",err,error,*998)

    IF(.NOT.ASSOCIATED(rowsFieldVariable)) CALL FlagError("Rows field variable is not associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(colsFieldVariable)) CALL FlagError("Columns field variable is not associated.",err,error,*998)
    IF(ALLOCATED(nodalMatrix%rowDofs)) CALL FlagError("Nodal matrix row dofs already allocated.",err,error,*998)
    IF(ALLOCATED(nodalMatrix%columnDofs)) CALL FlagError("Nodal matrix column dofs already allocated.",err,error,*998)
    IF(ALLOCATED(nodalMatrix%matrix)) CALL FlagError("Nodal matrix already allocated.",err,error,*998)
    
    nodalMatrix%maxNumberOfRows=rowsFieldVariable%maxNumberNodeInterpolationParameters*rowsFieldVariable%NUMBER_OF_COMPONENTS
    nodalMatrix%maxNumberOfColumns=colsFieldVariable%maxNumberNodeInterpolationParameters*colsFieldVariable%NUMBER_OF_COMPONENTS
    ALLOCATE(nodalMatrix%rowDofs(nodalMatrix%maxNumberOfRows),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate nodal matrix row dofs.",err,error,*999)
    ALLOCATE(nodalMatrix%columnDofs(nodalMatrix%maxNumberOfColumns),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate nodal matrix column dofs.",err,error,*999)
    ALLOCATE(nodalMatrix%matrix(nodalMatrix%maxNumberOfRows,nodalMatrix%maxNumberOfColumns),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate nodal matrix.",err,error,*999)
    
    EXITS("EquationsMatrices_NodalMatrixSetup")
    RETURN
999 CALL EquationsMatrices_NodalMatrixFinalise(nodalMatrix,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_NodalMatrixSetup",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NodalMatrixSetup

  !
  !================================================================================================================================
  !

  !>Sets up the nodal vector for the row field variables.
  SUBROUTINE EquationsMatrices_NodalVectorSetup(nodalVector,rowsFieldVariable,err,error,*)

    !Argument variables
    TYPE(NodalVectorType) :: nodalVector !<The nodal vector to setup
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMatrices_NodalVectorSetup",err,error,*998)

    IF(.NOT.ASSOCIATED(rowsFieldVariable)) CALL FlagError("Rows field variable is not associated.",err,error,*998)
    IF(ALLOCATED(nodalVector%rowDofs)) CALL FlagError("Nodal vector row dofs is already allocated.",err,error,*998)
    IF(ALLOCATED(nodalVector%vector)) CALL FlagError("Nodal vector vector already allocated.",err,error,*998)
        
    nodalVector%maxNumberOfRows=rowsFieldVariable%maxNumberNodeInterpolationParameters*rowsFieldVariable%NUMBER_OF_COMPONENTS
    ALLOCATE(nodalVector%rowDofs(nodalVector%maxNumberOfRows),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate nodal vector row dofs.",err,error,*999)
    ALLOCATE(nodalVector%vector(nodalVector%maxNumberOfRows),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate nodal vector vector.",err,error,*999)
    
    EXITS("EquationsMatrices_NodalVectorSetup")
    RETURN
999 CALL EquationsMatrices_NodalVectorFinalise(nodalVector,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_NodalVectorSetup",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NodalVectorSetup

  !
  !================================================================================================================================
  !

  !>Finalise the nodal calculation information and deallocate all memory
  SUBROUTINE EquationsMatrices_NodalFinalise(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<The equations matrices for which to finalise the nodals
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsJacobianType), POINTER :: jacobianMatrix
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    
    ENTERS("EquationsMatrices_NodalFinalise",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)
    
    dynamicMatrices=>vectorMatrices%dynamicMatrices
    IF(ASSOCIATED(dynamicMatrices)) THEN
      !Finalise the dynamic nodal matrices
      DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
        NULLIFY(equationsMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,equationsMatrix,err,error,*999)
        equationsMatrix%nodalMatrix%maxNumberOfRows=0
        equationsMatrix%nodalMatrix%maxNumberOfColumns=0
        IF(ALLOCATED(equationsMatrix%nodalMatrix%rowDofs)) DEALLOCATE(equationsMatrix%nodalMatrix%rowDofs)
        IF(ALLOCATED(equationsMatrix%nodalMatrix%columnDofs)) DEALLOCATE(equationsMatrix%nodalMatrix%columnDofs)
        IF(ALLOCATED(equationsMatrix%nodalMatrix%matrix)) DEALLOCATE(equationsMatrix%nodalMatrix%matrix)
      ENDDO !matrixIdx
    ENDIF
    linearMatrices=>vectorMatrices%linearMatrices
    IF(ASSOCIATED(linearMatrices)) THEN
      !Finalise the linear nodal matrices
      DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
        NULLIFY(equationsMatrix)
        CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,equationsMatrix,err,error,*999)
        equationsMatrix%nodalMatrix%maxNumberOfRows=0
        equationsMatrix%nodalMatrix%maxNumberOfColumns=0
        IF(ALLOCATED(equationsMatrix%nodalMatrix%rowDofs)) DEALLOCATE(equationsMatrix%nodalMatrix%rowDofs)
        IF(ALLOCATED(equationsMatrix%nodalMatrix%columnDofs)) DEALLOCATE(equationsMatrix%nodalMatrix%columnDofs)
        IF(ALLOCATED(equationsMatrix%nodalMatrix%matrix)) DEALLOCATE(equationsMatrix%nodalMatrix%matrix)
      ENDDO !matrixIdx
    ENDIF
    nonlinearMatrices=>vectorMatrices%nonlinearMatrices
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      DO matrixIdx=1,nonlinearMatrices%numberOfJacobians
        NULLIFY(jacobianMatrix)
        CALL EquationsMatricesNonlinear_JacobianMatrixGet(nonlinearMatrices,matrixIdx,jacobianMatrix,err,error,*999)
        jacobianMatrix%nodalJacobian%maxNumberOfRows=0
        jacobianMatrix%nodalJacobian%maxNumberOfColumns=0
        IF(ALLOCATED(jacobianMatrix%nodalJacobian%rowDofs)) DEALLOCATE(jacobianMatrix%nodalJacobian%rowDofs)
        IF(ALLOCATED(jacobianMatrix%nodalJacobian%columnDofs)) DEALLOCATE(jacobianMatrix%nodalJacobian%columnDofs)
        IF(ALLOCATED(jacobianMatrix%nodalJacobian%matrix)) DEALLOCATE(jacobianMatrix%nodalJacobian%matrix)
      ENDDO !matrixIdx
      nonlinearMatrices%nodalResidual%maxNumberOfRows=0
      IF(ALLOCATED(nonlinearMatrices%nodalResidual%rowDofs)) DEALLOCATE(nonlinearMatrices%nodalResidual%rowDofs)
      IF(ALLOCATED(nonlinearMatrices%nodalResidual%vector)) DEALLOCATE(nonlinearMatrices%nodalResidual%vector)
    ENDIF
    rhsVector=>vectorMatrices%rhsVector
    IF(ASSOCIATED(rhsVector)) THEN
      !Finalise the nodal vector
      rhsVector%nodalVector%maxNumberOfRows=0
      IF(ALLOCATED(rhsVector%nodalVector%rowDofs)) DEALLOCATE(rhsVector%nodalVector%rowDofs)
      IF(ALLOCATED(rhsVector%nodalVector%vector)) DEALLOCATE(rhsVector%nodalVector%vector)
    ENDIF
    sourceVector=>vectorMatrices%sourceVector
    IF(ASSOCIATED(sourceVector)) THEN
      !Finalise the nodal source vector
      sourceVector%nodalVector%maxNumberOfRows=0
      IF(ALLOCATED(sourceVector%nodalVector%rowDofs)) DEALLOCATE(sourceVector%nodalVector%rowDofs)
      IF(ALLOCATED(sourceVector%nodalVector%vector)) DEALLOCATE(sourceVector%nodalVector%vector)
    ENDIF
   
    EXITS("EquationsMatrices_NodalFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NodalFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the nodal matrix.
  SUBROUTINE EquationsMatrices_NodalMatrixInitialise(nodalMatrix,err,error,*)

    !Argument variables
    TYPE(NodalMatrixType) :: nodalMatrix !The nodal matrix to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMatrices_NodalMatrixInitialise",err,error,*999)

    nodalMatrix%equationsMatrixNumber=0
    nodalMatrix%numberOfRows=0
    nodalMatrix%numberOfColumns=0
    nodalMatrix%maxNumberOfRows=0
    nodalMatrix%maxNumberOfColumns=0
       
    EXITS("EquationsMatrices_NodalMatrixInitialise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalMatrixInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NodalMatrixInitialise

  !
  !================================================================================================================================
  !

  !>Finalise an nodal matrix and deallocate all memory
  SUBROUTINE EquationsMatrices_NodalMatrixFinalise(nodalMatrix,err,error,*)

    !Argument variables
    TYPE(NodalMatrixType):: nodalMatrix !<The nodal matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("EquationsMatrices_NodalMatrixFinalise",err,error,*999)

    IF(ALLOCATED(nodalMatrix%rowDofs)) DEALLOCATE(nodalMatrix%rowDofs)
    IF(ALLOCATED(nodalMatrix%columnDofs)) DEALLOCATE(nodalMatrix%columnDofs)
    IF(ALLOCATED(nodalMatrix%matrix)) DEALLOCATE(nodalMatrix%matrix)
    
    EXITS("EquationsMatrices_NodalMatrixFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalMatrixFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NodalMatrixFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the nodal vector
  SUBROUTINE EquationsMatrices_NodalVectorInitialise(nodalVector,err,error,*)

    !Argument variables
    TYPE(NodalVectorType) :: nodalVector !The nodal vector to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMatrices_NodalVectorInitialise",err,error,*999)

    nodalVector%numberOfRows=0
    nodalVector%maxNumberOfRows=0
       
    EXITS("EquationsMatrices_NodalVectorInitialise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalVectorInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NodalVectorInitialise
  
  !
  !================================================================================================================================
  !

  !>Finalise an nodal vector and deallocate all memory
  SUBROUTINE EquationsMatrices_NodalVectorFinalise(nodalVector,err,error,*)

    !Argument variables
    TYPE(NodalVectorType):: nodalVector !<The nodal vector to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("EquationsMatrices_NodalVectorFinalise",err,error,*999)

    IF(ALLOCATED(nodalVector%rowDofs)) DEALLOCATE(nodalVector%rowDofs)
    IF(ALLOCATED(nodalVector%vector)) DEALLOCATE(nodalVector%vector)
    
    EXITS("EquationsMatrices_NodalVectorFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalVectorFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NodalVectorFinalise

  !
  !================================================================================================================================
  !

  !>Adds the Jacobian matrices into the equations Jacobian.
  SUBROUTINE EquationsMatrices_JacobianNodeAdd(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: jacobianMatrixIdx
    TYPE(EquationsJacobianType), POINTER :: jacobianMatrix
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("EquationsMatrices_JacobianNodeAdd()")
#endif

    ENTERS("EquationsMatrices_JacobianNodeAdd",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not allocated.",err,error,*999)
    
    nonlinearMatrices=>vectorMatrices%nonlinearMatrices
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      DO jacobianMatrixIdx=1,nonlinearMatrices%numberOfJacobians
        NULLIFY(jacobianMatrix)
        CALL EquationsMatricesNonlinear_JacobianMatrixGet(nonlinearMatrices,jacobianMatrixIdx,jacobianMatrix,err,error,*999)
        IF(jacobianMatrix%updateJacobian) THEN
          !Add in Jacobian element matrices
          CALL DistributedMatrix_ValuesAdd(jacobianMatrix%jacobian,jacobianMatrix%nodalJacobian%rowDofs(1: &
            & jacobianMatrix%nodalJacobian%numberOfRows),jacobianMatrix%nodalJacobian%columnDofs(1: &
            & jacobianMatrix%nodalJacobian%numberOfColumns),jacobianMatrix%nodalJacobian%matrix(1: &
            & jacobianMatrix%nodalJacobian%numberOfRows,1:jacobianMatrix%nodalJacobian%numberOfColumns), &
            & err,error,*999)
        ENDIF
      ENDDO !jacobianMatrixIdx
    ENDIF
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("EquationsMatrices_JacobianNodeAdd()")
#endif
    
    EXITS("EquationsMatrices_JacobianNodeAdd")
    RETURN
999 ERRORSEXITS("EquationsMatrices_JacobianNodeAdd",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_JacobianNodeAdd

  !
  !================================================================================================================================
  !

  !>Finalise the element calculation information and deallocate all memory
  SUBROUTINE EquationsMatrices_ElementFinalise(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<The equations matrices for which to finalise the elements
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsJacobianType), POINTER :: jacobianMatrix
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    
    ENTERS("EquationsMatrices_ElementFinalise",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)
    
    dynamicMatrices=>vectorMatrices%dynamicMatrices
    IF(ASSOCIATED(dynamicMatrices)) THEN
      !Finalise the dynamic element matrices
      DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
        NULLIFY(equationsMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,equationsMatrix,err,error,*999)
        CALL EquationsMatrices_ElementMatrixFinalise(equationsMatrix%elementMatrix,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    linearMatrices=>vectorMatrices%linearMatrices
    IF(ASSOCIATED(linearMatrices)) THEN
      !Finalise the linear element matrices
      DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
        NULLIFY(equationsMatrix)
        CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,equationsMatrix,err,error,*999)
        CALL EquationsMatrices_ElementMatrixFinalise(equationsMatrix%elementMatrix,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    nonlinearMatrices=>vectorMatrices%nonlinearMatrices
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      DO matrixIdx=1,nonlinearMatrices%numberOfJacobians
        NULLIFY(jacobianMatrix)
        CALL EquationsMatricesNonlinear_JacobianMatrixGet(nonlinearMatrices,matrixIdx,jacobianMatrix,err,error,*999)
        CALL EquationsMatrices_ElementMatrixFinalise(jacobianMatrix%elementJacobian,err,error,*999)
      ENDDO !matrixIdx
      CALL EquationsMatrices_ElementVectorFinalise(nonlinearMatrices%elementResidual,err,error,*999)
    ENDIF
    rhsVector=>vectorMatrices%rhsVector
    IF(ASSOCIATED(rhsVector)) THEN
      !Finalise the element vector
      CALL EquationsMatrices_ElementVectorFinalise(rhsVector%elementVector,err,error,*999)
    ENDIF
    sourceVector=>vectorMatrices%sourceVector
    IF(ASSOCIATED(sourceVector)) THEN
      !Finalise the element source vector
      CALL EquationsMatrices_ElementVectorFinalise(sourceVector%elementVector,err,error,*999)
    ENDIF
    
    EXITS("EquationsMatrices_ElementFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_ElementFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_ElementFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the element calculation information for the equations matrices
  SUBROUTINE EquationsMatrices_ElementInitialise(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !The equations matrices to initialise the element information for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    INTEGER(INTG) :: rowsNumberOfElements,colsNumberOfElements !Number of elements in the row and col variables whose dofs are present in the element matrix
    TYPE(EquationsJacobianType), POINTER :: jacobianMatrix
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable,colFieldVariable
    
    ENTERS("EquationsMatrices_ElementInitialise",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)
    
    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*999)

    rowsNumberOfElements=1
    colsNumberOfElements=1
    dynamicMatrices=>vectorMatrices%dynamicMatrices
    IF(ASSOCIATED(dynamicMatrices)) THEN
      !Initialise the dynamic element matrices
      NULLIFY(dynamicMapping)
      CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
        NULLIFY(equationsMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,equationsMatrix,err,error,*999)        
        fieldVariable=>dynamicMapping%equationsMatrixToVarMaps(matrixIdx)%VARIABLE
        CALL EquationsMatrices_ElementMatrixSetup(equationsMatrix%elementMatrix,fieldVariable,fieldVariable, &
          & rowsNumberOfElements,colsNumberOfElements,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    linearMatrices=>vectorMatrices%linearMatrices
    IF(ASSOCIATED(linearMatrices)) THEN
      !Initialise the linear element matrices
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
        NULLIFY(equationsMatrix)
        CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,equationsMatrix,err,error,*999)        
        fieldVariable=>linearMapping%equationsMatrixToVarMaps(matrixIdx)%VARIABLE
        CALL EquationsMatrices_ElementMatrixSetup(equationsMatrix%elementMatrix,fieldVariable,fieldVariable, &
          & rowsNumberOfElements,colsNumberOfElements,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    nonlinearMatrices=>vectorMatrices%nonlinearMatrices
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      !Initialise the Jacobian element matrices
      NULLIFY(nonlinearMapping)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
      fieldVariable=>nonlinearMapping%jacobianToVarMap(1)%VARIABLE
      DO matrixIdx=1,nonlinearMatrices%numberOfJacobians
        NULLIFY(jacobianMatrix)
        CALL EquationsMatricesNonlinear_JacobianMatrixGet(nonlinearMatrices,matrixIdx,jacobianMatrix,err,error,*999)        
        colFieldVariable=>nonlinearMapping%jacobianToVarMap(matrixIdx)%VARIABLE
        CALL EquationsMatrices_ElementMatrixSetup(jacobianMatrix%elementJacobian,fieldVariable,colFieldVariable, &
          & rowsNumberOfElements,colsNumberOfElements,err,error,*999)
      ENDDO !matrixIdx
      !Use RHS variable for residual vector, otherwise first nonlinear variable if no RHS
      rhsMapping=>vectorMapping%rhsMapping
      IF(ASSOCIATED(rhsMapping)) THEN
        fieldVariable=>rhsMapping%rhsVariable
      ELSE
        fieldVariable=>nonlinearMapping%jacobianToVarMap(1)%variable
      ENDIF
      CALL EquationsMatrices_ElementVectorSetup(nonlinearMatrices%elementResidual,fieldVariable,err,error,*999)
      nonlinearMatrices%elementResidualCalculated=0
    ENDIF
    rhsVector=>vectorMatrices%rhsVector
    IF(ASSOCIATED(rhsVector)) THEN
      !Initialise the RHS element vector
      NULLIFY(rhsMapping)
      CALL EquationsMappingVector_RHSMappingGet(vectorMapping,rhsMapping,err,error,*999)
      fieldVariable=>rhsMapping%rhsVariable
      CALL EquationsMatrices_ElementVectorSetup(rhsVector%elementVector,fieldVariable,err,error,*999)
    ENDIF
    sourceVector=>vectorMatrices%sourceVector
    IF(ASSOCIATED(sourceVector)) THEN
      !Initialise the source element vector. Note that the number of rows in the source vector is taken, for now, from the RHS
      !vector
      IF(ASSOCIATED(rhsVector)) THEN
        !Initialise the RHS element vector
        NULLIFY(rhsMapping)
        CALL EquationsMappingVector_RHSMappingGet(vectorMapping,rhsMapping,err,error,*999)
        fieldVariable=>rhsMapping%rhsVariable
        CALL EquationsMatrices_ElementVectorSetup(sourceVector%elementVector,fieldVariable,err,error,*999)
      ELSE
        CALL FlagError("Not implemented.",err,error,*999)
      ENDIF
    ENDIF
   
    EXITS("EquationsMatrices_ElementInitialise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_ElementInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_ElementInitialise

  !
  !================================================================================================================================
  !

  !>Finalise a equations matrix and deallocate all memory
  SUBROUTINE EquationsMatrices_EquationsMatrixFinalise(equationsMatrix,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("EquationsMatrices_EquationsMatrixFinalise",err,error,*999)

    IF(ASSOCIATED(equationsMatrix)) THEN
      IF(ASSOCIATED(equationsMatrix%matrix)) CALL DistributedMatrix_Destroy(equationsMatrix%matrix,err,error,*999)
      CALL EquationsMatrices_ElementMatrixFinalise(equationsMatrix%elementMatrix,err,error,*999)
      CALL EquationsMatrices_NodalMatrixFinalise(equationsMatrix%nodalMatrix,err,error,*999)
      IF(ASSOCIATED(equationsMatrix%tempVector)) CALL DistributedVector_Destroy(equationsMatrix%tempVector,err,error,*999)
    ENDIF
    
    EXITS("EquationsMatrices_EquationsMatrixFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_EquationsMatrixFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_EquationsMatrixFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the dynamic equations matrix.
  SUBROUTINE EquationsMatrices_EquationsMatrixDynamicInitialise(dynamicMatrices,matrixNumber,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices !<A pointer to the dynamic matrices to initialise the dynamic equations matrix for
    INTEGER(INTG) :: matrixNumber !<The dynamic matrix number in the dynamic equations matrices to initialise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("EquationsMatrices_EquationsMatrixDynamicInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(dynamicMatrices)) CALL FlagError("Dynamic matrices is not associated.",err,error,*998)
    IF(matrixNumber<1.OR.matrixNumber>dynamicMatrices%numberOfDynamicMatrices) THEN
      localError="The specified dynamic matrix number of "//TRIM(NumberToVString(matrixNumber,"*",err,error))// &
        & " is invalid. The matrix number must be >= 1 and <= "// &
        & TRIM(NumberToVString(dynamicMatrices%numberOfDynamicMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    IF(.NOT.ALLOCATED(dynamicMatrices%matrices)) CALL FlagError("Dynamic matrices matrices is not allocated.",err,error,*998)
    IF(ASSOCIATED(dynamicMatrices%matrices(matrixNumber)%ptr)) THEN
      localError="Equations matrix for dynamic matrix number "//TRIM(NumberToVString(matrixNumber,"*",err,error))// &
        & " is already associated."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    NULLIFY(dynamicMapping)
    CALL EquationsMatricesDynamic_DynamicMappingGet(dynamicMatrices,dynamicMapping,err,error,*999)
    
    ALLOCATE(dynamicMatrices%matrices(matrixNumber)%ptr,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations matrix.",err,error,*999)
    equationsMatrix=>dynamicMatrices%matrices(matrixNumber)%ptr
    equationsMatrix%matrixNumber=matrixNumber
    equationsMatrix%dynamicMatrices=>dynamicMatrices
    NULLIFY(equationsMatrix%linearMatrices)
    equationsMatrix%storageType=MATRIX_BLOCK_STORAGE_TYPE
    equationsMatrix%structureType=EQUATIONS_MATRIX_NO_STRUCTURE
    equationsMatrix%lumped=.FALSE.
    equationsMatrix%updateMatrix=.TRUE.
    equationsMatrix%firstAssembly=.TRUE.
    equationsMatrix%numberOfColumns=dynamicMapping%equationsMatrixToVarMaps(matrixNumber)%numberOfColumns
    dynamicMapping%equationsMatrixToVarMaps(matrixNumber)%equationsMatrix=>equationsMatrix
    NULLIFY(equationsMatrix%matrix)
    CALL EquationsMatrices_ElementMatrixInitialise(equationsMatrix%elementMatrix,err,error,*999)
    CALL EquationsMatrices_NodalMatrixInitialise(equationsMatrix%nodalMatrix,err,error,*999)
    NULLIFY(equationsMatrix%tempVector)
    
    EXITS("EquationsMatrices_EquationsMatrixDynamicInitialise")
    RETURN
999 CALL EquationsMatrices_EquationsMatrixFinalise(dynamicMatrices%matrices(matrixNumber)%ptr,dummyErr,dummyError,*998)
998 ERRORS("EquationsMatrices_EquationsMatrixDynamicInitialise",err,error)
    EXITS("EquationsMatrices_EquationsMatrixDynamicInitialise")
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_EquationsMatrixDynamicInitialise

  !
  !================================================================================================================================
  !

  !>Initialise the linear equations matrix.
  SUBROUTINE EquationsMatrices_EquationsMatrixLinearInitialise(linearMatrices,matrixNumber,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices !<A pointer to the linear matrices to initialise the linear equations matrix for
    INTEGER(INTG) :: matrixNumber !<The linear matrix number in the linear equations matrices to initialise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("EquationsMatrices_EquationsMatrixLinearInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(linearMatrices)) CALL FlagError("Linear matrices is not associated.",err,error,*998)    
    IF(matrixNumber<1.OR.matrixNumber>linearMatrices%numberOfLinearMatrices) THEN
      localError="The specified linear matrix number of "//TRIM(NumberToVString(matrixNumber,"*",err,error))// &
        & " is invalid. The matrix number must be >= 1 and <= "// &
        & TRIM(NumberToVString(linearMatrices%numberOfLinearMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    IF(.NOT.ALLOCATED(linearMatrices%matrices)) CALL FlagError("Linear matrices matrices is not allocated.",err,error,*998)
    IF(ASSOCIATED(linearMatrices%matrices(matrixNumber)%ptr)) THEN
      localError="Equations matrix for linear matrix number "//TRIM(NumberToVString(matrixNumber,"*",err,error))// &
        & " is already associated."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    NULLIFY(linearMapping)
    CALL EquationsMatricesLinear_LinearMappingGet(linearMatrices,linearMapping,err,error,*999)
    
    ALLOCATE(linearMatrices%matrices(matrixNumber)%ptr,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations matrix.",err,error,*999)
    equationsMatrix=>linearMatrices%matrices(matrixNumber)%ptr
    equationsMatrix%matrixNumber=matrixNumber
    NULLIFY(equationsMatrix%dynamicMatrices)
    equationsMatrix%linearMatrices=>linearMatrices
    equationsMatrix%storageType=MATRIX_BLOCK_STORAGE_TYPE
    equationsMatrix%structureType=EQUATIONS_MATRIX_NO_STRUCTURE
    equationsMatrix%lumped=.FALSE.
    equationsMatrix%updateMatrix=.TRUE.
    equationsMatrix%firstAssembly=.TRUE.
    equationsMatrix%numberOfColumns=linearMapping%equationsMatrixToVarMaps(matrixNumber)%numberOfColumns
    linearMapping%equationsMatrixToVarMaps(matrixNumber)%equationsMatrix=>equationsMatrix
    NULLIFY(equationsMatrix%matrix)
    CALL EquationsMatrices_ElementMatrixInitialise(equationsMatrix%elementMatrix,err,error,*999)
    CALL EquationsMatrices_NodalMatrixInitialise(equationsMatrix%nodalMatrix,err,error,*999)
    NULLIFY(equationsMatrix%tempVector)
   
    EXITS("EquationsMatrices_EquationsMatrixLinearInitialise")
    RETURN
999 CALL EquationsMatrices_EquationsMatrixFinalise(linearMatrices%matrices(matrixNumber)%ptr,dummyErr,dummyError,*998)
998 ERRORS("EquationsMatrices_EquationsMatrixLinearInitialise",err,error)
    EXITS("EquationsMatrices_EquationsMatrixLinearInitialise")
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_EquationsMatrixLinearInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the vector equations matrices dynamic matrices and deallocates all memory
  SUBROUTINE EquationsMatrices_DynamicFinalise(dynamicMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices !<A pointer to the vector equation matrices dynamic matrices to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
     
    ENTERS("EquationsMatrices_DynamicFinalise",err,error,*999)

    IF(ASSOCIATED(dynamicMatrices)) THEN
      IF(ALLOCATED(dynamicMatrices%matrices)) THEN
        DO matrixIdx=1,SIZE(dynamicMatrices%matrices,1)
          CALL EquationsMatrices_EquationsMatrixFinalise(dynamicMatrices%matrices(matrixIdx)%ptr,err,error,*999)
        ENDDO !matrixIdx
        DEALLOCATE(dynamicMatrices%matrices)
      ENDIF
      IF(ASSOCIATED(dynamicMatrices%tempVector)) CALL DistributedVector_Destroy(dynamicMatrices%tempVector,err,error,*999)
      DEALLOCATE(dynamicMatrices)
    ENDIF
    
    EXITS("EquationsMatrices_DynamicFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_DynamicFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_DynamicFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises the vector equations matrices dynamic matrices
  SUBROUTINE EquationsMatrices_DynamicInitialise(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equation matrices to initialise the dynamic matrices for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,matrixIdx
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(VARYING_STRING) :: dummyError
     
    ENTERS("EquationsMatrices_DynamicInitialise",err,error,*998)
    
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*998)
    IF(ASSOCIATED(vectorMatrices%dynamicMatrices)) &
      & CALL FlagError("Equations matrices dynamic matrices is already associated.",err,error,*998)
    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*998)
    dynamicMapping=>vectorMapping%dynamicMapping

    IF(ASSOCIATED(dynamicMapping)) THEN
      ALLOCATE(vectorMatrices%dynamicMatrices,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations matrices dynamic matrices.",err,error,*999)
      vectorMatrices%dynamicMatrices%vectorMatrices=>vectorMatrices
      vectorMatrices%dynamicMatrices%numberOfDynamicMatrices=dynamicMapping%numberOfDynamicMatrices
      ALLOCATE(vectorMatrices%dynamicMatrices%matrices(dynamicMapping%numberOfDynamicMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations matrices dynamic matrices matrices.",err,error,*999)
      DO matrixIdx=1,dynamicMapping%numberOfDynamicMatrices
        NULLIFY(vectorMatrices%dynamicMatrices%matrices(matrixIdx)%ptr)
        CALL EquationsMatrices_EquationsMatrixDynamicInitialise(vectorMatrices%dynamicMatrices,matrixIdx,err,error,*999)
      ENDDO !matrixIdx
      NULLIFY(vectorMatrices%dynamicMatrices%tempVector)
    ENDIF
    
    EXITS("EquationsMatrices_DynamicInitialise")
    RETURN
999 CALL EquationsMatrices_DynamicFinalise(vectorMatrices%dynamicMatrices,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_DynamicInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_DynamicInitialise
  
  !
  !================================================================================================================================
  !

  !>Adds the Hessian elmental matrices into the equations Hessian.
  SUBROUTINE EquationsMatrices_HessianElementAdd(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: hessianMatrixIdx
    TYPE(EquationsHessianType), POINTER :: hessianMatrix
    TYPE(EquationsMatricesOptimisationType), POINTER :: optimisationMatrices
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMatrices_HessianElementAdd",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Equations matrices is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(vectorMatrices%optimisationMatrices)) &
      & CALL FlagError("Equations matrices optimisation matrices is not associated.",err,error,*999)

    optimisationMatrices=>vectorMatrices%optimisationMatrices
    DO hessianMatrixIdx=1,optimisationMatrices%numberOfHessians
      hessianMatrix=>optimisationMatrices%hessians(hessianMatrixIdx)%ptr
      IF(ASSOCIATED(hessianMatrix)) THEN
        IF(hessianMatrix%updateHessian) THEN
          !Add in Hessian element matrices
          CALL DistributedMatrix_ValuesAdd(hessianMatrix%hessian,hessianMatrix%elementHessian%rowDOFS( &
            & 1:hessianMatrix%elementHessian%numberOfRows),hessianMatrix%elementHessian%columnDOFS( &
            & 1:hessianMatrix%elementHessian%numberOfColumns),hessianMatrix%elementHessian%matrix(1: &
            & hessianMatrix%elementHessian%numberOfRows,1:hessianMatrix%elementHessian%numberOfColumns), &
            & err,error,*999)
        ENDIF
      ELSE
        localError="The Hessian matrix for Hessian matrix index "//TRIM(NumberToVString(hessianMatrixIdx,"*",err,error))// &
          & " is not associated."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDDO !hessianMatrixIdx
    
    EXITS("EquationsMatrices_HessianElementAdd")
    RETURN
999 ERRORSEXITS("EquationsMatrices_HessianElementAdd",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_HessianElementAdd

  !
  !================================================================================================================================
  !

  !>Outputs the equations Hessian matrices
  SUBROUTINE EquationsMatrices_HessianOutput(id,vectorMatrices,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the ouptut stream
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations Hessian matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: hessianMatrixIdx
    TYPE(EquationsHessianType), POINTER :: hessianMatrix
    TYPE(EquationsMatricesOptimisationType), POINTER :: optimisationMatrices
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsMatrices_HessianOutput",err,error,*999)
    
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Equations matrices is not associated.",err,error,*999)
    IF(.NOT.vectorMatrices%vectorMatricesFinished) &
      & CALL FlagError("Equations matrices have not been finished.",err,error,*999)
    IF(.NOT.ASSOCIATED(vectorMatrices%optimisationMatrices)) &
      & CALL FlagError("Equations matrices optimisation matrices is not associated.",err,error,*999)

    optimisationMatrices=>vectorMatrices%optimisationMatrices
    
    CALL WriteString(id,"",err,error,*999)
    CALL WriteString(id,"Hessian matrices:",err,error,*999)
    DO hessianMatrixIdx=1,optimisationMatrices%numberOfHessians
      hessianMatrix=>optimisationMatrices%hessians(hessianMatrixIdx)%ptr
      IF(ASSOCIATED(hessianMatrix)) THEN        
        CALL WriteStringValue(id,"Hessian matrix: ",hessianMatrixIdx,err,error,*999)
        CALL DistributedMatrix_Output(id,hessianMatrix%hessian,err,error,*999)
      ELSE
        localError="The Hessian matrix for Hessian matrix index "//TRIM(NumberToVString(hessianMatrixIdx,"*",err,error))// &
          & " is not associated."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDDO !hessianMatrixIdx
    
    EXITS("EquationsMatrices_HessianOutput")
    RETURN
999 ERRORSEXITS("EquationsMatrices_HessianOutput",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_HessianOutput
  
  !
  !================================================================================================================================
  !

  !>Adds the Jacobain matrices into the equations Jacobian.
  SUBROUTINE EquationsMatrices_JacobianElementAdd(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: jacobianMatrixIdx
    TYPE(EquationsJacobianType), POINTER :: jacobianMatrix
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("EquationsMatrices_JacobianElementAdd()")
#endif

    ENTERS("EquationsMatrices_JacobianElementAdd",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)

    DO jacobianMatrixIdx=1,nonlinearMatrices%numberOfJacobians
      NULLIFY(jacobianMatrix)
      CALL EquationsMatricesNonlinear_JacobianMatrixGet(nonlinearMatrices,jacobianMatrixIdx,jacobianMatrix,err,error,*999)
      IF(jacobianMatrix%updateJacobian) THEN
        !Add in Jacobian element matrices
        CALL DistributedMatrix_ValuesAdd(jacobianMatrix%jacobian,jacobianMatrix%elementJacobian%rowDOFS(1: &
          & jacobianMatrix%elementJacobian%numberOfRows),jacobianMatrix%elementJacobian%columnDOFS(1: &
          & jacobianMatrix%elementJacobian%numberOfColumns),jacobianMatrix%elementJacobian%matrix(1: &
          & jacobianMatrix%elementJacobian%numberOfRows,1:jacobianMatrix%elementJacobian%numberOfColumns), &
          & err,error,*999)
      ENDIF
    ENDDO !jacobianMatrixIdx
      
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("EquationsMatrices_JacobianElementAdd()")
#endif
    
    EXITS("EquationsMatrices_JacobianElementAdd")
    RETURN
999 ERRORSEXITS("EquationsMatrices_JacobianElementAdd",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_JacobianElementAdd

  !
  !================================================================================================================================
  !

  !>Outputs the equations Jacobian matrices
  SUBROUTINE EquationsMatrices_JacobianOutput(id,vectorMatrices,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the ouptut stream
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations Jacobian matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: jacobianMatrixIdx
    TYPE(EquationsJacobianType), POINTER :: jacobianMatrix
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    
    ENTERS("EquationsMatrices_JacobianOutput",err,error,*999)    
    
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)    
    IF(.NOT.vectorMatrices%vectorMatricesFinished) &
      & CALL FlagError("Vector equations matrices have not been finished.",err,error,*999)
    
    nonlinearMatrices=>vectorMatrices%nonlinearMatrices
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      CALL WriteString(id,"",err,error,*999)
      CALL WriteString(id,"Jacobian matrices:",err,error,*999)
      CALL WriteStringValue(id,"Number of Jacobian matrices = ",nonlinearMatrices%numberOfJacobians,err,error,*999)
      DO jacobianMatrixIdx=1,nonlinearMatrices%numberOfJacobians
        NULLIFY(jacobianMatrix)
        CALL EquationsMatricesNonlinear_JacobianMatrixGet(nonlinearMatrices,jacobianMatrixIdx,jacobianMatrix,err,error,*999)
        CALL WriteStringValue(id,"Jacobian matrix: ",jacobianMatrixIdx,err,error,*999)
        CALL DistributedMatrix_Output(id,jacobianMatrix%jacobian,err,error,*999)
      ENDDO !jacobianMatrixIdx
    ENDIF
    
    EXITS("EquationsMatrices_JacobianOutput")
    RETURN
999 ERRORSEXITS("EquationsMatrices_JacobianOutput",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_JacobianOutput
  
  !
  !================================================================================================================================
  !

  !>Sets the Jacobian calculation types of the residual variables
  SUBROUTINE EquationsMatrices_JacobianTypesSet(vectorMatrices,jacobianTypes,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices to set the Jacobian type for.
    INTEGER(INTG), INTENT(IN) :: jacobianTypes(:) !<jacobianTypes(matrixIdx). The Jacobian calculation type for the matrixIdx'th Jacobian matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsJacobianType), POINTER :: jacobianMatrix
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    INTEGER(INTG) :: matrixIdx,numberOfjacobians,jacobianType
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMatrices_JacobianTypesSet",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)
    IF(vectorMatrices%vectorMatricesFinished) CALL FlagError("Vector equations matrices have been finished.",err,error,*999)

    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
    numberOfJacobians=SIZE(jacobianTypes,1)
    IF(numberOfJacobians/=nonlinearMatrices%numberOfJacobians) THEN
      localError="Invalid number of Jacobian calculation types. The number of types "// &
        & TRIM(NumberToVString(numberOfJacobians,"*",err,error))//" should be "// &
        & TRIM(NumberToVString(nonlinearMatrices%numberOfJacobians,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
     
    DO matrixIdx=1,numberOfJacobians
      jacobianType=jacobianTypes(matrixIdx)
      NULLIFY(jacobianMatrix)
      CALL EquationsMatricesNonlinear_JacobianMatrixGet(nonlinearMatrices,matrixIdx,jacobianMatrix,err,error,*999)
      SELECT CASE(jacobianType)
      CASE(EQUATIONS_JACOBIAN_FINITE_DIFFERENCE_CALCULATED,EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED)
        jacobianMatrix%jacobianCalculationType=jacobianType
      CASE DEFAULT
        localError="The Jacobian calculation type of "//TRIM(NumberToVString(jacobianType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !matrixIdx

    EXITS("EquationsMatrices_JacobianTypesSet")
    RETURN
999 ERRORSEXITS("EquationsMatrices_JacobianTypesSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_JacobianTypesSet

  !
  !================================================================================================================================
  !

  !>Finalises the vector equations matrices linear matrices and deallocates all memory
  SUBROUTINE EquationsMatrices_LinearFinalise(linearMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices !<A pointer to the vector equation matrices linear matrices to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
     
    ENTERS("EquationsMatrices_LinearFinalise",err,error,*999)

    IF(ASSOCIATED(linearMatrices)) THEN
      IF(ALLOCATED(linearMatrices%matrices)) THEN
        DO matrixIdx=1,SIZE(linearMatrices%matrices,1)
          CALL EquationsMatrices_EquationsMatrixFinalise(linearMatrices%matrices(matrixIdx)%ptr,err,error,*999)
        ENDDO !matrixIdx
        DEALLOCATE(linearMatrices%matrices)
      ENDIF
      DEALLOCATE(linearMatrices)
    ENDIF
    
    EXITS("EquationsMatrices_LinearFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_LinearFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_LinearFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises the equations matrices linear matrices
  SUBROUTINE EquationsMatrices_LinearInitialise(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equation matrices to initialise the linear matrices for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,matrixIdx
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(VARYING_STRING) :: dummyError
     
    ENTERS("EquationsMatrices_LinearInitialise",err,error,*998)
    
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*998)
    IF(ASSOCIATED(vectorMatrices%linearMatrices)) &
      & CALL FlagError("Vector equations matrices linear matrices is already associated.",err,error,*998)
    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*999)
    linearMapping=>vectorMapping%linearMapping

    IF(ASSOCIATED(linearMapping)) THEN
      ALLOCATE(vectorMatrices%linearMatrices,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations matrices linear matrices.",err,error,*999)
      vectorMatrices%linearMatrices%vectorMatrices=>vectorMatrices
      vectorMatrices%linearMatrices%numberOfLinearMatrices=linearMapping%numberOfLinearMatrices
      ALLOCATE(vectorMatrices%linearMatrices%matrices(linearMapping%numberOfLinearMatrices),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations matrices linear matrices matrices.",err,error,*999)
      DO matrixIdx=1,linearMapping%numberOfLinearMatrices
        NULLIFY(vectorMatrices%linearMatrices%matrices(matrixIdx)%ptr)
        CALL EquationsMatrices_EquationsMatrixLinearInitialise(vectorMatrices%linearMatrices,matrixIdx,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    
    EXITS("EquationsMatrices_LinearInitialise")
    RETURN
999 CALL EquationsMatrices_LinearFinalise(vectorMatrices%linearMatrices,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_LinearInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_LinearInitialise
  
  !
  !================================================================================================================================
  !

  !>Finalises the equations matrices nonlinear matrices and deallocates all memory
  SUBROUTINE EquationsMatrices_NonlinearFinalise(nonlinearMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices !<A pointer to the equation matrices nonlinear matrices to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
     
    ENTERS("EquationsMatrices_NonlinearFinalise",err,error,*999)

    IF(ASSOCIATED(nonlinearMatrices)) THEN
      IF(ALLOCATED(nonlinearMatrices%jacobians)) THEN
        DO matrixIdx=1,nonlinearMatrices%numberOfJacobians
          CALL EquationsMatrices_JacobianFinalise(nonlinearMatrices%jacobians(matrixIdx)%ptr,err,error,*999)
        ENDDO !matrixidx
        DEALLOCATE(nonlinearMatrices%jacobians)
      ENDIF
      IF(ASSOCIATED(nonlinearMatrices%residual)) CALL DistributedVector_Destroy(nonlinearMatrices%residual,err,error,*999)
      CALL EquationsMatrices_ElementVectorFinalise(nonlinearMatrices%elementResidual,err,error,*999)
      CALL EquationsMatrices_NodalVectorFinalise(nonlinearMatrices%nodalResidual,err,error,*999)
      DEALLOCATE(nonlinearMatrices)
    ENDIF
    
    EXITS("EquationsMatrices_NonlinearFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NonlinearFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NonlinearFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises the equations matrices nonlinear matrices
  SUBROUTINE EquationsMatrices_NonlinearInitialise(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equation matrices to initialise the nonlinear matrices for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,dummyErr
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("EquationsMatrices_NonlinearInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*998)
    IF(ASSOCIATED(vectorMatrices%nonlinearMatrices)) &
      & CALL FlagError("Vector equations matrices nonlinear matrices is already associated.",err,error,*998)

    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*999)
   
    nonlinearMapping=>vectorMapping%nonlinearMapping
    IF(ASSOCIATED(nonlinearMapping)) THEN
      ALLOCATE(vectorMatrices%nonlinearMatrices,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations matrices nonlinear matrices.",err,error,*999)
      vectorMatrices%nonlinearMatrices%vectorMatrices=>vectorMatrices
      vectorMatrices%nonlinearMatrices%updateResidual=.TRUE.
      vectorMatrices%nonlinearMatrices%firstAssembly=.TRUE.
      NULLIFY(vectorMatrices%nonlinearMatrices%residual)
      CALL EquationsMatrices_ElementVectorInitialise(vectorMatrices%nonlinearMatrices%elementResidual,err,error,*999)
      CALL EquationsMatrices_NodalVectorInitialise(vectorMatrices%nonlinearMatrices%nodalResidual,err,error,*999)
      vectorMatrices%nonlinearMatrices%numberOfJacobians=nonlinearMapping%numberOfResidualVariables
      ALLOCATE(vectorMatrices%nonlinearMatrices%jacobians(nonlinearMapping%numberOfResidualVariables),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations matrices Jacobian matrices.",err,error,*999)
      DO matrixIdx=1,nonlinearMapping%numberOfResidualVariables
        NULLIFY(vectorMatrices%nonlinearMatrices%jacobians(matrixIdx)%ptr)
        CALL EquationsMatrices_JacobianInitialise(vectorMatrices%nonlinearMatrices,matrixIdx,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    
    EXITS("EquationsMatrices_NonlinearInitialise")
    RETURN
999 CALL EquationsMatrices_NonlinearFinalise(vectorMatrices%nonlinearMatrices,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_NonlinearInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NonlinearInitialise
  
  !
  !================================================================================================================================
  !

  !>Finalises the equations matrices optimisation matrices and deallocates all memory
  SUBROUTINE EquationsMatrices_OptimisationFinalise(optimisationMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesOptimisationType), POINTER :: optimisationMatrices !<A pointer to the equation matrices optimisation matrices to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
     
    ENTERS("EquationsMatrices_OptimisationFinalise",err,error,*999)

    IF(ASSOCIATED(optimisationMatrices)) THEN
      IF(ALLOCATED(optimisationMatrices%hessians)) THEN
        DO matrixIdx=1,SIZE(optimisationMatrices%hessians,1)
          CALL EquationsMatrices_HessianFinalise(optimisationMatrices%hessians(matrixIdx)%ptr,err,error,*999)
        ENDDO !matrixIdx
        DEALLOCATE(optimisationMatrices%hessians)
      ENDIF
      IF(ASSOCIATED(optimisationMatrices%gradient)) CALL DistributedVector_Destroy(optimisationMatrices%gradient,err,error,*999)
      CALL EquationsMatrices_ElementVectorFinalise(optimisationMatrices%elementGradient,err,error,*999)
      IF(ASSOCIATED(optimisationMatrices%constraints)) CALL DistributedVector_Destroy(optimisationMatrices%constraints, &
        & err,error,*999)
      CALL EquationsMatrices_ElementVectorFinalise(optimisationMatrices%elementConstraints,err,error,*999)
      IF(ASSOCIATED(optimisationMatrices%lowerBounds)) CALL DistributedVector_Destroy(optimisationMatrices%lowerBounds, &
        & err,error,*999)
      CALL EquationsMatrices_ElementVectorFinalise(optimisationMatrices%elementLowerBounds,err,error,*999)
      IF(ASSOCIATED(optimisationMatrices%upperBounds)) CALL DistributedVector_Destroy(optimisationMatrices%upperBounds, &
        & err,error,*999)
      CALL EquationsMatrices_ElementVectorFinalise(optimisationMatrices%elementUpperBounds,err,error,*999)
      IF(ASSOCIATED(optimisationMatrices%residual)) CALL DistributedVector_Destroy(optimisationMatrices%residual,err,error,*999)
      CALL EquationsMatrices_ElementVectorFinalise(optimisationMatrices%elementResidual,err,error,*999)
      DEALLOCATE(optimisationMatrices)
    ENDIF
    
    EXITS("EquationsMatrices_OptimisationFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_OptimisationFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_OptimisationFinalise
  
  !
  !================================================================================================================================
  !

  !>Outputs the equations matrices
  SUBROUTINE EquationsMatrices_VectorOutput(id,vectorMatrices,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the ouptut stream
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    
    ENTERS("EquationsMatrices_VectorOutput",err,error,*999)
    
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)
    IF(.NOT.vectorMatrices%vectorMatricesFinished) &
      & CALL FlagError("Vector equations matrices have not been finished.",err,error,*999)
    
    CALL WriteString(id,"",err,error,*999)
    CALL WriteString(id,"Equations matrices:",err,error,*999)
    dynamicMatrices=>vectorMatrices%dynamicMatrices
    IF(ASSOCIATED(dynamicMatrices)) THEN
      CALL WriteString(id,"Dynamic matrices:",err,error,*999)
      CALL WriteStringValue(id,"Number of dynamic matrices = ",dynamicMatrices%numberOfDynamicMatrices,err,error,*999)
      DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
        NULLIFY(equationsMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,equationsMatrix,err,error,*999)
        CALL WriteStringValue(id,"Dynamic matrix: ",matrixIdx,err,error,*999)
        CALL DistributedMatrix_Output(id,equationsMatrix%matrix,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    linearMatrices=>vectorMatrices%linearMatrices
    IF(ASSOCIATED(linearMatrices)) THEN
      CALL WriteString(id,"Linear matrices:",err,error,*999)
      CALL WriteStringValue(id,"Number of linear matrices = ",linearMatrices%numberOfLinearMatrices,err,error,*999)
      DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
        NULLIFY(equationsMatrix)
        CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,equationsMatrix,err,error,*999)
        CALL WriteStringValue(id,"Linear matrix: ",matrixIdx,err,error,*999)
        CALL DistributedMatrix_Output(id,equationsMatrix%matrix,err,error,*999)
      ENDDO !matrixIdx
    ENDIF
    nonlinearMatrices=>vectorMatrices%nonlinearMatrices
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      CALL WriteString(id,"Nonlinear vectors:",err,error,*999)
      IF(.NOT.ASSOCIATED(nonlinearMatrices%residual)) &
        & CALL FlagError("Nonlinear matrices residual is not associated.",err,error,*999)
      CALL WriteString(id,"Residual vector:",err,error,*999)
      CALL DistributedVector_Output(id,nonlinearMatrices%residual,err,error,*999)
    ENDIF
    rhsVector=>vectorMatrices%rhsVector
    IF(ASSOCIATED(rhsVector)) THEN
      CALL WriteString(id,"RHS vector:",err,error,*999)
      CALL DistributedVector_Output(id,rhsVector%vector,err,error,*999)
    ENDIF
    sourceVector=>vectorMatrices%sourceVector
    IF(ASSOCIATED(sourceVector)) THEN
      CALL WriteString(id,"Source vector:",err,error,*999)
      CALL DistributedVector_Output(id,sourceVector%vector,err,error,*999)
    ENDIF
    
    EXITS("EquationsMatrices_Output")
    RETURN
999 ERRORSEXITS("EquationsMatrices_Output",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_VectorOutput
  
  !
  !================================================================================================================================
  !

  !>Finalises the vector equations matrices RHS vector and deallocates all memory
  SUBROUTINE EquationsMatrices_RHSFinalise(rhsVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector !<A pointer to the vector equation matrices RHS vector to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
     
    ENTERS("EquationsMatrices_RHSFinalise",err,error,*999)

    IF(ASSOCIATED(rhsVector)) THEN
      IF(ASSOCIATED(rhsVector%vector)) CALL DistributedVector_Destroy(rhsVector%vector,err,error,*999)
      CALL EquationsMatrices_ElementVectorFinalise(rhsVector%elementVector,err,error,*999)
      CALL EquationsMatrices_NodalVectorFinalise(rhsVector%nodalVector,err,error,*999)
      DEALLOCATE(rhsVector)
    ENDIF      
     
    EXITS("EquationsMatrices_RHSFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_RHSFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_RHSFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises the vector equations matrices RHS vector
  SUBROUTINE EquationsMatrices_RHSInitialise(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equation matrices to initialise the rhs vector for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("EquationsMatrices_RHSInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*998)
    IF(ASSOCIATED(vectorMatrices%rhsVector)) &
      & CALL FlagError("Vector equations matrices RHS vector is already associated.",err,error,*998)
    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*998)
    
    rhsMapping=>vectorMapping%rhsMapping
    IF(ASSOCIATED(rhsMapping)) THEN
      ALLOCATE(vectorMatrices%rhsVector,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations matrices RHS vector.",err,error,*999)
      vectorMatrices%rhsVector%updateVector=.TRUE.
      vectorMatrices%rhsVector%firstAssembly=.TRUE.
      NULLIFY(vectorMatrices%rhsVector%vector)
      CALL EquationsMatrices_ElementVectorInitialise(vectorMatrices%rhsVector%elementVector,err,error,*999)
      CALL EquationsMatrices_NodalVectorInitialise(vectorMatrices%rhsVector%nodalVector,err,error,*999)
    ENDIF
    
    EXITS("EquationsMatrices_RHSInitialise")
    RETURN
999 CALL EquationsMatrices_RHSFinalise(vectorMatrices%rhsVector,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_RHSInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_RHSInitialise
  
   !
  !================================================================================================================================
  !

  !>Destroy the scalar equations matrices
  SUBROUTINE EquationsMatrices_ScalarDestroy(scalarMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesScalarType), POINTER :: scalarMatrices !<A pointer the scalar equations matrices to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMatrices_ScalarDestroy",err,error,*999)

    IF(.NOT.ASSOCIATED(scalarMatrices)) CALL FlagError("Scalar equations matrices is not associated",err,error,*999)
    
    CALL EquationsMatrices_ScalarFinalise(scalarMatrices,err,error,*999)
        
    EXITS("EquationsMatrices_ScalarDestroy")
    RETURN
999 ERRORSEXITS("EquationsMatrices_ScalarDestroy",err,error)    
    RETURN 1
   
  END SUBROUTINE EquationsMatrices_ScalarDestroy

  !
  !================================================================================================================================
  !

  !>Finalises the vector equations matrices source vector and deallocates all memory
  SUBROUTINE EquationsMatrices_SourceFinalise(sourceVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<A pointer to the vector equation matrices source vector to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
     
    ENTERS("EquationsMatrices_SourceFinalise",err,error,*999)

    IF(ASSOCIATED(sourceVector)) THEN
      IF(ASSOCIATED(sourceVector%vector)) CALL DistributedVector_Destroy(sourceVector%vector,err,error,*999)
      CALL EquationsMatrices_ElementVectorFinalise(sourceVector%elementVector,err,error,*999)
      CALL EquationsMatrices_NodalVectorFinalise(sourceVector%nodalVector,err,error,*999)
      DEALLOCATE(sourceVector)
    ENDIF      
     
    EXITS("EquationsMatrices_SourceFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_SourceFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_SourceFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises the equations matrices source vector
  SUBROUTINE EquationsMatrices_SourceInitialise(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equation matrices to initialise the source vector for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("EquationsMatrices_SourceInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*998)
    IF(ASSOCIATED(vectorMatrices%sourceVector)) &
      & CALL FlagError("Vector equations matrices source vector is already associated.",err,error,*998)
    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*998)
    
    sourceMapping=>vectorMapping%sourceMapping
    IF(ASSOCIATED(sourceMapping)) THEN
      ALLOCATE(vectorMatrices%sourceVector,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate equations matrices source vector.",err,error,*999)
      vectorMatrices%sourceVector%updateVector=.TRUE.
      vectorMatrices%sourceVector%firstAssembly=.TRUE.
      NULLIFY(vectorMatrices%sourceVector%vector)
      CALL EquationsMatrices_ElementVectorInitialise(vectorMatrices%sourceVector%elementVector,err,error,*999)
      CALL EquationsMatrices_NodalVectorInitialise(vectorMatrices%sourceVector%nodalVector,err,error,*999)
    ENDIF
    
    EXITS("EquationsMatrices_SourceInitialise")
    RETURN
999 CALL EquationsMatrices_SourceFinalise(vectorMatrices%sourceVector,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_SourceInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_SourceInitialise
  
  !
  !================================================================================================================================
  !

  !>Sets the lumping of the dynamic equations matrices
  SUBROUTINE EquationsMatrices_DynamicLumpingTypeSet(vectorMatrices,lumpingType,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector eqautions matrices
    INTEGER(INTG), INTENT(IN) :: lumpingType(:) !<lumpingType(matrixIdx). The lumping type for the matrixIdx'th dynamic equations matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsMatrices_DynamicLumpingTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)
    IF(vectorMatrices%vectorMatricesFinished) CALL FlagError("Vector equations matrices have already been finished.",err,error,*999)
    NULLIFY(dynamicMatrices)
    CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
    IF(SIZE(lumpingType,1)/=dynamicMatrices%numberOfDynamicMatrices) THEN
      localError="The size of the lumping type array ("//TRIM(NumberToVString(SIZE(lumpingType,1),"*",err,error))// &
        & ") is not equal to the number of dynamic matrices ("// &
        & TRIM(NumberToVString(dynamicMatrices%numberOfDynamicMatrices,"*",err,error))//")."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
      NULLIFY(equationsMatrix)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,equationsMatrix,err,error,*999)
      SELECT CASE(lumpingType(matrixIdx))
      CASE(EQUATIONS_MATRIX_UNLUMPED)
        equationsMatrix%lumped=.FALSE.
      CASE(EQUATIONS_MATRIX_LUMPED)
        equationsMatrix%lumped=.TRUE.        
      CASE DEFAULT
        localError="The specified lumping type of "//TRIM(NumberToVString(lumpingType(matrixIdx),"*",err,error))// &
          & " for dynamic matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !matrixIdx
    
    EXITS("EquationsMatrices_DynamicLumpingTypeSet")
    RETURN
999 ERRORSEXITS("EquationsMatrices_DynamicLumpingTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_DynamicLumpingTypeSet

  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of the dynamic equations matrices
  SUBROUTINE EquationsMatrices_DynamicStorageTypeSet(vectorMatrices,storageType,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(IN) :: storageType(:) !<storageType(matrixIdx). The storage type for the matrixIdx'th dynamic equations matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsMatrices_DynamicStorageTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)
    IF(vectorMatrices%vectorMatricesFinished) CALL FlagError("Vector equations matrices have already been finished.",err,error,*999)
    NULLIFY(dynamicMatrices)
    CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
    IF(SIZE(storageType,1)/=dynamicMatrices%numberOfDynamicMatrices) THEN
      localError="The size of the storage type array ("//TRIM(NumberToVString(SIZE(storageType,1),"*",err,error))// &
        & ") is not equal to the number of dynamic matrices ("// &
        & TRIM(NumberToVString(dynamicMatrices%numberOfDynamicMatrices,"*",err,error))//")."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
      NULLIFY(equationsMatrix)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,equationsMatrix,err,error,*999)
      SELECT CASE(storageType(matrixIdx))
      CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
        equationsMatrix%storageType=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
        equationsMatrix%storageType=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE        
      CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
        equationsMatrix%storageType=DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
        equationsMatrix%storageType=DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        equationsMatrix%storageType=DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
        equationsMatrix%storageType=DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
        equationsMatrix%storageType=DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE
      CASE DEFAULT
        localError="The specified storage type of "//TRIM(NumberToVString(storageType(matrixIdx),"*",err,error))// &
          & " for the dynamic matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !matrixIdx
    
    EXITS("EquationsMatrices_DynamicStorageTypeSet")
    RETURN
999 ERRORSEXITS("EquationsMatrices_DynamicStorageTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_DynamicStorageTypeSet

  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of the linear equations matrices
  SUBROUTINE EquationsMatrices_LinearStorageTypeSet(vectorMatrices,storageType,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(IN) :: storageType(:) !<storageType(matrixIdx). The storage type for the matrixIdx'th linear equations matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsMatrices_LinearStorageTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)
    IF(vectorMatrices%vectorMatricesFinished) CALL FlagError("Vector equations matrices have been finished.",err,error,*999)
    NULLIFY(linearMatrices)
    CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
    IF(SIZE(storageType,1)/=linearMatrices%numberOfLinearMatrices) THEN
      localError="The size of the storage type array ("//TRIM(NumberToVString(SIZE(storageType,1),"*",err,error))// &
        & ") is not equal to the number of linear matrices ("// &
        & TRIM(NumberToVString(linearMatrices%numberOfLinearMatrices,"*",err,error))//")."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
      NULLIFY(equationsMatrix)
      CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,equationsMatrix,err,error,*999)
      SELECT CASE(storageType(matrixIdx))
      CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
        equationsMatrix%storageType=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
        equationsMatrix%storageType=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE        
      CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
        equationsMatrix%storageType=DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
        equationsMatrix%storageType=DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        equationsMatrix%storageType=DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
        equationsMatrix%storageType=DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
        equationsMatrix%storageType=DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE
      CASE DEFAULT
        localError="The specified storage type of "//TRIM(NumberToVString(storageType(matrixIdx),"*",err,error))// &
          & " for the linear matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !matrixIdx
    
    EXITS("EquationsMatrices_LinearStorageTypeSet")
    RETURN
999 ERRORSEXITS("EquationsMatrices_LinearStorageTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_LinearStorageTypeSet

  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of the nonlinear (Jacobian) equations matrices
  SUBROUTINE EquationsMatrices_NonlinearStorageTypeSet0(vectorMatrices,storageType,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector eqautions matrices
    INTEGER(INTG), INTENT(IN) :: storageType(:) !<storageType(matrixIdx). The storage type for the matrixIdx'th Jacobian equations matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsJacobianType), POINTER :: jacobianMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsMatrices_NonlinearStorageTypeSet0",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)
    IF(vectorMatrices%vectorMatricesFinished) CALL FlagError("Vector equations matrices have been finished.",err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
    IF(SIZE(storageType,1)/=nonlinearMatrices%numberOfJacobians) THEN
      localError="The size of the storage type array ("//TRIM(NumberToVString(SIZE(storageType,1),"*",err,error))// &
        & ") is not equal to the number of Jacobian matrices ("// &
        & TRIM(NumberToVString(nonlinearMatrices%numberOfJacobians,"*",err,error))//")."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO matrixIdx=1,nonlinearMatrices%numberOfJacobians
      NULLIFY(jacobianMatrix)
      CALL EquationsMatricesNonlinear_JacobianMatrixGet(nonlinearMatrices,matrixIdx,jacobianMatrix,err,error,*999)
      SELECT CASE(storageType(matrixIdx))
      CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
        jacobianMatrix%storageType=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
        jacobianMatrix%storageType=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE        
      CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
        jacobianMatrix%storageType=DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
        jacobianMatrix%storageType=DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        jacobianMatrix%storageType=DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
        jacobianMatrix%storageType=DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE
      CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
        jacobianMatrix%storageType=DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE
      CASE DEFAULT
        localError="The specified storage type of "//TRIM(NumberToVString(storageType(matrixIdx),"*",err,error))// &
          & " for the Jacobian matrix is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !matrixIdx
    
    EXITS("EquationsMatrices_NonlinearStorageTypeSet0")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NonlinearStorageTypeSet0",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NonlinearStorageTypeSet0

  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of all nonlinear (Jacobian) equations matrices
  SUBROUTINE EquationsMatrices_NonlinearStorageTypeSet1(vectorMatrices,storageType,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(IN) :: storageType !<storageType. The storage type for all Jacobian equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG), ALLOCATABLE :: storageTypes(:)
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices

    ENTERS("EquationsMatrices_NonlinearStorageTypeSet1",err,error,*998)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*998)
    IF(vectorMatrices%vectorMatricesFinished) CALL FlagError("Vector equations matrices have been finished.",err,error,*998)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*998)

    ALLOCATE(storageTypes(nonlinearMatrices%numberOfJacobians),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate storage types.",err,error,*999)
    storageTypes=storageType
    CALL EquationsMatrices_NonlinearStorageTypeSet0(vectorMatrices,storageTypes,err,error,*999)
    DEALLOCATE(storageTypes)

    EXITS("EquationsMatrices_NonlinearStorageTypeSet1")
    RETURN
999 IF(ALLOCATED(storageTypes)) DEALLOCATE(storageTypes)
998 ERRORSEXITS("EquationsMatrices_NonlinearStorageTypeSet1",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NonlinearStorageTypeSet1

  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of the dynamic equations matrices
  SUBROUTINE EquationsMatrices_DynamicStructureTypeSet(vectorMatrices,structureType,err,error,*)
    
    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(IN) :: structureType(:) !<structureType(matrixIdx). The storage type for the matrixIdx'th dynamic equations matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMatrices_DynamicStructureTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)
    IF(vectorMatrices%vectorMatricesFinished) CALL FlagError("Vector equations matrices have already been finished.",err,error,*999)
    NULLIFY(dynamicMatrices)
    CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
    IF(SIZE(structureType,1)/=dynamicMatrices%numberOfDynamicMatrices) THEN
      localError="The size of the structure type array ("//TRIM(NumberToVString(SIZE(structureType,1),"*",err,error))// &
        & ") is not equal to the number of dynamic matrices ("// &
        & TRIM(NumberToVString(dynamicMatrices%numberOfDynamicMatrices,"*",err,error))//")."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
      NULLIFY(equationsMatrix)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,equationsMatrix,err,error,*999)
      SELECT CASE(structureType(matrixIdx))
      CASE(EQUATIONS_MATRIX_NO_STRUCTURE)
        equationsMatrix%structureType=EQUATIONS_MATRIX_NO_STRUCTURE
      CASE(EQUATIONS_MATRIX_FEM_STRUCTURE)
        equationsMatrix%structureType=EQUATIONS_MATRIX_FEM_STRUCTURE
      CASE(EQUATIONS_MATRIX_DIAGONAL_STRUCTURE)
        equationsMatrix%structureType=EQUATIONS_MATRIX_DIAGONAL_STRUCTURE
      CASE(EQUATIONS_MATRIX_NODAL_STRUCTURE)
        equationsMatrix%structureType=EQUATIONS_MATRIX_NODAL_STRUCTURE
      CASE DEFAULT
        localError="The specified strucutre type of "//TRIM(NumberToVString(structureType(matrixIdx),"*",err,error))// &
          & " for dynamic matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !matrixIdx
     
    EXITS("EquationsMatrices_DynamicStructureTypeSet")
    RETURN
999 ERRORSEXITS("EquationsMatrices_DynamicStructureTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_DynamicStructureTypeSet
  
  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of the linear equations matrices
  SUBROUTINE EquationsMatrices_LinearStructureTypeSet(vectorMatrices,structureType,err,error,*)
    
    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(IN) :: structureType(:) !<structureType(matrixIdx). The storage type for the matrixIdx'th linear equations matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMatrices_LinearStructureTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)
    IF(vectorMatrices%vectorMatricesFinished) CALL FlagError("Vector equations matrices have been finished.",err,error,*999)
    NULLIFY(linearMatrices)
    CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
    IF(SIZE(structureType,1)/=linearMatrices%numberOfLinearMatrices) THEN
      localError="The size of the structure type array ("//TRIM(NumberToVString(SIZE(structureType,1),"*",err,error))// &
        & ") is not equal to the number of linear matrices ("// &
        & TRIM(NumberToVString(linearMatrices%numberOfLinearMatrices,"*",err,error))//")."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
      NULLIFY(equationsMatrix)
      CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,equationsMatrix,err,error,*999)
      SELECT CASE(structureType(matrixIdx))
      CASE(EQUATIONS_MATRIX_NO_STRUCTURE)
        equationsMatrix%structureType=EQUATIONS_MATRIX_NO_STRUCTURE
      CASE(EQUATIONS_MATRIX_FEM_STRUCTURE)
        equationsMatrix%structureType=EQUATIONS_MATRIX_FEM_STRUCTURE
      CASE(EQUATIONS_MATRIX_DIAGONAL_STRUCTURE)
        equationsMatrix%structureType=EQUATIONS_MATRIX_DIAGONAL_STRUCTURE
      CASE(EQUATIONS_MATRIX_NODAL_STRUCTURE)
        equationsMatrix%structureType=EQUATIONS_MATRIX_NODAL_STRUCTURE
      CASE DEFAULT
        localError="The specified strucutre type of "//TRIM(NumberToVString(structureType(matrixIdx),"*",err,error))// &
          & " for linear matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !matrixIdx
    
    EXITS("EquationsMatrices_LinearStructureTypeSet")
    RETURN
999 ERRORSEXITS("EquationsMatrices_LinearStructureTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_LinearStructureTypeSet
  
  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of the nonlinear (Jacobian) equations matrices
  SUBROUTINE EquationsMatrices_NonlinearStructureTypeSet0(vectorMatrices,structureType,err,error,*)
    
    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(IN) :: structureType(:) !<structureType(matrixIdx). The structure type for the matrixIdx'th Jacobian equations matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsJacobianType), POINTER :: jacobianMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMatrices_NonlinearStructureTypeSet0",err,error,*999)

   IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)
    IF(vectorMatrices%vectorMatricesFinished) CALL FlagError("Vector equations matrices have been finished.",err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
    IF(SIZE(structureType,1)/=nonlinearMatrices%numberOfJacobians) THEN
      localError="The size of the structure type array ("//TRIM(NumberToVString(SIZE(structureType,1),"*",err,error))// &
        & ") is not equal to the number of Jacobian matrices ("// &
        & TRIM(NumberToVString(nonlinearMatrices%numberOfJacobians,"*",err,error))//")."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    DO matrixIdx=1,nonlinearMatrices%numberOfJacobians
      NULLIFY(jacobianMatrix)
      CALL EquationsMatricesNonlinear_JacobianMatrixGet(nonlinearMatrices,matrixIdx,jacobianMatrix,err,error,*999)
      SELECT CASE(structureType(matrixIdx))
      CASE(EQUATIONS_MATRIX_NO_STRUCTURE)
        jacobianMatrix%structureType=EQUATIONS_MATRIX_NO_STRUCTURE
      CASE(EQUATIONS_MATRIX_FEM_STRUCTURE)
        jacobianMatrix%structureType=EQUATIONS_MATRIX_FEM_STRUCTURE
      CASE(EQUATIONS_MATRIX_DIAGONAL_STRUCTURE)
        jacobianMatrix%structureType=EQUATIONS_MATRIX_DIAGONAL_STRUCTURE
      CASE(EQUATIONS_MATRIX_NODAL_STRUCTURE)
        jacobianMatrix%structureType=EQUATIONS_MATRIX_NODAL_STRUCTURE
      CASE DEFAULT
        localError="The specified strucutre type of "//TRIM(NumberToVString(structureType(matrixIdx),"*",err,error))// &
          & " for the Jacobian matrix is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !matrixIdx
    
    EXITS("EquationsMatrices_NonlinearStructureTypeSet0")
    RETURN
999 ERRORS("EquationsMatrices_NonlinearStructureTypeSet0",err,error)
    EXITS("EquationsMatrices_NonlinearStructureTypeSet0")
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NonlinearStructureTypeSet0

  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of all nonlinear (Jacobian) equations matrices
  SUBROUTINE EquationsMatrices_NonlinearStructureTypeSet1(vectorMatrices,structureType,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(IN) :: structureType !<The structure type for all Jacobian equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG), ALLOCATABLE :: structureTypeS(:)
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices

    ENTERS("EquationsMatrices_NonlinearStructureTypeSet1",err,error,*998)

    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*998)
    IF(vectorMatrices%vectorMatricesFinished) CALL FlagError("Vector equations matrices have been finished.",err,error,*998)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*998)

    ALLOCATE(structureTypes(nonlinearMatrices%numberOfJacobians),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate storage types.",err,error,*999)
    structureTypes=structureType
    CALL EquationsMatrices_NonlinearStructureTypeSet0(vectorMatrices,structureTypes,err,error,*999)
    DEALLOCATE(structureTypes)

    EXITS("EquationsMatrices_NonlinearStructureTypeSet1")
    RETURN
999 IF(ALLOCATED(structureTypes)) DEALLOCATE(structureTypes)
998 ERRORS("EquationsMatrices_NonlinearStructureTypeSet1",err,error)
    EXITS("EquationsMatrices_NonlinearStructureTypeSet1")
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NonlinearStructureTypeSet1

  !
  !================================================================================================================================
  !

  !>Finalise the scalar equations matrices and deallocate all memory.
  SUBROUTINE EquationsMatrices_ScalarFinalise(scalarMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesScalarType), POINTER :: scalarMatrices !<A pointer to the scalar equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("EquationsMatrices_ScalarFinalise",err,error,*999)

    IF(ASSOCIATED(scalarMatrices)) THEN
      DEALLOCATE(scalarMatrices)
    ENDIF
       
    EXITS("EquationsMatrices_ScalarFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_ScalarFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_ScalarFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the scalar equations matrices for the scalar equations.
  SUBROUTINE EquationsMatrices_ScalarInitialise(scalarEquations,err,error,*)
    
    !Argument variables
    TYPE(EquationsScalarType), POINTER :: scalarEquations !<A pointer to the scalar equations to initialise the scalar equations matrices for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(EquationsMappingScalarType), POINTER :: scalarMapping
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("EquationsMatrices_ScalarInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(scalarEquations)) CALL FlagError("Scalar equations is not associated.",err,error,*998)
    IF(ASSOCIATED(scalarEquations%scalarMatrices)) &
      & CALL FlagError("Scalar equations matrices is already associated for this equations.",err,error,*998)
    NULLIFY(scalarMapping)
    CALL EquationsScalar_ScalarMappingGet(scalarEquations,scalarMapping,err,error,*998)
    IF(.NOT.scalarMapping%scalarMappingFinished) CALL FlagError("Scalar equations mapping has not been finished.",err,error,*998)
    
    ALLOCATE(scalarEquations%scalarMatrices,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate scalar equations scalar equations matrices.",err,error,*999)
    scalarEquations%scalarMatrices%scalarEquations=>scalarEquations
    scalarEquations%scalarMatrices%scalarMatricesFinished=.FALSE.
    scalarEquations%scalarMatrices%scalarMapping=>scalarMapping
    NULLIFY(scalarEquations%scalarMatrices%solverMapping)
    NULLIFY(scalarEquations%scalarMatrices%functions)
    NULLIFY(scalarEquations%scalarMatrices%normMatrices)
    NULLIFY(scalarEquations%scalarMatrices%dotProductMatrices)
    NULLIFY(scalarEquations%scalarMatrices%quadraticMatrices)
       
    EXITS("EquationsMatrices_ScalarInitialise")
    RETURN
999 CALL EquationsMatrices_ScalarFinalise(scalarEquations%scalarMatrices,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_ScalarInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_ScalarInitialise

  !
  !================================================================================================================================
  !

  !>Finalise the vector equations matrices and deallocate all memory.
  SUBROUTINE EquationsMatrices_VectorFinalise(vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the vector equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("EquationsMatrices_VectorFinalise",err,error,*999)

    IF(ASSOCIATED(vectorMatrices)) THEN
      CALL EquationsMatrices_DynamicFinalise(vectorMatrices%dynamicMatrices,err,error,*999)
      CALL EquationsMatrices_LinearFinalise(vectorMatrices%linearMatrices,err,error,*999)
      CALL EquationsMatrices_NonlinearFinalise(vectorMatrices%nonlinearMatrices,err,error,*999)
      CALL EquationsMatrices_RHSFinalise(vectorMatrices%rhsVector,err,error,*999)      
      CALL EquationsMatrices_SourceFinalise(vectorMatrices%sourceVector,err,error,*999)      
      DEALLOCATE(vectorMatrices)
    ENDIF
       
    EXITS("EquationsMatrices_VectorFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_VectorFinalise",err,error)
    RETURN 1
  END SUBROUTINE EquationsMatrices_VectorFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the vector equations matrices for the vector equations.
  SUBROUTINE EquationsMatrices_VectorInitialise(vectorEquations,err,error,*)
    
     !Argument variables
    TYPE(EquationsVectorType), POINTER :: vectorEquations !<A pointer to the vector equations to initialise the vector equations matrices for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("EquationsMatrices_VectorInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is not associated.",err,error,*998)
    IF(ASSOCIATED(vectorEquations%vectorMatrices)) &
      & CALL FlagError("Vector equations matrices is already associated for this equations.",err,error,*998)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*998)
    IF(.NOT.vectorMapping%vectorMappingFinished) CALL FlagError("Vector equations mapping has not been finished.",err,error,*998)
    
    ALLOCATE(vectorEquations%vectorMatrices,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate vector equations vector equations matrices.",err,error,*999)
    vectorEquations%vectorMatrices%vectorEquations=>vectorEquations
    vectorEquations%vectorMatrices%vectorMatricesFinished=.FALSE.
    vectorEquations%vectorMatrices%vectorMapping=>vectorMapping
    NULLIFY(vectorEquations%vectorMatrices%solverMapping)
    vectorEquations%vectorMatrices%numberOfRows=vectorMapping%numberOfRows
    vectorEquations%vectorMatrices%totalNumberOfRows=vectorMapping%totalNumberOfRows
    vectorEquations%vectorMatrices%numberOfGlobalRows=vectorMapping%numberOfGlobalRows
    NULLIFY(vectorEquations%vectorMatrices%dynamicMatrices)
    NULLIFY(vectorEquations%vectorMatrices%linearMatrices)
    NULLIFY(vectorEquations%vectorMatrices%nonlinearMatrices)
    NULLIFY(vectorEquations%vectorMatrices%rhsVector)
    NULLIFY(vectorEquations%vectorMatrices%sourceVector)            
    CALL EquationsMatrices_DynamicInitialise(vectorEquations%vectorMatrices,err,error,*999)            
    CALL EquationsMatrices_LinearInitialise(vectorEquations%vectorMatrices,err,error,*999)            
    CALL EquationsMatrices_NonlinearInitialise(vectorEquations%vectorMatrices,err,error,*999)            
    CALL EquationsMatrices_RHSInitialise(vectorEquations%vectorMatrices,err,error,*999)            
    CALL EquationsMatrices_SourceInitialise(vectorEquations%vectorMatrices,err,error,*999)            
       
    EXITS("EquationsMatrices_VectorInitialise")
    RETURN
999 CALL EquationsMatrices_VectorFinalise(vectorEquations%vectorMatrices,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_VectorInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_VectorInitialise

  !
  !================================================================================================================================
  !

  !>Initialise the values of the equations matrices and vectors to the given value e.g., 0.0_DP
  SUBROUTINE EquationsMatrices_VectorValuesInitialise(vectorMatrices,selectionType,value,err,error,*)
    
    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations matrices to initialise the values for
    INTEGER(INTG), INTENT(IN) :: selectionType !<The selection of equations matrices to be initialised \see EquationsMatricesRoutines::SelectMatricesTypes,EquationsMatricesRoutines
    REAL(DP), INTENT(IN) :: value !<The value to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsJacobianType), POINTER :: jacobianMatrix
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    
    ENTERS("EquationsMatrices_VectorValuesInitialise",err,error,*999)
    
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector equations matrices is not associated.",err,error,*999)
    
    IF(selectionType==EQUATIONS_MATRICES_ALL.OR. &
      & selectionType==EQUATIONS_MATRICES_DYNAMIC_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_LINEAR_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_NONLINEAR_ONLY) THEN
      dynamicMatrices=>vectorMatrices%dynamicMatrices
      IF(ASSOCIATED(dynamicMatrices)) THEN
        DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
          NULLIFY(equationsMatrix)
          CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,equationsMatrix,err,error,*999)
          IF(equationsMatrix%updateMatrix) CALL DistributedMatrix_AllValuesSet(equationsMatrix%matrix,value,err,error,*999)
        ENDDO !matrixIdx
      ENDIF
    ENDIF
    IF(selectionType==EQUATIONS_MATRICES_ALL.OR. &
      & selectionType==EQUATIONS_MATRICES_DYNAMIC_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_LINEAR_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_NONLINEAR_ONLY) THEN
      linearMatrices=>vectorMatrices%linearMatrices
      IF(ASSOCIATED(linearMatrices)) THEN
        DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
          NULLIFY(equationsMatrix)
          CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,equationsMatrix,err,error,*999)
          IF(equationsMatrix%updateMatrix) CALL DistributedMatrix_AllValuesSet(equationsMatrix%matrix,value,err,error,*999)
       ENDDO !matrixIdx
      ENDIF
    ENDIF
    IF(selectionType==EQUATIONS_MATRICES_ALL.OR. &
      & selectionType==EQUATIONS_MATRICES_NONLINEAR_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_JACOBIAN_ONLY) THEN
      nonlinearMatrices=>vectorMatrices%nonlinearMatrices
      IF(ASSOCIATED(nonlinearMatrices)) THEN
        DO matrixIdx=1,nonlinearMatrices%numberOfJacobians
          NULLIFY(jacobianMatrix)
          CALL EquationsMatricesNonlinear_JacobianMatrixGet(nonlinearMatrices,matrixIdx,jacobianMatrix,err,error,*999)
          IF(jacobianMatrix%updateJacobian) CALL DistributedMatrix_AllValuesSet(jacobianMatrix%jacobian,VALUE,err,error,*999)
        ENDDO !matrixIdx
      ENDIF
    ENDIF
    IF(selectionType==EQUATIONS_MATRICES_ALL.OR. &
      & selectionType==EQUATIONS_MATRICES_NONLINEAR_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_RESIDUAL_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_RHS_RESIDUAL_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_RESIDUAL_SOURCE_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_VECTORS_ONLY) THEN
      nonlinearMatrices=>vectorMatrices%nonlinearMatrices
      IF(ASSOCIATED(nonlinearMatrices)) THEN
        IF(nonlinearMatrices%updateResidual) CALL DistributedVector_AllValuesSet(nonlinearMatrices%residual,value,err,error,*999)
      ENDIF
    ENDIF
    IF(selectionType==EQUATIONS_MATRICES_ALL.OR. &
      & selectionType==EQUATIONS_MATRICES_DYNAMIC_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_LINEAR_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_NONLINEAR_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_RHS_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_RHS_RESIDUAL_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_RHS_SOURCE_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_VECTORS_ONLY) THEN    
      rhsVector=>vectorMatrices%rhsVector
      IF(ASSOCIATED(rhsVector)) THEN
        IF(rhsVector%updateVector) CALL DistributedVector_AllValuesSet(rhsVector%vector,value,err,error,*999)
      ENDIF
    ENDIF
    IF(selectionType==EQUATIONS_MATRICES_ALL.OR. &
      & selectionType==EQUATIONS_MATRICES_DYNAMIC_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_LINEAR_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_NONLINEAR_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_SOURCE_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_RHS_SOURCE_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_RESIDUAL_SOURCE_ONLY.OR. &
      & selectionType==EQUATIONS_MATRICES_VECTORS_ONLY) THEN    
      sourceVector=>vectorMatrices%sourceVector
      IF(ASSOCIATED(sourceVector)) THEN
        IF(sourceVector%updateVector) CALL DistributedVector_AllValuesSet(sourceVector%vector,value,err,error,*999)
      ENDIF
    ENDIF
 
    EXITS("EquationsMatrices_VectorValuesInitialise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_VectorValuesInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_VectorValuesInitialise

  !
  !================================================================================================================================
  !

  !>Caclulates the matrix structure (sparsity) for an equations matrix.
  SUBROUTINE EquationsMatrix_StructureCalculate(equationsMatrix,numberOfNonZeros,rowIndices,columnIndices,list,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to calculate the structure for
    INTEGER(INTG), INTENT(OUT) :: numberOfNonZeros !<On return the number of non-zeros in the matrix
    INTEGER(INTG), POINTER :: rowIndices(:) !<On return a pointer to row location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: columnIndices(:) !<On return a pointer to the column location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    TYPE(LinkedList), POINTER :: list(:) 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,component,componentIdx,derivative,derivativeIdx,dofIdx,dummyErr,element,elementIdx, &
      & globalColumn,localColumn,localDOF,localDOFIdx,localNodeIdx,matrixNumber,node,node2,numberOfColumns, &
      & numberOfDerivatives,numberOfVersions,version,versionIdx
    INTEGER(INTG), ALLOCATABLE :: columns(:)
    REAL(DP) :: sparsity
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: dependentDofsDomainMapping
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: domainElements
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_DOF_TO_PARAM_MAP_TYPE), POINTER :: dependentDofsParamMapping
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: columnIndicesLists(:)
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("EquationsMatrix_StructureCalculate",err,error,*998)

    numberOfNonZeros=0
    
    IF(.NOT.ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is not associated.",err,error,*998)
    IF(ASSOCIATED(rowIndices)) CALL FlagError("Row indices is already associated.",err,error,*998)
    IF(ASSOCIATED(columnIndices)) CALL FlagError("Column indices is already associated.",err,error,*998)
    linearMatrices=>equationsMatrix%linearMatrices
    dynamicMatrices=>equationsMatrix%dynamicMatrices
    IF(.NOT.ASSOCIATED(dynamicMatrices).AND..NOT.ASSOCIATED(linearMatrices)) &
      & CALL FlagError("Either equations matrix dynamic or linear matrices is not associated.",err,error,*998)
    NULLIFY(vectorMatrices)
    IF(ASSOCIATED(dynamicMatrices)) THEN
      CALL EquationsMatricesDynamic_VectorMatricesGet(dynamicMatrices,vectorMatrices,err,error,*999)
    ELSE
      CALL EquationsMatricesLinear_VectorMatricesGet(linearMatrices,vectorMatrices,err,error,*999)
    ENDIF
    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*998)
    NULLIFY(vectorEquations)
    CALL EquationsMatricesVector_VectorEquationsGet(vectorMatrices,vectorEquations,err,error,*998)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*998)
    NULLIFY(equationsSet)
    CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*998)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*998)
    matrixNumber=equationsMatrix%matrixNumber
    NULLIFY(dynamicMapping)
    NULLIFY(linearMapping)
    IF(ASSOCIATED(dynamicMatrices)) THEN
      CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      fieldVariable=>dynamicMapping%equationsMatrixToVarMaps(matrixNumber)%variable
    ELSE
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      fieldVariable=>linearMapping%equationsMatrixToVarMaps(matrixNumber)%variable
    ENDIF
    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Dependent field variable is not associated.",err,error,*998)          
    dependentDofsDomainMapping=>fieldVariable%DOMAIN_MAPPING
    IF(.NOT.ASSOCIATED(dependentDofsDomainMapping)) &
      & CALL FlagError("Dependent dofs domain mapping is not associated.",err,error,*998)
    dependentDofsParamMapping=>fieldVariable%DOF_TO_PARAM_MAP
    IF(.NOT.ASSOCIATED(dependentDofsParamMapping)) &
      & CALL FlagError("Dependent dofs parameter mapping is not associated.",err,error,*998)
    
    SELECT CASE(equationsMatrix%structureType)
    CASE(EQUATIONS_MATRIX_NO_STRUCTURE)
      CALL FlagError("There is no structure to calculate for a matrix with no structure.",err,error,*998)
    CASE(EQUATIONS_MATRIX_FEM_STRUCTURE)
      SELECT CASE(equationsMatrix%storageType)
      CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        
        !Allocate lists
        ALLOCATE(columnIndicesLists(dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate column indices lists.",err,error,*999)
        !Allocate row indices
        ALLOCATE(rowIndices(dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL+1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate row indices.",err,error,*999)
        rowIndices(1)=1
        
        !First, loop over the rows and calculate the number of non-zeros
        numberOfNonZeros=0
        DO localDOFIdx=1,dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
          IF(dependentDofsParamMapping%DOF_TYPE(1,localDOFIdx)/=FIELD_NODE_DOF_TYPE) THEN
            localError="Local DOF number "//TRIM(NumberToVString(localDOFIdx,"*",err,error))//" is not a node based DOF."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          dofIdx=dependentDofsParamMapping%DOF_TYPE(2,localDOFIdx) !value for a particular field dof (localDOFIdx)
          node=dependentDofsParamMapping%NODE_DOF2PARAM_MAP(3,dofIdx) !node number of the field parameter
          component=dependentDofsParamMapping%NODE_DOF2PARAM_MAP(4,dofIdx) !component number of the field parameter
          domainNodes=>fieldVariable%components(component)%domain%topology%nodes          
          !Set up list
          NULLIFY(columnIndicesLists(localDOFIdx)%ptr)
          CALL List_CreateStart(columnIndicesLists(localDOFIdx)%ptr,err,error,*999)
          CALL List_DataTypeSet(columnIndicesLists(localDOFIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
          CALL List_InitialSizeSet(columnIndicesLists(localDOFIdx)%ptr,domainNodes%nodes(node)% &
            & NUMBER_OF_SURROUNDING_ELEMENTS*fieldVariable%components(component)% &
            & maxNumberElementInterpolationParameters,err,error,*999)
          CALL List_CreateFinish(columnIndicesLists(localDOFIdx)%ptr,err,error,*999)
          !Loop over all elements containing the dof
          DO elementIdx=1,domainNodes%nodes(node)%NUMBER_OF_SURROUNDING_ELEMENTS
            element=domainNodes%nodes(node)%SURROUNDING_ELEMENTS(elementIdx)
            DO componentIdx=1,fieldVariable%NUMBER_OF_COMPONENTS
              domainElements=>fieldVariable%components(componentIdx)%domain%topology%elements
              basis=>domainElements%elements(element)%basis
              DO localNodeIdx=1,basis%NUMBER_OF_NODES
                node2=domainElements%elements(element)%ELEMENT_NODES(localNodeIdx)
                DO derivativeIdx=1,basis%NUMBER_OF_DERIVATIVES(localNodeIdx)
                  derivative=domainElements%elements(element)%ELEMENT_DERIVATIVES(derivativeIdx,localNodeIdx)
                  version=domainElements%elements(element)%elementVersions(derivativeIdx,localNodeIdx)
                  !Find the local and global column and add the global column to the indices list
                  localColumn=fieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                    & nodes(node2)%derivatives(derivative)%versions(version)
                  globalColumn=fieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localColumn)
                  
                  CALL List_ItemAdd(columnIndicesLists(localDOFIdx)%ptr,globalColumn,err,error,*999)
                  
                ENDDO !derivative
              ENDDO !localNodeIdx
            ENDDO !componentIdx
          ENDDO !elementIdx
          CALL List_RemoveDuplicates(columnIndicesLists(localDOFIdx)%ptr,err,error,*999)
          CALL List_NumberOfItemsGet(columnIndicesLists(localDOFIdx)%ptr,numberOfColumns,err,error,*999)
          numberOfNonZeros=numberOfNonZeros+numberOfColumns
          rowIndices(localDOFIdx+1)=numberOfNonZeros+1
        ENDDO !localDOFIdx
        
      CASE DEFAULT
        localError="The matrix storage type of "//TRIM(NumberToVString(equationsMatrix%storageType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      
    CASE(EQUATIONS_MATRIX_NODAL_STRUCTURE)
      SELECT CASE(equationsMatrix%storageType)
      CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        !Allocate lists
        ALLOCATE(columnIndicesLists(dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate column indices lists.",err,error,*999)
        !Allocate row indices
        ALLOCATE(rowIndices(dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL+1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate row indices.",err,error,*999)
        rowIndices(1)=1
        
        !First, loop over the rows and calculate the number of non-zeros
        numberOfNonZeros=0
        DO localDofIdx=1,dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
          IF(dependentDofsParamMapping%DOF_TYPE(1,localDOFIdx)/=FIELD_NODE_DOF_TYPE) THEN
            localError="Local DOF number "//TRIM(NumberToVString(localDOFIdx,"*",err,error))//" is not a node based DOF."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          dofIdx=dependentDofsParamMapping%DOF_TYPE(2,localDofIdx)!value for a particular field dof (localDofIdx)
          node=dependentDofsParamMapping%NODE_DOF2PARAM_MAP(3,dofIdx)!node number (node) of the field parameter
          component=dependentDofsParamMapping%NODE_DOF2PARAM_MAP(4,dofIdx)!component number (component) of the field parameter
          domainNodes=>fieldVariable%components(component)%domain%topology%NODES
          
          !Set up list
          NULLIFY(columnIndicesLists(localDofIdx)%ptr)
          CALL List_CreateStart(columnIndicesLists(localDofIdx)%ptr,err,error,*999)
          CALL List_DataTypeSet(columnIndicesLists(localDofIdx)%ptr,LIST_INTG_TYPE,err,error,*999)          
          CALL List_InitialSizeSet(columnIndicesLists(localDofIdx)%ptr,fieldVariable%NUMBER_OF_COMPONENTS* &
            & fieldVariable%maxNumberElementInterpolationParameters,err,error,*999)          
          CALL List_CreateFinish(columnIndicesLists(localDofIdx)%ptr,err,error,*999)
          !Loop over all components, nodes, derivatives and versions
          DO componentIdx=1,fieldVariable%NUMBER_OF_COMPONENTS
            numberOfDerivatives=fieldVariable%components(componentIdx)%domain%topology%nodes%nodes(node)%NUMBER_OF_DERIVATIVES
            DO derivativeIdx=1,numberOfDerivatives
              numberOfVersions=fieldVariable%components(componentIdx)%domain%topology%nodes%nodes(node)% &
                & derivatives(derivativeIdx)%numberOfVersions
              DO versionIdx=1,numberOfVersions
                localColumn=fieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                  & nodes(node)%derivatives(derivativeIdx)%versions(versionIdx)
                globalColumn=fieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localColumn)
                
                CALL List_ItemAdd(columnIndicesLists(localDofIdx)%ptr,globalColumn,err,error,*999)
                
              ENDDO !versionIdx
            ENDDO !derivativeIdx
          ENDDO !componentIdx            
          CALL List_RemoveDuplicates(columnIndicesLists(localDofIdx)%ptr,err,error,*999)
          CALL List_NumberOfItemsGet(columnIndicesLists(localDofIdx)%ptr,numberOfColumns,err,error,*999)
          numberOfNonZeros=numberOfNonZeros+numberOfColumns
          rowIndices(localDofIdx+1)=numberOfNonZeros+1
        ENDDO !localDofIdx
        
      CASE DEFAULT
        localError="The matrix storage type of "//TRIM(NumberToVString(equationsMatrix%storageType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      
    CASE(EQUATIONS_MATRIX_DIAGONAL_STRUCTURE)
      CALL FlagError("There is not structure to calculate for a diagonal matrix.",err,error,*998)
    CASE DEFAULT
      localError="The matrix structure type of "//TRIM(NumberToVString(equationsMatrix%structureType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*998)
    END SELECT

    IF(numberOfNonZeros>0) THEN
      !Allocate and setup the column locations
      ALLOCATE(columnIndices(numberOfNonZeros),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate column indices.",err,error,*999)        
      ALLOCATE(list(dependentDofsDomainMapping%NUMBER_OF_GLOBAL),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate list.",err,error,*999)
      
      DO localDOFIdx=1,dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL          
        CALL List_DetachAndDestroy(columnIndicesLists(localDOFIdx)%ptr,numberOfColumns,columns,err,error,*999)        
        DO columnIdx=1,numberOfColumns
          !columns stores the list of nonzero column indices for each local row (localDOFIdx)
          columnIndices(rowIndices(localDOFIdx)+columnIdx-1)=columns(columnIdx)             
          !global to local columns
          IF(ASSOCIATED(linearMapping).OR.ASSOCIATED(dynamicMapping)) THEN
            IF(ASSOCIATED(dynamicMatrices)) THEN
              localColumn=vectorMatrices%vectorMapping%dynamicMapping &
                & %equationsMatrixToVarMaps(1)%columnDOFSMapping%global_to_local_map &
                & (columns(columnIdx))%LOCAL_NUMBER(1)
              localDOF = localColumn
              !Column to dof mapping?
              !localDOF=vectorMatrices%vectorMapping%dynamicMapping%equationsMatrixToVarMaps(1)%columnToDOFMap(localColumn)
            ELSE
              localColumn=vectorMatrices%vectorMapping%linearMapping%equationsMatrixToVarMaps(1)%columnDOFSMapping% &
                & global_To_Local_Map(columns(columnIdx))%LOCAL_NUMBER(1)
              localDOF = localColumn
            ENDIF
          ENDIF
          dofIdx=dependentDofsParamMapping%DOF_TYPE(2,localDOF)
          node=dependentDofsParamMapping%NODE_DOF2PARAM_MAP(3,dofIdx)
          component=dependentDofsParamMapping%NODE_DOF2PARAM_MAP(4,dofIdx)
          domainNodes=>fieldVariable%components(component)%domain%topology%NODES
          
          !Check whether boundary node    
          IF(domainNodes%nodes(node)%BOUNDARY_NODE) CALL LinkedList_Add(list(columns(columnIdx)),localDOFIdx,err,error,*999)
          
        ENDDO !columnIdx
        DEALLOCATE(columns)                                    
      ENDDO !localDOFIdx
    ENDIF
      
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Equations matrix structure:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Equations matrix number : ",matrixNumber,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ",dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL, &
        & err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ",dependentDofsDomainMapping%NUMBER_OF_GLOBAL, &
        & err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of non zeros = ",numberOfNonZeros,err,error,*999)
      IF(dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL*dependentDofsDomainMapping%NUMBER_OF_GLOBAL/=0) THEN
        sparsity=(1.0_DP-REAL(numberOfNonZeros,DP)/REAL(dependentDofsDomainMapping% &
          & TOTAL_NUMBER_OF_LOCAL*dependentDofsDomainMapping%NUMBER_OF_GLOBAL,DP))*100.0_DP
        CALL WriteStringFmtValue(DIAGNOSTIC_OUTPUT_TYPE,"  Sparsity (% of zeros) = ",sparsity,"F6.2",err,error,*999)
      ENDIF
      IF(numberOfNonZeros>0) THEN
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL+1,8,8, &
          & rowIndices,'("  Row indices    :",8(X,I13))','(18X,8(X,I13))',err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfNonZeros,8,8,columnIndices, &
          & '("  Column indices :",8(X,I13))','(18X,8(X,I13))', err,error,*999)
      ENDIF
    ENDIF
      
    EXITS("EquationsMatrix_StructureCalculate")
    RETURN
999 IF(ASSOCIATED(rowIndices)) DEALLOCATE(rowIndices)
    IF(ASSOCIATED(columnIndices)) DEALLOCATE(columnIndices)
    IF(ALLOCATED(columns)) DEALLOCATE(columns)
    IF(ALLOCATED(columnIndicesLists)) THEN
      DO localDofIdx=1,SIZE(columnIndicesLists,1)
        IF(ASSOCIATED(columnIndicesLists(localDofIdx)%ptr)) &
          & CALL List_Destroy(columnIndicesLists(localDofIdx)%ptr,dummyErr,dummyError,*998)
      ENDDO !localDofIdx
      DEALLOCATE(columnIndicesLists)
    ENDIF
998 ERRORSEXITS("EquationsMatrix_StructureCalculate",err,error)
    RETURN 1
  END SUBROUTINE EquationsMatrix_StructureCalculate

  !
  !================================================================================================================================
  !

  !>Caclulates the matrix structure (sparsity) for a Jacobian matrix.
  SUBROUTINE JacobianMatrix_StructureCalculate(jacobianMatrix,numberOfNonZeros,rowIndices,columnIndices,err,error,*)

    !Argument variables
    TYPE(EquationsJacobianType), POINTER :: jacobianMatrix !<A pointer to the Jacobian matrix to calculate the strucute for
    INTEGER(INTG), INTENT(OUT) :: numberOfNonZeros !<On return the number of non-zeros in the matrix
    INTEGER(INTG), POINTER :: rowIndices(:) !<On return a pointer to row location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: columnIndices(:) !<On return a pointer to the column location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,component,componentIdx,derivative,derivativeIdx,dofIdx,dummyErr,element,elementIdx,globalColumn, &
      & localColumn,localDOFIdx,localNodeIdx,matrixNumber,numberOfColumns,node,node2,numberOfDerivatives, &
      & numberOfVersions,version,versionIdx
    INTEGER(INTG), ALLOCATABLE :: columns(:)
    REAL(DP) :: sparsity
    TYPE(BASIS_TYPE), POINTER :: basis,basis2
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: dependentDofsDomainMapping,rowDofsDomainMapping
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: domainElements,domainElements2
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: NONlinearMatrices
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_DOF_TO_PARAM_MAP_TYPE), POINTER :: dependentDofsParamMapping,rowDofsParamMapping
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable,rowVariable
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: columnIndicesLists(:)
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("JacobianMatrix_StructureCalculate",err,error,*998)

    numberOfNonZeros=0
    IF(.NOT.ASSOCIATED(jacobianMatrix)) CALL FlagError("Jacobian matrix is not associated.",err,error,*998)
    IF(ASSOCIATED(rowIndices)) CALL FlagError("Row indices is already associated.",err,error,*998)
    IF(ASSOCIATED(columnIndices)) CALL FlagError("Column indices is already associated.",err,error,*998)
    matrixNumber=jacobianMatrix%jacobianNumber
    nonlinearMatrices=>jacobianMatrix%nonlinearMatrices
    IF(.NOT.ASSOCIATED(nonlinearMatrices)) CALL FlagError("Jacobian matrix nonlinear matrices is not associated.",err,error,*998)
    NULLIFY(vectorMatrices)
    CALL EquationsMatricesNonlinear_VectorMatricesGet(nonlinearMatrices,vectorMatrices,err,error,*998)
    NULLIFY(vectorMapping)
    CALL EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*998)
    NULLIFY(vectorEquations)
    CALL EquationsMatricesVector_VectorEquationsGet(vectorMatrices,vectorEquations,err,error,*998)
    NULLIFY(equations)
    CALL EquationsVector_EquationsGet(vectorEquations,equations,err,error,*998)
    NULLIFY(equationsSet)
    CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*998)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*998)
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*998)
    fieldVariable=>nonlinearMapping%jacobianToVarMap(matrixNumber)%VARIABLE
    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Dependent field variable is not associated.",err,error,*998)          
    dependentDofsDomainMapping=>fieldVariable%DOMAIN_MAPPING
    IF(.NOT.ASSOCIATED(dependentDofsDomainMapping)) &
      & CALL FlagError("Dependent dofs domain mapping is not associated.",err,error,*998)
    dependentDofsParamMapping=>fieldVariable%DOF_TO_PARAM_MAP
    IF(.NOT.ASSOCIATED(dependentDofsParamMapping)) &
      & CALL FlagError("Dependent dofs parameter mapping is not associated.",err,error,*998)
    !If RHS variable exists, use this for row DOFs, else use the first nonlinear variable
    IF(ASSOCIATED(vectorMapping%rhsMapping)) THEN
      rowVariable=>vectorMapping%rhsMapping%rhsVariable
    ELSE
      rowVariable=>nonlinearMapping%jacobianToVarMap(1)%VARIABLE
    ENDIF
    IF(.NOT.ASSOCIATED(rowVariable)) CALL FlagError("RHS or first nonlinear variable is not associated",err,error,*999)
    rowDofsDomainMapping=>rowVariable%DOMAIN_MAPPING
    rowDofsParamMapping=>rowVariable%DOF_TO_PARAM_MAP      
    IF(.NOT.ASSOCIATED(rowDofsDomainMapping)) CALL FlagError("Row dofs domain mapping is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(rowDofsParamMapping)) CALL FlagError("Row dofs parameter mapping is not associated.",err,error,*999)
    
    SELECT CASE(jacobianMatrix%structureType)
    CASE(EQUATIONS_MATRIX_NO_STRUCTURE)
      CALL FlagError("Not implemented.",err,error,*998)
    CASE(EQUATIONS_MATRIX_FEM_STRUCTURE)
      SELECT CASE(jacobianMatrix%storageType)
      CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        !Allocate lists
        ALLOCATE(columnIndicesLists(rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate column indices lists.",err,error,*999)
        !Allocate row indices
        ALLOCATE(rowIndices(rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL+1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate row indices.",err,error,*999)
        rowIndices(1)=1
        !First, loop over the rows and calculate the number of non-zeros
        numberOfNonZeros=0
        DO localDOFIdx=1,rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
          SELECT CASE(rowDofsParamMapping%DOF_TYPE(1,localDOFIdx))
          CASE(FIELD_CONSTANT_INTERPOLATION)
            CALL FlagError("Constant interpolation is not implemented yet.",err,error,*999)
          CASE(FIELD_NODE_DOF_TYPE)
            dofIdx=rowDofsParamMapping%DOF_TYPE(2,localDOFIdx)
            node=rowDofsParamMapping%NODE_DOF2PARAM_MAP(3,dofIdx) !node number
            component=rowDofsParamMapping%NODE_DOF2PARAM_MAP(4,dofIdx) !component number
            domainNodes=>rowVariable%components(component)%domain%topology%nodes
            !Set up list
            NULLIFY(columnIndicesLists(localDOFIdx)%ptr)
            CALL List_CreateStart(columnIndicesLists(localDOFIdx)%ptr,err,error,*999)
            CALL List_DataTypeSet(columnIndicesLists(localDOFIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
            CALL List_InitialSizeSet(columnIndicesLists(localDOFIdx)%ptr,domainNodes%nodes(node)% &
              & NUMBER_OF_SURROUNDING_ELEMENTS*rowVariable%components(component)% &
              & maxNumberElementInterpolationParameters,err,error,*999)
            CALL List_CreateFinish(columnIndicesLists(localDOFIdx)%ptr,err,error,*999)
            !Loop over all elements containing the dof
            DO elementIdx=1,domainNodes%nodes(node)%NUMBER_OF_SURROUNDING_ELEMENTS
              element=domainNodes%nodes(node)%SURROUNDING_ELEMENTS(elementIdx)
              DO componentIdx=1,fieldVariable%NUMBER_OF_COMPONENTS
                SELECT CASE(fieldVariable%components(componentIdx)%INTERPOLATION_TYPE)
                CASE(FIELD_CONSTANT_INTERPOLATION)
                  ! do nothing? this will probably never be encountered...?
                  CALL FlagError("Not implemented?",err,error,*999)
                CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                  localColumn=fieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%elements(element)
                  globalColumn=fieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localColumn)
                  CALL List_ItemAdd(columnIndicesLists(localDOFIdx)%ptr,globalColumn,err,error,*999)
                CASE(FIELD_NODE_BASED_INTERPOLATION)
                  domainElements=>fieldVariable%components(componentIdx)%domain%topology%elements
                  basis=>domainElements%elements(element)%basis
                  DO localNodeIdx=1,basis%NUMBER_OF_NODES
                    node2=domainElements%elements(element)%ELEMENT_NODES(localNodeIdx)
                    DO derivativeIdx=1,basis%NUMBER_OF_DERIVATIVES(localNodeIdx)
                      derivative=domainElements%elements(element)%ELEMENT_DERIVATIVES(derivativeIdx,localNodeIdx)
                      version=domainElements%elements(element)%elementVersions(derivativeIdx,localNodeIdx)
                      !Find the local and global column and add the global column to the indices list
                      localColumn=fieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP% &
                        & NODE_PARAM2DOF_MAP%NODES(node2)%derivatives(derivative)%versions(version)
                      globalColumn=fieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localColumn)
                      CALL List_ItemAdd(columnIndicesLists(localDOFIdx)%ptr,globalColumn,err,error,*999)
                    ENDDO !derivative
                  ENDDO !localNodeIdx
                CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                  CALL FlagError("Grid point based interpolation is not implemented yet.",err,error,*999)
                CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                  CALL FlagError("Gauss point based interpolation is not implemented yet.",err,error,*999)
                CASE DEFAULT
                  localError="The interpolation type of "// &
                    & TRIM(NumberToVString(fieldVariable%components(componentIdx)%INTERPOLATION_TYPE,"*",err,error))// &
                    & " for local DOF number "//TRIM(NumberToVString(localDOFIdx,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ENDDO !componentIdx
            ENDDO !elementIdx
            CALL List_RemoveDuplicates(columnIndicesLists(localDOFIdx)%ptr,err,error,*999)
            CALL List_NumberOfItemsGet(columnIndicesLists(localDOFIdx)%ptr,numberOfColumns, &
              & err,error,*999)
            numberOfNonZeros=numberOfNonZeros+numberOfColumns
            rowIndices(localDOFIdx+1)=numberOfNonZeros+1
          CASE(FIELD_ELEMENT_DOF_TYPE)
            ! row corresponds to a variable that's element-wisely interpolated
            dofIdx=rowDofsParamMapping%DOF_TYPE(2,localDOFIdx)          ! dofIdx = index in ELEMENT_DOF2PARAM_MAP
            element=rowDofsParamMapping%ELEMENT_DOF2PARAM_MAP(1,dofIdx)   ! current element (i.e. corresponds to current dof)
            component=rowDofsParamMapping%ELEMENT_DOF2PARAM_MAP(2,dofIdx)   ! current variable component
            domainElements=>rowVariable%components(component)%domain%topology%ELEMENTS
            basis=>domainElements%elements(element)%BASIS
            !Set up list
            NULLIFY(columnIndicesLists(localDOFIdx)%ptr)
            CALL List_CreateStart(columnIndicesLists(localDOFIdx)%ptr,err,error,*999)
            CALL List_DataTypeSet(columnIndicesLists(localDOFIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
            CALL List_InitialSizeSet(columnIndicesLists(localDOFIdx)%ptr, &
              & rowVariable%components(component)%maxNumberElementInterpolationParameters+1,err,error,*999) ! size = all nodal dofs + itself
            CALL List_CreateFinish(columnIndicesLists(localDOFIdx)%ptr,err,error,*999)
            DO componentIdx=1,fieldVariable%NUMBER_OF_COMPONENTS
              domainElements2=>fieldVariable%components(componentIdx)%domain%topology%elements
              basis2=>domainElements2%elements(element)%basis
              SELECT CASE(fieldVariable%components(componentIdx)%INTERPOLATION_TYPE)
              CASE(FIELD_CONSTANT_INTERPOLATION)
                CALL FlagError("Constant interpolation is not implemented yet.",err,error,*999)
              CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                ! it's assumed that element-based variables arne't directly coupled
                ! put a diagonal entry
                localColumn=fieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%elements(element)
                globalColumn=fieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localColumn)
                CALL List_ItemAdd(columnIndicesLists(localDOFIdx)%ptr,globalColumn,err,error,*999)
              CASE(FIELD_NODE_BASED_INTERPOLATION)
                ! loop over all nodes in the element (and dofs belonging to them)
                DO localNodeIdx=1,basis2%NUMBER_OF_NODES
                  node2=domainElements2%elements(element)%ELEMENT_NODES(localNodeIdx)
                  DO derivativeIdx=1,basis2%NUMBER_OF_DERIVATIVES(localNodeIdx)
                    derivative=domainElements2%elements(element)%ELEMENT_DERIVATIVES(derivativeIdx,localNodeIdx)
                    version=domainElements2%elements(element)%elementVersions(derivativeIdx,localNodeIdx)
                    !Find the local and global column and add the global column to the indices list
                    localColumn=fieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                      & nodes(node2)%derivatives(derivative)%versions(version)
                    globalColumn=fieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localColumn)
                    CALL List_ItemAdd(columnIndicesLists(localDOFIdx)%ptr,globalColumn,err,error,*999)
                  ENDDO !derivativeIdx
                ENDDO !localNodeIdx
              CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                CALL FlagError("Grid point based interpolation is not implemented yet.",err,error,*999)
              CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                CALL FlagError("Gauss point based interpolation is not implemented yet.",err,error,*999)
              CASE DEFAULT
                localError="Local dof number "//TRIM(NumberToVString(localDOFIdx,"*",err,error))// &
                  & " has invalid interpolation type."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDDO !componentIdx
            !Clean up the list
            CALL List_RemoveDuplicates(columnIndicesLists(localDOFIdx)%ptr,err,error,*999)
            CALL List_NumberOfItemsGet(columnIndicesLists(localDOFIdx)%ptr,numberOfColumns,err,error,*999)
            numberOfNonZeros=numberOfNonZeros+numberOfColumns
            rowIndices(localDOFIdx+1)=numberOfNonZeros+1
          CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
            CALL FlagError("Grid point based interpolation is not implemented yet.",err,error,*999)
          CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
            CALL FlagError("Gauss point based interpolation is not implemented yet.",err,error,*999)
          CASE DEFAULT
            localError="Local dof number "//TRIM(NumberToVString(localDOFIdx,"*",err,error))//" has an invalid type."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDDO !localDOFIdx
       CASE DEFAULT
        localError="The matrix storage type of "//TRIM(NumberToVString(jacobianMatrix%storageType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT      
    CASE(EQUATIONS_MATRIX_NODAL_STRUCTURE)
      SELECT CASE(jacobianMatrix%storageType)
      CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        !Allocate lists
        ALLOCATE(columnIndicesLists(rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate column indices lists.",err,error,*999)
        !Allocate row indices
        ALLOCATE(rowIndices(rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL+1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate row indices.",err,error,*999)
        rowIndices(1)=1
        !First, loop over the rows and calculate the number of non-zeros
        numberOfNonZeros=0
        DO localDofIdx=1,rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
          SELECT CASE(rowDofsParamMapping%DOF_TYPE(1,localDofIdx))
          CASE(FIELD_CONSTANT_INTERPOLATION)
            CALL FlagError("Constant interpolation is not implemented yet.",err,error,*999)
          CASE(FIELD_NODE_DOF_TYPE)
            dofIdx=dependentDofsParamMapping%DOF_TYPE(2,localDofIdx)!value for a particular field dof (localDofIdx)
            node=dependentDofsParamMapping%NODE_DOF2PARAM_MAP(3,dofIdx)!node number (node) of the field parameter
            component=dependentDofsParamMapping%NODE_DOF2PARAM_MAP(4,dofIdx)!component number (component) of the field parameter
            domainNodes=>fieldVariable%components(component)%domain%topology%NODES            
            !Set up list
            NULLIFY(columnIndicesLists(localDofIdx)%ptr)
            CALL List_CreateStart(columnIndicesLists(localDofIdx)%ptr,err,error,*999)
            CALL List_DataTypeSet(columnIndicesLists(localDofIdx)%ptr,LIST_INTG_TYPE,err,error,*999)            
            CALL List_InitialSizeSet(columnIndicesLists(localDofIdx)%ptr,fieldVariable%NUMBER_OF_COMPONENTS* &
              & fieldVariable%maxNumberElementInterpolationParameters,err,error,*999)            
            CALL List_CreateFinish(columnIndicesLists(localDofIdx)%ptr,err,error,*999)
            !Loop over all components,nodes,derivatives, and versions
            DO componentIdx=1,fieldVariable%NUMBER_OF_COMPONENTS
              SELECT CASE(fieldVariable%components(componentIdx)%INTERPOLATION_TYPE)
              CASE(FIELD_NODE_BASED_INTERPOLATION)
                numberOfDerivatives=fieldVariable%components(componentIdx)%domain%topology%nodes%nodes(node)%NUMBER_OF_DERIVATIVES
                DO derivativeIdx=1,numberOfDerivatives
                  numberOfVersions=fieldVariable%components(componentIdx)%domain%topology%nodes%nodes(node)% &
                    & derivatives(derivativeIdx)%numberOfVersions
                  DO versionIdx=1,numberOfVersions
                    localColumn=fieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%nodes(node)% &
                      & derivatives(derivativeIdx)%versions(versionIdx)
                    globalColumn=fieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localColumn)                    
                    CALL List_ItemAdd(columnIndicesLists(localDofIdx)%ptr,globalColumn,err,error,*999)                    
                  ENDDO !versionIdx
                ENDDO !derivativeIdx
              CASE DEFAULT
                localError="Local dof number "//TRIM(NumberToVString(localDofIdx,"*",err,error))//" has invalid interpolation type."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDDO !componentIdx            
            CALL List_RemoveDuplicates(columnIndicesLists(localDofIdx)%ptr,err,error,*999)
            CALL List_NumberOfItemsGet(columnIndicesLists(localDofIdx)%ptr,numberOfColumns,err,error,*999)
            numberOfNonZeros=numberOfNonZeros+numberOfColumns
            rowIndices(localDofIdx+1)=numberOfNonZeros+1
          CASE(FIELD_ELEMENT_DOF_TYPE)
            CALL FlagError("Element based interpolation is not implemented yet.",err,error,*999)
          CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
            CALL FlagError("Grid point based interpolation is not implemented yet.",err,error,*999)
          CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
            CALL FlagError("Gauss point based interpolation is not implemented yet.",err,error,*999)
          CASE DEFAULT
            localError="Local dof number "//TRIM(NumberToVString(localDofIdx,"*",err,error))//" has an invalid type."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDDO !localDofIdx
      CASE DEFAULT
        localError="The matrix storage type of "//TRIM(NumberToVString(jacobianMatrix%storageType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The matrix structure type of "//TRIM(NumberToVString(jacobianMatrix%structureType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*998)
    END SELECT

    IF(numberOfNonZeros>0) THEN
      !Allocate and setup the column locations
      ALLOCATE(columnIndices(numberOfNonZeros),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate column indices.",err,error,*999)
      DO localDOFIdx=1,rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
        CALL LIST_DETACH_AND_DESTROY(columnIndicesLists(localDOFIdx)%ptr,numberOfColumns,columns,err,error,*999)
        DO columnIdx=1,numberOfColumns
          columnIndices(rowIndices(localDOFIdx)+columnIdx-1)=columns(columnIdx)
        ENDDO !columnIdx
        DEALLOCATE(columns)
      ENDDO !localDOFIdx
    ENDIF
    
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Jacobian matrix structure:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ",dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL, &
        & err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ",dependentDofsDomainMapping%NUMBER_OF_GLOBAL, &
        & err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of non zeros = ",numberOfNonZeros,err,error,*999)
      IF(dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL*dependentDofsDomainMapping%NUMBER_OF_GLOBAL/=0) THEN
        sparsity=(1.0_DP-REAL(numberOfNonZeros,DP)/REAL(dependentDofsDomainMapping% &
          & TOTAL_NUMBER_OF_LOCAL*dependentDofsDomainMapping%NUMBER_OF_GLOBAL,DP))*100.0_DP
        CALL WriteStringFmtValue(DIAGNOSTIC_OUTPUT_TYPE,"  Sparsity (% of zeros) = ",sparsity,"F6.2",err,error,*999)
      ENDIF
      IF(numberOfNonZeros>0) THEN
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL+1,8,8, &
          & rowIndices,'("  Row indices    :",8(X,I13))','(18X,8(X,I13))',err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfNonZeros,8,8,columnIndices,&
          & '("  Column indices :",8(X,I13))','(18X,8(X,I13))', err,error,*999)
      ENDIF
    ENDIF
      
    EXITS("JacobianMatrix_StructureCalculate")
    RETURN
999 IF(ASSOCIATED(rowIndices)) DEALLOCATE(rowIndices)
    IF(ASSOCIATED(columnIndices)) DEALLOCATE(columnIndices)
    IF(ALLOCATED(columns)) DEALLOCATE(columns)
    IF(ALLOCATED(columnIndicesLists)) THEN
      DO localDofIdx=1,dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
        IF(ASSOCIATED(columnIndicesLists(localDofIdx)%ptr)) &
          & CALL List_Destroy(columnIndicesLists(localDofIdx)%ptr,dummyErr,dummyError,*998)
      ENDDO !localDofIdx
      DEALLOCATE(columnIndicesLists)
    ENDIF
998 ERRORSEXITS("JacobianMatrix_StructureCalculate",err,error)
    RETURN 1
  END SUBROUTINE JacobianMatrix_StructureCalculate

  !
  !================================================================================================================================
  !
 
END MODULE EquationsMatricesRoutines
