!> \file
!> \author Chris Bradley
!> \brief This module contains all equations matrices access method routines.
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

!> This module contains all equations matrices access method routines.
MODULE EquationsMatricesAccessRoutines
  
  USE BaseRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE Strings
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup EquationsMatricesRoutines_MatrixTypes EquationsMatricesRoutines::MatrixTypes
  !> \brief Equations matrix types
  !> \see EquationsMatricesRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRIX_LINEAR=1 !<The euations matrix is a linear matrix \see EquationsMatricesRoutines_MatrixTypes,EquationsMatricesRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRIX_DYNAMIC=2 !<The equations matrix is a dynamic matrix \see EquationsMatricesRoutines_MatrixTypes,EquationsMatricesRoutines
 !>@}
 
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
  !> \see SolverRoutines
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
  
  !> \addtogroup EquationsMatricesRoutines_VectorTemporalTypes EquationsMatricesRoutines::VectorTemporalTypes
  !> \brief The temporal types of equations matrices vectors
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_CURRENT_VECTOR=1 !<The current equations matrices vector
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_PREVIOUS_VECTOR=2 !<The previous equations matrices vector
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_PREVIOUS2_VECTOR=3 !<The second previous equations matrices vector
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_PREVIOUS3_VECTOR=4 !<The third previous equations matrices vector
  !>@}
 
  !Module types

  !Module variables

  !Interfaces

  PUBLIC EQUATIONS_MATRIX_LINEAR,EQUATIONS_MATRIX_DYNAMIC

  PUBLIC EQUATIONS_MATRIX_NO_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE, &
    & EQUATIONS_MATRIX_NODAL_STRUCTURE

  PUBLIC EQUATIONS_MATRIX_UNLUMPED,EQUATIONS_MATRIX_LUMPED

  PUBLIC EQUATIONS_MATRICES_SPARSE_MATRICES,EQUATIONS_MATRICES_FULL_MATRICES

  PUBLIC EQUATIONS_MATRICES_ALL,EQUATIONS_MATRICES_DYNAMIC_ONLY,EQUATIONS_MATRICES_LINEAR_ONLY, &
    & EQUATIONS_MATRICES_NONLINEAR_ONLY,EQUATIONS_MATRICES_JACOBIAN_ONLY,EQUATIONS_MATRICES_RESIDUAL_ONLY, &
    & EQUATIONS_MATRICES_RHS_ONLY,EQUATIONS_MATRICES_SOURCE_ONLY,EQUATIONS_MATRICES_RHS_RESIDUAL_ONLY, &
    & EQUATIONS_MATRICES_RHS_SOURCE_ONLY,EQUATIONS_MATRICES_RESIDUAL_SOURCE_ONLY,EQUATIONS_MATRICES_VECTORS_ONLY
  
  PUBLIC EQUATIONS_GRADIENT_FINITE_DIFFERENCE_CALCULATED,EQUATIONS_GRADIENT_ANALYTIC_CALCULATED

  PUBLIC EQUATIONS_JACOBIAN_FINITE_DIFFERENCE_CALCULATED,EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED

  PUBLIC EQUATIONS_HESSIAN_FINITE_DIFFERENCE_CALCULATED,EQUATIONS_HESSIAN_ANALYTIC_CALCULATED

  PUBLIC EQUATIONS_MATRICES_CURRENT_VECTOR,EQUATIONS_MATRICES_PREVIOUS_VECTOR,EQUATIONS_MATRICES_PREVIOUS2_VECTOR, &
    & EQUATIONS_MATRICES_PREVIOUS3_VECTOR

  PUBLIC ElementMatrix_Output
  
  PUBLIC ElementVector_Output
  
  PUBLIC EquationsMatricesDynamic_DampingMatrixGet
  
  PUBLIC EquationsMatricesDynamic_DynamicMappingGet
  
  PUBLIC EquationsMatricesDynamic_EquationsMatrixGet

  PUBLIC EquationsMatricesDynamic_MassMatrixGet
  
  PUBLIC EquationsMatricesDynamic_NumberOfDynamicMatricesGet

  PUBLIC EquationsMatricesDynamic_StiffnessMatrixGet
  
  PUBLIC EquationsMatricesDynamic_TempDistributedVectorExists
  
  PUBLIC EquationsMatricesDynamic_TempDistributedVectorGet
  
  PUBLIC EquationsMatricesDynamic_VectorMatricesGet

  PUBLIC EquationsMatricesLinear_EquationsMatrixGet
  
  PUBLIC EquationsMatricesLinear_LinearMappingGet

  PUBLIC EquationsMatricesLinear_NumberOfLinearMatricesGet

  PUBLIC EquationsMatricesLinear_TempDistributedVectorExists

  PUBLIC EquationsMatricesLinear_TempDistributedVectorGet

  PUBLIC EquationsMatricesLinear_VectorMatricesGet

  PUBLIC EquationsMatricesNonlinear_NonlinearMappingGet
  
  PUBLIC EquationsMatricesNonlinear_NumberOfResidualsGet

  PUBLIC EquationsMatricesNonlinear_ResidualVectorGet
  
  PUBLIC EquationsMatricesNonlinear_TempDistributedVectorExists

  PUBLIC EquationsMatricesNonlinear_TempDistributedVectorGet

  PUBLIC EquationsMatricesNonlinear_VectorMatricesGet

  PUBLIC EquationsMatricesOptimisation_HessianMatrixGet

  PUBLIC EquationsMatricesResidual_DistributedVectorExists

  PUBLIC EquationsMatricesResidual_DistributedVectorGet

  PUBLIC EquationsMatricesResidual_ElementVectorOutput

  PUBLIC EquationsMatricesResidual_FirstAssemblyGet

  PUBLIC EquationsMatricesResidual_CurrentDistributedVectorExists

  PUBLIC EquationsMatricesResidual_CurrentDistributedVectorGet

  PUBLIC EquationsMatricesResidual_PreviousDistributedVectorExists

  PUBLIC EquationsMatricesResidual_PreviousDistributedVectorGet

  PUBLIC EquationsMatricesResidual_Previous2DistributedVectorExists

  PUBLIC EquationsMatricesResidual_Previous2DistributedVectorGet

  PUBLIC EquationsMatricesResidual_Previous3DistributedVectorExists

  PUBLIC EquationsMatricesResidual_Previous3DistributedVectorGet

  PUBLIC EquationsMatricesResidual_JacobianMatrixGet

  PUBLIC EquationsMatricesResidual_NodalVectorOutput

  PUBLIC EquationsMatricesResidual_NonlinearMatricesGet

  PUBLIC EquationsMatricesResidual_NumberOfJacobiansGet

  PUBLIC EquationsMatricesResidual_ResidualMappingGet

  PUBLIC EquationsMatricesResidual_ResidualNumberGet

  PUBLIC EquationsMatricesResidual_VectorCoefficientGet

  PUBLIC EquationsMatricesResidual_UpdateVectorGet

  PUBLIC EquationsMatricesResidual_UpdateVectorSet

  PUBLIC EquationsMatricesRHS_DistributedVectorExists
  
  PUBLIC EquationsMatricesRHS_DistributedVectorGet

  PUBLIC EquationsMatricesRHS_ElementVectorOutput

  PUBLIC EquationsMatricesRHS_FirstAssemblyGet

  PUBLIC EquationsMatricesRHS_CurrentDistributedVectorExists
  
  PUBLIC EquationsMatricesRHS_CurrentDistributedVectorGet

  PUBLIC EquationsMatricesRHS_PreviousDistributedVectorExists
  
  PUBLIC EquationsMatricesRHS_PreviousDistributedVectorGet

  PUBLIC EquationsMatricesRHS_Previous2DistributedVectorExists
  
  PUBLIC EquationsMatricesRHS_Previous2DistributedVectorGet

  PUBLIC EquationsMatricesRHS_Previous3DistributedVectorExists
  
  PUBLIC EquationsMatricesRHS_Previous3DistributedVectorGet

  PUBLIC EquationsMatricesRHS_NodalVectorOutput

  PUBLIC EquationsMatricesRHS_VectorMatricesGet

  PUBLIC EquationsMatricesRHS_VectorCoefficientGet

  PUBLIC EquationsMatricesRHS_UpdateVectorGet

  PUBLIC EquationsMatricesRHS_UpdateVectorSet

  PUBLIC EquationsMatricesScalar_EquationsScalarGet

  PUBLIC EquationsMatricesSource_DistributedVectorExists

  PUBLIC EquationsMatricesSource_DistributedVectorGet

  PUBLIC EquationsMatricesSource_ElementVectorOutput

  PUBLIC EquationsMatricesSource_FirstAssemblyGet

  PUBLIC EquationsMatricesSource_CurrentDistributedVectorExists

  PUBLIC EquationsMatricesSource_CurrentDistributedVectorGet

  PUBLIC EquationsMatricesSource_PreviousDistributedVectorExists

  PUBLIC EquationsMatricesSource_PreviousDistributedVectorGet

  PUBLIC EquationsMatricesSource_Previous2DistributedVectorExists

  PUBLIC EquationsMatricesSource_Previous2DistributedVectorGet

  PUBLIC EquationsMatricesSource_Previous3DistributedVectorExists

  PUBLIC EquationsMatricesSource_Previous3DistributedVectorGet

  PUBLIC EquationsMatricesSource_NodalVectorOutput

  PUBLIC EquationsMatricesSource_SourceNumberGet

  PUBLIC EquationsMatricesSource_SourceVectorsGet

  PUBLIC EquationsMatricesSource_VectorCoefficientGet

  PUBLIC EquationsMatricesSource_UpdateVectorGet

  PUBLIC EquationsMatricesSource_UpdateVectorSet

  PUBLIC EquationsMatricesSources_NumberOfSourcesGet

  PUBLIC EquationsMatricesSources_SourceVectorGet

  PUBLIC EquationsMatricesSources_TempDistributedVectorExists

  PUBLIC EquationsMatricesSources_TempDistributedVectorGet

  PUBLIC EquationsMatricesSources_VectorMatricesGet

  PUBLIC EquationsMatricesVector_AssertIsFinished,EquationsMatricesVector_AssertNotFinished

  PUBLIC EquationsMatricesVector_DynamicMatricesExists

  PUBLIC EquationsMatricesVector_DynamicMatricesGet

  PUBLIC EquationsMatricesVector_LinearMatricesExists
  
  PUBLIC EquationsMatricesVector_LinearMatricesGet

  PUBLIC EquationsMatricesVector_NonlinearMatricesExists
  
  PUBLIC EquationsMatricesVector_NonlinearMatricesGet

  PUBLIC EquationsMatricesVector_OptimisationMatricesExists

  PUBLIC EquationsMatricesVector_OptimisationMatricesGet

  PUBLIC EquationsMatricesVector_NumberOfRowsGet

  PUBLIC EquationsMatricesVector_NumberOfGlobalRowsGet

  PUBLIC EquationsMatricesVector_RHSVectorExists
  
  PUBLIC EquationsMatricesVector_RHSVectorGet

  PUBLIC EquationsMatricesVector_SourceVectorsExists
  
  PUBLIC EquationsMatricesVector_SourceVectorsGet

  PUBLIC EquationsMatricesVector_TotalNumberOfRowsGet

  PUBLIC EquationsMatricesVector_VectorEquationsGet

  PUBLIC EquationsMatricesVector_VectorMappingGet

  PUBLIC EquationsMatrix_DistributedMatrixGet

  PUBLIC EquationsMatrix_DynamicMatricesExists

  PUBLIC EquationsMatrix_DynamicMatricesGet

  PUBLIC EquationsMatrix_ElementMatrixOutput

  PUBLIC EquationsMatrix_FirstAssemblyGet

  PUBLIC EquationsMatrix_LinearMatricesExists
  
  PUBLIC EquationsMatrix_LinearMatricesGet

  PUBLIC EquationsMatrix_LumpedFlagGet

  PUBLIC EquationsMatrix_MatrixCoefficientGet

  PUBLIC EquationsMatrix_MatrixNumberGet

  PUBLIC EquationsMatrix_NodalMatrixOutput

  PUBLIC EquationsMatrix_NumberOfColumnsGet
  
  PUBLIC EquationsMatrix_StorageTypeGet

  PUBLIC EquationsMatrix_StructureTypeGet

  PUBLIC EquationsMatrix_SymmetryFlagGet

  PUBLIC EquationsMatrix_TempDistributedVectorExists

  PUBLIC EquationsMatrix_TempDistributedVectorGet

  PUBLIC EquationsMatrix_UpdateMatrixGet

  PUBLIC EquationsMatrix_UpdateMatrixSet

  PUBLIC JacobianMatrix_CalculationTypeGet

  PUBLIC JacobianMatrix_DistributedMatrixGet

  PUBLIC JacobianMatrix_ElementMatrixOutput

  PUBLIC JacobianMatrix_FiniteDifferenceStepSizeGet

  PUBLIC JacobianMatrix_FirstAssemblyGet

  PUBLIC JacobianMatrix_MatrixCoefficientGet

  PUBLIC JacobianMatrix_MatrixNumberGet

  PUBLIC JacobianMatrix_NodalMatrixOutput

  PUBLIC JacobianMatrix_NumberOfColumnsGet
  
  PUBLIC JacobianMatrix_ResidualVectorGet

  PUBLIC JacobianMatrix_StorageTypeGet

  PUBLIC JacobianMatrix_StructureTypeGet

  PUBLIC JacobianMatrix_SymmetryFlagGet

  PUBLIC JacobianMatrix_UpdateMatrixGet

  PUBLIC JacobianMatrix_UpdateMatrixSet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Outputs the element matrix information.
  SUBROUTINE ElementMatrix_Output(id,elementMatrix,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the ouptut stream
    TYPE(ElementMatrixType), INTENT(IN) :: elementMatrix !<The element matrix to output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("ElementMatrix_Output",err,error,*999)

    CALL WriteStringValue(id,"  Number of rows = ",elementMatrix%numberOfRows,err,error,*999)
    CALL WriteStringValue(id,"  Number of columns = ",elementMatrix%numberOfColumns,err,error,*999)
    CALL WriteStringValue(id,"  Maximum number of rows = ",elementMatrix%maxNumberOfRows,err,error,*999)
    CALL WriteStringValue(id,"  Maximum number of columns = ",elementMatrix%maxNumberOfColumns,err,error,*999)
    CALL WriteStringVector(id,1,1,elementMatrix%numberOfRows,8,8,elementMatrix%rowDOFS, &
      & '("  Row dofs         :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
    CALL WriteStringVector(id,1,1,elementMatrix%numberOfColumns,8,8,elementMatrix%columnDOFS, &
      & '("  Column dofs      :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
    CALL WriteStringMatrix(id,1,1,elementMatrix%numberOfRows,1,1,elementMatrix%numberOfColumns,8,8, &
      & elementMatrix%matrix(1:elementMatrix%numberOfRows,1:elementMatrix%numberOfColumns),WRITE_STRING_MATRIX_NAME_AND_INDICES, &
      & '("  Matrix','(",I2,",:)','     :",8(X,E13.6))','(20X,8(X,E13.6))',err,error,*999)
    
    EXITS("ElementMatrix_Output")
    RETURN
999 ERRORSEXITS("ElementMatrix_Output",err,error)
    RETURN 1
    
  END SUBROUTINE ElementMatrix_Output
  
  !
  !================================================================================================================================
  !

  !>Outputs the element vector information .
  SUBROUTINE ElementVector_Output(id,elementVector,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the ouptut stream
    TYPE(ElementVectorType), INTENT(IN) :: elementVector !<The element vector to output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("ElementVector_Output",err,error,*999)

    CALL WriteStringValue(id,"  Number of rows = ",elementVector%numberOfRows,err,error,*999)
    CALL WriteStringValue(id,"  Maximum number of rows = ",elementVector%maxNumberOfRows,err,error,*999)
    CALL WriteStringVector(id,1,1,elementVector%numberOfRows,8,8,elementVector%rowDOFS, &
      & '("  Row dofs         :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
    CALL WriteStringVector(id,1,1,elementVector%numberOfRows,8,8,elementVector%vector, &
      & '("  Vector(:)        :",8(X,E13.6))','(20X,8(X,E13.6))',err,error,*999)
    
    EXITS("ElementVector_Output")
    RETURN
999 ERRORSEXITS("ElementVector_Output",err,error)
    RETURN 1
    
  END SUBROUTINE ElementVector_Output
  
  !
  !================================================================================================================================
  !

  !>Gets the dynamic damping matrix for dynamic matrices.
  SUBROUTINE EquationsMatricesDynamic_DampingMatrixGet(dynamicMatrices,dampingMatrix,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices !<A pointer to the dynamic matrices to get the damping matrix for
    TYPE(EquationsMatrixType), POINTER :: dampingMatrix !<On exit, a pointer to the damping matrix for the dynamic matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dampingMatrixNumber
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsMatricesDynamic_DampingMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dampingMatrix)) CALL FlagError("Damping matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dynamicMatrices)) CALL FlagError("Dynamic matrices is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(dynamicMatrices%vectorMatrices)) &
      & CALL FlagError("Vector matrices is not associated for the dynamic matrices.",err,error,*999)
    IF(.NOT.ALLOCATED(dynamicMatrices%matrices)) &
      & CALL FlagError("Matrices is not allocated for the dynamic matrices.",err,error,*999)
    IF(.NOT.ASSOCIATED(dynamicMatrices%vectorMatrices%vectorEquations)) &
      & CALL FlagError("Vector equations is not associated for the vector matrices.",err,error,*999)
    IF(.NOT.ASSOCIATED(dynamicMatrices%vectorMatrices%vectorEquations%vectorMapping)) &
      & CALL FlagError("Vector mapping is not associated for the vector equations.",err,error,*999)
    IF(.NOT.ASSOCIATED(dynamicMatrices%vectorMatrices%vectorEquations%vectorMapping%dynamicMapping)) &
      & CALL FlagError("Dynamic mapping is not associated for the vector mapping.",err,error,*999)
#endif    
    dampingMatrixNumber=dynamicMatrices%vectorMatrices%vectorEquations%vectorMapping%dynamicMapping%dampingMatrixNumber
#ifdef WITH_PRECHECKS
    IF(dampingMatrixNumber==0) CALL FlagError("The dynamic matrices have no damping matrix.",err,error,*999)
    IF(dampingMatrixNumber<0.OR.dampingMatrixNumber>dynamicMatrices%numberOfDynamicMatrices) THEN
      localError="The damping matrix number of "//TRIM(NumberToVString(dampingMatrixNumber,"*",err,error))// &
        & " is invalid. The damping matrix number should be >= 0 and <= "// &
        & TRIM(NumberToVString(dynamicMatrices%numberOfDynamicMatrices,"*",err,error))//"."
    ENDIF
#endif    

    dampingMatrix=>dynamicMatrices%matrices(dampingMatrixNumber)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dampingMatrix)) &
      & CALL FlagError("Damping matrix is not associated for the dynamic matrices.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesDynamic_DampingMatrixGet")
    RETURN
999 NULLIFY(dampingMatrix)
998 ERRORSEXITS("EquationsMatricesDynamic_DampingMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesDynamic_DampingMatrixGet

  !
  !================================================================================================================================
  !

  !>Gets the dynamic mapping for dynamic matrices.
  SUBROUTINE EquationsMatricesDynamic_DynamicMappingGet(dynamicMatrices,dynamicMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices !<A pointer to the dynamic matrices to get the dynamic mappipng for
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping !<On exit, a pointer to the dynamic mapping for the dynamic matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesDynamic_DynamicMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dynamicMapping)) CALL FlagError("Vector matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dynamicMatrices)) CALL FlagError("Dynamic matrices is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(dynamicMatrices%vectorMatrices)) &
      & CALL FlagError("Vector matrices is not associated for the dynamic matrices.",err,error,*999)
    IF(.NOT.ASSOCIATED(dynamicMatrices%vectorMatrices%vectorEquations)) &
      & CALL FlagError("Vector equations is not associated for the vector matrices.",err,error,*999)
    IF(.NOT.ASSOCIATED(dynamicMatrices%vectorMatrices%vectorEquations%vectorMapping)) &
      & CALL FlagError("Vector mapping is not associated for the vector equations.",err,error,*999)
#endif    

    dynamicMapping=>dynamicMatrices%vectorMatrices%vectorEquations%vectorMapping%dynamicMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dynamicMapping)) &
      & CALL FlagError("Dynamic mapping is not associated for the dynamic matrices.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesDynamic_DynamicMappingGet")
    RETURN
999 NULLIFY(dynamicMapping)
998 ERRORSEXITS("EquationsMatricesDynamic_DynamicMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesDynamic_DynamicMappingGet

  !
  !================================================================================================================================
  !

  !>Gets the specified equations matrix for dynamic matrices.
  SUBROUTINE EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,equationsMatrix,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices !<A pointer to the dynamic matrices to get the equations matrix for
    INTEGER(INTG), INTENT(IN) :: matrixIdx !<The matrix index of the equations matrix to get
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<On exit, a pointer to the equations matrix for the matrixIdx'th dynamic matrix. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsMatricesDynamic_EquationsMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dynamicMatrices)) CALL FlagError("Dynamic matrices is not associated.",err,error,*999)
    IF(matrixIdx<1.OR.matrixIdx>dynamicMatrices%numberOfDynamicMatrices) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is invalid. The matrix index must be >= 1 and <= "// &
        & TRIM(NumberToVString(dynamicMatrices%numberOfDynamicMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(dynamicMatrices%matrices)) CALL FlagError("Dynamic matrices matrices is not allocated.",err,error,*999)
#endif    
    
    equationsMatrix=>dynamicMatrices%matrices(matrixIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equationsMatrix)) THEN
      localError="Equations matrix for dynamic matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    EXITS("EquationsMatricesDynamic_EquationsMatrixGet")
    RETURN
999 NULLIFY(equationsMatrix)
998 ERRORSEXITS("EquationsMatricesDynamic_EquationsMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesDynamic_EquationsMatrixGet

  !
  !================================================================================================================================
  !

  !>Gets the dynamic mass matrix for dynamic matrices.
  SUBROUTINE EquationsMatricesDynamic_MassMatrixGet(dynamicMatrices,massMatrix,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices !<A pointer to the dynamic matrices to get the mass matrix for
    TYPE(EquationsMatrixType), POINTER :: massMatrix !<On exit, a pointer to the mass matrix for the dynamic matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: massMatrixNumber
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsMatricesDynamic_MassMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(massMatrix)) CALL FlagError("Mass matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dynamicMatrices)) CALL FlagError("Dynamic matrices is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(dynamicMatrices%vectorMatrices)) &
      & CALL FlagError("Vector matrices is not associated for the dynamic matrices.",err,error,*999)
    IF(.NOT.ALLOCATED(dynamicMatrices%matrices)) &
      & CALL FlagError("Matrices is not allocated for the dynamic matrices.",err,error,*999)
    IF(.NOT.ASSOCIATED(dynamicMatrices%vectorMatrices%vectorEquations)) &
      & CALL FlagError("Vector equations is not associated for the vector matrices.",err,error,*999)
    IF(.NOT.ASSOCIATED(dynamicMatrices%vectorMatrices%vectorEquations%vectorMapping)) &
      & CALL FlagError("Vector mapping is not associated for the vector equations.",err,error,*999)
    IF(.NOT.ASSOCIATED(dynamicMatrices%vectorMatrices%vectorEquations%vectorMapping%dynamicMapping)) &
      & CALL FlagError("Dynamic mapping is not associated for the vector mapping.",err,error,*999)
#endif    
    massMatrixNumber=dynamicMatrices%vectorMatrices%vectorEquations%vectorMapping%dynamicMapping%massMatrixNumber
#ifdef WITH_PRECHECKS
    IF(massMatrixNumber==0) CALL FlagError("The dynamic matrices have no mass matrix.",err,error,*999)
    IF(massMatrixNumber<0.OR.massMatrixNumber>dynamicMatrices%numberOfDynamicMatrices) THEN
      localError="The mass matrix number of "//TRIM(NumberToVString(massMatrixNumber,"*",err,error))// &
        & " is invalid. The mass matrix number should be >= 0 and <= "// &
        & TRIM(NumberToVString(dynamicMatrices%numberOfDynamicMatrices,"*",err,error))//"."
    ENDIF
#endif    

    massMatrix=>dynamicMatrices%matrices(massMatrixNumber)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(massMatrix)) &
      & CALL FlagError("Mass matrix is not associated for the dynamic matrices.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesDynamic_MassMatrixGet")
    RETURN
999 NULLIFY(massMatrix)
998 ERRORSEXITS("EquationsMatricesDynamic_MassMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesDynamic_MassMatrixGet

  !
  !================================================================================================================================
  !

  !>Returns the number of dynamic matrices for dynamic matrices.
  SUBROUTINE EquationsMatricesDynamic_NumberOfDynamicMatricesGet(dynamicMatrices,numberOfDynamicMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices !<A pointer to the dynamic matrices to get the number of dynamic matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfDynamicMatrices !<On exit, the number of dynamic matrices for the dynamic matrices.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesDynamic_NumberOfDynamicMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dynamicMatrices)) CALL FlagError("Dynamic matrices is not associated.",err,error,*999)
#endif    

    numberOfDynamicMatrices=dynamicMatrices%numberOfDynamicMatrices

    EXITS("EquationsMatricesDynamic_NumberOfDynamicMatricesGet")
    RETURN
999 ERRORSEXITS("EquationsMatricesDynamic_NumberOfDynamicMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesDynamic_NumberOfDynamicMatricesGet

  !
  !================================================================================================================================
  !

  !>Gets the dynamic stiffness matrix for dynamic matrices.
  SUBROUTINE EquationsMatricesDynamic_StiffnessMatrixGet(dynamicMatrices,stiffnessMatrix,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices !<A pointer to the dynamic matrices to get the stiffness matrix for
    TYPE(EquationsMatrixType), POINTER :: stiffnessMatrix !<On exit, a pointer to the stiffness matrix for the dynamic matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: stiffnessMatrixNumber
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsMatricesDynamic_StiffnessMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(stiffnessMatrix)) CALL FlagError("Stiffness matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dynamicMatrices)) CALL FlagError("Dynamic matrices is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(dynamicMatrices%vectorMatrices)) &
      & CALL FlagError("Vector matrices is not associated for the dynamic matrices.",err,error,*999)
    IF(.NOT.ALLOCATED(dynamicMatrices%matrices)) &
      & CALL FlagError("Matrices is not allocated for the dynamic matrices.",err,error,*999)
    IF(.NOT.ASSOCIATED(dynamicMatrices%vectorMatrices%vectorEquations)) &
      & CALL FlagError("Vector equations is not associated for the vector matrices.",err,error,*999)
    IF(.NOT.ASSOCIATED(dynamicMatrices%vectorMatrices%vectorEquations%vectorMapping)) &
      & CALL FlagError("Vector mapping is not associated for the vector equations.",err,error,*999)
    IF(.NOT.ASSOCIATED(dynamicMatrices%vectorMatrices%vectorEquations%vectorMapping%dynamicMapping)) &
      & CALL FlagError("Dynamic mapping is not associated for the vector mapping.",err,error,*999)
#endif    
    stiffnessMatrixNumber=dynamicMatrices%vectorMatrices%vectorEquations%vectorMapping%dynamicMapping%stiffnessMatrixNumber
#ifdef WITH_PRECHECKS
    IF(stiffnessMatrixNumber==0) CALL FlagError("The dynamic matrices have no stiffness matrix.",err,error,*999)
    IF(stiffnessMatrixNumber<0.OR.stiffnessMatrixNumber>dynamicMatrices%numberOfDynamicMatrices) THEN
      localError="The stiffness matrix number of "//TRIM(NumberToVString(stiffnessMatrixNumber,"*",err,error))// &
        & " is invalid. The stiffness matrix number should be >= 0 and <= "// &
        & TRIM(NumberToVString(dynamicMatrices%numberOfDynamicMatrices,"*",err,error))//"."
    ENDIF
#endif    

    stiffnessMatrix=>dynamicMatrices%matrices(stiffnessMatrixNumber)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(stiffnessMatrix)) &
      & CALL FlagError("Stiffness matrix is not associated for the dynamic matrices.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesDynamic_StiffnessMatrixGet")
    RETURN
999 NULLIFY(stiffnessMatrix)
998 ERRORSEXITS("EquationsMatricesDynamic_StiffnessMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesDynamic_StiffnessMatrixGet

  !
  !================================================================================================================================
  !

  !>Checks if the temporary distributed vector for dynamic matrices exists.
  SUBROUTINE EquationsMatricesDynamic_TempDistributedVectorExists(dynamicMatrices,tempDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices !<A pointer to the dynamic matrices to check the temp distributed vector exists for
    TYPE(DistributedVectorType), POINTER :: tempDistributedVector !<On exit, a pointer to the temp distributed vector for the dynamic matrices if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesDynamic_TempDistributedVectorExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(tempDistributedVector)) CALL FlagError("Temporary distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dynamicMatrices)) CALL FlagError("Dynamic matrices is not associated.",err,error,*999)
#endif    

    tempDistributedVector=>dynamicMatrices%tempVector

    EXITS("EquationsMatricesDynamic_TempDistributedVectorExists")
    RETURN
999 NULLIFY(tempDistributedVector)
998 ERRORSEXITS("EquationsMatricesDynamic_TempDistributedVectorExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesDynamic_TempDistributedVectorExists

  !
  !================================================================================================================================
  !

  !>Gets the temporary distributed vector for dynamic matrices.
  SUBROUTINE EquationsMatricesDynamic_TempDistributedVectorGet(dynamicMatrices,tempDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices !<A pointer to the dynamic matrices to get the temp distributed vector for
    TYPE(DistributedVectorType), POINTER :: tempDistributedVector !<On exit, a pointer to the temp distributed vector for the dynamic matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesDynamic_TempDistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(tempDistributedVector)) CALL FlagError("Temporary distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dynamicMatrices)) CALL FlagError("Dynamic matrices is not associated.",err,error,*999)
#endif    

    tempDistributedVector=>dynamicMatrices%tempVector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(tempDistributedVector)) &
      & CALL FlagError("Temporary distributed vector is not associated for the dynamic matrices.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesDynamic_TempDistributedVectorGet")
    RETURN
999 NULLIFY(tempDistributedVector)
998 ERRORSEXITS("EquationsMatricesDynamic_TempDistributedVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesDynamic_TempDistributedVectorGet

  !
  !================================================================================================================================
  !

  !>Gets the vector matrices for dynamic matrices.
  SUBROUTINE EquationsMatricesDynamic_VectorMatricesGet(dynamicMatrices,vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices !<A pointer to the dynamic matrices to get the vector matrices for
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<On exit, a pointer to the vector matrices for the dynamic matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesDynamic_VectorMatricesGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dynamicMatrices)) CALL FlagError("Dynamic matrices is not associated.",err,error,*999)
#endif    

    vectorMatrices=>dynamicMatrices%vectorMatrices

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(vectorMatrices)) &
      & CALL FlagError("Vector matrices is not associated for the dynamic matrices.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesDynamic_VectorMatricesGet")
    RETURN
999 NULLIFY(vectorMatrices)
998 ERRORSEXITS("EquationsMatricesDynamic_VectorMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesDynamic_VectorMatricesGet

  !
  !================================================================================================================================
  !

  !>Gets the specified equations matrix for linear matrices.
  SUBROUTINE EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,equationsMatrix,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices !<A pointer to the linear matrices to get the equations matrix for
    INTEGER(INTG), INTENT(IN) :: matrixIdx !<The matrix index of the equations matrix to get
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<On exit, a pointer to the equations matrix for the matrixIdx'th linear matrix. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsMatricesLinear_EquationsMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(linearMatrices)) CALL FlagError("Linear matrices is not associated.",err,error,*999)
    IF(matrixIdx<1.OR.matrixIdx>linearMatrices%numberOfLinearMatrices) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is invalid. The matrix index must be >= 1 and <= "// &
        & TRIM(NumberToVString(linearMatrices%numberOfLinearMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(linearMatrices%matrices)) CALL FlagError("Linear matrices matrices is not allocated.",err,error,*999)
#endif    
    
    equationsMatrix=>linearMatrices%matrices(matrixIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equationsMatrix)) THEN
      localError="Equations matrix for linear matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    EXITS("EquationsMatricesLinear_EquationsMatrixGet")
    RETURN
999 NULLIFY(equationsMatrix)
998 ERRORSEXITS("EquationsMatricesLinear_EquationsMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesLinear_EquationsMatrixGet

  !
  !================================================================================================================================
  !

  !>Gets the linear mapping for linear matrices.
  SUBROUTINE EquationsMatricesLinear_LinearMappingGet(linearMatrices,linearMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices !<A pointer to the linear matrices to get the linear mappipng for
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping !<On exit, a pointer to the linear mapping for the linear matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesLinear_LinearMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linearMapping)) CALL FlagError("Vector matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(linearMatrices)) CALL FlagError("Linear matrices is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(linearMatrices%vectorMatrices)) &
      & CALL FlagError("Vector matrices is not associated for the linear matrices.",err,error,*999)
    IF(.NOT.ASSOCIATED(linearMatrices%vectorMatrices%vectorEquations)) &
      & CALL FlagError("Vector equations is not associated for the vector matrices.",err,error,*999)
    IF(.NOT.ASSOCIATED(linearMatrices%vectorMatrices%vectorEquations%vectorMapping)) &
      & CALL FlagError("Vector mapping is not associated for the vector equations.",err,error,*999)
#endif    

    linearMapping=>linearMatrices%vectorMatrices%vectorEquations%vectorMapping%linearMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(linearMapping)) &
      & CALL FlagError("Linear mapping is not associated for the linear matrices.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesLinear_LinearMappingGet")
    RETURN
999 NULLIFY(linearMapping)
998 ERRORSEXITS("EquationsMatricesLinear_LinearMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesLinear_LinearMappingGet

  !
  !================================================================================================================================
  !

  !>Returns the number of linear matrices for linear matrices.
  SUBROUTINE EquationsMatricesLinear_NumberOfLinearMatricesGet(linearMatrices,numberOfLinearMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices !<A pointer to the linear matrices to get the number of linear matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfLinearMatrices !<On exit, the number of linear matrices for the linear matrices.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesLinear_NumberOfLinearMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(linearMatrices)) CALL FlagError("Linear matrices is not associated.",err,error,*999)
#endif    

    numberOfLinearMatrices=linearMatrices%numberOfLinearMatrices

    EXITS("EquationsMatricesLinear_NumberOfLinearMatricesGet")
    RETURN
999 ERRORSEXITS("EquationsMatricesLinear_NumberOfLinearMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesLinear_NumberOfLinearMatricesGet

  !
  !================================================================================================================================
  !

  !>Checks if the temporary distributed vector for linear matrices exists.
  SUBROUTINE EquationsMatricesLinear_TempDistributedVectorExists(linearMatrices,tempDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices !<A pointer to the linear matrices to check the temp distributed vector exists for
    TYPE(DistributedVectorType), POINTER :: tempDistributedVector !<On exit, a pointer to the temp distributed vector for the linear matrices if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesLinear_TempDistributedVectorExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(tempDistributedVector)) CALL FlagError("Temporary distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(linearMatrices)) CALL FlagError("Linear matrices is not associated.",err,error,*999)
#endif    

    tempDistributedVector=>linearMatrices%tempVector

    EXITS("EquationsMatricesLinear_TempDistributedVectorExists")
    RETURN
999 NULLIFY(tempDistributedVector)
998 ERRORSEXITS("EquationsMatricesLinear_TempDistributedVectorExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesLinear_TempDistributedVectorExists

  !
  !================================================================================================================================
  !

  !>Gets the temporary distributed vector for linear matrices.
  SUBROUTINE EquationsMatricesLinear_TempDistributedVectorGet(linearMatrices,tempDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices !<A pointer to the linear matrices to get the temp distributed vector for
    TYPE(DistributedVectorType), POINTER :: tempDistributedVector !<On exit, a pointer to the temp distributed vector for the linear matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesLinear_TempDistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(tempDistributedVector)) CALL FlagError("Temporary distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(linearMatrices)) CALL FlagError("Linear matrices is not associated.",err,error,*999)
#endif    

    tempDistributedVector=>linearMatrices%tempVector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(tempDistributedVector)) &
      & CALL FlagError("Temporary distributed vector is not associated for the linear matrices.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesLinear_TempDistributedVectorGet")
    RETURN
999 NULLIFY(tempDistributedVector)
998 ERRORSEXITS("EquationsMatricesLinear_TempDistributedVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesLinear_TempDistributedVectorGet

  !
  !================================================================================================================================
  !

  !>Gets the vector matrices for linear matrices.
  SUBROUTINE EquationsMatricesLinear_VectorMatricesGet(linearMatrices,vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices !<A pointer to the linear matrices to get the vector matrices for
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<On exit, a pointer to the vector matrices for the linear matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesLinear_VectorMatricesGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(linearMatrices)) CALL FlagError("Linear matrices is not associated.",err,error,*999)
#endif    

    vectorMatrices=>linearMatrices%vectorMatrices

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(vectorMatrices)) &
      & CALL FlagError("Vector matrices is not associated for the linear matrices.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesLinear_VectorMatricesGet")
    RETURN
999 NULLIFY(vectorMatrices)
998 ERRORSEXITS("EquationsMatricesLinear_VectorMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesLinear_VectorMatricesGet

  !
  !================================================================================================================================
  !

  !>Gets the nonlinear mapping for nonlinear matrices.
  SUBROUTINE EquationsMatricesNonlinear_NonlinearMappingGet(nonlinearMatrices,nonlinearMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices !<A pointer to the nonlinear matrices to get the nonlinear mappipng for
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping !<On exit, a pointer to the nonlinear mapping for the nonlinear matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesNonlinear_NonlinearMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(nonlinearMapping)) CALL FlagError("Vector matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(nonlinearMatrices)) CALL FlagError("Nonlinear matrices is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(nonlinearMatrices%vectorMatrices)) &
      & CALL FlagError("Vector matrices is not associated for the nonlinear matrices.",err,error,*999)
    IF(.NOT.ASSOCIATED(nonlinearMatrices%vectorMatrices%vectorEquations)) &
      & CALL FlagError("Vector equations is not associated for the vector matrices.",err,error,*999)
    IF(.NOT.ASSOCIATED(nonlinearMatrices%vectorMatrices%vectorEquations%vectorMapping)) &
      & CALL FlagError("Vector mapping is not associated for the vector equations.",err,error,*999)
#endif    

    nonlinearMapping=>nonlinearMatrices%vectorMatrices%vectorEquations%vectorMapping%nonlinearMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(nonlinearMapping)) &
      & CALL FlagError("Nonlinear mapping is not associated for the nonlinear matrices.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesNonlinear_NonlinearMappingGet")
    RETURN
999 NULLIFY(nonlinearMapping)
998 ERRORS("EquationsMatricesNonlinear_NonlinearMappingGet",err,error)
    EXITS("EquationsMatricesNonlinear_NonlinearMappingGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesNonlinear_NonlinearMappingGet

  !
  !================================================================================================================================
  !

  !>Returns the number of residuals for nonlinear matrices.
  SUBROUTINE EquationsMatricesNonlinear_NumberOfResidualsGet(nonlinearMatrices,numberOfResiduals,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices !<A pointer to the nonlinear matrices to get the number of residuals for
    INTEGER(INTG), INTENT(OUT) :: numberOfResiduals !<On exit, the number of residuals for the nonlinear matrices.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesNonlinear_NumberOfResidualsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(nonlinearMatrices)) CALL FlagError("Nonlinear matrices is not associated.",err,error,*999)
#endif    

    numberOfResiduals=nonlinearMatrices%numberOfResiduals

    EXITS("EquationsMatricesNonlinear_NumberOfResidualsGet")
    RETURN
999 ERRORSEXITS("EquationsMatricesNonlinear_NumberOfResidualsGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesNonlinear_NumberOfResidualsGet

  !
  !================================================================================================================================
  !

  !>Gets a residual vector for nonlinear matrices.
  SUBROUTINE EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices !<A pointer to the nonlinear matrices to get the residual vector for
    INTEGER(INTG), INTENT(IN) :: residualIdx !<The residual number of the residual vector to get
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<On exit, a pointer to the specified residual vector for the nonlinear matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsMatricesNonlinear_ResidualVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(residualVector)) CALL FlagError("Residual vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(nonlinearMatrices)) CALL FlagError("Nonlinear matrices is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(nonlinearMatrices%residuals)) CALL FlagError("Nonlinear matrices residuals is not allocated.",err,error,*999)
    IF(residualIdx<1.OR.residualIdx>nonlinearMatrices%numberOfResiduals) THEN
      localError="The specified residual number of "//TRIM(NumberToVString(residualIdx,"*",err,error))// &
        & " is invalid. The residual number should be >=1 and <= "// &
        & TRIM(NumberToVString(nonlinearMatrices%numberOfResiduals,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    residualVector=>nonlinearMatrices%residuals(residualIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(residualVector)) THEN
      localError="The residual vector for residual number "//TRIM(NumberToVString(residualIdx,"*",err,error))// &
        & " of the nonlinear matrices is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
    EXITS("EquationsMatricesNonlinear_ResidualVectorGet")
    RETURN
999 NULLIFY(residualVector)
998 ERRORS("EquationsMatricesNonlinear_ResidualVectorGet",err,error)
    EXITS("EquationsMatricesNonlinear_ResidualVectorGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesNonlinear_ResidualVectorGet

  !
  !================================================================================================================================
  !

  !>Checks if the temporary distributed vector for nonlinear matrices exists.
  SUBROUTINE EquationsMatricesNonlinear_TempDistributedVectorExists(nonlinearMatrices,tempDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices !<A pointer to the nonlinear matrices to check the temp distributed vector exists for
    TYPE(DistributedVectorType), POINTER :: tempDistributedVector !<On exit, a pointer to the temp distributed vector for the nonlinear matrices if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesNonlinear_TempDistributedVectorExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(tempDistributedVector)) CALL FlagError("Temporary distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(nonlinearMatrices)) CALL FlagError("Nonlinear matrices is not associated.",err,error,*999)
#endif    

    tempDistributedVector=>nonlinearMatrices%tempVector

    EXITS("EquationsMatricesNonlinear_TempDistributedVectorExists")
    RETURN
999 NULLIFY(tempDistributedVector)
998 ERRORSEXITS("EquationsMatricesNonlinear_TempDistributedVectorExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesNonlinear_TempDistributedVectorExists

  !
  !================================================================================================================================
  !

  !>Gets the temporary distributed vector for nonlinear matrices.
  SUBROUTINE EquationsMatricesNonlinear_TempDistributedVectorGet(nonlinearMatrices,tempDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices !<A pointer to the nonlinear matrices to get the temp distributed vector for
    TYPE(DistributedVectorType), POINTER :: tempDistributedVector !<On exit, a pointer to the temp distributed vector for the nonlinear matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesNonlinear_TempDistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(tempDistributedVector)) CALL FlagError("Temporary distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(nonlinearMatrices)) CALL FlagError("Nonlinear matrices is not associated.",err,error,*999)
#endif    

    tempDistributedVector=>nonlinearMatrices%tempVector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(tempDistributedVector)) &
      & CALL FlagError("Temporary distributed vector is not associated for the nonlinear matrices.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesNonlinear_TempDistributedVectorGet")
    RETURN
999 NULLIFY(tempDistributedVector)
998 ERRORSEXITS("EquationsMatricesNonlinear_TempDistributedVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesNonlinear_TempDistributedVectorGet

  !
  !================================================================================================================================
  !

  !>Gets the vector matrices for nonlinear matrices.
  SUBROUTINE EquationsMatricesNonlinear_VectorMatricesGet(nonlinearMatrices,vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices !<A pointer to the nonlinear matrices to get the vector matrices for
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<On exit, a pointer to the vector matrices for the nonlinear matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesNonlinear_VectorMatricesGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(nonlinearMatrices)) CALL FlagError("Nonlinear matrices is not associated.",err,error,*999)
#endif    

    vectorMatrices=>nonlinearMatrices%vectorMatrices

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(vectorMatrices)) &
      & CALL FlagError("Vector matrices is not associated for the nonlinear matrices.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesNonlinear_VectorMatricesGet")
    RETURN
999 NULLIFY(vectorMatrices)
998 ERRORS("EquationsMatricesNonlinear_VectorMatricesGet",err,error)
    EXITS("EquationsMatricesNonlinear_VectorMatricesGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesNonlinear_VectorMatricesGet

  !
  !================================================================================================================================
  !

  !>Gets the Hessian matrix for optimisation matrices.
  SUBROUTINE EquationsMatricesOptimisation_HessianMatrixGet(optimisationMatrices,hessianMatrixIdx,hessianMatrix,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesOptimisationType), POINTER :: optimisationMatrices !<A pointer to the optimisation matrices to get the Hessian matrix for
    INTEGER(INTG), INTENT(IN) :: hessianMatrixIdx !<The Hessian matrix index to get
    TYPE(HessianMatrixType), POINTER :: hessianMatrix !<On exit, a pointer to the Hessian matrix for the optimisation matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsMatricesOptimisation_HessianMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(hessianMatrix)) CALL FlagError("Hessian matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(optimisationMatrices)) CALL FlagError("Optimisation matrices is not associated.",err,error,*999)
    IF(hessianMatrixIdx<1.OR.hessianMatrixIdx>optimisationMatrices%numberOfHessians) THEN
      localError="The specified Hessian matrix index of "//TRIM(NumberToVString(hessianMatrixIdx,"*",err,error))// &
        & " is invalid. The Hessian matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(optimisationMatrices%numberOfHessians,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    IF(.NOT.ALLOCATED(optimisationMatrices%hessians)) &
      & CALL FlagError("The Hessians array is not allocated for the optimisation matrices.",err,error,*999)
#endif    

    hessianMatrix=>optimisationMatrices%hessians(hessianMatrixIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(hessianMatrix)) THEN
      localError="The Hessian matrix is not associated for Hessian matrix index "// &
        & TRIM(NumberToVString(hessianMatrixIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
    EXITS("EquationsMatricesOptimisation_HessianMatrixGet")
    RETURN
999 NULLIFY(hessianMatrix)
998 ERRORS("EquationsMatricesOptimisation_HessianMatrixGet",err,error)
    EXITS("EquationsMatricesOptimisation_HessianMatrixGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesOptimisation_HessianMatrixGet

  !
  !================================================================================================================================
  !

  !>Checks that the temporal type residual distributed vector exists for a equations matrices residual vector
  SUBROUTINE EquationsMatricesResidual_DistributedVectorExists(residualVector,temporalType,residualDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual to get the distributed vector for
    INTEGER(INTG), INTENT(IN) :: temporalType !<The temporal type of the vector to get \see EquationsMatricesRoutines_VectorTemporalTypes
    TYPE(DistributedVectorType), POINTER :: residualDistributedVector !<On exit, a pointer to the residual distributed vector for the equations matrices residual. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsMatricesResidual_DistributedVectorExists",err,error,*998)

    SELECT CASE(temporalType)
    CASE(EQUATIONS_MATRICES_CURRENT_VECTOR)
      CALL EquationsMatricesResidual_CurrentDistributedVectorExists(residualVector,residualDistributedVector,err,error,*998)
    CASE(EQUATIONS_MATRICES_PREVIOUS_VECTOR)
      CALL EquationsMatricesResidual_PreviousDistributedVectorExists(residualVector,residualDistributedVector,err,error,*998)
    CASE(EQUATIONS_MATRICES_PREVIOUS2_VECTOR)
      CALL EquationsMatricesResidual_Previous2DistributedVectorExists(residualVector,residualDistributedVector,err,error,*998)
    CASE(EQUATIONS_MATRICES_PREVIOUS3_VECTOR)
      CALL EquationsMatricesResidual_Previous3DistributedVectorExists(residualVector,residualDistributedVector,err,error,*998)
    CASE DEFAULT
      localError="The specified residual vector temporal type of "//TRIM(NumberToVString(temporalType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("EquationsMatricesResidual_DistributedVectorExists")
    RETURN
999 NULLIFY(residualDistributedVector)
998 ERRORSEXITS("EquationsMatricesResidual_DistributedVectorExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_DistributedVectorExists

  !
  !================================================================================================================================
  !

  !>Gets the residual distributed vector for a equations matrices residual vector
  SUBROUTINE EquationsMatricesResidual_DistributedVectorGet(residualVector,temporalType,residualDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual to get the distributed vector for
    INTEGER(INTG), INTENT(IN) :: temporalType !<The temporal type of the vector to get \see EquationsMatricesRoutines_VectorTemporalTypes
    TYPE(DistributedVectorType), POINTER :: residualDistributedVector !<On exit, a pointer to the residual distributed vector for the equations matrices residual. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsMatricesResidual_DistributedVectorGet",err,error,*998)

    SELECT CASE(temporalType)
    CASE(EQUATIONS_MATRICES_CURRENT_VECTOR)
      CALL EquationsMatricesResidual_CurrentDistributedVectorGet(residualVector,residualDistributedVector,err,error,*998)
    CASE(EQUATIONS_MATRICES_PREVIOUS_VECTOR)
      CALL EquationsMatricesResidual_PreviousDistributedVectorGet(residualVector,residualDistributedVector,err,error,*998)
    CASE(EQUATIONS_MATRICES_PREVIOUS2_VECTOR)
      CALL EquationsMatricesResidual_Previous2DistributedVectorGet(residualVector,residualDistributedVector,err,error,*998)
    CASE(EQUATIONS_MATRICES_PREVIOUS3_VECTOR)
      CALL EquationsMatricesResidual_Previous3DistributedVectorGet(residualVector,residualDistributedVector,err,error,*998)
    CASE DEFAULT
      localError="The specified residual vector temporal type of "//TRIM(NumberToVString(temporalType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("EquationsMatricesResidual_DistributedVectorGet")
    RETURN
999 NULLIFY(residualDistributedVector)
998 ERRORSEXITS("EquationsMatricesResidual_DistributedVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_DistributedVectorGet

  !
  !================================================================================================================================
  !

  !>Outputs the element vector for a equations matrices residual vector
  SUBROUTINE EquationsMatricesResidual_ElementVectorOutput(id,residualVector,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the ouptut stream
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual to output the element vector for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesResidual_ElementVectorOutput",err,error,*999)

    IF(.NOT.ASSOCIATED(residualVector)) CALL FlagError("Residual vector is not associated.",err,error,*999)

    CALL ElementVector_Output(id,residualVector%elementResidual,err,error,*999)
    
    EXITS("EquationsMatricesResidual_ElementVectorOutput")
    RETURN
999 ERRORSEXITS("EquationsMatricesResidual_ElementVectorOutput",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_ElementVectorOutput

  !
  !================================================================================================================================
  !

  !>Returns the first assembly flag for a residual vector.
  SUBROUTINE EquationsMatricesResidual_FirstAssemblyGet(residualVector,firstAssembly,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual vector to get the first assembly flag for
    LOGICAL, INTENT(OUT) :: firstAssembly !<On exit, the update flag for the residual vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesResidual_FirstAssemblyGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(residualVector)) CALL FlagError("Residual vector is not associated.",err,error,*999)
#endif    

    firstAssembly=residualVector%firstAssembly

    EXITS("EquationsMatricesResidual_FirstAssemblyGet")
    RETURN
999 ERRORSEXITS("EquationsMatricesResidual_FirstAssemblyGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_FirstAssemblyGet

  !
  !================================================================================================================================
  !

  !>Checks the current residual distributed vector exists for a equations matrices residual vector 
  SUBROUTINE EquationsMatricesResidual_CurrentDistributedVectorExists(residualVector,currentResidualDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual to check the currrent distributed vector exists
    TYPE(DistributedVectorType), POINTER :: currentResidualDistributedVector !<On exit, a pointer to the current residual distributed vector for the equations matrices residual if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesResidual_CurrentDistributedVectorExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(currentResidualDistributedVector)) &
      & CALL FlagError("The current residual distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(residualVector)) CALL FlagError("Residual vector is not associated.",err,error,*999)
#endif    

    currentResidualDistributedVector=>residualVector%residual
      
    EXITS("EquationsMatricesResidual_CurrentDistributedVectorExists")
    RETURN
999 NULLIFY(currentResidualDistributedVector)
998 ERRORS("EquationsMatricesResidual_CurrentDistributedVectorExists",err,error)
    EXITS("EquationsMatricesResidual_CurrentDistributedVectorExists")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_CurrentDistributedVectorExists

  !
  !================================================================================================================================
  !

  !>Gets the current residual distributed vector for a equations matrices residual vector
  SUBROUTINE EquationsMatricesResidual_CurrentDistributedVectorGet(residualVector,currentResidualDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual to get the currrent distributed vector for
    TYPE(DistributedVectorType), POINTER :: CurrentResidualDistributedVector !<On exit, a pointer to the current residual distributed vector for the equations matrices residual. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesResidual_CurrentDistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(currentResidualDistributedVector)) &
      & CALL FlagError("The current residual distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(residualVector)) CALL FlagError("Residual vector is not associated.",err,error,*999)
#endif    

    currentResidualDistributedVector=>residualVector%residual

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(currentResidualDistributedVector)) &
      & CALL FlagError("The current distributed vector is not associated for the residual vector.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesResidual_CurrentDistributedVectorGet")
    RETURN
999 NULLIFY(currentResidualDistributedVector)
998 ERRORS("EquationsMatricesResidual_CurrentDistributedVectorGet",err,error)
    EXITS("EquationsMatricesResidual_CurrentDistributedVectorGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_CurrentDistributedVectorGet

  !
  !================================================================================================================================
  !

  !>Checks if the previous residual distributed vector exists for a equations matrices residual vector
  SUBROUTINE EquationsMatricesResidual_PreviousDistributedVectorExists(residualVector,previousResidualDistributedVector, &
    & err,error,*)

    !Argument variables
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual to check if the previous distributed vector exists
    TYPE(DistributedVectorType), POINTER :: previousResidualDistributedVector !<On exit, a pointer to the previous residual distributed vector for the equations matrices residual if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesResidual_PreviousDistributedVectorExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(previousResidualDistributedVector)) &
      & CALL FlagError("Previous residual distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(residualVector)) CALL FlagError("Residual vector is not associated.",err,error,*999)
#endif    

    previousResidualDistributedVector=>residualVector%previousResidual
       
    EXITS("EquationsMatricesResidual_PreviousDistributedVectorExists")
    RETURN
999 NULLIFY(previousResidualDistributedVector)
998 ERRORS("EquationsMatricesResidual_PreviousDistributedVectorExists",err,error)
    EXITS("EquationsMatricesResidual_PreviousDistributedVectorExists")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_PreviousDistributedVectorExists

  !
  !================================================================================================================================
  !

  !>Gets the previous residual distributed vector for a equations matrices residual vector
  SUBROUTINE EquationsMatricesResidual_PreviousDistributedVectorGet(residualVector,previousResidualDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual to get the previous distributed vector for
    TYPE(DistributedVectorType), POINTER :: previousResidualDistributedVector !<On exit, a pointer to the previous residual distributed vector for the equations matrices residual. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesResidual_PreviousDistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(previousResidualDistributedVector)) &
      & CALL FlagError("Previous residual distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(residualVector)) CALL FlagError("Residual vector is not associated.",err,error,*999)
#endif    

    previousResidualDistributedVector=>residualVector%previousResidual

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(previousResidualDistributedVector)) &
      & CALL FlagError("The previous distributed vector is not associated for the residual vector.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesResidual_PreviousDistributedVectorGet")
    RETURN
999 NULLIFY(previousResidualDistributedVector)
998 ERRORS("EquationsMatricesResidual_PreviousDistributedVectorGet",err,error)
    EXITS("EquationsMatricesResidual_PreviousDistributedVectorGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_PreviousDistributedVectorGet

  !
  !================================================================================================================================
  !

  !>Checks if the second previous residual distributed vector exists for a equations matrices residual vector
  SUBROUTINE EquationsMatricesResidual_Previous2DistributedVectorExists(residualVector,previous2ResidualDistributedVector, &
    & err,error,*)

    !Argument variables
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual to check if the second previous distributed vector exists
    TYPE(DistributedVectorType), POINTER :: previous2ResidualDistributedVector !<On exit, a pointer to the second previous residual distributed vector for the equations matrices residual if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesResidual_Previous2DistributedVectorExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(previous2ResidualDistributedVector)) &
      & CALL FlagError("Second previous residual distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(residualVector)) CALL FlagError("Residual vector is not associated.",err,error,*999)
#endif    

    previous2ResidualDistributedVector=>residualVector%previous2Residual
       
    EXITS("EquationsMatricesResidual_Previous2DistributedVectorExists")
    RETURN
999 NULLIFY(previous2ResidualDistributedVector)
998 ERRORS("EquationsMatricesResidual_Previous2DistributedVectorExists",err,error)
    EXITS("EquationsMatricesResidual_Previous2DistributedVectorExists")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_Previous2DistributedVectorExists

  !
  !================================================================================================================================
  !

  !>Gets the second previous residual distributed vector for a equations matrices residual vector
  SUBROUTINE EquationsMatricesResidual_Previous2DistributedVectorGet(residualVector,previous2ResidualDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual to get the second previous distributed vector for
    TYPE(DistributedVectorType), POINTER :: previous2ResidualDistributedVector !<On exit, a pointer to the second previous residual distributed vector for the equations matrices residual. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesResidual_Previous2DistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(previous2ResidualDistributedVector)) &
      & CALL FlagError("Second previous residual distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(residualVector)) CALL FlagError("Residual vector is not associated.",err,error,*999)
#endif    

    previous2ResidualDistributedVector=>residualVector%previous2Residual

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(previous2ResidualDistributedVector)) &
      & CALL FlagError("The second previous distributed vector is not associated for the residual vector.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesResidual_Previous2DistributedVectorGet")
    RETURN
999 NULLIFY(previous2ResidualDistributedVector)
998 ERRORS("EquationsMatricesResidual_Previous2DistributedVectorGet",err,error)
    EXITS("EquationsMatricesResidual_Previous2DistributedVectorGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_Previous2DistributedVectorGet

  !
  !================================================================================================================================
  !

  !>Checks if the third previous residual distributed vector exists for a equations matrices residual vector
  SUBROUTINE EquationsMatricesResidual_Previous3DistributedVectorExists(residualVector,previous3ResidualDistributedVector, &
    & err,error,*)

    !Argument variables
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual to check if the third previous distributed vector exists
    TYPE(DistributedVectorType), POINTER :: previous3ResidualDistributedVector !<On exit, a pointer to the third  previous residual distributed vector for the equations matrices residual if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesResidual_Previous3DistributedVectorExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(previous3ResidualDistributedVector)) &
      & CALL FlagError("Third previous residual distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(residualVector)) CALL FlagError("Residual vector is not associated.",err,error,*999)
#endif    

    previous3ResidualDistributedVector=>residualVector%previous3Residual
       
    EXITS("EquationsMatricesResidual_Previous3DistributedVectorExists")
    RETURN
999 NULLIFY(previous3ResidualDistributedVector)
998 ERRORS("EquationsMatricesResidual_Previous3DistributedVectorExists",err,error)
    EXITS("EquationsMatricesResidual_Previous3DistributedVectorExists")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_Previous3DistributedVectorExists

  !
  !================================================================================================================================
  !

  !>Gets the third previous residual distributed vector for a equations matrices residual vector
  SUBROUTINE EquationsMatricesResidual_Previous3DistributedVectorGet(residualVector,previous3ResidualDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual to get the third previous distributed vector for
    TYPE(DistributedVectorType), POINTER :: previous3ResidualDistributedVector !<On exit, a pointer to the third previous residual distributed vector for the equations matrices residual. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesResidual_Previous3DistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(previous3ResidualDistributedVector)) &
      & CALL FlagError("The third previous residual distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(residualVector)) CALL FlagError("Residual vector is not associated.",err,error,*999)
#endif    

    previous3ResidualDistributedVector=>residualVector%previous3Residual

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(previous3ResidualDistributedVector)) &
      & CALL FlagError("The third previous distributed vector is not associated for the residual vector.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesResidual_Previous3DistributedVectorGet")
    RETURN
999 NULLIFY(previous3ResidualDistributedVector)
998 ERRORS("EquationsMatricesResidual_Previous3DistributedVectorGet",err,error)
    EXITS("EquationsMatricesResidual_Previous3DistributedVectorGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_Previous3DistributedVectorGet

  !
  !================================================================================================================================
  !

  !>Gets the specified Jacobian matrix for a residual vector.
  SUBROUTINE EquationsMatricesResidual_JacobianMatrixGet(residualVector,matrixIdx,JacobianMatrix,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual vector to get the Jacobian matrix for
    INTEGER(INTG), INTENT(IN) :: matrixIdx !<The matrix index of the Jacobian matrix to get
    TYPE(JacobianMatrixType), POINTER :: JacobianMatrix !<On exit, a pointer to the Jacobian matrix for the matrixIdx'th Jacobian matrix of the residual. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsMatricesResidual_JacobianMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(JacobianMatrix)) CALL FlagError("Jacobian matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(residualVector)) CALL FlagError("Residual vector is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(residualVector%jacobians)) CALL FlagError("Residual vecotr Jacobians is not allocated.",err,error,*999)
    IF(matrixIdx<1.OR.matrixIdx>residualVector%numberOfJacobians) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is invalid for residual vector number "//TRIM(NumberToVString(residualVector%residualNumber,"*",err,error))// &
        & ". The matrix index must be >= 1 and <= "//TRIM(NumberToVString(residualVector%numberOfJacobians,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    jacobianMatrix=>residualVector%jacobians(matrixIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(jacobianMatrix)) THEN
      localError="Jacobian matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is not associated for residual number "//TRIM(NumberToVString(residualVector%residualNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("EquationsMatricesResidual_JacobianMatrixGet")
    RETURN
999 NULLIFY(jacobianMatrix)
998 ERRORS("EquationsMatricesResidual_JacobianMatrixGet",err,error)
    EXITS("EquationsMatricesResidual_JacobianMatrixGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_JacobianMatrixGet

  !
  !================================================================================================================================
  !

  !>Outputs the nodal vector for a equations matrices residual vector
  SUBROUTINE EquationsMatricesResidual_NodalVectorOutput(id,residualVector,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the ouptut stream
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual to output the nodal vector for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesResidual_NodalVectorOutput",err,error,*999)

    IF(.NOT.ASSOCIATED(residualVector)) CALL FlagError("Residual vector is not associated.",err,error,*999)

    CALL NodalVector_Output(id,residualVector%nodalResidual,err,error,*999)
    
    EXITS("EquationsMatricesResidual_NodalVectorOutput")
    RETURN
999 ERRORSEXITS("EquationsMatricesResidual_NodalVectorOutput",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_NodalVectorOutput

  !
  !================================================================================================================================
  !

  !>Gets the nonlinear matrices for a residual vector.
  SUBROUTINE EquationsMatricesResidual_NonlinearMatricesGet(residualVector,nonlinearMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual vector to get the nonlinear matrices for
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices !<On exit, a pointer to the nonlinear matrices for the residual vector. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesResidual_NonlinearMatricesGet",err,error,*998)

#ifdef WITH_PRECHEKCS    
    IF(ASSOCIATED(nonlinearMatrices)) CALL FlagError("Nonlinear matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(residualVector)) CALL FlagError("Residual vector is not associated.",err,error,*999)
#endif    

    nonlinearMatrices=>residualVector%nonlinearMatrices

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(nonlinearMatrices)) &
      & CALL FlagError("Nonlinear matrices is not associated for the residual vector.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesResidual_NonlinearMatricesGet")
    RETURN
999 NULLIFY(nonlinearMatrices)
998 ERRORS("EquationsMatricesResidual_NonlinearMatricesGet",err,error)
    EXITS("EquationsMatricesResidual_NonlinearMatricesGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_NonlinearMatricesGet

  !
  !================================================================================================================================
  !

  !>Returns the number of Jacobians for a residual vector.
  SUBROUTINE EquationsMatricesResidual_NumberOfJacobiansGet(residualVector,numberOfJacobians,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual vector to get the number of Jacobianas for
    INTEGER(INTG), INTENT(OUT) :: numberOfJacobians !<On exit, the number of Jacobians for the residual vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesResidual_NumberOfJacobiansGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(residualVector)) CALL FlagError("Residual vector is not associated.",err,error,*999)
#endif    

    numberOfJacobians=residualVector%numberOfJacobians

    EXITS("EquationsMatricesResidual_NumberOfJacobiansGet")
    RETURN
999 ERRORSEXITS("EquationsMatricesResidual_NumberOfJacobiansGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_NumberOfJacobiansGet

  !
  !================================================================================================================================
  !

  !>Gets the residual mapping for a residual vector.
  SUBROUTINE EquationsMatricesResidual_ResidualMappingGet(residualVector,residualMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual vector to get the residual mappipng for
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping !<On exit, a pointer to the residual mapping for the residual vector. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesResidual_ResidualMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(residualMapping)) CALL FlagError("Residual mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(residualVector)) CALL FlagError("Residual vector is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(residualVector%nonlinearMatrices)) &
      & CALL FlagError("Nonlinear matrices is not associated for the residual vector.",err,error,*999)
    IF(.NOT.ASSOCIATED(residualVector%nonlinearMatrices%vectorMatrices)) &
      & CALL FlagError("Vector matrices is not associated for the nonlinear matrices.",err,error,*999)
    IF(.NOT.ASSOCIATED(residualVector%nonlinearMatrices%vectorMatrices%vectorEquations)) &
      & CALL FlagError("Vector equations is not associated for the vector matrices.",err,error,*999)
    IF(.NOT.ASSOCIATED(residualVector%nonlinearMatrices%vectorMatrices%vectorEquations%vectorMapping)) &
      & CALL FlagError("Vector mapping is not associated for the vector equations.",err,error,*999)
    IF(.NOT.ASSOCIATED(residualVector%nonlinearMatrices%vectorMatrices%vectorEquations%vectorMapping%nonlinearMapping)) &
      & CALL FlagError("Nonlinear mapping is not associated for the vector mapping.",err,error,*999)
    IF(.NOT.ALLOCATED(residualVector%nonlinearMatrices%vectorMatrices%vectorEquations%vectorMapping%nonlinearMapping%residuals)) &
      & CALL FlagError("Residuals is not allocated for the nonlinear mspping.",err,error,*999)
#endif    

    residualMapping=>residualVector%nonlinearMatrices%vectorMatrices%vectorEquations%vectorMapping%nonlinearMapping% &
      & residuals(residualVector%residualNumber)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(residualMapping)) &
      & CALL FlagError("Residual mapping is not associated for the residual vector.",err,error,*999)
#endif
       
    EXITS("EquationsMatricesResidual_ResidualMappingGet")
    RETURN
999 NULLIFY(residualMapping)
998 ERRORS("EquationsMatricesResidual_ResidualMappingGet",err,error)
    EXITS("EquationsMatricesResidual_ResidualMappingGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_ResidualMappingGet

  !
  !================================================================================================================================
  !

  !>Returns the residual number for a residual vector.
  SUBROUTINE EquationsMatricesResidual_ResidualNumberGet(residualVector,residualNumber,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual vector to get the residual number for
    INTEGER(INTG), INTENT(OUT) :: residualNumber !<On exit, the residual number for the residual vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesResidual_ResidualNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(residualVector)) CALL FlagError("Residual vector is not associated.",err,error,*999)
#endif    

    residualNumber=residualVector%residualNumber

    EXITS("EquationsMatricesResidual_ResidualNumberGet")
    RETURN
999 ERRORSEXITS("EquationsMatricesResidual_ResidualNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_ResidualNumberGet

  !
  !================================================================================================================================
  !

  !>Returns the update flag for a residual vector.
  SUBROUTINE EquationsMatricesResidual_UpdateVectorGet(residualVector,updateVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual vector to get the  update flag for
    LOGICAL, INTENT(OUT) :: updateVector !<On exit, the update flag for the residual vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesResidual_UpdateVectorGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(residualVector)) CALL FlagError("Residual vector is not associated.",err,error,*999)
#endif    

    updateVector=residualVector%updateResidual

    EXITS("EquationsMatricesResidual_UpdateVectorGet")
    RETURN
999 ERRORSEXITS("EquationsMatricesResidual_UpdateVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_UpdateVectorGet

  !
  !================================================================================================================================
  !

  !>Sets/changes the update flag for a residual vector.
  SUBROUTINE EquationsMatricesResidual_UpdateVectorSet(residualVector,updateVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual vector to set the update flag for
    LOGICAL, INTENT(IN) :: updateVector !<The update flag for the residual vector to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesResidual_UpdateVectorSet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(residualVector)) CALL FlagError("Residual vector is not associated.",err,error,*999)
#endif    

    residualVector%updateResidual=updateVector

    EXITS("EquationsMatricesResidual_UpdateVectorSet")
    RETURN
999 ERRORSEXITS("EquationsMatricesResidual_UpdateVectorSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_UpdateVectorSet

  !
  !================================================================================================================================
  !

  !>Returns the vector coefficient for a residual vector.
  SUBROUTINE EquationsMatricesResidual_VectorCoefficientGet(residualVector,vectorCoefficient,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual vector to get the vector coefficient for
    REAL(DP), INTENT(OUT) :: vectorCoefficient !<On exit, the vector coefficient for the residual vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesResidual_VectorCoefficientGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(residualVector)) CALL FlagError("Residual vector is not associated.",err,error,*999)
#endif    

    vectorCoefficient=residualVector%residualCoefficient

    EXITS("EquationsMatricesResidual_VectorCoefficientGet")
    RETURN
999 ERRORSEXITS("EquationsMatricesResidual_VectorCoefficientGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_VectorCoefficientGet

  !
  !================================================================================================================================
  !

  !>Checks that the temporal type rhs distributed vector exists for a equations matrices rhs vector
  SUBROUTINE EquationsMatricesRHS_DistributedVectorExists(rhsVector,temporalType,rhsDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector !<A pointer to the RHS to get the distributed vector for
    INTEGER(INTG), INTENT(IN) :: temporalType !<The temporal type of the vector to get \see EquationsMatricesRoutines_VectorTemporalTypes
    TYPE(DistributedVectorType), POINTER :: rhsDistributedVector !<On exit, a pointer to the RHS distributed vector for the equations matrices RHS. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsMatricesRHS_DistributedVectorExists",err,error,*998)

    SELECT CASE(temporalType)
    CASE(EQUATIONS_MATRICES_CURRENT_VECTOR)
      CALL EquationsMatricesRHS_CurrentDistributedVectorExists(rhsVector,rhsDistributedVector,err,error,*998)
    CASE(EQUATIONS_MATRICES_PREVIOUS_VECTOR)
      CALL EquationsMatricesRHS_PreviousDistributedVectorExists(rhsVector,rhsDistributedVector,err,error,*998)
    CASE(EQUATIONS_MATRICES_PREVIOUS2_VECTOR)
      CALL EquationsMatricesRHS_Previous2DistributedVectorExists(rhsVector,rhsDistributedVector,err,error,*998)
    CASE(EQUATIONS_MATRICES_PREVIOUS3_VECTOR)
      CALL EquationsMatricesRHS_Previous3DistributedVectorExists(rhsVector,rhsDistributedVector,err,error,*998)
    CASE DEFAULT
      localError="The specified RHS vector temporal type of "//TRIM(NumberToVString(temporalType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("EquationsMatricesRHS_DistributedVectorExists")
    RETURN
999 NULLIFY(rhsDistributedVector)
998 ERRORSEXITS("EquationsMatricesRHS_DistributedVectorExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesRHS_DistributedVectorExists

  !
  !================================================================================================================================
  !

  !>Gets the RHS distributed vector for a equations matrices RHS vector
  SUBROUTINE EquationsMatricesRHS_DistributedVectorGet(rhsVector,temporalType,rhsDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector !<A pointer to the RHS to get the distributed vector for
    INTEGER(INTG), INTENT(IN) :: temporalType !<The temporal type of the vector to get \see EquationsMatricesRoutines_VectorTemporalTypes
    TYPE(DistributedVectorType), POINTER :: rhsDistributedVector !<On exit, a pointer to the RHS distributed vector for the equations matrices RHS. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsMatricesRHS_DistributedVectorGet",err,error,*998)

    SELECT CASE(temporalType)
    CASE(EQUATIONS_MATRICES_CURRENT_VECTOR)
      CALL EquationsMatricesRHS_CurrentDistributedVectorGet(rhsVector,rhsDistributedVector,err,error,*998)
    CASE(EQUATIONS_MATRICES_PREVIOUS_VECTOR)
      CALL EquationsMatricesRHS_PreviousDistributedVectorGet(rhsVector,rhsDistributedVector,err,error,*998)
    CASE(EQUATIONS_MATRICES_PREVIOUS2_VECTOR)
      CALL EquationsMatricesRHS_Previous2DistributedVectorGet(rhsVector,rhsDistributedVector,err,error,*998)
    CASE(EQUATIONS_MATRICES_PREVIOUS3_VECTOR)
      CALL EquationsMatricesRHS_Previous3DistributedVectorGet(rhsVector,rhsDistributedVector,err,error,*998)
    CASE DEFAULT
      localError="The specified RHS vector temporal type of "//TRIM(NumberToVString(temporalType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("EquationsMatricesRHS_DistributedVectorGet")
    RETURN
999 NULLIFY(rhsDistributedVector)
998 ERRORSEXITS("EquationsMatricesRHS_DistributedVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesRHS_DistributedVectorGet

  !
  !================================================================================================================================
  !

  !>Outputs the element vector for a equations matrices RHS vector
  SUBROUTINE EquationsMatricesRHS_ElementVectorOutput(id,rhsVector,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the ouptut stream
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector !<A pointer to the RHS to output the element vector for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesRHS_ElementVectorOutput",err,error,*999)

    IF(.NOT.ASSOCIATED(rhsVector)) CALL FlagError("RHS vector is not associated.",err,error,*999)

    CALL ElementVector_Output(id,rhsVector%elementVector,err,error,*999)
    
    EXITS("EquationsMatricesRHS_ElementVectorOutput")
    RETURN
999 ERRORSEXITS("EquationsMatricesRHS_ElementVectorOutput",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesRHS_ElementVectorOutput

  !
  !================================================================================================================================
  !

  !>Returns the first assembly flag for a RHS vector.
  SUBROUTINE EquationsMatricesRHS_FirstAssemblyGet(rhsVector,firstAssembly,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector !<A pointer to the RHS vector to get the first assembly flag for
    LOGICAL, INTENT(OUT) :: firstAssembly !<On exit, the update flag for the RHS vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesRHS_FirstAssemblyGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(rhsVector)) CALL FlagError("RHS vector is not associated.",err,error,*999)
#endif    

    firstAssembly=rhsVector%firstAssembly

    EXITS("EquationsMatricesRHS_FirstAssemblyGet")
    RETURN
999 ERRORSEXITS("EquationsMatricesRHS_FirstAssemblyGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesRHS_FirstAssemblyGet

  !
  !================================================================================================================================
  !

  !>Checks if the current RHS distributed vector exists for a equations matrices RHS
  SUBROUTINE EquationsMatricesRHS_CurrentDistributedVectorExists(equationsMatricesRHS,currentRHSDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesRHSType), POINTER :: equationsMatricesRHS !<A pointer to the equations matrices RHS to check the current RHS distributed vector exists
    TYPE(DistributedVectorType), POINTER :: currentRHSDistributedVector !<On exit, a pointer to the current RHS distributed vector for the equations matrices RHS if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesRHS_CurrentDistributedVectorExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(currentRHSDistributedVector)) &
      & CALL FlagError("The current RHS distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatricesRHS)) CALL FlagError("Equations matrices RHS is not associated.",err,error,*999)
#endif    

    currentRHSDistributedVector=>equationsMatricesRHS%vector
       
    EXITS("EquationsMatricesRHS_CurrentDistributedVectorExists")
    RETURN
999 NULLIFY(currentRHSDistributedVector)
998 ERRORS("EquationsMatricesRHS_CurrentDistributedVectorExists",err,error)
    EXITS("EquationsMatricesRHS_CurrentDistributedVectorExists")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesRHS_CurrentDistributedVectorExists

  !
  !================================================================================================================================
  !

  !>Gets the current RHS distributed vector for a equations matrices RHS
  SUBROUTINE EquationsMatricesRHS_CurrentDistributedVectorGet(equationsMatricesRHS,currentRHSDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesRHSType), POINTER :: equationsMatricesRHS !<A pointer to the equations matrices RHS to get the current RHS distributed vector for
    TYPE(DistributedVectorType), POINTER :: currentRHSDistributedVector !<On exit, a pointer to the current RHS distributed vector for the equations matrices RHS. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesRHS_CurrentDistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(currentRHSDistributedVector)) &
      & CALL FlagError("The current RHS distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatricesRHS)) CALL FlagError("Equations matrices RHS is not associated.",err,error,*999)
#endif    

    currentRHSDistributedVector=>equationsMatricesRHS%vector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(currentRHSDistributedVector)) &
      & CALL FlagError("Equations matrices current RHS distributed vector is not associated.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesRHS_CurrentDistributedVectorGet")
    RETURN
999 NULLIFY(currentRHSDistributedVector)
998 ERRORS("EquationsMatricesRHS_CurrentDistributedVectorGet",err,error)
    EXITS("EquationsMatricesRHS_CurrentDistributedVectorGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesRHS_CurrentDistributedVectorGet

  !
  !================================================================================================================================
  !

  !>Checks if the previous RHS distributed vector for a equations matrices RHS exists.
  SUBROUTINE EquationsMatricesRHS_PreviousDistributedVectorExists(equationsMatricesRHS,previousRHSDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesRHSType), POINTER :: equationsMatricesRHS !<A pointer to the equations matrices RHS to check the previous RHS distributed vector exists
    TYPE(DistributedVectorType), POINTER :: previousRHSDistributedVector !<On exit, a pointer to the previous RHS distributed vector for the equations matrices RHS if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesRHS_PreviousDistributedVectorExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(previousRHSDistributedVector)) &
      & CALL FlagError("THe previous RHS distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatricesRHS)) CALL FlagError("Equations matrices RHS is not associated.",err,error,*999)
#endif    

    previousRHSDistributedVector=>equationsMatricesRHS%previousRHSVector
      
    EXITS("EquationsMatricesRHS_PreviousDistributedVectorExists")
    RETURN
999 NULLIFY(previousRHSDistributedVector)
998 ERRORSEXITS("EquationsMatricesRHS_PreviousDistributedVectorExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesRHS_PreviousDistributedVectorExists

  !
  !================================================================================================================================
  !

  !>Gets the previous RHS distributed vector for a equations matrices RHS
  SUBROUTINE EquationsMatricesRHS_PreviousDistributedVectorGet(equationsMatricesRHS,previousRHSDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesRHSType), POINTER :: equationsMatricesRHS !<A pointer to the equations matrices RHS to get the previous RHS distributed vector for
    TYPE(DistributedVectorType), POINTER :: previousRHSDistributedVector !<On exit, a pointer to the previous RHS distributed vector for the equations matrices RHS. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesRHS_PreviousDistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(previousRHSDistributedVector)) &
      & CALL FlagError("THe previous RHS distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatricesRHS)) CALL FlagError("Equations matrices RHS is not associated.",err,error,*999)
#endif    

    previousRHSDistributedVector=>equationsMatricesRHS%previousRHSVector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(previousRHSDistributedVector)) &
      & CALL FlagError("Equations matrices previous RHS distributed vector is not associated.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesRHS_PreviousDistributedVectorGet")
    RETURN
999 NULLIFY(previousRHSDistributedVector)
998 ERRORS("EquationsMatricesRHS_PreviousDistributedVectorGet",err,error)
    EXITS("EquationsMatricesRHS_PreviousDistributedVectorGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesRHS_PreviousDistributedVectorGet

  !
  !================================================================================================================================
  !

  !>Checks if the second previous RHS distributed vector for a equations matrices RHS exists.
  SUBROUTINE EquationsMatricesRHS_Previous2DistributedVectorExists(equationsMatricesRHS,previous2RHSDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesRHSType), POINTER :: equationsMatricesRHS !<A pointer to the equations matrices RHS to check the second previous RHS distributed vector exists
    TYPE(DistributedVectorType), POINTER :: previous2RHSDistributedVector !<On exit, a pointer to the second previous RHS distributed vector for the equations matrices RHS if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesRHS_Previous2DistributedVectorExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(previous2RHSDistributedVector)) &
      & CALL FlagError("THe second previous RHS distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatricesRHS)) CALL FlagError("Equations matrices RHS is not associated.",err,error,*999)
#endif    

    previous2RHSDistributedVector=>equationsMatricesRHS%previous2RHSVector
      
    EXITS("EquationsMatricesRHS_Previous2DistributedVectorExists")
    RETURN
999 NULLIFY(previous2RHSDistributedVector)
998 ERRORS("EquationsMatricesRHS_Previous2DistributedVectorExists",err,error)
    EXITS("EquationsMatricesRHS_Previous2DistributedVectorExists")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesRHS_Previous2DistributedVectorExists

  !
  !================================================================================================================================
  !

  !>Gets the second previous RHS distributed vector for a equations matrices RHS
  SUBROUTINE EquationsMatricesRHS_Previous2DistributedVectorGet(equationsMatricesRHS,previous2RHSDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesRHSType), POINTER :: equationsMatricesRHS !<A pointer to the equations matrices RHS to get the second previous RHS distributed vector for
    TYPE(DistributedVectorType), POINTER :: previous2RHSDistributedVector !<On exit, a pointer to the second previous RHS distributed vector for the equations matrices RHS. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesRHS_Previous2DistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(previous2RHSDistributedVector)) &
      & CALL FlagError("THe second previous RHS distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatricesRHS)) CALL FlagError("Equations matrices RHS is not associated.",err,error,*999)
#endif    

    previous2RHSDistributedVector=>equationsMatricesRHS%previous2RHSVector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(previous2RHSDistributedVector)) &
      & CALL FlagError("Equations matrices second previous RHS distributed vector is not associated.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesRHS_Previous2DistributedVectorGet")
    RETURN
999 NULLIFY(previous2RHSDistributedVector)
998 ERRORS("EquationsMatricesRHS_Previous2DistributedVectorGet",err,error)
    EXITS("EquationsMatricesRHS_Previous2DistributedVectorGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesRHS_Previous2DistributedVectorGet

  !
  !================================================================================================================================
  !

  !>Checks if the third previous RHS distributed vector for a equations matrices RHS exists.
  SUBROUTINE EquationsMatricesRHS_Previous3DistributedVectorExists(equationsMatricesRHS,previous3RHSDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesRHSType), POINTER :: equationsMatricesRHS !<A pointer to the equations matrices RHS to check the third previous RHS distributed vector exists
    TYPE(DistributedVectorType), POINTER :: previous3RHSDistributedVector !<On exit, a pointer to the third previous RHS distributed vector for the equations matrices RHS if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesRHS_Previous3DistributedVectorExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(previous3RHSDistributedVector)) &
      & CALL FlagError("THe third previous RHS distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatricesRHS)) CALL FlagError("Equations matrices RHS is not associated.",err,error,*999)
#endif    

    previous3RHSDistributedVector=>equationsMatricesRHS%previous3RHSVector
      
    EXITS("EquationsMatricesRHS_Previous3DistributedVectorExists")
    RETURN
999 NULLIFY(previous3RHSDistributedVector)
998 ERRORS("EquationsMatricesRHS_Previous3DistributedVectorExists",err,error)
    EXITS("EquationsMatricesRHS_Previous3DistributedVectorExists")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesRHS_Previous3DistributedVectorExists

  !
  !================================================================================================================================
  !

  !>Gets the third previous RHS distributed vector for a equations matrices RHS
  SUBROUTINE EquationsMatricesRHS_Previous3DistributedVectorGet(equationsMatricesRHS,previous3RHSDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesRHSType), POINTER :: equationsMatricesRHS !<A pointer to the equations matrices RHS to get the third previous RHS distributed vector for
    TYPE(DistributedVectorType), POINTER :: previous3RHSDistributedVector !<On exit, a pointer to the third previous RHS distributed vector for the equations matrices RHS. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesRHS_Previous3DistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(previous3RHSDistributedVector)) &
      & CALL FlagError("THe third previous RHS distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatricesRHS)) CALL FlagError("Equations matrices RHS is not associated.",err,error,*999)
#endif    

    previous3RHSDistributedVector=>equationsMatricesRHS%previous3RHSVector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(previous3RHSDistributedVector)) &
      & CALL FlagError("Equations matrices third previous RHS distributed vector is not associated.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesRHS_Previous3DistributedVectorGet")
    RETURN
999 NULLIFY(previous3RHSDistributedVector)
998 ERRORS("EquationsMatricesRHS_Previous3DistributedVectorGet",err,error)
    EXITS("EquationsMatricesRHS_Previous3DistributedVectorGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesRHS_Previous3DistributedVectorGet

  !
  !================================================================================================================================
  !

  !>Outputs the nodal vector for a equations matrices RHS vector
  SUBROUTINE EquationsMatricesRHS_NodalVectorOutput(id,rhsVector,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the ouptut stream
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector !<A pointer to the RHS to output the nodal vector for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesRHS_NodalVectorOutput",err,error,*999)

    IF(.NOT.ASSOCIATED(rhsVector)) CALL FlagError("RHS vector is not associated.",err,error,*999)

    CALL NodalVector_Output(id,rhsVector%nodalVector,err,error,*999)
    
    EXITS("EquationsMatricesRHS_NodalVectorOutput")
    RETURN
999 ERRORSEXITS("EquationsMatricesRHS_NodalVectorOutput",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesRHS_NodalVectorOutput

  !
  !================================================================================================================================
  !

  !>Returns the vector coefficient for a RHS vector.
  SUBROUTINE EquationsMatricesRHS_VectorCoefficientGet(rhsVector,vectorCoefficient,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector !<A pointer to the RHS vector to get the vector coefficient for
    REAL(DP), INTENT(OUT) :: vectorCoefficient !<On exit, the vector coefficient for the RHS vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesRHS_VectorCoefficientGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(rhsVector)) CALL FlagError("RHS vector is not associated.",err,error,*999)
#endif    

    vectorCoefficient=rhsVector%rhsCoefficient

    EXITS("EquationsMatricesRHS_VectorCoefficientGet")
    RETURN
999 ERRORSEXITS("EquationsMatricesRHS_VectorCoefficientGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesRHS_VectorCoefficientGet

  !
  !================================================================================================================================
  !

  !>Gets the vector matrices for a equations matrices RHS
  SUBROUTINE EquationsMatricesRHS_VectorMatricesGet(equationsMatricesRHS,vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesRHSType), POINTER :: equationsMatricesRHS !<A pointer to the equations matrices RHS to get the vector matrices for
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<On exit, a pointer to the vector matrices for the equations matrices RHS. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesRHS_VectorMatricesGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatricesRHS)) CALL FlagError("Equations matrices RHS is not associated.",err,error,*999)
#endif    

    vectorMatrices=>equationsMatricesRHS%vectorMatrices

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(vectorMatrices)) &
      & CALL FlagError("Equations matrices RHS vector matrices is not associated.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesRHS_VectorMatricesGet")
    RETURN
999 NULLIFY(vectorMatrices)
998 ERRORSEXITS("EquationsMatricesRHS_VectorMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesRHS_VectorMatricesGet

  !
  !================================================================================================================================
  !

  !>Returns the update flag for a RHS vector.
  SUBROUTINE EquationsMatricesRHS_UpdateVectorGet(rhsVector,updateVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector !<A pointer to the RHS vector to get the update flag for
    LOGICAL, INTENT(OUT) :: updateVector !<On exit, the update flag for the RHS vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesRHS_UpdateVectorGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(rhsVector)) CALL FlagError("RHS vector is not associated.",err,error,*999)
#endif    

    updateVector=rhsVector%updateVector

    EXITS("EquationsMatricesRHS_UpdateVectorGet")
    RETURN
999 ERRORSEXITS("EquationsMatricesRHS_UpdateVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesRHS_UpdateVectorGet

  !
  !================================================================================================================================
  !

  !>Sets/changes the update flag for a RHS vector.
  SUBROUTINE EquationsMatricesRHS_UpdateVectorSet(rhsVector,updateVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector !<A pointer to the RHS vector to set the update flag for
    LOGICAL, INTENT(IN) :: updateVector !<The update flag for the RHS vector to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesRHS_UpdateVectorSet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(rhsVector)) CALL FlagError("RHS vector is not associated.",err,error,*999)
#endif    

    rhsVector%updateVector=updateVector

    EXITS("EquationsMatricesRHS_UpdateVectorSet")
    RETURN
999 ERRORSEXITS("EquationsMatricesRHS_UpdateVectorSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesRHS_UpdateVectorSet

  !
  !================================================================================================================================
  !

  !>Gets the scalar equations for an scalar equations matrices.
  SUBROUTINE EquationsMatricesScalar_EquationsScalarGet(scalarMatrices,scalarEquations,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesScalarType), POINTER :: scalarMatrices !<A pointer to the equations scalar matrices to get the scalar equations for
    TYPE(EquationsScalarType), POINTER :: scalarEquations !<On exit, a pointer to the scalar equations in the specified scalar equations matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesScalar_EquationsScalarGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(scalarEquations)) CALL FlagError("Scalar equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(scalarMatrices)) CALL FlagError("Scalar matrices is not associated.",err,error,*999)
#endif    

    scalarEquations=>scalarMatrices%scalarEquations

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(scalarEquations)) &
      & CALL FlagError("Scalar equations is not associated for the scalar matrices.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesScalar_EquationsScalarGet")
    RETURN
999 NULLIFY(scalarEquations)
998 ERRORSEXITS("EquationsMatricesScalar_EquationsScalarGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesScalar_EquationsScalarGet

  !
  !================================================================================================================================
  !

  !>Checks that the temporal type source distributed vector exists for a equations matrices source vector
  SUBROUTINE EquationsMatricesSource_DistributedVectorExists(sourceVector,temporalType,sourceDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<A pointer to the source to get the distributed vector for
    INTEGER(INTG), INTENT(IN) :: temporalType !<The temporal type of the vector to get \see EquationsMatricesRoutines_VectorTemporalTypes
    TYPE(DistributedVectorType), POINTER :: sourceDistributedVector !<On exit, a pointer to the source distributed vector for the equations matrices source. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsMatricesSource_DistributedVectorExists",err,error,*998)

    SELECT CASE(temporalType)
    CASE(EQUATIONS_MATRICES_CURRENT_VECTOR)
      CALL EquationsMatricesSource_CurrentDistributedVectorExists(sourceVector,sourceDistributedVector,err,error,*998)
    CASE(EQUATIONS_MATRICES_PREVIOUS_VECTOR)
      CALL EquationsMatricesSource_PreviousDistributedVectorExists(sourceVector,sourceDistributedVector,err,error,*998)
    CASE(EQUATIONS_MATRICES_PREVIOUS2_VECTOR)
      CALL EquationsMatricesSource_Previous2DistributedVectorExists(sourceVector,sourceDistributedVector,err,error,*998)
    CASE(EQUATIONS_MATRICES_PREVIOUS3_VECTOR)
      CALL EquationsMatricesSource_Previous3DistributedVectorExists(sourceVector,sourceDistributedVector,err,error,*998)
    CASE DEFAULT
      localError="The specified source vector temporal type of "//TRIM(NumberToVString(temporalType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("EquationsMatricesSource_DistributedVectorExists")
    RETURN
999 NULLIFY(sourceDistributedVector)
998 ERRORSEXITS("EquationsMatricesSource_DistributedVectorExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSource_DistributedVectorExists

  !
  !================================================================================================================================
  !

  !>Gets the source distributed vector for a equations matrices source vector
  SUBROUTINE EquationsMatricesSource_DistributedVectorGet(sourceVector,temporalType,sourceDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<A pointer to the source to get the distributed vector for
    INTEGER(INTG), INTENT(IN) :: temporalType !<The temporal type of the vector to get \see EquationsMatricesRoutines_VectorTemporalTypes
    TYPE(DistributedVectorType), POINTER :: sourceDistributedVector !<On exit, a pointer to the source distributed vector for the equations matrices source. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsMatricesSource_DistributedVectorGet",err,error,*998)

    SELECT CASE(temporalType)
    CASE(EQUATIONS_MATRICES_CURRENT_VECTOR)
      CALL EquationsMatricesSource_CurrentDistributedVectorGet(sourceVector,sourceDistributedVector,err,error,*998)
    CASE(EQUATIONS_MATRICES_PREVIOUS_VECTOR)
      CALL EquationsMatricesSource_PreviousDistributedVectorGet(sourceVector,sourceDistributedVector,err,error,*998)
    CASE(EQUATIONS_MATRICES_PREVIOUS2_VECTOR)
      CALL EquationsMatricesSource_Previous2DistributedVectorGet(sourceVector,sourceDistributedVector,err,error,*998)
    CASE(EQUATIONS_MATRICES_PREVIOUS3_VECTOR)
      CALL EquationsMatricesSource_Previous3DistributedVectorGet(sourceVector,sourceDistributedVector,err,error,*998)
    CASE DEFAULT
      localError="The specified source vector temporal type of "//TRIM(NumberToVString(temporalType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("EquationsMatricesSource_DistributedVectorGet")
    RETURN
999 NULLIFY(sourceDistributedVector)
998 ERRORSEXITS("EquationsMatricesSource_DistributedVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSource_DistributedVectorGet

  !
  !================================================================================================================================
  !

  !>Outputs the element vector for a equations matrices source vector
  SUBROUTINE EquationsMatricesSource_ElementVectorOutput(id,sourceVector,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the ouptut stream
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<A pointer to the source to output the element vector for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesSource_ElementVectorOutput",err,error,*999)

    IF(.NOT.ASSOCIATED(sourceVector)) CALL FlagError("Source vector is not associated.",err,error,*999)

    CALL ElementVector_Output(id,sourceVector%elementVector,err,error,*999)
    
    EXITS("EquationsMatricesSource_ElementVectorOutput")
    RETURN
999 ERRORSEXITS("EquationsMatricesSource_ElementVectorOutput",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSource_ElementVectorOutput

  !
  !================================================================================================================================
  !

  !>Returns the first assembly flag for a source vector.
  SUBROUTINE EquationsMatricesSource_FirstAssemblyGet(sourceVector,firstAssembly,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<A pointer to the source vector to get the first assembly flag for
    LOGICAL, INTENT(OUT) :: firstAssembly !<On exit, the update flag for the source vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesSource_FirstAssemblyGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(sourceVector)) CALL FlagError("Source vector is not associated.",err,error,*999)
#endif    

    firstAssembly=sourceVector%firstAssembly

    EXITS("EquationsMatricesSource_FirstAssemblyGet")
    RETURN
999 ERRORSEXITS("EquationsMatricesSource_FirstAssemblyGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSource_FirstAssemblyGet

  !
  !================================================================================================================================
  !

  !>Checks the current source distributed vector exists for a equations matrices source vector
  SUBROUTINE EquationsMatricesSource_CurrentDistributedVectorExists(sourceVector,currentSourceDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<A pointer to the equations matrices source vector to check the current source distributed vector for
    TYPE(DistributedVectorType), POINTER :: currentSourceDistributedVector !<On exit, a pointer to the current source distributed vector for the equations matrices source if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesSource_CurrentDistributedVectorExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(currentSourceDistributedVector)) &
      & CALL FlagError("The current source distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(sourceVector)) CALL FlagError("Source vector is not associated.",err,error,*999)
#endif    

    currentSourceDistributedVector=>sourceVector%vector

    EXITS("EquationsMatricesSource_CurrentDistributedVectorExists")
    RETURN
999 NULLIFY(currentSourceDistributedVector)
998 ERRORS("EquationsMatricesSource_CurrentDistributedVectorExists",err,error)
    EXITS("EquationsMatricesSource_CurrentDistributedVectorExists")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSource_CurrentDistributedVectorExists

  !
  !================================================================================================================================
  !

  !>Gets the current source distributed vector for a equations matrices source vector
  SUBROUTINE EquationsMatricesSource_CurrentDistributedVectorGet(sourceVector,currentSourceDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<A pointer to the equations matrices source vector to get the current source distributed vector for
    TYPE(DistributedVectorType), POINTER :: currentSourceDistributedVector !<On exit, a pointer to the current source distributed vector for the equations matrices SOURCE. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesSource_CurrentDistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(currentSourceDistributedVector)) &
      & CALL FlagError("The current source distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(sourceVector)) CALL FlagError("Source vector is not associated.",err,error,*999)
#endif    

    currentSourceDistributedVector=>sourceVector%vector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(currentSourceDistributedVector)) &
      & CALL FlagError("The current source vector distributed vector is not associated.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesSource_CurrentDistributedVectorGet")
    RETURN
999 NULLIFY(currentSourceDistributedVector)
998 ERRORS("EquationsMatricesSource_CurrentDistributedVectorGet",err,error)
    EXITS("EquationsMatricesSource_CurrentDistributedVectorGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSource_CurrentDistributedVectorGet

  !
  !================================================================================================================================
  !

  !>Checks if the previous source distributed vector for a equations matrices source vector exists.
  SUBROUTINE EquationsMatricesSource_PreviousDistributedVectorExists(sourceVector,previousSourceDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<A pointer to the equations matrices source vector to check the previous source distributed vector exists for
    TYPE(DistributedVectorType), POINTER :: previousSourceDistributedVector !<On exit, a pointer to the previous source distributed vector for the equations matrices source if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesSource_PreviousDistributedVectorExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(previousSourceDistributedVector)) &
      & CALL FlagError("Previous source distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(sourceVector)) CALL FlagError("Source vector is not associated.",err,error,*999)
#endif    

    previousSourceDistributedVector=>sourceVector%previousSourceVector
       
    EXITS("EquationsMatricesSource_PreviousDistributedVectorExists")
    RETURN
999 NULLIFY(previousSourceDistributedVector)
998 ERRORS("EquationsMatricesSource_PreviousDistributedVectorExists",err,error)
    EXITS("EquationsMatricesSource_PreviousDistributedVectorExists")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSource_PreviousDistributedVectorExists

  !
  !================================================================================================================================
  !

  !>Gets the previous source distributed vector for a equations matrices source vector
  SUBROUTINE EquationsMatricesSource_PreviousDistributedVectorGet(sourceVector,previousSourceDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<A pointer to the equations matrices source vector to get the previous source distributed vector for
    TYPE(DistributedVectorType), POINTER :: previousSourceDistributedVector !<On exit, a pointer to the previous source distributed vector for the equations matrices SOURCE. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesSource_PreviousDistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(previousSourceDistributedVector)) &
      & CALL FlagError("Previous source distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(sourceVector)) CALL FlagError("Source vector is not associated.",err,error,*999)
#endif    

    previousSourceDistributedVector=>sourceVector%previousSourceVector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(previousSourceDistributedVector)) &
      & CALL FlagError("The previous source vector distributed vector is not associated.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesSource_PreviousDistributedVectorGet")
    RETURN
999 NULLIFY(previousSourceDistributedVector)
998 ERRORS("EquationsMatricesSource_PreviousDistributedVectorGet",err,error)
    EXITS("EquationsMatricesSource_PreviousDistributedVectorGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSource_PreviousDistributedVectorGet

  !
  !================================================================================================================================
  !

  !>Checks if the second previous source distributed vector for a equations matrices source vector exists.
  SUBROUTINE EquationsMatricesSource_Previous2DistributedVectorExists(sourceVector,previous2SourceDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<A pointer to the equations matrices source vector to check the second previous source distributed vector exists for
    TYPE(DistributedVectorType), POINTER :: previous2SourceDistributedVector !<On exit, a pointer to the second previous source distributed vector for the equations matrices source if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesSource_Previous2DistributedVectorExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(previous2SourceDistributedVector)) &
      & CALL FlagError("The second previous source distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(sourceVector)) CALL FlagError("Source vector is not associated.",err,error,*999)
#endif    

    previous2SourceDistributedVector=>sourceVector%previous2SourceVector
       
    EXITS("EquationsMatricesSource_Previous2DistributedVectorExists")
    RETURN
999 NULLIFY(previous2SourceDistributedVector)
998 ERRORS("EquationsMatricesSource_Previous2DistributedVectorExists",err,error)
    EXITS("EquationsMatricesSource_Previous2DistributedVectorExists")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSource_Previous2DistributedVectorExists

  !
  !================================================================================================================================
  !

  !>Gets the second previous source distributed vector for a equations matrices source vector
  SUBROUTINE EquationsMatricesSource_Previous2DistributedVectorGet(sourceVector,previous2SourceDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<A pointer to the equations matrices source vector to get the second previous source distributed vector for
    TYPE(DistributedVectorType), POINTER :: previous2SourceDistributedVector !<On exit, a pointer to the second previous source distributed vector for the equations matrices SOURCE. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesSource_Previous2DistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(previous2SourceDistributedVector)) &
      & CALL FlagError("The second previous source distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(sourceVector)) CALL FlagError("Source vector is not associated.",err,error,*999)
#endif    

    previous2SourceDistributedVector=>sourceVector%previous2SourceVector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(previous2SourceDistributedVector)) &
      & CALL FlagError("The second previous source vector distributed vector is not associated.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesSource_Previous2DistributedVectorGet")
    RETURN
999 NULLIFY(previous2SourceDistributedVector)
998 ERRORS("EquationsMatricesSource_Previous2DistributedVectorGet",err,error)
    EXITS("EquationsMatricesSource_Previous2DistributedVectorGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSource_Previous2DistributedVectorGet

  !
  !================================================================================================================================
  !

  !>Checks if the third previous source distributed vector for a equations matrices source vector exists.
  SUBROUTINE EquationsMatricesSource_Previous3DistributedVectorExists(sourceVector,previous3SourceDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<A pointer to the equations matrices source vector to check the third previous source distributed vector exists for
    TYPE(DistributedVectorType), POINTER :: previous3SourceDistributedVector !<On exit, a pointer to the third previous source distributed vector for the equations matrices source if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesSource_Previous3DistributedVectorExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(previous3SourceDistributedVector)) &
      & CALL FlagError("The third previous source distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(sourceVector)) CALL FlagError("Source vector is not associated.",err,error,*999)
#endif    

    previous3SourceDistributedVector=>sourceVector%previous3SourceVector
       
    EXITS("EquationsMatricesSource_Previous3DistributedVectorExists")
    RETURN
999 NULLIFY(previous3SourceDistributedVector)
998 ERRORS("EquationsMatricesSource_Previous3DistributedVectorExists",err,error)
    EXITS("EquationsMatricesSource_Previous3DistributedVectorExists")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSource_Previous3DistributedVectorExists

  !
  !================================================================================================================================
  !

  !>Gets the third previous source distributed vector for a equations matrices source vector
  SUBROUTINE EquationsMatricesSource_Previous3DistributedVectorGet(sourceVector,previous3SourceDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<A pointer to the equations matrices source vector to get the third previous source distributed vector for
    TYPE(DistributedVectorType), POINTER :: previous3SourceDistributedVector !<On exit, a pointer to the third previous source distributed vector for the equations matrices SOURCE. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesSource_Previous3DistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(previous3SourceDistributedVector)) &
      & CALL FlagError("The third previous source distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(sourceVector)) CALL FlagError("Source vector is not associated.",err,error,*999)
#endif    

    previous3SourceDistributedVector=>sourceVector%previous3SourceVector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(previous3SourceDistributedVector)) &
      & CALL FlagError("The third previous source vector distributed vector is not associated.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesSource_Previous3DistributedVectorGet")
    RETURN
999 NULLIFY(previous3SourceDistributedVector)
998 ERRORS("EquationsMatricesSource_Previous3DistributedVectorGet",err,error)
    EXITS("EquationsMatricesSource_Previous3DistributedVectorGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSource_Previous3DistributedVectorGet

  !
  !================================================================================================================================
  !

  !>Outputs the nodal vector for a equations matrices source vector
  SUBROUTINE EquationsMatricesSource_NodalVectorOutput(id,sourceVector,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the ouptut stream
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<A pointer to the source to output the element vector for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesSource_NodalVectorOutput",err,error,*999)

    IF(.NOT.ASSOCIATED(sourceVector)) CALL FlagError("Source vector is not associated.",err,error,*999)

    CALL NodalVector_Output(id,sourceVector%nodalVector,err,error,*999)
    
    EXITS("EquationsMatricesSource_NodalVectorOutput")
    RETURN
999 ERRORSEXITS("EquationsMatricesSource_NodalVectorOutput",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSource_NodalVectorOutput

  !
  !================================================================================================================================
  !

  !>Returns the source number for a source vector.
  SUBROUTINE EquationsMatricesSource_SourceNumberGet(sourceVector,sourceNumber,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<A pointer to the source vector to get the source number for
    INTEGER(INTG), INTENT(OUT) :: sourceNumber !<On exit, the source number for the source vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesSource_SourceNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(sourceVector)) CALL FlagError("Source vector is not associated.",err,error,*999)
#endif    

    sourceNumber=sourceVector%sourceNumber

    EXITS("EquationsMatricesSource_SourceNumberGet")
    RETURN
999 ERRORSEXITS("EquationsMatricesSource_SourceNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSource_SourceNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the source vectors for a equations matrices source vector
  SUBROUTINE EquationsMatricesSource_SourceVectorsGet(sourceVector,sourceVectors,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<A pointer to the equations matrices source to get the vector matrices for
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors !<On exit, a pointer to the source vectors for the equations matrices source vector. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesSource_SourceVectorsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(sourceVectors)) CALL FlagError("Source vectors is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(sourceVector)) CALL FlagError("Source vector is not associated.",err,error,*999)
#endif    

    sourceVectors=>sourceVector%sources

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(sourceVectors)) &
      & CALL FlagError("Source vector source vectors is not associated.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesSource_SourceVectorsGet")
    RETURN
999 NULLIFY(sourceVectors)
998 ERRORSEXITS("EquationsMatricesSource_SourceVectorsGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSource_SourceVectorsGet

  !
  !================================================================================================================================
  !

  !>Returns the vector coefficient for a source vector.
  SUBROUTINE EquationsMatricesSource_VectorCoefficientGet(sourceVector,vectorCoefficient,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<A pointer to the source vector to get the vector coefficient for
    REAL(DP), INTENT(OUT) :: vectorCoefficient !<On exit, the vector coefficient for the source vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesSource_VectorCoefficientGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(sourceVector)) CALL FlagError("Source vector is not associated.",err,error,*999)
#endif    

    vectorCoefficient=sourceVector%sourceCoefficient

    EXITS("EquationsMatricesSource_VectorCoefficientGet")
    RETURN
999 ERRORSEXITS("EquationsMatricesSource_VectorCoefficientGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSource_VectorCoefficientGet

  !
  !================================================================================================================================
  !

  !>Returns the update flag for a source vector.
  SUBROUTINE EquationsMatricesSource_UpdateVectorGet(sourceVector,updateVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<A pointer to the source vector to get the  update flag for
    LOGICAL, INTENT(OUT) :: updateVector !<On exit, the vector update flag for the source vector.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesSource_UpdateVectorGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(sourceVector)) CALL FlagError("Source vector is not associated.",err,error,*999)
#endif    

    updateVector=sourceVector%updateVector

    EXITS("EquationsMatricesSource_UpdateVectorGet")
    RETURN
999 ERRORSEXITS("EquationsMatricesSource_UpdateVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSource_UpdateVectorGet

  !
  !================================================================================================================================
  !

  !>Sets/changes the update flag for a source vector.
  SUBROUTINE EquationsMatricesSource_UpdateVectorSet(sourceVector,updateVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<A pointer to the source vector to set the update flag for
    LOGICAL, INTENT(IN) :: updateVector !<The vector update flag for the source vector to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesSource_UpdateVectorSet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(sourceVector)) CALL FlagError("Source vector is not associated.",err,error,*999)
#endif    

    sourceVector%updateVector=updateVector

    EXITS("EquationsMatricesSource_UpdateVectorSet")
    RETURN
999 ERRORSEXITS("EquationsMatricesSource_UpdateVectorSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSource_UpdateVectorSet

  !
  !================================================================================================================================
  !

  !>Returns the number of sources for a equations matrices sources.
  SUBROUTINE EquationsMatricesSources_NumberOfSourcesGet(sourceVectors,numberOfSources,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors !<A pointer to the source vectors to get the number of sources for
    INTEGER(INTG), INTENT(OUT) :: numberOfSources !<On exit, the number of sources for the source vectors.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesSources_NumberOfSourcesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(sourceVectors)) CALL FlagError("Source vectors is not associated.",err,error,*999)
#endif    

    numberOfSources=sourceVectors%numberOfSources

    EXITS("EquationsMatricesSources_NumberOfSourcesGet")
    RETURN
999 ERRORSEXITS("EquationsMatricesSources_NumberOfSourcesGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSources_NumberOfSourcesGet

  !
  !================================================================================================================================
  !

  !>Gets a source vector from the source vectors.
  SUBROUTINE EquationsMatricesSources_SourceVectorGet(sourceVectors,sourceIdx,sourceVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors !<A pointer to the source vectors to get the source vector for
    INTEGER(INTG), INTENT(IN) :: sourceIdx !<The source number of the source vector to get
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<On exit, a pointer to the specified source vector for the soure vectors. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsMatricesSources_SourceVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(sourceVector)) CALL FlagError("Source vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(sourceVectors)) CALL FlagError("Source vectors is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(sourceVectors%sources)) CALL FlagError("Source vectors sources is not allocated.",err,error,*999)
    IF(sourceIdx<1.OR.sourceIdx>sourceVectors%numberOfSources) THEN
      localError="The specified source index of "//TRIM(NumberToVString(sourceIdx,"*",err,error))// &
        & " is invalid. The source index should be >=1 and <= "// &
        & TRIM(NumberToVString(sourceVectors%numberOfSources,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    sourceVector=>sourceVectors%sources(sourceIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(sourceVector)) THEN
      localError="The source vector for source index "//TRIM(NumberToVString(sourceIdx,"*",err,error))// &
        & " of the source vectors is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
    EXITS("EquationsMatricesSources_SourceVectorGet")
    RETURN
999 NULLIFY(sourceVector)
998 ERRORS("EquationsMatricesSources_SourceVectorGet",err,error)
    EXITS("EquationsMatricesSources_SourceVectorGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSources_SourceVectorGet

  !
  !================================================================================================================================
  !

  !>Checks if the temporary distributed vector for source vectors exists.
  SUBROUTINE EquationsMatricesSources_TempDistributedVectorExists(sourceVectors,tempDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors !<A pointer to the source vectors to check the temp distributed vector exists for
    TYPE(DistributedVectorType), POINTER :: tempDistributedVector !<On exit, a pointer to the temp distributed vector for the source vectors if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesSources_TempDistributedVectorExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(tempDistributedVector)) CALL FlagError("Temporary distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(sourceVectors)) CALL FlagError("Source vectors is not associated.",err,error,*999)
#endif    

    tempDistributedVector=>sourceVectors%tempVector

    EXITS("EquationsMatricesSources_TempDistributedVectorExists")
    RETURN
999 NULLIFY(tempDistributedVector)
998 ERRORSEXITS("EquationsMatricesSources_TempDistributedVectorExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSources_TempDistributedVectorExists

  !
  !================================================================================================================================
  !

  !>Gets the temporary distributed vector for source vectors.
  SUBROUTINE EquationsMatricesSources_TempDistributedVectorGet(sourceVectors,tempDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors !<A pointer to the source vectors to get the temp distributed vector for
    TYPE(DistributedVectorType), POINTER :: tempDistributedVector !<On exit, a pointer to the temp distributed vector for the source vectors. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesSources_TempDistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(tempDistributedVector)) CALL FlagError("Temporary distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(sourceVectors)) CALL FlagError("Source vectors is not associated.",err,error,*999)
#endif    

    tempDistributedVector=>sourceVectors%tempVector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(tempDistributedVector)) &
      & CALL FlagError("Temporary distributed vector is not associated for the source vectors.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesSources_TempDistributedVectorGet")
    RETURN
999 NULLIFY(tempDistributedVector)
998 ERRORSEXITS("EquationsMatricesSources_TempDistributedVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSources_TempDistributedVectorGet

  !
  !================================================================================================================================
  !

  !>Gets the vector matrices for a equations matrices source vectors
  SUBROUTINE EquationsMatricesSources_VectorMatricesGet(sourceVectors,vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors !<A pointer to the equations matrices source vectors to get the vector matrices for
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<On exit, a pointer to the vector matrices for the equations matrices source vectors. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesSources_VectorMatricesGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(sourceVectors)) CALL FlagError("Source vectors is not associated.",err,error,*999)
#endif    

    vectorMatrices=>sourceVectors%vectorMatrices

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(vectorMatrices)) &
      & CALL FlagError("The source vectors vector matrices is not associated.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesSources_VectorMatricesGet")
    RETURN
999 NULLIFY(vectorMatrices)
998 ERRORSEXITS("EquationsMatricesSources_VectorMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSources_VectorMatricesGet

  !
  !=================================================================================================================================
  !

  !>Assert that a vector equations matrices has been finished
  SUBROUTINE EquationsMatricesVector_AssertIsFinished(vectorMatrices,err,error,*)

    !Argument Variables
    TYPE(EquationsMatricesVectorType), POINTER, INTENT(IN) :: vectorMatrices !<The vector matrices to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)
#endif    

    IF(.NOT.vectorMatrices%vectorMatricesFinished) &
      & CALL FlagError("Vector equations matrices has not been finished.",err,error,*999)
    
    EXITS("EquationsMatricesVector_AssertIsFinished")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a vector equations matrices has not been finished
  SUBROUTINE EquationsMatricesVector_AssertNotFinished(vectorMatrices,err,error,*)

    !Argument Variables
    TYPE(EquationsMatricesVectorType), POINTER, INTENT(IN) :: vectorMatrices !<The vector matrices to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)
#endif    

    IF(vectorMatrices%vectorMatricesFinished) &
      & CALL FlagError("Vector equations matrices has already been finished.",err,error,*999)
    
    EXITS("EquationsMatricesVector_AssertNotFinished")
    RETURN
999 ERRORSEXITS("EquationsMatricesVector_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_AssertNotFinished
  
  !
  !================================================================================================================================
  !

  !>Checks that the dynamic vector matrices exist for the vector matrices.
  SUBROUTINE EquationsMatricesVector_DynamicMatricesExists(vectorMatrices,dynamicMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations vector matrices to check the dynamic matrices exist for
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices !<On exit, a pointer to the dynamic Matrices in the specified vector equations matrices if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_DynamicMatricesExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dynamicMatrices)) CALL FlagError("Dynamic matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)
#endif    

    dynamicMatrices=>vectorMatrices%dynamicMatrices
       
    EXITS("EquationsMatricesVector_DynamicMatricesExists")
    RETURN
999 NULLIFY(dynamicMatrices)
998 ERRORSEXITS("EquationsMatricesVector_DynamicMatricesExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_DynamicMatricesExists

  !
  !================================================================================================================================
  !

  !>Gets the dynamic vector matrices for the vector matrices.
  SUBROUTINE EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations vector matrices to get the dynamic matrices for
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices !<On exit, a pointer to the dynamic Matrices in the specified vector equations matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_DynamicMatricesGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dynamicMatrices)) CALL FlagError("Dynamic matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)
#endif    

    dynamicMatrices=>vectorMatrices%dynamicMatrices

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dynamicMatrices)) CALL FlagError("Dynamic matrices is not associated for the vector matrices.", &
      & err,error,*999)
#endif    
       
    EXITS("EquationsMatricesVector_DynamicMatricesGet")
    RETURN
999 NULLIFY(dynamicMatrices)
998 ERRORSEXITS("EquationsMatricesVector_DynamicMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_DynamicMatricesGet

  !
  !================================================================================================================================
  !

  !>Checks that the linear vector matrices exists for the vector matrices.
  SUBROUTINE EquationsMatricesVector_LinearMatricesExists(vectorMatrices,linearMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations vector matrices to check the linear matrices exists for
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices !<On exit, a pointer to the linear Matrices in the specified vector equations matrices if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_LinearMatricesExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linearMatrices)) CALL FlagError("Linear matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)
#endif    

    linearMatrices=>vectorMatrices%linearMatrices

    EXITS("EquationsMatricesVector_LinearMatricesExists")
    RETURN
999 NULLIFY(linearMatrices)
998 ERRORSEXITS("EquationsMatricesVector_LinearMatricesExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_LinearMatricesExists

  !
  !================================================================================================================================
  !

  !>Gets the linear vector matrices for the vector matrices.
  SUBROUTINE EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations vector matrices to get the linear matrices for
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices !<On exit, a pointer to the linear Matrices in the specified vector equations matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_LinearMatricesGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linearMatrices)) CALL FlagError("Linear matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)
#endif    

    linearMatrices=>vectorMatrices%linearMatrices

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(linearMatrices)) CALL FlagError("Linear matrices is not associated for the vector matrices.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesVector_LinearMatricesGet")
    RETURN
999 NULLIFY(linearMatrices)
998 ERRORSEXITS("EquationsMatricesVector_LinearMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_LinearMatricesGet

  !
  !================================================================================================================================
  !

  !>Checks that the nonlinear vector matrices exists for the vector matrices.
  SUBROUTINE EquationsMatricesVector_NonlinearMatricesExists(vectorMatrices,nonlinearMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations vector matrices to checks the nonlinear matrices exists for
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices !<On exit, a pointer to the nonlinear Matrices in the specified vector equations matrices if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_NonlinearMatricesExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(nonlinearMatrices)) CALL FlagError("Nonlinear matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)
#endif    

    nonlinearMatrices=>vectorMatrices%nonlinearMatrices
      
    EXITS("EquationsMatricesVector_NonlinearMatricesExists")
    RETURN
999 NULLIFY(nonlinearMatrices)
998 ERRORS("EquationsMatricesVector_NonlinearMatricesExists",err,error)
    EXITS("EquationsMatricesVector_NonlinearMatricesExists")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_NonlinearMatricesExists

  !
  !================================================================================================================================
  !

  !>Gets the nonlinear vector matrices for the vector matrices.
  SUBROUTINE EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations vector matrices to get the nonlinear matrices for
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices !<On exit, a pointer to the nonlinear Matrices in the specified vector equations matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_NonlinearMatricesGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(nonlinearMatrices)) CALL FlagError("Nonlinear matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)
#endif    

    nonlinearMatrices=>vectorMatrices%nonlinearMatrices

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(nonlinearMatrices))  &
      & CALL FlagError("Nonlinear matrices is not associated for the vector matrices.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesVector_NonlinearMatricesGet")
    RETURN
999 NULLIFY(nonlinearMatrices)
998 ERRORS("EquationsMatricesVector_NonlinearMatricesGet",err,error)
    EXITS("EquationsMatricesVector_NonlinearMatricesGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_NonlinearMatricesGet

  !
  !================================================================================================================================
  !

  !>Checks if the optimisation vector matrices for the vector matrices exist.
  SUBROUTINE EquationsMatricesVector_OptimisationMatricesExists(vectorMatrices,optimisationMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations vector matrices to check the optimisation matrices for
    TYPE(EquationsMatricesOptimisationType), POINTER :: optimisationMatrices !<On exit, a pointer to the optimisation matrices in the specified vector equations matrices if they exist. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_OptimisationMatricesExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(optimisationMatrices)) CALL FlagError("Optimisation matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)
#endif    

    optimisationMatrices=>vectorMatrices%optimisationMatrices
       
    EXITS("EquationsMatricesVector_OptimisationMatricesExists")
    RETURN
999 NULLIFY(optimisationMatrices)
998 ERRORS("EquationsMatricesVector_OptimisationMatricesExists",err,error)
    EXITS("EquationsMatricesVector_OptimisationMatricesExists")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_OptimisationMatricesExists

  !
  !================================================================================================================================
  !

  !>Gets the optimisation vector matrices for the vector matrices.
  SUBROUTINE EquationsMatricesVector_OptimisationMatricesGet(vectorMatrices,optimisationMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations vector matrices to get the optimisation matrices for
    TYPE(EquationsMatricesOptimisationType), POINTER :: optimisationMatrices !<On exit, a pointer to the optimisation matrices in the specified vector equations matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_OptimisationMatricesGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(optimisationMatrices)) CALL FlagError("Optimisation matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)
#endif    

    optimisationMatrices=>vectorMatrices%optimisationMatrices

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(optimisationMatrices))  &
      & CALL FlagError("Optimisation matrices is not associated for the vector matrices.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesVector_OptimisationMatricesGet")
    RETURN
999 NULLIFY(optimisationMatrices)
998 ERRORS("EquationsMatricesVector_OptimisationMatricesGet",err,error)
    EXITS("EquationsMatricesVector_OptimisationMatricesGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_OptimisationMatricesGet

  !
  !================================================================================================================================
  !

  !>Gets the number of rows in the vector matrices.
  SUBROUTINE EquationsMatricesVector_NumberOfRowsGet(vectorMatrices,numberOfRows,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations vector matrices to get the number of rows for
    INTEGER(INTG), INTENT(OUT) :: numberOfRows !<On exit, the number of rows in the specified vector equations matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_NumberOfRowsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)
#endif    

    numberOfRows=vectorMatrices%numberOfRows
       
    EXITS("EquationsMatricesVector_NumberOfRowsGet")
    RETURN
999 ERRORS("EquationsMatricesVector_NumberOfRowsGet",err,error)
    EXITS("EquationsMatricesVector_NumberOfRowsGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_NumberOfRowsGet

  !
  !================================================================================================================================
  !

  !>Gets the number of global rows in the vector matrices.
  SUBROUTINE EquationsMatricesVector_NumberOfGlobalRowsGet(vectorMatrices,numberOfGlobalRows,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations vector matrices to get the number of global rows for
    INTEGER(INTG), INTENT(OUT) :: numberOfGlobalRows !<On exit, the number of global rows in the specified vector equations matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_NumberOfGlobalRowsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)
#endif    

    numberOfGlobalRows=vectorMatrices%numberOfGlobalRows
       
    EXITS("EquationsMatricesVector_NumberOfGlobalRowsGet")
    RETURN
999 ERRORS("EquationsMatricesVector_NumberOfGlobalRowsGet",err,error)
    EXITS("EquationsMatricesVector_NumberOfGlobalRowsGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_NumberOfGlobalRowsGet

  !
  !================================================================================================================================
  !

  !>Checks that the RHS vector exists for the vector matrices.
  SUBROUTINE EquationsMatricesVector_RHSVectorExists(vectorMatrices,equationsMatricesRHS,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations vector matrices to check the rhs vector exists for
    TYPE(EquationsMatricesRHSType), POINTER :: equationsMatricesRHS !<On exit, a pointer to the RHS vector in the specified vector equations matrices if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_RHSVectorExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(equationsMatricesRHS)) CALL FlagError("RHS vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)
#endif    

    equationsMatricesRHS=>vectorMatrices%rhsVector
       
    EXITS("EquationsMatricesVector_RHSVectorExists")
    RETURN
999 NULLIFY(equationsMatricesRHS)
998 ERRORSEXITS("EquationsMatricesVector_RHSVectorExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_RHSVectorExists

  !
  !================================================================================================================================
  !

  !>Gets the RHS vector for an vector matrices.
  SUBROUTINE EquationsMatricesVector_RHSVectorGet(vectorMatrices,equationsMatricesRHS,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations vector matrices to get the rhs vector exists for
    TYPE(EquationsMatricesRHSType), POINTER :: equationsMatricesRHS !<On exit, a pointer to the RHS vector in the specified vector equations matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_RHSVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(equationsMatricesRHS)) CALL FlagError("RHS vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)
#endif    

    equationsMatricesRHS=>vectorMatrices%rhsVector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equationsMatricesRHS)) CALL FlagError("RHS vector is not associated for the vector matrices.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesVector_RHSVectorGet")
    RETURN
999 NULLIFY(equationsMatricesRHS)
998 ERRORSEXITS("EquationsMatricesVector_RHSVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_RHSVectorGet

  !
  !================================================================================================================================
  !

  !>Checks that the source vectors exists for the vector matrices.
  SUBROUTINE EquationsMatricesVector_SourceVectorsExists(vectorMatrices,sourceVectors,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations vector matrices to check the source vectors exists for
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors !<On exit, a pointer to the source vectors in the specified vector equations matrices if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_SourceVectorsExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(sourceVectors)) CALL FlagError("Source vectors is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)
#endif    

    sourceVectors=>vectorMatrices%sourceVectors

   EXITS("EquationsMatricesVector_SourceVectorsExists")
    RETURN
999 NULLIFY(sourceVectors)
998 ERRORSEXITS("EquationsMatricesVector_SourceVectorsExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_SourceVectorsExists

  !
  !================================================================================================================================
  !

  !>Gets the source vectors for the vector matrices.
  SUBROUTINE EquationsMatricesVector_SourceVectorsGet(vectorMatrices,sourceVectors,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations vector matrices to get the source vectorS for
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors !<On exit, a pointer to the source vectors in the specified vector equations matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_SourceVectorsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(sourceVectors)) CALL FlagError("Source vectors is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)
#endif    

    sourceVectors=>vectorMatrices%sourceVectors

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(sourceVectors)) &
      & CALL FlagError("Source vectors is not associated for the vector matrices.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesVector_SourceVectorsGet")
    RETURN
999 NULLIFY(sourceVectors)
998 ERRORSEXITS("EquationsMatricesVector_SourceVectorsGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_SourceVectorsGet

  !
  !================================================================================================================================
  !

  !>Gets the total number of rows in the vector matrices.
  SUBROUTINE EquationsMatricesVector_TotalNumberOfRowsGet(vectorMatrices,totalNumberOfRows,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations vector matrices to get the total number of rows for
    INTEGER(INTG), INTENT(OUT) :: totalNumberOfRows !<On exit, the total number of rows in the specified vector equations matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_TotalNumberOfRowsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)
#endif    

    totalNumberOfRows=vectorMatrices%totalNumberOfRows
       
    EXITS("EquationsMatricesVector_TotalNumberOfRowsGet")
    RETURN
999 ERRORS("EquationsMatricesVector_TotalNumberOfRowsGet",err,error)
    EXITS("EquationsMatricesVector_TotalNumberOfRowsGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_TotalNumberOfRowsGet

  !
  !================================================================================================================================
  !

  !>Gets the vector equations for an vector equations matrices.
  SUBROUTINE EquationsMatricesVector_VectorEquationsGet(vectorMatrices,vectorEquations,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations vector matrices to get the vector equations for
    TYPE(EquationsVectorType), POINTER :: vectorEquations !<On exit, a pointer to the vector equations in the specified vector equations matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_VectorEquationsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)
#endif    

    vectorEquations=>vectorMatrices%vectorEquations

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is not associated for the vector matrices.", &
      & err,error,*999)
#endif    
       
    EXITS("EquationsMatricesVector_VectorEquationsGet")
    RETURN
999 NULLIFY(vectorEquations)
998 ERRORSEXITS("EquationsMatricesVector_VectorEquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_VectorEquationsGet

  !
  !================================================================================================================================
  !

  !>Gets the vector mapping for vector matrices.
  SUBROUTINE EquationsMatricesVector_VectorMappingGet(vectorMatrices,vectorMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations vector matrices to get the vector mapping for
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<On exit, a pointer to the vector mapping for the specified vector equations matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_VectorMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(vectorMatrices%vectorEquations)) &
      & CALL FlagError("Vector matrices vector equations is not associated.",err,error,*999)
#endif    
    
    vectorMapping=>vectorMatrices%vectorEquations%vectorMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated for the vector matrices.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesVector_VectorMappingGet")
    RETURN
999 NULLIFY(vectorMapping)
998 ERRORSEXITS("EquationsMatricesVector_VectorMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_VectorMappingGet

  !
  !================================================================================================================================
  !

  !>Gets the distributed matrix for a equations matrix
  SUBROUTINE EquationsMatrix_DistributedMatrixGet(equationsMatrix,distributedMatrix,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to get the distributed matrix for
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<On exit, a pointer to the distributed matrix for the equations matrix. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatrix_DistributedMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is not associated.",err,error,*999)
#endif    

    distributedMatrix=>equationsMatrix%matrix

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(distributedMatrix)) &
      & CALL FlagError("Equations matrix distributed matrix is not associated.",err,error,*999)
#endif    
       
    EXITS("EquationsMatrix_DistributedMatrixGet")
    RETURN
999 NULLIFY(distributedMatrix)
998 ERRORSEXITS("EquationsMatrix_DistributedMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrix_DistributedMatrixGet

  !
  !================================================================================================================================
  !

  !>Checks if the dynamic matrices for a equations matrix exists
  SUBROUTINE EquationsMatrix_DynamicMatricesExists(equationsMatrix,dynamicMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to check the dynamic matrices for
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices !<On exit, a pointer to the dynamic matrices for the equations matrix if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatrix_DynamicMatricesExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dynamicMatrices)) CALL FlagError("Dynamic matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is not associated.",err,error,*999)
#endif    

    dynamicMatrices=>equationsMatrix%dynamicMatrices

    EXITS("EquationsMatrix_DynamicMatricesExists")
    RETURN
999 NULLIFY(dynamicMatrices)
998 ERRORSEXITS("EquationsMatrix_DynamicMatricesExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrix_DynamicMatricesExists

  !
  !================================================================================================================================
  !

  !>Gets the dynamic matrices for a equations matrix
  SUBROUTINE EquationsMatrix_DynamicMatricesGet(equationsMatrix,dynamicMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to get the dynamic matrices for
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices !<On exit, a pointer to the dynamic matrices for the equations matrix. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatrix_DynamicMatricesGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dynamicMatrices)) CALL FlagError("Dynamic matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is not associated.",err,error,*999)
#endif    

    dynamicMatrices=>equationsMatrix%dynamicMatrices

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dynamicMatrices)) &
      & CALL FlagError("Equations matrix dynamic matrices is not associated.",err,error,*999)
#endif    
       
    EXITS("EquationsMatrix_DynamicMatricesGet")
    RETURN
999 NULLIFY(dynamicMatrices)
998 ERRORSEXITS("EquationsMatrix_DynamicMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrix_DynamicMatricesGet

  !
  !================================================================================================================================
  !

  !>Outputs the element matrix information for an equations matrix.
  SUBROUTINE EquationsMatrix_ElementMatrixOutput(id,equationsMatrix,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the ouptut stream
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to output the element matrix for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("EquationsMatrix_ElementMatrixOutput",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is not associated.",err,error,*999)

    CALL ElementMatrix_Output(id,equationsMatrix%elementMatrix,err,error,*999)
         
    EXITS("EquationsMatrix_ElementMatrixOutput")
    RETURN
999 ERRORSEXITS("EquationsMatrix_ElementMatrixOutput",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrix_ElementMatrixOutput
  
  !
  !================================================================================================================================
  !

  !>Gets the first assembly flag for an equations matrix
  SUBROUTINE EquationsMatrix_FirstAssemblyGet(equationsMatrix,firstAssembly,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to get the first assembly flag for
    LOGICAL, INTENT(OUT) :: firstAssembly !<On exit, the first assembly flag of the equations matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatrix_FirstAssemblyGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is not associated.",err,error,*999)
#endif    

    firstAssembly=equationsMatrix%firstAssembly
       
    EXITS("EquationsMatrix_FirstAssemblyGet")
    RETURN
999 ERRORSEXITS("EquationsMatrix_FirstAssemblyGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrix_FirstAssemblyGet

  !
  !================================================================================================================================
  !

  !>Checks if the linear matrices for a equations matrix exists
  SUBROUTINE EquationsMatrix_LinearMatricesExists(equationsMatrix,linearMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to checks the linear matrices for
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices !<On exit, a pointer to the linear matrices for the equations matrix if they exist. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatrix_LinearMatricesExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linearMatrices)) CALL FlagError("Linear matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is not associated.",err,error,*999)
#endif    

    linearMatrices=>equationsMatrix%linearMatrices
       
    EXITS("EquationsMatrix_LinearMatricesExists")
    RETURN
999 NULLIFY(linearMatrices)
998 ERRORSEXITS("EquationsMatrix_LinearMatricesExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrix_LinearMatricesExists

  !
  !================================================================================================================================
  !

  !>Gets the linear matrices for a equations matrix
  SUBROUTINE EquationsMatrix_LinearMatricesGet(equationsMatrix,linearMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to get the linear matrices for
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices !<On exit, a pointer to the linear matrices for the equations matrix. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatrix_LinearMatricesGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linearMatrices)) CALL FlagError("Linear matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is not associated.",err,error,*999)
#endif    

    linearMatrices=>equationsMatrix%linearMatrices

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(linearMatrices)) &
      & CALL FlagError("Equations matrix linear matrices is not associated.",err,error,*999)
#endif    
       
    EXITS("EquationsMatrix_LinearMatricesGet")
    RETURN
999 NULLIFY(linearMatrices)
998 ERRORSEXITS("EquationsMatrix_LinearMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrix_LinearMatricesGet

  !
  !================================================================================================================================
  !

  !>Gets the lumped flag for an equations matrix
  SUBROUTINE EquationsMatrix_LumpedFlagGet(equationsMatrix,lumpedFlag,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to get the lumped flag for
    LOGICAL, INTENT(OUT) :: lumpedFlag !<On exit, the lumped flag of the equations matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatrix_LumpedFlagGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is not associated.",err,error,*999)
#endif    

    lumpedFlag=equationsMatrix%lumped
       
    EXITS("EquationsMatrix_LumpedFlagGet")
    RETURN
999 ERRORSEXITS("EquationsMatrix_LumpedFlagGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrix_LumpedFlagGet

  !
  !================================================================================================================================
  !

  !>Gets the matrix coefficient for an equations matrix
  SUBROUTINE EquationsMatrix_MatrixCoefficientGet(equationsMatrix,matrixCoefficient,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to get the matrix coefficient for
    REAL(DP), INTENT(OUT) :: matrixCoefficient !<On exit, the matrix coefficient of the equations matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatrix_MatrixCoefficientGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is not associated.",err,error,*999)
#endif    

    matrixCoefficient=equationsMatrix%matrixCoefficient
       
    EXITS("EquationsMatrix_MatrixCoefficientGet")
    RETURN
999 ERRORSEXITS("EquationsMatrix_MatrixCoefficientGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrix_MatrixCoefficientGet

  !
  !================================================================================================================================
  !

  !>Gets the matrix number for an equations matrix
  SUBROUTINE EquationsMatrix_MatrixNumberGet(equationsMatrix,matrixNumber,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to get the matrix number for
    INTEGER(INTG), INTENT(OUT) :: matrixNumber !<On exit, the matrix number of the equations matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatrix_MatrixNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is not associated.",err,error,*999)
#endif    

    matrixNumber=equationsMatrix%matrixNumber
       
    EXITS("EquationsMatrix_MatrixNumberGet")
    RETURN
999 ERRORSEXITS("EquationsMatrix_MatrixNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrix_MatrixNumberGet

  !
  !================================================================================================================================
  !

  !>Outputs the nodal matrix information for an equations matrix.
  SUBROUTINE EquationsMatrix_NodalMatrixOutput(id,equationsMatrix,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the ouptut stream
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to output the nodal matrix for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("EquationsMatrix_NodalMatrixOutput",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is not associated.",err,error,*999)

    CALL NodalMatrix_Output(id,equationsMatrix%nodalMatrix,err,error,*999)
         
    EXITS("EquationsMatrix_NodalMatrixOutput")
    RETURN
999 ERRORSEXITS("EquationsMatrix_NodalMatrixOutput",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrix_NodalMatrixOutput
  
  !
  !================================================================================================================================
  !

  !>Gets the number of columns for an equations matrix
  SUBROUTINE EquationsMatrix_NumberOfColumnsGet(equationsMatrix,numberOfColumns,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to get the number of columns for
    INTEGER(INTG), INTENT(OUT) :: numberOfColumns !<On exit, the number of columns of the equations matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatrix_NumberOfColumnsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is not associated.",err,error,*999)
#endif    

    numberOfColumns=equationsMatrix%numberOfColumns
       
    EXITS("EquationsMatrix_NumberOfColumnsGet")
    RETURN
999 ERRORSEXITS("EquationsMatrix_NumberOfColumnsGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrix_NumberOfColumnsGet

  !
  !================================================================================================================================
  !

  !>Gets the storage type for an equations matrix
  SUBROUTINE EquationsMatrix_StorageTypeGet(equationsMatrix,storageType,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to get the storage type for
    INTEGER(INTG), INTENT(OUT) :: storageType !<On exit, the storage type of the equations matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatrix_StorageTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is not associated.",err,error,*999)
#endif    

    storageType=equationsMatrix%storageType
       
    EXITS("EquationsMatrix_StorageTypeGet")
    RETURN
999 ERRORSEXITS("EquationsMatrix_StorageTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrix_StorageTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the structure type for an equations matrix
  SUBROUTINE EquationsMatrix_StructureTypeGet(equationsMatrix,structureType,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to get the structure type for
    INTEGER(INTG), INTENT(OUT) :: structureType !<On exit, the structure type of the equations matrix. \see EquationsMatricesRoutines_EquationsMatrixStructureType
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatrix_StructureTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is not associated.",err,error,*999)
#endif    

    structureType=equationsMatrix%structureType
       
    EXITS("EquationsMatrix_StructureTypeGet")
    RETURN
999 ERRORSEXITS("EquationsMatrix_StructureTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrix_StructureTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the symmetry flag for an equations matrix
  SUBROUTINE EquationsMatrix_SymmetryFlagGet(equationsMatrix,symmetryFlag,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to get the symmetry flag for
    LOGICAL, INTENT(OUT) :: symmetryFlag !<On exit, the symmetry flag of the equations matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatrix_SymmetryFlagGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is not associated.",err,error,*999)
#endif    

    symmetryFlag=equationsMatrix%symmetric
       
    EXITS("EquationsMatrix_SymmetryFlagGet")
    RETURN
999 ERRORSEXITS("EquationsMatrix_SymmetryFlagGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrix_SymmetryFlagGet

  !
  !================================================================================================================================
  !

  !>Checks if the temp distributed vector for a equations matrix exists
  SUBROUTINE EquationsMatrix_TempDistributedVectorExists(equationsMatrix,tempDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to check the temp distributed vector existance for
    TYPE(DistributedVectorType), POINTER :: tempDistributedVector !<On exit, a pointer to the temp distributed vector for the equations matrix if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatrix_TempDistributedVectorExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(tempDistributedVector)) CALL FlagError("Temporary distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is not associated.",err,error,*999)
#endif    

    tempDistributedVector=>equationsMatrix%tempVector

    EXITS("EquationsMatrix_TempDistributedVectorExists")
    RETURN
999 NULLIFY(tempDistributedVector)
998 ERRORSEXITS("EquationsMatrix_TempDistributedVectorExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrix_TempDistributedVectorExists

  !
  !================================================================================================================================
  !

  !>Gets the temp distributed vector for a equations matrix
  SUBROUTINE EquationsMatrix_TempDistributedVectorGet(equationsMatrix,tempDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to get the temp distributed vector for
    TYPE(DistributedVectorType), POINTER :: tempDistributedVector !<On exit, a pointer to the temp distributed vector for the equations matrix. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatrix_TempDistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(tempDistributedVector)) CALL FlagError("Temporary distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is not associated.",err,error,*999)
#endif    

    tempDistributedVector=>equationsMatrix%tempVector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(tempDistributedVector)) &
      & CALL FlagError("Equations matrix temporary distributed vector is not associated.",err,error,*999)
#endif    
       
    EXITS("EquationsMatrix_TempDistributedVectorGet")
    RETURN
999 NULLIFY(tempDistributedVector)
998 ERRORSEXITS("EquationsMatrix_TempDistributedVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrix_TempDistributedVectorGet

  !
  !================================================================================================================================
  !

  !>Gets the update matrix flag for an equations matrix
  SUBROUTINE EquationsMatrix_UpdateMatrixGet(equationsMatrix,updateMatrix,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to get the update matrix flag for
    LOGICAL, INTENT(OUT) :: updateMatrix !<On exit, the update matrix flag of the equations matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatrix_UpdateMatrixGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is not associated.",err,error,*999)
#endif    

    updateMatrix=equationsMatrix%updateMatrix
       
    EXITS("EquationsMatrix_UpdateMatrixGet")
    RETURN
999 ERRORSEXITS("EquationsMatrix_UpdateMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrix_UpdateMatrixGet

  !
  !================================================================================================================================
  !

  !>Sets/changes the update matrix flag for an equations matrix
  SUBROUTINE EquationsMatrix_UpdateMatrixSet(equationsMatrix,updateMatrix,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix to set the update matrix flag for
    LOGICAL, INTENT(IN) :: updateMatrix !The update matrix flag of the equations matrix to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatrix_UpdateMatrixSet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is not associated.",err,error,*999)
#endif    

    equationsMatrix%updateMatrix=updateMatrix
       
    EXITS("EquationsMatrix_UpdateMatrixSet")
    RETURN
999 ERRORSEXITS("EquationsMatrix_UpdateMatrixSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrix_UpdateMatrixSet

  !
  !================================================================================================================================
  !

  !>Gets the calculation type for an Jacobian matrix
  SUBROUTINE JacobianMatrix_CalculationTypeGet(jacobianMatrix,calculationType,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix !<A pointer to the Jacobian matrix to get the calculation type for
    INTEGER(INTG), INTENT(OUT) :: calculationType !<On exit, the calculation type of the Jacobian matrix. \see EquationsMatricesRoutines_JacobianCalculationTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("JacobianMatrix_CalculationTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(jacobianMatrix)) CALL FlagError("Jacobian matrix is not associated.",err,error,*999)
#endif    

    calculationType=jacobianMatrix%jacobianCalculationType
       
    EXITS("JacobianMatrix_CalculationTypeGet")
    RETURN
999 ERRORSEXITS("JacobianMatrix_CalculationTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE JacobianMatrix_CalculationTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the distributed matrix for a Jacobian matrix
  SUBROUTINE JacobianMatrix_DistributedMatrixGet(jacobianMatrix,distributedMatrix,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix !<A pointer to the Jacobian matrix to get the distributed matrix for
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<On exit, a pointer to the distributed matrix for the Jacobian matrix. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("JacobianMatrix_DistributedMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(distributedMatrix)) CALL FlagError("Distributed matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(jacobianMatrix)) CALL FlagError("Jacobian matrix is not associated.",err,error,*999)
#endif    

    distributedMatrix=>jacobianMatrix%jacobian

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(distributedMatrix)) &
      & CALL FlagError("Jacobian matrix distributed matrix is not associated.",err,error,*999)
#endif    
       
    EXITS("JacobianMatrix_DistributedMatrixGet")
    RETURN
999 NULLIFY(distributedMatrix)
998 ERRORSEXITS("JacobianMatrix_DistributedMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE JacobianMatrix_DistributedMatrixGet

  !
  !================================================================================================================================
  !

  !>Outputs the element matrix information for an Jacobian matrix.
  SUBROUTINE JacobianMatrix_ElementMatrixOutput(id,jacobianMatrix,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the ouptut stream
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix !<A pointer to the Jacobian matrix to output the element matrix for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("JacobianMatrix_ElementMatrixOutput",err,error,*999)

    IF(.NOT.ASSOCIATED(jacobianMatrix)) CALL FlagError("Jacobian matrix is not associated.",err,error,*999)
     
    CALL ElementMatrix_Output(id,jacobianMatrix%elementJacobian,err,error,*999)
    
    EXITS("JacobianMatrix_ElementMatrixOutput")
    RETURN
999 ERRORSEXITS("JacobianMatrix_ElementMatrixOutput",err,error)
    RETURN 1
    
  END SUBROUTINE JacobianMatrix_ElementMatrixOutput
  
  !
  !================================================================================================================================
  !

  !>Gets the finite difference step size for an jacobian matrix
  SUBROUTINE JacobianMatrix_FiniteDifferenceStepSizeGet(jacobianMatrix,finiteDifferenceStepSize,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix !<A pointer to the jacobian matrix to get the finite difference step size for
    REAL(DP), INTENT(OUT) :: finiteDifferenceStepSize !<On exit, the finite difference step size  of the jacobian matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("JacobianMatrix_FiniteDifferenceStepSizeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(jacobianMatrix)) CALL FlagError("Jacobian matrix is not associated.",err,error,*999)
#endif    

    finiteDifferenceStepSize=jacobianMatrix%jacobianFiniteDifferenceStepSize
       
    EXITS("JacobianMatrix_FiniteDifferenceStepSizeGet")
    RETURN
999 ERRORSEXITS("JacobianMatrix_FiniteDifferenceStepSizeGet",err,error)
    RETURN 1
    
  END SUBROUTINE JacobianMatrix_FiniteDifferenceStepSizeGet

  !
  !================================================================================================================================
  !

  !>Gets the first assembly flag for an Jacobian matrix
  SUBROUTINE JacobianMatrix_FirstAssemblyGet(jacobianMatrix,firstAssembly,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix !<A pointer to the Jacobian matrix to get the first assembly flag for
    LOGICAL, INTENT(OUT) :: firstAssembly !<On exit, the first assembly flag of the Jacobian matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("JacobianMatrix_FirstAssemblyGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(jacobianMatrix)) CALL FlagError("Jacobian matrix is not associated.",err,error,*999)
#endif    

    firstAssembly=jacobianMatrix%firstAssembly
       
    EXITS("JacobianMatrix_FirstAssemblyGet")
    RETURN
999 ERRORSEXITS("JacobianMatrix_FirstAssemblyGet",err,error)
    RETURN 1
    
  END SUBROUTINE JacobianMatrix_FirstAssemblyGet

  !
  !================================================================================================================================
  !

  !>Gets the matrix coefficient for an jacobian matrix
  SUBROUTINE JacobianMatrix_MatrixCoefficientGet(jacobianMatrix,matrixCoefficient,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix !<A pointer to the jacobian matrix to get the matrix coefficient for
    REAL(DP), INTENT(OUT) :: matrixCoefficient !<On exit, the matrix coefficient of the jacobian matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("JacobianMatrix_MatrixCoefficientGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(jacobianMatrix)) CALL FlagError("Jacobian matrix is not associated.",err,error,*999)
#endif    

    matrixCoefficient=jacobianMatrix%jacobianCoefficient
       
    EXITS("JacobianMatrix_MatrixCoefficientGet")
    RETURN
999 ERRORSEXITS("JacobianMatrix_MatrixCoefficientGet",err,error)
    RETURN 1
    
  END SUBROUTINE JacobianMatrix_MatrixCoefficientGet

  !
  !================================================================================================================================
  !

  !>Gets the matrix number for an jacobian matrix
  SUBROUTINE JacobianMatrix_MatrixNumberGet(jacobianMatrix,matrixNumber,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix !<A pointer to the jacobian matrix to get the matrix number for
    INTEGER(INTG), INTENT(OUT) :: matrixNumber !<On exit, the matrix number of the jacobian matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("JacobianMatrix_MatrixNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(jacobianMatrix)) CALL FlagError("Jacobian matrix is not associated.",err,error,*999)
#endif    

    matrixNumber=jacobianMatrix%jacobianNumber
       
    EXITS("JacobianMatrix_MatrixNumberGet")
    RETURN
999 ERRORSEXITS("JacobianMatrix_MatrixNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE JacobianMatrix_MatrixNumberGet

  !
  !================================================================================================================================
  !

  !>Outputs the nodal matrix information for an Jacobian matrix.
  SUBROUTINE JacobianMatrix_NodalMatrixOutput(id,jacobianMatrix,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the ouptut stream
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix !<A pointer to the Jacobian matrix to output the nodal matrix for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("JacobianMatrix_NodalMatrixOutput",err,error,*999)

    IF(.NOT.ASSOCIATED(jacobianMatrix)) CALL FlagError("Jacobian matrix is not associated.",err,error,*999)
     
    CALL NodalMatrix_Output(id,jacobianMatrix%nodalJacobian,err,error,*999)
    
    EXITS("JacobianMatrix_NodalMatrixOutput")
    RETURN
999 ERRORSEXITS("JacobianMatrix_NodalMatrixOutput",err,error)
    RETURN 1
    
  END SUBROUTINE JacobianMatrix_NodalMatrixOutput
  
  !
  !================================================================================================================================
  !

  !>Gets the number of columns for an jacobian matrix
  SUBROUTINE JacobianMatrix_NumberOfColumnsGet(jacobianMatrix,numberOfColumns,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix !<A pointer to the jacobian matrix to get the number of columns for
    INTEGER(INTG), INTENT(OUT) :: numberOfColumns !<On exit, the number of columns of the jacobian matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("JacobianMatrix_NumberOfColumnsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(jacobianMatrix)) CALL FlagError("Jacobian matrix is not associated.",err,error,*999)
#endif    

    numberOfColumns=jacobianMatrix%numberOfColumns
       
    EXITS("JacobianMatrix_NumberOfColumnsGet")
    RETURN
999 ERRORSEXITS("JacobianMatrix_NumberOfColumnsGet",err,error)
    RETURN 1
    
  END SUBROUTINE JacobianMatrix_NumberOfColumnsGet

  !
  !================================================================================================================================
  !

  !>Gets the residual vector for a Jacobian matrix
  SUBROUTINE JacobianMatrix_ResidualVectorGet(jacobianMatrix,residualVector,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix !<A pointer to the Jacobian matrix to get the residual vector for
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<On exit, a pointer to the residual vectorfor the Jacobian matrix. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("JacobianMatrix_ResidualVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(residualVector)) CALL FlagError("Residual vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(jacobianMatrix)) CALL FlagError("Jacobian matrix is not associated.",err,error,*999)
#endif    

    residualVector=>jacobianMatrix%residualVector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(residualVector)) &
      & CALL FlagError("Jacobian matrix residual vector is not associated.",err,error,*999)
#endif    
       
    EXITS("JacobianMatrix_ResidualVectorGet")
    RETURN
999 NULLIFY(residualVector)
998 ERRORSEXITS("JacobianMatrix_ResidualVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE JacobianMatrix_ResidualVectorGet

  !
  !================================================================================================================================
  !

  !>Gets the storage type for an jacobian matrix
  SUBROUTINE JacobianMatrix_StorageTypeGet(jacobianMatrix,storageType,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix !<A pointer to the jacobian matrix to get the storage type for
    INTEGER(INTG), INTENT(OUT) :: storageType !<On exit, the storage type of the jacobian matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("JacobianMatrix_StorageTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(jacobianMatrix)) CALL FlagError("Jacobian matrix is not associated.",err,error,*999)
#endif    

    storageType=jacobianMatrix%storageType
       
    EXITS("JacobianMatrix_StorageTypeGet")
    RETURN
999 ERRORSEXITS("JacobianMatrix_StorageTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE JacobianMatrix_StorageTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the structure type for an jacobian matrix
  SUBROUTINE JacobianMatrix_StructureTypeGet(jacobianMatrix,structureType,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix !<A pointer to the jacobian matrix to get the structure type for
    INTEGER(INTG), INTENT(OUT) :: structureType !<On exit, the structure type of the jacobian matrix. \see JacobianMatricesRoutines_JacobianMatrixStructureType
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("JacobianMatrix_StructureTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(jacobianMatrix)) CALL FlagError("Jacobian matrix is not associated.",err,error,*999)
#endif    

    structureType=jacobianMatrix%structureType
       
    EXITS("JacobianMatrix_StructureTypeGet")
    RETURN
999 ERRORSEXITS("JacobianMatrix_StructureTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE JacobianMatrix_StructureTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the symmetry flag for an jacobian matrix
  SUBROUTINE JacobianMatrix_SymmetryFlagGet(jacobianMatrix,symmetryFlag,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix !<A pointer to the jacobian matrix to get the symmetry flag for
    LOGICAL, INTENT(OUT) :: symmetryFlag !<On exit, the symmetry flag of the jacobian matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("JacobianMatrix_SymmetryFlagGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(jacobianMatrix)) CALL FlagError("Jacobian matrix is not associated.",err,error,*999)
#endif    

    symmetryFlag=jacobianMatrix%symmetric
       
    EXITS("JacobianMatrix_SymmetryFlagGet")
    RETURN
999 ERRORSEXITS("JacobianMatrix_SymmetryFlagGet",err,error)
    RETURN 1
    
  END SUBROUTINE JacobianMatrix_SymmetryFlagGet

  !
  !================================================================================================================================
  !

  !>Gets the update matrix flag for an jacobian matrix
  SUBROUTINE JacobianMatrix_UpdateMatrixGet(jacobianMatrix,updateMatrix,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix !<A pointer to the jacobian matrix to get the update matrix flag for
    LOGICAL, INTENT(OUT) :: updateMatrix !<On exit, the update matrix flag of the jacobian matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("JacobianMatrix_UpdateMatrixGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(jacobianMatrix)) CALL FlagError("Jacobian matrix is not associated.",err,error,*999)
#endif    

    updateMatrix=jacobianMatrix%updateJacobian
       
    EXITS("JacobianMatrix_UpdateMatrixGet")
    RETURN
999 ERRORSEXITS("JacobianMatrix_UpdateMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE JacobianMatrix_UpdateMatrixGet

  !
  !================================================================================================================================
  !

  !>Sets/changes the update matrix flag for a Jacobian matrix
  SUBROUTINE JacobianMatrix_UpdateMatrixSet(jacobianMatrix,updateMatrix,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix !<A pointer to the Jacobian matrix to set the update matrix flag for
    LOGICAL, INTENT(IN) :: updateMatrix !<The update matrix flag of the Jacobian matrix to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("JacobianMatrix_UpdateMatrixSet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(jacobianMatrix)) CALL FlagError("Jacobian matrix is not associated.",err,error,*999)
#endif    

    jacobianMatrix%updateJacobian=updateMatrix
       
    EXITS("JacobianMatrix_UpdateMatrixSet")
    RETURN
999 ERRORSEXITS("JacobianMatrix_UpdateMatrixSet",err,error)
    RETURN 1
    
  END SUBROUTINE JacobianMatrix_UpdateMatrixSet

  !
  !================================================================================================================================
  !

  !>Outputs the nodal matrix information.
  SUBROUTINE NodalMatrix_Output(id,nodalMatrix,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the ouptut stream
    TYPE(NodalMatrixType), INTENT(IN) :: nodalMatrix !<The nodal matrix to output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("NodalMatrix_Output",err,error,*999)

    CALL WriteStringValue(id,"  Number of rows = ",nodalMatrix%numberOfRows,err,error,*999)
    CALL WriteStringValue(id,"  Number of columns = ",nodalMatrix%numberOfColumns,err,error,*999)
    CALL WriteStringValue(id,"  Maximum number of rows = ",nodalMatrix%maxNumberOfRows,err,error,*999)
    CALL WriteStringValue(id,"  Maximum number of columns = ",nodalMatrix%maxNumberOfColumns,err,error,*999)
    CALL WriteStringVector(id,1,1,nodalMatrix%numberOfRows,8,8,nodalMatrix%rowDOFS, &
      & '("  Row dofs         :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
    CALL WriteStringVector(id,1,1,nodalMatrix%numberOfColumns,8,8,nodalMatrix%columnDOFS, &
      & '("  Column dofs      :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
    CALL WriteStringMatrix(id,1,1,nodalMatrix%numberOfRows,1,1,nodalMatrix%numberOfColumns,8,8, &
      & nodalMatrix%matrix(1:nodalMatrix%numberOfRows,1:nodalMatrix%numberOfColumns),WRITE_STRING_MATRIX_NAME_AND_INDICES, &
      & '("  Matrix','(",I2,",:)','     :",8(X,E13.6))','(20X,8(X,E13.6))',err,error,*999)
    
    EXITS("NodalMatrix_Output")
    RETURN
999 ERRORSEXITS("NodalMatrix_Output",err,error)
    RETURN 1
    
  END SUBROUTINE NodalMatrix_Output
  
  !
  !================================================================================================================================
  !

  !>Outputs the nodal vector information .
  SUBROUTINE NodalVector_Output(id,nodalVector,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the ouptut stream
    TYPE(NodalVectorType), INTENT(IN) :: nodalVector !<The element vector to output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("NodalVector_Output",err,error,*999)

    CALL WriteStringValue(id,"  Number of rows = ",nodalVector%numberOfRows,err,error,*999)
    CALL WriteStringValue(id,"  Maximum number of rows = ",nodalVector%maxNumberOfRows,err,error,*999)
    CALL WriteStringVector(id,1,1,nodalVector%numberOfRows,8,8,nodalVector%rowDOFS, &
      & '("  Row dofs         :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
    CALL WriteStringVector(id,1,1,nodalVector%numberOfRows,8,8,nodalVector%vector, &
      & '("  Vector(:)        :",8(X,E13.6))','(20X,8(X,E13.6))',err,error,*999)
    
    EXITS("NodalVector_Output")
    RETURN
999 ERRORSEXITS("NodalVector_Output",err,error)
    RETURN 1
    
  END SUBROUTINE NodalVector_Output
  
  !
  !================================================================================================================================
  !

END MODULE EquationsMatricesAccessRoutines

