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
 
  !Module types

  !Module variables

  !Interfaces

  PUBLIC EQUATION_MATRIX_LINEAR,EQUATIONS_MATRIX_DYNAMIC

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

  PUBLIC EquationsMatricesDynamic_DynamicMappingGet
  
  PUBLIC EquationsMatricesDynamic_EquationsMatrixGet
  
  PUBLIC EquationsMatricesDynamic_VectorMatricesGet

  PUBLIC EquationsMatricesLinear_EquationsMatrixGet
  
  PUBLIC EquationsMatricesLinear_LinearMappingGet

  PUBLIC EquationsMatricesLinear_VectorMatricesGet

  PUBLIC EquationsMatricesNonlinear_NonlinearMappingGet
  
  PUBLIC EquationsMatricesNonlinear_VectorMatricesGet

  PUBLIC EquationsMatricesResidual_DistributedVectorGet

  PUBLIC EquationsMatriceResidual_JacobianMatrixGet

  PUBLIC EquationsMatriceResidual_NonlinearMatricesGet

  PUBLIC EquationsMatriceResidual_ResidualMappingGet

  PUBLIC EquationsMatricesRHS_DistributedVectorGet

  PUBLIC EquationsMatricesRHS_VectorMatricesGet

  PUBLIC EquationsMatricesScalar_EquationsScalarGet

  PUBLIC EquationsMatricesSource_DistributedVectorGet

  PUBLIC EquationsMatricesSource_SourceVectorsGet

  PUBLIC EquationsMatricesSources_SourceVectorGet

  PUBLIC EquationsMatricesSources_VectorMatricesGet

  PUBLIC EquationsMatricesVector_AssertIsFinish,EquationsMatricesVector_AssertNotFinished

  PUBLIC EquationsMatricesVector_DynamicMatricesExists

  PUBLIC EquationsMatricesVector_DynamicMatricesGet

  PUBLIC EquationsMatricesVector_LinearMatricesExists
  
  PUBLIC EquationsMatricesVector_LinearMatricesGet

  PUBLIC EquationsMatricesVector_NonlinearMatricesExists
  
  PUBLIC EquationsMatricesVector_NonlinearMatricesGet

  PUBLIC EquationsMatricesVector_RHSVectorExists
  
  PUBLIC EquationsMatricesVector_RHSVectorGet

  PUBLIC EquationsMatricesVector_SourceVectorsExists
  
  PUBLIC EquationsMatricesVector_SourceVectorsGet

  PUBLIC EquationsMatricesVector_VectorEquationsGet

  PUBLIC EquationsMatricesVector_VectorMappingGet

  PUBLIC EquationsMatrix_DistributedMatrixGet

  PUBLIC EquationsMatrix_DynamicMatricesExists

  PUBLIC EquationsMatrix_DynamicMatricesGet

  PUBLIC EquationsMatrix_LinearMatricesExists
  
  PUBLIC EquationsMatrix_LinearMatricesGet

  PUBLIC EquationsMatrix_TempDistributedVectorGet

  PUBLIC JacobianMatrix_DistributedMatrixGet

  PUBLIC JacobianMatrix_NonlinearMatricesGet


CONTAINS

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

  !>Gets the residual distributed vector for a equations matrices residual vector
  SUBROUTINE EquationsMatricesResidual_DistributedVectorGet(residualVector,residualDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer to the residual to get the  distributed vector for
    TYPE(DistributedVectorType), POINTER :: residualDistributedVector !<On exit, a pointer to the residual distributed vector for the equations matrices residual. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesResidual_DistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(residualDistributedVector)) CALL FlagError("Residual distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatricesRHS)) CALL FlagError("Equations matrices RHS is not associated.",err,error,*999)
#endif    

    residualDistributedVector=>residualVector%residual

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(residualDistributedVector)) &
      & CALL FlagError("The distributed vector is not associated for the residual vector.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesResidual_DistributedVectorGet")
    RETURN
999 NULLIFY(ressidualDistributedVector)
998 ERRORSEXITS("EquationsMatricesResidual_DistributedVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesResidual_DistributedVectorGet

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
        & ". The matrix index must be >= 1 and <= "//TRIM(NumberToVString(residualMatrices%numberOfJacobians,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    jacobianMatrix=>nonlinearMatrices%jacobians(matrixIdx)%ptr

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

  !>Gets the RHS distributed vector for a equations matrices RHS
  SUBROUTINE EquationsMatricesRHS_DistributedVectorGet(equationsMatricesRHS,rhsDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesRHSType), POINTER :: equationsMatricesRHS !<A pointer to the equations matrices RHS to get the RHS distributed vector for
    TYPE(DistributedVectorType), POINTER :: rhsDistributedVector !<On exit, a pointer to the RHS distributed vector for the equations matrices RHS. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesRHS_DistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(rhsDistributedVector)) CALL FlagError("RHS distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatricesRHS)) CALL FlagError("Equations matrices RHS is not associated.",err,error,*999)
#endif    

    rhsDistributedVector=>equationsMatricesRHS%vector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(rhsDistributedVector)) &
      & CALL FlagError("Equations matrices RHS distributed vector is not associated.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesRHS_DistributedVectorGet")
    RETURN
999 NULLIFY(rhsDistributedVector)
998 ERRORSEXITS("EquationsMatricesRHS_DistributedVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesRHS_DistributedVectorGet

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
999 NULLIFY(vectorMatriecs)
998 ERRORSEXITS("EquationsMatricesRHS_VectorMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesRHS_VectorMatricesGet

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

  !>Gets the source distributed vector for a equations matrices source vector
  SUBROUTINE EquationsMatricesSource_DistributedVectorGet(sourceVector,sourceDistributedVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<A pointer to the equations matrices source vector to get the source distributed vector for
    TYPE(DistributedVectorType), POINTER :: sourceDistributedVector !<On exit, a pointer to the source distributed vector for the equations matrices SOURCE. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesSource_DistributedVectorGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(sourceDistributedVector)) CALL FlagError("Source distributed vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(sourceVector)) CALL FlagError("Source vector is not associated.",err,error,*999)
#endif    

    sourceDistributedVector=>sourceVector%vector

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(sourceDistributedVector)) &
      & CALL FlagError("The source vector distributed vector is not associated.",err,error,*999)
#endif    
       
    EXITS("EquationsMatricesSource_DistributedVectorGet")
    RETURN
999 NULLIFY(sourceDistributedVector)
998 ERRORSEXITS("EquationsMatricesSource_DistributedVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesSource_DistributedVectorGet

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

    sourceVectors=>sourceVector%sourceVectors

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

  !>Gets a source vector from the source vectors.
  SUBROUTINE EquationsMatricesSources_SourceVectorGet(sourceVectors,sourceIdx,sourceVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors !<A pointer to the source vectors to get the source vector for
    INTEGER(INTG), INTENT(IN) :: sourcelIdx !<The source number of the source vector to get
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
      & CALL FlagError("Vector equations matrices has not been finished."
    
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

    sourceVectrs=>vectorMatrices%sourceVectors

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

END MODULE EquationsMatricesAccessRoutines

