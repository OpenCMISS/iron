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
  USE Kinds
  USE Strings
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC EquationsMatricesDynamic_DynamicMappingGet
  
  PUBLIC EquationsMatricesDynamic_EquationsMatrixGet
  
  PUBLIC EquationsMatricesDynamic_VectorMatricesGet

  PUBLIC EquationsMatricesLinear_EquationsMatrixGet
  
  PUBLIC EquationsMatricesLinear_LinearMappingGet

  PUBLIC EquationsMatricesLinear_VectorMatricesGet

  PUBLIC EquationsMatricesNonlinear_JacobianMatrixGet
  
  PUBLIC EquationsMatricesNonlinear_NonlinearMappingGet

  PUBLIC EquationsMatricesNonlinear_VectorMatricesGet

  PUBLIC EquationsMatricesScalar_EquationsScalarGet

  PUBLIC EquationsMatricesVector_DynamicMatricesGet
  
  PUBLIC EquationsMatricesVector_LinearMatricesGet
  
  PUBLIC EquationsMatricesVector_NonlinearMatricesGet
  
  PUBLIC EquationsMatricesVector_RHSVectorGet
  
  PUBLIC EquationsMatricesVector_SourceVectorGet

  PUBLIC EquationsMatricesVector_VectorEquationsGet

  PUBLIC EquationsMatricesVector_VectorMappingGet

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

    IF(ASSOCIATED(dynamicMapping)) CALL FlagError("Vector matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dynamicMatrices)) CALL FlagError("Dynamic matrices is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(dynamicMatrices%vectorMatrices)) &
      & CALL FlagError("Vector matrices is not associated for the dynamic matrices.",err,error,*999)
    IF(.NOT.ASSOCIATED(dynamicMatrices%vectorMatrices%vectorEquations)) &
      & CALL FlagError("Vector equations is not associated for the vector matrices.",err,error,*999)
    IF(.NOT.ASSOCIATED(dynamicMatrices%vectorMatrices%vectorEquations%vectorMapping)) &
      & CALL FlagError("Vector mapping is not associated for the vector equations.",err,error,*999)

    dynamicMapping=>dynamicMatrices%vectorMatrices%vectorEquations%vectorMapping%dynamicMapping
    IF(.NOT.ASSOCIATED(dynamicMapping)) &
      & CALL FlagError("Dynamic mapping is not associated for the dynamic matrices.",err,error,*999)
       
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
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsMatricesDynamic_EquationsMatrixGet",err,error,*998)

    IF(ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dynamicMatrices)) CALL FlagError("Dynamic matrices is not associated.",err,error,*999)
    IF(matrixIdx<1.OR.matrixIdx>dynamicMatrices%numberOfDynamicMatrices) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is invalid. The matrix index must be >= 1 and <= "// &
        & TRIM(NumberToVString(dynamicMatrices%numberOfDynamicMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(dynamicMatrices%matrices)) CALL FlagError("Dynamic matrices matrices is not allocated.",err,error,*999)
    
    equationsMatrix=>dynamicMatrices%matrices(matrixIdx)%ptr
    IF(.NOT.ASSOCIATED(equationsMatrix)) THEN
      localError="Equations matrix for dynamic matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
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

    IF(ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dynamicMatrices)) CALL FlagError("Dynamic matrices is not associated.",err,error,*999)

    vectorMatrices=>dynamicMatrices%vectorMatrices
    IF(.NOT.ASSOCIATED(vectorMatrices)) &
      & CALL FlagError("Vector matrices is not associated for the dynamic matrices.",err,error,*999)
       
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
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsMatricesLinear_EquationsMatrixGet",err,error,*998)

    IF(ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(linearMatrices)) CALL FlagError("Linear matrices is not associated.",err,error,*999)
    IF(matrixIdx<1.OR.matrixIdx>linearMatrices%numberOfLinearMatrices) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is invalid. The matrix index must be >= 1 and <= "// &
        & TRIM(NumberToVString(linearMatrices%numberOfLinearMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(linearMatrices%matrices)) CALL FlagError("Linear matrices matrices is not allocated.",err,error,*999)
    
    equationsMatrix=>linearMatrices%matrices(matrixIdx)%ptr
    IF(.NOT.ASSOCIATED(equationsMatrix)) THEN
      localError="Equations matrix for linear matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
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

    IF(ASSOCIATED(linearMapping)) CALL FlagError("Vector matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(linearMatrices)) CALL FlagError("Linear matrices is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(linearMatrices%vectorMatrices)) &
      & CALL FlagError("Vector matrices is not associated for the linear matrices.",err,error,*999)
    IF(.NOT.ASSOCIATED(linearMatrices%vectorMatrices%vectorEquations)) &
      & CALL FlagError("Vector equations is not associated for the vector matrices.",err,error,*999)
    IF(.NOT.ASSOCIATED(linearMatrices%vectorMatrices%vectorEquations%vectorMapping)) &
      & CALL FlagError("Vector mapping is not associated for the vector equations.",err,error,*999)

    linearMapping=>linearMatrices%vectorMatrices%vectorEquations%vectorMapping%linearMapping
    IF(.NOT.ASSOCIATED(linearMapping)) &
      & CALL FlagError("Linear mapping is not associated for the linear matrices.",err,error,*999)
       
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

    IF(ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(linearMatrices)) CALL FlagError("Linear matrices is not associated.",err,error,*999)

    vectorMatrices=>linearMatrices%vectorMatrices
    IF(.NOT.ASSOCIATED(vectorMatrices)) &
      & CALL FlagError("Vector matrices is not associated for the linear matrices.",err,error,*999)
       
    EXITS("EquationsMatricesLinear_VectorMatricesGet")
    RETURN
999 NULLIFY(vectorMatrices)
998 ERRORSEXITS("EquationsMatricesLinear_VectorMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesLinear_VectorMatricesGet

  !
  !================================================================================================================================
  !

  !>Gets the specified Jacobian matrix for nonlinear matrices.
  SUBROUTINE EquationsMatricesNonlinear_JacobianMatrixGet(nonlinearMatrices,matrixIdx,JacobianMatrix,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices !<A pointer to the nonlinear matrices to get the Jacobian matrix for
    INTEGER(INTG), INTENT(IN) :: matrixIdx !<The matrix index of the Jacobian matrix to get
    TYPE(EquationsJacobianType), POINTER :: JacobianMatrix !<On exit, a pointer to the Jacobian matrix for the matrixIdx'th nonlinear matrix. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsMatricesNonlinear_JacobianMatrixGet",err,error,*998)

    IF(ASSOCIATED(JacobianMatrix)) CALL FlagError("Jacobian matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(nonlinearMatrices)) CALL FlagError("Nonlinear matrices is not associated.",err,error,*999)
    IF(matrixIdx<1.OR.matrixIdx>nonlinearMatrices%numberOfJacobians) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is invalid. The matrix index must be >= 1 and <= "// &
        & TRIM(NumberToVString(nonlinearMatrices%numberOfJacobians,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(nonlinearMatrices%jacobians)) CALL FlagError("Nonlinear matrices Jacobians is not allocated.",err,error,*999)
    
    jacobianMatrix=>nonlinearMatrices%jacobians(matrixIdx)%ptr
    IF(.NOT.ASSOCIATED(jacobianMatrix)) THEN
      localError="Jacobian matrix for nonlinear matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("EquationsMatricesNonlinear_JacobianMatrixGet")
    RETURN
999 NULLIFY(jacobianMatrix)
998 ERRORS("EquationsMatricesNonlinear_JacobianMatrixGet",err,error)
    EXITS("EquationsMatricesNonlinear_JacobianMatrixGet")
    RETURN 1
    
  END SUBROUTINE EquationsMatricesNonlinear_JacobianMatrixGet

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

    IF(ASSOCIATED(nonlinearMapping)) CALL FlagError("Vector matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(nonlinearMatrices)) CALL FlagError("Nonlinear matrices is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(nonlinearMatrices%vectorMatrices)) &
      & CALL FlagError("Vector matrices is not associated for the nonlinear matrices.",err,error,*999)
    IF(.NOT.ASSOCIATED(nonlinearMatrices%vectorMatrices%vectorEquations)) &
      & CALL FlagError("Vector equations is not associated for the vector matrices.",err,error,*999)
    IF(.NOT.ASSOCIATED(nonlinearMatrices%vectorMatrices%vectorEquations%vectorMapping)) &
      & CALL FlagError("Vector mapping is not associated for the vector equations.",err,error,*999)

    nonlinearMapping=>nonlinearMatrices%vectorMatrices%vectorEquations%vectorMapping%nonlinearMapping
    IF(.NOT.ASSOCIATED(nonlinearMapping)) &
      & CALL FlagError("Nonlinear mapping is not associated for the nonlinear matrices.",err,error,*999)
       
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

  !>Gets the vector matrices for nonlinear matrices.
  SUBROUTINE EquationsMatricesNonlinear_VectorMatricesGet(nonlinearMatrices,vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices !<A pointer to the nonlinear matrices to get the vector matrices for
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<On exit, a pointer to the vector matrices for the nonlinear matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesNonlinear_VectorMatricesGet",err,error,*998)

    IF(ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(nonlinearMatrices)) CALL FlagError("Nonlinear matrices is not associated.",err,error,*999)

    vectorMatrices=>nonlinearMatrices%vectorMatrices
    IF(.NOT.ASSOCIATED(vectorMatrices)) &
      & CALL FlagError("Vector matrices is not associated for the nonlinear matrices.",err,error,*999)
       
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

  !>Gets the scalar equations for an scalar equations matrices.
  SUBROUTINE EquationsMatricesScalar_EquationsScalarGet(scalarMatrices,scalarEquations,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesScalarType), POINTER :: scalarMatrices !<A pointer to the equations scalar matrices to get the scalar equations for
    TYPE(EquationsScalarType), POINTER :: scalarEquations !<On exit, a pointer to the scalar equations in the specified scalar equations matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesScalar_EquationsScalarGet",err,error,*998)

    IF(ASSOCIATED(scalarEquations)) CALL FlagError("Scalar equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(scalarMatrices)) CALL FlagError("Scalar matrices is not associated.",err,error,*999)

    scalarEquations=>scalarMatrices%scalarEquations
    IF(.NOT.ASSOCIATED(scalarEquations)) &
      & CALL FlagError("Scalar equations is not associated for the scalar matrices.",err,error,*999)
       
    EXITS("EquationsMatricesScalar_EquationsScalarGet")
    RETURN
999 NULLIFY(scalarEquations)
998 ERRORSEXITS("EquationsMatricesScalar_EquationsScalarGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesScalar_EquationsScalarGet

  !
  !================================================================================================================================
  !

  !>Gets the dynamic vector matrices for vector matrices.
  SUBROUTINE EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations vector matrices to get the dynamic matrices for
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices !<On exit, a pointer to the dynamic Matrices in the specified vector equations matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_DynamicMatricesGet",err,error,*998)

    IF(ASSOCIATED(dynamicMatrices)) CALL FlagError("Dynamic matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)

    dynamicMatrices=>vectorMatrices%dynamicMatrices
    IF(.NOT.ASSOCIATED(dynamicMatrices)) CALL FlagError("Dynamic matrices is not associated for the vector matrices.", &
      & err,error,*999)
       
    EXITS("EquationsMatricesVector_DynamicMatricesGet")
    RETURN
999 NULLIFY(dynamicMatrices)
998 ERRORSEXITS("EquationsMatricesVector_DynamicMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_DynamicMatricesGet

  !
  !================================================================================================================================
  !

  !>Gets the linear vector matrices for an vector matrices.
  SUBROUTINE EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations vector matrices to get the linear matrices for
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices !<On exit, a pointer to the linear Matrices in the specified vector equations matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_LinearMatricesGet",err,error,*998)

    IF(ASSOCIATED(linearMatrices)) CALL FlagError("Linear matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)

    linearMatrices=>vectorMatrices%linearMatrices
    IF(.NOT.ASSOCIATED(linearMatrices)) CALL FlagError("Linear matrices is not associated for the vector matrices.",err,error,*999)
       
    EXITS("EquationsMatricesVector_LinearMatricesGet")
    RETURN
999 NULLIFY(linearMatrices)
998 ERRORSEXITS("EquationsMatricesVector_LinearMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_LinearMatricesGet

  !
  !================================================================================================================================
  !

  !>Gets the nonlinear vector matrices for an vector matrices.
  SUBROUTINE EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations vector matrices to get the nonlinear matrices for
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices !<On exit, a pointer to the nonlinear Matrices in the specified vector equations matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_NonlinearMatricesGet",err,error,*998)

    IF(ASSOCIATED(nonlinearMatrices)) CALL FlagError("Nonlinear matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)

    nonlinearMatrices=>vectorMatrices%nonlinearMatrices
    IF(.NOT.ASSOCIATED(nonlinearMatrices))  &
      & CALL FlagError("Nonlinear matrices is not associated for the vector matrices.",err,error,*999)
       
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

  !>Gets the rhs vector for an vector matrices.
  SUBROUTINE EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations vector matrices to get the rhs matrices for
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector !<On exit, a pointer to the RHS vector in the specified vector equations matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_RHSVectorGet",err,error,*998)

    IF(ASSOCIATED(rhsVector)) CALL FlagError("RHS vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)

    rhsVector=>vectorMatrices%rhsVector
    IF(.NOT.ASSOCIATED(rhsVector)) CALL FlagError("RHS vector is not associated for the vector matrices.",err,error,*999)
       
    EXITS("EquationsMatricesVector_RHSVectorGet")
    RETURN
999 NULLIFY(rhsVector)
998 ERRORSEXITS("EquationsMatricesVector_RHSVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_RHSVectorGet

  !
  !================================================================================================================================
  !

  !>Gets the source vector for an vector matrices.
  SUBROUTINE EquationsMatricesVector_SourceVectorGet(vectorMatrices,sourceVector,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations vector matrices to get the source vector for
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector !<On exit, a pointer to the source vector in the specified vector equations matrices. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMatricesVector_SourceVectorGet",err,error,*998)

    IF(ASSOCIATED(sourceVector)) CALL FlagError("Source vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)

    sourceVector=>vectorMatrices%sourceVector
    IF(.NOT.ASSOCIATED(sourceVector)) CALL FlagError("Source vector is not associated for the vector matrices.",err,error,*999)
       
    EXITS("EquationsMatricesVector_SourceVectorGet")
    RETURN
999 NULLIFY(sourceVector)
998 ERRORSEXITS("EquationsMatricesVector_SourceVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_SourceVectorGet

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

    IF(ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)

    vectorEquations=>vectorMatrices%vectorEquations
    IF(.NOT.ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is not associated for the vector matrices.", &
      & err,error,*999)
       
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

    IF(ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(vectorMatrices%vectorEquations)) &
      & CALL FlagError("Vector matrices vector equations is not associated.",err,error,*999)
    
    vectorMapping=>vectorMatrices%vectorEquations%vectorMapping
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated for the vector matrices.",err,error,*999)
       
    EXITS("EquationsMatricesVector_VectorMappingGet")
    RETURN
999 NULLIFY(vectorMapping)
998 ERRORSEXITS("EquationsMatricesVector_VectorMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatricesVector_VectorMappingGet

  !
  !================================================================================================================================
  !

END MODULE EquationsMatricesAccessRoutines

