!> \file
!> \author Ishani Roy & Sander Land
!> \brief This module handles all Monodomain equations routines which use Strang Splitting with hard-coded cell models.
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
!> Contributor(s): Ishani Roy, Sander Land, Chris Bradley
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

!>This module handles all Monodomain equations routines.
MODULE MonodomainEquationsRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BoundaryConditionsRoutines
  USE Constants
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE DecompositionAccessRoutines
  USE DistributedMatrixVector
  USE DistributedMatrixVectorAccessRoutines
  USE DomainMappings
  USE ElectrophysiologyCellRoutines
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingRoutines
  USE EquationsMappingAccessRoutines
  USE EquationsMatricesRoutines
  USE EquationsMatricesAccessRoutines
  USE EquationsSetAccessRoutines
  USE FIELD_IO_ROUTINES
  USE FieldRoutines
  USE FieldAccessRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE MatrixVector
  USE ProblemAccessRoutines
  USE RegionAccessRoutines
  USE Strings
  USE SolverRoutines
  USE SolverAccessRoutines
  USE SolverMappingAccessRoutines
  USE Timer
  USE Types

#include "macros.h"  


  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC Monodomain_PostLoop
 
  PUBLIC Monodomain_EquationsSetSetup
  
  PUBLIC Monodomain_FiniteElementCalculate

  PUBLIC Monodomain_EquationsSetSolutionMethodSet
  
  PUBLIC Monodomain_EquationsSetSpecificationSet

  PUBLIC Monodomain_ProblemSpecificationSet
  
  PUBLIC Monodomain_ProblemSetup
  
  PUBLIC Monodomain_PreSolve

  PUBLIC Monodomain_PostSolve
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Runs after each control loop iteration
  SUBROUTINE Monodomain_PostLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: currentIteration,equationsSetIdx,inputIteration,loopType,numberOfEquationsSets,outputIteration,outputType, &
      & regionUserNumber
    REAL(DP) :: currentTime,startTime,stopTime,timeIncrement
    TYPE(ControlLoopType), POINTER :: parentLoop
    TYPE(ControlLoopTimeType), POINTER :: timeLoop,parentTimeLoop
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldsType), POINTER :: regionFields
    TYPE(ProblemType), POINTER :: problem
    TYPE(RegionType), POINTER :: region   
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: filename,localError,method

    ENTERS("Monodomain_PostLoop",err,error,*999)

    CALL ControlLoop_OutputTypeGet(controlLoop,outputType,err,error,*999)
    IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
      CALL ControlLoop_TypeGet(controlLoop,loopType,err,error,*999)
      SELECT CASE(loopType)
      CASE(CONTROL_SIMPLE_TYPE)
        !Do nothing
      CASE(CONTROL_FIXED_LOOP_TYPE)
        !Do nothing
      CASE(CONTROL_TIME_LOOP_TYPE)
        !Export the dependent field for this time step
        CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
          & outputIteration,inputIteration,err,error,*999)
        NULLIFY(problem)
        CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
        !Get the solver
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)            
        NULLIFY(solver)
        CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
        !Loop over the equations sets associated with the solver
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
        DO equationsSetIdx=1,numberOfEquationsSets
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
          NULLIFY(region)
          CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
          CALL Region_UserNumberGet(region,regionUserNumber,err,error,*999)
          NULLIFY(dependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
          filename="Time_"//TRIM(NumberToVString(regionUserNumber,"*",err,error))// &
            & "_"//TRIM(NumberToVString(currentIteration,"*",err,error))           
          method="FORTRAN"
          IF(outputIteration/=0) THEN
            IF(MOD(currentIteration,outputIteration)==0) THEN
              NULLIFY(regionFields)
              CALL Region_FieldsGet(region,regionFields,err,error,*999)
              CALL FIELD_IO_NODES_EXPORT(regionFields,filename,method,err,error,*999)
            ENDIF
          ENDIF
        ENDDO !equationsSetIdx
      CASE(CONTROL_WHILE_LOOP_TYPE)
        !Do nothing
      CASE(CONTROL_LOAD_INCREMENT_LOOP_TYPE)
        !Do nothing
      CASE DEFAULT
        localError="The control loop type of "//TRIM(NumberToVString(loopType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
     
    EXITS("Monodomain_PostLoop")
    RETURN
999 ERRORSEXITS("Monodomain_PostLoop",err,error)
    RETURN 1
    
  END SUBROUTINE Monodomain_PostLoop
  
  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a monodomain equation set class.
  SUBROUTINE Monodomain_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSubtype,esType
    TYPE(VARYING_STRING) :: localError

    ENTERS("Monodomain_EquationsSetSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(SIZE(specification,1)<3) THEN
      localError="The size of the specified specification array of "// &
        & TRIM(NumberToVString(SIZE(specification,1),"*",err,error))//" is invalid. The size should be >= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    esType=specification(2)
    esSubtype=specification(3)
    SELECT CASE(esType)
    CASE(EQUATIONS_SET_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
      SELECT CASE(esSubtype)
      CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE, &
        & EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
        !ok
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVstring(esSubtype,"*",err,error))// &
          & " is not valid for a Monodomain equation type of a Strang splitting equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="Equations set equation type "//TRIM(NumberToVstring(esType,"*",err,error))// &
        & " is not valid for a monodomain equations set class."
    END SELECT
    !Set full specification
    IF(ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is already allocated.",err,error,*999)
    ALLOCATE(equationsSet%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
    equationsSet%specification(1:3)=[EQUATIONS_SET_BIOELECTRICS_CLASS,EQUATIONS_SET_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE, &
      & esSubtype]

    EXITS("Monodomain_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("Monodomain_EquationsSetSpecificationSet",err,error)
    EXITS("Monodomain_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Monodomain_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a monodomain problem class.
  SUBROUTINE Monodomain_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemType,problemSubtype

    ENTERS("Monodomain_ProblemSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(ALLOCATED(problem%specification)) CALL FlagError("Problem specification is already allocated.",err,error,*999)
    IF(SIZE(problemSpecification,1)<3) THEN
      localError="The size of the specified problem specification array of "// &
        & TRIM(NumberToVString(SIZE(problemSpecification,1),"*",err,error))// &
        & " is invalid. The size should be >= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    problemType=problemSpecification(2)
    SELECT CASE(problemType)
    CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
      problemSubtype=problemSpecification(3)
      SELECT CASE(problemSubtype)
      CASE(PROBLEM_MONODOMAIN_BUENOOROVIO_SUBTYPE, &
        & PROBLEM_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
        !ok
      CASE DEFAULT
        localError="Problem subtype "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
          & " is not valid for a Monodomain equation type of a Strang splitting  problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="Problem equation type "//TRIM(NumberToVstring(problemType,"*",err,error))// &
        & " is not valid for a monodomain problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    ALLOCATE(problem%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
    problem%specification(1:3)=[PROBLEM_BIOELECTRICS_CLASS,PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE,problemSubtype]
 
    EXITS("Monodomain_ProblemSpecificationSet")
    RETURN
999 ERRORSEXITS("Monodomain_ProblemSpecificationSet",err,error)
    RETURN 1
    
  END SUBROUTINE Monodomain_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a Monodomain equation finite element equations set.
  SUBROUTINE Monodomain_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) colsVariableType,columnComponentIdx,columnElementDOFIdx,columnElementParameterIdx,columnXiIdx,componentIdx, &
      & componentIdx2,esSpecification(3),gaussPointIdx,numberOfColumnElementParameters,numberOfColsComponents,numberOfDimensions, &
      & numberOfGauss,numberOfRowElementParameters,numberOfRowsComponents,numberOfXi,rowComponentIdx,rowElementDOFIdx, &
      & rowElementParameterIdx,rowsVariableType,rowXiIdx,scalingType,solutionMethod,xiIdx
    REAL(DP) :: columnPhi,columndPhidX(3),columndPhidXi(3),D(3,3),Df,Dt,dXidX(3,3),f(3),fnorm,gaussWeight,jacobian, &
      & jacobianGaussWeight,rowPhi,rowdPhidXi(3),rowdPhidX(3),sum
    LOGICAL :: update,updateDamping,updateMatrices,updateRHS,updateStiffness
    TYPE(BasisType), POINTER :: columnBasis,dependentBasis,geometricBasis,rowBasis
    TYPE(DecompositionType), POINTER :: dependentDecomposition,geometricDecomposition
    TYPE(DomainType), POINTER :: columnDomain,dependentDomain,geometricDomain,rowDomain
    TYPE(DomainElementsType), POINTER :: columnDomainElements,dependentDomainElements,geometricDomainElements,rowDomainElements
    TYPE(DomainTopologyType), POINTER :: columnDomainTopology,dependentDomainTopology,geometricDomainTopology,rowDomainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatrixType), POINTER :: dampingMatrix,stiffnessMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,fibreField,geometricField,independentField,materialsField
    TYPE(FieldInterpolationParametersType), POINTER :: colsInterpParameters,fibreInterpParameters,geometricInterpParameters, &
      & materialsInterpParameters,rowsInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: fibreInterpPoint,geometricInterpPoint,materialsInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: colsVariable,fibreVariable,geometricVariable,materialsVariable,rowsVariable
    TYPE(QuadratureSchemeType), POINTER :: columnQuadratureScheme,dependentQuadratureScheme,geometricQuadratureScheme, &
      & rowQuadratureScheme
    TYPE(VARYING_STRING) :: localError
     
    ENTERS("Monodomain_FiniteElementCalculate",err,error,*999)
    
    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE, &
      & EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Monodomain equation type of a Strang splitting equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(lhsMapping)
    CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
    NULLIFY(rowsVariable)
    CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(dynamicMapping)
    CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
    NULLIFY(colsVariable)
    CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,colsVariable,err,error,*999)
    NULLIFY(dynamicMatrices)
    CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
    NULLIFY(stiffnessMatrix)
    CALL EquationsMatricesDynamic_StiffnessMatrixGet(dynamicMatrices,stiffnessMatrix,err,error,*999)
    CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
    NULLIFY(dampingMatrix)
    CALL EquationsMatricesDynamic_DampingMatrixGet(dynamicMatrices,dampingMatrix,err,error,*999)      
    CALL EquationsMatrix_UpdateMatrixGet(dampingMatrix,updateDamping,err,error,*999)
    NULLIFY(rhsMapping)
    CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
    NULLIFY(rhsVector)
    updateRHS=.FALSE.
    IF(ASSOCIATED(rhsMapping)) THEN
      CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
      CALL EquationsMatricesRHS_UpdateVectorGet(rhsVector,updateRHS,err,error,*999)
    ENDIF

    updateMatrices=(updateStiffness.OR.updateDamping)
    update=(updateMatrices.OR.updateRHS)

    IF(update) THEN

      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)

       NULLIFY(geometricVariable)
      CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
      NULLIFY(geometricDecomposition)
      CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
      NULLIFY(geometricDomain)
      CALL Decomposition_DomainGet(geometricDecomposition,0,geometricDomain,err,error,*999)
      NULLIFY(geometricDomainTopology)
      CALL Domain_DomainTopologyGet(geometricDomain,geometricDomainTopology,err,error,*999)
      NULLIFY(geometricDomainElements)
      CALL DomainTopology_DomainElementsGet(geometricDomainTopology,geometricDomainElements,err,error,*999)
      NULLIFY(geometricBasis)
      CALL DomainElements_ElementBasisGet(geometricDomainElements,elementNumber,geometricBasis,err,error,*999)
      CALL Basis_NumberOfXiGet(geometricBasis,numberOfXi,err,error,*999)
      NULLIFY(geometricQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(geometricBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,geometricQuadratureScheme,err,error,*999)
      NULLIFY(geometricInterpParameters)
      CALL EquationsInterpolation_GeometricParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpParameters,err,error,*999)
      NULLIFY(geometricInterpPoint)
      CALL EquationsInterpolation_GeometricPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpPoint, &
        & err,error,*999)
      NULLIFY(geometricInterpPointMetrics)
      CALL EquationsInterpolation_GeometricPointMetricsGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpPointMetrics,err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,geometricInterpParameters,err,error,*999)

      NULLIFY(fibreField)
      CALL EquationsSet_FibreFieldExists(equationsSet,fibreField,err,error,*999)
      NULLIFY(fibreVariable)
      NULLIFY(fibreInterpParameters)
      NULLIFY(fibreInterpPoint)
      IF(ASSOCIATED(fibreField)) THEN
        CALL Field_VariableGet(fibreField,FIELD_U_VARIABLE_TYPE,fibreVariable,err,error,*999)
        CALL EquationsInterpolation_FibreParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,fibreInterpParameters, &
          & err,error,*999)
        CALL EquationsInterpolation_FibrePointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,fibreInterpPoint,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,fibreInterpParameters,err,error,*999)
      ENDIF

      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(dependentDecomposition)
      CALL Field_DecompositionGet(dependentField,dependentDecomposition,err,error,*999)
      NULLIFY(dependentDomain)
      CALL Decomposition_DomainGet(dependentDecomposition,0,dependentDomain,err,error,*999)
      NULLIFY(dependentDomainTopology)
      CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
      NULLIFY(dependentDomainElements)
      CALL DomainTopology_DomainElementsGet(dependentDomainTopology,dependentDomainElements,err,error,*999)
      NULLIFY(dependentBasis)
      CALL DomainElements_ElementBasisGet(dependentDomainElements,elementNumber,dependentBasis,err,error,*999)
      NULLIFY(dependentQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(dependentBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,dependentQuadratureScheme,err,error,*999)
      CALL BasisQuadratureScheme_NumberOfGaussGet(dependentQuadratureScheme,numberOfGauss,err,error,*999)
    
      NULLIFY(materialsField)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
      NULLIFY(materialsVariable)
      CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,materialsVariable,err,error,*999)
      CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & materialsInterpParameters,err,error,*999)
      CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,materialsInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsInterpParameters,err,error,*999)
  
      CALL FieldVariable_VariableTypeGet(rowsVariable,rowsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(rowsVariable,numberOfRowsComponents,err,error,*999)

      CALL FieldVariable_VariableTypeGet(colsVariable,colsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(colsVariable,numberOfColsComponents,err,error,*999)     

      DO gaussPointIdx=1,numberOfGauss
        !Get interpolated geometric and material interpolated point
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
        CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,materialsInterpPoint, &
          & err,error,*999)
        IF(ASSOCIATED(fibreVariable)) &
          & CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,fibreInterpPoint, &
          & err,error,*999)

        !Diffusion tensor  D  =  Dt I + (Df - Dt) f f^T  where Dt and Df are diffusivity/conductivity in fiber/transverse directions
        Df = materialsInterpPoint%values(2,NO_PART_DERIV) ! 2 = Df
        Dt = materialsInterpPoint%values(3,NO_PART_DERIV) ! 3 = Dt
        fnorm = 0.0
        DO componentIdx=1,numberOfDimensions
          f(componentIdx) = materialsInterpPoint%values(3+componentIdx,NO_PART_DERIV) ! 4,5[,6] = f
          fnorm = fnorm + f(componentIdx)*f(componentIdx)
        ENDDO !componentIdx
        !Normalize f, and fill in default for 0,0,0 -> 1,0,0
        fnorm = SQRT(fnorm)
        IF(fnorm < 1e-6) THEN
          f = [ 1.0, 0.0, 0.0 ] ! default
        ELSE
          f = f / fnorm
        ENDIF
        DO componentIdx=1,numberOfDimensions
          D(componentIdx,:)  = 0.0
          D(componentIdx,componentIdx) = Dt
          DO componentIdx2=1,numberOfDimensions
            D(componentIdx,componentIdx2) = D(componentIdx,componentIdx2) + (Df - Dt) * f(componentIdx) * f(componentIdx2)
          ENDDO !componentIdx2
        ENDDO !componentIdx
            
        !Calculate weight = det J * gauss pt weight
        CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
        CALL BasisQuadratureScheme_GaussWeightGet(geometricQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
        jacobianGaussWeight=jacobian*gaussWeight

        DO columnComponentIdx=1,numberOfColsComponents
          DO rowXiIdx=1,numberOfXi
            dXidX(rowXiIdx,columnComponentIdx)=geometricInterpPointMetrics%dXidX(rowXiIdx,columnComponentIdx)
          ENDDO !rowXiIdx
        ENDDO !columnXiIdx
       
        !Loop over field components
        rowElementDOFIdx=0          
        DO rowComponentIdx=1,numberOfRowsComponents
          NULLIFY(rowDomain)
          CALL FieldVariable_ComponentDomainGet(rowsVariable,rowComponentIdx,rowDomain,err,error,*999)
          NULLIFY(rowDomainTopology)
          CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
          NULLIFY(rowDomainElements)
          CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
          NULLIFY(rowBasis)
          CALL DomainElements_ElementBasisGet(rowDomainElements,elementNumber,rowBasis,err,error,*999)
          CALL Basis_QuadratureSchemeGet(rowBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,rowQuadratureScheme,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(rowBasis,numberOfRowElementParameters,err,error,*999)
          !Loop over element rows
          DO rowElementParameterIdx=1,numberOfRowElementParameters
            rowElementDOFIdx=rowElementDOFIdx+1
            CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
              & gaussPointIdx,rowPhi,err,error,*999)
            IF(updateMatrices) THEN
              rowdPhidX=0.0_DP
              DO rowXiIdx=1,numberOfXi
                CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx, &
                  & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(rowXiIdx),gaussPointIdx,rowdPhidXi(rowXiIdx), &
                  & err,error,*999)
                DO componentIdx=1,numberOfDimensions
                  rowdPhidX(componentIdx)=rowdPhidX(componentIdx)+rowdPhidXi(columnXiIdx)*dXidX(rowXiIdx,componentIdx)
                ENDDO !componentIdx
              ENDDO !columnXiiIdx
              columnElementDOFIdx=0
              !Loop over element columns. TODO: use symmetry?
              DO columnComponentIdx=1,numberOfColsComponents
                NULLIFY(columnDomain)
                CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,columnDomain,err,error,*999)
                NULLIFY(columnDomainTopology)
                CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
                NULLIFY(columnDomainElements)
                CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
                NULLIFY(columnBasis)
                CALL DomainElements_ElementBasisGet(columnDomainElements,elementNumber,columnBasis,err,error,*999)
                CALL Basis_QuadratureSchemeGet(columnBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,columnQuadratureScheme, &
                  & err,error,*999)
                CALL Basis_NumberOfElementParametersGet(columnBasis,numberOfColumnElementParameters,err,error,*999)
                DO columnElementParameterIdx=1,numberOfColumnElementParameters
                  columnElementDOFIdx=columnElementDOFIdx+1
                  CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme,columnElementParameterIdx, &
                    & NO_PART_DERIV,gaussPointIdx,columnPhi,err,error,*999)
                  IF(updateStiffness) THEN
                    columndPhidX=0.0_DP
                    DO columnXiIdx=1,numberOfXi
                      CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme,columnElementParameterIdx, &
                        & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(columnXiIdx),gaussPointIdx,columndPhidXi(columnXiIdx), &
                        & err,error,*999)
                      DO componentIdx=1,numberOfDimensions
                        columndPhidX(componentIdx)=columndPhidX(componentIdx)+ &
                          & columndPhidXi(columnXiIdx)*dXidX(columnXiIdx,componentIdx)
                      ENDDO !componentIdx
                    ENDDO !columnXiiIdx
                    sum=0.0_DP
                    DO componentIdx=1,numberOfDimensions
                      DO componentIdx2=1,numberOfDimensions
                        sum=sum+D(componentIdx,componentIdx2)*rowdPhidX(componentIdx)*columndPhidX(componentIdx2)
                      ENDDO !componentIdx2
                    ENDDO !componentIdx
                    stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                      & sum*jacobianGaussWeight ! Aij = int D_ij * dphi_m/dx_i * dphi_n/dx_j 
                  ENDIF !update stiffness        
                  IF(updateDamping) THEN 
                    dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                      & rowPhi*columnPhi*jacobianGaussWeight
                  ENDIF !update damping
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDIF !update matrices
            IF(updateRHS) rhsVector%elementVector%vector(rowElementDOFIdx)=0.0_DP
          ENDDO !rowElementParameterIdx
        ENDDO !rowComponentIdx
      ENDDO !gaussPointIdx

      CALL Field_ScalingTypeGet(dependentField,scalingType,err,error,*999)
      IF(scalingType/=FIELD_NO_SCALING) THEN
        NULLIFY(rowsInterpParameters)
        CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,rowsVariableType,rowsInterpParameters, &
          & err,error,*999)
        NULLIFY(colsInterpParameters)
        CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,colsVariableType,colsInterpParameters, &
          & err,error,*999)
        CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,rowsInterpParameters,err,error,*999)
        CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,colsInterpParameters,err,error,*999)
        !Loop over element rows
        rowElementDOFIdx=0          
        DO rowComponentIdx=1,numberOfRowsComponents
          NULLIFY(rowDomain)
          CALL FieldVariable_ComponentDomainGet(rowsVariable,rowComponentIdx,rowDomain,err,error,*999)
          NULLIFY(rowDomainTopology)
          CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
          NULLIFY(rowDomainElements)
          CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
          NULLIFY(rowBasis)
          CALL DomainElements_ElementBasisGet(rowDomainElements,elementNumber,rowBasis,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(rowBasis,numberOfRowElementParameters,err,error,*999)
          DO rowElementParameterIdx=1,numberOfRowElementParameters
            rowElementDOFIdx=rowElementDOFIdx+1                    
            IF(updateMatrices) THEN
              !Loop over element columns
              columnElementDOFIdx=0
              DO columnComponentIdx=1,numberOfColsComponents
                NULLIFY(columnDomain)
                CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,columnDomain,err,error,*999)
                NULLIFY(columnDomainTopology)
                CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
                NULLIFY(columnDomainElements)
                CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
                NULLIFY(columnBasis)
                CALL DomainElements_ElementBasisGet(columnDomainElements,elementNumber,columnBasis,err,error,*999)
                CALL Basis_NumberOfElementParametersGet(columnBasis,numberOfColumnElementParameters,err,error,*999)
                DO columnElementParameterIdx=1,numberOfColumnElementParameters
                  columnElementDOFIdx=columnElementDOFIdx+1
                  IF(updateStiffness) THEN
                    stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                      & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)* &
                      & colsInterpParameters%scaleFactors(columnElementParameterIdx,columnComponentIdx)
                  ENDIF !update stiffness
                  IF(updateDamping) THEN
                    dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                      & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)* &
                      & colsInterpParameters%scaleFactors(columnElementParameterIdx,columnComponentIdx)
                  ENDIF !update damping
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDIF !update matrices
            IF(updateRHS) THEN
              rhsVector%elementVector%vector(rowElementDOFIdx)=rhsVector%elementVector%vector(rowElementDOFIdx)* &
                & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
            ENDIF !update RHS
          ENDDO !rowElementParameterIdx
        ENDDO !rowComponentIdx
      ENDIF !scaling
    ENDIF !update

    EXITS("Monodomain_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("Monodomain_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE Monodomain_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Sets up the Monodomain equation type of a Strang splitting equations set class.
  SUBROUTINE Monodomain_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup a Monodomain equation on.
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(3)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Monodomain_EquationsSetSetup",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE, &
      & EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
      CALL Monodomain_EquationsSetSubtypeSetUP(equationsSet,equationsSetSetup,err,error,*999)        
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Monodomain equation type of a Strang splitting equation set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Monodomain_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("Monodomain_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Monodomain_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a Monodomain equation type of an Strang splitting equations set class.
  SUBROUTINE Monodomain_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(3)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Monodomain_EquationsSetSolutionMethodSet",err,error,*999)
    
    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE, &
      & EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)      
      SELECT CASE(solutionMethod)
      CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
        equationsSet%solutionMethod=EQUATIONS_SET_FEM_SOLUTION_METHOD
      CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The specified solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="Equations set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Monodomain equation type of monodomain Strang splitting equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Monodomain_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("Monodomain_EquationsSetSolutionMethodSet",err,error)
    EXITS("Monodomain_EquationsSetSolutionMethodSet")
    RETURN 1

  END SUBROUTINE Monodomain_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets up the Monodomain equation.
  SUBROUTINE Monodomain_EquationsSetSubtypeSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,esSpecification(3),equationsSetSubtype,geometricMeshComponent,geometricScalingType, &
      & numberOfDimensions,numberOfMaterialsComponents,numberOfComponents,solutionMethod,sparsityType
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsSetMaterialsType), POINTER :: equationsMaterials
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: geometricField
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Monodomain_EquationsSetSubtypeSetup",err,error,*999)
 
    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    equationsSetSubtype=esSpecification(3)
    SELECT CASE(equationsSetSubtype)
    CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE, &
      & EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="The equations set subtype of "//TRIM(NumberToVString(equationsSetSubtype,"*",err,error))// &
        & " does not equal a Monodomain equation subtype."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    NULLIFY(region)
    CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
    
    SELECT CASE(equationsSetSetup%setupType)
    CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
      !-----------------------------------------------------------------
      ! I n i t i a l   s e t u p
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL Monodomain_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        !Do nothing
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a Monodomain equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
      !-----------------------------------------------------------------
      ! G e o m e t r i c   f i e l d
      !-----------------------------------------------------------------
      !Do nothing
    CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
      !-----------------------------------------------------------------
      ! D e p e n d e n t   f i e l d
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
          NULLIFY(geometricField)
          CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
         !Create the auto created dependent field
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSet%dependent%dependentField,err,error,*999)
          CALL Field_TypeSetAndLock(equationsSet%dependent%dependentField,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_DependentTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DEPENDENT_TYPE,err,error,*999)
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          CALL Field_DecompositionSetAndLock(equationsSet%dependent%dependentField,geometricDecomposition,err,error,*999)
          CALL Field_GeometricFieldSetAndLock(equationsSet%dependent%dependentField,geometricField,err,error,*999)
          CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%dependentField,2,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsSet%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
            & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & FIELD_DP_TYPE,err,error,*999)           
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1,&
            & err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
            & err,error,*999)          
          !Default to the geometric interpolation setup
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
          CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
            & geometricMeshComponent,err,error,*999)
          CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
            & geometricMeshComponent,err,error,*999)
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
              & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
              & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            !Default the scaling to the geometric field scaling
            CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
            CALL Field_ScalingTypeSet(equationsSet%dependent%dependentField,geometricScalingType,err,error,*999)
          CASE DEFAULT
            localError="The solution method of "//TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))// &
              & " is invalid or not implemented"
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          !User specified field
          CALL FlagError("No user specified field supported!",err,error,*999)
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
          CALL Field_CreateFinish(equationsSet%dependent%dependentField,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a Monodomain equation"
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
      !-----------------------------------------------------------------
      ! I n d e p e n d e n t   f i e l d
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        IF(equationsSet%INDEPENDENT%independentFieldAutoCreated) THEN
          NULLIFY(geometricField)
          CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
          SELECT CASE(esSpecification(3))
          CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE)
            numberOfComponents = 4
          CASE(EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
            numberOfComponents = 19
          CASE DEFAULT
            localError="Equations set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
              & " is not valid for a Monodomain equation type of monodomain Strang splitting equations set class."
            CALL FlagError(localError,err,error,*999)
          END SELECT          
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSet%INDEPENDENT%independentField,err,error,*999)
          CALL Field_TypeSetAndLock(equationsSet%INDEPENDENT%independentField,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_DependentTypeSetAndLock(equationsSet%INDEPENDENT%independentField,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          CALL Field_DecompositionSetAndLock(equationsSet%INDEPENDENT%independentField,geometricDecomposition, &
                & err,error,*999)
          CALL Field_GeometricFieldSetAndLock(equationsSet%INDEPENDENT%independentField,geometricField,err,error,*999)
          CALL Field_NumberOfVariablesSetAndLock(equationsSet%INDEPENDENT%independentField,1,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsSet%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE,err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, &
            & err,error,*999)         
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE,&
            & numberOfComponents,err,error,*999)
          CALL Field_DOFOrderTypeSet(equationsSet%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE,&
            &  FIELD_CONTIGUOUS_COMPONENT_DOF_ORDER,err,error,*999) ! dofs continuous, so first + (x-1) is x'th component index
          !Default to the geometric interpolation setup
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
          CALL Field_ComponentMeshComponentSet(equationsSet%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE,1, &
            & geometricMeshComponent,err,error,*999)
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            CALL Field_ComponentInterpolationSetAndLock(equationsSet%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE,1, &
              & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            !Default the scaling to the geometric field scaling
            CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
            CALL Field_ScalingTypeSet(equationsSet%INDEPENDENT%independentField,geometricScalingType,err,error,*999)
          END SELECT
        ELSE
          ! user specified field
          CALL FlagError("No user specified field supported!",err,error,*999)
        ENDIF        
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsSet%INDEPENDENT%independentFieldAutoCreated) THEN
          CALL Field_CreateFinish(equationsSet%INDEPENDENT%independentField,err,error,*999)
          SELECT CASE(esSpecification(3))
          CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE)
            !Initialize to y0
            CALL Electrophysiology_BuenoOrovioInitialise(equationsSet%INDEPENDENT%independentField,err,error,*999) 
          CASE(EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
            !Initialize to y0
            CALL Electrophysiology_TenTusscher06Initialise(equationsSet%INDEPENDENT%independentField,err,error,*999)
          CASE DEFAULT
            localError="Equations set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
              & " is not valid for a Monodomain equation type of monodomain Strang splitting equations set class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a Monodomain equation"
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
      !-----------------------------------------------------------------
      ! M a t e r i a l s   f i e l d
      !-----------------------------------------------------------------
      NULLIFY(equationsMaterials)
      CALL EquationsSet_MaterialsGet(equationsSet,equationsMaterials,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Create the auto created materials field
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsMaterials%materialsField,err,error,*999)
          CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
          CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition,err,error,*999)
          CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,geometricField,err,error,*999)
          CALL Field_NumberOfVariablesSetAndLock(equationsMaterials%materialsField,1,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,err,error,*999)
          CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          !Materials field components are
          ! 1. activation  factor (usually 0.0 or 1.0)
          ! 2,3 for fiber/transverse conductivity   . defaults to constant interpolation 
          ! 4,5[,6] : fiber unit vector in dimension
          ! 7: out - activation times
          numberOfMaterialsComponents= 7 !numberOfDimensions + 3
          !Set the number of materials components
          CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
            & numberOfMaterialsComponents,err,error,*999)
          ! 1st = activation = node based, 2 3 diffusion constants
          CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,1, &
            & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
          CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,2, &
            & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
          CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,3, &
            & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
          ! 4 5 (6) fiber unit vector
          DO componentIdx=1,3 !numberOfDimensions
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,componentIdx,geometricMeshComponent, &
              & err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx+3, &
              & geometricMeshComponent,err,error,*999)
            CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx+3, &
              & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
          ENDDO !componentIdx
          CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
            & numberOfMaterialsComponents,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
          !Default the field scaling to that of the geometric field
          CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
          CALL Field_ScalingTypeSet(equationsMaterials%materialsField,geometricScalingType,err,error,*999)
        ELSE
          !User specified field
          CALL FlagError("No user specified field supported!",err,error,*999)
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Finish creating the materials field
          CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
          !Set the default values for the materials field
          !Materials field components are 1 for each dimension i.e., k in div(k.grad(u(x)))
          numberOfMaterialsComponents=numberOfDimensions                             
          !First set the k values to 1.0
          DO componentIdx=1,numberOfDimensions
            CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & componentIdx,1.0_DP,err,error,*999)
          ENDDO !componentIdx
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a Monodomain equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
      !-----------------------------------------------------------------
      ! S o u r c e   f i e l d
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        !Do nothing
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        !Do nothing
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a monodomain equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
      !-----------------------------------------------------------------
      ! E q u a t i o n s    t y p e
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
        NULLIFY(equations)
        CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
        CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
        CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
        SELECT CASE(solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          !Finish the equations creation
          NULLIFY(equations)
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          CALL Equations_CreateFinish(equations,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          !Create the equations mapping.
          NULLIFY(vectorMapping)
          CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
!!! Check this 
          !CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,1,err,error,*999)
          CALL EquationsMappingVector_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
          !CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
          CALL EquationsMappingVector_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
!!!! Check the above two lines
          CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
          CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
          !Create the equations matrices
          NULLIFY(vectorMatrices)
          CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
          CALL Equations_SparsityTypeGet(equations,sparsityType,err,error,*999)
          SELECT CASE(sparsityType)
          CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
            CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,MATRIX_BLOCK_STORAGE_TYPE,err,error,*999)
          CASE(EQUATIONS_MATRICES_SPARSE_MATRICES) 
            CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices, [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE, &
              & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
            CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE, &
              & EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)     
          CASE DEFAULT
            localError="The equations matrices sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
        CASE DEFAULT
          localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))// &
            & " is invalid or not implemented."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a Monodomain equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
        & " is invalid for a Monodomain equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Monodomain_EquationsSetSubtypeSetup")
    RETURN
999 ERRORSEXITS("Monodomain_EquationsSetSubtypeSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Monodomain_EquationsSetSubtypeSetup

  !
  !================================================================================================================================
  !
 
  !>Sets up the Monodomain solution.
  SUBROUTINE Monodomain_ProblemSetup(problem,problemSetup,err,error,*)
    
    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the solutions set to setup a Monodomain equation on.
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Monodomain_ProblemSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(2))
    CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
      CALL Monodomain_ProblemStrangSplittingSetup(problem,problemSetup,err,error,*999)
    CASE DEFAULT
      localError="Problem type "//TRIM(NumberToVString(pSpecification(2),"*",err,error))// &
        & " is not valid for a Monodomain equation Strang splitting problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Monodomain_ProblemSetup")
    RETURN
999 ERRORSEXITS("Monodomain_ProblemSetup",err,error)
    RETURN 1
  END SUBROUTINE Monodomain_ProblemSetup
  
!
  !================================================================================================================================
  !
 
  !>Sets up the Monodomain solution.
  SUBROUTINE Monodomain_ProblemStrangSplittingSetup(problem,problemSetup,err,error,*)
    
    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the solutions set to setup a Monodomain equation on.
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Monodomain_ProblemStrangSplittingSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_MONODOMAIN_BUENOOROVIO_SUBTYPE, &
      & PROBLEM_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
      CALL Monodomain_ProblemSubtypeSetup(problem,problemSetup,err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Monodomain equation type of a Strang splitting problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Monodomain_ProblemStrangSplittingSetup")
    RETURN
999 ERRORS("Monodomain_ProblemStrangSplittingSetup",err,error)
    EXITS("Monodomain_ProblemStrangSplittingSetup")
    RETURN 1
    
  END SUBROUTINE Monodomain_ProblemStrangSplittingSetup
  
  !
  !================================================================================================================================
  !

  !>Sets up the linear Monodomain equations solution.
  SUBROUTINE Monodomain_ProblemSubtypeSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to setup
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop,controlLoopRoot
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Monodomain_ProblemSubtypeSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_MONODOMAIN_BUENOOROVIO_SUBTYPE, &
      & PROBLEM_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    SELECT CASE(problemSetup%setupType)
      !
      ! Initial Setup
      !
    CASE(PROBLEM_SETUP_INITIAL_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Do nothing????
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Do nothing???
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a monodomain equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_CONTROL_TYPE)
      !
      ! Control loop setup
      !
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Set up a simple control loop
        NULLIFY(controlLoop)
        CALL ControlLoop_CreateStart(problem,controlLoop,err,error,*999)
        CALL ControlLoop_TypeSet(controlLoop,CONTROL_TIME_LOOP_TYPE,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Finish the control loops
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        NULLIFY(controlLoop)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        CALL ControlLoop_CreateFinish(controlLoop,err,error,*999)            
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a monodomain equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVERS_TYPE)
      !
      ! Solvers setup
      !
      !Get the control loop
      NULLIFY(controlLoopRoot)
      CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
      NULLIFY(controlLoop)
      CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Start the solvers creation
        NULLIFY(solvers)
        CALL Solvers_CreateStart(controlLoop,solvers,err,error,*999)
        CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
        !Set the solver to be a first order dynamic solver
        NULLIFY(solver)
        CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
        CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
        CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
        !Set solver defaults
        CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
        CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
        CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the solvers
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        !Finish the solvers creation
        CALL Solvers_CreateFinish(solvers,err,error,*999)
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a monodomain equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
      !
      ! Solver equations setup
      !
      NULLIFY(controlLoopRoot)
      CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
      NULLIFY(controlLoop)
      CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
      NULLIFY(solvers)
      CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
      NULLIFY(solver)
      CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Create the solver equations
        NULLIFY(solverEquations)
        CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
        CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
        CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
        CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the solver equations
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        !Finish the solver equations creation
        CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)             
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a monodomain equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
        & " is invalid for a Monodomain equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("Monodomain_ProblemSubtypeSetup")
    RETURN
999 ERRORSEXITS("Monodomain_ProblemSubtypeSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Monodomain_ProblemSubtypeSetup

  !
  !================================================================================================================================
  !

  !>Sets up the output type for a monodomain problem class.
  SUBROUTINE Monodomain_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField,independentField
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("Monodomain_PreSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(2))
    CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(independentField)
      CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)

      CALL Field_ParametersToFieldParametersCopy(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1, &
        & dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,err,error,*999)
      CALL Field_ParametersToFieldParametersCopy(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1, &
        & dependentField,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,1,err,error,*999) ! also to prev.

    CASE DEFAULT
      localError="Problem type "//TRIM(NumberToVString(pSpecification(2),"*",err,error))// &
        & " is not valid for a monodomain problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Monodomain_PreSolve")
    RETURN
999 ERRORSEXITS("Monodomain_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Monodomain_PreSolve

  !
  !================================================================================================================================
  !

  !>Sets up the output type for a monodomain problem class.
  SUBROUTINE Monodomain_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: nodeIdx,numberOfNodes,pSpecification(3)
    REAL(DP) :: currentTime,tmpa,tmpv,timeIncrement
    LOGICAL :: updateDamping,updateRHS,updateStiffness
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatrixType), POINTER :: dampingMatrix,stiffnessMatrix
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,materialsField,independentField,geometricField
    TYPE(FieldVariableType), POINTER :: dependentVariable,independentVariable,materialsVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: solver2
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("Monodomain_PostSolve",err,error,*999)
  
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(2))
    CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
      CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
      NULLIFY(solvers)
      CALL Solver_SolversGet(solver,solvers,err,error,*999)
      NULLIFY(solver2)
      CALL Solvers_SolverGet(solvers,1,solver2,err,error,*999)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver2,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
      NULLIFY(equations)
      CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
      NULLIFY(vectorEquations)
      CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
      NULLIFY(vectorMatrices)
      CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
      NULLIFY(dynamicMatrices)
      CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
      NULLIFY(stiffnessMatrix)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,1,stiffnessMatrix,err,error,*999)
      CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
      NULLIFY(dampingMatrix)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,2,dampingMatrix,err,error,*999)
      CALL EquationsMatrix_UpdateMatrixGet(dampingMatrix,updateStiffness,err,error,*999)
      NULLIFY(rhsVector)
      CALL EquationsMatricesVector_RHSVectorExists(vectorMatrices,rhsVector,err,error,*999)
      updateRHS=.FALSE.
      IF(ASSOCIATED(rhsVector)) CALL EquationsMatricesRHS_UpdateVectorGet(rhsVector,updateRHS,err,error,*999)

      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(dependentVariable)
      CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentVariable,err,error,*999)
      NULLIFY(materialsField)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
      NULLIFY(materialsVariable)
      CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,materialsVariable,err,error,*999)
      NULLIFY(independentField)
      CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
      NULLIFY(independentVariable)
      CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,independentVariable,err,error,*999)
      NULLIFY(domain)
      CALL FieldVariable_ComponentDomainGet(independentVariable,1,domain,err,error,*999)
      NULLIFY(domainTopology)
      CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
      NULLIFY(domainNodes)
      CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
      CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
      
      !Integrate  cell models
 
      CALL Field_ParametersToFieldParametersCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1, &
        & independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,err,error,*999) ! dependent -> independent

      SELECT CASE(pSpecification(3))
      CASE(PROBLEM_MONODOMAIN_BUENOOROVIO_SUBTYPE)
        !From t-dt to t
        CALL Electrophysiology_BuenoOrovioIntegrate(independentField,materialsField,currentTime-timeIncrement,currentTime, &
          & err,error,*999)
      CASE(PROBLEM_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
         !From t-dt to t
        CALL Electrophysiology_TenTusscher06Integrate(independentField,materialsField,currentTime-timeIncrement,currentTime, &
          & err,error,*999)
      CASE DEFAULT
        localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
          & " is not valid for a Monodomain equation type of a Strang splitting problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      
      DO nodeIdx=1,numberOfNodes
        !Default to version 1 of each derivative
        CALL FieldVariable_ParameterSetGetNode(independentVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,1,tmpv,err,error,*999)
        IF(tmpv > 0) THEN
          !Default to version 1 of each derivative
          CALL FieldVariable_ParameterSetGetNode(materialsVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,7,tmpa,err,error,*999) 
          IF(ABS(tmpa)<ZERO_TOLERANCE) THEN
            !Default to version 1 of each derivative
            CALL FieldVariable_ParameterSetUpdateNode(materialsVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,7,currentTime, &
              & err,error,*999)  
          ENDIF
        ENDIF
      ENDDO !nodeIdx

    CASE DEFAULT
      localError="Problem type "//TRIM(NumberToVString(pSpecification(2),"*",err,error))// &
        & " is not valid for a monodomain problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Monodomain_PostSolve")
    RETURN
999 ERRORSEXITS("Monodomain_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Monodomain_PostSolve

  !
  !================================================================================================================================
  !

END MODULE MonodomainEquationsRoutines
