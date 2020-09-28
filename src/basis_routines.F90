!> \file
!> \author Chris Bradley
!> \brief This module contains all basis function routines.
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
!> the terms of any one of the MPL, the GPL or the LGP  L.
!>

!> \defgroup OpenCMISS_Basis OpenCMISS::Iron::Basis
!> This module contains all basis function routines.
MODULE BasisRoutines

  USE BaseRoutines
  USE BasisAccessRoutines
  USE Constants
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE Strings
  USE Types

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  
  !!Module types
  
  !Module variables

  !Interfaces


  !>Evaluates the appropriate partial derivative index for the specificied basis function at a Xi location \see BasisRoutines
  INTERFACE Basis_EvaluateXi
    MODULE PROCEDURE Basis_EvaluateXiDP
  END INTERFACE Basis_EvaluateXi

  !>Interpolates the appropriate partial derivative index of the elements parameters for basis function at a Gauss point \see BasisRoutines
  INTERFACE Basis_InterpolateGauss
    MODULE PROCEDURE Basis_InterpolateGaussDP
  END INTERFACE Basis_InterpolateGauss
  
  !>Interpolates the appropriate partial derivative index of the elements parameters for basis function at Xi location \see BasisRoutines
  INTERFACE Basis_InterpolateXi
    MODULE PROCEDURE Basis_InterpolateXiDP
  END INTERFACE Basis_InterpolateXi

  !>Interpolates the requested partial derivative index(ices) of the element parameters for basis function at a face Gauss point \see BasisRoutines
  INTERFACE Basis_InterpolateLocalFaceGauss
    MODULE PROCEDURE Basis_InterpolateLocalFaceGaussDP
  END INTERFACE Basis_InterpolateLocalFaceGauss
 
  !>Evaluates a quadratic Hermite basis function
  INTERFACE Hermite_QuadraticEvaluate
    MODULE PROCEDURE Hermite_QuadraticEvaluateDP
  END INTERFACE Hermite_QuadraticEvaluate

  !>Evaluates a cubic Hermite basis function
  INTERFACE Hermite_CubicEvaluate
    MODULE PROCEDURE Hermite_CubicEvaluateDP
  END INTERFACE Hermite_CubicEvaluate

  !>Evaluates a linear Lagrange basis function
  INTERFACE Lagrange_LinearEvaluate
    MODULE PROCEDURE Lagrange_LinearEvaluateDP
  END INTERFACE Lagrange_LinearEvaluate

  !>Evaluates a quadratic Lagrange basis function
  INTERFACE Lagrange_QuadraticEvaluate
    MODULE PROCEDURE Lagrange_QuadraticEvaluateDP
  END INTERFACE Lagrange_QuadraticEvaluate

  !>Evaluates a cubic Lagrange basis function
  INTERFACE Lagrange_CubicEvaluate
    MODULE PROCEDURE Lagrange_CubicEvaluateDP
  END INTERFACE Lagrange_CubicEvaluate

  !>Evaluates a linear Simplex basis function
  INTERFACE Simplex_LinearEvaluate
    MODULE PROCEDURE Simplex_LinearEvaluateDP
  END INTERFACE Simplex_LinearEvaluate

  !>Evaluates a quadratic Simplex basis function
  INTERFACE Simplex_QuadraticEvaluate
    MODULE PROCEDURE Simplex_QuadraticEvaluateDP
  END INTERFACE Simplex_QuadraticEvaluate

  !>Evaluates a cubic Simplex basis function
  INTERFACE Simplex_CubicEvaluate
    MODULE PROCEDURE Simplex_CubicEvaluateDP
  END INTERFACE Simplex_CubicEvaluate

  !>Evaluates the Lagrange/Hermite/Fourier tensor product basis function for the given basis
  INTERFACE Basis_LHTPBasisEvaluate
    MODULE PROCEDURE Basis_LHTPBasisEvaluateDP
  END INTERFACE Basis_LHTPBasisEvaluate

  PUBLIC Basis_AreaToXiCoordinates

  PUBLIC Basis_BoundaryXiToXi,Basis_XiToBoundaryXi
  
  PUBLIC Basis_CollapsedXiSet

  PUBLIC Basis_CreateStart,Basis_CreateFinish

  PUBLIC Basis_Destroy
  
  PUBLIC Basis_EvaluateXi

  PUBLIC Basis_GaussPointsCalculate
  
  PUBLIC Basis_InterpolateGauss

  PUBLIC Basis_InterpolateLocalFaceGauss

  PUBLIC Basis_InterpolateXi

  PUBLIC Basis_LocalNodeXiCalculate

  PUBLIC Basis_InterpolationXiSet

  PUBLIC Basis_NumberOfXiSet

  PUBLIC Basis_QuadratureDestroy

  PUBLIC Basis_QuadratureLocalFaceGaussEvaluateSet
  
  PUBLIC Basis_QuadratureNumberOfGaussXiSet

  PUBLIC Basis_QuadratureOrderSet

  PUBLIC Basis_QuadratureTypeSet

  PUBLIC Basis_TypeSet

  PUBLIC Basis_XiToAreaCoordinates

  PUBLIC BasisFunctions_Finalise,BasisFunctions_Initialise

CONTAINS

  !
  !================================================================================================================================
  !

  !>Converts area coordinates to xi coordinates. \see BasisRoutines::Basis_XiToAreaCoordinates
  SUBROUTINE Basis_AreaToXiCoordinates(areaCoordinates,xiCoordinates,err,error,*)
    
    !Argument variables
    REAL(DP), INTENT(IN) :: areaCoordinates(:) !<The area coordinates to convert
    REAL(DP), INTENT(OUT) :: xiCoordinates(:) !<On return, the xi coordinates corresponding to the area coordinates
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Basis_AreaToXiCoordinates",err,error,*999)

    IF(SIZE(areaCoordinates,1)/=(SIZE(xiCoordinates,1)+1)) THEN
      localError="Invalid number of coordinates. The number of area coordinates of "// &
        & TRIM(NumberToVString(SIZE(areaCoordinates,1),"*",err,error))// &
        & " should be equal to the number of xi coordinates of "// &
        & TRIM(NumberToVString(SIZE(xiCoordinates,1),"*",err,error))//" plus one."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    SELECT CASE(SIZE(xiCoordinates,1))
    CASE(1)
      xiCoordinates(1)=1.0_DP-areaCoordinates(1)
    CASE(2)
      xiCoordinates(1)=1.0_DP-areaCoordinates(1)
      xiCoordinates(2)=1.0_DP-areaCoordinates(2)
    CASE(3)
      xiCoordinates(1)=1.0_DP-areaCoordinates(1)
      xiCoordinates(2)=1.0_DP-areaCoordinates(2)
      xiCoordinates(3)=1.0_DP-areaCoordinates(3)
    CASE DEFAULT
      localError="The number of xi coordinates of "//TRIM(NumberToVString(SIZE(xiCoordinates,1),"*",err,error))// &
        & " is invalid. The number must be >= 1 and <= 3."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Basis_AreaToXiCoordinates")
    RETURN
999 ERRORSEXITS("Basis_AreaToXiCoordinates",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_AreaToXiCoordinates

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a new basis \see BasisRoutines::Basis_CreateStart,OpenCMISS::Iron::cmfe_Basis_CreateFinish
  SUBROUTINE Basis_CreateFinish(basis,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: xiIdx,xiCoordIdx,localNodeIdx,localNodeIdx1,localNodeIdx2,localNodeIdx3,localNodeIdx4,elemParamIdx, &
      & localLineIdx,localFaceIdx,columnStart,columnStart2,columnStop,columnStop2
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Basis_CreateFinish",err,error,*999)

    CALL Basis_AssertNotFinished(basis,err,error,*999)
    
    SELECT CASE(basis%TYPE)
    CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
      CALL Basis_LHTPFamilyCreate(basis,err,error,*999)
    CASE(BASIS_SIMPLEX_TYPE)
      CALL Basis_SimplexFamilyCreate(basis,err,error,*999)
    CASE(BASIS_RADIAL_TYPE)
      CALL Basis_RadialFamilyCreate(basis,err,error,*999)
    CASE DEFAULT
      localError="Basis type "//TRIM(NumberToVString(basis%type,"*",err,error))//" is invalid or not implemented"
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    basis%basisFinished=.TRUE.
     
    IF(diagnostics1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Basis user number = ",basis%userNumber,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Basis family number = ",basis%familyNumber,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Basis global number = ",basis%globalNumber,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Basis type = ",basis%type,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Basis degenerate = ",basis%degenerate,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of Xi directions = ",basis%numberOfXi,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of Xi coordinates = ",basis%numberOfXiCoordinates,err,error,*999)
      
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%numberOfXiCoordinates,4,4,basis%interpolationType, &
        & '("  Interpolation type(:)  :",4(X,I2))','(25X,4(X,I2))',err,error,*999)      
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%numberOfXiCoordinates,4,4,basis%interpolationOrder, &
        & '("  Interpolation order(:) :",4(X,I2))','(26X,4(X,I2))',err,error,*999)      
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%numberOfXi,3,3,basis%collapsedXi, &
        & '("  Collapsed xi(:) :",3(X,I2))','(26X,3(X,I2))',err,error,*999)      
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of partial derivatives = ",basis%numberOfPartialDerivatives, &
        & err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of nodes = ",basis%numberOfNodes,err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%numberOfXiCoordinates,4,4,basis%numberOfNodesXiC, &
        & '("  Number of nodes(:)  :",4(X,I2))','(22X,4(X,I2))',err,error,*999)      
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%numberOfNodes,8,8,basis%numberOfDerivatives, &
        & '("  Number of derivatives(:) :",8(X,I2))','(28X,8(X,I2))',err,error,*999)      
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of element parameters = ",basis%numberOfElementParameters, &
        & err,error,*999)
! CPB 23/07/07 Doxygen may or may not like this line for some reason????      
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%numberOfNodes,8,8,basis%nodeAtCollapse, &
        & '("  Node at collapse(:) :",8(X,L1))','(23X,8(X,L1))',err,error,*999)      
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Node position index:",err,error,*999)
      DO xiCoordIdx=1,basis%numberOfXiCoordinates
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Xi coordinate = ",xiCoordIdx,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%numberOfNodes,16,16,basis%nodePositionIndex(:,xiCoordIdx), &
          & '("      Index(:)   :",16(X,I2))','(18X,16(X,I2))',err,error,*999)        
      ENDDO !xiCoordIdx
      
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Inverse node position index:",err,error,*999)
      SELECT CASE(basis%numberOfXiCoordinates)
      CASE(1)
        DO localNodeIdx1=1,basis%numberOfNodesXiC(1)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Xi coordinate = 1, Node position index = ",localNodeIdx1, &
            & err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Index = ",basis%nodePositionIndexInv(localNodeIdx1,1,1,1), &
            & err,error,*999)          
        ENDDO !localNodeIdx1
      CASE(2)
        DO localNodeIdx2=1,basis%numberOfNodesXiC(2)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Xi coordinate = 2, Node position index = ",localNodeIdx2, &
            & err,error,*999)
          DO localNodeIdx1=1,basis%numberOfNodesXiC(1)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Xi coordinate = 1, Node position index = ",localNodeIdx1, &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Index = ",basis%nodePositionIndexInv(localNodeIdx1, &
              & localNodeIdx2,1,1),err,error,*999)
          ENDDO !localNodeIdx1
        ENDDO !localNodeIdx2
      CASE(3)
        DO localNodeIdx3=1,basis%numberOfNodesXiC(3)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Xi coordinate = 3, Node position index = ",localNodeIdx3, &
            & err,error,*999)
          DO localNodeIdx2=1,basis%numberOfNodesXiC(2)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Xi coordinate = 2, Node position index = ",localNodeIdx2, &
              & err,error,*999)
            DO localNodeIdx1=1,basis%numberOfNodesXiC(1)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Xi coordinate = 1, Node position index = ",localNodeIdx1, &
                & err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Index = ",basis%nodePositionIndexInv(localNodeIdx1, &
                & localNodeIdx2,localNodeIdx3,1),err,error,*999)
            ENDDO !localNodeIdx1
          ENDDO !localNodeIdx2
        ENDDO !localNodeIdx3
      CASE(4)
        DO localNodeIdx4=1,basis%numberOfNodesXiC(4)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Xi coordinate = 4, Node position index = ",localNodeIdx4, &
            & err,error,*999)
          DO localNodeIdx3=1,basis%numberOfNodesXiC(3)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Xi coordinate = 3, Node position index = ",localNodeIdx3, &
              & err,error,*999)
            DO localNodeIdx2=1,basis%numberOfNodesXiC(2)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Xi coordinate = 2, Node position index = ",localNodeIdx2, &
                & err,error,*999)
              DO localNodeIdx1=1,basis%numberOfNodesXiC(1) 
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Xi coordinate = 1, Node position index = ",localNodeIdx1, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Index = " &
                  & ,basis%nodePositionIndexInv(localNodeIdx1,localNodeIdx2,localNodeIdx3,localNodeIdx4),err,error,*999)
              ENDDO !localNodeIdx1
            ENDDO !localNodeIdx2
          ENDDO !localNodeIdx3
        ENDDO !localNodeIdx4
      CASE DEFAULT
        CALL FlagError("Invalid number of xi coordinates",err,error,*999)
      END SELECT
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Derivative order index:",err,error,*999)
      DO xiIdx=1,basis%numberOfXi
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Xi = ",xiIdx,err,error,*999)
        DO localNodeIdx=1,basis%numberOfNodes
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Node = ",localNodeIdx,err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%numberOfDerivatives(localNodeIdx),8,8, &
            & basis%derivativeOrderIndex(:,localNodeIdx,xiIdx),'("        Index(:) :",8(X,I2))','(18X,8(X,I2))',err,error,*999)
        ENDDO !localNodeIdx
      ENDDO !xiIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Inverse derivative order index:",err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Element parameter index:",err,error,*999)
      DO localNodeIdx=1,basis%numberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Node = ",localNodeIdx,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%numberOfDerivatives(localNodeIdx),8,8, &
          & basis%elementParameterIndex(:,localNodeIdx),'("      Index(:)   :",8(X,I2))','(18X,8(X,I2))',err,error,*999)
      ENDDO !localNodeIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Inverse element parameter index:",err,error,*999)
      DO elemParamIdx=1,basis%numberOfElementParameters
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Element parameter = ",elemParamIdx,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,2,2, &
          & basis%elementParameterIndexInv(:,elemParamIdx),'("      Index(:)  :",2(X,I2))','(18X,2(X,I2))',err,error,*999)
      ENDDO !elemParamIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Partial derivative index:",err,error,*999)
      DO localNodeIdx=1,basis%numberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Node = ",localNodeIdx,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%numberOfDerivatives(localNodeIdx),8,8, &
          & basis%partialDerivativeIndex(:,localNodeIdx),'("      Index(:)   :",8(X,I2))','(18X,8(X,I2))',err,error,*999)
      ENDDO !localNodeIdx
      IF(basis%numberOfXi==3) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Local faces:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of local faces = ",basis%numberOfLocalFaces,err,error,*999)
        DO localFaceIdx=1,basis%numberOfLocalFaces
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Local face = ",localFaceIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of nodes in local face = ", &
            & basis%numberOfNodesInLocalFace(localFaceIdx),err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%numberOfNodesInLocalFace(localFaceIdx),4,4, &
            & basis%nodeNumbersInLocalFace(:,localFaceIdx),'("      Nodes in local face       :",4(X,I2))','(33X,4(X,I2))', &
            & err,error,*999)
          DO localNodeIdx=1,basis%numberOfNodesInLocalFace(localFaceIdx)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Local face node: ",localNodeIdx,err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%derivativeNumbersInLocalFace(0,localNodeIdx, &
              & localFaceIdx),4,4,basis%derivativeNumbersInLocalFace(1:,localNodeIdx,localFaceIdx), &
              & '("      Derivatives in local face :",4(X,I2))','(33X,4(X,I2))',err,error,*999)
          ENDDO !localNodeIdx
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%numberOfXiCoordinates-1,3,3, &
            & basis%localFaceXiDirections(:,localFaceIdx),'("      Local face xi directions  :",3(X,I2))','(33X,3(X,I2))', &
            & err,error,*999)
          CALL WriteStringFmtValue(DIAGNOSTIC_OUTPUT_TYPE,"      Local face xi normal      : ", &
            & basis%localFaceXiNormal(localFaceIdx),"I2",err,error,*999)
        ENDDO !localFaceIdx
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,2*basis%numberOfXiCoordinates+1,9,9, &
          & basis%xiNormalLocalFace(-basis%numberOfXiCoordinates:basis%numberOfXiCoordinates), &
          & '("    Xi normal local face :",9(X,I2))','(26X,9(X,I2))',err,error,*999)
      ENDIF
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Local lines:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of local lines = ",basis%numberOfLocalLines,err,error,*999)
      DO localLineIdx=1,basis%numberOfLocalLines
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Local line = ",localLineIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of nodes in local line = ", &
          & basis%numberOfNodesInLocalLine(localLineIdx),err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%numberOfNodesInLocalLine(localLineIdx),4,4, &
          & basis%nodeNumbersInLocalLine(:,localLineIdx),'("      Nodes in local line       :",4(X,I2))','(33X,4(X,I2))', &
          & err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%numberOfNodesInLocalLine(localLineIdx),4,4, &
          & basis%derivativeNumbersInLocalLine(:,localLineIdx),'("      Derivatives in local line :",4(X,I2))', &
          & '(33X,4(X,I2))',err,error,*999)
        CALL WriteStringFmtValue(DIAGNOSTIC_OUTPUT_TYPE,"      Local line xi direction   : ", &
          & basis%localLineXiDirection(localLineIdx),"I2",err,error,*999)
        IF(basis%numberOfXi>1) THEN
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%numberOfXi-1,2,2,basis%localLineXiNormals(:,localLineIdx), &
            &  '("      Local line xi normals     :",2(X,I2))','(31X,2(X,I2))',err,error,*999)
        ENDIF
      ENDDO !localLineIdx
      IF(basis%numberOfXi>=2) THEN
        IF(basis%numberOfXi==3) THEN
          columnStart=1
          columnStart2=-basis%numberOfXiCoordinates
          columnStop=2*basis%numberOfXiCoordinates+1
          columnStop2=basis%numberOfXiCoordinates
        ELSE
          columnStart=1
          columnStart2=1
          columnStop=1
          columnStop2=1
        ENDIF
        CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,2*basis%numberOfXiCoordinates+1, &
          & columnStart,1,columnStop,9,9,basis%xiNormalsLocalLine(-basis%numberOfXiCoordinates:basis%numberOfXiCoordinates, &
          & columnStart2:columnStop2),WRITE_STRING_MATRIX_NAME_AND_INDICES, &
          & '("    Xi normal local line','(",I2,",:)',':",9(X,I2))','(31X,9(X,I2))',err,error,*999)
      ENDIF
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of sub-bases = ",basis%numberOfSubBases,err,error,*999)
    ENDIF
    
    EXITS("Basis_CreateFinish")
    RETURN
999 ERRORSEXITS("Basis_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_CreateFinish

  !
  !================================================================================================================================
  !

  !>Starts the creation of a new basis \see BasisRoutines::Basis_CreateFinish,OpenCMISS::Iron::cmfe_Basis_CreateStart
  !>The default values of the BASIS attributes are:
  !>- TYPE: 1 (BASIS_LAGRANGE_HERMITE_TP_TYPE)
  !>- NUMBER_OF_XI: 3
  !>- INTERPOLATION_XI: (1,1,1) (BASIS_LINEAR_LAGRANGE_INTERPOLATIONs)
  !>- INTERPOLATION_TYPE: (1,1,1) (BASIS_LAGRANGE_INTERPOLATIONs)
  !>- INTERPOLATION_ORDER: (1,1,1) (BASIS_LINEAR_INTERPOLATION_ORDERs)
  !>- DEGENERATE: false
  !>- COLLAPSED_XI: (4,4,4) (BASIS_NOT_COLLAPSEDs)
  !>- QUADRATURE: 
  !>  - TYPE: 1 (BASIS_LAGRANGE_HERMITE_TP_TYPE)
  !>  - NUMBER_OF_GAUSS_XI: (2,2,2)
  !>  - GAUSS_ORDER: 0 
  SUBROUTINE Basis_CreateStart(userNumber,basisFunctions,basis,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the basis to start the creation of
    TYPE(BasisFunctionsType), POINTER :: basisFunctions !<The basis functions to create the basis for.
    TYPE(BasisType), POINTER :: basis !<On return, A pointer to the created basis. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: basisIdx,dummyErr
    TYPE(BasisType), POINTER :: newBasis
    TYPE(BasisPtrType), ALLOCATABLE :: newBases(:)
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("Basis_CreateStart",err,error,*999)

    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(basisFunctions)) CALL FlagError("Basis functions is not associated.",err,error,*999)
    
    !See if basis number has already been created
    CALL Basis_UserNumberFind(basisFunctions,userNumber,basis,err,error,*999)
    IF(ASSOCIATED(basis)) THEN
      localError="A basis with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " already exists."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    !Allocate new basis function and add it to the basis functions
    ALLOCATE(newBases(basisFunctions%numberOfBasisFunctions+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new bases.",err,error,*999)
    NULLIFY(newBasis)
    CALL Basis_Initialise(newBasis,err,error,*999)
    newBasis%basisFunctions=>basisFunctions
    DO basisIdx=1,basisFunctions%numberOfBasisFunctions
      newBases(basisIdx)%ptr=>basisFunctions%bases(basisIdx)%ptr
    ENDDO !basisIdx
    basisFunctions%numberOfBasisFunctions=basisFunctions%numberOfBasisFunctions+1
    newBases(basisFunctions%numberOfBasisFunctions)%ptr=>newBasis
    CALL MOVE_ALLOC(newBases,basisFunctions%bases)
    !Set the basis parameters
    newBasis%userNumber=userNumber
    newBasis%familyNumber=0
    newBasis%globalNumber=basisFunctions%numberOfBasisFunctions
    !Set the default basis parameters
    newBasis%type=BASIS_LAGRANGE_HERMITE_TP_TYPE
    newBasis%numberOfXi=3
    ALLOCATE(newBasis%interpolationXi(3),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate basis interpolation xi.",err,error,*999)
    newBasis%interpolationXi(1:3)=[BASIS_LINEAR_LAGRANGE_INTERPOLATION,BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
      & BASIS_LINEAR_LAGRANGE_INTERPOLATION]
    ALLOCATE(newBasis%collapsedXi(3),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate basis collapsed xi.",err,error,*999)
    newBasis%collapsedXi(1:3)=BASIS_NOT_COLLAPSED
    !Initialise the basis quadrature
    NULLIFY(newBasis%QUADRATURE%basis)
    CALL Basis_QuadratureInitialise(newBasis,err,error,*999)        
    basis=>newBasis
    
    EXITS("Basis_CreateStart")
    RETURN
999 IF(ASSOCIATED(newBasis)) CALL BASIS_DESTROY(newBasis,dummyErr,dummyError,*998)
998 IF(ALLOCATED(newBases)) DEALLOCATE(newBases)
    ERRORSEXITS("Basis_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_CreateStart

  !
  !================================================================================================================================
  !

  !>Destroys a basis identified by its basis user number \see BasisRoutines::BASIS_DESTROY_FAMILY,OpenCMISS::Iron::cmfe_Basis_Destroy
  RECURSIVE SUBROUTINE Basis_DestroyNumber(basisFunctions,userNumber,err,error,*)

    !Argument variables
    TYPE(BasisFunctionsType), POINTER :: basisFunctions !<The basis functions for the basis to destroy
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the basis to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Basis_DestroyNumber",err,error,*999)

    CALL Basis_FamilyDestroy(basisFunctions,userNumber,0,err,error,*999)
    
    EXITS("Basis_DestroyNumber")
    RETURN
999 ERRORSEXITS("Basis_DestroyNumber",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_DestroyNumber

  !
  !================================================================================================================================
  !

  !>Destroys a basis. \see BasisRoutines::BASIS_DESTROY_FAMILY,OpenCMISS::Iron::cmfe_Basis_Destroy
  RECURSIVE SUBROUTINE Basis_Destroy(basis,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: userNumber
    TYPE(BasisFunctionsType), POINTER :: basisFunctions
        
    ENTERS("Basis_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)

    NULLIFY(basisFunctions)
    CALL Basis_BasisFunctionsGet(basis,basisFunctions,err,error,*999)
    userNumber=basis%userNumber
    CALL Basis_FamilyDestroy(basisFunctions,userNumber,0,err,error,*999)
    
    EXITS("Basis_Destroy")
    RETURN
999 ERRORSEXITS("Basis_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_Destroy

  !
  !================================================================================================================================
  !

  !>Evaluates the appropriate partial derivative index at position xi for the basis for double precision arguments.
  !>Note for simplex basis functions the xi coordinates should exclude the last area coordinate.
  FUNCTION Basis_EvaluateXiDP(basis,elementParameterIndex,partialDerivativeIndex,xi,err,error)
  
    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: elementParameterIndex !<The element parameter index to evaluate i.e., the local basis index within the element basis.
    INTEGER(INTG), INTENT(IN) :: partialDerivativeIndex !<The partial derivative index to evaluate \see Constants_PartialDerivativeConstants
    REAL(DP), INTENT(IN) :: xi(:) !<The Xi position to evaluate the basis function at
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: Basis_EvaluateXiDP
    !Local Variables
    INTEGER(INTG) :: localNodeIdx,derivativeIdx
    REAL(DP) :: xil(SIZE(xi,1)+1)
    TYPE(VARYING_STRING) :: localError

    ENTERS("Basis_EvaluateXiDP",err,error,*999)
    
    Basis_EvaluateXiDP=0.0_DP
    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)    
    IF(elementParameterIndex<1.OR.elementParameterIndex>basis%numberOfElementParameters) THEN
      localError="The specified element parameter index of "// &
        & TRIM(NumberToVString(elementParameterIndex,"*",err,error))// &
        & " is invalid. The index must be > 0 and <= "// &
        & TRIM(NumberToVString(basis%numberOfElementParameters,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(basis%type)
    CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
      localNodeIdx=basis%elementParameterIndexInv(1,elementParameterIndex)
      derivativeIdx=basis%elementParameterIndexInv(2,elementParameterIndex)
      Basis_EvaluateXiDP=Basis_LHTPBasisEvaluate(basis,localNodeIdx,derivativeIdx,partialDerivativeIndex,xi,err,error)
      IF(err/=0) GOTO 999
    CASE(BASIS_SIMPLEX_TYPE)
      !Create the area coordinates from the xi coordinates
      CALL Basis_XiToAreaCoordinates(xi(1:SIZE(xi,1)),xil(1:SIZE(xi,1)+1),err,error,*999)
      localNodeIdx=basis%elementParameterIndexInv(1,elementParameterIndex)
      Basis_EvaluateXiDP=Basis_SimplexBasisEvaluate(basis,localNodeIdx,partialDerivativeIndex,xil,err,error)
      IF(err/=0) GOTO 999
    CASE(BASIS_RADIAL_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
     CASE DEFAULT
      localError="Basis type "//TRIM(NumberToVString(basis%TYPE,"*",err,error))//" is invalid or not implemented."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Basis_EvaluateXiDP")
    RETURN
999 ERRORSEXITS("Basis_EvaluateXiDP",err,error)
    RETURN
    
  END FUNCTION Basis_EvaluateXiDP
  
  !
  !================================================================================================================================
  !

  !>Converts a xi location on a boundary face or line to a full xi location.
  SUBROUTINE Basis_BoundaryXiToXi(basis,localLineFaceNumber,boundaryXi,fullXi,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to convert the boundary xi for
    INTEGER(INTG), INTENT(IN) :: localLineFaceNumber !<The local line/face number containing the boundary xi
    REAL(DP), INTENT(IN) :: boundaryXi(:) !<The boundary xi location to convert to the full xi location. Note the size of the boundary xi array will determine if we are dealing with a line xi or a face xi.
    REAL(DP), INTENT(OUT) :: fullXi(:) !<On exit, the equivalent full xi location of the boundary xi location
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: normalXi1,normalXi2,numberOfBoundaryXi,numberOfXi
    TYPE(VARYING_STRING) :: localError
        
    ENTERS("Basis_BoundaryXiToXi",err,error,*999)

    CALL Basis_AssertIsFinished(basis,err,error,*999)

    numberOfBoundaryXi=SIZE(boundaryXi,1)    
    numberOfXi=basis%numberOfXi
    IF(numberOfBoundaryXi>numberOfXi) THEN
      localError="The size of the boundary xi array of "//TRIM(NumberToVString(numberOfBoundaryXi,"*",err,error))// &
        & " is invalid. The size must be >= 1 and <= 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(fullXi,1)<numberOfXi) THEN
      localError="The size of the full xi array of "//TRIM(NumberToVString(SIZE(fullXi,1),"*",err,error))// &
        & " must be >= the number of basis xi directions of "//TRIM(NumberToVString(numberOfXi,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    IF(numberOfBoundaryXi==numberOfXi) THEN
      !Basis is of the same dimension as the boundary so just copy over
      fullXi(1:numberOfXi)=boundaryXi(1:numberOfXi)
    ELSE
      SELECT CASE(numberOfBoundaryXi)
      CASE(1)
        !On a line
        IF(localLineFaceNumber<1.OR.localLineFaceNumber>basis%numberOfLocalLines) THEN
          localError="The specified local line/face number of "// &
            & TRIM(NumberToVString(localLineFaceNumber,"*",err,error))// &
            & " is invalid. The local line number must be >=1 and <= "// &
            & TRIM(NumberToVString(basis%numberOfLocalLines,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        SELECT CASE(basis%TYPE)
        CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
          SELECT CASE(numberOfXi)
          CASE(2)
            !2D element
            normalXi1=basis%localLineXiNormals(1,localLineFaceNumber)
            SELECT CASE(normalXi1)
            CASE(-2)
              fullXi(1)=boundaryXi(1)
              fullXi(2)=0.0_DP
            CASE(-1)
              fullXi(1)=0.0_DP
              fullXi(2)=boundaryXi(1)
            CASE(1)
              fullXi(1)=1.0_DP
              fullXi(2)=boundaryXi(1)
            CASE(2)
              fullXi(1)=boundaryXi(1)
              fullXi(2)=1.0_DP
            CASE DEFAULT
              localError="The normal xi direction of "//TRIM(NumberToVString(normalXi1,"*",err,error))// &
                & " is invalid for local line number "//TRIM(NumberToVString(localLineFaceNumber,"*",err,error))// &
                & " on a Lagrange-Hermite basis with two xi directions."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(3)
            !3D element
            normalXi1=basis%localLineXiNormals(1,localLineFaceNumber)
            normalXi2=basis%localLineXiNormals(2,localLineFaceNumber)
            SELECT CASE(normalXi1)
            CASE(-3)
              SELECT CASE(normalXi2)
              CASE(-2)
                fullXi(1)=boundaryXi(1)
                fullXi(2)=0.0_DP
                fullXi(3)=0.0_DP
              CASE(-1)
                fullXi(1)=0.0_DP
                fullXi(2)=boundaryXi(1)
                fullXi(3)=0.0_DP
              CASE(1)
                fullXi(1)=1.0_DP
                fullXi(2)=boundaryXi(1)
                fullXi(3)=0.0_DP
              CASE(2)
                fullXi(1)=boundaryXi(1)
                fullXi(2)=1.0_DP
                fullXi(3)=0.0_DP
              CASE DEFAULT
                localError="The second normal xi direction of "//TRIM(NumberToVString(normalXi2,"*",err,error))// &
                  & " with a first normal xi direction of "//TRIM(NumberToVString(normalXi1,"*",err,error))// &
                  & " is invalid for local line number "//TRIM(NumberToVString(localLineFaceNumber,"*",err,error))// &
                  & " on a Lagrange-Hermite basis with three xi directions."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(-2)
              SELECT CASE(normalXi2)
              CASE(-3)
                fullXi(1)=boundaryXi(1)
                fullXi(2)=0.0_DP
                fullXi(3)=0.0_DP
              CASE(-1)
                fullXi(1)=0.0_DP
                fullXi(2)=0.0_DP
                fullXi(3)=boundaryXi(1)
              CASE(1)
                fullXi(1)=1.0_DP
                fullXi(2)=0.0_DP
                fullXi(3)=boundaryXi(1)
              CASE(3)
                fullXi(1)=boundaryXi(1)
                fullXi(2)=0.0_DP
                fullXi(3)=1.0_DP
              CASE DEFAULT
                localError="The second normal xi direction of "//TRIM(NumberToVString(normalXi2,"*",err,error))// &
                  & " with a first normal xi direction of "//TRIM(NumberToVString(normalXi1,"*",err,error))// &
                  & " is invalid for local line number "//TRIM(NumberToVString(localLineFaceNumber,"*",err,error))// &
                  & " on a Lagrange-Hermite basis with three xi directions."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(-1)
              SELECT CASE(normalXi2)
              CASE(-3)
                fullXi(1)=0.0_DP
                fullXi(2)=boundaryXi(1)
                fullXi(3)=0.0_DP
              CASE(-2)
                fullXi(1)=0.0_DP
                fullXi(2)=0.0_DP
                fullXi(3)=boundaryXi(1)
              CASE(2)
                fullXi(1)=0.0_DP
                fullXi(2)=1.0_DP
                fullXi(3)=boundaryXi(1)
              CASE(3)
                fullXi(1)=0.0_DP
                fullXi(2)=boundaryXi(1)
                fullXi(3)=1.0_DP
              CASE DEFAULT
                localError="The second normal xi direction of "//TRIM(NumberToVString(normalXi2,"*",err,error))// &
                  & " with a first normal xi direction of "//TRIM(NumberToVString(normalXi1,"*",err,error))// &
                  & " is invalid for local line number "//TRIM(NumberToVString(localLineFaceNumber,"*",err,error))// &
                  & " on a Lagrange-Hermite basis with three xi directions."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(1)
              SELECT CASE(normalXi2)
              CASE(-3)
                fullXi(1)=1.0_DP
                fullXi(2)=boundaryXi(1)
                fullXi(3)=0.0_DP
              CASE(-2)
                fullXi(1)=1.0_DP
                fullXi(2)=0.0_DP
                fullXi(3)=boundaryXi(1)
              CASE(2)
                fullXi(1)=1.0_DP
                fullXi(2)=1.0_DP
                fullXi(3)=boundaryXi(1)
              CASE(3)
                fullXi(1)=1.0_DP
                fullXi(2)=boundaryXi(1)
                fullXi(3)=1.0_DP
              CASE DEFAULT
                localError="The second normal xi direction of "//TRIM(NumberToVString(normalXi2,"*",err,error))// &
                  & " with a first normal xi direction of "//TRIM(NumberToVString(normalXi1,"*",err,error))// &
                  & " is invalid for local line number "//TRIM(NumberToVString(localLineFaceNumber,"*",err,error))// &
                  & " on a Lagrange-Hermite basis with three xi directions."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(2)
              SELECT CASE(normalXi2)
              CASE(-3)
                fullXi(1)=boundaryXi(1)
                fullXi(2)=1.0_DP
                fullXi(3)=0.0_DP
              CASE(-1)
                fullXi(1)=0.0_DP
                fullXi(2)=1.0_DP
                fullXi(3)=boundaryXi(1)
              CASE(1)
                fullXi(1)=1.0_DP
                fullXi(2)=1.0_DP
                fullXi(3)=boundaryXi(1)
              CASE(3)
                fullXi(1)=boundaryXi(1)
                fullXi(2)=1.0_DP
                fullXi(3)=1.0_DP
              CASE DEFAULT
                localError="The second normal xi direction of "//TRIM(NumberToVString(normalXi2,"*",err,error))// &
                  & " with a first normal xi direction of "//TRIM(NumberToVString(normalXi1,"*",err,error))// &
                  & " is invalid for local line number "//TRIM(NumberToVString(localLineFaceNumber,"*",err,error))// &
                  & " on a Lagrange-Hermite basis with three xi directions."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(3)
              SELECT CASE(normalXi2)
              CASE(-2)
                fullXi(1)=boundaryXi(1)
                fullXi(2)=0.0_DP
                fullXi(3)=1.0_DP
              CASE(-1)
                fullXi(1)=0.0_DP
                fullXi(2)=boundaryXi(1)
                fullXi(3)=1.0_DP
              CASE(1)
                fullXi(1)=1.0_DP
                fullXi(2)=boundaryXi(1)
                fullXi(3)=1.0_DP
              CASE(2)
                fullXi(1)=boundaryXi(1)
                fullXi(2)=1.0_DP
                fullXi(3)=1.0_DP
              CASE DEFAULT
                localError="The second normal xi direction of "//TRIM(NumberToVString(normalXi2,"*",err,error))// &
                  & " with a first normal xi direction of "//TRIM(NumberToVString(normalXi1,"*",err,error))// &
                  & " is invalid for local line number "//TRIM(NumberToVString(localLineFaceNumber,"*",err,error))// &
                  & " on a Lagrange-Hermite basis with three xi directions."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The normal xi direction of "//TRIM(NumberToVString(normalXi1,"*",err,error))// &
                & " is invalid for a local line "//TRIM(NumberToVString(localLineFaceNumber,"*",err,error))// &
                & " on a basis with two xi directions."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The number of xi directions of "//TRIM(NumberToVString(numberOfXi,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(BASIS_SIMPLEX_TYPE)
          SELECT CASE(numberOfXi)
          CASE(2)
            !2D element
            normalXi1=basis%localLineXiNormals(1,localLineFaceNumber)
            SELECT CASE(normalXi1)
            CASE(1)
              fullXi(1)=0.0_DP
              fullXi(2)=boundaryXi(1)
            CASE(2)
              fullXi(1)=boundaryXi(1)
              fullXi(2)=0.0_DP
            CASE(3)
              fullXi(1)=boundaryXi(1)
              fullXi(2)=1.0_DP-boundaryXi(1)
            CASE DEFAULT
              localError="The normal xi direction of "//TRIM(NumberToVString(normalXi1,"*",err,error))// &
                & " is invalid for a local line "//TRIM(NumberToVString(localLineFaceNumber,"*",err,error))// &
                & " on a simplex basis with two xi directions."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(3)
            !3D element
            normalXi1=basis%localLineXiNormals(1,localLineFaceNumber)
            normalXi2=basis%localLineXiNormals(2,localLineFaceNumber)
            SELECT CASE(normalXi1)
            CASE(1)
              SELECT CASE(normalXi2)
              CASE(2)
                fullXi(1)=0.0_DP
                fullXi(2)=0.0_DP
                fullXi(3)=boundaryXi(1)
              CASE(3)
                fullXi(1)=0.0_DP
                fullXi(2)=boundaryXi(1)
                fullXi(3)=0.0_DP
              CASE(4)
                fullXi(1)=0.0_DP
                fullXi(2)=boundaryXi(1)
                fullXi(3)=1.0_DP-boundaryXi(1)
              CASE DEFAULT
                localError="The second normal xi direction of "//TRIM(NumberToVString(normalXi2,"*",err,error))// &
                  & " with a first normal xi direction of "//TRIM(NumberToVString(normalXi1,"*",err,error))// &
                  & " is invalid for local line number "//TRIM(NumberToVString(localLineFaceNumber,"*",err,error))// &
                  & " on a simplex basis with three xi directions."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(2)
              SELECT CASE(normalXi2)
              CASE(1)
                fullXi(1)=0.0_DP
                fullXi(2)=0.0_DP
                fullXi(3)=boundaryXi(1)
              CASE(3)
                fullXi(1)=boundaryXi(1)
                fullXi(2)=0.0_DP
                fullXi(3)=0.0_DP
              CASE(4)
                fullXi(1)=boundaryXi(1)
                fullXi(2)=0.0_DP
                fullXi(3)=1.0_DP-boundaryXi(1)
              CASE DEFAULT
                localError="The second normal xi direction of "//TRIM(NumberToVString(normalXi2,"*",err,error))// &
                  & " with a first normal xi direction of "//TRIM(NumberToVString(normalXi1,"*",err,error))// &
                  & " is invalid for local line number "//TRIM(NumberToVString(localLineFaceNumber,"*",err,error))// &
                  & " on a simplex basis with three xi directions."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(3)
              SELECT CASE(normalXi2)
              CASE(1)
                fullXi(1)=0.0_DP
                fullXi(2)=boundaryXi(1)
                fullXi(3)=0.0_DP
              CASE(2)
                fullXi(1)=boundaryXi(1)
                fullXi(2)=0.0_DP
                fullXi(3)=0.0_DP
              CASE(4)
                fullXi(1)=boundaryXi(1)
                fullXi(2)=1.0_DP-boundaryXi(1)
                fullXi(3)=0.0_DP
              CASE DEFAULT
                localError="The second normal xi direction of "//TRIM(NumberToVString(normalXi2,"*",err,error))// &
                  & " with a first normal xi direction of "//TRIM(NumberToVString(normalXi1,"*",err,error))// &
                  & " is invalid for local line number "//TRIM(NumberToVString(localLineFaceNumber,"*",err,error))// &
                  & " on a simplex basis with three xi directions."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(4)
              SELECT CASE(normalXi2)
              CASE(1)
                fullXi(1)=0.0_DP
                fullXi(2)=boundaryXi(1)
                fullXi(3)=1.0_DP-boundaryXi(1)
              CASE(2)
                fullXi(1)=boundaryXi(1)
                fullXi(2)=0.0_DP
                fullXi(3)=1.0_DP-boundaryXi(1)
              CASE(3)
                fullXi(1)=boundaryXi(1)
                fullXi(2)=1.0_DP-boundaryXi(1)
                fullXi(3)=0.0_DP
              CASE DEFAULT
                localError="The second normal xi direction of "//TRIM(NumberToVString(normalXi2,"*",err,error))// &
                  & " with a first normal xi direction of "//TRIM(NumberToVString(normalXi1,"*",err,error))// &
                  & " is invalid for local line number "//TRIM(NumberToVString(localLineFaceNumber,"*",err,error))// &
                  & " on a simplex basis with three xi directions."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The normal xi direction of "//TRIM(NumberToVString(normalXi1,"*",err,error))// &
                & " is invalid for a local line "//TRIM(NumberToVString(localLineFaceNumber,"*",err,error))// &
                & " on a simplex basis with three xi directions."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The number of xi directions of "//TRIM(NumberToVString(numberOfXi,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(BASIS_SERENDIPITY_TYPE,BASIS_AUXILLIARY_TYPE,BASIS_B_SPLINE_TP_TYPE, &
          & BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE,BASIS_EXTENDED_LAGRANGE_TP_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The basis type of "//TRIM(NumberToVString(basis%TYPE,"*",err,error))// &
            & " for basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(2)
        !On a face
        IF(localLineFaceNumber<1.OR.localLineFaceNumber>basis%numberOfLocalFaces) THEN
          localError="The specified local line/face number of "// &
            & TRIM(NumberToVString(localLineFaceNumber,"*",err,error))// &
            & " is invalid. The local face number must be >=1 and <= "// &
            & TRIM(NumberToVString(basis%numberOfLocalFaces,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        normalXi1=basis%localFaceXiNormal(localLineFaceNumber)
        SELECT CASE(basis%TYPE)
        CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)                    
          SELECT CASE(normalXi1)
          CASE(-3)
            fullXi(1)=boundaryXi(1)
            fullXi(2)=boundaryXi(2)
            fullXi(3)=0.0_DP
          CASE(-2)
            fullXi(1)=boundaryXi(1)
            fullXi(2)=0.0_DP
            fullXi(3)=boundaryXi(2)
          CASE(-1)
            fullXi(1)=0.0_DP
            fullXi(2)=boundaryXi(1)
            fullXi(3)=boundaryXi(2)
          CASE(1)
            fullXi(1)=1.0_DP
            fullXi(2)=boundaryXi(1)
            fullXi(3)=boundaryXi(2)
          CASE(2)
            fullXi(1)=boundaryXi(1)
            fullXi(2)=1.0_DP
            fullXi(3)=boundaryXi(2)
          CASE(3)
            fullXi(1)=boundaryXi(1)
            fullXi(2)=boundaryXi(2)
            fullXi(3)=1.0_DP
          CASE DEFAULT
            localError="The normal xi direction of "//TRIM(NumberToVString(normalXi1,"*",err,error))// &
              & " is invalid for a local line "//TRIM(NumberToVString(localLineFaceNumber,"*",err,error))// &
              & " on a Lagrange-Hermite basis with three xi directions."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(BASIS_SIMPLEX_TYPE)
          SELECT CASE(normalXi1)
          CASE(1)
            fullXi(1)=0.0_DP
            fullXi(2)=boundaryXi(1)
            fullXi(3)=boundaryXi(2)
          CASE(2)
            fullXi(1)=boundaryXi(1)
            fullXi(2)=0.0_DP
            fullXi(3)=1-boundaryXi(1)-boundaryXi(2)
          CASE(3)
            fullXi(1)=boundaryXi(1)
            fullXi(2)=boundaryXi(2)
            fullXi(3)=0.0_DP
          CASE(4)
            fullXi(1)=boundaryXi(1)
            fullXi(2)=1.0_DP-boundaryXi(1)-boundaryXi(2)
            fullXi(3)=boundaryXi(2)
          CASE DEFAULT
            localError="The normal xi direction of "//TRIM(NumberToVString(normalXi1,"*",err,error))// &
              & " is invalid for a local line "//TRIM(NumberToVString(localLineFaceNumber,"*",err,error))// &
              & " on a simplex basis with three xi directions."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(BASIS_SERENDIPITY_TYPE,BASIS_AUXILLIARY_TYPE,BASIS_B_SPLINE_TP_TYPE, &
          & BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE,BASIS_EXTENDED_LAGRANGE_TP_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The basis type of "//TRIM(NumberToVString(basis%TYPE,"*",err,error))// &
            & " for basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="Invalid number of boundary xi directions. The number of boundary xi directions of "// &
          & TRIM(NumberToVString(numberOfBoundaryXi,"*",err,error))//" must be <= the number of basis xi directions of "// &
          & TRIM(NumberToVString(numberOfXi,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    
    EXITS("Basis_BoundaryXiToXi")
    RETURN
999 ERRORSEXITS("Basis_BoundaryXiToXi",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_BoundaryXiToXi

  !
  !================================================================================================================================
  !

  !>Converts a full xi location on a boundary face or line if the xi location is on a face or line.
  SUBROUTINE Basis_XiToBoundaryXi(basis,fullXi,boundaryXiType,localLineFaceNumber,boundaryXi,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to convert the full xi for
    REAL(DP), INTENT(IN) :: fullXi(:) !<The equivalent full xi location to find the boundary xi location
    INTEGER(INTG), INTENT(OUT) :: boundaryXiType !<On exit, the type of boundary xi location. If the fullXi location is not on the boundary the the xi type will be BASIS_NO_BOUNDARY_XI. If the fullXi location is on a line then the type will be BASIS_LINE_BOUNDARY_XI. If the fullXi location is on a face then the type will be BASIS_FACE_BOUNDARY_XI.
    INTEGER(INTG), INTENT(OUT) :: localLineFaceNumber !<On exit, the local line/face number containing the boundary xi. The the full xi location is not on a face or line the the localLineFaceNumber will be zero.
    REAL(DP), INTENT(OUT) :: boundaryXi(:) !<On exit, the equivalent boundary xi location for the full xi location. Note the size of the boundary xi array will determine if we are dealing with a line xi or a face xi.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: normalDirection1,normalDirection2,numberOfBoundaries,numberOfXi,numberOfXiCoordinates,onBoundary(4),xiIdx
    TYPE(VARYING_STRING) :: localError
        
    ENTERS("Basis_XiToBoundaryXi",err,error,*999)

    boundaryXiType=BASIS_NO_BOUNDARY_XI
    localLineFaceNumber=0
    boundaryXi=0.0_DP
    
    CALL Basis_AssertIsFinished(basis,err,error,*999)
    numberOfXi=basis%numberOfXi
    numberOfXiCoordinates=basis%numberOfXiCoordinates
    IF(SIZE(fullXi,1)<numberOfXi) THEN
      localError="The size of the full xi array of "//TRIM(NumberToVString(SIZE(fullXi,1),"*",err,error))// &
        & " must be >= the number of basis xi directions of "//TRIM(NumberToVString(numberOfXi,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    !Determine if we are near the boundary
    onBoundary=0
    numberOfBoundaries=0
    DO xiIdx=1,numberOfXi
      IF(ABS(fullXi(xiIdx))<ZERO_TOLERANCE) THEN
        onBoundary(xiIdx)=-1
        numberOfBoundaries=numberOfBoundaries+1
      ELSE IF(ABS(fullXi(xiIdx)-1.0_DP)<ZERO_TOLERANCE) THEN
        onBoundary(xiIdx)=+1
        numberOfBoundaries=numberOfBoundaries+1
      ENDIF
    ENDDO !xiIdx
    SELECT CASE(numberOfBoundaries)
    CASE(0) 
      !The xi location is not near the boundary. Return zero
      localLineFaceNumber=0
      boundaryXiType=BASIS_NO_BOUNDARY_XI
      boundaryXi=0.0_DP
    CASE(1)
      !The xi location is on a face in 3D or a line in 2D.
      IF(SIZE(boundaryXi,1)<(numberOfXi-1)) THEN
        localError="The size of the boundary xi array of "//TRIM(NumberToVString(SIZE(boundaryXi,1),"*",err,error))// &
          & " must be >= 2."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Find the normal
      normalDirection1=0
      DO xiIdx=1,numberOfXi
        IF(onBoundary(xiIdx)/=0) THEN
          normalDirection1=onBoundary(xiIdx)
          EXIT
        ENDIF
      ENDDO !xiIdx
      IF(normalDirection1==0) CALL FlagError("Could not find normal direction.",err,error,*999)
      SELECT CASE(basis%type)
      CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
        SELECT CASE(numberOfXi)
        CASE(2)
          localLineFaceNumber=basis%xiNormalsLocalLine(normalDirection1,1)
          boundaryXiType=BASIS_LINE_BOUNDARY_XI
          SELECT CASE(normalDirection1)
          CASE(-2)
            boundaryXi(1)=fullXi(1)
          CASE(-1)
            boundaryXi(1)=fullXi(2)
          CASE(1)
            boundaryXi(1)=fullXi(2)
          CASE(2)
            boundaryXi(1)=fullXi(1)
          CASE DEFAULT
            localError="The normal xi direction of "//TRIM(NumberToVString(normalDirection1,"*",err,error))// &
              & " is invalid on a Lagrange-Hermite basis with two xi directions."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(3)
          localLineFaceNumber=basis%xiNormalLocalFace(normalDirection1)
          boundaryXiType=BASIS_FACE_BOUNDARY_XI
          SELECT CASE(normalDirection1)
          CASE(-3)
            boundaryXi(1)=fullXi(2)
            boundaryXi(2)=fullXi(1)
          CASE(-2)
            boundaryXi(1)=fullXi(1)
            boundaryXi(2)=fullXi(3)
          CASE(-1)
            boundaryXi(1)=fullXi(3)
            boundaryXi(2)=fullXi(2)
          CASE(1)
            boundaryXi(1)=fullXi(2)
            boundaryXi(2)=fullXi(3)
          CASE(2)
            boundaryXi(1)=fullXi(3)
            boundaryXi(2)=fullXi(1)
          CASE(3)
            boundaryXi(1)=fullXi(1)
            boundaryXi(2)=fullXi(2)
          CASE DEFAULT
            localError="The normal xi direction of "//TRIM(NumberToVString(normalDirection1,"*",err,error))// &
              & " is invalid on a Lagrange-Hermite basis with three xi directions."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The number of xi directions of "//TRIM(NumberToVString(numberOfXi,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(BASIS_SIMPLEX_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(BASIS_SERENDIPITY_TYPE,BASIS_AUXILLIARY_TYPE,BASIS_B_SPLINE_TP_TYPE, &
        & BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE,BASIS_EXTENDED_LAGRANGE_TP_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The basis type of "//TRIM(NumberToVString(basis%type,"*",err,error))// &
          & " for basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(2)
      !The xi location is on a line in 3D
      IF(SIZE(boundaryXi,1)<(numberOfXi-1)) THEN
        localError="The size of the boundary xi array of "//TRIM(NumberToVString(SIZE(boundaryXi,1),"*",err,error))// &
          & " must be >= 2."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Find the normal
      normalDirection1=0
      normalDirection2=0
      DO xiIdx=1,numberOfXi
        IF(onBoundary(xiIdx)/=0) THEN
          IF(normalDirection1==0) THEN
            normalDirection1=onBoundary(xiIdx)
          ELSE
            normalDirection2=onBoundary(xiIdx)
            EXIT
          ENDIF
        ENDIF
      ENDDO !xiIdx
      IF(normalDirection1==0) CALL FlagError("Could not find normal direction 1.",err,error,*999)
      IF(normalDirection2==0) CALL FlagError("Could not find normal direction 2.",err,error,*999)
      SELECT CASE(basis%type)
      CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
        SELECT CASE(numberOfXi)
        CASE(3)
          localLineFaceNumber=basis%xiNormalsLocalLine(normalDirection1,normalDirection2)
          boundaryXiType=BASIS_LINE_BOUNDARY_XI
          SELECT CASE(normalDirection1)
          CASE(-3)
            SELECT CASE(normalDirection2)
            CASE(-2)
              boundaryXi(1)=fullXi(1)
            CASE(-1)
              boundaryXi(1)=fullXi(2)
            CASE(1)
              boundaryXi(1)=fullXi(2)
            CASE(2)
              boundaryXi(1)=fullXi(1)
            CASE DEFAULT
              localError="The second normal xi direction of "//TRIM(NumberToVString(normalDirection2,"*",err,error))// &
                & " is invalid on a Lagrange-Hermite basis with a first normal direction of "// &
                & TRIM(NumberToVString(normalDirection1,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            END SELECT            
          CASE(-2)
            SELECT CASE(normalDirection2)
            CASE(-3)
              boundaryXi(1)=fullXi(1)
            CASE(-1)
              boundaryXi(1)=fullXi(3)
            CASE(1)
              boundaryXi(1)=fullXi(3)
            CASE(3)
              boundaryXi(1)=fullXi(1)
            CASE DEFAULT
              localError="The second normal xi direction of "//TRIM(NumberToVString(normalDirection2,"*",err,error))// &
                & " is invalid on a Lagrange-Hermite basis with a first normal direction of "// &
                & TRIM(NumberToVString(normalDirection1,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            END SELECT            
          CASE(-1)
            SELECT CASE(normalDirection2)
            CASE(-3)
              boundaryXi(1)=fullXi(2)
            CASE(-2)
              boundaryXi(1)=fullXi(3)
            CASE(2)
              boundaryXi(1)=fullXi(3)
            CASE(3)
              boundaryXi(1)=fullXi(2)
            CASE DEFAULT
              localError="The second normal xi direction of "//TRIM(NumberToVString(normalDirection2,"*",err,error))// &
                & " is invalid on a Lagrange-Hermite basis with a first normal direction of "// &
                & TRIM(NumberToVString(normalDirection1,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            END SELECT            
          CASE(1)
            SELECT CASE(normalDirection2)
            CASE(-3)
              boundaryXi(1)=fullXi(2)
            CASE(-2)
              boundaryXi(1)=fullXi(3)
            CASE(2)
              boundaryXi(1)=fullXi(3)
            CASE(3)
              boundaryXi(1)=fullXi(2)
            CASE DEFAULT
              localError="The second normal xi direction of "//TRIM(NumberToVString(normalDirection2,"*",err,error))// &
                & " is invalid on a Lagrange-Hermite basis with a first normal direction of "// &
                & TRIM(NumberToVString(normalDirection1,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            END SELECT            
            boundaryXi(1)=fullXi(2)
          CASE(2)
            SELECT CASE(normalDirection2)
            CASE(-3)
              boundaryXi(1)=fullXi(1)
            CASE(-1)
              boundaryXi(1)=fullXi(3)
            CASE(1)
              boundaryXi(1)=fullXi(3)
            CASE(3)
              boundaryXi(1)=fullXi(1)
            CASE DEFAULT
              localError="The second normal xi direction of "//TRIM(NumberToVString(normalDirection2,"*",err,error))// &
                & " is invalid on a Lagrange-Hermite basis with a first normal direction of "// &
                & TRIM(NumberToVString(normalDirection1,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            END SELECT            
            boundaryXi(1)=fullXi(1)
          CASE(3)
            SELECT CASE(normalDirection2)
            CASE(-2)
              boundaryXi(1)=fullXi(1)
            CASE(-1)
              boundaryXi(1)=fullXi(2)
            CASE(1)
              boundaryXi(1)=fullXi(2)
            CASE(2)
              boundaryXi(1)=fullXi(1)
            CASE DEFAULT
              localError="The second normal xi direction of "//TRIM(NumberToVString(normalDirection2,"*",err,error))// &
                & " is invalid on a Lagrange-Hermite basis with three xi directions."
              CALL FlagError(localError,err,error,*999)
            END SELECT            
          CASE DEFAULT
            localError="The first normal xi direction of "//TRIM(NumberToVString(normalDirection1,"*",err,error))// &
              & " is invalid on a Lagrange-Hermite basis with three xi directions."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The number of xi directions of "//TRIM(NumberToVString(numberOfXi,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(BASIS_SIMPLEX_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(BASIS_SERENDIPITY_TYPE,BASIS_AUXILLIARY_TYPE,BASIS_B_SPLINE_TP_TYPE, &
        & BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE,BASIS_EXTENDED_LAGRANGE_TP_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The basis type of "//TRIM(NumberToVString(basis%type,"*",err,error))// &
          & " for basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The number of boundary xi locations of "//TRIM(NumberToVString(numberOfBoundaries,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
          
    EXITS("Basis_XiToBoundaryXi")
    RETURN
999 ERRORSEXITS("Basis_XiToBoundaryXi",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_XiToBoundaryXi

  !
  !================================================================================================================================
  !
  
  !>Calculates the gauss points and weights for a basis function of a particular order for all xi directions. The gauss points will be with respect to xi coordinates. 
  SUBROUTINE Basis_GaussPointsCalculate(basis,order,numberOfXi,numberOfGaussPoints,gaussPoints,gaussWeights,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: order !<The order of interpolation required
    INTEGER(INTG), INTENT(IN) :: numberOfXi !<The number of xi directions of the system in which to calculate gauss points (1D, 2D, 3D)
    INTEGER(INTG), INTENT(OUT) :: numberOfGaussPoints !<On return, the number of gauss points calculated
    REAL(DP), INTENT(OUT) :: gaussPoints(:,:) !<gaussPoints(xiCoordIdx,gaussPointIdx). On return, the calculated gauss point coordinates . Not these will be with respect to xi coordinates (even for simplex elements).
    REAL(DP), INTENT(OUT) :: gaussWeights(:) !<gaussWeights(gaussPointIdx). On return, gauss weight for particular gauss point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: gaussPointIdx,i,j,k,maxNumberOfGauss,numberOfGauss1,numberOfGauss2,numberOfGauss3,xiIdx
    REAL(DP) :: xi(3)
    REAL(DP), ALLOCATABLE :: x(:,:),w(:,:)
    TYPE(VARYING_STRING) :: localError

    ENTERS("Basis_GaussPointsCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(order<1) THEN
      localError="The specified order of "//TRIM(NumberToVString(order,"*",err,error))// &
        & " is invalid. The order must be between >= 1."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(gaussPoints,1)<numberOfXi) THEN
      localError="The size of the first dimension of the specified Gauss points array of "// &
        & TRIM(NumberToVString(SIZE(gaussPoints,1),"*",err,error))//" is too small. The size should be >= "// &
        & TRIM(NumberToVString(numberOfXi,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
        
    !current code assumes same order in each direction
    numberOfGauss1=MAX(CEILING((order+1.0_DP)/2.0_DP),1)
    SELECT CASE(numberOfXi)
    CASE(1)
      numberOfGauss2=1
      numberOfGauss3=1
    CASE(2)
      numberOfGauss2=numberOfGauss1
      numberOfGauss3=1
    CASE(3)
      numberOfGauss2=numberOfGauss1
      numberOfGauss3=numberOfGauss1
    CASE DEFAULT
      localError="The specified number of xi coordinates of " &
        & //TRIM(NumberToVString(numberOfXi,"*",err,error))//" is invalid. The number must be between 1 and 3."
      CALL FlagError(localError,err,error,*999)      
    END SELECT
    maxNumberOfGauss=numberOfGauss1*numberOfGauss2*numberOfGauss3

    SELECT CASE(basis%type)
    CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
      IF(SIZE(gaussPoints,2)<maxNumberOfGauss) THEN
        localError="The size of the second dimension of the specified Gauss points array of "// &
        & TRIM(NumberToVString(SIZE(gaussPoints,2),"*",err,error))//" is too small. The size should be >= "// &
        & TRIM(NumberToVString(maxNumberOfGauss,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(SIZE(gaussWeights,1)<maxNumberOfGauss) THEN
        localError="The size of the second dimension of the specified Gauss weights array of "// &
        & TRIM(NumberToVString(SIZE(gaussPoints,2),"*",err,error))//" is too small. The size should be >= "// &
        & TRIM(NumberToVString(maxNumberOfGauss,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Allocate arrays    
      ALLOCATE(w(numberOfGauss1,3),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate Gauss weights.",err,error,*999)
      ALLOCATE(x(numberOfGauss1,3),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate Gauss point coordinates.",err,error,*999)
      w=1.0_DP
      x=0.0_DP
      !Calculate the one-dimensional Gauss points
      DO xiIdx=1,numberOfXi
        CALL Gauss_Legendre(numberOfGauss1,0.0_DP,1.0_DP,x(:,xiIdx),w(:,xiIdx),err,error,*999)
      ENDDO !xiIdx
      !Form gauss point array for Lagrange-Hermite tensor product.
      numberOfGaussPoints=0
      DO k=1,numberOfGauss3
        DO j=1,numberOfGauss2
          DO i=1,numberOfGauss1
            numberOfGaussPoints=numberOfGaussPoints+1
            xi=[x(i,1),x(j,2),x(k,3)]
            gaussPoints(1:numberOfXi,numberOfGaussPoints)=xi(1:numberOfXi)
            gaussWeights(numberOfGaussPoints)=w(i,1)*w(j,2)*w(k,3)
          ENDDO !i
        ENDDO !j
      ENDDO !k
      DEALLOCATE(x)
      DEALLOCATE(w)
    CASE(BASIS_SIMPLEX_TYPE)      
      !Allocate arrays    
      ALLOCATE(w(maxNumberOfGauss,1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate Gauss weights.",err,error,*999)
      ALLOCATE(x(numberOfXi+1,maxNumberOfGauss),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate Gauss point coordinates.",err,error,*999)
      CALL Gauss_Simplex(order,numberOfXi+1,numberOfGaussPoints,x,w(:,1),err,error,*999)
      DO gaussPointIdx=1,numberOfGaussPoints
        CALL Basis_AreaToXiCoordinates(x(1:numberOfXi+1,gaussPointIdx),xi,err,error,*999)
        gaussPoints(1:numberOfXi,gaussPointIdx)=xi(1:numberOfXi)
        gaussWeights(gaussPointIdx)=w(gaussPointIdx,1)
      ENDDO !gaussPointIdx
      DEALLOCATE(x)
      DEALLOCATE(w)
    CASE DEFAULT
      localError="The specified basis type of "//TRIM(NumberToVString(basis%type,"*",err,error))// &
        & " is invalid or not implemented."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    EXITS("Basis_GaussPointsCalculate")
    RETURN
999 IF(ALLOCATED(x)) DEALLOCATE(x)
    IF(ALLOCATED(w)) DEALLOCATE(w)
    ERRORSEXITS("Basis_GaussPointsCalculate",err,error)
    RETURN 1

  END SUBROUTINE Basis_GaussPointsCalculate

  !
  !================================================================================================================================
  !
  
  !>Destroys a basis identified by its basis user number and family number. Called from the library visible routine Basis_Destroy
  !> \see BasisRoutines::Basis_Destroy
  RECURSIVE SUBROUTINE Basis_FamilyDestroy(basisFunctions,userNumber,familyNumber,err,error,*)

    !Argument variables
    TYPE(BasisFunctionsType), POINTER :: basisFunctions !<The basis functions with the basis to destroy
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the basis to destroy
    INTEGER(INTG), INTENT(IN) :: familyNumber !<The family number of the basis to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: count,basisIdx
    TYPE(BasisType), POINTER :: basis
    TYPE(BasisPtrType), ALLOCATABLE :: newSubBases(:)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Basis_FamilyDestroy",err,error,*999)

    NULLIFY(basis)
    CALL Basis_FamilyNumberFind(basisFunctions,userNumber,familyNumber,basis,err,error,*999)
    IF(.NOT.ASSOCIATED(basis)) THEN
      localError="The basis with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " and a family number of "//TRIM(NumberToVString(familyNumber,"*",err,error))//" does not exist."      
      CALL FlagError(localError,err,error,*999)
    ENDIF

!!NOTE: We have to find a pointer to the basis to destroy within this routine rather than passing in a pointer to a
!!DESTROY_BASIS_PTR type routine because we need to change basis%subBases of the parent basis and this would violate section
!!12.4.1.6 of the Fortran standard if the dummy basis pointer argument was associated with the subBases(x)%ptr actual
!!argument.
      
    IF(basis%numberOfSubBases==0) THEN
      !No more sub-bases so delete this instance
      IF(ASSOCIATED(basis%parentBasis)) THEN
        !Sub-basis function - delete this instance from the parentBasis
        IF(basis%parentBasis%numberOfSubBases>1) THEN
          !If the parent basis has more than one sub basis then remove this instance from its sub-bases list
          ALLOCATE(newSubBases(basis%parentBasis%numberOfSubBases-1),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate new sub-bases.",err,error,*999)
          count=0
          DO basisIdx=1,basis%parentBasis%numberOfSubBases
            IF(basis%parentBasis%subBases(basisIdx)%ptr%userNumber==basis%userNumber.AND. &
              & basis%parentBasis%subBases(basisIdx)%ptr%familyNumber/=basis%familyNumber) THEN
              count=count+1
              newSubBases(count)%ptr=>basis%parentBasis%subBases(basisIdx)%ptr
            ENDIF
          ENDDO !basisIdx
        ENDIF
        basis%parentBasis%numberOfSubBases=basis%parentBasis%numberOfSubBases-1
        CALL MOVE_ALLOC(newSubBases,basis%parentBasis%subBases)
      ELSE
        !Master basis function - delete this instance from basisFunctions
        IF(basisFunctions%numberOfBasisFunctions>1) THEN
          !If there is more than one basis defined then remove this instance from the basis functions
          ALLOCATE(newSubBases(basisFunctions%numberOfBasisFunctions-1),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate new sub-bases.",err,error,*999)
          count=0
          DO basisIdx=1,basisFunctions%numberOfBasisFunctions
            IF(basisFunctions%bases(basisIdx)%ptr%userNumber/=basis%userNumber.AND. &
              & basisFunctions%bases(basisIdx)%ptr%familyNumber==0) THEN
              count=count+1
              newSubBases(count)%ptr=>basisFunctions%bases(basisIdx)%ptr
            ENDIF
          ENDDO !basisIdx
        ENDIF
        CALL MOVE_ALLOC(newSubBases,basisFunctions%bases)
        basisFunctions%numberOfBasisFunctions=basisFunctions%numberOfBasisFunctions-1
      ENDIF

      CALL Basis_Finalise(basis,err,error,*999)
         
    ELSE
      !Recursively delete sub-bases first
      DO WHILE(basis%numberOfSubBases>0)
        CALL Basis_FamilyDestroy(basisFunctions,basis%subBases(1)%ptr%userNumber, &
          & basis%subBases(1)%ptr%familyNumber,err,error,*999)
      ENDDO
      !Now delete this instance
      CALL Basis_FamilyDestroy(basisFunctions,userNumber,familyNumber,err,error,*999)
    ENDIF
    
    EXITS("Basis_FamilyDestroy")
    RETURN
999 ERRORSEXITS("Basis_FamilyDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_FamilyDestroy

  !
  !================================================================================================================================
  !

  !>Finalises a basis and deallocates all memory.
  SUBROUTINE Basis_Finalise(basis,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Basis_Finalise",err,error,*999)

    IF(ASSOCIATED(basis)) THEN
      IF(ALLOCATED(basis%interpolationXi)) DEALLOCATE(basis%interpolationXi)
      IF(ALLOCATED(basis%interpolationType)) DEALLOCATE(basis%interpolationType)
      IF(ALLOCATED(basis%interpolationOrder)) DEALLOCATE(basis%interpolationOrder)
      IF(ALLOCATED(basis%collapsedXi)) DEALLOCATE(basis%collapsedXi)
      IF(ALLOCATED(basis%nodeAtCollapse)) DEALLOCATE(basis%nodeAtCollapse)
      CALL Basis_QuadratureFinalise(basis,err,error,*999)
      IF(ALLOCATED(basis%numberOfNodesXiC)) DEALLOCATE(basis%numberOfNodesXiC)
      IF(ALLOCATED(basis%numberOfDerivatives)) DEALLOCATE(basis%numberOfDerivatives)
      IF(ALLOCATED(basis%nodePositionIndex)) DEALLOCATE(basis%nodePositionIndex)
      IF(ALLOCATED(basis%nodePositionIndexInv)) DEALLOCATE(basis%nodePositionIndexInv)
      IF(ALLOCATED(basis%derivativeOrderIndex)) DEALLOCATE(basis%derivativeOrderIndex)
      IF(ALLOCATED(basis%derivativeOrderIndexInv)) DEALLOCATE(basis%derivativeOrderIndexInv)
      IF(ALLOCATED(basis%partialDerivativeIndex)) DEALLOCATE(basis%partialDerivativeIndex)
      IF(ALLOCATED(basis%elementParameterIndex)) DEALLOCATE(basis%elementParameterIndex)
      IF(ALLOCATED(basis%elementParameterIndexInv)) DEALLOCATE(basis%elementParameterIndexInv)
      IF(ALLOCATED(basis%localLineBasis)) DEALLOCATE(basis%localLineBasis)
      IF(ALLOCATED(basis%localLineXiDirection)) DEALLOCATE(basis%localLineXiDirection)
      IF(ALLOCATED(basis%localLineXiNormals)) DEALLOCATE(basis%localLineXiNormals)
      IF(ALLOCATED(basis%xiNormalsLocalLine)) DEALLOCATE(basis%xiNormalsLocalLine)
      IF(ALLOCATED(basis%numberOfNodesInLocalLine)) DEALLOCATE(basis%numberOfNodesInLocalLine)
      IF(ALLOCATED(basis%nodeNumbersInLocalLine)) DEALLOCATE(basis%nodeNumbersInLocalLine)
      IF(ALLOCATED(basis%derivativeNumbersInLocalLine)) DEALLOCATE(basis%derivativeNumbersInLocalLine)
      IF(ALLOCATED(basis%elementParametersInLocalLine)) DEALLOCATE(basis%elementParametersInLocalLine)
      IF(ALLOCATED(basis%localFaceBasis)) DEALLOCATE(basis%localFaceBasis)
      IF(ALLOCATED(basis%localFaceXiDirections)) DEALLOCATE(basis%localFaceXiDirections)
      IF(ALLOCATED(basis%localFaceXiNormal)) DEALLOCATE(basis%localFaceXiNormal)
      IF(ALLOCATED(basis%xiNormalLocalFace)) DEALLOCATE(basis%xiNormalLocalFace)
      IF(ALLOCATED(basis%numberOfNodesInLocalFace)) DEALLOCATE(basis%numberOfNodesInLocalFace)
      IF(ALLOCATED(basis%nodeNumbersInLocalFace)) DEALLOCATE(basis%nodeNumbersInLocalFace)
      IF(ALLOCATED(basis%derivativeNumbersInLocalFace)) DEALLOCATE(basis%derivativeNumbersInLocalFace)
      IF(ALLOCATED(basis%elementParametersInLocalFace)) DEALLOCATE(basis%elementParametersInLocalFace)
      IF(ALLOCATED(basis%lineBases)) DEALLOCATE(basis%lineBases)
      IF(ALLOCATED(basis%faceBases)) DEALLOCATE(basis%faceBases)
      IF(ALLOCATED(basis%subBases)) DEALLOCATE(basis%subBases)      
      DEALLOCATE(basis)
    ENDIF
   
    EXITS("Basis_Finalise")
    RETURN
999 ERRORSEXITS("Basis_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_Finalise

  !
  !================================================================================================================================
  !

  !>Allocates and initialises a basis.
  SUBROUTINE Basis_Initialise(basis,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to allocate and initialise. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Basis_Initialise",err,error,*998)

    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated.",err,error,*998)

    ALLOCATE(basis,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate basis.",err,error,*999)
    
    basis%userNumber=0
    basis%globalNumber=0
    basis%familyNumber=0
    NULLIFY(basis%basisFunctions)
    basis%basisFinished=.FALSE.
    basis%hermite=.FALSE.
    basis%type=0
    basis%numberOfXi=0
    basis%numberOfXiCoordinates=0
    basis%degenerate=.FALSE.
    basis%numberOfCollapsedXi=0
    basis%numberOfPartialDerivatives=0
    basis%numberOfNodes=0
    basis%numberOfElementParameters=0
    basis%maximumNumberOfDerivatives=0
    basis%numberOfLocalLines=0
    basis%numberOfLocalFaces=0
    basis%numberOfSubBases=0
    NULLIFY(basis%parentBasis)
   
    EXITS("Basis_Initialise")
    RETURN
999 CALL Basis_Finalise(basis,dummyErr,dummyError,*998)
998 ERRORSEXITS("Basis_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_Initialise

  !
  !================================================================================================================================
  ! 

  !>Interpolates the appropriate partial derivative index of the element parameters at a gauss point for the basis
  !>for double precision arguments. Note the interpolated value returned needs to be adjusted for the particular
  !!>coordinate system with CoordinateSystem_InterpolateAdjust. 
  FUNCTION Basis_InterpolateGaussDP(basis,partialDerivativeIndex,quadratureSchemeIdx,gaussPointNumber,elementParameters,err,error)
  
    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: partialDerivativeIndex !<The partial derivative index to interpolate \see Constants_PartialDerivativeConstants
    INTEGER(INTG), INTENT(IN) :: quadratureSchemeIdx !<The quadrature scheme to use \see BasisRoutines_QuadratureSchemes
    INTEGER(INTG), INTENT(IN) :: gaussPointNumber !<The Gauss point number in the scheme to interpolte
    REAL(DP), INTENT(IN) :: elementParameters(:) !<The element parameters to interpolate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: Basis_InterpolateGaussDP
    !Local Variables
    INTEGER(INTG) :: elementParameterIdx
    TYPE(QuadratureSchemeType), POINTER :: quadratureScheme
    TYPE(VARYING_STRING) :: localError

    ENTERS("Basis_InterpolateGaussDP",err,error,*999)
    
    Basis_InterpolateGaussDP=0.0_DP
    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated",err,error,*999)
    NULLIFY(quadratureScheme)
    CALL Basis_QuadratureSchemeGet(basis,quadratureSchemeIdx,quadratureScheme,err,error,*999)
    IF(gaussPointNumber<1.OR.gaussPointNumber>quadratureScheme%numberOfGauss) THEN
      localError="The specified Gauss point number of "//TRIM(NumberToVString(gaussPointNumber,"*",err,error))// &
        & " is invalid. The Gauss point number should be >= 1 and <= "// &
        & TRIM(NumberToVString(quadratureScheme%numberOfGauss,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(partialDerivativeIndex<1.OR.partialDerivativeIndex>basis%numberOfPartialDerivatives) THEN
      localError="The partial derivative index of "//TRIM(NumberToVString(partialDerivativeIndex,"*",err,error))// &
        & " is invalid. It must be between 1 and "// &
        & TRIM(NumberToVString(basis%numberOfPartialDerivatives,"*",err,error))
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO elementParameterIdx=1,basis%numberOfElementParameters
      Basis_InterpolateGaussDP=Basis_InterpolateGaussDP+ &
        & quadratureScheme%gaussBasisFunctions(elementParameterIdx,partialDerivativeIndex,gaussPointNumber)* &
        & elementParameters(elementParameterIdx)
    ENDDO !elementParameterIdx
   
    EXITS("Basis_InterpolateGaussDP")
    RETURN
999 ERRORSEXITS("Basis_InterpolateGaussDP",err,error)
    RETURN
    
  END FUNCTION Basis_InterpolateGaussDP

  !
  !================================================================================================================================
  !

  !>Interpolates the appropriate partial derivative index of the element local face parameters at a face gauss point for the basis
  !>for double precision arguments. Note the interpolated value returned needs to be adjusted for the particular
  !>coordinate system with CoordinateSystem_InterpolateAdjust. 
  FUNCTION Basis_InterpolateLocalFaceGaussDP(basis,partialDerivativeIndex,quadratureSchemeIdx,localFaceNumber,gaussPointNumber, &
    & faceParameters,err,error)
  
    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: partialDerivativeIndex !<The partial derivative index to interpolate \see Constants_PartialDerivativeConstants
    INTEGER(INTG), INTENT(IN) :: quadratureSchemeIdx !<The quadrature scheme to use \see BasisRoutines_QuadratureSchemes
    INTEGER(INTG), INTENT(IN) :: localFaceNumber !<The index number of the face to interpolate on
    INTEGER(INTG), INTENT(IN) :: gaussPointNumber !<The face Gauss point number in the scheme to interpolate
    REAL(DP), INTENT(IN) :: faceParameters(:) !<The face parameters to interpolate (in 3D coordinates)
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: Basis_InterpolateLocalFaceGaussDP
    !Local Variables
    INTEGER(INTG) :: elementParameterIdx
    TYPE(QuadratureSchemeType), POINTER :: quadratureScheme
    TYPE(VARYING_STRING) :: localError

    ENTERS("Basis_InterpolateLocalFaceGaussDP",err,error,*999)
    
    Basis_InterpolateLocalFaceGaussDP=0.0_DP
    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    NULLIFY(quadratureScheme)
    CALL Basis_QuadratureSchemeGet(basis,quadratureSchemeIdx,quadratureScheme,err,error,*999)    
    IF(basis%quadrature%evaluateFaceGauss) &
      & CALL FlagError("The face gauss interpolation scheme has not been created",err,error,*999)
    IF(localFaceNumber<1.OR.localFaceNumber>basis%numberOfLocalFaces) THEN
      localError="The specified local face number of "//TRIM(NumberToVString(localFaceNumber,"*",err,error))// &
        & " is invalid. The Gauss point number should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfLocalFaces,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(gaussPointNumber<1.OR.gaussPointNumber>quadratureScheme%numberOfFaceGauss(localFaceNumber)) THEN
      localError="The specified Gauss point number of "//TRIM(NumberToVString(gaussPointNumber,"*",err,error))// &
        & " is invalid. The Gauss point number should be >= 1 and <= "// &
        & TRIM(NumberToVString(quadratureScheme%numberOfFaceGauss(localFaceNumber),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(partialDerivativeIndex<1.OR.partialDerivativeIndex>basis%numberOfPartialDerivatives) THEN
      localError="The partial derivative index of "//TRIM(NumberToVString(partialDerivativeIndex,"*",err,error))// &
        & " is invalid. It must be between 1 and "// &
        & TRIM(NumberToVString(basis%numberOfPartialDerivatives,"*",err,error))
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO elementParameterIdx=1,basis%numberOfElementParameters
      Basis_InterpolateLocalFaceGaussDP=Basis_InterpolateLocalFaceGaussDP+ &
        & quadratureScheme%faceGaussBasisFunctions(elementParameterIdx,partialDerivativeIndex,gaussPointNumber,localFaceNumber)* &
        & faceParameters(elementParameterIdx)
    ENDDO !elementParameterIdx
    
    EXITS("Basis_InterpolateLocalFaceGaussDP")
    RETURN
999 ERRORSEXITS("Basis_InterpolateLocalFaceGaussDP",err,error)
    RETURN
    
  END FUNCTION Basis_InterpolateLocalFaceGaussDP

  !
  !================================================================================================================================
  !

  !>Interpolates the appropriate partial derivative index of the element parameters at position xi for the basis
  !>for double precision arguments. Note the interpolated value returned needs to be adjusted for the particular
  !>coordinate system with CoordinateSystem_InterpolateAdjust. Note for simplex basis functions the xi coordinates should
  !>exclude the last area coordinate.
  FUNCTION Basis_InterpolateXiDP(basis,partialDerivativeIndex,xi,elementParameters,err,error)
  
    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: partialDerivativeIndex !<The partial derivative index to interpolate \see Constants_PartialDerivativeConstants
    REAL(DP), INTENT(IN) :: xi(:) !<The Xi position to interpolate the basis function at
    REAL(DP), INTENT(IN) :: elementParameters(:) !<The element parameters to interpolate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: Basis_InterpolateXiDP
    !Local Variables
    INTEGER(INTG) :: localNodeIdx,derivativeIdx,elementParameterIdx
    REAL(DP) :: xil(SIZE(xi,1)+1)
    TYPE(VARYING_STRING) :: localError

    ENTERS("Basis_InterpolateXiDP",err,error,*999)
    
    Basis_InterpolateXiDP=0.0_DP
    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated",err,error,*999)
    
    SELECT CASE(basis%type)
    CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
      elementParameterIdx=0
      DO localNodeIdx=1,basis%numberOfNodes
        DO derivativeIdx=1,basis%numberOfDerivatives(localNodeIdx)
          elementParameterIdx=elementParameterIdx+1
          Basis_InterpolateXiDP=Basis_InterpolateXiDP+ &
            & Basis_LHTPBasisEvaluate(basis,localNodeIdx,derivativeIdx,partialDerivativeIndex,xi,err,error)* &
            & elementParameters(elementParameterIdx)
        ENDDO !derivativeIdx
      ENDDO !localNodeIdx
      IF(err/=0) GOTO 999
    CASE(BASIS_SIMPLEX_TYPE)
      !Create the area coordinates from the xi coordinates
      CALL Basis_XiToAreaCoordinates(xi(1:SIZE(xi,1)),xil(1:SIZE(xi,1)+1),err,error,*999)
      elementParameterIdx=0
      DO localNodeIdx=1,basis%numberOfNodes
        elementParameterIdx=elementParameterIdx+1
        Basis_InterpolateXiDP=Basis_InterpolateXiDP+ &
          & Basis_SimplexBasisEvaluate(basis,localNodeIdx,partialDerivativeIndex,xil,err,error)* &
          & elementParameters(elementParameterIdx)
      ENDDO !localNodeIdx
      IF(err/=0) GOTO 999
    CASE(BASIS_RADIAL_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
      IF(err/=0) GOTO 999
    CASE DEFAULT
      localError="Basis type "//TRIM(NumberToVString(basis%TYPE,"*",err,error))//" is invalid or not implemented"
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Basis_InterpolateXiDP")
    RETURN
999 ERRORSEXITS("Basis_InterpolateXiDP",err,error)
    RETURN
    
  END FUNCTION Basis_InterpolateXiDP
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the interpolation type in each xi directions for a basis identified by a pointer. \see OpenCMISS::Iron::cmfe_Basis_InterpolationXiSet
  SUBROUTINE Basis_InterpolationXiSet(basis,interpolationXi,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to set the interpolation xi
    INTEGER(INTG), INTENT(IN) :: interpolationXi(:) !<The interpolation xi parameters for each Xi direction \see BasisRoutines_InterpolationSpecifications
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: xiIdx,lastInterpolation
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Basis_InterpolationXiSet",err,error,*999)

    CALL Basis_AssertNotFinished(basis,err,error,*999)
    IF(SIZE(interpolationXi,1)<basis%numberOfXi) THEN
      localError="The size of the specified interpolation xi array of "// &
        & TRIM(NumberToVString(SIZE(interpolationXi,1),"*",err,error))//" does not match the number of xi directions of "// &
        & TRIM(NumberToVString(basis%numberOfXi,"*",err,error))//" for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    !Check the input values
    SELECT CASE(basis%type)
    CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
      DO xiIdx=1,basis%numberOfXi
        SELECT CASE(interpolationXi(xiIdx))
        CASE(BASIS_LINEAR_LAGRANGE_INTERPOLATION,BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,BASIS_CUBIC_LAGRANGE_INTERPOLATION, &
          & BASIS_CUBIC_HERMITE_INTERPOLATION,BASIS_QUADRATIC1_HERMITE_INTERPOLATION,BASIS_QUADRATIC2_HERMITE_INTERPOLATION)
          !Do nothing
        CASE DEFAULT
          localError="Interpolation xi value "//TRIM(NumberToVString(interpolationXi(xiIdx),"*",err,error))// &
            & " for xi direction "//TRIM(NumberToVString(xiIdx,"*",err,error))//" is invalid for a Lagrange-Hermite TP basis."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !xiIdx
    CASE(BASIS_SIMPLEX_TYPE)
      lastInterpolation=interpolationXi(1)
      DO xiIdx=1,basis%numberOfXi
        SELECT CASE(interpolationXi(xiIdx))
        CASE(BASIS_LINEAR_SIMPLEX_INTERPOLATION,BASIS_QUADRATIC_SIMPLEX_INTERPOLATION,BASIS_CUBIC_SIMPLEX_INTERPOLATION)
          IF(interpolationXi(xiIdx)/=lastInterpolation) THEN
            CALL FlagError("The interpolation xi value must be the same for all xi directions for a simplex basis.", &
              & err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="Interpolation xi value "//TRIM(NumberToVString(interpolationXi(xiIdx),"*",err,error))// &
            & " for xi direction "//TRIM(NumberToVString(xiIdx,"*",err,error))//" is invalid for a simplex basis."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !xiIdx
    CASE DEFAULT
      CALL FlagError("Invalid basis type or not implemented",err,error,*999)
    END SELECT
    
    !Set the interpolation xi
    basis%interpolationXi(1:basis%numberOfXi)=interpolationXi(1:basis%numberOfXi)
    
    EXITS("Basis_InterpolationXiSet")
    RETURN
999 ERRORSEXITS("Basis_InterpolationXiSet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_InterpolationXiSet

  !
  !================================================================================================================================
  !

  !>Creates and initialises a Lagrange-Hermite tensor product basis that has already been allocated Basis_CreateStart
  !> \see BasisRoutines::Basis_CreateStart
  SUBROUTINE Basis_LHTPBasisCreate(basis,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to create
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: maximumNumberOfNodes,numberOfDerivatives,xiIdx,xiIdx1,xiIdx2,xiIdx3,derivativeIdx,localNode,localLineNodeIdx, &
      & elementParameter,oldNumberOfDerivatives,position(4),collapsedPosition(3),maximumNodeExtent(3),collapsedXi(3), &
      & numberOfNodes,numberOfLocalLines,nodeCount,specialNodeCount,nodesInLine(4),numberOfLocalFaces,localFaceIdx, &
      & localNodeIdx,localNodeIdx1,localNodeIdx2,localNodeIdx3,directionIdx,localFaceDerivative,localNodeCount, &
      & localLineParameter,localFaceParameter
    LOGICAL, ALLOCATABLE :: nodeAtCollapse(:)
    LOGICAL :: atCollapse,collapsedFace,firstCollapsedPosition
    
    ENTERS("Basis_LHTPBasisCreate",err,error,*999)

    CALL  Basis_AssertIsLHTPBasis(basis,err,error,*999)
     
    basis%numberOfXiCoordinates=basis%numberOfXi
    basis%numberOfPartialDerivatives=basis%numberOfXiCoordinates**2+2
    ALLOCATE(basis%interpolationType(basis%numberOfXiCoordinates),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate basis interpolation type.",err,error,*999)
    ALLOCATE(basis%interpolationOrder(basis%numberOfXiCoordinates),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate basis interpolation order.",err,error,*999)
    ALLOCATE(basis%numberOfNodesXiC(basis%numberOfXiCoordinates),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate basis number of nodes xic.",err,error,*999)
    numberOfNodes=1
    maximumNumberOfNodes=0
    basis%degenerate=.FALSE.
    basis%numberOfCollapsedXi=0
    DO xiIdx=1,basis%numberOfXi
      !Set up the interpolation types, orders and number of nodes in each xi from the user specified interpolation xi.
      SELECT CASE(basis%interpolationXi(xiIdx))
      CASE(BASIS_LINEAR_LAGRANGE_INTERPOLATION)
        basis%interpolationType(xiIdx)=BASIS_LAGRANGE_INTERPOLATION
        basis%interpolationOrder(xiIdx)=BASIS_LINEAR_INTERPOLATION_ORDER
        basis%numberOfNodesXiC(xiIdx)=2            
      CASE(BASIS_QUADRATIC_LAGRANGE_INTERPOLATION)
        basis%interpolationType(xiIdx)=BASIS_LAGRANGE_INTERPOLATION
        basis%interpolationOrder(xiIdx)=BASIS_QUADRATIC_INTERPOLATION_ORDER
        basis%numberOfNodesXiC(xiIdx)=3
      CASE(BASIS_CUBIC_LAGRANGE_INTERPOLATION)
        basis%interpolationType(xiIdx)=BASIS_LAGRANGE_INTERPOLATION
        basis%interpolationOrder(xiIdx)=BASIS_CUBIC_INTERPOLATION_ORDER
        basis%numberOfNodesXiC(xiIdx)=4
      CASE(BASIS_CUBIC_HERMITE_INTERPOLATION)
        basis%interpolationType(xiIdx)=BASIS_HERMITE_INTERPOLATION
        basis%interpolationOrder(xiIdx)=BASIS_CUBIC_INTERPOLATION_ORDER
        basis%numberOfNodesXiC(xiIdx)=2
      CASE(BASIS_QUADRATIC1_HERMITE_INTERPOLATION)
        basis%interpolationType(xiIdx)=BASIS_HERMITE_INTERPOLATION
        basis%interpolationOrder(xiIdx)=BASIS_QUADRATIC1_INTERPOLATION_ORDER
        basis%numberOfNodesXiC(xiIdx)=2
      CASE(BASIS_QUADRATIC2_HERMITE_INTERPOLATION)
        basis%interpolationType(xiIdx)=BASIS_HERMITE_INTERPOLATION
        basis%interpolationOrder(xiIdx)=BASIS_QUADRATIC2_INTERPOLATION_ORDER
        basis%numberOfNodesXiC(xiIdx)=2
      CASE DEFAULT 
        CALL FlagError("Invalid interpolation type",err,error,*999)
      END SELECT
      IF(basis%collapsedXi(xiIdx)==BASIS_XI_COLLAPSED) THEN
        basis%numberOfCollapsedXi=basis%numberOfCollapsedXi+1
        collapsedXi(basis%numberOfCollapsedXi)=xiIdx
        basis%degenerate=.TRUE.
      ENDIF
      numberOfNodes=numberOfNodes*basis%numberOfNodesXiC(xiIdx)
      IF(basis%numberOfNodesXiC(xiIdx)>maximumNumberOfNodes) maximumNumberOfNodes=basis%numberOfNodesXiC(xiIdx)
    ENDDO !xiIdx
    !If a degenerate (collapsed) basis recalculate the number of nodes from the maximum posible number of nodes
    IF(basis%degenerate) THEN
      !Calculate the nodeAtCollapse array.
      ALLOCATE(nodeAtCollapse(numberOfNodes),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate at collapse",err,error,*999)
      position=1
      basis%numberOfNodes=0
      !Loop over the maximum number of nodes which is currently set for the basis
      DO localNodeIdx=1,numberOfNodes
        atCollapse=.FALSE.
        DO xiIdx=1,basis%numberOfXi
          IF(basis%collapsedXi(xiIdx)==BASIS_COLLAPSED_AT_XI0.AND.position(xiIdx)==1.OR. &
            & basis%collapsedXi(xiIdx)==BASIS_COLLAPSED_AT_XI1.AND.position(xiIdx)==basis%numberOfNodesXiC(xiIdx)) THEN
            atCollapse=.TRUE.
            firstCollapsedPosition=ALL(position(collapsedXi(1:basis%numberOfCollapsedXi))==1)
            EXIT
          ENDIF
        ENDDO !xiIdx
        IF(atCollapse) THEN
          IF(firstCollapsedPosition) THEN
            basis%numberOfNodes=basis%numberOfNodes+1
            nodeAtCollapse(basis%numberOfNodes)=.TRUE.
          ENDIF
        ELSE
          basis%numberOfNodes=basis%numberOfNodes+1
          nodeAtCollapse(basis%numberOfNodes)=.FALSE.
        ENDIF
        position(1)=position(1)+1
        DO xiIdx=1,basis%numberOfXi
          IF(position(xiIdx)>basis%numberOfNodesXiC(xiIdx)) THEN
            position(xiIdx)=1
            position(xiIdx+1)=position(xiIdx+1)+1
          ENDIF
        ENDDO !xiIdx
      ENDDO !localNodeIdx
      CALL MOVE_ALLOC(nodeAtCollapse,basis%nodeAtCollapse)
    ELSE        
      basis%numberOfNodes=numberOfNodes
      ALLOCATE(basis%nodeAtCollapse(basis%numberOfNodes),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate basis node at collapse.",err,error,*999)
      basis%nodeAtCollapse=.FALSE.
      collapsedXi(1)=1
    ENDIF

    ALLOCATE(basis%nodePositionIndex(basis%numberOfNodes,basis%numberOfXiCoordinates),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate basis node position index.",err,error,*999)
    SELECT CASE(basis%numberOfXiCoordinates)
    CASE(1)
      ALLOCATE(basis%nodePositionIndexInv(maximumNumberOfNodes,1,1,1),STAT=err)
    CASE(2)
      ALLOCATE(basis%nodePositionIndexInv(maximumNumberOfNodes,maximumNumberOfNodes,1,1),STAT=err)
    CASE(3)
      ALLOCATE(basis%nodePositionIndexInv(maximumNumberOfNodes,maximumNumberOfNodes,maximumNumberOfNodes,1),STAT=err)
    CASE DEFAULT
      CALL FlagError("Invalid number of xi coordinates.",err,error,*999)
    END SELECT
    IF(err/=0) CALL FlagError("Could not allocate node position index inverse.",err,error,*999)
    basis%nodePositionIndexInv=0
    
    !Determine the node position index and its inverse
    position=1
    collapsedPosition=1
    localNode=0
    firstCollapsedPosition=.TRUE.
    DO localNodeIdx1=1,numberOfNodes
      atCollapse=.FALSE.
      IF(basis%degenerate) THEN
        DO xiIdx=1,basis%numberOfXi
          IF(basis%collapsedXi(xiIdx)==BASIS_COLLAPSED_AT_XI0.AND.position(xiIdx)==1.OR. &
            & basis%collapsedXi(xiIdx)==BASIS_COLLAPSED_AT_XI1.AND.position(xiIdx)==basis%numberOfNodesXiC(xiIdx)) THEN 
            atCollapse=.TRUE.
            firstCollapsedPosition=ALL(position(collapsedXi(1:basis%numberOfCollapsedXi))==1)
            EXIT
          ENDIF
        ENDDO !xiIdx
      ENDIF
      IF(atCollapse) THEN
        IF(firstCollapsedPosition) THEN
          localNode=localNode+1
          basis%nodePositionIndex(localNode,1:basis%numberOfXi)=position(1:basis%numberOfXi)
          basis%nodePositionIndexInv(position(1),position(2),position(3),1)=localNode
        ELSE
          !The second node in the collapsed xi is set to the same node number as the first node in that xi direction.
          collapsedPosition(1:basis%numberOfXi)=position(1:basis%numberOfXi)
          collapsedPosition(collapsedXi(1:basis%numberOfCollapsedXi))=1
          basis%nodePositionIndexInv(position(1),position(2),position(3),1)= &
            & basis%nodePositionIndexInv(collapsedPosition(1),collapsedPosition(2),collapsedPosition(3),1)
        ENDIF
      ELSE
        localNode=localNode+1
        basis%nodePositionIndex(localNode,1:basis%numberOfXi)=position(1:basis%numberOfXi)
        basis%nodePositionIndexInv(position(1),position(2),position(3),1)=localNode
      ENDIF
      position(1)=position(1)+1
      DO xiIdx=1,basis%numberOfXi
        IF(position(xiIdx)>basis%numberOfNodesXiC(xiIdx)) THEN
          position(xiIdx)=1
          position(xiIdx+1)=position(xiIdx+1)+1
        ENDIF
      ENDDO !xiIdx
    ENDDO !localNodeIdx1
    !Calculate the maximum number of derivatives and the number of element parameters
    basis%maximumNumberOfDerivatives=-1
    basis%numberOfElementParameters=0
    DO localNodeIdx=1,basis%numberOfNodes
      numberOfDerivatives=1
      DO xiIdx=1,basis%numberOfXi
        IF((.NOT.basis%nodeAtCollapse(localNodeIdx).OR.basis%collapsedXi(xiIdx)==BASIS_NOT_COLLAPSED).AND. &
          & basis%interpolationType(xiIdx)==BASIS_HERMITE_INTERPOLATION.AND. &
          & (basis%interpolationOrder(xiIdx)==BASIS_CUBIC_INTERPOLATION_ORDER.OR. &
          & (basis%nodePositionIndex(localNodeIdx,xiIdx)==1.AND. &
          & basis%interpolationOrder(xiIdx)==BASIS_QUADRATIC2_INTERPOLATION_ORDER).OR. &
          & (basis%nodePositionIndex(localNodeIdx,xiIdx)==2.AND. &
          & basis%interpolationOrder(xiIdx)==BASIS_QUADRATIC1_INTERPOLATION_ORDER))) THEN
          !Derivative in this direction
          numberOfDerivatives=numberOfDerivatives*2
        ENDIF
      ENDDO !xiIdx
      basis%numberOfElementParameters=basis%numberOfElementParameters+numberOfDerivatives
      IF(numberOfDerivatives>basis%maximumNumberOfDerivatives) basis%maximumNumberOfDerivatives=numberOfDerivatives
    ENDDO !localNodeIdx
    !Now set up the number of derivatives and derivative order index
    ALLOCATE(basis%numberOfDerivatives(basis%numberOfNodes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate number of derivatives.",err,error,*999)
    ALLOCATE(basis%derivativeOrderIndex(basis%maximumNumberOfDerivatives,basis%numberOfNodes,basis%numberOfXi),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate derivative order index.",err,error,*999)
    ALLOCATE(basis%derivativeOrderIndexInv(FIRST_PART_DERIV,FIRST_PART_DERIV,FIRST_PART_DERIV,basis%numberOfNodes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate derivative order index inverse.",err,error,*999)
    ALLOCATE(basis%partialDerivativeIndex(basis%maximumNumberOfDerivatives,basis%numberOfNodes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate partial derivative index.",err,error,*999)
    ALLOCATE(basis%elementParameterIndex(basis%maximumNumberOfDerivatives,basis%numberOfNodes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate element parameter index.",err,error,*999)
    ALLOCATE(basis%elementParameterIndexInv(2,basis%numberOfElementParameters),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate element parameter index inverse.",err,error,*999)
    !Set the derivative order index and its inverse, the element parameter index and the partial derivative index.
    elementParameter=0
    basis%derivativeOrderIndex=0
    basis%derivativeOrderIndexInv=0
    DO localNodeIdx=1,basis%numberOfNodes
      basis%numberOfDerivatives(localNodeIdx)=1
      DO xiIdx1=1,basis%numberOfXi
        IF((.NOT.basis%nodeAtCollapse(localNodeIdx).OR.basis%collapsedXi(xiIdx1)==BASIS_NOT_COLLAPSED).AND. &
          & basis%interpolationType(xiIdx1)==BASIS_HERMITE_INTERPOLATION.AND. &
          & (basis%interpolationOrder(xiIdx1)==BASIS_CUBIC_INTERPOLATION_ORDER.OR. &
          & (basis%nodePositionIndex(localNodeIdx,xiIdx1)==1.AND. &
          & basis%interpolationOrder(xiIdx1)==BASIS_QUADRATIC2_INTERPOLATION_ORDER).OR. &
          & (basis%nodePositionIndex(localNodeIdx,xiIdx1)==2.AND. &
          & basis%interpolationOrder(xiIdx1)==BASIS_QUADRATIC1_INTERPOLATION_ORDER))) THEN
          oldNumberOfDerivatives=basis%numberOfDerivatives(localNodeIdx)
          basis%numberOfDerivatives(localNodeIdx)=basis%numberOfDerivatives(localNodeIdx)*2
          DO derivativeIdx=1,oldNumberOfDerivatives
            basis%derivativeOrderIndex(derivativeIdx,localNodeIdx,xiIdx1)=NO_PART_DERIV
            basis%derivativeOrderIndex(oldNumberOfDerivatives+derivativeIdx,localNodeIdx,xiIdx1)=FIRST_PART_DERIV
            DO xiIdx2=1,xiIdx1-1
              basis%derivativeOrderIndex(oldNumberOfDerivatives+derivativeIdx,localNodeIdx,xiIdx2)= &
                & basis%derivativeOrderIndex(derivativeIdx,localNodeIdx,xiIdx2)
            ENDDO !xiIdx2
          ENDDO !derivativeIdx
        ELSE
          DO derivativeIdx=1,basis%numberOfDerivatives(localNodeIdx)
            basis%derivativeOrderIndex(derivativeIdx,localNodeIdx,xiIdx1)=NO_PART_DERIV
          ENDDO !derivativeIdx
        ENDIF
      ENDDO !xiIdx1
      
      DO derivativeIdx=1,basis%numberOfDerivatives(localNodeIdx)
        elementParameter=elementParameter+1
        basis%elementParameterIndex(derivativeIdx,localNodeIdx)=elementParameter
        basis%elementParameterIndexInv(1,elementParameter)=localNodeIdx
        basis%elementParameterIndexInv(2,elementParameter)=derivativeIdx
        SELECT CASE(basis%numberOfXi)
        CASE(1)
          basis%derivativeOrderIndexInv(basis%derivativeOrderIndex(derivativeIdx,localNodeIdx,1),1,1,localNodeIdx)= &
            & derivativeIdx
          SELECT CASE(basis%derivativeOrderIndex(derivativeIdx,localNodeIdx,1))
          CASE(NO_PART_DERIV)
            basis%partialDerivativeIndex(derivativeIdx,localNodeIdx)=NO_PART_DERIV
          CASE(FIRST_PART_DERIV)
            basis%partialDerivativeIndex(derivativeIdx,localNodeIdx)=PART_DERIV_S1
          CASE DEFAULT
            CALL FlagError("Invalid derivative order index.",err,error,*999)
          END SELECT
        CASE(2)
          basis%derivativeOrderIndexInv(basis%derivativeOrderIndex(derivativeIdx,localNodeIdx,1), &
            & basis%derivativeOrderIndex(derivativeIdx,localNodeIdx,2),1,localNodeIdx)=derivativeIdx
          SELECT CASE(basis%derivativeOrderIndex(derivativeIdx,localNodeIdx,1))
          CASE(NO_PART_DERIV)
            SELECT CASE(basis%derivativeOrderIndex(derivativeIdx,localNodeIdx,2))
            CASE(NO_PART_DERIV)
              basis%partialDerivativeIndex(derivativeIdx,localNodeIdx)=NO_PART_DERIV
            CASE(FIRST_PART_DERIV)
              basis%partialDerivativeIndex(derivativeIdx,localNodeIdx)=PART_DERIV_S2
            CASE DEFAULT
              CALL FlagError("Invalid derivative order index.",err,error,*999)
            END SELECT
          CASE(FIRST_PART_DERIV)
            SELECT CASE(basis%derivativeOrderIndex(derivativeIdx,localNodeIdx,2))
            CASE(NO_PART_DERIV)
              basis%partialDerivativeIndex(derivativeIdx,localNodeIdx)=PART_DERIV_S1
            CASE(FIRST_PART_DERIV)
              basis%partialDerivativeIndex(derivativeIdx,localNodeIdx)=PART_DERIV_S1_S2
            CASE DEFAULT
              CALL FlagError("Invalid derivative order index.",err,error,*999)
            END SELECT
          CASE DEFAULT
            CALL FlagError("Invalid derivative order index.",err,error,*999)
          END SELECT
        CASE(3)
          basis%derivativeOrderIndexInv(basis%derivativeOrderIndex(derivativeIdx,localNodeIdx,1), &
            & basis%derivativeOrderIndex(derivativeIdx,localNodeIdx,2), &
            & basis%derivativeOrderIndex(derivativeIdx,localNodeIdx,3),localNodeIdx)=derivativeIdx
          SELECT CASE(basis%derivativeOrderIndex(derivativeIdx,localNodeIdx,1))
          CASE(NO_PART_DERIV)
            SELECT CASE(basis%derivativeOrderIndex(derivativeIdx,localNodeIdx,2))
            CASE(NO_PART_DERIV)
              SELECT CASE(basis%derivativeOrderIndex(derivativeIdx,localNodeIdx,3))
              CASE(NO_PART_DERIV)                
                basis%partialDerivativeIndex(derivativeIdx,localNodeIdx)=NO_PART_DERIV
              CASE(FIRST_PART_DERIV)
                basis%partialDerivativeIndex(derivativeIdx,localNodeIdx)=PART_DERIV_S3
              CASE DEFAULT
                CALL FlagError("Invalid derivative order index.",err,error,*999)
              END SELECT
            CASE(FIRST_PART_DERIV)
              SELECT CASE(basis%derivativeOrderIndex(derivativeIdx,localNodeIdx,3))
              CASE(NO_PART_DERIV)                
                basis%partialDerivativeIndex(derivativeIdx,localNodeIdx)=PART_DERIV_S2
              CASE(FIRST_PART_DERIV)
                basis%partialDerivativeIndex(derivativeIdx,localNodeIdx)=PART_DERIV_S2_S3
              CASE DEFAULT
                CALL FlagError("Invalid derivative order index.",err,error,*999)
              END SELECT
            CASE DEFAULT
              CALL FlagError("Invalid derivative order index.",err,error,*999)
            END SELECT
          CASE(FIRST_PART_DERIV)
            SELECT CASE(basis%derivativeOrderIndex(derivativeIdx,localNodeIdx,2))
            CASE(NO_PART_DERIV)
              SELECT CASE(basis%derivativeOrderIndex(derivativeIdx,localNodeIdx,3))
              CASE(NO_PART_DERIV)                
                basis%partialDerivativeIndex(derivativeIdx,localNodeIdx)=PART_DERIV_S1
              CASE(FIRST_PART_DERIV)
                basis%partialDerivativeIndex(derivativeIdx,localNodeIdx)=PART_DERIV_S1_S3
              CASE DEFAULT
                CALL FlagError("Invalid derivative order index.",err,error,*999)
              END SELECT
            CASE(FIRST_PART_DERIV)
              SELECT CASE(basis%derivativeOrderIndex(derivativeIdx,localNodeIdx,3))
              CASE(NO_PART_DERIV)                
                basis%partialDerivativeIndex(derivativeIdx,localNodeIdx)=PART_DERIV_S1_S2
              CASE(FIRST_PART_DERIV)
                basis%partialDerivativeIndex(derivativeIdx,localNodeIdx)=PART_DERIV_S1_S2_S3
              CASE DEFAULT
                CALL FlagError("Invalid derivative order index.",err,error,*999)
              END SELECT
            CASE DEFAULT
              CALL FlagError("Invalid derivative order index.",err,error,*999)
            END SELECT
          CASE DEFAULT
            CALL FlagError("Invalid derivative order index.",err,error,*999)
          END SELECT
        CASE DEFAULT
          CALL FlagError("Invalid number of xi direcions.",err,error,*999)
        END SELECT
      ENDDO !derivativeIdx
    ENDDO !localNodeIdx

    !Set up the line information
    SELECT CASE(basis%numberOfXi)
    CASE(1) !1 xi directions
      numberOfLocalLines=1
      basis%numberOfLocalLines=1
      ALLOCATE(basis%numberOfNodesInLocalLine(numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate number of nodes in local line.",err,error,*999)
      basis%numberOfNodesInLocalLine(1)=basis%numberOfNodesXiC(1)
      ALLOCATE(basis%localLineXiDirection(numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate local line xi direction.",err,error,*999)
      basis%localLineXiDirection(1)=1
      ALLOCATE(basis%nodeNumbersInLocalLine(basis%numberOfNodesXiC(1),numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate node numbers in local line.",err,error,*999)
      ALLOCATE(basis%derivativeNumbersInLocalLine(basis%numberOfNodesXiC(1),numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate derivative numbers in local line.",err,error,*999)
      basis%derivativeNumbersInLocalLine=NO_PART_DERIV
      ALLOCATE(basis%elementParametersInLocalLine(basis%numberOfNodesXiC(1)**2,numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate element parameters in local line.",err,error,*999)
      basis%elementParametersInLocalLine=1
      localLineParameter=0
      DO localNodeIdx2=1,basis%numberOfNodesXiC(1)
        DO localNodeIdx1=1,basis%numberOfNodes
          IF(basis%nodePositionIndex(localNodeIdx1,1)==localNodeIdx2) THEN
            basis%nodeNumbersInLocalLine(localNodeIdx2,1)=localNodeIdx1
            DO derivativeIdx=1,basis%numberOfDerivatives(localNodeIdx2)
              localLineParameter=localLineParameter+1
              basis%elementParametersInLocalLine(localLineParameter,1)=basis%elementParameterIndex( &
                & derivativeIdx,localNodeIdx1)
              IF(basis%derivativeOrderIndex(derivativeIdx,localNodeIdx2,1)==FIRST_PART_DERIV) THEN
                basis%derivativeNumbersInLocalLine(localNodeIdx2,1)=derivativeIdx
                EXIT
              ENDIF
            ENDDO !derivativeIdx
            EXIT
          ENDIF
        ENDDO !localNodeIdx2
      ENDDO !localNodeIdx1
    CASE(2) !2 xi directions
      !Determine the maximum node extent of the basis
      maximumNodeExtent(1)=MAXVAL(basis%nodePositionIndex(:,1))
      maximumNodeExtent(2)=MAXVAL(basis%nodePositionIndex(:,2))
      !Allocate and calculate the lines
      numberOfLocalLines=4-basis%numberOfCollapsedXi
      ALLOCATE(basis%numberOfNodesInLocalLine(numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate number of nodes in local line.",err,error,*999)
      basis%numberOfNodesInLocalLine=0
      ALLOCATE(basis%localLineXiDirection(numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate local line xi direction.",err,error,*999)
      ALLOCATE(BASIS%localLineXiNormals(1,numberOfLocalLines),STAT=ERR)
      IF(err/=0) CALL FlagError("Could not allocate local line xi normals.",err,error,*999)
      ALLOCATE(basis%xiNormalsLocalLine(-2:2,1),STAT=ERR)
      basis%xiNormalsLocalLine=0
      IF(err/=0) CALL FlagError("Could not allocate xi normals local line.",err,error,*999)
      ALLOCATE(basis%nodeNumbersInLocalLine(MAXVAL(basis%numberOfNodesXiC),numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate node numbers in local line",err,error,*999)
      basis%nodeNumbersInLocalLine=0
      ALLOCATE(basis%derivativeNumbersInLocalLine(MAXVAL(basis%numberOfNodesXiC),numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate derivative numbers in local line.",err,error,*999)
      basis%derivativeNumbersInLocalLine=NO_PART_DERIV
      ALLOCATE(basis%elementParametersInLocalLine(MAXVAL(basis%numberOfNodesXiC)**2,numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate element parameters in local line.",err,error,*999)
      basis%elementParametersInLocalLine=1
      !Find the lines
      basis%numberOfLocalLines=0
      DO xiIdx1=1,2
        xiIdx2=OTHER_XI_DIRECTIONS2(xiIdx1)
        !We are looking for lines in the xiIdx1 direction from the direction of xiIdx1=0
        !Loop over the element extremes in the xiIdx2 direction
        DO localNodeIdx2=1,maximumNodeExtent(xiIdx2),maximumNodeExtent(xiIdx2)-1
          nodeCount=0
          specialNodeCount=0
          nodesInLine=0
          DO localNodeIdx1=1,basis%numberOfNodes
            IF(basis%collapsedXi(xiIdx1)/=BASIS_NOT_COLLAPSED) THEN
              !The current xi direction, xiIdx1, is in a degenerate plane
              IF(basis%collapsedXi(xiIdx2)==BASIS_XI_COLLAPSED) THEN
                !The other xi direction is collapsed (must be the case)
                IF(basis%collapsedXi(xiIdx1)==BASIS_COLLAPSED_AT_XI0) THEN !Collapsed at the xi=0 end
                  IF(basis%nodePositionIndex(localNodeIdx1,xiIdx2)==localNodeIdx2.OR. &
                    & basis%nodePositionIndex(localNodeIdx1,xiIdx1)==1) THEN
                    nodeCount=nodeCount+1
                    nodesInLine(nodeCount)=localNodeIdx1
                  ENDIF
                ELSE !Collapsed at the xi=1 end
                  IF(basis%nodePositionIndex(localNodeIdx1,xiIdx2)==localNodeIdx2) THEN
                    nodeCount=nodeCount+1
                    nodesInLine(nodeCount)=localNodeIdx1
                  ELSE IF(basis%nodePositionIndex(localNodeIdx1,xiIdx1)==maximumNodeExtent(xiIdx1)) THEN
                    IF(xiIdx1<2) THEN !Special case - put the collapsed node at the end of the line
                      specialNodeCount=specialNodeCount+1
                      nodesInLine(maximumNodeExtent(xiIdx1))=localNodeIdx1
                    ELSE
                      nodeCount=nodeCount+1
                      nodesInLine(nodeCount)=localNodeIdx1
                    ENDIF
                  ENDIF
                ENDIF
              ELSE
                !The current xi direction must be collapsed
                IF(basis%nodePositionIndex(localNodeIdx1,xiIdx2)==localNodeIdx2) THEN
                  nodeCount=nodeCount+1
                  nodesInLine(nodeCount)=localNodeIdx1
                ENDIF
              ENDIF
            ELSE
              !The current xi direction, xiIdx1, is not involved in any collapsed (degenerate) planes
              IF(basis%nodePositionIndex(localNodeIdx1,xiIdx2)==localNodeIdx2) THEN
                nodeCount=nodeCount+1
                nodesInLine(nodeCount)=localNodeIdx1
              ENDIF
            ENDIF
          ENDDO !localNodeIdx1
          IF((nodeCount+specialNodeCount)>1) THEN !More than one node so it is a proper line 
            basis%numberOfLocalLines=basis%numberOfLocalLines+1
            basis%numberOfNodesInLocalLine(basis%numberOfLocalLines)=nodeCount+specialNodeCount
            basis%nodeNumbersInLocalLine(1:basis%numberOfNodesInLocalLine(basis%numberOfLocalLines), &
              & basis%numberOfLocalLines)=nodesInLine(1:basis%numberOfNodesInLocalLine(basis%numberOfLocalLines))
            localLineParameter=0
            DO localLineNodeIdx=1,basis%numberOfNodesInLocalLine(basis%numberOfLocalLines)
              DO derivativeIdx=1,basis%numberOfDerivatives(basis%nodeNumbersInLocalLine(localLineNodeIdx, &
                & basis%numberOfLocalLines))
                IF(basis%derivativeOrderIndex(derivativeIdx,basis%nodeNumbersInLocalLine(localLineNodeIdx, &
                  & basis%numberOfLocalLines),xiIdx2)==NO_PART_DERIV) THEN
                  localLineParameter=localLineParameter+1
                  basis%elementParametersInLocalLine(localLineParameter,basis%numberOfLocalLines)= &
                    & basis%elementParameterIndex(derivativeIdx,basis%nodeNumbersInLocalLine(localLineNodeIdx, &
                    & basis%numberOfLocalLines))
                  IF(basis%derivativeOrderIndex(derivativeIdx,basis%nodeNumbersInLocalLine(localLineNodeIdx, &
                    & basis%numberOfLocalLines),xiIdx1)==FIRST_PART_DERIV) THEN
                    basis%derivativeNumbersInLocalLine(localLineNodeIdx,basis%numberOfLocalLines)=derivativeIdx
                  ENDIF
                ENDIF
              ENDDO !derivativeIdx
            ENDDO !localLineNodeIdx
            basis%localLineXiDirection(basis%numberOfLocalLines)=xiIdx1
            IF(localNodeIdx2==1) THEN
              basis%localLineXiNormals(1,basis%numberOfLocalLines)=-xiIdx2
              basis%xiNormalsLocalLine(-xiIdx2,1)=basis%numberOfLocalLines
            ELSE
              basis%locallineXiNormals(1,basis%numberOfLocalLines)=xiIdx2
              basis%xiNormalsLocalLine(xiIdx2,1)=basis%numberOfLocalLines
            ENDIF
          ENDIF
        ENDDO !localNodeIdx2
      ENDDO !localNodeIdx1
    CASE(3) !3 xi directions
      !Determine the maximum node extent of the basis
      maximumNodeExtent(1)=MAXVAL(basis%nodePositionIndex(:,1))
      maximumNodeExtent(2)=MAXVAL(basis%nodePositionIndex(:,2))
      maximumNodeExtent(3)=MAXVAL(basis%nodePositionIndex(:,3))
      !Allocate and calculate the lines
      IF(basis%numberOfCollapsedXi==1) THEN
        numberOfLocalLines=9
        numberOfLocalFaces=5
      ELSE IF(basis%numberOfCollapsedXi==2) THEN
        numberOfLocalLines=8
        numberOfLocalFaces=5
      ELSE
        numberOfLocalLines=12
        numberOfLocalFaces=6
      ENDIF
      basis%numberOfLocalFaces=numberOfLocalFaces
      
      ALLOCATE(basis%localLineXiDirection(numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate local line xi direction.",err,error,*999)
      ALLOCATE(basis%localLineXiNormals(2,numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate local line xi normals.",err,error,*999)
      ALLOCATE(basis%xiNormalsLocalLine(-3:3,-3:3),STAT=ERR)
      IF(err/=0) CALL FlagError("Could not allocate xi normals local line.",err,error,*999)
      basis%xiNormalsLocalLine=0
      
      ALLOCATE(basis%numberOfNodesInLocalLine(numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate number of nodes in local line.",err,error,*999)
      basis%numberOfNodesInLocalLine=0
      
      ALLOCATE(basis%numberOfNodesInLocalFace(numberOfLocalFaces),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate number of nodes in local face.",err,error,*999)
      basis%numberOfNodesInLocalFace=0
      
      ALLOCATE(basis%nodeNumbersInLocalLine(MAXVAL(basis%numberOfNodesXiC),numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate node numbers in local line.",err,error,*999)
      basis%nodeNumbersInLocalLine=0
      
      ALLOCATE(basis%derivativeNumbersInLocalLine(MAXVAL(basis%numberOfNodesXiC),numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate derivative numbers in local line.",err,error,*999)
      basis%derivativeNumbersInLocalLine=NO_PART_DERIV
      
      ALLOCATE(basis%elementParametersInLocalLine(0:MAXVAL(basis%numberOfNodesXiC)**2,numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate element parameters in local line.",err,error,*999)
      basis%elementParametersInLocalLine=1
      basis%elementParametersInLocalLine(0,:)=1
      
      ALLOCATE(basis%localFaceXiDirections(2,numberOfLocalFaces),STAT=ERR)
      IF(err/=0) CALL FlagError("Could not allocate local face xi directions.",err,error,*999)
      ALLOCATE(basis%localFaceXiNormal(numberOfLocalFaces),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate local face xi direction.",err,error,*999)
      ALLOCATE(basis%xiNormalLocalFace(-3:3),STAT=ERR)
      IF(err/=0) CALL FlagError("Could not allocate xi normal local face.",err,error,*999)
      basis%xiNormalLocalFace=0
      
      ALLOCATE(basis%derivativeNumbersInLocalFace(0:basis%maximumNumberOfDerivatives, &
        & MAXVAL(basis%numberOfNodesXiC)**2,numberOfLocalFaces),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate derivative numbers in local face.",err,error,*999)
      basis%derivativeNumbersInLocalFace=NO_PART_DERIV
      basis%derivativeNumbersInLocalFace(0,:,:)=1
      
      ALLOCATE(basis%elementParametersInLocalFace(0:MAXVAL(basis%numberOfNodesXiC)**2* &
        & basis%maximumNumberOfDerivatives,numberOfLocalFaces),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate element parameters in local face.",err,error,*999)
      basis%elementParametersInLocalFace=1
      basis%elementParametersInLocalFace(0,:)=1
      
      ALLOCATE(basis%nodeNumbersInLocalFace(MAX(maximumNodeExtent(2)*maximumNodeExtent(3), &
        & maximumNodeExtent(3)*maximumNodeExtent(1),maximumNodeExtent(2)*maximumNodeExtent(1)), &
        & numberOfLocalFaces),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate node numbers in local face.",err,error,*999)
      basis%nodeNumbersInLocalFace=0
            
      !Find the lines and faces
      basis%numberOfLocalLines=0
      DO xiIdx1=1,3
        xiIdx2=OTHER_XI_DIRECTIONS3(xiIdx1,2,1)
        xiIdx3=OTHER_XI_DIRECTIONS3(xiIdx1,3,1)
        !We are looking for lines going in the xiIdx1 direction, starting from xiIdx1=0.
        DO localNodeIdx3=1,maximumNodeExtent(xiIdx3),maximumNodeExtent(xiIdx3)-1 
          DO localNodeIdx2=1,maximumNodeExtent(xiIdx2),maximumNodeExtent(xiIdx2)-1
            nodeCount=0
            specialNodeCount=0
            nodesInLine=0
            !Iterate over nodes in the line of interest
            DO localNodeIdx1=1,basis%numberOfNodes
              IF(basis%collapsedXi(xiIdx1)/=BASIS_NOT_COLLAPSED) THEN
                !The current xi direction, xiIdx1, is involved in a collapsed (degenerate) plane 
                IF(basis%collapsedXi(xiIdx2)==BASIS_XI_COLLAPSED.AND.basis%collapsedXi(xiIdx3)==BASIS_XI_COLLAPSED) THEN
                  !Both of the other two xi directions are collapsed
                  IF(basis%collapsedXi(xiIdx1)==BASIS_COLLAPSED_AT_XI0) THEN !Collapsed at the xi=0 end
                    IF((basis%nodePositionIndex(localNodeIdx1,xiIdx2)==localNodeIdx2.OR. &
                      & basis%nodePositionIndex(localNodeIdx1,xiIdx1)==1).AND. &
                      & (basis%nodePositionIndex(localNodeIdx1,xiIdx3)==localNodeIdx3.OR. &
                      & basis%nodePositionIndex(localNodeIdx1,xiIdx1)==1)) THEN
                      nodeCount=nodeCount+1
                      nodesInLine(nodeCount)=localNodeIdx1
                    ENDIF
                  ELSE !Collapsed at the xi=1 end
                    IF(basis%nodePositionIndex(localNodeIdx1,xiIdx2)==localNodeIdx2.AND. &
                      & basis%nodePositionIndex(localNodeIdx1,xiIdx3)==localNodeIdx3) THEN
                      nodeCount=nodeCount+1
                      nodesInLine(nodeCount)=localNodeIdx1
                    ELSE IF(basis%nodePositionIndex(localNodeIdx1,xiIdx1)==maximumNodeExtent(xiIdx1)) THEN
                      IF(xiIdx1<3) THEN !Special case - put the collapsed node at the end of the line
                        specialNodeCount=specialNodeCount+1
                        nodesInLine(maximumNodeExtent(xiIdx1))=localNodeIdx1
                      ELSE
                        nodeCount=nodeCount+1
                        nodesInLine(nodeCount)=localNodeIdx1
                      ENDIF
                    ENDIF
                  ENDIF
                ELSE
                  IF(basis%collapsedXi(xiIdx2)==BASIS_XI_COLLAPSED) THEN
                    !The other xiIdx2 xi direction is collapsed
                    IF(basis%collapsedXi(xiIdx1)==BASIS_COLLAPSED_AT_XI0) THEN !Collapsed at the xi=0 end
                      IF((basis%nodePositionIndex(localNodeIdx1,xiIdx2)==localNodeIdx2.OR. &
                        & basis%nodePositionIndex(localNodeIdx1,xiIdx1)==1).AND. &
                        & basis%nodePositionIndex(localNodeIdx1,xiIdx3)==localNodeIdx3) THEN
                        nodeCount=nodeCount+1
                        nodesInLine(nodeCount)=localNodeIdx1
                      ENDIF
                    ELSE IF(basis%collapsedXi(xiIdx1)==BASIS_COLLAPSED_AT_XI1) THEN !Collapsed at the xi=1 end
                      IF(basis%nodePositionIndex(localNodeIdx1,xiIdx2)==localNodeIdx2.AND. &
                        & basis%nodePositionIndex(localNodeIdx1,xiIdx3)==localNodeIdx3) THEN
                        nodeCount=nodeCount+1
                        nodesInLine(nodeCount)=localNodeIdx1
                      ELSE IF(basis%nodePositionIndex(localNodeIdx1,xiIdx1)==maximumNodeExtent(xiIdx1).AND. &
                        & basis%nodePositionIndex(localNodeIdx1,xiIdx3)==localNodeIdx3) THEN
                        IF(xiIdx1<xiIdx2) THEN !Special case - put the collapsed node at the end of the line
                          specialNodeCount=specialNodeCount+1
                          nodesInLine(maximumNodeExtent(xiIdx1))=localNodeIdx1
                        ELSE
                          nodeCount=nodeCount+1
                          nodesInLine(nodeCount)=localNodeIdx1
                        ENDIF
                      ENDIF
                    ELSE
                      !Not collapsed at a xi end
                      IF(basis%nodePositionIndex(localNodeIdx1,xiIdx2)==localNodeIdx2.AND. &
                        & basis%nodePositionIndex(localNodeIdx1,xiIdx3)==localNodeIdx3) THEN
                        nodeCount=nodeCount+1
                        nodesInLine(nodeCount)=localNodeIdx1
                      ENDIF
                    ENDIF
                  ELSE IF(basis%collapsedXi(xiIdx3)==BASIS_XI_COLLAPSED) THEN
                    !The other xiIdx3 xi direction is collapsed
                    IF(basis%collapsedXi(xiIdx1)==BASIS_COLLAPSED_AT_XI0) THEN !Collapsed at the xi=0 end
                      IF(basis%nodePositionIndex(localNodeIdx1,xiIdx2)==localNodeIdx2.AND. &
                        & (basis%nodePositionIndex(localNodeIdx1,xiIdx3)==localNodeIdx3.OR. &
                        & basis%nodePositionIndex(localNodeIdx1,xiIdx1)==1)) THEN
                        nodeCount=nodeCount+1
                        nodesInLine(nodeCount)=localNodeIdx1
                      ENDIF
                    ELSE IF(basis%collapsedXi(xiIdx1)==BASIS_COLLAPSED_AT_XI1) THEN !Collapsed at the xi=1 end
                      IF(basis%nodePositionIndex(localNodeIdx1,xiIdx2)==localNodeIdx2.AND. &
                        & basis%nodePositionIndex(localNodeIdx1,xiIdx3)==localNodeIdx3) THEN
                        nodeCount=nodeCount+1
                        nodesInLine(nodeCount)=localNodeIdx1
                      ELSE IF(basis%nodePositionIndex(localNodeIdx1,xiIdx1)==maximumNodeExtent(xiIdx1).AND. &
                        & basis%nodePositionIndex(localNodeIdx1,xiIdx2)==localNodeIdx2) THEN
                        IF(xiIdx1<xiIdx3) THEN !Special case - put the collapsed node at the end of the line
                          specialNodeCount=specialNodeCount+1
                          nodesInLine(maximumNodeExtent(xiIdx1))=localNodeIdx1
                        ELSE
                          nodeCount=nodeCount+1
                          nodesInLine(nodeCount)=localNodeIdx1
                        ENDIF
                      ENDIF
                    ELSE
                      !Not collapsed at a xi end
                      IF(basis%nodePositionIndex(localNodeIdx1,xiIdx2)==localNodeIdx2.AND. &
                        & basis%nodePositionIndex(localNodeIdx1,xiIdx3)==localNodeIdx3) THEN
                        nodeCount=nodeCount+1
                        nodesInLine(nodeCount)=localNodeIdx1
                      ENDIF
                    ENDIF
                  ELSE
                    !The current xi must be collapsed
                    IF(basis%nodePositionIndex(localNodeIdx1,xiIdx2)==localNodeIdx2.AND. &
                      & basis%nodePositionIndex(localNodeIdx1,xiIdx3)==localNodeIdx3) THEN
                      nodeCount=nodeCount+1
                      nodesInLine(nodeCount)=localNodeIdx1
                    ENDIF
                  ENDIF
                ENDIF
              ELSE
                !The current xi direction, xiIdx1, is not involved in any collapsed (degenerate) planes
                IF(basis%nodePositionIndex(localNodeIdx1,xiIdx2)==localNodeIdx2.AND. &
                  & basis%nodePositionIndex(localNodeIdx1,xiIdx3)==localNodeIdx3) THEN
                  nodeCount=nodeCount+1
                  nodesInLine(nodeCount)=localNodeIdx1
                ENDIF
              ENDIF
            ENDDO !localNodeIdx1
            IF((nodeCount+specialNodeCount)>1) THEN !More than one node so it is a proper line 
              basis%numberOfLocalLines=basis%numberOfLocalLines+1
              basis%numberOfNodesInLocalLine(basis%numberOfLocalLines)=nodeCount+specialNodeCount
              basis%nodeNumbersInLocalLine(1:basis%numberOfNodesInLocalLine(basis%numberOfLocalLines), &
                & basis%numberOfLocalLines)=nodesInLine(1:basis%numberOfNodesInLocalLine(basis%numberOfLocalLines))
              localLineParameter=0
              DO localLineNodeIdx=1,basis%numberOfNodesInLocalLine(basis%numberOfLocalLines)
                DO derivativeIdx=1,basis%numberOfDerivatives(basis%nodeNumbersInLocalLine(localLineNodeIdx, &
                  & basis%numberOfLocalLines))
                  IF(basis%derivativeOrderIndex(derivativeIdx,basis%nodeNumbersInLocalLine(localLineNodeIdx, &
                    & basis%numberOfLocalLines),xiIdx2)==NO_PART_DERIV.AND. &
                    & basis%derivativeOrderIndex(derivativeIdx,basis%nodeNumbersInLocalLine(localLineNodeIdx, &
                    & basis%numberOfLocalLines),xiIdx3)==NO_PART_DERIV) THEN
                    localLineParameter=localLineParameter+1
                    basis%elementParametersInLocalLine(localLineParameter,basis%numberOfLocalLines)= &
                      & basis%elementParameterIndex(derivativeIdx,basis%nodeNumbersInLocalLine(localLineNodeIdx, &
                      & basis%numberOfLocalLines))
                    IF(basis%derivativeOrderIndex(derivativeIdx,basis%nodeNumbersInLocalLine(localLineNodeIdx, &
                      & basis%numberOfLocalLines),xiIdx1)==FIRST_PART_DERIV) THEN
                      basis%derivativeNumbersInLocalLine(localLineNodeIdx,basis%numberOfLocalLines)=derivativeIdx
                    ENDIF
                  ENDIF
                ENDDO !derivativeIdx
              ENDDO !localLineNodeIdx
              basis%elementParametersInLocalLine(0,basis%numberOfLocalLines)=localLineParameter
              basis%localLineXiDirection(basis%numberOfLocalLines)=xiIdx1
              IF(localNodeIdx2==1) THEN
                basis%localLineXiNormals(1,basis%numberOfLocalLines)=-xiIdx2
                IF(localNodeIdx3==1) THEN
                  basis%localLineXiNormals(2,basis%numberOfLocalLines)=-xiIdx3
                  basis%xiNormalsLocalLine(-xiIdx2,-xiIdx3)=basis%numberOfLocalLines
                ELSE
                  basis%localLineXiNormals(2,basis%numberOfLocalLines)=xiIdx3
                  basis%xiNormalsLocalLine(-xiIdx2,xiIdx3)=basis%numberOfLocalLines
                ENDIF
              ELSE
                basis%localLineXiNormals(1,basis%numberOfLocalLines)=xiIdx2
                IF(localNodeIdx3==1) THEN
                  basis%localLineXiNormals(2,basis%numberOfLocalLines)=-xiIdx3
                  basis%xiNormalsLocalLine(xiIdx2,-xiIdx3)=basis%numberOfLocalLines
                ELSE
                  basis%localLineXiNormals(2,basis%numberOfLocalLines)=xiIdx3
                  basis%xiNormalsLocalLine(xiIdx2,xiIdx3)=basis%numberOfLocalLines
                ENDIF
              ENDIF
            ENDIF
          ENDDO !localNodeIdx2
        ENDDO !localNodeIdx3
      ENDDO !xiIdx1
      
      !Find the local nodes and derivatives in each face and the local face xi direction
      localFaceIdx=0
      !Loop over the -'ve and +'ve xi direction
      DO directionIdx=-1,1,2
        !Loop over the three xi directions
        DO xiIdx1=1,3
          !xiIdx1 is the +/- face normal direction. xiIdx2 and xiIdx3 are the xi directions in the face.
          xiIdx2=OTHER_XI_DIRECTIONS3(xiIdx1,2,1)
          xiIdx3=OTHER_XI_DIRECTIONS3(xiIdx1,3,1)

          IF(directionIdx==1) THEN
            !The +'ve xi direction
            localNodeIdx1=maximumNodeExtent(xiIdx1)
            !Compute if the face in the +xiIdx1 direction is collapsed.
            collapsedFace=basis%collapsedXi(xiIdx1)==BASIS_COLLAPSED_AT_XI1
          ELSE
            !The -'ve xi direction
            localNodeIdx1=1
            !Compute if the face in the -xiIdx1 direction is collapsed.
            collapsedFace=basis%collapsedXi(xiIdx1)==BASIS_COLLAPSED_AT_XI0
          ENDIF
          localNodeCount=0
          IF(.NOT.collapsedFace) THEN
            !If the face has not been collapsed
            localFaceIdx=localFaceIdx+1
            !Loop over the local nodes in the face
            DO localNodeIdx3=1,maximumNodeExtent(xiIdx2)
              DO localNodeIdx2=1,maximumNodeExtent(xiIdx3)
                IF(xiIdx1==1) THEN
                  localNodeIdx=basis%nodePositionIndexInv(localNodeIdx1,localNodeIdx2,localNodeIdx3,1)
                ELSE IF(xiIdx1==2) THEN
                  localNodeIdx=basis%nodePositionIndexInv(localNodeIdx2,localNodeIdx1,localNodeIdx3,1)
                ELSE
                  localNodeIdx=basis%nodePositionIndexInv(localNodeIdx2,localNodeIdx3,localNodeIdx1,1)
                ENDIF
                IF(ALL(basis%nodeNumbersInLocalFace(1:localNodeCount,localFaceIdx)/=localNodeIdx)) THEN
                  !The node hasn't been collapsed
                  localNodeCount=localNodeCount+1
                  basis%nodeNumbersInLocalFace(localNodeCount,localFaceIdx)=localNodeIdx
                ENDIF
              ENDDO !localNodeIdx3
            ENDDO !localNodexIdx2
            basis%numberOfNodesInLocalFace(localFaceIdx)=localNodeCount
            basis%localFaceXiDirections(1,localFaceIdx)=xiIdx2
            basis%localFaceXiDirections(2,localFaceIdx)=xiIdx3
            basis%localFaceXiNormal(localFaceIdx)=directionIdx*xiIdx1
            basis%xiNormalLocalFace(directionIdx*xiIdx1)=localFaceIdx
            !Compute derivatives and element parameters in the face
            localFaceParameter=0
            DO localNodeIdx=1,basis%numberOfNodesInLocalFace(localFaceIdx)
              localNode=basis%nodeNumbersInLocalFace(localNodeIdx,localFaceIdx)
              localFaceDerivative=0
              DO derivativeIdx=1,basis%numberOfDerivatives(localNode)
                IF(basis%derivativeOrderIndex(derivativeIdx,localNode,xiIdx1)==NO_PART_DERIV) THEN
                  localFaceParameter=localFaceParameter+1
                  localFaceDerivative=localFaceDerivative+1
                  basis%derivativeNumbersInLocalFace(localFaceDerivative,localNodeIdx,localFaceIdx)=derivativeIdx
                  basis%elementParametersInLocalFace(localFaceParameter,localFaceIdx)= &
                    & basis%elementParameterIndex(derivativeIdx,localNode)
                ENDIF
              ENDDO !derivativeIdx
              basis%derivativeNumbersInLocalFace(0,localNodeIdx,localFaceIdx)=localFaceDerivative
            ENDDO !localNodeIdx
            basis%elementParametersInLocalFace(0,localFaceIdx)=localFaceParameter
          ENDIF
        ENDDO !xiIdx1
      ENDDO !directionIdx
    CASE DEFAULT
      CALL FlagError("Invalid number of xi directions.",err,error,*999)
    END SELECT

    CALL Basis_QuadratureCreate(basis,err,error,*999)

    EXITS("Basis_LHTPBasisCreate")
    RETURN
999 IF(ALLOCATED(nodeAtCollapse)) DEALLOCATE(nodeAtCollapse)
    ERRORSEXITS("Basis_LHTPBasisCreate",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_LHTPBasisCreate

  !
  !
  !================================================================================================================================
  !
  
  !>Evaluates the double precision Lagrange/Hermite/Fourier tensor product basis function for the given basis.
  FUNCTION Basis_LHTPBasisEvaluateDP(basis,nodeNumber,derivativeNumber,partialDerivativeIndex,xi,err,error)
      
    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to evaluate 
    INTEGER(INTG), INTENT(IN) :: nodeNumber !<The local node number of the tensor product basis to evaluate
    INTEGER(INTG), INTENT(IN) :: derivativeNumber !<The local derivative number of the tensor product basis to evaluate
    INTEGER(INTG), INTENT(IN) :: partialDerivativeIndex  !<The partial derivative index to interpolate \see Constants_PartialDerivativeConstants
    REAL(DP), INTENT(IN) :: xi(:) !<The Xi position to evaluate the basis function at
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: Basis_LHTPBasisEvaluateDP !<On return the evaluated basis funtion.
    !Local variables
    INTEGER(INTG) :: xiIdx,localNodeIdx
    REAL(DP) :: sum
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Basis_LHTPBasisEvaluateDP",err,error,*999)
    
    Basis_LHTPBasisEvaluateDP=1.0_DP
    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(SIZE(xi,1)<basis%numberOfXi) THEN
      localError="The size of the supplied xi array of "//TRIM(NumberToVString(SIZE(xi,1),"*",err,error))// &
        & " is invalid. Basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))// &
        & " has "//TRIM(NumberToVString(basis%numberOfXi,"*",err,error))//" xi directions."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO xiIdx=1,basis%numberOfXi
      IF(basis%nodeAtCollapse(nodeNumber).AND.basis%collapsedXi(xiIdx)==BASIS_XI_COLLAPSED) THEN
        !We are at a collapsed node in the collapsed xi direction. Sum the basis functions in the collapsed xi direction.
        sum=0.0_DP
        SELECT CASE(basis%interpolationType(xiIdx))
        CASE(BASIS_LAGRANGE_INTERPOLATION)
          SELECT CASE(basis%interpolationOrder(xiIdx))
          CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
            DO localNodeIdx=1,2
              sum=sum+Lagrange_LinearEvaluate(localNodeIdx,PARTIAL_DERIVATIVE_INDEX(partialDerivativeIndex,xiIdx),xi(xiIdx), &
                & err,error)
            ENDDO !localNodeIdx
            IF(err/=0) GOTO 999
          CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
            DO localNodeIdx=1,3
              sum=sum+Lagrange_QuadraticEvaluate(localNodeIdx,PARTIAL_DERIVATIVE_INDEX(partialDerivativeIndex,xiIdx),xi(xiIdx), &
                & err,error)
            ENDDO !localNodeIdx
            IF(err/=0) GOTO 999
          CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
            DO localNodeIdx=1,4
              sum=sum+Lagrange_CubicEvaluate(localNodeIdx,PARTIAL_DERIVATIVE_INDEX(partialDerivativeIndex,xiIdx),xi(xiIdx), &
                & err,error)
            ENDDO !localNodeIdx
            IF(err/=0) GOTO 999
          CASE DEFAULT
            localError="Interpolation order value "//TRIM(NumberToVString(basis%interpolationOrder(xiIdx),"*",err,error))// &
              & " for xi direction "//TRIM(NumberToVString(xiIdx,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(BASIS_HERMITE_INTERPOLATION)
          SELECT CASE(basis%interpolationOrder(xiIdx))
          CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
            DO localNodeIdx=1,2
              sum=sum+Hermite_CubicEvaluate(localNodeIdx,basis%derivativeOrderIndex(derivativeNumber,nodeNumber,xiIdx), &
                & PARTIAL_DERIVATIVE_INDEX(partialDerivativeIndex,xiIdx),xi(xiIdx),err,error)
            ENDDO !localNodeIdx
            IF(err/=0) GOTO 999
          CASE DEFAULT
            localError="Interpolation order value "//TRIM(NumberToVString(basis%interpolationOrder(xiIdx),"*",err,error))// &
              & " for xi direction "//TRIM(NumberToVString(xiIdx,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="Interpolation type value "//TRIM(NumberToVString(basis%interpolationType(xiIdx),"*",err,error))// &
            & " for xi direction "//TRIM(NumberToVString(xiIdx,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        Basis_LHTPBasisEvaluateDP=Basis_LHTPBasisEvaluateDP*sum
      ELSE
        SELECT CASE(basis%interpolationType(xiIdx))
        CASE(BASIS_LAGRANGE_INTERPOLATION)
          SELECT CASE(basis%interpolationOrder(xiIdx))
          CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
            Basis_LHTPBasisEvaluateDP=Basis_LHTPBasisEvaluateDP* &
              & Lagrange_LinearEvaluate(basis%nodePositionIndex(nodeNumber,xiIdx), &
              & PARTIAL_DERIVATIVE_INDEX(partialDerivativeIndex,xiIdx),xi(xiIdx),err,error)
            IF(err/=0) GOTO 999
          CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
            Basis_LHTPBasisEvaluateDP=Basis_LHTPBasisEvaluateDP* &
              & Lagrange_QuadraticEvaluate(basis%nodePositionIndex(nodeNumber,xiIdx), &
              & PARTIAL_DERIVATIVE_INDEX(partialDerivativeIndex,xiIdx),xi(xiIdx),err,error)
            IF(err/=0) GOTO 999
          CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
            Basis_LHTPBasisEvaluateDP=Basis_LHTPBasisEvaluateDP* &
              & Lagrange_CubicEvaluate(basis%nodePositionIndex(nodeNumber,xiIdx), &
              & PARTIAL_DERIVATIVE_INDEX(partialDerivativeIndex,xiIdx),xi(xiIdx),err,error)
            IF(err/=0) GOTO 999
          CASE DEFAULT
            localError="Interpolation order value "//TRIM(NumberToVString(basis%interpolationOrder(xiIdx),"*",err,error))// &
              & " for xi direction "//TRIM(NumberToVString(xiIdx,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(BASIS_HERMITE_INTERPOLATION)
          SELECT CASE(basis%interpolationOrder(xiIdx))
          CASE(BASIS_QUADRATIC1_INTERPOLATION_ORDER)
            Basis_LHTPBasisEvaluateDP=Basis_LHTPBasisEvaluateDP* &
              & Hermite_QuadraticEvaluate(basis%nodePositionIndex(nodeNumber,xiIdx), &
              & basis%derivativeOrderIndex(derivativeNumber,nodeNumber,xiIdx), &
              & PARTIAL_DERIVATIVE_INDEX(partialDerivativeIndex,xiIdx),1,xi(xiIdx),err,error)
            IF(err/=0) GOTO 999
          CASE(BASIS_QUADRATIC2_INTERPOLATION_ORDER)
            Basis_LHTPBasisEvaluateDP=Basis_LHTPBasisEvaluateDP* &
              & Hermite_QuadraticEvaluate(basis%nodePositionIndex(nodeNumber,xiIdx), &
              & basis%derivativeOrderIndex(derivativeNumber,nodeNumber,xiIdx), &
              & PARTIAL_DERIVATIVE_INDEX(partialDerivativeIndex,xiIdx),2,xi(xiIdx),err,error)
            IF(err/=0) GOTO 999
          CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
            Basis_LHTPBasisEvaluateDP=Basis_LHTPBasisEvaluateDP* &
              & Hermite_CubicEvaluate(basis%nodePositionIndex(nodeNumber,xiIdx), &
              & basis%derivativeOrderIndex(derivativeNumber,nodeNumber,xiIdx), &
              & PARTIAL_DERIVATIVE_INDEX(partialDerivativeIndex,xiIdx),xi(xiIdx),err,error)
            IF(err/=0) GOTO 999
          CASE DEFAULT
            localError="Interpolation order value "//TRIM(NumberToVString(basis%interpolationOrder(xiIdx),"*",err,error))// &
              & " for xi direction "//TRIM(NumberToVString(xiIdx,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="Interpolation type value "//TRIM(NumberToVString(basis%interpolationType(xiIdx),"*",err,error))// &
            & " for xi direction "//TRIM(NumberToVString(xiIdx,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDIF
    ENDDO !xiIdx

    EXITS("Basis_LHTPBasisEvaluateDP")
    RETURN
999 ERRORSEXITS("Basis_LHTPBasisEvaluateDP",err,error)
    RETURN
    
  END FUNCTION Basis_LHTPBasisEvaluateDP

  !
  !================================================================================================================================
  !

  !>Creates and initialises a Lagrange-Hermite tensor product basis family that has already been allocated by Basis_CreateStart.
  !>\see BasisRoutines::Basis_LHTPCreate
  SUBROUTINE Basis_LHTPFamilyCreate(basis,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to create
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,faceXi(2),faceXi2(2),localFaceIdx,localLineIdx,xiIdx,xiIdx2
    LOGICAL :: lineBasisDone,faceBasisDone
    TYPE(BasisType), POINTER :: newSubBasis
    TYPE(VARYING_STRING) :: dummyError

    NULLIFY(newSubBasis)

    ENTERS("Basis_LHTPFamilyCreate",err,error,*999)

    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    
    !Create the main (parent) basis
    CALL Basis_LHTPBasisCreate(basis,err,error,*999)
    
    IF(basis%numberOfXi>1) THEN
      !Create the line bases as sub-basis types
      ALLOCATE(basis%lineBases(basis%numberOfXi),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate basis line bases.",err,error,*999)
      ALLOCATE(basis%localLineBasis(basis%numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate local line basis.",err,error,*999)
      DO xiIdx=1,basis%numberOfXi
        lineBasisDone=.FALSE.
        NULLIFY(newSubBasis)
        DO xiIdx2=1,xiIdx-1
          IF(basis%interpolationXi(xiIdx2)==basis%interpolationXi(xiIdx).AND. &
            basis%quadrature%numberOfGaussXi(xiIdx2)==basis%quadrature%numberOfGaussXi(xiIdx)) THEN
            lineBasisDone=.TRUE.
            EXIT
          ENDIF
        ENDDO !xiIdx2
        IF(lineBasisDone) THEN
          basis%lineBases(xiIdx)%ptr=>basis%lineBases(xiIdx2)%ptr
        ELSE
          !Create the new sub-basis
          CALL Basis_SubBasisCreate(basis,1,[xiIdx],newSubBasis,err,error,*999)
          !Fill in the basis information
          CALL Basis_LHTPBasisCreate(newSubBasis,err,error,*999)
          basis%lineBases(xiIdx)%ptr=>newSubBasis
        ENDIF
      ENDDO !xiIdx
      DO localLineIdx=1,basis%numberOfLocalLines
        xiIdx=basis%localLineXiDirection(localLineIdx)
        basis%localLineBasis(localLineIdx)%ptr=>basis%lineBases(xiIdx)%ptr
      ENDDO !localLineIdx
      IF(basis%numberOfXi>2) THEN
        !Create the face bases as sub-basis types
        ALLOCATE(basis%faceBases(basis%numberOfXi),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate basis face bases.",err,error,*999)
        ALLOCATE(basis%localFaceBasis(basis%numberOfLocalFaces),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate local face basis.",err,error,*999)
        DO xiIdx=1,basis%numberOfXi
          !Determine the face xi directions that lie in this xi direction
          faceXi(1)=OTHER_XI_DIRECTIONS3(xiIdx,2,1)
          faceXi(2)=OTHER_XI_DIRECTIONS3(xiIdx,3,1)
          faceBasisDone=.FALSE.
          NULLIFY(newSubBasis)
          DO xiIdx2=1,xiIdx-1
            !Determine the face xi directions that lie in this xi direction
            faceXi2(1)=OTHER_XI_DIRECTIONS3(xiIdx2,2,1)
            faceXi2(2)=OTHER_XI_DIRECTIONS3(xiIdx2,3,1)
            IF(basis%interpolationXi(faceXi2(1))==basis%interpolationXi(faceXi(1)).AND. &
              & basis%interpolationXi(faceXi2(2))==basis%interpolationXi(faceXi(2)).AND. &
              & basis%quadrature%numberOfGaussXi(faceXi2(1))==basis%quadrature%numberOfGaussXi(faceXi(1)).AND. &
              & basis%quadrature%numberOfGaussXi(faceXi2(2))==basis%quadrature%numberOfGaussXi(faceXi(2)).AND. &
              & basis%collapsedXi(faceXi2(1))==basis%collapsedXi(faceXi(1)).AND. &
              & basis%collapsedXi(faceXi2(2))==basis%collapsedXi(faceXi(2))) THEN              
              faceBasisDone=.TRUE.
              EXIT
            ENDIF
          ENDDO !xiIdx2
          IF(faceBasisDone) THEN
            basis%faceBases(xiIdx)%ptr=>basis%faceBases(xiIdx2)%ptr
          ELSE
            !Create the new sub-basis
            CALL Basis_SubBasisCreate(basis,2,[faceXi(1),faceXi(2)],newSubBasis,err,error,*999)
            !Fill in the basis information
            CALL Basis_LHTPBasisCreate(newSubBasis,err,error,*999)
            newSubBasis%lineBases(1)%ptr=>basis%lineBases(faceXi(1))%ptr
            newSubBasis%lineBases(2)%ptr=>basis%lineBases(faceXi(2))%ptr
            basis%faceBases(xiIdx)%ptr=>newSubBasis
            ALLOCATE(newSubBasis%localLineBasis(newSubBasis%numberOfLocalLines),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate new sub basis local line basis.",err,error,*999)
            DO localLineIdx=1,newSubBasis%numberOfLocalLines
              IF(newSubBasis%localLineXiDirection(localLineIdx)==1) THEN
                newSubBasis%localLineBasis(localLineIdx)%ptr=>newSubBasis%lineBases(1)%ptr
              ELSE IF(newSubBasis%localLineXiDirection(localLineIdx)==2) THEN
                newSubBasis%localLineBasis(localLineIdx)%ptr=>newSubBasis%lineBases(2)%ptr
              ENDIF
            ENDDO !localFaceIdx
          ENDIF
        ENDDO !xiIdx
        DO localFaceIdx=1,basis%numberOfLocalFaces
          xiIdx=basis%localFaceXiNormal(localFaceIdx)
          basis%localFaceBasis(localFaceIdx)%ptr=>basis%faceBases(ABS(xiIdx))%ptr
        ENDDO !localFaceIdx
      ENDIF
    ENDIF
    
    EXITS("Basis_LHTPFamilyCreate")
    RETURN
999 IF(ASSOCIATED(newSubBasis)) CALL Basis_FamilyDestroy(newSubBasis%basisFunctions,newSubBasis%userNumber, &
      & newSubBasis%familyNumber,dummyErr,dummyError,*998)
998 ERRORSEXITS("Basis_LHTPFamilyCreate",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_LHTPFamilyCreate

  !
  !================================================================================================================================
  !

  !>Creates and initialises a Radial basis family that has already been allocated by Basis_CreateStart.
  !>\see BasisRoutines::Basis_RadialCreate
  SUBROUTINE Basis_RadialFamilyCreate(basis,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to create
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Basis_RadialFamilyCreate",err,error,*999)

    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated",err,error,*999)
    
    !Create the main (parent) basis
    CALL FlagError("Not implemented.",err,error,*999)

    EXITS("Basis_RadialFamilyCreate")
    RETURN
999 ERRORSEXITS("Basis_RadialFamilyCreate",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_RadialFamilyCreate

  !
  !================================================================================================================================
  !

  !>Calculates the xi location of a local node in a basis.
  SUBROUTINE Basis_LocalNodeXiCalculate(basis,localNodeNumber,xi,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to calculate the xi for
    INTEGER(INTG), INTENT(IN) :: localNodeNumber !<The local node number to calculate the xi for
    REAL(DP), INTENT(OUT) :: xi(:) !<xi(xiIdx). On return, the xi position of the local node in the basis
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: xiIdx
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Basis_LocalNodeXiCalculate",err,error,*999)

    CALL Basis_AssertIsFinished(basis,err,error,*999)   
    IF(localNodeNumber<1.OR.localNodeNumber>basis%numberOfNodes) THEN
      localError="The specified local node number of "//TRIM(NumberToVString(localNodeNumber,"*",err,error))// &
        & " is invalid. The local node number must be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfNodes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF      
    IF(SIZE(xi,1)>=basis%numberOfXi) THEN
      localError="The size of the specified xic array of "//TRIM(NumberToVString(SIZE(xi,1),"*",err,error))// &
        & " is invalid. The size of the xi array must be >= "// &
        & TRIM(NumberToVString(basis%numberOfXi,"*",err,error))//"."            
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(basis%type)
    CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
      DO xiIdx=1,basis%numberOfXi
        xi(xiIdx)=REAL(basis%nodePositionIndex(localNodeNumber,xiIdx)-1,DP)/REAL(basis%numberOfNodesXiC(xiIdx)-1,DP)
      ENDDO !xiIdx
    CASE(BASIS_SIMPLEX_TYPE)
      DO xiIdx=1,basis%numberOfXi
        xi(xiIdx)=REAL(basis%numberOfNodesXiC(xiIdx)-basis%nodePositionIndex(localNodeNumber,xiIdx),DP)/ &
          & REAL(basis%numberOfNodesXiC(xiIdx)-1,DP)
      ENDDO !xiIdx
    CASE(BASIS_SERENDIPITY_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(BASIS_AUXILLIARY_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(BASIS_B_SPLINE_TP_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(BASIS_EXTENDED_LAGRANGE_TP_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The basis type of "//TRIM(NumberToVString(basis%type,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Basis_LocalNodeXiCalculate")
    RETURN
999 ERRORSEXITS("Basis_LocalNodeXiCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_LocalNodeXiCalculate

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of xi directions for a basis. \see OpenCMISS::Iron::cmfe_Basis_NumberOfXiSet
  SUBROUTINE Basis_NumberOfXiSet(basis,numberOfXi,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis function to set the number of xi for
    INTEGER(INTG), INTENT(IN) :: numberOfXi !<The number of Xi directions to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG), ALLOCATABLE :: newCollapsedXi(:),newInterpolationXi(:),newNumberOfGaussXi(:)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Basis_NumberOfXiSet",err,error,*999)

    CALL Basis_AssertNotFinished(basis,err,error,*999)
    IF(numberOfXi<1.OR.numberOfXi>3) THEN
      localError="The specified number of xi directions of "//TRIM(NumberToVString(numberOfXi,"*",err,error))// &
        & " is invalid. The number of xi directions must be >= 1 and <=3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    IF(basis%numberOfXi/=numberOfXi) THEN     
      !Reallocate the basis information arrays that depend on the number of xi directions
      ALLOCATE(newInterpolationXi(numberOfXi),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new interpolation xi.",err,error,*999)
      ALLOCATE(newCollapsedXi(numberOfXi),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new COLLAPSED xi.",err,error,*999)
      IF(numberOfXi>basis%numberOfXi) THEN
        newInterpolationXi(1:basis%numberOfXi)=basis%interpolationXi(1:basis%numberOfXi)
        newInterpolationXi(basis%numberOfXi+1:numberOfXi)=basis%interpolationXi(1)
        newCollapsedXi(1:basis%numberOfXi)=basis%collapsedXi(1:basis%numberOfXi)
        newCollapsedXi(basis%numberOfXi+1:numberOfXi)=basis%collapsedXi(1)
      ELSE
        newInterpolationXi(1:numberOfXi)=basis%interpolationXi(1:numberOfXi)
        newCollapsedXi(1:numberOfXi)=basis%collapsedXi(1:numberOfXi)
      ENDIF
      CALL MOVE_ALLOC(newInterpolationXi,basis%interpolationXi)
      CALL MOVE_ALLOC(newCollapsedXi,basis%collapsedXi)
      IF(ASSOCIATED(basis%quadrature%basis)) THEN
        ALLOCATE(newNumberOfGaussXi(numberOfXi),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new number of Gauss xi.",err,error,*999)
        IF(numberOfXi>basis%numberOfXi) THEN
          newNumberOfGaussXi(1:basis%numberOfXi)=basis%quadrature%numberOfGaussXi(1:basis%numberOfXi)
          newNumberOfGaussXi(basis%numberOfXi+1:numberOfXi)=basis%quadrature%numberOfGaussXi(1)
        ELSE
          newNumberOfGaussXi(1:numberOfXi)=basis%quadrature%numberOfGaussXi(1:numberOfXi)
        ENDIF
        CALL MOVE_ALLOC(newNumberOfGaussXi,basis%quadrature%numberOfGaussXi)        
      ENDIF
      basis%numberOfXi=numberOfXi
    ENDIF
    
    EXITS("Basis_NumberOfXiSet")
    RETURN
999 IF(ALLOCATED(newInterpolationXi)) DEALLOCATE(newInterpolationXi)
    IF(ALLOCATED(newCollapsedXi)) DEALLOCATE(newCollapsedXi)
    IF(ALLOCATED(newNumberOfGaussXi)) DEALLOCATE(newNumberOfGaussXi)
    ERRORSEXITS("Basis_NumberOfXiSet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_NumberOfXiSet
  
  !
  !================================================================================================================================
  !

  !>Creates the quadrature and quadrature schemes on a basis.
  SUBROUTINE Basis_QuadratureCreate(basis,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: schemeIdx,i,j,k,maximumNumberOfGauss,gaussPointIdx,xiIdx,derivativeIdx,localNodeIdx,elementParameterIdx, &
      & partialDerivativeIdx,numberOfGauss1,numberOfGauss2,numberOfGauss3
    REAL(DP) :: xi(3),gsx(4,20),gsw(20)
    REAL(DP), ALLOCATABLE :: positions(:,:),postitionsMatrix(:,:,:,:),weights(:,:)
    TYPE(QuadratureSchemeType), POINTER :: newScheme,scheme
    TYPE(QuadratureSchemePtrType), POINTER :: newSchemes(:)
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: maximumNumberOfFaceGauss,faceIdx,normal,faceXi(2),numberOfFaceXiCoordinates

    NULLIFY(newScheme)
    NULLIFY(newSchemes)

    ENTERS("Basis_QuadratureCreate",err,error,*999)

    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*998)
    IF(ASSOCIATED(basis%quadrature%schemes)) THEN
      localError="The quadrature schemes on basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))// &
        & " are already associated"
      CALL FlagError(localError,err,error,*998)
    ENDIF
     
!!TODO: \todo Sort this properly by having a create values cache.
    !Reset the basis quadrature - 
    !CALL Basis_QuadratureFinalise(basis,err,error,*999)      ! Kumar - I think this is not correct as it 
    !Initialise the basis quadrature                           !         resets the quadrature scheme already set.
    !CALL Basis_QuadratureInitialise(basis,err,error,*999)    ! 
    SELECT CASE(basis%quadrature%TYPE)
    CASE(BASIS_GAUSS_LEGENDRE_QUADRATURE)            
      !Allocate one scheme and add it to the list of schemes
      ALLOCATE(newScheme,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new quadrature scheme.",err,error,*999)
      newScheme%quadrature=>basis%quadrature
      basis%quadrature%numberOfSchemes=1
      ALLOCATE(newSchemes(basis%quadrature%numberOfSchemes),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new quadratures scheme.",err,error,*999)
      newSchemes(1)%ptr=>newScheme
      basis%quadrature%schemes=>newSchemes
      !Set up the quadrature scheme map
      ALLOCATE(basis%quadrature%quadratureSchemeMap(BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate quadrature scheme map.",err,error,*999)
      DO schemeIdx=1,BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES
        NULLIFY(basis%quadrature%quadratureSchemeMap(schemeIdx)%ptr)
      ENDDO !schemeIdx
      basis%quadrature%quadratureSchemeMap(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr=>newScheme
      !Set up the gauss point arrays
      newScheme%numberOfGauss=1
      maximumNumberOfGauss=-1
      DO xiIdx=1,basis%numberOfXi
        newScheme%numberOfGauss=newScheme%numberOfGauss*basis%quadrature%numberOfGaussXi(xiIdx)
        IF(basis%quadrature%numberOfGaussXi(xiIdx)>maximumNumberOfGauss) &
          & maximumNumberOfGauss=basis%quadrature%numberOfGaussXi(xiIdx)
      ENDDO !xiIdx
      ALLOCATE(newScheme%gaussPositions(basis%numberOfXiCoordinates,newScheme%numberOfGauss),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate Gauss positions.",err,error,*999)
      ALLOCATE(newScheme%gaussWeights(newScheme%numberOfGauss),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate Gauss weights.",err,error,*999)
      ALLOCATE(newScheme%gaussBasisFunctions(basis%numberOfElementParameters,basis%numberOfPartialDerivatives, &
        & newScheme%numberOfGauss),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate Gauss basis functions.",err,error,*999)
      ALLOCATE(weights(maximumNumberOfGauss,3),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate weights.",err,error,*999)
      ALLOCATE(positions(maximumNumberOfGauss,3),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate positions.",err,error,*999)
      ALLOCATE(postitionsMatrix(maximumNumberOfGauss,maximumNumberOfGauss,maximumNumberOfGauss,3),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate positions matrix.",err,error,*999)
      weights=1.0_DP
      positions=0.0_DP
      postitionsMatrix=0.0_DP
      DO xiIdx=1,basis%numberOfXi
        CALL Gauss_Legendre(basis%quadrature%numberOfGaussXi(xiIdx),0.0_DP,1.0_DP, &
          & positions(1:basis%quadrature%numberOfGaussXi(xiIdx),xiIdx), &
          & weights(1:basis%quadrature%numberOfGaussXi(xiIdx),xiIdx),err,error,*999)
      ENDDO !xiIdx
      SELECT CASE(basis%numberOfXi)
      CASE(1)
        numberOfGauss1=basis%quadrature%numberOfGaussXi(1)
        numberOfGauss2=1
        numberOfGauss3=1
      CASE(2)
        numberOfGauss1=basis%quadrature%numberOfGaussXi(1)
        numberOfGauss2=basis%quadrature%numberOfGaussXi(2)
        numberOfGauss3=1
      CASE(3)
        numberOfGauss1=basis%quadrature%numberOfGaussXi(1)
        numberOfGauss2=basis%quadrature%numberOfGaussXi(2)
        numberOfGauss3=basis%quadrature%numberOfGaussXi(3)
      CASE DEFAULT
        localError="The number of basis xi directions of "//TRIM(NumberToVString(basis%numberOfXi,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      DO k=1,numberOfGauss3
        DO j=1,numberOfGauss2
          DO i=1,numberOfGauss1
            postitionsMatrix(i,j,k,1)=positions(i,1)
            postitionsMatrix(i,j,k,2)=positions(j,2)
            postitionsMatrix(i,j,k,3)=positions(k,3)
            xi(1:basis%numberOfXi)=postitionsMatrix(i,j,k,1:basis%numberOfXi)
            gaussPointIdx=i+(j-1+(k-1)*numberOfGauss2)*numberOfGauss1
            newScheme%gaussWeights(gaussPointIdx)=weights(i,1)*weights(j,2)*weights(k,3)
            newScheme%gaussPositions(1:basis%numberOfXiCoordinates,gaussPointIdx)=xi(1:basis%numberOfXiCoordinates)
            elementParameterIdx=0
            DO localNodeIdx=1,basis%numberOfNodes
              DO derivativeIdx=1,basis%numberOfDerivatives(localNodeIdx)
                elementParameterIdx=elementParameterIdx+1
                DO partialDerivativeIdx=1,basis%numberOfPartialDerivatives
                  SELECT CASE(basis%type)
                  CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
                    newScheme%gaussBasisFunctions(elementParameterIdx,partialDerivativeIdx,gaussPointIdx)= &
                      & Basis_LHTPBasisEvaluate(basis,localNodeIdx,derivativeIdx,partialDerivativeIdx,xi,err,error)
                    IF(err/=0) GOTO 999                        
                  CASE DEFAULT
                    CALL FlagError("Not implemented.",err,error,*999)
                  END SELECT
                ENDDO !partialDerivativeIdx
              ENDDO !derivativeIdx
            ENDDO !localNodeIdx
          ENDDO !i
        ENDDO !j
      ENDDO !k
      !Create face quadrature scheme, if requested
      IF(basis%quadrature%evaluateFaceGauss) THEN
        IF(basis%numberOfXi/=3) &
          & CALL FlagError("Can not evaluate face quadrature schemes for a non three dimensional element.",err,error,*999)       
        !Find maximum number of face gauss points and allocate the arrays
        maximumNumberOfFaceGauss=PRODUCT(basis%quadrature%numberOfGaussXi(1:basis%numberOfXi))
        maximumNumberOfFaceGauss=maximumNumberOfFaceGauss/MINVAL(basis%quadrature%numberOfGaussXi(1:basis%numberOfXi))
        ALLOCATE(newScheme%numberOfFaceGauss(basis%numberOfLocalFaces),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate number of face gauss.",err,error,*999)
        ALLOCATE(newScheme%faceGaussPositions(basis%numberOfXiCoordinates,maximumNumberOfFaceGauss, &
          & basis%numberOfLocalFaces),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate face Gauss positions.",err,error,*999)
        ALLOCATE(newScheme%faceGaussWeights(maximumNumberOfFaceGauss,basis%numberOfLocalFaces),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate face Gauss weights.",err,error,*999)
        ALLOCATE(newScheme%faceGaussBasisFunctions(basis%numberOfElementParameters,basis%numberOfPartialDerivatives, &
          & maximumNumberOfFaceGauss,basis%numberOfLocalFaces),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate face Gauss basis function values array.",err,error,*999)
        !Zero them out just to be safe
        newScheme%faceGaussPositions=0.0_DP
        newScheme%faceGaussWeights=0.0_DP
        newScheme%faceGaussBasisFunctions=0.0_DP
        !Populate face_gauss_positions, weights, basis_fn
        DO faceIdx=1,basis%numberOfLocalFaces
          !What's the normal?
          normal=basis%localFaceXiNormal(faceIdx)
          IF(normal<0_INTG) THEN
            xi(ABS(normal))=0.0_DP
          ELSE
            xi(ABS(normal))=1.0_DP
          ENDIF
          normal=ABS(normal)
          faceXi=[OTHER_XI_DIRECTIONS3(normal,2,1), OTHER_XI_DIRECTIONS3(normal,3,1)]
          !How many gauss points are in this face?
          newScheme%numberOfFaceGauss(faceIdx)=PRODUCT(basis%quadrature%numberOfGaussXi(faceXi))
          gaussPointIdx=0_INTG
          DO j=1,basis%quadrature%numberOfGaussXi(faceXi(2))
            xi(faceXi(2))=positions(j,faceXi(2))
            DO i=1,basis%quadrature%numberOfGaussXi(faceXi(1))
              xi(faceXi(1))=positions(i,faceXi(1))
              gaussPointIdx=gaussPointIdx+1_INTG
              !Gauss point xi and weights first
              newScheme%faceGaussWeights(gaussPointIdx,faceIdx)=weights(i,faceXi(1))*weights(j,faceXi(2))
              newScheme%faceGaussPositions(1:3,gaussPointIdx,faceIdx)=xi(1:3)
              !Evaluate basis fn values at the Gauss points now
              elementParameterIdx=0
              DO localNodeIdx=1,basis%numberOfNodes
                DO derivativeIdx=1,basis%numberOfDerivatives(localNodeIdx)
                  elementParameterIdx=elementParameterIdx+1
                  DO partialDerivativeIdx=1,basis%numberOfPartialDerivatives
                    SELECT CASE(basis%TYPE)
                    CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
                      newScheme%faceGaussBasisFunctions(elementParameterIdx,partialDerivativeIdx,gaussPointIdx,faceIdx)= &
                        & Basis_LHTPBasisEvaluate(basis,localNodeIdx,derivativeIdx,partialDerivativeIdx,xi,err,error)
                      IF(err/=0) GOTO 999                        
                    CASE DEFAULT
                      CALL FlagError("Not implemented.",err,error,*999)
                    END SELECT
                  ENDDO !partialDerivativeIdx
                ENDDO !derivativeIdx
              ENDDO !localNodeIdx              
            ENDDO !i
          ENDDO !j
        ENDDO !faceIdx
      ENDIF
      !Clean up
      DEALLOCATE(weights)
      DEALLOCATE(positions)
      DEALLOCATE(postitionsMatrix)
    CASE(BASIS_GAUSS_LAGUERRE_QUADRATURE)
      CALL FlagError("Gauss Laguerre quadrature type not implemented.",err,error,*999)
    CASE(BASIS_GUASS_HERMITE_QUADRATURE)
      CALL FlagError("Gauss Hermite quadrature type not implemented.",err,error,*999)
    CASE(BASIS_GAUSS_SIMPLEX_QUADRATURE)
      !Allocate one scheme and add it to the list of schemes
      ALLOCATE(newScheme,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new quadrature scheme.",err,error,*999)
      newScheme%quadrature=>basis%quadrature
      basis%quadrature%numberOfSchemes=1
      ALLOCATE(newSchemes(basis%quadrature%numberOfSchemes),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new quadratures scheme.",err,error,*999)
      newSchemes(1)%ptr=>newScheme
      basis%quadrature%schemes=>newSchemes
      !Set up the quadrature scheme map
      ALLOCATE(basis%quadrature%quadratureSchemeMap(BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate quadrature scheme map.",err,error,*999)
      DO schemeIdx=1,BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES
        NULLIFY(basis%quadrature%quadratureSchemeMap(schemeIdx)%ptr)
      ENDDO !schemeIdx
      basis%quadrature%quadratureSchemeMap(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr=>newScheme
      !Set up the gauss point arrays
      CALL Gauss_Simplex(basis%quadrature%gaussOrder,basis%numberOfXiCoordinates,newScheme%numberOfGauss,gsx,gsw, &
        & err,error,*999)
      ALLOCATE(newScheme%gaussPositions(basis%numberOfXiCoordinates,newScheme%numberOfGauss),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate Gauss positions.",err,error,*999)
      ALLOCATE(newScheme%gaussWeights(newScheme%numberOfGauss),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate Gauss weights.",err,error,*999)
      ALLOCATE(newScheme%gaussBasisFunctions(basis%numberOfElementParameters,basis%numberOfPartialDerivatives, &
        & newScheme%numberOfGauss),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate Gauss basis functions.",err,error,*999)
      newScheme%gaussPositions(1:basis%numberOfXiCoordinates,1:newScheme%numberOfGauss)= &
        & gsx(1:basis%numberOfXiCoordinates,1:newScheme%numberOfGauss)
      newScheme%gaussWeights(1:newScheme%numberOfGauss)=gsw(1:newScheme%numberOfGauss)
      DO gaussPointIdx=1,newScheme%numberOfGauss
        elementParameterIdx=0
        DO localNodeIdx=1,basis%numberOfNodes
          DO derivativeIdx=1,basis%numberOfDerivatives(localNodeIdx)
            elementParameterIdx=elementParameterIdx+1
            DO partialDerivativeIdx=1,basis%numberOfPartialDerivatives
              SELECT CASE(basis%TYPE)
              CASE(BASIS_SIMPLEX_TYPE)
                !Gauss positions are in area coordinates so call the simplex basis evaluate directly
                newScheme%gaussBasisFunctions(elementParameterIdx,partialDerivativeIdx,gaussPointIdx)= &
                  & Basis_SimplexBasisEvaluate(basis,localNodeIdx,partialDerivativeIdx, &
                  & newScheme%gaussPositions(1:basis%numberOfXiCoordinates,gaussPointIdx),err,error)
                IF(err/=0) GOTO 999                        
              CASE DEFAULT
                CALL FlagError("Not implemented.",err,error,*999)
              END SELECT
            ENDDO !partialDerivativeIdx
          ENDDO !derivativeIdx
        ENDDO !localNodeIdx
      ENDDO !gaussPointIdx
      !Create face quadrature scheme, if requested
      IF(basis%quadrature%evaluateFaceGauss) THEN
        IF(basis%numberOfXi/=3) &
          & CALL FlagError("Cannot evaluate face quadrature schemes for a non three dimensional element.",err,error,*999)
        !Find maximum number of face gauss points and allocate the arrays
        maximumNumberOfFaceGauss=PRODUCT(basis%quadrature%numberOfGaussXi(1:basis%numberOfXi))
        maximumNumberOfFaceGauss=maximumNumberOfFaceGauss/MINVAL(basis%quadrature%numberOfGaussXi(1:basis%numberOfXi))
        ALLOCATE(newScheme%numberOfFaceGauss(basis%numberOfLocalFaces),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate number of face gauss.",err,error,*999)
        ALLOCATE(newScheme%faceGaussPositions(basis%numberOfXiCoordinates,maximumNumberOfFaceGauss, &
          & basis%numberOfLocalFaces),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate face Gauss positions.",err,error,*999)
        ALLOCATE(newScheme%faceGaussWeights(maximumNumberOfFaceGauss,basis%numberOfLocalFaces),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate face Gauss weights.",err,error,*999)
        ALLOCATE(newScheme%faceGaussBasisFunctions(basis%numberOfElementParameters,basis%numberOfPartialDerivatives, &
          & maximumNumberOfFaceGauss,basis%numberOfLocalFaces),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate face Gauss basis function values array.",err,error,*999)
        !Zero them out just to be safe
        newScheme%faceGaussPositions=0.0_DP
        newScheme%faceGaussWeights=0.0_DP
        newScheme%faceGaussBasisFunctions=0.0_DP
        !Populate face_gauss_positions, weights, basis_fn
        DO faceIdx=1,basis%numberOfLocalFaces
          !The number of face xi coordinates will be 3 for triangular face on a tet
          numberOfFaceXiCoordinates = basis%numberOfXi
          !Set up the gauss point arrays for the face
          CALL Gauss_Simplex(basis%quadrature%gaussOrder,numberOfFaceXiCoordinates, &
            & newScheme%numberOfFaceGauss(faceIdx),gsx,gsw,err,error,*999)
          IF(err/=0) CALL FlagError("Could not allocate Gauss basis functions.",err,error,*999)
          newScheme%faceGaussPositions(1:numberOfFaceXiCoordinates,1:newScheme%numberOfFaceGauss(faceIdx), &
            & faceIdx)=gsx(1:numberOfFaceXiCoordinates,1:newScheme%numberOfFaceGauss(faceIdx))
          newScheme%faceGaussWeights(1:newScheme%numberOfFaceGauss(faceIdx),faceIdx)= &
            & gsw(1:newScheme%numberOfFaceGauss(faceIdx))
          
          DO gaussPointIdx=1,newScheme%numberOfFaceGauss(faceIdx)
            elementParameterIdx=0
            DO localNodeIdx=1,basis%numberOfNodesInLocalFace(faceIdx)
              DO derivativeIdx=1,basis%numberOfDerivatives(localNodeIdx)
                elementParameterIdx=elementParameterIdx+1
                DO partialDerivativeIdx=1,basis%numberOfPartialDerivatives
                  SELECT CASE(basis%TYPE)
                  CASE(BASIS_SIMPLEX_TYPE)
                    newScheme%faceGaussBasisFunctions(elementParameterIdx,partialDerivativeIdx,gaussPointIdx,faceIdx)= &
                      & Basis_SimplexBasisEvaluate(basis,localNodeIdx,partialDerivativeIdx, &
                      & newScheme%faceGaussPositions(1:numberOfFaceXiCoordinates,gaussPointIdx,faceIdx),err,error)
                    IF(err/=0) GOTO 999                        
                  CASE DEFAULT
                    CALL FlagError("Not implemented",err,error,*999)
                  END SELECT
                ENDDO !partialDerivativeIdx
              ENDDO !derivativeIdx
            ENDDO !localNodeIdx
          ENDDO !gaussPointIdx          
        ENDDO !faceIdx
      ENDIF
    CASE DEFAULT
      localError="Quadrature type "//TRIM(NumberToVString(basis%quadrature%TYPE,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    IF(diagnostics1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Quadrature type = ",basis%quadrature%type,err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%numberOfXi,3,3,basis%quadrature%numberOfGaussXi, &
        & '("  Number of gauss points(xiIdx):",3(X,I2))','(22X,3(X,I2))',err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of quadrature schemes = ",basis%quadrature%numberOfSchemes, &
        & err,error,*999)
      DO schemeIdx=1,basis%quadrature%numberOfSchemes
        scheme=>basis%quadrature%schemes(schemeIdx)%ptr
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Scheme = ",schemeIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Total number of gauss points = ",scheme%numberOfGauss, &
          & err,error,*999)
        IF(diagnostics2) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      Gauss point positions and weights:",err,error,*999)
          DO gaussPointIdx=1,scheme%numberOfGauss
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Gauss point = ",gaussPointIdx,err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%numberOfXiCoordinates,3,3, &
              & scheme%gaussPositions(:,gaussPointIdx), &
              & '("          position(xiIdx)   :",3(X,F12.4))','(26X,3(X,F12.4))',err,error,*999)
            CALL WriteStringFmtValue(DIAGNOSTIC_OUTPUT_TYPE,"          WEIGHT         : ",scheme%gaussWeights(gaussPointIdx), &
              & "(F12.4)",err,error,*999)
          ENDDO !gaussPointIdx          
        ENDIF
        IF(diagnostics3) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      Basis functions evaluated at Gauss points:",err,error,*999)
          DO gaussPointIdx=1,scheme%numberOfGauss
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Gauss point = ",gaussPointIdx,err,error,*999)
            DO partialDerivativeIdx=1,basis%numberOfPartialDerivatives
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Partial derivative number = ",partialDerivativeIdx, &
                & err,error,*999)
              CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%numberOfElementParameters,4,4, &
                & scheme%gaussBasisFunctions(:,partialDerivativeIdx,gaussPointIdx), &
                & '("          BASIS FNS(elementParameterIdx)  :",4(X,F12.4))','(26X,4(X,F12.4))',err,error,*999)
            ENDDO !partialDerivativeIdx
          ENDDO !gaussPointIdx
        ENDIF
      ENDDO !schemeIdx
    ENDIF
    
    EXITS("Basis_QuadratureCreate")
    RETURN
999 IF(ASSOCIATED(newScheme)) THEN
      IF(ALLOCATED(newScheme%gaussPositions)) DEALLOCATE(newScheme%gaussPositions)
      IF(ALLOCATED(newScheme%gaussWeights)) DEALLOCATE(newScheme%gaussWeights)
      IF(ALLOCATED(newScheme%gaussBasisFunctions)) DEALLOCATE(newScheme%gaussBasisFunctions)
      DEALLOCATE(newScheme)     
    ENDIF
    IF(ALLOCATED(weights)) DEALLOCATE(weights)
    IF(ALLOCATED(positions)) DEALLOCATE(positions)
    IF(ALLOCATED(postitionsMatrix)) DEALLOCATE(postitionsMatrix)
    IF(ASSOCIATED(newSchemes)) DEALLOCATE(newSchemes)
    NULLIFY(basis%quadrature%schemes)
998 ERRORSEXITS("Basis_QuadratureCreate",err,error)    
    RETURN 1
   
  END SUBROUTINE Basis_QuadratureCreate
        
  !
  !================================================================================================================================
  !
  
  !>Destroys a quadrature on a given basis and deallocates all memory. \todo fix all this basis/quadrature into standard form.
  SUBROUTINE Basis_QuadratureDestroy(quadrature,err,error,*)

    !Argument variables
    TYPE(QuadratureType), POINTER :: quadrature !<A pointer to the quadrature
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Basis_QuadratureDestroy",err,error,*999)

    IF(.NOT.ASSOCIATED(quadrature)) CALL FlagError("Basis is not associated",err,error,*999)
    
    CALL FlagError("Not implemented.",err,error,*999)     
   
    EXITS("Basis_QuadratureDestroy")
    RETURN
999 ERRORSEXITS("Basis_QuadratureDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_QuadratureDestroy

  !
  !================================================================================================================================
  !
    
  !>Finalises a quadrature on a given basis and deallocates all memory
  SUBROUTINE Basis_QuadratureFinalise(basis,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: schemeIdx
    TYPE(QuadratureSchemeType), POINTER :: scheme

    ENTERS("Basis_QuadratureFinalise",err,error,*999)

    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)    
    IF(.NOT.ASSOCIATED(basis%quadrature%basis)) CALL FlagError("Basis quadrature basis is not associated.",err,error,*999)
    
    DO schemeIdx=1,basis%quadrature%numberOfSchemes
      scheme=>basis%quadrature%schemes(schemeIdx)%ptr
      !Destroy all scheme components
      IF (ASSOCIATED(scheme)) THEN
        IF(ALLOCATED(scheme%gaussPositions)) DEALLOCATE(scheme%gaussPositions)
        IF(ALLOCATED(scheme%gaussWeights)) DEALLOCATE(scheme%gaussWeights)
        IF(ALLOCATED(scheme%gaussBasisFunctions)) DEALLOCATE(scheme%gaussBasisFunctions)
        DEALLOCATE(scheme)
      ENDIF
    ENDDO !schemeIdx
    IF(ASSOCIATED(basis%quadrature%schemes)) DEALLOCATE(basis%quadrature%schemes)
    basis%quadrature%numberOfSchemes=0
    IF(ALLOCATED(basis%quadrature%numberOfGaussXi)) DEALLOCATE(basis%quadrature%numberOfGaussXi)
    NULLIFY(basis%quadrature%basis)
    basis%quadrature%TYPE=-1
    IF(ALLOCATED(basis%quadrature%quadratureSchemeMap)) DEALLOCATE(basis%quadrature%quadratureSchemeMap)
   
    EXITS("Basis_QuadratureFinalise")
    RETURN
999 ERRORSEXITS("Basis_QuadratureFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_QuadratureFinalise

  !
  !================================================================================================================================
  !

  !>Initialises a quadrature on the given basis.
  SUBROUTINE Basis_QuadratureInitialise(basis,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: xiIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("Basis_QuadratureInitialise",err,error,*999)

    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*998)
    IF(ASSOCIATED(basis%quadrature%basis)) THEN
      localError="Basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))// &
        & " already has a quadrature associated."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    
    SELECT CASE(basis%TYPE)
    CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
      !Set up a default Gauss Legendre quadrature 
      basis%quadrature%numberOfSchemes=0
      NULLIFY(basis%quadrature%schemes)
      basis%quadrature%basis=>basis        
      basis%quadrature%TYPE=BASIS_GAUSS_LEGENDRE_QUADRATURE
      !Set up a default number of Gauss points appropriate for the given interpolation order in each direction.
      ALLOCATE(basis%quadrature%numberOfGaussXi(basis%numberOfXi),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate number of Gauss in each xi direction.",err,error,*999)
      DO xiIdx=1,basis%numberOfXi
        SELECT CASE(basis%interpolationXi(xiIdx))
        CASE(BASIS_LINEAR_LAGRANGE_INTERPOLATION)
          basis%quadrature%numberOfGaussXi(xiIdx)=2
        CASE(BASIS_QUADRATIC_LAGRANGE_INTERPOLATION)
          basis%quadrature%numberOfGaussXi(xiIdx)=3
        CASE(BASIS_CUBIC_LAGRANGE_INTERPOLATION)
          basis%quadrature%numberOfGaussXi(xiIdx)=4
        CASE(BASIS_QUADRATIC1_HERMITE_INTERPOLATION,BASIS_QUADRATIC2_HERMITE_INTERPOLATION)
          basis%quadrature%numberOfGaussXi(xiIdx)=3
        CASE(BASIS_CUBIC_HERMITE_INTERPOLATION)
          basis%quadrature%numberOfGaussXi(xiIdx)=4
        CASE DEFAULT
          localError="Interpolation xi value "//TRIM(NumberToVString(basis%interpolationXi(xiIdx),"*",err,error))// &
            & " in xi direction "//TRIM(NumberToVString(xiIdx,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !xiIdx
    CASE(BASIS_SIMPLEX_TYPE)
      !Set up a default quadrature 
      basis%quadrature%numberOfSchemes=0
      NULLIFY(basis%quadrature%schemes)
      basis%quadrature%basis=>basis        
      basis%quadrature%TYPE=BASIS_GAUSS_SIMPLEX_QUADRATURE
      !Set up a default order appropriate for the given interpolation.
      ALLOCATE(basis%quadrature%numberOfGaussXi(basis%numberOfXi),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate number of Gauss in each xi direction.",err,error,*999)
!!TODO: \todo Set these to something more meaningfull!
      SELECT CASE(basis%interpolationXi(1))
      CASE(BASIS_LINEAR_SIMPLEX_INTERPOLATION)
        SELECT CASE(basis%numberOfXi)
        CASE(1)
          basis%quadrature%gaussOrder=2
        CASE(2)
          basis%quadrature%gaussOrder=3
        CASE(3)
          basis%quadrature%gaussOrder=3
        CASE DEFAULT
          localError="The number of xi directions of "//TRIM(NumberToVString(basis%numberOfXi,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(BASIS_QUADRATIC_SIMPLEX_INTERPOLATION)
        SELECT CASE(basis%numberOfXi)
        CASE(1)
          basis%quadrature%gaussOrder=3
        CASE(2)
          basis%quadrature%gaussOrder=3
        CASE(3)
          basis%quadrature%gaussOrder=5
        CASE DEFAULT
          localError="The number of xi directions "//TRIM(NumberToVString(basis%numberOfXi,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(BASIS_CUBIC_SIMPLEX_INTERPOLATION)
        SELECT CASE(basis%numberOfXi)
        CASE(1)
          basis%quadrature%gaussOrder=3
        CASE(2)
          basis%quadrature%gaussOrder=5
        CASE(3)
          basis%quadrature%gaussOrder=5
        CASE DEFAULT
          localError="The number of xi directions of "//TRIM(NumberToVString(basis%numberOfXi,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="Interpolation xi value "//TRIM(NumberToVString(basis%interpolationXi(1),"*",err,error))// &
          & " in xi direction 1 is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      basis%quadrature%numberOfGaussXi=basis%quadrature%gaussOrder
    CASE DEFAULT
      localError="Basis type value "//TRIM(NumberToVString(basis%interpolationXi(xiIdx),"*",err,error))// &
        & " is invalid or not implemented."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("Basis_QuadratureInitialise")
    RETURN
999 IF(ALLOCATED(basis%quadrature%numberOfGaussXi)) DEALLOCATE(basis%quadrature%numberOfGaussXi)
998 ERRORSEXITS("Basis_QuadratureInitialise",err,error)    
    RETURN 1
    
  END SUBROUTINE Basis_QuadratureInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of Gauss points in each xi direction on a basis quadrature identified by a pointer. \see OpenCMISS::Iron::cmfe_Basis_QuadratureNumberOfGaussSet
  SUBROUTINE Basis_QuadratureNumberOfGaussXiSet(basis,numberOfGaussXi,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: numberOfGaussXi(:) !numberOfGaussXi(xiIdx). <The number of Gauss in each Xi direction
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: xiIdx
    TYPE(VARYING_STRING) :: localError,localWarning
    
    ENTERS("Basis_QuadratureNumberOfGaussXiSet",err,error,*999)

    CALL Basis_AssertNotFinished(basis,err,error,*999)
    IF(.NOT.ASSOCIATED(basis%quadrature%basis)) CALL FlagError("Quadrature basis is not associated.",err,error,*999)          
    IF(SIZE(numberOfGaussXi,1)/=basis%numberOfXi) THEN
      localError="The size of the number of Gauss array ("// &
        & TRIM(NumberToVString(SIZE(numberOfGaussXi,1),"*",err,error))// &
        & ") does not match the number of xi directions ("// &
        & TRIM(NumberToVString(basis%numberOfXi,"*",err,error))//") for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(ANY(numberOfGaussXi<1)) CALL FlagError("Invalid number of Gauss values.",err,error,*999)
    
    basis%quadrature%numberOfGaussXi(1:basis%numberOfXi)=numberOfGaussXi(1:basis%numberOfXi)
    !Check the number of gauss points is sufficient for the interpolation order and flag a warning if not
    DO xiIdx=1,basis%numberOfXi
      SELECT CASE(basis%interpolationXi(xiIdx))
      CASE(BASIS_LINEAR_LAGRANGE_INTERPOLATION)
        IF(basis%quadrature%numberOfGaussXi(xiIdx)<2) THEN
          localWarning=TRIM(NumberToVString(basis%quadrature%numberOfGaussXi(xiIdx),"*",err,error))// &
            & " Gauss points are insufficient for linear Lagrange interpolation."
          CALL FlagWarning(localWarning,err,error,*999)
        ENDIF
      CASE(BASIS_QUADRATIC_LAGRANGE_INTERPOLATION)
        IF(basis%quadrature%numberOfGaussXi(xiIdx)<2) THEN
          localWarning=TRIM(NumberToVString(basis%quadrature%numberOfGaussXi(xiIdx),"*",err,error))//&
            & " Gauss points are insufficient for quadratic Lagrange interpolation."
          CALL FlagWarning(localWarning,err,error,*999)
        ENDIF
      CASE(BASIS_CUBIC_LAGRANGE_INTERPOLATION)
        IF(basis%quadrature%numberOfGaussXi(xiIdx)<3) THEN
          localWarning=TRIM(NumberToVString(basis%quadrature%numberOfGaussXi(xiIdx),"*",err,error))//&
            & " Gauss points are insufficient for cubic Lagrange interpolation."
          CALL FlagWarning(localWarning,err,error,*999)
        ENDIF
      CASE(BASIS_QUADRATIC1_HERMITE_INTERPOLATION,BASIS_QUADRATIC2_HERMITE_INTERPOLATION)
        IF(basis%quadrature%numberOfGaussXi(xiIdx)<2) THEN
          localWarning=TRIM(NumberToVString(basis%quadrature%numberOfGaussXi(xiIdx),"*",err,error))//&
            & " Gauss points are insufficient for quadratic Hermite interpolation."
          CALL FlagWarning(localWarning,err,error,*999)
        ENDIF
      CASE(BASIS_CUBIC_HERMITE_INTERPOLATION)
        IF(basis%quadrature%numberOfGaussXi(xiIdx)<3) THEN
          localWarning=TRIM(NumberToVString(basis%quadrature%numberOfGaussXi(xiIdx),"*",err,error))//&
            & " Gauss points are insufficient for cubic Hermite interpolation."
          CALL FlagWarning(localWarning,err,error,*999)
        ENDIF
      CASE(BASIS_LINEAR_SIMPLEX_INTERPOLATION)
        CALL FlagWarning("For simplex elements please set quadrature order rather than number of gauss points.",err,error,*999)
      CASE(BASIS_QUADRATIC_SIMPLEX_INTERPOLATION)
        CALL FlagWarning("For simplex elements please set quadrature order rather than number of gauss points.",err,error,*999)
      CASE DEFAULT
        localError="Interpolation xi value "//TRIM(NumberToVString(basis%interpolationXi(xiIdx),"*",err,error))// &
          & " is invalid for xi direction "//TRIM(NumberToVString(xiIdx,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !xiIdx
      
    EXITS("Basis_QuadratureNumberOfGaussXiSet")
    RETURN
999 ERRORSEXITS("Basis_QuadratureNumberOfGaussXiSet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_QuadratureNumberOfGaussXiSet
  !
  !================================================================================================================================
  !

  !>Sets/changes the order of a quadrature for a basis quadrature identified by a pointer. \see OpenCMISS::Iron::cmfe_Basis_QuadratureOrderSet
  SUBROUTINE Basis_QuadratureOrderSet(basis,order,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: order !<The quadrature order to be set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Basis_QuadratureOrderSet",err,error,*999)

    CALL Basis_AssertNotFinished(basis,err,error,*999)
    CALL Basis_AssertIsSimplexBasis(basis,err,error,*999)
    IF(ASSOCIATED(basis%quadrature%basis)) CALL FlagError("Quadrature basis is not associated.",err,error,*999)
    IF(order<1.OR.order>5) THEN
      localError="The specified order of "//TRIM(NumberToVString(order,"*",err,error))// &
        & " is invalid. The order must be >= 1 and <= 5."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    basis%quadrature%gaussOrder=order
      
    EXITS("Basis_QuadratureOrderSet")
    RETURN
999 ERRORSEXITS("Basis_QuadratureOrderSet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_QuadratureOrderSet

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the quadrature type on a basis identified by a pointer.
  SUBROUTINE Basis_QuadratureTypeSet(basis,quadratureType,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: quadratureType !<The quadrature type to be set \see BasisRoutines_QuadratureTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Basis_QuadratureTypeSet",err,error,*999)

    CALL Basis_AssertNotFinished(basis,err,error,*999)
    IF(.NOT.ASSOCIATED(basis%quadrature%basis)) THEN
      localError="The basis quadrature basis is not associated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(quadratureType)
    CASE(BASIS_GAUSS_LEGENDRE_QUADRATURE)
      basis%quadrature%TYPE=basis_GAUSS_LEGENDRE_QUADRATURE
    CASE(BASIS_GAUSS_LAGUERRE_QUADRATURE)
      basis%quadrature%TYPE=BASIS_GAUSS_LAGUERRE_QUADRATURE
      CALL FlagError("Gauss Laguerre quadrature is not implemented.",err,error,*999)
    CASE(BASIS_GUASS_HERMITE_QUADRATURE)
      basis%quadrature%TYPE=BASIS_GUASS_HERMITE_QUADRATURE
      CALL FlagError("Gauss Hermite quadrature is not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The specified quadrature type of "//TRIM(NumberToVString(quadratureType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
     
    EXITS("Basis_QuadratureTypeSet")
    RETURN
999 ERRORSEXITS("Basis_QuadratureTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_QuadratureTypeSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the local face Gauss evaluation flag on a basis
  SUBROUTINE Basis_QuadratureLocalFaceGaussEvaluateSet(basis,faceGaussEvaluate,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis
    LOGICAL, INTENT(IN) :: faceGaussEvaluate !<face Gauss evaluation flag
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Basis_QuadratureLocalFaceGaussEvaluateSet",err,error,*999)

    CALL Basis_AssertNotFinished(basis,err,error,*999)
    
    basis%quadrature%evaluateFaceGauss=faceGaussEvaluate
    
    EXITS("Basis_QuadratureLocalFaceGaussEvaluateSet")
    RETURN
999 ERRORSEXITS("Basis_QuadratureLocalFaceGaussEvaluateSet",err,error)
    RETURN 1

  END SUBROUTINE Basis_QuadratureLocalFaceGaussEvaluateSet

  !
  !================================================================================================================================
  !
  
  !>Creates and initialises a simplex basis that has already been allocated Basis_CreateStart
  !>\see BasisRoutines::Basis_CreateStart
  SUBROUTINE Basis_SimplexBasisCreate(basis,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: maximumNumberOfNodes,xiIdx,localNodeIdx,elementParametersIdx
    INTEGER(INTG), ALLOCATABLE :: nodesInFace(:)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Basis_SimplexBasisCreate",err,error,*999)

    CALL Basis_AssertIsSimplexBasis(basis,err,error,*999)
    
    basis%numberOfXiCoordinates=basis%numberOfXi+1 !Simplex bases have an additional area coordinate
    ALLOCATE(basis%interpolationType(basis%numberOfXiCoordinates),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate the basis interpolation type array.",err,error,*999)
    ALLOCATE(basis%interpolationOrder(basis%numberOfXiCoordinates),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate the basis interpolation order array.",err,error,*999)
    ALLOCATE(basis%numberOfNodesXiC(basis%numberOfXiCoordinates),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate the basis number of nodes xic array.",err,error,*999)
    basis%degenerate=.FALSE.
    basis%numberOfCollapsedXi=0
    SELECT CASE(basis%numberOfXi)
    CASE(1)
      basis%numberOfPartialDerivatives=3
      SELECT CASE(basis%interpolationXi(1))
      CASE(BASIS_LINEAR_SIMPLEX_INTERPOLATION)
        basis%interpolationType(1)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(1)=BASIS_LINEAR_INTERPOLATION_ORDER
        basis%interpolationType(2)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(2)=BASIS_LINEAR_INTERPOLATION_ORDER
        basis%numberOfNodesXiC(1)=2            
        basis%numberOfNodesXiC(2)=2            
        maximumNumberOfNodes=2
        basis%numberOfNodes=2
      CASE(BASIS_QUADRATIC_SIMPLEX_INTERPOLATION)
        basis%interpolationType(1)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(1)=BASIS_QUADRATIC_INTERPOLATION_ORDER
        basis%interpolationType(2)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(2)=BASIS_QUADRATIC_INTERPOLATION_ORDER
        basis%numberOfNodesXiC(1)=3
        basis%numberOfNodesXiC(2)=3
        maximumNumberOfNodes=3
        basis%numberOfNodes=3
      CASE(BASIS_CUBIC_SIMPLEX_INTERPOLATION)
        basis%interpolationType(1)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(1)=BASIS_CUBIC_INTERPOLATION_ORDER
        basis%interpolationType(2)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(2)=BASIS_CUBIC_INTERPOLATION_ORDER
        basis%numberOfNodesXiC(1)=4
        basis%numberOfNodesXiC(2)=4
        maximumNumberOfNodes=4
        basis%numberOfNodes=4
      CASE DEFAULT
        localError="The basis interpolation xi of "//TRIM(NumberToVString(basis%interpolationXi(1),"*",err,error))// &
          & " at xi index 1 is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(2)
      basis%numberOfPartialDerivatives=6
      SELECT CASE(basis%interpolationXi(2))
      CASE(BASIS_LINEAR_SIMPLEX_INTERPOLATION)
        basis%interpolationType(1)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(1)=BASIS_LINEAR_INTERPOLATION_ORDER
        basis%interpolationType(2)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(2)=BASIS_LINEAR_INTERPOLATION_ORDER
        basis%interpolationType(3)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(3)=BASIS_LINEAR_INTERPOLATION_ORDER
        basis%numberOfNodesXiC(1)=2            
        BASIS%numberOfNodesXiC(2)=2            
        basis%numberOfNodesXiC(3)=2            
        maximumNumberOfNodes=2
        BASIS%numberOfNodes=3
      CASE(BASIS_QUADRATIC_SIMPLEX_INTERPOLATION)
        basis%interpolationType(1)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(1)=BASIS_QUADRATIC_INTERPOLATION_ORDER
        basis%interpolationType(2)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(2)=BASIS_QUADRATIC_INTERPOLATION_ORDER
        basis%interpolationType(3)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(3)=BASIS_QUADRATIC_INTERPOLATION_ORDER
        basis%numberOfNodesXiC(1)=3
        basis%numberOfNodesXiC(2)=3
        basis%numberOfNodesXiC(3)=3
        maximumNumberOfNodes=3
        basis%numberOfNodes=6
      CASE(BASIS_CUBIC_SIMPLEX_INTERPOLATION)
        basis%interpolationType(1)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(1)=BASIS_CUBIC_INTERPOLATION_ORDER
        basis%interpolationType(2)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(2)=BASIS_CUBIC_INTERPOLATION_ORDER
        basis%interpolationType(3)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(3)=BASIS_CUBIC_INTERPOLATION_ORDER
        basis%numberOfNodesXiC(1)=4
        basis%numberOfNodesXiC(2)=4
        basis%numberOfNodesXiC(3)=4
        maximumNumberOfNodes=4
        basis%numberOfNodes=10
      CASE DEFAULT 
        localError="The basis interpolation xi of "//TRIM(NumberToVString(basis%interpolationXi(1),"*",err,error))// &
          & " at xi index 2 is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(3)
      basis%numberOfPartialDerivatives=11
      SELECT CASE(basis%interpolationXi(3))
      CASE(BASIS_LINEAR_SIMPLEX_INTERPOLATION)
        basis%interpolationType(1)=basis_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(1)=BASIS_LINEAR_INTERPOLATION_ORDER
        basis%interpolationType(2)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(2)=BASIS_LINEAR_INTERPOLATION_ORDER
        basis%interpolationType(3)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(3)=BASIS_LINEAR_INTERPOLATION_ORDER
        basis%interpolationType(4)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(4)=BASIS_LINEAR_INTERPOLATION_ORDER
        basis%numberOfNodesXiC(1)=2            
        BASIS%numberOfNodesXiC(2)=2            
        basis%numberOfNodesXiC(3)=2            
        BASIS%numberOfNodesXiC(4)=2            
        maximumNumberOfNodes=2
        basis%numberOfNodes=4
      CASE(BASIS_QUADRATIC_SIMPLEX_INTERPOLATION)
        basis%interpolationType(1)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(1)=BASIS_QUADRATIC_INTERPOLATION_ORDER
        basis%interpolationType(2)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(2)=BASIS_QUADRATIC_INTERPOLATION_ORDER
        basis%interpolationType(3)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(3)=BASIS_QUADRATIC_INTERPOLATION_ORDER
        basis%interpolationType(4)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(4)=BASIS_QUADRATIC_INTERPOLATION_ORDER
        basis%numberOfNodesXiC(1)=3
        basis%numberOfNodesXiC(2)=3
        basis%numberOfNodesXiC(3)=3
        basis%numberOfNodesXiC(4)=3
        maximumNumberOfNodes=3
        basis%numberOfNodes=10
      CASE(BASIS_CUBIC_SIMPLEX_INTERPOLATION)
        basis%interpolationType(1)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(1)=BASIS_CUBIC_INTERPOLATION_ORDER
        basis%interpolationType(2)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(2)=BASIS_CUBIC_INTERPOLATION_ORDER
        basis%interpolationType(3)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(3)=BASIS_CUBIC_INTERPOLATION_ORDER
        basis%interpolationType(4)=BASIS_SIMPLEX_INTERPOLATION
        basis%interpolationOrder(4)=BASIS_CUBIC_INTERPOLATION_ORDER
        basis%numberOfNodesXiC(1)=4
        basis%numberOfNodesXiC(2)=4
        basis%numberOfNodesXiC(3)=4
        basis%numberOfNodesXiC(4)=4
        maximumNumberOfNodes=4
        basis%numberOfNodes=20
      CASE DEFAULT 
        localError="The basis interpolation xi of "//TRIM(NumberToVString(basis%interpolationXi(1),"*",err,error))// &
          & " at xi index 2 is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The number of xi directions of "//TRIM(NumberToVString(basis%numberOfXi,"*",err,error))// &
        & " for basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))// &
        & " is invalid. The number of xi directions should be >= 1 and <= 3."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    ALLOCATE(basis%nodeAtCollapse(basis%numberOfNodes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate node at collapse.",err,error,*999)
    basis%nodeAtCollapse=.FALSE.
    
    ALLOCATE(basis%nodePositionIndex(basis%numberOfNodes,basis%numberOfXiCoordinates),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate node position index.",err,error,*999) 
    SELECT CASE(basis%numberOfXiCoordinates)
    CASE(2)
      ALLOCATE(basis%nodePositionIndexInv(maximumNumberOfNodes,maximumNumberOfNodes,1,1),STAT=err)
    CASE(3)
      ALLOCATE(basis%nodePositionIndexInv(maximumNumberOfNodes,maximumNumberOfNodes,maximumNumberOfNodes,1),STAT=err)
    CASE(4)
      ALLOCATE(basis%nodePositionIndexInv(maximumNumberOfNodes,maximumNumberOfNodes,maximumNumberOfNodes, &
        & maximumNumberOfNodes),STAT=err)
    CASE DEFAULT
      localError="The number of xi directions of "//TRIM(NumberToVString(basis%numberOfXiCoordinates,"*",err,error))// &
        & " for basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))// &
        & " is invalid. The number of xi directions should be >= 2 and <= 4 for simplex bases."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    IF(err/=0) CALL FlagError("Could not allocate inverse node position index.",err,error,*999)
    basis%nodePositionIndexInv=0
        
    !Determine the node position index and it's inverse
    SELECT CASE(basis%numberOfXi)
    CASE(1)
      SELECT CASE(basis%interpolationOrder(1))
      CASE(basis_LINEAR_INTERPOLATION_ORDER)
        !Node 1
        basis%nodePositionIndex(1,1)=2
        BASIS%nodePositionIndex(1,2)=1
        basis%nodePositionIndexInv(2,1,1,1)=1
        !Node 2
        basis%nodePositionIndex(2,1)=1
        basis%nodePositionIndex(2,2)=2
        basis%nodePositionIndexInv(1,2,1,1)=2
      CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
        !Node 1
        basis%nodePositionIndex(1,1)=3
        basis%nodePositionIndex(1,2)=1
        basis%nodePositionIndexInv(3,1,1,1)=1
        !Node 2
        basis%nodePositionIndex(2,1)=2
        basis%nodePositionIndex(2,2)=2
        basis%nodePositionIndexInv(2,2,1,1)=2
        !Node 3
        basis%nodePositionIndex(3,1)=1
        basis%nodePositionIndex(3,2)=3
        basis%nodePositionIndexInv(1,3,1,1)=3
      CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
        !Node 1
        basis%nodePositionIndex(1,1)=4
        basis%nodePositionIndex(1,2)=1
        basis%nodePositionIndexInv(4,1,1,1)=1
        !Node 2
        basis%nodePositionIndex(2,1)=3
        basis%nodePositionIndex(2,2)=2
        basis%nodePositionIndexInv(3,2,1,1)=2
        !Node 3
        basis%nodePositionIndex(3,1)=2
        basis%nodePositionIndex(3,2)=3
        basis%nodePositionIndexInv(2,3,1,1)=3
        !Node 4
        basis%nodePositionIndex(4,1)=1
        basis%nodePositionIndex(4,2)=4
        basis%nodePositionIndexInv(1,4,1,1)=4
      CASE DEFAULT
        localError="The basis interpolation order of "//TRIM(NumberToVString(basis%interpolationOrder(1),"*",err,error))// &
          & " at xi coordinate index 1 is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(2)
      SELECT CASE(basis%interpolationOrder(1))
      CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
        !Node 1
        basis%nodePositionIndex(1,1)=2
        basis%nodePositionIndex(1,2)=1
        basis%nodePositionIndex(1,3)=1
        basis%nodePositionIndexInv(2,1,1,1)=1
        !Node 2
        basis%nodePositionIndex(2,1)=1
        basis%nodePositionIndex(2,2)=2
        basis%nodePositionIndex(2,3)=1
        basis%nodePositionIndexInv(1,2,1,1)=2
        !Node 3
        basis%nodePositionIndex(3,1)=1
        basis%nodePositionIndex(3,2)=1
        basis%nodePositionIndex(3,3)=2
        basis%nodePositionIndexInv(1,1,2,1)=3
      CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
        !Node 1
        basis%nodePositionIndex(1,1)=3
        basis%nodePositionIndex(1,2)=1
        basis%nodePositionIndex(1,3)=1
        basis%nodePositionIndexInv(3,1,1,1)=1
        !Node 2
        basis%nodePositionIndex(2,1)=1
        basis%nodePositionIndex(2,2)=3
        basis%nodePositionIndex(2,3)=1
        basis%nodePositionIndexInv(1,3,1,1)=2
        !Node 3
        basis%nodePositionIndex(3,1)=1
        basis%nodePositionIndex(3,2)=1
        basis%nodePositionIndex(3,3)=3
        basis%nodePositionIndexInv(1,1,3,1)=3
        !Node 4
        basis%nodePositionIndex(4,1)=2
        basis%nodePositionIndex(4,2)=2
        basis%nodePositionIndex(4,3)=1
        basis%nodePositionIndexInv(2,2,1,1)=4
        !Node 5
        basis%nodePositionIndex(5,1)=1
        basis%nodePositionIndex(5,2)=2
        basis%nodePositionIndex(5,3)=2
        basis%nodePositionIndexInv(1,2,2,1)=5
        !Node 6
        basis%nodePositionIndex(6,1)=2
        basis%nodePositionIndex(6,2)=1
        basis%nodePositionIndex(6,3)=2
        basis%nodePositionIndexInv(2,1,2,1)=6
      CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
        !Node 1
        basis%nodePositionIndex(1,1)=4
        basis%nodePositionIndex(1,2)=1
        basis%nodePositionIndex(1,3)=1
        basis%nodePositionIndexInv(4,1,1,1)=1
        !Node 2
        basis%nodePositionIndex(2,1)=1
        basis%nodePositionIndex(2,2)=4
        basis%nodePositionIndex(2,3)=1
        basis%nodePositionIndexInv(1,4,1,1)=2
        !Node 3
        basis%nodePositionIndex(3,1)=1
        basis%nodePositionIndex(3,2)=1
        basis%nodePositionIndex(3,3)=4
        basis%nodePositionIndexInv(1,1,4,1)=3
        !Node 4
        basis%nodePositionIndex(4,1)=3
        basis%nodePositionIndex(4,2)=2
        basis%nodePositionIndex(4,3)=1
        basis%nodePositionIndexInv(3,2,1,1)=4
        !Node 5
        basis%nodePositionIndex(5,1)=2
        basis%nodePositionIndex(5,2)=3
        basis%nodePositionIndex(5,3)=1
        basis%nodePositionIndexInv(2,3,1,1)=5
        !Node 6
        basis%nodePositionIndex(6,1)=1
        basis%nodePositionIndex(6,2)=3
        basis%nodePositionIndex(6,3)=2
        basis%nodePositionIndexInv(1,3,2,1)=6
        !Node 7
        basis%nodePositionIndex(7,1)=1
        basis%nodePositionIndex(7,2)=2
        basis%nodePositionIndex(7,3)=3
        basis%nodePositionIndexInv(1,2,3,1)=7
        !Node 8
        basis%nodePositionIndex(8,1)=2
        basis%nodePositionIndex(8,2)=1
        basis%nodePositionIndex(8,3)=3
        basis%nodePositionIndexInv(2,1,3,1)=8
        !Node 9
        basis%nodePositionIndex(9,1)=3
        basis%nodePositionIndex(9,2)=1
        basis%nodePositionIndex(9,3)=2
        basis%nodePositionIndexInv(3,1,2,1)=9
        !Node 10
        basis%nodePositionIndex(10,1)=2
        basis%nodePositionIndex(10,2)=2
        basis%nodePositionIndex(10,3)=2
        basis%nodePositionIndexInv(2,2,2,1)=10
      CASE DEFAULT
        localError="The basis interpolation order of "//TRIM(NumberToVString(basis%interpolationOrder(1),"*",err,error))// &
          & " at xi coordinate index 1 is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(3)
      SELECT CASE(basis%interpolationOrder(1))
      CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
        !Node 1
        basis%nodePositionIndex(1,1)=2
        basis%nodePositionIndex(1,2)=1
        basis%nodePositionIndex(1,3)=1
        basis%nodePositionIndex(1,4)=1
        basis%nodePositionIndexInv(2,1,1,1)=1
        !Node 2
        basis%nodePositionIndex(2,1)=1
        basis%nodePositionIndex(2,2)=2
        basis%nodePositionIndex(2,3)=1
        basis%nodePositionIndex(2,4)=1
        basis%nodePositionIndexInv(1,2,1,1)=2
        !Node 3
        basis%nodePositionIndex(3,1)=1
        basis%nodePositionIndex(3,2)=1
        basis%nodePositionIndex(3,3)=2
        basis%nodePositionIndex(3,4)=1
        basis%nodePositionIndexInv(1,1,2,1)=3
        !Node 4
        basis%nodePositionIndex(4,1)=1
        basis%nodePositionIndex(4,2)=1
        basis%nodePositionIndex(4,3)=1
        basis%nodePositionIndex(4,4)=2
        basis%nodePositionIndexInv(1,1,1,2)=4
        
        ALLOCATE(nodesInFace(12),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate nodes in face.",err,error,*999) 
        nodesInFace(:)=[2,3,4,1,3,4,1,2,4,1,2,3] !12 Nodes
        
      CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
        !Node 1
        basis%nodePositionIndex(1,1)=3
        basis%nodePositionIndex(1,2)=1
        basis%nodePositionIndex(1,3)=1
        basis%nodePositionIndex(1,4)=1
        basis%nodePositionIndexInv(3,1,1,1)=1
        !Node 2
        basis%nodePositionIndex(2,1)=1
        basis%nodePositionIndex(2,2)=3
        basis%nodePositionIndex(2,3)=1
        basis%nodePositionIndex(2,4)=1
        basis%nodePositionIndexInv(1,3,1,1)=2
        !Node 3
        basis%nodePositionIndex(3,1)=1
        basis%nodePositionIndex(3,2)=1
        basis%nodePositionIndex(3,3)=3
        basis%nodePositionIndex(3,4)=1
        basis%nodePositionIndexInv(1,1,3,1)=3
        !Node 4
        basis%nodePositionIndex(4,1)=1
        basis%nodePositionIndex(4,2)=1
        basis%nodePositionIndex(4,3)=1
        basis%nodePositionIndex(4,4)=3
        basis%nodePositionIndexInv(1,1,1,3)=4
        !Node 5
        basis%nodePositionIndex(5,1)=2
        basis%nodePositionIndex(5,2)=2
        basis%nodePositionIndex(5,3)=1
        basis%nodePositionIndex(5,4)=1
        basis%nodePositionIndexInv(2,2,1,1)=5
        !Node 6
        basis%nodePositionIndex(6,1)=2
        basis%nodePositionIndex(6,2)=1
        basis%nodePositionIndex(6,3)=2
        basis%nodePositionIndex(6,4)=1
        basis%nodePositionIndexInv(2,1,2,1)=6
        !Node 7
        basis%nodePositionIndex(7,1)=2
        basis%nodePositionIndex(7,2)=1
        basis%nodePositionIndex(7,3)=1
        basis%nodePositionIndex(7,4)=2
        basis%nodePositionIndexInv(2,1,1,2)=7
        !Node 8
        basis%nodePositionIndex(8,1)=1
        basis%nodePositionIndex(8,2)=2
        basis%nodePositionIndex(8,3)=2
        basis%nodePositionIndex(8,4)=1
        basis%nodePositionIndexInv(1,2,2,1)=8
        !Node 9
        basis%nodePositionIndex(9,1)=1
        basis%nodePositionIndex(9,2)=1
        basis%nodePositionIndex(9,3)=2
        basis%nodePositionIndex(9,4)=2
        basis%nodePositionIndexInv(1,1,2,2)=9
        !Node 10
        basis%nodePositionIndex(10,1)=1
        basis%nodePositionIndex(10,2)=2
        basis%nodePositionIndex(10,3)=1
        basis%nodePositionIndex(10,4)=2
        basis%nodePositionIndexInv(1,2,1,2)=10
        
        ALLOCATE(nodesInFace(24),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate nodes in face.",err,error,*999) 
        nodesInFace(:)=[2,3,4,8,9,10,1,3,4,6,9,7,1,2,4,5,10,7,1,2,3,5,8,6] !24 Nodes

      CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
        !Node 1
        basis%nodePositionIndex(1,1)=4
        basis%nodePositionIndex(1,2)=1
        basis%nodePositionIndex(1,3)=1
        basis%nodePositionIndex(1,4)=1
        basis%nodePositionIndexInv(4,1,1,1)=1
        !Node 2
        basis%nodePositionIndex(2,1)=1
        basis%nodePositionIndex(2,2)=4
        basis%nodePositionIndex(2,3)=1
        basis%nodePositionIndex(2,4)=1
        basis%nodePositionIndexInv(1,4,1,1)=2
        !Node 3
        basis%nodePositionIndex(3,1)=1
        basis%nodePositionIndex(3,2)=1
        basis%nodePositionIndex(3,3)=4
        basis%nodePositionIndex(3,4)=1
        basis%nodePositionIndexInv(1,1,4,1)=3
        !Node 4
        basis%nodePositionIndex(4,1)=1
        basis%nodePositionIndex(4,2)=1
        basis%nodePositionIndex(4,3)=1
        basis%nodePositionIndex(4,4)=4
        basis%nodePositionIndexInv(1,1,1,4)=4
        !Node 5
        basis%nodePositionIndex(5,1)=3
        basis%nodePositionIndex(5,2)=2
        basis%nodePositionIndex(5,3)=1
        basis%nodePositionIndex(5,4)=1
        basis%nodePositionIndexInv(3,2,1,1)=5
        !Node 6
        basis%nodePositionIndex(6,1)=2
        basis%nodePositionIndex(6,2)=3
        basis%nodePositionIndex(6,3)=1
        basis%nodePositionIndex(6,4)=1
        basis%nodePositionIndexInv(2,3,1,1)=6
        !Node 7
        basis%nodePositionIndex(7,1)=3
        basis%nodePositionIndex(7,2)=1
        basis%nodePositionIndex(7,3)=2
        basis%nodePositionIndex(7,4)=1
        basis%nodePositionIndexInv(3,1,2,1)=7
        !Node 8
        basis%nodePositionIndex(8,1)=2
        basis%nodePositionIndex(8,2)=1
        basis%nodePositionIndex(8,3)=3
        basis%nodePositionIndex(8,4)=1
        basis%nodePositionIndexInv(2,1,3,1)=8
        !Node 9
        basis%nodePositionIndex(9,1)=3
        basis%nodePositionIndex(9,2)=1
        basis%nodePositionIndex(9,3)=1
        basis%nodePositionIndex(9,4)=2
        basis%nodePositionIndexInv(3,1,1,2)=9
        !Node 10
        basis%nodePositionIndex(10,1)=2
        basis%nodePositionIndex(10,2)=1
        basis%nodePositionIndex(10,3)=1
        basis%nodePositionIndex(10,4)=3
        basis%nodePositionIndexInv(2,1,1,3)=10
        !Node 11
        basis%nodePositionIndex(11,1)=1
        basis%nodePositionIndex(11,2)=3
        basis%nodePositionIndex(11,3)=2
        basis%nodePositionIndex(11,4)=1
        basis%nodePositionIndexInv(1,3,2,1)=11
        !Node 12
        basis%nodePositionIndex(12,1)=1
        basis%nodePositionIndex(12,2)=2
        basis%nodePositionIndex(12,3)=3
        basis%nodePositionIndex(12,4)=1
        basis%nodePositionIndexInv(1,2,3,1)=12
        !Node 13
        basis%nodePositionIndex(13,1)=1
        basis%nodePositionIndex(13,2)=1
        basis%nodePositionIndex(13,3)=3
        basis%nodePositionIndex(13,4)=2
        basis%nodePositionIndexInv(1,1,3,2)=13
        !Node 14
        basis%nodePositionIndex(14,1)=1
        basis%nodePositionIndex(14,2)=1
        basis%nodePositionIndex(14,3)=2
        basis%nodePositionIndex(14,4)=3
        basis%nodePositionIndexInv(1,1,2,3)=14
        !Node 15
        basis%nodePositionIndex(15,1)=1
        basis%nodePositionIndex(15,2)=3
        basis%nodePositionIndex(15,3)=1
        basis%nodePositionIndex(15,4)=2
        basis%nodePositionIndexInv(1,3,1,2)=15
        !Node 16
        basis%nodePositionIndex(16,1)=1
        basis%nodePositionIndex(16,2)=2
        basis%nodePositionIndex(16,3)=1
        basis%nodePositionIndex(16,4)=3
        basis%nodePositionIndexInv(1,2,1,3)=16
        !Node 17
        basis%nodePositionIndex(17,1)=2
        basis%nodePositionIndex(17,2)=2
        basis%nodePositionIndex(17,3)=2
        basis%nodePositionIndex(17,4)=1
        basis%nodePositionIndexInv(2,2,2,1)=17
        !Node 18
        basis%nodePositionIndex(18,1)=2
        basis%nodePositionIndex(18,2)=2
        basis%nodePositionIndex(18,3)=1
        basis%nodePositionIndex(18,4)=2
        basis%nodePositionIndexInv(2,2,1,2)=18
        !Node 19
        basis%nodePositionIndex(19,1)=2
        basis%nodePositionIndex(19,2)=1
        basis%nodePositionIndex(19,3)=2
        basis%nodePositionIndex(19,4)=2
        basis%nodePositionIndexInv(2,1,2,2)=19
        !Node 20
        basis%nodePositionIndex(20,1)=1
        basis%nodePositionIndex(20,2)=2
        basis%nodePositionIndex(20,3)=2
        basis%nodePositionIndex(20,4)=2
        basis%nodePositionIndexInv(1,2,2,2)=20
        
        ALLOCATE(nodesInFace(40),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate nodes in face.",err,error,*999) 
        nodesInFace(:)=[2,3,4,11,12,13,14,16,15,20,1,3,4,7,8,13,14,10,9,&
          &19,1,2,4,5,6,15,16,10,9,18,1,2,3,5,6,14,12,8,7,17] !40 nodes
        
      CASE DEFAULT
        localError="The basis interpolation order of "//TRIM(NumberToVString(basis%interpolationOrder(1),"*",err,error))// &
          & " at xi coordinate index 1 is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The number of xi directions of "//TRIM(NumberToVString(basis%numberOfXi,"*",err,error))// &
        & " for basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))// &
        & " is invalid. The number of xi directions should be >= 1 and <= 3."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    !Calculate the maximum number of derivatives (1 for simplex bases) and the number of element parameters
    basis%maximumNumberOfDerivatives=1
    basis%numberOfElementParameters=basis%numberOfNodes
    !Now set up the number of derivatives and derivative order index
    ALLOCATE(basis%numberOfDerivatives(basis%numberOfNodes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate number of derivatives.",err,error,*999)
    ALLOCATE(basis%derivativeOrderIndex(basis%maximumNumberOfDerivatives,basis%numberOfNodes,basis%numberOfXi),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate derivative order index.",err,error,*999)
    ALLOCATE(basis%derivativeOrderIndexInv(FIRST_PART_DERIV,FIRST_PART_DERIV,FIRST_PART_DERIV,basis%numberOfNodes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate inverse derivative order index.",err,error,*999)
    ALLOCATE(basis%partialDerivativeIndex(basis%maximumNumberOfDerivatives,basis%numberOfNodes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate partial derivative index.",err,error,*999)
    ALLOCATE(basis%elementParameterIndex(basis%maximumNumberOfDerivatives,basis%numberOfNodes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate element parameter index.",err,error,*999)
    ALLOCATE(basis%elementParameterIndexInv(2,basis%numberOfElementParameters),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate inverse element parameter index.",err,error,*999)
    !Set the derivative order index and its inverse, the element parameter index and the partial derivative index.
    elementParametersIdx=0
    basis%derivativeOrderIndexInv=0
    DO localNodeIdx=1,basis%numberOfNodes
      basis%numberOfDerivatives(localNodeIdx)=1
      DO xiIdx=1,basis%numberOfXi
        basis%derivativeOrderIndex(1,localNodeIdx,xiIdx)=1
      ENDDO !xiIdx
      elementParametersIdx=elementParametersIdx+1
      basis%elementParameterIndex(1,localNodeIdx)=elementParametersIdx
      basis%elementParameterIndexInv(1,elementParametersIdx)=localNodeIdx
      basis%elementParameterIndexInv(2,elementParametersIdx)=1
      basis%partialDerivativeIndex(1,localNodeIdx)=NO_PART_DERIV
      basis%derivativeOrderIndexInv(basis%derivativeOrderIndex(1,localNodeIdx,1),1,1,localNodeIdx)=1
    ENDDO !localNodeIdx
      
    !Set up the line and face information
    SELECT CASE(basis%numberOfXi)
    CASE(1)
      basis%numberOfLocalLines=1
      ALLOCATE(basis%numberOfNodesInLocalLine(basis%numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate number of nodes in local line.",err,error,*999)
      ALLOCATE(basis%localLineXiDirection(basis%numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate local line xi direction.",err,error,*999)
      basis%localLineXiDirection(1)=1
      ALLOCATE(basis%nodeNumbersInLocalLine(basis%numberOfNodesXiC(1),basis%numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate node numbers in local line.",err,error,*999)
      ALLOCATE(basis%derivativeNumbersInLocalLine(basis%numberOfNodesXiC(1),basis%numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate derivative numbers in local line.",err,error,*999)
      basis%derivativeNumbersInLocalLine=NO_PART_DERIV
      ALLOCATE(basis%elementParametersInLocalLine(basis%numberOfNodesXiC(1)**2,basis%numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate element parameters in local line.",err,error,*999)
      basis%elementParametersInLocalLine=1
      !Set the line values
      SELECT CASE(basis%interpolationOrder(1))
      CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
        !Line 1
        basis%numberOfNodesInLocalLine(1)=2
        basis%nodeNumbersInLocalLine(1,1)=1
        basis%nodeNumbersInLocalLine(2,1)=2
      CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
        !Line 1
        basis%numberOfNodesInLocalLine(1)=3
        basis%nodeNumbersInLocalLine(1,1)=1
        basis%nodeNumbersInLocalLine(2,1)=2
        basis%nodeNumbersInLocalLine(3,1)=3
      CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
        !Line 1
        basis%numberOfNodesInLocalLine(1)=4
        basis%nodeNumbersInLocalLine(1,1)=1
        basis%nodeNumbersInLocalLine(2,1)=2
        basis%nodeNumbersInLocalLine(3,1)=3
        basis%nodeNumbersInLocalLine(4,1)=4
      CASE DEFAULT 
        localError="The basis interpolation order of "//TRIM(NumberToVString(basis%interpolationOrder(1),"*",err,error))// &
          & " at xi coordinate index 1 is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(2)
      !Allocate and calculate the lines
      !Simplex hence three local lines
      basis%numberOfLocalLines=3
      ALLOCATE(basis%localLineXiDirection(basis%numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate local line xi direction",err,error,*999)
      ALLOCATE(basis%localLineXiNormals(1,basis%numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate local line xi normals.",err,error,*999)
      ALLOCATE(basis%xiNormalsLocalLine(-3:3,1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate xi normals local line.",err,error,*999)
      basis%xiNormalsLocalLine=0
      ALLOCATE(basis%numberOfNodesInLocalLine(basis%numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate number of nodes in local line.",err,error,*999)
      ALLOCATE(basis%nodeNumbersInLocalLine(MAXVAL(basis%numberOfNodesXiC),basis%numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate node numbers in local line.",err,error,*999)
      ALLOCATE(basis%derivativeNumbersInLocalLine(MAXVAL(basis%numberOfNodesXiC),basis%numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate derivative numbers in local line.",err,error,*999)
      basis%derivativeNumbersInLocalLine=NO_PART_DERIV
      ALLOCATE(basis%elementParametersInLocalLine(MAXVAL(basis%numberOfNodesXiC)**2,basis%numberOfLocalLines) &
        & ,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate element parameters in local line.",err,error,*999)
      basis%elementParametersInLocalLine=1
      !Set the line values
      SELECT CASE(basis%interpolationOrder(1))
      CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
        !Line 1
        basis%numberOfNodesInLocalLine(1)=2
        basis%nodeNumbersInLocalLine(1,1)=2
        basis%nodeNumbersInLocalLine(2,1)=3
        !Line 2
        basis%numberOfNodesInLocalLine(2)=2
        basis%nodeNumbersInLocalLine(1,2)=3
        basis%nodeNumbersInLocalLine(2,2)=1
        !Line 3
        basis%numberOfNodesInLocalLine(3)=2
        basis%nodeNumbersInLocalLine(1,3)=1
        basis%nodeNumbersInLocalLine(2,3)=2
      CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
        !Line 1
        basis%numberOfNodesInLocalLine(1)=3
        basis%nodeNumbersInLocalLine(1,1)=2
        basis%nodeNumbersInLocalLine(2,1)=5
        basis%nodeNumbersInLocalLine(3,1)=3
        !Line 2
        basis%numberOfNodesInLocalLine(2)=3
        basis%nodeNumbersInLocalLine(1,2)=3
        basis%nodeNumbersInLocalLine(2,2)=6
        basis%nodeNumbersInLocalLine(3,2)=2
        !Line 3
        basis%numberOfNodesInLocalLine(3)=3
        basis%nodeNumbersInLocalLine(1,3)=1
        basis%nodeNumbersInLocalLine(2,3)=4
        basis%nodeNumbersInLocalLine(3,3)=2
      CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
        !Line 1
        basis%numberOfNodesInLocalLine(1)=4
        basis%nodeNumbersInLocalLine(1,1)=2
        basis%nodeNumbersInLocalLine(2,1)=6
        basis%nodeNumbersInLocalLine(3,1)=7
        basis%nodeNumbersInLocalLine(4,1)=3
        !Line 2
        basis%numberOfNodesInLocalLine(2)=4
        basis%nodeNumbersInLocalLine(1,2)=3
        basis%nodeNumbersInLocalLine(2,2)=8
        basis%nodeNumbersInLocalLine(3,2)=9
        basis%nodeNumbersInLocalLine(4,2)=1
        !Line 3
        basis%numberOfNodesInLocalLine(3)=4
        basis%nodeNumbersInLocalLine(1,3)=1
        basis%nodeNumbersInLocalLine(2,3)=4
        basis%nodeNumbersInLocalLine(3,3)=5
        basis%nodeNumbersInLocalLine(4,3)=2
      CASE DEFAULT 
        localError="The basis interpolation order of "//TRIM(NumberToVString(basis%interpolationOrder(1),"*",err,error))// &
          & " at xi coordinate index 1 is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set line directions and normals
      !Line 1
      basis%localLineXiDirection(1)=1
      basis%localLineXiNormals(1,1)=3
      basis%xiNormalsLocalLine(3,1)=1
      !Line 2
      basis%localLineXiDirection(2)=2
      basis%localLineXiNormals(1,2)=2
      basis%xiNormalsLocalLine(2,1)=2
      !Line 3
      basis%localLineXiDirection(3)=3
      basis%localLineXiNormals(1,3)=1
      basis%xiNormalsLocalLine(1,1)=3
    CASE(3)
      basis%numberOfLocalLines=6
      basis%numberOfLocalFaces=4
      
      ALLOCATE(basis%localLineXiDirection(basis%numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate local line xi direction.",err,error,*999)
      ALLOCATE(basis%localLineXiNormals(2,basis%numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate local line xi normals.",err,error,*999)
      ALLOCATE(basis%xiNormalsLocalLine(-4:4,-4:4),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate xi normals local line.",err,error,*999)
      basis%xiNormalsLocalLine=0
      
      ALLOCATE(basis%numberOfNodesInLocalLine(basis%numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate number of nodes in local line.",err,error,*999)
      ALLOCATE(basis%numberOfNodesInLocalFace(basis%numberOfLocalFaces),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate number of nodes in local face.",err,error,*999)
      
      ALLOCATE(basis%nodeNumbersInLocalLine(MAXVAL(basis%numberOfNodesXiC),basis%numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate node numbers in local line.",err,error,*999)
      ALLOCATE(basis%derivativeNumbersInLocalLine(MAXVAL(basis%numberOfNodesXiC),basis%numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate derivative numbers in local line.",err,error,*999)
      basis%derivativeNumbersInLocalLine=NO_PART_DERIV
      ALLOCATE(basis%elementParametersInLocalLine(MAXVAL(basis%numberOfNodesXiC)**2,basis%numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate element parameters in local line.",err,error,*999)
      basis%elementParametersInLocalLine=1
      
      ALLOCATE(basis%localFaceXiDirections(3,basis%numberOfLocalFaces),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate local face xi directions.",err,error,*999)
      ALLOCATE(basis%localFaceXiNormal(basis%numberOfLocalFaces),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate local face xi normal.",err,error,*999)
      ALLOCATE(basis%xiNormalLocalFace(-4:4),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate xi normal local face.",err,error,*999)
      basis%xiNormalLocalFace=0
      
      SELECT CASE(basis%interpolationOrder(1))
      CASE(basis_LINEAR_INTERPOLATION_ORDER)
        ALLOCATE(basis%nodeNumbersInLocalFace(3,basis%numberOfLocalFaces),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate node numbers in local face.",err,error,*999)
        !\todo Number of local face node derivatives currenlty set to 1 (calculation of basis%derivativeNumbersInLocalFace for simplex elements has not been implemented yet)
        ALLOCATE(basis%derivativeNumbersInLocalFace(0:1,3,basis%numberOfLocalFaces),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate derivative numbers in local face.",err,error,*999)
      CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
        ALLOCATE(basis%nodeNumbersInLocalFace(6,basis%numberOfLocalFaces),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate node numbers in local face.",err,error,*999)
        !\todo Number of local face node derivatives currenlty set to 1 (calculation of basis%derivativeNumbersInLocalFace for simplex elements has not been implemented yet)
        ALLOCATE(basis%derivativeNumbersInLocalFace(0:1,6,basis%numberOfLocalFaces),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate derivative numbers in local face.",err,error,*999)
      CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
        ALLOCATE(basis%nodeNumbersInLocalFace(10,basis%numberOfLocalFaces),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate node numbers in local face.",err,error,*999)
        !\todo Number of local face node derivatives currenlty set to 1 (calculation of basis%derivativeNumbersInLocalFace for simplex elements has not been implemented yet)
        ALLOCATE(basis%derivativeNumbersInLocalFace(0:1,10,basis%numberOfLocalFaces),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate derivative numbers in local face.",err,error,*999)
      CASE DEFAULT
        localError="The basis interpolation order of "//TRIM(NumberToVString(basis%interpolationOrder(1),"*",err,error))// &
          & " at xi coordinate index 1 is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      basis%derivativeNumbersInLocalFace(1,:,:)=NO_PART_DERIV
      basis%derivativeNumbersInLocalFace(0,:,:)=1
          
      !Set the line and face values
      SELECT CASE(basis%interpolationOrder(1))
      CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
        !Line 1
        basis%numberOfNodesInLocalLine(1)=2
        basis%nodeNumbersInLocalLine(1,1)=1
        basis%nodeNumbersInLocalLine(2,1)=2
        !Line 2
        basis%numberOfNodesInLocalLine(2)=2
        basis%nodeNumbersInLocalLine(1,2)=1
        basis%nodeNumbersInLocalLine(2,2)=3
        !Line 3
        basis%numberOfNodesInLocalLine(3)=2
        basis%nodeNumbersInLocalLine(1,3)=1
        basis%nodeNumbersInLocalLine(2,3)=4
        !Line 4
        basis%numberOfNodesInLocalLine(4)=2
        basis%nodeNumbersInLocalLine(1,4)=2
        basis%nodeNumbersInLocalLine(2,4)=3
        !Line 5
        basis%numberOfNodesInLocalLine(5)=2
        basis%nodeNumbersInLocalLine(1,5)=2
        basis%nodeNumbersInLocalLine(2,5)=4
        !Line 6
        basis%numberOfNodesInLocalLine(6)=2
        basis%nodeNumbersInLocalLine(1,6)=3
        basis%nodeNumbersInLocalLine(2,6)=4
        !Face 1
        basis%numberOfNodesInLocalFace(1)=3
        basis%nodeNumbersInLocalFace(1,1)=2
        basis%nodeNumbersInLocalFace(2,1)=3
        basis%nodeNumbersInLocalFace(3,1)=4
        !Face 2
        basis%numberOfNodesInLocalFace(2)=3
        basis%nodeNumbersInLocalFace(1,2)=1
        basis%nodeNumbersInLocalFace(2,2)=4
        basis%nodeNumbersInLocalFace(3,2)=3
        !Face 3
        basis%numberOfNodesInLocalFace(3)=3
        basis%nodeNumbersInLocalFace(1,3)=1
        basis%nodeNumbersInLocalFace(2,3)=2
        basis%nodeNumbersInLocalFace(3,3)=4
        !Face 4
        basis%numberOfNodesInLocalFace(4)=3
        basis%nodeNumbersInLocalFace(1,4)=1
        basis%nodeNumbersInLocalFace(2,4)=3
        basis%nodeNumbersInLocalFace(3,4)=2
      CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
        !Line 1
        basis%numberOfNodesInLocalLine(1)=3
        basis%nodeNumbersInLocalLine(1,1)=1
        basis%nodeNumbersInLocalLine(2,1)=5
        basis%nodeNumbersInLocalLine(3,1)=2
        !Line 2
        basis%numberOfNodesInLocalLine(2)=3
        basis%nodeNumbersInLocalLine(1,2)=1
        basis%nodeNumbersInLocalLine(2,2)=6
        basis%nodeNumbersInLocalLine(3,2)=3
        !Line 3
        basis%numberOfNodesInLocalLine(3)=3
        basis%nodeNumbersInLocalLine(1,3)=1
        basis%nodeNumbersInLocalLine(2,3)=7
        basis%nodeNumbersInLocalLine(3,3)=4
        !Line 4
        basis%numberOfNodesInLocalLine(4)=3
        basis%nodeNumbersInLocalLine(1,4)=2
        basis%nodeNumbersInLocalLine(2,4)=8
        basis%nodeNumbersInLocalLine(3,4)=3
        !Line 5
        basis%numberOfNodesInLocalLine(5)=3
        basis%nodeNumbersInLocalLine(1,5)=2
        basis%nodeNumbersInLocalLine(2,5)=10
        basis%nodeNumbersInLocalLine(3,5)=4
        !Line 6
        basis%numberOfNodesInLocalLine(6)=3
        basis%nodeNumbersInLocalLine(1,6)=3
        basis%nodeNumbersInLocalLine(2,6)=9
        basis%nodeNumbersInLocalLine(3,6)=4
        !Face 1
        basis%numberOfNodesInLocalFace(1)=6
        basis%nodeNumbersInLocalFace(1,1)=2
        basis%nodeNumbersInLocalFace(2,1)=3
        basis%nodeNumbersInLocalFace(3,1)=4
        basis%nodeNumbersInLocalFace(4,1)=8
        basis%nodeNumbersInLocalFace(5,1)=9
        basis%nodeNumbersInLocalFace(6,1)=10
        !Face 2
        basis%numberOfNodesInLocalFace(2)=6
        basis%nodeNumbersInLocalFace(1,2)=1
        basis%nodeNumbersInLocalFace(2,2)=4
        basis%nodeNumbersInLocalFace(3,2)=3
        basis%nodeNumbersInLocalFace(4,2)=7
        basis%nodeNumbersInLocalFace(5,2)=9
        basis%nodeNumbersInLocalFace(6,2)=6
        !Face 3
        basis%numberOfNodesInLocalFace(3)=6
        basis%nodeNumbersInLocalFace(1,3)=1
        basis%nodeNumbersInLocalFace(2,3)=2
        basis%nodeNumbersInLocalFace(3,3)=4
        basis%nodeNumbersInLocalFace(4,3)=5
        basis%nodeNumbersInLocalFace(5,3)=10
        basis%nodeNumbersInLocalFace(6,3)=7
        !Face 4
        basis%numberOfNodesInLocalFace(4)=6
        basis%nodeNumbersInLocalFace(1,4)=1
        basis%nodeNumbersInLocalFace(2,4)=3
        basis%nodeNumbersInLocalFace(3,4)=2
        basis%nodeNumbersInLocalFace(4,4)=6
        basis%nodeNumbersInLocalFace(5,4)=8
        basis%nodeNumbersInLocalFace(6,4)=5
      CASE(basis_CUBIC_INTERPOLATION_ORDER)
        !Line 1
        basis%numberOfNodesInLocalLine(1)=4
        basis%nodeNumbersInLocalLine(1,1)=1
        basis%nodeNumbersInLocalLine(2,1)=5
        basis%nodeNumbersInLocalLine(3,1)=6
        basis%nodeNumbersInLocalLine(4,1)=2
        !Line 2
        basis%numberOfNodesInLocalLine(2)=4
        basis%nodeNumbersInLocalLine(1,2)=1
        basis%nodeNumbersInLocalLine(2,2)=7
        basis%nodeNumbersInLocalLine(3,2)=8
        basis%nodeNumbersInLocalLine(4,2)=3
        !Line 3
        basis%numberOfNodesInLocalLine(3)=4
        basis%nodeNumbersInLocalLine(1,3)=1
        basis%nodeNumbersInLocalLine(2,3)=9
        basis%nodeNumbersInLocalLine(3,3)=10
        basis%nodeNumbersInLocalLine(4,3)=4
        !Line 4
        basis%numberOfNodesInLocalLine(4)=4
        basis%nodeNumbersInLocalLine(1,4)=2
        basis%nodeNumbersInLocalLine(2,4)=11
        basis%nodeNumbersInLocalLine(3,4)=12
        basis%nodeNumbersInLocalLine(4,4)=3
        !Line 5
        basis%numberOfNodesInLocalLine(5)=4
        basis%nodeNumbersInLocalLine(1,5)=2
        basis%nodeNumbersInLocalLine(2,5)=15
        basis%nodeNumbersInLocalLine(3,5)=16
        basis%nodeNumbersInLocalLine(4,5)=4
        !Line 6
        basis%numberOfNodesInLocalLine(6)=4
        basis%nodeNumbersInLocalLine(1,6)=3
        basis%nodeNumbersInLocalLine(2,6)=13
        basis%nodeNumbersInLocalLine(3,6)=14
        basis%nodeNumbersInLocalLine(4,6)=4
        !Face 1
        basis%numberOfNodesInLocalFace(1)=10
        basis%nodeNumbersInLocalFace(1,1)=2
        basis%nodeNumbersInLocalFace(2,1)=3
        basis%nodeNumbersInLocalFace(3,1)=4
        basis%nodeNumbersInLocalFace(4,1)=11
        basis%nodeNumbersInLocalFace(5,1)=12
        basis%nodeNumbersInLocalFace(6,1)=13
        basis%nodeNumbersInLocalFace(7,1)=14
        basis%nodeNumbersInLocalFace(8,1)=16
        basis%nodeNumbersInLocalFace(9,1)=15
        basis%nodeNumbersInLocalFace(10,1)=20
        !Face 2
        basis%numberOfNodesInLocalFace(2)=10
        basis%nodeNumbersInLocalFace(1,2)=1
        basis%nodeNumbersInLocalFace(2,2)=4
        basis%nodeNumbersInLocalFace(3,2)=3
        basis%nodeNumbersInLocalFace(4,2)=9
        basis%nodeNumbersInLocalFace(5,2)=10
        basis%nodeNumbersInLocalFace(6,2)=14
        basis%nodeNumbersInLocalFace(7,2)=13
        basis%nodeNumbersInLocalFace(8,2)=8
        basis%nodeNumbersInLocalFace(9,2)=7
        basis%nodeNumbersInLocalFace(10,2)=19
        !Face 3
        basis%numberOfNodesInLocalFace(3)=10
        basis%nodeNumbersInLocalFace(1,3)=1
        basis%nodeNumbersInLocalFace(2,3)=2
        basis%nodeNumbersInLocalFace(3,3)=4
        basis%nodeNumbersInLocalFace(4,3)=5
        basis%nodeNumbersInLocalFace(5,3)=6
        basis%nodeNumbersInLocalFace(6,3)=15
        basis%nodeNumbersInLocalFace(7,3)=16
        basis%nodeNumbersInLocalFace(8,3)=10
        basis%nodeNumbersInLocalFace(9,3)=9
        basis%nodeNumbersInLocalFace(10,3)=18
        !Face 4
        basis%numberOfNodesInLocalFace(4)=10
        basis%nodeNumbersInLocalFace(1,4)=1
        basis%nodeNumbersInLocalFace(2,4)=3
        basis%nodeNumbersInLocalFace(3,4)=2
        basis%nodeNumbersInLocalFace(4,4)=7
        basis%nodeNumbersInLocalFace(5,4)=8
        basis%nodeNumbersInLocalFace(6,4)=12
        basis%nodeNumbersInLocalFace(7,4)=11
        basis%nodeNumbersInLocalFace(8,4)=6
        basis%nodeNumbersInLocalFace(9,4)=5
        basis%nodeNumbersInLocalFace(10,4)=17
      CASE DEFAULT
        localError="The basis interpolation order of "//TRIM(NumberToVString(basis%interpolationOrder(1),"*",err,error))// &
          & " at xi coordinate index 1 is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set line and face directions and normals
      !Line 1
      basis%localLineXiDirection(1)=1
      basis%localLineXiNormals(1,1)=3
      basis%localLineXiNormals(2,1)=4
      basis%xiNormalsLocalLine(3,4)=1
      basis%xiNormalsLocalLine(4,3)=1
      !Line 2
      basis%localLineXiDirection(2)=1
      basis%localLineXiNormals(1,2)=2
      basis%localLineXiNormals(2,2)=4
      basis%xiNormalsLocalLine(2,4)=2
      basis%xiNormalsLocalLine(4,2)=2
      !Line 3
      basis%localLineXiDirection(3)=1
      basis%localLineXiNormals(1,3)=2
      basis%localLineXiNormals(2,3)=3
      basis%xiNormalsLocalLine(2,3)=3
      basis%xiNormalsLocalLine(3,2)=3
      !Line 4
      basis%localLineXiDirection(4)=2
      basis%localLineXiNormals(1,4)=1
      basis%localLineXiNormals(2,4)=4
      basis%xiNormalsLocalLine(1,4)=4
      basis%xiNormalsLocalLine(4,1)=4
      !Line 5
      basis%localLineXiDirection(5)=2
      basis%localLineXiNormals(1,5)=1
      basis%localLineXiNormals(2,5)=3
      basis%xiNormalsLocalLine(1,3)=5
      basis%xiNormalsLocalLine(3,1)=5
      !Line 6
      basis%localLineXiDirection(6)=3
      basis%localLineXiNormals(1,6)=1
      basis%localLineXiNormals(2,6)=2
      basis%xiNormalsLocalLine(1,2)=6
      basis%xiNormalsLocalLine(2,1)=6
      !Face 1
      basis%localFaceXiDirections(1,1)=2
      basis%localFaceXiDirections(2,1)=3
      basis%localFaceXiDirections(3,1)=4          
      basis%localFaceXiNormal(1)=1
      basis%xiNormalLocalFace(1)=1
      !Face 2
      basis%localFaceXiDirections(1,2)=1
      basis%localFaceXiDirections(2,2)=3
      basis%localFaceXiDirections(3,2)=4          
      basis%localFaceXiNormal(2)=2
      basis%xiNormalLocalFace(2)=2
      !Face 3
      basis%localFaceXiDirections(1,3)=1
      basis%localFaceXiDirections(2,3)=2
      basis%localFaceXiDirections(3,3)=4          
      basis%localFaceXiNormal(3)=3 
      basis%xiNormalLocalFace(3)=3
      !Face 4
      basis%localFaceXiDirections(1,4)=1
      basis%localFaceXiDirections(2,4)=2
      basis%localFaceXiDirections(3,4)=3          
      basis%localFaceXiNormal(4)=4
      basis%xiNormalLocalFace(4)=4          
    CASE DEFAULT
     localError="The number of xi directions of "//TRIM(NumberToVString(basis%numberOfXi,"*",err,error))// &
        & " for basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))// &
        & " is invalid. The number of xi directions should be >= 1 and <= 3."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    CALL Basis_QuadratureCreate(basis,err,error,*999)
              
    EXITS("Basis_SimplexBasisCreate")
    RETURN
999 ERRORSEXITS("Basis_SimplexBasisCreate",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_SimplexBasisCreate

  !
  !================================================================================================================================
  !

  !>Creates and initialises a simplex basis family that has already been allocated by Basis_CreateStart.
  !> \see BasisRoutines::Basis_SimplexBasisCreate
  SUBROUTINE Basis_SimplexFamilyCreate(basis,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,faceXi(2),faceXi2(2),localFaceIdx,localLineIdx,xiIdx,xiIdx2
    LOGICAL :: lineBasisDone,faceBasisDone
    TYPE(BasisType), POINTER :: newSubBasis
    TYPE(VARYING_STRING) :: dummyError

    NULLIFY(newSubBasis)

    ENTERS("Basis_SimplexFamilyCreate",err,error,*999)

    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    
    !Create the main (parent) basis
    CALL Basis_SimplexBasisCreate(basis,err,error,*999)
    
    IF(basis%numberOfXi>1) THEN
      !Create the line bases as sub-basis types
      ALLOCATE(basis%lineBases(basis%numberOfXi),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate basis line bases.",err,error,*999)
      ALLOCATE(basis%localLineBasis(basis%numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate local line basis.",err,error,*999)
      DO xiIdx=1,basis%numberOfXi
        lineBasisDone=.FALSE.
        NULLIFY(newSubBasis)
        DO xiIdx2=1,xiIdx-1
          IF(basis%interpolationXi(xiIdx2)==basis%interpolationXi(xiIdx).AND. &
            basis%quadrature%numberOfGaussXi(xiIdx2)==basis%quadrature%numberOfGaussXi(xiIdx)) THEN
            lineBasisDone=.TRUE.
            EXIT
          ENDIF
        ENDDO !xiIdx2
        IF(lineBasisDone) THEN
          basis%lineBases(xiIdx)%ptr=>basis%lineBases(xiIdx2)%ptr
        ELSE
          !Create the new sub-basis
          CALL Basis_SubBasisCreate(basis,1,[xiIdx],newSubBasis,err,error,*999)
          !Fill in the basis information
          CALL Basis_SimplexBasisCreate(newSubBasis,err,error,*999)
          basis%lineBases(xiIdx)%ptr=>newSubBasis
        ENDIF
      ENDDO !localLineIdx
      DO localLineIdx=1,basis%numberOfLocalLines
!!\TODO Until we allow for tensor products of simplex bases each line basis will be the same and equal to the first one.
        basis%localLineBasis(localLineIdx)%ptr=>basis%lineBases(1)%ptr       
      ENDDO !localLineIdx
      IF(basis%numberOfXi>2) THEN
        !Create the face bases as sub-basis types
        ALLOCATE(basis%faceBases(basis%numberOfXi),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate basis face bases.",err,error,*999)
        ALLOCATE(basis%localFaceBasis(basis%numberOfLocalFaces),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate local face basis.",err,error,*999)
        DO xiIdx=1,basis%numberOfXi
          !Determine the face xi directions that lie in this xi direction
          faceXi(1)=OTHER_XI_DIRECTIONS3(xiIdx,2,1)
          faceXi(2)=OTHER_XI_DIRECTIONS3(xiIdx,3,1)
          faceBasisDone=.FALSE.
          NULLIFY(newSubBasis)
          DO xiIdx2=1,xiIdx-1
            !Determine the face xi directions that lie in this xi direction
            faceXi2(1)=OTHER_XI_DIRECTIONS3(xiIdx2,2,1)
            faceXi2(2)=OTHER_XI_DIRECTIONS3(xiIdx2,3,1)
            IF(basis%interpolationXi(faceXi2(1))==basis%interpolationXi(faceXi(1)).AND. &
              & basis%interpolationXi(faceXi2(2))==basis%interpolationXi(faceXi(2)).AND. &
              & basis%quadrature%numberOfGaussXi(faceXi2(1))==basis%quadrature%numberOfGaussXi(faceXi(1)).AND. &
              & basis%quadrature%numberOfGaussXi(faceXi2(2))==basis%quadrature%numberOfGaussXi(faceXi(2))) THEN
              faceBasisDone=.TRUE.
              EXIT
            ENDIF
          ENDDO !xiIdx2
          IF(faceBasisDone) THEN
            basis%faceBases(xiIdx)%ptr=>basis%faceBases(xiIdx2)%ptr
          ELSE
            !Create the new sub-basis
            CALL Basis_SubBasisCreate(basis,2,[faceXi(1),faceXi(2)],newSubBasis,err,error,*999)
            !Fill in the basis information
            CALL Basis_SimplexBasisCreate(newSubBasis,err,error,*999)
            newSubBasis%lineBases(1)%ptr=>basis%lineBases(faceXi(1))%ptr
            newSubBasis%lineBases(2)%ptr=>basis%lineBases(faceXi(2))%ptr
            basis%faceBases(xiIdx)%ptr=>newSubBasis
            ALLOCATE(newSubBasis%localLineBasis(newSubBasis%numberOfLocalLines),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate new sub basis local line basis.",err,error,*999)
             DO localLineIdx=1,newSubBasis%numberOfLocalLines
              newSubBasis%localLineBasis(localLineIdx)%ptr=>newSubBasis%lineBases(1)%ptr
            ENDDO !localFaceIdx
          ENDIF
        ENDDO !xiIdx
        DO localFaceIdx=1,basis%numberOfLocalFaces
!!\TODO Until we allow for tensor products of simplex bases each face basis will be the same and equal to the first one.
          basis%localFaceBasis(localFaceIdx)%ptr=>basis%faceBases(1)%ptr
        ENDDO !localFaceIdx
      ENDIF
    ENDIF
    
    EXITS("Basis_SimplexFamilyCreate")
    RETURN
999 IF(ASSOCIATED(newSubBasis)) CALL Basis_FamilyDestroy(newSubBasis%basisFunctions,newSubBasis%userNumber, &
      & newSubBasis%familyNumber,dummyErr,dummyError,*998)
998 ERRORSEXITS("Basis_SimplexFamilyCreate",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_SimplexFamilyCreate

  !
  !================================================================================================================================
  !

  !>Evaluates a simplex basis function and its derivatives with respect to external \f$\mathbf{\xi}\f$ coordinates.
  !>For Simplex line elements there are two area coordinates which are a function of \f$\xi_1\f$ : \f$L_1 = 1 - \xi_1\f$ and
  !>\f$L_2 = \xi_1 - 1\f$.The derivatives wrt to external coordinates are then given by \f$\frac{\partial\mathbf{N}}{\partial\xi_1}=
  !>\frac{\partial\mathbf(x)}{\partial L_2}-\frac{\partial \mathbf{N}}{\partial L_1}\f$ and \f$\frac{\partial^2\mathbf{N}}{
  !>\partial \xi_1^2} = \frac{\partial^2\mathbf{N}}{\partial L_1^2}-2\frac{\partial^2\mathbf{N}}{\partial L_1 \partial L_2}+
  !>\frac{\partial^2\mathbf{N}}{\partial L_2^2}\f$.
  !>For Simplex triangle elements there are three area coordinates which are a function of \f$\xi_1\f$ and
  !>\f$\xi_2\f$ : \f$L_1 = 1 - \xi_1\f$, \f$L_2 = 1 - \xi_2\f$ and \f$L_3=\xi_1 + \xi_2 - 1\f$. The derivatives wrt to external
  !>coordinates are then given by \f$\frac{\partial \mathbf{N}}{\partial\xi_1}=\frac{\partial\mathbf(N)}{\partial L_3}-
  !>\frac{\partial \mathbf{N}}{\partial L_1}\f$, \f$\frac{\partial \mathbf{N}}{\partial\xi_2}=\frac{\partial\mathbf(x)}{
  !>\partial L_3}-\frac{\partial \mathbf{N}}{\partial L_2}\f$, \f$\frac{\partial^2\mathbf{N}}{\partial \xi_1^2} =
  !>\frac{\partial^2\mathbf{N}}{\partial L_1^2}-2\frac{\partial^2\mathbf{N}}{\partial L_1 \partial L_3}+
  !>\frac{\partial^2\mathbf{N}}{\partial L_3^2}\f$, \f$\frac{\partial^2\mathbf{N}}{\partial \xi_2^2} =
  !>\frac{\partial^2\mathbf{N}}{\partial L_2^2}-2\frac{\partial^2\mathbf{N}}{\partial L_2 \partial L_3}+
  !>\frac{\partial^2\mathbf{N}}{\partial L_3^2}\f$ and \f$\frac{\partial^2\mathbf{N}}{\partial \xi_1 \partial \xi_2} =
  !>\frac{\partial^2\mathbf{N}}{\partial L_3^2}-\frac{\partial^2\mathbf{N}}{\partial L_1 \partial L_3}-
  !>\frac{\partial^2\mathbf{N}}{\partial L_2 \partial L_3}+\frac{\partial^2\mathbf{N}}{\partial L_1 \partial L_2}\f$.
  !>For Simplex tetrahedral elements there are four area coordinates which are a function of \f$\xi_1\f$,\f$\xi_2\f$ and
  !>\f$\xi_3\f$ : \f$L_1 = 1 - \xi_1\f$, \f$L_2 = 1 - \xi_2\f$, \f$L_3 = 1 - \xi_3\f$ and
  !>\f$L_4 = \xi_1 + \xi_2 + \xi_3 - 1\f$. The derivatives wrt to external coordinates are then given by
  !>\f$\frac{\partial \mathbf{N}}{\partial\xi_1}=\frac{\partial\mathbf(x)}{\partial L_4}-
  !>\frac{\partial \mathbf{N}}{\partial L_1}\f$,
  !>\f$\frac{\partial \mathbf{N}}{\partial\xi_2}=\frac{\partial\mathbf(x)}{\partial L_4}-
  !>\frac{\partial \mathbf{N}}{\partial L_2}\f$,
  !>\f$\frac{\partial \mathbf{N}}{\partial\xi_3}=\frac{\partial\mathbf(x)}{\partial L_4}-
  !>\frac{\partial \mathbf{N}}{\partial L_3}\f$,
  !>\f$\frac{\partial^2\mathbf{N}}{\partial \xi_1^2} = \frac{\partial^2\mathbf{N}}{\partial L_1^2}-
  !>2\frac{\partial^2\mathbf{N}}{\partial L_1 \partial L_4}+\frac{\partial^2\mathbf{N}}{\partial L_4^2}\f$,
  !>\f$\frac{\partial^2\mathbf{N}}{\partial \xi_2^2} = \frac{\partial^2\mathbf{N}}{\partial L_2^2}-
  !>2\frac{\partial^2\mathbf{N}}{\partial L_2 \partial L_4}+\frac{\partial^2\mathbf{N}}{\partial L_4^2}\f$
  !>\f$\frac{\partial^2\mathbf{N}}{\partial \xi_3^2} = \frac{\partial^2\mathbf{N}}{\partial L_3^2}-
  !>2\frac{\partial^2\mathbf{N}}{\partial L_3 \partial L_4}+\frac{\partial^2\mathbf{N}}{\partial L_4^2}\f$,
  !>\f$\frac{\partial^2\mathbf{N}}{\partial\xi_1\partial \xi_2}=\frac{\partial^2\mathbf{N}}{\partial L_4^2}-
  !>\frac{\partial^2\mathbf{N}}{\partial L_1 \partial L_4}-\frac{\partial^2\mathbf{N}}{\partial L_2 \partial L_4}+
  !>\frac{\partial^2\mathbf{N}}{\partial L_1 \partial L_2}\f$,
  !>\f$\frac{\partial^2\mathbf{N}}{\partial\xi_1\partial\xi_3}=\frac{\partial^2\mathbf{N}}{\partial L_4^2}-
  !>\frac{\partial^2\mathbf{N}}{\partial L_1 \partial L_4}-\frac{\partial^2\mathbf{N}}{\partial L_3 \partial L_4}+
  !>\frac{\partial^2\mathbf{N}}{\partial L_1 \partial L_3}\f$,
  !>\f$\frac{\partial^2\mathbf{N}}{\partial\xi_2\partial\xi_3}=\frac{\partial^2\mathbf{N}}{\partial L_4^2}-
  !>\frac{\partial^2\mathbf{N}}{\partial L_2 \partial L_4}-\frac{\partial^2\mathbf{N}}{\partial L_3 \partial L_4}+
  !>\frac{\partial^2\mathbf{N}}{\partial L_2 \partial L_3}\f$ and
  !>\f$\frac{\partial^3\mathbf{N}}{\partial \xi_1 \partial \xi_2 \partial \xi_3} = \frac{\partial^3\mathbf{N}}{\partial L_4^3}-
  !>\frac{\partial^3\mathbf{N}}{\partial L_1 \partial L_4^2}-\frac{\partial^3\mathbf{N}}{\partial L_2 \partial L_4^2}-
  !>\frac{\partial^3\mathbf{N}}{\partial L_3 \partial L_4^2}+\frac{\partial^3\mathbf{N}}{\partial L_1 \partial 2 \partial L_4}+
  !>\frac{\partial^3\mathbf{N}}{\partial L_1 \partial L_3 \partial L_4}+\frac{\partial^3\mathbf{N}}{\partial L_2 \partial L_3
  !>\partial L_4}-\frac{\partial^3\mathbf{N}}{\partial L_1 \partial L_2 \partial L_3}\f$.
  FUNCTION Basis_SimplexBasisEvaluate(basis,localNodeNumber,partialDerivativeIndex,xl,err,error)
    
    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis function to evaluate.
    INTEGER(INTG), INTENT(IN) :: localNodeNumber !<The node number defines the actual basis function to evaluate.
    INTEGER(INTG), INTENT(IN) :: partialDerivativeIndex !<The partial derivative index in Xi space of the basis to evaluate.
    REAL(DP), INTENT(IN) :: xl(:) !<xl(xiCoordIdx). The area coordinates to evaluate the basis function at.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: Basis_SimplexBasisEvaluate !<On return the evaluated basis function
    !Local variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Basis_SimplexBasisEvaluate",err,error,*999)
    
    Basis_SimplexBasisEvaluate=0.0_DP

    CALL Basis_AssertIsSimplexBasis(basis,err,error,*999)
    
    SELECT CASE(basis%numberOfXi)
    CASE(1)
      SELECT CASE(partialDerivativeIndex)
      CASE(NO_PART_DERIV)
        Basis_SimplexBasisEvaluate= &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,NO_PART_DERIV,xl,err,error)
        IF(err/=0) GOTO 999
      CASE(PART_DERIV_S1)
        Basis_SimplexBasisEvaluate= &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S2,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S1,xl,err,error)
        IF(err/=0) GOTO 999
      CASE(PART_DERIV_S1_S1)
        Basis_SimplexBasisEvaluate= &
          Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S1_S1,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & 2.0_DP*Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S1_S2,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate+ &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S2_S2,xl,err,error)
        IF(err/=0) GOTO 999
      CASE DEFAULT
        localError="The specified partial derivative index of "//TRIM(NumberToVString(partialDerivativeIndex,"*",err,error))// &
          & " is invalid for a Simplex line basis."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(2)
      SELECT CASE(partialDerivativeIndex)
      CASE(NO_PART_DERIV)
        Basis_SimplexBasisEvaluate= &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,NO_PART_DERIV,xl,err,error)
        IF(err/=0) GOTO 999
      CASE(PART_DERIV_S1)
        Basis_SimplexBasisEvaluate= &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S3,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S1,xl,err,error)
        IF(err/=0) GOTO 999
      CASE(PART_DERIV_S1_S1)
        Basis_SimplexBasisEvaluate= &
          Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S1_S1,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & 2.0_DP*Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S1_S3,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate+ &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S3_S3,xl,err,error)
        IF(err/=0) GOTO 999
      CASE(PART_DERIV_S2)
        Basis_SimplexBasisEvaluate= &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S3,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S2,xl,err,error)
        IF(err/=0) GOTO 999
     CASE(PART_DERIV_S2_S2)
        Basis_SimplexBasisEvaluate= &
          Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S2_S2,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & 2.0_DP*Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S2_S3,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate+ &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S3_S3,xl,err,error)
        IF(err/=0) GOTO 999
     CASE(PART_DERIV_S1_S2)
        Basis_SimplexBasisEvaluate= &
          Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S3_S3,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S1_S3,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S2_S3,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate+ &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S1_S2,xl,err,error)
        IF(err/=0) GOTO 999
      CASE DEFAULT
        localError="The specified partial derivative index of "//TRIM(NumberToVString(partialDerivativeIndex,"*",err,error))// &
          & " is invalid for a Simplex triangle basis."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(3)
      SELECT CASE(partialDerivativeIndex)
      CASE(NO_PART_DERIV)
        Basis_SimplexBasisEvaluate= &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,NO_PART_DERIV,xl,err,error)
        IF(err/=0) GOTO 999
      CASE(PART_DERIV_S1)
        Basis_SimplexBasisEvaluate= &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S4,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S1,xl,err,error)
        IF(err/=0) GOTO 999
      CASE(PART_DERIV_S1_S1)
        Basis_SimplexBasisEvaluate= &
          Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S1_S1,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & 2.0_DP*Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S1_S4,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate+ &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S4_S4,xl,err,error)
        IF(err/=0) GOTO 999
      CASE(PART_DERIV_S2)
        Basis_SimplexBasisEvaluate= &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S4,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S2,xl,err,error)
        IF(err/=0) GOTO 999
      CASE(PART_DERIV_S2_S2)
        Basis_SimplexBasisEvaluate= &
          Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S2_S2,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & 2.0_DP*Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S2_S4,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate+ &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S4_S4,xl,err,error)
        IF(err/=0) GOTO 999
      CASE(PART_DERIV_S1_S2)
        Basis_SimplexBasisEvaluate= &
          Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S4_S4,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S1_S4,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S2_S4,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate+ &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S1_S2,xl,err,error)
        IF(err/=0) GOTO 999
      CASE(PART_DERIV_S3)
        Basis_SimplexBasisEvaluate= &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S4,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S3,xl,err,error)
        IF(err/=0) GOTO 999
      CASE(PART_DERIV_S3_S3)
        Basis_SimplexBasisEvaluate= &
          Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S3_S3,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & 2.0_DP*Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S3_S4,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate+ &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S4_S4,xl,err,error)
        IF(err/=0) GOTO 999
      CASE(PART_DERIV_S1_S3)
        Basis_SimplexBasisEvaluate= &
          Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S4_S4,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S1_S4,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S3_S4,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate+ &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S1_S3,xl,err,error)
        IF(err/=0) GOTO 999
      CASE(PART_DERIV_S2_S3)
        Basis_SimplexBasisEvaluate= &
          Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S4_S4,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S2_S4,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S3_S4,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate+ &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S2_S3,xl,err,error)
        IF(err/=0) GOTO 999
      CASE(PART_DERIV_S1_S2_S3)
        Basis_SimplexBasisEvaluate= &
          Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S4_S4_S4,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S1_S4_S4,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S2_S4_S4,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S3_S4_S4,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate+ &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S1_S2_S4,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate+ &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S1_S3_S4,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate+ &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S2_S3_S4,xl,err,error)
        IF(err/=0) GOTO 999
        Basis_SimplexBasisEvaluate=Basis_SimplexBasisEvaluate- &
          & Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,PART_DERIV_S1_S2_S3,xl,err,error)
        IF(err/=0) GOTO 999
      CASE DEFAULT
        localError="The specified partial derivative index of "//TRIM(NumberToVString(partialDerivativeIndex,"*",err,error))// &
          & " is invalid for a Simplex tetrahedra basis."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="Invalid number of Xi coordinates. The number of xi coordinates for this basis is "// &
        & TRIM(NumberToVString(basis%numberOfXi,"*",err,error))//" which should be between 1 and 3."
      CALL FlagError(localError,err,error,*999)
    END SELECT
     
    EXITS("Basis_SimplexBasisEvaluate")
    RETURN
999 ERRORSEXITS("Basis_SimplexBasisEvaluate",err,error)
    RETURN
    
  END FUNCTION Basis_SimplexBasisEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates partial derivatives of a simplex basis function with respect to area coordinates.
  FUNCTION Basis_SimplexBasisDerivativeEvaluate(basis,localNodeNumber,partialDerivativeIndex,xl,err,error)
    
    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis function to evaluate.
    INTEGER(INTG), INTENT(IN) :: localNodeNumber !<The node number defines the actual basis function to evaluate.
    INTEGER(INTG), INTENT(IN) :: partialDerivativeIndex !<The partial derivative index in area coordinates of the basis to evaluate.
    REAL(DP), INTENT(IN) :: xl(:) !<xl(xiCoordIdx). The area coordinates to evaluate the basis function at.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: Basis_SimplexBasisDerivativeEvaluate !<On return the evaluated basis function
    !Local variables
    INTEGER(INTG) :: xiCoordIdx
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Basis_SimplexBasisDerivativeEvaluate",err,error,*999)
    
    Basis_SimplexBasisDerivativeEvaluate=1.0_DP
    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    
    DO xiCoordIdx=1,basis%numberOfXiCoordinates       
      SELECT CASE(basis%interpolationOrder(xiCoordIdx))
      CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
        Basis_SimplexBasisDerivativeEvaluate=Basis_SimplexBasisDerivativeEvaluate* &
          & Simplex_LinearEvaluate(basis%nodePositionIndex(localNodeNumber,xiCoordIdx), &
          & PARTIAL_DERIVATIVE_INDEX(partialDerivativeIndex,xiCoordIdx),xl(xiCoordIdx),err,error)
        IF(err/=0) GOTO 999
      CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
        Basis_SimplexBasisDerivativeEvaluate=Basis_SimplexBasisDerivativeEvaluate* &
          & Simplex_QuadraticEvaluate(basis%nodePositionIndex(localNodeNumber,xiCoordIdx), &
          & PARTIAL_DERIVATIVE_INDEX(partialDerivativeIndex,xiCoordIdx),xl(xiCoordIdx),err,error)
        IF(err/=0) GOTO 999
      CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
        Basis_SimplexBasisDerivativeEvaluate=Basis_SimplexBasisDerivativeEvaluate* &
          & Simplex_CubicEvaluate(basis%nodePositionIndex(localNodeNumber,xiCoordIdx), &
          & PARTIAL_DERIVATIVE_INDEX(partialDerivativeIndex,xiCoordIdx),xl(xiCoordIdx),err,error)
        IF(err/=0) GOTO 999
      CASE DEFAULT
        localError="Interpolation order value "//TRIM(NumberToVString(basis%interpolationOrder(xiCoordIdx),"*",err,error))// &
          & " for xi coordinate direction "//TRIM(NumberToVString(xiCoordIdx,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      IF(err/=0) GOTO 999
    ENDDO !xiCoordIdx
   
    EXITS("Basis_SimplexBasisDerivativeEvaluate")
    RETURN
999 ERRORSEXITS("Basis_SimplexBasisDerivativeEvaluate",err,error)
    RETURN
    
  END FUNCTION Basis_SimplexBasisDerivativeEvaluate

  !
  !================================================================================================================================
  !

  !>Creates a sub-basis on a parent basis.
  SUBROUTINE Basis_SubBasisCreate(parentBasis,numberOfXi,xiDirections,subBasis,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: parentBasis !<A pointer to the parent basis
    INTEGER(INTG), INTENT(IN) :: numberOfXi !<The number of Xi directions to create
    INTEGER(INTG), INTENT(IN) :: xiDirections(:) !<xiDirections(xiIdx). Gives the ii direction indices of the parent basis which are used to create the sub-basis
    TYPE(BasisType), POINTER :: subBasis !<On return, a pointer to the created sub-basis. The pointer must be NULL on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: basisIdx,xiIdx,numberCollapsed,numberEndCollapsed
    TYPE(BasisType), POINTER :: newSubBasis  
    TYPE(BasisPtrType), ALLOCATABLE :: newSubBases(:)
    TYPE(VARYING_STRING) :: localError
    
    NULLIFY(newSubBasis)
    
    ENTERS("Basis_SubBasisCreate",err,error,*999)

    IF(.NOT.ASSOCIATED(parentBasis)) CALL FlagError("Parent basis is not associated.",err,error,*999)
    IF(ASSOCIATED(subBasis)) CALL FlagError("The sub-basis is already associated.",err,error,*999)    
    IF(numberOfXi<=0.OR.numberOfXi>3) THEN
      localError="The number of xi directions specified of "//TRIM(NumberToVString(numberOfXi,"*",err,error))// &
        & " is invalid. The number of xi directions must be between 1 and 3."      
      CALL FlagError(localError,err,error,*999)
    ENDIF                
    IF(SIZE(xiDirections,1)/=numberOfXi) &
      & CALL FlagError("The size of the xi directions array must be the same as the number of xi directions",err,error,*999)
    IF(ANY(xiDirections<1).OR.ANY(xiDirections>3)) CALL FlagError("Invalid xi directions specified.",err,error,*999)

    CALL Basis_Initialise(newSubBasis,err,error,*999)
    newSubBasis%userNumber=parentBasis%userNumber
    newSubBasis%globalNumber=parentBasis%globalNumber
    newSubBasis%familyNumber=parentBasis%numberOfSubBases+1
    newSubBasis%basisFunctions=>parentBasis%basisFunctions
    newSubBasis%parentBasis=>parentBasis
    newSubBasis%numberOfXi=numberOfXi
    newSubBasis%type=parentBasis%type
    ALLOCATE(newSubBasis%interpolationXi(numberOfXi),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate sub-basis interpolation xi.",err,error,*999)
    ALLOCATE(newSubBasis%collapsedXi(numberOfXi),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate sub-basis collapsed xi.",err,error,*999)        
    numberCollapsed=0
    numberEndCollapsed=0
    DO xiIdx=1,numberOfXi
      newSubBasis%interpolationXi(xiIdx)=parentBasis%interpolationXi(xiDirections(xiIdx))
      newSubBasis%collapsedXi(xiIdx)=parentBasis%collapsedXi(xiDirections(xiIdx))
      IF(newSubBasis%collapsedXi(xiIdx)==BASIS_XI_COLLAPSED) THEN
        numberCollapsed=numberCollapsed+1
      ELSE IF(newSubBasis%collapsedXi(xiIdx)==BASIS_COLLAPSED_AT_XI0.OR.&
        & newSubBasis%collapsedXi(xiIdx)==BASIS_COLLAPSED_AT_XI1) THEN
        numberEndCollapsed=numberEndCollapsed+1
      ENDIF
    ENDDO !xiIdx
    IF(numberCollapsed==0.OR.numberEndCollapsed==0) newSubBasis%collapsedXi(1:numberOfXi)=BASIS_NOT_COLLAPSED
    NULLIFY(newSubBasis%quadrature%basis)
    CALL Basis_QuadratureInitialise(newSubBasis,err,error,*999)
    newSubBasis%quadrature%type=parentBasis%quadrature%type
    DO xiIdx=1,numberOfXi
      newSubBasis%quadrature%numberOfGaussXi(xiIdx)=parentBasis%quadrature%numberOfGaussXi(xiDirections(xiIdx))
    ENDDO !xiIdx
    newSubBasis%quadrature%gaussOrder=parentBasis%quadrature%gaussOrder
    newSubBasis%basisFinished=.TRUE.
    IF(numberOfXi>1) THEN
      ALLOCATE(newSubBasis%lineBases(numberOfXi),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate sub-basis line bases.",err,error,*999)
    ENDIF
    !Add the new sub-basis to the list of sub-bases in the parent basis
    ALLOCATE(newSubBases(parentBasis%numberOfSubBases+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new sub-bases.",err,error,*999)
    DO basisIdx=1,parentBasis%numberOfSubBases
      newSubBases(basisIdx)%ptr=>parentBasis%subBases(basisIdx)%ptr
    ENDDO !basisIdx
    newSubBases(parentBasis%numberOfSubBases+1)%ptr=>newSubBasis
    parentBasis%numberOfSubBases=parentBasis%numberOfSubBases+1
    CALL MOVE_ALLOC(newSubBases,parentBasis%subBases)
    subBasis=>newSubBasis
    
    EXITS("Basis_SubBasisCreate")
    RETURN
999 ERRORSEXITS("Basis_SubBasisCreate",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_SubBasisCreate
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the type for a basis is identified by a a pointer. \see OpenCMISS::Iron::cmfe_Basis_TypeSet
  SUBROUTINE Basis_TypeSet(basis,type,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to set
    INTEGER(INTG), INTENT(IN) :: type !<The type of the basis to be set. \see BasisRoutines_BasisTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Basis_TypeSet",err,error,*999)

    CALL Basis_AssertNotFinished(basis,err,error,*999)
    
    SELECT CASE(type)
    CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
      basis%type=BASIS_LAGRANGE_HERMITE_TP_TYPE
    CASE(BASIS_SIMPLEX_TYPE)
      !Reset the quadrature
      CALL Basis_QuadratureFinalise(basis,err,error,*999)
      !Change the default parameters for the old basis
      basis%type=BASIS_SIMPLEX_TYPE
      basis%interpolationXi(1:basis%numberOfXi)=BASIS_LINEAR_SIMPLEX_INTERPOLATION
      NULLIFY(basis%quadrature%basis)
      CALL Basis_QuadratureInitialise(basis,err,error,*999)
    CASE DEFAULT
      localError="Basis type "//TRIM(NumberToVString(type,"*",err,error))//" is invalid or not implemented."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Basis_TypeSet")
    RETURN
999 ERRORSEXITS("Basis_TypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_TypeSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the collapsed xi flags for a basis is identified by a a pointer. \see OpenCMISS::Iron::cmfe_Basis_CollapsedXiSet
  SUBROUTINE Basis_CollapsedXiSet(basis,collapsedXi,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: collapsedXi(:) !<collapsedXi(xiIdx). The collapse parameter for each Xi direction. \see BasisRoutinesS_XiCollapse
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: xiIdx1,xiIdx2,xiIdx3,numberCollapsed,collapsedXiDirection(3)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Basis_CollapsedXiSet",err,error,*999)

    CALL Basis_AssertNotFinished(basis,err,error,*999)
    CALL Basis_AssertIsLHTPBasis(basis,err,error,*999)    
    IF(basis%numberOfXi<=1) CALL FlagError("Can not collapse a basis with only 1 xi direction",err,error,*999) 
    IF(SIZE(collapsedXi,1)==basis%numberOfXi) THEN
      localError="The size of the xi collapsed array ("// &
        & TRIM(NumberToVString(SIZE(collapsedXi,1),"*",err,error))//") does not match the number of xi directions ("// &
        & TRIM(NumberToVString(basis%numberOfXi,"*",err,error))//") for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    numberCollapsed=0
    DO xiIdx1=1,basis%numberOfXi
      SELECT CASE(collapsedXi(xiIdx1))
      CASE(BASIS_XI_COLLAPSED)
        numberCollapsed=numberCollapsed+1
        collapsedXiDirection(numberCollapsed)=xiIdx1
      CASE(BASIS_COLLAPSED_AT_XI0,BASIS_COLLAPSED_AT_XI1,BASIS_NOT_COLLAPSED)
        !Do nothing
      CASE DEFAULT
        localError="Collapsed xi value "//TRIM(NumberToVString(collapsedXi(xiIdx1),"*",err,error))// &
          & " in xi direction "//TRIM(NumberToVString(xiIdx1,"*",err,error))//" is invalid"
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !xiIdx1
    IF(numberCollapsed>0) THEN
      IF(numberCollapsed<basis%numberOfXi) THEN
        IF(basis%numberOfXi==2) THEN
          !Two dimensional collapsed basis
          xiIdx1=collapsedXiDirection(1)
          xiIdx2=OTHER_XI_DIRECTIONS2(xiIdx1)
          IF(collapsedXi(xiIdx2)==BASIS_COLLAPSED_AT_XI0) THEN
            IF(basis%interpolationXi(xiIdx2)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
              & basis%interpolationXi(xiIdx2)=BASIS_QUADRATIC1_HERMITE_INTERPOLATION
          ELSE IF(collapsedXi(xiIdx2)==BASIS_COLLAPSED_AT_XI1) THEN
            IF(basis%interpolationXi(xiIdx2)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
              & basis%interpolationXi(xiIdx2)=BASIS_QUADRATIC2_HERMITE_INTERPOLATION
          ELSE
            localError="Invalid collapsing of a two dimensional basis. Xi direction "// &
              & TRIM(NumberToVString(xiIdx1,"*",err,error))//" is collapsed so xi direction "// &
              & TRIM(NumberToVString(xiIdx2,"*",err,error))//" must be collapsed at an end"
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          !Three dimensional collapsed basis
          IF(numberCollapsed==1) THEN
            !One collapse - wedge element
            xiIdx1=collapsedXiDirection(1)
            xiIdx2=OTHER_XI_DIRECTIONS3(xiIdx1,2,1)
            xiIdx3=OTHER_XI_DIRECTIONS3(xiIdx1,3,1)
            IF(collapsedXi(xiIdx2)==BASIS_NOT_COLLAPSED) THEN
              IF(collapsedXi(xiIdx3)==BASIS_COLLAPSED_AT_XI0) THEN
                IF(basis%interpolationXi(xiIdx3)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
                  & basis%interpolationXi(xiIdx3)=BASIS_QUADRATIC1_HERMITE_INTERPOLATION
              ELSE IF(collapsedXi(xiIdx3)==BASIS_COLLAPSED_AT_XI1) THEN
                IF(basis%interpolationXi(xiIdx3)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
                  & basis%interpolationXi(xiIdx3)=BASIS_QUADRATIC2_HERMITE_INTERPOLATION
              ELSE
                localError="Invalid collapsing of a three dimensional basis. Xi direction "// &
                  & TRIM(NumberToVString(xiIdx1,"*",err,error))//" is collapsed and xi direction "// &
                  & TRIM(NumberToVString(xiIdx2,"*",err,error))//" is not collapsed so xi direction "// &
                  & TRIM(NumberToVString(xiIdx3,"*",err,error))//" must be collapsed at an end"
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ELSE IF(collapsedXi(xiIdx3)==BASIS_NOT_COLLAPSED) THEN
              IF(collapsedXi(xiIdx2)==BASIS_COLLAPSED_AT_XI0) THEN
                IF(basis%interpolationXi(xiIdx2)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
                  & basis%interpolationXi(xiIdx2)=BASIS_QUADRATIC1_HERMITE_INTERPOLATION
              ELSE IF(collapsedXi(xiIdx2)==BASIS_COLLAPSED_AT_XI1) THEN
                IF(basis%interpolationXi(xiIdx2)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
                  & basis%interpolationXi(xiIdx2)=BASIS_QUADRATIC2_HERMITE_INTERPOLATION
              ELSE
                localError="Invalid collapsing of a three dimensional basis. Xi direction "// &
                  & TRIM(NumberToVString(xiIdx1,"*",err,error))//" is collapsed and xi direction "// &
                  & TRIM(NumberToVString(xiIdx3,"*",err,error))//" is not collapsed so xi direction "// &
                  & TRIM(NumberToVString(xiIdx2,"*",err,error))//" must be collapsed at an end"
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ELSE
              localError="Invalid collapsing of a three dimensional basis. Xi direction "// &
                & TRIM(NumberToVString(xiIdx1,"*",err,error))//" is collapsed so one of xi directions "// &
                & TRIM(NumberToVString(xiIdx2,"*",err,error))//" or "// &
                & TRIM(NumberToVString(xiIdx3,"*",err,error))//" must be collapsed at an end"
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            !Two collapses - pyramid element
            xiIdx1=collapsedXiDirection(1)
            xiIdx2=collapsedXiDirection(2)
            xiIdx3=OTHER_XI_DIRECTIONS3(xiIdx1,xiIdx2,2)
            IF(collapsedXi(xiIdx3)==BASIS_COLLAPSED_AT_XI0) THEN
              IF(basis%interpolationXi(xiIdx3)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
                & basis%interpolationXi(xiIdx3)=BASIS_QUADRATIC1_HERMITE_INTERPOLATION
            ELSE IF(collapsedXi(xiIdx3)==BASIS_COLLAPSED_AT_XI1) THEN
              IF(basis%interpolationXi(xiIdx3)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
                & basis%interpolationXi(xiIdx3)=BASIS_QUADRATIC2_HERMITE_INTERPOLATION
            ELSE
              localError="Invalid collapsing of a three dimensional basis. Xi directions "// &
                & TRIM(NumberToVString(xiIdx1,"*",err,error))//" and "// &
                & TRIM(NumberToVString(xiIdx2,"*",err,error))//" are collapsed so xi direction "// &
                & TRIM(NumberToVString(xiIdx3,"*",err,error))//" must be collapsed at an end"
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDIF
        ENDIF
      ELSE
        localError="Invalid collapsing of basis. The number of collapsed directions ("// &
          & TRIM(NumberToVString(numberCollapsed,"*",err,error))// &
          & ") must be less than the number of xi directions ("// &
          & TRIM(NumberToVString(basis%numberOfXi,"*",err,error))//")"
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      !No collapses in any xi direction - Reset interpolation_xi if necessary
      DO xiIdx1=1,basis%numberOfXi
        IF(basis%interpolationXi(xiIdx1)==BASIS_QUADRATIC1_HERMITE_INTERPOLATION.OR. &
          & basis%interpolationXi(xiIdx1)==BASIS_QUADRATIC2_HERMITE_INTERPOLATION) THEN
          basis%interpolationXi(xiIdx1)=BASIS_CUBIC_HERMITE_INTERPOLATION                  
        ENDIF
      ENDDO !xiIdx1
    ENDIF
    basis%collapsedXi(1:basis%numberOfXi)=collapsedXi(1:basis%numberOfXi)
     
    EXITS("Basis_CollapsedXiSet")
    RETURN
999 ERRORSEXITS("Basis_CollapsedXiSet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_CollapsedXiSet

  !
  !================================================================================================================================
  !

  !>Converts xi coordinates to area coordinates. \see BasisRoutiness::Basis_AreaToXiCoordinates
  SUBROUTINE Basis_XiToAreaCoordinates(xiCoordinates,areaCoordinates,err,error,*)
    
    !Argument variables
    REAL(DP), INTENT(IN) :: xiCoordinates(:) !<The xi coordinates to convert
    REAL(DP), INTENT(OUT) :: areaCoordinates(:) !<On return, the area coordinates corresponding to the xi coordinates
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Basis_XiToAreaCoordinates",err,error,*999)

    IF((SIZE(xiCoordinates,1)+1)/=SIZE(areaCoordinates,1)) THEN
      localError="Invalid number of coordinates. The number of xi coordinates of "// &
        & TRIM(NumberToVString(SIZE(xiCoordinates,1),"*",err,error))// &
        & " plus one should be equal to the number of area coordinates of "// &
        & TRIM(NumberToVString(SIZE(areaCoordinates,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    SELECT CASE(SIZE(xiCoordinates,1))
    CASE(1)
      areaCoordinates(1)=1.0_DP-xiCoordinates(1)
      areaCoordinates(2)=xiCoordinates(1)-1.0_DP
    CASE(2)
      areaCoordinates(1)=1.0_DP-xiCoordinates(1)
      areaCoordinates(2)=1.0_DP-xiCoordinates(2)
      areaCoordinates(3)=xiCoordinates(1)+xiCoordinates(2)-1.0_DP
    CASE(3)
      areaCoordinates(1)=1.0_DP-xiCoordinates(1)
      areaCoordinates(2)=1.0_DP-xiCoordinates(2)
      areaCoordinates(3)=1.0_DP-xiCoordinates(3)
      areaCoordinates(4)=xiCoordinates(1)+xiCoordinates(2)+xiCoordinates(3)-2.0_DP
    CASE DEFAULT
      localError="The number of xi coordinates of "//TRIM(NumberToVString(SIZE(xiCoordinates,1),"*",err,error))// &
        & " is invalid. The number must be >= 1 and <= 3."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Basis_XiToAreaCoordinates")
    RETURN
999 ERRORSEXITS("Basis_XiToAreaCoordinates",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_XiToAreaCoordinates

  !
  !================================================================================================================================
  !

  !>This routine calculates the weights and abscissae for a Gauss-Legendre quadrature scheme.
  !>\todo Fix this.
  SUBROUTINE Gauss_Legendre(n,alpha,beta,x,w,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: n !<The number of of Gauss points required.
    REAL(DP), INTENT(IN) :: alpha !<The lower limit of the integration scheme
    REAL(DP), INTENT(IN) :: beta !<The upper limit of the integration scheme
    REAL(DP), INTENT(OUT) :: x(:) !<x(gaussPointIdx). On exit the gaussPointIdx'th Gauss point location
    REAL(DP), INTENT(OUT) :: w(:) !<w(gaussPointIdx). On exit the gaussPointIdx'th Gauss point weight.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: gaussPointIdx
    REAL(DP) :: difference,t1,t2
    TYPE(VARYING_STRING) :: localError

    INTEGER(INTG) :: gaussStart(4) = [ 0,1,3,6 ]
    REAL(DP) :: gaussPointLocations(10),gaussPointWeights(10)
 
    gaussPointLocations = [ 0.500000000000000_DP, &
      &      (-1.0_DP/SQRT(3.0_DP)+1.0_DP)/2.0_DP,(+1.0_DP/SQRT(3.0_DP)+1.0_DP)/2.0_DP, &
      &      (-SQRT(0.6_DP)+1.0_DP)/2.0_DP, 0.5_DP, (+SQRT(0.6_DP)+1.0_DP)/2.0_DP, &
      &      (-SQRT((3.0_DP+2.0_DP*SQRT(6.0_DP/5.0_DP))/7.0_DP)+1.0_DP)/2.0_DP, &      
      &      (-SQRT((3.0_DP-2.0_DP*SQRT(6.0_DP/5.0_DP))/7.0_DP)+1.0_DP)/2.0_DP, &
      &      (+SQRT((3.0_DP-2.0_DP*SQRT(6.0_DP/5.0_DP))/7.0_DP)+1.0_DP)/2.0_DP, &
      &      (+SQRT((3.0_DP+2.0_DP*SQRT(6.0_DP/5.0_DP))/7.0_DP)+1.0_DP)/2.0_DP ]
    gaussPointWeights = [ 1.000000000000000_DP, &
      &      0.500000000000000_DP,0.500000000000000_DP, &
      &      2.5_DP/9.0_DP, 4.0_DP/9.0_DP, 2.5_DP/9.0_DP, &
      &      (18.0_DP-SQRT(30.0_DP))/72.0_DP, &
      &      (18.0_DP+SQRT(30.0_DP))/72.0_DP, &
      &      (18.0_DP+SQRT(30.0_DP))/72.0_DP, &
      &      (18.0_DP-SQRT(30.0_DP))/72.0_DP ]             
    
    ENTERS("Gauss_Legendre",err,error,*999)

    IF(n<1.OR.n>4) THEN
      localError="The specified number of Gauss points of "//TRIM(NumberToVString(n,"*",err,error))// &
        & " is invalid or not implemented."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(x,1)<n) THEN
      localError="The size of the x array of "//TRIM(NumberToVString(SIZE(x,1),"*",err,error))// &
        & " is too small for the specified number of Gauss points. The array needs to be of size "// &
        & TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(w,1)<n) THEN
      localError="The size of the w array of "//TRIM(NumberToVString(SIZE(w,1),"*",err,error))// &
        & " is too small for the specified number of Gauss points. The array needs to be of size "// &
        & TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO gaussPointIdx=1,n
      x(gaussPointIdx)=gaussPointLocations(gaussStart(n)+gaussPointIdx)
      w(gaussPointIdx)=gaussPointWeights(gaussStart(n)+gaussPointIdx)
    ENDDO !i
    
    IF(diagnostics1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of gauss points = ",n,err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,n,5,5,x,'("Gauss point locations :",5(X,F13.5))','(23X,5(X,F13.5))', &
        & err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,n,5,5,w,'("Gauss point weights   :",5(X,F13.5))','(23X,5(X,F13.5))', &
        & err,error,*999)
      IF(diagnostics2) THEN
        !Check by integrating y=x+1
        t1=0.0_DP !Numerical
        t2=0.0_DP !Analytic
        DO gaussPointIdx=1,n
          t1=t1+((x(gaussPointIdx)+1.0_DP)*w(gaussPointIdx))
        ENDDO !gaussPointIdx
        t2=(beta**2.0_DP/2.0_DP+beta)-(alpha**2.0_DP/2.0_DP-alpha)
        difference=ABS(t2-t1)
        CALL WriteStringFmtValue(DIAGNOSTIC_OUTPUT_TYPE,"Numerical Integration Test Difference = ",difference,'(F14.6)', &
          & err,error,*999)
      ENDIF
    ENDIF

    EXITS("Gauss_Legendre")
    RETURN
999 ERRORSEXITS("Gauss_Legendre",err,error)
    RETURN 1
    
  END SUBROUTINE Gauss_Legendre
  
  !
  !================================================================================================================================
  !

  !>This routine calculates the weights and abscissae for a Gauss quadrature scheme for simplex elements.
  !>
  !>Reference: Liu, Yen and Vinokur, Marcel. "Exact Integrations of Polynomials and Symmetric Quadrature Formulas
  !> over Arbitrary Polyhedral Grids", Journal of Computational Physics, 140:122-147 (1998).
  !>
  SUBROUTINE Gauss_Simplex(order,numberOfVertices,n,x,w,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: order !<The desired order of the scheme. Currently, the maximum order is 5.
    INTEGER(INTG), INTENT(IN) :: numberOfVertices !<The number of vertices. 2, 3 or 4 for lines, triangles or tetrahedra.
    INTEGER(INTG), INTENT(OUT) :: n !<On exit, the number of Gauss points 
    REAL(DP), INTENT(OUT) :: x(:,:) !<X(coordinateIdx,gaussPointIdx). On exit, the returned positions in area coordinates.
    REAL(DP), INTENT(OUT) :: w(:) !<W(gaussPointIdx). On exit, the returned weights.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: gaussIdx
    REAL(DP) :: alpha1,alpha2,beta,lambda,lc,l1Alpha1,l2Alpha1,l3Alpha1,l4Alpha1,l1Alpha2,l2Alpha2,l3Alpha2,l4Alpha2,l1Beta, &
      & l2Beta,l3Beta,l4Beta,wc,wAlpha1,wAlpha2,wBeta,aCosArg
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Gauss_Simplex",err,error,*999)
    
    IF(SIZE(x,1)<numberOfVertices) THEN
      localError="The first dimension of the X array is "//TRIM(NumberToVString(SIZE(x,1),"*",err,error))// &
        & " and it must be >= the number of vertices of "//TRIM(NumberToVString(numberOfVertices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    SELECT CASE(numberOfVertices)
    CASE(2)
      !Line
      SELECT CASE(order)
      CASE(1)
        n=1
        IF(SIZE(x,2)<n) THEN
          localError="The second dimension of the X array is "//TRIM(NumberToVString(SIZE(x,2),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(SIZE(w,1)<n) THEN
          localError="The first dimension of the W array is "//TRIM(NumberToVString(SIZE(w,1),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        !Gauss point 1
        x(1,1)=1.0_DP/2.0_DP
        x(2,1)=1.0_DP/2.0_DP
        w(1)=1.0_DP
      CASE(2)
        n=2
        IF(SIZE(x,2)<n) THEN
          localError="The second dimension of the X array is "//TRIM(NumberToVString(SIZE(x,2),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(SIZE(w,1)<n) THEN
          localError="The first dimension of the W array is "//TRIM(NumberToVString(SIZE(w,1),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        !Gauss point 1
        x(1,1)=(1.0_DP-1.0_DP/sqrt(3.0_DP))/2.0_DP
        x(2,1)=1.0_DP-x(1,1)
        w(1)=1.0_DP/2.0_DP
        !Gauss point 2
        x(1,2)=(1.0_DP+1.0_DP/SQRT(3.0_DP))/2.0_DP
        x(2,2)=1.0_DP-x(1,2)
        w(2)=1.0_DP/2.0_DP
      CASE(3)
        n=2
        IF(SIZE(x,2)<n) THEN
          localError="The second dimension of the X array is "//TRIM(NumberToVString(SIZE(x,2),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(SIZE(w,1)<n) THEN
          localError="The first dimension of the W array is "//TRIM(NumberToVString(SIZE(w,1),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        !Gauss point 1
        x(1,1)=(1.0_DP-1.0_DP/sqrt(3.0_DP))/2.0_DP
        x(2,1)=1.0_DP-x(1,1)
        w(1)=1.0_DP/2.0_DP
        !Gauss point 2
        x(1,2)=(1.0_DP+1.0_DP/SQRT(3.0_DP))/2.0_DP
        x(2,2)=1.0_DP-x(1,2)
        w(2)=1.0_DP/2.0_DP
      CASE(4)
        n=3
        IF(SIZE(x,2)<n) THEN
          localError="The second dimension of the X array is "//TRIM(NumberToVString(SIZE(x,2),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(SIZE(w,1)<n) THEN
          localError="The first dimension of the W array is "//TRIM(NumberToVString(SIZE(w,1),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        !Gauss point 1
        x(1,1)=(1.0_DP-SQRT(0.6_DP))/2.0_DP
        x(2,1)=1-x(1,1)
        w(1)=5.0_DP/18.0_DP
        !Gauss point 2
        x(1,2)=1.0_DP/2.0_DP
        x(2,2)=1.0_DP-x(1,2)
        w(2)=4.0_DP/9.0_DP
        !Gauss point 3
        x(1,3)=(1.0_DP+SQRT(0.6_DP))/2.0_DP
        x(2,3)=1.0_DP-x(1,3)
        w(3)=5.0_DP/18.0_DP
      CASE(5)
        n=3
        IF(SIZE(x,2)<n) THEN
          localError="The second dimension of the X array is "//TRIM(NumberToVString(SIZE(x,2),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(SIZE(w,1)<n) THEN
          localError="The first dimension of the W array is "//TRIM(NumberToVString(SIZE(w,1),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        !Gauss point 1
        x(1,1)=(1.0_DP-SQRT(0.6_DP))/2.0_DP
        x(2,1)=1.0_DP-x(1,1)
        w(1)=5.0_DP/18.0_DP
        !Gauss point 2
        x(1,2)=1.0_DP/2.0_DP
        x(2,2)=1.0_DP-x(1,2)
        w(2)=4.0_DP/9.0_DP
        !Gauss point 3
        x(1,3)=(1.0_DP+SQRT(0.6_DP))/2.0_DP
        x(2,3)=1.0_DP-x(1,3)
        w(3)=5.0_DP/18.0_DP
      CASE DEFAULT
        localError="The specified Gauss order of "//TRIM(NumberToVString(order,"*",err,error))// &
          & " is an invalid. You must specify an order between 1 and 5."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(3)
      !Triangle
      SELECT CASE(order)
      CASE(1)
        n=1
        IF(SIZE(x,2)<n) THEN
          localError="The second dimension of the x array is "//TRIM(NumberToVString(SIZE(x,2),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(SIZE(w,1)<n) THEN
          localError="The first dimension of the w array is "//TRIM(NumberToVString(SIZE(w,1),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        lC=1.0_DP/3.0_DP
        wC=1.0_DP
        !Gauss point 1
        x(1,1)=lC
        x(2,1)=lC
        x(3,1)=lC
        w(1)=wC/2.0_DP
      CASE(2)
        n=3
        IF(SIZE(x,2)<n) THEN
          localError="The second dimension of the x array is "//TRIM(NumberToVString(SIZE(x,2),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(SIZE(w,1)<n) THEN
          localError="The first dimension of the w array is "//TRIM(NumberToVString(SIZE(w,1),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        alpha1=-1.0_DP/2.0_DP
        wAlpha1=1.0_DP/3.0_DP
        l1Alpha1=(1.0_DP+2.0_DP*alpha1)/3.0_DP
        l2Alpha1=(1.0_DP-alpha1)/3.0_DP
        l3Alpha1=1.0_DP-l1Alpha1-l2Alpha1
        !Gauss point 1
        x(1,1)=l1Alpha1
        x(2,1)=l2Alpha1
        x(3,1)=l3Alpha1
        w(1)=wAlpha1/2.0_DP
        !Gauss point 2
        x(1,2)=l3Alpha1
        x(2,2)=l1Alpha1
        x(3,2)=l2Alpha1
        w(2)=wAlpha1/2.0_DP
        !Gauss point 3
        x(1,3)=l2Alpha1
        x(2,3)=l3Alpha1
        x(3,3)=l1Alpha1
        w(3)=wAlpha1/2.0_DP
      CASE(3)
        n=4
        IF(SIZE(x,2)<n) THEN
          localError="The second dimension of the x array is "//TRIM(NumberToVString(SIZE(x,2),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(SIZE(w,1)<n) THEN
          localError="The first dimension of the w array is "//TRIM(NumberToVString(SIZE(w,1),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        lC=1.0_DP/3.0_DP
        wC=-3.0_DP/4.0_DP
        alpha1=2.0_DP/5.0_DP
        wAlpha1=25.0_DP/48.0_DP
        l1Alpha1=(1.0_DP+2.0_DP*alpha1)/3.0_DP
        l2Alpha1=(1.0_DP-alpha1)/3.0_DP
        l3Alpha1=1.0_DP-l1Alpha1-l2Alpha1
        !Gauss point 1
        x(1,1)=lC
        x(2,1)=lC
        x(3,1)=lC
        w(1)=wC/2.0_DP
        !Gauss point 2
        x(1,2)=l1Alpha1
        x(2,2)=l2Alpha1
        x(3,2)=l3Alpha1
        w(2)=wAlpha1/2.0_DP
        !Gauss point 3
        x(1,3)=l3Alpha1
        x(2,3)=l1Alpha1
        x(3,3)=l2Alpha1
        w(3)=wAlpha1/2.0_DP
        !Gauss point 4
        x(1,4)=l2Alpha1
        x(2,4)=l3Alpha1
        x(3,4)=l1Alpha1
        w(4)=wAlpha1/2.0_DP
      CASE(4)
        n=6
        IF(SIZE(x,2)<n) THEN
          localError="The second dimension of the x array is "//TRIM(NumberToVString(SIZE(x,2),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(SIZE(w,1)<n) THEN
          localError="The first dimension of the w array is "//TRIM(NumberToVString(SIZE(w,1),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        alpha1=(-10.0_DP+5.0_DP*SQRT(10.0_DP)+SQRT(950.0_DP-220.0_DP*SQRT(10.0_DP)))/30.0_DP
        alpha2=(-10.0_DP+5.0_DP*SQRT(10.0_DP)-SQRT(950.0_DP-220.0_DP*SQRT(10.0_DP)))/30.0_DP
        wAlpha1=(5.0_DP*alpha2-2.0_DP)/(60.0_DP*alpha1*alpha1*(alpha2-alpha1))
        wAlpha2=(5.0_DP*alpha1-2.0_DP)/(60.0_DP*alpha2*alpha2*(alpha1-alpha2))
        l1Alpha1=(1.0_DP+2.0_DP*alpha1)/3.0_DP
        l2Alpha1=(1.0_DP-alpha1)/3.0_DP
        l3Alpha1=1.0_DP-l1Alpha1-l2Alpha1
        l1Alpha2=(1.0_DP+2.0_DP*alpha2)/3.0_DP
        l2Alpha2=(1.0_DP-alpha2)/3.0_DP
        l3Alpha2=1.0_DP-l1Alpha2-l2Alpha2
        !Gauss point 1
        x(1,1)=l1Alpha1
        x(2,1)=l2Alpha1
        x(3,1)=l3Alpha1
        w(1)=wAlpha1/2.0_DP 
        !Gauss point 2
        x(1,2)=l3Alpha1
        x(2,2)=l1Alpha1
        x(3,2)=l2Alpha1
        w(2)=wAlpha1/2.0_DP
        !Gauss point 3
        x(1,3)=l2Alpha1
        x(2,3)=l3Alpha1
        x(3,3)=l1Alpha1
        w(3)=wAlpha1/2.0_DP
        !Gauss point 4
        x(1,4)=l1Alpha2
        x(2,4)=l2Alpha2
        x(3,4)=l3Alpha2
        w(4)=wAlpha2/2.0_DP
        !Gauss point 5
        x(1,5)=l3Alpha2
        x(2,5)=l1Alpha2
        x(3,5)=l2Alpha2
        w(5)=wAlpha2/2.0_DP
        !Gauss point 6
        x(1,6)=l2Alpha2
        x(2,6)=l3Alpha2
        x(3,6)=l1Alpha2
        w(6)=wAlpha2/2.0_DP
      CASE(5)
        n=7
        IF(SIZE(x,2)<n) THEN
          localError="The second dimension of the x array is "//TRIM(NumberToVString(SIZE(x,2),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(SIZE(w,1)<n) THEN
          localError="The first dimension of the w array is "//TRIM(NumberToVString(SIZE(w,1),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        lC=1.0_DP/3.0_DP
        wC=9.0_DP/40.0_DP
        alpha1=(1.0_DP+SQRT(15.0_DP))/7.0_DP
        alpha2=(1.0_DP-SQRT(15.0_DP))/7.0_DP
        wAlpha1=(155.0_DP-SQRT(15.0_DP))/1200.0_DP
        wAlpha2=(155.0_DP+SQRT(15.0_DP))/1200.0_DP
        l1Alpha1=(1.0_DP+2.0_DP*alpha1)/3.0_DP
        l2Alpha1=(1.0_DP-alpha1)/3.0_DP
        l3Alpha1=1.0_DP-l1Alpha1-l2Alpha1
        l1Alpha2=(1.0_DP+2.0_DP*alpha2)/3.0_DP
        l2Alpha2=(1.0_DP-alpha2)/3.0_DP
        l3Alpha2=1.0_DP-l1Alpha2-l2Alpha2
        !Gauss point 1
        x(1,1)=lC
        x(2,1)=lC
        x(3,1)=lC
        w(1)=wC/2.0_DP
        !Gauss point 2
        x(1,2)=l1Alpha1
        x(2,2)=l2Alpha1
        x(3,2)=l3Alpha1
        w(2)=wAlpha1/2.0_DP
        !Gauss point 3
        x(1,3)=l3Alpha1
        x(2,3)=l1Alpha1
        x(3,3)=l2Alpha1
        w(3)=wAlpha1/2.0_DP
        !Gauss point 4
        x(1,4)=l2Alpha1
        x(2,4)=l3Alpha1
        x(3,4)=l1Alpha1
        w(4)=wAlpha1/2.0_DP
        !Gauss point 5
        x(1,5)=l1Alpha2
        x(2,5)=l2Alpha2
        x(3,5)=l3Alpha2
        w(5)=wAlpha2/2.0_DP
        !Gauss point 6
        x(1,6)=l3Alpha2
        x(2,6)=l1Alpha2
        x(3,6)=l2Alpha2
        w(6)=wAlpha2/2.0_DP
        !Gauss point 7
        x(1,7)=l2Alpha2
        x(2,7)=l3Alpha2
        x(3,7)=l1Alpha2
        w(7)=wAlpha2/2.0_DP
      CASE DEFAULT
        localError="The specified order of "//TRIM(NumberToVString(order,"*",err,error))// &
          & " is an invalid. You must have an order between 1 and 5."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(4)
      !Tetrahedra
      SELECT CASE(ORDER)
      CASE(1)
        n=1
        IF(SIZE(x,2)<n) THEN
          localError="The second dimension of the x array is "//TRIM(NumberToVString(SIZE(x,2),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(SIZE(w,1)<n) THEN
          localError="The first dimension of the w array is "//TRIM(NumberToVString(SIZE(w,1),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        lC=1.0_DP/4.0_DP
        wC=1.0_DP
        !Gauss point 1
        x(1,1)=lC
        x(2,1)=lC
        x(3,1)=lC
        w(1)=wC/6.0_DP 
      CASE(2)
        n=4
        IF(SIZE(x,2)<n) THEN
          localError="The second dimension of the x array is "//TRIM(NumberToVString(SIZE(x,2),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(SIZE(w,1)<n) THEN
          localError="The first dimension of the w array is "//TRIM(NumberToVString(SIZE(w,1),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        alpha1=1.0_DP/SQRT(5.0_DP)
        wAlpha1=1.0_DP/4.0_DP
        l1Alpha1=(1.0_DP+3.0_DP*alpha1)/4.0_DP
        l2Alpha1=(1.0_DP-alpha1)/4.0_DP
        l3Alpha1=(1.0_DP-alpha1)/4.0_DP
        l4Alpha1=1.0_DP-l1Alpha1-l2Alpha1-l3Alpha1
        !Gauss point 1
        x(1,1)=l1Alpha1
        x(2,1)=l2Alpha1
        x(3,1)=l3Alpha1
        x(4,1)=l4Alpha1
        w(1)=wAlpha1/6.0_DP 
        !Gauss point 2
        x(1,2)=l4Alpha1
        x(2,2)=l1Alpha1
        x(3,2)=l2Alpha1
        x(4,2)=l3Alpha1
        w(2)=wAlpha1/6.0_DP
        !Gauss point 3
        x(1,3)=l3Alpha1
        x(2,3)=l4Alpha1
        x(3,3)=l1Alpha1
        x(4,3)=l2Alpha1
        w(3)=wAlpha1/6.0_DP
        !Gauss point 4
        x(1,4)=l2Alpha1
        x(2,4)=l3Alpha1
        x(3,4)=l4Alpha1
        x(4,4)=l1Alpha1
        w(4)=wAlpha1/6.0_DP
      CASE(3)
        n=5
        IF(SIZE(x,2)<n) THEN
          localError="The second dimension of the x array is "//TRIM(NumberToVString(SIZE(x,2),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(SIZE(w,1)<n) THEN
          localError="The first dimension of the w array is "//TRIM(NumberToVString(SIZE(w,1),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        lC=1.0_DP/4.0_DP
        wC=-4.0_DP/5.0_DP
        alpha1=1.0_DP/3.0_DP
        wAlpha1=9.0_DP/20.0_DP
        l1Alpha1=(1.0_DP+3.0_DP*alpha1)/4.0_DP
        l2Alpha1=(1.0_DP-alpha1)/4.0_DP
        l3Alpha1=(1.0_DP-alpha1)/4.0_DP
        l4Alpha1=1.0_DP-l1Alpha1-l2Alpha1-l3Alpha1
        !Gauss point 1
        x(1,1)=lC
        x(2,1)=lC
        x(3,1)=lC
        x(4,1)=lC
        w(1)=wC/6.0_DP 
        !Gauss point 2
        x(1,2)=l1Alpha1
        x(2,2)=l2Alpha1
        x(3,2)=l3Alpha1
        x(4,2)=l4Alpha1
        w(2)=wAlpha1/6.0_DP
        !Gauss point 3
        x(1,3)=l4Alpha1
        x(2,3)=l1Alpha1
        x(3,3)=l2Alpha1
        x(4,3)=l3Alpha1
        w(3)=wAlpha1/6.0_DP
        !Gauss point 4
        x(1,4)=l3Alpha1
        x(2,4)=l4Alpha1
        x(3,4)=l1Alpha1
        x(4,4)=l2Alpha1
        w(4)=wAlpha1/6.0_DP
        !Gauss point 5
        x(1,5)=l2Alpha1
        x(2,5)=l3Alpha1
        x(3,5)=l4Alpha1
        x(4,5)=l1Alpha1
        w(5)=wAlpha1/6.0_DP
      CASE(4)
        n=11
        IF(SIZE(x,2)<n) THEN
          localError="The second dimension of the x array is "//TRIM(NumberToVString(SIZE(x,2),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(SIZE(w,1)<n) THEN
          localError="The first dimension of the w array is "//TRIM(NumberToVString(SIZE(w,1),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        lC=1.0_DP/4.0_DP
        wC=-148.0_DP/1875.0_DP
        alpha1=5.0_DP/7.0_DP
        beta=SQRT(70.0_DP)/28.0_DP
        wAlpha1=343.0_DP/7500.0_DP
        wBeta=56.0_DP/375.0_DP
        l1Alpha1=(1.0_DP+3.0_DP*alpha1)/4.0_DP
        l2Alpha1=(1.0_DP-alpha1)/4.0_DP
        l3Alpha1=(1.0_DP-alpha1)/4.0_DP
        l4Alpha1=1.0_DP-l1Alpha1-l2Alpha1-l3Alpha1
        l1Beta=(1.0_DP+2.0_DP*beta)/4.0_DP
        l2Beta=l1Beta
        l3Beta=(1.0_DP-2.0_DP*beta)/4.0_DP
        l4Beta=1.0_DP-l1Beta-l2Beta-l3Beta
        !Gauss point 1
        x(1,1)=lC
        x(2,1)=lC
        x(3,1)=lC
        x(4,1)=lC
        w(1)=wC/6.0_DP
        !Gauss point 2
        x(1,2)=l1Alpha1
        x(2,2)=l2Alpha1
        x(3,2)=l3Alpha1
        x(4,2)=l4Alpha1
        w(2)=wAlpha1/6.0_DP
        !Gauss point 3
        x(1,3)=l4Alpha1
        x(2,3)=l1Alpha1
        x(3,3)=l2Alpha1
        x(4,3)=l3Alpha1
        w(3)=wAlpha1/6.0_DP
        !Gauss point 4
        x(1,4)=l3Alpha1
        x(2,4)=l4Alpha1
        x(3,4)=l1Alpha1
        x(4,4)=l2Alpha1        
        w(4)=wAlpha1/6.0_DP
        !Gauss point 5
        x(1,5)=l2Alpha1
        x(2,5)=l3Alpha1
        x(3,5)=l4Alpha1
        x(4,5)=l1Alpha1
        w(5)=wAlpha1/6.0_DP
        !Gauss point 6
        x(1,6)=l1Beta
        x(2,6)=l2Beta
        x(3,6)=l3Beta
        x(4,6)=l4Beta
        w(6)=wBeta/6.0_DP
        !Gauss point 7
        x(1,7)=l1Beta
        x(2,7)=l3Beta
        x(3,7)=l2Beta
        x(4,7)=l4Beta
        w(7)=wBeta/6.0_DP
        !Gauss point 8
        x(1,8)=l1Beta
        x(2,8)=l3Beta
        x(3,8)=l4Beta
        x(4,8)=l2Beta
        w(8)=wBeta/6.0_DP
        !Gauss point 9
        x(1,9)=l3Beta
        x(2,9)=l1Beta
        x(3,9)=l2Beta
        x(4,9)=l4Beta
        w(9)=wBeta/6.0_DP
        !Gauss point 10
        x(1,10)=l3Beta
        x(2,10)=l1Beta
        x(3,10)=l4Beta
        x(4,10)=l2Beta
        w(10)=wBeta/6.0_DP
        !Gauss point 11
        x(1,11)=l3Beta
        x(2,11)=l4Beta
        x(3,11)=l1Beta
        x(4,11)=l2Beta
        w(11)=wBeta/6.0_DP
      CASE(5)
        n=14
        IF(SIZE(x,2)<n) THEN
          localError="The second dimension of the x array is "//TRIM(NumberToVString(SIZE(x,2),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(SIZE(w,1)<n) THEN
          localError="The first dimension of the w array is "//TRIM(NumberToVString(SIZE(w,1),"*",err,error))// &
            & " and it must be >= "//TRIM(NumberToVString(n,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        aCosArg=67.0_DP*SQRT(79.0_DP)/24964.0_DP
        lambda=4.0_DP/27.0_DP*(4.0_DP*SQRT(79.0_DP)*COS(((ACOS(aCosArg)+TWOPI)/3.0_DP))+71.0_DP)        
        alpha1=(SQRT(9.0_DP*lambda*lambda-248.0_DP*lambda+1680.0_DP)+28.0_DP-3.0_DP*lambda)/ &
          & (112.0_DP-10.0_DP*lambda)
        alpha2=(-1.0_DP*SQRT(9.0_DP*lambda*lambda-248.0_DP*lambda+1680.0_DP)+28.0_DP-3.0_DP*lambda)/ &
          & (112.0_DP-10.0_DP*lambda)
        beta=1.0_DP/SQRT(lambda)
        wAlpha1=((21.0_DP-lambda)*alpha2-7.0_DP)/(420.0_DP*alpha1*alpha1*(alpha2-alpha1))
        wAlpha2=((21.0_DP-lambda)*alpha1-7.0_DP)/(420.0_DP*alpha2*alpha2*(alpha1-alpha2))
        wBeta=lambda*lambda/840.0_DP
        l1Alpha1=(1.0_DP+3.0_DP*alpha1)/4.0_DP
        l2Alpha1=(1.0_DP-alpha1)/4.0_DP
        l3Alpha1=(1.0_DP-alpha1)/4.0_DP
        l4Alpha1=1.0_DP-l1Alpha1-l2Alpha1-l3Alpha1
        l1Alpha2=(1.0_DP+3.0_DP*alpha2)/4.0_DP
        l2Alpha2=(1.0_DP-alpha2)/4.0_DP
        l3Alpha2=(1.0_DP-alpha2)/4.0_DP
        l4Alpha2=1.0_DP-l1Alpha2-l2Alpha2-l3Alpha2
        l1Beta=(1.0_DP+2.0_DP*beta)/4.0_DP
        l2Beta=l1Beta
        l3Beta=(1.0_DP-2.0_DP*beta)/4.0_DP
        l4Beta=1.0_DP-l1Beta-l2Beta-l3Beta
        !Gauss point 1
        x(1,1)=l1Alpha1
        x(2,1)=l2Alpha1
        x(3,1)=l3Alpha1
        x(4,1)=l4Alpha1
        w(1)=wAlpha1/6.0_DP
        !Gauss point 2
        x(1,2)=l4Alpha1
        x(2,2)=l1Alpha1
        x(3,2)=l2Alpha1
        x(4,2)=l3Alpha1
        w(2)=wAlpha1/6.0_DP
        !Gauss point 3
        x(1,3)=l3Alpha1
        x(2,3)=l4Alpha1
        x(3,3)=l1Alpha1
        x(4,3)=l2Alpha1
        w(3)=wAlpha1/6.0_DP
        !Gauss point 4
        x(1,4)=l2Alpha1
        x(2,4)=l3Alpha1
        x(3,4)=l4Alpha1
        x(4,4)=l1Alpha1
        w(4)=wAlpha1/6.0_DP
        !Gauss point 5
        x(1,5)=l1Alpha2
        x(2,5)=l2Alpha2
        x(3,5)=l3Alpha2
        x(4,5)=l4Alpha2
        w(5)=wAlpha2/6.0_DP
        !Gauss point 6
        x(1,6)=l4Alpha2
        x(2,6)=l1Alpha2
        x(3,6)=l2Alpha2
        x(4,6)=l3Alpha2
        w(6)=wAlpha2/6.0_DP
        !Gauss point 7
        x(1,7)=l3Alpha2
        x(2,7)=l4Alpha2
        x(3,7)=l1Alpha2
        x(4,7)=l2Alpha2
        w(7)=wAlpha2/6.0_DP
        !Gauss point 8
        x(1,8)=l2Alpha2
        x(2,8)=l3Alpha2
        x(3,8)=l4Alpha2
        x(4,8)=l1Alpha2
        w(8)=wAlpha2/6.0_DP
        !Gauss point 9
        x(1,9)=l1Beta
        x(2,9)=l2Beta
        x(3,9)=l3Beta
        x(4,9)=l4Beta
        w(9)=wBeta/6.0_DP
        !Gauss point 10
        x(1,10)=l1Beta
        x(2,10)=l3Beta
        x(3,10)=l2Beta
        x(4,10)=l4Beta
        w(10)=wBeta/6.0_DP
        !Gauss point 11
        x(1,11)=l1Beta
        x(2,11)=l3Beta
        x(3,11)=l4Beta
        x(4,11)=l2Beta
        w(11)=wBeta/6.0_DP
        !Gauss point 12
        x(1,12)=l3Beta
        x(2,12)=l1Beta
        x(3,12)=l2Beta
        x(4,12)=l4Beta
        w(12)=wBeta/6.0_DP
        !Gauss point 13
        x(1,13)=l3Beta
        x(2,13)=l1Beta
        x(3,13)=l4Beta
        x(4,13)=l2Beta
        w(13)=wBeta/6.0_DP
        !Gauss point 14
        x(1,14)=l3Beta
        x(2,14)=l4Beta
        x(3,14)=l1Beta
        x(4,14)=l2Beta
        w(14)=wBeta/6.0_DP
      CASE DEFAULT
        localError="The specified order of "//TRIM(NumberToVString(order,"*",err,error))// &
          & " is an invalid. You must have an order between 1 and 5."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The specified number of vertices of "//TRIM(NumberToVString(numberOfVertices,"*",err,error))// &
        & " is an invalid. You must have between 2 and 4 vertices."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Simplex Gauss quadrature points:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of vertices = ",numberOfVertices,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Order = ",order,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of gauss points = ",n,err,error,*999)
      DO gaussIdx=1,n
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Gauss point ",gaussIdx,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,order,4,4,x(:,gaussIdx),'("        Location :",4(X,F13.5))', &
          & '(18X,4(X,F13.5))',err,error,*999)
        CALL WriteStringFmtValue(DIAGNOSTIC_OUTPUT_TYPE,"        Weight   : ",w(gaussIdx),"F13.5",err,error,*999)
      ENDDO !ng
      IF(diagnostics2) THEN
!!TODO: \todo add in integral check
      ENDIF
    ENDIF

    EXITS("Gauss_Simplex")
    RETURN
999 ERRORSEXITS("Gauss_Simplex",err,error)
    RETURN 1
    
  END SUBROUTINE Gauss_Simplex
  
  !
  !================================================================================================================================
  !

  !>Evaluates a 1D cubic Hermite basis function.
  FUNCTION Hermite_CubicEvaluateDP(localNodeIndex,nodeDerivativeIndex,partialDerivativeIndex,xi,err,error)
  
   !Argument variables
    INTEGER(INTG), INTENT(IN) :: localNodeIndex !<The local node number of the basis. Must be between 1 and 2.
    INTEGER(INTG), INTENT(IN) :: nodeDerivativeIndex !<The local derivative number of the basis. Must be between 1 and 2.
    INTEGER(INTG), INTENT(IN) :: partialDerivativeIndex !<The partial derivative index to evaluate.
    REAL(DP), INTENT(IN) :: xi !<xi(xiIdx). The Xi location to evaluate the basis at.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: Hermite_CubicEvaluateDP !<On exit the evaluated basis function.
    !Local variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Hermite_CubicEvaluateDP",err,error,*999)

    Hermite_CubicEvaluateDP=0.0_DP
    
    SELECT CASE(partialDerivativeIndex)
    CASE(NO_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        SELECT CASE(nodeDerivativeIndex)
        CASE(1)
          Hermite_CubicEvaluateDP=(2.0_DP*xi-3.0_DP)*xi*xi+1.0_DP ! 2xi^3-3xi^2+1
        CASE(2)
          Hermite_CubicEvaluateDP=((xi-2.0_DP)*xi+1.0_DP)*xi ! xi^3-2xi^2+xi
        CASE DEFAULT
          localError="The specified node derivative index of "//TRIM(NumberToVString(nodeDerivativeIndex,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)          
        END SELECT
      CASE(2)
        SELECT CASE(nodeDerivativeIndex)
        CASE(1)
          Hermite_CubicEvaluateDP=xi*xi*(3.0_DP-2.0_DP*xi) ! -2xi^3+3xi^2
        CASE(2)
          Hermite_CubicEvaluateDP=xi*xi*(xi-1.0_DP) ! xi^3-xi^2
        CASE DEFAULT
          localError="The specified node derivative index of "//TRIM(NumberToVString(nodeDerivativeIndex,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)          
        END SELECT
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE(FIRST_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        SELECT CASE(nodeDerivativeIndex)
        CASE(1)
          Hermite_CubicEvaluateDP=6.0_DP*xi*(xi-1.0_DP) ! 6xi^2-6xi
        CASE(2)
          Hermite_CubicEvaluateDP=(3.0_DP*xi-4.0_DP)*xi+1.0_DP ! 3xi^2-4xi+1
        CASE DEFAULT
          localError="The specified node derivative index of "//TRIM(NumberToVString(nodeDerivativeIndex,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)          
        END SELECT
      CASE(2)
        SELECT CASE(nodeDerivativeIndex)
        CASE(1)
          Hermite_CubicEvaluateDP=6.0_DP*xi*(1.0_DP-xi) ! -6xi^2+6xi
        CASE(2)
          Hermite_CubicEvaluateDP=xi*(3.0_DP*xi-2.0_DP) ! 3xi^2-2xi
        CASE DEFAULT
          localError="The specified node derivative index of "//TRIM(NumberToVString(nodeDerivativeIndex,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)          
        END SELECT
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE(SECOND_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        SELECT CASE(nodeDerivativeIndex)
        CASE(1)
          Hermite_CubicEvaluateDP=12.0_DP*xi-6.0_DP ! 12xi-6
        CASE(2)
          Hermite_CubicEvaluateDP=6.0_DP*xi-4.0_DP ! 6xi-4
        CASE DEFAULT
          localError="The specified node derivative index of "//TRIM(NumberToVString(nodeDerivativeIndex,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)          
        END SELECT
      CASE(2)
        SELECT CASE(nodeDerivativeIndex)
        CASE(1)
          Hermite_CubicEvaluateDP=6.0_DP-12.0_DP*xi ! -12xi+6
        CASE(2)
          Hermite_CubicEvaluateDP=6.0_DP*xi-2.0_DP ! 6xi-2
        CASE DEFAULT
          localError="The specified node derivative index of "//TRIM(NumberToVString(nodeDerivativeIndex,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)          
        END SELECT
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE DEFAULT
      localError="The specified partial derivative index of "//TRIM(NumberToVString(partialDerivativeIndex,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)          
    END SELECT

    EXITS("Hermite_CubicEvaluateDP")
    RETURN
999 ERRORSEXITS("Hermite_CubicEvaluateDP",err,error)
    RETURN
    
  END FUNCTION Hermite_CubicEvaluateDP
 
  !
  !================================================================================================================================
  !

  !>Evaluates a 1D quadratic Hermite basis function at position xi, and with the give localNodeIndex, nodeDerivativeIndex and partialDerivativeIndex. specialNodeIndex is the node with no derivative term.
  FUNCTION Hermite_QuadraticEvaluateDP(localNodeIndex,nodeDerivativeIndex,partialDerivativeIndex,specialNodeIndex,xi,err,error)
      
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: localNodeIndex !<The local node number of the basis. Must be between 1 and 2.
    INTEGER(INTG), INTENT(IN) :: nodeDerivativeIndex !<The local derivative number of the basis. Must be between 1 and 2.
    INTEGER(INTG), INTENT(IN) :: partialDerivativeIndex !<The partial derivative index to evaluate.
    INTEGER(INTG), INTENT(IN) :: specialNodeIndex !<The local node number with no derivative term.
    REAL(DP), INTENT(IN) :: xi !<The Xi location to evaluate the basis at.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: Hermite_QuadraticEvaluateDP
    !Local variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Hermite_QuadraticEvaluateDP",err,error,*999)
    
    Hermite_QuadraticEvaluateDP=0.0_DP
    
    SELECT CASE(specialNodeIndex)
    CASE(1)
      SELECT CASE(partialDerivativeIndex)
      CASE(NO_PART_DERIV)
        SELECT CASE(localNodeIndex)
        CASE(1)
          SELECT CASE(nodeDerivativeIndex)
          CASE(1)
            Hermite_QuadraticEvaluateDP=(xi-2.0_DP)*xi+1.0_DP ! xi^2-2xi+1
          CASE(2)
            Hermite_QuadraticEvaluateDP=0.0_DP ! 0
          CASE DEFAULT
            localError="The specified node derivative index of "//TRIM(NumberToVString(nodeDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)          
          END SELECT
        CASE(2)
          SELECT CASE(nodeDerivativeIndex)
          CASE(1)
            Hermite_QuadraticEvaluateDP=(2.0_DP-xi)*xi ! -xi^2+2xi
          CASE(2)
            Hermite_QuadraticEvaluateDP=(xi-1.0_DP)*xi ! xi^2-xi
          CASE DEFAULT
            localError="The specified node derivative index of "//TRIM(NumberToVString(nodeDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)          
          END SELECT
        CASE DEFAULT
          localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)          
        END SELECT
      CASE(FIRST_PART_DERIV)
        SELECT CASE(localNodeIndex)
        CASE(1)
          SELECT CASE(nodeDerivativeIndex)
          CASE(1)
            Hermite_QuadraticEvaluateDP=2.0_DP*xi-2.0_DP  ! 2xi-2
          CASE(2)
            Hermite_QuadraticEvaluateDP=0.0_DP ! 0
          CASE DEFAULT
            localError="The specified node derivative index of "//TRIM(NumberToVString(nodeDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)          
          END SELECT
        CASE(2)
          SELECT CASE(nodeDerivativeIndex)
          CASE(1)
            Hermite_QuadraticEvaluateDP=-2.0_DP*xi+2.0_DP ! -2xi+2
          CASE(2)
            Hermite_QuadraticEvaluateDP=2.0_DP*xi-1.0_DP ! 2xi-1
          CASE DEFAULT
            localError="The specified node derivative index of "//TRIM(NumberToVString(nodeDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)          
          END SELECT
        CASE DEFAULT
          localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)          
        END SELECT
      CASE(SECOND_PART_DERIV)
        SELECT CASE(localNodeIndex)
        CASE(1)
          SELECT CASE(nodeDerivativeIndex)
          CASE(1)
            Hermite_QuadraticEvaluateDP=2.0_DP ! 2
          CASE(2)
            Hermite_QuadraticEvaluateDP=0.0_DP ! 0
          CASE DEFAULT
            localError="The specified node derivative index of "//TRIM(NumberToVString(nodeDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)          
          END SELECT
        CASE(2)
          SELECT CASE(nodeDerivativeIndex)
          CASE(1)
            Hermite_QuadraticEvaluateDP=-2.0_DP ! -2
          CASE(2)
            Hermite_QuadraticEvaluateDP=2.0_DP ! 2
          CASE DEFAULT
            localError="The specified node derivative index of "//TRIM(NumberToVString(nodeDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)          
          END SELECT
        CASE DEFAULT
          localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)          
        END SELECT
      CASE DEFAULT
        localError="The specified partial derivative index of "//TRIM(NumberToVString(partialDerivativeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE(2)
      SELECT CASE(partialDerivativeIndex)
      CASE(NO_PART_DERIV)
        SELECT CASE(localNodeIndex)
        CASE(1)
          SELECT CASE(nodeDerivativeIndex)
          CASE(1)
            Hermite_QuadraticEvaluateDP=1.0_DP-xi*xi ! -xi^2+1
          CASE(2)
            Hermite_QuadraticEvaluateDP=xi*(1.0_DP-xi) ! -xi^2+xi
          CASE DEFAULT
            localError="The specified node derivative index of "//TRIM(NumberToVString(nodeDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)          
          END SELECT
        CASE(2)
          SELECT CASE(nodeDerivativeIndex)
          CASE(1)
            Hermite_QuadraticEvaluateDP=xi*xi ! xi^2
          CASE(2)
            Hermite_QuadraticEvaluateDP=0.0_DP ! 0
          CASE DEFAULT
            localError="The specified node derivative index of "//TRIM(NumberToVString(nodeDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)          
          END SELECT
        CASE DEFAULT
          localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)          
        END SELECT
      CASE(FIRST_PART_DERIV)
        SELECT CASE(localNodeIndex)
        CASE(1)
          SELECT CASE(nodeDerivativeIndex)
          CASE(1)
            Hermite_QuadraticEvaluateDP=-2.0_DP*xi ! -2xi
          CASE(2)
            Hermite_QuadraticEvaluateDP=1.0_DP-2.0_DP*xi ! -2xi+1
          CASE DEFAULT
            localError="The specified node derivative index of "//TRIM(NumberToVString(nodeDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)          
          END SELECT
        CASE(2)
          SELECT CASE(nodeDerivativeIndex)
          CASE(1)
            Hermite_QuadraticEvaluateDP=2.0_DP*xi ! 2xi
          CASE(2)
            Hermite_QuadraticEvaluateDP=0.0_DP ! 0
          CASE DEFAULT
            localError="The specified node derivative index of "//TRIM(NumberToVString(nodeDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)          
          END SELECT
        CASE DEFAULT
          localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)          
        END SELECT
      CASE(SECOND_PART_DERIV)
        SELECT CASE(localNodeIndex)
        CASE(1)
          SELECT CASE(nodeDerivativeIndex)
          CASE(1)
            Hermite_QuadraticEvaluateDP=-2.0_DP ! -2
          CASE(2)
            Hermite_QuadraticEvaluateDP=-2.0_DP ! -2
          CASE DEFAULT
            localError="The specified node derivative index of "//TRIM(NumberToVString(nodeDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)          
          END SELECT
        CASE(2)
          SELECT CASE(nodeDerivativeIndex)
          CASE(1)
            Hermite_QuadraticEvaluateDP=2.0_DP ! 2
          CASE(2)
            Hermite_QuadraticEvaluateDP=0.0_DP ! 0
          CASE DEFAULT
            localError="The specified node derivative index of "//TRIM(NumberToVString(nodeDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)          
          END SELECT
        CASE DEFAULT
          localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)          
        END SELECT
      CASE DEFAULT
        localError="The specified partial derivative index of "//TRIM(NumberToVString(partialDerivativeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE DEFAULT
        localError="The specified special node index of "//TRIM(NumberToVString(specialNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
    END SELECT

    EXITS("Hermite_QuadraticEvaluateDP")
    RETURN
999 ERRORSEXITS("Hermite_QuadraticEvaluateDP",err,error)
    RETURN
    
  END FUNCTION Hermite_QuadraticEvaluateDP

  !
  !================================================================================================================================
  !

  !>Evaluates a 1D cubic Lagrange basis function.
  FUNCTION Lagrange_CubicEvaluateDP(localNodeIndex,partialDerivativeIndex,xi,err,error)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: localNodeIndex !<The local node of the basis to evaluate. Must be between 1 and 4.
    INTEGER(INTG), INTENT(IN) :: partialDerivativeIndex !<The partial derivative to evaluate.
    REAL(DP), INTENT(IN) :: xi !<The Xi location to evaluate the basis at.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: Lagrange_CubicEvaluateDP !<On exit the evaluated basis function.
    !Local variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Lagrange_CubicEvaluateDP",err,error,*999)
    
    Lagrange_CubicEvaluateDP=0.0_DP
    
    SELECT CASE(partialDerivativeIndex)
    CASE(NO_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        Lagrange_CubicEvaluateDP=0.5_DP*(3.0_DP*xi-1.0_DP)*(3.0_DP*xi-2.0_DP)*(1.0_DP-xi) !
      CASE(2)
        Lagrange_CubicEvaluateDP=4.5_DP*xi*(3.0_DP*xi-2.0_DP)*(xi-1.0_DP) !
      CASE(3)
        Lagrange_CubicEvaluateDP=4.5_DP*xi*(3.0_DP*xi-1.0_DP)*(1.0_DP-xi) !
      CASE(4)
        Lagrange_CubicEvaluateDP=0.5_DP*xi*(3.0_DP*xi-1.0_DP)*(3.0_DP*xi-2.0_DP) !
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE(FIRST_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        Lagrange_CubicEvaluateDP=-13.5_DP*xi*xi+18.0_DP*xi-5.5_DP ! -13.5xi^2+18xi-5.5
      CASE(2)
        Lagrange_CubicEvaluateDP= 40.5_DP*xi*xi-45.0_DP*xi+9.0_DP ! 40.5xi^2-45xi+9
      CASE(3)
        Lagrange_CubicEvaluateDP=-40.5_DP*xi*xi+36.0_DP*xi-4.5_DP ! -40.5xi^2+36xi-4.5
      CASE(4)
        Lagrange_CubicEvaluateDP= 13.5_DP*xi*xi- 9.0_DP*xi+1.0_DP ! 13.5xi^2-9xi+1
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE(SECOND_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        Lagrange_CubicEvaluateDP=9.0_DP*(2.0_DP-3.0_DP*xi) ! 18-27xi
      CASE(2)
        Lagrange_CubicEvaluateDP=9.0_DP*(9.0_DP*xi-5.0_DP) ! 81xi-45
      CASE(3)
        Lagrange_CubicEvaluateDP=9.0_DP*(4.0_DP-9.0_DP*xi) ! 36-81xi
      CASE(4)
        Lagrange_CubicEvaluateDP=9.0_DP*(3.0_DP*xi-1.0_DP) ! 27xi-9
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE DEFAULT
      localError="The specified partial derivative index of "//TRIM(NumberToVString(partialDerivativeIndex,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)          
    END SELECT

    EXITS("Lagrange_CubicEvaluateDP")
    RETURN
999 ERRORSEXITS("Lagrange_CubicEvaluateDP",err,error)
    RETURN
    
  END FUNCTION Lagrange_CubicEvaluateDP

  !
  !================================================================================================================================
  !
  
  !> Evaluates a 1D linear Lagrange basis function.
  FUNCTION Lagrange_LinearEvaluateDP(localNodeIndex,partialDerivativeIndex,xi,err,error)
      
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: localNodeIndex !<The local node of the basis to evaluate. Must be between 1 and 2.
    INTEGER(INTG), INTENT(IN) :: partialDerivativeIndex !<The partial derivative to evaluate.
    REAL(DP), INTENT(IN) :: xi !<The Xi location to evaluate the basis at.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: Lagrange_LinearEvaluateDP !<On exit the evaluated basis function.
    !Local variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Lagrange_LinearEvaluateDP",err,error,*999)

    Lagrange_LinearEvaluateDP=0.0_DP
    SELECT CASE(partialDerivativeIndex)
    CASE(NO_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        Lagrange_LinearEvaluateDP=1.0_DP-xi ! 1-xi
      CASE(2)
        Lagrange_LinearEvaluateDP=xi !xi
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE(FIRST_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        Lagrange_LinearEvaluateDP=-1.0_DP ! -1
      CASE(2)
        Lagrange_LinearEvaluateDP=1.0_DP ! 1
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE(SECOND_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        Lagrange_LinearEvaluateDP=0.0_DP ! 0
      CASE(2)
        Lagrange_LinearEvaluateDP=0.0_DP ! 0
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE DEFAULT
      localError="The specified partial derivative index of "//TRIM(NumberToVString(partialDerivativeIndex,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)          
    END SELECT
    
    EXITS("Lagrange_LinearEvaluateDP")
    RETURN
999 ERRORSEXITS("Lagrange_LinearEvaluateDP",err,error)
    RETURN
    
  END FUNCTION Lagrange_LinearEvaluateDP

  !
  !================================================================================================================================
  !

  !>Evaluates a 1D quadratic Lagrange basis function.
  FUNCTION Lagrange_QuadraticEvaluateDP(localNodeIndex,partialDerivativeIndex,xi,err,error)
     
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: localNodeIndex !<The local node of the basis to evaluate. Must be between 1 and 3.
    INTEGER(INTG), INTENT(IN) :: partialDerivativeIndex !<The partial derivative to evaluate.
    REAL(DP), INTENT(IN) :: xi !<The Xi location to evaluate the basis at.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: Lagrange_QuadraticEvaluateDP !<On exit the evaluated basis function.
    !Local variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Lagrange_QuadraticEvaluateDP",err,error,*999)

    Lagrange_QuadraticEvaluateDP=0.0_DP
    SELECT CASE(partialDerivativeIndex)
    CASE(NO_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        Lagrange_QuadraticEvaluateDP=1.0_DP-3.0_DP*xi+2.0_DP*xi*xi ! 1-3xi+2xi^2
      CASE(2)
        Lagrange_QuadraticEvaluateDP=4.0_DP*xi*(1.0_DP-xi) ! 4xi-4xi^2
      CASE(3)
        Lagrange_QuadraticEvaluateDP=xi*(xi+xi-1.0_DP) ! 2xi^2-xi
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE(FIRST_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        Lagrange_QuadraticEvaluateDP=4.0_DP*xi-3.0_DP ! 4xi-3
      CASE(2)
        Lagrange_QuadraticEvaluateDP=4.0_DP-8.0_DP*xi ! 4-8xi
      CASE(3)
        Lagrange_QuadraticEvaluateDP=4.0_DP*xi-1.0_DP ! 4xi-1
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE(SECOND_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        Lagrange_QuadraticEvaluateDP=4.0_DP ! 4
      CASE(2)
        Lagrange_QuadraticEvaluateDP=-8.0_DP ! -8
      CASE(3)
        Lagrange_QuadraticEvaluateDP=4.0_DP ! 4
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE DEFAULT
      localError="The specified partial derivative index of "//TRIM(NumberToVString(partialDerivativeIndex,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)          
    END SELECT

    EXITS("Lagrange_QuadraticEvaluateDP")
    RETURN
999 ERRORSEXITS("Lagrange_QuadraticEvaluateDP",err,error)
    RETURN
    
  END FUNCTION Lagrange_QuadraticEvaluateDP

  !
  !================================================================================================================================
  !
  
  !>Evaluates a cubic simpelx basis function at a specificed area position and node index and with a given partial derivative index
  !>with respect to area coordinates for double precision arguments
  FUNCTION Simplex_CubicEvaluateDP(localNodeIndex,partialDerivativeIndex,xl,err,error)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: localNodeIndex !<The node index to evaluate
    INTEGER(INTG), INTENT(IN) :: partialDerivativeIndex !<The partial derivative index wrt area coordinates to evaluate
    REAL(DP), INTENT(IN) :: xl !<The area coordinate to evaluate at.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: Simplex_CubicEvaluateDP
    !Local variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Simplex_CubicEvaluateDP",err,error,*999)
    
    Simplex_CubicEvaluateDP=0.0_DP
        
    SELECT CASE(partialDerivativeIndex)
    CASE(NO_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        Simplex_CubicEvaluateDP=1.0 !1
      CASE(2)
        Simplex_CubicEvaluateDP=3.0_DP*xl !3L
      CASE(3)
        Simplex_CubicEvaluateDP=3.0_DP/2.0_DP*xl*(3.0_DP*xl-1.0_DP) !3/2.L(3L-1)
      CASE(4)
        Simplex_CubicEvaluateDP=0.5_DP*xl*(3.0_DP*xl-1.0_DP)*(3.0_DP*xl-2.0_DP) !1/2.L(3L-1)(3L-2)
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE(FIRST_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        Simplex_CubicEvaluateDP=0.0_DP !0
      CASE(2)
        Simplex_CubicEvaluateDP=3.0_DP !3
      CASE(3)
        Simplex_CubicEvaluateDP=3.0_DP/2.0_DP*(6.0_DP*xl-1) !3/2.(6L-1)
      CASE(4)
        Simplex_CubicEvaluateDP=13.5_DP*xl*xl-9.0_DP*xl+1.0_DP !27/2.L^2-9L+1
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE(SECOND_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        Simplex_CubicEvaluateDP=0.0_DP !0
      CASE(2)
        Simplex_CubicEvaluateDP=0.0_DP !0
      CASE(3)
        Simplex_CubicEvaluateDP=9.0_DP !9
      CASE(4)
        Simplex_CubicEvaluateDP=2.0_DP*xl-9.0_DP
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE(THIRD_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        Simplex_CubicEvaluateDP=0.0_DP !0
      CASE(2)
        Simplex_CubicEvaluateDP=0.0_DP !0
      CASE(3)
        Simplex_CubicEvaluateDP=0.0_DP !0
      CASE(4)
        Simplex_CubicEvaluateDP=2.0_DP !2
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE DEFAULT
      localError="The specified partial derivative index of "//TRIM(NumberToVString(partialDerivativeIndex,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)          
    END SELECT

    EXITS("Simplex_CubicEvaluateDP")
    RETURN
999 ERRORSEXITS("Simplex_CubicEvaluateDP",err,error)
    RETURN
    
  END FUNCTION Simplex_CubicEvaluateDP

  !
  !================================================================================================================================
  !
  
  !>Evaluates a linear simpelx basis function at a specificed area position and node index and with a given partial derivative index
  !>with respect to area coordinates for double precision arguments
  FUNCTION Simplex_LinearEvaluateDP(localNodeIndex,partialDerivativeIndex,xl,err,error)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: localNodeIndex !<The node index to evaluate
    INTEGER(INTG), INTENT(IN) :: partialDerivativeIndex !<The partial derivative index wrt area coordinates to evaluate
    REAL(DP), INTENT(IN) :: xl !<The area coordinate to evaluate at.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: Simplex_LinearEvaluateDP
    !Local variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Simplex_LinearEvaluateDP",err,error,*999)

    Simplex_LinearEvaluateDP=0.0_DP
    SELECT CASE(partialDerivativeIndex)
    CASE(NO_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        Simplex_LinearEvaluateDP=1.0 !1
      CASE(2)
        Simplex_LinearEvaluateDP=xl  !L
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE(FIRST_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        Simplex_LinearEvaluateDP=0.0_DP  !0
      CASE(2)
        Simplex_LinearEvaluateDP=1.0_DP !1
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE(SECOND_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        Simplex_LinearEvaluateDP=0.0_DP !0
      CASE(2)
        Simplex_LinearEvaluateDP=0.0_DP !0
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE(THIRD_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        Simplex_LinearEvaluateDP=0.0_DP !0
      CASE(2)
        Simplex_LinearEvaluateDP=0.0_DP !0
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE DEFAULT
      localError="The specified partial derivative index of "//TRIM(NumberToVString(partialDerivativeIndex,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)          
    END SELECT
    
    EXITS("Simplex_LinearEvaluateDP")
    RETURN
999 ERRORSEXITS("Simplex_LinearEvaluateDP",err,error)
    RETURN
    
  END FUNCTION Simplex_LinearEvaluateDP

  !
  !================================================================================================================================
  !
  
  !>Evaluates a quadratic simpelx basis function at a specificed area position and node index and with a given partial derivative index
  !>with respect to area coordinates for double precision arguments.
   FUNCTION Simplex_QuadraticEvaluateDP(localNodeIndex,partialDerivativeIndex,xl,err,error)
      
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: localNodeIndex !<The node index to evaluate
    INTEGER(INTG), INTENT(IN) :: partialDerivativeIndex !<The partial derivative index wrt area coordinates to evaluate
    REAL(DP), INTENT(IN) :: xl !<The area coordinate to evaluate at.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: Simplex_QuadraticEvaluateDP
    !Local variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Simplex_QuadraticEvaluateDP",err,error,*999)

    Simplex_QuadraticEvaluateDP=0.0_DP
    
    SELECT CASE(partialDerivativeIndex)
    CASE(NO_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        Simplex_QuadraticEvaluateDP=1.0_DP !1
      CASE(2)
        Simplex_QuadraticEvaluateDP=2.0_DP*xl !2L
      CASE(3)
        Simplex_QuadraticEvaluateDP=xl*(2.0_DP*xl-1.0_DP) !L(2L-1)
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE(FIRST_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        Simplex_QuadraticEvaluateDP=0.0_DP !0
      CASE(2)
        Simplex_QuadraticEvaluateDP=2.0_DP !4
      CASE(3)
        Simplex_QuadraticEvaluateDP=4.0_DP*xl-1.0_DP !4L-1
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE(SECOND_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        Simplex_QuadraticEvaluateDP=0.0_DP !0
      CASE(2)
        Simplex_QuadraticEvaluateDP=0.0_DP !0
      CASE(3)
        Simplex_QuadraticEvaluateDP=4.0_DP !4
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE(THIRD_PART_DERIV)
      SELECT CASE(localNodeIndex)
      CASE(1)
        Simplex_QuadraticEvaluateDP=0.0_DP !0
      CASE(2)
        Simplex_QuadraticEvaluateDP=0.0_DP !0
      CASE(3)
        Simplex_QuadraticEvaluateDP=0.0_DP !0
      CASE DEFAULT
        localError="The specified local node index of "//TRIM(NumberToVString(localNodeIndex,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)          
      END SELECT
    CASE DEFAULT
      localError="The specified partial derivative index of "//TRIM(NumberToVString(partialDerivativeIndex,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)          
    END SELECT

    EXITS("Simplex_QuadraticEvaluateDP")
    RETURN
999 ERRORSEXITS("Simplex_QuadraticEvaluateDP",err,error)
    RETURN
    
  END FUNCTION Simplex_QuadraticEvaluateDP

  !
  !================================================================================================================================
  !

  !>Finalises the basis functions for a context and deallocates all memory
  SUBROUTINE BasisFunctions_Finalise(basisFunctions,err,error,*)

    !Argument variables
    TYPE(BasisFunctionsType), POINTER :: basisFunctions !<A pointer to the basis functions to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("BasisFunctions_Finalise",err,error,*999)

    IF(ASSOCIATED(basisFunctions)) THEN
      !Destroy any created basis functions
      DO WHILE(basisFunctions%numberOfBasisFunctions>0)
        CALL Basis_Destroy(basisFunctions%bases(1)%ptr,err,error,*999)
      ENDDO !nb
      !Destroy basis functions and deallocated any memory allocated
      IF(ALLOCATED(basisFunctions%bases)) DEALLOCATE(basisFunctions%bases)
      DEALLOCATE(basisFunctions)
    ENDIF
    
    EXITS("BasisFunctions_Finalise")
    RETURN
999 ERRORSEXITS("BasisFunctions_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE BasisFunctions_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the basis functions for a context.
  SUBROUTINE BasisFunctions_Initialise(context,err,error,*)

    !Argument variables
    TYPE(ContextType), POINTER :: context !<The context to intialise the basis functions for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("BasisFunctions_Initialise",err,error,*999)

    IF(.NOT.ASSOCIATED(context)) CALL FlagError("Context is not associated.",err,error,*998)
    IF(ASSOCIATED(context%basisFunctions)) CALL FlagError("Context basis functions is already associated.",err,error,*998)

    ALLOCATE(context%basisFunctions,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate basis functions.",err,error,*999)
    !Initialise
    context%basisFunctions%context=>context    
    context%basisFunctions%numberOfBasisFunctions=0    
    
    EXITS("BasisFunctions_Initialise")
    RETURN
999 CALL BasisFunctions_Finalise(context%basisFunctions,dummyErr,dummyError,*998)
998 ERRORSEXITS("BasisFunctions_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE BasisFunctions_Initialise

  !
  !================================================================================================================================
  !

END MODULE BasisRoutines

