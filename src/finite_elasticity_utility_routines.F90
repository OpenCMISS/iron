!> \file
!> \author Chris Bradley
!> \brief This module handles utility routines for finite elasticity to allow access from other modules.
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

!>This module handles utility routines for finite elasticity to allow access from other modules.
MODULE FiniteElasticityUtilityRoutines

  USE BaseRoutines
  USE CoordinateSystemRoutines
  USE FieldAccessRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE Maths
  USE Strings
  USE Types
  
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

  PUBLIC FiniteElasticity_DeformationGradientTensorCalculate

CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the deformation gradient tensor at a given interpolated point
  SUBROUTINE FiniteElasticity_DeformationGradientTensorCalculate(dependentInterpPointMetrics,geometricInterpPointMetrics,&
    & fibreInterpPoint,F,J,err,error,*)

    !Argument variables
    TYPE(FieldInterpolatedPointMetricsType), POINTER, INTENT(IN) :: dependentInterpPointMetrics !<The interpolated point metrics of the deformed/spatial geometry
    TYPE(FieldInterpolatedPointMetricsType), POINTER, INTENT(IN) :: geometricInterpPointMetrics !<The interpolated point metrics of the undeformed/reference geometry
    TYPE(FieldInterpolatedPointType), POINTER, INTENT(IN) :: fibreInterpPoint !<A pointer to the fibre field interpolated at the point. If there is no fibre field defined this will be NULL.
    REAL(DP), INTENT(OUT) :: F(:,:) !<F(dependentCoordinateIdx,geometricCoordinateIdx). On return, the deformation gradient tensor F
    REAL(DP), INTENT(OUT) :: J !<On return, the Jacobian of the deformation i.e., determinant F
    INTEGER(INTG), INTENT(OUT) :: err   !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfXDimensions,numberOfXiDimensions,numberOfZDimensions
    REAL(DP) :: dNudX(3,3),dXdNu(3,3),dNudXi(3,3),dXidNu(3,3)
    TYPE(VARYING_STRING) :: localError

    ENTERS("FiniteElasticity_DeformationGradientTensorCalculate",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dependentInterpPointMetrics)) &
      & CALL FlagError("Dependent interpolated point metrics is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(geometricInterpPointMetrics)) &
      & CALL FlagError("Geometric interpolated point metrics is not associated.",err,error,*999)
    IF(SIZE(F,1)<dependentInterpPointMetrics%numberOfXDimensions) THEN
      localError="The size of the first index of F is too small. The size should be >= "// &
        & TRIM(NumberToVString(dependentInterpPointMetrics%numberOfXDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(F,2)<geometricInterpPointMetrics%numberOfXDimensions) THEN
      localError="The size of the second index of F is too small. The size should be >= "// &
        & TRIM(NumberToVString(geometricInterpPointMetrics%numberOfXDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    CALL FieldInterpolatedPointMetrics_NumberOfXDimensionsGet(geometricInterpPointMetrics,numberOfXDimensions,err,error,*999)
    CALL FieldInterpolatedPointMetrics_NumberOfXiDimensionsGet(geometricInterpPointMetrics,numberOfXiDimensions,err,error,*999)
    CALL FieldInterpolatedPointMetrics_NumberOfXDimensionsGet(dependentInterpPointMetrics,numberOfZDimensions,err,error,*999)

    CALL Identity(F,err,error,*999)
    
    CALL CoordinateSystem_MaterialSystemCalculate(geometricInterpPointMetrics,fibreInterpPoint, &
      & dNudXi(1:numberOfXDimensions,1:numberOfXiDimensions),dXidNu(1:numberOfXiDimensions,1:numberOfXDimensions), &
      & err,error,*999)
    !F = dZ/dNu = dZ/dXi * dXi/dNu  (deformation gradient tensor, F)
    CALL MatrixProduct(dependentInterpPointMetrics%dXdXi(1:numberOfZDimensions,1:numberOfXiDimensions), &
      & dXidNu(1:numberOfXiDimensions,1:numberOfXDimensions),F(1:numberOfZDimensions,1:numberOfXDimensions), &
      & err,error,*999)
    
    CALL Determinant(F(1:numberOfZDimensions,1:numberOfXDimensions),J,err,error,*999)

    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Calculated deformation gradient tensor:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of Z dimensions  = ",numberOfZDimensions,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of Xi dimensions = ",numberOfXiDimensions,err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Derivative of Xi wrt to Nu coordinates:",err,error,*999)
      CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfXiDimensions,1,1,numberOfXDimensions, &
        & numberOfXDimensions,numberOfXDimensions,dXidNu,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
        & '("    dXidNu','(",I1,",:)','   :",3(X,E13.6))','(19X,3(X,E13.6))',err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Deformation gradient tensor wrt Nu coordinates:",err,error,*999)
      CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfZDimensions,1,1,numberOfXDimensions, &
        & numberOfXDimensions,numberOfXDimensions,F,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
        & '("         F','(",I1,",:)','   :",3(X,E13.6))','(19X,3(X,E13.6))',err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Determinant F  = ",J,err,error,*999)
    ENDIF

    EXITS("FiniteElasticity_DeformationGradientTensorCalculate")
    RETURN
999 ERRORS("FiniteElasticity_GaussDeformationGradientTensorCalculate",err,error)
    EXITS("FiniteElasticity_GaussDeformationGradientTensorCalculate")
    RETURN 1

  END SUBROUTINE FiniteElasticity_DeformationGradientTensorCalculate

  !
  !================================================================================================================================
  !

  
END MODULE FiniteElasticityUtilityRoutines
