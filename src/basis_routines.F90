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
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> This module contains all basis function routines.
MODULE BasisRoutines

  USE BaseRoutines
  USE BasisAccessRoutines
  USE CONSTANTS
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE STRINGS
  USE TYPES

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup BASIS_ROUTINES_BasisTypes BASIS_ROUTINES::BasisTypes
  !> \brief Basis definition type parameters
  !> \todo Combine simplex and serendipity elements???
  !> \see BASIS_ROUTINES,OPENCMISS_BasisTypes
  !>@{ 
  INTEGER(INTG), PARAMETER :: BASIS_LAGRANGE_HERMITE_TP_TYPE=1 !<Lagrange-Hermite tensor product basis type \see BASIS_ROUTINES_BasisTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_SIMPLEX_TYPE=2 !<Simplex basis type \see BASIS_ROUTINES_BasisTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_SERENDIPITY_TYPE=3 !<Serendipity basis type \see BASIS_ROUTINES_BasisTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_AUXILLIARY_TYPE=4 !<Auxillary basis type \see BASIS_ROUTINES_BasisTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_B_SPLINE_TP_TYPE=5 !<B-spline basis type \see BASIS_ROUTINES_BasisTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE=6 !<Fourier-Lagrange tensor product basis type \see BASIS_ROUTINES_BasisTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_EXTENDED_LAGRANGE_TP_TYPE=7 !< Extendend Lagrange tensor product basis type \see BASIS_ROUTINES_BasisTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_RADIAL_TYPE=7 !< Radial basis typee \see BASIS_ROUTINES_BasisTypes,BASIS_ROUTINES
  !>@}

  !> \addtogroup BASIS_ROUTINES_InterpolationSpecifications BASIS_ROUTINES::InterpolationSpecifications
  !> \brief Interpolation specification parameters
  !> \see BASIS_ROUTINES,OPENCMISS_InterpolationSpecifications
  !>@{ 
  INTEGER(INTG), PARAMETER :: BASIS_LINEAR_LAGRANGE_INTERPOLATION=1 !<Linear Lagrange interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_QUADRATIC_LAGRANGE_INTERPOLATION=2 !<Quadratic Lagrange interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_CUBIC_LAGRANGE_INTERPOLATION=3 !<Cubic Lagrange interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_CUBIC_HERMITE_INTERPOLATION=4 !<Cubic Hermite interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_QUADRATIC1_HERMITE_INTERPOLATION=5 !<Quadratic Hermite (no derivative at xi=0) interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_QUADRATIC2_HERMITE_INTERPOLATION=6 !<Quadratic Hermite (no derivative at xi=1) interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_LINEAR_SIMPLEX_INTERPOLATION=7 !<Linear Simplex interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_QUADRATIC_SIMPLEX_INTERPOLATION=8 !<Quadratic Simplex interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_CUBIC_SIMPLEX_INTERPOLATION=9 !<Cubic Simplex interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_GAUSSIAN_RADIAL_INTERPOLATION=10 !<Gaussian Radial interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_MULTIQUARTIC_RADIAL_INTERPOLATION=11 !<Multiquartic Radial interpolation specification \see BASIS_ROUTINES_InterpolationSpecifications,BASIS_ROUTINES
  !>@}

  !> \addtogroup BASIS_ROUTINES_InterpolationTypes BASIS_ROUTINES::InterpolationTypes
  !> \brief Interpolation type parameters for a Xi direction
  !> \see BASIS_ROUTINES
  !>@{ 
  INTEGER(INTG), PARAMETER :: BASIS_LAGRANGE_INTERPOLATION=1 !<Lagrange interpolation \see BASIS_ROUTINES_InterpolationTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_HERMITE_INTERPOLATION=2 !<Hermite interpolation \see BASIS_ROUTINES_InterpolationTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_SIMPLEX_INTERPOLATION=3 !<Simplex interpolation \see BASIS_ROUTINES_InterpolationTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_SERENDIPITY_INTERPOLATION=4 !<Serendipity interpolation \see BASIS_ROUTINES_InterpolationTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_TRANSITION_INTERPOLATION=5 !<Transition interpolation \see BASIS_ROUTINES_InterpolationTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_SINGULAR_INTERPOLATION=6 !<Singular interpolation \see BASIS_ROUTINES_InterpolationTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_FOURIER_INTERPOLATION=7 !<Fourier interpolation \see BASIS_ROUTINES_InterpolationTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_RADIAL_INTERPOLATION=8 !<Radial interpolation \see BASIS_ROUTINES_InterpolationTypes,BASIS_ROUTINES
  !>@}
  
  !> \addtogroup BASIS_ROUTINES_InterpolationOrder BASIS_ROUTINES::InterpolationOrder
  !> \brief Interpolation order for a Xi direction
  !> \see BASIS_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: BASIS_LINEAR_INTERPOLATION_ORDER=1 !<Linear interpolation order \see BASIS_ROUTINES_InterpolationOrder,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_QUADRATIC_INTERPOLATION_ORDER=2 !<Quadratic interpolation order \see BASIS_ROUTINES_InterpolationOrder,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_CUBIC_INTERPOLATION_ORDER=3 !<Cubic interpolation order \see BASIS_ROUTINES_InterpolationOrder,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_QUADRATIC1_INTERPOLATION_ORDER=4 !<Quadratic (no derivative at xi=0) interpolation order \see BASIS_ROUTINES_InterpolationOrder,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_QUADRATIC2_INTERPOLATION_ORDER=5 !<Quadratic (no derivative at xi=1) interpolation order \see BASIS_ROUTINES_InterpolationOrder,BASIS_ROUTINES
  !>@}
  
  !> \addtogroup BASIS_ROUTINES_QuadratureTypes BASIS_ROUTINES::QuadratureTypes
  !> \brief Quadrature type parameters
  !> \see BASIS_ROUTINES,OPENCMISS_QuadratureTypes
  !>@{
  INTEGER(INTG), PARAMETER :: BASIS_GAUSS_LEGENDRE_QUADRATURE=1 !<Gauss-Legendre quadrature  \see BASIS_ROUTINES_QuadratureTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_GAUSS_LAGUERRE_QUADRATURE=2 !<Gauss-Laguerre quadrature  \see BASIS_ROUTINES_QuadratureTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_GUASS_HERMITE_QUADRATURE=3 !<Gauss-Hermite quadrature  \see BASIS_ROUTINES_QuadratureTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_ADAPTIVE_GAUSS_LEGENDRE_QUADRATURE=4 !<Adaptive Gauss-Legendre quadrature  \see BASIS_ROUTINES_QuadratureTypes,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_GAUSS_SIMPLEX_QUADRATURE=5 !<Gauss-Legendre for Simplex elements quadrature  \see BASIS_ROUTINES_QuadratureTypes,BASIS_ROUTINES
  !>@}

  !> \addtogroup BASIS_ROUTINES_XiCollapse BASIS_ROUTINES::XiCollapse
  !> \brief Xi collapse parameters
  !> \see BASIS_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: BASIS_XI_COLLAPSED=1 !<The Xi direction is collapsed \see BASIS_ROUTINES_XiCollapse,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_COLLAPSED_AT_XI0=2 !<The Xi direction at the xi=0 end of this Xi direction is collapsed \see BASIS_ROUTINES_XiCollapse,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_COLLAPSED_AT_XI1=3 !<The Xi direction at the xi=1 end of this Xi direction is collapsed \see BASIS_ROUTINES_XiCollapse,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_NOT_COLLAPSED=4 !<The Xi direction is not collapsed \see BASIS_ROUTINES_XiCollapse,BASIS_ROUTINES
  !>@}
  
  !> \addtogroup BASIS_ROUTINES_BoundaryXiType BASIS_ROUTINES::BoundaryXiType
  !> \brief The boundary xi types
  !> \see BASIS_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: BASIS_NO_BOUNDARY_XI=0 !<The xi is not on a boundary \see BASIS_ROUTINES_XiBoundaryType,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_LINE_BOUNDARY_XI=1 !<The xi is on a boundary line \see BASIS_ROUTINES_XiBoundaryType,BASIS_ROUTINES
  INTEGER(INTG), PARAMETER :: BASIS_FACE_BOUNDARY_XI=2 !<The xi is on a boundary face \see BASIS_ROUTINES_XiBoundaryType,BASIS_ROUTINES
  !>@}
  
  !!Module types
  ! 
  !!>Contains information on the defined basis functions
  !TYPE BASIS_FUNCTIONS_TYPE
  !  PRIVATE
  !  INTEGER(INTG) :: numberOfBasisFunctions !<The number of basis functions definegd
  !  TYPE(BASIS_PTR_TYPE), POINTER :: BASES(:) !<The array of pointers to the defined basis functions
  !END TYPE BASIS_FUNCTIONS_TYPE
  
  !Module variables

  !Interfaces

  INTERFACE Basis_CollapsedXiGet
    MODULE PROCEDURE BASIS_COLLAPSED_XI_GET
  END INTERFACE Basis_CollapsedXiGet
  
  !>Sets/changes the collapsed Xi flags for a basis.
  INTERFACE BASIS_COLLAPSED_XI_SET
    MODULE PROCEDURE BASIS_COLLAPSED_XI_SET_NUMBER
    MODULE PROCEDURE BASIS_COLLAPSED_XI_SET_PTR
  END INTERFACE BASIS_COLLAPSED_XI_SET

  !>Sets/changes the collapsed Xi flags for a basis.
  INTERFACE Basis_CollapsedXiSet
    MODULE PROCEDURE BASIS_COLLAPSED_XI_SET_NUMBER
    MODULE PROCEDURE BASIS_COLLAPSED_XI_SET_PTR
  END INTERFACE Basis_CollapsedXiSet

  INTERFACE BASIS_CREATE_START
    MODULE PROCEDURE Basis_CreateStart
  END INTERFACE BASIS_CREATE_START

  INTERFACE BASIS_CREATE_FINISH
    MODULE PROCEDURE Basis_CreateFinish
  END INTERFACE BASIS_CREATE_FINISH

  !>Evaluates the appropriate partial derivative index for the specificied basis function at a Xi location \see BASIS_ROUTINES
  INTERFACE BASIS_EVALUATE_XI
    MODULE PROCEDURE BASIS_EVALUATE_XI_DP
  END INTERFACE BASIS_EVALUATE_XI

  !>Evaluates the appropriate partial derivative index for the specificied basis function at a Xi location \see BASIS_ROUTINES
  INTERFACE Basis_EvaluateXi
    MODULE PROCEDURE BASIS_EVALUATE_XI_DP
  END INTERFACE Basis_EvaluateXi

  !>Interpolates the appropriate partial derivative index of the elements parameters for basis function at a Gauss point \see BASIS_ROUTINES
  INTERFACE BASIS_INTERPOLATE_GAUSS
    MODULE PROCEDURE BASIS_INTERPOLATE_GAUSS_DP
  END INTERFACE BASIS_INTERPOLATE_GAUSS
  
  !>Interpolates the appropriate partial derivative index of the elements parameters for basis function at a Gauss point \see BASIS_ROUTINES
  INTERFACE Basis_InterpolateGauss
    MODULE PROCEDURE BASIS_INTERPOLATE_GAUSS_DP
  END INTERFACE Basis_InterpolateGauss
  
  !>Interpolates the appropriate partial derivative index of the elements parameters for basis function at Xi location \see BASIS_ROUTINES
  INTERFACE BASIS_INTERPOLATE_XI
    MODULE PROCEDURE BASIS_INTERPOLATE_XI_DP
  END INTERFACE BASIS_INTERPOLATE_XI

   !>Interpolates the appropriate partial derivative index of the elements parameters for basis function at Xi location \see BASIS_ROUTINES
  INTERFACE Basis_InterpolateXi
    MODULE PROCEDURE BASIS_INTERPOLATE_XI_DP
  END INTERFACE Basis_InterpolateXi

  !>Interpolates the requested partial derivative index(ices) of the element parameters for basis function at a face Gauss point \see BASIS_ROUTINES
  INTERFACE BASIS_INTERPOLATE_LOCAL_FACE_GAUSS
    MODULE PROCEDURE BASIS_INTERPOLATE_LOCAL_FACE_GAUSS_DP
  END INTERFACE BASIS_INTERPOLATE_LOCAL_FACE_GAUSS

   !>Interpolates the requested partial derivative index(ices) of the element parameters for basis function at a face Gauss point \see BASIS_ROUTINES
  INTERFACE Basis_InterpolateLocalFaceGauss
    MODULE PROCEDURE BASIS_INTERPOLATE_LOCAL_FACE_GAUSS_DP
  END INTERFACE Basis_InterpolateLocalFaceGauss

  INTERFACE Basis_LocalNodeXiCalculate
    MODULE PROCEDURE BASIS_LOCAL_NODE_XI_CALCULATE
  END INTERFACE Basis_LocalNodeXiCalculate

  INTERFACE Basis_InterpolationXiGet
    MODULE PROCEDURE BASIS_INTERPOLATION_XI_GET
  END INTERFACE Basis_InterpolationXiGet
  
  !>Sets/changes the interpolation type in each Xi direction for a basis
  INTERFACE BASIS_INTERPOLATION_XI_SET
    MODULE PROCEDURE BASIS_INTERPOLATION_XI_SET_NUMBER
    MODULE PROCEDURE BASIS_INTERPOLATION_XI_SET_PTR
  END INTERFACE BASIS_INTERPOLATION_XI_SET

   !>Sets/changes the interpolation type in each Xi direction for a basis
  INTERFACE Basis_InterpolationXiSet
    MODULE PROCEDURE BASIS_INTERPOLATION_XI_SET_NUMBER
    MODULE PROCEDURE BASIS_INTERPOLATION_XI_SET_PTR
  END INTERFACE Basis_InterpolationXiSet

  INTERFACE Basis_NumberOfLocalNodesGet
    MODULE PROCEDURE BASIS_NUMBER_OF_LOCAL_NODES_GET
  END INTERFACE Basis_NumberOfLocalNodesGet

  INTERFACE BASIS_NUMBER_OF_XI_GET
    MODULE PROCEDURE Basis_NumberOfXiGet
  END INTERFACE BASIS_NUMBER_OF_XI_GET
  
  INTERFACE BASIS_NUMBER_OF_XI_SET
    MODULE PROCEDURE Basis_NumberOfXiSet
  END INTERFACE BASIS_NUMBER_OF_XI_SET

  INTERFACE Basis_QuadratureDestroy
    MODULE PROCEDURE BASIS_QUADRATURE_DESTROY
  END INTERFACE Basis_QuadratureDestroy

  INTERFACE Basis_QuadratureMultipleGaussXiGet
    MODULE PROCEDURE BASIS_QUADRATURE_MULTIPLE_GAUSS_XI_GET
  END INTERFACE Basis_QuadratureMultipleGaussXiGet
  
  INTERFACE Basis_QuadratureNumberOfGaussXiGet
    MODULE PROCEDURE BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_GET
  END INTERFACE Basis_QuadratureNUmberOfGaussXiGet

  INTERFACE Basis_QuadratureNumberOfGaussXiSet
    MODULE PROCEDURE BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET
  END INTERFACE Basis_QuadratureNUmberOfGaussXiSet

  INTERFACE Basis_QuadratureOrderGet
    MODULE PROCEDURE BASIS_QUADRATURE_ORDER_GET
  END INTERFACE Basis_QuadratureOrderGet
  
  !>Sets/changes the order of a quadrature for a basis quadrature.
  INTERFACE BASIS_QUADRATURE_ORDER_SET
    MODULE PROCEDURE BASIS_QUADRATURE_ORDER_SET_NUMBER
    MODULE PROCEDURE BASIS_QUADRATURE_ORDER_SET_PTR
  END INTERFACE BASIS_QUADRATURE_ORDER_SET

  !>Sets/changes the order of a quadrature for a basis quadrature.
  INTERFACE Basis_QuadratureOrderSet
    MODULE PROCEDURE BASIS_QUADRATURE_ORDER_SET_NUMBER
    MODULE PROCEDURE BASIS_QUADRATURE_ORDER_SET_PTR
  END INTERFACE Basis_QuadratureOrderSet

  INTERFACE Basis_QuadratureSingleGaussXiGet
    MODULE PROCEDURE BASIS_QUADRATURE_SINGLE_GAUSS_XI_GET
  END INTERFACE Basis_QuadratureSingleGaussXiGet

  INTERFACE Basis_QuadratureTypeGet
    MODULE PROCEDURE BASIS_QUADRATURE_TYPE_GET
  END INTERFACE Basis_QuadratureTypeGet
  
  !>Sets/changes the quadrature type for a basis
  INTERFACE BASIS_QUADRATURE_TYPE_SET
    MODULE PROCEDURE BASIS_QUADRATURE_TYPE_SET_NUMBER
    MODULE PROCEDURE BASIS_QUADRATURE_TYPE_SET_PTR
  END INTERFACE BASIS_QUADRATURE_TYPE_SET
  
  !>Sets/changes the quadrature type for a basis
  INTERFACE Basis_QuadratureTypeSet
    MODULE PROCEDURE BASIS_QUADRATURE_TYPE_SET_NUMBER
    MODULE PROCEDURE BASIS_QUADRATURE_TYPE_SET_PTR
  END INTERFACE Basis_QuadratureTypeSet

  INTERFACE Basis_TypeGet
    MODULE PROCEDURE BASIS_TYPE_GET
  END INTERFACE Basis_TypeGet
  
  !>Sets/changes the type for a basis.
  INTERFACE BASIS_TYPE_SET
    MODULE PROCEDURE BASIS_TYPE_SET_NUMBER
    MODULE PROCEDURE BASIS_TYPE_SET_PTR
  END INTERFACE BASIS_TYPE_SET

  !>Sets/changes the type for a basis.
  INTERFACE Basis_TypeSet
    MODULE PROCEDURE BASIS_TYPE_SET_NUMBER
    MODULE PROCEDURE BASIS_TYPE_SET_PTR
  END INTERFACE Basis_TypeSet

  !>Evaluates a linear Simplex basis function
  INTERFACE SIMPLEX_LINEAR_EVALUATE
    MODULE PROCEDURE SIMPLEX_LINEAR_EVALUATE_DP
  END INTERFACE SIMPLEX_LINEAR_EVALUATE

  !>Evaluates a quadratic Simplex basis function
  INTERFACE SIMPLEX_QUADRATIC_EVALUATE
    MODULE PROCEDURE SIMPLEX_QUADRATIC_EVALUATE_DP
  END INTERFACE SIMPLEX_QUADRATIC_EVALUATE

  !>Evaluates a cubic Simplex basis function
  INTERFACE SIMPLEX_CUBIC_EVALUATE
    MODULE PROCEDURE SIMPLEX_CUBIC_EVALUATE_DP
  END INTERFACE SIMPLEX_CUBIC_EVALUATE

  !>Evaluates the Lagrange/Hermite/Fourier tensor product basis function for the given basis
  INTERFACE BASIS_LHTP_BASIS_EVALUATE
    MODULE PROCEDURE BASIS_LHTP_BASIS_EVALUATE_DP
  END INTERFACE BASIS_LHTP_BASIS_EVALUATE

  !>Evaluates the Lagrange/Hermite/Fourier tensor product basis function for the given basis
  INTERFACE Basis_LHTPBasisEvaluate
    MODULE PROCEDURE BASIS_LHTP_BASIS_EVALUATE_DP
  END INTERFACE Basis_LHTPBasisEvaluate

  PUBLIC BASIS_LAGRANGE_HERMITE_TP_TYPE,BASIS_SIMPLEX_TYPE,BASIS_SERENDIPITY_TYPE,BASIS_AUXILLIARY_TYPE,BASIS_B_SPLINE_TP_TYPE, &
    & BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE,BASIS_EXTENDED_LAGRANGE_TP_TYPE, BASIS_RADIAL_TYPE

  PUBLIC BASIS_LINEAR_LAGRANGE_INTERPOLATION,BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,BASIS_CUBIC_LAGRANGE_INTERPOLATION, &
    & BASIS_CUBIC_HERMITE_INTERPOLATION,BASIS_QUADRATIC1_HERMITE_INTERPOLATION,BASIS_QUADRATIC2_HERMITE_INTERPOLATION, &
    & BASIS_LINEAR_SIMPLEX_INTERPOLATION,BASIS_QUADRATIC_SIMPLEX_INTERPOLATION,BASIS_CUBIC_SIMPLEX_INTERPOLATION, &
    & BASIS_GAUSSIAN_RADIAL_INTERPOLATION, BASIS_MULTIQUARTIC_RADIAL_INTERPOLATION

  PUBLIC BASIS_LINEAR_INTERPOLATION_ORDER,BASIS_QUADRATIC_INTERPOLATION_ORDER,BASIS_CUBIC_INTERPOLATION_ORDER, &
    & BASIS_QUADRATIC1_INTERPOLATION_ORDER,BASIS_QUADRATIC2_INTERPOLATION_ORDER
  
  PUBLIC BASIS_GAUSS_LEGENDRE_QUADRATURE,BASIS_GAUSS_LAGUERRE_QUADRATURE,BASIS_GUASS_HERMITE_QUADRATURE,&
    & BASIS_ADAPTIVE_GAUSS_LEGENDRE_QUADRATURE,BASIS_GAUSS_SIMPLEX_QUADRATURE

  PUBLIC BASIS_XI_COLLAPSED,BASIS_COLLAPSED_AT_XI0,BASIS_COLLAPSED_AT_XI1,BASIS_NOT_COLLAPSED

  PUBLIC BASIS_NO_BOUNDARY_XI,BASIS_LINE_BOUNDARY_XI,BASIS_FACE_BOUNDARY_XI

  PUBLIC Basis_AreaToXiCoordinates

  PUBLIC Basis_BoundaryXiToXi,Basis_XiToBoundaryXi
  
  PUBLIC BASIS_COLLAPSED_XI_SET

  PUBLIC Basis_CollapsedXiSet

  PUBLIC BASIS_EVALUATE_XI

  PUBLIC Basis_EvaluateXi

  PUBLIC Basis_GaussPointsCalculate
  
  PUBLIC BASIS_INTERPOLATE_GAUSS,BASIS_INTERPOLATE_XI,BASIS_INTERPOLATE_LOCAL_FACE_GAUSS

  PUBLIC Basis_InterpolateGauss,Basis_InterpolateXi,Basis_InterpolateLocalFaceGauss

  PUBLIC BASIS_LOCAL_NODE_XI_CALCULATE

  PUBLIC Basis_LocalNodeXiCalculate

  PUBLIC BASIS_NUMBER_OF_LOCAL_NODES_GET

  PUBLIC Basis_NumberOfLocalNodesGet

  PUBLIC BASIS_INTERPOLATION_XI_SET

  PUBLIC Basis_InterpolationXiSet

  PUBLIC BASIS_NUMBER_OF_XI_SET

  PUBLIC Basis_NumberOfXiSet

  PUBLIC BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET

  PUBLIC Basis_QuadratureNumberOfGaussXiSet

  PUBLIC BASIS_QUADRATURE_DESTROY,BASIS_QUADRATURE_ORDER_SET,BASIS_QUADRATURE_TYPE_SET

  PUBLIC Basis_QuadratureDestroy,Basis_QuadratureOrderSet,Basis_QuadratureTypeSet

  PUBLIC BASIS_TYPE_SET

  PUBLIC Basis_TypeSet

  PUBLIC Basis_QuadratureLocalFaceGaussEvaluateSet
  
  PUBLIC BASIS_CREATE_START,BASIS_CREATE_FINISH

  PUBLIC Basis_CreateStart,Basis_CreateFinish

  PUBLIC Basis_Destroy

  PUBLIC Bases_Finalise,Bases_Initialise

  PUBLIC BASIS_COLLAPSED_XI_GET

  PUBLIC Basis_CollapsedXiGet

  PUBLIC BASIS_INTERPOLATION_XI_GET

  PUBLIC Basis_InterpolationXiGet

  PUBLIC BASIS_NUMBER_OF_XI_GET

  PUBLIC Basis_NumberOfXiGet

  PUBLIC BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_GET

  PUBLIC Basis_QuadratureNumberOfGaussXiGet

!!TODO: These should be changed to operate on an array rather than have separate single/multiple routines
  PUBLIC BASIS_QUADRATURE_SINGLE_GAUSS_XI_GET, BASIS_QUADRATURE_MULTIPLE_GAUSS_XI_GET

  PUBLIC Basis_QuadratureSingleGaussXiGet,Basis_QuadratureMultipleGaussXiGet

  PUBLIC BASIS_QUADRATURE_ORDER_GET

  PUBLIC Basis_QuadratureOrderGet

  PUBLIC BASIS_QUADRATURE_TYPE_GET

  PUBLIC Basis_QuadratureTypeGet

  PUBLIC BASIS_TYPE_GET

  PUBLIC Basis_TypeGet

  PUBLIC Basis_XiToAreaCoordinates

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finalises the bases and deallocates all memory
  SUBROUTINE Bases_Finalise(err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Bases_Finalise",err,error,*999)

    !Destroy any created basis functions
    DO WHILE(basisFunctions%numberOfBasisFunctions>0)
      CALL BASIS_DESTROY(basisFunctions%bases(1)%ptr,err,error,*999)
    ENDDO !nb
    !Destroy basis functions and deallocated any memory allocated
    basisFunctions%numberOfBasisFunctions=0
    IF(ALLOCATED(basisFunctions%bases)) DEALLOCATE(basisFunctions%bases)
    
    EXITS("Bases_Finalise")
    RETURN
999 ERRORSEXITS("Bases_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Bases_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the bases.
  SUBROUTINE Bases_Initialise(err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Bases_Initialise",err,error,*999)

    basisFunctions%numberOfBasisFunctions=0    
    
    EXITS("Bases_Initialise")
    RETURN
999 ERRORSEXITS("Bases_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Bases_Initialise

  !
  !================================================================================================================================
  !

  !>Converts area coordinates to xi coordinates. \see BASIS_ROUTINES::Basis_XiToAreaCoordinates
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

  !>Finishes the creation of a new basis \see BASIS_ROUTINES::Basis_CreateStart,OpenCMISS::BasisCreateFinish
  SUBROUTINE Basis_CreateFinish(basis,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: basis !<A pointer to the basis to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: xiIdx,xiCoordIdx,localNodeIdx,localNodeIdx1,localNodeIdx2,localNodeIdx3,localNodeIdx4,elemParamIdx, &
      & localLineIdx,localFaceIdx,columnStart,columnStart2,columnStop,columnStop2
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Basis_CreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(basis%BASIS_FINISHED) CALL FlagError("Basis has already been finished.",err,error,*999)
    
    SELECT CASE(basis%type)
    CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
      CALL Basis_LHTPFamilyCreate(basis,err,error,*999)
    CASE(BASIS_SIMPLEX_TYPE)
      CALL Basis_SimplexFamilyCreate(basis,err,error,*999)
    CASE(BASIS_RADIAL_TYPE)
      CALL BASIS_RADIAL_FAMILY_CREATE(basis,err,error,*999)
    CASE DEFAULT
      localError="Basis type "//TRIM(NumberToVString(BASIS%TYPE,"*",err,error))//" is invalid or not implemented"
      CALL FlagError(localError,err,error,*999)
    END SELECT
    basis%BASIS_FINISHED=.TRUE.
     
    IF(diagnostics1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Basis user number = ",basis%USER_NUMBER,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Basis family number = ",basis%FAMILY_NUMBER,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Basis global number = ",basis%GLOBAL_NUMBER,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Basis type = ",basis%type,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Basis degenerate = ",basis%degenerate,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of Xi directions = ",basis%NUMBER_OF_XI,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of Xi coordinates = ",basis%NUMBER_OF_XI_COORDINATES,err,error,*999)
      
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%NUMBER_OF_XI_COORDINATES,4,4,basis%INTERPOLATION_TYPE, &
        & '("  Interpolation type(:)  :",4(X,I2))','(25X,4(X,I2))',err,error,*999)      
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%NUMBER_OF_XI_COORDINATES,4,4,basis%INTERPOLATION_ORDER, &
        & '("  Interpolation order(:) :",4(X,I2))','(26X,4(X,I2))',err,error,*999)      
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%NUMBER_OF_XI,3,3,basis%COLLAPSED_XI, &
        & '("  Collapsed xi(:) :",3(X,I2))','(26X,3(X,I2))',err,error,*999)      
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of partial derivatives = ",basis%NUMBER_OF_PARTIAL_DERIVATIVES, &
        & err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of nodes = ",basis%NUMBER_OF_NODES,err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%NUMBER_OF_XI_COORDINATES,4,4,basis%NUMBER_OF_NODES_XIC, &
        & '("  Number of nodes(:)  :",4(X,I2))','(22X,4(X,I2))',err,error,*999)      
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%NUMBER_OF_NODES,8,8,basis%NUMBER_OF_DERIVATIVES, &
        & '("  Number of derivatives(:) :",8(X,I2))','(28X,8(X,I2))',err,error,*999)      
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of element parameters = ",basis%NUMBER_OF_ELEMENT_PARAMETERS, &
        & err,error,*999)
! CPB 23/07/07 Doxygen may or may not like this line for some reason????      
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%NUMBER_OF_NODES,8,8,basis%NODE_AT_COLLAPSE, &
        & '("  Node at collapse(:) :",8(X,L1))','(23X,8(X,L1))',err,error,*999)      
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Node position index:",err,error,*999)
      DO xiCoordIdx=1,basis%NUMBER_OF_XI_COORDINATES
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Xi coordinate = ",xiCoordIdx,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%NUMBER_OF_NODES,16,16,basis%NODE_POSITION_INDEX(:,xiCoordIdx), &
          & '("      Index(:)   :",16(X,I2))','(18X,16(X,I2))',err,error,*999)        
      ENDDO !xiCoordIdx
      
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Inverse node position index:",err,error,*999)
      SELECT CASE(basis%NUMBER_OF_XI_COORDINATES)
      CASE(1)
        DO localNodeIdx1=1,basis%NUMBER_OF_NODES_XIC(1)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Xi coordinate = 1, Node position index = ",localNodeIdx1, &
            & err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Index = ",basis%NODE_POSITION_INDEX_INV(localNodeIdx1,1,1,1), &
            & err,error,*999)          
        ENDDO !localNodeIdx1
      CASE(2)
        DO localNodeIdx2=1,basis%NUMBER_OF_NODES_XIC(2)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Xi coordinate = 2, Node position index = ",localNodeIdx2, &
            & err,error,*999)
          DO localNodeIdx1=1,basis%NUMBER_OF_NODES_XIC(1)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Xi coordinate = 1, Node position index = ",localNodeIdx1, &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Index = ",basis%NODE_POSITION_INDEX_INV(localNodeIdx1, &
              & localNodeIdx2,1,1),err,error,*999)
          ENDDO !localNodeIdx1
        ENDDO !localNodeIdx2
      CASE(3)
        DO localNodeIdx3=1,basis%NUMBER_OF_NODES_XIC(3)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Xi coordinate = 3, Node position index = ",localNodeIdx3, &
            & err,error,*999)
          DO localNodeIdx2=1,basis%NUMBER_OF_NODES_XIC(2)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Xi coordinate = 2, Node position index = ",localNodeIdx2, &
              & err,error,*999)
            DO localNodeIdx1=1,basis%NUMBER_OF_NODES_XIC(1)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Xi coordinate = 1, Node position index = ",localNodeIdx1, &
                & err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Index = ",basis%NODE_POSITION_INDEX_INV(localNodeIdx1, &
                & localNodeIdx2,localNodeIdx3,1),err,error,*999)
            ENDDO !localNodeIdx1
          ENDDO !localNodeIdx2
        ENDDO !localNodeIdx3
      CASE(4)
        DO localNodeIdx4=1,basis%NUMBER_OF_NODES_XIC(4)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Xi coordinate = 4, Node position index = ",localNodeIdx4, &
            & err,error,*999)
          DO localNodeIdx3=1,basis%NUMBER_OF_NODES_XIC(3)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Xi coordinate = 3, Node position index = ",localNodeIdx3, &
              & err,error,*999)
            DO localNodeIdx2=1,basis%NUMBER_OF_NODES_XIC(2)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Xi coordinate = 2, Node position index = ",localNodeIdx2, &
                & err,error,*999)
              DO localNodeIdx1=1,basis%NUMBER_OF_NODES_XIC(1) 
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Xi coordinate = 1, Node position index = ",localNodeIdx1, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Index = " &
                  & ,basis%NODE_POSITION_INDEX_INV(localNodeIdx1,localNodeIdx2,localNodeIdx3,localNodeIdx4),err,error,*999)
              ENDDO !localNodeIdx1
            ENDDO !localNodeIdx2
          ENDDO !localNodeIdx3
        ENDDO !localNodeIdx4
      CASE DEFAULT
        CALL FlagError("Invalid number of xi coordinates",err,error,*999)
      END SELECT
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Derivative order index:",err,error,*999)
      DO xiIdx=1,basis%NUMBER_OF_XI
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Xi = ",xiIdx,err,error,*999)
        DO localNodeIdx=1,basis%NUMBER_OF_NODES
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Node = ",localNodeIdx,err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%NUMBER_OF_DERIVATIVES(localNodeIdx),8,8, &
            & basis%DERIVATIVE_ORDER_INDEX(:,localNodeIdx,xiIdx),'("        Index(:) :",8(X,I2))','(18X,8(X,I2))',err,error,*999)
        ENDDO !localNodeIdx
      ENDDO !xiIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Inverse derivative order index:",err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Element parameter index:",err,error,*999)
      DO localNodeIdx=1,basis%NUMBER_OF_NODES
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Node = ",localNodeIdx,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%NUMBER_OF_DERIVATIVES(localNodeIdx),8,8, &
          & basis%ELEMENT_PARAMETER_INDEX(:,localNodeIdx),'("      Index(:)   :",8(X,I2))','(18X,8(X,I2))',err,error,*999)
      ENDDO !localNodeIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Inverse element parameter index:",err,error,*999)
      DO elemParamIdx=1,basis%NUMBER_OF_ELEMENT_PARAMETERS
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Element parameter = ",elemParamIdx,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,2,2, &
          & basis%ELEMENT_PARAMETER_INDEX_INV(:,elemParamIdx),'("      Index(:)  :",2(X,I2))','(18X,2(X,I2))',err,error,*999)
      ENDDO !elemParamIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Partial derivative index:",err,error,*999)
      DO localNodeIdx=1,basis%NUMBER_OF_NODES
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Node = ",localNodeIdx,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%NUMBER_OF_DERIVATIVES(localNodeIdx),8,8, &
          & basis%PARTIAL_DERIVATIVE_INDEX(:,localNodeIdx),'("      Index(:)   :",8(X,I2))','(18X,8(X,I2))',err,error,*999)
      ENDDO !localNodeIdx
      IF(basis%NUMBER_OF_XI==3) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Local faces:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of local faces = ",basis%NUMBER_OF_LOCAL_FACES,err,error,*999)
        DO localFaceIdx=1,basis%NUMBER_OF_LOCAL_FACES
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Local face = ",localFaceIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of nodes in local face = ", &
            & basis%NUMBER_OF_NODES_IN_LOCAL_FACE(localFaceIdx),err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%NUMBER_OF_NODES_IN_LOCAL_FACE(localFaceIdx),4,4, &
            & basis%NODE_NUMBERS_IN_LOCAL_FACE(:,localFaceIdx),'("      Nodes in local face       :",4(X,I2))','(33X,4(X,I2))', &
            & err,error,*999)
          DO localNodeIdx=1,basis%NUMBER_OF_NODES_IN_LOCAL_FACE(localFaceIdx)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Local face node: ",localNodeIdx,err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(0,localNodeIdx, &
              & localFaceIdx),4,4,basis%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(1:,localNodeIdx,localFaceIdx), &
              & '("      Derivatives in local face :",4(X,I2))','(33X,4(X,I2))',err,error,*999)
          ENDDO !localNodeIdx
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%NUMBER_OF_XI_COORDINATES-1,3,3, &
            & basis%localFaceXiDirections(:,localFaceIdx),'("      Local face xi directions  :",3(X,I2))','(33X,3(X,I2))', &
            & err,error,*999)
          CALL WriteStringFmtValue(DIAGNOSTIC_OUTPUT_TYPE,"      Local face xi normal      : ", &
            & basis%localFaceXiNormal(localFaceIdx),"I2",err,error,*999)
        ENDDO !localFaceIdx
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,2*basis%NUMBER_OF_XI_COORDINATES+1,9,9, &
          & basis%xiNormalLocalFace(-basis%NUMBER_OF_XI_COORDINATES:basis%NUMBER_OF_XI_COORDINATES), &
          & '("    Xi normal local face :",9(X,I2))','(26X,9(X,I2))',err,error,*999)
      ENDIF
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Local lines:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of local lines = ",basis%NUMBER_OF_LOCAL_LINES,err,error,*999)
      DO localLineIdx=1,basis%NUMBER_OF_LOCAL_LINES
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Local line = ",localLineIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of nodes in local line = ", &
          & basis%NUMBER_OF_NODES_IN_LOCAL_LINE(localLineIdx),err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%NUMBER_OF_NODES_IN_LOCAL_LINE(localLineIdx),4,4, &
          & basis%NODE_NUMBERS_IN_LOCAL_LINE(:,localLineIdx),'("      Nodes in local line       :",4(X,I2))','(33X,4(X,I2))', &
          & err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%NUMBER_OF_NODES_IN_LOCAL_LINE(localLineIdx),4,4, &
          & basis%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(:,localLineIdx),'("      Derivatives in local line :",4(X,I2))', &
          & '(33X,4(X,I2))',err,error,*999)
        CALL WriteStringFmtValue(DIAGNOSTIC_OUTPUT_TYPE,"      Local line xi direction   : ", &
          & basis%localLineXiDirection(localLineIdx),"I2",err,error,*999)
        IF(basis%NUMBER_OF_XI>1) THEN
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%NUMBER_OF_XI-1,2,2,basis%localLineXiNormals(:,localLineIdx), &
            &  '("      Local line xi normals     :",2(X,I2))','(31X,2(X,I2))',err,error,*999)
        ENDIF
      ENDDO !localLineIdx
      IF(basis%NUMBER_OF_XI>=2) THEN
        IF(basis%NUMBER_OF_XI==3) THEN
          columnStart=1
          columnStart2=-basis%NUMBER_OF_XI_COORDINATES
          columnStop=2*basis%NUMBER_OF_XI_COORDINATES+1
          columnStop2=basis%NUMBER_OF_XI_COORDINATES
        ELSE
          columnStart=1
          columnStart2=1
          columnStop=1
          columnStop2=1
        ENDIF
        CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,2*basis%NUMBER_OF_XI_COORDINATES+1, &
          & columnStart,1,columnStop,9,9,basis%xiNormalsLocalLine(-basis%NUMBER_OF_XI_COORDINATES:basis%NUMBER_OF_XI_COORDINATES, &
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

  !>Starts the creation of a new basis 
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
  SUBROUTINE Basis_CreateStart(userNumber,basis,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the basis to start the creation of
    TYPE(BASIS_TYPE), POINTER :: basis !<On return, A pointer to the created basis. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: basisIdx
    TYPE(BASIS_TYPE), POINTER :: newBasis
    TYPE(BASIS_PTR_TYPE), ALLOCATABLE :: newBases(:)
    TYPE(VARYING_STRING) :: localError

    ENTERS("Basis_CreateStart",err,error,*999)

    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated",err,error,*999)
    
    !See if basis number has already been created
    CALL Basis_UserNumberFind(userNumber,basis,err,error,*999)
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
    DO basisIdx=1,basisFunctions%numberOfBasisFunctions
      newBases(basisIdx)%ptr=>basisFunctions%bases(basisIdx)%ptr
    ENDDO !basisIdx
    basisFunctions%numberOfBasisFunctions=basisFunctions%numberOfBasisFunctions+1
    newBases(basisFunctions%numberOfBasisFunctions)%ptr=>newBasis
    CALL MOVE_ALLOC(newBases,basisFunctions%bases)
    !Set the basis parameters
    newBasis%USER_NUMBER=userNumber
    newBasis%FAMILY_NUMBER=0
    newBasis%GLOBAL_NUMBER=basisFunctions%numberOfBasisFunctions
    !Set the default basis parameters
    newBasis%type=BASIS_LAGRANGE_HERMITE_TP_TYPE
    newBasis%NUMBER_OF_XI=3
    ALLOCATE(newBasis%INTERPOLATION_XI(3),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate basis interpolation xi.",err,error,*999)
    newBasis%INTERPOLATION_XI(1:3)=[BASIS_LINEAR_LAGRANGE_INTERPOLATION,BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
      & BASIS_LINEAR_LAGRANGE_INTERPOLATION]
    ALLOCATE(newBasis%COLLAPSED_XI(3),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate basis collapsed xi.",err,error,*999)
    newBasis%COLLAPSED_XI(1:3)=BASIS_NOT_COLLAPSED
    !Initialise the basis quadrature
    NULLIFY(newBasis%QUADRATURE%basis)
    CALL BASIS_QUADRATURE_INITIALISE(newBasis,err,error,*999)        
    basis=>newBasis
    
    EXITS("Basis_CreateStart")
    RETURN
999 IF(ASSOCIATED(newBasis)) CALL BASIS_DESTROY(newBasis,err,error,*998)
998 IF(ALLOCATED(newBases)) DEALLOCATE(newBases)
    ERRORSEXITS("Basis_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_CreateStart

  !
  !================================================================================================================================
  !

  !>Destroys a basis identified by its basis user number \see BASIS_ROUTINES::BASIS_DESTROY_FAMILY,OpenCMISS::Iron::cmfe_BasisDestroy
  RECURSIVE SUBROUTINE BASIS_DESTROY_NUMBER(USER_NUMBER,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the basis to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("BASIS_DESTROY_NUMBER",err,error,*999)

    CALL Basis_FamilyDestroy(USER_NUMBER,0,err,error,*999)
    
    EXITS("BASIS_DESTROY_NUMBER")
    RETURN
999 ERRORSEXITS("BASIS_DESTROY_NUMBER",err,error)
    RETURN 1
    
  END SUBROUTINE BASIS_DESTROY_NUMBER

  !
  !================================================================================================================================
  !

  !>Destroys a basis. \see BASIS_ROUTINES::BASIS_DESTROY_FAMILY,OpenCMISS::Iron::cmfe_BasisDestroy
  RECURSIVE SUBROUTINE BASIS_DESTROY(BASIS,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: USER_NUMBER
        
    ENTERS("BASIS_DESTROY",err,error,*999)

    IF(ASSOCIATED(BASIS)) THEN
      USER_NUMBER=BASIS%USER_NUMBER
      CALL Basis_FamilyDestroy(USER_NUMBER,0,err,error,*999)
      !NULLIFY(BASIS)
    ELSE
      CALL FlagError("Basis is not associated.",err,error,*999)
    ENDIF
    
    EXITS("BASIS_DESTROY")
    RETURN
999 ERRORSEXITS("BASIS_DESTROY",err,error)
    RETURN 1
    
  END SUBROUTINE BASIS_DESTROY

  !
  !================================================================================================================================
  !

  !>Evaluates the appropriate partial derivative index at position XI for the basis for double precision arguments.
  !>Note for simplex basis functions the XI coordinates should exclude the last area coordinate.
  FUNCTION BASIS_EVALUATE_XI_DP(BASIS,ELEMENT_PARAMETER_INDEX,PARTIAL_DERIV_INDEX,XI,err,error)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: ELEMENT_PARAMETER_INDEX !<The element parameter index to evaluate i.e., the local basis index within the element basis.
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIV_INDEX !<The partial derivative index to evaluate \see CONSTANTS_PartialDerivativeConstants
    REAL(DP), INTENT(IN) :: XI(:) !<The Xi position to evaluate the basis function at
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: BASIS_EVALUATE_XI_DP
    !Local Variables
    INTEGER(INTG) :: nn,nk
    REAL(DP) :: XIL(SIZE(XI,1)+1)
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("BASIS_EVALUATE_XI_DP",err,error,*999)
    
    BASIS_EVALUATE_XI_DP=0.0_DP
    IF(ASSOCIATED(BASIS)) THEN
      IF(ELEMENT_PARAMETER_INDEX>0.AND.ELEMENT_PARAMETER_INDEX<=BASIS%NUMBER_OF_ELEMENT_PARAMETERS) THEN
        SELECT CASE(BASIS%TYPE)
        CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
          nn=BASIS%ELEMENT_PARAMETER_INDEX_INV(1,ELEMENT_PARAMETER_INDEX)
          nk=BASIS%ELEMENT_PARAMETER_INDEX_INV(2,ELEMENT_PARAMETER_INDEX)
          BASIS_EVALUATE_XI_DP=BASIS_LHTP_BASIS_EVALUATE(BASIS,nn,nk,PARTIAL_DERIV_INDEX,XI,err,error)
          IF(ERR/=0) GOTO 999
        CASE(BASIS_SIMPLEX_TYPE)
          !Create the area coordinates from the xi coordinates
          CALL Basis_XiToAreaCoordinates(XI(1:SIZE(XI,1)),XIL(1:SIZE(XI,1)+1),err,error,*999)
          nn=BASIS%ELEMENT_PARAMETER_INDEX_INV(1,ELEMENT_PARAMETER_INDEX)
          BASIS_EVALUATE_XI_DP=BASIS_SIMPLEX_BASIS_EVALUATE(BASIS,nn,PARTIAL_DERIV_INDEX,XIL,err,error)
          IF(ERR/=0) GOTO 999
        CASE(BASIS_RADIAL_TYPE)


          IF(ERR/=0) GOTO 999
        CASE DEFAULT
          LOCAL_ERROR="Basis type "//TRIM(NumberToVString(BASIS%TYPE,"*",err,error))//" is invalid or not implemented."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The specified element parameter index of "// &
          & TRIM(NumberToVString(ELEMENT_PARAMETER_INDEX,"*",err,error))// &
          & " is invalid. The index must be > 0 and <= "// &
          & TRIM(NumberToVString(BASIS%NUMBER_OF_ELEMENT_PARAMETERS,"*",err,error))//"."
        CALL FlagError(LOCAL_ERROR,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated.",err,error,*999)
    ENDIF
    
    EXITS("BASIS_EVALUATE_XI_DP")
    RETURN
999 ERRORSEXITS("BASIS_EVALUATE_XI_DP",err,error)
    RETURN
    
  END FUNCTION BASIS_EVALUATE_XI_DP
  
  !
  !================================================================================================================================
  !

  !>Converts a xi location on a boundary face or line to a full xi location.
  SUBROUTINE Basis_BoundaryXiToXi(basis,localLineFaceNumber,boundaryXi,fullXi,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: basis !<A pointer to the basis to convert the boundary xi for
    INTEGER(INTG), INTENT(IN) :: localLineFaceNumber !<The local line/face number containing the boundary xi
    REAL(DP), INTENT(IN) :: boundaryXi(:) !<The boundary xi location to convert to the full xi location. Note the size of the boundary xi array will determine if we are dealing with a line xi or a face xi.
    REAL(DP), INTENT(OUT) :: fullXi(:) !<On exit, the equivalent full xi location of the boundary xi location
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: normalXi1,normalXi2,numberOfBoundaryXi,numberOfXi
    TYPE(VARYING_STRING) :: localError
        
    ENTERS("Basis_BoundaryXiToXi",err,error,*999)

    IF(ASSOCIATED(basis)) THEN
      IF(basis%BASIS_FINISHED) THEN
        numberOfBoundaryXi=SIZE(boundaryXi,1)
        numberOfXi=basis%NUMBER_OF_XI
        IF(numberOfBoundaryXi<=numberOfXi) THEN
          IF(SIZE(fullXi,1)>=numberOfXi) THEN
            IF(numberOfBoundaryXi==numberOfXi) THEN
              !Basis is of the same dimension as the boundary so just copy over
              fullXi(1:numberOfXi)=boundaryXi(1:numberOfXi)
            ELSE
              SELECT CASE(numberOfBoundaryXi)
              CASE(1)
                !On a line
                IF(localLineFaceNumber>=1.AND.localLineFaceNumber<=basis%NUMBER_OF_LOCAL_LINES) THEN
                  SELECT CASE(basis%type)
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
                      & " for basis number "//TRIM(NumberToVString(basis%USER_NUMBER,"*",err,error))//" is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                ELSE
                  localError="The specified local line/face number of "// &
                    & TRIM(NumberToVString(localLineFaceNumber,"*",err,error))// &
                    & " is invalid. The local line number must be >=1 and <= "// &
                    & TRIM(NumberToVString(basis%NUMBER_OF_LOCAL_LINES,"*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
              CASE(2)
                !On a face
                IF(localLineFaceNumber>=1.AND.localLineFaceNumber<=basis%NUMBER_OF_LOCAL_FACES) THEN
                  normalXi1=basis%localFaceXiNormal(localLineFaceNumber)
                  SELECT CASE(basis%type)
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
                      & " for basis number "//TRIM(NumberToVString(basis%USER_NUMBER,"*",err,error))//" is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                ELSE
                  localError="The specified local line/face number of "// &
                    & TRIM(NumberToVString(localLineFaceNumber,"*",err,error))// &
                    & " is invalid. The local face number must be >=1 and <= "// &
                    & TRIM(NumberToVString(basis%NUMBER_OF_LOCAL_FACES,"*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
              CASE DEFAULT
                localError="Invalid number of boundary xi directions. The number of boundary xi directions of "// &
                  & TRIM(NumberToVString(numberOfBoundaryXi,"*",err,error))//" must be <= the number of basis xi directions of "// &
                  & TRIM(NumberToVString(numberOfXi,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
          ELSE
            localError="The size of the full xi array of "//TRIM(NumberToVString(SIZE(fullXi,1),"*",err,error))// &
              & " must be >= the number of basis xi directions of "//TRIM(NumberToVString(numberOfXi,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          localError="The size of the boundary xi array of "//TRIM(NumberToVString(numberOfBoundaryXi,"*",err,error))// &
            & " is invalid. The size must be >= 1 and <= 2."
          CALL FlagError(localError,err,error,*999)          
        ENDIF
      ELSE
        CALL FlagError("Basis has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated.",err,error,*999)
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
    TYPE(BASIS_TYPE), POINTER :: basis !<A pointer to the basis to convert the full xi for
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
    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(.NOT.basis%BASIS_FINISHED) CALL FlagError("Basis has not been finished.",err,error,*999)
    numberOfXi=basis%NUMBER_OF_XI
    numberOfXiCoordinates=basis%NUMBER_OF_XI_COORDINATES
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
      SELECT CASE(basis%TYPE)
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
        localError="The basis type of "//TRIM(NumberToVString(basis%TYPE,"*",err,error))// &
          & " for basis number "//TRIM(NumberToVString(basis%USER_NUMBER,"*",err,error))//" is invalid."
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
      SELECT CASE(basis%TYPE)
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
        localError="The basis type of "//TRIM(NumberToVString(basis%TYPE,"*",err,error))// &
          & " for basis number "//TRIM(NumberToVString(basis%USER_NUMBER,"*",err,error))//" is invalid."
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
    TYPE(BASIS_TYPE), POINTER :: basis !<A pointer to the basis
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

    SELECT CASE(basis%TYPE)
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
  
  !>Destroys a basis identified by its basis user number and family number. Called from the library visible routine BASIS_DESTROY
  !> \see BASIS_ROUTINES::BASIS_DESTROY
  RECURSIVE SUBROUTINE Basis_FamilyDestroy(userNumber,familyNumber,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the basis to destroy
    INTEGER(INTG), INTENT(IN) :: familyNumber !<The family number of the basis to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: count,basisIdx
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(BASIS_PTR_TYPE), ALLOCATABLE :: newSubBases(:)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Basis_FamilyDestroy",err,error,*999)

    NULLIFY(basis)
    CALL Basis_FamilyNumberFind(userNumber,familyNumber,basis,err,error,*999)
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
            IF(basis%parentBasis%subBases(basisIdx)%ptr%USER_NUMBER==basis%USER_NUMBER.AND. &
              & basis%parentBasis%subBases(basisIdx)%ptr%FAMILY_NUMBER/=basis%FAMILY_NUMBER) THEN
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
            IF(basisFunctions%bases(basisIdx)%ptr%USER_NUMBER/=basis%USER_NUMBER.AND. &
              & basisFunctions%bases(basisIdx)%ptr%FAMILY_NUMBER==0) THEN
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
        CALL Basis_FamilyDestroy(basis%subBases(1)%ptr%USER_NUMBER,basis%subBases(1)%ptr%FAMILY_NUMBER,err,error,*999)
      ENDDO
      !Now delete this instance
      CALL Basis_FamilyDestroy(userNumber,familyNumber,err,error,*999)
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
    TYPE(BASIS_TYPE), POINTER :: basis !<A pointer to the basis to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Basis_Finalise",err,error,*999)

    IF(ASSOCIATED(basis)) THEN
      IF(ALLOCATED(basis%INTERPOLATION_XI)) DEALLOCATE(basis%INTERPOLATION_XI)
      IF(ALLOCATED(basis%INTERPOLATION_TYPE)) DEALLOCATE(basis%INTERPOLATION_TYPE)
      IF(ALLOCATED(basis%INTERPOLATION_ORDER)) DEALLOCATE(basis%INTERPOLATION_ORDER)
      IF(ALLOCATED(basis%COLLAPSED_XI)) DEALLOCATE(basis%COLLAPSED_XI)
      IF(ALLOCATED(basis%NODE_AT_COLLAPSE)) DEALLOCATE(basis%NODE_AT_COLLAPSE)
      CALL BASIS_QUADRATURE_FINALISE(basis,err,error,*999)
      IF(ALLOCATED(basis%NUMBER_OF_NODES_XIC)) DEALLOCATE(basis%NUMBER_OF_NODES_XIC)
      IF(ALLOCATED(basis%NUMBER_OF_DERIVATIVES)) DEALLOCATE(basis%NUMBER_OF_DERIVATIVES)
      IF(ALLOCATED(basis%NODE_POSITION_INDEX)) DEALLOCATE(basis%NODE_POSITION_INDEX)
      IF(ALLOCATED(basis%NODE_POSITION_INDEX_INV)) DEALLOCATE(basis%NODE_POSITION_INDEX_INV)
      IF(ALLOCATED(basis%DERIVATIVE_ORDER_INDEX)) DEALLOCATE(basis%DERIVATIVE_ORDER_INDEX)
      IF(ALLOCATED(basis%DERIVATIVE_ORDER_INDEX_INV)) DEALLOCATE(basis%DERIVATIVE_ORDER_INDEX_INV)
      IF(ALLOCATED(basis%PARTIAL_DERIVATIVE_INDEX)) DEALLOCATE(basis%PARTIAL_DERIVATIVE_INDEX)
      IF(ALLOCATED(basis%ELEMENT_PARAMETER_INDEX)) DEALLOCATE(basis%ELEMENT_PARAMETER_INDEX)
      IF(ALLOCATED(basis%ELEMENT_PARAMETER_INDEX_INV)) DEALLOCATE(basis%ELEMENT_PARAMETER_INDEX_INV)
      IF(ALLOCATED(basis%localLineBasis)) DEALLOCATE(basis%localLineBasis)
      IF(ALLOCATED(basis%localLineXiDirection)) DEALLOCATE(basis%localLineXiDirection)
      IF(ALLOCATED(basis%localLineXiNormals)) DEALLOCATE(basis%localLineXiNormals)
      IF(ALLOCATED(basis%xiNormalsLocalLine)) DEALLOCATE(basis%xiNormalsLocalLine)
      IF(ALLOCATED(basis%NUMBER_OF_NODES_IN_LOCAL_LINE)) DEALLOCATE(basis%NUMBER_OF_NODES_IN_LOCAL_LINE)
      IF(ALLOCATED(basis%NODE_NUMBERS_IN_LOCAL_LINE)) DEALLOCATE(basis%NODE_NUMBERS_IN_LOCAL_LINE)
      IF(ALLOCATED(basis%DERIVATIVE_NUMBERS_IN_LOCAL_LINE)) DEALLOCATE(basis%DERIVATIVE_NUMBERS_IN_LOCAL_LINE)
      IF(ALLOCATED(basis%ELEMENT_PARAMETERS_IN_LOCAL_LINE)) DEALLOCATE(basis%ELEMENT_PARAMETERS_IN_LOCAL_LINE)
      IF(ALLOCATED(basis%localFaceBasis)) DEALLOCATE(basis%localFaceBasis)
      IF(ALLOCATED(basis%localFaceXiDirections)) DEALLOCATE(basis%localFaceXiDirections)
      IF(ALLOCATED(basis%localFaceXiNormal)) DEALLOCATE(basis%localFaceXiNormal)
      IF(ALLOCATED(basis%xiNormalLocalFace)) DEALLOCATE(basis%xiNormalLocalFace)
      IF(ALLOCATED(basis%NUMBER_OF_NODES_IN_LOCAL_FACE)) DEALLOCATE(basis%NUMBER_OF_NODES_IN_LOCAL_FACE)
      IF(ALLOCATED(basis%NODE_NUMBERS_IN_LOCAL_FACE)) DEALLOCATE(basis%NODE_NUMBERS_IN_LOCAL_FACE)
      IF(ALLOCATED(basis%DERIVATIVE_NUMBERS_IN_LOCAL_FACE)) DEALLOCATE(basis%DERIVATIVE_NUMBERS_IN_LOCAL_FACE)
      IF(ALLOCATED(basis%ELEMENT_PARAMETERS_IN_LOCAL_FACE)) DEALLOCATE(basis%ELEMENT_PARAMETERS_IN_LOCAL_FACE)
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
    TYPE(BASIS_TYPE), POINTER :: basis !<A pointer to the basis to allocate and initialise. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Basis_Initialise",err,error,*998)

    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated.",err,error,*998)

    ALLOCATE(basis,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate basis.",err,error,*999)
    
    basis%USER_NUMBER=0
    basis%GLOBAL_NUMBER=0
    basis%FAMILY_NUMBER=0
    basis%BASIS_FINISHED=.FALSE.
    basis%hermite=.FALSE.
    basis%type=0
    basis%NUMBER_OF_XI=0
    basis%NUMBER_OF_XI_COORDINATES=0
    basis%degenerate=.FALSE.
    basis%NUMBER_OF_COLLAPSED_XI=0
    basis%NUMBER_OF_PARTIAL_DERIVATIVES=0
    basis%NUMBER_OF_NODES=0
    basis%NUMBER_OF_ELEMENT_PARAMETERS=0
    basis%MAXIMUM_NUMBER_OF_DERIVATIVES=0
    basis%NUMBER_OF_LOCAL_LINES=0
    basis%NUMBER_OF_LOCAL_FACES=0
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
  !!>coordinate system with COORDINATE_INTERPOLATE_ADJUST. 
  FUNCTION BASIS_INTERPOLATE_GAUSS_DP(BASIS,PARTIAL_DERIV_INDEX,QUADRATURE_SCHEME,GAUSS_POINT_NUMBER,ELEMENT_PARAMETERS,err,error)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIV_INDEX !<The partial derivative index to interpolate \see CONSTANTS_PartialDerivativeConstants
    INTEGER(INTG), INTENT(IN) :: QUADRATURE_SCHEME !<The quadrature scheme to use \see BASIS_ROUTINE_QuadratureSchemes
    INTEGER(INTG), INTENT(IN) :: GAUSS_POINT_NUMBER !<The Gauss point number in the scheme to interpolte
    REAL(DP), INTENT(IN) :: ELEMENT_PARAMETERS(:) !<The element parameters to interpolate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: BASIS_INTERPOLATE_GAUSS_DP
    !Local Variables
    INTEGER(INTG) :: ns
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: BASIS_QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("BASIS_INTERPOLATE_GAUSS_DP",err,error,*999)
    
    BASIS_INTERPOLATE_GAUSS_DP=0.0_DP
    IF(ASSOCIATED(BASIS)) THEN
      IF(QUADRATURE_SCHEME>0.AND.QUADRATURE_SCHEME<=BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES) THEN
        BASIS_QUADRATURE_SCHEME=>BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(QUADRATURE_SCHEME)%PTR
        IF(ASSOCIATED(BASIS_QUADRATURE_SCHEME)) THEN
          IF(GAUSS_POINT_NUMBER>0.AND.GAUSS_POINT_NUMBER<=BASIS_QUADRATURE_SCHEME%NUMBER_OF_GAUSS) THEN
            IF(PARTIAL_DERIV_INDEX>0.AND.PARTIAL_DERIV_INDEX<=BASIS%NUMBER_OF_PARTIAL_DERIVATIVES) THEN
              DO ns=1,BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                BASIS_INTERPOLATE_GAUSS_DP=BASIS_INTERPOLATE_GAUSS_DP+ &
                  & BASIS_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIV_INDEX,GAUSS_POINT_NUMBER)* &
                  & ELEMENT_PARAMETERS(ns)
              ENDDO !ns
            ELSE
              LOCAL_ERROR="The partial derivative index of "//TRIM(NumberToVString(PARTIAL_DERIV_INDEX,"*",err,error))// &
                & " is invalid. It must be between 1 and "// &
                & TRIM(NumberToVString(BASIS%NUMBER_OF_PARTIAL_DERIVATIVES,"*",err,error))
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            ENDIF
          ENDIF
        ELSE
          CALL FlagError("The quadrature scheme has not been created",err,error,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The quadrature scheme type number of "//TRIM(NumberToVString(QUADRATURE_SCHEME,"*",err,error))// &
          & " is invalid. It must be between 1 and "// &
          & TRIM(NumberToVString(BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES,"*",err,error))
        CALL FlagError(LOCAL_ERROR,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated",err,error,*999)
    ENDIF
    
    EXITS("BASIS_INTERPOLATE_GAUSS_DP")
    RETURN
999 ERRORSEXITS("BASIS_INTERPOLATE_GAUSS_DP",err,error)
    RETURN
  END FUNCTION BASIS_INTERPOLATE_GAUSS_DP

  !
  !================================================================================================================================
  !

  !>Interpolates the appropriate partial derivative index of the element local face parameters at a face gauss point for the basis
  !>for double precision arguments. Note the interpolated value returned needs to be adjusted for the particular
  !!>coordinate system with COORDINATE_INTERPOLATE_ADJUST. 
  FUNCTION BASIS_INTERPOLATE_LOCAL_FACE_GAUSS_DP(BASIS,PARTIAL_DERIV_INDEX,QUADRATURE_SCHEME, &
    & LOCAL_FACE_NUMBER,GAUSS_POINT_NUMBER,FACE_PARAMETERS,err,error)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIV_INDEX !<The partial derivative index to interpolate \see CONSTANTS_PartialDerivativeConstants
    INTEGER(INTG), INTENT(IN) :: QUADRATURE_SCHEME !<The quadrature scheme to use \see BASIS_ROUTINE_QuadratureSchemes
    INTEGER(INTG), INTENT(IN) :: LOCAL_FACE_NUMBER !<The index number of the face to interpolate on
    INTEGER(INTG), INTENT(IN) :: GAUSS_POINT_NUMBER !<The face Gauss point number in the scheme to interpolate
    REAL(DP), INTENT(IN) :: FACE_PARAMETERS(:) !<The face parameters to interpolate (in 3D coordinates)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: BASIS_INTERPOLATE_LOCAL_FACE_GAUSS_DP
    !Local Variables
    INTEGER(INTG) :: ns
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: BASIS_QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("BASIS_INTERPOLATE_LOCAL_FACE_GAUSS_DP",err,error,*999)
    
    BASIS_INTERPOLATE_LOCAL_FACE_GAUSS_DP=0.0_DP
    IF(ASSOCIATED(BASIS)) THEN
      IF(QUADRATURE_SCHEME>0.AND.QUADRATURE_SCHEME<=BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES) THEN
        BASIS_QUADRATURE_SCHEME=>BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(QUADRATURE_SCHEME)%PTR
        IF(ASSOCIATED(BASIS_QUADRATURE_SCHEME)) THEN
          IF(BASIS%QUADRATURE%EVALUATE_FACE_GAUSS) THEN !Alternartively, can check whether scheme's face arrays are allocated?
            IF(LOCAL_FACE_NUMBER>0.AND.LOCAL_FACE_NUMBER<=BASIS%NUMBER_OF_LOCAL_FACES) THEN
              IF(GAUSS_POINT_NUMBER>0.AND.GAUSS_POINT_NUMBER<=BASIS_QUADRATURE_SCHEME%NUMBER_OF_FACE_GAUSS(LOCAL_FACE_NUMBER)) THEN
                IF(PARTIAL_DERIV_INDEX>0.AND.PARTIAL_DERIV_INDEX<=BASIS%NUMBER_OF_PARTIAL_DERIVATIVES) THEN
                  DO ns=1,BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                    BASIS_INTERPOLATE_LOCAL_FACE_GAUSS_DP=BASIS_INTERPOLATE_LOCAL_FACE_GAUSS_DP+ &
                      & BASIS_QUADRATURE_SCHEME%FACE_GAUSS_BASIS_FNS(ns,PARTIAL_DERIV_INDEX,GAUSS_POINT_NUMBER,LOCAL_FACE_NUMBER)* &
                      & FACE_PARAMETERS(ns)
                  ENDDO !ns
                ELSE
                  LOCAL_ERROR="The partial derivative index of "//TRIM(NumberToVString(PARTIAL_DERIV_INDEX,"*",err,error))// &
                    & " is invalid. It must be between 1 and "// &
                    & TRIM(NumberToVString(BASIS%NUMBER_OF_PARTIAL_DERIVATIVES,"*",err,error))
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FlagError("The local face number index is invalid.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("The face gauss interpolation scheme has not been created",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("The quadrature scheme has not been created",err,error,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The quadrature scheme type number of "//TRIM(NumberToVString(QUADRATURE_SCHEME,"*",err,error))// &
          & " is invalid. It must be between 1 and "// &
          & TRIM(NumberToVString(BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES,"*",err,error))
        CALL FlagError(LOCAL_ERROR,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated",err,error,*999)
    ENDIF
    
    EXITS("BASIS_INTERPOLATE_LOCAL_FACE_GAUSS_DP")
    RETURN
999 ERRORSEXITS("BASIS_INTERPOLATE_LOCAL_FACE_GAUSS_DP",err,error)
    RETURN
  END FUNCTION BASIS_INTERPOLATE_LOCAL_FACE_GAUSS_DP

  !
  !================================================================================================================================
  !

  !>Interpolates the appropriate partial derivative index of the element parameters at position XI for the basis
  !>for double precision arguments. Note the interpolated value returned needs to be adjusted for the particular
  !>coordinate system with COORDINATE_INTERPOLATE_ADJUST. Note for simplex basis functions the XI coordinates should
  !>exclude the last area coordinate.
  FUNCTION BASIS_INTERPOLATE_XI_DP(BASIS,PARTIAL_DERIV_INDEX,XI,ELEMENT_PARAMETERS,err,error)
  
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIV_INDEX !<The partial derivative index to interpolate \see CONSTANTS_PartialDerivativeConstants
    REAL(DP), INTENT(IN) :: XI(:) !<The Xi position to interpolate the basis function at
    REAL(DP), INTENT(IN) :: ELEMENT_PARAMETERS(:) !<The element parameters to interpolate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: BASIS_INTERPOLATE_XI_DP
    !Local Variables
    INTEGER(INTG) :: nn,nk,ns
    REAL(DP) :: XIL(SIZE(XI,1)+1)
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("BASIS_INTERPOLATE_XI_DP",err,error,*999)
    
    BASIS_INTERPOLATE_XI_DP=0.0_DP
    IF(ASSOCIATED(BASIS)) THEN
      SELECT CASE(BASIS%TYPE)
      CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
        ns=0
        DO nn=1,BASIS%NUMBER_OF_NODES
          DO nk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
            ns=ns+1
            BASIS_INTERPOLATE_XI_DP=BASIS_INTERPOLATE_XI_DP+ &
              & BASIS_LHTP_BASIS_EVALUATE(BASIS,nn,nk,PARTIAL_DERIV_INDEX,XI,err,error)* &
              & ELEMENT_PARAMETERS(ns)
          ENDDO !nk
        ENDDO !nn
        IF(ERR/=0) GOTO 999
      CASE(BASIS_SIMPLEX_TYPE)
        !Create the area coordinates from the xi coordinates
        CALL Basis_XiToAreaCoordinates(XI(1:SIZE(XI,1)),XIL(1:SIZE(XI,1)+1),err,error,*999)
        ns=0
        DO nn=1,BASIS%NUMBER_OF_NODES
          ns=ns+1
          BASIS_INTERPOLATE_XI_DP=BASIS_INTERPOLATE_XI_DP+ &
            & BASIS_SIMPLEX_BASIS_EVALUATE(BASIS,nn,PARTIAL_DERIV_INDEX,XIL,err,error)* &
            & ELEMENT_PARAMETERS(ns)
        ENDDO !nn
        IF(ERR/=0) GOTO 999
      CASE(BASIS_RADIAL_TYPE)

        IF(ERR/=0) GOTO 999
      CASE DEFAULT
        LOCAL_ERROR="Basis type "//TRIM(NumberToVString(BASIS%TYPE,"*",err,error))//" is invalid or not implemented"
        CALL FlagError(LOCAL_ERROR,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Basis is not associated",err,error,*999)
    ENDIF
    
    EXITS("BASIS_INTERPOLATE_XI_DP")
    RETURN
999 ERRORSEXITS("BASIS_INTERPOLATE_XI_DP",err,error)
    RETURN
  END FUNCTION BASIS_INTERPOLATE_XI_DP
  
  !
  !================================================================================================================================
  !
  
  !>Gets/changes the interpolation type in each xi directions for a basis identified by a pointer.
  SUBROUTINE BASIS_INTERPOLATION_XI_GET(BASIS,INTERPOLATION_XI,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to get the interpolation xi
    INTEGER(INTG), INTENT(OUT) :: INTERPOLATION_XI(:) !<On return, the interpolation xi parameters for each Xi direction \see BASIS_ROUTINES_InterpolationSpecifications
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("BASIS_INTERPOLATION_XI_GET",err,error,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        IF(SIZE(INTERPOLATION_XI,1)>=SIZE(BASIS%INTERPOLATION_XI,1)) THEN
          INTERPOLATION_XI=BASIS%INTERPOLATION_XI
        ELSE
          LOCAL_ERROR="The size of INTERPOLATION_XI is too small. The supplied size is "// &
            & TRIM(NumberToVString(SIZE(INTERPOLATION_XI,1),"*",err,error))//" and it needs to be >= "// &
            & TRIM(NumberToVString(SIZE(BASIS%INTERPOLATION_XI,1),"*",err,error))//"."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Basis has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated.",err,error,*999)
    ENDIF
    
    EXITS("BASIS_INTERPOLATION_XI_GET")
    RETURN
999 ERRORSEXITS("BASIS_INTERPOLATION_XI_GET",err,error)
    RETURN
  END SUBROUTINE BASIS_INTERPOLATION_XI_GET
  

  !
  !================================================================================================================================
  !

  !>Sets/changes the interpolation type in each xi directions where the basis is identified by user number. \see OpenCMISS::Iron::cmfe_BasisInterpolationXiSet
  SUBROUTINE BASIS_INTERPOLATION_XI_SET_NUMBER(USER_NUMBER,INTERPOLATION_XI,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the basis to set the interpolation xi
    INTEGER(INTG), INTENT(IN) :: INTERPOLATION_XI(:) !<The interpolation xi parameters for each Xi direction \see BASIS_ROUTINES_InterpolationSpecifications
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    
    ENTERS("BASIS_INTERPOLATION_XI_SET_NUMBER",err,error,*999)

    CALL Basis_UserNumberFind(USER_NUMBER,BASIS,err,error,*999)
    CALL BASIS_INTERPOLATION_XI_SET(BASIS,INTERPOLATION_XI,err,error,*999)
    
    EXITS("BASIS_INTERPOLATION_XI_SET_NUMBER")
    RETURN
999 ERRORSEXITS("BASIS_INTERPOLATION_XI_SET_NUMBER",err,error)
    RETURN 1
  END SUBROUTINE BASIS_INTERPOLATION_XI_SET_NUMBER

  !
  !================================================================================================================================
  !

  !>Sets/changes the interpolation type in each xi directions for a basis identified by a pointer. \see OpenCMISS::Iron::cmfe_BasisInterpolationXiSet
  SUBROUTINE BASIS_INTERPOLATION_XI_SET_PTR(BASIS,INTERPOLATION_XI,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to set the interpolation xi
    INTEGER(INTG), INTENT(IN) :: INTERPOLATION_XI(:) !<The interpolation xi parameters for each Xi direction \see BASIS_ROUTINES_InterpolationSpecifications
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ni,LAST_INTERP
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("BASIS_INTERPOLATION_XI_SET_PTR",err,error,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        CALL FlagError("Basis has been finished",err,error,*999)
      ELSE
        IF(SIZE(INTERPOLATION_XI,1)==BASIS%NUMBER_OF_XI) THEN
          !Check the input values
          SELECT CASE(BASIS%TYPE)
          CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
            DO ni=1,BASIS%NUMBER_OF_XI
              SELECT CASE(INTERPOLATION_XI(ni))
              CASE(BASIS_LINEAR_LAGRANGE_INTERPOLATION,BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,BASIS_CUBIC_LAGRANGE_INTERPOLATION, &
                & BASIS_CUBIC_HERMITE_INTERPOLATION,BASIS_QUADRATIC1_HERMITE_INTERPOLATION,BASIS_QUADRATIC2_HERMITE_INTERPOLATION)
                !Do nothing
              CASE DEFAULT
                LOCAL_ERROR="Interpolation xi value "//TRIM(NumberToVString(INTERPOLATION_XI(ni),"*",err,error))// &
                  & " for xi direction "//TRIM(NumberToVString(ni,"*",err,error))//" is invalid for a Lagrange-Hermite TP basis."
                CALL FlagError(LOCAL_ERROR,err,error,*999)
              END SELECT
            ENDDO !ni
          CASE(BASIS_SIMPLEX_TYPE)
            LAST_INTERP=INTERPOLATION_XI(1)
            DO ni=1,BASIS%NUMBER_OF_XI
              SELECT CASE(INTERPOLATION_XI(ni))
              CASE(BASIS_LINEAR_SIMPLEX_INTERPOLATION,BASIS_QUADRATIC_SIMPLEX_INTERPOLATION,BASIS_CUBIC_SIMPLEX_INTERPOLATION)
                IF(INTERPOLATION_XI(ni)/=LAST_INTERP) THEN
                  CALL FlagError("The interpolation xi value must be the same for all xi directions for a simplex basis.", &
                    & err,error,*999)
                ENDIF
              CASE DEFAULT
                LOCAL_ERROR="Interpolation xi value "//TRIM(NumberToVString(INTERPOLATION_XI(ni),"*",err,error))// &
                  & " for xi direction "//TRIM(NumberToVString(ni,"*",err,error))//" is invalid for a simplex basis."
                CALL FlagError(LOCAL_ERROR,err,error,*999)
              END SELECT
            ENDDO !ni
          CASE DEFAULT
            CALL FlagError("Invalid basis type or not implemented",err,error,*999)
          END SELECT
          !Set the interpolation xi
          BASIS%INTERPOLATION_XI(1:BASIS%NUMBER_OF_XI)=INTERPOLATION_XI(1:BASIS%NUMBER_OF_XI)
        ELSE
          LOCAL_ERROR="The size of the interpolation xi array ("// &
            & TRIM(NumberToVString(SIZE(INTERPOLATION_XI,1),"*",err,error))//") does not match the number of xi directions ("// &
            & TRIM(NumberToVString(BASIS%NUMBER_OF_XI,"*",err,error))//") for basis number "// &
            & TRIM(NumberToVString(BASIS%USER_NUMBER,"*",err,error))//"."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated.",err,error,*999)
    ENDIF
    
    EXITS("BASIS_INTERPOLATION_XI_SET_PTR")
    RETURN
999 ERRORSEXITS("BASIS_INTERPOLATION_XI_SET_PTR",err,error)
    RETURN 1
  END SUBROUTINE BASIS_INTERPOLATION_XI_SET_PTR

  !
  !================================================================================================================================
  !

  !>Creates and initialises a Lagrange-Hermite tensor product basis that has already been allocated BASIS_CREATE_START
  !> \see BASIS_ROUTINES::BASIS_CREATE_START
  SUBROUTINE Basis_LHTPBasisCreate(basis,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: basis !<A pointer to the basis to create
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
    
    ENTERS("Basis_LHTPBasisCreate",ERR,error,*999)

    IF(ASSOCIATED(basis)) THEN
      IF(basis%type==BASIS_LAGRANGE_HERMITE_TP_TYPE) THEN
        basis%NUMBER_OF_XI_COORDINATES=basis%NUMBER_OF_XI
        basis%NUMBER_OF_PARTIAL_DERIVATIVES=basis%NUMBER_OF_XI_COORDINATES**2+2
        ALLOCATE(basis%INTERPOLATION_TYPE(basis%NUMBER_OF_XI_COORDINATES),STAT=err)
        IF(ERR/=0) CALL FlagError("Could not allocate basis interpolation type.",err,error,*999)
        ALLOCATE(basis%INTERPOLATION_ORDER(basis%NUMBER_OF_XI_COORDINATES),STAT=err)
        IF(ERR/=0) CALL FlagError("Could not allocate basis interpolation order.",err,error,*999)
        ALLOCATE(basis%NUMBER_OF_NODES_XIC(basis%NUMBER_OF_XI_COORDINATES),STAT=err)
        IF(ERR/=0) CALL FlagError("Could not allocate basis number of nodes xic.",err,error,*999)
        numberOfNodes=1
        maximumNumberOfNodes=0
        basis%degenerate=.FALSE.
        basis%NUMBER_OF_COLLAPSED_XI=0
        DO xiIdx=1,basis%NUMBER_OF_XI
          !Set up the interpolation types, orders and number of nodes in each xi from the user specified interpolation xi.
          SELECT CASE(basis%INTERPOLATION_XI(xiIdx))
          CASE(BASIS_LINEAR_LAGRANGE_INTERPOLATION)
            basis%INTERPOLATION_TYPE(xiIdx)=BASIS_LAGRANGE_INTERPOLATION
            basis%INTERPOLATION_ORDER(xiIdx)=BASIS_LINEAR_INTERPOLATION_ORDER
            basis%NUMBER_OF_NODES_XIC(xiIdx)=2            
          CASE(BASIS_QUADRATIC_LAGRANGE_INTERPOLATION)
            basis%INTERPOLATION_TYPE(xiIdx)=BASIS_LAGRANGE_INTERPOLATION
            basis%INTERPOLATION_ORDER(xiIdx)=BASIS_QUADRATIC_INTERPOLATION_ORDER
            basis%NUMBER_OF_NODES_XIC(xiIdx)=3
          CASE(BASIS_CUBIC_LAGRANGE_INTERPOLATION)
            basis%INTERPOLATION_TYPE(xiIdx)=BASIS_LAGRANGE_INTERPOLATION
            basis%INTERPOLATION_ORDER(xiIdx)=BASIS_CUBIC_INTERPOLATION_ORDER
            basis%NUMBER_OF_NODES_XIC(xiIdx)=4
          CASE(BASIS_CUBIC_HERMITE_INTERPOLATION)
            basis%INTERPOLATION_TYPE(xiIdx)=BASIS_HERMITE_INTERPOLATION
            basis%INTERPOLATION_ORDER(xiIdx)=BASIS_CUBIC_INTERPOLATION_ORDER
            basis%NUMBER_OF_NODES_XIC(xiIdx)=2
          CASE(BASIS_QUADRATIC1_HERMITE_INTERPOLATION)
            basis%INTERPOLATION_TYPE(xiIdx)=BASIS_HERMITE_INTERPOLATION
            basis%INTERPOLATION_ORDER(xiIdx)=BASIS_QUADRATIC1_INTERPOLATION_ORDER
            basis%NUMBER_OF_NODES_XIC(xiIdx)=2
          CASE(BASIS_QUADRATIC2_HERMITE_INTERPOLATION)
            basis%INTERPOLATION_TYPE(xiIdx)=BASIS_HERMITE_INTERPOLATION
            basis%INTERPOLATION_ORDER(xiIdx)=BASIS_QUADRATIC2_INTERPOLATION_ORDER
            basis%NUMBER_OF_NODES_XIC(xiIdx)=2
          CASE DEFAULT 
            CALL FlagError("Invalid interpolation type",err,error,*999)
          END SELECT
          IF(basis%COLLAPSED_XI(xiIdx)==BASIS_XI_COLLAPSED) THEN
            basis%NUMBER_OF_COLLAPSED_XI=basis%NUMBER_OF_COLLAPSED_XI+1
            collapsedXi(basis%NUMBER_OF_COLLAPSED_XI)=xiIdx
            basis%degenerate=.TRUE.
          ENDIF
          numberOfNodes=numberOfNodes*basis%NUMBER_OF_NODES_XIC(xiIdx)
          IF(basis%NUMBER_OF_NODES_XIC(xiIdx)>maximumNumberOfNodes) maximumNumberOfNodes=basis%NUMBER_OF_NODES_XIC(xiIdx)
        ENDDO !xiIdx
        !If a degenerate (collapsed) basis recalculate the number of nodes from the maximum posible number of nodes
        IF(basis%degenerate) THEN
          !Calculate the NODE_AT_COLLAPSE array.
          ALLOCATE(nodeAtCollapse(numberOfNodes),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate at collapse",err,error,*999)
          position=1
          basis%NUMBER_OF_NODES=0
          !Loop over the maximum number of nodes which is currently set for the basis
          DO localNodeIdx=1,numberOfNodes
            atCollapse=.FALSE.
            DO xiIdx=1,basis%NUMBER_OF_XI
              IF(basis%COLLAPSED_XI(xiIdx)==BASIS_COLLAPSED_AT_XI0.AND.position(xiIdx)==1.OR. &
                & basis%COLLAPSED_XI(xiIdx)==BASIS_COLLAPSED_AT_XI1.AND.position(xiIdx)==basis%NUMBER_OF_NODES_XIC(xiIdx)) &
                & THEN
                atCollapse=.TRUE.
                firstCollapsedPosition=ALL(position(collapsedXi(1:basis%NUMBER_OF_COLLAPSED_XI))==1)
                EXIT
              ENDIF
            ENDDO !xiIdx
            IF(atCollapse) THEN
              IF(firstCollapsedPosition) THEN
                basis%NUMBER_OF_NODES=basis%NUMBER_OF_NODES+1
                nodeAtCollapse(basis%NUMBER_OF_NODES)=.TRUE.
              ENDIF
            ELSE
              basis%NUMBER_OF_NODES=basis%NUMBER_OF_NODES+1
              nodeAtCollapse(basis%NUMBER_OF_NODES)=.FALSE.
            ENDIF
            position(1)=position(1)+1
            DO xiIdx=1,basis%NUMBER_OF_XI
              IF(position(xiIdx)>basis%NUMBER_OF_NODES_XIC(xiIdx)) THEN
                position(xiIdx)=1
                position(xiIdx+1)=position(xiIdx+1)+1
              ENDIF
            ENDDO !xiIdx
          ENDDO !localNodeIdx
          CALL MOVE_ALLOC(nodeAtCollapse,basis%NODE_AT_COLLAPSE)
        ELSE        
          basis%NUMBER_OF_NODES=numberOfNodes
          ALLOCATE(basis%NODE_AT_COLLAPSE(basis%NUMBER_OF_NODES),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate basis node at collapse.",err,error,*999)
          basis%NODE_AT_COLLAPSE=.FALSE.
          collapsedXi(1)=1
        ENDIF

        ALLOCATE(basis%NODE_POSITION_INDEX(basis%NUMBER_OF_NODES,basis%NUMBER_OF_XI_COORDINATES),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate basis node position index.",err,error,*999)
        SELECT CASE(basis%NUMBER_OF_XI_COORDINATES)
        CASE(1)
          ALLOCATE(basis%NODE_POSITION_INDEX_INV(maximumNumberOfNodes,1,1,1),STAT=err)
        CASE(2)
          ALLOCATE(basis%NODE_POSITION_INDEX_INV(maximumNumberOfNodes,maximumNumberOfNodes,1,1),STAT=err)
        CASE(3)
          ALLOCATE(basis%NODE_POSITION_INDEX_INV(maximumNumberOfNodes,maximumNumberOfNodes,maximumNumberOfNodes,1),STAT=err)
        CASE DEFAULT
          CALL FlagError("Invalid number of xi coordinates.",err,error,*999)
        END SELECT
        IF(ERR/=0) CALL FlagError("Could not allocate node position index inverse.",err,error,*999)
        basis%NODE_POSITION_INDEX_INV=0

        !Determine the node position index and its inverse
        position=1
        collapsedPosition=1
        localNode=0
        firstCollapsedPosition=.TRUE.
        DO localNodeIdx1=1,numberOfNodes
          atCollapse=.FALSE.
          IF(basis%degenerate) THEN
            DO xiIdx=1,basis%NUMBER_OF_XI
              IF(basis%COLLAPSED_XI(xiIdx)==BASIS_COLLAPSED_AT_XI0.AND.position(xiIdx)==1.OR. &
                & basis%COLLAPSED_XI(xiIdx)==BASIS_COLLAPSED_AT_XI1.AND.position(xiIdx)==basis%NUMBER_OF_NODES_XIC(xiIdx)) &
                & THEN 
                atCollapse=.TRUE.
                firstCollapsedPosition=ALL(position(collapsedXi(1:basis%NUMBER_OF_COLLAPSED_XI))==1)
                EXIT
              ENDIF
            ENDDO !xiIdx
          ENDIF
          IF(atCollapse) THEN
            IF(firstCollapsedPosition) THEN
              localNode=localNode+1
              basis%NODE_POSITION_INDEX(localNode,1:basis%NUMBER_OF_XI)=position(1:basis%NUMBER_OF_XI)
              basis%NODE_POSITION_INDEX_INV(position(1),position(2),position(3),1)=localNode
            ELSE
              !The second node in the collapsed xi is set to the same node number as the first node in that xi direction.
              collapsedPosition(1:basis%NUMBER_OF_XI)=position(1:basis%NUMBER_OF_XI)
              collapsedPosition(collapsedXi(1:basis%NUMBER_OF_COLLAPSED_XI))=1
              basis%NODE_POSITION_INDEX_INV(position(1),position(2),position(3),1)= &
                & basis%NODE_POSITION_INDEX_INV(collapsedPosition(1),collapsedPosition(2),collapsedPosition(3),1)
            ENDIF
          ELSE
            localNode=localNode+1
            basis%NODE_POSITION_INDEX(localNode,1:basis%NUMBER_OF_XI)=position(1:basis%NUMBER_OF_XI)
            basis%NODE_POSITION_INDEX_INV(position(1),position(2),position(3),1)=localNode
          ENDIF
          position(1)=position(1)+1
          DO xiIdx=1,basis%NUMBER_OF_XI
            IF(position(xiIdx)>basis%NUMBER_OF_NODES_XIC(xiIdx)) THEN
              position(xiIdx)=1
              position(xiIdx+1)=position(xiIdx+1)+1
            ENDIF
          ENDDO !xiIdx
        ENDDO !localNodeIdx1
        !Calculate the maximum number of derivatives and the number of element parameters
        basis%MAXIMUM_NUMBER_OF_DERIVATIVES=-1
        basis%NUMBER_OF_ELEMENT_PARAMETERS=0
        DO localNodeIdx=1,basis%NUMBER_OF_NODES
          numberOfDerivatives=1
          DO xiIdx=1,basis%NUMBER_OF_XI
            IF((.NOT.basis%NODE_AT_COLLAPSE(localNodeIdx).OR.basis%COLLAPSED_XI(xiIdx)==BASIS_NOT_COLLAPSED).AND. &
              & basis%INTERPOLATION_TYPE(xiIdx)==BASIS_HERMITE_INTERPOLATION.AND. &
              & (basis%INTERPOLATION_ORDER(xiIdx)==BASIS_CUBIC_INTERPOLATION_ORDER.OR. &
              & (basis%NODE_POSITION_INDEX(localNodeIdx,xiIdx)==1.AND. &
              & basis%INTERPOLATION_ORDER(xiIdx)==BASIS_QUADRATIC2_INTERPOLATION_ORDER).OR. &
              & (basis%NODE_POSITION_INDEX(localNodeIdx,xiIdx)==2.AND. &
              & basis%INTERPOLATION_ORDER(xiIdx)==BASIS_QUADRATIC1_INTERPOLATION_ORDER))) THEN
              !Derivative in this direction
              numberOfDerivatives=numberOfDerivatives*2
            ENDIF
          ENDDO !xiIdx
          basis%NUMBER_OF_ELEMENT_PARAMETERS=basis%NUMBER_OF_ELEMENT_PARAMETERS+numberOfDerivatives
          IF(numberOfDerivatives>basis%MAXIMUM_NUMBER_OF_DERIVATIVES) basis%MAXIMUM_NUMBER_OF_DERIVATIVES=numberOfDerivatives
        ENDDO !localNodeIdx
        !Now set up the number of derivatives and derivative order index
        ALLOCATE(basis%NUMBER_OF_DERIVATIVES(basis%NUMBER_OF_NODES),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate number of derivatives.",err,error,*999)
        ALLOCATE(basis%DERIVATIVE_ORDER_INDEX(basis%MAXIMUM_NUMBER_OF_DERIVATIVES,basis%NUMBER_OF_NODES, &
          & basis%NUMBER_OF_XI),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate derivative order index.",err,error,*999)
        ALLOCATE(basis%DERIVATIVE_ORDER_INDEX_INV(FIRST_PART_DERIV,FIRST_PART_DERIV,FIRST_PART_DERIV, &
          & basis%NUMBER_OF_NODES),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate derivative order index inverse.",err,error,*999)
        ALLOCATE(basis%PARTIAL_DERIVATIVE_INDEX(basis%MAXIMUM_NUMBER_OF_DERIVATIVES,basis%NUMBER_OF_NODES),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate partial derivative index.",err,error,*999)
        ALLOCATE(basis%ELEMENT_PARAMETER_INDEX(basis%MAXIMUM_NUMBER_OF_DERIVATIVES,basis%NUMBER_OF_NODES),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate element parameter index.",err,error,*999)
        ALLOCATE(basis%ELEMENT_PARAMETER_INDEX_INV(2,basis%NUMBER_OF_ELEMENT_PARAMETERS),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate element parameter index inverse.",err,error,*999)
        !Set the derivative order index and its inverse, the element parameter index and the partial derivative index.
        elementParameter=0
        basis%DERIVATIVE_ORDER_INDEX=0
        basis%DERIVATIVE_ORDER_INDEX_INV=0
        DO localNodeIdx=1,basis%NUMBER_OF_NODES
          basis%NUMBER_OF_DERIVATIVES(localNodeIdx)=1
          DO xiIdx1=1,basis%NUMBER_OF_XI
            IF((.NOT.basis%NODE_AT_COLLAPSE(localNodeIdx).OR.basis%COLLAPSED_XI(xiIdx1)==BASIS_NOT_COLLAPSED).AND. &
              & basis%INTERPOLATION_TYPE(xiIdx1)==BASIS_HERMITE_INTERPOLATION.AND. &
              & (basis%INTERPOLATION_ORDER(xiIdx1)==BASIS_CUBIC_INTERPOLATION_ORDER.OR. &
              & (basis%NODE_POSITION_INDEX(localNodeIdx,xiIdx1)==1.AND. &
              & basis%INTERPOLATION_ORDER(xiIdx1)==BASIS_QUADRATIC2_INTERPOLATION_ORDER).OR. &
              & (basis%NODE_POSITION_INDEX(localNodeIdx,xiIdx1)==2.AND. &
              & basis%INTERPOLATION_ORDER(xiIdx1)==BASIS_QUADRATIC1_INTERPOLATION_ORDER))) THEN
              oldNumberOfDerivatives=basis%NUMBER_OF_DERIVATIVES(localNodeIdx)
              basis%NUMBER_OF_DERIVATIVES(localNodeIdx)=basis%NUMBER_OF_DERIVATIVES(localNodeIdx)*2
              DO derivativeIdx=1,oldNumberOfDerivatives
                basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,localNodeIdx,xiIdx1)=NO_PART_DERIV
                basis%DERIVATIVE_ORDER_INDEX(oldNumberOfDerivatives+derivativeIdx,localNodeIdx,xiIdx1)=FIRST_PART_DERIV
                DO xiIdx2=1,xiIdx1-1
                  basis%DERIVATIVE_ORDER_INDEX(oldNumberOfDerivatives+derivativeIdx,localNodeIdx,xiIdx2)= &
                    & basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,localNodeIdx,xiIdx2)
                ENDDO !xiIdx2
              ENDDO !derivativeIdx
            ELSE
              DO derivativeIdx=1,basis%NUMBER_OF_DERIVATIVES(localNodeIdx)
                basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,localNodeIdx,xiIdx1)=NO_PART_DERIV
              ENDDO !derivativeIdx
            ENDIF
          ENDDO !xiIdx1

          DO derivativeIdx=1,basis%NUMBER_OF_DERIVATIVES(localNodeIdx)
            elementParameter=elementParameter+1
            basis%ELEMENT_PARAMETER_INDEX(derivativeIdx,localNodeIdx)=elementParameter
            basis%ELEMENT_PARAMETER_INDEX_INV(1,elementParameter)=localNodeIdx
            basis%ELEMENT_PARAMETER_INDEX_INV(2,elementParameter)=derivativeIdx
            SELECT CASE(basis%NUMBER_OF_XI)
            CASE(1)
              basis%DERIVATIVE_ORDER_INDEX_INV(basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,localNodeIdx,1),1,1,localNodeIdx)= &
                & derivativeIdx
              SELECT CASE(basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,localNodeIdx,1))
              CASE(NO_PART_DERIV)
                basis%PARTIAL_DERIVATIVE_INDEX(derivativeIdx,localNodeIdx)=NO_PART_DERIV
              CASE(FIRST_PART_DERIV)
                basis%PARTIAL_DERIVATIVE_INDEX(derivativeIdx,localNodeIdx)=PART_DERIV_S1
              CASE DEFAULT
                CALL FlagError("Invalid derivative order index.",err,error,*999)
              END SELECT
            CASE(2)
              basis%DERIVATIVE_ORDER_INDEX_INV(basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,localNodeIdx,1), &
                & basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,localNodeIdx,2),1,localNodeIdx)=derivativeIdx
              SELECT CASE(basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,localNodeIdx,1))
              CASE(NO_PART_DERIV)
                SELECT CASE(basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,localNodeIdx,2))
                CASE(NO_PART_DERIV)
                  basis%PARTIAL_DERIVATIVE_INDEX(derivativeIdx,localNodeIdx)=NO_PART_DERIV
                CASE(FIRST_PART_DERIV)
                  basis%PARTIAL_DERIVATIVE_INDEX(derivativeIdx,localNodeIdx)=PART_DERIV_S2
                CASE DEFAULT
                  CALL FlagError("Invalid derivative order index.",err,error,*999)
                END SELECT
              CASE(FIRST_PART_DERIV)
                SELECT CASE(basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,localNodeIdx,2))
                CASE(NO_PART_DERIV)
                  basis%PARTIAL_DERIVATIVE_INDEX(derivativeIdx,localNodeIdx)=PART_DERIV_S1
                CASE(FIRST_PART_DERIV)
                  basis%PARTIAL_DERIVATIVE_INDEX(derivativeIdx,localNodeIdx)=PART_DERIV_S1_S2
                CASE DEFAULT
                  CALL FlagError("Invalid derivative order index.",err,error,*999)
                END SELECT
              CASE DEFAULT
                CALL FlagError("Invalid derivative order index.",err,error,*999)
              END SELECT
            CASE(3)
              basis%DERIVATIVE_ORDER_INDEX_INV(basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,localNodeIdx,1), &
                & basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,localNodeIdx,2), &
                & basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,localNodeIdx,3),localNodeIdx)=derivativeIdx
              SELECT CASE(basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,localNodeIdx,1))
              CASE(NO_PART_DERIV)
                SELECT CASE(basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,localNodeIdx,2))
                CASE(NO_PART_DERIV)
                  SELECT CASE(basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,localNodeIdx,3))
                  CASE(NO_PART_DERIV)                
                    basis%PARTIAL_DERIVATIVE_INDEX(derivativeIdx,localNodeIdx)=NO_PART_DERIV
                  CASE(FIRST_PART_DERIV)
                    basis%PARTIAL_DERIVATIVE_INDEX(derivativeIdx,localNodeIdx)=PART_DERIV_S3
                  CASE DEFAULT
                    CALL FlagError("Invalid derivative order index.",err,error,*999)
                  END SELECT
                CASE(FIRST_PART_DERIV)
                  SELECT CASE(basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,localNodeIdx,3))
                  CASE(NO_PART_DERIV)                
                    basis%PARTIAL_DERIVATIVE_INDEX(derivativeIdx,localNodeIdx)=PART_DERIV_S2
                  CASE(FIRST_PART_DERIV)
                    basis%PARTIAL_DERIVATIVE_INDEX(derivativeIdx,localNodeIdx)=PART_DERIV_S2_S3
                  CASE DEFAULT
                    CALL FlagError("Invalid derivative order index.",err,error,*999)
                  END SELECT
                CASE DEFAULT
                  CALL FlagError("Invalid derivative order index.",err,error,*999)
                END SELECT
              CASE(FIRST_PART_DERIV)
                SELECT CASE(basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,localNodeIdx,2))
                CASE(NO_PART_DERIV)
                  SELECT CASE(basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,localNodeIdx,3))
                  CASE(NO_PART_DERIV)                
                    basis%PARTIAL_DERIVATIVE_INDEX(derivativeIdx,localNodeIdx)=PART_DERIV_S1
                  CASE(FIRST_PART_DERIV)
                    basis%PARTIAL_DERIVATIVE_INDEX(derivativeIdx,localNodeIdx)=PART_DERIV_S1_S3
                  CASE DEFAULT
                    CALL FlagError("Invalid derivative order index.",err,error,*999)
                  END SELECT
                CASE(FIRST_PART_DERIV)
                  SELECT CASE(basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,localNodeIdx,3))
                  CASE(NO_PART_DERIV)                
                    basis%PARTIAL_DERIVATIVE_INDEX(derivativeIdx,localNodeIdx)=PART_DERIV_S1_S2
                  CASE(FIRST_PART_DERIV)
                    basis%PARTIAL_DERIVATIVE_INDEX(derivativeIdx,localNodeIdx)=PART_DERIV_S1_S2_S3
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
        SELECT CASE(basis%NUMBER_OF_XI)
        CASE(1) !1 xi directions
          numberOfLocalLines=1
          basis%NUMBER_OF_LOCAL_LINES=1
          ALLOCATE(basis%NUMBER_OF_NODES_IN_LOCAL_LINE(numberOfLocalLines),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate number of nodes in local line.",err,error,*999)
          basis%NUMBER_OF_NODES_IN_LOCAL_LINE(1)=basis%NUMBER_OF_NODES_XIC(1)
          ALLOCATE(basis%localLineXiDirection(numberOfLocalLines),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate local line xi direction.",err,error,*999)
          basis%localLineXiDirection(1)=1
          ALLOCATE(basis%NODE_NUMBERS_IN_LOCAL_LINE(basis%NUMBER_OF_NODES_XIC(1),numberOfLocalLines),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate node numbers in local line.",err,error,*999)
          ALLOCATE(basis%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(basis%NUMBER_OF_NODES_XIC(1),numberOfLocalLines),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate derivative numbers in local line.",err,error,*999)
          basis%DERIVATIVE_NUMBERS_IN_LOCAL_LINE=NO_PART_DERIV
          ALLOCATE(basis%ELEMENT_PARAMETERS_IN_LOCAL_LINE(basis%NUMBER_OF_NODES_XIC(1)**2,numberOfLocalLines),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate element parameters in local line.",err,error,*999)
          basis%ELEMENT_PARAMETERS_IN_LOCAL_LINE=1
          localLineParameter=0
          DO localNodeIdx2=1,basis%NUMBER_OF_NODES_XIC(1)
            DO localNodeIdx1=1,basis%NUMBER_OF_NODES
              IF(basis%NODE_POSITION_INDEX(localNodeIdx1,1)==localNodeIdx2) THEN
                basis%NODE_NUMBERS_IN_LOCAL_LINE(localNodeIdx2,1)=localNodeIdx1
                DO derivativeIdx=1,basis%NUMBER_OF_DERIVATIVES(localNodeIdx2)
                  localLineParameter=localLineParameter+1
                  basis%ELEMENT_PARAMETERS_IN_LOCAL_LINE(localLineParameter,1)=basis%ELEMENT_PARAMETER_INDEX( &
                    & derivativeIdx,localNodeIdx1)
                  IF(basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,localNodeIdx2,1)==FIRST_PART_DERIV) THEN
                    basis%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(localNodeIdx2,1)=derivativeIdx
                    EXIT
                  ENDIF
                ENDDO !derivativeIdx
                EXIT
              ENDIF
            ENDDO !localNodeIdx2
          ENDDO !localNodeIdx1
        CASE(2) !2 xi directions
          !Determine the maximum node extent of the basis
          maximumNodeExtent(1)=MAXVAL(basis%NODE_POSITION_INDEX(:,1))
          maximumNodeExtent(2)=MAXVAL(basis%NODE_POSITION_INDEX(:,2))
          !Allocate and calculate the lines
          numberOfLocalLines=4-basis%NUMBER_OF_COLLAPSED_XI
          ALLOCATE(basis%NUMBER_OF_NODES_IN_LOCAL_LINE(numberOfLocalLines),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate number of nodes in local line.",err,error,*999)
          basis%NUMBER_OF_NODES_IN_LOCAL_LINE=0
          ALLOCATE(basis%localLineXiDirection(numberOfLocalLines),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate local line xi direction.",err,error,*999)
          ALLOCATE(BASIS%localLineXiNormals(1,numberOfLocalLines),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate local line xi normals.",err,error,*999)
          ALLOCATE(basis%xiNormalsLocalLine(-2:2,1),STAT=ERR)
          basis%xiNormalsLocalLine=0
          IF(ERR/=0) CALL FlagError("Could not allocate xi normals local line.",err,error,*999)
          ALLOCATE(basis%NODE_NUMBERS_IN_LOCAL_LINE(MAXVAL(basis%NUMBER_OF_NODES_XIC),numberOfLocalLines),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate node numbers in local line",err,error,*999)
          basis%NODE_NUMBERS_IN_LOCAL_LINE=0
          ALLOCATE(basis%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(MAXVAL(basis%NUMBER_OF_NODES_XIC),numberOfLocalLines),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate derivative numbers in local line.",err,error,*999)
          basis%DERIVATIVE_NUMBERS_IN_LOCAL_LINE=NO_PART_DERIV
          ALLOCATE(basis%ELEMENT_PARAMETERS_IN_LOCAL_LINE(MAXVAL(basis%NUMBER_OF_NODES_XIC)**2,numberOfLocalLines),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate element parameters in local line.",err,error,*999)
          basis%ELEMENT_PARAMETERS_IN_LOCAL_LINE=1
          !Find the lines
          basis%NUMBER_OF_LOCAL_LINES=0
          DO xiIdx1=1,2
            xiIdx2=OTHER_XI_DIRECTIONS2(xiIdx1)
            !We are looking for lines in the xiIdx1 direction from the direction of xiIdx1=0
            !Loop over the element extremes in the xiIdx2 direction
            DO localNodeIdx2=1,maximumNodeExtent(xiIdx2),maximumNodeExtent(xiIdx2)-1
              nodeCount=0
              specialNodeCount=0
              nodesInLine=0
              DO localNodeIdx1=1,basis%NUMBER_OF_NODES
                IF(basis%COLLAPSED_XI(xiIdx1)/=BASIS_NOT_COLLAPSED) THEN
                  !The current xi direction, xiIdx1, is in a degenerate plane
                  IF(basis%COLLAPSED_XI(xiIdx2)==BASIS_XI_COLLAPSED) THEN
                    !The other xi direction is collapsed (must be the case)
                    IF(basis%COLLAPSED_XI(xiIdx1)==BASIS_COLLAPSED_AT_XI0) THEN !Collapsed at the xi=0 end
                      IF(basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx2)==localNodeIdx2.OR. &
                        & basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx1)==1) THEN
                        nodeCount=nodeCount+1
                        nodesInLine(nodeCount)=localNodeIdx1
                      ENDIF
                    ELSE !Collapsed at the xi=1 end
                      IF(basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx2)==localNodeIdx2) THEN
                        nodeCount=nodeCount+1
                        nodesInLine(nodeCount)=localNodeIdx1
                      ELSE IF(basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx1)==maximumNodeExtent(xiIdx1)) THEN
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
                    IF(basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx2)==localNodeIdx2) THEN
                      nodeCount=nodeCount+1
                      nodesInLine(nodeCount)=localNodeIdx1
                    ENDIF
                  ENDIF
                ELSE
                  !The current xi direction, xiIdx1, is not involved in any collapsed (degenerate) planes
                  IF(basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx2)==localNodeIdx2) THEN
                    nodeCount=nodeCount+1
                    nodesInLine(nodeCount)=localNodeIdx1
                  ENDIF
                ENDIF
              ENDDO !nn1
              IF((nodeCount+specialNodeCount)>1) THEN !More than one node so it is a proper line 
                basis%NUMBER_OF_LOCAL_LINES=basis%NUMBER_OF_LOCAL_LINES+1
                basis%NUMBER_OF_NODES_IN_LOCAL_LINE(basis%NUMBER_OF_LOCAL_LINES)=nodeCount+specialNodeCount
                basis%NODE_NUMBERS_IN_LOCAL_LINE(1:basis%NUMBER_OF_NODES_IN_LOCAL_LINE(basis%NUMBER_OF_LOCAL_LINES), &
                  & basis%NUMBER_OF_LOCAL_LINES)=nodesInLine(1:basis%NUMBER_OF_NODES_IN_LOCAL_LINE(basis%NUMBER_OF_LOCAL_LINES))
                localLineParameter=0
                DO localLineNodeIdx=1,basis%NUMBER_OF_NODES_IN_LOCAL_LINE(basis%NUMBER_OF_LOCAL_LINES)
                  DO derivativeIdx=1,basis%NUMBER_OF_DERIVATIVES(basis%NODE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx, &
                    & basis%NUMBER_OF_LOCAL_LINES))
                    IF(basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,basis%NODE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx, &
                      & basis%NUMBER_OF_LOCAL_LINES),xiIdx2)==NO_PART_DERIV) THEN
                      localLineParameter=localLineParameter+1
                      basis%ELEMENT_PARAMETERS_IN_LOCAL_LINE(localLineParameter,basis%NUMBER_OF_LOCAL_LINES)= &
                        & basis%ELEMENT_PARAMETER_INDEX(derivativeIdx,basis%NODE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx, &
                        & basis%NUMBER_OF_LOCAL_LINES))
                      IF(basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,basis%NODE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx, &
                        & basis%NUMBER_OF_LOCAL_LINES),xiIdx1)==FIRST_PART_DERIV) THEN
                        basis%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx,basis%NUMBER_OF_LOCAL_LINES)=derivativeIdx
                      ENDIF
                    ENDIF
                  ENDDO !derivativeIdx
                ENDDO !localLineNodeIdx
                basis%localLineXiDirection(basis%NUMBER_OF_LOCAL_LINES)=xiIdx1
                IF(localNodeIdx2==1) THEN
                  basis%localLineXiNormals(1,basis%NUMBER_OF_LOCAL_LINES)=-xiIdx2
                  basis%xiNormalsLocalLine(-xiIdx2,1)=basis%NUMBER_OF_LOCAL_LINES
                ELSE
                  basis%locallineXiNormals(1,basis%NUMBER_OF_LOCAL_LINES)=xiIdx2
                  basis%xiNormalsLocalLine(xiIdx2,1)=basis%NUMBER_OF_LOCAL_LINES
                ENDIF
              ENDIF
            ENDDO !localNodeIdx2
          ENDDO !localNodeIdx1
        CASE(3) !3 xi directions
          !Determine the maximum node extent of the basis
          maximumNodeExtent(1)=MAXVAL(basis%NODE_POSITION_INDEX(:,1))
          maximumNodeExtent(2)=MAXVAL(basis%NODE_POSITION_INDEX(:,2))
          maximumNodeExtent(3)=MAXVAL(basis%NODE_POSITION_INDEX(:,3))
          !Allocate and calculate the lines
          IF(basis%NUMBER_OF_COLLAPSED_XI==1) THEN
            numberOfLocalLines=9
            numberOfLocalFaces=5
          ELSE IF(basis%NUMBER_OF_COLLAPSED_XI==2) THEN
            numberOfLocalLines=8
            numberOfLocalFaces=5
          ELSE
            numberOfLocalLines=12
            numberOfLocalFaces=6
          ENDIF
          basis%NUMBER_OF_LOCAL_FACES=numberOfLocalFaces

          ALLOCATE(basis%localLineXiDirection(numberOfLocalLines),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate local line xi direction.",err,error,*999)
          ALLOCATE(basis%localLineXiNormals(2,numberOfLocalLines),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate local line xi normals.",err,error,*999)
          ALLOCATE(basis%xiNormalsLocalLine(-3:3,-3:3),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate xi normals local line.",err,error,*999)
          basis%xiNormalsLocalLine=0
          
          ALLOCATE(basis%NUMBER_OF_NODES_IN_LOCAL_LINE(numberOfLocalLines),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate number of nodes in local line.",err,error,*999)
          basis%NUMBER_OF_NODES_IN_LOCAL_LINE=0

          ALLOCATE(basis%NUMBER_OF_NODES_IN_LOCAL_FACE(numberOfLocalFaces),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate number of nodes in local face.",err,error,*999)
          basis%NUMBER_OF_NODES_IN_LOCAL_FACE=0

          ALLOCATE(basis%NODE_NUMBERS_IN_LOCAL_LINE(MAXVAL(basis%NUMBER_OF_NODES_XIC),numberOfLocalLines),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate node numbers in local line.",err,error,*999)
          basis%NODE_NUMBERS_IN_LOCAL_LINE=0

          ALLOCATE(basis%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(MAXVAL(basis%NUMBER_OF_NODES_XIC),numberOfLocalLines),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate derivative numbers in local line.",err,error,*999)
          basis%DERIVATIVE_NUMBERS_IN_LOCAL_LINE=NO_PART_DERIV

          ALLOCATE(basis%ELEMENT_PARAMETERS_IN_LOCAL_LINE(MAXVAL(basis%NUMBER_OF_NODES_XIC)**2,numberOfLocalLines),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate element parameters in local line.",err,error,*999)
          basis%ELEMENT_PARAMETERS_IN_LOCAL_LINE=1

          ALLOCATE(BASIS%localFaceXiDirections(2,numberOfLocalFaces),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate local face xi directions.",err,error,*999)
          ALLOCATE(basis%localFaceXiNormal(numberOfLocalFaces),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate local face xi direction.",err,error,*999)
          ALLOCATE(basis%xiNormalLocalFace(-3:3),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate xi normal local face.",err,error,*999)
          basis%xiNormalLocalFace=0
          
          ALLOCATE(basis%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(0:basis%MAXIMUM_NUMBER_OF_DERIVATIVES, &
            & MAXVAL(basis%NUMBER_OF_NODES_XIC)**2,numberOfLocalFaces),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate derivative numbers in local face.",err,error,*999)
          basis%DERIVATIVE_NUMBERS_IN_LOCAL_FACE=NO_PART_DERIV
          basis%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(0,:,:)=1

          ALLOCATE(basis%ELEMENT_PARAMETERS_IN_LOCAL_FACE(MAXVAL(basis%NUMBER_OF_NODES_XIC)**2* &
            & basis%MAXIMUM_NUMBER_OF_DERIVATIVES,numberOfLocalFaces),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate element parameters in local face.",err,error,*999)
          basis%ELEMENT_PARAMETERS_IN_LOCAL_FACE=1

          ALLOCATE(basis%NODE_NUMBERS_IN_LOCAL_FACE(MAX(maximumNodeExtent(2)*maximumNodeExtent(3), &
            & maximumNodeExtent(3)*maximumNodeExtent(1),maximumNodeExtent(2)*maximumNodeExtent(1)), &
            & numberOfLocalFaces),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate node numbers in local face.",err,error,*999)
          basis%NODE_NUMBERS_IN_LOCAL_FACE=0
          

          !Find the lines and faces
          basis%NUMBER_OF_LOCAL_LINES=0
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
                DO localNodeIdx1=1,basis%NUMBER_OF_NODES
                  IF(basis%COLLAPSED_XI(xiIdx1)/=BASIS_NOT_COLLAPSED) THEN
                    !The current xi direction, xiIdx1, is involved in a collapsed (degenerate) plane 
                    IF(basis%COLLAPSED_XI(xiIdx2)==BASIS_XI_COLLAPSED.AND.basis%COLLAPSED_XI(xiIdx3)==BASIS_XI_COLLAPSED) THEN
                      !Both of the other two xi directions are collapsed
                      IF(basis%COLLAPSED_XI(xiIdx1)==BASIS_COLLAPSED_AT_XI0) THEN !Collapsed at the xi=0 end
                        IF((basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx2)==localNodeIdx2.OR. &
                          & basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx1)==1).AND. &
                          & (basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx3)==localNodeIdx3.OR. &
                          & basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx1)==1)) THEN
                          nodeCount=nodeCount+1
                          nodesInLine(nodeCount)=localNodeIdx1
                        ENDIF
                      ELSE !Collapsed at the xi=1 end
                        IF(basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx2)==localNodeIdx2.AND. &
                          & basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx3)==localNodeIdx3) THEN
                          nodeCount=nodeCount+1
                          nodesInLine(nodeCount)=localNodeIdx1
                        ELSE IF(basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx1)==maximumNodeExtent(xiIdx1)) THEN
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
                      IF(basis%COLLAPSED_XI(xiIdx2)==BASIS_XI_COLLAPSED) THEN
                        !The other xiIdx2 xi direction is collapsed
                        IF(basis%COLLAPSED_XI(xiIdx1)==BASIS_COLLAPSED_AT_XI0) THEN !Collapsed at the xi=0 end
                          IF((basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx2)==localNodeIdx2.OR. &
                            & basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx1)==1).AND. &
                            & basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx3)==localNodeIdx3) THEN
                            nodeCount=nodeCount+1
                            nodesInLine(nodeCount)=localNodeIdx1
                          ENDIF
                        ELSE IF(basis%COLLAPSED_XI(xiIdx1)==BASIS_COLLAPSED_AT_XI1) THEN !Collapsed at the xi=1 end
                          IF(basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx2)==localNodeIdx2.AND. &
                            & basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx3)==localNodeIdx3) THEN
                            nodeCount=nodeCount+1
                            nodesInLine(nodeCount)=localNodeIdx1
                          ELSE IF(basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx1)==maximumNodeExtent(xiIdx1).AND. &
                            & basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx3)==localNodeIdx3) THEN
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
                          IF(basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx2)==localNodeIdx2.AND. &
                            & basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx3)==localNodeIdx3) THEN
                            nodeCount=nodeCount+1
                            nodesInLine(nodeCount)=localNodeIdx1
                          ENDIF
                        ENDIF
                      ELSE IF(basis%COLLAPSED_XI(xiIdx3)==BASIS_XI_COLLAPSED) THEN
                        !The other xiIdx3 xi direction is collapsed
                        IF(basis%COLLAPSED_XI(xiIdx1)==BASIS_COLLAPSED_AT_XI0) THEN !Collapsed at the xi=0 end
                          IF(basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx2)==localNodeIdx2.AND. &
                            & (basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx3)==localNodeIdx3.OR. &
                            & basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx1)==1)) THEN
                            nodeCount=nodeCount+1
                            nodesInLine(nodeCount)=localNodeIdx1
                          ENDIF
                        ELSE IF(basis%COLLAPSED_XI(xiIdx1)==BASIS_COLLAPSED_AT_XI1) THEN !Collapsed at the xi=1 end
                          IF(basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx2)==localNodeIdx2.AND. &
                            & basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx3)==localNodeIdx3) THEN
                            nodeCount=nodeCount+1
                            nodesInLine(nodeCount)=localNodeIdx1
                          ELSE IF(basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx1)==maximumNodeExtent(xiIdx1).AND. &
                            & basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx2)==localNodeIdx2) THEN
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
                          IF(basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx2)==localNodeIdx2.AND. &
                            & basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx3)==localNodeIdx3) THEN
                            nodeCount=nodeCount+1
                            nodesInLine(nodeCount)=localNodeIdx1
                          ENDIF
                        ENDIF
                      ELSE
                        !The current xi must be collapsed
                        IF(basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx2)==localNodeIdx2.AND. &
                          & basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx3)==localNodeIdx3) THEN
                          nodeCount=nodeCount+1
                          nodesInLine(nodeCount)=localNodeIdx1
                        ENDIF
                      ENDIF
                    ENDIF
                  ELSE
                    !The current xi direction, xiIdx1, is not involved in any collapsed (degenerate) planes
                    IF(basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx2)==localNodeIdx2.AND. &
                      & basis%NODE_POSITION_INDEX(localNodeIdx1,xiIdx3)==localNodeIdx3) THEN
                      nodeCount=nodeCount+1
                      nodesInLine(nodeCount)=localNodeIdx1
                    ENDIF
                  ENDIF
                ENDDO !localNodeIdx1
                IF((nodeCount+specialNodeCount)>1) THEN !More than one node so it is a proper line 
                  basis%NUMBER_OF_LOCAL_LINES=basis%NUMBER_OF_LOCAL_LINES+1
                  basis%NUMBER_OF_NODES_IN_LOCAL_LINE(basis%NUMBER_OF_LOCAL_LINES)=nodeCount+specialNodeCount
                  basis%NODE_NUMBERS_IN_LOCAL_LINE(1:basis%NUMBER_OF_NODES_IN_LOCAL_LINE(basis%NUMBER_OF_LOCAL_LINES), &
                    & basis%NUMBER_OF_LOCAL_LINES)=nodesInLine(1:basis%NUMBER_OF_NODES_IN_LOCAL_LINE(basis%NUMBER_OF_LOCAL_LINES))
                  localLineParameter=0
                  DO localLineNodeIdx=1,basis%NUMBER_OF_NODES_IN_LOCAL_LINE(basis%NUMBER_OF_LOCAL_LINES)
                    DO derivativeIdx=1,basis%NUMBER_OF_DERIVATIVES(basis%NODE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx, &
                      & basis%NUMBER_OF_LOCAL_LINES))
                      IF(basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,basis%NODE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx, &
                        & basis%NUMBER_OF_LOCAL_LINES),xiIdx2)==NO_PART_DERIV.AND. &
                        & basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,basis%NODE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx, &
                        & basis%NUMBER_OF_LOCAL_LINES),xiIdx3)==NO_PART_DERIV) THEN
                        localLineParameter=localLineParameter+1
                        basis%ELEMENT_PARAMETERS_IN_LOCAL_LINE(localLineParameter,basis%NUMBER_OF_LOCAL_LINES)= &
                          & basis%ELEMENT_PARAMETER_INDEX(derivativeIdx,basis%NODE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx, &
                          & basis%NUMBER_OF_LOCAL_LINES))
                        IF(basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,basis%NODE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx, &
                          & basis%NUMBER_OF_LOCAL_LINES),xiIdx1)==FIRST_PART_DERIV) THEN
                          basis%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx,basis%NUMBER_OF_LOCAL_LINES)=derivativeIdx
                        ENDIF
                      ENDIF
                    ENDDO !derivativeIdx
                  ENDDO !localLineNodeIdx
                  basis%localLineXiDirection(basis%NUMBER_OF_LOCAL_LINES)=xiIdx1
                  IF(localNodeIdx2==1) THEN
                    basis%localLineXiNormals(1,basis%NUMBER_OF_LOCAL_LINES)=-xiIdx2
                    IF(localNodeIdx3==1) THEN
                      basis%localLineXiNormals(2,basis%NUMBER_OF_LOCAL_LINES)=-xiIdx3
                      basis%xiNormalsLocalLine(-xiIdx2,-xiIdx3)=basis%NUMBER_OF_LOCAL_LINES
                    ELSE
                      basis%localLineXiNormals(2,basis%NUMBER_OF_LOCAL_LINES)=xiIdx3
                      basis%xiNormalsLocalLine(-xiIdx2,xiIdx3)=basis%NUMBER_OF_LOCAL_LINES
                    ENDIF
                  ELSE
                    basis%localLineXiNormals(1,basis%NUMBER_OF_LOCAL_LINES)=xiIdx2
                    IF(localNodeIdx3==1) THEN
                      basis%localLineXiNormals(2,basis%NUMBER_OF_LOCAL_LINES)=-xiIdx3
                      basis%xiNormalsLocalLine(xiIdx2,-xiIdx3)=basis%NUMBER_OF_LOCAL_LINES
                    ELSE
                      basis%localLineXiNormals(2,basis%NUMBER_OF_LOCAL_LINES)=xiIdx3
                      basis%xiNormalsLocalLine(xiIdx2,xiIdx3)=basis%NUMBER_OF_LOCAL_LINES
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
                collapsedFace=basis%COLLAPSED_XI(xiIdx1)==BASIS_COLLAPSED_AT_XI1
              ELSE
                !The -'ve xi direction
                localNodeIdx1=1
                !Compute if the face in the -xiIdx1 direction is collapsed.
                collapsedFace=basis%COLLAPSED_XI(xiIdx1)==BASIS_COLLAPSED_AT_XI0
              ENDIF
              localNodeCount=0
              IF(.NOT.collapsedFace) THEN
                !If the face has not been collapsed
                localFaceIdx=localFaceIdx+1
                !Loop over the local nodes in the face
                DO localNodeIdx3=1,maximumNodeExtent(xiIdx2)
                  DO localNodeIdx2=1,maximumNodeExtent(xiIdx3)
                    IF(xiIdx1==1) THEN
                      localNodeIdx=basis%NODE_POSITION_INDEX_INV(localNodeIdx1,localNodeIdx2,localNodeIdx3,1)
                    ELSE IF(xiIdx1==2) THEN
                      localNodeIdx=basis%NODE_POSITION_INDEX_INV(localNodeIdx2,localNodeIdx1,localNodeIdx3,1)
                    ELSE
                      localNodeIdx=basis%NODE_POSITION_INDEX_INV(localNodeIdx2,localNodeIdx3,localNodeIdx1,1)
                    ENDIF
                    IF(ALL(basis%NODE_NUMBERS_IN_LOCAL_FACE(1:localNodeCount,localFaceIdx)/=localNodeIdx)) THEN
                      !The node hasn't been collapsed
                      localNodeCount=localNodeCount+1
                      basis%NODE_NUMBERS_IN_LOCAL_FACE(localNodeCount,localFaceIdx)=localNodeIdx
                    ENDIF
                  ENDDO !localNodeIdx3
                ENDDO !localNodexIdx2
                basis%NUMBER_OF_NODES_IN_LOCAL_FACE(localFaceIdx)=localNodeCount
                basis%localFaceXiDirections(1,localFaceIdx)=xiIdx2
                basis%localFaceXiDirections(2,localFaceIdx)=xiIdx3
                basis%localFaceXiNormal(localFaceIdx)=directionIdx*xiIdx1
                basis%xiNormalLocalFace(directionIdx*xiIdx1)=localFaceIdx
                !Compute derivatives and element parameters in the face
                localFaceParameter=0
                DO localNodeIdx=1,basis%NUMBER_OF_NODES_IN_LOCAL_FACE(localFaceIdx)
                  localNode=basis%NODE_NUMBERS_IN_LOCAL_FACE(localNodeIdx,localFaceIdx)
                  localFaceDerivative=0
                  DO derivativeIdx=1,basis%NUMBER_OF_DERIVATIVES(localNode)
                    IF(basis%DERIVATIVE_ORDER_INDEX(derivativeIdx,localNode,xiIdx1)==NO_PART_DERIV) THEN
                      localFaceParameter=localFaceParameter+1
                      localFaceDerivative=localFaceDerivative+1
                      basis%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(localFaceDerivative,localNodeIdx,localFaceIdx)=derivativeIdx
                      basis%ELEMENT_PARAMETERS_IN_LOCAL_FACE(localFaceParameter,localFaceIdx)= &
                        & basis%ELEMENT_PARAMETER_INDEX(derivativeIdx,localNode)
                    ENDIF
                  ENDDO !derivativeIdx
                  basis%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(0,localNodeIdx,localFaceIdx)=localFaceDerivative
                ENDDO !localNodeIdx
              ENDIF
            ENDDO !xiIdx1
          ENDDO !directionIdx
        CASE DEFAULT
          CALL FlagError("Invalid number of xi directions.",err,error,*999)
        END SELECT
        
        CALL BASIS_QUADRATURE_CREATE(basis,err,error,*999)

      ELSE
        CALL FlagError("Basis is not a Lagrange Hermite tensor product basis.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated.",err,error,*999)
    ENDIF
    
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
  
  !>Evaluates the double precision Lagrange/Hermite/Fourier tensor product basis function for the given BASIS.
  FUNCTION BASIS_LHTP_BASIS_EVALUATE_DP(BASIS,NODE_NUMBER,DERIVATIVE_NUMBER,PARTIAL_DERIV_INDEX,XI,err,error)
      
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to evaluate 
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBER !<The local node number of the tensor product basis to evaluate
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The local derivative number of the tensor product basis to evaluate
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIV_INDEX  !<The partial derivative index to interpolate \see CONSTANTS_PartialDerivativeConstants
    REAL(DP), INTENT(IN) :: XI(:) !<The Xi position to evaluate the basis function at
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: BASIS_LHTP_BASIS_EVALUATE_DP !<On return the evaluated basis funtion.
    !Local variables
    INTEGER(INTG) :: ni,nn
    REAL(DP) :: SUM
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("BASIS_LHTP_BASIS_EVALUATE_DP",err,error,*999)
    
    BASIS_LHTP_BASIS_EVALUATE_DP=1.0_DP
    IF(ASSOCIATED(BASIS)) THEN
      DO ni=1,BASIS%NUMBER_OF_XI
        IF(BASIS%NODE_AT_COLLAPSE(NODE_NUMBER).AND.BASIS%COLLAPSED_XI(ni)==BASIS_XI_COLLAPSED) THEN
          !We are at a collapsed node in the collapsed xi direction. Sum the basis functions in the collapsed xi direction.
          SUM=0.0_DP
          SELECT CASE(BASIS%INTERPOLATION_TYPE(ni))
          CASE(BASIS_LAGRANGE_INTERPOLATION)
            SELECT CASE(BASIS%INTERPOLATION_ORDER(ni))
            CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
              DO nn=1,2
                SUM=SUM+LAGRANGE_LINEAR_EVALUATE(nn,PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,ni),XI(ni),err,error)
              ENDDO !nn
            CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
              DO nn=1,3
                SUM=SUM+LAGRANGE_QUADRATIC_EVALUATE(nn,PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,ni),XI(ni),err,error)
              ENDDO !nn
            CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
              DO nn=1,4
                SUM=SUM+LAGRANGE_CUBIC_EVALUATE(nn,PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,ni),XI(ni),err,error)
              ENDDO !nn
            CASE DEFAULT
              LOCAL_ERROR="Interpolation order value "//TRIM(NumberToVString(BASIS%INTERPOLATION_ORDER(ni),"*",err,error))// &
                & " for xi direction "//TRIM(NumberToVString(ni,"*",err,error))//" is invalid"
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          CASE(BASIS_HERMITE_INTERPOLATION)
            SELECT CASE(BASIS%INTERPOLATION_ORDER(ni))
            CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
              DO nn=1,2
                SUM=SUM+HERMITE_CUBIC_EVALUATE(nn,BASIS%DERIVATIVE_ORDER_INDEX(DERIVATIVE_NUMBER,NODE_NUMBER,ni), &
                  & PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,ni),XI(ni),err,error)
              ENDDO !nn
            CASE DEFAULT
              LOCAL_ERROR="Interpolation order value "//TRIM(NumberToVString(BASIS%INTERPOLATION_ORDER(ni),"*",err,error))// &
                & " for xi direction "//TRIM(NumberToVString(ni,"*",err,error))//" is invalid"
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
            IF(ERR/=0) GOTO 999
          CASE DEFAULT
            LOCAL_ERROR="Interpolation type value "//TRIM(NumberToVString(BASIS%INTERPOLATION_TYPE(ni),"*",err,error))// &
              & " for xi direction "//TRIM(NumberToVString(ni,"*",err,error))//" is invalid"
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
          BASIS_LHTP_BASIS_EVALUATE_DP=BASIS_LHTP_BASIS_EVALUATE_DP*SUM
        ELSE
          SELECT CASE(BASIS%INTERPOLATION_TYPE(ni))
          CASE(BASIS_LAGRANGE_INTERPOLATION)
            SELECT CASE(BASIS%INTERPOLATION_ORDER(ni))
            CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
              BASIS_LHTP_BASIS_EVALUATE_DP=BASIS_LHTP_BASIS_EVALUATE_DP* &
                & LAGRANGE_LINEAR_EVALUATE(BASIS%NODE_POSITION_INDEX(NODE_NUMBER,ni), &
                & PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,ni),XI(ni),err,error)
            CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
              BASIS_LHTP_BASIS_EVALUATE_DP=BASIS_LHTP_BASIS_EVALUATE_DP* &
                & LAGRANGE_QUADRATIC_EVALUATE(BASIS%NODE_POSITION_INDEX(NODE_NUMBER,ni), &
                & PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,ni),XI(ni),err,error)
            CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
              BASIS_LHTP_BASIS_EVALUATE_DP=BASIS_LHTP_BASIS_EVALUATE_DP* &
                & LAGRANGE_CUBIC_EVALUATE(BASIS%NODE_POSITION_INDEX(NODE_NUMBER,ni), &
                & PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,ni),XI(ni),err,error)
            CASE DEFAULT
              LOCAL_ERROR="Interpolation order value "//TRIM(NumberToVString(BASIS%INTERPOLATION_ORDER(ni),"*",err,error))// &
                & " for xi direction "//TRIM(NumberToVString(ni,"*",err,error))//" is invalid"
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          CASE(BASIS_HERMITE_INTERPOLATION)
            SELECT CASE(BASIS%INTERPOLATION_ORDER(ni))
            CASE(BASIS_QUADRATIC1_INTERPOLATION_ORDER)
              BASIS_LHTP_BASIS_EVALUATE_DP=BASIS_LHTP_BASIS_EVALUATE_DP* &
                & HERMITE_QUADRATIC_EVALUATE(BASIS%NODE_POSITION_INDEX(NODE_NUMBER,ni), &
                & BASIS%DERIVATIVE_ORDER_INDEX(DERIVATIVE_NUMBER,NODE_NUMBER,ni), &
                & PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,ni),1,XI(ni),err,error)
            CASE(BASIS_QUADRATIC2_INTERPOLATION_ORDER)
              BASIS_LHTP_BASIS_EVALUATE_DP=BASIS_LHTP_BASIS_EVALUATE_DP* &
                & HERMITE_QUADRATIC_EVALUATE(BASIS%NODE_POSITION_INDEX(NODE_NUMBER,ni), &
                & BASIS%DERIVATIVE_ORDER_INDEX(DERIVATIVE_NUMBER,NODE_NUMBER,ni), &
                & PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,ni),2,XI(ni),err,error)
            CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
              BASIS_LHTP_BASIS_EVALUATE_DP=BASIS_LHTP_BASIS_EVALUATE_DP* &
                & HERMITE_CUBIC_EVALUATE(BASIS%NODE_POSITION_INDEX(NODE_NUMBER,ni), &
                & BASIS%DERIVATIVE_ORDER_INDEX(DERIVATIVE_NUMBER,NODE_NUMBER,ni), &
                & PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,ni),XI(ni),err,error)
            CASE DEFAULT
              LOCAL_ERROR="Interpolation order value "//TRIM(NumberToVString(BASIS%INTERPOLATION_ORDER(ni),"*",err,error))// &
                & " for xi direction "//TRIM(NumberToVString(ni,"*",err,error))//" is invalid"
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
            IF(ERR/=0) GOTO 999
          CASE DEFAULT
            LOCAL_ERROR="Interpolation type value "//TRIM(NumberToVString(BASIS%INTERPOLATION_TYPE(ni),"*",err,error))// &
              & " for xi direction "//TRIM(NumberToVString(ni,"*",err,error))//" is invalid"
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        ENDIF
      ENDDO !ni
    ELSE
      CALL FlagError("Basis is not associated",err,error,*999)
    ENDIF

    EXITS("BASIS_LHTP_BASIS_EVALUATE_DP")
    RETURN
999 ERRORSEXITS("BASIS_LHTP_BASIS_EVALUATE_DP",err,error)
    RETURN 
  END FUNCTION BASIS_LHTP_BASIS_EVALUATE_DP

  !
  !================================================================================================================================
  !

  !>Creates and initialises a Lagrange-Hermite tensor product basis family that has already been allocated by BASIS_CREATE_START.
  !>\see BASIS_ROUTINES::BASIS_LHTP_CREATE
  SUBROUTINE Basis_LHTPFamilyCreate(basis,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: basis !<A pointer to the basis to create
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,faceXi(2),faceXi2(2),localFaceIdx,localLineIdx,xiIdx,xiIdx2
    LOGICAL :: lineBasisDone,faceBasisDone
    TYPE(BASIS_TYPE), POINTER :: newSubBasis
    TYPE(VARYING_STRING) :: dummyError

    NULLIFY(newSubBasis)

    ENTERS("Basis_LHTPFamilyCreate",err,error,*999)

    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    
    !Create the main (parent) basis
    CALL Basis_LHTPBasisCreate(basis,err,error,*999)
    
    IF(basis%NUMBER_OF_XI>1) THEN
      !Create the line bases as sub-basis types
      ALLOCATE(basis%lineBases(basis%NUMBER_OF_XI),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate basis line bases.",err,error,*999)
      ALLOCATE(basis%localLineBasis(basis%NUMBER_OF_LOCAL_LINES),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate local line basis.",err,error,*999)
      DO xiIdx=1,basis%NUMBER_OF_XI
        lineBasisDone=.FALSE.
        NULLIFY(newSubBasis)
        DO xiIdx2=1,xiIdx-1
          IF(basis%INTERPOLATION_XI(xiIdx2)==basis%INTERPOLATION_XI(xiIdx).AND. &
            basis%QUADRATURE%NUMBER_OF_GAUSS_XI(xiIdx2)==basis%QUADRATURE%NUMBER_OF_GAUSS_XI(xiIdx)) THEN
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
      DO localLineIdx=1,basis%NUMBER_OF_LOCAL_LINES
        xiIdx=basis%localLineXiDirection(localLineIdx)
        basis%localLineBasis(localLineIdx)%ptr=>basis%lineBases(xiIdx)%ptr
      ENDDO !localLineIdx
      IF(basis%NUMBER_OF_XI>2) THEN
        !Create the face bases as sub-basis types
        ALLOCATE(basis%faceBases(basis%NUMBER_OF_XI),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate basis face bases.",err,error,*999)
        ALLOCATE(basis%localFaceBasis(basis%NUMBER_OF_LOCAL_FACES),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate local face basis.",err,error,*999)
        DO xiIdx=1,basis%NUMBER_OF_XI
          !Determine the face xi directions that lie in this xi direction
          faceXi(1)=OTHER_XI_DIRECTIONS3(xiIdx,2,1)
          faceXi(2)=OTHER_XI_DIRECTIONS3(xiIdx,3,1)
          faceBasisDone=.FALSE.
          NULLIFY(newSubBasis)
          DO xiIdx2=1,xiIdx-1
            !Determine the face xi directions that lie in this xi direction
            faceXi2(1)=OTHER_XI_DIRECTIONS3(xiIdx2,2,1)
            faceXi2(2)=OTHER_XI_DIRECTIONS3(xiIdx2,3,1)
            IF(basis%INTERPOLATION_XI(faceXi2(1))==basis%INTERPOLATION_XI(faceXi(1)).AND. &
              & basis%INTERPOLATION_XI(faceXi2(2))==basis%INTERPOLATION_XI(faceXi(2)).AND. &
              & basis%QUADRATURE%NUMBER_OF_GAUSS_XI(faceXi2(1))==basis%QUADRATURE%NUMBER_OF_GAUSS_XI(faceXi(1)).AND. &
              & basis%QUADRATURE%NUMBER_OF_GAUSS_XI(faceXi2(2))==basis%QUADRATURE%NUMBER_OF_GAUSS_XI(faceXi(2)).AND. &
              & basis%COLLAPSED_XI(faceXi2(1))==basis%COLLAPSED_XI(faceXi(1)).AND. &
              & basis%COLLAPSED_XI(faceXi2(2))==basis%COLLAPSED_XI(faceXi(2))) THEN              
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
            ALLOCATE(newSubBasis%localLineBasis(newSubBasis%NUMBER_OF_LOCAL_LINES),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate new sub basis local line basis.",err,error,*999)
            DO localLineIdx=1,newSubBasis%NUMBER_OF_LOCAL_LINES
              IF(newSubBasis%localLineXiDirection(localLineIdx)==1) THEN
                newSubBasis%localLineBasis(localLineIdx)%ptr=>newSubBasis%lineBases(1)%ptr
              ELSE IF(newSubBasis%localLineXiDirection(localLineIdx)==2) THEN
                newSubBasis%localLineBasis(localLineIdx)%ptr=>newSubBasis%lineBases(2)%ptr
              ENDIF
            ENDDO !localFaceIdx
          ENDIF
        ENDDO !xiIdx
        DO localFaceIdx=1,basis%NUMBER_OF_LOCAL_FACES
          xiIdx=basis%localFaceXiNormal(localFaceIdx)
          basis%localFaceBasis(localFaceIdx)%ptr=>basis%faceBases(ABS(xiIdx))%ptr
        ENDDO !localFaceIdx
      ENDIF
    ENDIF
    
    EXITS("Basis_LHTPFamilyCreate")
    RETURN
999 IF(ASSOCIATED(newSubBasis)) CALL Basis_FamilyDestroy(newSubBasis%USER_NUMBER,newSubBasis%FAMILY_NUMBER, &
      & dummyErr,dummyError,*998)
998 ERRORSEXITS("Basis_LHTPFamilyCreate",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_LHTPFamilyCreate

  !
  !================================================================================================================================
  !

  !>Creates and initialises a Radial basis family that has already been allocated by BASIS_CREATE_START.
  !>\see BASIS_ROUTINES::BASIS_RADIAL_CREATE
  SUBROUTINE BASIS_RADIAL_FAMILY_CREATE(BASIS,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to create
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("BASIS_RADIAL_FAMILY_CREATE",err,error,*999)

    IF(ASSOCIATED(BASIS)) THEN
      !Create the main (parent) basis
      CALL FlagError("Not implemented.",err,error,*999)
    ELSE
      CALL FlagError("Basis is not associated",err,error,*999)
    ENDIF

    EXITS("BASIS_RADIAL_FAMILY_CREATE")
    RETURN
999 ERRORSEXITS("BASIS_RADIAL_FAMILY_CREATE",err,error)
    RETURN 1
    
  END SUBROUTINE BASIS_RADIAL_FAMILY_CREATE

  !
  !================================================================================================================================
  !

  !>Calculates the xi location of a local node in a basis.
  SUBROUTINE BASIS_LOCAL_NODE_XI_CALCULATE(BASIS,LOCAL_NODE_NUMBER,XI,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to calculate the xi for
    INTEGER(INTG), INTENT(IN) :: LOCAL_NODE_NUMBER !<The local node number to calculate the xi for
    REAL(DP), INTENT(OUT) :: XI(:) !<XI(ni). On return, the xi position of the local node in the basis
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: xi_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("BASIS_LOCAL_NODE_XI_CALCULATE",err,error,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        IF(LOCAL_NODE_NUMBER>0.AND.LOCAL_NODE_NUMBER<=BASIS%NUMBER_OF_NODES) THEN
          IF(SIZE(XI,1)>=BASIS%NUMBER_OF_XI) THEN
            SELECT CASE(BASIS%TYPE)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
              DO xi_idx=1,BASIS%NUMBER_OF_XI
                XI(xi_idx)=REAL(BASIS%NODE_POSITION_INDEX(LOCAL_NODE_NUMBER,xi_idx)-1,DP)/ &
                  & REAL(BASIS%NUMBER_OF_NODES_XIC(xi_idx)-1,DP)
              ENDDO !xi_idx
            CASE(BASIS_SIMPLEX_TYPE)
              DO xi_idx=1,BASIS%NUMBER_OF_XI
                XI(xi_idx)=REAL(BASIS%NUMBER_OF_NODES_XIC(xi_idx)-BASIS%NODE_POSITION_INDEX(LOCAL_NODE_NUMBER,xi_idx),DP)/ &
                  & REAL(BASIS%NUMBER_OF_NODES_XIC(xi_idx)-1,DP)
              ENDDO !xi_idx
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
              LOCAL_ERROR="The basis type of "//TRIM(NumberToVString(BASIS%TYPE,"*",err,error))// &
                & " is invalid."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          ELSE
            LOCAL_ERROR="The size of the specified xic array of "//TRIM(NumberToVString(SIZE(XI,1),"*",err,error))// &
              & " is invalid. The size of the xi array must be >= "// &
              & TRIM(NumberToVString(BASIS%NUMBER_OF_XI,"*",err,error))//"."            
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified local node number of "//TRIM(NumberToVString(LOCAL_NODE_NUMBER,"*",err,error))// &
            & " is invalid. The local node number must be > 0 and <= "// &
            & TRIM(NumberToVString(BASIS%NUMBER_OF_NODES,"*",err,error))//"."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Basis has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated",err,error,*999)
    ENDIF
    
    EXITS("BASIS_LOCAL_NODE_XI_CALCULATE")
    RETURN
999 ERRORSEXITS("BASIS_LOCAL_NODE_XI_CALCULATE",err,error)
    RETURN 1
  END SUBROUTINE BASIS_LOCAL_NODE_XI_CALCULATE

  !
  !================================================================================================================================
  !

  !>Returns the number of local nodes in the specified basis \see OpenCMISS::Iron::cmfe_BasisNumberOfLocalNodesGet
  SUBROUTINE BASIS_NUMBER_OF_LOCAL_NODES_GET(BASIS,NUMBER_OF_LOCAL_NODES,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to get the number of nodes
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_LOCAL_NODES !<On return, the number of local nodes in the basis
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("BASIS_NUMBER_OF_LOCAL_NODES_GET",err,error,*999)

    IF(ASSOCIATED(BASIS)) THEN
      NUMBER_OF_LOCAL_NODES=BASIS%NUMBER_OF_NODES
    ELSE
      CALL FlagError("Basis is not associated",err,error,*999)
    ENDIF
    
    EXITS("BASIS_NUMBER_OF_LOCAL_NODES_GET")
    RETURN
999 ERRORSEXITS("BASIS_NUMBER_OF_LOCAL_NODES_GET",err,error)
    RETURN 1
  END SUBROUTINE BASIS_NUMBER_OF_LOCAL_NODES_GET

  !
  !================================================================================================================================
  !
  
  !>Gets the number of xi directions for a basis. \see OpenCMISS::Iron::cmfe_Basis_NumberOfXiGet
  SUBROUTINE Basis_NumberOfXiGet(basis,numberOfXi,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: basis !<A pointer to the basis function to get the number of xi for
    INTEGER(INTG), INTENT(OUT) :: numberOfXi !<On return, the number of Xi directions for the specified basis.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Basis_NumberOfXiGet",err,error,*999)

    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(.NOT.basis%BASIS_FINISHED) CALL FlagError("Basis has not been finished.",err,error,*999)
    
    numberOfXi=basis%NUMBER_OF_XI
   
    EXITS("Basis_NumberOfXiGet")
    RETURN
999 ERRORSEXITS("Basis_NumberOfXiGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_NumberOfXiGet

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of xi directions for a basis. \see OpenCMISS::Iron::cmfe_Basis_NumberOfXiSet
  SUBROUTINE Basis_NumberOfXiSet(basis,numberOfXi,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: basis !<A pointer to the basis function to set the number of xi for
    INTEGER(INTG), INTENT(IN) :: numberOfXi !<The number of Xi directions to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG), ALLOCATABLE :: newCollapsedXi(:),newInterpolationXi(:),newNumberOfGaussXi(:)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Basis_NumberOfXiSet",err,error,*999)

    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(basis%BASIS_FINISHED) CALL FlagError("Basis has already been finished.",err,error,*999)
    IF(numberOfXi<1.OR.numberOfXi>3) THEN
      localError="The specified number of xi directions of "//TRIM(NumberToVString(numberOfXi,"*",err,error))// &
        & " is invalid. The number of xi directions must be >= 1 and <=3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    IF(basis%NUMBER_OF_XI/=numberOfXi) THEN     
      !Reallocate the basis information arrays that depend on the number of xi directions
      ALLOCATE(newInterpolationXi(numberOfXi),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new interpolation xi.",err,error,*999)
      ALLOCATE(newCollapsedXi(numberOfXi),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new COLLAPSED xi.",err,error,*999)
      IF(numberOfXi>basis%NUMBER_OF_XI) THEN
        newInterpolationXi(1:basis%NUMBER_OF_XI)=basis%INTERPOLATION_XI(1:basis%NUMBER_OF_XI)
        newInterpolationXi(basis%NUMBER_OF_XI+1:numberOfXi)=basis%INTERPOLATION_XI(1)
        newCollapsedXi(1:basis%NUMBER_OF_XI)=basis%COLLAPSED_XI(1:basis%NUMBER_OF_XI)
        newCollapsedXi(basis%NUMBER_OF_XI+1:numberOfXi)=basis%COLLAPSED_XI(1)
      ELSE
        newInterpolationXi(1:numberOfXi)=basis%INTERPOLATION_XI(1:numberOfXi)
        newCollapsedXi(1:numberOfXi)=basis%COLLAPSED_Xi(1:numberOfXi)
      ENDIF
      CALL MOVE_ALLOC(newInterpolationXi,basis%INTERPOLATION_XI)
      CALL MOVE_ALLOC(newCollapsedXi,basis%COLLAPSED_XI)
      IF(ASSOCIATED(basis%quadrature%basis)) THEN
        ALLOCATE(newNumberOfGaussXi(numberOfXi),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new number of Gauss xi.",err,error,*999)
        IF(numberOfXi>basis%NUMBER_OF_XI) THEN
          newNumberOfGaussXi(1:basis%NUMBER_OF_XI)=basis%quadrature%NUMBER_OF_GAUSS_XI(1:basis%NUMBER_OF_XI)
          newNumberOfGaussXi(basis%NUMBER_OF_XI+1:numberOfXi)=basis%quadrature%NUMBER_OF_GAUSS_XI(1)
        ELSE
          newNumberOfGaussXi(1:numberOfXi)=basis%quadrature%NUMBER_OF_GAUSS_XI(1:numberOfXi)
        ENDIF
        CALL MOVE_ALLOC(newNumberOfGaussXi,basis%quadrature%NUMBER_OF_GAUSS_XI)        
      ENDIF
      basis%NUMBER_OF_XI=numberOfXi
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
  SUBROUTINE BASIS_QUADRATURE_CREATE(BASIS,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: scheme_idx,i,j,k,MAX_NUM_GAUSS,ng,ni,nk,nn,ns,nu,NUM_GAUSS_1,NUM_GAUSS_2,NUM_GAUSS_3
    REAL(DP) :: XI(3),GSX(4,20),GSW(20)
    REAL(DP), ALLOCATABLE :: POSITIONS(:,:),POSITIONS_MATRIX(:,:,:,:),WEIGHTS(:,:)
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: NEW_SCHEME,SCHEME
    TYPE(QUADRATURE_SCHEME_PTR_TYPE), POINTER :: NEW_SCHEMES(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: MAX_NUM_FACE_GAUSS,face_idx,NORMAL,FACE_XI(2),numberOfFaceXiCoordinates

    NULLIFY(NEW_SCHEME)
    NULLIFY(NEW_SCHEMES)

    ENTERS("BASIS_QUADRATURE_CREATE",err,error,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(ASSOCIATED(BASIS%QUADRATURE%SCHEMES)) THEN
        LOCAL_ERROR="The quadrature schemes on basis number "//TRIM(NumberToVString(BASIS%USER_NUMBER,"*",err,error))// &
          & " are already associated"
        CALL FlagError(LOCAL_ERROR,err,error,*998)
      ELSE
!!TODO: \todo Sort this properly by having a create values cache.
        !Reset the basis quadrature - 
        !CALL BASIS_QUADRATURE_FINALISE(BASIS,err,error,*999)      ! Kumar - I think this is not correct as it 
        !Initialise the basis quadrature                           !         resets the quadrature scheme already set.
        !CALL BASIS_QUADRATURE_INITIALISE(BASIS,err,error,*999)    ! 
        SELECT CASE(BASIS%QUADRATURE%TYPE)
        CASE(BASIS_GAUSS_LEGENDRE_QUADRATURE)            
          !Allocate one scheme and add it to the list of schemes
          ALLOCATE(NEW_SCHEME,STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate new quadrature scheme",err,error,*999)
          NEW_SCHEME%QUADRATURE=>BASIS%QUADRATURE
          BASIS%QUADRATURE%NUMBER_OF_SCHEMES=1
          ALLOCATE(NEW_SCHEMES(BASIS%QUADRATURE%NUMBER_OF_SCHEMES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate new quadratures scheme",err,error,*999)
          NEW_SCHEMES(1)%PTR=>NEW_SCHEME
          BASIS%QUADRATURE%SCHEMES=>NEW_SCHEMES
          !Set up the quadrature scheme map
          ALLOCATE(BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate quadrature scheme map",err,error,*999)
          DO scheme_idx=1,BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES
            NULLIFY(BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(scheme_idx)%PTR)
          ENDDO !scheme_idx
          BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR=>NEW_SCHEME
          !Set up the gauss point arrays
          NEW_SCHEME%NUMBER_OF_GAUSS=1
          MAX_NUM_GAUSS=-1
          DO ni=1,BASIS%NUMBER_OF_XI
            NEW_SCHEME%NUMBER_OF_GAUSS=NEW_SCHEME%NUMBER_OF_GAUSS*BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)
            IF(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)>MAX_NUM_GAUSS) MAX_NUM_GAUSS=BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)
          ENDDO !ni
          ALLOCATE(NEW_SCHEME%GAUSS_POSITIONS(BASIS%NUMBER_OF_XI_COORDINATES,NEW_SCHEME%NUMBER_OF_GAUSS),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate Gauss positions",err,error,*999)
          ALLOCATE(NEW_SCHEME%GAUSS_WEIGHTS(NEW_SCHEME%NUMBER_OF_GAUSS),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate Gauss weights",err,error,*999)
          ALLOCATE(NEW_SCHEME%GAUSS_BASIS_FNS(BASIS%NUMBER_OF_ELEMENT_PARAMETERS,BASIS%NUMBER_OF_PARTIAL_DERIVATIVES, &
            & NEW_SCHEME%NUMBER_OF_GAUSS),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate Gauss basis functions",err,error,*999)
          ALLOCATE(WEIGHTS(MAX_NUM_GAUSS,3),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate weights",err,error,*999)
          ALLOCATE(POSITIONS(MAX_NUM_GAUSS,3),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate positions",err,error,*999)
          ALLOCATE(POSITIONS_MATRIX(MAX_NUM_GAUSS,MAX_NUM_GAUSS,MAX_NUM_GAUSS,3),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate positions matrix",err,error,*999)
          WEIGHTS=1.0_DP
          POSITIONS=0.0_DP
          POSITIONS_MATRIX=0.0_DP
          DO ni=1,BASIS%NUMBER_OF_XI
            CALL Gauss_Legendre(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni),0.0_DP,1.0_DP, &
              & POSITIONS(1:BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni),ni), &
              & WEIGHTS(1:BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni),ni),err,error,*999)
          ENDDO !ni
          SELECT CASE(BASIS%NUMBER_OF_XI)
          CASE(1)
            NUM_GAUSS_1=BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(1)
            NUM_GAUSS_2=1
            NUM_GAUSS_3=1
          CASE(2)
            NUM_GAUSS_1=BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(1)
            NUM_GAUSS_2=BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(2)
            NUM_GAUSS_3=1
          CASE(3)
            NUM_GAUSS_1=BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(1)
            NUM_GAUSS_2=BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(2)
            NUM_GAUSS_3=BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(3)
          CASE DEFAULT
            CALL FlagError("Invalid number of xi directions",err,error,*999)
          END SELECT
          DO k=1,NUM_GAUSS_3
            DO j=1,NUM_GAUSS_2
              DO i=1,NUM_GAUSS_1
                POSITIONS_MATRIX(i,j,k,1)=POSITIONS(i,1)
                POSITIONS_MATRIX(i,j,k,2)=POSITIONS(j,2)
                POSITIONS_MATRIX(i,j,k,3)=POSITIONS(k,3)
                XI(1:BASIS%NUMBER_OF_XI)=POSITIONS_MATRIX(i,j,k,1:BASIS%NUMBER_OF_XI)
                ng=i+(j-1+(k-1)*NUM_GAUSS_2)*NUM_GAUSS_1
                NEW_SCHEME%GAUSS_WEIGHTS(ng)=WEIGHTS(i,1)*WEIGHTS(j,2)*WEIGHTS(k,3)
                NEW_SCHEME%GAUSS_POSITIONS(1:BASIS%NUMBER_OF_XI_COORDINATES,ng)=XI(1:BASIS%NUMBER_OF_XI_COORDINATES)
                ns=0
                DO nn=1,BASIS%NUMBER_OF_NODES
                  DO nk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                    ns=ns+1
                    DO nu=1,BASIS%NUMBER_OF_PARTIAL_DERIVATIVES
                      SELECT CASE(BASIS%TYPE)
                      CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
                        NEW_SCHEME%GAUSS_BASIS_FNS(ns,nu,ng)=BASIS_LHTP_BASIS_EVALUATE(BASIS,nn,nk,nu,XI,err,error)
                        IF(ERR/=0) GOTO 999                        
                      CASE DEFAULT
                        CALL FlagError("Not implemented",err,error,*999)
                      END SELECT
                    ENDDO !nu
                  ENDDO !nk
                ENDDO !nn
              ENDDO !i
            ENDDO !j
          ENDDO !k
          !Create face quadrature scheme, if requested
          IF(BASIS%QUADRATURE%EVALUATE_FACE_GAUSS) THEN
            IF(BASIS%NUMBER_OF_XI==3) THEN
              !Find maximum number of face gauss points and allocate the arrays
              MAX_NUM_FACE_GAUSS=PRODUCT(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(1:BASIS%NUMBER_OF_XI))
              MAX_NUM_FACE_GAUSS=MAX_NUM_FACE_GAUSS/MINVAL(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(1:BASIS%NUMBER_OF_XI))
              ALLOCATE(NEW_SCHEME%NUMBER_OF_FACE_GAUSS(BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate number of face gauss",err,error,*999)
              ALLOCATE(NEW_SCHEME%FACE_GAUSS_POSITIONS(BASIS%NUMBER_OF_XI_COORDINATES,MAX_NUM_FACE_GAUSS, &
                & BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate face Gauss positions",err,error,*999)
              ALLOCATE(NEW_SCHEME%FACE_GAUSS_WEIGHTS(MAX_NUM_FACE_GAUSS,BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate face Gauss weights",err,error,*999)
              ALLOCATE(NEW_SCHEME%FACE_GAUSS_BASIS_FNS(BASIS%NUMBER_OF_ELEMENT_PARAMETERS,BASIS%NUMBER_OF_PARTIAL_DERIVATIVES, &
                & MAX_NUM_FACE_GAUSS,BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate face Gauss basis function values array",err,error,*999)
              !Zero them out just to be safe
              NEW_SCHEME%FACE_GAUSS_POSITIONS=0.0_DP
              NEW_SCHEME%FACE_GAUSS_WEIGHTS=0.0_DP
              NEW_SCHEME%FACE_GAUSS_BASIS_FNS=0.0_DP
              !Populate face_gauss_positions, weights, basis_fn
              DO face_idx=1,BASIS%NUMBER_OF_LOCAL_FACES
                !What's the normal?
                NORMAL=BASIS%localFaceXiNormal(face_idx)
                IF(NORMAL<0_INTG) THEN
                  XI(ABS(NORMAL))=0.0_DP
                ELSE
                  XI(ABS(NORMAL))=1.0_DP
                ENDIF
                NORMAL=ABS(NORMAL)
                FACE_XI=[OTHER_XI_DIRECTIONS3(NORMAL,2,1), OTHER_XI_DIRECTIONS3(NORMAL,3,1)]
                !How many gauss points are in this face?
                NEW_SCHEME%NUMBER_OF_FACE_GAUSS(face_idx)=PRODUCT(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(FACE_XI))
                ng=0_INTG
                DO j=1,BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(FACE_XI(2))
                  XI(FACE_XI(2))=POSITIONS(j,FACE_XI(2))
                  DO i=1,BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(FACE_XI(1))
                    XI(FACE_XI(1))=POSITIONS(i,FACE_XI(1))
                    ng=ng+1_INTG
                    !Gauss point xi and weights first
                    NEW_SCHEME%FACE_GAUSS_WEIGHTS(ng,face_idx)=WEIGHTS(i,FACE_XI(1))*WEIGHTS(j,FACE_XI(2))
                    NEW_SCHEME%FACE_GAUSS_POSITIONS(1:3,ng,face_idx)=XI(1:3)
                    !Evaluate basis fn values at the Gauss points now
                    ns=0
                    DO nn=1,BASIS%NUMBER_OF_NODES
                      DO nk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                        ns=ns+1
                        DO nu=1,BASIS%NUMBER_OF_PARTIAL_DERIVATIVES
                          SELECT CASE(BASIS%TYPE)
                          CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
                            NEW_SCHEME%FACE_GAUSS_BASIS_FNS(ns,nu,ng,face_idx)= &
                              & BASIS_LHTP_BASIS_EVALUATE(BASIS,nn,nk,nu,XI,err,error)
                            IF(ERR/=0) GOTO 999                        
                          CASE DEFAULT
                            CALL FlagError("Not implemented",err,error,*999)
                          END SELECT
                        ENDDO !nu
                      ENDDO !nk
                    ENDDO !nn

                  ENDDO !i
                ENDDO !j
              ENDDO !face_idx
            ELSE
              CALL FlagError("Cannot evaluate face quadrature schemes for a non three dimensional element.",err,error,*999)
            ENDIF
          ENDIF
          !Clean up
          DEALLOCATE(WEIGHTS)
          DEALLOCATE(POSITIONS)
          DEALLOCATE(POSITIONS_MATRIX)
        CASE(BASIS_GAUSS_LAGUERRE_QUADRATURE)
          CALL FlagError("Gauss Laguerre quadrature type not implemented.",err,error,*999)
        CASE(BASIS_GUASS_HERMITE_QUADRATURE)
          CALL FlagError("Gauss Hermite quadrature type not implemented.",err,error,*999)
        CASE(BASIS_GAUSS_SIMPLEX_QUADRATURE)
          !Allocate one scheme and add it to the list of schemes
          ALLOCATE(NEW_SCHEME,STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate new quadrature scheme.",err,error,*999)
          NEW_SCHEME%QUADRATURE=>BASIS%QUADRATURE
          BASIS%QUADRATURE%NUMBER_OF_SCHEMES=1
          ALLOCATE(NEW_SCHEMES(BASIS%QUADRATURE%NUMBER_OF_SCHEMES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate new quadratures scheme.",err,error,*999)
          NEW_SCHEMES(1)%PTR=>NEW_SCHEME
          BASIS%QUADRATURE%SCHEMES=>NEW_SCHEMES
          !Set up the quadrature scheme map
          ALLOCATE(BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate quadrature scheme map.",err,error,*999)
          DO scheme_idx=1,BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES
            NULLIFY(BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(scheme_idx)%PTR)
          ENDDO !scheme_idx
          BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR=>NEW_SCHEME
          !Set up the gauss point arrays
          CALL Gauss_Simplex(BASIS%QUADRATURE%GAUSS_ORDER,BASIS%NUMBER_OF_XI_COORDINATES,NEW_SCHEME%NUMBER_OF_GAUSS,GSX,GSW, &
            & err,error,*999)
          ALLOCATE(NEW_SCHEME%GAUSS_POSITIONS(BASIS%NUMBER_OF_XI_COORDINATES,NEW_SCHEME%NUMBER_OF_GAUSS),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate Gauss positions.",err,error,*999)
          ALLOCATE(NEW_SCHEME%GAUSS_WEIGHTS(NEW_SCHEME%NUMBER_OF_GAUSS),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate Gauss weights.",err,error,*999)
          ALLOCATE(NEW_SCHEME%GAUSS_BASIS_FNS(BASIS%NUMBER_OF_ELEMENT_PARAMETERS,BASIS%NUMBER_OF_PARTIAL_DERIVATIVES, &
            & NEW_SCHEME%NUMBER_OF_GAUSS),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate Gauss basis functions.",err,error,*999)
          NEW_SCHEME%GAUSS_POSITIONS(1:BASIS%NUMBER_OF_XI_COORDINATES,1:NEW_SCHEME%NUMBER_OF_GAUSS)= &
            & GSX(1:BASIS%NUMBER_OF_XI_COORDINATES,1:NEW_SCHEME%NUMBER_OF_GAUSS)
          NEW_SCHEME%GAUSS_WEIGHTS(1:NEW_SCHEME%NUMBER_OF_GAUSS)=GSW(1:NEW_SCHEME%NUMBER_OF_GAUSS)
          DO ng=1,NEW_SCHEME%NUMBER_OF_GAUSS
            ns=0
            DO nn=1,BASIS%NUMBER_OF_NODES
              DO nk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                ns=ns+1
                DO nu=1,BASIS%NUMBER_OF_PARTIAL_DERIVATIVES
                  SELECT CASE(BASIS%TYPE)
                  CASE(BASIS_SIMPLEX_TYPE)
                    !Gauss positions are in area coordinates so call the simplex basis evaluate directly
                    NEW_SCHEME%GAUSS_BASIS_FNS(ns,nu,ng)= &
                      & BASIS_SIMPLEX_BASIS_EVALUATE(BASIS,nn,nu,NEW_SCHEME%GAUSS_POSITIONS(1:BASIS%NUMBER_OF_XI_COORDINATES,ng), &
                      & err,error)
                    IF(ERR/=0) GOTO 999                        
                  CASE DEFAULT
                    CALL FlagError("Not implemented.",err,error,*999)
                  END SELECT
                ENDDO !nu
              ENDDO !nk
            ENDDO !nn
          ENDDO !ng
          !Create face quadrature scheme, if requested
          IF(BASIS%QUADRATURE%EVALUATE_FACE_GAUSS) THEN
            IF(BASIS%NUMBER_OF_XI==3) THEN
              !Find maximum number of face gauss points and allocate the arrays
              MAX_NUM_FACE_GAUSS=PRODUCT(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(1:BASIS%NUMBER_OF_XI))
              MAX_NUM_FACE_GAUSS=MAX_NUM_FACE_GAUSS/MINVAL(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(1:BASIS%NUMBER_OF_XI))
              ALLOCATE(NEW_SCHEME%NUMBER_OF_FACE_GAUSS(BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate number of face gauss",err,error,*999)
              ALLOCATE(NEW_SCHEME%FACE_GAUSS_POSITIONS(BASIS%NUMBER_OF_XI_COORDINATES,MAX_NUM_FACE_GAUSS, &
                & BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate face Gauss positions",err,error,*999)
              ALLOCATE(NEW_SCHEME%FACE_GAUSS_WEIGHTS(MAX_NUM_FACE_GAUSS,BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate face Gauss weights",err,error,*999)
              ALLOCATE(NEW_SCHEME%FACE_GAUSS_BASIS_FNS(BASIS%NUMBER_OF_ELEMENT_PARAMETERS,BASIS%NUMBER_OF_PARTIAL_DERIVATIVES, &
                & MAX_NUM_FACE_GAUSS,BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate face Gauss basis function values array",err,error,*999)
              !Zero them out just to be safe
              NEW_SCHEME%FACE_GAUSS_POSITIONS=0.0_DP
              NEW_SCHEME%FACE_GAUSS_WEIGHTS=0.0_DP
              NEW_SCHEME%FACE_GAUSS_BASIS_FNS=0.0_DP
              !Populate face_gauss_positions, weights, basis_fn
              DO face_idx=1,BASIS%NUMBER_OF_LOCAL_FACES
                !The number of face xi coordinates will be 3 for triangular face on a tet
                numberOfFaceXiCoordinates = BASIS%NUMBER_OF_XI
                !Set up the gauss point arrays for the face
                CALL Gauss_Simplex(BASIS%QUADRATURE%GAUSS_ORDER,numberOfFaceXiCoordinates, &
                  & NEW_SCHEME%NUMBER_OF_FACE_GAUSS(face_idx),GSX,GSW,err,error,*999)
                IF(ERR/=0) CALL FlagError("Could not allocate Gauss basis functions",err,error,*999)
                NEW_SCHEME%FACE_GAUSS_POSITIONS(1:numberOfFaceXiCoordinates,1:NEW_SCHEME%NUMBER_OF_FACE_GAUSS(face_idx), &
                  & face_idx)=GSX(1:numberOfFaceXiCoordinates,1:NEW_SCHEME%NUMBER_OF_FACE_GAUSS(face_idx))
                NEW_SCHEME%FACE_GAUSS_WEIGHTS(1:NEW_SCHEME%NUMBER_OF_FACE_GAUSS(face_idx),face_idx)= &
                  & GSW(1:NEW_SCHEME%NUMBER_OF_FACE_GAUSS(face_idx))

                DO ng=1,NEW_SCHEME%NUMBER_OF_FACE_GAUSS(face_idx)
                  ns=0
                  DO nn=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(face_idx)
                    DO nk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                      ns=ns+1
                      DO nu=1,BASIS%NUMBER_OF_PARTIAL_DERIVATIVES
                        SELECT CASE(BASIS%TYPE)
                        CASE(BASIS_SIMPLEX_TYPE)
                          NEW_SCHEME%FACE_GAUSS_BASIS_FNS(ns,nu,ng,face_idx)= &
                            & BASIS_SIMPLEX_BASIS_EVALUATE(BASIS,nn,nu, &
                            & NEW_SCHEME%FACE_GAUSS_POSITIONS(1:numberOfFaceXiCoordinates,ng,face_idx),err,error)
                          IF(ERR/=0) GOTO 999                        
                        CASE DEFAULT
                          CALL FlagError("Not implemented",err,error,*999)
                        END SELECT
                      ENDDO !nu
                    ENDDO !nk
                  ENDDO !nn
                ENDDO !ng

              ENDDO !face_idx
            ELSE
              CALL FlagError("Cannot evaluate face quadrature schemes for a non three dimensional element.",err,error,*999)
            ENDIF
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="Quadrature type "//TRIM(NumberToVString(BASIS%QUADRATURE%TYPE,"*",err,error))//" is invalid."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated.",err,error,*998)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Quadrature type = ",BASIS%QUADRATURE%TYPE,err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,BASIS%NUMBER_OF_XI,3,3,BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI, &
        & '("  Number of gauss points(ni):",3(X,I2))','(22X,3(X,I2))',err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of quadrature schemes = ",BASIS%QUADRATURE%NUMBER_OF_SCHEMES, &
        & err,error,*999)
      DO scheme_idx=1,BASIS%QUADRATURE%NUMBER_OF_SCHEMES
        SCHEME=>BASIS%QUADRATURE%SCHEMES(scheme_idx)%PTR
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Scheme = ",scheme_idx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Total number of gauss points = ",SCHEME%NUMBER_OF_GAUSS, &
          & err,error,*999)
        IF(DIAGNOSTICS2) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      Gauss point positions and weights:",err,error,*999)
          DO ng=1,SCHEME%NUMBER_OF_GAUSS
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Gauss point = ",ng,err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,BASIS%NUMBER_OF_XI_COORDINATES,3,3,SCHEME%GAUSS_POSITIONS(:,ng), &
              & '("          position(ni)   :",3(X,F12.4))','(26X,3(X,F12.4))',err,error,*999)
            CALL WriteStringFmtValue(DIAGNOSTIC_OUTPUT_TYPE,"          WEIGHT         : ",SCHEME%GAUSS_WEIGHTS(ng), &
              & "(F12.4)",err,error,*999)
          ENDDO !ng          
        ENDIF
        IF(DIAGNOSTICS3) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"      Basis functions evaluated at Gauss points:",err,error,*999)
          DO ng=1,SCHEME%NUMBER_OF_GAUSS
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Gauss point = ",ng,err,error,*999)
            DO nu=1,BASIS%NUMBER_OF_PARTIAL_DERIVATIVES
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Partial derivative number = ",nu,err,error,*999)
              CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,BASIS%NUMBER_OF_ELEMENT_PARAMETERS,4,4, &
                & SCHEME%GAUSS_BASIS_FNS(:,nu,ng),'("          BASIS FNS(ns)  :",4(X,F12.4))','(26X,4(X,F12.4))',err,error,*999)
            ENDDO !nu
          ENDDO !ng
        ENDIF
      ENDDO !scheme_idx
    ENDIF
    
    EXITS("BASIS_QUADRATURE_CREATE")
    RETURN
999 IF(ASSOCIATED(NEW_SCHEME)) THEN
      IF(ALLOCATED(NEW_SCHEME%GAUSS_POSITIONS)) DEALLOCATE(NEW_SCHEME%GAUSS_POSITIONS)
      IF(ALLOCATED(NEW_SCHEME%GAUSS_WEIGHTS)) DEALLOCATE(NEW_SCHEME%GAUSS_WEIGHTS)
      IF(ALLOCATED(NEW_SCHEME%GAUSS_BASIS_FNS)) DEALLOCATE(NEW_SCHEME%GAUSS_BASIS_FNS)
      DEALLOCATE(NEW_SCHEME)     
    ENDIF
    IF(ALLOCATED(WEIGHTS)) DEALLOCATE(WEIGHTS)
    IF(ALLOCATED(POSITIONS)) DEALLOCATE(POSITIONS)
    IF(ALLOCATED(POSITIONS_MATRIX)) DEALLOCATE(POSITIONS_MATRIX)
    IF(ASSOCIATED(NEW_SCHEMES)) DEALLOCATE(NEW_SCHEMES)
    NULLIFY(BASIS%QUADRATURE%SCHEMES)
998 ERRORSEXITS("BASIS_QUADRATURE_CREATE",err,error)    
    RETURN 1
   
  END SUBROUTINE BASIS_QUADRATURE_CREATE
        
  !
  !================================================================================================================================
  !
  
  !>Destroys a quadrature on a given basis and deallocates all memory. \todo fix all this basis/quadrature into standard form.
  SUBROUTINE BASIS_QUADRATURE_DESTROY(QUADRATURE,err,error,*)

    !Argument variables
    TYPE(QUADRATURE_TYPE), POINTER :: QUADRATURE !<A pointer to the quadrature
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("BASIS_QUADRATURE_DESTROY",err,error,*999)

    IF(ASSOCIATED(QUADRATURE)) THEN
      CALL FlagError("Not implemented.",err,error,*999)     
    ELSE
      CALL FlagError("Basis is not associated",err,error,*999)
    ENDIF
    
    EXITS("BASIS_QUADRATURE_DESTROY")
    RETURN
999 ERRORSEXITS("BASIS_QUADRATURE_DESTROY",err,error)
    RETURN 1
  END SUBROUTINE BASIS_QUADRATURE_DESTROY

  !
  !================================================================================================================================
  !
    
  !>Finalises a quadrature on a given basis and deallocates all memory
  SUBROUTINE BASIS_QUADRATURE_FINALISE(BASIS,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: scheme_idx
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: SCHEME

    ENTERS("BASIS_QUADRATURE_FINALISE",err,error,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(ASSOCIATED(BASIS%QUADRATURE%BASIS)) THEN
        DO scheme_idx=1,BASIS%QUADRATURE%NUMBER_OF_SCHEMES
          SCHEME=>BASIS%QUADRATURE%SCHEMES(scheme_idx)%PTR
          !Destroy all scheme components
          IF (ASSOCIATED(SCHEME)) THEN
            IF(ALLOCATED(SCHEME%GAUSS_POSITIONS)) DEALLOCATE(SCHEME%GAUSS_POSITIONS)
            IF(ALLOCATED(SCHEME%GAUSS_WEIGHTS)) DEALLOCATE(SCHEME%GAUSS_WEIGHTS)
            IF(ALLOCATED(SCHEME%GAUSS_BASIS_FNS)) DEALLOCATE(SCHEME%GAUSS_BASIS_FNS)
            DEALLOCATE(SCHEME)
          ENDIF
        ENDDO !scheme_idx
        IF(ASSOCIATED(BASIS%QUADRATURE%SCHEMES)) DEALLOCATE(BASIS%QUADRATURE%SCHEMES)
        BASIS%QUADRATURE%NUMBER_OF_SCHEMES=0
        IF(ALLOCATED(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI)) DEALLOCATE(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI)
        NULLIFY(BASIS%QUADRATURE%BASIS)
        BASIS%QUADRATURE%TYPE=-1
        IF(ALLOCATED(BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP)) DEALLOCATE(BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP)
      ELSE
        CALL FlagError("Basis quadrature basis is not associated",err,error,*999)
      ENDIF      
    ELSE
      CALL FlagError("Basis is not associated",err,error,*999)
    ENDIF
    
    EXITS("BASIS_QUADRATURE_FINALISE")
    RETURN
999 ERRORSEXITS("BASIS_QUADRATURE_FINALISE",err,error)
    RETURN 1
  END SUBROUTINE BASIS_QUADRATURE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a quadrature on the given basis.
  SUBROUTINE BASIS_QUADRATURE_INITIALISE(BASIS,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ni
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("BASIS_QUADRATURE_INITIALISE",err,error,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(ASSOCIATED(BASIS%QUADRATURE%BASIS)) THEN
        LOCAL_ERROR="Basis number "//TRIM(NumberToVString(BASIS%USER_NUMBER,"*",err,error))// &
          & " already has a quadrature associated"
        CALL FlagError(LOCAL_ERROR,err,error,*998)
      ELSE
        SELECT CASE(BASIS%TYPE)
        CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
          !Set up a default Gauss Legendre quadrature 
          BASIS%QUADRATURE%NUMBER_OF_SCHEMES=0
          NULLIFY(BASIS%QUADRATURE%SCHEMES)
          BASIS%QUADRATURE%BASIS=>BASIS        
          BASIS%QUADRATURE%TYPE=BASIS_GAUSS_LEGENDRE_QUADRATURE
          !Set up a default number of Gauss points appropriate for the given interpolation order in each direction.
          ALLOCATE(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(BASIS%NUMBER_OF_XI),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate number of Gauss in each xi direction",err,error,*999)
          DO ni=1,BASIS%NUMBER_OF_XI
            SELECT CASE(BASIS%INTERPOLATION_XI(ni))
            CASE(BASIS_LINEAR_LAGRANGE_INTERPOLATION)
              BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)=2
            CASE(BASIS_QUADRATIC_LAGRANGE_INTERPOLATION)
              BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)=3
            CASE(BASIS_CUBIC_LAGRANGE_INTERPOLATION)
              BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)=4
            CASE(BASIS_QUADRATIC1_HERMITE_INTERPOLATION,BASIS_QUADRATIC2_HERMITE_INTERPOLATION)
              BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)=3
            CASE(BASIS_CUBIC_HERMITE_INTERPOLATION)
              BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)=4
            CASE DEFAULT
              LOCAL_ERROR="Interpolation xi value "//TRIM(NumberToVString(BASIS%INTERPOLATION_XI(ni),"*",err,error))// &
                & " in xi direction "//TRIM(NumberToVString(ni,"*",err,error))//" is invalid"
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          ENDDO !ni
        CASE(BASIS_SIMPLEX_TYPE)
         !Set up a default quadrature 
          BASIS%QUADRATURE%NUMBER_OF_SCHEMES=0
          NULLIFY(BASIS%QUADRATURE%SCHEMES)
          BASIS%QUADRATURE%BASIS=>BASIS        
          BASIS%QUADRATURE%TYPE=BASIS_GAUSS_SIMPLEX_QUADRATURE
          !Set up a default order appropriate for the given interpolation.
          ALLOCATE(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(BASIS%NUMBER_OF_XI),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate number of Gauss in each xi direction",err,error,*999)
!!TODO: \todo Set these to something more meaningfull!
          SELECT CASE(BASIS%INTERPOLATION_XI(1))
          CASE(BASIS_LINEAR_SIMPLEX_INTERPOLATION)
            SELECT CASE(BASIS%NUMBER_OF_XI)
            CASE(1)
              BASIS%QUADRATURE%GAUSS_ORDER=2
            CASE(2)
              BASIS%QUADRATURE%GAUSS_ORDER=3
            CASE(3)
              BASIS%QUADRATURE%GAUSS_ORDER=3
            CASE DEFAULT
              LOCAL_ERROR="The number of xi directions ("//TRIM(NumberToVString(BASIS%NUMBER_OF_XI,"*",err,error))// &
                & ") is invalid"
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          CASE(BASIS_QUADRATIC_SIMPLEX_INTERPOLATION)
            SELECT CASE(BASIS%NUMBER_OF_XI)
            CASE(1)
              BASIS%QUADRATURE%GAUSS_ORDER=3
            CASE(2)
              BASIS%QUADRATURE%GAUSS_ORDER=3
            CASE(3)
              BASIS%QUADRATURE%GAUSS_ORDER=5
            CASE DEFAULT
              LOCAL_ERROR="The number of xi directions ("//TRIM(NumberToVString(BASIS%NUMBER_OF_XI,"*",err,error))// &
                & ") is invalid"
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          CASE(BASIS_CUBIC_SIMPLEX_INTERPOLATION)
            SELECT CASE(BASIS%NUMBER_OF_XI)
            CASE(1)
              BASIS%QUADRATURE%GAUSS_ORDER=3
            CASE(2)
              BASIS%QUADRATURE%GAUSS_ORDER=5
            CASE(3)
              BASIS%QUADRATURE%GAUSS_ORDER=5
            CASE DEFAULT
              LOCAL_ERROR="The number of xi directions ("//TRIM(NumberToVString(BASIS%NUMBER_OF_XI,"*",err,error))// &
                & ") is invalid"
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="Interpolation xi value "//TRIM(NumberToVString(BASIS%INTERPOLATION_XI(1),"*",err,error))// &
              & " in xi direction 1 is invalid"
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
          BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI=BASIS%QUADRATURE%GAUSS_ORDER
        CASE DEFAULT
          LOCAL_ERROR="Basis type value "//TRIM(NumberToVString(BASIS%INTERPOLATION_XI(ni),"*",err,error))// &
            & " is invalid or not implemented"
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated",err,error,*998)
    ENDIF
    
    EXITS("BASIS_QUADRATURE_INITIALISE")
    RETURN
999 IF(ALLOCATED(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI)) DEALLOCATE(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI)
998 ERRORSEXITS("BASIS_QUADRATURE_INITIALISE",err,error)    
    RETURN 1
    
  END SUBROUTINE BASIS_QUADRATURE_INITIALISE

  !
  !================================================================================================================================
  !
  
  !>Get the number of Gauss points in each xi direction on a basis quadrature identified by a pointer. \see OpenCMISS::Iron::cmfe_BasisQuadratureNumberOfGaussXiGet
  SUBROUTINE BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_GET(BASIS,QUADRATURE_NUMBER_OF_GAUSS_XI,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: QUADRATURE_NUMBER_OF_GAUSS_XI(:) !<On return, the number of Gauss in each Xi direction
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_GET",err,error,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        IF(ASSOCIATED(BASIS%QUADRATURE%BASIS)) THEN
          IF(SIZE(QUADRATURE_NUMBER_OF_GAUSS_XI,1)>=SIZE(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI,1)) THEN
            QUADRATURE_NUMBER_OF_GAUSS_XI=BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI
          ELSE
            LOCAL_ERROR="The size of QUADRATURE_NUMBER_OF_GAUSS_XI is too small. The supplied size is "// &
              & TRIM(NumberToVString(SIZE(QUADRATURE_NUMBER_OF_GAUSS_XI,1),"*",err,error))//" and it needs to be >= "// &
              & TRIM(NumberToVString(SIZE(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI,1),"*",err,error))//"."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Quadrature basis is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Basis has not finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated.",err,error,*999)
    ENDIF
      
    EXITS("BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_GET")
    RETURN
999 ERRORSEXITS("BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_GET",err,error)
    RETURN
  END SUBROUTINE BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_GET
    
  !
  !================================================================================================================================
  !
  !================================================================================================================================
  !

  !>Sets/changes the number of Gauss points in each xi direction on a basis quadrature identified by a pointer. \see OpenCMISS::Iron::cmfe_BasisQuadratureNumberOfGaussSet
  SUBROUTINE BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET(BASIS,NUMBER_OF_GAUSS_XI,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_GAUSS_XI(:) !<The number of Gauss in each Xi direction
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ni
    TYPE(VARYING_STRING) :: LOCAL_ERROR,LOCAL_WARNING
    
    ENTERS("BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET",err,error,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        CALL FlagError("Basis has been finished.",err,error,*999)
      ELSE
        IF(ASSOCIATED(BASIS%QUADRATURE%BASIS)) THEN          
          IF(SIZE(NUMBER_OF_GAUSS_XI,1)==BASIS%NUMBER_OF_XI) THEN
            IF(ANY(NUMBER_OF_GAUSS_XI<1)) CALL FlagError("Invalid number of gauss values.",err,error,*999)
            BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(1:BASIS%NUMBER_OF_XI)=NUMBER_OF_GAUSS_XI(1:BASIS%NUMBER_OF_XI)
            !Check the number of gauss points is sufficient for the interpolation order and flag a warning if not
            DO ni=1,BASIS%NUMBER_OF_XI
              SELECT CASE(BASIS%INTERPOLATION_XI(ni))
              CASE(BASIS_LINEAR_LAGRANGE_INTERPOLATION)
                IF(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)<2) THEN
                  LOCAL_WARNING=TRIM(NumberToVString(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni),"*",err,error))// &
                    & " Gauss points are insufficient for linear Lagrange interpolation."
                  CALL FLAG_WARNING(LOCAL_WARNING,err,error,*999)
                ENDIF
              CASE(BASIS_QUADRATIC_LAGRANGE_INTERPOLATION)
                IF(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)<2) THEN
                  LOCAL_WARNING=TRIM(NumberToVString(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni),"*",err,error))//&
                    & " Gauss points are insufficient for quadratic Lagrange interpolation."
                  CALL FLAG_WARNING(LOCAL_WARNING,err,error,*999)
                ENDIF
              CASE(BASIS_CUBIC_LAGRANGE_INTERPOLATION)
                IF(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)<3) THEN
                  LOCAL_WARNING=TRIM(NumberToVString(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni),"*",err,error))//&
                    & " Gauss points are insufficient for cubic Lagrange interpolation."
                  CALL FLAG_WARNING(LOCAL_WARNING,err,error,*999)
                ENDIF
               CASE(BASIS_QUADRATIC1_HERMITE_INTERPOLATION,BASIS_QUADRATIC2_HERMITE_INTERPOLATION)
                IF(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)<2) THEN
                  LOCAL_WARNING=TRIM(NumberToVString(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni),"*",err,error))//&
                    & " Gauss points are insufficient for quadratic Hermite interpolation."
                  CALL FLAG_WARNING(LOCAL_WARNING,err,error,*999)
                ENDIF
              CASE(BASIS_CUBIC_HERMITE_INTERPOLATION)
                IF(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni)<3) THEN
                  LOCAL_WARNING=TRIM(NumberToVString(BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI(ni),"*",err,error))//&
                    & " Gauss points are insufficient for cubic Hermite interpolation."
                  CALL FLAG_WARNING(LOCAL_WARNING,err,error,*999)
                ENDIF
              CASE(BASIS_LINEAR_SIMPLEX_INTERPOLATION)
                LOCAL_WARNING="For simplex elements please set quadrature order rather than number of gauss points."
                CALL FLAG_WARNING(LOCAL_WARNING,err,error,*999)
              CASE(BASIS_QUADRATIC_SIMPLEX_INTERPOLATION)
                LOCAL_WARNING="For simplex elements please set quadrature order rather than number of gauss points."
                CALL FLAG_WARNING(LOCAL_WARNING,err,error,*999)
              CASE DEFAULT
                LOCAL_ERROR="Interpolation xi value "//TRIM(NumberToVString(BASIS%INTERPOLATION_XI(ni),"*",err,error))// &
                  & " is invalid for xi direction "//TRIM(NumberToVString(ni,"*",err,error))//"."
                CALL FlagError(LOCAL_ERROR,err,error,*999)
              END SELECT
            ENDDO !xi
          ELSE
            LOCAL_ERROR="The size of the number of Gauss array ("// &
              & TRIM(NumberToVString(SIZE(NUMBER_OF_GAUSS_XI,1),"*",err,error))// &
              & ") does not match the number of xi directions ("// &
              & TRIM(NumberToVString(BASIS%NUMBER_OF_XI,"*",err,error))//") for basis number "// &
              & TRIM(NumberToVString(BASIS%USER_NUMBER,"*",err,error))//"."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Quadrature basis is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated.",err,error,*999)
    ENDIF
      
    EXITS("BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET")
    RETURN
999 ERRORSEXITS("BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET",err,error)
    RETURN 1
  END SUBROUTINE BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET
  !
  !================================================================================================================================
  !
  
  !>Returns the xi positions of a Gauss point on a basis quadrature identified by a pointer. \see OpenCMISS::Iron::cmfe_BasisQuadratureGaussXiGet
  SUBROUTINE BASIS_QUADRATURE_SINGLE_GAUSS_XI_GET(BASIS,SCHEME,GAUSS_POINT,GAUSS_XI,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: SCHEME !<The quadrature scheme to return the Gauss points for
    INTEGER(INTG), INTENT(IN) :: GAUSS_POINT !<The Gauss point to return the element xi position for.
    REAL(DP), INTENT(OUT) :: GAUSS_XI(:) !<On return, GAUSS_XI(xi_direction) the xi position of the specified Gauss point for the specified quadrature scheme.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("BASIS_QUADRATURE_SINGLE_GAUSS_XI_GET",err,error,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        IF(ASSOCIATED(BASIS%QUADRATURE%BASIS)) THEN
          QUADRATURE_SCHEME=>BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(SCHEME)%PTR
          IF(ASSOCIATED(QUADRATURE_SCHEME)) THEN
            IF(SIZE(GAUSS_XI)==BASIS%NUMBER_OF_XI) THEN
              IF(GAUSS_POINT>0.AND.GAUSS_POINT<=QUADRATURE_SCHEME%NUMBER_OF_GAUSS) THEN
                GAUSS_XI(:)=QUADRATURE_SCHEME%GAUSS_POSITIONS(:,GAUSS_POINT)
              ELSE
                LOCAL_ERROR="The specified Gauss point number of "// & 
                  & TRIM(NumberToVString(GAUSS_POINT,"*",err,error))//" is invalid for the specified "// &
                  & "quadrature scheme of the specified element for this field which has "// &
                  & TRIM(NumberToVString(QUADRATURE_SCHEME%NUMBER_OF_GAUSS,"*",err,error))//" Gauss points."
                CALL FlagError(LOCAL_ERROR,err,error,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The number of xi values to return is invalid and needs to be "// &
                & TRIM(NumberToVString(BASIS%NUMBER_OF_XI,"*",err,error))//" for the specified basis."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("The specified quadrature scheme is not associated for this basis.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Quadrature basis is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Basis has not finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated.",err,error,*999)
    ENDIF
      
    EXITS("BASIS_QUADRATURE_SINGLE_GAUSS_XI_GET")
    RETURN
999 ERRORSEXITS("BASIS_QUADRATURE_SINGLE_GAUSS_XI_GET",err,error)
    RETURN
  END SUBROUTINE BASIS_QUADRATURE_SINGLE_GAUSS_XI_GET

  !
  !================================================================================================================================
  !
  
  !>Returns the xi positions of Gauss points on a basis quadrature identified by a pointer. If no Gauss points are specified then xi positions of all Gauss points are returned. \see OpenCMISS::Iron::cmfe_BasisQuadratureGaussXiGet
  SUBROUTINE BASIS_QUADRATURE_MULTIPLE_GAUSS_XI_GET(BASIS,SCHEME,GAUSS_POINTS,GAUSS_XI,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: SCHEME !<The quadrature scheme to return the Gauss points for
    INTEGER(INTG), INTENT(IN) :: GAUSS_POINTS(:) !<The Gauss points to return the element xi positions for.
    REAL(DP), INTENT(OUT) :: GAUSS_XI(:,:) !<On return, GAUSS_XI(xi_direction,Gauss_point) the Gauss xi positions for the specified quadrature scheme.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: Gauss_point
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("BASIS_QUADRATURE_MULTIPLE_GAUSS_XI_GET",err,error,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        IF(ASSOCIATED(BASIS%QUADRATURE%BASIS)) THEN
          QUADRATURE_SCHEME=>BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(SCHEME)%PTR
          IF(ASSOCIATED(QUADRATURE_SCHEME)) THEN
            IF(SIZE(GAUSS_XI,1)==BASIS%NUMBER_OF_XI) THEN
              IF(SIZE(GAUSS_POINTS)==0) THEN !Return all Gauss point xi locations.
                IF(SIZE(GAUSS_XI,2)==QUADRATURE_SCHEME%NUMBER_OF_GAUSS) THEN
                  GAUSS_XI=QUADRATURE_SCHEME%GAUSS_POSITIONS
                ELSE
                  LOCAL_ERROR="The number of Gauss Points to return the xi values for is invalid and needs to be "// &
                    & TRIM(NumberToVString(QUADRATURE_SCHEME%NUMBER_OF_GAUSS,"*",err,error))//"."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                ENDIF
              ELSE !Return only specified Gauss point xi locations.
                DO Gauss_point=1,SIZE(GAUSS_POINTS)
                  IF(GAUSS_POINTS(Gauss_point)>0.AND.GAUSS_POINTS(Gauss_point)<=QUADRATURE_SCHEME%NUMBER_OF_GAUSS) THEN
                    GAUSS_XI(:,Gauss_point)=QUADRATURE_SCHEME%GAUSS_POSITIONS(:,GAUSS_POINTS(Gauss_point))
                  ELSE
                    LOCAL_ERROR="The specified Gauss point number of "// & 
                      & TRIM(NumberToVString(GAUSS_POINTS(Gauss_point),"*",err,error))//" is invalid for the specified "// &
                      & "quadrature scheme of the specified element for this field which has "// &
                      & TRIM(NumberToVString(QUADRATURE_SCHEME%NUMBER_OF_GAUSS,"*",err,error))//" Gauss points."
                    CALL FlagError(LOCAL_ERROR,err,error,*999)
                  ENDIF
                ENDDO
              ENDIF
            ELSE
              LOCAL_ERROR="The number of xi values to return is invalid and needs to be "// &
                & TRIM(NumberToVString(BASIS%NUMBER_OF_XI,"*",err,error))//" for the specified basis."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("The specified quadrature scheme is not associated for this basis.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Quadrature basis is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Basis has not finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated.",err,error,*999)
    ENDIF
      
    EXITS("BASIS_QUADRATURE_MULTIPLE_GAUSS_XI_GET")
    RETURN
999 ERRORSEXITS("BASIS_QUADRATURE_MULTIPLE_GAUSS_XI_GET",err,error)
    RETURN
  END SUBROUTINE BASIS_QUADRATURE_MULTIPLE_GAUSS_XI_GET

  !
  !================================================================================================================================
  !

  !>Get the order of a quadrature for a basis quadrature identified by a pointer. \see OpenCMISS::Iron::cmfe_BasisQuadratureOrderGet
  SUBROUTINE BASIS_QUADRATURE_ORDER_GET(BASIS,QUADRATURE_ORDER,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: QUADRATURE_ORDER !<On return, the quadrature order for the specified basis.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("BASIS_QUADRATURE_ORDER_GET",err,error,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        IF(ASSOCIATED(BASIS%QUADRATURE%BASIS)) THEN
          QUADRATURE_ORDER=BASIS%QUADRATURE%GAUSS_ORDER
        ELSE
          CALL FlagError("Quadrature basis is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Basis has not finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated.",err,error,*999)
    ENDIF
      
    EXITS("BASIS_QUADRATURE_ORDER_GET")
    RETURN
999 ERRORSEXITS("BASIS_QUADRATURE_ORDER_GET",err,error)
    RETURN
  END SUBROUTINE BASIS_QUADRATURE_ORDER_GET

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the order of a quadrature for a basis quadrature identified by a user number.
  SUBROUTINE BASIS_QUADRATURE_ORDER_SET_NUMBER(USER_NUMBER,ORDER,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the basis
    INTEGER(INTG), INTENT(IN) :: ORDER !<The quadrature order to be set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    
    ENTERS("BASIS_QUADRATURE_ORDER_SET_NUMBER",err,error,*999)

    CALL Basis_UserNumberFind(USER_NUMBER,BASIS,err,error,*999)
    CALL BASIS_QUADRATURE_ORDER_SET(BASIS,ORDER,err,error,*999)
    
    EXITS("BASIS_QUADRATURE_ORDER_SET_NUMBER")
    RETURN
999 ERRORSEXITS("BASIS_QUADRATURE_ORDER_SET_NUMBER",err,error)
    RETURN 1
  END SUBROUTINE BASIS_QUADRATURE_ORDER_SET_NUMBER

  !
  !================================================================================================================================
  !

  !>Sets/changes the order of a quadrature for a basis quadrature identified by a pointer. \see OpenCMISS::Iron::cmfe_BasisQuadratureOrderSet
  SUBROUTINE BASIS_QUADRATURE_ORDER_SET_PTR(BASIS,ORDER,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: ORDER !<The quadrature order to be set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("BASIS_QUADRATURE_ORDER_SET_PTR",err,error,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        CALL FlagError("Basis has been finished",err,error,*999)
      ELSE
        IF(ASSOCIATED(BASIS%QUADRATURE%BASIS)) THEN
          IF(BASIS%TYPE==BASIS_SIMPLEX_TYPE) THEN !Relax this i.e., use this to set gauss points in each direction for LHTP's???
            IF(ORDER>1.AND.ORDER<=5) THEN
              BASIS%QUADRATURE%GAUSS_ORDER=ORDER
            ELSE
              LOCAL_ERROR="An order value of "//TRIM(NumberToVString(ORDER,"*",err,error))// &
                & " is invalid. You must specify and order between 1 and 5."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Can only set the quadrature order for simplex basis types.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Quadrature basis is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated.",err,error,*999)
    ENDIF
      
    EXITS("BASIS_QUADRATURE_ORDER_SET_PTR")
    RETURN
999 ERRORSEXITS("BASIS_QUADRATURE_ORDER_SET_PTR",err,error)
    RETURN 1
  END SUBROUTINE BASIS_QUADRATURE_ORDER_SET_PTR

  !
  !================================================================================================================================
  !
  
  !>get the quadrature type on a basis identified by a pointer. \see OpenCMISS::Iron::cmfe_BasisQuadratureTypeGet
  SUBROUTINE BASIS_QUADRATURE_TYPE_GET(BASIS,QUADRATURE_TYPE,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: QUADRATURE_TYPE !<On return, the quadrature type to be get \see BASIS_ROUTINES_QuadratureTypes
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("BASIS_QUADRATURE_TYPE_GET",err,error,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        IF(ASSOCIATED(BASIS%QUADRATURE%BASIS)) THEN
          QUADRATURE_TYPE=BASIS%QUADRATURE%TYPE
        ELSE
          CALL FlagError("Basis quadrature basis is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Basis has not finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated.",err,error,*999)
    ENDIF
    
    EXITS("BASIS_QUADRATURE_TYPE_GET")
    RETURN
999 ERRORSEXITS("BASIS_QUADRATURE_TYPE_GET",err,error)
    RETURN
  END SUBROUTINE BASIS_QUADRATURE_TYPE_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the quadrature type for a basis quadrature identified by a user number.
  SUBROUTINE BASIS_QUADRATURE_TYPE_SET_NUMBER(USER_NUMBER,TYPE,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the basis
    INTEGER(INTG), INTENT(IN) :: TYPE !<The quadrature type to be set \see BASIS_ROUTINES_QuadratureTypes
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    
    ENTERS("BASIS_QUADRATURE_TYPE_SET_NUMBER",err,error,*999)

    CALL Basis_UserNumberFind(USER_NUMBER,BASIS,err,error,*999)
    CALL BASIS_QUADRATURE_TYPE_SET_PTR(BASIS,TYPE,err,error,*999)
    
    EXITS("BASIS_QUADRATURE_TYPE_SET_NUMBER")
    RETURN
999 ERRORSEXITS("BASIS_QUADRATURE_TYPE_SET_NUMBER",err,error)
    RETURN 1
  END SUBROUTINE BASIS_QUADRATURE_TYPE_SET_NUMBER

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the quadrature type on a basis identified by a pointer.
  SUBROUTINE BASIS_QUADRATURE_TYPE_SET_PTR(BASIS,TYPE,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: TYPE !<The quadrature type to be set \see BASIS_ROUTINES_QuadratureTypes
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("BASIS_QUADRATURE_TYPE_SET_PTR",err,error,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        CALL FlagError("Basis has been finished.",err,error,*999)
      ELSE
        IF(ASSOCIATED(BASIS%QUADRATURE%BASIS)) THEN
          SELECT CASE(TYPE)
          CASE(BASIS_GAUSS_LEGENDRE_QUADRATURE)
            BASIS%QUADRATURE%TYPE=BASIS_GAUSS_LEGENDRE_QUADRATURE
          CASE(BASIS_GAUSS_LAGUERRE_QUADRATURE)
            BASIS%QUADRATURE%TYPE=BASIS_GAUSS_LAGUERRE_QUADRATURE
            CALL FlagError("Gauss Laguerre quadrature is not implemented.",err,error,*999)
          CASE(BASIS_GUASS_HERMITE_QUADRATURE)
            BASIS%QUADRATURE%TYPE=BASIS_GUASS_HERMITE_QUADRATURE
            CALL FlagError("Gauss Hermite quadrature is not implemented.",err,error,*999)
          CASE DEFAULT
            LOCAL_ERROR="Quadrature type "//TRIM(NumberToVString(TYPE,"*",err,error))//" is invalid."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Basis quadrature basis is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated.",err,error,*999)
    ENDIF
    
    EXITS("BASIS_QUADRATURE_TYPE_SET_PTR")
    RETURN
999 ERRORSEXITS("BASIS_QUADRATURE_TYPE_SET_PTR",err,error)
    RETURN 1
  END SUBROUTINE BASIS_QUADRATURE_TYPE_SET_PTR

  !
  !================================================================================================================================
  !

  !>Sets/changes the local face Gauss evaluation flag on a basis
  SUBROUTINE Basis_QuadratureLocalFaceGaussEvaluateSet(BASIS,FACE_GAUSS_EVALUATE,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    LOGICAL, INTENT(IN) :: FACE_GAUSS_EVALUATE !<face Gauss evaluation flag
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("Basis_QuadratureLocalFaceGaussEvaluateSet",err,error,*999)
    
    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        CALL FlagError("Basis has been finished.",err,error,*999)
      ELSE
        BASIS%QUADRATURE%EVALUATE_FACE_GAUSS=FACE_GAUSS_EVALUATE
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated.",err,error,*999)
    ENDIF
    
    EXITS("Basis_QuadratureLocalFaceGaussEvaluateSet")
    RETURN
999 ERRORSEXITS("Basis_QuadratureLocalFaceGaussEvaluateSet",err,error)
    RETURN 1

  END SUBROUTINE Basis_QuadratureLocalFaceGaussEvaluateSet

  !
  !================================================================================================================================
  !
  
  !>Creates and initialises a simplex basis that has already been allocated BASIS_CREATE_START
  !>\see BASIS_ROUTINES::BASIS_CREATE_START
  SUBROUTINE Basis_SimplexBasisCreate(BASIS,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: MAX_NUM_NODES,ni,nn,ns
    INTEGER(INTG), ALLOCATABLE :: NODES_IN_FACE(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("Basis_SimplexBasisCreate",err,error,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%TYPE==BASIS_SIMPLEX_TYPE) THEN
        BASIS%NUMBER_OF_XI_COORDINATES=BASIS%NUMBER_OF_XI+1 !Simplex bases have an additional area coordinate
        ALLOCATE(BASIS%INTERPOLATION_TYPE(BASIS%NUMBER_OF_XI_COORDINATES),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate INTERPOLATION_TYPE array.",err,error,*999)
        ALLOCATE(BASIS%INTERPOLATION_ORDER(BASIS%NUMBER_OF_XI_COORDINATES),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate INTERPOLATION_ORDER array.",err,error,*999)
        ALLOCATE(BASIS%NUMBER_OF_NODES_XIC(BASIS%NUMBER_OF_XI_COORDINATES),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate NUMBER_OF_NODES_XIC array.",err,error,*999)
        BASIS%DEGENERATE=.FALSE.
        BASIS%NUMBER_OF_COLLAPSED_XI=0
        SELECT CASE(BASIS%NUMBER_OF_XI)
        CASE(1)
          BASIS%NUMBER_OF_PARTIAL_DERIVATIVES=3
          SELECT CASE(BASIS%INTERPOLATION_XI(1))
          CASE(BASIS_LINEAR_SIMPLEX_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(1)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(1)=BASIS_LINEAR_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(2)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(2)=BASIS_LINEAR_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(1)=2            
            BASIS%NUMBER_OF_NODES_XIC(2)=2            
            MAX_NUM_NODES=2
            BASIS%NUMBER_OF_NODES=2
          CASE(BASIS_QUADRATIC_SIMPLEX_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(1)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(1)=BASIS_QUADRATIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(2)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(2)=BASIS_QUADRATIC_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(1)=3
            BASIS%NUMBER_OF_NODES_XIC(2)=3
            MAX_NUM_NODES=3
            BASIS%NUMBER_OF_NODES=3
          CASE(BASIS_CUBIC_SIMPLEX_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(1)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(1)=BASIS_CUBIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(2)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(2)=BASIS_CUBIC_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(1)=4
            BASIS%NUMBER_OF_NODES_XIC(2)=4
            MAX_NUM_NODES=4
            BASIS%NUMBER_OF_NODES=4
          CASE DEFAULT 
            CALL FlagError("Invalid interpolation type.",err,error,*999)
          END SELECT
        CASE(2)
          BASIS%NUMBER_OF_PARTIAL_DERIVATIVES=6
          SELECT CASE(BASIS%INTERPOLATION_XI(2))
          CASE(BASIS_LINEAR_SIMPLEX_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(1)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(1)=BASIS_LINEAR_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(2)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(2)=BASIS_LINEAR_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(3)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(3)=BASIS_LINEAR_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(1)=2            
            BASIS%NUMBER_OF_NODES_XIC(2)=2            
            BASIS%NUMBER_OF_NODES_XIC(3)=2            
            MAX_NUM_NODES=2
            BASIS%NUMBER_OF_NODES=3
          CASE(BASIS_QUADRATIC_SIMPLEX_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(1)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(1)=BASIS_QUADRATIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(2)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(2)=BASIS_QUADRATIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(3)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(3)=BASIS_QUADRATIC_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(1)=3
            BASIS%NUMBER_OF_NODES_XIC(2)=3
            BASIS%NUMBER_OF_NODES_XIC(3)=3
            MAX_NUM_NODES=3
            BASIS%NUMBER_OF_NODES=6
          CASE(BASIS_CUBIC_SIMPLEX_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(1)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(1)=BASIS_CUBIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(2)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(2)=BASIS_CUBIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(3)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(3)=BASIS_CUBIC_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(1)=4
            BASIS%NUMBER_OF_NODES_XIC(2)=4
            BASIS%NUMBER_OF_NODES_XIC(3)=4
            MAX_NUM_NODES=4
            BASIS%NUMBER_OF_NODES=10
          CASE DEFAULT 
            CALL FlagError("Invalid interpolation type.",err,error,*999)
          END SELECT
        CASE(3)
          BASIS%NUMBER_OF_PARTIAL_DERIVATIVES=11
          SELECT CASE(BASIS%INTERPOLATION_XI(3))
          CASE(BASIS_LINEAR_SIMPLEX_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(1)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(1)=BASIS_LINEAR_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(2)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(2)=BASIS_LINEAR_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(3)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(3)=BASIS_LINEAR_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(4)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(4)=BASIS_LINEAR_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(1)=2            
            BASIS%NUMBER_OF_NODES_XIC(2)=2            
            BASIS%NUMBER_OF_NODES_XIC(3)=2            
            BASIS%NUMBER_OF_NODES_XIC(4)=2            
            MAX_NUM_NODES=2
            BASIS%NUMBER_OF_NODES=4
          CASE(BASIS_QUADRATIC_SIMPLEX_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(1)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(1)=BASIS_QUADRATIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(2)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(2)=BASIS_QUADRATIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(3)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(3)=BASIS_QUADRATIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(4)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(4)=BASIS_QUADRATIC_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(1)=3
            BASIS%NUMBER_OF_NODES_XIC(2)=3
            BASIS%NUMBER_OF_NODES_XIC(3)=3
            BASIS%NUMBER_OF_NODES_XIC(4)=3
            MAX_NUM_NODES=3
            BASIS%NUMBER_OF_NODES=10
          CASE(BASIS_CUBIC_SIMPLEX_INTERPOLATION)
            BASIS%INTERPOLATION_TYPE(1)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(1)=BASIS_CUBIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(2)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(2)=BASIS_CUBIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(3)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(3)=BASIS_CUBIC_INTERPOLATION_ORDER
            BASIS%INTERPOLATION_TYPE(4)=BASIS_SIMPLEX_INTERPOLATION
            BASIS%INTERPOLATION_ORDER(4)=BASIS_CUBIC_INTERPOLATION_ORDER
            BASIS%NUMBER_OF_NODES_XIC(1)=4
            BASIS%NUMBER_OF_NODES_XIC(2)=4
            BASIS%NUMBER_OF_NODES_XIC(3)=4
            BASIS%NUMBER_OF_NODES_XIC(4)=4
            MAX_NUM_NODES=4
            BASIS%NUMBER_OF_NODES=20
          CASE DEFAULT 
            CALL FlagError("Invalid interpolation type.",err,error,*999)
          END SELECT
        CASE DEFAULT
          CALL FlagError("Invalid number of xi directions.",err,error,*999)
        END SELECT
        
        ALLOCATE(BASIS%NODE_AT_COLLAPSE(BASIS%NUMBER_OF_NODES),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate node at collapse.",err,error,*999)
        BASIS%NODE_AT_COLLAPSE=.FALSE.
        
        ALLOCATE(BASIS%NODE_POSITION_INDEX(BASIS%NUMBER_OF_NODES,BASIS%NUMBER_OF_XI_COORDINATES),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate NODE_POSITION_INDEX.",err,error,*999) 
        SELECT CASE(BASIS%NUMBER_OF_XI_COORDINATES)
        CASE(2)
          ALLOCATE(BASIS%NODE_POSITION_INDEX_INV(MAX_NUM_NODES,MAX_NUM_NODES,1,1),STAT=ERR)
        CASE(3)
          ALLOCATE(BASIS%NODE_POSITION_INDEX_INV(MAX_NUM_NODES,MAX_NUM_NODES,MAX_NUM_NODES,1),STAT=ERR)
        CASE(4)
          ALLOCATE(BASIS%NODE_POSITION_INDEX_INV(MAX_NUM_NODES,MAX_NUM_NODES,MAX_NUM_NODES,MAX_NUM_NODES),STAT=ERR)
        CASE DEFAULT
          CALL FlagError("Invalid number of coordinates.",err,error,*999)
        END SELECT
        IF(ERR/=0) CALL FlagError("Could not allocate NODE_POSITION_INDEX_INV.",err,error,*999)
        BASIS%NODE_POSITION_INDEX_INV=0
        
        !Determine the node position index and it's inverse
        SELECT CASE(BASIS%NUMBER_OF_XI)
        CASE(1)
          SELECT CASE(BASIS%INTERPOLATION_ORDER(1))
          CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
            !Node 1
            BASIS%NODE_POSITION_INDEX(1,1)=2
            BASIS%NODE_POSITION_INDEX(1,2)=1
            BASIS%NODE_POSITION_INDEX_INV(2,1,1,1)=1
            !Node 2
            BASIS%NODE_POSITION_INDEX(2,1)=1
            BASIS%NODE_POSITION_INDEX(2,2)=2
            BASIS%NODE_POSITION_INDEX_INV(1,2,1,1)=2
          CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
            !Node 1
            BASIS%NODE_POSITION_INDEX(1,1)=3
            BASIS%NODE_POSITION_INDEX(1,2)=1
            BASIS%NODE_POSITION_INDEX_INV(3,1,1,1)=1
            !Node 2
            BASIS%NODE_POSITION_INDEX(2,1)=2
            BASIS%NODE_POSITION_INDEX(2,2)=2
            BASIS%NODE_POSITION_INDEX_INV(2,2,1,1)=2
            !Node 3
            BASIS%NODE_POSITION_INDEX(3,1)=1
            BASIS%NODE_POSITION_INDEX(3,2)=3
            BASIS%NODE_POSITION_INDEX_INV(1,3,1,1)=3
          CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
            !Node 1
            BASIS%NODE_POSITION_INDEX(1,1)=4
            BASIS%NODE_POSITION_INDEX(1,2)=1
            BASIS%NODE_POSITION_INDEX_INV(4,1,1,1)=1
            !Node 2
            BASIS%NODE_POSITION_INDEX(2,1)=3
            BASIS%NODE_POSITION_INDEX(2,2)=2
            BASIS%NODE_POSITION_INDEX_INV(3,2,1,1)=2
            !Node 3
            BASIS%NODE_POSITION_INDEX(3,1)=2
            BASIS%NODE_POSITION_INDEX(3,2)=3
            BASIS%NODE_POSITION_INDEX_INV(2,3,1,1)=3
            !Node 4
            BASIS%NODE_POSITION_INDEX(4,1)=1
            BASIS%NODE_POSITION_INDEX(4,2)=4
            BASIS%NODE_POSITION_INDEX_INV(1,4,1,1)=4
          CASE DEFAULT
            CALL FlagError("Invalid interpolation order.",err,error,*999)
          END SELECT
        CASE(2)
          SELECT CASE(BASIS%INTERPOLATION_ORDER(1))
          CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
            !Node 1
            BASIS%NODE_POSITION_INDEX(1,1)=2
            BASIS%NODE_POSITION_INDEX(1,2)=1
            BASIS%NODE_POSITION_INDEX(1,3)=1
            BASIS%NODE_POSITION_INDEX_INV(2,1,1,1)=1
            !Node 2
            BASIS%NODE_POSITION_INDEX(2,1)=1
            BASIS%NODE_POSITION_INDEX(2,2)=2
            BASIS%NODE_POSITION_INDEX(2,3)=1
            BASIS%NODE_POSITION_INDEX_INV(1,2,1,1)=2
            !Node 3
            BASIS%NODE_POSITION_INDEX(3,1)=1
            BASIS%NODE_POSITION_INDEX(3,2)=1
            BASIS%NODE_POSITION_INDEX(3,3)=2
            BASIS%NODE_POSITION_INDEX_INV(1,1,2,1)=3
          CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
            !Node 1
            BASIS%NODE_POSITION_INDEX(1,1)=3
            BASIS%NODE_POSITION_INDEX(1,2)=1
            BASIS%NODE_POSITION_INDEX(1,3)=1
            BASIS%NODE_POSITION_INDEX_INV(3,1,1,1)=1
            !Node 2
            BASIS%NODE_POSITION_INDEX(2,1)=1
            BASIS%NODE_POSITION_INDEX(2,2)=3
            BASIS%NODE_POSITION_INDEX(2,3)=1
            BASIS%NODE_POSITION_INDEX_INV(1,3,1,1)=2
            !Node 3
            BASIS%NODE_POSITION_INDEX(3,1)=1
            BASIS%NODE_POSITION_INDEX(3,2)=1
            BASIS%NODE_POSITION_INDEX(3,3)=3
            BASIS%NODE_POSITION_INDEX_INV(1,1,3,1)=3
            !Node 4
            BASIS%NODE_POSITION_INDEX(4,1)=2
            BASIS%NODE_POSITION_INDEX(4,2)=2
            BASIS%NODE_POSITION_INDEX(4,3)=1
            BASIS%NODE_POSITION_INDEX_INV(2,2,1,1)=4
            !Node 5
            BASIS%NODE_POSITION_INDEX(5,1)=1
            BASIS%NODE_POSITION_INDEX(5,2)=2
            BASIS%NODE_POSITION_INDEX(5,3)=2
            BASIS%NODE_POSITION_INDEX_INV(1,2,2,1)=5
            !Node 6
            BASIS%NODE_POSITION_INDEX(6,1)=2
            BASIS%NODE_POSITION_INDEX(6,2)=1
            BASIS%NODE_POSITION_INDEX(6,3)=2
            BASIS%NODE_POSITION_INDEX_INV(2,1,2,1)=6
          CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
            !Node 1
            BASIS%NODE_POSITION_INDEX(1,1)=4
            BASIS%NODE_POSITION_INDEX(1,2)=1
            BASIS%NODE_POSITION_INDEX(1,3)=1
            BASIS%NODE_POSITION_INDEX_INV(4,1,1,1)=1
            !Node 2
            BASIS%NODE_POSITION_INDEX(2,1)=1
            BASIS%NODE_POSITION_INDEX(2,2)=4
            BASIS%NODE_POSITION_INDEX(2,3)=1
            BASIS%NODE_POSITION_INDEX_INV(1,4,1,1)=2
            !Node 3
            BASIS%NODE_POSITION_INDEX(3,1)=1
            BASIS%NODE_POSITION_INDEX(3,2)=1
            BASIS%NODE_POSITION_INDEX(3,3)=4
            BASIS%NODE_POSITION_INDEX_INV(1,1,4,1)=3
            !Node 4
            BASIS%NODE_POSITION_INDEX(4,1)=3
            BASIS%NODE_POSITION_INDEX(4,2)=2
            BASIS%NODE_POSITION_INDEX(4,3)=1
            BASIS%NODE_POSITION_INDEX_INV(3,2,1,1)=4
            !Node 5
            BASIS%NODE_POSITION_INDEX(5,1)=2
            BASIS%NODE_POSITION_INDEX(5,2)=3
            BASIS%NODE_POSITION_INDEX(5,3)=1
            BASIS%NODE_POSITION_INDEX_INV(2,3,1,1)=5
            !Node 6
            BASIS%NODE_POSITION_INDEX(6,1)=1
            BASIS%NODE_POSITION_INDEX(6,2)=3
            BASIS%NODE_POSITION_INDEX(6,3)=2
            BASIS%NODE_POSITION_INDEX_INV(1,3,2,1)=6
            !Node 7
            BASIS%NODE_POSITION_INDEX(7,1)=1
            BASIS%NODE_POSITION_INDEX(7,2)=2
            BASIS%NODE_POSITION_INDEX(7,3)=3
            BASIS%NODE_POSITION_INDEX_INV(1,2,3,1)=7
            !Node 8
            BASIS%NODE_POSITION_INDEX(8,1)=2
            BASIS%NODE_POSITION_INDEX(8,2)=1
            BASIS%NODE_POSITION_INDEX(8,3)=3
            BASIS%NODE_POSITION_INDEX_INV(2,1,3,1)=8
            !Node 9
            BASIS%NODE_POSITION_INDEX(9,1)=3
            BASIS%NODE_POSITION_INDEX(9,2)=1
            BASIS%NODE_POSITION_INDEX(9,3)=2
            BASIS%NODE_POSITION_INDEX_INV(3,1,2,1)=9
            !Node 10
            BASIS%NODE_POSITION_INDEX(10,1)=2
            BASIS%NODE_POSITION_INDEX(10,2)=2
            BASIS%NODE_POSITION_INDEX(10,3)=2
            BASIS%NODE_POSITION_INDEX_INV(2,2,2,1)=10
          CASE DEFAULT
            CALL FlagError("Invalid interpolation order",err,error,*999)
          END SELECT
        CASE(3)
          SELECT CASE(BASIS%INTERPOLATION_ORDER(1))
          CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
            !Node 1
            BASIS%NODE_POSITION_INDEX(1,1)=2
            BASIS%NODE_POSITION_INDEX(1,2)=1
            BASIS%NODE_POSITION_INDEX(1,3)=1
            BASIS%NODE_POSITION_INDEX(1,4)=1
            BASIS%NODE_POSITION_INDEX_INV(2,1,1,1)=1
            !Node 2
            BASIS%NODE_POSITION_INDEX(2,1)=1
            BASIS%NODE_POSITION_INDEX(2,2)=2
            BASIS%NODE_POSITION_INDEX(2,3)=1
            BASIS%NODE_POSITION_INDEX(2,4)=1
            BASIS%NODE_POSITION_INDEX_INV(1,2,1,1)=2
            !Node 3
            BASIS%NODE_POSITION_INDEX(3,1)=1
            BASIS%NODE_POSITION_INDEX(3,2)=1
            BASIS%NODE_POSITION_INDEX(3,3)=2
            BASIS%NODE_POSITION_INDEX(3,4)=1
            BASIS%NODE_POSITION_INDEX_INV(1,1,2,1)=3
            !Node 4
            BASIS%NODE_POSITION_INDEX(4,1)=1
            BASIS%NODE_POSITION_INDEX(4,2)=1
            BASIS%NODE_POSITION_INDEX(4,3)=1
            BASIS%NODE_POSITION_INDEX(4,4)=2
            BASIS%NODE_POSITION_INDEX_INV(1,1,1,2)=4

            ALLOCATE(NODES_IN_FACE(12),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate NODES_IN_FACE.",err,error,*999) 
            NODES_IN_FACE(:)=[2,3,4,1,3,4,1,2,4,1,2,3] !12 Nodes

          CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
            !Node 1
            BASIS%NODE_POSITION_INDEX(1,1)=3
            BASIS%NODE_POSITION_INDEX(1,2)=1
            BASIS%NODE_POSITION_INDEX(1,3)=1
            BASIS%NODE_POSITION_INDEX(1,4)=1
            BASIS%NODE_POSITION_INDEX_INV(3,1,1,1)=1
            !Node 2
            BASIS%NODE_POSITION_INDEX(2,1)=1
            BASIS%NODE_POSITION_INDEX(2,2)=3
            BASIS%NODE_POSITION_INDEX(2,3)=1
            BASIS%NODE_POSITION_INDEX(2,4)=1
            BASIS%NODE_POSITION_INDEX_INV(1,3,1,1)=2
            !Node 3
            BASIS%NODE_POSITION_INDEX(3,1)=1
            BASIS%NODE_POSITION_INDEX(3,2)=1
            BASIS%NODE_POSITION_INDEX(3,3)=3
            BASIS%NODE_POSITION_INDEX(3,4)=1
            BASIS%NODE_POSITION_INDEX_INV(1,1,3,1)=3
            !Node 4
            BASIS%NODE_POSITION_INDEX(4,1)=1
            BASIS%NODE_POSITION_INDEX(4,2)=1
            BASIS%NODE_POSITION_INDEX(4,3)=1
            BASIS%NODE_POSITION_INDEX(4,4)=3
            BASIS%NODE_POSITION_INDEX_INV(1,1,1,3)=4
            !Node 5
            BASIS%NODE_POSITION_INDEX(5,1)=2
            BASIS%NODE_POSITION_INDEX(5,2)=2
            BASIS%NODE_POSITION_INDEX(5,3)=1
            BASIS%NODE_POSITION_INDEX(5,4)=1
            BASIS%NODE_POSITION_INDEX_INV(2,2,1,1)=5
            !Node 6
            BASIS%NODE_POSITION_INDEX(6,1)=2
            BASIS%NODE_POSITION_INDEX(6,2)=1
            BASIS%NODE_POSITION_INDEX(6,3)=2
            BASIS%NODE_POSITION_INDEX(6,4)=1
            BASIS%NODE_POSITION_INDEX_INV(2,1,2,1)=6
            !Node 7
            BASIS%NODE_POSITION_INDEX(7,1)=2
            BASIS%NODE_POSITION_INDEX(7,2)=1
            BASIS%NODE_POSITION_INDEX(7,3)=1
            BASIS%NODE_POSITION_INDEX(7,4)=2
            BASIS%NODE_POSITION_INDEX_INV(2,1,1,2)=7
            !Node 8
            BASIS%NODE_POSITION_INDEX(8,1)=1
            BASIS%NODE_POSITION_INDEX(8,2)=2
            BASIS%NODE_POSITION_INDEX(8,3)=2
            BASIS%NODE_POSITION_INDEX(8,4)=1
            BASIS%NODE_POSITION_INDEX_INV(1,2,2,1)=8
            !Node 9
            BASIS%NODE_POSITION_INDEX(9,1)=1
            BASIS%NODE_POSITION_INDEX(9,2)=1
            BASIS%NODE_POSITION_INDEX(9,3)=2
            BASIS%NODE_POSITION_INDEX(9,4)=2
            BASIS%NODE_POSITION_INDEX_INV(1,1,2,2)=9
            !Node 10
            BASIS%NODE_POSITION_INDEX(10,1)=1
            BASIS%NODE_POSITION_INDEX(10,2)=2
            BASIS%NODE_POSITION_INDEX(10,3)=1
            BASIS%NODE_POSITION_INDEX(10,4)=2
            BASIS%NODE_POSITION_INDEX_INV(1,2,1,2)=10

            ALLOCATE(NODES_IN_FACE(24),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate NODES_IN_FACE.",err,error,*999) 
            NODES_IN_FACE(:)=[2,3,4,8,9,10,1,3,4,6,9,7,1,2,4,5,10,7,1,2,3,5,8,6] !24 Nodes

          CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
            !Node 1
            BASIS%NODE_POSITION_INDEX(1,1)=4
            BASIS%NODE_POSITION_INDEX(1,2)=1
            BASIS%NODE_POSITION_INDEX(1,3)=1
            BASIS%NODE_POSITION_INDEX(1,4)=1
            BASIS%NODE_POSITION_INDEX_INV(4,1,1,1)=1
            !Node 2
            BASIS%NODE_POSITION_INDEX(2,1)=1
            BASIS%NODE_POSITION_INDEX(2,2)=4
            BASIS%NODE_POSITION_INDEX(2,3)=1
            BASIS%NODE_POSITION_INDEX(2,4)=1
            BASIS%NODE_POSITION_INDEX_INV(1,4,1,1)=2
            !Node 3
            BASIS%NODE_POSITION_INDEX(3,1)=1
            BASIS%NODE_POSITION_INDEX(3,2)=1
            BASIS%NODE_POSITION_INDEX(3,3)=4
            BASIS%NODE_POSITION_INDEX(3,4)=1
            BASIS%NODE_POSITION_INDEX_INV(1,1,4,1)=3
            !Node 4
            BASIS%NODE_POSITION_INDEX(4,1)=1
            BASIS%NODE_POSITION_INDEX(4,2)=1
            BASIS%NODE_POSITION_INDEX(4,3)=1
            BASIS%NODE_POSITION_INDEX(4,4)=4
            BASIS%NODE_POSITION_INDEX_INV(1,1,1,4)=4
            !Node 5
            BASIS%NODE_POSITION_INDEX(5,1)=3
            BASIS%NODE_POSITION_INDEX(5,2)=2
            BASIS%NODE_POSITION_INDEX(5,3)=1
            BASIS%NODE_POSITION_INDEX(5,4)=1
            BASIS%NODE_POSITION_INDEX_INV(3,2,1,1)=5
            !Node 6
            BASIS%NODE_POSITION_INDEX(6,1)=2
            BASIS%NODE_POSITION_INDEX(6,2)=3
            BASIS%NODE_POSITION_INDEX(6,3)=1
            BASIS%NODE_POSITION_INDEX(6,4)=1
            BASIS%NODE_POSITION_INDEX_INV(2,3,1,1)=6
            !Node 7
            BASIS%NODE_POSITION_INDEX(7,1)=3
            BASIS%NODE_POSITION_INDEX(7,2)=1
            BASIS%NODE_POSITION_INDEX(7,3)=2
            BASIS%NODE_POSITION_INDEX(7,4)=1
            BASIS%NODE_POSITION_INDEX_INV(3,1,2,1)=7
            !Node 8
            BASIS%NODE_POSITION_INDEX(8,1)=2
            BASIS%NODE_POSITION_INDEX(8,2)=1
            BASIS%NODE_POSITION_INDEX(8,3)=3
            BASIS%NODE_POSITION_INDEX(8,4)=1
            BASIS%NODE_POSITION_INDEX_INV(2,1,3,1)=8
            !Node 9
            BASIS%NODE_POSITION_INDEX(9,1)=3
            BASIS%NODE_POSITION_INDEX(9,2)=1
            BASIS%NODE_POSITION_INDEX(9,3)=1
            BASIS%NODE_POSITION_INDEX(9,4)=2
            BASIS%NODE_POSITION_INDEX_INV(3,1,1,2)=9
            !Node 10
            BASIS%NODE_POSITION_INDEX(10,1)=2
            BASIS%NODE_POSITION_INDEX(10,2)=1
            BASIS%NODE_POSITION_INDEX(10,3)=1
            BASIS%NODE_POSITION_INDEX(10,4)=3
            BASIS%NODE_POSITION_INDEX_INV(2,1,1,3)=10
            !Node 11
            BASIS%NODE_POSITION_INDEX(11,1)=1
            BASIS%NODE_POSITION_INDEX(11,2)=3
            BASIS%NODE_POSITION_INDEX(11,3)=2
            BASIS%NODE_POSITION_INDEX(11,4)=1
            BASIS%NODE_POSITION_INDEX_INV(1,3,2,1)=11
            !Node 12
            BASIS%NODE_POSITION_INDEX(12,1)=1
            BASIS%NODE_POSITION_INDEX(12,2)=2
            BASIS%NODE_POSITION_INDEX(12,3)=3
            BASIS%NODE_POSITION_INDEX(12,4)=1
            BASIS%NODE_POSITION_INDEX_INV(1,2,3,1)=12
            !Node 13
            BASIS%NODE_POSITION_INDEX(13,1)=1
            BASIS%NODE_POSITION_INDEX(13,2)=1
            BASIS%NODE_POSITION_INDEX(13,3)=3
            BASIS%NODE_POSITION_INDEX(13,4)=2
            BASIS%NODE_POSITION_INDEX_INV(1,1,3,2)=13
            !Node 14
            BASIS%NODE_POSITION_INDEX(14,1)=1
            BASIS%NODE_POSITION_INDEX(14,2)=1
            BASIS%NODE_POSITION_INDEX(14,3)=2
            BASIS%NODE_POSITION_INDEX(14,4)=3
            BASIS%NODE_POSITION_INDEX_INV(1,1,2,3)=14
            !Node 15
            BASIS%NODE_POSITION_INDEX(15,1)=1
            BASIS%NODE_POSITION_INDEX(15,2)=3
            BASIS%NODE_POSITION_INDEX(15,3)=1
            BASIS%NODE_POSITION_INDEX(15,4)=2
            BASIS%NODE_POSITION_INDEX_INV(1,3,1,2)=15
            !Node 16
            BASIS%NODE_POSITION_INDEX(16,1)=1
            BASIS%NODE_POSITION_INDEX(16,2)=2
            BASIS%NODE_POSITION_INDEX(16,3)=1
            BASIS%NODE_POSITION_INDEX(16,4)=3
            BASIS%NODE_POSITION_INDEX_INV(1,2,1,3)=16
            !Node 17
            BASIS%NODE_POSITION_INDEX(17,1)=2
            BASIS%NODE_POSITION_INDEX(17,2)=2
            BASIS%NODE_POSITION_INDEX(17,3)=2
            BASIS%NODE_POSITION_INDEX(17,4)=1
            BASIS%NODE_POSITION_INDEX_INV(2,2,2,1)=17
            !Node 18
            BASIS%NODE_POSITION_INDEX(18,1)=2
            BASIS%NODE_POSITION_INDEX(18,2)=2
            BASIS%NODE_POSITION_INDEX(18,3)=1
            BASIS%NODE_POSITION_INDEX(18,4)=2
            BASIS%NODE_POSITION_INDEX_INV(2,2,1,2)=18
            !Node 19
            BASIS%NODE_POSITION_INDEX(19,1)=2
            BASIS%NODE_POSITION_INDEX(19,2)=1
            BASIS%NODE_POSITION_INDEX(19,3)=2
            BASIS%NODE_POSITION_INDEX(19,4)=2
            BASIS%NODE_POSITION_INDEX_INV(2,1,2,2)=19
            !Node 20
            BASIS%NODE_POSITION_INDEX(20,1)=1
            BASIS%NODE_POSITION_INDEX(20,2)=2
            BASIS%NODE_POSITION_INDEX(20,3)=2
            BASIS%NODE_POSITION_INDEX(20,4)=2
            BASIS%NODE_POSITION_INDEX_INV(1,2,2,2)=20

            ALLOCATE(NODES_IN_FACE(40),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate NODES_IN_FACE.",err,error,*999) 
            NODES_IN_FACE(:)=[2,3,4,11,12,13,14,16,15,20,1,3,4,7,8,13,14,10,9,&
                               &19,1,2,4,5,6,15,16,10,9,18,1,2,3,5,6,14,12,8,7,17] !40 nodes

          CASE DEFAULT
            CALL FlagError("Invalid interpolation order.",err,error,*999)
          END SELECT
        CASE DEFAULT
          CALL FlagError("Invalid number of xi directions.",err,error,*999)
        END SELECT
        !Calculate the maximum number of derivatives (1 for simplex bases) and the number of element parameters
        BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES=1
        BASIS%NUMBER_OF_ELEMENT_PARAMETERS=BASIS%NUMBER_OF_NODES
        !Now set up the number of derivatives and derivative order index
        ALLOCATE(BASIS%NUMBER_OF_DERIVATIVES(BASIS%NUMBER_OF_NODES),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate NUMBER_OF_DERIVATIVES.",err,error,*999)
        ALLOCATE(BASIS%DERIVATIVE_ORDER_INDEX(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES,BASIS%NUMBER_OF_NODES,BASIS%NUMBER_OF_XI), &
          & STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate DERIVATIVE_ORDER_INDEX.",err,error,*999)
        ALLOCATE(BASIS%DERIVATIVE_ORDER_INDEX_INV(FIRST_PART_DERIV,FIRST_PART_DERIV,FIRST_PART_DERIV,BASIS%NUMBER_OF_NODES), &
          & STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate DERIVATIVE_ORDER_INDEX_INV.",err,error,*999)
        ALLOCATE(BASIS%PARTIAL_DERIVATIVE_INDEX(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES,BASIS%NUMBER_OF_NODES),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate PARTIAL_DERIVATIVE_INDEX.",err,error,*999)
        ALLOCATE(BASIS%ELEMENT_PARAMETER_INDEX(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES,BASIS%NUMBER_OF_NODES),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate ELEMENT_PARAMETER_INDEX.",err,error,*999)
        ALLOCATE(BASIS%ELEMENT_PARAMETER_INDEX_INV(2,BASIS%NUMBER_OF_ELEMENT_PARAMETERS),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate ELEMENT_PARAMETER_INDEX_INV.",err,error,*999)
        !Set the derivative order index and its inverse, the element parameter index and the partial derivative index.
        ns=0
        BASIS%DERIVATIVE_ORDER_INDEX_INV=0
        DO nn=1,BASIS%NUMBER_OF_NODES
          BASIS%NUMBER_OF_DERIVATIVES(nn)=1
          DO ni=1,BASIS%NUMBER_OF_XI
            BASIS%DERIVATIVE_ORDER_INDEX(1,nn,ni)=1
          ENDDO !ni
          ns=ns+1
          BASIS%ELEMENT_PARAMETER_INDEX(1,nn)=ns
          BASIS%ELEMENT_PARAMETER_INDEX_INV(1,ns)=nn
          BASIS%ELEMENT_PARAMETER_INDEX_INV(2,ns)=1
          BASIS%PARTIAL_DERIVATIVE_INDEX(1,nn)=NO_PART_DERIV
          BASIS%DERIVATIVE_ORDER_INDEX_INV(BASIS%DERIVATIVE_ORDER_INDEX(1,nn,1),1,1,nn)=1
        ENDDO !nn
      
        !Set up the line and face information
        SELECT CASE(BASIS%NUMBER_OF_XI)
        CASE(1)
          BASIS%NUMBER_OF_LOCAL_LINES=1
          ALLOCATE(BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate number of nodes in local line.",err,error,*999)
          ALLOCATE(BASIS%localLineXiDirection(BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate local line xi direction.",err,error,*999)
          BASIS%localLineXiDirection(1)=1
          ALLOCATE(BASIS%NODE_NUMBERS_IN_LOCAL_LINE(BASIS%NUMBER_OF_NODES_XIC(1),BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate node numbers in local line.",err,error,*999)
          ALLOCATE(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(BASIS%NUMBER_OF_NODES_XIC(1),BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate derivative numbers in local line.",err,error,*999)
          BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE=NO_PART_DERIV
          ALLOCATE(basis%ELEMENT_PARAMETERS_IN_LOCAL_LINE(basis%NUMBER_OF_NODES_XIC(1)**2,BASIS%NUMBER_OF_LOCAL_LINES) &
            & ,STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate element parameters in local line.",err,error,*999)
          basis%ELEMENT_PARAMETERS_IN_LOCAL_LINE=1
          !Set the line values
          SELECT CASE(BASIS%INTERPOLATION_ORDER(1))
          CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
            !Line 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(1)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,1)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,1)=2
          CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
            !Line 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(1)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,1)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,1)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,1)=3
          CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
            !Line 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(1)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,1)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,1)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,1)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(4,1)=4
          CASE DEFAULT 
            LOCAL_ERROR="Interpolation order "//TRIM(NumberToVString(BASIS%INTERPOLATION_ORDER(1),"*", &
              & err,error))//" is invalid for a simplex basis type."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(2)
          !Allocate and calculate the lines
          !Simplex hence three local lines
          BASIS%NUMBER_OF_LOCAL_LINES=3
          ALLOCATE(BASIS%localLineXiDirection(BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate local line xi direction",err,error,*999)
          ALLOCATE(BASIS%localLineXiNormals(1,BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate local line xi normals.",err,error,*999)
          ALLOCATE(BASIS%xiNormalsLocalLine(-3:3,1),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate xi normals local line.",err,error,*999)
          basis%xiNormalsLocalLine=0
          ALLOCATE(BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate number of nodes in local line.",err,error,*999)
          ALLOCATE(BASIS%NODE_NUMBERS_IN_LOCAL_LINE(MAXVAL(BASIS%NUMBER_OF_NODES_XIC),BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate node numbers in local line.",err,error,*999)
          ALLOCATE(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(MAXVAL(BASIS%NUMBER_OF_NODES_XIC),BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate derivative numbers in local line.",err,error,*999)
          BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE=NO_PART_DERIV
          ALLOCATE(BASIS%ELEMENT_PARAMETERS_IN_LOCAL_LINE(MAXVAL(BASIS%NUMBER_OF_NODES_XIC)**2,BASIS%NUMBER_OF_LOCAL_LINES) &
            & ,STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate element parameters in local line.",err,error,*999)
          BASIS%ELEMENT_PARAMETERS_IN_LOCAL_LINE=1
          !Set the line values
          SELECT CASE(BASIS%INTERPOLATION_ORDER(1))
          CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
            !Line 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(1)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,1)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,1)=3
            !Line 2
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(2)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,2)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,2)=1
            !Line 3
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(3)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,3)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,3)=2
          CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
            !Line 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(1)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,1)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,1)=5
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,1)=3
            !Line 2
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(2)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,2)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,2)=6
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,2)=2
            !Line 3
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(3)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,3)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,3)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,3)=2
          CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
            !Line 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(1)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,1)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,1)=6
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,1)=7
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(4,1)=3
            !Line 2
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(2)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,2)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,2)=8
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,2)=9
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(4,2)=1
            !Line 3
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(3)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,3)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,3)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,3)=5
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(4,3)=2
          CASE DEFAULT 
            LOCAL_ERROR="Interpolation order "//TRIM(NumberToVString(BASIS%INTERPOLATION_ORDER(1),"*",err,error))// &
              & " is invalid for a simplex basis type."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
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
          BASIS%NUMBER_OF_LOCAL_LINES=6
          BASIS%NUMBER_OF_LOCAL_FACES=4

          ALLOCATE(BASIS%localLineXiDirection(BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate local line xi direction.",err,error,*999)
          ALLOCATE(BASIS%localLineXiNormals(2,BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate local line xi normals.",err,error,*999)
          ALLOCATE(BASIS%xiNormalsLocalLine(-4:4,-4:4),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate xi normals local line.",err,error,*999)
          basis%xiNormalsLocalLine=0

          ALLOCATE(BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate number of nodes in local line.",err,error,*999)
          ALLOCATE(BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate number of nodes in local face.",err,error,*999)

          ALLOCATE(BASIS%NODE_NUMBERS_IN_LOCAL_LINE(MAXVAL(BASIS%NUMBER_OF_NODES_XIC),BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate node numbers in local line.",err,error,*999)
          ALLOCATE(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(MAXVAL(BASIS%NUMBER_OF_NODES_XIC),BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate derivative numbers in local line.",err,error,*999)
          BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE=NO_PART_DERIV
          ALLOCATE(BASIS%ELEMENT_PARAMETERS_IN_LOCAL_LINE(MAXVAL(BASIS%NUMBER_OF_NODES_XIC)**2,BASIS%NUMBER_OF_LOCAL_LINES), &
            & STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate element parameters in local line.",err,error,*999)
          BASIS%ELEMENT_PARAMETERS_IN_LOCAL_LINE=1

          ALLOCATE(BASIS%localFaceXiDirections(3,BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate local face xi directions.",err,error,*999)
          ALLOCATE(BASIS%localFaceXiNormal(BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate local face xi normal.",err,error,*999)
          ALLOCATE(BASIS%xiNormalLocalFace(-4:4),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate xi normal local face.",err,error,*999)
          basis%xiNormalLocalFace=0
          
          SELECT CASE(BASIS%INTERPOLATION_ORDER(1))
          CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
            ALLOCATE(BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate node numbers in local face.",err,error,*999)
            !\todo Number of local face node derivatives currenlty set to 1 (calculation of BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE for simplex elements has not been implemented yet)
            ALLOCATE(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(0:1,3,BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate derivative numbers in local face.",err,error,*999)
          CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
            ALLOCATE(BASIS%NODE_NUMBERS_IN_LOCAL_FACE(6,BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate node numbers in local face.",err,error,*999)
            !\todo Number of local face node derivatives currenlty set to 1 (calculation of BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE for simplex elements has not been implemented yet)
            ALLOCATE(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(0:1,6,BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate derivative numbers in local face.",err,error,*999)
          CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
            ALLOCATE(BASIS%NODE_NUMBERS_IN_LOCAL_FACE(10,BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate node numbers in local face.",err,error,*999)
            !\todo Number of local face node derivatives currenlty set to 1 (calculation of BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE for simplex elements has not been implemented yet)
            ALLOCATE(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(0:1,10,BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate derivative numbers in local face.",err,error,*999)
          CASE DEFAULT
            LOCAL_ERROR="Interpolation order "//TRIM(NumberToVString(BASIS%INTERPOLATION_ORDER(1),"*",err,error))// &
              & " is invalid for a simplex basis type."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
          BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(1,:,:)=NO_PART_DERIV
          BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(0,:,:)=1
          
          !Set the line and face values
          SELECT CASE(BASIS%INTERPOLATION_ORDER(1))
          CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
            !Line 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(1)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,1)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,1)=2
            !Line 2
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(2)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,2)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,2)=3
            !Line 3
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(3)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,3)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,3)=4
            !Line 4
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(4)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,4)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,4)=3
            !Line 5
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(5)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,5)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,5)=4
            !Line 6
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(6)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,6)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,6)=4
            !Face 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(1)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,1)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,1)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,1)=4
            !Face 2
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(2)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,2)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,2)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,2)=3
            !Face 3
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(3)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,3)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,3)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,3)=4
            !Face 4
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(4)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,4)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,4)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,4)=2
          CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
            !Line 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(1)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,1)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,1)=5
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,1)=2
            !Line 2
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(2)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,2)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,2)=6
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,2)=3
            !Line 3
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(3)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,3)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,3)=7
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,3)=4
            !Line 4
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(4)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,4)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,4)=8
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,4)=3
            !Line 5
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(5)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,5)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,5)=10
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,5)=4
            !Line 6
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(6)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,6)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,6)=9
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,6)=4
            !Face 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(1)=6
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,1)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,1)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,1)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(4,1)=8
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(5,1)=9
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(6,1)=10
            !Face 2
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(2)=6
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,2)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,2)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,2)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(4,2)=7
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(5,2)=9
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(6,2)=6
            !Face 3
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(3)=6
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,3)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,3)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,3)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(4,3)=5
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(5,3)=10
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(6,3)=7
            !Face 4
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(4)=6
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,4)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,4)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,4)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(4,4)=6
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(5,4)=8
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(6,4)=5
           CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
            !Line 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(1)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,1)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,1)=5
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,1)=6
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(4,1)=2
            !Line 2
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(2)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,2)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,2)=7
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,2)=8
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(4,2)=3
            !Line 3
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(3)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,3)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,3)=9
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,3)=10
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(4,3)=4
            !Line 4
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(4)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,4)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,4)=11
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,4)=12
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(4,4)=3
            !Line 5
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(5)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,5)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,5)=15
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,5)=16
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(4,5)=4
            !Line 6
            BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(6)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(1,6)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(2,6)=13
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(3,6)=14
            BASIS%NODE_NUMBERS_IN_LOCAL_LINE(4,6)=4
            !Face 1
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(1)=10
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,1)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,1)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,1)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(4,1)=11
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(5,1)=12
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(6,1)=13
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(7,1)=14
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(8,1)=16
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(9,1)=15
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(10,1)=20
            !Face 2
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(2)=10
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,2)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,2)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,2)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(4,2)=9
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(5,2)=10
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(6,2)=14
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(7,2)=13
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(8,2)=8
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(9,2)=7
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(10,2)=19
            !Face 3
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(3)=10
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,3)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,3)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,3)=4
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(4,3)=5
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(5,3)=6
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(6,3)=15
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(7,3)=16
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(8,3)=10
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(9,3)=9
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(10,3)=18
            !Face 4
            BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(4)=10
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(1,4)=1
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(2,4)=3
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(3,4)=2
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(4,4)=7
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(5,4)=8
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(6,4)=12
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(7,4)=11
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(8,4)=6
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(9,4)=5
            BASIS%NODE_NUMBERS_IN_LOCAL_FACE(10,4)=17
          CASE DEFAULT
            LOCAL_ERROR="Interpolation order "//TRIM(NumberToVString(BASIS%INTERPOLATION_ORDER(1),"*",err,error))// &
              & " is invalid for a simplex basis type."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
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
          CALL FlagError("Invalid number of xi directions.",err,error,*999)
        END SELECT
        
        CALL BASIS_QUADRATURE_CREATE(BASIS,err,error,*999)
        
      ELSE
        CALL FlagError("Basis is not a simplex basis.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated.",err,error,*999)
    ENDIF
       
    EXITS("Basis_SimplexBasisCreate")
    RETURN
999 ERRORSEXITS("Basis_SimplexBasisCreate",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_SimplexBasisCreate

  !
  !================================================================================================================================
  !

  !>Creates and initialises a simplex basis family that has already been allocated by BASIS_CREATE_START.
  !> \see BASIS_ROUTINES::BASIS_SIMPLEX_BASIS_CREATE
  SUBROUTINE Basis_SimplexFamilyCreate(basis,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: basis !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,faceXi(2),faceXi2(2),localFaceIdx,localLineIdx,xiIdx,xiIdx2
    LOGICAL :: lineBasisDone,faceBasisDone
    TYPE(BASIS_TYPE), POINTER :: newSubBasis
    TYPE(VARYING_STRING) :: dummyError

    NULLIFY(newSubBasis)

    ENTERS("Basis_SimplexFamilyCreate",err,error,*999)

    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    
    !Create the main (parent) basis
    CALL Basis_SimplexBasisCreate(basis,err,error,*999)
    
    IF(basis%NUMBER_OF_XI>1) THEN
      !Create the line bases as sub-basis types
      ALLOCATE(basis%lineBases(basis%NUMBER_OF_XI),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate basis line bases.",err,error,*999)
      ALLOCATE(basis%localLineBasis(basis%NUMBER_OF_LOCAL_LINES),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate local line basis.",err,error,*999)
      DO xiIdx=1,basis%NUMBER_OF_XI
        lineBasisDone=.FALSE.
        NULLIFY(newSubBasis)
        DO xiIdx2=1,xiIdx-1
          IF(basis%INTERPOLATION_XI(xiIdx2)==basis%INTERPOLATION_XI(xiIdx).AND. &
            basis%QUADRATURE%NUMBER_OF_GAUSS_XI(xiIdx2)==basis%QUADRATURE%NUMBER_OF_GAUSS_XI(xiIdx)) THEN
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
      DO localLineIdx=1,basis%NUMBER_OF_LOCAL_LINES
!!\TODO Until we allow for tensor products of simplex bases each line basis will be the same and equal to the first one.
        basis%localLineBasis(localLineIdx)%ptr=>basis%lineBases(1)%ptr       
      ENDDO !localLineIdx
      IF(basis%NUMBER_OF_XI>2) THEN
        !Create the face bases as sub-basis types
        ALLOCATE(basis%faceBases(basis%NUMBER_OF_XI),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate basis face bases.",err,error,*999)
        ALLOCATE(basis%localFaceBasis(basis%NUMBER_OF_LOCAL_FACES),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate local face basis.",err,error,*999)
        DO xiIdx=1,basis%NUMBER_OF_XI
          !Determine the face xi directions that lie in this xi direction
          faceXi(1)=OTHER_XI_DIRECTIONS3(xiIdx,2,1)
          faceXi(2)=OTHER_XI_DIRECTIONS3(xiIdx,3,1)
          faceBasisDone=.FALSE.
          NULLIFY(newSubBasis)
          DO xiIdx2=1,xiIdx-1
            !Determine the face xi directions that lie in this xi direction
            faceXi2(1)=OTHER_XI_DIRECTIONS3(xiIdx2,2,1)
            faceXi2(2)=OTHER_XI_DIRECTIONS3(xiIdx2,3,1)
            IF(basis%INTERPOLATION_XI(faceXi2(1))==basis%INTERPOLATION_XI(faceXi(1)).AND. &
              & basis%INTERPOLATION_XI(faceXi2(2))==basis%INTERPOLATION_XI(faceXi(2)).AND. &
              & basis%QUADRATURE%NUMBER_OF_GAUSS_XI(faceXi2(1))==basis%QUADRATURE%NUMBER_OF_GAUSS_XI(faceXi(1)).AND. &
              & basis%QUADRATURE%NUMBER_OF_GAUSS_XI(faceXi2(2))==basis%QUADRATURE%NUMBER_OF_GAUSS_XI(faceXi(2))) THEN
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
            ALLOCATE(newSubBasis%localLineBasis(newSubBasis%NUMBER_OF_LOCAL_LINES),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate new sub basis local line basis.",err,error,*999)
             DO localLineIdx=1,newSubBasis%NUMBER_OF_LOCAL_LINES
              newSubBasis%localLineBasis(localLineIdx)%ptr=>newSubBasis%lineBases(1)%ptr
            ENDDO !localFaceIdx
          ENDIF
        ENDDO !xiIdx
        DO localFaceIdx=1,basis%NUMBER_OF_LOCAL_FACES
!!\TODO Until we allow for tensor products of simplex bases each face basis will be the same and equal to the first one.
          basis%localFaceBasis(localFaceIdx)%ptr=>basis%faceBases(1)%ptr
        ENDDO !localFaceIdx
      ENDIF
    ENDIF
    
    EXITS("Basis_SimplexFamilyCreate")
    RETURN
999 IF(ASSOCIATED(newSubBasis)) CALL Basis_FamilyDestroy(newSubBasis%USER_NUMBER,newSubBasis%FAMILY_NUMBER, &
      & dummyErr,dummyError,*998)
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
  FUNCTION BASIS_SIMPLEX_BASIS_EVALUATE(BASIS,NODE_NUMBER,PARTIAL_DERIV_INDEX,XL,err,error)
    
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis function to evaluate.
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBER !<The node number defines the actual basis function to evaluate.
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIV_INDEX !<The partial derivative index in Xi space of the basis to evaluate.
    REAL(DP), INTENT(IN) :: XL(:) !<XL(nic). The area coordinates to evaluate the basis function at.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: BASIS_SIMPLEX_BASIS_EVALUATE !<On return the evaluated basis function
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("BASIS_SIMPLEX_BASIS_EVALUATE",err,error,*999)
    
    BASIS_SIMPLEX_BASIS_EVALUATE=0.0_DP
    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%TYPE==BASIS_SIMPLEX_TYPE) THEN
        SELECT CASE(BASIS%NUMBER_OF_XI)
        CASE(1)
          SELECT CASE(PARTIAL_DERIV_INDEX)
          CASE(NO_PART_DERIV)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,NO_PART_DERIV,XL,err,error)
          CASE(PART_DERIV_S1)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1,XL,err,error)
          CASE(PART_DERIV_S1_S1)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S1,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & 2.0_DP*BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S2,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2_S2,XL,err,error)
          CASE DEFAULT
            LOCAL_ERROR="The specified partial derivative index of "//TRIM(NumberToVString(PARTIAL_DERIV_INDEX,"*",err,error))// &
              & " is invalid for a Simplex line basis."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(2)
          SELECT CASE(PARTIAL_DERIV_INDEX)
          CASE(NO_PART_DERIV)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,NO_PART_DERIV,XL,err,error)
          CASE(PART_DERIV_S1)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S3,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1,XL,err,error)
          CASE(PART_DERIV_S1_S1)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S1,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & 2.0_DP*BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S3,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S3_S3,XL,err,error)
          CASE(PART_DERIV_S2)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S3,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2,XL,err,error)
          CASE(PART_DERIV_S2_S2)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2_S2,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & 2.0_DP*BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2_S3,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S3_S3,XL,err,error)
          CASE(PART_DERIV_S1_S2)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S3_S3,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S3,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2_S3,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S2,XL,err,error)
          CASE DEFAULT
            LOCAL_ERROR="The specified partial derivative index of "//TRIM(NumberToVString(PARTIAL_DERIV_INDEX,"*",err,error))// &
              & " is invalid for a Simplex triangle basis."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(3)
          SELECT CASE(PARTIAL_DERIV_INDEX)
          CASE(NO_PART_DERIV)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,NO_PART_DERIV,XL,err,error)
          CASE(PART_DERIV_S1)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S4,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1,XL,err,error)
          CASE(PART_DERIV_S1_S1)
             BASIS_SIMPLEX_BASIS_EVALUATE= &
              BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S1,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & 2.0_DP*BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S4,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S4_S4,XL,err,error)
          CASE(PART_DERIV_S2)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S4,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2,XL,err,error)
          CASE(PART_DERIV_S2_S2)
             BASIS_SIMPLEX_BASIS_EVALUATE= &
              BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2_S2,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & 2.0_DP*BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2_S4,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S4_S4,XL,err,error)
          CASE(PART_DERIV_S1_S2)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S4_S4,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S4,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2_S4,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S2,XL,err,error)
          CASE(PART_DERIV_S3)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S4,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S3,XL,err,error)
          CASE(PART_DERIV_S3_S3)
             BASIS_SIMPLEX_BASIS_EVALUATE= &
               BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S3_S3,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & 2.0_DP*BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S3_S4,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S4_S4,XL,err,error)
          CASE(PART_DERIV_S1_S3)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S4_S4,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S4,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S3_S4,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S3,XL,err,error)
          CASE(PART_DERIV_S2_S3)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S4_S4,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2_S4,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S3_S4,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2_S3,XL,err,error)
          CASE(PART_DERIV_S1_S2_S3)
            BASIS_SIMPLEX_BASIS_EVALUATE= &
              BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S4_S4_S4,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S4_S4,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2_S4_S4,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S3_S4_S4,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S2_S4,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S3_S4,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE+ &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S2_S3_S4,XL,err,error)
            IF(ERR/=0) GOTO 999
            BASIS_SIMPLEX_BASIS_EVALUATE=BASIS_SIMPLEX_BASIS_EVALUATE- &
              & BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PART_DERIV_S1_S2_S3,XL,err,error)
          CASE DEFAULT
            LOCAL_ERROR="The specified partial derivative index of "//TRIM(NumberToVString(PARTIAL_DERIV_INDEX,"*",err,error))// &
              & " is invalid for a Simplex tetrahedra basis."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="Invalid number of Xi coordinates. The number of xi coordinates for this basis is "// &
            & TRIM(NumberToVString(BASIS%NUMBER_OF_XI,"*",err,error))//" which should be between 1 and 3."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
        IF(ERR/=0) GOTO 999
      ELSE
        CALL FlagError("Basis is not a simplex basis.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated.",err,error,*999)
    ENDIF
    
    EXITS("BASIS_SIMPLEX_BASIS_EVALUATE")
    RETURN
999 ERRORSEXITS("BASIS_SIMPLEX_BASIS_EVALUATE",err,error)
    RETURN 
  END FUNCTION BASIS_SIMPLEX_BASIS_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates partial derivatives of a simplex basis function with respect to area coordinates.
  FUNCTION BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE(BASIS,NODE_NUMBER,PARTIAL_DERIV_INDEX,XL,err,error)
    
    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis function to evaluate.
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBER !<The node number defines the actual basis function to evaluate.
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIV_INDEX !<The partial derivative index in area coordinates of the basis to evaluate.
    REAL(DP), INTENT(IN) :: XL(:) !<XL(nic). The area coordinates to evaluate the basis function at.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE !<On return the evaluated basis function
    !Local variables
    INTEGER(INTG) :: nic
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE",err,error,*999)
    
    BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE=1.0_DP
    IF(ASSOCIATED(BASIS)) THEN      
      DO nic=1,BASIS%NUMBER_OF_XI_COORDINATES       
        SELECT CASE(BASIS%INTERPOLATION_ORDER(nic))
        CASE(BASIS_LINEAR_INTERPOLATION_ORDER)
          BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE=BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE* &
            & SIMPLEX_LINEAR_EVALUATE(BASIS%NODE_POSITION_INDEX(NODE_NUMBER,nic), &
            & PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,nic),XL(nic),err,error)
        CASE(BASIS_QUADRATIC_INTERPOLATION_ORDER)
          BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE=BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE* &
            & SIMPLEX_QUADRATIC_EVALUATE(BASIS%NODE_POSITION_INDEX(NODE_NUMBER,nic), &
            & PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,nic),XL(nic),err,error)
        CASE(BASIS_CUBIC_INTERPOLATION_ORDER)
          BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE=BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE* &
            & SIMPLEX_CUBIC_EVALUATE(BASIS%NODE_POSITION_INDEX(NODE_NUMBER,nic), &
            & PARTIAL_DERIVATIVE_INDEX(PARTIAL_DERIV_INDEX,nic),XL(nic),err,error)
        CASE DEFAULT
          LOCAL_ERROR="Interpolation order value "//TRIM(NumberToVString(BASIS%INTERPOLATION_ORDER(nic),"*",err,error))// &
            & " for xi coordinate direction "//TRIM(NumberToVString(nic,"*",err,error))//" is invalid."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
        IF(ERR/=0) GOTO 999
      ENDDO !nic
    ELSE
      CALL FlagError("Basis is not associated.",err,error,*999)
    ENDIF
    
    EXITS("BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE")
    RETURN
999 ERRORSEXITS("BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE",err,error)
    RETURN 
  END FUNCTION BASIS_SIMPLEX_BASIS_DERIVATIVE_EVALUATE

  !
  !================================================================================================================================
  !

  !>Creates a sub-basis on a parent basis.
  SUBROUTINE Basis_SubBasisCreate(parentBasis,numberOfXi,xiDirections,subBasis,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: parentBasis !<A pointer to the parent basis
    INTEGER(INTG), INTENT(IN) :: numberOfXi !<The number of Xi directions to create
    INTEGER(INTG), INTENT(IN) :: xiDirections(:) !<xiDirections(xiIdx). Gives the ii direction indices of the parent basis which are used to create the sub-basis
    TYPE(BASIS_TYPE), POINTER :: subBasis !<On return, a pointer to the created sub-basis. The pointer must be NULL on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: basisIdx,xiIdx,numberCollapsed,numberEndCollapsed
    TYPE(BASIS_TYPE), POINTER :: newSubBasis  
    TYPE(BASIS_PTR_TYPE), ALLOCATABLE :: newSubBases(:)
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
    newSubBasis%USER_NUMBER=parentBasis%USER_NUMBER
    newSubBasis%GLOBAL_NUMBER=parentBasis%GLOBAL_NUMBER
    newSubBasis%FAMILY_NUMBER=parentBasis%numberOfSubBases+1
    newSubBasis%parentBasis=>parentBasis
    newSubBasis%NUMBER_OF_XI=numberOfXi
    newSubBasis%TYPE=parentBasis%TYPE
    ALLOCATE(newSubBasis%INTERPOLATION_XI(numberOfXi),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate sub-basis interpolation xi.",err,error,*999)
    ALLOCATE(newSubBasis%COLLAPSED_XI(numberOfXi),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate sub-basis collapsed xi.",err,error,*999)        
    numberCollapsed=0
    numberEndCollapsed=0
    DO xiIdx=1,numberOfXi
      newSubBasis%INTERPOLATION_XI(xiIdx)=parentBasis%INTERPOLATION_XI(xiDirections(xiIdx))
      newSubBasis%COLLAPSED_XI(xiIdx)=parentBasis%COLLAPSED_XI(xiDirections(xiIdx))
      IF(newSubBasis%COLLAPSED_XI(xiIdx)==BASIS_XI_COLLAPSED) THEN
        numberCollapsed=numberCollapsed+1
      ELSE IF(newSubBasis%COLLAPSED_XI(xiIdx)==BASIS_COLLAPSED_AT_XI0.OR.&
        & newSubBasis%COLLAPSED_XI(xiIdx)==BASIS_COLLAPSED_AT_XI1) THEN
        numberEndCollapsed=numberEndCollapsed+1
      ENDIF
    ENDDO !xiIdx
    IF(numberCollapsed==0.OR.numberEndCollapsed==0) newSubBasis%COLLAPSED_XI(1:numberOfXi)=BASIS_NOT_COLLAPSED
    NULLIFY(newSubBasis%quadrature%basis)
    CALL BASIS_QUADRATURE_INITIALISE(newSubBasis,err,error,*999)
    newSubBasis%quadrature%type=parentBasis%quadrature%type
    DO xiIdx=1,numberOfXi
      newSubBasis%quadrature%NUMBER_OF_GAUSS_XI(xiIdx)=parentBasis%QUADRATURE%NUMBER_OF_GAUSS_XI(xiDirections(xiIdx))
    ENDDO !xiIdx
    newSubBasis%quadrature%GAUSS_ORDER=parentBasis%quadrature%GAUSS_ORDER
    newSubBasis%BASIS_FINISHED=.TRUE.
    IF(numberOfXi>1) THEN
      ALLOCATE(newSubBasis%lineBases(numberOfXi),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate sub-basis line bases",err,error,*999)
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
  
  !>get the type for a basis is identified by a a pointer. \see OpenCMISS::Iron::cmfe_BasisTypeGet
  SUBROUTINE BASIS_TYPE_GET(BASIS,TYPE,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to get
    INTEGER(INTG), INTENT(OUT) :: TYPE !<On return, the type of the specified basis. \see BASIS_ROUTINES_BasisTypes
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("BASIS_TYPE_GET",err,error,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        TYPE=BASIS%TYPE
      ELSE
        CALL FlagError("Basis has not been finished yet",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated",err,error,*999)
    ENDIF
    
    EXITS("BASIS_TYPE_GET")
    RETURN
999 ERRORSEXITS("BASIS_TYPE_GET",err,error)
    RETURN
  END SUBROUTINE BASIS_TYPE_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the type for a basis is identified by a user number.
  SUBROUTINE BASIS_TYPE_SET_NUMBER(USER_NUMBER,TYPE,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the basis to set.
    INTEGER(INTG), INTENT(IN) :: TYPE !<The type of the basis to set \see BASIS_ROUTINES_BasisTypes
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    
    ENTERS("BASIS_TYPE_SET_NUMBER",err,error,*999)

    CALL Basis_UserNumberFind(USER_NUMBER,BASIS,err,error,*999)
    CALL BASIS_TYPE_SET_PTR(BASIS,TYPE,err,error,*999)
    
    EXITS("BASIS_TYPE_SET_NUMBER")
    RETURN
999 ERRORSEXITS("BASIS_TYPE_SET_NUMBER",err,error)
    RETURN 1
  END SUBROUTINE BASIS_TYPE_SET_NUMBER

  !
  !================================================================================================================================
  !

  !>Sets/changes the type for a basis is identified by a a pointer. \see OpenCMISS::Iron::cmfe_BasisTypeGet
  SUBROUTINE BASIS_TYPE_SET_PTR(BASIS,TYPE,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to set
    INTEGER(INTG), INTENT(IN) :: TYPE !<The type of the basis to be set. \see BASIS_ROUTINES_BasisTypes
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("BASIS_TYPE_SET_PTR",err,error,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        CALL FlagError("Basis has been finished",err,error,*999)
      ELSE
        SELECT CASE(TYPE)
        CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
          BASIS%TYPE=BASIS_LAGRANGE_HERMITE_TP_TYPE
        CASE(BASIS_SIMPLEX_TYPE)
          !Reset the quadrature
          CALL BASIS_QUADRATURE_FINALISE(BASIS,err,error,*999)
          !Change the default parameters for the old basis
          BASIS%TYPE=BASIS_SIMPLEX_TYPE
          BASIS%INTERPOLATION_XI(1:BASIS%NUMBER_OF_XI)=BASIS_LINEAR_SIMPLEX_INTERPOLATION
          NULLIFY(BASIS%QUADRATURE%BASIS)
          CALL BASIS_QUADRATURE_INITIALISE(BASIS,err,error,*999)
        CASE DEFAULT
          LOCAL_ERROR="Basis type "//TRIM(NumberToVString(TYPE,"*",err,error))//" is invalid or not implemented"
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated",err,error,*999)
    ENDIF
    
    EXITS("BASIS_TYPE_SET_PTR")
    RETURN
999 ERRORSEXITS("BASIS_TYPE_SET_PTR",err,error)
    RETURN 1
  END SUBROUTINE BASIS_TYPE_SET_PTR

  !
  !================================================================================================================================
  !
  
  !>Gets the collapsed xi flags for a basis is identified by a a pointer. \see OpenCMISS::Iron::cmfe_BasisCollapsedXiGet
  SUBROUTINE BASIS_COLLAPSED_XI_GET(BASIS,COLLAPSED_XI,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: COLLAPSED_XI(:) !<COLLAPSED_XI(ni). On return, the collapse parameter for each Xi direction. \see BASIS_ROUTINES_XiCollapse
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("BASIS_COLLAPSED_XI_GET",err,error,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        IF(SIZE(COLLAPSED_XI,1)>=SIZE(BASIS%COLLAPSED_XI)) THEN
          COLLAPSED_XI=BASIS%COLLAPSED_XI
        ELSE
          LOCAL_ERROR="The size of COLLAPSED_XI is too small. The supplied size is "// &
            & TRIM(NumberToVString(SIZE(COLLAPSED_XI,1),"*",err,error))//" and it needs to be >= "// &
            & TRIM(NumberToVString(SIZE(BASIS%COLLAPSED_XI,1),"*",err,error))//"."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Basis has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated.",err,error,*999)
    ENDIF
    
    EXITS("BASIS_COLLAPSED_XI_GET")
    RETURN
999 ERRORSEXITS("BASIS_COLLAPSED_XI_GET",err,error)
    RETURN
  END SUBROUTINE BASIS_COLLAPSED_XI_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the collapsed xi flags for a basis is identified by a user number.
  SUBROUTINE BASIS_COLLAPSED_XI_SET_NUMBER(USER_NUMBER,COLLAPSED_XI,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the basis to be set
    INTEGER(INTG), INTENT(IN) :: COLLAPSED_XI(:) !<COLLAPSED_XI(ni). The collapse parameter for each Xi direction. \see BASIS_ROUTINES_XiCollapse
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: BASIS
    
    ENTERS("BASIS_COLLAPSED_XI_SET_NUMBER",err,error,*999)

    CALL Basis_UserNumberFind(USER_NUMBER,BASIS,err,error,*999)
    CALL BASIS_COLLAPSED_XI_SET_PTR(BASIS,COLLAPSED_XI,err,error,*999)
    
    EXITS("BASIS_COLLAPSED_XI_SET_NUMBER")
    RETURN
999 ERRORSEXITS("BASIS_COLLAPSED_XI_SET_NUMBER",err,error)
    RETURN 1
  END SUBROUTINE BASIS_COLLAPSED_XI_SET_NUMBER

  !
  !================================================================================================================================
  !

  !>Sets/changes the collapsed xi flags for a basis is identified by a a pointer. \see OpenCMISS::Iron::cmfe_BasisCollapsedXiSet
  SUBROUTINE BASIS_COLLAPSED_XI_SET_PTR(BASIS,COLLAPSED_XI,err,error,*)

    !Argument variables
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: COLLAPSED_XI(:) !<COLLAPSED_XI(ni). The collapse parameter for each Xi direction. \see BASIS_ROUTINES_XiCollapse
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ni1,ni2,ni3,NUMBER_COLLAPSED,COLLAPSED_XI_DIR(3)
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("BASIS_COLLAPSED_XI_SET_PTR",err,error,*999)

    IF(ASSOCIATED(BASIS)) THEN
      IF(BASIS%BASIS_FINISHED) THEN
        CALL FlagError("Basis has been finished",err,error,*999)
      ELSE
        IF(BASIS%TYPE==BASIS_LAGRANGE_HERMITE_TP_TYPE) THEN
          IF(BASIS%NUMBER_OF_XI>1) THEN
            IF(SIZE(COLLAPSED_XI,1)==BASIS%NUMBER_OF_XI) THEN
              NUMBER_COLLAPSED=0
              DO ni1=1,BASIS%NUMBER_OF_XI
                SELECT CASE(COLLAPSED_XI(ni1))
                CASE(BASIS_XI_COLLAPSED)
                  NUMBER_COLLAPSED=NUMBER_COLLAPSED+1
                  COLLAPSED_XI_DIR(NUMBER_COLLAPSED)=ni1
                CASE(BASIS_COLLAPSED_AT_XI0,BASIS_COLLAPSED_AT_XI1,BASIS_NOT_COLLAPSED)
                  !Do nothing
                CASE DEFAULT
                  LOCAL_ERROR="Collapsed xi value "//TRIM(NumberToVString(COLLAPSED_XI(ni1),"*",err,error))// &
                    & " in xi direction "//TRIM(NumberToVString(ni1,"*",err,error))//" is invalid"
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              ENDDO !ni1
              IF(NUMBER_COLLAPSED>0) THEN
                IF(NUMBER_COLLAPSED<BASIS%NUMBER_OF_XI) THEN
                  IF(BASIS%NUMBER_OF_XI==2) THEN
                    !Two dimensional collapsed basis
                    ni1=COLLAPSED_XI_DIR(1)
                    ni2=OTHER_XI_DIRECTIONS2(ni1)
                    IF(COLLAPSED_XI(ni2)==BASIS_COLLAPSED_AT_XI0) THEN
                      IF(BASIS%INTERPOLATION_XI(ni2)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
                        & BASIS%INTERPOLATION_XI(ni2)=BASIS_QUADRATIC1_HERMITE_INTERPOLATION
                    ELSE IF(COLLAPSED_XI(ni2)==BASIS_COLLAPSED_AT_XI1) THEN
                      IF(BASIS%INTERPOLATION_XI(ni2)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
                        & BASIS%INTERPOLATION_XI(ni2)=BASIS_QUADRATIC2_HERMITE_INTERPOLATION
                    ELSE
                      LOCAL_ERROR="Invalid collapsing of a two dimensional basis. Xi direction "// &
                        & TRIM(NumberToVString(ni1,"*",err,error))//" is collapsed so xi direction "// &
                        & TRIM(NumberToVString(ni2,"*",err,error))//" must be collapsed at an end"
                      CALL FlagError(LOCAL_ERROR,err,error,*999)
                    ENDIF
                  ELSE
                    !Three dimensional collapsed basis
                    IF(NUMBER_COLLAPSED==1) THEN
                      !One collapse - wedge element
                      ni1=COLLAPSED_XI_DIR(1)
                      ni2=OTHER_XI_DIRECTIONS3(ni1,2,1)
                      ni3=OTHER_XI_DIRECTIONS3(ni1,3,1)
                      IF(COLLAPSED_XI(ni2)==BASIS_NOT_COLLAPSED) THEN
                        IF(COLLAPSED_XI(ni3)==BASIS_COLLAPSED_AT_XI0) THEN
                          IF(BASIS%INTERPOLATION_XI(ni3)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
                            & BASIS%INTERPOLATION_XI(ni3)=BASIS_QUADRATIC1_HERMITE_INTERPOLATION
                        ELSE IF(COLLAPSED_XI(ni3)==BASIS_COLLAPSED_AT_XI1) THEN
                          IF(BASIS%INTERPOLATION_XI(ni3)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
                            & BASIS%INTERPOLATION_XI(ni3)=BASIS_QUADRATIC2_HERMITE_INTERPOLATION
                        ELSE
                          LOCAL_ERROR="Invalid collapsing of a three dimensional basis. Xi direction "// &
                            & TRIM(NumberToVString(ni1,"*",err,error))//" is collapsed and xi direction "// &
                            & TRIM(NumberToVString(ni2,"*",err,error))//" is not collapsed so xi direction "// &
                            & TRIM(NumberToVString(ni3,"*",err,error))//" must be collapsed at an end"
                          CALL FlagError(LOCAL_ERROR,err,error,*999)
                        ENDIF
                      ELSE IF(COLLAPSED_XI(ni3)==BASIS_NOT_COLLAPSED) THEN
                        IF(COLLAPSED_XI(ni2)==BASIS_COLLAPSED_AT_XI0) THEN
                          IF(BASIS%INTERPOLATION_XI(ni2)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
                            & BASIS%INTERPOLATION_XI(ni2)=BASIS_QUADRATIC1_HERMITE_INTERPOLATION
                        ELSE IF(COLLAPSED_XI(ni2)==BASIS_COLLAPSED_AT_XI1) THEN
                          IF(BASIS%INTERPOLATION_XI(ni2)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
                            & BASIS%INTERPOLATION_XI(ni2)=BASIS_QUADRATIC2_HERMITE_INTERPOLATION
                        ELSE
                          LOCAL_ERROR="Invalid collapsing of a three dimensional basis. Xi direction "// &
                            & TRIM(NumberToVString(ni1,"*",err,error))//" is collapsed and xi direction "// &
                            & TRIM(NumberToVString(ni3,"*",err,error))//" is not collapsed so xi direction "// &
                            & TRIM(NumberToVString(ni2,"*",err,error))//" must be collapsed at an end"
                          CALL FlagError(LOCAL_ERROR,err,error,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="Invalid collapsing of a three dimensional basis. Xi direction "// &
                          & TRIM(NumberToVString(ni1,"*",err,error))//" is collapsed so one of xi directions "// &
                          & TRIM(NumberToVString(ni2,"*",err,error))//" or "// &
                          & TRIM(NumberToVString(ni3,"*",err,error))//" must be collapsed at an end"
                        CALL FlagError(LOCAL_ERROR,err,error,*999)
                      ENDIF
                    ELSE
                      !Two collapses - pyramid element
                      ni1=COLLAPSED_XI_DIR(1)
                      ni2=COLLAPSED_XI_DIR(2)
                      ni3=OTHER_XI_DIRECTIONS3(ni1,ni2,2)
                      IF(COLLAPSED_XI(ni3)==BASIS_COLLAPSED_AT_XI0) THEN
                        IF(BASIS%INTERPOLATION_XI(ni3)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
                          & BASIS%INTERPOLATION_XI(ni3)=BASIS_QUADRATIC1_HERMITE_INTERPOLATION
                      ELSE IF(COLLAPSED_XI(ni3)==BASIS_COLLAPSED_AT_XI1) THEN
                        IF(BASIS%INTERPOLATION_XI(ni3)==BASIS_CUBIC_HERMITE_INTERPOLATION) &
                          & BASIS%INTERPOLATION_XI(ni3)=BASIS_QUADRATIC2_HERMITE_INTERPOLATION
                      ELSE
                        LOCAL_ERROR="Invalid collapsing of a three dimensional basis. Xi directions "// &
                          & TRIM(NumberToVString(ni1,"*",err,error))//" and "// &
                          & TRIM(NumberToVString(ni2,"*",err,error))//" are collapsed so xi direction "// &
                          & TRIM(NumberToVString(ni3,"*",err,error))//" must be collapsed at an end"
                        CALL FlagError(LOCAL_ERROR,err,error,*999)
                      ENDIF
                    ENDIF
                  ENDIF
                ELSE
                  LOCAL_ERROR="Invalid collapsing of basis. The number of collapsed directions ("// &
                    & TRIM(NumberToVString(NUMBER_COLLAPSED,"*",err,error))// &
                    & ") must be less than the number of xi directions ("// &
                    & TRIM(NumberToVString(BASIS%NUMBER_OF_XI,"*",err,error))//")"
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                ENDIF
              ELSE
                !No collapses in any xi direction - Reset interpolation_xi if necessary
                DO ni1=1,BASIS%NUMBER_OF_XI
                  IF(BASIS%INTERPOLATION_XI(ni1)==BASIS_QUADRATIC1_HERMITE_INTERPOLATION.OR. &
                    & BASIS%INTERPOLATION_XI(ni1)==BASIS_QUADRATIC2_HERMITE_INTERPOLATION) THEN
                    BASIS%INTERPOLATION_XI(ni1)=BASIS_CUBIC_HERMITE_INTERPOLATION                  
                  ENDIF
                ENDDO
              ENDIF
              BASIS%COLLAPSED_XI(1:BASIS%NUMBER_OF_XI)=COLLAPSED_XI(1:BASIS%NUMBER_OF_XI)
            ELSE
              LOCAL_ERROR="The size of the xi collapsed array ("// &
                & TRIM(NumberToVString(SIZE(COLLAPSED_XI,1),"*",err,error))//") does not match the number of xi directions ("// &
                & TRIM(NumberToVString(BASIS%NUMBER_OF_XI,"*",err,error))//") for basis number "// &
                & TRIM(NumberToVString(BASIS%USER_NUMBER,"*",err,error))
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            ENDIF
          ELSE          
            CALL FlagError("Can not collapse a basis with only 1 xi direction",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Can only set collapsed xi directions for a Lagrange Hermite tensor product basis type",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Basis is not associated",err,error,*999)
    ENDIF
    
    EXITS("BASIS_COLLAPSED_XI_SET_PTR")
    RETURN
999 ERRORSEXITS("BASIS_COLLAPSED_XI_SET_PTR",err,error)
    RETURN 1
  END SUBROUTINE BASIS_COLLAPSED_XI_SET_PTR

  !
  !================================================================================================================================
  !

  !>Converts xi coordinates to area coordinates. \see BASIS_ROUTINES::Basis_AreaToXiCoordinates
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
        x(1,1)=1.0_DP/2.0
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
        wC=-9.0_DP/16.0
        alpha1=25.0_DP
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
  FUNCTION HERMITE_CUBIC_EVALUATE(NODE_INDEX,NODE_DERIVATIVE_INDEX,PARTIAL_DERIVATIVE_INDEX,XI,err,error)
  
   !Argument variables
    INTEGER(INTG), INTENT(IN) :: NODE_INDEX !<The local node number of the basis. Must be between 1 and 2.
    INTEGER(INTG), INTENT(IN) :: NODE_DERIVATIVE_INDEX !<The local derivative number of the basis. Must be between 1 and 2.
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIVATIVE_INDEX !<The partial derivative index to evaluate.
    REAL(DP), INTENT(IN) :: XI !<The Xi location to evaluate the basis at.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: HERMITE_CUBIC_EVALUATE !<On exit the evaluated basis function.
    !Local variables
    
    ENTERS("HERMITE_CUBIC_EVALUATE",err,error,*999)

    HERMITE_CUBIC_EVALUATE=0.0_DP
    SELECT CASE(PARTIAL_DERIVATIVE_INDEX)
    CASE(NO_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SELECT CASE(NODE_DERIVATIVE_INDEX)
        CASE(1)
          HERMITE_CUBIC_EVALUATE=(2.0_DP*XI-3.0_DP)*XI*XI+1.0_DP ! 2xi^3-3xi^2+1
        CASE(2)
          HERMITE_CUBIC_EVALUATE=((XI-2.0_DP)*XI+1.0_DP)*XI ! xi^3-2xi^2+xi
        CASE DEFAULT
          CALL FlagError("Invalid node derivative index",err,error,*999)
        END SELECT
      CASE(2)
        SELECT CASE(NODE_DERIVATIVE_INDEX)
        CASE(1)
          HERMITE_CUBIC_EVALUATE=XI*XI*(3.0_DP-2.0_DP*XI) ! -2xi^3+3xi^2
        CASE(2)
          HERMITE_CUBIC_EVALUATE=XI*XI*(XI-1.0_DP) ! xi^3-xi^2
        CASE DEFAULT
          CALL FlagError("Invalid node derivative index",err,error,*999)
        END SELECT
      CASE DEFAULT
        CALL FlagError("Invalid node index",err,error,*999)
      END SELECT
    CASE(FIRST_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SELECT CASE(NODE_DERIVATIVE_INDEX)
        CASE(1)
          HERMITE_CUBIC_EVALUATE=6.0_DP*XI*(XI-1.0_DP) ! 6xi^2-6xi
        CASE(2)
          HERMITE_CUBIC_EVALUATE=(3.0_DP*XI-4.0_DP)*XI+1.0_DP ! 3xi^2-4xi+1
        CASE DEFAULT
          CALL FlagError("Invalid node derivative index",err,error,*999)
        END SELECT
      CASE(2)
        SELECT CASE(NODE_DERIVATIVE_INDEX)
        CASE(1)
          HERMITE_CUBIC_EVALUATE=6.0_DP*XI*(1.0_DP-XI) ! -6xi^2+6xi
        CASE(2)
          HERMITE_CUBIC_EVALUATE=XI*(3.0_DP*XI-2.0_DP) ! 3xi^2-2xi
        CASE DEFAULT
          CALL FlagError("Invalid node derivative index",err,error,*999)
        END SELECT
      CASE DEFAULT
        CALL FlagError("Invalid node index",err,error,*999)
      END SELECT
    CASE(SECOND_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SELECT CASE(NODE_DERIVATIVE_INDEX)
        CASE(1)
          HERMITE_CUBIC_EVALUATE=12.0_DP*XI-6.0_DP ! 12xi-6
        CASE(2)
          HERMITE_CUBIC_EVALUATE=6.0_DP*XI-4.0_DP ! 6xi-4
        CASE DEFAULT
          CALL FlagError("Invalid node derivative index",err,error,*999)
        END SELECT
      CASE(2)
        SELECT CASE(NODE_DERIVATIVE_INDEX)
        CASE(1)
          HERMITE_CUBIC_EVALUATE=6.0_DP-12.0_DP*XI ! -12xi+6
        CASE(2)
          HERMITE_CUBIC_EVALUATE=6.0_DP*XI-2.0_DP ! 6xi-2
        CASE DEFAULT
          CALL FlagError("Invalid node derivative index",err,error,*999)
        END SELECT
      CASE DEFAULT
        CALL FlagError("Invalid node index",err,error,*999)
      END SELECT
    CASE DEFAULT
      CALL FlagError("Invalid partial derivative index",err,error,*999)
    END SELECT

    EXITS("HERMITE_CUBIC_EVALUATE")
    RETURN
999 ERRORSEXITS("HERMITE_CUBIC_EVALUATE",err,error)
    RETURN 
  END FUNCTION HERMITE_CUBIC_EVALUATE
 
  !
  !================================================================================================================================
  !
  
  !#### Generic-Function: HERMITE_QUADRATIC_EVALUATE
  !###  Description:
  !###    Evaluates a 1D quadratic Hermite basis function at position XI,and with the give NODE_INDEX, NODE_DERIVATIVE_INDEX and
  !###    PARTIAL_DERIVATIVE_INDEX. SPECIAL_NODE_INDEX is the node with no derivative term.
  !###  Child-functions: HERMITE_QUADRATIC_EVALUATE_DP

  !
  !================================================================================================================================
  !

  !>Evaluates a 1D quadratic Hermite basis function
  FUNCTION HERMITE_QUADRATIC_EVALUATE(NODE_INDEX,NODE_DERIVATIVE_INDEX,PARTIAL_DERIVATIVE_INDEX,SPECIAL_NODE_INDEX,XI,err,error)
      
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NODE_INDEX !<The local node number of the basis. Must be between 1 and 2.
    INTEGER(INTG), INTENT(IN) :: NODE_DERIVATIVE_INDEX !<The local derivative number of the basis. Must be between 1 and 2.
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIVATIVE_INDEX !<The partial derivative index to evaluate.
    INTEGER(INTG), INTENT(IN) :: SPECIAL_NODE_INDEX !<The local node number with no derivative term.
    REAL(DP), INTENT(IN) :: XI !<The Xi location to evaluate the basis at.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: HERMITE_QUADRATIC_EVALUATE
    !Local variables
    
    ENTERS("HERMITE_QUADRATIC_EVALUATE",err,error,*999)
    
    HERMITE_QUADRATIC_EVALUATE=0.0_DP
    SELECT CASE(SPECIAL_NODE_INDEX)
    CASE(1)
      SELECT CASE(PARTIAL_DERIVATIVE_INDEX)
      CASE(NO_PART_DERIV)
        SELECT CASE(NODE_INDEX)
        CASE(1)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=(XI-2.0_DP)*XI+1.0_DP ! xi^2-2xi+1
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=0.0_DP ! 0
          CASE DEFAULT
            CALL FlagError("Invalid node derivative index",err,error,*999)
          END SELECT
        CASE(2)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=(2.0_DP-XI)*XI ! -xi^2+2xi
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=(XI-1.0_DP)*XI ! xi^2-xi
          CASE DEFAULT
            CALL FlagError("Invalid node derivative index",err,error,*999)
          END SELECT
        CASE DEFAULT
          CALL FlagError("Invalid node index",err,error,*999)
        END SELECT
      CASE(FIRST_PART_DERIV)
        SELECT CASE(NODE_INDEX)
        CASE(1)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=2.0_DP*XI-2.0_DP  ! 2xi-2
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=0.0_DP ! 0
          CASE DEFAULT
            CALL FlagError("Invalid node derivative index",err,error,*999)
          END SELECT
        CASE(2)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=-2.0_DP*XI+2.0_DP ! -2xi+2
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=2.0_DP*XI-1.0_DP ! 2xi-1
          CASE DEFAULT
            CALL FlagError("Invalid node derivative index",err,error,*999)
          END SELECT
        CASE DEFAULT
          CALL FlagError("Invalid node index",err,error,*999)
        END SELECT
      CASE(SECOND_PART_DERIV)
        SELECT CASE(NODE_INDEX)
        CASE(1)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=2.0_DP ! 2
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=0.0_DP ! 0
          CASE DEFAULT
            CALL FlagError("Invalid node derivative index",err,error,*999)
          END SELECT
        CASE(2)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=-2.0_DP ! -2
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=2.0_DP ! 2
          CASE DEFAULT
            CALL FlagError("Invalid node derivative index",err,error,*999)
          END SELECT
        CASE DEFAULT
          CALL FlagError("Invalid node index",err,error,*999)
        END SELECT
      CASE DEFAULT
        CALL FlagError("Invalid partial derivative index",err,error,*999)
      END SELECT
    CASE(2)
      SELECT CASE(PARTIAL_DERIVATIVE_INDEX)
      CASE(NO_PART_DERIV)
        SELECT CASE(NODE_INDEX)
        CASE(1)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=1.0_DP-XI*XI ! -xi^2+1
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=XI*(1.0_DP-XI) ! -xi^2+xi
          CASE DEFAULT
            CALL FlagError("Invalid node derivative index",err,error,*999)
          END SELECT
        CASE(2)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=XI*XI ! xi^2
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=0.0_DP ! 0
          CASE DEFAULT
            CALL FlagError("Invalid node derivative index",err,error,*999)
          END SELECT
        CASE DEFAULT
          CALL FlagError("Invalid node index",err,error,*999)
        END SELECT
      CASE(FIRST_PART_DERIV)
        SELECT CASE(NODE_INDEX)
        CASE(1)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=-2.0_DP*XI ! -2xi
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=1.0_DP-2.0_DP*XI ! -2xi+1
          CASE DEFAULT
            CALL FlagError("Invalid node derivative index",err,error,*999)
          END SELECT
        CASE(2)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=2.0_DP*XI ! 2xi
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=0.0_DP ! 0
          CASE DEFAULT
            CALL FlagError("Invalid node derivative index",err,error,*999)
          END SELECT
        CASE DEFAULT
          CALL FlagError("Invalid node index",err,error,*999)
        END SELECT
      CASE(SECOND_PART_DERIV)
        SELECT CASE(NODE_INDEX)
        CASE(1)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=-2.0_DP ! -2
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=-2.0_DP ! -2
          CASE DEFAULT
            CALL FlagError("Invalid node derivative index",err,error,*999)
          END SELECT
        CASE(2)
          SELECT CASE(NODE_DERIVATIVE_INDEX)
          CASE(1)
            HERMITE_QUADRATIC_EVALUATE=2.0_DP ! 2
          CASE(2)
            HERMITE_QUADRATIC_EVALUATE=0.0_DP ! 0
          CASE DEFAULT
            CALL FlagError("Invalid node derivative index",err,error,*999)
          END SELECT
        CASE DEFAULT
          CALL FlagError("Invalid node index",err,error,*999)
        END SELECT
      CASE DEFAULT
        CALL FlagError("Invalid partial derivative index",err,error,*999)
      END SELECT
    CASE DEFAULT
      CALL FlagError("Invalid special node index",err,error,*999)
    END SELECT

    EXITS("HERMITE_QUADRATIC_EVALUATE")
    RETURN
999 ERRORSEXITS("HERMITE_QUADRATIC_EVALUATE",err,error)
    RETURN 
  END FUNCTION HERMITE_QUADRATIC_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates a 1D cubic Lagrange basis function.
  FUNCTION LAGRANGE_CUBIC_EVALUATE(NODE_INDEX,PARTIAL_DERIVATIVE_INDEX,XI,err,error)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NODE_INDEX !<The local node of the basis to evaluate. Must be between 1 and 4.
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIVATIVE_INDEX !<The partial derivative to evaluate.
    REAL(DP), INTENT(IN) :: XI !<The Xi location to evaluate the basis at.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: LAGRANGE_CUBIC_EVALUATE !<On exit the evaluated basis function.
    !Local variables
    
    ENTERS("LAGRANGE_CUBIC_EVALUATE",err,error,*999)
    
    LAGRANGE_CUBIC_EVALUATE=0.0_DP
    SELECT CASE(PARTIAL_DERIVATIVE_INDEX)
    CASE(NO_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        LAGRANGE_CUBIC_EVALUATE=0.5_DP*(3.0_DP*XI-1.0_DP)*(3.0_DP*XI-2.0_DP)*(1.0_DP-XI) !
      CASE(2)
        LAGRANGE_CUBIC_EVALUATE=4.5_DP*XI*(3.0_DP*XI-2.0_DP)*(XI-1.0_DP) !
      CASE(3)
        LAGRANGE_CUBIC_EVALUATE=4.5_DP*XI*(3.0_DP*XI-1.0_DP)*(1.0_DP-XI) !
      CASE(4)
        LAGRANGE_CUBIC_EVALUATE=0.5_DP*XI*(3.0_DP*XI-1.0_DP)*(3.0_DP*XI-2.0_DP) !
      CASE DEFAULT
        CALL FlagError("Invalid node index",err,error,*999)
      END SELECT
    CASE(FIRST_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        LAGRANGE_CUBIC_EVALUATE=-13.5_DP*XI*XI+18.0_DP*XI-5.5_DP ! -13.5xi^2+18xi-5.5
      CASE(2)
        LAGRANGE_CUBIC_EVALUATE= 40.5_DP*XI*XI-45.0_DP*XI+9.0_DP ! 40.5xi^2-45xi+9
      CASE(3)
        LAGRANGE_CUBIC_EVALUATE=-40.5_DP*XI*XI+36.0_DP*XI-4.5_DP ! -40.5xi^2+36xi-4.5
      CASE(4)
        LAGRANGE_CUBIC_EVALUATE= 13.5_DP*XI*XI- 9.0_DP*XI+1.0_DP ! 13.5xi^2-9xi+1
      CASE DEFAULT
        CALL FlagError("Invalid node index",err,error,*999)
      END SELECT
    CASE(SECOND_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        LAGRANGE_CUBIC_EVALUATE=9.0_DP*(2.0_DP-3.0_DP*XI) ! 18-27xi
      CASE(2)
        LAGRANGE_CUBIC_EVALUATE=9.0_DP*(9.0_DP*XI-5.0_DP) ! 81xi-45
      CASE(3)
        LAGRANGE_CUBIC_EVALUATE=9.0_DP*(4.0_DP-9.0_DP*XI) ! 36-81xi
      CASE(4)
        LAGRANGE_CUBIC_EVALUATE=9.0_DP*(3.0_DP*XI-1.0_DP) ! 27xi-9
      CASE DEFAULT
        CALL FlagError("Invalid node index",err,error,*999)
      END SELECT
    CASE DEFAULT
      CALL FlagError("Invalid partial derivative index",err,error,*999)
    END SELECT

    EXITS("LAGRANGE_CUBIC_EVALUATE")
    RETURN
999 ERRORSEXITS("LAGRANGE_CUBIC_EVALUATE",err,error)
    RETURN 
  END FUNCTION LAGRANGE_CUBIC_EVALUATE

  !
  !================================================================================================================================
  !
  
  !> Evaluates a 1D linear Lagrange basis function.
  FUNCTION LAGRANGE_LINEAR_EVALUATE(NODE_INDEX,PARTIAL_DERIVATIVE_INDEX,XI,err,error)
      
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NODE_INDEX !<The local node of the basis to evaluate. Must be between 1 and 2.
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIVATIVE_INDEX !<The partial derivative to evaluate.
    REAL(DP), INTENT(IN) :: XI !<The Xi location to evaluate the basis at.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: LAGRANGE_LINEAR_EVALUATE !<On exit the evaluated basis function.
    !Local variables
    
    ENTERS("LAGRANGE_LINEAR_EVALUATE",err,error,*999)

    LAGRANGE_LINEAR_EVALUATE=0.0_DP
    SELECT CASE(PARTIAL_DERIVATIVE_INDEX)
    CASE(NO_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        LAGRANGE_LINEAR_EVALUATE=1.0_DP-XI ! 1-xi
      CASE(2)
        LAGRANGE_LINEAR_EVALUATE=XI !xi
      CASE DEFAULT
        CALL FlagError("Invalid node index",err,error,*999)
      END SELECT
    CASE(FIRST_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        LAGRANGE_LINEAR_EVALUATE=-1.0_DP ! -1
      CASE(2)
        LAGRANGE_LINEAR_EVALUATE=1.0_DP ! 1
      CASE DEFAULT
        CALL FlagError("Invalid node index",err,error,*999)
      END SELECT
    CASE(SECOND_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        LAGRANGE_LINEAR_EVALUATE=0.0_DP ! 0
      CASE(2)
        LAGRANGE_LINEAR_EVALUATE=0.0_DP ! 0
      CASE DEFAULT
        CALL FlagError("Invalid node index",err,error,*999)
      END SELECT
    CASE DEFAULT
      CALL FlagError("Invalid partial derivative index",err,error,*999)
    END SELECT
    
    EXITS("LAGRANGE_LINEAR_EVALUATE")
    RETURN
999 ERRORSEXITS("LAGRANGE_LINEAR_EVALUATE",err,error)
    RETURN 
  END FUNCTION LAGRANGE_LINEAR_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates a 1D quadratic Lagrange basis function.
  FUNCTION LAGRANGE_QUADRATIC_EVALUATE(NODE_INDEX,PARTIAL_DERIVATIVE_INDEX,XI,err,error)
     
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NODE_INDEX !<The local node of the basis to evaluate. Must be between 1 and 3.
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIVATIVE_INDEX !<The partial derivative to evaluate.
    REAL(DP), INTENT(IN) :: XI !<The Xi location to evaluate the basis at.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: LAGRANGE_QUADRATIC_EVALUATE !<On exit the evaluated basis function.
    !Local variables
    
    ENTERS("LAGRANGE_QUADRATIC_EVALUATE",err,error,*999)

    LAGRANGE_QUADRATIC_EVALUATE=0.0_DP
    SELECT CASE(PARTIAL_DERIVATIVE_INDEX)
    CASE(NO_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        LAGRANGE_QUADRATIC_EVALUATE=1.0_DP-3.0_DP*XI+2.0_DP*XI*XI ! 1-3xi+2xi^2
      CASE(2)
        LAGRANGE_QUADRATIC_EVALUATE=4.0_DP*XI*(1.0_DP-XI) ! 4xi-4xi^2
      CASE(3)
        LAGRANGE_QUADRATIC_EVALUATE=XI*(XI+XI-1.0_DP) ! 2xi^2-xi
      CASE DEFAULT
        CALL FlagError("Invalid node index",err,error,*999)
      END SELECT
    CASE(FIRST_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        LAGRANGE_QUADRATIC_EVALUATE=4.0_DP*XI-3.0_DP ! 4xi-3
      CASE(2)
        LAGRANGE_QUADRATIC_EVALUATE=4.0_DP-8.0_DP*XI ! 4-8xi
      CASE(3)
        LAGRANGE_QUADRATIC_EVALUATE=4.0_DP*XI-1.0_DP ! 4xi-1
      CASE DEFAULT
        CALL FlagError("Invalid node index",err,error,*999)
      END SELECT
    CASE(SECOND_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        LAGRANGE_QUADRATIC_EVALUATE=4.0_DP ! 4
      CASE(2)
        LAGRANGE_QUADRATIC_EVALUATE=-8.0_DP ! -8
      CASE(3)
        LAGRANGE_QUADRATIC_EVALUATE=4.0_DP ! 4
      CASE DEFAULT
        CALL FlagError("Invalid node index",err,error,*999)
      END SELECT
    CASE DEFAULT
      CALL FlagError("Invalid partial derivative index",err,error,*999)
    END SELECT

    EXITS("LAGRANGE_QUADRATIC_EVALUATE")
    RETURN
999 ERRORSEXITS("LAGRANGE_QUADRATIC_EVALUATE",err,error)
    RETURN 
  END FUNCTION LAGRANGE_QUADRATIC_EVALUATE

  !
  !================================================================================================================================
  !
  
  !>Evaluates a cubic simpelx basis function at a specificed area position and node index and with a given partial derivative index
  !>with respect to area coordinates for double precision arguments
  FUNCTION SIMPLEX_CUBIC_EVALUATE_DP(NODE_INDEX,PARTIAL_DERIVATIVE_INDEX,XL,err,error)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NODE_INDEX !<The node index to evaluate
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIVATIVE_INDEX !<The partial derivative index wrt area coordinates to evaluate
    REAL(DP), INTENT(IN) :: XL !<The area coordinate to evaluate at.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: SIMPLEX_CUBIC_EVALUATE_DP
    !Local variables
    
    ENTERS("SIMPLEX_CUBIC_EVALUATE_DP",err,error,*999)
    
    SIMPLEX_CUBIC_EVALUATE_DP=0.0_DP
        
    SELECT CASE(PARTIAL_DERIVATIVE_INDEX)
    CASE(NO_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_CUBIC_EVALUATE_DP=1.0 !1
      CASE(2)
        SIMPLEX_CUBIC_EVALUATE_DP=3.0_DP*XL !3L
      CASE(3)
        SIMPLEX_CUBIC_EVALUATE_DP=3.0_DP/2.0_DP*XL*(3.0_DP*XL-1.0_DP) !3/2.L(3L-1)
      CASE(4)
        SIMPLEX_CUBIC_EVALUATE_DP=0.5_DP*XL*(3.0_DP*XL-1.0_DP)*(3.0_DP*XL-2.0_DP) !1/2.L(3L-1)(3L-2)
      CASE DEFAULT
        CALL FlagError("Invalid node index.",err,error,*999)
      END SELECT
    CASE(FIRST_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_CUBIC_EVALUATE_DP=0.0_DP !0
      CASE(2)
        SIMPLEX_CUBIC_EVALUATE_DP=3.0_DP !3
      CASE(3)
        SIMPLEX_CUBIC_EVALUATE_DP=3.0_DP/2.0_DP*(6.0_DP*XL-1) !3/2.(6L-1)
      CASE(4)
        SIMPLEX_CUBIC_EVALUATE_DP=13.5_DP*XL*XL-9.0_DP*XL+1.0_DP !27/2.L^2-9L+1
      CASE DEFAULT
        CALL FlagError("Invalid node index.",err,error,*999)
      END SELECT
    CASE(SECOND_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_CUBIC_EVALUATE_DP=0.0_DP !0
      CASE(2)
        SIMPLEX_CUBIC_EVALUATE_DP=0.0_DP !0
      CASE(3)
        SIMPLEX_CUBIC_EVALUATE_DP=9.0_DP !9
      CASE(4)
        SIMPLEX_CUBIC_EVALUATE_DP=2.0_DP*XL-9.0_DP
      CASE DEFAULT
        CALL FlagError("Invalid node index.",err,error,*999)
      END SELECT
    CASE(THIRD_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_CUBIC_EVALUATE_DP=0.0_DP !0
      CASE(2)
        SIMPLEX_CUBIC_EVALUATE_DP=0.0_DP !0
      CASE(3)
        SIMPLEX_CUBIC_EVALUATE_DP=0.0_DP !0
      CASE(4)
        SIMPLEX_CUBIC_EVALUATE_DP=2.0_DP !2
      CASE DEFAULT
        CALL FlagError("Invalid node index.",err,error,*999)
      END SELECT
    CASE DEFAULT
      CALL FlagError("Invalid partial derivative index.",err,error,*999)
    END SELECT

    EXITS("SIMPLEX_CUBIC_EVALUATE_DP")
    RETURN
999 ERRORSEXITS("SIMPLEX_CUBIC_EVALUATE_DP",err,error)
    RETURN 
  END FUNCTION SIMPLEX_CUBIC_EVALUATE_DP

  !
  !================================================================================================================================
  !
  
  !>Evaluates a linear simpelx basis function at a specificed area position and node index and with a given partial derivative index
  !>with respect to area coordinates for double precision arguments
  FUNCTION SIMPLEX_LINEAR_EVALUATE_DP(NODE_INDEX,PARTIAL_DERIVATIVE_INDEX,XL,err,error)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NODE_INDEX !<The node index to evaluate
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIVATIVE_INDEX !<The partial derivative index wrt area coordinates to evaluate
    REAL(DP), INTENT(IN) :: XL !<The area coordinate to evaluate at.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: SIMPLEX_LINEAR_EVALUATE_DP
    !Local variables
    
    ENTERS("SIMPLEX_LINEAR_EVALUATE_DP",err,error,*999)

    SIMPLEX_LINEAR_EVALUATE_DP=0.0_DP
    SELECT CASE(PARTIAL_DERIVATIVE_INDEX)
    CASE(NO_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_LINEAR_EVALUATE_DP=1.0 !1
      CASE(2)
        SIMPLEX_LINEAR_EVALUATE_DP=XL  !L
      CASE DEFAULT
        CALL FlagError("Invalid node index",err,error,*999)
      END SELECT
    CASE(FIRST_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_LINEAR_EVALUATE_DP=0.0_DP  !0
      CASE(2)
        SIMPLEX_LINEAR_EVALUATE_DP=1.0_DP !1
      CASE DEFAULT
        CALL FlagError("Invalid node index",err,error,*999)
      END SELECT
    CASE(SECOND_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_LINEAR_EVALUATE_DP=0.0_DP !0
      CASE(2)
        SIMPLEX_LINEAR_EVALUATE_DP=0.0_DP !0
      CASE DEFAULT
        CALL FlagError("Invalid node index",err,error,*999)
      END SELECT
    CASE(THIRD_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_LINEAR_EVALUATE_DP=0.0_DP !0
      CASE(2)
        SIMPLEX_LINEAR_EVALUATE_DP=0.0_DP !0
      CASE DEFAULT
        CALL FlagError("Invalid node index",err,error,*999)
      END SELECT
    CASE DEFAULT
      CALL FlagError("Invalid partial derivative index",err,error,*999)
    END SELECT
    
    EXITS("SIMPLEX_LINEAR_EVALUATE_DP")
    RETURN
999 ERRORSEXITS("SIMPLEX_LINEAR_EVALUATE_DP",err,error)
    RETURN 
  END FUNCTION SIMPLEX_LINEAR_EVALUATE_DP

  !
  !================================================================================================================================
  !
  
  !>Evaluates a quadratic simpelx basis function at a specificed area position and node index and with a given partial derivative index
  !>with respect to area coordinates for double precision arguments.
   FUNCTION SIMPLEX_QUADRATIC_EVALUATE_DP(NODE_INDEX,PARTIAL_DERIVATIVE_INDEX,XL,err,error)
      
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NODE_INDEX !<The node index to evaluate
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIVATIVE_INDEX !<The partial derivative index wrt area coordinates to evaluate
    REAL(DP), INTENT(IN) :: XL !<The area coordinate to evaluate at.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: SIMPLEX_QUADRATIC_EVALUATE_DP
    !Local variables
    
    ENTERS("SIMPLEX_QUADRATIC_EVALUATE_DP",err,error,*999)

    SIMPLEX_QUADRATIC_EVALUATE_DP=0.0_DP
    SELECT CASE(PARTIAL_DERIVATIVE_INDEX)
    CASE(NO_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_QUADRATIC_EVALUATE_DP=1.0_DP !1
      CASE(2)
        SIMPLEX_QUADRATIC_EVALUATE_DP=2.0_DP*XL !2L
      CASE(3)
        SIMPLEX_QUADRATIC_EVALUATE_DP=XL*(2.0_DP*XL-1.0_DP) !L(2L-1)
      CASE DEFAULT
        CALL FlagError("Invalid node index.",err,error,*999)
      END SELECT
    CASE(FIRST_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_QUADRATIC_EVALUATE_DP=0.0_DP !0
      CASE(2)
        SIMPLEX_QUADRATIC_EVALUATE_DP=2.0_DP !4
      CASE(3)
        SIMPLEX_QUADRATIC_EVALUATE_DP=4.0_DP*XL-1.0_DP !4L-1
      CASE DEFAULT
        CALL FlagError("Invalid node index",err,error,*999)
      END SELECT
    CASE(SECOND_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_QUADRATIC_EVALUATE_DP=0.0_DP !0
      CASE(2)
        SIMPLEX_QUADRATIC_EVALUATE_DP=0.0_DP !0
      CASE(3)
        SIMPLEX_QUADRATIC_EVALUATE_DP=4.0_DP !4
      CASE DEFAULT
        CALL FlagError("Invalid node index.",err,error,*999)
      END SELECT
    CASE(THIRD_PART_DERIV)
      SELECT CASE(NODE_INDEX)
      CASE(1)
        SIMPLEX_QUADRATIC_EVALUATE_DP=0.0_DP !0
      CASE(2)
        SIMPLEX_QUADRATIC_EVALUATE_DP=0.0_DP !0
      CASE(3)
        SIMPLEX_QUADRATIC_EVALUATE_DP=0.0_DP !0
      CASE DEFAULT
        CALL FlagError("Invalid node index.",err,error,*999)
      END SELECT
    CASE DEFAULT
      CALL FlagError("Invalid partial derivative index.",err,error,*999)
    END SELECT

    EXITS("SIMPLEX_QUADRATIC_EVALUATE_DP")
    RETURN
999 ERRORSEXITS("SIMPLEX_QUADRATIC_EVALUATE_DP",err,error)
    RETURN 
  END FUNCTION SIMPLEX_QUADRATIC_EVALUATE_DP

  !
  !================================================================================================================================
  !

END MODULE BasisRoutines

