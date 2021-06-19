!> \file
!> \author Chris Bradley
!> \brief This module contains all basis access method routines.
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

!> \addtogroup OpenCMISS_Basis OpenCMISS::Iron::Basis
!> This module contains all basis access method routines.
MODULE BasisAccessRoutines
  
  USE BaseRoutines
  USE Kinds
  USE ISO_VARYING_STRING
  USE Strings
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters


  !> \addtogroup BasisRoutines_BasisTypes BasisRoutines::BasisTypes
  !> \brief Basis definition type parameters
  !> \todo Combine simplex and serendipity elements???
  !> \see BasisRoutines,OpenCMISS_BasisTypes
  !>@{ 
  INTEGER(INTG), PARAMETER :: BASIS_LAGRANGE_HERMITE_TP_TYPE=1 !<Lagrange-Hermite tensor product basis type \see BasisRoutines_BasisTypes,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_SIMPLEX_TYPE=2 !<Simplex basis type \see BasisRoutines_BasisTypes,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_SERENDIPITY_TYPE=3 !<Serendipity basis type \see BasisRoutines_BasisTypes,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_AUXILLIARY_TYPE=4 !<Auxillary basis type \see BasisRoutines_BasisTypes,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_B_SPLINE_TP_TYPE=5 !<B-spline basis type \see BasisRoutines_BasisTypes,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE=6 !<Fourier-Lagrange tensor product basis type \see BasisRoutines_BasisTypes,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_EXTENDED_LAGRANGE_TP_TYPE=7 !< Extendend Lagrange tensor product basis type \see BasisRoutines_BasisTypes,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_RADIAL_TYPE=7 !< Radial basis typee \see BasisRoutines_BasisTypes,BasisRoutines
  !>@}

  !> \addtogroup BasisRoutines_InterpolationSpecifications BasisRoutines::InterpolationSpecifications
  !> \brief Interpolation specification parameters
  !> \see BasisRoutines,OpenCMISS_InterpolationSpecifications
  !>@{ 
  INTEGER(INTG), PARAMETER :: BASIS_LINEAR_LAGRANGE_INTERPOLATION=1 !<Linear Lagrange interpolation specification \see BasisRoutines_InterpolationSpecifications,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_QUADRATIC_LAGRANGE_INTERPOLATION=2 !<Quadratic Lagrange interpolation specification \see BasisRoutines_InterpolationSpecifications,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_CUBIC_LAGRANGE_INTERPOLATION=3 !<Cubic Lagrange interpolation specification \see BasisRoutines_InterpolationSpecifications,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_CUBIC_HERMITE_INTERPOLATION=4 !<Cubic Hermite interpolation specification \see BasisRoutines_InterpolationSpecifications,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_QUADRATIC1_HERMITE_INTERPOLATION=5 !<Quadratic Hermite (no derivative at xi=0) interpolation specification \see BasisRoutines_InterpolationSpecifications,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_QUADRATIC2_HERMITE_INTERPOLATION=6 !<Quadratic Hermite (no derivative at xi=1) interpolation specification \see BasisRoutines_InterpolationSpecifications,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_LINEAR_SIMPLEX_INTERPOLATION=7 !<Linear Simplex interpolation specification \see BasisRoutines_InterpolationSpecifications,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_QUADRATIC_SIMPLEX_INTERPOLATION=8 !<Quadratic Simplex interpolation specification \see BasisRoutines_InterpolationSpecifications,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_CUBIC_SIMPLEX_INTERPOLATION=9 !<Cubic Simplex interpolation specification \see BasisRoutines_InterpolationSpecifications,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_GAUSSIAN_RADIAL_INTERPOLATION=10 !<Gaussian Radial interpolation specification \see BasisRoutines_InterpolationSpecifications,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_MULTIQUARTIC_RADIAL_INTERPOLATION=11 !<Multiquartic Radial interpolation specification \see BasisRoutines_InterpolationSpecifications,BasisRoutines
  !>@}

  !> \addtogroup BasisRoutines_InterpolationTypes BasisRoutines::InterpolationTypes
  !> \brief Interpolation type parameters for a Xi direction
  !> \see BasisRoutines
  !>@{ 
  INTEGER(INTG), PARAMETER :: BASIS_LAGRANGE_INTERPOLATION=1 !<Lagrange interpolation \see BasisRoutines_InterpolationTypes,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_HERMITE_INTERPOLATION=2 !<Hermite interpolation \see BasisRoutines_InterpolationTypes,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_SIMPLEX_INTERPOLATION=3 !<Simplex interpolation \see BasisRoutines_InterpolationTypes,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_SERENDIPITY_INTERPOLATION=4 !<Serendipity interpolation \see BasisRoutines_InterpolationTypes,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_TRANSITION_INTERPOLATION=5 !<Transition interpolation \see BasisRoutines_InterpolationTypes,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_SINGULAR_INTERPOLATION=6 !<Singular interpolation \see BasisRoutines_InterpolationTypes,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_FOURIER_INTERPOLATION=7 !<Fourier interpolation \see BasisRoutines_InterpolationTypes,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_RADIAL_INTERPOLATION=8 !<Radial interpolation \see BasisRoutines_InterpolationTypes,BasisRoutines
  !>@}
  
  !> \addtogroup BasisRoutines_InterpolationOrder BasisRoutines::InterpolationOrder
  !> \brief Interpolation order for a Xi direction
  !> \see BasisRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: BASIS_LINEAR_INTERPOLATION_ORDER=1 !<Linear interpolation order \see BasisRoutines_InterpolationOrder,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_QUADRATIC_INTERPOLATION_ORDER=2 !<Quadratic interpolation order \see BasisRoutines_InterpolationOrder,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_CUBIC_INTERPOLATION_ORDER=3 !<Cubic interpolation order \see BasisRoutines_InterpolationOrder,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_QUADRATIC1_INTERPOLATION_ORDER=4 !<Quadratic (no derivative at xi=0) interpolation order \see BasisRoutines_InterpolationOrder,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_QUADRATIC2_INTERPOLATION_ORDER=5 !<Quadratic (no derivative at xi=1) interpolation order \see BasisRoutines_InterpolationOrder,BasisRoutines
  !>@}
  
  !> \addtogroup BasisRoutines_QuadratureTypes BasisRoutines::QuadratureTypes
  !> \brief Quadrature type parameters
  !> \see BasisRoutines,OpenCMISS_QuadratureTypes
  !>@{
  INTEGER(INTG), PARAMETER :: BASIS_GAUSS_LEGENDRE_QUADRATURE=1 !<Gauss-Legendre quadrature  \see BasisRoutines_QuadratureTypes,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_GAUSS_LAGUERRE_QUADRATURE=2 !<Gauss-Laguerre quadrature  \see BasisRoutines_QuadratureTypes,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_GUASS_HERMITE_QUADRATURE=3 !<Gauss-Hermite quadrature  \see BasisRoutines_QuadratureTypes,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_ADAPTIVE_GAUSS_LEGENDRE_QUADRATURE=4 !<Adaptive Gauss-Legendre quadrature  \see BasisRoutines_QuadratureTypes,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_GAUSS_SIMPLEX_QUADRATURE=5 !<Gauss-Legendre for Simplex elements quadrature  \see BasisRoutines_QuadratureTypes,BasisRoutines
  !>@}

  !> \addtogroup BasisRoutines_XiCollapse BasisRoutines::XiCollapse
  !> \brief Xi collapse parameters
  !> \see BasisRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: BASIS_XI_COLLAPSED=1 !<The Xi direction is collapsed \see BasisRoutines_XiCollapse,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_COLLAPSED_AT_XI0=2 !<The Xi direction at the xi=0 end of this Xi direction is collapsed \see BasisRoutines_XiCollapse,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_COLLAPSED_AT_XI1=3 !<The Xi direction at the xi=1 end of this Xi direction is collapsed \see BasisRoutines_XiCollapse,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_NOT_COLLAPSED=4 !<The Xi direction is not collapsed \see BasisRoutines_XiCollapse,BasisRoutines
  !>@}
  
  !> \addtogroup BasisRoutines_BoundaryXiType BasisRoutines::BoundaryXiType
  !> \brief The boundary xi types
  !> \see BasisRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: BASIS_NO_BOUNDARY_XI=0 !<The xi is not on a boundary \see BasisRoutines_XiBoundaryType,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_LINE_BOUNDARY_XI=1 !<The xi is on a boundary line \see BasisRoutines_XiBoundaryType,BasisRoutines
  INTEGER(INTG), PARAMETER :: BASIS_FACE_BOUNDARY_XI=2 !<The xi is on a boundary face \see BasisRoutines_XiBoundaryType,BasisRoutines
  !>@}
  
  !> \addtogroup Basis_QuadratureSchemes Basis::QuadratureSchemes  
  !> \brief Quadrature scheme parameters. NOTE: Quadratures schemes have not been implemented yet. For now you should just use the BASIS_DEFAULT_QUADRATURE_SCHEME.
  !> \see BasisRoutines,BasisAccessRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES=4 !<The number of currently defined quadrature schemes \see Basis_QuadratureScheme,BasisRoutines,BasisAccessRoutines
  INTEGER(INTG), PARAMETER :: BASIS_DEFAULT_QUADRATURE_SCHEME=1 !<Identifier for the default quadrature scheme \see Basis_QuadratureScheme,BasisRoutines,BasisAccessRoutines
  INTEGER(INTG), PARAMETER :: BASIS_LOW_QUADRATURE_SCHEME=2 !<Identifier for a low order quadrature scheme \see Basis_QuadratureScheme,BasisRoutines,BasisAccessRoutines
  INTEGER(INTG), PARAMETER :: BASIS_MID_QUADRATURE_SCHEME=3 !<Identifier for a mid order quadrature scheme \see Basis_QuadratureScheme,BasisRoutines,BasisAccessRoutines
  INTEGER(INTG), PARAMETER :: BASIS_HIGH_QUADRATURE_SCHEME=4 !<Identifier for a high order quadrature scheme \see Basis_QuadratureScheme,BasisRoutines,BasisAccessRoutines
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

  INTERFACE Basis_LineXiNormalsGet
    MODULE PROCEDURE Basis_LineXiNormalsGet0
    MODULE PROCEDURE Basis_LineXiNormalsGet1
  END INTERFACE Basis_LineXiNormalsGet
  
  INTERFACE Basis_QuadratureGaussXiGet
    MODULE PROCEDURE Basis_QuadratureGaussXiGet0
    MODULE PROCEDURE Basis_QuadratureGaussXiGet1
  END INTERFACE Basis_QuadratureGaussXiGet
  
  PUBLIC BASIS_LAGRANGE_HERMITE_TP_TYPE,BASIS_SIMPLEX_TYPE,BASIS_SERENDIPITY_TYPE,BASIS_AUXILLIARY_TYPE,BASIS_B_SPLINE_TP_TYPE, &
    & BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE,BASIS_EXTENDED_LAGRANGE_TP_TYPE, BASIS_RADIAL_TYPE

  PUBLIC BASIS_LINEAR_LAGRANGE_INTERPOLATION,BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,BASIS_CUBIC_LAGRANGE_INTERPOLATION, &
    & BASIS_CUBIC_HERMITE_INTERPOLATION,BASIS_QUADRATIC1_HERMITE_INTERPOLATION,BASIS_QUADRATIC2_HERMITE_INTERPOLATION, &
    & BASIS_LINEAR_SIMPLEX_INTERPOLATION,BASIS_QUADRATIC_SIMPLEX_INTERPOLATION,BASIS_CUBIC_SIMPLEX_INTERPOLATION, &
    & BASIS_GAUSSIAN_RADIAL_INTERPOLATION, BASIS_MULTIQUARTIC_RADIAL_INTERPOLATION

  PUBLIC BASIS_LAGRANGE_INTERPOLATION,BASIS_HERMITE_INTERPOLATION,BASIS_SIMPLEX_INTERPOLATION,BASIS_SERENDIPITY_INTERPOLATION, &
    & BASIS_TRANSITION_INTERPOLATION,BASIS_SINGULAR_INTERPOLATION,BASIS_FOURIER_INTERPOLATION,BASIS_RADIAL_INTERPOLATION

  PUBLIC BASIS_LINEAR_INTERPOLATION_ORDER,BASIS_QUADRATIC_INTERPOLATION_ORDER,BASIS_CUBIC_INTERPOLATION_ORDER, &
    & BASIS_QUADRATIC1_INTERPOLATION_ORDER,BASIS_QUADRATIC2_INTERPOLATION_ORDER
  
  PUBLIC BASIS_GAUSS_LEGENDRE_QUADRATURE,BASIS_GAUSS_LAGUERRE_QUADRATURE,BASIS_GUASS_HERMITE_QUADRATURE,&
    & BASIS_ADAPTIVE_GAUSS_LEGENDRE_QUADRATURE,BASIS_GAUSS_SIMPLEX_QUADRATURE

  PUBLIC BASIS_XI_COLLAPSED,BASIS_COLLAPSED_AT_XI0,BASIS_COLLAPSED_AT_XI1,BASIS_NOT_COLLAPSED

  PUBLIC BASIS_NO_BOUNDARY_XI,BASIS_LINE_BOUNDARY_XI,BASIS_FACE_BOUNDARY_XI

  PUBLIC BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES,BASIS_DEFAULT_QUADRATURE_SCHEME,BASIS_LOW_QUADRATURE_SCHEME, &
    & BASIS_MID_QUADRATURE_SCHEME,BASIS_HIGH_QUADRATURE_SCHEME

  PUBLIC Basis_AssertIsFinished,Basis_AssertNotFinished

  PUBLIC Basis_AssertIsLHTPBasis,Basis_AssertIsSimplexBasis

  PUBLIC Basis_BasisFunctionsGet

  PUBLIC Basis_CollapsedXiGet

  PUBLIC Basis_ContextGet

  PUBLIC Basis_ElementParameterGet

  PUBLIC Basis_FaceBasisGet

  PUBLIC Basis_FaceElementParameterGet

  PUBLIC Basis_FaceNodeNumberGet

  PUBLIC Basis_FaceNodeDerivativeNumberGet

  PUBLIC Basis_FaceNodeNumberOfDerivativesGet

  PUBLIC Basis_FaceNumberOfNodesGet

  PUBLIC Basis_FaceXiNormalGet
  
  PUBLIC Basis_FamilyNumberFind

  PUBLIC Basis_Get

  PUBLIC Basis_InterpolationOrderGet
  
  PUBLIC Basis_InterpolationTypeGet

  PUBLIC Basis_InterpolationXiGet

  PUBLIC Basis_LineBasisGet

  PUBLIC Basis_LineElementParameterGet

  PUBLIC Basis_LineNodeNumberGet

  PUBLIC Basis_LineNodeDerivativeNumberGet

  PUBLIC Basis_LineNumberOfNodesGet

  PUBLIC Basis_LineXiNormalsGet

  PUBLIC Basis_LocalFaceNumberGet

  PUBLIC Basis_LocalLineNumberGet

  PUBLIC Basis_MaximumNumberOfDerivativesGet

  PUBLIC Basis_NodeNumberOfDerivativesGet

  PUBLIC Basis_NumberOfElementParametersGet

  PUBLIC Basis_NumberOfLocalFacesGet

  PUBLIC Basis_NumberOfLocalLinesGet

  PUBLIC Basis_NumberOfLocalNodesGet

  PUBLIC Basis_NumberOfNodesXiCGet

  PUBLIC Basis_NumberOfXiGet

  PUBLIC Basis_PartialDerivativeGet

  PUBLIC Basis_QuadratureGaussXiGet

  PUBLIC Basis_QuadratureNumberOfGaussXiGet

  PUBLIC Basis_QuadratureOrderGet

  PUBLIC Basis_QuadratureSchemeGet

  PUBLIC Basis_QuadratureTypeGet

  PUBLIC Basis_TypeGet

  PUBLIC Basis_UserNumberFind

  PUBLIC Basis_UserNumberGet

  PUBLIC BasisQuadratureScheme_GaussBasisFunctionGet

  PUBLIC BasisQuadratureScheme_GaussPositionGet

  PUBLIC BasisQuadratureScheme_GaussWeightGet

  PUBLIC BasisQuadratureScheme_NumberOfGaussGet
  

CONTAINS

  !
  !=================================================================================================================================
  !

  !>Assert that a basis has been finished
  SUBROUTINE Basis_AssertIsFinished(basis,err,error,*)

    !Argument Variables
    TYPE(BasisType), POINTER, INTENT(INOUT) :: basis !<The basis to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Basis_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
#endif    

    IF(.NOT.basis%basisFinished) THEN
      localError="Basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))// &
        & " has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Basis_AssertIsFinished")
    RETURN
999 ERRORSEXITS("Basis_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a basis has not been finished
  SUBROUTINE Basis_AssertNotFinished(basis,err,error,*)

    !Argument Variables
    TYPE(BasisType), POINTER, INTENT(INOUT) :: basis !<The basis to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Basis_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
#endif    

    IF(basis%basisFinished) THEN
      localError="Basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))// &
        & " has already been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Basis_AssertNotFinished")
    RETURN
999 ERRORSEXITS("Basis_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_AssertNotFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a basis is a LHTP basis
  SUBROUTINE Basis_AssertIsLHTPBasis(basis,err,error,*)

    !Argument Variables
    TYPE(BasisType), POINTER, INTENT(INOUT) :: basis !<The basis to assert the LHTP status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Basis_AssertIsLHTPBasis",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
#endif    

    IF(basis%TYPE/=BASIS_LAGRANGE_HERMITE_TP_TYPE) THEN
      localError="Basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))// &
        & " is not a LHTP basis."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Basis_AssertIsLHTPBasis")
    RETURN
999 ERRORSEXITS("Basis_AssertIsLHTPBasis",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_AssertIsLHTPBasis

  !
  !=================================================================================================================================
  !

  !>Assert that a basis is a simplex basis
  SUBROUTINE Basis_AssertIsSimplexBasis(basis,err,error,*)

    !Argument Variables
    TYPE(BasisType), POINTER, INTENT(INOUT) :: basis !<The basis to assert the simplex status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Basis_AssertIsSimplexBasis",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
#endif    

    IF(basis%TYPE/=BASIS_SIMPLEX_TYPE) THEN
      localError="Basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))// &
        & " is not a simplex basis."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Basis_AssertIsSimplexBasis")
    RETURN
999 ERRORSEXITS("Basis_AssertIsSimplexBasis",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_AssertIsSimplexBasis

  !
  !================================================================================================================================
  !

  !>Returns the basis functionsfor the basis. 
  SUBROUTINE Basis_BasisFunctionsGet(basis,basisFunctions,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the quadrature scheme for
    TYPE(BasisFunctionsType), POINTER :: basisFunctions !<On return, the basis fucntions for the basis. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("Basis_BasisFunctionsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(basisFunctions)) CALL FlagError("Basis functions is already associated.",err,error,*998)
    CALL Basis_AssertIsFinished(basis,err,error,*998)
#endif    
    
    basisFunctions=>basis%basisFunctions

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(basisFunctions)) THEN
      localError="Basis functions is not associated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    EXITS("Basis_BasisFunctionsGet")
    RETURN
999 NULLIFY(basisFunctions)
998 ERRORSEXITS("Basis_BasisFunctionsGet",err,error)    
    RETURN 1
    
  END SUBROUTINE Basis_BasisFunctionsGet

  !
  !================================================================================================================================
  !
  
  !>Gets the collapsed xi flags for a basis is identified by a a pointer. \see OpenCMISS::Iron::cmfe_Basis_CollapsedXiGet
  SUBROUTINE Basis_CollapsedXiGet(basis,collapsedXi,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: collapsedXi(:) !<collapsedXi(xiIdx). On return, the collapse parameter for each Xi direction. \see BasisRoutines_XiCollapse
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("Basis_CollapsedXiGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    CALL Basis_AssertIsFinished(basis,err,error,*999)
    IF(SIZE(collapsedXi,1)<SIZE(basis%collapsedXi)) THEN
      localError="The size of collapsed xi is too small. The supplied size is "// &
        & TRIM(NumberToVString(SIZE(collapsedXi,1),"*",err,error))//" and it needs to be >= "// &
        & TRIM(NumberToVString(SIZE(basis%collapsedXi,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    collapsedXi(1:basis%numberOfXi)=basis%collapsedXi(1:basis%numberOfXi)
    
    EXITS("Basis_CollapsedXiGet")
    RETURN
999 ERRORSEXITS("Basis_CollapsedXiGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_CollapsedXiGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the context for a basis.
  SUBROUTINE Basis_ContextGet(basis,context,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the context for
    TYPE(ContextType), POINTER :: context !<On exit, a pointer to the context for the basis. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Basis_ContextGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(context)) CALL FlagError("Context is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(basis%basisFunctions)) THEN
      localError="Basis functions is not associated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    context=>basis%basisFunctions%context

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(context)) THEN
      localError="The context is not associated for the basis functions for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("Basis_ContextGet")
    RETURN
999 NULLIFY(context)
998 ERRORSEXITS("Basis_ContextGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_ContextGet

  !
  !================================================================================================================================
  !

  !>Returns the basis element parameter for a local node and derivative in a basis.
  SUBROUTINE Basis_ElementParameterGet(basis,localDerivativeIdx,localNodeIdx,elementParameterIdx,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the element parameter for
    INTEGER(INTG), INTENT(IN) :: localDerivativeIdx !<The local derivative index to get the element parameter for
    INTEGER(INTG), INTENT(IN) :: localNodeIdx !<The local node index to get the element parameter for
    INTEGER(INTG), INTENT(OUT) :: elementParameterIdx !<On exit, the element parameter index for the local node and derivative in the basis.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Basis_ElementParameterGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(localNodeIdx<1.OR.localNodeIdx>basis%numberOfNodes) THEN
      localError="The specified local node index of "//TRIM(NumberToVString(localNodeIdx,"*",err,error))// &
        & " is invalid. The local node index should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfNodes,"*",err,error))//" for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%numberOfDerivatives)) THEN
      localError="The number of derivatives array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localDerivativeIdx<1.OR.localDerivativeIdx>basis%numberOfDerivatives(localNodeIdx)) THEN
      localError="The specified local derivative index of "//TRIM(NumberToVString(localDerivativeIdx,"*",err,error))// &
        & " is invalid. The local derivative index should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfDerivatives(localNodeIdx),"*",err,error))// &
        & " for local node index "//TRIM(NumberToVString(localNodeIdx,"*",err,error))//" of basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%elementParameterIndex)) THEN
      localError="The element parameter index array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    elementParameterIdx=basis%elementParameterIndex(localDerivativeIdx,localNodeIdx)
    
    EXITS("Basis_ElementParameterGet")
    RETURN
999 ERRORSEXITS("Basis_ElementParameterGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_ElementParameterGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the basis for a local face in a basis.
  SUBROUTINE Basis_FaceBasisGet(basis,localFaceNumber,faceBasis,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the face basis for
    INTEGER(INTG), INTENT(IN) :: localFaceNumber !<The local face number to get the face basis for
    TYPE(BasisType), POINTER :: faceBasis !<On exit, a pointer to the basis for the specified local face. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Basis_FaceBasisGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(faceBasis)) CALL FlagError("Face Basis is already associated.",err,error,*998)
    IF(basis%numberOfXi<2) THEN
      localError="Can only get the face basis for a basis with > 1 xi directions. The specified basis number of "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//" has "//TRIM(NumberToVString(basis%numberOfXi,"*",err,error))// &
        & " xi directions."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(localFaceNumber<0.OR.localFaceNumber>basis%numberOfLocalFaces) THEN
      localError="The specified local face number of "//TRIM(NumberToVString(localFaceNumber,"*",err,error))// &
        & " is invalid. The local face number should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfLocalFaces,"*",err,error))//" for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%localFaceBasis)) THEN
      localError="The local face basis array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    faceBasis=>basis%localFaceBasis(localFaceNumber)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(faceBasis)) THEN
      localError="The face basis is not associated for local face number "// &
        & TRIM(NumberToVString(localFaceNumber,"*",err,error))//" of basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("Basis_FaceBasisGet")
    RETURN
999 NULLIFY(faceBasis)
998 ERRORSEXITS("Basis_FaceBasisGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_FaceBasisGet

  !
  !================================================================================================================================
  !

  !>Returns the basis element parameter for a face parameter in a local face in a basis.
  SUBROUTINE Basis_FaceElementParameterGet(basis,faceParameterIdx,localFaceNumber,elementParameterIdx,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the element parameter for
    INTEGER(INTG), INTENT(IN) :: faceParameterIdx !<The face parameter index to get the equivalent element parameter for
    INTEGER(INTG), INTENT(IN) :: localFaceNumber !<The local face number to get the element parameter for
    INTEGER(INTG), INTENT(OUT) :: elementParameterIdx !<On exit, the equivalent element parameter index for the face parameter index on the local face.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Basis_FaceElementParameterGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(basis%numberOfXi<2) THEN
      localError="Can only get the face element parameters for a basis with > 1 xi directions. The specified basis number of "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//" has "//TRIM(NumberToVString(basis%numberOfXi,"*",err,error))// &
        & " xi directions."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
    IF(localFaceNumber<0.OR.localFaceNumber>basis%numberOfLocalFaces) THEN
      localError="The specified local face number of "//TRIM(NumberToVString(localFaceNumber,"*",err,error))// &
        & " is invalid. The local face number should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfLocalFaces,"*",err,error))//" for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%elementParametersInLocalFace)) THEN
      localError="The element parameters in local face array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(faceParameterIdx<1.OR.faceParameterIdx>basis%elementParametersInLocalFace(0,localFaceNumber)) THEN
      localError="The specified face parameter index of "//TRIM(NumberToVString(faceParameterIdx,"*",err,error))// &
        & " is invalid for local face number "//TRIM(NumberToVString(localFaceNumber,"*",err,error))// &
        & " of basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))// &
        & ". The face parameter index should be >=1 and <= "// &
        & TRIM(NumberToVString(basis%elementParametersInLocalFace(0,localFaceNumber),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    elementParameterIdx=basis%elementParametersInLocalFace(faceParameterIdx,localFaceNumber)
    
    EXITS("Basis_FaceElementParameterGet")
    RETURN
999 ERRORSEXITS("Basis_FaceElementParameterGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_FaceElementParameterGet

  !
  !================================================================================================================================
  !

  !>Returns the element local node number for a face node number in a local face in a basis.
  SUBROUTINE Basis_FaceNodeNumberGet(basis,localFaceNodeIdx,localFaceNumber,localElementNodeIdx,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the node number for
    INTEGER(INTG), INTENT(IN) :: localFaceNodeIdx !<The face node index to get the equivalent element node index for
    INTEGER(INTG), INTENT(IN) :: localFaceNumber !<The local face number to get the element node index for
    INTEGER(INTG), INTENT(OUT) :: localElementNodeIdx !<On exit, the equivalent element node index for the face node index on the local face.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Basis_FaceNodeNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(basis%numberOfXi<2) THEN
      localError="Can only get the face node numbers for a basis with > 1 xi directions. The specified basis number of "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//" has "//TRIM(NumberToVString(basis%numberOfXi,"*",err,error))// &
        & " xi directions."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
    IF(localFaceNumber<0.OR.localFaceNumber>basis%numberOfLocalFaces) THEN
      localError="The specified local face number of "//TRIM(NumberToVString(localFaceNumber,"*",err,error))// &
        & " is invalid. The local face number should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfLocalFaces,"*",err,error))//" for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%nodeNumbersInLocalFace)) THEN
      localError="The node numbers in local face array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%numberOfNodesInLocalFace)) THEN
      localError="The number of nodes in local face array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localFaceNodeIdx<1.OR.localFaceNodeIdx>basis%numberOfNodesInLocalFace(localFaceNumber)) THEN
      localError="The specified local face node index of "//TRIM(NumberToVString(localFaceNodeIdx,"*",err,error))// &
        & " is invalid for local face number "//TRIM(NumberToVString(localFaceNumber,"*",err,error))// &
        & " of basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))// &
        & ". The local face node index should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfNodesInLocalFace(localFaceNumber),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    localElementNodeIdx=basis%nodeNumbersInLocalFace(localFaceNodeIdx,localFaceNumber)
    
    EXITS("Basis_FaceNodeNumberGet")
    RETURN
999 ERRORSEXITS("Basis_FaceNodeNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_FaceNodeNumberGet

  !
  !================================================================================================================================
  !

  !>Returns the element local derivative index for a face derivative and node index in a local face in a basis.
  SUBROUTINE Basis_FaceNodeDerivativeNumberGet(basis,localFaceDerivativeIdx,localFaceNodeIdx,localFaceNumber, &
    & localElementDerivativeIdx,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the element derivative index for
    INTEGER(INTG), INTENT(IN) :: localFaceDerivativeIdx !<The face derivative index to get the equivalent element derivative index for
    INTEGER(INTG), INTENT(IN) :: localFaceNodeIdx !<The face node index to get the element derivative index for
    INTEGER(INTG), INTENT(IN) :: localFaceNumber !<The local face number to get the element derivative index for
    INTEGER(INTG), INTENT(OUT) :: localElementDerivativeIdx !<On exit, the equivalent element derivative index for the face node and derivative index on the local face.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Basis_FaceNodeDerivativeNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(basis%numberOfXi<2) THEN
      localError="Can only get the face node numbers for a basis with > 1 xi directions. The specified basis number of "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//" has "//TRIM(NumberToVString(basis%numberOfXi,"*",err,error))// &
        & " xi directions."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
    IF(localFaceNumber<0.OR.localFaceNumber>basis%numberOfLocalFaces) THEN
      localError="The specified local face number of "//TRIM(NumberToVString(localFaceNumber,"*",err,error))// &
        & " is invalid. The local face number should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfLocalFaces,"*",err,error))//" for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%derivativeNumbersInLocalFace)) THEN
      localError="The derivative numbers in local face array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%numberOfNodesInLocalFace)) THEN
      localError="The number of nodes in local face array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localFaceNodeIdx<1.OR.localFaceNodeIdx>basis%numberOfNodesInLocalFace(localFaceNumber)) THEN
      localError="The specified local face node index of "//TRIM(NumberToVString(localFaceNodeIdx,"*",err,error))// &
        & " is invalid for local face number "//TRIM(NumberToVString(localFaceNumber,"*",err,error))// &
        & " of basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))// &
        & ". The local face node index should be >=1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfNodesInLocalFace(localFaceNumber),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localFaceDerivativeIdx<1.OR. &
      & localFaceDerivativeIdx>basis%derivativeNumbersInLocalFace(0,localFaceNodeIdx,localFaceNumber)) THEN
      localError="The specified local face derivative index of "//TRIM(NumberToVString(localFaceDerivativeIdx,"*",err,error))// &
        & " is invalid for local face node index "//TRIM(NumberToVString(localFaceNodeIdx,"*",err,error))// &
        & " of local face number "//TRIM(NumberToVString(localFaceNumber,"*",err,error))// &
        & " of basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))// &
        & ". The local face derivative index should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%derivativeNumbersInLocalFace(0,localFaceNodeIdx,localFaceNumber),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    localElementDerivativeIdx=basis%derivativeNumbersInLocalFace(localFaceDerivativeIdx,localFaceNodeIdx,localFaceNumber)
    
    EXITS("Basis_FaceNodeDerivativeNumberGet")
    RETURN
999 ERRORSEXITS("Basis_FaceNodeDerivativeNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_FaceNodeDerivativeNumberGet

  !
  !================================================================================================================================
  !

  !>Returns the number of face derivatives for a face  node index in a local face in a basis.
  SUBROUTINE Basis_FaceNodeNumberOfDerivativesGet(basis,localFaceNodeIdx,localFaceNumber,numberOfNodeDerivatives,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the number of node derivatives for
    INTEGER(INTG), INTENT(IN) :: localFaceNodeIdx !<The face node index to get the number of derivatives for
    INTEGER(INTG), INTENT(IN) :: localFaceNumber !<The local face number to get the number of derivatives for
    INTEGER(INTG), INTENT(OUT) :: numberOfNodeDerivatives !<On exit, the numbe of node derivatives for the face node index on the local face.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Basis_FaceNodeNumberOfDerivativesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(basis%numberOfXi<2) THEN
      localError="Can only get the face node numbers of derivatives for a basis with > 1 xi directions. The specified basis "// &
        & "number of "//TRIM(NumberToVString(basis%userNumber,"*",err,error))//" has "// &
        & TRIM(NumberToVString(basis%numberOfXi,"*",err,error))//" xi directions."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
    IF(localFaceNumber<0.OR.localFaceNumber>basis%numberOfLocalFaces) THEN
      localError="The specified local face number of "//TRIM(NumberToVString(localFaceNumber,"*",err,error))// &
        & " is invalid. The local face number should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfLocalFaces,"*",err,error))//" for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%derivativeNumbersInLocalFace)) THEN
      localError="The derivative numbers in local face array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%numberOfNodesInLocalFace)) THEN
      localError="The number of nodes in local face array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localFaceNodeIdx<1.OR.localFaceNodeIdx>basis%numberOfNodesInLocalFace(localFaceNumber)) THEN
      localError="The specified local face node index of "//TRIM(NumberToVString(localFaceNodeIdx,"*",err,error))// &
        & " is invalid for local face number "//TRIM(NumberToVString(localFaceNumber,"*",err,error))// &
        & " of basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))// &
        & ". The local face node index should be >=1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfNodesInLocalFace(localFaceNumber),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    numberOfNodeDerivatives=basis%derivativeNumbersInLocalFace(0,localFaceNodeIdx,localFaceNumber)
    
    EXITS("Basis_FaceNodeNumberOfDerivativesGet")
    RETURN
999 ERRORSEXITS("Basis_FaceNodeNumberOfDerivativesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_FaceNodeNumberOfDerivativesGet

  !
  !================================================================================================================================
  !

  !>Returns the number of nodes for a local face in a basis.
  SUBROUTINE Basis_FaceNumberOfNodesGet(basis,localFaceNumber,numberOfNodes,err,error,*)
    
    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the number of face nodes for
    INTEGER(INTG), INTENT(IN) :: localFaceNumber !<The local face number to get the number of nodes for
    INTEGER(INTG), INTENT(OUT) :: numberOfNodes !<On exit, the number of nodes for the local face.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Basis_FaceNumberOfNodesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(basis%numberOfXi<2) THEN
      localError="Can only get the face number of nodes for a basis with > 1 xi directions. The specified basis "// &
        & "number of "//TRIM(NumberToVString(basis%userNumber,"*",err,error))//" has "// &
        & TRIM(NumberToVString(basis%numberOfXi,"*",err,error))//" xi directions."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
    IF(localFaceNumber<0.OR.localFaceNumber>basis%numberOfLocalFaces) THEN
      localError="The specified local face number of "//TRIM(NumberToVString(localFaceNumber,"*",err,error))// &
        & " is invalid. The local face number should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfLocalFaces,"*",err,error))//" for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%numberOfNodesInLocalFace)) THEN
      localError="The number of nodes in local face array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    numberOfNodes=basis%numberOfNodesInLocalFace(localFaceNumber)
    
    EXITS("Basis_FaceNumberOfNodesGet")
    RETURN
999 ERRORSEXITS("Basis_FaceNumberOfNodesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_FaceNumberOfNodesGet

  !
  !================================================================================================================================
  !

  !>Returns the xi normal direction for a local face in a basis.
  SUBROUTINE Basis_FaceXiNormalGet(basis,localFaceNumber,faceXiNormal,err,error,*)
    
    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the xi normal direction for
    INTEGER(INTG), INTENT(IN) :: localFaceNumber !<The local face number to get the xi normal direction for
    INTEGER(INTG), INTENT(OUT) :: faceXiNormal !<On exit, the xi normal direction for the local face. \see Constants_ElementNormalXiDirections
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Basis_FaceXiNormalGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(basis%numberOfXi<2) THEN
      localError="Can only get the face xi normal for a basis with > 1 xi directions. The specified basis "// &
        & "number of "//TRIM(NumberToVString(basis%userNumber,"*",err,error))//" has "// &
        & TRIM(NumberToVString(basis%numberOfXi,"*",err,error))//" xi directions."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
    IF(localFaceNumber<0.OR.localFaceNumber>basis%numberOfLocalFaces) THEN
      localError="The specified local face number of "//TRIM(NumberToVString(localFaceNumber,"*",err,error))// &
        & " is invalid. The local face number should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfLocalFaces,"*",err,error))//" for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%localFaceXiNormal)) THEN
      localError="The local face xi normal array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    faceXiNormal=basis%localFaceXiNormal(localFaceNumber)
    
    EXITS("Basis_FaceXiNormalGet")
    RETURN
999 ERRORSEXITS("Basis_FaceXiNormalGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_FaceXiNormalGet

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the basis with the given user number and family number. If no basis with that
  !>number and family number exists then basis is returned nullified \see BasisAccessRoutines::Basis_UserNumberFind
  RECURSIVE SUBROUTINE Basis_FamilyNumberFind(basisFunctions,userNumber,familyNumber,basis,err,error,*)

    !Argument variables
    TYPE(BasisFunctionsType), POINTER :: basisFunctions !<The basis functions to find the user number. 
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the basis to find
    INTEGER(INTG), INTENT(IN) :: familyNumber !<The family number of the basis to find
    TYPE(BasisType), POINTER :: basis !<On exit, A pointer to the basis. If no basis with the specified user and family numbers can be found the pointer is not associated.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: basisIdx,subBasisIdx
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("Basis_FamilyNumberFind",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(basisFunctions)) CALL FlagError("Basis functions is not associated.",err,error,*999)
#endif    
    
    IF(ALLOCATED(basisFunctions%bases)) THEN
      BasisLoop: DO basisIdx=1,basisFunctions%numberOfBasisFunctions
#ifdef WITH_PRECHECKS        
        IF(.NOT.ASSOCIATED(basisFunctions%bases(basisIdx)%ptr)) THEN
          localError="The basis is not associated for basis index "//&
            & TRIM(NumberToVString(basisIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
#endif        
        IF(basisFunctions%bases(basisIdx)%ptr%userNumber==userNumber) THEN
          IF(familyNumber==0) THEN
            basis=>basisFunctions%bases(basisIdx)%ptr
            EXIT BasisLoop
          ELSE
!!TODO: \todo This only works for one level of sub-bases at the moment
            IF(ALLOCATED(basisFunctions%bases(basisIdx)%ptr%subBases)) THEN
              SubBasisLoop: DO subBasisIdx=1,basisFunctions%bases(basisIdx)%ptr%numberOfSubBases
#ifdef WITH_PRECHECKS                
                IF(.NOT.ASSOCIATED(basisFunctions%bases(basisIdx)%ptr%subBases(subBasisIdx)%ptr)) THEN
                  localError="The sub basis is not associated for sub-basis index "// &
                    & TRIM(NumberToVString(subBasisIdx,"*",err,error))//" of basis index "// &
                    & TRIM(NumberToVString(basisIdx,"*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
#endif                
                IF(basisFunctions%bases(basisIdx)%ptr%subBases(subBasisIdx)%ptr%familyNumber==familyNumber) THEN
                  basis=>basisFunctions%bases(basisIdx)%ptr%subBases(subBasisIdx)%ptr
                  EXIT BasisLoop
                ENDIF
              ENDDO SubBasisLoop !subBasisIdx
            ENDIF
          ENDIF
        ENDIF
      ENDDO BasisLoop !basisIdx
    ENDIF
        
    EXITS("Basis_FamilyNumberFind")
    RETURN
999 NULLIFY(basis)
998 ERRORSEXITS("Basis_FamilyNumberFind",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_FamilyNumberFind

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the basis with the given user number. 
  SUBROUTINE Basis_Get(basisFunctions,userNumber,basis,err,error,*)

    !Argument variables
    TYPE(BasisFunctionsType), POINTER :: basisFunctions !<The basis functions to get the user number for. 
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the basis to find
    TYPE(BasisType), POINTER :: basis !<On exit, a pointer to the basis with the specified user number if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("Basis_Get",err,error,*999)

    CALL Basis_UserNumberFind(basisFunctions,userNumber,basis,err,error,*999)
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(basis)) THEN
      localError="A basis with an user number of "//TRIM(NumberToVString(userNumber,"*",err,error))//" does not exist."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
  
    EXITS("Basis_Get")
    RETURN
999 ERRORSEXITS("Basis_Get",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_Get

  !
  !================================================================================================================================
  !
  
  !>Gets the interpolation order in each xic direction for a basis identified by a pointer.
  SUBROUTINE Basis_InterpolationOrderGet(basis,interpolationOrder,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the interpolation order for
    INTEGER(INTG), INTENT(OUT) :: interpolationOrder(:) !<interpolationOrder(xicIdx). On return, the interpolation order parameters for the xic'th direction
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("Basis_InterpolationOrderGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    CALL Basis_AssertIsFinished(basis,err,error,*999)
    IF(.NOT.ALLOCATED(basis%interpolationOrder)) CALL FlagError("The basis interpolation order is not allocated.",err,error,*999)
    IF(SIZE(interpolationOrder,1)<SIZE(basis%interpolationOrder,1)) THEN
      localError="The size of interpolation order is too small. The supplied size is "// &
        & TRIM(NumberToVString(SIZE(interpolationOrder,1),"*",err,error))//" and it needs to be >= "// &
        & TRIM(NumberToVString(SIZE(basis%interpolationOrder,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    interpolationOrder(1:SIZE(basis%interpolationOrder,1))=basis%interpolationOrder(1:SIZE(basis%interpolationOrder,1))
    
    EXITS("Basis_InterpolationOrderGet")
    RETURN
999 ERRORSEXITS("Basis_InterpolationOrderGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_InterpolationOrderGet

  !
  !================================================================================================================================
  !
  
  !>Gets the interpolation type in each xic direction for a basis identified by a pointer.
  SUBROUTINE Basis_InterpolationTypeGet(basis,interpolationType,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the interpolation type for
    INTEGER(INTG), INTENT(OUT) :: interpolationType(:) !<interpolationType(xicIdx). On return, the interpolation type parameters for the xic'th direction
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("Basis_InterpolationTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    CALL Basis_AssertIsFinished(basis,err,error,*999)
    IF(.NOT.ALLOCATED(basis%interpolationType)) CALL FlagError("The basis interpolation type is not allocated.",err,error,*999)
    IF(SIZE(interpolationType,1)<SIZE(basis%interpolationType,1)) THEN
      localError="The size of interpolation type is too small. The supplied size is "// &
        & TRIM(NumberToVString(SIZE(interpolationType,1),"*",err,error))//" and it needs to be >= "// &
        & TRIM(NumberToVString(SIZE(basis%interpolationType,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    interpolationType(1:SIZE(basis%interpolationType,1))=basis%interpolationType(1:SIZE(basis%interpolationType,1))
    
    EXITS("Basis_InterpolationTypeGet")
    RETURN
999 ERRORSEXITS("Basis_InterpolationTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_InterpolationTypeGet

  !
  !================================================================================================================================
  !
  
  !>Gets/changes the interpolation type in each xi direction for a basis identified by a pointer. \see OpenCMISS::Iron::cmfe_Basis_InterpolationXiGet
  SUBROUTINE Basis_InterpolationXiGet(basis,interpolationXi,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the interpolation xi
    INTEGER(INTG), INTENT(OUT) :: interpolationXi(:) !<interpolationXi(xiIdx). On return, the interpolation xi parameters for each xiIdx'th direction \see BasisRoutines_InterpolationSpecifications
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("Basis_InterpolationXiGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    CALL Basis_AssertIsFinished(basis,err,error,*999)
    IF(.NOT.ALLOCATED(basis%interpolationXi)) CALL FlagError("The basis interpolation xi is not allocated.",err,error,*999)
    IF(SIZE(interpolationXi,1)<SIZE(basis%interpolationXi,1)) THEN
      localError="The size of interpolation xi is too small. The supplied size is "// &
        & TRIM(NumberToVString(SIZE(interpolationXi,1),"*",err,error))//" and it needs to be >= "// &
        & TRIM(NumberToVString(SIZE(basis%interpolationXi,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    interpolationXi(1:SIZE(basis%interpolationXi,1))=basis%interpolationXi(1:SIZE(basis%interpolationXi,1))
    
    EXITS("Basis_InterpolationXiGet")
    RETURN
999 ERRORSEXITS("Basis_InterpolationXiGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_InterpolationXiGet  

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the basis for a local line in a basis.
  SUBROUTINE Basis_LineBasisGet(basis,localLineNumber,lineBasis,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the line basis for
    INTEGER(INTG), INTENT(IN) :: localLineNumber !<The local line number to get the line basis for
    TYPE(BasisType), POINTER :: lineBasis !<On exit, a pointer to the basis for the specified local line. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Basis_LineBasisGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(lineBasis)) CALL FlagError("Line Basis is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(localLineNumber<0.OR.localLineNumber>basis%numberOfLocalLines) THEN
      localError="The specified local line number of "//TRIM(NumberToVString(localLineNumber,"*",err,error))// &
        & " is invalid. The local line number should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfLocalLines,"*",err,error))//" for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%localLineBasis)) THEN
      localError="The local line basis array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    lineBasis=>basis%localLineBasis(localLineNumber)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(lineBasis)) THEN
      localError="The line basis is not associated for local line number "// &
        & TRIM(NumberToVString(localLineNumber,"*",err,error))//" of basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("Basis_LineBasisGet")
    RETURN
999 NULLIFY(lineBasis)
998 ERRORSEXITS("Basis_LineBasisGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_LineBasisGet

  !
  !================================================================================================================================
  !

  !>Returns the basis element parameter for a line parameter in a local line in a basis.
  SUBROUTINE Basis_LineElementParameterGet(basis,lineParameterIdx,localLineNumber,elementParameterIdx,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the element parameter for
    INTEGER(INTG), INTENT(IN) :: lineParameterIdx !<The line parameter index to get the equivalent element parameter for
    INTEGER(INTG), INTENT(IN) :: localLineNumber !<The local line number to get the element parameter for
    INTEGER(INTG), INTENT(OUT) :: elementParameterIdx !<On exit, the equivalent element parameter index for the line parameter index on the local line.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Basis_LineElementParameterGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(basis%numberOfXi<2) THEN
      localError="Can only get the line element parameters for a basis with > 1 xi directions. The specified basis number of "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//" has "//TRIM(NumberToVString(basis%numberOfXi,"*",err,error))// &
        & " xi directions."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
    IF(localLineNumber<0.OR.localLineNumber>basis%numberOfLocalLines) THEN
      localError="The specified local line number of "//TRIM(NumberToVString(localLineNumber,"*",err,error))// &
        & " is invalid. The local line number should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfLocalLines,"*",err,error))//" for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%elementParametersInLocalLine)) THEN
      localError="The element parameters in local line array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(lineParameterIdx<1.OR.lineParameterIdx>basis%elementParametersInLocalLine(0,localLineNumber)) THEN
      localError="The specified line parameter index of "//TRIM(NumberToVString(lineParameterIdx,"*",err,error))// &
        & " is invalid for local line number "//TRIM(NumberToVString(localLineNumber,"*",err,error))// &
        & " of basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))// &
        & ". The line parameter index should be >=1 and <= "// &
        & TRIM(NumberToVString(basis%elementParametersInLocalLine(0,localLineNumber),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    elementParameterIdx=basis%elementParametersInLocalLine(lineParameterIdx,localLineNumber)
    
    EXITS("Basis_LineElementParameterGet")
    RETURN
999 ERRORSEXITS("Basis_LineElementParameterGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_LineElementParameterGet

  !
  !================================================================================================================================
  !

  !>Returns the element local node number for a line node number in a local line in a basis.
  SUBROUTINE Basis_LineNodeNumberGet(basis,localLineNodeIdx,localLineNumber,localElementNodeIdx,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the node number for
    INTEGER(INTG), INTENT(IN) :: localLineNodeIdx !<The line node index to get the equivalent element node index for
    INTEGER(INTG), INTENT(IN) :: localLineNumber !<The local line number to get the element node index for
    INTEGER(INTG), INTENT(OUT) :: localElementNodeIdx !<On exit, the equivalent element node index for the line node index on the local line.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Basis_LineNodeNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(basis%numberOfXi<2) THEN
      localError="Can only get the line node numbers for a basis with > 1 xi directions. The specified basis number of "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//" has "//TRIM(NumberToVString(basis%numberOfXi,"*",err,error))// &
        & " xi directions."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
    IF(localLineNumber<0.OR.localLineNumber>basis%numberOfLocalLines) THEN
      localError="The specified local line number of "//TRIM(NumberToVString(localLineNumber,"*",err,error))// &
        & " is invalid. The local line number should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfLocalLines,"*",err,error))//" for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%nodeNumbersInLocalLine)) THEN
      localError="The node numbers in local line array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%numberOfNodesInLocalLine)) THEN
      localError="The number of nodes in local line array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localLineNodeIdx<1.OR.localLineNodeIdx>basis%numberOfNodesInLocalLine(localLineNumber)) THEN
      localError="The specified local line node index of "//TRIM(NumberToVString(localLineNodeIdx,"*",err,error))// &
        & " is invalid for local line number "//TRIM(NumberToVString(localLineNumber,"*",err,error))// &
        & " of basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))// &
        & ". The local line node index should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfNodesInLocalLine(localLineNumber),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    localElementNodeIdx=basis%nodeNumbersInLocalLine(localLineNodeIdx,localLineNumber)
    
    EXITS("Basis_LineNodeNumberGet")
    RETURN
999 ERRORSEXITS("Basis_LineNodeNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_LineNodeNumberGet

  !
  !================================================================================================================================
  !

  !>Returns the element local derivative index for a line derivative and node index in a local line in a basis.
  SUBROUTINE Basis_LineNodeDerivativeNumberGet(basis,localLineNodeIdx,localLineNumber,localLineDerivativeIdx,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the element derivative index for
    INTEGER(INTG), INTENT(IN) :: localLineNodeIdx !<The line node index to get the element derivative index for
    INTEGER(INTG), INTENT(IN) :: localLineNumber !<The local line number to get the element derivative index for
    INTEGER(INTG), INTENT(OUT) :: localLineDerivativeIdx !<On exit, the line derivative index 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Basis_LineNodeDerivativeNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(basis%numberOfXi<2) THEN
      localError="Can only get the line node numbers for a basis with > 1 xi directions. The specified basis number of "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//" has "//TRIM(NumberToVString(basis%numberOfXi,"*",err,error))// &
        & " xi directions."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
    IF(localLineNumber<0.OR.localLineNumber>basis%numberOfLocalLines) THEN
      localError="The specified local line number of "//TRIM(NumberToVString(localLineNumber,"*",err,error))// &
        & " is invalid. The local line number should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfLocalLines,"*",err,error))//" for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%derivativeNumbersInLocalLine)) THEN
      localError="The derivative numbers in local line array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%numberOfNodesInLocalLine)) THEN
      localError="The number of nodes in local line array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localLineNodeIdx<1.OR.localLineNodeIdx>basis%numberOfNodesInLocalLine(localLineNumber)) THEN
      localError="The specified local line node index of "//TRIM(NumberToVString(localLineNodeIdx,"*",err,error))// &
        & " is invalid for local line number "//TRIM(NumberToVString(localLineNumber,"*",err,error))// &
        & " of basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))// &
        & ". The local line node index should be >=1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfNodesInLocalLine(localLineNumber),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    localLineDerivativeIdx=basis%derivativeNumbersInLocalLine(localLineNodeIdx,localLineNumber)
    
    EXITS("Basis_LineNodeDerivativeNumberGet")
    RETURN
999 ERRORSEXITS("Basis_LineNodeDerivativeNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_LineNodeDerivativeNumberGet

  !
  !================================================================================================================================
  !

  !>Returns the number of nodes for a local line in a basis.
  SUBROUTINE Basis_LineNumberOfNodesGet(basis,localLineNumber,numberOfNodes,err,error,*)
    
    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the number of line nodes for
    INTEGER(INTG), INTENT(IN) :: localLineNumber !<The local line number to get the number of nodes for
    INTEGER(INTG), INTENT(OUT) :: numberOfNodes !<On exit, the number of nodes for the local line.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Basis_LineNumberOfNodesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(basis%numberOfXi<2) THEN
      localError="Can only get the line number of nodes for a basis with > 1 xi directions. The specified basis "// &
        & "number of "//TRIM(NumberToVString(basis%userNumber,"*",err,error))//" has "// &
        & TRIM(NumberToVString(basis%numberOfXi,"*",err,error))//" xi directions."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
    IF(localLineNumber<0.OR.localLineNumber>basis%numberOfLocalLines) THEN
      localError="The specified local line number of "//TRIM(NumberToVString(localLineNumber,"*",err,error))// &
        & " is invalid. The local line number should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfLocalLines,"*",err,error))//" for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%numberOfNodesInLocalLine)) THEN
      localError="The number of nodes in local line array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    numberOfNodes=basis%numberOfNodesInLocalLine(localLineNumber)
    
    EXITS("Basis_LineNumberOfNodesGet")
    RETURN
999 ERRORSEXITS("Basis_LineNumberOfNodesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_LineNumberOfNodesGet

  !
  !================================================================================================================================
  !

  !>Returns the xi normal direction for a local line in a 2D basis.
  SUBROUTINE Basis_LineXiNormalsGet0(basis,localLineNumber,lineXiNormal,err,error,*)
    
    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the xi normal direction for
    INTEGER(INTG), INTENT(IN) :: localLineNumber !<The local line number to get the xi normal direction for
    INTEGER(INTG), INTENT(OUT) :: lineXiNormal !<lineXiNormal. On exit, the xi normal direction for the local line. \see Constants_ElementNormalXiDirections
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: lineXiNormals(1)
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Basis_LineXiNormalsGet0",err,error,*999)

    CALL Basis_LineXiNormalsGet1(basis,localLineNumber,lineXiNormals,err,error,*999)    
    lineXiNormal=lineXiNormals(1)
    
    EXITS("Basis_LineXiNormalsGet0")
    RETURN
999 ERRORSEXITS("Basis_LineXiNormalsGet0",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_LineXiNormalsGet0

  !
  !================================================================================================================================
  !

  !>Returns the xi normal directions for a local line in a basis.
  SUBROUTINE Basis_LineXiNormalsGet1(basis,localLineNumber,lineXiNormals,err,error,*)
    
    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the xi normal direction for
    INTEGER(INTG), INTENT(IN) :: localLineNumber !<The local line number to get the xi normal direction for
    INTEGER(INTG), INTENT(OUT) :: lineXiNormals(:) !<lineXiNormals(xiCoordIdx). On exit, the xi normal direction for the local line. \see Constants_ElementNormalXiDirections
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Basis_LineXiNormalsGet1",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(basis%numberOfXi<2) THEN
      localError="Can only get the line xi normal for a basis with > 1 xi directions. The specified basis "// &
        & "number of "//TRIM(NumberToVString(basis%userNumber,"*",err,error))//" has "// &
        & TRIM(NumberToVString(basis%numberOfXi,"*",err,error))//" xi directions."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
    IF(localLineNumber<0.OR.localLineNumber>basis%numberOfLocalLines) THEN
      localError="The specified local line number of "//TRIM(NumberToVString(localLineNumber,"*",err,error))// &
        & " is invalid. The local line number should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfLocalLines,"*",err,error))//" for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%localLineXiNormals)) THEN
      localError="The local line xi normal array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(basis%localLineXiNormals,1)<(basis%numberOfXi-1)) THEN
      localError="The local line xi normal array is too small. The size of the xi normal array is "// &
        & TRIM(NumberToVString(SIZE(basis%localLineXiNormals,1),"*",err,error))//" and it needs to be >= "// &
        & TRIM(NumberToVString(basis%numberOfXi-1,"*",err,error))//" for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    lineXiNormals(1:basis%numberOfXi-1)=basis%localLineXiNormals(1:basis%numberOfXi-1,localLineNumber)
    
    EXITS("Basis_LineXiNormalsGet1")
    RETURN
999 ERRORSEXITS("Basis_LineXiNormalsGet1",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_LineXiNormalsGet1

  !
  !================================================================================================================================
  !

  !>Finds the local face number that corresponds to a normal xi direction for the basis. 
  SUBROUTINE Basis_LocalFaceNumberGet(basis,normalXiDirection,localFaceNumber,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the local face number for
    INTEGER(INTG), INTENT(IN) :: normalXiDirection !<The normal xi direction of the face.
    INTEGER(INTG), INTENT(OUT) :: localFaceNumber !<On exit, the local face number corresponding to the normal xi direction.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("Basis_LocalFaceNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    CALL Basis_AssertIsFinished(basis,err,error,*999)
    IF(normalXiDirection < -basis%numberOfXiCoordinates .OR. normalXiDirection > basis%numberOfXiCoordinates) THEN
      localError="The specified normal xi direction of "//TRIM(NumberToVString(normalXiDirection,"*",err,error))// &
        & " is invalid for basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))// &
        & ". The normal xi direction must be >= "// &
        & TRIM(NumberToVString(-basis%numberOfXiCoordinates,"*",err,error))// &
        & " and <= "//TRIM(NumberToVString(basis%numberOfXiCoordinates,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%xiNormalLocalFace)) THEN
      localError="The xi normal local face array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif   

    localFaceNumber=basis%xiNormalLocalFace(normalXiDirection)

#ifdef WITH_POSTCHECKS    
    IF(localFaceNumber<1.OR.localFaceNumber>basis%numberOfLocalFaces) THEN
      localError="The local face number of "//TRIM(NumberToVString(localFaceNumber,"*",err,error))//" for xi normal direction "// &
        & TRIM(NumberToVString(normalXiDirection,"*",err,error))//" of basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//" is invalid. The local face number should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfLocalFaces,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("Basis_LocalFaceNumberGet")
    RETURN
999 ERRORSEXITS("Basis_LocalFaceNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_LocalFaceNumberGet

  !
  !================================================================================================================================
  !

  !>Finds the local line number that corresponds to the normal xi directions for the basis. 
  SUBROUTINE Basis_LocalLineNumberGet(basis,normalXiDirections,localLineNumber,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the local face number for
    INTEGER(INTG), INTENT(IN) :: normalXiDirections(:) !<The normal xi directions of the line.
    INTEGER(INTG), INTENT(OUT) :: localLineNumber !<On exit, the local line number corresponding to the normal xi directions.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("Basis_LocalLineNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    CALL Basis_AssertIsFinished(basis,err,error,*999)
    IF(SIZE(normalXiDirections,1)<2) THEN
      localError="The specified number of normal xi directions of "// &
        & TRIM(NumberToVString(SIZE(normalXiDirections,1),"*",err,error))// &
        & " is invalid. There should be at least 2 normal xi directions."
      CALL FlagError(localError,err,error,*999)
    END IF
    IF(normalXiDirections(1) < -basis%numberOfXiCoordinates .OR. normalXiDirections(1) > basis%numberOfXiCoordinates) THEN
      localError="The first specified normal xi direction of "//TRIM(NumberToVString(normalXiDirections(1),"*",err,error))// &
        & " is invalid for basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))// &
        & ". The normal xi direction must be >= "//TRIM(NumberToVString(-basis%numberOfXiCoordinates,"*",err,error))// &
        & " and <= "//TRIM(NumberToVString(basis%numberOfXiCoordinates,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(normalXiDirections(2) < -basis%numberOfXiCoordinates .OR. normalXiDirections(2) > basis%numberOfXiCoordinates) THEN
      localError="The second specified normal xi direction of "//TRIM(NumberToVString(normalXiDirections(2),"*",err,error))// &
        & " is invalid for basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))// &
        & ". The normal xi direction must be >= "//TRIM(NumberToVString(-basis%numberOfXiCoordinates,"*",err,error))// &
        & " and <= "//TRIM(NumberToVString(basis%numberOfXiCoordinates,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%xiNormalsLocalLine)) THEN
      localError="The xi normals local line array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    localLineNumber=basis%xiNormalsLocalLine(normalXiDirections(1),normalXiDirections(2))

#ifdef WITH_POSTCHECKS    
    IF(localLineNumber<1.OR.localLineNumber>basis%numberOfLocalLines) THEN
      localError="The local line number of "//TRIM(NumberToVString(localLineNumber,"*",err,error))//" for xi normal directions "// &
        & TRIM(NumberToVString(normalXiDirections(1),"*",err,error))//","// &
        & TRIM(NumberToVString(normalXiDirections(2),"*",err,error))//" of basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//" is invalid. The local line number should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfLocalLines,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("Basis_LocalLineNumberGet")
    RETURN
999 ERRORSEXITS("Basis_LocalLineNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_LocalLineNumberGet

  !
  !================================================================================================================================
  !

  !>Returns the number of derivatives for a local nodes in the specified basis 
  SUBROUTINE Basis_NodeNumberOfDerivativesGet(basis,localNodeIdx,numberOfDerivatives,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the number of derivatives
    INTEGER(INTG), INTENT(IN) :: localNodeIdx !<The local node index to get the number of derivatives for
    INTEGER(INTG), INTENT(OUT) :: numberOfDerivatives !<On return, the number of derivatives at the local node in the basis
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("Basis_NodeNumberOfDerivativesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    CALL Basis_AssertIsFinished(basis,err,error,*999)
    IF(localNodeIdx<1.OR.localNodeIdx>basis%numberOfNodes) THEN
      localError="The specified local node index of "//TRIM(NumberToVString(localNodeIdx,"*",err,error))// &
        & " is invalid for basis number "//TRIM(NumberToVString(basis%userNumber,"*",err,error))// &
        & ". The local node index should be >= 1 and <= "//TRIM(NumberToVString(basis%numberOfNodes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)     
    ENDIF
    IF(.NOT.ALLOCATED(basis%numberOfDerivatives)) THEN
      localError="The number of derivatives array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
   
    numberOfDerivatives=basis%numberOfDerivatives(localNodeIdx)
    
    EXITS("Basis_NodeNumberOfDerivativesGet")
    RETURN
999 ERRORSEXITS("Basis_NodeNumberOfDerivativesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_NodeNumberOfDerivativesGet

  !
  !================================================================================================================================
  !

  !>Returns the maximum number of derivatives for the specified basis 
  SUBROUTINE Basis_MaximumNumberOfDerivativesGet(basis,maximumNumberOfDerivatives,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the maximum number of derivatives for
    INTEGER(INTG), INTENT(OUT) :: maximumNumberOfDerivatives !<On return, the maximum number of derivatives in the basis
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Basis_MaximumNumberOfDerivativesGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
#endif    
   
    maximumNumberOfDerivatives=basis%maximumNumberOfDerivatives
    
    EXITS("Basis_MaximumNumberOfDerivativesGet")
    RETURN
999 ERRORSEXITS("Basis_MaximumNumberOfDerivativesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_MaximumNumberOfDerivativesGet

  !
  !================================================================================================================================
  !

  !>Returns the number of element parameters in the specified basis
  SUBROUTINE Basis_NumberOfElementParametersGet(basis,numberOfElementParameters,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the number of element parameters for
    INTEGER(INTG), INTENT(OUT) :: numberOfElementParameters !<On return, the number of element parameters in the basis
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Basis_NumberOfElementParametersGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    CALL Basis_AssertIsFinished(basis,err,error,*999)
#endif    
   
    numberOfElementParameters=basis%numberOfElementParameters
    
    EXITS("Basis_NumberOfElementParametersGet")
    RETURN
999 ERRORSEXITS("Basis_NumberOfElementParametersGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_NumberOfElementParametersGet

  !
  !================================================================================================================================
  !

  !>Returns the number of local faces in the specified basis
  SUBROUTINE Basis_NumberOfLocalFacesGet(basis,numberOfLocalFaces,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the number of faces for
    INTEGER(INTG), INTENT(OUT) :: numberOfLocalFaces !<On return, the number of local faces in the basis
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Basis_NumberOfLocalFacesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    CALL Basis_AssertIsFinished(basis,err,error,*999)
#endif    
   
    numberOfLocalFaces=basis%numberOfLocalFaces
    
    EXITS("Basis_NumberOfLocalFacesGet")
    RETURN
999 ERRORSEXITS("Basis_NumberOfLocalFacesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_NumberOfLocalFacesGet

  !
  !================================================================================================================================
  !

  !>Returns the number of local lines in the specified basis
  SUBROUTINE Basis_NumberOfLocalLinesGet(basis,numberOfLocalLines,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the number of lines for
    INTEGER(INTG), INTENT(OUT) :: numberOfLocalLines !<On return, the number of local lines in the basis
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Basis_NumberOfLocalLinesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    CALL Basis_AssertIsFinished(basis,err,error,*999)
#endif    
   
    numberOfLocalLines=basis%numberOfLocalLines
    
    EXITS("Basis_NumberOfLocalLinesGet")
    RETURN
999 ERRORSEXITS("Basis_NumberOfLocalLinesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_NumberOfLocalLinesGet

  !
  !================================================================================================================================
  !

  !>Returns the number of local nodes in the specified basis \see OpenCMISS::Iron::cmfe_Basis_NumberOfLocalNodesGet
  SUBROUTINE Basis_NumberOfLocalNodesGet(basis,numberOfLocalNodes,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the number of nodes
    INTEGER(INTG), INTENT(OUT) :: numberOfLocalNodes !<On return, the number of local nodes in the basis
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Basis_NumberOfLocalNodesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    CALL Basis_AssertIsFinished(basis,err,error,*999)
#endif    
   
    numberOfLocalNodes=basis%numberOfNodes
    
    EXITS("Basis_NumberOfLocalNodesGet")
    RETURN
999 ERRORSEXITS("Basis_NumberOfLocalNodesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_NumberOfLocalNodesGet

  !
  !================================================================================================================================
  !

  !>Returns the number of nodes in each xi coordinate direction in the specified basis 
  SUBROUTINE Basis_NumberOfNodesXiCGet(basis,numberOfNodesXiC,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the number of nodes in each xi coordinate diretion
    INTEGER(INTG), INTENT(OUT) :: numberOfNodesXiC(:) !<numberOfNodesXiC(xicIdx). On return, the number of nodes in each xi coordinate direction in the basis
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("Basis_NumberOfNodesXiCGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    CALL Basis_AssertIsFinished(basis,err,error,*999)
    IF(SIZE(numberOfNodesXiC,1)<basis%numberOfXiCoordinates) THEN
      localError="The size of the specified number of nodes xi array is too small for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//". The size should be >= "// &
        & TRIM(NumberToVString(basis%numberOfXiCoordinates,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%numberOfNodesXic)) THEN
      localError="The number of nodes xi coordinate array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
   
    numberOfNodesXiC(1:basis%numberOfXiCoordinates)=basis%numberOfNodesXic(1:basis%numberOfXiCoordinates)
    
    EXITS("Basis_NumberOfNodesXiCGet")
    RETURN
999 ERRORSEXITS("Basis_NumberOfNodesXiCGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_NumberOfNodesXiCGet

  !
  !================================================================================================================================
  !
  
  !>Gets the number of xi directions for a basis. \see OpenCMISS::Iron::cmfe_Basis_NumberOfXiGet
  SUBROUTINE Basis_NumberOfXiGet(basis,numberOfXi,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis function to get the number of xi for
    INTEGER(INTG), INTENT(OUT) :: numberOfXi !<On return, the number of Xi directions for the specified basis.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Basis_NumberOfXiGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    CALL Basis_AssertIsFinished(basis,err,error,*999)
#endif    
    
    numberOfXi=basis%numberOfXi
   
    EXITS("Basis_NumberOfXiGet")
    RETURN
999 ERRORSEXITS("Basis_NumberOfXiGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_NumberOfXiGet

  !
  !================================================================================================================================
  !

  !>Returns the basis partial derivative index for a local node and derivative in a basis.
  SUBROUTINE Basis_PartialDerivativeGet(basis,localDerivativeIdx,localNodeIdx,partialDerivativeIdx,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the partial derivative for
    INTEGER(INTG), INTENT(IN) :: localDerivativeIdx !<The local derivative index to get the partial derivative for
    INTEGER(INTG), INTENT(IN) :: localNodeIdx !<The local node index to get the partial derivative for
    INTEGER(INTG), INTENT(OUT) :: partialDerivativeIdx !<On exit, the partial derivative index for the local node and derivative in the basis.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Basis_PartialDerivativeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(localNodeIdx<1.OR.localNodeIdx>basis%numberOfNodes) THEN
      localError="The specified local node index of "//TRIM(NumberToVString(localNodeIdx,"*",err,error))// &
        & " is invalid. The local node index should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfNodes,"*",err,error))//" for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%numberOfDerivatives)) THEN
      localError="The number of derivatives array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localDerivativeIdx<1.OR.localDerivativeIdx>basis%numberOfDerivatives(localNodeIdx)) THEN
      localError="The specified local derivative index of "//TRIM(NumberToVString(localDerivativeIdx,"*",err,error))// &
        & " is invalid. The local derivative index should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfDerivatives(localNodeIdx),"*",err,error))// &
        & " for local node index "//TRIM(NumberToVString(localNodeIdx,"*",err,error))//" of basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(basis%partialDerivativeIndex)) THEN
      localError="The partial derivative index array is not allocated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    partialDerivativeIdx=basis%partialDerivativeIndex(localDerivativeIdx,localNodeIdx)
    
    EXITS("Basis_PartialDerivativeGet")
    RETURN
999 ERRORSEXITS("Basis_PartialDerivativeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_PartialDerivativeGet

  !
  !================================================================================================================================
  !
  
  !>Returns the xi positions of a Gauss point on a basis quadrature identified by a pointer. \see OpenCMISS::Iron::cmfe_Basis_QuadratureGaussXiGet
  SUBROUTINE Basis_QuadratureGaussXiGet0(basis,quadratureSchemeIdx,gaussPointNumber,gaussXi,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: quadratureSchemeIdx !<The quadrature scheme to return the Gauss points for
    INTEGER(INTG), INTENT(IN) :: gaussPointNumber !<The Gauss point to return the element xi position for.
    REAL(DP), INTENT(OUT) :: gaussXi(:) !<On return, gaussXi(xiIdx) the xi position of the specified Gauss point for the specified quadrature scheme.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: gaussPointNumberArray(1)
    REAL(DP) :: gaussXiArray(SIZE(gaussXi,1),1)
    
    ENTERS("Basis_QuadratureGaussXiGet0",err,error,*999)

    gaussPointNumberArray(1)=gaussPointNumber
    CALL Basis_QuadratureGaussXiGet1(basis,quadratureSchemeIdx,gaussPointNumberArray,gaussXiArray,err,error,*999)
    gaussXi(:)=gaussXiArray(:,1)
      
    EXITS("Basis_QuadratureGaussXiGet0")
    RETURN
999 ERRORSEXITS("Basis_QuadratureGaussXiGet0",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_QuadratureGaussXiGet0

  !
  !================================================================================================================================
  !
  
  !>Returns the xi positions of Gauss points on a basis quadrature identified by a pointer. If no Gauss points are specified then xi positions of all Gauss points are returned. \see OpenCMISS::Iron::cmfe_Basis_QuadratureGaussXiGet
  SUBROUTINE Basis_QuadratureGaussXiGet1(basis,quadratureSchemeIdx,gaussPointNumbers,gaussXi,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis
    INTEGER(INTG), INTENT(IN) :: quadratureSchemeIdx !<The quadrature scheme to return the Gauss points for
    INTEGER(INTG), INTENT(IN) :: gaussPointNumbers(:) !<gaussPointNumbers(gaussPointIdx). The Gauss points to return the element xi positions for.
    REAL(DP), INTENT(OUT) :: gaussXi(:,:) !<On return, gaussXi(xiIdx,gaussPointIdx) the Gauss xi positions for the specified quadrature scheme.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: gaussPointIdx
    TYPE(QuadratureSchemeType), POINTER :: quadratureScheme
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("Basis_QuadratureGaussXiGet1",err,error,*999)

#ifdef WITH_PRECHECKS    
    CALL Basis_AssertIsFinished(basis,err,error,*999)
    IF(.NOT.ASSOCIATED(basis%quadrature%basis)) CALL FlagError("Quadrature basis is not associated.",err,error,*999)
    IF(SIZE(gaussXi,1)<basis%numberOfXi) THEN
      localError="The number of xi values to return of "//TRIM(NumberToVString(SIZE(gaussXi,1),"*",err,error))// &
        & " is invalid and needs to be >= "//TRIM(NumberToVString(basis%numberOfXi,"*",err,error))//" for the specified basis."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    NULLIFY(quadratureScheme)
    CALL Basis_QuadratureSchemeGet(basis,quadratureSchemeIdx,quadratureScheme,err,error,*999)    
    IF(SIZE(gaussPointNumbers,1)==0) THEN !Return all Gauss point xi locations.
#ifdef WITH_PRECHECKS      
      IF(SIZE(gaussXi,2)<quadratureScheme%numberOfGauss) THEN
        localError="The number of Gauss Points to return the xi values for of "// &
          & TRIM(NumberToVString(SIZE(gaussXi,2),"*",err,error))//" is invalid and needs to be >= "// &
          & TRIM(NumberToVString(quadratureScheme%numberOfGauss,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
#endif      
      gaussXi(1:basis%numberOfXi,quadratureScheme%numberOfGauss)= &
        & quadratureScheme%gaussPositions(1:basis%numberOfXi,quadratureScheme%numberOfGauss)
    ELSE !Return only specified Gauss point xi locations.
#ifdef WITH_PRECHECKS      
      IF(SIZE(gaussXi,2)<SIZE(gaussPointNumbers,1)) THEN
        localError="The number of Gauss Points to return the xi values for of "// &
          & TRIM(NumberToVString(SIZE(gaussXi,2),"*",err,error))//" is invalid and needs to be >= "// &
          & TRIM(NumberToVString(SIZE(gaussPointNumbers,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
#endif      
      DO gaussPointIdx=1,SIZE(gaussPointNumbers)
#ifdef WITH_PRECHECKS        
        IF(gaussPointNumbers(gaussPointIdx)<1.OR.gaussPointNumbers(gaussPointIdx)>quadratureScheme%numberOfGauss) THEN
          localError="The specified Gauss point number of "// &
            & TRIM(NumberToVString(gaussPointNumbers(gaussPointIdx),"*",err,error))// &
            & " at position index "//TRIM(NumberToVString(gaussPointIdx,"*",err,error))// &
            & " is invalid for quadrature scheme index "// &
            & TRIM(NumberToVString(quadratureSchemeIdx,"*",err,error))//" of basis number "// &
            & TRIM(NumberToVString(basis%userNumber,"*",err,error))//". The Gauss point number should be >=1 and <= "// &
            & TRIM(NumberToVString(quadratureScheme%numberOfGauss,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
#endif        
        gaussXi(1:basis%numberOfXi,gaussPointIdx)= &
          & quadratureScheme%gaussPositions(1:basis%numberOfXi,gaussPointNumbers(gaussPointIdx))
      ENDDO !gaussPointIdx
    ENDIF
      
    EXITS("Basis_QuadratureGaussXiGet1")
    RETURN
999 ERRORSEXITS("Basis_QuadratureGaussXiGet1",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_QuadratureGaussXiGet1

  !
  !================================================================================================================================
  !
  
  !>Get the number of Gauss points in each xi direction on a basis quadrature identified by a pointer. \see OpenCMISS::Iron::cmfe_Basis_QuadratureNumberOfGaussXiGet
  SUBROUTINE Basis_QuadratureNumberOfGaussXiGet(basis,quadratureNumberOfGaussXi,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: quadratureNumberOfGaussXi(:) !<On return, the number of Gauss in each Xi direction
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("Basis_QuadratureNumberOfGaussXiGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    CALL Basis_AssertIsFinished(basis,err,error,*999)
    IF(.NOT.ASSOCIATED(basis%quadrature%basis)) CALL FlagError("Quadrature basis is not associated.",err,error,*999)    
    IF(SIZE(quadratureNumberOfGaussXi,1)<SIZE(basis%quadrature%numberOfGaussXi,1)) THEN
      localError="The size of quadratureNumberOfGaussXi is too small. The supplied size is "// &
        & TRIM(NumberToVString(SIZE(quadratureNumberOfGaussXi,1),"*",err,error))//" and it needs to be >= "// &
        & TRIM(NumberToVString(SIZE(basis%quadrature%numberOfGaussXi,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    quadratureNumberOfGaussXi=basis%quadrature%numberOfGaussXi
     
    EXITS("Basis_QuadratureNumberOfGaussXiGet")
    RETURN
999 ERRORSEXITS("Basis_QuadratureNumberOfGaussXiGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_QuadratureNumberOfGaussXiGet
    
  !
  !================================================================================================================================
  !

  !>Get the order of a quadrature for a basis quadrature identified by a pointer. \see OpenCMISS::Iron::cmfe_Basis_QuadratureOrderGet
  SUBROUTINE Basis_QuadratureOrderGet(basis,quadratureOrder,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: quadratureOrder !<On return, the quadrature order for the specified basis.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Basis_QuadratureOrderGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    CALL Basis_AssertIsFinished(basis,err,error,*999)
    IF(.NOT.ASSOCIATED(basis%quadrature%basis)) CALL FlagError("Quadrature basis is not associated.",err,error,*999)
#endif    
    
    quadratureOrder=basis%quadrature%gaussOrder
      
    EXITS("Basis_QuadratureOrderGet")
    RETURN
999 ERRORSEXITS("Basis_QuadratureOrderGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_QuadratureOrderGet

  !
  !================================================================================================================================
  !

  !>Returns a quadrature scheme for the basis. 
  SUBROUTINE Basis_QuadratureSchemeGet(basis,quadratureSchemeIdx,quadratureScheme,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the quadrature scheme for
    INTEGER(INTG), INTENT(IN) :: quadratureSchemeIdx !<The index of the quadrature scheme to get. \see Basis_QuadratureScheme,BasisRoutines,BasisAccessRoutines
    TYPE(QuadratureSchemeType), POINTER :: quadratureScheme !<On exit, the basis quadrature scheme corresponding to the quadrature index. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("Basis_QuadratureSchemeGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(quadratureScheme)) CALL FlagError("Quadrature scheme is already associated.",err,error,*998)
    CALL Basis_AssertIsFinished(basis,err,error,*999)
    IF(quadratureSchemeIdx<1.OR.quadratureSchemeIdx>BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES) THEN
      localError="The specified quadrature scheme index of "//TRIM(NumberToVString(quadratureSchemeIdx,"*",err,error))// &
        & " is invalid. The quadrature scheme index should be >= 1 and <= "// &
        & TRIM(NumberToVString(BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    END IF
    IF(.NOT.ALLOCATED(basis%quadrature%quadratureSchemeMap)) &
      & CALL FlagError("Basis quadrature scheme map has not been allocated.",err,error,*999)
#endif    

    quadratureScheme=>basis%quadrature%quadratureSchemeMap(quadratureSchemeIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(quadratureScheme)) THEN
      localError="The quadrature scheme for the quadrature scheme index of "// &
        & TRIM(NumberToVString(quadratureSchemeIdx,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("Basis_QuadratureSchemeGet")
    RETURN
999 NULLIFY(quadratureScheme)
998 ERRORSEXITS("Basis_QuadratureSchemeGet",err,error)    
    RETURN 1
    
  END SUBROUTINE Basis_QuadratureSchemeGet

  !
  !================================================================================================================================
  !
  
  !>get the quadrature type on a basis identified by a pointer. \see OpenCMISS::Iron::cmfe_Basis_QuadratureTypeGet
  SUBROUTINE Basis_QuadratureTypeGet(basis,quadratureType,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis
    INTEGER(INTG), INTENT(OUT) :: quadratureType !<On return, the quadrature type to be get \see BasisRoutines_QuadratureTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("Basis_QuadratureTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    CALL Basis_AssertIsFinished(basis,err,error,*999)
    IF(.NOT.ASSOCIATED(basis%quadrature%basis)) THEN
      localError="The basis quadrature basis is not associated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    quadratureType=basis%quadrature%type
   
    EXITS("Basis_QuadratureTypeGet")
    RETURN
999 ERRORSEXITS("Basis_QuadratureTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_QuadratureTypeGet

  !
  !================================================================================================================================
  !
  
  !>get the type for a basis is identified by a a pointer. \see OpenCMISS::Iron::cmfe_Basis_TypeGet
  SUBROUTINE Basis_TypeGet(basis,type,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get
    INTEGER(INTG), INTENT(OUT) :: type !<On return, the type of the specified basis. \see BasisRoutines_BasisTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Basis_TypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    CALL Basis_AssertIsFinished(basis,err,error,*999)
#endif    
    
    type=basis%type
    
    EXITS("Basis_TypeGet")
    RETURN
999 ERRORSEXITS("Basis_TypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_TypeGet

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to a basis with the given user number. If no basis with that number exists basis is left nullified.
  SUBROUTINE Basis_UserNumberFind(basisFunctions,userNumber,basis,err,error,*)

    !Argument variables
    TYPE(BasisFunctionsType), POINTER :: basisFunctions !<The basis functions to find the user number for.
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the basis to find.
    TYPE(BasisType), POINTER :: basis !<On exit, a pointer to the found basis. If no basis with the given user number exists the pointer is NULL.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Basis_UserNumberFind",err,error,*999)
    
    CALL Basis_FamilyNumberFind(basisFunctions,userNumber,0,basis,err,error,*999)
  
    EXITS("Basis_UserNumberFind")
    RETURN
999 ERRORSEXITS("Basis_UserNumberFind",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_UserNumberFind

  !
  !================================================================================================================================
  !

  !>Returns the user number for a basis.
  SUBROUTINE Basis_UserNumberGet(basis,userNumber,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the user number for
    INTEGER(INTG), INTENT(OUT) :: userNumber !<On exit, the user number of the basis.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Basis_UserNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
#endif    

    userNumber=basis%userNumber
  
    EXITS("Basis_UserNumberGet")
    RETURN
999 ERRORSEXITS("Basis_UserNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_UserNumberGet

  !
  !================================================================================================================================
  !
  
  !>Returns the Gauss basis function for Gauss point in a basis quadrature scheme.  
  SUBROUTINE BasisQuadratureScheme_GaussBasisFunctionGet(quadratureScheme,elementParameterIdx,partialDerivativeIdx, &
    & gaussIdx,gaussBasisFunction,err,error,*)

    !Argument variables
    TYPE(QuadratureSchemeType), POINTER :: quadratureScheme !<The basis quadrature scheme to get the Gauss position for.
    INTEGER(INTG), INTENT(IN) :: elementParameterIdx !<The element parameter index to get the Gauss basis function for.
    INTEGER(INTG), INTENT(IN) :: partialDerivativeIdx !<The partial derivative index to get the Gauss basis function for.
    INTEGER(INTG), INTENT(IN) :: gaussIdx !<The Gauss point index to get the Gauss position for.
    REAL(DP), INTENT(OUT) :: gaussBasisFunction !<On exit, the specified Gauss basis function for the Gauss point in the quadrature scheme.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("BasisQuadratureScheme_GaussBasisFunctionGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(quadratureScheme)) CALL FlagError("Quadrature scheme is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(quadratureScheme%gaussBasisFunctions)) &
      & CALL FLagError("Gauss basis functions are not allocated for the quadrature scheme.",err,error,*999)
    IF(elementParameterIdx<1.OR.elementParameterIdx>SIZE(quadratureScheme%gaussBasisFunctions,1)) THEN
      localError="The specified element parameter index of "//TRIM(NumberToVString(elementParameterIdx,"*",err,error))// &
        & " is invalid. The element parameter index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(quadratureScheme%gaussBasisFunctions,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(partialDerivativeIdx<1.OR.partialDerivativeIdx>SIZE(quadratureScheme%gaussBasisFunctions,2)) THEN
      localError="The specified partial derivative index of "//TRIM(NumberToVString(partialDerivativeIdx,"*",err,error))// &
        & " is invalid. The partial derivative index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(quadratureScheme%gaussBasisFunctions,2),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(gaussIdx<1.OR.gaussIdx>quadratureScheme%numberOfGauss) THEN
      localError="The specified Gauss point index of "//TRIM(NumberToVString(gaussIdx,"*",err,error))// &
        & " is invalid. The Gauss point index should be >= 1 and <= "// &
        & TRIM(NumberToVString(quadratureScheme%numberOfGauss,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    gaussBasisFunction=quadratureScheme%gaussBasisFunctions(elementParameterIdx,partialDerivativeIdx,gaussIdx)
    
    EXITS("BasisQuadratureScheme_GaussBasisFunctionGet")
    RETURN
999 ERRORSEXITS("BasisQuadratureScheme_GaussBasisFunctionGet",err,error)    
    RETURN 1
    
  END SUBROUTINE BasisQuadratureScheme_GaussBasisFunctionGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the Gauss position for Gauss point in a basis quadrature scheme.  
  SUBROUTINE BasisQuadratureScheme_GaussPositionGet(quadratureScheme,gaussIdx,gaussPosition,err,error,*)

    !Argument variables
    TYPE(QuadratureSchemeType), POINTER :: quadratureScheme !<The basis quadrature scheme to get the Gauss position for.
    INTEGER(INTG), INTENT(IN) :: gaussIdx !<The Gauss point index to get the Gauss position for.
    REAL(DP), INTENT(OUT) :: gaussPosition(:) !<gaussPosition(xiCoordinateIdx). On exit, the Gauss position for the Gauss point in the quadrature scheme.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("BasisQuadratureScheme_GaussPositionGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(quadratureScheme)) CALL FlagError("Quadrature scheme is not associated.",err,error,*999)
    IF(gaussIdx<1.OR.gaussIdx>quadratureScheme%numberOfGauss) THEN
      localError="The specified Gauss point index of "//TRIM(NumberToVString(gaussIdx,"*",err,error))// &
        & " is invalid. The Gauss point index should be >= 1 and <= "// &
        & TRIM(NumberToVString(quadratureScheme%numberOfGauss,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(quadratureScheme%gaussPositions)) &
      & CALL FLagError("Gauss positions are not allocated for the quadrature scheme.",err,error,*999)
    IF(SIZE(gaussPosition,1)<SIZE(quadratureScheme%gaussPositions,1)) THEN
      localError="The size of the specified Gauss position array of "// &
        & TRIM(NumberToVString(SIZE(gaussPosition,1),"*",err,error))//" is too small. The size should be >= "// &
        & TRIM(NumberToVString(SIZE(quadratureScheme%gaussPositions,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    gaussPosition(1:SIZE(quadratureScheme%gaussPositions,1))= &
      & quadratureScheme%gaussPositions(1:SIZE(quadratureScheme%gaussPositions,1),gaussIdx)
    
    EXITS("BasisQuadratureScheme_GaussPositionGet")
    RETURN
999 ERRORSEXITS("BasisQuadratureScheme_GaussPositionGet",err,error)    
    RETURN 1
    
  END SUBROUTINE BasisQuadratureScheme_GaussPositionGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the Gauss weight for Gauss point in a basis quadrature scheme.  
  SUBROUTINE BasisQuadratureScheme_GaussWeightGet(quadratureScheme,gaussIdx,gaussWeight,err,error,*)

    !Argument variables
    TYPE(QuadratureSchemeType), POINTER :: quadratureScheme !<The basis quadrature scheme to get the Gauss weight for.
    INTEGER(INTG), INTENT(IN) :: gaussIdx !<The Gauss point index to get the Gauss weight for.
    REAL(DP), INTENT(OUT) :: gaussWeight !<On exit, the Gauss wieght for the Gauss point in the quadrature scheme.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("BasisQuadratureScheme_GaussWeightGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(quadratureScheme)) CALL FlagError("Quadrature scheme is not associated.",err,error,*999)
    IF(gaussIdx<1.OR.gaussIdx>quadratureScheme%numberOfGauss) THEN
      localError="The specified Gauss point index of "//TRIM(NumberToVString(gaussIdx,"*",err,error))// &
        & " is invalid. The Gauss point index should be >= 1 and <= "// &
        & TRIM(NumberToVString(quadratureScheme%numberOfGauss,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(quadratureScheme%gaussWeights)) &
      & CALL FLagError("Gauss weights are not allocated for the quadrature scheme.",err,error,*999)
#endif    

    gaussWeight=quadratureScheme%gaussWeights(gaussIdx)
    
    EXITS("BasisQuadratureScheme_GaussWeightGet")
    RETURN
999 ERRORSEXITS("BasisQuadratureScheme_GaussWeightGet",err,error)    
    RETURN 1
    
  END SUBROUTINE BasisQuadratureScheme_GaussWeightGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of Gauss for a basis quadrature scheme.  
  SUBROUTINE BasisQuadratureScheme_NumberOfGaussGet(quadratureScheme,numberOfGauss,err,error,*)

    !Argument variables
    TYPE(QuadratureSchemeType), POINTER :: quadratureScheme !<The basis quadrature scheme to get the number of Gauss points for.
    INTEGER(INTG), INTENT(OUT) :: numberOfGauss !<On exit, the number of Gauss points in the quadrature scheme.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("BasisQuadratureScheme_NumberOfGaussGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(quadratureScheme)) CALL FlagError("Quadrature scheme is not associated.",err,error,*999)
#endif    

    numberOfGauss=quadratureScheme%numberOfGauss
    
    EXITS("BasisQuadratureScheme_NumberOfGaussGet")
    RETURN
999 ERRORSEXITS("BasisQuadratureScheme_NumberOfGaussGet",err,error)    
    RETURN 1
    
  END SUBROUTINE BasisQuadratureScheme_NumberOfGaussGet
  
  !
  !================================================================================================================================
  !

END MODULE BasisAccessRoutines
