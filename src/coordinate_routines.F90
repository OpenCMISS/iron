!> \file
!> \author Chris Bradley
!> \brief This module contains all coordinate transformation and support routines.
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
!> Contributor(s): Kumar Mithraratne
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

!> This module contains all coordinate transformation and support routines.
MODULE COORDINATE_ROUTINES

  USE BaseRoutines
  USE Constants
  USE CoordinateSystemAccessRoutines
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Kinds
  USE Maths
  USE Strings
  USE Types
  
#include "macros.h"

  IMPLICIT NONE

!!TODO: Should all of the get/set/create/destroy routines be accessed via pointers???

  PRIVATE

  !Module parameters

  !> \addtogroup COORINDATE_ROUTINES_CoordinateSystemTypes COORDINATE_ROUTINES::CoordinateSystemTypes
  !> \see COORDINATE_ROUTINES
  !> Coordinate system type parameters
  !>@{ 
  INTEGER(INTG), PARAMETER :: COORDINATE_RECTANGULAR_CARTESIAN_TYPE=1 !<Rectangular Cartesian coordinate system type \see COORDINATE_ROUTINES_CoordinateSystemTypes,COORDINATE_ROUTINES
  INTEGER(INTG), PARAMETER :: COORDINATE_CYLINDRICAL_POLAR_TYPE=2 !<Cylindrical polar coordinate system type \see COORDINATE_ROUTINES_CoordinateSystemTypes,COORDINATE_ROUTINES
  INTEGER(INTG), PARAMETER :: COORDINATE_SPHERICAL_POLAR_TYPE=3 !<Spherical polar coordinate system type \see COORDINATE_ROUTINES_CoordinateSystemTypes,COORDINATE_ROUTINES
  INTEGER(INTG), PARAMETER :: COORDINATE_PROLATE_SPHEROIDAL_TYPE=4 !<Prolate spheroidal coordinate system type \see COORDINATE_ROUTINES_CoordinateSystemTypes,COORDINATE_ROUTINES
  INTEGER(INTG), PARAMETER :: COORDINATE_OBLATE_SPHEROIDAL_TYPE=5 !<Oblate spheroidal coordinate system type \see COORDINATE_ROUTINES_CoordinateSystemTypes,COORDINATE_ROUTINES
  !>@}

  !> \addtogroup COORDINATE_ROUTINES_RadialInterpolations COORDINATE_ROUTINES::RadialInterpolations
  !> \see COORDINATE_ROUTINES
  !> \brief The type of radial interpolation for polar coordinate systems
  !>@{
  INTEGER(INTG), PARAMETER :: COORDINATE_NO_RADIAL_INTERPOLATION_TYPE=0 !<No radial interpolation \see COORDINATE_ROUTINES_RadialInterpolations,COORDINATE_ROUTINES
  INTEGER(INTG), PARAMETER :: COORDINATE_RADIAL_INTERPOLATION_TYPE=1 !<r radial interpolation \see COORDINATE_ROUTINES_RadialInterpolations,COORDINATE_ROUTINES
  INTEGER(INTG), PARAMETER :: COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE=2 !<r^2 radial interpolation \see COORDINATE_ROUTINES_RadialInterpolations,COORDINATE_ROUTINES
  INTEGER(INTG), PARAMETER :: COORDINATE_RADIAL_CUBED_INTERPOLATION_TYPE=3 !<r^3 radial interpolation \see COORDINATE_ROUTINES_RadialInterpolations,COORDINATE_ROUTINES
  !>@}
  
  !> \addtogroup COORDINATE_ROUTINES_JacobianType COORDINATE_ROUTINES::JacobianType
  !> \see COORDINATE_ROUTINES
  !> \brief The type of Jacobian to return when coordinate metrics are calculated.
  !>@{
  INTEGER(INTG), PARAMETER :: COORDINATE_JACOBIAN_NO_TYPE=0 !<No Jacobian \see COORDINATE_ROUTINES_JacobianTypes,COORDINATE_ROUTINES
  INTEGER(INTG), PARAMETER :: COORDINATE_JACOBIAN_LINE_TYPE=1 !<Line type Jacobian \see COORDINATE_ROUTINES_JacobianTypes,COORDINATE_ROUTINES
  INTEGER(INTG), PARAMETER :: COORDINATE_JACOBIAN_AREA_TYPE=2 !<Area type Jacobian \see COORDINATE_ROUTINES_JacobianTypes,COORDINATE_ROUTINES
  INTEGER(INTG), PARAMETER :: COORDINATE_JACOBIAN_VOLUME_TYPE=3 !<Volume type Jacobian \see COORDINATE_ROUTINES_JacobianTypes,COORDINATE_ROUTINES
  !>@}
  
  !Module types
  
  !Module variables

  CHARACTER(LEN=21) :: COORDINATE_SYSTEM_TYPE_STRING(5) = &
    & [ "Rectangular Cartesian",&
    &    "Cylindrical Polar    ", &
    &    "Spherical Polar      ", &
    &    "Prolate Spheroidal   ", &
    &    "Oblate Spheroidal    " ]

  !Interfaces

  !>COORDINATE_CONVERT_FROM_RC performs a coordinate transformation from a rectangular cartesian coordinate at the point with
  !>coordinate Z(:) to the returned point with coordinate X(:) in the coordinate system identified by COORDINATE_SYSTEM.
  INTERFACE COORDINATE_CONVERT_FROM_RC
    MODULE PROCEDURE COORDINATE_CONVERT_FROM_RC_DP
    MODULE PROCEDURE COORDINATE_CONVERT_FROM_RC_SP
  END INTERFACE !COORDINATE_CONVERT_FROM_RC

  !>COORDINATE_CONVERT_TO_RC performs a coordinate transformation from a coordinate system identified by COORDINATE_SYSTEM
  !>at the point X(:) to the returned point Z(:) in rectangular cartesian coordinates.
  INTERFACE COORDINATE_CONVERT_TO_RC
    MODULE PROCEDURE COORDINATE_CONVERT_TO_RC_DP
    MODULE PROCEDURE COORDINATE_CONVERT_TO_RC_SP
  END INTERFACE !COORDINATE_CONVERT_TO_RC

  !>Calculates the difference (or delta) between two points in a coordinate system. Discontinuities for polar coordinate
  !>systems are accounted for
  INTERFACE COORDINATE_DELTA_CALCULATE
    MODULE PROCEDURE COORDINATE_DELTA_CALCULATE_DP
    !MODULE PROCEDURE COORDINATE_DELTA_CALCULATE_SP
  END INTERFACE !COORDINATE_DELTA_CALCULATE

  !>Calculates DX(:)/DZ(I) at X, where Z(I) are rectangular cartesian and X(:) are curvilinear coordinates defined by COORDINATE_SYSTEM. \todo CHANGE NAME TO SOMETHING MORE MEANINGFULL?
  INTERFACE DXZ
    MODULE PROCEDURE DXZ_DP
    !MODULE PROCEDURE DXZ_SP
  END INTERFACE !Dxz

  !!TODO:: CHANGE NAME TO SOMETHING MORE MEANINGFULL?
  INTERFACE D2ZX
    MODULE PROCEDURE D2ZX_DP
    !MODULE PROCEDURE D2ZX_SP
  END INTERFACE !D2ZX

  !!TODO:: CHANGE NAME TO SOMETHING MORE MEANINGFULL?
  INTERFACE DZX
    MODULE PROCEDURE DZX_DP
    !MODULE PROCEDURE DZX_SP
  END INTERFACE !DZX

  INTERFACE COORDINATE_DERIVATIVE_CONVERT_TO_RC
    MODULE PROCEDURE COORDINATE_DERIVATIVE_CONVERT_TO_RC_DP
    MODULE PROCEDURE COORDINATE_DERIVATIVE_CONVERT_TO_RC_SP
  END INTERFACE !COORDINATE_DERIVATIVE_CONVERT_TO_RC

  PUBLIC COORDINATE_RECTANGULAR_CARTESIAN_TYPE,COORDINATE_CYLINDRICAL_POLAR_TYPE,COORDINATE_SPHERICAL_POLAR_TYPE, &
    & COORDINATE_PROLATE_SPHEROIDAL_TYPE,COORDINATE_OBLATE_SPHEROIDAL_TYPE

  PUBLIC COORDINATE_NO_RADIAL_INTERPOLATION_TYPE,COORDINATE_RADIAL_INTERPOLATION_TYPE, &
    & COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE,COORDINATE_RADIAL_CUBED_INTERPOLATION_TYPE
  
  PUBLIC COORDINATE_JACOBIAN_NO_TYPE,COORDINATE_JACOBIAN_LINE_TYPE,COORDINATE_JACOBIAN_AREA_TYPE,COORDINATE_JACOBIAN_VOLUME_TYPE
  
  PUBLIC COORDINATE_SYSTEM_TYPE_STRING

  PUBLIC COORDINATE_CONVERT_FROM_RC,COORDINATE_CONVERT_TO_RC,COORDINATE_DELTA_CALCULATE,COORDINATE_DERIVATIVE_NORM, &
    & COORDINATE_INTERPOLATION_ADJUST,COORDINATE_INTERPOLATION_PARAMETERS_ADJUST,COORDINATE_METRICS_CALCULATE
  
  PUBLIC COORDINATE_SYSTEM_DIMENSION_GET,COORDINATE_SYSTEM_DIMENSION_SET

  PUBLIC COORDINATE_SYSTEM_FOCUS_GET,COORDINATE_SYSTEM_FOCUS_SET

  PUBLIC Coordinates_RadialInterpolationTypeGet,Coordinates_RadialInterpolationTypeSet

  PUBLIC COORDINATE_SYSTEM_TYPE_GET,COORDINATE_SYSTEM_TYPE_SET

  PUBLIC COORDINATE_SYSTEM_ORIGIN_GET,COORDINATE_SYSTEM_ORIGIN_SET
  
  PUBLIC COORDINATE_SYSTEM_ORIENTATION_GET,COORDINATE_SYSTEM_ORIENTATION_SET

  PUBLIC COORDINATE_SYSTEM_CREATE_START,COORDINATE_SYSTEM_CREATE_FINISH

  PUBLIC COORDINATE_SYSTEM_DESTROY

  PUBLIC COORDINATE_DERIVATIVE_CONVERT_TO_RC
  
  PUBLIC Coordinates_MaterialSystemCalculate

  PUBLIC CoordinateSystems_Initialise,CoordinateSystems_Finalise

CONTAINS

  !
  !================================================================================================================================
  !

  !>Performs a coordinate transformation from a rectangular cartesian coordinate at the point with coordinate Z(:) to the returned point with coordinate X(:) in the coordinate system identified by COORDINATE_SYSTEM for double precision coordinates.
  FUNCTION COORDINATE_CONVERT_FROM_RC_DP(COORDINATE_SYSTEM,Z,err,error)
  
    !Argument variables
    TYPE(CoordinateSystemType), INTENT(IN) :: COORDINATE_SYSTEM !<The coordinate system to perform the conversion on
    REAL(DP), INTENT(IN) :: Z(:) !<The rectangular cartesian coordiantes to convert
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: COORDINATE_CONVERT_FROM_RC_DP(SIZE(Z,1))
    !Local variables
    REAL(DP) :: A1,A2,A3,A4,A5,A6,A7,A8,A9,FOCUS
    
    ENTERS("COORDINATE_CONVERT_FROM_RC_DP",err,error,*999)

    COORDINATE_CONVERT_FROM_RC_DP=0.0_DP

    IF(SIZE(Z,1)<COORDINATE_SYSTEM%numberOfDimensions) &
      & CALL FlagError("Size of Z is less than the number of dimensions.",err,error,*999)
    
    SELECT CASE(COORDINATE_SYSTEM%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      COORDINATE_CONVERT_FROM_RC_DP(1:COORDINATE_SYSTEM%numberOfDimensions)=Z(1:COORDINATE_SYSTEM%numberOfDimensions)
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      SELECT CASE(COORDINATE_SYSTEM%numberOfDimensions)
      CASE(2)
        COORDINATE_CONVERT_FROM_RC_DP(1)=SQRT(Z(1)**2+Z(2)**2)
        COORDINATE_CONVERT_FROM_RC_DP(2)=ATAN2(Z(1),Z(2))
      CASE(3)
        COORDINATE_CONVERT_FROM_RC_DP(1)=SQRT(Z(1)**2+Z(2)**2)
        COORDINATE_CONVERT_FROM_RC_DP(2)=ATAN2(Z(1),Z(2))
        COORDINATE_CONVERT_FROM_RC_DP(3)=Z(3)
      CASE DEFAULT
        CALL FlagError("Invalid number of coordinates.",err,error,*999)
      END SELECT
      IF(COORDINATE_CONVERT_FROM_RC_DP(2)<0.0_DP) &
        & COORDINATE_CONVERT_FROM_RC_DP(2)=COORDINATE_CONVERT_FROM_RC_DP(2)+2.0_DP*PI !reference coordinate 0->2*pi
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(COORDINATE_SYSTEM%numberOfDimensions==3) THEN
        COORDINATE_CONVERT_FROM_RC_DP(1)=SQRT(Z(1)**2+Z(2)**2+Z(3)**2)
        IF(ABS(Z(1))>=ZERO_TOLERANCE.OR.ABS(Z(2))>=ZERO_TOLERANCE) THEN
          COORDINATE_CONVERT_FROM_RC_DP(2)=ATAN2(Z(2),Z(1))
        ELSE
          COORDINATE_CONVERT_FROM_RC_DP(2)=0.0_DP
        ENDIF
        A1=SQRT(Z(1)**2+Z(2)**2)
        IF(ABS(Z(3))>=ZERO_TOLERANCE.OR.ABS(A1)>=ZERO_TOLERANCE) THEN
          COORDINATE_CONVERT_FROM_RC_DP(3)=ATAN2(Z(3),A1)
        ELSE
          COORDINATE_CONVERT_FROM_RC_DP(3)=0.0_DP
        ENDIF
      ELSE
        CALL FlagError("Invalid number of coordinates.",err,error,*999)
      ENDIF
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      FOCUS=COORDINATE_SYSTEM%FOCUS
      IF(COORDINATE_SYSTEM%numberOfDimensions==3) THEN
        A1=Z(1)**2+Z(2)**2+Z(3)**2-FOCUS**2
        A2=SQRT(A1**2+4.0_DP*(FOCUS**2)*(Z(2)**2+Z(3)**2))
        A3=2.0_DP*FOCUS**2
        A4=MAX((A2+A1)/A3,0.0_DP)
        A5=MAX((A2-A1)/A3,0.0_DP)
        A6=SQRT(A4)
        A7=MIN(SQRT(A5),1.0_DP)
        IF(ABS(A7)<=1.0_DP) THEN
          A8=ASIN(A7)
        ELSE
          A8=0.0_DP
          CALL FlagWarning("Put A8=0 since ABS(A8)>1.",err,error,*999)
        ENDIF
        IF((ABS(Z(3))<ZERO_TOLERANCE).OR.(ABS(A6)<ZERO_TOLERANCE).OR.(ABS(A7)<ZERO_TOLERANCE)) THEN
          A9=0.0_DP
        ELSE
          IF(ABS(A6*A7)>0.0_DP) THEN
            A9=Z(3)/(FOCUS*A6*A7)
          ELSE
            A9=0.0_DP
            CALL FlagWarning("Put A9=0 since A6*A7=0.",err,error,*999)
          ENDIF
          IF(A9>=1.0_DP) THEN
            A9=PI/2.0_DP
          ELSE IF(A9<=-1.0_DP) THEN
            A9=-PI/2.0_DP
          ELSE
            A9=ASIN(A9)
          ENDIF
        ENDIF
        COORDINATE_CONVERT_FROM_RC_DP(1)=LOG(A6+SQRT(A4+1.0_DP))
        IF(ABS(Z(1))>=ZERO_TOLERANCE) THEN
          COORDINATE_CONVERT_FROM_RC_DP(2)=A8
        ELSE
          COORDINATE_CONVERT_FROM_RC_DP(2)=PI-A8
        ENDIF
        IF(ABS(Z(2))>ZERO_TOLERANCE) THEN
          COORDINATE_CONVERT_FROM_RC_DP(3)=MOD(A9+2.0_DP*PI,2.0_DP*PI)
        ELSE
          COORDINATE_CONVERT_FROM_RC_DP(3)=PI-A9
        ENDIF
      ELSE
        CALL FlagError("Invalid number of coordinates.",err,error,*999)
      ENDIF
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      CALL FlagError("Invalid coordinate type.",err,error,*999)
    END SELECT

    EXITS("COORDINATE_CONVERT_FROM_RC_DP")
    RETURN
999 ERRORSEXITS("COORDINATE_CONVERT_FROM_RC_DP",err,error)
    RETURN 
  END FUNCTION COORDINATE_CONVERT_FROM_RC_DP
  
  !
  !================================================================================================================================
  !

  !>COORDINATE_CONVERT_FROM_RC_SP performs a coordinate transformation from a rectangular cartesian coordinate at the
  !>point with coordinate Z(:) to the returned point with coordinate X(:) in the coordinate system identified by
  !>COORDINATE_SYSTEM for single precision coordinates.
  FUNCTION COORDINATE_CONVERT_FROM_RC_SP(COORDINATE_SYSTEM,Z,err,error)
  
    !Argument variables
    TYPE(CoordinateSystemType), INTENT(IN) :: COORDINATE_SYSTEM !<The coordinate system to convert from RC to
    REAL(SP), INTENT(IN) :: Z(:) !<The coordinate to convert from rectangular cartesian
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(SP) :: COORDINATE_CONVERT_FROM_RC_SP(SIZE(Z,1))
    !Local variables
    REAL(SP) :: A1,A2,A3,A4,A5,A6,A7,A8,A9,FOCUS
    
    ENTERS("COORDINATE_CONVERT_FROM_RC_SP",err,error,*999)

    COORDINATE_CONVERT_FROM_RC_SP=0.0_SP
    
    IF(SIZE(Z,1)<COORDINATE_SYSTEM%numberOfDimensions) &
      & CALL FlagError("Size of Z is less than the number of dimensions.",err,error,*999)
    
    SELECT CASE(COORDINATE_SYSTEM%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      COORDINATE_CONVERT_FROM_RC_SP(1:COORDINATE_SYSTEM%numberOfDimensions)=Z(1:COORDINATE_SYSTEM%numberOfDimensions)
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      SELECT CASE(COORDINATE_SYSTEM%numberOfDimensions)
      CASE(2)
        COORDINATE_CONVERT_FROM_RC_SP(1)=SQRT(Z(1)**2+Z(2)**2)
        COORDINATE_CONVERT_FROM_RC_SP(2)=ATAN2(Z(1),Z(2))
      CASE(3)
        COORDINATE_CONVERT_FROM_RC_SP(1)=SQRT(Z(1)**2+Z(2)**2)
        COORDINATE_CONVERT_FROM_RC_SP(2)=ATAN2(Z(1),Z(2))
        COORDINATE_CONVERT_FROM_RC_SP(3)=Z(3)
      CASE DEFAULT
        CALL FlagError("Invalid number of coordinates.",err,error,*999)
      END SELECT
      IF(COORDINATE_CONVERT_FROM_RC_SP(2)<0.0_SP)  &
        & COORDINATE_CONVERT_FROM_RC_SP(2)=COORDINATE_CONVERT_FROM_RC_SP(2)+2.0_SP*REAL(PI,SP) !reference coordinate 0->2*pi
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(COORDINATE_SYSTEM%numberOfDimensions==3) THEN
        COORDINATE_CONVERT_FROM_RC_SP(1)=SQRT(Z(1)**2+Z(2)**2+Z(3)**2)
        IF(ABS(Z(1))>=ZERO_TOLERANCE_SP.OR.ABS(Z(2))>ZERO_TOLERANCE_SP) THEN
          COORDINATE_CONVERT_FROM_RC_SP(2)=ATAN2(Z(2),Z(1))
        ELSE
          COORDINATE_CONVERT_FROM_RC_SP(2)=0.0_SP
        ENDIF
        A1=SQRT(Z(1)**2+Z(2)**2)
        IF(ABS(Z(3))>=ZERO_TOLERANCE_SP.OR.ABS(A1)>ZERO_TOLERANCE_SP) THEN
          COORDINATE_CONVERT_FROM_RC_SP(3)=ATAN2(Z(3),A1)
        ELSE
          COORDINATE_CONVERT_FROM_RC_SP(3)=0.0_SP
        ENDIF
      ELSE
        CALL FlagError("Invalid number of coordinates.",err,error,*999)
      ENDIF
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      FOCUS=REAL(COORDINATE_SYSTEM%FOCUS,SP)
      IF(COORDINATE_SYSTEM%numberOfDimensions==3) THEN
        A1=Z(1)**2+Z(2)**2+Z(3)**2-FOCUS**2
        A2=SQRT(A1**2+4.0_SP*(FOCUS**2)*(Z(2)**2+Z(3)**2))
        A3=2.0_SP*FOCUS**2
        A4=MAX((A2+A1)/A3,0.0_SP)
        A5=MAX((A2-A1)/A3,0.0_SP)
        A6=SQRT(A4)
        A7=MIN(SQRT(A5),1.0_SP)
        IF(ABS(A7)<=1.0_SP) THEN
          A8=ASIN(A7)
        ELSE
          A8=0.0_SP
          CALL FlagWarning("Put A8=0 since ABS(A8)>1.",err,error,*999)
        ENDIF
        IF((ABS(Z(3))<ZERO_TOLERANCE_SP).OR.(ABS(A6)<ZERO_TOLERANCE_SP).OR.(ABS(A7)<ZERO_TOLERANCE_SP)) THEN
          A9=0.0_SP
        ELSE
          IF(ABS(A6*A7)>ZERO_TOLERANCE_SP) THEN
            A9=Z(3)/(FOCUS*A6*A7)
          ELSE
            A9=0.0_SP
            CALL FlagWarning("Put A9=0 since A6*A7=0.",err,error,*999)
          ENDIF
          IF(A9>=1.0_SP) THEN
            A9=REAL(PI,SP)/2.0_SP
          ELSE IF(A9<=-1.0_SP) THEN
            A9=-REAL(PI,SP)/2.0_SP
          ELSE
            A9=ASIN(A9)
          ENDIF
        ENDIF
        COORDINATE_CONVERT_FROM_RC_SP(1)=LOG(A6+SQRT(A4+1.0_SP))
        IF(ABS(Z(1))>=ZERO_TOLERANCE_SP) THEN
          COORDINATE_CONVERT_FROM_RC_SP(2)=A8
        ELSE
          COORDINATE_CONVERT_FROM_RC_SP(2)=REAL(PI,SP)-A8
        ENDIF
        IF(ABS(Z(2))>=ZERO_TOLERANCE_SP) THEN
          COORDINATE_CONVERT_FROM_RC_SP(3)=MOD(A9+2.0_SP*REAL(PI,SP),2.0_SP*&
            & REAL(PI,SP))
        ELSE
          COORDINATE_CONVERT_FROM_RC_SP(3)=REAL(PI,SP)-A9
        ENDIF
      ELSE
        CALL FlagError("Invalid number of coordinates.",err,error,*999)
      ENDIF
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      CALL FlagError("Invalid coordinate type.",err,error,*999)
    END SELECT

    EXITS("COORDINATE_CONVERT_FROM_RC_SP")
    RETURN
999 ERRORSEXITS("COORDINATE_CONVERT_FROM_RC_SP",err,error)
    RETURN 
  END FUNCTION COORDINATE_CONVERT_FROM_RC_SP

  !
  !================================================================================================================================
  !

  !>COORDINATE_CONVERT_TO_RC_DP performs a coordinate transformation from a coordinate system identified by
  !>COORDINATE_SYSTEM at the point X(:) to the returned point Z(:) in rectangular cartesian coordinates for
  !>double precision coordinates.
  FUNCTION COORDINATE_CONVERT_TO_RC_DP(COORDINATE_SYSTEM,X,err,error)
  
    !Argument variables
    TYPE(CoordinateSystemType), INTENT(IN) :: COORDINATE_SYSTEM !<The coordinate system to convert to rectangular cartesian
    REAL(DP), INTENT(IN) :: X(:) !<The coordiante to convert
    INTEGER(INTG), INTENT(OUT) :: err !<The error coode
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: COORDINATE_CONVERT_TO_RC_DP(SIZE(X,1))
    !Local variables
    REAL(DP) :: FOCUS
    
    ENTERS("COORDINATE_CONVERT_TO_RC_DP",err,error,*999)
    
    COORDINATE_CONVERT_TO_RC_DP=0.0_DP

    IF(SIZE(X,1)<COORDINATE_SYSTEM%numberOfDimensions) &
      & CALL FlagError("Size of X is less than the number of dimensions.",err,error,*999)

    SELECT CASE(COORDINATE_SYSTEM%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      COORDINATE_CONVERT_TO_RC_DP(1:COORDINATE_SYSTEM%numberOfDimensions)=X(1:COORDINATE_SYSTEM%numberOfDimensions)
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      SELECT CASE(COORDINATE_SYSTEM%numberOfDimensions)
      CASE(2)
        COORDINATE_CONVERT_TO_RC_DP(1)=X(1)*COS(X(2))
        COORDINATE_CONVERT_TO_RC_DP(2)=X(1)*SIN(X(2))
      CASE(3)
        COORDINATE_CONVERT_TO_RC_DP(1)=X(1)*COS(X(2))
        COORDINATE_CONVERT_TO_RC_DP(2)=X(1)*SIN(X(2))
        COORDINATE_CONVERT_TO_RC_DP(3)=X(3)
      CASE DEFAULT
        CALL FlagError("Invalid number of coordinates.",err,error,*999)
      END SELECT
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(COORDINATE_SYSTEM%numberOfDimensions==3) THEN  
        COORDINATE_CONVERT_TO_RC_DP(1)=X(1)*COS(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_DP(2)=X(1)*SIN(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_DP(3)=X(1)*SIN(X(3))
      ELSE
        CALL FlagError("Invalid number of coordinates.",err,error,*999)
      ENDIF
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      IF(COORDINATE_SYSTEM%numberOfDimensions==3) THEN
        FOCUS=COORDINATE_SYSTEM%FOCUS
        COORDINATE_CONVERT_TO_RC_DP(1)=FOCUS*COSH(X(1))*COS(X(2))
        COORDINATE_CONVERT_TO_RC_DP(2)=FOCUS*SINH(X(1))*SIN(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_DP(3)=FOCUS*SINH(X(1))*SIN(X(2))*SIN(X(3))
      ELSE
        CALL FlagError("Invalid number of coordinates.",err,error,*999)
      ENDIF
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      IF(COORDINATE_SYSTEM%numberOfDimensions==3) THEN
        FOCUS=COORDINATE_SYSTEM%FOCUS
        COORDINATE_CONVERT_TO_RC_DP(1)=FOCUS*COSH(X(1))*COS(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_DP(2)=FOCUS*SINH(X(1))*SIN(X(2))
        COORDINATE_CONVERT_TO_RC_DP(3)=FOCUS*COSH(X(1))*COS(X(2))*SIN(X(3))
      ELSE
        CALL FlagError("Invalid number of coordinates.",err,error,*999)
      ENDIF
    CASE DEFAULT
      CALL FlagError("Invalid coordinate type.",err,error,*999)
    END SELECT

    EXITS("COORDINATE_CONVERT_TO_RC_DP")
    RETURN
999 ERRORSEXITS("COORDINATE_CONVERT_TO_RC_DP",err,error)
    RETURN 
  END FUNCTION COORDINATE_CONVERT_TO_RC_DP

  !
  !================================================================================================================================
  !

  !>COORDINATE_CONVERT_TO_RC_SP performs a coordinate transformation from a coordinate system identified by
  !>COORDINATE_SYSTEM at the point X(:) to the returned point Z(:) in rectangular cartesian coordinates for
  !>single precision coordinates.
  FUNCTION COORDINATE_CONVERT_TO_RC_SP(COORDINATE_SYSTEM,X,err,error)
  
    !Argument variables
    TYPE(CoordinateSystemType), INTENT(IN) :: COORDINATE_SYSTEM !<The coordinate system to convert to rectangular cartesian
    REAL(SP), INTENT(IN) :: X(:) !<The coordinate to convert
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(SP) :: COORDINATE_CONVERT_TO_RC_SP(SIZE(X,1))
    !Local variables
    REAL(SP) :: FOCUS
    
    ENTERS("COORDINATE_CONVERT_TO_RC_SP",err,error,*999)
    
    COORDINATE_CONVERT_TO_RC_SP=0.0_SP

    IF(SIZE(X,1)<COORDINATE_SYSTEM%numberOfDimensions) &
      & CALL FlagError("Size of X is less than the number of dimensions.",err,error,*999)

    SELECT CASE(COORDINATE_SYSTEM%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      COORDINATE_CONVERT_TO_RC_SP(1:COORDINATE_SYSTEM%numberOfDimensions)=X(1:COORDINATE_SYSTEM%numberOfDimensions)
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      SELECT CASE(COORDINATE_SYSTEM%numberOfDimensions)
      CASE(2)
        COORDINATE_CONVERT_TO_RC_SP(1)=X(1)*COS(X(2))
        COORDINATE_CONVERT_TO_RC_SP(2)=X(1)*SIN(X(2))
      CASE(3)
        COORDINATE_CONVERT_TO_RC_SP(1)=X(1)*COS(X(2))
        COORDINATE_CONVERT_TO_RC_SP(2)=X(1)*SIN(X(2))
        COORDINATE_CONVERT_TO_RC_SP(3)=X(3)
      CASE DEFAULT
        CALL FlagError("Invalid number of coordinates.",err,error,*999)
      END SELECT
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(COORDINATE_SYSTEM%numberOfDimensions==3) THEN  
        COORDINATE_CONVERT_TO_RC_SP(1)=X(1)*COS(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_SP(2)=X(1)*SIN(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_SP(3)=X(1)*SIN(X(3))
      ELSE
        CALL FlagError("Invalid number of coordinates.",err,error,*999)
      ENDIF
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      IF(COORDINATE_SYSTEM%numberOfDimensions==3) THEN
        FOCUS=REAL(COORDINATE_SYSTEM%FOCUS,SP)
        COORDINATE_CONVERT_TO_RC_SP(1)=FOCUS*COSH(X(1))*COS(X(2))
        COORDINATE_CONVERT_TO_RC_SP(2)=FOCUS*SINH(X(1))*SIN(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_SP(3)=FOCUS*SINH(X(1))*SIN(X(2))*SIN(X(3))
      ELSE
        CALL FlagError("Invalid number of coordinates.",err,error,*999)
      ENDIF
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      IF(COORDINATE_SYSTEM%numberOfDimensions==3) THEN
        FOCUS=REAL(COORDINATE_SYSTEM%FOCUS,SP)
        COORDINATE_CONVERT_TO_RC_SP(1)=FOCUS*COSH(X(1))*COS(X(2))*COS(X(3))
        COORDINATE_CONVERT_TO_RC_SP(2)=FOCUS*SINH(X(1))*SIN(X(2))
        COORDINATE_CONVERT_TO_RC_SP(3)=FOCUS*COSH(X(1))*COS(X(2))*SIN(X(3))
      ELSE
        CALL FlagError("Invalid number of coordinates.",err,error,*999)
      ENDIF
    CASE DEFAULT
      CALL FlagError("Invalid coordinate type.",err,error,*999)
    END SELECT

    EXITS("COORDINATE_CONVERT_TO_RC_SP")
    RETURN
999 ERRORSEXITS("COORDINATE_CONVERT_TO_RC_SP",err,error)
    RETURN 
  END FUNCTION COORDINATE_CONVERT_TO_RC_SP

  !
  !================================================================================================================================
  !

  !>Calculates the difference (or detlta) between the point X and the point Y i.e., Y-X, in the given coordinate system.
  !>0->2Pi discontinuities with polar coordinates are accounted for.
  FUNCTION COORDINATE_DELTA_CALCULATE_DP(COORDINATE_SYSTEM,X,Y,err,error)
  
    !Argument variables
    TYPE(CoordinateSystemType), INTENT(IN) :: COORDINATE_SYSTEM !<The coordinate system to calculate the delta for
    REAL(DP), INTENT(IN) :: X(:) !<The first coordinate
    REAL(DP), INTENT(IN) :: Y(:) !<The second coordinate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: COORDINATE_DELTA_CALCULATE_DP(SIZE(X,1))
    !Local variables

    ENTERS("COORDINATE_DELTA_CALCULATE_DP",err,error,*999)

    COORDINATE_DELTA_CALCULATE_DP=0.0_DP

    IF(SIZE(X,1)<COORDINATE_SYSTEM%numberOfDimensions) &
      & CALL FlagError("Size of X is less than the number of dimensions.",err,error,*999)

    IF(SIZE(X,1)/=SIZE(Y,1)) &
      & CALL FlagError("Size of X is different to the size of Y.",err,error,*999)
   
    COORDINATE_DELTA_CALCULATE_DP(1:COORDINATE_SYSTEM%numberOfDimensions)=Y(1:COORDINATE_SYSTEM%numberOfDimensions)- &
      & X(1:COORDINATE_SYSTEM%numberOfDimensions)
    SELECT CASE(COORDINATE_SYSTEM%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      !Do nothing
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      CALL FlagError("Invalid coordinate type.",err,error,*999)
    END SELECT

    EXITS("COORDINATE_DELTA_CALCULATE_DP")
    RETURN
999 ERRORSEXITS("COORDINATE_DELTA_CALCULATE_DP",err,error)
    RETURN 
  END FUNCTION COORDINATE_DELTA_CALCULATE_DP

  !
  !================================================================================================================================
  !

  !>Calculates the covariant metric tensor GL(i,j), the contravariant metric tensor GU(i,J), the Jacobian and derivative of the interpolated coordinate system (XI_i) with respect to the given coordinate (X_j) system (DXI_DX) at a point (X - normally a Gauss point). Old cmiss name: XGMG
  SUBROUTINE COORDINATE_METRICS_CALCULATE(COORDINATE_SYSTEM,jacobianType,METRICS,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to calculate the metrics for
    INTEGER(INTG), INTENT(IN) :: jacobianType !<The type of Jacobian to calculate \see COORDINATE_ROUTINES_JacobianTypes,COORDINATE_ROUTINES
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: METRICS !<A pointer to the metrics to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: mi,ni,nu
    REAL(DP) :: DET_GL,DET_DX_DXI,DX_DXI2(3),DX_DXI3(3),FF,G1,G3,LENGTH,MU,R,RC,RCRC,RR,SCALE
    TYPE(FieldInterpolatedPointType), POINTER :: INTERPOLATED_POINT
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    NULLIFY(INTERPOLATED_POINT)

    ENTERS("COORDINATE_METRICS_CALCULATE",err,error,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(ASSOCIATED(METRICS)) THEN
        INTERPOLATED_POINT=>METRICS%interpolatedPoint
        IF(ASSOCIATED(INTERPOLATED_POINT)) THEN
          IF(INTERPOLATED_POINT%partialDerivativeType>=FIRST_PART_DERIV) THEN

            SELECT CASE(METRICS%numberOfXiDimensions)
            CASE(1)
              !Calculate the derivatives of X with respect to XI
              nu=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(1)
              METRICS%dXdXi(1:METRICS%numberOfXDimensions,1)=INTERPOLATED_POINT%VALUES(1:METRICS%numberOfXDimensions,nu)
              !Initialise the covariant metric tensor to the identity matrix
              METRICS%GL(1,1)=1.0_DP
            CASE(2)
              !Calculate the derivatives of X with respect to XI
              nu=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(1)
              METRICS%dXdXi(1:METRICS%numberOfXDimensions,1)=INTERPOLATED_POINT%VALUES(1:METRICS%numberOfXDimensions,nu)
              nu=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(2)
              METRICS%dXdXi(1:METRICS%numberOfXDimensions,2)=INTERPOLATED_POINT%VALUES(1:METRICS%numberOfXDimensions,nu)
              !Initialise the covariant metric tensor to the identity matrix
              METRICS%GL(1,1)=1.0_DP
              METRICS%GL(1,2)=0.0_DP
              METRICS%GL(2,1)=0.0_DP
              METRICS%GL(2,2)=1.0_DP
            CASE(3)
              !Calculate the derivatives of X with respect to XI
              nu=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(1)
              METRICS%dXdXi(1:METRICS%numberOfXDimensions,1)=INTERPOLATED_POINT%VALUES(1:METRICS%numberOfXDimensions,nu)
              nu=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(2)
              METRICS%dXdXi(1:METRICS%numberOfXDimensions,2)=INTERPOLATED_POINT%VALUES(1:METRICS%numberOfXDimensions,nu)
              nu=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(3)
              METRICS%dXdXi(1:METRICS%numberOfXDimensions,3)=INTERPOLATED_POINT%VALUES(1:METRICS%numberOfXDimensions,nu)
              !Initialise the covariant metric tensor to the identity matrix
              METRICS%GL(1,1)=1.0_DP
              METRICS%GL(1,2)=0.0_DP
              METRICS%GL(1,3)=0.0_DP
              METRICS%GL(2,1)=0.0_DP
              METRICS%GL(2,2)=1.0_DP
              METRICS%GL(2,3)=0.0_DP
              METRICS%GL(3,1)=0.0_DP
              METRICS%GL(3,2)=0.0_DP
              METRICS%GL(3,3)=1.0_DP
            CASE DEFAULT
              CALL FlagError("Not implemented.",err,error,*999)
            END SELECT
                        
            !Calculate the covariant metric tensor GL(i,j)
            SELECT CASE(COORDINATE_SYSTEM%TYPE)
            CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
              SELECT CASE(METRICS%numberOfXDimensions)
              CASE(1)
                DO mi=1,METRICS%numberOfXiDimensions
                  DO ni=1,METRICS%numberOfXiDimensions
                    METRICS%GL(mi,ni)=METRICS%dXdXi(1,mi)*METRICS%dXdXi(1,ni)                    
                  ENDDO !ni
                ENDDO !mi
              CASE(2)
                DO mi=1,METRICS%numberOfXiDimensions
                  DO ni=1,METRICS%numberOfXiDimensions
                    METRICS%GL(mi,ni)=METRICS%dXdXi(1,mi)*METRICS%dXdXi(1,ni)+ &
                      & METRICS%dXdXi(2,mi)*METRICS%dXdXi(2,ni)
                  ENDDO !ni
                ENDDO !mi
              CASE(3)
                DO mi=1,METRICS%numberOfXiDimensions
                  DO ni=1,METRICS%numberOfXiDimensions
                    METRICS%GL(mi,ni)=METRICS%dXdXi(1,mi)*METRICS%dXdXi(1,ni)+ &
                      & METRICS%dXdXi(2,mi)*METRICS%dXdXi(2,ni)+ &
                      & METRICS%dXdXi(3,mi)*METRICS%dXdXi(3,ni)
                  ENDDO !ni
                ENDDO !mi
              CASE DEFAULT
                LOCAL_ERROR=TRIM(NUMBER_TO_VSTRING(METRICS%numberOfXDimensions,"*",err,error))// &
                  & " is an invalid number of dimensions for a rectangular cartesian coordinate system."
                CALL FlagError(LOCAL_ERROR,err,error,*999)
              END SELECT
            CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
              R=INTERPOLATED_POINT%VALUES(1,1)
              RR=R*R
              IF(METRICS%numberOfXDimensions==2) THEN
                DO mi=1,METRICS%numberOfXiDimensions
                  DO ni=1,METRICS%numberOfXiDimensions
                    METRICS%GL(mi,ni)=METRICS%dXdXi(1,mi)*METRICS%dXdXi(1,ni)+RR*METRICS%dXdXi(2,mi)*METRICS%dXdXi(2,ni)
                  ENDDO !ni
                ENDDO !mi
              ELSE IF(METRICS%numberOfXDimensions==3) THEN
                DO mi=1,METRICS%numberOfXiDimensions
                  DO ni=1,METRICS%numberOfXiDimensions
                    METRICS%GL(mi,ni)=METRICS%dXdXi(1,mi)*METRICS%dXdXi(1,ni)+RR*METRICS%dXdXi(2,mi)*METRICS%dXdXi(2,ni)+ &
                      & METRICS%dXdXi(3,mi)*METRICS%dXdXi(3,ni)
                  ENDDO !ni
                ENDDO !mi
              ELSE
                LOCAL_ERROR=TRIM(NUMBER_TO_VSTRING(METRICS%numberOfXDimensions,"*",err,error))// &
                  & " is an invalid number of dimensions for a cylindrical polar coordinate system."
                CALL FlagError(LOCAL_ERROR,err,error,*999)
              ENDIF
            CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
              R=INTERPOLATED_POINT%VALUES(1,1)
              RR=R*R
              RC=R*COS(INTERPOLATED_POINT%VALUES(3,1))
              RCRC=RC*RC          
              DO mi=1,METRICS%numberOfXiDimensions
                DO ni=1,METRICS%numberOfXiDimensions
                  METRICS%GL(mi,ni)=METRICS%dXdXi(1,mi)*METRICS%dXdXi(1,ni)+RCRC*METRICS%dXdXi(2,mi)*METRICS%dXdXi(2,ni)+ &
                    & RR*METRICS%dXdXi(3,mi)*METRICS%dXdXi(3,ni)
                ENDDO !ni
              ENDDO !mi
            CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
              IF(ABS(INTERPOLATED_POINT%VALUES(2,1))<ZERO_TOLERANCE) THEN
                CALL FlagWarning("Mu is zero.",err,error,*999)
              ELSE
                FF=COORDINATE_SYSTEM%FOCUS*COORDINATE_SYSTEM%FOCUS
                R=INTERPOLATED_POINT%VALUES(1,1)
                MU=INTERPOLATED_POINT%VALUES(2,1)
                G1=FF*(SINH(R)*SINH(R)+SIN(MU)*SIN(MU))
                G3=FF*SINH(R)*SINH(R)*SIN(MU)*SIN(MU)
                IF(METRICS%numberOfXDimensions==2) THEN
                  DO mi=1,METRICS%numberOfXiDimensions
                    DO ni=1,METRICS%numberOfXiDimensions
                      METRICS%GL(mi,ni)=G1*(METRICS%dXdXi(1,mi)*METRICS%dXdXi(1,ni)+METRICS%dXdXi(2,mi)*METRICS%dXdXi(2,ni))
                    ENDDO !ni
                  ENDDO !mi
                ELSE IF(METRICS%numberOfXDimensions==3) THEN
                  DO mi=1,METRICS%numberOfXiDimensions
                    DO ni=1,METRICS%numberOfXiDimensions
                      METRICS%GL(mi,ni)=G1*(METRICS%dXdXi(1,mi)*METRICS%dXdXi(1,ni)+METRICS%dXdXi(2,mi)*METRICS%dXdXi(2,ni))+ &
                        & G3*METRICS%dXdXi(3,mi)*METRICS%dXdXi(3,ni)
                    ENDDO !ni
                  ENDDO !mi
                ELSE
                  LOCAL_ERROR=TRIM(NUMBER_TO_VSTRING(METRICS%numberOfXDimensions,"*",err,error))// &
                    & " is an invalid number of dimensions for a prolate spheroidal coordinate system."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                ENDIF
              ENDIF
            CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              LOCAL_ERROR="The coordinate system type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM%TYPE,"*",err,error))// &
                & " is invalid."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
            
            !Calcualte the contravariant metric tensor
            CALL INVERT(METRICS%GL(1:METRICS%numberOfXiDimensions,1:METRICS%numberOfXiDimensions), &
              & METRICS%GU(1:METRICS%numberOfXiDimensions,1:METRICS%numberOfXiDimensions),DET_GL, &
              & err,error,*999)

            !Calculate the Jacobian
            SELECT CASE(jacobianType)
            CASE(COORDINATE_JACOBIAN_NO_TYPE)
              METRICS%jacobian=0.0
              METRICS%jacobianType=COORDINATE_JACOBIAN_NO_TYPE
            CASE(COORDINATE_JACOBIAN_LINE_TYPE)
              METRICS%jacobian=SQRT(ABS(METRICS%GL(1,1)))
              METRICS%jacobianType=COORDINATE_JACOBIAN_LINE_TYPE
            CASE(COORDINATE_JACOBIAN_AREA_TYPE)
              IF(METRICS%numberOfXiDimensions==3) THEN
                METRICS%jacobian=SQRT(ABS(DET_GL*METRICS%GU(3,3)))
              ELSE
                METRICS%jacobian=SQRT(ABS(DET_GL))
              ENDIF
              METRICS%jacobianType=COORDINATE_JACOBIAN_AREA_TYPE
            CASE(COORDINATE_JACOBIAN_VOLUME_TYPE)
              METRICS%jacobian=SQRT(ABS(DET_GL))
              METRICS%jacobianType=COORDINATE_JACOBIAN_VOLUME_TYPE
            CASE DEFAULT
              LOCAL_ERROR="The Jacobian type of "//TRIM(NUMBER_TO_VSTRING(jacobianType,"*",err,error))// &
                & " is invalid."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
            
            !Calculate the derivatives of Xi with respect to X - DXI_DX
            IF(METRICS%numberOfXiDimensions==METRICS%numberOfXDimensions) THEN
              CALL INVERT(METRICS%dXdXi,METRICS%dXidX,DET_DX_DXI,err,error,*999)
            ELSE
              !We have a line or a surface embedded in a higher dimensional space
              SELECT CASE(METRICS%numberOfXiDimensions)
              CASE(1)
                !Line in space
                SELECT CASE(METRICS%numberOfXDimensions)
                CASE(2)
                  IF(INTERPOLATED_POINT%partialDerivativeType>FIRST_PART_DERIV) THEN
                    !We have curvature information. Form the frenet vector frame.
                    !Calculate the normal vector from the normalised second derivative of the position vector.
                    nu=PARTIAL_DERIVATIVE_SECOND_DERIVATIVE_MAP(1)
                    CALL Normalise(INTERPOLATED_POINT%VALUES(1:2,nu),DX_DXI2,err,error,*999)
                  ELSE
                    !No curvature information but obtain other normal frenet vector by rotating tangent vector 90 deg.
                    DX_DXI2(1)=-1.0_DP*METRICS%dXdXi(2,1)
                    DX_DXI2(2)=METRICS%dXdXi(1,1)                    
                  ENDIF
                  DET_DX_DXI=METRICS%dXdXi(1,1)*DX_DXI2(2)-METRICS%dXdXi(2,1)*DX_DXI2(1)
                  IF(ABS(DET_DX_DXI)>ZERO_TOLERANCE) THEN
                    METRICS%dXidX(1,1)=DX_DXI2(2)/DET_DX_DXI
                    METRICS%dXidX(1,2)=-1.0_DP*DX_DXI2(1)/DET_DX_DXI
                    !Normalise to ensure that g^11=g^1.g^1
                    CALL L2Norm(METRICS%dXidX(1,1:2),LENGTH,err,error,*999)
                    SCALE=SQRT(ABS(METRICS%GU(1,1)))/LENGTH
                    METRICS%dXidX(1,1:2)=SCALE*METRICS%dXidX(1,1:2)
                  ELSE
                    CALL FlagWarning("Zero determinant. Unable to obtain dxi/dx.",err,error,*999)
                    METRICS%dXidX=0.0_DP                    
                  ENDIF
                CASE(3)
                  IF(INTERPOLATED_POINT%partialDerivativeType>FIRST_PART_DERIV) THEN
                    !We have curvature information. Form the frenet vector frame.
                    !Calculate the normal vector from the normalised second derivative of the position vector.
                    nu=PARTIAL_DERIVATIVE_SECOND_DERIVATIVE_MAP(1)
                    CALL Normalise(INTERPOLATED_POINT%VALUES(1:3,nu),DX_DXI2,err,error,*999)
                    !Calculate the bi-normal vector from the normalised cross product of the tangent and normal vectors
                    CALL NormaliseCrossProduct(METRICS%dXdXi(1:3,1),DX_DXI2,DX_DXI3,err,error,*999)
                    DET_DX_DXI=METRICS%dXdXi(1,1)*(DX_DXI2(2)*DX_DXI3(3)-DX_DXI2(3)*DX_DXI3(2))+ &
                      & DX_DXI2(1)*(METRICS%dXdXi(3,1)*DX_DXI3(2)-DX_DXI3(3)*METRICS%dXdXi(2,1))+ &
                      & DX_DXI3(1)*(METRICS%dXdXi(2,1)*DX_DXI2(3)-METRICS%dXdXi(3,1)*DX_DXI2(2))
                    IF(ABS(DET_DX_DXI)>ZERO_TOLERANCE) THEN
                      METRICS%dXidX(1,1)=(DX_DXI3(3)*DX_DXI2(2)-DX_DXI2(3)*DX_DXI3(2))/DET_DX_DXI
                      METRICS%dXidX(1,2)=-1.0_DP*(DX_DXI3(3)*DX_DXI2(1)-DX_DXI2(3)*DX_DXI3(1))/DET_DX_DXI
                      METRICS%dXidX(1,3)=(DX_DXI3(2)*DX_DXI2(1)-DX_DXI2(2)*DX_DXI3(1))/DET_DX_DXI
                      !Normalise to ensure that g^11=g^1.g^1
                      CALL L2Norm(METRICS%dXidX(1,1:3),LENGTH,err,error,*999)
                      SCALE=SQRT(ABS(METRICS%GU(1,1)))/LENGTH
                      METRICS%dXidX(1,1:3)=SCALE*METRICS%dXdXi(1,1:3)
                    ELSE
                      CALL FlagWarning("Zero determinant. Unable to obtain dxi/dx.",err,error,*999)
                      METRICS%dXidX=0.0_DP                    
                    ENDIF
                  ELSE
                    METRICS%dXidX(1,1)=METRICS%dXdXi(1,1)
                    METRICS%dXidX(1,2)=METRICS%dXdXi(2,1)
                    METRICS%dXidX(1,3)=METRICS%dXdXi(3,1)
                    !Normalise to ensure that g^11=g^1.g^1
                    CALL L2Norm(METRICS%dXidX(1,1:3),LENGTH,err,error,*999)
                    SCALE=SQRT(ABS(METRICS%GU(1,1)))/LENGTH
                    METRICS%dXidX(1,1:3)=SCALE*METRICS%dXidX(1,1:3)
                  ENDIF
                CASE DEFAULT
                  CALL FlagError("Invalid embedding of a line in space.",err,error,*999)
                END SELECT
              CASE(2)
                !Surface in space
                IF(METRICS%numberOfXDimensions==3) THEN
                  !Surface in 3D space.
                  !Calculate the covariant vectors g^1 and g^2. These are calculated as follows:
                  !First define g_3=g_1 x g_2, and then define g^1=((g_2 x g_3)_b)/DET_GL and g^2=((g_3 x g_1)_b)/DET_GL. 
                  !The _b means lowering the index with the metric tensor of the curvilinear coordinate system.
                  !This way we have a consistent set of covariant and covariant vectors, i.e.  <g_M,g^N>=delta_M^N.
                  METRICS%dXidX(1,1:3)=(METRICS%GL(2,2)*METRICS%dXdXi(1:3,1)-METRICS%GL(1,2)*METRICS%dXdXi(1:3,2))/DET_GL
                  METRICS%dXidX(2,1:3)=(METRICS%GL(1,1)*METRICS%dXdXi(1:3,2)-METRICS%GL(2,1)*METRICS%dXdXi(1:3,1))/DET_GL
                  SELECT CASE(COORDINATE_SYSTEM%TYPE)
                  CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
                    !Do nothing
                  CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
                    R=INTERPOLATED_POINT%VALUES(1,1)
                    RR=R*R
                    METRICS%dXidX(1:2,2)=METRICS%dXidX(1:2,2)*RR
                  CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
                    R=INTERPOLATED_POINT%VALUES(1,1)
                    RR=R*R
                    RC=R*COS(INTERPOLATED_POINT%VALUES(3,1))
                    RCRC=RC*RC          
                    METRICS%dXidX(1:2,2)=METRICS%dXidX(1:2,2)*RCRC
                    METRICS%dXidX(1:2,3)=METRICS%dXidX(1:2,3)*RR
                  CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
                    IF(ABS(INTERPOLATED_POINT%VALUES(2,1))<ZERO_TOLERANCE) THEN
                      CALL FlagWarning("Mu is zero.",err,error,*999)
                    ELSE
                      FF=COORDINATE_SYSTEM%FOCUS*COORDINATE_SYSTEM%FOCUS
                      R=INTERPOLATED_POINT%VALUES(1,1)
                      MU=INTERPOLATED_POINT%VALUES(2,1)
                      G1=FF*(SINH(R)*SINH(R)+SIN(MU)*SIN(MU))
                      G3=FF*SINH(R)*SINH(R)*SIN(MU)*SIN(MU)
                      METRICS%dXidX(1:2,1)=METRICS%dXidX(1:2,1)*G1
                      METRICS%dXidX(1:2,2)=METRICS%dXidX(1:2,2)*G1
                      METRICS%dXidX(1:2,3)=METRICS%dXidX(1:2,3)*G3
                    ENDIF
                  CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE DEFAULT
                    LOCAL_ERROR="The coordinate system type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM%TYPE,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(LOCAL_ERROR,err,error,*999)
                  END SELECT
                ELSE
                  CALL FlagError("Invalid embedding of a surface in space.",err,error,*999)
                ENDIF
              CASE DEFAULT
                CALL FlagError("Invalid embedding in space.",err,error,*999)
              END SELECT
            ENDIF
          ELSE
            CALL FlagError("Metrics interpolated point has not been interpolated to include first derivatives.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Metrics interpolated point is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Metrics is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Coordinate system is not associated.",err,error,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Coordinate system metrics:",err,error,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Coordinate system type = ",TRIM(COORDINATE_SYSTEM_TYPE_STRING( &
        & COORDINATE_SYSTEM%TYPE)),err,error,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of X dimensions = ",METRICS%numberOfXDimensions,err,error,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of Xi dimensions = ",METRICS%numberOfXiDimensions,err,error,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Location of metrics:",err,error,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,METRICS%numberOfXDimensions,3,3,INTERPOLATED_POINT%VALUES(:,1), &
        & '("    X           :",3(X,E13.6))','(17X,3(X,E13.6))',err,error,*999)      
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Derivative of X wrt Xi:",err,error,*999)
      CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,METRICS%numberOfXDimensions,1,1,METRICS%numberOfXiDimensions, &
        & 3,3,METRICS%dXdXi,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    dX_dXi','(",I1,",:)',' :",3(X,E13.6))', &
        & '(17X,3(X,E13.6))',err,error,*999)
      IF(METRICS%numberOfXDimensions/=METRICS%numberOfXiDimensions) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Constructed derivative of X wrt Xi:",err,error,*999)
        SELECT CASE(METRICS%numberOfXiDimensions)
        CASE(1)
          !Line in space
          SELECT CASE(METRICS%numberOfXDimensions)
          CASE(2)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,METRICS%numberOfXDimensions,3,3,DX_DXI2, &
              & '("    dX_dXi(:,2) :",3(X,E13.6))','(17X,3(X,E13.6))',err,error,*999)      
          CASE(3)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,METRICS%numberOfXDimensions,3,3,DX_DXI2, &
              & '("    dX_dXi(:,2) :",3(X,E13.6))','(17X,3(X,E13.6))',err,error,*999)      
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,METRICS%numberOfXDimensions,3,3,DX_DXI3, &
              & '("    dX_dXi(:,3) :",3(X,E13.6))','(17X,3(X,E13.6))',err,error,*999)      
          CASE DEFAULT
            CALL FlagError("Invalid embedding of a line in space.",err,error,*999)
          END SELECT
        CASE(2)
          !Surface in space
          SELECT CASE(METRICS%numberOfXDimensions)
          CASE(3)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,METRICS%numberOfXDimensions,3,3,DX_DXI3, &
              & '("    dX_dXi(:,3) :",3(X,E13.6))','(17X,3(X,E13.6))',err,error,*999)      
          CASE DEFAULT
            CALL FlagError("Invalid embedding of a surface in space.",err,error,*999)
          END SELECT
        CASE DEFAULT
          CALL FlagError("Invalid embedding in space.",err,error,*999)
        END SELECT
      ENDIF
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  det dX_dXi    = ",DET_DX_DXI,err,error,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Derivative of Xi wrt X:",err,error,*999)
      CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,METRICS%numberOfXiDimensions,1,1,METRICS%numberOfXDimensions, &
        & 3,3,METRICS%dXidX,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    dXi_dX','(",I1,",:)',' :",3(X,E13.6))', &
        & '(17X,3(X,E13.6))',err,error,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Covariant metric tensor:",err,error,*999)
      CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,METRICS%numberOfXiDimensions,1,1,METRICS%numberOfXiDimensions, &
        & 3,3,METRICS%GL,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    GL','(",I1,",:)','     :",3(X,E13.6))','(17X,3(X,E13.6))', &
        & err,error,*999)      
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Contravariant metric tensor:",err,error,*999)
      CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,METRICS%numberOfXiDimensions,1,1,METRICS%numberOfXiDimensions, &
        & 3,3,METRICS%GU,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    GU','(",I1,",:)','     :",3(X,E13.6))','(17X,3(X,E13.6))', &
        & err,error,*999)      
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Jacobian type = ",METRICS%jacobianType,err,error,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Jacobian      = ",METRICS%jacobian,err,error,*999)
    ENDIF
    
    EXITS("COORDINATE_METRICS_CALCULATE")
    RETURN
999 ERRORSEXITS("COORDINATE_METRICS_CALCULATE",err,error)
    RETURN 1
  END SUBROUTINE COORDINATE_METRICS_CALCULATE

  !
  !================================================================================================================================
  !

  !>Calculates the normal vector, N, at the point X. IF REVERSE is true the reversed normal is returned. Old-cmiss-name: NORMAL
  SUBROUTINE COORDINATE_SYSTEM_NORMAL_CALCULATE(COORDINATE_SYSTEM,REVERSE,X,N,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: COORDINATE_SYSTEM !<The coordinate system to calculate the normal for
    LOGICAL, INTENT(IN) :: REVERSE !<If .TRUE. the reversed normal is returned.
    REAL(DP), INTENT(IN) :: X(:,:) !<The coordinate and it's derivatives to calcualte the normal at
    REAL(DP), INTENT(OUT) :: N(3) !<On exit, the normal vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfXDimensions,d_s1,d_s2,d2_s1
    REAL(DP) :: LENGTH,R,TANGENT1(3),TANGENT2(3)
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("COORDINATE_SYSTEM_NORMAL_CALCULATE",err,error,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%coordinateSystemFinished) THEN
        CALL FlagError("Coordinate system has been finished.",err,error,*999)
      ELSE
        numberOfXDimensions=COORDINATE_SYSTEM%numberOfDimensions
        d_s1=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(1)
        d_s2=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(2)
        d2_s1=PARTIAL_DERIVATIVE_SECOND_DERIVATIVE_MAP(1)
        SELECT CASE(COORDINATE_SYSTEM%TYPE)
        CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
          IF(numberOfXDimensions==2) THEN
            TANGENT1(1)=X(1,d_s1)
            TANGENT1(2)=X(2,d_s1)
          ELSE IF(numberOfXDimensions==3) THEN
            TANGENT1(1)=X(1,d_s1)
            TANGENT1(2)=X(2,d_s1)
            TANGENT1(3)=X(3,d_s1)
            TANGENT2(1)=X(1,d_s2)
            TANGENT2(2)=X(2,d_s2)
            TANGENT2(3)=X(3,d_s2)
          ELSE
            LOCAL_ERROR=TRIM(NUMBER_TO_VSTRING(numberOfXDimensions,"*",err,error))// &
              & " is an invalid number of dimensions to calculate a normal from in a rectangular cartesian coordinate system."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          ENDIF
        CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
          R=X(1,1)
          IF(numberOfXDimensions==2) THEN
            TANGENT1(1)=X(1,d_s1)*COS(X(1,d_s1))-R*SIN(X(1,d_s1))*X(2,d_s1)
            TANGENT1(2)=X(2,d_s1)*SIN(X(1,d_s1))+R*COS(X(1,d_s1))*X(2,d_s1)
          ELSE IF(numberOfXDimensions==3) THEN
            TANGENT1(1)=X(1,d_s1)*COS(X(1,d_s1))-R*SIN(X(1,d_s1))*X(2,d_s1)
            TANGENT1(2)=X(2,d_s1)*SIN(X(1,d_s1))+R*COS(X(1,d_s1))*X(2,d_s1)
            TANGENT1(3)=X(3,d_s1)
            TANGENT2(1)=X(1,d_s2)*COS(X(1,d_s1))-R*SIN(X(1,d_s1))*X(2,d_s2)
            TANGENT2(2)=X(1,d_s2)*SIN(X(1,d_s1))+R*COS(X(1,d_s1))*X(2,d_s2)
            TANGENT2(3)=X(3,d_s2)
           ELSE
            LOCAL_ERROR=TRIM(NUMBER_TO_VSTRING(numberOfXDimensions,"*",err,error))// &
              & " is an invalid number of dimensions to calculate a normal from in a rectangular cartesian coordinate system."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          ENDIF          
        CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
          R=X(1,1)
          TANGENT1(1)=X(1,d_s1)*COS(X(1,d2_s1))*COS(X(1,d_s1))- &
            &                 R*SIN(X(1,d2_s1))*COS(X(1,d_s1))*X(3,d_s1)- &
            &                 R*COS(X(1,d2_s1))*SIN(X(1,d_s1))*X(2,d_s1)
          TANGENT1(2)=X(1,d_s1)*COS(X(1,d2_s1))*SIN(X(1,d_s1))- &
            &                 R*SIN(X(1,d2_s1))*SIN(X(1,d_s1))*X(3,d_s1)+ &
            &                 R*COS(X(1,d2_s1))*COS(X(1,d_s1))*X(2,d_s1)
          TANGENT1(3)=X(1,d_s1)*SIN(X(1,d2_s1))+R*COS(X(1,d2_s1))*X(3,d_s1)
          TANGENT2(1)=X(1,d_s2)*COS(X(1,d2_s1))*COS(X(1,d_s1))- &
            &                 R*SIN(X(1,d2_s1))*COS(X(1,d_s1))*X(3,d_s2)- &
            &                 R*COS(X(1,d2_s1))*SIN(X(1,d_s1))*X(2,d_s2)
          TANGENT2(2)=X(1,d_s2)*COS(X(1,d2_s1))*SIN(X(1,d_s1))- &
            &                 R*SIN(X(1,d2_s1))*SIN(X(1,d_s1))*X(3,d_s2)+ &
            &                 R*COS(X(1,d2_s1))*COS(X(1,d_s1))*X(2,d_s2)
          TANGENT2(3)=X(1,d_s2)*SIN(X(1,d2_s1))+R*COS(X(1,d2_s1))*X(3,d_s2)
        CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          LOCAL_ERROR="The coordinate system type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM%TYPE,"*",err,error))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
        IF(numberOfXDimensions==2) THEN
          N(1)=-TANGENT1(2)
          N(2)=TANGENT1(1)
          LENGTH=SQRT(N(1)*N(1)+N(2)*N(2))
          IF(ABS(LENGTH)<ZERO_TOLERANCE) CALL FlagError("Zero normal vector length.",err,error,*999)
          IF(REVERSE) THEN
            N(1)=-N(1)/LENGTH
            N(2)=-N(2)/LENGTH
          ELSE            
            N(1)=N(1)/LENGTH
            N(2)=N(2)/LENGTH
          ENDIF
        ELSE
          N(1)=TANGENT1(2)*TANGENT2(3)-TANGENT1(3)*TANGENT2(2)
          N(2)=TANGENT1(3)*TANGENT2(1)-TANGENT1(1)*TANGENT2(3)
          N(3)=TANGENT1(1)*TANGENT2(2)-TANGENT1(2)*TANGENT2(1)
          LENGTH=SQRT(N(1)*N(1)+N(2)*N(2)+N(3)*N(3))
          IF(REVERSE) THEN
            N(1)=-N(1)/LENGTH
            N(2)=-N(2)/LENGTH
            N(3)=-N(3)/LENGTH
          ELSE            
            N(1)=N(1)/LENGTH
            N(2)=N(2)/LENGTH
            N(3)=N(3)/LENGTH
          ENDIF
        ENDIF        
      ENDIF
    ELSE
      CALL FlagError("Coordinate system is not associated.",err,error,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Coordinate system metrics:",err,error,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Coordinate system type = ",TRIM(COORDINATE_SYSTEM_TYPE_STRING( &
        & COORDINATE_SYSTEM%TYPE)),err,error,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of X dimensions = ",numberOfXDimensions,err,error,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfXDimensions,3,3,X(:,1),'("  X         :",3(X,E13.6))', &
        & '(13X,3(X,E13.6))',err,error,*999)      
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfXDimensions,3,3,TANGENT1,'("  Tangent 1 :",3(X,E13.6))', &
        & '(13X,3(X,E13.6))',err,error,*999)
      IF(numberOfXDimensions==3) THEN
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfXDimensions,3,3,TANGENT2,'("  Tangent 2 :",3(X,E13.6))', &
          & '(13X,3(X,E13.6))',err,error,*999)
      ENDIF
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfXDimensions,3,3,N,'("  Normal    :",3(X,E13.6))', &
        & '(13X,3(X,E13.6))',err,error,*999)            
    ENDIF
    
    EXITS("COORDINATE_SYSTEM_NORMAL_CALCULATE")
    RETURN
999 ERRORSEXITS("COORDINATE_SYSTEM_NORMAL_CALCULATE",err,error)
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_NORMAL_CALCULATE

   !
  !================================================================================================================================
  !

  !>Finalises a coordinate system and deallocates all memory. 
  SUBROUTINE COORDINATE_SYSTEM_FINALISE(COORDINATE_SYSTEM,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("COORDINATE_SYSTEM_FINALISE",err,error,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      DEALLOCATE(COORDINATE_SYSTEM)
    ENDIF
   
    EXITS("COORDINATE_SYSTEM_FINALISE")
    RETURN
999 ERRORSEXITS("COORDINATE_SYSTEM_FINALISE",err,error)
    RETURN 1

  END SUBROUTINE COORDINATE_SYSTEM_FINALISE

  !
  !================================================================================================================================
  !

  !>Returns the coordinate system focus. 
  SUBROUTINE COORDINATE_SYSTEM_FOCUS_GET(COORDINATE_SYSTEM,FOCUS,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to get the focus for
    REAL(DP), INTENT(OUT) :: FOCUS !<On return, the focus of the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("COORDINATE_SYSTEM_FOCUS_GET",err,error,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%coordinateSystemFinished) THEN
        SELECT CASE(COORDINATE_SYSTEM%TYPE)
        CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE,COORDINATE_OBLATE_SPHEROIDAL_TYPE)
          FOCUS=COORDINATE_SYSTEM%FOCUS
        CASE DEFAULT
          CALL FlagError("No focus defined for this coordinate system type.",err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Coordinate system has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Coordinate system is not associated.",err,error,*999)
    ENDIF
    
    EXITS("COORDINATE_SYSTEM_FOCUS_GET")
    RETURN
999 ERRORSEXITS("COORDINATE_SYSTEM_FOCUS_GET",err,error)
    RETURN 1
    
  END SUBROUTINE COORDINATE_SYSTEM_FOCUS_GET

  !
  !================================================================================================================================
  !


  !>Gets the coordinate system radial interpolation type. 
  SUBROUTINE Coordinates_RadialInterpolationTypeGet(coordinateSystem,radialInterpolationType,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<The coordinate system to get the radial interpolation for
    INTEGER(INTG), INTENT(OUT) :: radialInterpolationType !<On return, the radial interpolation type for the coordinate system \see COORDINATE_ROUTINES_RadialInterpolations,COORDINATE_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Coordinates_RadialInterpolationTypeGet",err,error,*999)

    IF(ASSOCIATED(coordinateSystem)) THEN
      IF(coordinateSystem%coordinateSystemFinished) THEN
        SELECT CASE(coordinateSystem%TYPE)
        CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE,COORDINATE_SPHERICAL_POLAR_TYPE)
          radialInterpolationType=coordinateSystem%radialInterpolationType
        CASE DEFAULT
          CALL FlagError("No radial interpolation type defined for this coordinate system interpolation.",err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Coordinate system has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Coordinate system is not associated.",err,error,*999)
    ENDIF
    
    EXITS("Coordinates_RadialInterpolationTypeGet")
    RETURN
999 ERRORSEXITS("Coordinates_RadialInterpolationTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Coordinates_RadialInterpolationTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the coordinate system type. 
  SUBROUTINE COORDINATE_SYSTEM_TYPE_GET(COORDINATE_SYSTEM,SYSTEM_TYPE,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to get the type for
    INTEGER(INTG), INTENT(OUT) :: SYSTEM_TYPE !<On return, the type for the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("COORDINATE_SYSTEM_TYPE_GET",err,error,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%coordinateSystemFinished) THEN
        SYSTEM_TYPE=COORDINATE_SYSTEM%TYPE
      ELSE
        CALL FlagError("Coordinate system has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Coordinate system is not associated.",err,error,*999)
    ENDIF
    
    EXITS("COORDINATE_SYSTEM_TYPE_GET")
    RETURN
999 ERRORSEXITS("COORDINATE_SYSTEM_TYPE_GET",err,error)
    RETURN 1

  END SUBROUTINE COORDINATE_SYSTEM_TYPE_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the dimension of the coordinate system. \see OPENCMISS::CMISSCoordinateSystemDimensionSet
  SUBROUTINE COORDINATE_SYSTEM_DIMENSION_SET(COORDINATE_SYSTEM,DIMENSION,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: COORDINATE_SYSTEM !<A pointer the coordinate system to set the dimension for
    INTEGER(INTG), INTENT(IN) :: DIMENSION !<The dimension to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("COORDINATE_SYSTEM_DIMENSION_SET",err,error,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%coordinateSystemFinished) THEN
        CALL FlagError("Coordinate system has been finished.",err,error,*999)
      ELSE
        SELECT CASE(COORDINATE_SYSTEM%TYPE)
        CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
          IF(DIMENSION>=1.AND.DIMENSION<=3) THEN
            COORDINATE_SYSTEM%numberOfDimensions=DIMENSION
          ELSE
            CALL FlagError("Invalid number of dimensions.",err,error,*999)
          ENDIF
        CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
          IF(DIMENSION>=2.AND.DIMENSION<=3) THEN
            COORDINATE_SYSTEM%numberOfDimensions=DIMENSION
          ELSE
            CALL FlagError("Invalid number of dimensions.",err,error,*999)
          ENDIF
        CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
          IF(DIMENSION==3) THEN
            COORDINATE_SYSTEM%numberOfDimensions=DIMENSION
          ELSE
            CALL FlagError("Invalid number of dimensions.",err,error,*999)
          ENDIF
        CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
          IF(DIMENSION==3) THEN
            COORDINATE_SYSTEM%numberOfDimensions=DIMENSION
          ELSE
            CALL FlagError("Invalid number of dimensions.",err,error,*999)
          ENDIF
        CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
          IF(DIMENSION==3) THEN
            COORDINATE_SYSTEM%numberOfDimensions=DIMENSION
          ELSE
            CALL FlagError("Invalid number of dimensions.",err,error,*999)
          ENDIF
        CASE DEFAULT
          CALL FlagError("Invalid coordinate system type.",err,error,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Coordinate system is not associated.",err,error,*999)
    ENDIF
    
    EXITS("COORDINATE_SYSTEM_DIMENSION_SET")
    RETURN
999 ERRORSEXITS("COORDINATE_SYSTEM_DIMENSION_SET",err,error)
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_DIMENSION_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the focus of a coordinate system. \see OPENCMISS::CMISSCoordinateSystemFocusSet
  SUBROUTINE COORDINATE_SYSTEM_FOCUS_SET(COORDINATE_SYSTEM,FOCUS,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to set the focus for
    REAL(DP), INTENT(IN) :: FOCUS !<The focus to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("COORDINATE_SYSTEM_FOCUS_SET",err,error,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%coordinateSystemFinished) THEN
        CALL FlagError("Coordinate system has been finished.",err,error,*999)
      ELSE
        SELECT CASE(COORDINATE_SYSTEM%TYPE)
        CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
          IF(FOCUS>ZERO_TOLERANCE) THEN
            COORDINATE_SYSTEM%FOCUS=FOCUS
          ELSE
            CALL FlagError("Focus is less than zero.",err,error,*999)
          ENDIF
        CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
          IF(FOCUS>ZERO_TOLERANCE) THEN
            COORDINATE_SYSTEM%FOCUS=FOCUS
          ELSE
            CALL FlagError("Focus is less than zero.",err,error,*999)
          ENDIF
        CASE DEFAULT
          CALL FlagError("Invalid coordinate system type.",err,error,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Coordinate system is not associated.",err,error,*999)
    ENDIF
      
    EXITS("COORDINATE_SYSTEM_FOCUS_SET")
    RETURN
999 ERRORSEXITS("COORDINATE_SYSTEM_FOCUS_SET",err,error)
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_FOCUS_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the radial interpolation type of a coordinate system. \see OPENCMISS::CMISSCoordinateSystemRadialInterpolationTypeSet
  SUBROUTINE Coordinates_RadialInterpolationTypeSet(coordinateSystem,radialInterpolationType,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<The coordinate system to set the interpolation type for
    INTEGER(INTG), INTENT(IN) :: radialInterpolationType !<The interpolation type to set \see COORDINATE_ROUTINES_RadialInterpolations,COORDINATE_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Coordinates_RadialInterpolationTypeSet",err,error,*999)

    IF(ASSOCIATED(coordinateSystem)) THEN
      IF(coordinateSystem%coordinateSystemFinished) THEN
        CALL FlagError("Coordinate system has been finished.",err,error,*999)
      ELSE
        SELECT CASE(coordinateSystem%TYPE)
        CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
          SELECT CASE(radialInterpolationType)
          CASE(COORDINATE_NO_RADIAL_INTERPOLATION_TYPE)
            coordinateSystem%radialInterpolationType=COORDINATE_NO_RADIAL_INTERPOLATION_TYPE
          CASE DEFAULT
            localERROR="The radial interpolation type of "//TRIM(NUMBER_TO_VSTRING(radialInterpolationType,"*",err,error))// &
              & " is invalid for a rectangular cartesian coordinate system."
            CALL FlagError(localERROR,err,error,*999)
          END SELECT
        CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE,COORDINATE_SPHERICAL_POLAR_TYPE)
          SELECT CASE(radialInterpolationType)
          CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
            coordinateSystem%radialInterpolationType=COORDINATE_RADIAL_INTERPOLATION_TYPE
          CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
            coordinateSystem%radialInterpolationType=COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE
          CASE DEFAULT
            localERROR="The radial interpolation type of "//TRIM(NUMBER_TO_VSTRING(radialInterpolationType,"*",err,error))// &
              & " is invalid for a cylindrical/spherical coordinate system."
            CALL FlagError(localERROR,err,error,*999)
          END SELECT
        CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
          SELECT CASE(radialInterpolationType)
          CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
            coordinateSystem%radialInterpolationType=COORDINATE_RADIAL_INTERPOLATION_TYPE
          CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
            coordinateSystem%radialInterpolationType=COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE
          CASE(COORDINATE_RADIAL_CUBED_INTERPOLATION_TYPE)
            coordinateSystem%radialInterpolationType=COORDINATE_RADIAL_CUBED_INTERPOLATION_TYPE
          CASE DEFAULT
            localERROR="The radial interpolation type of "//TRIM(NumberToVString(radialInterpolationType,"*",err,error))// &
              & " is invalid for a prolate spheroidal coordinate system."
            CALL FlagError(localERROR,err,error,*999)
          END SELECT
        CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          CALL FlagError("Invalid coordinate system type.",err,error,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Coordinate system is not associated.",err,error,*999)
    ENDIF
    
    EXITS("Coordinates_RadialInterpolationTypeSet")
    RETURN
999 ERRORSEXITS("Coordinates_RadialInterpolationTypeSet",err,error)
    RETURN 1
  END SUBROUTINE Coordinates_RadialInterpolationTypeSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of a coordinate system. \see OPENCMISS::CMISSCoordinateSystemTypeSet
  SUBROUTINE COORDINATE_SYSTEM_TYPE_SET(COORDINATE_SYSTEM,TYPE,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: COORDINATE_SYSTEM !<The coordinate system to set the type for
    INTEGER(INTG), INTENT(IN) :: TYPE !<The coordinate system type to set \see COORDINATE_ROUTINES_CoordinateSystemTypes,COORDINATE_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("COORDINATE_SYSTEM_TYPE_SET",err,error,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%coordinateSystemFinished) THEN
        CALL FlagError("Coordinate system has been finished.",err,error,*999)
      ELSE
        SELECT CASE(TYPE)
        CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
          COORDINATE_SYSTEM%TYPE=COORDINATE_RECTANGULAR_CARTESIAN_TYPE
        CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
          COORDINATE_SYSTEM%TYPE=COORDINATE_CYLINDRICAL_POLAR_TYPE
        CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
          COORDINATE_SYSTEM%TYPE=COORDINATE_SPHERICAL_POLAR_TYPE
        CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
          COORDINATE_SYSTEM%TYPE=COORDINATE_PROLATE_SPHEROIDAL_TYPE
        CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
          COORDINATE_SYSTEM%TYPE=COORDINATE_OBLATE_SPHEROIDAL_TYPE
        CASE DEFAULT
          CALL FlagError("Invalid coordinate system type.",err,error,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Coordinate system is not associated.",err,error,*999)
    ENDIF
    
    EXITS("COORDINATE_SYSTEM_TYPE_SET")
    RETURN
999 ERRORSEXITS("COORDINATE_SYSTEM_TYPE_SET",err,error)
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Returns the origin of a coordinate system. \see OPENCMISS::CMISSCoordinateSystemOriginGet
  SUBROUTINE COORDINATE_SYSTEM_ORIGIN_GET(COORDINATE_SYSTEM,ORIGIN,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to get the origin for
    REAL(DP), INTENT(OUT) :: ORIGIN(:) !<On return, the origin of the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("COORDINATE_SYSTEM_ORIGIN_GET",err,error,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%coordinateSystemFinished) THEN
        IF(SIZE(ORIGIN)>=3) THEN
          ORIGIN(1:3)=COORDINATE_SYSTEM%ORIGIN
        ELSE
          CALL FlagError("The origin must have >= 3 components.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Coordinate system has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Coordinate system is not associated.",err,error,*999)
    ENDIF
    
    EXITS("COORDINATE_SYSTEM_ORIGIN_GET")
    RETURN
999 ERRORSEXITS("COORDINATE_SYSTEM_ORIGIN_GET",err,error)
    RETURN 1
    
  END SUBROUTINE COORDINATE_SYSTEM_ORIGIN_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the origin of a coordinate system. \see OPENCMISS::CMISSCoordinateSystemOriginSet
  SUBROUTINE COORDINATE_SYSTEM_ORIGIN_SET(COORDINATE_SYSTEM,ORIGIN,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to set the origin for
    REAL(DP), INTENT(IN) :: ORIGIN(:) !<The origin to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("COORDINATE_SYSTEM_ORIGIN_SET",err,error,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%coordinateSystemFinished) THEN
        CALL FlagError("Coordinate system has been finished.",err,error,*999)
      ELSE
        IF(SIZE(ORIGIN)==3) THEN
          COORDINATE_SYSTEM%ORIGIN=ORIGIN
        ELSE
          CALL FlagError("The origin must have exactly 3 components.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Coordinate system is not associated.",err,error,*999)
    ENDIF
    
    EXITS("COORDINATE_SYSTEM_ORIGIN_SET")
    RETURN
999 ERRORSEXITS("COORDINATE_SYSTEM_ORIGIN_SET",err,error)
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_ORIGIN_SET

  !
  !================================================================================================================================
  !

  !>Returns the orientation of a coordinate system. \see OPENCMISS::CMISSCoordinateSystemOrientationSets/changesGet
  SUBROUTINE COORDINATE_SYSTEM_ORIENTATION_GET(COORDINATE_SYSTEM,ORIENTATION,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to get the orientation for
    REAL(DP), INTENT(OUT) :: ORIENTATION(:,:) !<On return, the orientation of the coordinate system
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("COORDINATE_SYSTEM_ORIENTATION_GET",err,error,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%coordinateSystemFinished) THEN
        IF(SIZE(ORIENTATION,1)>=3.AND.SIZE(ORIENTATION,2)>=3) THEN
          ORIENTATION(1:3,1:3)=COORDINATE_SYSTEM%ORIENTATION
        ELSE
          CALL FlagError("The orientation matrix must have >= 3x3 components.",err,error,*999)
        ENDIF
      ELSE
         CALL FlagError("Coordinate system has not been finished.",err,error,*999)
       ENDIF
    ELSE
      CALL FlagError("Coordinate system is not associated.",err,error,*999)
    ENDIF
    
    EXITS("COORDINATE_SYSTEM_ORIENTATION_GET")
    RETURN
999 ERRORSEXITS("COORDINATE_SYSTEM_ORIENTATION_GET",err,error)
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_ORIENTATION_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the orientation of a coordinate system. \see OPENCMISS::CMISSCoordinateSystemOrientationSet
  SUBROUTINE COORDINATE_SYSTEM_ORIENTATION_SET(COORDINATE_SYSTEM,ORIENTATION,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to set the orientation for
    REAL(DP), INTENT(IN) :: ORIENTATION(:,:) !<The orientation to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("COORDINATE_SYSTEM_ORIENTATION_SET",err,error,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%coordinateSystemFinished) THEN
        CALL FlagError("Coordinate system has been finished.",err,error,*999)
      ELSE
        IF(SIZE(ORIENTATION,1)==3.AND.SIZE(ORIENTATION,2)==3) THEN
!!TODO: \todo Check orientation matrix vectors are orthogonal to each other etc.
          COORDINATE_SYSTEM%ORIENTATION=ORIENTATION
        ELSE
          CALL FlagError("The orientation matrix must have exactly 3x3 components.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Coordinate system is not associated.",err,error,*999)
    ENDIF
    
    EXITS("COORDINATE_SYSTEM_ORIENTATION_SET")
    RETURN
999 ERRORSEXITS("COORDINATE_SYSTEM_ORIENTATION_SET",err,error)
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_ORIENTATION_SET
  
  !
  !================================================================================================================================
  !

  !>Starts the creation of and initialises a new coordinate system. \see OPENCMISS::CMISSCoordinateSystemCreateStart
  !>The default values of the COORDINATE_SYSTEM's attributes are:
  !>- TYPE: 1 (COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
  !>- radialInterpolationType: 0 (COORDINATE_NO_RADIAL_INTERPOLATION_TYPE)
  !>- Dimensions: 3
  !>- Focus: 1.0
  !>- Origin: (0.0,0.0,0.0)
  !>- Oritention: ((1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0))
  SUBROUTINE COORDINATE_SYSTEM_CREATE_START(USER_NUMBER,coordinateSystems,COORDINATE_SYSTEM,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number for the created coordinate system
    TYPE(CoordinateSystemsType), POINTER :: coordinateSystems !<A pointer to the coordinate systems to create the coordinate system for.
    TYPE(CoordinateSystemType), POINTER :: COORDINATE_SYSTEM !<On exit, a pointer to the created coordinate system. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coord_system_idx
    TYPE(CoordinateSystemType), POINTER :: NEW_COORDINATE_SYSTEM
    TYPE(CoordinateSystemPtrType), POINTER :: NEW_COORDINATE_SYSTEMS(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    NULLIFY(NEW_COORDINATE_SYSTEM)
    NULLIFY(NEW_COORDINATE_SYSTEMS)

    ENTERS("COORDINATE_SYSTEM_CREATE_START",err,error,*998)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) CALL FlagError("Coordinate system is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(coordinateSystems)) CALL FlagError("Coordinate systems is not associated.",err,error,*999)
    
    NULLIFY(NEW_COORDINATE_SYSTEM)
    CALL CoordinateSystem_UserNumberFind(coordinateSystems,USER_NUMBER,NEW_COORDINATE_SYSTEM,err,error,*999)
    IF(ASSOCIATED(NEW_COORDINATE_SYSTEM)) THEN
      LOCAL_ERROR="Coordinate system number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",err,error))// &
        & " has already been created."
      CALL FlagError(LOCAL_ERROR,err,error,*998)
    ELSE
      IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
        CALL FlagError("Coordinate system is already associated.",err,error,*999)
      ELSE
        NULLIFY(NEW_COORDINATE_SYSTEM)
        ALLOCATE(NEW_COORDINATE_SYSTEM,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new coordinate system.",err,error,*999)
      
        NEW_COORDINATE_SYSTEM%userNumber=USER_NUMBER
        NEW_COORDINATE_SYSTEM%coordinateSystems=>coordinateSystems
        NEW_COORDINATE_SYSTEM%coordinateSystemFinished=.FALSE.
        NEW_COORDINATE_SYSTEM%TYPE=COORDINATE_RECTANGULAR_CARTESIAN_TYPE
        NEW_COORDINATE_SYSTEM%radialInterpolationType=COORDINATE_NO_RADIAL_INTERPOLATION_TYPE
        NEW_COORDINATE_SYSTEM%numberOfDimensions=3
        NEW_COORDINATE_SYSTEM%FOCUS=1.0_DP    
        NEW_COORDINATE_SYSTEM%ORIGIN=[0.0_DP,0.0_DP,0.0_DP]
        NEW_COORDINATE_SYSTEM%ORIENTATION=RESHAPE(&
          & [1.0_DP,0.0_DP,0.0_DP, &
          &   0.0_DP,1.0_DP,0.0_DP, &
          &   0.0_DP,0.0_DP,1.0_DP], &
          & [3,3])
        
        ALLOCATE(NEW_COORDINATE_SYSTEMS(coordinateSystems%numberOfCoordinateSystems+1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new coordinate systems.",err,error,*999)
        DO coord_system_idx=1,coordinateSystems%numberOfCoordinateSystems
          NEW_COORDINATE_SYSTEMS(coord_system_idx)%ptr=>coordinateSystems%coordinateSystems(coord_system_idx)%ptr
        ENDDO !coord_system_idx
        NEW_COORDINATE_SYSTEMS(coordinateSystems%numberOfCoordinateSystems+1)%ptr=>NEW_COORDINATE_SYSTEM
        DEALLOCATE(coordinateSystems%coordinateSystems)
        coordinateSystems%coordinateSystems=>NEW_COORDINATE_SYSTEMS
        coordinateSystems%numberOfCoordinateSystems=coordinateSystems%numberOfCoordinateSystems+1
        
        COORDINATE_SYSTEM=>NEW_COORDINATE_SYSTEM
      ENDIF
    ENDIF
        
    EXITS("COORDINATE_SYSTEM_CREATE_START")
    RETURN
999 IF(ASSOCIATED(NEW_COORDINATE_SYSTEM)) DEALLOCATE(NEW_COORDINATE_SYSTEM)
    IF(ASSOCIATED(NEW_COORDINATE_SYSTEMS)) DEALLOCATE(NEW_COORDINATE_SYSTEMS)
    NULLIFY(COORDINATE_SYSTEM)
998 ERRORSEXITS("COORDINATE_SYSTEM_CREATE_START",err,error)
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_CREATE_START

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a new coordinate system.
  SUBROUTINE COORDINATE_SYSTEM_CREATE_FINISH(COORDINATE_SYSTEM,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to finish
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coord_system_idx
    TYPE(CoordinateSystemsType), POINTER :: coordinateSystems

    ENTERS("COORDINATE_SYSTEM_CREATE_FINISH",err,error,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      COORDINATE_SYSTEM%coordinateSystemFinished=.TRUE.
    ELSE
      CALL FlagError("Coordinate system is not associated.",err,error,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      NULLIFY(coordinateSystems)
      CALL CoordinateSystem_CoordinateSystemsGet(COORDINATE_SYSTEM,coordinateSystems,err,error,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of coordinate systems = ", &
        & coordinateSystems%numberOfCoordinateSystems,err,error,*999)
      DO coord_system_idx=1,coordinateSystems%numberOfCoordinateSystems
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Coordinate system : ",coord_system_idx,err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number = ", &
          & coordinateSystems%coordinateSystems(coord_system_idx)%ptr%userNumber,err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Type = ", &
          & COORDINATE_SYSTEM_TYPE_STRING(coordinateSystems%coordinateSystems(coord_system_idx)%ptr%TYPE),err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of dimensions = ", &
          & coordinateSystems%coordinateSystems(coord_system_idx)%ptr%numberOfDimensions,err,error,*999)
      ENDDO !coord_system_idx
    ENDIF
    
    EXITS("COORDINATE_SYSTEM_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("COORDINATE_SYSTEM_CREATE_FINISH",err,error)
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Destroys a coordinate system.
  SUBROUTINE COORDINATE_SYSTEM_DESTROY(COORDINATE_SYSTEM,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coord_system_no,new_coord_system_no
    LOGICAL :: FOUND
    TYPE(CoordinateSystemPtrType), POINTER :: NEW_COORDINATE_SYSTEMS(:)
    TYPE(CoordinateSystemsType), POINTER :: coordinateSystems

    ENTERS("COORDINATE_SYSTEM_DESTROY",err,error,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(COORDINATE_SYSTEM%userNumber==0) THEN
        CALL FlagError("Cannot destroy the world coordinate system.",err,error,*999)
      ELSE
        NULLIFY(coordinateSystems)
        CALL CoordinateSystem_CoordinateSystemsGet(COORDINATE_SYSTEM,coordinateSystems,err,error,*999)
        FOUND=.FALSE.
        new_coord_system_no=0
        ALLOCATE(NEW_COORDINATE_SYSTEMS(coordinateSystems%numberOfCoordinateSystems-1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new coordianate systems.",err,error,*999)
        DO coord_system_no=1,coordinateSystems%numberOfCoordinateSystems
          IF(coordinateSystems%coordinateSystems(coord_system_no)%ptr%userNumber==COORDINATE_SYSTEM%userNumber) THEN
            FOUND=.TRUE.
          ELSE
            new_coord_system_no=new_coord_system_no+1
            NEW_COORDINATE_SYSTEMS(new_coord_system_no)%ptr=>coordinateSystems%coordinateSystems(coord_system_no)%ptr
          ENDIF
        ENDDO !coord_system_no
        IF(FOUND) THEN
          CALL COORDINATE_SYSTEM_FINALISE(COORDINATE_SYSTEM,err,error,*999)
          DEALLOCATE(coordinateSystems%coordinateSystems)
          coordinateSystems%coordinateSystems=>NEW_COORDINATE_SYSTEMS
          coordinateSystems%numberOfCoordinateSystems=coordinateSystems%numberOfCoordinateSystems-1
        ELSE
          DEALLOCATE(NEW_COORDINATE_SYSTEMS)
          CALL FlagError("Coordinate system number to destroy does not exist.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Coordinate system is not associated.",err,error,*999)
    ENDIF
      
    EXITS("COORDINATE_SYSTEM_DESTROY")
    RETURN
999 ERRORSEXITS("COORDINATE_SYSTEM_DESTROY",err,error)
    RETURN 1
  END SUBROUTINE COORDINATE_SYSTEM_DESTROY

  !
  !================================================================================================================================
  !
  
  !>Calculates DX(:)/DZ(I) at X, where Z(I) are rectangular cartesian and X(:) are curvilinear coordinates defined by COORDINATE_SYSTEM for double precision coordinates.
  FUNCTION DXZ_DP(COORDINATE_SYSTEM,I,X,err,error)
  
    !Argument variables
    TYPE(CoordinateSystemType), INTENT(IN) :: COORDINATE_SYSTEM
    INTEGER(INTG), INTENT(IN) :: I
    REAL(DP), INTENT(IN) :: X(:)
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Function variable
    REAL(DP) :: DXZ_DP(SIZE(X,1))
    !Local variables
    REAL(DP) :: RD,FOCUS

    ENTERS("DXZ_DP",err,error,*999)

    DXZ_DP=0.0_DP

    IF(SIZE(X,1)<COORDINATE_SYSTEM%numberOfDimensions) &
      & CALL FlagError("Size of X is less than the number of dimensions", &
      & err,error,*999)
   
    SELECT CASE(COORDINATE_SYSTEM%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      IF(I>0.AND.I<=COORDINATE_SYSTEM%numberOfDimensions) THEN
        DXZ_DP(I)=1.0_DP
      ELSE
        CALL FlagError("Invalid i value",err,error,*999)
      ENDIF
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      SELECT CASE(COORDINATE_SYSTEM%numberOfDimensions)
      CASE(2)
        SELECT CASE(I)
        CASE(1)
          DXZ_DP(1)=COS(X(2))
          DXZ_DP(2)=-SIN(X(2))/X(1)
        CASE(2)
          DXZ_DP(1)=SIN(X(2))
          DXZ_DP(2)=COS(X(2))/X(1)
        CASE DEFAULT
          CALL FlagError("Invalid i value",err,error,*999)
        END SELECT
      CASE(3)
        SELECT CASE(I)
        CASE(1)
          DXZ_DP(1)=COS(X(2))
          DXZ_DP(2)=-SIN(X(2))/X(1)
          DXZ_DP(3)=0.0_DP
        CASE(2)
          DXZ_DP(1)=SIN(X(2))
          DXZ_DP(2)=COS(X(2))/X(1)
          DXZ_DP(3)=0.0_DP
        CASE(3)
          DXZ_DP(1)=0.0_DP
          DXZ_DP(2)=0.0_DP
          DXZ_DP(3)=1.0_DP
        CASE DEFAULT
          CALL FlagError("Invalid i value",err,error,*999)
        END SELECT
      CASE DEFAULT
        CALL FlagError("Invalid number of coordinates",err,error,*999)
      END SELECT
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(COORDINATE_SYSTEM%numberOfDimensions==3) THEN
        SELECT CASE(I)
        CASE(1)
          DXZ_DP(1)=COS(X(2))*COS(X(3))
          DXZ_DP(2)=-SIN(X(2))/(X(1)*COS(X(3)))
          DXZ_DP(3)=-COS(X(2))*SIN(X(3))/X(1)
        CASE(2)
          DXZ_DP(1)=SIN(X(2))*COS(X(3))
          DXZ_DP(2)=COS(X(2))/(X(1)*COS(X(3)))
          DXZ_DP(3)=-SIN(X(2))*SIN(X(3))/X(1)
        CASE(3)
          DXZ_DP(1)=SIN(X(3))
          DXZ_DP(2)=0.0_DP
          DXZ_DP(3)=COS(X(3))/X(1)
        CASE DEFAULT
          CALL FlagError("Invalid i value",err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Invalid number of coordinates",err,error,*999)
      ENDIF
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      IF(COORDINATE_SYSTEM%numberOfDimensions==3) THEN
        FOCUS=COORDINATE_SYSTEM%FOCUS
        RD=FOCUS*(COSH(X(1))*COSH(X(1))-COS(X(2))*COS(X(2)))
        SELECT CASE(I)
        CASE(1)
          DXZ_DP(1)=SINH(X(1))*COS(X(2))/RD
          DXZ_DP(2)=-COSH(X(1))*SIN(X(2))/RD
          DXZ_DP(3)=0.0_DP
        CASE(2)
          DXZ_DP(1)=COSH(X(1))*SIN(X(2))*COS(X(3))/RD
          DXZ_DP(2)=SINH(X(1))*COS(X(2))*COS(X(3))/RD
          DXZ_DP(3)=-SIN(X(3))/(FOCUS*SINH(X(1))*SIN(X(2)))
        CASE(3)
          DXZ_DP(1)=COSH(X(1))*SIN(X(2))*SIN(X(3))/RD
          DXZ_DP(2)=SINH(X(1))*COS(X(2))*SIN(X(3))/RD
          DXZ_DP(3)=COS(X(3))/(FOCUS*SINH(X(1))*SIN(X(2)))
        CASE DEFAULT
          CALL FlagError("Invalid i value",err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Invalid number of coordinates",err,error,*999)
      ENDIF
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      CALL FlagError("Not implemented",err,error,*999)
    CASE DEFAULT
      CALL FlagError("Invalid coordinate type",err,error,*999)
    END SELECT

    EXITS("DXZ_DP")
    RETURN
999 ERRORSEXITS("DXZ_DP",err,error)
    RETURN 
  END FUNCTION DXZ_DP

  !
  !================================================================================================================================
  !

  !#### Generic-Function: D2ZX
  !###  Description:
  !###    Calculates D2Z(:)/DX(I)DX(J) at X(:), where Z(:) are rectangalar
  !###    Cartesian and X(I) and X(J) are curvilinear coordinates defined by 
  !###    COORDINATE_SYSTEM.
  !###  Child-functions: D2ZX_SP,D2ZX_SP

  !
  !================================================================================================================================
  !
  
  FUNCTION D2ZX_DP(COORDINATE_SYSTEM,I,J,X,err,error)
  
    !#### Function: D2ZX_DP
    !###  Type: REAL(DP)(SIZE(X,1))
    !###  Description:
    !###    Calculates D2Z(:)/DX(I)DX(J) at X(:), where Z(:) are rectangalar
    !###    Cartesian and X(I) and X(J) are curvilinear coordinates defined by 
    !###    COORDINATE_SYSTEM.
    !###  Parent-function: D2ZX
    
    !Argument variables
    TYPE(CoordinateSystemType), INTENT(IN) :: COORDINATE_SYSTEM
    INTEGER(INTG), INTENT(IN) :: I,J
    REAL(DP), INTENT(IN) :: X(:)
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Function variable
    REAL(DP) :: D2ZX_DP(SIZE(X,1))
    !Local variables
    REAL(DP) :: FOCUS

    ENTERS("D2ZX_DP",err,error,*999)

    D2ZX_DP=0.0_DP

    IF(SIZE(X,1)<COORDINATE_SYSTEM%numberOfDimensions) &
      & CALL FlagError("Size of X is less than the number of dimensions", &
      & err,error,*999)
   
    SELECT CASE(COORDINATE_SYSTEM%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      D2ZX_DP(1:COORDINATE_SYSTEM%numberOfDimensions)=0.0_DP
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      SELECT CASE(COORDINATE_SYSTEM%numberOfDimensions)
      CASE(2)
        SELECT CASE(I)
        CASE(1)
          SELECT CASE(J)
          CASE(1)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=0.0_DP
          CASE(2)
            D2ZX_DP(1)=-SIN(X(2))
            D2ZX_DP(2)=COS(X(2))
          CASE DEFAULT
            CALL FlagError("Invalid j value",err,error,*999)
          END SELECT
        CASE(2)
          SELECT CASE(J)
          CASE(1)
            D2ZX_DP(1)=-SIN(X(2))
            D2ZX_DP(2)=COS(X(2))
          CASE(2)
            D2ZX_DP(1)=-X(1)*COS(X(2))
            D2ZX_DP(2)=-X(1)*SIN(X(2))
          CASE DEFAULT
            CALL FlagError("Invalid j value",err,error,*999)
          END SELECT
        CASE DEFAULT
          CALL FlagError("Invalid i value",err,error,*999)
        END SELECT
      CASE(3)
        SELECT CASE(I)
        CASE(1)
          SELECT CASE(J)
          CASE(1)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=0.0_DP
            D2ZX_DP(3)=0.0_DP
          CASE(2)
            D2ZX_DP(1)=-SIN(X(2))
            D2ZX_DP(2)=COS(X(2))
            D2ZX_DP(3)=0.0_DP
          CASE(3)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=0.0_DP
            D2ZX_DP(3)=0.0_DP
          CASE DEFAULT
            CALL FlagError("Invalid j value",err,error,*999)
          END SELECT
        CASE(2)
          SELECT CASE(J)
          CASE(1)
            D2ZX_DP(1)=-SIN(X(2))
            D2ZX_DP(2)=COS(X(2))
            D2ZX_DP(3)=0.0_DP
          CASE(2)
            D2ZX_DP(1)=-X(1)*COS(X(2))
            D2ZX_DP(2)=-X(1)*SIN(X(2))
            D2ZX_DP(3)=0.0_DP
          CASE(3)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=0.0_DP
            D2ZX_DP(3)=0.0_DP
          CASE DEFAULT
            CALL FlagError("Invalid j value",err,error,*999)
          END SELECT
        CASE(3)
          SELECT CASE(J)
          CASE(1)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=0.0_DP
            D2ZX_DP(3)=0.0_DP
          CASE(2)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=0.0_DP
            D2ZX_DP(3)=0.0_DP
          CASE(3)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=0.0_DP
            D2ZX_DP(3)=0.0_DP
          CASE DEFAULT
            CALL FlagError("Invalid j value",err,error,*999)
          END SELECT
        CASE DEFAULT
          CALL FlagError("Invalid i value",err,error,*999)
        END SELECT
      CASE DEFAULT
        CALL FlagError("Invalid number of coordinates",err,error,*999)
      END SELECT
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(COORDINATE_SYSTEM%numberOfDimensions==3) THEN
        SELECT CASE(I)
        CASE(1)
          SELECT CASE(J)
          CASE(1)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=0.0_DP
            D2ZX_DP(3)=0.0_DP
          CASE(2)
            D2ZX_DP(1)=-SIN(X(2))*COS(X(3))
            D2ZX_DP(2)=COS(X(2))*COS(X(3))
            D2ZX_DP(3)=0.0_DP
          CASE(3)
            D2ZX_DP(1)=-COS(X(2))*SIN(X(3))
            D2ZX_DP(2)=-SIN(X(2))*SIN(X(3))
            D2ZX_DP(3)=COS(X(3))
          CASE DEFAULT
            CALL FlagError("Invalid j value",err,error,*999)
          END SELECT
        CASE(2)
          SELECT CASE(J)
          CASE(1)
            D2ZX_DP(1)=-SIN(X(2))*COS(X(3))
            D2ZX_DP(2)=COS(X(2))*COS(X(3))
            D2ZX_DP(3)=0.0_DP
          CASE(2)
            D2ZX_DP(1)=-X(1)*COS(X(2))*COS(X(3))
            D2ZX_DP(2)=-X(1)*SIN(X(2))*COS(X(3))
            D2ZX_DP(3)=0.0_DP
          CASE(3)
            D2ZX_DP(1)=X(1)*SIN(X(2))*SIN(X(3))
            D2ZX_DP(2)=-X(1)*COS(X(2))*SIN(X(3))
            D2ZX_DP(3)=0.0_DP
          CASE DEFAULT
            CALL FlagError("Invalid j value",err,error,*999)
          END SELECT
        CASE(3)
          SELECT CASE(J)
          CASE(1)
            D2ZX_DP(1)=-COS(X(2))*SIN(X(3))
            D2ZX_DP(2)=-SIN(X(2))*SIN(X(3))
            D2ZX_DP(3)=COS(X(3))
          CASE(2)
            D2ZX_DP(1)=X(1)*SIN(X(2))*SIN(X(3))
            D2ZX_DP(2)=-X(1)*COS(X(2))*SIN(X(3))
            D2ZX_DP(3)=0.0_DP
          CASE(3)
            D2ZX_DP(1)=-X(1)*COS(X(2))*COS(X(3))
            D2ZX_DP(2)=-X(1)*SIN(X(2))*COS(X(3))
            D2ZX_DP(3)=-X(1)*SIN(X(3))
          CASE DEFAULT
            CALL FlagError("Invalid j value",err,error,*999)
          END SELECT
        CASE DEFAULT
          CALL FlagError("Invalid i value",err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Invalid number of coordinates",err,error,*999)
      ENDIF
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      IF(COORDINATE_SYSTEM%numberOfDimensions==3) THEN
        FOCUS=COORDINATE_SYSTEM%FOCUS
        SELECT CASE(I)
        CASE(1)
          SELECT CASE(J)
          CASE(1)
            D2ZX_DP(1)=FOCUS*COSH(X(1))*COS(X(2))          
            D2ZX_DP(2)=FOCUS*SINH(X(1))*SIN(X(2))*COS(X(3))
            D2ZX_DP(3)=FOCUS*SINH(X(1))*SIN(X(2))*SIN(X(3))
          CASE(2)
            D2ZX_DP(1)=-FOCUS*SINH(X(1))*SIN(X(2))
            D2ZX_DP(2)=FOCUS*COSH(X(1))*COS(X(2))*COS(X(3))
            D2ZX_DP(3)=FOCUS*COSH(X(1))*COS(X(2))*SIN(X(3))
          CASE(3)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=-FOCUS*COSH(X(1))*SIN(X(2))*SIN(X(3))
            D2ZX_DP(3)=FOCUS*COSH(X(1))*SIN(X(2))*COS(X(3))
          CASE DEFAULT
            CALL FlagError("Invalid j value",err,error,*999)
          END SELECT
        CASE(2)
          SELECT CASE(J)
          CASE(1)
            D2ZX_DP(1)=-FOCUS*SINH(X(1))*SIN(X(2))
            D2ZX_DP(2)=FOCUS*COSH(X(1))*COS(X(2))*COS(X(3))
            D2ZX_DP(3)=FOCUS*COSH(X(1))*COS(X(2))*SIN(X(3))
          CASE(2)
            D2ZX_DP(1)=-FOCUS*COSH(X(1))*COS(X(2))
            D2ZX_DP(2)=-FOCUS*SINH(X(1))*SIN(X(2))*COS(X(3))
            D2ZX_DP(3)=-FOCUS*SINH(X(1))*SIN(X(2))*SIN(X(3))
          CASE(3)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=-FOCUS*SINH(X(1))*COS(X(2))*SIN(X(3))
            D2ZX_DP(3)=FOCUS*SINH(X(1))*COS(X(2))*COS(X(3))
          CASE DEFAULT
            CALL FlagError("Invalid j value",err,error,*999)
          END SELECT
        CASE(3)
          SELECT CASE(J)
          CASE(1)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=-FOCUS*COSH(X(1))*SIN(X(2))*SIN(X(3))
            D2ZX_DP(3)=FOCUS*COSH(X(1))*SIN(X(2))*COS(X(3))
          CASE(2)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=-FOCUS*SINH(X(1))*COS(X(2))*SIN(X(3))
            D2ZX_DP(3)=FOCUS*SINH(X(1))*COS(X(2))*COS(X(3))
          CASE(3)
            D2ZX_DP(1)=0.0_DP
            D2ZX_DP(2)=-FOCUS*SINH(X(1))*SIN(X(2))*COS(X(3))
            D2ZX_DP(3)=-FOCUS*SINH(X(1))*SIN(X(2))*SIN(X(3))
          CASE DEFAULT
            CALL FlagError("Invalid j value",err,error,*999)
          END SELECT
        CASE DEFAULT
          CALL FlagError("Invalid i value",err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Invalid number of coordinates",err,error,*999)
      ENDIF
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      CALL FlagError("Not implemented",err,error,*999)
    CASE DEFAULT
      CALL FlagError("Invalid coordinate type",err,error,*999)
    END SELECT

    EXITS("D2ZX_DP")
    RETURN
999 ERRORSEXITS("D2ZX_DP",err,error)
    RETURN 
  END FUNCTION D2ZX_DP

  !
  !================================================================================================================================
  !
  
  !#### Generic-Function: DZX
  !###  Description:
  !###    Calculates DZ(:)/DX(J) at X, where Z(:) are rectangalar
  !###    Cartesian and X(J) are curvilinear coordinates defined by 
  !###    COORDINATE_SYSTEM.
  !###  Child-functions: DZX_DP,DZX_SP

  !
  !================================================================================================================================
  !
  
  FUNCTION DZX_DP(COORDINATE_SYSTEM,I,X,err,error)
  
    !#### Function: DZX_DP
    !###  Type: REAL(DP)(SIZE(X,1))
    !###  Description:
    !###    Calculates DZ(:)/DX(I) at X, where Z(:) are rectangalar
    !###    Cartesian and X(I) are curvilinear coordinates defined by 
    !###    COORDINATE_SYSTEM.
    !###  Parent-function: DZX
    
    !Argument variables
    TYPE(CoordinateSystemType), INTENT(IN) :: COORDINATE_SYSTEM
    INTEGER(INTG), INTENT(IN) :: I
    REAL(DP), INTENT(IN) :: X(:)
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Function variable
    REAL(DP) :: DZX_DP(SIZE(X,1))
    !Local variables
    REAL(DP) :: FOCUS

    ENTERS("DZX_DP",err,error,*999)

    DZX_DP=0.0_DP

    IF(SIZE(X,1)<COORDINATE_SYSTEM%numberOfDimensions) &
      & CALL FlagError("Size of X is less than the number of dimensions", &
      & err,error,*999)
   
   SELECT CASE(COORDINATE_SYSTEM%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      IF(I>0.AND.I<=COORDINATE_SYSTEM%numberOfDimensions) THEN
        DZX_DP(I)=1.0_DP
      ELSE
        CALL FlagError("Invalid i value",err,error,*999)
      ENDIF
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      SELECT CASE(COORDINATE_SYSTEM%numberOfDimensions)
      CASE(2)
        SELECT CASE(I)
        CASE(1)
          DZX_DP(1)=COS(X(2))
          DZX_DP(2)=SIN(X(2))
        CASE(2)
          DZX_DP(1)=-X(1)*SIN(X(2))
          DZX_DP(2)=X(1)*COS(X(2))
        CASE DEFAULT
          CALL FlagError("Invalid i value",err,error,*999)
        END SELECT
      CASE(3)
        SELECT CASE(I)
        CASE(1)
          DZX_DP(1)=COS(X(2))
          DZX_DP(2)=SIN(X(2))
          DZX_DP(3)=0.0_DP
        CASE(2)
          DZX_DP(1)=-X(1)*SIN(X(2))
          DZX_DP(2)=X(1)*COS(X(2))
          DZX_DP(3)=0.0_DP
        CASE(3)
          DZX_DP(1)=0.0_DP
          DZX_DP(2)=0.0_DP
          DZX_DP(3)=1.0_DP
        CASE DEFAULT
          CALL FlagError("Invalid i value",err,error,*999)
        END SELECT
      CASE DEFAULT
        CALL FlagError("Invalid number of coordinates",err,error,*999)
      END SELECT
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(COORDINATE_SYSTEM%numberOfDimensions==3) THEN
        SELECT CASE(I)
        CASE(1)
          DZX_DP(1)=COS(X(2))*COS(X(3))
          DZX_DP(2)=SIN(X(2))*COS(X(3))
          DZX_DP(3)=SIN(X(3))
        CASE(2)
          DZX_DP(1)=-X(1)*SIN(X(2))*COS(X(3))
          DZX_DP(2)=X(1)*COS(X(2))*COS(X(3))
          DZX_DP(3)=0.0_DP
        CASE(3)
          DZX_DP(1)=-X(1)*COS(X(2))*SIN(X(3))
          DZX_DP(2)=-X(1)*SIN(X(2))*SIN(X(3))
          DZX_DP(3)=X(1)*COS(X(3))
        CASE DEFAULT
          CALL FlagError("Invalid i value",err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Invalid number of coordinates",err,error,*999)
      ENDIF
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      IF(COORDINATE_SYSTEM%numberOfDimensions==3) THEN
        FOCUS=COORDINATE_SYSTEM%FOCUS
        SELECT CASE(I)
        CASE(1)
          DZX_DP(1)=FOCUS*SINH(X(1))*COS(X(2))
          DZX_DP(2)=FOCUS*COSH(X(1))*SIN(X(2))*COS(X(3))
          DZX_DP(3)=FOCUS*COSH(X(1))*SIN(X(2))*SIN(X(3))
        CASE(2)
          DZX_DP(1)=-FOCUS*COSH(X(1))*SIN(X(2))
          DZX_DP(2)=FOCUS*SINH(X(1))*COS(X(2))*COS(X(3))
          DZX_DP(3)=FOCUS*SINH(X(1))*COS(X(2))*SIN(X(3))
        CASE(3)
          DZX_DP(1)=0.0_DP
          DZX_DP(2)=-FOCUS*SINH(X(1))*SIN(X(2))*SIN(X(3))
          DZX_DP(3)=FOCUS*SINH(X(1))*SIN(X(2))*COS(X(3))
        CASE DEFAULT
          CALL FlagError("Invalid i value",err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Invalid number of coordinates",err,error,*999)
      ENDIF
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      IF(COORDINATE_SYSTEM%numberOfDimensions==3) THEN
        FOCUS=COORDINATE_SYSTEM%FOCUS
        SELECT CASE(I)
        CASE(1)
          DZX_DP(1)=FOCUS*SINH(X(1))*COS(X(2))*COS(X(3))
          DZX_DP(2)=FOCUS*COSH(X(1))*SIN(X(2))
          DZX_DP(3)=FOCUS*SINH(X(1))*COS(X(2))*SIN(X(3))
        CASE(2)
          DZX_DP(1)=-FOCUS*COSH(X(1))*SIN(X(2))*COS(X(3))
          DZX_DP(2)=FOCUS*SINH(X(1))*COS(X(2))
          DZX_DP(3)=-FOCUS*COSH(X(1))*SIN(X(2))*SIN(X(3))
        CASE(3)
          DZX_DP(1)=-FOCUS*COSH(X(1))*COS(X(2))*SIN(X(3))
          DZX_DP(2)=0.0_DP
          DZX_DP(3)=FOCUS*COSH(X(1))*COS(X(2))*COS(X(3))
        CASE DEFAULT
          CALL FlagError("Invalid i value",err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Invalid number of coordinates",err,error,*999)
      ENDIF
    CASE DEFAULT
      CALL FlagError("Invalid coordinate type",err,error,*999)
    END SELECT

    EXITS("DZX_DP")
    RETURN
999 ERRORSEXITS("DZX_DP",err,error)
    RETURN 
  END FUNCTION DZX_DP

  !
  !================================================================================================================================
  !

  !!TODO:: CHANGE THIS TO A FUNCTION
  
  !#### Generic-subroutine: COORDINATE_DERIVATIVE_CONVERT_TO_RC
  !###  Description:
  !###    COORDINATE_DERIVATIVE_CONVERT_TO_RC performs a coordinate transformation from a
  !###    coordinate system identified by COORDINATE at the point X
  !###    with coordinates/derivatives X(nj,nu) to the point with 
  !###    coordinates/derivatives Z(nj) in rectangular cartesian coordinates
  !###  Child-subroutines: COORDINATE_DERIVATIVE_CONVERT_TO_RC_DP,COORDINATE_DERIVATIVE_CONVERT_TO_RC_SP

  !
  !================================================================================================================================
  !
  
  SUBROUTINE COORDINATE_DERIVATIVE_CONVERT_TO_RC_DP(COORDINATE_SYSTEM,PART_DERIV_TYPE,X,Z,&
    & err,error,*)
  
    !#### Subroutine: COORDINATE_DERIVATIVE_CONVERT_TO_RC_DP
    !###  Description:
    !###    COORDINATE_DERIVATIVE_CONVERT_TO_RC_DP performs a coordinate transformation from a
    !###    coordinate system identified by COORDINATE_SYSTEM at the point X
    !###    with coordinates/derivatives X(nj,nu) to the point with 
    !###    coordinates/derivatives Z(nj) in rectangular cartesian coordinates
    !###    for double precision coordinates.
    !###  Parent-function: COORDINATE_DERIVATIVE_CONVERT_TO_RC
    
    !Argument variables
    TYPE(CoordinateSystemType), INTENT(IN) :: COORDINATE_SYSTEM
    INTEGER(INTG), INTENT(IN) :: PART_DERIV_TYPE
    REAL(DP), INTENT(IN) :: X(:,:)
    REAL(DP), INTENT(OUT) :: Z(:)
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local variables
    REAL(DP) :: FOCUS
    
    ENTERS("COORDINATE_DERIVATIVE_CONVERT_TO_RC_DP",err,error,*999)

!!TODO: change all second index X(:,?) numbers to their apropriate constant
!!as defined in constants e.g. X(1,2) == X(1,PART_DERIV_S1)
    
    IF(SIZE(X,1)<COORDINATE_SYSTEM%numberOfDimensions) &
      & CALL FlagError("Size of X is less than the number of dimensions", &
      & err,error,*999)
    
    IF(SIZE(X,1)==SIZE(Z,1)) THEN
      SELECT CASE(COORDINATE_SYSTEM%TYPE)
      CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
        IF(SIZE(X,2)>=PART_DERIV_TYPE) THEN
          Z=X(:,PART_DERIV_TYPE)
        ELSE
          CALL FlagError("Invalid derivative type",err,error,*999)
        ENDIF
      CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
        SELECT CASE(PART_DERIV_TYPE)
        CASE(NO_PART_DERIV)
          Z=COORDINATE_CONVERT_TO_RC(COORDINATE_SYSTEM,X(:,1),err,error)
          IF(err/=0) GOTO 999
        CASE(PART_DERIV_S1)
          IF(SIZE(X,2)>=2) THEN
            SELECT CASE(COORDINATE_SYSTEM%numberOfDimensions)
            CASE(2)
              Z(1)=COS(X(2,1))*X(1,2)-X(1,1)*SIN(X(2,1))*X(2,2) !d(x)/d(s1)
              Z(2)=SIN(X(2,1))*X(1,2)+X(1,1)*COS(X(2,1))*X(2,2) !d(y)/d(s1)
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,2)-X(1,1)*SIN(X(2,1))*X(2,2) !d(x)/d(s1)
              Z(2)=SIN(X(2,1))*X(1,2)+X(1,1)*COS(X(2,1))*X(2,2) !d(y)/d(s1)
              Z(3)=X(3,2) !d(z)/d(s1)
            CASE DEFAULT
              CALL FlagError("Invalid number of coordinates",err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Not enough X derivatives supplied",err,error,*999)
          ENDIF
        CASE(PART_DERIV_S1_S1)
          CALL FlagError("Not implemented",err,error,*999)
        CASE(PART_DERIV_S2)
          IF(SIZE(X,2)>=4) THEN
            SELECT CASE(COORDINATE_SYSTEM%numberOfDimensions)
            CASE(2)
              Z(1)=COS(X(2,1))*X(1,4)-X(1,1)*SIN(X(2,1))*X(2,4) !d(x)/d(s2)
              Z(2)=SIN(X(2,1))*X(1,4)+X(1,1)*COS(X(2,1))*X(2,4) !d(y)/d(s2)
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,4)-X(1,1)*SIN(X(2,1))*X(2,4) !d(x)/d(s2)
              Z(2)=SIN(X(2,1))*X(1,4)+X(1,1)*COS(X(2,1))*X(2,4) !d(y)/d(s2)
              Z(3)=X(3,4) !d(z)/d(s2)
            CASE DEFAULT
              CALL FlagError("Invalid number of coordinates",err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Not enough X derivatives supplied",err,error,*999)
          ENDIF
        CASE(PART_DERIV_S2_S2)
          CALL FlagError("Not implemented",err,error,*999)
        CASE(PART_DERIV_S1_S2)
          IF(SIZE(X,2)>=6) THEN
            SELECT CASE(SIZE(X,1))
            CASE(2)
              Z(1)=COS(X(2,1))*X(1,6)-X(1,2)*SIN(X(2,1))*X(2,4)-&
                & (X(1,4)*SIN(X(2,1))*X(2,2)+X(1,1)*COS(X(2,1))*&
                & X(2,2)*X(2,4)+X(1,1)*SIN(X(2,1))*X(2,6)) !d2(x)/d(s1)d(s2)
              Z(2)=SIN(X(2,1))*X(1,6)+X(1,2)*COS(X(2,1))*X(2,4)+&
                & (X(1,4)*COS(X(2,1))*X(2,2)-X(1,1)*SIN(X(2,1))*&
                & X(2,2)*X(2,4)+X(1,1)*COS(X(2,1))*X(2,6)) !d2(y)/d(s1)d(s2)
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,6)-X(1,2)*SIN(X(2,1))*X(2,4)-&
                & (X(1,4)*SIN(X(2,1))*X(2,2)+X(1,1)*COS(X(2,1))*&
                & X(2,2)*X(2,4)+X(1,1)*SIN(X(2,1))*X(2,6)) !d2(x)/d(s1)d(s2)
              Z(2)=SIN(X(2,1))*X(1,6)+X(1,2)*COS(X(2,1))*X(2,4)+&
                & (X(1,4)*COS(X(2,1))*X(2,2)-X(1,1)*SIN(X(2,1))*&
                & X(2,2)*X(2,4)+X(1,1)*COS(X(2,1))*X(2,6)) !d2(y)/d(s1)d(s2)
              Z(3)=X(3,6) !d2(z)/d(s1)d(s2)
            CASE DEFAULT
              CALL FlagError("Invalid number of coordinates",err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Not enough X derivatives supplied",err,error,*999)
          ENDIF
        CASE(PART_DERIV_S3)
          IF(SIZE(X,2)>=7) THEN
            SELECT CASE(COORDINATE_SYSTEM%numberOfDimensions)
            CASE(2)
              Z(1)=0.0_DP
              Z(2)=0.0_DP
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,7)-X(1,1)*SIN(X(2,1))*X(2,7) !d(x)/d(s3)
              Z(2)=SIN(X(2,1))*X(1,7)+X(1,1)*COS(X(2,1))*X(2,7) !d(y)/d(s3)
              Z(3)=X(3,7) !d(z)/d(s3)
            CASE DEFAULT
              CALL FlagError("Invalid number of coordinates",err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Not enough X derivatives supplied",err,error,*999)
          ENDIF
        CASE(PART_DERIV_S3_S3)
          CALL FlagError("Not implemented",err,error,*999)
        CASE(PART_DERIV_S1_S3)
          IF(SIZE(X,2)>=9) THEN
            SELECT CASE(COORDINATE_SYSTEM%numberOfDimensions)
            CASE(2)
              Z(1)=0.0_DP
              Z(2)=0.0_DP
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,9)-X(1,2)*SIN(X(2,1))*X(2,7)-&
                & (X(1,7)*SIN(X(2,1))*X(2,2)+X(1,1)*COS(X(2,1))*&
                & X(2,2)*X(2,7)+X(1,1)*SIN(X(2,1))*X(2,9)) !d2(x)/d(s1)d(s3)
              Z(2)=SIN(X(2,1))*X(1,9)+X(1,2)*COS(X(2,1))*X(2,7)+&
                & (X(1,7)*COS(X(2,1))*X(2,2)-X(1,1)*SIN(X(2,1))*&
                & X(2,2)*X(2,7)+X(1,1)*COS(X(2,1))*X(2,9)) !d2(y)/d(s1)d(s3)
              Z(3)=X(3,9) !d2(z)/d(s1)d(s3)
            CASE DEFAULT
              CALL FlagError("Invalid number of coordinates",err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Not enough X derivatives supplied",err,error,*999)
          ENDIF
        CASE(PART_DERIV_S2_S3)
          IF(SIZE(X,2)>=10) THEN
            SELECT CASE(COORDINATE_SYSTEM%numberOfDimensions)
            CASE(2)
              Z(1)=0.0_DP
              Z(2)=0.0_DP
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,10)-X(1,4)*SIN(X(2,1))*X(2,7)-&
                & (X(1,7)*SIN(X(2,1))*X(2,4)+X(1,1)*COS(X(2,1))*&
                & X(2,4)*X(2,7)+X(1,1)*SIN(X(2,1))*X(2,10)) !d2(x)/d(s2)d(s3)
              Z(2)=SIN(X(2,1))*X(1,10)+X(1,4)*COS(X(2,1))*X(2,7)+&
                & (X(1,7)*COS(X(2,1))*X(2,4)-X(1,1)*SIN(X(2,1))*&
                & X(2,4)*X(2,7)+X(1,1)*COS(X(2,1))*X(2,10)) !d2(y)/d(s2)d(s3)
              Z(3)=X(3,10) !d2(z)/d(s2)d(s3)
            CASE DEFAULT
              CALL FlagError("Invalid number of coordinates",err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Not enough X derivatives supplied",err,error,*999)
          ENDIF
        CASE(PART_DERIV_S1_S2_S3)
          IF(SIZE(X,2)>=11) THEN
            SELECT CASE(COORDINATE_SYSTEM%numberOfDimensions)
            CASE(2)
              Z(1)=0.0_DP
              Z(2)=0.0_DP
            CASE(3)  
              Z(1)=-COS(X(2,1))*X(2,2)*X(2,4)*X(1,7)-&
                & SIN(X(2,1))*(X(2,6)*X(1,7)+X(2,4)*X(1,9))-&
                & SIN(X(2,1))*X(2,2)*X(1,10)+COS(X(2,1))*X(1,11)-&
                & COS(X(2,1))*X(2,2)*X(1,4)*X(2,7)-& 
                & SIN(X(2,1))*(X(1,6)*X(2,7)+X(1,4)*X(2,9))+&
                & (-COS(X(2,1))*X(1,2)+SIN(X(2,1))*X(1,1)*X(2,2))*&
                & X(2,4)*X(2,7)-X(1,1)*COS(X(2,1))*(X(2,6)*X(2,7)+&
                & X(2,4)*X(2,9))+(-SIN(X(2,1))*X(1,2)-X(1,1)*& 
                & COS(X(2,1))*X(2,2))*X(2,10)-X(1,1)*&
                & SIN(X(2,1))*X(2,11) !d3(x)/d(s1)d(s2)d(s3)
              Z(2)=-SIN(X(2,1))*X(2,2)*X(2,4)*X(1,7)+&
                & COS(X(2,1))*(X(2,6)*X(1,7)+X(2,4)*X(1,9))+&
                & COS(X(2,1))*X(2,2)*X(1,10)+SIN(X(2,1))*X(1,11)-&
                & SIN(X(2,1))*X(2,2)*X(1,4)*X(2,7)+&
                & COS(X(2,1))*(X(1,6)*X(2,7)+X(1,4)*X(2,9))+&
                & (-SIN(X(2,1))*X(1,2)-X(1,1)*COS(X(2,1))*X(2,2))*&
                & X(2,4)*X(2,7)-X(1,1)*SIN(X(2,1))*(X(2,6)*X(2,7)+&
                & X(2,4)*X(2,9))+(COS(X(2,1))*X(1,2)-X(1,1)*&
                & SIN(X(2,1))*X(2,2))*X(2,10)+X(1,1)*&
                & COS(X(2,1))*X(2,11) !d3(y)/d(s1)d(s2)d(s3)
              Z(3)=X(3,11) !d3(z)/d(s1)d(s2)d(s3)
            CASE DEFAULT
              CALL FlagError("Invalid number of coordinates",err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Not enough X derivatives supplied",err,error,*999)
          ENDIF
        CASE DEFAULT
          CALL FlagError("Invalid partial derivative type",err,error,*999)
        END SELECT
      CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
        IF(COORDINATE_SYSTEM%numberOfDimensions==3) THEN
          SELECT CASE(PART_DERIV_TYPE)
          CASE(NO_PART_DERIV)
            Z=COORDINATE_CONVERT_TO_RC(COORDINATE_SYSTEM,X(:,1),err,error)
            IF(err/=0) GOTO 999
          CASE(PART_DERIV_S1,PART_DERIV_S1_S1,PART_DERIV_S2,PART_DERIV_S2_S2,&
            & PART_DERIV_S1_S2,PART_DERIV_S3,PART_DERIV_S3_S3,&
            & PART_DERIV_S1_S3,PART_DERIV_S2_S3,PART_DERIV_S1_S2_S3)
            CALL FlagError("Not implemented",err,error,*999)
          CASE DEFAULT
            CALL FlagError("Invalid partial derivative type",err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Invalid number of coordinates",err,error,*999)
        ENDIF
      CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
        IF(COORDINATE_SYSTEM%numberOfDimensions==3) THEN
          FOCUS=COORDINATE_SYSTEM%FOCUS
          SELECT CASE(PART_DERIV_TYPE)
          CASE(NO_PART_DERIV)
            Z=COORDINATE_CONVERT_TO_RC(COORDINATE_SYSTEM,X(:,1),err,error)
            IF(err/=0) GOTO 999
          CASE(PART_DERIV_S1)
            IF(SIZE(X,2)>=2) THEN
              Z(1)=FOCUS*SINH(X(1,1))*COS(X(2,1))*X(1,2)-&
                & FOCUS*COSH(X(1,1))*SIN(X(2,1))*X(2,2) !d(x)/d(s1)
              Z(2)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,2) !d(y)/d(s1)
              Z(3)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,2) !d(z)/d(s1)
            ELSE
              CALL FlagError("Invalid derivative type",err,error,*999)
            ENDIF
          CASE(PART_DERIV_S1_S1)
            CALL FlagError("Not implemented",err,error,*999)
          CASE(PART_DERIV_S2)
            IF(SIZE(X,2)>=4) THEN
              Z(1)=FOCUS*SINH(X(1,1))*COS(X(2,1))*X(1,4)-&
                & FOCUS*COSH(X(1,1))*SIN(X(2,1))*X(2,4) !d(x)/d(s2)
              Z(2)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,4)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,4)-&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,4) !d(y)/d(s2)
              Z(3)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,4)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,4)+&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,4) !d(z)/d(s2)
            ELSE
              CALL FlagError("Not enough X derivatives supplied",&
                & err,error,*999)
            ENDIF
          CASE(PART_DERIV_S2_S2)
            CALL FlagError("Not implemented",err,error,*999)
          CASE(PART_DERIV_S1_S2)
            IF(SIZE(X,2)>=6) THEN
              Z(1)=FOCUS*(SINH(X(1,1))*COS(X(2,1))*X(1,6)+&
                & X(1,2)*(COSH(X(1,1))*COS(X(2,1))*X(1,4)-&
                & SINH(X(1,1))*SIN(X(2,1))*X(2,4))-&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,6)-&
                & X(2,2)*(SINH(X(1,1))*SIN(X(2,1))*X(1,4)+&
                & COSH(X(1,1))*COS(X(2,1))*X(2,4))) !d2(x)/d(s1)d(s2)
              Z(2)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,6)+&
                & X(1,2)*(SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(1,4)+COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & X(2,4)-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,4))+SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & X(2,6)+X(2,2)*(COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & X(1,4)-SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,4)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,4))-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,6)-&
                & X(3,2)*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(1,4)+SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & X(2,4)+SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,4))) !d2(y)/d(s1)d(s2)
              Z(3)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,6)+&
                & X(1,2)*(SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(1,4)+COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & X(2,4)+COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,4))+SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,6)+&
                & X(2,2)*(COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & X(1,4)-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(2,4)+SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & X(3,4))+SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,6)+&
                & X(3,2)*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(1,4)+SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & X(2,4)-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,4))) !d2(z)/d(s1)d(s2)
            ELSE
              CALL FlagError("Not enough X derivatives supplied",&
                & err,error,*999)
            ENDIF
          CASE(PART_DERIV_S3)
            IF(SIZE(X,2)>=7) THEN
              Z(1)=FOCUS*SINH(X(1,1))*COS(X(2,1))*X(1,7)-&
                & FOCUS*COSH(X(1,1))*SIN(X(2,1))*X(2,7) !d(x)/d(s3)
              Z(2)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,7)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,7)-&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,7) !d(y)/d(s3)
              Z(3)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,7)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,7) !d(z)/d(s3)
            ELSE
              CALL FlagError("Not enough X derivatives supplied",&
                & err,error,*999)
            ENDIF
          CASE(PART_DERIV_S3_S3)
            CALL FlagError("Not implemented",err,error,*999)
          CASE(PART_DERIV_S1_S3)
            IF(SIZE(X,2)>=9) THEN
              Z(1)=FOCUS*(SINH(X(1,1))*COS(X(2,1))*X(1,9)+&
                & X(1,2)*(COSH(X(1,1))*COS(X(2,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*X(2,7))-&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,9)-&
                & X(2,2)*(SINH(X(1,1))*SIN(X(2,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*X(2,7))) !d2(x)/d(s1)d(s3)
              Z(2)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,9)+&
                & X(1,2)*(SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,7)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,9)+&
                & X(2,2)*(COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,7)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,7))-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,9)-&
                & X(3,2)*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,7)+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))* X(2,7)+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,7))) !d2(y)/d(s1)d(s3) 
              Z(3)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,9)+&
                & X(1,2)*(SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,9)+&
                & X(2,2)*(COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,9)+&
                & X(3,2)*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,7)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,7))) !d2(z)/d(s1)d(s3)
            ELSE
              CALL FlagError("Not enough X derivatives supplied",&
                & err,error,*999)
            ENDIF
          CASE(PART_DERIV_S2_S3)
            IF(SIZE(X,2)>=10) THEN
              Z(1)=FOCUS*(SINH(X(1,1))*COS(X(2,1))*X(1,10)+&
                & X(1,4)*(COSH(X(1,1))*COS(X(2,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*X(2,7))-&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,10)-&
                & X(2,4)*(SINH(X(1,1))*SIN(X(2,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*X(2,7))) !d2(x)/d(s2)d(s3)
              Z(2)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,10)+&
                & X(1,4)*(SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,7)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,10)+&
                & X(2,4)*(COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,7)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,7))-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,10)-&
                & X(3,4)*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,7)+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,7))) !d2(y)/d(s2)d(s3)
              Z(3)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,10)+&
                & X(1,4)*(SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,10)+&
                & X(2,4)*(COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,10)+&
                & X(3,4)*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,7)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,7))) !d2(z)/d(s2)d(s3)
            ELSE
              CALL FlagError("Not enough X derivatives supplied",&
                & err,error,*999)
            ENDIF
          CASE(PART_DERIV_S1_S2_S3)
            IF(SIZE(X,2)>=11) THEN
              Z(1)=FOCUS*((SINH(X(1,1))*COS(X(2,1))*X(1,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,2))*X(1,4)*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*(X(1,6)*X(1,7)+X(1,4)*X(1,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*X(2,2))*X(2,4)*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*(X(2,6)*X(1,7)+X(2,4)*X(1,9))+&
                & (COSH(X(1,1))*COS(X(2,1))*X(1,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*X(2,2))*X(1,10)+&
                & SINH(X(1,1))*COS(X(2,1))*X(1,11)+&
                & (-COSH(X(1,1))*SIN(X(2,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*X(2,2))*X(1,4)*X(2,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*(X(1,6)*X(2,7)+X(1,4)*X(2,9))+&
                & (-SINH(X(1,1))*COS(X(2,1))*X(1,2)+&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,2))*X(2,4)*X(2,7)-&
                & COSH(X(1,1))*COS(X(2,1))*(X(2,6)*X(2,7)+X(2,4)*X(2,9))+&
                & (-SINH(X(1,1))*SIN(X(2,1))*X(1,2)-&
                & COSH(X(1,1))*COS(X(2,1))*X(2,2))*X(2,10)-&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,11)) !d3(x)/d(s1)d(s2)d(s3)
              Z(2)=FOCUS*((COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(1,4)*X(1,7)+SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & (X(1,6)*X(1,7)+X(1,4)*X(1,9))+&
                & (SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(2,4)*X(1,7)+COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & (X(2,6)*X(1,7)+X(2,4)*X(1,9))+&
                & (-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(3,4)*X(1,7)-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & (X(3,6)*X(1,7)+X(3,4)*X(1,9))+&
                & (SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*X(1,10)+&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,11))+&
                & FOCUS*((SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(1,4)*X(2,7)+COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & (X(1,6)*X(2,7)+X(1,4)*X(2,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)+&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(2,4)*X(2,7)-SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & (X(2,6)*X(2,7)+X(2,4)*X(2,9))+&
                & (-COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,2)+&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(3,4)*X(2,7)-SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & (X(3,6)*X(2,7)+X(3,4)*X(2,9))+&
                & (COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,2))*X(2,10)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,11))+&
                & FOCUS*((-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(1,4)*X(3,7)-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & (X(1,6)*X(3,7)+X(1,4)*X(3,9))+&
                & (-COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,2)+&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(2,4)*X(3,7)-SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & (X(2,6)*X(3,7)+X(2,4)*X(3,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)+&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(3,4)*X(3,7)-SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & (X(3,6)*X(3,7)+X(3,4)*X(3,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*X(3,10)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,11)) !d3(y)/d(s1)d(s2)d(s3)
              Z(3)=FOCUS*((COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(1,4)*X(1,7)+SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & (X(1,6)*X(1,7)+X(1,4)*X(1,9))+&
                & (SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(2,4)*X(1,7)+COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & (X(2,6)*X(1,7)+X(2,4)*X(1,9))+&
                & (SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(3,4)*X(1,7)+COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & (X(3,6)*X(1,7)+X(3,4)*X(1,9))+&
                & (SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)+&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*X(1,10)+&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,11))+&
                & FOCUS*((SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(1,4)*X(2,7)+COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & (X(1,6)*X(2,7)+X(1,4)*X(2,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(2,4)*X(2,7)-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & (X(2,6)*X(2,7)+X(2,4)*X(2,9))+&
                & (COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(3,4)*X(2,7)+SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & (X(3,6)*X(2,7)+X(3,4)*X(2,9))+&
                & (COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,2))*X(2,10)+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,11))+&
                & FOCUS*((SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(1,4)*X(3,7)+COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & (X(1,6)*X(3,7)+X(1,4)*X(3,9))+&
                & (COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(2,4)*X(3,7)+SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & (X(2,6)*X(3,7)+X(2,4)*X(3,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(3,4)*X(3,7)-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & (X(3,6)*X(3,7)+X(3,4)*X(3,9))+&
                & (COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*X(3,10)+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,11)) !d3(z)/d(s1)d(s2)d(s3)
            ELSE
              CALL FlagError("Not enough X derivatives supplied",&
                & err,error,*999)
            ENDIF
          CASE DEFAULT
            CALL FlagError("Invalid partial derivative type",err,error,*999)
          END SELECT
        ENDIF
      CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
        IF(COORDINATE_SYSTEM%numberOfDimensions==3) THEN
          FOCUS=COORDINATE_SYSTEM%FOCUS
          SELECT CASE(PART_DERIV_TYPE)
          CASE(NO_PART_DERIV)
            Z=COORDINATE_CONVERT_TO_RC(COORDINATE_SYSTEM,X(:,1),err,error)
            IF(err/=0) GOTO 999
          CASE(PART_DERIV_S1,PART_DERIV_S1_S1,PART_DERIV_S2,PART_DERIV_S2_S2,&
            & PART_DERIV_S1_S2,PART_DERIV_S3,PART_DERIV_S3_S3,&
            & PART_DERIV_S1_S3,PART_DERIV_S2_S3,PART_DERIV_S1_S2_S3)
            CALL FlagError("Not implemented",err,error,*999)
          CASE DEFAULT
            CALL FlagError("Invalid partial derivative type",err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Invalid number of coordinates",err,error,*999)
        ENDIF
      CASE DEFAULT
        CALL FlagError("Invalid coordinate type",err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("The sizes of the vectors X and Z do not match",&
        & err,error,*999)
    ENDIF

    EXITS("COORDINATE_DERIVATIVE_CONVERT_TO_RC_DP")
    RETURN
999 ERRORSEXITS("COORDINATE_DERIVATIVE_CONVERT_TO_RC_DP",err,error)
    RETURN 1
  END SUBROUTINE COORDINATE_DERIVATIVE_CONVERT_TO_RC_DP

  !
  !================================================================================================================================
  !
  
  SUBROUTINE COORDINATE_DERIVATIVE_CONVERT_TO_RC_SP(COORDINATE_SYSTEM,PART_DERIV_TYPE,X,Z,&
    & err,error,*)
  
    !#### Subroutine: COORDINATE_DERIVATIVE_CONVERT_TO_RC_SP
    !###  Description:
    !###    COORDINATE_DERIVATIVE_CONVERT_TO_RC_SP performs a coordinate transformation from a
    !###    coordinate system identified by COORDINATE_SYSTEM at the point X
    !###    with coordinates/derivatives X(nj,nu) to the point with 
    !###    coordinates/derivatives Z(nj) in rectangular cartesian coordinates
    !###    for single precision coordinates.
    !###  Parent-function: COORDINATE_DERIVATIVE_CONVERT_TO_RC
    
    !Argument variables
    TYPE(CoordinateSystemType), INTENT(IN) :: COORDINATE_SYSTEM
    INTEGER(INTG), INTENT(IN) :: PART_DERIV_TYPE
    REAL(SP), INTENT(IN) :: X(:,:)
    REAL(SP), INTENT(OUT) :: Z(:)
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local variables
    REAL(SP) :: FOCUS
    
    ENTERS("COORDINATE_DERIVATIVE_CONVERT_TO_RC_SP",err,error,*999)
    
!!TODO: change all second index X(:,?) numbers to their apropriate constant
!!as defined in constants e.g. X(1,2) == X(1,PART_DERIV_S1)

    IF(SIZE(X,1)<COORDINATE_SYSTEM%numberOfDimensions) &
      & CALL FlagError("Size of X is less than the number of dimensions", &
      & err,error,*999)
    
    IF(SIZE(X,1)==SIZE(Z,1)) THEN
      SELECT CASE(COORDINATE_SYSTEM%TYPE)
      CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
        IF(SIZE(X,2)>=PART_DERIV_TYPE) THEN
          Z=X(:,PART_DERIV_TYPE)
        ELSE
          CALL FlagError("Invalid partial derivative type",err,error,*999)
        ENDIF
      CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
        SELECT CASE(PART_DERIV_TYPE)
        CASE(NO_PART_DERIV)
          Z=COORDINATE_CONVERT_TO_RC(COORDINATE_SYSTEM,X(:,1),err,error)
          IF(err/=0) GOTO 999
        CASE(PART_DERIV_S1)
          IF(SIZE(X,2)>=2) THEN
            SELECT CASE(COORDINATE_SYSTEM%numberOfDimensions)
            CASE(2)
              Z(1)=COS(X(2,1))*X(1,2)-X(1,1)*SIN(X(2,1))*X(2,2) !d(x)/d(s1)
              Z(2)=SIN(X(2,1))*X(1,2)+X(1,1)*COS(X(2,1))*X(2,2) !d(y)/d(s1)
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,2)-X(1,1)*SIN(X(2,1))*X(2,2) !d(x)/d(s1)
              Z(2)=SIN(X(2,1))*X(1,2)+X(1,1)*COS(X(2,1))*X(2,2) !d(y)/d(s1)
              Z(3)=X(3,2) !d(z)/d(s1)
            CASE DEFAULT
              CALL FlagError("Invalid number of coordinates",err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Not enough X derivatives supplied",err,error,*999)
          ENDIF
        CASE(PART_DERIV_S1_S1)
          CALL FlagError("Not implemented",err,error,*999)
        CASE(PART_DERIV_S2)
          IF(SIZE(X,2)>=4) THEN
            SELECT CASE(COORDINATE_SYSTEM%numberOfDimensions)
            CASE(2)
              Z(1)=COS(X(2,1))*X(1,4)-X(1,1)*SIN(X(2,1))*X(2,4) !d(x)/d(s2)
              Z(2)=SIN(X(2,1))*X(1,4)+X(1,1)*COS(X(2,1))*X(2,4) !d(y)/d(s2)
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,4)-X(1,1)*SIN(X(2,1))*X(2,4) !d(x)/d(s2)
              Z(2)=SIN(X(2,1))*X(1,4)+X(1,1)*COS(X(2,1))*X(2,4) !d(y)/d(s2)
              Z(3)=X(3,4) !d(z)/d(s2)
            CASE DEFAULT
              CALL FlagError("Invalid number of coordinates",err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Not enough X derivatives supplied",err,error,*999)
          ENDIF
        CASE(PART_DERIV_S2_S2)
          CALL FlagError("Not implemented",err,error,*999)
        CASE(PART_DERIV_S1_S2)
          IF(SIZE(X,2)>=6) THEN
            SELECT CASE(COORDINATE_SYSTEM%numberOfDimensions)
            CASE(2)
              Z(1)=COS(X(2,1))*X(1,6)-X(1,2)*SIN(X(2,1))*X(2,4)-&
                & (X(1,4)*SIN(X(2,1))*X(2,2)+X(1,1)*COS(X(2,1))*&
                & X(2,2)*X(2,4)+X(1,1)*SIN(X(2,1))*X(2,6)) !d2(x)/d(s1)d(s2)
              Z(2)=SIN(X(2,1))*X(1,6)+X(1,2)*COS(X(2,1))*X(2,4)+&
                & (X(1,4)*COS(X(2,1))*X(2,2)-X(1,1)*SIN(X(2,1))*&
                & X(2,2)*X(2,4)+X(1,1)*COS(X(2,1))*X(2,6)) !d2(y)/d(s1)d(s2)
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,6)-X(1,2)*SIN(X(2,1))*X(2,4)-&
                & (X(1,4)*SIN(X(2,1))*X(2,2)+X(1,1)*COS(X(2,1))*&
                & X(2,2)*X(2,4)+X(1,1)*SIN(X(2,1))*X(2,6)) !d2(x)/d(s1)d(s2)
              Z(2)=SIN(X(2,1))*X(1,6)+X(1,2)*COS(X(2,1))*X(2,4)+&
                & (X(1,4)*COS(X(2,1))*X(2,2)-X(1,1)*SIN(X(2,1))*&
                & X(2,2)*X(2,4)+X(1,1)*COS(X(2,1))*X(2,6)) !d2(y)/d(s1)d(s2)
              Z(3)=X(3,6) !d2(z)/d(s1)d(s2)
            CASE DEFAULT
              CALL FlagError("Invalid number of coordinates",err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Not enough X derivatives supplied",err,error,*999)
          ENDIF
        CASE(PART_DERIV_S3)
          IF(SIZE(X,2)>=7) THEN
            SELECT CASE(COORDINATE_SYSTEM%numberOfDimensions)
            CASE(2)
              Z(1)=0.0_SP
              Z(2)=0.0_SP
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,7)-X(1,1)*SIN(X(2,1))*X(2,7) !d(x)/d(s3)
              Z(2)=SIN(X(2,1))*X(1,7)+X(1,1)*COS(X(2,1))*X(2,7) !d(y)/d(s3)
              Z(3)=X(3,7) !d(z)/d(s3)
            CASE DEFAULT
              CALL FlagError("Invalid number of coordinates",err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Not enough X derivatives supplied",err,error,*999)
          ENDIF
        CASE(PART_DERIV_S3_S3)
          CALL FlagError("Not implemented",err,error,*999)
        CASE(PART_DERIV_S1_S3)
          IF(SIZE(X,2)>=9) THEN
            SELECT CASE(COORDINATE_SYSTEM%numberOfDimensions)
            CASE(2)
              Z(1)=0.0_SP
              Z(2)=0.0_SP
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,9)-X(1,2)*SIN(X(2,1))*X(2,7)-&
                & (X(1,7)*SIN(X(2,1))*X(2,2)+X(1,1)*COS(X(2,1))*&
                & X(2,2)*X(2,7)+X(1,1)*SIN(X(2,1))*X(2,9)) !d2(x)/d(s1)d(s3)
              Z(2)=SIN(X(2,1))*X(1,9)+X(1,2)*COS(X(2,1))*X(2,7)+&
                & (X(1,7)*COS(X(2,1))*X(2,2)-X(1,1)*SIN(X(2,1))*&
                & X(2,2)*X(2,7)+X(1,1)*COS(X(2,1))*X(2,9)) !d2(y)/d(s1)d(s3)
              Z(3)=X(3,9) !d2(z)/d(s1)d(s3)
            CASE DEFAULT
              CALL FlagError("Invalid number of coordinates",err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Not enough X derivatives supplied",err,error,*999)
          ENDIF
        CASE(PART_DERIV_S2_S3)
          IF(SIZE(X,2)>=10) THEN
            SELECT CASE(COORDINATE_SYSTEM%numberOfDimensions)
            CASE(2)
              Z(1)=0.0_SP
              Z(2)=0.0_SP
            CASE(3)
              Z(1)=COS(X(2,1))*X(1,10)-X(1,4)*SIN(X(2,1))*X(2,7)-&
                & (X(1,7)*SIN(X(2,1))*X(2,4)+X(1,1)*COS(X(2,1))*&
                & X(2,4)*X(2,7)+X(1,1)*SIN(X(2,1))*X(2,10)) !d2(x)/d(s2)d(s3)
              Z(2)=SIN(X(2,1))*X(1,10)+X(1,4)*COS(X(2,1))*X(2,7)+&
                & (X(1,7)*COS(X(2,1))*X(2,4)-X(1,1)*SIN(X(2,1))*&
                & X(2,4)*X(2,7)+X(1,1)*COS(X(2,1))*X(2,10)) !d2(y)/d(s2)d(s3)
              Z(3)=X(3,10) !d2(z)/d(s2)d(s3)
            CASE DEFAULT
              CALL FlagError("Invalid number of coordinates",err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Not enough X derivatives supplied",err,error,*999)
          ENDIF
        CASE(PART_DERIV_S1_S2_S3)
          IF(SIZE(X,2)>=11) THEN
            SELECT CASE(COORDINATE_SYSTEM%numberOfDimensions)
            CASE(2)
              Z(1)=0.0_SP
              Z(2)=0.0_SP
            CASE(3)  
              Z(1)=-COS(X(2,1))*X(2,2)*X(2,4)*X(1,7)-&
                & SIN(X(2,1))*(X(2,6)*X(1,7)+X(2,4)*X(1,9))-&
                & SIN(X(2,1))*X(2,2)*X(1,10)+COS(X(2,1))*X(1,11)-&
                & COS(X(2,1))*X(2,2)*X(1,4)*X(2,7)-& 
                & SIN(X(2,1))*(X(1,6)*X(2,7)+X(1,4)*X(2,9))+&
                & (-COS(X(2,1))*X(1,2)+SIN(X(2,1))*X(1,1)*X(2,2))*&
                & X(2,4)*X(2,7)-X(1,1)*COS(X(2,1))*(X(2,6)*X(2,7)+&
                & X(2,4)*X(2,9))+(-SIN(X(2,1))*X(1,2)-X(1,1)*& 
                & COS(X(2,1))*X(2,2))*X(2,10)-X(1,1)*&
                & SIN(X(2,1))*X(2,11) !d3(x)/d(s1)d(s2)d(s3)
              Z(2)=-SIN(X(2,1))*X(2,2)*X(2,4)*X(1,7)+&
                & COS(X(2,1))*(X(2,6)*X(1,7)+X(2,4)*X(1,9))+&
                & COS(X(2,1))*X(2,2)*X(1,10)+SIN(X(2,1))*X(1,11)-&
                & SIN(X(2,1))*X(2,2)*X(1,4)*X(2,7)+&
                & COS(X(2,1))*(X(1,6)*X(2,7)+X(1,4)*X(2,9))+&
                & (-SIN(X(2,1))*X(1,2)-X(1,1)*COS(X(2,1))*X(2,2))*&
                & X(2,4)*X(2,7)-X(1,1)*SIN(X(2,1))*(X(2,6)*X(2,7)+&
                & X(2,4)*X(2,9))+(COS(X(2,1))*X(1,2)-X(1,1)*&
                & SIN(X(2,1))*X(2,2))*X(2,10)+X(1,1)*&
                & COS(X(2,1))*X(2,11) !d3(y)/d(s1)d(s2)d(s3)
              Z(3)=X(3,11) !d3(z)/d(s1)d(s2)d(s3)
            CASE DEFAULT
              CALL FlagError("Invalid number of coordinates",err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Not enough X derivatives supplied",err,error,*999)
          ENDIF
        CASE DEFAULT
          CALL FlagError("Invalid partial derivative type",err,error,*999)
        END SELECT
      CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
        IF(COORDINATE_SYSTEM%numberOfDimensions==3) THEN
          SELECT CASE(PART_DERIV_TYPE)
          CASE(NO_PART_DERIV)
            Z=COORDINATE_CONVERT_TO_RC(COORDINATE_SYSTEM,X(:,1),err,error)
            IF(err/=0) GOTO 999
          CASE(PART_DERIV_S1,PART_DERIV_S1_S1,PART_DERIV_S2,PART_DERIV_S2_S2,&
            & PART_DERIV_S1_S2,PART_DERIV_S3,PART_DERIV_S3_S3,&
            & PART_DERIV_S1_S3,PART_DERIV_S2_S3,PART_DERIV_S1_S2_S3)
            CALL FlagError("Not implemented",err,error,*999)
          CASE DEFAULT
            CALL FlagError("Invalid partial derivative type",err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Invalid number of coordinates",err,error,*999)
        ENDIF
      CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
        IF(COORDINATE_SYSTEM%numberOfDimensions==3) THEN
          FOCUS=REAL(COORDINATE_SYSTEM%FOCUS,SP)
          SELECT CASE(PART_DERIV_TYPE)
          CASE(NO_PART_DERIV)
            Z=COORDINATE_CONVERT_TO_RC(COORDINATE_SYSTEM,X(:,1),err,error)
            IF(err/=0) GOTO 999
          CASE(PART_DERIV_S1)
            IF(SIZE(X,2)>=2) THEN
              Z(1)=FOCUS*SINH(X(1,1))*COS(X(2,1))*X(1,2)-&
                & FOCUS*COSH(X(1,1))*SIN(X(2,1))*X(2,2) !d(x)/d(s1)
              Z(2)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,2) !d(y)/d(s1)
              Z(3)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,2) !d(z)/d(s1)
            ELSE
              CALL FlagError("Not enough X derivatives supplied", &
                & err,error,*999)
            ENDIF
          CASE(PART_DERIV_S1_S1)
            CALL FlagError("Not implemented",err,error,*999)
          CASE(PART_DERIV_S2)
            IF(SIZE(X,2)>=4) THEN
              Z(1)=FOCUS*SINH(X(1,1))*COS(X(2,1))*X(1,4)-&
                & FOCUS*COSH(X(1,1))*SIN(X(2,1))*X(2,4) !d(x)/d(s2)
              Z(2)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,4)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,4)-&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,4) !d(y)/d(s2)
              Z(3)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,4)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,4)+&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,4) !d(z)/d(s2)
            ELSE
              CALL FlagError("Not enough X derivatives supplied",&
                & err,error,*999)
            ENDIF
          CASE(PART_DERIV_S2_S2)
            CALL FlagError("Not implemented",err,error,*999)
          CASE(PART_DERIV_S1_S2)
            IF(SIZE(X,2)>=6) THEN
              Z(1)=FOCUS*(SINH(X(1,1))*COS(X(2,1))*X(1,6)+&
                & X(1,2)*(COSH(X(1,1))*COS(X(2,1))*X(1,4)-&
                & SINH(X(1,1))*SIN(X(2,1))*X(2,4))-&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,6)-&
                & X(2,2)*(SINH(X(1,1))*SIN(X(2,1))*X(1,4)+&
                & COSH(X(1,1))*COS(X(2,1))*X(2,4))) !d2(x)/d(s1)d(s2)
              Z(2)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,6)+&
                & X(1,2)*(SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(1,4)+COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & X(2,4)-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,4))+SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & X(2,6)+X(2,2)*(COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & X(1,4)-SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,4)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,4))-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,6)-&
                & X(3,2)*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(1,4)+SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & X(2,4)+SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,4))) !d2(y)/d(s1)d(s2)
              Z(3)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,6)+&
                & X(1,2)*(SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(1,4)+COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & X(2,4)+COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,4))+SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,6)+&
                & X(2,2)*(COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & X(1,4)-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(2,4)+SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & X(3,4))+SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,6)+&
                & X(3,2)*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(1,4)+SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & X(2,4)-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,4))) !d2(z)/d(s1)d(s2)
            ELSE
              CALL FlagError("Not enough X derivatives supplied",&
                & err,error,*999)
            ENDIF
          CASE(PART_DERIV_S3)
            IF(SIZE(X,2)>=7) THEN
              Z(1)=FOCUS*SINH(X(1,1))*COS(X(2,1))*X(1,7)-&
                & FOCUS*COSH(X(1,1))*SIN(X(2,1))*X(2,7) !d(x)/d(s3)
              Z(2)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,7)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,7)-&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,7) !d(y)/d(s3)
              Z(3)=FOCUS*COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,7)+&
                & FOCUS*SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & FOCUS*SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,7) !d(z)/d(s3)
            ELSE
              CALL FlagError("Not enough X derivatives supplied",&
                & err,error,*999)
            ENDIF
          CASE(PART_DERIV_S3_S3)
            CALL FlagError("Not implemented",err,error,*999)
          CASE(PART_DERIV_S1_S3)
            IF(SIZE(X,2)>=9) THEN
              Z(1)=FOCUS*(SINH(X(1,1))*COS(X(2,1))*X(1,9)+&
                & X(1,2)*(COSH(X(1,1))*COS(X(2,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*X(2,7))-&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,9)-&
                & X(2,2)*(SINH(X(1,1))*SIN(X(2,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*X(2,7))) !d2(x)/d(s1)d(s3)
              Z(2)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,9)+&
                & X(1,2)*(SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,7)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,9)+&
                & X(2,2)*(COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,7)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,7))-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,9)-&
                & X(3,2)*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,7)+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))* X(2,7)+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,7))) !d2(y)/d(s1)d(s3) 
              Z(3)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,9)+&
                & X(1,2)*(SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,9)+&
                & X(2,2)*(COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,9)+&
                & X(3,2)*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,7)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,7))) !d2(z)/d(s1)d(s3)
            ELSE
              CALL FlagError("Not enough X derivatives supplied",&
                & err,error,*999)
            ENDIF
          CASE(PART_DERIV_S2_S3)
            IF(SIZE(X,2)>=10) THEN
              Z(1)=FOCUS*(SINH(X(1,1))*COS(X(2,1))*X(1,10)+&
                & X(1,4)*(COSH(X(1,1))*COS(X(2,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*X(2,7))-&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,10)-&
                & X(2,4)*(SINH(X(1,1))*SIN(X(2,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*X(2,7))) !d2(x)/d(s2)d(s3)
              Z(2)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,10)+&
                & X(1,4)*(SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,7)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,10)+&
                & X(2,4)*(COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,7)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,7))-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,10)-&
                & X(3,4)*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,7)+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,7))) !d2(y)/d(s2)d(s3)
              Z(3)=FOCUS*(COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,10)+&
                & X(1,4)*(SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,10)+&
                & X(2,4)*(COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,7)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,7))+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,10)+&
                & X(3,4)*(COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,7)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,7))) !d2(z)/d(s2)d(s3)
            ELSE
              CALL FlagError("Not enough X derivatives supplied",&
                & err,error,*999)
            ENDIF
          CASE(PART_DERIV_S1_S2_S3)
            IF(SIZE(X,2)>=11) THEN
              Z(1)=FOCUS*((SINH(X(1,1))*COS(X(2,1))*X(1,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,2))*X(1,4)*X(1,7)+&
                & COSH(X(1,1))*COS(X(2,1))*(X(1,6)*X(1,7)+X(1,4)*X(1,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*X(2,2))*X(2,4)*X(1,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*(X(2,6)*X(1,7)+X(2,4)*X(1,9))+&
                & (COSH(X(1,1))*COS(X(2,1))*X(1,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*X(2,2))*X(1,10)+&
                & SINH(X(1,1))*COS(X(2,1))*X(1,11)+&
                & (-COSH(X(1,1))*SIN(X(2,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*X(2,2))*X(1,4)*X(2,7)-&
                & SINH(X(1,1))*SIN(X(2,1))*(X(1,6)*X(2,7)+X(1,4)*X(2,9))+&
                & (-SINH(X(1,1))*COS(X(2,1))*X(1,2)+&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,2))*X(2,4)*X(2,7)-&
                & COSH(X(1,1))*COS(X(2,1))*(X(2,6)*X(2,7)+X(2,4)*X(2,9))+&
                & (-SINH(X(1,1))*SIN(X(2,1))*X(1,2)-&
                & COSH(X(1,1))*COS(X(2,1))*X(2,2))*X(2,10)-&
                & COSH(X(1,1))*SIN(X(2,1))*X(2,11)) !d3(x)/d(s1)d(s2)d(s3)
              Z(2)=FOCUS*((COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(1,4)*X(1,7)+SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & (X(1,6)*X(1,7)+X(1,4)*X(1,9))+&
                & (SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(2,4)*X(1,7)+COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & (X(2,6)*X(1,7)+X(2,4)*X(1,9))+&
                & (-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(3,4)*X(1,7)-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & (X(3,6)*X(1,7)+X(3,4)*X(1,9))+&
                & (SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*X(1,10)+&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,11))+&
                & FOCUS*((SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(1,4)*X(2,7)+COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & (X(1,6)*X(2,7)+X(1,4)*X(2,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)+&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(2,4)*X(2,7)-SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & (X(2,6)*X(2,7)+X(2,4)*X(2,9))+&
                & (-COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,2)+&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(3,4)*X(2,7)-SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & (X(3,6)*X(2,7)+X(3,4)*X(2,9))+&
                & (COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,2))*X(2,10)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,11))+&
                & FOCUS*((-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(1,4)*X(3,7)-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & (X(1,6)*X(3,7)+X(1,4)*X(3,9))+&
                & (-COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,2)+&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(2,4)*X(3,7)-SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & (X(2,6)*X(3,7)+X(2,4)*X(3,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)+&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(3,4)*X(3,7)-SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & (X(3,6)*X(3,7)+X(3,4)*X(3,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*X(3,10)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & X(3,11)) !d3(y)/d(s1)d(s2)d(s3)
              Z(3)=FOCUS*((COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(1,4)*X(1,7)+SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & (X(1,6)*X(1,7)+X(1,4)*X(1,9))+&
                & (SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(2,4)*X(1,7)+COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & (X(2,6)*X(1,7)+X(2,4)*X(1,9))+&
                & (SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(3,4)*X(1,7)+COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & (X(3,6)*X(1,7)+X(3,4)*X(1,9))+&
                & (SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)+&
                & COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*X(1,10)+&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,11))+&
                & FOCUS*((SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(1,4)*X(2,7)+COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*&
                & (X(1,6)*X(2,7)+X(1,4)*X(2,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(2,4)*X(2,7)-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & (X(2,6)*X(2,7)+X(2,4)*X(2,9))+&
                & (COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(3,4)*X(2,7)+SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & (X(3,6)*X(2,7)+X(3,4)*X(2,9))+&
                & (COSH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(2,2)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(3,2))*X(2,10)+&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,11))+&
                & FOCUS*((SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(1,4)*X(3,7)+COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & (X(1,6)*X(3,7)+X(1,4)*X(3,9))+&
                & (COSH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(3,2))*&
                & X(2,4)*X(3,7)+SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*&
                & (X(2,6)*X(3,7)+X(2,4)*X(3,9))+&
                & (-COSH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(1,2)-&
                & SINH(X(1,1))*COS(X(2,1))*SIN(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(3,2))*&
                & X(3,4)*X(3,7)-SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*&
                & (X(3,6)*X(3,7)+X(3,4)*X(3,9))+&
                & (COSH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*X(1,2)+&
                & SINH(X(1,1))*COS(X(2,1))*COS(X(3,1))*X(2,2)-&
                & SINH(X(1,1))*SIN(X(2,1))*SIN(X(3,1))*X(3,2))*X(3,10)+&
                & SINH(X(1,1))*SIN(X(2,1))*COS(X(3,1))*&
                & X(3,11)) !d3(z)/d(s1)d(s2)d(s3)
            ELSE
              CALL FlagError("Not enough X derivatives supplied",&
                & err,error,*999)
            ENDIF
          CASE DEFAULT
            CALL FlagError("Invalid partial derivative type",err,error,*999)
          END SELECT
        ENDIF
      CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
        IF(COORDINATE_SYSTEM%numberOfDimensions==3) THEN
          SELECT CASE(PART_DERIV_TYPE)
          CASE(NO_PART_DERIV)
            Z=COORDINATE_CONVERT_TO_RC(COORDINATE_SYSTEM,X(:,1),err,error)
            IF(err/=0) GOTO 999
          CASE(PART_DERIV_S1,PART_DERIV_S1_S1,PART_DERIV_S2,PART_DERIV_S2_S2,&
            & PART_DERIV_S1_S2,PART_DERIV_S3,PART_DERIV_S3_S3,&
            & PART_DERIV_S1_S3,PART_DERIV_S2_S3,PART_DERIV_S1_S2_S3)
            CALL FlagError("Not implemented",err,error,*999)
          CASE DEFAULT
            CALL FlagError("Invalid partial derivative type",err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Invalid number of coordinates",err,error,*999)
        ENDIF
      CASE DEFAULT
        CALL FlagError("Invalid coordinate type",err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("The sizes of the vectors X and Z do not match",&
        & err,error,*999)
    ENDIF

    EXITS("COORDINATE_DERIVATIVE_CONVERT_TO_RC_SP")
    RETURN
999 ERRORSEXITS("COORDINATE_DERIVATIVE_CONVERT_TO_RC_SP",err,error)
    RETURN 1
  END SUBROUTINE COORDINATE_DERIVATIVE_CONVERT_TO_RC_SP

  !
  !================================================================================================================================
  !

  !>Calculates the norm of a derivative in a coordinate system identified by COORDINATE_SYSTEM at the given interpolated
  !>point and returns the value in NORM for single precision coordinates. PART_DERIV_INDEX is used to select the
  !>appropriate partial derivative (i.e., wrt S1, S2 or S3) to normalise.
  SUBROUTINE COORDINATE_DERIVATIVE_NORM(COORDINATE_SYSTEM,PART_DERIV_INDEX,INTERPOLATED_POINT,DERIV_NORM,err,error,*)
  
    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to calculate the derivative norm for
    INTEGER(INTG), INTENT(IN) :: PART_DERIV_INDEX !<The partial derivative index to select the direction to normalise
    TYPE(FieldInterpolatedPointType), POINTER :: INTERPOLATED_POINT !<A pointer to the interpolated point 
    REAL(DP), INTENT(OUT) :: DERIV_NORM !<On exit, the derivative norm of the coordinate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: component_idx,NUMBER_OF_COMPONENTS
    REAL(DP) :: FOCUS,SL,SM
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("COORDINATE_DERIVATIVE_NORM",err,error,*999)

    DERIV_NORM=0.0_DP
    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(ASSOCIATED(INTERPOLATED_POINT)) THEN
        IF(INTERPOLATED_POINT%partialDerivativeType>=FIRST_PART_DERIV) THEN
          IF(PART_DERIV_INDEX<=INTERPOLATED_POINT%maximumPartialDerivativeIndex) THEN
            NUMBER_OF_COMPONENTS=SIZE(INTERPOLATED_POINT%VALUES,1)
            SELECT CASE(PART_DERIV_INDEX)
            CASE(PART_DERIV_S1,PART_DERIV_S2,PART_DERIV_S3)
              SELECT CASE(COORDINATE_SYSTEM%TYPE)
              CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
                DO component_idx=1,NUMBER_OF_COMPONENTS
                  DERIV_NORM=DERIV_NORM+INTERPOLATED_POINT%VALUES(component_idx,PART_DERIV_INDEX)**2
                ENDDO !component_index
              CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
                IF(NUMBER_OF_COMPONENTS==2) THEN
                  DERIV_NORM=INTERPOLATED_POINT%VALUES(1,PART_DERIV_INDEX)**2+(INTERPOLATED_POINT%VALUES(1,1)* &
                    & INTERPOLATED_POINT%VALUES(2,PART_DERIV_INDEX))**2
                ELSE IF(NUMBER_OF_COMPONENTS==3) THEN
                  DERIV_NORM=INTERPOLATED_POINT%VALUES(1,PART_DERIV_INDEX)**2+(INTERPOLATED_POINT%VALUES(1,1)* &
                    & INTERPOLATED_POINT%VALUES(2,PART_DERIV_INDEX))**2+INTERPOLATED_POINT%VALUES(3,PART_DERIV_INDEX)**2
                ELSE
                  LOCAL_ERROR="The number of components for the interpolated point of "// &
                    & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_COMPONENTS,"*",err,error))// &
                    & " is invalid for a cylindrical polar coordinate system."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                ENDIF
              CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
                DERIV_NORM=INTERPOLATED_POINT%VALUES(1,PART_DERIV_INDEX)**2+(INTERPOLATED_POINT%VALUES(1,1)* &
                  & INTERPOLATED_POINT%VALUES(2,PART_DERIV_INDEX)*COS(INTERPOLATED_POINT%VALUES(3,1)))**2+ &
                  & (INTERPOLATED_POINT%VALUES(1,1)*INTERPOLATED_POINT%VALUES(3,PART_DERIV_INDEX))**2
              CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
                FOCUS=COORDINATE_SYSTEM%FOCUS
                SL=SINH(INTERPOLATED_POINT%VALUES(1,1))
                SM=SIN(INTERPOLATED_POINT%VALUES(2,1))
                DERIV_NORM=FOCUS*FOCUS*((SL*SL+SM*SM)*(INTERPOLATED_POINT%VALUES(1,PART_DERIV_INDEX)**2+ &
                  & INTERPOLATED_POINT%VALUES(2,PART_DERIV_INDEX))**2)+(SL*SM*INTERPOLATED_POINT%VALUES(3,PART_DERIV_INDEX))**2
              CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE DEFAULT
                LOCAL_ERROR="The coordinate system type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM%TYPE,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(LOCAL_ERROR,err,error,*999)
              END SELECT
              DERIV_NORM=SQRT(DERIV_NORM)
            CASE DEFAULT
              LOCAL_ERROR="The partial derivative index of "//TRIM(NUMBER_TO_VSTRING(PART_DERIV_INDEX,"*",err,error))// &
                & " is invalid."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          ELSE
            LOCAL_ERROR="The partial derivative index of "//TRIM(NUMBER_TO_VSTRING(PART_DERIV_INDEX,"*",err,error))// &
              & " is invalid. The interpolated point has a maximum number of partial derivatives of "// &
              & TRIM(NUMBER_TO_VSTRING(INTERPOLATED_POINT%maximumPartialDerivativeIndex,"*",err,error))//"."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("The point has not been interpolated to include first derivative values.",err,error,*999)
        ENDIF          
      ELSE
        CALL FlagError("Interpolated point is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Coordinate system is not associated.",err,error,*999)
    ENDIF
    
    EXITS("COORDINATE_DERIVATIVE_NORM")
    RETURN
999 ERRORSEXITS("COORDINATE_DERIVATIVE_NORM",err,error)
    RETURN 1
  END SUBROUTINE COORDINATE_DERIVATIVE_NORM

  !
  !================================================================================================================================
  !

  !>Adjusts the interpolation for non-rectangular cartesian coordinate systems.
  SUBROUTINE COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,PARTIAL_DERIVATIVE_INDEX,VALUE,err,error,*)
  
    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to adjust
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIVATIVE_INDEX !<The partial derivative index to adjust
    REAL(DP), INTENT(INOUT) :: VALUE !<On entry, the coordinate value to adjust. On exit, the adjusted value.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    REAL(DP) :: COSHX,CSS,D,DES,FOCUS,R,SS,SINHX,THETA
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("COORDINATE_INTERPOLATION_ADJUST",err,error,*999)

    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      SELECT CASE(COORDINATE_SYSTEM%TYPE)
      CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
        !Do nothing
      CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
        SELECT CASE(COORDINATE_SYSTEM%radialInterpolationType)
        CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
          !Do nothing
        CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
          R=SQRT(VALUE)
          IF(PARTIAL_DERIVATIVE_INDEX==NO_PART_DERIV) THEN
            VALUE=R
          ELSE
            VALUE=VALUE/(2.0_DP*R)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The radial interpolation type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM% &
            & radialInterpolationType,"*",err,error))//" is invalid for a cylindrical coordinate system."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
        SELECT CASE(COORDINATE_SYSTEM%radialInterpolationType)
        CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
          !Do nothing
        CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
          R=VALUE**(1.0_DP/3.0_DP)
          IF(PARTIAL_DERIVATIVE_INDEX==NO_PART_DERIV) THEN
            VALUE=R
          ELSE
            VALUE=VALUE/(3.0_DP*R*R)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The radial interpolation type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM% &
            & radialInterpolationType,"*",err,error))//" is invalid for a cylindrical/spherical coordinate system."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
        SELECT CASE(COORDINATE_SYSTEM%radialInterpolationType)
        CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
          !Do nothing
        CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
          FOCUS=COORDINATE_SYSTEM%FOCUS
          SS=VALUE/(FOCUS*FOCUS)
          SINHX=SQRT(SS)
          COSHX=SQRT(1.0_DP+SS)
          IF(PARTIAL_DERIVATIVE_INDEX==NO_PART_DERIV) THEN
            VALUE=LOG(SINHX+COSHX)
          ELSE
            VALUE=VALUE/(2.0_DP*FOCUS*FOCUS*SINHX*COSHX)
          ENDIF
        CASE(COORDINATE_RADIAL_CUBED_INTERPOLATION_TYPE)
          FOCUS=COORDINATE_SYSTEM%FOCUS
          CSS=VALUE/(FOCUS**3.0_DP)
          DES=CSS*CSS-4.0_DP/27.0_DP
          IF(DES>0.0_DP) THEN
            D=((CSS+SQRT(DES))/2.0_DP)**(1.0_DP/3.0_DP)
            COSHX=D+1.0_DP/(3.0_DP*D)
          ELSE
            THETA=ACOS(CSS*SQRT(27.0_DP)/2.0_DP)
            COSHX=2.0_DP/SQRT(3.0_DP)*COS(THETA/3.0_DP)
          ENDIF
          SINHX=SQRT(ABS(COSHX*COSHX-1.0_DP))
          IF(PARTIAL_DERIVATIVE_INDEX==NO_PART_DERIV) THEN
            VALUE=LOG(SINHX+COSHX)
          ELSE
            VALUE=VALUE/((3.0_DP*COSHX*COSHX-1.0_DP)*SINHX)/(FOCUS**3.0_DP)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The radial interpolation type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM% &
            & radialInterpolationType,"*",err,error))//" is invalid for a prolate spheroidal coordinate system."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        LOCAL_ERROR="The coordinate system type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM%TYPE,"*",err,error))// &
          & " is invalid."
        CALL FlagError(LOCAL_ERROR,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Coordinate system is not associated.",err,error,*999)
    ENDIF
      
    EXITS("COORDINATE_INTERPOLATION_ADJUST")
    RETURN
999 ERRORSEXITS("COORDINATE_INTERPOLATION_ADJUST",err,error)
    RETURN 1
  END SUBROUTINE COORDINATE_INTERPOLATION_ADJUST

  !
  !================================================================================================================================
  !

  !>Adjusts the interpolation parameters for non-rectangular cartesian coordinate systems.
  SUBROUTINE COORDINATE_INTERPOLATION_PARAMETERS_ADJUST(COORDINATE_SYSTEM,INTERPOLATION_PARAMETERS,err,error,*)
  
    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: COORDINATE_SYSTEM !<A pointer to the coordinate system to adjust
    TYPE(FieldInterpolationParametersType), POINTER :: INTERPOLATION_PARAMETERS !<A pointer to the interpolation parameters to adjust
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("COORDINATE_INTERPOLATION_PARAMETERS_ADJUST",err,error,*999)

!!TODO: Tidy up element parameters for non-rc coordinate systems. See bottom of XPXE and ZPZE.
    
    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      IF(ASSOCIATED(INTERPOLATION_PARAMETERS)) THEN
        SELECT CASE(COORDINATE_SYSTEM%TYPE)
        CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
          !Do nothing
        CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
          SELECT CASE(COORDINATE_SYSTEM%radialInterpolationType)
          CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
            !Do nothing
          CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
          CASE DEFAULT
            LOCAL_ERROR="The radial interpolation type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM% &
              & radialInterpolationType,"*",err,error))//" is invalid for a cylindrical coordinate system."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
          CALL FlagError("Not implemented",err,error,*999)
        CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
          SELECT CASE(COORDINATE_SYSTEM%radialInterpolationType)
          CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
            !Do nothing
          CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
          CASE DEFAULT
            LOCAL_ERROR="The radial interpolation type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM% &
              & radialInterpolationType,"*",err,error))//" is invalid for a spherical coordinate system."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
          SELECT CASE(COORDINATE_SYSTEM%radialInterpolationType)
          CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
            !Do nothing
          CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
          CASE(COORDINATE_RADIAL_CUBED_INTERPOLATION_TYPE)
          CASE DEFAULT
            LOCAL_ERROR="The radial interpolation type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM% &
              & radialInterpolationType,"*",err,error))//" is invalid for a prolate spheroidal coordinate system."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          LOCAL_ERROR="The coordinate system type of "//TRIM(NUMBER_TO_VSTRING(COORDINATE_SYSTEM%TYPE,"*",err,error))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Interpolation parameters is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Coordinate system is not associated.",err,error,*999)
    ENDIF
      
    EXITS("COORDINATE_INTERPOLATION_PARAMETERS_ADJUST")
    RETURN
999 ERRORSEXITS("COORDINATE_INTERPOLATION_PARAMETERS_ADJUST",err,error)
    RETURN 1
  END SUBROUTINE COORDINATE_INTERPOLATION_PARAMETERS_ADJUST

  !
  !================================================================================================================================
  !
 
  !>Calculates the tensor to get from material coordinate system, nu, to local coordinate system, xi.
  SUBROUTINE Coordinates_MaterialSystemCalculate(geometricInterpPointMetrics,fibreInterpPoint,dNudX,dXdNu,dNudXi,dXidNu,err,error,*)
  
    !Argument variables
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics !<The geometric interpolation point metrics at the point to calculate the material coordinate system from.
    TYPE(FieldInterpolatedPointType), POINTER :: fibreInterpPoint !<The fibre interpolation point at the point to calculate the material coordinate system from
    REAL(DP), INTENT(OUT) :: dNudX(:,:) !<dNudX(nuIdx,xIdx). On return, the tensor to transform from the material system to the geometric coordinate system
    REAL(DP), INTENT(OUT) :: dXdNu(:,:) !<dXdNu(xIdx,nuIdx). On return, the tensor to transform from the geometric coordinate system to the material coordinate system
    REAL(DP), INTENT(OUT) :: dNudXi(:,:) !<dNudXi(nuIdx,xiIdx). On return, the tensor to transform from the material system to the xi coordinate system
    REAL(DP), INTENT(OUT) :: dXidNu(:,:) !<dXidNu(xiIdx,nuIdx). On return, the tensor to transform from the xi coordinate system to the material coordinate system
    INTEGER(INTG), INTENT(OUT) :: err   !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) ::  error   !<The error string
    !Local variables
    INTEGER(INTG) :: numberOfXDimensions,numberOfXiDimensions,numberOfNuDimensions,xiIdx
    REAL(DP) :: dNudXiTemp(3,3),Jnuxi,JNuX
    TYPE(VARYING_STRING) :: localError 
     
    ENTERS("Coordinates_MaterialSystemCalculate",err,error,*999)
    
    IF(ASSOCIATED(geometricInterpPointMetrics)) THEN
      
      numberOfXDimensions=geometricInterpPointMetrics%numberOfXDimensions
      numberOfXiDimensions=geometricInterpPointMetrics%numberOfXiDimensions
      
      !Calculate dX/dNu and its inverse dNu/dX (same as transpose due to orthogonality)
      
      !The fibre interpolated point might not be used for isotropic constitutive relations
      IF(ASSOCIATED(fibreInterpPoint)) THEN
        !We have a fibre field
        numberOfNuDimensions=SIZE(fibreInterpPoint%values,1)
        SELECT CASE(numberOfXDimensions)
        CASE(1)
          dXdNu(1,1)=1.0_DP
        CASE(2)
          CALL Coordinates_MaterialSystemCalculatedXdNu2D(geometricInterpPointMetrics,fibreInterpPoint%values(1: &
            & numberOfNuDimensions,1),dXdNu(1:numberOfXDimensions,1:numberOfXDimensions),err,error,*999)
        CASE(3)
          CALL Coordinates_MaterialSystemCalculatedXdNu3D(geometricInterpPointMetrics,fibreInterpPoint%values(1: &
            & numberOfNuDimensions,1),dXdNu(1:numberOfXDimensions,1:numberOfXDimensions),err,error,*999)
        CASE DEFAULT
          localError="The number of dimensions in the geometric interpolated point of "// &
            & TRIM(NumberToVString(numberOfXDimensions,"*",err,error))// &
            & " is invalid. The number of dimensions must be >= 1 and <= 3."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        !No fibre field
        numberOfNuDimensions=0
        DO xiIdx=1,numberOfXiDimensions
          dNudXiTemp(1:numberOfXDimensions,xiIdx)=geometricInterpPointMetrics%dXdXi(1:numberOfXDimensions,xiIdx)
          CALL Normalise(dNudXiTemp(1:numberOfXDimensions,xiIdx),dXdNu(1:numberOfXDimensions,xiIdx),err,error,*999)
        ENDDO !xiIdx
      ENDIF
      !Calculate dNu/dX the inverse of dX/dNu (same as transpose due to orthogonality)
      CALL MatrixTranspose(dXdNu(1:numberOfXDimensions,1:numberOfXDimensions),dNudX(1:numberOfXDimensions,1: &
        & numberOfXDimensions),err,error,*999)
      !Calculate dNu/dXi = dNu/dX * dX/dXi and its inverse dXi/dNu
      CALL MatrixProduct(dNudX(1:numberOfXDimensions,1:numberOfXDimensions), &
        & geometricInterpPointMetrics%dXdXi(1:numberOfXDimensions,1:numberOfXiDimensions), &
        & dNudXiTemp(1:numberOfXDimensions,1:numberOfXiDimensions),err,error,*999)
      !Setup dNudXi
      CALL IdentityMatrix(dNudXi,err,error,*999)
      dNudXi(1:numberOfXDimensions,1:numberOfXiDimensions)=dNudXiTemp(1:numberOfXDimensions,1:numberOfXiDimensions)

      IF(numberOfXDimensions==numberOfXiDimensions) THEN
        CALL Invert(dNudXi(1:numberOfXDimensions,1:numberOfXiDimensions),dXidNu(1:numberOfXiDimensions,1:numberOfXDimensions), &
          & JNuXi,err,error,*999)
      ELSE
        CALL FlagError("Not implemented",err,error,*999)
      ENDIF

      IF(DIAGNOSTICS1) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Calculated material coordinate system:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of X dimensions  = ",numberOfXDimensions,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of Xi dimensions = ",numberOfXiDimensions,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of Nu dimensions = ",numberOfNuDimensions,err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Derivative of X wrt Nu:",err,error,*999)
        CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfXDimensions,1,1,numberOfXDimensions, &
          & numberOfXDimensions,numberOfXDimensions,dXdNu,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
          & '("    dX_dNu','(",I1,",:)','  :",3(X,E13.6))','(18X,3(X,E13.6))',err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Derivative of Nu wrt X:",err,error,*999)
        CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfXDimensions,1,1,numberOfXDimensions, &
          & numberOfXDimensions,numberOfXDimensions,dNudX,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
          & '("    dNu_dX','(",I1,",:)','  :",3(X,E13.6))','(18X,3(X,E13.6))',err,error,*999)
        CALL Determinant(dNudX(1:numberOfXDimensions,1:numberOfXDimensions),JNuX,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Determinant dNu_dX, JNuX = ", JNuX,err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Derivative of Nu wrt Xi:",err,error,*999)
        CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfXDimensions,1,1,numberOfXiDimensions, &
          & numberOfXiDimensions,numberOfXiDimensions,dNudXi,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
          & '("    dNu_dXi','(",I1,",:)',' :",3(X,E13.6))','(18X,3(X,E13.6))',err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Derivative of Xi wrt Nu:",err,error,*999)
        CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfXiDimensions,1,1,numberOfXDimensions, &
          & numberOfXDimensions,numberOfXDimensions,dXidNu,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
          & '("    dXi_dNu','(",I1,",:)',' :",3(X,E13.6))','(18X,3(X,E13.6))',err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Determinant dNu_dXi, JNuXi = ",JNuXi,err,error,*999)
      ENDIF
      
    ELSE
      CALL FlagError("Geometric interpolated point metrics is not associated.",err,error,*999)
    ENDIF    
    
    EXITS("Coordinates_MaterialSystemCalculate")
    RETURN
999 ERRORSEXITS("Coordinates_MaterialSystemCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE Coordinates_MaterialSystemCalculate
  
  !
  !================================================================================================================================
  !

  !>Calculates transformation between spatial CS and rotated reference orthogonal material CS in 2D space
  SUBROUTINE Coordinates_MaterialSystemCalculatedXdNu2D(geometricInterpPointMetrics,angle,dXdNu,err,error,*)

    !Argument variables
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics !<The geometric interpolated point metrics at the point to calculate dXdNu at. 
    REAL(DP), INTENT(IN) :: angle(:) !<angle(fibreIdx). The fibre angle (in radians) 
    REAL(DP), INTENT(OUT) :: dXdNu(:,:) !<dXdNu(coordinateIdx,coordinateIdx). On exit, the dXdNu tensor.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: det,dXdNuR(2,2),R(2,2),f(2),g(2)

    ENTERS("Coordinates_MaterialSystemCalculatedXdNu2D",err,error,*999)

    IF(ASSOCIATED(geometricInterpPointMetrics)) THEN
    
      !First calculate reference material CS
      
      !Reference material direction 1.
      f(1:2) = [ geometricInterpPointMetrics%dXdXi(1,1),geometricInterpPointMetrics%dXdXi(2,1) ]
      
      !Compute (normalised) vector orthogonal to material direction 1 to form material direction 2
      g(1:2) = [ -1.0_DP*f(2),f(1) ]
      
      CALL Normalise(f(1:2),dXdNuR(1:2,1),err,error,*999)
      CALL Normalise(g(1:2),dXdNuR(1:2,2),err,error,*999)
      
      !Rotate by multiply with rotation matrix
      R(:,1) = [ COS(angle(1)),-1.0_DP*SIN(angle(1)) ]
      R(:,2) = [ SIN(angle(1)),COS(angle(1)) ]
      
      CALL MatrixProduct(R,dXdNuR,dXdNu,err,error,*999)
      
      CALL Normalise(dXdNu(1:2,1),dXdNu(1:2,1),err,error,*999)
      CALL Normalise(dXdNu(1:2,2),dXdNu(1:2,2),err,error,*999)

      IF(diagnostics1) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"2D material system calculation:",err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Reference material directions:",err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,2,2,f,'("    f              :",2(X,E13.6))','(20X,2(X,E13.6))', &
          & err,error,*999)      
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,2,2,g,'("    g              :",2(X,E13.6))','(20X,2(X,E13.6))', &
          & err,error,*999)      
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Derivative of X wrt Nu (reference):",err,error,*999)
        CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,1,1,2,2,2,dXdNuR,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
          & '("      dX_dNuR','(",I1,",:)',' :",2(X,E13.6))','(20X,2(X,E13.6))',err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Fibre calculation:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Fibre angle = ",angle(1),err,error,*999)
        IF(diagnostics2) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Rotation matrix:",err,error,*999)
          CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,1,1,2,2,2,R,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
            & '("      R','(",I1,",:)','       :",2(X,E13.6))','(20X,2(X,E13.6))',err,error,*999)
        ENDIF
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Derivative of X wrt Nu (material):",err,error,*999)
        CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,1,1,2,2,2,dXdNu,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
          & '("      dX_dNu','(",I1,",:)','   :",2(X,E13.6))','(20X,2(X,E13.6))',err,error,*999)
        CALL Determinant(dXdNu,det,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Determinant dX_dNu = ",det,err,error,*999)
      ENDIF
      
    ELSE
      CALL FlagError("Geometry interpolated point metrics is not associated.",err,error,*999)
    ENDIF
        
    EXITS("Coordinates_MaterialSystemCalculatedXdNu2D")
    RETURN
999 ERRORSEXITS("Coordinates_MaterialSystemCalculatedXdNu2D",err,error)
    RETURN 1
    
  END SUBROUTINE Coordinates_MaterialSystemCalculatedXdNu2D

  !
  !================================================================================================================================
  !

  !>Calculates transformation between spatial CS and rotated reference orthogonal material CS in 3D space
  SUBROUTINE Coordinates_MaterialSystemCalculatedXdNu3D(geometricInterpPointMetrics,angle,dXdNu,err,error,*)

    !Argument variables
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics !<The geometric interpolated point metrics at the point to calculate dXdNu at. 
    REAL(DP), INTENT(IN) :: angle(:) !<angles(fibreIdx). The fibre, imbrication and sheet (in radians) 
    REAL(DP), INTENT(OUT) :: dXdNu(:,:) !<dXdNu(coordinateIdx,coordinateIdx). On exit, the dXdNu tensor.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: angles(3),det,dXdNu2(3,3),dXdNu3(3,3),dXdNuR(3,3),f(3),g(3),h(3), &
      & Raf(3,3),Rbf(3,3),Rai(3,3),Rbi(3,3),Ras(3,3),Rbs(3,3)
    
    ENTERS("Coordinates_MaterialSystemCalculatedXdNu3D",err,error,*999)
    
    IF(ASSOCIATED(geometricInterpPointMetrics)) THEN

      !Calculate fibre, imbrication and sheet angles (allowing for missing angles)
      angles=0.0_DP
      angles(1:SIZE(angle,1))=angle(1:SIZE(angle,1))
      
      !First calculate reference material CS
      f(1:3)=geometricInterpPointMetrics%dXdXi(1:3,1) !reference material direction 1.
      CALL CrossProduct(geometricInterpPointMetrics%dXdXi(1:3,1),geometricInterpPointMetrics%dXdXi(1:3,2),h,err,error,*999) !reference material direction 3.    
      CALL CrossProduct(h,f,g,err,error,*999) !reference material direction 2.      

      CALL Normalise(f,dXdNuR(1:3,1),err,error,*999)
      CALL Normalise(g,dXdNuR(1:3,2),err,error,*999)
      CALL Normalise(h,dXdNuR(1:3,3),err,error,*999)
      
      IF(diagnostics1) THEN
      ENDIF

      !FIBRE ANGLE(alpha) - angles(1) 
      !In order to rotate reference material CS by alpha(fibre angle) in anti-clockwise  
      !direction about its axis 3, following steps are performed.
      !(a) first align reference material direction 3 with Z(spatial) axis by rotating the ref material CS. 
      !(b) then rotate the aligned material CS by alpha about Z axis in anti-clockwise direction
      !(c) apply the inverse of step(a) to the CS in (b)
      !It can be shown that steps (a),(b) and (c) are equivalent to post-multiplying
      !rotation in (a) by rotation in (b). i.e. Ra*Rb  
      
      !The normalised reference material CS contains the transformation(rotation) between 
      !the spatial CS -> reference material CS. i.e. Raf
      Raf=dXdNuR
        
      !Initialise rotation matrix Rbf
      CALL IdentityMatrix(Rbf,err,error,*999)
      !Populate rotation matrix Rbf about axis 3 (Z)
      Rbf(1,1)=COS(angles(1))
      Rbf(1,2)=-1.0_DP*SIN(angles(1))
      Rbf(2,1)=SIN(angles(1))
      Rbf(2,2)=COS(angles(1))
        
      CALL MatrixProduct(Raf,Rbf,dXdNu3,err,error,*999)  

      !IMBRICATION ANGLE (beta) - angles(2)     
      !In order to rotate alpha-rotated material CS by beta(imbrication angle) in anti-clockwise  
      !direction about its new axis 2, following steps are performed.
      !(a) first align new material direction 2 with Y(spatial) axis by rotating the new material CS. 
      !(b) then rotate the aligned CS by beta about Y axis in anti-clockwise direction
      !(c) apply the inverse of step(a) to the CS in (b)
      !As mentioned above, (a),(b) and (c) are equivalent to post-multiplying
      !rotation in (a) by rotation in (b). i.e. Rai*Rbi  
      
      !dXdNu3 contains the transformation(rotation) between 
      !the spatial CS -> alpha-rotated reference material CS. i.e. Rai
      Rai=dXdNu3
      !Initialise rotation matrix Rbi
      CALL IdentityMatrix(Rbi,err,error,*999)
      !Populate rotation matrix Rbi about axis 2 (Y). Note the sign change
      Rbi(1,1)=COS(angles(2))
      Rbi(1,3)=SIN(angles(2))
      Rbi(3,1)=-1.0_DP*SIN(angles(2))
      Rbi(3,3)=COS(angles(2))
        
      CALL MatrixProduct(Rai,Rbi,dXdNu2,err,error,*999)  

      !SHEET ANGLE (gamma) - angles(3)    
      !In order to rotate alpha-beta-rotated material CS by gama(sheet angle) in anti-clockwise  
      !direction about its new axis 1, following steps are performed.
      !(a) first align new material direction 1 with X(spatial) axis by rotating the new material CS. 
      !(b) then rotate the aligned CS by gama about X axis in anti-clockwise direction
      !(c) apply the inverse of step(a) to the CS in (b)
      !Again steps (a),(b) and (c) are equivalent to post-multiplying
      !rotation in (a) by rotation in (b). i.e. Ras*Rbs  
      
      !dXdNu2 contains the transformation(rotation) between 
      !the spatial CS -> alpha-beta-rotated reference material CS. i.e. Ras
      Ras=dXdNu2
      !Initialise rotation matrix Rbs
      CALL IdentityMatrix(Rbs,err,error,*999)
      !Populate rotation matrix Rbs about axis 1 (X). 
      Rbs(2,2)=COS(angles(3))
      Rbs(2,3)=-1.0_DP*SIN(angles(3))
      Rbs(3,2)=SIN(angles(3))
      Rbs(3,3)=COS(angles(3))
      
      CALL MatrixProduct(Ras,Rbs,dXdNu,err,error,*999)  

      IF(diagnostics1) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"3D material system calculation:",err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Reference material directions:",err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,3,3,f,'("    f              :",3(X,E13.6))','(20X,3(X,E13.6))', &
          & err,error,*999)      
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,3,3,g,'("    g              :",3(X,E13.6))','(20X,3(X,E13.6))', &
          & err,error,*999)      
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,3,3,h,'("    h              :",3(X,E13.6))','(20X,3(X,E13.6))', &
          & err,error,*999)      
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Derivative of X wrt Nu (reference):",err,error,*999)
        CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3,3,3,dXdNuR,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
          & '("      dX_dNuR','(",I1,",:)',' :",3(X,E13.6))','(20X,3(X,E13.6))',err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Fibre calculation:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Fibre angle = ",angles(1),err,error,*999)
        IF(diagnostics2) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Rotation matrix A:",err,error,*999)
          CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3,3,3,Raf,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
            & '("      Ra','(",I1,",:)','      :",3(X,E13.6))','(20X,3(X,E13.6))',err,error,*999)
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Rotation matrix B:",err,error,*999)
          CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3,3,3,Rbf,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
            & '("      Rb','(",I1,",:)','      :",3(X,E13.6))','(20X,3(X,E13.6))',err,error,*999)
        ENDIF
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Derivative of X wrt Nu (after alpha rotation):",err,error,*999)
        CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3,3,3,dXdNu3,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
          & '("      dX_dNu3','(",I1,",:)',' :",3(X,E13.6))','(20X,3(X,E13.6))',err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Imbrication calculation:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Imbrication angle = ",angles(2),err,error,*999)
        IF(diagnostics2) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Rotation matrix A:",err,error,*999)
          CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3,3,3,Rai,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
            & '("      Ra','(",I1,",:)','      :",3(X,E13.6))','(20X,3(X,E13.6))',err,error,*999)
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Rotation matrix B:",err,error,*999)
          CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3,3,3,Rbi,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
            & '("      Rb','(",I1,",:)','      :",3(X,E13.6))','(20X,3(X,E13.6))',err,error,*999)
        ENDIF
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Derivative of X wrt Nu (after alpha-beta rotation):",err,error,*999)
        CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3,3,3,dXdNu2,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
          & '("      dX_dNu2','(",I1,",:)',' :",3(X,E13.6))','(20X,3(X,E13.6))',err,error,*999)
         CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Sheet calculation:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Sheet angle = ",angles(3),err,error,*999)
        IF(diagnostics2) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Rotation matrix A:",err,error,*999)
          CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3,3,3,Ras,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
            & '("      Ra','(",I1,",:)','      :",3(X,E13.6))','(20X,3(X,E13.6))',err,error,*999)
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Rotation matrix B:",err,error,*999)
          CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3,3,3,Rbs,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
            & '("      Rb','(",I1,",:)','      :",3(X,E13.6))','(20X,3(X,E13.6))',err,error,*999)
        ENDIF
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Derivative of X wrt Nu (material):",err,error,*999)
        CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3,3,3,dXdNu,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
          & '("      dX_dNu','(",I1,",:)','  :",3(X,E13.6))','(20X,3(X,E13.6))',err,error,*999)
        CALL Determinant(dXdNu,det,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Determinant dX_dNu = ",det,err,error,*999)
      ENDIF
            
    ELSE
      CALL FlagError("Geometry interpolated point metrics is not associated.",err,error,*999)
    ENDIF
    
    EXITS("Coordinates_MaterialSystemCalculatedXdNu3D")
    RETURN
999 ERRORSEXITS("Coordinates_MaterialSystemCalculatedXdNu3D",err,error)
    RETURN 1
    
  END SUBROUTINE Coordinates_MaterialSystemCalculatedXdNu3D

  !
  !================================================================================================================================
  !

  !>Finalises the coordinate systems and destroys all coordinate systems.
  SUBROUTINE CoordinateSystems_Finalise(coordinateSystems,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemsType), POINTER :: coordinateSystems !<A pointer to the coordinate systems to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coordinateSystemIdx
    
    ENTERS("CoordinateSystems_Finalise",err,error,*999)

    IF(ASSOCIATED(coordinateSystems)) THEN
      DO coordinateSystemIdx=1,coordinateSystems%numberOfCoordinateSystems
        CALL COORDINATE_SYSTEM_FINALISE(coordinateSystems%coordinateSystems(coordinateSystemIdx)%ptr,err,error,*999)
      ENDDO !coordinateSystemIdx
      DEALLOCATE(coordinateSystems%coordinateSystems)
      DEALLOCATE(coordinateSystems)
    ENDIF
  
    EXITS("CoordinateSystems_Finalise")
    RETURN
999 ERRORSEXITS("CoordinateSystems_Finalise",err,error)
    RETURN 1
  END SUBROUTINE CoordinateSystems_Finalise

  !
  !================================================================================================================================
  !
   
  !>Initialises the coordinate systems and creates the world coordinate system for a context.
  SUBROUTINE CoordinateSystems_Initialise(context,err,error,*)

    !Argument variables
    TYPE(ContextType), POINTER :: context !<The context to initialise the the coordinate systems for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("CoordinateSystems_Initialise",err,error,*998)

    IF(.NOT.ASSOCIATED(context)) CALL FlagError("Context is not associated.",err,error,*998)
    IF(ASSOCIATED(context%coordinateSystems)) &
      & CALL FlagError("Context coordinate systems is already associated.",err,error,*998)
    
    !Allocate the coordinate systems
    ALLOCATE(context%coordinateSystems,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate coordinate systems.",err,error,*999)
    !Initialise
    context%coordinateSystems%context=>context
    context%coordinateSystems%numberOfCoordinateSystems=0
    NULLIFY(context%coordinateSystems%coordinateSystems)
    !Add in a world coordinate system
    ALLOCATE(context%coordinateSystems%coordinateSystems(1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate list of coordinate systems.",err,error,*999)
    !Create the default RC World cooordinate system
    ALLOCATE(context%coordinateSystems%coordinateSystems(1)%ptr,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate world coordinate system.",err,error,*999)
    context%coordinateSystems%coordinateSystems(1)%ptr%userNumber=0
    context%coordinateSystems%coordinateSystems(1)%ptr%type=COORDINATE_RECTANGULAR_CARTESIAN_TYPE
    context%coordinateSystems%coordinateSystems(1)%ptr%numberOfDimensions=3
    context%coordinateSystems%coordinateSystems(1)%ptr%focus=1.0_DP    
    context%coordinateSystems%coordinateSystems(1)%ptr%origin=[0.0_DP,0.0_DP,0.0_DP]
    context%coordinateSystems%coordinateSystems(1)%ptr%orientation=RESHAPE(&
      & [1.0_DP,0.0_DP,0.0_DP, &
      &  0.0_DP,1.0_DP,0.0_DP, &
      &  0.0_DP,0.0_DP,1.0_DP], &
      & [3,3])    
    context%coordinateSystems%coordinateSystems(1)%ptr%coordinateSystemFinished=.TRUE.
    context%coordinateSystems%numberOfCoordinateSystems=1
   
    EXITS("CoordinateSystems_Initialise")
    RETURN
999 CALL CoordinateSystems_Finalise(context%coordinateSystems,dummyErr,dummyError,*998)
998 ERRORSEXITS("CoordinateSystems_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystems_Initialise

  !
  !================================================================================================================================
  !
  
END MODULE COORDINATE_ROUTINES

