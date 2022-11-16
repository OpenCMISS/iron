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

!> \defgroup OpenCMISS_CoordinateSystem OpenCMISS::Iron::CoordinateSystem
!> This module contains all coordinate transformation and support routines.
MODULE CoordinateSystemRoutines

  USE BaseRoutines
  USE Constants
  USE CoordinateSystemAccessRoutines
  USE FieldAccessRoutines
  USE InputOutput
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

  !Module types
  
  !Module variables

  !Interfaces

  !>Coordinates_ConvertFromRC performs a coordinate transformation from a rectangular cartesian coordinate at the point with
  !>coordinate Z(:) to the returned point with coordinate X(:) in the coordinate system identified by COORDINATE_SYSTEM.
  INTERFACE CoordinateSystem_ConvertFromRC
    MODULE PROCEDURE CoordinateSystem_ConvertFromRCDP
    MODULE PROCEDURE CoordinateSystem_ConvertFromRCSP
  END INTERFACE CoordinateSystem_ConvertFromRC

  !>CoordinateSystem_ConvertToRC performs a coordinate transformation from a coordinate system identified by COORDINATE_SYSTEM
  !>at the point X(:) to the returned point Z(:) in rectangular cartesian coordinates.
  INTERFACE CoordinateSystem_ConvertToRC
    MODULE PROCEDURE CoordinateSystem_ConvertToRCDP
    MODULE PROCEDURE CoordinateSystem_ConvertToRCSP
  END INTERFACE CoordinateSystem_ConvertToRC

  !>Calculates the difference (or delta) between two points in a coordinate system. Discontinuities for polar coordinate
  !>systems are accounted for
  INTERFACE CoordinateSystem_DeltaCalculate
    MODULE PROCEDURE CoordinateSystem_DeltaCalculateDP
    !MODULE PROCEDURE CoordinateSystem_DeltaCalculateSP
  END INTERFACE CoordinateSystem_DeltaCalculate

  !>Calculates dx(:)/dz(i) at x, where z(i) are rectangular cartesian and x(:) are curvilinear coordinates defined by COORDINATE_SYSTEM. \todo CHANGE NAME TO SOMETHING MORE MEANINGFULL?
  INTERFACE dXdZ
    MODULE PROCEDURE dXdZDP
    !MODULE PROCEDURE dXdZSP
  END INTERFACE dXdZ

  !!TODO:: CHANGE NAME TO SOMETHING MORE MEANINGFULL?
  INTERFACE d2ZdX1dX2
    MODULE PROCEDURE d2ZdX1dX2DP
    !MODULE PROCEDURE d2ZdX1dX2SP
  END INTERFACE d2ZdX1dX2

  !!TODO:: CHANGE NAME TO SOMETHING MORE MEANINGFULL?
  INTERFACE dZdX
    MODULE PROCEDURE dZdXDP
    !MODULE PROCEDURE dZdXSP
  END INTERFACE dZdX

  INTERFACE CoordinateSystem_DerivativeConvertToRC
    MODULE PROCEDURE CoordinateSystem_DerivativeConvertToRCDP
    MODULE PROCEDURE CoordinateSystem_DerivativeConvertToRCSP
  END INTERFACE CoordinateSystem_DerivativeConvertToRC

  PUBLIC CoordinateSystem_ConvertFromRC
  
  PUBLIC CoordinateSystem_ConvertToRC

  PUBLIC CoordinateSystem_CreateStart,CoordinateSystem_CreateFinish

  PUBLIC CoordinateSystem_DeltaCalculate

  PUBLIC CoordinateSystem_DerivativeConvertToRC
  
  PUBLIC CoordinateSystem_DerivativeNorm
  
  PUBLIC CoordinateSystem_Destroy
  
  PUBLIC CoordinateSystem_DimensionSet
  
  PUBLIC CoordinateSystem_FocusSet
  
  PUBLIC CoordinateSystem_InterpolationAdjust

  PUBLIC CoordinateSystem_InterpolationParametersAdjust

  PUBLIC CoordinateSystem_MaterialSystemCalculate

  PUBLIC CoordinateSystem_MaterialTransformSymTensor2

  PUBLIC CoordinateSystem_MetricsCalculate
  
  PUBLIC CoordinateSystem_RadialInterpolationTypeSet

  PUBLIC CoordinateSystem_TypeSet

  PUBLIC CoordinateSystem_OriginSet
  
  PUBLIC CoordinateSystem_OrientationSet
  
  PUBLIC CoordinateSystems_Initialise,CoordinateSystems_Finalise

CONTAINS

  !
  !================================================================================================================================
  !

  !>Performs a coordinate transformation from a rectangular cartesian coordinate at the point with coordinate z(:) to the returned point with coordinate X(:) in the coordinate system identified by coordinateSystem for double precision coordinates.
  FUNCTION CoordinateSystem_ConvertFromRCDP(coordinateSystem,z,err,error)
  
    !Argument variables
    TYPE(CoordinateSystemType), INTENT(IN) :: coordinateSystem !<The coordinate system to perform the conversion on
    REAL(DP), INTENT(IN) :: z(:) !<The rectangular cartesian coordiantes to convert
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: CoordinateSystem_ConvertFromRCDP(SIZE(z,1))
    !Local variables
    REAL(DP) :: a1,a2,a3,a4,a5,a6,a7,a8,a9,focus
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("CoordinateSystem_ConvertFromRCDP",err,error,*999)

    CoordinateSystem_ConvertFromRCDP=0.0_DP

    IF(SIZE(z,1)<coordinateSystem%numberOfDimensions) THEN
      localError="The size of z of "//TRIM(NumberToVString(SIZE(z,1),"*",err,error))// &
        & " is less than the number of coordinate system dimensions of "// &
        & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(coordinateSystem%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      CoordinateSystem_ConvertFromRCDP(1:coordinateSystem%numberOfDimensions)=z(1:coordinateSystem%numberOfDimensions)
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      SELECT CASE(coordinateSystem%numberOfDimensions)
      CASE(2)
        CoordinateSystem_ConvertFromRCDP(1)=SQRT(z(1)**2+z(2)**2)
        CoordinateSystem_ConvertFromRCDP(2)=ATAN2(z(1),z(2))
      CASE(3)
        CoordinateSystem_ConvertFromRCDP(1)=SQRT(z(1)**2+z(2)**2)
        CoordinateSystem_ConvertFromRCDP(2)=ATAN2(z(1),z(2))
        CoordinateSystem_ConvertFromRCDP(3)=z(3)
      CASE DEFAULT
        localError="The number of coordinates of "//TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))// &
          & " is invalid for a cylindrical polar coordinate system."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      IF(CoordinateSystem_ConvertFromRCDP(2)<0.0_DP) &
        & CoordinateSystem_ConvertFromRCDP(2)=CoordinateSystem_ConvertFromRCDP(2)+2.0_DP*PI !reference coordinate 0->2*pi
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(coordinateSystem%numberOfDimensions/=3) THEN
        localError="The number of coordinates of "//TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))// &
          & " is invalid for a spherical polar coordinate system."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      CoordinateSystem_ConvertFromRCDP(1)=SQRT(z(1)**2+z(2)**2+z(3)**2)
      IF(ABS(z(1))>=ZERO_TOLERANCE.OR.ABS(z(2))>=ZERO_TOLERANCE) THEN
        CoordinateSystem_ConvertFromRCDP(2)=ATAN2(z(2),z(1))
      ELSE
        CoordinateSystem_ConvertFromRCDP(2)=0.0_DP
      ENDIF
      a1=SQRT(z(1)**2+z(2)**2)
      IF(ABS(z(3))>=ZERO_TOLERANCE.OR.ABS(a1)>=ZERO_TOLERANCE) THEN
        CoordinateSystem_ConvertFromRCDP(3)=ATAN2(z(3),a1)
      ELSE
        CoordinateSystem_ConvertFromRCDP(3)=0.0_DP
      ENDIF
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      focus=coordinateSystem%focus
      IF(coordinateSystem%numberOfDimensions/=3) THEN
        localError="The number of coordinates of "//TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))// &
          & " is invalid for a prolate spheroidal coordinate system."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      a1=z(1)**2+z(2)**2+z(3)**2-focus**2
      a2=SQRT(a1**2+4.0_DP*(focus**2)*(z(2)**2+z(3)**2))
      a3=2.0_DP*focus**2
      a4=MAX((a2+a1)/a3,0.0_DP)
      a5=MAX((a2-a1)/a3,0.0_DP)
      a6=SQRT(a4)
      a7=MIN(SQRT(a5),1.0_DP)
      IF(ABS(a7)<=1.0_DP) THEN
        a8=ASIN(a7)
      ELSE
        a8=0.0_DP
        CALL FlagWarning("Put a8=0 since ABS(a8)>1.",err,error,*999)
      ENDIF
      IF((ABS(z(3))<ZERO_TOLERANCE).OR.(ABS(a6)<ZERO_TOLERANCE).OR.(ABS(a7)<ZERO_TOLERANCE)) THEN
        a9=0.0_DP
      ELSE
        IF(ABS(a6*a7)>0.0_DP) THEN
          a9=z(3)/(focus*a6*a7)
        ELSE
          a9=0.0_DP
          CALL FlagWarning("Put a9=0 since a6*a7=0.",err,error,*999)
        ENDIF
        IF(a9>=1.0_DP) THEN
          a9=PI/2.0_DP
        ELSE IF(a9<=-1.0_DP) THEN
          a9=-PI/2.0_DP
        ELSE
          a9=ASIN(a9)
        ENDIF
      ENDIF
      CoordinateSystem_ConvertFromRCDP(1)=LOG(a6+SQRT(a4+1.0_DP))
      IF(ABS(z(1))>=ZERO_TOLERANCE) THEN
        CoordinateSystem_ConvertFromRCDP(2)=a8
      ELSE
        CoordinateSystem_ConvertFromRCDP(2)=PI-a8
      ENDIF
      IF(ABS(z(2))>ZERO_TOLERANCE) THEN
        CoordinateSystem_ConvertFromRCDP(3)=MOD(a9+2.0_DP*PI,2.0_DP*PI)
      ELSE
        CoordinateSystem_ConvertFromRCDP(3)=PI-a9
      ENDIF
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The coordinate type of "//TRIM(NumberToVString(coordinateSystem%TYPE,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("CoordinateSystem_ConvertFromRCDP")
    RETURN
999 ERRORSEXITS("CoordinateSystem_ConvertFromRCDP",err,error)
    RETURN
    
  END FUNCTION CoordinateSystem_ConvertFromRCDP
  
  !
  !================================================================================================================================
  !

  !>CoordinateSystem_ConvertFromRCSP performs a coordinate transformation from a rectangular cartesian coordinate at the
  !>point with coordinate z(:) to the returned point with coordinate x(:) in the coordinate system identified by
  !>coordinateSystem for single precision coordinates.
  FUNCTION CoordinateSystem_ConvertFromRCSP(coordinateSystem,z,err,error)
  
    !Argument variables
    TYPE(CoordinateSystemType), INTENT(IN) :: coordinateSystem !<The coordinate system to convert from RC to
    REAL(SP), INTENT(IN) :: z(:) !<The coordinate to convert from rectangular cartesian
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(SP) :: CoordinateSystem_ConvertFromRCSP(SIZE(z,1))
    !Local variables
    REAL(SP) :: a1,a2,a3,a4,a5,a6,a7,a8,a9,focus
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("CoordinateSystem_ConvertFromRCSP",err,error,*999)

    CoordinateSystem_ConvertFromRCSP=0.0_SP
    
    IF(SIZE(z,1)<coordinateSystem%numberOfDimensions) THEN
      localError="The size of z of "//TRIM(NumberToVString(SIZE(z,1),"*",err,error))// &
        & " is less than the number of coordinate system dimensions of "// &
        & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(coordinateSystem%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      CoordinateSystem_ConvertFromRCSP(1:coordinateSystem%numberOfDimensions)=z(1:coordinateSystem%numberOfDimensions)
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      SELECT CASE(coordinateSystem%numberOfDimensions)
      CASE(2)
        CoordinateSystem_ConvertFromRCSP(1)=SQRT(z(1)**2+z(2)**2)
        CoordinateSystem_ConvertFromRCSP(2)=ATAN2(z(1),z(2))
      CASE(3)
        CoordinateSystem_ConvertFromRCSP(1)=SQRT(z(1)**2+z(2)**2)
        CoordinateSystem_ConvertFromRCSP(2)=ATAN2(z(1),z(2))
        CoordinateSystem_ConvertFromRCSP(3)=z(3)
      CASE DEFAULT
        localError="The number of coordinates of "//TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))// &
          & " is invalid for a cylindrical polar coordinate system."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      IF(CoordinateSystem_ConvertFromRCSP(2)<0.0_SP)  &
        & CoordinateSystem_ConvertFromRCSP(2)=CoordinateSystem_ConvertFromRCSP(2)+2.0_SP*REAL(PI,SP) !reference coordinate 0->2*pi
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(coordinateSystem%numberOfDimensions/=3) THEN
        localError="The number of coordinates of "//TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))// &
          & " is invalid for a spherical polar coordinate system."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      CoordinateSystem_ConvertFromRCSP(1)=SQRT(z(1)**2+z(2)**2+z(3)**2)
      IF(ABS(z(1))>=ZERO_TOLERANCE_SP.OR.ABS(z(2))>ZERO_TOLERANCE_SP) THEN
        CoordinateSystem_ConvertFromRCSP(2)=ATAN2(z(2),z(1))
      ELSE
        CoordinateSystem_ConvertFromRCSP(2)=0.0_SP
      ENDIF
      a1=SQRT(z(1)**2+z(2)**2)
      IF(ABS(z(3))>=ZERO_TOLERANCE_SP.OR.ABS(a1)>ZERO_TOLERANCE_SP) THEN
        CoordinateSystem_ConvertFromRCSP(3)=ATAN2(z(3),a1)
      ELSE
        CoordinateSystem_ConvertFromRCSP(3)=0.0_SP
      ENDIF
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      focus=REAL(coordinateSystem%focus,SP)
      IF(coordinateSystem%numberOfDimensions/=3) THEN
        localError="The number of coordinates of "//TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))// &
          & " is invalid for a prolate spheroidal coordinate system."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      a1=z(1)**2+z(2)**2+z(3)**2-focus**2
      a2=SQRT(a1**2+4.0_SP*(focus**2)*(z(2)**2+z(3)**2))
      a3=2.0_SP*focus**2
      a4=MAX((a2+a1)/a3,0.0_SP)
      a5=MAX((a2-a1)/a3,0.0_SP)
      a6=SQRT(a4)
      a7=MIN(SQRT(a5),1.0_SP)
      IF(ABS(a7)<=1.0_SP) THEN
        a8=ASIN(a7)
      ELSE
        a8=0.0_SP
        CALL FlagWarning("Put a8=0 since ABS(a8)>1.",err,error,*999)
      ENDIF
      IF((ABS(z(3))<ZERO_TOLERANCE_SP).OR.(ABS(a6)<ZERO_TOLERANCE_SP).OR.(ABS(a7)<ZERO_TOLERANCE_SP)) THEN
        a9=0.0_SP
      ELSE
        IF(ABS(a6*a7)>ZERO_TOLERANCE_SP) THEN
          a9=z(3)/(focus*a6*a7)
        ELSE
          a9=0.0_SP
          CALL FlagWarning("Put a9=0 since a6*a7=0.",err,error,*999)
        ENDIF
        IF(a9>=1.0_SP) THEN
          a9=REAL(PI,SP)/2.0_SP
        ELSE IF(a9<=-1.0_SP) THEN
          a9=-REAL(PI,SP)/2.0_SP
        ELSE
          a9=ASIN(a9)
        ENDIF
      ENDIF
      CoordinateSystem_ConvertFromRCSP(1)=LOG(a6+SQRT(a4+1.0_SP))
      IF(ABS(z(1))>=ZERO_TOLERANCE_SP) THEN
        CoordinateSystem_ConvertFromRCSP(2)=a8
      ELSE
        CoordinateSystem_ConvertFromRCSP(2)=REAL(PI,SP)-a8
      ENDIF
      IF(ABS(z(2))>=ZERO_TOLERANCE_SP) THEN
        CoordinateSystem_ConvertFromRCSP(3)=MOD(a9+2.0_SP*REAL(PI,SP),2.0_SP*&
          & REAL(PI,SP))
      ELSE
        CoordinateSystem_ConvertFromRCSP(3)=REAL(PI,SP)-a9
      ENDIF
   CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The coordinate type of "//TRIM(NumberToVString(coordinateSystem%TYPE,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("CoordinateSystem_ConvertFromRCSP")
    RETURN
999 ERRORSEXITS("CoordinateSystem_ConvertFromRCSP",err,error)
    RETURN
    
  END FUNCTION CoordinateSystem_ConvertFromRCSP

  !
  !================================================================================================================================
  !

  !>CoordinateSystem_ConvertToRCDP performs a coordinate transformation from a coordinate system identified by
  !>coordinateSystem at the point x(:) to the returned point z(:) in rectangular cartesian coordinates for
  !>double precision coordinates.
  FUNCTION CoordinateSystem_ConvertToRCDP(coordinateSystem,x,err,error)
  
    !Argument variables
    TYPE(CoordinateSystemType), INTENT(IN) :: coordinateSystem !<The coordinate system to convert to rectangular cartesian
    REAL(DP), INTENT(IN) :: x(:) !<The coordiante to convert
    INTEGER(INTG), INTENT(OUT) :: err !<The error coode
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: CoordinateSystem_ConvertToRCDP(SIZE(x,1))
    !Local variables
    REAL(DP) :: focus
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("CoordinateSystem_ConvertToRCDP",err,error,*999)
    
    CoordinateSystem_ConvertToRCDP=0.0_DP

    IF(SIZE(x,1)<coordinateSystem%numberOfDimensions) THEN
      localError="The size of x of "//TRIM(NumberToVString(SIZE(x,1),"*",err,error))// &
        & " is less than the number of coordinate system dimensions of "// &
        & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
 
    SELECT CASE(coordinateSystem%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      CoordinateSystem_ConvertToRCDP(1:coordinateSystem%numberOfDimensions)=x(1:coordinateSystem%numberOfDimensions)
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      SELECT CASE(coordinateSystem%numberOfDimensions)
      CASE(2)
        CoordinateSystem_ConvertToRCDP(1)=x(1)*COS(x(2))
        CoordinateSystem_ConvertToRCDP(2)=x(1)*SIN(x(2))
      CASE(3)
        CoordinateSystem_ConvertToRCDP(1)=x(1)*COS(x(2))
        CoordinateSystem_ConvertToRCDP(2)=x(1)*SIN(x(2))
        CoordinateSystem_ConvertToRCDP(3)=x(3)
      CASE DEFAULT
        localError="The number of coordinates of "//TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))// &
          & " is invalid for a cylindrical polar coordinate system."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(coordinateSystem%numberOfDimensions/=3) THEN
        localError="The number of coordinates of "//TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))// &
          & " is invalid for a spherical polar coordinate system."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      CoordinateSystem_ConvertToRCDP(1)=x(1)*COS(x(2))*COS(x(3))
      CoordinateSystem_ConvertToRCDP(2)=x(1)*SIN(x(2))*COS(x(3))
      CoordinateSystem_ConvertToRCDP(3)=x(1)*SIN(x(3))
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      IF(coordinateSystem%numberOfDimensions/=3) THEN
        localError="The number of coordinates of "//TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))// &
          & " is invalid for a prolate spheroidal coordinate system."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      focus=coordinateSystem%focus
      CoordinateSystem_ConvertToRCDP(1)=focus*COSH(x(1))*COS(x(2))
      CoordinateSystem_ConvertToRCDP(2)=focus*SINH(x(1))*SIN(x(2))*COS(x(3))
      CoordinateSystem_ConvertToRCDP(3)=focus*SINH(x(1))*SIN(x(2))*SIN(x(3))
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      IF(coordinateSystem%numberOfDimensions/=3) THEN
        localError="The number of coordinates of "//TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))// &
          & " is invalid for a oblate spheroidal coordinate system."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      focus=coordinateSystem%focus
      CoordinateSystem_ConvertToRCDP(1)=focus*COSH(x(1))*COS(x(2))*COS(x(3))
      CoordinateSystem_ConvertToRCDP(2)=focus*SINH(x(1))*SIN(x(2))
      CoordinateSystem_ConvertToRCDP(3)=focus*COSH(x(1))*COS(x(2))*SIN(x(3))
    CASE DEFAULT
      localError="The coordinate type of "//TRIM(NumberToVString(coordinateSystem%TYPE,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("CoordinateSystem_ConvertToRCDP")
    RETURN
999 ERRORSEXITS("CoordinateSystem_ConvertToRCDP",err,error)
    RETURN
    
  END FUNCTION CoordinateSystem_ConvertToRCDP

  !
  !================================================================================================================================
  !

  !>CoordinateSystem_ConvertToRCSP performs a coordinate transformation from a coordinate system identified by
  !>coordinateSystem at the point x(:) to the returned point z(:) in rectangular cartesian coordinates for
  !>single precision coordinates.
  FUNCTION CoordinateSystem_ConvertToRCSP(coordinateSystem,x,err,error)
  
    !Argument variables
    TYPE(CoordinateSystemType), INTENT(IN) :: coordinateSystem !<The coordinate system to convert to rectangular cartesian
    REAL(SP), INTENT(IN) :: x(:) !<The coordinate to convert
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(SP) :: CoordinateSystem_ConvertToRCSP(SIZE(x,1))
    !Local variables
    REAL(SP) :: focus
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("CoordinateSystem_ConvertToRCSP",err,error,*999)
    
    CoordinateSystem_ConvertToRCSP=0.0_SP

    IF(SIZE(x,1)<coordinateSystem%numberOfDimensions) THEN
     localError="The size of x of "//TRIM(NumberToVString(SIZE(x,1),"*",err,error))// &
        & " is less than the number of coordinate system dimensions of "// &
        & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    SELECT CASE(coordinateSystem%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      CoordinateSystem_ConvertToRCSP(1:coordinateSystem%numberOfDimensions)=x(1:coordinateSystem%numberOfDimensions)
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      SELECT CASE(coordinateSystem%numberOfDimensions)
      CASE(2)
        CoordinateSystem_ConvertToRCSP(1)=x(1)*COS(x(2))
        CoordinateSystem_ConvertToRCSP(2)=x(1)*SIN(x(2))
      CASE(3)
        CoordinateSystem_ConvertToRCSP(1)=x(1)*COS(x(2))
        CoordinateSystem_ConvertToRCSP(2)=x(1)*SIN(x(2))
        CoordinateSystem_ConvertToRCSP(3)=x(3)
      CASE DEFAULT
        localError="The number of coordinates of "//TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))// &
          & " is invalid for a cylindrical polar coordinate system."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(coordinateSystem%numberOfDimensions/=3) THEN  
        localError="The number of coordinates of "//TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))// &
          & " is invalid for a spherical polar coordinate system."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      CoordinateSystem_ConvertToRCSP(1)=x(1)*COS(x(2))*COS(x(3))
      CoordinateSystem_ConvertToRCSP(2)=x(1)*SIN(x(2))*COS(x(3))
      CoordinateSystem_ConvertToRCSP(3)=x(1)*SIN(x(3))
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      IF(coordinateSystem%numberOfDimensions/=3) THEN
        localError="The number of coordinates of "//TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))// &
          & " is invalid for a prolate spheroidal coordinate system."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      focus=REAL(coordinateSystem%focus,SP)
      CoordinateSystem_ConvertToRCSP(1)=focus*COSH(x(1))*COS(x(2))
      CoordinateSystem_ConvertToRCSP(2)=focus*SINH(x(1))*SIN(x(2))*COS(x(3))
      CoordinateSystem_ConvertToRCSP(3)=focus*SINH(x(1))*SIN(x(2))*SIN(x(3))
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      IF(coordinateSystem%numberOfDimensions/=3) THEN
        localError="The number of coordinates of "//TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))// &
          & " is invalid for a oblate spheroidal coordinate system."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      focus=REAL(coordinateSystem%focus,SP)
      CoordinateSystem_ConvertToRCSP(1)=focus*COSH(x(1))*COS(x(2))*COS(x(3))
      CoordinateSystem_ConvertToRCSP(2)=focus*SINH(x(1))*SIN(x(2))
      CoordinateSystem_ConvertToRCSP(3)=focus*COSH(x(1))*COS(x(2))*SIN(x(3))
    CASE DEFAULT
      localError="The coordinate type of "//TRIM(NumberToVString(coordinateSystem%TYPE,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("CoordinateSystem_ConvertToRCSP")
    RETURN
999 ERRORSEXITS("CoordinateSystem_ConvertToRCSP",err,error)
    RETURN
    
  END FUNCTION CoordinateSystem_ConvertToRCSP

  !
  !================================================================================================================================
  !

  !>Calculates the difference (or detlta) between the point x and the point y i.e., y-x, in the given coordinate system.
  !>0->2Pi discontinuities with polar coordinates are accounted for.
  FUNCTION CoordinateSystem_DeltaCalculateDP(coordinateSystem,x,y,err,error)
  
    !Argument variables
    TYPE(CoordinateSystemType), INTENT(IN) :: coordinateSystem !<The coordinate system to calculate the delta for
    REAL(DP), INTENT(IN) :: x(:) !<The first coordinate
    REAL(DP), INTENT(IN) :: y(:) !<The second coordinate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: CoordinateSystem_DeltaCalculateDP(SIZE(x,1))
    !Local variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("CoordinateSystem_DeltaCalculateDP",err,error,*999)

    CoordinateSystem_DeltaCalculateDP=0.0_DP

    IF(SIZE(x,1)<coordinateSystem%numberOfDimensions) THEN
      localError="The size of x of "//TRIM(NumberToVString(SIZE(x,1),"*",err,error))// &
        & " is less than the number of coordinate system dimensions of "// &
        & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(x,1)/=SIZE(y,1)) THEN
      localError="The size of x of "//TRIM(NumberToVString(SIZE(x,1),"*",err,error))// &
        & " does not matach the size of y of "//TRIM(NumberToVString(SIZE(y,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
   
    CoordinateSystem_DeltaCalculateDP(1:coordinateSystem%numberOfDimensions)=y(1:coordinateSystem%numberOfDimensions)- &
      & x(1:coordinateSystem%numberOfDimensions)
    SELECT CASE(coordinateSystem%type)
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
      localError="The coordinate type of "//TRIM(NumberToVString(coordinateSystem%TYPE,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("CoordinateSystem_DeltaCalculateDP")
    RETURN
999 ERRORSEXITS("CoordinateSystem_DeltaCalculateDP",err,error)
    RETURN
    
  END FUNCTION CoordinateSystem_DeltaCalculateDP

  !
  !================================================================================================================================
  !

  !>Calculates the covariant metric tensor gl(i,j), the contravariant metric tensor gu(i,j), the Jacobian and derivative of the interpolated coordinate system (xi_i) with respect to the given coordinate (x_j) system (dXidX) at a point (x - normally a Gauss point). Old cmiss name: XGMG
  SUBROUTINE CoordinateSystem_MetricsCalculate(coordinateSystem,jacobianType,metrics,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<A pointer to the coordinate system to calculate the metrics for
    INTEGER(INTG), INTENT(IN) :: jacobianType !<The type of Jacobian to calculate \see CoordinateRoutines_JacobianTypes,CoordinateRoutines
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: metrics !<A pointer to the metrics to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: xiIdx1,xiIdx2,partialDerivativeIndex
    REAL(DP) :: detGL,detdXdXi,dXdXi2(3),dXdXi3(3),ff,g1,g3,length,mu,r,rc,rcrc,rr,scale
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint
    TYPE(VARYING_STRING) :: localError

     ENTERS("CoordinateSystem_MetricsCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(metrics)) CALL FlagError("Metrics is not associated.",err,error,*999)

    NULLIFY(interpolatedPoint)
    CALL FieldInterpolatedPointMetrics_InterpolatedPointGet(metrics,interpolatedPoint,err,error,*999)
    
    IF(interpolatedPoint%partialDerivativeType<FIRST_PART_DERIV) THEN
      CALL FlagError("Metrics interpolated point has not been interpolated to include first derivatives.",err,error,*999)
    ENDIF

    SELECT CASE(metrics%numberOfXiDimensions)
    CASE(1)
      !Calculate the derivatives of x with respect to xi
      partialDerivativeIndex=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(1)
      metrics%dXdXi(1:metrics%numberOfXDimensions,1)= &
        & interpolatedPoint%values(1:metrics%numberOfXDimensions,partialDerivativeIndex)
      !Initialise the covariant metric tensor to the identity matrix
      metrics%gl(1,1)=1.0_DP
    CASE(2)
      !Calculate the derivatives of x with respect to xi
      partialDerivativeIndex=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(1)
      metrics%dXdXi(1:metrics%numberOfXDimensions,1)= &
        & interpolatedPoint%values(1:metrics%numberOfXDimensions,partialDerivativeIndex)
      partialDerivativeIndex=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(2)
      metrics%dXdXi(1:metrics%numberOfXDimensions,2)= &
        & interpolatedPoint%values(1:metrics%numberOfXDimensions,partialDerivativeIndex)
      !Initialise the covariant metric tensor to the identity matrix
      metrics%gl(1,1)=1.0_DP
      metrics%gl(1,2)=0.0_DP
      metrics%gl(2,1)=0.0_DP
      metrics%gl(2,2)=1.0_DP
    CASE(3)
      !Calculate the derivatives of x with respect to xi
      partialDerivativeIndex=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(1)
      metrics%dXdXi(1:metrics%numberOfXDimensions,1)= &
        & interpolatedPoint%values(1:metrics%numberOfXDimensions,partialDerivativeIndex)
      partialDerivativeIndex=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(2)
      metrics%dXdXi(1:metrics%numberOfXDimensions,2)= &
        & interpolatedPoint%values(1:metrics%numberOfXDimensions,partialDerivativeIndex)
      partialDerivativeIndex=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(3)
      metrics%dXdXi(1:metrics%numberOfXDimensions,3)= &
        & interpolatedPoint%values(1:metrics%numberOfXDimensions,partialDerivativeIndex)
      !Initialise the covariant metric tensor to the identity matrix
      metrics%gl(1,1)=1.0_DP
      metrics%gl(1,2)=0.0_DP
      metrics%gl(1,3)=0.0_DP
      metrics%gl(2,1)=0.0_DP
      metrics%gl(2,2)=1.0_DP
      metrics%gl(2,3)=0.0_DP
      metrics%gl(3,1)=0.0_DP
      metrics%gl(3,2)=0.0_DP
      metrics%gl(3,3)=1.0_DP
    CASE DEFAULT
      localError="The number of coordinates of "//TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    !Calculate the covariant metric tensor gl(i,j)
    SELECT CASE(coordinateSystem%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      SELECT CASE(metrics%numberOfXDimensions)
      CASE(1)
        DO xiIdx1=1,metrics%numberOfXiDimensions
          DO xiIdx2=1,metrics%numberOfXiDimensions
            metrics%gl(xiIdx1,xiIdx2)=metrics%dXdXi(1,xiIdx1)*metrics%dXdXi(1,xiIdx2)                    
          ENDDO !xiIdx2
        ENDDO !xiIdx1
      CASE(2)
        DO xiIdx1=1,metrics%numberOfXiDimensions
          DO xiIdx2=1,metrics%numberOfXiDimensions
            metrics%gl(xiIdx1,xiIdx2)=metrics%dXdXi(1,xiIdx1)*metrics%dXdXi(1,xiIdx2)+ &
              & metrics%dXdXi(2,xiIdx1)*metrics%dXdXi(2,xiIdx2)
          ENDDO !xiIdx2
        ENDDO !xiIdx1
      CASE(3)
        DO xiIdx1=1,metrics%numberOfXiDimensions
          DO xiIdx2=1,metrics%numberOfXiDimensions
            metrics%gl(xiIdx1,xiIdx2)=metrics%dXdXi(1,xiIdx1)*metrics%dXdXi(1,xiIdx2)+ &
              & metrics%dXdXi(2,xiIdx1)*metrics%dXdXi(2,xiIdx2)+ &
              & metrics%dXdXi(3,xiIdx1)*metrics%dXdXi(3,xiIdx2)
          ENDDO !xiIdx2
        ENDDO !xiIdx1
      CASE DEFAULT
        localError="The number of coordinates of "//TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))// &
          & " is invalid for a rectangular cartesian coordinate system."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      r=interpolatedPoint%values(1,1)
      rr=r*r
      IF(metrics%numberOfXDimensions==2) THEN
        DO xiIdx1=1,metrics%numberOfXiDimensions
          DO xiIdx2=1,metrics%numberOfXiDimensions
            metrics%gl(xiIdx1,xiIdx2)= &
              & metrics%dXdXi(1,xiIdx1)*metrics%dXdXi(1,xiIdx2)+ &
              & rr*metrics%dXdXi(2,xiIdx1)*metrics%dXdXi(2,xiIdx2)
          ENDDO !xiIdx2
        ENDDO !xiIdx1
      ELSE IF(metrics%numberOfXDimensions==3) THEN
        DO xiIdx1=1,metrics%numberOfXiDimensions
          DO xiIdx2=1,metrics%numberOfXiDimensions
            metrics%gl(xiIdx1,xiIdx2)= &
              & metrics%dXdXi(1,xiIdx1)*metrics%dXdXi(1,xiIdx2)+ &
              & rr*metrics%dXdXi(2,xiIdx1)*metrics%dXdXi(2,xiIdx2)+ &
              & metrics%dXdXi(3,xiIdx1)*metrics%dXdXi(3,xiIdx2)
          ENDDO !xiIdx2
        ENDDO !xiIdx1
      ELSE
        localError="The number of coordinates of "//TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))// &
          & " is invalid for a cylindrical polar coordinate system."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      r=interpolatedPoint%values(1,1)
      rr=r*r
      rc=r*COS(interpolatedPoint%values(3,1))
      rcrc=rc*rc          
      DO xiIdx1=1,metrics%numberOfXiDimensions
        DO xiIdx2=1,metrics%numberOfXiDimensions
          metrics%gl(xiIdx1,xiIdx2)=metrics%dXdXi(1,xiIdx1)*metrics%dXdXi(1,xiIdx2)+ &
            & rcrc*metrics%dXdXi(2,xiIdx1)*metrics%dXdXi(2,xiIdx2)+ &
            & rr*metrics%dXdXi(3,xiIdx1)*metrics%dXdXi(3,xiIdx2)
        ENDDO !xiIdx2
      ENDDO !xiIdx1
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      IF(ABS(interpolatedPoint%values(2,1))<ZERO_TOLERANCE) THEN
        CALL FlagWarning("Mu is zero.",err,error,*999)
      ELSE
        ff=coordinateSystem%focus*coordinateSystem%focus
        r=interpolatedPoint%values(1,1)
        mu=interpolatedPoint%values(2,1)
        g1=ff*(SINH(r)*SINH(r)+SIN(mu)*SIN(mu))
        g3=ff*SINH(r)*SINH(r)*SIN(mu)*SIN(mu)
        IF(metrics%numberOfXDimensions==2) THEN
          DO xiIdx1=1,metrics%numberOfXiDimensions
            DO xiIdx2=1,metrics%numberOfXiDimensions
              metrics%gl(xiIdx1,xiIdx2)=g1*(metrics%dXdXi(1,xiIdx1)*metrics%dXdXi(1,xiIdx2)+ &
                & metrics%dXdXi(2,xiIdx1)*metrics%dXdXi(2,xiIdx2))
            ENDDO !xiIdx2
          ENDDO !xiIdx1
        ELSE IF(metrics%numberOfXDimensions==3) THEN
          DO xiIdx1=1,metrics%numberOfXiDimensions
            DO xiIdx2=1,metrics%numberOfXiDimensions
              metrics%gl(xiIdx1,xiIdx2)=g1*(metrics%dXdXi(1,xiIdx1)*metrics%dXdXi(1,xiIdx2)+ &
                & metrics%dXdXi(2,xiIdx1)*metrics%dXdXi(2,xiIdx2))+ &
                & g3*metrics%dXdXi(3,xiIdx1)*metrics%dXdXi(3,xiIdx2)
            ENDDO !xiIdx2
          ENDDO !xiIdx1
        ELSE
          localError="The number of coordinates of "//TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))// &
            & " is invalid for a prolate spheroidal coordinate system."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The coordinate type of "//TRIM(NumberToVString(coordinateSystem%type,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    !Calcualte the contravariant metric tensor
    CALL Invert(metrics%gl(1:metrics%numberOfXiDimensions,1:metrics%numberOfXiDimensions), &
      & metrics%gu(1:metrics%numberOfXiDimensions,1:metrics%numberOfXiDimensions),detGL, &
      & err,error,*999)
    
    !Calculate the Jacobian
    SELECT CASE(jacobianType)
    CASE(COORDINATE_JACOBIAN_NO_TYPE)
      metrics%jacobian=0.0
      metrics%jacobianType=COORDINATE_JACOBIAN_NO_TYPE
    CASE(COORDINATE_JACOBIAN_LINE_TYPE)
      metrics%jacobian=SQRT(ABS(metrics%gl(1,1)))
      metrics%jacobianType=COORDINATE_JACOBIAN_LINE_TYPE
    CASE(COORDINATE_JACOBIAN_AREA_TYPE)
      IF(metrics%numberOfXiDimensions==3) THEN
        metrics%jacobian=SQRT(ABS(detGL*metrics%gu(3,3)))
      ELSE
        metrics%jacobian=SQRT(ABS(detGL))
      ENDIF
      metrics%jacobianType=COORDINATE_JACOBIAN_AREA_TYPE
    CASE(COORDINATE_JACOBIAN_VOLUME_TYPE)
      metrics%jacobian=SQRT(ABS(detGL))
      metrics%jacobianType=COORDINATE_JACOBIAN_VOLUME_TYPE
    CASE DEFAULT
      localError="The Jacobian type of "//TRIM(NumberToVString(jacobianType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    !Calculate the derivatives of xi with respect to x - dXidX
    IF(metrics%numberOfXiDimensions==metrics%numberOfXDimensions) THEN
      CALL Invert(metrics%dXdXi,metrics%dXidX,detdXdXi,err,error,*999)
    ELSE
      !We have a line or a surface embedded in a higher dimensional space
      SELECT CASE(metrics%numberOfXiDimensions)
      CASE(1)
        !Line in space
        SELECT CASE(metrics%numberOfXDimensions)
        CASE(2)
          IF(interpolatedPoint%partialDerivativeType>FIRST_PART_DERIV) THEN
            !We have curvature information. Form the frenet vector frame.
            !Calculate the normal vector from the normalised second derivative of the position vector.
            partialDerivativeIndex=PARTIAL_DERIVATIVE_SECOND_DERIVATIVE_MAP(1)
            CALL Normalise(interpolatedPoint%values(1:2,partialDerivativeIndex),dXdXi2,err,error,*999)
          ELSE
            !No curvature information but obtain other normal frenet vector by rotating tangent vector 90 deg.
            dXdXi2(1)=-1.0_DP*metrics%dXdXi(2,1)
            dXdXi2(2)=metrics%dXdXi(1,1)                    
          ENDIF
          detdXdXi=metrics%dXdXi(1,1)*dXdXi2(2)-metrics%dXdXi(2,1)*dXdXi2(1)
          IF(ABS(detdXdXi)>ZERO_TOLERANCE) THEN
            metrics%dXidX(1,1)=dXdXi2(2)/detdXdXi
            metrics%dXidX(1,2)=-1.0_DP*dXdXi2(1)/detdXdXi
            !Normalise to ensure that g^11=g^1.g^1
            CALL L2Norm(metrics%dXidX(1,1:2),length,err,error,*999)
            scale=SQRT(ABS(metrics%gu(1,1)))/length
            metrics%dXidX(1,1:2)=scale*metrics%dXidX(1,1:2)
          ELSE
            CALL FlagWarning("Zero determinant. Unable to obtain dxi/dx.",err,error,*999)
            metrics%dXidX=0.0_DP                    
          ENDIF
        CASE(3)
          IF(interpolatedPoint%partialDerivativeType>FIRST_PART_DERIV) THEN
            !We have curvature information. Form the frenet vector frame.
            !Calculate the normal vector from the normalised second derivative of the position vector.
            partialDerivativeIndex=PARTIAL_DERIVATIVE_SECOND_DERIVATIVE_MAP(1)
            CALL Normalise(interpolatedPoint%values(1:3,partialDerivativeIndex),dXdXi2,err,error,*999)
            !Calculate the bi-normal vector from the normalised cross product of the tangent and normal vectors
            CALL NormaliseCrossProduct(metrics%dXdXi(1:3,1),dXdXi2,dXdXi3,err,error,*999)
            detdXdXi=metrics%dXdXi(1,1)*(dXdXi2(2)*dXdXi3(3)-dXdXi2(3)*dXdXi3(2))+ &
              & dXdXi2(1)*(metrics%dXdXi(3,1)*dXdXi3(2)-dXdXi3(3)*metrics%dXdXi(2,1))+ &
              & dXdXi3(1)*(metrics%dXdXi(2,1)*dXdXi2(3)-metrics%dXdXi(3,1)*dXdXi2(2))
            IF(ABS(detdXdXi)>ZERO_TOLERANCE) THEN
              metrics%dXidX(1,1)=(dXdXi3(3)*dXdXi2(2)-dXdXi2(3)*dXdXi3(2))/detdXdXi
              metrics%dXidX(1,2)=-1.0_DP*(dXdXi3(3)*dXdXi2(1)-dXdXi2(3)*dXdXi3(1))/detdXdXi
              metrics%dXidX(1,3)=(dXdXi3(2)*dXdXi2(1)-dXdXi2(2)*dXdXi3(1))/detdXdXi
              !Normalise to ensure that g^11=g^1.g^1
              CALL L2Norm(metrics%dXidX(1,1:3),length,err,error,*999)
              scale=SQRT(ABS(metrics%gu(1,1)))/length
              metrics%dXidX(1,1:3)=scale*metrics%dXdXi(1,1:3)
            ELSE
              CALL FlagWarning("Zero determinant. Unable to obtain dxi/dx.",err,error,*999)
              metrics%dXidX=0.0_DP                    
            ENDIF
          ELSE
            metrics%dXidX(1,1)=metrics%dXdXi(1,1)
            metrics%dXidX(1,2)=metrics%dXdXi(2,1)
            metrics%dXidX(1,3)=metrics%dXdXi(3,1)
            !Normalise to ensure that g^11=g^1.g^1
            CALL L2Norm(metrics%dXidX(1,1:3),length,err,error,*999)
            scale=SQRT(ABS(metrics%gu(1,1)))/length
            metrics%dXidX(1,1:3)=scale*metrics%dXidX(1,1:3)
          ENDIF
        CASE DEFAULT
          CALL FlagError("Invalid embedding of a line in space.",err,error,*999)
        END SELECT
      CASE(2)
        !Surface in space
        IF(metrics%numberOfXDimensions==3) THEN
          !Surface in 3D space.
          !Calculate the covariant vectors g^1 and g^2. These are calculated as follows:
          !First define g_3=g_1 x g_2, and then define g^1=((g_2 x g_3)_b)/detGL and g^2=((g_3 x g_1)_b)/detGL. 
          !The _b means lowering the index with the metric tensor of the curvilinear coordinate system.
          !This way we have a consistent set of covariant and covariant vectors, i.e.  <g_M,g^N>=delta_M^N.
          metrics%dXidX(1,1:3)=(metrics%gl(2,2)*metrics%dXdXi(1:3,1)-metrics%gl(1,2)*metrics%dXdXi(1:3,2))/detGL
          metrics%dXidX(2,1:3)=(metrics%gl(1,1)*metrics%dXdXi(1:3,2)-metrics%gl(2,1)*metrics%dXdXi(1:3,1))/detGL
          SELECT CASE(coordinateSystem%type)
          CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
            !Do nothing
          CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
            r=interpolatedPoint%values(1,1)
            rr=r*r
            metrics%dXidX(1:2,2)=metrics%dXidX(1:2,2)*rr
          CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
            r=interpolatedPoint%values(1,1)
            rr=r*r
            rc=R*COS(interpolatedPoint%values(3,1))
            rcrc=rc*rc          
            metrics%dXidX(1:2,2)=metrics%dXidX(1:2,2)*rcrc
            metrics%dXidX(1:2,3)=metrics%dXidX(1:2,3)*rr
          CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
            IF(ABS(interpolatedPoint%values(2,1))<ZERO_TOLERANCE) THEN
              CALL FlagWarning("Mu is zero.",err,error,*999)
            ELSE
              ff=coordinateSystem%focus*coordinateSystem%focus
              r=interpolatedPoint%values(1,1)
              mu=interpolatedPoint%values(2,1)
              g1=ff*(SINH(r)*SINH(r)+SIN(mu)*SIN(mu))
              g3=ff*SINH(r)*SINH(r)*SIN(mu)*SIN(mu)
              metrics%dXidX(1:2,1)=metrics%dXidX(1:2,1)*g1
              metrics%dXidX(1:2,2)=metrics%dXidX(1:2,2)*g1
              metrics%dXidX(1:2,3)=metrics%dXidX(1:2,3)*g3
            ENDIF
          CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The coordinate system type of "//TRIM(NumberToVString(coordinateSystem%type,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Invalid embedding of a surface in space.",err,error,*999)
        ENDIF
      CASE DEFAULT
        CALL FlagError("Invalid embedding in space.",err,error,*999)
      END SELECT
    ENDIF
  
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Coordinate system metrics:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Coordinate system type = ",TRIM(COORDINATE_SYSTEM_TYPE_STRING( &
        & coordinateSystem%type)),err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of X dimensions = ",metrics%numberOfXDimensions,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of Xi dimensions = ",metrics%numberOfXiDimensions,err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Location of metrics:",err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,metrics%numberOfXDimensions,3,3,interpolatedPoint%values(:,1), &
        & '("    X           :",3(X,E13.6))','(17X,3(X,E13.6))',err,error,*999)      
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Derivative of X wrt Xi:",err,error,*999)
      CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,metrics%numberOfXDimensions,1,1,metrics%numberOfXiDimensions, &
        & 3,3,metrics%dXdXi,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    dX_dXi','(",I1,",:)',' :",3(X,E13.6))', &
        & '(17X,3(X,E13.6))',err,error,*999)
      IF(metrics%numberOfXDimensions/=metrics%numberOfXiDimensions) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Constructed derivative of X wrt Xi:",err,error,*999)
        SELECT CASE(metrics%numberOfXiDimensions)
        CASE(1)
          !Line in space
          SELECT CASE(metrics%numberOfXDimensions)
          CASE(2)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,metrics%numberOfXDimensions,3,3,dXdXi2, &
              & '("    dX_dXi(:,2) :",3(X,E13.6))','(17X,3(X,E13.6))',err,error,*999)      
          CASE(3)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,metrics%numberOfXDimensions,3,3,dXdXi2, &
              & '("    dX_dXi(:,2) :",3(X,E13.6))','(17X,3(X,E13.6))',err,error,*999)      
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,metrics%numberOfXDimensions,3,3,dXdXi3, &
              & '("    dX_dXi(:,3) :",3(X,E13.6))','(17X,3(X,E13.6))',err,error,*999)      
          CASE DEFAULT
            CALL FlagError("Invalid embedding of a line in space.",err,error,*999)
          END SELECT
        CASE(2)
          !Surface in space
          SELECT CASE(metrics%numberOfXDimensions)
          CASE(3)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,metrics%numberOfXDimensions,3,3,dXdXi3, &
              & '("    dX_dXi(:,3) :",3(X,E13.6))','(17X,3(X,E13.6))',err,error,*999)      
          CASE DEFAULT
            CALL FlagError("Invalid embedding of a surface in space.",err,error,*999)
          END SELECT
        CASE DEFAULT
          CALL FlagError("Invalid embedding in space.",err,error,*999)
        END SELECT
      ENDIF
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  det dX_dXi    = ",detdXdXi,err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Derivative of Xi wrt X:",err,error,*999)
      CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,metrics%numberOfXiDimensions,1,1,metrics%numberOfXDimensions, &
        & 3,3,metrics%dXidX,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    dXi_dX','(",I1,",:)',' :",3(X,E13.6))', &
        & '(17X,3(X,E13.6))',err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Covariant metric tensor:",err,error,*999)
      CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,metrics%numberOfXiDimensions,1,1,metrics%numberOfXiDimensions, &
        & 3,3,metrics%gl,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    GL','(",I1,",:)','     :",3(X,E13.6))','(17X,3(X,E13.6))', &
        & err,error,*999)      
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Contravariant metric tensor:",err,error,*999)
      CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,metrics%numberOfXiDimensions,1,1,metrics%numberOfXiDimensions, &
        & 3,3,metrics%gu,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    GU','(",I1,",:)','     :",3(X,E13.6))','(17X,3(X,E13.6))', &
        & err,error,*999)      
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Jacobian type = ",metrics%jacobianType,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Jacobian      = ",metrics%jacobian,err,error,*999)
    ENDIF
    
    EXITS("CoordinateSystem_MetricsCalculate")
    RETURN
999 ERRORSEXITS("CoordinateSystem_MetricsCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_MetricsCalculate

  !
  !================================================================================================================================
  !

  !>Calculates the normal vector, n, at the point x. IF reverse is true the reversed normal is returned. Old-cmiss-name: NORMAL
  SUBROUTINE CoordinateSystem_NormalCalculate(coordinateSystem,reverse,x,n,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<The coordinate system to calculate the normal for
    LOGICAL, INTENT(IN) :: reverse !<If .TRUE. the reversed normal is returned.
    REAL(DP), INTENT(IN) :: x(:,:) !<The coordinate and it's derivatives to calcualte the normal at
    REAL(DP), INTENT(OUT) :: n(:) !<On exit, the normal vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfXDimensions,dS1,dS2,d2S1
    REAL(DP) :: length,r,tangent1(3),tangent2(3)
    TYPE(VARYING_STRING) :: localError

    ENTERS("CoordinateSystem_NormalCalculate",err,error,*999)

    CALL CoordinateSystem_AssertIsFinished(coordinateSystem,err,error,*999)
    
    numberOfXDimensions=coordinateSystem%numberOfDimensions    
    dS1=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(1)
    dS2=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(2)
    d2S1=PARTIAL_DERIVATIVE_SECOND_DERIVATIVE_MAP(1)
    SELECT CASE(coordinateSystem%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      IF(numberOfXDimensions==2) THEN
        tangent1(1)=x(1,dS1)
        tangent1(2)=x(2,dS1)
      ELSE IF(numberOfXDimensions==3) THEN
        tangent1(1)=x(1,dS1)
        tangent1(2)=x(2,dS1)
        tangent1(3)=x(3,dS1)
        tangent2(1)=x(1,dS2)
        tangent2(2)=x(2,dS2)
        tangent2(3)=x(3,dS2)
      ELSE
        localError="The number of coordinates of "//TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))// &
          & " is invalid to calculate a normal in a rectangular cartesian coordinate system."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      r=x(1,1)
      IF(numberOfXDimensions==2) THEN
        tangent1(1)=x(1,dS1)*COS(x(1,dS1))-r*SIN(x(1,dS1))*x(2,dS1)
        tangent1(2)=x(2,dS1)*SIN(x(1,dS1))+r*COS(x(1,dS1))*x(2,dS1)
      ELSE IF(numberOfXDimensions==3) THEN
        tangent1(1)=x(1,dS1)*COS(x(1,dS1))-r*SIN(x(1,dS1))*x(2,dS1)
        tangent1(2)=x(2,dS1)*SIN(x(1,dS1))+r*COS(x(1,dS1))*x(2,dS1)
        tangent1(3)=x(3,dS1)
        tangent2(1)=x(1,dS2)*COS(x(1,dS1))-r*SIN(x(1,dS1))*x(2,dS2)
        tangent2(2)=x(1,dS2)*SIN(x(1,dS1))+r*COS(x(1,dS1))*x(2,dS2)
        tangent2(3)=x(3,dS2)
      ELSE
        localError="The number of coordinates of "//TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))// &
          & " is invalid to calculate a normal in a cylindrical polar coordinate system."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      r=x(1,1)
      tangent1(1)=x(1,dS1)*COS(x(1,d2S1))*COS(x(1,dS1))- &
        &                 r*SIN(x(1,d2S1))*COS(x(1,dS1))*x(3,dS1)- &
        &                 r*COS(x(1,d2S1))*SIN(x(1,dS1))*x(2,dS1)
      tangent1(2)=x(1,dS1)*COS(x(1,d2S1))*SIN(x(1,dS1))- &
        &                 r*SIN(x(1,d2S1))*SIN(x(1,dS1))*x(3,dS1)+ &
        &                 r*COS(x(1,d2S1))*COS(x(1,dS1))*x(2,dS1)
      tangent1(3)=x(1,dS1)*SIN(x(1,d2S1))+R*COS(x(1,d2S1))*x(3,dS1)
      tangent2(1)=x(1,dS2)*COS(x(1,d2S1))*COS(x(1,dS1))- &
        &                 r*SIN(x(1,d2S1))*COS(x(1,dS1))*x(3,dS2)- &
        &                 r*COS(x(1,d2S1))*SIN(x(1,dS1))*x(2,dS2)
      tangent2(2)=x(1,dS2)*COS(x(1,d2S1))*SIN(x(1,dS1))- &
        &                 r*SIN(x(1,d2S1))*SIN(x(1,dS1))*x(3,dS2)+ &
        &                 r*COS(x(1,d2S1))*COS(x(1,dS1))*x(2,dS2)
      tangent2(3)=x(1,dS2)*SIN(x(1,d2S1))+r*COS(x(1,d2S1))*x(3,dS2)
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The coordinate system type of "//TRIM(NumberToVString(coordinateSystem%type,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    IF(numberOfXDimensions==2) THEN
      n(1)=-tangent1(2)
      n(2)=tangent1(1)
      length=SQRT(n(1)*n(1)+n(2)*n(2))
      IF(ABS(length)<ZERO_TOLERANCE) CALL FlagError("Zero normal vector length.",err,error,*999)
      IF(reverse) THEN
        n(1)=-n(1)/length
        n(2)=-n(2)/length
      ELSE            
        n(1)=n(1)/length
        n(2)=n(2)/length
      ENDIF
    ELSE
      n(1)=tangent1(2)*tangent2(3)-tangent1(3)*tangent2(2)
      n(2)=tangent1(3)*tangent2(1)-tangent1(1)*tangent2(3)
      n(3)=tangent1(1)*tangent2(2)-tangent1(2)*tangent2(1)
      length=SQRT(n(1)*n(1)+n(2)*n(2)+n(3)*n(3))
      IF(reverse) THEN
        n(1)=-n(1)/length
        n(2)=-n(2)/length
        n(3)=-n(3)/length
      ELSE            
        n(1)=n(1)/length
        n(2)=n(2)/length
        n(3)=n(3)/length
      ENDIF
    ENDIF

    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Coordinate system metrics:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Coordinate system type = ",TRIM(COORDINATE_SYSTEM_TYPE_STRING( &
        & coordinateSystem%type)),err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of X dimensions = ",numberOfXDimensions,err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfXDimensions,3,3,x(:,1),'("  X         :",3(X,E13.6))', &
        & '(13X,3(X,E13.6))',err,error,*999)      
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfXDimensions,3,3,tangent1,'("  Tangent 1 :",3(X,E13.6))', &
        & '(13X,3(X,E13.6))',err,error,*999)
      IF(numberOfXDimensions==3) THEN
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfXDimensions,3,3,tangent2,'("  Tangent 2 :",3(X,E13.6))', &
          & '(13X,3(X,E13.6))',err,error,*999)
      ENDIF
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfXDimensions,3,3,N,'("  Normal    :",3(X,E13.6))', &
        & '(13X,3(X,E13.6))',err,error,*999)            
    ENDIF
    
    EXITS("CoordinateSystem_NormalCalculate")
    RETURN
999 ERRORSEXITS("CoordinateSystem_NormalCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_NormalCalculate

   !
  !================================================================================================================================
  !

  !>Finalises a coordinate system and deallocates all memory. 
  SUBROUTINE CoordinateSystem_Finalise(coordinateSystem,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<A pointer to the coordinate system to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CoordinateSystem_Finalise",err,error,*999)

    IF(ASSOCIATED(coordinateSystem)) THEN
      DEALLOCATE(coordinateSystem)
    ENDIF
   
    EXITS("CoordinateSystem_Finalise")
    RETURN
999 ERRORSEXITS("CoordinateSystem_Finalise",err,error)
    RETURN 1

  END SUBROUTINE CoordinateSystem_Finalise

  !
  !================================================================================================================================
  !

  !>Sets/changes the dimension of the coordinate system. \see OpenCMISS::Iron::cmfe_CoordinateSystem_DimensionSet
  SUBROUTINE CoordinateSystem_DimensionSet(coordinateSystem,dimension,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<A pointer the coordinate system to set the dimension for
    INTEGER(INTG), INTENT(IN) :: dimension !<The dimension to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("CoordinateSystem_DimensionSet",err,error,*999)

    CALL CoordinateSystem_AssertNotFinished(coordinateSystem,err,error,*999)
    
    SELECT CASE(coordinateSystem%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      IF(dimension<1.OR.dimension>3) THEN
        localError="The specified dimension of "//TRIM(NumberToVString(dimension,"*",err,error))// &
          & " is invalid for a rectangular cartesian coordinate system."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      coordinateSystem%numberOfDimensions=dimension
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      IF(dimension<2.OR.dimension>3) THEN
        localError="The specified dimension of "//TRIM(NumberToVString(dimension,"*",err,error))// &
          & " is invalid for a cylindrical polar coordinate system."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      coordinateSystem%numberOfDimensions=dimension
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(dimension/=3) THEN
        localError="The specified dimension of "//TRIM(NumberToVString(dimension,"*",err,error))// &
          & " is invalid for a spherical polar coordinate system."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      coordinateSystem%numberOfDimensions=dimension
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      IF(dimension/=3) THEN
        localError="The specified dimension of "//TRIM(NumberToVString(dimension,"*",err,error))// &
          & " is invalid for a prolate spheriodal coordinate system."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      coordinateSystem%numberOfDimensions=dimension
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      IF(dimension/=3) THEN
        localError="The specified dimension of "//TRIM(NumberToVString(dimension,"*",err,error))// &
          & " is invalid for a oblate spheriodal coordinate system."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      coordinateSystem%numberOfDimensions=dimension
    CASE DEFAULT
      localError="The coordinate type of "//TRIM(NumberToVString(coordinateSystem%type,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("CoordinateSystem_DimensionSet")
    RETURN
999 ERRORSEXITS("CoordinateSystem_DimensionSet",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_DimensionSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the focus of a coordinate system. \see OpenCMISS::Iron::cmfe_CoordinateSystem_FocusSet
  SUBROUTINE CoordinateSystem_FocusSet(coordinateSystem,focus,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<A pointer to the coordinate system to set the focus for
    REAL(DP), INTENT(IN) :: focus !<The focus to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("CoordinateSystem_FocusSet",err,error,*999)

    CALL CoordinateSystem_AssertNotFinished(coordinateSystem,err,error,*999)
    IF(focus<=ZERO_TOLERANCE) THEN
      localError="The specified focus of "//TRIM(NumberToVString(focus,"*",err,error))// &
        & " is invalid. The focus must be greater than zero."
      CALL FlagError(localError,err,error,*999)
    ENDIF
   
    SELECT CASE(coordinateSystem%type)
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE,COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      coordinateSystem%focus=focus
    CASE DEFAULT
      localError="The focus can not be defined for coordinate system type "// &
        & TRIM(NumberToVString(coordinateSystem%type,"*",err,error))
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    EXITS("CoordinateSystem_FocusSet")
    RETURN
999 ERRORSEXITS("CoordinateSystem_FocusSet",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_FocusSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the radial interpolation type of a coordinate system. \see OpenCMISS::Iron::cmfe_CoordinateSystem_RadialInterpolationTypeSet
  SUBROUTINE CoordinateSystem_RadialInterpolationTypeSet(coordinateSystem,radialInterpolationType,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<The coordinate system to set the interpolation type for
    INTEGER(INTG), INTENT(IN) :: radialInterpolationType !<The interpolation type to set \see CoordinateRoutines_RadialInterpolations,CoordinateRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("CoordinateSystem_RadialInterpolationTypeSet",err,error,*999)

    CALL CoordinateSystem_AssertNotFinished(coordinateSystem,err,error,*999)

    SELECT CASE(coordinateSystem%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      SELECT CASE(radialInterpolationType)
      CASE(COORDINATE_NO_RADIAL_INTERPOLATION_TYPE)
        coordinateSystem%radialInterpolationType=COORDINATE_NO_RADIAL_INTERPOLATION_TYPE
      CASE DEFAULT
        localError="The radial interpolation type of "//TRIM(NumberToVString(radialInterpolationType,"*",err,error))// &
          & " is invalid for a rectangular cartesian coordinate system."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE,COORDINATE_SPHERICAL_POLAR_TYPE)
      SELECT CASE(radialInterpolationType)
      CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
        coordinateSystem%radialInterpolationType=COORDINATE_RADIAL_INTERPOLATION_TYPE
      CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
        coordinateSystem%radialInterpolationType=COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE
      CASE DEFAULT
        localError="The radial interpolation type of "//TRIM(NumberToVString(radialInterpolationType,"*",err,error))// &
          & " is invalid for a cylindrical/spherical coordinate system."
        CALL FlagError(localError,err,error,*999)
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
        localError="The radial interpolation type of "//TRIM(NumberToVString(radialInterpolationType,"*",err,error))// &
          & " is invalid for a prolate spheroidal coordinate system."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The coordinate type of "//TRIM(NumberToVString(coordinateSystem%type,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("CoordinateSystem_RadialInterpolationTypeSet")
    RETURN
999 ERRORSEXITS("CoordinateSystem_RadialInterpolationTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_RadialInterpolationTypeSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of a coordinate system. \see OpenCMISS::Iron::cmfe_CoordinateSystem_TypeSet
  SUBROUTINE CoordinateSystem_TypeSet(coordinateSystem,systemType,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<The coordinate system to set the type for
    INTEGER(INTG), INTENT(IN) :: systemType !<The coordinate system type to set \see CoordinateRoutines_CoordinateSystemTypes,CoordinateRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("CoordinateSystem_TypeSet",err,error,*999)

    CALL CoordinateSystem_AssertNotFinished(coordinateSystem,err,error,*999)
    
    SELECT CASE(systemType)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      coordinateSystem%type=COORDINATE_RECTANGULAR_CARTESIAN_TYPE
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      coordinateSystem%type=COORDINATE_CYLINDRICAL_POLAR_TYPE
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      coordinateSystem%type=COORDINATE_SPHERICAL_POLAR_TYPE
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      coordinateSystem%type=COORDINATE_PROLATE_SPHEROIDAL_TYPE
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      coordinateSystem%type=COORDINATE_OBLATE_SPHEROIDAL_TYPE
    CASE DEFAULT
      localError="The specified coordinate system type of "// &
        & TRIM(NumberToVString(systemType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("CoordinateSystem_TypeSet")
    RETURN
999 ERRORSEXITS("CoordinateSystem_TypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_TypeSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the origin of a coordinate system. \see OpenCMISS::Iron::cmfe_CoordinateSystem_OriginSet
  SUBROUTINE CoordinateSystem_OriginSet(coordinateSystem,origin,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<A pointer to the coordinate system to set the origin for
    REAL(DP), INTENT(IN) :: origin(:) !<The origin to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("CoordinateSystem_OriginSet",err,error,*999)

    CALL CoordinateSystem_AssertNotFinished(coordinateSystem,err,error,*999)
    IF(SIZE(origin,1)/=3) THEN
      localError="The size of the specified orgin array of "//TRIM(NumberToVString(SIZE(origin,1),"*",err,error))// &
        & " is invalid. The size must be 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    coordinateSystem%origin=origin
    
    EXITS("CoordinateSystem_OriginSet")
    RETURN
999 ERRORSEXITS("CoordinateSystem_OriginSet",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_OriginSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the orientation of a coordinate system. \see OpenCMISS::Iron::cmfe_CoordinateSystem_OrientationSet
  SUBROUTINE CoordinateSystem_OrientationSet(coordinateSystem,orientation,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<A pointer to the coordinate system to set the orientation for
    REAL(DP), INTENT(IN) :: orientation(:,:) !<The orientation to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("CoordinateSystem_OrientationSet",err,error,*999)

    CALL CoordinateSystem_AssertNotFinished(coordinateSystem,err,error,*999)    
    IF(SIZE(orientation,1)/=3.OR.SIZE(orientation,2)/=3) THEN
      localError="The size of the specified orientation array is "//TRIM(NumberToVString(SIZE(orientation,1),"*",err,error))// &
        & "x"//TRIM(NumberToVString(SIZE(orientation,2),"*",err,error))//" and it must be 3x3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
!!TODO: \todo Check orientation matrix vectors are orthogonal to each other etc.
    coordinateSystem%orientation=orientation
    
    EXITS("CoordinateSystem_OrientationSet")
    RETURN
999 ERRORSEXITS("CoordinateSystem_OrientationSet",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_OrientationSet
  
  !
  !================================================================================================================================
  !

  !>Starts the creation of and initialises a new coordinate system. \see OpenCMISS::Iron::cmfe_CoordinateSystem_CreateStart
  !>The default values of the coordinateSystem's attributes are:
  !>- TYPE: 1 (COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
  !>- radialInterpolationType: 0 (COORDINATE_NO_RADIAL_INTERPOLATION_TYPE)
  !>- Dimensions: 3
  !>- Focus: 1.0
  !>- Origin: (0.0,0.0,0.0)
  !>- Oritention: ((1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0))
  SUBROUTINE CoordinateSystem_CreateStart(userNumber,coordinateSystems,coordinateSystem,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number for the created coordinate system
    TYPE(CoordinateSystemsType), POINTER :: coordinateSystems !<A pointer to the coordinate systems to create the coordinate system for.
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<On exit, a pointer to the created coordinate system. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coordinateSystemIdx
    TYPE(CoordinateSystemType), POINTER :: newCoordinateSystem
    TYPE(CoordinateSystemPtrType), POINTER :: newCoordinateSystems(:)
    TYPE(VARYING_STRING) :: localError

    NULLIFY(newCoordinateSystem)
    NULLIFY(newCoordinateSystems)

    ENTERS("CoordinateSystem_CreateStart",err,error,*998)

    IF(ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(coordinateSystems)) CALL FlagError("Coordinate systems is not associated.",err,error,*999)
   
    NULLIFY(newCoordinateSystem)
    CALL CoordinateSystem_UserNumberFind(coordinateSystems,userNumber,newCoordinateSystem,err,error,*999)
    IF(ASSOCIATED(newCoordinateSystem)) THEN
      localError="Coordinate system number "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " has already been created."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    
    NULLIFY(newCoordinateSystem)
    ALLOCATE(newCoordinateSystem,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new coordinate system.",err,error,*999)
      
    newCoordinateSystem%userNumber=userNumber
    newCoordinateSystem%coordinateSystems=>coordinateSystems
    newCoordinateSystem%coordinateSystemFinished=.FALSE.
    newCoordinateSystem%TYPE=COORDINATE_RECTANGULAR_CARTESIAN_TYPE
    newCoordinateSystem%radialInterpolationType=COORDINATE_NO_RADIAL_INTERPOLATION_TYPE
    newCoordinateSystem%numberOfDimensions=3
    newCoordinateSystem%focus=1.0_DP    
    newCoordinateSystem%origin=[0.0_DP,0.0_DP,0.0_DP]
    newCoordinateSystem%orientation=RESHAPE(&
      & [1.0_DP,0.0_DP,0.0_DP, &
      &   0.0_DP,1.0_DP,0.0_DP, &
      &   0.0_DP,0.0_DP,1.0_DP], &
      & [3,3])
    
    ALLOCATE(newCoordinateSystems(coordinateSystems%numberOfCoordinateSystems+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new coordinate systems.",err,error,*999)
    DO coordinateSystemIdx=1,coordinateSystems%numberOfCoordinateSystems
      newCoordinateSystems(coordinateSystemIdx)%ptr=>coordinateSystems%coordinateSystems(coordinateSystemIdx)%ptr
    ENDDO !coordinateSystemIdx
    newCoordinateSystems(coordinateSystems%numberOfCoordinateSystems+1)%ptr=>newCoordinateSystem
    DEALLOCATE(coordinateSystems%coordinateSystems)
    coordinateSystems%coordinateSystems=>newCoordinateSystems
    coordinateSystems%numberOfCoordinateSystems=coordinateSystems%numberOfCoordinateSystems+1
    
    coordinateSystem=>newCoordinateSystem
       
    EXITS("CoordinateSystem_CreateStart")
    RETURN
999 IF(ASSOCIATED(newCoordinateSystem)) DEALLOCATE(newCoordinateSystem)
    IF(ASSOCIATED(newCoordinateSystems)) DEALLOCATE(newCoordinateSystems)
    NULLIFY(coordinateSystem)
998 ERRORSEXITS("CoordinateSystem_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_CreateStart

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a new coordinate system. \see OpenCMISS::Iron::cmfe_CoordinateSystem_CreateFinish
  SUBROUTINE CoordinateSystem_CreateFinish(coordinateSystem,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<A pointer to the coordinate system to finish
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coordinateSystemIdx
    TYPE(CoordinateSystemsType), POINTER :: coordinateSystems

    ENTERS("CoordinateSystem_CreateFinish",err,error,*999)

    CALL CoordinateSystem_AssertNotFinished(coordinateSystem,err,error,*999)
      
    coordinateSystem%coordinateSystemFinished=.TRUE.
     
    IF(diagnostics1) THEN
      NULLIFY(coordinateSystems)
      CALL CoordinateSystem_CoordinateSystemsGet(coordinateSystem,coordinateSystems,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of coordinate systems = ", &
        & coordinateSystems%numberOfCoordinateSystems,err,error,*999)
      DO coordinateSystemIdx=1,coordinateSystems%numberOfCoordinateSystems
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Coordinate system : ",coordinateSystemIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number = ", &
          & coordinateSystems%coordinateSystems(coordinateSystemIdx)%ptr%userNumber,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Type = ", &
          & COORDINATE_SYSTEM_TYPE_STRING(coordinateSystems%coordinateSystems(coordinateSystemIdx)%ptr%type),err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of dimensions = ", &
          & coordinateSystems%coordinateSystems(coordinateSystemIdx)%ptr%numberOfDimensions,err,error,*999)
      ENDDO !coordinateSystemIdx
    ENDIF
    
    EXITS("CoordinateSystem_CreateFinish")
    RETURN
999 ERRORSEXITS("CoordinateSystem_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_CreateFinish

  !
  !================================================================================================================================
  !

  !>Destroys a coordinate system. \see OpenCMISS::Iron::cmfe_CoordinateSystem_Destroy
  SUBROUTINE CoordinateSystem_Destroy(coordinateSystem,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<A pointer to the coordinate system to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coordinateSystemNumber,newCoordinateSystemNumber
    LOGICAL :: found
    TYPE(CoordinateSystemPtrType), POINTER :: newCoordinateSystems(:)
    TYPE(CoordinateSystemsType), POINTER :: coordinateSystems
    TYPE(VARYING_STRING) :: localError

    ENTERS("CoordinateSystem_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is not associated.",err,error,*999)
    IF(coordinateSystem%userNumber==0) CALL FlagError("Cannot destroy the world coordinate system.",err,error,*999)
    
    NULLIFY(coordinateSystems)
    CALL CoordinateSystem_CoordinateSystemsGet(coordinateSystem,coordinateSystems,err,error,*999)
    found=.FALSE.
    newCoordinateSystemNumber=0
    ALLOCATE(newCoordinateSystems(coordinateSystems%numberOfCoordinateSystems-1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new coordianate systems.",err,error,*999)
    DO coordinateSystemNumber=1,coordinateSystems%numberOfCoordinateSystems
      IF(coordinateSystems%coordinateSystems(coordinateSystemNumber)%ptr%userNumber==coordinateSystem%userNumber) THEN
        found=.TRUE.
      ELSE
        newCoordinateSystemNumber=newCoordinateSystemNumber+1
        newCoordinateSystems(newCoordinateSystemNumber)%ptr=>coordinateSystems%coordinateSystems(coordinateSystemNumber)%ptr
      ENDIF
    ENDDO !coordinateSystemNumber
    IF(found) THEN
      CALL CoordinateSystem_Finalise(coordinateSystem,err,error,*999)
      DEALLOCATE(coordinateSystems%coordinateSystems)
      coordinateSystems%coordinateSystems=>newCoordinateSystems
      coordinateSystems%numberOfCoordinateSystems=coordinateSystems%numberOfCoordinateSystems-1
    ELSE
      DEALLOCATE(newCoordinateSystems)
      localError="The coordinate system number "//TRIM(NumberToVString(coordinateSystem%userNumber,"*",err,error))// &
        & " to destroy does not exist."
      CALL FlagError(localError,err,error,*999)
     ENDIF
       
    EXITS("CoordinateSystem_Destroy")
    RETURN
999 ERRORSEXITS("CoordinateSystem_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_Destroy

  !
  !================================================================================================================================
  !
  
  !>Calculates dx(:)/dz(i) at x, where z(i) are rectangular cartesian and x(:) are curvilinear coordinates defined by coordinateSystem for double precision coordinates.
  FUNCTION dXdZDP(coordinateSystem,coordinateIdx,x,err,error)
  
    !Argument variables
    TYPE(CoordinateSystemType), INTENT(IN) :: coordinateSystem
    INTEGER(INTG), INTENT(IN) :: coordinateIdx
    REAL(DP), INTENT(IN) :: x(:)
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Function variable
    REAL(DP) :: dXdZDP(SIZE(x,1))
    !Local variables
    REAL(DP) :: rd,focus
    TYPE(VARYING_STRING) :: localError

    ENTERS("dXdZDP",err,error,*999)

    dXdZDP=0.0_DP

    IF(SIZE(x,1)<coordinateSystem%numberOfDimensions) THEN
      localError="The size of x of "//TRIM(NumberToVString(SIZE(x,1),"*",err,error))// &
        & " is less than the coordinate system dimension of "// &
        & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(coordinateIdx<1.OR.coordinateIdx>coordinateSystem%numberOfDimensions) THEN
      localError="The specified coordinate index of "//TRIM(NumberToVString(coordinateIdx,"*",err,error))// &
        & " is invalid. The coordinate index must be >= 1 and <= "// &
        & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(coordinateSystem%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      dXdZDP(coordinateIdx)=1.0_DP
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      SELECT CASE(coordinateSystem%numberOfDimensions)
      CASE(2)
        SELECT CASE(coordinateIdx)
        CASE(1)
          dXdZDP(1)=COS(x(2))
          dXdZDP(2)=-SIN(x(2))/x(1)
        CASE(2)
          dXdZDP(1)=SIN(x(2))
          dXdZDP(2)=COS(x(2))/x(1)
        CASE DEFAULT
          CALL FlagError("Invalid coordinate index.",err,error,*999)
        END SELECT
      CASE(3)
        SELECT CASE(coordinateIdx)
        CASE(1)
          dXdZDP(1)=COS(x(2))
          dXdZDP(2)=-SIN(x(2))/x(1)
          dXdZDP(3)=0.0_DP
        CASE(2)
          dXdZDP(1)=SIN(x(2))
          dXdZDP(2)=COS(x(2))/x(1)
          dXdZDP(3)=0.0_DP
        CASE(3)
          dXdZDP(1)=0.0_DP
          dXdZDP(2)=0.0_DP
          dXdZDP(3)=1.0_DP
        CASE DEFAULT
          CALL FlagError("Invalid coordinate index.",err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The number of coordinate system dimensions of "// &
          & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(coordinateSystem%numberOfDimensions==3) THEN
        SELECT CASE(coordinateIdx)
        CASE(1)
          dXdZDP(1)=COS(x(2))*COS(x(3))
          dXdZDP(2)=-SIN(x(2))/(x(1)*COS(x(3)))
          dXdZDP(3)=-COS(x(2))*SIN(x(3))/x(1)
        CASE(2)
          dXdZDP(1)=SIN(x(2))*COS(x(3))
          dXdZDP(2)=COS(x(2))/(x(1)*COS(x(3)))
          dXdZDP(3)=-SIN(x(2))*SIN(x(3))/x(1)
        CASE(3)
          dXdZDP(1)=SIN(x(3))
          dXdZDP(2)=0.0_DP
          dXdZDP(3)=COS(x(3))/x(1)
        CASE DEFAULT
          CALL FlagError("Invalid coordinate index.",err,error,*999)
        END SELECT
      ELSE
        localError="The number of coordinate system dimensions of "// &
          & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      IF(coordinateSystem%numberOfDimensions==3) THEN
        focus=coordinateSystem%focus
        rd=focus*(COSH(x(1))*COSH(x(1))-COS(x(2))*COS(x(2)))
        SELECT CASE(coordinateIdx)
        CASE(1)
          dXdZDP(1)=SINH(x(1))*COS(x(2))/rd
          dXdZDP(2)=-COSH(x(1))*SIN(x(2))/rd
          dXdZDP(3)=0.0_DP
        CASE(2)
          dXdZDP(1)=COSH(x(1))*SIN(x(2))*COS(x(3))/rd
          dXdZDP(2)=SINH(x(1))*COS(x(2))*COS(x(3))/rd
          dXdZDP(3)=-SIN(x(3))/(focus*SINH(x(1))*SIN(x(2)))
        CASE(3)
          dXdZDP(1)=COSH(x(1))*SIN(x(2))*SIN(x(3))/rd
          dXdZDP(2)=SINH(x(1))*COS(x(2))*SIN(x(3))/rd
          dXdZDP(3)=COS(x(3))/(focus*SINH(x(1))*SIN(x(2)))
        CASE DEFAULT
          CALL FlagError("Invalid coordinate index.",err,error,*999)
        END SELECT
      ELSE
        localError="The number of coordinate system dimensions of "// &
          & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      CALL FlagError("Not implemented",err,error,*999)
    CASE DEFAULT
      localError="The specified coordinate system type of "// &
        & TRIM(NumberToVString(coordinateSystem%TYPE,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("dXdZDP")
    RETURN
999 ERRORSEXITS("dXdZDP",err,error)
    RETURN
    
  END FUNCTION dXdZDP

  !
  !================================================================================================================================
  !

  !>Calculates d2z(:)/dx(i)dx(j) at x(:), where z(:) are rectangalar Cartesian and x(i) and x(j) are curvilinear coordinates defined by coordinateSystem.
  FUNCTION d2ZdX1dX2DP(coordinateSystem,coordinateIdx1,coordinateIdx2,x,err,error)
  
   
    !Argument variables
    TYPE(CoordinateSystemType), INTENT(IN) :: coordinateSystem
    INTEGER(INTG), INTENT(IN) :: coordinateIdx1
    INTEGER(INTG), INTENT(IN) :: coordinateIdx2
    REAL(DP), INTENT(IN) :: x(:)
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Function variable
    REAL(DP) :: d2ZdX1dX2DP(SIZE(x,1))
    !Local variables
    REAL(DP) :: focus
    TYPE(VARYING_STRING) :: localError

    ENTERS("d2ZdX1dX2DP",err,error,*999)

    d2ZdX1dX2DP=0.0_DP

    IF(SIZE(x,1)<coordinateSystem%numberOfDimensions) THEN
      localError="The size of x of "//TRIM(NumberToVString(SIZE(x,1),"*",err,error))// &
        & " is less than the coordinate system dimension of "// &
        & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(coordinateIdx1<1.OR.coordinateIdx1>coordinateSystem%numberOfDimensions) THEN
      localError="The specified coordinate index 1 of "//TRIM(NumberToVString(coordinateIdx1,"*",err,error))// &
        & " is invalid. The coordinate index must be >= 1 and <= "// &
        & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(coordinateIdx2<1.OR.coordinateIdx2>coordinateSystem%numberOfDimensions) THEN
      localError="The specified coordinate index 2 of "//TRIM(NumberToVString(coordinateIdx2,"*",err,error))// &
        & " is invalid. The coordinate index must be >= 1 and <= "// &
        & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
   
    SELECT CASE(coordinateSystem%type)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      d2ZdX1dX2DP(1:coordinateSystem%numberOfDimensions)=0.0_DP
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      SELECT CASE(coordinateSystem%numberOfDimensions)
      CASE(2)
        SELECT CASE(coordinateIdx1)
        CASE(1)
          SELECT CASE(coordinateIdx2)
          CASE(1)
            d2ZdX1dX2DP(1)=0.0_DP
            d2ZdX1dX2DP(2)=0.0_DP
          CASE(2)
            d2ZdX1dX2DP(1)=-SIN(x(2))
            d2ZdX1dX2DP(2)=COS(x(2))
          CASE DEFAULT
            CALL FlagError("Invalid coordinate index 2.",err,error,*999)
          END SELECT
        CASE(2)
          SELECT CASE(coordinateIdx2)
          CASE(1)
            d2ZdX1dX2DP(1)=-SIN(x(2))
            d2ZdX1dX2DP(2)=COS(x(2))
          CASE(2)
            d2ZdX1dX2DP(1)=-x(1)*COS(x(2))
            d2ZdX1dX2DP(2)=-x(1)*SIN(x(2))
          CASE DEFAULT
            CALL FlagError("Invalid coordinate index 2.",err,error,*999)
          END SELECT
        CASE DEFAULT
          CALL FlagError("Invalid coordinate index 1.",err,error,*999)
        END SELECT
      CASE(3)
        SELECT CASE(coordinateIdx1)
        CASE(1)
          SELECT CASE(coordinateIdx2)
          CASE(1)
            d2ZdX1dX2DP(1)=0.0_DP
            d2ZdX1dX2DP(2)=0.0_DP
            d2ZdX1dX2DP(3)=0.0_DP
          CASE(2)
            d2ZdX1dX2DP(1)=-SIN(x(2))
            d2ZdX1dX2DP(2)=COS(x(2))
            d2ZdX1dX2DP(3)=0.0_DP
          CASE(3)
            d2ZdX1dX2DP(1)=0.0_DP
            d2ZdX1dX2DP(2)=0.0_DP
            d2ZdX1dX2DP(3)=0.0_DP
          CASE DEFAULT
            CALL FlagError("Invalid coordinate index 2.",err,error,*999)
          END SELECT
        CASE(2)
          SELECT CASE(coordinateIdx2)
          CASE(1)
            d2ZdX1dX2DP(1)=-SIN(x(2))
            d2ZdX1dX2DP(2)=COS(x(2))
            d2ZdX1dX2DP(3)=0.0_DP
          CASE(2)
            d2ZdX1dX2DP(1)=-x(1)*COS(x(2))
            d2ZdX1dX2DP(2)=-x(1)*SIN(x(2))
            d2ZdX1dX2DP(3)=0.0_DP
          CASE(3)
            d2ZdX1dX2DP(1)=0.0_DP
            d2ZdX1dX2DP(2)=0.0_DP
            d2ZdX1dX2DP(3)=0.0_DP
          CASE DEFAULT
            CALL FlagError("Invalid coordinate index 2.",err,error,*999)
          END SELECT
        CASE(3)
          SELECT CASE(coordinateIdx2)
          CASE(1)
            d2ZdX1dX2DP(1)=0.0_DP
            d2ZdX1dX2DP(2)=0.0_DP
            d2ZdX1dX2DP(3)=0.0_DP
          CASE(2)
            d2ZdX1dX2DP(1)=0.0_DP
            d2ZdX1dX2DP(2)=0.0_DP
            d2ZdX1dX2DP(3)=0.0_DP
          CASE(3)
            d2ZdX1dX2DP(1)=0.0_DP
            d2ZdX1dX2DP(2)=0.0_DP
            d2ZdX1dX2DP(3)=0.0_DP
          CASE DEFAULT
            CALL FlagError("Invalid coordinate index 2.",err,error,*999)
          END SELECT
        CASE DEFAULT
          CALL FlagError("Invalid coordinate index 1.",err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The number of coordinate system dimensions of "// &
          & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(coordinateSystem%numberOfDimensions==3) THEN
        SELECT CASE(coordinateIdx1)
        CASE(1)
          SELECT CASE(coordinateIdx2)
          CASE(1)
            d2ZdX1dX2DP(1)=0.0_DP
            d2ZdX1dX2DP(2)=0.0_DP
            d2ZdX1dX2DP(3)=0.0_DP
          CASE(2)
            d2ZdX1dX2DP(1)=-SIN(x(2))*COS(x(3))
            d2ZdX1dX2DP(2)=COS(x(2))*COS(x(3))
            d2ZdX1dX2DP(3)=0.0_DP
          CASE(3)
            d2ZdX1dX2DP(1)=-COS(x(2))*SIN(x(3))
            d2ZdX1dX2DP(2)=-SIN(x(2))*SIN(x(3))
            d2ZdX1dX2DP(3)=COS(x(3))
          CASE DEFAULT
            CALL FlagError("Invalid coordinate index 2.",err,error,*999)
          END SELECT
        CASE(2)
          SELECT CASE(coordinateIdx2)
          CASE(1)
            d2ZdX1dX2DP(1)=-SIN(x(2))*COS(x(3))
            d2ZdX1dX2DP(2)=COS(x(2))*COS(x(3))
            d2ZdX1dX2DP(3)=0.0_DP
          CASE(2)
            d2ZdX1dX2DP(1)=-x(1)*COS(x(2))*COS(x(3))
            d2ZdX1dX2DP(2)=-x(1)*SIN(x(2))*COS(x(3))
            d2ZdX1dX2DP(3)=0.0_DP
          CASE(3)
            d2ZdX1dX2DP(1)=x(1)*SIN(x(2))*SIN(x(3))
            d2ZdX1dX2DP(2)=-x(1)*COS(x(2))*SIN(x(3))
            d2ZdX1dX2DP(3)=0.0_DP
          CASE DEFAULT
            CALL FlagError("Invalid coordinate index 2.",err,error,*999)
          END SELECT
        CASE(3)
          SELECT CASE(coordinateIdx2)
          CASE(1)
            d2ZdX1dX2DP(1)=-COS(x(2))*SIN(x(3))
            d2ZdX1dX2DP(2)=-SIN(x(2))*SIN(x(3))
            d2ZdX1dX2DP(3)=COS(x(3))
          CASE(2)
            d2ZdX1dX2DP(1)=x(1)*SIN(x(2))*SIN(x(3))
            d2ZdX1dX2DP(2)=-x(1)*COS(x(2))*SIN(x(3))
            d2ZdX1dX2DP(3)=0.0_DP
          CASE(3)
            d2ZdX1dX2DP(1)=-x(1)*COS(x(2))*COS(x(3))
            d2ZdX1dX2DP(2)=-x(1)*SIN(x(2))*COS(x(3))
            d2ZdX1dX2DP(3)=-x(1)*SIN(x(3))
          CASE DEFAULT
            CALL FlagError("Invalid coordinate index 2.",err,error,*999)
          END SELECT
        CASE DEFAULT
          CALL FlagError("Invalid coordinate index 1.",err,error,*999)
        END SELECT
      ELSE
        localError="The number of coordinate system dimensions of "// &
          & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      IF(coordinateSystem%numberOfDimensions==3) THEN
        focus=coordinateSystem%focus
        SELECT CASE(coordinateIdx1)
        CASE(1)
          SELECT CASE(coordinateIdx2)
          CASE(1)
            d2ZdX1dX2DP(1)=focus*COSH(x(1))*COS(x(2))          
            d2ZdX1dX2DP(2)=focus*SINH(x(1))*SIN(x(2))*COS(x(3))
            d2ZdX1dX2DP(3)=focus*SINH(x(1))*SIN(x(2))*SIN(x(3))
          CASE(2)
            d2ZdX1dX2DP(1)=-focus*SINH(x(1))*SIN(x(2))
            d2ZdX1dX2DP(2)=focus*COSH(x(1))*COS(x(2))*COS(x(3))
            d2ZdX1dX2DP(3)=focus*COSH(x(1))*COS(x(2))*SIN(x(3))
          CASE(3)
            d2ZdX1dX2DP(1)=0.0_DP
            d2ZdX1dX2DP(2)=-focus*COSH(x(1))*SIN(x(2))*SIN(x(3))
            d2ZdX1dX2DP(3)=focus*COSH(x(1))*SIN(x(2))*COS(x(3))
          CASE DEFAULT
            CALL FlagError("Invalid coordinate index 2.",err,error,*999)
          END SELECT
        CASE(2)
          SELECT CASE(coordinateIdx2)
          CASE(1)
            d2ZdX1dX2DP(1)=-focus*SINH(x(1))*SIN(x(2))
            d2ZdX1dX2DP(2)=focus*COSH(x(1))*COS(x(2))*COS(x(3))
            d2ZdX1dX2DP(3)=focus*COSH(x(1))*COS(x(2))*SIN(x(3))
          CASE(2)
            d2ZdX1dX2DP(1)=-focus*COSH(x(1))*COS(x(2))
            d2ZdX1dX2DP(2)=-focus*SINH(x(1))*SIN(x(2))*COS(x(3))
            d2ZdX1dX2DP(3)=-focus*SINH(x(1))*SIN(x(2))*SIN(x(3))
          CASE(3)
            d2ZdX1dX2DP(1)=0.0_DP
            d2ZdX1dX2DP(2)=-focus*SINH(x(1))*COS(x(2))*SIN(x(3))
            d2ZdX1dX2DP(3)=focus*SINH(x(1))*COS(x(2))*COS(x(3))
          CASE DEFAULT
            CALL FlagError("Invalid coordinate index 2.",err,error,*999)
          END SELECT
        CASE(3)
          SELECT CASE(coordinateIdx2)
          CASE(1)
            d2ZdX1dX2DP(1)=0.0_DP
            d2ZdX1dX2DP(2)=-focus*COSH(x(1))*SIN(x(2))*SIN(x(3))
            d2ZdX1dX2DP(3)=focus*COSH(x(1))*SIN(x(2))*COS(x(3))
          CASE(2)
            d2ZdX1dX2DP(1)=0.0_DP
            d2ZdX1dX2DP(2)=-focus*SINH(x(1))*COS(x(2))*SIN(x(3))
            d2ZdX1dX2DP(3)=focus*SINH(x(1))*COS(x(2))*COS(x(3))
          CASE(3)
            d2ZdX1dX2DP(1)=0.0_DP
            d2ZdX1dX2DP(2)=-focus*SINH(x(1))*SIN(x(2))*COS(x(3))
            d2ZdX1dX2DP(3)=-focus*SINH(x(1))*SIN(x(2))*SIN(x(3))
          CASE DEFAULT
            CALL FlagError("Invalid coordinate index 2.",err,error,*999)
          END SELECT
        CASE DEFAULT
          CALL FlagError("Invalid coordinate index 1.",err,error,*999)
        END SELECT
      ELSE
        localError="The number of coordinate system dimensions of "// &
          & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      CALL FlagError("Not implemented",err,error,*999)
    CASE DEFAULT
      localError="The specified coordinate system type of "// &
        & TRIM(NumberToVString(coordinateSystem%type,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("d2ZdX1dX2DP")
    RETURN
999 ERRORSEXITS("d2ZdX1dX2DP",err,error)
    RETURN
    
  END FUNCTION d2ZdX1dX2DP

  !
  !================================================================================================================================
  !

  !>Calculates dZ(:)/dX(i) at x, where z(:) are rectangalar Cartesian and x(i) are curvilinear coordinates defined by coordinateSystem.
  FUNCTION dZdXDP(coordinateSystem,coordinateIdx,x,err,error)
     
    !Argument variables
    TYPE(CoordinateSystemType), INTENT(IN) :: coordinateSystem
    INTEGER(INTG), INTENT(IN) :: coordinateIdx
    REAL(DP), INTENT(IN) :: x(:)
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Function variable
    REAL(DP) :: dZdXDP(SIZE(x,1))
    !Local variables
    REAL(DP) :: focus
    TYPE(VARYING_STRING) :: localError

    ENTERS("dZdXDP",err,error,*999)

    dZdXDP=0.0_DP

   IF(SIZE(x,1)<coordinateSystem%numberOfDimensions) THEN
      localError="The size of x of "//TRIM(NumberToVString(SIZE(x,1),"*",err,error))// &
        & " is less than the coordinate system dimension of "// &
        & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(coordinateIdx<1.OR.coordinateIdx>coordinateSystem%numberOfDimensions) THEN
      localError="The specified coordinate index of "//TRIM(NumberToVString(coordinateIdx,"*",err,error))// &
        & " is invalid. The coordinate index must be >= 1 and <= "// &
        & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
   
   SELECT CASE(coordinateSystem%type)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      dZdXDP(coordinateIdx)=1.0_DP
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      SELECT CASE(coordinateSystem%numberOfDimensions)
      CASE(2)
        SELECT CASE(coordinateIdx)
        CASE(1)
          dZdXDP(1)=COS(x(2))
          dZdXDP(2)=SIN(x(2))
        CASE(2)
          dZdXDP(1)=-x(1)*SIN(x(2))
          dZdXDP(2)=x(1)*COS(x(2))
        CASE DEFAULT
          CALL FlagError("Invalid coordinate index.",err,error,*999)
        END SELECT
      CASE(3)
        SELECT CASE(coordinateIdx)
        CASE(1)
          dZdXDP(1)=COS(x(2))
          dZdXDP(2)=SIN(x(2))
          dZdXDP(3)=0.0_DP
        CASE(2)
          dZdXDP(1)=-x(1)*SIN(x(2))
          dZdXDP(2)=x(1)*COS(x(2))
          dZdXDP(3)=0.0_DP
        CASE(3)
          dZdXDP(1)=0.0_DP
          dZdXDP(2)=0.0_DP
          dZdXDP(3)=1.0_DP
        CASE DEFAULT
          CALL FlagError("Invalid coordinate index.",err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The number of coordinate system dimensions of "// &
          & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(coordinateSystem%numberOfDimensions==3) THEN
        SELECT CASE(coordinateIdx)
        CASE(1)
          dZdXDP(1)=COS(x(2))*COS(x(3))
          dZdXDP(2)=SIN(x(2))*COS(x(3))
          dZdXDP(3)=SIN(x(3))
        CASE(2)
          dZdXDP(1)=-x(1)*SIN(x(2))*COS(x(3))
          dZdXDP(2)=x(1)*COS(x(2))*COS(x(3))
          dZdXDP(3)=0.0_DP
        CASE(3)
          dZdXDP(1)=-x(1)*COS(x(2))*SIN(x(3))
          dZdXDP(2)=-x(1)*SIN(x(2))*SIN(x(3))
          dZdXDP(3)=x(1)*COS(x(3))
        CASE DEFAULT
          CALL FlagError("Invalid coordinate index.",err,error,*999)
        END SELECT
      ELSE
        localError="The number of coordinate system dimensions of "// &
          & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      IF(coordinateSystem%numberOfDimensions==3) THEN
        focus=coordinateSystem%focus
        SELECT CASE(coordinateIdx)
        CASE(1)
          dZdXDP(1)=focus*SINH(x(1))*COS(x(2))
          dZdXDP(2)=focus*COSH(x(1))*SIN(x(2))*COS(x(3))
          dZdXDP(3)=focus*COSH(x(1))*SIN(x(2))*SIN(x(3))
        CASE(2)
          dZdXDP(1)=-focus*COSH(x(1))*SIN(x(2))
          dZdXDP(2)=focus*SINH(x(1))*COS(x(2))*COS(x(3))
          dZdXDP(3)=focus*SINH(x(1))*COS(x(2))*SIN(x(3))
        CASE(3)
          dZdXDP(1)=0.0_DP
          dZdXDP(2)=-focus*SINH(x(1))*SIN(x(2))*SIN(x(3))
          dZdXDP(3)=focus*SINH(x(1))*SIN(x(2))*COS(x(3))
        CASE DEFAULT
          CALL FlagError("Invalid coordinate index.",err,error,*999)
        END SELECT
      ELSE
        localError="The number of coordinate system dimensions of "// &
          & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      IF(coordinateSystem%numberOfDimensions==3) THEN
        focus=coordinateSystem%focus
        SELECT CASE(coordinateIdx)
        CASE(1)
          dZdXDP(1)=focus*SINH(x(1))*COS(x(2))*COS(x(3))
          dZdXDP(2)=focus*COSH(x(1))*SIN(x(2))
          dZdXDP(3)=focus*SINH(x(1))*COS(x(2))*SIN(x(3))
        CASE(2)
          dZdXDP(1)=-focus*COSH(x(1))*SIN(x(2))*COS(x(3))
          dZdXDP(2)=focus*SINH(x(1))*COS(x(2))
          dZdXDP(3)=-focus*COSH(x(1))*SIN(x(2))*SIN(x(3))
        CASE(3)
          dZdXDP(1)=-focus*COSH(x(1))*COS(x(2))*SIN(x(3))
          dZdXDP(2)=0.0_DP
          dZdXDP(3)=focus*COSH(x(1))*COS(x(2))*COS(x(3))
        CASE DEFAULT
          CALL FlagError("Invalid coordinate index.",err,error,*999)
        END SELECT
      ELSE
        localError="The number of coordinate system dimensions of "// &
          & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The specified coordinate system type of "// &
        & TRIM(NumberToVString(coordinateSystem%type,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("dZdXDP")
    RETURN
999 ERRORSEXITS("dZdXDP",err,error)
    RETURN
    
  END FUNCTION dZdXDP

  !
  !================================================================================================================================
  !

  !>CoordinateSystem_DerivativeConvertToRC performs a coordinate transformation from a coordinate system identified by coordinateSystem at the point x with coordinates/derivatives x(coordinateIdx,partialDerivativeType) to the point with coordinates/derivatives z(coordinateIdx) in rectangular cartesian coordinates.
  SUBROUTINE CoordinateSystem_DerivativeConvertToRCDP(coordinateSystem,partialDerivativeType,x,z,err,error,*)
  
    !Argument variables
    TYPE(CoordinateSystemType), INTENT(IN) :: coordinateSystem
    INTEGER(INTG), INTENT(IN) :: partialDerivativeType
    REAL(DP), INTENT(IN) :: x(:,:)
    REAL(DP), INTENT(OUT) :: z(:)
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local variables
    REAL(DP) :: focus
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("CoordinateSystem_DerivativeConvertToRCDP",err,error,*999)
    
    IF(SIZE(x,1)<coordinateSystem%numberOfDimensions) THEN
      localError="The size of x of "//TRIM(NumberToVString(SIZE(x,1),"*",err,error))// &
        & " is less than the coordinate system dimension of "// &
        & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(x,1)/=SIZE(z,1)) THEN
      localError="The size of x of "//TRIM(NumberToVString(SIZE(x,1),"*",err,error))// &
        & " does not match the size of z of "//TRIM(NumberToVString(SIZE(z,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(x,2)<partialDerivativeType) THEN
      localError="The size of second index of x of "//TRIM(NumberToVString(SIZE(x,1),"*",err,error))// &
        & " is less than the partial derivative type of "//TRIM(NumberToVString(partialDerivativeType,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(partialDerivativeType<1.OR.partialDerivativeType>MAXIMUM_PARTIAL_DERIV_NUMBER) THEN
      localError="The specified partial derivative type of "//TRIM(NumberToVString(partialDerivativeType,"*",err,error))// &
        & " is invalid. The partial derivative type should be >= 1 and <= "// &
        & TRIM(NumberToVString(MAXIMUM_PARTIAL_DERIV_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    SELECT CASE(coordinateSystem%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      z=x(:,partialDerivativeType)
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      SELECT CASE(partialDerivativeType)
      CASE(NO_PART_DERIV)
        z=CoordinateSystem_ConvertToRC(coordinateSystem,x(:,NO_PART_DERIV),err,error)
        IF(err/=0) GOTO 999
      CASE(PART_DERIV_S1)
        SELECT CASE(coordinateSystem%numberOfDimensions)
        CASE(2)
          z(1)=COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)-x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1) !d(x)/d(s1)
          z(2)=SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1) !d(y)/d(s1)
        CASE(3)
          z(1)=COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)-x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1) !d(x)/d(s1)
          z(2)=SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1) !d(y)/d(s1)
          z(3)=x(3,PART_DERIV_S1) !d(z)/d(s1)
        CASE DEFAULT
          localError="The number of coordinate system dimensions of "// &
            & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PART_DERIV_S1_S1)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(PART_DERIV_S2)
        SELECT CASE(coordinateSystem%numberOfDimensions)
        CASE(2)
          z(1)=COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S2)-x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2) !d(x)/d(s2)
          z(2)=SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S2)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2) !d(y)/d(s2)
        CASE(3)
          z(1)=COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S2)-x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2) !d(x)/d(s2)
          z(2)=SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S2)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2) !d(y)/d(s2)
          z(3)=x(3,PART_DERIV_S2) !d(z)/d(s2)
        CASE DEFAULT
          localError="The number of coordinate system dimensions of "// &
            & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PART_DERIV_S2_S2)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(PART_DERIV_S1_S2)
        SELECT CASE(SIZE(x,1))
        CASE(2)
          z(1)=COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2)-x(1,PART_DERIV_S1)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2)-&
            & (x(1,PART_DERIV_S2)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*&
            & x(2,PART_DERIV_S1)*x(2,PART_DERIV_S2)+x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))* &
            & x(2,PART_DERIV_S1_S2)) !d2(x)/d(s1)d(s2)
          z(2)=SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2)+x(1,PART_DERIV_S1)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2)+&
            & (x(1,PART_DERIV_S2)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)-x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))*&
            & x(2,PART_DERIV_S1)*x(2,PART_DERIV_S2)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))* &
            & x(2,PART_DERIV_S1_S2)) !d2(y)/d(s1)d(s2)
        CASE(3)
          z(1)=COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2)-x(1,PART_DERIV_S1)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2)-&
            & (x(1,PART_DERIV_S2)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*&
            & x(2,PART_DERIV_S1)*x(2,PART_DERIV_S2)+x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))* &
            & x(2,PART_DERIV_S1_S2)) !d2(x)/d(s1)d(s2)
          z(2)=SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2)+x(1,PART_DERIV_S1)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2)+&
            & (x(1,PART_DERIV_S2)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)-x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))*&
            & x(2,PART_DERIV_S1)*x(2,PART_DERIV_S2)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))* &
            & x(2,PART_DERIV_S1_S2)) !d2(y)/d(s1)d(s2)
          z(3)=x(3,PART_DERIV_S1_S2) !d2(z)/d(s1)d(s2)
        CASE DEFAULT
          localError="The number of coordinate system dimensions of "// &
            & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PART_DERIV_S3)
        SELECT CASE(coordinateSystem%numberOfDimensions)
        CASE(2)
          z(1)=0.0_DP
          z(2)=0.0_DP
        CASE(3)
          z(1)=COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S3)-x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S3) !d(x)/d(s3)
          z(2)=SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S3)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S3) !d(y)/d(s3)
          z(3)=x(3,PART_DERIV_S3) !d(z)/d(s3)
        CASE DEFAULT
          localError="The number of coordinate system dimensions of "// &
            & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PART_DERIV_S3_S3)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(PART_DERIV_S1_S3)
        SELECT CASE(coordinateSystem%numberOfDimensions)
        CASE(2)
          z(1)=0.0_DP
          z(2)=0.0_DP
        CASE(3)
          z(1)=COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1_S3)-x(1,PART_DERIV_S1)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S3)-&
            & (x(1,PART_DERIV_S3)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*&
            & x(2,PART_DERIV_S1)*x(2,PART_DERIV_S3)+x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))* &
            & x(2,PART_DERIV_S1_S3)) !d2(x)/d(s1)d(s3)
          z(2)=SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1_S3)+x(1,PART_DERIV_S1)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S3)+&
            & (x(1,PART_DERIV_S3)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)-x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))*&
            & x(2,PART_DERIV_S1)*x(2,PART_DERIV_S3)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))* &
            & x(2,PART_DERIV_S1_S3)) !d2(y)/d(s1)d(s3)
          z(3)=x(3,PART_DERIV_S1_S3) !d2(z)/d(s1)d(s3)
        CASE DEFAULT
          localError="The number of coordinate system dimensions of "// &
            & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PART_DERIV_S2_S3)
        SELECT CASE(coordinateSystem%numberOfDimensions)
        CASE(2)
          z(1)=0.0_DP
          z(2)=0.0_DP
        CASE(3)
          z(1)=COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S2_S3)-x(1,PART_DERIV_S2)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S3)-&
            & (x(1,PART_DERIV_S3)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*&
            & x(2,PART_DERIV_S2)*x(2,PART_DERIV_S3)+x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))* &
            & x(2,PART_DERIV_S2_S3)) !d2(x)/d(s2)d(s3)
          z(2)=SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S2_S3)+x(1,PART_DERIV_S2)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S3)+&
            & (x(1,PART_DERIV_S3)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2)-x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))*&
            & x(2,PART_DERIV_S2)*x(2,PART_DERIV_S3)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))* &
            & x(2,PART_DERIV_S2_S3)) !d2(y)/d(s2)d(s3)
          z(3)=x(3,PART_DERIV_S2_S3) !d2(z)/d(s2)d(s3)
        CASE DEFAULT
          localError="The number of coordinate system dimensions of "// &
            & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PART_DERIV_S1_S2_S3)
        SELECT CASE(coordinateSystem%numberOfDimensions)
        CASE(2)
          z(1)=0.0_DP
          z(2)=0.0_DP
        CASE(3)  
          z(1)=-COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)*x(2,PART_DERIV_S2)*x(1,PART_DERIV_S3)-&
            & SIN(x(2,NO_PART_DERIV))*(x(2,PART_DERIV_S1_S2)*x(1,PART_DERIV_S3)+x(2,PART_DERIV_S2)*x(1,PART_DERIV_S1_S3))-&
            & SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)*x(1,PART_DERIV_S2_S3)+COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2_S3)-&
            & COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)*x(1,PART_DERIV_S2)*x(2,PART_DERIV_S3)-& 
            & SIN(x(2,NO_PART_DERIV))*(x(1,PART_DERIV_S1_S2)*x(2,PART_DERIV_S3)+x(1,PART_DERIV_S2)*x(2,PART_DERIV_S1_S3))+&
            & (-COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)+SIN(x(2,NO_PART_DERIV))*x(1,NO_PART_DERIV)*x(2,PART_DERIV_S1))*&
            & x(2,PART_DERIV_S2)*x(2,PART_DERIV_S3)-x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*(x(2,PART_DERIV_S1_S2)* &
            & x(2,PART_DERIV_S3)+x(2,PART_DERIV_S2)*x(2,PART_DERIV_S1_S3))+(-SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)- &
            & x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1))*x(2,PART_DERIV_S2_S3)-x(1,NO_PART_DERIV)*&
            & SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1_S2_S3) !d3(x)/d(s1)d(s2)d(s3)
          z(2)=-SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)*x(2,PART_DERIV_S2)*x(1,PART_DERIV_S3)+&
            & COS(x(2,NO_PART_DERIV))*(x(2,PART_DERIV_S1_S2)*x(1,PART_DERIV_S3)+x(2,PART_DERIV_S2)*x(1,PART_DERIV_S1_S3))+&
            & COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)*x(1,PART_DERIV_S2_S3)+SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2_S3)-&
            & SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)*x(1,PART_DERIV_S2)*x(2,PART_DERIV_S3)+&
            & COS(x(2,NO_PART_DERIV))*(x(1,PART_DERIV_S1_S2)*x(2,PART_DERIV_S3)+x(1,PART_DERIV_S2)*x(2,PART_DERIV_S1_S3))+&
            & (-SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)-x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1))*&
            & x(2,PART_DERIV_S2)*x(2,PART_DERIV_S3)-x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))*(x(2,PART_DERIV_S1_S2)* &
            & x(2,PART_DERIV_S3)+x(2,PART_DERIV_S2)*x(2,PART_DERIV_S1_S3))+(COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)- &
            & x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1))*x(2,PART_DERIV_S2_S3)+x(1,NO_PART_DERIV)*&
            & COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1_S2_S3) !d3(y)/d(s1)d(s2)d(s3)
          z(3)=x(3,PART_DERIV_S1_S2_S3) !d3(z)/d(s1)d(s2)d(s3)
        CASE DEFAULT
          localError="The number of coordinate system dimensions of "// &
            & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The specified partial derivative type of "// &
          & TRIM(NumberToVString(partialDerivativeType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(coordinateSystem%numberOfDimensions/=3) THEN
        localError="The number of coordinate system dimensions of "// &
          & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      SELECT CASE(partialDerivativeType)
      CASE(NO_PART_DERIV)
        z=CoordinateSystem_ConvertToRC(coordinateSystem,x(:,NO_PART_DERIV),err,error)
        IF(err/=0) GOTO 999
      CASE(PART_DERIV_S1,PART_DERIV_S1_S1,PART_DERIV_S2,PART_DERIV_S2_S2,&
        & PART_DERIV_S1_S2,PART_DERIV_S3,PART_DERIV_S3_S3,&
        & PART_DERIV_S1_S3,PART_DERIV_S2_S3,PART_DERIV_S1_S2_S3)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The specified partial derivative type of "// &
          & TRIM(NumberToVString(partialDerivativeType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      IF(coordinateSystem%numberOfDimensions/=3) THEN
        localError="The number of coordinate system dimensions of "// &
          & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      focus=coordinateSystem%focus
      SELECT CASE(partialDerivativeType)
      CASE(NO_PART_DERIV)
        z=CoordinateSystem_ConvertToRC(coordinateSystem,x(:,NO_PART_DERIV),err,error)
        IF(err/=0) GOTO 999
      CASE(PART_DERIV_S1)
        z(1)=focus*SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & focus*COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1) !d(x)/d(s1)
        z(2)=focus*COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & focus*SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & focus*SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S1) !d(y)/d(s1)
        z(3)=focus*COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & focus*SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)+&
          & focus*SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S1) !d(z)/d(s1)
      CASE(PART_DERIV_S1_S1)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(PART_DERIV_S2)
        z(1)=focus*SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S2)-&
          & focus*COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2) !d(x)/d(s2)
        z(2)=focus*COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S2)+&
          & focus*SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S2)-&
          & focus*SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S2) !d(y)/d(s2)
        z(3)=focus*COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S2)+&
          & focus*SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S2)+&
          & focus*SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S2) !d(z)/d(s2)
      CASE(PART_DERIV_S2_S2)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(PART_DERIV_S1_S2)
        z(1)=focus*(SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2)+&
          & x(1,PART_DERIV_S1)*(COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S2)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2))-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1_S2)-&
          & x(2,PART_DERIV_S1)*(SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S2)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2))) !d2(x)/d(s1)d(s2)
        z(2)=focus*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2)+&
          & x(1,PART_DERIV_S1)*(SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(1,PART_DERIV_S2)+COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(2,PART_DERIV_S2)-COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S2))+SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(2,PART_DERIV_S1_S2)+x(2,PART_DERIV_S1)*(COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(1,PART_DERIV_S2)-SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S2)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S2))-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1_S2)-&
          & x(3,PART_DERIV_S1)*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(1,PART_DERIV_S2)+SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(2,PART_DERIV_S2)+SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S2))) !d2(y)/d(s1)d(s2)
        z(3)=focus*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2)+&
          & x(1,PART_DERIV_S1)*(SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(1,PART_DERIV_S2)+COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(2,PART_DERIV_S2)+COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S2))+SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1_S2)+&
          & x(2,PART_DERIV_S1)*(COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(1,PART_DERIV_S2)-SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(2,PART_DERIV_S2)+SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S2))+SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1_S2)+&
          & x(3,PART_DERIV_S1)*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(1,PART_DERIV_S2)+SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(2,PART_DERIV_S2)-SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S2))) !d2(z)/d(s1)d(s2)
      CASE(PART_DERIV_S3)
        z(1)=focus*SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S3)-&
          & focus*COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S3) !d(x)/d(s3)
        z(2)=focus*COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & focus*SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)-&
          & focus*SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S3) !d(y)/d(s3)
        z(3)=focus*COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & focus*SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)+&
          & focus*SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S3) !d(z)/d(s3)
      CASE(PART_DERIV_S3_S3)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(PART_DERIV_S1_S3)
        z(1)=focus*(SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1_S3)+&
          & x(1,PART_DERIV_S1)*(COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S3))-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1_S3)-&
          & x(2,PART_DERIV_S1)*(SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S3))) !d2(x)/d(s1)d(s3)
        z(2)=focus*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1_S3)+&
          & x(1,PART_DERIV_S1)*(SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S3))+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1_S3)+&
          & x(2,PART_DERIV_S1)*(COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S3))-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1_S3)-&
          & x(3,PART_DERIV_S1)*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))* x(2,PART_DERIV_S3)+&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S3))) !d2(y)/d(s1)d(s3) 
        z(3)=focus*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1_S3)+&
          & x(1,PART_DERIV_S1)*(SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)+&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S3))+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1_S3)+&
          & x(2,PART_DERIV_S1)*(COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S3))+&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1_S3)+&
          & x(3,PART_DERIV_S1)*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S3))) !d2(z)/d(s1)d(s3)
      CASE(PART_DERIV_S2_S3)
        z(1)=focus*(SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S2_S3)+&
          & x(1,PART_DERIV_S2)*(COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S3))-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2_S3)-&
          & x(2,PART_DERIV_S2)*(SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S3))) !d2(x)/d(s2)d(s3)
        z(2)=focus*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S2_S3)+&
          & x(1,PART_DERIV_S2)*(SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S3))+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S2_S3)+&
          & x(2,PART_DERIV_S2)*(COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S3))-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S2_S3)-&
          & x(3,PART_DERIV_S2)*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)+&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S3))) !d2(y)/d(s2)d(s3)
        z(3)=focus*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S2_S3)+&
          & x(1,PART_DERIV_S2)*(SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)+&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S3))+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S2_S3)+&
          & x(2,PART_DERIV_S2)*(COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S3))+&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S2_S3)+&
          & x(3,PART_DERIV_S2)*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S3))) !d2(z)/d(s2)d(s3)
      CASE(PART_DERIV_S1_S2_S3)
        z(1)=focus*((SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1))*x(1,PART_DERIV_S2)*x(1,PART_DERIV_S3)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*(x(1,PART_DERIV_S1_S2)*x(1,PART_DERIV_S3)+x(1,PART_DERIV_S2)* &
          & x(1,PART_DERIV_S1_S3))+(-COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1))*x(2,PART_DERIV_S2)*x(1,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*(x(2,PART_DERIV_S1_S2)*x(1,PART_DERIV_S3)+x(2,PART_DERIV_S2)* &
          & x(1,PART_DERIV_S1_S3))+(COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1))*x(1,PART_DERIV_S2_S3)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2_S3)+&
          & (-COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1))*x(1,PART_DERIV_S2)*x(2,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*(x(1,PART_DERIV_S1_S2)*x(2,PART_DERIV_S3)+x(1,PART_DERIV_S2)* &
          & x(2,PART_DERIV_S1_S3))+(-SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1))*x(2,PART_DERIV_S2)*x(2,PART_DERIV_S3)-&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*(x(2,PART_DERIV_S1_S2)*x(2,PART_DERIV_S3)+x(2,PART_DERIV_S2)* &
          & x(2,PART_DERIV_S1_S3))+(-SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1))*x(2,PART_DERIV_S2_S3)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1_S2_S3)) !d3(x)/d(s1)d(s2)d(s3)
        z(2)=focus*((COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(1,PART_DERIV_S2)*x(1,PART_DERIV_S3)+SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & (x(1,PART_DERIV_S1_S2)*x(1,PART_DERIV_S3)+x(1,PART_DERIV_S2)*x(1,PART_DERIV_S1_S3))+&
          & (SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(2,PART_DERIV_S2)*x(1,PART_DERIV_S3)+COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & (x(2,PART_DERIV_S1_S2)*x(1,PART_DERIV_S3)+x(2,PART_DERIV_S2)*x(1,PART_DERIV_S1_S3))+&
          & (-SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(3,PART_DERIV_S2)*x(1,PART_DERIV_S3)-COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & (x(3,PART_DERIV_S1_S2)*x(1,PART_DERIV_S3)+x(3,PART_DERIV_S2)*x(1,PART_DERIV_S1_S3))+&
          & (SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*x(1,PART_DERIV_S2_S3)+&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2_S3))+&
          & focus*((SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(1,PART_DERIV_S2)*x(2,PART_DERIV_S3)+COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & (x(1,PART_DERIV_S1_S2)*x(2,PART_DERIV_S3)+x(1,PART_DERIV_S2)*x(2,PART_DERIV_S1_S3))+&
          & (-COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)+&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(2,PART_DERIV_S2)*x(2,PART_DERIV_S3)-SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & (x(2,PART_DERIV_S1_S2)*x(2,PART_DERIV_S3)+x(2,PART_DERIV_S2)*x(2,PART_DERIV_S1_S3))+&
          & (-COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(3,PART_DERIV_S2)*x(2,PART_DERIV_S3)-SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & (x(3,PART_DERIV_S1_S2)*x(2,PART_DERIV_S3)+x(3,PART_DERIV_S2)*x(2,PART_DERIV_S1_S3))+&
          & (COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*x(2,PART_DERIV_S2_S3)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1_S2_S3))+&
          & focus*((-SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(1,PART_DERIV_S2)*x(3,PART_DERIV_S3)-COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & (x(1,PART_DERIV_S1_S2)*x(3,PART_DERIV_S3)+x(1,PART_DERIV_S2)*x(3,PART_DERIV_S1_S3))+&
          & (-COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(2,PART_DERIV_S2)*x(3,PART_DERIV_S3)-SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & (x(2,PART_DERIV_S1_S2)*x(3,PART_DERIV_S3)+x(2,PART_DERIV_S2)*x(3,PART_DERIV_S1_S3))+&
          & (-COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)+&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(3,PART_DERIV_S2)*x(3,PART_DERIV_S3)-SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & (x(3,PART_DERIV_S1_S2)*x(3,PART_DERIV_S3)+x(3,PART_DERIV_S2)*x(3,PART_DERIV_S1_S3))+&
          & (-COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*x(3,PART_DERIV_S2_S3)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S1_S2_S3)) !d3(y)/d(s1)d(s2)d(s3)
        z(3)=focus*((COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)+&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(1,PART_DERIV_S2)*x(1,PART_DERIV_S3)+SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & (x(1,PART_DERIV_S1_S2)*x(1,PART_DERIV_S3)+x(1,PART_DERIV_S2)*x(1,PART_DERIV_S1_S3))+&
          & (SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(2,PART_DERIV_S2)*x(1,PART_DERIV_S3)+COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & (x(2,PART_DERIV_S1_S2)*x(1,PART_DERIV_S3)+x(2,PART_DERIV_S2)*x(1,PART_DERIV_S1_S3))+&
          & (SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(3,PART_DERIV_S2)*x(1,PART_DERIV_S3)+COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & (x(3,PART_DERIV_S1_S2)*x(1,PART_DERIV_S3)+x(3,PART_DERIV_S2)*x(1,PART_DERIV_S1_S3))+&
          & (SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)+&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*x(1,PART_DERIV_S2_S3)+&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2_S3))+&
          & focus*((SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(1,PART_DERIV_S2)*x(2,PART_DERIV_S3)+COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & (x(1,PART_DERIV_S1_S2)*x(2,PART_DERIV_S3)+x(1,PART_DERIV_S2)*x(2,PART_DERIV_S1_S3))+&
          & (-COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(2,PART_DERIV_S2)*x(2,PART_DERIV_S3)-SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & (x(2,PART_DERIV_S1_S2)*x(2,PART_DERIV_S3)+x(2,PART_DERIV_S2)*x(2,PART_DERIV_S1_S3))+&
          & (COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(3,PART_DERIV_S2)*x(2,PART_DERIV_S3)+SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & (x(3,PART_DERIV_S1_S2)*x(2,PART_DERIV_S3)+x(3,PART_DERIV_S2)*x(2,PART_DERIV_S1_S3))+&
          & (COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*x(2,PART_DERIV_S2_S3)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1_S2_S3))+&
          & focus*((SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(1,PART_DERIV_S2)*x(3,PART_DERIV_S3)+COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & (x(1,PART_DERIV_S1_S2)*x(3,PART_DERIV_S3)+x(1,PART_DERIV_S2)*x(3,PART_DERIV_S1_S3))+&
          & (COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(2,PART_DERIV_S2)*x(3,PART_DERIV_S3)+SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & (x(2,PART_DERIV_S1_S2)*x(3,PART_DERIV_S3)+x(2,PART_DERIV_S2)*x(3,PART_DERIV_S1_S3))+&
          & (-COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(3,PART_DERIV_S2)*x(3,PART_DERIV_S3)-SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & (x(3,PART_DERIV_S1_S2)*x(3,PART_DERIV_S3)+x(3,PART_DERIV_S2)*x(3,PART_DERIV_S1_S3))+&
          & (COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*x(3,PART_DERIV_S2_S3)+&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S1_S2_S3)) !d3(z)/d(s1)d(s2)d(s3)
      CASE DEFAULT
        localError="The specified partial derivative type of "// &
          & TRIM(NumberToVString(partialDerivativeType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      IF(coordinateSystem%numberOfDimensions/=3) THEN
        localError="The number of coordinate system dimensions of "// &
          & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      focus=coordinateSystem%focus
      SELECT CASE(partialDerivativeType)
      CASE(NO_PART_DERIV)
        z=CoordinateSystem_ConvertToRC(coordinateSystem,x(:,NO_PART_DERIV),err,error)
        IF(err/=0) GOTO 999
      CASE(PART_DERIV_S1,PART_DERIV_S1_S1,PART_DERIV_S2,PART_DERIV_S2_S2,&
        & PART_DERIV_S1_S2,PART_DERIV_S3,PART_DERIV_S3_S3,&
        & PART_DERIV_S1_S3,PART_DERIV_S2_S3,PART_DERIV_S1_S2_S3)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The specified partial derivative type of "// &
          & TRIM(NumberToVString(partialDerivativeType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The specified coordinate system type of "// &
        & TRIM(NumberToVString(coordinateSystem%TYPE,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
 
    EXITS("CoordinateSystem_DerivativeConvertToRCDP")
    RETURN
999 ERRORSEXITS("CoordinateSystem_DerivativeConvertToRCDP",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_DerivativeConvertToRCDP

  !
  !================================================================================================================================
  !

  !>CoordinateSystem_DerivativeConvertToRCSP performs a coordinate transformation from a coordinate system identified by coordinateSystem at the point x with coordinates/derivatives x(coordinateIdx,partialDerivativeIdx) to the point with coordinates/derivatives z(coordinateIdx) in rectangular cartesian coordinates for single precision coordinates.
  SUBROUTINE CoordinateSystem_DerivativeConvertToRCSP(coordinateSystem,partialDerivativeType,x,z,err,error,*)
  
   !Argument variables
    TYPE(CoordinateSystemType), INTENT(IN) :: coordinateSystem
    INTEGER(INTG), INTENT(IN) :: partialDerivativeType
    REAL(SP), INTENT(IN) :: x(:,:)
    REAL(SP), INTENT(OUT) :: z(:)
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local variables
    REAL(SP) :: focus
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("CoordinateSystem_DerivativeConvertToRCSP",err,error,*999)

    IF(SIZE(x,1)<coordinateSystem%numberOfDimensions) THEN
      localError="The size of x of "//TRIM(NumberToVString(SIZE(x,1),"*",err,error))// &
        & " is less than the coordinate system dimension of "// &
        & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(x,1)/=SIZE(z,1)) THEN
      localError="The size of x of "//TRIM(NumberToVString(SIZE(x,1),"*",err,error))// &
        & " does not match the size of z of "//TRIM(NumberToVString(SIZE(z,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(x,2)<partialDerivativeType) THEN
      localError="The size of second index of x of "//TRIM(NumberToVString(SIZE(x,1),"*",err,error))// &
        & " is less than the partial derivative type of "//TRIM(NumberToVString(partialDerivativeType,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(partialDerivativeType<1.OR.partialDerivativeType>MAXIMUM_PARTIAL_DERIV_NUMBER) THEN
      localError="The specified partial derivative type of "//TRIM(NumberToVString(partialDerivativeType,"*",err,error))// &
        & " is invalid. The partial derivative type should be >= 1 and <= "// &
        & TRIM(NumberToVString(MAXIMUM_PARTIAL_DERIV_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    SELECT CASE(coordinateSystem%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      z=x(:,partialDerivativeType)
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      SELECT CASE(partialDerivativeType)
      CASE(NO_PART_DERIV)
        z=CoordinateSystem_ConvertToRC(coordinateSystem,x(:,NO_PART_DERIV),err,error)
        IF(err/=0) GOTO 999
      CASE(PART_DERIV_S1)
        SELECT CASE(coordinateSystem%numberOfDimensions)
        CASE(2)
          z(1)=COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)-x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1) !d(x)/d(s1)
          z(2)=SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1) !d(y)/d(s1)
        CASE(3)
          z(1)=COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)-x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1) !d(x)/d(s1)
          z(2)=SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1) !d(y)/d(s1)
          z(3)=x(3,PART_DERIV_S1) !d(z)/d(s1)
        CASE DEFAULT
          localError="The number of coordinate system dimensions of "// &
            & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PART_DERIV_S1_S1)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(PART_DERIV_S2)
        SELECT CASE(coordinateSystem%numberOfDimensions)
        CASE(2)
          z(1)=COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S2)-x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2) !d(x)/d(s2)
          z(2)=SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S2)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2) !d(y)/d(s2)
        CASE(3)
          z(1)=COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S2)-x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2) !d(x)/d(s2)
          z(2)=SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S2)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2) !d(y)/d(s2)
          z(3)=x(3,PART_DERIV_S2) !d(z)/d(s2)
        CASE DEFAULT
          localError="The number of coordinate system dimensions of "// &
            & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PART_DERIV_S2_S2)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(PART_DERIV_S1_S2)
        SELECT CASE(SIZE(x,1))
        CASE(2)
          z(1)=COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2)-x(1,PART_DERIV_S1)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2)-&
            & (x(1,PART_DERIV_S2)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*&
            & x(2,PART_DERIV_S1)*x(2,PART_DERIV_S2)+x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))* &
            & x(2,PART_DERIV_S1_S2)) !d2(x)/d(s1)d(s2)
          z(2)=SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2)+x(1,PART_DERIV_S1)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2)+&
            & (x(1,PART_DERIV_S2)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)-x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))*&
            & x(2,PART_DERIV_S1)*x(2,PART_DERIV_S2)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))* &
            & x(2,PART_DERIV_S1_S2)) !d2(y)/d(s1)d(s2)
        CASE(3)
          z(1)=COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2)-x(1,PART_DERIV_S1)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2)-&
            & (x(1,PART_DERIV_S2)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*&
            & x(2,PART_DERIV_S1)*x(2,PART_DERIV_S2)+x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))* &
            & x(2,PART_DERIV_S1_S2)) !d2(x)/d(s1)d(s2)
          z(2)=SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2)+x(1,PART_DERIV_S1)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2)+&
            & (x(1,PART_DERIV_S2)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)-x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))*&
            & x(2,PART_DERIV_S1)*x(2,PART_DERIV_S2)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))* &
            & x(2,PART_DERIV_S1_S2)) !d2(y)/d(s1)d(s2)
          z(3)=x(3,PART_DERIV_S1_S2) !d2(z)/d(s1)d(s2)
        CASE DEFAULT
          localError="The number of coordinate system dimensions of "// &
            & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PART_DERIV_S3)
        SELECT CASE(coordinateSystem%numberOfDimensions)
        CASE(2)
          z(1)=0.0_SP
          z(2)=0.0_SP
        CASE(3)
          z(1)=COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S3)-x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S3) !d(x)/d(s3)
          z(2)=SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S3)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S3) !d(y)/d(s3)
          z(3)=x(3,PART_DERIV_S3) !d(z)/d(s3)
        CASE DEFAULT
          localError="The number of coordinate system dimensions of "// &
            & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PART_DERIV_S3_S3)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(PART_DERIV_S1_S3)
        SELECT CASE(coordinateSystem%numberOfDimensions)
        CASE(2)
          z(1)=0.0_SP
          z(2)=0.0_SP
        CASE(3)
          z(1)=COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1_S3)-x(1,PART_DERIV_S1)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S3)-&
            & (x(1,PART_DERIV_S3)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*&
            & x(2,PART_DERIV_S1)*x(2,PART_DERIV_S3)+x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))* &
            & x(2,PART_DERIV_S1_S3)) !d2(x)/d(s1)d(s3)
          z(2)=SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1_S3)+x(1,PART_DERIV_S1)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S3)+&
            & (x(1,PART_DERIV_S3)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)-x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))*&
            & x(2,PART_DERIV_S1)*x(2,PART_DERIV_S3)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))* &
            & x(2,PART_DERIV_S1_S3)) !d2(y)/d(s1)d(s3)
          z(3)=x(3,PART_DERIV_S1_S3) !d2(z)/d(s1)d(s3)
        CASE DEFAULT
          localError="The number of coordinate system dimensions of "// &
            & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PART_DERIV_S2_S3)
        SELECT CASE(coordinateSystem%numberOfDimensions)
        CASE(2)
          z(1)=0.0_SP
          z(2)=0.0_SP
        CASE(3)
          z(1)=COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S2_S3)-x(1,PART_DERIV_S2)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S3)-&
            & (x(1,PART_DERIV_S3)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*&
            & x(2,PART_DERIV_S2)*x(2,PART_DERIV_S3)+x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))* &
            & x(2,PART_DERIV_S2_S3)) !d2(x)/d(s2)d(s3)
          z(2)=SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S2_S3)+x(1,PART_DERIV_S2)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S3)+&
            & (x(1,PART_DERIV_S3)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2)-x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))*&
            & x(2,PART_DERIV_S2)*x(2,PART_DERIV_S3)+x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))* &
            & x(2,PART_DERIV_S2_S3)) !d2(y)/d(s2)d(s3)
          z(3)=x(3,PART_DERIV_S2_S3) !d2(z)/d(s2)d(s3)
        CASE DEFAULT
          localError="The number of coordinate system dimensions of "// &
            & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PART_DERIV_S1_S2_S3)
        SELECT CASE(coordinateSystem%numberOfDimensions)
        CASE(2)
          z(1)=0.0_SP
          z(2)=0.0_SP
        CASE(3)  
          z(1)=-COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)*x(2,PART_DERIV_S2)*x(1,PART_DERIV_S3)-&
            & SIN(x(2,NO_PART_DERIV))*(x(2,PART_DERIV_S1_S2)*x(1,PART_DERIV_S3)+x(2,PART_DERIV_S2)*x(1,PART_DERIV_S1_S3))-&
            & SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)*x(1,PART_DERIV_S2_S3)+COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2_S3)-&
            & COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)*x(1,PART_DERIV_S2)*x(2,PART_DERIV_S3)-& 
            & SIN(x(2,NO_PART_DERIV))*(x(1,PART_DERIV_S1_S2)*x(2,PART_DERIV_S3)+x(1,PART_DERIV_S2)*x(2,PART_DERIV_S1_S3))+&
            & (-COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)+SIN(x(2,NO_PART_DERIV))*x(1,NO_PART_DERIV)*x(2,PART_DERIV_S1))*&
            & x(2,PART_DERIV_S2)*x(2,PART_DERIV_S3)-x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*(x(2,PART_DERIV_S1_S2)* &
            & x(2,PART_DERIV_S3)+x(2,PART_DERIV_S2)*x(2,PART_DERIV_S1_S3))+(-SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)- &
            & x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1))*x(2,PART_DERIV_S2_S3)-x(1,NO_PART_DERIV)*&
            & SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1_S2_S3) !d3(x)/d(s1)d(s2)d(s3)
          z(2)=-SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)*x(2,PART_DERIV_S2)*x(1,PART_DERIV_S3)+&
            & COS(x(2,NO_PART_DERIV))*(x(2,PART_DERIV_S1_S2)*x(1,PART_DERIV_S3)+x(2,PART_DERIV_S2)*x(1,PART_DERIV_S1_S3))+&
            & COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)*x(1,PART_DERIV_S2_S3)+SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2_S3)-&
            & SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1)*x(1,PART_DERIV_S2)*x(2,PART_DERIV_S3)+&
            & COS(x(2,NO_PART_DERIV))*(x(1,PART_DERIV_S1_S2)*x(2,PART_DERIV_S3)+x(1,PART_DERIV_S2)*x(2,PART_DERIV_S1_S3))+&
            & (-SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)-x(1,NO_PART_DERIV)*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1))*&
            & x(2,PART_DERIV_S2)*x(2,PART_DERIV_S3)-x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))*(x(2,PART_DERIV_S1_S2)* &
            & x(2,PART_DERIV_S3)+x(2,PART_DERIV_S2)*x(2,PART_DERIV_S1_S3))+(COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)- &
            & x(1,NO_PART_DERIV)*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1))*x(2,PART_DERIV_S2_S3)+x(1,NO_PART_DERIV)*&
            & COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1_S2_S3) !d3(y)/d(s1)d(s2)d(s3)
          z(3)=x(3,PART_DERIV_S1_S2_S3) !d3(z)/d(s1)d(s2)d(s3)
        CASE DEFAULT
          localError="The number of coordinate system dimensions of "// &
            & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The specified partial derivative type of "// &
          & TRIM(NumberToVString(partialDerivativeType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      IF(coordinateSystem%numberOfDimensions/=3) THEN
        localError="The number of coordinate system dimensions of "// &
          & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      SELECT CASE(partialDerivativeType)
      CASE(NO_PART_DERIV)
        z=CoordinateSystem_ConvertToRC(coordinateSystem,x(:,NO_PART_DERIV),err,error)
        IF(err/=0) GOTO 999
      CASE(PART_DERIV_S1,PART_DERIV_S1_S1,PART_DERIV_S2,PART_DERIV_S2_S2,&
        & PART_DERIV_S1_S2,PART_DERIV_S3,PART_DERIV_S3_S3,&
        & PART_DERIV_S1_S3,PART_DERIV_S2_S3,PART_DERIV_S1_S2_S3)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The specified partial derivative type of "// &
          & TRIM(NumberToVString(partialDerivativeType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      IF(coordinateSystem%numberOfDimensions/=3) THEN
        localError="The number of coordinate system dimensions of "// &
          & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      focus=REAL(coordinateSystem%focus,SP)
      SELECT CASE(partialDerivativeType)
      CASE(NO_PART_DERIV)
        z=CoordinateSystem_ConvertToRC(coordinateSystem,x(:,NO_PART_DERIV),err,error)
        IF(err/=0) GOTO 999
      CASE(PART_DERIV_S1)
        z(1)=focus*SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & focus*COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1) !d(x)/d(s1)
        z(2)=focus*COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & focus*SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & focus*SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S1) !d(y)/d(s1)
        z(3)=focus*COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & focus*SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)+&
          & focus*SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S1) !d(z)/d(s1)
      CASE(PART_DERIV_S1_S1)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(PART_DERIV_S2)
        z(1)=focus*SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S2)-&
          & focus*COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2) !d(x)/d(s2)
        z(2)=focus*COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S2)+&
          & focus*SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S2)-&
          & focus*SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S2) !d(y)/d(s2)
        z(3)=focus*COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S2)+&
          & focus*SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S2)+&
          & focus*SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S2) !d(z)/d(s2)
      CASE(PART_DERIV_S2_S2)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(PART_DERIV_S1_S2)
        z(1)=focus*(SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2)+&
          & x(1,PART_DERIV_S1)*(COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S2)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2))-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1_S2)-&
          & x(2,PART_DERIV_S1)*(SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S2)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2))) !d2(x)/d(s1)d(s2)
        z(2)=focus*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2)+&
          & x(1,PART_DERIV_S1)*(SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(1,PART_DERIV_S2)+COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(2,PART_DERIV_S2)-COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S2))+SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(2,PART_DERIV_S1_S2)+x(2,PART_DERIV_S1)*(COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(1,PART_DERIV_S2)-SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S2)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S2))-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1_S2)-&
          & x(3,PART_DERIV_S1)*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(1,PART_DERIV_S2)+SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(2,PART_DERIV_S2)+SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S2))) !d2(y)/d(s1)d(s2)
        z(3)=focus*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2)+&
          & x(1,PART_DERIV_S1)*(SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(1,PART_DERIV_S2)+COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(2,PART_DERIV_S2)+COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S2))+SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1_S2)+&
          & x(2,PART_DERIV_S1)*(COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(1,PART_DERIV_S2)-SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(2,PART_DERIV_S2)+SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S2))+SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1_S2)+&
          & x(3,PART_DERIV_S1)*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(1,PART_DERIV_S2)+SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(2,PART_DERIV_S2)-SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S2))) !d2(z)/d(s1)d(s2)
      CASE(PART_DERIV_S3)
        z(1)=focus*SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S3)-&
          & focus*COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S3) !d(x)/d(s3)
        z(2)=focus*COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & focus*SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)-&
          & focus*SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S3) !d(y)/d(s3)
        z(3)=focus*COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & focus*SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)+&
          & focus*SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S3) !d(z)/d(s3)
      CASE(PART_DERIV_S3_S3)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(PART_DERIV_S1_S3)
        z(1)=focus*(SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1_S3)+&
          & x(1,PART_DERIV_S1)*(COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S3))-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1_S3)-&
          & x(2,PART_DERIV_S1)*(SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S3))) !d2(x)/d(s1)d(s3)
        z(2)=focus*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1_S3)+&
          & x(1,PART_DERIV_S1)*(SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S3))+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1_S3)+&
          & x(2,PART_DERIV_S1)*(COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S3))-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1_S3)-&
          & x(3,PART_DERIV_S1)*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))* x(2,PART_DERIV_S3)+&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S3))) !d2(y)/d(s1)d(s3) 
        z(3)=focus*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1_S3)+&
          & x(1,PART_DERIV_S1)*(SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)+&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S3))+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1_S3)+&
          & x(2,PART_DERIV_S1)*(COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S3))+&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1_S3)+&
          & x(3,PART_DERIV_S1)*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S3))) !d2(z)/d(s1)d(s3)
      CASE(PART_DERIV_S2_S3)
        z(1)=focus*(SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S2_S3)+&
          & x(1,PART_DERIV_S2)*(COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S3))-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S2_S3)-&
          & x(2,PART_DERIV_S2)*(SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S3))) !d2(x)/d(s2)d(s3)
        z(2)=focus*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S2_S3)+&
          & x(1,PART_DERIV_S2)*(SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S3))+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S2_S3)+&
          & x(2,PART_DERIV_S2)*(COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S3))-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S2_S3)-&
          & x(3,PART_DERIV_S2)*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)+&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S3))) !d2(y)/d(s2)d(s3)
        z(3)=focus*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S2_S3)+&
          & x(1,PART_DERIV_S2)*(SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)+&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S3))+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S2_S3)+&
          & x(2,PART_DERIV_S2)*(COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S3))+&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S2_S3)+&
          & x(3,PART_DERIV_S2)*(COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S3)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S3))) !d2(z)/d(s2)d(s3)
      CASE(PART_DERIV_S1_S2_S3)
        z(1)=focus*((SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1))*x(1,PART_DERIV_S2)*x(1,PART_DERIV_S3)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*(x(1,PART_DERIV_S1_S2)*x(1,PART_DERIV_S3)+x(1,PART_DERIV_S2)* &
          & x(1,PART_DERIV_S1_S3))+(-COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1))*x(2,PART_DERIV_S2)*x(1,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*(x(2,PART_DERIV_S1_S2)*x(1,PART_DERIV_S3)+x(2,PART_DERIV_S2)* &
          & x(1,PART_DERIV_S1_S3))+(COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1))*x(1,PART_DERIV_S2_S3)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2_S3)+&
          & (-COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1))*x(1,PART_DERIV_S2)*x(2,PART_DERIV_S3)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*(x(1,PART_DERIV_S1_S2)*x(2,PART_DERIV_S3)+x(1,PART_DERIV_S2)* &
          & x(2,PART_DERIV_S1_S3))+(-SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1))*x(2,PART_DERIV_S2)*x(2,PART_DERIV_S3)-&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*(x(2,PART_DERIV_S1_S2)*x(2,PART_DERIV_S3)+x(2,PART_DERIV_S2)* &
          & x(2,PART_DERIV_S1_S3))+(-SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1))*x(2,PART_DERIV_S2_S3)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*x(2,PART_DERIV_S1_S2_S3)) !d3(x)/d(s1)d(s2)d(s3)
        z(2)=focus*((COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(1,PART_DERIV_S2)*x(1,PART_DERIV_S3)+SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & (x(1,PART_DERIV_S1_S2)*x(1,PART_DERIV_S3)+x(1,PART_DERIV_S2)*x(1,PART_DERIV_S1_S3))+&
          & (SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(2,PART_DERIV_S2)*x(1,PART_DERIV_S3)+COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & (x(2,PART_DERIV_S1_S2)*x(1,PART_DERIV_S3)+x(2,PART_DERIV_S2)*x(1,PART_DERIV_S1_S3))+&
          & (-SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(3,PART_DERIV_S2)*x(1,PART_DERIV_S3)-COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & (x(3,PART_DERIV_S1_S2)*x(1,PART_DERIV_S3)+x(3,PART_DERIV_S2)*x(1,PART_DERIV_S1_S3))+&
          & (SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*x(1,PART_DERIV_S2_S3)+&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2_S3))+&
          & focus*((SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(1,PART_DERIV_S2)*x(2,PART_DERIV_S3)+COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & (x(1,PART_DERIV_S1_S2)*x(2,PART_DERIV_S3)+x(1,PART_DERIV_S2)*x(2,PART_DERIV_S1_S3))+&
          & (-COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)+&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(2,PART_DERIV_S2)*x(2,PART_DERIV_S3)-SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & (x(2,PART_DERIV_S1_S2)*x(2,PART_DERIV_S3)+x(2,PART_DERIV_S2)*x(2,PART_DERIV_S1_S3))+&
          & (-COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(3,PART_DERIV_S2)*x(2,PART_DERIV_S3)-SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & (x(3,PART_DERIV_S1_S2)*x(2,PART_DERIV_S3)+x(3,PART_DERIV_S2)*x(2,PART_DERIV_S1_S3))+&
          & (COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*x(2,PART_DERIV_S2_S3)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1_S2_S3))+&
          & focus*((-SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(1,PART_DERIV_S2)*x(3,PART_DERIV_S3)-COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & (x(1,PART_DERIV_S1_S2)*x(3,PART_DERIV_S3)+x(1,PART_DERIV_S2)*x(3,PART_DERIV_S1_S3))+&
          & (-COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(2,PART_DERIV_S2)*x(3,PART_DERIV_S3)-SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & (x(2,PART_DERIV_S1_S2)*x(3,PART_DERIV_S3)+x(2,PART_DERIV_S2)*x(3,PART_DERIV_S1_S3))+&
          & (-COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)+&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(3,PART_DERIV_S2)*x(3,PART_DERIV_S3)-SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & (x(3,PART_DERIV_S1_S2)*x(3,PART_DERIV_S3)+x(3,PART_DERIV_S2)*x(3,PART_DERIV_S1_S3))+&
          & (-COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*x(3,PART_DERIV_S2_S3)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S1_S2_S3)) !d3(y)/d(s1)d(s2)d(s3)
        z(3)=focus*((COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)+&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(1,PART_DERIV_S2)*x(1,PART_DERIV_S3)+SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & (x(1,PART_DERIV_S1_S2)*x(1,PART_DERIV_S3)+x(1,PART_DERIV_S2)*x(1,PART_DERIV_S1_S3))+&
          & (SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(2,PART_DERIV_S2)*x(1,PART_DERIV_S3)+COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & (x(2,PART_DERIV_S1_S2)*x(1,PART_DERIV_S3)+x(2,PART_DERIV_S2)*x(1,PART_DERIV_S1_S3))+&
          & (SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(3,PART_DERIV_S2)*x(1,PART_DERIV_S3)+COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & (x(3,PART_DERIV_S1_S2)*x(1,PART_DERIV_S3)+x(3,PART_DERIV_S2)*x(1,PART_DERIV_S1_S3))+&
          & (SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)+&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*x(1,PART_DERIV_S2_S3)+&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1_S2_S3))+&
          & focus*((SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(1,PART_DERIV_S2)*x(2,PART_DERIV_S3)+COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & (x(1,PART_DERIV_S1_S2)*x(2,PART_DERIV_S3)+x(1,PART_DERIV_S2)*x(2,PART_DERIV_S1_S3))+&
          & (-COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(2,PART_DERIV_S2)*x(2,PART_DERIV_S3)-SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & (x(2,PART_DERIV_S1_S2)*x(2,PART_DERIV_S3)+x(2,PART_DERIV_S2)*x(2,PART_DERIV_S1_S3))+&
          & (COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(3,PART_DERIV_S2)*x(2,PART_DERIV_S3)+SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & (x(3,PART_DERIV_S1_S2)*x(2,PART_DERIV_S3)+x(3,PART_DERIV_S2)*x(2,PART_DERIV_S1_S3))+&
          & (COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*x(2,PART_DERIV_S2_S3)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1_S2_S3))+&
          & focus*((SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(1,PART_DERIV_S2)*x(3,PART_DERIV_S3)+COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & (x(1,PART_DERIV_S1_S2)*x(3,PART_DERIV_S3)+x(1,PART_DERIV_S2)*x(3,PART_DERIV_S1_S3))+&
          & (COSH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(2,PART_DERIV_S2)*x(3,PART_DERIV_S3)+SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & (x(2,PART_DERIV_S1_S2)*x(3,PART_DERIV_S3)+x(2,PART_DERIV_S2)*x(3,PART_DERIV_S1_S3))+&
          & (-COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*&
          & x(3,PART_DERIV_S2)*x(3,PART_DERIV_S3)-SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*&
          & (x(3,PART_DERIV_S1_S2)*x(3,PART_DERIV_S3)+x(3,PART_DERIV_S2)*x(3,PART_DERIV_S1_S3))+&
          & (COSH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(1,PART_DERIV_S1)+&
          & SINH(x(1,NO_PART_DERIV))*COS(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*x(2,PART_DERIV_S1)-&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*SIN(x(3,NO_PART_DERIV))*x(3,PART_DERIV_S1))*x(3,PART_DERIV_S2_S3)+&
          & SINH(x(1,NO_PART_DERIV))*SIN(x(2,NO_PART_DERIV))*COS(x(3,NO_PART_DERIV))*&
          & x(3,PART_DERIV_S1_S2_S3)) !d3(z)/d(s1)d(s2)d(s3)
      CASE DEFAULT
        localError="The specified partial derivative type of "// &
          & TRIM(NumberToVString(partialDerivativeType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      IF(coordinateSystem%numberOfDimensions/=3) THEN
        localError="The number of coordinate system dimensions of "// &
          & TRIM(NumberToVString(coordinateSystem%numberOfDimensions,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      focus=REAL(coordinateSystem%focus,SP)
      SELECT CASE(partialDerivativeType)
      CASE(NO_PART_DERIV)
        z=CoordinateSystem_ConvertToRC(coordinateSystem,x(:,NO_PART_DERIV),err,error)
        IF(err/=0) GOTO 999
      CASE(PART_DERIV_S1,PART_DERIV_S1_S1,PART_DERIV_S2,PART_DERIV_S2_S2,&
        & PART_DERIV_S1_S2,PART_DERIV_S3,PART_DERIV_S3_S3,&
        & PART_DERIV_S1_S3,PART_DERIV_S2_S3,PART_DERIV_S1_S2_S3)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The specified partial derivative type of "// &
          & TRIM(NumberToVString(partialDerivativeType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The specified coordinate system type of "// &
        & TRIM(NumberToVString(coordinateSystem%TYPE,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("CoordinateSystem_DerivativeConvertToRCSP")
    RETURN
999 ERRORSEXITS("CoordinateSystem_DerivativeConvertToRCSP",err,error)
    RETURN 1
  END SUBROUTINE CoordinateSystem_DerivativeConvertToRCSP

  !
  !================================================================================================================================
  !

  !>Calculates the norm of a derivative in a coordinate system identified by coordinateSystem at the given interpolated
  !>point and returns the value in NORM for single precision coordinates. partialDerivativeIndex is used to select the
  !>appropriate partial derivative (i.e., wrt S1, S2 or S3) to normalise.
  SUBROUTINE CoordinateSystem_DerivativeNorm(coordinateSystem,partialDerivativeIndex,interpolatedPoint,derivativeNorm,err,error,*)
  
    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<A pointer to the coordinate system to calculate the derivative norm for
    INTEGER(INTG), INTENT(IN) :: partialDerivativeIndex !<The partial derivative index to select the direction to normalise
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint !<A pointer to the interpolated point 
    REAL(DP), INTENT(OUT) :: derivativeNorm !<On exit, the derivative norm of the coordinate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: componentIdx,numberOfComponents
    REAL(DP) :: focus,sl,sm
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("CoordinateSystem_DerivativeNorm",err,error,*999)

    derivativeNorm=0.0_DP
    
    IF(.NOT.ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(interpolatedPoint)) CALL FlagError("Interpolated point is not associated.",err,error,*999)
    IF(interpolatedPoint%partialDerivativeType<FIRST_PART_DERIV) &
      & CALL FlagError("The point has not been interpolated to include first derivative values.",err,error,*999)      
    IF(partialDerivativeIndex>interpolatedPoint%maximumPartialDerivativeIndex) THEN
      localError="The partial derivative index of "//TRIM(NumberToVString(partialDerivativeIndex,"*",err,error))// &
        & " is invalid. The interpolated point has a maximum number of partial derivatives of "// &
        & TRIM(NumberToVString(interpolatedPoint%maximumPartialDerivativeIndex,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    END IF
    
    numberOfComponents=SIZE(interpolatedPoint%values,1)
    SELECT CASE(partialDerivativeIndex)
    CASE(PART_DERIV_S1,PART_DERIV_S2,PART_DERIV_S3)
      SELECT CASE(coordinateSystem%TYPE)
      CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
        DO componentIdx=1,numberOfComponents
          derivativeNorm=derivativeNorm+interpolatedPoint%values(componentIdx,partialDerivativeIndex)**2
        ENDDO !component_index
      CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
        IF(numberOfComponents==2) THEN
          derivativeNorm=interpolatedPoint%values(1,partialDerivativeIndex)**2+(interpolatedPoint%values(1,1)* &
            & interpolatedPoint%values(2,partialDerivativeIndex))**2
        ELSE IF(numberOfComponents==3) THEN
          derivativeNorm=interpolatedPoint%values(1,partialDerivativeIndex)**2+(interpolatedPoint%values(1,1)* &
            & interpolatedPoint%values(2,partialDerivativeIndex))**2+interpolatedPoint%values(3,partialDerivativeIndex)**2
        ELSE
          localError="The number of components for the interpolated point of "// &
            & TRIM(NumberToVString(numberOfComponents,"*",err,error))// &
            & " is invalid for a cylindrical polar coordinate system."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
        derivativeNorm=interpolatedPoint%values(1,partialDerivativeIndex)**2+(interpolatedPoint%values(1,1)* &
          & interpolatedPoint%values(2,partialDerivativeIndex)*COS(interpolatedPoint%values(3,1)))**2+ &
          & (interpolatedPoint%values(1,1)*interpolatedPoint%values(3,partialDerivativeIndex))**2
      CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
        focus=coordinateSystem%focus
        sl=SINH(interpolatedPoint%values(1,1))
        sm=SIN(interpolatedPoint%values(2,1))
        derivativeNorm=focus*focus*((sl*sl+sm*sm)*(interpolatedPoint%values(1,partialDerivativeIndex)**2+ &
          & interpolatedPoint%values(2,partialDerivativeIndex))**2)+(sl*sm*interpolatedPoint%values(3,partialDerivativeIndex))**2
      CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The coordinate system type of "//TRIM(NumberToVString(coordinateSystem%type,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      derivativeNorm=SQRT(derivativeNorm)
    CASE DEFAULT
      localError="The partial derivative index of "//TRIM(NumberToVString(partialDerivativeIndex,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("CoordinateSystem_DerivativeNorm")
    RETURN
999 ERRORSEXITS("CoordinateSystem_DerivativeNorm",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_DerivativeNorm

  !
  !================================================================================================================================
  !

  !>Adjusts the interpolation for non-rectangular cartesian coordinate systems.
  SUBROUTINE CoordinateSystem_InterpolationAdjust(coordinateSystem,partialDerivativeIndex,value,err,error,*)
  
    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<A pointer to the coordinate system to adjust
    INTEGER(INTG), INTENT(IN) :: partialDerivativeIndex !<The partial derivative index to adjust
    REAL(DP), INTENT(INOUT) :: value !<On entry, the coordinate value to adjust. On exit, the adjusted value.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    REAL(DP) :: coshx,css,d,des,focus,r,ss,sinhx,theta
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("CoordinateSystem_InterpolationAdjust",err,error,*999)

    IF(.NOT.ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is not associated.",err,error,*999)
    
    SELECT CASE(coordinateSystem%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      !Do nothing
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      SELECT CASE(coordinateSystem%radialInterpolationType)
      CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
        !Do nothing
      CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
        r=SQRT(value)
        IF(partialDerivativeIndex==NO_PART_DERIV) THEN
          value=r
        ELSE
          value=value/(2.0_DP*R)
        ENDIF
      CASE DEFAULT
        localError="The radial interpolation type of "//TRIM(NumberToVString(coordinateSystem% &
          & radialInterpolationType,"*",err,error))//" is invalid for a cylindrical coordinate system."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      SELECT CASE(coordinateSystem%radialInterpolationType)
      CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
        !Do nothing
      CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
        r=value**(1.0_DP/3.0_DP)
        IF(partialDerivativeIndex==NO_PART_DERIV) THEN
          value=r
        ELSE
          value=value/(3.0_DP*r*r)
        ENDIF
      CASE DEFAULT
        localError="The radial interpolation type of "//TRIM(NumberToVString(coordinateSystem% &
          & radialInterpolationType,"*",err,error))//" is invalid for a cylindrical/spherical coordinate system."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      SELECT CASE(coordinateSystem%radialInterpolationType)
      CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
        !Do nothing
      CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
        focus=coordinateSystem%focus
        ss=value/(focus*focus)
        sinhx=SQRT(ss)
        coshx=SQRT(1.0_DP+ss)
        IF(partialDerivativeIndex==NO_PART_DERIV) THEN
          value=LOG(sinhx+coshx)
        ELSE
          value=value/(2.0_DP*focus*focus*sinhx*coshx)
        ENDIF
      CASE(COORDINATE_RADIAL_CUBED_INTERPOLATION_TYPE)
        focus=coordinateSystem%focus
        css=value/(focus**3.0_DP)
        des=css*css-4.0_DP/27.0_DP
        IF(des>0.0_DP) THEN
          d=((css+SQRT(des))/2.0_DP)**(1.0_DP/3.0_DP)
          coshx=d+1.0_DP/(3.0_DP*D)
        ELSE
          theta=ACOS(css*SQRT(27.0_DP)/2.0_DP)
          coshx=2.0_DP/SQRT(3.0_DP)*COS(theta/3.0_DP)
        ENDIF
        sinhx=SQRT(ABS(coshx*coshx-1.0_DP))
        IF(partialDerivativeIndex==NO_PART_DERIV) THEN
          value=LOG(sinhx+coshx)
        ELSE
          value=value/((3.0_DP*coshx*coshx-1.0_DP)*sinhx)/(focus**3.0_DP)
        ENDIF
      CASE DEFAULT
        localError="The radial interpolation type of "//TRIM(NumberToVString(coordinateSystem% &
          & radialInterpolationType,"*",err,error))//" is invalid for a prolate spheroidal coordinate system."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The coordinate system type of "//TRIM(NumberToVString(coordinateSystem%type,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
     
    EXITS("CoordinateSystem_InterpolationAdjust")
    RETURN
999 ERRORSEXITS("CoordinateSystem_InterpolationAdjust",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_InterpolationAdjust

  !
  !================================================================================================================================
  !

  !>Adjusts the interpolation parameters for non-rectangular cartesian coordinate systems.
  SUBROUTINE CoordinateSystem_InterpolationParametersAdjust(coordinateSystem,interpolationParameters,err,error,*)
  
    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<A pointer to the coordinate system to adjust
    TYPE(FieldInterpolationParametersType), POINTER :: interpolationParameters !<A pointer to the interpolation parameters to adjust
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("CoordinateSystem_InterpolationParametersAdjust",err,error,*999)

!!TODO: Tidy up element parameters for non-rc coordinate systems. See bottom of XPXE and ZPZE.
    
    IF(.NOT.ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(interpolationParameters)) CALL FlagError("Interpolation parameters is not associated.",err,error,*999)
    
    SELECT CASE(coordinateSystem%TYPE)
    CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
      !Do nothing
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE)
      SELECT CASE(coordinateSystem%radialInterpolationType)
      CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
        !Do nothing
      CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
      CASE DEFAULT
        localError="The radial interpolation type of "//TRIM(NumberToVString(coordinateSystem% &
          & radialInterpolationType,"*",err,error))//" is invalid for a cylindrical coordinate system."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      CALL FlagError("Not implemented",err,error,*999)
    CASE(COORDINATE_SPHERICAL_POLAR_TYPE)
      SELECT CASE(coordinateSystem%radialInterpolationType)
      CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
        !Do nothing
      CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
      CASE DEFAULT
        localError="The radial interpolation type of "//TRIM(NumberToVString(coordinateSystem% &
          & radialInterpolationType,"*",err,error))//" is invalid for a spherical coordinate system."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE)
      SELECT CASE(coordinateSystem%radialInterpolationType)
      CASE(COORDINATE_RADIAL_INTERPOLATION_TYPE)
        !Do nothing
      CASE(COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE)
      CASE(COORDINATE_RADIAL_CUBED_INTERPOLATION_TYPE)
      CASE DEFAULT
        localError="The radial interpolation type of "//TRIM(NumberToVString(coordinateSystem% &
          & radialInterpolationType,"*",err,error))//" is invalid for a prolate spheroidal coordinate system."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The coordinate system type of "//TRIM(NumberToVString(coordinateSystem%type,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    EXITS("CoordinateSystem_InterpolationParametersAdjust")
    RETURN
999 ERRORS("CoordinateSystem_InterpolationParametersAdjust",err,error)
    EXITS("CoordinateSystem_InterpolationParametersAdjust")
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_InterpolationParametersAdjust

  !
  !================================================================================================================================
  !
 
  !>Calculates the tensor to get from material coordinate system, nu, to local coordinate system, xi.
  SUBROUTINE CoordinateSystem_MaterialSystemCalculate(geometricInterpPointMetrics,fibreInterpPoint,dNudX,dXdNu,dNudXi,dXidNu, &
    & err,error,*)
  
    !Argument variables
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics !<The geometric interpolation point metrics at the point to calculate the material coordinate system from.
    TYPE(FieldInterpolatedPointType), POINTER :: fibreInterpPoint !<The fibre interpolation point at the point to calculate the material coordinate system from
    REAL(DP), INTENT(OUT) :: dNudx(:,:) !<dNudx(nuIdx,xIdx). On return, the tensor to transform from the material system to the geometric coordinate system
    REAL(DP), INTENT(OUT) :: dXdNu(:,:) !<dXdNu(xIdx,nuIdx). On return, the tensor to transform from the geometric coordinate system to the material coordinate system
    REAL(DP), INTENT(OUT) :: dNudXi(:,:) !<dNudXi(nuIdx,xiIdx). On return, the tensor to transform from the material system to the xi coordinate system
    REAL(DP), INTENT(OUT) :: dXidNu(:,:) !<dXidNu(xiIdx,nuIdx). On return, the tensor to transform from the xi coordinate system to the material coordinate system
    INTEGER(INTG), INTENT(OUT) :: err   !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) ::  error   !<The error string
    !Local variables
    INTEGER(INTG) :: numberOfXDimensions,numberOfXiDimensions,numberOfNuDimensions
    REAL(DP) :: dNudXiTemp(3,3),Jnuxi,JNuX
    TYPE(VARYING_STRING) :: localError 
     
    ENTERS("CoordinateSystem_MaterialSystemCalculate",err,error,*999)
    
    IF(.NOT.ASSOCIATED(geometricInterpPointMetrics)) &
      & CALL FlagError("Geometric interpolated point metrics is not associated.",err,error,*999)
      
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
        CALL CoordinateSystem_MaterialSystemCalculatedXdNu2D(geometricInterpPointMetrics,fibreInterpPoint%values(1: &
          & numberOfNuDimensions,1),dXdNu(1:numberOfXDimensions,1:numberOfXDimensions),err,error,*999)
      CASE(3)
        CALL CoordinateSystem_MaterialSystemCalculatedXdNu3D(geometricInterpPointMetrics,fibreInterpPoint%values(1: &
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
      CALL IdentityMatrix(dXdNu,err,error,*999)
    ENDIF
    !Calculate dNu/dX the inverse of dX/dNu (same as transpose due to orthogonality)
    CALL MatrixTranspose(dXdNu(1:numberOfXDimensions,1:numberOfXDimensions),dNudX(1:numberOfXDimensions,1: &
      & numberOfXDimensions),err,error,*999)
    !Calculate dNu/dXi = dNu/dX * dX/dXi and its inverse dXi/dNu
    CALL MatrixProduct(dNudx(1:numberOfXDimensions,1:numberOfXDimensions), &
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

    IF(diagnostics1) THEN
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
      CALL Determinant(dNudx(1:numberOfXDimensions,1:numberOfXDimensions),JNuX,err,error,*999)
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
   
    EXITS("CoordinateSystem_MaterialSystemCalculate")
    RETURN
999 ERRORSEXITS("CoordinateSystem_MaterialSystemCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_MaterialSystemCalculate
  
  !
  !================================================================================================================================
  !

  !>Calculates transformation between spatial CS and rotated reference orthogonal material CS in 2D space
  SUBROUTINE CoordinateSystem_MaterialSystemCalculatedXdNu2D(geometricInterpPointMetrics,angle,dXdNu,err,error,*)

    !Argument variables
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics !<The geometric interpolated point metrics at the point to calculate dXdNu at. 
    REAL(DP), INTENT(IN) :: angle(:) !<angle(fibreIdx). The fibre angle (in radians) 
    REAL(DP), INTENT(OUT) :: dXdNu(:,:) !<dXdNu(coordinateIdx,coordinateIdx). On exit, the dXdNu tensor.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: det,dXdNuR(2,2),R(2,2),f(2),g(2)

    ENTERS("CoordinateSystem_MaterialSystemCalculatedXdNu2D",err,error,*999)

    IF(.NOT.ASSOCIATED(geometricInterpPointMetrics)) &
      & CALL FlagError("Geometry interpolated point metrics is not associated.",err,error,*999)
    
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
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,2,2,f,'("    f              :",2(x,E13.6))','(20X,2(x,E13.6))', &
        & err,error,*999)      
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,2,2,g,'("    g              :",2(x,E13.6))','(20X,2(x,E13.6))', &
        & err,error,*999)      
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Derivative of X wrt Nu (reference):",err,error,*999)
      CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,1,1,2,2,2,dXdNuR,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
        & '("      dX_dNuR','(",I1,",:)',' :",2(x,E13.6))','(20X,2(x,E13.6))',err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Fibre calculation:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Fibre angle = ",angle(1),err,error,*999)
      IF(diagnostics2) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Rotation matrix:",err,error,*999)
        CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,1,1,2,2,2,R,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
          & '("      R','(",I1,",:)','       :",2(x,E13.6))','(20X,2(x,E13.6))',err,error,*999)
      ENDIF
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Derivative of X wrt Nu (material):",err,error,*999)
      CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,1,1,2,2,2,dXdNu,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
        & '("      dX_dNu','(",I1,",:)','   :",2(x,E13.6))','(20X,2(x,E13.6))',err,error,*999)
      CALL Determinant(dXdNu,det,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Determinant dX_dNu = ",det,err,error,*999)
    ENDIF
        
    EXITS("CoordinateSystem_MaterialSystemCalculatedXdNu2D")
    RETURN
999 ERRORS("CoordinateSystem_MaterialSystemCalculatedXdNu2D",err,error)
    EXITS("CoordinateSystem_MaterialSystemCalculatedXdNu2D")
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_MaterialSystemCalculatedXdNu2D

  !
  !================================================================================================================================
  !

  !>Calculates transformation between spatial CS and rotated reference orthogonal material CS in 3D space
  SUBROUTINE CoordinateSystem_MaterialSystemCalculatedXdNu3D(geometricInterpPointMetrics,angle,dXdNu,err,error,*)

    !Argument variables
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics !<The geometric interpolated point metrics at the point to calculate dXdNu at. 
    REAL(DP), INTENT(IN) :: angle(:) !<angles(fibreIdx). The fibre, imbrication and sheet (in radians) 
    REAL(DP), INTENT(OUT) :: dXdNu(:,:) !<dXdNu(coordinateIdx,coordinateIdx). On exit, the dXdNu tensor.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: angles(3),det,dXdNu2(3,3),dXdNu3(3,3),dXdNuR(3,3),f(3),g(3),h(3), &
      & Raf(3,3),Rbf(3,3),Rai(3,3),Rbi(3,3),Ras(3,3),Rbs(3,3)
    
    ENTERS("CoordinateSystem_MaterialSystemCalculatedXdNu3D",err,error,*999)
    
    IF(.NOT.ASSOCIATED(geometricInterpPointMetrics)) &
      & CALL FlagError("Geometry interpolated point metrics is not associated.",err,error,*999)

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
    !(a) first align reference material direction 3 with z(spatial) axis by rotating the ref material CS. 
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
    !(a) first align new material direction 1 with x(spatial) axis by rotating the new material CS. 
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
   
    EXITS("CoordinateSystem_MaterialSystemCalculatedXdNu3D")
    RETURN
999 ERRORS("CoordinateSystem_MaterialSystemCalculatedXdNu3D",err,error)
    EXITS("CoordinateSystem_MaterialSystemCalculatedXdNu3D")
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_MaterialSystemCalculatedXdNu3D

  !
  !================================================================================================================================
  !

  !>Transforms a symmetric second order tensor in materials coordinates (e.g., conductivity tensor) to a tensor in xi coordinates
  SUBROUTINE CoordinateSystem_MaterialTransformSymTensor2(geometricInterpPointMetrics,fibreInterpPoint,symmetricMaterialTensor, &
    & transformedTensor,err,error,*)

    !Argument variables
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics !<The geometric interpolated point metrics at the point to tranformed tensor.
    TYPE(FieldInterpolatedPointType), POINTER :: fibreInterpPoint !<The fibre interpolated point (if it exists).
    
    REAL(DP), INTENT(IN) :: symmetricMaterialTensor(:) !<symmetricMaterialTensor(voigtIdx). The original symmetric material tensor values (in Voigt form) to transform.
    REAL(DP), INTENT(OUT) :: transformedTensor(:,:) !<transformedTensor(xiCoordinateIdx,xiCoordinateIdx). On exit, the tensor transformed from material coordinates.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfDimensions,numberOfXi
    REAL(DP) :: dNudX(3,3),dNudXi(3,3),dXdNu(3,3),dXidNu(3,3),materialTensor(3,3),tempTensor(3,3)
    TYPE(VARYING_STRING) :: localError

    ENTERS("CoordinateSystem_MaterialTransformSymTensor2",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(geometricInterpPointMetrics)) &
      & CALL FlagError("Geometry interpolated point metrics is not associated.",err,error,*999)
    IF(geometricInterpPointMetrics%numberOfXDimensions/=geometricInterpPointMetrics%numberOfXiDimensions) THEN
      localError="A different number of x and xi dimensions is not implemented. The number of x dimensions is "// &
        & TRIM(NumberToVString(geometricInterpPointMetrics%numberOfXDimensions,"*",err,error))// &
        & " and the number of xi dimensions is "// &
        & TRIM(NumberToVString(geometricInterpPointMetrics%numberOfXiDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(symmetricMaterialTensor,1)<NUMBER_OF_VOIGT(geometricInterpPointMetrics%numberOfXDimensions)) THEN
      localError="The size of the specified symmetric tensor of "// &
        & TRIM(NumberToVString(SIZE(symmetricMaterialTensor,1),"*",err,error))// &
        & " is too small. The size should be >= "// &
        & TRIM(NumberToVString(NUMBER_OF_VOIGT(geometricInterpPointMetrics%numberOfXDimensions),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(transformedTensor,1)<geometricInterpPointMetrics%numberOfXiDimensions) THEN
      localError="The first dimension of the specified transformed tensor of "// &
        & TRIM(NumberToVString(SIZE(transformedTensor,1),"*",err,error))//" is too small. The size should be >= "// &
        & TRIM(NumberToVString(geometricInterpPointMetrics%numberOfXiDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(transformedTensor,2)<geometricInterpPointMetrics%numberOfXiDimensions) THEN
      localError="The second dimension of the specified transformed tensor of "// &
        & TRIM(NumberToVString(SIZE(transformedTensor,2),"*",err,error))//" is too small. The size should be >= "// &
        & TRIM(NumberToVString(geometricInterpPointMetrics%numberOfXiDimensions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    numberOfDimensions=geometricInterpPointMetrics%numberOfXDimensions
    numberOfXi=geometricInterpPointMetrics%numberOfXiDimensions
    
    !Calculate the untransformed tensor in material coordinates
    materialTensor(1:numberOfDimensions,1:numberOfDimensions)=0.0_DP
    SELECT CASE(numberOfDimensions)
    CASE(1)
      materialTensor(1,1)=symmetricMaterialTensor(TENSOR_TO_VOIGT1(1,1))
    CASE(2)
      materialTensor(1,1)=symmetricMaterialTensor(TENSOR_TO_VOIGT2(1,1))
      materialTensor(1,2)=symmetricMaterialTensor(TENSOR_TO_VOIGT2(1,2))
      materialTensor(2,1)=symmetricMaterialTensor(TENSOR_TO_VOIGT2(2,1))
      materialTensor(2,2)=symmetricMaterialTensor(TENSOR_TO_VOIGT2(2,2))
    CASE(3)
      materialTensor(1,1)=symmetricMaterialTensor(TENSOR_TO_VOIGT3(1,1))
      materialTensor(1,2)=symmetricMaterialTensor(TENSOR_TO_VOIGT3(1,2))
      materialTensor(1,3)=symmetricMaterialTensor(TENSOR_TO_VOIGT3(1,3))
      materialTensor(2,1)=symmetricMaterialTensor(TENSOR_TO_VOIGT3(2,1))
      materialTensor(2,2)=symmetricMaterialTensor(TENSOR_TO_VOIGT3(2,2))
      materialTensor(2,3)=symmetricMaterialTensor(TENSOR_TO_VOIGT3(2,3))
      materialTensor(3,1)=symmetricMaterialTensor(TENSOR_TO_VOIGT3(3,1))
      materialTensor(3,2)=symmetricMaterialTensor(TENSOR_TO_VOIGT3(3,2))
      materialTensor(3,3)=symmetricMaterialTensor(TENSOR_TO_VOIGT3(3,3))
    CASE DEFAULT
      localError="The number of dimensions of "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    !Calculate material coordinates
    dNudX(1:numberOfDimensions,1:numberOfDimensions)=0.0_RP
    dXdNu(1:numberOfDimensions,1:numberOfDimensions)=0.0_RP
    dNudXi(1:numberOfDimensions,1:numberOfXi)=0.0_RP
    dXidNu(1:numberOfXi,1:numberOfDimensions)=0.0_RP
    CALL CoordinateSystem_MaterialSystemCalculate(geometricInterpPointMetrics,fibreInterpPoint, &
      & dNudX(1:numberOfDimensions,1:numberOfDimensions),dXdNu(1:numberOfDimensions,1:numberOfDimensions), &
      & dNudXi(1:numberOfDimensions,1:numberOfXi),dXidNu(1:numberOfXi,1:numberOfDimensions), &
      & err,error,*999)

    !Rotate tensor from material coordinates
    tempTensor(1:numberOfXi,1:numberOfDimensions)=0.0_DP
    CALL MatrixProduct(dXidNu(1:numberOfXi,1:numberOfDimensions),materialTensor(1:numberOfDimensions,1:numberOfDimensions), &
      & tempTensor(1:numberOfXi,1:numberOfDimensions),err,error,*999)
    transformedTensor(1:numberOfXi,1:numberOfXi)=0.0_RP
    CALL MatrixProduct(tempTensor(1:numberOfXi,1:numberOfDimensions),dNudXi(1:numberOfDimensions,1:numberOfXi), &
      & transformedTensor(1:numberOfXi,1:numberOfXi),err,error,*999)
    
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Symmetric material tensor transformation:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of X  = ",numberOfDimensions,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of Xi = ",numberOfXi,err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Material tensor:",err,error,*999)
      CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfDimensions,1,1,numberOfDimensions,numberOfDimensions, &
        & numberOfDimensions,materialTensor,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
        & '("    MT','(",I1,",:)',' :",3(X,E13.6))','(12X,3(X,E13.6))',err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Transformed tensor:",err,error,*999)
      CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfXi,1,1,numberOfXi,numberOfXi, &
        & numberOfXi,transformedTensor,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
        & '("    TT','(",I1,",:)','   :",3(X,E13.6))','(12X,3(X,E13.6))',err,error,*999)
    ENDIF
        
    EXITS("CoordinateSystem_MaterialTransformSymTensor2")
    RETURN
999 ERRORS("CoordinateSystem_MaterialTransformSymTensor2",err,error)
    EXITS("CoordinateSystem_MaterialTransformSymTensor2")
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_MaterialTransformSymTensor2

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
        CALL CoordinateSystem_Finalise(coordinateSystems%coordinateSystems(coordinateSystemIdx)%ptr,err,error,*999)
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
  
END MODULE CoordinateSystemRoutines

