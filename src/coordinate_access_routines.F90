!> \file
!> \author Chris Bradley
!> \brief This module contains all coordinate system access method routines.
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

!> \addtogroup OpenCMISS_CoordinateSystem OpenCMISS::Iron::CoordinateSystem
!> This module contains all coordinate system access method routines.
MODULE CoordinateSystemAccessRoutines
  
  USE BaseRoutines
  USE Kinds
  USE ISO_VARYING_STRING
  USE Strings
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup OpenCMISS_CoordinateSystemConstants OpenCMISS::Iron::CoordinateSystem::Constants
  !> \see CoordinateRoutines
  !> Coordinate system type parameters
  !>@{ 
  !> \addtogroup CoordinateSystem_CoordinateSystemTypes OpenCMISS::Iron::CoordinateSystem::Constants::CoordinateSystemTypes
  !> \see CoordinateRoutines
  !> Coordinate system type parameters
  !>@{ 
  INTEGER(INTG), PARAMETER :: COORDINATE_RECTANGULAR_CARTESIAN_TYPE=1 !<Rectangular Cartesian coordinate system type \see CoordinateRoutines_CoordinateSystemTypes,CoordinateRoutines
  INTEGER(INTG), PARAMETER :: COORDINATE_CYLINDRICAL_POLAR_TYPE=2 !<Cylindrical polar coordinate system type \see CoordinateRoutines_CoordinateSystemTypes,CoordinateRoutines
  INTEGER(INTG), PARAMETER :: COORDINATE_SPHERICAL_POLAR_TYPE=3 !<Spherical polar coordinate system type \see CoordinateRoutines_CoordinateSystemTypes,CoordinateRoutines
  INTEGER(INTG), PARAMETER :: COORDINATE_PROLATE_SPHEROIDAL_TYPE=4 !<Prolate spheroidal coordinate system type \see CoordinateRoutines_CoordinateSystemTypes,CoordinateRoutines
  INTEGER(INTG), PARAMETER :: COORDINATE_OBLATE_SPHEROIDAL_TYPE=5 !<Oblate spheroidal coordinate system type \see CoordinateRoutines_CoordinateSystemTypes,CoordinateRoutines
  !>@}

  !> \addtogroup CoordinateRoutines_RadialInterpolations  OpenCMISS::Iron::CoordinateSystem::Constants::RadialInterpolations
  !> \see CoordinateRoutines
  !> \brief The type of radial interpolation for polar coordinate systems
  !>@{
  INTEGER(INTG), PARAMETER :: COORDINATE_NO_RADIAL_INTERPOLATION_TYPE=0 !<No radial interpolation \see CoordinateRoutines_RadialInterpolations,CoordinateRoutines
  INTEGER(INTG), PARAMETER :: COORDINATE_RADIAL_INTERPOLATION_TYPE=1 !<r radial interpolation \see CoordinateRoutines_RadialInterpolations,CoordinateRoutines
  INTEGER(INTG), PARAMETER :: COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE=2 !<r^2 radial interpolation \see CoordinateRoutines_RadialInterpolations,CoordinateRoutines
  INTEGER(INTG), PARAMETER :: COORDINATE_RADIAL_CUBED_INTERPOLATION_TYPE=3 !<r^3 radial interpolation \see CoordinateRoutines_RadialInterpolations,CoordinateRoutines
  !>@}
  
  !> \addtogroup CoordinateRoutines_JacobianType OpenCMISS::Iron::CoordinateSystem::Constants::JacobianType
  !> \see CoordinateRoutines
  !> \brief The type of Jacobian to return when coordinate metrics are calculated.
  !>@{
  INTEGER(INTG), PARAMETER :: COORDINATE_JACOBIAN_NO_TYPE=0 !<No Jacobian \see CoordinateRoutines_JacobianTypes,CoordinateRoutines
  INTEGER(INTG), PARAMETER :: COORDINATE_JACOBIAN_LINE_TYPE=1 !<Line type Jacobian \see CoordinateRoutines_JacobianTypes,CoordinateRoutines
  INTEGER(INTG), PARAMETER :: COORDINATE_JACOBIAN_AREA_TYPE=2 !<Area type Jacobian \see CoordinateRoutines_JacobianTypes,CoordinateRoutines
  INTEGER(INTG), PARAMETER :: COORDINATE_JACOBIAN_VOLUME_TYPE=3 !<Volume type Jacobian \see CoordinateRoutines_JacobianTypes,CoordinateRoutines
  !>@}
  !>@}
  
  !Module types

  !Module variables
  
  CHARACTER(LEN=21) :: COORDINATE_SYSTEM_TYPE_STRING(5) = &
    & [ "Rectangular Cartesian", &
    &    "Cylindrical Polar    ", &
    &    "Spherical Polar      ", &
    &    "Prolate Spheroidal   ", &
    &    "Oblate Spheroidal    " ]

  !Interfaces

  PUBLIC COORDINATE_RECTANGULAR_CARTESIAN_TYPE,COORDINATE_CYLINDRICAL_POLAR_TYPE,COORDINATE_SPHERICAL_POLAR_TYPE, &
    & COORDINATE_PROLATE_SPHEROIDAL_TYPE,COORDINATE_OBLATE_SPHEROIDAL_TYPE

  PUBLIC COORDINATE_NO_RADIAL_INTERPOLATION_TYPE,COORDINATE_RADIAL_INTERPOLATION_TYPE, &
    & COORDINATE_RADIAL_SQUARED_INTERPOLATION_TYPE,COORDINATE_RADIAL_CUBED_INTERPOLATION_TYPE
  
  PUBLIC COORDINATE_JACOBIAN_NO_TYPE,COORDINATE_JACOBIAN_LINE_TYPE,COORDINATE_JACOBIAN_AREA_TYPE,COORDINATE_JACOBIAN_VOLUME_TYPE

  PUBLIC COORDINATE_SYSTEM_TYPE_STRING
  
  PUBLIC CoordinateSystem_AssertIsFinished,CoordinateSystem_AssertNotFinished

  PUBLIC CoordinateSystem_CoordinateSystemsGet
  
  PUBLIC CoordinateSystem_DimensionGet

  PUBLIC CoordinateSystem_FocusGet

  PUBLIC CoordinateSystem_Get

  PUBLIC CoordinateSystem_OriginGet

  PUBLIC CoordinateSystem_OrientationGet

  PUBLIC CoordinateSystem_RadialInterpolationTypeGet

  PUBLIC CoordinateSystem_TypeGet
  
  PUBLIC CoordinateSystem_UserNumberFind  

  PUBLIC CoordinateSystems_WorldCoordinateSystemGet
  
CONTAINS

  !
  !=================================================================================================================================
  !

  !>Assert that a coordinate system has been finished
  SUBROUTINE CoordinateSystem_AssertIsFinished(coordinateSystem,err,error,*)

    !Argument Variables
    TYPE(CoordinateSystemType), POINTER, INTENT(INOUT) :: coordinateSystem !<The coordinate system to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("CoordinateSystem_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is not associated.",err,error,*999)
#endif    

    IF(.NOT.coordinateSystem%coordinateSystemFinished) THEN
      localError="Coordinate system number "//TRIM(NumberToVString(coordinateSystem%userNumber,"*",err,error))// &
        & " has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("CoordinateSystem_AssertIsFinished")
    RETURN
999 ERRORSEXITS("CoordinateSystem_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a coordinate system has not been finished
  SUBROUTINE CoordinateSystem_AssertNotFinished(coordinateSystem,err,error,*)

    !Argument Variables
    TYPE(CoordinateSystemType), POINTER, INTENT(INOUT) :: coordinateSystem !<The coordinate system to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("CoordinateSystem_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is not associated.",err,error,*999)
#endif    

    IF(coordinateSystem%coordinateSystemFinished) THEN
      localError="Coordinate system number "//TRIM(NumberToVString(coordinateSystem%userNumber,"*",err,error))// &
        & " has already been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("CoordinateSystem_AssertNotFinished")
    RETURN
999 ERRORSEXITS("CoordinateSystem_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Returns the coordinate systems for the coordinate system. 
  SUBROUTINE CoordinateSystem_CoordinateSystemsGet(coordinateSystem,coordinateSystems,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER, INTENT(IN) :: coordinateSystem !<A pointer to the coordinate system to get the coordinate systems for
    TYPE(CoordinateSystemsType), POINTER, INTENT(OUT) :: coordinateSystems !<On return, a pointer to the coordinate systems for the coordinate system. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CoordinateSystem_CoordinateSystemsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(coordinateSystems)) CALL FlagError("Coordinate systems is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is not associated.",err,error,*999)
#endif    

    coordinateSystems=>coordinateSystem%coordinateSystems

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(coordinateSystems)) THEN
      localError="Coordinate systems is not associated for coordinate system number "// &
        & TRIM(NumberToVString(coordinateSystem%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("CoordinateSystem_CoordinateSystemsGet")
    RETURN
998 NULLIFY(coordinateSystems)
999 ERRORS("CoordinateSystem_CoordinateSystemsGet",err,error)
    EXITS("CoordinateSystem_CoordinateSystemsGet")
    RETURN 1

  END SUBROUTINE CoordinateSystem_CoordinateSystemsGet

  !
  !================================================================================================================================
  !

  !>Gets the coordinate system dimension. 
  SUBROUTINE CoordinateSystem_DimensionGet(coordinateSystem,numberOfDimensions,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<A pointer to the coordinate system to get the dimension for
    INTEGER(INTG), INTENT(OUT) :: numberOfDimensions !<On return, the number of dimensions in the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CoordinateSystem_DimensionGet",err,error,*999)

    CALL CoordinateSystem_AssertIsFinished(coordinateSystem,err,error,*999)    

    numberOfDimensions=coordinateSystem%numberOfDimensions
    
    EXITS("CoordinateSystem_DimensionGet")
    RETURN
999 ERRORSEXITS("CoordinateSystem_DimensionGet",err,error)
    RETURN 1

  END SUBROUTINE CoordinateSystem_DimensionGet

  !
  !================================================================================================================================
  !

  !>Returns the coordinate system focus. 
  SUBROUTINE CoordinateSystem_FocusGet(coordinateSystem,focus,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<A pointer to the coordinate system to get the focus for
    REAL(DP), INTENT(OUT) :: focus !<On return, the focus of the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("CoordinateSystem_FocusGet",err,error,*999)

    CALL CoordinateSystem_AssertIsFinished(coordinateSystem,err,error,*999)
    
    SELECT CASE(coordinateSystem%type)
    CASE(COORDINATE_PROLATE_SPHEROIDAL_TYPE,COORDINATE_OBLATE_SPHEROIDAL_TYPE)
      focus=coordinateSystem%focus
    CASE DEFAULT
      localError="Focus is not defined for coordinate system type "// &
        & TRIM(NumberToVString(coordinateSystem%TYPE,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("CoordinateSystem_FocusGet")
    RETURN
999 ERRORSEXITS("CoordinateSystem_FocusGet",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_FocusGet

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the coordinate system with the given user number. 
  SUBROUTINE CoordinateSystem_Get(coordinateSystems,userNumber,coordinateSystem,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemsType), POINTER :: coordinateSystems !<A pointer to the coordinate systems to get the user number for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the coordinate system to get
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<On exit, a pointer to the coordinate system with the specified user number if it exists. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("CoordinateSystem_Get",err,error,*999)

    CALL CoordinateSystem_UserNumberFind(coordinateSystems,userNumber,coordinateSystem,err,error,*999)

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(coordinateSystem)) THEN
      localError="A coordinate system with an user number of "//TRIM(NumberToVString(userNumber,"*",err,error))//" does not exist."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
  
    EXITS("CoordinateSystem_Get")
    RETURN
999 ERRORSEXITS("CoordinateSystem_Get",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_Get

  !
  !================================================================================================================================
  !

  !>Returns the origin of a coordinate system. \see OpenCMISS::Iron::cmfe_CoordinateSystem_OriginGet
  SUBROUTINE CoordinateSystem_OriginGet(coordinateSystem,origin,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<A pointer to the coordinate system to get the origin for
    REAL(DP), INTENT(OUT) :: origin(:) !<On return, the origin of the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CoordinateSystem_OriginGet",err,error,*999)

    CALL CoordinateSystem_AssertIsFinished(coordinateSystem,err,error,*999)
#ifdef WITH_PRECHECKS    
    IF(SIZE(origin,1)<3) THEN
      localError="The size of the specified origin array is "// &
        & TRIM(NumberToVString(SIZE(origin,1),"*",err,error))//" is invalid. The size must be >= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    origin(1:3)=coordinateSystem%origin
     
    EXITS("CoordinateSystem_OriginGet")
    RETURN
999 ERRORSEXITS("CoordinateSystem_OriginGet",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_OriginGet

  !
  !================================================================================================================================
  !

  !>Returns the orientation of a coordinate system. \see OpenCMISS::Iron::cmfe_CoordinateSystem_OrientationGet
  SUBROUTINE CoordinateSystem_OrientationGet(coordinateSystem,orientation,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<A pointer to the coordinate system to get the orientation for
    REAL(DP), INTENT(OUT) :: orientation(:,:) !<On return, the orientation of the coordinate system
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CoordinateSystem_OrientationGet",err,error,*999)

    CALL CoordinateSystem_AssertIsFinished(coordinateSystem,err,error,*999)
#ifdef WITH_PRECHECKS    
    IF(SIZE(orientation,1)<3.OR.SIZE(orientation,2)<3) THEN
      localError="The size of the specified orientation array is "//TRIM(NumberToVString(SIZE(orientation,1),"*",err,error))// &
        & "x"//TRIM(NumberToVString(SIZE(orientation,2),"*",err,error))//" and it must be >= 3x3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    orientation(1:3,1:3)=coordinateSystem%orientation
   
    EXITS("CoordinateSystem_OrientationGet")
    RETURN
999 ERRORSEXITS("CoordinateSystem_OrientationGet",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_OrientationGet
  
  !
  !================================================================================================================================
  !

  !>Gets the coordinate system radial interpolation type. 
  SUBROUTINE CoordinateSystem_RadialInterpolationTypeGet(coordinateSystem,radialInterpolationType,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<The coordinate system to get the radial interpolation for
    INTEGER(INTG), INTENT(OUT) :: radialInterpolationType !<On return, the radial interpolation type for the coordinate system \see CoordinateRoutines_RadialInterpolations,CoordinateRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("CoordinateSystem_RadialInterpolationTypeGet",err,error,*999)

    CALL CoordinateSystem_AssertIsFinished(coordinateSystem,err,error,*999)
    
    SELECT CASE(coordinateSystem%type)
    CASE(COORDINATE_CYLINDRICAL_POLAR_TYPE,COORDINATE_SPHERICAL_POLAR_TYPE)
      radialInterpolationType=coordinateSystem%radialInterpolationType
    CASE DEFAULT
      localError="The radial interpolation type is not defined for coordinate system type "// &
        & TRIM(NumberToVString(coordinateSystem%TYPE,"*",err,error))
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("CoordinateSystem_RadialInterpolationTypeGet")
    RETURN
999 ERRORSEXITS("CoordinateSystem_RadialInterpolationTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_RadialInterpolationTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the coordinate system type. 
  SUBROUTINE CoordinateSystem_TypeGet(coordinateSystem,systemType,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<A pointer to the coordinate system to get the type for
    INTEGER(INTG), INTENT(OUT) :: systemType !<On return, the type for the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CoordinateSystem_TypeGet",err,error,*999)

    CALL CoordinateSystem_AssertIsFinished(coordinateSystem,err,error,*999)
    
    systemType=coordinateSystem%type
     
    EXITS("CoordinateSystem_TypeGet")
    RETURN
999 ERRORSEXITS("CoordinateSystem_TypeGet",err,error)
    RETURN 1

  END SUBROUTINE CoordinateSystem_TypeGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the coordinate system identified by a user number. If a coordinate system with that number is not
  !>found then coordinate system is set to NULL.
  SUBROUTINE CoordinateSystem_UserNumberFind(coordinateSystems,userNumber,coordinateSystem,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemsType), POINTER :: coordinateSystems !<The coordinate systems to find the user number for.
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the coordinate system to find.
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<On exit, a pointer to the coordinate system with the specified user number if it exists. If no coordinate system has the specified user number the pointer is returned as NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coordinateSystemIdx
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("CoordinateSystem_UserNumberFind",ERR,ERROR,*999)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(coordinateSystems)) CALL FlagError("Coordinate systems is not associated.",err,error,*999)
#endif    
   
    NULLIFY(coordinateSystem)
    IF(ASSOCIATED(coordinateSystems%coordinateSystems)) THEN
      DO coordinateSystemIdx=1,coordinateSystems%numberOfCoordinateSystems
#ifdef WITH_PRECHECKS        
        IF(.NOT.ASSOCIATED(coordinateSystems%coordinateSystems(coordinateSystemIdx)%ptr)) THEN
          localError="The coordinate system pointer in coordinate systems is not associated for coordinate system index "// &
            & TRIM(NumberToVString(coordinateSystemIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
#endif        
        IF(coordinateSystems%coordinateSystems(coordinateSystemIdx)%ptr%userNumber==userNumber) THEN
          coordinateSystem=>coordinateSystems%coordinateSystems(coordinateSystemIdx)%ptr
          EXIT
        ENDIF
      ENDDO !coordinateSystemIdx
    ENDIF
    
    EXITS("CoordinateSystem_UserNumberFind")
    RETURN
999 ERRORS("CoordinateSystem_UserNumberFind",err,error)
    EXITS("CoordinateSystem_UserNumberFind")
    RETURN 1
    
  END SUBROUTINE CoordinateSystem_UserNumberFind

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the world coordinate system for a coordinate systems.
  SUBROUTINE CoordinateSystems_WorldCoordinateSystemGet(coordinateSystems,worldCoordinateSystem,err,error,*)

    !Argument variables
    TYPE(CoordinateSystemsType), POINTER :: coordinateSystems !<A pointer to the coordinate systems to get the world coordinate system for
    TYPE(CoordinateSystemType), POINTER :: worldCoordinateSystem !<On exit, a pointer to the world coordinate system for the coordinate systems. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("CoordinateSystems_WorldCoordinateSystemGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(worldCoordinateSystem)) CALL FlagError("World coordinate system is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(coordinateSystems)) CALL FlagError("Coordinate systems is not associated.",err,error,*998)
#endif    

    CALL CoordinateSystem_UserNumberFind(coordinateSystems,0,worldCoordinateSystem,err,error,*999)

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(worldCoordinateSystem)) &
      & CALL FlagError("Coordinate systems world coordinate system is not associated.",err,error,*999)
#endif    
       
    EXITS("CoordinateSystems_WorldCoordinateSystemGet")
    RETURN
999 NULLIFY(worldCoordinateSystem)
998 ERRORSEXITS("CoordinateSystems_WorldCoordinateSystemGet",err,error)
    RETURN 1
    
  END SUBROUTINE CoordinateSystems_WorldCoordinateSystemGet
  
  !
  !================================================================================================================================
  !
  
END MODULE CoordinateSystemAccessRoutines
