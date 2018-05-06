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

  INTERFACE BASIS_FAMILY_NUMBER_FIND
    MODULE PROCEDURE Basis_FamilyNumberFind
  END INTERFACE BASIS_FAMILY_NUMBER_FIND

  INTERFACE BASIS_USER_NUMBER_FIND
    MODULE PROCEDURE Basis_UserNumberFind
  END INTERFACE BASIS_USER_NUMBER_FIND

  PUBLIC BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES,BASIS_DEFAULT_QUADRATURE_SCHEME,BASIS_LOW_QUADRATURE_SCHEME, &
    & BASIS_MID_QUADRATURE_SCHEME,BASIS_HIGH_QUADRATURE_SCHEME

  PUBLIC Basis_BasisFunctionsGet

  PUBLIC Basis_ContextGet
  
  PUBLIC Basis_FamilyNumberFind

  PUBLIC BASIS_FAMILY_NUMBER_FIND

  PUBLIC Basis_Get

  PUBLIC Basis_LocalFaceNumberGet

  PUBLIC Basis_LocalLineNumberGet

  PUBLIC Basis_QuadratureSchemeGet

  PUBLIC Basis_UserNumberFind

  PUBLIC BASIS_USER_NUMBER_FIND

  PUBLIC Basis_UserNumberGet

CONTAINS

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
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Basis_BasisFunctionsGet",err,error,*998)

    IF(ASSOCIATED(basisFunctions)) CALL FlagError("Basis functions is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)    
    IF(.NOT.basis%basisFinished) CALL FlagError("Basis has not been finished.",err,error,*999)
    
    basisFunctions=>basis%basisFunctions
    IF(.NOT.ASSOCIATED(basisFunctions)) THEN
      localError="Basis functions is not associated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    EXITS("Basis_BasisFunctionsGet")
    RETURN
999 NULLIFY(basisFunctions)
998 ERRORSEXITS("Basis_BasisFunctionsGet",err,error)    
    RETURN 1
    
  END SUBROUTINE Basis_BasisFunctionsGet

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
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Basis_ContextGet",err,error,*998)

    IF(ASSOCIATED(context)) CALL FlagError("Context is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)

    IF(.NOT.ASSOCIATED(basis%basisFunctions)) THEN
      localError="Basis functions is not associated for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    context=>basis%basisFunctions%context
    IF(.NOT.ASSOCIATED(context)) THEN
      localError="The context is not associated for the basis functions for basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Basis_ContextGet")
    RETURN
999 NULLIFY(context)
998 ERRORSEXITS("Basis_ContextGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_ContextGet

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
   TYPE(VARYING_STRING) :: localError

    ENTERS("Basis_FamilyNumberFind",err,error,*999)
    
    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(basisFunctions)) CALL FlagError("Basis functions is not associated.",err,error,*999)
    
    NULLIFY(basis)
    IF(ALLOCATED(basisFunctions%bases)) THEN
      BasisLoop: DO basisIdx=1,basisFunctions%numberOfBasisFunctions
        IF(ASSOCIATED(basisFunctions%bases(basisIdx)%ptr)) THEN
          IF(basisFunctions%bases(basisIdx)%ptr%userNumber==userNumber) THEN
            IF(familyNumber==0) THEN
              basis=>basisFunctions%bases(basisIdx)%ptr
              EXIT BasisLoop
            ELSE
!!TODO: \todo This only works for one level of sub-bases at the moment
              IF(ALLOCATED(basisFunctions%bases(basisIdx)%ptr%subBases)) THEN
                SubBasisLoop: DO subBasisIdx=1,basisFunctions%bases(basisIdx)%ptr%numberOfSubBases
                  IF(ASSOCIATED(basisFunctions%bases(basisIdx)%ptr%subBases(subBasisIdx)%ptr)) THEN
                    IF(basisFunctions%bases(basisIdx)%ptr%subBases(subBasisIdx)%ptr%familyNumber==familyNumber) THEN
                      basis=>basisFunctions%bases(basisIdx)%ptr%subBases(subBasisIdx)%ptr
                      EXIT BasisLoop
                    ENDIF
                  ELSE
                    localError="The sub basis is not associated for sub-basis index "// &
                      & TRIM(NumberToVString(subBasisIdx,"*",err,error))//" of basis index "// &
                      & TRIM(NumberToVString(basisIdx,"*",err,error))//"."
                    CALL FlagError(localError,err,error,*999)
                  ENDIF
                ENDDO SubBasisLoop !subBasisIdx
              ENDIF
            ENDIF
          ENDIF
        ELSE
          localError="The basis is not associated for basis index "//&
            & TRIM(NumberToVString(basisIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO BasisLoop !basisIdx
    ENDIF
        
    EXITS("Basis_FamilyNumberFind")
    RETURN
999 ERRORSEXITS("Basis_FamilyNumberFind",err,error)
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
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Basis_Get",err,error,*999)

    CALL Basis_UserNumberFind(basisFunctions,userNumber,basis,err,error,*999)
    IF(.NOT.ASSOCIATED(basis)) THEN
      localError="A basis with an user number of "//TRIM(NumberToVString(userNumber,"*",err,error))//" does not exist."
      CALL FlagError(localError,err,error,*999)
    ENDIF
  
    EXITS("Basis_Get")
    RETURN
999 ERRORSEXITS("Basis_Get",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_Get

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
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Basis_LocalFaceNumberGet",err,error,*999)

    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(.NOT.basis%basisFinished) CALL FlagError("Basis has not been finished.",err,error,*999)
    IF(normalXiDirection < -basis%numberOfXiCoordinates .OR. normalXiDirection > basis%numberOfXiCoordinates) THEN
      localError="The specified normal xi direction of "//TRIM(NumberToVString(normalXiDirection,"*",err,error))// &
        & " is invalid the normal xi direction must be >= "// &
        & TRIM(NumberToVString(-basis%numberOfXiCoordinates,"*",err,error))// &
        & " and <= "//TRIM(NumberToVString(basis%numberOfXiCoordinates,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    localFaceNumber=basis%xiNormalLocalFace(normalXiDirection)
    IF(localFaceNumber<1.OR.localFaceNumber>basis%numberOfLocalFaces) THEN
      localError="The local face number of "//TRIM(NumberToVString(localFaceNumber,"*",err,error))//" for xi normal direction "// &
        & TRIM(NumberToVString(normalXiDirection,"*",err,error))//" of basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//" is invalid. The local face number should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfLocalFaces,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
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
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Basis_LocalLineNumberGet",err,error,*999)

    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(.NOT.basis%basisFinished) CALL FlagError("Basis has not been finished.",err,error,*999)
    IF(.NOT.SIZE(normalXiDirections,1)==2) THEN
      localError="The specified number of normal xi directions of "// &
        & TRIM(NumberToVString(SIZE(normalXiDirections,1),"*",err,error))//" is invalid. There should be 2 normal xi directions."
      CALL FlagError(localError,err,error,*999)
    END IF
    IF(normalXiDirections(1) < -basis%numberOfXiCoordinates .OR. normalXiDirections(1) > basis%numberOfXiCoordinates) THEN
      localError="The first specified normal xi direction of "//TRIM(NumberToVString(normalXiDirections(1),"*",err,error))// &
        & " is invalid the normal xi direction must be >= "// &
        & TRIM(NumberToVString(-basis%numberOfXiCoordinates,"*",err,error))// &
        & " and <= "//TRIM(NumberToVString(basis%numberOfXiCoordinates,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(normalXiDirections(2) < -basis%numberOfXiCoordinates .OR. normalXiDirections(2) > basis%numberOfXiCoordinates) THEN
      localError="The second specified normal xi direction of "//TRIM(NumberToVString(normalXiDirections(2),"*",err,error))// &
        & " is invalid the normal xi direction must be >= "// &
        & TRIM(NumberToVString(-basis%numberOfXiCoordinates,"*",err,error))// &
        & " and <= "//TRIM(NumberToVString(basis%numberOfXiCoordinates,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    localLineNumber=basis%xiNormalsLocalLine(normalXiDirections(1),normalXiDirections(2))
    IF(localLineNumber<1.OR.localLineNumber>basis%numberOfLocalLines) THEN
      localError="The local line number of "//TRIM(NumberToVString(localLineNumber,"*",err,error))//" for xi normal directions "// &
        & TRIM(NumberToVString(normalXiDirections(1),"*",err,error))//","// &
        & TRIM(NumberToVString(normalXiDirections(2),"*",err,error))//" of basis number "// &
        & TRIM(NumberToVString(basis%userNumber,"*",err,error))//" is invalid. The local line number should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfLocalLines,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    EXITS("Basis_LocalLineNumberGet")
    RETURN
999 ERRORSEXITS("Basis_LocalLineNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_LocalLineNumberGet

  !
  !================================================================================================================================
  !

  !>Returns a quadrature scheme for the basis. 
  SUBROUTINE Basis_QuadratureSchemeGet(basis,quadratureSchemeIdx,quadratureScheme,err,error,*)

    !Argument variables
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to get the quadrature scheme for
    INTEGER(INTG), INTENT(IN) :: quadratureSchemeIdx !<The index of the quadrature scheme to get. \see Basis_QuadratureScheme,BasisRoutines,BasisAccessRoutines
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: quadratureScheme !<On exit, the basis quadrature scheme corresponding to the quadrature index. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Basis_QuadratureSchemeGet",err,error,*998)

    IF(ASSOCIATED(quadratureScheme)) CALL FlagError("Quadrature scheme is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(.NOT.basis%basisFinished) CALL FlagError("Basis has not been finished.",err,error,*999)
    IF(quadratureSchemeIdx<1.OR.quadratureSchemeIdx>BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES) THEN
      localError="The specified quadrature scheme index of "//TRIM(NumberToVString(quadratureSchemeIdx,"*",err,error))// &
        & " is invalid. The quadrature scheme index should be >= 1 and <= "// &
        & TRIM(NumberToVString(BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    END IF
    IF(.NOT.ALLOCATED(basis%quadrature%QUADRATURE_SCHEME_MAP)) &
      & CALL FlagError("Basis quadrature scheme map has not been allocated.",err,error,*999)

    quadratureScheme=>basis%quadrature%QUADRATURE_SCHEME_MAP(quadratureSchemeIdx)%ptr
    IF(.NOT.ASSOCIATED(quadratureScheme)) THEN
      localError="The quadrature scheme for the quadrature scheme index of "// &
        & TRIM(NumberToVString(quadratureSchemeIdx,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    EXITS("Basis_QuadratureSchemeGet")
    RETURN
999 NULLIFY(quadratureScheme)
998 ERRORSEXITS("Basis_QuadratureSchemeGet",err,error)    
    RETURN 1
    
  END SUBROUTINE Basis_QuadratureSchemeGet

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

    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)

    userNumber=basis%userNumber
  
    EXITS("Basis_UserNumberGet")
    RETURN
999 ERRORSEXITS("Basis_UserNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE Basis_UserNumberGet

  !
  !================================================================================================================================
  !

END MODULE BasisAccessRoutines
