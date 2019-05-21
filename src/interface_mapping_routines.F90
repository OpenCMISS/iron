!> \file
!> \author Chris Bradley
!> \brief This module contains all interface mapping routines.
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

!>This module contains all interface mapping routines.
MODULE INTERFACE_MAPPING_ROUTINES

  USE BaseRoutines
  USE FIELD_ROUTINES
  USE FieldAccessRoutines
  USE INPUT_OUTPUT
  USE INTERFACE_CONDITIONS_CONSTANTS
  USE InterfaceEquationsAccessRoutines
  USE InterfaceMappingAccessRoutines
  USE InterfaceMatricesAccessRoutines
  USE ISO_VARYING_STRING
  USE KINDS
  USE STRINGS
  USE TYPES

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  PUBLIC INTERFACE_MAPPING_CREATE_FINISH,INTERFACE_MAPPING_CREATE_START

  PUBLIC INTERFACE_MAPPING_DESTROY

  PUBLIC InterfaceMapping_LagrangeVariableSet

  PUBLIC INTERFACE_MAPPING_MATRICES_COEFFS_SET

  PUBLIC InterfaceMapping_MatricesColumnMeshIndicesSet,InterfaceMapping_MatricesRowMeshIndicesSet

  PUBLIC INTERFACE_MAPPING_MATRICES_NUMBER_SET

  PUBLIC INTERFACE_MAPPING_MATRICES_TRANSPOSE_SET

  PUBLIC INTERFACE_MAPPING_RHS_COEFF_SET

  PUBLIC INTERFACE_MAPPING_RHS_VARIABLE_TYPE_SET

CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates interface mapping
  SUBROUTINE INTERFACE_MAPPING_CALCULATE(INTERFACE_MAPPING,ERR,ERROR,*)
    
    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping to calculate the mapping for
    INTEGER(INTG), INTENT(OUT) ::       ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) ::    ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: column_idx,dof_idx,matrix_idx,mesh_idx,variable_idx,number_of_interface_matrices
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FieldType), POINTER :: LAGRANGE_FIELD
    TYPE(FieldVariableType), POINTER :: FIELD_VARIABLE,lagrangeVariable
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_DEPENDENT_TYPE), POINTER :: INTERFACE_DEPENDENT
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_LAGRANGE_TYPE), POINTER :: LAGRANGE
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(InterfaceMappingRHSType), POINTER :: rhsMapping
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("INTERFACE_MAPPING_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      createValuesCache=>INTERFACE_MAPPING%createValuesCache
      IF(ASSOCIATED(createValuesCache)) THEN
        INTERFACE_EQUATIONS=>INTERFACE_MAPPING%interfaceEquations
        IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
          INTERFACE_CONDITION=>INTERFACE_EQUATIONS%INTERFACE_CONDITION
          SELECT CASE(INTERFACE_CONDITION%METHOD)
          CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
            LAGRANGE=>INTERFACE_CONDITION%LAGRANGE
            IF(ASSOCIATED(LAGRANGE)) THEN
              INTERFACE_DEPENDENT=>INTERFACE_CONDITION%DEPENDENT
              IF(ASSOCIATED(INTERFACE_DEPENDENT)) THEN
                !Set the Lagrange variable information
                LAGRANGE_FIELD=>LAGRANGE%LAGRANGE_FIELD
                NULLIFY(lagrangeVariable)
                CALL Field_VariableGet(LAGRANGE_FIELD,createValuesCache%lagrangeVariableType,lagrangeVariable, &
                  & ERR,ERROR,*999)
                INTERFACE_MAPPING%lagrangeVariableType=createValuesCache%lagrangeVariableType
                INTERFACE_MAPPING%lagrangeVariable=>lagrangeVariable
                !Set the number of columns in the interface matrices
                INTERFACE_MAPPING%numberOfColumns=lagrangeVariable%numberOfDofs
                INTERFACE_MAPPING%totalNumberOfColumns=lagrangeVariable%totalNumberOfDofs
                INTERFACE_MAPPING%numberOfGlobalColumns=lagrangeVariable%numberOfGlobalDofs
                !Set the column dofs mapping
                INTERFACE_MAPPING%columnDOFSMapping=>lagrangeVariable%domainMapping
                ALLOCATE(INTERFACE_MAPPING%lagrangeDOFToColumnMap(lagrangeVariable%totalNumberOfDofs),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate Lagrange dof to column map.",ERR,ERROR,*999)
                !1-1 mapping for now
                DO dof_idx=1,lagrangeVariable%totalNumberOfDofs
                  column_idx=lagrangeVariable%domainMapping%localToGlobalMap(dof_idx)
                  INTERFACE_MAPPING%lagrangeDOFToColumnMap(dof_idx)=column_idx
                ENDDO
                !Set the number of interface matrices
                INTERFACE_MAPPING%numberOfInterfaceMatrices=createValuesCache%numberOfInterfaceMatrices
                ALLOCATE(INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(INTERFACE_MAPPING%numberOfInterfaceMatrices), &
                  & STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate interface matrix rows to variable maps.",ERR,ERROR,*999)
                !Loop over the interface matrices and calculate the row mappings
                !The pointers below have been checked for association above.
                SELECT CASE(INTERFACE_CONDITION%METHOD)
                CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                  number_of_interface_matrices=INTERFACE_MAPPING%numberOfInterfaceMatrices
                CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                  !Number of interface matrices whose rows/columns are related to Dependent/Lagrange variables and not Lagrange/Lagrange variables (last interface matrix is Lagrange/Lagrange (Penalty matrix)
                  number_of_interface_matrices=INTERFACE_MAPPING%numberOfInterfaceMatrices-1 
                ENDSELECT
                DO matrix_idx=1,number_of_interface_matrices
                  !Initialise and setup the interface matrix
                  CALL InterfaceMapping_MatrixToVarMapInitialise(INTERFACE_MAPPING,matrix_idx,ERR,ERROR,*999)
                  mesh_idx=createValuesCache%matrixRowFieldVariableIndices(matrix_idx)
                  NULLIFY(EQUATIONS_SET)
                  NULLIFY(FIELD_VARIABLE)
                  DO variable_idx=1,INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES
                    IF(INTERFACE_DEPENDENT%VARIABLE_MESH_INDICES(variable_idx)==mesh_idx) THEN
                      EQUATIONS_SET=>INTERFACE_DEPENDENT%EQUATIONS_SETS(variable_idx)%PTR
                      FIELD_VARIABLE=>INTERFACE_DEPENDENT%fieldVariables(variable_idx)%PTR
                      EXIT
                    ENDIF
                  ENDDO !variable_idx
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                      INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%equationsSet=>EQUATIONS_SET
                      INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%variableType=FIELD_VARIABLE%variableType
                      INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%variable=>FIELD_VARIABLE
                      INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%meshIndex=mesh_idx
                      INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%matrixCoefficient=INTERFACE_MAPPING% &
                        & createValuesCache%matrixCoefficients(matrix_idx)
                      INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%hasTranspose=INTERFACE_MAPPING% &
                        & createValuesCache%hasTranspose(matrix_idx)
                       !Set the number of rows
                      INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%numberOfRows=FIELD_VARIABLE%numberOfDofs
                      INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%totalNumberOfRows= &
                        & FIELD_VARIABLE%totalNumberOfDofs
                      INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%numberOfGlobalRows= &
                        & FIELD_VARIABLE%numberOfGlobalDofs
                      !Set the row mapping
                      INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%rowDOFsMapping=> &
                        & FIELD_VARIABLE%domainMapping
                      ALLOCATE(INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%variableDOFToRowMap( &
                        & FIELD_VARIABLE%totalNumberOfDofs),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate variable dof to row map.",ERR,ERROR,*999)
                      !1-1 mapping for now
                      DO dof_idx=1,FIELD_VARIABLE%totalNumberOfDofs
                        INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%variableDOFToRowMap(dof_idx)=dof_idx
                      ENDDO !dof_idx
                    ELSE
                      LOCAL_ERROR="Dependent variable for mesh index "//TRIM(NUMBER_TO_VSTRING(mesh_idx,"*",ERR,ERROR))// &
                        & " could not be found."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="Equations set for mesh index "//TRIM(NUMBER_TO_VSTRING(mesh_idx,"*",ERR,ERROR))// &
                      & " could not be found."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDDO !matrix_idx

                !The pointers below have been checked for association above.
                SELECT CASE(INTERFACE_CONDITION%METHOD)
                CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                  !Sets up the Lagrange-(Penalty) interface matrix mapping and calculate the row mappings
                  matrix_idx = INTERFACE_MAPPING%numberOfInterfaceMatrices !last of the interface matrices
                  !Initialise and setup the interface matrix
                  CALL InterfaceMapping_MatrixToVarMapInitialise(INTERFACE_MAPPING,matrix_idx,ERR,ERROR,*999)
                  mesh_idx=createValuesCache%matrixRowFieldVariableIndices(matrix_idx)
                  NULLIFY(lagrangeVariable)
                  CALL Field_VariableGet(LAGRANGE_FIELD,createValuesCache%lagrangeVariableType,lagrangeVariable, &
                    & ERR,ERROR,*999)
                  NULLIFY(INTERFACE_EQUATIONS)
                  NULLIFY(FIELD_VARIABLE)
                  FIELD_VARIABLE=>lagrangeVariable
                  INTERFACE_EQUATIONS=>INTERFACE_CONDITION%interfaceEquations
                  IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
                    IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                      INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%interfaceEquations=>INTERFACE_EQUATIONS
                      INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%variableType=FIELD_VARIABLE%variableType
                      INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%VARIABLE=>FIELD_VARIABLE
                      INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%meshIndex=mesh_idx
                      INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%matrixCoefficient=INTERFACE_MAPPING% &
                        & createValuesCache%matrixCoefficients(matrix_idx)
                      INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%hasTranspose=INTERFACE_MAPPING% &
                        & createValuesCache%hasTranspose(matrix_idx)
                        !Set the number of rows
                      INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%numberOfRows=FIELD_VARIABLE%numberOfDofs
                      INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%totalNumberOfRows= &
                        & FIELD_VARIABLE%totalNumberOfDofs
                      INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%numberOfGlobalRows= &
                        & FIELD_VARIABLE%numberOfGlobalDofs
                      !Set the row mapping
                      INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%rowDOFsMapping=> &
                        & FIELD_VARIABLE%domainMapping
                      ALLOCATE(INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%variableDOFToRowMap( &
                        & FIELD_VARIABLE%totalNumberOfDofs),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate variable dof to row map.",ERR,ERROR,*999)
                      !1-1 mapping for now
                      DO dof_idx=1,FIELD_VARIABLE%totalNumberOfDofs
                        INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%variableDOFToRowMap(dof_idx)=dof_idx
                      ENDDO !dof_idx
                    ELSE
                      LOCAL_ERROR="Lagrange variable for mesh index "//TRIM(NUMBER_TO_VSTRING(mesh_idx,"*",ERR,ERROR))// &
                        & " could not be found."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="Interface Equations for mesh index "//TRIM(NUMBER_TO_VSTRING(mesh_idx,"*",ERR,ERROR))// &
                      & " could not be found."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDSELECT

                !Calculate RHS mappings
                IF(createValuesCache%rhsLagrangeVariableType/=0) THEN
                  CALL INTERFACE_MAPPING_RHS_MAPPING_INITIALISE(INTERFACE_MAPPING,ERR,ERROR,*999)
                  rhsMapping=>INTERFACE_MAPPING%rhsMapping
                  IF(ASSOCIATED(rhsMapping)) THEN
                    rhsMapping%rhsVariableType=createValuesCache%rhsLagrangeVariableType
                    lagrangeVariable=>LAGRANGE_FIELD%variableTypeMap(createValuesCache%rhsLagrangeVariableType)%PTR
                    rhsMapping%rhsVariable=>lagrangeVariable
                    rhsMapping%rhsVariableMapping=>lagrangeVariable%domainMapping
                    rhsMapping%rhsCoefficient=createValuesCache%rhsCoefficient
                    !Allocate and set up the row mappings
                    ALLOCATE(rhsMapping%rhsDOFToInterfaceRowMap(lagrangeVariable%totalNumberOfDofs),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate rhs dof to interface row map.",ERR,ERROR,*999)
                    ALLOCATE(rhsMapping%interfaceRowToRHSDOFMap(INTERFACE_MAPPING%totalNumberOfColumns),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate interface row to dof map.",ERR,ERROR,*999)
                    DO dof_idx=1,lagrangeVariable%totalNumberOfDofs
                      !1-1 mapping for now
                      column_idx=dof_idx
                      rhsMapping%rhsDOFToInterfaceRowMap(dof_idx)=column_idx
                    ENDDO !dof_idx
                    DO column_idx=1,INTERFACE_MAPPING%totalNumberOfColumns
                      !1-1 mapping for now
                      dof_idx=column_idx
                      rhsMapping%interfaceRowToRHSDOFMap(column_idx)=dof_idx
                    ENDDO !column_idx
                  ELSE
                    CALL FlagError("RHS mapping is not associated.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ELSE
                CALL FlagError("Interface condition dependent is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Interface condition Lagrange is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The interface condition method of "// &
              & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Interface equations interface condition is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Interface mapping create values cache is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface mapping is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("INTERFACE_MAPPING_CALCULATE")
    RETURN
999 ERRORSEXITS("INTERFACE_MAPPING_CALCULATE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE INTERFACE_MAPPING_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finishes the creation of interface mapping
  SUBROUTINE INTERFACE_MAPPING_CREATE_FINISH(INTERFACE_MAPPING,ERR,ERROR,*)
    
    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping to finish the creation of.
    INTEGER(INTG), INTENT(OUT) ::       ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) ::    ERROR !<The error string
    !Local Variables

    ENTERS("INTERFACE_MAPPING_CREATE_FINISH",ERR,ERROR,*999)

    CALL InterfaceMapping_AssertNotFinished(INTERFACE_MAPPING,ERR,ERROR,*999)

    !Calculate the equations mapping and clean up
    CALL INTERFACE_MAPPING_CALCULATE(INTERFACE_MAPPING,ERR,ERROR,*999)
    CALL InterfaceMapping_CreateValuesCacheFinalise(INTERFACE_MAPPING%createValuesCache,ERR,ERROR,*999)
    INTERFACE_MAPPING%interfaceMappingFinished=.TRUE.
    
    EXITS("INTERFACE_MAPPING_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("INTERFACE_MAPPING_CREATE_FINISH",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE INTERFACE_MAPPING_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the process of creating the interface mapping for interface equations.
  SUBROUTINE INTERFACE_MAPPING_CREATE_START(INTERFACE_EQUATIONS,INTERFACE_MAPPING,ERR,ERROR,*)
    
    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to create the mapping for.
    TYPE(InterfaceMappingType), POINTER :: INTERFACE_MAPPING !<On exit, a pointer to the created interface mapping. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("INTERFACE_MAPPING_CREATE_START",ERR,ERROR,*999)

    CALL InterfaceEquations_AssertIsFinished(INTERFACE_EQUATIONS,ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      CALL FlagError("Interface mapping is already associated.",ERR,ERROR,*999)
    ELSE
      NULLIFY(INTERFACE_MAPPING)
      CALL INTERFACE_MAPPING_INITIALISE(INTERFACE_EQUATIONS,ERR,ERROR,*999)
      INTERFACE_MAPPING=>INTERFACE_EQUATIONS%interfaceMapping
    ENDIF
   
    EXITS("INTERFACE_MAPPING_CREATE_START")
    RETURN
999 ERRORSEXITS("INTERFACE_MAPPING_CREATE_START",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE INTERFACE_MAPPING_CREATE_START

  !
  !================================================================================================================================
  !

  !>Finalises an interface mapping create values cache and deallocates all memory
  SUBROUTINE InterfaceMapping_CreateValuesCacheFinalise(createValuesCache,ERR,ERROR,*)

    !Argument variables
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the create values cache to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("InterfaceMapping_CreateValuesCacheFinalise",ERR,ERROR,*999)

    IF(ASSOCIATED(createValuesCache)) THEN
      IF(ALLOCATED(createValuesCache%matrixCoefficients)) DEALLOCATE(createValuesCache%matrixCoefficients)
      IF(ALLOCATED(createValuesCache%hasTranspose)) DEALLOCATE(createValuesCache%hasTranspose)
      IF(ALLOCATED(createValuesCache%matrixRowFieldVariableIndices))  &
        & DEALLOCATE(createValuesCache%matrixRowFieldVariableIndices)
      IF(ALLOCATED(createValuesCache%matrixColFieldVariableIndices)) &
        & DEALLOCATE(createValuesCache%matrixColFieldVariableIndices)
      DEALLOCATE(createValuesCache)
    ENDIF
       
    EXITS("InterfaceMapping_CreateValuesCacheFinalise")
    RETURN
999 ERRORSEXITS("InterfaceMapping_CreateValuesCacheFinalise",ERR,ERROR)
    RETURN 1
  END SUBROUTINE InterfaceMapping_CreateValuesCacheFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an interface mapping create values cache 
  SUBROUTINE InterfaceMapping_CreateValuesCacheInitialise(INTERFACE_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping to create the create values cache for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,variable_idx,variable_type_idx,variable_type_idx2
    TYPE(FieldType), POINTER :: LAGRANGE_FIELD
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_DEPENDENT_TYPE), POINTER :: INTERFACE_DEPENDENT
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_LAGRANGE_TYPE), POINTER :: LAGRANGE
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    ENTERS("InterfaceMapping_CreateValuesCacheInitialise",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      IF(ASSOCIATED(INTERFACE_MAPPING%createValuesCache)) THEN
        CALL FlagError("Interface mapping create values cache is already associated.",ERR,ERROR,*998)
      ELSE
        INTERFACE_EQUATIONS=>INTERFACE_MAPPING%interfaceEquations
        IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
          INTERFACE_CONDITION=>INTERFACE_EQUATIONS%INTERFACE_CONDITION
          IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
            !Allocate and initialise the create values cache
            ALLOCATE(INTERFACE_MAPPING%createValuesCache,STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate interface mapping create values cache.",ERR,ERROR,*999)
            INTERFACE_MAPPING%createValuesCache%numberOfInterfaceMatrices=0
            INTERFACE_MAPPING%createValuesCache%lagrangeVariableType=0
            INTERFACE_MAPPING%createValuesCache%rhsLagrangeVariableType=0
            INTERFACE_MAPPING%createValuesCache%rhsCoefficient=0.0_DP
            !Set the default interface mapping in the create values cache
            !First calculate how many interface matrices we have and set the variable types
            SELECT CASE(INTERFACE_CONDITION%METHOD)
            CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
              LAGRANGE=>INTERFACE_CONDITION%LAGRANGE
              IF(ASSOCIATED(LAGRANGE)) THEN
                LAGRANGE_FIELD=>LAGRANGE%LAGRANGE_FIELD
                IF(ASSOCIATED(LAGRANGE_FIELD)) THEN
                  INTERFACE_DEPENDENT=>INTERFACE_CONDITION%DEPENDENT
                  IF(ASSOCIATED(INTERFACE_DEPENDENT)) THEN
                    !The pointers below have been checked for association above.
                    SELECT CASE(INTERFACE_CONDITION%METHOD)
                    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                      !Default the number of interface matrices to the number of added dependent variables
                      INTERFACE_MAPPING%createValuesCache%numberOfInterfaceMatrices= &
                      INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES
                    CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                      !Default the number of interface matrices to the number of added dependent variables plus a single Lagrange variable
                      INTERFACE_MAPPING%createValuesCache%numberOfInterfaceMatrices= &
                      INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES+1
                    END SELECT
                    !Default the Lagrange variable to the first Lagrange variable
                    INTERFACE_MAPPING%createValuesCache%lagrangeVariableType=0
                    DO variable_type_idx=1,FIELD_NUMBER_OF_VARIABLE_TYPES
                      IF(ASSOCIATED(LAGRANGE_FIELD%variableTypeMap(variable_type_idx)%PTR)) THEN
                        INTERFACE_MAPPING%createValuesCache%lagrangeVariableType=variable_type_idx
                        EXIT
                      ENDIF
                    ENDDO !variable_type_idx
                    IF(INTERFACE_MAPPING%createValuesCache%lagrangeVariableType==0) &
                      & CALL FlagError("Could not find a Lagrange variable type in the Lagrange field.",ERR,ERROR,*999)
                    !Default the RHS Lagrange variable to the second Lagrange variable
                    DO variable_type_idx2=variable_type_idx+1,FIELD_NUMBER_OF_VARIABLE_TYPES
                      IF(ASSOCIATED(LAGRANGE_FIELD%variableTypeMap(variable_type_idx2)%PTR)) THEN
                        INTERFACE_MAPPING%createValuesCache%rhsLagrangeVariableType=variable_type_idx2
                        EXIT
                      ENDIF
                    ENDDO !variable_type_idx2
                    IF(INTERFACE_MAPPING%createValuesCache%rhsLagrangeVariableType==0) &
                      & CALL FlagError("Could not find a RHS Lagrange variable type in the Lagrange field.",ERR,ERROR,*999)
                    ALLOCATE(INTERFACE_MAPPING%createValuesCache%matrixCoefficients(INTERFACE_MAPPING% &
                      & createValuesCache%numberOfInterfaceMatrices),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate create values cache matrix coefficients.",ERR,ERROR,*999)
                    !Default the interface matrices coefficients to add.
                    INTERFACE_MAPPING%createValuesCache%matrixCoefficients=1.0_DP
                    INTERFACE_MAPPING%createValuesCache%rhsCoefficient=1.0_DP
                    ALLOCATE(INTERFACE_MAPPING%createValuesCache%hasTranspose(INTERFACE_MAPPING% &
                      & createValuesCache%numberOfInterfaceMatrices),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate create values cache has transpose.",ERR,ERROR,*999)
                    !Default the interface matrices to all have a transpose
                    INTERFACE_MAPPING%createValuesCache%hasTranspose=.TRUE.
                    ALLOCATE(INTERFACE_MAPPING%createValuesCache%matrixRowFieldVariableIndices(INTERFACE_MAPPING% &
                      & createValuesCache%numberOfInterfaceMatrices),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate create values cache matrix row field variable indexes.", &
                      & ERR,ERROR,*999)
                    !Default the interface matrices to be in mesh index order.
                    DO variable_idx=1,INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES
                      INTERFACE_MAPPING%createValuesCache%matrixRowFieldVariableIndices(variable_idx)=variable_idx
                    ENDDO !variable_idx
                    !The pointers below have been checked for association above.
                    SELECT CASE(INTERFACE_CONDITION%METHOD)
                    CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                      !Default the interface matrix (Penalty) to have no transpose
                      INTERFACE_MAPPING%createValuesCache%hasTranspose(INTERFACE_MAPPING% &
                        & createValuesCache%numberOfInterfaceMatrices)=.FALSE.
                      !Default the interface matrices to be in mesh index order (and set Penalty matrix (last interface matrix)to be the first Lagrange variable).
                      INTERFACE_MAPPING%createValuesCache%matrixRowFieldVariableIndices(INTERFACE_DEPENDENT% &
                        & NUMBER_OF_DEPENDENT_VARIABLES+1)=1
                    END SELECT
                  ELSE
                    CALL FlagError("Interface condition depdendent is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Interface condition Lagrange field is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Interface condition Lagrange is not associated.",ERR,ERROR,*999)
              ENDIF
            CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The interface equations method of "// &
                & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FlagError("Interface equations interface condition is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Interface mapping interface equations is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface mapping is not associated.",ERR,ERROR,*998)
    ENDIF
       
    EXITS("InterfaceMapping_CreateValuesCacheInitialise")
    RETURN
999 CALL InterfaceMapping_CreateValuesCacheFinalise(INTERFACE_MAPPING%createValuesCache,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORS("InterfaceMapping_CreateValuesCacheInitialise",ERR,ERROR)
    EXITS("InterfaceMapping_CreateValuesCacheInitialise")
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_CreateValuesCacheInitialise

  !
  !================================================================================================================================
  !

  !>Destroys an interface mapping.
  SUBROUTINE INTERFACE_MAPPING_DESTROY(INTERFACE_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: INTERFACE_MAPPING !<A pointer the interface mapping to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("INTERFACE_MAPPING_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      CALL INTERFACE_MAPPING_FINALISE(INTERFACE_MAPPING,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Equations mapping is not associated.",ERR,ERROR,*999)
    ENDIF
        
    EXITS("INTERFACE_MAPPING_DESTROY")
    RETURN
999 ERRORSEXITS("INTERFACE_MAPPING_DESTROY",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE INTERFACE_MAPPING_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises the interface mapping and deallocates all memory.
  SUBROUTINE INTERFACE_MAPPING_FINALISE(INTERFACE_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx

    ENTERS("INTERFACE_MAPPING_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      IF(ALLOCATED(INTERFACE_MAPPING%lagrangeDOFToColumnMap)) DEALLOCATE(INTERFACE_MAPPING%lagrangeDOFToColumnMap)
      IF(ALLOCATED(INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps)) THEN
        DO matrix_idx=1,SIZE(INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps,1)
          CALL InterfaceMapping_MatrixToVarMapFinalise(INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx), &
            & ERR,ERROR,*999)
        ENDDO !matrix_idx
        DEALLOCATE(INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps)
      ENDIF
      CALL INTERFACE_MAPPING_RHS_MAPPING_FINALISE(INTERFACE_MAPPING%rhsMapping,ERR,ERROR,*999)
      CALL InterfaceMapping_CreateValuesCacheFinalise(INTERFACE_MAPPING%createValuesCache,ERR,ERROR,*999)
      DEALLOCATE(INTERFACE_MAPPING)
    ENDIF
       
    EXITS("INTERFACE_MAPPING_FINALISE")
    RETURN
999 ERRORSEXITS("INTERFACE_MAPPING_FINALISE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE INTERFACE_MAPPING_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the interface mapping and deallocates all memory.
  SUBROUTINE INTERFACE_MAPPING_INITIALISE(INTERFACE_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to initialise the interface mapping for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS("INTERFACE_MAPPING_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(ASSOCIATED(INTERFACE_EQUATIONS%interfaceMapping)) THEN
        CALL FlagError("Interface mapping is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(INTERFACE_EQUATIONS%interfaceMapping,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate interface equations interface mapping.",ERR,ERROR,*999)
        INTERFACE_EQUATIONS%interfaceMapping%interfaceEquations=>INTERFACE_EQUATIONS
        INTERFACE_EQUATIONS%interfaceMapping%interfaceMappingFinished=.FALSE.
        INTERFACE_EQUATIONS%interfaceMapping%lagrangeVariableType=0
        NULLIFY(INTERFACE_EQUATIONS%interfaceMapping%lagrangeVariable)
        INTERFACE_EQUATIONS%interfaceMapping%numberOfColumns=0
        INTERFACE_EQUATIONS%interfaceMapping%totalNumberOfColumns=0
        INTERFACE_EQUATIONS%interfaceMapping%numberOfGlobalColumns=0
        NULLIFY(INTERFACE_EQUATIONS%interfaceMapping%columnDOFSMapping)
        INTERFACE_EQUATIONS%interfaceMapping%numberOfInterfaceMatrices=0
        NULLIFY(INTERFACE_EQUATIONS%interfaceMapping%rhsMapping)
        NULLIFY(INTERFACE_EQUATIONS%interfaceMapping%createValuesCache)
        CALL InterfaceMapping_CreateValuesCacheInitialise(INTERFACE_EQUATIONS%interfaceMapping,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface equations is not associated.",ERR,ERROR,*998)
    ENDIF
       
    EXITS("INTERFACE_MAPPING_INITIALISE")
    RETURN
999 CALL INTERFACE_MAPPING_FINALISE(INTERFACE_EQUATIONS%interfaceMapping,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("INTERFACE_MAPPING_INITIALISE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE INTERFACE_MAPPING_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the Lagrange variable type for the interface mapping. 
  SUBROUTINE InterfaceMapping_LagrangeVariableSet(INTERFACE_MAPPING,lagrangeVariableType,ERR,ERROR,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping.
    INTEGER(INTG), INTENT(IN) :: lagrangeVariableType !<The Lagrange variable type to set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FieldType), POINTER :: LAGRANGE_FIELD
    TYPE(FieldVariableType), POINTER :: lagrangeVariable
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_LAGRANGE_TYPE), POINTER :: LAGRANGE
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("InterfaceMapping_LagrangeVariableSet",ERR,ERROR,*999)

    CALL InterfaceMapping_AssertNotFinished(INTERFACE_MAPPING,ERR,ERROR,*999)
    createValuesCache=>INTERFACE_MAPPING%createValuesCache
    IF(ASSOCIATED(createValuesCache)) THEN
      INTERFACE_EQUATIONS=>INTERFACE_MAPPING%interfaceEquations
      IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
        INTERFACE_CONDITION=>INTERFACE_EQUATIONS%INTERFACE_CONDITION
        IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
          SELECT CASE(INTERFACE_CONDITION%METHOD)
          CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
            LAGRANGE=>INTERFACE_CONDITION%LAGRANGE
            IF(ASSOCIATED(LAGRANGE)) THEN
              IF(LAGRANGE%LAGRANGE_FINISHED) THEN
                LAGRANGE_FIELD=>LAGRANGE%LAGRANGE_FIELD
                NULLIFY(lagrangeVariable)
                CALL Field_VariableGet(LAGRANGE_FIELD,lagrangeVariableType,lagrangeVariable,ERR,ERROR,*999)
                createValuesCache%lagrangeVariableType=lagrangeVariableType
              ELSE
                CALL FlagError("Interface condition Lagrange field has not been finished.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Interface condition Lagrange is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The interface condition method of "// &
              & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Interface equations interface condition is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Interface mapping interface equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface mapping create values cache is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("InterfaceMapping_LagrangeVariableSet")
    RETURN
999 ERRORSEXITS("InterfaceMapping_LagrangeVariableSet",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_LagrangeVariableSet

  !
  !================================================================================================================================
  !

  !>Finalises an interface matrix to variable map and deallocates all memory.
  SUBROUTINE InterfaceMapping_MatrixToVarMapFinalise(INTERFACE_MATRIX_TO_VAR_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(InterfaceMatrixToVarMapType) :: INTERFACE_MATRIX_TO_VAR_MAP !<The interface matrix to var map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
  
    ENTERS("InterfaceMapping_MatrixToVarMapFinalise",ERR,ERROR,*999)

    IF(ALLOCATED(INTERFACE_MATRIX_TO_VAR_MAP%variableDOFToRowMap)) &
      & DEALLOCATE(INTERFACE_MATRIX_TO_VAR_MAP%variableDOFToRowMap)
    
    EXITS("InterfaceMapping_MatrixToVarMapFinalise")
    RETURN
999 ERRORSEXITS("InterfaceMapping_MatrixToVarMapFinalise",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_MatrixToVarMapFinalise

  !
  !================================================================================================================================
  !

  !>Initialises an interface matrix to variable map.
  SUBROUTINE InterfaceMapping_MatrixToVarMapInitialise(INTERFACE_MAPPING,matrix_idx,ERR,ERROR,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping to initialise the matrix to variable map for a given matrix index.
    INTEGER(INTG), INTENT(IN) :: matrix_idx !<The matrix index to intialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
  
    ENTERS("InterfaceMapping_MatrixToVarMapInitialise",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      IF(ALLOCATED(INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps)) THEN
        IF(matrix_idx>0.AND.matrix_idx<=INTERFACE_MAPPING%numberOfInterfaceMatrices) THEN
          INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%matrixNumber=matrix_idx
          NULLIFY(INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%interfaceMatrix)
          NULLIFY(INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%equationsSet)
          INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%variableType=0
          NULLIFY(INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%variable)
          INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%meshIndex=0
          INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%matrixCoefficient=0.0_DP
          INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%hasTranspose=.FALSE.
          INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%numberOfRows=0
          INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%totalNumberOfRows=0
          INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%numberOfGlobalRows=0
          NULLIFY(INTERFACE_MAPPING%interfaceMatrixRowsToVarMaps(matrix_idx)%rowDOFsMapping)          
        ELSE
          LOCAL_ERROR="The specified matrix index of "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
            & " is invalid. The index must be > 0 and <= "// &
            & TRIM(NUMBER_TO_VSTRING(INTERFACE_MAPPING%numberOfInterfaceMatrices,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Interface mapping matrix rows to var maps is not allocated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface mapping is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("InterfaceMapping_MatrixToVarMapInitialise")
    RETURN
999 ERRORSEXITS("InterfaceMapping_MatrixToVarMapInitialise",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_MatrixToVarMapInitialise

  !
  !================================================================================================================================
  !

  !>Sets the coefficients for the interface matrices. 
  SUBROUTINE INTERFACE_MAPPING_MATRICES_COEFFS_SET(INTERFACE_MAPPING,matrixCoefficients,ERR,ERROR,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping.
    REAL(DP), INTENT(IN) :: matrixCoefficients(:) !<The interface matrix coefficients
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("INTERFACE_MAPPING_MATRICES_COEFFS_SET",ERR,ERROR,*999)

    CALL InterfaceMapping_AssertNotFinished(INTERFACE_MAPPING,ERR,ERROR,*999)
    
    createValuesCache=>INTERFACE_MAPPING%createValuesCache
    IF(ASSOCIATED(createValuesCache)) THEN          
      INTERFACE_EQUATIONS=>INTERFACE_MAPPING%interfaceEquations
      IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
        INTERFACE_CONDITION=>INTERFACE_EQUATIONS%INTERFACE_CONDITION
        IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
          SELECT CASE(INTERFACE_CONDITION%METHOD)
          CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
            !Check that the number of supplied coefficients matches the number of interface matrices
            IF(SIZE(matrixCoefficients,1)==createValuesCache%numberOfInterfaceMatrices) THEN
              createValuesCache%matrixCoefficients(1:createValuesCache%numberOfInterfaceMatrices)= &
                & matrixCoefficients(1:createValuesCache%numberOfInterfaceMatrices)
            ELSE
              LOCAL_ERROR="Invalid size of matrix coefficeints. The size of the supplied array ("// &
                & TRIM(NUMBER_TO_VSTRING(SIZE(matrixCoefficients,1),"*",ERR,ERROR))// &
                & ") must match the number of interface matrices ("// &
                & TRIM(NUMBER_TO_VSTRING(createValuesCache%numberOfInterfaceMatrices,"*",ERR,ERROR))//")."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The interface condition method of "// &
              & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Interface equations interface condition is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Interface mapping interface equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface mapping create values cache is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("INTERFACE_MAPPING_MATRICES_COEFFS_SET")
    RETURN
999 ERRORSEXITS("INTERFACE_MAPPING_MATRICES_COEFFS_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE INTERFACE_MAPPING_MATRICES_COEFFS_SET

  !
  !================================================================================================================================
  !

  !>Sets the column mesh indices for the interface matrices. 
  SUBROUTINE InterfaceMapping_MatricesColumnMeshIndicesSet(INTERFACE_MAPPING,COLUMN_MESH_INDICES,ERR,ERROR,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping.
    INTEGER(INTG), INTENT(IN) :: COLUMN_MESH_INDICES(:) !<The interface matrix column mesh indices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("InterfaceMapping_MatricesColumnMeshIndicesSet",ERR,ERROR,*999)

    CALL InterfaceMapping_AssertNotFinished(INTERFACE_MAPPING,ERR,ERROR,*999)
    
    createValuesCache=>INTERFACE_MAPPING%createValuesCache
    IF(ASSOCIATED(createValuesCache)) THEN
      INTERFACE_EQUATIONS=>INTERFACE_MAPPING%interfaceEquations
      IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
        INTERFACE_CONDITION=>INTERFACE_EQUATIONS%INTERFACE_CONDITION
        IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
          SELECT CASE(INTERFACE_CONDITION%METHOD)
          CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
            CALL FlagError("Can not set the column mesh indices when using the Lagrange multipliers "// &
              "interface condition method.",ERR,ERROR,*999)
          CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(INTERFACE_CONDITION_PENALTY_METHOD)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The interface condition method of "// &
              & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Interface equations interface condition is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Interface mapping interface equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface mapping create values cache is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("InterfaceMapping_MatricesColumnMeshIndicesSet")
    RETURN
999 ERRORS("InterfaceMapping_MatricesColumnMeshIndicesSet",ERR,ERROR)
    EXITS("InterfaceMapping_MatricesColumnMeshIndicesSet")
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_MatricesColumnMeshIndicesSet

  !
  !================================================================================================================================
  !

  !>Sets the row mesh indices for the interface matrices. 
  SUBROUTINE InterfaceMapping_MatricesRowMeshIndicesSet(INTERFACE_MAPPING,ROW_MESH_INDICES,ERR,ERROR,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping.
    INTEGER(INTG), INTENT(IN) :: ROW_MESH_INDICES(:) !<The interface matrix mesh indices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: mesh_idx,mesh_idx2,mesh_idx3
    LOGICAL :: FOUND
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_DEPENDENT_TYPE), POINTER :: INTERFACE_DEPENDENT
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("InterfaceMapping_MatricesRowMeshIndicesSet",ERR,ERROR,*999)

    CALL InterfaceMapping_AssertNotFinished(INTERFACE_MAPPING,ERR,ERROR,*999)
    
    createValuesCache=>INTERFACE_MAPPING%createValuesCache
    IF(ASSOCIATED(createValuesCache)) THEN
      INTERFACE_EQUATIONS=>INTERFACE_MAPPING%interfaceEquations
      IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
        INTERFACE_CONDITION=>INTERFACE_EQUATIONS%INTERFACE_CONDITION
        IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
          SELECT CASE(INTERFACE_CONDITION%METHOD)
          CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
            !Check the size of the mesh indicies array
            IF(SIZE(ROW_MESH_INDICES,1)==createValuesCache%numberOfInterfaceMatrices) THEN
              !Check that mesh indices are valid.
              INTERFACE_DEPENDENT=>INTERFACE_CONDITION%DEPENDENT
              IF(ASSOCIATED(INTERFACE_DEPENDENT)) THEN
                DO mesh_idx=1,createValuesCache%numberOfInterfaceMatrices
                  FOUND=.FALSE.
                  DO mesh_idx2=1,INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES
                    IF(ROW_MESH_INDICES(mesh_idx)==INTERFACE_DEPENDENT%VARIABLE_MESH_INDICES(mesh_idx2)) THEN
                      FOUND=.TRUE.
                      EXIT
                    ENDIF
                  ENDDO !mesh_idx2
                  IF(FOUND) THEN
                    !Check that the mesh index has not been repeated.
                    DO mesh_idx3=mesh_idx+1,createValuesCache%numberOfInterfaceMatrices
                      IF(ROW_MESH_INDICES(mesh_idx)==ROW_MESH_INDICES(mesh_idx3)) THEN
                        LOCAL_ERROR="The supplied mesh index of "// &
                          & TRIM(NUMBER_TO_VSTRING(ROW_MESH_INDICES(mesh_idx),"*",ERR,ERROR))// &
                          & " at position "//TRIM(NUMBER_TO_VSTRING(mesh_idx,"*",ERR,ERROR))// &
                          & " has been repeated at position "//TRIM(NUMBER_TO_VSTRING(mesh_idx3,"*",ERR,ERROR))//"."
                        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ENDDO !mesh_idx3
                    !Set the mesh indices
                    createValuesCache%matrixRowFieldVariableIndices(1:createValuesCache%numberOfInterfaceMatrices)= &
                      & ROW_MESH_INDICES(1:createValuesCache%numberOfInterfaceMatrices)
                  ELSE
                    LOCAL_ERROR="The supplied mesh index of "// &
                      & TRIM(NUMBER_TO_VSTRING(ROW_MESH_INDICES(mesh_idx),"*",ERR,ERROR))// &
                      & " at position "//TRIM(NUMBER_TO_VSTRING(mesh_idx,"*",ERR,ERROR))// &
                      & " has not been added as a dependent variable to the interface condition."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDDO !mesh_idx
              ELSE
                CALL FlagError("Interface condition dependent is not assocaited.",ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="Invalid size of mesh indices. The size of the supplied array ("// &
                & TRIM(NUMBER_TO_VSTRING(SIZE(ROW_MESH_INDICES,1),"*",ERR,ERROR))// &
                & ") must match the number of interface matrices ("// &
                & TRIM(NUMBER_TO_VSTRING(createValuesCache%numberOfInterfaceMatrices,"*",ERR,ERROR))//")."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The interface condition method of "// &
              & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Interface equations interface condition is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Interface mapping interface equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface mapping create values cache is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("InterfaceMapping_MatricesRowMeshIndicesSet")
    RETURN
999 ERRORSEXITS("InterfaceMapping_MatricesRowMeshIndicesSet",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_MatricesRowMeshIndicesSet

  !
  !================================================================================================================================
  !

  !>Sets the number of interface matrices for an interface mapping.
  SUBROUTINE INTERFACE_MAPPING_MATRICES_NUMBER_SET(INTERFACE_MAPPING,numberOfInterfaceMatrices,ERR,ERROR,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping to set the number of linear matrices for.
    INTEGER(INTG), INTENT(IN) :: numberOfInterfaceMatrices !<The number of interface matrices to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx,matrix_idx2,variable_idx,number_of_dependent_variables
    INTEGER(INTG), ALLOCATABLE :: OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES(:)
    LOGICAL :: FOUND
    LOGICAL, ALLOCATABLE :: OLD_MATRIX_TRANSPOSE(:)
    REAL(DP), ALLOCATABLE :: OLD_MATRIX_COEFFICIENTS(:)
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_DEPENDENT_TYPE), POINTER :: INTERFACE_DEPENDENT
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("INTERFACE_MAPPING_MATRICES_NUMBER_SET",ERR,ERROR,*999)

    CALL InterfaceMapping_AssertNotFinished(INTERFACE_MAPPING,ERR,ERROR,*999)

    createValuesCache=>INTERFACE_MAPPING%createValuesCache
    IF(ASSOCIATED(createValuesCache)) THEN
      INTERFACE_EQUATIONS=>INTERFACE_MAPPING%interfaceEquations
      IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
        INTERFACE_CONDITION=>INTERFACE_EQUATIONS%INTERFACE_CONDITION
        IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
          SELECT CASE(INTERFACE_CONDITION%METHOD)
          CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
            !Check the number of interface matrices
            IF(numberOfInterfaceMatrices>0) THEN
              INTERFACE_DEPENDENT=>INTERFACE_CONDITION%DEPENDENT
              IF(ASSOCIATED(INTERFACE_DEPENDENT)) THEN
                SELECT CASE(INTERFACE_CONDITION%METHOD)
                CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                  number_of_dependent_variables=INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES
                CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                  number_of_dependent_variables=INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES+1
                END SELECT
                IF(numberOfInterfaceMatrices<=number_of_dependent_variables) THEN
                  !If we need to reallocate and reset all the create values cache arrays and change the number of matrices
                  IF(numberOfInterfaceMatrices/=createValuesCache%numberOfInterfaceMatrices) THEN
                    ALLOCATE(OLD_MATRIX_COEFFICIENTS(createValuesCache%numberOfInterfaceMatrices),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate old matrix coefficients.",ERR,ERROR,*999)
                    ALLOCATE(OLD_MATRIX_TRANSPOSE(createValuesCache%numberOfInterfaceMatrices),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate old matrix transpose.",ERR,ERROR,*999)
                    ALLOCATE(OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES(createValuesCache%numberOfInterfaceMatrices),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate old matrix row field indexes.",ERR,ERROR,*999)
                    OLD_MATRIX_COEFFICIENTS(1:createValuesCache%numberOfInterfaceMatrices)= &
                      createValuesCache%matrixCoefficients(1:createValuesCache%numberOfInterfaceMatrices)
                    OLD_MATRIX_TRANSPOSE(1:createValuesCache%numberOfInterfaceMatrices)= &
                      & createValuesCache%hasTranspose(1:createValuesCache%numberOfInterfaceMatrices)
                    OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES(1:createValuesCache%numberOfInterfaceMatrices)= &
                      & createValuesCache%matrixRowFieldVariableIndices(1:createValuesCache% &
                      & numberOfInterfaceMatrices)
                    IF(ALLOCATED(createValuesCache%matrixCoefficients)) &
                      & DEALLOCATE(createValuesCache%matrixCoefficients)
                    IF(ALLOCATED(createValuesCache%hasTranspose)) &
                      & DEALLOCATE(createValuesCache%hasTranspose)
                    IF(ALLOCATED(createValuesCache%matrixRowFieldVariableIndices)) &
                      & DEALLOCATE(createValuesCache%matrixRowFieldVariableIndices)
                    ALLOCATE(createValuesCache%matrixCoefficients(numberOfInterfaceMatrices),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate matrix coefficients.",ERR,ERROR,*999)
                    ALLOCATE(createValuesCache%hasTranspose(numberOfInterfaceMatrices),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate matrix tranpose.",ERR,ERROR,*999)
                    ALLOCATE(createValuesCache%matrixRowFieldVariableIndices(numberOfInterfaceMatrices),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate matrix row field variable indexes.",ERR,ERROR,*999)
                    IF(numberOfInterfaceMatrices>createValuesCache%numberOfInterfaceMatrices) THEN
                      createValuesCache%matrixCoefficients(1:createValuesCache%numberOfInterfaceMatrices)= &
                        & OLD_MATRIX_COEFFICIENTS(1:createValuesCache%numberOfInterfaceMatrices)
                      createValuesCache%matrixCoefficients(createValuesCache%numberOfInterfaceMatrices+1: &
                        & numberOfInterfaceMatrices)=1.0_DP
                      createValuesCache%hasTranspose(1:createValuesCache%numberOfInterfaceMatrices)= &
                        & OLD_MATRIX_TRANSPOSE(1:createValuesCache%numberOfInterfaceMatrices)
                      createValuesCache%hasTranspose(createValuesCache%numberOfInterfaceMatrices+1: &
                        & numberOfInterfaceMatrices)=.TRUE.
                      createValuesCache%matrixRowFieldVariableIndices(1:createValuesCache% &
                        & numberOfInterfaceMatrices)=OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES(1:createValuesCache% &
                        & numberOfInterfaceMatrices)
                      !Loop through in mesh index order and set the default matrix to variable map to be in mesh index order
                      DO matrix_idx=createValuesCache%numberOfInterfaceMatrices+1,numberOfInterfaceMatrices
                        createValuesCache%matrixRowFieldVariableIndices(matrix_idx)=0
                        DO variable_idx=1,INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES
                          FOUND=.FALSE.
                          DO matrix_idx2=1,createValuesCache%numberOfInterfaceMatrices
                            IF(INTERFACE_DEPENDENT%VARIABLE_MESH_INDICES(variable_idx)==createValuesCache% &
                              matrixRowFieldVariableIndices(matrix_idx2)) THEN
                              FOUND=.TRUE.
                              EXIT
                            ENDIF
                          ENDDO !matrix_idx2
                          IF(.NOT.FOUND) THEN
                            createValuesCache%matrixRowFieldVariableIndices(matrix_idx)=INTERFACE_DEPENDENT% &
                              & VARIABLE_MESH_INDICES(variable_idx)
                          ENDIF
                        ENDDO !variable_idx2
                        IF(createValuesCache%matrixRowFieldVariableIndices(matrix_idx)==0) THEN
                          LOCAL_ERROR="Could not map an interface mesh index for interface matrix "// &
                            & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//"."
                          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                      ENDDO !matrix_idx
                    ELSE
                      createValuesCache%matrixCoefficients(1:numberOfInterfaceMatrices)= &
                        & OLD_MATRIX_COEFFICIENTS(1:numberOfInterfaceMatrices)
                      createValuesCache%hasTranspose(1:numberOfInterfaceMatrices)= &
                        & OLD_MATRIX_TRANSPOSE(1:numberOfInterfaceMatrices)
                      createValuesCache%matrixRowFieldVariableIndices(1:numberOfInterfaceMatrices)= &
                        & OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES(1:numberOfInterfaceMatrices)
                    ENDIF
                    IF(ALLOCATED(OLD_MATRIX_COEFFICIENTS)) DEALLOCATE(OLD_MATRIX_COEFFICIENTS)
                    IF(ALLOCATED(OLD_MATRIX_TRANSPOSE)) DEALLOCATE(OLD_MATRIX_TRANSPOSE)
                    IF(ALLOCATED(OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES)) DEALLOCATE(OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The specified number of interface matrices of "// &
                    & TRIM(NUMBER_TO_VSTRING(numberOfInterfaceMatrices,"*",ERR,ERROR))// &
                    & " is invalid. The number must be <= the number of added dependent variables of "// &
                    & TRIM(NUMBER_TO_VSTRING(number_of_dependent_variables,"*",ERR,ERROR))//"."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Interface condition dependent is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The specified number of interface matrices of "// &
                & TRIM(NUMBER_TO_VSTRING(numberOfInterfaceMatrices,"*",ERR,ERROR))// &
                & " is invalid. The number must be > 0."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The interface condition method of "// &
              & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Interface equations interface condition is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Interface mapping interface equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface mapping create values cache is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("INTERFACE_MAPPING_MATRICES_NUMBER_SET")
    RETURN
999 IF(ALLOCATED(OLD_MATRIX_COEFFICIENTS)) DEALLOCATE(OLD_MATRIX_COEFFICIENTS)
    IF(ALLOCATED(OLD_MATRIX_TRANSPOSE)) DEALLOCATE(OLD_MATRIX_TRANSPOSE)
    IF(ALLOCATED(OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES)) DEALLOCATE(OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES)
    ERRORSEXITS("INTERFACE_MAPPING_MATRICES_NUMBER_SET",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE INTERFACE_MAPPING_MATRICES_NUMBER_SET

  !
  !================================================================================================================================
  !

  !>Sets the transpose flag for the interface matrices. 
  SUBROUTINE INTERFACE_MAPPING_MATRICES_TRANSPOSE_SET(INTERFACE_MAPPING,MATRIX_TRANSPOSE,ERR,ERROR,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping.
    LOGICAL, INTENT(IN) :: MATRIX_TRANSPOSE(:) !<MATRIX_TRANSPOSE(matrix_idx). The interface matrix transpose flag for the matrix_idx'th interface matrix.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("INTERFACE_MAPPING_MATRICES_TRANSPOSE_SET",ERR,ERROR,*999)

    CALL InterfaceMapping_AssertNotFinished(INTERFACE_MAPPING,ERR,ERROR,*999)
    
    createValuesCache=>INTERFACE_MAPPING%createValuesCache
    IF(ASSOCIATED(createValuesCache)) THEN
      INTERFACE_EQUATIONS=>INTERFACE_MAPPING%interfaceEquations
      IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
        INTERFACE_CONDITION=>INTERFACE_EQUATIONS%INTERFACE_CONDITION
        IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
          SELECT CASE(INTERFACE_CONDITION%METHOD)
          CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
            !Check that the number of supplied coefficients matches the number of interface matrices
            IF(SIZE(MATRIX_TRANSPOSE,1)==createValuesCache%numberOfInterfaceMatrices) THEN
              createValuesCache%hasTranspose(1:createValuesCache%numberOfInterfaceMatrices)= &
                MATRIX_TRANSPOSE(1:createValuesCache%numberOfInterfaceMatrices)
            ELSE
              LOCAL_ERROR="Invalid size of matrix tranpose. The size of the supplied array ("// &
                & TRIM(NUMBER_TO_VSTRING(SIZE(MATRIX_TRANSPOSE,1),"*",ERR,ERROR))// &
                & ") must match the number of interface matrices ("// &
                & TRIM(NUMBER_TO_VSTRING(createValuesCache%numberOfInterfaceMatrices,"*",ERR,ERROR))//")."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The interface condition method of "// &
              & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Interface equations interface condition is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Interface mapping interface equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface mapping create values cache is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("INTERFACE_MAPPING_MATRICES_TRANSPOSE_SET")
    RETURN
999 ERRORSEXITS("INTERFACE_MAPPING_MATRICES_TRANSPOSE_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE INTERFACE_MAPPING_MATRICES_TRANSPOSE_SET

  !
  !================================================================================================================================
  !

  !>Sets the coefficient applied to the interface RHS vector.
  SUBROUTINE INTERFACE_MAPPING_RHS_COEFF_SET(INTERFACE_MAPPING,rhsCoefficient,ERR,ERROR,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping to set the RHS coefficent for
    REAL(DP), INTENT(IN) :: rhsCoefficient !<The coefficient applied to the interface RHS vector.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("INTERFACE_MAPPING_RHS_COEFF_SET",ERR,ERROR,*999)

    CALL InterfaceMapping_AssertNotFinished(INTERFACE_MAPPING,ERR,ERROR,*999)
   
    IF(ASSOCIATED(INTERFACE_MAPPING%createValuesCache)) THEN
      IF(INTERFACE_MAPPING%createValuesCache%rhsLagrangeVariableType/=0) THEN
        INTERFACE_MAPPING%createValuesCache%rhsCoefficient=rhsCoefficient
      ELSE
        CALL FlagError("The interface mapping RHS Lagrange variable type has not been set.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface mapping create values cache is not associated",ERR,ERROR,*999)
    ENDIF
       
    EXITS("INTERFACE_MAPPING_RHS_COEFF_SET")
    RETURN
999 ERRORSEXITS("INTERFACE_MAPPING_RHS_COEFF_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE INTERFACE_MAPPING_RHS_COEFF_SET
  
  !
  !================================================================================================================================
  !

  !>Finalises the interface mapping RHS mapping and deallocates all memory
  SUBROUTINE INTERFACE_MAPPING_RHS_MAPPING_FINALISE(rhsMapping,ERR,ERROR,*)

    !Argument variables
    TYPE(InterfaceMappingRHSType), POINTER :: rhsMapping !<A pointer to the RHS mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    ENTERS("INTERFACE_MAPPING_RHS_MAPPING_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(rhsMapping)) THEN
      IF(ALLOCATED(rhsMapping%rhsDOFToInterfaceRowMap)) DEALLOCATE(rhsMapping%rhsDOFToInterfaceRowMap)
      IF(ALLOCATED(rhsMapping%interfaceRowToRHSDOFMap)) DEALLOCATE(rhsMapping%interfaceRowToRHSDOFMap)
      DEALLOCATE(rhsMapping)
    ENDIF
       
    EXITS("INTERFACE_MAPPING_RHS_MAPPING_FINALISE")
    RETURN
999 ERRORSEXITS("INTERFACE_MAPPING_RHS_MAPPING_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE INTERFACE_MAPPING_RHS_MAPPING_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the interface mapping RHS mapping
  SUBROUTINE INTERFACE_MAPPING_RHS_MAPPING_INITIALISE(INTERFACE_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping to initialise the RHS mapping for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS("INTERFACE_MAPPING_RHS_MAPPING_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      IF(ASSOCIATED(INTERFACE_MAPPING%rhsMapping)) THEN
        CALL FlagError("Interface mapping RHS mapping is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(INTERFACE_MAPPING%rhsMapping,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate interface mapping RHS mapping.",ERR,ERROR,*999)
        INTERFACE_MAPPING%rhsMapping%interfaceMapping=>INTERFACE_MAPPING        
        INTERFACE_MAPPING%rhsMapping%rhsVariableType=0
        NULLIFY(INTERFACE_MAPPING%rhsMapping%rhsVariable)
        NULLIFY(INTERFACE_MAPPING%rhsMapping%rhsVariableMapping)
        INTERFACE_MAPPING%rhsMapping%rhsCoefficient=1.0_DP
      ENDIF
    ELSE
      CALL FlagError("Interface mapping is not associated.",ERR,ERROR,*998)
    ENDIF
       
    EXITS("INTERFACE_MAPPING_RHS_MAPPING_INITIALISE")
    RETURN
999 CALL INTERFACE_MAPPING_RHS_MAPPING_FINALISE(INTERFACE_MAPPING%rhsMapping,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("INTERFACE_MAPPING_RHS_MAPPING_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE INTERFACE_MAPPING_RHS_MAPPING_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the mapping between a Lagrange field variable and the interface rhs vector.
  SUBROUTINE INTERFACE_MAPPING_RHS_VARIABLE_TYPE_SET(INTERFACE_MAPPING,rhsVariableType,ERR,ERROR,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping to set the RHS variable type for.
    INTEGER(INTG), INTENT(IN) :: rhsVariableType !<The variable type associated with the interface rhs vector. If the interface condition does not have a rhs vector then the variable type on input should be zero.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_LAGRANGE_TYPE), POINTER :: INTERFACE_LAGRANGE
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(FieldType), POINTER :: LAGRANGE_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("INTERFACE_MAPPING_RHS_VARIABLE_TYPE_SET",ERR,ERROR,*999)

    CALL InterfaceMapping_AssertNotFinished(INTERFACE_MAPPING,ERR,ERROR,*999)

    createValuesCache=>INTERFACE_MAPPING%createValuesCache
    IF(ASSOCIATED(createValuesCache)) THEN
      IF(rhsVariableType==0) THEN
        createValuesCache%rhsLagrangeVariableType=0
      ELSE
        INTERFACE_EQUATIONS=>INTERFACE_MAPPING%interfaceEquations
        IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
          INTERFACE_CONDITION=>INTERFACE_EQUATIONS%INTERFACE_CONDITION
          IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
            SELECT CASE(INTERFACE_CONDITION%METHOD)
            CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
              INTERFACE_LAGRANGE=>INTERFACE_CONDITION%LAGRANGE
              IF(ASSOCIATED(INTERFACE_LAGRANGE)) THEN
                LAGRANGE_FIELD=>INTERFACE_LAGRANGE%LAGRANGE_FIELD
                IF(ASSOCIATED(LAGRANGE_FIELD)) THEN
                  !Check the RHS variable type is not being by the interface matrices
                  IF(createValuesCache%lagrangeVariableType==rhsVariableType) THEN
                    LOCAL_ERROR="The specified RHS variable type of "// &
                      & TRIM(NUMBER_TO_VSTRING(rhsVariableType,"*",ERR,ERROR))// &
                      & " is the same as the Lagrange variable type for the interface matrices."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                  !Check the RHS variable number is defined on the Lagrange field
                  IF(rhsVariableType>=1.AND.rhsVariableType<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                    IF(ASSOCIATED(LAGRANGE_FIELD%variableTypeMap(rhsVariableType)%PTR)) THEN
                      createValuesCache%rhsLagrangeVariableType=rhsVariableType
                    ELSE
                      LOCAL_ERROR="The specified RHS variable type of "// &
                        & TRIM(NUMBER_TO_VSTRING(rhsVariableType,"*",ERR,ERROR))// &
                        & " is not defined on the Lagrange field."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="The specified RHS variable type of "// &
                      & TRIM(NUMBER_TO_VSTRING(rhsVariableType,"*",ERR,ERROR))// &
                      & " is invalid. The number must either be zero or >= 1 and <= "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Lagrange field is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Interface Lagrange is not associated.",ERR,ERROR,*999)
              ENDIF
            CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The interface condition method of "// &
                & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FlagError("Interface equations interface condition is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Interface mapping interface equations is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface mapping create values cache is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("INTERFACE_MAPPING_RHS_VARIABLE_TYPE_SET")
    RETURN
999 ERRORSEXITS("INTERFACE_MAPPING_RHS_VARIABLE_TYPE_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE INTERFACE_MAPPING_RHS_VARIABLE_TYPE_SET
  
  !
  !================================================================================================================================
  !

END MODULE INTERFACE_MAPPING_ROUTINES
