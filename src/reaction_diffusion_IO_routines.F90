!> \file
!> \author Vijay Rajagopal
!> \brief This module handles some mesh/parameter input routines and cmgui output routines for reaction diffusion
!> routines and should be eventually replaces by field_IO_routines.F90
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
!> Contributor(s): Vijay Rajagopal
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

!> Temporary IO routines for fluid mechanics

MODULE REACTION_DIFFUSION_IO_ROUTINES

 USE BaseRoutines
 USE ComputationRoutines
 USE ComputationAccessRoutines
 USE ContextAccessRoutines
 USE EquationsSetConstants
 USE FieldRoutines
 USE FieldAccessRoutines
 USE InputOutput
 USE ISO_VARYING_STRING
 USE Kinds
 USE MeshRoutines
 USE RegionAccessRoutines
 USE Types

#ifndef NOMPIMOD
  USE MPI
#endif

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  PUBLIC REACTION_DIFFUSION_IO_WRITE_CMGUI

CONTAINS

  ! OK
  !================================================================================================================================
  !

  !> Writes solution into cmgui formats exelem and exnode.
  SUBROUTINE REACTION_DIFFUSION_IO_WRITE_CMGUI(REGION, EQUATIONS_SET_GLOBAL_NUMBER, NAME, exportExelem, ERR, ERROR,*)

    !Argument variables
    TYPE(RegionType), INTENT(IN), POINTER :: REGION !<A pointer to the region to get the coordinate system for
    CHARACTER(30),INTENT(IN) :: NAME !<the prefix name of file.
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_GLOBAL_NUMBER
    LOGICAL, INTENT(IN) :: exportExelem !<A flag to indicate whether to write an exelem file.
    INTEGER(INTG) :: ERR !<The error code    
    TYPE(VARYING_STRING):: ERROR !<The error string
    !Local Variables
    TYPE(ComputationEnvironmentType), POINTER :: computationEnvironment
    TYPE(ContextType), POINTER :: context
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET
    TYPE(DomainType), POINTER :: COMPUTATION_DOMAIN
    TYPE(FieldType), POINTER :: SOURCE_FIELD
    REAL(DP) :: NodeXValue,NodeYValue,NodeZValue,NodeUValue
    INTEGER(INTG):: myWorldComputationNodeNumber,NumberOfOutputFields,NumberOfDimensions,NumberOfElements,NumberOfNodes
    INTEGER(INTG):: NumberOfVariableComponents,NumberOfSourceComponents,I,J,K,ValueIndex,NODE_GLOBAL_NUMBER
    INTEGER(INTG) :: NodesInMeshComponent,BasisType,MaxNodesPerElement,NumberOfFieldComponents(3),ELEMENT_GLOBAL_NUMBER
    INTEGER(INTG) :: NODE_LOCAL_NUMBER,numberOfWorldComputationNodes
    INTEGER(INTG),ALLOCATABLE :: ElementNodes(:,:),SimplexOutputHelp(:)
    REAL(DP), ALLOCATABLE :: ElementNodesScales(:,:)
    LOGICAL :: OUTPUT_SOURCE
    TYPE(VARYING_STRING) :: FILENAME !<the prefix name of file.
    CHARACTER(50) :: INTG_STRING,INTG_STRING2


    ENTERS("REACTION_DIFFUSION_IO_WRITE_CMGUI",ERR,ERROR,*999)

    NULLIFY(context)
    CALL Region_ContextGet(region,context,err,error,*999)
    NULLIFY(computationEnvironment)
    CALL Context_ComputationEnvironmentGet(context,computationEnvironment,err,error,*999)

    CALL ComputationEnvironment_NumberOfWorldNodesGet(computationEnvironment,numberOfWorldComputationNodes,err,error,*999)
    CALL ComputationEnvironment_WorldNodeNumberGet(computationEnvironment,myWorldComputationNodeNumber,err,error,*999)

    EQUATIONS_SET => region%equationsSets%equationsSets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr
    NULLIFY(SOURCE_FIELD)
    COMPUTATION_DOMAIN=>REGION%MESHES%MESHES(1) & 
      & %ptr%DECOMPOSITIONS%DECOMPOSITIONS(1)%ptr%DOMAIN(1)%ptr

    NumberOfDimensions = COMPUTATION_DOMAIN%numberOfDimensions
    NumberOfNodes = COMPUTATION_DOMAIN%TOPOLOGY%NODES%numberOfNodes
    NodesInMeshComponent = REGION%meshes%meshes(1)%ptr%topology(1)%ptr%nodes%numberOfNodes
    NumberOfElements = COMPUTATION_DOMAIN%TOPOLOGY%ELEMENTS%numberOfElements
    NumberOfVariableComponents=region%equationsSets%equationsSets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependentField% &
      & variables(1)%numberOfComponents
    NumberOfOutputFields=2
    !determine if there is a source field
    OUTPUT_SOURCE = .FALSE.
    IF( (EQUATIONS_SET%SPECIFICATION(1)==EQUATIONS_SET_CLASSICAL_FIELD_CLASS) &
      & .AND.(EQUATIONS_SET%SPECIFICATION(2)==EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE) &
        & .AND.(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) )THEN
          SOURCE_FIELD=>region%equationsSets%equationsSets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%source%sourceField
          IF( ASSOCIATED(SOURCE_FIELD) ) OUTPUT_SOURCE = .FALSE. !currently set to false to rethink how source is accessed for output
     END IF

    IF( OUTPUT_SOURCE ) THEN
      NumberOfSourceComponents=region%equationsSets%equationsSets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%source%sourceField% &
        & variables(1)%numberOfComponents
      NumberOfOutputFields = NumberOfOutputFields + 1
    !  CALL Field_InterpolationParametersInitialise(SOURCE_FIELD,SOURCE_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
    !  CALL Field_InterpolatedPointsInitialise(SOURCE_INTERPOLATION_PARAMETERS,SOURCE_INTERPOLATED_POINT,ERR,ERROR,*999)
    END IF

    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Writing Nodes...",ERR,ERROR,*999)

    FILENAME="./output/"//TRIM(NAME)//".exnode"

    OPEN(UNIT=myWorldComputationNodeNumber, FILE=CHAR(FILENAME),STATUS='unknown')
    ! WRITING HEADER INFORMATION
    WRITE(myWorldComputationNodeNumber,*) 'Group name: Cell'
    WRITE(INTG_STRING,'(I0)') NumberOfOutputFields
    WRITE(myWorldComputationNodeNumber,*) '#Fields=',TRIM(INTG_STRING)

    ValueIndex=1
    WRITE(INTG_STRING,'(I0)') NumberOfDimensions
    WRITE(myWorldComputationNodeNumber,*) &
      & ' 1) coordinates,  coordinate, rectangular cartesian, #Components=',TRIM(INTG_STRING)
    DO I=1,NumberOfDimensions
      IF(I==1) THEN
        WRITE(INTG_STRING,'(I0)') ValueIndex
        WRITE(myWorldComputationNodeNumber,*) '   x.  Value index= ',TRIM(INTG_STRING),', #Derivatives= 0'
      ELSE IF(I==2) THEN
        WRITE(INTG_STRING,'(I0)') ValueIndex
        WRITE(myWorldComputationNodeNumber,*) '   y.  Value index= ',TRIM(INTG_STRING),', #Derivatives= 0'
      ELSE
        WRITE(INTG_STRING,'(I0)') ValueIndex
        WRITE(myWorldComputationNodeNumber,*) '   z.  Value index= ',TRIM(INTG_STRING),', #Derivatives= 0'
      END IF
      ValueIndex=ValueIndex+1
    END DO

    WRITE(INTG_STRING,'(I0)') NumberOfVariableComponents
    WRITE(myWorldComputationNodeNumber,*) ' 2) dependent, field, rectangular cartesian, #Components=', &
      & TRIM(INTG_STRING)

    DO I=1,NumberOfVariableComponents
      WRITE(INTG_STRING,'(I0)') ValueIndex
      WRITE(INTG_STRING2,'(I0)') I
      WRITE(myWorldComputationNodeNumber,*)  '  ',TRIM(INTG_STRING2),'. Value index= ',TRIM(INTG_STRING), &
        & ', #Derivatives= 0'
      ValueIndex=ValueIndex+1
    END DO

    IF( OUTPUT_SOURCE ) THEN !Watch out that no numbering conflict occurs with Analytic: 4.)
      WRITE(INTG_STRING,'(I0)') NumberOfSourceComponents
      WRITE(myWorldComputationNodeNumber,*) ' 3) source, field, rectangular cartesian, #Components=', &
        & TRIM(INTG_STRING)
      DO I=1,NumberOfSourceComponents
        WRITE(INTG_STRING,'(I0)') ValueIndex
        WRITE(INTG_STRING2,'(I0)') I
        WRITE(myWorldComputationNodeNumber,*)  '   ',TRIM(INTG_STRING2),'.  Value index= ', &
          & TRIM(INTG_STRING),', #Derivatives= 0'
        ValueIndex=ValueIndex+1
      END DO
    END IF

    !WRITE OUT NODE VALUES
    DO I = 1,NumberOfNodes
      NODE_GLOBAL_NUMBER = COMPUTATION_DOMAIN%TOPOLOGY%NODES%NODES(I)%globalNumber
      NodeXValue = region%equationsSets%equationsSets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%geometry%geometricField%variables(1) &
        & %parameterSets%parameterSets(1)%ptr%parameters%cmiss%dataDP(I)
      IF(NumberOfDimensions==2 .OR. NumberOfDimensions==3) THEN
        NodeYValue = region%equationsSets%equationsSets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%geometry%geometricField%variables(1) &
          & %parameterSets%parameterSets(1)%ptr%parameters%cmiss%dataDP(I+NumberOfNodes)
      ENDIF
      IF(NumberOfDimensions==3) THEN
        NodeZValue = region%equationsSets%equationsSets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%geometry%geometricField%variables(1) &
          & %parameterSets%parameterSets(1)%ptr%parameters%cmiss%dataDP(I+(2*NumberOfNodes))
      ENDIF
      NodeUValue=region%equationsSets%equationsSets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependentField% &
        & variables(1)%parameterSets%parameterSets(1)%ptr%parameters%cmiss%dataDP(I)

      WRITE(myWorldComputationNodeNumber,*) ' Node: ',NODE_GLOBAL_NUMBER
      WRITE(myWorldComputationNodeNumber,'("    ", es25.16 )')NodeXValue

      IF(NumberOfDimensions==2 .OR. NumberOfDimensions==3) THEN
        WRITE(myWorldComputationNodeNumber,'("    ", es25.16 )')NodeYValue
      END IF

      IF(NumberOfDimensions==3) THEN
        WRITE(myWorldComputationNodeNumber,'("    ", es25.16 )')NodeZValue
      END IF
      WRITE(myWorldComputationNodeNumber,'("    ", es25.16 )')NodeUValue

      IF( (EQUATIONS_SET%SPECIFICATION(1)==EQUATIONS_SET_CLASSICAL_FIELD_CLASS) &
        & .AND.(EQUATIONS_SET%SPECIFICATION(2)==EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE) &
          & .AND.(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) )THEN
          !source field
          IF( OUTPUT_SOURCE ) THEN
            !NodeSourceValue = SOURCE_INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,1)
            !WRITE(myWorldComputationNodeNumber,'("    ", es25.16 )')NodeSourceValue
          END IF
      END IF
    END DO !nodes I
    CLOSE(myWorldComputationNodeNumber)

    !OUTPUT ELEMENTS IN CURRENT DOMAIN
    MaxNodesPerElement=COMPUTATION_DOMAIN%TOPOLOGY%ELEMENTS%ELEMENTS(1)%basis%numberOfElementParameters
    BasisType = 1
    IF(NumberOfDimensions==2) THEN
      IF(MaxNodesPerElement==4.OR.MaxNodesPerElement==9.OR.MaxNodesPerElement==16) THEN
        BasisType=1
      ELSEIF(MaxNodesPerElement==3.OR.MaxNodesPerElement==6.OR.MaxNodesPerElement==10) THEN
        BasisType=2
      ENDIF
    ELSEIF(NumberOfDimensions==3) THEN
      BasisType=region%equationsSets%equationsSets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%equations% &
        & interpolation%geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr%bases(1)%ptr%type
    ENDIF

    IF(exportExelem) THEN
      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Writing Elements...",ERR,ERROR,*999)
      FILENAME="./output/"//TRIM(NAME)//".exelem"
      OPEN(UNIT=myWorldComputationNodeNumber, FILE=CHAR(FILENAME),STATUS='unknown')
      WRITE(myWorldComputationNodeNumber,*) 'World name: Cell'
      IF (BasisType==1) THEN !lagrange basis in 1 and 2D
        WRITE(INTG_STRING,'(I0)') NumberOfDimensions
        WRITE(myWorldComputationNodeNumber,*) 'Shape.  Dimension= ',TRIM(INTG_STRING)
        WRITE(myWorldComputationNodeNumber,*) '#Scale factor sets= 1'
        IF(NumberOfDimensions==1) THEN
          WRITE(INTG_STRING,'(I0)') MaxNodesPerElement
          WRITE(myWorldComputationNodeNumber,*) 'q.Lagrange, #Scale factors=',TRIM(INTG_STRING)
        ELSE IF (NumberOfDimensions==2) THEN
          IF(MaxNodesPerElement==4) THEN

            WRITE(INTG_STRING,'(I0)') MaxNodesPerElement
            WRITE(myWorldComputationNodeNumber,*) &
              & 'l.Lagrange*l.Lagrange, #Scale factors=',TRIM(INTG_STRING) !linear lagrange
          ELSE IF(MaxNodesPerElement==9) THEN
            WRITE(INTG_STRING,'(I0)') MaxNodesPerElement
            WRITE(myWorldComputationNodeNumber,*) &
              & 'q.Lagrange*q.Lagrange, #Scale factors=',TRIM(INTG_STRING) !quadratic lagrange
          ELSE IF(MaxNodesPerElement==16) THEN
            WRITE(INTG_STRING,'(I0)') MaxNodesPerElement
            WRITE(myWorldComputationNodeNumber,*) &
              & 'c.Lagrange*c.Lagrange, #Scale factors=',TRIM(INTG_STRING) !cubic lagrange
          END IF
        ELSE !three dimensions
          IF(MaxNodesPerElement==8) THEN
            WRITE(INTG_STRING,'(I0)') MaxNodesPerElement
            WRITE(myWorldComputationNodeNumber,*) &
              & 'l.Lagrange*l.Lagrange*l.Lagrange, #Scale factors=',TRIM(INTG_STRING)
          ELSE IF(MaxNodesPerElement==27) THEN
            WRITE(INTG_STRING,'(I0)') MaxNodesPerElement
            WRITE(myWorldComputationNodeNumber,*) &
              & 'q.Lagrange*q.Lagrange*q.Lagrange, #Scale factors=',TRIM(INTG_STRING)
          ELSE IF(MaxNodesPerElement==64) THEN
            WRITE(INTG_STRING,'(I0)') MaxNodesPerElement
            WRITE(myWorldComputationNodeNumber,*) &
              & 'c.Lagrange*c.Lagrange*c.Lagrange, #Scale factors=',TRIM(INTG_STRING)
          END IF
        END IF
      ELSEIF(BasisType==2) THEN
        IF(NumberOfDimensions==2) THEN
          WRITE(myWorldComputationNodeNumber,*) 'Shape.  Dimension=', &
            & NumberOfDimensions,', simplex(2)*simplex'
          IF(MaxNodesPerElement==3) THEN
            WRITE(myWorldComputationNodeNumber,*) '#Scale factor sets= 1'
            WRITE(INTG_STRING,'(I0)') MaxNodesPerElement
            WRITE(myWorldComputationNodeNumber,*)  &
              & ' l.simplex(2)*l.simplex, #Scale factors= ', TRIM(INTG_STRING)
          ELSE IF(MaxNodesPerElement==6) THEN
            WRITE(myWorldComputationNodeNumber,*) '#Scale factor sets= 1'
            WRITE(INTG_STRING,'(I0)') MaxNodesPerElement
            WRITE(myWorldComputationNodeNumber,*) &
              & ' l.simplex(2)*l.simplex, #Scale factors= ', TRIM(INTG_STRING)
          ELSE IF (MaxNodesPerElement== 10 ) THEN
            WRITE(myWorldComputationNodeNumber,*) '#Scale factor sets= 1'
            WRITE(INTG_STRING,'(I0)') MaxNodesPerElement
            WRITE(myWorldComputationNodeNumber,*) &
              & ' q.simplex(2)*q.simplex, #Scale factors= ', TRIM(INTG_STRING)
          ENDIF
        ELSE IF(NumberOfDimensions==3) THEN
          WRITE(INTG_STRING2,'(I0)') NumberOfDimensions
          WRITE(myWorldComputationNodeNumber,*) &
            & 'Shape.  Dimension=',TRIM(INTG_STRING2),', simplex(2;3)*simplex*simplex'
          IF(MaxNodesPerElement==4) THEN
            WRITE(myWorldComputationNodeNumber,*) &
              & '#Scale factor sets= 1'
            WRITE(INTG_STRING,'(I0)') MaxNodesPerElement
            WRITE(myWorldComputationNodeNumber,*) &
              & ' l.simplex(2;3)*l.simplex*l.simplex, #Scale factors= ', TRIM(INTG_STRING)
          ELSE IF (MaxNodesPerElement== 10 ) THEN
            WRITE(myWorldComputationNodeNumber,*) '#Scale factor sets= 1'
            WRITE(INTG_STRING,'(I0)') MaxNodesPerElement
            WRITE(myWorldComputationNodeNumber,*) &
              & ' q.simplex(2;3)*q.simplex*q.simplex, #Scale factors= ', TRIM(INTG_STRING)
          ELSE IF(MaxNodesPerElement==20) THEN
            WRITE(myWorldComputationNodeNumber,*) '#Scale factor sets= 1'
            WRITE(INTG_STRING,'(I0)') MaxNodesPerElement
            WRITE(myWorldComputationNodeNumber,*) &
              & ' q.simplex(2;3)*q.simplex*q.simplex, #Scale factors= ', TRIM(INTG_STRING)
          ENDIF
        ELSE
          WRITE(myWorldComputationNodeNumber,*) '#Scale factor sets= 0'
        END IF

      END IF
      WRITE(INTG_STRING,'(I0)') MaxNodesPerElement
      WRITE(myWorldComputationNodeNumber,*) '#Nodes= ',TRIM(INTG_STRING)
      WRITE(INTG_STRING,'(I0)') NumberOfOutputFields
      WRITE(myWorldComputationNodeNumber,*) '#Fields= ',TRIM(INTG_STRING)
      NumberOfFieldComponents(1) = NumberOfDimensions
      NumberOfFieldComponents(2) = NumberOfVariableComponents
      NumberOfFieldComponents(3) = NumberOfSourceComponents
      DO I=1,NumberOfOutputFields
        IF(I==1)THEN
          WRITE(INTG_STRING,'(I0)') NumberOfDimensions
          WRITE(myWorldComputationNodeNumber,*) &
            & ' 1) coordinates,  coordinate, rectangular cartesian, #Components= ',TRIM(INTG_STRING)
        ELSE IF(I==2) THEN
          WRITE(INTG_STRING,'(I0)') NumberOfVariableComponents
          WRITE(myWorldComputationNodeNumber,*) &
          & ' 2) dependent,  field,  rectangular cartesian, #Components= ',TRIM(INTG_STRING)
        ELSE IF(I==3) THEN
          WRITE(INTG_STRING,'(I0)') NumberOfSourceComponents
          WRITE(myWorldComputationNodeNumber,*) &
            & ' 3) source,  field,  rectangular cartesian, #Components= ',TRIM(INTG_STRING)
        END IF

        DO J=1,NumberOfFieldComponents(I)
          IF(NumberOfDimensions==1) THEN
            IF(I==1)THEN
              IF(J==1) THEN
                  WRITE(myWorldComputationNodeNumber,*)'   x.   l.Lagrange, no modify, standard node based.'
              ELSE IF(J==2) THEN
                  WRITE(myWorldComputationNodeNumber,*)'   y.   l.Lagrange, no modify, standard node based.'
              ELSE IF(J==3) THEN
                  WRITE(myWorldComputationNodeNumber,*)'   z.   l.Lagrange, no modify, standard node based.'
              END IF
            ELSE
              WRITE(myWorldComputationNodeNumber,*) &
                & '   ',J,'.   l.Lagrange, no modify, standard node based.'
            END IF
          ELSE IF(NumberOfDimensions==2) THEN
            IF(I==1)THEN
              IF(J==1) THEN
                IF(MaxNodesPerElement==4)THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   x.   l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==9) THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   x.   q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==16)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   x.   c.Lagrange*c.Lagrange, no modify, standard node based.'

                ELSE IF(MaxNodesPerElement==3)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   x.  l.simplex(2)*l.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==6)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   x.  q.simplex(2)*q.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==10)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   x.  c.simplex(2)*c.simplex, no modify, standard node based.'
                END IF
              ELSE IF(J==2) THEN
                IF(MaxNodesPerElement==4) THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   y.   l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==9)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   y.   q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==16)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   y.   c.Lagrange*c.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==3)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   y.  l.simplex(2)*l.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==6)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   y.  q.simplex(2)*q.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==10)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                  & '   y.  c.simplex(2)*c.simplex, no modify, standard node based.'
                END IF
              ELSE IF(J==3) THEN
                IF(MaxNodesPerElement==4) THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   z.   l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==9)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   z.   q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==16)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   z.   c.Lagrange*c.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==3)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   z.  l.simplex(2)*l.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==6)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   z.  q.simplex(2)*q.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==10)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   z.  c.simplex(2)*c.simplex, no modify, standard node based.'
                END IF
              END IF
            ELSE
                IF(MaxNodesPerElement==4) THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',J,'.   l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==9)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',J,'.   q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==16)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',J,'.   c.Lagrange*c.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==3)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',J,'.  l.simplex(2)*l.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==6)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',J,'.  q.simplex(2)*q.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==10)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',J,'.  c.simplex(2)*c.simplex, no modify, standard node based.'
                END IF
            END IF
          ELSE IF(NumberOfDimensions==3) THEN
            IF(I==1)THEN
              IF(J==1) THEN
                IF(MaxNodesPerElement==8) THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   x.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==27)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   x.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==64)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   x.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==4)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   x.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==10)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   x.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==20)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   x.  c.simplex(2;3)*c.simplex*c.simplex, no modify, standard node based.'
                END IF
              ELSE IF(J==2) THEN
                IF(MaxNodesPerElement==8) THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   y.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==27)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   y.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==64)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   y.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==4)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   y.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==10)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   y.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==20)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   y.  c.simplex(2;3)*c.simplex*c.simplex, no modify, standard node based.'
                END IF
              ELSE IF(J==3) THEN
                IF(MaxNodesPerElement==8) THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   z.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==27)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   z.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==64)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   z.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==4)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   z.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==10)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   z.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==20)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   z.  c.simplex(2;3)*c.simplex*c.simplex, no modify, standard node based.'
                END IF
              END IF
            ELSE
                IF(MaxNodesPerElement==8) THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',J,'.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==27)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',J,'.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==64)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',J,'.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==4)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',J,'.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==10)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',J,'.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.'
                ELSE IF(MaxNodesPerElement==20)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',J,'.  c.simplex(2;3)*c.simplex*c.simplex, no modify, standard node based.'
                END IF
            END IF
          END IF
          WRITE(INTG_STRING,'(I0)') MaxNodesPerElement
          WRITE(myWorldComputationNodeNumber,*) '   #Nodes= ',TRIM(INTG_STRING)

          DO K = 1,MaxNodesPerElement
            WRITE(INTG_STRING,'(I0)') K
            WRITE(myWorldComputationNodeNumber,*) '    ',TRIM(INTG_STRING),'.  #Values=1'
            WRITE(myWorldComputationNodeNumber,*) '     Value indices:     1'
            WRITE(myWorldComputationNodeNumber,*) '     Scale factor indices:   ',TRIM(INTG_STRING)
          END DO
        END DO !J loop
      END DO !I loop
      IF(.NOT.ALLOCATED(ElementNodes)) ALLOCATE(ElementNodes(NumberOfElements,MaxNodesPerElement))
      IF(.NOT.ALLOCATED(ElementNodesScales)) ALLOCATE(ElementNodesScales(NumberOfElements,MaxNodesPerElement))
      DO I=1,NumberOfElements
        ELEMENT_GLOBAL_NUMBER=COMPUTATION_DOMAIN%DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(K)%globalNumber
        DO J=1,MaxNodesPerElement
          NODE_LOCAL_NUMBER=COMPUTATION_DOMAIN%TOPOLOGY%ELEMENTS%ELEMENTS(I)%elementNodes(J)
          NODE_GLOBAL_NUMBER=COMPUTATION_DOMAIN%MAPPINGS%NODES%localToGlobalMap(NODE_LOCAL_NUMBER)
          ElementNodes(I,J)=NODE_GLOBAL_NUMBER
          ElementNodesScales(I,J)=1.0000000000000000E+00
        END DO
      END DO


      DO K=1,NumberOfElements
        ELEMENT_GLOBAL_NUMBER=COMPUTATION_DOMAIN%DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(K)%globalNumber
        IF (BasisType==1) THEN
          WRITE(INTG_STRING,'(I0)') ELEMENT_GLOBAL_NUMBER
          WRITE(myWorldComputationNodeNumber,*) 'Element:     ', TRIM(INTG_STRING),' 0  0'
          WRITE(myWorldComputationNodeNumber,*) '   Nodes:'
          WRITE(myWorldComputationNodeNumber,*) '   ', ElementNodes(K,1:MaxNodesPerElement)
          WRITE(myWorldComputationNodeNumber,*) '   Scale factors:'
          WRITE(myWorldComputationNodeNumber,*) '   ',ElementNodesScales(K,1:MaxNodesPerElement)

        ELSEIF(BasisType==2) THEN
          IF(.NOT.ALLOCATED(SimplexOutputHelp)) ALLOCATE(SimplexOutputHelp(MaxNodesPerElement))
          IF(NumberOfDimensions==2)THEN
            SimplexOutputHelp(1)=ElementNodes(K,1)
            SimplexOutputHelp(2)=ElementNodes(K,2)
            SimplexOutputHelp(3)=ElementNodes(K,3)
          ELSE IF(NumberOfDimensions==3) THEN
            SimplexOutputHelp(1)=ElementNodes(K,1)
            SimplexOutputHelp(2)=ElementNodes(K,4)
            SimplexOutputHelp(3)=ElementNodes(K,2)
            SimplexOutputHelp(4)=ElementNodes(K,3)
          END IF
          WRITE(INTG_STRING,'(I0)') ELEMENT_GLOBAL_NUMBER
          WRITE(myWorldComputationNodeNumber,*) 'Element:     ', TRIM(INTG_STRING),' 0  0'
          WRITE(myWorldComputationNodeNumber,*) '   Nodes:'
          WRITE(myWorldComputationNodeNumber,*) '   ', SimplexOutputHelp
          WRITE(myWorldComputationNodeNumber,*) '   Scale factors:'
          WRITE(myWorldComputationNodeNumber,*) '   ',ElementNodesScales(K,1:MaxNodesPerElement)
        END IF
      ENDDO
      CLOSE(myWorldComputationNodeNumber)
    END IF !exportExelem flag check 

    EXITS("REACTION_DIFFUSION_IO_WRITE_CMGUI")
    RETURN
999 ERRORSEXITS("REACTION_DIFFUSION_IO_WRITE_CMGUI",ERR,ERROR)
    RETURN 1

  END SUBROUTINE REACTION_DIFFUSION_IO_WRITE_CMGUI

  !================================================================================================================================
  !


END MODULE REACTION_DIFFUSION_IO_ROUTINES
