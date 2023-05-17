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
!> Contributor(s): Vijay Rajagopal, Chris Bradley
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

MODULE ReactionDiffusionIORoutines

 USE BaseRoutines
 USE BasisAccessRoutines
 USE ComputationRoutines
 USE ComputationAccessRoutines
 USE ContextAccessRoutines
 USE DecompositionAccessRoutines
 USE DomainMappings
 USE EquationsSetAccessRoutines
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

  PUBLIC ReactionDiffusion_IOWriteCMGUI

CONTAINS

  ! OK
  !================================================================================================================================
  !

  !> Writes solution into cmgui formats exelem and exnode.
  SUBROUTINE ReactionDiffusion_IOWriteCMGUI(region,equationsSetUserNumber,name,exportExelem,err,error,*)

    !Argument variables
    TYPE(RegionType), INTENT(IN), POINTER :: region !<A pointer to the region to get the coordinate system for
    CHARACTER(30),INTENT(IN) :: name !<the prefix name of file.
    INTEGER(INTG), INTENT(IN) :: equationsSetUserNumber !<The user number of the equations set. TODO:: WHY NOT PASS equations set?
    LOGICAL, INTENT(IN) :: exportExelem !<A flag to indicate whether to write an exelem file.
    INTEGER(INTG) :: err !<The error code    
    TYPE(VARYING_STRING):: error !<The error string
    !Local Variables
    INTEGER(INTG) :: basisShapeType,componentIdx,dependentVariableType,dimensionIdx,dofIndex,elementIdx,esSpecification(3), &
      & fieldIdx,globalElementNumber,globalNodeNumber,localNodeNumber,maxNodesPerElement,myWorldComputationNodeNumber,nodeIdx, &
      & numberOfDimensions,numberOfElements,numberOfFieldComponents(3),numberOfNodes,numberOfOutputFields, &
      & numberOfSourceComponents,numberOfVariableComponents,numberOfWorldComputationNodes,valueIndex
    INTEGER(INTG), ALLOCATABLE :: elementNodes(:,:),simplexOutputHelp(:)
    REAL(DP) :: nodeXValue,nodeYValue,nodeZValue,nodeUValue
    REAL(DP), ALLOCATABLE :: elementNodesScales(:,:)
    REAL(DP), POINTER :: dependentParameterData(:),geometricParameterData(:),sourceParameterData(:)
    LOGICAL :: outputSource
    CHARACTER(50) :: intgString,intgString2
    TYPE(BasisType), POINTER :: elementBasis
    TYPE(ComputationEnvironmentType), POINTER :: computationEnvironment
    TYPE(ContextType), POINTER :: context
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainMappingType), POINTER :: nodesMapping
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(FieldType), POINTER :: dependentField,geometricField,sourceField
    TYPE(FieldVariableType), POINTER :: dependentVariable,geometricVariable,sourceVariable
    TYPE(VARYING_STRING) :: filename !<the prefix name of file.
    
    ENTERS("ReactionDiffusion_IOWriteCMGUI",err,error,*996)

    NULLIFY(equationsSet)
    CALL Region_EquationsSetGet(region,equationsSetUserNumber,equationsSet,err,error,*996)
    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*996)
      
    IF((esSpecification(1)==EQUATIONS_SET_CLASSICAL_FIELD_CLASS).AND. &
      & (esSpecification(2)==EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE)) THEN

      NULLIFY(geometricParameterData)
      NULLIFY(dependentParameterData)
      NULLIFY(sourceParameterData)
      
      NULLIFY(context)
      CALL Region_ContextGet(region,context,err,error,*999)
      NULLIFY(computationEnvironment)
      CALL Context_ComputationEnvironmentGet(context,computationEnvironment,err,error,*999)

      CALL ComputationEnvironment_NumberOfWorldNodesGet(computationEnvironment,numberOfWorldComputationNodes,err,error,*999)
      CALL ComputationEnvironment_WorldNodeNumberGet(computationEnvironment,myWorldComputationNodeNumber,err,error,*999)

      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(geometricVariable)
      CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
      CALL FieldVariable_ParameterSetDataGet(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameterData,err,error,*999)
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(dependentVariable)
      CALL Field_VariableIndexGet(dependentField,1,dependentVariable,dependentVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfVariableComponents,err,error,*999)
      CALL FieldVariable_ParameterSetDataGet(dependentVariable,FIELD_VALUES_SET_TYPE,dependentParameterData,err,error,*999)
      NULLIFY(decomposition)
      CALL Field_DecompositionGet(geometricfield,decomposition,err,error,*999)      
      CALL Decomposition_NumberOfDimensionsGet(decomposition,numberOfDimensions,err,error,*999)
      NULLIFY(decompositionTopology)
      CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
      NULLIFY(decompositionElements)
      CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
      NULLIFY(domain)
      CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
      NULLIFY(domainTopology)
      CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
      NULLIFY(domainNodes)
      CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
      CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
      NULLIFY(domainElements)
      CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
      CALL DomainElements_NumberOfElementsGet(domainElements,numberOfElements,err,error,*999)
      NULLIFY(domainMappings)
      CALL Domain_DomainMappingsGet(domain,domainMappings,err,error,*999)
      NULLIFY(nodesMapping)
      CALL DomainMappings_NodesMappingGet(domainMappings,nodesMapping,err,error,*999)

      numberOfOutputFields=2
      NULLIFY(sourceField)
      NULLIFY(sourceVariable)
      !Determine if there is a source field
      outputSource=.FALSE.
      IF(esSpecification(3)==EQUATIONS_SET_CELLML_NOSPLIT_LINEAR_GEN_REACT_DIFF_SUBTYPE) THEN
        CALL EquationsSet_SourceFieldExists(equationsSet,sourceField,err,error,*999)
        IF(ASSOCIATED(sourceField)) outputSource=.FALSE. !currently set to false to rethink how source is accessed for output
      ENDIF

      IF(outputSource) THEN
        CALL Field_VariableGet(sourceField,FIELD_U_VARIABLE_TYPE,sourceVariable,err,error,*999)
        CALL FieldVariable_NumberOfComponentsGet(sourceVariable,numberOfSourceComponents,err,error,*999)
        CALL FieldVariable_ParameterSetDataGet(sourceVariable,FIELD_VALUES_SET_TYPE,sourceParameterData,err,error,*999)
        numberOfOutputFields=numberOfOutputFields+1
      ENDIF

      CALL WriteString(GENERAL_OUTPUT_TYPE,"Writing Nodes...",err,error,*999)

      filename="./output/"//TRIM(name)//".exnode"

      OPEN(UNIT=myWorldComputationNodeNumber, FILE=CHAR(filename),STATUS='unknown')
      !Write header information
      WRITE(myWorldComputationNodeNumber,*) 'Group name: Cell'
      WRITE(intgString,'(I0)') numberOfOutputFields
      WRITE(myWorldComputationNodeNumber,*) '#Fields=',TRIM(intgString)

      valueIndex=1
      WRITE(intgString,'(I0)') numberOfDimensions
      WRITE(myWorldComputationNodeNumber,*) &
        & ' 1) coordinates,  coordinate, rectangular cartesian, #Components=',TRIM(intgString)
      DO dimensionIdx=1,numberOfDimensions
        IF(dimensionIdx==1) THEN
          WRITE(intgString,'(I0)') valueIndex
          WRITE(myWorldComputationNodeNumber,*) '   x.  Value index= ',TRIM(intgString),', #Derivatives= 0'
        ELSE IF(dimensionIdx==2) THEN
          WRITE(intgString,'(I0)') valueIndex
          WRITE(myWorldComputationNodeNumber,*) '   y.  Value index= ',TRIM(intgString),', #Derivatives= 0'
        ELSE
          WRITE(intgString,'(I0)') valueIndex
          WRITE(myWorldComputationNodeNumber,*) '   z.  Value index= ',TRIM(intgString),', #Derivatives= 0'
        END IF
        valueIndex=valueIndex+1
      ENDDO !dimensionsIdx

      WRITE(intgString,'(I0)') numberOfVariableComponents
      WRITE(myWorldComputationNodeNumber,*) ' 2) dependent, field, rectangular cartesian, #Components=', &
        & TRIM(intgString)

      DO componentIdx=1,numberOfVariableComponents
        WRITE(intgString,'(I0)') valueIndex
        WRITE(intgString2,'(I0)') componentIdx
        WRITE(myWorldComputationNodeNumber,*)  '  ',TRIM(intgString2),'. Value index= ',TRIM(intgString), &
          & ', #Derivatives= 0'
        valueIndex=valueIndex+1
      ENDDO !componentIdx

      IF(outputSource) THEN !Watch out that no numbering conflict occurs with Analytic: 4.)
        WRITE(intgString,'(I0)') numberOfSourceComponents
        WRITE(myWorldComputationNodeNumber,*) ' 3) source, field, rectangular cartesian, #Components=', &
          & TRIM(intgString)
        DO componentIdx=1,numberOfSourceComponents
          WRITE(intgString,'(I0)') valueIndex
          WRITE(intgString2,'(I0)') componentIdx
          WRITE(myWorldComputationNodeNumber,*)  '   ',TRIM(intgString2),'.  Value index= ', &
            & TRIM(intgString),', #Derivatives= 0'
          valueIndex=valueIndex+1
        ENDDO !componentIdx
      ENDIF

      !Write out node values
      DO nodeIdx=1,numberOfNodes
        CALL DomainNodes_NodeGlobalNumberGet(domainNodes,nodeIdx,globalNodeNumber,err,error,*999)
        CALL FieldVariable_LocalNodeDOFGet(geometricVariable,1,1,nodeIdx,1,dofIndex,err,error,*999)
        nodeXValue=geometricParameterData(dofIndex)
        IF(numberOfDimensions>=2) THEN
          CALL FieldVariable_LocalNodeDOFGet(geometricVariable,1,1,nodeIdx,2,dofIndex,err,error,*999)
          nodeYValue=geometricParameterData(dofIndex)
          IF(numberOfDimensions==3) THEN
            CALL FieldVariable_LocalNodeDOFGet(geometricVariable,1,1,nodeIdx,3,dofIndex,err,error,*999)
            nodeZValue=geometricParameterData(dofIndex)
          ENDIF
        ENDIF
        CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,1,nodeIdx,1,dofIndex,err,error,*999)
        nodeUValue=dependentParameterData(dofIndex)

        WRITE(myWorldComputationNodeNumber,*) ' Node: ',globalNodeNumber
        WRITE(myWorldComputationNodeNumber,'("    ",ES25.16)') nodeXValue

        IF(numberOfDimensions>=2) THEN
          WRITE(myWorldComputationNodeNumber,'("    ",ES25.16)') nodeYValue
          IF(numberOfDimensions==3) THEN
            WRITE(myWorldComputationNodeNumber,'("    ",ES25.16)') nodeZValue
          ENDIF
        ENDIF
        WRITE(myWorldComputationNodeNumber,'("    ",ES25.16)') nodeUValue

        IF(esSpecification(3)==EQUATIONS_SET_CELLML_NOSPLIT_LINEAR_GEN_REACT_DIFF_SUBTYPE) THEN
          !Source field
          IF(outputSource) THEN
            !CALL FieldVariable_LocalNodeDOFGet(sourceVariable,1,1,nodeIdx,1,dofIndex,err,error,*999)
            !nodeSourceValue=sourceParameterData(dofIndex)
            !WRITE(myWorldComputationNodeNumber,'("    ",ES25.16)') nodeSourceValue
          ENDIF
        ENDIF
      ENDDO !nodeIdx
      CLOSE(myWorldComputationNodeNumber)

      !Output elements in the current domainOUTPUT ELEMENTS IN CURRENT DOMAIN
      NULLIFY(elementBasis)
      CALL DomainElements_ElementBasisGet(domainElements,1,elementBasis,err,error,*999)
      CALL Basis_NumberOfElementParametersGet(elementBasis,maxNodesPerElement,err,error,*999)
      CALL Basis_TypeGet(elementBasis,basisShapeType,err,error,*999)

      IF(exportExelem) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Writing Elements...",err,error,*999)
        filename="./output/"//TRIM(name)//".exelem"
        OPEN(UNIT=myWorldComputationNodeNumber, FILE=CHAR(filename),STATUS='unknown')
        WRITE(myWorldComputationNodeNumber,*) 'World name: Cell'
        IF (basisShapeType==BASIS_LAGRANGE_HERMITE_TP_TYPE) THEN !lagrange basis in 1 and 2D
          WRITE(intgString,'(I0)') numberOfDimensions
          WRITE(myWorldComputationNodeNumber,*) 'Shape.  Dimension= ',TRIM(intgString)
          WRITE(myWorldComputationNodeNumber,*) '#Scale factor sets= 1'
          IF(numberOfDimensions==1) THEN
            WRITE(intgString,'(I0)') maxNodesPerElement
            WRITE(myWorldComputationNodeNumber,*) 'q.Lagrange, #Scale factors=',TRIM(intgString)
          ELSE IF(numberOfDimensions==2) THEN
            IF(maxNodesPerElement==4) THEN              
              WRITE(intgString,'(I0)') maxNodesPerElement
              WRITE(myWorldComputationNodeNumber,*) &
                & 'l.Lagrange*l.Lagrange, #Scale factors=',TRIM(intgString) !linear lagrange
            ELSE IF(maxNodesPerElement==9) THEN
              WRITE(intgString,'(I0)') maxNodesPerElement
              WRITE(myWorldComputationNodeNumber,*) &
                & 'q.Lagrange*q.Lagrange, #Scale factors=',TRIM(intgString) !quadratic lagrange
            ELSE IF(maxNodesPerElement==16) THEN
              WRITE(intgString,'(I0)') maxNodesPerElement
              WRITE(myWorldComputationNodeNumber,*) &
                & 'c.Lagrange*c.Lagrange, #Scale factors=',TRIM(intgString) !cubic lagrange
            ENDIF
          ELSE !three dimensions
            IF(maxNodesPerElement==8) THEN
              WRITE(intgString,'(I0)') maxNodesPerElement
              WRITE(myWorldComputationNodeNumber,*) &
                & 'l.Lagrange*l.Lagrange*l.Lagrange, #Scale factors=',TRIM(intgString)
            ELSE IF(maxNodesPerElement==27) THEN
              WRITE(intgString,'(I0)') maxNodesPerElement
              WRITE(myWorldComputationNodeNumber,*) &
                & 'q.Lagrange*q.Lagrange*q.Lagrange, #Scale factors=',TRIM(intgString)
            ELSE IF(maxNodesPerElement==64) THEN
              WRITE(intgString,'(I0)') maxNodesPerElement
              WRITE(myWorldComputationNodeNumber,*) &
                & 'c.Lagrange*c.Lagrange*c.Lagrange, #Scale factors=',TRIM(intgString)
            ENDIF
          ENDIF
        ELSE IF(basisShapeType==BASIS_SIMPLEX_TYPE) THEN
          IF(numberOfDimensions==2) THEN
            WRITE(myWorldComputationNodeNumber,*) 'Shape.  Dimension=', &
              & numberOfDimensions,', simplex(2)*simplex'
            IF(maxNodesPerElement==3) THEN
              WRITE(myWorldComputationNodeNumber,*) '#Scale factor sets= 1'
              WRITE(intgString,'(I0)') maxNodesPerElement
              WRITE(myWorldComputationNodeNumber,*)  &
                & ' l.simplex(2)*l.simplex, #Scale factors= ', TRIM(intgString)
            ELSE IF(maxNodesPerElement==6) THEN
              WRITE(myWorldComputationNodeNumber,*) '#Scale factor sets= 1'
              WRITE(intgString,'(I0)') maxNodesPerElement
              WRITE(myWorldComputationNodeNumber,*) &
                & ' l.simplex(2)*l.simplex, #Scale factors= ', TRIM(intgString)
            ELSE IF (maxNodesPerElement== 10 ) THEN
              WRITE(myWorldComputationNodeNumber,*) '#Scale factor sets= 1'
              WRITE(intgString,'(I0)') maxNodesPerElement
              WRITE(myWorldComputationNodeNumber,*) &
                & ' q.simplex(2)*q.simplex, #Scale factors= ', TRIM(intgString)
            ENDIF
          ELSE IF(numberOfDimensions==3) THEN
            WRITE(intgString2,'(I0)') numberOfDimensions
            WRITE(myWorldComputationNodeNumber,*) &
              & 'Shape.  Dimension=',TRIM(intgString2),', simplex(2;3)*simplex*simplex'
            IF(maxNodesPerElement==4) THEN
              WRITE(myWorldComputationNodeNumber,*) &
                & '#Scale factor sets= 1'
              WRITE(intgString,'(I0)') maxNodesPerElement
              WRITE(myWorldComputationNodeNumber,*) &
                & ' l.simplex(2;3)*l.simplex*l.simplex, #Scale factors= ', TRIM(intgString)
            ELSE IF (maxNodesPerElement== 10 ) THEN
              WRITE(myWorldComputationNodeNumber,*) '#Scale factor sets= 1'
              WRITE(intgString,'(I0)') maxNodesPerElement
              WRITE(myWorldComputationNodeNumber,*) &
                & ' q.simplex(2;3)*q.simplex*q.simplex, #Scale factors= ', TRIM(intgString)
            ELSE IF(maxNodesPerElement==20) THEN
              WRITE(myWorldComputationNodeNumber,*) '#Scale factor sets= 1'
              WRITE(intgString,'(I0)') maxNodesPerElement
              WRITE(myWorldComputationNodeNumber,*) &
                & ' q.simplex(2;3)*q.simplex*q.simplex, #Scale factors= ', TRIM(intgString)
            ENDIF
          ELSE
            WRITE(myWorldComputationNodeNumber,*) '#Scale factor sets= 0'
          ENDIF
        ENDIF

        WRITE(intgString,'(I0)') maxNodesPerElement
        WRITE(myWorldComputationNodeNumber,*) '#Nodes= ',TRIM(intgString)
        WRITE(intgString,'(I0)') numberOfOutputFields
        WRITE(myWorldComputationNodeNumber,*) '#Fields= ',TRIM(intgString)
        numberOfFieldComponents(1) = numberOfDimensions
        numberOfFieldComponents(2) = numberOfVariableComponents
        numberOfFieldComponents(3) = numberOfSourceComponents
        DO fieldIdx=1,numberOfOutputFields
          IF(fieldIdx==1)THEN
            WRITE(intgString,'(I0)') numberOfDimensions
            WRITE(myWorldComputationNodeNumber,*) &
              & ' 1) coordinates,  coordinate, rectangular cartesian, #Components= ',TRIM(intgString)
          ELSE IF(fieldIdx==2) THEN
            WRITE(intgString,'(I0)') numberOfVariableComponents
            WRITE(myWorldComputationNodeNumber,*) &
              & ' 2) dependent,  field,  rectangular cartesian, #Components= ',TRIM(intgString)
          ELSE IF(fieldIdx==3) THEN
            WRITE(intgString,'(I0)') numberOfSourceComponents
            WRITE(myWorldComputationNodeNumber,*) &
              & ' 3) source,  field,  rectangular cartesian, #Components= ',TRIM(intgString)
          ENDIF

          DO componentIdx=1,numberOfFieldComponents(fieldIdx)
            IF(numberOfDimensions==1) THEN
              IF(fieldIdx==1)THEN
                IF(componentIdx==1) THEN
                  WRITE(myWorldComputationNodeNumber,*)'   x.   l.Lagrange, no modify, standard node based.'
                ELSE IF(componentIdx==2) THEN
                  WRITE(myWorldComputationNodeNumber,*)'   y.   l.Lagrange, no modify, standard node based.'
                ELSE IF(componentIdx==3) THEN
                  WRITE(myWorldComputationNodeNumber,*)'   z.   l.Lagrange, no modify, standard node based.'
                ENDIF
              ELSE
                WRITE(myWorldComputationNodeNumber,*) &
                  & '   ',componentIdx,'.   l.Lagrange, no modify, standard node based.'
              ENDIF
            ELSE IF(numberOfDimensions==2) THEN
              IF(fieldIdx==1)THEN
                IF(componentIdx==1) THEN
                  IF(maxNodesPerElement==4)THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   x.   l.Lagrange*l.Lagrange, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==9) THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   x.   q.Lagrange*q.Lagrange, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==16)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   x.   c.Lagrange*c.Lagrange, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==3)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   x.  l.simplex(2)*l.simplex, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==6)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   x.  q.simplex(2)*q.simplex, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==10)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   x.  c.simplex(2)*c.simplex, no modify, standard node based.'
                  ENDIF
                ELSE IF(componentIdx==2) THEN
                  IF(maxNodesPerElement==4) THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   y.   l.Lagrange*l.Lagrange, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==9)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   y.   q.Lagrange*q.Lagrange, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==16)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   y.   c.Lagrange*c.Lagrange, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==3)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   y.  l.simplex(2)*l.simplex, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==6)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   y.  q.simplex(2)*q.simplex, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==10)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   y.  c.simplex(2)*c.simplex, no modify, standard node based.'
                  ENDIF
                ELSE IF(componentIdx==3) THEN
                  IF(maxNodesPerElement==4) THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   z.   l.Lagrange*l.Lagrange, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==9)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   z.   q.Lagrange*q.Lagrange, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==16)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   z.   c.Lagrange*c.Lagrange, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==3)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   z.  l.simplex(2)*l.simplex, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==6)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   z.  q.simplex(2)*q.simplex, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==10)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   z.  c.simplex(2)*c.simplex, no modify, standard node based.'
                  ENDIF
                ENDIF
              ELSE
                IF(maxNodesPerElement==4) THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',componentIdx,'.   l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(maxNodesPerElement==9)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',componentIdx,'.   q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(maxNodesPerElement==16)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',componentIdx,'.   c.Lagrange*c.Lagrange, no modify, standard node based.'
                ELSE IF(maxNodesPerElement==3)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',componentIdx,'.  l.simplex(2)*l.simplex, no modify, standard node based.'
                ELSE IF(maxNodesPerElement==6)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',componentIdx,'.  q.simplex(2)*q.simplex, no modify, standard node based.'
                ELSE IF(maxNodesPerElement==10)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',componentIdx,'.  c.simplex(2)*c.simplex, no modify, standard node based.'
                ENDIF
              ENDIF
            ELSE IF(numberOfDimensions==3) THEN
              IF(fieldIdx==1)THEN
                IF(componentIdx==1) THEN
                  IF(maxNodesPerElement==8) THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   x.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==27)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   x.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==64)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   x.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==4)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   x.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==10)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   x.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==20)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   x.  c.simplex(2;3)*c.simplex*c.simplex, no modify, standard node based.'
                  ENDIF
                ELSE IF(componentIdx==2) THEN
                  IF(maxNodesPerElement==8) THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   y.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==27)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   y.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==64)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   y.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==4)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   y.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==10)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   y.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==20)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   y.  c.simplex(2;3)*c.simplex*c.simplex, no modify, standard node based.'
                  ENDIF
                ELSE IF(componentIdx==3) THEN
                  IF(maxNodesPerElement==8) THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   z.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==27)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   z.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==64)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   z.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==4)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   z.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==10)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   z.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.'
                  ELSE IF(maxNodesPerElement==20)  THEN
                    WRITE(myWorldComputationNodeNumber,*) &
                      & '   z.  c.simplex(2;3)*c.simplex*c.simplex, no modify, standard node based.'
                  ENDIF
                ENDIF
              ELSE
                IF(maxNodesPerElement==8) THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',componentIdx,'.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
                ELSE IF(maxNodesPerElement==27)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',componentIdx,'.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
                ELSE IF(maxNodesPerElement==64)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',componentIdx,'.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'
                ELSE IF(maxNodesPerElement==4)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',componentIdx,'.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.'
                ELSE IF(maxNodesPerElement==10)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',componentIdx,'.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.'
                ELSE IF(maxNodesPerElement==20)  THEN
                  WRITE(myWorldComputationNodeNumber,*) &
                    & '   ',componentIdx,'.  c.simplex(2;3)*c.simplex*c.simplex, no modify, standard node based.'
                ENDIF
              ENDIF
            ENDIF
            WRITE(intgString,'(I0)') maxNodesPerElement
            WRITE(myWorldComputationNodeNumber,*) '   #Nodes= ',TRIM(intgString)

            DO nodeIdx=1,maxNodesPerElement
              WRITE(intgString,'(I0)') nodeIdx
              WRITE(myWorldComputationNodeNumber,*) '    ',TRIM(intgString),'.  #Values=1'
              WRITE(myWorldComputationNodeNumber,*) '     Value indices:     1'
              WRITE(myWorldComputationNodeNumber,*) '     Scale factor indices:   ',TRIM(intgString)
            ENDDO !nodeIdx
          ENDDO !componentIdx
        ENDDO !fieldIdx

        IF(.NOT.ALLOCATED(elementNodes)) ALLOCATE(elementNodes(numberOfElements,maxNodesPerElement))
        IF(.NOT.ALLOCATED(elementNodesScales)) ALLOCATE(elementNodesScales(numberOfElements,maxNodesPerElement))

        DO elementIdx=1,numberOfElements
          DO nodeIdx=1,maxNodesPerElement
            CALL DomainElements_ElementNodeGet(domainElements,nodeIdx,elementIdx,localNodeNumber,err,error,*999)
            CALL DomainMapping_LocalToGlobalGet(nodesMapping,localNodeNumber,globalNodeNumber,err,error,*999)
            elementNodes(elementIdx,nodeIdx)=globalNodeNumber
            elementNodesScales(elementIdx,nodeIdx)=1.0000000000000000E+00
          ENDDO !nodeIdx
        ENDDO !elementIdx

        DO elementIdx=1,numberOfElements
          CALL DecompositionElements_ElementGlobalNumberGet(decompositionelements,elementIdx,globalElementNumber,err,error,*999)
          IF (basisShapeType==BASIS_LAGRANGE_HERMITE_TP_TYPE) THEN
            WRITE(intgString,'(I0)') globalElementNumber
            WRITE(myWorldComputationNodeNumber,*) 'Element:     ', TRIM(intgString),' 0  0'
            WRITE(myWorldComputationNodeNumber,*) '   Nodes:'
            WRITE(myWorldComputationNodeNumber,*) '   ', elementNodes(elementIdx,1:maxNodesPerElement)
            WRITE(myWorldComputationNodeNumber,*) '   Scale factors:'
            WRITE(myWorldComputationNodeNumber,*) '   ',elementNodesScales(elementIdx,1:maxNodesPerElement)
          ELSEIF(basisShapeType==BASIS_SIMPLEX_TYPE) THEN
            IF(.NOT.ALLOCATED(simplexOutputHelp)) ALLOCATE(simplexOutputHelp(maxNodesPerElement))
            IF(numberOfDimensions==2)THEN
              simplexOutputHelp(1)=elementNodes(elementIdx,1)
              simplexOutputHelp(2)=elementNodes(elementIdx,2)
              simplexOutputHelp(3)=elementNodes(elementIdx,3)
            ELSE IF(numberOfDimensions==3) THEN
              simplexOutputHelp(1)=elementNodes(elementIdx,1)
              simplexOutputHelp(2)=elementNodes(elementIdx,4)
              simplexOutputHelp(3)=elementNodes(elementIdx,2)
              simplexOutputHelp(4)=elementNodes(elementIdx,3)
            ENDIF
            WRITE(intgString,'(I0)') globalElementNumber
            WRITE(myWorldComputationNodeNumber,*) 'Element:     ',TRIM(intgString),' 0  0'
            WRITE(myWorldComputationNodeNumber,*) '   Nodes:'
            WRITE(myWorldComputationNodeNumber,*) '   ',simplexOutputHelp
            WRITE(myWorldComputationNodeNumber,*) '   Scale factors:'
            WRITE(myWorldComputationNodeNumber,*) '   ',elementNodesScales(elementIdx,1:maxNodesPerElement)
          ENDIF
        ENDDO !elementIdx
        CLOSE(myWorldComputationNodeNumber)
      ENDIF !exportExelem flag check
      IF(ALLOCATED(elementNodes)) DEALLOCATE(elementNodes)
      IF(ALLOCATED(elementNodesScales)) DEALLOCATE(elementNodesScales)
      IF(ALLOCATED(simplexOutputHelp)) DEALLOCATE(simplexOutputHelp)
      IF(outputSource) THEN
        IF(ASSOCIATED(sourceParameterData)) &
          & CALL FieldVariable_ParameterSetDataRestore(sourceVariable,FIELD_VALUES_SET_TYPE,sourceParameterData,err,error,*999)
      ENDIF
      IF(ASSOCIATED(dependentParameterData)) &
        & CALL FieldVariable_ParameterSetDataRestore(dependentVariable,FIELD_VALUES_SET_TYPE,dependentParameterData,err,error,*999)
      IF(ASSOCIATED(geometricParameterData)) &
        & CALL FieldVariable_ParameterSetDataRestore(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameterData,err,error,*999)
    ENDIF !Reaction-diffusion equations set

    EXITS("ReactionDiffusion_IOWriteCMGUI")
    RETURN
999 IF(ALLOCATED(elementNodes)) DEALLOCATE(elementNodes)
    IF(ALLOCATED(elementNodesScales)) DEALLOCATE(elementNodesScales)
    IF(ALLOCATED(simplexOutputHelp)) DEALLOCATE(simplexOutputHelp)
    IF(ASSOCIATED(sourceParameterData)) &
      & CALL FieldVariable_ParameterSetDataRestore(sourceVariable,FIELD_VALUES_SET_TYPE,sourceParameterData,err,error,*998)
998 IF(ASSOCIATED(dependentParameterData)) &
      & CALL FieldVariable_ParameterSetDataRestore(dependentVariable,FIELD_VALUES_SET_TYPE,dependentParameterData,err,error,*997)
997 IF(ASSOCIATED(geometricParameterData)) &
      & CALL FieldVariable_ParameterSetDataRestore(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameterData,err,error,*996)
996 ERRORSEXITS("ReactionDiffusion_IOWriteCMGUI",err,error)
    RETURN 1

  END SUBROUTINE ReactionDiffusion_IOWriteCMGUI

  !
  !===============================================================================================================================  
  !


END MODULE ReactionDiffusionIORoutines
