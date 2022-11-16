!> \file
!> \author Ting Yu
!> \brief This module handles all analytic analysis routines.
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
!> Contributor(s): Ting Yu, Chris Bradley
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

!> \defgroup OpenCMISS_AnalyticAnalysis OpenCMISS::Iron::AnalyticAnalysis
!> This module handles all analytic analysis routines.
MODULE AnalyticAnalysisRoutines

  USE BasisRoutines
  USE BasisAccessRoutines
  USE CmissMPI
  USE ComputationRoutines
  USE ComputationAccessRoutines
  USE Constants
  USE DecompositionAccessRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE MatrixVector
  USE MeshAccessRoutines
#ifndef NOMPIMOD
  USE MPI
#endif
  USE Strings
  USE Timer
  USE Types

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  !Module parameters

  !> \addtogroup AnalyticAnalysis_Constants AnalyticAnalysis::Constants
  !>@{
  !> \addtogroup AnalyticAnalysisRoutines_ErrorTypes AnalyticAnalysis::Constants::ErrorTypes
  !> \brief errors definition type parameters
  !> \see AnalyticAnalysisRoutines,OPENCMISS_ErrorTypes
  !>@{
  INTEGER(INTG), PARAMETER :: ANALYTIC_ABSOLUTE_ERROR_TYPE=1 !<The absolute type \see AnalyticAnalysisRoutines_ErrorTypes,AnalyticAnalysisRoutines
  INTEGER(INTG), PARAMETER :: ANALYTIC_PERCENTAGE_ERROR_TYPE=2 !<The percentage type \see AnalyticAnalysisRoutines_ErrorTypes,AnalyticAnalysisRoutines
  INTEGER(INTG), PARAMETER :: ANALYTIC_RELATIVE_ERROR_TYPE=3 !<The relative type \see AnalyticAnalysisRoutines_ErrorTypes,AnalyticAnalysisRoutines
  !>@}
  !>@}

  !Module types

  !Module variables

  !Interfaces

  PUBLIC ANALYTIC_ABSOLUTE_ERROR_TYPE,ANALYTIC_PERCENTAGE_ERROR_TYPE,ANALYTIC_RELATIVE_ERROR_TYPE

  PUBLIC AnalyticAnalysis_Output
  
  PUBLIC AnalyticAnalysis_AbsoluteErrorGetNode,AnalyticAnalysis_PercentageErrorGetNode, &
    & AnalyticAnalysis_RelativeErrorGetNode,AnalyticAnalysis_RMSErrorGetNode

  PUBLIC AnalyticAnalysis_AbsoluteErrorGetElement,AnalyticAnalysis_PercentageErrorGetElement, &
    & AnalyticAnalysis_RelativeErrorGetElement,AnalyticAnalysis_RMSErrorGetElement

  PUBLIC AnalyticAnalysis_AbsoluteErrorGetConstant,AnalyticAnalysis_PercentageErrorGetConstant, &
    & AnalyticAnalysis_RelativeErrorGetConstant
   
  PUBLIC AnalyticAnalysis_IntegralNumericalValueGet,AnalyticAnalysis_IntegralAnalyticValueGet, &
    & AnalyticAnalysis_IntegralPercentageErrorGet,AnalyticAnalysis_IntegralAbsoluteErrorGet, &
    & AnalyticAnalysis_IntegralRelativeErrorGet,AnalyticAnalysis_IntegralNIDNumericalValueGet, &
    & AnalyticAnalysis_IntegralNIDErrorGet
    
CONTAINS  

  !
  !=================================================================================================================================
  !  

  !>Output the analytic error analysis for a dependent field compared to the analytic values parameter set. \see OpenCMISS::Iron::cmfe_AnalyticAnalytisOutput
  SUBROUTINE AnalyticAnalysis_Output(field,filename,err,error,*)
  
    !Argument variables 
    TYPE(FieldType), INTENT(IN), POINTER :: field !<A pointer to the dependent field to calculate the analytic error analysis for
    CHARACTER(LEN=*) :: filename !<If not empty, the filename to output the analytic analysis to. If empty, the analysis will be output to the standard output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,derivativeIdx,elementIdx,ghostNumber(8),groupCommunicator,interpolationType,localDOF, &
      & mpiIerror,myGroupComputationNodeNumber,nodeIdx,number(MAXIMUM_GLOBAL_DERIV_NUMBER),numberOfComponents, &
      & numberOfGroupComputationNodes,numberOfDerivatives,numberOfElements,numberOfNodes,numberOfVariables,outputID, &
      & totalNumberOfElements,totalNumberOfNodes,variableIdx,variableType
    REAL(DP) :: ghostRMSErrorPer(MAXIMUM_GLOBAL_DERIV_NUMBER),ghostRMSErrorAbs(MAXIMUM_GLOBAL_DERIV_NUMBER), &
      & ghostRMSErrorRel(MAXIMUM_GLOBAL_DERIV_NUMBER),rmsErrorPer(MAXIMUM_GLOBAL_DERIV_NUMBER), &
      & rmsErrorAbs(MAXIMUM_GLOBAL_DERIV_NUMBER),rmsErrorRel(MAXIMUM_GLOBAL_DERIV_NUMBER),values(5)
    REAL(DP), POINTER :: analyticValues(:),numericalValues(:)
    REAL(DP), ALLOCATABLE :: integralErrors(:,:),ghostIntegralErrors(:,:)
    CHARACTER(LEN=40) :: firstFormat
    CHARACTER(LEN=MAXSTRLEN) :: fName
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: componentLabel,localError,localString
    TYPE(WorkGroupType), POINTER :: workGroup

    NULLIFY(analyticValues)
    NULLIFY(numericalValues)
    
    ENTERS("AnalyticAnalysis_Output",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    CALL Field_AssertIsDependent(field,err,error,*999)
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(field,decomposition,err,error,*999)
    NULLIFY(workGroup)
    CALL Decomposition_WorkGroupGet(decomposition,workGroup,err,error,*999)
    CALL WorkGroup_GroupCommunicatorGet(workGroup,groupCommunicator,err,error,*999)
    CALL WorkGroup_NumberOfGroupNodesGet(workGroup,numberOfGroupComputationNodes,err,error,*999)
    CALL WorkGroup_GroupNodeNumberGet(workGroup,myGroupComputationNodeNumber,err,error,*999)
    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    IF(LEN_TRIM(filename)>=1) THEN
!!TODO \todo have more general ascii file mechanism
      IF(numberOfGroupComputationNodes>1) THEN
        WRITE(fName,'(A,".opanal.",I0)') filename(1:LEN_TRIM(filename)),myGroupComputationNodeNumber
      ELSE
        fName=filename(1:LEN_TRIM(filename))//".opanal"
      ENDIF
      outputID=IO1_FILE_UNIT
      OPEN(UNIT=outputID,FILE=fName(1:LEN_TRIM(fName)),STATUS="REPLACE",FORM="FORMATTED",IOSTAT=err)
      IF(err/=0) CALL FlagError("Error opening analysis output file.",err,error,*999)            
    ELSE
      outputID=GENERAL_OUTPUT_TYPE
    ENDIF
    CALL WriteString(outputID,"Analytic error analysis:",err,error,*999)
    CALL WriteString(outputID,"",err,error,*999)
    localString="Field "//TRIM(NumberToVString(FIELD%userNumber,"*",err,error))//" : "//FIELD%LABEL
    IF(err/=0) GOTO 999
    CALL WriteString(outputID,localString,err,error,*999)
    !Loop over the variables
    CALL Field_NumberOfVariablesGet(field,numberOfVariables,err,error,*999)
    DO variableIdx=1,numberOfVariables
      NULLIFY(fieldVariable)
      CALL Field_VariableIndexGet(field,variableIdx,fieldVariable,variableType,err,error,*999)
      CALL WriteString(outputID,"",err,error,*999)
      localString="Variable "//TRIM(NumberToVString(variableType,"*",err,error))//" : "//fieldVariable%variableLabel
      IF(err/=0) GOTO 999
      CALL WriteString(outputID,localString,err,error,*999)
      CALL WriteString(outputID,"",err,error,*999)
      !Get the dependent and analytic parameter sets
      NULLIFY(numericalValues)
      CALL Field_ParameterSetDataGet(field,variableType,FIELD_VALUES_SET_TYPE,numericalValues,err,error,*999)
      NULLIFY(analyticValues)
      CALL Field_ParameterSetDataGet(field,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,analyticValues,err,error,*999)
      !Loop over the components
      CALL FieldVariable_NumberOfComponentsGet(fieldVariable,numberOfComponents,err,error,*999)
      DO componentIdx=1,numberOfComponents
        NULLIFY(domain)
        CALL FieldVariable_ComponentDomainGet(fieldVariable,componentIdx,domain,err,error,*999)
        NULLIFY(domainTopology)
        CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
        CALL FieldVariable_ComponentLabelGet(fieldVariable,componentIdx,componentLabel,err,error,*999)
        localString="Component "//TRIM(NumberToVString(componentIdx,"*",err,error))//" : "//componentLabel
        IF(err/=0) GOTO 999
        CALL WriteString(outputID,localString,err,error,*999)
        CALL WriteString(outputID,"",err,error,*999)
        CALL FieldVariable_ComponentInterpolationGet(fieldVariable,componentIdx,interpolationType,err,error,*999)
        SELECT CASE(interpolationType)
        CASE(FIELD_CONSTANT_INTERPOLATION)
          CALL WriteString(outputID,"Constant errors:",err,error,*999)
          localString="                       Numerical      Analytic       % error  Absolute err  Relative err"
          CALL WriteString(outputID,localString,err,error,*999)
          CALL FieldVariable_ConstantDOFGet(fieldVariable,componentIdx,localDOF,err,error,*999)
          values(1)=numericalValues(localDOF)
          values(2)=analyticValues(localDOF)
          values(3)=AnalyticAnalysis_PercentageError(values(1),values(2))
          values(4)=AnalyticAnalysis_AbsoluteError(values(1),values(2))
          values(5)=AnalyticAnalysis_RelativeError(values(1),values(2))
          CALL WriteStringVector(outputID,1,1,5,5,5,values,"(20X,5(2X,E12.5))","(20X,5(2X,E12.5))",err,error,*999)
        CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
          NULLIFY(domainElements)
          CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
          NULLIFY(decompositionElements)
          CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
          number=0
          rmsErrorPer=0.0_DP
          rmsErrorAbs=0.0_DP
          rmsErrorRel=0.0_DP
          ghostNumber=0
          ghostRMSErrorPer=0.0_DP
          ghostRMSErrorAbs=0.0_DP
          ghostRMSErrorRel=0.0_DP
          CALL WriteString(outputID,"Element errors:",err,error,*999)
          localString="  Element#             Numerical      Analytic       % error  Absolute err  Relative err"
          CALL WriteString(outputID,localString,err,error,*999)
          CALL DomainElements_NumberOfElementsGet(domainElements,numberOfElements,err,error,*999)
          CALL DomainElements_TotalNumberOfElementsGet(domainElements,totalNumberOfElements,err,error,*999)
          DO elementIdx=1,numberOfElements
            CALL FieldVariable_LocalElementDOFGet(fieldVariable,elementIdx,componentIdx,localDOF,err,error,*999)
            values(1)=numericalValues(localDOF)
            values(2)=analyticValues(localDOF)
            values(3)=AnalyticAnalysis_PercentageError(values(1),values(2))
            values(4)=AnalyticAnalysis_AbsoluteError(values(1),values(2))
            values(5)=AnalyticAnalysis_RelativeError(values(1),values(2))
            !Accumlate the RMS errors
            number(1)=number(1)+1
            rmsErrorPer(1)=rmsErrorPer(1)+values(3)*values(3)
            rmsErrorAbs(1)=rmsErrorAbs(1)+values(4)*values(4)
            rmsErrorRel(1)=rmsErrorRel(1)+values(5)*values(5)
            WRITE(firstFormat,"(A,I10,A)") "('",decompositionElements%elements(elementIdx)%userNumber,"',20X,3(2X,E12.5))"
            CALL WriteStringVector(outputID,1,1,5,5,5,VALUES,firstFormat,"(20X,5(2X,E12.5))",err,error,*999)
          ENDDO !elementIdx
          DO elementIdx=numberOfElements+1,totalNumberOfElements
            CALL FieldVariable_LocalElementDOFGet(fieldVariable,elementIdx,componentIdx,localDOF,err,error,*999)
            values(1)=numericalValues(localDOF)
            values(2)=analyticValues(localDOF)
            values(3)=AnalyticAnalysis_PercentageError(values(1),values(2))
            values(4)=AnalyticAnalysis_AbsoluteError(values(1),values(2))
            values(5)=AnalyticAnalysis_RelativeError(values(1),values(2))
            !Accumlate the RMS errors
            ghostNumber(1)=ghostNumber(1)+1
            ghostRMSErrorPer(1)=ghostRMSErrorPer(1)+values(3)*values(3)
            ghostRMSErrorAbs(1)=ghostRMSErrorAbs(1)+values(4)*values(4)
            ghostRMSErrorRel(1)=ghostRMSErrorRel(1)+values(5)*values(5)
            WRITE(firstFormat,"(A,I10,A)") "('",decompositionElements%elements(elementIdx)%userNumber,"',20X,3(2X,E12.5))"
            CALL WriteStringVector(outputID,1,1,5,5,5,VALUES,firstFormat,"(20X,5(2X,E12.5))",err,error,*999)
          ENDDO !elementIdx
          !Output RMS errors                  
          CALL WriteString(outputID,"",err,error,*999)
          IF(number(1)>0) THEN
            IF(numberOfGroupComputationNodes>1) THEN
              !Local elements only
              CALL WriteString(outputID,"Local RMS errors:",err,error,*999)
              localString="                                                     % error  Absolute err  Relative err"
              CALL WriteString(outputID,localString,err,error,*999)
              values(1)=SQRT(rmsErrorPer(derivativeIdx)/number(derivativeIdx))
              values(2)=SQRT(rmsErrorAbs(derivativeIdx)/number(derivativeIdx))
              values(3)=SQRT(rmsErrorRel(derivativeIdx)/number(derivativeIdx))
              CALL WriteStringVector(outputID,1,1,3,3,3,VALUES,"(46X,3(2X,E12.5))","(46X,3(2X,E12.5))",err,error,*999)
              !Local and ghost nodes
              CALL WriteString(outputID,"Local + Ghost RMS errors:",err,error,*999)
              localString="                                                     % error  Absolute err  Relative err"
              CALL WriteString(outputID,localString,err,error,*999)
              values(1)=SQRT((rmsErrorPer(1)+ghostRMSErrorPer(1))/(number(1)+ghostNumber(1)))
              values(2)=SQRT((rmsErrorAbs(1)+ghostRMSErrorAbs(1))/(number(1)+ghostNumber(1)))
              values(3)=SQRT((rmsErrorRel(1)+ghostRMSErrorRel(1))/(number(1)+ghostNumber(1)))
              CALL WriteStringVector(outputID,1,1,3,3,3,VALUES,"(46X,3(2X,E12.5))","(46X,3(2X,E12.5))",err,error,*999)
              !Global RMS values
              !Collect the values across the ranks
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,NUMBER,1,MPI_INTEGER,MPI_SUM,groupCommunicator,mpiIerror)
              CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIerror,err,error,*999)
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,rmsErrorPer,1,MPI_DOUBLE_PRECISION,MPI_SUM,groupCommunicator,mpiIerror)
              CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIerror,err,error,*999)
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,rmsErrorAbs,1,MPI_DOUBLE_PRECISION,MPI_SUM,groupCommunicator,mpiIerror)
              CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIerror,err,error,*999)
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,rmsErrorRel,1,MPI_DOUBLE_PRECISION,MPI_SUM,groupCommunicator,mpiIerror)
              CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIerror,err,error,*999)
              CALL WriteString(outputID,"Global RMS errors:",err,error,*999)
              localString="                                                     % error  Absolute err  Relative err"
              CALL WriteString(outputID,localString,err,error,*999)
              values(1)=SQRT(rmsErrorPer(1)/number(1))
              values(2)=SQRT(rmsErrorAbs(1)/number(1))
              values(3)=SQRT(rmsErrorRel(1)/number(1))
              CALL WriteStringVector(outputID,1,1,3,3,3,VALUES,"(46X,3(2X,E12.5))","(46X,3(2X,E12.5))",err,error,*999)
            ELSE
              CALL WriteString(outputID,"RMS errors:",err,error,*999)
              localString="                                                     % error  Absolute err  Relative err"
              CALL WriteString(outputID,localString,err,error,*999)
              values(1)=SQRT(rmsErrorPer(derivativeIdx)/number(derivativeIdx))
              values(2)=SQRT(rmsErrorAbs(derivativeIdx)/number(derivativeIdx))
              values(3)=SQRT(rmsErrorRel(derivativeIdx)/number(derivativeIdx))
              CALL WriteStringVector(outputID,1,1,3,3,3,VALUES,"(46X,3(2X,E12.5))","(46X,3(2X,E12.5))",err,error,*999)
            ENDIF
          ENDIF
        CASE(FIELD_NODE_BASED_INTERPOLATION)
          NULLIFY(domainNodes)
          CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
          number=0
          rmsErrorPer=0.0_DP
          rmsErrorAbs=0.0_DP
          rmsErrorRel=0.0_DP
          ghostNumber=0
          ghostRMSErrorPer=0.0_DP
          ghostRMSErrorAbs=0.0_DP
          ghostRMSErrorRel=0.0_DP
          CALL WriteString(outputID,"Nodal errors:",err,error,*999)
          localString="     Node#  Deriv#     Numerical      Analytic       % error  Absolute err  Relative err"
          CALL WriteString(outputID,localString,err,error,*999)
          CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
          CALL DomainNodes_TotalNumberOfNodesGet(domainNodes,totalNumberOfNodes,err,error,*999)
          DO nodeIdx=1,numberOfNodes
            CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfDerivatives,err,error,*999)
            DO derivativeIdx=1,numberOfDerivatives
              CALL FieldVariable_LocalNodeDOFGet(fieldVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOF,err,error,*999)
              values(1)=numericalValues(localDOF)
              values(2)=analyticValues(localDOF)
              values(3)=AnalyticAnalysis_PercentageError(values(1),values(2))
              values(4)=AnalyticAnalysis_AbsoluteError(values(1),values(2))
              values(5)=AnalyticAnalysis_RelativeError(values(1),values(2))
              !Accumlate the RMS errors
              number(derivativeIdx)=number(derivativeIdx)+1
              rmsErrorPer(derivativeIdx)=rmsErrorPer(derivativeIdx)+values(3)*values(3)
              rmsErrorAbs(derivativeIdx)=rmsErrorAbs(derivativeIdx)+values(4)*values(4)
              rmsErrorRel(derivativeIdx)=rmsErrorRel(derivativeIdx)+values(5)*values(5)
              IF(derivativeIdx==1) THEN
                WRITE(firstFormat,"(A,I10,A,I6,A)") "('",domainNodes%nodes(nodeIdx)%userNumber,"',2X,'", &
                  & derivativeIdx,"',5(2X,E12.5))"
              ELSE
                WRITE(firstFormat,"(A,I6,A)") "(12X,'",derivativeIdx,"',5(2X,E12.5))"
              ENDIF
              CALL WriteStringVector(outputID,1,1,5,5,5,VALUES,firstFormat,"(20X,5(2X,E12.5))",err,error,*999)
            ENDDO !derivativeIdx
          ENDDO !nodeIdx
          DO nodeIdx=numberOfNodes+1,totalNumberOfNodes
            CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfDerivatives,err,error,*999)
            DO derivativeIdx=1,numberOfDerivatives
              CALL FieldVariable_LocalNodeDOFGet(fieldVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOF,err,error,*999)
              values(1)=numericalValues(localDOF)
              values(2)=analyticValues(localDOF)
              values(3)=AnalyticAnalysis_PercentageError(values(1),values(2))
              values(4)=AnalyticAnalysis_AbsoluteError(values(1),values(2))
              values(5)=AnalyticAnalysis_RelativeError(values(1),values(2))
              !Accumlate the RMS errors
              ghostNumber(derivativeIdx)=ghostNumber(derivativeIdx)+1
              ghostRMSErrorPer(derivativeIdx)=ghostRMSErrorPer(derivativeIdx)+values(3)*values(3)
              ghostRMSErrorAbs(derivativeIdx)=ghostRMSErrorAbs(derivativeIdx)+values(4)*values(4)
              ghostRMSErrorRel(derivativeIdx)=ghostRMSErrorRel(derivativeIdx)+values(5)*values(5)
              IF(derivativeIdx==1) THEN
                WRITE(firstFormat,"(A,I10,A,I6,A)") "('",domainNodes%nodes(nodeIdx)%userNumber, &
                  & "',2X,'",derivativeIdx,"',5(2X,E12.5))"
              ELSE
                WRITE(firstFormat,"(A,I6,A)") "(12X,'",derivativeIdx,"',5(2X,E12.5))"
              ENDIF
              CALL WriteStringVector(outputID,1,1,5,5,5,VALUES,firstFormat,"(20X,5(2X,E12.5))",err,error,*999)
            ENDDO !derivativeIdx
          ENDDO !nodeIdx
          !Output RMS errors                  
          CALL WriteString(outputID,"",err,error,*999)
          IF(numberOfGroupComputationNodes>1) THEN
            IF(ANY(number>0)) THEN
              !Local nodes only
              CALL WriteString(outputID,"Local RMS errors:",err,error,*999)
              localString="            Deriv#                                   % error  Absolute err  Relative err"
              CALL WriteString(outputID,localString,err,error,*999)
              DO derivativeIdx=1,MAXIMUM_GLOBAL_DERIV_NUMBER
                IF(number(derivativeIdx)>0) THEN
                  values(1)=SQRT(rmsErrorPer(derivativeIdx)/number(derivativeIdx))
                  values(2)=SQRT(rmsErrorAbs(derivativeIdx)/number(derivativeIdx))
                  values(3)=SQRT(rmsErrorRel(derivativeIdx)/number(derivativeIdx))
                  WRITE(firstFormat,"(A,I6,A)") "(12X,'",derivativeIdx,"',28X,3(2X,E12.5))"
                  CALL WriteStringVector(outputID,1,1,3,3,3,VALUES,firstFormat,"(46X,3(2X,E12.5))",err,error,*999)
                ENDIF
              ENDDO !derivativeIdx
              !Local and ghost nodes
              CALL WriteString(outputID,"Local + Ghost RMS errors:",err,error,*999)
              localString="            Deriv#                                   % error  Absolute err  Relative err"
              CALL WriteString(outputID,localString,err,error,*999)
              DO derivativeIdx=1,MAXIMUM_GLOBAL_DERIV_NUMBER
                IF(number(derivativeIdx)>0) THEN
                  values(1)=SQRT((rmsErrorPer(derivativeIdx)+ghostRMSErrorPer(derivativeIdx))/ &
                    & (number(derivativeIdx)+ghostNumber(derivativeIdx)))
                  values(2)=SQRT((rmsErrorAbs(derivativeIdx)+ghostRMSErrorAbs(derivativeIdx))/ &
                    & (number(derivativeIdx)+ghostNumber(derivativeIdx)))
                  values(3)=SQRT((rmsErrorRel(derivativeIdx)+ghostRMSErrorRel(derivativeIdx))/ &
                    & (number(derivativeIdx)+ghostNumber(derivativeIdx)))
                  WRITE(firstFormat,"(A,I6,A)") "(12X,'",derivativeIdx,"',28X,3(2X,E12.5))"
                  CALL WriteStringVector(outputID,1,1,3,3,3,VALUES,firstFormat,"(46X,3(2X,E12.5))",err,error,*999)
                ENDIF
              ENDDO !derivativeIdx
              !Global RMS values
              !Collect the values across the ranks
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,NUMBER,MAXIMUM_GLOBAL_DERIV_NUMBER,MPI_INTEGER,MPI_SUM,groupCommunicator,mpiIerror)
              CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIerror,err,error,*999)
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,rmsErrorPer,MAXIMUM_GLOBAL_DERIV_NUMBER,MPI_DOUBLE_PRECISION,MPI_SUM, &
                & groupCommunicator,mpiIerror)
              CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIerror,err,error,*999)
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,rmsErrorAbs,MAXIMUM_GLOBAL_DERIV_NUMBER,MPI_DOUBLE_PRECISION,MPI_SUM, &
                & groupCommunicator,mpiIerror)
              CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIerror,err,error,*999)
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,rmsErrorRel,MAXIMUM_GLOBAL_DERIV_NUMBER,MPI_DOUBLE_PRECISION,MPI_SUM, &
                & groupCommunicator,mpiIerror)
              CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIerror,err,error,*999)
              CALL WriteString(outputID,"Global RMS errors:",err,error,*999)
              localString="            Deriv#                                   % error  Absolute err  Relative err"
              CALL WriteString(outputID,localString,err,error,*999)
              DO derivativeIdx=1,MAXIMUM_GLOBAL_DERIV_NUMBER
                IF(number(derivativeIdx)>0) THEN
                  values(1)=SQRT(rmsErrorPer(derivativeIdx)/number(derivativeIdx))
                  values(2)=SQRT(rmsErrorAbs(derivativeIdx)/number(derivativeIdx))
                  values(3)=SQRT(rmsErrorRel(derivativeIdx)/number(derivativeIdx))
                  WRITE(firstFormat,"(A,I6,A)") "(12X,'",derivativeIdx,"',28X,3(2X,E12.5))"
                  CALL WriteStringVector(outputID,1,1,3,3,3,VALUES,firstFormat,"(46X,3(2X,E12.5))",err,error,*999)
                ENDIF
              ENDDO !derivativeIdx
            ENDIF
          ELSE
            IF(ANY(number>0)) THEN
              CALL WriteString(outputID,"RMS errors:",err,error,*999)
              localString="            Deriv#                                   % error  Absolute err  Relative err"
              CALL WriteString(outputID,localString,err,error,*999)
              DO derivativeIdx=1,MAXIMUM_GLOBAL_DERIV_NUMBER
                IF(number(derivativeIdx)>0) THEN
                  values(1)=SQRT(rmsErrorPer(derivativeIdx)/number(derivativeIdx))
                  values(2)=SQRT(rmsErrorAbs(derivativeIdx)/number(derivativeIdx))
                  values(3)=SQRT(rmsErrorRel(derivativeIdx)/number(derivativeIdx))
                  WRITE(firstFormat,"(A,I6,A)") "(12X,'",derivativeIdx,"',28X,3(2X,E12.5))"
                  CALL WriteStringVector(outputID,1,1,3,3,3,VALUES,firstFormat,"(46X,3(2X,E12.5))",err,error,*999)
                ENDIF
              ENDDO !derivativeIdx
            ENDIF
          ENDIF
        CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The interpolation type of "// &
            & TRIM(NumberToVString(interpolationType,"*",err,error))//" for component number "// &
            & TRIM(NumberToVString(componentIdx,"*",err,error))//" of variable type "// &
            & TRIM(NumberToVString(variableType,"*",err,error))//" of field number "// &
            & TRIM(NumberToVString(field%userNumber,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        CALL WriteString(outputID,"",err,error,*999)
      ENDDO !componentIdx
      !Restore the dependent and analytic parameter sets
      CALL Field_ParameterSetDataRestore(FIELD,variableType,FIELD_VALUES_SET_TYPE,numericalValues,err,error,*999)
      CALL Field_ParameterSetDataRestore(FIELD,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,analyticValues,err,error,*999)
      !Allocated the integral errors
      ALLOCATE(integralErrors(6,numberOfComponents),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate integral errors.",err,error,*999)
      ALLOCATE(ghostIntegralErrors(6,numberOfComponents),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate ghost integral errors.",err,error,*999)
      CALL AnalyticAnalysis_IntegralErrors(fieldVariable,integralErrors,ghostIntegralErrors,err,error,*999)
      IF(numberOfGroupComputationNodes>1) THEN
        CALL WriteString(outputID,"Local Integral errors:",err,error,*999)
        localString="Component#             Numerical      Analytic       % error  Absolute err  Relative err"
        CALL WriteString(outputID,localString,err,error,*999)
        DO componentIdx=1,numberOfComponents
          values(1)=integralErrors(1,componentIdx)
          values(2)=integralErrors(3,componentIdx)
          values(3)=AnalyticAnalysis_PercentageError(values(1),values(2))
          values(4)=AnalyticAnalysis_AbsoluteError(values(1),values(2))
          values(5)=AnalyticAnalysis_RelativeError(values(1),values(2))
          WRITE(firstFormat,"(A,I10,A)") "('",componentIdx,"',2X,'Intg  ',5(2X,E12.5))"
          CALL WriteStringVector(outputID,1,1,5,5,5,VALUES,firstFormat,"(20X,5(2X,E12.5))",err,error,*999)
          values(1)=integralErrors(2,componentIdx)
          values(2)=integralErrors(4,componentIdx)
          values(3)=AnalyticAnalysis_PercentageError(values(1),values(2))
          values(4)=AnalyticAnalysis_AbsoluteError(values(1),values(2))
          values(5)=AnalyticAnalysis_RelativeError(values(1),values(2))
          WRITE(firstFormat,"(A)") "(12X,'Intg^2',5(2X,E12.5))"
          CALL WriteStringVector(outputID,1,1,5,5,5,VALUES,firstFormat,"(20X,5(2X,E12.5))",err,error,*999)
        ENDDO !componentIdx
        localString="Component#             Numerical      Analytic           NID        NID(%)"
        CALL WriteString(outputID,localString,err,error,*999)
        DO componentIdx=1,numberOfComponents
          values(1)=integralErrors(5,componentIdx)
          values(2)=integralErrors(3,componentIdx)
          values(3)=AnalyticAnalysis_NIDError(values(1),values(2))
          values(4)=values(3)*100.0_DP
          WRITE(firstFormat,"(A,I10,A)") "('",componentIdx,"',2X,'Diff  ',4(2X,E12.5))"
          CALL WriteStringVector(outputID,1,1,4,4,4,VALUES,firstFormat,"(20X,4(2X,E12.5))",err,error,*999)
          values(1)=integralErrors(6,componentIdx)
          values(2)=integralErrors(4,componentIdx)
          values(3)=AnalyticAnalysis_NIDError(values(1),values(2))
          values(4)=values(3)*100.0_DP
          WRITE(firstFormat,"(A)") "(12X,'Diff^2',4(2X,E12.5))"
          CALL WriteStringVector(outputID,1,1,4,4,4,VALUES,firstFormat,"(20X,4(2X,E12.5))",err,error,*999)
        ENDDO !componentIdx
        CALL WriteString(outputID,"Local + Ghost Integral errors:",err,error,*999)
        localString="Component#             Numerical      Analytic       % error  Absolute err  Relative err"
        CALL WriteString(outputID,localString,err,error,*999)
        DO componentIdx=1,numberOfComponents
          values(1)=integralErrors(1,componentIdx)+ghostIntegralErrors(1,componentIdx)
          values(2)=integralErrors(3,componentIdx)+ghostIntegralErrors(3,componentIdx)
          values(3)=AnalyticAnalysis_PercentageError(values(1),values(2))
          values(4)=AnalyticAnalysis_AbsoluteError(values(1),values(2))
          values(5)=AnalyticAnalysis_RelativeError(values(1),values(2))
          WRITE(firstFormat,"(A,I10,A)") "('",componentIdx,"',2X,'Intg  ',5(2X,E12.5))"
          CALL WriteStringVector(outputID,1,1,5,5,5,VALUES,firstFormat,"(20X,5(2X,E12.5))",err,error,*999)
          values(1)=integralErrors(2,componentIdx)
          values(2)=integralErrors(4,componentIdx)
          values(3)=AnalyticAnalysis_PercentageError(values(1),values(2))
          values(4)=AnalyticAnalysis_AbsoluteError(values(1),values(2))
          values(5)=AnalyticAnalysis_RelativeError(values(1),values(2))
          WRITE(firstFormat,"(A)") "(12X,'Intg^2',5(2X,E12.5))"
          CALL WriteStringVector(outputID,1,1,5,5,5,VALUES,firstFormat,"(20X,5(2X,E12.5))",err,error,*999)
        ENDDO !componentIdx
        localString="Component#             Numerical      Analytic           NID        NID(%)"
        CALL WriteString(outputID,localString,err,error,*999)
        DO componentIdx=1,numberOfComponents
          values(1)=integralErrors(5,componentIdx)+ghostIntegralErrors(5,componentIdx)
          values(2)=integralErrors(3,componentIdx)+ghostIntegralErrors(3,componentIdx)
          values(3)=AnalyticAnalysis_NIDError(values(1),values(2))
          values(4)=values(3)*100.0_DP
          WRITE(firstFormat,"(A,I10,A)") "('",componentIdx,"',2X,'Diff  ',4(2X,E12.5))"
          CALL WriteStringVector(outputID,1,1,4,4,4,VALUES,firstFormat,"(20X,4(2X,E12.5))",err,error,*999)
          values(1)=integralErrors(6,componentIdx)+ghostIntegralErrors(6,componentIdx)
          values(2)=integralErrors(4,componentIdx)+ghostIntegralErrors(4,componentIdx)
          values(3)=AnalyticAnalysis_NIDError(values(1),values(2))
          values(4)=values(3)*100.0_DP
          WRITE(firstFormat,"(A)") "(12X,'Diff^2',4(2X,E12.5))"
          CALL WriteStringVector(outputID,1,1,4,4,4,VALUES,firstFormat,"(20X,4(2X,E12.5))",err,error,*999)
        ENDDO !componentIdx
        !Collect the values across the ranks
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,integralErrors,6*numberOfComponents,MPI_DOUBLE_PRECISION, &
          & MPI_SUM,groupCommunicator,mpiIerror)
        CALL WriteString(outputID,"Global Integral errors:",err,error,*999)
        localString="Component#             Numerical      Analytic       % error  Absolute err  Relative err"
        CALL WriteString(outputID,localString,err,error,*999)
        DO componentIdx=1,numberOfComponents
          CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIerror,err,error,*999)
          values(1)=integralErrors(1,componentIdx)
          values(2)=integralErrors(3,componentIdx)
          values(4)=AnalyticAnalysis_AbsoluteError(values(1),values(2))
          values(5)=AnalyticAnalysis_RelativeError(values(1),values(2))
          WRITE(firstFormat,"(A,I10,A)") "('",componentIdx,"',2X,'Intg  ',5(2X,E12.5))"
          CALL WriteStringVector(outputID,1,1,5,5,5,VALUES,firstFormat,"(20X,5(2X,E12.5))",err,error,*999)
          values(1)=integralErrors(2,componentIdx)
          values(2)=integralErrors(4,componentIdx)
          values(3)=AnalyticAnalysis_PercentageError(values(1),values(2))
          values(4)=AnalyticAnalysis_AbsoluteError(values(1),values(2))
          values(5)=AnalyticAnalysis_RelativeError(values(1),values(2))
          WRITE(firstFormat,"(A)") "(12X,'Intg^2',5(2X,E12.5))"
          CALL WriteStringVector(outputID,1,1,5,5,5,VALUES,firstFormat,"(20X,5(2X,E12.5))",err,error,*999) 
        ENDDO !componentIdx
        localString="Component#             Numerical      Analytic           NID        NID(%)"
        CALL WriteString(outputID,localString,err,error,*999)
        DO componentIdx=1,numberOfComponents
          values(1)=integralErrors(5,componentIdx)
          values(2)=integralErrors(3,componentIdx)
          values(3)=AnalyticAnalysis_NIDError(values(1),values(2))
          values(4)=values(3)*100.0_DP
          WRITE(firstFormat,"(A,I10,A)") "('",componentIdx,"',2X,'Diff  ',4(2X,E12.5))"
          CALL WriteStringVector(outputID,1,1,4,4,4,VALUES,firstFormat,"(20X,4(2X,E12.5))",err,error,*999)
          values(1)=integralErrors(6,componentIdx)
          values(2)=integralErrors(4,componentIdx)
          values(3)=AnalyticAnalysis_NIDError(values(1),values(2))
          values(4)=values(3)*100.0_DP
          WRITE(firstFormat,"(A)") "(12X,'Diff^2',4(2X,E12.5))"
          CALL WriteStringVector(outputID,1,1,4,4,4,VALUES,firstFormat,"(20X,4(2X,E12.5))",err,error,*999)
        ENDDO !componentIdx
      ELSE
        CALL WriteString(outputID,"Integral errors:",err,error,*999)
        localString="Component#             Numerical      Analytic       % error  Absolute err  Relative err"
        CALL WriteString(outputID,localString,err,error,*999)
        DO componentIdx=1,numberOfComponents
          values(1)=integralErrors(1,componentIdx)
          values(2)=integralErrors(3,componentIdx)
          values(3)=AnalyticAnalysis_PercentageError(values(1),values(2))
          values(4)=AnalyticAnalysis_AbsoluteError(values(1),values(2))
          values(5)=AnalyticAnalysis_RelativeError(values(1),values(2))
          WRITE(firstFormat,"(A,I10,A)") "('",componentIdx,"',2X,'Intg  ',5(2X,E12.5))"
          CALL WriteStringVector(outputID,1,1,5,5,5,VALUES,firstFormat,"(20X,5(2X,E12.5))",err,error,*999)
          values(1)=integralErrors(2,componentIdx)
          values(2)=integralErrors(4,componentIdx)
          values(3)=AnalyticAnalysis_PercentageError(values(1),values(2))
          values(4)=AnalyticAnalysis_AbsoluteError(values(1),values(2))
          values(5)=AnalyticAnalysis_RelativeError(values(1),values(2))
          WRITE(firstFormat,"(A)") "(12X,'Intg^2',5(2X,E12.5))"
          CALL WriteStringVector(outputID,1,1,5,5,5,VALUES,firstFormat,"(20X,5(2X,E12.5))",err,error,*999)
        ENDDO !componentIdx
        localString="Component#             Numerical      Analytic           NID        NID(%)"
        CALL WriteString(outputID,localString,err,error,*999)
        DO componentIdx=1,numberOfComponents
          values(1)=integralErrors(5,componentIdx)
          values(2)=integralErrors(3,componentIdx)
          values(3)=AnalyticAnalysis_NIDError(values(1),values(2))
          values(4)=values(3)*100.0_DP
          WRITE(firstFormat,"(A,I10,A)") "('",componentIdx,"',2X,'Diff  ',4(2X,E12.5))"
          CALL WriteStringVector(outputID,1,1,4,4,4,VALUES,firstFormat,"(20X,4(2X,E12.5))",err,error,*999)
          values(1)=integralErrors(6,componentIdx)
          values(2)=integralErrors(4,componentIdx)
          values(3)=AnalyticAnalysis_NIDError(values(1),values(2))
          values(4)=values(3)*100.0_DP
          WRITE(firstFormat,"(A)") "(12X,'Diff^2',4(2X,E12.5))"
          CALL WriteStringVector(outputID,1,1,4,4,4,VALUES,firstFormat,"(20X,4(2X,E12.5))",err,error,*999)
        ENDDO !componentIdx
        CALL WriteString(outputID,"",err,error,*999)
      ENDIF
      IF(ALLOCATED(integralErrors)) DEALLOCATE(integralErrors)
      IF(ALLOCATED(ghostIntegralErrors)) DEALLOCATE(ghostIntegralErrors)
    ENDDO !variableIdx
    IF(LEN_TRIM(filename)>=1) THEN
      CLOSE(UNIT=outputID,IOSTAT=err)
      IF(err/=0) CALL FlagError("Error closing analysis output file.",err,error,*999)
    ENDIF
        
    EXITS("AnalyticAnalysis_Output")
    RETURN
999 IF(ALLOCATED(integralErrors)) DEALLOCATE(integralErrors)
    IF(ALLOCATED(ghostIntegralErrors)) DEALLOCATE(ghostIntegralErrors)
    ERRORSEXITS("AnalyticAnalysis_Output",err,error)
    RETURN 1
    
  END SUBROUTINE AnalyticAnalysis_Output

  !
  !================================================================================================================================
  !  

  !>Calculates the absolute error between a numeric value and an analytic value.
  PURE FUNCTION AnalyticAnalysis_AbsoluteError(numericValue,analyticValue)

    !Argument variables
    REAL(DP), INTENT(IN) :: numericValue !<The numerical value for the error calculation
    REAL(DP), INTENT(IN) :: analyticValue !<The analytic value for the error calculation
    !Function result
    REAL(DP) :: AnalyticAnalysis_AbsoluteError

    AnalyticAnalysis_AbsoluteError=ABS(analyticValue-numericValue)

  END FUNCTION AnalyticAnalysis_AbsoluteError
  
  !
  !================================================================================================================================
  !  

  !>Calculates the percentage error between a numeric value and an analytic value.
  PURE FUNCTION AnalyticAnalysis_PercentageError(numericValue,analyticValue)

    !Argument variables
    REAL(DP), INTENT(IN) :: numericValue !<The numerical value for the error calculation
    REAL(DP), INTENT(IN) :: analyticValue !<The analytic value for the error calculation
    !Function result
    REAL(DP) :: AnalyticAnalysis_PercentageError

    IF(ABS(analyticValue)>ZERO_TOLERANCE) THEN
      AnalyticAnalysis_PercentageError=(analyticValue-numericValue)/analyticValue*100.0_DP
    ELSE
      AnalyticAnalysis_PercentageError=0.0_DP
    ENDIF

  END FUNCTION AnalyticAnalysis_PercentageError
  
  !
  !================================================================================================================================
  !  

  !>Calculates the relative error between a numeric value and an analytic value.
  PURE FUNCTION AnalyticAnalysis_RelativeError(numericValue,analyticValue)

    !Argument variables
    REAL(DP), INTENT(IN) :: numericValue !<The numerical value for the error calculation
    REAL(DP), INTENT(IN) :: analyticValue !<The analytic value for the error calculation
    !Function result
    REAL(DP) :: AnalyticAnalysis_RelativeError

    IF(ABS(1.0_DP+analyticValue)>ZERO_TOLERANCE) THEN
      AnalyticAnalysis_RelativeError=(analyticValue-numericValue)/(1.0_DP+analyticValue)
    ELSE
      AnalyticAnalysis_RelativeError=0.0_DP
    ENDIF

  END FUNCTION AnalyticAnalysis_RelativeError
  
  !
  !================================================================================================================================
  !  

  !>Calculates the Normalised Integral Difference (NID) error with a numeric value and an analytic value.
  PURE FUNCTION AnalyticAnalysis_NIDError(numericValue,analyticValue)

    !Argument variables
    REAL(DP), INTENT(IN) :: numericValue !<The numerical value for the error calculation
    REAL(DP), INTENT(IN) :: analyticValue !<The analytic value for the error calculation
    !Function result
    REAL(DP) :: AnalyticAnalysis_NIDError

    IF(ABS(analyticValue)>ZERO_TOLERANCE) THEN
      AnalyticAnalysis_NIDError=(analyticValue-numericValue)/analyticValue
    ELSE
      AnalyticAnalysis_NIDError=(analyticValue-numericValue)/(1.0_DP+analyticValue)
    ENDIF

  END FUNCTION AnalyticAnalysis_NIDError
  
  !
  !================================================================================================================================
  !  

  !>Calculates the relative error between a numeric value and an analytic value.
  SUBROUTINE AnalyticAnalysis_IntegralErrors(fieldVariable,integralErrors,ghostIntegralErrors,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to calculate the integral errors for
    REAL(DP), INTENT(OUT) :: integralErrors(:,:) !<integralErrors(6,componentIdx). On exit, the integral errors for the local elements
    REAL(DP), INTENT(OUT) :: ghostIntegralErrors(:,:) !<ghostIntegralErrors(6,componentIdx). On exit, the integral errors for the ghost elements
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG), PARAMETER :: MAX_NUMBER_OF_COMPONENTS=99
    INTEGER(INTG) :: componentIdx,elementIdx,elementParameterIdx,gaussPointIdx,numberOfComponents, &
      & numberOfElementParameters(MAX_NUMBER_OF_COMPONENTS),numberOfElements,numberOfGauss,numberOfXi,scalingType, &
      & totalNumberOfElements,variableType
    REAL(DP) :: analyticIntegral,gaussWeight,jacobian,jacobianGaussWeight,numericalIntegral,phi
    TYPE(BasisType), POINTER :: basis,geometricBasis
    TYPE(BasisPtrType) :: dependentBasis(MAX_NUMBER_OF_COMPONENTS)
    TYPE(DecompositionType), POINTER :: dependentDecomposition,geometricDecomposition
    TYPE(DomainType), POINTER :: domain,dependentDomain,geometricDomain
    TYPE(DomainElementsType), POINTER :: domainElements,dependentDomainElements,geometricDomainElements
    TYPE(DomainTopologyType), POINTER :: domainTopology,dependentDomainTopology,geometricDomainTopology
    TYPE(FieldType), POINTER :: dependentField,geometricField
    TYPE(FieldInterpolatedPointType), POINTER :: geometricInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldInterpolationParametersType), POINTER :: analyticInterpParameters,geometricInterpParameters,numericalInterpParameters
    TYPE(FieldVariableType), POINTER :: geometricVariable
    TYPE(QuadratureSchemePtrType) :: quadratureScheme(MAX_NUMBER_OF_COMPONENTS)
    TYPE(VARYING_STRING) :: localError

    ENTERS("AnalyticAnalysis_IntegralErrors",err,error,*999)
    
    integralErrors=0.0_DP
    ghostIntegralErrors=0.0_DP
     
    CALL FieldVariable_NumberOfComponentsGet(fieldVariable,numberOfComponents,err,error,*999)
    
#ifdef WITH_PRECHECKS    
    IF(SIZE(integralErrors,1)<6.OR.SIZE(integralErrors,2)<numberOfComponents) THEN
      localError="Invalid size for integralErrors. The size is ("// &
        & TRIM(NumberToVString(SIZE(integralErrors,1),"*",err,error))//","// &
        & TRIM(NumberToVString(SIZE(integralErrors,2),"*",err,error))//") and it needs to be at least (6,"// &
        & TRIM(NumberToVString(numberOfComponents,"*",err,error))//")."
      CALL FlagError(localError,err,error,*999)
    ENDIF      
    IF(SIZE(ghostIntegralErrors,1)<6.OR.SIZE(ghostIntegralErrors,2)<numberOfComponents) THEN
      localError="Invalid size for ghostIntegralErrors. The size is ("// &
        & TRIM(NumberToVString(SIZE(ghostIntegralErrors,1),"*",err,error))//","// &
        & TRIM(NumberToVString(SIZE(ghostIntegralErrors,2),"*",err,error))//") and it needs to be at least (6,"// &
        & TRIM(NumberToVString(numberOfComponents,"*",err,error))//")."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    CALL FieldVariable_VariableTypeGet(fieldVariable,variableType,err,error,*999)
    NULLIFY(dependentField)
    CALL FieldVariable_FieldGet(fieldVariable,dependentField,err,error,*999)
    CALL Field_ScalingTypeGet(dependentField,scalingType,err,error,*999)
    NULLIFY(dependentDecomposition)
    CALL Field_DecompositionGet(dependentField,dependentDecomposition,err,error,*999)
    NULLIFY(geometricField)
    CALL Field_GeometricFieldGet(dependentField,geometricField,err,error,*999)
    NULLIFY(geometricDecomposition)
    CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
    NULLIFY(geometricVariable)
    CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
    NULLIFY(geometricInterpParameters)
    CALL FieldVariable_InterpolationParameterInitialise(geometricVariable,geometricInterpParameters,err,error,*999)
    NULLIFY(numericalInterpParameters)
    CALL FieldVariable_InterpolationParameterInitialise(fieldVariable,numericalInterpParameters,err,error,*999)
    NULLIFY(analyticInterpParameters)    
    CALL FieldVariable_InterpolationParameterInitialise(fieldVariable,analyticInterpParameters,err,error,*999)
    NULLIFY(geometricInterpPoint)
    CALL Field_InterpolatedPointInitialise(geometricInterpParameters,geometricInterpPoint,err,error,*999)
    NULLIFY(geometricInterpPointMetrics)
    CALL Field_InterpolatedPointMetricsInitialise(geometricInterpPoint,geometricInterpPointMetrics,err,error,*999)
    NULLIFY(dependentDomain)
    CALL Decomposition_DomainGet(dependentDecomposition,0,dependentDomain,err,error,*999)
    NULLIFY(dependentDomainTopology)
    CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
    NULLIFY(dependentDomainElements)
    CALL DomainTopology_DomainElementsGet(dependentDomainTopology,dependentDomainElements,err,error,*999)
    NULLIFY(geometricDomain)
    CALL Decomposition_DomainGet(geometricDecomposition,0,geometricDomain,err,error,*999)
    NULLIFY(geometricDomainTopology)
    CALL Domain_DomainTopologyGet(geometricDomain,geometricDomainTopology,err,error,*999)
    NULLIFY(geometricDomainElements)
    CALL DomainTopology_DomainElementsGet(geometricDomainTopology,geometricDomainElements,err,error,*999)
    CALL DomainElements_NumberOfElementsGet(dependentDomainElements,numberOfElements,err,error,*999)
    CALL DomainElements_TotalNumberOfElementsGet(dependentDomainElements,totalNumberOfElements,err,error,*999)

    !Cache component bases and quadrature schemes to avoid repeated calculations
    IF(numberOfComponents>MAX_NUMBER_OF_COMPONENTS) THEN
      localError="The number of components of "//TRIM(NumberToVString(numberOfComponents,"*",err,error))// &
        & " is greater than the maximum number of components of "//TRIM(NumberToVString(MAX_NUMBER_OF_COMPONENTS,"*", &
        & err,error))//". Increase MAX_NUMBER_OF_COMPONENTS."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    DO elementIdx=1,numberOfElements
      NULLIFY(geometricBasis)
      CALL DomainElements_ElementBasisGet(geometricDomainElements,elementIdx,geometricBasis,err,error,*999)
      CALL Basis_NumberOfXiGet(geometricBasis,numberOfXi,err,error,*999)
      DO componentIdx=1,numberOfComponents
        NULLIFY(dependentDomain)
        CALL FieldVariable_ComponentDomainGet(fieldVariable,componentIdx,dependentDomain,err,error,*999)
        NULLIFY(dependentDomainTopology)
        CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
        NULLIFY(dependentDomainElements)
        CALL DomainTopology_DomainElementsGet(dependentDomainTopology,dependentDomainElements,err,error,*999)
        NULLIFY(dependentBasis(componentIdx)%ptr)
        CALL DomainElements_ElementBasisGet(dependentDomainElements,elementIdx,dependentBasis(componentIdx)%ptr,err,error,*999)
        CALL Basis_NumberOfElementParametersGet(dependentBasis(componentIdx)%ptr,numberOfElementParameters(componentIdx), &
          & err,error,*999)
        NULLIFY(quadratureScheme(componentIdx)%ptr)
        CALL Basis_QuadratureSchemeGet(dependentBasis(componentIdx)%ptr,BASIS_DEFAULT_QUADRATURE_SCHEME, &
          & quadratureScheme(componentIdx)%ptr,err,error,*999)
      ENDDO !componentIdx
      CALL BasisQuadratureScheme_NumberOfGaussGet(quadratureScheme(1)%ptr,numberOfGauss,err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx,geometricInterpParameters,err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx,numericalInterpParameters,err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_ANALYTIC_VALUES_SET_TYPE,elementIdx,analyticInterpParameters, &
        & err,error,*999)
      DO gaussPointIdx=1,numberOfGauss
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
        CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
        CALL BasisQuadratureScheme_GaussWeightGet(quadratureScheme(1)%ptr,gaussPointIdx,gaussWeight,err,error,*999)
        jacobianGaussWeight=jacobian*gaussWeight
        DO componentIdx=1,numberOfComponents
          numericalIntegral=0.0_DP
          analyticIntegral=0.0_DP
          DO elementParameterIdx=1,numberOfElementParameters(componentIdx)
            CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureScheme(componentIdx)%ptr,elementParameterIdx,NO_PART_DERIV, &
              & gaussPointIdx,phi,err,error,*999)
            IF(scalingType==FIELD_NO_SCALING) THEN             
              numericalIntegral=numericalIntegral+phi*numericalInterpParameters%parameters(elementParameterIdx,componentIdx)
              analyticIntegral=analyticIntegral+phi*analyticInterpParameters%parameters(elementParameterIdx,componentIdx)
            ELSE
              numericalIntegral=numericalIntegral+phi*numericalInterpParameters%parameters(elementParameterIdx,componentIdx)* &
                & numericalInterpParameters%scaleFactors(elementParameterIdx,componentIdx)
              analyticIntegral=analyticIntegral+phi*analyticInterpParameters%parameters(elementParameterIdx,componentIdx)* &
                & analyticInterpParameters%scaleFactors(elementParameterIdx,componentIdx)
            ENDIF
          ENDDO !elementParameterIdx
          integralErrors(1,componentIdx)=integralErrors(1,componentIdx)+numericalIntegral*jacobianGaussWeight
          integralErrors(2,componentIdx)=integralErrors(2,componentIdx)+numericalIntegral**2*jacobianGaussWeight
          integralErrors(3,componentIdx)=integralErrors(3,componentIdx)+analyticIntegral*jacobianGaussWeight
          integralErrors(4,componentIdx)=integralErrors(4,componentIdx)+analyticIntegral**2*jacobianGaussWeight
          integralErrors(5,componentIdx)=integralErrors(5,componentIdx)+(analyticIntegral-numericalIntegral)*jacobianGaussWeight
          integralErrors(6,componentIdx)=integralErrors(6,componentIdx)+(analyticIntegral-numericalIntegral)**2*jacobianGaussWeight
        ENDDO !componentIdx
      ENDDO !gaussPointIdx
    ENDDO !elementIdx
    DO elementIdx=numberOfElements+1,totalNumberOfElements
      NULLIFY(geometricBasis)
      CALL DomainElements_ElementBasisGet(geometricDomainElements,elementIdx,geometricBasis,err,error,*999)
      CALL Basis_NumberOfXiGet(geometricBasis,numberOfXi,err,error,*999)
      DO componentIdx=1,numberOfComponents
        NULLIFY(dependentDomain)
        CALL FieldVariable_ComponentDomainGet(fieldVariable,componentIdx,dependentDomain,err,error,*999)
        NULLIFY(dependentDomainTopology)
        CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
        NULLIFY(dependentDomainElements)
        CALL DomainTopology_DomainElementsGet(dependentDomainTopology,dependentDomainElements,err,error,*999)
        NULLIFY(dependentBasis(componentIdx)%ptr)
        CALL DomainElements_ElementBasisGet(dependentDomainElements,elementIdx,dependentBasis(componentIdx)%ptr,err,error,*999)
        CALL Basis_NumberOfElementParametersGet(dependentBasis(componentIdx)%ptr,numberOfElementParameters(componentIdx), &
          & err,error,*999)
        NULLIFY(quadratureScheme(componentIdx)%ptr)
        CALL Basis_QuadratureSchemeGet(dependentBasis(componentIdx)%ptr,BASIS_DEFAULT_QUADRATURE_SCHEME, &
          & quadratureScheme(componentIdx)%ptr,err,error,*999)
      ENDDO !componentIdx
      CALL BasisQuadratureScheme_NumberOfGaussGet(quadratureScheme(1)%ptr,numberOfGauss,err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx,geometricInterpParameters,err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx,numericalInterpParameters,err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_ANALYTIC_VALUES_SET_TYPE,elementIdx,analyticInterpParameters, &
        & err,error,*999)     
      DO gaussPointIdx=1,numberOfGauss
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
        CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
        CALL BasisQuadratureScheme_GaussWeightGet(quadratureScheme(1)%ptr,gaussPointIdx,gaussWeight,err,error,*999)
        jacobianGaussWeight=jacobian*gaussWeight
        DO componentIdx=1,numberOfComponents
          numericalIntegral=0.0_DP
          analyticIntegral=0.0_DP
          DO elementParameterIdx=1,numberOfElementParameters(componentIdx)
            CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureScheme(componentIdx)%ptr,elementParameterIdx,NO_PART_DERIV, &
              & gaussPointIdx,phi,err,error,*999)
            IF(scalingType==FIELD_NO_SCALING) THEN             
              numericalIntegral=numericalIntegral+phi*numericalInterpParameters%parameters(elementParameterIdx,componentIdx)
              analyticIntegral=analyticIntegral+phi*analyticInterpParameters%parameters(elementParameterIdx,componentIdx)
            ELSE
              numericalIntegral=numericalIntegral+phi*numericalInterpParameters%parameters(elementParameterIdx,componentIdx)* &
                & numericalInterpParameters%scaleFactors(elementParameterIdx,componentIdx)
              analyticIntegral=analyticIntegral+phi*analyticInterpParameters%parameters(elementParameterIdx,componentIdx)* &
                & analyticInterpParameters%scaleFactors(elementParameterIdx,componentIdx)
            ENDIF
          ENDDO !elementParameterIdx
          ghostIntegralErrors(1,componentIdx)=ghostIntegralErrors(1,componentIdx)+numericalIntegral*jacobianGaussWeight
          ghostIntegralErrors(2,componentIdx)=ghostIntegralErrors(2,componentIdx)+numericalIntegral**2*jacobianGaussWeight
          ghostIntegralErrors(3,componentIdx)=ghostIntegralErrors(3,componentIdx)+analyticIntegral*jacobianGaussWeight
          ghostIntegralErrors(4,componentIdx)=ghostIntegralErrors(4,componentIdx)+analyticIntegral**2*jacobianGaussWeight
          ghostIntegralErrors(5,componentIdx)=ghostIntegralErrors(5,componentIdx)+ &
            & (analyticIntegral-numericalIntegral)*jacobianGaussWeight
          ghostIntegralErrors(6,componentIdx)=ghostIntegralErrors(6,componentIdx)+ &
            & (analyticIntegral-numericalIntegral)**2*jacobianGaussWeight
        ENDDO !componentIdx
      ENDDO !gaussPointIdx
    ENDDO !elementIdx
    CALL Field_InterpolatedPointMetricsFinalise(geometricInterpPointMetrics,err,error,*999)
    CALL Field_InterpolatedPointFinalise(geometricInterpPoint,err,error,*999)
    CALL FieldVariable_InterpolationParameterFinalise(analyticInterpParameters,err,error,*999)
    CALL FieldVariable_InterpolationParameterFinalise(numericalInterpParameters,err,error,*999)
    CALL FieldVariable_InterpolationParameterFinalise(geometricInterpParameters,err,error,*999)

    EXITS("AnalyticAnalysis_IntegralErrors")
    RETURN
999 ERRORSEXITS("AnalyticAnalysis_IntegralErrors",err,error)
    RETURN 1
    
  END SUBROUTINE AnalyticAnalysis_IntegralErrors

  !
  !================================================================================================================================
  !

  !>Get integral absolute error value for the field
  SUBROUTINE AnalyticAnalysis_IntegralAbsoluteErrorGet(field,variableType,componentNumber,integralError, &
    & ghostIntegralError,err,error,*)

    !Argument variables   
    TYPE(FieldType), POINTER :: field !<the field.
    INTEGER(INTG), INTENT(IN) :: variableType !<variable type
    INTEGER(INTG), INTENT(IN) :: componentNumber !<component number
    REAL(DP), INTENT(OUT) :: integralError(:) !<On return, the integral numerical value for local elements
    REAL(DP), INTENT(OUT) :: ghostIntegralError(:) !<On return, the integral numerical for global elements
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfComponents
    REAL(DP), ALLOCATABLE :: ghostIntegralErrors(:,:),integralErrors(:,:)
    TYPE(FieldVariableType), POINTER :: fieldVariable
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("AnalyticAnalysis_IntegralAbsoluteErrorGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(SIZE(integralError,1)<2) THEN
      localError="The size of the integral error array of "//TRIM(NumberToVString(SIZE(integralError,1),"*",err,error))// &
        & " is too small. The size of the array must be >= 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(ghostIntegralError,1)<2) THEN
      localError="The size of the ghost integral error array of "// &
        & TRIM(NumberToVString(SIZE(ghostIntegralError,1),"*",err,error))//" is too small. The size of the array must be >= 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL FieldVariable_NumberOfComponentsGet(fieldVariable,numberOfComponents,err,error,*999)
    ALLOCATE(integralErrors(6,numberOfComponents),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate integral errors.",err,error,*999)
    ALLOCATE(ghostIntegralErrors(6,numberOfComponents),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate ghost integral errors.",err,error,*999)
    CALL AnalyticAnalysis_IntegralErrors(fieldVariable,integralErrors,ghostIntegralErrors,err,error,*999)
    integralError(1)=AnalyticAnalysis_AbsoluteError(integralErrors(1,componentNumber),integralErrors(3,componentNumber))
    integralError(2)=AnalyticAnalysis_AbsoluteError(integralErrors(2,componentNumber),integralErrors(4,componentNumber))
    ghostIntegralError(1)=AnalyticAnalysis_AbsoluteError(ghostIntegralErrors(1,componentNumber), &
      & ghostIntegralErrors(3,componentNumber))
    ghostIntegralError(2)=AnalyticAnalysis_AbsoluteError(ghostIntegralErrors(2,componentNumber), &
      & ghostIntegralErrors(4,componentNumber))
    
    DEALLOCATE(integralErrors)
    DEALLOCATE(ghostIntegralErrors)
    
    EXITS("AnalyticAnalysis_IntegralAbsoluteErrorGet")
    RETURN
999 IF(ALLOCATED(integralErrors)) DEALLOCATE(integralErrors)
    IF(ALLOCATED(ghostIntegralErrors)) DEALLOCATE(ghostIntegralErrors)    
    ERRORSEXITS("AnalyticAnalysis_IntegralAbsoluteErrorGet",err,error)
    RETURN 1
    
  END SUBROUTINE AnalyticAnalysis_IntegralAbsoluteErrorGet
  
   !
  !================================================================================================================================
  !

  !>Get integral analytic value for the field TODO should we use analytical formula to calculate the integration?
  SUBROUTINE AnalyticAnalysis_IntegralAnalyticValueGet(field,variableType,componentNumber,integralError, &
    & ghostIntegralError,err,error,*)
  
    !Argument variables
    TYPE(FieldType), POINTER :: field !<the field.
    INTEGER(INTG), INTENT(IN) :: variableType !<variable type
    INTEGER(INTG), INTENT(IN) :: componentNumber !<component number
    REAL(DP), INTENT(OUT) :: integralError(:) !<On return, the integral numerical value for local elements
    REAL(DP), INTENT(OUT) :: ghostIntegralError(:) !<On return, the integral numerical for global elements
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfComponents
    REAL(DP), ALLOCATABLE ::  ghostIntegralErrors(:,:),integralErrors(:,:)
    TYPE(FieldVariableType), POINTER :: fieldVariable
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("AnalyticAnalysis_IntegralAnalyticValueGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(SIZE(integralError,1)<2) THEN
      localError="The size of the integral error array of "//TRIM(NumberToVString(SIZE(integralError,1),"*",err,error))// &
        & " is too small. The size of the array must be >= 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(ghostIntegralError,1)<2) THEN
      localError="The size of the ghost integral error array of "// &
        & TRIM(NumberToVString(SIZE(ghostIntegralError,1),"*",err,error))//" is too small. The size of the array must be >= 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL FieldVariable_NumberOfComponentsGet(fieldVariable,numberOfComponents,err,error,*999)
    ALLOCATE(integralErrors(6,numberOfComponents),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate integral errors.",err,error,*999)
    ALLOCATE(ghostIntegralErrors(6,numberOfComponents),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate ghost integral errors.",err,error,*999)
    CALL AnalyticAnalysis_IntegralErrors(fieldVariable,integralErrors,ghostIntegralErrors,err,error,*999)
    integralError(1)=integralErrors(3,componentNumber)
    integralError(2)=integralErrors(4,componentNumber)
    ghostIntegralError(1)=ghostIntegralErrors(3,componentNumber)
    ghostIntegralError(2)=ghostIntegralErrors(4,componentNumber)
    
    DEALLOCATE(integralErrors)
    DEALLOCATE(ghostIntegralErrors)
    
    EXITS("AnalyticAnalysis_IntegralAnalyticValueGet")
    RETURN
999 IF(ALLOCATED(integralErrors)) DEALLOCATE(integralErrors)
    IF(ALLOCATED(ghostIntegralErrors)) DEALLOCATE(ghostIntegralErrors)    
    ERRORSEXITS("AnalyticAnalysis_IntegralAnalyticValueGet",err,error)
    RETURN 1
    
  END SUBROUTINE AnalyticAnalysis_IntegralAnalyticValueGet
  
  !
  !================================================================================================================================
  !

  !>Get integral numerical value for the field, TODO check integral calculation
  SUBROUTINE AnalyticAnalysis_IntegralNumericalValueGet(field,variableType,componentNumber,integralError, &
    & ghostIntegralError,err,error,*)
  
    !Argument variables   
    TYPE(FieldType), POINTER :: field !<the field.
    INTEGER(INTG), INTENT(IN) :: variableType !<variable type
    INTEGER(INTG), INTENT(IN) :: componentNumber !<component number
    REAL(DP), INTENT(OUT) :: integralError(:) !<On return, the integral numerical value for local elements
    REAL(DP), INTENT(OUT) :: ghostIntegralError(:) !<On return, the integral numerical for global elements
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfComponents
    REAL(DP), ALLOCATABLE ::  ghostIntegralErrors(:,:),integralErrors(:,:)
    TYPE(FieldVariableType), POINTER :: fieldVariable
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("AnalyticAnalysis_IntegralNumericalValueGet",err,error,*999)       

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(SIZE(integralError,1)<2) THEN
      localError="The size of the integral error array of "//TRIM(NumberToVString(SIZE(integralError,1),"*",err,error))// &
        & " is too small. The size of the array must be >= 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(ghostIntegralError,1)<2) THEN
      localError="The size of the ghost integral error array of "// &
        & TRIM(NumberToVString(SIZE(ghostIntegralError,1),"*",err,error))//" is too small. The size of the array must be >= 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL FieldVariable_NumberOfComponentsGet(fieldVariable,numberOfComponents,err,error,*999)
    ALLOCATE(integralErrors(6,numberOfComponents),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate integral errors.",err,error,*999)
    ALLOCATE(ghostIntegralErrors(6,numberOfComponents),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate ghost integral errors.",err,error,*999)
    CALL AnalyticAnalysis_IntegralErrors(fieldVariable,integralErrors,ghostIntegralErrors,err,error,*999)
    CALL AnalyticAnalysis_IntegralErrors(fieldVariable,integralErrors,ghostIntegralErrors,err,error,*999)
    integralError(1)=integralErrors(1,componentNumber)
    integralError(2)=integralErrors(2,componentNumber)
    ghostIntegralError(1)=ghostIntegralErrors(1,componentNumber)
    ghostIntegralError(2)=ghostIntegralErrors(2,componentNumber)

    DEALLOCATE(integralErrors)
    DEALLOCATE(ghostIntegralErrors)
    
    EXITS("AnalyticAnalysis_IntegralNumericalValueGet")
    RETURN
999 IF(ALLOCATED(integralErrors)) DEALLOCATE(integralErrors)
    IF(ALLOCATED(ghostIntegralErrors)) DEALLOCATE(ghostIntegralErrors)    
    ERRORSEXITS("AnalyticAnalysis_IntegralNumericalValueGet",err,error)
    RETURN 1
    
  END SUBROUTINE AnalyticAnalysis_IntegralNumericalValueGet
  
  !
  !================================================================================================================================
  !

  !>Get integral nid numerical value for the field, TODO check integral calculation
  SUBROUTINE AnalyticAnalysis_IntegralNIDNumericalValueGet(field,variableType,componentNumber,integralError, &
    & ghostIntegralError,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<the field.
    INTEGER(INTG), INTENT(IN) :: variableType !<variable type
    INTEGER(INTG), INTENT(IN) :: componentNumber !<component number
    REAL(DP), INTENT(OUT) :: integralError(:) !<On return, the integral numerical value for local elements
    REAL(DP), INTENT(OUT) :: ghostIntegralError(:) !<On return, the integral numerical for global elements
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfComponents
    REAL(DP), ALLOCATABLE ::  ghostIntegralErrors(:,:),integralErrors(:,:)
    TYPE(FieldVariableType), POINTER :: fieldVariable
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("AnalyticAnalysis_IntegralNIDNumericalValueGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(SIZE(integralError,1)<2) THEN
      localError="The size of the integral error array of "//TRIM(NumberToVString(SIZE(integralError,1),"*",err,error))// &
        & " is too small. The size of the array must be >= 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(ghostIntegralError,1)<2) THEN
      localError="The size of the ghost integral error array of "// &
        & TRIM(NumberToVString(SIZE(ghostIntegralError,1),"*",err,error))//" is too small. The size of the array must be >= 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL FieldVariable_NumberOfComponentsGet(fieldVariable,numberOfComponents,err,error,*999)
    ALLOCATE(integralErrors(6,numberOfComponents),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate integral errors.",err,error,*999)
    ALLOCATE(ghostIntegralErrors(6,numberOfComponents),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate ghost integral errors.",err,error,*999)
    CALL AnalyticAnalysis_IntegralErrors(fieldVariable,integralErrors,ghostIntegralErrors,err,error,*999)
    integralError(1)=integralErrors(5,componentNumber)
    integralError(2)=integralErrors(6,componentNumber)
    ghostIntegralError(1)=ghostIntegralErrors(5,componentNumber)
    ghostIntegralError(2)=ghostIntegralErrors(6,componentNumber)

    DEALLOCATE(integralErrors)
    DEALLOCATE(ghostIntegralErrors)
    
    EXITS("AnalyticAnalysis_IntegralNIDNumericalValueGet")
    RETURN
999 IF(ALLOCATED(integralErrors)) DEALLOCATE(integralErrors)
    IF(ALLOCATED(ghostIntegralErrors)) DEALLOCATE(ghostIntegralErrors)    
    ERRORS("AnalyticAnalysis_IntegralNIDNumericalValueGet",err,error)
    EXITS("AnalyticAnalysis_IntegralNIDNumericalValueGet")
    RETURN 1
    
  END SUBROUTINE AnalyticAnalysis_IntegralNIDNumericalValueGet

  !
  !================================================================================================================================
  !

  !>Get integral nid error value for the field
  SUBROUTINE AnalyticAnalysis_IntegralNIDErrorGet(field,variableType,componentNumber,integralError, &
    & ghostIntegralError,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<the field.
    INTEGER(INTG), INTENT(IN) :: variableType !<variable type
    INTEGER(INTG), INTENT(IN) :: componentNumber !<component number
    REAL(DP), INTENT(OUT) :: integralError(:) !<On return, the integral numerical value for local elements
    REAL(DP), INTENT(OUT) :: ghostIntegralError(:) !<On return, the integral numerical for global elements
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfComponents
    REAL(DP), ALLOCATABLE ::  ghostIntegralErrors(:,:),integralErrors(:,:)
    TYPE(FieldVariableType), POINTER :: fieldVariable
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("AnalyticAnalysis_IntegralNIDErrorGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(SIZE(integralError,1)<2) THEN
      localError="The size of the integral error array of "//TRIM(NumberToVString(SIZE(integralError,1),"*",err,error))// &
        & " is too small. The size of the array must be >= 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(ghostIntegralError,1)<2) THEN
      localError="The size of the ghost integral error array of "// &
        & TRIM(NumberToVString(SIZE(ghostIntegralError,1),"*",err,error))//" is too small. The size of the array must be >= 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL FieldVariable_NumberOfComponentsGet(fieldVariable,numberOfComponents,err,error,*999)
    ALLOCATE(integralErrors(6,numberOfComponents),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate integral errors.",err,error,*999)
    ALLOCATE(ghostIntegralErrors(6,numberOfComponents),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate ghost integral errors.",err,error,*999)
    CALL AnalyticAnalysis_IntegralErrors(fieldVariable,integralErrors,ghostIntegralErrors,err,error,*999)
    integralError(1)=AnalyticAnalysis_NIDError(integralErrors(5,componentNumber),integralErrors(3,componentNumber))
    integralError(2)=AnalyticAnalysis_NIDError(integralErrors(6,componentNumber),integralErrors(4,componentNumber))
    ghostIntegralError(1)=AnalyticAnalysis_NIDError(ghostIntegralErrors(5,componentNumber), &
      & ghostIntegralErrors(3,componentNumber))
    ghostIntegralError(2)=AnalyticAnalysis_NIDError(ghostIntegralErrors(6,componentNumber), &
      & ghostIntegralErrors(4,componentNumber))
 
    DEALLOCATE(integralErrors)
    DEALLOCATE(ghostIntegralErrors)
    
    EXITS("AnalyticAnalysis_IntegralNIDErrorGet")
    RETURN
999 IF(ALLOCATED(integralErrors)) DEALLOCATE(integralErrors)
    IF(ALLOCATED(ghostIntegralErrors)) DEALLOCATE(ghostIntegralErrors)    
    ERRORSEXITS("AnalyticAnalysis_IntegralNIDErrorGet",err,error)
    RETURN 1
    
  END SUBROUTINE AnalyticAnalysis_IntegralNIDErrorGet

  !
  !================================================================================================================================
  !

  !>Get integral percentage error value for the field
  SUBROUTINE AnalyticAnalysis_IntegralPercentageErrorGet(field,variableType,componentNumber,integralError, &
    & ghostIntegralError,err,error,*)

    !Argument variables   
    TYPE(FieldType), POINTER :: field !<the field.
    INTEGER(INTG), INTENT(IN) :: variableType !<variable type
    INTEGER(INTG), INTENT(IN) :: componentNumber !<component number
    REAL(DP), INTENT(OUT) :: integralError(:) !<On return, the integral numerical value for local elements
    REAL(DP), INTENT(OUT) :: ghostIntegralError(:) !<On return, the integral numerical for global elements
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfComponents
    REAL(DP), ALLOCATABLE ::  ghostIntegralErrors(:,:),integralErrors(:,:)
    TYPE(FieldVariableType), POINTER :: fieldVariable
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("AnalyticAnalysis_IntegralPercentageErrorGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(SIZE(integralError,1)<2) THEN
      localError="The size of the integral error array of "//TRIM(NumberToVString(SIZE(integralError,1),"*",err,error))// &
        & " is too small. The size of the array must be >= 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(ghostIntegralError,1)<2) THEN
      localError="The size of the ghost integral error array of "// &
        & TRIM(NumberToVString(SIZE(ghostIntegralError,1),"*",err,error))//" is too small. The size of the array must be >= 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL FieldVariable_NumberOfComponentsGet(fieldVariable,numberOfComponents,err,error,*999)
    ALLOCATE(integralErrors(6,numberOfComponents),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate integral errors.",err,error,*999)
    ALLOCATE(ghostIntegralErrors(6,numberOfComponents),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate ghost integral errors.",err,error,*999)
    CALL AnalyticAnalysis_IntegralErrors(fieldVariable,integralErrors,ghostIntegralErrors,err,error,*999)
    integralError(1)=AnalyticAnalysis_PercentageError(integralErrors(1,componentNumber),integralErrors(3,componentNumber))
    integralError(2)=AnalyticAnalysis_PercentageError(integralErrors(2,componentNumber),integralErrors(4,componentNumber))
    ghostIntegralError(1)=AnalyticAnalysis_PercentageError(ghostIntegralErrors(1,componentNumber), &
      & ghostIntegralErrors(3,componentNumber))
    ghostIntegralError(2)=AnalyticAnalysis_PercentageError(ghostIntegralErrors(2,componentNumber), &
      & ghostIntegralErrors(4,componentNumber))
    
    DEALLOCATE(integralErrors)
    DEALLOCATE(ghostIntegralErrors)
    
    EXITS("AnalyticAnalysis_IntegralPercentageErrorGet")
    RETURN
999 IF(ALLOCATED(integralErrors)) DEALLOCATE(integralErrors)
    IF(ALLOCATED(ghostIntegralErrors)) DEALLOCATE(ghostIntegralErrors)    
    ERRORSEXITS("AnalyticAnalysis_IntegralPercentageErrorGet",err,error)
    RETURN 1
    
  END SUBROUTINE AnalyticAnalysis_IntegralPercentageErrorGet
  
   !
  !================================================================================================================================
  !

  !>Get integral relative error value for the field.
  SUBROUTINE AnalyticAnalysis_IntegralRelativeErrorGet(field,variableType,componentNumber,integralError, &
    & ghostIntegralError,err,error,*)
  
    !Argument variables   
    TYPE(FieldType), POINTER :: field !<the field.
    INTEGER(INTG), INTENT(IN) :: variableType !<variable type
    INTEGER(INTG), INTENT(IN) :: componentNumber !<component number
    REAL(DP), INTENT(OUT) :: integralError(:) !<On return, the integral numerical value for local elements
    REAL(DP), INTENT(OUT) :: ghostIntegralError(:) !<On return, the integral numerical for global elements
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfComponents
    REAL(DP), ALLOCATABLE ::  ghostIntegralErrors(:,:),integralErrors(:,:)
    TYPE(FieldVariableType), POINTER :: fieldVariable
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("AnalyticAnalysis_IntegralRelativeErrorGet",err,error,*999)       

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(SIZE(integralError,1)<2) THEN
      localError="The size of the integral error array of "//TRIM(NumberToVString(SIZE(integralError,1),"*",err,error))// &
        & " is too small. The size of the array must be >= 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(ghostIntegralError,1)<2) THEN
      localError="The size of the ghost integral error array of "// &
        & TRIM(NumberToVString(SIZE(ghostIntegralError,1),"*",err,error))//" is too small. The size of the array must be >= 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL FieldVariable_NumberOfComponentsGet(fieldVariable,numberOfComponents,err,error,*999)
    ALLOCATE(integralErrors(6,numberOfComponents),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate integral errors.",err,error,*999)
    ALLOCATE(ghostIntegralErrors(6,numberOfComponents),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate ghost integral errors.",err,error,*999)
    CALL AnalyticAnalysis_IntegralErrors(fieldVariable,integralErrors,ghostIntegralErrors,err,error,*999)
    integralError(1)=AnalyticAnalysis_RelativeError(integralErrors(1,componentNumber),integralErrors(3,componentNumber))
    integralError(2)=AnalyticAnalysis_RelativeError(integralErrors(2,componentNumber),integralErrors(4,componentNumber))
    ghostIntegralError(1)=AnalyticAnalysis_RelativeError(ghostIntegralErrors(1,componentNumber), &
      & ghostIntegralErrors(3,componentNumber))
    ghostIntegralError(2)=AnalyticAnalysis_RelativeError(ghostIntegralErrors(2,componentNumber), &
      & ghostIntegralErrors(4,componentNumber))
   
    DEALLOCATE(integralErrors)
    DEALLOCATE(ghostIntegralErrors)
    
    EXITS("AnalyticAnalysis_IntegralRelativeErrorGet")
    RETURN
999 IF(ALLOCATED(integralErrors)) DEALLOCATE(integralErrors)
    IF(ALLOCATED(ghostIntegralErrors)) DEALLOCATE(ghostIntegralErrors)    
    ERRORSEXITS("AnalyticAnalysis_IntegralRelativeErrorGet",err,error)
    RETURN 1
    
  END SUBROUTINE AnalyticAnalysis_IntegralRelativeErrorGet

  !
  !================================================================================================================================
  !

  !>Get absolute error value for the node
  SUBROUTINE AnalyticAnalysis_AbsoluteErrorGetNode(field,variableType,versionNumber,derivativeNumber,userNodeNumber, &
    & componentNumber,value,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<the field.
    INTEGER(INTG), INTENT(IN) :: versionNumber !<derivative version number
    INTEGER(INTG), INTENT(IN) :: derivativeNumber !<derivative number
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<node number
    INTEGER(INTG), INTENT(IN) :: componentNumber !<component number
    INTEGER(INTG), INTENT(IN) :: variableType !<variable type
    REAL(DP), INTENT(OUT) :: VALUE !<On return, the absolute error
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: numericalValue,analyticValue

    ENTERS("AnalyticAnalysis_AbsoluteErrorGetNode",err,error,*999)

    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)

    CALL Field_ParameterSetGetNode(field,variableType,FIELD_VALUES_SET_TYPE,versionNumber,derivativeNumber, &
        & userNodeNumber,componentNumber,numericalValue,err,error,*999)
    CALL Field_ParameterSetGetNode(field,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,versionNumber,derivativeNumber, &
      & userNodeNumber,componentNumber,analyticValue,err,error,*999)
    value=AnalyticAnalysis_AbsoluteError(numericalValue,analyticValue)
 
    EXITS("AnalyticAnalysis_AbsoluteErrorGetNode")
    RETURN
999 ERRORSEXITS("AnalyticAnalysis_AbsoluteErrorGetNode",err,error)
    RETURN 1
    
  END SUBROUTINE AnalyticAnalysis_AbsoluteErrorGetNode

  !
  !================================================================================================================================
  !

  !>Get percentage error value for the node
  SUBROUTINE AnalyticAnalysis_PercentageErrorGetNode(field,variableType,versionNumber,derivativeNumber,userNodeNumber, &
    & componentNumber,value,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<the field.
    INTEGER(INTG), INTENT(IN) :: versionNumber !<derivative version number
    INTEGER(INTG), INTENT(IN) :: derivativeNumber !<derivative number
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<node number
    INTEGER(INTG), INTENT(IN) :: componentNumber !<component number
    INTEGER(INTG), INTENT(IN) :: variableType !<variable type
    REAL(DP), INTENT(OUT) :: VALUE !<On return, the percentage error
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: numericalValue, analyticValue

    ENTERS("AnalyticAnalysis_PercentageErrorGetNode",err,error,*999)

    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    
    CALL Field_ParameterSetGetNode(field,variableType,FIELD_VALUES_SET_TYPE,versionNumber,derivativeNumber, &
      & userNodeNumber,componentNumber,numericalValue,err,error,*999)
    CALL Field_ParameterSetGetNode(field,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,versionNumber,derivativeNumber, &
      & userNodeNumber,componentNumber,analyticValue,err,error,*999)
    value=AnalyticAnalysis_PercentageError(numericalValue,analyticValue)

    EXITS("AnalyticAnalysis_PercentageErrorGetNode")
    RETURN
999 ERRORSEXITS("AnalyticAnalysis_PercentageErrorGetNode",err,error)
    RETURN 1
    
  END SUBROUTINE AnalyticAnalysis_PercentageErrorGetNode

  !
  !================================================================================================================================
  !

  !>Get relative error value for the node
  SUBROUTINE AnalyticAnalysis_RelativeErrorGetNode(field,variableType,versionNumber,derivativeNumber,userNodeNumber, &
    & componentNumber,value,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<the field.
    INTEGER(INTG), INTENT(IN) :: versionNumber !<derivative version number
    INTEGER(INTG), INTENT(IN) :: derivativeNumber !<derivative number
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<node number
    INTEGER(INTG), INTENT(IN) :: componentNumber !<component number
    INTEGER(INTG), INTENT(IN) :: variableType !<variable type
    REAL(DP), INTENT(OUT) :: VALUE !<On return, the relative error
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: numericalValue,analyticValue

    ENTERS("AnalyticAnalysis_RelativeErrorGetNode",err,error,*999)

    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    
    CALL Field_ParameterSetGetNode(field,variableType,FIELD_VALUES_SET_TYPE,versionNumber,derivativeNumber, &
      & userNodeNumber,componentNumber,numericalValue,err,error,*999)
    CALL Field_ParameterSetGetNode(field,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,versionNumber,derivativeNumber, &
      & userNodeNumber,componentNumber,analyticValue,err,error,*999)
    value=AnalyticAnalysis_RelativeError(numericalValue,analyticValue)

    EXITS("AnalyticAnalysis_RelativeErrorGetNode")
    RETURN
999 ERRORSEXITS("AnalyticAnalysis_RelativeErrorGetNode",err,error)
    RETURN 1
    
  END SUBROUTINE AnalyticAnalysis_RelativeErrorGetNode

  !
  !================================================================================================================================
  !

  !>Get absolute error value for an element
  SUBROUTINE AnalyticAnalysis_AbsoluteErrorGetElement(field,variableType,userElementNumber,componentNumber,value,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<the field.
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<node number
    INTEGER(INTG), INTENT(IN) :: componentNumber !<component number
    INTEGER(INTG), INTENT(IN) :: variableType !<variable type
    REAL(DP), INTENT(OUT) :: value !<On return, the absolute error
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: numericalValue,analyticValue

    ENTERS("AnalyticAnalysis_AbsoluteErrorGetElement",err,error,*999)

    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    
    CALL Field_ParameterSetGetElement(field,variableType,FIELD_VALUES_SET_TYPE,userElementNumber,componentNumber, &
      & numericalValue,err,error,*999)
    CALL Field_ParameterSetGetElement(field,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,userElementNumber, &
      & componentNumber,analyticValue,err,error,*999)
    value=AnalyticAnalysis_AbsoluteError(numericalValue,analyticValue)

    EXITS("AnalyticAnalysis_AbsoluteErrorGetElement")
    RETURN
999 ERRORSEXITS("AnalyticAnalysis_AbsoluteErrorGetElement",err,error)
    RETURN 1
    
  END SUBROUTINE AnalyticAnalysis_AbsoluteErrorGetElement

  !
  !================================================================================================================================
  !

  !>Get percentage error value for the node
  SUBROUTINE AnalyticAnalysis_PercentageErrorGetElement(field,variableType,userElementNumber,componentNumber,value,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<the field.
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<node number
    INTEGER(INTG), INTENT(IN) :: componentNumber !<component number
    INTEGER(INTG), INTENT(IN) :: variableType !<variable type
    REAL(DP), INTENT(OUT) :: value !<On return, the percentage error
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: numericalValue,analyticValue

    ENTERS("AnalyticAnalysis_PercentageErrorGetElement",err,error,*999)

    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    
    CALL Field_ParameterSetGetElement(field,variableType,FIELD_VALUES_SET_TYPE,userElementNumber,componentNumber, &
      & numericalValue,err,error,*999)
    CALL Field_ParameterSetGetElement(field,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,userElementNumber, &
      & componentNumber,analyticValue,err,error,*999)
    value=AnalyticAnalysis_PercentageError(numericalValue,analyticValue)

    EXITS("AnalyticAnalysis_PercentageErrorGetElement")
    RETURN
999 ERRORSEXITS("AnalyticAnalysis_PercentageErrorGetElement",err,error)
    RETURN 1
    
  END SUBROUTINE AnalyticAnalysis_PercentageErrorGetElement


  !
  !================================================================================================================================
  !

  !>Get relative error value for an element
  SUBROUTINE AnalyticAnalysis_RelativeErrorGetElement(field,variableType,userElementNumber,componentNumber,value,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<the field.
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<node number
    INTEGER(INTG), INTENT(IN) :: componentNumber !<component number
    INTEGER(INTG), INTENT(IN) :: variableType !<variable type
    REAL(DP), INTENT(OUT) :: value !<On return, the relative error
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: numericalValue,analyticValue

    ENTERS("AnalyticAnalysis_RelativeErrorGetElement",err,error,*999)

    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    
    CALL Field_ParameterSetGetElement(field,variableType,FIELD_VALUES_SET_TYPE,userElementNumber,componentNumber, &
      & numericalValue,err,error,*999)
    CALL Field_ParameterSetGetElement(field,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,userElementNumber, &
      & componentNumber,analyticValue,err,error,*999)
    value=AnalyticAnalysis_RelativeError(numericalValue,analyticValue)
 
    EXITS("AnalyticAnalysis_RelativeErrorGetElement")
    RETURN
999 ERRORSEXITS("AnalyticAnalysis_RelativeErrorGetElement",err,error)
    RETURN 1
    
  END SUBROUTINE AnalyticAnalysis_RelativeErrorGetElement

  !
  !================================================================================================================================
  !

  !>Get absolute error value for the node
  SUBROUTINE AnalyticAnalysis_AbsoluteErrorGetConstant(field,variableType,componentNumber,value,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<the field.
    INTEGER(INTG), INTENT(IN) :: componentNumber !<component number
    INTEGER(INTG), INTENT(IN) :: variableType !<variable type
    REAL(DP), INTENT(OUT) :: value !<On return, the absolute error
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: numericalValue,analyticValue

    ENTERS("AnalyticAnalysis_AbsoluteErrorGetConstant",err,error,*999)

    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    
    CALL Field_ParameterSetGetConstant(field,variableType,FIELD_VALUES_SET_TYPE,componentNumber,numericalValue,err,error, &
      & *999)
    CALL Field_ParameterSetGetConstant(field,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,componentNumber,analyticValue, &
      & err,error,*999)
    value=AnalyticAnalysis_AbsoluteError(numericalValue,analyticValue)

    EXITS("AnalyticAnalysis_AbsoluteErrorGetConstant")
    RETURN
999 ERRORSEXITS("AnalyticAnalysis_AbsoluteErrorGetConstant",err,error)
    RETURN 1
    
  END SUBROUTINE AnalyticAnalysis_AbsoluteErrorGetConstant

  !
  !================================================================================================================================
  !

  !>Get percentage error value for a constant
  SUBROUTINE AnalyticAnalysis_PercentageErrorGetConstant(field,variableType,componentNumber,value,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<the field.
    INTEGER(INTG), INTENT(IN) :: componentNumber !<component number
    INTEGER(INTG), INTENT(IN) :: variableType !<variable type
    REAL(DP), INTENT(OUT) :: value !<On return, the percentage error
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: numericalValue,analyticValue

    ENTERS("AnalyticAnalysis_PercentageErrorGetConstant",err,error,*999)

    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    
    CALL Field_ParameterSetGetConstant(field,variableType,FIELD_VALUES_SET_TYPE,componentNumber,numericalValue,err,error,*999)
    CALL Field_ParameterSetGetConstant(field,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,componentNumber,analyticValue, &
      & err,error,*999)
    value=AnalyticAnalysis_PercentageError(numericalValue,analyticValue)
 
    EXITS("AnalyticAnalysis_PercentageErrorGetConstant")
    RETURN
999 ERRORSEXITS("AnalyticAnalysis_PercentageErrorGetConstant",err,error)
    RETURN 1
    
  END SUBROUTINE AnalyticAnalysis_PercentageErrorGetConstant


  !
  !================================================================================================================================
  !

  !>Get relative error value for a constant
  SUBROUTINE AnalyticAnalysis_RelativeErrorGetConstant(field,variableType,componentNumber,value,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<the field.
    INTEGER(INTG), INTENT(IN) :: componentNumber !<component number
    INTEGER(INTG), INTENT(IN) :: variableType !<variable type
    REAL(DP), INTENT(OUT) :: value !<On return, the relative error
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: numericalValue,analyticValue

    ENTERS("AnalyticAnalysis_RelativeErrorGetConstant",err,error,*999)

    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    
    CALL Field_ParameterSetGetConstant(field,variableType,FIELD_VALUES_SET_TYPE,componentNumber,numericalValue,err,error, &
      & *999)
    CALL Field_ParameterSetGetConstant(field,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,componentNumber,analyticValue, &
      & err,error,*999)
    value=AnalyticAnalysis_RelativeError(numericalValue,analyticValue)

    EXITS("AnalyticAnalysis_RelativeErrorGetConstant")
    RETURN
999 ERRORSEXITS("AnalyticAnalysis_RelativeErrorGetConstant",err,error)
    RETURN 1
    
  END SUBROUTINE AnalyticAnalysis_RelativeErrorGetConstant

  !
  !================================================================================================================================
  !

  !>Get rms error value for the field
  SUBROUTINE AnalyticAnalysis_RMSErrorGetNode(field,variableType,componentNumber,errorType,localRMS,localGhostRMS, &
    & globalRMS,err,error,*)
  
    !Argument variables   
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the RMS error for.
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the RMS error for \see \see FieldRoutines_VariableTypes,FieldRoutines
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The field variable component number to get the RMS error for. 
    INTEGER(INTG), INTENT(IN) :: errorType !<The error type to get the RMS error for \see AnalyticAnalysisRoutines_ErrorTypes,AnalyticAnalysisRoutines
    REAL(DP), INTENT(OUT) :: localRMS(:) !<localRMS(derivativeIdx). On return, the local RMS error
    REAL(DP), INTENT(OUT) :: localGhostRMS(:) !<localGhostRMS(derivativeIdx). On return, the local + ghost RMS error
    REAL(DP), INTENT(OUT) :: globalRMS(:) !<globalRMS(derivativeIdx). On return, the global RMS error
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: derivativeIdx,maximumNumberOfDerivatives,nodeIdx,numberOfDerivatives,numberOfNodes,totalNumberOfNodes
    REAL(DP) :: errorValue
    INTEGER(INTG) :: ghostNumber(MAXIMUM_GLOBAL_DERIV_NUMBER),number(MAXIMUM_GLOBAL_DERIV_NUMBER),mpiIerror, &
      & numberOfGroupComputationNodes,myGroupComputationNodeNumber,groupCommunicator
    REAL(DP) :: rmsError(MAXIMUM_GLOBAL_DERIV_NUMBER),ghostRMSError(MAXIMUM_GLOBAL_DERIV_NUMBER)
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError
    TYPE(WorkGroupType), POINTER :: workGroup
        
    ENTERS("AnalyticAnalysis_RMSErrorGetNode",err,error,*999)
   
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(field,decomposition,err,error,*999)
    NULLIFY(workGroup)
    CALL Decomposition_WorkGroupGet(decomposition,workGroup,err,error,*999)
    CALL WorkGroup_GroupCommunicatorGet(workGroup,groupCommunicator,err,error,*999)
    CALL WorkGroup_NumberOfGroupNodesGet(workGroup,numberOfGroupComputationNodes,err,error,*999)
    CALL WorkGroup_GroupNodeNumberGet(workGroup,myGroupComputationNodeNumber,err,error,*999)
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    NULLIFY(domain)
    CALL FieldVariable_ComponentDomainGet(fieldVariable,componentNumber,domain,err,error,*999)
    NULLIFY(domainTopology)
    CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(domainNodes)
    CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
    CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
    CALL DomainNodes_TotalNumberOfNodesGet(domainNodes,totalNumberOfNodes,err,error,*999)
    CALL DomainNodes_MaximumNumberOfDerivativesGet(domainNodes,maximumNumberOfDerivatives,err,error,*999)
    IF(SIZE(localRMS,1)<maximumNumberOfDerivatives) THEN
      localError="The size of the local RMS array of "//TRIM(NumberToVString(SIZE(localRMS,1),"*",err,error))// &
        & " is too small. The size of the array must be >= "//TRIM(NumberToVString(maximumNumberOfDerivatives,"*",err,error))// &
        & ", the maximum number of derivatives for the nodes."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(localGhostRMS,1)<maximumNumberOfDerivatives) THEN
      localError="The size of the local ghost RMS array of "//TRIM(NumberToVString(SIZE(localGhostRMS,1),"*",err,error))// &
        & " is too small. The size of the array must be >= "//TRIM(NumberToVString(maximumNumberOfDerivatives,"*",err,error))// &
        & ", the maximum number of derivatives for the nodes."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(globalRMS,1)<maximumNumberOfDerivatives) THEN
      localError="The size of the global RMS array of "//TRIM(NumberToVString(SIZE(globalRMS,1),"*",err,error))// &
        & " is too small. The size of the array must be >= "//TRIM(NumberToVString(maximumNumberOfDerivatives,"*",err,error))// &
        & ", the maximum number of derivatives for the nodes."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    number=0
    rmsError=0.0_DP
    ghostNumber=0
    ghostRMSError=0.0_DP
    DO nodeIdx=1,numberOfNodes
      CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfDerivatives,err,error,*999)
      DO derivativeIdx=1,numberOfDerivatives
        SELECT CASE(errorType)
        CASE(ANALYTIC_ABSOLUTE_ERROR_TYPE)
          !Default to version 1 of each node derivative
          CALL AnalyticAnalysis_AbsoluteErrorGetNode(field,variableType,1,derivativeIdx,nodeIdx,componentNumber, &
            & errorValue,err,error,*999)
        CASE(ANALYTIC_PERCENTAGE_ERROR_TYPE)
          !Default to version 1 of each node derivative
          CALL AnalyticAnalysis_PercentageErrorGetNode(field,variableType,1,derivativeIdx,nodeIdx,componentNumber, &
            & errorValue,err,error,*999)
        CASE(ANALYTIC_RELATIVE_ERROR_TYPE)
          !Default to version 1 of each node derivative
          CALL AnalyticAnalysis_RelativeErrorGetNode(field,variableType,1,derivativeIdx,nodeIdx,componentNumber, &
            & errorValue,err,error,*999)
        CASE DEFAULT
          localError="The specified error type of "//TRIM(NumberToVString(errorType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        !Accumlate the RMS errors
        number(derivativeIdx)=number(derivativeIdx)+1
        rmsError(derivativeIdx)=rmsError(derivativeIdx)+errorValue*errorValue
      ENDDO !derivativeIdx
    ENDDO !nodeIdx
    DO nodeIdx=numberOfNodes+1,totalNumberOfNodes
      CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfDerivatives,err,error,*999)
      DO derivativeIdx=1,numberOfDerivatives
        SELECT CASE(errorType)
        CASE(ANALYTIC_ABSOLUTE_ERROR_TYPE)
          !Default to version 1 of each node derivative
          CALL AnalyticAnalysis_AbsoluteErrorGetNode(field,variableType,1,derivativeIdx,nodeIdx,componentNumber, &
            & errorValue,err,error,*999)
        CASE(ANALYTIC_PERCENTAGE_ERROR_TYPE)
          !Default to version 1 of each node derivative
          CALL AnalyticAnalysis_PercentageErrorGetNode(field,variableType,1,derivativeIdx,nodeIdx,componentNumber, &
            & errorValue,err,error,*999)
        CASE(ANALYTIC_RELATIVE_ERROR_TYPE)
          !Default to version 1 of each node derivative
          CALL AnalyticAnalysis_RelativeErrorGetNode(field,variableType,1,derivativeIdx,nodeIdx,componentNumber, &
            & errorValue,err,error,*999)
        CASE DEFAULT
          localError="The specified error type of "//TRIM(NumberToVString(errorType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        !Accumlate the RMS errors
        ghostNumber(derivativeIdx)=ghostNumber(derivativeIdx)+1
        ghostRMSError(derivativeIdx)=ghostRMSError(derivativeIdx)+errorValue*errorValue
      ENDDO !derivativeIdx
    ENDDO !nodeIdx

    IF(numberOfGroupComputationNodes>1) THEN
      IF(ANY(number>0)) THEN
        DO derivativeIdx=1,maximumNumberOfDerivatives
          IF(number(derivativeIdx)>0) localRMS(derivativeIdx)=SQRT(rmsError(derivativeIdx)/number(derivativeIdx))
        ENDDO !derivativeIdx
        DO derivativeIdx=1,maximumNumberOfDerivatives
          IF(number(derivativeIdx)>0) localGhostRMS(derivativeIdx)= &
            & SQRT((rmsError(derivativeIdx)+ghostRMSError(derivativeIdx))/(number(derivativeIdx)+ghostNumber(derivativeIdx)))
        ENDDO !derivativeIdx
        !Global RMS values
        !Collect the values across the ranks
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,number,maximumNumberOfDerivatives,MPI_INTEGER,MPI_SUM,groupCommunicator,mpiIerror)
        CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIerror,err,error,*999)
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,rmsError,maximumNumberOfDerivatives,MPI_DOUBLE_PRECISION,MPI_SUM,groupCommunicator, &
          & mpiIerror)
        CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIerror,err,error,*999)
        DO derivativeIdx=1,maximumNumberOfDerivatives
          IF(number(derivativeIdx)>0) globalRMS(derivativeIdx)=SQRT(rmsError(derivativeIdx)/number(derivativeIdx))
        ENDDO !derivativeIdx
      ENDIF
    ELSE
      IF(ANY(number>0)) THEN
        DO derivativeIdx=1,maximumNumberOfDerivatives
          IF(number(derivativeIdx)>0) THEN
            localRMS(derivativeIdx)=SQRT(rmsError(derivativeIdx)/number(derivativeIdx))
            globalRMS(derivativeIdx)=localRMS(derivativeIdx)
          ENDIF
        ENDDO !derivativeIdx
      ENDIF
    ENDIF
    
    EXITS("AnalyticAnalysis_RMSErrorGetNode")
    RETURN
999 ERRORSEXITS("AnalyticAnalysis_RMSErrorGetNode",err,error)
    RETURN 1
    
  END SUBROUTINE AnalyticAnalysis_RMSErrorGetNode

  !
  !================================================================================================================================
  !

  !>Get rms error value for the field
  SUBROUTINE AnalyticAnalysis_RMSErrorGetElement(field,variableType,componentNumber,errorType,localRMS,localGhostRMS, &
    & globalRMS,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the RMS error for.
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the RMS error for \see FieldRoutines_VariableTypes,FieldRoutines
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The field variable component number to get the RMS error for. 
    INTEGER(INTG), INTENT(IN) :: errorType !<The error type to get the RMS error for \see AnalyticAnalysisRoutines_ErrorTypes,AnalyticAnalysisRoutines
    REAL(DP), INTENT(OUT) :: localRMS !<On return, the local RMS error
    REAL(DP), INTENT(OUT) :: localGhostRMS !<On return, the local + ghost RMS error
    REAL(DP), INTENT(OUT) :: globalRMS !<On return, the global RMS error
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx
    INTEGER(INTG) :: ghostNumber,number,mpiIerror,numberOfGroupComputationNodes,myGroupComputationNodeNumber,groupCommunicator
    REAL(DP) :: errorValue,rmsError,ghostRMSError
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError
    TYPE(WorkGroupType), POINTER :: workGroup

    ENTERS("AnalyticAnalysis_RMSErrorGetElement",err,error,*999)

    NULLIFY(decomposition)
    CALL Field_DecompositionGet(field,decomposition,err,error,*999)
    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
    NULLIFY(workGroup)
    CALL Decomposition_WorkGroupGet(decomposition,workGroup,err,error,*999)
    CALL WorkGroup_GroupCommunicatorGet(workGroup,groupCommunicator,err,error,*999)
    CALL WorkGroup_NumberOfGroupNodesGet(workGroup,numberOfGroupComputationNodes,err,error,*999)
    CALL WorkGroup_GroupNodeNumberGet(workGroup,myGroupComputationNodeNumber,err,error,*999)
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    NULLIFY(domain)
    CALL FieldVariable_ComponentDomainGet(fieldVariable,componentNumber,domain,err,error,*999)
    NULLIFY(domainTopology)
    CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(domainElements)
    CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
    number=0
    rmsError=0.0_DP
    ghostNumber=0
    ghostRMSError=0.0_DP
    DO elementIdx=1,domainElements%numberOfElements
      SELECT CASE(errorType)
      CASE(ANALYTIC_ABSOLUTE_ERROR_TYPE)
        CALL AnalyticAnalysis_AbsoluteErrorGetElement(field,variableType,elementIdx,componentNumber,errorValue, &
          & err,error,*999)
      CASE(ANALYTIC_PERCENTAGE_ERROR_TYPE)
        CALL AnalyticAnalysis_PercentageErrorGetElement(field,variableType,elementIdx,componentNumber, &
          & errorValue,err,error,*999)
      CASE(ANALYTIC_RELATIVE_ERROR_TYPE)
        CALL AnalyticAnalysis_RelativeErrorGetElement(field,variableType,elementIdx,componentNumber,errorValue, &
          & err,error,*999)
      CASE DEFAULT
        localError="The specified error type of "//TRIM(NumberToVString(errorType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      number=number+1
      rmsError=rmsError+errorValue*errorValue
    ENDDO !elementIdx
    DO elementIdx=domainElements%numberOfElements+1,domainElements%totalNumberOfElements
      SELECT CASE(errorType)
      CASE(ANALYTIC_ABSOLUTE_ERROR_TYPE)
        CALL AnalyticAnalysis_AbsoluteErrorGetElement(field,variableType,elementIdx,componentNumber,errorValue, &
          & err,error,*999)
      CASE(ANALYTIC_PERCENTAGE_ERROR_TYPE)
        CALL AnalyticAnalysis_PercentageErrorGetElement(field,variableType,elementIdx,componentNumber, &
          & errorValue,err,error,*999)
      CASE(ANALYTIC_RELATIVE_ERROR_TYPE)
        CALL AnalyticAnalysis_RelativeErrorGetElement(field,variableType,elementIdx,componentNumber,errorValue, &
          & err,error,*999)
      CASE DEFAULT
        localError="The specified error type of "//TRIM(NumberToVString(errorType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      ghostNumber=ghostNumber+1
      ghostRMSError=ghostRMSError+errorValue*errorValue
    ENDDO !elementIdx
    IF(number>0) THEN
      IF(numberOfGroupComputationNodes>1) THEN
        !Local elements only
        localRMS=SQRT(rmsError/number)
        !Local and ghost elements
        localGhostRMS=SQRT((rmsError+ghostRMSError)/(number+ghostNumber))
        !Global RMS values
        !Collect the values across the ranks
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,number,1,MPI_INTEGER,MPI_SUM,groupCommunicator,mpiIerror)
        CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIerror,err,error,*999)
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,rmsError,1,MPI_DOUBLE_PRECISION,MPI_SUM,groupCommunicator,mpiIerror)
        CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIerror,err,error,*999)
        globalRMS=SQRT(rmsError/number)
      ENDIF
    ENDIF

    EXITS("AnalyticAnalysis_RMSErrorGetElement")
    RETURN
999 ERRORSEXITS("AnalyticAnalysis_RMSErrorGetElement",err,error)
    RETURN 1
  END SUBROUTINE AnalyticAnalysis_RMSErrorGetElement

  !
  !================================================================================================================================
  !    
 
END MODULE AnalyticAnalysisRoutines
