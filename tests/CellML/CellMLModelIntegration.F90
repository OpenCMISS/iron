!> \file
!> \author David Nickerson
!> \brief This is an example program to integrate a CellML model using OpenCMISS.
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
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
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

!> \example cellml/model-integration/Fortran/FortranExample.F90
!! Example program to integrate a CellML model using OpenCMISS.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/Bioelectrics/Monodomain/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/Bioelectrics/Monodomain/build-gnu'>Linux GNU Build</a>
!!
!<

!> Main program
PROGRAM CellMLIntegrationFortranExample

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !Test program parameters

  REAL(CMISSRP), PARAMETER :: HEIGHT=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=1.0_CMISSRP

  INTEGER(CMISSIntg), PARAMETER :: contextUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: coordinateSystemUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: regionUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: basisUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: generatedMeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: meshUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: decompositionUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: decomposerUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: geometricFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: equationsSetFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: dependentFieldUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: materialsFieldUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: cellMLUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: cellMLModelsFieldUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: cellMLStateFieldUserNumber=15
  INTEGER(CMISSIntg), PARAMETER :: cellMLIntermediateFieldUserNumber=16
  INTEGER(CMISSIntg), PARAMETER :: cellMLParametersFieldUserNumber=17
  INTEGER(CMISSIntg), PARAMETER :: equationsSetUserNumber=18
  INTEGER(CMISSIntg), PARAMETER :: problemUserNumber=19

  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: numberOfArguments,argumentLength,status
  CHARACTER(LEN=255) :: commandArgument,cellmlFile
  LOGICAL :: fileExist

  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS

  LOGICAL :: EXPORT_FIELD

  INTEGER(CMISSIntg) :: n98ModelIndex

  INTEGER(CMISSIntg) :: gNacomponent,stimcomponent,node_idx

  REAL(CMISSRP) :: X,Y,DISTANCE,gNa_VALUE
  
  INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_ELEMENTS=1
  INTEGER(CMISSIntg) :: OUTPUT_FREQUENCY = 1
  REAL(CMISSRP), PARAMETER :: STIM_VALUE = 100.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: STIM_STOP = 0.10_CMISSRP
  REAL(CMISSRP) :: TIME_STOP = 1.50_CMISSRP
  REAL(CMISSRP), PARAMETER :: ODE_TIME_STEP = 0.00001_CMISSRP
  REAL(CMISSRP) :: PDE_TIME_STEP = 0.001_CMISSRP
  REAL(CMISSRP), PARAMETER :: CONDUCTIVITY = 0.1_CMISSRP

  !CMISS variables

  TYPE(cmfe_BasisType) :: basis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CellMLType) :: CellML
  TYPE(cmfe_CellMLEquationsType) :: CellMLEquations
  TYPE(cmfe_ComputationEnvironmentType) :: ComputationEnvironment
  TYPE(cmfe_ContextType) :: context
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_DecomposerType) :: Decomposer
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,EquationsSetField,DependentField,MaterialsField
  TYPE(cmfe_FieldType) :: CellMLModelsField,CellMLStateField,CellMLIntermediateField,CellMLParametersField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh  
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations
  TYPE(cmfe_WorkGroupType) :: worldWorkGroup

  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: numberOfComputationNodes,ComputationNodeNumber
  INTEGER(CMISSIntg) :: decompositionIndex,equationsSetIndex,cellMLIndex
  INTEGER(CMISSIntg) :: FirstNodeNumber,LastNodeNumber
  INTEGER(CMISSIntg) :: FirstNodeDomain,LastNodeDomain,NodeDomain
  INTEGER(CMISSIntg) :: err


  !Process command line arguments before getting started.
  numberOfArguments = COMMAND_ARGUMENT_COUNT()
  !We must at least have a CellML File specified
  IF (numberOfArguments >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1,commandArgument,argumentLength,status)
    cellmlFile = adjustl(commandArgument)
    WRITE(*, '("CellML File: ", A)') cellmlFile
    inquire(file=cellmlFile, exist=fileExist)
    if (.not. fileExist) then
      write(*, '(">>ERROR: File does not exist")')
      stop
    endif
  ELSE
    WRITE(*,'(">>USAGE: ",A)') "FortranExample <CellML Model URL>"
    STOP
  ENDIF

!  IF(numberOfArguments >= 3) THEN
!    CALL GET_COMMAND_ARGUMENT(1,commandArgument,argumentLength,status)
!    !IF(status>0) CALL HANDLE_ERROR("error for command argument 1.")
!    READ(commandArgument(1:argumentLength),*) PDE_TIME_STEP
!    WRITE(*, '("PDE Step Size: ", E14.7)') PDE_TIME_STEP
!    CALL GET_COMMAND_ARGUMENT(2,commandArgument,argumentLength,status)
!    READ(commandArgument(1:argumentLength),*) TIME_STOP
!    WRITE(*, '("Stop Time: ", E14.7)') TIME_STOP
!    CALL GET_COMMAND_ARGUMENT(3,commandArgument,argumentLength,status)
!    READ(commandArgument(1:argumentLength),*) OUTPUT_FREQUENCY
!    WRITE(*, '("Output Frequency: ", I10)') OUTPUT_FREQUENCY
!  ELSE
!    !If there are not enough arguments die horribly
!    WRITE(*,'(">>USAGE: ",A)') "MonodomainExample <PDE step size> <stop time> <output frequency> <CellML Model URL>"
!    STOP
!  ENDIF

  !Intialise OpenCMISS
  CALL cmfe_Initialise(err)
  !Trap errors
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,err)
  !Create a context
  CALL cmfe_Context_Initialise(context,err)  
  CALL cmfe_Context_Create(contextUserNumber,context,err)  
  CALL cmfe_Region_Initialise(worldRegion,err)
  CALL cmfe_Context_WorldRegionGet(context,worldRegion,err)
  
  !Get the computation nodes information
  CALL cmfe_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL cmfe_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL cmfe_WorkGroup_Initialise(worldWorkGroup,err)
  CALL cmfe_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL cmfe_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationNodes,err)
  CALL cmfe_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationNodeNumber,err)

!  IF (numberOfComputationNodes .gt. 2)
!    WRITE(*,'(">>NOTE: ",A)') "It doesn't make any sense to use more than 2 computation nodes for this example?"
!    STOP
!  ENDIF

  !CALL cmfe_OutputSetOn("Monodomain",err)
    
  NUMBER_GLOBAL_X_ELEMENTS=NUMBER_OF_ELEMENTS
  NUMBER_GLOBAL_Y_ELEMENTS=0
  NUMBER_GLOBAL_Z_ELEMENTS=0
  
  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,err)
  CALL cmfe_CoordinateSystem_CreateStart(coordinateSystemUserNumber,context,CoordinateSystem,err)
  !Set the coordinate system to be 1D
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,1,err)
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,err)

  !Start the creation of the region
  CALL cmfe_Region_Initialise(Region,err)
  CALL cmfe_Region_CreateStart(regionUserNumber,WorldRegion,Region,err)
  !Set the regions coordinate system to the RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,err)
  !Set the region label
  CALL cmfe_Region_LabelSet(Region,"Region",err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region,err)

  !Start the creation of a basis (default is trilinear lagrange)
  CALL cmfe_Basis_Initialise(basis,err)
  CALL cmfe_Basis_CreateStart(basisUserNumber,context,basis,err)
  !Set the basis to be a bilinear Lagrange basis
  CALL cmfe_Basis_NumberOfXiSet(basis,1,err)
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(basis,err)

  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,basis,err)   
  !Define the mesh on the region
  CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH],err)
  CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS],err)
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(Mesh,err)
  CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,err)

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(decomposition,err)
  CALL cmfe_Decomposition_CreateStart(decompositionUserNumber,mesh,decomposition,err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(decomposition,err)

  !Decompose
  CALL cmfe_Decomposer_Initialise(decomposer,err)
  CALL cmfe_Decomposer_CreateStart(decomposerUserNumber,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL cmfe_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL cmfe_Decomposer_CreateFinish(decomposer,err)
  
  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(GeometricField,err)
  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,err)
  !Set the decomposition to use
  CALL cmfe_Field_DecompositionSet(GeometricField,Decomposition,err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,err)
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,err)
        
  !Create the equations_set
  CALL cmfe_EquationsSet_Initialise(EquationsSet,err)
  CALL cmfe_Field_Initialise(EquationsSetField,err)
  !Set the equations set to be a Monodomain equations set - but we won't actually use it?
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_BIOELECTRICS_CLASS, &
    & CMFE_EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE,CMFE_EQUATIONS_SET_NO_SUBTYPE],EquationsSetFieldUserNumber,EquationsSetField, &
    & EquationsSet,err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,err)

  !Create the equations set dependent field variables
  CALL cmfe_Field_Initialise(DependentField,err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,err)
  
  !Create the equations set materials field variables
  CALL cmfe_Field_Initialise(MaterialsField,err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,err)
  
  !Set Am
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
    & 193.6_CMISSRP,err)
  !Set Cm
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
    & 0.014651_CMISSRP,err)
  !Set conductivity
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
    & CONDUCTIVITY,err)

  !Create the CellML environment
  CALL cmfe_CellML_Initialise(CellML,err)
  CALL cmfe_CellML_CreateStart(CellMLUserNumber,Region,CellML,err)
  !Import a Noble 1998 model from a file
  CALL cmfe_CellML_ModelImport(CellML,cellmlFile,n98ModelIndex,err)
  CALL cmfe_CellML_VariableSetAsKnown(CellML,n98ModelIndex,"fast_sodium_current/g_Na ",err)
  CALL cmfe_CellML_VariableSetAsKnown(CellML,n98ModelIndex,"membrane/IStim",err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_K1",err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_to",err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_K",err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_K_ATP",err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_Ca_L_K",err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_b_K",err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_NaK",err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_Na",err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_b_Na",err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_Ca_L_Na",err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_NaCa",err)
  !Finish the CellML environment
  CALL cmfe_CellML_CreateFinish(CellML,err)

  !Start the creation of CellML <--> OpenCMISS field maps
  CALL cmfe_CellML_FieldMapsCreateStart(CellML,err)
  !Now we can set up the field variable component <--> CellML model variable mappings.
  !Map Vm
  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & n98ModelIndex,"membrane/V",CMFE_FIELD_VALUES_SET_TYPE,err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,n98ModelIndex,"membrane/V",CMFE_FIELD_VALUES_SET_TYPE, &
    & DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,err)
  !Finish the creation of CellML <--> OpenCMISS field maps
  CALL cmfe_CellML_FieldMapsCreateFinish(CellML,err)

  !todo - get vm initialial value.
  CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
    & -92.5_CMISSRP,err)
  
  !Start the creation of the CellML models field
  CALL cmfe_Field_Initialise(CellMLModelsField,err)
  CALL cmfe_CellML_ModelsFieldCreateStart(CellML,CellMLModelsFieldUserNumber,CellMLModelsField,err)
  !Finish the creation of the CellML models field
  CALL cmfe_CellML_ModelsFieldCreateFinish(CellML,err)

  !Start the creation of the CellML state field
  CALL cmfe_Field_Initialise(CellMLStateField,err)
  CALL cmfe_CellML_StateFieldCreateStart(CellML,CellMLStateFieldUserNumber,CellMLStateField,err)
  !Finish the creation of the CellML state field
  CALL cmfe_CellML_StateFieldCreateFinish(CellML,err)

  !Start the creation of the CellML intermediate field
  CALL cmfe_Field_Initialise(CellMLIntermediateField,err)
  CALL cmfe_CellML_IntermediateFieldCreateStart(CellML,CellMLIntermediateFieldUserNumber,CellMLIntermediateField,err)
  !Finish the creation of the CellML intermediate field
  CALL cmfe_CellML_IntermediateFieldCreateFinish(CellML,err)
  
  !Start the creation of CellML parameters field
  CALL cmfe_Field_Initialise(CellMLParametersField,err)
  CALL cmfe_CellML_ParametersFieldCreateStart(CellML,CellMLParametersFieldUserNumber,CellMLParametersField,err)
  !Finish the creation of CellML parameters
  CALL cmfe_CellML_ParametersFieldCreateFinish(CellML,err)
  
  !Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_TIMING_OUTPUT,err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_MATRIX_OUTPUT,err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,err)

  !Find the domains of the first and last nodes
  FirstNodeNumber=1
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)
  ELSE
    LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
  ENDIF
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,FirstNodeNumber,1,FirstNodeDomain,err)
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,LastNodeNumber,1,LastNodeDomain,err)

  CALL cmfe_CellML_FieldComponentGet(CellML,n98ModelIndex,CMFE_CELLML_PARAMETERS_FIELD,"membrane/IStim",stimcomponent,err)
  !Set the Stimulus all nodes?
  DO node_idx=1,NUMBER_OF_ELEMENTS+1
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,err)
    IF(NodeDomain==ComputationNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(CellMLParametersField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
        & node_idx,stimcomponent,STIM_VALUE,err)
    ENDIF
  ENDDO

!  !Set up the g_Na gradient
!  CALL cmfe_CellML_FieldComponentGet(CellML,n98ModelIndex,CMFE_CELLML_PARAMETERS_FIELD,"fast_sodium_current/g_Na", &
!    & gNacomponent,err)
!  !Loop over the nodes
!  DO node_idx=1,LastNodeNumber
!    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,err)
!    IF(NodeDomain==ComputationNodeNumber) THEN
!      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,node_idx,1, &
!        & X,err)
!      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,node_idx,2, &
!        & Y,err)
!      DISTANCE=SQRT(X**2+Y**2)/SQRT(2.0_CMISSRP)
!      gNa_VALUE=2.0_CMISSRP*(DISTANCE+0.5_CMISSRP)*385.5e-3_CMISSRP
!      CALL cmfe_Field_ParameterSetUpdateNode(CellMLParametersField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
!        & node_idx,gNacomponent,gNa_VALUE,err)
!    ENDIF
!  ENDDO

  !Start the creation of a problem.
  CALL cmfe_Problem_Initialise(Problem,err)
  !Set the problem to be a standard Monodomain problem
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,context,[CMFE_PROBLEM_BIOELECTRICS_CLASS,CMFE_PROBLEM_MONODOMAIN_EQUATION_TYPE, &
    & CMFE_PROBLEM_MONODOMAIN_GUDUNOV_SPLIT_SUBTYPE],Problem,err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,err)

  !Start the creation of the problem control loop
  !Loop in time for STIM_STOP with the Stimulus applied.
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,err)
  !Get the control loop
  CALL cmfe_ControlLoop_Initialise(ControlLoop,err)
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,err)
  !Set the times
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,0.0_CMISSRP,STIM_STOP,PDE_TIME_STEP,err)
  !Set the output
  CALL cmfe_ControlLoop_OutputTypeSet(ControlLoop,CMFE_CONTROL_LOOP_TIMING_OUTPUT,err)
  !Set the output frequency (0 for no output, n for output every n time steps)
  CALL cmfe_ControlLoop_TimeOutputSet(ControlLoop,OUTPUT_FREQUENCY,err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,err)

  !Start the creation of the problem solvers
  CALL cmfe_Problem_SolversCreateStart(Problem,err)
  !Get the first (DAE) solver
  CALL cmfe_Solver_Initialise(Solver,err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,err)
  !Set the DAE time step
  CALL cmfe_Solver_DAETimeStepSet(Solver,ODE_TIME_STEP,err)
  !CALL cmfe_Solver_DAESolverTypeSet(Solver,CMFE_SOLVER_DAE_EXTERNAL,err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_NO_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_TIMING_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_SOLVER_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_MATRIX_OUTPUT,err)
  !Get the second (Parabolic) solver
  CALL cmfe_Solver_Initialise(Solver,err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,2,Solver,err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_NO_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_TIMING_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_SOLVER_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_MATRIX_OUTPUT,err)
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(Problem,err)

  !Start the creation of the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateStart(Problem,err)
  !Get the first solver
  !Get the CellML equations
  CALL cmfe_Solver_Initialise(Solver,err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,err)
  CALL cmfe_CellMLEquations_Initialise(CellMLEquations,err)
  CALL cmfe_Solver_CellMLEquationsGet(Solver,CellMLEquations,err)
  !Add in the CellML environement
  CALL cmfe_CellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,err)
  !Finish the creation of the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateFinish(Problem,err)

  !Start the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,err)
  !Get the second solver
  !Get the solver equations
  CALL cmfe_Solver_Initialise(Solver,err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,2,Solver,err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,err)
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_FULL_MATRICES,err)
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,err)

  !Start the creation of the equations set boundary conditions
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,err)
 !Set the first node to 0.0 and the last node to 1.0
  IF(FirstNodeDomain==ComputationNodeNumber) THEN
   !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,FirstNodeNumber,1, &
   !  & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
  ENDIF
  IF(LastNodeDomain==ComputationNodeNumber) THEN
   !CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,LastNodeNumber,1, &
   !  & CMFE_BOUNDARY_CONDITION_FIXED,1.0_CMISSRP,err)
  ENDIF
 !Finish the creation of the equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,err)

  !Solve the problem for the first STIM_STOP
  CALL cmfe_Problem_Solve(Problem,err)

  !Now turn the stimulus off
  !Set the Stimulus at node 1
  DO node_idx=1,NUMBER_OF_ELEMENTS+1
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,err)
    IF(NodeDomain==ComputationNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(CellMLParametersField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
        & node_idx,stimcomponent,0.0_CMISSRP,err)
    ENDIF
  ENDDO !node_idx

  !Set the time loop from STIM_STOP to TIME_STOP
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,STIM_STOP,TIME_STOP,PDE_TIME_STEP,err)

  !Solve the problem for the next period
  CALL cmfe_Problem_Solve(Problem,err)

  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL cmfe_Fields_Initialise(Fields,err)
    CALL cmfe_Fields_Create(Region,Fields,err)
    CALL cmfe_Fields_NodesExport(Fields,"MonodomainExample","FORTRAN",err)
    CALL cmfe_Fields_ElementsExport(Fields,"MonodomainExample","FORTRAN",err)
    CALL cmfe_Fields_Finalise(Fields,err)
  ENDIF

  !Destroy the context
  CALL cmfe_Context_Destroy(context,err)
  !Finialise OpenCMISS
  CALL cmfe_Finalise(err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM CellMLIntegrationFortranExample
