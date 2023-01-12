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
!> Contributor(s): David Nickerson, Chris Bradley
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

!> Main program
PROGRAM CellMLIntegrationFortranExample

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  !Test program parameters

  REAL(CMISSRP), PARAMETER :: HEIGHT=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=1.0_CMISSRP

  INTEGER(CMISSIntg), PARAMETER :: CONTEXT_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=2
  INTEGER(CMISSIntg), PARAMETER :: REGION_USER_NUMBER=3
  INTEGER(CMISSIntg), PARAMETER :: BASIS_USER_NUMBER=4
  INTEGER(CMISSIntg), PARAMETER :: GENERATED_MESH_USER_NUMBER=5
  INTEGER(CMISSIntg), PARAMETER :: MESH_USER_NUMBER=6
  INTEGER(CMISSIntg), PARAMETER :: DECOMPOSITION_USER_NUMBER=7
  INTEGER(CMISSIntg), PARAMETER :: DECOMPOSER_USER_NUMBER=8
  INTEGER(CMISSIntg), PARAMETER :: GEOMETRIC_FIELD_USER_NUMBER=9
  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER=10
  INTEGER(CMISSIntg), PARAMETER :: DEPENDENT_FIELD_USER_NUMBER=11
  INTEGER(CMISSIntg), PARAMETER :: MATERIALS_FIELD_USER_NUMBER=12
  INTEGER(CMISSIntg), PARAMETER :: CELLML_USER_NUMBER=13
  INTEGER(CMISSIntg), PARAMETER :: CELLML_MODELS_FIELD_USER_NUMBER=14
  INTEGER(CMISSIntg), PARAMETER :: CELLML_STATE_FIELD_USER_NUMBER=15
  INTEGER(CMISSIntg), PARAMETER :: CELLML_INTERMEDIATE_FIELD_USER_NUMBER=16
  INTEGER(CMISSIntg), PARAMETER :: CELLML_PARAMETERS_FIELD_USER_NUMBER=17
  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=18
  INTEGER(CMISSIntg), PARAMETER :: PROBLEM_USER_NUMBER=19

  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: numberOfArguments,argumentLength,status
  CHARACTER(LEN=255) :: commandArgument,cellmlFile
  LOGICAL :: fileExist

  INTEGER(CMISSIntg) :: numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements

  LOGICAL :: exportField

  INTEGER(CMISSIntg) :: n98ModelIndex

  INTEGER(CMISSIntg) :: gNacomponent,stimComponent,nodeIdx

  REAL(CMISSRP) :: X,Y,distance,gNaValue
  
  INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_ELEMENTS=1
  
  INTEGER(CMISSIntg) :: outputFrequency = 1
  REAL(CMISSRP), PARAMETER :: STIM_VALUE = 100.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: STIM_STOP = 0.10_CMISSRP
  REAL(CMISSRP), PARAMETER :: TIME_STOP = 1.50_CMISSRP
  REAL(CMISSRP), PARAMETER :: ODE_TIME_STEP = 0.00001_CMISSRP
  REAL(CMISSRP), PARAMETER :: PDE_TIME_STEP = 0.001_CMISSRP
  REAL(CMISSRP), PARAMETER :: CONDUCTIVITY = 0.1_CMISSRP

  !CMISS variables

  TYPE(cmfe_BasisType) :: basis
  TYPE(cmfe_BoundaryConditionsType) :: boundaryConditions
  TYPE(cmfe_CellMLType) :: cellML
  TYPE(cmfe_CellMLEquationsType) :: cellMLEquations
  TYPE(cmfe_ComputationEnvironmentType) :: computationEnvironment
  TYPE(cmfe_ContextType) :: context
  TYPE(cmfe_ControlLoopType) :: controlLoop
  TYPE(cmfe_CoordinateSystemType) :: coordinateSystem
  TYPE(cmfe_DecompositionType) :: decomposition
  TYPE(cmfe_DecomposerType) :: decomposer
  TYPE(cmfe_EquationsType) :: equations
  TYPE(cmfe_EquationsSetType) :: equationsSet
  TYPE(cmfe_FieldType) :: geometricField,equationsSetField,dependentField,materialsField
  TYPE(cmfe_FieldType) :: cellMLModelsField,cellMLStateField,cellMLIntermediateField,cellMLParametersField
  TYPE(cmfe_FieldsType) :: fields
  TYPE(cmfe_GeneratedMeshType) :: generatedMesh  
  TYPE(cmfe_MeshType) :: mesh
  TYPE(cmfe_ProblemType) :: problem
  TYPE(cmfe_RegionType) :: region,worldRegion
  TYPE(cmfe_SolverType) :: solver
  TYPE(cmfe_SolverEquationsType) :: solverEquations
  TYPE(cmfe_WorkGroupType) :: worldWorkGroup

  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: numberOfComputationNodes,computationNodeNumber
  INTEGER(CMISSIntg) :: decompositionIndex,equationsSetIndex,cellMLIndex
  INTEGER(CMISSIntg) :: firstNodeNumber,lastNodeNumber
  INTEGER(CMISSIntg) :: firstNodeDomain,lastNodeDomain,nodeDomain
  INTEGER(CMISSIntg) :: err


  !Process command line arguments before getting started.
  numberOfArguments = COMMAND_ARGUMENT_COUNT()
  !We must at least have a CellML File specified
  IF(numberOfArguments >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1,commandArgument,argumentLength,status)
    cellmlFile = ADJUSTL(commandArgument)
    WRITE(*,'("CellML File: ", A)') cellmlFile
    INQUIRE(FILE=cellmlFile,EXIST=fileExist)
    IF(.NOT.fileExist) THEN
      WRITE(*,'(">>ERROR: File does not exist.")')
      STOP
    ENDIF
  ELSE
    WRITE(*,'(">>USAGE: ",A)') "FortranExample <CellML Model URL>"
    STOP
  ENDIF

  !Intialise OpenCMISS
  CALL cmfe_Initialise(err)
  !Trap errors
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,err)
  !Create a context
  CALL cmfe_Context_Initialise(context,err)  
  CALL cmfe_Context_Create(CONTEXT_USER_NUMBER,context,err)  
  CALL cmfe_Region_Initialise(worldRegion,err)
  CALL cmfe_Context_WorldRegionGet(context,worldRegion,err)
  
  !Get the computation nodes information
  CALL cmfe_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL cmfe_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL cmfe_WorkGroup_Initialise(worldWorkGroup,err)
  CALL cmfe_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL cmfe_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationNodes,err)
  CALL cmfe_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationNodeNumber,err)

  !CALL cmfe_OutputSetOn("Monodomain",err)
    
  numberOfGlobalXElements=NUMBER_OF_ELEMENTS
  numberOfGlobalYElements=0
  numberOfGlobalZElements=0
  
  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL cmfe_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
  !Set the coordinate system to be 1D
  CALL cmfe_CoordinateSystem_DimensionSet(coordinateSystem,1,err)
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(coordinateSystem,err)

  !Start the creation of the region
  CALL cmfe_Region_Initialise(region,err)
  CALL cmfe_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
  !Set the regions coordinate system to the RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(region,coordinateSystem,err)
  !Set the region label
  CALL cmfe_Region_LabelSet(region,"Region",err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(region,err)

  !Start the creation of a basis (default is trilinear lagrange)
  CALL cmfe_Basis_Initialise(basis,err)
  CALL cmfe_Basis_CreateStart(BASIS_USER_NUMBER,context,basis,err)
  !Set the basis to be a bilinear Lagrange basis
  CALL cmfe_Basis_NumberOfXiSet(basis,1,err)
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(basis,err)

  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(generatedMesh,err)
  CALL cmfe_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(generatedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(generatedMesh,basis,err)   
  !Define the mesh on the region
  CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH],err)
  CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements],err)
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(mesh,err)
  CALL cmfe_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh,err)

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(decomposition,err)
  CALL cmfe_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(decomposition,err)

  !Decompose
  CALL cmfe_Decomposer_Initialise(decomposer,err)
  CALL cmfe_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL cmfe_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL cmfe_Decomposer_CreateFinish(decomposer,err)
  
  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(geometricField,err)
  CALL cmfe_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region,geometricField,err)
  !Set the decomposition to use
  CALL cmfe_Field_DecompositionSet(geometricField,decomposition,err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,err)
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(geometricField,err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)
        
  !Create the equations_set
  CALL cmfe_EquationsSet_Initialise(equationsSet,err)
  CALL cmfe_Field_Initialise(equationsSetField,err)
  !Set the equations set to be a Monodomain equations set - but we won't actually use it?
  CALL cmfe_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,geometricField,[CMFE_EQUATIONS_SET_BIOELECTRICS_CLASS, &
    & CMFE_EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE,CMFE_EQUATIONS_SET_NO_SUBTYPE],EQUATIONS_SET_FIELD_USER_NUMBER, &
    & equationsSetField,equationsSet,err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(equationsSet,err)

  !Create the equations set dependent field variables
  CALL cmfe_Field_Initialise(dependentField,err)
  CALL cmfe_EquationsSet_DependentCreateStart(equationsSet,DEPENDENT_FIELD_USER_NUMBER,dependentField,err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(equationsSet,err)
  
  !Create the equations set materials field variables
  CALL cmfe_Field_Initialise(materialsField,err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(equationsSet,MATERIALS_FIELD_USER_NUMBER,materialsField,err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(equationsSet,err)
  
  !Set Am
  CALL cmfe_Field_ComponentValuesInitialise(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
    & 193.6_CMISSRP,err)
  !Set Cm
  CALL cmfe_Field_ComponentValuesInitialise(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
    & 0.014651_CMISSRP,err)
  !Set conductivity
  CALL cmfe_Field_ComponentValuesInitialise(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
    & CONDUCTIVITY,err)

  !Create the CellML environment
  CALL cmfe_CellML_Initialise(cellML,err)
  CALL cmfe_CellML_CreateStart(CELLML_USER_NUMBER,region,cellML,err)
  !Import a Noble 1998 model from a file
  CALL cmfe_CellML_ModelImport(cellML,cellmlFile,n98ModelIndex,err)
  CALL cmfe_CellML_VariableSetAsKnown(cellML,n98ModelIndex,"fast_sodium_current/g_Na ",err)
  CALL cmfe_CellML_VariableSetAsKnown(cellML,n98ModelIndex,"membrane/IStim",err)
  CALL cmfe_CellML_VariableSetAsWanted(cellML,n98ModelIndex,"membrane/i_K1",err)
  CALL cmfe_CellML_VariableSetAsWanted(cellML,n98ModelIndex,"membrane/i_to",err)
  CALL cmfe_CellML_VariableSetAsWanted(cellML,n98ModelIndex,"membrane/i_K",err)
  CALL cmfe_CellML_VariableSetAsWanted(cellML,n98ModelIndex,"membrane/i_K_ATP",err)
  CALL cmfe_CellML_VariableSetAsWanted(cellML,n98ModelIndex,"membrane/i_Ca_L_K",err)
  CALL cmfe_CellML_VariableSetAsWanted(cellML,n98ModelIndex,"membrane/i_b_K",err)
  CALL cmfe_CellML_VariableSetAsWanted(cellML,n98ModelIndex,"membrane/i_NaK",err)
  CALL cmfe_CellML_VariableSetAsWanted(cellML,n98ModelIndex,"membrane/i_Na",err)
  CALL cmfe_CellML_VariableSetAsWanted(cellML,n98ModelIndex,"membrane/i_b_Na",err)
  CALL cmfe_CellML_VariableSetAsWanted(cellML,n98ModelIndex,"membrane/i_Ca_L_Na",err)
  CALL cmfe_CellML_VariableSetAsWanted(cellML,n98ModelIndex,"membrane/i_NaCa",err)
  !Finish the CellML environment
  CALL cmfe_CellML_CreateFinish(cellML,err)

  !Start the creation of CellML <--> OpenCMISS field maps
  CALL cmfe_CellML_FieldMapsCreateStart(cellML,err)
  !Now we can set up the field variable component <--> CellML model variable mappings.
  !Map Vm
  CALL cmfe_CellML_CreateFieldToCellMLMap(cellML,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & n98ModelIndex,"membrane/V",CMFE_FIELD_VALUES_SET_TYPE,err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(cellML,n98ModelIndex,"membrane/V",CMFE_FIELD_VALUES_SET_TYPE, &
    & dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,err)
  !Finish the creation of CellML <--> OpenCMISS field maps
  CALL cmfe_CellML_FieldMapsCreateFinish(cellML,err)

  !todo - get vm initialial value.
  CALL cmfe_Field_ComponentValuesInitialise(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
    & -92.5_CMISSRP,err)
  
  !Start the creation of the CellML models field
  CALL cmfe_Field_Initialise(cellMLModelsField,err)
  CALL cmfe_CellML_ModelsFieldCreateStart(cellML,CELLML_MODELS_FIELD_USER_NUMBER,cellMLModelsField,err)
  !Finish the creation of the CellML models field
  CALL cmfe_CellML_ModelsFieldCreateFinish(cellML,err)

  !Start the creation of the CellML state field
  CALL cmfe_Field_Initialise(cellMLStateField,err)
  CALL cmfe_CellML_StateFieldCreateStart(cellML,CELLML_STATE_FIELD_USER_NUMBER,cellMLStateField,err)
  !Finish the creation of the CellML state field
  CALL cmfe_CellML_StateFieldCreateFinish(cellML,err)

  !Start the creation of the CellML intermediate field
  CALL cmfe_Field_Initialise(cellMLIntermediateField,err)
  CALL cmfe_CellML_IntermediateFieldCreateStart(cellML,CELLML_INTERMEDIATE_FIELD_USER_NUMBER,cellMLIntermediateField,err)
  !Finish the creation of the CellML intermediate field
  CALL cmfe_CellML_IntermediateFieldCreateFinish(cellML,err)
  
  !Start the creation of CellML parameters field
  CALL cmfe_Field_Initialise(cellMLParametersField,err)
  CALL cmfe_CellML_ParametersFieldCreateStart(cellML,CELLML_PARAMETERS_FIELD_USER_NUMBER,cellMLParametersField,err)
  !Finish the creation of CellML parameters
  CALL cmfe_CellML_ParametersFieldCreateFinish(cellML,err)
  
  !Create the equations set equations
  CALL cmfe_Equations_Initialise(equations,err)
  CALL cmfe_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(equations,CMFE_EQUATIONS_SPARSE_MATRICES,err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_NO_OUTPUT,err)
  !CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_TIMING_OUTPUT,err)
  !CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_MATRIX_OUTPUT,err)
  !CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(equationsSet,err)

  !Find the domains of the first and last nodes
  firstNodeNumber=1
  IF(numberOfGlobalZElements==0) THEN
    lastNodeNumber=(numberOfGlobalXElements+1)*(numberOfGlobalYElements+1)
  ELSE
    lastNodeNumber=(numberOfGlobalXElements+1)*(numberOfGlobalYElements+1)*(numberOfGlobalZElements+1)
  ENDIF
  CALL cmfe_Decomposition_NodeDomainGet(decomposition,firstNodeNumber,1,firstNodeDomain,err)
  CALL cmfe_Decomposition_NodeDomainGet(decomposition,lastNodeNumber,1,lastNodeDomain,err)

  CALL cmfe_CellML_FieldComponentGet(cellML,n98ModelIndex,CMFE_CELLML_PARAMETERS_FIELD,"membrane/IStim",stimComponent,err)
  !Set the Stimulus all nodes?
  DO nodeIdx=1,NUMBER_OF_ELEMENTS+1
    CALL cmfe_Decomposition_NodeDomainGet(decomposition,nodeIdx,1,nodeDomain,err)
    IF(nodeDomain==computationNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(cellMLParametersField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
        & nodeIdx,stimComponent,STIM_VALUE,err)
    ENDIF
  ENDDO

!  !Set up the g_Na gradient
!  CALL cmfe_CellML_FieldComponentGet(cellML,n98ModelIndex,CMFE_CELLML_PARAMETERS_FIELD,"fast_sodium_current/g_Na", &
!    & gNacomponent,err)
!  !Loop over the nodes
!  DO nodeIdx=1,lastNodeNumber
!    CALL cmfe_Decomposition_NodeDomainGet(decomposition,nodeIdx,1,nodeDomain,err)
!    IF(nodeDomain==computationNodeNumber) THEN
!      CALL cmfe_Field_ParameterSetGetNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,nodeIdx,1, &
!        & X,err)
!      CALL cmfe_Field_ParameterSetGetNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,nodeIdx,2, &
!        & Y,err)
!      distance=SQRT(X**2+Y**2)/SQRT(2.0_CMISSRP)
!      gNaValue=2.0_CMISSRP*(distance+0.5_CMISSRP)*385.5e-3_CMISSRP
!      CALL cmfe_Field_ParameterSetUpdateNode(cellMLParametersField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
!        & nodeIdx,gNacomponent,gNaValue,err)
!    ENDIF
!  ENDDO

  !Start the creation of a problem.
  CALL cmfe_Problem_Initialise(problem,err)
  !Set the problem to be a standard Monodomain problem
  CALL cmfe_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[CMFE_PROBLEM_BIOELECTRICS_CLASS, &
    & CMFE_PROBLEM_MONODOMAIN_EQUATION_TYPE,CMFE_PROBLEM_MONODOMAIN_GUDUNOV_SPLIT_SUBTYPE],problem,err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(problem,err)

  !Start the creation of the problem control loop
  !Loop in time for STIM_STOP with the Stimulus applied.
  CALL cmfe_Problem_ControlLoopCreateStart(problem,err)
  !Get the control loop
  CALL cmfe_ControlLoop_Initialise(controlLoop,err)
  CALL cmfe_Problem_ControlLoopGet(problem,CMFE_CONTROL_LOOP_NODE,controlLoop,err)
  !Set the times
  CALL cmfe_ControlLoop_TimesSet(controlLoop,0.0_CMISSRP,STIM_STOP,PDE_TIME_STEP,err)
  !Set the output
  CALL cmfe_ControlLoop_OutputTypeSet(controlLoop,CMFE_CONTROL_LOOP_TIMING_OUTPUT,err)
  !Set the output frequency (0 for no output, n for output every n time steps)
  CALL cmfe_ControlLoop_TimeOutputSet(controlLoop,outputFrequency,err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(problem,err)

  !Start the creation of the problem solvers
  CALL cmfe_Problem_SolversCreateStart(problem,err)
  !Get the first (DAE) solver
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
  !Set the DAE time step
  CALL cmfe_Solver_DAETimeStepSet(solver,ODE_TIME_STEP,err)
  !CALL cmfe_Solver_DAESolverTypeSet(solver,CMFE_SOLVER_DAE_EXTERNAL,err)
  CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_NO_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_PROGRESS_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_TIMING_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_SOLVER_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_MATRIX_OUTPUT,err)
  !Get the second (Parabolic) solver
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,2,solver,err)
  CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_NO_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_PROGRESS_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_TIMING_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_SOLVER_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_MATRIX_OUTPUT,err)
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(problem,err)

  !Start the creation of the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateStart(problem,err)
  !Get the first solver
  !Get the CellML equations
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
  CALL cmfe_CellMLEquations_Initialise(cellMLEquations,err)
  CALL cmfe_Solver_CellMLEquationsGet(solver,cellMLEquations,err)
  !Add in the CellML environement
  CALL cmfe_CellMLEquations_CellMLAdd(cellMLEquations,cellML,cellMLIndex,err)
  !Finish the creation of the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateFinish(problem,err)

  !Start the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateStart(problem,err)
  !Get the second solver
  !Get the solver equations
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,2,solver,err)
  CALL cmfe_SolverEquations_Initialise(solverEquations,err)
  CALL cmfe_Solver_SolverEquationsGet(solver,solverEquations,err)
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(solverEquations,CMFE_SOLVER_SPARSE_MATRICES,err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(solverEquations,CMFE_SOLVER_FULL_MATRICES,err)
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(problem,err)

  !Start the creation of the equations set boundary conditions
  CALL cmfe_BoundaryConditions_Initialise(boundaryConditions,err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)
 !Set the first node to 0.0 and the last node to 1.0
  IF(firstNodeDomain==computationNodeNumber) THEN
   !CALL cmfe_BoundaryConditions_SetNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,firstNodeNumber,1, &
   !  & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
  ENDIF
  IF(lastNodeDomain==computationNodeNumber) THEN
   !CALL cmfe_BoundaryConditions_SetNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,lastNodeNumber,1, &
   !  & CMFE_BOUNDARY_CONDITION_FIXED,1.0_CMISSRP,err)
  ENDIF
 !Finish the creation of the equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

  !Solve the problem for the first STIM_STOP
  CALL cmfe_Problem_Solve(problem,err)

  !Now turn the stimulus off
  !Set the Stimulus at node 1
  DO nodeIdx=1,NUMBER_OF_ELEMENTS+1
    CALL cmfe_Decomposition_NodeDomainGet(decomposition,nodeIdx,1,nodeDomain,err)
    IF(nodeDomain==computationNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(cellMLParametersField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1, &
        & nodeIdx,stimComponent,0.0_CMISSRP,err)
    ENDIF
  ENDDO !nodeIdx

  !Set the time loop from STIM_STOP to TIME_STOP
  CALL cmfe_ControlLoop_TimesSet(controlLoop,STIM_STOP,TIME_STOP,PDE_TIME_STEP,err)

  !Solve the problem for the next period
  CALL cmfe_Problem_Solve(problem,err)

  exportField=.TRUE.
  IF(exportField) THEN
    CALL cmfe_Fields_Initialise(fields,err)
    CALL cmfe_Fields_Create(region,fields,err)
    CALL cmfe_Fields_NodesExport(fields,"MonodomainExample","FORTRAN",err)
    CALL cmfe_Fields_ElementsExport(fields,"MonodomainExample","FORTRAN",err)
    CALL cmfe_Fields_Finalise(fields,err)
  ENDIF

  !Destroy the context
  CALL cmfe_Context_Destroy(context,err)
  !Finialise OpenCMISS
  CALL cmfe_Finalise(err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM CellMLIntegrationFortranExample
