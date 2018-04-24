!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a 1D diffusion equation and compare it to the analytic solution using OpenCMISS calls.
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

!> \example ClassicalField/Diffusion/Analytic1DDiffusion/src/Analytic1DDiffusionExample.f90
!! Example program to solve a 1D diffusion equation using OpenCMISS calls.
!!
!! \htmlinclude ClassicalField/Diffusion/Analytic1DDiffusion/history.html
!<

!> Main program
PROGRAM Analytic1DDiffusionExample

  USE OpenCMISS
  USE OpenCMISS_Iron

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(CMISSRP), PARAMETER ::    PI=3.141592653589793238462643383279502884197_CMISSRP

  INTEGER(CMISSIntg), PARAMETER :: NUMBER_GLOBAL_X_ELEMENTS=6
  REAL(CMISSRP), PARAMETER :: LENGTH=3.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: END_TIME=0.1_CMISSRP
  REAL(CMISSRP), PARAMETER :: TIME_STEP=0.01_CMISSRP
  REAL(CMISSRP), PARAMETER :: A=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: B=PI/2.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: C=0.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: K=1.0_CMISSRP
  
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: AnalyticFieldUserNumber=13
  !Program types
  
  !Program variables
  
  !CMISS variables

  TYPE(cmfe_BasisType) :: Basis
  TYPE(cmfe_ContextType) :: context
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,DependentField,EquationsSetField,MaterialsField,AnalyticField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh  
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver, LinearSolver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: err
  
#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

  !Intialise OpenCMISS
  CALL cmfe_Context_Initialise(context,err)
  CALL cmfe_Initialise(context,Err)

  CALL cmfe_Region_Initialise(worldRegion,err)
  CALL cmfe_Context_WorldRegionGet(context,worldRegion,err)
  
  CALL cmfe_errorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,err)

  CALL cmfe_OutputSetOn("Diffusion1DAnalytic",err)
  
  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,context,CoordinateSystem,err)
  !Set the coordinate system to be 1D
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,1,err)
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,err)
  
  !Start the creation of the region
  CALL cmfe_Region_Initialise(Region,err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,err)
  !Label the Region
  CALL cmfe_Region_LabelSet(Region,"Region",err)
  !Set the regions coordinate system to the 1D RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region,err)
  
  !Start the creation of a basis (default is trilinear lagrange)
  CALL cmfe_Basis_Initialise(Basis,err)
  CALL cmfe_Basis_CreateStart(BasisUserNumber,context,Basis,err)
  !Set the basis to be a linear Lagrange basis
  CALL cmfe_Basis_NumberOfXiSet(Basis,1,err)
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(BASIS,err)
  
  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,err)
  !Set up a regular mesh
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,Basis,err)   
  !Define the mesh on the region
  CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[LENGTH],err)
  CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS],err)
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(Mesh,err)
  CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,err)
  
  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,1,err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition,err)
  
  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(GeometricField,err)
  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,err)
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,err)
       
  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,err)
  
  !Create the equations_set
  CALL cmfe_EquationsSet_Initialise(EquationsSet,err)
  CALL cmfe_Field_Initialise(EquationsSetField,err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE],EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,err)

  !Create the equations set dependent field variables
  CALL cmfe_Field_Initialise(DependentField,err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,err)

  !Create the equations set material field variables
  CALL cmfe_Field_Initialise(MaterialsField,err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,err)
  !Set the conductivity
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,K,err)
 
  !Create the equations set analytic field variables
  CALL cmfe_Field_Initialise(AnalyticField,err)
  CALL cmfe_EquationsSet_AnalyticCreateStart(EquationsSet,CMFE_EQUATIONS_SET_DIFFUSION_EQUATION_ONE_DIM_1, & 
    & AnalyticFieldUserNumber,AnalyticField,err)
  !Finish the equations set analytic field variables
  CALL cmfe_EquationsSet_AnalyticCreateFinish(EquationsSet,err)
  !Set the analytic field parameters
  !Set the multiplicative constant. 
  CALL cmfe_Field_ComponentValuesInitialise(AnalyticField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,A,err)
  !Set the phase. 
  CALL cmfe_Field_ComponentValuesInitialise(AnalyticField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,B,err)
  !Set the offset. 
  CALL cmfe_Field_ComponentValuesInitialise(AnalyticField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,C,err)
  !Set the length to be the length of the mesh
  CALL cmfe_Field_ComponentValuesInitialise(AnalyticField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,LENGTH,err)
  
  !Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,err)
  !Set the equations set output
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_TIMING_OUTPUT,err)
  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_MATRIX_OUTPUT,err)
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,err)

  !Create the problem
  CALL cmfe_Problem_Initialise(Problem,err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,context,[CMFE_PROBLEM_CLASSICAL_FIELD_CLASS, &
    & CMFE_PROBLEM_DIFFUSION_EQUATION_TYPE,CMFE_PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE],Problem,err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,err)

  !Create the problem control
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,err)
  CALL cmfe_ControlLoop_Initialise(ControlLoop,err)
  !Get the control loop
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,err)
  !Set the times
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,0.0_CMISSRP,END_TIME,TIME_STEP,err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,err)

  !Start the creation of the problem solvers
  CALL cmfe_Solver_Initialise(Solver,err)
  CALL cmfe_Solver_Initialise(LinearSolver,err)
  CALL cmfe_Problem_SolversCreateStart(Problem,err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_NO_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_TIMING_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_SOLVER_OUTPUT,err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_MATRIX_OUTPUT,err)
  CALL cmfe_Solver_DynamicLinearSolverGet(Solver,LinearSolver,err)
  CALL cmfe_Solver_OutputTypeSet(LinearSolver,CMFE_SOLVER_PROGRESS_OUTPUT,err)
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(Problem,err)

  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(Solver,err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,err)
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,err)
  !Get the solve equations
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,err)
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_FULL_MATRICES,err)  
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,err)

  !Create the equations set boundary conditions
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,err)
  CALL cmfe_SolverEquations_BoundaryConditionsAnalytic(SolverEquations,err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,err)

  !Solve the problem
  CALL cmfe_Problem_Solve(Problem,err)

  !Output Analytic analysis
  !CALL cmfe_EquationsSet_AnalyticTimeSet(EquationsSet,END_TIME,err)
  !CALL cmfe_EquationsSet_AnalyticEvaluate(EquationsSet,err)
  CALL cmfe_AnalyticAnalysis_Output(DependentField,"Diffusion1DAnalytic",err)

  !Output fields
  CALL cmfe_Fields_Initialise(Fields,err)
  CALL cmfe_Fields_Create(Region,Fields,err)
  CALL cmfe_Fields_NodesExport(Fields,"Diffusion1DAnalytic","FORTRAN",err)
  CALL cmfe_Fields_ElementsExport(Fields,"Diffusion1DAnalytic","FORTRAN",err)
  CALL cmfe_Fields_Finalise(Fields,err)

  !Finalise and quit
  CALL cmfe_Finalise(context,err)
  WRITE(*,'(A)') "Program successfully completed."

  STOP
  
END PROGRAM Analytic1DDiffusionExample
