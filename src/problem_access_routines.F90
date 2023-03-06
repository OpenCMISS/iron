!> \file
!> \author Chris Bradley
!> \brief This module contains all problem access method routines.
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

!> This module contains all problem access method routines.
MODULE ProblemAccessRoutines
  
  USE BaseRoutines
  USE ControlLoopAccessRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE SolverAccessRoutines
  USE SolverMappingAccessRoutines
  USE Strings
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Problem Classes
  INTEGER(INTG), PARAMETER :: PROBLEM_NO_CLASS=0  
  INTEGER(INTG), PARAMETER :: PROBLEM_ELASTICITY_CLASS=1
  INTEGER(INTG), PARAMETER :: PROBLEM_FLUID_MECHANICS_CLASS=2
  INTEGER(INTG), PARAMETER :: PROBLEM_ELECTROMAGNETICS_CLASS=3
  INTEGER(INTG), PARAMETER :: PROBLEM_CLASSICAL_FIELD_CLASS=4  
  INTEGER(INTG), PARAMETER :: PROBLEM_BIOELECTRICS_CLASS=5
  INTEGER(INTG), PARAMETER :: PROBLEM_MODAL_CLASS=6
  INTEGER(INTG), PARAMETER :: PROBLEM_FITTING_CLASS=7
  INTEGER(INTG), PARAMETER :: PROBLEM_OPTIMISATION_CLASS=8
  INTEGER(INTG), PARAMETER :: PROBLEM_MULTI_PHYSICS_CLASS=9

  !Problem types
  INTEGER(INTG), PARAMETER :: PROBLEM_NO_TYPE=0
  !Elasticity class
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_ELASTICITY_TYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_FINITE_ELASTICITY_TYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_ELASTICITY_CONTACT_TYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE=4
  !Fluid mechanics class
  INTEGER(INTG), PARAMETER :: PROBLEM_STOKES_EQUATION_TYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_NAVIER_STOKES_EQUATION_TYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_DARCY_EQUATION_TYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_DARCY_PRESSURE_EQUATION_TYPE=4
  INTEGER(INTG), PARAMETER :: PROBLEM_POISEUILLE_EQUATION_TYPE=5
  INTEGER(INTG), PARAMETER :: PROBLEM_BURGERS_EQUATION_TYPE=6
  !Electromagnetics class
  INTEGER(INTG), PARAMETER :: PROBLEM_ELECTROSTATIC_TYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_MAGNETOSTATIC_TYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_MAXWELLS_EQUATIONS_TYPE=3
  !Classical field class
  INTEGER(INTG), PARAMETER :: PROBLEM_LAPLACE_EQUATION_TYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_POISSON_EQUATION_TYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_HELMHOLTZ_EQUATION_TYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_WAVE_EQUATION_TYPE=4
  INTEGER(INTG), PARAMETER :: PROBLEM_ADVECTION_EQUATION_TYPE=5
  INTEGER(INTG), PARAMETER :: PROBLEM_DIFFUSION_EQUATION_TYPE=6
  INTEGER(INTG), PARAMETER :: PROBLEM_ADVECTION_DIFFUSION_EQUATION_TYPE=7
  INTEGER(INTG), PARAMETER :: PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE=8
  INTEGER(INTG), PARAMETER :: PROBLEM_BIHARMONIC_EQUATION_TYPE=9
  INTEGER(INTG), PARAMETER :: PROBLEM_FITTING_TYPE=10
  INTEGER(INTG), PARAMETER :: PROBLEM_HJ_EQUATION_TYPE=11
  
  !Bioelectrics class
  INTEGER(INTG), PARAMETER :: PROBLEM_MONODOMAIN_EQUATION_TYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_BIDOMAIN_EQUATION_TYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE=3
  !Modal class
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_ELASTIC_MODAL_TYPE=1
  !Multi physics class
  INTEGER(INTG), PARAMETER :: PROBLEM_FINITE_ELASTICITY_DARCY_TYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_FINITE_ELASTICITY_STOKES_TYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_DIFFUSION_DIFFUSION_TYPE=4
  INTEGER(INTG), PARAMETER :: PROBLEM_DIFFUSION_ADVECTION_DIFFUSION_TYPE=5
  INTEGER(INTG), PARAMETER :: PROBLEM_MULTI_COMPARTMENT_TRANSPORT_TYPE=6 !<Problem type for the multi-compartment coupled transport, comprising either/or/both advection-diffusion & diffusion 
  INTEGER(INTG), PARAMETER :: PROBLEM_FINITE_ELASTICITY_FLUID_PRESSURE_TYPE=7
  INTEGER(INTG), PARAMETER :: PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE=8

  !Problem subtypes
  INTEGER(INTG), PARAMETER :: PROBLEM_NO_SUBTYPE=0
  !Elasticity class
  !  Linear elasticity
  INTEGER(INTG), PARAMETER :: PROBLEM_STOKES_DAMPING_LINEAR_ELASTICITY_SUBTYPE=1
  !  Finite elasticity
  INTEGER(INTG), PARAMETER :: PROBLEM_STATIC_FINITE_ELASTICITY_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_QUASISTATIC_FINITE_ELASTICITY_WITH_GROWTH_SUBTYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_DYNAMIC_FINITE_ELASTICITY_SUBTYPE=4
  INTEGER(INTG), PARAMETER :: PROBLEM_FINITE_ELASTICITY_WITH_ACTIVE_SUBTYPE=5
  INTEGER(INTG), PARAMETER :: PROBLEM_FINITE_ELASTICITY_WITH_CELLML_SUBTYPE=6
  INTEGER(INTG), PARAMETER :: PROBLEM_FINITE_ELASTICITY_WITH_GROWTH_CELLML_SUBTYPE=7
  INTEGER(INTG), PARAMETER :: PROBLEM_MULTISCALE_FINITE_ELASTICITY_SUBTYPE=8
  ! Linear elasticity subject to contact constraint
  INTEGER(INTG), PARAMETER :: PROBLEM_LE_CONTACT_TRANSFORM_REPROJECT_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_LE_CONTACT_TRANSFORM_SUBTYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_LE_CONTACT_REPROJECT_SUBTYPE=3
  ! Finite elasticity subject to contact constraint
  INTEGER(INTG), PARAMETER :: PROBLEM_FE_CONTACT_TRANSFORM_REPROJECT_SUBTYPE=4
  INTEGER(INTG), PARAMETER :: PROBLEM_FE_CONTACT_TRANSFORM_SUBTYPE=5
  INTEGER(INTG), PARAMETER :: PROBLEM_FE_CONTACT_REPROJECT_SUBTYPE=6

  !Fluid mechanics class
  !  Stokes equations
  INTEGER(INTG), PARAMETER :: PROBLEM_STATIC_STOKES_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_LAPLACE_STOKES_SUBTYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_TRANSIENT_STOKES_SUBTYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_OPTIMISED_STOKES_SUBTYPE=4
  INTEGER(INTG), PARAMETER :: PROBLEM_ALE_STOKES_SUBTYPE=5
  !  Navier-Stokes equations
  INTEGER(INTG), PARAMETER :: PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_OPTIMISED_NAVIER_STOKES_SUBTYPE=4
  INTEGER(INTG), PARAMETER :: PROBLEM_ALE_NAVIER_STOKES_SUBTYPE=5
  INTEGER(INTG), PARAMETER :: PROBLEM_QUASISTATIC_NAVIER_STOKES_SUBTYPE=6
  INTEGER(INTG), PARAMETER :: PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE=7
  INTEGER(INTG), PARAMETER :: PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE=8
  INTEGER(INTG), PARAMETER :: PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE=9
  INTEGER(INTG), PARAMETER :: PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE=10
  INTEGER(INTG), PARAMETER :: PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE=11
  INTEGER(INTG), PARAMETER :: PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE=12
  INTEGER(INTG), PARAMETER :: PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE=13
  INTEGER(INTG), PARAMETER :: PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE=14
  INTEGER(INTG), PARAMETER :: PROBLEM_CONSTITUTIVE_RBS_NAVIER_STOKES_SUBTYPE=15
  INTEGER(INTG), PARAMETER :: PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE=16
  !  Darcy equation
  INTEGER(INTG), PARAMETER :: PROBLEM_STANDARD_DARCY_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_QUASISTATIC_DARCY_SUBTYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_ALE_DARCY_SUBTYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_TRANSIENT_DARCY_SUBTYPE=4
  !  Poiseuille equation
  INTEGER(INTG), PARAMETER :: PROBLEM_STATIC_POISEUILLE_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_DYNAMIC_POISEUILLE_SUBTYPE=2
  !  Burgers equation
  INTEGER(INTG), PARAMETER :: PROBLEM_STATIC_BURGERS_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_DYNAMIC_BURGERS_SUBTYPE=2
  !Electromagnetics class
  !Classical field class
  !  Laplace equation
  INTEGER(INTG), PARAMETER :: PROBLEM_STANDARD_LAPLACE_SUBTYPE=1
  !  Hamilton-Jacobi equation
  INTEGER(INTG), PARAMETER :: PROBLEM_STANDARD_HJ_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_GENERALISED_HJ_SUBTYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_MOVING_MESH_HJ_SUBTYPE=3
  !  Poisson equation
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE=4
  INTEGER(INTG), PARAMETER :: PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE=5
  INTEGER(INTG), PARAMETER :: PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE=6
  INTEGER(INTG), PARAMETER :: PROBLEM_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE=7
  !  Helmholtz equation
  INTEGER(INTG), PARAMETER :: PROBLEM_STANDARD_HELMHOLTZ_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_GENERALISED_HELMHOLTZ_SUBTYPE=2
  !  Wave equation
  !  Diffusion equation
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_DIFFUSION_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_NONLINEAR_DIFFUSION_SUBTYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_ALE_DIFFUSION_SUBTYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_NONLINEAR_ALE_DIFFUSION_SUBTYPE=4
  !  Reaction-diffusion equation
  INTEGER(INTG), PARAMETER :: PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_CELLML_REAC_EVAL_REAC_DIFF_NO_SPLIT_SUBTYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_REAC_DIFF_NO_SPLIT_SUBTYPE=3
  !  Advection-diffusion equation
  INTEGER(INTG), PARAMETER :: PROBLEM_GENERALISED_ADVEC_DIFF_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_NONLINEAR_SOURCE_ADVEC_DIFF_SUBTYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE=4
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE=5
  INTEGER(INTG), PARAMETER :: PROBLEM_NONLINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE=6
  INTEGER(INTG), PARAMETER :: PROBLEM_ADVECTION_SUBTYPE=7
  !Subtypes for steady-state advection-diffusion equation
  INTEGER(INTG), PARAMETER :: PROBLEM_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE=8
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE=9
  INTEGER(INTG), PARAMETER :: PROBLEM_NONLINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE=10

  !Bioelectric class
  !  Monodomain equation
  INTEGER(INTG), PARAMETER :: PROBLEM_MONODOMAIN_GUDUNOV_SPLIT_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_MONODOMAIN_STRANG_SPLIT_SUBTYPE=2   
  INTEGER(INTG), PARAMETER :: PROBLEM_MONODOMAIN_DIRECT_MODEL_SUBTYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_MONODOMAIN_BUENOOROVIO_SUBTYPE=4
  INTEGER(INTG), PARAMETER :: PROBLEM_MONODOMAIN_TENTUSSCHER06_SUBTYPE=5
  !  Bidomain equation
  INTEGER(INTG), PARAMETER :: PROBLEM_BIDOMAIN_GUDUNOV_SPLIT_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_BIDOMAIN_STRANG_SPLIT_SUBTYPE=2

  !Fitting class
  !  Galerkin projection
  INTEGER(INTG), PARAMETER :: PROBLEM_STATIC_LINEAR_FITTING_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_STATIC_NONLINEAR_FITTING_SUBTYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_QUASISTATIC_LINEAR_FITTING_SUBTYPE=3  
  INTEGER(INTG), PARAMETER :: PROBLEM_QUASISTATIC_NONLINEAR_FITTING_SUBTYPE=4
  INTEGER(INTG), PARAMETER :: PROBLEM_DIV_FREE_VELOCITY_FITTING_SUBTYPE=5
  
  !Multi physics (subtype numbers must be different from Darcy ones)
  INTEGER(INTG), PARAMETER :: PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE=101
  INTEGER(INTG), PARAMETER :: PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE=103
  INTEGER(INTG), PARAMETER :: PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE=104
  INTEGER(INTG), PARAMETER :: PROBLEM_COUPLED_DIFFUSION_DIFFUSION_SUBTYPE=111
  INTEGER(INTG), PARAMETER :: PROBLEM_COUPLED_ALE_DIFFUSION_DIFFUSION_SUBTYPE=112
  INTEGER(INTG), PARAMETER :: PROBLEM_COUPLED_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE=121
  INTEGER(INTG), PARAMETER :: PROBLEM_COUPLED_ALE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE=122
  INTEGER(INTG), PARAMETER :: PROBLEM_STANDARD_MULTI_COMPARTMENT_TRANSPORT_SUBTYPE=131
  INTEGER(INTG), PARAMETER :: PROBLEM_STANDARD_ALE_MULTI_COMPARTMENT_TRANSPORT_SUBTYPE=132    
  INTEGER(INTG), PARAMETER :: PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE=133
  INTEGER(INTG), PARAMETER :: PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE=141
  INTEGER(INTG), PARAMETER :: PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE=142
  INTEGER(INTG), PARAMETER :: PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE=143
  INTEGER(INTG), PARAMETER :: PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE=144
  INTEGER(INTG), PARAMETER :: PROBLEM_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE=145
  INTEGER(INTG), PARAMETER :: PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE=151
  INTEGER(INTG), PARAMETER :: PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE=152
  INTEGER(INTG), PARAMETER :: PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE=153
  INTEGER(INTG), PARAMETER :: PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE=154
  INTEGER(INTG), PARAMETER :: PROBLEM_DYNAMIC_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE=155
  INTEGER(INTG), PARAMETER :: PROBLEM_DYNAMIC_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE=156

  !> \addtogroup Problem_SetupTypes Problem::Constants::SetupTypes
  !> \brief Setup type parameters
  !> \see Problem
  !>@{
  INTEGER(INTG), PARAMETER :: PROBLEM_SETUP_INITIAL_TYPE=1 !<Initial setup for a problem. \see Problem_SetupTypes,Problem
  INTEGER(INTG), PARAMETER :: PROBLEM_SETUP_CONTROL_TYPE=2 !<Solver setup for a problem. \see Problem_SetupTypes,Problem
  INTEGER(INTG), PARAMETER :: PROBLEM_SETUP_SOLVERS_TYPE=3 !<Solver setup for a problem. \see Problem_SetupTypes,Problem
  INTEGER(INTG), PARAMETER :: PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE=4 !<Solver equations setup for a problem. \see Problem_SetupTypes,Problem
  INTEGER(INTG), PARAMETER :: PROBLEM_SETUP_CELLML_EQUATIONS_TYPE=5 !<CellML equations setup for a problem. \see Problem_SetupTypes,Problem
  !>@}
  
  !> \addtogroup Problem_SetupActionTypes Problem::Constants::SetupActionTypes
  !> \brief Setup action type parameters
  !> \see Problem
  !>@{
  INTEGER(INTG), PARAMETER :: PROBLEM_SETUP_START_ACTION=1 !<Start setup action. \see Problem_SetupActionTypes,Problem
  INTEGER(INTG), PARAMETER :: PROBLEM_SETUP_FINISH_ACTION=2 !<Finish setup action. \see Problem_SetupActionTypes,Problem
  !>@}

  !> \addtogroup Problem_LinearityTypes Problem::Constants::LinearityTypes
  !> \brief Setup type parameters
  !> \see Problem
  !>@{
  INTEGER(INTG), PARAMETER :: PROBLEM_SOLVER_LINEAR=1 !<Linear problem. \see Problem_LinearityTypes,Problem
  INTEGER(INTG), PARAMETER :: PROBLEM_SOLVER_NONLINEAR=2 !<Nonlinear problem. \see Problem_LinearityTypes,Problem
  !>@}
  
  !Module types

  !Module variables
  
  !Interfaces

  INTERFACE Problem_CellMLEquationsGet
    MODULE PROCEDURE Problem_CellMLEquationsGet0
    MODULE PROCEDURE Problem_CellMLEquationsGet1
  END INTERFACE Problem_CellMLEquationsGet

  INTERFACE Problem_ControlLoopGet
    MODULE PROCEDURE Problem_ControlLoopGet0
    MODULE PROCEDURE Problem_ControlLoopGet1
  END INTERFACE Problem_ControlLoopGet

  INTERFACE Problem_SolverGet
    MODULE PROCEDURE Problem_SolverGet0
    MODULE PROCEDURE Problem_SolverGet1
  END INTERFACE Problem_SolverGet
  
  INTERFACE Problem_SolverEquationsGet
    MODULE PROCEDURE Problem_SolverEquationsGet0
    MODULE PROCEDURE Problem_SolverEquationsGet1
  END INTERFACE Problem_SolverEquationsGet

  PUBLIC PROBLEM_NO_CLASS,PROBLEM_ELASTICITY_CLASS,PROBLEM_FLUID_MECHANICS_CLASS,PROBLEM_ELECTROMAGNETICS_CLASS, &
    & PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_BIOELECTRICS_CLASS,PROBLEM_MODAL_CLASS,PROBLEM_FITTING_CLASS, &
    & PROBLEM_OPTIMISATION_CLASS,PROBLEM_MULTI_PHYSICS_CLASS

  PUBLIC PROBLEM_NO_TYPE

  PUBLIC PROBLEM_LINEAR_ELASTICITY_TYPE,PROBLEM_FINITE_ELASTICITY_TYPE,PROBLEM_LINEAR_ELASTICITY_CONTACT_TYPE, &
    & PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE

  PUBLIC PROBLEM_STOKES_EQUATION_TYPE,PROBLEM_NAVIER_STOKES_EQUATION_TYPE,PROBLEM_DARCY_EQUATION_TYPE, &
    & PROBLEM_DARCY_PRESSURE_EQUATION_TYPE,PROBLEM_POISEUILLE_EQUATION_TYPE,PROBLEM_BURGERS_EQUATION_TYPE

  PUBLIC PROBLEM_ELECTROSTATIC_TYPE,PROBLEM_MAGNETOSTATIC_TYPE,PROBLEM_MAXWELLS_EQUATIONS_TYPE

  PUBLIC PROBLEM_LAPLACE_EQUATION_TYPE,PROBLEM_POISSON_EQUATION_TYPE,PROBLEM_HELMHOLTZ_EQUATION_TYPE, &
    & PROBLEM_WAVE_EQUATION_TYPE,PROBLEM_ADVECTION_EQUATION_TYPE,PROBLEM_DIFFUSION_EQUATION_TYPE, &
    & PROBLEM_ADVECTION_DIFFUSION_EQUATION_TYPE,PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE, &
    & PROBLEM_BIHARMONIC_EQUATION_TYPE,PROBLEM_FITTING_TYPE,PROBLEM_HJ_EQUATION_TYPE

  PUBLIC PROBLEM_MONODOMAIN_EQUATION_TYPE,PROBLEM_BIDOMAIN_EQUATION_TYPE,PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE, &
    & PROBLEM_MONODOMAIN_DIRECT_MODEL_SUBTYPE,PROBLEM_MONODOMAIN_BUENOOROVIO_SUBTYPE,PROBLEM_MONODOMAIN_TENTUSSCHER06_SUBTYPE

  PUBLIC PROBLEM_LINEAR_ELASTIC_MODAL_TYPE

  PUBLIC PROBLEM_FINITE_ELASTICITY_DARCY_TYPE,PROBLEM_FINITE_ELASTICITY_STOKES_TYPE,PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE, &
    & PROBLEM_DIFFUSION_DIFFUSION_TYPE,PROBLEM_DIFFUSION_ADVECTION_DIFFUSION_TYPE,PROBLEM_MULTI_COMPARTMENT_TRANSPORT_TYPE, &
    & PROBLEM_FINITE_ELASTICITY_FLUID_PRESSURE_TYPE,PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE

  PUBLIC PROBLEM_NO_SUBTYPE

  PUBLIC PROBLEM_STOKES_DAMPING_LINEAR_ELASTICITY_SUBTYPE

  PUBLIC PROBLEM_STATIC_FINITE_ELASTICITY_SUBTYPE,PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE, &
    & PROBLEM_QUASISTATIC_FINITE_ELASTICITY_WITH_GROWTH_SUBTYPE,PROBLEM_DYNAMIC_FINITE_ELASTICITY_SUBTYPE, &
    & PROBLEM_FINITE_ELASTICITY_WITH_ACTIVE_SUBTYPE,PROBLEM_FINITE_ELASTICITY_WITH_CELLML_SUBTYPE, &
    & PROBLEM_FINITE_ELASTICITY_WITH_GROWTH_CELLML_SUBTYPE,PROBLEM_MULTISCALE_FINITE_ELASTICITY_SUBTYPE

  PUBLIC PROBLEM_LE_CONTACT_TRANSFORM_REPROJECT_SUBTYPE,PROBLEM_LE_CONTACT_TRANSFORM_SUBTYPE,PROBLEM_LE_CONTACT_REPROJECT_SUBTYPE

  PUBLIC PROBLEM_FE_CONTACT_TRANSFORM_REPROJECT_SUBTYPE,PROBLEM_FE_CONTACT_TRANSFORM_SUBTYPE,PROBLEM_FE_CONTACT_REPROJECT_SUBTYPE

  PUBLIC PROBLEM_STATIC_STOKES_SUBTYPE,PROBLEM_LAPLACE_STOKES_SUBTYPE,PROBLEM_TRANSIENT_STOKES_SUBTYPE, &
    & PROBLEM_OPTIMISED_STOKES_SUBTYPE,PROBLEM_ALE_STOKES_SUBTYPE

  PUBLIC PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE,PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE,PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
    & PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE,PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE, &
    & PROBLEM_OPTIMISED_NAVIER_STOKES_SUBTYPE,PROBLEM_ALE_NAVIER_STOKES_SUBTYPE, &
    & PROBLEM_QUASISTATIC_NAVIER_STOKES_SUBTYPE,PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
    & PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE,PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
    & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE,PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
    & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE,PROBLEM_CONSTITUTIVE_RBS_NAVIER_STOKES_SUBTYPE, &
    & PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE

  PUBLIC PROBLEM_STANDARD_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_DARCY_SUBTYPE,PROBLEM_ALE_DARCY_SUBTYPE, &
    & PROBLEM_TRANSIENT_DARCY_SUBTYPE

  PUBLIC PROBLEM_STATIC_POISEUILLE_SUBTYPE,PROBLEM_DYNAMIC_POISEUILLE_SUBTYPE

  PUBLIC PROBLEM_STATIC_BURGERS_SUBTYPE,PROBLEM_DYNAMIC_BURGERS_SUBTYPE

  PUBLIC PROBLEM_STANDARD_LAPLACE_SUBTYPE

  PUBLIC PROBLEM_STANDARD_HJ_SUBTYPE,PROBLEM_GENERALISED_HJ_SUBTYPE,PROBLEM_MOVING_MESH_HJ_SUBTYPE

  PUBLIC PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE,PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE, &
    & PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE,PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
    & PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE,PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE,PROBLEM_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE

  PUBLIC PROBLEM_STANDARD_HELMHOLTZ_SUBTYPE,PROBLEM_GENERALISED_HELMHOLTZ_SUBTYPE

  PUBLIC PROBLEM_LINEAR_DIFFUSION_SUBTYPE,PROBLEM_NONLINEAR_DIFFUSION_SUBTYPE, &
    & PROBLEM_LINEAR_ALE_DIFFUSION_SUBTYPE,PROBLEM_NONLINEAR_ALE_DIFFUSION_SUBTYPE

  PUBLIC PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE,PROBLEM_CELLML_REAC_EVAL_REAC_DIFF_NO_SPLIT_SUBTYPE, &
    & PROBLEM_REAC_DIFF_NO_SPLIT_SUBTYPE

  PUBLIC PROBLEM_GENERALISED_ADVEC_DIFF_SUBTYPE,PROBLEM_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
    & PROBLEM_NONLINEAR_SOURCE_ADVEC_DIFF_SUBTYPE,PROBLEM_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE, &
    & PROBLEM_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE,PROBLEM_NONLINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
    & PROBLEM_ADVECTION_SUBTYPE,PROBLEM_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE, &
    & PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE,PROBLEM_NONLINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE

  PUBLIC PROBLEM_MONODOMAIN_GUDUNOV_SPLIT_SUBTYPE,PROBLEM_MONODOMAIN_STRANG_SPLIT_SUBTYPE

  PUBLIC PROBLEM_BIDOMAIN_GUDUNOV_SPLIT_SUBTYPE,PROBLEM_BIDOMAIN_STRANG_SPLIT_SUBTYPE

  PUBLIC PROBLEM_STATIC_LINEAR_FITTING_SUBTYPE,PROBLEM_STATIC_NONLINEAR_FITTING_SUBTYPE, &
    & PROBLEM_QUASISTATIC_LINEAR_FITTING_SUBTYPE,PROBLEM_QUASISTATIC_NONLINEAR_FITTING_SUBTYPE, &
    & PROBLEM_DIV_FREE_VELOCITY_FITTING_SUBTYPE

  PUBLIC PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE, &
    & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE, &
    & PROBLEM_COUPLED_DIFFUSION_DIFFUSION_SUBTYPE,PROBLEM_COUPLED_ALE_DIFFUSION_DIFFUSION_SUBTYPE, &
    & PROBLEM_COUPLED_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE,PROBLEM_COUPLED_ALE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE, &
    & PROBLEM_STANDARD_MULTI_COMPARTMENT_TRANSPORT_SUBTYPE,PROBLEM_STANDARD_ALE_MULTI_COMPARTMENT_TRANSPORT_SUBTYPE, &
    & PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE,PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE, &
    & PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE,PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
    & PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE,PROBLEM_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE, &
    & PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE,PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE, &
    & PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE,PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE, &
    & PROBLEM_DYNAMIC_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE,PROBLEM_DYNAMIC_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE

  PUBLIC PROBLEM_SETUP_INITIAL_TYPE,PROBLEM_SETUP_CONTROL_TYPE,PROBLEM_SETUP_SOLVERS_TYPE,PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE, &
    & PROBLEM_SETUP_CELLML_EQUATIONS_TYPE

  PUBLIC PROBLEM_SETUP_START_ACTION,PROBLEM_SETUP_FINISH_ACTION

  PUBLIC PROBLEM_SOLVER_LINEAR,PROBLEM_SOLVER_NONLINEAR  

  PUBLIC Problem_AssertIsFinished,Problem_AssertNotFinished

  PUBLIC Problem_CellMLEquationsGet

  PUBLIC Problem_ContextGet

  PUBLIC Problem_ControlLoopGet

  PUBLIC Problem_ControlLoopRootGet

  PUBLIC Problem_Get

  PUBLIC Problem_ProblemsGet

  PUBLIC Problem_SolverGet

  PUBLIC Problem_SolverEquationsGet

  PUBLIC Problem_SpecificationGet

  PUBLIC Problem_SpecificationSizeGet

  PUBLIC Problem_UserNumberFind

  PUBLIC Problems_ContextGet

  PUBLIC Problems_ProblemGet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Assert that a problem has been finished
  SUBROUTINE Problem_AssertIsFinished(problem,err,error,*)

    !Argument Variables
    TYPE(ProblemType), POINTER, INTENT(IN) :: problem !<The problem to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Problem_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
#endif    

    IF(.NOT.problem%problemFinished) THEN
      localError="Problem number "//TRIM(NumberToVString(problem%userNumber,"*",err,error))//" has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Problem_AssertIsFinished")
    RETURN
999 ERRORSEXITS("Problem_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_AssertIsFinished

  !
  !================================================================================================================================
  !

  !>Assert that a problem has not been finished
  SUBROUTINE Problem_AssertNotFinished(problem,err,error,*)

    !Argument Variables
    TYPE(ProblemType), POINTER, INTENT(IN) :: problem !<The problem to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Problem_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
#endif    

    IF(problem%problemFinished) THEN
      localError="Problem number "//TRIM(NumberToVString(problem%userNumber,"*",err,error))//" has already been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Problem_AssertNotFinished")
    RETURN
999 ERRORSEXITS("Problem_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the CellML equations defined with a solver. \see OpenCMISS::Iron::cmfe_Problem_CellMLEquationsGet
  SUBROUTINE Problem_CellMLEquationsGet0(problem,controlLoopIdentifier,solverIndex,cellMLEquations,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to get solver CellML equations for
    INTEGER(INTG), INTENT(IN) :: controlLoopIdentifier !<The control loop identifier to get the solver CellML equations for
    INTEGER(INTG), INTENT(IN) :: solverIndex !<The solver index in the solvers to get the solver CellML equations for
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations !<On exit, a pointer to the specified solver CellML equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Problem_CellMLEquationsGet0",err,error,*999)

    CALL Problem_CellMLEquationsGet1(Problem,[controlLoopIdentifier],solverIndex,cellMLEquations,err,error,*999)
    
    EXITS("Problem_CellMLEquationsGet0")
    RETURN
999 ERRORSEXITS("Problem_CellMLEquationsGet0",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_CellMLEquationsGet0

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solver CellML equations defined with a solver. \see OPENCMISS::CMISSProblemSolverCellMLEquationsGet
  SUBROUTINE Problem_CellMLEquationsGet1(problem,controlLoopIdentifiers,solverIndex,cellMLEquations,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to get the CellML equations for
    INTEGER(INTG), INTENT(IN) :: controlLoopIdentifiers(:) !<The control loop identifier to get the CellML equations for
    INTEGER(INTG), INTENT(IN) :: solverIndex !<The solver index in the solvers to get the CellML equations for
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations !<On exit, a pointer to the specified CellML equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: controlLoop,controlLoopRoot
    TYPE(SolverType), POINTER :: solver
    TYPE(SolversType), POINTER :: solvers

    ENTERS("Problem_CellMLEquationsGet1",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLEquations)) CALL FlagError("CellML equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
#endif    
 
    NULLIFY(controlLoopRoot)
    CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
    NULLIFY(controlLoop)
    CALL ControlLoop_Get(controlLoopRoot,controlLoopIdentifiers,controlLoop,err,error,*999)
    NULLIFY(solvers)
    CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
    NULLIFY(solver)
    CALL Solvers_SolverGet(solvers,solverIndex,solver,err,error,*999)
    CALL Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*999)
    
    EXITS("Problem_CellMLEquationsGet1")
    RETURN
999 NULLIFY(cellMLEquations)
998 ERRORSEXITS("Problem_CellMLEquationsGet1",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_CellMLEquationsGet1
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the context for a problem.
  SUBROUTINE Problem_ContextGet(problem,context,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to get the context for
    TYPE(ContextType), POINTER :: context !<On exit, a pointer to the context for the problem. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Problem_ContextGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(context)) CALL FlagError("Context is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(problem%problems)) THEN
      localError="Problems is not associated for problem number "// &
        & TRIM(NumberToVString(problem%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    context=>problem%problems%context

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(context)) THEN
      localError="The context is not associated for the problems for problem number "// &
        & TRIM(NumberToVString(problem%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("Problem_ContextGet")
    RETURN
999 NULLIFY(context)
998 ERRORSEXITS("Problem_ContextGet",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_ContextGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the control loop for a problem. \see OpenCMISS::Iron::cmfe_Problem_ControlLoopGet
  SUBROUTINE Problem_ControlLoopGet0(problem,controlLoopIdentifier,controlLoop,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to get the control loop for.
    INTEGER(INTG), INTENT(IN) :: controlLoopIdentifier !<The control loop identifier
    TYPE(ControlLoopType), POINTER :: controlLoop !<On return, a pointer to the control loop. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Problem_ControlLoopGet0",err,error,*999)

    CALL Problem_ControlLoopGet1(problem,[controlLoopIdentifier],controlLoop,err,error,*999) 
       
    EXITS("Problem_ControlLoopGet0")
    RETURN
999 ERRORSEXITS("Problem_ControlLoopGet0",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_ControlLoopGet0
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the control_loop for a problem. \see OpenCMISS::Iron::cmfe_Problem_ControlLoopGet
  SUBROUTINE Problem_ControlLoopGet1(problem,controlLoopIdentifiers,controlLoop,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to get the control loop for.
    INTEGER(INTG), INTENT(IN) :: controlLoopIdentifiers(:) !<The control loop identifier.
    TYPE(ControlLoopType), POINTER :: controlLoop !<On return, a pointer to the control loop. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Problem_ControlLoopGet1",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(controlLoop)) CALL FlagError("Control loop is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(problem%controlLoop)) CALL FlagError("Problem control loop is not associated.",err,error,*999)
#endif    

    CALL ControlLoop_Get(problem%controlLoop,controlLoopIdentifiers,controlLoop,err,error,*999)
    
    EXITS("Problem_ControlLoopGet1")
    RETURN
999 NULLIFY(controlLoop)
998 ERRORSEXITS("Problem_ControlLoopGet1",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_ControlLoopGet1
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the root control loop for a problem. 
  SUBROUTINE Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<a pointer to the problem to get the control loop root for.
    TYPE(ControlLoopType), POINTER :: controlLoopRoot !<On return, a pointer to the control loop root for the problem. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("Problem_ControlLoopRootGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(controlLoopRoot)) CALL FlagError("Control loop root is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
#endif    

    controlLoopRoot=>problem%controlLoop

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(controlLoopRoot)) THEN
      localError="The problem control loop root is not associated for problem number "// &
        & TRIM(NumberToVString(problem%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
    EXITS("Problem_ControlLoopRootGet")
    RETURN
999 NULLIFY(controlLoopRoot)
998 ERRORSEXITS("Problem_ControlLoopRootGet",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_ControlLoopRootGet

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the problem with the given user number. 
  SUBROUTINE Problem_Get(problems,userNumber,problem,err,error,*)

    !Argument variables
    TYPE(ProblemsType), POINTER :: problems !<The problems to get the user number for.
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the problem to find
    TYPE(ProblemType), POINTER :: problem !<On exit, a pointer to the problem with the specified user number if it exists. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("Problem_Get",err,error,*999)

    CALL Problem_UserNumberFind(problems,userNumber,problem,err,error,*999)
#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(problem)) THEN
      localError="A problem with an user number of "//TRIM(NumberToVString(userNumber,"*",err,error))//" does not exist."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
  
    EXITS("Problem_Get")
    RETURN
999 ERRORSEXITS("Problem_Get",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_Get

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the problems for a problem. 
  SUBROUTINE Problem_ProblemsGet(problem,problems,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<a pointer to the problem to get the problems for.
    TYPE(ProblemsType), POINTER :: problems !<On return, a pointer to the problems for the problem. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("Problem_ProblemsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(problems)) CALL FlagError("Problems is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
#endif    

    problems=>problem%problems

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(problems)) THEN
      localError="The problem problems is not associated for problem number "// &
        & TRIM(NumberToVString(problem%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
    EXITS("Problem_ProblemsGet")
    RETURN
999 NULLIFY(problems)
998 ERRORSEXITS("Problem_ProblemsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_ProblemsGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solver for a problem control loop. \see OpenCMISS::Iron::cmfe_Problem_SolverGet
  SUBROUTINE Problem_SolverGet0(problem,controlLoopIdentifier,solverIndex,solver,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to get the solver for.
    INTEGER(INTG), INTENT(IN) :: controlLoopIdentifier !<The control loop identifier
    INTEGER(INTG), INTENT(IN) :: solverIndex !<The solver index to get the solver for.
    TYPE(SolverType), POINTER :: solver !<On return, a pointer to the solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Problem_SolverGet0",err,error,*999)

    CALL Problem_SolverGet1(problem,[controlLoopIdentifier],solverIndex,solver,err,error,*999) 
       
    EXITS("Problem_SolverGet0")
    RETURN
999 ERRORSEXITS("Problem_SolverGet0",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverGet0
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solver for a problem control loop. \see OpenCMISS::Iron::cmfe_Problem_SolverGet
  SUBROUTINE Problem_SolverGet1(problem,controlLoopIdentifiers,solverIndex,solver,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to get the solver for.
    INTEGER(INTG), INTENT(IN) :: controlLoopIdentifiers(:) !<The control loop identifier to get the solver for.
    INTEGER(INTG), INTENT(IN) :: solverIndex !<The solver index to get the solver for.
    TYPE(SolverType), POINTER :: solver !<On return, a pointer to the solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: controlLoop,controlLoopRoot
    TYPE(SolversType), POINTER :: solvers
 
    ENTERS("Problem_SolverGet1",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solver)) CALL FlagError("Solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
#endif    


    NULLIFY(controlLoopRoot)
    CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
    NULLIFY(controlLoop)
    CALL ControlLoop_Get(controlLoopRoot,controlLoopIdentifiers,controlLoop,err,error,*999)
    NULLIFY(solvers)
    CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
    CALL Solvers_SolverGet(solvers,solverIndex,solver,err,error,*999)
     
    EXITS("Problem_SolverGet1")
    RETURN
999 NULLIFY(solver)
998 ERRORSEXITS("Problem_SolverGet1",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverGet1
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to a solver equations defined with a solver. \see OpenCMISS::Iron::cmfe_Problem_SolverEquationsGet
  SUBROUTINE Problem_SolverEquationsGet0(problem,controlLoopIdentifier,solverIndex,solverEquations,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to get solver equations for
    INTEGER(INTG), INTENT(IN) :: controlLoopIdentifier !<The control loop identifier to get the solver equations for
    INTEGER(INTG), INTENT(IN) :: solverIndex !<The solver index in the solvers to get the solver equations for
    TYPE(SolverEquationsType), POINTER :: solverEquations !<On exit, a pointer to the specified solver equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Problem_SolverEquationsGet0",err,error,*999)

    CALL Problem_SolverEquationsGet1(problem,[controlLoopIdentifier],solverIndex,solverEquations,err,error,*999)
    
    EXITS("Problem_SolverEquationsGet0")
    RETURN
999 ERRORSEXITS("Problem_SolverEquationsGet0",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverEquationsGet0

  !
  !================================================================================================================================
  !

  !>Returns a pointer to a solver equations defined with a solver. \see OpenCMISS::Iron::cmfe_Problem_SolverEquationsGet
  SUBROUTINE Problem_SolverEquationsGet1(problem,controlLoopIdentifiers,solverIndex,solverEquations,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to get solver equations for
    INTEGER(INTG), INTENT(IN) :: controlLoopIdentifiers(:) !<controlLoopIdentifiers(identifierIdx). The control loop identifiers to get the solver equations for
    INTEGER(INTG), INTENT(IN) :: solverIndex !<The solver index in the solvers to get the solver equations for
    TYPE(SolverEquationsType), POINTER :: solverEquations !<On exit, a pointer to the specified solver equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: controlLoop,controlLoopRoot
    TYPE(SolverType), POINTER :: solver
    TYPE(SolversType), POINTER :: solvers

    ENTERS("Problem_SolverEquationsGet1",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
#endif    

    NULLIFY(controlLoopRoot)
    CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
    NULLIFY(controlLoop)
    CALL ControlLoop_Get(controlLoopRoot,controlLoopIdentifiers,controlLoop,err,error,*999)
    NULLIFY(solvers)
    CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
    NULLIFY(solver)
    CALL Solvers_SolverGet(solvers,solverIndex,solver,err,error,*999)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        
    EXITS("Problem_SolverEquationsGet1")
    RETURN
999 NULLIFY(solverEquations)
998 ERRORSEXITS("Problem_SolverEquationsGet1",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverEquationsGet1
  
  !
  !================================================================================================================================
  !

  !>Returns the problem specification i.e., problem class, type and subtype for an problem. \see OpenCMISS::Iron::cmfe_Problem_SpecificationGet
  SUBROUTINE Problem_SpecificationGet(problem,minSpecificationLength,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to get the specification for
    INTEGER(INTG), INTENT(IN) :: minSpecificationLength !<The minimum number of problem specification identifiers to return. If the minSpecificationLength is zero then all specification identifiers are required. 
    INTEGER(INTG), INTENT(INOUT) :: problemSpecification(:) !<On return, The problem specifcation array. Must be allocated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: specificationLength
    TYPE(VARYING_STRING) :: localError

    ENTERS("Problem_SpecificationGet",err,error,*999)
   
#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) THEN
      localError="The specification array has not been allocated for problem number "// &
        & TRIM(NumberToVString(problem%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    IF(minSpecificationLength==0) THEN
      specificationLength=problem%specificationLength
    ELSE
      IF(minSpecificationLength>problem%specificationLength) THEN
        localError="The requested minimum number of specification parameters of "// &
          & TRIM(NumberToVString(minSpecificationLength,"*",err,error))//" for problem number "// &
          & TRIM(NumberToVString(problem%userNumber,"*",err,error))//" is too large. The problem only has "// &
          & TRIM(NumberToVString(problem%specificationLength,"*",err,error))//" specification parameters."
        CALL FlagError(localError,err,error,*999)
      ELSE
        specificationLength=minSpecificationLength
      ENDIF
    ENDIF
#ifdef WITH_PRECHECKS
    IF(SIZE(problemSpecification,1)<specificationLength) THEN
      localError="The problem specification array size is "// &
        & TRIM(NumberToVstring(SIZE(problemSpecification,1),"*",err,error))// &
        & " and it needs to be >= "//TRIM(NumberToVstring(specificationLength,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    problemSpecification(1:specificationLength)=problem%specification(1:specificationLength)

    EXITS("Problem_SpecificationGet")
    RETURN
999 ERRORSEXITS("Problem_SpecificationGet",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SpecificationGet

  !
  !================================================================================================================================
  !

  !>Gets the size of the problem specification array for a problem identified by a pointer. \see OpenCMISS::Iron::cmfe_Problem_SpecificationSizeGet
  SUBROUTINE Problem_SpecificationSizeGet(problem,specificationSize,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to get the specification for.
    INTEGER(INTG), INTENT(OUT) :: specificationSize !<On return, the size of the problem specifcation array.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Problem_SpecificationSizeGet",err,error,*999)

    specificationSize=0
#ifdef WITH_PRECHECKS
    CALL Problem_AssertIsFinished(problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
#endif    
    
    specificationSize=SIZE(problem%specification,1)

    EXITS("Problem_SpecificationSizeGet")
    RETURN
999 ERRORSEXITS("Problem_SpecificationSizeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SpecificationSizeGet

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the problem identified by a user number. If no problem with that user number exists problem is left nullified.
  SUBROUTINE Problem_UserNumberFind(problems,userNumber,problem,err,error,*)

    !Argument variables
    TYPE(ProblemsType), POINTER :: problems !<The problems to find the user number for.
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to find.
    TYPE(ProblemType), POINTER :: problem !<On return, a pointer to the problem with the given user number. If no problem with that user number exists then the pointer is returned as NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemIdx
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("Problem_UserNumberFind",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(problem)) CALL FlagError("Problem is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(problems)) CALL FlagError("Problems is not associated.",err,error,*999)
#endif    
   
    IF(ASSOCIATED(problems%problems)) THEN
      DO problemIdx=1,problems%numberOfProblems
#ifdef WITH_PRECHECKS        
        IF(.NOT.ASSOCIATED(problems%problems(problemIdx)%ptr)) THEN
          localError="The problem pointer in problems is not associated for problem index "// &
            & TRIM(NumberToVString(problemIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
#endif        
        IF(problems%problems(problemIdx)%ptr%userNumber==userNumber) THEN
          problem=>problems%problems(problemIdx)%ptr
          EXIT
        ENDIF
      ENDDO !problemIdx
    ENDIF
    
    EXITS("Problem_UserNumberFind")
    RETURN
999 NULLIFY(problem)
998 ERRORSEXITS("Problem_UserNumberFind",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_UserNumberFind

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the context for a problems.
  SUBROUTINE Problems_ContextGet(problems,context,err,error,*)

    !Argument variables
    TYPE(ProblemsType), POINTER :: problems !<A pointer to the problems to get the context for
    TYPE(ContextType), POINTER :: context !<On exit, a pointer to the context for the problems. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
  
    ENTERS("Problems_ContextGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(context)) CALL FlagError("Context is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(problems)) CALL FlagError("Problems is not associated.",err,error,*999)
#endif    

    context=>problems%context

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(context)) CALL FlagError("The context is not associated for the problems.",err,error,*999)
#endif    
    
    EXITS("Problems_ContextGet")
    RETURN
999 NULLIFY(context)
998 ERRORSEXITS("Problems_ContextGet",err,error)
    RETURN 1
    
  END SUBROUTINE Problems_ContextGet

  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the specified problem in the list of problems.
  SUBROUTINE Problems_ProblemGet(problems,problemIndex,problem,err,error,*)

    !Argument variables
    TYPE(ProblemsType), POINTER :: problems !<A pointer to the problems to get the problem for
    INTEGER(INTG), INTENT(IN) :: problemIndex !<The specified problem to get
    TYPE(ProblemType), POINTER :: problem !<On exit, a pointer to the specified problem. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Problems_ProblemGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(problem)) CALL FlagError("Problem is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(problems)) CALL FlagError("Problems is not associated.",err,error,*999)
    IF(problemIndex<=0.OR.problemIndex>problems%numberOfProblems) THEN
      localError="The specified problem index of "//TRIM(NumberToVString(problemIndex,"*",err,error))// &
        & " is invalid. The problem index must be >= 1 and <= "// &
        & TRIM(NumberToVString(problems%numberOfProblems,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ASSOCIATED(problems%problems)) CALL FlagError("Problems problems is not associated.",err,error,*999)
#endif    
      
    problem=>problems%problems(problemIndex)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(problem)) THEN
      localError="The problems problem is not associated for problem index "// &
        & TRIM(NumberToVString(problemIndex,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
    EXITS("Problems_ProblemGet")
    RETURN
999 NULLIFY(problem)
998 ERRORSEXITS("Problems_ProblemGet",err,error)
    RETURN 1
    
  END SUBROUTINE Problems_ProblemGet

  !
  !================================================================================================================================
  !
 
END MODULE ProblemAccessRoutines
