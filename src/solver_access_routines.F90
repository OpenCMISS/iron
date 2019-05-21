!> \file
!> \author Chris Bradley
!> \brief This module contains all solver access method routines.
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

!> This module contains all solver access method routines.
MODULE SolverAccessRoutines
  
  USE BaseRoutines
  USE Constants
  USE ISO_VARYING_STRING
  USE Kinds
  USE Strings
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup SOLVER_ROUTINES_SolverTypes SOLVER_ROUTINES::SolverTypes
  !> \brief The types of solver
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_NUMBER_OF_SOLVER_TYPES=9 !<Number of different solver types possible \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_LINEAR_TYPE=1 !<A linear solver \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_NONLINEAR_TYPE=2 !<A nonliXnear solver  \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_TYPE=3 !<A dynamic solver \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_TYPE=4 !<A differential-algebraic equation solver \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_EIGENPROBLEM_TYPE=5 !<A eigenproblem solver \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_TYPE=6 !<An optimiser solver \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_CELLML_EVALUATOR_TYPE=7 !<A CellML evaluation solver \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_STATE_ITERATION_TYPE=8 !<An state iteration solver \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_GEOMETRIC_TRANSFORMATION_TYPE=9 !<An geometric transformation solver \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_SolverLibraries SOLVER_ROUTINES::SolverLibraries
  !> \brief The types of solver libraries
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_CMISS_LIBRARY=LIBRARY_CMISS_TYPE !<CMISS (internal) solver library \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_PETSC_LIBRARY=LIBRARY_PETSC_TYPE !<PETSc solver library \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_MUMPS_LIBRARY=LIBRARY_MUMPS_TYPE !<MUMPS solver library \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_SUPERLU_LIBRARY=LIBRARY_SUPERLU_TYPE !<SuperLU solver library \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_SPOOLES_LIBRARY=LIBRARY_SPOOLES_TYPE !<Spooles solver library \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_UMFPACK_LIBRARY=LIBRARY_UMFPACK_TYPE !<UMFPACK solver library \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_LUSOL_LIBRARY=LIBRARY_LUSOL_TYPE !<LUSOL solver library \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ESSL_LIBRARY=LIBRARY_ESSL_TYPE !<ESSL solver library \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_LAPACK_LIBRARY=LIBRARY_LAPACK_TYPE !<LAPACK solver library \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_HYPRE_LIBRARY=LIBRARY_HYPRE_TYPE !<Hypre solver library \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_PASTIX_LIBRARY=LIBRARY_PASTIX_TYPE !<PaStiX solver library \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
   !>@}

  !> \addtogroup SOLVER_ROUTINES_LinearSolverTypes SOLVER_ROUTINES::LinearSolverTypes
  !> \brief The types of linear solvers
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_LINEAR_DIRECT_SOLVE_TYPE=1 !<Direct linear solver type \see SOLVER_ROUTINES_LinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE=2 !<Iterative linear solver type \see SOLVER_ROUTINES_LinearSolverTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_DirectLinearSolverTypes SOLVER_ROUTINES::DirectLinearSolverTypes
  !> \brief The types of direct linear solvers
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_DIRECT_LU=1 !<LU direct linear solver \see SOLVER_ROUTINES_DirectLinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DIRECT_CHOLESKY=2 !<Cholesky direct linear solver \see SOLVER_ROUTINES_DirectLinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DIRECT_SVD=3 !<SVD direct linear solver \see SOLVER_ROUTINES_DirectLinearSolverTypes,SOLVER_ROUTINES
  !>@}
  
  !> \addtogroup SOLVER_ROUTINES_IterativeLinearSolverTypes SOLVER_ROUTINES::IterativeLinearSolverTypes
  !> \brief The types of iterative linear solvers
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_RICHARDSON=1 !<Richardson iterative solver type \see SOLVER_ROUTINES_IterativeLinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_CHEBYSHEV=2 !<Chebyshev iterative solver type \see SOLVER_ROUTINES_IterativeLinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_CONJUGATE_GRADIENT=3 !<Conjugate gradient iterative solver type \see SOLVER_ROUTINES_IterativeLinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_BICONJUGATE_GRADIENT=4 !<Bi-conjugate gradient iterative solver type \see SOLVER_ROUTINES_IterativeLinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_GMRES=5 !<Generalised minimum residual iterative solver type \see SOLVER_ROUTINES_IterativeLinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_BiCGSTAB=6 !<Stabalised bi-conjugate gradient iterative solver type \see SOLVER_ROUTINES_IterativeLinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_CONJGRAD_SQUARED=7 !<Conjugate gradient squared iterative solver type \see SOLVER_ROUTINES_IterativeLinearSolverTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_IterativePreconditionerTypes SOLVER_ROUTINES::IterativePreconditionerTypes
  !> \brief The types of iterative preconditioners
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_NO_PRECONDITIONER=0 !<No preconditioner type \see SOLVER_ROUTINES_IterativePreconditionerTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_JACOBI_PRECONDITIONER=1 !<Jacobi preconditioner type \see SOLVER_ROUTINES_IterativePreconditionerTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_BLOCK_JACOBI_PRECONDITIONER=2 !<Iterative block Jacobi preconditioner type \see SOLVER_ROUTINES_IterativePreconditionerTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_SOR_PRECONDITIONER=3 !<Successive over relaxation preconditioner type \see SOLVER_ROUTINES_IterativePreconditionerTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_INCOMPLETE_CHOLESKY_PRECONDITIONER=4 !<Incomplete Cholesky preconditioner type \see SOLVER_ROUTINES_IterativePreconditionerTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER=5 !<Incomplete LU preconditioner type \see SOLVER_ROUTINES_IterativePreconditionerTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_ADDITIVE_SCHWARZ_PRECONDITIONER=6 !<Additive Schwrz preconditioner type \see SOLVER_ROUTINES_IterativePreconditionerTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_NonlinearSolverTypes SOLVER_ROUTINES::NonlinearSolverTypes
  !> \brief The types of nonlinear solvers
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_NONLINEAR_NEWTON=1 !<Newton nonlinear solver type \see SOLVER_ROUTINES_NonlinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_NONLINEAR_BFGS_INVERSE=2 !<BFGS inverse nonlinear solver type \see SOLVER_ROUTINES_NonlinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_NONLINEAR_SQP=3 !<Sequential Quadratic Program nonlinear solver type \see SOLVER_ROUTINES_NonlinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_NONLINEAR_QUASI_NEWTON=4 !<Sequential Quasi-Newton nonlinear solver type \see SOLVER_ROUTINES_NonlinearSolverTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_QuasiNewtonSolverTypes SOLVER_ROUTINES::QuasiNewtonSolverTypes
  !> \brief The types of nonlinear Quasi-Newton solvers
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_LINESEARCH=1 !<Quasi-Newton line search nonlinear solver type \see SOLVER_ROUTINES_QuasiNewtonSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_TRUSTREGION=2 !<Quasi-Newton trust region nonlinear solver type \see SOLVER_ROUTINES_QuasiNewtonSolverTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_QuasiNewtonTypes SOLVER_ROUTINES::QuasiNewtonTypes
  !> \brief The nonlinear Quasi-Newton types
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_LBFGS=1 !<LBFGS Quasi-Newton type \see SOLVER_ROUTINES_QuasiNewtonTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_GOODBROYDEN=2 !<"Good" Broyden Quasi-Newton type \see SOLVER_ROUTINES_QuasiNewtonTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_BADBROYDEN=3 !<"Bad" Broyden Quasi-Newton type \see SOLVER_ROUTINES_QuasiNewtonTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_QuasiNewtonLineSearchTypes SOLVER_ROUTINES::QuasiNewtonLineSearchTypes
  !> \brief The types line search techniques for Quasi-Newton line search nonlinear solvers
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_LINESEARCH_BASIC=1 !<Simple damping line search. \see SOLVER_ROUTINES_QuasiNewtonLineSearchTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_LINESEARCH_L2=2 !<Secant line search over the L2 norm of the function  \see SOLVER_ROUTINES_QuasiNewtonLineSearchTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_LINESEARCH_CP=3 !<Critical point secant line search \see SOLVER_ROUTINES_QuasiNewtonLineSearchTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_QuasiNewtonRestartTypes SOLVER_ROUTINES::QuasiNewtonRestartTypes
  !> \brief The nonlinear Quasi-Newton restart types
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_RESTART_NONE=1 !<Never restart \see SOLVER_ROUTINES_QuasiNewtonRestartTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_RESTART_POWELL=2 !<Restart based upon descent criteria \see SOLVER_ROUTINES_QuasiNewtonRestartTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_RESTART_PERIODIC=3 !<Restart after a fixed number of iterations \see SOLVER_ROUTINES_QuasiNewtonRestartTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_QuasiNewtonScaleTypes SOLVER_ROUTINES::QuasiNewtonScaleTypes
  !> \brief The nonlinear Quasi-Newton scale types
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_SCALE_NONE=1 !<Don't scale the problem \see SOLVER_ROUTINES_QuasiNewtonScaleTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_SCALE_SHANNO=2 !<Use Shanno scaling \see SOLVER_ROUTINES_QuasiNewtonScaleTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_SCALE_LINESEARCH=3 !<Scale based upon line search lambda \see SOLVER_ROUTINES_QuasiNewtonScaleTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_SCALE_JACOBIAN=4 !<Scale by inverting a previously computed Jacobian \see SOLVER_ROUTINES_QuasiNewtonScaleTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_NewtonSolverTypes SOLVER_ROUTINES::NewtonSolverTypes
  !> \brief The types of nonlinear Newton solvers
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_LINESEARCH=1 !<Newton line search nonlinear solver type \see SOLVER_ROUTINES_NewtonSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_TRUSTREGION=2 !<Newton trust region nonlinear solver type \see SOLVER_ROUTINES_NewtonSolverTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_NewtonLineSearchTypes SOLVER_ROUTINES::NewtonLineSearchTypes
  !> \brief The types line search techniques for Newton line search nonlinear solvers
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_LINESEARCH_NONORMS=1 !<No norms line search for Newton line search nonlinear solves \see SOLVER_ROUTINES_NewtonLineSearchTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_LINESEARCH_LINEAR=2 !<Linear search for Newton line search nonlinear solves \see SOLVER_ROUTINES_NewtonLineSearchTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_LINESEARCH_QUADRATIC=3 !<Quadratic search for Newton line search nonlinear solves \see SOLVER_ROUTINES_NewtonLineSearchTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_LINESEARCH_CUBIC=4!<Cubic search for Newton line search nonlinear solves \see SOLVER_ROUTINES_NewtonLineSearchTypes,SOLVER_ROUTINES
  !>@}
  
  !> \addtogroup SOLVER_ROUTINES_JacobianCalculationTypes SOLVER_ROUTINES::JacobianCalculationTypes
  !> \brief The Jacobian calculation types for a nonlinear solver 
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_JACOBIAN_NOT_CALCULATED=1 !<The Jacobian values will not be calculated for the nonlinear equations set \see SOLVER_ROUTINES_JacobianCalculationTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED=2 !<The Jacobian values will be calculated analytically for the nonlinear equations set \see SOLVER_ROUTINES_JacobianCalculationTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_JACOBIAN_FD_CALCULATED=3 !<The Jacobian values will be calculated using finite differences for the nonlinear equations set \see SOLVER_ROUTINES_JacobianCalculationTypes,SOLVER_ROUTINES
  !>@} 
  
  !> \addtogroup SOLVER_ROUTINES_NewtonConvergenceTestTypes SOLVER_ROUTINES::NewtonConvergenceTestTypes
  !> \brief The convergence test types for a nonlinear solver 
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_CONVERGENCE_PETSC_DEFAULT=1 !<Petsc default convergence test \see SOLVER_ROUTINES_NewtonConvergenceTestTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_CONVERGENCE_ENERGY_NORM=2 !<Energy norm convergence test \see SOLVER_ROUTINES_NewtonConvergenceTestTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_CONVERGENCE_DIFFERENTIATED_RATIO=3 !<Sum of differentiated ratios of unconstrained to constrained residuals convergence test \see SOLVER_ROUTINES_NewtonConvergenceTestTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_DynamicOrderTypes SOLVER_ROUTINES::DynamicOrderTypes
  !> \brief The order types for a dynamic solver 
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_FIRST_ORDER=1 !<Dynamic solver has first order terms \see SOLVER_ROUTINES_DynamicOrderTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_SECOND_ORDER=2 !<Dynamic solver has second order terms \see SOLVER_ROUTINES_DynamicOrderTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_DynamicLinearityTypes SOLVER_ROUTINES::DynamicLinearityTypes
  !> \brief The time linearity types for a dynamic solver 
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_LINEAR=1 !<Dynamic solver has linear terms \see SOLVER_ROUTINES_DynamicLinearityTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_NONLINEAR=2 !<Dynamic solver has nonlinear terms \see SOLVER_ROUTINES_DynamicLinearityTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_DynamicDegreeTypes SOLVER_ROUTINES::DynamicDegreeTypes
  !> \brief The time interpolation polynomial degree types for a dynamic solver 
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_FIRST_DEGREE=1 !<Dynamic solver uses a first degree polynomial for time interpolation \see SOLVER_ROUTINES_DynamicDegreeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_SECOND_DEGREE=2 !<Dynamic solver uses a second degree polynomial for time interpolation \see SOLVER_ROUTINES_DynamicDegreeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_THIRD_DEGREE=3 !<Dynamic solver uses a third degree polynomial for time interpolation \see SOLVER_ROUTINES_DynamicDegreeTypes,SOLVER_ROUTINES
  !>@}    
  
  !> \addtogroup SOLVER_ROUTINES_DynamicSchemeTypes SOLVER_ROUTINES::DynamicSchemeTypes
  !> \brief The types of dynamic solver scheme
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_EULER_SCHEME=1 !<Euler (explicit) dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_BACKWARD_EULER_SCHEME=2 !<Backward Euler (implicit) dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME=3 !<Crank-Nicolson dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_GALERKIN_SCHEME=4 !<Galerkin dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_ZLAMAL_SCHEME=5 !<Zlamal dynamic solver \see SOLVER_ROUTINES_DynamicorderTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_SECOND_DEGREE_GEAR_SCHEME=6 !<2nd degree Gear dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_SECOND_DEGREE_LINIGER1_SCHEME=7 !<1st 2nd degree Liniger dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_SECOND_DEGREE_LINIGER2_SCHEME=8 !<2nd 2nd degree Liniger dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_NEWMARK1_SCHEME=9 !<1st Newmark dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_NEWMARK2_SCHEME=10 !<2nd Newmark dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_NEWMARK3_SCHEME=11 !<3rd Newmark dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_THIRD_DEGREE_GEAR_SCHEME=12 !<3rd degree Gear dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_THIRD_DEGREE_LINIGER1_SCHEME=13 !<1st 3rd degree Liniger dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_THIRD_DEGREE_LINIGER2_SCHEME=14 !<2nd 3rd degree Liniger dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_HOUBOLT_SCHEME=15 !<Houbolt dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_WILSON_SCHEME=16 !<Wilson dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_BOSSAK_NEWMARK1_SCHEME=17 !<1st Bossak-Newmark dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_BOSSAK_NEWMARK2_SCHEME=18 !<2nd Bossak-Newmark dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_HILBERT_HUGHES_TAYLOR1_SCHEME=19 !<1st Hilbert-Hughes-Taylor dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_HILBERT_HUGHES_TAYLOR2_SCHEME=20 !<1st Hilbert-Hughes-Taylor dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_USER_DEFINED_SCHEME=21 !<User specified degree and theta dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  !>@}
  
  !> \addtogroup SOLVER_ROUTINES_DynamicStartupTypes SOLVER_ROUTINES::DynamicStartupTypes
  !> \brief The dynamic solver previous values startup type
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_PREVIOUS_STARTUP_TYPE=1 !<Dynamic solver previous values startup type \see SOLVER_ROUTINES_DynamicStartupTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_PREVIOUS_RESTART_TYPE=2 !<Dynamic solver previous values startup type type \see SOLVER_ROUTINES_DynamicStartupTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_DAETypes SOLVER_ROUTINES::DAETypes
  !> \brief The type of differential-algebraic equation
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_DIFFERENTIAL_ONLY=0 !<Differential equations only \see SOLVER_ROUTINES_DAETypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_INDEX_1=1 !<Index 1 differential-algebraic equation \see SOLVER_ROUTINES_DAETypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_INDEX_2=2 !<Index 2 differential-algebraic equation \see SOLVER_ROUTINES_DAETypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_INDEX_3=3 !<Index 3 differential-algebraic equation \see SOLVER_ROUTINES_DAETypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_DAESolverTypes SOLVER_ROUTINES::DAESolverTypes
  !> \brief The differential-algebraic equation solver types for a differential-algebraic equation solver 
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_EULER=1 !<Euler differential-algebraic equation solver \see SOLVER_ROUTINES_DAESolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_CRANK_NICOLSON=2 !<Crank-Nicolson differential-algebraic equation solver \see SOLVER_ROUTINES_DAESolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_RUNGE_KUTTA=3 !<Runge-Kutta differential-algebraic equation solver \see SOLVER_ROUTINES_DAESolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_ADAMS_MOULTON=4 !<Adams-Moulton differential-algebraic equation solver \see SOLVER_ROUTINES_DAESolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_BDF=5 !<General BDF differential-algebraic equation solver \see SOLVER_ROUTINES_DAESolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_RUSH_LARSON=6 !<Rush-Larson differential-algebraic equation solver \see SOLVER_ROUTINES_DAESolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_EXTERNAL=7 !<External (e.g., CellML generated) differential-algebraic equation solver \see SOLVER_ROUTINES_DAESolverTypes,SOLVER_ROUTINES
  
  !>@}

  !> \addtogroup SOLVER_ROUTINES_EulerDAESolverTypes SOLVER_ROUTINES::EulerDAESolverTypes
  !> \brief The Euler solver types for a differential-algebriac equation solver 
  !> \see SOLVER_ROUTINES_DAESolverTypes,SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_EULER_FORWARD=1 !<Forward Euler differential equation solver \see SOLVER_ROUTINES_EulerDAESolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_EULER_BACKWARD=2 !<Backward Euler differential equation solver \see SOLVER_ROUTINES_EulerDAESolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_EULER_IMPROVED=3 !<Improved Euler differential equation solver \see SOLVER_ROUTINES_EulerDAESolverTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_OptimiserObjectiveTypes SOLVER_ROUTINES::OptimiserObjectiveTypes
  !> \brief The types of objective for an optimisation problem
  !> \see SOLVER_ROUTINES_OptimiserObjectiveTypes,SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_NO_OBJECTIVE=1 !<Optimisation problem has no objective (feasibility problem) \see SOLVER_ROUTINES_OptimiserObjectiveTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_ONE_OBJECTIVE=2 !<Optimisation problem has one objective (standard problem) \see SOLVER_ROUTINES_OptimiserObjectiveTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_MANY_OBJECTIVE=3 !<Optimisation problem has many objectives (multi-objective problem) \see SOLVER_ROUTINES_OptimiserObjectiveTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_OptimiserVariableTypes SOLVER_ROUTINES::OptimiserVariableTypes
  !> \brief The types of variables for an optimisation problem
  !> \see SOLVER_ROUTINES_OptimiserVariablesTypes,SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_CONTINUOUS_VARIABLES=1 !<Optimisation problem has continuous variables \see SOLVER_ROUTINES_OptimiserVariableTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_DISCRETE_VARIABLES=2 !<Optimisation problem has discrete variables \see SOLVER_ROUTINES_OptimiserVariableTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_OptimiserCertaintyTypes SOLVER_ROUTINES::OptimiserCertaintyTypes
  !> \brief The types of certainty for an optimisation problem
  !> \see SOLVER_ROUTINES_OptimiserCertaintyTypes,SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_DETERMINISTIC_CERTAINTY=1 !<Optimisation problem is determanistic \see SOLVER_ROUTINES_OptimiserCertaintyTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_STOCHASTIC_CERTAINTY=2 !<Optimisation problem is stochastic \see SOLVER_ROUTINES_OptimiserCertaintyTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_OptimiserConstraintTypes SOLVER_ROUTINES::OptimiserConstraintTypes
  !> \brief The types of constraints for an optimisation problem
  !> \see SOLVER_ROUTINES_OptimiserConstraintTypes,SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_UNCONSTRAINED=1 !<Unconstrained optimisation problem \see SOLVER_ROUTINES_OptimiserConstraintTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_BOUND_CONSTRAINED=2 !<Optimisation problem with bounds on the variables \see SOLVER_ROUTINES_OptimiserConstraintTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_LINEAR_CONSTRAINTS=3 !<Optimisation with linear constraints \see SOLVER_ROUTINES_OptimiserConstraintTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_NONLINEAR_CONSTRAINTS=4 !<Optimisation with non-linear constraints \see SOLVER_ROUTINES_OptimiserConstraintTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_PDE_CONSTRAINTS=5 !<PDE constrained optimisation \see SOLVER_ROUTINES_OptimiserConstraintTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_OptimiserSolverTypes SOLVER_ROUTINES::OptimiserSolverTypes
  !> \brief The types of solver for an optimisation problem
  !> \see SOLVER_ROUTINES_OptimiserSolverTypes,SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_NONLINEAR_LEAST_SQUARES=1 !<Non-linear least squares optimisation problem \see SOLVER_ROUTINES_OptimiserSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_LINEAR_PROGRAMMING=1 !<Linear programming optimisation problem \see SOLVER_ROUTINES_OptimiserSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_QUADRATIC_PROGRAMMING=1 !<Quadratic programming optimisation problem \see SOLVER_ROUTINES_OptimiserSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_INTEGER_PROGRAMMING=1 !<Integer programming optimisation problem \see SOLVER_ROUTINES_OptimiserSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_LINEAR_COMPLEMENTARITY=1 !<Linear complementarity optimisation problem \see SOLVER_ROUTINES_OptimiserSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_NONLINEAR_COMPLEMENTARITY=1 !<Nonlinear complementarity optimisation problem \see SOLVER_ROUTINES_OptimiserSolverTypes,SOLVER_ROUTINES
  !>@}
  
  !> \addtogroup SOLVER_ROUTINES_OptimiserGradientCalculationTypes SOLVER_ROUTINES::OptimiserGradientCalculationTypes
  !> \brief The gradient calculation types for an optimiser solver 
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_GRADIENT_NOT_CALCULATED=1 !<The gradient values will not be calculated for the optimiser equations set \see SOLVER_ROUTINES_OptimiserGradientCalculationTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_GRADIENT_EQUATIONS_CALCULATED=2 !<The gradient values will be calculated analytically for the optimiser equations set \see SOLVER_ROUTINES_OptimiserGradientCalculationTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_GRADIENT_FD_CALCULATED=3 !<The gradient values will be calculated using finite differences for the optimiser equations set \see SOLVER_ROUTINES_OptimiserGradientCalculationTypes,SOLVER_ROUTINES
  !>@} 
  
  !> \addtogroup SOLVER_ROUTINES_OptimiserHessianCalculationTypes SOLVER_ROUTINES::OptimiserHessianCalculationTypes
  !> \brief The Hessian calculation types for an optimiser solver 
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_HESSIAN_NOT_CALCULATED=1 !<The Hessian values will not be calculated for the optimiser equations set \see SOLVER_ROUTINES_OptimiserHessianCalculationTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_HESSIAN_EQUATIONS_CALCULATED=2 !<The Hessian values will be calculated analytically for the optimiser equations set \see SOLVER_ROUTINES_OptimiserHessianCalculationTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_HESSIAN_FD_CALCULATED=3 !<The Hessian values will be calculated using finite differences for the optimiser equations set \see SOLVER_ROUTINES_OptimiserHessianCalculationTypes,SOLVER_ROUTINES
  !>@} 
  
  !> \addtogroup SOLVER_ROUTINES_SolutionInitialiseTypes SOLVER_ROUTINES::SolutionInitialiseTypes
  !> \brief The types of solution initialisation
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_SOLUTION_INITIALISE_ZERO=0 !<Initialise the solution by zeroing it before a solve \see SOLVER_ROUTINES_SolutionInitialiseTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_SOLUTION_INITIALISE_CURRENT_FIELD=1 !<Initialise the solution by copying in the current dependent field values \see SOLVER_ROUTINES_SolutionInitialiseTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_SOLUTION_INITIALISE_NO_CHANGE=2 !<Do not change the solution before a solve \see SOLVER_ROUTINES_SolutionInitialiseTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_OutputTypes SOLVER_ROUTINES::OutputTypes
  !> \brief The types of output
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_NO_OUTPUT=0 !<No output from the solver routines \see SOLVER_ROUTINES_OutputTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_MONITOR_OUTPUT=1 !<Monitor output from solver routines \see SOLVER_ROUTINES_OutputTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_PROGRESS_OUTPUT=2 !<Progress output from solver routines \see SOLVER_ROUTINES_OutputTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_TIMING_OUTPUT=3 !<Timing output from the solver routines plus below \see SOLVER_ROUTINES_OutputTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_SOLVER_OUTPUT=4 !<Solver specific output from the solver routines plus below \see SOLVER_ROUTINES_OutputTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_MATRIX_OUTPUT=5 !<Solver matrices output from the solver routines plus below\see SOLVER_ROUTINES_OutputTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_SparsityTypes SOLVER_ROUTINES::SparsityTypes
  !> \brief The types of sparse solver matrices
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_SPARSE_MATRICES=1 !<Use sparse solver matrices \see SOLVER_ROUTINES_SparsityTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_FULL_MATRICES=2 !<Use fully populated solver matrices \see SOLVER_ROUTINES_SparsityTypes,SOLVER_ROUTINES
  !>@}

  
  !> \addtogroup SOLVER_ROUTINES_SymmetryTypes SOLVER_ROUTINES::SymmetryTypes
  !> \brief The types of symmetry in the solver matrices
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_SYMMETRIC_MATRICES=1 !<Use symmetric solver matrices \see SOLVER_ROUTINES_SymmetryTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_UNSYMMETRIC_MATRICES=2 !<Use unsymmetric solver matrices \see SOLVER_ROUTINES_SymmetryTypes,SOLVER_ROUTINES
  !>@}

 
  !> \addtogroup PROBLEM_CONSTANTS_EquationsLinearityTypes PROBLEM_CONSTANTS::EquationsLinearityTypes
  !> \brief The solver equations linearity types 
  !> \see PROBLEM_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_LINEAR=1 !<Solver equations are linear \see PROBLEM_CONSTANTS_EquationLinearityTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_NONLINEAR=2 !<Solver equations are nonlinear \see PROBLEM_CONSTANTS_EquationLinearityTypes,PROBLEM_CONSTANTS
  !>@}

  !> \addtogroup PROBLEM_CONSTANTS_EquationsTimeDependenceTypes PROBLEM_CONSTANTS::EquationsTimeDependenceTypes
  !> \brief The solver equations time dependence types 
  !> \see PROBLEM_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_STATIC=1 !<Solver equations are static \see PROBLEM_CONSTANTS_EquationTimeDependenceTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_QUASISTATIC=2 !<Solver equations are quasistatic \see PROBLEM_CONSTANTS_EquationTimeDependenceTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC=3 !<Solver equations are first order dynamic \see PROBLEM_CONSTANTS_EquationTimeDependenceTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_SECOND_ORDER_DYNAMIC=4 !<Solver equations are second order dynamic \see PROBLEM_CONSTANTS_EquationTimeDependenceTypes,PROBLEM_CONSTANTS
  !>@}

  
  !> \addtogroup PROBLEM_CONSTANTS_CellMLEquationsLinearityTypes OpenCMISS::Iron::CellMLEquationsLinearityTypes
  !> \brief The CellML equations linearity types 
  !> \see PROBLEM_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: CELLML_EQUATIONS_LINEAR=1 !<CellML equations are linear \see PROBLEM_CONSTANTS_CellMLEquationLinearityTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: CELLML_EQUATIONS_NONLINEAR=2 !<CellML equations are nonlinear \see PROBLEM_CONSTANTS_CellMLEquationLinearityTypes,PROBLEM_CONSTANTS
  !>@}

  !> \addtogroup PROBLEM_CONSTANTS_CellMLEquationsTimeDependenceTypes OpenCMISS:Iron::CellMLEquationsTimeDependenceTypes
  !> \brief The CellML equations time dependence types 
  !> \see PROBLEM_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: CELLML_EQUATIONS_STATIC=1 !<CellML equations are static \see PROBLEM_CONSTANTS_CellMLEquationTimeDependenceTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: CELLML_EQUATIONS_QUASISTATIC=2 !<CellML equations are quasistatic \see PROBLEM_CONSTANTS_CellMLEquationTimeDependenceTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: CELLML_EQUATIONS_DYNAMIC=3 !<CellML equations are dynamic \see PROBLEM_CONSTANTS_CellMLEquationTimeDependenceTypes,PROBLEM_CONSTANTS
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

  INTERFACE SOLVER_CELLML_EQUATIONS_GET
    MODULE PROCEDURE Solver_CellMLEquationsGet
  END INTERFACE SOLVER_CELLML_EQUATIONS_GET

  INTERFACE SOLVER_SOLVER_EQUATIONS_GET
    MODULE PROCEDURE Solver_SolverEquationsGet
  END INTERFACE SOLVER_SOLVER_EQUATIONS_GET

  INTERFACE SOLVERS_SOLVER_GET
    MODULE PROCEDURE Solvers_SolverGet
  END INTERFACE SOLVERS_SOLVER_GET
  
  INTERFACE SOLVER_EQUATIONS_BOUNDARY_CONDITIONS_GET
    MODULE PROCEDURE SolverEquations_BoundaryConditionsGet
  END INTERFACE SOLVER_EQUATIONS_BOUNDARY_CONDITIONS_GET

  PUBLIC SOLVER_NUMBER_OF_SOLVER_TYPES
  
  PUBLIC SOLVER_LINEAR_TYPE,SOLVER_NONLINEAR_TYPE,SOLVER_DYNAMIC_TYPE,SOLVER_DAE_TYPE,SOLVER_EIGENPROBLEM_TYPE, &
    & SOLVER_OPTIMISER_TYPE,SOLVER_CELLML_EVALUATOR_TYPE,SOLVER_GEOMETRIC_TRANSFORMATION_TYPE

  PUBLIC SOLVER_CMISS_LIBRARY,SOLVER_PETSC_LIBRARY,SOLVER_MUMPS_LIBRARY,SOLVER_SUPERLU_LIBRARY,SOLVER_SPOOLES_LIBRARY, &
    & SOLVER_UMFPACK_LIBRARY,SOLVER_LUSOL_LIBRARY,SOLVER_ESSL_LIBRARY,SOLVER_LAPACK_LIBRARY,SOLVER_HYPRE_LIBRARY, &
    & SOLVER_PASTIX_LIBRARY

  PUBLIC SOLVER_LINEAR_DIRECT_SOLVE_TYPE,SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE
 
  PUBLIC SOLVER_DIRECT_LU,SOLVER_DIRECT_CHOLESKY,SOLVER_DIRECT_SVD

  PUBLIC SOLVER_ITERATIVE_RICHARDSON,SOLVER_ITERATIVE_CONJUGATE_GRADIENT,SOLVER_ITERATIVE_CHEBYSHEV, &
    & SOLVER_ITERATIVE_BICONJUGATE_GRADIENT,SOLVER_ITERATIVE_GMRES,SOLVER_ITERATIVE_BiCGSTAB,SOLVER_ITERATIVE_CONJGRAD_SQUARED
  
  PUBLIC SOLVER_ITERATIVE_NO_PRECONDITIONER,SOLVER_ITERATIVE_JACOBI_PRECONDITIONER,SOLVER_ITERATIVE_BLOCK_JACOBI_PRECONDITIONER, &
    & SOLVER_ITERATIVE_SOR_PRECONDITIONER,SOLVER_ITERATIVE_INCOMPLETE_CHOLESKY_PRECONDITIONER, &
    & SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER,SOLVER_ITERATIVE_ADDITIVE_SCHWARZ_PRECONDITIONER

  PUBLIC SOLVER_NONLINEAR_NEWTON,SOLVER_NONLINEAR_BFGS_INVERSE,SOLVER_NONLINEAR_SQP,SOLVER_NONLINEAR_QUASI_NEWTON

  PUBLIC SOLVER_NEWTON_LINESEARCH,SOLVER_NEWTON_TRUSTREGION

  PUBLIC SOLVER_QUASI_NEWTON_LBFGS,SOLVER_QUASI_NEWTON_GOODBROYDEN, &
    SOLVER_QUASI_NEWTON_BADBROYDEN

  PUBLIC SOLVER_QUASI_NEWTON_LINESEARCH,SOLVER_QUASI_NEWTON_TRUSTREGION
  
  PUBLIC SOLVER_QUASI_NEWTON_LINESEARCH_BASIC,SOLVER_QUASI_NEWTON_LINESEARCH_L2, &
    & SOLVER_QUASI_NEWTON_LINESEARCH_CP

  PUBLIC SOLVER_QUASI_NEWTON_RESTART_NONE,SOLVER_QUASI_NEWTON_RESTART_POWELL, &
    & SOLVER_QUASI_NEWTON_RESTART_PERIODIC

  PUBLIC SOLVER_QUASI_NEWTON_SCALE_NONE,SOLVER_QUASI_NEWTON_SCALE_SHANNO, &
    & SOLVER_QUASI_NEWTON_SCALE_LINESEARCH,SOLVER_QUASI_NEWTON_SCALE_JACOBIAN

  PUBLIC SOLVER_NEWTON_LINESEARCH_NONORMS,SOLVER_NEWTON_LINESEARCH_LINEAR,SOLVER_NEWTON_LINESEARCH_QUADRATIC, &
    & SOLVER_NEWTON_LINESEARCH_CUBIC

  PUBLIC SOLVER_NEWTON_JACOBIAN_NOT_CALCULATED,SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED, &
    & SOLVER_NEWTON_JACOBIAN_FD_CALCULATED

  PUBLIC SOLVER_NEWTON_CONVERGENCE_PETSC_DEFAULT,SOLVER_NEWTON_CONVERGENCE_ENERGY_NORM, &
    & SOLVER_NEWTON_CONVERGENCE_DIFFERENTIATED_RATIO
  
  PUBLIC SOLVER_DYNAMIC_LINEAR,SOLVER_DYNAMIC_NONLINEAR

  PUBLIC SOLVER_DYNAMIC_FIRST_ORDER,SOLVER_DYNAMIC_SECOND_ORDER

  PUBLIC SOLVER_DYNAMIC_FIRST_DEGREE,SOLVER_DYNAMIC_SECOND_DEGREE,SOLVER_DYNAMIC_THIRD_DEGREE

  PUBLIC SOLVER_DYNAMIC_EULER_SCHEME,SOLVER_DYNAMIC_BACKWARD_EULER_SCHEME,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME, &
    & SOLVER_DYNAMIC_GALERKIN_SCHEME,SOLVER_DYNAMIC_ZLAMAL_SCHEME,SOLVER_DYNAMIC_SECOND_DEGREE_GEAR_SCHEME, &
    & SOLVER_DYNAMIC_SECOND_DEGREE_LINIGER1_SCHEME,SOLVER_DYNAMIC_SECOND_DEGREE_LINIGER2_SCHEME, &
    & SOLVER_DYNAMIC_NEWMARK1_SCHEME,SOLVER_DYNAMIC_NEWMARK2_SCHEME,SOLVER_DYNAMIC_NEWMARK3_SCHEME, &
    & SOLVER_DYNAMIC_THIRD_DEGREE_GEAR_SCHEME,SOLVER_DYNAMIC_THIRD_DEGREE_LINIGER1_SCHEME, &
    & SOLVER_DYNAMIC_THIRD_DEGREE_LINIGER2_SCHEME,SOLVER_DYNAMIC_HOUBOLT_SCHEME,SOLVER_DYNAMIC_WILSON_SCHEME, &
    & SOLVER_DYNAMIC_BOSSAK_NEWMARK1_SCHEME,SOLVER_DYNAMIC_BOSSAK_NEWMARK2_SCHEME, &
    & SOLVER_DYNAMIC_HILBERT_HUGHES_TAYLOR1_SCHEME,SOLVER_DYNAMIC_HILBERT_HUGHES_TAYLOR2_SCHEME, &
    & SOLVER_DYNAMIC_USER_DEFINED_SCHEME

  PUBLIC SOLVER_DYNAMIC_PREVIOUS_STARTUP_TYPE,SOLVER_DYNAMIC_PREVIOUS_RESTART_TYPE

  PUBLIC SOLVER_DAE_DIFFERENTIAL_ONLY,SOLVER_DAE_INDEX_1,SOLVER_DAE_INDEX_2,SOLVER_DAE_INDEX_3

  PUBLIC SOLVER_DAE_EULER,SOLVER_DAE_CRANK_NICOLSON,SOLVER_DAE_RUNGE_KUTTA,SOLVER_DAE_ADAMS_MOULTON,SOLVER_DAE_BDF, &
    & SOLVER_DAE_RUSH_LARSON,SOLVER_DAE_EXTERNAL

  PUBLIC SOLVER_DAE_EULER_FORWARD,SOLVER_DAE_EULER_BACKWARD,SOLVER_DAE_EULER_IMPROVED

  PUBLIC SOLVER_OPTIMISER_NO_OBJECTIVE,SOLVER_OPTIMISER_ONE_OBJECTIVE,SOLVER_OPTIMISER_MANY_OBJECTIVE

  PUBLIC SOLVER_OPTIMISER_CONTINUOUS_VARIABLES,SOLVER_OPTIMISER_DISCRETE_VARIABLES

  PUBLIC SOLVER_OPTIMISER_DETERMINISTIC_CERTAINTY,SOLVER_OPTIMISER_STOCHASTIC_CERTAINTY
  
  PUBLIC SOLVER_OPTIMISER_UNCONSTRAINED,SOLVER_OPTIMISER_BOUND_CONSTRAINED,SOLVER_OPTIMISER_LINEAR_CONSTRAINTS, &
    & SOLVER_OPTIMISER_NONLINEAR_CONSTRAINTS,SOLVER_OPTIMISER_PDE_CONSTRAINTS

  PUBLIC SOLVER_OPTIMISER_GRADIENT_NOT_CALCULATED,SOLVER_OPTIMISER_GRADIENT_EQUATIONS_CALCULATED, &
    & SOLVER_OPTIMISER_GRADIENT_FD_CALCULATED

  PUBLIC SOLVER_OPTIMISER_HESSIAN_NOT_CALCULATED,SOLVER_OPTIMISER_HESSIAN_EQUATIONS_CALCULATED, &
    & SOLVER_OPTIMISER_HESSIAN_FD_CALCULATED

  PUBLIC SOLVER_SOLUTION_INITIALISE_ZERO,SOLVER_SOLUTION_INITIALISE_CURRENT_FIELD,SOLVER_SOLUTION_INITIALISE_NO_CHANGE
  
  PUBLIC SOLVER_NO_OUTPUT,SOLVER_MONITOR_OUTPUT,SOLVER_PROGRESS_OUTPUT,SOLVER_TIMING_OUTPUT,SOLVER_SOLVER_OUTPUT, &
    & SOLVER_MATRIX_OUTPUT
  
  PUBLIC SOLVER_SPARSE_MATRICES,SOLVER_FULL_MATRICES

  PUBLIC SOLVER_SYMMETRIC_MATRICES,SOLVER_UNSYMMETRIC_MATRICES

  PUBLIC SOLVER_EQUATIONS_LINEAR,SOLVER_EQUATIONS_NONLINEAR

  PUBLIC SOLVER_EQUATIONS_STATIC,SOLVER_EQUATIONS_QUASISTATIC,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC, &
    & SOLVER_EQUATIONS_SECOND_ORDER_DYNAMIC

  PUBLIC CELLML_EQUATIONS_LINEAR,CELLML_EQUATIONS_NONLINEAR

  PUBLIC CELLML_EQUATIONS_STATIC,CELLML_EQUATIONS_QUASISTATIC,CELLML_EQUATIONS_DYNAMIC
  
  PUBLIC CellMLEquations_SolverGet

  PUBLIC Solver_AssertIsFinished,Solver_AssertNotFinished

  PUBLIC Solver_AssertIsDynamic,Solver_AssertIsLinear,Solver_AssertIsNonlinear,Solver_AssertIsOptimiser
  
  PUBLIC Solver_CellMLEquationsGet

  PUBLIC SOLVER_CELLML_EQUATIONS_GET

  PUBLIC Solver_ControlLoopGet

  PUBLIC Solver_DynamicSolverGet
  
  PUBLIC Solver_LinearSolverGet
  
  PUBLIC Solver_NonlinearSolverGet
  
  PUBLIC Solver_OptimiserSolverGet
  
  PUBLIC Solver_SolverEquationsGet

  PUBLIC SOLVER_SOLVER_EQUATIONS_GET

  PUBLIC Solver_LinkingSolverGet

  PUBLIC Solver_SolversGet

  PUBLIC Solver_WorkGroupGet

  PUBLIC Solvers_ControlLoopGet

  PUBLIC Solvers_SolverGet

  PUBLIC SOLVERS_SOLVER_GET

  PUBLIC SolverDynamic_LinearSolverGet

  PUBLIC SolverDynamic_NonlinearSolverGet

  PUBLIC SolverDynamic_SolverGet

  PUBLIC SolverEquations_AssertIsFinished,SolverEquations_AssertNotFinished

  PUBLIC SolverEquations_BoundaryConditionsGet

  PUBLIC SOLVER_EQUATIONS_BOUNDARY_CONDITIONS_GET
  
  PUBLIC SolverEquations_SolverGet

  PUBLIC SolverEquations_SolverMappingGet

  PUBLIC SolverEquations_SolverMatricesGet

  PUBLIC SolverLinear_SolverGet
  
  PUBLIC SolverNonlinear_LinearSolverGet
  
  PUBLIC SolverNonlinear_NewtonSolverGet
  
  PUBLIC SolverNonlinear_QuasiNewtonSolverGet
  
  PUBLIC SolverNonlinear_SolverGet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solver for a CellML equations.
  SUBROUTINE CellMLEquations_SolverGet(cellMLEquations,solver,err,error,*)

    TYPE(CELLML_EQUATIONS_TYPE), POINTER :: cellMLEquations !<A pointer to the CellML equations to get the solver for
    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<On exit, a pointer to the solver for the specified CellML equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("CellMLEquations_SolverGet",err,error,*998)

    IF(ASSOCIATED(solver)) CALL FlagError("Solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellMLEquations)) CALL FlagError("CellML equations is not associated.",err,error,*999)

    solver=>cellMLEquations%solver
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("CellML equations solver is not associated.",err,error,*999)
      
    EXITS("CellMLEquations_SolverGet")
    RETURN
999 NULLIFY(solver)
998 ERRORSEXITS("CellMLEquations_SolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLEquations_SolverGet
  
  !
  !================================================================================================================================
  !

  !>Assert that a solver has been finished
  SUBROUTINE Solver_AssertIsFinished(solver,err,error,*)

    !Argument Variables
    TYPE(SOLVER_TYPE), POINTER, INTENT(IN) :: solver !<The solver to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Solver_AssertIsFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)

    IF(.NOT.solver%SOLVER_FINISHED) CALL FlagError("Solver has not been finished.",err,error,*999)
    
    EXITS("Solver_AssertIsFinished")
    RETURN
999 ERRORSEXITS("Solver_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_AssertIsFinished

  !
  !================================================================================================================================
  !

  !>Assert that a solver has not been finished
  SUBROUTINE Solver_AssertNotFinished(solver,err,error,*)

    !Argument Variables
    TYPE(SOLVER_TYPE), POINTER, INTENT(IN) :: solver !<The solver to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_AssertNotFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)

    IF(solver%SOLVER_FINISHED) CALL FlagError("Solver has already been finished.",err,error,*999)
    
    EXITS("Solver_AssertNotFinished")
    RETURN
999 ERRORSEXITS("Solver_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Assert that a solver is a dynamic solver
  SUBROUTINE Solver_AssertIsDynamic(solver,err,error,*)

    !Argument Variables
    TYPE(SOLVER_TYPE), POINTER, INTENT(INOUT) :: solver !<The solver to assert the dynamic solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Solver_AssertIsDynamic",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("solver is not associated.",err,error,*999)

    IF(solver%SOLVE_TYPE/=SOLVER_DYNAMIC_TYPE) THEN
      localError="The solver type of "//TRIM(NumberToVString(solver%SOLVE_TYPE,"*",err,error))// &
        & " does not correspond to the required dynamic solve type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Solver_AssertIsDynamic")
    RETURN
999 ERRORSEXITS("Solver_AssertIsDynamic",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_AssertIsDynamic

  !
  !================================================================================================================================
  !

  !>Assert that a solver is a linear solver
  SUBROUTINE Solver_AssertIsLinear(solver,err,error,*)

    !Argument Variables
    TYPE(SOLVER_TYPE), POINTER, INTENT(INOUT) :: solver !<The solver to assert the linear solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Solver_AssertIsLinear",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("solver is not associated.",err,error,*999)

    IF(solver%SOLVE_TYPE/=SOLVER_LINEAR_TYPE) THEN
      localError="The solver type of "//TRIM(NumberToVString(solver%SOLVE_TYPE,"*",err,error))// &
        & " does not correspond to the required linear solve type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Solver_AssertIsLinear")
    RETURN
999 ERRORSEXITS("Solver_AssertIsLinear",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_AssertIsLinear

  !
  !================================================================================================================================
  !

  !>Assert that a solver is a nonlinear solver
  SUBROUTINE Solver_AssertIsNonlinear(solver,err,error,*)

    !Argument Variables
    TYPE(SOLVER_TYPE), POINTER, INTENT(INOUT) :: solver !<The solver to assert the nonlinear solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Solver_AssertIsNonlinear",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("solver is not associated.",err,error,*999)

    IF(solver%SOLVE_TYPE/=SOLVER_NONLINEAR_TYPE) THEN
      localError="The solver type of "//TRIM(NumberToVString(solver%SOLVE_TYPE,"*",err,error))// &
        & " does not correspond to the required nonlinear solve type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Solver_AssertIsNonlinear")
    RETURN
999 ERRORSEXITS("Solver_AssertIsNonlinear",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_AssertIsNonlinear

  !
  !================================================================================================================================
  !

  !>Assert that a solver is an optimiser solver
  SUBROUTINE Solver_AssertIsOptimiser(solver,err,error,*)

    !Argument Variables
    TYPE(SOLVER_TYPE), POINTER, INTENT(INOUT) :: solver !<The solver to assert the optimiser solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Solver_AssertIsOptimiser",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("solver is not associated.",err,error,*999)

    IF(solver%SOLVE_TYPE/=SOLVER_OPTIMISER_TYPE) THEN
      localError="The solver type of "//TRIM(NumberToVString(solver%SOLVE_TYPE,"*",err,error))// &
        & " does not correspond to the required optimiser solve type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Solver_AssertIsOptimiser")
    RETURN
999 ERRORSEXITS("Solver_AssertIsOptimiser",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_AssertIsOptimiser

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the CellML equations for a solver. \see OpenCMISS::Iron::cmfe_Solver_CellMLEquationsGet
  SUBROUTINE Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver to get the CellML equations for
    TYPE(CELLML_EQUATIONS_TYPE), POINTER :: cellMLEquations !<On exit, a pointer to the specified CellML equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_CellMLEquationsGet",err,error,*998)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)
    !IF(.NOT.solver%SOLVER_FINISHED) CALL FlagError("Solver has not been finished.",err,error,*998)
    IF(ASSOCIATED(cellMLEquations)) CALL FlagError("CellML equations is already associated.",err,error,*998)

    cellMLEquations=>solver%CELLML_EQUATIONS
    IF(.NOT.ASSOCIATED(cellMLEquations)) CALL FlagError("Solver CellML equations is not associated.",err,error,*999)
      
    EXITS("Solver_CellMLEquationsGet")
    RETURN
999 NULLIFY(cellMLEquations)
998 ERRORSEXITS("Solver_CellMLEquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_CellMLEquationsGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the control loop for a solver.
  SUBROUTINE Solver_ControlLoopGet(solver,controlLoop,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver to get the control loop for
    TYPE(ControlLoopType), POINTER :: controlLoop !<On exit, a pointer to the control loop for the specified solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(SOLVERS_TYPE), POINTER :: solvers
 
    ENTERS("Solver_ControlLoopGet",err,error,*998)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)
    IF(.NOT.solver%SOLVER_FINISHED) CALL FlagError("Solver has not been finished.",err,error,*998)
    IF(ASSOCIATED(controlLoop)) CALL FlagError("Control loop is already associated.",err,error,*998)

    solvers=>solver%solvers
    IF(.NOT.ASSOCIATED(solvers)) THEN
      IF(ASSOCIATED(solver%linking_solver)) THEN
        solvers=>solver%linking_solver%solvers
        IF(.NOT.ASSOCIATED(solvers)) CALL FlagError("Solver linking solver solvers is not associated.",err,error,*999)
      ELSE
        CALL FlagError("Solver solvers is not associated.",err,error,*999)
      ENDIF
    ENDIF
    controlLoop=>solvers%controlLoop
    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Solvers control loop is not associated.",err,error,*999)
    
    EXITS("Solver_ControlLoopGet")
    RETURN
999 NULLIFY(controlLoop)
998 ERRORSEXITS("Solver_ControlLoopGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_ControlLoopGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the dynamic solver for a solver.
  SUBROUTINE Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver to get the dynamic solver for
    TYPE(DYNAMIC_SOLVER_TYPE), POINTER :: dynamicSolver !<On exit, a pointer to the dynamic solver for the specified solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_DynamicSolverGet",err,error,*998)

    IF(ASSOCIATED(dynamicSolver)) CALL FlagError("Dynamic solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)

    dynamicSolver=>solver%DYNAMIC_SOLVER
    IF(.NOT.ASSOCIATED(dynamicSolver)) CALL FlagError("Solver dynamic solver is not associated.",err,error,*999)
    
    EXITS("Solver_DynamicSolverGet")
    RETURN
999 NULLIFY(dynamicSolver)
998 ERRORSEXITS("Solver_DynamicSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_DynamicSolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the linear solver for a solver.
  SUBROUTINE Solver_LinearSolverGet(solver,linearSolver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver to get the linear solver for
    TYPE(LINEAR_SOLVER_TYPE), POINTER :: linearSolver !<On exit, a pointer to the linear solver for the specified solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_LinearSolverGet",err,error,*998)

    IF(ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)

    linearSolver=>solver%LINEAR_SOLVER
    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("Solver linear solver is not associated.",err,error,*999)
    
    EXITS("Solver_LinearSolverGet")
    RETURN
999 NULLIFY(linearSolver)
998 ERRORSEXITS("Solver_LinearSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_LinearSolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the nonlinear solver for a solver.
  SUBROUTINE Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver to get the nonlinear solver for
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: nonlinearSolver !<On exit, a pointer to the nonlinear solver for the specified solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_NonlinearSolverGet",err,error,*998)

    IF(ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)

    nonlinearSolver=>solver%NONLINEAR_SOLVER
    IF(.NOT.ASSOCIATED(nonlinearSolver)) CALL FlagError("Solver nonlinear solver is not associated.",err,error,*999)
    
    EXITS("Solver_NonlinearSolverGet")
    RETURN
999 NULLIFY(nonlinearSolver)
998 ERRORSEXITS("Solver_NonlinearSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_NonlinearSolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the optimiser solver for a solver.
  SUBROUTINE Solver_OptimiserSolverGet(solver,optimiserSolver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver to get the optimiser solver for
    TYPE(OptimiserSolverType), POINTER :: optimiserSolver !<On exit, a pointer to the optimiser solver for the specified solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_OptimiserSolverGet",err,error,*998)

    IF(ASSOCIATED(optimiserSolver)) CALL FlagError("Optimiser solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)

    optimiserSolver=>solver%optimiserSolver
    IF(.NOT.ASSOCIATED(optimiserSolver)) CALL FlagError("Solver optimiser solver is not associated.",err,error,*999)
    
    EXITS("Solver_OptimiserSolverGet")
    RETURN
999 NULLIFY(optimiserSolver)
998 ERRORSEXITS("Solver_OptimiserSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_OptimiserSolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solver equations for a solver. \see OpenCMISS::Iron::cmfe_Solver_SolverEquationsGet
  SUBROUTINE Solver_SolverEquationsGet(solver,solverEquations,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver to get the solver equations for
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<On exit, a pointer to the specified solver equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_SolverEquationsGet",err,error,*998)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)
    IF(.NOT.solver%SOLVER_FINISHED) CALL FlagError("Solver has not been finished.",err,error,*998)
    IF(ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is already associated.",err,error,*998)

    solverEquations=>solver%SOLVER_EQUATIONS
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver solver equations is not associated.",err,error,*999)
      
    EXITS("Solver_SolverEquationsGet")
    RETURN
999 NULLIFY(solverEquations)
998 ERRORSEXITS("Solver_SolverEquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_SolverEquationsGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the linking solver for a solver.
  SUBROUTINE Solver_LinkingSolverGet(solver,linkingSolver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver to get the solvers for.
    TYPE(SOLVER_TYPE), POINTER :: linkingSolver !<On exit, A pointer to the linking solver for the solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_LinkingSolverGet",err,error,*998)

    IF(ASSOCIATED(linkingSolver)) CALL FlagError("Linking solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
      
    linkingSolver=>solver%LINKING_SOLVER
    IF(.NOT.ASSOCIATED(linkingSolver)) CALL FlagError("The solver linking solver is not associated.",err,error,*999)
       
    EXITS("Solver_LinkingSolverGet")
    RETURN
999 NULLIFY(linkingSolver)
998 ERRORSEXITS("Solver_LinkingSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_LinkingSolverGet

  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the solvers for a solver.
  SUBROUTINE Solver_SolversGet(solver,solvers,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver to get the solvers for.
    TYPE(SOLVERS_TYPE), POINTER :: solvers !<On exit, A pointer to the solvers for the solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: testSolver
 
    ENTERS("Solver_SolversGet",err,error,*998)

    IF(ASSOCIATED(solvers)) CALL FlagError("Solvers is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)

    testSolver=>solver
    DO WHILE(ASSOCIATED(testSolver%LINKING_SOLVER))
      testSolver=>testSolver%LINKING_SOLVER
    ENDDO
    solvers=>testSolver%solvers
    IF(.NOT.ASSOCIATED(solvers)) CALL FlagError("The solver solvers is not associated.",err,error,*999)
       
    EXITS("Solver_SolversGet")
    RETURN
999 NULLIFY(solvers)
998 ERRORSEXITS("Solver_SolversGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_SolversGet

  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the work group for a solver. FOR NOW JUST RETURN THE PROBLEM WORK GROUP.
  SUBROUTINE Solver_WorkGroupGet(solver,workGroup,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver to get the work group for.
    TYPE(WorkGroupType), POINTER :: workGroup !<On exit, A pointer to the work group for the solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(SOLVERS_TYPE), POINTER :: solvers
 
    ENTERS("Solver_WorkGroupGet",err,error,*998)

    IF(ASSOCIATED(workGroup)) CALL FlagError("Work group is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)

    NULLIFY(solvers)
    CALL Solver_SolversGet(solver,solvers,err,error,*999)
    IF(.NOT.ASSOCIATED(solvers)) CALL FlagError("Solver solvers is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(solvers%controlLoop)) &
      & CALL FlagError("Solver solvers control loop is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(solvers%controlLoop%problem)) &
      & CALL FlagError("Solver solvers control loop problem is not associated.",err,error,*999)

    workGroup=>solvers%controlLoop%problem%workGroup
    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("The solver work group is not associated.",err,error,*999)
       
    EXITS("Solver_WorkGroupGet")
    RETURN
999 NULLIFY(workGroup)
998 ERRORSEXITS("Solver_WorkGroupGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_WorkGroupGet

  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the control loop for a solvers.
  SUBROUTINE Solvers_ControlLoopGet(solvers,controlLoop,err,error,*)

    !Argument variables
    TYPE(SOLVERS_TYPE), POINTER :: solvers !<A pointer to the solvers to get the control loop for.
    TYPE(ControlLoopType), POINTER :: controlLoop !<On exit, A pointer to the control loop for the solvers. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solvers_ControlLoopGet",err,error,*998)

    IF(ASSOCIATED(controlLoop)) CALL FlagError("Control loop is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solvers)) CALL FlagError("Solvers is not associated.",err,error,*999)
      
    controlLoop=>solvers%controlLoop
    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("The solvers control loop is not associated.",err,error,*999)
       
    EXITS("Solvers_ControlLoopGet")
    RETURN
999 NULLIFY(controlLoop)
998 ERRORSEXITS("Solvers_ControlLoopGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solvers_ControlLoopGet

  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the specified solver in the list of solvers.
  SUBROUTINE Solvers_SolverGet(solvers,solverIndex,solver,err,error,*)

    !Argument variables
    TYPE(SOLVERS_TYPE), POINTER :: solvers !<A pointer to the solvers to get the solver for
    INTEGER(INTG), INTENT(IN) :: solverIndex !<The specified solver to get
    TYPE(SOLVER_TYPE), POINTER :: solver !<On exit, a pointer to the specified solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Solvers_SolverGet",err,error,*998)

    IF(.NOT.ASSOCIATED(solvers)) CALL FlagError("Solvers is not associated.",err,error,*998)
    IF(ASSOCIATED(solver)) CALL FlagError("Solver is already associated.",err,error,*998)
    IF(solverIndex<=0.OR.solverIndex>solvers%NUMBER_OF_SOLVERS) THEN
      localError="The specified solver index of "//TRIM(NumberToVString(solverIndex,"*",err,error))// &
        & " is invalid. The solver index must be >= 1 and <= "// &
        & TRIM(NumberToVString(solvers%NUMBER_OF_SOLVERS,"*",err,error))//"."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    IF(.NOT.ALLOCATED(solvers%solvers)) CALL FlagError("Solvers solvers is not associated.",err,error,*998)
      
    solver=>solvers%solvers(solverIndex)%ptr
    IF(.NOT.ASSOCIATED(solver)) THEN
      localError="The solvers solver is not associated for solver index "// &
        & TRIM(NumberToVString(solverIndex,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("Solvers_SolverGet")
    RETURN
999 NULLIFY(solver)
998 ERRORSEXITS("Solvers_SolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solvers_SolverGet

  !
  !================================================================================================================================
  !

  !>Gets the linear solver for a dynamic solver. 
  SUBROUTINE SolverDynamic_LinearSolverGet(dynamicSolver,linearSolver,err,error,*)

    !Argument variables
    TYPE(DYNAMIC_SOLVER_TYPE), POINTER :: dynamicSolver !<A pointer to the dynamic solver to get the linear solver for
    TYPE(SOLVER_TYPE), POINTER :: linearSolver !<On exit, a pointer to the linear solver for the specified dynamic solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverDynamic_LinearSolverGet",err,error,*998)

    IF(ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dynamicSolver)) CALL FlagError("Dynamic solver is not associated.",err,error,*999)

    linearSolver=>dynamicSolver%LINEAR_SOLVER
    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("Dynamic solver linear solver is not associated.",err,error,*999)
 
    EXITS("SolverDynamic_LinearSolverGet")
    RETURN
998 NULLIFY(linearSolver)
999 ERRORSEXITS("SolverDynamic_LinearSolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverDynamic_LinearSolverGet
     
  !
  !================================================================================================================================
  !

  !>Gets the nonlinear solver for a dynamic solver. 
  SUBROUTINE SolverDynamic_NonlinearSolverGet(dynamicSolver,nonlinearSolver,err,error,*)

    !Argument variables
    TYPE(DYNAMIC_SOLVER_TYPE), POINTER :: dynamicSolver !<A pointer to the dynamic solver to get the nonlinear solver for
    TYPE(SOLVER_TYPE), POINTER :: nonlinearSolver !<On exit, a pointer to the nonlinear solver for the specified dynamic solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverDynamic_NonlinearSolverGet",err,error,*998)

    IF(ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dynamicSolver)) CALL FlagError("Dynamic solver is not associated.",err,error,*999)

    nonlinearSolver=>dynamicSolver%NONLINEAR_SOLVER
    IF(.NOT.ASSOCIATED(nonlinearSolver)) CALL FlagError("Dynamic solver nonlinear solver is not associated.",err,error,*999)
 
    EXITS("SolverDynamic_NonlinearSolverGet")
    RETURN
998 NULLIFY(nonlinearSolver)
999 ERRORSEXITS("SolverDynamic_NonlinearSolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverDynamic_NonlinearSolverGet
     
  !
  !================================================================================================================================
  !

  !>Gets the solver for a dynamic solver. 
  SUBROUTINE SolverDynamic_SolverGet(dynamicSolver,solver,err,error,*)

    !Argument variables
    TYPE(DYNAMIC_SOLVER_TYPE), POINTER :: dynamicSolver !<A pointer to the dynamic solver to get the solver for
    TYPE(SOLVER_TYPE), POINTER :: solver !<On exit, a pointer to the solver for the specified dynamic solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverDynamic_SolverGet",err,error,*998)

    IF(ASSOCIATED(solver)) CALL FlagError("Solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dynamicSOlver)) CALL FlagError("Dynamic solver is not associated.",err,error,*999)

    solver=>dynamicSolver%solver
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Dynamic solver solver is not associated.",err,error,*999)
 
    EXITS("SolverDynamic_SolverGet")
    RETURN
998 NULLIFY(solver)
999 ERRORSEXITS("SolverDynamic_SolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverDynamic_SolverGet
     
  !
  !================================================================================================================================
  !

  !>Assert that a solver equations has been finished
  SUBROUTINE SolverEquations_AssertIsFinished(solverEquations,err,error,*)

    !Argument Variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER, INTENT(IN) :: solverEquations !<The solver equations to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverEquations_AssertIsFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)

    IF(.NOT.solverEquations%SOLVER_EQUATIONS_FINISHED) CALL FlagError("Solver equations has not been finished.",err,error,*999)
    
    EXITS("SolverEquations_AssertIsFinished")
    RETURN
999 ERRORSEXITS("SolverEquations_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE SolverEquations_AssertIsFinished

  !
  !================================================================================================================================
  !

  !>Assert that a solver equations has not been finished
  SUBROUTINE SolverEquations_AssertNotFinished(solverEquations,err,error,*)

    !Argument Variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER, INTENT(IN) :: solverEquations !<The solver equations to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverEquations_AssertNotFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)

    IF(solverEquations%SOLVER_EQUATIONS_FINISHED) CALL FlagError("Solver equations has already been finished.",err,error,*999)
    
    EXITS("SolverEquations_AssertNotFinished")
    RETURN
999 ERRORSEXITS("SolverEquations_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE SolverEquations_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Gets the boundary conditions for solver equations. \see OpenCMISS::Iron::cmfe_SolverEquations_BoundaryConditionsGet
  SUBROUTINE SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<A pointer to the solver equations to get the boundary conditions for
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: boundaryConditions !<On exit, a pointer to the boundary conditions for the specified solver equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEquations_BoundaryConditionsGet",err,error,*999)

    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
    IF(.NOT.solverEquations%SOLVER_EQUATIONS_FINISHED) CALL FlagError("Solver equations has not been finished.",err,error,*999)
    IF(ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is already associated.",err,error,*999)

    boundaryConditions=>solverEquations%BOUNDARY_CONDITIONS
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Solver equations boundary conditions is not associated.", &
      & err,error,*999)
 
    EXITS("SolverEquations_BoundaryConditionsGet")
    RETURN
999 ERRORSEXITS("SolverEquations_BoundaryConditionsGet",err,error)
    RETURN 1

  END SUBROUTINE SolverEquations_BoundaryConditionsGet
     
  !
  !================================================================================================================================
  !

  !>Gets the solver for solver equations. 
  SUBROUTINE SolverEquations_SolverGet(solverEquations,solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<A pointer to the solver equations to get the solver for
    TYPE(SOLVER_TYPE), POINTER :: solver !<On exit, a pointer to the solver for the specified solver equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEquations_SolverGet",err,error,*998)

    IF(ASSOCIATED(solver)) CALL FlagError("Solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)

    solver=>solverEquations%solver
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver equations solver is not associated.",err,error,*999)
 
    EXITS("SolverEquations_SolverGet")
    RETURN
998 NULLIFY(solver)
999 ERRORSEXITS("SolverEquations_SolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverEquations_SolverGet
     
  !
  !================================================================================================================================
  !

  !>Gets the solver mapping for solver equations. 
  SUBROUTINE SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<A pointer to the solver equations to get the solver mapping for
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping !<On exit, a pointer to the solver mapping for the specified solver equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEquations_SolverMappingGet",err,error,*998)

    IF(ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
    IF(.NOT.solverEquations%SOLVER_EQUATIONS_FINISHED) CALL FlagError("Solver equations has not been finished.",err,error,*999)

    solverMapping=>solverEquations%solverMapping
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver equations solver mapping is not associated.",err,error,*999)
 
    EXITS("SolverEquations_SolverMappingGet")
    RETURN
998 NULLIFY(solverMapping)
999 ERRORSEXITS("SolverEquations_SolverMappingGet",err,error)
    RETURN 1

  END SUBROUTINE SolverEquations_SolverMappingGet
     
  !
  !================================================================================================================================
  !

  !>Gets the solver matrices for solver equations. 
  SUBROUTINE SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<A pointer to the solver equations to get the solver matrices for
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: solverMatrices !<On exit, a pointer to the solver matrices for the specified solver equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEquations_SolverMatricesGet",err,error,*998)

    IF(ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
    IF(.NOT.solverEquations%SOLVER_EQUATIONS_FINISHED) CALL FlagError("Solver equations has not been finished.",err,error,*999)

    solverMatrices=>solverEquations%SOLVER_MATRICES
    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver equations solver matrices is not associated.",err,error,*999)
 
    EXITS("SolverEquations_SolverMatricesGet")
    RETURN
998 NULLIFY(solverMatrices)
999 ERRORSEXITS("SolverEquations_SolverMatricesGet",err,error)
    RETURN 1

  END SUBROUTINE SolverEquations_SolverMatricesGet
     
  !
  !================================================================================================================================
  !

  !>Gets the solver for a linear solver. 
  SUBROUTINE SolverLinear_SolverGet(linearSolver,solver,err,error,*)

    !Argument variables
    TYPE(LINEAR_SOLVER_TYPE), POINTER :: linearSolver !<A pointer to the linear solver to get the solver for
    TYPE(SOLVER_TYPE), POINTER :: solver !<On exit, a pointer to the solver for the specified linear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverLinear_SolverGet",err,error,*998)

    IF(ASSOCIATED(solver)) CALL FlagError("Solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is not associated.",err,error,*999)

    solver=>linearSolver%solver
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Linear solver solver is not associated.",err,error,*999)
 
    EXITS("SolverLinear_SolverGet")
    RETURN
998 NULLIFY(solver)
999 ERRORSEXITS("SolverLinear_SolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverLinear_SolverGet
     
  !
  !================================================================================================================================
  !

  !>Gets the linear solver for a nonlinear solver. 
  SUBROUTINE SolverNonlinear_LinearSolverGet(nonlinearSolver,linearSolver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: nonlinearSolver !<A pointer to the nonlinear solver to get the linear solver for
    TYPE(SOLVER_TYPE), POINTER :: linearSolver !<On exit, a pointer to the linear solver for the specified nonlinear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: newtonSolver
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: nonlinSolver
    TYPE(QUASI_NEWTON_SOLVER_TYPE), POINTER :: quasiNewtonSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverNonlinear_SolverLinearSolverGet",err,error,*998)

    IF(ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is not associated.",err,error,*999)

    nonlinSolver=>nonlinearSolver%NONLINEAR_SOLVER
    IF(.NOT.ASSOCIATED(nonlinSolver)) CALL FlagError("Nonlinear solver nonlinear solver is not associated.",err,error,*999)
    SELECT CASE(nonlinSolver%NONLINEAR_SOLVE_TYPE)
    CASE(SOLVER_NONLINEAR_NEWTON)
      newtonSolver=>nonlinSolver%NEWTON_SOLVER
      IF(.NOT.ASSOCIATED(newtonSolver)) CALL FlagError("Nonlinear solver Newton solver is not associated.",err,error,*999)
      linearSolver=>newtonSolver%LINEAR_SOLVER
    CASE(SOLVER_NONLINEAR_BFGS_INVERSE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_NONLINEAR_SQP)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
      quasiNewtonSolver=>nonlinSolver%QUASI_NEWTON_SOLVER
      IF(.NOT.ASSOCIATED(quasiNewtonSolver)) &
        & CALL FlagError("Nonlinear solver quasi-Newton solver is not associated.",err,error,*999)
      linearSolver=>quasiNewtonSolver%LINEAR_SOLVER
    CASE DEFAULT
      localError="The nonlinear solver type of "//TRIM(NumberToVString(nonlinSolver%NONLINEAR_SOLVE_TYPE,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("Nonlinear solver linear solver is not associated.",err,error,*999)
 
    EXITS("SolverNonlinear_LinearSolverGet")
    RETURN
998 NULLIFY(linearSolver)
999 ERRORSEXITS("SolverNonlinear_LinearSolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverNonlinear_LinearSolverGet
     
  !
  !================================================================================================================================
  !

  !>Gets the newton solver for a nonlinear solver. 
  SUBROUTINE SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*)

    !Argument variables
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: nonlinearSolver !<A pointer to the nonlinear solver to get the Newton solver for
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: newtonSolver !<On exit, a pointer to the newton solver for the specified nonlinear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinear_NewtonSolverGet",err,error,*998)

    IF(ASSOCIATED(newtonSolver)) CALL FlagError("Newton solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is not associated.",err,error,*999)

    newtonSolver=>nonlinearSolver%NEWTON_SOLVER
    IF(.NOT.ASSOCIATED(newtonSolver)) CALL FlagError("Nonlinear solver Newton solver is not associated.",err,error,*999)
 
    EXITS("SolverNonlinear_NewtonSolverGet")
    RETURN
998 NULLIFY(newtonSolver)
999 ERRORSEXITS("SolverNonlinear_NewtonSolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverNonlinear_NewtonSolverGet
     
  !
  !================================================================================================================================
  !

  !>Gets the quasi Newton solver for a nonlinear solver. 
  SUBROUTINE SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*)

    !Argument variables
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: nonlinearSolver !<A pointer to the nonlinear solver to get the quasi Newton solver for
    TYPE(QUASI_NEWTON_SOLVER_TYPE), POINTER :: quasiNewtonSolver !<On exit, a pointer to the quasi Newton solver for the specified nonlinear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinear_QuasiNewtonSolverGet",err,error,*998)

    IF(ASSOCIATED(quasiNewtonSolver)) CALL FlagError("Quasi Newton solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is not associated.",err,error,*999)

    quasiNewtonSolver=>nonlinearSolver%QUASI_NEWTON_SOLVER
    IF(.NOT.ASSOCIATED(quasiNewtonSolver)) CALL FlagError("Nonlinear solver quasi Newton solver is not associated.",err,error,*999)
 
    EXITS("SolverNonlinear_QuasiNewtonSolverGet")
    RETURN
998 NULLIFY(quasiNewtonSolver)
999 ERRORSEXITS("SolverNonlinear_QuasiNewtonSolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverNonlinear_QuasiNewtonSolverGet
     
  !
  !================================================================================================================================
  !

  !>Gets the solver for a non-linear solver. 
  SUBROUTINE SolverNonlinear_SolverGet(nonlinearSolver,solver,err,error,*)

    !Argument variables
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: nonlinearSolver !<A pointer to the nonlinear solver to get the solver for
    TYPE(SOLVER_TYPE), POINTER :: solver !<On exit, a pointer to the solver for the specified nonlinear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinear_SolverGet",err,error,*998)

    IF(ASSOCIATED(solver)) CALL FlagError("Solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is not associated.",err,error,*999)

    solver=>nonlinearSolver%solver
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Nonlinear solver solver is not associated.",err,error,*999)
 
    EXITS("SolverNonlinear_SolverGet")
    RETURN
998 NULLIFY(solver)
999 ERRORSEXITS("SolverNonlinear_SolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverNonlinear_SolverGet
     
  !
  !================================================================================================================================
  !

END MODULE SolverAccessRoutines
