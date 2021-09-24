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

  !> \addtogroup SolverRoutines_SolverTypes SolverRoutines::SolverTypes
  !> \brief The types of solver
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_NUMBER_OF_SOLVER_TYPES=9 !<Number of different solver types possible \see SolverRoutines_SolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_LINEAR_TYPE=1 !<A linear solver \see SolverRoutines_SolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_NONLINEAR_TYPE=2 !<A nonliXnear solver  \see SolverRoutines_SolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_TYPE=3 !<A dynamic solver \see SolverRoutines_SolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_TYPE=4 !<A differential-algebraic equation solver \see SolverRoutines_SolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_EIGENPROBLEM_TYPE=5 !<A eigenproblem solver \see SolverRoutines_SolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_TYPE=6 !<An optimiser solver \see SolverRoutines_SolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_CELLML_EVALUATOR_TYPE=7 !<A CellML evaluation solver \see SolverRoutines_SolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_STATE_ITERATION_TYPE=8 !<An state iteration solver \see SolverRoutines_SolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_GEOMETRIC_TRANSFORMATION_TYPE=9 !<An geometric transformation solver \see SolverRoutines_SolverTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_SolverLibraries SolverRoutines::SolverLibraries
  !> \brief The types of solver libraries
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_CMISS_LIBRARY=LIBRARY_CMISS_TYPE !<CMISS (internal) solver library \see SolverRoutines_SolverLibraries,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_PETSC_LIBRARY=LIBRARY_PETSC_TYPE !<PETSc solver library \see SolverRoutines_SolverLibraries,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_MUMPS_LIBRARY=LIBRARY_MUMPS_TYPE !<MUMPS solver library \see SolverRoutines_SolverLibraries,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_SUPERLU_LIBRARY=LIBRARY_SUPERLU_TYPE !<SuperLU solver library \see SolverRoutines_SolverLibraries,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_SPOOLES_LIBRARY=LIBRARY_SPOOLES_TYPE !<Spooles solver library \see SolverRoutines_SolverLibraries,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_UMFPACK_LIBRARY=LIBRARY_UMFPACK_TYPE !<UMFPACK solver library \see SolverRoutines_SolverLibraries,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_LUSOL_LIBRARY=LIBRARY_LUSOL_TYPE !<LUSOL solver library \see SolverRoutines_SolverLibraries,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_ESSL_LIBRARY=LIBRARY_ESSL_TYPE !<ESSL solver library \see SolverRoutines_SolverLibraries,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_LAPACK_LIBRARY=LIBRARY_LAPACK_TYPE !<LAPACK solver library \see SolverRoutines_SolverLibraries,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_HYPRE_LIBRARY=LIBRARY_HYPRE_TYPE !<Hypre solver library \see SolverRoutines_SolverLibraries,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_PASTIX_LIBRARY=LIBRARY_PASTIX_TYPE !<PaStiX solver library \see SolverRoutines_SolverLibraries,SolverRoutines
   !>@}

  !> \addtogroup SolverRoutines_LinearSolverTypes SolverRoutines::LinearSolverTypes
  !> \brief The types of linear solvers
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_LINEAR_DIRECT_SOLVE_TYPE=1 !<Direct linear solver type \see SolverRoutines_LinearSolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE=2 !<Iterative linear solver type \see SolverRoutines_LinearSolverTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_DirectLinearSolverTypes SolverRoutines::DirectLinearSolverTypes
  !> \brief The types of direct linear solvers
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_DIRECT_LU=1 !<LU direct linear solver \see SolverRoutines_DirectLinearSolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DIRECT_CHOLESKY=2 !<Cholesky direct linear solver \see SolverRoutines_DirectLinearSolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DIRECT_SVD=3 !<SVD direct linear solver \see SolverRoutines_DirectLinearSolverTypes,SolverRoutines
  !>@}
  
  !> \addtogroup SolverRoutines_IterativeLinearSolverTypes SolverRoutines::IterativeLinearSolverTypes
  !> \brief The types of iterative linear solvers
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_RICHARDSON=1 !<Richardson iterative solver type \see SolverRoutines_IterativeLinearSolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_CHEBYSHEV=2 !<Chebyshev iterative solver type \see SolverRoutines_IterativeLinearSolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_CONJUGATE_GRADIENT=3 !<Conjugate gradient iterative solver type \see SolverRoutines_IterativeLinearSolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_BICONJUGATE_GRADIENT=4 !<Bi-conjugate gradient iterative solver type \see SolverRoutines_IterativeLinearSolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_GMRES=5 !<Generalised minimum residual iterative solver type \see SolverRoutines_IterativeLinearSolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_BiCGSTAB=6 !<Stabalised bi-conjugate gradient iterative solver type \see SolverRoutines_IterativeLinearSolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_CONJGRAD_SQUARED=7 !<Conjugate gradient squared iterative solver type \see SolverRoutines_IterativeLinearSolverTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_IterativePreconditionerTypes SolverRoutines::IterativePreconditionerTypes
  !> \brief The types of iterative preconditioners
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_NO_PRECONDITIONER=0 !<No preconditioner type \see SolverRoutines_IterativePreconditionerTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_JACOBI_PRECONDITIONER=1 !<Jacobi preconditioner type \see SolverRoutines_IterativePreconditionerTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_BLOCK_JACOBI_PRECONDITIONER=2 !<Iterative block Jacobi preconditioner type \see SolverRoutines_IterativePreconditionerTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_SOR_PRECONDITIONER=3 !<Successive over relaxation preconditioner type \see SolverRoutines_IterativePreconditionerTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_INCOMPLETE_CHOLESKY_PRECONDITIONER=4 !<Incomplete Cholesky preconditioner type \see SolverRoutines_IterativePreconditionerTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER=5 !<Incomplete LU preconditioner type \see SolverRoutines_IterativePreconditionerTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_ADDITIVE_SCHWARZ_PRECONDITIONER=6 !<Additive Schwrz preconditioner type \see SolverRoutines_IterativePreconditionerTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_NonlinearSolverTypes SolverRoutines::NonlinearSolverTypes
  !> \brief The types of nonlinear solvers
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_NONLINEAR_NEWTON=1 !<Newton nonlinear solver type \see SolverRoutines_NonlinearSolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_NONLINEAR_BFGS_INVERSE=2 !<BFGS inverse nonlinear solver type \see SolverRoutines_NonlinearSolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_NONLINEAR_SQP=3 !<Sequential Quadratic Program nonlinear solver type \see SolverRoutines_NonlinearSolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_NONLINEAR_QUASI_NEWTON=4 !<Sequential Quasi-Newton nonlinear solver type \see SolverRoutines_NonlinearSolverTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_QuasiNewtonSolverTypes SolverRoutines::QuasiNewtonSolverTypes
  !> \brief The types of nonlinear Quasi-Newton solvers
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_LINESEARCH=1 !<Quasi-Newton line search nonlinear solver type \see SolverRoutines_QuasiNewtonSolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_TRUSTREGION=2 !<Quasi-Newton trust region nonlinear solver type \see SolverRoutines_QuasiNewtonSolverTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_QuasiNewtonTypes SolverRoutines::QuasiNewtonTypes
  !> \brief The nonlinear Quasi-Newton types
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_LBFGS=1 !<LBFGS Quasi-Newton type \see SolverRoutines_QuasiNewtonTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_GOODBROYDEN=2 !<"Good" Broyden Quasi-Newton type \see SolverRoutines_QuasiNewtonTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_BADBROYDEN=3 !<"Bad" Broyden Quasi-Newton type \see SolverRoutines_QuasiNewtonTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_QuasiNewtonLineSearchTypes SolverRoutines::QuasiNewtonLineSearchTypes
  !> \brief The types line search techniques for Quasi-Newton line search nonlinear solvers
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_LINESEARCH_BASIC=1 !<Simple damping line search. \see SolverRoutines_QuasiNewtonLineSearchTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_LINESEARCH_L2=2 !<Secant line search over the L2 norm of the function  \see SolverRoutines_QuasiNewtonLineSearchTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_LINESEARCH_CP=3 !<Critical point secant line search \see SolverRoutines_QuasiNewtonLineSearchTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_QuasiNewtonRestartTypes SolverRoutines::QuasiNewtonRestartTypes
  !> \brief The nonlinear Quasi-Newton restart types
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_RESTART_NONE=1 !<Never restart \see SolverRoutines_QuasiNewtonRestartTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_RESTART_POWELL=2 !<Restart based upon descent criteria \see SolverRoutines_QuasiNewtonRestartTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_RESTART_PERIODIC=3 !<Restart after a fixed number of iterations \see SolverRoutines_QuasiNewtonRestartTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_QuasiNewtonScaleTypes SolverRoutines::QuasiNewtonScaleTypes
  !> \brief The nonlinear Quasi-Newton scale types
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_SCALE_NONE=1 !<Don't scale the problem \see SolverRoutines_QuasiNewtonScaleTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_SCALE_SHANNO=2 !<Use Shanno scaling \see SolverRoutines_QuasiNewtonScaleTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_SCALE_LINESEARCH=3 !<Scale based upon line search lambda \see SolverRoutines_QuasiNewtonScaleTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_QUASI_NEWTON_SCALE_JACOBIAN=4 !<Scale by inverting a previously computed Jacobian \see SolverRoutines_QuasiNewtonScaleTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_NewtonSolverTypes SolverRoutines::NewtonSolverTypes
  !> \brief The types of nonlinear Newton solvers
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_LINESEARCH=1 !<Newton line search nonlinear solver type \see SolverRoutines_NewtonSolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_TRUSTREGION=2 !<Newton trust region nonlinear solver type \see SolverRoutines_NewtonSolverTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_NewtonLineSearchTypes SolverRoutines::NewtonLineSearchTypes
  !> \brief The types line search techniques for Newton line search nonlinear solvers
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_LINESEARCH_NONORMS=1 !<No norms line search for Newton line search nonlinear solves \see SolverRoutines_NewtonLineSearchTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_LINESEARCH_LINEAR=2 !<Linear search for Newton line search nonlinear solves \see SolverRoutines_NewtonLineSearchTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_LINESEARCH_QUADRATIC=3 !<Quadratic search for Newton line search nonlinear solves \see SolverRoutines_NewtonLineSearchTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_LINESEARCH_CUBIC=4!<Cubic search for Newton line search nonlinear solves \see SolverRoutines_NewtonLineSearchTypes,SolverRoutines
  !>@}
  
  !> \addtogroup SolverRoutines_JacobianCalculationTypes SolverRoutines::JacobianCalculationTypes
  !> \brief The Jacobian calculation types for a nonlinear solver 
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_JACOBIAN_NOT_CALCULATED=1 !<The Jacobian values will not be calculated for the nonlinear equations set \see SolverRoutines_JacobianCalculationTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED=2 !<The Jacobian values will be calculated analytically for the nonlinear equations set \see SolverRoutines_JacobianCalculationTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_JACOBIAN_FD_CALCULATED=3 !<The Jacobian values will be calculated using finite differences for the nonlinear equations set \see SolverRoutines_JacobianCalculationTypes,SolverRoutines
  !>@} 
  
  !> \addtogroup SolverRoutines_NewtonConvergenceTestTypes SolverRoutines::NewtonConvergenceTestTypes
  !> \brief The convergence test types for a nonlinear solver 
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_CONVERGENCE_PETSC_DEFAULT=1 !<Petsc default convergence test \see SolverRoutines_NewtonConvergenceTestTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_CONVERGENCE_ENERGY_NORM=2 !<Energy norm convergence test \see SolverRoutines_NewtonConvergenceTestTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_CONVERGENCE_DIFFERENTIATED_RATIO=3 !<Sum of differentiated ratios of unconstrained to constrained residuals convergence test \see SolverRoutines_NewtonConvergenceTestTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_DynamicOrderTypes SolverRoutines::DynamicOrderTypes
  !> \brief The order types for a dynamic solver 
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_FIRST_ORDER=1 !<Dynamic solver has first order terms \see SolverRoutines_DynamicOrderTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_SECOND_ORDER=2 !<Dynamic solver has second order terms \see SolverRoutines_DynamicOrderTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_DynamicLinearityTypes SolverRoutines::DynamicLinearityTypes
  !> \brief The time linearity types for a dynamic solver 
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_LINEAR=1 !<Dynamic solver has linear terms \see SolverRoutines_DynamicLinearityTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_NONLINEAR=2 !<Dynamic solver has nonlinear terms \see SolverRoutines_DynamicLinearityTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_DynamicDegreeTypes SolverRoutines::DynamicDegreeTypes
  !> \brief The time interpolation polynomial degree types for a dynamic solver 
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_FIRST_DEGREE=1 !<Dynamic solver uses a first degree polynomial for time interpolation \see SolverRoutines_DynamicDegreeTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_SECOND_DEGREE=2 !<Dynamic solver uses a second degree polynomial for time interpolation \see SolverRoutines_DynamicDegreeTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_THIRD_DEGREE=3 !<Dynamic solver uses a third degree polynomial for time interpolation \see SolverRoutines_DynamicDegreeTypes,SolverRoutines
  !>@}    
  
  !> \addtogroup SolverRoutines_DynamicSchemeTypes SolverRoutines::DynamicSchemeTypes
  !> \brief The types of dynamic solver scheme
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_EULER_SCHEME=1 !<Euler (explicit) dynamic solver \see SolverRoutines_DynamicSchemeTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_BACKWARD_EULER_SCHEME=2 !<Backward Euler (implicit) dynamic solver \see SolverRoutines_DynamicSchemeTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME=3 !<Crank-Nicolson dynamic solver \see SolverRoutines_DynamicSchemeTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_GALERKIN_SCHEME=4 !<Galerkin dynamic solver \see SolverRoutines_DynamicSchemeTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_ZLAMAL_SCHEME=5 !<Zlamal dynamic solver \see SolverRoutines_DynamicorderTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_SECOND_DEGREE_GEAR_SCHEME=6 !<2nd degree Gear dynamic solver \see SolverRoutines_DynamicSchemeTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_SECOND_DEGREE_LINIGER1_SCHEME=7 !<1st 2nd degree Liniger dynamic solver \see SolverRoutines_DynamicSchemeTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_SECOND_DEGREE_LINIGER2_SCHEME=8 !<2nd 2nd degree Liniger dynamic solver \see SolverRoutines_DynamicSchemeTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_NEWMARK1_SCHEME=9 !<1st Newmark dynamic solver \see SolverRoutines_DynamicSchemeTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_NEWMARK2_SCHEME=10 !<2nd Newmark dynamic solver \see SolverRoutines_DynamicSchemeTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_NEWMARK3_SCHEME=11 !<3rd Newmark dynamic solver \see SolverRoutines_DynamicSchemeTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_THIRD_DEGREE_GEAR_SCHEME=12 !<3rd degree Gear dynamic solver \see SolverRoutines_DynamicSchemeTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_THIRD_DEGREE_LINIGER1_SCHEME=13 !<1st 3rd degree Liniger dynamic solver \see SolverRoutines_DynamicSchemeTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_THIRD_DEGREE_LINIGER2_SCHEME=14 !<2nd 3rd degree Liniger dynamic solver \see SolverRoutines_DynamicSchemeTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_HOUBOLT_SCHEME=15 !<Houbolt dynamic solver \see SolverRoutines_DynamicSchemeTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_WILSON_SCHEME=16 !<Wilson dynamic solver \see SolverRoutines_DynamicSchemeTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_BOSSAK_NEWMARK1_SCHEME=17 !<1st Bossak-Newmark dynamic solver \see SolverRoutines_DynamicSchemeTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_BOSSAK_NEWMARK2_SCHEME=18 !<2nd Bossak-Newmark dynamic solver \see SolverRoutines_DynamicSchemeTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_HILBERT_HUGHES_TAYLOR1_SCHEME=19 !<1st Hilbert-Hughes-Taylor dynamic solver \see SolverRoutines_DynamicSchemeTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_HILBERT_HUGHES_TAYLOR2_SCHEME=20 !<1st Hilbert-Hughes-Taylor dynamic solver \see SolverRoutines_DynamicSchemeTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_USER_DEFINED_SCHEME=21 !<User specified degree and theta dynamic solver \see SolverRoutines_DynamicSchemeTypes,SolverRoutines
  !>@}
  
  !> \addtogroup SolverRoutines_DynamicStartupTypes SolverRoutines::DynamicStartupTypes
  !> \brief The dynamic solver previous values startup type
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_PREVIOUS_STARTUP_TYPE=1 !<Dynamic solver previous values startup type \see SolverRoutines_DynamicStartupTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_PREVIOUS_RESTART_TYPE=2 !<Dynamic solver previous values startup type type \see SolverRoutines_DynamicStartupTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_DAETypes SolverRoutines::DAETypes
  !> \brief The type of differential-algebraic equation
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_DIFFERENTIAL_ONLY=0 !<Differential equations only \see SolverRoutines_DAETypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_INDEX_1=1 !<Index 1 differential-algebraic equation \see SolverRoutines_DAETypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_INDEX_2=2 !<Index 2 differential-algebraic equation \see SolverRoutines_DAETypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_INDEX_3=3 !<Index 3 differential-algebraic equation \see SolverRoutines_DAETypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_DAESolverTypes SolverRoutines::DAESolverTypes
  !> \brief The differential-algebraic equation solver types for a differential-algebraic equation solver 
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_EULER=1 !<Euler differential-algebraic equation solver \see SolverRoutines_DAESolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_CRANK_NICOLSON=2 !<Crank-Nicolson differential-algebraic equation solver \see SolverRoutines_DAESolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_RUNGE_KUTTA=3 !<Runge-Kutta differential-algebraic equation solver \see SolverRoutines_DAESolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_ADAMS_MOULTON=4 !<Adams-Moulton differential-algebraic equation solver \see SolverRoutines_DAESolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_BDF=5 !<General BDF differential-algebraic equation solver \see SolverRoutines_DAESolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_RUSH_LARSON=6 !<Rush-Larson differential-algebraic equation solver \see SolverRoutines_DAESolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_EXTERNAL=7 !<External (e.g., CellML generated) differential-algebraic equation solver \see SolverRoutines_DAESolverTypes,SolverRoutines
  
  !>@}

  !> \addtogroup SolverRoutines_EulerDAESolverTypes SolverRoutines::EulerDAESolverTypes
  !> \brief The Euler solver types for a differential-algebriac equation solver 
  !> \see SolverRoutines_DAESolverTypes,SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_EULER_FORWARD=1 !<Forward Euler differential equation solver \see SolverRoutines_EulerDAESolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_EULER_BACKWARD=2 !<Backward Euler differential equation solver \see SolverRoutines_EulerDAESolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_DAE_EULER_IMPROVED=3 !<Improved Euler differential equation solver \see SolverRoutines_EulerDAESolverTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_OptimiserObjectiveTypes SolverRoutines::OptimiserObjectiveTypes
  !> \brief The types of objective for an optimisation problem
  !> \see SolverRoutines_OptimiserObjectiveTypes,SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_NO_OBJECTIVE=1 !<Optimisation problem has no objective (feasibility problem) \see SolverRoutines_OptimiserObjectiveTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_ONE_OBJECTIVE=2 !<Optimisation problem has one objective (standard problem) \see SolverRoutines_OptimiserObjectiveTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_MANY_OBJECTIVE=3 !<Optimisation problem has many objectives (multi-objective problem) \see SolverRoutines_OptimiserObjectiveTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_OptimiserVariableTypes SolverRoutines::OptimiserVariableTypes
  !> \brief The types of variables for an optimisation problem
  !> \see SolverRoutines_OptimiserVariablesTypes,SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_CONTINUOUS_VARIABLES=1 !<Optimisation problem has continuous variables \see SolverRoutines_OptimiserVariableTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_DISCRETE_VARIABLES=2 !<Optimisation problem has discrete variables \see SolverRoutines_OptimiserVariableTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_OptimiserCertaintyTypes SolverRoutines::OptimiserCertaintyTypes
  !> \brief The types of certainty for an optimisation problem
  !> \see SolverRoutines_OptimiserCertaintyTypes,SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_DETERMINISTIC_CERTAINTY=1 !<Optimisation problem is determanistic \see SolverRoutines_OptimiserCertaintyTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_STOCHASTIC_CERTAINTY=2 !<Optimisation problem is stochastic \see SolverRoutines_OptimiserCertaintyTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_OptimiserConstraintTypes SolverRoutines::OptimiserConstraintTypes
  !> \brief The types of constraints for an optimisation problem
  !> \see SolverRoutines_OptimiserConstraintTypes,SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_UNCONSTRAINED=1 !<Unconstrained optimisation problem \see SolverRoutines_OptimiserConstraintTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_BOUND_CONSTRAINED=2 !<Optimisation problem with bounds on the variables \see SolverRoutines_OptimiserConstraintTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_LINEAR_CONSTRAINTS=3 !<Optimisation with linear constraints \see SolverRoutines_OptimiserConstraintTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_NONLINEAR_CONSTRAINTS=4 !<Optimisation with non-linear constraints \see SolverRoutines_OptimiserConstraintTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_PDE_CONSTRAINTS=5 !<PDE constrained optimisation \see SolverRoutines_OptimiserConstraintTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_OptimiserSolverTypes SolverRoutines::OptimiserSolverTypes
  !> \brief The types of solver for an optimisation problem
  !> \see SolverRoutines_OptimiserSolverTypes,SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_NONLINEAR_LEAST_SQUARES=1 !<Non-linear least squares optimisation problem \see SolverRoutines_OptimiserSolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_LINEAR_PROGRAMMING=1 !<Linear programming optimisation problem \see SolverRoutines_OptimiserSolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_QUADRATIC_PROGRAMMING=1 !<Quadratic programming optimisation problem \see SolverRoutines_OptimiserSolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_INTEGER_PROGRAMMING=1 !<Integer programming optimisation problem \see SolverRoutines_OptimiserSolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_LINEAR_COMPLEMENTARITY=1 !<Linear complementarity optimisation problem \see SolverRoutines_OptimiserSolverTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_NONLINEAR_COMPLEMENTARITY=1 !<Nonlinear complementarity optimisation problem \see SolverRoutines_OptimiserSolverTypes,SolverRoutines
  !>@}
  
  !> \addtogroup SolverRoutines_OptimiserGradientCalculationTypes SolverRoutines::OptimiserGradientCalculationTypes
  !> \brief The gradient calculation types for an optimiser solver 
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_GRADIENT_NOT_CALCULATED=1 !<The gradient values will not be calculated for the optimiser equations set \see SolverRoutines_OptimiserGradientCalculationTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_GRADIENT_EQUATIONS_CALCULATED=2 !<The gradient values will be calculated analytically for the optimiser equations set \see SolverRoutines_OptimiserGradientCalculationTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_GRADIENT_FD_CALCULATED=3 !<The gradient values will be calculated using finite differences for the optimiser equations set \see SolverRoutines_OptimiserGradientCalculationTypes,SolverRoutines
  !>@} 
  
  !> \addtogroup SolverRoutines_OptimiserHessianCalculationTypes SolverRoutines::OptimiserHessianCalculationTypes
  !> \brief The Hessian calculation types for an optimiser solver 
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_HESSIAN_NOT_CALCULATED=1 !<The Hessian values will not be calculated for the optimiser equations set \see SolverRoutines_OptimiserHessianCalculationTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_HESSIAN_EQUATIONS_CALCULATED=2 !<The Hessian values will be calculated analytically for the optimiser equations set \see SolverRoutines_OptimiserHessianCalculationTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_OPTIMISER_HESSIAN_FD_CALCULATED=3 !<The Hessian values will be calculated using finite differences for the optimiser equations set \see SolverRoutines_OptimiserHessianCalculationTypes,SolverRoutines
  !>@} 
  
  !> \addtogroup SolverRoutines_SolutionInitialiseTypes SolverRoutines::SolutionInitialiseTypes
  !> \brief The types of solution initialisation
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_SOLUTION_INITIALISE_ZERO=0 !<Initialise the solution by zeroing it before a solve \see SolverRoutines_SolutionInitialiseTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_SOLUTION_INITIALISE_CURRENT_FIELD=1 !<Initialise the solution by copying in the current dependent field values \see SolverRoutines_SolutionInitialiseTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_SOLUTION_INITIALISE_NO_CHANGE=2 !<Do not change the solution before a solve \see SolverRoutines_SolutionInitialiseTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_OutputTypes SolverRoutines::OutputTypes
  !> \brief The types of output
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_NO_OUTPUT=0 !<No output from the solver routines \see SolverRoutines_OutputTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_MONITOR_OUTPUT=1 !<Monitor output from solver routines \see SolverRoutines_OutputTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_PROGRESS_OUTPUT=2 !<Progress output from solver routines \see SolverRoutines_OutputTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_TIMING_OUTPUT=3 !<Timing output from the solver routines plus below \see SolverRoutines_OutputTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_SOLVER_OUTPUT=4 !<Solver specific output from the solver routines plus below \see SolverRoutines_OutputTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_MATRIX_OUTPUT=5 !<Solver matrices output from the solver routines plus below\see SolverRoutines_OutputTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_EquationsSparsityTypes SolverRoutines::SparsityTypes
  !> \brief The types of sparse solver matrices
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_SPARSE_MATRICES=1 !<Use sparse solver matrices \see SolverRoutines_EqutionsSparsityTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_FULL_MATRICES=2 !<Use fully populated solver matrices \see SolverRoutines_EquationsSparsityTypes,SolverRoutines
  !>@}

  
  !> \addtogroup SolverRoutines_SymmetryTypes SolverRoutines::SymmetryTypes
  !> \brief The types of symmetry in the solver matrices
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_SYMMETRIC_MATRICES=1 !<Use symmetric solver matrices \see SolverRoutines_SymmetryTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_UNSYMMETRIC_MATRICES=2 !<Use unsymmetric solver matrices \see SolverRoutines_SymmetryTypes,SolverRoutines
  !>@}

 
  !> \addtogroup SolverRoutines_EquationsLinearityTypes SolverRoutines::EquationsLinearityTypes
  !> \brief The solver equations linearity types 
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_LINEAR=1 !<Solver equations are linear \see SolverRoutines_EquationLinearityTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_NONLINEAR=2 !<Solver equations are nonlinear \see SolverRoutines_EquationLinearityTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_EquationsTimeDependenceTypes SolverRoutines::EquationsTimeDependenceTypes
  !> \brief The solver equations time dependence types 
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_STATIC=1 !<Solver equations are static \see SolverRoutines_EquationTimeDependenceTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_QUASISTATIC=2 !<Solver equations are quasistatic \see SolverRoutines_EquationTimeDependenceTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC=3 !<Solver equations are first order dynamic \see SolverRoutines_EquationTimeDependenceTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_SECOND_ORDER_DYNAMIC=4 !<Solver equations are second order dynamic \see SolverRoutines_EquationTimeDependenceTypes,SolverRoutines
  !>@}

  
  !> \addtogroup SolverRoutines_CellMLEquationsLinearityTypes OpenCMISS::Iron::CellMLEquationsLinearityTypes
  !> \brief The CellML equations linearity types 
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: CELLML_EQUATIONS_LINEAR=1 !<CellML equations are linear \see SolverRoutines_CellMLEquationLinearityTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: CELLML_EQUATIONS_NONLINEAR=2 !<CellML equations are nonlinear \see SolverRoutines_CellMLEquationLinearityTypes,SolverRoutines
  !>@}

  !> \addtogroup SolverRoutines_CellMLEquationsTimeDependenceTypes OpenCMISS:Iron::CellMLEquationsTimeDependenceTypes
  !> \brief The CellML equations time dependence types 
  !> \see SolverRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: CELLML_EQUATIONS_STATIC=1 !<CellML equations are static \see SolverRoutines_CellMLEquationTimeDependenceTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: CELLML_EQUATIONS_QUASISTATIC=2 !<CellML equations are quasistatic \see SolverRoutines_CellMLEquationTimeDependenceTypes,SolverRoutines
  INTEGER(INTG), PARAMETER :: CELLML_EQUATIONS_DYNAMIC=3 !<CellML equations are dynamic \see SolverRoutines_CellMLEquationTimeDependenceTypes,SolverRoutines
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

  INTERFACE Solver_LabelGet
    MODULE PROCEDURE Solver_LabelGetC
    MODULE PROCEDURE Solver_LabelGetVS
  END INTERFACE Solver_LabelGet
  
  INTERFACE SolverDynamic_ThetaGet
    MODULE PROCEDURE SolverDynamic_ThetaGet0
    MODULE PROCEDURE SolverDynamic_ThetaGet1
  END INTERFACE SolverDynamic_ThetaGet
  
  PUBLIC SOLVER_NUMBER_OF_SOLVER_TYPES
  
  PUBLIC SOLVER_LINEAR_TYPE,SOLVER_NONLINEAR_TYPE,SOLVER_DYNAMIC_TYPE,SOLVER_DAE_TYPE,SOLVER_EIGENPROBLEM_TYPE, &
    & SOLVER_OPTIMISER_TYPE,SOLVER_CELLML_EVALUATOR_TYPE,SOLVER_STATE_ITERATION_TYPE,SOLVER_GEOMETRIC_TRANSFORMATION_TYPE

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
  
  PUBLIC CellMLEquations_AssertIsFinished,CellMLEquations_AssertNotFinished

  PUBLIC CellMLEquations_CellMLGet

  PUBLIC CellMLEquations_LinearityTypeGet

  PUBLIC CellMLEquations_NumberOfCellMLEnvironmentsGet

  PUBLIC CellMLEquations_SolverGet

  PUBLIC CellMLEquations_TimeDependenceTypeGet

  PUBLIC CellMLEquations_TimeGet

  PUBLIC Solver_AssertIsFinished,Solver_AssertNotFinished

  PUBLIC Solver_AssertIsDAE,Solver_AssertIsDynamic,Solver_AssertIsGeometricTransformation,Solver_AssertIsLinear, &
    & Solver_AssertIsNonlinear,Solver_AssertIsOptimiser
  
  PUBLIC Solver_CellMLEquationsExists

  PUBLIC Solver_CellMLEquationsGet

  PUBLIC Solver_ControlLoopGet

  PUBLIC Solver_DAEEulerSolverTypeGet
  
  PUBLIC Solver_DAESolverExists

  PUBLIC Solver_DAESolverGet

  PUBLIC Solver_DAESolverTypeGet

  PUBLIC Solver_DynamicDegreeGet

  PUBLIC Solver_DynamicLinearityTypeGet
  
  PUBLIC Solver_DynamicLinkedLinearSolverGet
  
  PUBLIC Solver_DynamicLinkedNonlinearSolverGet

  PUBLIC Solver_DynamicRestartGet
  
  PUBLIC Solver_DynamicSolverExists

  PUBLIC Solver_DynamicSolverGet

  PUBLIC Solver_DynamicSolverInitialisedGet

  PUBLIC Solver_GeometricTransformationSolverGet

  PUBLIC Solver_GlobalNumberGet

  PUBLIC Solver_LabelGet

  PUBLIC Solver_LibraryTypeGet
  
  PUBLIC Solver_LinearSolverExists
  
  PUBLIC Solver_LinearSolverGet
  
  PUBLIC Solver_LinkingSolverExists

  PUBLIC Solver_LinkingSolverGet

  PUBLIC Solver_NewtonLinkedCellMLSolverGet
  
  PUBLIC Solver_NewtonLinkedLinearSolverGet
  
  PUBLIC Solver_NonlinearSolverExists
  
  PUBLIC Solver_NonlinearSolverGet
  
  PUBLIC Solver_OptimiserSolverExists

  PUBLIC Solver_OptimiserSolverGet

  PUBLIC Solver_OutputTypeGet

  PUBLIC Solver_QuasiNewtonLinkedCellMLSolverGet
  
  PUBLIC Solver_QuasiNewtonLinkedLinearSolverGet

  PUBLIC Solver_SolverEquationsExists

  PUBLIC Solver_SolverEquationsGet

  PUBLIC Solver_SolverFinishedGet

  PUBLIC Solver_SolverSetupGet,Solver_SolverSetupSet

  PUBLIC Solver_SolversGet

  PUBLIC Solver_SolverTypeGet

  PUBLIC Solver_WorkGroupGet

  PUBLIC Solvers_AssertIsFinished,Solvers_AssertNotFinished

  PUBLIC Solvers_ControlLoopGet

  PUBLIC Solvers_NumberOfSolversGet

  PUBLIC Solvers_SolverGet

  PUBLIC SolverCellMLEvaluator_LibraryTypeGet

  PUBLIC SolverCellMLEvaluator_SolverGet

  PUBLIC SolverDAE_AssertIsAdamsMoulton,SolverDAE_AssertIsBDF,SolverDAE_AssertIsCrankNicolson,SolverDAE_AssertIsEuler, &
    & SolverDAE_AssertIsExternal,SolverDAE_AssertIsRungeKutta,SolverDAE_AssertIsRushLarson

  PUBLIC SolverDAE_AdamsMoultonSolverGet

  PUBLIC SolverDAE_BDFSolverGet

  PUBLIC SolverDAE_CrankNicolsonSolverGet

  PUBLIC SolverDAE_EulerSolverGet

  PUBLIC SolverDAE_ExternalSolverGet

  PUBLIC SolverDAE_LibraryTypeGet

  PUBLIC SolverDAE_RungeKuttaSolverGet

  PUBLIC SolverDAE_RushLarsonSolverGet

  PUBLIC SolverDAE_SolverGet

  PUBLIC SolverDAE_SolverTypeGet

  PUBLIC SolverDAEBDF_DAESolverGet

  PUBLIC SolverDAEEuler_BackwardEulerSolverGet
  
  PUBLIC SolverDAEEuler_ForwardEulerSolverGet

  PUBLIC SolverDAEEuler_ImprovedEulerSolverGet

  PUBLIC SolverDAEEuler_LibraryTypeGet

  PUBLIC SolverDAEEuler_DAESolverGet

  PUBLIC SolverDAEEulerBackward_EulerSolverGet
  
  PUBLIC SolverDAEEulerForward_EulerSolverGet
  
  PUBLIC SolverDAEEulerImproved_EulerSolverGet
  
  PUBLIC SolverDAEExternal_DAESolverGet

  PUBLIC SolverDynamic_DegreeGet

  PUBLIC SolverDynamic_LibraryTypeGet

  PUBLIC SolverDynamic_LinkedLinearSolverGet

  PUBLIC SolverDynamic_LinkedNonlinearSolverGet

  PUBLIC SolverDynamic_OrderGet

  PUBLIC SolverDynamic_RestartGet

  PUBLIC SolverDynamic_SolverGet

  PUBLIC SolverDynamic_SolverInitialisedGet

  PUBLIC SolverDynamic_ThetaGet

  PUBLIC SolverDynamic_TimesGet

  PUBLIC SolverEigenproblem_LibraryTypeGet

  PUBLIC SolverEigenproblem_SolverGet

  PUBLIC SolverEquations_AssertIsFinished,SolverEquations_AssertNotFinished

  PUBLIC SolverEquations_AssertIsLinear,SolverEquations_AssertIsNonlinear

  PUBLIC SolverEquations_AssertIsDynamic,SolverEquations_AssertIsStatic

  PUBLIC SolverEquations_BoundaryConditionsGet

  PUBLIC SolverEquations_LinearityTypeGet

  PUBLIC SolverEquations_SolverGet
    
  PUBLIC SolverEquations_SolverMappingExists
  
  PUBLIC SolverEquations_SolverMappingGet

  PUBLIC SolverEquations_SolverMatricesExists
  
  PUBLIC SolverEquations_SolverMatricesGet

  PUBLIC SolverEquations_SparsityTypeGet

  PUBLIC SolverEquations_SymmetryTypeGet

  PUBLIC SolverEquations_TimeDependenceTypeGet

  PUBLIC SolverGeometricTransformation_ArbitraryPathGet
  
  PUBLIC SolverGeometricTransformation_FieldGet
  
  PUBLIC SolverGeometricTransformation_FieldVariableGet
  
  PUBLIC SolverGeometricTransformation_NumberOfIncrementsGet
  
  PUBLIC SolverGeometricTransformation_SolverGet

  PUBLIC SolverLinear_AssertIsDirect,SolverLinear_AssertIsIterative

  PUBLIC SolverLinear_DirectSolverGet

  PUBLIC SolverLinear_IterativeSolverGet

  PUBLIC SolverLinear_LibraryTypeGet

  PUBLIC SolverLinear_SolverGet

  PUBLIC SolverLinear_SolverTypeGet

  PUBLIC SolverLinearDirect_LibraryTypeGet

  PUBLIC SolverLinearDirect_LinearSolverGet

  PUBLIC SolverLinearIterative_LibraryTypeGet

  PUBLIC SolverLinearIterative_LinearSolverGet

  PUBLIC SolverNonlinear_AssertIsNewton,SolverNonlinear_AssertIsQuasiNewton

  PUBLIC SolverNonlinear_LibraryTypeGet
  
  PUBLIC SolverNonlinear_LinearSolverGet
  
  PUBLIC SolverNonlinear_NewtonSolverGet
  
  PUBLIC SolverNonlinear_QuasiNewtonSolverGet
  
  PUBLIC SolverNonlinear_SolverGet

  PUBLIC SolverNonlinear_SolverTypeGet

  PUBLIC SolverNonlinearNewton_AssertIsLinesearch,SolverNonlinearNewton_AssertIsTrustregion

  PUBLIC SolverNonlinearNewton_ConvergenceTestGet

  PUBLIC SolverNonlinearNewton_LibraryTypeGet

  PUBLIC SolverNonlinearNewton_LinesearchSolverGet

  PUBLIC SolverNonlinearNewton_LinkedCellMLSolverExists

  PUBLIC SolverNonlinearNewton_LinkedCellMLSolverGet

  PUBLIC SolverNonlinearNewton_LinkedLinearSolverGet

  PUBLIC SolverNonlinearNewton_NonlinearSolverGet

  PUBLIC SolverNonlinearNewton_SolverTypeGet

  PUBLIC SolverNonlinearNewton_TrustregionSolverGet
 
  PUBLIC SolverNonlinearNewtonLinesearch_NewtonSolverGet  

  PUBLIC SolverNonlinearNewtonTrustregion_NewtonSolverGet  

  PUBLIC SolverNonlinearQuasiNewton_AssertIsLinesearch,SolverNonlinearQuasiNewton_AssertIsTrustregion

  PUBLIC SolverNonlinearQuasiNewton_ConvergenceTestGet

  PUBLIC SolverNonlinearQuasiNewton_LibraryTypeGet

  PUBLIC SolverNonlinearQuasiNewton_LinesearchSolverGet

  PUBLIC SolverNonlinearQuasiNewton_LinkedCellMLSolverExists

  PUBLIC SolverNonlinearQuasiNewton_LinkedCellMLSolverGet

  PUBLIC SolverNonlinearQuasiNewton_LinkedLinearSolverGet

  PUBLIC SolverNonlinearQuasiNewton_NonlinearSolverGet

  PUBLIC SolverNonlinearQuasiNewton_SolverTypeGet

  PUBLIC SolverNonlinearQuasiNewton_TrustregionSolverGet

  PUBLIC SolverNonlinearQuasiNewtonLinesearch_QuasiNewtonSolverGet  

  PUBLIC SolverNonlinearQuasiNewtonTrustregion_QuasiNewtonSolverGet

  PUBLIC SolverOptimiser_CertaintyTypeGet

  PUBLIC SolverOptimiser_ConstraintTypeGet

  PUBLIC SolverOptimiser_LibraryTypeGet

  PUBLIC SolverOptimiser_ObjectiveTypeGet

  PUBLIC SolverOptimiser_SolverGet

  PUBLIC SolverOptimiser_VariableTypeGet
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Assert that a CellML equations has been finished
  SUBROUTINE CellMLEquations_AssertIsFinished(cellMLEquations,err,error,*)

    !Argument Variables
    TYPE(CellMLEquationsType), POINTER, INTENT(IN) :: cellMLEquations !<The CellML equations to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("CellMLEquations_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellMLEquations)) CALL FlagError("CellML equations is not associated.",err,error,*999)
#endif    

    IF(.NOT.cellMLEquations%cellMLEquationsFinished) CALL FlagError("CellML equations has not been finished.",err,error,*999)
    
    EXITS("CellMLEquations_AssertIsFinished")
    RETURN
999 ERRORSEXITS("CellMLEquations_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLEquations_AssertIsFinished

  !
  !================================================================================================================================
  !

  !>Assert that a CellML equations has not been finished
  SUBROUTINE CellMLEquations_AssertNotFinished(cellMLEquations,err,error,*)

    !Argument Variables
    TYPE(CellMLEquationsType), POINTER, INTENT(IN) :: cellMLEquations !<The CellML equations to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("CellMLEquations_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellMLEquations)) CALL FlagError("CellML equations is not associated.",err,error,*999)
#endif    

    IF(cellMLEquations%cellMLEquationsFinished) CALL FlagError("CellML equations has already been finished.",err,error,*999)
    
    EXITS("CellMLEquations_AssertNotFinished")
    RETURN
999 ERRORSEXITS("CellMLEquations_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLEquations_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Returns a pointer to a CellML environment for a CellML equations.
  SUBROUTINE CellMLEquations_CellMLGet(cellMLEquations,cellMLIdx,cellML,err,error,*)

    TYPE(CellMLEquationsType), POINTER :: cellMLEquations !<A pointer to the CellML equations to get the CellML environment for
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: cellMLIdx !<The index of the CellML environment to get.
    TYPE(CellMLType), POINTER :: cellML !<On exit, a pointer to the CellML environments for the specified CellML equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("CellMLEquations_CellMLGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellML)) CALL FlagError("CellML environment is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellMLEquations)) CALL FlagError("CellML equations is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(cellMLEquations%cellMLEnvironments)) &
      & CALL FlagError("The CellML environments is not allocated for the CellML equations.",err,error,*999)
    IF(cellMLIdx<1.OR.cellMLIdx>cellMLEquations%numberOfCellMLEnvironments) THEN
      localError="The specified CellML index of "//TRIM(NumberToVString(cellMLIdx,"*",err,error))// &
        & " is invalid. The CellML index should be >= 1 and <= "// &
        & TRIM(NumberToVString(cellMLEquations%numberOfCellMLEnvironments,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    cellML=>cellMLEquations%cellMLEnvironments(cellMLIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(cellML)) THEN
      localError="The CellML environment is not associated for the specified CellML index of "// &
        & TRIM(NumberToVString(cellMLIdx,"*",err,error))//" of the CellML environments for the CellML equations."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("CellMLEquations_CellMLGet")
    RETURN
999 NULLIFY(cellML)
998 ERRORSEXITS("CellMLEquations_CellMLGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLEquations_CellMLGet
  
  !
  !================================================================================================================================
  !

  !>Returns the linearity type for CellML equations \see OpenCMISS::Iron::cmfe_CellMLEquations_LinearityTypeGet
  SUBROUTINE CellMLEquations_LinearityTypeGet(cellMLEquations,linearityType,err,error,*)

    !Argument variables
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations !<A pointer the CellML equations to get the linearity type for
    INTEGER(INTG), INTENT(OUT) :: linearityType !<On exit, the type of linearity of the CellML equations \see SolverRoutines_CellMLEquationLinearityTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("CellMLEquations_LinearityTypeGet",err,error,*999)

    CALL CellMLEquations_AssertIsFinished(cellMLEquations,err,error,*999)
   
    linearityType=cellMLEquations%linearity
   
    EXITS("CellMLEquations_LinearityTypeGet")
    RETURN
999 ERRORSEXITS("CellMLEquations_LinearityTypeGet",err,error)
    RETURN 1
   
  END SUBROUTINE CellMLEquations_LinearityTypeGet
        
  !
  !================================================================================================================================
  !

  !>Returns the number of CellML environments for CellML equations
  SUBROUTINE CellMLEquations_NumberOfCellMLEnvironmentsGet(cellMLEquations,numberOfCellMLEnvironments,err,error,*)

    !Argument variables
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations !<A pointer the CellML equations to get the number of CellML environments for
    INTEGER(INTG), INTENT(OUT) :: numberOfCellMLEnvironments !<On exit, the number of CellML environments for the CellML equations.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("CellMLEquations_NumberOfCellMLEnvironmentsGet",err,error,*999)

    CALL CellMLEquations_AssertIsFinished(cellMLEquations,err,error,*999)
   
    numberOfCellMLEnvironments=cellMLEquations%numberOfCellMLEnvironments
   
    EXITS("CellMLEquations_NumberOfCellMLEnvironmentsGet")
    RETURN
999 ERRORSEXITS("CellMLEquations_NumberOfCellMLEnvironmentsGet",err,error)
    RETURN 1
   
  END SUBROUTINE CellMLEquations_NumberOfCellMLEnvironmentsGet
        
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solver for a CellML equations.
  SUBROUTINE CellMLEquations_SolverGet(cellMLEquations,solver,err,error,*)

    TYPE(CellMLEquationsType), POINTER :: cellMLEquations !<A pointer to the CellML equations to get the solver for
    !Argument variables
    TYPE(SolverType), POINTER :: solver !<On exit, a pointer to the solver for the specified CellML equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("CellMLEquations_SolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solver)) CALL FlagError("Solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellMLEquations)) CALL FlagError("CellML equations is not associated.",err,error,*999)
#endif    

    solver=>cellMLEquations%solver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("CellML equations solver is not associated.",err,error,*999)
#endif    
      
    EXITS("CellMLEquations_SolverGet")
    RETURN
999 NULLIFY(solver)
998 ERRORSEXITS("CellMLEquations_SolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLEquations_SolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns the time dependence type for CellML equations \see OpenCMISS::Iron::cmfe_CellMLEquations_TimeDependenceTypeGet
  SUBROUTINE CellMLEquations_TimeDependenceTypeGet(cellMLEquations,timeDependenceType,err,error,*)

    !Argument variables
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations !<A pointer the CellML equations to get the time dependence type for
    INTEGER(INTG), INTENT(OUT) :: timeDependenceType !<On exit, the type of time dependence of the CellML equations \see SolverRoutines_CellMLEquationTimeDependenceTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("CellMLEquations_TimeDependenceTypeGet",err,error,*999)

    CALL CellMLEquations_AssertIsFinished(cellMLEquations,err,error,*999)
   
    timeDependenceType=cellMLEquations%timeDependence
   
    EXITS("CellMLEquations_TimeDependenceTypeGet")
    RETURN
999 ERRORSEXITS("CellMLEquations_TimeDependenceTypeGet",err,error)
    RETURN 1
   
  END SUBROUTINE CellMLEquations_TimeDependenceTypeGet
        
  !
  !================================================================================================================================
  !

  !>Returns the curent time for CellML equations.
  SUBROUTINE CellMLEquations_TimeGet(cellMLEquations,time,err,error,*)

    !Argument variables
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations !<A pointer the CellML equations to get the time for.
    REAL(DP), INTENT(OUT) :: time !<On exit, the time for the CellML equations
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("CellMLEquations_TimeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellMLEquations)) CALL FlagError("CellML evaluator solver is not associated.",err,error,*999)
#endif    
    
    time=cellMLEquations%currentTime
   
    EXITS("CellMLEquations_TimeGet")
    RETURN
999 ERRORSEXITS("CellMLEquations_TimeGet",err,error)
    RETURN 1
   
  END SUBROUTINE CellMLEquations_TimeGet

  !
  !================================================================================================================================
  !

  !>Assert that a solver has been finished
  SUBROUTINE Solver_AssertIsFinished(solver,err,error,*)

    !Argument Variables
    TYPE(SolverType), POINTER, INTENT(IN) :: solver !<The solver to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Solver_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    

    IF(.NOT.solver%solverFinished) CALL FlagError("Solver has not been finished.",err,error,*999)
    
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
    TYPE(SolverType), POINTER, INTENT(IN) :: solver !<The solver to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    

    IF(solver%solverFinished) CALL FlagError("Solver has already been finished.",err,error,*999)
    
    EXITS("Solver_AssertNotFinished")
    RETURN
999 ERRORSEXITS("Solver_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Assert that a solver is a DAE solver
  SUBROUTINE Solver_AssertIsDAE(solver,err,error,*)

    !Argument Variables
    TYPE(SolverType), POINTER, INTENT(INOUT) :: solver !<The solver to assert the DAE solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Solver_AssertIsDAE",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    

    IF(solver%solveType/=SOLVER_DAE_TYPE) THEN
      localError="The solver type of "//TRIM(NumberToVString(solver%solveType,"*",err,error))// &
        & " does not correspond to the required differential-algebraic equations solve type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Solver_AssertIsDAE")
    RETURN
999 ERRORSEXITS("Solver_AssertIsDAE",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_AssertIsDAE

  !
  !================================================================================================================================
  !

  !>Assert that a solver is a dynamic solver
  SUBROUTINE Solver_AssertIsDynamic(solver,err,error,*)

    !Argument Variables
    TYPE(SolverType), POINTER, INTENT(INOUT) :: solver !<The solver to assert the dynamic solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Solver_AssertIsDynamic",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    

    IF(solver%solveType/=SOLVER_DYNAMIC_TYPE) THEN
      localError="The solver type of "//TRIM(NumberToVString(solver%solveType,"*",err,error))// &
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

  !>Assert that a solver is a geometric transformation solver
  SUBROUTINE Solver_AssertIsGeometricTransformation(solver,err,error,*)

    !Argument Variables
    TYPE(SolverType), POINTER, INTENT(INOUT) :: solver !<The solver to assert the geometric transformation solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Solver_AssertIsGeometricTransformation",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    

    IF(solver%solveType/=SOLVER_GEOMETRIC_TRANSFORMATION_TYPE) THEN
      localError="The solver type of "//TRIM(NumberToVString(solver%solveType,"*",err,error))// &
        & " does not correspond to the required geometric transformation solve type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Solver_AssertIsGeometricTransformation")
    RETURN
999 ERRORSEXITS("Solver_AssertIsGeometricTransformation",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_AssertIsGeometricTransformation

  !
  !================================================================================================================================
  !

  !>Assert that a solver is a linear solver
  SUBROUTINE Solver_AssertIsLinear(solver,err,error,*)

    !Argument Variables
    TYPE(SolverType), POINTER, INTENT(INOUT) :: solver !<The solver to assert the linear solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Solver_AssertIsLinear",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    

    IF(solver%solveType/=SOLVER_LINEAR_TYPE) THEN
      localError="The solver type of "//TRIM(NumberToVString(solver%solveType,"*",err,error))// &
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
    TYPE(SolverType), POINTER, INTENT(INOUT) :: solver !<The solver to assert the nonlinear solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Solver_AssertIsNonlinear",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    

    IF(solver%solveType/=SOLVER_NONLINEAR_TYPE) THEN
      localError="The solver type of "//TRIM(NumberToVString(solver%solveType,"*",err,error))// &
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
    TYPE(SolverType), POINTER, INTENT(INOUT) :: solver !<The solver to assert the optimiser solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Solver_AssertIsOptimiser",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    

    IF(solver%solveType/=SOLVER_OPTIMISER_TYPE) THEN
      localError="The solver type of "//TRIM(NumberToVString(solver%solveType,"*",err,error))// &
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

  !>Checks if the CellML equations for a solver exist.
  SUBROUTINE Solver_CellMLEquationsExists(solver,cellMLEquations,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to check the CellML equations for
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations !<On exit, a pointer to the specified CellML equations if they exist, null if they don;t. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_CellMLEquationsExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLEquations)) CALL FlagError("CellML equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif
    
    cellMLEquations=>solver%cellMLEquations
      
    EXITS("Solver_CellMLEquationsExists")
    RETURN
999 NULLIFY(cellMLEquations)
998 ERRORSEXITS("Solver_CellMLEquationsExists",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_CellMLEquationsExists
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the CellML equations for a solver. \see OpenCMISS::Iron::cmfe_Solver_CellMLEquationsGet
  SUBROUTINE Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to get the CellML equations for
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations !<On exit, a pointer to the specified CellML equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_CellMLEquationsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLEquations)) CALL FlagError("CellML equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif
    
    cellMLEquations=>solver%cellMLEquations

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(cellMLEquations)) CALL FlagError("Solver CellML equations is not associated.",err,error,*999)
#endif    
      
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
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to get the control loop for
    TYPE(ControlLoopType), POINTER :: controlLoop !<On exit, a pointer to the control loop for the specified solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(SolversType), POINTER :: solvers
 
    ENTERS("Solver_ControlLoopGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(controlLoop)) CALL FlagError("Control loop is already associated.",err,error,*998)
#endif 
    CALL Solver_AssertIsFinished(solver,err,error,*999)

    solvers=>solver%solvers
    IF(.NOT.ASSOCIATED(solvers)) THEN
#ifdef WITH_PRECHECKS      
      IF(.NOT.ASSOCIATED(solver%linkingSolver)) CALL FlagError("Solver solvers is not associated.",err,error,*999)
#endif      
      solvers=>solver%linkingSolver%solvers
#ifdef WITH_PRECHECKS
      IF(.NOT.ASSOCIATED(solvers)) CALL FlagError("Solver linking solver solvers is not associated.",err,error,*999)
#endif      
    ENDIF
    controlLoop=>solvers%controlLoop
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Solvers control loop is not associated.",err,error,*999)
#endif    
    
    EXITS("Solver_ControlLoopGet")
    RETURN
999 NULLIFY(controlLoop)
998 ERRORSEXITS("Solver_ControlLoopGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_ControlLoopGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the solve type for an Euler differential-algebraic equation solver.
  SUBROUTINE Solver_DAEEulerSolverTypeGet(solver,daeEulerType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the Euler differential equation solver to get type for 
    INTEGER(INTG), INTENT(OUT) :: daeEulerType !<On return, the type of Euler solver for the Euler differential-algebraic equation to set \see SolverRoutines_EulerDAESolverTypes,SolverRoutines.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DAESolverType), POINTER :: daeSolver
    TYPE(EulerDAESolverType), POINTER :: eulerSolver
     
    ENTERS("Solver_DAEEulerSolverTypeGet",err,error,*999)

    CALL Solver_AssertIsFinished(solver,err,error,*999)
    CALL Solver_AssertIsDAE(solver,err,error,*999)
    NULLIFY(daeSolver)
    CALL Solver_DAESolverGet(solver,daeSolver,err,error,*999)
    CALL SolverDAE_AssertIsEuler(daeSolver,err,error,*999)
    NULLIFY(eulerSolver)
    CALL SolverDAE_EulerSolverGet(daeSolver,eulerSolver,err,error,*999)
    daeEulerType=eulerSolver%eulerType
         
    EXITS("Solver_DAEEulerSolverTypeGet")
    RETURN
999 ERRORSEXITS("Solver_DAEEulerSolverTypeGet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_DAEEulerSolverTypeGet

  !
  !================================================================================================================================
  !

  !>Checks if the DAE solver for a solver exists.
  SUBROUTINE Solver_DAESolverExists(solver,daeSolver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to check the DAE solver for
    TYPE(DAESolverType), POINTER :: daeSolver !<On exit, a pointer to the DAE solver for the specified solver if it exists. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_DAESolverExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(daeSolver)) CALL FlagError("DAE solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)
#endif    

    daeSolver=>solver%DAESolver
    
    EXITS("Solver_DAESolverExists")
    RETURN
999 NULLIFY(daeSolver)
998 ERRORSEXITS("Solver_DAESolverExists",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_DAESolverExists
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the DAE solver for a solver.
  SUBROUTINE Solver_DAESolverGet(solver,daeSolver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to get the DAE solver for
    TYPE(DAESolverType), POINTER :: daeSolver !<On exit, a pointer to the DAE solver for the specified solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_DAESolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(daeSolver)) CALL FlagError("DAE solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)
#endif    

    daeSolver=>solver%DAESolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("Solver DAE solver is not associated.",err,error,*999)
#endif    
    
    EXITS("Solver_DAESolverGet")
    RETURN
999 NULLIFY(daeSolver)
998 ERRORSEXITS("Solver_DAESolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_DAESolverGet
  
  !
  !================================================================================================================================
  !

  !>Assert that a DAE solver is an Adams Moulton solver
  SUBROUTINE SolverDAE_AssertIsAdamsMoulton(daeSolver,err,error,*)

    !Argument Variables
    TYPE(DAESolverType), POINTER, INTENT(INOUT) :: daeSolver !<The DAE solver to assert the Adams Moulton solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables    
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverDAE_AssertIsAdamsMoulton",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("DAE solver is not associated.",err,error,*999)
#endif    

    IF(daeSolver%daeSolveType/=SOLVER_DAE_ADAMS_MOULTON) THEN
      localError="The DAE solver type of "//TRIM(NumberToVString(daeSolver%daeSolveType,"*",err,error))// &
        & " does not correspond to the required Adams-Moulton differential-algebraic equations solve type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("SolverDAE_AssertIsAdamsMoulton")
    RETURN
999 ERRORSEXITS("SolverDAE_AssertIsAdamsMoulton",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_AssertIsAdamsMoulton

  !
  !================================================================================================================================
  !

  !>Assert that a DAE solver is a BDF solver
  SUBROUTINE SolverDAE_AssertIsBDF(daeSolver,err,error,*)

    !Argument Variables
    TYPE(DAESolverType), POINTER, INTENT(INOUT) :: daeSolver !<The DAE solver to assert the BDF solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverDAE_AssertIsBDF",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("DAE solver is not associated.",err,error,*999)
#endif    

    IF(daeSolver%daeSolveType/=SOLVER_DAE_BDF) THEN
      localError="The DAE solver type of "//TRIM(NumberToVString(daeSolver%daeSolveType,"*",err,error))// &
        & " does not correspond to the required BDF differential-algebraic equations solve type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("SolverDAE_AssertIsBDF")
    RETURN
999 ERRORSEXITS("SolverDAE_AssertIsBDF",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_AssertIsBDF

  !
  !================================================================================================================================
  !

  !>Assert that a DAE solver is a Crank-Nicolson solver
  SUBROUTINE SolverDAE_AssertIsCrankNicolson(daeSolver,err,error,*)

    !Argument Variables
    TYPE(DAESolverType), POINTER, INTENT(INOUT) :: daeSolver !<The DAE solver to assert the Crank-Nicolson solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverDAE_AssertIsCrankNicolson",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("DAE solver is not associated.",err,error,*999)
#endif    

    IF(daeSolver%daeSolveType/=SOLVER_DAE_CRANK_NICOLSON) THEN
      localError="The DAE solver type of "//TRIM(NumberToVString(daeSolver%daeSolveType,"*",err,error))// &
        & " does not correspond to the required Crank-Nicolson differential-algebraic equations solve type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("SolverDAE_AssertIsCrankNicolson")
    RETURN
999 ERRORSEXITS("SolverDAE_AssertIsCrankNicolson",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_AssertIsCrankNicolson

  !
  !================================================================================================================================
  !

  !>Assert that a DAE solver is an Euler solver
  SUBROUTINE SolverDAE_AssertIsEuler(daeSolver,err,error,*)

    !Argument Variables
    TYPE(DAESolverType), POINTER, INTENT(INOUT) :: daeSolver !<The DAE solver to assert the Euler solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverDAE_AssertIsEuler",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("DAE solver is not associated.",err,error,*999)
#endif    

    IF(daeSolver%daeSolveType/=SOLVER_DAE_EULER) THEN
      localError="The DAE solver type of "//TRIM(NumberToVString(daeSolver%daeSolveType,"*",err,error))// &
        & " does not correspond to the required Euler differential-algebraic equations solve type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("SolverDAE_AssertIsEuler")
    RETURN
999 ERRORSEXITS("SolverDAE_AssertIsEuler",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_AssertIsEuler

  !
  !================================================================================================================================
  !

  !>Assert that a DAE solver is an external solver
  SUBROUTINE SolverDAE_AssertIsExternal(daeSolver,err,error,*)

    !Argument Variables
    TYPE(DAESolverType), POINTER, INTENT(INOUT) :: daeSolver !<The DAE solver to assert the external solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverDAE_AssertIsExternal",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("DAE solver is not associated.",err,error,*999)
#endif    

    IF(daeSolver%daeSolveType/=SOLVER_DAE_EXTERNAL) THEN
      localError="The DAE solver type of "//TRIM(NumberToVString(daeSolver%daeSolveType,"*",err,error))// &
        & " does not correspond to the required external differential-algebraic equations solve type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("SolverDAE_AssertIsExternal")
    RETURN
999 ERRORSEXITS("SolverDAE_AssertIsExternal",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_AssertIsExternal

  !
  !================================================================================================================================
  !

  !>Assert that a DAE solver is an Runge-Kutta solver
  SUBROUTINE SolverDAE_AssertIsRungeKutta(daeSolver,err,error,*)

    !Argument Variables
    TYPE(DAESolverType), POINTER, INTENT(INOUT) :: daeSolver !<The DAE solver to assert the Runge-Kutta solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverDAE_AssertIsRungeKutta",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("DAE solver is not associated.",err,error,*999)
#endif    

    IF(daeSolver%daeSolveType/=SOLVER_DAE_RUNGE_KUTTA) THEN
      localError="The DAE solver type of "//TRIM(NumberToVString(daeSolver%daeSolveType,"*",err,error))// &
        & " does not correspond to the required Runge-Kutta differential-algebraic equations solve type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("SolverDAE_AssertIsRungeKutta")
    RETURN
999 ERRORSEXITS("SolverDAE_AssertIsRungeKutta",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_AssertIsRungeKutta

  !
  !================================================================================================================================
  !

  !>Assert that a DAE solver is an Rush-Larson solver
  SUBROUTINE SolverDAE_AssertIsRushLarson(daeSolver,err,error,*)

    !Argument Variables
    TYPE(DAESolverType), POINTER, INTENT(INOUT) :: daeSolver !<The DAE solver to assert the Rush-Larson solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverDAE_AssertIsRushLarson",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("DAE solver is not associated.",err,error,*999)
#endif    

    IF(daeSolver%daeSolveType/=SOLVER_DAE_RUSH_LARSON) THEN
      localError="The DAE solver type of "//TRIM(NumberToVString(daeSolver%daeSolveType,"*",err,error))// &
        & " does not correspond to the required Rush-Larson differential-algebraic equations solve type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("SolverDAE_AssertIsRushLarson")
    RETURN
999 ERRORSEXITS("SolverDAE_AssertIsRushLarson",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_AssertIsRushLarson

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the Adams-Moulton solver for a DAE solver.
  SUBROUTINE SolverDAE_AdamsMoultonSolverGet(daeSolver,adamsMoultonSolver,err,error,*)

    !Argument variables
    TYPE(DAESolverType), POINTER :: daeSolver !<A pointer to the DAE solver to get the Adams-Moulton solver for
    TYPE(AdamsMoultonDAESolverType), POINTER :: adamsMoultonSolver !<On exit, a pointer to the Adams-Moulton solver for the specified DAE solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverDAE_AdamsMoultonSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(adamsMoultonSolver)) CALL FlagError("Adams-Moulton solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("DAE solver is not associated.",err,error,*999)
#endif    

    adamsMoultonSolver=>daeSolver%adamsMoultonSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(adamsMoultonSolver)) CALL FlagError("DAE Solver Adams-Moulton solver is not associated.",err,error,*999)
#endif    
    
    EXITS("SolverDAE_AdamsMoultonSolverGet")
    RETURN
999 NULLIFY(adamsMoultonSolver)
998 ERRORSEXITS("SolverDAE_AdamsMoultonSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_AdamsMoultonSolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the BDF solver for a DAE solver.
  SUBROUTINE SolverDAE_BDFSolverGet(daeSolver,bdfSolver,err,error,*)

    !Argument variables
    TYPE(DAESolverType), POINTER :: daeSolver !<A pointer to the DAE solver to get the BDF solver for
    TYPE(BDFDAESolverType), POINTER :: bdfSolver !<On exit, a pointer to the BDF solver for the specified DAE solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverDAE_BDFSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(bdfSolver)) CALL FlagError("BDF solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("DAE solver is not associated.",err,error,*999)
#endif    

    bdfSolver=>daeSolver%bdfSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(bdfSolver)) CALL FlagError("DAE Solver BDF solver is not associated.",err,error,*999)
#endif    
    
    EXITS("SolverDAE_BDFSolverGet")
    RETURN
999 NULLIFY(bdfSolver)
998 ERRORSEXITS("SolverDAE_BDFSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_BDFSolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the Crank-Nicolson solver for a DAE solver.
  SUBROUTINE SolverDAE_CrankNicolsonSolverGet(daeSolver,crankNicolsonSolver,err,error,*)

    !Argument variables
    TYPE(DAESolverType), POINTER :: daeSolver !<A pointer to the DAE solver to get the Crank-Nicolson solver for
    TYPE(CrankNicolsonDAESolverType), POINTER :: crankNicolsonSolver !<On exit, a pointer to the Crank-Nicolson solver for the specified DAE solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverDAE_CrankNicolsonSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(crankNicolsonSolver)) CALL FlagError("Crank-Nicolson solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("DAE solver is not associated.",err,error,*999)
#endif    

    crankNicolsonSolver=>daeSolver%crankNicolsonSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(crankNicolsonSolver)) CALL FlagError("DAE Solver Crank-Nicolson solver is not associated.",err,error,*999)
#endif    
    
    EXITS("SolverDAE_CrankNicolsonSolverGet")
    RETURN
999 NULLIFY(crankNicolsonSolver)
998 ERRORSEXITS("SolverDAE_CrankNicolsonSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_CrankNicolsonSolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the Euler solver for a DAE solver.
  SUBROUTINE SolverDAE_EulerSolverGet(daeSolver,eulerSolver,err,error,*)

    !Argument variables
    TYPE(DAESolverType), POINTER :: daeSolver !<A pointer to the DAE solver to get the Euler solver for
    TYPE(EulerDAESolverType), POINTER :: eulerSolver !<On exit, a pointer to the Euler solver for the specified DAE solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverDAE_EulerSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(eulerSolver)) CALL FlagError("Euler solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("DAE solver is not associated.",err,error,*999)
#endif    

    eulerSolver=>daeSolver%eulerSolver

#ifdef WTIH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(eulerSolver)) CALL FlagError("DAE Solver Euler solver is not associated.",err,error,*999)
#endif    
    
    EXITS("SolverDAE_EulerSolverGet")
    RETURN
999 NULLIFY(eulerSolver)
998 ERRORSEXITS("SolverDAE_EulerSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_EulerSolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the external solver for a DAE solver.
  SUBROUTINE SolverDAE_ExternalSolverGet(daeSolver,externalSolver,err,error,*)

    !Argument variables
    TYPE(DAESolverType), POINTER :: daeSolver !<A pointer to the DAE solver to get the Euler solver for
    TYPE(ExternalDAESolverType), POINTER :: externalSolver !<On exit, a pointer to the external solver for the specified DAE solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverDAE_ExternalSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(externalSolver)) CALL FlagError("External solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("DAE solver is not associated.",err,error,*999)
#endif    

    externalSolver=>daeSolver%externalSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(externalSolver)) CALL FlagError("DAE Solver external solver is not associated.",err,error,*999)
#endif    
    
    EXITS("SolverDAE_ExternalSolverGet")
    RETURN
999 NULLIFY(externalSolver)
998 ERRORSEXITS("SolverDAE_ExternalSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_ExternalSolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns the type of library to use for a differential-algebraic equation solver
  SUBROUTINE SolverDAE_LibraryTypeGet(daeSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(DAESolverType), POINTER :: daeSolver !<A pointer the differential-algebraic equation solver to get the library type for
    INTEGER(INTG), INTENT(OUT) :: solverLibraryType !<On return, the type of library used for the differential-algebraic equation solver \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(AdamsMoultonDAESolverType), POINTER :: adamsMoultonSolver
    TYPE(BDFDAESolverType), POINTER :: bdfSolver
    TYPE(CrankNicolsonDAESolverType), POINTER :: crankNicolsonSolver
    TYPE(EulerDAESolverType), POINTER :: eulerSolver
    TYPE(RungeKuttaDAESolverType), POINTER :: rungeKuttaSolver
    TYPE(RushLarsonDAESolverType), POINTER :: rushLarsonSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverDAE_LibraryTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("DAE solver is not associated.",err,error,*999)
#endif    
    
    SELECT CASE(daeSolver%daeSolveType)
    CASE(SOLVER_DAE_EULER)
      NULLIFY(eulerSolver)
      CALL SolverDAE_EulerSolverGet(daeSolver,eulerSolver,err,error,*999)
      CALL SolverDAEEuler_LibraryTypeGet(eulerSolver,solverLibraryType,err,error,*999)
    CASE(SOLVER_DAE_CRANK_NICOLSON)
      NULLIFY(crankNicolsonSolver)
      CALL SolverDAE_CrankNicolsonSolverGet(daeSolver,crankNicolsonSolver,err,error,*999)
      solverLibraryType=crankNicolsonSolver%solverLibrary
    CASE(SOLVER_DAE_RUNGE_KUTTA)
      NULLIFY(rungeKuttaSolver)
      CALL SolverDAE_RungeKuttaSolverGet(daeSolver,rungeKuttaSolver,err,error,*999)
      solverLibraryType=rungeKuttaSolver%solverLibrary
    CASE(SOLVER_DAE_ADAMS_MOULTON)
      NULLIFY(adamsMoultonSolver)
      CALL SolverDAE_AdamsMoultonSolverGet(daeSolver,adamsMoultonSolver,err,error,*999)
      solverLibraryType=adamsMoultonSolver%solverLibrary
    CASE(SOLVER_DAE_BDF)
      NULLIFY(bdfSolver)
      CALL SolverDAE_BDFSolverGet(daeSolver,bdfSolver,err,error,*999)
      solverLibraryType=bdfSolver%solverLibrary
    CASE(SOLVER_DAE_RUSH_LARSON)
      NULLIFY(rushLarsonSolver)
      CALL SolverDAE_RushLarsonSolverGet(daeSolver,rushLarsonSolver,err,error,*999)
      solverLibraryType=rushLarsonSolver%solverLibrary
    CASE(SOLVER_DAE_EXTERNAL)
      CALL FlagError("Can not get the solver library for an external differntial-algebraic equations solver.",err,error,*999)
    CASE DEFAULT
      localError="The differential-algebraic equations solver type of "// &
        & TRIM(NumberToVString(daeSolver%daeSolveType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("SolverDAE_LibraryTypeGet")
    RETURN
999 ERRORSEXITS("SolverDAE_LibraryTypeGet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAE_LibraryTypeGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the Runge-Kutta solver for a DAE solver.
  SUBROUTINE SolverDAE_RungeKuttaSolverGet(daeSolver,rungeKuttaSolver,err,error,*)

    !Argument variables
    TYPE(DAESolverType), POINTER :: daeSolver !<A pointer to the DAE solver to get the Runge-Kutta solver for
    TYPE(RungeKuttaDAESolverType), POINTER :: rungeKuttaSolver !<On exit, a pointer to the Runge-Kutta solver for the specified DAE solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverDAE_RungeKuttaSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(rungeKuttaSolver)) CALL FlagError("Runge-Kutta solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("DAE solver is not associated.",err,error,*999)
#endif    

    rungeKuttaSolver=>daeSolver%rungeKuttaSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(rungeKuttaSolver)) CALL FlagError("DAE Solver Runge-Kutta solver is not associated.",err,error,*999)
#endif    
    
    EXITS("SolverDAE_RungeKuttaSolverGet")
    RETURN
999 NULLIFY(rungeKuttaSolver)
998 ERRORSEXITS("SolverDAE_RungeKuttaSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_RungeKuttaSolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the Rush-Larson solver for a DAE solver.
  SUBROUTINE SolverDAE_RushLarsonSolverGet(daeSolver,rushLarsonSolver,err,error,*)

    !Argument variables
    TYPE(DAESolverType), POINTER :: daeSolver !<A pointer to the DAE solver to get the Rush-Larson solver for
    TYPE(RushLarsonDAESolverType), POINTER :: rushLarsonSolver !<On exit, a pointer to the Rush-Larson solver for the specified DAE solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverDAE_RushLarsonSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(rushLarsonSolver)) CALL FlagError("Rush-Larson solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("DAE solver is not associated.",err,error,*999)
#endif    

    rushLarsonSolver=>daeSolver%rushLarsonSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(rushLarsonSolver)) CALL FlagError("DAE Solver Rush-Larson solver is not associated.",err,error,*999)
#endif    
    
    EXITS("SolverDAE_RushLarsonSolverGet")
    RETURN
999 NULLIFY(rushLarsonSolver)
998 ERRORSEXITS("SolverDAE_RushLarsonSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_RushLarsonSolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solver for a DAE solver.
  SUBROUTINE SolverDAE_SolverGet(daeSolver,solver,err,error,*)

    !Argument variables
    TYPE(DAESolverType), POINTER :: daeSolver !<A pointer to the DAE solver to get the solver for
    TYPE(SolverType), POINTER :: solver !<On exit, a pointer to the solver for the specified DAE solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverDAE_SolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solver)) CALL FlagError("Solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("DAE solver is not associated.",err,error,*999)
#endif    

    solver=>daeSolver%solver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("DAE Solver solver is not associated.",err,error,*999)
#endif    
    
    EXITS("SolverDAE_SolverGet")
    RETURN
999 NULLIFY(solver)
998 ERRORSEXITS("SolverDAE_SolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_SolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns the solve type of a DAE solver.
  SUBROUTINE SolverDAE_SolverTypeGet(daeSolver,daeSolveType,err,error,*)

    !Argument variables
    TYPE(DAESolverType), POINTER :: daeSolver !<A pointer to the DAE solver to get the solve type for
    INTEGER(INTG), INTENT(OUT) :: daeSolveType !<On return, the DAE solver type. \see SolverRoutines_DAESolverTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverDAE_SolverTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("DAE solver is not associated.",err,error,*999)
#endif    

    daeSolveType=daeSolver%daeSolveType
    
    EXITS("SolverDAE_SolverTypeGet")
    RETURN
999 ERRORSEXITS("SolverDAE_SolverTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_SolverTypeGet

  !
  !================================================================================================================================
  !
  
  !>Returns the solve type for an differential-algebraic equation solver.
  SUBROUTINE Solver_DAESolverTypeGet(solver,daeSolveType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to get the differential-algebraic equation solver type for 
    INTEGER(INTG), INTENT(OUT) :: daeSolveType !<On return, the type of solver for the differential-algebraic equation to set \see SolverRoutines_DAESolverTypes,SolverRoutines.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DAESolverType), POINTER :: daeSolver
     
    ENTERS("Solver_DAESolverTypeGet",err,error,*999)

    CALL Solver_AssertIsFinished(solver,err,error,*999)
    CALL Solver_AssertIsDAE(solver,err,error,*999)
    NULLIFY(daeSolver)
    CALL Solver_DAESolverGet(solver,daeSolver,err,error,*999)
    
    daeSolveType=daeSolver%daeSolveType
         
    EXITS("Solver_DAESolverTypeGet")
    RETURN
999 ERRORSEXITS("Solver_DAESolverTypeGet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_DAESolverTypeGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the DAE solver for an BDF DAE solver.
  SUBROUTINE SolverDAEBDF_DAESolverGet(bdfSolver,daeSolver,err,error,*)

    !Argument variables
    TYPE(BDFDAESolverType), POINTER :: bdfSolver !<A pointer to the BDF solver to get the DAE solver for
    TYPE(DAESolverType), POINTER :: daeSolver !<On exit, a pointer to the DAE solver for the specified BDF solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverDAEBDF_DAESolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(daeSolver)) CALL FlagError("DAE solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(bdfSolver)) CALL FlagError("BDF solver is not associated.",err,error,*999)
#endif    

    daeSolver=>bdfSolver%daeSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("BDF Solver DAE solver is not associated.",err,error,*999)
#endif    
    
    EXITS("SolverDAEBDF_DAESolverGet")
    RETURN
999 NULLIFY(daeSolver)
998 ERRORSEXITS("SolverDAEBDF_DAESolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAEBDF_DAESolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the backward Euler solver for an Euler DAE solver.
  SUBROUTINE SolverDAEEuler_BackwardEulerSolverGet(eulerSolver,backwardEulerSolver,err,error,*)

    !Argument variables
    TYPE(EulerDAESolverType), POINTER :: eulerSolver !<A pointer to the Euler solver to get the backward Euler solver for
    TYPE(BackwardEulerDAESolverType), POINTER :: backwardEulerSolver !<On exit, a pointer to the backward Euler solver for the specified Euler solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverDAEEuler_BackwardEulerSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(backwardEulerSolver)) CALL FlagError("Backward Euler solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(eulerSolver)) CALL FlagError("Euler solver is not associated.",err,error,*999)
#endif    

    backwardEulerSolver=>eulerSolver%backwardEulerSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(backwardEulerSolver)) CALL FlagError("Euler Solver backward Euler solver is not associated.",err,error,*999)
#endif    
    
    EXITS("SolverDAEEuler_BackwardEulerSolverGet")
    RETURN
999 NULLIFY(backwardEulerSolver)
998 ERRORSEXITS("SolverDAEEuler_BackwardEulerSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAEEuler_BackwardEulerSolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the forward Euler solver for an Euler DAE solver.
  SUBROUTINE SolverDAEEuler_ForwardEulerSolverGet(eulerSolver,forwardEulerSolver,err,error,*)

    !Argument variables
    TYPE(EulerDAESolverType), POINTER :: eulerSolver !<A pointer to the Euler solver to get the forward Euler solver for
    TYPE(ForwardEulerDAESolverType), POINTER :: forwardEulerSolver !<On exit, a pointer to the forward Euler solver for the specified Euler solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverDAEEuler_ForwardEulerSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(forwardEulerSolver)) CALL FlagError("Forward Euler solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(eulerSolver)) CALL FlagError("Euler solver is not associated.",err,error,*999)
#endif    

    forwardEulerSolver=>eulerSolver%forwardEulerSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(forwardEulerSolver)) CALL FlagError("Euler Solver forward Euler solver is not associated.",err,error,*999)
#endif
    
    EXITS("SolverDAEEuler_ForwardEulerSolverGet")
    RETURN
999 NULLIFY(forwardEulerSolver)
998 ERRORSEXITS("SolverDAEEuler_ForwardEulerSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAEEuler_ForwardEulerSolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the improved Euler solver for an Euler DAE solver.
  SUBROUTINE SolverDAEEuler_ImprovedEulerSolverGet(eulerSolver,improvedEulerSolver,err,error,*)

    !Argument variables
    TYPE(EulerDAESolverType), POINTER :: eulerSolver !<A pointer to the Euler solver to get the improved Euler solver for
    TYPE(ImprovedEulerDAESolverType), POINTER :: improvedEulerSolver !<On exit, a pointer to the improved Euler solver for the specified Euler solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverDAEEuler_ImprovedEulerSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(improvedEulerSolver)) CALL FlagError("Improved Euler solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(eulerSolver)) CALL FlagError("Euler solver is not associated.",err,error,*998)
#endif    

    improvedEulerSolver=>eulerSolver%improvedEulerSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(improvedEulerSolver)) CALL FlagError("Euler Solver improved Euler solver is not associated.",err,error,*999)
#endif    
    
    EXITS("SolverDAEEuler_ImprovedEulerSolverGet")
    RETURN
999 NULLIFY(improvedEulerSolver)
998 ERRORSEXITS("SolverDAEEuler_ImprovedEulerSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAEEuler_ImprovedEulerSolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns the type of library to use for an Euler differential-algebraic equation solver
  SUBROUTINE SolverDAEEuler_LibraryTypeGet(eulerSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(EulerDAESolverType), POINTER :: eulerSolver !<A pointer the differential-algebraic equation Euler solver to get the library type for
    INTEGER(INTG), INTENT(OUT) :: solverLibraryType !<On return, the type of library used for the differential-algebraic equation Euler solver \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(BackwardEulerDAESolverType), POINTER :: backwardEulerSolver
    TYPE(ForwardEulerDAESolverType), POINTER :: forwardEulerSolver
    TYPE(ImprovedEulerDAESolverType), POINTER :: improvedEulerSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverDAEEuler_LibraryTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(eulerSolver)) CALL FlagError("Euler DAE solver is not associated.",err,error,*999)
#endif    
    
    SELECT CASE(eulerSolver%eulerType)
    CASE(SOLVER_DAE_EULER_BACKWARD)
      NULLIFY(backwardEulerSolver)
      CALL SolverDAEEuler_BackwardEulerSolverGet(eulerSolver,backwardEulerSolver,err,error,*999)
      solverLibraryType=backwardEulerSolver%solverLibrary
    CASE(SOLVER_DAE_EULER_FORWARD)
      NULLIFY(forwardEulerSolver)
      CALL SolverDAEEuler_ForwardEulerSolverGet(eulerSolver,forwardEulerSolver,err,error,*999)
      solverLibraryType=forwardEulerSolver%solverLibrary
    CASE(SOLVER_DAE_EULER_IMPROVED)
      NULLIFY(improvedEulerSolver)
      CALL SolverDAEEuler_ImprovedEulerSolverGet(eulerSolver,improvedEulerSolver,err,error,*999)
      solverLibraryType=improvedEulerSolver%solverLibrary
    CASE DEFAULT
      localError="The Euler differential-algebraic equations solver type of "// &
        & TRIM(NumberToVString(eulerSolver%eulerType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("SolverDAEEuler_LibraryTypeGet")
    RETURN
999 ERRORSEXITS("SolverDAEEuler_LibraryTypeGet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAEEuler_LibraryTypeGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the DAE solver for an Euler DAE solver.
  SUBROUTINE SolverDAEEuler_DAESolverGet(eulerSolver,daeSolver,err,error,*)

    !Argument variables
    TYPE(EulerDAESolverType), POINTER :: eulerSolver !<A pointer to the Euler solver to get the DAE solver for
    TYPE(DAESolverType), POINTER :: daeSolver !<On exit, a pointer to the DAE solver for the specified Euler solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverDAEEuler_DAESolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(daeSolver)) CALL FlagError("DAE solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(eulerSolver)) CALL FlagError("Euler solver is not associated.",err,error,*999)
#endif    

    daeSolver=>eulerSolver%daeSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("Euler Solver DAE solver is not associated.",err,error,*999)
#endif    
    
    EXITS("SolverDAEEuler_DAESolverGet")
    RETURN
999 NULLIFY(daeSolver)
998 ERRORSEXITS("SolverDAEEuler_DAESolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAEEuler_DAESolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the Euler solver for a backward Euler DAE solver.
  SUBROUTINE SolverDAEEulerBackward_EulerSolverGet(backwardEulerSolver,eulerSolver,err,error,*)

    !Argument variables
    TYPE(BackwardEulerDAESolverType), POINTER :: backwardEulerSolver !<A pointer to the backward Euler solver to get the Euler solver for
    TYPE(EulerDAESolverType), POINTER :: eulerSolver !<On exit, a pointer to the Euler solver for the specified backward Euler solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverDAEEulerBackward_EulerSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(eulerSolver)) CALL FlagError("Euler solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(backwardEulerSolver)) CALL FlagError("Backward Euler solver is not associated.",err,error,*999)
#endif    

    eulerSolver=>backwardEulerSolver%eulerDAESolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(eulerSolver)) CALL FlagError("Backward Euler Solver Euler solver is not associated.",err,error,*999)
#endif    
    
    EXITS("SolverDAEEulerBackward_EulerSolverGet")
    RETURN
999 NULLIFY(eulerSolver)
998 ERRORSEXITS("SolverDAEEulerBackward_EulerSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAEEulerBackward_EulerSolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the Euler solver for a forward Euler DAE solver.
  SUBROUTINE SolverDAEEulerForward_EulerSolverGet(forwardEulerSolver,eulerSolver,err,error,*)

    !Argument variables
    TYPE(ForwardEulerDAESolverType), POINTER :: forwardEulerSolver !<A pointer to the forward Euler solver to get the Euler solver for
    TYPE(EulerDAESolverType), POINTER :: eulerSolver !<On exit, a pointer to the Euler solver for the specified forward Euler solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverDAEEulerForward_EulerSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(eulerSolver)) CALL FlagError("Euler solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(forwardEulerSolver)) CALL FlagError("Forward Euler solver is not associated.",err,error,*999)
#endif    

    eulerSolver=>forwardEulerSolver%eulerDAESolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(eulerSolver)) CALL FlagError("Forward Euler Solver Euler solver is not associated.",err,error,*999)
#endif    
    
    EXITS("SolverDAEEulerForward_EulerSolverGet")
    RETURN
999 NULLIFY(eulerSolver)
998 ERRORSEXITS("SolverDAEEulerForward_EulerSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAEEulerForward_EulerSolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the Euler solver for a improved Euler DAE solver.
  SUBROUTINE SolverDAEEulerImproved_EulerSolverGet(improvedEulerSolver,eulerSolver,err,error,*)

    !Argument variables
    TYPE(ImprovedEulerDAESolverType), POINTER :: improvedEulerSolver !<A pointer to the improved Euler solver to get the Euler solver for
    TYPE(EulerDAESolverType), POINTER :: eulerSolver !<On exit, a pointer to the Euler solver for the specified improved Euler solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverDAEEulerImproved_EulerSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(eulerSolver)) CALL FlagError("Euler solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(improvedEulerSolver)) CALL FlagError("Improved Euler solver is not associated.",err,error,*999)
#endif    

    eulerSolver=>improvedEulerSolver%eulerDAESolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(eulerSolver)) CALL FlagError("Improved Euler Solver Euler solver is not associated.",err,error,*999)
#endif    
    
    EXITS("SolverDAEEulerImproved_EulerSolverGet")
    RETURN
999 NULLIFY(eulerSolver)
998 ERRORSEXITS("SolverDAEEulerImproved_EulerSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAEEulerImproved_EulerSolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the DAE solver for an external DAE solver.
  SUBROUTINE SolverDAEExternal_DAESolverGet(externalSolver,daeSolver,err,error,*)

    !Argument variables
    TYPE(ExternalDAESolverType), POINTER :: externalSolver !<A pointer to the external solver to get the DAE solver for
    TYPE(DAESolverType), POINTER :: daeSolver !<On exit, a pointer to the DAE solver for the specified external solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverDAEEuler_DAESolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(daeSolver)) CALL FlagError("DAE solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(externalSolver)) CALL FlagError("External solver is not associated.",err,error,*999)
#endif    

    daeSolver=>externalSolver%daeSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("External solver DAE solver is not associated.",err,error,*999)
#endif    
    
    EXITS("SolverDAEExternal_DAESolverGet")
    RETURN
999 NULLIFY(daeSolver)
998 ERRORSEXITS("SolverDAEExternal_DAESolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAEExternal_DAESolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns the degree of the polynomial used to interpolate time for a dynamic solver.
  SUBROUTINE Solver_DynamicDegreeGet(solver,degree,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the dynamic solver to get the degree for
    INTEGER(INTG), INTENT(OUT) :: degree !<On return, the degree of the polynomial used for time interpolation in a dynamic solver \see SolverRoutines_DynamicDegreeTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    
    ENTERS("Solver_DynamicDegreeGet",err,error,*999)

    CALL Solver_AssertIsFinished(solver,err,error,*999)
    CALL Solver_AssertIsDynamic(solver,err,error,*999)
    NULLIFY(dynamicSolver)
    CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)
    
    degree=dynamicSolver%degree
    
    EXITS("Solver_DynamicDegreeGet")
    RETURN
999 ERRORSEXITS("Solver_DynamicDegreeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_DynamicDegreeGet

  !
  !================================================================================================================================
  !

  !>Returns the linearity type for the dynamic solver. \see OpenCMISS::Iron::cmfe_Solver_DynamicLinearityTypeGet
  SUBROUTINE Solver_DynamicLinearityTypeGet(solver,linearityType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to get the dynamic linearity type for 
    INTEGER(INTG), INTENT(OUT) :: linearityType !<On return, the type of linearity \see SolverRoutines_EquationLinearityTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver !<A pointer the dynamic solver to finalise
  
    ENTERS("Solver_DynamicLinearityTypeGet",err,error,*999)

    CALL Solver_AssertIsFinished(solver,err,error,*999)
    CALL Solver_AssertIsDynamic(solver,err,error,*999)
    NULLIFY(dynamicSolver)
    CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)
    
    linearityType=dynamicSolver%linearity
   
    EXITS("Solver_DynamicLinearityTypeGet")
    RETURN
999 ERRORSEXITS("Solver_DynamicLinearityTypeGet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_DynamicLinearityTypeGet

  !
  !================================================================================================================================
  !

  !>Returns the restart value for a dynamic solver.
  SUBROUTINE Solver_DynamicRestartGet(solver,restart,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the dynamic solver to get the degree for
    LOGICAL, INTENT(OUT) :: restart !<On return, the restart value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    
    ENTERS("Solver_DynamicRestartGet",err,error,*999)

    CALL Solver_AssertIsFinished(solver,err,error,*999)
    CALL Solver_AssertIsDynamic(solver,err,error,*999)
    NULLIFY(dynamicSolver)
    CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)
    
    restart=dynamicSolver%restart
    
    EXITS("Solver_DynamicRestartGet")
    RETURN
999 ERRORSEXITS("Solver_DynamicRestartGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_DynamicRestartGet

  !
  !================================================================================================================================
  !

  !>Checks if the dynamic solver for a solver exists.
  SUBROUTINE Solver_DynamicSolverExists(solver,dynamicSolver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to check the dynamic solver for
    TYPE(DynamicSolverType), POINTER :: dynamicSolver !<On exit, a pointer to the dynamic solver for the specified solver if it exists. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_DynamicSolverExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dynamicSolver)) CALL FlagError("Dynamic solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif
    
    dynamicSolver=>solver%dynamicSolver
    
    EXITS("Solver_DynamicSolverExists")
    RETURN
999 NULLIFY(dynamicSolver)
998 ERRORSEXITS("Solver_DynamicSolverExists",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_DynamicSolverExists
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the dynamic solver for a solver.
  SUBROUTINE Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to get the dynamic solver for
    TYPE(DynamicSolverType), POINTER :: dynamicSolver !<On exit, a pointer to the dynamic solver for the specified solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_DynamicSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dynamicSolver)) CALL FlagError("Dynamic solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif
    
    dynamicSolver=>solver%dynamicSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dynamicSolver)) CALL FlagError("Solver dynamic solver is not associated.",err,error,*999)
#endif    
    
    EXITS("Solver_DynamicSolverGet")
    RETURN
999 NULLIFY(dynamicSolver)
998 ERRORSEXITS("Solver_DynamicSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_DynamicSolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns the solver initialised value for a dynamic solver.
  SUBROUTINE Solver_DynamicSolverInitialisedGet(solver,solverInitialised,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the dynamic solver to get the solver initialised flag for
    LOGICAL, INTENT(OUT) :: solverInitialised !<On return, the solver initialised value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    
    ENTERS("Solver_DynamicSolverInitialisedGet",err,error,*999)

    CALL Solver_AssertIsFinished(solver,err,error,*999)
    CALL Solver_AssertIsDynamic(solver,err,error,*999)
    NULLIFY(dynamicSolver)
    CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)
    
    solverInitialised=dynamicSolver%solverInitialised
    
    EXITS("Solver_DynamicSolverInitialisedGet")
    RETURN
999 ERRORSEXITS("Solver_DynamicSolverInitialisedGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_DynamicSolverInitialisedGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the linked linear solver for a solvers dynamic solver. \see OpenCMISS::Iron::cmfe_Solver_DynamicLinearSolverGet
  SUBROUTINE Solver_DynamicLinkedLinearSolverGet(solver,linearSolver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to get the dynamic linear solver for
    TYPE(SolverType), POINTER :: linearSolver !<On exit, a pointer to the linked linear solver for the specified solvers dynamic solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
 
    ENTERS("Solver_DynamicLinkedLinearSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)
    CALL Solver_AssertIsDynamic(solver,err,error,*999)
    dynamicSolver=>solver%dynamicSolver
    IF(.NOT.ASSOCIATED(dynamicSolver)) CALL FlagError("Solver dynamic solver is not associated.",err,error,*999)
#endif    
    
    linearSolver=>solver%dynamicSolver%linearSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("The dynamic solver linked linear solver is not associated",err,error,*999)
#endif    
    
    EXITS("Solver_DynamicLinkedLinearSolverGet")
    RETURN
999 NULLIFY(linearSolver)
998 ERRORSEXITS("Solver_DynamicLinkedLinearSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_DynamicLinkedLinearSolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the linked nonlinear solver for a solvers dynamic solver. \see OpenCMISS::Iron::cmfe_Solver_DynamicNonlinearSolverGet
  SUBROUTINE Solver_DynamicLinkedNonlinearSolverGet(solver,nonlinearSolver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to get the dynamic linked nonlinear solver for
    TYPE(SolverType), POINTER :: nonlinearSolver !<On exit, a pointer to the linked nonlinear solver for the specified solvers dynamic solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
 
    ENTERS("Solver_DynamicLinkedNonlinearSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    CALL Solver_AssertIsDynamic(solver,err,error,*999)
    dynamicSolver=>solver%dynamicSolver
    IF(.NOT.ASSOCIATED(dynamicSolver)) CALL FlagError("Solver dynamic solver is not associated.",err,error,*999)
#endif    
    
    nonlinearSolver=>solver%dynamicSolver%nonlinearSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(nonlinearSolver)) CALL FlagError("The dynamic solver nonlinear solver is not associated",err,error,*999)
#endif    
    
    EXITS("Solver_DynamicLinkedNonlinearSolverGet")
    RETURN
999 NULLIFY(nonlinearSolver)
998 ERRORSEXITS("Solver_DynamicLinkedNonlinearSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_DynamicLinkedNonlinearSolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the geometric transformation solver for a solver.
  SUBROUTINE Solver_GeometricTransformationSolverGet(solver,geometricTransformationSolver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to get the geometric transformation solver for
    TYPE(GeometricTransformationSolverType), POINTER :: geometricTransformationSolver !<On exit, a pointer to the geometric transformation solver for the specified solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_GeometricTransformationSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(geometricTransformationSolver)) &
      & CALL FlagError("Geometric transformation solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    

    geometricTransformationSolver=>solver%geometricTransformationSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(geometricTransformationSolver)) &
      & CALL FlagError("Solver geometric transformation solver is not associated.",err,error,*999)
#endif    
    
    EXITS("Solver_GeometricTransformationSolverGet")
    RETURN
999 NULLIFY(geometricTransformationSolver)
998 ERRORSEXITS("Solver_GeometricTransformationSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_GeometricTransformationSolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns the global number of a solver.
  SUBROUTINE Solver_GlobalNumberGet(solver,globalNumber,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to get the global number for
    INTEGER(INTG), INTENT(OUT) :: globalNumber !<On return, the solver global number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Solver_GlobalNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    

    globalNumber=solver%globalNumber
    
    EXITS("Solver_GlobalNumberGet")
    RETURN
999 ERRORSEXITS("Solver_GlobalNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_GlobalNumberGet

  !
  !================================================================================================================================
  !

  !>Returns the label of a solver. \see OpenCMISS::Iron::cmfe_Solver_LabelGet
  SUBROUTINE Solver_LabelGetC(solver,label,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to get the label for
    CHARACTER(LEN=*), INTENT(OUT) :: label !<On return, the solver label.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: cLength,vsLength

    ENTERS("Solver_LabelGetC",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    
    
    cLength=LEN(label)
    vsLength=LEN_TRIM(solver%label)
    IF(cLength>vsLength) THEN
      label=CHAR(solver%label,vsLength)
    ELSE
      label=CHAR(solver%label,cLength)
    ENDIF
    
    EXITS("Solver_LabelGetC")
    RETURN
999 ERRORSEXITS("Solver_LabelGetC",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_LabelGetC

   !
  !================================================================================================================================
  !

  !>Returns the label of a solver. \see OpenCMISS::Iron::cmfe_Solver_LabelGet
  SUBROUTINE Solver_LabelGetVS(solver,label,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to get the label for
    TYPE(VARYING_STRING), INTENT(OUT) :: label !<On return, the solver label.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Solver_LabelGetVS",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    
    
    label=VAR_STR(CHAR(solver%label))
    
    EXITS("Solver_LabelGetVS")
    RETURN
999 ERRORSEXITS("Solver_LabelGetVS",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_LabelGetVS

  !
  !================================================================================================================================
  !

  !>Gets the type of library to use for the solver \see OpenCMISS::Iron::cmfe_Solver_LibraryTypeGet
  SUBROUTINE Solver_LibraryTypeGet(solver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to get the library type of
    INTEGER(INTG), INTENT(OUT) :: solverLibraryType !<On exit, the type of library used for the solver \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Solver_LibraryTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    
    
    SELECT CASE(solver%solveType)
    CASE(SOLVER_LINEAR_TYPE)
      CALL SolverLinear_LibraryTypeGet(solver%linearSolver,solverLibraryType,err,error,*999)
    CASE(SOLVER_NONLINEAR_TYPE)
      CALL SolverNonlinear_LibraryTypeGet(solver%nonlinearSolver,solverLibraryType,err,error,*999)
    CASE(SOLVER_DYNAMIC_TYPE)
      CALL SolverDynamic_LibraryTypeGet(solver%dynamicSolver,solverLibraryType,err,error,*999)
    CASE(SOLVER_DAE_TYPE)
      CALL SolverDAE_LibraryTypeGet(solver%daeSolver,solverLibraryType,err,error,*999)
    CASE(SOLVER_EIGENPROBLEM_TYPE)
      CALL SolverEigenproblem_LibraryTypeGet(solver%eigenproblemSolver,solverLibraryType,err,error,*999)
    CASE(SOLVER_OPTIMISER_TYPE)
      CALL SolverOptimiser_LibraryTypeGet(solver%optimiserSolver,solverLibraryType,err,error,*999)
    CASE(SOLVER_CELLML_EVALUATOR_TYPE)
      CALL SolverCellMLEvaluator_LibraryTypeGet(solver%cellMLEvaluatorSolver,solverLibraryType,err,error,*999)
    CASE DEFAULT
      localError="The solver type of "//TRIM(NumberToVString(solver%solveType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("Solver_LibraryTypeGet")
    RETURN
999 ERRORSEXITS("Solver_LibraryTypeGet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_LibraryTypeGet
  
  !
  !================================================================================================================================
  !

  !>Checks if the linear solver for a solver exists.
  SUBROUTINE Solver_LinearSolverExists(solver,linearSolver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to check the linear solver for
    TYPE(LinearSolverType), POINTER :: linearSolver !<On exit, a pointer to the linear solver for the specified solver if it exists. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_LinearSolverExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    

    linearSolver=>solver%linearSolver
   
    EXITS("Solver_LinearSolverExists")
    RETURN
999 NULLIFY(linearSolver)
998 ERRORSEXITS("Solver_LinearSolverExists",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_LinearSolverExists
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the linear solver for a solver.
  SUBROUTINE Solver_LinearSolverGet(solver,linearSolver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to get the linear solver for
    TYPE(LinearSolverType), POINTER :: linearSolver !<On exit, a pointer to the linear solver for the specified solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_LinearSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    

    linearSolver=>solver%linearSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("Solver linear solver is not associated.",err,error,*999)
#endif    
    
    EXITS("Solver_LinearSolverGet")
    RETURN
999 NULLIFY(linearSolver)
998 ERRORSEXITS("Solver_LinearSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_LinearSolverGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the CellML solver associated with a Newton solver  \see OpenCMISS::Iron::cmfe_Solver_NewtonCellMLSolverGetSet
  SUBROUTINE Solver_NewtonLinkedCellMLSolverGet(solver,cellMLSolver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the Newton solver to get the linear solver for
    TYPE(SolverType), POINTER :: cellMLSolver !<On exit, a pointer the linear solver linked to the Newton solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(SolverType), POINTER :: linkedLinearSolver,linkedNonlinearSolver

    ENTERS("Solver_QuasiNewtonLinkedCellMLSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLSolver)) CALL FlagError("CellML solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    
    NULLIFY(nonlinearSolver)
    IF(solver%solveType==SOLVER_NONLINEAR_TYPE) THEN
      CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    ELSE IF(solver%solveType==SOLVER_DYNAMIC_TYPE) THEN
      NULLIFY(dynamicSolver)
      CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)
      NULLIFY(linkedNonlinearSolver)
      CALL SolverDynamic_LinkedNonlinearSolverGet(dynamicSolver,linkedNonlinearSolver,err,error,*999)
      CALL Solver_NonlinearSolverGet(linkedNonlinearSolver,nonlinearSolver,err,error,*999)
    ELSE
      CALL FlagError("The specified solver is not a nonlinear or dynamic nonlinear solver.",err,error,*999)
    ENDIF
    CALL SolverNonlinear_AssertIsNewton(nonlinearSolver,err,error,*999)
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
#endif
    
    IF(solver%solveType==SOLVER_NONLINEAR_TYPE) THEN
      cellMLSolver=>solver%nonlinearSolver%newtonSolver%cellMLEvaluatorSolver
    ELSE
      cellMLSolver=>solver%dynamicSolver%nonlinearSolver%nonlinearSolver%newtonSolver%cellMLEvaluatorSolver
    ENDIF

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(cellMLSolver)) CALL FlagError("Newton solver CellML solver is not associated.",err,error,*999)
#endif    
    
    EXITS("Solver_NewtonLinkedCellMLSolverGet")
    RETURN
999 NULLIFY(cellMLSolver)
998 ERRORSEXITS("Solver_NewtonLinkedCellMLSolverGet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_NewtonLinkedCellMLSolverGet

  !
  !================================================================================================================================
  !

  !>Returns the linear solver associated with a Newton solver \see OpenCMISS::Iron::cmfe_Solver_NewtonLinearSolverGet
  SUBROUTINE Solver_NewtonLinkedLinearSolverGet(solver,linearSolver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the Newton solver to get the linear solver for
    TYPE(SolverType), POINTER :: linearSolver !<On exit, a pointer the linear solver linked to the Newton solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver

    ENTERS("Solver_NewtonLinkedLinearSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)        
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsNewton(nonlinearSolver,err,error,*999)
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
#endif    
    
    linearSolver=>solver%nonlinearSolver%newtonSolver%linearSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("Newton solver linear solver is not associated.",err,error,*999)
#endif    
    
    EXITS("Solver_NewtonLinkedLinearSolverGet")
    RETURN
999 NULLIFY(linearSolver)
998 ERRORSEXITS("Solver_NewtonLinkedLinearSolverGet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_NewtonLinkedLinearSolverGet

  !
  !================================================================================================================================
  !

  !>Checks if the nonlinear solver for a solver exists.
  SUBROUTINE Solver_NonlinearSolverExists(solver,nonlinearSolver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to checks the nonlinear solver for
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver !<On exit, a pointer to the nonlinear solver for the specified solver if it exists. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_NonlinearSolverExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    

    nonlinearSolver=>solver%nonlinearSolver
    
    EXITS("Solver_NonlinearSolverExists")
    RETURN
999 NULLIFY(nonlinearSolver)
998 ERRORSEXITS("Solver_NonlinearSolverExists",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_NonlinearSolverExists
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the nonlinear solver for a solver.
  SUBROUTINE Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to get the nonlinear solver for
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver !<On exit, a pointer to the nonlinear solver for the specified solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_NonlinearSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    

    nonlinearSolver=>solver%nonlinearSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(nonlinearSolver)) CALL FlagError("Solver nonlinear solver is not associated.",err,error,*999)
#endif    
    
    EXITS("Solver_NonlinearSolverGet")
    RETURN
999 NULLIFY(nonlinearSolver)
998 ERRORSEXITS("Solver_NonlinearSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_NonlinearSolverGet
  
  !
  !================================================================================================================================
  !

  !>Checks if the optimiser solver for a solver exists.
  SUBROUTINE Solver_OptimiserSolverExists(solver,optimiserSolver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to check the optimiser solver for
    TYPE(OptimiserSolverType), POINTER :: optimiserSolver !<On exit, a pointer to the optimiser solver for the specified solver if it exists. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_OptimiserSolverExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(optimiserSolver)) CALL FlagError("Optimiser solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    

    optimiserSolver=>solver%optimiserSolver
   
    EXITS("Solver_OptimiserSolverExists")
    RETURN
999 NULLIFY(optimiserSolver)
998 ERRORSEXITS("Solver_OptimiserSolverExists",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_OptimiserSolverExists
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the optimiser solver for a solver.
  SUBROUTINE Solver_OptimiserSolverGet(solver,optimiserSolver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to get the optimiser solver for
    TYPE(OptimiserSolverType), POINTER :: optimiserSolver !<On exit, a pointer to the optimiser solver for the specified solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_OptimiserSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(optimiserSolver)) CALL FlagError("Optimiser solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    

    optimiserSolver=>solver%optimiserSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(optimiserSolver)) CALL FlagError("Solver optimiser solver is not associated.",err,error,*999)
#endif
    
    EXITS("Solver_OptimiserSolverGet")
    RETURN
999 NULLIFY(optimiserSolver)
998 ERRORSEXITS("Solver_OptimiserSolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_OptimiserSolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns the output type of a solver.
  SUBROUTINE Solver_OutputTypeGet(solver,outputType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to get the output type for
    INTEGER(INTG), INTENT(OUT) :: outputType !<On return, the solver output type. \see SolverRoutines_OutputTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Solver_OutputTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    

    outputType=solver%outputType
    
    EXITS("Solver_OutputTypeGet")
    RETURN
999 ERRORSEXITS("Solver_OutputTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_OutputTypeGet

  !
  !================================================================================================================================
  !
  
  !>Returns the CellML solver associated with a Quasi-Newton solver  \see OpenCMISS::Iron::cmfe_Solver_QuasiNewtonCellMLSolverGetSet
  SUBROUTINE Solver_QuasiNewtonLinkedCellMLSolverGet(solver,cellMLSolver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the Quasi-Newton solver to get the linear solver for
    TYPE(SolverType), POINTER :: cellMLSolver !<On exit, a pointer the linear solver linked to the Quasi-Newton solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(SolverType), POINTER :: linkedLinearSolver,linkedNonlinearSolver

    ENTERS("Solver_QuasiNewtonLinkedCellMLSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLSolver)) CALL FlagError("CellML solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    
    NULLIFY(nonlinearSolver)
    IF(solver%solveType==SOLVER_NONLINEAR_TYPE) THEN
      CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    ELSE IF(solver%solveType==SOLVER_DYNAMIC_TYPE) THEN
      NULLIFY(dynamicSolver)
      CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)
      NULLIFY(linkedNonlinearSolver)
      CALL SolverDynamic_LinkedNonlinearSolverGet(dynamicSolver,linkedNonlinearSolver,err,error,*999)
      CALL Solver_NonlinearSolverGet(linkedNonlinearSolver,nonlinearSolver,err,error,*999)
    ELSE
      CALL FlagError("The specified solver is not a nonlinear or dynamic nonlinear solver.",err,error,*999)
    ENDIF
    CALL SolverNonlinear_AssertIsQuasiNewton(nonlinearSolver,err,error,*999)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
#endif
    
    IF(solver%solveType==SOLVER_NONLINEAR_TYPE) THEN
      cellMLSolver=>solver%nonlinearSolver%quasiNewtonSolver%cellMLEvaluatorSolver
    ELSE
      cellMLSolver=>solver%dynamicSolver%nonlinearSolver%nonlinearSolver%quasiNewtonSolver%cellMLEvaluatorSolver
    ENDIF

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(cellMLSolver)) CALL FlagError("Quasi-Newton solver CellML solver is not associated.",err,error,*999)
#endif    
    
    EXITS("Solver_QuasiNewtonLinkedCellMLSolverGet")
    RETURN
999 NULLIFY(cellMLSolver)
998 ERRORSEXITS("Solver_QuasiNewtonLinkedCellMLSolverGet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_QuasiNewtonLinkedCellMLSolverGet

  !
  !================================================================================================================================
  !

  !>Returns the linear solver associated with a Quasi-Newton solver \see OpenCMISS::Iron::cmfe_Solver_QuasiNewtonLinearSolverGet
  SUBROUTINE Solver_QuasiNewtonLinkedLinearSolverGet(solver,linearSolver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the Quasi-Newton solver to get the linear solver for
    TYPE(SolverType), POINTER :: linearSolver !<On exit, a pointer the linear solver linked to the Quasi-Newton solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver

    ENTERS("Solver_QuasiNewtonLinkedLinearSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)        
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsQuasiNewton(nonlinearSolver,err,error,*999)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
#endif    
    
    linearSolver=>solver%nonlinearSolver%quasiNewtonSolver%linearSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("Quasi-Newton solver linear solver is not associated.",err,error,*999)
#endif    
    
    EXITS("Solver_QuasiNewtonLinkedLinearSolverGet")
    RETURN
999 NULLIFY(linearSolver)
998 ERRORSEXITS("Solver_QuasiNewtonLinkedLinearSolverGet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_QuasiNewtonLinkedLinearSolverGet

  !
  !================================================================================================================================
  !

  !>Checks if the solver equations for a solver exist. 
  SUBROUTINE Solver_SolverEquationsExists(solver,solverEquations,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver check the solver equations for
    TYPE(SolverEquationsType), POINTER :: solverEquations !<On exit, a pointer to the specified solver equations if they exists, null if not. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_SolverEquationsExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    

    solverEquations=>solver%solverEquations
    
    EXITS("Solver_SolverEquationsExists")
    RETURN
999 NULLIFY(solverEquations)
998 ERRORSEXITS("Solver_SolverEquationsExists",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_SolverEquationsExists
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solver equations for a solver. \see OpenCMISS::Iron::cmfe_Solver_SolverEquationsGet
  SUBROUTINE Solver_SolverEquationsGet(solver,solverEquations,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to get the solver equations for
    TYPE(SolverEquationsType), POINTER :: solverEquations !<On exit, a pointer to the specified solver equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_SolverEquationsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is already associated.",err,error,*998)
    CALL Solver_AssertIsFinished(solver,err,error,*999)
#endif    

    solverEquations=>solver%solverEquations

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver solver equations is not associated.",err,error,*999)
#endif    
      
    EXITS("Solver_SolverEquationsGet")
    RETURN
999 NULLIFY(solverEquations)
998 ERRORSEXITS("Solver_SolverEquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_SolverEquationsGet
  
  !
  !================================================================================================================================
  !

  !>Returns the finished status for a solver.
  SUBROUTINE Solver_SolverFinishedGet(solver,solverFinished,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to get the solver finished status for
    LOGICAL, INTENT(OUT) :: solverFinished !<On exit, the solver finished status.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_SolverFinishedGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    

    solverFinished=solver%solverFinished
     
    EXITS("Solver_SolverFinishedGet")
    RETURN
999 ERRORSEXITS("Solver_SolverFinishedGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_SolverFinishedGet
  
  !
  !================================================================================================================================
  !

  !>Returns the setup status for a solver.
  SUBROUTINE Solver_SolverSetupGet(solver,solverSetup,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to get the solver setup status for
    LOGICAL, INTENT(OUT) :: solverSetup !<On exit, the solver setup status.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_SolverSetupGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    

    solverSetup=solver%solverSetup
     
    EXITS("Solver_SolverSetupGet")
    RETURN
999 ERRORSEXITS("Solver_SolverSetupGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_SolverSetupGet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the setup status for a solver.
  SUBROUTINE Solver_SolverSetupSet(solver,solverSetup,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to get the solver setup status for
    LOGICAL, INTENT(IN) :: solverSetup !<The solver setup status to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_SolverSetupSet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    

    solver%solverSetup=solverSetup
     
    EXITS("Solver_SolverSetupSet")
    RETURN
999 ERRORSEXITS("Solver_SolverSetupSet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_SolverSetupSet
  
  !
  !================================================================================================================================
  !
  
  !>Checks a solvers linking solver and returns a pointer to the linking solver for a solver if it exists
  SUBROUTINE Solver_LinkingSolverExists(solver,linkingSolver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to get the solvers for.
    TYPE(SolverType), POINTER :: linkingSolver !<On exit, a pointer to the linking solver for the solver if it has one. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_LinkingSolverExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linkingSolver)) CALL FlagError("Linking solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    
      
    linkingSolver=>solver%linkingSolver
       
    EXITS("Solver_LinkingSolverExists")
    RETURN
999 NULLIFY(linkingSolver)
998 ERRORSEXITS("Solver_LinkingSolverExists",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_LinkingSolverExists

  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the linking solver for a solver.
  SUBROUTINE Solver_LinkingSolverGet(solver,linkingSolver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to get the solvers for.
    TYPE(SolverType), POINTER :: linkingSolver !<On exit, A pointer to the linking solver for the solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_LinkingSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linkingSolver)) CALL FlagError("Linking solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    
      
    linkingSolver=>solver%linkingSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(linkingSolver)) CALL FlagError("The solver linking solver is not associated.",err,error,*999)
#endif    
       
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
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to get the solvers for.
    TYPE(SolversType), POINTER :: solvers !<On exit, A pointer to the solvers for the solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(SolverType), POINTER :: testSolver
 
    ENTERS("Solver_SolversGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solvers)) CALL FlagError("Solvers is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    

    testSolver=>solver
    DO WHILE(ASSOCIATED(testSolver%linkingSolver))
      testSolver=>testSolver%linkingSolver
    ENDDO
    solvers=>testSolver%solvers

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solvers)) CALL FlagError("The solver solvers is not associated.",err,error,*999)
#endif    
       
    EXITS("Solver_SolversGet")
    RETURN
999 NULLIFY(solvers)
998 ERRORSEXITS("Solver_SolversGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_SolversGet

  !
  !================================================================================================================================
  !

  !>Returns the solve type of a solver.
  SUBROUTINE Solver_SolverTypeGet(solver,solveType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to get the solve type for
    INTEGER(INTG), INTENT(OUT) :: solveType !<On return, the solver type. \see SolverRoutines_SolverTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Solver_SolverTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    

    solveType=solver%solveType
    
    EXITS("Solver_SolverTypeGet")
    RETURN
999 ERRORSEXITS("Solver_SolverTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_SolverTypeGet

  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the work group for a solver. FOR NOW JUST RETURN THE PROBLEM WORK GROUP.
  SUBROUTINE Solver_WorkGroupGet(solver,workGroup,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to get the work group for.
    TYPE(WorkGroupType), POINTER :: workGroup !<On exit, A pointer to the work group for the solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(SolversType), POINTER :: solvers
 
    ENTERS("Solver_WorkGroupGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(workGroup)) CALL FlagError("Work group is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)

    NULLIFY(solvers)
    CALL Solver_SolversGet(solver,solvers,err,error,*999)
    IF(.NOT.ASSOCIATED(solvers)) CALL FlagError("Solver solvers is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(solvers%controlLoop)) &
      & CALL FlagError("Solver solvers control loop is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(solvers%controlLoop%problem)) &
      & CALL FlagError("Solver solvers control loop problem is not associated.",err,error,*999)
#endif    

    workGroup=>solvers%controlLoop%problem%workGroup

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("The solver work group is not associated.",err,error,*999)
#endif    
       
    EXITS("Solver_WorkGroupGet")
    RETURN
999 NULLIFY(workGroup)
998 ERRORSEXITS("Solver_WorkGroupGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_WorkGroupGet

  !
  !================================================================================================================================
  !

  !>Assert that a solvers has been finished
  SUBROUTINE Solvers_AssertIsFinished(solvers,err,error,*)

    !Argument Variables
    TYPE(SolversType), POINTER, INTENT(IN) :: solvers !<The solvers to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Solvers_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solvers)) CALL FlagError("Solvers is not associated.",err,error,*999)
#endif    

    IF(.NOT.solvers%solversFinished) CALL FlagError("Solvers has not been finished.",err,error,*999)
    
    EXITS("Solvers_AssertIsFinished")
    RETURN
999 ERRORSEXITS("Solvers_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Solvers_AssertIsFinished

  !
  !================================================================================================================================
  !

  !>Assert that a solvers has not been finished
  SUBROUTINE Solvers_AssertNotFinished(solvers,err,error,*)

    !Argument Variables
    TYPE(SolversType), POINTER, INTENT(IN) :: solvers !<The solvers to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solvers_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solvers)) CALL FlagError("Solvers is not associated.",err,error,*999)
#endif    

    IF(solvers%solversFinished) CALL FlagError("Solvers has already been finished.",err,error,*999)
    
    EXITS("Solvers_AssertNotFinished")
    RETURN
999 ERRORSEXITS("Solvers_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Solvers_AssertNotFinished

  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the control loop for a solvers.
  SUBROUTINE Solvers_ControlLoopGet(solvers,controlLoop,err,error,*)

    !Argument variables
    TYPE(SolversType), POINTER :: solvers !<A pointer to the solvers to get the control loop for.
    TYPE(ControlLoopType), POINTER :: controlLoop !<On exit, A pointer to the control loop for the solvers. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solvers_ControlLoopGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(controlLoop)) CALL FlagError("Control loop is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solvers)) CALL FlagError("Solvers is not associated.",err,error,*999)
#endif    
      
    controlLoop=>solvers%controlLoop

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("The solvers control loop is not associated.",err,error,*999)
#endif    
       
    EXITS("Solvers_ControlLoopGet")
    RETURN
999 NULLIFY(controlLoop)
998 ERRORSEXITS("Solvers_ControlLoopGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solvers_ControlLoopGet

  !
  !================================================================================================================================
  !
  
  !>Returns the number of solvers for a solvers.
  SUBROUTINE Solvers_NumberOfSolversGet(solvers,numberOfSolvers,err,error,*)

    !Argument variables
    TYPE(SolversType), POINTER :: solvers !<A pointer to the solvers to get the number of solvers for.
    INTEGER(INTG), INTENT(OUT) :: numberOfSolvers !<On exit, the number of solvers for the solvers.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solvers_NumberOfSolversGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solvers)) CALL FlagError("Solvers is not associated.",err,error,*999)
#endif    
      
    numberOfSolvers=solvers%numberOfSolvers
       
    EXITS("Solvers_NumberOfSolversGet")
    RETURN
999 ERRORSEXITS("Solvers_NumberOfSolversGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solvers_NumberOfSolversGet

  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the specified solver in the list of solvers.
  SUBROUTINE Solvers_SolverGet(solvers,solverIndex,solver,err,error,*)

    !Argument variables
    TYPE(SolversType), POINTER :: solvers !<A pointer to the solvers to get the solver for
    INTEGER(INTG), INTENT(IN) :: solverIndex !<The specified solver to get
    TYPE(SolverType), POINTER :: solver !<On exit, a pointer to the specified solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Solvers_SolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solvers)) CALL FlagError("Solvers is not associated.",err,error,*998)
    IF(ASSOCIATED(solver)) CALL FlagError("Solver is already associated.",err,error,*998)
    IF(solverIndex<=0.OR.solverIndex>solvers%numberOfSolvers) THEN
      localError="The specified solver index of "//TRIM(NumberToVString(solverIndex,"*",err,error))// &
        & " is invalid. The solver index must be >= 1 and <= "// &
        & TRIM(NumberToVString(solvers%numberOfSolvers,"*",err,error))//"."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    IF(.NOT.ALLOCATED(solvers%solvers)) CALL FlagError("Solvers solvers is not associated.",err,error,*998)
#endif    
      
    solver=>solvers%solvers(solverIndex)%ptr

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solver)) THEN
      localError="The solvers solver is not associated for solver index "// &
        & TRIM(NumberToVString(solverIndex,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
    EXITS("Solvers_SolverGet")
    RETURN
999 NULLIFY(solver)
998 ERRORSEXITS("Solvers_SolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solvers_SolverGet

  !
  !================================================================================================================================
  !

  !>Returns the type of library to use for a CellML evaluator solver.
  SUBROUTINE SolverCellMLEvaluator_LibraryTypeGet(cellMLEvaluatorSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(CellMLEvaluatorSolverType), POINTER :: cellMLEvaluatorSolver !<A pointer the CellML evaluator solver to get the library type for.
    INTEGER(INTG), INTENT(OUT) :: solverLibraryType !<On exit, the type of library used for the CellML evaluator solver \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverCellMLEvaluator_LibraryTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellMLEvaluatorSolver)) CALL FlagError("CellML evaluator solver is not associated.",err,error,*999)
#endif    
    
    solverLibraryType=cellMLEvaluatorSolver%solverLibrary
   
    EXITS("SolverCellMLEvaluator_LibraryTypeGet")
    RETURN
999 ERRORSEXITS("SolverCellMLEvaluator_LibraryTypeGet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverCellMLEvaluator_LibraryTypeGet

  !
  !================================================================================================================================
  !

  !>Returns the solver for a CellML evaluator solver.
  SUBROUTINE SolverCellMLEvaluator_SolverGet(cellMLEvaluatorSolver,solver,err,error,*)

    !Argument variables
    TYPE(CellMLEvaluatorSolverType), POINTER :: cellMLEvaluatorSolver !<A pointer the CellML evaluator solver to get the solver for.
    TYPE(SolverType), POINTER :: SOLVER !<On exit, a pointer to the solver for the CellML evaluator solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverCellMLEvaluator_SolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solver)) CALL FlagError("Solver is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(cellMLEvaluatorSolver)) CALL FlagError("CellML evaluator solver is not associated.",err,error,*999)
#endif    
    
    solver=>cellMLEvaluatorSolver%solver

#ifdef WITH_POSTCHECKS
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("The solver is not associated for the CellML evaluator solver.",err,error,*999)
#endif    
   
    EXITS("SolverCellMLEvaluator_SolverGet")
    RETURN
999 NULLIFY(solver)
998 ERRORSEXITS("SolverCellMLEvaluator_SolverGet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverCellMLEvaluator_SolverGet

  !
  !================================================================================================================================
  !

  !>Returns the degree of the polynomial used to interpolate time for a dynamic solver.
  SUBROUTINE SolverDynamic_DegreeGet(dynamicSolver,degree,err,error,*)

    !Argument variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver !<A pointer to the dynamic solver to get the degree for
    INTEGER(INTG), INTENT(OUT) :: degree !<On return, the degree of the polynomial used for time interpolation in a dynamic solver \see SolverRoutines_DynamicDegreeTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverDynamic_DegreeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dynamicSolver)) CALL FlagError("Dynamic solver is not associated.",err,error,*999)
#endif
    
    degree=dynamicSolver%degree
    
    EXITS("SolverDynamic_DegreeGet")
    RETURN
999 ERRORSEXITS("SolverDynamic_DegreeGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDynamic_DegreeGet

  !
  !================================================================================================================================
  !

  !>Returns the type of library to use for a dynamic solver.
  SUBROUTINE SolverDynamic_LibraryTypeGet(dynamicSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver !<A pointer the dynamic solver to get the library type for.
    INTEGER(INTG), INTENT(OUT) :: solverLibraryType !<On exit, the type of library used for the dynamic solver \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverDynamic_LibraryTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(dynamicSolver)) CALL FlagError("Dynamic solver is not associated.",err,error,*999)
    
    solverLibraryType=dynamicSolver%solverLibrary
    
    EXITS("SolverDynamic_LibraryTypeGet")
    RETURN
999 ERRORSEXITS("SolverDynamic_LibraryTypeGet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDynamic_LibraryTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the linear solver for a dynamic solver. 
  SUBROUTINE SolverDynamic_LinkedLinearSolverGet(dynamicSolver,linearSolver,err,error,*)

    !Argument variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver !<A pointer to the dynamic solver to get the linear solver for
    TYPE(SolverType), POINTER :: linearSolver !<On exit, a pointer to the linear solver for the specified dynamic solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverDynamic_LinkedLinearSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dynamicSolver)) CALL FlagError("Dynamic solver is not associated.",err,error,*999)
#endif    

    linearSolver=>dynamicSolver%linearSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("Dynamic solver linked linear solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverDynamic_LinkedLinearSolverGet")
    RETURN
998 NULLIFY(linearSolver)
999 ERRORSEXITS("SolverDynamic_LinkedLinearSolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverDynamic_LinkedLinearSolverGet
     
  !
  !================================================================================================================================
  !

  !>Gets the linked nonlinear solver for a dynamic solver. 
  SUBROUTINE SolverDynamic_LinkedNonlinearSolverGet(dynamicSolver,nonlinearSolver,err,error,*)

    !Argument variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver !<A pointer to the dynamic solver to get the linked nonlinear solver for
    TYPE(SolverType), POINTER :: nonlinearSolver !<On exit, a pointer to the linked nonlinear solver for the specified dynamic solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverDynamic_LinkedNonlinearSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dynamicSolver)) CALL FlagError("Dynamic solver is not associated.",err,error,*999)
#endif    

    nonlinearSolver=>dynamicSolver%nonlinearSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(nonlinearSolver)) CALL FlagError("Dynamic solver linked nonlinear solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverDynamic_LinkedNonlinearSolverGet")
    RETURN
998 NULLIFY(nonlinearSolver)
999 ERRORSEXITS("SolverDynamic_LinkedNonlinearSolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverDynamic_LinkedNonlinearSolverGet
     
  !
  !================================================================================================================================
  !

  !>Returns the order of the dynamic solver.
  SUBROUTINE SolverDynamic_OrderGet(dynamicSolver,order,err,error,*)

    !Argument variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver !<A pointer to the dynamic solver to get the order for
    INTEGER(INTG), INTENT(OUT) :: order !<On return, the order of the dynamic solver \see SolverRoutines_DynamicOrderTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverDynamic_OrderGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dynamicSolver)) CALL FlagError("Dynamic solver is not associated.",err,error,*999)
#endif
    
    order=dynamicSolver%order
    
    EXITS("SolverDynamic_OrderGet")
    RETURN
999 ERRORSEXITS("SolverDynamic_OrderGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDynamic_OrderGet

  !
  !================================================================================================================================
  !

  !>Returns the restart status of the dynamic solver.
  SUBROUTINE SolverDynamic_RestartGet(dynamicSolver,restart,err,error,*)

    !Argument variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver !<A pointer to the dynamic solver to get the restart status for
    LOGICAL, INTENT(OUT) :: restart !<On return, the restart status of the dynamic solver 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverDynamic_RestartGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dynamicSolver)) CALL FlagError("Dynamic solver is not associated.",err,error,*999)
#endif
    
    restart=dynamicSolver%restart
    
    EXITS("SolverDynamic_RestartGet")
    RETURN
999 ERRORSEXITS("SolverDynamic_RestartGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDynamic_RestartGet

  !
  !================================================================================================================================
  !

  !>Gets the solver for a dynamic solver. 
  SUBROUTINE SolverDynamic_SolverGet(dynamicSolver,solver,err,error,*)

    !Argument variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver !<A pointer to the dynamic solver to get the solver for
    TYPE(SolverType), POINTER :: solver !<On exit, a pointer to the solver for the specified dynamic solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverDynamic_SolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solver)) CALL FlagError("Solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dynamicSOlver)) CALL FlagError("Dynamic solver is not associated.",err,error,*999)
#endif    

    solver=>dynamicSolver%solver

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Dynamic solver solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverDynamic_SolverGet")
    RETURN
998 NULLIFY(solver)
999 ERRORSEXITS("SolverDynamic_SolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverDynamic_SolverGet
     
  !
  !================================================================================================================================
  !

  !>Returns the solver initialised status of the dynamic solver.
  SUBROUTINE SolverDynamic_SolverInitialisedGet(dynamicSolver,solverInitialised,err,error,*)

    !Argument variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver !<A pointer to the dynamic solver to get the solver initialised status for
    LOGICAL, INTENT(OUT) :: solverInitialised !<On return, the solver initialised status of the dynamic solver 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverDynamic_SolverInitialisedGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dynamicSolver)) CALL FlagError("Dynamic solver is not associated.",err,error,*999)
#endif
    
    solverInitialised=dynamicSolver%solverInitialised
    
    EXITS("SolverDynamic_SolverInitialisedGet")
    RETURN
999 ERRORSEXITS("SolverDynamic_SolverInitialisedGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDynamic_SolverInitialisedGet

  !
  !================================================================================================================================
  !

  !>Returns the theta value of a dynamic solver.
  SUBROUTINE SolverDynamic_ThetaGet0(dynamicSolver,theta,err,error,*)

    !Argument variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver !<A pointer to the dynamic solver to get the theta for
    REAL(DP), INTENT(OUT) :: theta !<On return, the theta value of the dynamic solver 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: thetas(1)
    
    ENTERS("SolverDynamic_ThetaGet0",err,error,*999)

    CALL SolverDynamic_ThetaGet1(dynamicSolver,thetas,err,error,*999)
    theta=thetas(1)
     
    EXITS("SolverDynamic_ThetaGet0")
    RETURN
999 ERRORSEXITS("SolverDynamic_ThetaGet0",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDynamic_ThetaGet0

  !
  !================================================================================================================================
  !

  !>Returns the theta values of a dynamic solver.
  SUBROUTINE SolverDynamic_ThetaGet1(dynamicSolver,thetas,err,error,*)

    !Argument variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver !<A pointer to the dynamic solver to get the thetas for
    REAL(DP), INTENT(OUT) :: thetas(:) !<thetas(thetaIdx). On return, the theta value of the dynamic solver 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("SolverDynamic_ThetaGet1",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(dynamicSolver)) CALL FlagError("Dynamic solver is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(dynamicSolver%theta)) CALL FlagError("Dynamic solver theta is not allocated.",err,error,*999)
    IF(SIZE(thetas,1)<SIZE(dynamicSolver%theta,1)) THEN
      localError="The size of the specified thetas array of "//TRIM(NumberToVString(SIZE(thetas,1),"*",err,error))// &
        & " is invalid. The size must be >= "//TRIM(NumberToVString(SIZE(dynamicSolver%theta,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    thetas(1:SIZE(dynamicSolver%theta,1))=dynamicSolver%theta(1:SIZE(dynamicSolver%theta,1))
     
    EXITS("SolverDynamic_ThetaGet1")
    RETURN
999 ERRORSEXITS("SolverDynamic_ThetaGet1",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDynamic_ThetaGet1

  !
  !================================================================================================================================
  !

  !>Returns the current times for a dynamic solver.
  SUBROUTINE SolverDynamic_TimesGet(dynamicSolver,currentTime,timeIncrement,err,error,*)

    !Argument variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver !<A pointer to the dynamic solver to get the solver times for
    REAL(DP), INTENT(OUT) :: currentTime !<On return, the current time of the dynamic solver 
    REAL(DP), INTENT(OUT) :: timeIncrement !<On return, the current time increment of the dynamic solver 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverDynamic_TimesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dynamicSolver)) CALL FlagError("Dynamic solver is not associated.",err,error,*999)
#endif

    currentTime=dynamicSolver%currentTime
    timeIncrement=dynamicSolver%timeIncrement
    
    EXITS("SolverDynamic_TimesGet")
    RETURN
999 ERRORSEXITS("SolverDynamic_TimesGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDynamic_TimesGet

  !
  !================================================================================================================================
  !

  !>Returns the type of library to use for an eigenproblem solver.
  SUBROUTINE SolverEigenproblem_LibraryTypeGet(eigenproblemSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(EigenproblemSolverType), POINTER :: eigenproblemSolver !<A pointer the eigenproblem solver to get the library type for.
    INTEGER(INTG), INTENT(OUT) :: solverLibraryType !<On exit, the type of library used for the eigenproblem solver \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverEigenproblem_LibraryTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(eigenproblemSolver)) CALL FlagError("Eigenproblem solver is not associated.",err,error,*999)
#endif    
    
    solverLibraryType=eigenproblemSolver%solverLibrary
     
    EXITS("SolverEigenproblem_LibraryTypeGet")
    RETURN
999 ERRORSEXITS("SolverEigenproblem_LibraryTypeGet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverEigenproblem_LibraryTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the solver for a eigenproblem solver. 
  SUBROUTINE SolverEigenproblem_SolverGet(eigenproblemSolver,solver,err,error,*)

    !Argument variables
    TYPE(EigenproblemSolverType), POINTER :: eigenproblemSolver !<A pointer to the eigenproblem solver to get the solver for
    TYPE(SolverType), POINTER :: solver !<On exit, a pointer to the solver for the specified eigenproblem solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEigenproblem_SolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solver)) CALL FlagError("Solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(eigenproblemSOlver)) CALL FlagError("Eigenproblem solver is not associated.",err,error,*999)
#endif    

    solver=>eigenproblemSolver%solver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Eigenproblem solver solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverEigenproblem_SolverGet")
    RETURN
998 NULLIFY(solver)
999 ERRORSEXITS("SolverEigenproblem_SolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverEigenproblem_SolverGet
     
  !
  !================================================================================================================================
  !

  !>Assert that a solver equations has been finished
  SUBROUTINE SolverEquations_AssertIsFinished(solverEquations,err,error,*)

    !Argument Variables
    TYPE(SolverEquationsType), POINTER, INTENT(IN) :: solverEquations !<The solver equations to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverEquations_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
#endif    

    IF(.NOT.solverEquations%solverEquationsFinished) CALL FlagError("Solver equations has not been finished.",err,error,*999)
    
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
    TYPE(SolverEquationsType), POINTER, INTENT(IN) :: solverEquations !<The solver equations to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverEquations_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
#endif    

    IF(solverEquations%solverEquationsFinished) CALL FlagError("Solver equations has already been finished.",err,error,*999)
    
    EXITS("SolverEquations_AssertNotFinished")
    RETURN
999 ERRORSEXITS("SolverEquations_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE SolverEquations_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Assert that a solver equations is linear 
  SUBROUTINE SolverEquations_AssertIsLinear(solverEquations,err,error,*)

    !Argument Variables
    TYPE(SolverEquationsType), POINTER, INTENT(INOUT) :: solverEquations !<The solver equations to assert the linearity for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverEquations_AssertIsLinear",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
#endif    

    IF(solverEquations%linearity/=SOLVER_EQUATIONS_LINEAR) THEN
      localError="The solver equations linearity type of "//TRIM(NumberToVString(solverEquations%linearity,"*",err,error))// &
        & " does not correspond to the required linear solver equations type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("SolverEquations_AssertIsLinear")
    RETURN
999 ERRORSEXITS("SolverEquations_AssertIsLinear",err,error)
    RETURN 1
    
  END SUBROUTINE SolverEquations_AssertIsLinear

  !
  !================================================================================================================================
  !

  !>Assert that a solver equations is nonlinear 
  SUBROUTINE SolverEquations_AssertIsNonlinear(solverEquations,err,error,*)

    !Argument Variables
    TYPE(SolverEquationsType), POINTER, INTENT(INOUT) :: solverEquations !<The solver equations to assert the nonlinearity for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverEquations_AssertIsNonlinear",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
#endif    

    IF(solverEquations%linearity/=SOLVER_EQUATIONS_NONLINEAR) THEN
      localError="The solver equations linearity type of "//TRIM(NumberToVString(solverEquations%linearity,"*",err,error))// &
        & " does not correspond to the required nonlinear solver equations type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("SolverEquations_AssertIsNonlinear")
    RETURN
999 ERRORSEXITS("SolverEquations_AssertIsNonlinear",err,error)
    RETURN 1
    
  END SUBROUTINE SolverEquations_AssertIsNonlinear
  
  !
  !================================================================================================================================
  !

  !>Assert that solver equations time dependence is dynamic
  SUBROUTINE SolverEquations_AssertIsDynamic(solverEquations,err,error,*)

    !Argument Variables
    TYPE(SolverEquationsType), POINTER, INTENT(INOUT) :: solverEquations !<The solver equations to assert the dynamic time dependence for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverEquations_AssertIsDynamic",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
#endif    

    IF(solverEquations%timeDependence/=SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC.AND. &
      & solverEquations%timeDependence/=SOLVER_EQUATIONS_SECOND_ORDER_DYNAMIC) THEN
      localError="The solver equations time dependence type of "// &
        & TRIM(NumberToVString(solverEquations%timeDependence,"*",err,error))// &
        & " does not correspond to first or second order dynamic equations."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("SolverEquations_AssertIsDynamic")
    RETURN
999 ERRORSEXITS("SolverEquations_AssertIsDynamic",err,error)
    RETURN 1
    
  END SUBROUTINE SolverEquations_AssertIsDynamic

  !
  !================================================================================================================================
  !

  !>Assert that solver equations time dependence is static
  SUBROUTINE SolverEquations_AssertIsStatic(solverEquations,err,error,*)

    !Argument Variables
    TYPE(SolverEquationsType), POINTER, INTENT(INOUT) :: solverEquations !<The solver equations to assert the static time dependence for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverEquations_AssertIsStatic",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
#endif    

    IF(solverEquations%timeDependence/=SOLVER_EQUATIONS_STATIC.AND. &
      & solverEquations%timeDependence/=SOLVER_EQUATIONS_QUASISTATIC) THEN
      localError="The solver equations time dependence type of "// &
        & TRIM(NumberToVString(solverEquations%timeDependence,"*",err,error))// &
        & " does not correspond to static or quasistatic equations."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("SolverEquations_AssertIsStatic")
    RETURN
999 ERRORSEXITS("SolverEquations_AssertIsStatic",err,error)
    RETURN 1
    
  END SUBROUTINE SolverEquations_AssertIsStatic

  !
  !================================================================================================================================
  !

  !>Gets the boundary conditions for solver equations. \see OpenCMISS::Iron::cmfe_SolverEquations_BoundaryConditionsGet
  SUBROUTINE SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to get the boundary conditions for
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<On exit, a pointer to the boundary conditions for the specified solver equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEquations_BoundaryConditionsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is already associated.",err,error,*998)
    CALL SolverEquations_AssertIsFinished(solverEquations,err,error,*999)
#endif    

    boundaryConditions=>solverEquations%boundaryConditions

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Solver equations boundary conditions is not associated.", &
      & err,error,*999)
#endif    
 
    EXITS("SolverEquations_BoundaryConditionsGet")
    RETURN
999 NULLIFY(boundaryConditions)
998 ERRORSEXITS("SolverEquations_BoundaryConditionsGet",err,error)
    RETURN 1

  END SUBROUTINE SolverEquations_BoundaryConditionsGet
     
  !
  !================================================================================================================================
  !

  !>Gets the linearity type for solver equations. 
  SUBROUTINE SolverEquations_LinearityTypeGet(solverEquations,linearityType,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to get the linearity type for
    INTEGER(INTG), INTENT(OUT) :: linearityType !<On exit, the solver equations linearity type \see SolverRoutines_EquationsLinearityTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEquations_LinearityTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
#endif    

    linearityType=solverEquations%linearity
 
    EXITS("SolverEquations_LinearityTypeGet")
    RETURN
999 ERRORSEXITS("SolverEquations_LinearityTypeGet",err,error)
    RETURN 1

  END SUBROUTINE SolverEquations_LinearityTypeGet
     
  !
  !================================================================================================================================
  !

  !>Gets the solver for solver equations. 
  SUBROUTINE SolverEquations_SolverGet(solverEquations,solver,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to get the solver for
    TYPE(SolverType), POINTER :: solver !<On exit, a pointer to the solver for the specified solver equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEquations_SolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solver)) CALL FlagError("Solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
#endif    

    solver=>solverEquations%solver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver equations solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverEquations_SolverGet")
    RETURN
998 NULLIFY(solver)
999 ERRORSEXITS("SolverEquations_SolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverEquations_SolverGet
     
  !
  !================================================================================================================================
  !

  !>Checks the solver mapping for solver equations exists. 
  SUBROUTINE SolverEquations_SolverMappingExists(solverEquations,solverMapping,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to check the solver mapping for
    TYPE(SolverMappingType), POINTER :: solverMapping !<On exit, a pointer to the solver mapping for the specified solver equations if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEquations_SolverMappingExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
#endif    

    solverMapping=>solverEquations%solverMapping
    
    EXITS("SolverEquations_SolverMappingExists")
    RETURN
998 NULLIFY(solverMapping)
999 ERRORSEXITS("SolverEquations_SolverMappingExists",err,error)
    RETURN 1

  END SUBROUTINE SolverEquations_SolverMappingExists
     
  !
  !================================================================================================================================
  !

  !>Gets the solver mapping for solver equations. 
  SUBROUTINE SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to get the solver mapping for
    TYPE(SolverMappingType), POINTER :: solverMapping !<On exit, a pointer to the solver mapping for the specified solver equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEquations_SolverMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
#endif    

    solverMapping=>solverEquations%solverMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver equations solver mapping is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverEquations_SolverMappingGet")
    RETURN
998 NULLIFY(solverMapping)
999 ERRORSEXITS("SolverEquations_SolverMappingGet",err,error)
    RETURN 1

  END SUBROUTINE SolverEquations_SolverMappingGet
     
  !
  !================================================================================================================================
  !

  !>Checks the solver matrices for solver equations exists. 
  SUBROUTINE SolverEquations_SolverMatricesExists(solverEquations,solverMatrices,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to check the solver matrices for
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<On exit, a pointer to the solver matrices for the specified solver equations if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEquations_SolverMatricesExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is already associated.",err,error,*998)
    CALL SolverEquations_AssertIsFinished(solverEquations,err,error,*999)
#endif    
 
    solverMatrices=>solverEquations%solverMatrices
 
    EXITS("SolverEquations_SolverMatricesExists")
    RETURN
998 NULLIFY(solverMatrices)
999 ERRORSEXITS("SolverEquations_SolverMatricesExists",err,error)
    RETURN 1

  END SUBROUTINE SolverEquations_SolverMatricesExists
     
  !
  !================================================================================================================================
  !

  !>Gets the solver matrices for solver equations. 
  SUBROUTINE SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to get the solver matrices for
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<On exit, a pointer to the solver matrices for the specified solver equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEquations_SolverMatricesGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is already associated.",err,error,*998)
    CALL SolverEquations_AssertIsFinished(solverEquations,err,error,*999)
#endif    
 
    solverMatrices=>solverEquations%solverMatrices

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver equations solver matrices is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverEquations_SolverMatricesGet")
    RETURN
998 NULLIFY(solverMatrices)
999 ERRORSEXITS("SolverEquations_SolverMatricesGet",err,error)
    RETURN 1

  END SUBROUTINE SolverEquations_SolverMatricesGet
     
  !
  !================================================================================================================================
  !

  !>Gets the sparsity type for solver equations. 
  SUBROUTINE SolverEquations_SparsityTypeGet(solverEquations,sparsityType,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to get the sparsity type for
    INTEGER(INTG), INTENT(OUT) :: sparsityType !<On exit, the solver equations sparsity type \see SolverRoutines_EquationsSparsityTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEquations_SparsityTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
#endif    

    sparsityType=solverEquations%sparsityType
 
    EXITS("SolverEquations_SparsityTypeGet")
    RETURN
999 ERRORSEXITS("SolverEquations_SparsityTypeGet",err,error)
    RETURN 1

  END SUBROUTINE SolverEquations_SparsityTypeGet
     
  !
  !================================================================================================================================
  !

  !>Gets the symmetry type for solver equations. 
  SUBROUTINE SolverEquations_SymmetryTypeGet(solverEquations,symmetryType,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to get the symmetry type for
    INTEGER(INTG), INTENT(OUT) :: symmetryType !<On exit, the solver equations symmetry type \see SolverRoutines_EquationsSymmetryTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEquations_SymmetryTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
#endif    

    symmetryType=solverEquations%symmetryType
 
    EXITS("SolverEquations_SymmetryTypeGet")
    RETURN
999 ERRORSEXITS("SolverEquations_SymmetryTypeGet",err,error)
    RETURN 1

  END SUBROUTINE SolverEquations_SymmetryTypeGet
     
  !
  !================================================================================================================================
  !

  !>Gets the time dependence type for solver equations. 
  SUBROUTINE SolverEquations_TimeDependenceTypeGet(solverEquations,timeDependenceType,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to get the time dependence type for
    INTEGER(INTG), INTENT(OUT) :: timeDependenceType !<On exit, the solver equations time dependence type \see SolverRoutines_EquationsTimeDependenceTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEquations_TimeDependenceTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
#endif    

    timeDependenceType=solverEquations%timeDependence
 
    EXITS("SolverEquations_TimeDependenceTypeGet")
    RETURN
999 ERRORSEXITS("SolverEquations_TimeDependenceTypeGet",err,error)
    RETURN 1

  END SUBROUTINE SolverEquations_TimeDependenceTypeGet
     
  !
  !================================================================================================================================
  !

  !>Gets the arbitrary path for a geometric transformation solver. 
  SUBROUTINE SolverGeometricTransformation_ArbitraryPathGet(geometricTransformationSolver,arbitraryPath,err,error,*)

    !Argument variables
    TYPE(GeometricTransformationSolverType), POINTER :: geometricTransformationSolver !<A pointer to the geometric transformation solver to get the arbitrary path for
    LOGICAL, INTENT(OUT) :: arbitraryPath !<On exit, the arbitrary path flag for the specified geometric transformation solver.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverGeometricTransformation_ArbitraryPathGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(geometricTransformationSolver)) &
      & CALL FlagError("Geometric transformation solver is not associated.",err,error,*999)
#endif    

    arbitraryPath=geometricTransformationSolver%arbitraryPath
 
    EXITS("SolverGeometricTransformation_ArbitraryPathGet")
    RETURN
999 ERRORSEXITS("SolverGeometricTransformation_ArbitraryPathGet",err,error)
    RETURN 1

  END SUBROUTINE SolverGeometricTransformation_ArbitraryPathGet
     
  !
  !================================================================================================================================
  !

  !>Gets the field and variable type for a geometric transformation solver. 
  SUBROUTINE SolverGeometricTransformation_FieldGet(geometricTransformationSolver,field,variableType,err,error,*)

    !Argument variables
    TYPE(GeometricTransformationSolverType), POINTER :: geometricTransformationSolver !<A pointer to the geometric transformation solver to get the field for
    TYPE(FieldType), POINTER :: field !<On exit, a pointer to the field for the specified geometric transformation solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: variableType !<On exit, the field variable type for the specified geometric transformation solver.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverGeometricTransformation_FieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(field)) CALL FlagError("Field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(geometricTransformationSolver)) &
      & CALL FlagError("Geometric transformation solver is not associated.",err,error,*999)
#endif    

    field=>geometricTransformationSolver%field
    variableType=geometricTransformationSolver%fieldVariableType

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Geometric transformation solver field is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverGeometricTransformation_FieldGet")
    RETURN
999 NULLIFY(field)
998 ERRORSEXITS("SolverGeometricTransformation_FieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverGeometricTransformation_FieldGet
     
  !
  !================================================================================================================================
  !

  !>Gets the field variable for a geometric transformation solver. 
  SUBROUTINE SolverGeometricTransformation_FieldVariableGet(geometricTransformationSolver,fieldVariable,err,error,*)

    !Argument variables
    TYPE(GeometricTransformationSolverType), POINTER :: geometricTransformationSolver !<A pointer to the geometric transformation solver to get the field variablefor
    TYPE(FieldVariableType), POINTER :: fieldVariable !<On exit, a pointer to the field variable for the specified geometric transformation solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverGeometricTransformation_FieldVariableGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(geometricTransformationSolver)) &
      & CALL FlagError("Geometric transformation solver is not associated.",err,error,*999)
#endif    

    fieldVariable=>geometricTransformationSolver%fieldVariable

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(fieldVariable)) &
      & CALL FlagError("Geometric transformation solver field variable is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverGeometricTransformation_FieldVariableGet")
    RETURN
999 NULLIFY(fieldVariable)
998 ERRORSEXITS("SolverGeometricTransformation_FieldVariableGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverGeometricTransformation_FieldVariableGet
     
  !
  !================================================================================================================================
  !

  !>Gets the number of increments for a geometric transformation solver. 
  SUBROUTINE SolverGeometricTransformation_NumberOfIncrementsGet(geometricTransformationSolver,numberOfIncrements,err,error,*)

    !Argument variables
    TYPE(GeometricTransformationSolverType), POINTER :: geometricTransformationSolver !<A pointer to the geometric transformation solver to get the number of increments for
    INTEGER(INTG), INTENT(OUT) :: numberOfIncrements !<On exit, the number of increments for the specified geometric transformation solver.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverGeometricTransformation_NumberOfIncrementsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(geometricTransformationSolver)) &
      & CALL FlagError("Geometric transformation solver is not associated.",err,error,*999)
#endif    

    numberOfIncrements=geometricTransformationSolver%numberOfIncrements
 
    EXITS("SolverGeometricTransformation_NumberOfIncrementsGet")
    RETURN
999 ERRORSEXITS("SolverGeometricTransformation_NumberOfIncrementsGet",err,error)
    RETURN 1

  END SUBROUTINE SolverGeometricTransformation_NumberOfIncrementsGet
     
  !
  !================================================================================================================================
  !

  !>Gets the solver for a geometric transformation solver. 
  SUBROUTINE SolverGeometricTransformation_SolverGet(geometricTransformationSolver,solver,err,error,*)

    !Argument variables
    TYPE(GeometricTransformationSolverType), POINTER :: geometricTransformationSolver !<A pointer to the geometric transformation solver to get the solver for
    TYPE(SolverType), POINTER :: solver !<On exit, a pointer to the solver for the specified geometric transformation solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverGeometricTransformation_SolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solver)) CALL FlagError("Solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(geometricTransformationSolver)) &
      & CALL FlagError("Geometric transformation solver is not associated.",err,error,*999)
#endif    

    solver=>geometricTransformationSolver%solver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Geometric transformation solver solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverGeometricTransformation_SolverGet")
    RETURN
999 NULLIFY(solver)
998 ERRORSEXITS("SolverGeometricTransformation_SolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverGeometricTransformation_SolverGet
     
  !
  !================================================================================================================================
  !

  !>Assert that a linear solver is a direct linear solver
  SUBROUTINE SolverLinear_AssertIsDirect(linearSolver,err,error,*)

    !Argument Variables
    TYPE(LinearSolverType), POINTER, INTENT(INOUT) :: linearSolver !<The linear solver to assert the direct solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverLinear_AssertIsDirect",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is not associated.",err,error,*999)
#endif    

    IF(linearSolver%linearSolveType/=SOLVER_LINEAR_DIRECT_SOLVE_TYPE) THEN
      localError="The linear solver type of "//TRIM(NumberToVString(linearSolver%linearSolveType,"*",err,error))// &
        & " does not correspond to the required linear direct solver type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("SolverLinear_AssertIsDirect")
    RETURN
999 ERRORSEXITS("SolverLinear_AssertIsDirect",err,error)
    RETURN 1
    
  END SUBROUTINE SolverLinear_AssertIsDirect

  !
  !================================================================================================================================
  !

  !>Assert that a linear solver is a iterative linear solver
  SUBROUTINE SolverLinear_AssertIsIterative(linearSolver,err,error,*)

    !Argument Variables
    TYPE(LinearSolverType), POINTER, INTENT(INOUT) :: linearSolver !<The linear solver to assert the iterative solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverLinear_AssertIsIterative",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is not associated.",err,error,*999)
#endif    

    IF(linearSolver%linearSolveType/=SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE) THEN
      localError="The linear solver type of "//TRIM(NumberToVString(linearSolver%linearSolveType,"*",err,error))// &
        & " does not correspond to the required linear iterative solver type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("SolverLinear_AssertIsIterative")
    RETURN
999 ERRORSEXITS("SolverLinear_AssertIsIterative",err,error)
    RETURN 1
    
  END SUBROUTINE SolverLinear_AssertIsIterative

  !
  !================================================================================================================================
  !

  !>Gets the direct solver for a linear solver. 
  SUBROUTINE SolverLinear_DirectSolverGet(linearSolver,directSolver,err,error,*)

    !Argument variables
    TYPE(LinearSolverType), POINTER :: linearSolver !<A pointer to the linear solver to get the direct solver for
    TYPE(LinearDirectSolverType), POINTER :: directSolver !<On exit, a pointer to the direct solver for the specified linear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverLinear_DirectSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(directSolver)) CALL FlagError("Direct solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is not associated.",err,error,*999)
#endif    

    directSolver=>linearSolver%directSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(directSolver)) CALL FlagError("Linear solver direct solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverLinear_DirectSolverGet")
    RETURN
999 NULLIFY(directSolver)
998 ERRORSEXITS("SolverLinear_DirectSolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverLinear_DirectSolverGet
     
  !
  !================================================================================================================================
  !

  !>Gets the iterative solver for a linear solver. 
  SUBROUTINE SolverLinear_IterativeSolverGet(linearSolver,iterativeSolver,err,error,*)

    !Argument variables
    TYPE(LinearSolverType), POINTER :: linearSolver !<A pointer to the linear solver to get the iterative solver for
    TYPE(LinearIterativeSolverType), POINTER :: iterativeSolver !<On exit, a pointer to the iterative solver for the specified linear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverLinear_IterativeSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(iterativeSolver)) CALL FlagError("Iterative solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is not associated.",err,error,*999)
#endif    

    iterativeSolver=>linearSolver%iterativeSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(iterativeSolver)) CALL FlagError("Linear solver iterative solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverLinear_IterativeSolverGet")
    RETURN
999 NULLIFY(iterativeSolver)
998 ERRORSEXITS("SolverLinear_IterativeSolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverLinear_IterativeSolverGet
     
  !
  !================================================================================================================================
  !

  !>Returns the type of library to use for a linear solver.
  SUBROUTINE SolverLinear_LibraryTypeGet(linearSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(LinearSolverType), POINTER :: linearSolver !<A pointer the linear solver to get the library type for.
     INTEGER(INTG), INTENT(OUT) :: solverLibraryType !<On exit, the type of library used for the linear solver \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(LinearDirectSolverType), POINTER :: directSolver
    TYPE(LinearIterativeSolverType), POINTER :: iterativeSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverLinear_LibraryTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is not associated.",err,error,*999)
#endif    
    
    SELECT CASE(linearSolver%linearSolveType)
    CASE(SOLVER_LINEAR_DIRECT_SOLVE_TYPE)
      NULLIFY(directSolver)
      CALL SolverLinear_DirectSolverGet(linearSolver,directSolver,err,error,*999)
      CALL SolverLinearDirect_LibraryTypeGet(directSolver,solverLibraryType,err,error,*999)
    CASE(SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE)
      NULLIFY(iterativeSolver)
      CALL SolverLinear_IterativeSolverGet(linearSolver,iterativeSolver,err,error,*999)
      CALL SolverLinearIterative_LibraryTypeGet(iterativeSolver,solverLibraryType,err,error,*999)
    CASE DEFAULT
      localError="The linear solver type of "//TRIM(NumberToVString(linearSolver%linearSolveType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("SolverLinear_LibraryTypeGet")
    RETURN
999 ERRORSEXITS("SolverLinear_LibraryTypeGet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverLinear_LibraryTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the solver for a linear solver. 
  SUBROUTINE SolverLinear_SolverGet(linearSolver,solver,err,error,*)

    !Argument variables
    TYPE(LinearSolverType), POINTER :: linearSolver !<A pointer to the linear solver to get the solver for
    TYPE(SolverType), POINTER :: solver !<On exit, a pointer to the solver for the specified linear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverLinear_SolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solver)) CALL FlagError("Solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is not associated.",err,error,*999)
#endif    

    solver=>linearSolver%solver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Linear solver solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverLinear_SolverGet")
    RETURN
999 NULLIFY(solver)
998 ERRORSEXITS("SolverLinear_SolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverLinear_SolverGet
     
  !
  !================================================================================================================================
  !

  !>Returns the solve type of a linear solver.
  SUBROUTINE SolverLinear_SolverTypeGet(linearSolver,linearSolveType,err,error,*)

    !Argument variables
    TYPE(LinearSolverType), POINTER :: linearSolver !<A pointer to the linear solver to get the solve type for
    INTEGER(INTG), INTENT(OUT) :: linearSolveType !<On return, the linear solver type. \see SolverRoutines_LinearSolverTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverLinear_SolverTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is not associated.",err,error,*999)
#endif    

    linearSolveType=linearSolver%linearSolveType
    
    EXITS("SolverLinear_SolverTypeGet")
    RETURN
999 ERRORSEXITS("SolverLinear_SolverTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverLinear_SolverTypeGet

  !
  !================================================================================================================================
  !

  !>Returns the type of library to use for a direct linear solver.
  SUBROUTINE SolverLinearDirect_LibraryTypeGet(directSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(LinearDirectSolverType), POINTER :: directSolver !<A pointer the direct linear solver to get the library type for.
    INTEGER(INTG), INTENT(OUT) :: solverLibraryType !<On exit, the type of library used for the direct linear solver \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverLinearDirect_LibraryTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(directSolver)) CALL FlagError("Direct linear solver is not associated.",err,error,*999)
#endif    
    
    SELECT CASE(directSolver%directSolverType)
    CASE(SOLVER_DIRECT_LU)
      solverLibraryType=directSolver%solverLibrary
    CASE(SOLVER_DIRECT_CHOLESKY)
      solverLibraryType=directSolver%solverLibrary
    CASE(SOLVER_DIRECT_SVD)
      solverLibraryType=directSolver%solverLibrary
    CASE DEFAULT
      localError="The direct linear solver type of "//TRIM(NumberToVString(directSolver%directSolverType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("SolverLinearDirect_LibraryTypeGet")
    RETURN
999 ERRORSEXITS("SolverLinearDirect_LibraryTypeGet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverLinearDirect_LibraryTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the linear solver for a direct linear solver. 
  SUBROUTINE SolverLinearDirect_LinearSolverGet(directSolver,linearSolver,err,error,*)

    !Argument variables
    TYPE(LinearDirectSolverType), POINTER :: directSolver !<A pointer to the direct solver to get the linear solver for
    TYPE(LinearSolverType), POINTER :: linearSolver !<On exit, a pointer to the linear solver for the specified direct solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverLinearDirect_LinearSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(directSolver)) CALL FlagError("Direct linear solver is not associated.",err,error,*999)
#endif    

    linearSolver=>directSolver%linearSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("Linear direct solver linear solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverLinearDirect_LinearSolverGet")
    RETURN
999 NULLIFY(linearSolver)
998 ERRORSEXITS("SolverLinearDirect_LinearSolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverLinearDirect_LinearSolverGet
     
  !
  !================================================================================================================================
  !

  !>Returns the type of library to use for an iterative linear solver.
  SUBROUTINE SolverLinearIterative_LibraryTypeGet(iterativeSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(LinearIterativeSolverType), POINTER :: iterativeSolver !<A pointer the iterative linear solver to get the library type for.
    INTEGER(INTG), INTENT(OUT) :: solverLibraryType !<On exit, the type of library used for the iterative linear solver \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverLinearIterative_LibraryTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(iterativeSolver)) CALL FlagError("Iterative linear solver is not associated.",err,error,*999)
#endif    
    
    SELECT CASE(iterativeSolver%iterativeSolverType)
    CASE(SOLVER_ITERATIVE_RICHARDSON)
      solverLibraryType=iterativeSolver%solverLibrary
    CASE(SOLVER_ITERATIVE_CHEBYSHEV)
      solverLibraryType=iterativeSolver%solverLibrary
    CASE(SOLVER_ITERATIVE_CONJUGATE_GRADIENT)
      solverLibraryType=iterativeSolver%solverLibrary
    CASE(SOLVER_ITERATIVE_GMRES)
      solverLibraryType=iterativeSolver%solverLibrary
    CASE(SOLVER_ITERATIVE_BiCGSTAB)
      solverLibraryType=iterativeSolver%solverLibrary
    CASE(SOLVER_ITERATIVE_CONJGRAD_SQUARED)
      solverLibraryType=iterativeSolver%solverLibrary
    CASE DEFAULT
      localError="The iterative linear solver type of "// &
        & TRIM(NumberToVString(iterativeSolver%iterativeSolverType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("SolverLinearIterative_LibraryTypeGet")
    RETURN
999 ERRORSEXITS("SolverLinearIterative_LibraryTypeGet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverLinearIterative_LibraryTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the linear solver for a iterative linear solver. 
  SUBROUTINE SolverLinearIterative_LinearSolverGet(iterativeSolver,linearSolver,err,error,*)

    !Argument variables
    TYPE(LinearIterativeSolverType), POINTER :: iterativeSolver !<A pointer to the iterative solver to get the linear solver for
    TYPE(LinearSolverType), POINTER :: linearSolver !<On exit, a pointer to the linear solver for the specified iterative solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverLinearIterative_LinearSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(iterativeSolver)) CALL FlagError("Iterative linear solver is not associated.",err,error,*999)
#endif    

    linearSolver=>iterativeSolver%linearSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("Linear iterative solver linear solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverLinearIterative_LinearSolverGet")
    RETURN
999 NULLIFY(linearSolver)
998 ERRORSEXITS("SolverLinearIterative_LinearSolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverLinearIterative_LinearSolverGet
     
  !
  !================================================================================================================================
  !

  !>Assert that a nonlinear solver is a Newton solver
  SUBROUTINE SolverNonlinear_AssertIsNewton(nonlinearSolver,err,error,*)

    !Argument Variables
    TYPE(NonlinearSolverType), POINTER, INTENT(INOUT) :: nonlinearSolver !<The nonlinear solver to assert the Newton solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverNonlinear_AssertIsNewton",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is not associated.",err,error,*999)
#endif    

    IF(nonlinearSolver%nonlinearSolveType/=SOLVER_NONLINEAR_NEWTON) THEN
      localError="The nonlinear solver type of "//TRIM(NumberToVString(nonlinearSolver%nonlinearSolveType,"*",err,error))// &
        & " does not correspond to the required Newton nonlinear solve type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("SolverNonlinear_AssertIsNewton")
    RETURN
999 ERRORSEXITS("SolverNonlinear_AssertIsNewton",err,error)
    RETURN 1
    
  END SUBROUTINE SolverNonlinear_AssertIsNewton

  !
  !================================================================================================================================
  !

  !>Assert that a nonlinear solver is a quasi Newton solver
  SUBROUTINE SolverNonlinear_AssertIsQuasiNewton(nonlinearSolver,err,error,*)

    !Argument Variables
    TYPE(NonlinearSolverType), POINTER, INTENT(INOUT) :: nonlinearSolver !<The nonlinear solver to assert the quasi Newton solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverNonlinear_AssertIsQuasiNewton",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is not associated.",err,error,*999)
#endif    

    IF(nonlinearSolver%nonlinearSolveType/=SOLVER_NONLINEAR_QUASI_NEWTON) THEN
      localError="The nonlinear solver type of "//TRIM(NumberToVString(nonlinearSolver%nonlinearSolveType,"*",err,error))// &
        & " does not correspond to the required quasi Newton nonlinear solve type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("SolverNonlinear_AssertIsQuasiNewton")
    RETURN
999 ERRORSEXITS("SolverNonlinear_AssertIsQuasiNewton",err,error)
    RETURN 1
    
  END SUBROUTINE SolverNonlinear_AssertIsQuasiNewton

  !
  !================================================================================================================================
  !

  !>Returns the type of library to use for a nonlinear solver.
  SUBROUTINE SolverNonlinear_LibraryTypeGet(nonlinearSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver !<A pointer the nonlinear solver to get the library type for.
    INTEGER(INTG), INTENT(OUT) :: solverLibraryType !<On exit, the type of library used for the nonlinear solver \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverNonlinear_LibraryTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is not associated.",err,error,*999)
#endif    
    
    SELECT CASE(nonlinearSolver%nonlinearSolveType)
    CASE(SOLVER_NONLINEAR_NEWTON)
      NULLIFY(newtonSolver)
      CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
      CALL SolverNonlinearNewton_LibraryTypeGet(newtonSolver,solverLibraryType,err,error,*999)
    CASE(SOLVER_NONLINEAR_BFGS_INVERSE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_NONLINEAR_SQP)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
      NULLIFY(quasiNewtonSolver)
      CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
      CALL SolverNonlinearQuasiNewton_LibraryTypeGet(quasiNewtonSolver,solverLibraryType,err,error,*999)
   CASE DEFAULT
      localError="The nonlinear solver type of "// &
        & TRIM(NumberToVString(nonlinearSolver%nonlinearSolveType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("SolverNonlinear_LibraryTypeGet")
    RETURN
999 ERRORSEXITS("SolverNonlinear_LibraryTypeGet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverNonlinear_LibraryTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the linear solver for a nonlinear solver. 
  SUBROUTINE SolverNonlinear_LinearSolverGet(nonlinearSolver,linearSolver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: nonlinearSolver !<A pointer to the nonlinear solver to get the linear solver for
    TYPE(SolverType), POINTER :: linearSolver !<On exit, a pointer to the linear solver for the specified nonlinear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverNonlinear_SolverLinearSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(nonlinearSolver%nonlinearSolver)) &
      & CALL FlagError("Nonlinear solver nonlinear solver is not associated.",err,error,*999)
#endif    

    SELECT CASE(nonlinearSolver%nonlinearSolver%nonlinearSolveType)
    CASE(SOLVER_NONLINEAR_NEWTON)
#ifdef WITH_PRECHECKS      
      IF(.NOT.ASSOCIATED(nonlinearSolver%nonlinearSolver%newtonSolver)) &
        & CALL FlagError("Nonlinear solver Newton solver is not associated.",err,error,*999)
#endif      
      linearSolver=>nonlinearSolver%nonlinearSolver%newtonSolver%linearSolver
    CASE(SOLVER_NONLINEAR_BFGS_INVERSE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_NONLINEAR_SQP)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
#ifdef WITH_PRECHECKS      
      IF(.NOT.ASSOCIATED(nonlinearSolver%nonlinearSolver%quasiNewtonSolver)) &
        & CALL FlagError("Nonlinear solver quasi-Newton solver is not associated.",err,error,*999)
#endif
      linearSolver=>nonlinearSolver%nonlinearSolver%quasiNewtonSolver%linearSolver
    CASE DEFAULT
      localError="The nonlinear solver type of "// &
        & TRIM(NumberToVString(nonlinearSolver%nonlinearSolver%nonlinearSolveType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("Nonlinear solver linear solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverNonlinear_LinearSolverGet")
    RETURN
999 NULLIFY(linearSolver)
998 ERRORSEXITS("SolverNonlinear_LinearSolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverNonlinear_LinearSolverGet
     
  !
  !================================================================================================================================
  !

  !>Gets the newton solver for a nonlinear solver. 
  SUBROUTINE SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*)

    !Argument variables
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver !<A pointer to the nonlinear solver to get the Newton solver for
    TYPE(NewtonSolverType), POINTER :: newtonSolver !<On exit, a pointer to the newton solver for the specified nonlinear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinear_NewtonSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(newtonSolver)) CALL FlagError("Newton solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is not associated.",err,error,*999)
#endif    

    newtonSolver=>nonlinearSolver%newtonSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(newtonSolver)) CALL FlagError("Nonlinear solver Newton solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverNonlinear_NewtonSolverGet")
    RETURN
999 NULLIFY(newtonSolver)
998 ERRORSEXITS("SolverNonlinear_NewtonSolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverNonlinear_NewtonSolverGet
     
  !
  !================================================================================================================================
  !

  !>Gets the quasi Newton solver for a nonlinear solver. 
  SUBROUTINE SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*)

    !Argument variables
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver !<A pointer to the nonlinear solver to get the quasi Newton solver for
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver !<On exit, a pointer to the quasi Newton solver for the specified nonlinear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinear_QuasiNewtonSolverGet",err,error,*998)
    
#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(quasiNewtonSolver)) CALL FlagError("Quasi Newton solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is not associated.",err,error,*999)
#endif    

    quasiNewtonSolver=>nonlinearSolver%quasiNewtonSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(quasiNewtonSolver)) CALL FlagError("Nonlinear solver quasi Newton solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverNonlinear_QuasiNewtonSolverGet")
    RETURN
999 NULLIFY(quasiNewtonSolver)
998 ERRORSEXITS("SolverNonlinear_QuasiNewtonSolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverNonlinear_QuasiNewtonSolverGet
     
  !
  !================================================================================================================================
  !

  !>Gets the solver for a non-linear solver. 
  SUBROUTINE SolverNonlinear_SolverGet(nonlinearSolver,solver,err,error,*)

    !Argument variables
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver !<A pointer to the nonlinear solver to get the solver for
    TYPE(SolverType), POINTER :: solver !<On exit, a pointer to the solver for the specified nonlinear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinear_SolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solver)) CALL FlagError("Solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is not associated.",err,error,*999)
#endif    

    solver=>nonlinearSolver%solver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Nonlinear solver solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverNonlinear_SolverGet")
    RETURN
999 NULLIFY(solver)
998 ERRORSEXITS("SolverNonlinear_SolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverNonlinear_SolverGet
     
  !
  !================================================================================================================================
  !

  !>Returns the solve type of a nonlinear solver.
  SUBROUTINE SolverNonlinear_SolverTypeGet(nonlinearSolver,nonlinearSolveType,err,error,*)

    !Argument variables
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver !<A pointer to the nonlinear solver to get the solve type for
    INTEGER(INTG), INTENT(OUT) :: nonlinearSolveType !<On return, the nonlinear solver type. \see SolverRoutines_NonlinearSolverTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinear_SolverTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is not associated.",err,error,*999)
#endif    

    nonlinearSolveType=nonlinearSolver%nonlinearSolveType
    
    EXITS("SolverNonlinear_SolverTypeGet")
    RETURN
999 ERRORSEXITS("SolverNonlinear_SolverTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverNonlinear_SolverTypeGet

  !
  !================================================================================================================================
  !

  !>Assert that a nonlinear Newton solver is a linesearch solver
  SUBROUTINE SolverNonlinearNewton_AssertIsLinesearch(newtonSolver,err,error,*)

    !Argument Variables
    TYPE(NewtonSolverType), POINTER, INTENT(INOUT) :: newtonSolver !<The nonlinear Newton solver to assert the linesearch solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverNonlinearNewton_AssertIsLinesearch",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(newtonSolver)) CALL FlagError("Nonlinear Newwton solver is not associated.",err,error,*999)
#endif    

    IF(newtonSolver%newtonSolveType/=SOLVER_NEWTON_LINESEARCH) THEN
      localError="The nonlinear Newton solver type of "//TRIM(NumberToVString(newtonSolver%newtonSolveType,"*",err,error))// &
        & " does not correspond to the required Newton linesearch nonlinear solve type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("SolverNonlinearNewton_AssertIsLinesearch")
    RETURN
999 ERRORSEXITS("SolverNonlinearNewton_AssertIsLinesearch",err,error)
    RETURN 1
    
  END SUBROUTINE SolverNonlinearNewton_AssertIsLinesearch

  !
  !================================================================================================================================
  !

  !>Assert that a nonlinear Newton solver is a trustregion solver
  SUBROUTINE SolverNonlinearNewton_AssertIsTrustregion(newtonSolver,err,error,*)

    !Argument Variables
    TYPE(NewtonSolverType), POINTER, INTENT(INOUT) :: newtonSolver !<The nonlinear Newton solver to assert the trustregion solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverNonlinearNewton_AssertIsTrustregion",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(newtonSolver)) CALL FlagError("Nonlinear Newwton solver is not associated.",err,error,*999)
#endif    

    IF(newtonSolver%newtonSolveType/=SOLVER_NEWTON_TRUSTREGION) THEN
      localError="The nonlinear Newton solver type of "//TRIM(NumberToVString(newtonSolver%newtonSolveType,"*",err,error))// &
        & " does not correspond to the required Newton trustregion nonlinear solve type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("SolverNonlinearNewton_AssertIsTrustregion")
    RETURN
999 ERRORSEXITS("SolverNonlinearNewton_AssertIsTrustregion",err,error)
    RETURN 1
    
  END SUBROUTINE SolverNonlinearNewton_AssertIsTrustregion

  !
  !================================================================================================================================
  !

  !>Gets the convergence test for a Newton solver. 
  SUBROUTINE SolverNonlinearNewton_ConvergenceTestGet(newtonSolver,convergenceTest,err,error,*)

    !Argument variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver !<A pointer to the Newton nonlinear solver to get the convergence test for
    TYPE(NewtonSolverConvergenceTestType), POINTER :: convergenceTest !<On exit, a pointer to the convergence test for the specified Newton nonlinear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinearNewton_ConvergenceTestGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(convergenceTest)) CALL FlagError("Convergence test is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(newtonSolver)) CALL FlagError("Newton nonlinear solver is not associated.",err,error,*999)
#endif    

    convergenceTest=>newtonSolver%convergenceTest

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(convergenceTest)) &
      & CALL FlagError("Newton nonlinear solver convergence test is not associated.",err,error,*999)
#endif
 
    EXITS("SolverNonlinearNewton_ConvergenceTestGet")
    RETURN
999 NULLIFY(convergenceTest)
998 ERRORS("SolverNonlinearNewton_ConvergenceTestGet",err,error)
    EXITS("SolverNonlinearNewton_ConvergenceTestGet")
    RETURN 1

  END SUBROUTINE SolverNonlinearNewton_ConvergenceTestGet
     
  !
  !================================================================================================================================
  !

  !>Returns the type of library to use for a Newton solver.
  SUBROUTINE SolverNonlinearNewton_LibraryTypeGet(newtonSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver !<A pointer the Newton solver to get the library type for.
    INTEGER(INTG), INTENT(OUT) :: solverLibraryType !<On exit, the type of library used for the Newton solver \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(NewtonLinesearchSolverType), POINTER :: linesearchSolver
    TYPE(NewtonTrustregionSolverType), POINTER :: trustregionSolver
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverNonlinearNewton_LibraryTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(newtonSolver)) CALL FlagError("Newton solver is not associated.",err,error,*999)
#endif    
    
    SELECT CASE(newtonSolver%newtonSolveType)
    CASE(SOLVER_NEWTON_LINESEARCH)
#ifdef WITH_PRECHECKS      
      NULLIFY(linesearchSolver)
      CALL SolverNonlinearNewton_LinesearchSolverGet(newtonSolver,linesearchSolver,err,error,*999)
#endif      
      solverLibraryType=newtonSolver%linesearchSolver%solverLibrary
    CASE(SOLVER_NEWTON_TRUSTREGION)
#ifdef WITH_PRECHECKS      
      NULLIFY(trustregionSolver)
      CALL SolverNonlinearNewton_TrustregionSolverGet(newtonSolver,trustregionSolver,err,error,*999)
#endif      
      solverLibraryType=newtonSolver%trustregionSolver%solverLibrary
    CASE DEFAULT
      localError="The Newton solver type of "// &
        & TRIM(NumberToVString(newtonSolver%newtonSolveType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("SolverNonlinearNewton_LibraryTypeGet")
    RETURN
999 ERRORSEXITS("SolverNonlinearNewton_LibraryTypeGet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverNonlinearNewton_LibraryTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the linesearch solver for a Newton solver. 
  SUBROUTINE SolverNonlinearNewton_LinesearchSolverGet(newtonSolver,linesearchSolver,err,error,*)

    !Argument variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver !<A pointer to the Newton nonlinear solver to get the linesearch solver for
    TYPE(NewtonLinesearchSolverType), POINTER :: linesearchSolver !<On exit, a pointer to the linesearch solver for the specified Newton nonlinear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinearNewton_LinesearchSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linesearchSolver)) CALL FlagError("Linesearch solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(newtonSolver)) CALL FlagError("Newton nonlinear solver is not associated.",err,error,*999)
#endif    

    linesearchSolver=>newtonSolver%linesearchSolver

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(linesearchSolver)) &
      & CALL FlagError("Newton nonlinear solver linesearch solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverNonlinearNewton_LinesearchSolverGet")
    RETURN
999 NULLIFY(linesearchSolver)
998 ERRORS("SolverNonlinearNewton_LinesearchSolverGet",err,error)
    EXITS("SolverNonlinearNewton_LinesearchSolverGet")
    RETURN 1

  END SUBROUTINE SolverNonlinearNewton_LinesearchSolverGet
     
  !
  !================================================================================================================================
  !

  !>Checks if the linked CellML solver for a Newton solver exists. 
  SUBROUTINE SolverNonlinearNewton_LinkedCellMLSolverExists(newtonSolver,cellMLSolver,err,error,*)

    !Argument variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver !<A pointer to the Newton nonlinear solver to check the linked CellML solver for
    TYPE(SolverType), POINTER :: cellMLSolver !<On exit, a pointer to the linked CellML solver for the specified Newton nonlinear solver if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinearNewton_LinkedCellMLSolverExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLSolver)) CALL FlagError("CellML solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(newtonSolver)) CALL FlagError("Newton nonlinear solver is not associated.",err,error,*999)
#endif    

    cellMLSolver=>newtonSolver%cellMLEvaluatorSolver

    EXITS("SolverNonlinearNewton_LinkedCellMLSolverExists")
    RETURN
999 NULLIFY(cellMLSolver)
998 ERRORS("SolverNonlinearNewton_LinkedCellMLSolverExists",err,error)
    EXITS("SolverNonlinearNewton_LinkedCellMLSolverExists")
    RETURN 1

  END SUBROUTINE SolverNonlinearNewton_LinkedCellMLSolverExists
     
  !
  !================================================================================================================================
  !

  !>Gets the linked CellML solver for a Newton solver. 
  SUBROUTINE SolverNonlinearNewton_LinkedCellMLSolverGet(newtonSolver,cellMLSolver,err,error,*)

    !Argument variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver !<A pointer to the Newton nonlinear solver to get the linked CellML solver for
    TYPE(SolverType), POINTER :: cellMLSolver !<On exit, a pointer to the linked CellML solver for the specified Newton nonlinear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinearNewton_LinkedCellMLSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLSolver)) CALL FlagError("CellML solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(newtonSolver)) CALL FlagError("Newton nonlinear solver is not associated.",err,error,*999)
#endif    

    cellMLSolver=>newtonSolver%cellMLEvaluatorSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(cellMLSolver)) &
      & CALL FlagError("Newton nonlinear solver linked CellML solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverNonlinearNewton_LinkedCellMLSolverGet")
    RETURN
999 NULLIFY(cellMLSolver)
998 ERRORS("SolverNonlinearNewton_LinkedCellMLSolverGet",err,error)
    EXITS("SolverNonlinearNewton_LinkedCellMLSolverGet")
    RETURN 1

  END SUBROUTINE SolverNonlinearNewton_LinkedCellMLSolverGet
     
  !
  !================================================================================================================================
  !

  !>Gets the linked linear solver for a Newton solver. 
  SUBROUTINE SolverNonlinearNewton_LinkedLinearSolverGet(newtonSolver,linearSolver,err,error,*)

    !Argument variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver !<A pointer to the Newton nonlinear solver to get the linked linear solver for
    TYPE(SolverType), POINTER :: linearSolver !<On exit, a pointer to the linked linear solver for the specified Newton nonlinear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinearNewton_LinkedLinearSolverGet",err,error,*998)

    IF(ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(newtonSolver)) CALL FlagError("Newton nonlinear solver is not associated.",err,error,*999)

    linearSolver=>newtonSolver%linearSolver
    IF(.NOT.ASSOCIATED(linearSolver)) &
      & CALL FlagError("Newton nonlinear solver linked linear solver is not associated.",err,error,*999)
 
    EXITS("SolverNonlinearNewton_LinkedLinearSolverGet")
    RETURN
999 NULLIFY(linearSolver)
998 ERRORS("SolverNonlinearNewton_LinkedLinearSolverGet",err,error)
    EXITS("SolverNonlinearNewton_LinkedLinearSolverGet")
    RETURN 1

  END SUBROUTINE SolverNonlinearNewton_LinkedLinearSolverGet
     
  !
  !================================================================================================================================
  !

  !>Gets the nonlinear solver for a Newton solver. 
  SUBROUTINE SolverNonlinearNewton_NonlinearSolverGet(newtonSolver,nonlinearSolver,err,error,*)

    !Argument variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver !<A pointer to the Newton nonlinear solver to get the nonlinear solver for
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver !<On exit, a pointer to the nonlinear solver for the specified Newton nonlinear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinearNewton_NonlinearSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(newtonSolver)) CALL FlagError("Newton nonlinear solver is not associated.",err,error,*999)
#endif    

    nonlinearSolver=>newtonSolver%nonlinearSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(nonlinearSolver)) &
      & CALL FlagError("Newton nonlinear solver nonlinear solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverNonlinearNewton_NonlinearSolverGet")
    RETURN
999 NULLIFY(nonlinearSolver)
998 ERRORS("SolverNonlinearNewton_NonlinearSolverGet",err,error)
    EXITS("SolverNonlinearNewton_NonlinearSolverGet")
    RETURN 1

  END SUBROUTINE SolverNonlinearNewton_NonlinearSolverGet
     
  !
  !================================================================================================================================
  !

  !>Returns the solve type of a Newton nonlinear solver.
  SUBROUTINE SolverNonlinearNewton_SolverTypeGet(newtonNonlinearSolver,newtonNonlinearSolveType,err,error,*)

    !Argument variables
    TYPE(NewtonSolverType), POINTER :: newtonNonlinearSolver !<A pointer to the Newton nonlinear solver to get the solve type for
    INTEGER(INTG), INTENT(OUT) :: newtonNonlinearSolveType !<On return, the Newton nonlinear solver type. \see SolverRoutines_NewtonSolverTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinearNewton_SolverTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(newtonNonlinearSolver)) CALL FlagError("Newton nonlinear solver is not associated.",err,error,*999)
#endif    

    newtonNonlinearSolveType=newtonNonlinearSolver%newtonSolveType
    
    EXITS("SolverNonlinearNewton_SolverTypeGet")
    RETURN
999 ERRORSEXITS("SolverNonlinearNewton_SolverTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverNonlinearNewton_SolverTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the trustregion solver for a Newton solver. 
  SUBROUTINE SolverNonlinearNewton_TrustregionSolverGet(newtonSolver,trustregionSolver,err,error,*)

    !Argument variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver !<A pointer to the Newton nonlinear solver to get the trustregion solver for
    TYPE(NewtonTrustregionSolverType), POINTER :: trustregionSolver !<On exit, a pointer to the trustregion solver for the specified Newton nonlinear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinearNewton_TrustregionSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(trustregionSolver)) CALL FlagError("Trustregion solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(newtonSolver)) CALL FlagError("Newton nonlinear solver is not associated.",err,error,*999)
#endif    

    trustregionSolver=>newtonSolver%trustregionSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(trustregionSolver)) &
      & CALL FlagError("Newton nonlinear solver trustregion solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverNonlinearNewton_TrustregionSolverGet")
    RETURN
999 NULLIFY(trustregionSolver)
998 ERRORS("SolverNonlinearNewton_TrustregionSolverGet",err,error)
    EXITS("SolverNonlinearNewton_TrustregionSolverGet")
    RETURN 1

  END SUBROUTINE SolverNonlinearNewton_TrustregionSolverGet
     
  !
  !================================================================================================================================
  !

  !>Gets the Newton solver for a Newton linesearch solver. 
  SUBROUTINE SolverNonlinearNewtonLinesearch_NewtonSolverGet(linesearchSolver,newtonSolver,err,error,*)

    !Argument variables
    TYPE(NewtonLinesearchSolverType), POINTER :: linesearchSolver !<A pointer to the Newton nonlinear linesearch solver to get the Newton solver for
    TYPE(NewtonSolverType), POINTER :: newtonSolver !<On exit, a pointer to the Newton solver for the specified Newton linesearch solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinearNewtonLinesearch_NewtonSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(newtonSolver)) CALL FlagError("Newton solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(linesearchSolver)) CALL FlagError("Linesearch solver is not associated.",err,error,*999)
#endif    

    newtonSolver=>linesearchSolver%newtonSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(newtonSolver)) &
      & CALL FlagError("Newton nonlinear linesearch solver Newton solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverNonlinearNewtonLinesearch_NewtonSolverGet")
    RETURN
999 NULLIFY(newtonSolver)
998 ERRORS("SolverNonlinearNewtonLinesearch_NewtonSolverGet",err,error)
    EXITS("SolverNonlinearNewtonLinesearch_NewtonSolverGet")
    RETURN 1

  END SUBROUTINE SolverNonlinearNewtonLinesearch_NewtonSolverGet
     
  !
  !================================================================================================================================
  !

  !>Gets the Newton solver for a Newton trustregion solver. 
  SUBROUTINE SolverNonlinearNewtonTrustregion_NewtonSolverGet(trustregionSolver,newtonSolver,err,error,*)

    !Argument variables
    TYPE(NewtonTrustregionSolverType), POINTER :: trustregionSolver !<A pointer to the Newton nonlinear trustregion solver to get the Newton solver for
    TYPE(NewtonSolverType), POINTER :: newtonSolver !<On exit, a pointer to the Newton solver for the specified Newton trustregion solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinearNewtonTrustregion_NewtonSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(newtonSolver)) CALL FlagError("Newton solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(trustregionSolver)) CALL FlagError("Trustregion solver is not associated.",err,error,*999)
#endif    

    newtonSolver=>trustregionSolver%newtonSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(newtonSolver)) &
      & CALL FlagError("Newton nonlinear trustregion solver Newton solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverNonlinearNewtonTrustregion_NewtonSolverGet")
    RETURN
999 NULLIFY(newtonSolver)
998 ERRORS("SolverNonlinearNewtonTrustregion_NewtonSolverGet",err,error)
    EXITS("SolverNonlinearNewtonTrustregion_NewtonSolverGet")
    RETURN 1

  END SUBROUTINE SolverNonlinearNewtonTrustregion_NewtonSolverGet
     
  !
  !================================================================================================================================
  !

  !>Assert that a nonlinear quasi Newton solver is a linesearch solver
  SUBROUTINE SolverNonlinearQuasiNewton_AssertIsLinesearch(quasiNewtonSolver,err,error,*)

    !Argument Variables
    TYPE(QuasiNewtonSolverType), POINTER, INTENT(INOUT) :: quasiNewtonSolver !<The nonlinear quasi Newton solver to assert the linesearch solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverNonlinearQuasiNewton_AssertIsLinesearch",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(quasiNewtonSolver)) CALL FlagError("Nonlinear quasi Newwton solver is not associated.",err,error,*999)
#endif
    
    IF(quasiNewtonSolver%quasiNewtonSolveType/=SOLVER_QUASI_NEWTON_LINESEARCH) THEN
      localError="The nonlinear quasi Newton solver type of "// &
        & TRIM(NumberToVString(quasiNewtonSolver%quasiNewtonSolveType,"*",err,error))// &
        & " does not correspond to the required quasi Newton linesearch nonlinear solve type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("SolverNonlinearQuasiNewton_AssertIsLinesearch")
    RETURN
999 ERRORSEXITS("SolverNonlinearQuasiNewton_AssertIsLinesearch",err,error)
    RETURN 1
    
  END SUBROUTINE SolverNonlinearQuasiNewton_AssertIsLinesearch

  !
  !================================================================================================================================
  !

  !>Assert that a nonlinear quasi Newton solver is a trustregion solver
  SUBROUTINE SolverNonlinearQuasiNewton_AssertIsTrustregion(quasiNewtonSolver,err,error,*)

    !Argument Variables
    TYPE(QuasiNewtonSolverType), POINTER, INTENT(INOUT) :: quasiNewtonSolver !<The nonlinear quasi Newton solver to assert the trustregion solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverNonlinearQuasiNewton_AssertIsTrustregion",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(quasiNewtonSolver)) CALL FlagError("Nonlinear quasi Newwton solver is not associated.",err,error,*999)
#endif    

    IF(quasiNewtonSolver%quasiNewtonSolveType/=SOLVER_QUASI_NEWTON_TRUSTREGION) THEN
      localError="The nonlinear quasi Newton solver type of "// &
        & TRIM(NumberToVString(quasiNewtonSolver%quasiNewtonSolveType,"*",err,error))// &
        & " does not correspond to the required quasi Newton trustregion nonlinear solve type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("SolverNonlinearQuasiNewton_AssertIsTrustregion")
    RETURN
999 ERRORSEXITS("SolverNonlinearQuasiNewton_AssertIsTrustregion",err,error)
    RETURN 1
    
  END SUBROUTINE SolverNonlinearQuasiNewton_AssertIsTrustregion

  !
  !================================================================================================================================
  !

  !>Gets the convergence test for a quasi Newton solver. 
  SUBROUTINE SolverNonlinearQuasiNewton_ConvergenceTestGet(quasiNewtonSolver,convergenceTest,err,error,*)

    !Argument variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver !<A pointer to the quasi Newton nonlinear solver to get the convergence test for
    TYPE(NewtonSolverConvergenceTestType), POINTER :: convergenceTest !<On exit, a pointer to the convergence test for the specified quasi Newton nonlinear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinearQuasiNewton_ConvergenceTestGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(convergenceTest)) CALL FlagError("Convergence test is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(quasiNewtonSolver)) CALL FlagError("Quasi Newton nonlinear solver is not associated.",err,error,*999)
#endif    

    convergenceTest=>quasiNewtonSolver%convergenceTest

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(convergenceTest)) &
      & CALL FlagError("Quasi Newton nonlinear solver convergence test is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverNonlinearQuasiNewton_ConvergenceTestGet")
    RETURN
999 NULLIFY(convergenceTest)
998 ERRORS("SolverNonlinearQuasiNewton_ConvergenceTestGet",err,error)
    EXITS("SolverNonlinearQuasiNewton_ConvergenceTestGet")
    RETURN 1

  END SUBROUTINE SolverNonlinearQuasiNewton_ConvergenceTestGet
     
  !
  !================================================================================================================================
  !

  !>Returns the type of library to use for a Quasi-Newton solver.
  SUBROUTINE SolverNonlinearQuasiNewton_LibraryTypeGet(quasiNewtonSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver !<A pointer the Quasi-Newton solver to get the library type for.
    INTEGER(INTG), INTENT(OUT) :: solverLibraryType !<On exit, the type of library used for the Quasi-Newton solver \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(QuasiNewtonLinesearchSolverType), POINTER :: linesearchSolver
    TYPE(QuasiNewtonTrustregionSolverType), POINTER :: trustregionSolver
#endif    
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverNonlinearQuasiNewton_LibraryTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(quasiNewtonSolver))  CALL FlagError("Quasi-Newton solver is not associated.",err,error,*999)
#endif    
    
    SELECT CASE(quasiNewtonSolver%quasiNewtonSolveType)
    CASE(SOLVER_QUASI_NEWTON_LINESEARCH)
#ifdef WITH_PRECHECKS      
      NULLIFY(linesearchSolver)
      CALL SolverNonlinearQuasiNewton_LinesearchSolverGet(quasiNewtonSolver,linesearchSolver,err,error,*999)
#endif      
      solverLibraryType=quasiNewtonSolver%linesearchSolver%solverLibrary
    CASE(SOLVER_QUASI_NEWTON_TRUSTREGION)
#ifdef WITH_PRECHECKS      
      NULLIFY(trustregionSolver)
      CALL SolverNonlinearQuasiNewton_TrustRegionSolverGet(quasiNewtonSolver,trustregionSolver,err,error,*999)
#endif      
      solverLibraryType=quasiNewtonSolver%trustregionSolver%solverLibrary
    CASE DEFAULT
      localError="The Quasi-Newton solver type of "// &
        & TRIM(NumberToVString(quasiNewtonSolver%quasiNewtonSolveType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("SolverNonlinearQuasiNewton_LibraryTypeGet")
    RETURN
999 ERRORS("SolverNonlinearQuasiNewton_LibraryTypeGet",err,error)
    EXITS("SolverNonlinearQuasiNewton_LibraryTypeGet")
    RETURN 1
   
  END SUBROUTINE SolverNonlinearQuasiNewton_LibraryTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the linesearch solver for a quasi Newton solver. 
  SUBROUTINE SolverNonlinearQuasiNewton_LinesearchSolverGet(quasiNewtonSolver,linesearchSolver,err,error,*)

    !Argument variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver !<A pointer to the quasi Newton nonlinear solver to get the linesearch solver for
    TYPE(QuasiNewtonLinesearchSolverType), POINTER :: linesearchSolver !<On exit, a pointer to the linesearch solver for the specified quasi Newton nonlinear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinearQuasiNewton_LinesearchSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linesearchSolver)) CALL FlagError("Linesearch solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(quasiNewtonSolver)) CALL FlagError("Quasi Newton nonlinear solver is not associated.",err,error,*999)
#endif    

    linesearchSolver=>quasiNewtonSolver%linesearchSolver

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(linesearchSolver)) &
      & CALL FlagError("Quasi Newton nonlinear solver linesearch solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverNonlinearQuasiNewton_LinesearchSolverGet")
    RETURN
999 NULLIFY(linesearchSolver)
998 ERRORS("SolverNonlinearQuasiNewton_LinesearchSolverGet",err,error)
    EXITS("SolverNonlinearQuasiNewton_LinesearchSolverGet")
    RETURN 1

  END SUBROUTINE SolverNonlinearQuasiNewton_LinesearchSolverGet
     
  !
  !================================================================================================================================
  !

  !>Checkes if the linked CellML solver for a quasi Newton solver exists. 
  SUBROUTINE SolverNonlinearQuasiNewton_LinkedCellMLSolverExists(quasiNewtonSolver,cellMLSolver,err,error,*)

    !Argument variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver !<A pointer to the quasi Newton nonlinear solver to check the existstance of the linked CellML solver for
    TYPE(SolverType), POINTER :: cellMLSolver !<On exit, a pointer to the linked CellML solver for the specified quasi Newton nonlinear solver if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinearQuasiNewton_LinkedCellMLSolverExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLSolver)) CALL FlagError("CellML solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(quasiNewtonSolver)) CALL FlagError("Quasi Newton nonlinear solver is not associated.",err,error,*999)
#endif    

    cellMLSolver=>quasiNewtonSolver%cellMLEvaluatorSolver
 
    EXITS("SolverNonlinearQuasiNewton_LinkedCellMLSolverExists")
    RETURN
999 NULLIFY(cellMLSolver)
998 ERRORS("SolverNonlinearQuasiNewton_LinkedCellMLSolverExists",err,error)
    EXITS("SolverNonlinearQuasiNewton_LinkedCellMLSolverExists")
    RETURN 1

  END SUBROUTINE SolverNonlinearQuasiNewton_LinkedCellMLSolverExists
     
  !
  !================================================================================================================================
  !

  !>Gets the linked CellML solver for a quasi Newton solver. 
  SUBROUTINE SolverNonlinearQuasiNewton_LinkedCellMLSolverGet(quasiNewtonSolver,cellMLSolver,err,error,*)

    !Argument variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver !<A pointer to the quasi Newton nonlinear solver to get the linked CellML solver for
    TYPE(SolverType), POINTER :: cellMLSolver !<On exit, a pointer to the linked CellML solver for the specified quasi Newton nonlinear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinearQuasiNewton_LinkedCellMLSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLSolver)) CALL FlagError("CellML solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(quasiNewtonSolver)) CALL FlagError("Quasi Newton nonlinear solver is not associated.",err,error,*999)
#endif    

    cellMLSolver=>quasiNewtonSolver%cellMLEvaluatorSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(cellMLSolver)) &
      & CALL FlagError("Quasi Newton nonlinear solver linked CellML solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverNonlinearQuasiNewton_LinkedCellMLSolverGet")
    RETURN
999 NULLIFY(cellMLSolver)
998 ERRORS("SolverNonlinearQuasiNewton_LinkedCellMLSolverGet",err,error)
    EXITS("SolverNonlinearQuasiNewton_LinkedCellMLSolverGet")
    RETURN 1

  END SUBROUTINE SolverNonlinearQuasiNewton_LinkedCellMLSolverGet
     
  !
  !================================================================================================================================
  !

  !>Gets the linked linear solver for a quasi Newton solver. 
  SUBROUTINE SolverNonlinearQuasiNewton_LinkedLinearSolverGet(quasiNewtonSolver,linearSolver,err,error,*)

    !Argument variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver !<A pointer to the quasi Newton nonlinear solver to get the linked linear solver for
    TYPE(SolverType), POINTER :: linearSolver !<On exit, a pointer to the linked linear solver for the specified quasi Newton nonlinear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinearQuasiNewton_LinkedLinearSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(quasiNewtonSolver)) CALL FlagError("Quasi Newton nonlinear solver is not associated.",err,error,*999)
#endif    

    linearSolver=>quasiNewtonSolver%linearSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(linearSolver)) &
      & CALL FlagError("Quasi Newton nonlinear solver linked linear solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverNonlinearQuasiNewton_LinkedLinearSolverGet")
    RETURN
999 NULLIFY(linearSolver)
998 ERRORS("SolverNonlinearQuasiNewton_LinkedLinearSolverGet",err,error)
    EXITS("SolverNonlinearQuasiNewton_LinkedLinearSolverGet")
    RETURN 1

  END SUBROUTINE SolverNonlinearQuasiNewton_LinkedLinearSolverGet
     
  !
  !================================================================================================================================
  !

  !>Gets the nonlinear solver for a quasi Newton solver. 
  SUBROUTINE SolverNonlinearQuasiNewton_NonlinearSolverGet(quasiNewtonSolver,nonlinearSolver,err,error,*)

    !Argument variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver !<A pointer to the quasi Newton nonlinear solver to get the nonlinear solver for
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver !<On exit, a pointer to the nonlinear solver for the specified quasi Newton nonlinear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinearQuasiNewton_NonlinearSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(quasiNewtonSolver)) CALL FlagError("Quasi Newton nonlinear solver is not associated.",err,error,*999)
#endif    

    nonlinearSolver=>quasiNewtonSolver%nonlinearSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(nonlinearSolver)) &
      & CALL FlagError("Quasi Newton nonlinear solver nonlinear solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverNonlinearQuasiNewton_NonlinearSolverGet")
    RETURN
999 NULLIFY(nonlinearSolver)
998 ERRORS("SolverNonlinearQuasiNewton_NonlinearSolverGet",err,error)
    EXITS("SolverNonlinearQuasiNewton_NonlinearSolverGet")
    RETURN 1

  END SUBROUTINE SolverNonlinearQuasiNewton_NonlinearSolverGet
     
  !
  !================================================================================================================================
  !

  !>Returns the solve type of a quasi Newton nonlinear solver.
  SUBROUTINE SolverNonlinearQuasiNewton_SolverTypeGet(quasiNewtonNonlinearSolver,quasiNewtonNonlinearSolveType,err,error,*)

    !Argument variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonNonlinearSolver !<A pointer to the quasi Newton nonlinear solver to get the solve type for
    INTEGER(INTG), INTENT(OUT) :: quasiNewtonNonlinearSolveType !<On return, the quasi Newton nonlinear solver type. \see SolverRoutines_QuasiNewtonSolverTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinearQuasiNewton_SolverTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(quasiNewtonNonlinearSolver)) &
      & CALL FlagError("Quasi Newton nonlinear solver is not associated.",err,error,*999)
#endif    

    quasiNewtonNonlinearSolveType=quasiNewtonNonlinearSolver%quasiNewtonSolveType
    
    EXITS("SolverNonlinearQuasiNewton_SolverTypeGet")
    RETURN
999 ERRORSEXITS("SolverNonlinearQuasiNewton_SolverTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverNonlinearQuasiNewton_SolverTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the trustregion solver for a quasi Newton solver. 
  SUBROUTINE SolverNonlinearQuasiNewton_TrustregionSolverGet(quasiNewtonSolver,trustregionSolver,err,error,*)

    !Argument variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver !<A pointer to the quasi Newton nonlinear solver to get the trustregion solver for
    TYPE(QuasiNewtonTrustregionSolverType), POINTER :: trustregionSolver !<On exit, a pointer to the trustregion solver for the specified quasi Newton nonlinear solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinearQuasiNewton_TrustregionSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(trustregionSolver)) CALL FlagError("Trustregion solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(quasiNewtonSolver)) CALL FlagError("Quasi Newton nonlinear solver is not associated.",err,error,*999)
#endif    

    trustregionSolver=>quasiNewtonSolver%trustregionSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(trustregionSolver)) &
      & CALL FlagError("Quasi Newton nonlinear solver trustregion solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverNonlinearQuasiNewton_TrustregionSolverGet")
    RETURN
999 NULLIFY(trustregionSolver)
998 ERRORS("SolverNonlinearQuasiNewton_TrustregionSolverGet",err,error)
    EXITS("SolverNonlinearQuasiNewton_TrustregionSolverGet")
    RETURN 1

  END SUBROUTINE SolverNonlinearQuasiNewton_TrustregionSolverGet
     
  !
  !================================================================================================================================
  !

  !>Gets the quasi Newton solver for a quasi Newton linesearch solver. 
  SUBROUTINE SolverNonlinearQuasiNewtonLinesearch_QuasiNewtonSolverGet(linesearchSolver,quasiNewtonSolver,err,error,*)

    !Argument variables
    TYPE(QuasiNewtonLinesearchSolverType), POINTER :: linesearchSolver !<A pointer to the quasi Newton nonlinear linesearch solver to get the quasi Newton solver for
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver !<On exit, a pointer to the quasi Newton solver for the specified quasi Newton nonlinear linesearch solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinearQuasiNewtonLinesearch_QuasiNewtonSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(quasiNewtonSolver)) CALL FlagError("Quasi Newton solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(linesearchSolver)) CALL FlagError("Linesearch solver is not associated.",err,error,*999)
#endif    

    quasiNewtonSolver=>linesearchSolver%quasiNewtonSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(quasiNewtonSolver)) &
      & CALL FlagError("Quasi Newton nonlinear linesearch solver quasi Newton solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverNonlinearQuasiNewtonLinesearch_QuasiNewtonSolverGet")
    RETURN
999 NULLIFY(quasiNewtonSolver)
998 ERRORS("SolverNonlinearQuasiNewtonLinesearch_QuasiNewtonSolverGet",err,error)
    EXITS("SolverNonlinearQuasiNewtonLinesearch_QuasiNewtonSolverGet")
    RETURN 1

  END SUBROUTINE SolverNonlinearQuasiNewtonLinesearch_QuasiNewtonSolverGet
     
  !
  !================================================================================================================================
  !

  !>Gets the quasi Newton solver for a quasi Newton trustregion solver. 
  SUBROUTINE SolverNonlinearQuasiNewtonTrustregion_QuasiNewtonSolverGet(trustregionSolver,quasiNewtonSolver,err,error,*)

    !Argument variables
    TYPE(QuasiNewtonTrustregionSolverType), POINTER :: trustregionSolver !<A pointer to the quasi Newton nonlinear trustregion solver to get the quasi Newton solver for
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver !<On exit, a pointer to the quasi Newton solver for the specified quasi Newton nonlinear trustregion solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinearQuasiNewtonTrustregion_QuasiNewtonSolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(quasiNewtonSolver)) CALL FlagError("Quasi Newton solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(trustregionSolver)) CALL FlagError("Trustregion solver is not associated.",err,error,*999)
#endif    

    quasiNewtonSolver=>trustregionSolver%quasiNewtonSolver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(quasiNewtonSolver)) &
      & CALL FlagError("Quasi Newton nonlinear trustregion solver quasi Newton solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverNonlinearQuasiNewtonTrustregion_QuasiNewtonSolverGet")
    RETURN
999 NULLIFY(quasiNewtonSolver)
998 ERRORS("SolverNonlinearQuasiNewtonTrustregion_QuasiNewtonSolverGet",err,error)
    EXITS("SolverNonlinearQuasiNewtonTrustregion_QuasiNewtonSolverGet")
    RETURN 1

  END SUBROUTINE SolverNonlinearQuasiNewtonTrustregion_QuasiNewtonSolverGet
     
  !
  !================================================================================================================================
  !

  !>Returns the type of certainty for an optimiser solver.
  SUBROUTINE SolverOptimiser_CertaintyTypeGet(optimiserSolver,solverCertaintyType,err,error,*)

    !Argument variables
    TYPE(OptimiserSolverType), POINTER :: optimiserSolver !<A pointer the optimiser solver to get the certainty type for.
    INTEGER(INTG), INTENT(OUT) :: solverCertaintyType !<On exit, the type of certainty for the optimiser solver \see SolverRoutines_OptimiserCertaintyTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverOptimiser_CertaintyTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(optimiserSolver)) CALL FlagError("Optimiser solver is not associated.",err,error,*999)
#endif    
    
    solverCertaintyType=optimiserSolver%certaintyType
    
    EXITS("SolverOptimiser_CertaintyTypeGet")
    RETURN
999 ERRORSEXITS("SolverOptimiser_CertaintyTypeGet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverOptimiser_CertaintyTypeGet

  !
  !================================================================================================================================
  !

  !>Returns the type of constraints for an optimiser solver.
  SUBROUTINE SolverOptimiser_ConstraintTypeGet(optimiserSolver,solverConstraintType,err,error,*)

    !Argument variables
    TYPE(OptimiserSolverType), POINTER :: optimiserSolver !<A pointer the optimiser solver to get the constraint type for.
    INTEGER(INTG), INTENT(OUT) :: solverConstraintType !<On exit, the type of constraint for the optimiser solver \see SolverRoutines_OptimiserConstraintTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverOptimiser_ConstraintTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(optimiserSolver)) CALL FlagError("Optimiser solver is not associated.",err,error,*999)
#endif    
    
    solverConstraintType=optimiserSolver%constraintType
    
    EXITS("SolverOptimiser_ConstraintTypeGet")
    RETURN
999 ERRORSEXITS("SolverOptimiser_ConstraintTypeGet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverOptimiser_ConstraintTypeGet

  !
  !================================================================================================================================
  !

  !>Returns the type of library to use for an optimiser solver.
  SUBROUTINE SolverOptimiser_LibraryTypeGet(optimiserSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(OptimiserSolverType), POINTER :: optimiserSolver !<A pointer the optimiser solver to get the library type for.
    INTEGER(INTG), INTENT(OUT) :: solverLibraryType !<On exit, the type of library used for the optimiser solver \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverOptimiser_LibraryTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(optimiserSolver)) CALL FlagError("Optimiser solver is not associated.",err,error,*999)
#endif    
    
    solverLibraryType=optimiserSolver%solverLibrary
    
    EXITS("SolverOptimiser_LibraryTypeGet")
    RETURN
999 ERRORSEXITS("SolverOptimiser_LibraryTypeGet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverOptimiser_LibraryTypeGet

  !
  !================================================================================================================================
  !

  !>Returns the type of objective for an optimiser solver.
  SUBROUTINE SolverOptimiser_ObjectiveTypeGet(optimiserSolver,solverObjectiveType,err,error,*)

    !Argument variables
    TYPE(OptimiserSolverType), POINTER :: optimiserSolver !<A pointer the optimiser solver to get the objective type for.
    INTEGER(INTG), INTENT(OUT) :: solverObjectiveType !<On exit, the type of objective for the optimiser solver \see SolverRoutines_OptimiserObjectiveTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverOptimiser_ObjectiveTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(optimiserSolver)) CALL FlagError("Optimiser solver is not associated.",err,error,*999)
#endif    
    
    solverObjectiveType=optimiserSolver%objectiveType
    
    EXITS("SolverOptimiser_ObjectiveTypeGet")
    RETURN
999 ERRORSEXITS("SolverOptimiser_ObjectiveTypeGet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverOptimiser_ObjectiveTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the solver for a optimiser solver. 
  SUBROUTINE SolverOptimiser_SolverGet(optimiserSolver,solver,err,error,*)

    !Argument variables
    TYPE(OptimiserSolverType), POINTER :: optimiserSolver !<A pointer to the optimiser solver to get the solver for
    TYPE(SolverType), POINTER :: solver !<On exit, a pointer to the solver for the specified optimiser solver. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverOptimiser_SolverGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solver)) CALL FlagError("Solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(optimiserSolver)) CALL FlagError("Optimiser solver is not associated.",err,error,*999)
#endif    

    solver=>optimiserSolver%solver

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Optimiser solver solver is not associated.",err,error,*999)
#endif    
 
    EXITS("SolverOptimiser_SolverGet")
    RETURN
999 NULLIFY(solver)
998 ERRORSEXITS("SolverOptimiser_SolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverOptimiser_SolverGet
     
  !
  !================================================================================================================================
  !

  !>Returns the type of variable for an optimiser solver.
  SUBROUTINE SolverOptimiser_VariableTypeGet(optimiserSolver,solverVariableType,err,error,*)

    !Argument variables
    TYPE(OptimiserSolverType), POINTER :: optimiserSolver !<A pointer the optimiser solver to get the variable type for.
    INTEGER(INTG), INTENT(OUT) :: solverVariableType !<On exit, the type of variable for the optimiser solver \see SolverRoutines_OptimiserVariableTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverOptimiser_VariableTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(optimiserSolver)) CALL FlagError("Optimiser solver is not associated.",err,error,*999)
#endif    
    
    solverVariableType=optimiserSolver%variableType
    
    EXITS("SolverOptimiser_VariableTypeGet")
    RETURN
999 ERRORSEXITS("SolverOptimiser_VariableTypeGet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverOptimiser_VariableTypeGet

  !
  !================================================================================================================================
  !

END MODULE SolverAccessRoutines
