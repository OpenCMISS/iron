!> \file
!> \author Chris Bradley
!> \brief This module contains all field access method routines.
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

!> This module contains all field access method routines.
MODULE FieldAccessRoutines
  
  USE BaseRoutines
  USE DecompositionAccessRoutines
  USE Kinds
  USE ISO_VARYING_STRING
  USE Strings
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup FieldRoutines_DependentTypes FieldRoutines::DependentTypes
  !> \brief Depedent field parameter types
  !> \see FieldRoutines,OPENCMISS_FieldDependentTypes
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_INDEPENDENT_TYPE=1 !<Independent field type \see FieldRoutines_DependentTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DEPENDENT_TYPE=2 !<Dependent field type \see FieldRoutines_DependentTypes,FieldRoutines
  !>@}

  !> \addtogroup FieldRoutines_DimensionTypes FieldRoutines::DimensionTypes
  !> \brief Field dimension parameter types
  !> \see FieldRoutines,OPENCMISS_FieldDimensionTypes
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_SCALAR_DIMENSION_TYPE=1 !<Scalar field \see FieldRoutines_DimensionTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_VECTOR_DIMENSION_TYPE=2 !<Vector field \see FieldRoutines_DimensionTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_TENSOR_DIMENSION_TYPE=3 !<Tensor field \see FieldRoutines_DimensionTypes,FieldRoutines
  !>@}

  !> \addtogroup FieldRoutines_FieldTypes FieldRoutines::FieldTypes
  !> \brief Field type parameters
  !> \see FieldRoutines,OPENCMISS_FieldTypes
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_GEOMETRIC_TYPE=1 !<Geometric field \see FieldRoutines_FieldTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_FIBRE_TYPE=2 !<Fibre field \see FieldRoutines_FieldTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_GENERAL_TYPE=3 !<General field \see FieldRoutines_FieldTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_MATERIAL_TYPE=4 !<Material field \see FieldRoutines_FieldTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_GEOMETRIC_GENERAL_TYPE=5 !<Geometric general field \see FieldRoutines_FieldTypes,FieldRoutines
  !>@}

  !> \addtogroup FieldRoutines_DofTypes FieldRoutines::DofTypes
  !> \brief Field dof type parameters
  !> \see FieldRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_CONSTANT_DOF_TYPE=1 !<The dof is from a field variable component with constant interpolation \see FieldRoutines_DofTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_ELEMENT_DOF_TYPE=2 !<The dof is from a field variable component with element based interpolation \see FieldRoutines_DofTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_NODE_DOF_TYPE=3 !<The dof is from a field variable component with node based interpolation \see FieldRoutines_DofTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_GRID_POINT_DOF_TYPE=4 !<The dof is from a field variable component with grid point based interpolation \see FieldRoutines_DofTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_GAUSS_POINT_DOF_TYPE=5 !<The dof is from a field variable component with Gauss point based interpolation \see FieldRoutines_DofTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DATA_POINT_DOF_TYPE=6 !<The dof is from a field variable component with Gauss point based interpolation \see FieldRoutines_DofTypes,FieldRoutines
  !>@}
  !> \addtogroup FieldRoutines_DOFOrderTypes FieldRoutines::DOFOrderTypes
  !> \brief Field DOF order types
  !> \see FieldRoutines,OPENCMISS_FieldDOFOrderTypes
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_SEPARATED_COMPONENT_DOF_ORDER=1 !<Field variable component dofs are not contiguous \see FieldRoutines_DOFOrderTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_CONTIGUOUS_COMPONENT_DOF_ORDER=2 !<Field variable component dofs are contiguous \see FieldRoutines_DOFOrderTypes,FieldRoutines
  !>@}

  !> \addtogroup FieldRoutines_ScalingTypes FieldRoutines::ScalingTypes
  !> \brief Field scaling type parameters
  !> \see FieldRoutines,OPENCMISS_FieldScalingTypes
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_NO_SCALING=0 !<The field is not scaled \see FieldRoutines_ScalingTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_UNIT_SCALING=1 !<The field has unit scaling \see FieldRoutines_ScalingTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_ARC_LENGTH_SCALING=2 !<The field has arc length scaling \see FieldRoutines_ScalingTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_ARITHMETIC_MEAN_SCALING=3 !<The field has arithmetic mean of the arc length scaling \see FieldRoutines_ScalingTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_GEOMETRIC_MEAN_SCALING=4 !<The field has geometric mean of the arc length scaling \see FieldRoutines_ScalingTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_HARMONIC_MEAN_SCALING=5 !<The field has harmonic mean of the arc length scaling \see FieldRoutines_ScalingTypes,FieldRoutines
  !>@}

  !> \addtogroup FieldRoutines_InterpolationComponentsTypes FieldRoutines::InterpolationComponentsTypes
  !> \brief Field interpolation components types
  !> \see FieldRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_ALL_COMPONENTS_TYPE=1 !<The field is interpolated for all components \see FieldRoutines_InterpolationComponentsTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_GEOMETRIC_COMPONENTS_TYPE=2 !<The field is interpolated for geometric components \see FieldRoutines_InterpolationComponentsTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_NONGEOMETRIC_COMPONENTS_TYPE=3 !<The field is interpolated for non-geometric components \see FieldRoutines_InterpolationComponentsTypes,FieldRoutines
  !>@}
  
  !> \addtogroup FieldRoutines_VariableTypes FieldRoutines::VariableTypes
  !> \brief Field variable type parameters
  !> \see FieldRoutines,OPENCMISS_FieldVariableTypes
  !> \todo sort out variable access routines so that you are always accessing by variable type rather than variable number.
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_NUMBER_OF_VARIABLE_TYPES=49 !<Number of different field variable types possible \see FieldRoutines_VariableTypes,FieldRoutines 
  INTEGER(INTG), PARAMETER :: FIELD_NUMBER_OF_VARIABLE_SUBTYPES=4 !<The number of variants of a particular variable - currently 4. U, delUdelN, delUdelT,del2UdelT2 \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_U_VARIABLE_TYPE=1 !<Standard variable type i.e., u \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELUDELN_VARIABLE_TYPE=2 !<Normal derivative variable type i.e., du/dn \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELUDELT_VARIABLE_TYPE=3 !<First time derivative variable type i.e., du/dt \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DEL2UDELT2_VARIABLE_TYPE=4 !<Second time derivative variable type i.e., d^2u/dt^2 \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_V_VARIABLE_TYPE=5 !<Second standard variable type i.e., v \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELVDELN_VARIABLE_TYPE=6 !<Second normal derivative variable type i.e., dv/dn \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELVDELT_VARIABLE_TYPE=7 !<First time derivative variable type i.e., du/dt \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DEL2VDELT2_VARIABLE_TYPE=8 !<Second time derivative variable type i.e., d^2u/dt^2 \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_W_VARIABLE_TYPE=9 !<Third standard variable type i.e., w \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_U1_VARIABLE_TYPE=10 !<Third standard variable type i.e., v \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELU1DELN_VARIABLE_TYPE=11 !<Third normal derivative variable type i.e., dv/dn \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELU1DELT_VARIABLE_TYPE=12 !<Third time derivative variable type i.e., du/dt \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DEL2U1DELT2_VARIABLE_TYPE=13 !<Third time derivative variable type i.e., d^2u/dt^2 \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_U2_VARIABLE_TYPE=14 !<Fourth standard variable type i.e., v \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELU2DELN_VARIABLE_TYPE=15 !<Fourth normal derivative variable type i.e., dv/dn \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELU2DELT_VARIABLE_TYPE=16 !<Fourth time derivative variable type i.e., du/dt \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DEL2U2DELT2_VARIABLE_TYPE=17 !<Fourth time derivative variable type i.e., d^2u/dt^2 \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_U3_VARIABLE_TYPE=18 !<Fifth standard variable type i.e., v \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELU3DELN_VARIABLE_TYPE=19 !<Fifth normal derivative variable type i.e., dv/dn \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELU3DELT_VARIABLE_TYPE=20 !<Fifth time derivative variable type i.e., du/dt \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DEL2U3DELT2_VARIABLE_TYPE=21 !<Fifth time derivative variable type i.e., d^2u/dt^2 \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_U4_VARIABLE_TYPE=22 !<Sixth standard variable type i.e., v \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELU4DELN_VARIABLE_TYPE=23 !<Sixth normal derivative variable type i.e., dv/dn \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELU4DELT_VARIABLE_TYPE=24 !<Sixth time derivative variable type i.e., du/dt \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DEL2U4DELT2_VARIABLE_TYPE=25 !<Sixth time derivative variable type i.e., d^2u/dt^2 \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_U5_VARIABLE_TYPE=26 !<Seventh standard variable type i.e., v \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELU5DELN_VARIABLE_TYPE=27 !<Seventh normal derivative variable type i.e., dv/dn \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELU5DELT_VARIABLE_TYPE=28 !<Seventh time derivative variable type i.e., du/dt \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DEL2U5DELT2_VARIABLE_TYPE=29 !<Seventh time derivative variable type i.e., d^2u/dt^2 \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_U6_VARIABLE_TYPE=30 !<Eighth standard variable type i.e., v \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELU6DELN_VARIABLE_TYPE=31 !<Eighth normal derivative variable type i.e., dv/dn \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELU6DELT_VARIABLE_TYPE=32 !<Eighth time derivative variable type i.e., du/dt \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DEL2U6DELT2_VARIABLE_TYPE=33 !<Eighth time derivative variable type i.e., d^2u/dt^2 \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_U7_VARIABLE_TYPE=34 !<Ninth standard variable type i.e., v \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELU7DELN_VARIABLE_TYPE=35 !<Ninth normal derivative variable type i.e., dv/dn \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELU7DELT_VARIABLE_TYPE=36 !<Ninth time derivative variable type i.e., du/dt \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DEL2U7DELT2_VARIABLE_TYPE=37 !<Ninth time derivative variable type i.e., d^2u/dt^2 \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_U8_VARIABLE_TYPE=38 !<Tenth standard variable type i.e., v \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELU8DELN_VARIABLE_TYPE=39 !<Tenth normal derivative variable type i.e., dv/dn \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELU8DELT_VARIABLE_TYPE=40 !<Tenth time derivative variable type i.e., du/dt \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DEL2U8DELT2_VARIABLE_TYPE=41 !<Tenth time derivative variable type i.e., d^2u/dt^2 \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_U9_VARIABLE_TYPE=42 !<Eleventh standard variable type i.e., v \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELU9DELN_VARIABLE_TYPE=43 !<Eleventh normal derivative variable type i.e., dv/dn \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELU9DELT_VARIABLE_TYPE=44 !<Eleventh time derivative variable type i.e., du/dt \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DEL2U9DELT2_VARIABLE_TYPE=45 !<Eleventh time derivative variable type i.e., d^2u/dt^2 \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_U10_VARIABLE_TYPE=46 !<Twelfth standard variable type i.e., v \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELU10DELN_VARIABLE_TYPE=47 !<Twelfth normal derivative variable type i.e., dv/dn \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DELU10DELT_VARIABLE_TYPE=48 !<Twelfth time derivative variable type i.e., du/dt \see FieldRoutines_VariableTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DEL2U10DELT2_VARIABLE_TYPE=49 !<Twelfth time derivative variable type i.e., d^2u/dt^2 \see FieldRoutines_VariableTypes,FieldRoutines
  !>@}
  
  !> \addtogroup FieldRoutines_ParameterSetTypes FieldRoutines::ParameterSetTypes
  !> \brief Field parameter set type parameters \todo make program defined constants negative?
  !> \see FieldRoutines,OPENCMISS_FieldParameterSetTypes
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_NUMBER_OF_SET_TYPES=99 !<The maximum number of different parameter sets for a field \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_VALUES_SET_TYPE=1 !<The parameter set corresponding to the field values (at time T+DT for dynamic problems) \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_BOUNDARY_CONDITIONS_SET_TYPE=2 !<The parameter set corresponding to the field boundary conditions \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_INITIAL_VALUES_SET_TYPE=3 !<The parameter set corresponding to the field initial values \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_INCREMENTAL_VALUES_SET_TYPE=4 !<The parameter set corresponding to the field incremental values \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_ANALYTIC_VALUES_SET_TYPE=5 !<The parameter set corresponding to the analytic field values \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS_VALUES_SET_TYPE=6 !<The parameter set corresponding to the previous field values (at time T) \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS2_VALUES_SET_TYPE=7 !<The parameter set corresponding to the previous field values (at time T-DT) \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS3_VALUES_SET_TYPE=8 !<The parameter set corresponding to the previous field values (at time T-DT) \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_NEXT_VALUES_SET_TYPE=9 !<The parameter set corresponding to the next field values (at time T+dT) \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_MEAN_PREDICTED_DISPLACEMENT_SET_TYPE=10 !<The parameter set corresponding to the mean predicited values (at time T+DT) \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_VELOCITY_VALUES_SET_TYPE=11 !<The parameter set corresponding to the velocity values (at time T+DT) \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_INITIAL_VELOCITY_SET_TYPE=12 !<The parameter set corresponding to the initial velocity values for dynamic problems. This is also the previous velocity values \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS_VELOCITY_SET_TYPE=12 !<The parameter set corresponding to the previous velocity values (at time T). This is also the initial velocity values for dynamic problems. \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_MEAN_PREDICTED_VELOCITY_SET_TYPE=13 !<The parameter set corresponding to the mean predicited velocity values (at time T+DT) \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_ACCELERATION_VALUES_SET_TYPE=14 !<The parameter set corresponding to the acceleration values (at time T+DT) \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_INITIAL_ACCELERATION_SET_TYPE=15 !<The parameter set corresponding to the initial acceleration values for dynamic problems. This is also the previous accelearation values \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS_ACCELERATION_SET_TYPE=15 !<The parameter set corresponding to the previous acceleration values (at time T).This is also the initial acceleration values for dynamic problems. \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_MEAN_PREDICTED_ACCELERATION_SET_TYPE=16 !<The parameter set corresponding to the mean predicted acceleration values (at time T+DT) \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_PREDICTED_DISPLACEMENT_SET_TYPE=17 !<The parameter set corresponding to the predicted values (at time T+DT) \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_PREDICTED_VELOCITY_SET_TYPE=18 !<The parameter set corresponding to the predicted velocity values (at time T+DT) \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_PREDICTED_ACCELERATION_SET_TYPE=19 !<The parameter set corresponding to the predicted acceleration values (at time T+DT) \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_RESIDUAL_SET_TYPE=20 !<The parameter set corresponding to the evaluated residual values (at time T+DT) \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS_RESIDUAL_SET_TYPE=21 !<The parameter set corresponding to the residual values evaluated previously (at time T) \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS2_RESIDUAL_SET_TYPE=22 !<The parameter set corresponding to the residual values evaluated previously (at time T-DT) \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS3_RESIDUAL_SET_TYPE=23 !<The parameter set corresponding to the residual values evaluated previously (at time T-2*DT) \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_MESH_DISPLACEMENT_SET_TYPE=24 !<The parameter set corresponding to the mesh displacement values for ALE \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_MESH_VELOCITY_SET_TYPE=25 !<The parameter set corresponding to the mesh velocity values for ALE \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_BOUNDARY_SET_TYPE=26 !<The parameter set corresponding to the mesh velocity values for ALE \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_INPUT_DATA1_SET_TYPE=27 !<The parameter set corresponding to a input field \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_INPUT_DATA2_SET_TYPE=28 !<The parameter set corresponding to a input field \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_INPUT_DATA3_SET_TYPE=29 !<The parameter set corresponding to a input field \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_INPUT_VEL1_SET_TYPE=30 !<The parameter set corresponding to a input field (PPE)\see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_INPUT_VEL2_SET_TYPE=31 !<The parameter set corresponding to a input field (PPE)\see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_INPUT_VEL3_SET_TYPE=32 !<The parameter set corresponding to a input field (PPE)\see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_INPUT_LABEL_SET_TYPE=33 !<The parameter set corresponding to a input field (PPE)\see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_PRESSURE_VALUES_SET_TYPE=34 !<The parameter set corresponding to the surface pressure values (at time T+DT). \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS_PRESSURE_SET_TYPE=35 !<The parameter set corresponding to the previous surface pressure values (at previous increment step). \see FieldRoutines_ParameterSetTypes,FieldRoutines 
  INTEGER(INTG), PARAMETER :: FIELD_RELATIVE_VELOCITY_SET_TYPE=36 !<The parameter set corresponding to the relative velocity values for ALE \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_NEGATIVE_MESH_VELOCITY_SET_TYPE=37 !<The parameter set corresponding to the NEGATIVE mesh velocity values for ALE \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE=38 !<The parameter set corresponding to the previous iteration field values (at iteration n) \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE=39 !<The parameter set corresponding to the impermeable flag field values \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_INTEGRATED_NEUMANN_SET_TYPE=40 !<Stores integrated Neumann values calculated from Neumann point values \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_UPWIND_VALUES_SET_TYPE=41 !<Stores upwind values associated with a field. \see FieldRoutines_ParameterSetTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS_UPWIND_VALUES_SET_TYPE=42 !<Stores upwind values associated with a field from previous timestep. \see FieldRoutines_ParameterSetTypes,FieldRoutines
  !>@}

  !> \addtogroup FieldRoutines_InterpolationTypes FieldRoutines::InterpolationTypes
  !> \brief Field interpolation parameters
  !> \see FieldRoutines,OPENCMISS_FieldInterpolationTypes
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_CONSTANT_INTERPOLATION=1 !<Constant interpolation. One parameter for the field \see FieldRoutines_InterpolationTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_ELEMENT_BASED_INTERPOLATION=2 !<Element based interpolation. Parameters are different in each element \see FieldRoutines_InterpolationTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_NODE_BASED_INTERPOLATION=3 !<Node based interpolation. Parameters are nodal based and a basis function is used \see FieldRoutines_InterpolationTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_GRID_POINT_BASED_INTERPOLATION=4 !<Grid point based interpolation. Parameters are different at each grid point \see FieldRoutines_InterpolationTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_GAUSS_POINT_BASED_INTERPOLATION=5 !<Gauss point based interpolation. Parameters are different at each Gauss point \see FieldRoutines_InterpolationTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DATA_POINT_BASED_INTERPOLATION=6 !<data point based interpolation. Parameters are different at each data point \see FieldRoutines_InterpolationTypes,FieldRoutines
  !>@}
   
  !> \addtogroup FieldRoutines_DataTypes FieldRoutines::DataTypes
  !> \brief Field data types
  !> \see FieldRoutines,OPENCMISS_FieldDataTypes
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_INTG_TYPE=1 !<Integer field data type \see FieldRoutines_DataTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_SP_TYPE=2 !<Single precision real field data type \see FieldRoutines_DataTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_DP_TYPE=3 !<Double precision real field data type \see FieldRoutines_DataTypes,FieldRoutines
  INTEGER(INTG), PARAMETER :: FIELD_L_TYPE=4 !<Logical field data type \see FieldRoutines_DataTypes,FieldRoutines
  !>@}
   
  
  !Module types

  !Module variables
  
  !Interfaces

  !>Gets the label for a field variable component.
  INTERFACE Field_ComponentLabelGet
    MODULE PROCEDURE Field_ComponentLabelGetC
    MODULE PROCEDURE Field_ComponentLabelGetVS
  END INTERFACE Field_ComponentLabelGet
  
  !>Gets the label for a field.
  INTERFACE Field_LabelGet
    MODULE PROCEDURE Field_LabelGetC
    MODULE PROCEDURE Field_LabelGetVS
  END INTERFACE Field_LabelGet
 
  !>Finds and returns a field identified by a user number. If no field  with that number exists field is left nullified.
  INTERFACE Field_UserNumberFind
    MODULE PROCEDURE Field_UserNumberFindInterface
    MODULE PROCEDURE Field_UserNumberFindRegion
  END INTERFACE Field_UserNumberFind
  
  !>Gets the label for a field variable.
  INTERFACE Field_VariableLabelGet
    MODULE PROCEDURE Field_VariableLabelGetC
    MODULE PROCEDURE Field_VariableLabelGetVS
  END INTERFACE Field_VariableLabelGet
 
  !>Gets the label for a field variable.
  INTERFACE FieldVariable_LabelGet
    MODULE PROCEDURE FieldVariable_LabelGetC
    MODULE PROCEDURE FieldVariable_LabelGetVS
  END INTERFACE FieldVariable_LabelGet
  
   !>Gets the label for a field variable component.
  INTERFACE FieldVariable_ComponentLabelGet
    MODULE PROCEDURE FieldVariable_ComponentLabelGetC
    MODULE PROCEDURE FieldVariable_ComponentLabelGetVS
  END INTERFACE FieldVariable_ComponentLabelGet
  
  PUBLIC FIELD_INDEPENDENT_TYPE,FIELD_DEPENDENT_TYPE

  PUBLIC FIELD_SCALAR_DIMENSION_TYPE,FIELD_VECTOR_DIMENSION_TYPE,FIELD_TENSOR_DIMENSION_TYPE

  PUBLIC FIELD_GEOMETRIC_TYPE,FIELD_FIBRE_TYPE,FIELD_GENERAL_TYPE,FIELD_MATERIAL_TYPE,FIELD_GEOMETRIC_GENERAL_TYPE
  PUBLIC FIELD_CONSTANT_DOF_TYPE,FIELD_ELEMENT_DOF_TYPE,FIELD_NODE_DOF_TYPE,FIELD_GRID_POINT_DOF_TYPE,FIELD_GAUSS_POINT_DOF_TYPE, &
    & FIELD_DATA_POINT_DOF_TYPE

  PUBLIC FIELD_SEPARATED_COMPONENT_DOF_ORDER,FIELD_CONTIGUOUS_COMPONENT_DOF_ORDER

  PUBLIC FIELD_NO_SCALING,FIELD_UNIT_SCALING,FIELD_ARC_LENGTH_SCALING,FIELD_HARMONIC_MEAN_SCALING,FIELD_ARITHMETIC_MEAN_SCALING, &
    & FIELD_GEOMETRIC_MEAN_SCALING

  PUBLIC FIELD_ALL_COMPONENTS_TYPE,FIELD_GEOMETRIC_COMPONENTS_TYPE,FIELD_NONGEOMETRIC_COMPONENTS_TYPE
  
  PUBLIC FIELD_NUMBER_OF_VARIABLE_TYPES,FIELD_NUMBER_OF_VARIABLE_SUBTYPES,FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
    & FIELD_DELUDELT_VARIABLE_TYPE,FIELD_DEL2UDELT2_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE, &
    & FIELD_DELVDELT_VARIABLE_TYPE,FIELD_DEL2VDELT2_VARIABLE_TYPE,&
    & FIELD_W_VARIABLE_TYPE, &
    & FIELD_U1_VARIABLE_TYPE,FIELD_DELU1DELN_VARIABLE_TYPE,FIELD_DELU1DELT_VARIABLE_TYPE,FIELD_DEL2U1DELT2_VARIABLE_TYPE,&
    & FIELD_U2_VARIABLE_TYPE,FIELD_DELU2DELN_VARIABLE_TYPE,FIELD_DELU2DELT_VARIABLE_TYPE,FIELD_DEL2U2DELT2_VARIABLE_TYPE,&
    & FIELD_U3_VARIABLE_TYPE,FIELD_DELU3DELN_VARIABLE_TYPE,FIELD_DELU3DELT_VARIABLE_TYPE,FIELD_DEL2U3DELT2_VARIABLE_TYPE,&
    & FIELD_U4_VARIABLE_TYPE,FIELD_DELU4DELN_VARIABLE_TYPE,FIELD_DELU4DELT_VARIABLE_TYPE,FIELD_DEL2U4DELT2_VARIABLE_TYPE,&
    & FIELD_U5_VARIABLE_TYPE,FIELD_DELU5DELN_VARIABLE_TYPE,FIELD_DELU5DELT_VARIABLE_TYPE,FIELD_DEL2U5DELT2_VARIABLE_TYPE,&
    & FIELD_U6_VARIABLE_TYPE,FIELD_DELU6DELN_VARIABLE_TYPE,FIELD_DELU6DELT_VARIABLE_TYPE,FIELD_DEL2U6DELT2_VARIABLE_TYPE,&
    & FIELD_U7_VARIABLE_TYPE,FIELD_DELU7DELN_VARIABLE_TYPE,FIELD_DELU7DELT_VARIABLE_TYPE,FIELD_DEL2U7DELT2_VARIABLE_TYPE,&
    & FIELD_U8_VARIABLE_TYPE,FIELD_DELU8DELN_VARIABLE_TYPE,FIELD_DELU8DELT_VARIABLE_TYPE,FIELD_DEL2U8DELT2_VARIABLE_TYPE,&
    & FIELD_U9_VARIABLE_TYPE,FIELD_DELU9DELN_VARIABLE_TYPE,FIELD_DELU9DELT_VARIABLE_TYPE,FIELD_DEL2U9DELT2_VARIABLE_TYPE,&
    & FIELD_U10_VARIABLE_TYPE,FIELD_DELU10DELN_VARIABLE_TYPE,FIELD_DELU10DELT_VARIABLE_TYPE,FIELD_DEL2U10DELT2_VARIABLE_TYPE

  PUBLIC FIELD_NUMBER_OF_SET_TYPES,FIELD_VALUES_SET_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,FIELD_INITIAL_VALUES_SET_TYPE, &
    & FIELD_INCREMENTAL_VALUES_SET_TYPE,FIELD_ANALYTIC_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE, &
    & FIELD_PREVIOUS2_VALUES_SET_TYPE,FIELD_PREVIOUS3_VALUES_SET_TYPE,FIELD_NEXT_VALUES_SET_TYPE, &
    & FIELD_MEAN_PREDICTED_DISPLACEMENT_SET_TYPE,FIELD_VELOCITY_VALUES_SET_TYPE,FIELD_INITIAL_VELOCITY_SET_TYPE, &
    & FIELD_PREVIOUS_VELOCITY_SET_TYPE,FIELD_MEAN_PREDICTED_VELOCITY_SET_TYPE,FIELD_ACCELERATION_VALUES_SET_TYPE, &
    & FIELD_INITIAL_ACCELERATION_SET_TYPE,FIELD_PREVIOUS_ACCELERATION_SET_TYPE,FIELD_MEAN_PREDICTED_ACCELERATION_SET_TYPE, &
    & FIELD_PREDICTED_DISPLACEMENT_SET_TYPE,FIELD_PREDICTED_VELOCITY_SET_TYPE,FIELD_PREDICTED_ACCELERATION_SET_TYPE, &
    & FIELD_RESIDUAL_SET_TYPE,FIELD_PREVIOUS_RESIDUAL_SET_TYPE,FIELD_PREVIOUS2_RESIDUAL_SET_TYPE, &
    & FIELD_PREVIOUS3_RESIDUAL_SET_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,FIELD_MESH_VELOCITY_SET_TYPE, &
    & FIELD_BOUNDARY_SET_TYPE,FIELD_INPUT_DATA1_SET_TYPE,FIELD_INPUT_DATA2_SET_TYPE,FIELD_INPUT_DATA3_SET_TYPE, &
    & FIELD_PRESSURE_VALUES_SET_TYPE,FIELD_PREVIOUS_PRESSURE_SET_TYPE,FIELD_RELATIVE_VELOCITY_SET_TYPE, &
    & FIELD_NEGATIVE_MESH_VELOCITY_SET_TYPE,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,FIELD_INPUT_VEL1_SET_TYPE, &
    & FIELD_INPUT_VEL2_SET_TYPE,FIELD_INPUT_VEL3_SET_TYPE,FIELD_INPUT_LABEL_SET_TYPE,FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE, &
    & FIELD_INTEGRATED_NEUMANN_SET_TYPE,FIELD_UPWIND_VALUES_SET_TYPE,FIELD_PREVIOUS_UPWIND_VALUES_SET_TYPE

  PUBLIC FIELD_CONSTANT_INTERPOLATION,FIELD_ELEMENT_BASED_INTERPOLATION,FIELD_NODE_BASED_INTERPOLATION, &
    & FIELD_GRID_POINT_BASED_INTERPOLATION,FIELD_GAUSS_POINT_BASED_INTERPOLATION,FIELD_DATA_POINT_BASED_INTERPOLATION

  PUBLIC FIELD_INTG_TYPE,FIELD_SP_TYPE,FIELD_DP_TYPE,FIELD_L_TYPE
  
  PUBLIC Field_AssertIsFinished,Field_AssertNotFinished

  PUBLIC Field_ComponentDOFGetConstant

  PUBLIC Field_ComponentDOFGetUserDataPoint

  PUBLIC Field_ComponentDOFGetUserElement

  PUBLIC Field_ComponentDOFGetUserNode

  PUBLIC Field_ComponentInterpolationCheck
  
  PUBLIC Field_ComponentInterpolationGet
 
  PUBLIC Field_ComponentLabelGet

  PUBLIC Field_ComponentMeshComponentCheck

  PUBLIC Field_ComponentMeshComponentGet

  PUBLIC Field_CoordinateSystemGet  

  PUBLIC Field_CreateValuesCacheGet

  PUBLIC Field_DataProjectionGet

  PUBLIC Field_DataTypeCheck

  PUBLIC Field_DataTypeGet

  PUBLIC Field_DecompositionGet

  PUBLIC Field_DependentTypeCheck

  PUBLIC Field_DependentTypeGet

  PUBLIC Field_DimensionCheck

  PUBLIC Field_DimensionGet

  PUBLIC Field_DOFOrderTypeCheck

  PUBLIC Field_DOFOrderTypeGet

  PUBLIC Field_FieldsGet

  PUBLIC Field_GeometricFieldGet

  PUBLIC Field_GeometricGeneralFieldGet

  PUBLIC Field_GeometricParametersGet

  PUBLIC Field_InterfaceGet

  PUBLIC Field_IsInterfaceField

  PUBLIC Field_IsRegionField

  PUBLIC Field_LabelGet

  PUBLIC Field_NumberOfComponentsCheck

  PUBLIC Field_NumberOfComponentsGet

  PUBLIC Field_NumberOfVariablesCheck

  PUBLIC Field_NumberOfVariablesGet

  PUBLIC Field_ParameterSetGet

  PUBLIC Field_RegionGet
  
  PUBLIC Field_ScaleFactorsVectorGet

  PUBLIC Field_ScalingTypeCheck

  PUBLIC Field_ScalingTypeGet

  PUBLIC Field_TypeCheck

  PUBLIC Field_TypeGet

  PUBLIC Field_UserNumberFind

  PUBLIC Field_VariableCheck
  
  PUBLIC Field_VariableGet

  PUBLIC Field_VariableIndexGet

  PUBLIC Field_VariableLabelGet

  PUBLIC Field_VariableTypeCheck

  PUBLIC Field_VariableTypesCheck

  PUBLIC Field_VariableTypesGet

  PUBLIC FieldInterpolatedPoint_InterpolationParametersGet

  PUBLIC FieldInterpolatedPointMetrics_InterpolatedPointGet

  PUBLIC FieldInterpolationParameters_FieldVariableGet

  PUBLIC FieldParameterSet_ParametersGet

  PUBLIC FieldPhysicalPoint_FieldInterpolatedPointGet

  PUBLIC FieldPhysicalPoint_GeometricInterpolatedPointGet

  PUBLIC FieldVariable_AssertIsINTGData,FieldVariable_AssertIsSPData,FieldVariable_AssertIsDPData,FieldVariable_AssertIsLData
  
  PUBLIC FieldVariable_AssertComponentNumberOK
  
  PUBLIC FieldVariable_ComponentDomainGet

  PUBLIC FieldVariable_ComponentInterpolationCheck

  PUBLIC FieldVariable_ComponentInterpolationGet

  PUBLIC FieldVariable_ComponentLabelGet

  PUBLIC FieldVariable_ComponentMeshComponentCheck

  PUBLIC FieldVariable_ComponentMeshComponentGet
  
  PUBLIC FieldVariable_DataTypeCheck

  PUBLIC FieldVariable_DataTypeGet

  PUBLIC FieldVariable_DimensionCheck

  PUBLIC FieldVariable_DimensionGet
  
  PUBLIC FieldVariable_DOFOrderTypeCheck

  PUBLIC FieldVariable_DOFOrderTypeGet

  PUBLIC FieldVariable_DomainGet

  PUBLIC FieldVariable_DomainMappingGet

  PUBLIC FieldVariable_FieldGet

  PUBLIC FieldVariable_ConstantDOFGet
  
  PUBLIC FieldVariable_LocalElementDOFGet

  PUBLIC FieldVariable_UserElementDOFGet

  PUBLIC FieldVariable_LocalNodeDOFGet

  PUBLIC FieldVariable_UserNodeDOFGet

  PUBLIC FieldVariable_LocalGaussDOFGet

  PUBLIC FieldVariable_UserGaussDOFGet

  PUBLIC FieldVariable_LocalDataPointDOFGet

  PUBLIC FieldVariable_UserDataPointDOFGet

  PUBLIC FieldVariable_LocalElementDataDOFGet

  PUBLIC FieldVariable_UserElementDataDOFGet

  PUBLIC FieldVariable_LabelGet

  PUBLIC FieldVariable_NumberOfComponentsCheck

  PUBLIC FieldVariable_NumberOfComponentsGet

  PUBLIC FieldVariable_ParameterSetCheck

  PUBLIC FieldVariable_ParameterSetGet

  PUBLIC FieldVariable_VariableTypeGet
  
  PUBLIC Fields_RegionGet

CONTAINS

  !
  !=================================================================================================================================
  !

  !>Assert that a field has been finished
  SUBROUTINE Field_AssertIsFinished(field,err,error,*)

    !Argument Variables
    TYPE(FieldType), POINTER, INTENT(IN) :: field !<The field to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Field_AssertIsFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)

    IF(.NOT.field%fieldFinished) THEN
      localError="Field number "//TRIM(NumberToVString(field%userNumber,"*",err,error))// &
        & " has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Field_AssertIsFinished")
    RETURN
999 ERRORSEXITS("Field_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Field_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a field has not been finished
  SUBROUTINE Field_AssertNotFinished(field,err,error,*)

    !Argument Variables
    TYPE(FieldType), POINTER, INTENT(IN) :: field !<The field to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Field_AssertNotFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)

    IF(field%fieldFinished) THEN
      localError="Field number "//TRIM(NumberToVString(field%userNumber,"*",err,error))// &
        & " has already been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Field_AssertNotFinished")
    RETURN
999 ERRORSEXITS("Field_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Field_AssertNotFinished
  
  !
  !================================================================================================================================
  !

  !>Returns the dof numbers for a field variable component that corresponds to the specified constant
  SUBROUTINE Field_ComponentDOFGetConstant(field,variableType,componentNumber,localDOF,globalDOF,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the dof for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the dof for \see FieldRoutines_VariableTypes,FieldRoutines
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The field component number to get the dof for
    INTEGER(INTG), INTENT(OUT) :: localDOF !<On exit, the local dof corresponding to the constant
    INTEGER(INTG), INTENT(OUT) :: globalDOF !<On exit, the global dof corresponding to the constant
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DomainMappingType), POINTER :: domainMapping
    TYPE(FieldVariableType), POINTER :: fieldVariable

    ENTERS("Field_ComponentDOFGetConstant",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL FieldVariable_ConstantDOFGet(fieldVariable,componentNumber,localDOF,err,error,*999)
    NULLIFY(domainMapping)
    CALL FieldVariable_DomainMappingGet(fieldVariable,domainMapping,err,error,*999)
    globalDOF=domainMapping%localToGlobalMap(localDOF)

    EXITS("Field_ComponentDOFGetConstant")
    RETURN
999 ERRORSEXITS("Field_ComponentDOFGetConstant",err,error)
    RETURN 1

  END SUBROUTINE Field_ComponentDOFGetConstant

  !
  !================================================================================================================================
  !

  !>Returns the dof numbers for a field component that corresponds to the specified user data point.
  SUBROUTINE Field_ComponentDOFGetUserDataPoint(field,variableType,userDataPointNumber,componentNumber,localDof, &
    & globalDof,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the dof for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the dof for \see FieldRoutines_VariableTypes,FieldRoutines
    INTEGER(INTG), INTENT(IN) :: userDataPointNumber !<The user data point number to get the dof for
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The field component number to get the dof for
    INTEGER(INTG), INTENT(OUT) :: localDof !<On exit, the local dof corresponding to the user data point
    INTEGER(INTG), INTENT(OUT) :: globalDof !<On exit, the global dof corresponding to the user data point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    LOGICAL :: ghostDof
    TYPE(DomainMappingType), POINTER :: domainMapping
    TYPE(FieldVariableType), POINTER :: fieldVariable

    ENTERS("Field_ComponentDOFGetUserDataPoint",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL FieldVariable_UserDataPointDofGet(fieldVariable,userDataPointNumber,componentNumber,localDof,ghostDOf,err,error,*999)
    NULLIFY(domainMapping)
    CALL FieldVariable_DomainMappingGet(fieldVariable,domainMapping,err,error,*999)
    globalDof=domainMapping%localToGlobalMap(localDof)

    EXITS("Field_ComponentDOFGetUserDataPoint")
    RETURN
999 ERRORSEXITS("Field_ComponentDOFGetUserDataPoint",err,error)
    RETURN 1

  END SUBROUTINE Field_ComponentDOFGetUserDataPoint

  !
  !================================================================================================================================
  !

  !>Returns the dof numbers for a field component that corresponds to the specified user element.
  SUBROUTINE Field_ComponentDOFGetUserElement(field,variableType,userElementNumber,componentNumber,localDOF, &
    & globalDOF,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the dof for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the dof for \see FieldRoutines_VariableTypes,FieldRoutines
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to get the dof for
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The field component number to get the dof for
    INTEGER(INTG), INTENT(OUT) :: localDOF !<On exit, the local dof corresponding to the user element
    INTEGER(INTG), INTENT(OUT) :: globalDOF !<On exit, the global dof corresponding to the user element
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    LOGICAL :: ghostDOF
    TYPE(DomainMappingType), POINTER :: domainMapping
    TYPE(FieldVariableType), POINTER :: fieldVariable

    ENTERS("Field_ComponentDOFGetUserElement",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL FieldVariable_UserElementDofGet(fieldVariable,userElementNumber,componentNumber,localDOF,ghostDOF,err,error,*999)
    NULLIFY(domainMapping)
    CALL FieldVariable_DomainMappingGet(fieldVariable,domainMapping,err,error,*999)
    globalDOF=domainMapping%localToGlobalMap(localDOF)

    EXITS("Field_ComponentDOFGetUserElement")
    RETURN
999 ERRORSEXITS("Field_ComponentDOFGetUserElement",err,error)
    RETURN 1

  END SUBROUTINE Field_ComponentDOFGetUserElement

  !
  !================================================================================================================================
  !
  
  !>Returns the dof numbers for a field component that corresponds to the specified user node and derivative.
  SUBROUTINE Field_ComponentDOFGetUserNode(field,variableType,versionNumber,derivativeNumber,userNodeNumber, & 
    & componentNumber,localDOF,globalDOF,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the dof for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the dof for \see FieldRoutines_VariableTypes,FieldRoutines
    INTEGER(INTG), INTENT(IN) :: versionNumber !<The version number to get the dof for
    INTEGER(INTG), INTENT(IN) :: derivativeNumber !<The derivative number to get the dof for
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number to get the dof for
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The field component number to get the dof for
    INTEGER(INTG), INTENT(OUT) :: localDOF !<On exit, the local dof corresponding to the user node
    INTEGER(INTG), INTENT(OUT) :: globalDOF !<On exit, the global dof corresponding to the user node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    LOGICAL :: ghostDOF
    TYPE(DomainMappingType), POINTER :: domainMapping
    TYPE(FieldVariableType), POINTER :: fieldVariable

    ENTERS("Field_ComponentDOFGetUserNode",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL FieldVariable_UserNodeDofGet(fieldVariable,versionNumber,derivativeNumber,userNodeNumber,componentNumber, &
      & localDOF,ghostDOF,err,error,*999)
    NULLIFY(domainMapping)
    CALL FieldVariable_DomainMappingGet(fieldVariable,domainMapping,err,error,*999)
    globalDOF=domainMapping%localToGlobalMap(localDOF)

    EXITS("Field_ComponentDOFGetUserNode")
    RETURN
999 ERRORSEXITS("Field_ComponentDOFGetUserNode",err,error)
    RETURN 1

  END SUBROUTINE Field_ComponentDOFGetUserNode

  !
  !================================================================================================================================
  !

  !>Checks the interpolation type for a field variable component.
  SUBROUTINE Field_ComponentInterpolationCheck(field,variableType,componentNumber,interpolationType,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to check the interpolation for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type of the field variable component to check \see FieldRoutines_VariableTypes,FieldRoutines
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The field component number of the field variable component to check
    INTEGER(INTG), INTENT(IN) :: interpolationType !<The interpolation type of the field variable component to check \see FieldRoutines_VariableTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldVariableType), POINTER :: fieldVariable

    ENTERS("Field_ComponentInterpolationCheck",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL FieldVariable_ComponentInterpolationCheck(fieldVariable,componentNumber,interpolationType,err,error,*999)    
 
    EXITS("Field_ComponentInterpolationCheck")
    RETURN
999 ERRORSEXITS("Field_ComponentInterpolationCheck",err,error)
    RETURN 1

  END SUBROUTINE Field_ComponentInterpolationCheck

  !
  !================================================================================================================================
  !

  !>Gets the interpolation type for a field variable component identified by a pointer. \see OpenCMISS::Iron::cmfe_Field_ComponentInterpolationGet
  SUBROUTINE Field_ComponentInterpolationGet(field,variableType,componentNumber,interpolationType,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the interpolation type for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type of the field variable component to get \see FieldRoutines_VariableTypes,FieldRoutines
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The field component number of the field variable component to set
    INTEGER(INTG), INTENT(OUT) :: interpolationType !<On return, the interpolation type of the field variable component \see FieldRoutines_InterpolationTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldVariableType), POINTER :: fieldVariable

    ENTERS("Field_ComponentInterpolationGet",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL FieldVariable_ComponentInterpolationGet(fieldVariable,componentNumber,interpolationType,err,error,*999)   

    EXITS("Field_ComponentInterpolationGet")
    RETURN
999 ERRORSEXITS("Field_ComponentInterpolationGet",err,error)
    RETURN 1
    
  END SUBROUTINE Field_ComponentInterpolationGet

  !
  !================================================================================================================================
  !

  !>Gets the label for a field variable component for character labels. \see OpenCMISS::Iron::cmfe_Field_ComponentLabelGet
  SUBROUTINE Field_ComponentLabelGetC(field,variableType,componentNumber,label,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the label for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type \see FieldRoutines_VariableTypes,FieldRoutines
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The field component number to get
    CHARACTER(LEN=*), INTENT(OUT) :: label !<On return, the field variable label
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldVariableType), POINTER :: fieldVariable

    ENTERS("Field_ComponentLabelGetC",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL FieldVariable_ComponentLabelGet(fieldVariable,componentNumber,label,err,error,*999)

    EXITS("Field_ComponentLabelGetC")
    RETURN
999 ERRORSEXITS("Field_ComponentLabelGetC",err,error)
    RETURN 1

  END SUBROUTINE Field_ComponentLabelGetC

  !
  !================================================================================================================================
  !

  !>Gets the label for a field variable component for varying string labels. \see OpenCMISS::Iron::cmfe_Field_ComponentLabelGet
  SUBROUTINE Field_ComponentLabelGetVS(field,variableType,componentNumber,label,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the label for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type \see FieldRoutines_VariableTypes,FieldRoutines
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The field component number to get
    TYPE(VARYING_STRING), INTENT(OUT) :: label !<On return, the field variable label
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldVariableType), POINTER :: fieldVariable

    ENTERS("Field_ComponentLabelGetVS",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL FieldVariable_ComponentLabelGet(fieldVariable,componentNumber,label,err,error,*999)

    EXITS("Field_ComponentLabelGetVS")
    RETURN
999 ERRORSEXITS("Field_ComponentLabelGetVS",err,error)
    RETURN 1

  END SUBROUTINE Field_ComponentLabelGetVS

  !
  !================================================================================================================================
  !

  !>Checks the mesh component number type for a field variable component.
  SUBROUTINE Field_ComponentMeshComponentCheck(field,variableType,componentNumber,meshComponent,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to check the mesh component for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type of the field variable component to check \see FieldRoutines_VariableTypes,FieldRoutines
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The field component number of the field variable component to check
    INTEGER(INTG), INTENT(IN) :: meshComponent !<The mesh component to check for the specified field variable component
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldVariableType), POINTER :: fieldVariable

    ENTERS("Field_ComponentMeshComponentCheck",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL FieldVariable_ComponentMeshComponentCheck(fieldVariable,componentNumber,meshComponent,err,error,*999)    
 
    EXITS("Field_ComponentMeshComponentCheck")
    RETURN
999 ERRORSEXITS("Field_ComponentMeshComponentCheck",err,error)
    RETURN 1

  END SUBROUTINE Field_ComponentMeshComponentCheck

  !
  !================================================================================================================================
  !

  !>Gets the mesh component number for a field variable component identified by a pointer. \see OpenCMISS::Iron::cmfe_Field_ComponentMeshComponentnGet
  SUBROUTINE Field_ComponentMeshComponentGet(field,variableType,componentNumber,meshComponent,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the mesh component type for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type of the field variable component to get \see FieldRoutines_VariableTypes,FieldRoutines
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The field component number of the field variable component to set
    INTEGER(INTG), INTENT(OUT) :: meshComponent !<On return, the mesh component for the specified field variable component
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldVariableType), POINTER :: fieldVariable

    ENTERS("Field_ComponentMeshComponentGet",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL FieldVariable_ComponentMeshComponentGet(fieldVariable,componentNumber,meshComponent,err,error,*999)   

    EXITS("Field_ComponentMeshComponentGet")
    RETURN
999 ERRORSEXITS("Field_ComponentMeshComponentGet",err,error)
    RETURN 1
    
  END SUBROUTINE Field_ComponentMeshComponentGet

  !
  !================================================================================================================================
  !

  !>Returns the coordinate system for a field accounting for regions and interfaces
  SUBROUTINE Field_CoordinateSystemGet(field,coordinateSystem,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER, INTENT(IN) :: field !<A pointer to the field to get the coordinate system for
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem!<On return, the field coordinate system. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(InterfaceType), POINTER :: interface
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_CoordinateSystemGet",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is already associated.",err,error,*999)

    NULLIFY(coordinateSystem)
    NULLIFY(interface)
    region=>field%region
    IF(ASSOCIATED(region)) THEN
      coordinateSystem=>region%coordinateSystem
      IF(.NOT.ASSOCIATED(coordinateSystem)) THEN
        localError="The coordinate system is not associated for field number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))//" of region number "// &
          & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      interface=>field%interface
      IF(ASSOCIATED(interface)) THEN
        coordinateSystem=>interface%coordinateSystem
        IF(.NOT.ASSOCIATED(coordinateSystem)) THEN
          localError="The coordinate system is not associated for field number "// &
            & TRIM(NumberToVString(field%userNumber,"*",err,error))//" of interface number "// &
            & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        localError="A region or interface is not associated for field number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF

    EXITS("Field_CoordinateSystemGet")
    RETURN
999 ERRORSEXITS("Field_CoordinateSystemGet",err,error)
    RETURN 1
    
  END SUBROUTINE Field_CoordinateSystemGet

  !
  !================================================================================================================================
  !

  !>Gets the create values cache from a field.
  SUBROUTINE Field_CreateValuesCacheGet(field,createValuesCache,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<The field to get the create values cache for.
    TYPE(FieldCreateValuesCacheType), POINTER :: createValuesCache !<On exit, a pointer to the create values cache for the field. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Field_CreateValuesCacheGet",err,error,*998)

    !Check input arguments
    IF(ASSOCIATED(createValuesCache)) CALL FlagError("Create values cache is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)

    createValuesCache=>field%createValuesCache
    IF(.NOT.ASSOCIATED(createValuesCache)) THEN
      localError="Create values cache is not associated for field "//TRIM(NumberToVString(field%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Field_CreateValuesCacheGet")
    RETURN
999 NULLIFY(createValuesCache)
998 ERRORSEXITS("Field_CreateValuesCacheGet",err,error)
    RETURN 1
    
    
  END SUBROUTINE Field_CreateValuesCacheGet

  !
  !================================================================================================================================
  !

  !>Gets a data project from a field.
  SUBROUTINE Field_DataProjectionGet(field,dataProjection,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<The field to get the data projection for.
    TYPE(DataProjectionType), POINTER :: dataProjection !<On exit, a pointer to the data projection for the field. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Field_DataProjectionGet",err,error,*998)

    !Check input arguments
    IF(ASSOCIATED(dataProjection)) CALL FlagError("Data projection is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)

    dataProjection=>field%dataProjection
    IF(.NOT.ASSOCIATED(dataProjection)) THEN
      localError="Data projection  is not associated for field "//TRIM(NumberToVString(field%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Field_DataProjectionGet")
    RETURN
999 NULLIFY(dataProjection)
998 ERRORSEXITS("Field_DataProjectionGet",err,error)
    RETURN 1
    
  END SUBROUTINE Field_DataProjectionGet

  !
  !================================================================================================================================
  !

  !>Checks the data type for a field variable.
  SUBROUTINE Field_DataTypeCheck(field,variableType,dataType,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to check the data type for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to check \see FieldRoutines_VariableTypes,FieldRoutines
    INTEGER(INTG), INTENT(IN) :: dataType !<The data type of the field variable to check \see FieldRoutines_DataTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldVariableType), POINTER :: fieldVariable

    ENTERS("Field_DataTypeCheck",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL FieldVariable_DataTypeCheck(fieldVariable,dataType,err,error,*999)

    EXITS("Field_DataTypeCheck")
    RETURN
999 ERRORSEXITS("Field_DataTypeCheck",err,error)
    RETURN 1

  END SUBROUTINE Field_DataTypeCheck

  !
  !================================================================================================================================
  !

  !>Gets the data type for a field variable. \see OpenCMISS::Iron::cmfe_Field_DataTypeGet
  SUBROUTINE Field_DataTypeGet(field,variableType,dataType,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the data type for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get \see FieldRoutines_VariableTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: dataType !<On return, the data type of the field variable \see FieldRoutines_DataTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldVariableType), POINTER :: fieldVariable

    ENTERS("Field_DataTypeGet",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)

    dataType=fieldVariable%dataType

    EXITS("Field_DataTypeGet")
    RETURN
999 ERRORSEXITS("Field_DataTypeGet",err,error)
    RETURN 1

  END SUBROUTINE Field_DataTypeGet

  !
  !================================================================================================================================
  !

  !>Gets a decomposition from a field. \see OpenCMISS::Iron::cmfe_Field_DecompositionGet
  SUBROUTINE Field_DecompositionGet(field,decomposition,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<The field to get the decomposition for.
    TYPE(DecompositionType), POINTER :: decomposition !<On exit, a pointer to the decomposition for the field. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Field_DecompositionGet",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(ASSOCIATED(decomposition)) CALL FlagError("Decomposition is already associated.",err,error,*999)

    !Get the field decomposition
    decomposition=>field%decomposition
    !Check field decomposition is associated.
    IF(.NOT.ASSOCIATED(decomposition)) THEN
      localError="Field decomposition is not associated for field number "// &
        & TRIM(NumberToVString(field%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Field_DecompositionGet")
    RETURN
999 ERRORSEXITS("Field_DecompositionGet",err,error)
    RETURN 1
    
  END SUBROUTINE Field_DecompositionGet

  !
  !================================================================================================================================
  !

  !>Checks the dependent type for a field.
  SUBROUTINE Field_DependentTypeCheck(field,dependentType,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to check the dependent type for
    INTEGER(INTG), INTENT(IN) :: dependentType !<The dependent type to check \see FieldRoutines_DependentTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_DependentTypeCheck",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    
    SELECT CASE(dependentType)
    CASE(FIELD_INDEPENDENT_TYPE)
      IF(field%dependentType/=FIELD_INDEPENDENT_TYPE) THEN
        localError="Invalid dependent type. The dependent type of field number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))//" is "// &
          & TRIM(NumberToVString(field%dependentType,"*",err,error))// &
          & " which is not an independent field."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(FIELD_DEPENDENT_TYPE)
      IF(field%dependentType/=FIELD_DEPENDENT_TYPE) THEN
        localError="Invalid dependent type. The dependent type of field number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))//" is "// &
          & TRIM(NumberToVString(field%dependentType,"*",err,error))// &
          & " which is not a dependent field."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The specified dependent type of "//TRIM(NumberToVString(dependentType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Field_DependentTypeCheck")
    RETURN
999 ERRORSEXITS("Field_DependentTypeCheck",err,error)
    RETURN 1

  END SUBROUTINE Field_DependentTypeCheck

  !
  !================================================================================================================================
  !

  !>Gets the dependent type for a field. \see OpenCMISS::Iron::cmfe_Field_DependentTypeGet
  SUBROUTINE Field_DependentTypeGet(field,dependentType,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the dependent type for
    INTEGER(INTG), INTENT(OUT) :: dependentType !<On return, the dependent type to get \see FieldRoutines_DependentTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Field_DependentTypeGet",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)

    dependentType=field%dependentType

    EXITS("Field_DependentTypeGet")
    RETURN
999 ERRORSEXITS("Field_DependentTypeGet",err,error)
    RETURN 1

  END SUBROUTINE Field_DependentTypeGet

  !
  !================================================================================================================================
  !

  !>Checks the field dimension for a field variable.
  SUBROUTINE Field_DimensionCheck(field,variableType,dimensionType,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to check the dimension for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to check \see FieldRoutines_VariableTypes,FieldRoutines 
    INTEGER(INTG), INTENT(IN) :: dimensionType !<The field dimension to check \see FieldRoutines_DimensionTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldVariableType), POINTER :: fieldVariable

    ENTERS("Field_DimensionCheck",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL FieldVariable_DimensionCheck(fieldVariable,dimensionType,err,error,*999)

    EXITS("Field_DimensionCheck")
    RETURN
999 ERRORSEXITS("Field_DimensionCheck",err,error)
    RETURN 1

  END SUBROUTINE Field_DimensionCheck

  !
  !================================================================================================================================
  !

  !>Gets the field dimension for a field variable. \see OpenCMISS::Iron::cmfe_Field_DimensionGet
  SUBROUTINE Field_DimensionGet(field,variableType,dimension,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the dimension for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type \see FieldRoutines_VariableTypes,FieldRoutines 
    INTEGER(INTG), INTENT(OUT) :: dimension !<On return, the field dimension to get \see FieldRoutines_DimensionTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldVariableType), POINTER :: fieldVariable

    ENTERS("Field_DimensionGet",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)

    dimension=fieldVariable%dimension

    EXITS("Field_DimensionGet")
    RETURN
999 ERRORSEXITS("Field_DimensionGet",err,error)
    RETURN 1

  END SUBROUTINE Field_DimensionGet

  !
  !================================================================================================================================
  !

  !>Checks the DOF order type for a field variable.
  SUBROUTINE Field_DOFOrderTypeCheck(field,variableType,dofOrderType,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to check the DOF order type for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to check \see FieldRoutines_VariableTypes,FieldRoutines
    INTEGER(INTG), INTENT(IN) :: dofOrderType !<The DOF order type of the field variable to check \see FieldRoutines_DOFOrderTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldVariableType), POINTER :: fieldVariable

    ENTERS("Field_DOFOrderTypeCheck",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL FieldVariable_DOFOrderTypeCheck(fieldVariable,dofOrderType,err,error,*999)

    EXITS("Field_DOFOrderTypeCheck")
    RETURN
999 ERRORSEXITS("Field_DOFOrderTypeCheck",err,error)
    RETURN 1

  END SUBROUTINE Field_DOFOrderTypeCheck

  !
  !================================================================================================================================
  !

  !>Gets the DOF order type for a field variable. \see OpenCMISS::Iron::cmfe_Field_DOFOrderTypeGet
  SUBROUTINE Field_DOFOrderTypeGet(field,variableType,dofOrderType,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the DOF order type for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get \see FieldRoutines_VariableTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: dofOrderType !<On return, the DOF order type of the field variable \see FieldRoutines_DOFOrderTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldVariableType), POINTER :: fieldVariable

    ENTERS("Field_DOFOrderTypeGet",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)

    dofOrderType=fieldVariable%dofOrderType

    EXITS("Field_DOFOrderTypeGet")
    RETURN
999 ERRORSEXITS("Field_DOFOrderTypeGet",err,error)
    RETURN 1

  END SUBROUTINE Field_DOFOrderTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the fields for a field.
  SUBROUTINE Field_FieldsGet(field,fields,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the fields for
    TYPE(FieldsType), POINTER :: fields !<On return, a pointer to the field fields. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_FieldsGet",err,error,*998)

    IF(ASSOCIATED(fields)) CALL FlagError("Fields is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
     
    fields=>field%fields
    IF(.NOT.ASSOCIATED(fields)) THEN
      localError="The fields for field number "//TRIM(NumberToVString(field%userNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("Field_FieldsGet")
    RETURN
999 NULLIFY(fields)
998 ERRORSEXITS("Field_FieldsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Field_FieldsGet

  !
  !================================================================================================================================
  !

  !>Gets the geometric field for a field identified by a pointer. \see OpenCMISS::Iron::cmfe_Field_GeometricFieldGet
  SUBROUTINE Field_GeometricFieldGet(field,geometricField,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the geometric field for
    TYPE(FieldType), POINTER :: geometricField !<On return, a pointer to the geometric field. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_GeometricFieldGet",err,error,*998)

    IF(ASSOCIATED(geometricField)) CALL FlagError("Geometric field is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
     
    geometricField=>field%geometricField
    IF(.NOT.ASSOCIATED(geometricField)) THEN
      localError="The geometric field for field number "//TRIM(NumberToVString(field%userNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("Field_GeometricFieldGet")
    RETURN
999 NULLIFY(geometricField)
998 ERRORSEXITS("Field_GeometricFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE Field_GeometricFieldGet

  !
  !================================================================================================================================
  !

  !>Gets a geometric general field for a field if there is any (eg. the dependent field for a finite elasticity equation),
  !>otherwise the normal geometric field is returned if present. If no geometric field is found then an error is raised.
  SUBROUTINE Field_GeometricGeneralFieldGet(field,geometricField,generalFound,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER, INTENT(IN) :: field !<A pointer to the field to get the geometric field for
    TYPE(FieldType), POINTER, INTENT(OUT) :: geometricField !<On return, a pointer to the geometric field. Must not be associated on entry.
    LOGICAL, INTENT(OUT) :: generalFound !<On return, true if we found a geometric general field, otherwise false.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: fieldIdx
    TYPE(FieldType), POINTER :: otherField
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_GeometricGeneralFieldGet",err,error,*999)

    IF(ASSOCIATED(geometricField)) CALL FlagError("Geometric field is already associated.",err,error,*999)
    CALL Field_AssertIsFinished(field,err,error,*999)
    IF(.NOT.ASSOCIATED(field%fields)) CALL FlagError("Field fields are not associated.",err,error,*999)

    generalFound=.FALSE.
    ! Find the geometric general field associated with this field
    DO fieldIdx=1,field%fields%numberOfFields
      otherField=>field%fields%fields(fieldIdx)%ptr
      IF(.NOT.ASSOCIATED(otherField)) THEN
        localError="Field index "//TRIM(NumberToVString(fieldIdx,"*",err,error))//" is not associated."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(otherField%TYPE==FIELD_GEOMETRIC_GENERAL_TYPE) THEN
        geometricField=>otherField
        generalFound=.TRUE.
      ENDIF
    ENDDO !fieldIdx

    IF(.NOT.generalFound) THEN
      ! We couldn't find a geometric general field. Just return the undeformed geometric field.
      IF(.NOT.ASSOCIATED(field%geometricField)) &
        & CALL FlagError("Geometric general field not found and geometric field is not associated.",err,error,*999)
      geometricField=>field%geometricField
    END IF

    EXITS("Field_GeometricGeneralFieldGet")
    RETURN
999 ERRORSEXITS("Field_GeometricGeneralFieldGet",err,error)
    RETURN 1

  END SUBROUTINE Field_GeometricGeneralFieldGet

  !
  !================================================================================================================================
  !

  !>Gets the geometric parameters for a field. 
  SUBROUTINE Field_GeometricParametersGet(field,geometricParameters,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the geometric field for
    TYPE(FieldGeometricParametersType), POINTER :: geometricParameters !<On return, a pointer to the geometric parameters. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_GeometricParametersGet",err,error,*998)

    IF(ASSOCIATED(geometricParameters)) CALL FlagError("Geometric parameters is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
     
    geometricParameters=>field%geometricFieldParameters
    IF(.NOT.ASSOCIATED(geometricParameters)) THEN
      localError="The geometric field parameters for field number "//TRIM(NumberToVString(field%userNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("Field_GeometricParametersGet")
    RETURN
999 NULLIFY(geometricParameters)
998 ERRORSEXITS("Field_GeometricParametersGet",err,error)
    RETURN 1
    
  END SUBROUTINE Field_GeometricParametersGet

  !
  !================================================================================================================================
  !

  !>Returns the interface for a field
  SUBROUTINE Field_InterfaceGet(field,interface,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the interface for
    TYPE(InterfaceType), POINTER :: interface !<On return, the field interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_InterfaceGet",err,error,*998)

    !Check input arguments
    IF(ASSOCIATED(interface)) CALL FlagError("Interface is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)

    INTERFACE=>field%INTERFACE
    IF(.NOT.ASSOCIATED(interface)) THEN
      localError="The interface for field number "// &
        & TRIM(NumberToVString(field%userNumber,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
 
    EXITS("Field_InterfaceGet")
    RETURN
999 NULLIFY(interface)
998 ERRORSEXITS("Field_InterfaceGet",err,error)
    RETURN 1
    
  END SUBROUTINE Field_InterfaceGet

  !
  !================================================================================================================================
  !

  !>Determines if the given field is an interface field or not. 
  SUBROUTINE Field_IsInterfaceField(field,interfaceField,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to determine if it is an interface field or not.
    LOGICAL :: interfaceField !<On exit, .TRUE. if the given field is in an interface region, .FALSE. if not. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Field_IsInterfaceField",err,error,*999)

    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)

    interfaceField = ASSOCIATED(field%interface)
    
    EXITS("Field_IsInterfaceField")
    RETURN
999 ERRORSEXITS("Field_IsInterfaceField",err,error)
    RETURN 1
    
  END SUBROUTINE Field_IsInterfaceField

  !
  !================================================================================================================================
  !

  !>Determines if the given field is a region field or not. 
  SUBROUTINE Field_IsRegionField(field,regionField,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to determine if it is a region field or not.
    LOGICAL :: regionField !<On exit, .TRUE. if the given field is in a region, .FALSE. if not. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Field_IsRegionField",err,error,*999)

    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)

    regionField = ASSOCIATED(field%region)
    
    EXITS("Field_IsRegionField")
    RETURN
999 ERRORSEXITS("Field_IsRegionField",err,error)
    RETURN 1
    
  END SUBROUTINE Field_IsRegionField

  !
  !================================================================================================================================
  !

  !>Gets the field label for a field for character labels. \see OpenCMISS::Iron::cmfe_Field_LabelGet
  SUBROUTINE Field_LabelGetC(field,label,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the label for
    CHARACTER(LEN=*), INTENT(OUT) :: label !<On return, the field label for the specified field
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: cLength,vsLength

    ENTERS("Field_LabelGetC",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    cLength=LEN(label)
    vsLength=LEN_TRIM(field%label)
    IF(cLength>vsLength) THEN
      label=CHAR(LEN_TRIM(field%label))
    ELSE
      label=CHAR(field%label,cLength)
    ENDIF

    EXITS("Field_LabelGetC")
    RETURN
999 ERRORSEXITS("Field_LabelGetC",err,error)
    RETURN 1

  END SUBROUTINE Field_LabelGetC

  !
  !================================================================================================================================
  !

  !>Gets the field label for a field for varying string labels. \see OpenCMISS::Iron::cmfe_Field_LabelGet
  SUBROUTINE Field_LabelGetVS(field,label,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the label for
    TYPE(VARYING_STRING), INTENT(OUT) :: label !<On return, the field label for the specified field
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Field_LabelGetVS",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    
    label=field%label

    EXITS("Field_LabelGetVS")
    RETURN
999 ERRORSEXITS("Field_LabelGetVS",err,error)
    RETURN 1
    
  END SUBROUTINE Field_LabelGetVS

  !
  !================================================================================================================================
  !

  !>Checks the number of field components for a field variable.
  SUBROUTINE Field_NumberOfComponentsCheck(field,variableType,numberOfComponents,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to check the number of components
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to check \see FieldRoutines_VariableTypes,FieldRoutines 
    INTEGER(INTG), INTENT(IN) :: numberOfComponents !The number of components in the field variable to check
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_NumberOfComponentsCheck",err,error,*999)

    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    
    IF(field%fieldFinished) THEN
      NULLIFY(fieldVariable)
      CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
      IF(fieldVariable%numberOfComponents/=numberOfComponents) THEN
        localError="Invalid number of components. The number components for variable type "// &
          & TRIM(NumberToVString(variableType,"*",err,error))//" of field number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))//" is "// &
          & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))// &
          & " which does not correspond to the specified number of components of "// &
          & TRIM(NumberToVString(numberOFComponents,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      !Field has not been finished so check the create values cache.
      NULLIFY(createValuesCache)
      CALL Field_CreateValuesCacheGet(field,createValuesCache,err,error,*999)
      IF(.NOT.ALLOCATED(createValuesCache%numberOfComponents)) THEN
        localError="The create values cache number of components is not allocated for field number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(createValuesCache%numberOfComponents(variableType)/=numberOfComponents) THEN
        fieldVariable=>field%variableTypeMap(variableType)%ptr
        IF(ASSOCIATED(fieldVariable)) THEN
          localError="Invalid number of components. The number components for variable type "// &
            & TRIM(NumberToVString(variableType,"*",err,error))//" of field number "// &
            & TRIM(NumberToVString(field%userNumber,"*",err,error))//" is "// &
            & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))// &
            & " which does not correspond to the specified number of components of "// &
            & TRIM(NumberToVString(numberOFComponents,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ELSE
          localError="Invalid number of components. The number components for variable type "// &
            & TRIM(NumberToVString(variableType,"*",err,error))//" of field number "// &
            & TRIM(NumberToVString(field%userNumber,"*",err,error))// &
            & " does not correspond to the specified number of components of "// &
            & TRIM(NumberToVString(numberOFComponents,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    ENDIF

    EXITS("Field_NumberOfComponentsCheck")
    RETURN
999 ERRORSEXITS("Field_NumberOfComponentsCheck",err,error)
    RETURN 1

  END SUBROUTINE Field_NumberOfComponentsCheck

  !
  !================================================================================================================================
  !

  !>Gets the number of field components for a field variable. \see OpenCMISS::Iron::cmfe_Field_NumberOfComponentsGet
  SUBROUTINE Field_NumberOfComponentsGet(field,variableType,numberOfComponents,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the number of components
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get \see FieldRoutines_VariableTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: numberOfComponents !<On return, the number of components in the field variable
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_NumberOfComponentsGet",err,error,*999)

    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    
    IF(field%fieldFinished) THEN
      NULLIFY(fieldVariable)
      CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
      numberOfComponents=fieldVariable%numberOfComponents
    ELSE
      !Field has not been finished so check the create values cache.
      NULLIFY(createValuesCache)
      CALL Field_CreateValuesCacheGet(field,createValuesCache,err,error,*999)
      IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
        localError="The supplied variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
          & " is invalid. The field variable type must be > 1 and <= "// &
          & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(.NOT.ALLOCATED(createValuesCache%numberOfComponents)) THEN
        localError="The create values cache number of components is not allocated for field number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      numberOfComponents=createValuesCache%numberOfComponents(variableType)
    ENDIF

    EXITS("Field_NumberOfComponentsGet")
    RETURN
999 ERRORSEXITS("Field_NumberOfComponentsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Field_NumberOfComponentsGet

  !
  !================================================================================================================================
  !

  !>Checks the number of variables for a field.
  SUBROUTINE Field_NumberOfVariablesCheck(field,numberOfVariables,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to check the number of variables for
    INTEGER(INTG), INTENT(IN) :: numberOfVariables !<The number of variables in the specified field to check
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_NumberOfVariablesCheck",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    IF(field%numberOfVariables/=numberOfVariables) THEN
      localError="Invalid number of variables. The number of variables for field number "// &
        & TRIM(NumberToVString(field%userNumber,"*",err,error))//" is "// &
        & TRIM(NumberToVString(field%numberOfVariables,"*",err,error))// &
        & " which does correspond to the specified number of variables of "// &
        & TRIM(NumberToVString(numberOfVariables,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("Field_NumberOfVariablesCheck")
    RETURN
999 ERRORSEXITS("Field_NumberOfVariablesCheck",err,error)
    RETURN 1

  END SUBROUTINE Field_NumberOfVariablesCheck

  !
  !================================================================================================================================
  !

  !>Gets the number of variables for a field. \see OpenCMISS::Iron::cmfe_Field_NumberOfVariablesGet
  SUBROUTINE Field_NumberOfVariablesGet(field,numberOfVariables,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the number of variables for
    INTEGER(INTG), INTENT(OUT) :: numberOfVariables !<On return, the number of variables in the specified field
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Field_NumberOfVariablesGet",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    
    numberOfVariables=field%numberOfVariables
 
    EXITS("Field_NumberOfVariablesGet")
    RETURN
999 ERRORSEXITS("Field_NumberOfVariablesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Field_NumberOfVariablesGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the specified parameter set for the field variable.
  SUBROUTINE Field_ParameterSetGet(field,variableType,fieldSetType,parameterSet,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the parameter set for
    INTEGER(INTG), INTENT(IN) :: variableType!<The field variable type to get the parameter set for \see FieldRoutines_VariableTypes,FieldRoutines
    INTEGER(INTG), INTENT(IN) :: fieldSetType !<The field parameter set identifier to get \see FieldRoutines_ParameterSetTypes,FieldRoutines
    TYPE(FieldParameterSetType), POINTER :: parameterSet !<On return, a pointer to the specified parameter set. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldVariableType), POINTER :: fieldVariable

    ENTERS("Field_ParameterSetGet",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL FieldVariable_ParameterSetGet(fieldVariable,fieldSetType,parameterSet,err,error,*999)

    EXITS("Field_ParameterSetGet")
    RETURN
999 ERRORSEXITS("Field_ParameterSetGet",err,error)
    RETURN 1

  END SUBROUTINE Field_ParameterSetGet

  !
  !================================================================================================================================
  !

  !>Returns the region for a field accounting for regions and interfaces
  SUBROUTINE Field_RegionGet(field,region,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the region for
    TYPE(RegionType), POINTER :: region !<On return, the fields region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(InterfaceType), POINTER :: interface
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_RegionGet",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(ASSOCIATED(region)) CALL FlagError("Region is already associated.",err,error,*999)
        
    NULLIFY(region)
    NULLIFY(interface)
    region=>field%region
    IF(.NOT.ASSOCIATED(region)) THEN          
      INTERFACE=>field%INTERFACE
      IF(ASSOCIATED(INTERFACE)) THEN
        IF(ASSOCIATED(interface%parentRegion)) THEN
          region=>interface%parentRegion     
        ELSE
          localError="The parent region is not associated for field number "// &
            & TRIM(NumberToVString(field%userNumber,"*",err,error))//" of interface number "// &
            & TRIM(NumberToVString(interface%userNumber,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        localError="A region or interface is not associated for field number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF
    
    EXITS("Field_RegionGet")
    RETURN
999 ERRORSEXITS("Field_RegionGet",err,error)
    RETURN 1
    
  END SUBROUTINE Field_RegionGet

  !
  !================================================================================================================================
  !

  !>Returns the scale factors vector for a field scaling index
  SUBROUTINE Field_ScaleFactorsVectorGet(field,scalingIndex,scaleFactorsVector,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the scale factors vector for
    INTEGER(INTG), INTENT(IN) :: scalingIndex !<The scaling index of the scale factors vector to get
    TYPE(DistributedVectorType), POINTER :: scaleFactorsVector !<On return, the field scale factor vector for the scaling index. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_ScaleFactorsVectorGet",err,error,*998)

    IF(ASSOCIATED(scaleFactorsVector)) CALL FlagError("Scale factors vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(scalingIndex<1.OR.scalingIndex>field%scalings%numberOfScalingIndices) THEN
      localError="The specified scaling index of "//TRIM(NumberToVString(scalingIndex,"*",err,error))// &
        & " is invalid for field number "//TRIM(NumberToVString(field%userNumber,"*",err,error))// &
        & ". The scaling index should be >= 1 and <= "// &
        & TRIM(NumberToVString(field%scalings%numberOfScalingIndices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    IF(.NOT.ALLOCATED(field%scalings%scalings)) THEN
      localError="The field scalings is not allocated for field number "// &
        & TRIM(NumberToVString(field%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    scaleFactorsVector=>field%scalings%scalings(scalingIndex)%scaleFactors
    IF(.NOT.ASSOCIATED(scaleFactorsVector)) THEN
      localError="The scale factors vector is not associated for scaling index number "// &
        & TRIM(NumberToVString(scalingIndex,"*",err,error))//" of field number "// &
        & TRIM(NumberToVString(field%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    
    EXITS("Field_ScaleFactorsVectorGet")
    RETURN
999 NULLIFY(scaleFactorsVector)
998 ERRORSEXITS("Field_ScaleFactorsVectorGet",err,error)
    RETURN 1
    
  END SUBROUTINE Field_ScaleFactorsVectorGet

  !
  !================================================================================================================================
  !

  !>Checks the scaling type for a field.
  SUBROUTINE Field_ScalingTypeCheck(field,scalingType,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to check the scaling type for
    INTEGER(INTG), INTENT(IN) :: scalingType !<The scaling type for the specified field to check \see FieldRoutines_ScalingTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_ScalingTypeCheck",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    
    SELECT CASE(scalingType)
    CASE(FIELD_NO_SCALING)
      IF(field%scalings%scalingType/=FIELD_NO_SCALING) THEN
        localError="Invalid scaling type. The scaling type for field number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))//" is "// &
          & TRIM(NumberToVString(field%scalings%scalingType,"*",err,error))// &
          & " which is not no scaling."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(FIELD_UNIT_SCALING)
      IF(field%scalings%scalingType/=FIELD_UNIT_SCALING) THEN
        localError="Invalid scaling type. The scaling type for field number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))//" is "// &
          & TRIM(NumberToVString(field%scalings%scalingType,"*",err,error))// &
          & " which is not unit scaling."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(FIELD_ARC_LENGTH_SCALING)
      IF(field%scalings%scalingType/=FIELD_ARC_LENGTH_SCALING) THEN
        localError="Invalid scaling type. The scaling type for field number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))//" is "// &
          & TRIM(NumberToVString(field%scalings%scalingType,"*",err,error))// &
          & " which is not arc length scaling."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(FIELD_ARITHMETIC_MEAN_SCALING)
      IF(field%scalings%scalingType/=FIELD_ARITHMETIC_MEAN_SCALING) THEN
        localError="Invalid scaling type. The scaling type for field number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))//" is "// &
          & TRIM(NumberToVString(field%scalings%scalingType,"*",err,error))// &          
          & " which is not arithmetic mean scaling."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(FIELD_GEOMETRIC_MEAN_SCALING)
      IF(field%scalings%scalingType/=FIELD_GEOMETRIC_MEAN_SCALING) THEN
        localError="Invalid scaling type. The scaling type for field number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))//" is "// &
          & TRIM(NumberToVString(field%scalings%scalingType,"*",err,error))// &          
          & " which is not geometric mean scaling."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(FIELD_HARMONIC_MEAN_SCALING)
      IF(field%scalings%scalingType/=FIELD_HARMONIC_MEAN_SCALING) THEN
        localError="Invalid scaling type. The scaling type for field number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))//" is "// &
          & TRIM(NumberToVString(field%scalings%scalingType,"*",err,error))// &
          & " which is not harmonic mean scaling."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The specified scaling type of "//TRIM(NumberToVString(scalingType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Field_ScalingTypeCheck")
    RETURN
999 ERRORSEXITS("Field_ScalingTypeCheck",err,error)
    RETURN 1
    
  END SUBROUTINE Field_ScalingTypeCheck

  !
  !================================================================================================================================
  !

  !>Gets the scaling type for a field. \see OpenCMISS::Iron::cmfe_Field_ScalingTypeGet
  SUBROUTINE Field_ScalingTypeGet(field,scalingType,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the scaling type for
    INTEGER(INTG), INTENT(OUT) :: scalingType !<On return, the scaling type for the specified field to get \see FieldRoutines_ScalingTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Field_ScalingTypeGet",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    
    scalingType=field%scalings%scalingType
 
    EXITS("Field_ScalingTypeGet")
    RETURN
999 ERRORSEXITS("Field_ScalingTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Field_ScalingTypeGet

  !
  !================================================================================================================================
  !

  !>Checks the field type for a field.
  SUBROUTINE Field_TypeCheck(field,type,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to check the type for
    INTEGER(INTG), INTENT(IN) :: type !<The field type to check \see FieldRoutines_FieldTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_TypeCheck",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    
    SELECT CASE(type)
    CASE(FIELD_GEOMETRIC_TYPE)
      IF(field%type/=FIELD_GEOMETRIC_TYPE) THEN
        localError="Invalid field type. The field type for field number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))//" is "// &
          & TRIM(NumberToVString(field%type,"*",err,error))// &
          & " which is not a geometric field."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(FIELD_FIBRE_TYPE)
      IF(field%type/=FIELD_FIBRE_TYPE) THEN
        localError="Invalid field type. The field type for field number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))//" is "// &
          & TRIM(NumberToVString(field%type,"*",err,error))// &
          & " which is not a fibre field."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(FIELD_GENERAL_TYPE)
      IF(field%type/=FIELD_GENERAL_TYPE) THEN
        localError="Invalid field type. The field type for field number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))//" is "// &
          & TRIM(NumberToVString(field%type,"*",err,error))// &
          & " which is not a general field."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(FIELD_MATERIAL_TYPE)
      IF(field%type/=FIELD_MATERIAL_TYPE) THEN
        localError="Invalid field type. The field type for field number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))//" is "// &
          & TRIM(NumberToVString(field%type,"*",err,error))// &
          & " which is not a material field."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(FIELD_GEOMETRIC_GENERAL_TYPE)
      IF(field%type/=FIELD_GEOMETRIC_GENERAL_TYPE) THEN
        localError="Invalid field type. The field type for field number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))//" is "// &
          & TRIM(NumberToVString(field%type,"*",err,error))// &
          & " which is not a geometric general field."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The specified field type of "//TRIM(NumberToVString(type,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)          
    END SELECT

    EXITS("Field_TypeCheck")
    RETURN
999 ERRORSEXITS("Field_TypeCheck",err,error)
    RETURN 1
    
  END SUBROUTINE Field_TypeCheck

  !
  !================================================================================================================================
  !

  !>Gets the field type for a field. \see OpenCMISS::Iron::cmfe_FieldTypeGet
  SUBROUTINE Field_TypeGet(field,type,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the type for
    INTEGER(INTG), INTENT(OUT) :: type !<On return, the field type for the specified field \see FieldRoutines_FieldTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Field_TypeGet",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    
    type=field%type

    EXITS("Field_TypeGet")
    RETURN
999 ERRORSEXITS("Field_TypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Field_TypeGet

  !
  !================================================================================================================================
  !

  !>Finds and returns in field a pointer to the field identified by a user number in the given list of fields. If no field with that user number exists field is left nullified.
  SUBROUTINE Field_UserNumberFindGeneric(userNumber,fields,field,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The field user number to find
    TYPE(FieldsType), POINTER :: fields !<The list of fields containing the field
    TYPE(FieldType), POINTER :: field !<On return, a pointer to the field with the given user number. If no field with that user number exists in the list of fields the field is null. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: fieldIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_UserNumberFindGeneric",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(fields)) CALL FlagError("Fields is not associated.",err,error,*999)
    IF(ASSOCIATED(field)) CALL FlagError("Field is already associated.",err,error,*999)
   
    !Get the field from the user number
    NULLIFY(field)
    IF(ASSOCIATED(fields%fields)) THEN
      DO fieldIdx=1,fields%numberOfFields
        IF(ASSOCIATED(fields%fields(fieldIdx)%ptr)) THEN
          IF(fields%fields(fieldIdx)%ptr%userNumber==userNumber) THEN
            field=>fields%fields(fieldIdx)%ptr
            EXIT
          ENDIF
        ELSE
          localError="The field pointer in fields is not associated for field index "// &
            & TRIM(NumberToVString(fieldIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !fieldIdx
    ENDIF
      
    EXITS("Field_UserNumberFindGeneric")
    RETURN
999 ERRORSEXITS("Field_UserNumberFindGeneric",err,error)
    RETURN 1
    
  END SUBROUTINE Field_UserNumberFindGeneric

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the field identified by user number in a given interface. If no field with that user number exists a null pointer is returned.
  SUBROUTINE Field_UserNumberFindInterface(userNumber,interface,field,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The field user number to find
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface containing the field
    TYPE(FieldType), POINTER :: field !<On exit, a pointer to the field with the given user number. If no field with that user number exists in the interface the field is null. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Field_UserNumberFindInterface",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*999)

    CALL Field_UserNumberFindGeneric(userNumber,interface%fields,field,err,error,*999)
  
    EXITS("Field_UserNumberFindInterface")
    RETURN
999 ERRORSEXITS("Field_UserNumberFindInterface",err,error)
    RETURN 1
    
  END SUBROUTINE Field_UserNumberFindInterface

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the field identified by user number in a given region. If no field with that user number exists a null pointer is returned.
  SUBROUTINE Field_UserNumberFindRegion(userNumber,region,field,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The field user number to find
    TYPE(RegionType), POINTER :: region !<A pointer to the region containing the field
    TYPE(FieldType), POINTER :: field !<On exit, a pointer to the field with the given user number. If no field with that user number exists in the region the field is null. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Field_UserNumberFindRegion",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)

    CALL Field_UserNumberFindGeneric(userNumber,region%fields,field,err,error,*999)
  
    EXITS("Field_UserNumberFindRegion")
    RETURN
999 ERRORSEXITS("Field_UserNumberFindRegion",err,error)
    RETURN 1
    
  END SUBROUTINE Field_UserNumberFindRegion

  !
  !================================================================================================================================
  !

  !>Gets the label for a field variable for character labels. \see OpenCMISS::Iron::cmfe_Field_VariableLabelGet
  SUBROUTINE Field_VariableLabelGetC(field,variableType,label,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the label for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type \see FieldRoutines_VariableTypes,FieldRoutines
    CHARACTER(LEN=*), INTENT(OUT) :: label !<On return, the field variable label
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldVariableType), POINTER :: fieldVariable

    ENTERS("Field_VariableLabelGetC",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL FieldVariable_LabelGet(fieldVariable,label,err,error,*999)

    EXITS("Field_VariableLabelGetC")
    RETURN
999 ERRORSEXITS("Field_VariableLabelGetC",err,error)
    RETURN 1

  END SUBROUTINE Field_VariableLabelGetC

  !
  !================================================================================================================================
  !

  !>Gets the label for a field variable for varying string labels. \see OpenCMISS::Iron::cmfe_Field_VariableLabelGet
  SUBROUTINE Field_VariableLabelGetVS(field,variableType,label,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the label for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type \see FieldRoutines_VariableTypes,FieldRoutines
    TYPE(VARYING_STRING), INTENT(OUT) :: label !<On return, the field variable label
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldVariableType), POINTER :: fieldVariable

    ENTERS("Field_VariableLabelGetVS",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL FieldVariable_LabelGet(fieldVariable,label,err,error,*999)

    EXITS("Field_VariableLabelGetVS")
    RETURN
999 ERRORSEXITS("Field_VariableLabelGetVS",err,error)
    RETURN 1

  END SUBROUTINE Field_VariableLabelGetVS

  !
  !================================================================================================================================
  !

  !>Returns a pointer to interpolation parameters for a field interpolated point
  SUBROUTINE FieldInterpolatedPoint_InterpolationParametersGet(interpolatedPoint,interpolationParameters,err,error,*)

    !Argument variables
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint !<A pointer to the interpolated point to get the interpolation parameters for.
    TYPE(FieldInterpolationParametersType), POINTER :: interpolationParameters !<On exit, a pointer to the interpolation parameters for the interpolated point. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("FieldInterpolatedPoint_InterpolationParametersGet",err,error,*998)

    IF(ASSOCIATED(interpolationParameters)) CALL FlagError("Interpolation parameters is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interpolatedPoint)) CALL FlagError("Field interpolated point is not associated.",err,error,*999)

    interpolationParameters=>interpolatedPoint%interpolationParameters
    IF(.NOT.ASSOCIATED(interpolationParameters)) &
      & CALL FlagError("Interpolation parameters is not associated for the interpolated point.",err,error,*999)

    EXITS("FieldInterpolatedPoint_InterpolationParametersGet")
    RETURN
999 NULLIFY(interpolationParameters)
998 ERRORS("FieldInterpolatedPoint_InterpolationParametersGet",err,error)
    EXITS("FieldInterpolatedPoint_InterpolationParametersGet")
    RETURN 1
    
  END SUBROUTINE FieldInterpolatedPoint_InterpolationParametersGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to interpolated point for a field interpolated point metrics
  SUBROUTINE FieldInterpolatedPointMetrics_InterpolatedPointGet(interpolatedPointMetrics,interpolatedPoint,err,error,*)

    !Argument variables
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: interpolatedPointMetrics !<A pointer to the interpolated point metrics to get the interpolated point for.
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint !<On exit, a pointer to the interpolated point for the interpolated point metrics. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("FieldInterpolatedPointMetrics_InterpolatedPointGet",err,error,*998)

    IF(ASSOCIATED(interpolatedPoint)) CALL FlagError("Interpolated point is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interpolatedPointMetrics)) &
      & CALL FlagError("Field interpolated point metrics is not associated.",err,error,*999)

    interpolatedPoint=>interpolatedPointMetrics%interpolatedPoint
    IF(.NOT.ASSOCIATED(interpolatedPoint)) &
      & CALL FlagError("Interpolated point is not associated for the interpolated point metrics.",err,error,*999)

    EXITS("FieldInterpolatedPointMetrics_InterpolatedPointGet")
    RETURN
999 NULLIFY(interpolatedPoint)
998 ERRORS("FieldInterpolatedPointMetrics_InterpolatedPointGet",err,error)
    EXITS("FieldInterpolatedPointMetrics_InterpolatedPointGet")
    RETURN 1
    
  END SUBROUTINE FieldInterpolatedPointMetrics_InterpolatedPointGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to a field variable for a field interpolation parameters
  SUBROUTINE FieldInterpolationParameters_FieldVariableGet(interpolationParameters,fieldVariable,err,error,*)

    !Argument variables
    TYPE(FieldInterpolationParametersType), POINTER :: interpolationParameters !<A pointer to the interpolation parameters to get the field variable for.
    TYPE(FieldVariableType), POINTER :: fieldVariable !<On exit, a pointer to the field variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("FieldInterpolationParameters_FieldVariableGet",err,error,*998)

    IF(ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interpolationParameters)) CALL FlagError("Field interpolation parameters is not associated.",err,error,*999)

    fieldVariable=>interpolationParameters%fieldVariable
    IF(.NOT.ASSOCIATED(fieldVariable)) &
      & CALL FlagError("The field variable is not associated for the interpolation parameters.",err,error,*999)

    EXITS("FieldInterpolationParameters_FieldVariableGet")
    RETURN
999 NULLIFY(fieldVariable)
998 ERRORS("FieldInterpolationParameters_FieldVariableGet",err,error)
    EXITS("FieldInterpolationParameters_FieldVariableGet")
    RETURN 1
    
  END SUBROUTINE FieldInterpolationParameters_FieldVariableGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the parameters distributed vector for a field parameter set.
  SUBROUTINE FieldParameterSet_ParametersGet(parameterSet,parameters,err,error,*)

    !Argument variables
    TYPE(FieldParameterSetType), POINTER :: parameterSet !<A pointer to the parameter set to get the parameters vector for.
    TYPE(DistributedVectorType), POINTER :: parameters !<On exit, a pointer to the parameters vector for the parameter set. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("FieldParameterSet_ParametersGet",err,error,*998)

    IF(ASSOCIATED(parameters)) CALL FlagError("Parameters is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(parameterSet)) CALL FlagError("Field parameter set is not associated.",err,error,*999)

    parameters=>parameterSet%parameters
    IF(.NOT.ASSOCIATED(parameters)) &
      & CALL FlagError("Parameters is not associated for the parameter set.",err,error,*999)

    EXITS("FieldParameterSet_ParametersGet")
    RETURN
999 NULLIFY(parameters)
998 ERRORSEXITS("FieldParameterSet_ParametersGet",err,error)
    RETURN 1
    
  END SUBROUTINE FieldParameterSet_ParametersGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the field interpolated point for a field physical point
  SUBROUTINE FieldPhysicalPoint_FieldInterpolatedPointGet(physicalPoint,fieldInterpolatedPoint,err,error,*)

    !Argument variables
    TYPE(FieldPhysicalPointType), POINTER :: physicalPoint !<A pointer to the physical point to get the field interpolated point for.
    TYPE(FieldInterpolatedPointType), POINTER :: fieldInterpolatedPoint !<On exit, a pointer to the field interpolated point for the physical point. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("FieldPhysicalPoint_FieldInterpolatedPointGet",err,error,*998)

    IF(ASSOCIATED(fieldInterpolatedPoint)) CALL FlagError("Field interpolated point is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(physicalPoint)) CALL FlagError("Field physical point is not associated.",err,error,*999)

    fieldInterpolatedPoint=>physicalPoint%fieldInterpolatedPoint
    IF(.NOT.ASSOCIATED(fieldInterpolatedPoint)) &
      & CALL FlagError("Field interpolated point is not associated for the physical point.",err,error,*999)

    EXITS("FieldPhysicalPoint_FieldInterpolatedPointGet")
    RETURN
999 NULLIFY(fieldInterpolatedPoint)
998 ERRORS("FieldPhysicalPoint_FieldInterpolatedPointGet",err,error)
    EXITS("FieldPhysicalPoint_FieldInterpolatedPointGet")
    RETURN 1
    
  END SUBROUTINE FieldPhysicalPoint_FieldInterpolatedPointGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the geometric interpolated point for a field physical point
  SUBROUTINE FieldPhysicalPoint_GeometricInterpolatedPointGet(physicalPoint,geometricInterpolatedPoint,err,error,*)

    !Argument variables
    TYPE(FieldPhysicalPointType), POINTER :: physicalPoint !<A pointer to the physical point to get the geometric interpolated point for.
    TYPE(FieldInterpolatedPointType), POINTER :: geometricInterpolatedPoint !<On exit, a pointer to the geometric interpolated point for the physical point. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("FieldPhysicalPoint_GeometricInterpolatedPointGet",err,error,*998)

    IF(ASSOCIATED(geometricInterpolatedPoint)) CALL FlagError("Geometric interpolated point is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(physicalPoint)) CALL FlagError("Field physical point is not associated.",err,error,*999)

    geometricInterpolatedPoint=>physicalPoint%geometricInterpolatedPoint
    IF(.NOT.ASSOCIATED(geometricInterpolatedPoint)) &
      & CALL FlagError("Geometric interpolated point is not associated for the physical point.",err,error,*999)

    EXITS("FieldPhysicalPoint_GeometricInterpolatedPointGet")
    RETURN
999 NULLIFY(geometricInterpolatedPoint)
998 ERRORS("FieldPhysicalPoint_GeometricInterpolatedPointGet",err,error)
    EXITS("FieldPhysicalPoint_GeometricInterpolatedPointGet")
    RETURN 1
    
  END SUBROUTINE FieldPhysicalPoint_GeometricInterpolatedPointGet

  !
  !=================================================================================================================================
  !

  !>Assert that a field variable has integer data type
  SUBROUTINE FieldVariable_AssertIsINTGData(fieldVariable,err,error,*)

    !Argument Variables
    TYPE(FieldVariableType), POINTER, INTENT(INOUT) :: fieldVariable !<The field variable to assert the integer data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("FieldVariable_AssertIsINTGData",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)

    IF(fieldVariable%dataType/=FIELD_INTG_TYPE) THEN
      localError="The field variable data type of "//TRIM(NumberToVString(fieldVariable%dataType,"*",err,error))// &
        & " does not correspond to the required integer data type for field variable type "// &
        & TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("FieldVariable_AssertIsINTGData")
    RETURN
999 ERRORSEXITS("FieldVariable_AssertIsINTGData",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_AssertIsINTGData

  !
  !=================================================================================================================================
  !

  !>Assert that a field variable has single precision real data type
  SUBROUTINE FieldVariable_AssertIsSPData(fieldVariable,err,error,*)

    !Argument Variables
    TYPE(FieldVariableType), POINTER, INTENT(INOUT) :: fieldVariable !<The field variable to assert the single precision data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("FieldVariable_AssertIsSPData",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)

    IF(fieldVariable%dataType/=FIELD_SP_TYPE) THEN
      localError="The field variable data type of "//TRIM(NumberToVString(fieldVariable%dataType,"*",err,error))// &
        & " does not correspond to the required single precision real data type for field variable type "// &
        & TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("FieldVariable_AssertIsSPData")
    RETURN
999 ERRORSEXITS("FieldVariable_AssertIsSPData",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_AssertIsSPData

  !
  !=================================================================================================================================
  !

  !>Assert that a field variable has double precision real data type
  SUBROUTINE FieldVariable_AssertIsDPData(fieldVariable,err,error,*)

    !Argument Variables
    TYPE(FieldVariableType), POINTER, INTENT(INOUT) :: fieldVariable !<The field variable to assert the double precision data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("FieldVariable_AssertIsDPData",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)

    IF(fieldVariable%dataType/=FIELD_DP_TYPE) THEN
      localError="The field variable data type of "//TRIM(NumberToVString(fieldVariable%dataType,"*",err,error))// &
        & " does not correspond to the required double precision real data type for field variable type "// &
        & TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("FieldVariable_AssertIsDPData")
    RETURN
999 ERRORSEXITS("FieldVariable_AssertIsDPData",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_AssertIsDPData

  !
  !=================================================================================================================================
  !

  !>Assert that a field variable has logical data type
  SUBROUTINE FieldVariable_AssertIsLData(fieldVariable,err,error,*)

    !Argument Variables
    TYPE(FieldVariableType), POINTER, INTENT(INOUT) :: fieldVariable !<The field variable to assert the logical data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("FieldVariable_AssertIsLData",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)

    IF(fieldVariable%dataType/=FIELD_L_TYPE) THEN
      localError="The field variable data type of "//TRIM(NumberToVString(fieldVariable%dataType,"*",err,error))// &
        & " does not correspond to the required logical data type for field variable type "// &
        & TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("FieldVariable_AssertIsLData")
    RETURN
999 ERRORSEXITS("FieldVariable_AssertIsLData",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_AssertIsLData

  !
  !================================================================================================================================
  !

  !>Asserts that the component number of a field variable is valid.
  SUBROUTINE FieldVariable_AssertComponentNumberOK(fieldVariable,componentNumber,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to assert the component number for
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number to assert
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_AssertComponentNumberOK",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    
    IF(componentNumber<1.OR.componentNumber>fieldVariable%numberOfComponents) THEN
      localError="The specified component number of "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " is invalid for field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &        
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//". The field variable component must be >= 1 and <= "// &
        & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("FieldVariable_AssertComponentNumberOK")
    RETURN
999 ERRORSEXITS("FieldVariable_AssertComponentNumberOK",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_AssertComponentNumberOK

  !
  !================================================================================================================================
  !

  !>Checks if a field has a field variable of a specified type
  SUBROUTINE Field_VariableCheck(field,variableType,fieldVariable,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to check the variable for.
    INTEGER(INTG), INTENT(IN) :: variableType !<The type of field variable to check. \see FieldRoutines_VariableTypes,FieldRoutines
    TYPE(FieldVariableType), POINTER :: fieldVariable !<On exit, a pointer to the field variable if it exists. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_VariableCheck",err,error,*998)

    IF(ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(variableType<0.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The field variable type must be between 1 and "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    fieldVariable=>field%variableTypeMap(variableType)%ptr

    EXITS("Field_VariableCheck")
    RETURN
999 NULLIFY(fieldVariable)
998 ERRORSEXITS("Field_VariableCheck",err,error)
    RETURN 1
    
  END SUBROUTINE Field_VariableCheck

  !
  !================================================================================================================================
  !

  !>Returns a pointer to a field variable of a specified type
  SUBROUTINE Field_VariableGet(field,variableType,fieldVariable,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the variable for.
    INTEGER(INTG), INTENT(IN) :: variableType !<The type of field variable to set. \see FieldRoutines_VariableTypes,FieldRoutines
    TYPE(FieldVariableType), POINTER :: fieldVariable !<On exit, a pointer to the field variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_VariableGet",err,error,*998)

    IF(ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(variableType<0.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The field variable type must be between 1 and "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    fieldVariable=>field%variableTypeMap(variableType)%ptr
    IF(.NOT.ASSOCIATED(fieldVariable)) THEN
      localError="The field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " has not been defined on field number "//TRIM(NumberToVString(field%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("Field_VariableGet")
    RETURN
999 NULLIFY(fieldVariable)
998 ERRORSEXITS("Field_VariableGet",err,error)
    RETURN 1
    
  END SUBROUTINE Field_VariableGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to a field variable of a specified index
  SUBROUTINE Field_VariableIndexGet(field,variableIdx,fieldVariable,variableType,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the variable for.
    INTEGER(INTG), INTENT(IN) :: variableIdx !<The index of the field variable to set.
    TYPE(FieldVariableType), POINTER :: fieldVariable !<On exit, a pointer to the field variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: variableType !<On exit, the variable type of the specified variable index.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_VariableIndexGet",err,error,*998)

    IF(ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(field%variables)) CALL FlagError("Field variables is not allocated.",err,error,*999)
    IF(variableIdx<0.OR.variableIdx>SIZE(field%variables,1)) THEN
      localError="The specified field variable index of "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is invalid. The field variable index must be between 1 and "// &
        & TRIM(NumberToVString(SIZE(field%variables,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    variableType=field%variables(variableIdx)%variableType
    fieldVariable=>field%variableTypeMap(variableType)%ptr
    
    IF(.NOT.ASSOCIATED(fieldVariable)) THEN
      localError="The field variable index of "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is not associated for field number "//TRIM(NumberToVString(field%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("Field_VariableIndexGet")
    RETURN
999 NULLIFY(fieldVariable)
998 ERRORSEXITS("Field_VariableIndexGet",err,error)
    RETURN 1
    
  END SUBROUTINE Field_VariableIndexGet

  !
  !================================================================================================================================
  !

  !>Checks the field contains the given field variable type.
  SUBROUTINE Field_VariableTypeCheck(field,variableType,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to check the variable type for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to check for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: variableIdx
    LOGICAL :: variableFound
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_VariableTypeCheck",err,error,*999)

    variableFound=.FALSE.
    CALL Field_AssertIsFinished(field,err,error,*999)
    DO variableIdx=1,field%numberOfVariables
      IF(field%variables(variableIdx)%variableType==variableType) THEN
        variableFound=.TRUE.
        EXIT
      ENDIF
    ENDDO !variableIdx
    IF(.NOT.variableFound) THEN
      localError="Field number "//TRIM(NumberToVString(field%userNumber,"*",err,error))// &
        & " does not contain a variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("Field_VariableTypeCheck")
    RETURN
999 ERRORSEXITS("Field_VariableTypeCheck",err,error)
    RETURN 1
    
  END SUBROUTINE Field_VariableTypeCheck

  !
  !================================================================================================================================
  !

  !>Checks the field variable types for a field.
  SUBROUTINE Field_VariableTypesCheck(field,variableTypes,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to check the variable types for
    INTEGER(INTG), INTENT(IN) :: variableTypes(:) !<variableTypes(variableIdx). The field variable type for the variableIdx'th field variable to check
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: variableIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_VariableTypesCheck",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    IF(SIZE(variableTypes,1)<field%numberOfVariables) THEN
      localError="Invalid variable types. The size of the specified variable types array is "// &
        & TRIM(NumberToVString(SIZE(variableTypes,1),"*",err,error))//" and it must be >= "// &
        & TRIM(NumberToVString(field%numberOfVariables,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    DO variableIdx=1,field%numberOfVariables
      IF(variableTypes(variableIdx)<1.OR.variableTypes(variableIdx)>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
        localError="The specified variable type of "//TRIM(NumberToVString(variableTypes(variableIdx),"*",err,error))// &
          & " at position number "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
          & " is invalid. The variable type must be >= 1 and <= "// &
          & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(field%variables(variableIdx)%variableType/=variableTypes(variableIdx)) THEN
        localError="Invalid variable type. The variable type for variable index number "// &
          & TRIM(NumberToVString(variableIdx,"*",err,error))//" of field number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))//" is "// &
          & TRIM(NumberToVString(field%variables(variableIdx)%variableType,"*",err,error))// &
          & " which is does correspond to the specified variable_type of "// &
          & TRIM(NumberToVString(variableTypes(variableIdx),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDDO !variableIdx

    EXITS("Field_VariableTypesCheck")
    RETURN
999 ERRORSEXITS("Field_VariableTypesCheck",err,error)
    RETURN 1
    
  END SUBROUTINE Field_VariableTypesCheck

  !
  !================================================================================================================================
  !

  !>Gets the field variable types for a field. \see OpenCMISS::Iron::cmfe_Field_VariableTypesGet
  SUBROUTINE Field_VariableTypesGet(field,variableTypes,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the variable types for
    INTEGER(INTG), INTENT(OUT) :: variableTypes(:) !<variableTypes(variableIdx). On return, the field variable type variableIdx'th field variable
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: variableIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_VariableTypesGet",err,error,*999)

    CALL Field_AssertIsFinished(field,err,error,*999)
    IF(SIZE(variableTypes,1)<field%numberOfVariables) THEN
      localError="Invalid variable types. The size of the specified variable types array is "// &
        & TRIM(NumberToVString(SIZE(variableTypes,1),"*",err,error))//" and it must be >= "// &
        & TRIM(NumberToVString(field%numberOfVariables,"*",err,error))//" for field number "// &
        & TRIM(NumberToVString(field%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    variableTypes=0
    DO variableIdx=1,field%numberOfVariables
      variableTypes(variableIdx)=field%variables(variableIdx)%variableType
    ENDDO !variableIdx

    EXITS("Field_VariableTypesGet")
    RETURN
999 ERRORSEXITS("Field_VariableTypesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Field_VariableTypesGet

  !
  !================================================================================================================================
  !

  !>Returns the domain of the specified field variable component
  SUBROUTINE FieldVariable_ComponentDomainGet(fieldVariable,componentIdx,domain,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the domain for
    INTEGER(INTG), INTENT(IN) :: componentIdx !<The component index of the field variable to get the domain for. 
    TYPE(DomainType), POINTER :: domain  !<On exit, the domain for field variable component. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_ComponentDomainGet",err,error,*998)

    IF(ASSOCIATED(domain)) CALL FlagError("Domain is already assocaiated.",err,error,*998)
    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    IF(componentIdx<1.OR.componentIdx>fieldVariable%numberOfComponents) THEN
      localError="The specified component number of "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
        & " is invalid for field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &        
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//". The field variable component must be >= 1 and <= "// &
        & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components)) THEN
      localError="Field variable components have not been allocated"// &
        & " for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &        
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    domain=>fieldVariable%components(componentIdx)%domain
    IF(.NOT.ASSOCIATED(domain)) THEN
      localError="Domain is not associated"// &
        & " for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &        
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
 
    EXITS("FieldVariable_ComponentDomainGet")
    RETURN
999 NULLIFY(domain)
998 ERRORSEXITS("FieldVariable_ComponentDomainGet",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_ComponentDomainGet

  !
  !================================================================================================================================
  !

  !>Checks the interpolation type for a field variable component.
  SUBROUTINE FieldVariable_ComponentInterpolationCheck(fieldVariable,componentNumber,interpolationType,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to check the interpolation for
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The field component number of the field variable component to check
    INTEGER(INTG), INTENT(IN) :: interpolationType !<The interpolation type of the field variable component to check \see FieldRoutines_VariableTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_ComponentInterpolationCheck",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(fieldVariable%components)) THEN
      localError="The components are not allocated for variable type "// &
        & TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) &
        localError=localError//" of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(componentNumber<1.OR.componentNumber>fieldVariable%numberOfComponents) THEN
      localError="The specified omponent number of "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " is invalid for variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) &
        & localError=localError//" of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//". The component number should be >= 1 and <= "// &
        & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(interpolationType)
    CASE(FIELD_CONSTANT_INTERPOLATION)
      IF(fieldVariable%components(componentNumber)%interpolationType/=FIELD_CONSTANT_INTERPOLATION) THEN
        localError="Invalid interpolation type. The interpolation type for component number "// &
          & TRIM(NumberToVString(componentNumber,"*",err,error))// &
          & " of variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
        IF(ASSOCIATED(fieldVariable%field)) &
          & localError=localError// " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
        localError=localError//" is "// &
          & TRIM(NumberToVString(fieldVariable%components(componentNumber)%interpolationType,"*",err,error))// &
          & " which is not constant interpolation."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
      IF(fieldVariable%components(componentNumber)%interpolationType/=FIELD_ELEMENT_BASED_INTERPOLATION) THEN
        localError="Invalid interpolation type. The interpolation type for component number "// &
          & TRIM(NumberToVString(componentNumber,"*",err,error))// &
          & " of variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
        IF(ASSOCIATED(fieldVariable%field)) &
          & localError=localError// " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
        localError=localError//" is "// &
          & TRIM(NumberToVString(fieldVariable%components(componentNumber)%interpolationType,"*",err,error))// &
          & " which is not element based interpolation."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(FIELD_NODE_BASED_INTERPOLATION)
      IF(fieldVariable%components(componentNumber)%interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
        localError="Invalid interpolation type. The interpolation type for component number "// &
          & TRIM(NumberToVString(componentNumber,"*",err,error))// &
          & " of variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
        IF(ASSOCIATED(fieldVariable%field)) &
          & localError=localError// " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
        localError=localError//" is "// &
          & TRIM(NumberToVString(fieldVariable%components(componentNumber)%interpolationType,"*",err,error))// &
          & " which is not node based interpolation."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
      IF(fieldVariable%components(componentNumber)%interpolationType/=FIELD_GRID_POINT_BASED_INTERPOLATION) THEN
        localError="Invalid interpolation type. The interpolation type for component number "// &
          & TRIM(NumberToVString(componentNumber,"*",err,error))// &
          & " of variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
        IF(ASSOCIATED(fieldVariable%field)) &
          & localError=localError// " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
        localError=localError//" is "// &
          & TRIM(NumberToVString(fieldVariable%components(componentNumber)%interpolationType,"*",err,error))// &
          & " which is not grid point based interpolation."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
      IF(fieldVariable%components(componentNumber)%interpolationType/=FIELD_GAUSS_POINT_BASED_INTERPOLATION) THEN
        localError="Invalid interpolation type. The interpolation type for component number "// &
          & TRIM(NumberToVString(componentNumber,"*",err,error))// &
          & " of variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
        IF(ASSOCIATED(fieldVariable%field)) &
          & localError=localError// " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
        localError=localError//" is "// &
          & TRIM(NumberToVString(fieldVariable%components(componentNumber)%interpolationType,"*",err,error))// &
          & " which is not Gauss point based interpolation."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
      IF(fieldVariable%components(componentNumber)%interpolationType/=FIELD_DATA_POINT_BASED_INTERPOLATION) THEN
        localError="Invalid interpolation type. The interpolation type for component number "// &
          & TRIM(NumberToVString(componentNumber,"*",err,error))// &
          & " of variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
        IF(ASSOCIATED(fieldVariable%field)) &
          & localError=localError// " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
        localError=localError//" is "// &
          & TRIM(NumberToVString(fieldVariable%components(componentNumber)%interpolationType,"*",err,error))// &
          & " which is not data point based interpolation."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The specified interpolation type of "//TRIM(NumberToVString(interpolationType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("FieldVariable_ComponentInterpolationCheck")
    RETURN
999 ERRORSEXITS("FieldVariable_ComponentInterpolationCheck",err,error)
    RETURN 1

  END SUBROUTINE FieldVariable_ComponentInterpolationCheck

   !
  !================================================================================================================================
  !

  !>Returns the interpolation type of the specified field variable component
  SUBROUTINE FieldVariable_ComponentInterpolationGet(fieldVariable,componentNumber,interpolationType,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the interpolation type for
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number of the field variable to get the interpolation type for. 
    INTEGER(INTG), INTENT(OUT) :: interpolationType  !<On exit, the interpolation type for field variable component. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_ComponentInterpolationGet",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(fieldVariable%components)) THEN
      localError="The components are not allocated for variable type "// &
        & TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) &
        localError=localError//" of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(componentNumber<1.OR.componentNumber>fieldVariable%numberOfComponents) THEN
      localError="The specified omponent number of "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " is invalid for variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) &
        & localError=localError//" of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//". The component number should be >= 1 and <= "// &
        & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    interpolationType=fieldVariable%components(componentNumber)%interpolationType

    EXITS("FieldVariable_ComponentInterpolationGet")
    RETURN
999 ERRORSEXITS("FieldVariable_ComponentInterpolationGet",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_ComponentInterpolationGet

  !
  !================================================================================================================================
  !

  !>Gets the label for a field variable component for character labels. \see OpenCMISS::Iron::cmfe_Field_ComponentLabelGet
  SUBROUTINE FieldVariable_ComponentLabelGetC(fieldVariable,componentNumber,label,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the label for
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The field component number to get the label for
    CHARACTER(LEN=*), INTENT(OUT) :: label !<On return, the field variable label
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: cLength,vsLength
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_ComponentLabelGetC",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(fieldVariable%components)) &
      & CALL FlagError("Field variable components has not been allocated.",err,error,*999)
    IF(componentNumber<1.OR.componentNumber>fieldVariable%numberOfComponents) THEN
      localError="The specified component number of "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " is invalid for variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) &
        & localError=localError//" of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//". The number of components should be >= 1 and <= "// &
        & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    cLength=LEN(label)
    vsLength=LEN_TRIM(fieldVariable%components(componentNumber)%componentLabel)
    IF(cLength>vsLength) THEN
      label=CHAR(LEN_TRIM(fieldVariable%components(componentNumber)%componentLabel))
    ELSE
      label=CHAR(fieldVariable%components(componentNumber)%componentLabel,cLength)
    ENDIF

    EXITS("FieldVariable_ComponentLabelGetC")
    RETURN
999 ERRORSEXITS("FieldVariable_ComponentLabelGetC",err,error)
    RETURN 1

  END SUBROUTINE FieldVariable_ComponentLabelGetC

  !
  !================================================================================================================================
  !

  !>Gets the label for a field variable component for varying string labels.
  SUBROUTINE FieldVariable_ComponentLabelGetVS(fieldVariable,componentNumber,label,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the label for
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The field component number to get the label for
    TYPE(VARYING_STRING), INTENT(OUT) :: label !<On return, the field variable label
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_ComponentLabelGetVS",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(fieldVariable%components)) &
      & CALL FlagError("Field variable components has not been allocated.",err,error,*999)
    IF(componentNumber<1.OR.componentNumber>fieldVariable%numberOfComponents) THEN
      localError="The specified component number of "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " is invalid for variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) &
        & localError=localError//" of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//". The number of components should be >= 1 and <= "// &
        & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    label=fieldVariable%components(componentNumber)%componentLabel
 
    EXITS("FieldVariable_ComponentLabelGetVS")
    RETURN
999 ERRORSEXITS("FieldVariable_ComponentLabelGetVS",err,error)
    RETURN 1

  END SUBROUTINE FieldVariable_ComponentLabelGetVS

  !
  !================================================================================================================================
  !

  !>Check the mesh component number for a field variable component.
  SUBROUTINE FieldVariable_ComponentMeshComponentCheck(fieldVariable,componentNumber,meshComponent,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to check the mesh component for
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The field component number to check the field variable component for
    INTEGER(INTG), INTENT(IN) :: meshComponent !<The mesh component to check for the specified field variable component
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_ComponentMeshComponentCheck",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    
    IF(.NOT.ALLOCATED(fieldVariable%components)) THEN
      localError="The components are not allocated for variable type "// &
        & TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) &
        localError=localError//" of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(componentNumber<1.OR.componentNumber>fieldVariable%numberOfComponents) THEN
      localError="The specified omponent number of "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " is invalid for variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) &
        & localError=localError//" of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//". The component number should be >= 1 and <= "// &
        & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    IF(fieldVariable%components(componentNumber)%meshComponentNumber/=meshComponent) THEN
      localError="Invalid mesh component number. The mesh component number for component number "// &
        & TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) &
        & localError=localError//" of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//" is "// &
        & TRIM(NumberToVString(fieldVariable%components(componentNumber)%interpolationType,"*",err,error))// &
        & " which is does correspond to the specified mesh component number of "// &
        & TRIM(NumberToVString(meshComponent,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("FieldVariable_ComponentMeshComponentCheck")
    RETURN
999 ERRORSEXITS("FieldVariable_ComponentMeshComponentCheck",err,error)
    RETURN 1

  END SUBROUTINE FieldVariable_ComponentMeshComponentCheck

  !
  !================================================================================================================================
  !

  !>Returns the mesh component number of the specified field variable component
  SUBROUTINE FieldVariable_ComponentMeshComponentGet(fieldVariable,componentIdx,meshComponentNumber,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the mesh component number for
    INTEGER(INTG), INTENT(IN) :: componentIdx !<The component index of the field variable to get the mesh component number for. 
    INTEGER(INTG), INTENT(OUT) :: meshComponentNumber  !<On exit, the mesh component number for field variable component. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_ComponentMeshComponentGet",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(fieldVariable%components)) THEN
      localError="The components are not allocated for variable type "// &
        & TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) &
        localError=localError//" of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(componentIdx<1.OR.componentIdx>fieldVariable%numberOfComponents) THEN
      localError="The specified component number of "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
        & " is invalid for variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) &
        & localError=localError//" of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//". The number of components should be >= 1 and <= "// &
        & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
     
    meshComponentNumber=fieldVariable%components(componentIdx)%meshComponentNumber

    EXITS("FieldVariable_ComponentMeshComponentGet")
    RETURN
999 ERRORS("FieldVariable_ComponentMeshComponentGet",err,error)
    EXITS("FieldVariable_ComponentMeshComponentGet")
    RETURN 1
    
  END SUBROUTINE FieldVariable_ComponentMeshComponentGet

  !
  !================================================================================================================================
  !

  !>Checks the data type for a field variable.
  SUBROUTINE FieldVariable_DataTypeCheck(fieldVariable,dataType,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to check the data type for
    INTEGER(INTG), INTENT(IN) :: dataType !<The data type of the field variable to check \see FieldRoutines_DataTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_DataTypeCheck",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    
    SELECT CASE(dataType)              
    CASE(FIELD_INTG_TYPE)
      IF(fieldVariable%dataType/=FIELD_INTG_TYPE) THEN
        localError="Invalid data type. The data type for variable type "// &
          & TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
        IF(ASSOCIATED(fieldVariable%field)) &
          & localError=localError//" of field number "// &
          & TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
        localError=localError//" is "//TRIM(NumberToVString(fieldVariable%dataType,"*",err,error))// &
          & " which is not an integer data type."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(FIELD_SP_TYPE)
      IF(fieldVariable%dataType/=FIELD_SP_TYPE) THEN
        localError="Invalid data type. The data type for variable type "// &
          & TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
        IF(ASSOCIATED(fieldVariable%field)) &
          & localError=localError//" of field number "// &
          & TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
        localError=localError//" is "//TRIM(NumberToVString(fieldVariable%dataType,"*",err,error))// &
          & " which is not a single precision data type."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(FIELD_DP_TYPE)
      IF(fieldVariable%dataType/=FIELD_DP_TYPE) THEN
        localError="Invalid data type. The data type for variable type "// &
          & TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
        IF(ASSOCIATED(fieldVariable%field)) &
          & localError=localError//" of field number "// &
          & TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
        localError=localError//" is "//TRIM(NumberToVString(fieldVariable%dataType,"*",err,error))// &
          & " which is not a double precision data type."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(FIELD_L_TYPE)
      IF(fieldVariable%dataType/=FIELD_L_TYPE) THEN
        localError="Invalid data type. The data type for variable type "// &
          & TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
        IF(ASSOCIATED(fieldVariable%field)) &
          & localError=localError//" of field number "// &
          & TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
        localError=localError//" is "//TRIM(NumberToVString(fieldVariable%dataType,"*",err,error))// &
          & " which is not a logical data type."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The specified data type of "//TRIM(NumberToVString(dataType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("FieldVariable_DataTypeCheck")
    RETURN
999 ERRORSEXITS("FieldVariable_DataTypeCheck",err,error)
    RETURN 1

  END SUBROUTINE FieldVariable_DataTypeCheck

  !
  !================================================================================================================================
  !

  !>Gets the data type for a field variable. 
  SUBROUTINE FieldVariable_DataTypeGet(fieldVariable,dataType,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the data type for
     INTEGER(INTG), INTENT(OUT) :: dataType !<On return, the data type of the field variable \see FieldRoutines_DataTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("FieldVariable_DataTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)

    dataType=fieldVariable%dataType

    EXITS("FieldVariable_DataTypeGet")
    RETURN
999 ERRORSEXITS("FieldVariable_DataTypeGet",err,error)
    RETURN 1

  END SUBROUTINE FieldVariable_DataTypeGet

  !
  !================================================================================================================================
  !

  !>Checks the dimension for a field variable.
  SUBROUTINE FieldVariable_DimensionCheck(fieldVariable,dimensionType,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to check the dimension for
    INTEGER(INTG), INTENT(IN) :: dimensionType !<The field dimension to check \see FieldRoutines_DimensionTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_DimensionCheck",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    
    SELECT CASE(dimensionType)
    CASE(FIELD_SCALAR_DIMENSION_TYPE)
      IF(fieldVariable%dimension/=FIELD_SCALAR_DIMENSION_TYPE) THEN
        localError="Invalid dimension type. The dimension type for variable type "// &
          & TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
        IF(ASSOCIATED(fieldVariable%field)) &
          & localError=localError//" of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
        localError=localError//" is "//TRIM(NumberToVString(fieldVariable%dimension,"*",err,error))// &
          & " which is not a scalar field."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(FIELD_VECTOR_DIMENSION_TYPE)
      IF(fieldVariable%dimension/=FIELD_VECTOR_DIMENSION_TYPE) THEN
        localError="Invalid dimension type. The dimension type for variable type "// &
          & TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
        IF(ASSOCIATED(fieldVariable%field)) &
          & localError=localError//" of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
        localError=localError//" is "//TRIM(NumberToVString(fieldVariable%dimension,"*",err,error))// &
          & " which is not a vector field."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(FIELD_TENSOR_DIMENSION_TYPE) 
      IF(fieldVariable%dimension/=FIELD_TENSOR_DIMENSION_TYPE) THEN
        localError="Invalid dimension type. The dimension type for variable type "// &
          & TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
        IF(ASSOCIATED(fieldVariable%field)) &
          & localError=localError//" of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
        localError=localError//" is "//TRIM(NumberToVString(fieldVariable%dimension,"*",err,error))// &
         & " which is not a tensor field."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The specified dimension type of "//TRIM(NumberToVString(dimensionType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("FieldVariable_DimensionCheck")
    RETURN
999 ERRORSEXITS("FieldVariable_DimensionCheck",err,error)
    RETURN 1

  END SUBROUTINE FieldVariable_DimensionCheck

  !
  !================================================================================================================================
  !

  !>Gets the dimension for a field variable. 
  SUBROUTINE FieldVariable_DimensionGet(fieldVariable,dimension,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the dimension for
    INTEGER(INTG), INTENT(OUT) :: dimension !<On return, the field dimension to get \see FieldRoutines_DimensionTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("FieldVariable_DimensionGet",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)

    dimension=fieldVariable%dimension

    EXITS("FieldVariable_DimensionGet")
    RETURN
999 ERRORSEXITS("FieldVariable_DimensionGet",err,error)
    RETURN 1

  END SUBROUTINE FieldVariable_DimensionGet

  !
  !================================================================================================================================
  !

  !>Checks the DOF order type for a field variable.
  SUBROUTINE FieldVariable_DOFOrderTypeCheck(fieldVariable,dofOrderType,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to check the DOF order type for
    INTEGER(INTG), INTENT(IN) :: dofOrderType !<The DOF order type of the field variable to check \see FieldRoutines_DOFOrderTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_DOFOrderTypeCheck",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    
    SELECT CASE(dofOrderType)              
    CASE(FIELD_SEPARATED_COMPONENT_DOF_ORDER)
      IF(fieldVariable%dofOrderType/=FIELD_SEPARATED_COMPONENT_DOF_ORDER) THEN
        localError="Invalid DOF order type. The DOF order type for variable type "// &
          & TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
        IF(ASSOCIATED(fieldVariable%field)) &
          & localError=localError//" of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
        localError=localError//" is "//TRIM(NumberToVString(fieldVariable%dataType,"*",err,error))// &
          & " which is not a separated component DOF order type."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(FIELD_CONTIGUOUS_COMPONENT_DOF_ORDER)
      IF(fieldVariable%dofOrderType/=FIELD_CONTIGUOUS_COMPONENT_DOF_ORDER) THEN
        localError="Invalid DOF order type. The DOF order type for variable type "// &
          & TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
        IF(ASSOCIATED(fieldVariable%field)) &
          & localError=localError//" of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
        localError=localError//" is "//TRIM(NumberToVString(fieldVariable%dataType,"*",err,error))// &
          & " which is not a contiguous component DOF order type."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The specified DOF order type of "//TRIM(NumberToVString(dofOrderType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("FieldVariable_DOFOrderTypeCheck")
    RETURN
999 ERRORSEXITS("FieldVariable_DOFOrderTypeCheck",err,error)
    RETURN 1

  END SUBROUTINE FieldVariable_DOFOrderTypeCheck

  !
  !================================================================================================================================
  !

  !>Gets the DOF order type for a field variable.
  SUBROUTINE FieldVariable_DOFOrderTypeGet(fieldVariable,dofOrderType,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the DOF order type for
    INTEGER(INTG), INTENT(OUT) :: dofOrderType !<On return, the DOF order type of the field variable \see FieldRoutines_DOFOrderTypes,FieldRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("FieldVariable_DOFOrderTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)

    dofOrderType=fieldVariable%dofOrderType

    EXITS("FieldVariable_DOFOrderTypeGet")
    RETURN
999 ERRORSEXITS("FieldVariable_DOFOrderTypeGet",err,error)
    RETURN 1

  END SUBROUTINE FieldVariable_DOFOrderTypeGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the domain of the specified field variable component
  SUBROUTINE FieldVariable_DomainGet(fieldVariable,componentIdx,domain,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the domain for
    INTEGER(INTG), INTENT(IN) :: componentIdx !<The component index of the field variable to get the domain for. If 0 then the component used to decompose the domain is used. 
    TYPE(DomainType), POINTER :: domain  !<On exit, a pointer to domain for the field variable component. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: domainMeshComponent
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_DomainGet",err,error,*998)

    IF(ASSOCIATED(domain)) CALL FlagError("Domain is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(componentIdx == 0) THEN
      IF(.NOT.ASSOCIATED(fieldVariable%field)) CALL FlagError("Field variable field is not associated.",err,error,*999)
      IF(.NOT.ASSOCIATED(fieldVariable%field%decomposition)) &
        & CALL FlagError("Field variable field is not associated.",err,error,*999)
      IF(.NOT.ALLOCATED(fieldVariable%field%decomposition%domain)) &
        & CALL FlagError("Decomposition domain is not allocated.",err,error,*999)
      domainMeshComponent=fieldVariable%field%decomposition%meshComponentNumber
      IF(domainMeshComponent<1.OR.domainMeshComponent>fieldVariable%field%decomposition%numberOfComponents) THEN
        localError="The domain mesh component of "//TRIM(NumberToVString(domainMeshComponent,"*",err,error))// &
          & " is invalid. The mesh component must be >= 1 and <= "// &
          & TRIM(NumberToVString(fieldVariable%field%decomposition%numberOfComponents,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      domain=>fieldVariable%field%decomposition%domain(domainMeshComponent)%ptr
    ELSE
      IF(componentIdx<1.OR.componentIdx>fieldVariable%numberOfComponents) THEN
        localError="The specified field variable component of "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
          & " is invalid. The field variable component must be >= 1 and <= "// &
          & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(.NOT.ALLOCATED(fieldVariable%components)) &
        & CALL FlagError("Field variable components has not been allocated.",err,error,*999)
      domain=>fieldVariable%components(componentIdx)%domain
    ENDIF
    IF(.NOT.ASSOCIATED(domain)) THEN
      localError="The field variable domain is not associated for component number "// &
        & TRIM(NumberToVString(componentIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("FieldVariable_DomainGet")
    RETURN
999 NULLIFY(domain)
998 ERRORSEXITS("FieldVariable_DomainGet",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_DomainGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the domain mapping of the specified field variable
  SUBROUTINE FieldVariable_DomainMappingGet(fieldVariable,domainMapping,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the domain mapping for
    TYPE(DomainMappingType), POINTER :: domainMapping  !<On exit, a pointer to domain mapping for the field variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_DomainMappingGet",err,error,*998)

    IF(ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    
    domainMapping=>fieldVariable%domainMapping
    IF(.NOT.ASSOCIATED(domainMapping)) THEN
      localError="Domain mapping is not associated"// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    EXITS("FieldVariable_DomainMappingGet")
    RETURN
999 NULLIFY(domainMapping)
998 ERRORSEXITS("FieldVariable_DomainMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_DomainMappingGet

  !
  !================================================================================================================================
  !

  !>Returns the DOF number for a constant DOF
  SUBROUTINE FieldVariable_ConstantDOFGet(fieldVariable,componentNumber,dofNumber,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the DOF for
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(OUT) :: dofNumber  !<On exit, the dof number for field variable component. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_ConstantDOFGet",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    IF(componentNumber<1.OR.componentNumber>fieldVariable%numberOfComponents) THEN
      localError="The specified field variable component number of "// &
        & TRIM(NumberToVString(componentNumber,"*",err,error))//" is invalid. The component number should be >= 1 and <= "// &
        & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))// &
        & " for field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components)) THEN
      localError="Components is not allocated"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(fieldVariable%components(componentNumber)%interpolationType/=FIELD_CONSTANT_INTERPOLATION) THEN
      localError="Cannot index by constant"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      SELECT CASE(fieldVariable%components(componentNumber)%interpolationType)
      CASE(FIELD_CONSTANT_INTERPOLATION)
        !?
      CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
        localError=localError//" which has element based interpolation."            
      CASE(FIELD_NODE_BASED_INTERPOLATION)
        localError=localError//" which has node based interpolation."
      CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
        localError=localError//" which has grid point based interpolation."
      CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
        localError=localError//" which has Gauss point based interpolation."
      CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
        localError=localError//" which has data point based interpolation."
      CASE DEFAULT
        localError=localError//" which has unknown interpolation."
      END SELECT
      CALL FlagError(localError,err,error,*999)
    ENDIF      

    dofNumber=fieldVariable%components(componentNumber)%paramToDOFMap%constantParam2DOFMap

    EXITS("FieldVariable_ConstantDOFGet")
    RETURN
999 ERRORSEXITS("FieldVariable_ConstantDOFGet",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_ConstantDOFGet

  !
  !================================================================================================================================
  !

  !>Returns the DOF number for a local element DOF
  SUBROUTINE FieldVariable_LocalElementDOFGet(fieldVariable,localElementNumber,componentNumber,dofNumber,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the DOF for
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The local element number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(OUT) :: dofNumber  !<On exit, the dof number for field variable component. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DomainType), POINTER :: domain
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_LocalElementDOFGet",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    IF(componentNumber<1.OR.componentNumber>fieldVariable%numberOfComponents) THEN
      localError="The specified field variable component number of "// &
        & TRIM(NumberToVString(componentNumber,"*",err,error))//" is invalid. The component number should be >= 1 and <= "// &
        & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))// &
        & " for field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components)) THEN
      localError="Components is not allocated"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(fieldVariable%components(componentNumber)%interpolationType/=FIELD_ELEMENT_BASED_INTERPOLATION) THEN
      localError="Cannot index by element"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      SELECT CASE(fieldVariable%components(componentNumber)%interpolationType)
      CASE(FIELD_CONSTANT_INTERPOLATION)
        localError=localError//" which has constant interpolation."            
      CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
        !?
      CASE(FIELD_NODE_BASED_INTERPOLATION)
        localError=localError//" which has node based interpolation."
      CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
        localError=localError//" which has grid point based interpolation."
      CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
        localError=localError//" which has Gauss point based interpolation."
      CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
        localError=localError//" which has data point based interpolation."
      CASE DEFAULT
        localError=localError//" which has unknown interpolation."
      END SELECT
      CALL FlagError(localError,err,error,*999)
    ENDIF      
    domain=>fieldVariable%components(componentNumber)%domain
    IF(.NOT.ASSOCIATED(domain)) THEN
      localError="Domain is not associated"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(decomposition)
    CALL Domain_DecompositionGet(domain,decomposition,err,error,*999)
    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)    
    IF(localElementNumber<1.OR.localElementNumber>decompositionElements%totalNumberOfElements) THEN
      localError="The specified local element number of "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " does not exist in the domain"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
       IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
       localError=localError//" which has "//TRIM(NumberToVString(decompositionElements%totalNumberOfElements,"*",err,error))// &
         & " local elements."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components(componentNumber)%paramToDOFMap%elementParam2DOFMap%elements)) THEN
      localError="Element parameter to dof map elements is not allocated "// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
       IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
         & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
       localError=localError//"."
       CALL FlagError(localError,err,error,*999)
    ENDIF

    dofNumber=fieldVariable%components(componentNumber)%paramToDOFMap%elementParam2DOFMap%elements(localElementNumber)

    EXITS("FieldVariable_LocalElementDOFGet")
    RETURN
999 ERRORSEXITS("FieldVariable_LocalElementDOFGet",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_LocalElementDOFGet

  !
  !================================================================================================================================
  !

  !>Returns the DOF number for a user element DOF
  SUBROUTINE FieldVariable_UserElementDOFGet(fieldVariable,userElementNumber,componentNumber,dofNumber,ghostDof,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the DOF for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(OUT) :: dofNumber  !<On exit, the dof number for field variable component. 
    LOGICAL, INTENT(OUT) :: ghostDOF !<On exit, is .TRUE. if the dof corresponds to a ghost node, .FALSE. if not.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: localElementNumber
    LOGICAL :: userElementExists
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DomainType), POINTER :: domain
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_UserElementDOFGet",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    IF(componentNumber<1.OR.componentNumber>fieldVariable%numberOfComponents) THEN
      localError="The specified field variable component number of "// &
        & TRIM(NumberToVString(componentNumber,"*",err,error))//" is invalid. The component number should be >= 1 and <= "// &
        & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))// &
        & " for field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components)) THEN
      localError="Components is not allocated"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(fieldVariable%components(componentNumber)%interpolationType/=FIELD_ELEMENT_BASED_INTERPOLATION) THEN
      localError="Cannot index by element"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      SELECT CASE(fieldVariable%components(componentNumber)%interpolationType)
      CASE(FIELD_CONSTANT_INTERPOLATION)
        localError=localError//" which has constant interpolation."            
      CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
        !?
      CASE(FIELD_NODE_BASED_INTERPOLATION)
        localError=localError//" which has node based interpolation."
      CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
        localError=localError//" which has grid point based interpolation."
      CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
        localError=localError//" which has Gauss point based interpolation."
      CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
        localError=localError//" which has data point based interpolation."
      CASE DEFAULT
        localError=localError//" which has unknown interpolation."
      END SELECT
      CALL FlagError(localError,err,error,*999)
    ENDIF      
    domain=>fieldVariable%components(componentNumber)%domain
    IF(.NOT.ASSOCIATED(domain)) THEN
      localError="Domain is not associated"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(decomposition)
    CALL Domain_DecompositionGet(domain,decomposition,err,error,*999)
    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)    
    CALL DecompositionElements_ElementCheckExists(decompositionElements,userElementNumber,userElementExists,localElementNumber, &
      & ghostDOF,err,error,*999)
    IF(.NOT.userElementExists) THEN
      localError="The specified user element number of "//TRIM(NumberToVString(userElementNumber,"*",err,error))// &
        & " does not exist in the domain"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localElementNumber<1.OR.localElementNumber>decompositionElements%totalNumberOfElements) THEN
      localError="The specified local element number of "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " does not exist in the domain"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
       IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
       localError=localError//" which has "//TRIM(NumberToVString(decompositionElements%totalNumberOfElements,"*",err,error))// &
         & " local elements."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components(componentNumber)%paramToDOFMap%elementParam2DOFMap%elements)) THEN
      localError="Element parameter to dof map elements is not allocated "// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
       IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
         & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
       localError=localError//"."
       CALL FlagError(localError,err,error,*999)
    ENDIF

    dofNumber=fieldVariable%components(componentNumber)%paramToDOFMap%elementParam2DOFMap%elements(localElementNumber)

    EXITS("FieldVariable_UserElementDOFGet")
    RETURN
999 ERRORSEXITS("FieldVariable_UserElementDOFGet",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_UserElementDOFGet

  !
  !================================================================================================================================
  !

  !>Returns the DOF number for a local node DOF
  SUBROUTINE FieldVariable_LocalNodeDOFGet(fieldVariable,versionNumber,derivativeNumber,localNodeNumber,componentNumber, &
    & dofNumber,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the DOF for
    INTEGER(INTG), INTENT(IN) :: versionNumber !<The version number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(IN) :: derivativeNumber !<The derivative number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(IN) :: localNodeNumber !<The local node number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(OUT) :: dofNumber  !<On exit, the dof number for field variable component. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_LocalNodeDOFGet",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    IF(componentNumber<1.OR.componentNumber>fieldVariable%numberOfComponents) THEN
      localError="The specified field variable component number of "// &
        & TRIM(NumberToVString(componentNumber,"*",err,error))//" is invalid. The component number should be >= 1 and <= "// &
        & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))// &
        & " for field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components)) THEN
      localError="Components is not allocated"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(fieldVariable%components(componentNumber)%interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
      localError="Cannot index by node"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      SELECT CASE(fieldVariable%components(componentNumber)%interpolationType)
      CASE(FIELD_CONSTANT_INTERPOLATION)
        localError=localError//" which has constant interpolation."            
      CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
        localError=localError//" which has element based interpolation."
      CASE(FIELD_NODE_BASED_INTERPOLATION)
        !?
      CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
        localError=localError//" which has grid point based interpolation."
      CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
        localError=localError//" which has Gauss point based interpolation."
      CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
        localError=localError//" which has data point based interpolation."
      CASE DEFAULT
        localError=localError//" which has unknown interpolation."
      END SELECT
      CALL FlagError(localError,err,error,*999)
    ENDIF      
    domain=>fieldVariable%components(componentNumber)%domain
    IF(.NOT.ASSOCIATED(domain)) THEN
      localError="Domain is not associated"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(domainTopology)
    CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(domainNodes)
    CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)    
    IF(localNodeNumber<1.OR.localNodeNumber>domainNodes%totalNumberOfNodes) THEN
      localError="The specified local node number of "//TRIM(NumberToVString(localNodeNumber,"*",err,error))// &
        & " does not exist in the domain"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
       IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//" which has "//TRIM(NumberToVString(domainNodes%totalNumberOfNodes,"*",err,error))//" local nodes."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainNodes%nodes)) THEN
      localError="Domain nodes nodes is not allocated "// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
       IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
         & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
       localError=localError//"."
       CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(derivativeNumber<1.OR.derivativeNumber>domainNodes%nodes(localNodeNumber)%numberOfDerivatives) THEN
      localError="Derivative number "//TRIM(NumberToVString(derivativeNumber,"*",err,error))//" is invalid"// &
        & " for local node number "//TRIM(NumberToVString(localNodeNumber,"*",err,error))// &
        & " of field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//" which has "// &
        & TRIM(NumberToVString(domainNodes%nodes(localNodeNumber)%numberOfDerivatives,"*",err,error))//" derivatives."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainNodes%nodes(localNodeNumber)%derivatives)) THEN
      localError="Domain nodes derivatives is not allocated "// &
        & " for local node number "//TRIM(NumberToVString(localNodeNumber,"*",err,error))// &
        & " of field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(versionNumber<1.OR.versionNumber>domainNodes%nodes(localNodeNumber)%derivatives(derivativeNumber)%numberOfVersions) THEN
      localError="Version number "//TRIM(NumberToVString(versionNumber,"*",err,error))//" is invalid"// &
        & " for derivative number "//TRIM(NumberToVString(derivativeNumber,"*",err,error))// &
        & " of local node number "//TRIM(NumberToVString(localNodeNumber,"*",err,error))// &
        & " of field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldvariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//" which has a maximum of "//TRIM(NumberToVString(domainNodes%nodes(localNodeNumber)% &
        & derivatives(derivativeNumber)%numberOfVersions,"*",err,error))//" versions "// &
        & "(NOTE: version numbers are indexed directly from the value the user specifies during "// &
        & "element creation and no record is kept of the total number of versions the user sets."// &
        & "The maximum version number the user sets defines the total number of versions allocated)."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components(componentNumber)%paramToDOFMap%nodeParam2DOFMap%nodes)) THEN
      localError="Node parameter to dof map nodes is not allocated "// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
       IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
         & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
       localError=localError//"."
       CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components(componentNumber)%paramToDOFMap%nodeParam2DOFMap%nodes(localNodeNumber)% &
      & derivatives)) THEN
      localError="Node parameter to dof map nodes derivatives is not allocated "// &
        & " for local node number "//TRIM(NumberToVString(localNodeNumber,"*",err,error))// &
        & " of field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
       IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
         & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
       localError=localError//"."
       CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components(componentNumber)%paramToDOFMap%nodeParam2DOFMap%nodes(localNodeNumber)% &
      & derivatives(derivativeNumber)%versions)) THEN
      localError="Node parameter to dof map nodes derivatives versions is not allocated "// &
        & " for derivative number "//TRIM(NumberToVString(derivativeNumber,"*",err,error))// &
        & " of local node number "//TRIM(NumberToVString(localNodeNumber,"*",err,error))// &
        & " of field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
       IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
         & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
       localError=localError//"."
       CALL FlagError(localError,err,error,*999)
    ENDIF

    dofNumber=fieldVariable%components(componentNumber)%paramToDOFMap%nodeParam2DOFMap%nodes(localNodeNumber)% &
      & derivatives(derivativeNumber)%versions(versionNumber)

    EXITS("FieldVariable_LocalNodeDOFGet")
    RETURN
999 ERRORSEXITS("FieldVariable_LocalNodeDOFGet",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_LocalNodeDOFGet

  !
  !================================================================================================================================
  !

  !>Returns the DOF number for a user node DOF
  SUBROUTINE FieldVariable_UserNodeDOFGet(fieldVariable,versionNumber,derivativeNumber,userNodeNumber,componentNumber, &
    & dofNumber,ghostDof,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the DOF for
    INTEGER(INTG), INTENT(IN) :: versionNumber !<The version number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(IN) :: derivativeNumber !<The derivative number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(OUT) :: dofNumber  !<On exit, the dof number for field variable component.
    LOGICAL, INTENT(OUT) :: ghostDOF !<On exit, is .TRUE. if the dof corresponds to a ghost node, .FALSE. if not.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: localNodeNumber
    LOGICAL :: userNodeExists
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_UserNodeDOFGet",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    IF(componentNumber<1.OR.componentNumber>fieldVariable%numberOfComponents) THEN
      localError="The specified field variable component number of "// &
        & TRIM(NumberToVString(componentNumber,"*",err,error))//" is invalid. The component number should be >= 1 and <= "// &
        & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))// &
        & " for field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components)) THEN
      localError="Components is not allocated for component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(fieldVariable%components(componentNumber)%interpolationType/=FIELD_NODE_BASED_INTERPOLATION) THEN
      localError="Cannot index by node for component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      SELECT CASE(fieldVariable%components(componentNumber)%interpolationType)
      CASE(FIELD_CONSTANT_INTERPOLATION)
        localError=localError//" which has constant interpolation."            
      CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
        localError=localError//" which has element based interpolation."
      CASE(FIELD_NODE_BASED_INTERPOLATION)
        !?
      CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
        localError=localError//" which has grid point based interpolation."
      CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
        localError=localError//" which has Gauss point based interpolation."
      CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
        localError=localError//" which has data point based interpolation."
      CASE DEFAULT
        localError=localError//" which has unknown interpolation."
      END SELECT
      CALL FlagError(localError,err,error,*999)
    ENDIF      
    domain=>fieldVariable%components(componentNumber)%domain
    IF(.NOT.ASSOCIATED(domain)) THEN
      localError="Domain is not associated"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(domainTopology)
    CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(domainNodes)
    CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)    
    CALL DomainNodes_NodeCheckExists(domainNodes,userNodeNumber,userNodeExists,localNodeNumber,ghostDOF,err,error,*999)
    IF(.NOT.userNodeExists) THEN
      localError="The specified user node number of "//TRIM(NumberToVString(userNodeNumber,"*",err,error))// &
        & " does not exist in the domain"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localNodeNumber<1.OR.localNodeNumber>domainNodes%totalNumberOfNodes) THEN
      localError="The specified local node number of "//TRIM(NumberToVString(localNodeNumber,"*",err,error))// &
        & " does not exist in the domain"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
       IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//" which has "//TRIM(NumberToVString(domainNodes%totalNumberOfNodes,"*",err,error))//" local nodes."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainNodes%nodes)) THEN
      localError="Domain nodes nodes is not allocated "// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
       IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
         & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
       localError=localError//"."
       CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(derivativeNumber<1.OR.derivativeNumber>domainNodes%nodes(localNodeNumber)%numberOfDerivatives) THEN
      localError="Derivative number "//TRIM(NumberToVString(derivativeNumber,"*",err,error))//" is invalid"// &
        & " for local node number "//TRIM(NumberToVString(localNodeNumber,"*",err,error))// &
        & " of field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//" which has "// &
        & TRIM(NumberToVString(domainNodes%nodes(localNodeNumber)%numberOfDerivatives,"*",err,error))//" derivatives."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(domainNodes%nodes(localNodeNumber)%derivatives)) THEN
      localError="Domain nodes derivatives is not allocated "// &
        & " for local node number "//TRIM(NumberToVString(localNodeNumber,"*",err,error))// &
        & " of field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(versionNumber<1.OR.versionNumber>domainNodes%nodes(localNodeNumber)%derivatives(derivativeNumber)%numberOfVersions) THEN
      localError="Version number "//TRIM(NumberToVString(versionNumber,"*",err,error))//" is invalid"// &
        & " for derivative number "//TRIM(NumberToVString(derivativeNumber,"*",err,error))// &
        & " of local node number "//TRIM(NumberToVString(localNodeNumber,"*",err,error))// &
        & " of field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldvariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//" which has a maximum of "//TRIM(NumberToVString(domainNodes%nodes(localNodeNumber)% &
        & derivatives(derivativeNumber)%numberOfVersions,"*",err,error))//" versions "// &
        & "(NOTE: version numbers are indexed directly from the value the user specifies during "// &
        & "element creation and no record is kept of the total number of versions the user sets."// &
        & "The maximum version number the user sets defines the total number of versions allocated)."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components(componentNumber)%paramToDOFMap%nodeParam2DOFMap%nodes)) THEN
      localError="Node parameter to dof map nodes is not allocated "// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
       IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
         & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
       localError=localError//"."
       CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components(componentNumber)%paramToDOFMap%nodeParam2DOFMap%nodes(localNodeNumber)% &
      & derivatives)) THEN
      localError="Node parameter to dof map nodes derivatives is not allocated "// &
        & " for local node number "//TRIM(NumberToVString(localNodeNumber,"*",err,error))// &
        & " of field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
       IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
         & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
       localError=localError//"."
       CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components(componentNumber)%paramToDOFMap%nodeParam2DOFMap%nodes(localNodeNumber)% &
      & derivatives(derivativeNumber)%versions)) THEN
      localError="Node parameter to dof map nodes derivatives versions is not allocated "// &
        & " for derivative number "//TRIM(NumberToVString(derivativeNumber,"*",err,error))// &
        & " of local node number "//TRIM(NumberToVString(localNodeNumber,"*",err,error))// &
        & " of field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
       IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
         & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
       localError=localError//"."
       CALL FlagError(localError,err,error,*999)
    ENDIF

    dofNumber=fieldVariable%components(componentNumber)%paramToDOFMap%nodeParam2DOFMap%nodes(localNodeNumber)% &
      & derivatives(derivativeNumber)%versions(versionNumber)

    EXITS("FieldVariable_UserNodeDOFGet")
    RETURN
999 ERRORSEXITS("FieldVariable_UserNodeDOFGet",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_UserNodeDOFGet

  !
  !================================================================================================================================
  !

  !>Returns the DOF number for a local Gauss DOF
  SUBROUTINE FieldVariable_LocalGaussDOFGet(fieldVariable,gaussPointNumber,localElementNumber,componentNumber,dofNumber, &
    & err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the DOF for
    INTEGER(INTG), INTENT(IN) :: gaussPointNumber !<The Gauss point number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The local element number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(OUT) :: dofNumber  !<On exit, the dof number for field variable component. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DomainType), POINTER :: domain
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_LocalGaussDOFGet",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    IF(componentNumber<1.OR.componentNumber>fieldVariable%numberOfComponents) THEN
      localError="The specified field variable component number of "// &
        & TRIM(NumberToVString(componentNumber,"*",err,error))//" is invalid. The component number should be >= 1 and <= "// &
        & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))// &
        & " for field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components)) THEN
      localError="Components is not allocated"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(fieldVariable%components(componentNumber)%interpolationType/=FIELD_GAUSS_POINT_BASED_INTERPOLATION) THEN
      localError="Cannot index by Gauss point"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      SELECT CASE(fieldVariable%components(componentNumber)%interpolationType)
      CASE(FIELD_CONSTANT_INTERPOLATION)
        localError=localError//" which has constant interpolation."            
      CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
        localError=localError//" which has element based interpolation."            
      CASE(FIELD_NODE_BASED_INTERPOLATION)
        localError=localError//" which has node based interpolation."
      CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
        localError=localError//" which has grid point based interpolation."
      CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
        !?
      CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
        localError=localError//" which has data point based interpolation."
      CASE DEFAULT
        localError=localError//" which has unknown interpolation."
      END SELECT
      CALL FlagError(localError,err,error,*999)
    ENDIF      
    domain=>fieldVariable%components(componentNumber)%domain
    IF(.NOT.ASSOCIATED(domain)) THEN
      localError="Domain is not associated"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(decomposition)
    CALL Domain_DecompositionGet(domain,decomposition,err,error,*999)
    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)    
    IF(localElementNumber<1.OR.localElementNumber>decompositionElements%totalNumberOfElements) THEN
      localError="The specified local element number of "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " does not exist in the domain"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
       IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
       localError=localError//" which has "//TRIM(NumberToVString(decompositionElements%totalNumberOfElements,"*",err,error))// &
         & " local elements."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components(componentNumber)%paramToDOFMap%gaussPointParam2DOFMap%gaussPoints)) THEN
      localError="Gauss point parameter to dof map gauss points is not allocated "// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
         & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(gaussPointNumber<1.OR.gaussPointNumber>SIZE(fieldVariable%components(componentNumber)%paramToDOFMap% &
      & gaussPointParam2DOFMap%gaussPoints,1)) THEN
      localError="The specified Gauss point number of "//TRIM(NumberToVString(gaussPointNumber,"*",err,error))// &
        & " is invalid. The Gauss point number must be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(fieldVariable%components(componentNumber)%paramToDOFMap%gaussPointParam2DOFMap% &
        & gaussPoints,1),"*",err,error))// &
        & " for local element number of "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " of field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    dofNumber=fieldVariable%components(componentNumber)%paramToDOFMap%gaussPointParam2DOFMap% &
      & gaussPoints(gaussPointNumber,localElementNumber)

    EXITS("FieldVariable_LocalGaussDOFGet")
    RETURN
999 ERRORSEXITS("FieldVariable_LocalGaussDOFGet",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_LocalGaussDOFGet

  !
  !================================================================================================================================
  !

  !>Returns the DOF number for a user Gauss point DOF
  SUBROUTINE FieldVariable_UserGaussDOFGet(fieldVariable,gaussPointNumber,userElementNumber,componentNumber,dofNumber,ghostDof, &
    & err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the DOF for
    INTEGER(INTG), INTENT(IN) :: gaussPointNumber !<The Gauss point number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(OUT) :: dofNumber  !<On exit, the dof number for field variable component. 
    LOGICAL, INTENT(OUT) :: ghostDOF !<On exit, is .TRUE. if the dof corresponds to a ghost node, .FALSE. if not.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: localElementNumber
    LOGICAL :: userElementExists
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DomainType), POINTER :: domain
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_UserGaussDOFGet",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    IF(componentNumber<1.OR.componentNumber>fieldVariable%numberOfComponents) THEN
      localError="The specified field variable component number of "// &
        & TRIM(NumberToVString(componentNumber,"*",err,error))//" is invalid. The component number should be >= 1 and <= "// &
        & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))// &
        & " for field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components)) THEN
      localError="Components is not allocated"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(fieldVariable%components(componentNumber)%interpolationType/=FIELD_GAUSS_POINT_BASED_INTERPOLATION) THEN
      localError="Cannot index by Gauss point"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      SELECT CASE(fieldVariable%components(componentNumber)%interpolationType)
      CASE(FIELD_CONSTANT_INTERPOLATION)
        localError=localError//" which has constant interpolation."            
      CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
        localError=localError//" which has element based interpolation."            
      CASE(FIELD_NODE_BASED_INTERPOLATION)
        localError=localError//" which has node based interpolation."
      CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
        localError=localError//" which has grid point based interpolation."
      CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
        !?
      CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
        localError=localError//" which has data point based interpolation."
      CASE DEFAULT
        localError=localError//" which has unknown interpolation."
      END SELECT
      CALL FlagError(localError,err,error,*999)
    ENDIF      
    domain=>fieldVariable%components(componentNumber)%domain
    IF(.NOT.ASSOCIATED(domain)) THEN
      localError="Domain is not associated"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(decomposition)
    CALL Domain_DecompositionGet(domain,decomposition,err,error,*999)
    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)    
    CALL DecompositionElements_ElementCheckExists(decompositionElements,userElementNumber,userElementExists,localElementNumber, &
      & ghostDOF,err,error,*999)
    IF(.NOT.userElementExists) THEN
      localError="The specified user element number of "//TRIM(NumberToVString(userElementNumber,"*",err,error))// &
        & " does not exist in the domain"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localElementNumber<1.OR.localElementNumber>decompositionElements%totalNumberOfElements) THEN
      localError="The specified local element number of "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " does not exist in the domain"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
       IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
       localError=localError//" which has "//TRIM(NumberToVString(decompositionElements%totalNumberOfElements,"*",err,error))// &
         & " local elements."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components(componentNumber)%paramToDOFMap%gaussPointParam2DOFMap%gaussPoints)) THEN
      localError="Gauss point parameter to dof map gauss points is not allocated "// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
         & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(gaussPointNumber<1.OR.gaussPointNumber>SIZE(fieldVariable%components(componentNumber)%paramToDOFMap% &
      & gaussPointParam2DOFMap%gaussPoints,1)) THEN
      localError="The specified Gauss point number of "//TRIM(NumberToVString(gaussPointNumber,"*",err,error))// &
        & " is invalid. The Gauss point number must be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(fieldVariable%components(componentNumber)%paramToDOFMap%gaussPointParam2DOFMap% &
        & gaussPoints,1),"*",err,error))// &
        & " for local element number of "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " of field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    dofNumber=fieldVariable%components(componentNumber)%paramToDOFMap%gaussPointParam2DOFMap% &
      & gaussPoints(gaussPointNumber,localElementNumber)

    EXITS("FieldVariable_UserGaussDOFGet")
    RETURN
999 ERRORSEXITS("FieldVariable_UserGaussDOFGet",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_UserGaussDOFGet

  !
  !================================================================================================================================
  !

  !>Returns the DOF number for a local data point DOF
  SUBROUTINE FieldVariable_LocalDataPointDOFGet(fieldVariable,localDataPointNumber,componentNumber,dofNumber,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the DOF for
    INTEGER(INTG), INTENT(IN) :: localDataPointNumber !<The local data point number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(OUT) :: dofNumber  !<On exit, the dof number for field variable component. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DecompositionDataPointsType), POINTER :: decompositionDataPoints
    TYPE(DomainType), POINTER :: domain
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_LocalDataPointDOFGet",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    IF(componentNumber<1.OR.componentNumber>fieldVariable%numberOfComponents) THEN
      localError="The specified field variable component number of "// &
        & TRIM(NumberToVString(componentNumber,"*",err,error))//" is invalid. The component number should be >= 1 and <= "// &
        & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))// &
        & " for field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components)) THEN
      localError="Components is not allocated"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(fieldVariable%components(componentNumber)%interpolationType/=FIELD_DATA_POINT_BASED_INTERPOLATION) THEN
      localError="Cannot index by data point"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      SELECT CASE(fieldVariable%components(componentNumber)%interpolationType)
      CASE(FIELD_CONSTANT_INTERPOLATION)
        localError=localError//" which has constant interpolation."            
      CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
        localError=localError//" which has element based interpolation."            
      CASE(FIELD_NODE_BASED_INTERPOLATION)
        localError=localError//" which has node based interpolation."
      CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
        localError=localError//" which has grid point based interpolation."
      CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
        localError=localError//" which has Gauss point based interpolation."
      CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
        !?
      CASE DEFAULT
        localError=localError//" which has unknown interpolation."
      END SELECT
      CALL FlagError(localError,err,error,*999)
    ENDIF      
    domain=>fieldVariable%components(componentNumber)%domain
    IF(.NOT.ASSOCIATED(domain)) THEN
      localError="Domain is not associated"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(decomposition)
    CALL Domain_DecompositionGet(domain,decomposition,err,error,*999)
    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(decompositionDataPoints)
    CALL DecompositionTopology_DecompositionDataPointsGet(decompositionTopology,decompositionDataPoints,err,error,*999)    
    IF(localDataPointNumber<1.OR.localDataPointNumber>decompositionDataPoints%totalNumberOfDataPoints) THEN
      localError="The specified local data point number of "//TRIM(NumberToVString(localDataPointNumber,"*",err,error))// &
        & " does not exist in the domain"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
       IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
       localError=localError//" which has "// &
         & TRIM(NumberToVString(decompositionDataPoints%totalNumberOfDataPoints,"*",err,error))// &
         & " local data points."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components(componentNumber)%paramToDOFMap%dataPointParam2DOFMap%dataPoints)) THEN
      localError="Dat point parameter to dof map data points is not allocated "// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
         & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    dofNumber=fieldVariable%components(componentNumber)%paramToDOFMap%dataPointParam2DOFMap%dataPoints(localDataPointNumber)

    EXITS("FieldVariable_LocalDataPointDOFGet")
    RETURN
999 ERRORSEXITS("FieldVariable_LocalDataPointDOFGet",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_LocalDataPointDOFGet

  !
  !================================================================================================================================
  !

  !>Returns the DOF number for a user data point DOF
  SUBROUTINE FieldVariable_UserDataPointDOFGet(fieldVariable,userDataPointNumber,componentNumber,dofNumber,ghostDof,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the DOF for
    INTEGER(INTG), INTENT(IN) :: userDataPointNumber !<The user data point number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(OUT) :: dofNumber  !<On exit, the dof number for field variable component. 
    LOGICAL, INTENT(OUT) :: ghostDOF !<On exit, is .TRUE. if the dof corresponds to a ghost DOF, .FALSE. if not.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: localDataPointNumber
    LOGICAL :: userDataPointExists
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DecompositionDataPointsType), POINTER :: decompositionDataPoints
    TYPE(DomainType), POINTER :: domain
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_UserDataPointDOFGet",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    IF(componentNumber<1.OR.componentNumber>fieldVariable%numberOfComponents) THEN
      localError="The specified field variable component number of "// &
        & TRIM(NumberToVString(componentNumber,"*",err,error))//" is invalid. The component number should be >= 1 and <= "// &
        & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))// &
        & " for field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components)) THEN
      localError="Components is not allocated"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(fieldVariable%components(componentNumber)%interpolationType/=FIELD_DATA_POINT_BASED_INTERPOLATION) THEN
      localError="Cannot index by data point"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      SELECT CASE(fieldVariable%components(componentNumber)%interpolationType)
      CASE(FIELD_CONSTANT_INTERPOLATION)
        localError=localError//" which has constant interpolation."            
      CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
        localError=localError//" which has element based interpolation."            
      CASE(FIELD_NODE_BASED_INTERPOLATION)
        localError=localError//" which has node based interpolation."
      CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
        localError=localError//" which has grid point based interpolation."
      CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
        localError=localError//" which has Gauss point based interpolation."
      CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
        !?
      CASE DEFAULT
        localError=localError//" which has unknown interpolation."
      END SELECT
      CALL FlagError(localError,err,error,*999)
    ENDIF      
    domain=>fieldVariable%components(componentNumber)%domain
    IF(.NOT.ASSOCIATED(domain)) THEN
      localError="Domain is not associated"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(decomposition)
    CALL Domain_DecompositionGet(domain,decomposition,err,error,*999)
    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(decompositionDataPoints)
    CALL DecompositionTopology_DecompositionDataPointsGet(decompositionTopology,decompositionDataPoints,err,error,*999)    
    CALL DecompositionDataPoints_DataPointCheckExists(decompositionDataPoints,userDataPointNumber,userDataPointExists, &
      & localDataPointNumber,ghostDOF,err,error,*999)
    IF(.NOT.userDataPointExists) THEN
      localError="The specified user data point number of "//TRIM(NumberToVString(userDataPointNumber,"*",err,error))// &
        & " does not exist in the domain"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localDataPointNumber<1.OR.localDataPointNumber>decompositionDataPoints%totalNumberOfDataPoints) THEN
      localError="The specified local data point number of "//TRIM(NumberToVString(localDataPointNumber,"*",err,error))// &
        & " does not exist in the domain"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
       IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
       localError=localError//" which has "// &
         & TRIM(NumberToVString(decompositionDataPoints%totalNumberOfDataPoints,"*",err,error))// &
         & " local data points."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components(componentNumber)%paramToDOFMap%dataPointParam2DOFMap%dataPoints)) THEN
      localError="Data point parameter to dof map data points is not allocated "// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
         & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    dofNumber=fieldVariable%COMPONENTS(componentNumber)%paramToDOFMap%dataPointParam2DOFMap%dataPoints(localDataPointNumber)  

    EXITS("FieldVariable_UserDataPointDOFGet")
    RETURN
999 ERRORSEXITS("FieldVariable_UserDataPointDOFGet",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_UserDataPointDOFGet

  !
  !================================================================================================================================
  !

  !>Returns the DOF number for a local element data point DOF
  SUBROUTINE FieldVariable_LocalElementDataDOFGet(fieldVariable,localElementNumber,elementDataPointIndex,componentNumber, &
    & dofNumber,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the DOF for
    INTEGER(INTG), INTENT(IN) :: localElementNumber !<The local element number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(IN) :: elementDataPointIndex !<The index of the data point projected onto the element to get the DOF for. 
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(OUT) :: dofNumber  !<On exit, the dof number for field variable component. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: localDataPointNumber
    TYPE(DataProjectionType), POINTER :: dataProjection
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DecompositionDataPointsType), POINTER :: decompositionDataPoints
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DomainType), POINTER :: domain
    TYPE(FieldType), POINTER :: field
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_LocalElementDataDOFGet",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    IF(componentNumber<1.OR.componentNumber>fieldVariable%numberOfComponents) THEN
      localError="The specified field variable component number of "// &
        & TRIM(NumberToVString(componentNumber,"*",err,error))//" is invalid. The component number should be >= 1 and <= "// &
        & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))// &
        & " for field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components)) THEN
      localError="Components is not allocated"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(fieldVariable%components(componentNumber)%interpolationType/=FIELD_DATA_POINT_BASED_INTERPOLATION) THEN
      localError="Cannot index by data point"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      SELECT CASE(fieldVariable%components(componentNumber)%interpolationType)
      CASE(FIELD_CONSTANT_INTERPOLATION)
        localError=localError//" which has constant interpolation."            
      CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
        localError=localError//" which has element based interpolation."            
      CASE(FIELD_NODE_BASED_INTERPOLATION)
        localError=localError//" which has node based interpolation."
      CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
        localError=localError//" which has grid point based interpolation."
      CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
        localError=localError//" which has Gauss point based interpolation."
      CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
        !?
      CASE DEFAULT
        localError=localError//" which has unknown interpolation."
      END SELECT
      CALL FlagError(localError,err,error,*999)
    ENDIF
    domain=>fieldVariable%components(componentNumber)%domain
    IF(.NOT.ASSOCIATED(domain)) THEN
      localError="Domain is not associated"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(decomposition)
    CALL Domain_DecompositionGet(domain,decomposition,err,error,*999)
    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)    
    IF(localElementNumber<1.OR.localElementNumber>decompositionElements%totalNumberOfElements) THEN
      localError="The specified local element number of "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " does not exist in the domain"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
       IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
       localError=localError//" which has "//TRIM(NumberToVString(decompositionElements%totalNumberOfElements,"*",err,error))// &
         & " local elements."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(field)
    CALL FieldVariable_FieldGet(fieldVariable,field,err,error,*999)
    dataProjection=>field%dataProjection
    IF(.NOT.ASSOCIATED(dataProjection)) THEN
      localError="Data projection is not associated"
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " for field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(decompositionDataPoints)
    CALL DecompositionTopology_DecompositionDataPointsGet(decompositionTopology,decompositionDataPoints,err,error,*999)    
    IF(elementDataPointIndex<1.OR.elementDataPointIndex>decompositionDataPoints%elementDataPoints(localElementNumber)% &
      & numberOfProjectedData) THEN
      localError="The specified element data point index of "//TRIM(NumberToVString(elementDataPointIndex,"*",err,error))// &
        & " does not exist in local element number "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " of field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " for field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//" which has "//TRIM(NumberToVString(decompositionDataPoints%elementDataPoints(localElementNumber)% &
        & numberOfProjectedData,"*",err,error))//" projected data points."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    localDataPointNumber=decompositionDataPoints%elementDataPoints(localElementNumber)%dataIndices(elementDataPointIndex)% &
      & localNumber
    IF(localDataPointNumber<1.OR.localDataPointNumber>decompositionDataPoints%totalNumberOfDataPoints) THEN
      localError="The specified local data point number of "//TRIM(NumberToVString(localDataPointNumber,"*",err,error))// &
        & " does not exist in the domain"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
       IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
       localError=localError//" which has "// &
         & TRIM(NumberToVString(decompositionDataPoints%totalNumberOfDataPoints,"*",err,error))// &
         & " local data points."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components(componentNumber)%paramToDOFMap%dataPointParam2DOFMap%dataPoints)) THEN
      localError="Dat point parameter to dof map data points is not allocated "// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
         & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    dofNumber=fieldVariable%components(componentNumber)%paramToDOFMap%dataPointParam2DOFMap%dataPoints(localDataPointNumber)

    EXITS("FieldVariable_LocalElementDataDOFGet")
    RETURN
999 ERRORSEXITS("FieldVariable_LocalElementDataDOFGet",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_LocalElementDataDOFGet

  !
  !================================================================================================================================
  !

  !>Returns the DOF number for a local user element data point DOF
  SUBROUTINE FieldVariable_UserElementDataDOFGet(fieldVariable,userElementNumber,elementDataPointIndex,componentNumber, &
    & dofNumber,ghostDof,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the DOF for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(IN) :: elementDataPointIndex !<The index of the data point projected onto the element to get the DOF for. 
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number of the field variable to get the DOF for. 
    INTEGER(INTG), INTENT(OUT) :: dofNumber  !<On exit, the dof number for field variable component. 
    LOGICAL, INTENT(OUT) :: ghostDOF !<On exit, is .TRUE. if the dof corresponds to a ghost DOF, .FALSE. if not.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: localDataPointNumber,localElementNumber
    LOGICAL :: userElementExists
    TYPE(DataProjectionType), POINTER :: dataProjection
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DecompositionDataPointsType), POINTER :: decompositionDataPoints
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DomainType), POINTER :: domain
    TYPE(FieldType), POINTER :: field
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_UserElementDataDOFGet",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    IF(componentNumber<1.OR.componentNumber>fieldVariable%numberOfComponents) THEN
      localError="The specified field variable component number of "// &
        & TRIM(NumberToVString(componentNumber,"*",err,error))//" is invalid. The component number should be >= 1 and <= "// &
        & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))// &
        & " for field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components)) THEN
      localError="Components is not allocated"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(fieldVariable%components(componentNumber)%interpolationType/=FIELD_DATA_POINT_BASED_INTERPOLATION) THEN
      localError="Cannot index by data point"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      SELECT CASE(fieldVariable%components(componentNumber)%interpolationType)
      CASE(FIELD_CONSTANT_INTERPOLATION)
        localError=localError//" which has constant interpolation."            
      CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
        localError=localError//" which has element based interpolation."            
      CASE(FIELD_NODE_BASED_INTERPOLATION)
        localError=localError//" which has node based interpolation."
      CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
        localError=localError//" which has grid point based interpolation."
      CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
        localError=localError//" which has Gauss point based interpolation."
      CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
        !?
      CASE DEFAULT
        localError=localError//" which has unknown interpolation."
      END SELECT
      CALL FlagError(localError,err,error,*999)
    ENDIF
    domain=>fieldVariable%components(componentNumber)%domain
    IF(.NOT.ASSOCIATED(domain)) THEN
      localError="Domain is not associated"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(decomposition)
    CALL Domain_DecompositionGet(domain,decomposition,err,error,*999)
    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)    
    CALL DecompositionElements_ElementCheckExists(decompositionElements,userElementNumber,userElementExists,localElementNumber, &
      & ghostDOF,err,error,*999)
    IF(.NOT.userElementExists) THEN
      localError="The specified user element number of "//TRIM(NumberToVString(userElementNumber,"*",err,error))// &
        & " does not exist in the domain"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(localElementNumber<1.OR.localElementNumber>decompositionElements%totalNumberOfElements) THEN
      localError="The specified local element number of "//TRIM(NumberToVString(localElementNumber,"*",err,error))// &
        & " does not exist in the domain"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
       IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
       localError=localError//" which has "//TRIM(NumberToVString(decompositionElements%totalNumberOfElements,"*",err,error))// &
         & " local elements."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(field)
    CALL FieldVariable_FieldGet(fieldVariable,field,err,error,*999)
    dataProjection=>field%dataProjection
    IF(.NOT.ASSOCIATED(dataProjection)) THEN
      localError="Data projection is not associated"
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " for field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(decompositionDataPoints)
    CALL DecompositionTopology_DecompositionDataPointsGet(decompositionTopology,decompositionDataPoints,err,error,*999)    
    IF(elementDataPointIndex<1.OR.elementDataPointIndex>decompositionDataPoints%elementDataPoints(localElementNumber)% &
      & numberOfProjectedData) THEN
      localError="The specified element data point index of "//TRIM(NumberToVString(elementDataPointIndex,"*",err,error))// &
        & " does not exist in local element number "//TRIM(NumberToVString(localElementNumber,"*",err,error))// & 
        & " of field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " for field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//" which has "//TRIM(NumberToVString(decompositionDataPoints%elementDataPoints(localElementNumber)% &
        & numberOfProjectedData,"*",err,error))//" projected data points."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    localDataPointNumber=decompositionDataPoints%elementDataPoints(localElementNumber)%dataIndices(elementDataPointIndex)% &
      & localNumber
    IF(localDataPointNumber<1.OR.localDataPointNumber>decompositionDataPoints%totalNumberOfDataPoints) THEN
      localError="The specified local data point number of "//TRIM(NumberToVString(localDataPointNumber,"*",err,error))// &
        & " does not exist in the domain"// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
       IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
        & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
       localError=localError//" which has "// &
         & TRIM(NumberToVString(decompositionDataPoints%totalNumberOfDataPoints,"*",err,error))// &
         & " local data points."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components(componentNumber)%paramToDOFMap%dataPointParam2DOFMap%dataPoints)) THEN
      localError="Dat point parameter to dof map data points is not allocated "// &
        & " for field component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
        & " of field variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
      IF(ASSOCIATED(fieldVariable%field)) localError=localError// &
         & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    dofNumber=fieldVariable%components(componentNumber)%paramToDOFMap%dataPointParam2DOFMap%dataPoints(localDataPointNumber)

    EXITS("FieldVariable_UserElementDataDOFGet")
    RETURN
999 ERRORSEXITS("FieldVariable_UserElementDataDOFGet",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_UserElementDataDOFGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the field of the specified field variable component
  SUBROUTINE FieldVariable_FieldGet(fieldVariable,field,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the field for
    TYPE(FieldType), POINTER :: field  !<On exit, a pointer to field for the field variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("FieldVariable_FieldGet",err,error,*998)

    IF(ASSOCIATED(field)) CALL FlagError("Field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    
    field=>fieldVariable%field
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("The field variable field is not associated.",err,error,*999)

    EXITS("FieldVariable_FieldGet")
    RETURN
999 NULLIFY(field)
998 ERRORSEXITS("FieldVariable_FieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_FieldGet

  !
  !================================================================================================================================
  !

  !>Gets the label for a field variable for character labels.
  SUBROUTINE FieldVariable_LabelGetC(fieldVariable,label,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the label for
    CHARACTER(LEN=*), INTENT(OUT) :: label !<On return, the field variable label
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: cLength,vsLength

    ENTERS("FieldVariable_LabelGetC",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
     
    cLength=LEN(label)
    vsLength=LEN_TRIM(fieldVariable%variableLabel)
    IF(cLength>vsLength) THEN
      label=CHAR(LEN_TRIM(fieldVariable%variableLabel))
    ELSE
      label=CHAR(fieldVariable%variableLabel,cLength)
    ENDIF

    EXITS("FieldVariable_LabelGetC")
    RETURN
999 ERRORSEXITS("FieldVariable_LabelGetC",err,error)
    RETURN 1

  END SUBROUTINE FieldVariable_LabelGetC

  !
  !================================================================================================================================
  !

  !>Gets the label for a field variable for varying string labels.
  SUBROUTINE FieldVariable_LabelGetVS(fieldVariable,label,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the label for
    TYPE(VARYING_STRING), INTENT(OUT) :: label !<On return, the field variable label
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("FieldVariable_LabelGetVS",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
   
    label=fieldVariable%variableLabel
 
    EXITS("FieldVariable_LabelGetVS")
    RETURN
999 ERRORSEXITS("FieldVariable_LabelGetVS",err,error)
    RETURN 1

  END SUBROUTINE FieldVariable_LabelGetVS

  !
  !================================================================================================================================
  !

  !>Checks the number of field components for a field variable.
  SUBROUTINE FieldVariable_NumberOfComponentsCheck(fieldVariable,numberOfComponents,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to check the number of components
    INTEGER(INTG), INTENT(IN) :: numberOfComponents !The number of components in the field variable to check
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldType), POINTER :: field
    TYPE(FieldCreateValuesCacheType), POINTER :: createValuesCache
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_NumberOfComponentsCheck",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)

    NULLIFY(field)
    CALL FieldVariable_FieldGet(fieldVariable,field,err,error,*999)
    IF(field%fieldFinished) THEN
      IF(fieldVariable%numberOfComponents/=numberOfComponents) THEN
        localError="Invalid number of components. The number components for variable type "// &
          & TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))//" of field number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))//" is "// &
          & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))// &
          & " which does not correspond to the specified number of components of "// &
          & TRIM(NumberToVString(numberOfComponents,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      !Field has not been finished so check the create values cache.
      NULLIFY(createValuesCache)
      CALL Field_CreateValuesCacheGet(field,createValuesCache,err,error,*999)
      IF(.NOT.ALLOCATED(createValuesCache%numberOfComponents)) THEN
        localError="The create values cache number of components is not allocated for field number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(createValuesCache%numberOfComponents(fieldVariable%variableType)/=numberOfComponents) THEN
        localError="Invalid number of components. The number components for variable type "// &
          & TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))//" of field number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))//" is "// &
          & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))// &
          & " which does not correspond to the specified number of components of "// &
          & TRIM(NumberToVString(numberOfComponents,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF

    EXITS("FieldVariable_NumberOfComponentsCheck")
    RETURN
999 ERRORSEXITS("FieldVariable_NumberOfComponentsCheck",err,error)
    RETURN 1

  END SUBROUTINE FieldVariable_NumberOfComponentsCheck

  !
  !================================================================================================================================
  !

  !>Returns the number of components for a field variable
  SUBROUTINE FieldVariable_NumberOfComponentsGet(fieldVariable,numberOfComponents,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the number of components for
    INTEGER(INTG), INTENT(OUT) :: numberOfComponents !<On return, the number of components for the field variable.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("FieldVariable_NumberOfComponentsGet",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
       
    numberOfComponents=fieldVariable%numberOfComponents
    
    EXITS("FieldVariable_NumberOfComponentsGet")
    RETURN
999 ERRORSEXITS("FieldVariable_NumberOfComponentsGet",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_NumberOfComponentsGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to a field parameter set of a specified type if such a parameter set exists
  SUBROUTINE FieldVariable_ParameterSetCheck(fieldVariable,parameterSetType,parameterSet,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to check the parameter set for.
    INTEGER(INTG), INTENT(IN) :: parameterSetType !<The type of parameter set to check. \see FieldRoutines_ParameterSetTypes,FieldRoutines
    TYPE(FieldParameterSetType), POINTER :: parameterSet !<On exit, a pointer to the field parameter set. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_ParameterSetCheck",err,error,*998)

    IF(ASSOCIATED(parameterSet)) CALL FlagError("Field parameter set is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    IF(parameterSetType<0.OR.parameterSetType>FIELD_NUMBER_OF_SET_TYPES) THEN
      localError="The specified parameter set type of "//TRIM(NumberToVString(parameterSetType,"*",err,error))// &
        & " is invalid. The parameter set type must be between 1 and "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_SET_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    parameterSet=>fieldVariable%parameterSets%setType(parameterSetType)%ptr

    EXITS("FieldVariable_ParameterSetCheck")
    RETURN
999 NULLIFY(parameterSet)
998 ERRORSEXITS("FieldVariable_ParameterSetCheck",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_ParameterSetCheck

  !
  !================================================================================================================================
  !

  !>Returns a pointer to a field parameter set of a specified type
  SUBROUTINE FieldVariable_ParameterSetGet(fieldVariable,parameterSetType,parameterSet,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the parameter set for.
    INTEGER(INTG), INTENT(IN) :: parameterSetType !<The type of parameter set to get. \see FieldRoutines_ParameterSetTypes,FieldRoutines
    TYPE(FieldParameterSetType), POINTER :: parameterSet !<On exit, a pointer to the field parameter set. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_ParameterSetGet",err,error,*998)

    CALL FieldVariable_ParameterSetCheck(fieldVariable,parameterSetType,parameterSet,err,error,*999)
    IF(.NOT.ASSOCIATED(parameterSet)) THEN
      IF(ASSOCIATED(fieldVariable%field)) THEN
        localError="The parameter set type of "//TRIM(NumberToVString(parameterSetType,"*",err,error))// &
          & " has not been defined on field number "//TRIM(NumberToVString(fieldVariable%field%userNumber, &
          & "*",err,error))//"."
      ELSE
        localError="The parameter set type of "//TRIM(NumberToVString(parameterSetType,"*",err,error))// &
          & " has not been defined on the field variable."
      ENDIF
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("FieldVariable_ParameterSetGet")
    RETURN
999 NULLIFY(parameterSet)
998 ERRORSEXITS("FieldVariable_ParameterSetGet",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_ParameterSetGet

  !
  !================================================================================================================================
  !

  !>Returns the variable type for a field variable
  SUBROUTINE FieldVariable_VariableTypeGet(fieldVariable,variableType,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the variable type for
    INTEGER(INTG), INTENT(OUT) :: variableType !<On return, the variable type for the field variable.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("FieldVariable_VariableTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
       
    variableType=fieldVariable%variableType
    
    EXITS("FieldVariable_VariableTypeGet")
    RETURN
999 ERRORSEXITS("FieldVariable_VariableTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_VariableTypeGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the region for a fields accounting for region and interfaces
  SUBROUTINE Fields_RegionGet(fields,region,err,error,*)

    !Argument variables
    TYPE(FieldsType), POINTER :: fields !<A pointer to the field variable to get the number of components for
    TYPE(RegionType), POINTER :: region !<On return, a pointer to the region for the fields. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(InterfaceType), POINTER :: interface
    TYPE(VARYING_STRING) :: localError

    ENTERS("Fields_RegionGet",err,error,*998)

    IF(ASSOCIATED(region)) CALL FlagError("Region is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(fields)) CALL FlagError("Fields is not associated.",err,error,*999)
       
    NULLIFY(region)
    NULLIFY(interface)
    region=>fields%region
    IF(.NOT.ASSOCIATED(region)) THEN          
      interface=>fields%interface
      IF(ASSOCIATED(interface)) THEN
        IF(ASSOCIATED(interface%parentRegion)) THEN
          region=>interface%parentRegion     
        ELSE
          localError="The parent region is not associated for interface number "// &
            & TRIM(NumberToVString(interface%userNumber,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("A region or interface is not associated for the fields.",err,error,*999)
      ENDIF
    ENDIF
    
    EXITS("Fields_RegionGet")
    RETURN
999 NULLIFY(region)
998 ERRORSEXITS("Fields_RegionGet",err,error)
    RETURN 1
    
  END SUBROUTINE Fields_RegionGet

  !
  !================================================================================================================================
  !

END MODULE FieldAccessRoutines
