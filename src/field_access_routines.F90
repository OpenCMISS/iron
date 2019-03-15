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

  !> \addtogroup FIELD_ROUTINES_VariableTypes FIELD_ROUTINES::VariableTypes
  !> \brief Field variable type parameters
  !> \see FIELD_ROUTINES,OPENCMISS_FieldVariableTypes
  !> \todo sort out variable access routines so that you are always accessing by variable type rather than variable number.
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_NUMBER_OF_VARIABLE_TYPES=49 !<Number of different field variable types possible \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
  INTEGER(INTG), PARAMETER :: FIELD_NUMBER_OF_VARIABLE_SUBTYPES=4 !<The number of variants of a particular variable - currently 4. U, delUdelN, delUdelT,del2UdelT2 \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_U_VARIABLE_TYPE=1 !<Standard variable type i.e., u \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELUDELN_VARIABLE_TYPE=2 !<Normal derivative variable type i.e., du/dn \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELUDELT_VARIABLE_TYPE=3 !<First time derivative variable type i.e., du/dt \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DEL2UDELT2_VARIABLE_TYPE=4 !<Second time derivative variable type i.e., d^2u/dt^2 \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_V_VARIABLE_TYPE=5 !<Second standard variable type i.e., v \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELVDELN_VARIABLE_TYPE=6 !<Second normal derivative variable type i.e., dv/dn \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELVDELT_VARIABLE_TYPE=7 !<First time derivative variable type i.e., du/dt \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DEL2VDELT2_VARIABLE_TYPE=8 !<Second time derivative variable type i.e., d^2u/dt^2 \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_W_VARIABLE_TYPE=9 !<Third standard variable type i.e., w \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_U1_VARIABLE_TYPE=10 !<Third standard variable type i.e., v \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELU1DELN_VARIABLE_TYPE=11 !<Third normal derivative variable type i.e., dv/dn \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELU1DELT_VARIABLE_TYPE=12 !<Third time derivative variable type i.e., du/dt \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DEL2U1DELT2_VARIABLE_TYPE=13 !<Third time derivative variable type i.e., d^2u/dt^2 \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_U2_VARIABLE_TYPE=14 !<Fourth standard variable type i.e., v \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELU2DELN_VARIABLE_TYPE=15 !<Fourth normal derivative variable type i.e., dv/dn \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELU2DELT_VARIABLE_TYPE=16 !<Fourth time derivative variable type i.e., du/dt \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DEL2U2DELT2_VARIABLE_TYPE=17 !<Fourth time derivative variable type i.e., d^2u/dt^2 \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_U3_VARIABLE_TYPE=18 !<Fifth standard variable type i.e., v \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELU3DELN_VARIABLE_TYPE=19 !<Fifth normal derivative variable type i.e., dv/dn \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELU3DELT_VARIABLE_TYPE=20 !<Fifth time derivative variable type i.e., du/dt \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DEL2U3DELT2_VARIABLE_TYPE=21 !<Fifth time derivative variable type i.e., d^2u/dt^2 \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_U4_VARIABLE_TYPE=22 !<Sixth standard variable type i.e., v \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELU4DELN_VARIABLE_TYPE=23 !<Sixth normal derivative variable type i.e., dv/dn \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELU4DELT_VARIABLE_TYPE=24 !<Sixth time derivative variable type i.e., du/dt \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DEL2U4DELT2_VARIABLE_TYPE=25 !<Sixth time derivative variable type i.e., d^2u/dt^2 \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_U5_VARIABLE_TYPE=26 !<Seventh standard variable type i.e., v \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELU5DELN_VARIABLE_TYPE=27 !<Seventh normal derivative variable type i.e., dv/dn \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELU5DELT_VARIABLE_TYPE=28 !<Seventh time derivative variable type i.e., du/dt \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DEL2U5DELT2_VARIABLE_TYPE=29 !<Seventh time derivative variable type i.e., d^2u/dt^2 \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_U6_VARIABLE_TYPE=30 !<Eighth standard variable type i.e., v \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELU6DELN_VARIABLE_TYPE=31 !<Eighth normal derivative variable type i.e., dv/dn \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELU6DELT_VARIABLE_TYPE=32 !<Eighth time derivative variable type i.e., du/dt \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DEL2U6DELT2_VARIABLE_TYPE=33 !<Eighth time derivative variable type i.e., d^2u/dt^2 \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_U7_VARIABLE_TYPE=34 !<Ninth standard variable type i.e., v \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELU7DELN_VARIABLE_TYPE=35 !<Ninth normal derivative variable type i.e., dv/dn \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELU7DELT_VARIABLE_TYPE=36 !<Ninth time derivative variable type i.e., du/dt \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DEL2U7DELT2_VARIABLE_TYPE=37 !<Ninth time derivative variable type i.e., d^2u/dt^2 \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_U8_VARIABLE_TYPE=38 !<Tenth standard variable type i.e., v \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELU8DELN_VARIABLE_TYPE=39 !<Tenth normal derivative variable type i.e., dv/dn \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELU8DELT_VARIABLE_TYPE=40 !<Tenth time derivative variable type i.e., du/dt \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DEL2U8DELT2_VARIABLE_TYPE=41 !<Tenth time derivative variable type i.e., d^2u/dt^2 \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_U9_VARIABLE_TYPE=42 !<Eleventh standard variable type i.e., v \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELU9DELN_VARIABLE_TYPE=43 !<Eleventh normal derivative variable type i.e., dv/dn \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELU9DELT_VARIABLE_TYPE=44 !<Eleventh time derivative variable type i.e., du/dt \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DEL2U9DELT2_VARIABLE_TYPE=45 !<Eleventh time derivative variable type i.e., d^2u/dt^2 \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_U10_VARIABLE_TYPE=46 !<Twelfth standard variable type i.e., v \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELU10DELN_VARIABLE_TYPE=47 !<Twelfth normal derivative variable type i.e., dv/dn \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELU10DELT_VARIABLE_TYPE=48 !<Twelfth time derivative variable type i.e., du/dt \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DEL2U10DELT2_VARIABLE_TYPE=49 !<Twelfth time derivative variable type i.e., d^2u/dt^2 \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  !>@}
  
  !> \addtogroup FIELD_ROUTINES_ParameterSetTypes FIELD_ROUTINES::ParameterSetTypes
  !> \brief Field parameter set type parameters \todo make program defined constants negative?
  !> \see FIELD_ROUTINES,OPENCMISS_FieldParameterSetTypes
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_NUMBER_OF_SET_TYPES=99 !<The maximum number of different parameter sets for a field \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_VALUES_SET_TYPE=1 !<The parameter set corresponding to the field values (at time T+DT for dynamic problems) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_BOUNDARY_CONDITIONS_SET_TYPE=2 !<The parameter set corresponding to the field boundary conditions \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INITIAL_VALUES_SET_TYPE=3 !<The parameter set corresponding to the field initial values \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INCREMENTAL_VALUES_SET_TYPE=4 !<The parameter set corresponding to the field incremental values \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_ANALYTIC_VALUES_SET_TYPE=5 !<The parameter set corresponding to the analytic field values \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS_VALUES_SET_TYPE=6 !<The parameter set corresponding to the previous field values (at time T) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS2_VALUES_SET_TYPE=7 !<The parameter set corresponding to the previous field values (at time T-DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS3_VALUES_SET_TYPE=8 !<The parameter set corresponding to the previous field values (at time T-DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_NEXT_VALUES_SET_TYPE=9 !<The parameter set corresponding to the next field values (at time T+dT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_MEAN_PREDICTED_DISPLACEMENT_SET_TYPE=10 !<The parameter set corresponding to the mean predicited values (at time T+DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_VELOCITY_VALUES_SET_TYPE=11 !<The parameter set corresponding to the velocity values (at time T+DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INITIAL_VELOCITY_SET_TYPE=12 !<The parameter set corresponding to the initial velocity values for dynamic problems. This is also the previous velocity values \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS_VELOCITY_SET_TYPE=12 !<The parameter set corresponding to the previous velocity values (at time T). This is also the initial velocity values for dynamic problems. \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_MEAN_PREDICTED_VELOCITY_SET_TYPE=13 !<The parameter set corresponding to the mean predicited velocity values (at time T+DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_ACCELERATION_VALUES_SET_TYPE=14 !<The parameter set corresponding to the acceleration values (at time T+DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INITIAL_ACCELERATION_SET_TYPE=15 !<The parameter set corresponding to the initial acceleration values for dynamic problems. This is also the previous accelearation values \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS_ACCELERATION_SET_TYPE=15 !<The parameter set corresponding to the previous acceleration values (at time T).This is also the initial acceleration values for dynamic problems. \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_MEAN_PREDICTED_ACCELERATION_SET_TYPE=16 !<The parameter set corresponding to the mean predicted acceleration values (at time T+DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREDICTED_DISPLACEMENT_SET_TYPE=17 !<The parameter set corresponding to the predicted values (at time T+DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREDICTED_VELOCITY_SET_TYPE=18 !<The parameter set corresponding to the predicted velocity values (at time T+DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREDICTED_ACCELERATION_SET_TYPE=19 !<The parameter set corresponding to the predicted acceleration values (at time T+DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_RESIDUAL_SET_TYPE=20 !<The parameter set corresponding to the evaluated residual values (at time T+DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS_RESIDUAL_SET_TYPE=21 !<The parameter set corresponding to the residual values evaluated previously (at time T) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS2_RESIDUAL_SET_TYPE=22 !<The parameter set corresponding to the residual values evaluated previously (at time T-DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS3_RESIDUAL_SET_TYPE=23 !<The parameter set corresponding to the residual values evaluated previously (at time T-2*DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_MESH_DISPLACEMENT_SET_TYPE=24 !<The parameter set corresponding to the mesh displacement values for ALE \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_MESH_VELOCITY_SET_TYPE=25 !<The parameter set corresponding to the mesh velocity values for ALE \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_BOUNDARY_SET_TYPE=26 !<The parameter set corresponding to the mesh velocity values for ALE \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INPUT_DATA1_SET_TYPE=27 !<The parameter set corresponding to a input field \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INPUT_DATA2_SET_TYPE=28 !<The parameter set corresponding to a input field \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INPUT_DATA3_SET_TYPE=29 !<The parameter set corresponding to a input field \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INPUT_VEL1_SET_TYPE=30 !<The parameter set corresponding to a input field (PPE)\see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INPUT_VEL2_SET_TYPE=31 !<The parameter set corresponding to a input field (PPE)\see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INPUT_VEL3_SET_TYPE=32 !<The parameter set corresponding to a input field (PPE)\see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INPUT_LABEL_SET_TYPE=33 !<The parameter set corresponding to a input field (PPE)\see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PRESSURE_VALUES_SET_TYPE=34 !<The parameter set corresponding to the surface pressure values (at time T+DT). \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS_PRESSURE_SET_TYPE=35 !<The parameter set corresponding to the previous surface pressure values (at previous increment step). \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES 
  INTEGER(INTG), PARAMETER :: FIELD_RELATIVE_VELOCITY_SET_TYPE=36 !<The parameter set corresponding to the relative velocity values for ALE \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_NEGATIVE_MESH_VELOCITY_SET_TYPE=37 !<The parameter set corresponding to the NEGATIVE mesh velocity values for ALE \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE=38 !<The parameter set corresponding to the previous iteration field values (at iteration n) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE=39 !<The parameter set corresponding to the impermeable flag field values \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INTEGRATED_NEUMANN_SET_TYPE=40 !<Stores integrated Neumann values calculated from Neumann point values \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_UPWIND_VALUES_SET_TYPE=41 !<Stores upwind values associated with a field. \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS_UPWIND_VALUES_SET_TYPE=42 !<Stores upwind values associated with a field from previous timestep. \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  !>@}

  !> \addtogroup FIELD_ROUTINES_InterpolationTypes FIELD_ROUTINES::InterpolationTypes
  !> \brief Field interpolation parameters
  !> \see FIELD_ROUTINES,OPENCMISS_FieldInterpolationTypes
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_CONSTANT_INTERPOLATION=1 !<Constant interpolation. One parameter for the field \see FIELD_ROUTINES_InterpolationTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_ELEMENT_BASED_INTERPOLATION=2 !<Element based interpolation. Parameters are different in each element \see FIELD_ROUTINES_InterpolationTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_NODE_BASED_INTERPOLATION=3 !<Node based interpolation. Parameters are nodal based and a basis function is used \see FIELD_ROUTINES_InterpolationTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_GRID_POINT_BASED_INTERPOLATION=4 !<Grid point based interpolation. Parameters are different at each grid point \see FIELD_ROUTINES_InterpolationTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_GAUSS_POINT_BASED_INTERPOLATION=5 !<Gauss point based interpolation. Parameters are different at each Gauss point \see FIELD_ROUTINES_InterpolationTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DATA_POINT_BASED_INTERPOLATION=6 !<data point based interpolation. Parameters are different at each data point \see FIELD_ROUTINES_InterpolationTypes,FIELD_ROUTINES
  !>@}
   
  !> \addtogroup FIELD_ROUTINES_DataTypes FIELD_ROUTINES::DataTypes
  !> \brief Field data types
  !> \see FIELD_ROUTINES,OPENCMISS_FieldDataTypes
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_INTG_TYPE=1 !<Integer field data type \see FIELD_ROUTINES_DataTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_SP_TYPE=2 !<Single precision real field data type \see FIELD_ROUTINES_DataTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DP_TYPE=3 !<Double precision real field data type \see FIELD_ROUTINES_DataTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_L_TYPE=4 !<Logical field data type \see FIELD_ROUTINES_DataTypes,FIELD_ROUTINES
  !>@}
   
  
  !Module types

  !Module variables
  
  !Interfaces
  
  INTERFACE FIELD_COORDINATE_SYSTEM_GET    
    MODULE PROCEDURE Field_CoordinateSystemGet
  END INTERFACE FIELD_COORDINATE_SYSTEM_GET

  INTERFACE FIELD_GEOMETRIC_FIELD_GET
    MODULE PROCEDURE Field_GeometricFieldGet
  END INTERFACE FIELD_GEOMETRIC_FIELD_GET

  INTERFACE FIELD_MESH_DECOMPOSITION_GET
    MODULE PROCEDURE Field_DecompositionGet
  END INTERFACE FIELD_MESH_DECOMPOSITION_GET

  INTERFACE FIELD_NUMBER_OF_VARIABLES_GET
    MODULE PROCEDURE Field_NumberOfVariablesGet
  END INTERFACE FIELD_NUMBER_OF_VARIABLES_GET

  INTERFACE FIELD_PARAMETER_SET_VECTOR_GET
    MODULE PROCEDURE FieldParameterSet_ParametersGet
  END INTERFACE FIELD_PARAMETER_SET_VECTOR_GET
  
  INTERFACE Field_ParameterSetVectorGet
    MODULE PROCEDURE FieldParameterSet_ParametersGet
  END INTERFACE Field_ParameterSetVectorGet

  INTERFACE FIELD_REGION_GET
    MODULE PROCEDURE Field_RegionGet
  END INTERFACE FIELD_REGION_GET
  
  !>Finds and returns a field identified by a user number. If no field  with that number exists field is left nullified.
  INTERFACE Field_UserNumberFind
    MODULE PROCEDURE Field_UserNumberFindInterface
    MODULE PROCEDURE Field_UserNumberFindRegion
  END INTERFACE Field_UserNumberFind
  
  !>Finds and returns a field identified by a user number. If no field  with that number exists field is left nullified.
  INTERFACE FIELD_USER_NUMBER_FIND
    MODULE PROCEDURE Field_UserNumberFindInterface
    MODULE PROCEDURE Field_UserNumberFindRegion
  END INTERFACE FIELD_USER_NUMBER_FIND

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
  
  PUBLIC Field_ComponentInterpolationGet
 
  PUBLIC Field_CoordinateSystemGet  

  PUBLIC FIELD_COORDINATE_SYSTEM_GET

  PUBLIC Field_CreateValuesCacheGet

  PUBLIC Field_DataProjectionGet

  PUBLIC Field_DecompositionGet

  PUBLIC Field_GeometricFieldGet

  PUBLIC FIELD_GEOMETRIC_FIELD_GET

  PUBLIC FIELD_MESH_DECOMPOSITION_GET

  PUBLIC Field_NumberOfComponentsGet

  PUBLIC Field_NumberOfVariablesGet

  PUBLIC FIELD_NUMBER_OF_VARIABLES_GET

  PUBLIC Field_RegionGet
  
  PUBLIC FIELD_REGION_GET

  PUBLIC Field_ScaleFactorsVectorGet

  PUBLIC Field_UserNumberFind

  PUBLIC FIELD_USER_NUMBER_FIND

  PUBLIC Field_VariableGet

  PUBLIC FieldInterpolatedPoint_InterpolationParametersGet

  PUBLIC FieldInterpolatedPointMetrics_InterpolatedPointGet

  PUBLIC FieldInterpolationParameters_FieldVariableGet

  PUBLIC FieldParameterSet_ParametersGet

  PUBLIC FieldVariable_AssertIsINTGData,FieldVariable_AssertIsSPData,FieldVariable_AssertIsDPData,FieldVariable_AssertIsLData
  
  PUBLIC FieldVariable_AssertComponentNumberOK
  
  PUBLIC FieldVariable_ComponentDomainGet

  PUBLIC FieldVariable_ComponentInterpolationGet

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

  PUBLIC FieldVariable_NumberOfComponentsGet

  PUBLIC FieldVariable_ParameterSetCheck

  PUBLIC FieldVariable_ParameterSetGet
  
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

  !>Gets the interpolation type for a field variable component identified by a pointer. \see OpenCMISS::Iron::cmfe_FieldComponentInterpolationGet
  SUBROUTINE Field_ComponentInterpolationGet(field,variableType,componentNumber,interpolationType,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the interpolation type for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type of the field variable component to get \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The field component number of the field variable component to set
    INTEGER(INTG), INTENT(OUT) :: interpolationType !<On return, the interpolation type of the field variable component \see FIELD_ROUTINES_InterpolationTypes,FIELD_ROUTINES
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

  !>Gets the number of field components for a field variable. \see OpenCMISS::Iron::cmfe_Field_NumberOfComponentsGet
  SUBROUTINE Field_NumberOfComponentsGet(field,variableType,numberOfComponents,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the number of components
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
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

  !>Returns a pointer to a field variable of a specified type
  SUBROUTINE Field_VariableGet(field,variableType,fieldVariable,err,error,*)

    !Argument variables
    TYPE(FieldType), POINTER :: field !<A pointer to the field to get the variable for.
    INTEGER(INTG), INTENT(IN) :: variableType !<The type of field variable to set. \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
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

  !>Returns the interpolation type of the specified field variable component
  SUBROUTINE FieldVariable_ComponentInterpolationGet(fieldVariable,componentIdx,interpolationType,err,error,*)

    !Argument variables
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the interpolation type for
    INTEGER(INTG), INTENT(IN) :: componentIdx !<The component index of the field variable to get the interpolation type for. 
    INTEGER(INTG), INTENT(OUT) :: interpolationType  !<On exit, the interpolation type for field variable component. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_ComponentInterpolationGet",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(componentIdx<1.OR.componentIdx>fieldVariable%numberOfComponents) THEN
      IF(ASSOCIATED(fieldVariable%field)) THEN
        localError="The specified component number of "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
          & " is invalid for variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))// &
          & " of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))// &
          & ". The field variable component must be >= 1 and <= "// &
          & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))//"."
      ELSE
        localError="The specified component number of "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
          & " is invalid for variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))// &
          & ". The field variable component must be >= 1 and <= "// &
          & TRIM(NumberToVString(fieldVariable%numberOfComponents,"*",err,error))//"."
      ENDIF
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components)) &
      & CALL FlagError("Field variable components has not been allocated.",err,error,*999)
    
    interpolationType=fieldVariable%components(componentIdx)%interpolationType

    EXITS("FieldVariable_ComponentInterpolationGet")
    RETURN
999 ERRORSEXITS("FieldVariable_ComponentInterpolationGet",err,error)
    RETURN 1
    
  END SUBROUTINE FieldVariable_ComponentInterpolationGet

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

    dofNumber=fieldVariable%COMPONENTS(componentNumber)%paramToDOFMap%dataPointParam2DOFMap%dataPoints(localDataPointNumber)

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
    INTEGER(INTG), INTENT(IN) :: parameterSetType !<The type of parameter set to check. \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
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
    INTEGER(INTG), INTENT(IN) :: parameterSetType !<The type of parameter set to get. \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
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
