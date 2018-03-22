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
  USE Kinds
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
  INTEGER(INTG), PARAMETER :: FIELD_NEXT_VALUES_SET_TYPE=7 !<The parameter set corresponding to the next field values (at time T+dT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_MEAN_PREDICTED_DISPLACEMENT_SET_TYPE=8 !<The parameter set corresponding to the mean predicited values (at time T+DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_VELOCITY_VALUES_SET_TYPE=9 !<The parameter set corresponding to the velocity values (at time T+DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INITIAL_VELOCITY_SET_TYPE=10 !<The parameter set corresponding to the initial velocity values for dynamic problems. This is also the previous velocity values \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS_VELOCITY_SET_TYPE=10 !<The parameter set corresponding to the previous velocity values (at time T). This is also the initial velocity values for dynamic problems. \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_MEAN_PREDICTED_VELOCITY_SET_TYPE=11 !<The parameter set corresponding to the mean predicited velocity values (at time T+DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_ACCELERATION_VALUES_SET_TYPE=12 !<The parameter set corresponding to the acceleration values (at time T+DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INITIAL_ACCELERATION_SET_TYPE=13 !<The parameter set corresponding to the initial acceleration values for dynamic problems. This is also the previous accelearation values \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS_ACCELERATION_SET_TYPE=13 !<The parameter set corresponding to the previous acceleration values (at time T).This is also the initial acceleration values for dynamic problems. \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_MEAN_PREDICTED_ACCELERATION_SET_TYPE=14 !<The parameter set corresponding to the mean predicted acceleration values (at time T+DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREDICTED_DISPLACEMENT_SET_TYPE=15 !<The parameter set corresponding to the predicted values (at time T+DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREDICTED_VELOCITY_SET_TYPE=16 !<The parameter set corresponding to the predicted velocity values (at time T+DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREDICTED_ACCELERATION_SET_TYPE=17 !<The parameter set corresponding to the predicted acceleration values (at time T+DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_RESIDUAL_SET_TYPE=18 !<The parameter set corresponding to the evaluated residual values (at time T) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS_RESIDUAL_SET_TYPE=19 !<The parameter set corresponding to the residual values evaluated previously (at time T-DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_MESH_DISPLACEMENT_SET_TYPE=20 !<The parameter set corresponding to the mesh displacement values for ALE \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_MESH_VELOCITY_SET_TYPE=21 !<The parameter set corresponding to the mesh velocity values for ALE \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_BOUNDARY_SET_TYPE=22 !<The parameter set corresponding to the mesh velocity values for ALE \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INPUT_DATA1_SET_TYPE=23 !<The parameter set corresponding to a input field \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INPUT_DATA2_SET_TYPE=24 !<The parameter set corresponding to a input field \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INPUT_DATA3_SET_TYPE=25 !<The parameter set corresponding to a input field \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INPUT_VEL1_SET_TYPE=26 !<The parameter set corresponding to a input field (PPE)\see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INPUT_VEL2_SET_TYPE=27 !<The parameter set corresponding to a input field (PPE)\see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INPUT_VEL3_SET_TYPE=28 !<The parameter set corresponding to a input field (PPE)\see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INPUT_LABEL_SET_TYPE=29 !<The parameter set corresponding to a input field (PPE)\see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PRESSURE_VALUES_SET_TYPE=30 !<The parameter set corresponding to the surface pressure values (at time T+DT). \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS_PRESSURE_SET_TYPE=31 !<The parameter set corresponding to the previous surface pressure values (at previous increment step). \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES 
  INTEGER(INTG), PARAMETER :: FIELD_RELATIVE_VELOCITY_SET_TYPE=32 !<The parameter set corresponding to the relative velocity values for ALE \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_NEGATIVE_MESH_VELOCITY_SET_TYPE=33 !<The parameter set corresponding to the NEGATIVE mesh velocity values for ALE \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE=34 !<The parameter set corresponding to the previous iteration field values (at iteration n) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE=35 !<The parameter set corresponding to the impermeable flag field values \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INTEGRATED_NEUMANN_SET_TYPE=36 !<Stores integrated Neumann values calculated from Neumann point values \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_UPWIND_VALUES_SET_TYPE=37 !<Stores upwind values associated with a field. \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS_UPWIND_VALUES_SET_TYPE=38 !<Stores upwind values associated with a field from previous timestep. \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  !>@}

  !Module types

  !Module variables

  !Interfaces
  
  INTERFACE FIELD_COORDINATE_SYSTEM_GET    
    MODULE PROCEDURE Field_CoordinateSystemGet
  END INTERFACE FIELD_COORDINATE_SYSTEM_GET
  
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
    & FIELD_INCREMENTAL_VALUES_SET_TYPE,FIELD_ANALYTIC_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,FIELD_NEXT_VALUES_SET_TYPE, &
    & FIELD_MEAN_PREDICTED_DISPLACEMENT_SET_TYPE,FIELD_VELOCITY_VALUES_SET_TYPE,FIELD_INITIAL_VELOCITY_SET_TYPE, &
    & FIELD_PREVIOUS_VELOCITY_SET_TYPE,FIELD_MEAN_PREDICTED_VELOCITY_SET_TYPE,FIELD_ACCELERATION_VALUES_SET_TYPE, &
    & FIELD_INITIAL_ACCELERATION_SET_TYPE,FIELD_PREVIOUS_ACCELERATION_SET_TYPE,FIELD_MEAN_PREDICTED_ACCELERATION_SET_TYPE, &
    & FIELD_PREDICTED_DISPLACEMENT_SET_TYPE,FIELD_PREDICTED_VELOCITY_SET_TYPE,FIELD_PREDICTED_ACCELERATION_SET_TYPE, &
    & FIELD_RESIDUAL_SET_TYPE,FIELD_PREVIOUS_RESIDUAL_SET_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,FIELD_MESH_VELOCITY_SET_TYPE, &
    & FIELD_BOUNDARY_SET_TYPE,FIELD_INPUT_DATA1_SET_TYPE,FIELD_INPUT_DATA2_SET_TYPE,FIELD_INPUT_DATA3_SET_TYPE, &
    & FIELD_PRESSURE_VALUES_SET_TYPE,FIELD_PREVIOUS_PRESSURE_SET_TYPE,FIELD_RELATIVE_VELOCITY_SET_TYPE, &
    & FIELD_NEGATIVE_MESH_VELOCITY_SET_TYPE,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,FIELD_INPUT_VEL1_SET_TYPE, &
    & FIELD_INPUT_VEL2_SET_TYPE,FIELD_INPUT_VEL3_SET_TYPE,FIELD_INPUT_LABEL_SET_TYPE,FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE, &
    & FIELD_INTEGRATED_NEUMANN_SET_TYPE,FIELD_UPWIND_VALUES_SET_TYPE,FIELD_PREVIOUS_UPWIND_VALUES_SET_TYPE

  PUBLIC Field_ComponentInterpolationGet
 
  PUBLIC Field_CoordinateSystemGet  

  PUBLIC FIELD_COORDINATE_SYSTEM_GET

  PUBLIC Field_DataProjectionGet

  PUBLIC Field_DecompositionGet

  PUBLIC Field_RegionGet
  
  PUBLIC FIELD_REGION_GET

  PUBLIC Field_UserNumberFind

  PUBLIC FIELD_USER_NUMBER_FIND

  PUBLIC Field_VariableGet

  PUBLIC FieldVariable_ComponentInterpolationGet

  PUBLIC FieldVariable_DomainGet

  PUBLIC FieldVariable_FieldGet

  PUBLIC FieldVariable_ParameterSetGet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Gets the interpolation type for a field variable component identified by a pointer. \see OpenCMISS::Iron::cmfe_FieldComponentInterpolationGet
  SUBROUTINE Field_ComponentInterpolationGet(field,variableType,componentNumber,interpolationType,err,error,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: field !<A pointer to the field to get the interpolation type for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type of the field variable component to get \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The field component number of the field variable component to set
    INTEGER(INTG), INTENT(OUT) :: interpolationType !<On return, the interpolation type of the field variable component \see FIELD_ROUTINES_InterpolationTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_ComponentInterpolationGet",err,error,*999)

    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(.NOT.field%FIELD_FINISHED) THEN
      localError="Field number "//TRIM(NumberToVString(field%USER_NUMBER,"*",err,error))// &
        & " has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
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
    TYPE(FIELD_TYPE), POINTER :: field !<A pointer to the field to get the coordinate system for
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem!<On return, the field coordinate system. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(INTERFACE_TYPE), POINTER :: interface
    TYPE(REGION_TYPE), POINTER :: region
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_CoordinateSystemGet",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is already associated.",err,error,*999)

    NULLIFY(coordinateSystem)
    NULLIFY(interface)
    region=>field%region
    IF(ASSOCIATED(region)) THEN
      coordinateSystem=>region%COORDINATE_SYSTEM
      IF(.NOT.ASSOCIATED(coordinateSystem)) THEN
        localError="The coordinate system is not associated for field number "// &
          & TRIM(NumberToVString(field%USER_NUMBER,"*",err,error))//" of region number "// &
          & TRIM(NumberToVString(region%USER_NUMBER,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      interface=>field%interface
      IF(ASSOCIATED(interface)) THEN
        coordinateSystem=>interface%COORDINATE_SYSTEM
        IF(.NOT.ASSOCIATED(coordinateSystem)) THEN
          localError="The coordinate system is not associated for field number "// &
            & TRIM(NumberToVString(field%USER_NUMBER,"*",err,error))//" of interface number "// &
            & TRIM(NumberToVString(region%USER_NUMBER,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        localError="A region or interface is not associated for field number "// &
          & TRIM(NumberToVString(field%USER_NUMBER,"*",err,error))//"."
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

  !>Gets a data project from a field.
  SUBROUTINE Field_DataProjectionGet(field,dataProjection,err,error,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: field !<The field to get the data projection for.
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
      localError="Data projection  is not associated for field "//TRIM(NumberToVString(field%USER_NUMBER,"*",err,error))//"."
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

  !>Gets a decomposition from a field.
  SUBROUTINE Field_DecompositionGet(field,decomposition,err,error,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: field !<The field to get the decomposition for.
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition !<On exit, a pointer to the decomposition for the field. Must not be associated on entry.
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
      localError="Field decomposition is not associated for decomposition "// &
        & TRIM(NumberToVString(field%USER_NUMBER,"*",err,error))//"."
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

  !>Returns the region for a field accounting for regions and interfaces
  SUBROUTINE Field_RegionGet(field,region,err,error,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: field !<A pointer to the field to get the region for
    TYPE(REGION_TYPE), POINTER :: region !<On return, the fields region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(INTERFACE_TYPE), POINTER :: interface
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
        IF(ASSOCIATED(interface%PARENT_REGION)) THEN
          region=>interface%PARENT_REGION     
        ELSE
          localError="The parent region is not associated for field number "// &
            & TRIM(NumberToVString(field%USER_NUMBER,"*",err,error))//" of interface number "// &
            & TRIM(NumberToVString(interface%USER_NUMBER,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        localError="A region or interface is not associated for field number "// &
          & TRIM(NumberToVString(field%USER_NUMBER,"*",err,error))//"."
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

  !>Finds and returns in field a pointer to the field identified by a user number in the given list of fields. If no field with that user number exists field is left nullified.
  SUBROUTINE Field_UserNumberFindGeneric(userNumber,fields,field,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The field user number to find
    TYPE(FIELDS_TYPE), POINTER :: fields !<The list of fields containing the field
    TYPE(FIELD_TYPE), POINTER :: field !<On return, a pointer to the field with the given user number. If no field with that user number exists in the list of fields the field is null. Must not be associated on entry.
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
      DO fieldIdx=1,fields%NUMBER_OF_FIELDS
        IF(ASSOCIATED(fields%fields(fieldIdx)%ptr)) THEN
          IF(fields%fields(fieldIdx)%ptr%USER_NUMBER==userNumber) THEN
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
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface containing the field
    TYPE(FIELD_TYPE), POINTER :: field !<On exit, a pointer to the field with the given user number. If no field with that user number exists in the interface the field is null. Must not be associated on entry.
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
    TYPE(REGION_TYPE), POINTER :: region !<A pointer to the region containing the field
    TYPE(FIELD_TYPE), POINTER :: field !<On exit, a pointer to the field with the given user number. If no field with that user number exists in the region the field is null. Must not be associated on entry.
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
  !
  !================================================================================================================================
  !

  !>Returns a pointer to a field variable of a specified type
  SUBROUTINE Field_VariableGet(field,variableType,fieldVariable,err,error,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: field !<A pointer to the field to get the variable for.
    INTEGER(INTG), INTENT(IN) :: variableType !<The type of field variable to set. \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable !<On exit, a pointer to the field variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Field_VariableGet",err,error,*998)

    IF(ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(field)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(.NOT.field%FIELD_FINISHED) THEN
      localError="Field number "//TRIM(NumberToVString(field%USER_NUMBER,"*",err,error))//" has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(variableType<0.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The field variable type must be between 1 and "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    fieldVariable=>field%VARIABLE_TYPE_MAP(variableType)%ptr
    IF(.NOT.ASSOCIATED(fieldVariable)) THEN
      localError="The field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " has not been defined on field number "//TRIM(NumberToVString(field%USER_NUMBER,"*",err,error))//"."
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

  !>Returns the interpolation type of the specified field variable component
  SUBROUTINE FieldVariable_ComponentInterpolationGet(fieldVariable,componentIdx,interpolationType,err,error,*)

    !Argument variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable !<A pointer to the field variable to get the interpolation type for
    INTEGER(INTG), INTENT(IN) :: componentIdx !<The component index of the field variable to get the interpolation type for. 
    INTEGER(INTG), INTENT(OUT) :: interpolationType  !<On exit, the interpolation type for field variable component. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_ComponentInterpolationGet",err,error,*999)

    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(componentIdx<1.OR.componentIdx>fieldVariable%NUMBER_OF_COMPONENTS) THEN
      IF(ASSOCIATED(fieldVariable%field)) THEN
        localError="The specified component number of "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
          & " is invalid for variable type "//TRIM(NumberToVString(fieldVariable%VARIABLE_TYPE,"*",err,error))// &
          & " of field number "//TRIM(NumberToVString(fieldVariable%field%USER_NUMBER,"*",err,error))// &
          & ". The field variable component must be >= 1 and <= "// &
          & TRIM(NumberToVString(fieldVariable%NUMBER_OF_COMPONENTS,"*",err,error))//"."
      ELSE
        localError="The specified component number of "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
          & " is invalid for variable type "//TRIM(NumberToVString(fieldVariable%VARIABLE_TYPE,"*",err,error))// &
          & ". The field variable component must be >= 1 and <= "// &
          & TRIM(NumberToVString(fieldVariable%NUMBER_OF_COMPONENTS,"*",err,error))//"."
      ENDIF
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(fieldVariable%components)) &
      & CALL FlagError("Field variable components has not been allocated.",err,error,*999)
    
    interpolationType=fieldVariable%components(componentIdx)%INTERPOLATION_TYPE

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
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable !<A pointer to the field variable to get the domain for
    INTEGER(INTG), INTENT(IN) :: componentIdx !<The component index of the field variable to get the domain for. If 0 then the component used to decompose the domain is used. 
    TYPE(DOMAIN_TYPE), POINTER :: domain  !<On exit, a pointer to domain for the field variable component. Must not be associated on entry.
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
      IF(.NOT.ASSOCIATED(fieldVariable%field%decomposition%domain)) &
        & CALL FlagError("Decomposition domain is not associated.",err,error,*999)
      domainMeshComponent=fieldVariable%field%decomposition%MESH_COMPONENT_NUMBER
      IF(domainMeshComponent<1.OR.domainMeshComponent>fieldVariable%field%decomposition%numberOfComponents) THEN
        localError="The domain mesh component of "//TRIM(NumberToVString(domainMeshComponent,"*",err,error))// &
          & " is invalid. The mesh component must be >= 1 and <= "// &
          & TRIM(NumberToVString(fieldVariable%field%decomposition%numberOfComponents,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      domain=>fieldVariable%field%decomposition%domain(domainMeshComponent)%ptr
    ELSE
      IF(componentIdx<1.OR.componentIdx>fieldVariable%NUMBER_OF_COMPONENTS) THEN
        localError="The specified field variable component of "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
          & " is invalid. The field variable component must be >= 1 and <= "// &
          & TRIM(NumberToVString(fieldVariable%NUMBER_OF_COMPONENTS,"*",err,error))//"."
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

  !>Returns a pointer to the field of the specified field variable component
  SUBROUTINE FieldVariable_FieldGet(fieldVariable,field,err,error,*)

    !Argument variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable !<A pointer to the field variable to get the field for
    TYPE(FIELD_TYPE), POINTER :: field  !<On exit, a pointer to field for the field variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("FieldVariable_FieldGet",err,error,*998)

    IF(ASSOCIATED(field)) CALL FlagError("Field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field is not associated.",err,error,*999)
    
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

  !>Returns a pointer to a field parameter set of a specified type
  SUBROUTINE FieldVariable_ParameterSetGet(fieldVariable,parameterSetType,parameterSet,err,error,*)

    !Argument variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable !<A pointer to the field variable to get the parameter set for.
    INTEGER(INTG), INTENT(IN) :: parameterSetType !<The type of parameter set to get. \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: parameterSet !<On exit, a pointer to the field parameter set. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("FieldVariable_ParameterSetGet",err,error,*998)

    IF(ASSOCIATED(parameterSet)) CALL FlagError("Field parameter set is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(parameterSetType<0.OR.parameterSetType>FIELD_NUMBER_OF_SET_TYPES) THEN
      localError="The specified parameter set type of "//TRIM(NumberToVString(parameterSetType,"*",err,error))// &
        & " is invalid. The parameter set type must be between 1 and "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_SET_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    parameterSet=>fieldVariable%PARAMETER_SETS%SET_TYPE(parameterSetType)%ptr
    IF(.NOT.ASSOCIATED(parameterSet)) THEN
      IF(ASSOCIATED(fieldVariable%field)) THEN
        localError="The parameter set type of "//TRIM(NumberToVString(parameterSetType,"*",err,error))// &
          & " has not been defined on field number "//TRIM(NumberToVString(fieldVariable%field%USER_NUMBER, &
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

END MODULE FieldAccessRoutines
