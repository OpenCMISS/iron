!> \file
!> \author Richard Christie
!> \brief Test code for FieldML API arguments.
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
!> The Original Code is OpenCMISS-Iron FieldML Arguments Test:
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand. Portions created by the University of Auckland are
!> Copyright (C) 2016 by the University of Auckland. All Rights Reserved.
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

MODULE IRON_TEST_FIELDML_ARGUMENTS

  USE OpenCMISS
  USE OpenCMISS_Iron
  ! USE ISO_C_BINDING
  USE IRON_TEST_FRAMEWORK

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: TestFieldMLArguments

CONTAINS

  SUBROUTINE TestFieldMLArguments(worldRegion)
    TYPE(cmfe_RegionType), INTENT(IN) :: worldRegion

  CALL TestFieldMLInvalidInfo(worldRegion)
  END SUBROUTINE TestFieldMLArguments

  ! Call FieldML routines with initialised but not created fieldmlInfo
  SUBROUTINE TestFieldMLInvalidInfo(worldRegion)
    TYPE(cmfe_RegionType), INTENT(IN) :: worldRegion
    ! local variables
    TYPE(cmfe_RegionType) :: region
    ! local variables
    TYPE(cmfe_CoordinateSystemType) :: coordinateSystem
    TYPE(cmfe_FieldType) :: geometricField
    TYPE(cmfe_NodesType) :: nodes
    TYPE(cmfe_FieldMLIOType) :: fieldmlInfo
    INTEGER(CMISSIntg) :: originalErrorHandlingMode
    INTEGER(CMISSIntg) :: err

    err = 0
    CALL BEGIN_TEST("invalid info")

    ! Don't want to trap errors; this tests that they fail gracefully
    CALL cmfe_ErrorHandlingModeGet(originalErrorHandlingMode, err)
    CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_OUTPUT_ERROR, err)

    CALL cmfe_FieldMLIO_Initialise(fieldmlInfo, err) 

    CALL cmfe_CoordinateSystem_Initialise(coordinateSystem, err)
    CALL cmfe_FieldML_InputCoordinateSystemCreateStart( fieldmlInfo, "dummy field name", coordinateSystem, &
      & AUTO_USER_NUMBER(), err )
    CALL EXPECT_EQ("cmfe_FieldML_InputCoordinateSystemCreateStart result with null FieldML info argument", 1, err)

    CALL cmfe_Region_Initialise(region, err)
    CALL cmfe_Nodes_Initialise(nodes, err)
    CALL cmfe_FieldML_InputNodesCreateStart(fieldmlInfo, "dummy nodes argument name", region, nodes, err)
    CALL EXPECT_EQ("cmfe_FieldML_InputNodesCreateStart result with null FieldML info argument", 1, err)

    CALL cmfe_Field_Initialise(geometricField, err)
    CALL cmfe_FieldML_OutputAddField(fieldmlInfo, "coordinates", "PLAIN_TEXT", geometricField, &
      & CMFE_FIELD_U_VARIABLE_TYPE, CMFE_FIELD_VALUES_SET_TYPE, err)
    CALL EXPECT_EQ("cmfe_FieldML_OutputAddField result with null FieldML info argument", 1, err)

    CALL cmfe_FieldMLIO_Finalise(fieldmlInfo, err)
    CALL EXPECT_EQ("cmfe_FieldMLIO_Finalise result with null FieldML info argument", 1, err)

    CALL cmfe_ErrorHandlingModeSet(originalErrorHandlingMode, err)

    CALL END_TEST()
  END SUBROUTINE TestFieldMLInvalidInfo

END MODULE IRON_TEST_FIELDML_ARGUMENTS
