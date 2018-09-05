!> \file
!> \author Richard Christie
!> \brief Test code for Iron FieldML I/O.
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
!> The Original Code is OpenCMISS-Iron FieldML I/O Test:
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

PROGRAM IRON_TEST_FIELDML_IO

  USE OpenCMISS
  USE OpenCMISS_Iron
  USE IRON_TEST_FRAMEWORK
  USE IRON_TEST_FIELDML_CUBE
  USE IRON_TEST_FIELDML_ARGUMENTS

  IMPLICIT NONE

  ! CMISS variables
  TYPE(cmfe_CoordinateSystemType) :: worldCoordinateSystem
  TYPE(cmfe_RegionType) :: worldRegion

  ! Generic CMISS variables

  INTEGER(CMISSIntg) :: numberOfComputationalNodes, computationalNodeNumber
  INTEGER(CMISSIntg) :: err

  CALL INITIALISE_TESTS()

  ! Initialise OpenCMISS-Iron

  CALL cmfe_Initialise(worldCoordinateSystem, worldRegion, err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR, err)

  ! Get computational nodes information

  CALL cmfe_ComputationalNumberOfNodesGet(numberOfComputationalNodes, err)
  CALL cmfe_ComputationalNodeNumberGet(computationalNodeNumber, err)

  CALL TestFieldMLIOCube(worldRegion)
  CALL TestFieldMLArguments(worldRegion)

  CALL cmfe_Finalise(err)

  CALL FINALISE_TESTS(err)

  IF (err /= 0) THEN
    WRITE(*,'(A)') "FieldML I/O tests failed."
    STOP 1
  ENDIF

  WRITE(*,'(A)') "FieldML I/O tests successfully completed."
  STOP

END PROGRAM IRON_TEST_FIELDML_IO
