!> \file
!> \author Richard Christie
!> \brief Test code for FieldML cube I/O.
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
!> The Original Code is OpenCMISS-Iron FieldML I/O Cube Test:
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

MODULE IRON_TEST_FIELDML_CUBE

  USE OpenCMISS
  USE OpenCMISS_Iron
  USE ISO_C_BINDING
  USE IRON_TEST_FRAMEWORK

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: TestFieldMLIOCube

CONTAINS

  SUBROUTINE TestFieldMLIOCube(worldRegion)
    TYPE(cmfe_RegionType), INTENT(IN) :: worldRegion
    ! local variables
    TYPE(cmfe_RegionType) :: region
    TYPE(cmfe_MeshType) :: mesh
    TYPE(cmfe_FieldType) :: geometricField
    INTEGER(CMISSIntg) :: err

    err = 0
    CALL BEGIN_TEST("cube")

    ! Initial FieldML file is using Zinc-like naming conventions
    ! Files exported from Zinc are currently not yet readable in Iron due to the
    ! additional derivative/version mappings for element field parameters, which
    ! are present even for Lagrange bases (and set to defaults of 1).
    CALL ReadCube(worldRegion, "input/cube.fieldml", &
      & region, mesh, geometricField, &
      & geometricFieldName = "coordinates", &
      & geometricFieldNodeParametersName = "nodes.coordinates", &
      & nodesArgumentName = "nodes.argument", &
      & basisEvaluatorName = "mesh3d.trilinearLagrange", & ! will become "mesh3d.interpolation1" if exported from Zinc
      & meshArgumentName = "mesh3d.argument", &
      & meshComponentTemplateName = "mesh3d.template1")
    CALL CheckCube(mesh, geometricField)

    CALL WriteCube(mesh, geometricField, "cube", "cube.fieldml")

    CALL cmfe_Field_Destroy(geometricField, err)
    CALL cmfe_Mesh_Destroy(mesh, err)
    CALL cmfe_Region_Destroy(region, err)

    ! File is re-read with Iron naming conventions
    CALL ReadCube(worldRegion, "cube.fieldml", &
      & region, mesh, geometricField, &
      & geometricFieldName = "coordinates", &
      & geometricFieldNodeParametersName = "coordinates.dofs.node", &
      & nodesArgumentName = "cube.nodes.argument", &
      & basisEvaluatorName = "cube.component1trilinearLagrange_3.evaluator", &
      & meshArgumentName = "cube.mesh.argument", &
      & meshComponentTemplateName = "cube.component1.template")
    CALL CheckCube(mesh, geometricField)

    CALL cmfe_Field_Destroy(geometricField, err)
    CALL cmfe_Mesh_Destroy(mesh, err)
    CALL cmfe_Region_Destroy(region, err)

    CALL END_TEST()
  END SUBROUTINE TestFieldMLIOCube

  SUBROUTINE ReadCube(worldRegion, inputFilename, region, mesh, geometricField, &
      & geometricFieldName, geometricFieldNodeParametersName, nodesArgumentName, &
      & basisEvaluatorName, meshArgumentName, meshComponentTemplateName)
    TYPE(cmfe_RegionType), INTENT(IN) :: worldRegion
    CHARACTER(KIND=C_CHAR,LEN=*), INTENT(IN) :: inputFilename
    CHARACTER(KIND=C_CHAR,LEN=*), INTENT(IN) :: geometricFieldName
    CHARACTER(KIND=C_CHAR,LEN=*), INTENT(IN) :: geometricFieldNodeParametersName
    CHARACTER(KIND=C_CHAR,LEN=*), INTENT(IN) :: nodesArgumentName
    CHARACTER(KIND=C_CHAR,LEN=*), INTENT(IN) :: basisEvaluatorName
    CHARACTER(KIND=C_CHAR,LEN=*), INTENT(IN) :: meshArgumentName
    CHARACTER(KIND=C_CHAR,LEN=*), INTENT(IN) :: meshComponentTemplateName
    TYPE(cmfe_RegionType), INTENT(OUT) :: region
    TYPE(cmfe_MeshType), INTENT(OUT) :: mesh
    TYPE(cmfe_FieldType), INTENT(OUT) :: geometricField
    ! local variables
    TYPE(cmfe_BasisType) :: basis
    TYPE(cmfe_CoordinateSystemType) :: coordinateSystem
    TYPE(cmfe_DecompositionType) :: decomposition
    TYPE(cmfe_NodesType) :: nodes
    TYPE(cmfe_FieldMLIOType) :: fieldmlInfo
    INTEGER(CMISSIntg) :: numberOfComputationalNodes, computationalNodeNumber
    INTEGER(CMISSIntg) :: err

    err = 0

    ! Get computational nodes information

    CALL cmfe_ComputationalNumberOfNodesGet(numberOfComputationalNodes, err)
    CALL cmfe_ComputationalNodeNumberGet(computationalNodeNumber, err)

    ! Initialise FieldML and parse input file

    CALL cmfe_FieldMLIO_Initialise(fieldmlInfo, err) 
    CALL cmfe_FieldML_InputCreateFromFile(inputFilename, fieldmlInfo, err)

    ! Define Coordinate System from value type of coordinates field evaluator

    CALL cmfe_CoordinateSystem_Initialise(coordinateSystem, err)
    CALL cmfe_FieldML_InputCoordinateSystemCreateStart( fieldmlInfo, geometricFieldName, coordinateSystem, &
      & AUTO_USER_NUMBER(), err )
    CALL cmfe_CoordinateSystem_CreateFinish(coordinateSystem, err)
    ! CALL cmfe_CoordinateSystem_DimensionGet(coordinateSystem, coordinateCount, err )

    ! Define region

    CALL cmfe_Region_Initialise(region, err)
    CALL cmfe_Region_CreateStart(AUTO_USER_NUMBER(), worldRegion, region, err)
    CALL cmfe_Region_LabelSet(region, "Test", err)
    CALL cmfe_Region_CoordinateSystemSet(region, coordinateSystem, err)
    CALL cmfe_Region_CreateFinish(region, err)

    ! Define nodes from FieldML ensemble argument

    CALL cmfe_Nodes_Initialise(nodes, err)
    CALL cmfe_FieldML_InputNodesCreateStart(fieldmlInfo, nodesArgumentName, region, nodes, err)
    CALL cmfe_Nodes_CreateFinish(nodes, err)

    ! Define bases from FieldML evaluator (referencing interpolator)

    CALL cmfe_Basis_Initialise(basis, err)
    CALL cmfe_FieldML_InputBasisCreateStart(fieldmlInfo, basisEvaluatorName, AUTO_USER_NUMBER(), basis, err)
    CALL cmfe_Basis_CreateFinish(basis, err)

    ! Define mesh from FieldML mesh argument

    CALL cmfe_Mesh_Initialise(mesh, err)
    CALL cmfe_FieldML_InputMeshCreateStart(fieldmlInfo, meshArgumentName, mesh, AUTO_USER_NUMBER(), region, err)
    CALL cmfe_Mesh_NumberOfComponentsSet(mesh, 1, err)
    ! Define Field template
    CALL cmfe_FieldML_InputCreateMeshComponent(fieldmlInfo, mesh, 1, meshComponentTemplateName, err )
    CALL cmfe_Mesh_CreateFinish(mesh, err)

    ! Create Domain decomposition

    CALL cmfe_Decomposition_Initialise(decomposition, err)
    CALL cmfe_Decomposition_CreateStart(AUTO_USER_NUMBER(), mesh, decomposition, err)
    CALL cmfe_Decomposition_TypeSet(decomposition, CMFE_DECOMPOSITION_ALL_TYPE, err)
    CALL cmfe_Decomposition_NumberOfDomainsSet(decomposition, numberOfComputationalNodes, err)
    CALL cmfe_Decomposition_CreateFinish(decomposition, err)

    ! Define Geometric Field

    CALL cmfe_Field_Initialise(geometricField, err)
    CALL cmfe_FieldML_InputFieldCreateStart(fieldmlInfo, region, decomposition, AUTO_USER_NUMBER(), geometricField, &
      & CMFE_FIELD_U_VARIABLE_TYPE, geometricFieldName, err)
    CALL cmfe_Field_CreateFinish( geometricField, err)
    CALL cmfe_FieldML_InputFieldParametersUpdate(fieldmlInfo, geometricField, geometricFieldNodeParametersName, &
      & CMFE_FIELD_U_VARIABLE_TYPE, CMFE_FIELD_VALUES_SET_TYPE, err)
    CALL cmfe_Field_ParameterSetUpdateStart(geometricField, CMFE_FIELD_U_VARIABLE_TYPE, CMFE_FIELD_VALUES_SET_TYPE, err)
    CALL cmfe_Field_ParameterSetUpdateFinish(geometricField, CMFE_FIELD_U_VARIABLE_TYPE, CMFE_FIELD_VALUES_SET_TYPE, err)

    CALL cmfe_FieldMLIO_Finalise(fieldmlInfo, err)
  END SUBROUTINE ReadCube

  SUBROUTINE WriteCube(mesh, geometricField, Basename, OutputFileName)
    TYPE(cmfe_MeshType), INTENT(IN) :: mesh
    TYPE(cmfe_FieldType), INTENT(IN) :: geometricField
    CHARACTER(KIND=C_CHAR,LEN=*), INTENT(IN) :: Basename !< Name prefixed to all objects output, separated by .
    CHARACTER(KIND=C_CHAR,LEN=*), INTENT(IN) :: OutputFileName
    ! parameters
    CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: OutputDirectory = ""
    CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: DataFormat = "PLAIN_TEXT"
    ! local variables
    TYPE(cmfe_FieldMLIOType) :: fieldmlInfo
    INTEGER(CMISSIntg) :: err

    err = 0
    CALL cmfe_FieldMLIO_Initialise(fieldmlInfo, err ) 
    CALL cmfe_FieldML_OutputCreate(mesh, OutputDirectory, Basename, DataFormat, fieldmlInfo, err)
    CALL cmfe_FieldML_OutputAddField(fieldmlInfo, "coordinates", DataFormat, geometricField, &
      & CMFE_FIELD_U_VARIABLE_TYPE, CMFE_FIELD_VALUES_SET_TYPE, err)
    CALL cmfe_FieldML_OutputWrite(fieldmlInfo, OutputFileName, err)
    CALL cmfe_FieldMLIO_Finalise(fieldmlInfo, err)
  END SUBROUTINE WriteCube

  SUBROUTINE CheckCube(mesh, geometricField)
    TYPE(cmfe_MeshType), INTENT(IN) :: mesh
    TYPE(cmfe_FieldType), INTENT(IN) :: geometricField
    ! local variables
    REAL(CMISSRP), PARAMETER, DIMENSION(3, 4) :: xi = &
      & RESHAPE([ 0.5, 0.5, 0.5,  0.0, 0.1, 0.2,  0.0, 0.0, 0.0,  0.667, 0.123, 0.456], [3, 4])
    INTEGER(CMISSIntg), PARAMETER :: elementUserNumber = 1
    REAL(CMISSRP), DIMENSION(3) :: values
    REAL(CMISSRP), PARAMETER :: tolerance = MERGE(1.0E-6_CMISSRP, 1.0E-14_CMISSRP, KIND(0.0_CMISSRP) == KIND(0.0_CMISSSP))
    INTEGER :: i, p
    TYPE(cmfe_MeshNodesType) :: meshNodes
    INTEGER(CMISSIntg) :: elementCount, err, meshComponentCount, nodeCount

    err = 0
    CALL cmfe_Mesh_NumberOfElementsGet(mesh, elementCount, err)
    CALL EXPECT_EQ("element count", 1, elementCount)
    CALL cmfe_Mesh_NumberOfComponentsGet(mesh, meshComponentCount, err)
    CALL EXPECT_EQ("mesh component count", 1, meshComponentCount)
    CALL cmfe_MeshNodes_Initialise(meshNodes, err)
    CALL cmfe_Mesh_NodesGet(mesh, 1, meshNodes, err)
    CALL cmfe_MeshNodes_NumberOfNodesGet(meshNodes, nodeCount, err)
    CALL EXPECT_EQ("node count", 8, nodeCount)
    DO p = 1, 4
      CALL cmfe_Field_ParameterSetInterpolateXi(geometricField, CMFE_FIELD_U_VARIABLE_TYPE, CMFE_FIELD_VALUES_SET_TYPE, &
        & 1, elementUserNumber, xi(:,p), values, err)
      DO i = 1, 3
        CALL EXPECT_NEAR("coordinates", xi(i,p), values(i), tolerance)
      ENDDO
    ENDDO
  END SUBROUTINE CheckCube

END MODULE IRON_TEST_FIELDML_CUBE
