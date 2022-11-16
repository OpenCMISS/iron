!> \file
!> \author Chris Bradley
!> \brief This module contains all type definitions in order to avoid cyclic module references.
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
!> Contributor(s): Chris Bradley
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
!#### Index: ne
!###  Description:
!###    Index label for a element.
!#### Index: ng
!###  Description:
!###    Index label for a gauss point.
!#### Index: ni
!###  Description:
!###    Index label for a xi direction.
!#### Index: nk
!###  Description:
!###    Index label for a derivative with respect to the global directions.
!#### Index: nn
!###  Description:
!###    Index for a local node within an element.
!#### Index: np
!###  Description:
!###    Index for a node.
!#### Index: ns
!###  Description:
!###    Index for a element parameter within an element.
!#### Index: nu
!###  Description:
!###    Index for a partial derivative.

!> This module contains all type definitions in order to avoid cyclic module references.
MODULE Types

  USE CMISSPETScTypes, ONLY : PetscISColoringType,PetscKspType,PetscMatType,PetscMatColoringType,PetscMatFDColoringType, &
    & PetscPCType,PetscSnesType,PetscSnesLineSearchType,PetscTaoType,PetscVecType
  USE Constants
  USE Kinds
  USE ISO_C_BINDING
  USE ISO_VARYING_STRING
#ifndef NOMPIMOD
  USE MPI
#endif
  USE Trees
  use linkedlist_routines

  IMPLICIT NONE

  PRIVATE

#ifdef NOMPIMOD
#include "mpif.h"
#endif
  
  
  !  
  !================================================================================================================================
  !
  ! Base types
  !

  TYPE RealDPPtrType
    REAL(DP), POINTER :: ptr(:)
  END TYPE RealDPPtrType
  
  TYPE IntegerIntgPtrType
    INTEGER(INTG), POINTER :: ptr(:)
  END TYPE IntegerIntgPtrType

   TYPE IntegerCIntAllocType
    INTEGER(C_INT), ALLOCATABLE :: array(:)
  END TYPE IntegerCIntAllocType

  PUBLIC RealDPPtrType,IntegerIntgPtrType,IntegerCIntAllocType
  
  !  
  !================================================================================================================================
  !
  ! List types
  !

  !>Buffer type to allow arrays of pointers to a list
  TYPE ListPtrType
    TYPE(ListType), POINTER :: ptr !<The pointer to the list
  END TYPE ListPtrType

  !>Contains information on a list
  TYPE ListType
    LOGICAL :: mutable !<If true, list entries can be changed once they have been added. False by default.
    LOGICAL :: listFinished !<Is .TRUE. if the list has finished being created, .FALSE. if not.
    INTEGER(INTG) :: numberInList !<The number of items currently in the list
    INTEGER(INTG) :: dataDimension !<The dimension of the data being stored
    INTEGER(INTG) :: initialSize !<The size of the list when it was initially created.
    INTEGER(INTG) :: size !<The current size of the list.
    INTEGER(INTG) :: dataType !<The data type of the list \see Lists_DataType
    INTEGER(INTG) :: keyDimension !<The key dimension number i.e., the dimension index used for indexing and sorting
    INTEGER(INTG) :: sortOrder !<The ordering to be used when sorting the list \see Lists_SortingOrder
    INTEGER(INTG) :: sortMethod !<The sorting method to be used when sorting the list \see Lists_SortingMethod
    INTEGER(INTG) :: searchMethod !<The search method to be used when searching the list \see Lists_SearchingMethod
    INTEGER(INTG), ALLOCATABLE :: listIntg(:) !<The integer data (dimension = 1) for integer lists. 
    INTEGER(INTG), ALLOCATABLE :: listIntg2(:,:) !<The integer data (dimension > 1) for integer lists. 
    REAL(SP), ALLOCATABLE :: listSP(:) !<The single precision data (dimension = 1)for single precision real lists. 
    REAL(SP), ALLOCATABLE :: listSP2(:,:) !<The single precision data (dimension > 1) for single precision real lists. 
    REAL(DP), ALLOCATABLE :: listDP(:) !<The double precision data (dimension = 1)for double precision real lists. 
    REAL(DP), ALLOCATABLE :: listDP2(:,:) !<The double precision data (dimension > 1) for double precision real lists. 
  END TYPE ListType

  PUBLIC ListPtrType,ListType
    
  !
  !================================================================================================================================
  !
  ! Quadrature types
  !

  !>Contains information for a particular quadrature scheme. \see OpenCMISS::Iron::cmfe_QuadratureSchemeType \todo Also evaluate the product of the basis functions at gauss points for speed???
  TYPE QuadratureSchemeType
    INTEGER(INTG) :: globalNumber !<The global number of the quadrature scheme in the list of quadrature schemes for a particular quadrature.
    TYPE(QuadratureType), POINTER :: quadrature !<The pointer back to the quadrature for a particular quadrature scheme
    INTEGER(INTG) :: numberOfGauss !<The number of gauss points for the quadrature scheme.
    REAL(DP), ALLOCATABLE :: gaussPositions(:,:) !<gaussPositions(xiCoordIdx,gaussPointIdx). The positions in the xiCoordIdx'th xi coordinate of Gauss point gaussPointIdx. Old CMISS name XIG(ni,ng,nb).
    REAL(DP), ALLOCATABLE :: gaussWeights(:) !<gaussWeights(gaussPointIdx). The weight applied to Gauss point gaussPointIdx. Old CMISS name WG(ng,nb).
    REAL(DP), ALLOCATABLE :: gaussBasisFunctions(:,:,:) !<gaussBasisFunctions(elementParameterIdx,partialDerivativeIdx,gaussPointIdx). The value of the basis functions evaluated at Gauss point gaussPointIdx for the partialDerivativeIdx'th derivative of the basis function associated with the elementParameterIdx'th element parameter. Old CMISS name PG(ns,nu,ng,nb)
    !Quadrature information at faces
    INTEGER(INTG), ALLOCATABLE :: numberOfFaceGauss(:) !<numberOfFaceGauss(localFaceIdx) number of gauss points in each local face of the element
    REAL(DP), ALLOCATABLE :: faceGaussBasisFunctions(:,:,:,:) !<faceGaussBasisFunctions(elementParameterIdx,partialDerivativeIdx,gaussPointIdx,localFaceIdx)
    REAL(DP), ALLOCATABLE :: faceGaussPositions(:,:,:) !<faceGaussPositions(xiCoordIdx,gaussPointIdx,localFaceIdx)
    REAL(DP), ALLOCATABLE :: faceGaussWeights(:,:) !<faceGaussWeights(gaussPointIdx,localFaceIdx)
  END TYPE QuadratureSchemeType

  !>A buffer type to allow for an array of pointers to a QuadratureSchemeType \see Types::QuadratureSchemeType
  TYPE QuadratureSchemePtrType
    TYPE(QuadratureSchemeType), POINTER :: ptr !<A pointer to the quadrature scheme
  END TYPE QuadratureSchemePtrType

  !>Contains information on the quadrature to be used for integrating a basis. \see OpenCMISS::Iron::cmfe_QuadratureType
  TYPE QuadratureType
    INTEGER(INTG) :: type !<The type of the quadrature \see BasisRoutines_QuadratureTypes
    TYPE(BasisType), POINTER :: basis !<The pointer back to the basis
    INTEGER(INTG), ALLOCATABLE :: numberOfGaussXi(:) !<numberOfGaussXi(xiIdx). For standard Gauss schemes the number of Gauss points to be used in the xiIdx'th xi direction.
    INTEGER(INTG) :: gaussOrder !<For simplex Gauss schemes the order of the Quadrature scheme i.e., the order/dimension of the polynomial that can be integrated.
    TYPE(QuadratureSchemePtrType), ALLOCATABLE :: quadratureSchemeMap(:) !<quadratureSchemeMap(schemeIdx). The pointer map to the defined quadrature schemes. The size of array is given by BasisRoutines::BASIS_NUMBER_OF_QUADRATURE_SCHEME_TYPES. If the quadrature scheme is not defined for the particular type then the array element is NULL. \see BasisRoutines_QuadratureSchemes.
    INTEGER(INTG) :: numberOfSchemes !<The number of quadrature schemes defined for this quadrature
    TYPE(QuadratureSchemePtrType), POINTER :: schemes(:) !<schemes(schemeIdx). The array of pointers to the quadrature schemes defined for the basis. scheme_idx must be between 1 and QuadratureType::NUMBER_OF_SCHEMES.
    LOGICAL :: evaluateFaceGauss=.FALSE. !! \todo should this be here??
  END TYPE QuadratureType

  PUBLIC QuadratureSchemeType,QuadratureSchemePtrType,QuadratureType

  !
  !================================================================================================================================
  !
  ! Basis types
  !

  !> A buffer type to allow for an array of pointers to a BasisType.
  TYPE BasisPtrType
    TYPE(BasisType), POINTER :: ptr !<The pointer to the basis.
  END TYPE BasisPtrType

  !> Contains all information about a basis .
  TYPE BasisType
    !\todo Add in different sub types for the different types of bases???
    INTEGER(INTG) :: userNumber !<The user defined identifier for the basis. The user number must be unique.
    INTEGER(INTG) :: globalNumber !<The global number for the basis i.e., the position indentifier for the list of bases defined.
    INTEGER(INTG) :: familyNumber !<The family number for the basis. A basis has a number of sub-bases attached which make a basis family. The main parent basis is the basis defined by the user and it will have a family number of 0. The sub-bases of the parent basis will be the line and face bases that make up the basis. These will have different family numbers.
    TYPE(BasisFunctionsType), POINTER :: basisFunctions !<A pointer back to the basis functions for the basis.
    LOGICAL :: basisFinished !<Is .TRUE. if the basis has finished being created, .FALSE. if not.
    LOGICAL :: hermite !<Is .TRUE. if the basis is a hermite basis, .FALSE. if not.
    INTEGER(INTG) :: type !< The type of basis \see BasisRoutines_BasisTypes 
    INTEGER(INTG) :: numberOfXi !<The number of xi directions for the basis.
    INTEGER(INTG) :: numberOfXiCoordinates !<The number of xi coordinate directions for the basis. For Lagrange Hermite tensor product basis functions this is equal to the number of Xi directions. For simplex basis functions this is equal to the number of Xi directions + 1
    INTEGER(INTG), ALLOCATABLE :: interpolationXi(:) !<interpolationXi(xiIdx). The interpolation specification used in the xiIdx'th Xi direction \see BasisRoutines_InterpolationSpecifications
    INTEGER(INTG), ALLOCATABLE :: interpolationType(:) !<interpolationType(xicIdx). The interpolation type in the xiIdx'th Xi coordinate direction. Old CMISS name IBT(1,xiIdx,basisIdx) \see BasisRoutines_InterpolationTypes
    INTEGER(INTG), ALLOCATABLE :: interpolationOrder(:)!<interpolationOrder(xicIdx). The interpolation order in the xiIdx'th Xi coordinate direction. Old CMISS name IBT(2,xiIdx,basisIdx) \see BasisRoutines_InterpolationOrder 
    !Degenerate information
    LOGICAL :: degenerate !<Is .TRUE. if the basis is a degenerate basis (i.e., has collapsed nodes), .FALSE. if not.
    INTEGER(INTG), ALLOCATABLE :: collapsedXi(:) !<collapsedXi(xiIdx). The collapsed state of the xiIdx'th direction. collapsedXi can be either XI_COLLAPSED, COLLAPSED_AT_XI0, COLLAPSED_AT_XI1 or NOT_COLLAPSED depending on whether or not the xiIdx'th direction is collapsed, has a perpendicular Xi collapsed at the xi=0 end of the xiIdx'th direction, has a perpendicular xi collapsed at the xi=1 of the xiIdx'th direction or is not collapsed. NOTE: in old cmiss the value IBT(1,xiIdx) = 5 or 6 was set for the xiIdx that was collapsed. The perpendicular line/face was then stored in IBT(3,xiIdx). For this the quadratic1 and quadratic2 type interpolation types are set on the perpendicular xi direction and the xiIdx direction that is collapsed will have collapsedXi(xiIdx) set to XI_COLLAPSED. BE CAREFUL WITH THIS WHEN TRANSLATING OLD CMISS CODE. Old CMISS name IBT(1,xiIdx) ???? \see BasisRoutines_XiCollapse
    INTEGER(INTG) :: numberOfCollapsedXi !<The number of xi directions in the basis that are collapsed.
    LOGICAL, ALLOCATABLE :: nodeAtCollapse(:) !<nodeAtCollapse(localNodeIdx). Is .TRUE. if the localNodeIdx'th node of the basis is at a collapse, .FALSE. if not.
    !Quadrature
    TYPE(QuadratureType) :: quadrature !<The quadrature schemes for the basis.
    INTEGER(INTG) :: numberOfPartialDerivatives !<The number of partial derivatives for the basis. Old CMISS name NUT(basisFamilyIdx)
    INTEGER(INTG) :: numberOfNodes !<The number of local nodes in the basis. Old CMISS name NNT(basisFamilyIdx)
    !\todo
    INTEGER(INTG), ALLOCATABLE :: numberOfNodesXic(:) !<numberOfNodesXic(xiCoordIdx). The number of local nodes in the xiCoordIdx'th coordinate in the basis. Old CMISS name IBT(2,xiIdx,basisIdx).
    INTEGER(INTG) :: numberOfElementParameters  !<The number of element parameters in the basis. Old CMISS name NST(basisFamilyIdx). 
    INTEGER(INTG) :: maximumNumberOfDerivatives !<The maximum number of derivatives at any node in the basis. Old CMISS name NKT(0,basisFamilyIdx)
    INTEGER(INTG), ALLOCATABLE :: numberOfDerivatives(:) !<numberOfDerivatives(localNodeIdx). The number of derivatives at the localNodeIdx'th node in the basis. Old CMISS name NKT(localNodeIdx,basisFamilyIdx).
    INTEGER(INTG), ALLOCATABLE :: nodePositionIndex(:,:) !<nodePositionIndex(localNodeIdx,xiCoordIdx). The index of the node position for the localNodeIdx'th local node in the xiCoordiIdx'th coordinate. For Lagrange-Hermite tensor product basis functions: The number of coordinates equals the number of xi directions. Thus if nodePositionIndex(localNodeIdx,:)=1,2,2 then local node localNodeIdx is the first node in the xi(coord)=1 direction, the second node in the xi(coord)=2 direction and the second node in the xi(coord)=3 direction; For simplex basis functions: The number of coordinates equals the number of xi directions plus one. The index specifies the inverse distance away from the corner/end of that area coordinate. Thus if an element has quadratic interpolation the index will range from 3 (closest to the corner/end of the element that the area coordinate has the value 1) to 1 (closest to the corners/end of the element that the area coordinate has the value 0). If M is the order of the element then in nodePositionIndex(localNodeIdx,:)=1,1,M then that node is apex for the third area coordinate. In general the index values will add up to M+number of xi directions+1 (i.e., subtract one from the indicies to get the standard simplex coordinates. Old CMISS name INP(localNodeIdx,xiIdx,basisIdx). \see Types::BasisType::nodePositionIndexInv.
    INTEGER(INTG), ALLOCATABLE :: nodePositionIndexInv(:,:,:,:) !<nodePositionIndexInv(localNodeCoord1,localNodeCoord2,localNodeCoord3,localNodeCoord4). The inverse of the node position index for the basis. The nodePositionIndexInv gives the local node number for the node that has node position indices of localNodeCoord1 in the 1st xi(coord) direction, localNodeCoord2 in the 2nd xi(coord) direction, localNodeCoord3 in the 3rd xi(coord) direction and localNodeCoord4 in the 4th xi(coord) direction. NOTE: that if the basis has less than 4 xi(coord) direction the position index is 1. Old CMISS name NNB(inp1,inp2,inp3,nbf). \see Types::BasisType::NodePositionIndex.
    INTEGER(INTG), ALLOCATABLE :: derivativeOrderIndex(:,:,:) !<derivativeOrderIndex(derivativeIdx,localNodeIdx,0:xiIdx). The index of the derivative order for the derivativeIdx'th derivative of the localNodeIdx'th node in the xiIdx'th direction of the basis. The derivative index is NO_PART_DERIV for zeroth order, FIRST_PART_DERIV for the first order and SECOND_PART_DERIV for the second order derivative. Thus a derivativeOrderIndex(derivativeIdx,localNodeIdx,1..) of {NO_PART_DERIV,FIRST_PART_DERIV,NO_PART_DERIV} indicates that the derivativeIdx'th derivative of the localNodeIdx'th node of the basis is the first derivative with respect to the s2 direction. Old CMISS name IDO(derivativeIdx,localNodeIdx,1:xiIdx,basisFamilyIdx). \see Types::BasisType::derivativeOrderIndexInv,CONSTANTS_PartialDerivativeConstants
    INTEGER(INTG), ALLOCATABLE :: derivativeOrderIndexInv(:,:,:,:) !<derivativeOrderIndexInv(partialDerivIdx1,partialDerivIdx2,partialDerivIdx3,localNodeIdx). The inverse of the derivative order index for the nn'th local node of the basis. derivativeOrderIndexInv gives the derivative number for the partialDerivIdx1 partial derivative in the 1st xi direction, the partialDerivIdx2 partial derivative in the 2nd xi direction and the partialDerivIdx3 partial derivative in the 3rd xi direction. NOTE: local node localNodeIdx does not carry any derivatives of the requested partial derivative type then derivativeOrderIndexInv will return 0. If the basis has less than 3 xi directions then the partial derivative index is 1. \see Types::BasisType::derivativeOrderIndex
    INTEGER(INTG), ALLOCATABLE :: partialDerivativeIndex(:,:) !<partialDerivativeIndex(derivativeIdx,localNodeIdx). Gives the partial derivative number (partialDerivIdx) of the derivativeIdx'th derivative of the localNodeIdx'th local node for the basis. Old CMISS name IDO(derivativeIdx,localNodeIdx,0,basisFamilyIdx).
    INTEGER(INTG), ALLOCATABLE :: elementParameterIndex(:,:) !<elementParameterIndex(derivativeIdx,localNodeIdx). Gives the element parameter number (elementParamIdx) of the derivativeIdx'th derivative of the localNodeIdx'th local node for the basis. Old CMISS name NSB(derivativeIdx,localNodeIdx,basisFamilyIdx).
    INTEGER(INTG), ALLOCATABLE :: elementParameterIndexInv(:,:) !<elementParameterIndexInv(1..2,elementParamIdx). Gives the inverse of the element parameter index. elementParameterIndexInv(1,elementParamIdx) gives the local node number corresponding to the elementParamIdx'th element parameter. elementParameterIndexInv(2,elementParamIdx) gives the local derivative number corresponding to the elementParamIdx'th element parameter.
    !Line information
    INTEGER(INTG) :: numberOfLocalLines !<The number of local lines in the basis.
    TYPE(BasisPtrType), ALLOCATABLE :: localLineBasis(:) !<localLineBasis(localLineIdx). localLineBasis(localLineIdx)%ptr is a pointer to the sub-basis used for the localLineIdx'th local line. 
    INTEGER(INTG), ALLOCATABLE :: localLineXiDirection(:) !<localLineXiDirection(localLineIdx). The xi direction of the localLineIdx'th local line for the basis.
    INTEGER(INTG), ALLOCATABLE :: localLineXiNormals(:,:) !<localLineXiNormals(xiCoordIdx,localLineIdx). The xiCoordIdx'th xi directions that are "normal" to the localLineIdx'th local line. There are number of xi coordinate directions - 2 normal xi directions. Not allocated for bases with only 1 xi direction. Note: Normals are always outward.
    INTEGER(INTG), ALLOCATABLE :: xiNormalsLocalLine(:,:) !<xiNormalLocalLine(xiCoordIdx1,xiCoordIdx2). The local line number corresponding to the intersection of the xiCoordIdx1'th and xiCoordIdx2'th normal direction. xiCoordIdx1 and xiCoordIdx2 vary from -numberOfXiCoordinates to + numberOfXiCoordinates. For bases with 2 xi directions xiCoordIdx2 is only 1. Not allocated for bases with only 1 xi direction.
    INTEGER(INTG), ALLOCATABLE :: numberOfNodesInLocalLine(:) !<numberOfNodesInLocalLine(localLineIdx). The the number of nodes in the localLineIdx'th local line for the basis. Old CMISS name NNL(0,localLineIdx,basisIdx).
    INTEGER(INTG), ALLOCATABLE :: nodeNumbersInLocalLine(:,:) !<nodeNumbersInLocalLine(localLineNodeIdx,localLineIdx). The local node numbers (localNodeIdx) for the localLineNodeIdx'th line node in the localLineIdx'th local line for the basis. Old CMISS name NNL(1..,localLineIdx,basisIdx).
    INTEGER(INTG), ALLOCATABLE :: derivativeNumbersInLocalLine(:,:) !<derivativeNumbersInLocalLine(localLineNodeIdx,localLineIdx). The derivative numbers (derivativeIdx) for the localLineNodeIdx'th line node in the localLineIdx'th local line for the basis.
    INTEGER(INTG), ALLOCATABLE :: elementParametersInLocalLine(:,:) !<elementParametersInLocalLine(lineParameterIdx,elementLineIdx). The local element parameter for the lineParameterIdx'th line parameter in the elementLineIdx'th local line for the basis. The number of line parameters in the elementLineIdx'th local line is given by elementParametersInLocalLine(0,elementLineIdx).
    !Face information
    INTEGER(INTG) :: numberOfLocalFaces !<The number of local faces in the basis.
    TYPE(BasisPtrType), ALLOCATABLE :: localFaceBasis(:) !<localFaceBasis(localFaceIdx). localFaceBasis(localFaceIdx)%ptr is a pointer to the sub-basis used for the localFaceIdx'th local face. 
    INTEGER(INTG), ALLOCATABLE :: localFaceXiDirections(:,:) !<localFaceXiDirections(xiCoordIdx,localFaceIdx). The xiCoordIdx'th Xi direction in the localFaceIdx'th local face for the basis. There are numberOfXiCoords - 1 face xi directions. Not allocated for bases with only 2 xi directions. 
    INTEGER(INTG), ALLOCATABLE :: localFaceXiNormal(:) !localFaceXiNormal(localFaceIdx). The xi coordinate direction normal to the localFaceIdx'th local face for the basis. Not allocated for bases with only 2 xi directions. Note: Normals are always outward. \see Constants_ElementNormalXiDirections
    INTEGER(INTG), ALLOCATABLE :: xiNormalLocalFace(:) !xiNormalLocalFace(xiCoordIdx). The local face number of that corresponds to the xiCoordIdx'th normal direction. xiCoordIdx varies from -numberOfXiCoord to + numberOfXiCoord. Not allocated for bases with only 2 xi directions. 
    INTEGER(INTG), ALLOCATABLE :: numberOfNodesInLocalFace(:) !<numberOfNodesInLocalFace(localFaceIdx). The the number of nodes in the localFaceIdx'th local face for the basis. Old CMISS name NNL(0,localFaceIdx,basisIdx).
    INTEGER(INTG), ALLOCATABLE :: nodeNumbersInLocalFace(:,:) !<nodeNumbersInLocalFace(localFaceNodeIdx,localFaceIdx). The local element node numbers (localNodeIdx) for the localFaceNodeIdx'th face node in the localFaceIdx'th local face for the basis. Old CMISS name NNL(1..,localFaceIdx,basisIdx).
    INTEGER(INTG), ALLOCATABLE :: derivativeNumbersInLocalFace(:,:,:) !<derivativeNumbersInLocalFace(0:derivativeIdx,localFaceNodeIdx,localFaceIdx). The element derivative numbers for the derivativeIdx'th face derivative's of the localFaceNodeIdx'th face node in the localFaceIdx'th local face for the basis. The number of derivatives at the localFaceNodeIdx'th face node in the localFaceIdx'th local face is given by derivativeNumbersInLocalFace(0,localFaceNodeIdx,localFaceIdx).
    INTEGER(INTG), ALLOCATABLE :: elementParametersInLocalFace(:,:) !<elementParametersInLocalFace(0:faceParameterIdx,localFaceIdx). The local element parameter for the faceParameterIdx'th face parameter in the localFaceIdx'th local face for the basis. The number of face parameters in the localFaceIdx'th local face is given by elementParametersInLocalFace(0,localFaceIdx).
    !Sub-basis information
    TYPE(BasisPtrType), ALLOCATABLE :: lineBases(:) !<lineBases(xiIdx). lineBases(xiIdx)%ptr is the pointer to the basis for the xiIdx'th xi direction for the basis. Only allocated if the number of xi directions > 1.
    TYPE(BasisPtrType), ALLOCATABLE :: faceBases(:) !<faceBases(xicIdx). faceBases(xicIdx)%ptr pointer to the basis for the xic'th direction of the face normal for the basis. Only allocated if the number of xi directions > 2.
    INTEGER(INTG) :: numberOfSubBases !<The number of sub-bases (lines, faces) for the basis.
    TYPE(BasisPtrType), ALLOCATABLE :: subBases(:) !<subBases(subBasisIdx). The pointer to the subBasisIdx'th sub-basis for the basis.
    TYPE(BasisType), POINTER :: parentBasis !<The pointer to the parent basis for the basis. NOTE: that if the basis is not a sub-basis of another basis this pointer will be NULL. 
  END TYPE BasisType

  !>Contains information on the defined basis functions
  TYPE BasisFunctionsType
    TYPE(ContextType), POINTER :: context !<A pointer to the context for the basis functions
    INTEGER(INTG) :: numberOfBasisFunctions !<The number of basis functions defined
    TYPE(BasisPtrType), ALLOCATABLE :: bases(:) !<The array of pointers to the defined basis functions
  END TYPE BasisFunctionsType

  PUBLIC BasisPtrType,BasisType,BasisFunctionsType
  
  !
  !================================================================================================================================
  !
  ! Coordinate system types
  !

  !>Contains information on a coordinate system. \todo Have a list of coordinate systems and have a pointer in the CoordinateSystemType back to the regions that use them.
  TYPE CoordinateSystemType
    INTEGER(INTG) :: userNumber !<The user defined identifier for the coordinate. The user number must be unique.
    TYPE(CoordinateSystemsType), POINTER :: coordinateSystems !<A pointer back to the coordinate systems for the coordinate system
    LOGICAL :: coordinateSystemFinished !<Is .TRUE. if the coordinate system has finished being created, .FALSE. if not.
    INTEGER(INTG) :: type !<The type of coordinate system. Old CMISS name ITYP10(nr). \see CoordinateRoutines_CoordinateSystemTypes
    INTEGER(INTG) :: radialInterpolationType !<The type of radial interpolation type for non-rectangular cartesian systems. Old CMISS name JTYP10(nr). \see CoordinateRoutines_RadialInterpolations
    INTEGER(INTG) :: numberOfDimensions !<The number of dimensions for the coordinate system. Old CMISS name NJT.
    REAL(DP) :: focus !<The focus of the coordinate system for a prolate-spheriodal coordinate system.
    REAL(DP) :: origin(3) !<origin(coordinateIdx). The coordinateIdx'th component of the origin of the coordinate system wrt the global coordinate system. NOTE: maybe this should be wrt to the parent regions coordinate system - this would then go into the REGION type.
    REAL(DP) :: orientation(3,3) !<ORIENTATION(coordinateIdx1,coordinateIdx2). The orientation matrix for the orientation of the coordinate system wrt to the global coordinate system. NOTE: maybe this should be wrt to the parent regions coordinate system - this would then go into the RegionType.
  END TYPE CoordinateSystemType

  !>A buffer type to allow for an array of pointers to a CoordinateSystemType.
  TYPE CoordinateSystemPtrType
    TYPE(CoordinateSystemType), POINTER :: ptr !<A pointer to the coordinate system
  END TYPE CoordinateSystemPtrType

  !>Contains information on the list of coordinate systems
  TYPE CoordinateSystemsType
    TYPE(ContextType), POINTER :: context !<A pointer to the context for the list of coordinate systems
    INTEGER(INTG) :: numberOfCoordinateSystems !<The number of coordinate systems defined
    TYPE(CoordinateSystemPtrType), POINTER :: coordinateSystems(:) !<coordinateSystems(coordinateSystemIdx). The coordinateSystemIdx'th coordinate system.
  END TYPE CoordinateSystemsType

  PUBLIC CoordinateSystemType,CoordinateSystemPtrType,CoordinateSystemsType
  
  !
  !================================================================================================================================
  !
  ! Data projection types
  !
  
  !>Contains information about a data projection result.
  TYPE DataProjectionResultType
    INTEGER(INTG) :: userDataPointNumber !<The user number of the data point to which the projection result corresponds to.
    INTEGER(INTG) :: globalDataPointNumber !<The global number of the data point to which the projection results cooresponds to.
    REAL(DP) :: distance !<The distances between the data point and the projection. Assigned only if dataPointsProjected is .TRUE.
    INTEGER(INTG) :: elementNumber !<The local element number of the mesh the data point projects onto. Assigned only if dataPointsProjected is .TRUE.
    INTEGER(INTG) :: globalElementNumber !<The global element number of the mesh the data point projects onto. Temp whilst we still have global data points. Assigned only if dataPointsProjected is .TRUE.
    INTEGER(INTG) :: elementLineFaceNumber !<The element line/face of the mesh the data point projects onto. Assigned only if dataPointsProjected is .TRUE. and DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE or DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE is chosen    
    INTEGER(INTG) :: exitTag !<The exit tag of the data projection. Assigned only if dataPointsProjected is .TRUE. \See DataProtectionRoutines,DataProjectionRoutines_DataProjectionTypes 
    REAL(DP), ALLOCATABLE :: xi(:) !<The xi coordinate of the projection. Assigned only if dataPointsProjected is .TRUE.
    REAL(DP), ALLOCATABLE :: elementXi(:) !<The xi coordinates in the element of the projection i.e., the full xi rather than the projected face or line xi. Assigned only if dataPointsProjected is .TRUE.
    REAL(DP), ALLOCATABLE :: projectionVector(:) !<The projection vector from data point to the projected point. 
  END TYPE DataProjectionResultType

  !>Contains information on projection candidates
  TYPE DataProjectionCandidateType
    INTEGER(INTG), ALLOCATABLE :: candidateElementNumbers(:) !<candidateElementNumbers(candidateElementIdx). The candidate element local numbers
    INTEGER(INTG), ALLOCATABLE :: localFaceLineNumbers(:) !<localFaceLineNumbers(candidateElementIdx). The user specified corresponding element face/line numbers for the candidate elements
  END TYPE DataProjectionCandidateType

  !>Contains information on a data point projection
  TYPE DataProjectionType
    INTEGER(INTG) :: globalNumber !<The global number of data projection. 
    INTEGER(INTG) :: userNumber !<The user defined number of data projection. 
    TYPE(VARYING_STRING) :: label !<A string label for the data projection.
    LOGICAL :: dataProjectionFinished !<Is .TRUE. if the data projection has finished being created, .FALSE. if not.
    LOGICAL :: dataProjectionProjected !<Is .TRUE. if the data projection have been projected, .FALSE. if not.
    TYPE(DataPointsType), POINTER :: dataPoints !<The pointer to the data points for this data projection.
    TYPE(DataProjectionsType), POINTER :: dataProjections !<A pointer back to the data projections
    TYPE(FieldType), POINTER :: projectionField !<The pointer to the geometric/dependent field for this data projection.
    TYPE(FieldVariableType), POINTER :: projectionVariable !<A point to the geometric/dependent field variable for this data projection
    INTEGER(INTG) :: projectionVariableType !<The variable type of the geometric/dependent field for this data projection.
    INTEGER(INTG) :: projectionSetType !<The parameter set type of the geometric/dependent field for this data projection. 
    TYPE(DecompositionType), POINTER :: decomposition !<The pointer to the decomposition where data points are projected
    INTEGER(INTG) :: numberOfCoordinates !<The number of coordinates of this data projection.
    INTEGER(INTG) :: numberOfElementXi !<The number of xi of the mesh, ie. the mesh dimension
    INTEGER(INTG) :: numberOfXi !<The number of projection xi.
    INTEGER(INTG) :: projectionType !<type of projection to perform. \See DataProjectionRoutines_DataProjectionTypes
    REAL(DP) :: maximumIterationUpdate !<The maximum xi update allowed at each newton iteration, analogous to maximum trust region size in the trust region model approach.
    INTEGER(INTG) :: maximumNumberOfIterations !<The maximum number of iterations
    INTEGER(INTG) :: numberOfClosestElements !<The number of closest elements to perform full projection on. The algorithm first find the distance of the data point to each elements base on starting xi, full projection is only performed on the first few elements sorted by the distance
    REAL(DP) :: absoluteTolerance !<The absolute tolerance of the iteration update
    REAL(DP) :: relativeTolerance !<The relative tolerance of the iteration update
    REAL(DP), ALLOCATABLE :: startingXi(:) !<The starting value of the element xi
    INTEGER(INTG) :: maxNumberOfCandidates !<The maximum number of projection candidate elements.
    TYPE(DataProjectionCandidateType), ALLOCATABLE :: dataProjectionCandidates(:) !<projectionCandidates(dataIdx). The projection candidates for the dataIdx'th data point. The 0'th index contains the default projection candidates which can then be overridden for specific data points.
    TYPE(DataProjectionResultType), ALLOCATABLE :: dataProjectionResults(:) !<dataProjectionResults(dataIdx). The data projection results for the dataIdx'th data point.
    REAL(DP) :: rmsError !<The RMS error for the data projection.
    REAL(DP) :: maximumError !<The maximum error for the data projection.
    INTEGER(INTG) :: maximumErrorDataPoint !<The global data point number where the maximum error occurs.
    REAL(DP) :: minimumError !<The minimum error for the data projection.
    INTEGER(INTG) :: minimumErrorDataPoint !<The global data point number where the minimum error occurs.
  END TYPE DataProjectionType

  !>A buffer type to allow for an array of pointers to a DataProjectionType.
  TYPE DataProjectionPtrType
    TYPE(DataProjectionType), POINTER :: ptr !<The pointer to the data projection.
  END TYPE DataProjectionPtrType

  !>Contains information on the data point projectiosn defined on data points
  TYPE DataProjectionsType
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer back to the data points
    INTEGER(INTG) :: numberOfDataProjections !<The number of data projections defined.
    TYPE(DataProjectionPtrType), ALLOCATABLE :: dataProjections(:) !<dataProjections(projectionIdx). A pointer to the projection_idx'th data projection for the data points.
  END TYPE DataProjectionsType

  PUBLIC DataProjectionResultType,DataProjectionCandidateType,DataProjectionType,DataProjectionPtrType,DataProjectionsType

  !
  !================================================================================================================================
  !
  ! Data point types
  !

  !>Contains information about a data point.
  TYPE DataPointType
    INTEGER(INTG) :: globalNumber !<The global number of data point. 
    INTEGER(INTG) :: userNumber !<The user defined number of data point. 
    TYPE(VARYING_STRING) :: label !<A string label for the data point.
    REAL(DP), ALLOCATABLE :: position(:) !position(coordinateIdx). Values of the data point specifying the spatial position in the region, has the size of region dimension the data point belongs to.
    REAL(DP), ALLOCATABLE :: weights(:) !<weights(coordinateIdx). Weights of the data point, has the size of region dimension the data point belongs to.
  END TYPE DataPointType

  !>Contains information on the data points defined on a region. \see OpenCMISS::Iron::cmfe_DataPointsType
  TYPE DataPointsType
    INTEGER(INTG) :: globalNumber !<The global number of data points
    INTEGER(INTG) :: userNumber !<The user number of the data points
    TYPE(VARYING_STRING) :: label !<A string label for the data points.
    TYPE(DataPointSetsType), POINTER :: dataPointSets !<The pointer to the data point sets containing these data points.
    TYPE(RegionType), POINTER :: region !<A pointer to the region containing the data points. If the data points are in an interface rather than a region then this pointer will be NULL and the interface pointer should be used.
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface containing the data points. If the data points are in a region rather than an interface then this pointer will be NULL and the interface pointer should be used.
    LOGICAL :: dataPointsFinished !<Is .TRUE. if the data points have finished being created, .FALSE. if not.
    INTEGER(INTG) :: numberOfDimensions !<The number of dimensions for the data points.
    INTEGER(INTG) :: numberOfDataPoints !<The number of data points defined on the region/interface.
    TYPE(DataPointType), ALLOCATABLE :: dataPoints(:) !<dataPoints(dataPointIdx). The data point information for the dataPointIdx'th global data point.
    TYPE(TreeType), POINTER :: dataPointsTree !<The tree for user to global data point mapping.
    TYPE(DataProjectionsType), POINTER :: dataProjections !<dataProjections(projectionIdx). A pointer to the projectionIdx'th data projection for the data points.
  END TYPE DataPointsType

  !>A buffer type to allow for an array of pointers to a DataPointsType.
  TYPE DataPointsPtrType
    TYPE(DataPointsType), POINTER :: ptr !<The pointer to the data points.
  END TYPE DataPointsPtrType

  !>Contains information on the data point sets defined on a region.
  TYPE DataPointSetsType
    TYPE(RegionType), POINTER :: region !<A pointer to the region containg the data points. If the data points are in an interface rather than a region then this pointer will be NULL and the interface pointer should be used.
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface containg the data points. If the data points are in a region rather than an interface then this pointer will be NULL and the region pointer should be used.
    INTEGER(INTG) :: numberOfDataPointSets !<The number of data point sets defined on the region.
    TYPE(DataPointsPtrType), ALLOCATABLE :: dataPointSets(:) !<dataPointSets(setIdx). The array of pointers to the data points.
    TYPE(TreeType), POINTER :: dataPointSetsTree !<The tree for user to global data point sets mapping.
  END TYPE DataPointSetsType

  PUBLIC DataPointType,DataPointsType,DataPointsPtrType,DataPointSetsType

  !
  !================================================================================================================================
  !
  ! Node types
  !

  !>Contains information about a node.
  TYPE NodeType
    INTEGER(INTG) :: globalNumber !<The global number of the node.
    INTEGER(INTG) :: userNumber !<The user defined number of the node.
    TYPE(VARYING_STRING) :: label !<A string label for the node
  END TYPE NodeType

  !>Contains information on the nodes defined on a region. \see OpenCMISS::Iron::cmfe_NodesType
  TYPE NodesType
    TYPE(RegionType), POINTER :: region !<A pointer to the region containing the nodes. If the nodes are in an interface rather than a region then this pointer will be NULL and the interface pointer should be used.
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface containing the nodes. If the nodes are in a region rather than an interface then this pointer will be NULL and the region pointer should be used.
    LOGICAL :: nodesFinished !<Is .TRUE. if the nodes have finished being created, .FALSE. if not.
    INTEGER(INTG) :: numberOfNodes!<The number of nodes defined on the region.
    TYPE(NodeType), ALLOCATABLE :: nodes(:) !<nodes(nodesIdx). The nodal information for the nodesIdx'th global node.
    TYPE(TreeType), POINTER :: nodesTree !<The tree for user to global node mapping
  END TYPE NodesType

  PUBLIC NodeType,NodesType

  !
  !================================================================================================================================
  !
  ! Mesh types
  !

  !>Contains information on the dofs for a mesh.
  TYPE MeshDofsType
    TYPE(MeshTopologyType), POINTER :: meshTopology !<The pointer to the mesh component topology for the dofs information.
    INTEGER(INTG) :: numberOfDofs !<The number of dofs in the mesh.
  END TYPE MeshDofsType

  !>Contains information on the mesh adjacent elements for a xi coordinate
  TYPE MeshAdjacentElementType
    INTEGER(INTG) :: numberOfAdjacentElements !<The number of adjacent elements for the xi coordinate
    INTEGER(INTG), ALLOCATABLE :: adjacentElements(:) !<The global element numbers of the elements adjacent to this element for the xi coordinate
  END TYPE MeshAdjacentElementType
  
  !>Contains the information for an element in a mesh.
  TYPE MeshElementType
    INTEGER(INTG) :: globalNumber !<The global element number in the mesh.
    INTEGER(INTG) :: userNumber !<The corresponding user number for the element.
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis function for the element.
    INTEGER(INTG), ALLOCATABLE :: meshElementNodes(:) !<meshElementNodes(localNodeIdx). The mesh node number in the mesh of the localNodeIdx'th local node in the element. Old CMISS name NPNE(nn,nbf,ne).
    INTEGER(INTG), ALLOCATABLE :: globalElementNodes(:) !<globalElementNodes(localNodeIdx). The global node number in the mesh of the localNodeIdx'th local node in the element. Old CMISS name NPNE(nn,nbf,ne).
    INTEGER(INTG), ALLOCATABLE :: userElementNodeVersions(:,:) !< userElementNodeVersions(derivativeIdx,localNodeIdx).  The version number for the localNodeIdx'th user node's derivativeIdx'th derivative. Size of array dependent on the maximum number of derivatives for the basis of the specified element.
    INTEGER(INTG), ALLOCATABLE :: userElementNodes(:) !<userElementNodes(localNodeIdx). The user node number in the mesh of the localNodeIdx'th local node in the element. Old CMISS name NPNE(nn,nbf,ne).
    TYPE(MeshAdjacentElementType), ALLOCATABLE :: adjacentElements(:) !<adjacentElements(-nic:nic). The adjacent elements information in the nic'th xi coordinate direction. Note that -nic gives the adjacent element before the element in the nic'th direction and +nic gives the adjacent element after the element in the nic'th direction. The ni=0 index will give the information on the current element. Old CMISS name NXI(-ni:ni,0:nei,ne).
    !INTEGER(INTG), ALLOCATABLE :: numberOfAdjacentElements(:) !<numberOfAdjacentElements(-ni:ni). The number of elements adjacent to this element in the ni'th xi direction. Note that -ni gives the adjacent element before the element in the ni'th direction and +ni gives the adjacent element after the element in the ni'th direction. The ni=0 index should be 1 for the current element. Old CMISS name NXI(-ni:ni,0:nei,ne).
    !INTEGER(INTG), ALLOCATABLE :: adjacentElements(:,:) !<adjacentElements(nei,-ni:ni). The local element numbers of the elements adjacent to this element in the ni'th xi direction. Note that -ni gives the adjacent elements before the element in the ni'th direction and +ni gives the adjacent elements after the element in the ni'th direction. The ni=0 index should give the current element number. Old CMISS name NXI(-ni:ni,0:nei,ne)
    LOGICAL :: boundaryElement !<Is .TRUE. if the mesh element is on the boundary of a mesh, .FALSE. if not.
  END TYPE MeshElementType

  !>Contains the information for the elements of a mesh.
  TYPE MeshElementsType
    TYPE(MeshTopologyType), POINTER :: meshTopology !<The pointer to the mesh component topology for the elements information.
    INTEGER(INTG) :: numberOfElements !< The number of elements in the mesh.
    LOGICAL :: elementsFinished !<Is .TRUE. if the mesh elements have finished being created, .FALSE. if not.
    TYPE(MeshElementType), ALLOCATABLE :: elements(:) !<elements(elementIdx). The array of information for the elements of this mesh. elements(elementIdx) contains the information for the elementIdx'th global element of the mesh. 
    TYPE(TreeType), POINTER :: elementsTree !<A tree mapping the mesh global element number to the mesh user element number.
  END TYPE MeshElementsType

  !>Contains the information for a node derivative of a mesh.
  TYPE MeshNodeDerivativeType
    INTEGER(INTG) :: numberOfVersions !The number of global versions at the node for the mesh.
    INTEGER(INTG), ALLOCATABLE :: userVersionNumbers(:) !userVersionNumbers(versionIdx). The user version numbers for the versionIdx'th version for the node.
    INTEGER(INTG), ALLOCATABLE :: dofIndex(:) !The global dof version index (nv) in the domain of the nk'th global derivative for the node.
    INTEGER(INTG) :: globalDerivativeIndex !The global derivative index of the nk'th global derivative for the node.
    INTEGER(INTG) :: partialDerivativeIndex !The partial derivative index (nu) of the nk'th global derivative for the node. Old CMISS name NUNK(nk,nj,np).
  END TYPE MeshNodeDerivativeType

  !>Contains the topology information for a global node of a mesh.
  TYPE MeshNodeType
    INTEGER(INTG) :: meshNumber !<The mesh node number in the mesh.
    INTEGER(INTG) :: globalNumber !<The global node number in the mesh.
    INTEGER(INTG) :: userNumber !<The corresponding user number for the node.
    INTEGER(INTG) :: numberOfDerivatives !<The number of global derivatives at the node for the mesh. Old CMISS name NKT(nj,np).
    TYPE(MeshNodeDerivativeType), ALLOCATABLE :: derivatives(:) !<derivatives(derivativeIdx). Contains information on the derivativeIdx'th derivative of the node.
    INTEGER(INTG) :: numberOfSurroundingElements !<The number of elements surrounding the node in the mesh. Old CMISS name NENP(np,0,0:nr).
    INTEGER(INTG), ALLOCATABLE :: surroundingElements(:) !<surroudingElements(localElementIdx). The global element number of the localElementIdx'th element that is surrounding the node. Old CMISS name NENP(np,nep,0:nr). 
    LOGICAL :: boundaryNode !<Is .TRUE. if the mesh node is on the boundary of the mesh, .FALSE. if not.
  END TYPE MeshNodeType

  !>Contains the information for the nodes of a mesh.
  TYPE MeshNodesType
    TYPE(MeshTopologyType), POINTER :: meshTopology !<The pointer to the mesh component topology for the nodes information.
    INTEGER(INTG) :: numberOfnodes !<The number of nodes in the mesh.
    TYPE(MeshNodeType), ALLOCATABLE :: nodes(:) !<nodes(nodeIdx). The pointer to the array of topology information for the nodes of the mesh. node(nodeIdx) contains the topological information for the nodeIdx'th global node of the mesh. 
    TYPE(TreeType), POINTER :: nodesTree !<A tree mapping the mesh global number to the region nodes global number.
  END TYPE MeshNodesType
  
  TYPE MeshElementDataPointType
    INTEGER(INTG) :: userNumber !<User number of the projected data point
    INTEGER(INTG) :: globalNumber !<Global number of the data point, sequence is according to element number sequence  
  END TYPE MeshElementDataPointType
  
  !>Contains information on the projected data points on an element
  TYPE MeshElementDataPointsType
    INTEGER(INTG) :: numberOfProjectedData !<Number of projected data on this element
    INTEGER(INTG) :: elementNumber !<The mesh global element number 
    TYPE(MeshElementDataPointType), ALLOCATABLE :: dataIndices(:) !<dataIndices(elementDatePointIdx). The global and user number of this data point
  END TYPE MeshElementDataPointsType
  
  !>Contains information of the projected data point 
  TYPE MeshDataPointType 
    INTEGER(INTG) :: userNumber !<User number of the projected data point
    INTEGER(INTG) :: globalNumber !<Global number of the data point, sequence is according to element number sequence 
    INTEGER(INTG) :: elementNumber !<The global element number which the data point is projected on 
  END TYPE MeshDataPointType 
  
  !<Contains the information for the data points of a mesh
  TYPE MeshDataPointsType
    TYPE(MeshTopologyType), POINTER :: meshTopology !<The pointer to the mesh component topology for the data information.
    INTEGER(INTG) :: totalNumberOfProjectedData !<Number of projected data in this mesh
    TYPE(MeshDataPointType), ALLOCATABLE :: dataPoints(:) !<dataPoints(dataPointIdx). Information of the projected data points
    TYPE(MeshElementDataPointsType), ALLOCATABLE :: elementDataPoints(:) !<elementDataPoints(elementIdx). Information of the projected data on the elements 
  END TYPE MeshDataPointsType

  !>Contains information on the (global) topology of a mesh.
  TYPE MeshTopologyType
    TYPE(MeshType), POINTER :: mesh !<Pointer to the parent mesh.
    INTEGER(INTG) :: meshComponentNumber !<The mesh component number for this mesh topology.
    TYPE(MeshNodesType), POINTER :: nodes !<Pointer to the nodes within the mesh topology.
    TYPE(MeshElementsType), POINTER :: elements !<Pointer to the elements within the mesh topology.
    TYPE(MeshDofsType), POINTER :: dofs !<Pointer to the dofs within the mesh topology.
    TYPE(MeshDataPointsType), POINTER :: dataPoints !<Pointer to the data points within the mesh topology
  END TYPE MeshTopologyType

  !>A buffer type to allow for an array of pointers to a MeshTopologyType.
  TYPE MeshTopologyPtrType
    TYPE(MeshTopologytype), POINTER :: ptr !<The pointer to the mesh topology.
  END TYPE MeshTopologyPtrType

  !>Embedded mesh types
  TYPE EMBEDDING_XI_TYPE
    INTEGER(INTG) :: NUMBER_OF_NODES !<Number of nodes embedded in this element
    INTEGER(INTG), ALLOCATABLE :: NODE_NUMBERS(:) !<NODE_NUMBERS(node_idx) Node numbers in child mesh for the node_idx'th embedded node in this element
    REAL(DP), ALLOCATABLE :: XI_COORDS(:,:) !<XI_COORDS(:,node_idx) Xi coordinates of the node_idx'th embedded node this element
  END TYPE EMBEDDING_XI_TYPE

  TYPE EMBEDDING_GAUSSPOINT_TYPE
    INTEGER(INTG) :: ELEMENT_NUMBER !<Element number in child mesh
    REAL(DP), ALLOCATABLE :: CHILD_XI_COORD(:) !<Xi coord in this element
    REAL(DP), ALLOCATABLE :: PARENT_XI_COORD(:) !<Xi coordinates in parent element, not really needed but can be useful
  END TYPE EMBEDDING_GAUSSPOINT_TYPE

  TYPE MESH_EMBEDDING_TYPE
    TYPE(MeshType), POINTER :: CHILD_MESH, PARENT_MESH !<The parent and child mesh
    TYPE(EMBEDDING_XI_TYPE), ALLOCATABLE :: CHILD_NODE_XI_POSITION(:)            !<Location of embedded nodes for each element in the parent mesh
    TYPE(EMBEDDING_GAUSSPOINT_TYPE), ALLOCATABLE :: GAUSS_POINT_XI_POSITION(:,:) !<GAUSS_POINT_XI_POSITION(gauss_idx,element_idx) Location of the gauss_idx'th Gauss point of the element_idx'th parent mesh in the child mesh
  END TYPE MESH_EMBEDDING_TYPE

  !>Contains information on a mesh defined on a region. \see OpenCMISS::Iron::cmfe_MeshType
  TYPE MeshType
    INTEGER(INTG) :: userNumber !<The user number of the mesh. The user number must be unique.
    INTEGER(INTG) :: globalNumber !<The corresponding global number for the mesh.
    LOGICAL :: meshFinished !<Is .TRUE. if the mesh has finished being created, .FALSE. if not.
    TYPE(MeshesType), POINTER :: meshes !<A pointer to the meshes for this mesh.
    TYPE(RegionType), POINTER :: region !<A pointer to the region containing the mesh. If the mesh is in an interface rather than a region then this pointer will be NULL and the interface pointer should be used.
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface containing the mesh. If the mesh is in a region rather than an interface then this pointer will be NULL and the interface pointer should be used.
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh generate this mesh.
    INTEGER(INTG) :: numberOfDimensions !<The number of dimensions (Xi directions) for this mesh.
    INTEGER(INTG) :: numberOfComponents !<The number of mesh components in this mesh.
    LOGICAL :: MESH_EMBEDDED !<Is .TRUE. if the mesh is embedded in another mesh, .FALSE. if not.
    TYPE(MeshType), POINTER :: EMBEDDING_MESH !<If this mesh is embedded the pointer to the mesh that this mesh is embedded in. IF the mesh is not embedded the pointer is NULL.
    INTEGER(INTG) :: NUMBER_OF_EMBEDDED_MESHES !<The number of meshes that are embedded in this mesh.
    TYPE(MeshPtrType), POINTER :: EMBEDDED_MESHES(:) !<EMBEDDED_MESHES(mesh_idx). A pointer to the mesh_idx'th mesh that is embedded in this mesh.
    INTEGER(INTG) :: numberOfElements
    TYPE(MeshTopologyPtrType), ALLOCATABLE :: topology(:) !<topology(meshComponentIdx). A pointer to the topology meshComponentIdx'th mesh component.
    TYPE(DecompositionsType), POINTER :: decompositions !<A pointer to the decompositions for this mesh.
    LOGICAL :: surroundingElementsCalculate !<Boolean flag to determine whether surrounding elements should be calculated.
  END TYPE MeshType

  !>A buffer type to allow for an array of pointers to a MeshType.
  TYPE MeshPtrType
    TYPE(MeshType), POINTER :: ptr !<The pointer to the mesh. 
  END TYPE MeshPtrType

  !>Contains information on the meshes defined on a region.
  TYPE MeshesType
    TYPE(RegionType), POINTER :: region !<A pointer to the region containg the meshes. If the meshes are in an interface rather than a region then this pointer will be NULL and the interface pointer should be used.
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface containg the meshes. If the meshes are in a region rather than an interface then this pointer will be NULL and the region pointer should be used.
    INTEGER(INTG) :: numberOfMeshes !<The number of meshes defined on the region.
    TYPE(MeshPtrType), ALLOCATABLE :: meshes(:) !<meshes(meshesIdx). The array of pointers to the meshes.
  END TYPE MeshesType

  PUBLIC MeshDofsType

  PUBLIC MeshAdjacentElementType,MeshElementType,MeshElementsType

  PUBLIC MeshNodeDerivativeType,MeshNodeType,MeshNodesType

  PUBLIC MeshElementDataPointsType,MeshDataPointType,MeshDataPointsType

  PUBLIC MeshTopologytype,MeshTopologyPtrType

  PUBLIC EMBEDDING_XI_TYPE,EMBEDDING_GAUSSPOINT_TYPE,MESH_EMBEDDING_TYPE

  PUBLIC MeshType,MeshPtrType,MeshesType

  !
  !================================================================================================================================
  !
  ! Generated Mesh types
  !

  !>Contains information on a generated regular mesh
  TYPE GeneratedMeshRegularType
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh
    TYPE(BasisPtrType), ALLOCATABLE :: bases(:) !<The pointers to the bases used in the regular mesh. TODO: move up a level and combine
    INTEGER(INTG) :: coordinateDimension !<The number of coordinates for the regular mesh.
    INTEGER(INTG) :: meshDimension !<The dimension/number of Xi directions of the regular mesh.
    REAL(DP), ALLOCATABLE :: origin(:) !<origin(coordinateIdx). The position of the origin (first) corner of the regular mesh
    REAL(DP), ALLOCATABLE :: maximumExtent(:) !<maximumExtent(coordinateIdx). The extent/size in each coordinateIdx'th direction of the regular mesh.
    REAL(DP), ALLOCATABLE :: baseVectors(:,:) !<baseVectors(coordinateIdx,xiIdx). The base vectors indicating the geometrinc direction of the xiIdx mesh coordinate.
    INTEGER(INTG), ALLOCATABLE :: numberOfElementsXi(:) !<numberOfElementsXi(xiIdx). The number of elements in the xiIdx'th Xi direction for the mesh.
  END TYPE GeneratedMeshRegularType

  !>Contains information of a generated cylinder mesh
  !>Allows only a 3D cylinder mesh with xi directions (r,theta,z)
  TYPE GeneratedMeshCylinderType
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh.
    REAL(DP), ALLOCATABLE :: origin(:) !<origin(coordinateIdx). The position of the origin (centre) of lower face of cylinder mesh.
    REAL(DP), ALLOCATABLE :: cylinderExtent(:) !<cylinderExtent(coordinateIdx). The size of inner & outer radii and height of cylinder.
    INTEGER(INTG) :: meshDimension !<The dimension/number of Xi directions of the cylinder mesh.
    INTEGER(INTG), ALLOCATABLE :: numberOfElementsXi(:) !<NUMBER_OF_ELEMENTS(xiIdx). The number of elements in radial, circumferential and axial directions
    TYPE(BasisPtrType), ALLOCATABLE :: bases(:) !<The pointers to the bases used in the cylinder mesh. TODO: move up a level and combine
    LOGICAL :: appendLinearComponent=.FALSE. !<True when two mesh components are needed
 END TYPE GeneratedMeshCylinderType
 
  !>Contains information of a generated ellipsoid mesh
  !>Allows only a 3D ellipsoid mesh
 TYPE GeneratedMeshEllipsoidType
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<A pointer to the generated mesh.
    REAL(DP), ALLOCATABLE :: origin(:) !<origin(coordinateIdx). The position of the origin (centre) of lower face of ellipsoid mesh.
    REAL(DP), ALLOCATABLE :: ellipsoidExtent(:) !<ellipsoidExtent(coordinateIdx). The size of long axis, short axis, wall thickness and cut off angle of ellipsoid.
    INTEGER(INTG) :: meshDimension !<The dimension/number of Xi directions of the ellipsoid mesh.
    INTEGER(INTG), ALLOCATABLE :: numberOfElementsXi(:) !<numberOfElementsXi(xiIdx). The number of elements in circumferential, longitudinal and transmural directions
    TYPE(BasisPtrType), ALLOCATABLE :: bases(:) !<The pointers to the bases used in the ellipsoid mesh. TODO: move up a level and combine
    LOGICAL :: appendLinearComponent=.FALSE. !<True when two mesh components are needed 
END TYPE GeneratedMeshEllipsoidType

  !>Contains information on a generated mesh. \see OpenCMISS::Iron::cmfe_GeneratedMeshType
  TYPE GeneratedMeshType
    INTEGER(INTG) :: userNumber !<The user number of the generated mesh. The user number must be unique.
    INTEGER(INTG) :: globalNumber !<The corresponding global number for the generated mesh.
    TYPE(GeneratedMeshesType), POINTER :: generatedMeshes !<Is .TRUE. if the generated mesh has finished being created, .FALSE. if not.
    LOGICAL :: generatedMeshFinished !<Is .TRUE. if the generated mesh has finished being created, .FALSE. if not.
    TYPE(RegionType), POINTER :: region !<A pointer to the region containing the generated mesh. If the generated mesh is in an interface rather than a region then this pointer will be NULL and the interface pointer should be used.
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface containing the generated mesh. If the generated mesh is in a region rather than an interface then this pointer will be NULL and the interface pointer should be used.
    INTEGER(INTG) :: generatedType !<The type of the generated mesh. \see GeneratedMeshRoutines_GeneratedMeshTypes,GeneratedMeshRoutines
    TYPE(GeneratedMeshRegularType), POINTER :: regularMesh !<A pointer to the regular generated mesh information if the generated mesh is a regular mesh, NULL if not.
    TYPE(GeneratedMeshCylinderType), POINTER :: cylinderMesh !<A pointer to the cylinder generate mesh information if the generated mesh is a cylinder mesh, NULL if not.
    TYPE(GeneratedMeshEllipsoidType), POINTER :: ellipsoidMesh !<A pointer to the ellipsoid generate mesh information if the generated mesh is a ellipsoid mesh, NULL if not.
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh that has been generated.
  END TYPE GeneratedMeshType
  
  !>A buffer type to allow for an array of pointers to a GeneratedMeshType.
  TYPE GeneratedMeshPtrType
    TYPE(GeneratedMeshType), POINTER :: ptr !<The pointer to the generated mesh.
  END TYPE GeneratedMeshPtrType
       
  !>Contains information on the generated meshes defined on a region.
  TYPE GeneratedMeshesType
    TYPE(RegionType), POINTER :: region !<A pointer to the region containing the generated meshes. If the generated meshes are in an interface rather than a region then this pointer will be NULL and the interface pointer should be used.
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface containing the generated meshes. If the generated meshes are in a region rather than an interface then this pointer will be NULL and the interface pointer should be used.
    INTEGER(INTG) :: numberOfGeneratedMeshes !<The number of generated meshes defined.
    TYPE(GeneratedMeshPtrType), ALLOCATABLE :: generatedMeshes(:) !<The array of pointers to the generated meshes.
  END TYPE GeneratedMeshesType

  PUBLIC GeneratedMeshRegularType,GeneratedMeshCylinderType,GeneratedMeshEllipsoidType

  PUBLIC GeneratedMeshType,GeneratedMeshPtrType,GeneratedMeshesType
  
  !
  !================================================================================================================================
  !
  ! Domain types
  !
  
  !>Contains information on the degrees-of-freedom (dofs) for a domain.
  TYPE DomainDOFsType
    TYPE(DomainTopologyType), POINTER :: domainTopology !<A pointer back to the domain topology.
    INTEGER(INTG) :: numberOfDofs !<The number of degrees-of-freedom (excluding ghost dofs) in the domain.
    INTEGER(INTG) :: totalNumberOfDofs !<The total number of degrees-of-freedom (including ghost dofs) in the domain.
    INTEGER(INTG) :: numberOfGlobalDofs !<The number of global degrees-of-freedom in the domains.
    INTEGER(INTG), ALLOCATABLE :: dofIndex(:,:) !<dofIndex(i,ny). The index for the ny'th degree-of-freedom. When i=1 dofIndex will give the global derivative version number (version_idx) associated with the dof. When i=2 dofIndex will give the global derivative number (derivative_idx) associated with the dof. When i=3 dofIndex will give the local node number (node_idx) associated with the dof.
  END TYPE DomainDOFsType

  !>Contains the information for a line in a domain.
  TYPE DomainLineType
    INTEGER(INTG) :: number !<The line number in the domain.
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis function for the line.
    INTEGER(INTG), ALLOCATABLE :: nodesInLine(:) !<nodesInLine(localNodeIdx). The local node number in the domain of the localNodeIdx'th local node in the line. Old CMISS name NPL(2..5,nj,nl).
    INTEGER(INTG), ALLOCATABLE :: derivativesInLine(:,:,:) !<derivativesInLine(i,localDerivativeIdx,localNodeIdx). When i=1 derivativesInLine will give the global derivative number of the localDerivativeIdx'th local derivative for the localNodeIdx'th local node in the line, When i=2 derivativesInLine will give the global derivative version number of the localDerivativeIdx'th local derivative for the localNodeIdx'th local node in the line. Old CMISS name NPL(4..5,nj,nl).
    LOGICAL :: boundaryLine !<Is .TRUE. if the line is on the boundary of the mesh for the domain, .FALSE. if not.
    INTEGER(INTG) :: elementNumber !<The element number of the element on which the line is on
  END TYPE DomainLineType

  !>A buffer type to allow for an array of pointers to a DomainLineType
  TYPE DomainLinePtrType
    TYPE(DomainLineType), POINTER :: ptr !<A pointer to the domain line.
  END TYPE DomainLinePtrType

  !>Contains the topology information for the lines of a domain.
  TYPE DomainLinesType
    TYPE(DomainTopologyType), POINTER :: domainTopology !<The pointer back to the domain topology.
    INTEGER(INTG) :: numberOfLines !<The number of lines in this domain topology.
    TYPE(DomainLineType), ALLOCATABLE :: lines(:) !<lines(lineIdx). The pointer to the array of topology information for the lines of this domain. lines(lineIdx) contains the topological information for the lineIdx'th local line of the domain.
  END TYPE DomainLinesType

  !>Contains the information for a face in a domain.
  TYPE DomainFaceType
    INTEGER(INTG) :: number !<The face number in the domain.
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis function for the face.
    INTEGER(INTG), ALLOCATABLE :: nodesInFace(:) !<nodesInFace(localFaceNodeIdx). The local node number in the domain of the localFaceNodeIdx'th local node in the face. Old CMISS name NPNF(nn,nbf).
    INTEGER(INTG), ALLOCATABLE :: derivativesInFace(:,:,:) !<derivativesInFace(i,localDerivativeIdx,localNodeIdx). When i=1 derivativesInFace will give the global derivative number of the local derivativeIdx'th local derivative for the localNodeIdx'th local node in the face, When i=2 derivativesInFace will give the global derivative version number of the local derivativeIdx'th local derivative for the localNodeIdx'th local node in the face.
    LOGICAL :: boundaryFace !<Is .TRUE. if the face is on the boundary of the mesh for the domain, .FALSE. if not.
    INTEGER(INTG) :: elementNumber !<The element number of the element on which the face is on
  END TYPE DomainFaceType

  !>A buffer type to allow for an array of pointers to a DomainFaceType.
  TYPE DomainFacePtrType
    TYPE(DomainFaceType), POINTER :: ptr !<The pointer to the domain face.
  END TYPE DomainFacePtrType

  !>Contains the topology information for the faces of a domain.
  TYPE DomainFacesType
    TYPE(DomainTopologyType), POINTER :: domainTopology !<The pointer back to the domain topology.
    INTEGER(INTG) :: numberOfFaces !<The number of faces in this domain topology.
    TYPE(DomainFaceType), ALLOCATABLE :: faces(:) !<faces(faceIdx). The pointer to the array of topology information for the faces of this domain. facecs(faceIdx) contains the topological information for the faceIdx'th local face of the domain.
  END TYPE DomainFacesType

  !>Contains the information for an element in a domain.
  TYPE DomainElementType
    INTEGER(INTG) :: number !<The local element number in the domain.
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis function for the element.
    INTEGER(INTG), ALLOCATABLE :: elementNodes(:) !<elementNodes(localElementNodeIdx). The local node number in the domain of the localElementNodeIdx'th local node in the element. Old CMISS name NPNE(nn,nbf,ne).
    INTEGER(INTG), ALLOCATABLE :: elementDerivatives(:,:) !<elementDerivatives(localDerivativeIdx,localElementNodeIdx).  Will give the global derivative number of the local_derivative_idx'th local derivative for the local_element_node_idx'th local node in the element. Old CMISS name NKJE(nk,nn,nj,ne).
    INTEGER(INTG), ALLOCATABLE :: elementVersions(:,:) !<elementVersions(localDerivativeIdx,localElementNodeIdx). will give the version number of the localDerivativeIdx'th local derivative for the localElementNodeIdx'th local node in the element. Old CMISS name NVJE(nn,nbf,nj,ne).
  END TYPE DomainElementType
  
  !>Contains the topology information for the elements of a domain.
  TYPE DomainElementsType
    TYPE(DomainTopologyType), POINTER :: domainTopology !<The pointer back to the domain topology.
    INTEGER(INTG) :: numberOfElements !<The number of elements (excluding ghost elements) in this domain topology.
    INTEGER(INTG) :: totalNumberOfElements !<The total number of elements (including ghost elements) in this domain topology.
    INTEGER(INTG) :: numberOfGlobalElements !<The number of global elements in this domain topology.
    TYPE(DomainElementType), ALLOCATABLE :: elements(:) !<elements(elementIdx). The pointer to the array of topology information for the elements of this domain. elements(elementIdx) contains the topological information for the elementIdx'th local elements of the domain.
    INTEGER(INTG) :: maximumNumberOfElementParameters !<The maximum number of element parameters (ns) for all the elements in the domain.
  END TYPE DomainElementsType

  !>Contains the topology information for a local node derivative of a domain.
  TYPE DomainNodeDerivativeType
    INTEGER(INTG) :: numberOfVersions !<The number of global versions at the node for the mesh.
    INTEGER(INTG), ALLOCATABLE :: userVersionNumbers(:) !<userVersionNumbers(derivativeIdx). The user version index of the derivativeIdx'th global derivative for the node.
    INTEGER(INTG), ALLOCATABLE :: dofIndex(:) !<dofIndex(derivativeIdx). The local dof derivative version index in the domain of the derivativeIdx'th global derivative for the node.
    INTEGER(INTG) :: globalDerivativeIndex !<The global derivative index of the nk'th global derivative for the node.
    INTEGER(INTG) :: partialDerivativeIndex !<The partial derivative index (nu) of the nk'th global derivative for the node. Old CMISS name NUNK(nk,nj,np).
  END TYPE DomainNodeDerivativeType

  !>Contains the topology information for a local node of a domain.
  TYPE DomainNodeType
    INTEGER(INTG) :: localNumber !<The local node number in the domain.
    INTEGER(INTG) :: meshNumber !<The corresponding global node number in the mesh of the local node number in the domain i.e., the mesh node number.
    INTEGER(INTG) :: globalNumber !<The corresponding global number for the node i.e., the node number in the list of nodes for the region.
    INTEGER(INTG) :: userNumber !<The corresponding user number for the node.
    INTEGER(INTG) :: numberOfDerivatives !<The number of global derivatives at the node for the domain. Old CMISS name NKT(nj,np)
    TYPE(DomainNodeDerivativeType), ALLOCATABLE :: derivatives(:) !<derivatives(derivativeIdx). Derivative information for derivativeIdx'th derivative at the node
    INTEGER(INTG) :: numberOfSurroundingElements !<The number of elements surrounding the node in the domain. Old CMISS name NENP(np,0,0:nr).
    INTEGER(INTG), ALLOCATABLE :: surroundingElements(:) !<surroundingElements(surroundingElementIdx). The local element number of the surroundingElementIdx'th element that is surrounding the node. Old CMISS name NENP(np,nep,0:nr).
    INTEGER(INTG) :: numberOfNodeLines !<The number of lines surrounding the node in the domain.
    INTEGER(INTG), ALLOCATABLE :: nodeLines(:) !<nodeLines(lineIdx). The local line number of the lineIdx'th line that is surrounding the node.
    INTEGER(INTG) :: numberOfNodeFaces !<The number of faces surrounding the node in the domain.
    INTEGER(INTG), ALLOCATABLE :: nodeFaces(:) !<nodeFaces(faceIdx). The local face number of the faceIdx'th face that is surrounding the node.
    LOGICAL :: boundaryNode !<Is .TRUE. if the node is on the boundary of the mesh for the domain, .FALSE. if not.
  END TYPE DomainNodeType

  !>Contains the topology information for the nodes of a domain
  TYPE DomainNodesType
    TYPE(DomainTopologyType), POINTER :: domainTopology !<The pointer back to the domain topology
    INTEGER(INTG) :: numberOfNodes !<The number of nodes (excluding ghost nodes) in this domain topology.
    INTEGER(INTG) :: totalNumberOfNodes !<The total number of nodes (including ghost nodes) in this domain topology.
    INTEGER(INTG) :: numberOfGlobalNodes !<The number of global nodes in this domain topology.
    INTEGER(INTG) :: maximumNumberOfDerivatives !<The maximum number of derivatives over the nodes in this domain topology.
    TYPE(DomainNodeType), ALLOCATABLE :: nodes(:) !<nodes(nodeIdx). The pointer to the array of topology information for the nodes of this domain. nodes(nodeIdx) contains the topological information for the nodeIdx'th local node of the domain.
    TYPE(TreeType), POINTER :: nodesTree !<A tree mapping the domain local number to the region nodes user number.
  END TYPE DomainNodesType

  !>Contains the topology information for a domain
  TYPE DomainTopologyType
    TYPE(DomainType), POINTER :: domain !<The pointer to the domain for this topology information.
    TYPE(DomainNodesType), POINTER :: nodes !<The pointer to the topology information for the nodes of this domain.
    TYPE(DomainDOFsType), POINTER :: dofs !<The pointer to the topology information for the dofs of this domain. /todo think about getting rid of domain dofs
    TYPE(DomainElementsType), POINTER :: elements !<The pointer to the topology information for the elements of this domain.
    TYPE(DomainFacesType), POINTER :: faces !<The pointer to the topology information for the faces of this domain.
    TYPE(DomainLinesType), POINTER :: lines !<The pointer to the topology information for the lines of this domain.
  END TYPE DomainTopologyType

  PUBLIC DomainDOFsType

  PUBLIC DomainLineType,DomainLinePtrType,DomainLinesType

  PUBLIC DomainFaceType,DomainFacePtrType,DomainFacesType

  PUBLIC DomainElementType,DomainElementsType

  PUBLIC DomainNodeDerivativeType,DomainNodeType,DomainNodesType

  PUBLIC DomainTopologyType  

  !
  !================================================================================================================================
  !
  ! Distributed matrix vector types
  !

  !>Contains the information for an adjacent domain for transfering the ghost data of a distributed vector to/from the
  !>current domain.
  TYPE DistributedVectorTransferType
    TYPE(DistributedVectorCMISSType), POINTER :: cmissVector !<The pointer to the CMISS distributed vector object for this transfer information.
    INTEGER(INTG) :: dataType !<The data type of the distributed vector. This is "inherited" from the distributed vector.
    INTEGER(INTG) :: sendBufferSize !<The size of the buffer to send distributed vector data from the current domain to the adjacent domain.
    INTEGER(INTG) :: receiveBufferSize !<The size of the buffer to receive distributed vector data from the adjacent domain to the current domain.
    INTEGER(INTG) :: sendTagNumber !<The MPI tag number for the data sending from the current domain to the adjacent domain. It is calculated as an offset from the base tag number of the distribued vector.
    INTEGER(INTG) :: receiveTagNumber !<The MPI tag number for the data receiving from the adjacent domain to the current domain. It is calculated as an offset from the base tag number of the distribued vector.
    INTEGER(INTG) :: mpiSendRequest !<The MPI request pointer for sending data from the current domain to the adjacent domain.
    INTEGER(INTG) :: mpiReceiveRequest !<The MPI request pointer for sending data from the adjacent domain to the current domain.
    INTEGER(INTG), ALLOCATABLE :: sendBufferIntg(:) !<The integer buffer for sending the distributed integer vector data from the current domain to the adjacent domain.
    REAL(DP), ALLOCATABLE :: sendBufferDP(:) !<The double precision real buffer for sending the distributed real vector data from the current domain to the adjacent domain.
    REAL(SP), ALLOCATABLE :: sendBufferSP(:) !<The single precision real buffer for sending the distributed real vector data from the current domain to the adjacent domain.
    LOGICAL, ALLOCATABLE :: sendBufferL(:) !<The logical buffer for sending the distributed logical vector data from the current domain to the adjacent domain.
    INTEGER(INTG), ALLOCATABLE :: receiveBufferIntg(:) !<The integer buffer for receiving the distributed integer vector data from the adjacent domain to the current domain.
    REAL(DP), ALLOCATABLE :: receiveBufferDP(:) !<The double precision real buffer for receiving the distributed real vector data from the adjacent domain to the current domain.
    REAL(SP), ALLOCATABLE :: receiveBufferSP(:) !<The single precision real buffer for receiving the distributed real vector data from the adjacent domain to the current domain.
    LOGICAL, ALLOCATABLE :: receiveBufferL(:) !<The logical buffer for receiving the distributed logical vector data from the adjacent domain to the current domain.  
  END TYPE DistributedVectorTransferType

  !>Contains information for a CMISS distributed vector
  TYPE DistributedVectorCMISSType
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG) :: baseTagNumber !<The base number for the MPI tag numbers that will be used to communicate the distributed vector data amongst the domains. The base tag number can be thought of as the identification number for the distributed vector object.
    INTEGER(INTG) :: n !<The size of the distributed vector
    INTEGER(INTG) :: dataSize !<The size of the distributed vector that is held locally by the domain.
    INTEGER(INTG), ALLOCATABLE :: dataIntg(:) !<dataIntg(i). The integer data for an integer distributed vector. The i'th component contains the data for the i'th local number of distributed vector data on the domain. 
    REAL(DP), ALLOCATABLE :: dataDP(:) !<dataDP(i). The real data for a double precision real distributed vector. The i'th component contains the data for the i'th local number of distributed vector data on the domain. 
    REAL(SP), ALLOCATABLE :: dataSP(:) !<dataSP(i). The real data for a single precision real distributed vector. The i'th component contains the data for the i'th local number of distributed vector data on the domain. 
    LOGICAL, ALLOCATABLE :: dataL(:) !<dataL(i). The logical data for a logical distributed vector. The i'th component contains the data for the i'th local number of distributed vector data on the domain.  
    TYPE(DistributedVectorTransferType), ALLOCATABLE :: transfers(:) !<transfers(adjacentDomainIdx). The transfer information for the adjacentDomainIdx'th adjacent domain to this domain. 
  END TYPE DistributedVectorCMISSType

  !>Contains information for a PETSc distributed vector
  TYPE DistributedVectorPETScType
    TYPE(DistributedVectorType), POINTER :: distributedVector !<A pointer to the distributed vector
    INTEGER(INTG) :: n !<The number of local components in the vector
    INTEGER(INTG) :: globalN !<The number of global components in the vector
    INTEGER(INTG) :: dataSize  !<The size of the distributed vector that is held locally by the domain.
    INTEGER(INTG), ALLOCATABLE :: globalNumbers(:) !<globalNumbers(i). The PETSc global number corresponding to the i'th local number.
    LOGICAL :: useOverrideVector !<Is .TRUE. if the override vector is used instead of the standard vector
    TYPE(PetscVecType) :: vector !<The PETSc vector
    TYPE(PetscVecType) :: overrideVector !<The PETSc override vector
  END TYPE DistributedVectorPETScType
  
  !>Contains the information for a vector that is distributed across a number of domains.
  TYPE DistributedVectorType
    LOGICAL :: vectorFinished !<!<Is .TRUE. if the distributed vector has finished being created, .FALSE. if not.
    INTEGER(INTG) :: libraryType !<The format of the distributed vector \see DistributedMatrixVector_LibraryTypes,DistributedMatrixVector
    INTEGER(INTG) :: ghostingType !<The ghosting type \see DistributedMatrixVector_GhostingTypes,DistributedMatrixVector
    TYPE(DomainMappingType), POINTER :: domainMapping !<The pointer for the domain mapping that identifies how the vector is distributed amongst the domains.
    INTEGER(INTG) :: dataType !<The type of data for the distributed vector \see DistributedMatrixVector_DataTypes 
    TYPE(DistributedVectorCMISSType), POINTER :: cmiss !<A pointer to the CMISS distributed vector information
    TYPE(DistributedVectorPETScType), POINTER :: petsc !<A pointer to the PETSc distributed vector information
  END TYPE DistributedVectorType

  !>Contains information for a CMISS distributed matrix
  TYPE DistributedMatrixCMISSType
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG) :: baseTagNumber !<The base number for the MPI tag numbers that will be used to communicate the distributed matrix data amongst the domains. The base tag number can be thought of as the identification number for the distributed matrix object.
    INTEGER(INTG) :: globalM !<The number of global rows in the matrix
    TYPE(MatrixType), POINTER :: matrix !<A pointer to the matrix to store the rows corresponding to this domain.
  END TYPE DistributedMatrixCMISSType

  !>Contains information for a PETSc distributed matrix
  TYPE DistributedMatrixPETScType
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix !<A pointer to the distributed matrix
    INTEGER(INTG) :: m !<The number of local rows in the PETSc matrix
    INTEGER(INTG) :: n !<The number of local columns in the PETSc matrix
    INTEGER(INTG) :: globalM !<The number of global rows in the PETSc matrix
    INTEGER(INTG) :: globalN !<The number of global columns in the PETSc matrix
    INTEGER(INTG) :: storageType !<The storage type (sparsity) of the PETSc matrix
    INTEGER(INTG) :: symmetryType !<The symmetry type of the PETSc matrix \see DistributedMatrixVector_SymmetryTypes
    INTEGER(INTG) :: transposeType !<The transpose type of the matrix \see DistributedMatrixVector_TransposeTypes
    INTEGER(INTG) :: numberOfNonZeros !<The number of non-zeros in the PETSc matrix
    INTEGER(INTG) :: dataSize !<The size of the allocated data in the PETSc matrix
    INTEGER(INTG) :: maximumColumnIndicesPerRow!<The maximum number of column indicies for the rows.
    INTEGER(INTG), ALLOCATABLE :: diagonalNumberOfNonZeros(:) !<diagonalNumberOfNonZeros(i). The number of non-zeros in the diagonal part of the the i'th row
    INTEGER(INTG), ALLOCATABLE :: offdiagonalNumberOfNonZeros(:) !<offdiagonalNumberOfNonZeros(i). The number of non-zeros in the off diagonal part of the the i'th row
    INTEGER(INTG), ALLOCATABLE :: rowIndices(:) !<rowIndices(i). The row indices for the matrix.
    INTEGER(INTG), ALLOCATABLE :: columnIndices(:) !<columnIndices(i). The column indices for the matrix.
    TYPE(LINKEDLIST),POINTER :: list(:) !< \todo Comment
    INTEGER(INTG), ALLOCATABLE :: globalRowNumbers(:) !<globalRowNumbers(i). The PETSc global row number corresponding to the i'th local row number.
    REAL(DP), POINTER :: dataDP(:) !<dataDP(i). The real data for the matrix. \todo Is this used???
    LOGICAL :: useOverrideMatrix !<Is .TRUE. if the override matrix is to be used instead of the standard matrix
    TYPE(PetscMatType) :: matrix !<The PETSc matrix
    TYPE(PetscMatType) :: overrideMatrix !<The PETSc override matrix
  END TYPE DistributedMatrixPETScType
  
  !>Contains the information for a matrix that is distributed across a number of domains.
  TYPE DistributedMatrixType
    LOGICAL :: matrixFinished !<Is .TRUE. if the distributed matrix has finished being created, .FALSE. if not.
    INTEGER(INTG) :: libraryType !<The library of the distributed matrix \see DistributedMatrixVector_LibraryTypes,DistributedMatrixVector
    INTEGER(INTG) :: ghostingType !<The ghosting type \see DistributedMatrixVector_GhostingTypes,DistributedMatrixVector
    TYPE(DomainMappingType), POINTER :: rowDomainMapping !<The pointer for the domain mapping that identifies how the matrix rows are distributed amongst the domains.
    TYPE(DomainMappingType), POINTER :: columnDomainMapping !<The pointer for the domain mapping that identifies how the matrix columns are distributed amongst the domains.
    INTEGER(INTG) :: dataType !<The type of data for the distributed matrix \see DistributedMatrixVector_DataTypes
    TYPE(DistributedMatrixCMISSType), POINTER :: cmiss !<A pointer to the CMISS distributed matrix information
    TYPE(DistributedMatrixPETScType), POINTER :: petsc !<A pointer to the PETSc distributed matrix information
  END TYPE DistributedMatrixType

  PUBLIC DistributedVectorTransferType,DistributedVectorCMISSType,DistributedVectorPETScType,DistributedVectorType

  PUBLIC DistributedMatrixCMISSType,DistributedMatrixPETScType,DistributedMatrixType

  !
  !================================================================================================================================
  !
  ! Matrix vector types
  !

  !>Contains information for a vector
  TYPE VectorType
    INTEGER(INTG) :: ID !<The ID of the vector.
    LOGICAL :: vectorFinished !<Is .TRUE. if the vector has finished being created, .FALSE. if not.
    INTEGER(INTG) :: n !<The length of the vector
    INTEGER(INTG) :: dataType !<The data type of the vector \see MatrixVector_DataTypes 
    INTEGER(INTG) :: size !<The size of the data array of the vector
    INTEGER(INTG), ALLOCATABLE :: dataIntg(:) !<dataIntg(i). The integer data for an integer vector. The i'th component contains the data for the i'th component vector data on the domain. 
    REAL(SP), ALLOCATABLE :: dataSP(:) !<dataSP(i). The real data for a single precision real vector. The i'th component contains the data for the i'th component vector data on the domain. 
    REAL(DP), ALLOCATABLE :: dataDP(:) !<dataDP(i). The real data for a double precision real vector. The i'th component contains the data for the i'th component vector data on the domain. 
    LOGICAL, ALLOCATABLE :: dataL(:) !<dataL(i). The real for a logical vector. The i'th component contains the data for the i'th component vector data on the domain. 
  END TYPE VectorType

  !>Contains information for a matrix
  TYPE MatrixType
    INTEGER(INTG) :: ID !<The ID of the matrix
    LOGICAL :: matrixFinished !<Is .TRUE. if the matrix has finished being created, .FALSE. if not.
    INTEGER(INTG) :: m !<The number of rows in the matrix
    INTEGER(INTG) :: n !<The number of columns in the matrix
    INTEGER(INTG) :: maxM !<The maximum number of columns in the matrix storage
    INTEGER(INTG) :: maxN !<The maximum number of rows in the matrix storage
    INTEGER(INTG) :: dataType !<The data type of the matrix  \see MatrixVector_DataTypes 
    INTEGER(INTG) :: storageType !<The storage type of the matrix \see MatrixVector_StorageTypes 
    INTEGER(INTG) :: symmetryType !<The symmetry type of the matrix \see MatrixVector_SymmetryTypes
    INTEGER(INTG) :: transposeType !<The transpose type of the matrix \see MatrixVector_TransposeTypes
    INTEGER(INTG) :: numberOfNonZeros !<The number of non-zero elements in the matrix 
    INTEGER(INTG) :: size !<The size of the data arrays
    INTEGER(INTG) :: maximumColumnIndicesPerRow !<The maximum number of column indicies for the rows.
    INTEGER(INTG) :: maximumRowIndicesPerColumn !<If the matrix transpose is required the maximum number of row indicies for the columns.
    INTEGER(INTG) :: blockSize !<The block size for blocked storage types.
    INTEGER(INTG) :: numberOfBlocks !<The number of non-zero blocks
    INTEGER(INTG) :: numberOfRowBlocks !<The number of blocks in the row direction
    INTEGER(INTG) :: numberOfColumnBlocks !<The number of blocks in the column direction.
    INTEGER(INTG), ALLOCATABLE :: rowIndices(:) !<rowIndices(i). The row indices for the matrix storage scheme. \see MatrixVector_MatrixStorageStructures
    INTEGER(INTG), ALLOCATABLE :: columnIndices(:) !<columnIndices(i). The column indices for the matrix storage scheme. \see MatrixVector_MatrixStorageStructures
    INTEGER(INTG), ALLOCATABLE :: transposeRowsColumns(:) !<transposeRowsColumns(i). If the matrix has the transpose set for a number of columns then transposeRowsColumns(i) stores the i'th row of column in the transpose. \see MatrixVector_TransposeTypes
    INTEGER(INTG), ALLOCATABLE :: rowIndicesT(:) !<rowIndicesT(i). The transpose row indices for the matrix storage scheme. \see MatrixVector_MatrixStorageStructures
    INTEGER(INTG), ALLOCATABLE :: columnIndicesT(:) !<columnIndicesT(i). The transpose column indices for the matrix storage scheme. \see MatrixVector_MatrixStorageStructures
    INTEGER(INTG), ALLOCATABLE :: transposeDataSwivel(:) !transposeDataSwivel(i). The mapping from the transpose non-zero location to the corresponding non-zero location in the non-tranposed data array.
    INTEGER(INTG), ALLOCATABLE :: rowColStart(:) !<rowColStart(i). If the matrix transpose is required then rowColStart(i) stores the index into the rowColList which is the start of the i'th row or column. \see MatrixVector_MatrixStorageStructures \see MatrixVector_TransposeTypes
    INTEGER(INTG), ALLOCATABLE :: rowColList(:) !<rowColList(i). If the matrix transpose is required then the rowColList(j) stores the next entry in the row or column. If the entry is negative k then that indicates that is the end of the kth row or column. \see MatrixVector_MatrixStorageStructures \see MatrixVector_TransposeTypes
    TYPE(LINKEDLIST),POINTER :: list(:) !\todo Comment
    INTEGER(INTG), ALLOCATABLE :: dataIntg(:) !<dataIntg(i). The integer data for an integer matrix. The i'th component contains the data for the i'th matrix data stored on the domain.
    REAL(SP), ALLOCATABLE :: dataSP(:) !<dataSP(i). The real data for a single precision matrix. The i'th component contains the data for the i'th matrix data stored on the domain.
    REAL(DP), ALLOCATABLE :: dataDP(:) !<dataDP(i). The real data for a double precision matrix. The i'th component contains the data for the i'th matrix data stored on the domain.
    LOGICAL, ALLOCATABLE :: dataL(:) !<dataL(i). The logical data for a logical matrix. The i'th component contains the data for the i'th matrix data stored on the domain.
  END TYPE MatrixType

  !>Type to allow for coupling of equations/interface rows/columns to a solver row/column
  TYPE MatrixRowColCouplingType
    INTEGER(INTG) :: numberOfRowCols !<The number of rows/columns in the coupling
    INTEGER(INTG), ALLOCATABLE :: rowCols(:) !<rowCols(i). The row/column for the i'th row/column that is mapped.
    REAL(DP), ALLOCATABLE :: couplingCoefficients(:) !<couplingCoefficients(i). The coupling coefficient for the i'th row/column that is mapped.
  END TYPE MatrixRowColCouplingType
  
  PUBLIC VectorType

  PUBLIC MatrixType

  PUBLIC MatrixRowColCouplingType
  
  !
  !================================================================================================================================
  !
  ! Domain mapping types
  !
  
  !>Contains the information on an adjacent domain to a domain in a domain mapping. 
  TYPE DomainAdjacentDomainType
    INTEGER(INTG) :: domainNumber !<The number of the domain that is adjacent to a domain in a mapping.
    INTEGER(INTG) :: numberOfSendGhosts !<The number of ghost locals in the current domain that are to be sent to this domain.
    INTEGER(INTG) :: numberOfReceiveGhosts !<The number of ghost locals in the current domain that are to be received from this domain.
    INTEGER(INTG), ALLOCATABLE :: localGhostSendIndices(:) !<localGhostSendIndices(i). The local numbers of the ghosts in the current domain that are to be sent to this domain.
    INTEGER(INTG), ALLOCATABLE :: localGhostReceiveIndices(:) !<localGhostReceiveIndices(i). The local numbers of the ghosts in the current domain that are to be received from this domain.
  END TYPE DomainAdjacentDomainType
  
  !>Contains the local information for a global mapping number for a domain mapping.
  TYPE DomainGlobalMappingType
    INTEGER(INTG) :: numberOfDomains !<The number of domains that the global number is mapped to a local number in.
    INTEGER(INTG), ALLOCATABLE :: localNumber(:) !<localNumber(domainIdx). The mapped local number for the domain_Ix'th domain for the global number.
    INTEGER(INTG), ALLOCATABLE :: domainNumber(:) !<domainNumber(domainIdx). The domain number for the domainIdx'th domain for which the global number is mapped to a local number
    INTEGER(INTG), ALLOCATABLE :: localType(:) !<localType(domainIdx). The type of local for the domain_idx'th domain for which the global number is mapped to a local number. The types depend on whether the mapped local number in the domainIdx'th domain is an internal, boundary or ghost local number. \see DomainMappings_DomainType
  END TYPE DomainGlobalMappingType

  !>Contains information on the domain mappings (i.e., local and global numberings).
  TYPE DomainMappingType
    TYPE(WorkGroupType), POINTER :: workGroup !<A pointer to the work group for the domain mapping
    INTEGER(INTG) :: numberOfLocal !<The number of local numbers in the domain excluding ghost numbers
    INTEGER(INTG) :: totalNumberOfLocal !<The total number of local numbers in the domain including ghost numbers.
    INTEGER(INTG), ALLOCATABLE :: numberOfDomainLocal(:) !<numberOfDomainLocal(domainNumber). The number of locals for domainNumber'th domain. NOTE: the domain_no goes from 0 to the number of domains-1.
    INTEGER(INTG), ALLOCATABLE :: numberOfDomainGhost(:) !<numberOfDomainGhost(domainNumber). The total number of ghosts for domainNumber'th domain. NOTE: the domain_no goes from 0 to the number of domains-1.
    INTEGER(INTG) :: numberOfGlobal !<The number of global numbers for this mapping.
    !INTEGER(INTG) :: numberOfDomains !<The number of domains in this mapping.
    INTEGER(INTG) :: numberOfInternal !<The number of internal numbers in this mapping.
    INTEGER(INTG) :: numberOfBoundary !<The number of boundary numbers in this mapping.
    INTEGER(INTG) :: numberOfGhost !<The number of ghost numbers in this mapping.
    INTEGER(INTG) :: internalStart !<The start postition in the domainList for the list of internal numbers
    INTEGER(INTG) :: internalFinish !<The finish postition in the domainList for the list of internal numbers
    INTEGER(INTG) :: boundaryStart !<The start postition in the domainList for the list of boundary numbers
    INTEGER(INTG) :: boundaryFinish !<The finish postition in the domainList for the list of boundary numbers
    INTEGER(INTG) :: ghostStart !<The start postition in the domainList for the list of ghost numbers
    INTEGER(INTG) :: ghostFinish !<The finish postition in the domainList for the list of ghost numbers
    INTEGER(INTG), ALLOCATABLE :: domainList(:) !<domainList(i). The list of local numbers grouped so that the internal numbers are from internalStart to internalFinish, the boundary numbers are from boundaryStart to boundaryFinish and the ghost numbers are from ghostStart to ghostFinish
    INTEGER(INTG), ALLOCATABLE :: localToGlobalMap(:) !<localToGlobalMap(i). The global number for the i'th local number for the mapping.
    TYPE(DomainGlobalMappingType), ALLOCATABLE :: globalToLocalMap(:) !<globalToLocalMap(i). The local information for the i'th global number for the mapping.
    INTEGER(INTG) :: numberOfAdjacentDomains !<The number of domains that are adjacent to this domain in the mapping.
    INTEGER(INTG), ALLOCATABLE :: adjacentDomainsPtr(:) !<adjacentDomainsPtr(domainNumber). The pointer to the list of adjacent domains for domain_no. adjacentDomainsPtr(domain_no) gives the starting position in adjacentDomainsList for the first adjacent domain number for domain number domain_no. adjacentDomainsPtr(domain_no+1) gives the last+1 position in adjacentDomainsList for the last adjacent domain number for domain number domain_no. NOTE: the index for adjacentDomainsPtr varies from 0 to the number of domains.
    INTEGER(INTG), ALLOCATABLE :: adjacentDomainsList(:) !<adjacentDomainsList(i). The list of adjacent domains for each domain. The start and end positions for the list for domain number domainNumber are given by adjacentDomainPtr(domainNumber) and adjacentDomainPtr(domainNumber+1)-1 respectively.
    TYPE(DomainAdjacentDomainType), ALLOCATABLE :: adjacentDomains(:) !<adjacentDomains(adjacentDomainIdx). The adjacent domain information for the adjacentDomainIdx'th adjacent domain to this domain. 
  END TYPE DomainMappingType

  !>Contains information on the domain decomposition mappings.
  TYPE DomainMappingsType
    TYPE(DomainType), POINTER :: domain !<A pointer to the domain decomposition.
    TYPE(DomainMappingType), POINTER :: elements !<Pointer to the element mappings for the domain decomposition.
    TYPE(DomainMappingType), POINTER :: nodes !<Pointer to the node mappings for the domain decomposition.
    TYPE(DomainMappingType), POINTER :: dofs !<Pointer to the dof mappings for the domain decomposition.
  END TYPE DomainMappingsType
  
  !>Contains information for a domain.
  TYPE DomainType
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the domain decomposition for this domain.
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh for this domain.
    TYPE(RegionType), POINTER :: region !<A pointer to the region that this domain is in. This is "inherited" from the mesh region
    INTEGER(INTG) :: meshComponentNumber !<The mesh component number of the mesh which this domain was decomposed from.        
    INTEGER(INTG) :: numberOfDimensions !<The number of dimensions for this domain. This is "inherited" from the mesh.
    INTEGER(INTG), ALLOCATABLE :: nodeDomain(:) !<nodeDomain(nodeIdx). The domain number that the nodeIdx'th global node is in for the domain decomposition. Note: the domain numbers start at 0 and go up to the numberOfDomains-1.
    TYPE(DomainMappingsType), POINTER :: mappings !<Pointer to the mappings for the domain  e.g., global to local and local to global maps
    TYPE(DomainTopologyType), POINTER :: topology !<Pointer to the topology for the domain.
  END TYPE DomainType

  !>A buffer type to allow for an array of pointers to a DOMAIN_TYPE.
  TYPE DomainPtrType
    TYPE(DomainType), POINTER :: ptr !<The pointer to the domain.
  END TYPE DomainPtrType

  PUBLIC DomainAdjacentDomainType

  PUBLIC DomainGlobalMappingType,DomainMappingType,DomainMappingsType

  PUBLIC DomainType,DomainPtrType

  !
  !================================================================================================================================
  !
  ! Decomposition types
  !
  
  !>Contains the information for a line in a decomposition.
  TYPE DecompositionLineType
    INTEGER(INTG) :: number !<The line number in the decomposition.
    INTEGER(INTG) :: xiDirection!<The Xi direction of the line. Old CMISS name NPL(1,0,nl)
    INTEGER(INTG) :: numberOfSurroundingElements !<The number of elements that surround (use) this line.
    INTEGER(INTG), ALLOCATABLE :: surroundingElements(:) !<surroundingElements(nel). The local element number of the nel'th element that surrounds (uses) this line. 
    INTEGER(INTG), ALLOCATABLE :: elementLines(:) !<elementLines(nel). The local arc number of the nel'th element that surrounds (uses) this line.
    INTEGER(INTG) :: adjacentLines(0:1) !<adjacentLines(0:1). The line number of adjacent lines. adjacentLines(0) is the line number adjacent in the -xi direction. adjacentLines(1) is the line number adjacent in the +xi direction. Old CMISS name NPL(2..3,0,nl).
    LOGICAL :: boundaryLine !<Is .TRUE. if the line is on the boundary of the mesh for the domain, .FALSE. if not.
  END TYPE DecompositionLineType

  !>Contains the topology information for the lines of a decomposition.
  TYPE DecompositionLinesType
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<The pointer back to the decomposition topology.
    INTEGER(INTG) :: numberOfLines !<The number of lines in this decomposition topology.
    TYPE(DecompositionLineType), ALLOCATABLE :: lines(:) !<(lineIdx). The pointer to the array of topology information for the lines of this decomposition. lines(lineIdx) contains the topological information for the lineIdx'th local line of the decomposition.
  END TYPE DecompositionLinesType

  !>Contains the information for a face in a decomposition.
  TYPE DecompositionFaceType
    INTEGER(INTG) :: number !<The face number in the decomposition.
    INTEGER(INTG) :: xiNormalDirection !<The Xi direction of the face, the direction of the normal to the face
    INTEGER(INTG) :: numberOfSurroundingElements !<The number of elements that surround (use) this face.
    INTEGER(INTG), ALLOCATABLE :: surroundingElements(:) !<surroundingElements(nel). The local element number of the nel'th element that surrounds (uses) this face. 
    INTEGER(INTG), ALLOCATABLE :: elementFaces(:) !<elementFaces(nel). The local arc number of the nel'th element that surrounds (uses) this face.
!    INTEGER(INTG) :: adjacentFaces(0:1) !<adjacentFaces(0:1). The face number of adjacent faces. adjacentFaces(0) is the face number adjacent in the -xi direction. adjacentFaces(1) is the face number adjacent in the +xi direction. Old CMISS name NPL(2..3,0,nl).
    LOGICAL :: boundaryFace !<Is .TRUE. if the face is on the boundary of the mesh for the domain, .FALSE. if not.
    INTEGER(INTG) :: elementNumber !<The element number of the element on which the face is on
  END TYPE DecompositionFaceType

  !>Contains the topology information for the faces of a decomposition.
  TYPE DecompositionFacesType
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<The pointer back to the decomposition topology.
    INTEGER(INTG) :: numberOfFaces !<The number of faces in this decomposition topology.
    TYPE(DecompositionFaceType), ALLOCATABLE :: faces(:) !<faces(faceIdx). The pointer to the array of topology information for the faces of this decomposition. faces(faceIdx) contains the topological information for the faceIdx'th local face of the decomposition.
  END TYPE DecompositionFacesType

  !>Contains information on the decomposition adjacent elements for a xi coordinate
  TYPE DecompositionAdjacentElementType
    INTEGER(INTG) :: numberOfAdjacentElements !<The number of adjacent elements for the xi coordinate
    INTEGER(INTG), ALLOCATABLE :: adjacentElements(:) !<The local element numbers of the elements adjacent to this element for the xi coordinate
  END TYPE DecompositionAdjacentElementType
  
  !>Contains the information for an element in a decomposition.
  TYPE DecompositionElementType
    INTEGER(INTG) :: localNumber !<The local element number in the decomposition.
    INTEGER(INTG) :: globalNumber !<The corresponding global element number in the mesh of the local element number in the decomposition.
    INTEGER(INTG) :: userNumber !<The corresponding user number for the element.
    TYPE(DecompositionAdjacentElementType), ALLOCATABLE :: adjacentElements(:) !<adjacentElements(-xiCoordinateIdx:xiCoordinateIdx). The adjacent elements information in the xiCoordinateIdx'th xi coordinate direction. Note that -xiCoordinateIdx gives the adjacent element before the element in the xiCoordinateIdx'th direction and +xiCoordinateIdx gives the adjacent element after the element in the xiCoordinateIdx'th direction. The xiIdx=0 index will give the information on the current element. Old CMISS name NXI(-xiIdx:xiIdx,0:nei,ne).
    !\todo INTEGER(INTG), ALLOCATABLE :: numberOfAdjacentElements(:) !<numberOfAdjacentElements(-xiIdx:xiIdx). The number of elements adjacent to this element in the xiIdx'th xi direction. Note that -xiIdx gives the adjacent element before the element in the xiIdx'th direction and +xiIdx gives the adjacent element after the element in the xiIdx'th direction. The xiIdx=0 index should be 1 for the current element. Old CMISS name NXI(-xiIdx:xiIdx,0:nei,ne).
    !\todo INTEGER(INTG), ALLOCATABLE :: adjacentElements(:,:) !<adjacentElements(nei,-xiIdx:xiIdx). The local element numbers of the elements adjacent to this element in the xiIdx'th xi direction. Note that -xiIdx gives the adjacent elements before the element in the xiIdx'th direction and +xiIdx gives the adjacent elements after the element in the xiIdx'th direction. The xiIdx=0 index should give the current element number. Old CMISS name NXI(-xiIdx:xiIdx,0:nei,ne).
    INTEGER(INTG), ALLOCATABLE :: elementLines(:) !<elementLines(elementLineIdx). The local decomposition line number corresponding to the elementLineIdx'th local line of the element. Old CMISS name NLL(nae,ne). 
    INTEGER(INTG), ALLOCATABLE :: elementFaces(:) !<elementFaces(elementFaceIdx). The local decomposition face number corresponding to the elementFaceIdx'th local face of the element. Old CMISS name NLL(nae,ne). 
    LOGICAL :: boundaryElement !<Is .TRUE. if the element is on the boundary of the mesh for the domain, .FALSE. if not.
  END TYPE DecompositionElementType

  !>Contains the topology information for the elements of a decomposition.
  TYPE DecompositionElementsType
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<The pointer back to the decomposition topology.
    INTEGER(INTG) :: numberOfElements !<The number of elements excluding ghost elements in this decomposition topology.
    INTEGER(INTG) :: totalNumberOfElements !<The total number of elements in this decomposition topology.
    INTEGER(INTG) :: numberOfGlobalElements !<The number of global elements in this decomposition topology.
    TYPE(DecompositionElementType), ALLOCATABLE :: elements(:) !<elements(elementIdx). The pointer to the array of topology information for the elements of this decomposition. elements(elementIdx) contains the topological information for the elementIdx'th local element of the decomposition. 
    TYPE(TreeType), POINTER :: elementsTree !<A tree mapping the decomposition local element number to the decomposition user element number.
  END TYPE DecompositionElementsType
  
  !>Contains data point information
  TYPE DecompositionElementDataPointType
    INTEGER(INTG) :: userNumber !<User number of the projected data point
    INTEGER(INTG) :: globalNumber !<Global number of the data point, sequence is according to element number sequence  
    INTEGER(INTG) :: localNumber !<local number of the data point
  END TYPE DecompositionElementDataPointType
  
  !>Contains information on the projected data points on an element, for decomposition since data points on elements go with the elements
  TYPE DecompositionElementDataPointsType
    INTEGER(INTG) :: numberOfProjectedData !<Number of projected data on this element
    INTEGER(INTG) :: globalElementNumber !<The global number of this element (element index is local number and element number is global number)
    TYPE(DecompositionElementDataPointType), ALLOCATABLE :: dataIndices(:) !<dataIndices(dataPointIdx). The global,local and user number of this data point
  END TYPE DecompositionElementDataPointsType
  
  !>Contains data point decompostion topology   
  TYPE DecompositionDataPointsType
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<The pointer back to the decomposition topology.
    INTEGER(INTG) :: numberOfDataPoints !<The number of data points excluding ghost data points in this decomposition topology.
    INTEGER(INTG) :: totalNumberOfDataPoints !<Number of projected data in this decomposition topology.
    INTEGER(INTG) :: numberOfGlobalDataPoints !<Number of global projected data 
    INTEGER(INTG), ALLOCATABLE :: numberOfDomainLocal(:) !<numberOfDomainLocal(compDomainIdx). Number of local data points in each computation domain 
    INTEGER(INTG), ALLOCATABLE :: numberOfDomainGhost(:) !<numberOfDomainGhost(compDomainIdx). Number of ghost data points in each computation domain 
    INTEGER(INTG), ALLOCATABLE :: numberOfElementDataPoints(:) !<numberOfElementDataPoints(elementIdx). Number of data points in each global element
    TYPE(DecompositionElementDataPointsType), ALLOCATABLE :: elementDataPoints(:) !<elementDataPoints(elementIdx). Information of the projected data on the elements for decomposition of data points
    TYPE(TreeType), POINTER :: dataPointsTree !<A tree mapping the domain local number to the region data point user number.
  END TYPE DecompositionDataPointsType

   !>Contains the topology information for a decomposition
  TYPE DecompositionTopologyType
    TYPE(DecompositionType), POINTER :: decomposition !<The pointer to the decomposition for this topology information.
    TYPE(DecompositionElementsType), POINTER :: elements !<The pointer to the topology information for the elements of this decomposition.
    TYPE(DecompositionLinesType), POINTER :: lines !<The pointer to the topology information for the lines of this decomposition.
    TYPE(DecompositionFacesType), POINTER :: faces !<The pointer to the topology information for the faces of this decomposition.
    TYPE(DecompositionDataPointsType), POINTER :: dataPoints !<The pointer to the topology information for the data of this decomposition.
  END TYPE DecompositionTopologyType

  !>Contains information on the decomposition connectivity for a given coupled mesh element
  TYPE DecompositionElementConnectivityType
    INTEGER(INTG) :: coupledElementNumber !< The coupled number to define the connectivity for.
    REAL(DP), ALLOCATABLE :: xi(:,:,:) !<xi(xiIdx,decompositionComponent,elementParameterIdx) The xiIdx'th xi of a coupled element to be copuled to a decompositionsh's decompositionMeshComponentIdx'th decomposition components's elementParameterIdx'th element parameter. !\todo the XI array needs to be restructured to efficiently utilize memory when coupling bodies with 2xi directions to bodies with 3xi directions using an interface.
    INTEGER(INTG) :: connectedLineFace !<The coupled element line/face number to be connected to the decomposition.
  END TYPE DecompositionElementConnectivityType

  !>A buffer type to allow for arrays of pointers to DecompositionElementConnectivityType
  TYPE DecompositionElementConnectivityPtrType
    TYPE(DecompositionElementConnectivityType), POINTER :: ptr !<A pointer to DecompositionElementConnectivityType
  END TYPE DecompositionElementConnectivityPtrType

  !>Contains information on the coupling between decompositions in an interface
  TYPE DecompositionConnectivityType
    TYPE(DecompositionCouplingType), POINTER :: decompositionCoupling !<A pointer back to the decomposition coupling
    TYPE(BasisType), POINTER :: basis !<A pointer to the inteface basis
    INTEGER(INTG) :: numberOfInterfaceElements !<The number of elements in the interface
    INTEGER(INTG) :: totalNumberOfInterfaceElements !<The total number of elements in the interface
    INTEGER(INTG) :: numberOfGlobalInterfaceElements !<The global number of elements in the interface
    INTEGER(INTG) :: numberOfCoupledDecompositions !<The number of coupled decompositions in the interface
    TYPE(DecompositionElementConnectivityPtrType), ALLOCATABLE :: elementConnectivity(:,:) !<elementConnectivity(elementIdx,coupledDecompositionIdx) !<The decomposition connectivity for a given interface element
    INTEGER(INTG), ALLOCATABLE :: coupledNodes(:,:) !<coupledNodes(coupledMeshIdx,interfaceNodeIdx). Coupled nodes numbers
  END TYPE DecompositionConnectivityType
  
  !>Contains information on the mesh decomposition. \see OpenCMISS::Iron::cmfe_DecompositionType
  TYPE DecompositionType
    INTEGER(INTG) :: userNumber!<The user defined identifier for the domain decomposition. The user number must be unique.
    INTEGER(INTG) :: globalNumber !<The global number of the domain decomposition in the list of domain decompositions for a particular mesh.
    LOGICAL :: decompositionFinished !<Is .TRUE. if the decomposition has finished being created, .FALSE. if not.
    TYPE(DecompositionsType), POINTER :: decompositions !<A pointer to the decompositions for this decomposition.
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh for this decomposition.
    TYPE(RegionType), POINTER :: region !<A pointer to the region containing the mesh. If the mesh is in an interface rather than a region then this pointer will be NULL and the interface pointer should be used.
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface containing the mesh. If the mesh is in a region rather than an interface then this pointer will be NULL and the interface pointer should be used.
    TYPE(DecomposerType), POINTER :: decomposer !<A pointer to the decomposer for this decomposition.
    INTEGER(INTG) :: numberOfDimensions !<The number of dimensions (Xi directions) for this decomposition/mesh.
    INTEGER(INTG) :: numberOfComponents !<The number of mesh components in this decomposition/mesh.
    INTEGER(INTG) :: meshComponentNumber !<The component number (index) of the mesh component that this decomposition belongs to (i.e., was generated from).
    INTEGER(INTG) :: domainDecompositionType !<The type of the domain decomposition \see MeshRoutines_DecompositionTypes.
    TYPE(WorkGroupType), POINTER :: workGroup !<The work group to use for this decomposition
    INTEGER(INTG) :: numberOfElements !<The number of elements in the decomposition
    INTEGER(INTG), ALLOCATABLE :: elementDomain(:) !<elementDomain(elementIdx). The domain number that the elementIdx'th global element is in for the decomposition. Note: the domain numbers start at 0 and go up to the numberOfDomains-1.
    TYPE(DecompositionTopologyType), POINTER :: topology !<A pointer to the topology for this decomposition.
    TYPE(DomainPtrType), ALLOCATABLE :: domain(:) !<domain(meshComponentIdx). A pointer to the domain for mesh component for the domain associated with the computation node.
    TYPE(DecompositionCouplingType), POINTER :: decompositionCoupling !<A pointer to the decomposition coupling if any.
    LOGICAL :: calculateFaces !<Boolean flag to determine whether faces should be calculated
    LOGICAL :: calculateLines !<Boolean flag to determine whether lines should be calculated
  END TYPE DecompositionType

  !>A buffer type to allow for an array of pointers to a DecompositionType.
  TYPE DecompositionPtrType
    TYPE(DecompositionType), POINTER :: ptr !<The pointer to the domain decomposition. 
  END TYPE DecompositionPtrType

  !>Contains information on coupling decompositions (across interfaces)
  TYPE DecompositionCouplingType
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer back to the decomposition
    INTEGER(INTG) :: numberOfCoupledDecompositions !<The number of decompositions coupled to this decomposition
    TYPE(DecompositionPtrType), ALLOCATABLE :: coupledDecompositions(:) !<coupledDecompositions(coupledDecompositionIdx). The array of pointers to the coupled decompositions
    TYPE(DecompositionConnectivityType), POINTER :: decompositionConnectivity !<A pointer to information on the decompositon node element connectivity
  END TYPE DecompositionCouplingType

  !>Contains information on the domain decompositions defined on a mesh.
  TYPE DecompositionsType
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh.
    INTEGER(INTG) :: numberOfDecompositions !<The number of decompositions defined on the mesh.
    TYPE(DecompositionPtrType), ALLOCATABLE :: decompositions(:) !<decompositions(decompositionIdx). The array of pointers to the domain decompositions.
  END TYPE DecompositionsType

  !>Contains information on a link between nodes in a decomposer graph
  TYPE DecomposerGraphLinkType
    TYPE(DecomposerGraphNodeType), POINTER :: parentGraphNode !<A pointer back to the parent node for the link
    TYPE(DecomposerGraphNodeType), POINTER :: linkedGraphNode !<A pointer forward to the linked node for this link
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the interface decomposition for this link
  END TYPE DecomposerGraphLinkType

  !>A buffer type to allow for an array of pointers to DecomposerGraphLinkType
  TYPE DecomposerGraphLinkPtrType
    TYPE(DecomposerGraphLinkType), POINTER :: ptr !<The pointer to the DecomposerGraphLinkType
  END TYPE DecomposerGraphLinkPtrType

  !>Contains information on a node in a decomposer graph
  TYPE DecomposerGraphNodeType
    TYPE(DecomposerGraphType), POINTER :: decomposerGraph !<A pointer back to the decomposer graph
    TYPE(DecompositionType), POINTER :: decomposition !<The decomposition for the graph node
    LOGICAL :: rootNode !<Is .TRUE. if the node is a root of the graph, .FALSE. if not. 
    INTEGER(INTG) :: numberOfGraphLinks !<The number of links to other decomposer graph nodes
    TYPE(DecomposerGraphLinkPtrType), ALLOCATABLE :: graphLinks(:) !<graphLinks(graphLinkIdx). graphLinks(graphLinkIdx)%ptr points to the graphLinkIdx'th link to another node in the decomposition graph.
    INTEGER(INTG) :: currentLinkIndex !<The current link index whilst the decomposer graph is being traversed.
    TYPE(DecomposerGraphNodeType), POINTER :: previousNode !<The previous node whilst the decomposer graph is being traversed.
    INTEGER(INTG) :: elementOffset !<The element number offset for this node in the supergraph.
    INTEGER(INTG) :: vertexOffset !<The vertex number offset for this node in the supergraph.
  END TYPE DecomposerGraphNodeType

  !>A buffer type to allow for an array of pointers to DecomposerGraphNodeType
  TYPE DecomposerGraphNodePtrType
    TYPE(DecomposerGraphNodeType), POINTER :: ptr !<The pointer to the DecomposerGraphNodeType
  END TYPE DecomposerGraphNodePtrType

  !>Contains information on the decomposer graph of connected decompositions to decompose.
  TYPE DecomposerGraphType
    TYPE(DecomposerType), POINTER :: decomposer !<A pointer back to the decomposer.
    INTEGER(INTG) :: numberOfGraphRoots !<The number of graph root pointers. A graph root points to a graph of connected decompositions
    TYPE(DecomposerGraphNodePtrType), ALLOCATABLE :: graphRoots(:) !<graphRoots(graphRootIdx). graphRoots(graphRootIdx)%ptr points to the root of the graphRootIdx'th connected decompositions graph 
  END TYPE DecomposerGraphType

  !>Contains information of a decomposer.
  TYPE DecomposerType
    INTEGER(INTG) :: userNumber !<The user number of the region decomposer
    TYPE(DecomposersType), POINTER :: decomposers !<A pointer back to the region decomposers
    TYPE(RegionType), POINTER :: region !<A pointer back to the region for the decomposer
    LOGICAL :: decomposerFinished !<Is .TRUE. if the region decomposer has been finished
    INTEGER(INTG) :: numberOfDecompositions !<The number of decompositions in the region decomposer
    TYPE(DecompositionPtrType), ALLOCATABLE :: decompositions(:) !<decompositions(decompositionIdx)%ptr is the pointer to the decompositionIdx'th decomposition in the decomposer
    TYPE(DecomposerGraphType), POINTER :: decomposerGraph !<A pointer to the decomposer graph of decompositions.
    TYPE(WorkGroupType), POINTER :: workGroup
    INTEGER(INTG) :: outputType !<The output type for the decomposer \see DecompositionRoutines_DecomposerOutputTypes,DecompositionRoutines
    INTEGER(INTG) :: numberOfPartitions !<The number of partitions for the decomposer
    INTEGER(INTG) :: numberOfEdgesCut !<The number of edges cut in the decomposition
  END TYPE DecomposerType

  !>A buffer type to allow for an array of pointers to DecomposerType
  TYPE DecomposerPtrType
    TYPE(DecomposerType), POINTER :: ptr !<A pointer to the region decomposer.
  END TYPE DecomposerPtrType

  !>Contains information on region decomposers defined on a region.
  TYPE DecomposersType
    TYPE(RegionType), POINTER :: region !<A pointer back to the region for the region decomposers.
    INTEGER(INTG) :: numberOfDecomposers !<The number of region decomposers.
    TYPE(DecomposerPtrType), ALLOCATABLE :: decomposers(:) !<decomposers(decomposeridx). The array of pointers to the region decomposers. 
  END TYPE DecomposersType

  PUBLIC DecompositionLineType,DecompositionLinesType

  PUBLIC DecompositionFaceType,DecompositionFacesType

  PUBLIC DecompositionAdjacentElementType,DecompositionElementType,DecompositionElementsType

  PUBLIC DecompositionElementDataPointType,DecompositionElementDataPointsType,DecompositionDataPointsType

  PUBLIC DecompositionTopologyType

  PUBLIC DecompositionCouplingType,DecompositionConnectivityType

  PUBLIC DecompositionElementConnectivityType,DecompositionElementConnectivityPtrType

  PUBLIC DecompositionType,DecompositionPtrType,DecompositionsType

  PUBLIC DecomposerGraphLinkType,DecomposerGraphLinkPtrType

  PUBLIC DecomposerGraphNodeType,DecomposerGraphNodePtrType

  PUBLIC DecomposerGraphType

  PUBLIC DecomposersType,DecomposerPtrType,DecomposerType

  !
  !================================================================================================================================
  !
  ! Field types
  !

  !>Contains information on a physical point in a field.
  TYPE FieldPhysicalPointType
    TYPE(FieldInterpolatedPointType), POINTER :: fieldInterpolatedPoint !<A pointer to the interpolated point of the field that is to be interpolated.
    TYPE(FieldInterpolatedPointType), POINTER :: geometricInterpolatedPoint !<A pointer to the interpolated point of the geometric field that is to be interpolated.
    INTEGER(INTG) :: physicalDerivativeType !<The type of the physical derivatives that have been interpolated. \see CONSTANTS_PhysicalDerivativeConstants
    REAL(DP), ALLOCATABLE :: values(:) !<values(componentIdx). The physical field component values.
  END TYPE FieldPhysicalPointType

  !>Buffer type to allow for arrays of pointers to FieldPhysicalPointType
  TYPE FieldPhysicalPointPtrType
    TYPE(FieldPhysicalPointType), POINTER :: ptr
  END TYPE FieldPhysicalPointPtrType
  
  !> Contains the interpolated point coordinate metrics. Old CMISS name GL,GU,RG.
  TYPE FieldInterpolatedPointMetricsType
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint !<A pointer to the interpolated point.
    INTEGER(INTG) :: numberOfXDimensions !<The number of X dimensions.
    INTEGER(INTG) :: numberOfXiDimensions !<The number of Xi dimensions
    REAL(DP), ALLOCATABLE :: gl(:,:) !<GL(mXiIdx,nXiIdx). Covariant metric tensor. Old CMISS name GL.
    REAL(DP), ALLOCATABLE :: gu(:,:) !<GU(mXiIdx,nXiIdx). Contravariant metric tensor. Old CMISS name GU.
    REAL(DP), ALLOCATABLE :: dXdXi(:,:) !<dXdXi(coordinateIdx,xiIdx). Rate of change of the X coordinate system wrt the Xi coordinate system.
    REAL(DP), ALLOCATABLE :: dXidX(:,:) !<dXidX(xiIdx,coordinateIdx). Rate of change of the Xi coordinate system wrt the X coordinate system. 
    REAL(DP) :: jacobian !<The Jacobian of the Xi to X coordinate system transformation. Old CMISS name RG.
    INTEGER(INTG) :: jacobianType !<The type of Jacobian. \see CoordinateRoutines_JacobianType
  END TYPE FieldInterpolatedPointMetricsType

  TYPE FieldInterpolatedPointMetricsPtrType
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: ptr
  END TYPE FieldInterpolatedPointMetricsPtrType
  
  !>Contains the interpolated value (and the derivatives wrt xi) of a field at a point. Old CMISS name XG.
  TYPE FieldInterpolatedPointType
    TYPE(FieldInterpolationParametersType), POINTER :: interpolationParameters !<A pointer to the interpolation parameters of the field that is to be interpolated.
    INTEGER(INTG) :: maximumPartialDerivativeIndex !<The maximum number of partial derivatives that have been allocated for the values component.
    INTEGER(INTG) :: partialDerivativeType !<The type of the partial derivatives that have been interpolated. partialDerivativeType can be either NO_PART_DERIV, FIRST_PART_DERIV or SECOND_PART_DERIV depending on whether just the field value, the field value and all first derivatives (including cross derivatives) or the first value and all first and second derivatives have been interpolated.
    REAL(DP), ALLOCATABLE :: values(:,:) !<values(componentIdx,partialDerivativeIdx). The interpolated field components and their partial derivatives.
  END TYPE FieldInterpolatedPointType

  TYPE FieldInterpolatedPointPtrType
    TYPE(FieldInterpolatedPointType), POINTER :: ptr
  END TYPE FieldInterpolatedPointPtrType
  
  !>Contains the parameters required to interpolate a field variable within an element. Old CMISS name XE
  TYPE FieldInterpolationParametersType
    TYPE(FieldType), POINTER :: field !<A pointer to the field to be interpolated.
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to be interpolated.
    INTEGER(INTG) :: numberOfXi !<The number of xi directions for the interpolation parameters.
    TYPE(BasisPtrType), ALLOCATABLE :: bases(:) !<bases(componentIdx). An array to hold a pointer to the basis (if any) used for interpolating the componentIdx'th component of the field variable.
    INTEGER(INTG), ALLOCATABLE :: numberOfParameters(:) !<numberOfParameters(componentIdx). The number of interpolation parameters used for interpolating the componentIdx'th component of the field variable.
    REAL(DP), ALLOCATABLE :: parameters(:,:) !<parameters(elementParameterIdx,componentIdx). The elementParameterIdx'th interpolation parameter used for interpolating the componentIdx'th component of the field variable.
    REAL(DP), ALLOCATABLE :: scaleFactors(:,:) !<scaleFactors(elementParameterIdx,componentIdx). The scale factors used for scaling the componentIdx'th component of the field variable.
  END TYPE FieldInterpolationParametersType

  TYPE FieldInterpolationParametersPtrType
    TYPE(FieldInterpolationParametersType), POINTER :: ptr
  END TYPE FieldInterpolationParametersPtrType
  
  !>Contains the geometric parameters (lines, faces, volumes etc.) for a geometric field decomposition.
  TYPE FieldGeometricParametersType
    INTEGER(INTG) :: numberOfLines !<The number of lines in the field. 
    INTEGER(INTG) :: numberOfAreas !<The number of areas in the field. 
    INTEGER(INTG) :: numberOfVolumes !<The number of volumes in the field. Inherited from the field decomposition.
    REAL(DP), ALLOCATABLE :: lengths(:) !<lengths(lineIdx). The length of the lineIdx'th line in the field decomposition.
    REAL(DP), ALLOCATABLE :: areas(:) !<areas(faceIdx). The area of the faceIdx'th face in the field decomposition.
    REAL(DP), ALLOCATABLE :: volumes(:) !<volumes(elementIdx). The volume of the elementIdx'th element in the field decomposition.
    INTEGER(INTG) :: numberOfFieldsUsing !<The number of fields that use these geometric parameters for their scaling. 
    TYPE(FieldPtrType), POINTER :: fieldsUsing(:) !< fieldsUsing(fieldIdx). A pointer to the fieldIdx'th field that uses these geometric parameters for its scaling.
  END TYPE FieldGeometricParametersType

  !>A type to hold the scale factors for the appropriate mesh component of a field. 
  TYPE FieldScalingType
    INTEGER(INTG) :: meshComponentNumber !<The mesh component number of a field variable component that the scaling factors are associated with.
    INTEGER(INTG) :: maximumNumberOfDerivatives !<The maximum number of derivatives in the mesh component. 
    INTEGER(INTG) :: maxNumberOfElementParameters !<The maximum number of element parameters in the mesh component. 
    TYPE(DistributedVectorType), POINTER :: scaleFactors !<scaleFactors(nk,np). The scale factor that is applied to the nk'th derivative of the np'th node of the mesh component. \todo  Make scale factors nodally based for now. Will have to revert to element based and extended to be a matrix to allow for a global derivative to be mapped onto many different element derivatives at the points that closes meshes or for inconsistent xi directions
  END TYPE FieldScalingType

  !>A type to hold the field scalings for the field.
  TYPE FieldScalingsType
    INTEGER(INTG) :: scalingType !<The type of scaling that is applied to the field. \see FieldRoutines_ScalingTypes
    INTEGER(INTG) :: numberOfScalingIndices !<The number of scaling indices (or sets of scale factors) for the field. In general there will be a set of scale factors (or a scaling index) for each different mesh component that is used by the field variable components.
    TYPE(FieldScalingType), ALLOCATABLE :: scalings(:) !<scalings(scalingIdx). The scaling factors for the scalingIdx'th set of scaling factors. 
  END TYPE FieldScalingsType

  !> A type to hold the mapping from field dof numbers to field parameters (nodes, elements, etc)
  TYPE FieldDOFToParamMapType
    INTEGER(INTG) :: numberOfDOFs !<The number of degrees-of-freedom for the field.
    INTEGER(INTG), ALLOCATABLE :: DOFType(:,:) !<DOFType(i=1..2,ny). The parameter type of the ny'th dof. When i=1 the DOFType is the type of the dof, i.e., 1=constant field dof, 2=element based field dof, 3=node based field dof, 4=point based field dof. When i=2 the DOFType gives the nyy'th number of the different field types and is used to index the xxxDOF2ParamMap arrays.
    INTEGER(INTG) :: numberOfConstantDOFs !<The number of constant degrees-of-freedom in the field dofs.
    INTEGER(INTG) :: numberOfElementDOFs !<The number of element based degrees-of-freedom in the field dofs.
    INTEGER(INTG) :: numberOfNodeDOFs !<The number of node based degrees-of-freedom in the field dofs.
    INTEGER(INTG) :: numberOfGridPointDOFs !<The number of grid point based degrees-of-freedom in the field dofs.
    INTEGER(INTG) :: numberOfGaussPointDOFs !<The number of Gauss point based degrees-of-freedom in the field dofs.
    INTEGER(INTG) :: numberOfDataPointDOFs !<The number of data point based degrees-of-freedom in the field dofs.
    INTEGER(INTG), ALLOCATABLE :: constantDOF2ParamMap(:) !<constantDOF2ParamMap(nyy). The mapping from constant field dofs to field parameters for the nyy'th constant field dof. The DOF2ParamMap gives the component number (nh) of the field parameter. The nyy value for a particular field dof (ny) is given by the DOFType component of this type.
    INTEGER(INTG), ALLOCATABLE :: elementDOF2ParamMap(:,:) !<elementDOF2ParamMap(i=1..2,nyy). The mapping from element based field dofs to field parameters for the nyy'th constant field dof. When i=1 the DOF2ParamMap gives the element number (ne) of the field parameter. When i=2 the DOF2ParamMap gives the component number (nh) of the field parameter. The nyy value for a particular field dof (ny) is given by the DOFType component of this type.
    INTEGER(INTG), ALLOCATABLE :: nodeDOF2ParamMap(:,:) !<nodeDOF2ParamMap(i=1..4,nyy). The mapping from node based field dofs to field parameters for the nyy'th constant field dof. When i=1 the DOF2ParamMap gives the version number (version_idx) of the field parameter. When i=2 the DOF2ParamMap gives the derivative number (derivative_idx) of the field parameter. When i=3 the DOF2ParamMap gives the node number (node_idx) of the field parameter. When i=4 the DOF2ParamMap gives the component number (componentIdx) of the field parameter. The nyy value for a particular field dof (ny) is given by the DOFType component of this type.
    INTEGER(INTG), ALLOCATABLE :: gridPointDOF2ParamMap(:,:) !<gridPointDOF2ParamMap(i=1..2,nyy). The mapping from grid point based field dofs to field parameters for the nyy'th grid point field dof. When i=1 the DOF2ParamMap gives the grid point number (nq) of the field parameter. When i=2 the DOF2ParamMap gives the component number (nh) of the field parameter. The nyy value for a particular field dof (ny) is given by the DOFType component of this type.  
    INTEGER(INTG), ALLOCATABLE :: gaussPointDOF2ParamMap(:,:) !<gaussPointDOF2ParamMap(i=1..3,nyy). The mapping from Gauss point based field dofs to field parameters for the nyy'th grid point field dof. When i=1 the DOF2ParamMap gives the Gauss point number (ng) of the field parameter. When i=2 the DOF2ParamMap gives the element number (ne) of the field parameter. When i=3 the DOF2ParamMap gives the component number (nh) of the field parameter. The nyy value for a particular field dof (ny) is given by the DOFType component of this type.  
    INTEGER(INTG), ALLOCATABLE :: dataPointDOF2ParamMap(:,:)!<dataPointDOF2ParamMap(i=1..3,nyy). The mapping from data point based field dofs to field parameters for the nyy'th data point field dof. When i=1 the DOF2ParamMap gives the data point number (np) of the field parameter. When i=2 the DOF2ParamMap gives the element number (ne) of the field parameter. When i=3 the DOF2ParamMap gives the component number (nh) of the field parameter. The nyy value for a particular field dof (ny) is given by the DOFType component of this type.  
  END TYPE FieldDOFToParamMapType

  !>A type to hold the mapping from a field node derivative's versions to field dof numbers for a particular field variable component.
  TYPE FieldNodeParamToDOFMapDerivativeType
    INTEGER(INTG) :: numberOfVersions !<The number of versions for the node derivative parameters of this field variable component.
    INTEGER(INTG), ALLOCATABLE :: versions(:) !<nodeParam2DOFMap%nodes%(nodeIdx)%derivatives(derivativeIdx)%versions(versionIdx). The field variable dof number of the nodeIdx'th node's derivativeIdx'th derivative's versionIdx'th version.
  END TYPE FieldNodeParamToDOFMapDerivativeType

  !>A type to hold the mapping from a field node's derivative to field dof numbers for a particular field variable component.
  TYPE FieldNodeParamToDOFMapNodeType
    INTEGER(INTG) :: numberOfDerivatives !<The number of derivatives for the node parameter of this field variable component.
    TYPE(FieldNodeParamToDOFMapDerivativeType), ALLOCATABLE :: derivatives(:) !<derivatives(?). The mapping from field node derivative parameter to a dof
  END TYPE FieldNodeParamToDOFMapNodeType

  !>A type to hold the mapping from field nodes to field dof numbers for a particular field variable component.
  TYPE FieldNodeParamToDOFMapType
    INTEGER(INTG) :: numberOfNodeParameters !<The number of node based field parameters for this field variable component.
    TYPE(FieldNodeParamToDOFMapNodeType), ALLOCATABLE :: nodes(:)  !<nodes(?). The mapping from field node parameter to a dof
  END TYPE FieldNodeParamToDOFMapType

  !>A type to hold the mapping from field elements to field dof numbers for a particular field variable component.
  TYPE FieldElementParamToDOFMapType
    INTEGER(INTG) :: numberOfElementParameters !<The number of element based field parameters for this field variable component.
    INTEGER(INTG), ALLOCATABLE :: elements(:) !<elementParam2DOFMap%elements(elementIdx). The field variable dof number of the elementIdx'th element based parameter for this field variable component. \todo Allow for multiple element parameters per element.
  END TYPE FieldElementParamToDOFMapType

  !>A type to hold the mapping from field grid points to field dof numbers for a particular field variable component.
  TYPE FieldGridPointParamToDOFMapType
    INTEGER(INTG) :: numberOfGridPointParameters !<The number of grid point based field parameters for this field variable component.
    INTEGER(INTG), ALLOCATABLE :: gridPoints(:) !<gridPointParam2DOFMap%gridPoints(gridPointIdx). The field variable dof number of the gridPointIdx'th point based parameter for this field variable component. 
  END TYPE FieldGridPointParamToDOFMapType

  !>A type to hold the mapping from field Gauss points to field dof numbers for a particular field variable component.
  TYPE FieldGaussPointParamToDOFMapType
    INTEGER(INTG) :: numberOfGridPointParameters !<The number of Gauss point based field parameters for this field variable component.
    INTEGER(INTG), ALLOCATABLE :: gaussPoints(:,:) !<gaussPointParam2DOFMap%gaussPoints(gaussPointIdx,elementIdx). The field variable dof number of the gaussPointIdx'th Gauss point in the elementIdx'th element based parameter for this field variable component. 
  END TYPE FieldGaussPointParamToDOFMapType
  
  !>A type to hold the mapping from field data points to field dof numbers for a particular field variable component.
  
  TYPE FieldDataPointParamToDOFMapType
    INTEGER(INTG) :: numberOfDataPointParameters !<The number of data point based field parameters for this field variable component.
    INTEGER(INTG), ALLOCATABLE :: dataPoints(:) !<dataPointParam2DOFMap%dataPoints(dataPointIdx). The field variable dof number of the dataPointIdx'th data point based parameter for this field variable component. 
  END TYPE FieldDataPointParamToDOFMapType

  !>A type to hold the mapping from field parameters (nodes, elements, etc) to field dof numbers for a particular field variable component.
  TYPE FieldParamToDOFMapType
    INTEGER(INTG) :: numberOfConstantParameters !<The number of constant field parameters for this field variable component. Note: this is currently always 1 but is included for completeness and to allow for multiple constants per field variable component in the future.
    INTEGER(INTG) :: constantParam2DOFMap !<The field variable dof number of the constant parameter for this field variable component.
    TYPE(FieldElementParamToDOFMapType) :: elementParam2DOFMap !> A type to hold the mapping from field element parameters to field dof numbers
    TYPE(FieldNodeParamToDOFMapType) :: nodeParam2DOFMap !> A type to hold the mapping from field node parameters to field dof numbers
    TYPE(FieldGridPointParamToDOFMapType) :: gridPointParam2DOFMap !> A type to hold the mapping from grid point element parameters to field dof numbers
    TYPE(FieldGaussPointParamToDOFMapType) :: gaussPointParam2DOFMap !> A type to hold the mapping from field gauss point parameters to field dof numbers
    TYPE(FieldDataPointParamToDOFMapType) :: dataPointParam2DOFMap !> A type to hold the mapping from field data point parameters to field dof numbers
  END TYPE FieldParamToDOFMapType

  !>Contains information for a component of a field variable.
  TYPE FieldVariableComponentType
    INTEGER(INTG) :: componentNumber !<The number of the field variable component.
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable for this component.
    TYPE(VARYING_STRING) :: componentLabel !<The label for the field variable component
    INTEGER(INTG) :: interpolationType !<The interpolation type of the field variable component \see FieldRoutines_InterpolationTypes
    INTEGER(INTG) :: meshComponentNumber !<The mesh component of the field decomposition for this field variable component.
    INTEGER(INTG) :: scalingIndex !<The index into the defined field scalings for this field variable component.
    TYPE(DomainType), POINTER :: domain !<A pointer to the domain of the field decomposition for this field variable component.
    INTEGER(INTG) :: maxNumberElementInterpolationParameters !<The maximum number of interpolation parameters in an element for a field variable component. 
    INTEGER(INTG) :: maxNumberNodeInterpolationParameters !<The maximum number of interpolation parameters in an element for a field variable component. 
    TYPE(FieldParamToDOFMapType) :: paramToDOFMap !<The mapping of the field parameters to the field dofs for this field variable component.
  END TYPE FieldVariableComponentType

  !>A type to hold the parameter sets for a field.
  TYPE FieldParameterSetType
    INTEGER(INTG) :: setIndex !<The global set index (from 1 to the Types::FieldParameterSetsType::numberOfParameterSets) that this parameter set corresponds to.
    INTEGER(INTG) :: setType !<The user set type (index) (from 1 to FieldRoutines::FIELD_NUMBER_OF_SET_TYPES) that this parameter set \see FieldRoutines_ParameterSetTypes
  !###      corresponds to.
    TYPE(DistributedVectorType), POINTER :: parameters !<A pointer to the distributed vector that contains the field parameters for this field parameter set.
  END TYPE FieldParameterSetType
  
  !>A buffer type to allow for an array of pointers to a FieldParameterSetType.
  TYPE FieldParameterSetPtrType
    TYPE(FieldParameterSetType), POINTER :: ptr !<The pointer to the field parameter set. 
  END TYPE FieldParameterSetPtrType

  !>A type to store the parameter sets for a field.
  TYPE FieldParameterSetsType    
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable that these parameter sets are defined on.
    INTEGER(INTG) :: numberOfParameterSets !<The number of parameter sets that are currently defined on the field.
    TYPE(FieldParameterSetPtrType), POINTER :: setType(:) !<setType(setTypeIdx). A pointer to an array of pointers to the field set types. setType(setTypeIdx)%ptr is a pointer to the parameter set type for the setTypeIdx'th parameter set. setTypeIdx can vary from 1 to FieldRoutines::FIELD_NUMBER_OF_SET_TYPES. The value of the pointer will be NULL if the parameter set corresponding to the setTypeIdx'th parameter set has not yet been created for the field.
    TYPE(FieldParameterSetPtrType), POINTER :: parameterSets(:) !<parameterSets(setTypeIdx). \todo change to allocatable. A pointer to an array of pointers to the parameter sets that have been created on the field. parameterSets(setTypeIdx)%ptr is a pointer to the parameter set type for the setTypeIdx'th parameter set that has been created. setTypeIdx can vary from 1 to the number of parameter set types that have currently been created for the field i.e., Types::FieldParameterSetsType::numberOfParameterSets.
  END TYPE FieldParameterSetsType

  !>Contains information for a field variable defined on a field.
  TYPE FieldVariableType
    INTEGER(INTG) :: variableNumber !<The number of the field variable
    INTEGER(INTG) :: variableType !<The type of the field variable. \see FieldRoutines_VariableTypes
    TYPE(VARYING_STRING) :: variableLabel !<The label for the variable
    TYPE(FieldType), POINTER :: field !<A pointer to the field for this field variable.
    TYPE(RegionType), POINTER :: region !<A pointer to the region for this field variable.
    INTEGER(INTG) :: dimension !<The dimension of the field variable. \see FieldRoutines_DimensionTypes
    INTEGER(INTG) :: dataType !<The data type of the field variable.  \see FieldRoutines_DataTypes,FieldRoutines
    INTEGER(INTG) :: dofOrderType !<The order of the DOF's in the field variable \see FieldRoutines_DOFOrderTypes,FieldRoutines
    INTEGER(INTG) :: maxNumberElementInterpolationParameters !<The maximum number of interpolation parameters in an element for a field variable. 
    INTEGER(INTG) :: maxNumberNodeInterpolationParameters !<The maximum number of interpolation parameters in an element for a field variable. 
    INTEGER(INTG) :: numberOfDofs !<Number of local degress of freedom for this field variable (excluding ghosted dofs). Old CMISS name NYNR(0,0,nc,nr,nx).
    INTEGER(INTG) :: totalNumberOfDofs !<Number of local degrees of freedom for this field variable (including ghosted dofs). Old CMISS name NYNR(0,0,nc,nr,nx).
    INTEGER(INTG) :: numberOfGlobalDofs !<Number of global degrees of freedom for this field variable. Old CMISS name NYNR(0,0,nc,nr,nx).
    TYPE(DomainMappingType), POINTER :: domainMapping !<Domain mapping for this variable. 
    INTEGER(INTG) :: numberOfComponents !<The number of components in the field variable.
    TYPE(FieldVariableComponentType), ALLOCATABLE :: components(:) !<components(componentIdx). The array of field variable components.
    TYPE(FieldDOFToParamMapType) :: dofToParamMap !<The mappings for the field dofs to the field parameters
    TYPE(FieldParameterSetsType) :: parameterSets !<The parameter sets for the field variable
  END TYPE FieldVariableType
  
  !>A buffer type to allow for an array of pointers to a FieldVariableType.
  TYPE FieldVariablePtrType
    TYPE(FieldVariableType), POINTER :: ptr !<The pointer to the field variable. 
  END TYPE FieldVariablePtrType

  !>A type to store a list of field variables
  TYPE FieldVariablesListType
    INTEGER(INTG) :: numberOfVariables !<The number of variables is a list
    TYPE(FieldVariablePtrType), ALLOCATABLE :: variables(:) !<variables(variableIdx). variables(variableIdx)%ptr is a pointer to the variableIdx'th field variable in the list.
  END TYPE FieldVariablesListType

  !>A type to temporarily hold (cache) the user modifiable values which are used to create a field. 
  TYPE FieldCreateValuesCacheType
    LOGICAL :: labelLocked !<Is .TRUE. if the field label has been locked, .FALSE. if not.
    LOGICAL :: decompositionLocked !<Is .TRUE. if the field decomposition has been locked, .FALSE. if not.
    LOGICAL :: DataProjectionLocked !<Is .TRUE. if the field data projection has been locked, .FALSE. if not.
    LOGICAL :: dependentTypeLocked !<Is .TRUE. if the field dependent type has been locked, .FALSE. if not.
    LOGICAL :: numberOfVariablesLocked !<Is .TRUE. if the number of field variables has been locked, .FALSE. if not.
    LOGICAL :: geometricFieldLocked !<Is .TRUE. if the geometric field has been locked, .FALSE. if not.
    LOGICAL :: scalingTypeLocked !<Is .TRUE. if the scaling type has been locked, .FALSE. if not.        
    LOGICAL :: typeLocked !<Is .TRUE. if the field type has been locked, .FALSE. if not.        
    INTEGER(INTG), ALLOCATABLE :: variableTypes(:) !<variableTypes(variableIdx). The cache of the variable type for the given variableIdx of the field. \see FieldRoutines_VariableTypes
    LOGICAL :: variableTypesLocked !<Is .TRUE. if the variable types have been locked, .FALSE. if not.
    TYPE(VARYING_STRING), ALLOCATABLE :: variableLabels(:) !<variableLabels(variableTypeIdx). The variable label for the variableTypeIdx'th variable type of the field.
    LOGICAL, ALLOCATABLE :: variableLabelsLocked(:) !<variableLabelsLocked(variableTypeIdx). Is .TRUE. if the variable label for the variableTypeIdx'th variable type has been locked, .FALSE. if not.
    INTEGER(INTG), ALLOCATABLE :: dimension(:) !<dimension(variableTypeIdx). The cache of the variable dimension for the variableTypeIdx'th variable type of the field. \see FieldRoutines_DimensionTypes
    LOGICAL, ALLOCATABLE :: dimensionLocked(:) !<dimensionLocked(variableTypeIdx). Is .TRUE. if the dimension for the variableTypeIdx'th variable type has been locked, .FALSE. if not.
    INTEGER(INTG), ALLOCATABLE :: dataTypes(:) !<dataTypes(variableTypeIdx). The cache of the variable data type for the variableTypeIdx'th variable type of the field. \see FieldRoutines_DataTypes
    LOGICAL, ALLOCATABLE :: dataTypesLocked(:) !<dataTypesLocked(variableTypeIdx). Is .TRUE. if the data type for the variableTypeIdx'th variable type has been locked, .FALSE. if not.
    INTEGER(INTG), ALLOCATABLE :: dofOrderTypes(:) !<dofOrderTypes(variableTypeIdx). The cache of the variable dof order type for the variableTypeIdx'th variable type of the field. \see FieldRoutines_DataTypes
    LOGICAL, ALLOCATABLE :: dofOrderTypesLocked(:) !<dofOrderTypesLocked(variableTypeIdx). Is .TRUE. if the dof order type for the variableTypeIdx'th variable type has been locked, .FALSE. if not.
    INTEGER(INTG), ALLOCATABLE :: numberOfComponents(:) !<numberOfComponents(variableTypeIdx). The number of components in the field for the variableTypeIdx'th field variable type.
    LOGICAL, ALLOCATABLE :: numberOfComponentsLocked(:) !<numberOfComponentsLocked(variableTypeIdx). Is .TRUE. if the number of components has been locked for the variableTypeIdx'th variable type, .FALSE. if not.
    TYPE(VARYING_STRING), ALLOCATABLE :: componentLabels(:,:) !<componentLabels(componentIdx,variableTypeIdx). The cache of the component label for the given component and variable type of the field.
    LOGICAL, ALLOCATABLE :: componentLabelsLocked(:,:) !<componentLabelsLocked(componentIdx,variableTypeIdx). Is .TRUE. if the component label of the componentIdx'th component of the variableTypeIdx'th varible type has been locked, .FALSE. if not.
    INTEGER(INTG), ALLOCATABLE :: interpolationTypes(:,:) !<interpolationTypes(componentIdx,variableTypeIdx). The cache of the interpolation type for the given component and variable type of the field. \see FieldRoutines_InterpolationTypes
    LOGICAL, ALLOCATABLE :: interpolationTypesLocked(:,:) !<interpolationTypesLocked(componentIdx,variableTypeIdx). Is .TRUE. if the interpolation type of the componentIdx'th component of the variableTypeIdx'th varible type has been locked, .FALSE. if not.
    INTEGER(INTG), ALLOCATABLE :: meshComponentNumber(:,:) !<meshComponentNumber(componentIdx,varaible_type_idx). The cache of the mesh component number for the given component and variable type of the field.
    LOGICAL, ALLOCATABLE :: meshComponentNumberLocked(:,:) !<meshComponentNumberLocked(componentIdx,variableTypeIdx). Is .TRUE. if the mesh component number of the componentIdx'th component of the variableTypeIdx'th varible type has been locked, .FALSE. if not.
  END TYPE FieldCreateValuesCacheType

  !>Contains information for a field defined on a region. \see OpenCMISS::Iron::cmfe_FieldType
  TYPE FieldType
    INTEGER(INTG) :: globalNumber !<The global number of the field in the list of fields for a region.
    INTEGER(INTG) :: userNumber !<The user defined identifier for the field. The user number must be unique.
    TYPE(VARYING_STRING) :: label !<The label for the field
    LOGICAL :: fieldFinished !<Is .TRUE. if the field has finished being created, .FALSE. if not.
    TYPE(FieldsType), POINTER :: fields !<A pointer to the fields for this region.
    TYPE(RegionType), POINTER :: region !<A pointer to the region containing the field. If the field are in an interface rather than a region then this pointer will be NULL and the interface pointer should be used.
    TYPE(InterfaceType), POINTER :: interface!<A pointer to the interface containing the field. If the field are in a region rather than an interface then this pointer will be NULL and the interface pointer should be used.
    INTEGER(INTG) :: type !<The type of the field. NOTE: this should be a field variable attribute as you may have a, say, geometric field variable and a general field variable bundled together in the same field. \see FieldRoutines_FieldTypes
    INTEGER(INTG) :: dependentType !<The dependent type of the field. \see FieldRoutines_DependentTypes
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition of the mesh for which the field is defined on.
    INTEGER(INTG) :: numberOfVariables !<The number of variable types in the field. Old CMISS name NCT(nr,nx)
    TYPE(FieldVariablePtrType), ALLOCATABLE :: variableTypeMap(:) !<variableTypeMap(variableIdx). The map from the available field variable types to the field variable types that are defined for the field. variableIdx varies from 1 to FieldRoutines::FIELD_NUMBER_OF_VARIABLE_TYPES. If the particular field variable type has not been defined on the field then the variableTypeMap will be NULL. \see FieldRoutines_VariableTypes
    TYPE(FieldVariableType), ALLOCATABLE :: variables(:) !<variables(variableIdx). The array of field variables. 
    TYPE(FieldScalingsType) :: scalings !<The scaling parameters for the field
    TYPE(FieldType), POINTER :: geometricField !<A pointer to the geometric field that this field uses. If the field itself is a geometric field then this will be a pointer back to itself.
    TYPE(FieldGeometricParametersType), POINTER :: geometricFieldParameters !<If the field is a geometric field the pointer to the geometric parameters (lines, areas, volumes etc.). If the field is not a geometric field the pointer is NULL.
    TYPE(FieldCreateValuesCacheType), POINTER :: createValuesCache !<The create values cache for the field.
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection that this field uses.
  END TYPE FieldType

  !>A buffer type to allow for an array of pointers to a FieldType.
  TYPE FieldPtrType
    TYPE(FieldType), POINTER :: ptr !<The pointer to the field.  
  END TYPE FieldPtrType

  !>Contains information on the fields defined on a region.
  TYPE FieldsType
    TYPE(RegionType), POINTER :: region !<A pointer to the region containing the fields. If the fields are in an interface rather than a region then this pointer will be NULL and the interface pointer should be used.
    TYPE(InterfaceType), POINTER :: INTERFACE !<A pointer to the interface containing the fields. If the fields are in a region rather than an interface then this pointer will be NULL and the interface pointer should be used.
    INTEGER(INTG) :: numberOfFields !<The number of fields defined on the region.
    TYPE(FieldPtrType), POINTER :: fields(:) !<fields(fieldsIdx). The array of pointers to the fields.
  END TYPE FieldsType

  PUBLIC FieldPhysicalPointType,FieldPhysicalPointPtrType

  PUBLIC FieldInterpolatedPointMetricsType,FieldInterpolatedPointMetricsPtrType

  PUBLIC FieldInterpolatedPointType,FieldInterpolatedPointPtrType

  PUBLIC FieldInterpolationParametersType,FieldInterpolationParametersPtrType

  PUBLIC FieldGeometricParametersType

  PUBLIC FieldScalingType,FieldScalingsType

  PUBLIC FieldDOFToParamMapType

  PUBLIC FieldNodeParamToDOFMapDerivativeType,FieldNodeParamToDOFMapNodeType,FieldNodeParamToDOFMapType

  PUBLIC FieldElementParamToDOFMapType

  PUBLIC FieldGridPointParamToDOFMapType

  PUBLIC FieldGaussPointParamToDOFMapType

  PUBLIC FieldDataPointParamToDOFMapType

  PUBLIC FieldParamToDOFMapType

  PUBLIC FieldVariableComponentType

  PUBLIC FieldParameterSetType,FieldParameterSetPtrType,FieldParameterSetsType

  PUBLIC FieldVariableType,FieldVariablePtrType

  PUBLIC FieldVariablesListType

  PUBLIC FieldCreateValuesCacheType

  PUBLIC FieldType,FieldPtrType,FieldsType

  !
  !================================================================================================================================
  !
  ! Equations matrices types
  !
  
  !>Contains information for an element matrix.
  TYPE ElementMatrixType
    INTEGER(INTG) :: equationsMatrixNumber !<The equations matrix number that this element matrix belongs to.
    INTEGER(INTG) :: structureType !<The structure type of the element matrix. \see EquationsMatricesRoutines_EquationsMatrixStructureTypes,EquationsMatricesRoutines
    INTEGER(INTG) :: numberOfRows !<The current number of rows in the element matrix.
    INTEGER(INTG) :: numberOfColumns !<The current number of columns in the element matrix.
    INTEGER(INTG) :: maxNumberOfRows !<The maximum (allocated) number of rows in the element matrix.
    INTEGER(INTG) :: maxNumberOfColumns !<The maximu (allocated) number of columns in the element matrix.
    INTEGER(INTG), ALLOCATABLE :: rowDOFS(:) !<rowDOFS(i). The equations row that the i'th row of the element matrix belongs to.
    INTEGER(INTG), ALLOCATABLE :: columnDOFS(:) !<columnDOFS(j). The equations column that the j'th column of the element matrix bleongs to.
    REAL(DP), ALLOCATABLE :: matrix(:,:) !<matrix(i,j). The vlaue of the i'th row and the j'th column of the element matrix.
  END TYPE ElementMatrixType

  !>Contains information for an element vector.
  TYPE ElementVectorType
    INTEGER(INTG) :: numberOfRows !<The current number of rows in the element vector
    INTEGER(INTG) :: maxNumberOfRows !<The maximum (allocated) number of rows in the element vecotr
    INTEGER(INTG), ALLOCATABLE :: rowDOFS(:) !<rowDOFS(i). The equations row that the i'th row of the element vector belongs to
    REAL(DP), ALLOCATABLE :: vector(:) !<vector(i). The value of the i'th row of the element vector
  END TYPE ElementVectorType

  !>Contains information for an nodal matrix.
  TYPE NodalMatrixType
    INTEGER(INTG) :: equationsMatrixNumber !<The equations matrix number that this nodal matrix belongs to.
    INTEGER(INTG) :: structureType !<The structure type of the nodal matrix. \see EquationsMatricesRoutines_EquationsMatrixStructureTypes,EquationsMatricesRoutines
    INTEGER(INTG) :: numberOfRows !<The current number of rows in the nodal matrix.
    INTEGER(INTG) :: numberOfColumns !<The current number of columns in the nodal matrix.
    INTEGER(INTG) :: maxNumberOfRows !<The maximum (allocated) number of rows in the nodal matrix.
    INTEGER(INTG) :: maxNumberOfColumns !<The maximum (allocated) number of columns in the nodal matrix.
    INTEGER(INTG), ALLOCATABLE :: rowDOFS(:) !<rowDOFS(i). The equations row that the i'th row of the nodal matrix belongs to.
    INTEGER(INTG), ALLOCATABLE :: columnDOFS(:) !<columnDOFS(j). The equations column that the j'th column of the nodal matrix bleongs to.
    REAL(DP), ALLOCATABLE :: matrix(:,:) !<matrix(i,j). The vlaue of the i'th row and the j'th column of the nodal matrix.
  END TYPE NodalMatrixType

  !>Contains information for an nodal vector.
  TYPE NodalVectorType
    INTEGER(INTG) :: numberOfRows !<The current number of rows in the nodal vector
    INTEGER(INTG) :: maxNumberOfRows !<The maximum (allocated) number of rows in the nodal vecotr
    INTEGER(INTG), ALLOCATABLE :: rowDOFS(:) !<rowDOFS(i). The equations row that the i'th row of the nodal vector belongs to
    REAL(DP), ALLOCATABLE :: vector(:) !<vector(i). The value of the i'th row of the nodal vector
  END TYPE NodalVectorType

  !>Contains information about an equations matrix.
  TYPE EquationsMatrixType
    INTEGER(INTG) :: matrixNumber !<The number of the equations matrix
    INTEGER(INTG) :: matrixMatricesType !<The type of matrices for the matrix \see EquationsMatricesRoutines_MatrixTypes,EquationsMatricesRoutines
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices !<A pointer to the dynamic equations matrices for the dynamic equation matrix.
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices !<A pointer to the linear equations matrices for the linear equation matrix.
    REAL(DP) :: matrixCoefficient !<The coefficient applied to the equations matrix
    INTEGER(INTG) :: storageType !<The storage (sparsity) type for this matrix 
    INTEGER(INTG) :: structureType !<The structure (sparsity) type for this matrix \see EquationsMatricesRoutines_EquationsMatrixStructureTypes,EquationsMatricesRoutines
    LOGICAL :: lumped !<Is .TRUE. if the equations matrix is lumped
    LOGICAL :: symmetric !<Is .TRUE. if the equations matrix is symmetric
    INTEGER(INTG) :: numberOfColumns !<The number of columns in this equations matrix
    LOGICAL :: updateMatrix !<Is .TRUE. if this equations matrix is to be updated
    LOGICAL :: firstAssembly !<Is .TRUE. if this equations matrix has not been assembled
    TYPE(DistributedMatrixType), POINTER :: matrix !<A pointer to the distributed equations matrix data
    TYPE(ElementMatrixType) :: elementMatrix !<The element matrix for this equations matrix
    TYPE(NodalMatrixType) :: nodalMatrix !<The nodal matrix for this equations matrix
    TYPE(DistributedVectorType), POINTER :: tempVector !<Temporary vector used for assembly. 
  END TYPE EquationsMatrixType

  !>A buffer type to allow for an array of pointers to a EquationsMatrixType \see Types::EquationsMatrixType.
  TYPE EquationsMatrixPtrType
    TYPE(EquationsMatrixType), POINTER :: ptr !<A pointer to the equations matrix.
  END TYPE EquationsMatrixPtrType
 
  !>Contains information on the Jacobian matrix for nonlinear problems
  TYPE JacobianMatrixType
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector !<A pointer back to the residual vector for this Jacobian
    INTEGER(INTG) :: jacobianNumber !<The equations Jacobian matrix number
    REAL(DP) :: jacobianCoefficient !<The coefficient applied to the Jacobian matrix
    INTEGER(INTG) :: storageType !<The storage (sparsity) type for this matrix
    INTEGER(INTG) :: structureType !<The structure (sparsity) type for this matrix
    LOGICAL :: symmetric !<Is .TRUE. if the Jacobian matrix is symmetric
    INTEGER(INTG) :: numberOfColumns !<The number of columns in this global matrix
    LOGICAL :: updateJacobian !<Is .TRUE. if this Jacobian matrix is to be updated
    TYPE(DistributedMatrixType), POINTER :: jacobian !<A pointer to the distributed jacobian matrix data
    LOGICAL :: firstAssembly !<Is .TRUE. if this Jacobian matrix has not been assembled
    TYPE(ElementMatrixType) :: elementJacobian !<The element matrix for this Jacobian matrix. This is not used if the Jacobian is not supplied.
    TYPE(NodalMatrixType) :: nodalJacobian !<The nodal matrix for this Jacobian matrix. This is not used if the Jacobian is not supplied.
    INTEGER(INTG) :: jacobianCalculationType !<The calculation type (analytic of finite difference) of the Jacobian. \see EquationsMatricesRoutines_JacobianCalculationTypes
    REAL(DP) :: jacobianFiniteDifferenceStepSize !<The finite difference step size used for calculating the Jacobian.
  END TYPE JacobianMatrixType

  !>A buffer type to allow for an array of pointers to a JacobianMatrixType \see Types::JacobianMatrixType.
  TYPE JacobianMatrixPtrType
    TYPE(JacobianMatrixType), POINTER :: ptr
  END TYPE JacobianMatrixPtrType

  !>Contains information on the Hessian matrix for optimisation problems
  TYPE HessianMatrixType
    INTEGER(INTG) :: hessianNumber !<The equations Hessian matrix number
    TYPE(EquationsMatricesOptimisationType), POINTER :: optimisationMatrices !<A pointer back to the optimisation matrices for this Hessian
    INTEGER(INTG) :: storageType !<The storage (sparsity) type for this matrix
    INTEGER(INTG) :: structureType!<The structure (sparsity) type for this matrix
    LOGICAL :: symmetric !<Is .TRUE. if the Hessian matrix is symmetric
    INTEGER(INTG) :: numberOfColumns !<The number of columns in this global matrix
    LOGICAL :: updateHessian !<Is .TRUE. if this Hessian matrix is to be updated
    TYPE(DistributedMatrixType), POINTER :: hessian !<A pointer to the distributed Hessian matrix data
    LOGICAL :: firstAssembly !<Is .TRUE. if this Hessian matrix has not been assembled
    TYPE(ElementMatrixType) :: elementHessian !<The element matrix for this Hessian matrix. This is not used if the Hessian is not supplied.
    INTEGER(INTG) :: hessianCalculationType !<The calculation type (analytic of finite difference) of the Hessian.
  END TYPE HessianMatrixType

  !>A buffer type to allow for an array of pointers to a HessianMatrixType \see Types::HessianMatrixType
  TYPE HessianMatrixPtrType
    TYPE(HessianMatrixType), POINTER :: ptr
  END TYPE HessianMatrixPtrType

  !>Contains information on functions for scalar equations
  TYPE EquationsMatricesFunctionType
    TYPE(EquationsMatricesScalarType), POINTER :: scalarMatrices !<A pointer back to the scalar matrices
  END TYPE EquationsMatricesFunctionType
  
  !>Contains information on norms for scalar equations
  TYPE EquationsMatricesNormType
    TYPE(EquationsMatricesScalarType), POINTER :: scalarMatrices !<A pointer back to the scalar matrices
  END TYPE EquationsMatricesNormType
  
  !>Contains information on dot products for scalar equations
  TYPE EquationsMatricesDotProductType
    TYPE(EquationsMatricesScalarType), POINTER :: scalarMatrices !<A pointer back to the scalar matrices
  END TYPE EquationsMatricesDotProductType
  
  !>Contains information on quadratic matrix vector forms for scalar equations
  TYPE EquationsMatricesQuadraticType
    TYPE(EquationsMatricesScalarType), POINTER :: scalarMatrices !<A pointer back to the scalar matrices
  END TYPE EquationsMatricesQuadraticType
  
  !>Contains information of the dynamic matrices for equations matrices
  TYPE EquationsMatricesDynamicType
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer back to the vector equations matrices.
    INTEGER(INTG) :: numberOfDynamicMatrices !<The number of dynamic equations matrices defined for the equations set.
    TYPE(EquationsMatrixPtrType), ALLOCATABLE :: matrices(:) !<matrices(matrixIdx)%ptr contains the information on the matrixIdx'th dynamic equations matrix.
    TYPE(DistributedVectorType), POINTER :: tempVector !<Temporary vector used for assembly. 
  END TYPE EquationsMatricesDynamicType

  !>Contains information of the linear matrices for equations matrices
  TYPE EquationsMatricesLinearType
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer back to the vector equations matrices.
    INTEGER(INTG) :: numberOfLinearMatrices !<The number of linear equations matrices defined for the equations set.
    TYPE(EquationsMatrixPtrType), ALLOCATABLE :: matrices(:) !<matrices(matrixIdx)%ptr contains the information on the matrixIdx'th linear equations matrix.
    TYPE(DistributedVectorType), POINTER :: tempVector !<Temporary vector used for assembly. 
  END TYPE EquationsMatricesLinearType

  !>Contains information on the residual vectors for nonlinear matrices
  TYPE EquationsMatricesResidualType
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices !<A pointer back to the nonlinear matrices.
    INTEGER(INTG) :: residualNumber !<The residual number for this residual
    REAL(DP) :: residualCoefficient !<The coefficient applied to the residual coefficient
    LOGICAL :: updateResidual !<Is .TRUE. if the equations residual vector is to be updated
    LOGICAL :: firstAssembly !<Is .TRUE. if this residual vector has not been assembled
    TYPE(DistributedVectorType), POINTER :: residual !<A pointer to the distributed residual vector for nonlinear equations
    TYPE(DistributedVectorType), POINTER :: previousResidual !<A point to the previous distributed residual vector.
    TYPE(DistributedVectorType), POINTER :: previous2Residual !<A point to the second previous distributed residual vector.
    TYPE(DistributedVectorType), POINTER :: previous3Residual !<A point to the third previous distributed residual vector.
    TYPE(ElementVectorType) :: elementResidual !<The element residual information for nonlinear equations. Old CMISS name RE1
    TYPE(NodalVectorType) :: nodalResidual !<The nodal residual information for nonlinear equations.
    INTEGER(INTG) :: nodalResidualCalculated !<The number of the nodal the residual is calculated for, or zero if it isn't calculated
    INTEGER(INTG) :: elementResidualCalculated !<The number of the element the residual is calculated for, or zero if it isn't calculated
    INTEGER(INTG) :: numberOfJacobians !<The number of Jacobian matrices for this residual.
    TYPE(JacobianMatrixPtrType), ALLOCATABLE :: jacobians(:) !<jacobians(matrixIdx)%ptr is a pointer to the matrixIdx'th Jacobian matrix for this residual     
  END TYPE EquationsMatricesResidualType

  !>A buffer type to allow for arrays of pointers to a EquationsMatricesResidualType \see Types::EquationsMatricesResidualType
  TYPE EquationsMatricesResidualPtrType
    TYPE(EquationsMatricesResidualType), POINTER :: ptr !<A pointer to the EquationsMatricesResidualType
  END TYPE EquationsMatricesResidualPtrType

  !>Contains information of the nolinear matrices and vectors for equations matrices
  TYPE EquationsMatricesNonlinearType
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer back to the vector equations matrices.
    INTEGER(INTG) :: numberOfResiduals !<The number of residuals for this nonlinear matrices
    TYPE(EquationsMatricesResidualPtrType), ALLOCATABLE :: residuals(:) !<residuals(residualIdx). residuals(residualIdx)%ptr is a pointer to the residualIdx'th residual for this nonlinear matrices
    TYPE(DistributedVectorType), POINTER :: tempVector !<Temporary vector used for assembly. 
  END TYPE EquationsMatricesNonlinearType

  !>Contains information of the optimisation matrices and vectors for equations matrices
  TYPE EquationsMatricesOptimisationType
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer back to the equations matrices.
    REAL(DP) :: objective !<The value of the objective
    INTEGER(INTG) :: numberOfHessians !<The number of Hessian matrices for the equations set.
    TYPE(HessianMatrixPtrType), ALLOCATABLE :: hessians(:) !<hessians(matrixIdx)%ptr is a pointer to the matrixIdx'th Hessian matrix for optimisation equations
    LOGICAL :: updateGradient !<Is .TRUE. if the equations gradient vector is to be updated
    TYPE(DistributedVectorType), POINTER :: gradient !<A pointer to the distributed gradient vector for optimisation equations
    TYPE(ElementVectorType) :: elementGradient !<The element gradient information for optimisation equations.
    LOGICAL :: updateConstraints !<Is .TRUE. if the equations constraints vector is to be updated
    TYPE(DistributedVectorType), POINTER :: constraints !<A pointer to the distributed constraints vector for optimisation equations
    TYPE(ElementVectorType) :: elementConstraints !<The element constraints information for optimisation equations.
    LOGICAL :: updateBounds !<Is .TRUE. if the equations bounds vectors are to be updated
    TYPE(DistributedVectorType), POINTER :: lowerBounds !<A pointer to the distributed lower bounds vector for optimisation equations
    TYPE(DistributedVectorType), POINTER :: upperBounds !<A pointer to the distributed upper bounds vector for optimisation equations
    TYPE(ElementVectorType) :: elementLowerBounds !<The element lower bounds information for optimisation equations.
    TYPE(ElementVectorType) :: elementUpperBounds !<The element upper bounds information for optimisation equations.
    LOGICAL :: updateResidual !<Is .TRUE. if the equations residual vector is to be updated
    LOGICAL :: firstAssembly !<Is .TRUE. if this residual vector has not been assembled
    TYPE(DistributedVectorType), POINTER :: residual !<A pointer to the distributed residual vector for optimisation equations
    TYPE(ElementVectorType) :: elementResidual !<The element residual information for optimisation equations. Old CMISS name RE1
  END TYPE EquationsMatricesOptimisationType

  !>Contains information of the RHS vector for equations matrices
  TYPE EquationsMatricesRHSType
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer back to the vector equations matrices.
    REAL(DP) :: rhsCoefficient !<The coefficient applied to the RHS vector
    LOGICAL :: updateVector !<Is .TRUE. if the equations rhs vector is to be updated
    LOGICAL :: firstAssembly !<Is .TRUE. if this rhs vector has not been assembled
    TYPE(DistributedVectorType), POINTER :: vector !<A pointer to the distributed global rhs vector data \todo rename this RHS_VECTOR
    TYPE(DistributedVectorType), POINTER :: previousRHSVector !<A pointer to the previous distributed RHS vector data
    TYPE(DistributedVectorType), POINTER :: previous2RHSVector !<A pointer to the second previous distributed RHS vector data
    TYPE(DistributedVectorType), POINTER :: previous3RHSVector !<A pointer to the third previous distributed RHS vector data
    TYPE(ElementVectorType) :: elementVector !<The element rhs information
    TYPE(NodalVectorType) :: nodalVector !<The nodal rhs information
  END TYPE EquationsMatricesRHSType
  
  !>Contains information of the source vector for equations matrices
  TYPE EquationsMatricesSourceType
    TYPE(EquationsMatricesSourcesType), POINTER :: sources !<A pointer back to the sources.
    INTEGER(INTG) :: sourceNumber !<The source number for the sources in the vector equations matrices
    REAL(DP) :: sourceCoefficient !<The coefficient applied to the source vector
    LOGICAL :: updateVector !<Is .TRUE. if the equations rhs vector is to be updated
    LOGICAL :: firstAssembly !<Is .TRUE. if this source vector has not been assembled
    TYPE(DistributedVectorType), POINTER :: vector !<A pointer to the distributed source vector data \todo rename this SOURCE_VECTOR
    TYPE(DistributedVectorType), POINTER :: previousSourceVector !<A pointer to the previous distributed source vector data
    TYPE(DistributedVectorType), POINTER :: previous2SourceVector !<A pointer to the second previous distributed source vector data
    TYPE(DistributedVectorType), POINTER :: previous3SourceVector !<A pointer to the third previous distributed source vector data
    TYPE(ElementVectorType) :: elementVector !<The element source information
    TYPE(NodalVectorType) :: nodalVector !<The nodal source information
  END TYPE EquationsMatricesSourceType

  !>A buffer type to allow for arrays of pointer to EquationsMatricesSourceType
  TYPE EquationsMatricesSourcePtrType
    TYPE(EquationsMatricesSourceType), POINTER :: ptr !<A pointer to the EquationsMatricesSourceType
  END TYPE EquationsMatricesSourcePtrType
  
  !>Contains information on the source vectors for the vector equations matrices.
  TYPE EquationsMatricesSourcesType
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer back to the vector equations matrices
    INTEGER(INTG) :: numberOfSources !<The number of sources in the vector equations. 
    TYPE(EquationsMatricesSourcePtrType), ALLOCATABLE :: sources(:) !<sources(sourceIdx). sources(sourceIdx)%ptr is a pointer to the sourceIdx'th source.
    TYPE(DistributedVectorType), POINTER :: tempVector !<Temporary vector used for assembly. 
  END TYPE EquationsMatricesSourcesType
  
  !>Contains information on the scalar equations matrices, vectors and scalars
  TYPE EquationsMatricesScalarType
    TYPE(EquationsScalarType), POINTER :: scalarEquations !<A pointer back to the scalar equations
    LOGICAL :: scalarMatricesFinished !<Is .TRUE. if the scalar equations matrices have finished being created, .FALSE. if not.
    TYPE(EquationsMappingScalarType), POINTER :: scalarMapping !<A pointer to the scalar mapping for the scalar equations matrices.
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping for the equations matrices
    TYPE(EquationsMatricesFunctionType), POINTER :: functions !<A pointer to the functions information for the scalar equations matrices
    TYPE(EquationsMatricesNormType), POINTER :: normMatrices !<A pointer to the norm matrices and vectors for the scalar equations matrices
    TYPE(EquationsMatricesDotProductType), POINTER :: dotProductMatrices !<A pointer to the dot product matrices and vectors for the scalar equations matrices
    TYPE(EquationsMatricesQuadraticType), POINTER :: quadraticMatrices !<A pointer to the quadratic matrices and vectors for the scalar equations matrices
  END TYPE EquationsMatricesScalarType
  
  !>Contains information on the vector equations matrices and vectors
  TYPE EquationsMatricesVectorType
    TYPE(EquationsVectorType), POINTER :: vectorEquations !<A pointer back to the vector equations
    LOGICAL :: vectorMatricesFinished !<Is .TRUE. if the vector equations matrices have finished being created, .FALSE. if not.
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the vector mapping for the vector equations matrices.
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping for the equations matrices
    INTEGER(INTG) :: numberOfRows !<The number of local rows (excluding ghost rows) in the distributed vector equations matrices and vectors
    INTEGER(INTG) :: totalNumberOfRows !<The number of local rows (including ghost rows) in the distributed vector equations matrices and vectors
    INTEGER(INTG) :: numberOfGlobalRows !<The number of global rows in the distributed vector equations matrices and vectors
    !Equations matrices components
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices!<A pointer to the dynamic matrices information for the vector equations matrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices !<A pointer to the linear matrices information for the vector equations matrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices !<A pointer to the nonlinear matrices and vectors information for the vector equations matrices
    TYPE(EquationsMatricesOptimisationType), POINTER :: optimisationMatrices !<A pointer to the optimisation matrices and vectors information for the equations matrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector !<A pointer to the RHS vector information for the vector equations matrices
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors !<A pointer to the source vectors information for the vector equations matrices
  END TYPE EquationsMatricesVectorType

  PUBLIC ElementMatrixType,ElementVectorType

  PUBLIC NodalMatrixType,NodalVectorType

  PUBLIC EquationsMatrixType,EquationsMatrixPtrType

  PUBLIC JacobianMatrixType,JacobianMatrixPtrType

  PUBLIC HessianMatrixType,HessianMatrixPtrType

  PUBLIC EquationsMatricesFunctionType

  PUBLIC EquationsMatricesNormType

  PUBLIC EquationsMatricesDotProductType

  PUBLIC EquationsMatricesQuadraticType

  PUBLIC EquationsMatricesDynamicType,EquationsMatricesLinearType,EquationsMatricesNonlinearType

  PUBLIC EquationsMatricesOptimisationType

  PUBLIC EquationsMatricesResidualType,EquationsMatricesResidualPtrType
  
  PUBLIC EquationsMatricesRHSType

  PUBLIC EquationsMatricesSourceType,EquationsMatricesSourcePtrType

  PUBLIC EquationsMatricesSourcesType

  PUBLIC EquationsMatricesScalarType

  PUBLIC EquationsMatricesVectorType
  
  !
  !================================================================================================================================
  !
  ! Equations mapping types
  !
  
   !>Contains the information about the mapping of a variable DOF to an equations matrix column
  TYPE VarToEquationsColumnMapType
    INTEGER(INTG), ALLOCATABLE :: columnDOF(:) !<columnDOF(dofIdx). The equations column number for this equations matrix that the dofIdx'th variable DOF is mapped to.  
  END TYPE VarToEquationsColumnMapType

  !>Contains the mapping for a dependent variable type to the equations matrices
  TYPE VarToEquationsMatricesMapType
    INTEGER(INTG) :: variableIndex !<The variable index for this variable to equations matrices map
    INTEGER(INTG) :: variableType !<The variable type for this variable to equations matrices map
    TYPE(FieldVariableType), POINTER :: variable !<A pointer to the field variable for this variable to equations matrices map
    INTEGER(INTG) :: numberOfEquationsMatrices !<The number of equations matrices (linear or dynamic) this variable type is mapped to. If the number is -1 the variable is mapped to the RHS vector. If the number is zero then this variable type is not involved in the equations set and the rest of the type is not allocated.
    INTEGER(INTG), ALLOCATABLE :: equationsMatrixNumbers(:) !<equationsMatrixNumbers(i). The equations matrix number for the i'th matrix that this variable type is mapped to.
    TYPE(varToEquationsColumnMapType), ALLOCATABLE :: dofToColumnsMaps(:) !<dofToColumnsMaps(i). The variable dof to equations columns for the i'th equations matrix.
  END TYPE VarToEquationsMatricesMapType

  !>A buffer type of pointers to VarToEquationsMatricesMapType \see Types::VarToEquationsMatricesMapType
  TYPE VarToEquationsMatricesMapPtrType
    TYPE(VarToEquationsMatricesMapType), POINTER :: ptr !<A pointer to a VarToEquationsMatricesMapType
  END TYPE VarToEquationsMatricesMapPtrType

  !>Contains information for mapping an equations matrix to a field variable.
  TYPE EquationsMatrixToVarMapType
    INTEGER(INTG) :: matrixNumber !<The equations matrix number
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix
    INTEGER(INTG) :: variableIndex !<The variable index in the list of linear variables for this equation matrix.
    INTEGER(INTG) :: variableType !<The dependent variable type mapped to this equations matrix
    TYPE(FieldVariableType), POINTER :: variable !<A pointer to the field variable that is mapped to this equations matrix
    INTEGER(INTG) :: numberOfColumns !<The number of columns in this equations matrix.
    REAL(DP) :: matrixCoefficient !<The multiplicative coefficent for the matrix in the equation set
    INTEGER(INTG), POINTER :: columnToDOFMap(:) !<columnToDOFMap(columnIdx). The variable DOF that the columnIdx'th column of this equations matrix is mapped to.
    TYPE(DomainMappingType), POINTER :: columnDOFsMapping !<A pointer to the column dofs domain mapping for the matrix variable
  END TYPE EquationsMatrixToVarMapType

  !>A buffer type of pointers to a EquationsMatrixToVarMapType \see Types::EquationsMatrixToVarMapType
  TYPE EquationsMatrixToVarMapPtrType
    TYPE(EquationsMatrixToVarMapType), POINTER :: ptr !<A pointer to a EquationsMatrixToVarMapType
  END TYPE EquationsMatrixToVarMapPtrType
  
  !>Contains information for function mapping in the scalar equations mapping
  TYPE EquationsMappingFunctionType
    INTEGER(INTG) :: functionNumber !<The number of the function in the functions scalar mapping
 END TYPE EquationsMappingFunctionType
    
  !>Contains information for mapping field variables to functions in the scalar equations mapping
  TYPE EquationsMappingFunctionsType
    TYPE(EquationsMappingScalarType), POINTER :: scalarMapping !<A pointer to the equations scalar mapping
    INTEGER(INTG) :: numberOfFunctions !<The number of functions in the mapping
    TYPE(EquationsMappingFunctionType), ALLOCATABLE :: functions(:) !<functions(functionIdx). Information on the functionIdx'th function mapping.
  END TYPE EquationsMappingFunctionsType
    
  !>Contains information for mapping field variables to a norm i.e., ||x|| in the scalar equations mapping
  TYPE EquationsMappingNormType
    INTEGER(INTG) :: normNumber !<The number of the norm in the norms scalar mapping
    TYPE(FieldVariableType), POINTER :: normVariable !<A pointer to the field variable for this norm mapping
    REAL(DP) :: normCoefficient !<The multiplicative coefficient applied to the norm.
  END TYPE EquationsMappingNormType
    
  !>Contains information for mapping field variables to norms i.e., ||x|| in the scalar equations mapping
  TYPE EquationsMappingNormsType
    TYPE(EquationsMappingScalarType), POINTER :: scalarMapping !<A pointer to the equations scalar mapping
    INTEGER(INTG) :: numberOfNorms !<The number of norms in the scalar mapping
    TYPE(EquationsMappingNormType), ALLOCATABLE :: norms(:) !<norms(normIdx). Information on the normIdx'th norm mapping.
  END TYPE EquationsMappingNormsType
    
  !>Contains information for mapping field variables to a dot product i.e., x^T.y in the scalar equations mapping
  TYPE EquationsMappingDotProductType
    INTEGER(INTG) :: dotProductNumber !<The number of dot product in the dot products scalar mapping
    TYPE(FieldVariablePtrType) :: dotProductVariables(2) !<dotProductVaraibles(variableIdx). The variableIdx'th field variable in the dot product. For x^T.y the first variable is x and the second variable is y.
    REAL(DP) :: dotProductCoefficient !<The multiplicative coefficient applied to the dot product.
  END TYPE EquationsMappingDotProductType
    
   !>Contains information for mapping field variables to dot products i.e., x^T.y in the equations set of the mapping
  TYPE EquationsMappingDotProductsType
    TYPE(EquationsMappingScalarType), POINTER :: scalarMapping !<A pointer to the equations scalar mapping
    INTEGER(INTG) :: numberOfDotProducts !<The number of dot products in the scalar mapping
    TYPE(EquationsMappingDotProductType), ALLOCATABLE :: dotProducts(:) !<dotProducts(dotProductIdx). Information on the dotProductIdx'th dot product mapping.
  END TYPE EquationsMappingDotProductsType
    
  !>Contains information for mapping field variables to a quadratic form i.e., x^T.A.y in the scalar equations mapping
  TYPE EquationsMappingQuadraticType
    INTEGER(INTG) :: quadraticNumber !<The number of dot product in the dot products scalar mapping    
    TYPE(FieldVariablePtrType) :: quadraticVariables(2) !<dotProductVaraibles(variableIdx). The variableIdx'th field variable in the dot product. For x^T.A.y the first variable is x and the second variable is y.
    REAL(DP) :: quadraticCoefficient !<The multiplicative coefficient applied to the quadratic form.
  END TYPE EquationsMappingQuadraticType
    
  !>Contains information for mapping field variables to quadratic forms i.e., x^T.A.y in the scalar equations mapping
  TYPE EquationsMappingQuadraticsType
    TYPE(EquationsMappingScalarType), POINTER :: scalarMapping !<A pointer to the equations scalar mapping
    INTEGER(INTG) :: numberOfQuadratics !<The number of quadratics in the scalar mapping
    TYPE(EquationsMappingQuadraticType), ALLOCATABLE :: quadratics(:) !<quadratics(quadraticIdx). Information on the quadraticIdx'th quadratic mapping.
  END TYPE EquationsMappingQuadraticsType
    
  !>Contains information for mapping field variables to the dynamic matrices in the vector equations mapping
  TYPE EquationsMappingDynamicType
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping
    INTEGER(INTG) :: numberOfDynamicMatrices !<The number of dynamic equations matrices in this mapping
    INTEGER(INTG) :: stiffnessMatrixNumber !<The matrix number of the dynamic stiffness matrix. 0 if there is no dynamic stiffness matrix
    INTEGER(INTG) :: dampingMatrixNumber !<The matrix number of the dynamic damping matrix. 0 if there is no dynamic damping matrix
    INTEGER(INTG) :: massMatrixNumber !<The matrix number of the dynamic mass matrix. 0 if there is no dynamic mass matrix
    INTEGER(INTG) :: dynamicVariableType !<The variable type involved in the equations matrix mapping.
    TYPE(FieldVariableType), POINTER :: dynamicVariable !<A pointer to the variable that is mapped to the dynamic matrices.
    TYPE(VarToEquationsMatricesMapType), POINTER :: varToEquationsMatricesMap !<The equations matrices mapping for the dynamic variable.
    TYPE(EquationsMatrixToVarMapPtrType), ALLOCATABLE :: equationsMatrixToVarMaps(:) !<equationsMatrixToVarMaps(matrixIdx). The mappings for the matrixIdx'th equations matrix.
    INTEGER(INTG), POINTER :: equationsRowToVariableDOFMaps(:) !<equationsRowToVariableDOFMaps(rowIdx). The row mappings for the rowIdx'th row of the equations matrices to the dynamic variable.
  END TYPE EquationsMappingDynamicType

  !>Contains information for mapping field variables to the linear matrices in the vector equations mapping
  TYPE EquationsMappingLinearType
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping
    INTEGER(INTG) :: numberOfLinearMatrices !<The number of linear equations matrices in this mapping
    INTEGER(INTG) :: numberOfLinearVariables !<The number of independent variables in the linear mapping
    INTEGER(INTG), ALLOCATABLE :: linearVariableTypes(:) !<linearVariableTypes(variableIdx). The variable type for the variableIdx'th independent linear variable.
    INTEGER(INTG), ALLOCATABLE :: linearVariableTypesMap(:) !<linearVariableTypesMap(variableTypeIdx). The variable index for the variableTypeIdx'th variable type. 
    TYPE(FieldVariablePtrType), ALLOCATABLE :: linearVariables(:) !<linearVariables(variableIdx). linearVariables(variableIdx)%ptr is the pointer to the linear field variable for the variableIdx'th independent linear variable.
    TYPE(VarToEquationsMatricesMapPtrType), ALLOCATABLE :: varToEquationsMatricesMaps(:) !<varToEquationsMatricesMaps(variableIdx). The equations matrices mapping for the variableIdx'th variable.
    TYPE(EquationsMatrixToVarMapPtrType), ALLOCATABLE :: equationsMatrixToVarMaps(:) !<equationsMatrixToVarMaps(matrixIdx). The mappings for the matrixIdx'th equations matrix.
    INTEGER(INTG), ALLOCATABLE :: equationsRowToVariableDOFMaps(:,:) !<equationsRowToVariableDOFMaps(rowIdx,variableIdx). The row mappings for the rowIdx'th row of the equations matrices to the variableIdx'th variable.
  END TYPE EquationsMappingLinearType

  !>Contains the mapping from the Jacobian back to the nonlinear residual variables.
  TYPE JacobianMatrixToVarMapType
    INTEGER(INTG) :: jacobianNumber !<The equations Jacobian matrix number
    INTEGER(INTG) :: variableType !<The dependent variable type mapped to this equations matrix
    TYPE(FieldVariableType), POINTER :: variable !<A pointer to the field variable that is mapped to this equations matrix
    TYPE(JacobianMatrixType), POINTER :: jacobian !<A pointer to the equations matrix for this variable
    INTEGER(INTG) :: numberOfColumns !<The number of columns in this equations matrix.
    REAL(DP) :: jacobianCoefficient !<The multiplicative coefficent for the matrix in the equation set
    INTEGER(INTG), POINTER :: equationsColumnToDOFVariableMap(:) !<equationsColumnToDOFVariableMap(columnIdx). The variable DOF that the columnIdx'th column of this equations matrix is mapped to.
    TYPE(DomainMappingType), POINTER :: columnDOFSMapping !<A pointer to the column dofs domain mapping for the matrix variable
  END TYPE JacobianMatrixToVarMapType

  !>A buffer type of pointers to an JacobianMatrixToVarMapPtrType \see Types::JacobianMatrixToVarMapPtrType
  TYPE JacobianMatrixToVarMapPtrType
    TYPE(JacobianMatrixToVarMapType), POINTER :: ptr !<A pointer to the JacobianMatrixToVarMapPtrType
  END TYPE JacobianMatrixToVarMapPtrType

  !>Contains the mapping for a dependent variable type to the nonlinear Jacobian matrix
  TYPE VarToJacobianMatrixMapType
    INTEGER(INTG) :: jacobianNumber !<The equations Jacobian matrix number
    INTEGER(INTG) :: variableType !<The variable type for this variable to equations matrices map
    TYPE(FieldVariableType), POINTER :: variable !<A pointer to the field variable for this variable to equations matrices map
    INTEGER(INTG), ALLOCATABLE :: dofToColumnsMap(:) !<dofToColumnsMap(dofIdx). The Jacobian column number for dofIdx'th variable dof
  END TYPE VarToJacobianMatrixMapType

  !>A buffer type of pointer to a VarToJacobianMatrixMapType \see Types::VarToJacobianMatrixMapType
  TYPE VarToJacobianMatrixMapPtrType
    TYPE(VarToJacobianMatrixMapType), POINTER :: ptr !<A pointer to the VarToJacobianMatrixMapType
  END TYPE VarToJacobianMatrixMapPtrType

  !>Contains information on the equations mapping for a residual 
  TYPE EquationsMappingResidualType
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping !<A pointer back to the nonlinear mapping for this residual
    INTEGER(INTG) :: residualNumber !<The residual number in the nonlinear mapping
    INTEGER(INTG) :: numberOfVariables !<The number of residual variables in this mapping
    INTEGER(INTG) :: numberOfJacobianMatrices !<The number of nonlinear Jacobian matrices in this mapping
    INTEGER(INTG), ALLOCATABLE :: variableTypes(:) !<variableTypes(variableIdx). The variable type for the variableIdx'th independent residual variable.
    INTEGER(INTG), ALLOCATABLE :: variableTypesMap(:) !<variableTypesMap(variableTypeIdx). The residual variable index for the variableTypeIdx'th variable type. 
    TYPE(FieldVariablePtrType), ALLOCATABLE :: variables(:) !<variables(variableIdx). variables(variableId)%ptr is a pointer to the variableIdx'th residual variable.
    TYPE(VarToJacobianMatrixMapPtrType), ALLOCATABLE :: varToJacobianMatrixMaps(:) !<varToJacobianMap(variableIdx). varToJacobianMap(variableIdx)%ptr is a pointer to the mapping from the residual variable to the Jacobain matrix for the variableIdx'th residual variable.
    TYPE(JacobianMatrixToVarMapPtrType), ALLOCATABLE :: jacobianMatrixToVarMaps(:) !<jacobianMatrixToVarMaps(jacobianIdx). <jacobianMatrixToVarMaps(jacobianIdx)%ptr is a pointer to the mapping from the Jacobian matrix to the residual variables for the jacobianIdx'th Jacobian.
    REAL(DP) :: residualCoefficient !<The multiplicative coefficient applied to the residual vector
    INTEGER(INTG), POINTER :: equationsRowToResidualDOFMap(:) !<equationsRowToResidualDOFMap(rowIdx). The mapping from the rowIdx'th row of the equations to the residual variable DOFs.
  END TYPE EquationsMappingResidualType

  !>A buffer type for pointers to a EquationsMappingResidualType \see Types::EquationsMappingResidualType
  TYPE EquationsMappingResidualPtrType
    TYPE(EquationsMappingResidualType), POINTER :: ptr !<A pointer to a EquationsMappingResidualType
  END TYPE EquationsMappingResidualPtrType
  
  !>Contains information on the equations mapping for nonlinear matrices i.e., how a field variable is mapped to residual
  !>vectors, and how the field variables are mapped to the rows and columns of the associated Jacobian matrices of the
  !>vector equations mapping. There may be multiple residual variables with a Jacobian matrix for each variable
  TYPE EquationsMappingNonlinearType
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping
    INTEGER(INTG) :: numberOfResiduals !<The number of residuals in this mapping.
    INTEGER(INTG) :: numberOfNonlinearVariables !<The number of independent variables in the nonlinear mapping
    INTEGER(INTG), ALLOCATABLE :: nonlinearVariableTypes(:) !<nonlinearVariableTypes(variableIdx). The variable type for the variableIdx'th independent nonlinear variable.
    INTEGER(INTG), ALLOCATABLE :: nonlinearVariableTypesMap(:) !<linearVariableTypesMap(variableTypeIdx). The variable index for the variableTypeIdx'th variable type. 
    TYPE(FieldVariablePtrType), ALLOCATABLE :: nonlinearVariables(:) !<nonlinearVariables(variableIdx). A pointer to the nonlinear field variable for the variableIdx'th independent nonlinear variable.
!!TODO: Should have a map from the nonlinear variable to the residuals that variable is mapped to.
    TYPE(EquationsMappingResidualPtrType), ALLOCATABLE :: residuals(:) !<residuals(residualIdx). residuals(residualIdx)%ptr is a pointer to the residualIdx'th residual.
  END TYPE EquationsMappingNonlinearType

  !>Contains information on the equations mapping for a RHS i.e., how a field variable is mapped to the RHS vector for
  !>the vector equations mapping.
  TYPE EquationsMappingRHSType
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping
    INTEGER(INTG) :: rhsVariableType !<The variable type number mapped to the RHS vector
    TYPE(FieldVariableType), POINTER :: rhsVariable !<A pointer to the variable that is mapped to the RHS vector
    TYPE(DomainMappingType), POINTER :: rhsVariableMapping !<A pointer to the RHS variable domain mapping
    REAL(DP) :: rhsCoefficient !<The multiplicative coefficient applied to the RHS vector
    INTEGER(INTG), POINTER :: rhsDOFToEquationsRowMap(:) !<rhsDOFToEquationsRowMap(rhsDOFIdx). The mapping from the rhsDofIdx'th RHS DOF in the rhs variable to the equations row.   
    INTEGER(INTG), POINTER :: equationsRowToRHSDOFMap(:) !<equationsRowToRHSDOFMap(rowIdx). The mapping from the rowIdx'th row of the equations to the RHS DOF.   
  END TYPE EquationsMappingRHSType

  !>Contains information on the equations mapping for a source vector i.e., how a field variable is mapped to the source vector for
  !>the vector equation mapping.
  TYPE EquationsMappingSourceType
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping !<A pointer to the sources mapping
    INTEGER(INTG) :: sourceNumber !<The source number for the sources mapping
    INTEGER(INTG) :: sourceVariableType !<The variable type number mapped from the source vector
    TYPE(FieldVariableType), POINTER :: sourceVariable !<A pointer to the source variable 
    TYPE(DomainMappingType), POINTER :: sourceVariableMapping !<A pointer to the domain mapping for the source variable.
    REAL(DP) :: sourceCoefficient !<The multiplicative coefficient applied to the source vector
    INTEGER(INTG), POINTER :: sourceDOFToEquationsRowMap(:) !<sourceDOFToEquationsRowMap(sourceDofIdx). The mapping from the sourceDofIdx'th source dof in the source variable to the equations row.   
    INTEGER(INTG), POINTER :: equationsRowToSourceDOFMap(:) !<equationsRowToSourceDOFMap(rowIdx). The mapping from the rowIdx'th row of the equations to the source dof.
  END TYPE EquationsMappingSourceType

  !>A buffer type for pointers to a EquationsMappingSourceType \see Types::EquationsMappingSouceType
  TYPE EquationsMappingSourcePtrType
    TYPE(EquationsMappingSourceType), POINTER :: ptr !<A pointer to the EquationsMappingSourceType
  END TYPE EquationsMappingSourcePtrType

  !>Contains information on the sources mapping 
  TYPE EquationsMappingSourcesType
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer back to the vector equations mapping
    INTEGER(INTG) :: numberOfSources !<The number of sources in the vector equations mapping
    TYPE(EquationsMappingSourcePtrType), ALLOCATABLE :: sources(:) !<sources(sourceIdx). sources(sourceIdx)%ptr is a pointer to the source mapping for the sourceIdx'th source vector.
  END TYPE EquationsMappingSourcesType

  !>Contains information on the create values cache for the scalar equations mapping. Because we do not want to allocate and
  !deallocate large data structures as the equations mapping options are changed between create start and create finish we
  !cache the important information and the allocate and process the data structures at create finish.
  TYPE EquationsMappingScalarCreateValuesCacheType
  END TYPE EquationsMappingScalarCreateValuesCacheType

  !>Contains information on the create values cache for the vector equations mapping. Because we do not want to allocate and
  !deallocate large data structures as the equations mapping options are changed between create start and create finish we
  !cache the important information and the allocate and process the data structures at create finish.
  TYPE EquationsMappingVectorCreateValuesCacheType
    INTEGER(INTG) :: numberOfDynamicMatrices !<The number of dynamic matrices in the equations mapping
    INTEGER(INTG) :: dynamicStiffnessMatrixNumber !<The dynamic matrix number corresponding to the dynamic stiffness matrix
    INTEGER(INTG) :: dynamicDampingMatrixNumber !<The dynamic matrix number corresponding to the dynamic damping matrix
    INTEGER(INTG) :: dynamicMassMatrixNumber!<The dynamic matrix number corresponding to the dynamic mass matrix
    INTEGER(INTG) :: dynamicVariableType !<The dependent variable type mapped to the dynamic equations matrices.
    REAL(DP), ALLOCATABLE :: dynamicMatrixCoefficients(:) !<dynamicsMatrixCoefficients(matrixIdx). The coefficient of the matrixIdx'th dynamic matrix in the equations set.
    INTEGER(INTG) :: numberOfLinearMatrices !<The number of linear matrices in the equations mapping
    INTEGER(INTG), ALLOCATABLE :: linearMatrixVariableTypes(:) !<linearMatrixVariableTypes(matrixIdx). The dependent variable type mapped to the matrixIdx'th linear equations matrix.
    REAL(DP), ALLOCATABLE :: linearMatrixCoefficients(:) !<linearMatrixCoefficients(matrixIdx). The coefficient of the matrixIdx'th linear matrix in the equations set.
    INTEGER(INTG) :: numberOfResiduals !<The number of residuals for the nonlinear equations
    INTEGER(INTG), ALLOCATABLE :: numberOfResidualVariables(:) !<numberOfResidualVariables(residualIdx). The number of residual variables for the residualIdx'th residual in the nonlinear equations.
    INTEGER(INTG), ALLOCATABLE :: residualVariableTypes(:,:) !<residualVariableTypes(variableIdx,residualIdx). The type of the variableIdx'th variable for the residualIdx'th residual.
    REAL(DP), ALLOCATABLE :: residualCoefficients(:) !<residualCoefficients(residualIdx). The coefficient multiplying the residualIdx'th residual vector.
    INTEGER(INTG) :: rhsVariableType !<The dependent variable type mapped to the rhs vector.
    REAL(DP) :: rhsCoefficient !<The coefficient multiplying the RHS vector.
    INTEGER(INTG) :: numberOfSources !<The number of sources in the vector equations mapping
    INTEGER(INTG), ALLOCATABLE :: sourceVariableTypes(:) !<sourceVariableTypes(sourceIdx). The source variable type mapped to the sourceIdx'th source vector
    REAL(DP), ALLOCATABLE :: sourceCoefficients(:) !<sourceCoefficients(sourceIdx). The coefficient multiplying the sourceIdx'th source vector.
  END TYPE EquationsMappingVectorCreateValuesCacheType

  !>Contains information on the vector equations LHS mapping i.e., how a field variable is mapped to the vector equations rows
  !>for this equations mapping.
  TYPE EquationsMappingLHSType
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping
    INTEGER(INTG) :: lhsVariableType !<The variable type number mapped to the RHS vector
    TYPE(FieldVariableType), POINTER :: lhsVariable !<A pointer to the variable that is mapped to the RHS vector
    INTEGER(INTG) :: numberOfRows !<The number of local rows (excluding ghost rows) in the equations 
    INTEGER(INTG) :: totalNumberOfRows !<The number of local rows (including ghost rows) in the equations 
    INTEGER(INTG) :: numberOfGlobalRows !<The number of global rows in the equations 
    TYPE(DomainMappingType), POINTER :: rowDOFsMapping !<A pointer to the row variable domain mapping
    INTEGER(INTG), POINTER :: lhsDOFToEquationsRowMap(:) !<lhsDOFToEquationsRowMap(dofIdx). The mapping from the dofIdx'th LHS dof to the vector equations row.   
    INTEGER(INTG), POINTER :: equationsRowToLHSDOFMap(:) !<equationsRowToLHSDOFMap(rowIdx). The mapping from the rowIdx'th row of the equations to the LHS dof.   
  END TYPE EquationsMappingLHSType

  !>Contains information on the mapping of field variables for a scalar equation
  TYPE EquationsMappingScalarType
    TYPE(EquationsScalarType), POINTER :: scalarEquations !<A pointer to the scalar equations for the mapping.
    LOGICAL :: scalarMappingFinished !<Is .TRUE. if the scalar mapping has been finished. .FALSE. if not.
    TYPE(EquationsMatricesScalarType), POINTER :: scalarMatrices !<A pointer to the equations scalar matrices associated with this scalar equations mapping.
    TYPE(EquationsMappingFunctionsType), POINTER :: functionMappings !<A pointer to the equations mapping for functions
    TYPE(EquationsMappingNormsType), POINTER :: normMappings !<A pointer to the equations mapping for vector norms
    TYPE(EquationsMappingDotProductsType), POINTER :: dotProductMappings !<A pointer to the equations mapping for vector dot products i.e., x^T.y
    TYPE(EquationsMappingQuadraticsType), POINTER :: quadraticMappings !<A pointer to the equations mapping for the quadratic matrices i.e., x^T.A.y
    TYPE(EquationsMappingScalarCreateValuesCacheType), POINTER :: createValuesCache !<The create values cache for the scalar equations mapping
   END TYPE EquationsMappingScalarType
  
  !>Contains information on the mapping of field variables for vector equations
  TYPE EquationsMappingVectorType
    TYPE(EquationsVectorType), POINTER :: vectorEquations !<A pointer to the vector equations for the mapping
    LOGICAL :: vectorMappingFinished !<Is .TRUE. if the vector mapping has been finished. .FALSE. if not.
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the equations vector matrices associated with this vector equations mapping.
    !Row mappings
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping !<A pointer to the LHS i.e., vector equations rows, mapping.
    !Equations mapping components
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping !<A pointer to the equations mapping for dynamic matrices i.e., M.a + C.v + K.u
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping !<A pointer to the equations mapping for the linear matrices i.e., A.x
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping !<A pointer to the equations mapping for the nonlinear matrices and vectors i.e., g(x,y)
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping !<A pointer to the equations mapping for the source vector i.e., s
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping !<A pointer to the equations mapping for the RHS vector i.e., b
    !Create values cache
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache !<The create values cache for the vector equations mapping
  END TYPE EquationsMappingVectorType

  PUBLIC VarToEquationsColumnMapType

  PUBLIC VarToEquationsMatricesMapType,VarToEquationsMatricesMapPtrType

  PUBLIC EquationsMatrixToVarMapType,EquationsMatrixToVarMapPtrType

  PUBLIC EquationsMappingFunctionType,EquationsMappingFunctionsType

  PUBLIC EquationsMappingNormType,EquationsMappingNormsType

  PUBLIC EquationsMappingDotProductType,EquationsMappingDotProductsType

  PUBLIC EquationsMappingQuadraticType,EquationsMappingQuadraticsType

  PUBLIC EquationsMappingDynamicType

  PUBLIC EquationsMappingLinearType

  PUBLIC JacobianMatrixToVarMapType,JacobianMatrixToVarMapPtrType

  PUBLIC VarToJacobianMatrixMapType,VarToJacobianMatrixMapPtrType

  PUBLIC EquationsMappingNonlinearType

  PUBLIC EquationsMappingResidualType,EquationsMappingResidualPtrType

  PUBLIC EquationsMappingRHSType

  PUBLIC EquationsMappingSourceType

  PUBLIC EquationsMappingSourcesType

  PUBLIC EquationsMappingScalarCreateValuesCacheType

  PUBLIC EquationsMappingVectorCreateValuesCacheType

  PUBLIC EquationsMappingLHSType

  PUBLIC EquationsMappingScalarType

  PUBLIC EquationsMappingVectorType
  
  !
  !================================================================================================================================
  !
  ! Equations types
  !
  
  !>Contains information on the interpolation for the equations
  TYPE EquationsInterpolationType
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations
    TYPE(FieldType), POINTER :: geometricField !<A pointer to the geometric field for the equations.
    TYPE(FieldType), POINTER :: fibreField !<A pointer to the fibre field for the equations (if one is defined).
    TYPE(FieldType), POINTER :: dependentField !<A pointer to the dependent field for the equations 
    TYPE(FieldType), POINTER :: independentField !<A pointer to the independent field for the equations 
    TYPE(FieldType), POINTER :: materialsField !<A pointer to the material field for the equations (if one is defined).
    TYPE(FieldType), POINTER :: sourceField !<A pointer to the source field for the equations (if one is defined).
    TYPE(FieldInterpolationParametersPtrType), POINTER :: geometricInterpParameters(:) !<GEOMETRIC_INTERP_PARAMETERS(variableTypeIdx). A pointer to the variableTypeIdx'th geometric interpolation parameters for the equations.
    TYPE(FieldInterpolationParametersPtrType), POINTER :: fibreInterpParameters(:) !<FIBRE_INTERP_PARAMETERS(variableTypeIdx). A pointer to the fibre interpolation parameters for the equations (if a fibre field is defined). 
    TYPE(FieldInterpolationParametersPtrType), POINTER :: dependentInterpParameters(:) !<DEPENDENT_INTERP_PARAMETERS(variableTypeIdx). A pointer to the variableTypeIdx'th dependent interpolation parameters for the equations. 
    TYPE(FieldInterpolationParametersPtrType), POINTER :: prevDependentInterpParameters(:) !<prevDependentInterpParameters(fieldVariableType). A pointer to the previous fieldVariableType'th dependent interpolation parameters for the equations. Only allocated for non-static equations.
    TYPE(FieldInterpolationParametersPtrType), POINTER :: independentInterpParameters(:) !<INDEPENDENT_INTERP_PARAMETERS(variableTypeIdx). A pointer to the variableTypeIdx'th independent interpolation parameters for the equations. 
    TYPE(FieldInterpolationParametersPtrType), POINTER :: materialsInterpParameters(:) !<materialsInterpParameters(variableTypeIdx). A pointer to the variableTypeIdx'th material interpolation parameters for the equations (if a material field is defined). 
    TYPE(FieldInterpolationParametersPtrType), POINTER :: sourceInterpParameters(:) !<sourceInterpParameters(variableTypeIdx). A pointer to the variableTypeIdx'th source interpolation parameters for the equations (if a source field is defined). 
    TYPE(FieldInterpolatedPointPtrType), POINTER :: geometricInterpPoint(:) !<geometricInterpPoint(variableTypeIdx). A pointer to the variableTypeIdx'th geometric interpolated point information for the equations. 
    TYPE(FieldInterpolatedPointPtrType), POINTER :: fibreInterpPoint(:) !<fibreInterpPoint(variableTypeIdx). A pointer to the variableTypeIdx'th fibre interpolated point information for the equations (if a fibre field is defined). 
    TYPE(FieldInterpolatedPointPtrType), POINTER :: dependentInterpPoint(:) !<dependentInterpPoint(variableTypeIdx). A pointer to the variableTypeIdx'th dependent interpolated point information for the equations. 
    TYPE(FieldInterpolatedPointPtrType), POINTER :: prevDependentInterpPoint(:) !<prevDependentInterpPoint(fieldVariableType). A pointer to the previous fieldVariableType'th dependent interpolated point information for the equations. Only allocated for non-static equations.
    TYPE(FieldInterpolatedPointPtrType), POINTER :: independentInterpPoint(:) !<independentInterpPoint(variableTypeIdx). A pointer to the variableTypeIdx'th independent interpolated point information for the equations. 
    TYPE(FieldInterpolatedPointPtrType), POINTER :: materialsInterpPoint(:) !<materialsInterpPoint(variableTypeIdx). A pointer to the variableTypeIdx'th material interpolated point information for the equations (if a material field is defined). 
    TYPE(FieldInterpolatedPointPtrType), POINTER :: sourceInterpPoint(:) !<sourceInterpPoint(variableTypeIdx). A pointer to the variableTypeIdx'th source interpolated point information for the equations (if a source field is defined).
    TYPE(FieldPhysicalPointPtrType), POINTER :: dependentPhysicalPoint(:) !<dependentPhysicalPoint(variableTypeIdx). A pointer to the variableTypeIdx'th dependent physical point information for the equations. 
    TYPE(FieldInterpolatedPointMetricsPtrType), POINTER :: dependentInterpPointMetrics(:) !<dependentInterpPointMetrics(variableTypeIdx). A pointer to the variableTypeIdx'th dependent interpolated point metrics information 
    TYPE(FieldInterpolatedPointMetricsPtrType), POINTER :: prevDependentInterpPointMetrics(:) !<prevDependentInterpPointMetrics(fieldVariableType). A pointer to the previous fieldVariableType'th dependent interpolated point metrics information. Only allocated for non-static equations. 
    TYPE(FieldInterpolatedPointMetricsPtrType), POINTER :: independentInterpPointMetrics(:) !<independentInterpPointMetrics(variableTypeIdx). A pointer to the variableTypeIdx'th independent interpolated point metrics information 
    TYPE(FieldInterpolatedPointMetricsPtrType), POINTER :: geometricInterpPointMetrics(:) !<geometricInterpPointMetrics(variableTypeIdx). A pointer to the variableTypeIdx'th geometric interpolated point metrics information 
    TYPE(FieldInterpolatedPointMetricsPtrType), POINTER :: fibreInterpPointMetrics(:) !<fibreInterpPointMetrics(variableTypeIdx). A pointer to the variableTypeIdx'th fibre interpolated point metrics information 
  END TYPE EquationsInterpolationType

  !>Contains information about scalar equations (i.e., a single equation row).
  TYPE EquationsScalarType
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations
    TYPE(EquationsMappingScalarType), POINTER :: scalarMapping !<A pointer to the mapping for the scalar equation
    TYPE(EquationsMatricesScalarType), POINTER :: scalarMatrices !<A pointer to the matrices, vectors and scalars for the scalar equation
  END TYPE EquationsScalarType

  !>Contains information about vector equations (i.e., a number of equation rows). 
  TYPE EquationsVectorType
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the mapping for the vector equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<A pointer to the matrices, vectors for the vector equations
  END TYPE EquationsVectorType
  
  !>Contains information about the equations in an equations set. \see OpenCMISS::Iron::cmfe_EquationsType
  TYPE EquationsType
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations_set
    LOGICAL :: equationsFinished !<Is .TRUE. if the equations have been finished, .FALSE. if not.
    INTEGER(INTG) :: equationType !<The equations type \see EquationsRoutines_EquationTypes,EquationsRoutines
    INTEGER(INTG) :: equalityType !<The equations equality type \see EquationsRoutines_EquationEqualityTypes,EquationsRoutines
    INTEGER(INTG) :: linearity !<The equations linearity type \see EquationsSetRoutines_LinearityTypes,EquationsSetRoutines
    INTEGER(INTG) :: timeDependence !<The equations time dependence type \see EquationsSetRoutines_TimeDependenceTypes,EquationsSetRoutines
    INTEGER(INTG) :: outputType !<The output type for the equations \see EquationsRoutines_EquationsOutputTypes,EquationsRoutines
    INTEGER(INTG) :: sparsityType !<The sparsity type for the equation matrices of the equations \see EquationsRoutines_EquationsSparsityTypes,EquationsRoutines
    INTEGER(INTG) :: lumpingType !<The lumping type for the equation matrices of the equations \see EquationsRoutines_EquationsLumpingTypes,EquationsRoutines
    TYPE(EquationsInterpolationType), POINTER :: interpolation !<A pointer to the interpolation information used in the equations.
    TYPE(EquationsScalarType), POINTER :: scalarEquations !<A pointer to the scalar equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations !<A pointer to the vector equations
  END TYPE EquationsType

  TYPE EquationsPtrType
    TYPE(EquationsType), POINTER :: ptr
  END TYPE EquationsPtrType

  PUBLIC EquationsInterpolationType

  PUBLIC EquationsScalarType

  PUBLIC EquationsVectorType

  PUBLIC EquationsType,EquationsPtrType

  !
  !================================================================================================================================
  !
  ! Boundary conditions types
  !

  !>Contains information on the boundary conditions for the row variables
  TYPE BoundaryConditionsRowVariableType
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions for this boundary conditions variable
    TYPE(FieldVariableType), POINTER :: lhsVariable !<A pointer to the row field variable for this boundary condition variable
    TYPE(FieldVariableType), POINTER :: rhsVariable !<A pointer to the RHS field variable for this LHS variable. Temporary for now until BC restructure.
    INTEGER(INTG) :: numberOfDOFs !<The number of DOFs in the row variable
    INTEGER(INTG) :: totalNumberOfDOFs !<The total number of DOFs in the row variable
    INTEGER(INTG) :: numberOfGlobalDOFs !<The number of global DOFs in the row variable. TEMP WHILST WE SWITCH FROM GLOBAL DOFS.
    INTEGER(INTG), ALLOCATABLE :: rowConditionTypes(:) !<rowConditionTypes(dofIdx). The row DOF boundary condition type. \see BoundaryConditionsRoutines_RowConditionsTypes
  END TYPE BoundaryConditionsRowVariableType

  !>A buffer type to allow for an array of pointers to a BoundaryConditionsRowVariableType \see Types::BoundaryConditionsRowVariableType
  TYPE BoundaryConditionsRowVariablePtrType
    TYPE(BoundaryConditionsRowVariableType), POINTER :: ptr !<A pointer to the boundary conditions row variable
  END TYPE BoundaryConditionsRowVariablePtrType

  !>Contains information on the boundary conditions for a dependent field variable
  TYPE BoundaryConditionsVariableType
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions for this boundary conditions variable
    INTEGER(INTG) :: variableType !<The type of variable for this variable boundary conditions
    TYPE(FieldVariableType), POINTER :: variable !<A pointer to the field variable for this boundary condition variable
    INTEGER(INTG), ALLOCATABLE :: DOFTypes(:) !<DOFTypes(dofIdx). The general boundary condition type (eg. fixed or free) of the dofIdx'th dof in the dependent field variable. \see  BoundaryConditionsRoutines_DOFTypes
    INTEGER(INTG), ALLOCATABLE :: conditionTypes(:) !<conditionTypes(dofIdx). The specific boundary condition type (eg. incremented pressure) of the dofIdx'th dof of the dependent field variable, which might be specific to an equation set. The solver routines should not need to use this array, and should only need the DOFTypes array. \see OpenCMISS_BoundaryConditionsDOFTypes,OpenCMISS
    INTEGER(INTG) :: numberOfDirichletConditions !<Stores the number of dirichlet conditions associated with this variable
    TYPE(BoundaryConditionsDirichletType), POINTER :: dirichletBoundaryConditions !<A pointer to the dirichlet boundary condition type for this boundary condition variable
    TYPE(BoundaryConditionsNeumannType), POINTER :: neumannBoundaryConditions !<A pointer to the Neumann boundary condition type for this boundary condition variable
    TYPE(BoundaryConditionsPressureIncrementedType), POINTER :: pressureIncrementedBoundaryConditions !<A pointer to the pressure incremented condition type for this boundary condition variable
    INTEGER(INTG), ALLOCATABLE :: dofCounts(:) !<dofCounts(CONDITION_TYPE): The number of DOFs that have a CONDITION_TYPE boundary condition set
    LOGICAL, ALLOCATABLE :: parameterSetRequired(:) !<parameterSetRequired(PARAMETER_SET) is true if any boundary condition has been set that requires the PARAMETER_SET field parameter set
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints !<A pointer to the linear DOF constraints structure.
  END TYPE BoundaryConditionsVariableType

  !>A buffer type to allow for an array of pointers to a BoundaryConditionsVariableType \see Types::BoundaryConditionsVariableType
  TYPE BoundaryConditionsVariablePtrType
    TYPE(BoundaryConditionsVariableType), POINTER :: ptr !<A pointer to the boundary conditions variable
  END TYPE BoundaryConditionsVariablePtrType

  !>Contains information on the boundary conditions for the solver equations. \see OpenCMISS::Iron::cmfe_BoundaryConditionsType
  TYPE BoundaryConditionsType
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations.
    LOGICAL :: boundaryConditionsFinished !<Is .TRUE. if the boundary conditions for the equations set has finished being created, .FALSE. if not.
    INTEGER(INTG) :: numberOfBoundaryConditionsRowVariables !<The number of boundary conditions row variables
    TYPE(BoundaryConditionsRowVariablePtrType), ALLOCATABLE :: boundaryConditionsRowVariables(:) !<boundaryConditionsRowVariables(variableIdx). boundaryConditionsRowVariables(variableIdx)%ptr is the pointer to the variableIdx'th boundary conditions row variable. variableIdx ranges from 1 to numberOfBoundaryConditionsRowVariables    
    INTEGER(INTG) :: numberOfBoundaryConditionsVariables !<The number of boundary conditions variables
    TYPE(BoundaryConditionsVariablePtrType), ALLOCATABLE :: boundaryConditionsVariables(:) !<boundaryConditionsVariables(variableIdx). boundaryConditionsVariables(variableIdx)%ptr is the pointer to the variableIdx'th boundary conditions variable. variableIdx ranges from 1 to numberOfBoundaryConditionsVariables
    INTEGER(INTG) :: neumannMatrixSparsity !<The sparsity type of the Neumann integration matrices. \see SolverRoutines_EquationsSparsityTypes,SolverRoutines
  END TYPE BoundaryConditionsType

  !>A buffer type to allow for an array of pointers to a BoundaryConditionsSparsityIndicesType \see Types::BoundaryConditionsSparsityIndicesType
  TYPE BoundaryConditionsSparsityIndicesPtrType
    TYPE(BoundaryConditionsSparsityIndicesType), POINTER :: ptr !<A pointer to the boundary conditions sparsity indices type
  END TYPE BoundaryConditionsSparsityIndicesPtrType

  !> Contains information on dofs with associated dirichlet conditions and corresponding non-zero elements in the equations matrices
  TYPE BoundaryConditionsDirichletType
    INTEGER(INTG), ALLOCATABLE :: dirichletDOFIndices(:)  !<dirichletDOFIndices(dofIdx). Stores the dofIdx of the dofs which are subject to a dirichlet boundary condition \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
    TYPE(BoundaryConditionsSparsityIndicesPtrType), ALLOCATABLE :: linearSparsityIndices(:,:) !<linearSparsityIndices(equationsSetIdx,equationsMatrixIdx). Stores the indices of the non-zero elements of the equationsSetIdx'th equation set and equationsMatrixIdx'th linear equation matrix in the columns corresponding to the dofs which are subject to a dirichlet boundary condition
    TYPE(BoundaryConditionsSparsityIndicesPtrType), ALLOCATABLE :: dynamicSparsityIndices(:,:) !<dynamicSparsityIndices(equationsSetIdx,equationsMatrixIdx). Stores the indices of the non-zero elements of the equationsSetIdx'th equation set and equationsMatrixIdx'th dynamic equation matrix in the columns corresponding to the dofs which are subject to a dirichlet boundary condition
  END TYPE BoundaryConditionsDirichletType

  !> Contains information on indices of non-zero elements with associated dirichlet conditions
  !> Indices stored in compressed column format without a values array
  TYPE BoundaryConditionsSparsityIndicesType
    INTEGER(INTG), ALLOCATABLE :: sparseRowIndices(:) !<sparseRowIndices(sparseColumnIndices(columnIdx)). Between sparseColumnIndices(columnIdx) and sparseColumnIndices(columnIdx+1)-1 are the row indices of non-zero elements of the 'columnIdx'th column
    INTEGER(INTG), ALLOCATABLE :: sparseColumnIndices(:) !<sparseColumnIndices(columnIdx). Between sparseColumnIndices(columnIdx) and sparseColumnIndices(columnIdx+1)-1 are the row indices of non-zero elements of the 'columnIdx'th column
  END TYPE BoundaryConditionsSparsityIndicesType

  !>Contains information used to integrate Neumann boundary conditions
  TYPE BoundaryConditionsNeumannType
    INTEGER(INTG), ALLOCATABLE :: setDofs(:) !<setDofs(neumann_idx): the global dof for the neumann_idx'th Neumann condition
    TYPE(DistributedMatrixType), POINTER :: integrationMatrix !<The N matrix that multiples the point values vector q to give the integrated values f. Number of rows equals number of local dofs, and number of columns equals number of set point DOFs.
    TYPE(DistributedVectorType), POINTER :: pointValues !<The vector of set point values q
    TYPE(DomainMappingType), POINTER :: pointDofMapping !<The domain mapping for DOFs with Neumann point conditions set.
  END TYPE BoundaryConditionsNeumannType

  !>Contains information on dofs associated with pressure incremented conditions
  TYPE BoundaryConditionsPressureIncrementedType
    INTEGER(INTG), ALLOCATABLE :: pressureIncrementedDOFIndices(:)  !<pressureIncrementedDOFIndices(idx). Stores the dofIdx of the dofs which are subject to a pressure incremented boundary condition \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  END TYPE BoundaryConditionsPressureIncrementedType

  !>Describes the value of a DOF as a linear combination of other DOFs.
  TYPE BoundaryConditionsDOFConstraintType
    INTEGER(INTG) :: numberOfDOFs !<The number of other DOFs that contribute to the value of this DOF.
    INTEGER(INTG) :: globalDOF !<The global DOF number of this DOF.
    INTEGER(INTG), ALLOCATABLE :: dofs(:) !<dofs(dofIdx) is the dofIdx'th DOF in the linear combination.
    REAL(DP), ALLOCATABLE :: coefficients(:) !<coefficients(dofIdx) is the linear proportion of the dofIdx'th dof in the value for this DOF.
  END TYPE BoundaryConditionsDOFConstraintType

  !>A pointer to a linear DOF constraint.
  TYPE BoundaryConditionsDOFConstraintPtrType
    TYPE(BoundaryConditionsDOFConstraintType), POINTER :: ptr
  END TYPE BoundaryConditionsDOFConstraintPtrType

  !>The coupled equations DOF information for the DOF constraints.
  !>The BoundaryConditionsDofConstraintType describes how an
  !>equations DOF is a linear combination of other equations DOFs.
  !>This data structure describes a solver row or column that is mapped to multiple
  !>equations rows or columns and is used to help build the solver mapping.
  !>The first equation row/column is used as the owner of the solver row/column.
  TYPE BoundaryConditionsCoupledDOFsType
    INTEGER(INTG) :: numberOfDOFs !<The number of equations rows or columns that are mapped to this solver DOF.
    INTEGER(INTG), ALLOCATABLE :: globalDOFs(:) !<globalDOFs(dofIdx) is the dofIdx'th global DOF mapped to the row/column.
    INTEGER(INTG), ALLOCATABLE :: localDOFs(:) !<localDOFs(dofIdx) is the dofIdx'th local DOF mapped to the row/column.
    REAL(DP), ALLOCATABLE :: coefficients(:) !<coefficients(dofIdx) is the linear proportion of the dofIdx'th dof in the mapping for this row or column.
  END TYPE BoundaryConditionsCoupledDOFsType

  !>A pointer to the coupled equations DOF information for the DOF constraints.
  TYPE BoundaryConditionsCoupledDOFsPtrType
    TYPE(BoundaryConditionsCoupledDOFsType), POINTER :: ptr
  END TYPE BoundaryConditionsCoupledDOFsPtrType

  !>Describes linear constraints between solver DOFs in the solver mapping.
  TYPE BoundaryConditionsDOFConstraintsType
    INTEGER(INTG) :: numberOfConstraints !<The number of DOF constraints.
    INTEGER(INTG) :: numberOfDOFs !<The number of global DOFs.
    TYPE(BoundaryConditionsDOFConstraintPtrType), ALLOCATABLE :: constraints(:) !<constraints(constraintIdx)%ptr is a pointer to the dof constraint for the constraintIdx'th constraint.
    TYPE(BoundaryConditionsCoupledDOFsPtrType), ALLOCATABLE :: dofCouplings(:) !<dofCouplings(dofIdx)%ptr is a pointer to the coupled DOF information for the solver row/column corresponding to the dofIdx'th equations DOF.
  END TYPE BoundaryConditionsDOFConstraintsType

  PUBLIC BoundaryConditionsRowVariableType,BoundaryConditionsRowVariablePtrType

  PUBLIC BoundaryConditionsVariableType,BoundaryConditionsVariablePtrType

  PUBLIC BoundaryConditionsType

  PUBLIC BoundaryConditionsSparsityIndicesType,BoundaryConditionsSparsityIndicesPtrType

  PUBLIC BoundaryConditionsDirichletType

  PUBLIC BoundaryConditionsNeumannType

  PUBLIC BoundaryConditionsPressureIncrementedType

  PUBLIC BoundaryConditionsDOFConstraintType,BoundaryConditionsDOFConstraintPtrType

  PUBLIC BoundaryConditionsCoupledDOFsType,BoundaryConditionsCoupledDOFsPtrType

  PUBLIC BoundaryConditionsDOFConstraintsType

  !
  !================================================================================================================================
  !
  ! Equations set types
  !

  !>Contains information on the setup information for an equations set
  TYPE EquationsSetSetupType
    INTEGER(INTG) :: setupType !<The setup type for the equations set setup \see EquationsSetRoutines_SetupTypes,EquationsSetRoutines
    INTEGER(INTG) :: actionType !<The action type for the equations set setup \see EquationsSetRoutines_SetupActionTypes,EquationsSetRoutines
    INTEGER(INTG) :: fieldUserNumber !<The user number for the field for the equations set setup.
    TYPE(FieldType), POINTER :: field !<A pointer to the field for the equations set setup.
    INTEGER(INTG) :: analyticFunctionType !<The analytic function type to use.
  END TYPE EquationsSetSetupType  
  
  !>Contains information on the geometry for an equations set
  TYPE EquationsSetGeometryType
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set.
    TYPE(FieldType), POINTER :: geometricField !<The geometric field for this equations set.
    TYPE(FieldType), POINTER :: fibreField !<The fibre field for this equations set if one is defined. If no fibre field is defined the pointer is NULL.
  END TYPE EquationsSetGeometryType

  TYPE EquationsSetMaterialsType
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set.
    LOGICAL :: materialsFinished !<Is .TRUE. if the materials for the equations set has finished being created, .FALSE. if not.
    LOGICAL :: materialsFieldAutoCreated !<Is .TRUE. if the materials field has been auto created, .FALSE. if not.
    TYPE(FieldType), POINTER :: materialsField !<A pointer to the materials field for the equations set if one is defined. If no material field is defined the pointer is NULL.
  END TYPE EquationsSetMaterialsType

  !>Contains information on the dependent variables for the equations set.
  TYPE EquationsSetDependentType
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set.
    LOGICAL :: dependentFinished !<Is .TRUE. if the dependent variables for the equations set has finished being created, .FALSE. if not.
    LOGICAL :: dependentFieldAutoCreated !<Is .TRUE. if the dependent field has been auto created, .FALSE. if not.
    TYPE(FieldType), POINTER :: dependentField !<A pointer to the dependent field for the equations set.
  END TYPE EquationsSetDependentType

  !>Contains information on the derived variables for the equations set, eg. stress or strain
  TYPE EquationsSetDerivedType
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer back to the equations set.
    LOGICAL :: derivedFinished !<Is .TRUE. if the derived variables for the equations set have finished being created, .FALSE. if not.
    LOGICAL :: derivedFieldAutoCreated !<Is .TRUE. if the derived field has been or will be auto created by the equations set setup, .FALSE. if it was already created.
    TYPE(FieldType), POINTER :: derivedField !<A pointer to the derived field for the equations set.
    INTEGER(INTG) :: numberOfVariables !<The number of derived variables used.
    INTEGER(INTG), ALLOCATABLE :: variableTypes(:) !<variableTypes(derivedType) is zero if the derived type is not used, otherwise it is the field variable type in the derived field for the derived variable type
  END TYPE EquationsSetDerivedType

  !>Contains information on the independent variables for the equations set.
  TYPE EquationsSetIndependentType
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set.
    LOGICAL :: independentFinished !<Is .TRUE. if the independent variables for the equations set has finished being created, .FALSE. if not.
    LOGICAL :: independentFieldAutoCreated !<Is .TRUE. if the independent field has been auto created, .FALSE. if not.
    TYPE(FieldType), POINTER :: independentField !<A pointer to the independent field for the equations set.
  END TYPE EquationsSetIndependentType

  !>Contains information on the source for the equations set.
  TYPE EquationsSetSourceType
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set.
    LOGICAL :: sourceFinished !<Is .TRUE. if the source for the equations set has finished being created, .FALSE. if not.
    LOGICAL :: sourceFieldAutoCreated !<Is .TRUE. if the source field has been auto created, .FALSE. if not.
    TYPE(FieldType), POINTER :: sourceField !<A pointer to the source field for the equations set if one is defined. If no source is defined the pointer is NULL.
  END TYPE EquationsSetSourceType
  
  !>Contains information on the analytic setup for the equations set.
  TYPE EquationsSetAnalyticType
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set.
    INTEGER(INTG) :: analyticFunctionType !<The analytic function identifier
    LOGICAL :: analyticFinished !<Is .TRUE. if the analytic setup for the problem has finished being created, .FALSE. if not.
    LOGICAL :: analyticFieldAutoCreated !<Is .TRUE. if the analytic field has been auto created, .FALSE. if not.
    TYPE(FieldType), POINTER :: analyticField !<A pointer to the analytic field for the equations set if one is defined. If no source is defined the pointer is NULL.
    REAL(DP) :: analyticTime !<The time value to use for analytic evaluations.
    REAL(DP) :: analyticUserParams(20)  !<A small array that can be used to hold various parameters often required in analytic problems. \todo should this be allocated?
  END TYPE EquationsSetAnalyticType

  TYPE EquationsSetEquationsFieldType
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set.
    LOGICAL :: equationsSetFieldFinished !<Is .TRUE. if the equations set field for the equations set has finished being created, .FALSE. if not.
    LOGICAL :: equationsSetFieldAutoCreated !<Is .TRUE. if the equations set field has been auto created, .FALSE. if not.
    TYPE(FieldType), POINTER :: equationsSetField !<A pointer to the equations set field for the equations set.
  END TYPE EquationsSetEquationsFieldType

  !>Contains information on an equations set. \see OpenCMISS::Iron::cmfe_EquationsSetType
  TYPE EquationsSetType
    INTEGER(INTG) :: userNumber !<The user identifying number of the equations set
    INTEGER(INTG) :: globalNumber !<The global index of the equations set in the region.
    LOGICAL :: equationsSetFinished !<Is .TRUE. if the equations set have finished being created, .FALSE. if not.
    TYPE(EquationsSetsType), POINTER :: equationsSets !<A pointer back to the equations sets
    TYPE(VARYING_STRING) :: label !<A user defined label for the equations set.
    TYPE(RegionType), POINTER :: region !<A pointer back to the region containing the equations set.
    INTEGER(INTG), ALLOCATABLE :: specification(:) !<The equations set specification array, eg. [class, type, subtype], although there can be more or fewer identifiers. Unused identifiers are set to zero.
    INTEGER(INTG) :: specificationLength !<The number of non-zero identifiers in the specification.
    REAL(DP) :: currentTime !<The current time for the equations set
    REAL(DP) :: deltaTime !<The current time increment for the equations set
    INTEGER(INTG) :: outputType !<The output type for the equations set \see EquationsSetRoutines_OutputTypes,EquationsSetRoutines
    INTEGER(INTG) :: solutionMethod !<The solution method for the equations set \see EquationsRoutines_SolutionMethods 
    TYPE(EquationsSetGeometryType) :: geometry !<The geometry information for the equations set.
    TYPE(EquationsSetMaterialsType), POINTER :: materials !<A pointer to the materials information for the equations set.
    TYPE(EquationsSetSourceType), POINTER :: source !<A pointer to the source information for the equations set.
    TYPE(EquationsSetDependentType) :: dependent !<The depedent variable information for the equations set.
    TYPE(EquationsSetIndependentType), POINTER :: independent !<A pointer to the indepedent field information for the equations set.
    TYPE(EquationsSetAnalyticType), POINTER :: analytic !<A pointer to the analytic setup information for the equations set.
    TYPE(EquationsSetDerivedType), POINTER :: derived !<A pointer to the derived field information.
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations information for the equations set
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary condition information for the equations set. THIS SHOULD BE CHANGED AS EQUATIONS SETS SHOULDN'T HAVE BOUNDARY CONDITIONS
    TYPE(EquationsSetEquationsFieldType), POINTER :: equationsField !<A pointer to the equations field for the equations set.
  END TYPE EquationsSetType
  
  !>A buffer type to allow for an array of pointers to a EquationsSetType \see Types::EquationsSetType
  TYPE EquationsSetPtrType
    TYPE(EquationsSetType), POINTER :: ptr !<A pointer to the equations set.
  END TYPE EquationsSetPtrType
  
  TYPE EquationsSetsType
    TYPE(RegionType) , POINTER :: region !<A pointer to the region containing the equations sets
    INTEGER(INTG) :: numberOfEquationsSets!<The number of equations sets
    TYPE(EquationsSetPtrType), POINTER :: equationsSets(:) !<EquationsSets(equationsSetIdx). EquationsSets(equationsSetIdx)%ptr is the pointer to the equationsSetIdx'th equations set.
  END TYPE EquationsSetsType

  PUBLIC EquationsSetSetupType

  PUBLIC EquationsSetGeometryType

  PUBLIC EquationsSetMaterialsType

  PUBLIC EquationsSetDependentType

  PUBLIC EquationsSetDerivedType

  PUBLIC EquationsSetIndependentType

  PUBLIC EquationsSetSourceType

  PUBLIC EquationsSetAnalyticType

  PUBLIC EquationsSetEquationsFieldType

  PUBLIC EquationsSetType,EquationsSetPtrType,EquationsSetsType
  
  !
  !================================================================================================================================
  !
  ! Interface types
  
  !>Contains information about an interface matrix.
  TYPE InterfaceMatrixType
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer to the interface matrices for the interface matrix.
    INTEGER(INTG) :: matrixNumber !<The number of the interface matrix
    REAL(DP) :: matrixCoefficient !<The coefficient applied to the matrix
    REAL(DP) :: transposeMatrixCoefficient !<The coefficient applied to the transpose matrix
    INTEGER(INTG) :: storageType !<The storage (sparsity) type for this matrix
    INTEGER(INTG) :: structureType !<The structure (sparsity) type for this matrix
    INTEGER(INTG) :: numberOfRows !<The number of rows in this interface matrix
    INTEGER(INTG) :: totalNumberOfRows !<The number of rows in this interface matrix
    INTEGER(INTG) :: interfaceMatrixTimeDependenceType !<Determines where the interface matrix is mapped to
    INTEGER(INTG) :: interfaceMatrixTransposeTimeDependenceType !<Determines where the transpose of the interface matrix is mapped to
    LOGICAL :: updateMatrix !<Is .TRUE. if this interface matrix is to be updated
    LOGICAL :: firstAssembly !<Is .TRUE. if this interface matrix has not been assembled
    LOGICAL :: hasTranspose !<Is .TRUE. if this interface matrix has has transpose
    TYPE(DistributedMatrixType), POINTER :: matrix !<A pointer to the distributed interface matrix data
    TYPE(DistributedMatrixType), POINTER :: matrixTranspose !<A pointer to the distributed interface matrix transpose data
    TYPE(DistributedVectorType), POINTER :: tempVector !<Temporary vector used for assembly. 
    TYPE(DistributedVectorType), POINTER :: tempTransposeVector !<Temporary vector used for assembly. 
    TYPE(ElementMatrixType) :: elementMatrix !<The element matrix for this interface matrix
  END TYPE InterfaceMatrixType

  !>A buffer type to allow for an array of pointers to a InterfaceMatrixType \see Types::InterfaceMatrixType.
  TYPE InterfaceMatrixPtrType
    TYPE(InterfaceMatrixType), POINTER :: ptr !<A pointer to the interface matrix.
  END TYPE InterfaceMatrixPtrType

  !>Contains information of the RHS vector for interface matrices
  TYPE InterfaceRHSType
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer back to the interface matrices.
    REAL(DP) :: rhsCoefficient !<The coefficient applied to the RHS vetor
    LOGICAL :: updateVector !<Is .TRUE. if the interface rhs vector is to be updated
    LOGICAL :: firstAssembly !<Is .TRUE. if this rhs vector has not been assembled
    TYPE(DistributedVectorType), POINTER :: rhsVector !<A pointer to the distributed RHS vector data 
    TYPE(ElementVectorType) :: elementVector !<The element rhs information
  END TYPE InterfaceRHSType
  
  !>Contains information on the interface matrices 
  TYPE InterfaceMatricesType
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer back to the interface equations
    LOGICAL :: interfaceMatricesFinished !<Is .TRUE. if the interface  matrices have finished being created, .FALSE. if not.
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface equations mapping for the interface equations matrices.
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping for the interface equations matrices
    INTEGER(INTG) :: numberOfColumns !<The number of local columns in the interface matrices
    INTEGER(INTG) :: totalNumberOfColumns !<The total number of local columns in the interface matrices
    INTEGER(INTG) :: numberOfGlobalColumns !<The number of global columns in the interface matrices
    INTEGER(INTG) :: numberOfInterfaceMatrices !<The number of interfaces matrices defined for the interface condition.
    TYPE(InterfaceMatrixPtrType), ALLOCATABLE :: matrices(:) !<matrices(matrixIdx)%ptr contains the information on the matrixIdx'th  interface matrix.
    TYPE(InterfaceRHSType), POINTER :: rhsVector !<A pointer to the RHS vector information for the interface matrices.
  END TYPE InterfaceMatricesType

  !>Contains information on interface variable mapping for an interface matrix.
  TYPE InterfaceMatrixToVarMapType
    INTEGER(INTG) :: matrixNumber !<The interface matrix number
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix !<A pointer to the interface matrix
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set containing the dependent variable that is mapped to this interface matrix.
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface condition containing the Lagrange variable that is mapped to this interface matrix.
    INTEGER(INTG) :: variableType !<The dependent variable type mapped to this interface matrix
    TYPE(FieldVariableType), POINTER :: variable !<A pointer to the field variable that is mapped to this interface matrix
    INTEGER(INTG) :: meshIndex !<The mesh index for the matrix in the interface.
    LOGICAL :: hasTranspose !<.TRUE. if the interface matrix has a tranpose, .FALSE. if not.  
    REAL(DP) :: matrixCoefficient !<The multiplicative coefficent for the matrix.    
    REAL(DP) :: transposeMatrixCoefficient !<The multiplicative coefficent for the matrix transpose.    
    INTEGER(INTG) :: numberOfRows !<The number of rows  in this interface matrix.
    INTEGER(INTG) :: totalNumberOfRows !<The total number of rows in this interface matrix.
    INTEGER(INTG) :: numberOfGlobalRows !<The global number of rows in this interface matrix.
    TYPE(DomainMappingType), POINTER :: rowDOFsMapping !<A pointer to the domain mapping for the row dofs.
    INTEGER(INTG), POINTER :: variableDOFToRowMap(:) !<variableDOFToRowMap(dofIdx). The mapping from the dofIdx'th variable dof to the rows on the interface matrix
  END TYPE InterfaceMatrixToVarMapType

  !>A buffer type to allow for arrays of pointers to InterfaceMatrixToVarMapType \see Types::InterfaceMatrixToVarMapType
  TYPE InterfaceMatrixToVarMapPtrType
    TYPE(InterfaceMatrixToVarMapType), POINTER :: ptr !<A pointer to a InterfaceMatrixToVarMapType
  END TYPE InterfaceMatrixToVarMapPtrType
  
  TYPE InterfaceMappingRHSType
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer back to the interface mapping
    INTEGER(INTG) :: rhsVariableType !<The variable type number mapped to the RHS vector
    TYPE(FieldVariableType), POINTER :: rhsVariable !<A pointer to the variable that is mapped to the RHS vector
    TYPE(DomainMappingType), POINTER :: rhsVariableMapping !<A pointer to the RHS variable domain mapping
    REAL(DP) :: rhsCoefficient !<The multiplicative coefficient applied to the RHS vector
    INTEGER(INTG), POINTER :: rhsDOFToInterfaceRowMap(:) !<rhsDOFToInterfaceRowMap(rhsDOFIdx). The mapping from the rhsDOFIdx'th RHS dof in the rhs variable to the interface row.   
    INTEGER(INTG), POINTER :: interfaceRowToRHSDOFMap(:) !<interfaceRowToRHSDOFMap(rowIdx). The mapping from the rowIdx'th row of the interface to the RHS dof.   
  END TYPE InterfaceMappingRHSType
  
  TYPE InterfaceMappingCreateValuesCacheType
    INTEGER(INTG) :: numberOfInterfaceMatrices !<Cache of the number of interface matrices
    INTEGER(INTG) :: lagrangeVariableType !<The variable type of the Lagrange field.
    REAL(DP), ALLOCATABLE :: matrixCoefficients(:) !<matrixCoefficients(matrixIdx). The matrix cooefficient for the matrixIdx'th interface matrix.
    REAL(DP), ALLOCATABLE :: transposeMatrixCoefficients(:) !<transposeMatrixCoefficients(matrixIdx). The matrix cooefficient for the matrixIdx'th transpose interface matrix.
    LOGICAL, ALLOCATABLE :: hasTranspose(:) !<hasTranspose(matrixIdx). .TRUE. if the matrixIdx'th interface matrix has an tranpose, .FALSE. if not.
    INTEGER(INTG), ALLOCATABLE :: matrixRowFieldVariableIndices(:) !<matrixRowFieldVariableIndices(variableIdx). The field variable index that are mapped to the the interface matrix rows.
    INTEGER(INTG), ALLOCATABLE :: matrixColFieldVariableIndices(:) !<matrixColFieldVariableIndices(variableIdx). The field variable index that are mapped to the the interface matrix columns.
    INTEGER(INTG) :: rhsLagrangeVariableType !<The Lagrange variable type mapped to the rhs vector
    REAL(DP) :: rhsCoefficient !<The coefficient multiplying the RHS vector.
  END TYPE InterfaceMappingCreateValuesCacheType
  
  !>Contains information on an interface mapping. TODO: Generalise to non-Lagrange multipler mappings
  TYPE InterfaceMappingType  
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations  !<A pointer to the interface equations for this interface mapping
    LOGICAL :: interfaceMappingFinished !<Is .TRUE. if the interface mapping has finished being created, .FALSE. if not.
    INTEGER(INTG) :: lagrangeVariableType !<The variable type of the mapped Lagrange field variable.
    TYPE(FieldVariableType), POINTER :: lagrangeVariable !<A pointer to the variable that is mapped to the Lagrange multiplier field variable.
    INTEGER(INTG) :: numberOfColumns !<The number of columns in the interface mapping
    INTEGER(INTG) :: totalNumberOfColumns !<The total number of columns in the interface mapping
    INTEGER(INTG) :: numberOfGlobalColumns !<The global number of columns in the interface mapping
    TYPE(DomainMappingType), POINTER :: columnDOFSMapping !<A pointer to the domain mapping for the columns
    INTEGER(INTG), POINTER :: lagrangeDOFToColumnMap(:) !<lagrangeDOFToColumnMap(dofIdx). The mapping from the dofIdx'th Lagrange dof to the interface matrices column
    INTEGER(INTG) :: numberOfInterfaceMatrices !<The number of interface matrices that the mapping is set up for.
    TYPE(InterfaceMatrixToVarMapPtrType), ALLOCATABLE :: interfaceMatrixToVarMaps(:) !<interfaceMatrixToVarMaps(interfaceMatrixIdx). Information for the interface matrix to dependent variable maps for the interfaceMatrixIdx'th interface matrix.
    TYPE(InterfaceMappingRHSType), POINTER :: rhsMapping !<A pointer to the interface mapping for the RHS vector.
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache
  END TYPE InterfaceMappingType

  !>Contains information about the interpolation for a parameter set in interface equations
  TYPE InterfaceEquationsInterpolationSetType
    TYPE(FieldInterpolationParametersPtrType), POINTER :: interpolationParameters(:) !<interpolationParameters(fieldVariableType). A pointer to the field_variable_type'th field interpolation parameters.
    TYPE(FieldInterpolatedPointPtrType), POINTER :: interpolatedPoint(:) !<interpolatedPoint(fieldVariableType). A pointer to the field_variable_type'th field interpolated point. 
    TYPE(FieldInterpolatedPointMetricsPtrType), POINTER :: interpolatedPointMetrics(:) !<interpolatedPointMetrics(fieldVariableType). A pointer to the field_variable_type'th field interpolated point metrics.
  END TYPE InterfaceEquationsInterpolationSetType

  !>Contains information about the interpolation for a domain (interface or coupled mesh) in the interface equations
  TYPE InterfaceEquationsDomainInterpolationType
    TYPE(InterfaceEquationsInterpolationType), POINTER :: interpolation !<A pointer to the interpolation information used in the interface equations.
    TYPE(FieldType), POINTER :: geometricField !<A pointer to the geometric field for the domain
    INTEGER(INTG) :: numberOfGeometricInterpolationSets !<The number of geometric interpolation sets in the domain
    TYPE(InterfaceEquationsInterpolationSetType), ALLOCATABLE :: geometricInterpolation(:) !<geometricInterpolation(interpolationSetIdx). The geometric interpolation information for the interpolationSetIdx'th interpolation set.
    TYPE(FieldType), POINTER :: dependentField !<A pointer to the dependent field for the domain
    INTEGER(INTG) :: numberOfDependentInterpolationSets !<The number of dependent interpolation sets in the domain
    TYPE(InterfaceEquationsInterpolationSetType), ALLOCATABLE :: dependentInterpolation(:) !<dependentInterpolation(interpolationSetIdx). The dependent interpolation information for the interpolationSetIdx'th interpolation set.
    TYPE(FieldType), POINTER :: penaltyField !<A pointer to the penalty field for the domain
    INTEGER(INTG) :: numberOfPenaltyInterpolationSets !<The number of penalty interpolation sets in the domain
    TYPE(InterfaceEquationsInterpolationSetType), ALLOCATABLE :: penaltyInterpolation(:) !<penaltyInterpolation(interpolationSetIdx). The penalty interpolation information for the interpolationSetIdx'th interpolation set.
  END TYPE InterfaceEquationsDomainInterpolationType
  
  !>Contains information on the interpolation for the interface equations
  TYPE InterfaceEquationsInterpolationType
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the equations
    TYPE(InterfaceEquationsDomainInterpolationType) :: interfaceInterpolation !<The interpolation information for interpolating on the interface.
    TYPE(InterfaceEquationsDomainInterpolationType), ALLOCATABLE :: variableInterpolation(:) !<variableInterpolation(variableIdx). The interpolation for the variableIdx'th field variable in the interface condition.
  END TYPE InterfaceEquationsInterpolationType

  !>Contains information about the interface equations for an interface condition. 
  TYPE InterfaceEquationsType
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition
    LOGICAL :: interfaceEquationsFinished !<Is .TRUE. if the interface equations have finished being created, .FALSE. if not.
    INTEGER(INTG) :: outputType !<The output type for the interface equations \see InterfaceEquations_OutputTypes,InterfaceEquations
    INTEGER(INTG) :: sparsityType !<The sparsity type for the interface equation matrices of the interface equations \see InterfaceEquations_SparsityTypes,InterfaceEquations
    INTEGER(INTG) :: linearity !<The interface equations linearity type \see InterfaceEquations_LinearityTypes,InterfaceEquations
    INTEGER(INTG) :: timeDependence !<The interface equations time dependence type \see InterfaceEquations_TimeDependenceTypes,InterfaceEquations
    TYPE(InterfaceEquationsInterpolationType), POINTER :: interpolation !<A pointer to the interpolation information used in the interface equations.
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface equations mapping for the interface.
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<A pointer to the interface equations matrices and vectors used for the interface equations.
  END TYPE InterfaceEquationsType
  
  !>Contains information on the geometry for an interface condition
  TYPE InterfaceGeometryType
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition.
    TYPE(FieldType), POINTER :: geometricField !<The geometric field for this equations set.
  END TYPE InterfaceGeometryType

  !>Contains information about the penalty field information for an interface condition. 
  TYPE InterfacePenaltyType
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition
    LOGICAL :: penaltyFinished !<Is .TRUE. if the interface penalty field has finished being created, .FALSE. if not.
    LOGICAL :: penaltyFieldAutoCreated !<Is .TRUE. if the penalty field has been auto created, .FALSE. if not.
    TYPE(FieldType), POINTER :: penaltyField !<A pointer to the penalty field.
  END TYPE InterfacePenaltyType

  !>Contains information about the Lagrange field information for an interface condition. 
  TYPE InterfaceLagrangeType
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition
    LOGICAL :: lagrangeFinished !<Is .TRUE. if the interface Lagrange field has finished being created, .FALSE. if not.
    LOGICAL :: lagrangeFieldAutoCreated !<Is .TRUE. if the Lagrange field has been auto created, .FALSE. if not.
    TYPE(FieldType), POINTER :: lagrangeField !<A pointer to the lagrange field.
    INTEGER(INTG) :: numberOfComponents !<The number of components in the Lagrange field.
  END TYPE InterfaceLagrangeType

  !>Contains information about the dependent field information for an interface condition. 
  TYPE InterfaceDependentType
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition
    INTEGER(INTG) :: numberOfDependentVariables !<The number of dependent variables in the interface condition
    TYPE(EquationsSetPtrType), POINTER :: equationsSets(:) !<equationsSets(variableIdx). The pointer to the equations set containing the dependent variable for the variableIdx'th added dependent variable.
    TYPE(FieldVariablePtrType), POINTER :: fieldVariables(:) !<fieldVariables(variableIdx). The pointer to the variableIdx'th dependent variable in the interface condition.
    INTEGER(INTG), POINTER :: variableMeshIndices(:) !<variableMeshIndices(variableIdx). The mesh index of the variableIdx'th dependent variable in the interface condition.
  END TYPE InterfaceDependentType

  !>Contains information for the interface condition data.
  TYPE InterfaceConditionType
    INTEGER(INTG) :: userNumber !<The user identifying number of the interface condition. Must be unique.
    INTEGER(INTG) :: globalNumber !<The global index of the interface condition in the interface conditions.
    LOGICAL :: interfaceConditionFinished !<Is .TRUE. ifand where  the interfaand where and where ce condition has finished being created, .FALSE. if not.
    TYPE(InterfaceConditionsType), POINTER :: interfaceConditions !<A pointer back to the interface conditions.
    TYPE(VARYING_STRING) :: label !<A user defined label for the interface condition.
    TYPE(InterfaceType), POINTER :: INTERFACE !<A pointer back to the interface.
    INTEGER(INTG) :: outputType !<The output type for the interface \see InterfaceCondition_OutputTypes,InterfaceCondition
    INTEGER(INTG) :: method !<An integer which denotes the interface condition method. \see InterfaceCondition_Methods,InterfaceCondition
    INTEGER(INTG) :: operator !<An integer which denotes the type of interface operator. \see InterfaceCondition_Operator,InterfaceCondition
    INTEGER(INTG) :: integrationType !<An integer which denotes the integration type. \see InterfaceCondition_IntegrationType,InterfaceCondition
    TYPE(InterfaceGeometryType) :: geometry !<The geometry information for the interface condition.
    TYPE(InterfacePenaltyType), POINTER :: penalty !<A pointer to the interface condition penalty information if there are any for this interface condition.
    TYPE(InterfaceLagrangeType), POINTER :: lagrange !<A pointer to the interface condition Lagrange multipler information if there are any for this interface condition.
    TYPE(InterfaceDependentType), POINTER :: dependent !<A pointer to the interface condition dependent field information if there is any for this interface condition.
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations if there are any for this interface condition.
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary condition information for this interface condition.
  END TYPE InterfaceConditionType

  !>A buffer type to allow for an array of pointers to a InterfaceConditionType.
  TYPE InterfaceConditionPtrType
    TYPE(InterfaceConditionType), POINTER :: ptr !<The pointer to the interface condition.
  END TYPE InterfaceConditionPtrType
  
  !>Contains information for interface region specific data that is not of 'region' importance. <<>>
  TYPE InterfaceConditionsType
    TYPE(InterfaceType), POINTER :: INTERFACE !<A pointer back to the interface containing the interface conditions.
    INTEGER(INTG) :: numberOfInterfaceConditions !<The number of interface conditions
    TYPE(InterfaceConditionPtrType), POINTER :: interfaceConditions(:) !<interfaceConditions(interfaceConditionIdx). A pointer to the interfaceConditionIdx'th interface condition.
  END TYPE InterfaceConditionsType
  
  !>Contains information on the mesh connectivity for a given coupled mesh element
  TYPE InterfaceElementConnectivityType
    INTEGER(INTG) :: coupledElementNumber !< The coupled number to define the connectivity for.
    REAL(DP), ALLOCATABLE :: xi(:,:,:) !<xi(xiIdx,meshComponent,elementParameterIdx) The xiIdx'th xi of a coupled element to be copuled to an interface mesh's interfaceMeshComponentIdx'th interface mesh components's elementParameterIdx'th element parameter. !\todo the XI array needs to be restructured to efficiently utilize memory when coupling bodies with 2xi directions to bodies with 3xi directions using an interface.
    INTEGER(INTG) :: connectedLineFace !<The coupled  element line/face number to be connected to the interface mesh.
  END TYPE InterfaceElementConnectivityType

  !>A buffer type to allow arrays of pointers to InterfaceElementConnectivityType
  TYPE InterfaceElementConnectivityPtrType
    TYPE(InterfaceElementConnectivityType), POINTER :: ptr !<A pointer to the InterfaceElementConnectivityType
  END TYPE InterfaceElementConnectivityPtrType

  !>Contains information on the coupling between meshes in an interface
  TYPE InterfaceMeshConnectivityType
    TYPE(InterfaceType), POINTER :: interface !<A pointer back to the interface for the coupled mesh connectivity
    TYPE(MeshType), POINTER :: interfaceMesh !<A pointer to the inteface mesh
    TYPE(BasisType), POINTER :: basis !<A pointer to the inteface mesh basis
    LOGICAL :: meshConnectivityFinished !<Is .TRUE. if the coupled mesh connectivity has finished being created, .FALSE. if not.
    INTEGER(INTG) :: numberOfInterfaceElements !<The number of elements in the interface
    INTEGER(INTG) :: numberOfInterfaceNodes !<The number of nodes in the interface
    INTEGER(INTG) :: numberOfCoupledMeshes !<The number of coupled meshes in the interface
    TYPE(InterfaceElementConnectivityPtrType), ALLOCATABLE :: elementConnectivity(:,:) !<elementConnectivity(elementIdx,coupledMeshIdx) !<The mesh connectivity for a given interface mesh element
    INTEGER(INTG), ALLOCATABLE :: coupledNodes(:,:) !<coupledNodes(coupledMeshIdx,interfaceNodeIdx). Coupled nodes numbers
  END TYPE InterfaceMeshConnectivityType

  !>Contains information on the coupling between decompositions in an interface
  TYPE InterfaceDecompositionConnectivityType
    TYPE(InterfaceType), POINTER :: interface !<A pointer back to the interface for the coupled mesh connectivity
    TYPE(DecompositionType), POINTER :: interfaceDecomposition !<A pointer to the inteface decomposition
    TYPE(BasisType), POINTER :: basis !<A pointer to the inteface mesh basis
    LOGICAL :: decompositionConnectivityFinished !<Is .TRUE. if the coupled decomposition connectivity has finished being created, .FALSE. if not.
    INTEGER(INTG) :: numberOfInterfaceElements !<The number of elements in the interface
    INTEGER(INTG) :: totalNumberOfInterfaceElements !<The total number (including ghosts) of elements in the interface
    INTEGER(INTG) :: numberOfGlobalInterfaceElements !<The number of global interface elements in the interface
    INTEGER(INTG) :: numberOfCoupledDecompositions !<The number of coupled decompositions in the interface
    TYPE(InterfaceElementConnectivityPtrType), ALLOCATABLE :: elementConnectivity(:,:) !<elementConnectivity(elementIdx,coupledMeshIdx) !<The mesh connectivity for a given interface mesh element
    INTEGER(INTG), ALLOCATABLE :: coupledNodes(:,:) !<coupledNodes(coupledMeshIdx,interfaceNodeIdx). Coupled nodes numbers
  END TYPE InterfaceDecompositionConnectivityType

  !>Contains information on a data connectivity point 
  TYPE InterfacePointConnectivityType
    INTEGER(INTG) :: coupledElementNumber !<The element number this point is connected to in the coupled mesh
    INTEGER(INTG) :: elementLineFaceNumber !<The local connected face/line number in the coupled mesh
    REAL(DP), ALLOCATABLE :: xi(:) !<xi(xiIdx). The full xi location the data point is connected to in this coupled mesh
    REAL(DP), ALLOCATABLE :: reducedXi(:) !<reducedXi(xiIdx). The reduced (face/line) xi location the data point is connected to in this coupled mesh
  END TYPE InterfacePointConnectivityType

  !>A buffer type for arrays of pointers to InterfacePointConnectivityType
  TYPE InterfacePointConnectivityPtrType
    TYPE(InterfacePointConnectivityType), POINTER :: ptr !<A pointer to the InterfacePointConnectivityType
  END TYPE InterfacePointConnectivityPtrType
  
  !Contains information on coupled mesh elements that are connected to each interface element.
  TYPE InterfaceCoupledElementsType
    INTEGER(INTG) :: numberOfCoupledElements
    INTEGER(INTG), ALLOCATABLE :: elementNumbers(:) !<elementNumbers(elementIdx). The global/local(if updated after decomposition) numbers of the coupled mesh elements that are connected to this interface element.
  END TYPE InterfaceCoupledElementsType
  
  !>Contains information on the data point coupling/points connectivity between meshes in the an interface
  TYPE InterfacePointsConnectivityType
    TYPE(InterfaceType), POINTER :: INTERFACE !<A pointer back to the interface for the coupled mesh connectivity
    TYPE(MeshType), POINTER :: interfaceMesh !<A pointer to the interface mesh where the xi locations of data points are defined
    LOGICAL :: pointsConnectivityFinished !<Is .TRUE. if the data points connectivity has finished being created, .FALSE. if not.
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points defined on the interface for the connectivity
    TYPE(InterfacePointConnectivityPtrType), ALLOCATABLE :: pointsConnectivity(:,:) !<pointsConnectivity(dataPointIndex,coupledMeshIdx). A pointer to the points connectivity information for each data point in each coupled mesh. 
    TYPE(InterfaceCoupledElementsType), ALLOCATABLE :: coupledElements(:,:) !<coupledElements(interfaceElementIdx,coupledMeshIdx). The coupled mesh elements that are connected to each interface element.
    INTEGER(INTG), ALLOCATABLE :: maxNumberOfCoupledElements(:) !<maxNumberOfCoupledElements(coupledMeshIdx). The maximum number of coupled elements to an interface element in coupledMeshIdx'th mesh
  END TYPE InterfacePointsConnectivityType
 
  !>Contains information for the interface data.
  TYPE InterfaceType
    INTEGER :: userNumber !<The user defined identifier for the interface. The user number must be unique.
    INTEGER :: globalNumber !<The global number of the interface in the list of interfaces for a particular parent region.
    LOGICAL :: interfaceFinished !<Is .TRUE. if the interface has finished being created, .FALSE. if not.
    TYPE(VARYING_STRING) :: label !<A user defined label for the region.
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<A pointer to the coordinate system used by the interface. 
    TYPE(InterfacesType), POINTER :: interfaces !<A pointer back to the parent interfaces
    TYPE(RegionType), POINTER :: parentRegion !<A point to the parent region containing the interface.
    INTEGER(INTG) :: numberOfCoupledMeshes !<The number of coupled meshes in the interface.
    TYPE(MeshPtrType), ALLOCATABLE :: coupledMeshes(:) !<coupledMeshes(meshIdx). coupledMeshes(meshIdx)%ptr is the pointer to the meshIdx'th mesh involved in the interface.
    TYPE(InterfaceDecompositionConnectivityType), POINTER :: decompositionConnectivity !<A pointer to the decomposition connectivity the interface.
    TYPE(InterfaceMeshConnectivityType), POINTER :: meshConnectivity !<A pointer to the meshes connectivity the interface.
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity !<A pointer to the points connectivity the interface.
    TYPE(DataPointSetsType), POINTER :: dataPointSets  !<A pointer to the data points defined in an interface.
    TYPE(NodesType), POINTER :: nodes !<A pointer to the nodes in an interface
    TYPE(MeshesType), POINTER :: meshes !<A pointer to the mesh in an interface.
    TYPE(GeneratedMeshesType), POINTER :: generatedMeshes !<A pointer to the generated meshes in an interface.
    TYPE(FieldsType), POINTER :: fields !<A pointer to the fields defined over an interface.
    TYPE(InterfaceConditionsType), POINTER :: interfaceConditions !<The pointer to the interface conditions for this interface
  END TYPE InterfaceType

  !>A buffer type to allow for an array of pointers to a InterfaceType.
  TYPE InterfacePtrType
    TYPE(InterfaceType), POINTER :: ptr !<The pointer to the interface.
  END TYPE InterfacePtrType
  
  !>Contains information for interfaces on a parent region.
  TYPE InterfacesType
    TYPE(RegionType), POINTER :: parentRegion !<A pointer back to the parent region containing the interfaces.
    INTEGER(INTG) :: numberOfInterfaces !<The number of interfaces
    TYPE(InterfacePtrType), ALLOCATABLE :: interfaces(:) !<interfaces(interfaceIdx). A pointer to the interfaceIdx'th interface.
  END TYPE InterfacesType

  PUBLIC InterfaceMatrixType,InterfaceMatrixPtrType

  PUBLIC InterfaceMatrixToVarMapType

  PUBLIC InterfaceMatricesType

  PUBLIC InterfaceRHSType

  PUBLIC InterfaceMappingRHSType

  PUBLIC InterfaceMappingCreateValuesCacheType

  PUBLIC InterfaceMappingType

  PUBLIC InterfaceEquationsInterpolationSetType

  PUBLIC InterfaceEquationsDomainInterpolationType

  PUBLIC InterfaceEquationsInterpolationType

  PUBLIC InterfaceEquationsType

  PUBLIC InterfaceGeometryType

  PUBLIC InterfacePenaltyType

  PUBLIC InterfaceLagrangeType

  PUBLIC InterfaceDependentType

  PUBLIC InterfaceConditionType,InterfaceConditionPtrType,InterfaceConditionsType

  PUBLIC InterfaceDecompositionConnectivityType
  
  PUBLIC InterfaceElementConnectivityType,InterfaceElementConnectivityPtrType

  PUBLIC InterfaceMeshConnectivityType

  PUBLIC InterfacePointConnectivityType,InterfacePointsConnectivityType,InterfacePointConnectivityPtrType

  PUBLIC InterfaceCoupledElementsType

  PUBLIC InterfaceType,InterfacePtrType,InterfacesType
  
  !
  !================================================================================================================================
  !
  ! CellML types

  !>This type is a wrapper for the C_PTR which references the actual CellML model definition object.
  TYPE CellMLModelType
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML environment.
    INTEGER(INTG) :: globalNumber !< The global number of this CellML model within the parent CellML environment.
    TYPE(VARYING_STRING) :: modelId !<The ID of the model.
    TYPE(C_PTR) :: ptr !< The handle for the actual C++ CellML model definition object
    INTEGER(INTG) :: numberOfState !<The number of state variables in the CellML model.
    TYPE(VARYING_STRING), ALLOCATABLE :: stateVariableId(:) !<stateVariableId(stateVariableIdx). The ID for the stateVariableIdx'th state variable.
    INTEGER(INTG) :: numberOfIntermediate !<The number of intermediate variables in the CellML model.
    TYPE(VARYING_STRING), ALLOCATABLE :: intermediateVariableId(:) !<intermediateVariableId(intermediateVariableIdx). The ID for the intermediateVariableIdx'th intermediate variable.
    INTEGER(INTG) :: numberOfParameters !<The number of parameters in the CellML model.
    TYPE(VARYING_STRING), ALLOCATABLE :: parameterVariableId(:) !<parameterVariableId(parameterVariableIdx). The ID for the parameter_variable_idx'th parameter variable.
 END TYPE CellMLModelType

  !>A buffer type to allow for an array of pointers to a CellMLModelType
  TYPE CellMLModelPtrType
    TYPE(CellMLModelType), POINTER :: ptr
  END TYPE CellMLModelPtrType

  !>Contains information on the models field for a CellML environment
  TYPE CellMLModelsFieldType
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML environment 
    LOGICAL :: modelsFieldFinished  !<Is .TRUE. if the models field has finished being created, .FALSE. if not.
    LOGICAL :: modelsFieldAutoCreated !<Is .TRUE. if the models field has been auto created, .FALSE. if not.
    TYPE(FieldType), POINTER :: modelsField !<A pointer to the models field
    INTEGER(INTG) :: onlyOneModelIndex !<If only one model is used in the models field for the CellML environment then this will be equal to the model index. It will be zero otherwise.
  END TYPE CellMLModelsFieldType
  
  !>Contains information on the state field for a CellML environment
  TYPE CellMLStateFieldType
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML environment 
    LOGICAL :: stateFieldFinished  !<Is .TRUE. if the state field has finished being created, .FALSE. if not.
    LOGICAL :: stateFieldAutoCreated !<Is .TRUE. if the state field has been auto created, .FALSE. if not.
    TYPE(FieldType), POINTER :: stateField !<A pointer to the state field
  END TYPE CellMLStateFieldType
  
  !>Contains information on the intermediate field for a CellML environment
  TYPE CellMLIntermediateFieldType
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML environment 
    LOGICAL :: intermediateFieldFinished  !<Is .TRUE. if the intermediate field has finished being created, .FALSE. if not.
    LOGICAL :: intermediateFieldAutoCreated!<Is .TRUE. if the intermediate field has been auto created, .FALSE. if not.
    TYPE(FieldType), POINTER :: intermediateField !<A pointer to the intermediate field
  END TYPE CellMLIntermediateFieldType
  
  !>Contains information on the parameters field for a CellML environment
  TYPE CellMLParametersFieldType
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML environment 
    LOGICAL :: parametersFieldFinished  !<Is .TRUE. if the parameters field has finished being created, .FALSE. if not.
    LOGICAL :: parametersFieldAutoCreated !<Is .TRUE. if the parameters field has been auto created, .FALSE. if not.
    TYPE(FieldType), POINTER :: parametersField !<A pointer to the parameters field
  END TYPE CellMLParametersFieldType
 
  !> Contains information on the solver, cellml, dof etc. for which cellml equations are to be evaluated by petsc
  TYPE CellMLPETScContextType
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    TYPE(CellMLType), POINTER :: cellml !<A pointer to the CellML environment
    INTEGER(INTG) :: dofIdx !<The current DOF to be evaluated
    REAL(DP), POINTER :: rates(:) !<A pointer to the temporary rates array
    INTEGER(INTG), ALLOCATABLE :: ratesIndices(:) !<The PETSc array indices for the rates
  END TYPE CellMLPETScContextType
  
  !>Contains information on the mapping between CellML fields and OpenCMISS fields and vise versa.
  TYPE CellMLModelMapType
    INTEGER(INTG) :: cellMLMapType !<The direction of the mapping. \see CellML_FieldMappingTypes,CmissCellML
    TYPE(FieldType), POINTER :: field !<A pointer to/from the field being mapped.
    INTEGER(INTG) :: variableType !<The field variable type being mapped.
    INTEGER(INTG) :: componentNumber !<The field variable component number being mapped.
    INTEGER(INTG) :: fieldParameterSet !<The field variable parameter set being mapped.
    TYPE(VARYING_STRING) :: variableId !<The variable ID of the CellML variable being mapped.
    INTEGER(INTG) :: cellMLFieldType !<The type of CellML field variable being mapped. \see CellML_FieldTypes,CmissCellML
    INTEGER(INTG) :: cellMLVariableNumber !<The CellML variable component number being mapped.
    INTEGER(INTG) :: cellMLParameterSet !<The CellML variable parameter set being mapped.
  END TYPE CellMLModelMapType

  !>Buffer type to allow an array of pointers to CellMLModelMapType
  TYPE CellMLModelMapPtrType
    TYPE(CellMLModelMapType), POINTER :: ptr !<A pointer to the CellMLModelMapType
  END TYPE CellMLModelMapPtrType

  !>Contains information on the maps between a CellML model and external OpenCMISS fields.
  TYPE CellMLModelMapsType
    TYPE(CellMLFieldMapsType), POINTER :: cellMLFieldMaps !<A pointer back to the field maps
    INTEGER(INTG) :: numberOfFieldsMappedTo !<The number of OpenCMISS fields mapped to this CellML model's variables.
    TYPE(CellMLModelMapPtrType), ALLOCATABLE :: fieldsMappedTo(:) !<fieldsMappedTo(mapIdx). The mapIdx'th field mapping for OpenCMISS fields that the CellML model maps to.
    INTEGER(INTG) :: numberOfFieldsMappedFrom !<The number of CellML variable fields mapped to OpenCMISS fields.
    TYPE(CellMLModelMapPtrType), ALLOCATABLE :: fieldsMappedFrom(:) !<fieldsMappedFrom(mapIdx). The mapIdx'th field mapping for OpenCMISS fields that the CellML model maps from.
  END TYPE CellMLModelMapsType

  !>Buffer type to allow arrays of pointer to CellMLModelMapsType
  TYPE CellMLModelMapsPtrType
    TYPE(CellMLModelMapsType), POINTER :: ptr !<A pointer to the CellMLModelMapsType
  END TYPE CellMLModelMapsPtrType
  
  !>Contains information on the maps between CellML and external OpenCMISS fields.
  TYPE CellMLFieldMapsType
    TYPE(CellMLType), POINTER :: cellML !<A pointer back to the CellML environment.
    LOGICAL :: cellMLFieldMapsFinished !<Is .TRUE. if the CellML maps have finished being created, .FALSE. if not.
    TYPE(FieldType), POINTER :: sourceGeometricField !<The source geometric field for the CellML environment.
    TYPE(FieldVariableType), POINTER :: sourceFieldVariable !<The source field variable for the CellML environment.
    TYPE(DomainType), POINTER :: sourceFieldDomain !<The source field domain for the CellML environment.
    INTEGER(INTG) :: sourceFieldInterpolationType !<The source field interpolation type for the CellML environment.
    TYPE(CellMLModelMapsPtrType), ALLOCATABLE :: modelMaps(:) !<modelMaps(modelIdx). Contains information on the maps between the model_idx'th CellML model and external OpenCMISS fields.
  END TYPE CellMLFieldMapsType
  
  !>Contains information for a CellML environment.
  TYPE CellMLType
    TYPE(RegionType), POINTER :: region !<A pointer to the region containing this CellML environment.
    INTEGER(INTG) :: globalNumber !<The global number of the CellML environment in the list of environments for a field.
    INTEGER(INTG) :: userNumber !<The user defined identifier for the CellML environment. The user number must be unique.
    TYPE(CellMLEnvironmentsType), POINTER :: environments !<A pointer back to the CellML environments.
    LOGICAL :: cellMLFinished !<Is .TRUE. if the environment has finished being created, .FALSE. if not.
    INTEGER(INTG) :: numberOfModels !< The number of models defined in the CellML environment
    INTEGER(INTG) :: maximumNumberOfState !<The maximum number of state variables across all models.
    INTEGER(INTG) :: maximumNumberOfParameters !<The maximum number of parameters variables across all models.
    INTEGER(INTG) :: maximumNumberOfIntermediate !<The maximum number of intermediate variables across all models.
    TYPE(CellMLModelPtrType), ALLOCATABLE :: models(:) !< models(modelIdx). The array of pointers to the models.
    TYPE(CellMLFieldMapsType), POINTER :: fieldMaps !<A pointer the information on CellML<-->Field maps.
    TYPE(CellMLModelsFieldType), POINTER :: modelsField !<A pointer to the models field information.
    TYPE(CellMLStateFieldType), POINTER :: stateField !<A pointer to the state field information.
    TYPE(CellMLIntermediateFieldType), POINTER :: intermediateField !<A pointer to the intermediate field information.
    TYPE(CellMLParametersFieldType), POINTER :: parametersField !<A pointer to the parameters field information.
    LOGICAL :: cellMLGenerated !<Is .TRUE. if the CellML environment has finished being generated, .FALSE. if not.
  END TYPE CellMLType

  !> A buffer type to allow for an array of pointers to a CellMLType.
  TYPE CellMLPtrType
    TYPE(CellMLType), POINTER :: ptr !< The pointer to the CellML environment.
  END TYPE CellMLPtrType

  !>Contains information on the CellML environments defined.
  TYPE CellMLEnvironmentsType
    TYPE(RegionType), POINTER :: region !<A pointer to the region containing the CellML environments
    INTEGER(INTG) :: numberOfEnvironments !<The number of environments defined.
    TYPE(CellMLPtrType), ALLOCATABLE :: environments(:) !<The array of pointers to the CellML environments.
  END TYPE CellMLEnvironmentsType

  PUBLIC CellMLModelType,CellMLModelPtrType

  PUBLIC CellMLModelsFieldType

  PUBLIC CellMLStateFieldType

  PUBLIC CellMLIntermediateFieldType

  PUBLIC CellMLParametersFieldType

  PUBLIC CellMLPETScContextType

  PUBLIC CellMLModelMapType,CellMLModelMapPtrType,CellMLModelMapsType,CellMLModelMapsPtrType

  PUBLIC CellMLFieldMapsType

  PUBLIC CellMLType,CellMLPtrType

  PUBLIC CellMLEnvironmentsType

  !
  !================================================================================================================================
  !
  ! Solver matrices types
  !
  
  !>Contains information on the solver matrix
  TYPE SolverMatrixType
    INTEGER(INTG) :: matrixNumber !<The number of the solver matrix
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices for this solver matrix
    LOGICAL :: updateMatrix !<Is .TRUE. if the solver matrix is to be updated
    INTEGER(INTG) :: storageType !<The storage type for the solver matrix.
    INTEGER(INTG) :: symmetryType !<The solver matrix symmetry type.
    INTEGER(INTG) :: numberOfColumns !<The number of columns in the distributed solver matrix
    TYPE(DistributedVectorType), POINTER :: solverVector !<A pointer to the distributed solver vector associated with the matrix
    TYPE(DistributedMatrixType), POINTER :: matrix !<A pointer to the distributed solver matrix data
  END TYPE SolverMatrixType

  !>A buffer type to allow for an array of pointers to a SolverMatrixType \see TYPES:SOLUTION_MATRIX_TYPE
  TYPE SolverMatrixPtrType
    TYPE(SolverMatrixType), POINTER :: ptr !<The pointer to the solver matrix.
  END TYPE SolverMatrixPtrType
  
  !>Contains information on the solver matrices and rhs vector
  TYPE SolverMatricesType
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver
    LOGICAL :: solverMatricesFinished !<Is .TRUE. if the solver matrices have finished being created, .FALSE. if not.
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping for these solver matrices
    INTEGER(INTG) :: numberOfRows !<The number of (local) rows in the distributed solution matrix for this computation node
    INTEGER(INTG) :: numberOfGlobalRows !<The number of global rows in the distributed solution matrix
    INTEGER(INTG) :: solverLibraryType !<The library type of the solver for the solver matrices \see SolverRoutines_SolverLibraries
    INTEGER(INTG) :: matrixLibraryType !<The library type for the solver matrices \see DistributedMatrixVector_LibraryTypes
    !Linear matrices
    INTEGER(INTG) :: numberOfMatrices !<The number of solver matrices defined for the problem
    TYPE(SolverMatrixPtrType), ALLOCATABLE :: matrices(:) !<matrices(matrixIdx)%ptr contains the information on the matrixIdx'th solver matrix
    !Nonlinear matrices and vectors
    LOGICAL :: updateResidual !<Is .TRUE. if the residual vector is to be updated
    TYPE(DistributedVectorType), POINTER :: residual !<A pointer to the distributed residual vector for nonlinear problems
    !Optimiser matrices and vectors
    LOGICAL :: updateGradient !<Is .TRUE. if the gradient vector is to be updated
    TYPE(DistributedVectorType), POINTER :: gradient !<A pointer to the distributed gradient vector for optimisation problems
    !Right hand side vector
    LOGICAL :: updateRHSVector !<Is .TRUE. if the RHS vector is to be updated
    TYPE(DistributedVectorType), POINTER :: rhsVector !<A pointer to the distributed RHS vector for the solver matrices
  END TYPE SolverMatricesType

  PUBLIC SolverMatrixType,SolverMatrixPtrType

  PUBLIC SolverMatricesType

  !
  !================================================================================================================================
  !
  ! Solver equations types
  !

  !>Contains information about the solver equations for a solver. \see OpenCMISS::Iron::cmfe_SolverEquationsType
  TYPE SolverEquationsType
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    LOGICAL :: solverEquationsFinished !<Is .TRUE. if the solver equations have finished being created, .FALSE. if not.
    INTEGER(INTG) :: linearity !<The linearity type of the solver equations
    INTEGER(INTG) :: timeDependence !<The time dependence type of the solver equations
    INTEGER(INTG) :: sparsityType !<The type of sparsity to use in the solver matrices \see SolverRoutines_EquationsSparsityTypes,SolverRoutines
    INTEGER(INTG) :: symmetryType !<The type of symmetry to use in the solver matrices \see SolverRoutines_SymmetryTypes,SolverRoutines
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping
    TYPE(SolverMatricesType), POINTER :: solverMatrices !<A pointer to the solver matrices for the problem
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary condition information for the solver equations.
  END TYPE SolverEquationsType

  PUBLIC SolverEquationsType

  !
  !================================================================================================================================
  !
  ! CellML equations types
  !

  !>Contains information about the CellML equations for a solver.
  TYPE CellMLEquationsType
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    LOGICAL :: cellMLEquationsFinished !<Is .TRUE. if the CellML equations have finished being created, .FALSE. if not.
    INTEGER(INTG) :: linearity !<The linearity type of the CellML equations
    INTEGER(INTG) :: timeDependence !<The time dependence type of the CellML equations
    REAL(DP) :: currentTime !<The current time to evaluate the equations at
    INTEGER(INTG) :: numberOfCellMLEnvironments !<The number of CellML environments in the equations
    TYPE(CellMLPtrType), ALLOCATABLE :: cellMLEnvironments(:) !<cellMLEnvironments(cellmlIdx). The array of pointers to the CellML environments for these CellML equations.     
  END TYPE CellMLEquationsType

  PUBLIC CellMLEquationsType
  
  !
  !================================================================================================================================
  !
  ! Solver types
  !
  
  !>Contains information for a dynamic solver
  TYPE DynamicSolverType
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG) :: solverLibrary !<The library type for the dynamic solver \see SolverRoutines_SolverLibraries,SolverRoutines
    LOGICAL :: solverInitialised !<Is .TRUE. if the solver has been initialised, .FALSE. if not.
    INTEGER(INTG) :: numberOfSolves !<The number of dynamic solves for this solver. 
    INTEGER(INTG) :: linearity !<The linearity type of the dynamic solver \see SolverRoutines_DynamicLinearityTypes,SolverRoutines
    INTEGER(INTG) :: order !<The order of the dynamic solve \see SolverRoutines_DynamicOrderTypes,SolverRoutines
    INTEGER(INTG) :: degree !<The degree of the time interpolation polynomial \see SolverRoutines_DynamicDegreeTypes,SolverRoutines
    INTEGER(INTG) :: scheme !<The dyanamic solver scheme \see SolverRoutines_DynamicSchemeTypes,SolverRoutines
    INTEGER(INTG) :: startupType !<The previous values startup type i.e., what prevous values to use during startup \see SolverRoutines_DynamicStartupTypes,SolverRoutines
    REAL(DP), ALLOCATABLE :: theta(:) !<theta(degree_idx). The theta value for the degree_idx'th polynomial in the dynamic solver
    LOGICAL :: explicit !<Is .TRUE. if the dynamic scheme is an explicit scheme, .FALSE. if not.
    LOGICAL :: restart !<Is .TRUE. if the dynamic scheme is to be restarted (i.e., recalculate values at the current time step), .FALSE. if not.    
    LOGICAL :: ale !<Is .TRUE. if the dynamic scheme is an ALE scheme, .FALSE. if not.
    REAL(DP) :: currentTime !<The current time value for the dynamic solver.
    REAL(DP) :: timeIncrement !<The time increment for the dynamic solver to solver for.
    TYPE(SolverType), POINTER :: linearSolver !<A pointer to the linked linear solver
    TYPE(SolverType), POINTER :: nonlinearSolver !<A pointer to the linked nonlinear solver
  END TYPE DynamicSolverType
  
   !>Contains information for an forward Euler differential-algebraic equation solver
  TYPE ForwardEulerDAESolverType
    TYPE(EulerDAESolverType), POINTER :: eulerDAESolver !<A pointer to the differential-algebraic solver
    INTEGER(INTG) :: solverLibrary !<The library type for the forward Euler differential-algebraic equation solver \see SolverRoutines_SolverLibraries,SolverRoutines
  END TYPE ForwardEulerDAESolverType

  !>Contains information for an backward Euler differential-algebraic equation solver
  TYPE BackwardEulerDAESolverType
    TYPE(EulerDAESolverType), POINTER :: eulerDAESolver !<A pointer to the differential-algebraic solver
    INTEGER(INTG) :: solverLibrary !<The library type for the backward Euler differential-algebraic equation solver \see SolverRoutines_SolverLibraries,SolverRoutines
  END TYPE BackwardEulerDAESolverType

  !>Contains information for an improved Euler differential-algebraic equation solver
  TYPE ImprovedEulerDAESolverType
    TYPE(EulerDAESolverType), POINTER :: eulerDAESolver !<A pointer to the differential-algebraic solver
    INTEGER(INTG) :: solverLibrary !<The library type for the improved Euler differential-algebraic equation solver \see SolverRoutines_SolverLibraries,SolverRoutines
  END TYPE ImprovedEulerDAESolverType
  
  !>Contains information for an Euler differential-algebraic equation solver
  TYPE EulerDAESolverType
    TYPE(DAESolverType), POINTER :: DAESolver !<A pointer to the differential-algebraic solver
    INTEGER(INTG) :: eulerType !<The type of Euler differential-algebraic equation solver \see SolverRoutines_EulerDESolverTypes,SolverRoutines
    TYPE(ForwardEulerDAESolverType), POINTER :: forwardEulerSolver !<A pointer to the forward Euler solver information
    TYPE(BackwardEulerDAESolverType), POINTER :: backwardEulerSolver !<A pointer to the backward Euler solver information
    TYPE(ImprovedEulerDAESolverType), POINTER :: improvedEulerSolver !<A pointer to the improved Euler solver information
    INTEGER(INTG) :: solverLibrary !<The library type for the Euler differential equation solver \see SolverRoutines_SolverLibraries,SolverRoutines
  END TYPE EulerDAESolverType

  !>Contains information for a Crank-Nicholson differential-algebraic equation solver
  TYPE CrankNicolsonDAESolverType
    TYPE(DAESolverType), POINTER :: DAESolver !<A pointer to the differential-algebraic solver
    INTEGER(INTG) :: solverLibrary !<The library type for the Crank-Nicholson differential-algebraic equation solver \see SolverRoutines_SolverLibraries,SolverRoutines
  END TYPE CrankNicolsonDAESolverType
  
  !>Contains information for a Runge-Kutta differential-algebraic equation solver
  TYPE RungeKuttaDAESolverType
    TYPE(DAESolverType), POINTER :: DAESolver !<A pointer to the differential-algebraic solver
    INTEGER(INTG) :: solverLibrary !<The library type for the Runge-Kutta differential-algebraic equation solver \see SolverRoutines_SolverLibraries,SolverRoutines
  END TYPE RungeKuttaDAESolverType
  
  !>Contains information for an Adams-Moulton differential-algebraic equation solver
  TYPE AdamsMoultonDAESolverType
    TYPE(DAESolverType), POINTER :: DAESolver !<A pointer to the differential-algebraic solver
    INTEGER(INTG) :: solverLibrary !<The library type for the Adams-Moulton differential-algebraic equation solver \see SolverRoutines_SolverLibraries,SolverRoutines
  END TYPE AdamsMoultonDAESolverType
  
  !>Contains information for a BDF differential-algebraic equation solver
  TYPE BDFDAESolverType
    TYPE(DAESolverType), POINTER :: DAESolver !<A pointer to the differential-algebraic solver
    INTEGER(INTG) :: solverLibrary !<The library type for the BDF differential-algebraic equation solver \see SolverRoutines_SolverLibraries,SolverRoutines
  END TYPE BDFDAESolverType
  
  !>Contains information for a Rush-Larson differential-algebraic equation solver
  TYPE RushLarsonDAESolverType
    TYPE(DAESolverType), POINTER :: DAESolver !<A pointer to the differential-algebraic solver
    INTEGER(INTG) :: solverLibrary !<The library type for the Rush-Larson differential-algebraic equation solver \see SolverRoutines_SolverLibraries,SolverRoutines
  END TYPE RushLarsonDAESolverType
  
  !>Contains information for an external differential-algebraic equation solver
  TYPE ExternalDAESolverType
    TYPE(DAESolverType), POINTER :: DAESolver !<A pointer to the differential-algebraic solver
  END TYPE ExternalDAESolverType
  
  !>Contains information for an differential-algebraic equation solver
  TYPE DAESolverType
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG) :: daeType !<The differential-algebraic equation type \see SolverRoutines_DAETypes,SolverRoutines
    INTEGER(INTG) :: daeSolveType !<The solve type for the differential-algebraic equation solver \see SolverRoutines_DAESolveTypes,SolverRoutines
    REAL(DP) :: startTime !<The start time to integrate from
    REAL(DP) :: endTime !<The end time to integrate to
    REAL(DP) :: initialStep !<The (initial) time step
    TYPE(EulerDAESolverType), POINTER :: eulerSolver !<A pointer to information for an Euler solver
    TYPE(CrankNicolsonDAESolverType), POINTER :: crankNicolsonSolver !<A pointer to information for a Crank-Nicholson solver
    TYPE(RungeKuttaDAESolverType), POINTER :: rungeKuttaSolver !<A pointer to information for a Runge-Kutta solver
    TYPE(AdamsMoultonDAESolverType), POINTER :: adamsMoultonSolver !<A pointer to information for an Adams-Moulton solver
    TYPE(BDFDAESolverType), POINTER :: bdfSolver !<A pointer to information for a BDF solver
    TYPE(RushLarsonDAESolverType), POINTER :: rushLarsonSolver !<A pointer to information for a Rush-Larson solver
    TYPE(ExternalDAESolverType), POINTER :: externalSolver !<A pointer to information for an external solver
  END TYPE DAESolverType
  
  !>Contains information for a direct linear solver
  TYPE LinearDirectSolverType
    TYPE(LinearSolverType), POINTER :: linearSolver !<A pointer to the linear solver
    INTEGER(INTG) :: solverLibrary !<The library type for the linear direct solver \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG) :: directSolverType !<The type of direct linear solver
    TYPE(PetscPCType) :: pc !<The PETSc preconditioner object
    TYPE(PetscKspType) :: ksp !<The PETSc solver object
  END TYPE LinearDirectSolverType

  !>Contains information for an iterative linear solver
  TYPE LinearIterativeSolverType
    TYPE(LinearSolverType), POINTER :: linearSolver !<A pointer to the linear solver
    INTEGER(INTG) :: solverLibrary !<The library type for the linear solver \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG) :: iterativeSolverType !<The type of iterative solver \see SolverRoutines_IterativeLinearSolverTypes,SolverRoutines    
    INTEGER(INTG) :: iterativePreconditionerType !<The type of iterative preconditioner \see SolverRoutines_IterativePreconditionerTypes,SolverRoutines
    INTEGER(INTG) :: solutionInitialiseType !<The type of solution vector initialisation \see SolverRoutines_SolutionInitialiseTypes,SolverRoutines
    INTEGER(INTG) :: maximumNumberOfIterations !<The maximum number of iterations
    REAL(DP) :: relativeTolerance !<The relative tolerance between the rhs and residual norm
    REAL(DP) :: absoluteTolerance !<The absolute tolerance of the residual norm
    REAL(DP) :: divergenceTolerance !<The absolute tolerance of the residual norm
    INTEGER(INTG) :: gmresRestart !<The GMRES restart iterations size
    TYPE(PetscPCType) :: pc !<The PETSc preconditioner object
    TYPE(PetscKspType) :: ksp !<The PETSc solver object
  END TYPE LinearIterativeSolverType
  
  !>Contains information for a linear solver
  TYPE LinearSolverType
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG) :: linearSolveType !<The type of linear solver \see SolverRoutines_LinearSolverTypes,SolverRoutines
    LOGICAL :: linkedNewtonPetSCSolver !<Is .TRUE. if this linear solver has been linked from a PETSc Newton solver, .FALSE. if not.
    TYPE(LinearDirectSolverType), POINTER :: directSolver !<A pointer to the direct solver information
    TYPE(LinearIterativeSolverType), POINTER :: iterativeSolver !<A pointer to the iterative solver information
  END TYPE LinearSolverType

  !>Contains information for a Newton line search nonlinear solver
  TYPE NewtonLinesearchSolverType
    TYPE(NewtonSolverType), POINTER :: newtonSolver !<A pointer to the Newton solver
    INTEGER(INTG) :: solverLibrary !<The library type for the Newton linesearch solver \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG) :: linesearchType !<The line search type \see SolverRoutines_NonlinearLineSearchTypes,SolverRoutines
    REAL(DP) :: linesearchAlpha !<The line search alpha
    REAL(DP) :: linesearchMaxstep !<The line search maximum step
    REAL(DP) :: linesearchStepTolerance !<The line search step tolerance
    TYPE(PetscMatColoringType) :: jacobianMatColoring !<The Jacobian matrix colouring
    TYPE(PetscISColoringType) :: jacobianISColoring !<The Jacobian matrix index set colouring
    TYPE(PetscMatFDColoringType) :: jacobianMatFDColoring !<The Jacobian matrix finite difference colouring
    TYPE(PetscSnesType) :: snes !<The PETSc nonlinear solver object
    TYPE(PetscSnesLineSearchType) :: snesLineSearch !<The PETSc SNES line search object
    LOGICAL :: linesearchMonitorOutput !<Flag to determine whether to enable/disable linesearch monitor output.
  END TYPE NewtonLinesearchSolverType
  
  !>Contains information for a Newton trust region nonlinear solver
  TYPE NewtonTrustregionSolverType
    TYPE(NewtonSolverType), POINTER :: newtonSolver !<A pointer to the Newton solver 
    INTEGER(INTG) :: solverLibrary !<The library type for the nonlinear solver \see SolverRoutines_SolverLibraries,SolverRoutines
    REAL(DP) :: trustregionTolerance !<The trust region tolerance
    REAL(DP) :: trustregionDelta0 !<The trust region delta0
    TYPE(PetscSnesType) :: snes !<The PETSc nonlinear solver object
  END TYPE NewtonTrustregionSolverType

  !>Contains information about the convergence test for a newton solver
  TYPE NewtonSolverConvergenceTestType
    REAL(DP) :: energyFirstIter !<The energy for the first iteration
    REAL(DP) :: normalisedEnergy !<The normalized energy for the subsequent iterations
  END TYPE NewtonSolverConvergenceTestType

  !>Contains information for a Newton nonlinear solver
  TYPE NewtonSolverType
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver !<A pointer to the nonlinear solver
    INTEGER(INTG) :: newtonSolveType !<The type of Newton solver
    INTEGER(INTG) :: solutionInitialiseType !<The type of solution vector initialisation \see SolverRoutines_SolutionInitialiseTypes,SolverRoutines
    INTEGER(INTG) :: totalNumberOfFunctionEvaluations !<The number of function evaluations performed by the nonlinear solver
    INTEGER(INTG) :: totalNumberOfJacobianEvaluations !<The number of Jacobian evaluations performed by the nonlinear solver
    INTEGER(INTG) :: maximumNumberOfIterations !<The maximum number of iterations
    INTEGER(INTG) :: maximumNumberOfFunctionEvaluations !<The maximum number of function evaluations
    INTEGER(INTG) :: jacobianCalculationType !<The type of calculation used for the Jacobian \see SolverRoutines_JacobianCalculationTypes,SolverRoutines
    INTEGER(INTG) :: convergenceTestType !<The type of convergence test \see SolverRoutines_NewtonConvergenceTestTypes,SolverRoutines
    REAL(DP) :: absoluteTolerance !<The tolerance between the absolute decrease between the solution norm and the initial guess
    REAL(DP) :: relativeTolerance !<The tolerance between the relative decrease between the solution norm and the initial guess
    REAL(DP) :: solutionTolerance !<The tolerance of the change in the norm of the solution
    TYPE(NewtonSolverConvergenceTestType), POINTER :: convergenceTest !<A pointer to the newton solver convergence test 
    TYPE(NewtonLinesearchSolverType), POINTER :: linesearchSolver !<A pointer to the Newton line search solver information
    TYPE(NewtonTrustregionSolverType), POINTER :: trustregionSolver !<A pointer to the Newton trust region solver information
    TYPE(SolverType), POINTER :: linearSolver !<A pointer to the linked linear solver
    TYPE(SolverType), POINTER :: cellMLEvaluatorSolver !<A pointer to the linked CellML solver
  END TYPE NewtonSolverType

  !>Contains information for a Quasi-Newton line search nonlinear solver
  TYPE QuasiNewtonLinesearchSolverType
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver !<A pointer to the Quasi-Newton solver
    INTEGER(INTG) :: solverLibrary !<The library type for the Quasi-Newton linesearch solver \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG) :: linesearchType !<The line search type \see SolverRoutines_NonlinearLineSearchTypes,SolverRoutines
    REAL(DP) :: linesearchMaxstep !<The line search maximum step
    REAL(DP) :: linesearchStepTolerance !<The line search step tolerance
    TYPE(PetscMatColoringType) :: jacobianMatColoring !<The Jacobian matrix colouring
    TYPE(PetscISColoringType) :: jacobianISColoring !<The Jacobian matrix index set colouring
    TYPE(PetscMatFDColoringType) :: jacobianMatFDColoring !<The Jacobian matrix finite difference colouring
    TYPE(PetscSnesType) :: snes !<The PETSc nonlinear solver object
    TYPE(PetscSnesLineSearchType) :: snesLineSearch !<The PETSc SNES line search object
    LOGICAL :: linesearchMonitorOutput !<Flag to determine whether to enable/disable linesearch monitor output.
  END TYPE QuasiNewtonLinesearchSolverType
  
  !>Contains information for a Quasi-Newton trust region nonlinear solver
  TYPE QuasiNewtonTrustregionSolverType
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver !<A pointer to the Quasi-Newton solver 
    INTEGER(INTG) :: solverLibrary !<The library type for the nonlinear solver \see SolverRoutines_SolverLibraries,SolverRoutines
    REAL(DP) :: trustregionTolerance !<The trust region tolerance
    REAL(DP) :: trustregionDelta0 !<The trust region delta0
    TYPE(PetscSnesType) :: snes !<The PETSc nonlinear solver object
  END TYPE QuasiNewtonTrustregionSolverType
  
  !>Contains information for a Quasi-Newton nonlinear solver
  TYPE QuasiNewtonSolverType
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver !<A pointer to the nonlinear solver
    INTEGER(INTG) :: quasiNewtonSolveType !<The type of Quasi-Newton solver
    INTEGER(INTG) :: quasiNewtonType !<The type of Quasi-Newton variant
    INTEGER(INTG) :: restartType !<The restart type of the Quasi-Newton solver
    INTEGER(INTG) :: restart !<Number of saved states Quasi-Newton solver.
    INTEGER(INTG) :: scaleType !<The scaling type of the Quasi-Newton solver
    INTEGER(INTG) :: solutionInitialiseType !<The type of solution vector initialisation \see SolverRoutines_SolutionInitialiseTypes,SolverRoutines
    INTEGER(INTG) :: totalNumberOfFunctionEvaluations !<The number of function evaluations performed by the nonlinear solver
    INTEGER(INTG) :: totalNumberOfJacobianEvaluations !<The number of Jacobian evaluations performed by the nonlinear solver
    INTEGER(INTG) :: maximumNumberOfIterations !<The maximum number of iterations
    INTEGER(INTG) :: maximumNumberOfFunctionEvaluations !<The maximum number of function evaluations
    INTEGER(INTG) :: jacobianCalculationType !<The type of calculation used for the Jacobian \see SolverRoutines_JacobianCalculationTypes,SolverRoutines
    INTEGER(INTG) :: convergenceTestType !<The type of convergence test \see SolverRoutines_NewtonConvergenceTestTypes,SolverRoutines
    REAL(DP) :: absoluteTolerance !<The tolerance between the absolute decrease between the solution norm and the initial guess
    REAL(DP) :: relativeTolerance !<The tolerance between the relative decrease between the solution norm and the initial guess
    REAL(DP) :: solutionTolerance !<The tolerance of the change in the norm of the solution
    TYPE(NewtonSolverConvergenceTestType), POINTER :: convergenceTest !<A pointer to the (Quasi-)Newton solver convergence test 
    TYPE(QuasiNewtonLinesearchSolverType), POINTER :: linesearchSolver !<A pointer to the Quasi-Newton line search solver information
    TYPE(QuasiNewtonTrustregionSolverType), POINTER :: trustregionSolver !<A pointer to the Quasi-Newton trust region solver information
    TYPE(SolverType), POINTER :: linearSolver !<A pointer to the linked linear solver
    TYPE(SolverType), POINTER :: cellMLEvaluatorSolver !<A pointer to the linked CellML solver
  END TYPE QuasiNewtonSolverType

  !>Contains information for a nonlinear solver
  TYPE NonlinearSolverType
    TYPE(SolverType), POINTER :: solver !<A pointer to the problem_solver
    INTEGER(INTG) :: nonlinearSolveType !<The type of nonlinear solver \see SolverRoutines_NonlinearSolverTypes,SolverRoutines
    TYPE(NewtonSolverType), POINTER :: newtonSolver !<A pointer to the Newton solver information
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver !<A pointer to the Quasi-Newton solver information
  END TYPE NonlinearSolverType
  
  !>Contains information for an eigenproblem solver
  TYPE EigenproblemSolverType
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG) :: solverLibrary !<The library type for the eigenproblem solver \see SolverRoutines_SolverLibraries,SolverRoutines
  END TYPE EigenproblemSolverType
  
  !>Contains information for an optimiser solver
  TYPE OptimiserSolverType
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG) :: solverLibrary !<The library type for the optimiser solver \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG) :: objectiveType !The type of objective for the optimiser solver \see SolverRoutines_OptimiserObjectiveTypes,SolverRoutines
    INTEGER(INTG) :: constraintType !The type of constraints for the optimiser solver \see SolverRoutines_OptimiserConstraintTypes,SolverRoutines
    INTEGER(INTG) :: certaintyType !The type of certainty for the optimiser solver \see SolverRoutines_OptimiserCertaintyTypes,SolverRoutines
    INTEGER(INTG) :: variableType !The type of variable type for the optimiser solver \see SolverRoutines_OptimiserVariableTypes,SolverRoutines
    INTEGER(INTG) :: gradientCalculationType !The type of calculation to obtain the gradient for the optimiser solver \see SolverRoutines_OptimiserGradientCalculationTypes,SolverRoutines
    INTEGER(INTG) :: hessianCalculationType !The type of calculation to obtain the Hessian for the optimiser solver \see SolverRoutines_OptimiserHessianCalculationTypes,SolverRoutines
    INTEGER(INTG) :: totalNumberOfObjectiveEvaluations !<The total number of objective evaluations for the optimiser solver.
    TYPE(PetscTaoType) :: tao !<The PETSc tao optimiser object
  END TYPE OptimiserSolverType

  !>Contains information for a CellML evaluation solver
  TYPE CellMLEvaluatorSolverType
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG) :: solverLibrary !<The library type for the CellML evaluation solver \see SolverRoutines_SolverLibraries,SolverRoutines
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML environment for the solver
  END TYPE CellMLEvaluatorSolverType
  
  !>Contains information for a geometric transformation solver
  TYPE GeometricTransformationSolverType
    TYPE(SolverType), POINTER :: solver !<A pointer to the problem_solver
    LOGICAL :: arbitraryPath !<.TRUE. if the transformation has an arbitrary path, .FALSE. if it's uni-directional(default)
    INTEGER(INTG) :: numberOfIncrements !<The number of increments used to apply the transformation.
    REAL(DP), ALLOCATABLE :: scalings(:) !scaling(loadIncrementIdx), the scaling factors for each load increment, apply the full transformation in 1 load increment if unallocated. Only allocated if there are multiple load steps and if the transformation is uni-directional.
    REAL(DP), ALLOCATABLE :: transformationMatrices(:,:,:) !<transformationMatrices(spatialCoord+1,spatialCoord+1,incrementIdx). 4x4 matrices for 3D transformation, 3x3 for 2D transformation
    TYPE(FieldType), POINTER :: field !<fields to which the geometric transformations are applied 
    INTEGER(INTG) :: fieldVariableType !<The field variable type index to transform
    TYPE(FieldVariableType), POINTER :: fieldVariable !<The field variable to which the geometric transformations are applied.
  END TYPE GeometricTransformationSolverType

   !>A buffer type to allow for an array of pointers to a SolverType \see Types::SolverType
  TYPE SolverPtrType
    TYPE(SolverType), POINTER :: ptr
  END TYPE SolverPtrType

 !>Contains information on the type of solver to be used. \see OpenCMISS::Iron::cmfe_SolverType
  TYPE SolverType
    TYPE(SolversType), POINTER :: solvers !<A pointer to the control loop solvers. Note that if this is a linked solver this will be NULL and solvers should be accessed through the linking solver.
    INTEGER(INTG) :: globalNumber !<The global number of the solver in the list of solvers
    TYPE(SolverType), POINTER :: linkingSolver !<A pointer to any solver that is linking to this solver
    INTEGER(INTG) :: numberOfLinkedSolvers !<The number of linked solvers
    TYPE(SolverPtrType), ALLOCATABLE :: linkedSolvers(:) !<linkedSolvers(linkedSolverIdx). A pointer to the linked solvers
    TYPE(SolverPtrType), ALLOCATABLE :: linkedSolverTypeMap(:) !<linkedSolverTypeMap(linkedSolverIdx). The map from the available linked solver types to the linked solver types that are defined for the solver. linked_solver_idx varies from 1 to SolverRoutines::SOLVER_NUMBER_OF_SOLVER_TYPES. If the particular linked solver type has not been defined on the solver then the linkedSolverTypeMap will be NULL. \see SolverRoutines_SolverTypes
    LOGICAL :: solverFinished !<Is .TRUE. if the solver has finished being created, .FALSE. if not.
    LOGICAL :: solverSetup !<Is .TRUE. if the solver has been setup for the first time, .FALSE. if not.
    TYPE(VARYING_STRING) :: label !<A user defined label for the solver.     
    INTEGER(INTG) :: outputType !<The type of output required \see SolverRoutines_OutputTypes,SolverRoutines    
    INTEGER(INTG) :: solveType !<The type of the solver \see SolverRoutines_SolverTypes,SolverRoutines    
    TYPE(LinearSolverType), POINTER :: linearSolver !<A pointer to the linear solver information
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver !<A pointer to the nonlinear solver information
    TYPE(DynamicSolverType), POINTER :: dynamicSolver !<A pointer to the dynamic solver information
    TYPE(DAESolverType), POINTER :: DAESolver !<A pointer to the differential-algebraic equation solver information    
    TYPE(EigenproblemSolverType), POINTER :: eigenproblemSolver !<A pointer to the eigenproblem solver information
    TYPE(OptimiserSolverType), POINTER :: optimiserSolver !<A pointer to the optimiser solver information
    TYPE(CellMLEvaluatorSolverType), POINTER :: cellMLEvaluatorSolver !<A pointer to the CellML solver information
    TYPE(GeometricTransformationSolverType), POINTER :: geometricTransformationSolver !<A pointer to the geometric transformation solver information
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations !<A pointer to the CellML equations
    TYPE(WorkGroupType), POINTER :: workGroup !<A pointer to the computation work group for the solver
  END TYPE SolverType

  !>Contains information on the solvers to be used in a control loop
  TYPE SolversType
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the parent control loop
    LOGICAL :: solversFinished !<Is .TRUE. if the solvers have finished being created, .FALSE. if not.
    INTEGER(INTG) :: numberOfSolvers !<The number of solvers
    TYPE(SolverPtrType), ALLOCATABLE :: solvers(:) !<A pointer to the solvers information for the problem.
  END TYPE SolversType

  PUBLIC DynamicSolverType

  PUBLIC ForwardEulerDAESolverType,BackwardEulerDAESolverType,ImprovedEulerDAESolverType,EulerDAESolverType

  PUBLIC CrankNicolsonDAESolverType,RungeKuttaDAESolverType,AdamsMoultonDAESolverType,BDFDAESolverType

  PUBLIC RushLarsonDAESolverType,ExternalDAESolverType,DAESolverType

  PUBLIC LinearDirectSolverType,LinearIterativeSolverType,LinearSolverType

  PUBLIC NewtonLinesearchSolverType,NewtonTrustregionSolverType,NewtonSolverType

  PUBLIC NewtonSolverConvergenceTestType

  PUBLIC QuasiNewtonLinesearchSolverType,QuasiNewtonTrustregionSolverType,QuasiNewtonSolverType

  PUBLIC NonlinearSolverType

  PUBLIC EigenproblemSolverType

  PUBLIC OptimiserSolverType

  PUBLIC CellMLEvaluatorSolverType

  PUBLIC GeometricTransformationSolverType

  PUBLIC SolverType,SolverPtrType,SolversType
  
  !
  !================================================================================================================================
  !
  ! Solver mapping types
  !

  !>Information on an equations matrix to solver matrix map
  TYPE EquationsMatrixToSolverMatrixMapType
    INTEGER(INTG) :: equationsMatrixType !<The type of equations matrix \see SolverMappingRoutines_EquationsMatrixTypes,SolverMappingRoutines
    INTEGER(INTG) :: equationsMatrixNumber !<The equations matrix number being mapped.
    INTEGER(INTG) :: solverMatrixNumber !<The solver matrix number being mapped.
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<A pointer to the equations matrix being mapped.
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<A pointer to the solver matrix being mapped.
    TYPE(MatrixRowColCouplingType), POINTER :: equationsColToSolverColsMap(:) !<equationsColToSolverColsMap(columnIdx). The mapping from the columnIdx'th column of this equations matrix to the solver matrix columns. Note for interface matrices this will actually be for the transposed interface matrix i.e., for the interface matrix rows.
  END TYPE EquationsMatrixToSolverMatrixMapType

  !>A buffer type to allow for an array of pointers to a EquationsMatrixToSolverMatrixMapType \see Types::EquationsMatrixToSolverMatrixMapType
  TYPE EquationsMatrixToSolverMatrixPtrType
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: ptr !<A pointer to the equations to solver maps
  END TYPE EquationsMatrixToSolverMatrixPtrType
  
  TYPE InterfaceMatrixToSolverMatrixMapType
    INTEGER(INTG) :: interfaceMatrixNumber !<The interface matrix number being mapped.
    INTEGER(INTG) :: solverMatrixNumber !<The solver matrix number being mapped.
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix !<A pointer to the interface matrix being mapped.
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<A pointer to the solver matrix being mapped.
    TYPE(MatrixRowColCouplingType), POINTER :: interfaceRowToSolverColsMap(:) !<interfaceRowToSolverColsMap(rowIdx). The mapping from the rowIdx'th row of this interface matrix to the solver matrix columns.
  END TYPE InterfaceMatrixToSolverMatrixMapType

  !>A buffer type to allow for an array of pointers to a InterfaceMatrixToSolverMatrixMapType \see Types::InterfaceMatrixToSolverMatrixMapType
  TYPE InterfaceMatrixToSolverMatrixMapPtrType
    TYPE(InterfaceMatrixToSolverMatrixMapType), POINTER :: ptr !<A pointer to the interface to solver maps
  END TYPE InterfaceMatrixToSolverMatrixMapPtrType
  
  TYPE JacobianMatrixToSolverMatrixMapType
    INTEGER(INTG) :: jacobianMatrixNumber !<The Jacobian matrix number being mapped.
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix !<A pointer to the Jacobian matrix being mapped.
    INTEGER(INTG) :: solverMatrixNumber !<The solver matrix number being mapped.
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<A pointer to the solver matrix being mapped.
    TYPE(MatrixRowColCouplingType), POINTER :: jacobianColToSolverColsMap(:) !<jacobianColToSolverColsMap(columnIdx). The mapping from the columnIdx'th column of the Jacobian matrix to the solver matrix columns.
  END TYPE JacobianMatrixToSolverMatrixMapType

  TYPE JacobianMatrixToSolverMatrixMapPtrType
    TYPE(JacobianMatrixToSolverMatrixMapType), POINTER :: ptr !<A pointer to the Jacobian to solver map
  END TYPE JacobianMatrixToSolverMatrixMapPtrType

  !>Contains information on the mappings between field variable dofs inequations and the solver matrix columns (solver dofs) 
  TYPE VariableDOFToSolverDOFsMapType
    INTEGER(INTG), ALLOCATABLE :: dofNumbers(:) !<dofNumbers(variableDOFIdx). The solver DOF number that the variableDOFIdx'th variable DOF is mapped to.
    REAL(DP), ALLOCATABLE :: couplingCoefficients(:) !<couplingCoefficients(variableDOFIdx). The multiplicative constant for the mapping between the variableDOFIdx'th variable DOF and the solver DOF.
    REAL(DP), ALLOCATABLE :: additiveConstants(:) !<additiveConstants(variableDOFIdx). The additive constant for the mapping between the variableDOFIdx'th variable dof and the solver DOF.
  END TYPE VariableDOFToSolverDOFsMapType

  !>A buffer type to the VariableDOFToSolverDOFsMapType \see Types::VariableDOFToSolverDOFsMapType
  TYPE VariableDOFToSolverDOFsMapPtrType
    TYPE(VariableDOFToSolverDOFsMapType), POINTER :: ptr !<A pointer to the varDOFToSolverDOFsMap
  END TYPE VariableDOFToSolverDOFsMapPtrType

  TYPE SMInterfaceConditionToEquationsSetMapType
    INTEGER(INTG) :: interfaceConditionIndex
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    INTEGER(INTG) :: interfaceMatrixNumber
  END TYPE SMInterfaceConditionToEquationsSetMapType
  
  !>Contains information on the equations to solver matrix mappings when indexing by solver matrix number
  TYPE EquationsMatricesToSolverMatrixMapType
    INTEGER(INTG) :: solverMatrixNumber !<The number of the solver matrix for these mappings

    INTEGER(INTG) :: numberOfVariables !<The number of variables mapped to this solver matrix
    INTEGER(INTG), ALLOCATABLE :: variableTypes(:) !<variableTypes(variableIdx). The variable type for the variableIdx'th variable that is mapped to the solver matrix.
    TYPE(FieldVariablePtrType), ALLOCATABLE :: variables(:) !<variables(variableIdx). variables(variableIdx)%ptr points to the variableIdx'th variable that is mapped to the solver matrix.
    TYPE(VariableDOFToSolverDOFsMapPtrType), ALLOCATABLE :: varDOFToSolverDOFsMaps(:) !<varDOFToSolverDOFsMaps(variableIdx)%ptr. The mappings from the variable dofs to the solver dofs for the variableIdx'th variable in the equations set that is mapped to the solver matrix.
    INTEGER(INTG) :: numberOfDynamicMatrices !<The number of dynamic matrices mapped to this solver matrix
    TYPE(EquationsMatrixToSolverMatrixPtrType), ALLOCATABLE :: dynamicMatrixToSolverMatrixMaps(:) !<dynamicMatrixToSolverMatrixMaps(dynamicMatrixIdx). The maps from the dynamicMatrixIdx'th dynamic matrix to solver matrix
    INTEGER(INTG) :: numberOfLinearMatrices !<The number of linear matrices mapped to this solver matrix
    TYPE(EquationsMatrixToSolverMatrixPtrType), ALLOCATABLE :: linearMatrixToSolverMatrixMaps(:) !<linearMatrixToSolverMatrixMaps(linerMatrixIdx). The maps from the linearMatrixIdx'th linear equations matrix to solver matrix
    INTEGER(INTG) :: numberOfJacobianMatrices !<The number of nonlinear equations Jacobian matrices mapped to this solver matrix
    TYPE(JacobianMatrixToSolverMatrixMapPtrType), ALLOCATABLE :: jacobianMatrixToSolverMatrixMaps(:) !<jacobianMatrixToSolverMatrixMaps(jacobianMatrixIdx). The map from the jacobianMatrixIdx'th Jacobian matrix to the solver matrix
  END TYPE EquationsMatricesToSolverMatrixMapType

  !>A buffer type to allow for an array of pointers to EquationsMatricesToSolverMatrixMapType
  TYPE EquationsMatricesToSolverMatrixMapPtrType
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: ptr !<A pointer to EquationsMatricesToSolverMatrixMapType
  END TYPE EquationsMatricesToSolverMatrixMapPtrType
  
  !>Contains information on the equations to solver matrix mappings when indexing by equations matrix number
  TYPE EquationsMatrixToSolverMatricesMapType
    INTEGER(INTG) :: equationsMatrixNumber !<The number of the equations matrix for these mappings
    INTEGER(INTG) :: numberOfSolverMatrices !<The number of solver matrices mapped to this equations matrix
    TYPE(EquationsMatrixToSolverMatrixPtrType), ALLOCATABLE :: equationsMatrixToSolverMatrixMaps(:) !<equationsMatrixToSolverMatrixMaps(solverMatrixIdx). The maps from the equation matrix to the solverMatrixIdx'th solver matrix
  END TYPE EquationsMatrixToSolverMatricesMapType

  !>A buffer type to allow for an array of pointers to EquationsMatrixToSolverMatricesMapType
  TYPE EquationsMatrixToSolverMatricesMapPtrType
    TYPE(EquationsMatrixToSolverMatricesMapType), POINTER :: ptr !<A pointer to EquationsMatrixToSolverMatricesMapType
  END TYPE EquationsMatrixToSolverMatricesMapPtrType
  
  !>Contains information on the mappings from an equations set to the solver matrices
  TYPE EquationsSetToSolverMatricesMapType
    INTEGER(INTG) :: equationsSetIndex !<The index of the equations set for these mappings
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mappings
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations in this equations set
    INTEGER(INTG) :: numberOfInterfaceConditions !<The number of interface conditions affecting this equations set.
    TYPE(SMInterfaceConditionToEquationsSetMapType), ALLOCATABLE :: interfaceConditionToEquationsSetMaps(:) !<interfaceConditionToEquationsSetMaps(interfaceConditionIdx). Information on the interfaceConditionIdx'th interface condition affecting this equations set
    INTEGER(INTG) :: numberOfSolverMatrices !<The number of solver matrices that matrices in this equations set is mapped to
    TYPE(EquationsMatricesToSolverMatrixMapPtrType), ALLOCATABLE :: equationsMatricesToSolverMatrixMaps(:) !<equationsMatricesToSolverMatrixMaps(solverMatrixIdx). The mappings from the equations matrices in this equation set to the solverMatrixIdx'th solver_matrix
    INTEGER(INTG) :: numberOfEquationsMatrices !<The number of equations matrices in this equations set mapped to the solver matrices
    TYPE(EquationsMatrixToSolverMatricesMapPtrType), ALLOCATABLE :: equationsMatrixToSolverMatricesMaps(:) !<equationsMatrixToSolverMatricesMaps(equationsMatrixIdx). The mappings from the equation_matrix_idx'th equations matrix in this equation set to the solver_matrices.
    INTEGER(INTG) :: numberOfJacobianMatrices !<The number of Jacobian matrices in this equatons set mapped to the solver matrices
    TYPE(JacobianMatrixToSolverMatrixMapPtrType), ALLOCATABLE :: jacobianMatrixToSolverMatrixMaps(:) !<jacobianMatrixToSolverMatrixMaps(jacobianMatrixIdx). The mappings from the jacobianMatrixIdx'th Jacobian matrix in this equation set to the solver matrices map.
    TYPE(MatrixRowColCouplingType), POINTER :: equationsRowToSolverRowsMap(:) !<equationsRowToSolverRowsMap(equationsRowIdx). The mappings from the equationsRowIdx'th equations row to the solver matrices rows.
  END TYPE EquationsSetToSolverMatricesMapType

  !>A buffer type to EquationsSetToSolverMatricesMapType \see Types::EquationsSetToSolverMatricesMapType
  TYPE EquationsSetToSolverMatricesMapPtrType
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: ptr !<A pointer to the equationsSetToSolverMatricesMap
  END TYPE EquationsSetToSolverMatricesMapPtrType

  !>Contains information on the interface to solver matrix mappings when indexing by solver matrix number
  TYPE InterfaceMatricesToSolverMatrixMapType
    INTEGER(INTG) :: solverMatrixNumber !<The number of the solver matrix for these mappings
    INTEGER(INTG) :: lagrangeVariableType !<LThe variable type for the Lagrange variable that is mapped to the solver matrix.
    TYPE(FieldVariableType), POINTER :: lagrangeVariable !<A pointer to the Lagrange variable that is mapped to the solver matrix.
    TYPE(VariableDOFToSolverDOFsMapType), POINTER :: lagrangeVarDOFToSolverDOFsMap !<The mappings from the Lagrange variable dofs to the solver dofs.
    INTEGER(INTG) :: numberOfDependentVariables !<The number of dependent variables mapped to this solver matrix
    INTEGER(INTG), ALLOCATABLE :: dependentVariableTypes(:) !<dependentVariableTypes(variableIdx). The variable type for the variableIdx'th dependent variable that is mapped to the solver matrix.
    TYPE(FieldVariablePtrType), ALLOCATABLE :: dependentVariables(:) !<dependentVariables(variableIdx). depndentVariables(variableIdx)%ptr points to the variableIdx'th dependent variable that is mapped to the solver matrix.
    TYPE(VariableDOFToSolverDOFsMapPtrType), ALLOCATABLE :: dependentVarDOFToSolverDOFsMaps(:) !<dependentVarDOFToSolverDOFsMaps(interfaceMatrixIdx). The mappings from the dependent variable dofs to the solver dofs for the interfaceMatrixIdx'th interface matrix dependent variable in the interface condition that is mapped to the solver matrix.
    INTEGER(INTG) :: numberOfInterfaceMatrices !<The number of interface matrices mapped to this solver matrix
    TYPE(InterfaceMatrixToSolverMatrixMapPtrType), ALLOCATABLE :: interfaceMatrixToSolverMatrixMaps(:) !<interfaceMatrixToSolverMatrixMaps(interfaceMatrixIdx). The maps from the interfaceMatrixIdx'th interface matrix to solver matrix
    TYPE(MatrixRowColCouplingType), POINTER :: interfaceColToSolverColsMap(:) !<interfaceColToSolverColsMap(columnIdx). The mapping from the columnIdx'th column of this interface matrix to the solver matrix columns.
  END TYPE InterfaceMatricesToSolverMatrixMapType

  !>A buffer type to allow for arrays of pointers to InterfaceMatricesToSolverMatrixMapType
  TYPE InterfaceMatricesToSolverMatrixMapPtrType
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: ptr !<A pointer to InterfaceMatricesToSolverMatrixMapType
  END TYPE InterfaceMatricesToSolverMatrixMapPtrType
  
  !>Contains information on the interface to solver matrix mappings when indexing by interface matrix number
  TYPE InterfaceMatrixToSolverMatricesMapType
    INTEGER(INTG) :: interfaceMatrixNumber !<The number of the interface  matrix for these mappings
    INTEGER(INTG) :: numberOfSolverMatrices !<The number of solver matrices mapped to this interface matrix
    TYPE(InterfaceMatrixToSolverMatrixMapPtrType), ALLOCATABLE :: interfaceMatrixToSolverMatrixMaps(:) !<interfaceMatrixToSolverMatrixMaps(solverMatrixIdx). The maps from the interface matrix to the solverMatrixIdx'th solver matrix
   TYPE(MatrixRowColCouplingType), POINTER :: interfaceRowToSolverRowsMap(:) !<interfaceRowToSolverRowsMap(interfaceRowIdx). The mapping from the interfaceRowIdx'th interface matrix to a solver row.
  END TYPE InterfaceMatrixToSolverMatricesMapType

  !>A buffer type to allow for arrays of pointers to InterfaceMatrixToSolverMatricesMapType
  TYPE InterfaceMatrixToSolverMatricesMapPtrType
    TYPE(InterfaceMatrixToSolverMatricesMapType), POINTER :: ptr !<A pointer to InterfaceMatrixToSolverMatricesMapType
  END TYPE InterfaceMatrixToSolverMatricesMapPtrType
  
  TYPE SMEquationsSetToInterfaceConditionMapType
    INTEGER(INTG) :: equationsSetIndex
    TYPE(EquationsSetType), POINTER :: equationsSet
    INTEGER(INTG) :: interfaceMatrixIndex
  END TYPE SMEquationsSetToInterfaceConditionMapType
  
  !>Contains information on the mappings from an interface condition to the solver matrices
  TYPE InterfaceConditionToSolverMatricesMapType
    INTEGER(INTG) :: interfaceConditionIndex !<The index of the interface condition for these mappings
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mappings
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations in this interface condition
    INTEGER(INTG) :: numberOfEquationsSets !<The number of equations sets that the interface condition affects
    TYPE(SMEquationsSetToInterfaceConditionMapType), ALLOCATABLE :: equationsSetToInterfaceConditionMaps(:) !<equationsSetToInterfaceConditionMaps(equationsSetIdx). The equations set information of the equationsSetIdx'th equations set that the interface condition affects.
    INTEGER(INTG) :: numberOfSolverMatrices !<The number of solver matrices this interface condition is mapped to
    TYPE(InterfaceMatricesToSolverMatrixMapPtrType), ALLOCATABLE :: interfaceMatricesToSolverMatrixMaps(:) !<interfaceMatricesToSolverMatrixMaps(solverMatrixIdx). The mappings from the interface matrices in this interface condition to the solverMatrixIdx'th solver_matrix
    INTEGER(INTG) :: numberOfInterfaceMatrices !<The number of interface matrices in this interface condition that are mapped to solver matrices
    TYPE(InterfaceMatrixToSolverMatricesMapPtrType), ALLOCATABLE :: interfaceMatrixToSolverMatricesMaps(:) !<interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx). The mappings from the interfaceMatrixIdx'th interface matrix in this interface condition to the solver_matrices.
    TYPE(MatrixRowColCouplingType), POINTER :: interfaceColToSolverRowsMap(:) !<interfaceColToSolverRowsMap(interfaceColumnIdx). The mapping from the interfaceColumnIdx'th interface column to solver rows.
  END TYPE InterfaceConditionToSolverMatricesMapType

  !>A buffer type to InterfaceConditionToSolverMatricesMapType \see Types::InterfaceConditionToSolverMatricesMapType
  TYPE InterfaceConditionToSolverMatricesMapPtrType
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: ptr
  END TYPE InterfaceConditionToSolverMatricesMapPtrType
  
 !>Contains information about the mapping from a solver matrix column to dynamic equations matrices and variables
  TYPE SolverColToDynamicEquationsMapType
    INTEGER(INTG) :: numberOfDynamicMatrices !<The number of dynamic matrices the solver column is mapped to in this equations set
    INTEGER(INTG), ALLOCATABLE :: equationsMatrixNumbers(:) !<equationsMatrixNumbers(i). The i'th equations matrix number in the equations that the solver column is mapped to
    INTEGER(INTG), ALLOCATABLE :: equationsColumnNumbers(:) !<equationsColumnNumbers(i). The i'th equations column number in the equation set the solver column is mapped to.
    REAL(DP), ALLOCATABLE :: couplingCoefficients(:) !<couplingCoefficients(i). The i'th coupling coefficient for solver column mapping
  END TYPE SolverColToDynamicEquationsMapType

  !>A buffer type to allow for arrays of pointers to SolverColToDynamicEquationsMapType
  TYPE SolverColToDynamicEquationsMapPtrType
    TYPE(SolverColToDynamicEquationsMapType), POINTER :: ptr !<The pointer to SolverColToDynamicEquationsMapType
  END TYPE SolverColToDynamicEquationsMapPtrType

  !>Contains information about the mapping from a solver matrix column to linear equations matrices and variables
  TYPE SolverColToLinearEquationsMapType
    INTEGER(INTG) :: numberOfLinearMatrices !<The number of linear matrices the solver column is mapped to in this equations set
    INTEGER(INTG), ALLOCATABLE :: equationsMatrixNumbers(:) !<equationsMatrixNumbers(i). The i'th equations matrix number in the equations that the solver column is mapped to
    INTEGER(INTG), ALLOCATABLE :: equationsColumnNumbers(:) !<equationsColumnNumbers(i). The i'th equations column number in the equation set the solver column is mapped to.
    REAL(DP), ALLOCATABLE :: couplingCoefficients(:) !<couplingCoefficients(i). The i'th coupling coefficient for solver column mapping
  END TYPE SolverColToLinearEquationsMapType

  !>A buffer type to allow for arrays of pointers to SolverColToLinearEquationsMapType
  TYPE SolverColToLinearEquationsMapPtrType
    TYPE(SolverColToLinearEquationsMapType), POINTER :: ptr !<The pointer to SolverColToLinearEquationsMapType
  END TYPE SolverColToLinearEquationsMapPtrType

  !>Contains information about the mapping from a solver matrix column to nonlinear equations matrices and variables
  TYPE SolverColToNonlinearEquationsMapType
    INTEGER(INTG) :: numberOfJacobianMatrices !<The number of Jacobian matrices the solver column is mapped to in this equations set
    INTEGER(INTG), ALLOCATABLE :: jacobianMatrixNumbers(:) !<jacobianMatrixNumbers(i). The i'th Jacobian matrix number in the equations that the solver column is mapped to
    INTEGER(INTG), ALLOCATABLE :: jacobianColumnNumbers(:) !<jacobianColumnNumbers(i). The i'th Jacobian column number in the equation set the solver column is mapped to.
    REAL(DP), ALLOCATABLE :: couplingCoefficients(:) !<couplingCoefficients(i). The i'th coupling coefficient for solver column mapping
  END TYPE SolverColToNonlinearEquationsMapType

  !>A buffer type to allow for arrays of pointers to SolverColToNonlinearEquationsMapType
  TYPE SolverColToNonlinearEquationsMapPtrType
    TYPE(SolverColToNonlinearEquationsMapType), POINTER :: ptr !<The pointer to SolverColToNonlinearEquationsMapType
  END TYPE SolverColToNonlinearEquationsMapPtrType

  !>Contains information about mapping the solver DOF to the field variable DOFs in the equations.
  TYPE SolverDOFToVariableDOFsMapType
    INTEGER(INTG) :: numberOfEquationDOFs !<The number of equation DOFs this solver DOF is mapped to
    INTEGER(INTG), ALLOCATABLE :: equationTypes(:) !<equationTypes(equationDOFIdx). The type of equation of the equationDOFIdx'th dof that the solver dof is mapped to. \see SolverMappingRoutines_EquationsTypes
    INTEGER(INTG), ALLOCATABLE :: equationIndices(:) !<equationIndices(equationDOFIdx). The index of either the equations set or interface condition of the equationDOFIdx'th DOF that this solver DOF is mapped to.
    TYPE(FieldVariablePtrType), ALLOCATABLE :: variable(:) !<variable(equationDOFIdx)%ptr is a pointer to the field variable that the solver dof is mapped to in the equationDOFIdx'th equation
    INTEGER(INTG), ALLOCATABLE :: variableDOF(:) !<variableDOF(equationDOFIdx). The variable dof number that the solver dof is mapped to in the equationDOFIdx'th equation
    REAL(DP), ALLOCATABLE :: variableCoefficient(:) !<variableCoefficient(equationDOFIdx). The mulitplicative coefficient for the mapping to the equationDOFIdx'th equations
    REAL(DP), ALLOCATABLE :: additiveConstant(:) !<additiveConstant(equationDOFIdx). The additive constant for the mapping to the equationDOFIdx'th equations
  END TYPE SolverDOFToVariableDOFsMapType

  !A buffer type to allow for an array of pointers to SolverDOFToVariableDOFsMapType
  TYPE SolverDOFToVariableDOFsMapPtrType
    TYPE(SolverDOFToVariableDOFsMapType), POINTER :: ptr !<A pointer to SolverDOFToVariableDOFsMapType
  END TYPE SolverDOFToVariableDOFsMapPtrType    
  
  !>Contains information about the mappings from a solver matrix to the equations in an equations set
  TYPE SolverMatrixToEquationsSetMapType
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations in the equations set that these columns map to.
    LOGICAL :: haveDynamic !<Is .TRUE. if there are any dynamic equations in the map.
    LOGICAL :: haveLinear !<Is .TRUE. if there are any linear equations in the map.
    LOGICAL :: haveNonlinear !<Is .TRUE. if there are any nonlinear equations in the map.
    TYPE(SolverColToDynamicEquationsMapPtrType), ALLOCATABLE :: solverColToDynamicEquationsMaps(:) !<solverColToDynamicEquationsMaps(columnIdx). The mappings from the columnIdx'th column of the solver matrix to the dynamic equations in the equations set.
    TYPE(SolverColToLinearEquationsMapPtrType), ALLOCATABLE :: solverColToLinearEquationsMaps(:) !<solverColToLinearEquationsMaps(columnIdx). The mappings from the columnIdx'th column of the solver matrix to the linear equations in the equations set.
    TYPE(SolverColToNonlinearEquationsMapPtrType), ALLOCATABLE :: solverColToNonlinearEquationsMaps(:) !<solverColToLinearEquationsMaps(columnIdx). The mappings from the columnIdx'th column of the solver matrix to the nonlinear equations in the equations set.
  END TYPE SolverMatrixToEquationsSetMapType

  !>A buffer type to allow for arrays of pointers to SolverMatrixToEquationsSetMapType
  TYPE SolverMatrixToEquationsSetMapPtrType
    TYPE(SolverMatrixToEquationsSetMapType), POINTER :: ptr !<A pointer to SolverMatrixToEquationsSetMapType
  END TYPE SolverMatrixToEquationsSetMapPtrType
  
  !>Contains information about the mapping from a solver matrix column to interface equations matrices and variables
  TYPE SolverColToInterfaceEquationsMapType
    INTEGER(INTG) :: numberOfInterfaceMatrices !<The number of interface matrices the solver column is mapped to in this interface condition
    INTEGER(INTG), ALLOCATABLE :: interfaceMatrixNumbers(:) !<interfaceMatrixNumbers(i). The i'th interface matrix number in the interface equations that the solver column is mapped to
    INTEGER(INTG), ALLOCATABLE :: interfaceColumnNumbers(:) !<interfaceColumnNumbers(i). The i'th interface interface column number in the interface condition the solver column is mapped to.
    REAL(DP), ALLOCATABLE :: couplingCoefficients(:) !<couplingCoefficients(i). The i'th coupling coefficient for solver column mapping
  END TYPE SolverColToInterfaceEquationsMapType

  !>A buffer type to allow for an array of pointers to SolverColToInterfaceEquationsMapType
  TYPE SolverColToInterfaceEquationsMapPtrType
    TYPE(SolverColToInterfaceEquationsMapType), POINTER :: ptr !<A pointer to SolverColToInterfaceEquationsMapType
  END TYPE SolverColToInterfaceEquationsMapPtrType

  !>Contains information about the mappings from a solver matrix to the equations in an equations set
  TYPE SolverMatrixToInterfaceConditionMapType
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations in the interface conditionthat these columns map to.
    TYPE(SolverColToInterfaceEquationsMapPtrType), ALLOCATABLE :: solverColToInterfaceEquationsMaps(:) !<solverColToInterfaceEquationsMaps(columnIdx). The mappings from the columnIdx'th column of the solver matrix to the interface equations in the interface condition.
  END TYPE SolverMatrixToInterfaceConditionMapType

  !>A buffer type to allow for an array of pointers to SolverMatrixToInterfaceConditionMapType
  TYPE SolverMatrixToInterfaceConditionMapPtrType
    TYPE(SolverMatrixToInterfaceConditionMapType), POINTER :: ptr !<A pointer to SolverMatrixToInterfaceConditionMapType
  END TYPE SolverMatrixToInterfaceConditionMapPtrType
  
  !>Contains information on the mappings from a solver columns to equations 
  TYPE SolverColsToEquationsMapType
    INTEGER(INTG) :: solverMatrixNumber !<The number of this solver matrix
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping for this solver matrix mapping
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<A pointer to the solver matrix being mappind
    INTEGER(INTG) :: numberOfColumns !<The number of columns in this distributed solver matrix
    TYPE(SolverMatrixToEquationsSetMapType), ALLOCATABLE :: solverColToEquationsSetMaps(:) !<solverColToEquationsSetMaps(equationsSetIdx). The solver columns to equations matrix maps for the equationsSetIdx'th equations set.
    TYPE(SolverMatrixToInterfaceConditionMapType), ALLOCATABLE :: solverColToInterfaceConditionMaps(:) !<solverColToInterfaceConditionMaps(interfaceConditionIdx). The solver columns to interface matrix maps for the interfaceConditionIdx'th interface condition.
    INTEGER(INTG) :: numberOfDofs !<The number of local dofs (excluding ghost values) in the solver vector associated with this solver matrix
    INTEGER(INTG) :: totalNumberOfDofs !<The number of local dofs (including ghost values) in the solver vector associated with this solver matrix
    INTEGER(INTG) :: numberOfGlobalDofs !<The number of global dofs in the solver vector associated with this solver matrix.
!!TODO: should this be index by solver dof rather than column???
    TYPE(SolverDOFToVariableDOFsMapType), ALLOCATABLE :: solverDOFToVariableDOFsMap(:) !<solverDOFToVariableDOFsMap(dofIdx). The mappings from the dofIdx'th local solver dof to the field variables in the equations set.
    TYPE(DomainMappingType), POINTER :: columnDOFSMapping !<The domain mapping for solver matrix column dofs
  END TYPE SolverColsToEquationsMapType
  
  !>A buffer type to allow for arrays of pointers to a SolverColToEquationsMapType
  TYPE SolverColsToEquationsMapPtrType
    TYPE(SolverColsToEquationsMapType), POINTER :: ptr !<A pointer to SolverColToEquationsMapType
  END TYPE SolverColsToEquationsMapPtrType

  !>Contains information on the mappings from a solver row to the equations.
  TYPE SolverRowToEquationsMapType
    INTEGER(INTG) :: numberOfEquationsSetRows !<The number of equations set rows the solver row is mapped to. If the rows are interface rows then this will be zero.
    INTEGER(INTG) :: interfaceConditionIndex !<The interface condition index that the solver row is mapped to. If the rows are from an equations set then this will be zero.
    INTEGER(INTG), ALLOCATABLE :: equationsIndex(:) !<equationsIndex(rowColIdx). If the rows are equations set rows then this is the index of the equations set that the rowColIdx'th equations set row belongs to. If the rows are interface rows this will not be allocated.
    INTEGER(INTG), ALLOCATABLE :: rowColNumber(:) !<rowColNumber(rowColIdx). If the rows are equations set rows the rowColIdx'th equations row number that the solver row is mapped to. If the rows are interface rows then the rowColIdx'th column number that the solver row is mapped to.
    REAL(DP), ALLOCATABLE :: couplingCoefficients(:) !<couplingCoefficients(rowColIdx). The rowColIdx'th coupling coefficient for the solver row to equations mapping
  END TYPE SolverRowToEquationsMapType

  !>A buffer type to allow for arrays of pointers to a SolverRowToEquationsMapType
  TYPE SolverRowToEquationsMapPtrType
    TYPE(SolverRowToEquationsMapType), POINTER :: ptr !<A pointer to SolverRowToEquationsMapType
  END TYPE SolverRowToEquationsMapPtrType

  !>Contains information about the cached create values for a solver mapping
  TYPE SolverMappingCreateValuesCacheType
    TYPE(ListPtrType), POINTER :: equationsVariableList(:) !<equationsVariableList(solverMatrixIdx). The list of equations set variables in the solver mapping.
    TYPE(ListType), POINTER :: equationsRHSVariableList !<The list of equations RHS variables in the solver mapping.
    INTEGER, ALLOCATABLE :: dynamicVariableType(:) !<dynamicVariableType(equationsSetIdx). The variable type that is mapped to the dynamic matrices for the equationsSetIdx'th equations set.
    INTEGER(INTG), ALLOCATABLE :: matrixVariableTypes(:,:,:) !<matrixVariableTypes(0:variableType,equationsSetIdx,matrixIdx). The list of matrix variable types in the equationsSetIdx'th equations set for the matrixIdx'th solver matrix. matrixVariableTypes(0,equationsSetIdx,matrixIdx) is the number of variable types in the equationsSetIdx'th equations set mapped to the matrixIdx'th solver matrix and matrixVariableTypes(1..,equationsSetIdx,matrixIdx) is the list of the variable types in the equations set.
    INTEGER(INTG), ALLOCATABLE :: residualVariableTypes(:,:,:) !<residualVariableTypes(0:variableType,residualIdx,equationsSetIdx). The list of residual variable types in the residualIdx'th residual in the equationsSetIdx'th equations set. residualVariableTypes(0,residualIdx,equationsSetIdx) is the number of variable types in the residualIdx'th residual in equationsSetIdx'th equations set and residualVariableTypes(1..,residualIdx,equationsSetIdx) is the list of the variable types in the residual in the equations set.
    INTEGER(INTG), ALLOCATABLE :: rhsVariableType(:) !<rhsVariableType(equationsSetIdx). The variable type that is mapped to the solution RHS for the equationsSetIdx'th equations set
    INTEGER, ALLOCATABLE :: sourceVariableTypes(:,:,:) !<sourceVariableTypes(variableType,sourceIdx,equationsSetIdx). The source variable type that is mapped to the sourceIdx'th source for the equationsSetIdx'th equations set.
    TYPE(ListPtrType), POINTER :: interfaceVariableList(:) !<interfaceVariableList(solverMatrixIdx). The list of interface condition variables in the solver mapping for the solver_matrix idx'th solver matrix.
    TYPE(ListPtrType), POINTER :: interfaceIndices(:) !<interfaceIndices(equationsSetIdx). The list of interface condition indices in the equationsSetIdx'th equations set.
  END TYPE SolverMappingCreateValuesCacheType

  !>Contains information on a field variable in a solver matrix
  TYPE SolverMappingVariableType
    INTEGER(INTG) :: variableIdx !<The index of this variable in the solver mapping variables
    TYPE(FieldVariableType), POINTER :: variable !<The field variable in the solver matrix
    INTEGER(INTG) :: variableType !<The field variable type in the solver matrix
    INTEGER(INTG) :: numberOfEquations !<The number of equations (set or interface) involving this field variable
    INTEGER(INTG), ALLOCATABLE :: equationTypes(:) !<equationTypes(equationsIdx). The type of the equationsIdx'th equation the variable is involved in. \see SolverMappingRoutines_EquationsTypes
    INTEGER(INTG), ALLOCATABLE :: equationIndices(:) !<equationIndices(equationsIdx). The equations set\interface condition index for the equationsIdx'th equation the variable is involved in. 
    INTEGER(INTG), ALLOCATABLE :: numberOfMatrices(:) !<numberOfMatrices(equationsIdx). The number of equations matrices that are mapped to the equationsIdx'th equation the variable is involved in.
    INTEGER(INTG), ALLOCATABLE :: matrixNumbers(:,:) !<matrixIndices(matrixIdx,equationsIdx). The matrix number of the matrixIdx'th matrix in the equationsIdx'th equation the variable is involved in. 
  END TYPE SolverMappingVariableType

  !>A buffer type to a SolverMappingVariableType \see Types::SolverMappingVariableType  
  TYPE SolverMappingVariablePtrType
    TYPE(SolverMappingVariableType), POINTER :: ptr !<A pointer to SolverMappingVariableType
  END TYPE SolverMappingVariablePtrType

  !>Contains information on the variables involved in a solver matrix
  TYPE SolverMappingVariablesType
    INTEGER(INTG) :: numberOfVariables !<The number of variables involved in the solver matrix.
    TYPE(SolverMappingVariablePtrType), ALLOCATABLE :: variables(:) !<variables(variableIdx)%ptr. A pointer to the variable information for the variableIdx'th variable involved in the solver matrix.
  END TYPE SolverMappingVariablesType

  !>Describes the coupled rows or columns in the solver mapping
  TYPE SolverMappingDofCouplingsType
    INTEGER(INTG) :: numberOfCouplings !<The number of couplings in the list.
    INTEGER(INTG) :: capacity !<The allocated length of the couplings array.
    TYPE(BoundaryConditionsCoupledDofsPtrType), ALLOCATABLE :: dofCouplings(:) !<dofCouplings(couplingIdx) is a pointer to the coupled DOF information for the couplingIdx'th coupling in the solver rows or columns.
  END TYPE SolverMappingDofCouplingsType

  TYPE SolverMatrixToEquationsMapType
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer back to the solver mapping
    INTEGER(INTG) :: solverMatrixNumber !<The number of this matrix to equations map
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<A pointer to the solver matrix being mapped
    TYPE(SolverMappingVariablesType), POINTER :: variablesList !<The list of variable for the solverMatrixIdx'th solver matrix.
    INTEGER(INTG) :: numberOfColumns !<The number of columns in this distributed solver matrix
    INTEGER(INTG) :: numberOfDOFs !<The number of local dofs (excluding ghost values) in the solver vector associated with this solver matrix
    INTEGER(INTG) :: totalNumberOfDOFs !<The number of local dofs (including ghost values) in the solver vector associated with this solver matrix
    INTEGER(INTG) :: numberOfGlobalDOFs !<The number of global dofs in the solver vector associated with this solver matrix.
    INTEGER(INTG) :: numberOfEquationsSets !<The number of equations sets mapped to this solver matrix
    TYPE(SolverMatrixToEquationsSetMapPtrType), ALLOCATABLE :: solverMatrixToEquationsSetMaps(:) !<solverMatrixToEquationsSetMaps(equationsSetIdx). The solver matrix to equations set maps for the equationsSetIdx'th equations set.
    INTEGER(INTG) :: numberOfInterfaceConditions !<The number of interface conditions mapped to this sovler matrix
    TYPE(SolverMatrixToInterfaceConditionMapPtrType), ALLOCATABLE :: solverMatrixToInterfaceConditionMaps(:) !<solverMatrixToInterfaceConditionMaps(interfaceConditionIdx). The solver matrix to interface condition maps for the interfaceConditionIdx'th interface condition.
    TYPE(SolverDOFToVariableDOFsMapPtrType), ALLOCATABLE :: solverDOFToVariableDOFsMaps(:) !<solverDOFToVariableDOFsMaps(dofIdx). A pointer to the mappings from the dofIdx'th solver dof to the field variables in the equations.
    TYPE(DomainMappingType), POINTER :: columnDOFSMapping !<The domain mapping for solver matrix column dofs
   END TYPE SolverMatrixToEquationsMapType

  !>A buffer type to a SolverMatrixToEquationsMapType \see Types::SolverMatrixToEquationsMapType
  TYPE SolverMatrixToEquationsMapPtrType
    TYPE(SolverMatrixToEquationsMapType), POINTER :: ptr !<A pointer to the solverMatrixToEquationsMap
  END TYPE SolverMatrixToEquationsMapPtrType

  !>Contains information on the solver mapping between the global equation sets and the solver matrices.
  TYPE SolverMappingType
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations for this mapping.
    LOGICAL :: solverMappingFinished !<Is .TRUE. if the solution mapping has finished being created, .FALSE. if not.
    INTEGER(INTG) :: numberOfSolverMatrices !<The number of solution matrices in this mapping.
    INTEGER(INTG) :: numberOfRows !<The number of (local) rows in the solver matrices
    INTEGER(INTG) :: numberOfGlobalRows !<The number of global rows in the solver matrices
    INTEGER(INTG) :: numberOfEquationsSets !<The number of equations sets in the solution mapping.
    TYPE(EquationsSetPtrType), ALLOCATABLE :: equationsSets(:) !<The list of equations sets that are in this solution mapping
    TYPE(EquationsSetToSolverMatricesMapPtrType), ALLOCATABLE :: equationsSetToSolverMatricesMaps(:) !<equationsSetToSolverMatricesMaps(equationsSetIdx). The mapping from the equationsSetIdx'th equations set to the solver matrices.
    INTEGER(INTG) :: numberOfInterfaceConditions !<The number of interface conditions in the solution mapping.
    TYPE(InterfaceConditionPtrType), ALLOCATABLE :: interfaceConditions(:) !<The list of interface conditions that are in this
    TYPE(InterfaceConditionToSolverMatricesMapPtrType), ALLOCATABLE :: interfaceConditionToSolverMatricesMaps(:) !<interfaceConditionToSolverMatricesMaps(interfaceConditionIdx). The mapping from the interfaceConditionIdx'th interface condition to the solver matrices.
    !!TODO: This should probably be solverMatricesToEquationsMaps
    TYPE(SolverMatrixToEquationsMapPtrType), ALLOCATABLE :: solverMatricesToEquationsMaps(:) !<solverMatricesToEquationsMaps(solverMatrixIx). The mapping from the solverMatrixIdx'th solver matrix to the equations (equations set and interface condition) matrices
    !TYPE(SolverMappingVariablesType), ALLOCATABLE :: variablesList(:) !<variablesList(solverMatrixIdx). The list of variable for the solverMatrixIdx'th solver matrix.
    !TYPE(SolverColsToEquationsMapsPtrType), ALLOCATABLE :: solverColsToEquationsColsMaps(:) !<solverColsToEquationsColsMap(solverMatrixIdx). The mapping from the solverMatrixIdx'th solver matrix columns to the equations columns. 
    TYPE(SolverRowToEquationsMapPtrType), ALLOCATABLE :: solverRowToEquationsMaps(:) !<solverRowToEquationsMaps(localRowIdx). The mappings from the local_row_idx'th solver row to the equations set rows.
    !LOGICAL :: HAVE_JACOBIAN !<Is .TRUE. if the Jacobian exists for nonlinear problems.
    TYPE(DomainMappingType), POINTER :: rowDOFsMapping !<The domain mapping for the solver rows.
    TYPE(SolverMappingVariablesType), POINTER :: rhsVariablesList !<The list of variables for the RHS vector
    TYPE(SolverMappingCreateValuesCacheType), POINTER :: createValuesCache !<The create values cache for the solver mapping
  END TYPE SolverMappingType

  PUBLIC EquationsMatrixToSolverMatrixMapType,EquationsMatrixToSolverMatrixPtrType

  PUBLIC InterfaceMatrixToSolverMatrixMapType,InterfaceMatrixToSolverMatrixMapPtrType

  PUBLIC JacobianMatrixToSolverMatrixMapType,JacobianMatrixToSolverMatrixMapPtrType

  PUBLIC VariableDOFToSolverDOFsMapType,VariableDOFToSolverDOFsMapPtrType

  PUBLIC SMInterfaceConditionToEquationsSetMapType

  PUBLIC EquationsMatricesToSolverMatrixMapType,EquationsMatricesToSolverMatrixMapPtrType

  PUBLIC EquationsMatrixToSolverMatricesMapType,EquationsMatrixToSolverMatricesMapPtrType

  PUBLIC EquationsSetToSolverMatricesMapType,EquationsSetToSolverMatricesMapPtrType

  PUBLIC InterfaceMatricesToSolverMatrixMapType,InterfaceMatricesToSolverMatrixMapPtrType

  PUBLIC InterfaceMatrixToSolverMatricesMapType,InterfaceMatrixToSolverMatricesMapPtrType

  PUBLIC SMEquationsSetToInterfaceConditionMapType

  PUBLIC InterfaceConditionToSolverMatricesMapType,InterfaceConditionToSolverMatricesMapPtrType

  PUBLIC SolverColToDynamicEquationsMapType,SolverColToDynamicEquationsMapPtrType

  PUBLIC SolverColToLinearEquationsMapType,SolverColToLinearEquationsMapPtrType

  PUBLIC SolverColToNonlinearEquationsMapType,SolverColToNonlinearEquationsMapPtrType

  PUBLIC SolverColToInterfaceEquationsMapType,SolverColToInterfaceEquationsMapPtrType

  PUBLIC SolverDOFToVariableDOFsMapType,SolverDOFToVariableDOFsMapPtrType

  PUBLIC SolverMatrixToEquationsSetMapType,SolverMatrixToInterfaceConditionMapType

  PUBLIC SolverMatrixToEquationsSetMapPtrType,SolverMatrixToInterfaceConditionMapPtrType

  !PUBLIC SolverColsToEquationsMapType,SolverColsToEquationsMapPtrType

  PUBLIC SolverRowToEquationsMapType,SolverRowToEquationsMapPtrType

  PUBLIC SolverMappingCreateValuesCacheType

  PUBLIC SolverMappingVariableType,SolverMappingVariablesType

  PUBLIC SolverMappingDofCouplingsType

  PUBLIC SolverMatrixToEquationsMapType,SolverMatrixToEquationsMapPtrType

  PUBLIC SolverMappingType
  
  !
  !================================================================================================================================
  !
  ! History types
  !

  !>Contains information about a history file for a control loop. \see OpenCMISS::Iron::cmfe_HistoryType
  TYPE HISTORY_TYPE
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop for the history file
    LOGICAL :: HISTORY_FINISHED !<Is .TRUE. if the history file has finished being created, .FALSE. if not.
    INTEGER(INTG) :: FILE_FORMAT !<The format of the history file \see HISTORY_ROUTINES_FileFormatTypes,HISTORY_ROUTINES
    TYPE(VARYING_STRING) :: FILENAME !<The file name of the history file
    INTEGER(INTG) :: UNIT_NUMBER !<The unit number of the history file.
    TYPE(FieldType), POINTER :: FIELD !<A pointer to the field that will be read/written in the history file
  END TYPE HISTORY_TYPE

  PUBLIC HISTORY_TYPE
  
  !
  !================================================================================================================================
  !
  ! Control types
  !
  
  !>Contains information on a simple (execute once) control loop
  TYPE ControlLoopSimpleType
    TYPE(ControlLoopType), POINTER :: controlLoop
  END TYPE ControlLoopSimpleType

  !>Contains information on a fixed iteration control loop
  TYPE ControlLoopFixedType
    TYPE(ControlLoopType), POINTER :: controlLoop
    INTEGER(INTG) :: iterationNumber !<The iteration number for this fixed loop progress
    INTEGER(INTG) :: startIteration !<The start iteration for the fixed loop
    INTEGER(INTG) :: stopIteration !<The stop iteration for the fixed loop
    INTEGER(INTG) :: iterationIncrement !<The iteration increment for the fixed loop
    INTEGER(INTG) :: outputNumber !<The frequency of output
    INTEGER(INTG) :: inputNumber !<The frequency of input
  END TYPE ControlLoopFixedType

  !>Contains information on a time iteration control loop
  TYPE ControlLoopTimeType
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer back to the control loop
    INTEGER(INTG) :: iterationNumber !<The iteration number for this time loop progress
    INTEGER(INTG) :: globalIterationNumber !<The total number of iterations for all times the loop has executed
    INTEGER(INTG) :: outputNumber !<The frequency of output
    INTEGER(INTG) :: inputNumber !<The frequency of input
    REAL(DP) :: currentTime !<The current time in the time loop
    REAL(DP) :: startTime !<The start time for the time loop
    REAL(DP) :: stopTime !<The stop time for the time loop
    REAL(DP) :: timeIncrement !<The time increment
    INTEGER(INTG) :: numberOfIterations !<The total number of iterations for this loop, if 0 it will be computed from time span and increment
  END TYPE ControlLoopTimeType

  !>Contains information on a do-while control loop
  TYPE ControlLoopWhileType
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer back to the control loop
    INTEGER(INTG) :: iterationNumber !<The iteration number for this while loop progress
    INTEGER(INTG) :: maximumNumberOfIterations !<The maximum number of iterations for this while loop
    REAL(DP) :: absoluteTolerance !<The absolute tolerance for the while loop continue test
    REAL(DP) :: relativeTolerance !<The relative tolerance for the while loop contintue test
    LOGICAL :: continueLoop !<If .TRUE. the while loop will continue, if .FALSE. the loop will finish
  END TYPE ControlLoopWhileType

  !>Contains information on a load-increment control loop
  TYPE ControlLoopLoadIncrementType
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer back to the control loop
    INTEGER(INTG) :: iterationNumber !<The iteration number for this load increment loop progress
    INTEGER(INTG) :: maximumNumberOfIterations !<The maximum number of iterations for this load increment loop
    INTEGER(INTG) :: outputNumber !<The frequency of output
    INTEGER(INTG) :: inputNumber !<The frequency of input
  END TYPE ControlLoopLoadIncrementType

  !>Contains information about a dependent field variable involved in a control loop solver.
  TYPE ControlLoopFieldVariableType
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable
    INTEGER(INTG) :: timeDependence !<The time dependence type of the field variable in the control loop solvers
    INTEGER(INTG) :: linearity !<The linearity type of the field variable in the control loop solvers.
  END TYPE ControlLoopFieldVariableType

  !>Contains information on the list of dependent field variables involved in the control loop solvers.
  TYPE ControlLoopFieldVariablesType
    INTEGER(INTG) :: numberOfFieldVariables !<The number of field variables in the control loop
    TYPE(ControlLoopFieldVariableType), ALLOCATABLE :: fieldVariables(:) !<fieldVariables(fieldVariableIdx). The information for the fieldVariableIdx'th field variable.
  END TYPE ControlLoopFieldVariablesType

  !>A buffer type to allow for an array of pointers to a ControlLoopType \see Types::ControlLoopType
  TYPE ControlLoopPtrType
    TYPE(ControlLoopType), POINTER :: ptr !<The pointer to the control loop
  END TYPE ControlLoopPtrType

  !>Contains information on a control loop. \see OpenCMISS::Iron::cmfe_ControlLoopType
  TYPE ControlLoopType
    TYPE(ProblemType), POINTER :: problem !<A pointer back to the problem for the control loop
    TYPE(ControlLoopType), POINTER :: parentLoop !<A pointer back to the parent control loop if this is a sub loop
    LOGICAL :: controlLoopFinished !<Is .TRUE. if the problem control has finished being created, .FALSE. if not.
    TYPE(VARYING_STRING) :: label !<A user defined label for the control loop.
    
    INTEGER(INTG) :: loopType !<The type of control loop \see ControlLoopRoutines_ControlTypes,ControlLoopRoutines
    INTEGER(INTG) :: controlLoopLevel !<The level of the control loop
    INTEGER(INTG) :: subLoopIndex !<The position of this loop in the sub loops of the parent loop if this is a sub loop.
    INTEGER(INTG) :: outputType !<The output type of the control loop \see ControlLoopRoutines_OutputTypes,ControlLoopRoutines
    
    TYPE(ControlLoopSimpleType), POINTER :: simpleLoop !<A pointer to the simple loop information
    TYPE(ControlLoopFixedType), POINTER :: fixedLoop !<A pointer to the fixed loop information
    TYPE(ControlLoopTimeType), POINTER :: timeLoop !<A pointer to the time loop information
    TYPE(ControlLoopWhileType), POINTER :: whileLoop !<A pointer to the while loop information
    TYPE(ControlLoopLoadIncrementType), POINTER :: loadIncrementLoop !<A pointer to the load increment loop information
    
    INTEGER(INTG) :: numberOfSubLoops !<The number of control loops below this loop
    TYPE(ControlLoopPtrType), ALLOCATABLE :: subLoops(:) !<A array of pointers to the loops below this loop.

    TYPE(ControlLoopFieldVariablesType), POINTER :: fieldVariables !<A pointer to the field variables information for this control loop.
    TYPE(SolversType), POINTER :: solvers !<A pointer to the solvers for this control loop
    TYPE(HISTORY_TYPE), POINTER :: history !<A pointer to the history file for this control loop.
  END TYPE ControlLoopType

  PUBLIC ControlLoopSimpleType,ControlLoopFixedType,ControlLoopTimeType,ControlLoopWhileType,ControlLoopLoadIncrementType

  PUBLIC ControlLoopFieldVariableType,ControlLoopFieldVariablesType

  PUBLIC ControlLoopType,ControlLoopPtrType
  
  !
  !================================================================================================================================
  !
  ! Problem types
  !

  TYPE ProblemSetupType
    INTEGER(INTG) :: setupType !<The setup type \see ProblemRoutines_SetupTypes,ProblemRoutines
    INTEGER(INTG) :: actionType !<The action type \see ProblemRoutines_SetupActionTypes,ProblemRoutines
  END TYPE ProblemSetupType
  
  !>Contains information for a problem. \see OpenCMISS::Iron::cmfe_ProblemType
  TYPE ProblemType
    INTEGER(INTG) :: userNumber !<The user defined identifier for the problem. The user number must be unique.
    INTEGER(INTG) :: globalNumber !<The global number of the problem in the list of problems.
    LOGICAL :: problemFinished !<Is .TRUE. if the problem has finished being created, .FALSE. if not.
    TYPE(ProblemsType), POINTER :: problems !<A pointer to the problems for this problem.
    INTEGER(INTG), ALLOCATABLE :: specification(:) !<The problem specification array, eg. [class, type, subtype], although there can be more or fewer identifiers. Unused identifiers are set to zero.
    INTEGER(INTG) :: specificationLength !<The number of non-zero specification identitiers
    TYPE(WorkGroupType), POINTER :: workGroup !<The work group to use for the problem.
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop information for the problem.
  END TYPE ProblemType
  
  !>A buffer type to allow for an array of pointers to a ProblemType \see Types:ProblemType
  TYPE ProblemPtrType
    TYPE(ProblemType), POINTER :: ptr !<The pointer to the problem.
  END TYPE ProblemPtrType
       
  !>Contains information on the problems defined.
  TYPE ProblemsType
    TYPE(ContextType), POINTER :: context !<A pointer to the context for the problems
    INTEGER(INTG) :: numberOfProblems !<The number of problems defined.
    TYPE(ProblemPtrType), POINTER :: problems(:) !<The array of pointers to the problems.
  END TYPE ProblemsType

  PUBLIC ProblemSetupType

  PUBLIC ProblemType,ProblemPtrType,ProblemsType

  !
  !================================================================================================================================
  !
  ! Region types
  
  !>A buffer type to allow for an array of pointers to a RegionType.
  TYPE RegionPtrType
    TYPE(RegionType), POINTER :: ptr !<The pointer to the region.
  END TYPE RegionPtrType
     
  !>Contains information for a region. \see OpenCMISS::Iron::cmfe_RegionType
  TYPE RegionType 
    INTEGER(INTG) :: userNumber !<The user defined identifier for the region. The user number must be unique.
    TYPE(RegionsType), POINTER :: regions !<A pointer back to the regions for this region.
    LOGICAL :: regionFinished !<Is .TRUE. if the region has finished being created, .FALSE. if not.
    TYPE(VARYING_STRING) :: label !<A user defined label for the region.
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<A pointer to the coordinate system used by the region.
    TYPE(DataPointSetsType), POINTER :: dataPointSets !<A pointer to the data point sets defined on the region.          
    TYPE(NodesType), POINTER :: nodes !<A pointer to the nodes defined on the region.
    TYPE(MeshesType), POINTER :: meshes !<A pointer to the meshes defined on the region.
    TYPE(GeneratedMeshesType), POINTER :: generatedMeshes !<A pointer to the generated meshes defined on the region.
    TYPE(DecomposersType), POINTER :: decomposers !<A pointer to the decomposers defined on the region.
    TYPE(FieldsType), POINTER :: fields !<A pointer to the fields defined on the region.
    TYPE(EquationsSetsType), POINTER :: equationsSets !<A pointer to the equation sets defined on the region.
    TYPE(CellMLEnvironmentsType), POINTER :: cellMLEnvironments !<A pointer to the CellML environments for the region.
    TYPE(RegionType), POINTER :: parentRegion !<A pointer to the parent region for the region. If the region has no parent region then it is the global (world) region and parentRegion is NULL.
    INTEGER(INTG) :: numberOfSubRegions !<The number of sub-regions defined for the region.
    TYPE(RegionPtrType), POINTER :: subRegions(:) !<An array of pointers to the sub-regions defined on the region. \todo make this allocatable
    TYPE(InterfacesType), POINTER :: interfaces !<A pointer to the interfaces defined on the region.
  END TYPE RegionType

  !>Contains information about the regions
  TYPE RegionsType
    TYPE(ContextType), POINTER :: context !<A pointer to the context for the regions
    TYPE(RegionType), POINTER :: worldRegion !<A pointer to the world region
  END TYPE RegionsType

  PUBLIC RegionType,RegionPtrType,RegionsType
  
  !    
  !================================================================================================================================
  !
  ! Computation types
  
  !>Pointer type to WorkGroupType
  TYPE WorkGroupPtrType
    TYPE(WorkGroupType), POINTER :: ptr
  END TYPE WorkGroupPtrType
  
  !>Contains information on logical working groups
  TYPE WorkGroupType
    INTEGER(INTG) :: userNumber !<The user number of the work group
    LOGICAL :: workGroupFinished !<Is .TRUE. if the work group has been finished. .FALSE. if not. 
    TYPE(WorkGroupType), POINTER:: parentWorkGroup !<Parent of this working groups
    TYPE(VARYING_STRING) :: label !<The label of the work group
    INTEGER(INTG) :: numberOfGroupComputationNodes !<The size of the total computational nodes belonging to this group
    INTEGER(INTG), ALLOCATABLE :: worldRanks(:) !<worldRanks(rankIdx). The rank in the world communicator corresponding to the rankIdx'th group rank. 
    INTEGER(INTG) :: numberOfAvailableRanks !<The number of available ranks for this work group. Note that the numberOfAvaiableRanks plus the sum of the number of ranks in one sub group level down should be the number of nodes in the parent work group.
    INTEGER(INTG), ALLOCATABLE :: availableRanks(:) !<availableRanks(rankIdx). The list of available ranks for this work group.
    INTEGER(INTG) :: numberOfSubGroups !<The number of sub work groups
    TYPE(WorkGroupPtrType), ALLOCATABLE:: subGroups(:) !<subGroups(subgg365GroupIdx). A pointer to the subGroupIdx'th sub work group.
    TYPE(ComputationEnvironmentType), POINTER :: computationEnvironment !<A pointer to the computational environment
    INTEGER(INTG) :: mpiGroupCommunicator !<The MPI communicator for this work group
    INTEGER(INTG) :: mpiGroup !<The MPI communicator for this work group
    INTEGER(INTG) :: myGroupComputationNodeNumber !<The rank number in the group communicator
    INTEGER(INTG) :: myWorldComputationNodeNumber !<The rank number in the world communicator
  END TYPE WorkGroupType

  !>Contains information on a cache heirarchy
  TYPE ComputationCacheType
    INTEGER(INTG) :: numberOfLevels !<The number of levels in the cache hierarchy
    INTEGER(INTG), ALLOCATABLE :: size(:) !<size(levelIdx). The size of the levelIdx'th cache level.
  END TYPE ComputationCacheType

  !>Contains information on a computation node containing a number of processors
  TYPE ComputationNodeType
    INTEGER(INTG) :: numberOfProcessors !<The number of processors for this computation node
    INTEGER(INTG) :: rank !<The MPI rank of this computation node in the world communicator
    TYPE(ComputationCacheType) :: cache !<Information about the caches of this computational node (not currently used).
    INTEGER(INTG) :: nodeNameLength !<The length of the name of the computation node
    CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME) :: nodeName !<The name of the computation node
  END TYPE ComputationNodeType

  !>Contains information on the MPI type to transfer information about a computation node
  TYPE MPIComputationNodeType
    INTEGER(INTG) :: mpiType !<The MPI data type
    INTEGER(INTG) :: numberOfBlocks !<The number of blocks in the MPI data type. This will be equal to 4.
    INTEGER(INTG) :: blockLengths(4) !<The length of each block.
    INTEGER(INTG) :: types(4) !<The data types of each block.
    INTEGER(MPI_ADDRESS_KIND) :: displacements(4) !<The address displacements to each block.
  END TYPE MPIComputationNodeType

  !>Contains information on the computation environment the program is running in.
  TYPE ComputationEnvironmentType
    TYPE(ContextType), POINTER :: context !<A pointer back to the context for the compuation environment
    INTEGER(INTG) :: mpiVersion !<The version of MPI that we are running with
    INTEGER(INTG) :: mpiSubVersion !<The sub-version of MPI that we are running with
    INTEGER(INTG) :: mpiCommWorld !<The clone of the MPI world communicator for OpenCMISS
    INTEGER(INTG) :: mpiGroupWorld !<The group of the cloned MPI world communicator for OpenCMISS
    INTEGER(INTG) :: numberOfWorldComputationNodes !<The number of computation nodes in the world communicator
    INTEGER(INTG) :: myWorldComputationNodeNumber !<The rank of the running process in the world communicator
    TYPE(ComputationNodeType), ALLOCATABLE :: computationNodes(:) !<computationNodes(node_idx). Contains information on the node_idx'th computation node.
    TYPE(MPIComputationNodeType) :: mpiComputationNode !<The MPI data type information to transfer the computation node information.
    TYPE(WorkGroupType), POINTER :: worldWorkGroup !<A pointer to the work group corresponding to the world communicator
  END TYPE ComputationEnvironmentType
  
  PUBLIC WorkGroupType,WorkGroupPtrType

  PUBLIC ComputationCacheType

  PUBLIC ComputationNodeType,MPIComputationNodeType

  PUBLIC ComputationEnvironmentType
  
  !    
  !================================================================================================================================
  !
  ! Context type

  !>Contains information on the OpenCMISS context that contains all base level objects
  TYPE ContextType
    INTEGER(INTG) :: userNumber !<The context user number
    TYPE(ContextsType), POINTER :: contexts !<A pointer to the contexts for the context.
    INTEGER(INTG), ALLOCATABLE :: cmissRandomSeeds(:) !<The current random seeds for OpenCMISS
    TYPE(BasisFunctionsType), POINTER :: basisFunctions !<A pointer to the list of basis functions
    TYPE(ComputationEnvironmentType), POINTER :: computationEnvironment !<A pointer to the computtional environment
    TYPE(ProblemsType), POINTER :: problems !<A pointer to the list of problems
    TYPE(CoordinateSystemsType), POINTER :: coordinateSystems !<A pointer to the list of coordinate systems for the context
    TYPE(RegionsType), POINTER :: regions !<A pointer to the regions for the context
  END TYPE ContextType

  !>A buffer type to allow for an array of pointers to a ContexType \see OpenCMISS::Types::ContexTypes
  TYPE ContextPtrType
    TYPE(ContextType), POINTER :: ptr !<A pointer to the context
  END TYPE ContextPtrType

  !Having a problem with Heap corruption so for now just allocate a fixed number of contexts until corruption found.
  INTEGER(INTG), PARAMETER :: MAXIMUM_NUMBER_OF_CONTEXTS=99
  
  !>Contains information on the OpenCMISS contexts. 
  TYPE ContextsType
    INTEGER(INTG) :: numberOfContexts !<The number of contexts
    !Set aside a fixed number of contexts rather than allocate to avoid Heap corruption for now
    !TYPE(ContextPtrType), ALLOCATABLE :: contexts(:) !<contexts(contexIdx)%ptr is a pointer to the contextIdx'th context.
    TYPE(ContextPtrType) :: contexts(MAXIMUM_NUMBER_OF_CONTEXTS) !<contexts(contexIdx)%ptr is a pointer to the contextIdx'th context.
  END TYPE ContextsType

  PUBLIC MAXIMUM_NUMBER_OF_CONTEXTS
  
  PUBLIC ContextType,ContextPtrType,ContextsType
  
  !  
  !================================================================================================================================
  !
  

END MODULE Types

