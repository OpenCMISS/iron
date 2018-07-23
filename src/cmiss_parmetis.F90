!> \file
!> \author Chris Bradley
!> \brief This module is a CMISS buffer module to the ParMETIS library.
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

!> This module is a CMISS buffer module to the ParMETIS library.
MODULE CmissParMETIS
  
  USE BaseRoutines
  USE Kinds
  USE ISO_VARYING_STRING
  
#include "macros.h"

  IMPLICIT NONE
 
  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  INTERFACE

    FUNCTION ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, ncon, nparts, tpwgts, ubvec, &
      & options, edgecut, part, comm)
#ifdef WIN32
      !DEC$ ATTRIBUTES C, reference, alias:'_ParMETIS_V3_PartKway' :: ParMETIS_V3_PartKway
#endif      
      USE Kinds
      INTEGER(INTG) :: vtxdist(*)
      INTEGER(INTG) :: xadj(*)
      INTEGER(INTG) :: adjncy(*)
      INTEGER(INTG) :: vwgt(*)
      INTEGER(INTG) :: adjwgt(*)
      INTEGER(INTG) :: wgtflag
      INTEGER(INTG) :: numflag
      INTEGER(INTG) :: ncon
      INTEGER(INTG) :: nparts
      REAL(DP) :: tpwgts(*)
      REAL(DP) :: ubvec(*)
      INTEGER(INTG) :: options(*)
      INTEGER(INTG) :: edgecut
      INTEGER(INTG) :: part(*)
      INTEGER(INTG) :: comm

      INTEGER(INTG) :: ParMETIS_V3_PartKway
    END FUNCTION ParMETIS_V3_PartKWay

    FUNCTION ParMETIS_V3_PartMeshKway(elmdist, eptr, eind, elmwgt, wgtflag, numflag, ncon, ncommonnodes, nparts, tpwgts, &
      & ubvec, options, edgecut, part, comm)
#ifdef WIN32
      !DEC$ ATTRIBUTES C, reference, alias:'_ParMETIS_V3_PartMeshKway' :: ParMETIS_V3_PartMeshKway
#endif      
      USE Kinds
      INTEGER(INTG) :: elmdist(*)
      INTEGER(INTG) :: eptr(*)
      INTEGER(INTG) :: eind(*)
      INTEGER(INTG) :: elmwgt(*)
      INTEGER(INTG) :: wgtflag
      INTEGER(INTG) :: numflag
      INTEGER(INTG) :: ncon
      INTEGER(INTG) :: ncommonnodes
      INTEGER(INTG) :: nparts
      REAL(DP) :: tpwgts(*)
      REAL(DP) :: ubvec(*)
      INTEGER(INTG) :: options(*)
      INTEGER(INTG) :: edgecut
      INTEGER(INTG) :: part(*)
      INTEGER(INTG) :: comm

      INTEGER(INTG) :: ParMETIS_V3_PartMeshKway
    END FUNCTION ParMETIS_V3_PartMeshKway
    
  END INTERFACE

  PUBLIC ParMETIS_PartKWay

  PUBLIC ParMETIS_PartMeshKWay
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Buffer routine to the ParMetis ParMETIS_V3_PartKway routine.
  SUBROUTINE ParMETIS_PartKWay(vertexDistance,xadj,adjncy,vertexWeight,adjWeight,weightFlag,numFlag,nCon, &
    & numberParts,tpWeights,ubVec,options,numberEdgesCut,partition,communicator,err,error,*)

    !Argument Variables
    INTEGER(INTG), INTENT(IN) :: vertexDistance(:)
    INTEGER(INTG), INTENT(IN) :: xadj(:)
    INTEGER(INTG), INTENT(IN) :: adjncy(:)
    INTEGER(INTG), INTENT(IN) :: vertexWeight(:)
    INTEGER(INTG), INTENT(IN) :: adjWeight(:)
    INTEGER(INTG), INTENT(IN) :: weightFlag
    INTEGER(INTG), INTENT(IN) :: numFlag
    INTEGER(INTG), INTENT(IN) :: nCon
    INTEGER(INTG), INTENT(IN) :: numberParts
    REAL(DP), INTENT(IN) :: tpWeights(:)
    REAL(DP), INTENT(IN) :: ubVec(:)
    INTEGER(INTG), INTENT(IN) :: options(:)
    INTEGER(INTG), INTENT(OUT) :: numberEdgesCut
    INTEGER(INTG), INTENT(OUT) :: partition(:)
    INTEGER(INTG), INTENT(IN) :: communicator
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    INTEGER(INTG) :: status

    ENTERS("ParMETIS_PartKWay",err,error,*999)

    status=ParMETIS_V3_PartKway(vertexDistance,xadj,adjncy,vertexWeight,adjWeight,weightFlag,numFlag,nCon, &
      & numberParts,tpWeights,ubVec,options,numberEdgesCut,partition,communicator)
    
    IF(status/=1) CALL FlagError("ParMetis error in ParMETIS_V3_PartKway",err,error,*999)
    
    EXITS("ParMETIS_PartKWay")
    RETURN
999 ERRORSEXITS("ParMETIS_PartKWay",err,error)
    RETURN 1
    
  END SUBROUTINE ParMETIS_PartKWay

  !
  !================================================================================================================================
  !

  !>Buffer routine to the ParMetis ParMETIS_V3_PartMeshKway routine.
  SUBROUTINE ParMETIS_PartMeshKWay(elementDistance,elementPtr,elementIndex,elementWeight,weightFlag,numFlag,nCon, &
    & numberCommonNodes,numberParts,tpWeights,ubVec,options,numberEdgesCut,partition,communicator,err,error,*)

    !Argument Variables
    INTEGER(INTG), INTENT(IN) :: elementDistance(:)
    INTEGER(INTG), INTENT(IN) :: elementPtr(:)
    INTEGER(INTG), INTENT(IN) :: elementIndex(:)
    INTEGER(INTG), INTENT(IN) :: elementWeight(:)
    INTEGER(INTG), INTENT(IN) :: weightFlag
    INTEGER(INTG), INTENT(IN) :: numFlag
    INTEGER(INTG), INTENT(IN) :: nCon
    INTEGER(INTG), INTENT(IN) :: numberCommonNodes
    INTEGER(INTG), INTENT(IN) :: numberParts
    REAL(DP), INTENT(IN) :: tpWeights(:)
    REAL(DP), INTENT(IN) :: ubVec(:)
    INTEGER(INTG), INTENT(IN) :: options(:)
    INTEGER(INTG), INTENT(OUT) :: numberEdgesCut
    INTEGER(INTG), INTENT(OUT) :: partition(:)
    INTEGER(INTG), INTENT(IN) :: communicator
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    INTEGER(INTG) :: status

    ENTERS("ParMETIS_PartMeshKWay",err,error,*999)
    
    status=ParMETIS_V3_PartMeshKway(elementDistance,elementPtr,elementIndex,elementWeight,weightFlag,numflag,nCon, &
      & numberCommonNodes,numberParts,tpWeights,ubVec,options,numberEdgesCut,partition,communicator)
    
    IF(status/=1) CALL FlagError("ParMetis error in ParMETIS_V3_PartMeshKway",ERR,ERROR,*999)
    
    EXITS("ParMETIS_PartMeshKWay")
    RETURN
999 ERRORSEXITS("ParMETIS_PartMeshKWay",err,error)
    RETURN 1
    
  END SUBROUTINE ParMETIS_PartMeshKWay

  !
  !================================================================================================================================
  !
    
END MODULE CmissParMETIS
