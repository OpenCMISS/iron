!> \file
!> \author Chris Bradley
!> \brief This module contains all node access method routines.
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

!>This module contains all node access method routines.
MODULE NodeAccessRoutines
  
  USE BaseRoutines
  USE Kinds
  USE ISO_VARYING_STRING
  USE Strings
  USE Trees
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  INTERFACE Node_CheckExists
    MODULE PROCEDURE Nodes_NodeCheckExists
  END INTERFACE Node_CheckExists

  PUBLIC Nodes_NodeCheckExists

  PUBLIC Node_CheckExists

CONTAINS
  
  !
  !================================================================================================================================
  !

  !>Checks that a user node number is defined on the specified region.
  SUBROUTINE Nodes_NodeCheckExists(nodes,userNumber,nodeExists,globalNumber,err,error,*)

    !Argument variables
    TYPE(NodesType), POINTER :: nodes !<A pointer to the nodes to check
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user node number to check if it exists
    LOGICAL, INTENT(OUT) :: nodeExists !<On exit, is .TRUE. if the node user number exists in the region, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: globalNumber !<On exit, if the node exists the global number corresponding to the user node number. If the node does not exist then global number will be 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(TreeNodeType), POINTER :: treeNode
   
    ENTERS("Nodes_NodeCheckExists",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(nodes)) CALL FlagError("Nodes is not associated.",err,error,*999)
#endif    

    nodeExists=.FALSE.
    globalNumber=0
    NULLIFY(treeNode)
    CALL Tree_Search(nodes%nodesTree,userNumber,treeNode,err,error,*999)
    IF(ASSOCIATED(treeNode)) THEN
      CALL Tree_NodeValueGet(nodes%nodesTree,treeNode,globalNumber,err,error,*999)
      nodeExists=.TRUE.
    ENDIF

    EXITS("Nodes_NodeCheckExists")
    RETURN
999 ERRORSEXITS("Nodes_NodeCheckExists",err,error)
    RETURN 1
    
  END SUBROUTINE Nodes_NodeCheckExists
  
  !
  !================================================================================================================================
  !

END MODULE NodeAccessRoutines
