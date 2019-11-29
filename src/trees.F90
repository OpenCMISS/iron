!> \file
!> \author Chris Bradley
!> \brief Implements trees of base types.
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

!> Implements trees of base types.
MODULE Trees

  USE BaseRoutines
  USE Constants
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Kinds
  USE Strings

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup OpenCMISS_TreeConstants OpenCMISS::Iron::Trees::Constants
  !> \brief Trees constants.
  !>@{
  !> \addtogroup Trees_TreeNodeColourTypes Trees::TreeNodeColourTypes
  !> \brief The colour of the tree nodes
  !> \see Trees
  !>@{
  INTEGER(INTG), PARAMETER :: TREE_BLACK_NODE=0 !<The black colour type for a tree node \see Trees_TreeNodeColourTypes,Trees
  INTEGER(INTG), PARAMETER :: TREE_RED_NODE=1 !<The red colour type for a tree node \see Trees_TreeNodeColourTypes,Trees
  !>@}
  !> \addtogroup Trees_TreeNodeInsertStatus Trees::TreeNodeInsertStatus
  !> \brief The insert status for tree nodes
  !> \see Trees
  !>@{
  INTEGER(INTG), PARAMETER :: TREE_NODE_INSERT_SUCESSFUL=1 !<Successful insert status \see Trees_TreeNodeInsertStatus,Trees
  INTEGER(INTG), PARAMETER :: TREE_NODE_DUPLICATE_KEY=2 !<Duplicate key found for those trees that do not allow duplicate keys \see Trees_TreeNodeInsertStatus,Trees
  !>@}
  !> \addtogroup Trees_TreeInsertTypes Trees::TreeInsertTypes
  !> \brief The insert type for a tree
  !> \see Trees
  !>@{
  INTEGER(INTG), PARAMETER :: TREE_DUPLICATES_ALLOWED_TYPE=1 !<Duplicate keys allowed tree type \see Trees_TreeInsertTypes,Trees
  INTEGER(INTG), PARAMETER :: TREE_NO_DUPLICATES_ALLOWED=2 !<No duplicate keys allowed tree type \see Trees_TreeInsertTypes,Trees
  !>@}
  !> \addtogroup Trees_DataType Trees::DataTypes
  !> \brief Data type parameters for a tree.
  !> \see Trees
  !>@{
  INTEGER(INTG), PARAMETER :: TREE_INTG_TYPE=INTEGER_TYPE !<Integer data type for tree values \see Trees_DataType,Trees
  INTEGER(INTG), PARAMETER :: TREE_DP_TYPE=DOUBLE_REAL_TYPE !<Double precision real data type for tree values \see Trees_DataType,Trees
  !>@}
  !>@}
  
  !Module types

  !>Contains information for a node in a binary search tree
  TYPE TreeNodeType
    PRIVATE
    INTEGER(INTG) :: key !<The key for the tree node
    INTEGER(INTG), ALLOCATABLE :: valuesIntg(:) !<The integer values stored at the tree node
    REAL(DP), ALLOCATABLE :: valuesDP(:) !<The double precision real values stored at the tree node
    INTEGER(INTG) :: colour !<The colour of the node for the red-black tree
    TYPE(TreeNodeType), POINTER :: left !<The pointer to the left daughter tree node if any
    TYPE(TreeNodeType), POINTER :: right !<The pointer to the right daughter tree node if any
    TYPE(TreeNodeType), POINTER :: parent !<The pointer to the parent tree node
  END TYPE TreeNodeType

  !>Contains information on a Red-Black binary search tree
  TYPE TreeType
    PRIVATE
    LOGICAL :: treeFinished !<Is .TRUE. if the tree has finished being created, .FALSE. if not.
    INTEGER(INTG) :: insertType !<The insert type for duplicate keys for the tree \see Trees_InsertTypes,Trees
    INTEGER(INTG) :: numberInTree !<The number of items currently in the tree
    INTEGER(INTG) :: maximumNumberOfValues !<The maximum number of values in the tree nodes
    INTEGER(INTG) :: dataType !<The data type of the tree values \see Trees_DataTypes,Trees
    INTEGER(INTG) :: valueInitialiseIntg !<The integer value to initialise the tree node values to.
    REAL(DP) :: valueInitialiseDP !<The double precision real value to initialise the tree node values to.    
    TYPE(TreeNodeType), POINTER :: root !<The pointer to the root of the tree
    TYPE(TreeNodeType), POINTER :: nil !<The pointer to the nil of the tree
  END TYPE TreeType

  !>A buffer type to allow for arrays of pointers to tree types
  TYPE TreePtrType
    TYPE(TreeType), POINTER :: ptr !<The pointer to the tree
  END TYPE TreePtrType

  !Module variables
  
  !Interfaces

  INTERFACE Tree_Detach
    MODULE PROCEDURE Tree_DetachIntg1
    MODULE PROCEDURE Tree_DetachIntg2
    MODULE PROCEDURE Tree_DetachDP1
    MODULE PROCEDURE Tree_DetachDP2
  END INTERFACE Tree_Detach
  
  INTERFACE Tree_DetachAndDestroy
    MODULE PROCEDURE Tree_DetachAndDestroyIntg1
    MODULE PROCEDURE Tree_DetachAndDestroyIntg2
    MODULE PROCEDURE Tree_DetachAndDestroyDP1
    MODULE PROCEDURE Tree_DetachAndDestroyDP2
  END INTERFACE Tree_DetachAndDestroy
  
  INTERFACE Tree_DetachInOrder
    MODULE PROCEDURE Tree_DetachInOrderIntg1
    MODULE PROCEDURE Tree_DetachInOrderIntg2
    MODULE PROCEDURE Tree_DetachInOrderDP1
    MODULE PROCEDURE Tree_DetachInOrderDP2
  END INTERFACE Tree_DetachInOrder
  
  INTERFACE Tree_ItemInsert
    MODULE PROCEDURE Tree_ItemInsertIntg0
    MODULE PROCEDURE Tree_ItemInsertIntg1
    MODULE PROCEDURE Tree_ItemInsertDP0
    MODULE PROCEDURE Tree_ItemInsertDP1
  END INTERFACE Tree_ItemInsert
  
  INTERFACE Tree_NodeValueGet
    MODULE PROCEDURE Tree_NodeValueGetIntg0
    MODULE PROCEDURE Tree_NodeValueGetIntg1
    MODULE PROCEDURE Tree_NodeValueGetDP0
    MODULE PROCEDURE Tree_NodeValueGetDP1
  END INTERFACE Tree_NodeValueGet

  INTERFACE Tree_NodeValueSet
    MODULE PROCEDURE Tree_NodeValueSetIntg0
    MODULE PROCEDURE Tree_NodeValueSetIntg1
    MODULE PROCEDURE Tree_NodeValueSetDP0
    MODULE PROCEDURE Tree_NodeValueSetDP1
  END INTERFACE Tree_NodeValueSet
  
  INTERFACE Tree_ValueInitialiseGet
    MODULE PROCEDURE Tree_ValueInitialiseGetIntg
    MODULE PROCEDURE Tree_ValueInitialiseGetDP
  END INTERFACE Tree_ValueInitialiseGet
  
  INTERFACE Tree_ValueInitialiseSet
    MODULE PROCEDURE Tree_ValueInitialiseSetIntg
    MODULE PROCEDURE Tree_ValueInitialiseSetDP
  END INTERFACE Tree_ValueInitialiseSet
  
  PUBLIC TreeType,TreeNodeType,TreePtrType

  PUBLIC TREE_NODE_INSERT_SUCESSFUL,TREE_NODE_DUPLICATE_KEY

  PUBLIC TREE_DUPLICATES_ALLOWED_TYPE,TREE_NO_DUPLICATES_ALLOWED

  PUBLIC TREE_INTG_TYPE,TREE_DP_TYPE
  
  PUBLIC Tree_CreateFinish,Tree_CreateStart

  PUBLIC Tree_DataTypeGet,Tree_DataTypeSet

  PUBLIC Tree_Destroy

  PUBLIC Tree_Detach,Tree_DetachAndDestroy

  PUBLIC Tree_InsertTypeSet

  PUBLIC Tree_ItemDelete,Tree_ItemInsert

  PUBLIC Tree_NodeKeyGet

  PUBLIC Tree_NodeValueGet,Tree_NodeValueSet

  PUBLIC Tree_Output

  PUBLIC Tree_Search

  PUBLIC Tree_ValueInitialiseGet,Tree_ValueInitialiseSet

CONTAINS
  
  !
  !=================================================================================================================================
  !

  !>Assert that a tree has been finished
  SUBROUTINE Tree_AssertIsFinished(tree,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER, INTENT(IN) :: tree !<The tree to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Tree_AssertIsFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(tree)) CALL FlagError("Tree is not associated.",err,error,*999)

    IF(.NOT.tree%treeFinished) CALL FlagError("Tree has not been finished.",err,error,*999)
    
    EXITS("Tree_AssertIsFinished")
    RETURN
999 ERRORSEXITS("Tree_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a tree has not been finished
  SUBROUTINE Tree_AssertNotFinished(tree,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER, INTENT(IN) :: tree !<The tree to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Tree_AssertNotFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(tree)) CALL FlagError("Tree is not associated.",err,error,*999)

    IF(tree%treeFinished) CALL FlagError("Tree has already been finished.",err,error,*999)
    
    EXITS("Tree_AssertNotFinished")
    RETURN
999 ERRORSEXITS("Tree_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Assert that a tree has an integer data type
  SUBROUTINE Tree_AssertIsIntgData(tree,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER, INTENT(IN) :: tree !<The tree to assert the integer data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Tree_AssertIsIntgData",err,error,*999)

    IF(.NOT.ASSOCIATED(tree)) CALL FlagError("Tree is not associated.",err,error,*999)

    IF(tree%dataType/=TREE_INTG_TYPE) THEN
      localError="The tree data type of "//TRIM(NumberToVString(tree%dataType,"*",err,error))// &
        & " does not correspond to the required integer data type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Tree_AssertIsIntgData")
    RETURN
999 ERRORSEXITS("Tree_AssertIsIntgData",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_AssertIsIntgData

  !
  !================================================================================================================================
  !

  !>Assert that a tree has a double precision real data type
  SUBROUTINE Tree_AssertIsDPData(tree,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER, INTENT(IN) :: tree !<The tree to assert the double precision real data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Tree_AssertIsDPData",err,error,*999)

    IF(.NOT.ASSOCIATED(tree)) CALL FlagError("Tree is not associated.",err,error,*999)

    IF(tree%dataType/=TREE_DP_TYPE) THEN
      localError="The tree data type of "//TRIM(NumberToVString(tree%dataType,"*",err,error))// &
        & " does not correspond to the required double precision real data type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Tree_AssertIsDPData")
    RETURN
999 ERRORSEXITS("Tree_AssertIsDPData",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_AssertIsDPData

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a tree created with Tree_CreateStart \see{Trees::Tree_CreateStart}.
  SUBROUTINE Tree_CreateFinish(tree,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree to finish
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    
    ENTERS("Tree_CreateFinish",err,error,*998)

    CALL Tree_AssertNotFinished(tree,err,error,*998)

    !Allocate the nil tree node
    ALLOCATE(tree%nil,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate NIL tree node.",err,error,*999)
    CALL Tree_NodeInitialise(tree,tree%nil,err,error,*999)
    tree%nil%key=-99999999 !Set it to something identifiable for debugging
    tree%nil%left=>tree%nil
    tree%nil%right=>tree%nil
    tree%nil%parent=>tree%nil
    !Set the root tree node to NIL        
    tree%root=>tree%nil
    !Finish the tree creation
    tree%treeFinished=.TRUE.

    EXITS("Tree_CreateFinish")
    RETURN
999 CALL Tree_Finalise(tree,err,error,*998)
998 ERRORSEXITS("Tree_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_CreateFinish

  !
  !================================================================================================================================
  !

  !>Starts the creation of a tree and returns a pointer to the created tree \see{Trees::Tree_CreateFinish}.
  SUBROUTINE Tree_CreateStart(tree,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree to create. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variable

    ENTERS("Tree_CreateStart",err,error,*998)

    IF(ASSOCIATED(tree)) CALL FlagError("Tree is already associated.",err,error,*998)
    
    ALLOCATE(tree,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate tree.",err,error,*999)
    CALL Tree_Initialise(tree,err,error,*999)
    !Set Defaults
    tree%insertType=TREE_DUPLICATES_ALLOWED_TYPE

    EXITS("Tree_CreateStart")
    RETURN
999 CALL Tree_Finalise(tree,err,error,*998)
998 ERRORSEXITS("Tree_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_CreateStart

  !
  !================================================================================================================================
  !

  !>Gets the data type for a tree.
  SUBROUTINE Tree_DataTypeGet(tree,dataType,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER, INTENT(IN) :: tree !<A pointer to the tree 
    INTEGER(INTG), INTENT(OUT) :: dataType !<On return, the data type of the tree \see Trees_DataTypes,Trees
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Tree_DataTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(tree)) CALL FlagError("Tree is not associated.",err,error,*999)

    dataType=tree%dataType

    EXITS("Tree_DataTypeGet")
    RETURN
999 ERRORSEXITS("Tree_DataTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_DataTypeGet

  !
  !================================================================================================================================
  !

  !>Sets/changes the data type for a tree.
  SUBROUTINE Tree_DataTypeSet(tree,dataType,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER, INTENT(IN) :: tree !<A pointer to the tree 
    INTEGER(INTG), INTENT(IN) :: dataType !<The data type of the tree to set \see Trees_DataTypes,Trees
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Tree_DataTypeSet",err,error,*999)

    CALL Tree_AssertNotFinished(tree,err,error,*999)
    
    SELECT CASE(dataType)
    CASE(TREE_INTG_TYPE)
      tree%dataType=TREE_INTG_TYPE
    CASE(TREE_DP_TYPE)
      tree%dataType=TREE_DP_TYPE
    CASE DEFAULT
      localError="The data type of "//TRIM(NumberToVString(dataType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Tree_DataTypeSet")
    RETURN
999 ERRORSEXITS("Tree_DataTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_DataTypeSet

  !
  !================================================================================================================================
  !

  !>Destroys a tree
  SUBROUTINE Tree_Destroy(tree,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    
    ENTERS("Tree_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(tree)) CALL FlagError("Tree is not associated.",err,error,*999)
    
    CALL Tree_Finalise(tree,err,error,*999)

    EXITS("Tree_Destroy")
    RETURN
999 ERRORSEXITS("Tree_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_Destroy

  !
  !================================================================================================================================
  !

  !>Detaches the integer tree values and returns them as a pointer to an array
  SUBROUTINE Tree_DetachIntg1(tree,numberInTree,treeValues,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree to detach
    INTEGER(INTG), INTENT(OUT) :: numberInTree !<On exit, the number in the array that has been detached
    INTEGER(INTG), POINTER :: treeValues(:) !<On exit, a pointer to the dettached tree values. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Tree_DetachIntg1",err,error,*998)

    CALL Tree_AssertIsFinished(tree,err,error,*998)
    CALL Tree_AssertIsIntgData(tree,err,error,*998)
    IF(ASSOCIATED(treeValues)) CALL FlagError("Tree values is already associated.",err,error,*998)
    IF(tree%maximumNumberOfValues/=1) THEN
      localError="The maximum number of values in the tree nodes is "// &
        & TRIM(NumberToVString(tree%maximumNumberOfValues,"*",err,error))// &
        & " which does not match the rank of one of the supplied values pointer."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    NULLIFY(treeValues)
    ALLOCATE(treeValues(tree%numberInTree),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate tree values.",err,error,*999)
    treeValues=tree%valueInitialiseIntg
    numberInTree=0
    CALL Tree_DetachInOrderIntg1(tree,tree%root,numberInTree,treeValues,err,error,*999)

    EXITS("Tree_DetachIntg1")
    RETURN
999 IF(ASSOCIATED(treeValues)) DEALLOCATE(treeValues)
    numberInTree=0
998 ERRORSEXITS("Tree_DetachIntg1",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_DetachIntg1

  !
  !================================================================================================================================
  !

  !>Detaches the integer tree values and returns them as a pointer to the an array
  SUBROUTINE Tree_DetachIntg2(tree,numberInTree,treeValues,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree to detach
    INTEGER(INTG), INTENT(OUT) :: numberInTree !<On exit, the number in the array that has been detached
    INTEGER(INTG), POINTER :: treeValues(:,:) !<On exit, a pointer to the dettached tree values. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    
    ENTERS("Tree_DetachIntg2",err,error,*998)

    CALL Tree_AssertIsFinished(tree,err,error,*998)
    CALL Tree_AssertIsIntgData(tree,err,error,*998)
    IF(ASSOCIATED(treeValues)) CALL FlagError("Tree values is already associated.",err,error,*998)
    
    NULLIFY(treeValues)
    ALLOCATE(treeValues(tree%numberInTree,tree%maximumNumberOfValues),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate tree values.",err,error,*999)
    treeValues=tree%valueInitialiseIntg
    numberInTree=0
    CALL Tree_DetachInOrderIntg2(tree,tree%root,numberInTree,treeValues,err,error,*999)

    EXITS("Tree_DetachIntg2")
    RETURN
999 IF(ASSOCIATED(treeValues)) DEALLOCATE(treeValues)
    numberInTree=0
998 ERRORSEXITS("Tree_DetachIntg2",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_DetachIntg2

  !
  !================================================================================================================================
  !

  !>Detaches the double precision real tree values and returns them as a pointer to an array
  SUBROUTINE Tree_DetachDP1(tree,numberInTree,treeValues,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree to detach
    INTEGER(INTG), INTENT(OUT) :: numberInTree !<On exit, the number in the array that has been detached
    REAL(DP), POINTER :: treeValues(:) !<On exit, a pointer to the dettached tree values. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Tree_DetachDP1",err,error,*998)

    CALL Tree_AssertIsFinished(tree,err,error,*998)
    CALL Tree_AssertIsDPData(tree,err,error,*998)
    IF(ASSOCIATED(treeValues)) CALL FlagError("Tree values is already associated.",err,error,*998)
    IF(tree%maximumNumberOfValues/=1) THEN
      localError="The maximum number of values in the tree nodes is "// &
        & TRIM(NumberToVString(tree%maximumNumberOfValues,"*",err,error))// &
        & " which does not match the rank of one of the supplied values pointer."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    NULLIFY(treeValues)
    ALLOCATE(treeValues(tree%numberInTree),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate tree values.",err,error,*999)
    treeValues=tree%valueInitialiseDP
    numberInTree=0
    CALL Tree_DetachInOrderDP1(tree,tree%root,numberInTree,treeValues,err,error,*999)

    EXITS("Tree_DetachDP1")
    RETURN
999 IF(ASSOCIATED(treeValues)) DEALLOCATE(treeValues)
    numberInTree=0
998 ERRORSEXITS("Tree_DetachDP1",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_DetachDP1

  !
  !================================================================================================================================
  !

  !>Detaches the double precision real tree values and returns them as a pointer to the an array
  SUBROUTINE Tree_DetachDP2(tree,numberInTree,treeValues,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree to detach
    INTEGER(INTG), INTENT(OUT) :: numberInTree !<On exit, the number in the array that has been detached
    REAL(DP), POINTER :: treeValues(:,:) !<On exit, a pointer to the dettached tree values. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    
    ENTERS("Tree_DetachDP2",err,error,*998)

    CALL Tree_AssertIsFinished(tree,err,error,*998)
    CALL Tree_AssertIsDPData(tree,err,error,*998)
    IF(ASSOCIATED(treeValues)) CALL FlagError("Tree values is already associated.",err,error,*998)
    
    NULLIFY(treeValues)
    ALLOCATE(treeValues(tree%numberInTree,tree%maximumNumberOfValues),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate tree values.",err,error,*999)
    treeValues=tree%valueInitialiseDP
    numberInTree=0
    CALL Tree_DetachInOrderDP2(tree,tree%root,numberInTree,treeValues,err,error,*999)

    EXITS("Tree_DetachDP2")
    RETURN
999 IF(ASSOCIATED(treeValues)) DEALLOCATE(treeValues)
    numberInTree=0
998 ERRORSEXITS("Tree_DetachDP2",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_DetachDP2

  !
  !================================================================================================================================
  !

  !>Detaches the integer tree values and returns them as a pointer to the an array and then destroys the tree
  SUBROUTINE Tree_DetachAndDestroyIntg1(tree,numberInTree,treeValues,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree to detach and destroy
    INTEGER(INTG), INTENT(OUT) :: numberInTree !<On exit, the number in the array that has been detached
    INTEGER(INTG), POINTER :: treeValues(:) !<On exit, a pointer to the dettached tree values. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Tree_DetachAndDestroyIntg1",err,error,*998)

    CALL Tree_AssertIsFinished(tree,err,error,*998)
    CALL Tree_AssertIsIntgData(tree,err,error,*998)
    IF(ASSOCIATED(treeValues)) CALL FlagError("Tree values is already associated.",err,error,*998)
    IF(tree%maximumNumberOfValues/=1) THEN
      localError="The maximum number of values in the tree nodes is "// &
        & TRIM(NumberToVString(tree%maximumNumberOfValues,"*",err,error))// &
        & " which does not match the rank of one of the supplied values pointer."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    NULLIFY(treeValues)
    ALLOCATE(treeValues(tree%numberInTree),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate tree values.",err,error,*999)
    treeValues=tree%valueInitialiseIntg
    numberInTree=0
    CALL Tree_DetachInOrderIntg1(tree,tree%root,numberInTree,treeValues,err,error,*999)
    CALL Tree_Finalise(tree,err,error,*999)

    EXITS("Tree_DetachAndDestroyIntg1")
    RETURN
999 IF(ASSOCIATED(treeValues)) DEALLOCATE(treeValues)
    numberInTree=0
998 ERRORSEXITS("Tree_DetachAndDestroyIntg1",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_DetachAndDestroyIntg1

  !
  !================================================================================================================================
  !

  !>Detaches the integer tree values and returns them as a pointer to the an array and then destroys the tree
  SUBROUTINE Tree_DetachAndDestroyIntg2(tree,numberInTree,treeValues,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree to detach and destroy
    INTEGER(INTG), INTENT(OUT) :: numberInTree !<On exit, the number in the array that has been detached
    INTEGER(INTG), POINTER :: treeValues(:,:) !<On exit, a pointer to the dettached tree values. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    
    ENTERS("Tree_DetachAndDestroyIntg2",err,error,*998)

    CALL Tree_AssertIsFinished(tree,err,error,*998)
    CALL Tree_AssertIsIntgData(tree,err,error,*998)
    IF(ASSOCIATED(treeValues)) CALL FlagError("Tree values is already associated.",err,error,*998)
    
    NULLIFY(treeValues)
    ALLOCATE(treeValues(tree%numberInTree,tree%maximumNumberOfValues),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate tree values.",err,error,*999)
    treeValues=tree%valueInitialiseIntg
    numberInTree=0
    CALL Tree_DetachInOrderIntg2(tree,tree%root,numberInTree,treeValues,err,error,*999)
    CALL Tree_Finalise(tree,err,error,*999)

    EXITS("Tree_DetachAndDestroyIntg2")
    RETURN
999 IF(ASSOCIATED(treeValues)) DEALLOCATE(treeValues)
    numberInTree=0
998 ERRORSEXITS("Tree_DetachAndDestroyIntg2",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_DetachAndDestroyIntg2

  !
  !================================================================================================================================
  !

  !>Detaches the double precision real tree values and returns them as a pointer to the an array and then destroys the tree
  SUBROUTINE Tree_DetachAndDestroyDP1(tree,numberInTree,treeValues,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree to detach and destroy
    INTEGER(INTG), INTENT(OUT) :: numberInTree !<On exit, the number in the array that has been detached
    REAL(DP), POINTER :: treeValues(:) !<On exit, a pointer to the dettached tree values. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Tree_DetachAndDestroyDP1",err,error,*998)

    CALL Tree_AssertIsFinished(tree,err,error,*998)
    CALL Tree_AssertIsDPData(tree,err,error,*998)
    IF(ASSOCIATED(treeValues)) CALL FlagError("Tree values is already associated.",err,error,*998)
    IF(tree%maximumNumberOfValues/=1) THEN
      localError="The maximum number of values in the tree nodes is "// &
        & TRIM(NumberToVString(tree%maximumNumberOfValues,"*",err,error))// &
        & " which does not match the rank of one of the supplied values pointer."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    NULLIFY(treeValues)
    ALLOCATE(treeValues(tree%numberInTree),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate tree values.",err,error,*999)
    treeValues=tree%valueInitialiseDP
    numberInTree=0
    CALL Tree_DetachInOrderDP1(tree,tree%root,numberInTree,treeValues,err,error,*999)
    CALL Tree_Finalise(tree,err,error,*999)

    EXITS("Tree_DetachAndDestroyDP1")
    RETURN
999 IF(ASSOCIATED(treeValues)) DEALLOCATE(treeValues)
    numberInTree=0
998 ERRORSEXITS("Tree_DetachAndDestroyDP1",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_DetachAndDestroyDP1

  !
  !================================================================================================================================
  !

  !>Detaches the double precision real tree values and returns them as a pointer to the an array and then destroys the tree
  SUBROUTINE Tree_DetachAndDestroyDP2(tree,numberInTree,treeValues,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree to detach and destroy
    INTEGER(INTG), INTENT(OUT) :: numberInTree !<On exit, the number in the array that has been detached
    REAL(DP), POINTER :: treeValues(:,:) !<On exit, a pointer to the dettached tree values. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    
    ENTERS("Tree_DetachAndDestroyDP2",err,error,*998)

    CALL Tree_AssertIsFinished(tree,err,error,*998)
    CALL Tree_AssertIsDPData(tree,err,error,*998)
    IF(ASSOCIATED(treeValues)) CALL FlagError("Tree values is already associated.",err,error,*998)
    
    NULLIFY(treeValues)
    ALLOCATE(treeValues(tree%numberInTree,tree%maximumNumberOfValues),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate tree values.",err,error,*999)
    treeValues=tree%valueInitialiseDP
    numberInTree=0
    CALL Tree_DetachInOrderDP2(tree,tree%root,numberInTree,treeValues,err,error,*999)
    CALL Tree_Finalise(tree,err,error,*999)

    EXITS("Tree_DetachAndDestroyDP2")
    RETURN
999 IF(ASSOCIATED(treeValues)) DEALLOCATE(treeValues)
    numberInTree=0
998 ERRORSEXITS("Tree_DetachAndDestroyDP2",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_DetachAndDestroyDP2

  !
  !================================================================================================================================
  !

  !>Detaches the integer tree values in order from the specified tree node and adds them to the tree values array
  RECURSIVE SUBROUTINE Tree_DetachInOrderIntg1(tree,x,count,treeValues,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree to detach
    TYPE(TreeNodeType), POINTER :: x !<A pointer to the specified tree node to detach from
    INTEGER(INTG), INTENT(INOUT) :: count !<The current number in the detached tree values array
    INTEGER(INTG), INTENT(INOUT) :: treeValues(:) !<The current detached tree values array
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Tree_DetachInOrderIntg1",err,error,*999)

    CALL Tree_AssertIsIntgData(tree,err,error,*999)
    
    IF(.NOT.ASSOCIATED(x,tree%nil)) THEN
      CALL Tree_DetachInOrderIntg1(tree,x%left,count,treeValues,err,error,*999)
      count=count+1
      IF(count<=SIZE(treeValues,1)) THEN
        IF(.NOT.ALLOCATED(x%valuesIntg)) THEN
          localError="Tree node values is not allocated for key "//TRIM(NumberToVString(x%key,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        treeValues(count)=x%valuesIntg(1)
      ELSE
        localError="The current count of the tree values of "//TRIM(NumberToVString(count,"*",err,error))// &
          & " is greater than the size of the tree values array of "// &
          & TRIM(NumberToVString(SIZE(treeValues,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      CALL Tree_DetachInOrderIntg1(tree,x%right,count,treeValues,err,error,*999)
    ENDIF

    EXITS("Tree_DetachInOrderIntg1")
    RETURN
999 ERRORSEXITS("Tree_DetachInOrderIntg1",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_DetachInOrderIntg1

  !
  !================================================================================================================================
  !

  !>Detaches the integer tree values in order from the specified tree node and adds them to the tree values array
  RECURSIVE SUBROUTINE Tree_DetachInOrderIntg2(tree,x,count,treeValues,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree to detach
    TYPE(TreeNodeType), POINTER :: x !<A pointer to the specified tree node to detach from
    INTEGER(INTG), INTENT(INOUT) :: count !<The current number in the detached tree values array
    INTEGER(INTG), INTENT(INOUT) :: treeValues(:,:) !<The current detached tree values array
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    INTEGER(INTG) :: valueIdx
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Tree_DetachInOrderIntg2",err,error,*999)

    CALL Tree_AssertIsIntgData(tree,err,error,*999)
    
    IF(.NOT.ASSOCIATED(x,tree%nil)) THEN
      CALL Tree_DetachInOrderIntg2(tree,x%left,count,treeValues,err,error,*999)
      count=count+1
      IF(count<=SIZE(treeValues,1)) THEN
        IF(.NOT.ALLOCATED(x%valuesIntg)) THEN
          localError="Tree node values is not allocated for key "//TRIM(NumberToVString(x%key,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(SIZE(x%valuesIntg,1)>SIZE(treeValues,2)) THEN
          localError="The second dimension of the supplied tree values array of "// &
            & TRIM(NumberToVString(SIZE(treeValues,2),"*",err,error))// &
            & " is too small. The size of the second dimension must be >= "// &
            & TRIM(NumberToVString(SIZE(x%valuesIntg,1),"*",err,error))// &
            & " for key number "//TRIM(NumberToVString(x%key,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        DO valueIdx=1,SIZE(x%valuesIntg,1)
          treeValues(count,valueIdx)=x%valuesIntg(valueIdx)
        ENDDO !valueIdx
      ELSE
        localError="The current count of the tree values of "//TRIM(NumberToVString(count,"*",err,error))// &
          & " is greater than the size of the tree values array of "// &
          & TRIM(NumberToVString(SIZE(treeValues,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      CALL Tree_DetachInOrderIntg2(tree,x%right,count,treeValues,err,error,*999)
    ENDIF

    EXITS("Tree_DetachInOrderIntg2")
    RETURN
999 ERRORSEXITS("Tree_DetachInOrderIntg2",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_DetachInOrderIntg2

  !
  !================================================================================================================================
  !

  !>Detaches the double precision real tree values in order from the specified tree node and adds them to the tree values array
  RECURSIVE SUBROUTINE Tree_DetachInOrderDP1(tree,x,count,treeValues,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree to detach
    TYPE(TreeNodeType), POINTER :: x !<A pointer to the specified tree node to detach from
    INTEGER(INTG), INTENT(INOUT) :: count !<The current number in the detached tree values array
    REAL(DP), INTENT(INOUT) :: treeValues(:) !<The current detached tree values array
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Tree_DetachInOrderDP1",err,error,*999)

    CALL Tree_AssertIsDPData(tree,err,error,*999)
    
    IF(.NOT.ASSOCIATED(x,tree%nil)) THEN
      CALL Tree_DetachInOrderDP1(tree,x%left,count,treeValues,err,error,*999)
      count=count+1
      IF(count<=SIZE(treeValues,1)) THEN
        IF(.NOT.ALLOCATED(x%valuesDP)) THEN
          localError="Tree node values is not allocated for key "//TRIM(NumberToVString(x%key,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        treeValues(count)=x%valuesDP(1)
      ELSE
        localError="The current count of the tree values of "//TRIM(NumberToVString(count,"*",err,error))// &
          & " is greater than the size of the tree values array of "// &
          & TRIM(NumberToVString(SIZE(treeValues,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      CALL Tree_DetachInOrderDP1(tree,x%right,count,treeValues,err,error,*999)
    ENDIF

    EXITS("Tree_DetachInOrderDP1")
    RETURN
999 ERRORSEXITS("Tree_DetachInOrderDP1",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_DetachInOrderDP1

  !
  !================================================================================================================================
  !

  !>Detaches the double precision real tree values in order from the specified tree node and adds them to the tree values array
  RECURSIVE SUBROUTINE Tree_DetachInOrderDP2(tree,x,count,treeValues,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree to detach
    TYPE(TreeNodeType), POINTER :: x !<A pointer to the specified tree node to detach from
    INTEGER(INTG), INTENT(INOUT) :: count !<The current number in the detached tree values array
    REAL(DP), INTENT(INOUT) :: treeValues(:,:) !<The current detached tree values array
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    INTEGER(INTG) :: valueIdx
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Tree_DetachInOrderDP2",err,error,*999)

    CALL Tree_AssertIsDPData(tree,err,error,*999)
    
    IF(.NOT.ASSOCIATED(x,tree%nil)) THEN
      CALL Tree_DetachInOrderDP2(tree,x%left,count,treeValues,err,error,*999)
      count=count+1
      IF(count<=SIZE(treeValues,1)) THEN
        IF(.NOT.ALLOCATED(x%valuesDP)) THEN
          localError="Tree node values is not allocated for key "//TRIM(NumberToVString(x%key,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(SIZE(x%valuesDP,1)>SIZE(treeValues,2)) THEN
          localError="The second dimension of the supplied tree values array of "// &
            & TRIM(NumberToVString(SIZE(treeValues,2),"*",err,error))// &
            & " is too small. The size of the second dimension must be >= "// &
            & TRIM(NumberToVString(SIZE(x%valuesDP,1),"*",err,error))// &
            & " for key number "//TRIM(NumberToVString(x%key,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        DO valueIdx=1,SIZE(x%valuesDP,1)
          treeValues(count,valueIdx)=x%valuesDP(valueIdx)
        ENDDO !valueIdx
      ELSE
        localError="The current count of the tree values of "//TRIM(NumberToVString(count,"*",err,error))// &
          & " is greater than the size of the tree values array of "// &
          & TRIM(NumberToVString(SIZE(treeValues,1),"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      CALL Tree_DetachInOrderDP2(tree,x%right,count,treeValues,err,error,*999)
    ENDIF

    EXITS("Tree_DetachInOrderDP2")
    RETURN
999 ERRORSEXITS("Tree_DetachInOrderDP2",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_DetachInOrderDP2

  !
  !================================================================================================================================
  !

  !>Finalises a tree and deallocates all memory
  SUBROUTINE Tree_Finalise(tree,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    
    ENTERS("Tree_Finalise",err,error,*999)

    IF(ASSOCIATED(tree)) THEN
      CALL Tree_NodeFinalise(tree,tree%root,err,error,*999)
      IF(ASSOCIATED(tree%nil)) DEALLOCATE(tree%nil)
      DEALLOCATE(tree)
    ENDIF

    EXITS("Tree_Finalise")
    RETURN
999 ERRORSEXITS("Tree_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises a tree
  SUBROUTINE Tree_Initialise(tree,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    
    ENTERS("Tree_Initialise",err,error,*999)

    IF(.NOT.ASSOCIATED(tree)) CALL FlagError("Tree is not associated.",err,error,*999)
    
    tree%treeFinished=.FALSE.
    tree%insertType=0
    tree%numberInTree=0
    tree%maximumNumberOfValues=0
    tree%dataType=TREE_INTG_TYPE
    tree%valueInitialiseIntg=0
    tree%valueInitialiseDP=0.0_DP
    NULLIFY(tree%root)
    NULLIFY(tree%nil)

    EXITS("Tree_Initialise")
    RETURN
999 ERRORSEXITS("Tree_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_Initialise

  !
  !================================================================================================================================
  !

  !>Gets the insert type for a tree
  SUBROUTINE Tree_InsertTypeGet(tree,insertType,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree
    INTEGER(INTG), INTENT(OUT) :: insertType !<On return, the tree insert type \see Trees_TreeInsertTypes,Trees::TreeInsertTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Tree_InsertTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(tree)) CALL FlagError("Tree is not associated.",err,error,*999)

    insertType=tree%insertType

    EXITS("Tree_InsertTypeGet")
    RETURN
999 ERRORSEXITS("Tree_InsertTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_InsertTypeGet

  !
  !================================================================================================================================
  !

  !>Sets/changes the insert type for a tree
  SUBROUTINE Tree_InsertTypeSet(tree,insertType,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree
    INTEGER(INTG), INTENT(IN) :: insertType !<The insert type to set \see Trees_TreeInsertTypes,Trees::TreeInsertTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Tree_InsertTypeSet",err,error,*999)

    CALL Tree_AssertNotFinished(tree,err,error,*999)
    
    SELECT CASE(insertType)
    CASE(TREE_DUPLICATES_ALLOWED_TYPE)
      tree%insertType=TREE_DUPLICATES_ALLOWED_TYPE
    CASE(TREE_NO_DUPLICATES_ALLOWED)
      tree%insertType=TREE_NO_DUPLICATES_ALLOWED
    CASE DEFAULT
      localError="The specified insert type of "//TRIM(NumberToVString(insertType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Tree_InsertTypeSet")
    RETURN
999 ERRORSEXITS("Tree_InsertTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_InsertTypeSet

  !
  !================================================================================================================================
  !

  !>Deletes a tree node specified by a key from a tree 
  SUBROUTINE Tree_ItemDelete(tree,key,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the Red-Black tree to delete from
    INTEGER(INTG), INTENT(IN) :: key !<A pointer to the tree node to delete
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    INTEGER(INTG) :: compareValue
    TYPE(TreeNodeType), POINTER :: u,v,w,x,y,z
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Tree_ItemDelete",err,error,*999)

    CALL Tree_AssertIsFinished(tree,err,error,*999)
    
    !Try and find the key to delete
    z=>tree%root
    IF(.NOT.ASSOCIATED(z,tree%nil)) THEN
      compareValue=z%key-key
      DO WHILE(compareValue/=0)
        IF(compareValue>0) THEN !z%key > key
          z=>z%left
        ELSE !z%key < key
          z=>z%right
        ENDIF
        IF(ASSOCIATED(z,tree%nil)) THEN
          EXIT
        ELSE
          compareValue=z%key-key
        ENDIF
      ENDDO
      IF(compareValue==0) THEN
        !Found the key so delete it
        IF(ASSOCIATED(z%left,tree%nil).OR.ASSOCIATED(z%right,tree%nil)) THEN
          y=>z
        ELSE
          y=>Tree_Successor(tree,z,err,error)
          IF(err/=0) GOTO 999
        ENDIF
        IF(.NOT.ASSOCIATED(y%left,tree%nil)) THEN
          x=>y%left
        ELSE
          x=>y%right
        ENDIF
        x%parent=>y%parent
        IF(ASSOCIATED(y%parent,tree%nil)) THEN
          tree%root=>x
        ELSE
          IF(ASSOCIATED(y,y%parent%left)) THEN
            y%parent%left=>x
          ELSE
            y%parent%right=>x
          ENDIF
        ENDIF
        IF(y%colour==TREE_BLACK_NODE) THEN
          !Fixup the delete to ensure the tree has red black properties
          !Note: Due to Fortran restrictions on aliasing pointers in dummy arguments we need to do the fixup and rotations
          !inside this routine rather than call fixup and rotate left and rotate right subroutines.
          DO WHILE(.NOT.ASSOCIATED(x,tree%root).AND.x%colour==TREE_BLACK_NODE)
            IF(ASSOCIATED(x,x%parent%left)) THEN
              w=>x%parent%right
              IF(w%colour==TREE_RED_NODE) THEN
                w%colour=TREE_BLACK_NODE
                x%parent%colour=TREE_RED_NODE
                !Rotate left on x->parent
                u=>x%parent
                v=>u%right
                u%right=>v%left
                IF(.NOT.ASSOCIATED(v%left,tree%nil)) v%left%parent=>u
                v%parent=>u%parent
                IF(ASSOCIATED(u%parent,tree%nil)) THEN
                  tree%root=>v
                ELSE
                  IF(ASSOCIATED(u,u%parent%left)) THEN
                    u%parent%left=>v
                  ELSE
                    u%parent%right=>v
                  ENDIF
                ENDIF
                v%left=>u
                u%parent=>v
                w=>x%parent%right
              ENDIF
              IF(w%left%colour==TREE_BLACK_NODE.AND.w%right%colour==TREE_BLACK_NODE) THEN
                w%colour=TREE_RED_NODE
                x=>x%parent
              ELSE
                IF(w%right%colour==TREE_BLACK_NODE) THEN
                  w%left%colour=TREE_BLACK_NODE
                  w%colour=TREE_RED_NODE
                  !Rotate right on w
                  u=>w
                  v=>u%left
                  u%left=>v%right
                  IF(.NOT.ASSOCIATED(v%right,tree%nil)) v%right%parent=>u
                  v%parent=>u%parent
                  IF(ASSOCIATED(v%parent,tree%nil)) THEN
                    tree%root=>v
                  ELSE
                    IF(ASSOCIATED(u,u%parent%right)) THEN
                      u%parent%right=>v
                    ELSE
                      u%parent%left=>v
                    ENDIF
                  ENDIF
                  v%right=>u
                  u%parent=>v
                  w=>x%parent%right
                ENDIF
                w%colour=x%parent%colour
                x%parent%colour=TREE_BLACK_NODE
                w%right%colour=TREE_BLACK_NODE
                !Rotate left on x->parent
                u=>x%parent
                v=>u%right
                u%right=>v%left
                IF(.NOT.ASSOCIATED(v%left,tree%nil)) v%left%parent=>u
                v%parent=>u%parent
                IF(ASSOCIATED(u%parent,tree%nil)) THEN
                  tree%root=>v
                ELSE
                  IF(ASSOCIATED(u,u%parent%left)) THEN
                    u%parent%left=>v
                  ELSE
                    u%parent%right=>v
                  ENDIF
                ENDIF
                v%left=>u
                u%parent=>v
                x=>tree%root
              ENDIF
            ELSE
              w=>x%parent%left
              IF(w%colour==TREE_RED_NODE) THEN
                w%colour=TREE_BLACK_NODE
                x%parent%colour=TREE_RED_NODE
                !Rotate right on x->parent
                u=>x%parent
                v=>u%left
                u%left=>v%right
                IF(.NOT.ASSOCIATED(v%right,tree%nil)) v%right%parent=>u
                v%parent=>u%parent
                IF(ASSOCIATED(v%parent,tree%nil)) THEN
                  tree%root=>v
                ELSE
                  IF(ASSOCIATED(u,u%parent%right)) THEN
                    u%parent%right=>v
                  ELSE
                    u%parent%left=>v
                  ENDIF
                ENDIF
                v%right=>u
                u%parent=>v
                w=>x%parent%left
              ENDIF
              IF(w%right%colour==TREE_BLACK_NODE.AND.w%left%colour==TREE_BLACK_NODE) THEN
                w%colour=TREE_RED_NODE
                x=>x%parent
              ELSE
                IF(w%left%colour==TREE_BLACK_NODE) THEN
                  w%right%colour=TREE_BLACK_NODE
                  w%colour=TREE_RED_NODE
                  !Rotate left on w
                  u=>w
                  v=>u%right
                  u%right=>v%left
                  IF(.NOT.ASSOCIATED(v%left,tree%nil)) v%left%parent=>U
                  v%parent=>u%parent
                  IF(ASSOCIATED(u%parent,tree%nil)) THEN
                    tree%root=>v
                  ELSE
                    IF(ASSOCIATED(u,u%parent%left)) THEN
                      u%parent%left=>v
                    ELSE
                      u%parent%right=>v
                    ENDIF
                  ENDIF
                  v%left=>u
                  u%parent=>v
                  w=>x%parent%left
                ENDIF
                w%colour=x%parent%colour
                x%parent%colour=TREE_BLACK_NODE
                w%left%colour=TREE_BLACK_NODE
                !Rotate right on x->parent
                u=>x%parent
                v=>u%left
                u%left=>v%right
                IF(.NOT.ASSOCIATED(v%right,tree%nil)) v%right%parent=>u
                v%parent=>u%parent
                IF(ASSOCIATED(v%parent,tree%nil)) THEN
                  tree%root=>v
                ELSE
                  IF(ASSOCIATED(U,u%parent%right)) THEN
                    u%parent%right=>v
                  ELSE
                    u%parent%left=>v
                  ENDIF
                ENDIF
                v%right=>u
                u%parent=>v
                x=>tree%root
              ENDIF
            ENDIF
          ENDDO
          x%colour=TREE_BLACK_NODE
        ENDIF
        IF(.NOT.ASSOCIATED(Y,Z)) THEN
          y%left=>z%left
          y%right=>z%right
          y%parent=>z%parent
          y%colour=z%colour
          z%left%parent=>y
          z%right%parent=>y
          IF(ASSOCIATED(Z,z%parent%left)) THEN
            z%parent%left=>y              
          ELSE
            z%parent%right=>y
          ENDIF
        ENDIF
        DEALLOCATE(z)
        tree%numberInTree=tree%numberInTree-1
      ELSE
        localError="Could not find the key "//TRIM(NumberToVString(key,"*",err,error))//" in the tree."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("The tree root is NIL. Can not delete the key.",err,error,*999)
    ENDIF

    EXITS("Tree_ItemDelete")
    RETURN
999 ERRORSEXITS("Tree_ItemDelete",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_ItemDelete

  !
  !================================================================================================================================
  !

  !>Inserts an integer tree node into a red-black tree 
  SUBROUTINE Tree_ItemInsertIntg0(tree,key,VALUE,insertStatus,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the Red-Black tree to insert into
    INTEGER(INTG), INTENT(IN) :: key !<The key to insert
    INTEGER(INTG), INTENT(IN) :: value !<The value to insert
    INTEGER(INTG), INTENT(OUT) :: insertStatus !<On exit, the status of the insert \see Trees_TreeNodeInsertStatus,Trees::TreeNodeInsertStatus
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    
    ENTERS("Tree_ItemInsertIntg0",err,error,*999)

    CALL Tree_ItemInsertIntg1(tree,key,[value],insertStatus,err,error,*999)

    EXITS("Tree_ItemInsertIntg0")
    RETURN
999 ERRORSEXITS("Tree_ItemInsertIntg0",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_ItemInsertIntg0
  
  !
  !================================================================================================================================
  !

  !>Inserts an integer tree node into a red-black tree 
  SUBROUTINE Tree_ItemInsertIntg1(tree,key,values,insertStatus,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the Red-Black tree to insert into
    INTEGER(INTG), INTENT(IN) :: key !<The key to insert
    INTEGER(INTG), INTENT(IN) :: values(:) !<The value to insert
    INTEGER(INTG), INTENT(OUT) :: insertStatus !<On exit, the status of the insert \see Trees_TreeNodeInsertStatus,Trees::TreeNodeInsertStatus
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    INTEGER(INTG), ALLOCATABLE :: newValues(:)
    LOGICAL :: duplicateKey
    TYPE(TreeNodeType), POINTER :: newTreeNode,x,y,z

    NULLIFY(newTreeNode)
    
    ENTERS("Tree_ItemInsertIntg1",err,error,*998)

    CALL Tree_AssertIsFinished(tree,err,error,*999)
    CALL Tree_AssertIsIntgData(tree,err,error,*999)
    
    !Find the position to insert
    y=>tree%nil
    x=>tree%root
    duplicateKey=.FALSE.
    DO WHILE(.NOT.ASSOCIATED(x,tree%nil))
      y=>x
      duplicateKey=(key==x%key)
      IF(duplicateKey) THEN
        EXIT
      ELSE IF(key<x%key) THEN
        x=>x%left
      ELSE
        x=>x%right
      ENDIF
    ENDDO
    IF(duplicateKey) THEN
      IF(tree%insertType==TREE_NO_DUPLICATES_ALLOWED) THEN
        insertStatus=TREE_NODE_DUPLICATE_KEY
      ELSE
        !Append the values to the current values
        ALLOCATE(newValues(SIZE(x%valuesIntg,1)+SIZE(values,1)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new values.",err,error,*999)
        newValues(1:SIZE(x%valuesIntg,1))=x%valuesIntg(1:SIZE(x%valuesIntg,1))
        newValues(SIZE(x%valuesIntg,1)+1:SIZE(x%valuesIntg,1)+SIZE(values,1))=values(1:SIZE(values,1))
        CALL MOVE_ALLOC(newValues,x%valuesIntg)
        IF((SIZE(x%valuesIntg,1)+SIZE(values,1))>tree%maximumNumberOfValues) &
          & tree%maximumNumberOfValues=SIZE(x%valuesIntg,1)+SIZE(values,1)
        insertStatus=TREE_NODE_INSERT_SUCESSFUL
      ENDIF
    ELSE
      !Allocate the new tree node and set its key and value
      ALLOCATE(newTreeNode,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new tree node.",err,error,*999)
      CALL Tree_NodeInitialise(tree,newTreeNode,err,error,*999)
      newTreeNode%key=key
      ALLOCATE(newTreeNode%valuesIntg(SIZE(values,1)),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new tree node values.",err,error,*999)      
      newTreeNode%valuesIntg(1:SIZE(values,1))=values(1:SIZE(values,1))
      IF(SIZE(values,1)>tree%maximumNumberOfValues) tree%maximumNumberOfValues=SIZE(values,1)
      !Insert the new tree node into the tree
      newTreeNode%colour=TREE_RED_NODE
      newTreeNode%left=>tree%nil
      newTreeNode%right=>tree%nil
      newTreeNode%parent=>y
      IF(ASSOCIATED(y,tree%nil)) THEN
        tree%root=>newTreeNode
      ELSE
        IF(newTreeNode%key<y%key) THEN
          y%left=>newTreeNode
        ELSE
          y%right=>newTreeNode
        ENDIF
      ENDIF
      !Fix up the tree to keep red-black properties
      !Note: Due to Fortran restrictions on aliasing pointers in dummy arguments we need to do the fixup and rotations
      !inside this routine rather than call fixup and rotate left and rotate right subroutines.
      z=>newTreeNode
      DO WHILE(z%parent%colour==TREE_RED_NODE)
        IF(ASSOCIATED(z%parent,z%parent%parent%left)) THEN
          y=>z%parent%parent%right
          IF(y%colour==TREE_RED_NODE) THEN
            z%parent%colour=TREE_BLACK_NODE
            y%colour=TREE_BLACK_NODE
            z%parent%parent%colour=TREE_RED_NODE
            z=>z%parent%parent
          ELSE
            IF(ASSOCIATED(z,z%parent%right)) THEN
              z=>z%parent
              !Rotate the tree left at z
              x=>z                  
              y=>x%right
              x%right=>y%left
              IF(.NOT.ASSOCIATED(y%left,tree%nil)) y%left%parent=>x
              y%parent=>x%parent
              IF(ASSOCIATED(x%parent,tree%nil)) THEN
                tree%root=>y
              ELSE
                IF(ASSOCIATED(x,x%parent%left)) THEN
                  x%parent%left=>y
                ELSE
                  x%parent%right=>y
                ENDIF
              ENDIF
              y%left=>x
              x%parent=>y
            ENDIF
            z%parent%colour=TREE_BLACK_NODE
            z%parent%parent%colour=TREE_RED_NODE
            !Rotate the tree right at z->parent->parent
            x=>z%parent%parent
            y=>x%left
            x%left=>y%right
            IF(.NOT.ASSOCIATED(y%right,tree%nil)) y%right%parent=>x
            y%parent=>x%parent
            IF(ASSOCIATED(x%parent,tree%nil)) THEN
              tree%root=>y
            ELSE
              IF(ASSOCIATED(x,x%parent%right)) THEN
                x%parent%right=>y
              ELSE
                x%parent%left=>y
              ENDIF
            ENDIF
            y%right=>x
            x%parent=>y
          ENDIF
        ELSE
          y=>z%parent%parent%left
          IF(y%colour==TREE_RED_NODE) THEN
            z%parent%colour=TREE_BLACK_NODE
            y%colour=TREE_BLACK_NODE
            z%parent%parent%colour=TREE_RED_NODE
            z=>z%parent%parent
          ELSE
            IF(ASSOCIATED(z,z%parent%left)) THEN
              z=>z%parent
              x=>z
              !Rotate the tree right at z
              y=>x%left
              x%left=>y%right
              IF(.NOT.ASSOCIATED(y%right,tree%nil)) y%right%parent=>x
              y%parent=>x%parent
              IF(ASSOCIATED(x%parent,tree%nil)) THEN
                tree%root=>y
              ELSE
                IF(ASSOCIATED(x,x%parent%right)) THEN
                  x%parent%right=>y
                ELSE
                  x%parent%left=>y
                ENDIF
              ENDIF
              y%right=>x
              x%parent=>y
            ENDIF
            z%parent%colour=TREE_BLACK_NODE
            z%parent%parent%colour=TREE_RED_NODE
            !Rotate the tree left at z->parent->parent
            x=>z%parent%parent
            y=>x%right
            x%right=>y%left
            IF(.NOT.ASSOCIATED(y%left,tree%nil)) y%left%parent=>x
            y%parent=>x%parent
            IF(ASSOCIATED(x%parent,tree%nil)) THEN
              tree%root=>y
            ELSE
              IF(ASSOCIATED(X,x%parent%left)) THEN
                x%parent%left=>y
              ELSE
                x%parent%right=>y
              ENDIF
            ENDIF
            y%left=>x
            x%parent=>y
          ENDIF
        ENDIF
      ENDDO
      tree%root%colour=TREE_BLACK_NODE
      !Increment the number in the tree and indicate a successful insertion
      tree%numberInTree=tree%numberInTree+1
      insertStatus=TREE_NODE_INSERT_SUCESSFUL
    ENDIF

    EXITS("Tree_ItemInsertIntg1")
    RETURN
999 IF(ASSOCIATED(newTreeNode)) DEALLOCATE(newTreeNode)
    IF(ALLOCATED(newValues)) DEALLOCATE(newValues)
998 ERRORSEXITS("Tree_ItemInsertIntg1",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_ItemInsertIntg1

  !
  !================================================================================================================================
  !

  !>Inserts a double precision real tree node into a red-black tree 
  SUBROUTINE Tree_ItemInsertDP0(tree,key,VALUE,insertStatus,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the Red-Black tree to insert into
    INTEGER(INTG), INTENT(IN) :: key !<The key to insert
    REAL(DP), INTENT(IN) :: value !<The value to insert
    INTEGER(INTG), INTENT(OUT) :: insertStatus !<On exit, the status of the insert \see Trees_TreeNodeInsertStatus,Trees::TreeNodeInsertStatus
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    
    ENTERS("Tree_ItemInsertDP0",err,error,*999)

    CALL Tree_ItemInsertDP1(tree,key,[value],insertStatus,err,error,*999)

    EXITS("Tree_ItemInsertDP0")
    RETURN
999 ERRORSEXITS("Tree_ItemInsertDP0",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_ItemInsertDP0
  
  !
  !================================================================================================================================
  !

  !>Inserts a double precision real tree node into a red-black tree 
  SUBROUTINE Tree_ItemInsertDP1(tree,key,values,insertStatus,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the Red-Black tree to insert into
    INTEGER(INTG), INTENT(IN) :: key !<The key to insert
    REAL(DP), INTENT(IN) :: values(:) !<The value to insert
    INTEGER(INTG), INTENT(OUT) :: insertStatus !<On exit, the status of the insert \see Trees_TreeNodeInsertStatus,Trees::TreeNodeInsertStatus
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    REAL(DP), ALLOCATABLE :: newValues(:)
    LOGICAL :: duplicateKey
    TYPE(TreeNodeType), POINTER :: newTreeNode,x,y,z

    NULLIFY(newTreeNode)
    
    ENTERS("Tree_ItemInsertDP1",err,error,*998)

    CALL Tree_AssertIsFinished(tree,err,error,*999)
    CALL Tree_AssertIsDPData(tree,err,error,*999)
    
    !Find the position to insert
    y=>tree%nil
    x=>tree%root
    duplicateKey=.FALSE.
    DO WHILE(.NOT.ASSOCIATED(x,tree%nil))
      y=>x
      duplicateKey=(key==x%key)
      IF(duplicateKey) THEN
        EXIT
      ELSE IF(key<x%key) THEN
        x=>x%left
      ELSE
        x=>x%right
      ENDIF
    ENDDO
    IF(duplicateKey) THEN
      IF(tree%insertType==TREE_NO_DUPLICATES_ALLOWED) THEN
        insertStatus=TREE_NODE_DUPLICATE_KEY
      ELSE
        !Append the values to the current values
        ALLOCATE(newValues(SIZE(x%valuesDP,1)+SIZE(values,1)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new values.",err,error,*999)
        newValues(1:SIZE(x%valuesDP,1))=x%valuesDP(1:SIZE(x%valuesDP,1))
        newValues(SIZE(x%valuesDP,1)+1:SIZE(x%valuesDP,1)+SIZE(values,1))=values(1:SIZE(values,1))
        CALL MOVE_ALLOC(newValues,x%valuesDP)
        IF((SIZE(x%valuesDP,1)+SIZE(values,1))>tree%maximumNumberOfValues) &
          & tree%maximumNumberOfValues=SIZE(x%valuesDP,1)+SIZE(values,1)
        insertStatus=TREE_NODE_INSERT_SUCESSFUL
      ENDIF
    ELSE
      !Allocate the new tree node and set its key and value
      ALLOCATE(newTreeNode,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new tree node.",err,error,*999)
      CALL Tree_NodeInitialise(tree,newTreeNode,err,error,*999)
      newTreeNode%key=key
      ALLOCATE(newTreeNode%valuesDP(SIZE(values,1)),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new tree node values.",err,error,*999)      
      newTreeNode%valuesDP(1:SIZE(values,1))=values(1:SIZE(values,1))
      IF(SIZE(values,1)>tree%maximumNumberOfValues) tree%maximumNumberOfValues=SIZE(values,1)
      !Insert the new tree node into the tree
      newTreeNode%colour=TREE_RED_NODE
      newTreeNode%left=>tree%nil
      newTreeNode%right=>tree%nil
      newTreeNode%parent=>y
      IF(ASSOCIATED(y,tree%nil)) THEN
        tree%root=>newTreeNode
      ELSE
        IF(newTreeNode%key<y%key) THEN
          y%left=>newTreeNode
        ELSE
          y%right=>newTreeNode
        ENDIF
      ENDIF
      !Fix up the tree to keep red-black properties
      !Note: Due to Fortran restrictions on aliasing pointers in dummy arguments we need to do the fixup and rotations
      !inside this routine rather than call fixup and rotate left and rotate right subroutines.
      z=>newTreeNode
      DO WHILE(z%parent%colour==TREE_RED_NODE)
        IF(ASSOCIATED(z%parent,z%parent%parent%left)) THEN
          y=>z%parent%parent%right
          IF(y%colour==TREE_RED_NODE) THEN
            z%parent%colour=TREE_BLACK_NODE
            y%colour=TREE_BLACK_NODE
            z%parent%parent%colour=TREE_RED_NODE
            z=>z%parent%parent
          ELSE
            IF(ASSOCIATED(z,z%parent%right)) THEN
              z=>z%parent
              !Rotate the tree left at z
              x=>z                  
              y=>x%right
              x%right=>y%left
              IF(.NOT.ASSOCIATED(y%left,tree%nil)) y%left%parent=>x
              y%parent=>x%parent
              IF(ASSOCIATED(x%parent,tree%nil)) THEN
                tree%root=>y
              ELSE
                IF(ASSOCIATED(x,x%parent%left)) THEN
                  x%parent%left=>y
                ELSE
                  x%parent%right=>y
                ENDIF
              ENDIF
              y%left=>x
              x%parent=>y
            ENDIF
            z%parent%colour=TREE_BLACK_NODE
            z%parent%parent%colour=TREE_RED_NODE
            !Rotate the tree right at z->parent->parent
            x=>z%parent%parent
            y=>x%left
            x%left=>y%right
            IF(.NOT.ASSOCIATED(y%right,tree%nil)) y%right%parent=>x
            y%parent=>x%parent
            IF(ASSOCIATED(x%parent,tree%nil)) THEN
              tree%root=>y
            ELSE
              IF(ASSOCIATED(x,x%parent%right)) THEN
                x%parent%right=>y
              ELSE
                x%parent%left=>y
              ENDIF
            ENDIF
            y%right=>x
            x%parent=>y
          ENDIF
        ELSE
          y=>z%parent%parent%left
          IF(y%colour==TREE_RED_NODE) THEN
            z%parent%colour=TREE_BLACK_NODE
            y%colour=TREE_BLACK_NODE
            z%parent%parent%colour=TREE_RED_NODE
            z=>z%parent%parent
          ELSE
            IF(ASSOCIATED(z,z%parent%left)) THEN
              z=>z%parent
              x=>z
              !Rotate the tree right at z
              y=>x%left
              x%left=>y%right
              IF(.NOT.ASSOCIATED(y%right,tree%nil)) y%right%parent=>x
              y%parent=>x%parent
              IF(ASSOCIATED(x%parent,tree%nil)) THEN
                tree%root=>y
              ELSE
                IF(ASSOCIATED(x,x%parent%right)) THEN
                  x%parent%right=>y
                ELSE
                  x%parent%left=>y
                ENDIF
              ENDIF
              y%right=>x
              x%parent=>y
            ENDIF
            z%parent%colour=TREE_BLACK_NODE
            z%parent%parent%colour=TREE_RED_NODE
            !Rotate the tree left at z->parent->parent
            x=>z%parent%parent
            y=>x%right
            x%right=>y%left
            IF(.NOT.ASSOCIATED(y%left,tree%nil)) y%left%parent=>x
            y%parent=>x%parent
            IF(ASSOCIATED(x%parent,tree%nil)) THEN
              tree%root=>y
            ELSE
              IF(ASSOCIATED(X,x%parent%left)) THEN
                x%parent%left=>y
              ELSE
                x%parent%right=>y
              ENDIF
            ENDIF
            y%left=>x
            x%parent=>y
          ENDIF
        ENDIF
      ENDDO
      tree%root%colour=TREE_BLACK_NODE
      !Increment the number in the tree and indicate a successful insertion
      tree%numberInTree=tree%numberInTree+1
      insertStatus=TREE_NODE_INSERT_SUCESSFUL
    ENDIF

    EXITS("Tree_ItemInsertDP1")
    RETURN
999 IF(ASSOCIATED(newTreeNode)) DEALLOCATE(newTreeNode)
    IF(ALLOCATED(newValues)) DEALLOCATE(newValues)
998 ERRORSEXITS("Tree_ItemInsertDP1",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_ItemInsertDP1

  !
  !================================================================================================================================
  !

  !>Finalises a tree node and deallocates all memory
  RECURSIVE SUBROUTINE Tree_NodeFinalise(tree,treeNode,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree containing the tree node to finalise
    TYPE(TreeNodeType), POINTER :: treeNode !<A pointer to the tree node to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    
    ENTERS("Tree_NodeFinalise",err,error,*999)

    IF(.NOT.ASSOCIATED(tree)) CALL FlagError("Tree is not associated.",err,error,*999)
    
    IF(.NOT.ASSOCIATED(treeNode,tree%nil)) THEN
      CALL Tree_NodeFinalise(tree,treeNode%left,err,error,*999)
      CALL Tree_NodeFinalise(tree,treeNode%right,err,error,*999)
      IF(ALLOCATED(treeNode%valuesIntg)) DEALLOCATE(treeNode%valuesIntg)
      IF(ALLOCATED(treeNode%valuesDP)) DEALLOCATE(treeNode%valuesDP)
      DEALLOCATE(treeNode)
    ENDIF

    EXITS("Tree_NodeFinalise")
    RETURN
999 ERRORSEXITS("Tree_NodeFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_NodeFinalise

  !
  !================================================================================================================================
  !

  !>Initialises a tree node
  SUBROUTINE Tree_NodeInitialise(tree,treeNode,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree containing the tree node to initialise
    TYPE(TreeNodeType), POINTER :: treeNode !<A pointer to the tree node to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    
    ENTERS("Tree_NodeInitialise",err,error,*999)

    IF(.NOT.ASSOCIATED(tree)) CALL FlagError("Tree is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(treeNode)) CALL FlagError("Tree node is not associated.",err,error,*999)
    
    treeNode%key=0
    treeNode%colour=TREE_BLACK_NODE
    NULLIFY(treeNode%left)
    NULLIFY(treeNode%right)
    NULLIFY(treeNode%parent)

    EXITS("Tree_NodeInitialise")
    RETURN
999 ERRORSEXITS("Tree_NodeInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_NodeInitialise

  !
  !================================================================================================================================
  !

  !>Gets the key at a specified tree node
  SUBROUTINE Tree_NodeKeyGet(tree,treeNode,key,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree containing the tree node 
    TYPE(TreeNodeType), POINTER :: treeNode !<A pointer to the tree node to get the key of
    INTEGER(INTG), INTENT(OUT) :: key !<On exit, the key of the specified tree node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    
    ENTERS("Tree_NodeKeyGet",err,error,*999)

    CALL Tree_AssertIsFinished(tree,err,error,*999)
    IF(.NOT.ASSOCIATED(treeNode)) CALL FlagError("Tree node is not associated.",err,error,*999)
    
    key=treeNode%key

    EXITS("Tree_NodeKeyGet")
    RETURN
999 ERRORSEXITS("Tree_NodeKeyGet",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_NodeKeyGet

  !
  !================================================================================================================================
  !

  !>Gets the integer value at a specified tree node
  SUBROUTINE Tree_NodeValueGetIntg0(tree,treeNode,VALUE,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree containing the tree node 
    TYPE(TreeNodeType), POINTER :: treeNode !<A pointer to the tree node to get the value of
    INTEGER(INTG), INTENT(OUT) :: value !<On exit, the value of the specified tree node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    INTEGER(INTG) :: values(1)
    
    ENTERS("Tree_NodeValueGetIntg0",err,error,*999)

    CALL Tree_NodeValueGetIntg1(tree,treeNode,values,err,error,*999)
    value=values(1)

    EXITS("Tree_NodeValueGetIntg0")
    RETURN
999 ERRORSEXITS("Tree_NodeValueGetIntg0",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_NodeValueGetIntg0

  !
  !================================================================================================================================
  !

  !>Gets the integer values at a specified tree node
  SUBROUTINE Tree_NodeValueGetIntg1(tree,treeNode,values,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree containing the tree node 
    TYPE(TreeNodeType), POINTER :: treeNode !<A pointer to the tree node to get the value of
    INTEGER(INTG), INTENT(OUT) :: values(:) !<On exit, the value of the specified tree node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Tree_NodeValueGetIntg1",err,error,*999)

    CALL Tree_AssertIsFinished(tree,err,error,*999)
    CALL Tree_AssertIsIntgData(tree,err,error,*999)
    IF(.NOT.ASSOCIATED(treeNode)) CALL FlagError("Tree node is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(treeNode%valuesIntg)) THEN
      localError="The tree node values is not allocated for tree node key "// &
        & TRIM(NumberToVString(treeNode%key,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(values,1)<SIZE(treeNode%valuesIntg,1)) THEN
      localError="The size of the supplied values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))// &
        & " is too small. The size should be > "//TRIM(NumberToVString(SIZE(treeNode%valuesIntg,1),"*",err,error))// &
        & " for tree node key "//TRIM(NumberToVString(treeNode%key,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    values(1:SIZE(treeNode%valuesIntg,1))=treeNode%valuesIntg(1:SIZE(treeNode%valuesIntg,1))

    EXITS("Tree_NodeValueGetIntg1")
    RETURN
999 ERRORSEXITS("Tree_NodeValueGetIntg1",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_NodeValueGetIntg1

  !
  !================================================================================================================================
  !

  !>Gets the double precision real value at a specified tree node
  SUBROUTINE Tree_NodeValueGetDP0(tree,treeNode,VALUE,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree containing the tree node 
    TYPE(TreeNodeType), POINTER :: treeNode !<A pointer to the tree node to get the value of
    REAL(DP), INTENT(OUT) :: value !<On exit, the value of the specified tree node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    REAL(DP) :: values(1)
    
    ENTERS("Tree_NodeValueGetDP0",err,error,*999)

    CALL Tree_NodeValueGetDP1(tree,treeNode,values,err,error,*999)
    value=values(1)

    EXITS("Tree_NodeValueGetDP0")
    RETURN
999 ERRORSEXITS("Tree_NodeValueGetDP0",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_NodeValueGetDP0

  !
  !================================================================================================================================
  !

  !>Gets the double precision real values at a specified tree node
  SUBROUTINE Tree_NodeValueGetDP1(tree,treeNode,values,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree containing the tree node 
    TYPE(TreeNodeType), POINTER :: treeNode !<A pointer to the tree node to get the value of
    REAL(DP), INTENT(OUT) :: values(:) !<On exit, the value of the specified tree node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Tree_NodeValueGetDP1",err,error,*999)

    CALL Tree_AssertIsFinished(tree,err,error,*999)
    CALL Tree_AssertIsDPData(tree,err,error,*999)
    IF(.NOT.ASSOCIATED(treeNode)) CALL FlagError("Tree node is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(treeNode%valuesDP)) THEN
      localError="The tree node values is not allocated for tree node key "// &
        & TRIM(NumberToVString(treeNode%key,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(values,1)<SIZE(treeNode%valuesDP,1)) THEN
      localError="The size of the supplied values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))// &
        & " is too small. The size should be > "//TRIM(NumberToVString(SIZE(treeNode%valuesDP,1),"*",err,error))// &
        & " for tree node key "//TRIM(NumberToVString(treeNode%key,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    values(1:SIZE(treeNode%valuesDP,1))=treeNode%valuesDP(1:SIZE(treeNode%valuesDP,1))

    EXITS("Tree_NodeValueGetDP1")
    RETURN
999 ERRORSEXITS("Tree_NodeValueGetDP1",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_NodeValueGetDP1

  !
  !================================================================================================================================
  !

  !>Sets the integer value at a specified tree node
  SUBROUTINE Tree_NodeValueSetIntg0(tree,treeNode,VALUE,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree containing the tree node 
    TYPE(TreeNodeType), POINTER :: treeNode !<A pointer to the tree node to set the value of
    INTEGER(INTG), INTENT(IN) :: value !<The value of the specified tree node to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    
    ENTERS("Tree_NodeValueSetIntg0",err,error,*999)

    CALL Tree_NodeValueSetIntg1(tree,treeNode,[value],err,error,*999)
    
    EXITS("Tree_NodeValueSetIntg0")
    RETURN
999 ERRORSEXITS("Tree_NodeValueSetIntg0",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_NodeValueSetIntg0

  !
  !================================================================================================================================
  !

  !>Sets the integer values at a specified tree node
  SUBROUTINE Tree_NodeValueSetIntg1(tree,treeNode,values,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree containing the tree node 
    TYPE(TreeNodeType), POINTER :: treeNode !<A pointer to the tree node to set the value of
    INTEGER(INTG), INTENT(IN) :: values(:) !<The values of the specified tree node to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    INTEGER(INTG), ALLOCATABLE :: newValues(:)
    
    ENTERS("Tree_NodeValueSetIntg1",err,error,*999)

    CALL Tree_AssertIsFinished(tree,err,error,*999)
    CALL Tree_AssertIsIntgData(tree,err,error,*999)
    IF(.NOT.ASSOCIATED(treeNode)) CALL FlagError("Tree node is not associated.",err,error,*999)

    ALLOCATE(newValues(SIZE(values,1)),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new values.",err,error,*999)
    newValues(1:SIZE(values,1))=values(1:SIZE(values,1))
    CALL MOVE_ALLOC(newValues,treeNode%valuesIntg)
    IF(SIZE(values,1)>tree%maximumNumberOfValues) tree%maximumNumberOfValues=SIZE(values,1)

    EXITS("Tree_NodeValueSetIntg1")
    RETURN
999 IF(ALLOCATED(newValues)) DEALLOCATE(newValues)
    ERRORSEXITS("Tree_NodeValueSetIntg1",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_NodeValueSetIntg1

  !
  !================================================================================================================================
  !

  !>Sets the double precision real value at a specified tree node
  SUBROUTINE Tree_NodeValueSetDP0(tree,treeNode,VALUE,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree containing the tree node 
    TYPE(TreeNodeType), POINTER :: treeNode !<A pointer to the tree node to set the value of
    REAL(DP), INTENT(IN) :: value !<The value of the specified tree node to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    
    ENTERS("Tree_NodeValueSetDP0",err,error,*999)

    CALL Tree_NodeValueSetDP1(tree,treeNode,[value],err,error,*999)
    
    EXITS("Tree_NodeValueSetDP0")
    RETURN
999 ERRORSEXITS("Tree_NodeValueSetDP0",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_NodeValueSetDP0

  !
  !================================================================================================================================
  !

  !>Sets the double precision real values at a specified tree node
  SUBROUTINE Tree_NodeValueSetDP1(tree,treeNode,values,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree containing the tree node 
    TYPE(TreeNodeType), POINTER :: treeNode !<A pointer to the tree node to set the value of
    REAL(DP), INTENT(IN) :: values(:) !<The values of the specified tree node to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    REAL(DP), ALLOCATABLE :: newValues(:)
    
    ENTERS("Tree_NodeValueSetDP1",err,error,*999)

    CALL Tree_AssertIsFinished(tree,err,error,*999)
    CALL Tree_AssertIsDPData(tree,err,error,*999)
    IF(.NOT.ASSOCIATED(treeNode)) CALL FlagError("Tree node is not associated.",err,error,*999)

    ALLOCATE(newValues(SIZE(values,1)),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new values.",err,error,*999)
    newValues(1:SIZE(values,1))=values(1:SIZE(values,1))
    CALL MOVE_ALLOC(newValues,treeNode%valuesDP)
    IF(SIZE(values,1)>tree%maximumNumberOfValues) tree%maximumNumberOfValues=SIZE(values,1)

    EXITS("Tree_NodeValueSetDP1")
    RETURN
999 IF(ALLOCATED(newValues)) DEALLOCATE(newValues)
    ERRORSEXITS("Tree_NodeValueSetDP1",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_NodeValueSetDP1

  !
  !================================================================================================================================
  !

  !>Outputs a tree to the specified output stream id
  SUBROUTINE Tree_Output(id,tree,err,error,*)

    !Argument Variables
    INTEGER(INTG), INTENT(IN) :: id !<The id of the output stream
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree to search
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    
    ENTERS("Tree_Output",err,error,*999)

    CALL Tree_AssertIsFinished(tree,err,error,*999)
    
    CALL WriteString(id,"Tree:",err,error,*999)
    CALL WriteStringValue(id,"Number of tree nodes = ",tree%numberInTree,err,error,*999)
    CALL WriteStringValue(id,"Maximum number of values = ",tree%maximumNumberOfValues,err,error,*999)
    CALL WriteStringValue(id,"Tree insert type = ",tree%insertType,err,error,*999)
    CALL WriteStringValue(id,"Tree data type = ",tree%dataType,err,error,*999)
    CALL Tree_OutputInOrder(id,tree,tree%root,err,error,*999)

    EXITS("Tree_Output")
    RETURN
999 ERRORSEXITS("Tree_Output",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_Output

  !
  !================================================================================================================================
  !

  !>Outputs a tree in order to the specified output stream id from the specified tree node
  RECURSIVE SUBROUTINE Tree_OutputInOrder(id,tree,x,err,error,*)

    !Argument Variables
    INTEGER(INTG), INTENT(IN) :: id !<The id of the output stream
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree to search
    TYPE(TreeNodeType), POINTER :: x !<A pointer to the tree node to output from
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Tree_OutputInOrder",err,error,*999)

    IF(.NOT.ASSOCIATED(tree)) CALL FlagError("Tree is not associated.",err,error,*999)
    
    IF(.NOT.ASSOCIATED(x,tree%nil)) THEN
      !Output the left subtree first
      CALL Tree_OutputInOrder(id,tree,x%left,err,error,*999)
      !Now output the information for this node
      CALL WriteString(id,"  Tree Node:",err,error,*999)
      CALL WriteStringValue(id,"    Key = ",x%key,err,error,*999)
      SELECT CASE(tree%dataType)
      CASE(TREE_INTG_TYPE)
        IF(ALLOCATED(x%valuesIntg)) CALL WriteStringValue(id,"    Value = ",x%valuesIntg(1),err,error,*999)
      CASE(TREE_DP_TYPE)
        IF(ALLOCATED(x%valuesDP)) CALL WriteStringValue(id,"    Value = ",x%valuesDP(1),err,error,*999)
      CASE DEFAULT
        localError="The tree data type of "//TRIM(NumberToVString(tree%dataType,"*",err,error))//" is invalid or not implemented."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      CALL WriteStringValue(id,"    Colour = ",x%colour,err,error,*999)
      IF(ASSOCIATED(x%left,tree%nil)) THEN
        CALL WriteString(id,"    Left Key = NIL",err,error,*999)
      ELSE
        CALL WriteStringValue(id,"    Left Key = ",x%left%key,err,error,*999)
      ENDIF
      IF(ASSOCIATED(x%right,tree%nil)) THEN
        CALL WriteString(id,"    Right Key = NIL",err,error,*999)
      ELSE
        CALL WriteStringValue(id,"    Right Key = ",x%right%key,err,error,*999)
      ENDIF
      IF(ASSOCIATED(x%parent,tree%nil)) THEN
        CALL WriteString(id,"    Parent Key = NIL",err,error,*999)
      ELSE
        CALL WriteStringValue(id,"    Parent Key = ",x%parent%key,err,error,*999)
      ENDIF
      !Output the right subtree last
      CALL Tree_OutputInOrder(id,tree,x%right,err,error,*999)
    ENDIF

    EXITS("Tree_OutputInOrder")
    RETURN
999 ERRORSEXITS("Tree_OutputInOrder",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_OutputInOrder

  !
  !================================================================================================================================
  !

  !>Returns the predeccessor of a tree at a specified tree node
  FUNCTION Tree_Predecessor(tree,x,err,error)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the Red-Black tree to find the predecessor of
    TYPE(TreeNodeType), POINTER :: x !<A pointer to the tree node to return the predecessor of
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Function variable
    TYPE(TreeNodeType), POINTER :: Tree_Predecessor !<On return the pointer to the predecessor of X or NIL if no predecessor exits
    !Local Variables
    TYPE(TreeNodeType), POINTER :: y

    NULLIFY(Tree_Predecessor)
    
    ENTERS("Tree_Predecessor",err,error,*999)

    IF(.NOT.ASSOCIATED(tree)) CALL FlagError("Tree is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(x)) CALL FlagError("Tree node x is not associated.",err,error,*999)
    
    y=>x%left
    IF(ASSOCIATED(y,tree%nil)) THEN
      DO WHILE(.NOT.ASSOCIATED(y%right,tree%nil))
        y=>y%right
      ENDDO
      Tree_Predecessor=>y
    ELSE
      y=>x%parent
      DO WHILE(ASSOCIATED(x,y%left))
        IF(ASSOCIATED(y,tree%root)) THEN
          Tree_Predecessor=>tree%nil
          EXIT
        ELSE
          x=>y
          y=>y%parent
        ENDIF
      ENDDO
      IF(.NOT.ASSOCIATED(Tree_Predecessor)) Tree_Predecessor=>y
    ENDIF

    EXITS("Tree_Predecessor")
    RETURN
999 ERRORSEXITS("Tree_Predecessor",err,error)
    RETURN
    
  END FUNCTION Tree_Predecessor

  !
  !================================================================================================================================
  !

  !>Searches a tree to see if it contains a key
  SUBROUTINE Tree_Search(tree,key,x,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree to search
    INTEGER(INTG), INTENT(IN) :: key !<The key to search for
    TYPE(TreeNodeType), POINTER :: x !<On return a pointer to the tree node containing the key. If the key does not exist NULL is returned
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    INTEGER(INTG) :: compareValue
    TYPE(TreeNodeType), POINTER :: y
    
    ENTERS("Tree_Search",err,error,*999)

    CALL Tree_AssertIsFinished(tree,err,error,*999)
    IF(ASSOCIATED(x)) CALL FlagError("The tree node x is already associated.",err,error,*999)
    
    NULLIFY(x)
    y=>tree%root
    IF(.NOT.ASSOCIATED(y,tree%nil)) THEN
      compareValue=y%key-key
      DO WHILE(compareValue/=0)
        IF(compareValue>0) THEN !y%key > key
          y=>y%left
        ELSE !y%key < key
          y=>y%right
        ENDIF
        IF(ASSOCIATED(y,tree%nil)) THEN
          EXIT
        ELSE
          compareValue=y%key-key
        ENDIF
      ENDDO
      IF(compareValue==0) x=>y
    ENDIF

    EXITS("Tree_Search")
    RETURN
999 ERRORSEXITS("Tree_Search",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_Search

  !
  !================================================================================================================================
  !

  !>Returns the successor of a tree at a specified tree node
  FUNCTION Tree_Successor(tree,x,err,error)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the Red-Black tree to find the successor of
    TYPE(TreeNodeType), POINTER :: x !<A pointer to the tree node to return the successor of
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Function variable
    TYPE(TreeNodeType), POINTER :: tree_Successor !<On return the pointer to the successor of X or NIL if no sucessor exits
    !Local Variables
    TYPE(TreeNodeType), POINTER :: y

    NULLIFY(Tree_Successor)
    
    ENTERS("Tree_Successor",err,error,*999)

    IF(.NOT.ASSOCIATED(tree)) CALL FlagError("Tree is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(x)) CALL FlagError("Tree node x is not associated.",err,error,*999)
    
    y=>x%right
    IF(ASSOCIATED(y,tree%nil)) THEN
      DO WHILE(.NOT.ASSOCIATED(y%left,tree%nil))
        y=>y%left
      ENDDO
      Tree_Successor=>y
      EXITS("Tree_Successor")
      RETURN
    ELSE
      y=>x%parent
      DO WHILE(ASSOCIATED(x,y%right))
        x=>y
        y=>y%parent
      ENDDO
      IF(ASSOCIATED(y,tree%root)) THEN
        Tree_Successor=>tree%nil
      ELSE
        Tree_Successor=>y
      ENDIF
    ENDIF

    EXITS("Tree_Successor")
    RETURN
999 ERRORSEXITS("Tree_Successor",err,error)
    RETURN
    
  END FUNCTION Tree_Successor

  !
  !================================================================================================================================
  !

  !>Returns values initalise integer value for a tree
  SUBROUTINE Tree_ValueInitialiseGetIntg(tree,valueInitialise,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree
    INTEGER(INTG), INTENT(OUT) :: valueInitialise !<On return, the integer value initialise for the tree.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    
    ENTERS("Tree_ValueInitialiseGetIntg",err,error,*999)

    CALL Tree_AssertNotFinished(tree,err,error,*999)
    CALL Tree_AssertIsIntgData(tree,err,error,*999)
    
    valueInitialise=tree%valueInitialiseIntg

    EXITS("Tree_ValueInitaliseGetIntg")
    RETURN
999 ERRORSEXITS("Tree_ValueInitialiseGetIntg",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_ValueInitialiseGetIntg
  
  !
  !================================================================================================================================
  !

  !>Returns the values initalise double precision real value for a tree
  SUBROUTINE Tree_ValueInitialiseGetDP(tree,valueInitialise,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree
    REAL(DP), INTENT(OUT) :: valueInitialise !<On return, the double precision real value initialise for the tree.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    
    ENTERS("Tree_ValueInitialiseGetDP",err,error,*999)

    CALL Tree_AssertNotFinished(tree,err,error,*999)
    CALL Tree_AssertIsDPData(tree,err,error,*999)
    
    valueInitialise=tree%valueInitialiseDP

    EXITS("Tree_ValueInitaliseSetGP")
    RETURN
999 ERRORSEXITS("Tree_ValueInitialiseSetGP",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_ValueInitialiseGetDP
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the values initalise integer value for a tree
  SUBROUTINE Tree_ValueInitialiseSetIntg(tree,valueInitialise,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree
    INTEGER(INTG), INTENT(IN) :: valueInitialise !<The value initialise for the tree to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    
    ENTERS("Tree_ValueInitialiseSetIntg",err,error,*999)

    CALL Tree_AssertNotFinished(tree,err,error,*999)
    CALL Tree_AssertIsIntgData(tree,err,error,*999)
    
    tree%valueInitialiseIntg=valueInitialise

    EXITS("Tree_ValueInitaliseSetIntg")
    RETURN
999 ERRORSEXITS("Tree_ValueInitialiseSetIntg",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_ValueInitialiseSetIntg
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the values initalise double precision real value for a tree
  SUBROUTINE Tree_ValueInitialiseSetDP(tree,valueInitialise,err,error,*)

    !Argument Variables
    TYPE(TreeType), POINTER :: tree !<A pointer to the tree
    REAL(DP), INTENT(IN) :: valueInitialise !<The value initialise for the tree to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    
    ENTERS("Tree_ValueInitialiseSetDP",err,error,*999)

    CALL Tree_AssertNotFinished(tree,err,error,*999)
    CALL Tree_AssertIsDPData(tree,err,error,*999)
    
    tree%valueInitialiseDP=valueInitialise

    EXITS("Tree_ValueInitaliseSetDP")
    RETURN
999 ERRORSEXITS("Tree_ValueInitialiseSetDP",err,error)
    RETURN 1
    
  END SUBROUTINE Tree_ValueInitialiseSetDP
  
  !  
  !================================================================================================================================
  !


END MODULE Trees
