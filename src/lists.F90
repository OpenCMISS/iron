!> \file
!> \author Chris Bradley
!> \brief Implements lists of base types.
!> \todo Fix up and have this module use the sorting module for sorts
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

!> Implements lists of base types.
MODULE Lists

  USE BaseRoutines
  USE Constants
  USE ISO_VARYING_STRING
  USE Kinds
  USE Strings
  USE Types

#include "macros.h"  
  
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup OpenCMISS_ListsConstants OpenCMISS::Iron::Lists::Constants
  !> \brief Lists constants.
  !>@{
  !> \addtogroup Lists_DataType Lists::DataTypes
  !> \brief Data type parameters for a list.
  !> \see Lists
  !>@{
  INTEGER(INTG), PARAMETER :: LIST_INTG_TYPE=INTEGER_TYPE !<Integer data type for a list \see Lists_DataType,Lists
  INTEGER(INTG), PARAMETER :: LIST_SP_TYPE=SINGLE_REAL_TYPE !<Single precision real data type for a list \see Lists_DataTypes,Lists
  INTEGER(INTG), PARAMETER :: LIST_DP_TYPE=DOUBLE_REAL_TYPE !<Double precision real data type for a list \see Lists_DataTypes,Lists
  !>@}
  
  !> \addtogroup Lists_SortingOrder Lists::SortingOrder
  !> \brief Sorting order parameters for a list.
  !> \see Lists
  !>@{
  INTEGER(INTG), PARAMETER :: LIST_UNSORTED_TYPE=1 !<Unsorted list type \see Lists_SortingOrder,Lists
  INTEGER(INTG), PARAMETER :: LIST_SORT_ASCENDING_TYPE=2 !<Ascending order for sort \see Lists_SortingOrder,Lists
  INTEGER(INTG), PARAMETER :: LIST_SORT_DESCENDING_TYPE=3 !<Descending order for sort \see Lists_SortingOrder,Lists
  !>@}


  !\todo Change lists sorting to use sorting module
  
  !> \addtogroup Lists_SortingMethod Lists::SortingMethod
  !> \brief Sorting method parameters for a list.
  !> \see Lists
  !>@{
  INTEGER(INTG), PARAMETER :: LIST_BUBBLE_SORT_METHOD=1 !<Bubble sort method \see Lists_SortingMethod,Lists
  INTEGER(INTG), PARAMETER :: LIST_SHELL_SORT_METHOD=2 !<Shell sort method \see Lists_SortingMethod,Lists
  INTEGER(INTG), PARAMETER :: LIST_HEAP_SORT_METHOD=3 !<Heap sort method \see Lists_SortingMethod,Lists
  !>@}

  !> \addtogroup Lists_SearchingMethod Lists::SeachingMethod
  !> \brief Searching method parameters for a list.
  !> \see Lists
  !>@{
  INTEGER(INTG), PARAMETER :: LIST_LINEAR_SEARCH_METHOD=1 !<Bubble sort method \see Lists_SeachingMethod,Lists
  INTEGER(INTG), PARAMETER :: LIST_BINARY_SEARCH_METHOD=2 !<Shell sort method \see Lists_SeachingMethod,Lists
  !>@}
  !>@}

  !Module types

  !Module variables
  
  !Interfaces

  INTERFACE LIST_CREATE_FINISH
    MODULE PROCEDURE List_CreateFinish
  END INTERFACE LIST_CREATE_FINISH
  
  INTERFACE LIST_CREATE_START
    MODULE PROCEDURE List_CreateStart
  END INTERFACE LIST_CREATE_START
  
  INTERFACE LIST_DATA_DIMENSION_SET
    MODULE PROCEDURE List_DataDimensionSet
  END INTERFACE LIST_DATA_DIMENSION_SET
  
  INTERFACE LIST_DATA_TYPE_SET
    MODULE PROCEDURE List_DataTypeSet
  END INTERFACE LIST_DATA_TYPE_SET
  
  !>Detaches the list values from a list and returns them as a pointer to a array of base type before destroying the list \see Lists.
  INTERFACE LIST_DETACH_AND_DESTROY
    MODULE PROCEDURE List_DetachAndDestroyIntg1
    MODULE PROCEDURE List_DetachAndDestroyIntg2
    MODULE PROCEDURE List_DetachAndDestroySP1
    MODULE PROCEDURE List_DetachAndDestroySP2
    MODULE PROCEDURE List_DetachAndDestroyDP1
    MODULE PROCEDURE List_DetachAndDestroyDP2
  END INTERFACE LIST_DETACH_AND_DESTROY

  !>Detaches the list values from a list and returns them as a pointer to a array of base type before destroying the list \see Lists.
  INTERFACE List_DetachAndDestroy
    MODULE PROCEDURE List_DetachAndDestroyIntg1
    MODULE PROCEDURE List_DetachAndDestroyIntg2
    MODULE PROCEDURE List_DetachAndDestroySP1
    MODULE PROCEDURE List_DetachAndDestroySP2
    MODULE PROCEDURE List_DetachAndDestroyDP1
    MODULE PROCEDURE List_DetachAndDestroyDP2
  END INTERFACE List_DetachAndDestroy

  !>Adds an item to the end of a list \see Lists.
  INTERFACE LIST_ITEM_ADD
    MODULE PROCEDURE List_ItemAddIntg1
    MODULE PROCEDURE List_ItemAddIntg2
    MODULE PROCEDURE List_ItemAddSP1
    MODULE PROCEDURE List_ItemAddSP2
    MODULE PROCEDURE List_ItemAddDP1
    MODULE PROCEDURE List_ItemAddDP2
  END INTERFACE LIST_ITEM_ADD
  
  !>Adds an item to the end of a list \see Lists.
  INTERFACE List_ItemAdd
    MODULE PROCEDURE List_ItemAddIntg1
    MODULE PROCEDURE List_ItemAddIntg2
    MODULE PROCEDURE List_ItemAddSP1
    MODULE PROCEDURE List_ItemAddSP2
    MODULE PROCEDURE List_ItemAddDP1
    MODULE PROCEDURE List_ItemAddDP2
  END INTERFACE List_ItemAdd
  
  INTERFACE LIST_ITEM_DELETE
    MODULE PROCEDURE List_ItemDelete
  END INTERFACE LIST_ITEM_DELETE
  
  !>Sets an item in the list \see Lists.
  INTERFACE LIST_ITEM_SET
    MODULE PROCEDURE List_ItemSetIntg1
    MODULE PROCEDURE List_ItemSetIntg2
    MODULE PROCEDURE List_ItemSetSP1
    MODULE PROCEDURE List_ItemSetSP2
    MODULE PROCEDURE List_ItemSetDP1
    MODULE PROCEDURE List_ItemSetDP2
  END INTERFACE LIST_ITEM_SET
  
  !>Sets an item in the list \see Lists.
  INTERFACE List_ItemSet
    MODULE PROCEDURE List_ItemSetIntg1
    MODULE PROCEDURE List_ItemSetIntg2
    MODULE PROCEDURE List_ItemSetSP1
    MODULE PROCEDURE List_ItemSetSP2
    MODULE PROCEDURE List_ItemSetDP1
    MODULE PROCEDURE List_ItemSetDP2
  END INTERFACE List_ItemSet
  
  !>Returns an item in a list at a specififed position. \see Lists.
  INTERFACE LIST_ITEM_GET
    MODULE PROCEDURE List_ItemGetIntg1
    MODULE PROCEDURE List_ItemGetIntg2
    MODULE PROCEDURE List_ItemGetSP1
    MODULE PROCEDURE List_ItemGetSP2
    MODULE PROCEDURE List_ItemGetDP1
    MODULE PROCEDURE List_ItemGetDP2
  END INTERFACE LIST_ITEM_GET

  !>Returns an item in a list at a specififed position. \see Lists.
  INTERFACE List_ItemGet
    MODULE PROCEDURE List_ItemGetIntg1
    MODULE PROCEDURE List_ItemGetIntg2
    MODULE PROCEDURE List_ItemGetSP1
    MODULE PROCEDURE List_ItemGetSP2
    MODULE PROCEDURE List_ItemGetDP1
    MODULE PROCEDURE List_ItemGetDP2
  END INTERFACE List_ItemGet

  !>Determines if an item is in a list and returns the position of the item \see Lists.
  INTERFACE LIST_ITEM_IN_LIST
    MODULE PROCEDURE List_ItemInListIntg1
    MODULE PROCEDURE List_ItemInListIntg2
    MODULE PROCEDURE List_ItemInListSP1
    MODULE PROCEDURE List_ItemInListSP2
    MODULE PROCEDURE List_ItemInListDP1
    MODULE PROCEDURE List_ItemInListDP2
  END INTERFACE LIST_ITEM_IN_LIST

  !>Determines if an item is in a list and returns the position of the item \see Lists.
  INTERFACE List_ItemInList
    MODULE PROCEDURE List_ItemInListIntg1
    MODULE PROCEDURE List_ItemInListIntg2
    MODULE PROCEDURE List_ItemInListSP1
    MODULE PROCEDURE List_ItemInListSP2
    MODULE PROCEDURE List_ItemInListDP1
    MODULE PROCEDURE List_ItemInListDP2
  END INTERFACE List_ItemInList

  INTERFACE LIST_INITIAL_SIZE_SET
    MODULE PROCEDURE List_InitialSizeSet
  END INTERFACE LIST_INITIAL_SIZE_SET
  
  INTERFACE LIST_KEY_DIMENSION_SET
    MODULE PROCEDURE List_KeyDimensionSet
  END INTERFACE LIST_KEY_DIMENSION_SET
  
  INTERFACE LIST_MUTABLE_SET
    MODULE PROCEDURE List_MutableSet
  END INTERFACE LIST_MUTABLE_SET
  
  INTERFACE LIST_NUMBER_OF_ITEMS_GET
    MODULE PROCEDURE List_NumberOfItemsGet
  END INTERFACE LIST_NUMBER_OF_ITEMS_GET
  
  INTERFACE LIST_REMOVE_DUPLICATES
    MODULE PROCEDURE List_RemoveDuplicates
  END INTERFACE LIST_REMOVE_DUPLICATES
  
  !>Searches a list for a given value and returns the position in the list if the value exists \see Lists.
  INTERFACE List_Search
    MODULE PROCEDURE List_SearchListIntg1
    MODULE PROCEDURE List_SearchListIntg2
    MODULE PROCEDURE List_SearchListSP1
    MODULE PROCEDURE List_SearchListSP2
    MODULE PROCEDURE List_SearchListDP1
    MODULE PROCEDURE List_SearchListDP2
    MODULE PROCEDURE List_SearchIntg1Array
    MODULE PROCEDURE List_SearchIntg2Array
    MODULE PROCEDURE List_SearchSP1Array
    MODULE PROCEDURE List_SearchSP2Array
    MODULE PROCEDURE List_SearchDP1Array
    MODULE PROCEDURE List_SearchDP2Array
  END INTERFACE List_Search

  !>Searches a list using the binary search method.
  INTERFACE List_SearchBinary
    MODULE PROCEDURE List_SearchBinaryIntg1Array
    MODULE PROCEDURE List_SearchBinaryIntg2Array
    MODULE PROCEDURE List_SearchBinarySP1Array
    MODULE PROCEDURE List_SearchBinarySP2Array
    MODULE PROCEDURE List_SearchBinaryDP1Array
    MODULE PROCEDURE List_SearchBinaryDP2Array
  END INTERFACE List_SearchBinary
  
  !>Searches a list using the linear search method.
  INTERFACE LIST_SEARCH_LINEAR
    MODULE PROCEDURE List_SearchLinearIntg1Array
    MODULE PROCEDURE List_SearchLinearIntg2Array
    MODULE PROCEDURE List_SearchLinearSP1Array
    MODULE PROCEDURE List_SearchLinearSP2Array
    MODULE PROCEDURE List_SearchLinearDP1Array
    MODULE PROCEDURE List_SearchLinearDP2Array
  END INTERFACE LIST_SEARCH_LINEAR

  !>Searches a list using the linear search method.
  INTERFACE List_SearchLinear
    MODULE PROCEDURE List_SearchLinearIntg1Array
    MODULE PROCEDURE List_SearchLinearIntg2Array
    MODULE PROCEDURE List_SearchLinearSP1Array
    MODULE PROCEDURE List_SearchLinearSP2Array
    MODULE PROCEDURE List_SearchLinearDP1Array
    MODULE PROCEDURE List_SearchLinearDP2Array
  END INTERFACE List_SearchLinear

  !>Sorts a list into ascending order.
  INTERFACE List_Sort
    MODULE PROCEDURE List_SortList
    MODULE PROCEDURE List_SortIntg1Array
    MODULE PROCEDURE List_SortIntg2Array
    MODULE PROCEDURE List_SortSP1Array
    MODULE PROCEDURE List_SortSP2Array
    MODULE PROCEDURE List_SortDP1Array
    MODULE PROCEDURE List_SortDP2Array
  END INTERFACE List_Sort

  !>Sorts a list into assending order using the bubble sort method.
  INTERFACE LIST_SORT_BUBBLE
    MODULE PROCEDURE List_SortBubbleIntg1Array
    MODULE PROCEDURE List_SortBubbleIntg2Array
    MODULE PROCEDURE List_SortBubbleSP1Array
    MODULE PROCEDURE List_SortBubbleSP2Array
    MODULE PROCEDURE List_SortBubbleDP1Array
    MODULE PROCEDURE List_SortBubbleDP2Array
  END INTERFACE LIST_SORT_BUBBLE

  !>Sorts a list into assending order using the bubble sort method.
  INTERFACE List_SortBubble
    MODULE PROCEDURE List_SortBubbleIntg1Array
    MODULE PROCEDURE List_SortBubbleIntg2Array
    MODULE PROCEDURE List_SortBubbleSP1Array
    MODULE PROCEDURE List_SortBubbleSP2Array
    MODULE PROCEDURE List_SortBubbleDP1Array
    MODULE PROCEDURE List_SortBubbleDP2Array
  END INTERFACE List_SortBubble

  !>Sorts a list into assending order using the heap sort method.
  INTERFACE LIST_SORT_HEAP
    MODULE PROCEDURE List_SortHeapIntg1Array
    MODULE PROCEDURE List_SortHeapIntg2Array
    MODULE PROCEDURE List_SortHeapSP1Array
    MODULE PROCEDURE List_SortHeapSP2Array
    MODULE PROCEDURE List_SortHeapDP1Array
    MODULE PROCEDURE List_SortHeapDP2Array
  END INTERFACE LIST_SORT_HEAP

  !>Sorts a list into assending order using the heap sort method.
  INTERFACE List_SortHeap
    MODULE PROCEDURE List_SortHeapIntg1Array
    MODULE PROCEDURE List_SortHeapIntg2Array
    MODULE PROCEDURE List_SortHeapSP1Array
    MODULE PROCEDURE List_SortHeapSP2Array
    MODULE PROCEDURE List_SortHeapDP1Array
    MODULE PROCEDURE List_SortHeapDP2Array
  END INTERFACE List_SortHeap

  !>Sorts a list into either assending or descending order using the shell sort method.
  INTERFACE LIST_SORT_SHELL
    MODULE PROCEDURE List_SortShellIntg1Array
    MODULE PROCEDURE List_SortShellIntg2Array
    MODULE PROCEDURE List_SortShellSP1Array
    MODULE PROCEDURE List_SortShellSP2Array
    MODULE PROCEDURE List_SortShellDP1Array
    MODULE PROCEDURE List_SortShellDP2Array
  END INTERFACE LIST_SORT_SHELL

  !>Sorts a list into either assending or descending order using the shell sort method.
  INTERFACE List_SortShell
    MODULE PROCEDURE List_SortShellIntg1Array
    MODULE PROCEDURE List_SortShellIntg2Array
    MODULE PROCEDURE List_SortShellSP1Array
    MODULE PROCEDURE List_SortShellSP2Array
    MODULE PROCEDURE List_SortShellDP1Array
    MODULE PROCEDURE List_SortShellDP2Array
  END INTERFACE List_SortShell

  !>Calculates the intersection of two arrays
  INTERFACE List_Itersection
    MODULE PROCEDURE List_IntersectionIntgArray
  END INTERFACE List_Itersection

  !>Checks whether an array is a subset of another array
  INTERFACE LIST_SUBSET_OF
    MODULE PROCEDURE List_SubsetOfIntgArray
  END INTERFACE LIST_SUBSET_OF

  !>Checks whether an array is a subset of another array
  INTERFACE List_SubsetOf
    MODULE PROCEDURE List_SubsetOfIntgArray
  END INTERFACE List_SubsetOf

  PUBLIC LIST_INTG_TYPE,LIST_SP_TYPE,LIST_DP_TYPE

  PUBLIC LIST_BUBBLE_SORT_METHOD,LIST_SHELL_SORT_METHOD,LIST_HEAP_SORT_METHOD

  PUBLIC LIST_LINEAR_SEARCH_METHOD,LIST_BINARY_SEARCH_METHOD
  
  PUBLIC List_AppendList
  
  PUBLIC List_ClearItems

  PUBLIC LIST_CREATE_FINISH,LIST_CREATE_START

  PUBLIC List_CreateFinish,List_CreateStart
  
  PUBLIC LIST_DATA_DIMENSION_SET

  PUBLIC List_DataDimensionSet
  
  PUBLIC LIST_DATA_TYPE_SET

  PUBLIC List_DataTypeSet

  PUBLIC LIST_DETACH_AND_DESTROY

  PUBLIC List_Destroy,List_DetachAndDestroy

  PUBLIC LIST_ITEM_ADD

  PUBLIC List_ItemAdd
  
  PUBLIC LIST_ITEM_DELETE

  PUBLIC List_ItemDelete

  PUBLIC LIST_ITEM_GET

  PUBLIC List_ItemGet

  PUBLIC LIST_ITEM_IN_LIST

  PUBLIC List_ItemInList

  PUBLIC LIST_ITEM_SET

  PUBLIC List_ItemSet
  
  PUBLIC LIST_INITIAL_SIZE_SET

  PUBLIC List_InitialSizeSet

  PUBLIC LIST_KEY_DIMENSION_SET

  PUBLIC List_KeyDimensionSet

  PUBLIC LIST_MUTABLE_SET

  PUBLIC List_MutableSet

  PUBLIC LIST_NUMBER_OF_ITEMS_GET

  PUBLIC List_NumberOfItemsGet

  PUBLIC LIST_REMOVE_DUPLICATES

  PUBLIC List_RemoveDuplicates

  PUBLIC List_SearchingMethodSet
  
  PUBLIC List_Search

  PUBLIC LIST_SEARCH_LINEAR

  PUBLIC List_SearchBinary,List_SearchLinear
  
  PUBLIC List_SortingMethodSet

  PUBLIC List_Sort
  
  PUBLIC LIST_SORT_BUBBLE,LIST_SORT_HEAP,LIST_SORT_SHELL
  
  PUBLIC List_SortBubble,List_SortHeap,List_SortShell

  PUBLIC List_Itersection

  PUBLIC LIST_SUBSET_OF
  
  PUBLIC List_SubsetOf

CONTAINS

  !
  !================================================================================================================================
  !

  !>Assert that a list has been finished
  SUBROUTINE List_AssertIsFinished(list,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The list to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("List_AssertIsFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(list)) CALL FlagError("List is not associated.",err,error,*999)

    IF(.NOT.list%listFinished) CALL FlagError("List has not been finished.",err,error,*999)
    
    EXITS("List_AssertIsFinished")
    RETURN
999 ERRORSEXITS("List_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE List_AssertIsFinished

  !
  !================================================================================================================================
  !

  !>Assert that a list has not been finished
  SUBROUTINE List_AssertNotFinished(list,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The list to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("List_AssertNotFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(list)) CALL FlagError("List is not associated.",err,error,*999)

    IF(list%listFinished) CALL FlagError("List has already been finished.",err,error,*999)
    
    EXITS("List_AssertNotFinished")
    RETURN
999 ERRORSEXITS("List_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE List_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Assert that a list is mutable
  SUBROUTINE List_AssertIsMutable(list,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The list to assert the mutable status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("List_AssertIsMutable",err,error,*999)

    IF(.NOT.ASSOCIATED(list)) CALL FlagError("List is not associated.",err,error,*999)

    IF(.NOT.list%mutable) CALL FlagError("List is not mutable.",err,error,*999)
    
    EXITS("List_AssertIsMutable")
    RETURN
999 ERRORSEXITS("List_AssertIsMutable",err,error)
    RETURN 1
    
  END SUBROUTINE List_AssertIsMutable

  !
  !================================================================================================================================
  !

  !>Assert that a list is not mutable
  SUBROUTINE List_AssertNotMutable(list,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The list to assert the mutable status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("List_AssertNotMutable",err,error,*999)

    IF(.NOT.ASSOCIATED(list)) CALL FlagError("List is not associated.",err,error,*999)

    IF(list%mutable) CALL FlagError("List is mutable.",err,error,*999)
    
    EXITS("List_AssertNotMutable")
    RETURN
999 ERRORSEXITS("List_AssertNotMutable",err,error)
    RETURN 1
    
  END SUBROUTINE List_AssertNotMutable

  !
  !================================================================================================================================
  !

  !>Assert that a list has an integer data type
  SUBROUTINE List_AssertIsIntgData(list,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The list to assert the integer data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("List_AssertIsIntgData",err,error,*999)

    IF(.NOT.ASSOCIATED(list)) CALL FlagError("List is not associated.",err,error,*999)

    IF(list%dataType/=LIST_INTG_TYPE) THEN
      localError="The list data type of "//TRIM(NumberToVString(list%dataType,"*",err,error))// &
        & " does not correspond to the required integer data type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("List_AssertIsIntgData")
    RETURN
999 ERRORSEXITS("List_AssertIsIntgData",err,error)
    RETURN 1
    
  END SUBROUTINE List_AssertIsIntgData

  !
  !================================================================================================================================
  !

  !>Assert that a list has a single precision real data type
  SUBROUTINE List_AssertIsSPData(list,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The list to assert the single precision real data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("List_AssertIsSPData",err,error,*999)

    IF(.NOT.ASSOCIATED(list)) CALL FlagError("List is not associated.",err,error,*999)

    IF(list%dataType/=LIST_SP_TYPE) THEN
      localError="The list data type of "//TRIM(NumberToVString(list%dataType,"*",err,error))// &
        & " does not correspond to the required single precision real data type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("List_AssertIsSPData")
    RETURN
999 ERRORSEXITS("List_AssertIsSPData",err,error)
    RETURN 1
    
  END SUBROUTINE List_AssertIsSPData

  !
  !================================================================================================================================
  !

  !>Assert that a list has a double precision real data type
  SUBROUTINE List_AssertIsDPData(list,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The list to assert the double precision real data type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("List_AssertIsDPData",err,error,*999)

    IF(.NOT.ASSOCIATED(list)) CALL FlagError("List is not associated.",err,error,*999)

    IF(list%dataType/=LIST_DP_TYPE) THEN
      localError="The list data type of "//TRIM(NumberToVString(list%dataType,"*",err,error))// &
        & " does not correspond to the required double precision real data type."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("List_AssertIsDPData")
    RETURN
999 ERRORSEXITS("List_AssertIsDPData",err,error)
    RETURN 1
    
  END SUBROUTINE List_AssertIsDPData

 !
  !================================================================================================================================
  !

  !>Finishes the creation of a list created with List_CreateStart \see Lists::List_CreateStart.
  SUBROUTINE List_CreateFinish(list,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(INOUT) :: list !<A pointer to the list to finish
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("List_CreateFinish",err,error,*998)

    CALL List_AssertNotFinished(list,err,error,*999)

    !Allocate the list
    IF(list%dataDimension==1) THEN
      SELECT CASE(list%dataType)
      CASE(LIST_INTG_TYPE)
        ALLOCATE(list%listIntg(list%initialSize),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate list integer data.",err,error,*999)
      CASE(LIST_SP_TYPE)
        ALLOCATE(list%listSP(list%initialSize),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate list single precision data.",err,error,*999)
      CASE(LIST_DP_TYPE)
        ALLOCATE(list%listDP(list%initialSize),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate list double precision data.",err,error,*999)
      CASE DEFAULT
        localError="The list data type of "//TRIM(NumberToVString(list%dataType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      SELECT CASE(list%dataType)
      CASE(LIST_INTG_TYPE)
        ALLOCATE(list%listIntg2(list%dataDimension,list%initialSize),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate list integer data.",err,error,*999)
      CASE(LIST_SP_TYPE)
        ALLOCATE(list%listSP2(list%dataDimension,list%initialSize),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate list single precision data.",err,error,*999)
      CASE(LIST_DP_TYPE)
        ALLOCATE(list%listDP2(list%dataDimension,list%initialSize),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate list double precision data.",err,error,*999)
      CASE DEFAULT
        localError="The list data type of "//TRIM(NumberToVString(list%dataType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    list%size=list%initialSize
    list%listFinished=.TRUE.

    EXITS("List_CreateFinish")
    RETURN
999 CALL List_Finalise(list,dummyErr,dummyError,*998)
998 ERRORSEXITS("List_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE List_CreateFinish

  !
  !================================================================================================================================
  !

  !>Starts the creation of a list and returns a pointer to the created list \see Lists::List_CreateFinish.
  SUBROUTINE List_CreateStart(list,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(INOUT) :: list !<On exit, pointer to the list to create. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables

    ENTERS("List_CreateStart",err,error,*999)

    CALL List_Initialise(list,err,error,*999)
    
    EXITS("List_CreateStart")
    RETURN
999 ERRORSEXITS("List_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE List_CreateStart

  !
  !================================================================================================================================
  !

  !>Sets/changes the data dimension for a list.
  SUBROUTINE List_DataDimensionSet(list,dataDimension,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<A pointer to the list 
    INTEGER(INTG), INTENT(IN) :: dataDimension !<The data dimension of the list to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("List_DataDimensionSet",err,error,*999)

    CALL List_AssertNotFinished(list,err,error,*999)
    IF(dataDimension<=0) THEN
      localError="The specified data dimension of "//TRIM(NumberToVString(dataDimension,"*",err,error))// &
        & " is invalid. The dimension must be > 0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    list%dataDimension=dataDimension

    EXITS("List_DataDimensionSet")
    RETURN
999 ERRORSEXITS("List_DataDimensionSet",err,error)
    RETURN 1
    
  END SUBROUTINE List_DataDimensionSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the data dimension for a list.
  SUBROUTINE List_MutableSet(list,mutable,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<A pointer to the list 
    LOGICAL, INTENT(IN) :: mutable !<The mutability of the list to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("List_MutableSet",err,error,*999)

    CALL List_AssertNotFinished(list,err,error,*999)
    
    list%mutable = mutable
 
    EXITS("List_MutableSet")
    RETURN
999 ERRORSEXITS("List_MutableSet",err,error)
    RETURN 1
    
  END SUBROUTINE List_MutableSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the data type for a list.
  SUBROUTINE List_DataTypeSet(list,dataType,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<A pointer to the list 
    INTEGER(INTG), INTENT(IN) :: dataType !<The data type of the list to set \see Lists_DataTypes,Lists
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("List_DataTypeSet",err,error,*999)

    CALL List_AssertNotFinished(list,err,error,*999)
    
    SELECT CASE(dataType)
    CASE(LIST_INTG_TYPE)
      list%dataType=LIST_INTG_TYPE
    CASE(LIST_SP_TYPE)
      list%dataType=LIST_SP_TYPE
    CASE(LIST_DP_TYPE)
      list%dataType=LIST_DP_TYPE
    CASE DEFAULT
      localError="The data type of "//TRIM(NumberToVString(dataType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("List_DataTypeSet")
    RETURN
999 ERRORSEXITS("List_DataTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE List_DataTypeSet

  !
  !================================================================================================================================
  !

  !>Destroys a list.
  SUBROUTINE List_Destroy(list,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(INOUT) :: list !<A pointer to the list to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("List_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(list)) CALL FlagError("List is not associated.",err,error,*999)
    
    CALL List_Finalise(list,err,error,*999)

    EXITS("List_Destroy")
    RETURN
999 ERRORSEXITS("List_Destroy",err,error)
    RETURN 1
  END SUBROUTINE List_Destroy

  !
  !================================================================================================================================
  !

  !>Finalises a list and deallocates all memory.
  SUBROUTINE List_Finalise(list,err,error,*)    

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(INOUT) :: list !<A pointer to the list to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("List_Finalise",err,error,*999)

    IF(ASSOCIATED(list)) THEN
      IF(ALLOCATED(list%listIntg)) DEALLOCATE(list%listIntg)
      IF(ALLOCATED(list%listIntg2)) DEALLOCATE(list%listIntg2)
      IF(ALLOCATED(list%listSP)) DEALLOCATE(list%listSP)
      IF(ALLOCATED(list%listSP2)) DEALLOCATE(list%listSP2)
      IF(ALLOCATED(list%listDP)) DEALLOCATE(list%listDP)
      IF(ALLOCATED(list%listDP2)) DEALLOCATE(list%listDP2)
      DEALLOCATE(list)
    ENDIF

    EXITS("List_Finalise")
    RETURN
999 ERRORSEXITS("List_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE List_Finalise

  !
  !================================================================================================================================
  !

  !>Appends a list to the end of this list
  SUBROUTINE List_AppendList(list,appendedList,err,error,*)
   !Argument Variables
    TYPE(ListType), POINTER, INTENT(INOUT) :: list !<A pointer to the list
    TYPE(ListType), POINTER, INTENT(IN) :: appendedList !<The list to append
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: newSize
    INTEGER(INTG), ALLOCATABLE :: newListIntg(:)
    REAL(SP), ALLOCATABLE :: newListSP(:)
    REAL(DP), ALLOCATABLE :: newListDP(:)
    TYPE(VARYING_STRING) :: localError

    ENTERS("List_AppendList",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsFinished(appendedList,err,error,*999)

    IF(list%dataType/=appendedList%dataType) THEN
      localError="The list data type of "//TRIM(NumberToVString(list%dataType,"*",err,error))// &
        & " does not match the data type of the list to append."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(list%dataDimension/=appendedList%dataDimension) THEN
      localError="Invalid data dimension. The list to append has data dimension of "// &
        & TRIM(NumberToVString(appendedList%dataDimension,"*",err,error))//" and the list data dimension is "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    SELECT CASE(list%dataDimension)
    CASE(1)
      SELECT CASE(list%dataType)
      CASE(LIST_INTG_TYPE)
        IF(list%size<list%numberInList+appendedList%numberInList) THEN
          !Reallocate
          newSize=MAX(2*list%numberInList,list%numberInList+appendedList%numberInList*2)
          ALLOCATE(newListIntg(newSize),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate new list.",err,ERROR,*999)
          newListIntg(1:list%numberInList)=list%listIntg(1:list%numberInList)
          CALL MOVE_ALLOC(newListIntg,list%listIntg)
          list%size=newSize
        END IF
        list%listIntg(list%numberInList+1:list%numberInList+appendedList%numberInList)= &
          & appendedList%listIntg(1:appendedList%numberInList)
        list%numberInList=list%numberInList+appendedList%numberInList
      CASE(LIST_SP_TYPE)
        IF(list%size<list%numberInList+appendedList%numberInList) THEN
          !Reallocate
          newSize=MAX(2*list%numberInList,list%numberInList+appendedList%numberInList*2)
          ALLOCATE(newListSP(newSize),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate new list.",err,ERROR,*999)
          newListSP(1:list%numberInList)=list%listSP(1:list%numberInList)
          CALL MOVE_ALLOC(newListSP,list%listSP)
          list%size=newSize
        END IF
        list%listSP(list%numberInList+1:list%numberInList+appendedList%numberInList)= &
          & appendedList%listSP(1:appendedList%numberInList)
        list%numberInList=list%numberInList+appendedList%numberInList
      CASE(LIST_DP_TYPE)
        IF(list%size<list%numberInList+appendedList%numberInList) THEN
          !Reallocate
          newSize=MAX(2*list%numberInList,list%numberInList+appendedList%numberInList*2)
          ALLOCATE(newListDP(newSize),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate new list.",err,ERROR,*999)
          newListDP(1:list%numberInList)=list%listDP(1:list%numberInList)
          CALL MOVE_ALLOC(newListDP,list%listDP)
          list%size=newSize
        END IF
        list%listDP(list%numberInList+1:list%numberInList+appendedList%numberInList)= &
          & appendedList%listDP(1:appendedList%numberInList)
        list%numberInList=list%numberInList+appendedList%numberInList
      CASE DEFAULT
        CALL FlagError("The list data type of "//TRIM(NumberToVString(list%dataType,"*",err,error))// &
          & " is invalid.",err,error,*999)
      END SELECT
    CASE DEFAULT
      CALL FlagError("Dimensions > 1 not implemented for appended to a list.",err,error,*999)
    END SELECT

    EXITS("List_AppendList")
    RETURN
999 IF(ALLOCATED(newListIntg)) DEALLOCATE(newListIntg)
    IF(ALLOCATED(newListSP)) DEALLOCATE(newListSP)
    IF(ALLOCATED(newListDP)) DEALLOCATE(newListDP)
    ERRORSEXITS("List_AppendList",err,error)
    RETURN 1
    
  END SUBROUTINE List_AppendList

  !
  !================================================================================================================================
  !

  !>Clears all the items from a list
  SUBROUTINE List_ClearItems(list,err,error,*)
   !Argument Variables
    TYPE(ListType), POINTER, INTENT(INOUT) :: list !<A pointer to the list
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("List_ClearItems",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsMutable(list,err,error,*999)
    
    list%numberInList=0

    EXITS("List_ClearItems")
    RETURN
999 ERRORSEXITS("List_ClearItems",err,error)
    RETURN 1
    
  END SUBROUTINE List_ClearItems
  !
  !================================================================================================================================
  !

  !>Initialises a list and all its components
  SUBROUTINE List_Initialise(list,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(INOUT) :: list !<A pointer to the list to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError    

    ENTERS("List_Initialise",err,error,*998)

    IF(ASSOCIATED(list)) CALL FlagError("List is already associated.",err,error,*998)
    
    ALLOCATE(list,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate list.",err,error,*999)
    list%listFinished=.FALSE.
    list%mutable=.FALSE.
    list%numberInList=0
    list%dataDimension=1
    list%initialSize=10
    list%size=0
    list%dataType=LIST_INTG_TYPE
    list%keyDimension=1
    list%sortOrder=LIST_SORT_ASCENDING_TYPE
    list%sortMethod=LIST_HEAP_SORT_METHOD
    list%searchMethod=LIST_LINEAR_SEARCH_METHOD

    EXITS("List_Initialise")
    RETURN
999 CALL List_Finalise(list,dummyErr,dummyError,*998)
998 ERRORSEXITS("List_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE List_Initialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the initial size for a list
  SUBROUTINE List_InitialSizeSet(list,initialSize,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: initialSize !<The initial size of the list to set. Must be greater than zero.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("List_InitialSizeSet",err,error,*999)

    CALL List_AssertNotFinished(list,err,error,*999)
    IF(initialSize<=0) THEN
      localError="The initial size of "//TRIM(NumberToVString(initialSize,"*",err,error))// &
        & " is invalid. The size must be > 0."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    list%initialSize=initialSize

    EXITS("List_InitialSizeSet")
    RETURN
999 ERRORSEXITS("List_InitialSizeSet",err,error)
    RETURN 1
    
  END SUBROUTINE List_InitialSizeSet

  !
  !================================================================================================================================
  !

  !>Adds an item to the end of an integer list of data dimension 1. 
  SUBROUTINE List_ItemAddIntg1(list,item,err,error,*)
   !Argument Variables
    TYPE(ListType), POINTER, INTENT(INOUT) :: list !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: item !<The item to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: newSize
    INTEGER(INTG), ALLOCATABLE :: newList(:)
    TYPE(VARYING_STRING) :: localError

    ENTERS("List_ItemAddIntg1",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsIntgData(list,err,error,*999)
    IF(list%dataDimension/=1) THEN
      localError="Invalid data dimension. The supplied data dimension is 1 and the list data dimension is "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    IF(list%numberInList==list%size) THEN
      !Reallocate
      newSize=MAX(2*list%numberInList,1)
      ALLOCATE(newList(newSize),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new list.",err,error,*999)
      newList(1:list%numberInList)=list%listIntg(1:list%numberInList)
      CALL MOVE_ALLOC(newList,list%listIntg)
      list%size=newSize
    ENDIF
    list%listIntg(list%numberInList+1)=item
    list%numberInList=list%numberInList+1
    
    EXITS("List_ItemAddIntg1")
    RETURN
999 IF(ALLOCATED(newList)) DEALLOCATE(newList)
    ERRORSEXITS("List_ItemAddIntg1",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemAddIntg1
  
  !
  !================================================================================================================================
  !

  !>Adds an item to the end of an integer list of data dimension > 1. 
  SUBROUTINE List_ItemAddIntg2(list,item,err,error,*)
   !Argument Variables
    TYPE(ListType), POINTER, INTENT(INOUT) :: list !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: item(:) !<The item to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: newSize
    INTEGER(INTG), ALLOCATABLE :: newList(:,:)
    TYPE(VARYING_STRING) :: localError

    ENTERS("List_ItemAddIntg2",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsIntgData(list,err,error,*999)
    IF(list%dataDimension/=SIZE(item,1)) THEN
      localError="Invalid data dimension. The supplied data dimension is "// &
        & TRIM(NumberToVString(SIZE(item,1),"*",err,error))//" and the list data dimension is "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    IF(list%numberInList==list%size) THEN
      !Reallocate
      newSize=MAX(2*list%numberInList,1)
      ALLOCATE(newList(list%dataDimension,newSize),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new list.",err,error,*999)
      newList(:,1:list%numberInList)=list%listIntg2(:,1:list%numberInList)
      CALL MOVE_ALLOC(newList,list%listIntg2)
      list%size=newSize
    ENDIF
    list%listIntg2(:,list%numberInList+1)=item
    list%numberInList=list%numberInList+1
   
    EXITS("List_ItemAddIntg2")
    RETURN
999 IF(ALLOCATED(newList)) DEALLOCATE(newList)
    ERRORSEXITS("List_ItemAddIntg2",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemAddIntg2
  
  !
  !================================================================================================================================
  !

  !>Adds an item to the end of a single precision real list of data dimension 1. 
  SUBROUTINE List_ItemAddSP1(list,item,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(INOUT) :: list !<A pointer to the list
    REAL(SP), INTENT(IN) :: item !<The item to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: newSize
    REAL(SP), ALLOCATABLE :: newList(:)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_ItemAddSP1",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsSPData(list,err,error,*999)
    IF(list%dataDimension/=1) THEN
      localError="Invalid data dimension. The supplied data dimension is 1 and the list data dimension is "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    IF(list%numberInList==list%size) THEN
      !Reallocate
      newSize=MAX(2*list%numberInList,1)
      ALLOCATE(newList(newSize),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new list.",err,error,*999)
      newList(1:list%numberInList)=list%listSP(1:list%numberInList)
      CALL MOVE_ALLOC(newList,list%listSP)
      list%size=newSize
    ENDIF
    list%listSP(list%numberInList+1)=item
    list%numberInList=list%numberInList+1
    
    EXITS("List_ItemAddSP1")
    RETURN
999 IF(ALLOCATED(newList)) DEALLOCATE(newList)
    ERRORSEXITS("List_ItemAddSP1",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemAddSP1
  
  !
  !================================================================================================================================
  !

  !>Adds an item to the end of a single precision real list of data dimension > 1. 
  SUBROUTINE List_ItemAddSP2(list,ITEM,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(INOUT) :: list !<A pointer to the list
    REAL(SP), INTENT(IN) :: ITEM(:) !<The item to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: newSize
    REAL(SP), ALLOCATABLE :: newList(:,:)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_ItemAddSP2",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsSPData(list,err,error,*999)
    IF(list%dataDimension/=SIZE(item,1)) THEN
      localError="Invalid data dimension. The supplied data dimension is "// &
        & TRIM(NumberToVString(SIZE(item,1),"*",err,error))//" and the list data dimension is "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(list%numberInList==list%size) THEN
      !Reallocate
      newSize=MAX(2*list%numberInList,1)
      ALLOCATE(newList(list%dataDimension,newSize),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new list.",err,error,*999)
      newList(:,1:list%numberInList)=list%listSP2(:,1:list%numberInList)
      CALL MOVE_ALLOC(newList,list%listSP2)
      list%size=newSize
    ENDIF
    list%listSP2(:,list%numberInList+1)=item
    list%numberInList=list%numberInList+1
    
    EXITS("List_ItemAddSP2")
    RETURN
999 IF(ALLOCATED(newList)) DEALLOCATE(newList)
    ERRORSEXITS("List_ItemAddSP2",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemAddSP2
  
  !
  !================================================================================================================================
  !

  !>Adds an item to the end of a double precision real list of data dimension 1.
  SUBROUTINE List_ItemAddDP1(list,item,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(INOUT) :: list !<A pointer to the list
    REAL(DP), INTENT(IN) :: item !<The item to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: newSize
    REAL(DP), ALLOCATABLE :: newList(:)
    TYPE(VARYING_STRING) :: localError

    ENTERS("List_ItemAddDP1",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsDPData(list,err,error,*999)
    IF(list%dataDimension/=1) THEN
      localError="Invalid data dimension. The supplied data dimension is 1 and the list data dimension is "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    IF(list%numberInList==list%size) THEN
      !Reallocate
      newSize=MAX(2*list%numberInList,1)
      ALLOCATE(newList(newSize),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new list.",err,error,*999)
      newList(1:list%numberInList)=list%listDP(1:list%numberInList)
      CALL MOVE_ALLOC(newList,list%listDP)
      list%size=newSize
    ENDIF
    list%listDP(list%numberInList+1)=item
    list%numberInList=list%numberInList+1
    
    EXITS("List_ItemAddDP1")
    RETURN
999 IF(ALLOCATED(newList)) DEALLOCATE(newList)
    ERRORSEXITS("List_ItemAddDP1",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemAddDP1
  
  !
  !================================================================================================================================
  !

  !>Adds an item to the end of a double precision real list of data dimension > 1.
  SUBROUTINE List_ItemAddDP2(list,item,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(INOUT) :: list !<A pointer to the list
    REAL(DP), INTENT(IN) :: item(:) !<The item to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: newSize
    REAL(DP), ALLOCATABLE :: newList(:,:)
    TYPE(VARYING_STRING) :: localError

    ENTERS("List_ItemAddDP2",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsDPData(list,err,error,*999)
    IF(list%dataDimension/=SIZE(item,1)) THEN
      localError="Invalid data dimension. The supplied data dimension is "// &
        & TRIM(NumberToVString(SIZE(item,1),"*",err,error))//" and the list data dimension is "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    IF(list%numberInList==list%size) THEN
      !Reallocate
      newSize=MAX(2*list%numberInList,1)
      ALLOCATE(newList(list%dataDimension,newSize),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new list.",err,error,*999)
      newList(:,1:list%numberInList)=list%listDP2(:,1:list%numberInList)
      CALL MOVE_ALLOC(newList,list%listDP2)
      list%size=newSize
    ENDIF
    list%listDP2(:,list%numberInList+1)=item
    list%numberInList=list%numberInList+1
    
    EXITS("List_ItemAddDP2")
    RETURN
999 IF(ALLOCATED(newList)) DEALLOCATE(newList)
    ERRORSEXITS("List_ItemAddDP2",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemAddDP2
  
  !
  !================================================================================================================================
  !

  !>Sets an item in an integer list of data dimension 1. 
  SUBROUTINE List_ItemSetIntg1(list,listItem,item,err,error,*)
   !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: listItem !<The index of the item to set
    INTEGER(INTG), INTENT(IN) :: item !<The item to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("List_ItemSetIntg1",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsIntgData(list,err,error,*999)
    CALL List_AssertIsMutable(list,err,error,*999)
    IF(list%dataDimension/=1) THEN
      localError="Invalid data dimension. The supplied data dimension is 1 and the list data dimension is "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(listItem<=0.OR.listItem>list%numberInList) THEN
      localError="Invalid list index. The supplied index is "// &
        & TRIM(NumberToVString(listItem,"*",err,error))//" and that list entry count is"// &
        & TRIM(NumberToVString(list%numberInList,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    list%listIntg(listItem)=item
    
    EXITS("List_ItemSetIntg1")
    RETURN
999 ERRORSEXITS("List_ItemSetIntg1",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemSetIntg1
  
  !
  !================================================================================================================================
  !

  !>Set an item in an integer list of data dimension > 1. 
  SUBROUTINE List_ItemSetIntg2(list,listItem,item,err,error,*)
   !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: listItem !<The index of the item to set
    INTEGER(INTG), INTENT(IN) :: item(:) !<The item to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("List_ItemSetIntg2",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsIntgData(list,err,error,*999)
    CALL List_AssertIsMutable(list,err,error,*999)
    IF(list%dataDimension/=SIZE(item,1)) THEN
      localError="Invalid data dimension. The supplied data dimension is "// &
        & TRIM(NumberToVString(SIZE(item,1),"*",err,error))//" and the list data dimension is "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(listItem<=0.OR.listItem>list%numberInList) THEN
      localError="Invalid list index. The supplied index is "// &
        & TRIM(NumberToVString(listItem,"*",err,error))//" and that list entry count is"// &
        & TRIM(NumberToVString(list%numberInList,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    list%listIntg2(:,listItem)=item
    
    EXITS("List_ItemSetIntg2")
    RETURN
999 ERRORSEXITS("List_ItemSetIntg2",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemSetIntg2
  
  !
  !================================================================================================================================
  !

  !>Sets an item in a single precision real list of data dimension 1. 
  SUBROUTINE List_ItemSetSP1(list,listItem,item,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: listItem !<The index of the item to set
    REAL(SP), INTENT(IN) :: item !<The item to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_ItemSetSP1",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsSPData(list,err,error,*999)
    CALL List_AssertIsMutable(list,err,error,*999)
    IF(list%dataDimension/=1) THEN
      localError="Invalid data dimension. The supplied data dimension is 1 and the list data dimension is "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(listItem<=0.OR.listItem>list%numberInList) THEN
      localError="Invalid list index. The supplied index is "// &
        & TRIM(NumberToVString(listItem,"*",err,error))//" and that list entry count is"// &
        & TRIM(NumberToVString(list%numberInList,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    list%listSP(listItem)=item
    
    EXITS("List_ItemSetSP1")
    RETURN
999 ERRORSEXITS("List_ItemSetSP1",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemSetSP1
  
  !
  !================================================================================================================================
  !

  !>Sets an item in a single precision real list of data dimension > 1. 
  SUBROUTINE List_ItemSetSP2(list,listItem,item,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: listItem !<The index of the item to set
    REAL(SP), INTENT(IN) :: item(:) !<The item to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_ItemSetSP2",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsSPData(list,err,error,*999)
    CALL List_AssertIsMutable(list,err,error,*999)
    IF(list%dataDimension/=SIZE(item,1)) THEN
      localError="Invalid data dimension. The supplied data dimension is "// &
        & TRIM(NumberToVString(SIZE(item,1),"*",err,error))//" and the list data dimension is "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(listItem<=0.OR.listItem>list%numberInList) THEN
      localError="Invalid list index. The supplied index is "// &
        & TRIM(NumberToVString(listItem,"*",err,error))//" and that list entry count is"// &
        & TRIM(NumberToVString(list%numberInList,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    list%listSP2(:,listItem)=item
    
    EXITS("List_ItemSetSP2")
    RETURN
999 ERRORSEXITS("List_ItemSetSP2",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemSetSP2
  
  !
  !================================================================================================================================
  !

  !>Sets an item in a double precision real list of data dimension 1.
  SUBROUTINE List_ItemSetDP1(list,listItem,item,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: listItem !<The index of the item to set. 
    REAL(DP), INTENT(IN) :: item !<The item to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("List_ItemSetDP1",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsDPData(list,err,error,*999)
    CALL List_AssertIsMutable(list,err,error,*999)
    IF(list%dataDimension/=1) THEN
      localError="Invalid data dimension. The supplied data dimension is 1 and the list data dimension is "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(listItem<=0.OR.listItem>list%numberInList) THEN
      localError="Invalid list index. The supplied index is "// &
        & TRIM(NumberToVString(listItem,"*",err,error))//" and that list entry count is"// &
        & TRIM(NumberToVString(list%numberInList,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    list%listDP(listItem)=item
    
    EXITS("List_ItemSetDP1")
    RETURN
999 ERRORSEXITS("List_ItemSetDP1",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemSetDP1
  
  !
  !================================================================================================================================
  !

  !>Sets an item in a double precision real list of data dimension > 1.
  SUBROUTINE List_ItemSetDP2(list,listItem,item,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: listItem !<The index of the item to set.
    REAL(DP), INTENT(IN) :: item(:) !<The item to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("List_ItemSetDP2",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsDPData(list,err,error,*999)
    CALL List_AssertIsMutable(list,err,error,*999)
    IF(list%dataDimension/=SIZE(item,1)) THEN
      localError="Invalid data dimension. The supplied data dimension is "// &
        & TRIM(NumberToVString(SIZE(item,1),"*",err,error))//" and the list data dimension is "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(listItem<=0.OR.listItem>list%numberInList) THEN
      localError="The specified list item position of "//TRIM(NumberToVString(listItem,"*",err,error))// &
        & " is invalid. The list item position must be > 0 and <= "// &
        & TRIM(NumberToVString(list%numberInList,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    list%listDP2(:,listItem)=item
    
    EXITS("List_ItemSetDP2")
    RETURN
999 ERRORSEXITS("List_ItemSetDP2",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemSetDP2
  
  !
  !================================================================================================================================
  !

  !>Returns the item in a list at position listItem in the given integer list. 
  SUBROUTINE List_ItemGetIntg1(list,listItem,item,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The pointer to the list
    INTEGER(INTG), INTENT(IN) :: listItem !<The position of the item to get
    INTEGER(INTG), INTENT(OUT) :: item !<On exit, the item at the specified position
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_ItemGetIntg1",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsIntgData(list,err,error,*999)
    IF(list%dataDimension/=1) THEN
      localError="Invalid data dimension. The supplied data dimension is 1 and the list data dimension is "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(listItem<=0.OR.listItem>list%numberInList) THEN
      localError="Invalid list index. The supplied index is "// &
        & TRIM(NumberToVString(listItem,"*",err,error))//" and that list entry count is"// &
        & TRIM(NumberToVString(list%numberInList,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    item=list%listIntg(listItem)

    EXITS("List_ItemGetIntg1")
    RETURN
999 ERRORSEXITS("List_ItemGetIntg1",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemGetIntg1
  
  !
  !================================================================================================================================
  !

  !>Returns the item in a list at position listItem in the given integer list. 
  SUBROUTINE List_ItemGetIntg2(list,listItem,item,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The pointer to the list
    INTEGER(INTG), INTENT(IN) :: listItem !<The position of the item to get
    INTEGER(INTG), INTENT(OUT) :: item(:) !<On exit, the item at the specified position
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_ItemGetIntg2",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsIntgData(list,err,error,*999)
    IF(list%dataDimension/=SIZE(item,1)) THEN
      localError="Invalid item dimension. The supplied item has dimension "// &
        & TRIM(NumberToVString(SIZE(item,1),"*",err,error))//" and the list is of dimesion "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(listItem<=0.OR.listItem>list%numberInList) THEN
      localError="The specified list item position of "//TRIM(NumberToVString(listItem,"*",err,error))// &
        & " is invalid. The list item position must be > 0 and <= "// &
        & TRIM(NumberToVString(list%numberInList,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    item=list%listIntg2(:,listItem)

    EXITS("List_ItemGetIntg2")
    RETURN
999 ERRORSEXITS("List_ItemGetIntg2",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemGetIntg2
  
  !
  !================================================================================================================================
  !

  !>Returns the item in a list at position listItem in the given single precision list. 
  SUBROUTINE List_ItemGetSP1(list,listItem,item,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The pointer to the list
    INTEGER(INTG), INTENT(IN) :: listItem !<The position of the item to get
    REAL(SP), INTENT(OUT) :: item !<On exit, the item at the specified position.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_ItemGetSP1",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsSPData(list,err,error,*999)
    IF(list%dataDimension/=1) THEN
      localError="Invalid data dimension. The supplied data dimension is 1 and the list data dimension is "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(listItem<=0.OR.listItem>list%numberInList) THEN
      localError="Invalid list index. The supplied index is "// &
        & TRIM(NumberToVString(listItem,"*",err,error))//" and that list entry count is"// &
        & TRIM(NumberToVString(list%numberInList,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    item=list%listSP(listItem)

    EXITS("List_ItemGetSP1")
    RETURN
999 ERRORSEXITS("List_ItemGetSP1",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemGetSP1
  
  !
  !================================================================================================================================
  !

  !>Returns the item in a list at position listItem in the given single precision list. 
  SUBROUTINE List_ItemGetSP2(list,listItem,item,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The pointer to the list
    INTEGER(INTG), INTENT(IN) :: listItem !<The position of the item to get
    REAL(SP), INTENT(OUT) :: item(:) !<On exit, the item at the specified position
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_ItemGetSP2",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsSPData(list,err,error,*999)
    IF(list%dataDimension/=SIZE(item,1)) THEN
      localError="Invalid item dimension. The supplied item has dimension "// &
        & TRIM(NumberToVString(SIZE(item,1),"*",err,error))//" and the list is of dimesion "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(listItem<=0.OR.listItem>list%numberInList) THEN
      localError="The specified list item position of "//TRIM(NumberToVString(listItem,"*",err,error))// &
        & " is invalid. The list item position must be > 0 and <= "// &
        & TRIM(NumberToVString(list%numberInList,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    item=list%listSP2(:,listItem)

    EXITS("List_ItemGetSP2")
    RETURN
999 ERRORSEXITS("List_ItemGetSP2",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemGetSP2
  
  !
  !================================================================================================================================
  !

  !>Returns the item in a list at position listItem in the given double precision list. 
  SUBROUTINE List_ItemGetDP1(list,listItem,item,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The pointer to the list
    INTEGER(INTG), INTENT(IN) :: listItem !<The position of the item to get
    REAL(DP), INTENT(OUT) :: item !<On exit, the item at the specified position.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_ItemGetDP1",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsDPData(list,err,error,*999)
    IF(list%dataDimension/=1) THEN
      localError="Invalid data dimension. The supplied data dimension is 1 and the list data dimension is "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(listItem<=0.OR.listItem>list%numberInList) THEN
      localError="Invalid list index. The supplied index is "// &
        & TRIM(NumberToVString(listItem,"*",err,error))//" and that list entry count is"// &
        & TRIM(NumberToVString(list%numberInList,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    item=list%listDP(listItem)
 
    EXITS("List_ItemGetDP1")
    RETURN
999 ERRORSEXITS("List_ItemGetDP1",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemGetDP1
  
  !
  !================================================================================================================================
  !

  !>Returns the item in a list at position listItem in the given double precision list. 
  SUBROUTINE List_ItemGetDP2(list,listItem,item,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The pointer to the list
    INTEGER(INTG), INTENT(IN) :: listItem !<The position of the item to get
    REAL(DP), INTENT(OUT) :: item(:) !<On exit, the item at the specified position
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_ItemGetDP2",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsDPData(list,err,error,*999)
    IF(list%dataDimension/=SIZE(item,1)) THEN
      localError="Invalid item dimension. The supplied item has dimension "// &
        & TRIM(NumberToVString(SIZE(item,1),"*",err,error))//" and the list is of dimesion "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(listItem<=0.OR.listItem>list%numberInList) THEN
      localError="The specified list item position of "//TRIM(NumberToVString(listItem,"*",err,error))// &
        & " is invalid. The list item position must be > 0 and <= "// &
        & TRIM(NumberToVString(list%numberInList,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    item=list%listDP2(:,listItem)

    EXITS("List_ItemGetDP2")
    RETURN
999 ERRORSEXITS("List_ItemGetDP2",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemGetDP2
  
  !
  !================================================================================================================================
  !

  !>Determines if item is in the given integer list. If it is listItem is the index in the list. If not listItem is 0.
  SUBROUTINE List_ItemInListIntg1(list,item,listItem,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The pointer to the list
    INTEGER(INTG), INTENT(IN) :: item !<The item to find.
    INTEGER(INTG), INTENT(OUT) :: listItem !<On exit, the position of the item in the list. If the item does not exist then the value of 0 is returned.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    
    ENTERS("List_ItemInListIntg1",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsIntgData(list,err,error,*999)
    
    CALL List_Search(list,item,listItem,err,error,*999)

    EXITS("List_ItemInListIntg1")
    RETURN
999 ERRORSEXITS("List_ItemInListIntg1",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemInListIntg1
  
  !
  !================================================================================================================================
  !

  !>Determines if item is in the given integer list. If it is listItem is the index in the list. If not listItem is 0.
  SUBROUTINE List_ItemInListIntg2(list,items,listItem,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The pointer to the list
    INTEGER(INTG), INTENT(IN) :: items(:) !<The item to find.
    INTEGER(INTG), INTENT(OUT) :: listItem !<On exit, the position of the item in the list. If the item does not exist then the value of 0 is returned.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    
    ENTERS("List_ItemInListIntg2",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsIntgData(list,err,error,*999)
    
    CALL List_Search(list,items,listItem,err,error,*999)

    EXITS("List_ItemInListIntg2")
    RETURN
999 ERRORSEXITS("List_ItemInListIntg2",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemInListIntg2
  
  !
  !================================================================================================================================
  !

  !> Determines if item is in the given single precision real list. If it is listItem is the index in the list. If not
  !> listItem is 0.
  SUBROUTINE List_ItemInListSP1(list,item,listItem,err,error,*)

    !Argument Variables    
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The pointer to the list
    REAL(SP), INTENT(IN) :: item !<The item to find.
    INTEGER(INTG), INTENT(OUT) :: listItem !<On exit, the position of the item in the list. If the item does not exist then the value of 0 is returned.     
    INTEGER(INTG), INTENT(OUT) :: err !<The error code    
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("List_ItemInListSP1",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsSPData(list,err,error,*999)
    
    CALL List_Search(list,item,listItem,err,error,*999)

    EXITS("List_ItemInListSP1")
    RETURN
999 ERRORSEXITS("List_ItemInListSP1",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemInListSP1
  
  !
  !================================================================================================================================
  !

  !> Determines if item is in the given single precision real list. If it is listItem is the index in the list. If not
  !> listItem is 0.
  SUBROUTINE List_ItemInListSP2(list,items,listItem,err,error,*)

    !Argument Variables    
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The pointer to the list
    REAL(SP), INTENT(IN) :: items(:) !<The item to find.
    INTEGER(INTG), INTENT(OUT) :: listItem !<On exit, the position of the item in the list. If the item does not exist then the value of 0 is returned.     
    INTEGER(INTG), INTENT(OUT) :: err !<The error code    
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("List_ItemInListSP2",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsSPData(list,err,error,*999)
    
    CALL List_Search(list,items,listItem,err,error,*999)

    EXITS("List_ItemInListSP2")
    RETURN
999 ERRORSEXITS("List_ItemInListSP2",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemInListSP2
  
  !
  !================================================================================================================================
  !

  !> Determines if item is in the given double precision real list. If it is listItem is the index in the list. If not
  !> listItem is 0.
  SUBROUTINE List_ItemInListDP1(list,item,listItem,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The pointer to the list
    REAL(DP), INTENT(IN) :: item  !<The item to find.
    INTEGER(INTG), INTENT(OUT) :: listItem !<On exit, the position of the item in the list. If the item does not exist then the value of 0 is returned.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("List_ItemInListDP1",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsDPData(list,err,error,*999)
    
    CALL List_Search(list,item,listItem,err,error,*999)

    EXITS("List_ItemInListDP1")
    RETURN
999 ERRORSEXITS("List_ItemInListDP1",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemInListDP1

  !
  !================================================================================================================================
  !

  !> Determines if item is in the given double precision real list. If it is listItem is the index in the list. If not
  !> listItem is 0.
  SUBROUTINE List_ItemInListDP2(list,items,listItem,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The pointer to the list
    REAL(DP), INTENT(IN) :: items(:)  !<The item to find.
    INTEGER(INTG), INTENT(OUT) :: listItem !<On exit, the position of the item in the list. If the item does not exist then the value of 0 is returned.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("List_ItemInListDP2",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsDPData(list,err,error,*999)
    
    CALL List_Search(list,items,listItem,err,error,*999)

    EXITS("List_ItemInListDP2")
    RETURN
999 ERRORSEXITS("List_ItemInListDP2",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemInListDP2

  !
  !================================================================================================================================
  !
  
  !>Deletes the item given by the listItem index from the given list.
  SUBROUTINE List_ItemDelete(list,listItem,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The pointer to the list
    INTEGER(INTG), INTENT(IN) :: listItem !<The position in the list to delete.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("List_ItemDelete",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsIntgData(list,err,error,*999)
    IF(listItem<=0.OR.listItem>list%numberInList) THEN
      localError="The specified list item position of "//TRIM(NumberToVString(listItem,"*",err,error))// &
        & " is invalid. The list item position must be > 0 and <= "// &
        & TRIM(NumberToVString(list%numberInList,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    IF(list%dataDimension==1) THEN
      SELECT CASE(list%dataType)
      CASE(LIST_INTG_TYPE)
        list%listIntg(1:listItem-1)=list%listIntg(1:listItem-1)
        list%listIntg(listItem:list%numberInList-1)=list%listIntg(listItem+1:list%numberInList)
      CASE(LIST_SP_TYPE)
        list%listSP(1:listItem-1)=list%listSP(1:listItem-1)
        list%listSP(listItem:list%numberInList-1)=list%listSP(listItem+1:list%numberInList)
      CASE(LIST_DP_TYPE)
        list%listDP(1:listItem-1)=list%listDP(1:listItem-1)
        list%listDP(listItem:list%numberInList-1)=list%listDP(listItem+1:list%numberInList)
      CASE DEFAULT
        localError="The list data type of "//TRIM(NumberToVString(list%dataType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      SELECT CASE(list%dataType)
      CASE(LIST_INTG_TYPE)
        list%listIntg2(:,1:listItem-1)=list%listIntg2(:,1:listItem-1)
        list%listIntg2(:,listItem:list%numberInList-1)=list%listIntg2(:,listItem+1:list%numberInList)
      CASE(LIST_SP_TYPE)
        list%listSP2(:,1:listItem-1)=list%listSP2(:,1:listItem-1)
        list%listSP2(:,listItem:list%numberInList-1)=list%listSP2(:,listItem+1:list%numberInList)
      CASE(LIST_DP_TYPE)
        list%listDP2(:,1:listItem-1)=list%listDP2(:,1:listItem-1)
        list%listDP2(:,listItem:list%numberInList-1)=list%listDP2(:,listItem+1:list%numberInList)
      CASE DEFAULT
        localError="The list data type of "//TRIM(NumberToVString(list%dataType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    list%numberInList=list%numberInList-1

    EXITS("List_ItemDelete")
    RETURN
999 ERRORSEXITS("List_ItemDelete",err,error)
    RETURN 1
    
  END SUBROUTINE List_ItemDelete
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the key dimension (i.e., the dimension for searching and sorting) for a list
  SUBROUTINE List_KeyDimensionSet(list,keyDimension,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: keyDimension !<The key dimension to set. Must be greater than zero and <= the data dimension.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("List_KeyDimensionSet",err,error,*999)

    IF(.NOT.ASSOCIATED(list)) CALL FlagError("List is not associated.",err,error,*999)
    IF(keyDimension<=0.OR.keyDimension>list%dataDimension) THEN
      localError="The specified key dimension of "//TRIM(NumberToVString(keyDimension,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    list%keyDimension=keyDimension

    EXITS("List_KeyDimensionSet")
    RETURN
999 ERRORSEXITS("List_KeyDimensionSet",err,error)
    RETURN 1
    
  END SUBROUTINE List_KeyDimensionSet

  !
  !================================================================================================================================
  !

  !>Gets the current number of items in a list
  SUBROUTINE List_NumberOfItemsGet(list,numberOfItems,err,error,*)
      
    !Argument variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<A pointer to the list 
    INTEGER(INTG), INTENT(OUT) :: numberOfItems !<On exit, the current number of items in the list
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    
    ENTERS("List_NumberOfItemsGet",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    
    numberOfItems=list%numberInList
   
    EXITS("List_NumberOfItemsGet")
    RETURN
999 ERRORSEXITS("List_NumberOfItemsGet",err,error)
    RETURN 1
    
  END SUBROUTINE List_NumberOfItemsGet
  
  !
  !================================================================================================================================
  !

  !>Detaches the list values from an integer list of data dimension 1 and returns them as an array of base type
  !>before destroying the list. The listValues array must not be allocated on entry. It is up to the user to then deallocate
  !>the returned list memory.
  SUBROUTINE List_DetachAndDestroyIntg1(list,numberInList,listValues,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(INOUT) :: list !<The pointer to the list
    INTEGER(INTG), INTENT(OUT) :: numberInList !<On exit, the number in the list that has been detached.
    INTEGER(INTG), ALLOCATABLE, INTENT(INOUT) :: listValues(:) !<On exit, the detached list. Must not be allocated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("List_DetachAndDestroyIntg1",err,error,*999)

    IF(ALLOCATED(listValues)) CALL FlagError("List values is already allocated.",err,error,*999)
    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsIntgData(list,err,error,*999)
    IF(list%dataDimension/=1) THEN
      localError="Invalid data dimension. The supplied data dimension is 1 and the list data dimension is "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    numberInList=list%numberInList
    !Note this will return more memory as the list will be bigger. Maybe copy to an array the correct size?
    CALL MOVE_ALLOC(list%listIntg,listValues)
    CALL List_Finalise(list,err,error,*999)
    
    EXITS("List_DetachAndDestroyIntg1")
    RETURN
999 ERRORSEXITS("List_DetachAndDestroyIntg1",err,error)
    RETURN 1
    
  END SUBROUTINE List_DetachAndDestroyIntg1

  !
  !================================================================================================================================
  !

  !>Detaches the list values from an integer list of data dimension > 1 and returns them as an array of base type
  !>before destroying the list. The listValues array must not be allocated on entry. It is up to the user to then deallocate
  !>the returned list memory.
  SUBROUTINE List_DetachAndDestroyIntg2(list,numberInList,listValues,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(INOUT) :: list !<The pointer to the list
    INTEGER(INTG), INTENT(OUT) :: numberInList !<On exit, the number in the list that has been detached.
    INTEGER(INTG), ALLOCATABLE, INTENT(INOUT) :: listValues(:,:) !<On exit, the detached list. Must not be allocated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("List_DetachAndDestroyIntg2",err,error,*999)

    IF(ALLOCATED(listValues)) CALL FlagError("List values is already allocated.",err,error,*999)
    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsIntgData(list,err,error,*999)
    IF(list%dataDimension==1) THEN
      CALL FlagError("Invalid data dimension. The supplied data dimension is > 1 and the list data dimension is 1.", &
        & err,error,*999)
    ENDIF
    
    numberInList=list%numberInList
    !Note this will return more memory as the list will be bigger. Maybe copy to an array the correct size?
    CALL MOVE_ALLOC(list%listIntg2,listValues)
    CALL List_Finalise(list,err,error,*999)
    
    EXITS("List_DetachAndDestroyIntg2")
    RETURN
999 ERRORSEXITS("List_DetachAndDestroyIntg2",err,error)
    RETURN 1
    
  END SUBROUTINE List_DetachAndDestroyIntg2

  !
  !================================================================================================================================
  !

  !>Detaches the list values from a single precision real list of data dimension 1 and returns them as an array
  !>of base type before destroying the list. The listValues array must not be allocated on entry. It is up to the user to
  !>then deallocate the returned list memory.
  SUBROUTINE List_DetachAndDestroySP1(list,numberInList,listValues,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(INOUT) :: list !<The pointer to the list
    INTEGER(INTG), INTENT(OUT) :: numberInList !<On exit, the number in the list that has been detached.
    REAL(SP), ALLOCATABLE, INTENT(INOUT) :: listValues(:) !<On exit, the detached list. Must not be allocated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("List_DetachAndDestroySP1",err,error,*999)

    IF(ALLOCATED(listValues)) CALL FlagError("List values is already allocated.",err,error,*999)
    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsSPData(list,err,error,*999)
    IF(list%dataDimension/=1) THEN
      localError="Invalid data dimension. The supplied data dimension is 1 and the list data dimension is "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    numberInList=list%numberInList
    !Note this will return more memory as the list will be bigger. Maybe copy to an array the correct size?
    CALL MOVE_ALLOC(list%listSP,listValues)
    CALL List_Finalise(list,err,error,*999)
    
    EXITS("List_DetachAndDestroySP1")
    RETURN
999 ERRORSEXITS("List_DetachAndDestroySP1",err,error)
    RETURN 1
    
  END SUBROUTINE List_DetachAndDestroySP1
  
  !
  !================================================================================================================================
  !
  
  !>Detaches the list values from a single precision real list of data dimension > 1 and returns them as an array
  !>of base type before destroying the list. The listValues array must not be allocated on entry. It is up to the user to
  !>then deallocate the returned list memory.
  SUBROUTINE List_DetachAndDestroySP2(list,numberInList,listValues,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(INOUT) :: list !<The pointer to the list
    INTEGER(INTG), INTENT(OUT) :: numberInList !<On exit, the number in the list that has been detached.
    REAL(SP), ALLOCATABLE, INTENT(INOUT) :: listValues(:,:) !<On exit, the detached list. Must not be allocated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("List_DetachAndDestroySP2",err,error,*999)

    IF(ALLOCATED(listValues)) CALL FlagError("List values is already allocated.",err,error,*999)
    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsSPData(list,err,error,*999)
    IF(list%dataDimension/=1) THEN
      localError="Invalid data dimension. The supplied data dimension is 1 and the list data dimension is "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    numberInList=list%numberInList
    !Note this will return more memory as the list will be bigger. Maybe copy to an array the correct size?
    CALL MOVE_ALLOC(list%listSP2,listValues)
    CALL List_Finalise(list,err,error,*999)
    
    EXITS("List_DetachAndDestroySP2")
    RETURN
999 ERRORSEXITS("List_DetachAndDestroySP2",err,error)
    RETURN 1
    
  END SUBROUTINE List_DetachAndDestroySP2

  !
  !================================================================================================================================
  !
  
  !>Detaches the list values from a double precision real list of data dimension 1 and returns them as an array
  !>of base type before destroying the list. The listValues array must not be allocated on entry. It is up to the user
  !>to then deallocate the returned list memory.
  SUBROUTINE List_DetachAndDestroyDP1(list,numberInList,listValues,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(INOUT) :: list !<The pointer to the list
    INTEGER(INTG), INTENT(OUT) :: numberInList !<On exit, the number in the list that has been detached.
    REAL(DP), ALLOCATABLE, INTENT(INOUT) :: listValues(:) !<On exit, the detached list. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("List_DetachAndDestroyDP1",err,error,*999)

    IF(ALLOCATED(listValues)) CALL FlagError("List values is already allocated.",err,error,*999)
    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsDPData(list,err,error,*999)
    IF(list%dataDimension/=1) THEN
      localError="Invalid data dimension. The supplied data dimension is 1 and the list data dimension is "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    numberInList=list%numberInList
    !Note this will return more memory as the list will be bigger. Maybe copy to an array the correct size?
    CALL MOVE_ALLOC(list%listDP,listValues)
    CALL List_Finalise(list,err,error,*999)
    
    EXITS("List_DetachAndDestroyDP1")
    RETURN
999 ERRORSEXITS("List_DetachAndDestroyDP1",err,error)
    RETURN 1
    
  END SUBROUTINE List_DetachAndDestroyDP1

  !
  !================================================================================================================================
  !
  
  !>Detaches the list values from a double precision real list of data dimension > 1 and returns them as an array
  !>of base type before destroying the list. The listValues array must not be allocated on entry. It is up to the user
  !>to then deallocate the returned list memory.
  SUBROUTINE List_DetachAndDestroyDP2(list,numberInList,listValues,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(INOUT) :: list !<The pointer to the list
    INTEGER(INTG), INTENT(OUT) :: numberInList !<On exit, the number in the list that has been detached.
    REAL(DP), ALLOCATABLE, INTENT(INOUT) :: listValues(:,:) !<On exit, the detached list. Must not be allocated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("List_DetachAndDestroyDP2",err,error,*999)

    IF(ALLOCATED(listValues)) CALL FlagError("List values is already allocated.",err,error,*999)
    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsDPData(list,err,error,*999)
    IF(list%dataDimension==1) THEN
      CALL FlagError("Invalid data dimension. The supplied data dimension is > 1 and the list data dimension is 1.", &
        & err,error,*999)
    ENDIF
    
    numberInList=list%numberInList
    !Note this will return more memory as the list will be bigger. Maybe copy to an array the correct size?
    CALL MOVE_ALLOC(list%listDP2,listValues)
    CALL List_Finalise(list,err,error,*999)
    
    EXITS("List_DetachAndDestroyDP2")
    RETURN
999 ERRORSEXITS("List_DetachAndDestroyDP2",err,error)
    RETURN 1
    
  END SUBROUTINE List_DetachAndDestroyDP2

  !
  !================================================================================================================================
  !

  !>Removes duplicate entries from a list. A side effect of this is that the list is sorted.
  SUBROUTINE List_RemoveDuplicates(list,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(INOUT) :: list !<The pointer to the list
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i,j,numberRemoved
    LOGICAL :: sameValue
    TYPE(VARYING_STRING) :: localError

    ENTERS("List_RemoveDuplicates",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    
    IF(list%numberInList>0) THEN
      IF(list%dataDimension==1) THEN
        SELECT CASE(list%dataType)
        CASE(LIST_INTG_TYPE)              
          CALL List_Sort(list%listIntg(1:list%numberInList),err,error,*999)
          i=1
          DO WHILE(i<=list%numberInList)
            !Find the extent of duplicate values if any
            j=i+1
            sameValue=.TRUE.
            DO WHILE(j<=list%numberInList.AND.sameValue)
              IF(list%listIntg(j)==list%listIntg(i)) THEN
                j=j+1
              ELSE
                sameValue=.FALSE.
              ENDIF
            ENDDO !j
            IF(j>i+1.OR.sameValue) THEN
              !We have duplicates so remove them
              IF(sameValue) THEN
                !Duplicates to the end of the list so just set the number in the list
                list%numberInList=i
              ELSE
                numberRemoved=j-i-1
                list%listIntg(i+1:list%numberInList-numberRemoved)=list%listIntg(j:list%numberInList)
                list%numberInList=list%numberInList-numberRemoved
              ENDIF
            ENDIF
            i=i+1
          ENDDO !i
        CASE(LIST_SP_TYPE)
          CALL List_Sort(list%listSP(1:list%numberInList),err,error,*999)
          i=1
          DO WHILE(i<=list%numberInList)
            !Find the extent of duplicate values if any
            j=i+1
            sameValue=.TRUE.
            DO WHILE(j<=list%numberInList.AND.sameValue)
              IF(ABS(list%listSP(j)-list%listSP(i))<=ZERO_TOLERANCE_SP) THEN
                j=j+1
              ELSE
                sameValue=.FALSE.
              ENDIF
            ENDDO !j
            IF(j>i+1.OR.sameValue) THEN
              !We have duplicates so remove them
              IF(sameValue) THEN
                !Duplicates to the end of the list so just set the number in the list
                list%numberInList=i
              ELSE
                numberRemoved=j-i-1
                list%listSP(i+1:list%numberInList-numberRemoved)=list%listSP(j:list%numberInList)
                list%numberInList=list%numberInList-numberRemoved
              ENDIF
            ENDIF
            i=i+1
          ENDDO !i
        CASE(LIST_DP_TYPE)
          CALL List_Sort(list%listDP(1:list%numberInList),err,error,*999)
          i=1
          DO WHILE(i<=list%numberInList)
            !Find the extent of duplicate values if any
            j=i+1
            sameValue=.TRUE.
            DO WHILE(j<=list%numberInList.AND.sameValue)
              IF(ABS(list%listDP(j)-list%listDP(i))<=ZERO_TOLERANCE) THEN
                j=j+1
              ELSE
                sameValue=.FALSE.
              ENDIF
            ENDDO !j
            IF(j>i+1.OR.sameValue) THEN
              !We have duplicates so remove them
              IF(sameValue) THEN
                !Duplicates to the end of the list so just set the number in the list
                list%numberInList=i
              ELSE
                numberRemoved=j-i-1
                list%listDP(i+1:list%numberInList-numberRemoved)=list%listDP(j:list%numberInList)
                list%numberInList=list%numberInList-numberRemoved
              ENDIF
            ENDIF
            i=i+1
          ENDDO !i
        CASE DEFAULT
          localError="The list data type of "//TRIM(NumberToVString(list%dataType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        SELECT CASE(list%dataType)
        CASE(LIST_INTG_TYPE)              
          CALL List_Sort(list%listIntg2(:,1:list%numberInList),list%keyDimension,err,error,*999)
          i=1
          DO WHILE(i<=list%numberInList)
            !Find the extent of duplicate values if any
            j=i+1
            sameValue=.TRUE.
            DO WHILE(j<=list%numberInList.AND.sameValue)
              IF(list%listIntg2(list%keyDimension,j)==list%listIntg2(list%keyDimension,i)) THEN
                j=j+1
              ELSE
                sameValue=.FALSE.
              ENDIF
            ENDDO !j
            IF(j>i+1.OR.sameValue) THEN
              !We have duplicates so remove them
              IF(sameValue) THEN
                !Duplicates to the end of the list so just set the number in the list
                list%numberInList=i
              ELSE
                numberRemoved=j-i-1
                list%listIntg2(:,i+1:list%numberInList-numberRemoved)=list%listIntg2(:,j:list%numberInList)
                list%numberInList=list%numberInList-numberRemoved
              ENDIF
            ENDIF
            i=i+1
          ENDDO !i
        CASE(LIST_SP_TYPE)
          CALL List_Sort(list%listSP2(:,1:list%numberInList),list%keyDimension,err,error,*999)
          i=1
          DO WHILE(i<=list%numberInList)
            !Find the extent of duplicate values if any
            j=i+1
            sameValue=.TRUE.
            DO WHILE(j<=list%numberInList.AND.sameValue)
              IF(ABS(list%listSP2(list%keyDimension,j)-list%listSP2(list%keyDimension,i))<=ZERO_TOLERANCE_SP) THEN
                j=j+1
              ELSE
                sameValue=.FALSE.
              ENDIF
            ENDDO !j
            IF(j>i+1.OR.sameValue) THEN
              !We have duplicates so remove them
              IF(sameValue) THEN
                !Duplicates to the end of the list so just set the number in the list
                list%numberInList=i
              ELSE
                numberRemoved=j-i-1
                list%listSP2(:,i+1:list%numberInList-numberRemoved)=list%listSP2(:,j:list%numberInList)
                list%numberInList=list%numberInList-numberRemoved
              ENDIF
            ENDIF
            i=i+1
          ENDDO !i
        CASE(LIST_DP_TYPE)
          CALL List_Sort(list%listDP2(:,1:list%numberInList),list%keyDimension,err,error,*999)
          i=1
          DO WHILE(i<=list%numberInList)
            !Find the extent of duplicate values if any
            j=i+1
            sameValue=.TRUE.
            DO WHILE(j<=list%numberInList.AND.sameValue)
              IF(ABS(list%listDP2(list%keyDimension,j)-list%listDP2(list%keyDimension,i))<=ZERO_TOLERANCE) THEN
                j=j+1
              ELSE
                sameValue=.FALSE.
              ENDIF
            ENDDO !j
            IF(j>i+1.OR.sameValue) THEN
              !We have duplicates so remove them
              IF(sameValue) THEN
                !Duplicates to the end of the list so just set the number in the list
                list%numberInList=i
              ELSE
                numberRemoved=j-i-1
                list%listDP2(:,i+1:list%numberInList-numberRemoved)=list%listDP2(:,j:list%numberInList)
                list%numberInList=list%numberInList-numberRemoved
              ENDIF
            ENDIF
            i=i+1
          ENDDO !i
        CASE DEFAULT
          localError="The list data type of "//TRIM(NumberToVString(list%dataType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDIF
    ENDIF
  
    EXITS("List_RemoveDuplicates")
    RETURN
999 ERRORSEXITS("List_RemoveDuplicates",err,error)
    RETURN 1
    
  END SUBROUTINE List_RemoveDuplicates
 
  !
  !================================================================================================================================
  !

  !>Sets/changes the searching method for a list
  SUBROUTINE List_SearchingMethodSet(list,searchingMethod,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<A pointer to the list to set the searching method for
    INTEGER(INTG), INTENT(IN) :: searchingMethod !<The searching method to set. See \see Lists_SeachingMethod,Lists
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("List_SearchingMethodSet",err,error,*999)

    CALL List_AssertNotFinished(list,err,error,*999)

    SELECT CASE(searchingMethod)
    CASE(LIST_LINEAR_SEARCH_METHOD)
      list%searchMethod=LIST_LINEAR_SEARCH_METHOD
    CASE(LIST_BINARY_SEARCH_METHOD)
      list%searchMethod=LIST_BINARY_SEARCH_METHOD
    CASE DEFAULT
      localError="The specified list searching method of "//TRIM(NumberToVString(searchingMethod,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("List_SearchingMethodSet")
    RETURN
999 ERRORSEXITS("List_SearchingMethodSet",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchingMethodSet

  !
  !================================================================================================================================
  !

  !>Seaches an integer list to see if it contains the value and returns the position in the list of that value. If value is not found the position will be zero.
  SUBROUTINE List_SearchListIntg1(list,value,position,err,error,*)
  
    !Argument variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The list to search
    INTEGER(INTG), INTENT(IN) :: value !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On return, the position of value in the list. If value is not in the list then position will be zero.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_SearchListIntg1",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsIntgData(list,err,error,*999)
    IF(list%dataDimension/=1) THEN
      localError="Invalid data dimension. The supplied value has dimension 1 but the supplied list has data dimension "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
   
    SELECT CASE(list%searchMethod)
    CASE(LIST_LINEAR_SEARCH_METHOD)
      CALL List_SearchLinearIntg1Array(list%listIntg(1:list%numberInList),value,position,err,error,*999)
    CASE(LIST_BINARY_SEARCH_METHOD)
      CALL List_SearchBinaryIntg1Array(list%listIntg(1:list%numberInList),value,position,err,error,*999)
    CASE DEFAULT
      localError="The list search method of "//TRIM(NumberToVString(list%searchMethod,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("List_SearchListIntg1")
    RETURN
999 ERRORSEXITS("List_SearchListIntg1",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchListIntg1
  
  !
  !================================================================================================================================
  !

  !>Seaches an integer list to see if it contains the value and returns the position in the list of that value. If value is not found the position will be zero.
  SUBROUTINE List_SearchListIntg2(list,values,position,err,error,*)
  
    !Argument variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The list to search
    INTEGER(INTG), INTENT(IN) :: values(:) !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On return, the position of value in the list. If value is not in the list then position will be zero.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_SearchListIntg2",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsIntgData(list,err,error,*999)
    IF(list%dataDimension==1) THEN
      localError="Invalid data dimension. The supplied list has data dimension 1 but the supplied value has data dimension "// &
        & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
   
    SELECT CASE(list%searchMethod)
    CASE(LIST_LINEAR_SEARCH_METHOD)
      CALL List_SearchLinearIntg2Array(list%listIntg2(:,1:list%numberInList),list%keyDimension,values,position,err,error,*999)
    CASE(LIST_BINARY_SEARCH_METHOD)
      CALL List_SearchBinaryIntg2Array(list%listIntg2(:,1:list%numberInList),list%keyDimension,values,position,err,error,*999)
    CASE DEFAULT
      localError="The list search method of "//TRIM(NumberToVString(list%searchMethod,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("List_SearchListIntg2")
    RETURN
999 ERRORSEXITS("List_SearchListIntg2",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchListIntg2
  
  !
  !================================================================================================================================
  !

  !>Seaches an single precision real list to see if it contains the value and returns the position in the list of that value. If value is not found the position will be zero.
  SUBROUTINE List_SearchListSP1(list,value,position,err,error,*)
  
    !Argument variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The list to search
    REAL(SP), INTENT(IN) :: value !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On return, the position of value in the list. If value is not in the list then position will be zero.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_SearchListSP1",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsSPData(list,err,error,*999)
    IF(list%dataDimension/=1) THEN
      localError="Invalid data dimension. The supplied value has dimension 1 but the supplied list has data dimension "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(list%searchMethod)
    CASE(LIST_LINEAR_SEARCH_METHOD)
      CALL List_SearchLinearSP1Array(list%listSP(1:list%numberInList),value,position,err,error,*999)
    CASE(LIST_BINARY_SEARCH_METHOD)
      CALL List_SearchBinarySP1Array(list%listSP(1:list%numberInList),value,position,err,error,*999)
    CASE DEFAULT
      localError="The list search method of "//TRIM(NumberToVString(list%searchMethod,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("List_SearchListSP1")
    RETURN
999 ERRORSEXITS("List_SearchListSP1",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchListSP1
  
  !
  !================================================================================================================================
  !

  !>Seaches an single precision real list to see if it contains the value and returns the position in the list of that value. If value is not found the position will be zero.
  SUBROUTINE List_SearchListSP2(list,values,position,err,error,*)
  
    !Argument variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The list to search
    REAL(SP), INTENT(IN) :: values(:) !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On return, the position of value in the list. If value is not in the list then position will be zero.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_SearchListSP2",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsSPData(list,err,error,*999)
    IF(list%dataDimension==1) THEN
      localError="Invalid data dimension. The supplied list has data dimension 1 but the supplied value has data dimension "// &
        & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(list%searchMethod)
    CASE(LIST_LINEAR_SEARCH_METHOD)
      CALL List_SearchLinearSP2Array(list%listSP2(:,1:list%numberInList),list%keyDimension,values,position,err,error,*999)
    CASE(LIST_BINARY_SEARCH_METHOD)
      CALL List_SearchBinarySP2Array(list%listSP2(:,1:list%numberInList),list%keyDimension,values,position,err,error,*999)
    CASE DEFAULT
      localError="The list search method of "//TRIM(NumberToVString(list%searchMethod,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("List_SearchListSP2")
    RETURN
999 ERRORSEXITS("List_SearchListSP2",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchListSP2
  
  !
  !================================================================================================================================
  !

  !>Seaches an double precision real list to see if it contains the value and returns the position in the list of that value. If value is not found the position will be zero.
  SUBROUTINE List_SearchListDP1(list,value,position,err,error,*)
  
    !Argument variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The list to search
    REAL(DP), INTENT(IN) :: value !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On return, the position of value in the list. If value is not in the list then position will be zero.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_SearchListDP1",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsDPData(list,err,error,*999)
     IF(list%dataDimension/=1) THEN
      localError="Invalid data dimension. The supplied value has dimension 1 but the supplied list has data dimension "// &
        & TRIM(NumberToVString(list%dataDimension,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
   
    SELECT CASE(list%searchMethod)
    CASE(LIST_LINEAR_SEARCH_METHOD)
      CALL List_SearchLinearDP1Array(list%listDP(1:list%numberInList),value,position,err,error,*999)
    CASE(LIST_BINARY_SEARCH_METHOD)
      CALL List_SearchBinaryDP1Array(list%listDP(1:list%numberInList),value,position,err,error,*999)
    CASE DEFAULT
      localError="The list search method of "//TRIM(NumberToVString(list%searchMethod,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("List_SearchListDP1")
    RETURN
999 ERRORSEXITS("List_SearchListDP1",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchListDP1
  
  !
  !================================================================================================================================
  !

  !>Seaches an double precision real list to see if it contains the value and returns the position in the list of that value. If value is not found the position will be zero.
  SUBROUTINE List_SearchListDP2(list,values,position,err,error,*)
  
    !Argument variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<The list to search
    REAL(DP), INTENT(IN) :: values(:) !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On return, the position of value in the list. If value is not in the list then position will be zero.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_SearchListDP2",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    CALL List_AssertIsDPData(list,err,error,*999)
    IF(list%dataDimension==1) THEN
      localError="Invalid data dimension. The supplied list has data dimension 1 but the supplied value has data dimension "// &
        & TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
   
    SELECT CASE(list%searchMethod)
    CASE(LIST_LINEAR_SEARCH_METHOD)
      CALL List_SearchLinearDP2Array(list%listDP2(:,1:list%numberInList),list%keyDimension,values,position,err,error,*999)
    CASE(LIST_BINARY_SEARCH_METHOD)
      CALL List_SearchBinaryDP2Array(list%listDP2(:,1:list%numberInList),list%keyDimension,values,position,err,error,*999)
    CASE DEFAULT
      localError="The list search method of "//TRIM(NumberToVString(list%searchMethod,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("List_SearchListDP2")
    RETURN
999 ERRORSEXITS("List_SearchListDP2",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchListDP2
  
  !
  !================================================================================================================================
  !

  !>Searches an integer array list A for value. If the search is successful position contains the index of the position of value in the list otherwise position is zero.
  SUBROUTINE List_SearchIntg1Array(a,value,position,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: a(:) !<The list to search
    INTEGER(INTG), INTENT(IN) :: value !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    
    ENTERS("List_SearchIntg1Array",err,error,*999)

    !Default search method is a linear search
    CALL List_SearchLinear(a,value,position,err,error,*999)
    
    EXITS("List_SearchIntg1Array")
    RETURN
999 ERRORSEXITS("List_SearchIntg1Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchIntg1Array
  
  !
  !================================================================================================================================
  !

  !>Searches an integer array list A for value. If the search is successful position contains the index of the position of value in the list otherwise position is zero.
  SUBROUTINE List_SearchIntg2Array(a,keyDimension,values,position,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: a(:,:) !<The list to search
    INTEGER(INTG), INTENT(IN) :: keyDimension !<The key dimension to search
    INTEGER(INTG), INTENT(IN) :: values(:) !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    
    ENTERS("List_SearchIntg2Array",err,error,*999)

    !Default search method is a linear search
    CALL List_SearchLinear(a,keyDimension,values,position,err,error,*999)
    
    EXITS("List_SearchIntg2Array")
    RETURN
999 ERRORSEXITS("List_SearchIntg2Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchIntg2Array
  
  !
  !================================================================================================================================
  !

  !>Searches a single precision real array list a for value. If the search is successful position contains the index of the position of value in the list otherwise position is zero.
  SUBROUTINE List_SearchSP1Array(a,value,position,err,error,*)
  
   !Argument variables
    REAL(SP), INTENT(IN) :: a(:) !<The list to search
    REAL(SP), INTENT(IN) :: value !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
   
    ENTERS("List_SearchSP1Array",err,error,*999)

    !Default search method is a linear search
    CALL List_SearchLinear(a,value,position,err,error,*999)    
    
    EXITS("List_SearchSP1Array")
    RETURN
999 ERRORSEXITS("List_SearchSP1Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchSP1Array
  
  !
  !================================================================================================================================
  !

  !>Searches a single precision real array list a for value. If the search is successful position contains the index of the position of value in the list otherwise position is zero.
  SUBROUTINE List_SearchSP2Array(a,keyDimension,values,position,err,error,*)
  
   !Argument variables
    REAL(SP), INTENT(IN) :: a(:,:) !<The list to search
    INTEGER(INTG), INTENT(IN) :: keyDimension !<The key dimension to search
    REAL(SP), INTENT(IN) :: values(:) !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
   
    ENTERS("List_SearchSP2Array",err,error,*999)

    !Default search method is a linear search
    CALL List_SearchLinear(a,keyDimension,values,position,err,error,*999)    
    
    EXITS("List_SearchSP2Array")
    RETURN
999 ERRORSEXITS("List_SearchSP2Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchSP2Array
  
  !
  !================================================================================================================================
  !

  !>Searches a double precision real array list a for value. If the search is successful position contains the index of the position of value in the list otherwise position is zero.
  SUBROUTINE List_SearchDP1Array(a,value,position,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: a(:) !<The list to search
    REAL(DP), INTENT(IN) :: value !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    
    ENTERS("List_SearchDP1Array",err,error,*999)

    !Default search method is a linear search
    CALL List_SearchLinear(a,value,position,err,error,*999)    
    
    EXITS("List_SearchDP1Array")
    RETURN
999 ERRORSEXITS("List_SearchDP1Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchDP1Array
  
  !
  !================================================================================================================================
  !

  !>Searches a double precision real array list a for value. If the search is successful position contains the index of the position of value in the list otherwise position is zero.
  SUBROUTINE List_SearchDP2Array(a,keyDimension,values,position,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: a(:,:) !<The list to search
    INTEGER(INTG), INTENT(IN) :: keyDimension !<The key dimension to search
    REAL(DP), INTENT(IN) :: values(:) !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    
    ENTERS("List_SearchDP2Array",err,error,*999)

    !Default search method is a linear search
    CALL List_SearchLinear(a,keyDimension,values,position,err,error,*999)    
    
    EXITS("List_SearchDP2Array")
    RETURN
999 ERRORSEXITS("List_SearchDP2Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchDP2Array
  
  !
  !================================================================================================================================
  !

  !>Searches an integer array list a for value using the binary search method. If the search is successful position contains the index of the position of value in the list otherwise position is zero. This routine will assume that the list a is sorted into ascending order
  SUBROUTINE List_SearchBinaryIntg1Array(a,value,position,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: a(:) !<The list to search
    INTEGER(INTG), INTENT(IN) :: value !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: lowLimit,upLimit,midPoint
    
    ENTERS("List_SearchBinaryIntg1Array",err,error,*999)

    position=0
    lowLimit=1
    upLimit=SIZE(a,1)
    IF(a(lowLimit)==VALUE) THEN
      position=lowLimit
    ELSE IF(a(upLimit)==VALUE) THEN
      position=upLimit
    ELSE
      DO WHILE(lowLimit<=upLimit)
        midPoint=(upLimit+lowLimit)/2
        IF(a(midPoint)==VALUE) THEN
          position=midPoint
          EXIT
        ELSE IF(a(midPoint)>value) THEN
          upLimit=midPoint
        ELSE
          lowLimit=midPoint
        ENDIF
      ENDDO
    ENDIF
     
    EXITS("List_SearchBinaryIntg1Array")
    RETURN
999 ERRORSEXITS("List_SearchBinaryIntg1Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchBinaryIntg1Array
  
  !
  !================================================================================================================================
  !

  !>Searches an integer array list a for value using the binary search method. If the search is successful position contains the index of the position of value in the list otherwise position is zero. This routine will assume that the list a is sorted into ascending order
  SUBROUTINE List_SearchBinaryIntg2Array(a,keyDimension,values,position,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: a(:,:) !<The list to search
    INTEGER(INTG), INTENT(IN) :: keyDimension !<The key dimension to search
    INTEGER(INTG), INTENT(IN) :: values(:) !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: lowLimit,upLimit,midPoint
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_SearchBinaryIntg2Array",err,error,*999)

    IF(keyDimension<=0.OR.keyDimension>SIZE(a,1)) THEN
      localError="The specified key dimension of "//TRIM(NumberToVString(keyDimension,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(a,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(a,1)/=SIZE(values,1)) THEN
      localError="The size of the array data dimension of "//TRIM(NumberToVString(SIZE(a,1),"*",err,error))// &
        & " does not match the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    position=0
    lowLimit=1
    upLimit=SIZE(a,2)
    IF(a(keyDimension,lowLimit)==values(keyDimension)) THEN
      position=lowLimit
    ELSE IF(a(keyDimension,upLimit)==values(keyDimension)) THEN
      position=upLimit
    ELSE
      DO WHILE(lowLimit<=upLimit)
        midPoint=(upLimit+lowLimit)/2
        IF(a(keyDimension,midPoint)==values(keyDimension)) THEN
          position=midPoint
          EXIT
        ELSE IF(a(keyDimension,midPoint)>values(keyDimension)) THEN
          upLimit=midPoint
        ELSE
          lowLimit=midPoint
        ENDIF
      ENDDO
    ENDIF
     
    EXITS("List_SearchBinaryIntg2Array")
    RETURN
999 ERRORSEXITS("List_SearchBinaryIntg2Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchBinaryIntg2Array
  
  !
  !================================================================================================================================
  !

  !>Searches a single precision real array list a for value using the binary search method. If the search is successful position contains the index of the position of value in the list otherwise position is zero. This routine will assume that the list a is sorted into ascending order
  SUBROUTINE List_SearchBinarySP1Array(a,value,position,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: a(:) !<The list to search
    REAL(SP), INTENT(IN) :: value !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: lowLimit,upLimit,midPoint
    
    ENTERS("List_SearchBinarySP1Array",err,error,*999)

    position=0
    lowLimit=1
    upLimit=SIZE(a,1)
    IF(ABS(a(lowLimit)-VALUE)<ZERO_TOLERANCE_SP) THEN
      position=lowLimit
    ELSE IF(ABS(a(upLimit)-VALUE)<ZERO_TOLERANCE_SP) THEN
      position=upLimit
    ELSE
      DO WHILE(lowLimit<=upLimit)
        midPoint=(upLimit+lowLimit)/2
        IF(ABS(a(midPoint)-VALUE)<ZERO_TOLERANCE_SP) THEN
          position=midPoint
          EXIT
        ELSE IF(a(midPoint)>value) THEN
          upLimit=midPoint
        ELSE
          lowLimit=midPoint
        ENDIF
      ENDDO
    ENDIF
     
    EXITS("List_SearchBinarySP1Array")
    RETURN
999 ERRORSEXITS("List_SearchBinarySP1Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchBinarySP1Array
  
  !
  !================================================================================================================================
  !

  !>Searches a single precision real array list a for value using the binary search method. If the search is successful position contains the index of the position of value in the list otherwise position is zero. This routine will assume that the list a is sorted into ascending order
  SUBROUTINE List_SearchBinarySP2Array(a,keyDimension,values,position,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: a(:,:) !<The list to search
    INTEGER(INTG), INTENT(IN) :: keyDimension !<The key dimension to search
    REAL(SP), INTENT(IN) :: values(:) !<The values to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: lowLimit,upLimit,midPoint
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_SearchBinarySP2Array",err,error,*999)

    IF(keyDimension<=0.OR.keyDimension>SIZE(a,1)) THEN
      localError="The specified key dimension of "//TRIM(NumberToVString(keyDimension,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(a,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(a,1)/=SIZE(values,1)) THEN
      localError="The size of the array data dimension of "//TRIM(NumberToVString(SIZE(a,1),"*",err,error))// &
        & " does not match the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    position=0
    lowLimit=1
    upLimit=SIZE(a,2)
    IF(ABS(a(keyDimension,lowLimit)-values(keyDimension))<ZERO_TOLERANCE_SP) THEN
      position=lowLimit
    ELSE IF(ABS(a(keyDimension,upLimit)-values(keyDimension))<ZERO_TOLERANCE_SP) THEN
      position=upLimit
    ELSE
      DO WHILE(lowLimit<=upLimit)
        midPoint=(upLimit+lowLimit)/2
        IF(ABS(a(keyDimension,midPoint)-values(keyDimension))<ZERO_TOLERANCE_SP) THEN
          position=midPoint
          EXIT
        ELSE IF(a(keyDimension,midPoint)>values(keyDimension)) THEN
          upLimit=midPoint
        ELSE
          lowLimit=midPoint
        ENDIF
      ENDDO
    ENDIF
     
    EXITS("List_SearchBinarySP2Array")
    RETURN
999 ERRORSEXITS("List_SearchBinarySP2Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchBinarySP2Array
  
  !
  !================================================================================================================================
  !

  !>Searches a double precision real array list a for value using the binary search method. If the search is successful position contains the index of the position of value in the list otherwise position is zero. This routine will assume that the list a is sorted into ascending order
  SUBROUTINE List_SearchBinaryDP1Array(a,value,position,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: a(:) !<The list to search
    REAL(DP), INTENT(IN) :: value !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: lowLimit,upLimit,midPoint
    
    ENTERS("List_SearchBinaryDP1Array",err,error,*999)

    position=0
    lowLimit=1
    upLimit=SIZE(a,1)
    IF(ABS(a(lowLimit)-VALUE)<ZERO_TOLERANCE_DP) THEN
      position=lowLimit
    ELSE IF(ABS(a(upLimit)-VALUE)<ZERO_TOLERANCE_DP) THEN
      position=upLimit
    ELSE
      DO WHILE(lowLimit<=upLimit)
        midPoint=(upLimit+lowLimit)/2
        IF(ABS(a(midPoint)-VALUE)<ZERO_TOLERANCE_DP) THEN
          position=midPoint
          EXIT
        ELSE IF(a(midPoint)>value) THEN
          upLimit=midPoint
        ELSE
          lowLimit=midPoint
        ENDIF
      ENDDO
    ENDIF
      
    EXITS("List_SearchBinaryDP1Array")
    RETURN
999 ERRORSEXITS("List_SearchBinaryDP1Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchBinaryDP1Array
  
  !
  !================================================================================================================================
  !

  !>Searches a double precision real array list a for value using the binary search method. If the search is successful position contains the index of the position of value in the list otherwise position is zero. This routine will assume that the list a is sorted into ascending order
  SUBROUTINE List_SearchBinaryDP2Array(a,keyDimension,values,position,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: a(:,:) !<The list to search
    INTEGER(INTG), INTENT(IN) :: keyDimension !<The key dimension to search
    REAL(DP), INTENT(IN) :: values(:) !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: lowLimit,upLimit,midPoint
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_SearchBinaryDP2Array",err,error,*999)

    IF(keyDimension<=0.OR.keyDimension>SIZE(a,1)) THEN
      localError="The specified key dimension of "//TRIM(NumberToVString(keyDimension,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(a,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(a,1)/=SIZE(values,1)) THEN
      localError="The size of the array data dimension of "//TRIM(NumberToVString(SIZE(a,1),"*",err,error))// &
        & " does not match the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    position=0
    lowLimit=1
    upLimit=SIZE(a,2)
    IF(ABS(a(keyDimension,lowLimit)-values(keyDimension))<ZERO_TOLERANCE_DP) THEN
      position=lowLimit
    ELSE IF(ABS(a(keyDimension,upLimit)-values(keyDimension))<ZERO_TOLERANCE_DP) THEN
      position=upLimit
    ELSE
      DO WHILE(lowLimit<=upLimit)
        midPoint=(upLimit+lowLimit)/2
        IF(ABS(a(keyDimension,midPoint)-values(keyDimension))<ZERO_TOLERANCE_DP) THEN
          position=midPoint
          EXIT
        ELSE IF(a(keyDimension,midPoint)>values(keyDimension)) THEN
          upLimit=midPoint
        ELSE
          lowLimit=midPoint
        ENDIF
      ENDDO
    ENDIF
     
    EXITS("List_SearchBinaryDP2Array")
    RETURN
999 ERRORSEXITS("List_SearchBinaryDP2Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchBinaryDP2Array
  
  !
  !================================================================================================================================
  !

  !>Searches an integer array list A for value using the linear search method. If the search is successful position contains the index of the position of value in the list otherwise position is zero.
  SUBROUTINE List_SearchLinearIntg1Array(a,value,position,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: a(:) !<The list to search
    INTEGER(INTG), INTENT(IN) :: value !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i
    LOGICAL :: found
    
    ENTERS("List_SearchLinearIntg1Array",err,error,*999)

    found=.FALSE.
    i=1
    DO WHILE(i<=SIZE(a,1).AND..NOT.found)
      IF(a(i)==value) THEN
        found=.TRUE.
      ELSE
        i=i+1
      ENDIF
    ENDDO
    IF(found) THEN
      position=i
    ELSE
      position=0
    ENDIF
    
    EXITS("List_SearchLinearIntg1Array")
    RETURN
999 ERRORSEXITS("List_SearchLinearIntg1Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchLinearIntg1Array
  
  !
  !================================================================================================================================
  !

  !>Searches an integer array list A for value using the linear search method. If the search is successful position contains the index of the position of value in the list otherwise position is zero.
  SUBROUTINE List_SearchLinearIntg2Array(a,keyDimension,values,position,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: a(:,:) !<The list to search
    INTEGER(INTG), INTENT(IN) :: keyDimension !<The key dimension to search
    INTEGER(INTG), INTENT(IN) :: values(:) !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i
    LOGICAL :: found
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_SearchLinearIntg2Array",err,error,*999)

    IF(keyDimension<=0.OR.keyDimension>SIZE(a,1)) THEN
      localError="The specified key dimension of "//TRIM(NumberToVString(keyDimension,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(a,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(a,1)/=SIZE(values,1)) THEN
      localError="The size of the array data dimension of "//TRIM(NumberToVString(SIZE(a,1),"*",err,error))// &
        & " does not match the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    found=.FALSE.
    i=1
    DO WHILE(i<=SIZE(a,2).AND..NOT.found)
      IF(a(keyDimension,i)==values(keyDimension)) THEN
        found=.TRUE.
      ELSE
        i=i+1
      ENDIF
    ENDDO  
    IF(found) THEN
      position=i
    ELSE
      position=0
    ENDIF
    
    EXITS("List_SearchLinearIntg2Array")
    RETURN
999 ERRORSEXITS("List_SearchLinearIntg2Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchLinearIntg2Array
  
  !
  !================================================================================================================================
  !

  !>Searches a single precision real array list a for value using the linear search method. If the search is successful position contains the index of the position of value in the list otherwise position is zero.
  SUBROUTINE List_SearchLinearSP1Array(A,value,position,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: a(:) !<The list to search
    REAL(SP), INTENT(IN) :: value !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i
    LOGICAL :: found
    
    ENTERS("List_SearchLinearSP1Array",err,error,*999)

    found=.FALSE.
    i=1
    DO WHILE(i<=SIZE(a,1).AND..NOT.found)
      IF(ABS(a(i)-value)<ZERO_TOLERANCE_SP) THEN
        found=.TRUE.
      ELSE
        i=i+1
      ENDIF
    ENDDO
    IF(found) THEN
      position=i
    ELSE
      position=0
    ENDIF
    
    EXITS("List_SearchLinearSP1Array")
    RETURN
999 ERRORSEXITS("List_SearchLinearSP1Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchLinearSP1Array
  
  !
  !================================================================================================================================
  !

  !>Searches a single precision real array list a for value using the linear search method. If the search is successful position contains the index of the position of value in the list otherwise position is zero.
  SUBROUTINE List_SearchLinearSP2Array(a,keyDimension,values,position,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: a(:,:) !<The list to search
    INTEGER(INTG), INTENT(IN) :: keyDimension !<The key dimension to search
    REAL(SP), INTENT(IN) :: values(:) !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i
    LOGICAL :: found
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_SearchLinearSP2Array",err,error,*999)

    IF(keyDimension<=0.OR.keyDimension>SIZE(a,1)) THEN
      localError="The specified key dimension of "//TRIM(NumberToVString(keyDimension,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(a,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(a,1)/=SIZE(values,1)) THEN
      localError="The size of the array data dimension of "//TRIM(NumberToVString(SIZE(a,1),"*",err,error))// &
        & " does not match the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    found=.FALSE.
    i=1
    DO WHILE(i<=SIZE(a,2).AND..NOT.found)
      IF(ABS(a(keyDimension,i)-values(keyDimension))<ZERO_TOLERANCE_SP) THEN
        found=.TRUE.
      ELSE
        i=i+1
      ENDIF
    ENDDO
    IF(found) THEN
      position=i
    ELSE
      position=0
    ENDIF
    
    EXITS("List_SearchLinearSP2Array")
    RETURN
999 ERRORSEXITS("List_SearchLinearSP2Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchLinearSP2Array
  
  !
  !================================================================================================================================
  !

  !>Searches a double precision real array list A for value using the linear search method. If the search is successful position contains the index of the position of value in the list otherwise position is zero.
  SUBROUTINE List_SearchLinearDP1Array(a,value,position,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: a(:) !<The list to search
    REAL(DP), INTENT(IN) :: value !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i
    LOGICAL :: found
    
    ENTERS("List_SearchLinearDP1Array",err,error,*999)

    found=.FALSE.
    i=1
    DO WHILE(i<=SIZE(a,1).AND..NOT.found)
      IF(ABS(a(i)-value)<ZERO_TOLERANCE_DP) THEN
        found=.TRUE.
      ELSE
        i=i+1
      ENDIF
    ENDDO
    IF(found) THEN
      position=i
    ELSE
      position=0
    ENDIF
    
    EXITS("List_SearchLinearDP1Array")
    RETURN
999 ERRORSEXITS("List_SearchLinearDP1Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchLinearDP1Array
  
  !
  !================================================================================================================================
  !

  !>Searches a double precision real array list a for value using the linear search method. If the search is successful position contains the index of the position of value in the list otherwise position is zero.
  SUBROUTINE List_SearchLinearDP2Array(a,keyDimension,values,position,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: a(:,:) !<The list to search
    INTEGER(INTG), INTENT(IN) :: keyDimension !<The key dimension to search
    REAL(DP), INTENT(IN) :: values(:) !<The value to search for
    INTEGER(INTG), INTENT(OUT) :: position !<On exit, the position of value in the list. If value does not exist in the list the returned position is zero. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i
    LOGICAL :: found
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_SearchLinearDP2Array",err,error,*999)

    IF(keyDimension<=0.OR.keyDimension>SIZE(a,1)) THEN
      localError="The specified key dimension of "//TRIM(NumberToVString(keyDimension,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(a,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(a,1)/=SIZE(values,1)) THEN
      localError="The size of the array data dimension of "//TRIM(NumberToVString(SIZE(a,1),"*",err,error))// &
        & " does not match the size of the values array of "//TRIM(NumberToVString(SIZE(values,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
   
    found=.FALSE.
    i=1
    DO WHILE(i<=SIZE(a,2).AND..NOT.found)
      IF(ABS(a(keyDimension,i)-values(keyDimension))<ZERO_TOLERANCE_DP) THEN
        found=.TRUE.
      ELSE
        i=i+1
      ENDIF
    ENDDO
    IF(found) THEN
      position=i
    ELSE
      position=0
    ENDIF
    
    EXITS("List_SearchLinearDP2Array")
    RETURN
999 ERRORSEXITS("List_SearchLinearDP2Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SearchLinearDP2Array
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the sorting method for a list
  SUBROUTINE List_SortingMethodSet(list,sortingMethod,err,error,*)

    !Argument Variables
    TYPE(ListType), POINTER, INTENT(IN) :: list !<A pointer to the list to set the sorting method for
    INTEGER(INTG), INTENT(IN) :: sortingMethod !<The searching method to set. See \see Lists_SortingMethod,Lists
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("List_SortingMethodSet",err,error,*999)

    CALL List_AssertNotFinished(list,err,error,*999)

    SELECT CASE(sortingMethod)
    CASE(LIST_BUBBLE_SORT_METHOD)
      list%sortMethod=LIST_BUBBLE_SORT_METHOD
    CASE(LIST_SHELL_SORT_METHOD)
      list%sortMethod=LIST_SHELL_SORT_METHOD
    CASE(LIST_HEAP_SORT_METHOD)
      list%sortMethod=LIST_HEAP_SORT_METHOD
    CASE DEFAULT
      localError="The specified list sorting method of "//TRIM(NumberToVString(sortingMethod,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("List_SortingMethodSet")
    RETURN
999 ERRORSEXITS("List_SortingMethodSet",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortingMethodSet

  !
  !================================================================================================================================
  !

  !>Sorts a list of into ascending order.
  SUBROUTINE List_SortList(list,err,error,*)
  
    !Argument variables
    TYPE(ListType), POINTER, INTENT(INOUT) :: list !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_SortList",err,error,*999)

    CALL List_AssertIsFinished(list,err,error,*999)
    
    SELECT CASE(list%sortMethod)
    CASE(LIST_BUBBLE_SORT_METHOD)
      IF(list%dataDimension==1) THEN
        SELECT CASE(list%dataType)
        CASE(LIST_INTG_TYPE)
          CALL List_SortBubbleIntg1Array(list%listIntg(1:list%numberInList),err,error,*999)
        CASE(LIST_SP_TYPE)
          CALL List_SortBubbleSP1Array(list%listSP(1:list%numberInList),err,error,*999)              
        CASE(LIST_DP_TYPE)
          CALL List_SortBubbleDP1Array(list%listDP(1:list%numberInList),err,error,*999) 
        CASE DEFAULT
          localError="The list data type of "//TRIM(NumberToVString(list%dataType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        SELECT CASE(list%dataType)
        CASE(LIST_INTG_TYPE)
          CALL List_SortBubbleIntg2Array(list%listIntg2(:,1:list%numberInList),list%keyDimension,err,error,*999)
        CASE(LIST_SP_TYPE)
          CALL List_SortBubbleSP2Array(list%listSP2(:,1:list%numberInList),list%keyDimension,err,error,*999)              
        CASE(LIST_DP_TYPE)
          CALL List_SortBubbleDP2Array(list%listDP2(:,1:list%numberInList),list%keyDimension,err,error,*999)
        CASE DEFAULT
          localError="The list data type of "//TRIM(NumberToVString(list%dataType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDIF
    CASE(LIST_SHELL_SORT_METHOD)
      IF(list%dataDimension==1) THEN
        SELECT CASE(list%dataType)
        CASE(LIST_INTG_TYPE)
          CALL List_SortShellIntg1Array(list%listIntg(1:list%numberInList),err,error,*999)
        CASE(LIST_SP_TYPE)
          CALL List_SortShellSP1Array(list%listSP(1:list%numberInList),err,error,*999)              
        CASE(LIST_DP_TYPE)
          CALL List_SortShellDP1Array(list%listDP(1:list%numberInList),err,error,*999)
        CASE DEFAULT
          localError="The list data type of "//TRIM(NumberToVString(list%dataType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        SELECT CASE(list%dataType)
        CASE(LIST_INTG_TYPE)
          CALL List_SortShellIntg2Array(list%listIntg2(:,1:list%numberInList),list%keyDimension,err,error,*999)
        CASE(LIST_SP_TYPE)
          CALL List_SortShellSP2Array(list%listSP2(:,1:list%numberInList),list%keyDimension,err,error,*999)              
        CASE(LIST_DP_TYPE)
          CALL List_SortShellDP2Array(list%listDP2(:,1:list%numberInList),list%keyDimension,err,error,*999)
        CASE DEFAULT
          localError="The list data type of "//TRIM(NumberToVString(list%dataType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDIF
    CASE(LIST_HEAP_SORT_METHOD)
      IF(list%dataDimension==1) THEN
        SELECT CASE(list%dataType)
        CASE(LIST_INTG_TYPE)
          CALL List_SortHeapIntg1Array(list%listIntg(1:list%numberInList),err,error,*999)
        CASE(LIST_SP_TYPE)
          CALL List_SortHeapSP1Array(list%listSP(1:list%numberInList),err,error,*999)              
        CASE(LIST_DP_TYPE)
          CALL List_SortHeapDP1Array(list%listDP(1:list%numberInList),err,error,*999)                            
        CASE DEFAULT
          localError="The list data type of "//TRIM(NumberToVString(list%dataType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        SELECT CASE(list%dataType)
        CASE(LIST_INTG_TYPE)
          CALL List_SortHeapIntg2Array(list%listIntg2(:,1:list%numberInList),list%keyDimension,err,error,*999)
        CASE(LIST_SP_TYPE)
          CALL List_SortHeapSP2Array(list%listSP2(:,1:list%numberInList),list%keyDimension,err,error,*999)              
        CASE(LIST_DP_TYPE)
          CALL List_SortHeapDP2Array(list%listDP2(:,1:list%numberInList),list%keyDimension,err,error,*999) 
        CASE DEFAULT
          localError="The list data type of "//TRIM(NumberToVString(list%dataType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDIF
    CASE DEFAULT
      localError="The list sort method of "//TRIM(NumberToVString(list%sortMethod,"*",err,error))//" is invlaid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("List_SortList")
    RETURN
999 ERRORSEXITS("List_SortList",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortList
  
  !
  !================================================================================================================================
  !

  !>Sorts an integer array list of data dimension 1 into ascending order.
  SUBROUTINE List_SortIntg1Array(a,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: a(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    
    ENTERS("List_SortIntg1Array",err,error,*999)

    !Default sort method is a heap sort
    CALL List_SortHeap(a,err,error,*999)    

    EXITS("List_SortIntg1Array")
    RETURN
999 ERRORSEXITS("List_SortIntg1Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortIntg1Array
  
  !
  !================================================================================================================================
  !

  !>Sorts an integer array list of data dimension > 1 into ascending order.
  SUBROUTINE List_SortIntg2Array(a,keyDimension,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: a(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: keyDimension !<The key dimension to sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    
    ENTERS("List_SortIntg2Array",err,error,*999)

    !Default sort method is a heap sort
    CALL List_SortHeap(a,keyDimension,err,error,*999)    

    EXITS("List_SortIntg2Array")
    RETURN
999 ERRORSEXITS("List_SortIntg2Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortIntg2Array
  
  !
  !================================================================================================================================
  !

  !>Sorts an single precision array list of data dimension 1 into ascending order.
  SUBROUTINE List_SortSP1Array(a,err,error,*)
      
    !Argument variables
    REAL(SP), INTENT(INOUT) :: a(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    
    ENTERS("List_SortSP1Array",err,error,*999)

    !Default sort method is a heap sort
    CALL List_SortHeap(a,err,error,*999)    

    EXITS("List_SortSP1Array")
    RETURN
999 ERRORSEXITS("List_SortSP1Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortSP1Array
  
  !
  !================================================================================================================================
  !

  !>Sorts an single precision array list of data dimension > 1 into ascending order.
  SUBROUTINE List_SortSP2Array(a,keyDimension,err,error,*)
      
    !Argument variables
    REAL(SP), INTENT(INOUT) :: a(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: keyDimension !<The key dimension to sort the list on.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    
    ENTERS("List_SortSP2Array",err,error,*999)

    !Default sort method is a heap sort
    CALL List_SortHeap(a,keyDimension,err,error,*999)    

    EXITS("List_SortSP2Array")
    RETURN
999 ERRORSEXITS("List_SortSP2Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortSP2Array
  
  !
  !================================================================================================================================
  !

  !>Sorts an double precision array list of data dimension 1 into ascending order.
  SUBROUTINE List_SortDP1Array(a,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(INOUT) :: a(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
     
    ENTERS("List_SortDP1Array",err,error,*999)

    !Default sort method is a heap sort
    CALL List_SortHeap(a,err,error,*999)    

    EXITS("List_SortDP1Array")
    RETURN
999 ERRORSEXITS("List_SortDP1Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortDP1Array
  
  !
  !================================================================================================================================
  !

  !>Sorts an double precision array list of data dimension > 1 into ascending order.
  SUBROUTINE List_SortDP2Array(a,keyDimension,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(INOUT) :: a(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: keyDimension !<The key dimension to sort on.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
     
    ENTERS("List_SortDP2Array",err,error,*999)

    !Default sort method is a heap sort
    CALL List_SortHeap(a,keyDimension,err,error,*999)    

    EXITS("List_SortDP2Array")
    RETURN
999 ERRORSEXITS("List_SortDP2Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortDP2Array
  
  !
  !================================================================================================================================
  !

  !>Performs a bubble sort on an integer array of data dimension 1 list
  SUBROUTINE List_SortBubbleIntg1Array(a,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: a(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: flag,i,j,k,value
    
    ENTERS("List_SortBubbleIntg1Array",err,error,*999)

    IF(SIZE(a,1)>1) THEN
      flag=SIZE(a,1)
      DO i=1,SIZE(a,1)
        k=flag-1
        flag=0
        DO j=1,k
          IF(a(j)>a(j+1)) THEN
            value=a(j)
            a(j)=a(j+1)
            a(j+1)=value
            flag=j
          ENDIF
        ENDDO !j
        IF(flag==0) EXIT
      ENDDO !i
    ENDIF

    EXITS("List_SortBubbleIntg1Array")
    RETURN
999 ERRORSEXITS("List_SortBubbleIntg1Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortBubbleIntg1Array
  
  !
  !================================================================================================================================
  !

  !>Performs a bubble sort on an integer array of data dimension > 1 list
  SUBROUTINE List_SortBubbleIntg2Array(a,keyDimension,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: a(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: keyDimension !<The key dimension of A to do the sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: flag,i,j,k,value(SIZE(a,1))
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_SortBubbleIntg2Array",err,error,*999)

    IF(keyDimension<=0.OR.keyDimension>SIZE(a,1)) THEN
      localError="The specified key dimension of "//TRIM(NumberToVString(keyDimension,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(A,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    IF(SIZE(a,2)>1) THEN
      flag=SIZE(a,2)
      DO i=1,SIZE(a,2)
        k=flag-1
        flag=0
        DO j=1,k
          IF(a(keyDimension,j)>a(keyDimension,j+1)) THEN
            value=a(:,j)
            a(:,j)=a(:,j+1)
            a(:,j+1)=value
            flag=j
          ENDIF
        ENDDO !j
        IF(flag==0) EXIT
      ENDDO !i
    ENDIF
    
    EXITS("List_SortBubbleIntg2Array")
    RETURN
999 ERRORSEXITS("List_SortBubbleIntg2Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortBubbleIntg2Array
  
  !
  !================================================================================================================================
  !

  !>BUBBLE_SORT_SP performs a bubble sort on a single precision array of data dimension 1 list
  SUBROUTINE List_SortBubbleSP1Array(A,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(INOUT) :: A(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: flag,i,j,k
    REAL(SP) :: value
    
    ENTERS("List_SortBubbleSP1Array",err,error,*999)

    IF(SIZE(A,1)>1) THEN
      flag=SIZE(A,1)
      DO i=1,SIZE(A,1)
        k=flag-1
        flag=0
        DO j=1,k
          IF(A(j)>A(j+1)) THEN
            value=A(j)
            A(j)=A(j+1)
            A(j+1)=value
            flag=j
          ENDIF
        ENDDO
        IF(flag==0) EXIT
      ENDDO
    ENDIF

    EXITS("List_SortBubbleSP1Array")
    RETURN
999 ERRORSEXITS("List_SortBubbleSP1Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortBubbleSP1Array
  
   !
  !================================================================================================================================
  !

  !>Performs a bubble sort on a single precision array of data dimension > 1 list
  SUBROUTINE List_SortBubbleSP2Array(a,keyDimension,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(INOUT) :: a(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: keyDimension !<The key dimension of A to do the sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: flag,i,j,k
    REAL(SP) :: value(SIZE(a,1))
    TYPE(VARYING_STRING) :: localError
     
    ENTERS("List_SortBubbleSP2Array",err,error,*999)

    IF(keyDimension<=0.OR.keyDimension>SIZE(a,1)) THEN
       localError="The specified key dimension of "//TRIM(NumberToVString(keyDimension,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(a,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    IF(SIZE(a,2)>1) THEN
      flag=SIZE(a,2)
      DO i=1,SIZE(a,2)
        k=flag-1
        flag=0
        DO j=1,k
          IF(a(keyDimension,j)>a(keyDimension,j+1)) THEN
            value=a(:,j)
            a(:,j)=a(:,j+1)
            a(:,j+1)=value
            flag=j
          ENDIF
        ENDDO !j
        IF(flag==0) EXIT
      ENDDO !i
    ENDIF

    EXITS("List_SortBubbleSP2Array")
    RETURN
999 ERRORSEXITS("List_SortBubbleSP2Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortBubbleSP2Array
  
  !
  !================================================================================================================================
  !

  !>Performs a bubble sort on a double precision of data dimension 1 list
  SUBROUTINE List_SortBubbleDP1Array(a,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(INOUT) :: a(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: flag,i,j,k
    REAL(DP) :: value
    
    ENTERS("List_SortBubbleDP1Array",err,error,*999)

    IF(SIZE(a,1)>1) THEN
      flag=SIZE(a,1)
      DO i=1,SIZE(a,1)
        k=flag-1
        flag=0
        DO j=1,k
          IF(a(j)>a(j+1)) THEN
            value=a(j)
            a(j)=a(j+1)
            a(j+1)=value
            flag=j
          ENDIF
        ENDDO !j
        IF(flag==0) EXIT
      ENDDO !i
    ENDIF

    EXITS("List_SortBubbleDP1Array")
    RETURN
999 ERRORSEXITS("List_SortBubbleDP1Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortBubbleDP1Array
  
  !
  !================================================================================================================================
  !

  !>Performs a bubble sort on a double precision of data dimension > 1 list
  SUBROUTINE List_SortBubbleDP2Array(a,keyDimension,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(INOUT) :: a(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: keyDimension !<The key dimension of A to do the sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: flag,i,j,k
    REAL(DP) :: value(SIZE(A,1))
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_SortBubbleDP2Array",err,error,*999)

    IF(keyDimension<=0.OR.keyDimension>SIZE(a,1)) THEN
      localError="The specified key dimension of "//TRIM(NumberToVString(keyDimension,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(a,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    IF(SIZE(a,2)>1) THEN
      flag=SIZE(a,2)
      DO i=1,SIZE(a,2)
        k=flag-1
        flag=0
        DO j=1,k
          IF(a(keyDimension,j)>a(keyDimension,j+1)) THEN
            value=a(:,j)
            a(:,j)=a(:,j+1)
            a(:,j+1)=value
            flag=j
          ENDIF
        ENDDO !j
        IF(flag==0) EXIT
      ENDDO !i
    ENDIF

    EXITS("List_SortBubbleDP2Array")
    RETURN
999 ERRORSEXITS("List_SortBubbleDP2Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortBubbleDP2Array
  
  !
  !================================================================================================================================
  !

  !>Sorts an integer array of data dimension 1 list into assending order using the heap sort method.
  SUBROUTINE List_SortHeapIntg1Array(a,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: a(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,iValue,j,l,value
    
    ENTERS("List_SortHeapIntg1Array",err,error,*999)

    IF(SIZE(a,1)>1) THEN      
      l=SIZE(a,1)/2+1
      iValue=SIZE(a,1)
      DO 
        IF(l>1) THEN
          l=l-1
          value=a(l)
        ELSE
          value=a(iValue)
          a(iValue)=a(1)
          iValue=iValue-1
          IF(iValue==1) THEN
            a(1)=value
            EXIT
          ENDIF
        ENDIF
        i=l
        j=l+l
        DO WHILE(j<=iValue)
          IF(j<iValue) THEN
            IF(a(j)<a(j+1)) j=j+1
          ENDIF
          IF(value<a(j)) THEN
            a(i)=a(j)
            i=j
            j=j+j
          ELSE
            j=iValue+1
          ENDIF
        ENDDO
        a(i)=value
      ENDDO
    ENDIF

    EXITS("List_SortHeapIntg1Array")
    RETURN
999 ERRORSEXITS("List_SortHeapIntg1Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortHeapIntg1Array
  
  !
  !================================================================================================================================
  !

  !>Sorts an integer array of data dimension > 1 list into assending order using the heap sort method.
  SUBROUTINE List_SortHeapIntg2Array(a,keyDimension,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: a(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: keyDimension !<The key dimension of A to do the sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,iValue,j,l,value(SIZE(a,1))
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_SortHeapIntg2Array",err,error,*999)
    
    IF(keyDimension<=0.OR.keyDimension>SIZE(a,1)) THEN
      localError="The specified key dimension of "//TRIM(NumberToVString(keyDimension,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(a,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    IF(SIZE(a,2)>1) THEN      
      l=SIZE(a,2)/2+1
      iValue=SIZE(a,2)
      DO 
        IF(l>1) THEN
          l=l-1
          value=a(:,l)
        ELSE
          value=a(:,iValue)
          a(:,iValue)=a(:,1)
          iValue=iValue-1
          IF(iValue==1) THEN
            a(:,1)=value
            EXIT
          ENDIF
        ENDIF
        i=l
        j=l+l
        DO WHILE(j<=iValue)
          IF(j<iValue) THEN
            IF(a(keyDimension,j)<a(keyDimension,j+1)) j=j+1
          ENDIF
          IF(value(keyDimension)<a(keyDimension,j)) THEN
            a(:,i)=a(:,j)
            i=j
            j=j+j
          ELSE
            j=iValue+1
          ENDIF
        ENDDO
        a(:,i)=value
      ENDDO
    ENDIF

    EXITS("List_SortHeapIntg2Array")
    RETURN
999 ERRORSEXITS("List_SortHeapIntg2Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortHeapIntg2Array
  
  !
  !================================================================================================================================
  !

  !>Sorts a real single precision array of data dimension 1 list into assending order using the heap sort method.
  SUBROUTINE List_SortHeapSP1Array(a,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(INOUT) :: a(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,iValue,j,l
    REAL(SP) :: value
    
    ENTERS("List_SortHeapSP1Array",err,error,*999)

    IF(SIZE(a,1)>1) THEN      
      l=SIZE(a,1)/2+1
      iValue=SIZE(a,1)
      DO 
        IF(l>1) THEN
          l=l-1
          value=a(l)
        ELSE
          value=a(iValue)
          a(iValue)=a(1)
          iValue=iValue-1
          IF(iValue==1) THEN
            a(1)=value
            EXIT
          ENDIF
        ENDIF
        i=l
        j=l+l
        DO WHILE(j<=iValue)
          IF(j<iValue) THEN
            IF(a(j)<a(j+1)) j=j+1
          ENDIF
          IF(value<a(j)) THEN
            a(i)=a(j)
            i=j
            j=j+j
          ELSE
            j=iValue+1
          ENDIF
        ENDDO
        a(i)=value
      ENDDO
    ENDIF

    EXITS("List_SortHeapSP1Array")
    RETURN
999 ERRORSEXITS("List_SortHeapSP1Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortHeapSP1Array
  
  !
  !================================================================================================================================
  !

  !>Sorts a real single precision array of data dimension > 1 list into assending order using the heap sort method.
  SUBROUTINE List_SortHeapSP2Array(a,keyDimension,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(INOUT) :: a(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: keyDimension !<The key dimension of A to do the sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,iValue,j,l
    REAL(SP) :: value(SIZE(a,1))
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_SortHeapSP2Array",err,error,*999)

    IF(keyDimension<=0.OR.keyDimension<=SIZE(A,1)) THEN
      localError="The specified key dimension of "//TRIM(NumberToVString(keyDimension,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(a,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    IF(SIZE(a,2)>1) THEN      
      l=SIZE(a,2)/2+1
      iValue=SIZE(a,2)
      DO 
        IF(l>1) THEN
          l=l-1
          value=a(:,l)
        ELSE
          value=a(:,iValue)
          a(:,iValue)=a(:,1)
          iValue=iValue-1
          IF(iValue==1) THEN
            a(:,1)=value
            EXIT
          ENDIF
        ENDIF
        i=l
        j=l+l
        DO WHILE(j<=iValue)
          IF(j<iValue) THEN
            IF(a(keyDimension,j)<a(keyDimension,j+1)) j=j+1
          ENDIF
          IF(value(keyDimension)<a(keyDimension,j)) THEN
            a(:,i)=a(:,j)
            i=j
            j=j+j
          ELSE
            j=iValue+1
          ENDIF
        ENDDO
        a(:,i)=value
      ENDDO
    ENDIF

    EXITS("List_SortHeapSP2Array")
    RETURN
999 ERRORSEXITS("List_SortHeapSP2Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortHeapSP2Array
  
  !
  !================================================================================================================================
  !
  
  !>Sorts a real double precision array of data dimension 1 list into assending order using the heap sort method.
  SUBROUTINE List_SortHeapDP1Array(a,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(INOUT) :: a(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,iValue,j,l
    REAL(DP) :: value
    
    ENTERS("List_SortHeapDP1Array",err,error,*999)

    IF(SIZE(a,1)>1) THEN      
      l=SIZE(a,1)/2+1
      iValue=SIZE(a,1)
      DO 
        IF(l>1) THEN
          l=l-1
          value=a(l)
        ELSE
          value=a(iValue)
          a(iValue)=a(1)
          iValue=iValue-1
          IF(iValue==1) THEN
            a(1)=value
            EXIT
          ENDIF
        ENDIF
        i=l
        j=l+l
        DO WHILE(j<=iValue)
          IF(j<iValue) THEN
            IF(a(j)<a(j+1)) j=j+1
          ENDIF
          IF(value<a(j)) THEN
            a(i)=a(j)
            i=j
            j=j+j
          ELSE
            j=iValue+1
          ENDIF
        ENDDO
        a(i)=value
      ENDDO
    ENDIF

    EXITS("List_SortHeapDP1Array")
    RETURN
999 ERRORSEXITS("List_SortHeapDP1Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortHeapDP1Array
  
  !
  !================================================================================================================================
  !
  
  !>Sorts a real double precision array of data dimension > 1 list into assending order using the heap sort method.
  SUBROUTINE List_SortHeapDP2Array(a,keyDimension,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(INOUT) :: a(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: keyDimension !<The key dimension of A to do the sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,iValue,j,l
    REAL(DP) :: value(SIZE(a,1))
    TYPE(VARYING_STRING) :: localError
   
    ENTERS("List_SortHeapDP2Array",err,error,*999)

    IF(keyDimension<=0.OR.keyDimension>SIZE(a,1)) THEN
      localError="The specified key dimension of "//TRIM(NumberToVString(keyDimension,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(a,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    IF(SIZE(a,2)>1) THEN      
      l=SIZE(a,2)/2+1
      iValue=SIZE(a,2)
      DO 
        IF(l>1) THEN
          l=l-1
          value=a(:,l)
        ELSE
          value=a(:,iValue)
          a(:,iValue)=a(:,1)
          iValue=iValue-1
          IF(iValue==1) THEN
            a(:,1)=value
            EXIT
          ENDIF
        ENDIF
        i=l
        j=l+l
        DO WHILE(j<=iValue)
          IF(j<iValue) THEN
            IF(a(keyDimension,j)<a(keyDimension,j+1)) j=j+1
          ENDIF
          IF(value(keyDimension)<a(keyDimension,j)) THEN
            a(:,i)=a(:,j)
            i=j
            j=j+j
          ELSE
            j=iValue+1
          ENDIF
        ENDDO
        a(:,i)=value
      ENDDO
    ENDIF

    EXITS("List_SortHeapDP2Array")
    RETURN
999 ERRORSEXITS("List_SortHeapDP2Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortHeapDP2Array
  
  !
  !================================================================================================================================
  !

  !>Sorts an integer array of data dimension 1 list into either assending or descending order using the shell sort method.
  SUBROUTINE List_SortShellIntg1Array(a,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: a(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,inc,j,value
    
    ENTERS("List_SortShellIntg1Array",err,error,*999)

    inc=4
    DO WHILE(inc<=SIZE(a,1))
      inc=3*inc+1
    ENDDO
    DO WHILE(inc>1)
      inc=inc/3
      DO i=inc+1,SIZE(a,1)
        value=A(i)
        j=i
        DO WHILE(a(j-inc)>value)
          a(j)=a(j-inc)
          j=j-inc
          IF(j<=inc) EXIT
        ENDDO
        a(j)=value
      ENDDO !i
    ENDDO

    EXITS("List_SortShellIntg1Array")
    RETURN
999 ERRORSEXITS("List_SortShellIntg1Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortShellIntg1Array
  
  !
  !================================================================================================================================
  !

  !>Sorts an integer array of data dimension > 1 list into either assending or descending order using the shell sort method.
  SUBROUTINE List_SortShellIntg2Array(a,keyDimension,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: a(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: keyDimension !<The key dimension of A to do the sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,inc,j,value(SIZE(a,1))
    TYPE(VARYING_STRING) :: localError
   
    ENTERS("List_SortShellIntg2Array",err,error,*999)

    IF(keyDimension<=0.OR.keyDimension>SIZE(a,1)) THEN
      localError="The specified key dimension of "//TRIM(NumberToVString(keyDimension,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(a,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    inc=4
    DO WHILE(inc<=SIZE(a,2))
      inc=3*inc+1
    ENDDO
    DO WHILE(inc>1)
      inc=inc/3
      DO i=inc+1,SIZE(a,2)
        value=a(:,i)
        j=i
        DO WHILE(a(keyDimension,j-inc)>value(keyDimension))
          a(:,j)=a(:,j-inc)
          j=j-inc
          IF(j<=inc) EXIT
        ENDDO
        a(:,j)=value
      ENDDO !i
    ENDDO

    EXITS("List_SortShellIntg2Array")
    RETURN
999 ERRORSEXITS("List_SortShellIntg2Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortShellIntg2Array
  
  !
  !================================================================================================================================
  !

  !>Sorts a real single precision array of data dimension 1 list into either assending or descending order using the shell
  !>sort method.
  SUBROUTINE List_SortShellSP1Array(a,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(INOUT) :: a(:) !<The list to sort
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,inc,j
    REAL(SP) :: value
    
    ENTERS("List_SortShellSP1Array",err,error,*999)

    inc=4
    DO WHILE(inc<=SIZE(a,1))
      inc=3*inc+1
    ENDDO
    DO WHILE(inc>1)
      inc=inc/3
      DO i=inc+1,SIZE(a,1)
        value=a(i)
        j=i
        DO WHILE(a(j-inc)>value)
          a(j)=a(j-inc)
          j=j-inc
          IF(j<=inc) EXIT
        ENDDO
        a(j)=value
      ENDDO !i
    ENDDO

    EXITS("List_SortShellSP1Array")
    RETURN
999 ERRORSEXITS("List_SortShellSP1Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortShellSP1Array
  
  !
  !================================================================================================================================
  !

  !>Sorts a real single precision array of data dimension > 1 list into either assending or descending order using the shell
  !>sort method.
  SUBROUTINE List_SortShellSP2Array(a,keyDimension,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(INOUT) :: A(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: keyDimension !<The key dimension of A to do the sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,inc,j
    REAL(SP) :: value(SIZE(a,1))
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_SortShellSP2Array",err,error,*999)

    IF(keyDimension<=0.OR.keyDimension>SIZE(a,1)) THEN
      localError="The specified key dimension of "//TRIM(NumberToVString(keyDimension,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(a,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    inc=4
    DO WHILE(inc<=SIZE(a,2))
      inc=3*inc+1
    ENDDO
    DO WHILE(inc>1)
      inc=inc/3
      DO i=inc+1,SIZE(a,2)
        value=a(:,i)
        j=i
        DO WHILE(a(keyDimension,j-inc)>value(keyDimension))
          a(:,j)=a(:,j-inc)
          j=j-inc
          IF(j<=inc) EXIT
        ENDDO
        a(:,j)=VALUE
      ENDDO !i
    ENDDO

    EXITS("List_SortShellSP2Array")
    RETURN
999 ERRORSEXITS("List_SortShellSP2Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortShellSP2Array
  
  !
  !================================================================================================================================
  !

  !>Sorts a real double precision array of data dimension 1 list into either assending or descending order using the shell
  !>sort method.
  SUBROUTINE List_SortShellDP1Array(a,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(INOUT) :: a(:)
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,inc,j
    REAL(DP) :: value
    
    ENTERS("List_SortShellDP1Array",err,error,*999)

    inc=4
    DO WHILE(inc<=SIZE(a,1))
      inc=3*inc+1
    ENDDO
    DO WHILE(inc>1)
      inc=inc/3
      DO i=inc+1,SIZE(a,1)
        value=a(i)
        j=i
        DO WHILE(a(j-inc)>value)
          a(j)=a(j-inc)
          j=j-inc
          IF(j<=inc) EXIT
        ENDDO
        a(j)=value
      ENDDO !i
    ENDDO

    EXITS("List_SortShellDP1Array")
    RETURN
999 ERRORSEXITS("List_SortShellDP1Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortShellDP1Array
  
  !
  !================================================================================================================================
  !

  !>Sorts a real double precision array of data dimension 2 list into either assending or descending order using the shell
  !>sort method.
  SUBROUTINE List_SortShellDP2Array(a,keyDimension,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(INOUT) :: A(:,:) !<The list to sort
    INTEGER(INTG), INTENT(IN) :: keyDimension !<The key dimension of A to do the sort on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,inc,j
    REAL(DP) :: value(SIZE(a,1))
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("List_SortShellDP2Array",err,error,*999)

    IF(keyDimension<=0.OR.keyDimension>SIZE(a,1)) THEN
      localError="The specified key dimension of "//TRIM(NumberToVString(keyDimension,"*",err,error))// &
        & " is invalid. The key dimension must be > 0 and <= "//TRIM(NumberToVString(SIZE(a,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    inc=4
    DO WHILE(inc<=SIZE(a,2))
      inc=3*inc+1
    ENDDO
    DO WHILE(inc>1)
      inc=inc/3
      DO i=inc+1,SIZE(a,2)
        value=a(:,i)
        j=i
        DO WHILE(a(keyDimension,j-inc)>value(keyDimension))
          a(:,j)=a(:,j-inc)
          j=j-inc
          IF(j<=inc) EXIT
        ENDDO
        a(:,j)=value
      ENDDO !i
    ENDDO

    EXITS("List_SortShellDP2Array")
    RETURN
999 ERRORSEXITS("List_SortShellDP2Array",err,error)
    RETURN 1
    
  END SUBROUTINE List_SortShellDP2Array
  
  !
  !================================================================================================================================
  !

  !>Finds the intersection of two sets (arrays), leaving the original arrays intact
  SUBROUTINE List_IntersectionIntgArray(a,b,c,err,error,*)
    
    ! Argument variables
    INTEGER(INTG), INTENT(IN), TARGET :: a(:)   !<One of the two arrays to find the intersection for
    INTEGER(INTG), INTENT(IN), TARGET :: b(:)   !<Other array to find the intersection for
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: c(:) !<On exit, contains the list of common elements of the arrays. Must not be allocated on entry.
    INTEGER(INTG), INTENT(OUT) :: err          !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    ! Local variables
    INTEGER(INTG) :: sizeShorter,sizeLonger
    INTEGER(INTG) :: i,j,start,numberOfMatches
    INTEGER(INTG), POINTER :: longer(:),shorter(:)
    INTEGER(INTG), ALLOCATABLE :: matches(:)
    INTEGER(INTG), ALLOCATABLE :: longArray(:),shortArray(:)   !<copies, if needed
    
    ENTERS("List_IntersectionIntgArray",err,error,*999)

    !If the lists are small, it's probably easier to directly compare: O(n^2)
    !but if they're big, sort first then compare: O(n log n)*2 + 2*O(n)

    IF(ALLOCATED(c)) CALL FlagError("Output C array is already allocated.",err,error,*999)
    
    !Start finding the intersection
    NULLIFY(longer)
    NULLIFY(shorter)
    !It's quicker to compare shorter array elements to longer ones
    IF(SIZE(a)>SIZE(b)) THEN
      longer=>a
      shorter=>b
    ELSE
      longer=>b
      shorter=>a
    ENDIF
    sizeShorter=SIZE(shorter)
    sizeLonger=SIZE(longer)
    ALLOCATE(matches(sizeShorter),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate matches array.",err,error,*999)
    numberOfMatches=0
    
    !Long or short lists?
    IF(sizeLonger*sizeShorter<=1E4) THEN  ! a rather arbitrary cutoff...
      !'short' lists - begin comparing straight away
      DO i=1,sizeShorter
        DO j=1,sizeLonger
          IF(shorter(i)==longer(j)) THEN
            numberOfMatches=numberOfMatches+1
            matches(numberOfMatches)=shorter(i)
          ENDIF
        ENDDO !j
      ENDDO !i
    ELSE
      !'long' lists - make copies of the arrays
      ALLOCATE(longArray(sizeLonger),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate long array.",err,error,*999)
      ALLOCATE(shortArray(sizeShorter),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate short array.",err,error,*999)
      longArray(1:sizeLonger)=longer(1:sizeLonger)
      shortArray(1:sizeShorter)=shorter(1:sizeShorter)
      !Sort both arrays
      CALL List_Sort(longArray,err,error,*999)
      CALL List_Sort(shortArray,err,error,*999)
      ! compare now
      start=1
      DO i=1,sizeShorter
        DO j=start,sizeLonger
          IF(longArray(j)==shortArray(i)) THEN
            numberOfMatches=numberOfMatches+1
            matches(numberOfMatches)=shortArray(i)
            start=MIN(j+1,sizeLonger)
            EXIT
          ELSE IF(longArray(j)>shortArray(i)) THEN
            !Can start here next time
            start=MAX(j-1,1)
            EXIT
          ENDIF
        ENDDO
      ENDDO
      DEALLOCATE(longArray)
      DEALLOCATE(shortArray)
    ENDIF !Long or short lists
    !Cut the array down to size
    ALLOCATE(c(numberOfMatches),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate c array.",err,error,*999)
    c(1:numberOfMatches)=matches(1:numberOfMatches)
    DEALLOCATE(matches)

    EXITS("List_IntersectionIntgArray")
    RETURN
999 IF(ALLOCATED(matches)) DEALLOCATE(matches)
    IF(ALLOCATED(longArray)) DEALLOCATE(longArray)
    IF(ALLOCATED(shortArray)) DEALLOCATE(shortArray)
     ERRORSEXITS("List_IntersectionIntgArray",err,error)
    RETURN 1

  END SUBROUTINE List_IntersectionIntgArray

  !
  !================================================================================================================================
  !

  !>Finds out whether array A is a subset of array B
  SUBROUTINE List_SubsetOfIntgArray(a,b,subset,err,error,*)
    ! Argument variables
    INTEGER(INTG), INTENT(IN) :: a(:) !<Supposed subset (to test for)
    INTEGER(INTG), INTENT(IN) :: b(:) !<Supposed superset
    LOGICAL, INTENT(OUT) :: subset !<On exit, .TRUE. if a is a subset of b, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    ! Logical variables
    INTEGER(INTG) :: sizeA,sizeB,i,j,start,sizeReduce
    INTEGER(INTG), ALLOCATABLE :: aSorted(:),bSorted(:)

    ENTERS("List_SubsetOfIntgArray",err,error,*999)

    sizeA=SIZE(a)
    sizeB=SIZE(b)
    subset=.FALSE.
    
    IF(sizeA<=sizeB) THEN
      sizeReduce=0
      DO i=1,sizeA
        IF(a(i)==0) sizeReduce=sizeReduce+1
      ENDDO !i
      sizeA=sizeA-sizeReduce
      sizeReduce=0
      DO i=1,sizeB
        IF(b(i)==0) sizeReduce=sizeReduce+1
      ENDDO !i
      sizeB=sizeB-sizeReduce

      !Short of long arrays?
      IF(sizeA*sizeB<=1E4) THEN
        !'short' arrays - just compare without sorting
        iLoop: DO i=1,sizeA
          jLoop: DO j=1,sizeB
            IF(a(i)==b(j)) THEN
              EXIT jLoop
            ELSE IF(j==sizeB) THEN
              EXIT iLoop
            ENDIF
          ENDDO jLoop
          IF(I==sizeA) subset=.TRUE.
        ENDDO iLoop
      ELSE
        !'long' arrays - sort first
        ALLOCATE(aSorted(sizeA),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate a sorted array.",err,error,*999)
        ALLOCATE(bSorted(sizeB),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate b sorted array.",err,error,*999)
        aSorted(1:sizeA)=a(1:sizeA)
        bSorted(1:sizeB)=b(1:sizeB)
        CALL List_Sort(aSorted,err,error,*999)
        CALL List_Sort(bSorted,err,error,*999)
        start=1
        iLoop2: DO i=1,sizeA
          jLoop2: DO j=1,sizeB
            IF(a(i)==b(j)) THEN
              start=MIN(j+1,sizeB)
              EXIT jLoop2
            ELSE IF(a(i)<b(j)) THEN
              EXIT iLoop2
            ENDIF
          ENDDO jLoop2
          IF(i==sizeA) subset=.TRUE.
        ENDDO iLoop2
        DEALLOCATE(aSorted)
        DEALLOCATE(bSorted)
      ENDIF
    ENDIF

    EXITS("List_SubsetOfIntgArray")
    RETURN
999 IF(ALLOCATED(aSorted)) DEALLOCATE(aSorted)
    IF(ALLOCATED(bSorted)) DEALLOCATE(bSorted)
    ERRORSEXITS("List_SubsetOfIntgArray",err,error)
    RETURN 1

  END SUBROUTINE List_SubsetOfIntgArray

  !
  !================================================================================================================================
  !

END MODULE Lists
