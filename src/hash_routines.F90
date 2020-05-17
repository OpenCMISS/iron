! Implements hash tables for static sets (i.e. do not change) as described in
! Fredman, Komlos and Szemeredi 1984
! using (arrays of) Lists that allow for INTG or REAL (implemented SP) of different sizes as data content (values)

MODULE HashRoutines

  USE BaseRoutines
  USE INPUT_OUTPUT
  USE CONSTANTS
  USE KINDS
  USE ISO_VARYING_STRING
  USE STRINGS 
  USE ComputationEnvironment
  USE Types
  USE Lists

#include "macros.h"  

  IMPLICIT NONE
  
  PRIVATE  ! by default

  !Module parameters
  INTEGER(INTG), PARAMETER :: HASH_INITIAL_ADDITIONAL_VALUES_SIZE=5 !<Initial number of additional lists of values 

  ! TYPEs

  ! Items for Linked Lists with 2 entries
  TYPE HashLinkedListItemType
    PRIVATE
    INTEGER(INTG) :: key    
    INTEGER(INTG) :: value
    TYPE(HashLinkedListItemType),POINTER :: next => NULL()
  END TYPE

  ! Linked Lists
  TYPE HashLinkedListType
    PRIVATE
    TYPE(HashLinkedListItemType),POINTER :: root => NULL()
    TYPE(HashLinkedListItemType),POINTER :: last => NULL()
  END TYPE

  TYPE HashLinkedListPtrType
    PRIVATE
    TYPE(HashLinkedListType), POINTER :: ptr
  END TYPE HashLinkedListPtrType

  ! Array of Linked Lists
  TYPE HashListArrayType
     PRIVATE
     TYPE(HashLinkedListPtrType), ALLOCATABLE :: vec(:)
     INTEGER(INTG)                            :: vecLen
     LOGICAL(INTG)                            :: isFinished
  END TYPE HashListArrayType

  ! The Hash Table
  ! TYPE HashTableType
  ! in types.f90

  INTERFACE HashLinkedList_Add
    MODULE PROCEDURE HashLinkedList_AddData
    MODULE PROCEDURE HashLinkedList_AddList
  END INTERFACE HashLinkedList_Add

  INTERFACE HashTable_ListValuesSet
    MODULE PROCEDURE HashTable_ListRealSpValuesSet
    MODULE PROCEDURE HashTable_ListIntValuesSet
  END INTERFACE HashTable_ListValuesSet

  INTERFACE HashTable_ValuesSetAndInsert
    MODULE PROCEDURE HashTable_IntValuesSetAndInsert
    MODULE PROCEDURE HashTable_RealSpValuesSetAndInsert
  END INTERFACE HashTable_ValuesSetAndInsert

  INTERFACE HashTable_GetValue
    MODULE PROCEDURE HashTable_GetRealSpValue
    MODULE PROCEDURE HashTable_GetIntValue
    MODULE PROCEDURE HashTable_GetListArrayIntValue
    MODULE PROCEDURE HashTable_GetListArrayRealSpValue
  END INTERFACE HashTable_GetValue

  INTERFACE HashTable_AdditionalValuesSet
    MODULE PROCEDURE HashTable_AdditionalIntValuesSet
    MODULE PROCEDURE HashTable_AdditionalRealSpValuesSet
  END INTERFACE HashTable_AdditionalValuesSet

  ! public subs
  ! Array of lists (keep private)
  !PUBLIC :: HashListArray_CreateStart, HashListArray_CreateFinish
  ! Table create
  PUBLIC :: HashTable_ValuesSetAndInsert, HashTable_CreateStart, HashTable_CreateFinish, HashTable_AdditionalValuesSet
  ! Table get values
  PUBLIC :: HashTable_GetKey, HashTable_GetValue

CONTAINS

!==================================================

  !> initialises or adds a piece of data to list
  SUBROUTINE HashLinkedList_AddData(list,key, value, isDuplicateKey, err,error,*)

    !Argument Variables
    TYPE(HashLinkedListType), POINTER :: list

    INTEGER(INTG),INTENT(IN) :: key, value

    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    LOGICAL, INTENT (OUT) :: isDuplicateKey

    !Local Variables
    TYPE(HashLinkedListItemType),POINTER :: current => NULL()

    ENTERS("HashLinkedList_AddData",err,error,*999)

    isDuplicateKey = .FALSE.

    ! Allocate the list if not done already
    IF(.NOT. ASSOCIATED(list)) THEN
        ALLOCATE(list, STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate list",err,error,*999)
    END IF    

    IF(ASSOCIATED(list%root)) THEN
      current => list%root
       DO
         ! key is there already
         IF(current%key==key) THEN
         isDuplicateKey = .TRUE.
         EXIT ! do nothing
         END IF 

         ! tail is reached
         IF(.NOT. ASSOCIATED(current%next)) THEN
            ALLOCATE(current%next)
            current%next%key   = key
            current%next%value = value
            list%last => current%next
         EXIT
         END IF
          
         ! proceed
         current => current%next 

        END DO  
       
     ELSE     ! list is empty: create root 

      ALLOCATE(list%root)
      list%root%key   = key
      list%root%value = value
      list%last => list%root
    END IF

    EXITS("HashLinkedList_AddData")
    RETURN
999 ERRORSEXITS("HashLinkedList_AddData",err,error)
    RETURN 1

  END SUBROUTINE HashLinkedList_AddData

!==================================================

  !> Returns length of list
  SUBROUTINE HashLinkedList_Size(list,n,err,error,*)

    ! argument variables
    TYPE(HashLinkedListType), POINTER :: list
    INTEGER(INTG), INTENT(OUT) :: n
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
 
    ! local variables
    TYPE(HashLinkedListItemType),POINTER :: current => NULL()
    LOGICAL :: isEmpty

    ENTERS("HashLinkedList_Size", err, error, *999)

    ! Check if list is allocated
    ! Return zero-size array also if list is empty
    CALL HashLinkedList_IsEmpty(list, isEmpty, err, error, *999)
    IF((.NOT. ASSOCIATED(list)).OR. isEmpty) THEN
        n = 0
    ELSE

    ! Traversing to find size
     current => list%root
     n=1
     DO
       IF(ASSOCIATED(current%next)) THEN
         n=n+1
         current => current%next
       ELSE
         EXIT
       END IF
     END DO
    
    END IF


    EXITS("HashLinkedList_Size")
    RETURN
999 ERRORSEXITS("HashLinkedList_Size",err,error)
    RETURN 1

  END SUBROUTINE HashLinkedList_Size

!==================================================

  SUBROUTINE HashLinkedList_GetData(list,key, value, isFound, err, error, *)

    ! argument variables
    TYPE(HashLinkedListType), POINTER :: list
    INTEGER(INTG),INTENT(IN)   :: key
    INTEGER(INTG),INTENT(OUT)  :: value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.

    ! local variables
    TYPE(HashLinkedListItemType),POINTER :: current => NULL()
    LOGICAL :: isFound

    ENTERS("HashLinkedList_GetData", err, error, *999)
 
    isFound = .FALSE.
   
    IF(ASSOCIATED(list%root)) THEN
    
    current => list%root
    DO
         ! key is found 
         IF(current%key==key)THEN
          value = current%value 
           
          isFound = .TRUE.

         END IF 

         ! tail is reached
         IF(.NOT. ASSOCIATED(current%next)) EXIT
          
         ! otherwise proceed
         current => current%next 
    END DO
  
    END IF
   
    EXITS("HashLinkedList_GetData")
    RETURN
999 ERRORSEXITS("HashLinkedList_GetData",err,error)
    RETURN 1


   END SUBROUTINE HashLinkedList_GetData

!==================================================

  !> adds all data from one list (addList) to another (list)
  SUBROUTINE HashLinkedList_AddList(list,addList,err,error, *)

    ! argument variables
    TYPE(HashLinkedListType), POINTER :: list
    TYPE(HashLinkedListType), POINTER :: addList
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    ! local variables
    LOGICAL :: isFound, isEmpty
    TYPE(HashLinkedListItemType),POINTER :: current => NULL()

    ENTERS("HashLinkedList_AddList", err, error, *999)

    CALL HashLinkedList_IsEmpty(addList, isEmpty, err, error, *999)
    IF(.NOT. isEmpty) THEN

    current => addList%root
    DO
      CALL HashLinkedList_AddData(list,current%key, current%value, isFound, err,error,*999)
      IF(ASSOCIATED(current%next)) THEN
        current => current%next
      ELSE
        EXIT
      END IF
    END DO

    END IF

    EXITS("HashLinkedList_AddList")
    RETURN
999 ERRORSEXITS("HashLinkedList_AddList",err,error)
    RETURN 1

  END SUBROUTINE HashLinkedList_AddList

!==================================================

  !> removes the first item from list and returns its value in data
  SUBROUTINE HashLinkedList_RemoveFirst(list,key,value,err,error,*)
 
    ! Argument Variables
    TYPE(HashLinkedListType), POINTER:: list
    INTEGER(INTG),INTENT(OUT) :: key, value
    INTEGER(INTG), INTENT(OUT)        :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.

    !Local Variables
    TYPE(HashLinkedListItemType),POINTER :: next

    ENTERS("HashLinkedList_RemoveFirst", err, error, *999)
    
    IF(ASSOCIATED(list%root)) THEN
      key   = list%root%key
      value = list%root%value
      next => list%root%next
      DEALLOCATE(list%root)
      list%root => next
      IF(ASSOCIATED(list%root)) THEN
        IF(.NOT.ASSOCIATED(list%root%next)) list%last => list%root  ! only one left
      ELSE
        list%last => NULL()
      END IF
    ELSE
      IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE, &
        & ">>> warning: Linked list is empty and first item cannot be removed.",err,error,*999) 
      ! write(*,*) ">>> warning: Linked list is empty and first item cannot be removed."
      END IF
    END IF


    EXITS("HashLinkedList_RemoveFirst")
    RETURN
999 ERRORSEXITS("HashLinkedList_RemoveFirst",err,error)
    RETURN

  END SUBROUTINE HashLinkedList_RemoveFirst

!==================================================

  !> removes the last item from list and returns its value in data
  SUBROUTINE HashLinkedList_RemoveLast(list,key,value,err,error,*)
    
    ! argument variables
    TYPE(HashLinkedListType), POINTER :: list
    INTEGER(INTG),INTENT(OUT) :: key,value
    INTEGER(INTG), INTENT(OUT)        :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.

    ! local variables
    TYPE(HashLinkedListItemType),POINTER :: current => NULL()

    ENTERS("HashLinkedList_RemoveLast", err, error, *999)

    IF(.NOT.ASSOCIATED(list%root)) THEN
      IF(diagnostics1) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE, &
          & ">>> warning: Linked list is empty and last item cannot be removed.",err,error,*999) 
      END IF
      !write(*,*) ">>> warning: Linked list is empty and last item cannot be removed."
      GOTO 987
    END IF
    current => list%root

    DO
      IF(ASSOCIATED(current%next)) THEN
        IF(ASSOCIATED(current%next%next)) THEN
          current => current%next
        ELSE
          ! next one is the last one
          key = current%next%key
          value = current%next%value
          DEALLOCATE(current%next)
          current%next => NULL()
          list%last => current
          EXIT
        END IF
      ELSE
        ! there must be only one item in the list(?)!
        key   = current%key
        value = current%value
        DEALLOCATE(list%root)
        list%root => NULL()
        list%last => NULL()
        EXIT
      END IF
    END DO

987 EXITS("HashLinkedList_RemoveLast")
    RETURN
999 ERRORSEXITS("HashLinkedList_RemoveLast",err,error)
    RETURN

  END SUBROUTINE HashLinkedList_RemoveLast

!==================================================

  !> will delete and deallocate all items
  SUBROUTINE HashLinkedList_Destroy(list,err,error, *)

    ! argument variables
    TYPE(HashLinkedListType), POINTER :: list
    INTEGER(INTG), INTENT(OUT)        :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    ! local variables
    TYPE(HashLinkedListItemType), POINTER :: current =>NULL()
    TYPE(HashLinkedListItemType), POINTER :: next=>NULL()

    ENTERS("HashLinkedList_Destroy", err, error, *999)

    IF(.NOT.ASSOCIATED(list%root)) GOTO 987

    current => list%root
    DO
      IF(ASSOCIATED(current%next)) THEN
        next => current%next
        DEALLOCATE(current)
        current => next
      ELSE
        DEALLOCATE(current)
        EXIT
      END IF
    END DO
    
    list%root => NULL()
    list%last => NULL()
    
    ! deallocate memory
    DEALLOCATE(list)

987 EXITS("HashLinkedList_Destroy")
    RETURN
999 ERRORSEXITS("HashLinkedList_Destroy",err,error)
    RETURN

  END SUBROUTINE HashLinkedList_Destroy

!==================================================

  !> returns true if the list is empty
  SUBROUTINE HashLinkedList_IsEmpty(list, isEmpty, err, error, *)

    ! Argument Variables
    TYPE(HashLinkedListType), POINTER :: list
    LOGICAL, INTENT(OUT)              :: isEmpty
    INTEGER(INTG), INTENT(OUT)        :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.

    ENTERS("HashLinkedList_IsEmpty", err, error, *999)

    isEmpty = .TRUE.

    IF(ASSOCIATED(list)) THEN
      IF(ASSOCIATED(list%root)) isEmpty = .FALSE.
    END IF

    EXITS("HashLinkedList_IsEmpty")
    RETURN
999 ERRORSEXITS("HashLinkedList_IsEmpty",err,error)
    RETURN

  END SUBROUTINE HashLinkedList_IsEmpty

!==================================================

  !> copies out the data to an allocatable array of size 2xn
  SUBROUTINE HashLinkedList_ListToArray(list,array,err,error,*)

    ! argument variables
    TYPE(HashLinkedListType),POINTER:: list

    INTEGER(INTG),ALLOCATABLE,INTENT(OUT) :: array(:,:)
    INTEGER(INTG), INTENT(OUT)        :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.

    ! local variables
    INTEGER(INTG) :: i,n
    LOGICAL       :: isEmpty
    TYPE(HashLinkedListItemType),POINTER :: current => NULL()

    ENTERS("HashLinkedList_ListToArray", err, error, *999)

    ! return zero-size array if list is empty
    CALL HashLinkedList_IsEmpty(list, isEmpty, err, error, *999)
    IF(isEmpty) THEN
      ALLOCATE(array(0,0), STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate list array",err,error,*999)
      GOTO 987
    END IF

    ! first traversing to find size
    current => list%root
    n=1
    DO
      IF(ASSOCIATED(current%next)) THEN
        n=n+1
        current => current%next
      ELSE
        EXIT
      END IF
    END DO

    ! copy to array
    IF(ALLOCATEd(array)) DEALLOCATE(array)
    ALLOCATE(array(2,n),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate list array",err,error,*999)
    current => list%root
    DO i=1,n
      array(1,i)=current%key
      array(2,i)=current%value

      current => current%next
    END DO

987 EXITS("HashLinkedList_ListToArray")
    RETURN
999 ERRORSEXITS("HashLinkedList_ListToArray",err,error)
    RETURN

  END SUBROUTINE HashLinkedList_ListToArray

!==================================================

  !>Starts the creation of a list array and returns a pointer to the created list
  SUBROUTINE HashListArray_CreateStart(listArray,err,error,*)

    !Argument Variables 
    TYPE (HashListArrayType), POINTER :: listArray !<A pointer to the list array to create. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    INTEGER(INTG)        :: dummyErr !<The error code.
    TYPE(VARYING_STRING) :: dummyError !<The error string.

    ENTERS("HashListArray_CreateStart",err,error,*999)

    IF(ASSOCIATED(listArray)) THEN
      CALL FlagError("List Array is already associated",err,error,*998)
    ELSE
      ALLOCATE(listArray,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate list array",err,error,*999)
    ! CALL HashListArrayType_Initialise(tbl,err,error,*999)
      !Set Defaults
      listArray%vecLen=10 ! default number of elements in array
      listArray%isFinished = .FALSE.
    END IF

    EXITS("HashListArray_CreateStart")
    RETURN
999 CALL HashListArray_Finalise(listArray,dummyErr,dummyError,*998)
998 ERRORSEXITS("HashListArray_CreateStart",err,error)
    RETURN 1

  END SUBROUTINE HashListArray_CreateStart

!==================================================

  !>Finishes the creation of a List Array created with HashListArray_CreateStart 
  SUBROUTINE HashListArray_CreateFinish(listArray,err,error,*)

    !Argument Variables 
    TYPE (HashListArrayType), POINTER :: listArray !<A pointer to the list array to create. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    INTEGER(INTG) :: i

    ENTERS("HashListArray_CreateFinish",err,error,*999)

    IF(ASSOCIATED(listArray)) THEN
       IF(listArray%isFinished) THEN
         CALL FlagError("List Array is already finished",err,error,*999)
       ELSE  
         ! Allocate lists
         ALLOCATE(listArray%vec(0:listArray%vecLen-1), STAT=err)
         IF(err/=0) CALL FlagError("Could not allocate list array",err,error,*999)
         ! Nullify each list pointer
         DO i = 0, listArray%vecLen-1
            NULLIFY(listArray%vec(i)%ptr)
         END DO

         listArray%isFinished = .TRUE.
       END IF
    ELSE
       CALL FlagError("List Array is not associated",err,error,*999)
    END IF

    EXITS("HashListArray_CreateFinish")
    RETURN
999 ERRORSEXITS("HashListArray_CreateFinish",err,error)
    RETURN 1

  END SUBROUTINE HashListArray_CreateFinish

!==================================================

  !>Finalises a list array and deallocates all memory.
  SUBROUTINE HashListArray_Finalise(listArray,err,error,*)

    !Argument Variables
    TYPE (HashListArrayType), POINTER :: listArray !<The list array to finalise
    INTEGER(INTG), INTENT(OUT)    :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("HashListArray_Finalise",err,error,*999)
 
   IF(ASSOCIATED(listArray)) THEN
      IF(ALLOCATED(listArray%vec)) DEALLOCATE(listArray%vec)
      DEALLOCATE(listArray)
   END IF

    EXITS("HashListArray_Finalise")
    RETURN
999 ERRORSEXITS("HashListArray_Finalise",err,error)
    RETURN 1
  END SUBROUTINE HashListArray_Finalise


  !>Finishes the creation of a table created with HashTable_CreateStart
  SUBROUTINE HashTable_CreateFinish(table,err,error,*)

    !Argument Variables
    TYPE(HashTableType), POINTER :: table !<A pointer to the table to finish
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables

!    not needed    
!    TYPE(VARYING_STRING) :: dummyError !<The error string.
!    INTEGER(INTG)        :: dummyErr

    ENTERS("HashTable_CreateFinish",err,error,*998)

    IF(ASSOCIATED(table)) THEN
      IF(table%hashTableFinished) THEN
        CALL FlagError("Table is already finished",err,error,*998)
      ELSE
         ! Allocate pointers
         ! No pointers to allocate
         ! Finish the tree creation
        table%hashTableFinished=.TRUE.
      END IF
    ELSE
      CALL FlagError("Table is not associated",err,error,*998)
    END IF

    EXITS("HashTable_CreateFinish")
    RETURN
!999 CALL HashTable_Finalise(table,dummyErr,dummyError,*998)
! not needed
998 ERRORSEXITS("HashTable_CreateFinish",err,error)
    RETURN 1
  END SUBROUTINE HashTable_CreateFinish

!==================================================

  !>Starts the creation of a hash table and returns a pointer to the created table.
  SUBROUTINE HashTable_CreateStart(table,err,error,*)

    !Argument Variables
    TYPE(HashTableType), POINTER :: table !<A pointer to the table to create. Must NOT be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variable

    TYPE(VARYING_STRING) :: dummyError !<The error string.
    INTEGER(INTG)        :: dummyErr

    ENTERS("HashTable_CreateStart",err,error,*998)

    IF(ASSOCIATED(table)) THEN
      CALL FlagError("Table is already associated",err,error,*998)
    ELSE
      ALLOCATE(table,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate table",err,error,*999)
      CALL HashTable_Initialise(table,err,error,*999)
    END IF

    EXITS("HashTable_CreateStart")
    RETURN
999 CALL HashTable_Finalise(table,dummyErr,dummyError,*998)
998 ERRORSEXITS("HashTable_CreateStart",err,error)
    RETURN 1
  END SUBROUTINE HashTable_CreateStart

!==================================================

 !>Initialises a hash table
  SUBROUTINE HashTable_Initialise(table,err,error,*)

    !Argument Variables
    TYPE(HashTableType), POINTER :: table !<A pointer to the table to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local Variables
    
    ENTERS("HashTable_Initialise",err,error,*999)

    IF(ASSOCIATED(table)) THEN
      table%hashTableFinished=.FALSE.
      table%n=0
      table%p=0
      table%listSVal=>NULL()
      IF(ALLOCATED(table%vecTKey)) DEALLOCATE(table%vecTKey)
      IF(ALLOCATED(table%vecTVal)) DEALLOCATE(table%vecTVal)
      IF(ALLOCATED(table%vecSKey)) DEALLOCATE(table%vecSKey)
      IF(ALLOCATED(table%arrayOfListSVal)) DEALLOCATE(table%arrayOfListSVal)
    ELSE
      CALL FlagError("Table is not associated",err,error,*999)
    END IF

    EXITS("HashTable_Initialise")
    RETURN
999 ERRORSEXITS("HashTable_Initialise",err,error)
    RETURN 1
  END SUBROUTINE HashTable_Initialise

!==================================================

! Queries for a value q and returns index
  SUBROUTINE HashTable_GetKey(table, q, indexFound, isFound, err, error, *) 

    !Argument variables
    TYPE(HashTableType), POINTER :: table
    INTEGER(INTG), INTENT(IN)    :: q          ! queried value
    LOGICAL, INTENT(OUT)         :: isFound    ! if object has been found
    INTEGER(INTG), INTENT(OUT)   :: indexFound ! and corresponding index

    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    INTEGER(INTG), INTENT(OUT) :: err

    ! Local variables    
    INTEGER(INTG) :: j, blockInd

    ENTERS("HashTable_GetKey",err,error,*999)

    isFound = .FALSE.
    indexFound = 0 !Returns meaningless value (0) if q not found: Always consider flag found
    ! copy procedure as in paper
    j = MyAlgModule(MyAlgModule(table%vecTKey(0)*q,table%p),table%n)

    blockInd = table%vecTKey(j)
    IF(blockInd /= 0) THEN ! if block is non empty
       indexFound = MyAlgModule(MyAlgModule(table%vecTKey(blockInd+1)*q,table%p),table%vecTKey(blockInd)**2)+1+blockInd
       IF(table%vecTKey(indexFound)==q) isFound = .TRUE.
    END IF

    EXITS("HashTable_GetKey")
    RETURN
999 ERRORSEXITS("HashTable_GetKey",err,error)
    RETURN 1

  END SUBROUTINE HashTable_GetKey

!==================================================

  ! Returns the INT value(s) from the array of lists at key location indexFound
  ! The list number(s) are stored in the array listNum
  SUBROUTINE HashTable_GetListArrayIntValue(table, indexFound, value, listNum, err,error,*) 

    ! Argument variables
    TYPE(HashTableType), POINTER :: table    

    INTEGER(INTG), INTENT(IN)  :: indexFound 
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: value(:,:) ! Allow for different sizes

    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: listNum(:) ! The lists where the values are found

    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    INTEGER(INTG), INTENT(OUT) :: err

    ! Local variables
    INTEGER(INTG) :: indexHash
    INTEGER(INTG), ALLOCATABLE :: listData(:)
    INTEGER(INTG) :: listTot, listAssociated, i, j

    ENTERS("HashTable_GetListArrayIntValue",err,error,*999)

    IF(ASSOCIATED(table)) THEN
      ! First get the index in the table
      indexHash = table%vecTVal(indexFound)

      ! Check if array of lists has been created
      IF(ALLOCATED(table%arrayOfListSVal)) THEN
        ! First list is the original list (listSVal). ListData has the dimension of the number of arrays.
        IF(ASSOCIATED(table%arrayOfListSVal(1)%ptr)) THEN
          ALLOCATE(listData(size(table%arrayOfListSVal,1)), STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate list data.",err,error,*999)  
        END IF 

        listTot = 0        ! The number of lists with INT values
        listAssociated = 0 ! The total number of ASSOCIATED lists

        DO i = 1,size(table%arrayOfListSVal,1)
          IF(ASSOCIATED(table%arrayOfListSVal(i)%ptr)) THEN
            IF(table%arrayOfListSVal(i)%ptr%DATA_TYPE==LIST_INTG_TYPE) THEN
              listData(i)     = table%arrayOfListSVal(i)%ptr%DATA_DIMENSION ! store the dimension of INT list
              listTot         = listTot +1 ! increase counter of INT lists
            ELSE  
              listData(i) = 0 ! append dimension 0
            END IF
            listAssociated = listAssociated +1 ! increase counter of ASSOCIATED lists
          ELSE
            EXIT ! Stop at first not-associated list
          END IF
        END DO

        IF (listTot == 0) THEN
          CALL FlagError("Could not find any data of the desired data type!",err,error,*999)
        END IF

        ALLOCATE(listNum(listTot), STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate listnum.",err,error,*999)  

        ! value is an array #list X max dimension
        ALLOCATE(value(listTot,maxval(listData(1:listAssociated),1)), STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate value.",err,error,*999)
        value = 0

        j = 1
        ! Loop over associated lists
        DO i = 1,listAssociated
          SELECT CASE(listData(i)) ! Dimension of data in the list
          CASE(0) ! Dimension 0 for not-INT list
             ! Do nothing
          CASE(1:)
            CALL List_ItemGet(table%arrayOfListSVal(i)%ptr, indexHash, value(j,1:listData(i)),err,error,*999)
            listNum(j) = i 
            j = j+1
          CASE DEFAULT
            CALL FlagError("Invalid dimension of list data!",err,error,*999)
          END SELECT
        END DO
       
      ELSE
        CALL FlagError("No extra values have been added",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Table is not associated",err,error,*999)
    END IF

    EXITS("HashTable_GetListArrayIntValue")
    RETURN
999 IF (ALLOCATED(listData)) DEALLOCATE(listData)
    ERRORSEXITS("HashTable_GetListArrayIntValue",err,error)
    RETURN 1

  END SUBROUTINE HashTable_GetListArrayIntValue

!==================================================

  ! Returns the REAL SP value(s) from the array of lists at key location indexFound
  ! The list number(s) are stored in the array listNum
  SUBROUTINE HashTable_GetListArrayRealSpValue(table, indexFound, value, listNum, err,error,*) 

    ! Argument variables
    TYPE(HashTableType), POINTER :: table    

    INTEGER(INTG), INTENT(IN)  :: indexFound 
    REAL(SP), ALLOCATABLE, INTENT(OUT) :: value(:,:) ! Allow for different sizes

    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: listNum(:) ! The lists where the values are found

    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    INTEGER(INTG), INTENT(OUT) :: err

    ! Local variables
    INTEGER(INTG) :: indexHash
    INTEGER(INTG), ALLOCATABLE :: listData(:)
    INTEGER(INTG) :: listTot, listAssociated, i, j

    ENTERS("HashTable_GetListArrayRealSpValue",err,error,*999)

    IF(ASSOCIATED(table)) THEN
      ! First get the index in the table
      indexHash = table%vecTVal(indexFound)

      ! Check if array of lists has been created
      IF(ALLOCATED(table%arrayOfListSVal)) THEN
        ! First list is the original list (listSVal). ListData has the dimension of the number of arrays.
        IF(ASSOCIATED(table%arrayOfListSVal(1)%ptr)) THEN
          ALLOCATE(listData(size(table%arrayOfListSVal,1)), STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate list data.",err,error,*999)  
        END IF 

        listTot = 0        ! The number of lists with REAL SP values
        listAssociated = 0 ! The total number of ASSOCIATED lists

        DO i = 1,size(table%arrayOfListSVal,1)
          IF(ASSOCIATED(table%arrayOfListSVal(i)%ptr)) THEN
            IF(table%arrayOfListSVal(i)%ptr%DATA_TYPE==LIST_SP_TYPE) THEN
              listData(i)     = table%arrayOfListSVal(i)%ptr%DATA_DIMENSION ! store the dimension of REAL SP list
              listTot         = listTot +1 ! increase counter of REAL SP lists
            ELSE  
              listData(i) = 0 ! append dimension 0
            END IF
            listAssociated = listAssociated +1 ! increase counter of ASSOCIATED lists
          ELSE
            EXIT ! Stop at first not-associated list
          END IF
        END DO

        IF (listTot == 0) THEN
          CALL FlagError("Could not find any data of the desired data type!",err,error,*999)
        END IF

        ALLOCATE(listNum(listTot), STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate listnum.",err,error,*999)  

        ! value is an array #list X max dimension
        ALLOCATE(value(listTot,maxval(listData(1:listAssociated),1)), STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate value.",err,error,*999)
        value = 0

        j = 1
        ! Loop over associated lists
        DO i = 1,listAssociated
          SELECT CASE(listData(i)) ! Dimension of data in the list
          CASE(0) ! Dimension 0 for not-REAL SP list
             ! Do nothing
          CASE(1:)
            CALL List_ItemGet(table%arrayOfListSVal(i)%ptr, indexHash, value(j,1:listData(i)),err,error,*999)
            listNum(j) = i 
            j = j+1
          CASE DEFAULT
            CALL FlagError("Invalid dimension of list data!",err,error,*999)
          END SELECT
        END DO
       
      ELSE
        CALL FlagError("No extra values have been added",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Table is not associated",err,error,*999)
    END IF

    EXITS("HashTable_GetListArrayRealSpValue")
    RETURN
999 IF (ALLOCATED(listData)) DEALLOCATE(listData)
    ERRORSEXITS("HashTable_GetListArrayRealSpValue",err,error)
    RETURN 1

  END SUBROUTINE HashTable_GetListArrayRealSpValue


!==================================================
  ! Returns the SP value at index found
  SUBROUTINE HashTable_GetRealSpValue(table, indexFound, value, err,error,*) 

    ! Argument variables
    TYPE(HashTableType), POINTER :: table    

    INTEGER(INTG), INTENT(IN)  :: indexFound 
    REAL(SP), ALLOCATABLE, INTENT(OUT) :: value(:) ! Allow for different sizes

    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    INTEGER(INTG), INTENT(OUT) :: err

    ! Local variables
    INTEGER(INTG) :: indexHash, listDataDimension
    REAL(SP) :: value1D

    ENTERS("HashTable_GetRealSpValue",err,error,*999)

    IF(ASSOCIATED(table)) THEN

         ! Use list (original list SVal) and the index stored in TVal
         IF(ASSOCIATED(table%listSVal)) THEN

           ! First get the index
           indexHash = table%vecTVal(indexFound)
           listDataDimension = table%listSVal%DATA_DIMENSION

           ALLOCATE(value(listDataDimension), STAT=err)
           IF(err/=0) CALL FlagError("Could not allocate value",err,error,*999)  

           SELECT CASE(listDataDimension)
           CASE(1)
             CALL List_ItemGet(table%listSVal, indexHash, value1D,err,error,*999) 
             value(1) = value1D
           CASE(2:)
             CALL List_ItemGet(table%listSVal, indexHash, value,err,error,*999)
           CASE DEFAULT
             CALL FlagError("Invalid dimension of list data!",err,error,*999)
           END SELECT
         ELSE
           CALL FlagError("Values are not there!",err,error,*999)
         END IF
    ELSE
       CALL FlagError("Table is not associated",err,error,*999)
    END IF
    
    EXITS("HashTable_GetRealSpValue")
    RETURN
999 ERRORSEXITS("HashTable_GetRealSpValue",err,error)
    RETURN 1

  END SUBROUTINE HashTable_GetRealSpValue

!==================================================
  ! Returns the INT value at index found
  SUBROUTINE HashTable_GetIntValue(table, indexFound, value, err,error,*) 

    ! Argument variables
    TYPE(HashTableType), POINTER :: table    

    INTEGER(INTG), INTENT(IN)  :: indexFound 
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: value(:) ! Allow for different sizes

    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    INTEGER(INTG), INTENT(OUT) :: err

    ! Local variables
    INTEGER(INTG)  :: indexHash, listDataDimension
    INTEGER(INTG)  :: value1D

    ENTERS("HashTable_GetRealSpValue",err,error,*999)

    IF(ASSOCIATED(table)) THEN

         ! Use list (original list SVal) and the index stored in TVal
         IF(ASSOCIATED(table%listSVal)) THEN

           ! First get the index
           indexHash = table%vecTVal(indexFound)
           listDataDimension = table%listSVal%DATA_DIMENSION

           ALLOCATE(value(listDataDimension), STAT=err)
           IF(err/=0) CALL FlagError("Could not allocate value",err,error,*999)  

           SELECT CASE(listDataDimension)
           CASE(1)
             CALL List_ItemGet(table%listSVal, indexHash, value1D,err,error,*999) 
             value(1) = value1D
           CASE(2:)
             CALL List_ItemGet(table%listSVal, indexHash, value,err,error,*999)
           CASE DEFAULT
             CALL FlagError("Invalid dimension of list data!",err,error,*999)
           END SELECT
         ELSE
           CALL FlagError("Values are not there!",err,error,*999)
         END IF
    ELSE
       CALL FlagError("Table is not associated",err,error,*999)
    END IF
    
    EXITS("HashTable_GetIntValue")
    RETURN
999 ERRORSEXITS("HashTable_GetIntValue",err,error)
    RETURN 1

  END SUBROUTINE HashTable_GetIntValue

!==================================================

  ! Build the list of INT values according to input
  SUBROUTINE HashTable_ListIntValuesSet (table, vecSVal, addElements, err,error,*)

    ! Argument variables
    TYPE(HashTableType), POINTER :: table    

    INTEGER(INTG), INTENT(IN) :: vecSVal(:,:)  ! array of values

    LOGICAL, INTENT (IN) :: addElements

    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    INTEGER(INTG), INTENT(OUT) :: err

    ! Local variables
    ! related to Lists
    INTEGER(INTG) :: listSizeX, listSizeY, i
    INTEGER(INTG) :: vecS1DVal(size(vecSVal,2)) ! in case only one row of values => 1D List

    ENTERS("HashTable_ListIntValuesSet",err,error,*999)

    ! The size of the input: the values correspond to the rows of the input array
    ! Dimension of data
    listSizeY = size(vecSVal,1) 
    ! Number of entries (must correspond to number of keys, will be checked on insert)
    listSizeX = size(vecSVal,2)

    CALL HashTable_ListInitialise (table, LIST_INTG_TYPE, listsizeX, listSizeY, addElements, err,error,*999)

    ! One-row case
    IF (listSizeY==1) vecS1DVal = vecSVal(1,:) 

    ! Add items (might be adding to an already filled list)
    ! Resize before adding NEW items (otherwise ItemAdd would add unnecessary space)
    IF (addElements) CALL List_Resize(table%listSVal, listSizeX, err, error, *999) 
    DO i = 1,listSizeX
        IF (listSizeY==1) THEN
          CALL List_ItemAdd(table%listSVal,vecS1DVal(i), err, error, *999)
        ELSE
          CALL List_ItemAdd(table%listSVal,vecSVal(:,i), err, error, *999)
        END IF
    END DO

    EXITS("HashTable_ListIntValuesSet")
    RETURN
999 ERRORSEXITS("HashTable_ListIntValuesSet",err,error)
    RETURN 1 

END SUBROUTINE

!==================================================


  ! Build the list of SP values according to input
  SUBROUTINE HashTable_ListRealSpValuesSet (table, vecSVal, addElements, err,error,*)

    ! Argument variables
    TYPE(HashTableType), POINTER :: table    

    REAL(SP), INTENT(IN) :: vecSVal(:,:)  ! array of values

    LOGICAL, INTENT (IN) :: addElements

    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    INTEGER(INTG), INTENT(OUT) :: err

    ! Local variables
    ! related to Lists
    INTEGER(INTG) :: listSizeX, listSizeY, i
    REAL(SP) :: vecS1DVal(size(vecSVal,2)) ! in case only one row of values => 1D List

    ENTERS("HashTable_ListRealSpValuesSet",err,error,*999)

    ! The size of the input: the values correspond to the rows of the input array
    ! Dimension of data
    listSizeY = size(vecSVal,1) 
    ! Number of entries (must correspond to number of keys, will be checked on insert)
    listSizeX = size(vecSVal,2)

    CALL HashTable_ListInitialise (table, LIST_SP_TYPE, listsizeX, listSizeY, addElements, err,error,*999)

    ! One-row case
    IF (listSizeY==1) vecS1DVal = vecSVal(1,:) 

    ! Add items (can be adding to an already filled list)
    ! Resize before adding NEW items (otherwise ItemAdd would add unnecessary space)
    IF (addElements) CALL List_Resize(table%listSVal, listSizeX, err, error, *999) 
    DO i = 1,listSizeX
        IF (listSizeY==1) THEN
         CALL List_ItemAdd(table%listSVal,vecS1DVal(i), err, error, *999)
        ELSE
         CALL List_ItemAdd(table%listSVal,vecSVal(:,i), err, error, *999)
        END IF
    END DO

    EXITS("HashTable_ListRealSpValuesSet")
    RETURN
999 ERRORSEXITS("HashTable_ListRealSpValuesSet",err,error)
    RETURN 1 

END SUBROUTINE

!==================================================

  ! Initialise the list of values independent of data type
  SUBROUTINE HashTable_ListInitialise (table, dataType, listsizeX, listSizeY, addElements, err,error,*)

    ! Argument variables
    TYPE(HashTableType), POINTER :: table    
    INTEGER(INTG), INTENT(IN)    :: dataType, listsizeX, listSizeY
    LOGICAL, INTENT(IN)          :: addElements

    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    INTEGER(INTG), INTENT(OUT) :: err


    ENTERS("HashTable_ListInitialise",err,error,*999)

    IF(ASSOCIATED(table)) THEN
      IF(table%hashTableFinished) THEN
       ! Do nothing 
      ELSE
        CALL FlagError("The table has not been finished.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Table not associated.",err,error,*999)
    END IF

    ! Check if extra values are there
    IF (.NOT. addElements) THEN
      ! Deallocate extra values if creating a new table (and go on)
      IF (ALLOCATED(table%arrayOfListSVal)) DEALLOCATE(table%arrayOfListSVal)
    ELSE
      IF (ALLOCATED(table%arrayOfListSVal)) THEN
        DEALLOCATE(table%arrayOfListSVal)
        CALL FlagError("Cannot add new elements to a list with extra values! (FIX?)",err,error,*999)
      END IF
    END IF
    

    IF (.NOT. addElements .AND. ASSOCIATED(table%listSVal)) CALL List_Destroy(table%listSVal, err, error, *999) ! destroy the existing list 
    ! If not already existing, create a new one
    IF (.NOT. ASSOCIATED(table%listSVal)) THEN
    ! Pointers nullified in hash initialise!
      CALL List_CreateStart(table%listSVal, err, error, *999)
!   Defaults in initialise:
!      LIST%LIST_FINISHED=.FALSE.
!      LIST%MUTABLE=.FALSE.
!      LIST%NUMBER_IN_LIST=0
!      LIST%DATA_DIMENSION=1
!      LIST%INITIAL_SIZE=10
!      LIST%SIZE=0
!      LIST%DATA_TYPE=LIST_INTG_TYPE
!      LIST%KEY_DIMENSION=1
!      LIST%SORT_ORDER=LIST_SORT_ASCENDING_TYPE
!      LIST%SORT_METHOD=LIST_HEAP_SORT_METHOD

    ! Let us create a list of intg of dim listSizeY
      CALL List_DataDimensionSet(table%listSVal,listSizeY, err, error, *999)
    ! Initial size listSizeX
      CALL List_InitialSizeSet(table%listSVal, listSizeX, err, error, *999)
      CALL List_DataTypeSet(table%listSVal, dataType, err, error, *999) ! Ok default
      CALL List_CreateFinish(table%listSVal, err, error, *999)
    END IF

    EXITS("HashTable_ListInitialise")
    RETURN
999 ERRORSEXITS("HashTable_ListInitialise",err,error)
    RETURN 1 

END SUBROUTINE

!==================================================

  ! Build the list of additional INT values according to input
  SUBROUTINE HashTable_AdditionalIntValuesSet (table, vecSVal, dataType, err,error,*)

    ! Argument variables
    TYPE(HashTableType), POINTER :: table    

    INTEGER(INTG), INTENT(IN) :: vecSVal(:,:)  ! array of values
    INTEGER(INTG) :: dataType  ! data type

    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    INTEGER(INTG), INTENT(OUT) :: err

    ! Local variables
    ! related to Lists
    INTEGER(INTG) :: listSizeX, listSizeY, listNum, j
    INTEGER(INTG) :: vecS1DVal(size(vecSVal,2)) ! in case only one row of values => 1D List

    ENTERS("HashTable_AdditionalIntValuesSet",err,error,*999)

    ! The size of the input: the values correspond to the rows of the input array
    ! Dimension of data
    listSizeY = size(vecSVal,1) 
    ! Number of entries (must correspond to number of keys, will be checked on insert)
    listSizeX = size(vecSVal,2)

    ! Init new list and return number in the list array
    CALL HashTable_AdditionalValuesInit (table, dataType, listSizeX, listSizeY, listNum, err,error,*999)

    IF (listSizeY==1) vecS1DVal = vecSVal(1,:) 

    ! Add items 
    DO j = 1,listSizeX
      IF (listSizeY==1) THEN
        CALL List_ItemAdd(table%arrayOfListSVal(listNum)%ptr,vecS1DVal(j), err, error, *999)
      ELSE
        CALL List_ItemAdd(table%arrayOfListSVal(listNum)%ptr,vecSVal(:,j), err, error, *999)
      END IF
    END DO

    EXITS("HashTable_AdditionalIntValuesSet")
    RETURN
999 ERRORSEXITS("HashTable_AdditionalIntValuesSet",err,error)
    RETURN 1 

END SUBROUTINE

!==================================================

  ! Build the list of additional REAL SP values according to input
  SUBROUTINE HashTable_AdditionalRealSpValuesSet (table, vecSVal, dataType, err,error,*)

    ! Argument variables
    TYPE(HashTableType), POINTER :: table    

    REAL(SP), INTENT(IN) :: vecSVal(:,:)  ! array of values
    INTEGER(INTG) :: dataType  ! data type

    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    INTEGER(INTG), INTENT(OUT) :: err

    ! Local variables
    ! related to Lists
    INTEGER(INTG) :: listSizeX, listSizeY, listNum, j
    REAL(SP) :: vecS1DVal(size(vecSVal,2)) ! in case only one row of values => 1D List

    ENTERS("HashTable_AdditionalRealSpValuesSet",err,error,*999)

    ! The size of the input: the values correspond to the rows of the input array
    ! Dimension of data
    listSizeY = size(vecSVal,1) 
    ! Number of entries (must correspond to number of keys, will be checked on insert)
    listSizeX = size(vecSVal,2)

    ! Init new list and return number in the list array
    CALL HashTable_AdditionalValuesInit (table, dataType, listSizeX, listSizeY, listNum, err,error,*999)

    IF (listSizeY==1) vecS1DVal = vecSVal(1,:) 

    ! Add items 
    DO j = 1,listSizeX
      IF (listSizeY==1) THEN
        CALL List_ItemAdd(table%arrayOfListSVal(listNum)%ptr,vecS1DVal(j), err, error, *999)
      ELSE
        CALL List_ItemAdd(table%arrayOfListSVal(listNum)%ptr,vecSVal(:,j), err, error, *999)
      END IF
    END DO

    EXITS("HashTable_AdditionalRealSpValuesSet")
    RETURN
999 ERRORSEXITS("HashTable_AdditionalRealSpValuesSet",err,error)
    RETURN 1 

END SUBROUTINE

!==================================================

  ! Initialise additional list of values (returns newly allocated list)
  SUBROUTINE HashTable_AdditionalValuesInit (table, dataType, listSizeX, listSizeY, listNum, err,error,*)

    ! Argument variables
    TYPE(HashTableType), POINTER :: table    

    INTEGER(INTG), INTENT(IN)  :: dataType, listSizeX, listSizeY
    INTEGER(INTG), INTENT(OUT) :: listNum

    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    INTEGER(INTG), INTENT(OUT) :: err

    ! Local variables
    INTEGER(INTG) :: i, sizeInit

    ENTERS("HashTable_AdditionalValuesInit",err,error,*999)
   
   ! New array of lists if none allocated yet
    IF (.NOT. ALLOCATED(table%arrayOfListSVal)) THEN

      ! Allocate a sufficient dimension of (empty) lists
      ALLOCATE(table%arrayOfListSVal(HASH_INITIAL_ADDITIONAL_VALUES_SIZE),STAT=err)   
      DO i=1,HASH_INITIAL_ADDITIONAL_VALUES_SIZE
        table%arrayOfListSVal(i)%ptr=>NULL() ! Define the ptr to list as not associated
      END DO
    
      CALL List_NumberOfItemsGet(table%listSVal, sizeInit, err, error, *999)
 
      IF (sizeInit /= listSizeX) THEN
        DEALLOCATE(table%arrayOfListSVal)
        CALL FlagError("Adding values of different size than keys in table!!", err, error,*999)
      END IF

      ! First list in array is the one already provided in listSval (at table creation)
      ! No need to initialise / finish new pointer. Just allocate.
      ALLOCATE(table%arrayOfListSVal(1)%ptr,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate array list 1.", err, error,*999)
      table%arrayOfListSVal(1)%ptr = table%listSval 

    END IF

    ! Allocate an additional list (look for the first free spot)
    listNum = 2
    DO
      IF (size(table%arrayOfListSVal)>=listNum) THEN
        IF (.NOT. ASSOCIATED(table%arrayOfListSVal(listNum)%ptr)) THEN
          CALL List_CreateStart(table%arrayOfListSVal(listNum)%ptr, err, error, *999)
!   Defaults in initialise:
!      LIST%LIST_FINISHED=.FALSE.
!      LIST%MUTABLE=.FALSE.
!      LIST%NUMBER_IN_LIST=0
!      LIST%DATA_DIMENSION=1
!      LIST%INITIAL_SIZE=10
!      LIST%SIZE=0
!      LIST%DATA_TYPE=LIST_INTG_TYPE
!      LIST%KEY_DIMENSION=1
!      LIST%SORT_ORDER=LIST_SORT_ASCENDING_TYPE
!      LIST%SORT_METHOD=LIST_HEAP_SORT_METHOD

       ! Let us create a list of intg of dim listSizeY
          CALL List_DataDimensionSet(table%arrayOfListSVal(listNum)%ptr, listSizeY, err, error, *999)
          CALL List_DataTypeSet     (table%arrayOfListSVal(listNum)%ptr, dataType, err, error, *999) 
       ! Initial size listSizeX
          CALL List_InitialSizeSet  (table%arrayOfListSVal(listNum)%ptr, listSizeX, err, error, *999)
          CALL List_CreateFinish    (table%arrayOfListSVal(listNum)%ptr, err, error, *999)
          EXIT
        END IF
        listNum = listNum+1
      ELSE
        CALL FlagError("Number of additional lists exceeds initial allocation!",err,error,*999)
      END IF      
    END DO

    EXITS("HashTable_AdditionalValuesInit")
    RETURN
999 ERRORSEXITS("HashTable_AdditionalValuesInit",err,error)
    RETURN 1 

END SUBROUTINE

!==================================================

    ! Create the list of values and call the hash table create routine 
    SUBROUTINE HashTable_IntValuesSetAndInsert(table,  vecSKey, vecSVal, addElements, err,error,*)
    ! Argument variables
    TYPE(HashTableType), POINTER :: table    

    INTEGER(INTG), INTENT(IN) :: vecSKey(:)    ! array of keys
    INTEGER(INTG), INTENT(IN) :: vecSVal(:,:)  ! array of values

    LOGICAL, INTENT(IN) :: addElements ! adding to existing or new table

    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    INTEGER(INTG), INTENT(OUT) :: err

    ENTERS("HashTable_IntMultValuesSetAndInsert",err,error,*999)

    ! Create the list of values (if empty)
    CALL HashTable_ListValuesSet(table, vecSVal, addElements, err,error,*999)
    ! Compute the table based on the key array
    CALL HashTable_Insert(table, vecSKey, addElements, err,error,*999)

    EXITS("HashTable_IntMultValuesSetAndInsert")
    RETURN
999 ERRORSEXITS("HashTable_IntMultValuesSetAndInsert",err,error)
    RETURN 1 

    END SUBROUTINE

!==================================================

    ! Create the list of values and call the hash table create routine 
    SUBROUTINE HashTable_RealSpValuesSetAndInsert(table,  vecSKey, vecSVal, addElements, err,error,*)
    ! Argument variables
    TYPE(HashTableType), POINTER :: table    

    INTEGER(INTG), INTENT(IN) :: vecSKey (:)  ! array of keys
    REAL(SP), INTENT(IN) :: vecSVal(:,:)      ! array of values

    LOGICAL, INTENT(IN) :: addElements ! adding to existing or new table

    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    INTEGER(INTG), INTENT(OUT) :: err

    ENTERS("HashTable_RealSpMultValuesSetAndInsert",err,error,*999)

    ! Create the list of values (if empty)
    CALL HashTable_ListValuesSet(table, vecSVal, addElements, err,error,*999)
    ! Compute the table based on the key array
    CALL HashTable_Insert(table, vecSKey, addElements, err,error,*999)

    EXITS("HashTable_RealSpMultValuesSetAndInsert")
    RETURN
999 ERRORSEXITS("HashTable_RealSpMultValuesSetAndInsert",err,error)
    RETURN 1 

    END SUBROUTINE


!==================================================


  ! Build the table according to Fredman 1984: This is the most important subroutine
   SUBROUTINE HashTable_Insert(table, vecSKey, addElements, err,error,*)
    ! Argument variables
    TYPE(HashTableType), POINTER :: table    

    INTEGER(INTG), INTENT(IN) :: vecSKey(:)    ! array of keys

    LOGICAL, INTENT(IN) :: addElements ! adding to existing or new table

    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    INTEGER(INTG), INTENT(OUT) :: err

    ! Local variables
    ! related to subset W
    TYPE (HashListArrayType), POINTER  :: subsetW=>NULL() ! W as a pointer to an array of lists
    TYPE(HashLinkedListItemType) :: subsetWListItem
    INTEGER(INTG) :: cardSubsetW, storeCardSubsetW, checkSumW, checkSquaredSumW, countFilledSubsetW, indexSubsetW
    ! related to array T
    INTEGER(INTG) :: indexTLoop, sizeOfT, prevIndexT
    INTEGER(INTG), ALLOCATABLE :: indexTArray (:), tmpKey(:)
    ! related to Lists
    INTEGER(INTG) :: listSize

    TYPE(VARYING_STRING) :: dummyError !<The error string.
    INTEGER(INTG)        :: dummyErr
    ! generic
    INTEGER(INTG) :: i,k,j
    LOGICAL :: isFound = .FALSE.

    ENTERS("HashTable_Insert",err,error,*999)

    IF(addElements) THEN ! add to existing table
     
      IF(.NOT. ASSOCIATED(table))        CALL FlagError("Table should be already created.",err,error,*999)
      IF(.NOT. table%hashTableFinished)  CALL FlagError("Table should be already initialised.",err,error,*999) 
      IF(ALLOCATED(table%vecTKey)) DEALLOCATE(table%vecTKey)
      IF(ALLOCATED(table%vecTVal)) DEALLOCATE(table%vecTVal)

      IF(.NOT. ALLOCATED(table%vecSKey)) CALL FlagError("Cannot put new keys in an empty table.",err,error,*999)
      IF(.NOT. ASSOCIATED(table%listSVal)) CALL FlagError("Cannot put new values in an empty table.",err,error,*999)

      table%n = size(table%vecSKey,1)  ! old size of table

      ALLOCATE(tmpKey(table%n+size(vecSKey,1)), STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate tmp in table.",err,error,*999)
      tmpKey(1:table%n)  = table%vecSKey
      tmpKey(table%n+1:) = vecSKey(1:)
      CALL MOVE_ALLOC (tmpKey, table%vecSKey)

      table%n = size(table%vecSKey,1) ! new size of table

    ELSE ! new table

      IF(ASSOCIATED(table)) THEN
        IF(table%hashTableFinished) THEN
      !   DEALLOCATE already existing memory in table
          IF(ALLOCATED(table%vecTKey)) DEALLOCATE(table%vecTKey)
          IF(ALLOCATED(table%vecTVal)) DEALLOCATE(table%vecTVal)
          IF(ALLOCATED(table%vecSKey)) DEALLOCATE(table%vecSKey)
        ELSE
          CALL FlagError("The table has not been finished.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Table not associated.",err,error,*999)
      END IF

      table%n = size(vecSKey,1) 

      ALLOCATE(table%vecSKey(table%n), STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate vecSKey in table",err,error,*999)
      table%vecSKey(1:table%n) = vecSKey(1:table%n)

    END IF


    ! Check that size of list corresponds to key size
     IF (ASSOCIATED(table%listSVal)) THEN
       CALL List_NumberOfItemsGet(table%listSVal, listSize,err, error, *999)
       IF(table%n /= listSize)  CALL FlagError("Keys and list of values have different size!",err,error,*999)
     END IF       
  

    ! Compute p as the first prime number larger than biggest key
    table%p = maxval(table%vecSKey,1)+1
    SELECT CASE(table%p)
    CASE (1)
    CASE (2)
    CASE DEFAULT
     i=2
     DO !i = 2, table%p-1
        IF(MOD(table%p,i) == 0)  THEN
          table%p = table%p+1
          i = 2
          CYCLE
        END IF
        i = i+1
        IF(i == table%p) EXIT
     END DO
    END SELECT 

    ! Initialise w as an array of size n
    CALL HashListArray_CreateStart(subsetW, err,error,*999)
    subsetW%vecLen = table%n
    CALL HashListArray_CreateFinish(subsetW, err,error,*999)

    ! Compute the W_j
    k=2 ! use 2 as start as in paper (but 1 also works for paper case!!)

    DO

     DO i=1,table%n

      indexSubsetW = MyAlgModule(MyAlgModule((k*table%vecSKey(i)),table%p),table%n)
      !PRINT *, indexSubsetW
      ! CALL HashLinkedList_AddData(subsetW%vec(indexSubsetW-1)%ptr,table%vecSKey(i),table%vecSVal(i), isFound, err, error, *999) ! numbered from 0!
      ! Change! Store INDEX i instead of value!
      CALL HashLinkedList_AddData(subsetW%vec(indexSubsetW-1)%ptr,table%vecSKey(i), i, isFound, err, error, *999) ! numbered from 0!
      IF(isFound) CALL FlagError("Duplicate key, error!",err,error,*999) 
         !PRINT *, "Success with ", indexSubsetW
     END DO

     checkSquaredSumW =0
     checkSumW   = 0
     countFilledSubsetW  = 0
     DO i=1,table%n
      CALL HashLinkedList_Size(subsetW%vec(i-1)%ptr, cardSubsetW, err, error,*999) 
      IF(cardSubsetW /= 0) THEN
      checkSquaredSumW = checkSquaredSumW + cardSubsetW**2
      checkSumW    = checkSumW    + cardSubsetW
      countFilledSubsetW   = countFilledSubsetW +1
      END IF
     END DO

    IF(checkSquaredSumW<3*table%n) THEN ! k is correct
     EXIT
    ELSE
     k = k+1 ! try another k
     IF(k>table%p) CALL FlagError("Suitable k could not be found!",err,error,*999) 
    END IF

   END DO

   ! Allocate T: T(0)=k, T(1),...,T(n) access indices, T(*) |W_j| T(**) k' T(***) = W_j
   sizeOfT = 1+table%n+2*countFilledSubsetW+checkSquaredSumW
   ALLOCATE (table%vecTKey(0:(sizeOfT-1)),STAT=err)
   IF(err/=0) CALL FlagError("Could not allocate vecTKey in table.",err,error,*999)
   ALLOCATE (table%vecTVal(0:(sizeOfT-1)),STAT=err)
   IF(err/=0) CALL FlagError("Could not allocate vecTVal in table.",err,error,*999)
   table%vecTKey = 0
   table%vecTVal = 0
   table%vecTKey(0) = k
   IF(sizeOfT<=table%n+1) CALL FlagError("No keys to be assigned!",err,error,*999) 

   prevIndexT = table%n

   DO i=1,table%n

     ! Cardinality of block W_i 
     CALL HashLinkedList_Size(subsetW%vec(i-1)%ptr, storeCardSubsetW, err,error,*999)

     IF(storeCardSubsetW == 0) THEN
       table%vecTKey(i) = 0 ! Block is empty
     ELSE

       IF(prevIndexT == table%n) THEN
        table%vecTKey(i) = table%n+1 ! First non-zero entry is just index that follows n
       ELSE 
        ! Compute card of 
        ! CALL HashLinkedList_Size(subsetW%vec(prevIndexT-1), cardSubsetW, err)    
        ! replaced with cardSubsetW = storeCardSubsetW below!!!!   
        table%vecTKey(i) = table%vecTKey(prevIndexT)+2+cardSubsetW**2    
       END IF

       table%vecTKey(table%vecTKey(i))   = storeCardSubsetW ! First  element of each block: card of W
       table%vecTKey(table%vecTKey(i)+1) = 1            ! Second element of each block: k' (set to 1)

       ALLOCATE (indexTArray(storeCardSubsetW), STAT=err)  ! indexTArray contains the indices where keys will be stored in T for each block    
       IF(err/=0) CALL FlagError("Could not allocate indexTArray in table.",err,error,*999)
       subsetWListItem = subsetW%vec(i-1)%ptr%root            ! take the first item in the block
       j = 1 ! init j for do loop

       DO 
         ! function given by corollary 2
         ! gives LOCAL index
         ! +2 translation
         indexTArray(j) =  MyAlgModule(MyAlgModule(table%vecTKey(table%vecTKey(i)+1)*subsetWListItem%key,table%p) &
                             & ,storeCardSubsetW**2)+2
         IF((j==1).OR. any(indexTArray(1:(j-1))/=indexTArray(j))) THEN ! insert items until conflict

            table%vecTKey(indexTArray(j)+table%vecTKey(i)-1) = subsetWListItem%key 
            table%vecTVal(indexTArray(j)+table%vecTKey(i)-1) = subsetWListItem%value

            IF(ASSOCIATED(subsetWListItem%next)) subsetWListItem = subsetWListItem%next
            
            IF(j==storeCardSubsetW) EXIT ! everything filled correctly
            j = j+1  ! next key
            
         ELSE ! conflict: increase k' and start over
             table%vecTKey(table%vecTKey(i)+1) =table%vecTKey(table%vecTKey(i)+1)+1 

             DO indexTLoop =1,size(indexTArray,1) 
                table%vecTKey(indexTArray(indexTLoop) + table%vecTKey(i)-1) = 0   ! reset to zero what already filled (keys)
                table%vecTVal(indexTArray(indexTLoop) + table%vecTKey(i)-1) = 0   ! reset to zero what already filled (values)
             END DO 

             indexTArray = 0                       ! reset to zero indices
             j = 1                              ! reset j
             subsetWListItem = subsetW%vec(i-1)%ptr%root  ! go back to list start
         END IF
 
       END DO
       
       !PRINT *, "IND IN T",  indexTArray
       DEALLOCATE(indexTArray)

       prevIndexT = i
       cardSubsetW = storeCardSubsetW
     END IF

   END DO

   ! Destroy the list array of Ws
   CALL HashListArray_Finalise(subsetW, dummyErr,dummyError,*998)

   ! DEALLOCATE local memory
   IF(ALLOCATED(indexTArray)) DEALLOCATE(indexTArray)

   ! Diagnostics instead of print
   IF(diagnostics1) THEN

     CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE, "The number of keys n:    ", table%n, err, error, *999)

     CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE, "The size of the array T: ", size(table%vecTKey,1), err, error, *999)

     CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE, "The value of k = T(0):   ", table%vecTkey(0), err, error, *999)

     CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE, 2,1, table%n+1, 4,4, table%vecTKey, & 
       & '("T(1) to T(n) (indices to T blocks):", 4(x,I8))','(35x,4(x,I8))', err, error, *999)

     CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE, table%n+2,1, size(table%vecTKey,1), 4,4, table%vecTKey, & 
       & '("T(n+1) to T(end):                  ", 4(x,I8))','(35x,4(x,I8))', err, error, *999)

!       print *, "T vector n    ", table%vecTKey 
  !>Writes the given INTEGER VECTOR to the given output stream specified by ID. The FIRST_FORMAT is the format initially used, followed by the REPEAT_FORMAT which is repeated as many times as necessary. NUMBER_FIRST is the number of data items in the FIRST_FORMAT and NUMBER_REPEAT is the number of data items in the REPEAT_FORMAT. FIRST_IDX and LAST_IDX are the extents of the data and DELTA is the NUMBER of indices to skip for each index.
!          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_LINE%BASIS%NUMBER_OF_NODES,4,4,DOMAIN_LINE%NODES_IN_LINE, &
!           & '("        Nodes in line        :",4(X,I8))','(30X,4(X,I8))',err,error,*999)

!  SUBROUTINE WRITE_STRING_VECTOR_DP(ID,FIRST_IDX,DELTA,LAST_IDX,                               NUMBER_FIRST,NUMBER_REPEAT,VECTOR,FIRST_FORMAT,REPEAT_FORMAT, &
!    & err,error,*)

       !print *, "n    ", table%n
       !print *, "T vector n    ", table%vecTKey (1:table%n)
       !print *, "T vector rest ", table%vecTKey (table%n+1:)
       !print *, "T value vector rest ", table%vecTVal (table%n+1:)
       !print *, "of size ", size(table%vecTKey)
  END IF

    EXITS("HashTable_Insert")
    RETURN
999 CALL HashTable_Finalise(Table,dummyErr,dummyError,*998)
998 ERRORSEXITS("HashTable_Insert",err,error)
    RETURN 1

  END SUBROUTINE HashTable_Insert

!==================================================

  !>Finalises a hash table and deallocates all memory.
  SUBROUTINE HashTable_Finalise(table,err,error,*)    

    !Argument Variables
    TYPE (HashTableType), POINTER :: table !<The table to finalise
    INTEGER(INTG), INTENT(OUT)    :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("HashTable_Finalise",err,error,*999)
 
   IF(ASSOCIATED(table)) THEN
      IF(ALLOCATED(table%vecTKey)) DEALLOCATE(table%vecTKey)
      IF(ALLOCATED(table%vecTVal)) DEALLOCATE(table%vecTVal)
      IF(ALLOCATED(table%vecSKey)) DEALLOCATE(table%vecSKey)
      DEALLOCATE(table)
   END IF

    EXITS("HashTable_Finalise")
    RETURN
999 ERRORSEXITS("HashTable_Finalise",err,error)
    RETURN 1

  END SUBROUTINE HashTable_Finalise

!==================================================

! Computes modulo as in Fredman84 (returns divisor instead of 0)
FUNCTION MyAlgModule(a,b) RESULT(c)

    INTEGER(INTG), INTENT(IN)  :: a,b
    INTEGER(INTG) :: c

    c = MOD(a,b)
    IF(c==0) c = b

END FUNCTION MyAlgModule


END MODULE HashRoutines
