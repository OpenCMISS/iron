!> HLinked List for 2xinteger data TYPE
! + Hash functions

MODULE hash_routines

  USE BaseRoutines
  USE INPUT_OUTPUT
  USE CONSTANTS
  USE KINDS
  USE ISO_VARYING_STRING
  USE STRINGS 
  USE ComputationEnvironment

#include "macros.h"  

  IMPLICIT NONE
  
! delete me!
!  INTEGER, PARAMETER :: tbl_size = 500 ! a "big" default size for the table

  PRIVATE  ! by default

  ! TYPEs

  ! Items for Linked Lists with 2 entries
  TYPE HashLinkedListItemType
    integer(INTG) :: key    
    integer(INTG) :: value
    TYPE(HashLinkedListItemType),pointer :: next => NULL()
  END TYPE

  ! Linked Lists
  TYPE HashLinkedListType
    TYPE(HashLinkedListItemType),pointer :: root => NULL()
    TYPE(HashLinkedListItemType),pointer :: last => NULL()
  END TYPE

  ! Array of Linked Lists
  TYPE HashListArrayType
     TYPE(HashLinkedListType), ALLOCATABLE :: vec(:)
     INTEGER(INTG)                         :: vec_len
     LOGICAL(INTG)                         :: is_finished
  END TYPE HashListArrayType

  ! The Hash Table
  TYPE HashTableType
    PRIVATE
    LOGICAL :: HashTable_Finished !<Is .TRUE. if the table has finished being created, .FALSE. if not.

    INTEGER(INTG), DIMENSION(:), ALLOCATABLE :: T_key   ! The hash vector of keys
    INTEGER(INTG), DIMENSION(:), ALLOCATABLE :: T_val   ! The hash vector of values

    INTEGER(INTG), DIMENSION(:), ALLOCATABLE :: S_key   ! The original vector of keys
    INTEGER(INTG), DIMENSION(:), ALLOCATABLE :: S_val   ! The original vector of values

    INTEGER(INTG) :: N !<The number of items currently in the table (number of keys)
    INTEGER(INTG) :: P !<The prime number required from the table algorithm
  END TYPE HashTableType

  INTERFACE HashLinkedListType_Add
    MODULE PROCEDURE HashLinkedListType_Add_Data
    MODULE PROCEDURE HashLinkedListType_Add_List
  END INTERFACE HashLinkedListType_Add

  INTERFACE HashTable_Insert
    MODULE PROCEDURE HashTable_PutAll
    MODULE PROCEDURE HashTable_PutOne
  END INTERFACE HashTable_Insert

  ! public TYPEs
  PUBLIC :: HashTableType

  ! public subs
  PUBLIC :: HashTable_Insert, HashTable_Get, HashTable_CreateStart, HashTable_CreateFinish, HashTable_GetValue


CONTAINS

! -------------------------------------------------------------------

  !> initialises or adds a piece of data to list
  SUBROUTINE HashLinkedListType_Add_Data(list,key, value, is_duplicate_key, ERR)!,ERROR)

    TYPE(HashLinkedListType),INTENT(INOUT) :: list
    INTEGER(INTG),INTENT(IN) :: key, value
    INTEGER(INTG), INTENT(OUT) :: ERR
    LOGICAL, INTENT (OUT) :: is_duplicate_key
    ! CHARACTER(LEN=*), INTENT(OUT) :: ERROR
    ! local variables
    TYPE(HashLinkedListItemType),pointer :: current

    ! ENTERS("HashLinkedListType_Add_Data",ERR,ERROR,*999)

    is_duplicate_key = .FALSE.

    if (associated(list%root)) then
      current => list%root
       do
         ! key is there already
         IF (current%key==key) then
         is_duplicate_key = .TRUE.
         EXIT ! do nothing
         END if 

         ! tail is reached
         IF (.NOT. associated(current%next)) THEN
            allocate(current%next)
            current%next%key   = key
            current%next%value = value
            list%last => current%next
         EXIT
         END IF
          
         ! proceed
         current => current%next 

        END DO  
       
     ELSE     ! list is empty: create root 

      allocate(list%root)
      list%root%key   = key
      list%root%value = value
      list%last => list%root
    END IF

    !EXITS("HashLinkedListType_Add_Data")
    RETURN
!999 ERRORSEXITS("HashLinkedListType_Add_Data",ERR,ERROR)
!    RETURN 1

  END SUBROUTINE HashLinkedListType_Add_Data


  SUBROUTINE HashLinkedListType_Get_Data(list,key, value, ERR)
!,ERROR)

    TYPE(HashLinkedListType),INTENT(INOUT) :: list
    INTEGER(INTG),INTENT(IN)   :: key
    INTEGER(INTG),INTENT(OUT)  :: value
    INTEGER(INTG), INTENT(OUT) :: ERR

    ! local variables
    TYPE(HashLinkedListItemType),pointer :: current

    ERR = 1000 ! set to 0 if a value is found
    
    if (associated(list%root)) then
      current => list%root
       do
         ! key is found 
         IF (current%key==key)THEN
          value = current%value 
           
          ERR = 0

         END IF 

         ! tail is reached
         IF (.NOT. associated(current%next)) THEN
         ! value not found!!!
         EXIT
         END IF
          
         ! otherwise proceed
         current => current%next 

        END DO  
    END IF
   ! if list is empty: do nothing 
   
   END SUBROUTINE HashLinkedListType_Get_Data
! -------------------------------------------------------------------

  !> adds all data from one list (addlist) to another (list)
  Subroutine HashLinkedListType_Add_List(list,addlist,ERR,ERROR)
    TYPE(HashLinkedListType),intent(inout) :: list
    TYPE(HashLinkedListType),intent(in) :: addlist
    INTEGER(INTG), INTENT(OUT) :: ERR
    CHARACTER(LEN=*), INTENT(OUT) :: ERROR
    LOGICAL :: is_found
    ! local variables
    TYPE(HashLinkedListItemType),pointer :: current

    if (HashLinkedListType_is_Empty(addlist)) return

    current => addlist%root
    do
      call HashLinkedListType_Add_Data(list,current%key, current%value, is_found, ERR)
     !,ERROR)
      !,*999)
      if (associated(current%next)) then
        current => current%next
      else
        exit
      ENDif
    ENDdo

!    EXITS("HashLinkedListType_Add_List")
    RETURN
!999 ERRORSEXITS("HashLinkedListType_Add_List",ERR,ERROR)
!    RETURN 1

  END Subroutine HashLinkedListType_Add_List

! -------------------------------------------------------------------

  !> removes the first item from list and returns its value in data
  Subroutine HashLinkedListType_Remove_First(list,key,value,ERR,ERROR,*)
    TYPE(HashLinkedListType),intent(inout) :: list
    integer(INTG),intent(out) :: key, value
    INTEGER(INTG), INTENT(OUT) :: ERR
    CHARACTER(LEN=*), INTENT(OUT) :: ERROR
    ! local variables
    TYPE(HashLinkedListItemType),pointer :: next
    
    IF (associated(list%root)) THEN
      key   = list%root%key
      value = list%root%value
      next => list%root%next
      deallocate(list%root)
      list%root => next
      IF (associated(list%root)) THEN
        if (.not.associated(list%root%next)) list%last => list%root  ! only one left
      ELSE
        list%last => NULL()
      ENDIF
    ELSE
      write(*,*) ">>> warning: HLinked list is empty and cannot remove first item"
    ENDIF

  END Subroutine HashLinkedListType_Remove_First

! -------------------------------------------------------------------

  !> removes the last item from list and returns its value in data
  Subroutine HashLinkedListType_Remove_Last(list,key,value,ERR,ERROR,*)
    TYPE(HashLinkedListType),intent(inout) :: list
    integer(INTG),intent(out) :: key,value
    INTEGER(INTG), INTENT(OUT) :: ERR
    CHARACTER(LEN=*), INTENT(OUT) :: ERROR
    ! local variables
    TYPE(HashLinkedListItemType),pointer :: current

    if (.not.associated(list%root)) then
      write(*,*) ">>> warning: HLinked list is empty and cannot remove last item"
      return
    ENDif
    current => list%root

    do
      if (associated(current%next)) then
        if (associated(current%next%next)) then
          current => current%next
        else
          ! next one is the last one
          key = current%next%key
          value = current%next%value
          deallocate(current%next)
          current%next => NULL()
          list%last => current
          exit
        ENDif
      else
        ! there must be only one item in the list(?)!
        key   = current%key
        value = current%value
        deallocate(list%root)
        list%root => NULL()
        list%last => NULL()
        exit
      ENDif
    ENDdo

  END Subroutine HashLinkedListType_Remove_Last

! -------------------------------------------------------------------

  !> will delete and deallocate all items
  Subroutine HashLinkedListType_Destroy(list,ERR)
!,ERROR)

    TYPE(HashLinkedListType), INTENT(inout) :: list
    INTEGER(INTG), INTENT(OUT) :: ERR
 !   CHARACTER(LEN=*), INTENT(OUT) :: ERROR
    ! local variables
    TYPE(HashLinkedListItemType), POINTER :: current,next

    if (.not.associated(list%root)) return

    current => list%root
    do
      if (associated(current%next)) then
        next => current%next
        deallocate(current)
        current => next
      else
        deallocate(current)
        exit
      ENDif
    ENDdo
    
    list%root => NULL()
    list%last => NULL()

  END SUBROUTINE HashLinkedListType_Destroy

! -------------------------------------------------------------------

  !> returns true if the list is empty
  Function HashLinkedListType_is_Empty(list)
    TYPE(HashLinkedListType),intent(in) :: list
    logical :: HashLinkedListType_is_Empty

    HashLinkedListType_is_Empty = .true.
    if (associated(list%root)) HashLinkedListType_is_Empty = .false.
    
  END Function HashLinkedListType_is_Empty

! -------------------------------------------------------------------

  !> Returns length of list
  SUBROUTINE HashLinkedListType_Size(list,n,ERR)
    TYPE(HashLinkedListType),INTENT(IN) :: list
    INTEGER(INTG), INTENT(OUT) :: n,ERR
 
    ! local variables
    integer(INTG) :: i
    TYPE(HashLinkedListItemType),pointer :: current

    ! return zero-size array if list is empty
    IF (HashLinkedListType_is_Empty(list)) THEN
      n=0
      RETURN
    END IF

    ! Traversing to find size
    current => list%root
    n=1
    do
      if (associated(current%next)) then
        n=n+1
        current => current%next
      else
        exit
      END IF
    END DO

  END SUBROUTINE HashLinkedListType_Size

  !> copies out the data to an allocatable array of size 2xn
  Subroutine HashLinkedListType_to_Array(list,array,ERR,ERROR,*)
    TYPE(HashLinkedListType),intent(in) :: list
    integer(INTG),allocatable,intent(out) :: array(:,:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    CHARACTER(LEN=*), INTENT(OUT) :: ERROR
    ! local variables
    integer(INTG) :: i,n
    TYPE(HashLinkedListItemType),pointer :: current

    ! return zero-size array if list is empty
    if (HashLinkedListType_is_Empty(list)) then
      allocate(array(0,0))
      return
    ENDif

    ! first traversing to find size
    current => list%root
    n=1
    do
      if (associated(current%next)) then
        n=n+1
        current => current%next
      else
        exit
      ENDif
    ENDdo

    ! copy to array
    if (allocated(array)) deallocate(array)
    allocate(array(2,n),stat=err)
    !IF (ERR/=0) CALL ...
    current => list%root
    do i=1,n
      array(1,i)=current%key
      array(2,i)=current%value

      current => current%next
    ENDdo

  END Subroutine HashLinkedListType_to_Array

  !>Initialises a list and all its components
  SUBROUTINE HashListArrayType_Initialise(tbl,ERR,ERROR,*)
 
    TYPE (HashListArrayType), POINTER :: tbl
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR    

    ENTERS("HashListArrayType_Initialise",ERR,ERROR,*998)

    IF(ASSOCIATED(tbl)) THEN
      CALL FlagError("List is already associated.",ERR,ERROR,*998)
    ELSE
      ALLOCATE(tbl,STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate list.",ERR,ERROR,*999)
      tbl%is_finished=.FALSE.
      tbl%vec_len=10 ! default number of elements in array
    ENDIF

    EXITS("HashListArrayType_Initialise")
    RETURN
999 CALL HashListArrayType_Finalise(tbl,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("HashListArrayType_Initialise",ERR,ERROR)
    RETURN 1

  END SUBROUTINE HashListArrayType_Initialise

  !>Starts the creation of a list and returns a pointer to the created list
  SUBROUTINE HashListArrayType_CreateStart(tbl,ERR,ERROR,*)
 
    TYPE (HashListArrayType), POINTER :: tbl
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables

    ENTERS("HashListArrayType_CreateStart",ERR,ERROR,*999)


    CALL HashListArrayType_Initialise(tbl,ERR,ERROR,*999)
   
    tbl%is_init = .TRUE.

    EXITS("HashListArrayType_CreateStart")
    RETURN
999 ERRORSEXITS("HashListArrayType_CreateStart",ERR,ERROR)
    RETURN 1

  END SUBROUTINE HashListArrayType_CreateStart


  !>Finishes the creation of a table created with HashTable_CreateStart
  SUBROUTINE HashTable_CreateFinish(Table,ERR,ERROR,*)

    !Argument Variables
    TYPE(HashTableType), POINTER :: Table !<A pointer to the table to finish
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    
    TYPE(VARYING_STRING) :: dummyERROR !<The error string.
    INTEGER(INTG)        :: dummyERR

    ENTERS("HashTable_CreateFinish",ERR,ERROR,*998)

    IF(ASSOCIATED(Table)) THEN
      IF(Table%HashTable_Finished) THEN
        CALL FlagError("Table is already finished",ERR,ERROR,*998)
      ELSE
         ! Allocate pointers
         ! No pointers to allocate
         ! Finish the tree creation
        Table%HashTable_Finished=.TRUE.
      ENDIF
    ELSE
      CALL FlagError("Table is not associated",ERR,ERROR,*998)
    ENDIF

    EXITS("HashTable_CreateFinish")
    RETURN
999 CALL HashTable_Finalise(Table,dummyERR,dummyERROR,*998)
998 ERRORSEXITS("HashTable_CreateFinish",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HashTable_CreateFinish

  !
  !================================================================================================================================
  !

  !>Starts the creation of a hash table and returns a pointer to the created table.
  SUBROUTINE HashTable_CreateStart(Table,ERR,ERROR,*)

    !Argument Variables
    TYPE(HashTableType), POINTER :: Table !<A pointer to the table to create. Must NOT be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variable

    TYPE(VARYING_STRING) :: dummyERROR !<The error string.
    INTEGER(INTG)        :: dummyERR

    ENTERS("HashTable_CreateStart",ERR,ERROR,*998)

    IF(ASSOCIATED(Table)) THEN
      CALL FlagError("Table is already associated",ERR,ERROR,*998)
    ELSE
      ALLOCATE(Table,STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate table",ERR,ERROR,*999)
      CALL HashTable_Initialise(Table,ERR,ERROR,*999)
    ENDIF

    EXITS("HashTable_CreateStart")
    RETURN
999 CALL HashTable_Finalise(Table,dummyERR,dummyERROR,*998)
998 ERRORSEXITS("HashTable_CreateStart",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HashTable_CreateStart

 !>Initialises a tree
  SUBROUTINE HashTable_Initialise(Table,ERR,ERROR,*)

    !Argument Variables
    TYPE(HashTableType), POINTER :: Table !<A pointer to the table to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    
    ENTERS("HashTable_Initialise",ERR,ERROR,*999)

    IF(ASSOCIATED(Table)) THEN
      Table%HashTable_Finished=.FALSE.
      Table%n=0
      Table%p=0
    ELSE
      CALL FlagError("Table is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("HashTable_Initialise")
    RETURN
999 ERRORSEXITS("HashTable_Initialise",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HashTable_Initialise

  !
  !================================================================================================================================
  !

! queries for a value q and returns index
  SUBROUTINE HashTable_Get(Table, q, index_found, found, ERR, ERROR, *) 

    TYPE(HashTableType), POINTER :: Table
!   INTEGER(INTG), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: T_key ! the vector of keys
!   INTEGER(INTG), INTENT(IN)  :: q, n, p                         ! queried value, size, and prime number
    INTEGER(INTG), INTENT(IN)  :: q ! queried value
    INTEGER(INTG), INTENT(OUT) :: index_found
    LOGICAL, INTENT(OUT) :: found ! if object has been found

    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    INTEGER(INTG), INTENT(OUT) :: ERR
    
    INTEGER(INTG) :: j,block_ind

    ENTERS("HashTable_Get",ERR,ERROR,*999)

    found = .FALSE.
    index_found = 0 ! change??? returns k if used inappropriately!!!
    ! copy procedure as in paper
    j = my_alg_module(my_alg_module(Table%T_key(0)*q,Table%p),Table%n)

    block_ind = Table%T_key(j)
    IF (block_ind /= 0) THEN ! if block is non empty
       index_found = my_alg_module(my_alg_module(Table%T_key(block_ind+1)*q,Table%p),Table%T_key(block_ind)**2)+1+block_ind
       IF (Table%T_key(index_found)==q) found = .TRUE.
    END IF

    EXITS("HashTable_Get")
    RETURN
999 ERRORSEXITS("HashTable_Get",ERR,ERROR)
    RETURN 1

  END SUBROUTINE HashTable_Get

  ! Returns the value at index found
  SUBROUTINE HashTable_GetValue(Table, index_found, value, ERR,ERROR,*) 

    TYPE(HashTableType), POINTER :: Table    

    INTEGER(INTG), INTENT(IN)  :: index_found 
    INTEGER(INTG), INTENT(OUT) :: value

    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    INTEGER(INTG), INTENT(OUT) :: ERR


    ENTERS("HashTable_GetValue",ERR,ERROR,*999)
 
    IF (ASSOCIATED(Table)) THEN
       IF (ALLOCATED(Table%T_val)) THEN
         value = Table%T_val(index_found)
       ELSE
         CALL FlagError("Values are not there!",ERR,ERROR,*999)
       END IF
    ELSE
       CALL FlagError("Table is not associated",ERR,ERROR,*999)
    END IF
    
    EXITS("HashTable_GetValue")
    RETURN
999 ERRORSEXITS("HashTable_GetValue",ERR,ERROR)
    RETURN 1

  END SUBROUTINE HashTable_GetValue


  ! Build the table according to Fredman 1984
  SUBROUTINE HashTable_PutAll(Table,  S_key, S_val, ERR,ERROR,*)

    TYPE(HashTableType), POINTER :: Table    

    TYPE (HashListArrayType), POINTER  :: W_tbl ! W as a pointer to an array of lists
    INTEGER(INTG), INTENT(IN) :: S_key(:)       ! array of keys
    INTEGER(INTG), INTENT(IN) :: S_val(:)

    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    INTEGER(INTG), INTENT(OUT) :: ERR

    TYPE(VARYING_STRING) :: dummyERROR !<The error string.
    INTEGER(INTG)        :: dummyERR

    INTEGER(INTG) :: i,k,j, card_w, check_sum, check_sum_sq, count_fill, dim_t, prev_ind, index_w, card_w_store
    INTEGER(INTG), DIMENSION(:), ALLOCATABLE :: ind_in_T
    TYPE(HashLinkedListItemType) :: W_list_item

    LOGICAL :: is_found

    ENTERS("HashTable_PutAll",ERR,ERROR,*999)

    is_found = .FALSE.    

    ! Deallocate existing memory in table
   IF(ASSOCIATED(Table)) THEN
      IF(ALLOCATED(Table%T_key)) DEALLOCATE(Table%T_key)
      IF(ALLOCATED(Table%T_val)) DEALLOCATE(Table%T_val)
      IF(ALLOCATED(Table%S_key)) DEALLOCATE(Table%S_key)
      IF(ALLOCATED(Table%S_val)) DEALLOCATE(Table%S_val)
   END IF


    ! Store S_key and S_val for future additions
!    IF (PRESENT(S_key).AND.PRESENT(S_val)) THEN
      Table%n = size(S_key,1) 
      IF (Table%n /= size(S_val,1))  CALL FlagError("Keys and values have different size!",ERR,ERROR,*999)
      ALLOCATE(Table%S_key(Table%n), STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate S_key in table",ERR,ERROR,*999)
      ALLOCATE(Table%S_val(Table%n), STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate S_val in table",ERR,ERROR,*999)
      Table%S_key = S_key
      Table%S_val = S_val
 !   ELSE ! Make sure they are already there!
 !     IF (.NOT. ALLOCATED(Table%S_key)) CALL FlagError("Keys are not present!",ERR,ERROR,*999) 
 !     IF (.NOT. ALLOCATED(Table%S_val)) CALL FlagError("Values are not present!",ERR,ERROR,*999) 
 !     Table%n = size(Table%S_key,1)
 !   END IF

    ! Compute p as the first prime number larger than biggest key
    Table%p = maxval(Table%S_key,1)+1
    SELECT CASE(Table%p)
    CASE (1)
    CASE (2)
    CASE DEFAULT
     i=2
     DO !i = 2, Table%p-1
        IF (MOD(Table%p,i) == 0)  THEN
          Table%p = Table%p+1
          i = 2
          CYCLE
        END IF
        i = i+1
        IF (i == Table%p) EXIT
     END DO
    END SELECT 

    ! Initialise w as an array of size n
    ! HERE!!!
    CALL HTable_Init(W_tbl, ERR,ERROR,*999)
!    CALL HTable_Init(W_tbl, Table%n, ,ERR,ERROR,*999))
!    IF (ALLOCATED(tbl%vec)) DEALLOCATE(tbl%vec)
!    IF (PRESENT(tbl_len)) THEN
!      ALLOCATE(tbl%vec(0:tbl_len-1))
!       tbl%vec_len = tbl_len
!    ELSE
!       ALLOCATE(tbl%vec(0:tbl_size-1))
!       tbl%vec_len = tbl_size
!    END IF



    ! Compute the W_j
    k=2 ! use 2 as start as in paper (but 1 also works for paper case!!)

    DO

     DO i=1,Table%n

      index_w = my_alg_module(my_alg_module((k*Table%S_key(i)),Table%p),Table%n)
      !PRINT *, index_w
      CALL HashLinkedListType_Add_Data(W_tbl%vec(index_w-1),Table%S_key(i),Table%S_val(i), is_found, ERR) ! numbered from 0!
      IF (is_found) CALL FlagError("Duplicate key, error!",ERR,ERROR,*999) 
         !PRINT *, "Success with ", index_w
     END DO

     check_sum_sq =0
     check_sum   = 0
     count_fill  = 0
     DO i=1,Table%n
      CALL HashLinkedListType_Size(W_tbl%vec(i-1), card_w, ERR)
      IF (card_w /= 0) THEN
      check_sum_sq = check_sum_sq + card_w**2
      check_sum    = check_sum    + card_w
      count_fill   = count_fill +1
      END IF
     END DO

    IF (check_sum_sq<3*Table%n) THEN ! k is correct
     EXIT
    ELSE
     k = k+1 ! try another k
     IF (k>Table%p) CALL FlagError("Suitable k could not be found!",ERR,ERROR,*999) 
    END IF

   END DO

   ! Allocate T: T(0)=k, T(1),...,T(n) access indices, T(*) |W_j| T(**) k' T(***) = W_j
   dim_T = 1+Table%n+2*count_fill+check_sum_sq
   ALLOCATE (Table%T_key(0:(dim_T-1)),STAT=ERR)
   IF(ERR/=0) CALL FlagError("Could not allocate T_key in table.",ERR,ERROR,*999)
   ALLOCATE (Table%T_val(0:(dim_T-1)),STAT=ERR)
   IF(ERR/=0) CALL FlagError("Could not allocate T_val in table.",ERR,ERROR,*999)
   Table%T_key = 0
   Table%T_val = 0
   Table%T_key(0) = k
   IF (dim_T<=Table%n+1) CALL FlagError("No keys to be assigned!",ERR,ERROR,*999) 

   prev_ind = Table%n

   DO i=1,Table%n

     ! Cardinality of block W_i 
     CALL HashLinkedListType_Size(W_tbl%vec(i-1), card_w_store, ERR)

     IF (card_w_store == 0) THEN
       Table%T_key(i) = 0 ! Block is empty
     ELSE

       IF (prev_ind == Table%n) THEN
        Table%T_key(i) = Table%n+1 ! First non-zero entry is just index that follows n
       ELSE 
        ! Compute card of 
        ! CALL HashLinkedListType_Size(W_tbl%vec(prev_ind-1), card_w, ERR)    
        ! replaced with card_w = card_w_store below!!!!   
        Table%T_key(i) = Table%T_key(prev_ind)+2+card_w**2    
       END IF

       Table%T_key(Table%T_key(i))   = card_w_store ! First  element of each block: card of W
       Table%T_key(Table%T_key(i)+1) = 1            ! Second element of each block: k' (set to 1)

       ALLOCATE (ind_in_T(card_w_store), STAT=ERR)  ! ind_in_T contains the indices where keys will be stored in T for each block    
       IF(ERR/=0) CALL FlagError("Could not allocate ind_in_T in table.",ERR,ERROR,*999)
       W_list_item = W_tbl%vec(i-1)%root            ! take the first item in the block
       j = 1 ! init j for do loop

       DO 
         ! function given by corollary 2
         ! gives LOCAL index
         ! +2 translation
         ind_in_T(j) =  my_alg_module(my_alg_module(Table%T_key(Table%T_key(i)+1)*W_list_item%key,Table%p),card_w_store**2)+2
         IF ((j==1).OR. any(ind_in_T(1:(j-1))/=ind_in_T(j))) THEN ! insert items until conflict

            Table%T_key(ind_in_T(j)+Table%T_key(i)-1) = W_list_item%key 
            Table%T_val(ind_in_T(j)+Table%T_key(i)-1) = W_list_item%value

            IF (associated(W_list_item%next)) W_list_item = W_list_item%next
            
            IF (j==card_w_store) EXIT ! everything filled correctly
            j = j+1  ! next key
            
         ELSE ! conflict: increase k' and start over
             Table%T_key(Table%T_key(i)+1) =Table%T_key(Table%T_key(i)+1)+1 
             Table%T_key(ind_in_T + Table%T_key(i)-1) = 0   ! reset to zero what already filled (keys)
             Table%T_val(ind_in_T + Table%T_key(i)-1) = 0   ! reset to zero what already filled (values)
             ind_in_T = 0                       ! reset to zero indices
             j = 1                              ! reset j
             W_list_item = W_tbl%vec(i-1)%root  ! go back to list start
         END IF
 
       END DO
       
       !PRINT *, "IND IN T",  ind_in_T
       DEALLOCATE(ind_in_T)

       prev_ind = i
       card_w = card_w_store
     END IF

   END DO


   IF (ComputationalEnvironment_NodeNumberGet(ERR,ERROR)==0) THEN
   print *, "T vector n    ", Table%T_key (1:Table%n)
   print *, "T vector rest ", Table%T_key (Table%n+1:)
   print *, "T value vector rest ", Table%T_val (Table%n+1:)
   print *, "of size ", size(Table%T_key)
   END IF

   ! insert data in the list corresponding to hash index
   ! CALL HashLinkedListType_Add_Data(tbl%vec(hash),key,val, ERR)
    EXITS("HashTable_PutAll")
    RETURN
999 CALL HashTable_Finalise(Table,dummyERR,dummyERROR,*998)
998 ERRORSEXITS("HashTable_PutAll",ERR,ERROR)
    RETURN 1

  END SUBROUTINE HashTable_PutAll

  SUBROUTINE HashTable_PutOne(Table, S_key_one, S_val_one, ERR,ERROR,*)

    TYPE (HashTableType), POINTER :: Table
    INTEGER(INTG)        :: length
    INTEGER(INTG), ALLOCATABLE  :: tmp_key(:), tmp_val(:) 


    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    
    INTEGER(INTG), INTENT(IN) :: S_key_one, S_val_one

    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS("HashTable_PutOne",ERR,ERROR,*999)

    IF (.NOT. table%HashTable_Finished)  CALL FlagError("Table should be already initialised.",ERR,ERROR,*999) 
    IF (ALLOCATED(Table%T_key)) DEALLOCATE(Table%T_key)
    IF (ALLOCATED(Table%T_val)) DEALLOCATE(Table%T_val)

    IF(.NOT. ALLOCATED(Table%S_key)) CALL FlagError("Cannot put one new element in an empty table.",ERR,ERROR,*999)
    IF(.NOT. ALLOCATED(Table%S_val)) CALL FlagError("Cannot put one new element in an empty table.",ERR,ERROR,*999)

    length = size(Table%S_key)

    !Table%S_key(length+1) = S_key_one
    !Table%S_val(length+1) = S_val_one
    ! not allowed!!! Then:

    ALLOCATE(tmp_key(length+1), STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate tmp in table.",ERR,ERROR,*999)

    tmp_key(1:length) = Table%S_key
    tmp_key(length+1) = S_key_one
!    CALL MOVE_ALLOC(tmp, Table%S_key)

    ALLOCATE(tmp_val(length+1), STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate tmp in table.",ERR,ERROR,*999)

    tmp_val(1:length) = Table%S_val
    tmp_val(length+1) = S_val_one
!    CALL MOVE_ALLOC(tmp, Table%S_val)

    ! Recompute the table.
    CALL HashTable_PutAll(Table, tmp_key, tmp_val, ERR,ERROR,*999)

    EXITS("HashTable_PutOne")
    RETURN
999 CALL HashTable_Finalise(Table,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("HashTable_PutOne",ERR,ERROR)
    RETURN 1

  END SUBROUTINE HashTable_PutOne

  !>Finalises a hash table and deallocates all memory.
  SUBROUTINE HashTable_Finalise(Table,ERR,ERROR,*)    

    !Argument Variables
    TYPE (HashTableType), POINTER :: Table !<The table to finalise
    INTEGER(INTG), INTENT(OUT)    :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("HashTable_Finalise",ERR,ERROR,*999)
 
   IF(ASSOCIATED(Table)) THEN
      IF(ALLOCATED(Table%T_key)) DEALLOCATE(Table%T_key)
      IF(ALLOCATED(Table%T_val)) DEALLOCATE(Table%T_val)
      IF(ALLOCATED(Table%S_key)) DEALLOCATE(Table%S_key)
      IF(ALLOCATED(Table%S_val)) DEALLOCATE(Table%S_val)
      DEALLOCATE(Table)
   END IF

    EXITS("HashTable_Finalise")
    RETURN
999 ERRORSEXITS("HashTable_Finalise",ERR,ERROR)
    RETURN 1

  END SUBROUTINE HashTable_Finalise



FUNCTION my_alg_module(a,b) RESULT(c)

    INTEGER(INTG), INTENT(IN)  :: a,b
    INTEGER(INTG)  :: c

    c = MOD(a,b)
    IF (c==0) c = b

END FUNCTION my_alg_module


END MODULE hash_routines
