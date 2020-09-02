!> Only for integer data type for now
MODULE LinkedList_routines

#include "macros.h"  

  USE BaseRoutines
  USE Constants
  USE Kinds
  USE ISO_VARYING_STRING

  IMPLICIT NONE
  
  PRIVATE  ! by default

  ! types
  TYPE LinkedListItem
    INTEGER(INTG) :: data
    TYPE(LinkedListItem),pointer :: next => NULL()
  END TYPE

  TYPE LinkedList
    TYPE(LinkedListItem),pointer :: root => NULL()
    TYPE(LinkedListItem),pointer :: last => NULL()
  END TYPE

  INTERFACE LinkedList_Add
    MODULE PROCEDURE LinkedList_Add_Data
    MODULE PROCEDURE LinkedList_Add_List
  END INTERFACE

  ! public types
  PUBLIC :: LinkedListItem,LinkedList

  ! public subs
  PUBLIC :: LinkedList_Add,LinkedList_Destroy,LinkedList_Remove_First,LinkedList_Remove_Last
  PUBLIC :: LinkedList_is_Empty,LinkedList_to_Array

contains

! -------------------------------------------------------------------

  !> initialises or adds a piece of data to list
  SUBROUTINE LinkedList_Add_Data(list,data,err,error,*)

    TYPE(LinkedList),INTENT(INOUT) :: list
    INTEGER(INTG),INTENT(IN) :: data
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    ! local variables
    TYPE(LinkedListItem), POINTER :: current

    ENTERS("LinkedList_Add_Data",err,error,*999)

    if (associated(list%root)) then
      ! add to the tail end (for now)
      current => list%last
      allocate(current%next)
      current%next%data = data
      list%last => current%next
    else
      allocate(list%root)
      list%root%data = data
      list%last => list%root
    endif

    EXITS("LinkedList_Add_Data")
    RETURN
999 ERRORSEXITS("LinkedList_Add_Data",err,error)
    RETURN 1

  END SUBROUTINE LinkedList_Add_Data


! -------------------------------------------------------------------

  !> adds all data from one list to another
  SUBROUTINE LinkedList_Add_List(list,addlist,err,error,*)
    TYPE(LinkedList),intent(inout) :: list
    TYPE(LinkedList),intent(in) :: addlist
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    ! local variables
    TYPE(LinkedListItem), POINTER :: current

    ENTERS("LinkedList_Add_List",err,error,*999)
    
    if (LinkedList_is_Empty(addlist)) return

    current => addlist%root
    do
      call LinkedList_Add_Data(list,current%data,err,error,*999)
      if (associated(current%next)) then
        current => current%next
      else
        exit
      endif
    enddo

    EXITS("LinkedList_Add_List")
    RETURN
999 ERRORSEXITS("LinkedList_Add_List",err,error)
    RETURN 1

  END SUBROUTINE LinkedList_Add_List

! -------------------------------------------------------------------

  !> removes the first item from list and returns its value in data
  SUBROUTINE LinkedList_Remove_First(list,data,err,error,*)
    TYPE(LinkedList),intent(inout) :: list
    INTEGER(INTG),intent(out) :: data
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    ! local variables
    TYPE(LinkedListItem), POINTER :: next
    
    ENTERS("LinkedList_Remove_First",err,error,*999)
    
    if (associated(list%root)) then
      data = list%root%data
      next => list%root%next
      deallocate(list%root)
      list%root => next
      if (associated(list%root)) then
        if (.not.associated(list%root%next)) list%last => list%root  ! only one left
      else
        list%last => NULL()
      endif
    else
      CALL FlagWarning("Linked list is empty and cannot remove first item.",err,error,*999)
    endif

    EXITS("LinkedList_Remove_First")
    RETURN
999 ERRORSEXITS("LinkedList_Remove_First",err,error)
    RETURN 1
    
  END SUBROUTINE LinkedList_Remove_First

! -------------------------------------------------------------------

  !> removes the first item from list and returns its value in data
  SUBROUTINE LinkedList_Remove_Last(list,data,err,error,*)
    TYPE(LinkedList),intent(inout) :: list
    INTEGER(INTG),intent(out) :: data
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    ! local variables
    TYPE(LinkedListItem), POINTER :: current

    ENTERS("LinkedList_Remove_Last",err,error,*999)
    
    if (.not.associated(list%root)) then
      CALL FlagWarning("linked list is empty and cannot remove last item.",err,error,*999)
      return
    endif
    current => list%root

    do
      if (associated(current%next)) then
        if (associated(current%next%next)) then
          current => current%next
        else
          ! next one is the last one
          data = current%next%data
          deallocate(current%next)
          current%next => NULL()
          list%last => current
          exit
        endif
      else
        ! there must be only one item in the list?
        data = current%data
        deallocate(list%root)
        list%root => NULL()
        list%last => NULL()
        exit
      endif
    enddo

    EXITS("LinkedList_Remove_Last")
    RETURN
999 ERRORSEXITS("LinkedList_Remove_Last",err,error)
    RETURN 1
    
  END SUBROUTINE LinkedList_Remove_Last

! -------------------------------------------------------------------

  !> will delete and deallocate all items
  SUBROUTINE LinkedList_Destroy(list,err,error,*)

    TYPE(LinkedList), INTENT(inout) :: list
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    ! local variables
    TYPE(LinkedListItem), POINTER :: current,next

    ENTERS("LinkedList_Destroy",err,error,*999)
    
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
      endif
    enddo
    
    list%root => NULL()
    list%last => NULL()

    EXITS("LinkedList_Destroy")
    RETURN
999 ERRORSEXITS("LinkedList_Destroy",err,error)
    RETURN 1
    
   END SUBROUTINE LinkedList_Destroy

! -------------------------------------------------------------------

  !> returns true if the list is empty
  Function LinkedList_is_Empty(list)
    TYPE(LinkedList),intent(in) :: list
    logical :: LinkedList_is_Empty

    LinkedList_is_Empty = .true.
    if (associated(list%root)) LinkedList_is_Empty = .false.
    
  End Function LinkedList_is_Empty

! -------------------------------------------------------------------

  !> copies out the data to an allocatable array
  SUBROUTINE LinkedList_to_Array(list,array,err,error,*)
    TYPE(LinkedList),intent(in) :: list
    INTEGER(INTG),allocatable,intent(out) :: array(:)
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    ! local variables
    INTEGER(INTG) :: i,n
    TYPE(LinkedListItem), POINTER :: current

    ENTERS("LinkedList_to_Array",err,error,*999)
    
    ! return zero-size array if list is empty
    if (LinkedList_is_Empty(list)) then
      allocate(array(0))
      return
    endif

    ! first traversing to find size
    current => list%root
    n=1
    do
      if (associated(current%next)) then
        n=n+1
        current => current%next
      else
        exit
      endif
    enddo

    ! copy to array
    if (allocated(array)) deallocate(array)
    allocate(array(n),stat=err)
    !IF (ERR/=0) CALL ...
    current => list%root
    do i=1,n
      array(i)=current%data
      current => current%next
    enddo

    EXITS("LinkedList_to_Array")
    RETURN
999 ERRORSEXITS("LinkedList_to_Array",err,error)
    RETURN 1
  END SUBROUTINE LinkedList_to_Array

End Module LinkedList_routines
