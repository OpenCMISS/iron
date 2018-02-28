!> \file
!> \author Chris Bradley
!> \brief This module contains all procedures for sorting. NOTE: THE ROUTINES IN THIS MODULE HAVE NOT BEEN TESTED!!!
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

!>This module contains all procedures for sorting. NOTE: THE ROUTINES IN THIS MODULE HAVE NOT BEEN TESTED!!!
MODULE Sorting

  USE BaseRoutines
  USE Constants
  USE Kinds
  USE ISO_VARYING_STRING

#include "macros.h"  
  
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Interfaces

  !>Sorts a list into assending order using the bubble sort method, returning the sorting index.
  INTERFACE Sorting_BubbleIndexSort
    MODULE PROCEDURE Sorting_BubbleIndexSortIntg
    MODULE PROCEDURE Sorting_BubbleIndexSortSP
    MODULE PROCEDURE Sorting_BubbleIndexSortDP
  END INTERFACE Sorting_BubbleIndexSort

  !>Sorts a list into assending order using the bubble sort method.
  INTERFACE Sorting_BubbleSort
    MODULE PROCEDURE Sorting_BubbleSortIntg
    MODULE PROCEDURE Sorting_BubbleSortSP
    MODULE PROCEDURE Sorting_BubbleSortDP
  END INTERFACE Sorting_BubbleSort

  !>Sorts a list into assending order using the heap sort method.
  INTERFACE Sorting_HeapSort
    MODULE PROCEDURE Sorting_HeapSortIntg
    MODULE PROCEDURE Sorting_HeapSortSP
    MODULE PROCEDURE Sorting_HeapSortDP
  END INTERFACE Sorting_HeapSort

  !>Sorts a list into assending order using the shell sort method.
  INTERFACE Sorting_ShellSort
    MODULE PROCEDURE Sorting_ShellSortIntg
    MODULE PROCEDURE Sorting_ShellSortSP
    MODULE PROCEDURE Sorting_ShellSortDP
  END INTERFACE Sorting_ShellSort

  PUBLIC Sorting_BubbleIndexSort

  PUBLIC Sorting_BubbleSort

  PUBLIC Sorting_HeapSort

  PUBLIC Sorting_ShellSort

CONTAINS

  !
  !================================================================================================================================
  !

  !>Sorts an integer list into assending order using the bubble sort method, returning the sorting index.
  SUBROUTINE Sorting_BubbleIndexSortIntg(a,indices,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: a(:) !<a(index). On entry the unsorted list of numbers. On exit, the sorted list of numbers in ascending order. 
    INTEGER(INTG), INTENT(OUT) :: indices(:) !<indices(index). On exit, the list of sorted indices      
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: flag,i,j,k,value,indexValue
    
    ENTERS("Sorting_BubbleIndexSortIntg",err,error,*999)

    IF(SIZE(indices,1)==SIZE(a,1)) THEN
      indices(1)=1  
      IF(SIZE(a,1)>1) THEN
        flag=SIZE(a,1)
        DO i=1,SIZE(a,1)
          k=flag-1
          flag=0
          DO j=1,k
            IF(i==1) indices(j+1)=j+1          
            IF(a(j)>a(j+1)) THEN
              value=a(j)
              a(j)=a(j+1)
              a(j+1)=value
              indexValue=indices(j)
              indices(j)=indices(j+1)
              indices(j+1)=indexValue              
              flag=j
            ENDIF
          ENDDO !j
          IF(flag==0) EXIT
        ENDDO !i
      ENDIF
    ELSE
      CALL FlagError("The size of the list to sort and indices arrays does not match.",err,error,*999)
    ENDIF      

    EXITS("Sorting_BubbleIndexSortIntg")
    RETURN
999 ERRORSEXITS("Sorting_BubbleIndexSortIntg",err,error)
    RETURN 1
    
  END SUBROUTINE Sorting_BubbleIndexSortIntg
  
  !
  !================================================================================================================================
  !
  !>Sorts a single precision list into assending order using the bubble sort method, returning the sorting index.
  SUBROUTINE Sorting_BubbleIndexSortSP(a,indices,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(INOUT) :: a(:) !<a(index). On entry the unsorted list of numbers. On exit, the sorted list of numbers in ascending order. 
    INTEGER(INTG), INTENT(OUT) :: indices(:) !<indices(index). On exit, the list of sorted indices      
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: flag,i,j,k,indexValue
    REAL(SP) :: value
    
    ENTERS("Sorting_BubbleIndexSortSP",err,error,*999)

    IF(SIZE(indices,1)==SIZE(a,1)) THEN
      indices(1)=1  
      IF(SIZE(a,1)>1) THEN
        flag=SIZE(a,1)
        DO i=1,SIZE(a,1)
          k=flag-1
          flag=0
          DO j=1,k
            IF(i==1) indices(j+1)=j+1          
            IF(a(j)>a(j+1)) THEN
              value=a(j)
              a(j)=a(j+1)
              a(j+1)=value
              indexValue=indices(j)
              indices(j)=indices(j+1)
              indices(j+1)=indexValue              
              flag=j
            ENDIF
          ENDDO !j
          IF(flag==0) EXIT
        ENDDO !i
      ENDIF
    ELSE
      CALL FlagError("The size of the list to sort and indices arrays does not match.",err,error,*999)
    ENDIF      

    EXITS("Sorting_BubbleIndexSortSP")
    RETURN
999 ERRORSEXITS("Sorting_BubbleIndexSortSP",err,error)
    RETURN 1
    
  END SUBROUTINE Sorting_BubbleIndexSortSP
  
  !
  !================================================================================================================================
  !
  !>Sorts a double precision list into assending order using the bubble sort method, returning the sorting index.
  SUBROUTINE Sorting_BubbleIndexSortDP(a,indices,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(INOUT) :: a(:) !<a(index). On entry the unsorted list of numbers. On exit, the sorted list of numbers in ascending order. 
    INTEGER(INTG), INTENT(OUT) :: indices(:) !<indices(index). On exit, the list of sorted indices      
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: flag,i,j,k,indexValue
    REAL(DP) :: VALUE
    
    ENTERS("Sorting_BubbleIndexSortDP",err,error,*999)

    IF(SIZE(indices,1)==SIZE(a,1)) THEN
      indices(1)=1  
      IF(SIZE(a,1)>1) THEN
        flag=SIZE(a,1)
        DO i=1,SIZE(a,1)
          k=flag-1
          flag=0
          DO j=1,k
            IF(i==1) indices(j+1)=j+1          
            IF(a(j)>a(j+1)) THEN
              value=a(j)
              a(j)=a(j+1)
              a(j+1)=value
              indexValue=indices(j)
              indices(j)=indices(j+1)
              indices(j+1)=indexValue              
              flag=j
            ENDIF
          ENDDO !j
          IF(flag==0) EXIT
        ENDDO !i
      ENDIF
    ELSE
      CALL FlagError("The size of the list to sort and indices arrays does not match.",err,error,*999)
    ENDIF      

    EXITS("Sorting_BubbleIndexSortDP")
    RETURN
999 ERRORSEXITS("Sorting_BubbleIndexSortDP",err,error)
    RETURN 1
    
  END SUBROUTINE Sorting_BubbleIndexSortDP
  
  !
  !================================================================================================================================
  !

  !>Sorts an integer list into assending order using the bubble sort method.
  SUBROUTINE Sorting_BubbleSortIntg(a,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: a(:) !<a(index). On entry the unsorted list of numbers. On exit, the sorted list of numbers in ascending order. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: flag,i,j,k,value
    
    ENTERS("Sorting_BubbleSortIntg",err,error,*999)

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

    EXITS("Sorting_BubbleSortIntg")
    RETURN
999 ERRORSEXITS("Sorting_BubbleSortIntg",err,error)
    RETURN 1
    
  END SUBROUTINE Sorting_BubbleSortIntg
  
  !
  !================================================================================================================================
  !

  !>Sorts a single precision list into assending order using the bubble sort method.
  SUBROUTINE Sorting_BubbleSortSP(a,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(INOUT) :: a(:) !<a(index). On entry the unsorted list of numbers. On exit, the sorted list of numbers in ascending order. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: flag,i,j,k
    REAL(SP) :: value
    
    ENTERS("Sorting_BubbleSortSP",err,error,*999)

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

    EXITS("Sorting_BubbleSortSP")
    RETURN
999 ERRORSEXITS("Sorting_BubbleSortSP",err,error)
    RETURN 1
    
  END SUBROUTINE Sorting_BubbleSortSP
  
  !
  !================================================================================================================================
  !

  !>Sorts a double precision list into assending order using the bubble sort method.
  SUBROUTINE Sorting_BubbleSortDP(a,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(INOUT) :: a(:) !<a(index). On entry the unsorted list of numbers. On exit, the sorted list of numbers in ascending order. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: flag,i,j,k
    REAL(DP) :: VALUE
    
    ENTERS("Sorting_BubbleSortDP",err,error,*999)

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

    EXITS("Sorting_BubbleSortDP")
    RETURN
999 ERRORSEXITS("Sorting_BubbleSortDP",err,error)
    RETURN 1
    
  END SUBROUTINE Sorting_BubbleSortDP
  
  !
  !================================================================================================================================
  !
  
  !>Sorts an integer list into assending order using the heap sort method.
  SUBROUTINE Sorting_HeapSortIntg(a,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: a(:) !<a(index). On entry the unsorted list of numbers. On exit, the sorted list of numbers in ascending order. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,indexValue,j,k,VALUE
    
    ENTERS("Sorting_HeapSortIntg",err,error,*999)

    IF(SIZE(a,1)>1) THEN      
      k=SIZE(a,1)/2+1
      indexValue=SIZE(a,1)
      DO 
        IF(k>1) THEN
          k=k-1
          value=a(k)
        ELSE
          value=a(indexValue)
          a(indexValue)=a(1)
          indexValue=indexValue-1
          IF(indexValue==1) THEN
            a(1)=value
            EXIT
          ENDIF
        ENDIF
        i=k
        j=k+k
        DO WHILE(j<=indexValue)
          IF(j<indexValue) THEN
            IF(a(j)<a(j+1)) j=j+1
          ENDIF
          IF(value<a(j)) THEN
            a(i)=a(j)
            i=j
            j=j+j
          ELSE
            j=indexValue+1
          ENDIF
        ENDDO
        a(i)=value
      ENDDO
    ENDIF

    EXITS("Sorting_HeapSortIntg")
    RETURN
999 ERRORSEXITS("Sorting_HeapSortIntg",err,error)
    RETURN 1
    
  END SUBROUTINE Sorting_HeapSortIntg
  
  !
  !================================================================================================================================
  !
  
  !>Sorts a single precision list into assending order using the heap sort method.
  SUBROUTINE Sorting_HeapSortSP(a,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(INOUT) :: a(:) !<a(index). On entry the unsorted list of numbers. On exit, the sorted list of numbers in ascending order. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,indexValue,j,k
    REAL(SP) :: value
    
    ENTERS("Sorting_HeapSortSP",err,error,*999)

    IF(SIZE(a,1)>1) THEN      
      k=SIZE(a,1)/2+1
      indexValue=SIZE(a,1)
      DO 
        IF(k>1) THEN
          k=k-1
          value=a(k)
        ELSE
          value=a(indexValue)
          a(indexValue)=a(1)
          indexValue=indexValue-1
          IF(indexValue==1) THEN
            a(1)=value
            EXIT
          ENDIF
        ENDIF
        i=k
        j=k+k
        DO WHILE(j<=indexValue)
          IF(j<indexValue) THEN
            IF(a(j)<a(j+1)) j=j+1
          ENDIF
          IF(value<a(j)) THEN
            a(i)=a(j)
            i=j
            j=j+j
          ELSE
            j=indexValue+1
          ENDIF
        ENDDO
        a(i)=value
      ENDDO
    ENDIF

    EXITS("Sorting_HeapSortSP")
    RETURN
999 ERRORSEXITS("Sorting_HeapSortSP",err,error)
    RETURN 1
    
  END SUBROUTINE Sorting_HeapSortSP
  
  !
  !================================================================================================================================
  !
  
  !>Sorts a double precision list into assending order using the heap sort method.
  SUBROUTINE Sorting_HeapSortDP(a,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(INOUT) :: a(:) !<a(index). On entry the unsorted list of numbers. On exit, the sorted list of numbers in ascending order. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,indexValue,j,k
    REAL(DP) :: value
    
    ENTERS("Sorting_HeapSortDP",err,error,*999)

    IF(SIZE(a,1)>1) THEN      
      k=SIZE(a,1)/2+1
      indexValue=SIZE(a,1)
      DO 
        IF(k>1) THEN
          k=k-1
          value=a(k)
        ELSE
          value=a(indexValue)
          a(indexValue)=a(1)
          indexValue=indexValue-1
          IF(indexValue==1) THEN
            a(1)=value
            EXIT
          ENDIF
        ENDIF
        i=k
        j=k+k
        DO WHILE(j<=indexValue)
          IF(j<indexValue) THEN
            IF(a(j)<a(j+1)) j=j+1
          ENDIF
          IF(value<a(j)) THEN
            a(i)=a(j)
            i=j
            j=j+j
          ELSE
            j=indexValue+1
          ENDIF
        ENDDO
        a(i)=value
      ENDDO
    ENDIF

    EXITS("Sorting_HeapSortDP")
    RETURN
999 ERRORSEXITS("Sorting_HeapSortDP",err,error)
    RETURN 1
    
  END SUBROUTINE Sorting_HeapSortDP
  
  !
  !================================================================================================================================
  !

  !>Sorts an integer list into assending order using the shell sort method.
  SUBROUTINE Sorting_ShellSortIntg(a,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: a(:) !<a(index). On entry the unsorted list of numbers. On exit, the sorted list of numbers in ascending order. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,increment,j,value
    
    ENTERS("Sorting_ShellSortIntg",err,error,*999)

    increment=4
    DO WHILE(increment<=SIZE(a,1))
      increment=3*increment+1
    ENDDO
    DO WHILE(increment>1)
      increment=increment/3
      DO i=increment+1,SIZE(a,1)
        value=a(i)
        j=i
        DO WHILE(a(j-increment)>value)
          a(j)=a(j-increment)
          j=j-increment
          IF(j<=increment) EXIT
        ENDDO
        a(j)=value
      ENDDO !i
    ENDDO

    EXITS("Sorting_ShellSortIntg")
    RETURN
999 ERRORSEXITS("Sorting_ShellSortIntg",err,error)
    RETURN 1
    
  END SUBROUTINE Sorting_ShellSortIntg
  
  !
  !================================================================================================================================
  !

  !>Sorts a single precision list into assending order using the shell sort method.
  SUBROUTINE Sorting_ShellSortSP(a,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(INOUT) :: a(:) !<a(index). On entry the unsorted list of numbers. On exit, the sorted list of numbers in ascending order. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,increment,j
    REAL(SP) :: value
    
    ENTERS("Sorting_ShellSortSP",err,error,*999)

    increment=4
    DO WHILE(increment<=SIZE(a,1))
      increment=3*increment+1
    ENDDO
    DO WHILE(increment>1)
      increment=increment/3
      DO i=increment+1,SIZE(a,1)
        value=a(i)
        j=i
        DO WHILE(a(j-increment)>value)
          a(j)=a(j-increment)
          j=j-increment
          IF(j<=increment) EXIT
        ENDDO
        a(j)=value
      ENDDO !i
    ENDDO

    EXITS("Sorting_ShellSortSP")
    RETURN
999 ERRORSEXITS("Sorting_ShellSortSP",err,error)
    RETURN 1
    
  END SUBROUTINE Sorting_ShellSortSP
  
  !
  !================================================================================================================================
  !

  !>Sorts a double precision list into assending order using the shell sort method.
  SUBROUTINE Sorting_ShellSortDP(a,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(INOUT) :: a(:) !<a(index). On entry the unsorted list of numbers. On exit, the sorted list of numbers in ascending order. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,increment,j
    REAL(DP) :: value
    
    ENTERS("Sorting_ShellSortDP",err,error,*999)

    increment=4
    DO WHILE(increment<=SIZE(a,1))
      increment=3*increment+1
    ENDDO
    DO WHILE(increment>1)
      increment=increment/3
      DO i=increment+1,SIZE(a,1)
        value=a(i)
        j=i
        DO WHILE(a(j-increment)>value)
          a(j)=a(j-increment)
          j=j-increment
          IF(j<=increment) EXIT
        ENDDO
        a(j)=value
      ENDDO !i
    ENDDO

    EXITS("Sorting_ShellSortDP")
    RETURN
999 ERRORSEXITS("Sorting_ShellSortDP",err,error)
    RETURN 1
    
  END SUBROUTINE Sorting_ShellSortDP
  
  !
  !================================================================================================================================
  !
  
END MODULE Sorting
