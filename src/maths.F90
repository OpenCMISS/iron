!> \file
!> \author Chris Bradley
!> \brief This module contains all mathematics support routines.
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
!> Contributor(s): Chris Bradley, Kumar Mithraratne
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

!> This module contains all mathematics support routines.
MODULE Maths

  USE BaseRoutines
  USE Constants
  USE Kinds
  USE ISO_VARYING_STRING
  USE Strings

#include "macros.h"  
  
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Interfaces

  !>Returns hyperbolic cotangent of an argument
  INTERFACE Coth
    MODULE PROCEDURE CothSP
    MODULE PROCEDURE CothDP
  END INTERFACE Coth 

  !>Calculates the vector cross product of two vectors
  INTERFACE CrossProduct
    MODULE PROCEDURE CrossProductIntg
    MODULE PROCEDURE CrossProductSP
    MODULE PROCEDURE CrossProductDP
  END INTERFACE CrossProduct

  !>Returns hyperbolic cosecant of an argument
  INTERFACE Csch
    MODULE PROCEDURE CschSP
    MODULE PROCEDURE CschDP
  END INTERFACE Csch

  !>Calculates the the vector cross product of a x b in c and the n derivatives, dc, of the vector cross product given the derivatives da and db of a and b
  INTERFACE dCrossProduct
    MODULE PROCEDURE dCrossProductIntg
    MODULE PROCEDURE dCrossProductSP
    MODULE PROCEDURE dCrossProductDP
  END INTERFACE dCrossProduct

  !>Decomposes a matrix/tensor into its spherical and deviatoric parts
  INTERFACE DecomposeSphericalDeviatoric
    MODULE PROCEDURE DecomposeSphericalDeviatoricMatrix2SP
    MODULE PROCEDURE DecomposeSphericalDeviatoricMatrix2DP
  END INTERFACE DecomposeSphericalDeviatoric
  
  !>Decomposes a matrix/tensor into its symmetric and skew-symmetric parts
  INTERFACE DecomposeSymmetricSkew
    MODULE PROCEDURE DecomposeSymmetricSkewMatrix2SP
    MODULE PROCEDURE DecomposeSymmetricSkewMatrix2DP
  END INTERFACE DecomposeSymmetricSkew
  
  !>Returns the determinant of a matrix
  INTERFACE Determinant
    MODULE PROCEDURE DeterminantFullIntg
    MODULE PROCEDURE DeterminantFullSP
    MODULE PROCEDURE DeterminantFullDP
  END INTERFACE Determinant
        
  !>Calculates and returns the dot product between matrices and vectors i.e., C=A.B
  INTERFACE DotProduct
    MODULE PROCEDURE DotProductMatrix2Matrix2SP
    MODULE PROCEDURE DotProductMatrix2Matrix2DP
    MODULE PROCEDURE DotProductMatrixVectorSP
    MODULE PROCEDURE DotProductMatrixVectorDP
    MODULE PROCEDURE DotProductVectorVectorSP
    MODULE PROCEDURE DotProductVectorVectorDP
  END INTERFACE DotProduct
  
  !>Calculates and returns the dot product between a transposed matrices and vectors i.e., C=A.B^T
  INTERFACE DotProductTranspose
    MODULE PROCEDURE DotProductTransposeMatrix2Matrix2SP
    MODULE PROCEDURE DotProductTransposeMatrix2Matrix2DP
  END INTERFACE DotProductTranspose
  
  !>Calculates and returns the dot product between a transposed matrices and vectors i.e., C=A^T.B
  INTERFACE DotTransposeProduct
    MODULE PROCEDURE DotTransposeProductMatrix2Matrix2SP
    MODULE PROCEDURE DotTransposeProductMatrix2Matrix2DP
    MODULE PROCEDURE DotTransposeProductMatrixVectorSP
    MODULE PROCEDURE DotTransposeProductMatrixVectorDP
  END INTERFACE DotTransposeProduct
  
  !>Calculates and returns the double dot product between matrices and vectors i.e., C=A:B
  INTERFACE DoubleDotProduct
    MODULE PROCEDURE DoubleDotProductMatrix2Matrix2SP
    MODULE PROCEDURE DoubleDotProductMatrix2Matrix2DP
    MODULE PROCEDURE DoubleDotProductMatrix2Matrix4SP
    MODULE PROCEDURE DoubleDotProductMatrix2Matrix4DP
    MODULE PROCEDURE DoubleDotProductMatrix4Matrix2SP
    MODULE PROCEDURE DoubleDotProductMatrix4Matrix2DP
  END INTERFACE DoubleDotProduct

  !>Calculates the elliptic integral of the second kind - E(m).
  INTERFACE Edp
    MODULE PROCEDURE EdpSP
    MODULE PROCEDURE EdpDP
  END INTERFACE Edp

  !>Returns the eigenvalues of a matrix.
  INTERFACE Eigenvalue
    MODULE PROCEDURE EigenvalueFullSP
    MODULE PROCEDURE EigenvalueFullDP
  END INTERFACE Eigenvalue

  !>Returns the eigenvectors of a matrix.
  INTERFACE Eigenvector
    MODULE PROCEDURE EigenvectorFullSP
    MODULE PROCEDURE EigenvectorFullDP
  END INTERFACE Eigenvector

  !>Calculates the modified Bessel function of the first kind of order 0 using the approximation of Abromowitz and Stegun.
  INTERFACE I0
    MODULE PROCEDURE I0DP
    MODULE PROCEDURE I0SP
  END INTERFACE I0

  !>Calculates the modified Bessel function of the first kind of order 1 using the approximation of Abromowitz and Stegun.
  INTERFACE I1
    MODULE PROCEDURE I1DP
    MODULE PROCEDURE I1SP
  END INTERFACE I1

  !>Returns the invaiants of a matrix.
  INTERFACE Invariants
    MODULE PROCEDURE InvariantsOne2SP
    MODULE PROCEDURE InvariantsOne2DP
    MODULE PROCEDURE InvariantsTwo2SP
    MODULE PROCEDURE InvariantsTwo2DP
    MODULE PROCEDURE InvariantsThree2SP
    MODULE PROCEDURE InvariantsThree2DP
  END INTERFACE Invariants

  !>Returns the first invaiant of a matrix.
  INTERFACE InvariantsOne
    MODULE PROCEDURE InvariantsOne2SP
    MODULE PROCEDURE InvariantsOne2DP
  END INTERFACE InvariantsOne

  !>Returns the first three invaiants of a matrix.
  INTERFACE InvariantsThree
    MODULE PROCEDURE InvariantsThree2SP
    MODULE PROCEDURE InvariantsThree2DP
  END INTERFACE InvariantsThree

  !>Returns the first two invaiants of a matrix.
  INTERFACE InvariantsTwo
    MODULE PROCEDURE InvariantsTwo2SP
    MODULE PROCEDURE InvariantsTwo2DP
  END INTERFACE InvariantsTwo

  !>Returns the inverse of a matrix.
  INTERFACE Invert
    MODULE PROCEDURE InvertFullSP
    MODULE PROCEDURE InvertFullDP
  END INTERFACE Invert

  !>Returns the inverse of a transposed matrix.
  INTERFACE InvertTranspose
    MODULE PROCEDURE InvertTransposeFullSP
    MODULE PROCEDURE InvertTransposeFullDP
  END INTERFACE InvertTranspose

  !>Calculates the modified Bessel function of the second kind of order 0 using the approximation of Abromowitz and Stegun.
  INTERFACE K0
    MODULE PROCEDURE K0DP
    MODULE PROCEDURE K0SP
  END INTERFACE K0

  !>Calculates the modified Bessel function of the second kind of order 1 using the approximation of Abromowitz and Stegun.
  INTERFACE K1
    MODULE PROCEDURE K1DP
    MODULE PROCEDURE K1SP
  END INTERFACE K1

  !>Calculates the elliptic integral of the first kind - K(m).
  INTERFACE Kdp
    MODULE PROCEDURE KdpDP
    MODULE PROCEDURE KdpSP
  END INTERFACE Kdp

  !>Returns the identity matrix.
  INTERFACE IdentityMatrix
    MODULE PROCEDURE IdentityMatrixSP
    MODULE PROCEDURE IdentityMatrixDP
  END INTERFACE IdentityMatrix

  !>Returns the L1 norm of a matrix or vector.
  INTERFACE L1Norm
    MODULE PROCEDURE L1NormMatrixSP
    MODULE PROCEDURE L1NormMatrixDP
    MODULE PROCEDURE L1NormVectorSP
    MODULE PROCEDURE L1NormVectorDP
  END INTERFACE L1Norm
    
  !>Returns the L2 norm of a matrix or vector.
  INTERFACE L2Norm
    MODULE PROCEDURE L2NormMatrixSP
    MODULE PROCEDURE L2NormMatrixDP
    MODULE PROCEDURE L2NormVectorSP
    MODULE PROCEDURE L2NormVectorDP
  END INTERFACE L2Norm
    
  !>Returns the Linf norm of a matrix or vector.
  INTERFACE LInfNorm
    MODULE PROCEDURE LInfNormMatrixSP
    MODULE PROCEDURE LInfNormMatrixDP
    MODULE PROCEDURE LInfNormVectorSP
    MODULE PROCEDURE LInfNormVectorDP
  END INTERFACE LInfNorm
    
  !>Calculates the Macaulay Bracket of a number
  INTERFACE MacaulayBracket
    MODULE PROCEDURE MacaulayBracketSP
    MODULE PROCEDURE MacaulayBracketDP
  END INTERFACE MacaulayBracket
  
  !>Calculates and returns the matrix product between matrices and vectors i.e., C=A.B
  INTERFACE MatrixProduct
    MODULE PROCEDURE DotProductMatrix2Matrix2SP
    MODULE PROCEDURE DotProductMatrix2Matrix2DP
  END INTERFACE MatrixProduct
  
  !>Calculates and returns the matrix product between a matrix and a vector i.e., C=A.b
  INTERFACE MatrixVectorProduct
    MODULE PROCEDURE DotProductMatrixVectorSP
    MODULE PROCEDURE DotProductMatrixVectorDP
  END INTERFACE MatrixVectorProduct
  
  !>Calculates and returns the matrix product between a transposed matrices and vectors i.e., C=A.B^T
  INTERFACE MatrixProductTranspose
    MODULE PROCEDURE DotProductTransposeMatrix2Matrix2SP
    MODULE PROCEDURE DotProductTransposeMatrix2Matrix2DP
  END INTERFACE MatrixProductTranspose
  
  !>Calculates and returns the matrix product between a transposed matrices and vectors i.e., C=A^T.B
  INTERFACE MatrixTransposeProduct
    MODULE PROCEDURE DotTransposeProductMatrix2Matrix2SP
    MODULE PROCEDURE DotTransposeProductMatrix2Matrix2DP
  END INTERFACE MatrixTransposeProduct
  
  !>Calculates and returns the matrix product between a transposed matrix and a vector i.e., C=A^T.b
  INTERFACE MatrixTransposeVectorProduct
    MODULE PROCEDURE DotTransposeProductMatrixVectorSP
    MODULE PROCEDURE DotTransposeProductMatrixVectorDP
  END INTERFACE MatrixTransposeVectorProduct
  
  !>Normalises matrices and vectors
  INTERFACE Normalise
    MODULE PROCEDURE NormaliseMatrixSP
    MODULE PROCEDURE NormaliseMatrixDP
    MODULE PROCEDURE NormaliseVectorSP
    MODULE PROCEDURE NormaliseVectorDP
  END INTERFACE Normalise

  !>Calculates the normalised vector cross product of two vectors
  INTERFACE NormaliseCrossProduct
    MODULE PROCEDURE NormaliseCrossProductSP
    MODULE PROCEDURE NormaliseCrossProductDP
  END INTERFACE NormaliseCrossProduct

  !>Returns hyperbolic secant of an argument
  INTERFACE Sech
    MODULE PROCEDURE SechSP
    MODULE PROCEDURE SechDP
  END INTERFACE Sech

  !>Solves a small linear system Ax=b.
  INTERFACE SolveEigenproblem
    MODULE PROCEDURE SolveEigenproblemSP
    MODULE PROCEDURE SolveEigenproblemDP
  END INTERFACE SolveEigenproblem

  !>Solves a small linear system Ax=b.
  INTERFACE SolveSmallLinearSystem
    MODULE PROCEDURE SolveSmallLinearSystemSP
    MODULE PROCEDURE SolveSmallLinearSystemDP
  END INTERFACE SolveSmallLinearSystem

  !>Calculates and returns the tensor product between matrices and vectors.
  INTERFACE TensorProduct
    MODULE PROCEDURE TensorProductMatrix2Matrix2SP
    MODULE PROCEDURE TensorProductMatrix2Matrix2DP
    MODULE PROCEDURE TensorProductVectorVectorSP
    MODULE PROCEDURE TensorProductVectorVectorDP
  END INTERFACE TensorProduct
  
  !>Converts a tensor to a Voigt representation
  INTERFACE TensorToVoigt
    MODULE PROCEDURE TensorToVoigt2SP
    MODULE PROCEDURE TensorToVoigt2DP
    MODULE PROCEDURE TensorToVoigt4SP
    MODULE PROCEDURE TensorToVoigt4DP
  END INTERFACE TensorToVoigt

  !>Performs a transformation on a tensor in one coordinate system to the tensor in another coordinate system.
  INTERFACE TensorTransform
    MODULE PROCEDURE TensorTransformVectorDP0
    MODULE PROCEDURE TensorTransformVectorDP1
    MODULE PROCEDURE TensorTransformTensorTwoDP
    MODULE PROCEDURE TensorTransformTensorFourDP
  END INTERFACE TensorTransform
  
  !>Calculates and returns the trace of a matrix
  INTERFACE Trace
    MODULE PROCEDURE TraceSP
    MODULE PROCEDURE TraceDP
    MODULE PROCEDURE TraceIntg
  END INTERFACE Trace
  
  !>Returns the transpose of a matrix A in A^T.
  INTERFACE MatrixTranspose
    MODULE PROCEDURE TransposeMatrixSP
    MODULE PROCEDURE TransposeMatrixDP
  END INTERFACE MatrixTranspose

  !>Calculates and returns the unimodular form of a matrix/tensor
  INTERFACE Unimodular
    MODULE PROCEDURE UnimodularMatrixSP
    MODULE PROCEDURE UnimodularMatrixDP
  END INTERFACE Unimodular
  
  !>Converts a Voigt representation to a tensor
  INTERFACE VoigtToTensor
    MODULE PROCEDURE VoigtToTensor2SP0
    MODULE PROCEDURE VoigtToTensor2SP1
    MODULE PROCEDURE VoigtToTensor2DP0
    MODULE PROCEDURE VoigtToTensor2DP1
    MODULE PROCEDURE VoigtToTensor4SP
    MODULE PROCEDURE VoigtToTensor4DP
  END INTERFACE VoigtToTensor

  !>Zeros a matrix
  INTERFACE ZeroMatrix
    MODULE PROCEDURE ZeroMatrixSP
    MODULE PROCEDURE ZeroMatrixDP
  END INTERFACE ZeroMatrix

  PUBLIC Coth

  PUBLIC CrossProduct

  PUBLIC Csch
  
  PUBLIC dCrossProduct

  PUBLIC DecomposeSphericalDeviatoric

  PUBLIC DecomposeSymmetricSkew

  PUBLIC DotProduct

  PUBLIC DotProductTranspose

  PUBLIC DotTransposeProduct

  PUBLIC DoubleDotProduct

  PUBLIC Determinant

  PUBLIC Eigenvalue,Eigenvector

  PUBLIC I0,I1

  PUBLIC IdentityMatrix

  PUBLIC Invariants

  PUBLIC InvariantsOne

  PUBLIC InvariantsTwo

  PUBLIC InvariantsThree

  PUBLIC Invert

  PUBLIC InvertTranspose

  PUBLIC K0,K1,Kdp

  PUBLIC L1Norm

  PUBLIC L2Norm

  PUBLIC LInfNorm

  PUBLIC MacaulayBracket

  PUBLIC MatrixProduct

  PUBLIC MatrixProductTranspose

  PUBLIC MatrixTranspose

  PUBLIC MatrixTransposeProduct

  PUBLIC MatrixTransposeVectorProduct

  PUBLIC MatrixVectorProduct

  PUBLIC Normalise

  PUBLIC NormaliseCrossProduct

  PUBLIC Sech

  PUBLIC SolveEigenproblem

  PUBLIC SolveSmallLinearSystem

  PUBLIC spline_cubic_set,s3_fs,spline_cubic_val

  PUBLIC TensorProduct

  PUBLIC TensorToVoigt

  PUBLIC TensorTransform

  PUBLIC Trace

  PUBLIC Unimodular

  PUBLIC VoigtToTensor

  PUBLIC ZeroMatrix
  
CONTAINS

  !
  !================================================================================================================================
  !
  
  !>Calculates single precision hyperbolic cotangent function
  PURE FUNCTION CothSP(a)

    !Argument variables
    REAL(SP), INTENT(IN) :: a !<argument to perform coth() on
    !Function variable
    REAL(SP) :: CothSP
    
    CothSP=COSH(a)/SINH(a)

    RETURN
    
  END FUNCTION CothSP

  !
  !================================================================================================================================
  !

  !> Calculates double precision hyperbolic cotangent function
  PURE FUNCTION CothDP(a)

    !Argument variables
    REAL(DP), INTENT(IN) :: a !<argument to perform coth() on
    !Function variable
    REAL(DP) :: CothDP

    CothDP=COSH(a)/SINH(a)

    RETURN
    
  END FUNCTION CothDP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the vector cross-product of the integer vectors a x b in c.
  SUBROUTINE CrossProductIntg(a,b,c,err,error,*)
      
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: a(:) !<The first vector in the cross product
    INTEGER(INTG), INTENT(IN) :: b(:) !<The second vector in the cross product
    INTEGER(INTG), INTENT(OUT) :: c(:) !<On exit, the cross product of the first and second vectors
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    
    ENTERS("CrossProductIntg",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(a,1)/=SIZE(b,1)) CALL FlagError("The vectors a and b are not the same size.",err,error,*999)
    IF(SIZE(c,1)==3) CALL FlagError("The vector c is not the correct size.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(a,1))
    CASE(3)
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
    CASE DEFAULT
      CALL FlagError("Invalid vector size.",err,error,*999)
    END SELECT
 
    EXITS("CrossProductIntg")
    RETURN
999 ERRORSEXITS("CrossProductIntg",err,error)
    RETURN 1
    
  END SUBROUTINE CrossProductIntg
  
  !
  !================================================================================================================================
  !

  !>Calculates and returns the vector cross-product of the single precision vectors axb in c.
  SUBROUTINE CrossProductSP(a,b,c,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: a(:) !<The first vector in the cross product
    REAL(SP), INTENT(IN) :: b(:) !<The second vector in the cross product
    REAL(SP), INTENT(OUT) :: c(:) !<On exit, the cross product of the first and second vectors
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    
    ENTERS("CrossProductSP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(a,1)/=SIZE(b,1)) CALL FlagError("The vectors a and b are not the same size.",err,error,*999)
    IF(SIZE(c,1)/=3) CALL FlagError("The vector c is not the correct size.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1))
    CASE(3)
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
    CASE DEFAULT
      CALL FlagError("Invalid vector size.",err,error,*999)
    END SELECT

    EXITS("CrossProductSP")
    RETURN
999 ERRORSEXITS("CrossProductSP",err,error)
    RETURN 1
    
  END SUBROUTINE CrossProductSP
  
  !
  !================================================================================================================================
  !

  !>Calculates and returns the vector cross-product of the double precision vectors axb in c.
  SUBROUTINE CrossProductDP(a,b,c,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: a(:) !<The first vector in the cross product
    REAL(DP), INTENT(IN) :: b(:) !<The second vector in the cross product
    REAL(DP), INTENT(OUT) :: c(:) !<On exit, the cross product of the first and second vectors
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    
    ENTERS("CrossProductDP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(a,1)/=SIZE(b,1)) CALL FlagError("The vectors a and b are not the same size.",err,error,*999)
    IF(SIZE(c,1)/=3) CALL FlagError("The vector c is not the correct size.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(a,1))
    CASE(3)
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
    CASE DEFAULT
      CALL FlagError("Invalid vector size.",err,error,*999)
    END SELECT

    EXITS("CrossProductDP")
    RETURN
999 ERRORSEXITS("CrossProductDP",err,error)
    RETURN 1
    
  END SUBROUTINE CrossProductDP
  
  !
  !================================================================================================================================
  !
  
  !>Calculates single precision hyperbolic cosecant function
  PURE FUNCTION CschSP(a)

    !Argument variables
    REAL(SP), INTENT(IN) :: a !<argument to perform csch() on
    !Function variable
    REAL(SP) :: CschSP
    
    CschSP=1.0/SINH(a)

    RETURN
    
  END FUNCTION CschSP

  !
  !================================================================================================================================
  !

  !>Calculates double precision hyperbolic cosecant function
  PURE FUNCTION CschDP(a)

    !Argument variables
    REAL(DP), INTENT(IN) :: a !<argument to perform csch() on
    !Function variable
    REAL(DP) :: CschDP

    CschDP=1.0_DP/SINH(a)

    RETURN
    
  END FUNCTION CschDP

  !
  !================================================================================================================================
  !
  
  !>Calculates the the vector cross product of axb in c and the n derivatives, dc, of the vector cross product given the 
  !>derivatives da and db of a and b for integer vectors.
  SUBROUTINE dCrossProductIntg(n,a,b,c,da,db,dc,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: n !<The number of derivatives
    INTEGER(INTG), INTENT(IN) :: a(:) !<The a vector
    INTEGER(INTG), INTENT(IN) :: b(:) !<The b vector
    INTEGER(INTG), INTENT(OUT) :: c(:) !<On exit, the cross product of a x b
    INTEGER(INTG), INTENT(IN) :: da(:,:) !<The n derivatives of a
    INTEGER(INTG), INTENT(IN) :: db(:,:) !<The n derivatives of b
    INTEGER(INTG), INTENT(OUT) :: dc(:,:) !<On exit, the derivatives of c
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: derivIdx
    
    ENTERS("dCrossProductIntg",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(da,1)/=SIZE(db,1).OR.SIZE(a,1)/=SIZE(da,1).OR.SIZE(b,1)/=SIZE(db,1)) &
      & CALL FlagError("The vectors for da and db are not the same size.",err,error,*999)
    IF(SIZE(da,2)<n.OR.SIZE(db,2)<n) CALL FlagError("The number of derivative vectors is too small.",err,error,*999)
    IF(SIZE(c,1)/=3) CALL FlagError("The vector c is not the correct size.",err,error,*999)
#endif    
    
    CALL CrossProduct(a,b,c,err,error,*999)
    SELECT CASE(SIZE(a,1))
    CASE(3)
      DO derivIdx=1,n
        dc(1,derivIdx)=da(2,derivIdx)*b(3)-da(3,derivIdx)*b(2)+a(2)*db(3,derivIdx)-a(3)*db(2,derivIdx)
        dc(2,derivIdx)=da(3,derivIdx)*b(1)-da(1,derivIdx)*b(3)+a(3)*db(1,derivIdx)-a(1)*db(3,derivIdx)
        dc(3,derivIdx)=da(1,derivIdx)*b(2)-da(2,derivIdx)*b(1)+a(1)*db(2,derivIdx)-a(2)*db(1,derivIdx)
      ENDDO !derivIdx
    CASE DEFAULT
      CALL FlagError("Invalid vector size.",err,error,*999)
    END SELECT

    EXITS("dCrossProductIntg")
    RETURN
999 ERRORSEXITS("dCrossProductIntg",err,error)
    RETURN 1
    
  END SUBROUTINE dCrossProductIntg
  
  !
  !================================================================================================================================
  !
  
  !>Calculates the the vector cross product of axb in c and the n derivatives, dc, of the vector cross product given the 
  !>derivatives da and db of a and b for single precision vectors.
  SUBROUTINE dCrossProductSP(n,a,b,c,da,db,dc,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: n !<The number of derivatives
    REAL(SP), INTENT(IN) :: a(:) !<The a vector
    REAL(SP), INTENT(IN) :: b(:) !<The b vector
    REAL(SP), INTENT(OUT) :: c(:) !<On exit, the cross product of a x b
    REAL(SP), INTENT(IN) :: da(:,:) !<The n derivatives of a
    REAL(SP), INTENT(IN) :: db(:,:) !<The n derivatives of b
    REAL(SP), INTENT(OUT) :: dc(:,:) !<On exit, the derivatives of c
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: derivIdx
    
    ENTERS("dCrossProductSP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(da,1)/=SIZE(db,1).OR.SIZE(a,1)/=SIZE(da,1).OR.SIZE(b,1)/=SIZE(db,1)) &
      & CALL FlagError("The vectors for da and db are not the same size.",err,error,*999)
    IF(SIZE(da,2)<n.OR.SIZE(db,2)<n) CALL FlagError("The number of derivative vectors is too small.",err,error,*999)
    IF(SIZE(c,1)/=3) CALL FlagError("The vector c is not the correct size.",err,error,*999)
#endif    
    
    CALL CrossProduct(a,b,c,err,error,*999)
    SELECT CASE(SIZE(a,1))
    CASE(3)
      DO derivIdx=1,n
        dc(1,derivIdx)=da(2,derivIdx)*b(3)-da(3,derivIdx)*b(2)+a(2)*db(3,derivIdx)-a(3)*db(2,derivIdx)
        dc(2,derivIdx)=da(3,derivIdx)*b(1)-da(1,derivIdx)*b(3)+a(3)*db(1,derivIdx)-a(1)*db(3,derivIdx)
        dc(3,derivIdx)=da(1,derivIdx)*b(2)-da(2,derivIdx)*b(1)+a(1)*db(2,derivIdx)-a(2)*db(1,derivIdx)
      ENDDO !derivIdx
    CASE DEFAULT
      CALL FlagError("Invalid vector size.",err,error,*999)
    END SELECT

    EXITS("dCrossProductSP")
    RETURN
999 ERRORSEXITS("dCrossProductSP",err,error)
    RETURN 1
    
  END SUBROUTINE dCrossProductSP
  
  !
  !================================================================================================================================
  !

  !>Calculates the the vector cross product of axb in c and the n derivatives, dc, of the vector cross product given the 
  !>derivatives da and db of a and b for double precision vectors.
  SUBROUTINE dCrossProductDP(n,a,b,c,da,db,dc,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: n !<The number of derivatives
    REAL(DP), INTENT(IN) :: a(:) !<The a vector
    REAL(DP), INTENT(IN) :: b(:) !<The b vector
    REAL(DP), INTENT(OUT) :: c(:) !<On exit, the cross product of a x b
    REAL(DP), INTENT(IN) :: da(:,:) !<The n derivatives of a
    REAL(DP), INTENT(IN) :: db(:,:) !<The n derivatives of b
    REAL(DP), INTENT(OUT) :: dc(:,:) !<On exit, the derivatives of c
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: derivIdx
    
    ENTERS("dCrossProductDP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(da,1)/=SIZE(db,1).OR.SIZE(a,1)/=SIZE(da,1).OR.SIZE(b,1)/=SIZE(db,1)) &
      & CALL FlagError("The vectors for da and db are not the same size.",err,error,*999)
    IF(SIZE(da,2)<n.OR.SIZE(db,2)<n) CALL FlagError("The number of derivative vectors is too small.",err,error,*999)
    IF(SIZE(c,1)/=3) CALL FlagError("The vector c is not the correct size.",err,error,*999)
#endif    

    CALL CrossProduct(a,b,c,err,error,*999)
    SELECT CASE(SIZE(a,1))
    CASE(3)
      DO derivIdx=1,n
        dc(1,derivIdx)=da(2,derivIdx)*b(3)-da(3,derivIdx)*b(2)+a(2)*db(3,derivIdx)-a(3)*db(2,derivIdx)
        dc(2,derivIdx)=da(3,derivIdx)*b(1)-da(1,derivIdx)*b(3)+a(3)*db(1,derivIdx)-a(1)*db(3,derivIdx)
        dc(3,derivIdx)=da(1,derivIdx)*b(2)-da(2,derivIdx)*b(1)+a(1)*db(2,derivIdx)-a(2)*db(1,derivIdx)
      ENDDO !derivIdx
    CASE DEFAULT
      CALL FlagError("Invalid vector size.",err,error,*999)
    END SELECT

    EXITS("dCrossProductDP")
    RETURN
999 ERRORSEXITS("dCrossProductDP",err,error)
    RETURN 1
    
  END SUBROUTINE dCrossProductDP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the spherical and deviatoric parts of a single precision matrix2 matrix A i.e., sphA in B and devA in C.
  SUBROUTINE DecomposeSphericalDeviatoricMatrix2SP(A,B,C,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The matrix A to decompose
    REAL(SP), INTENT(OUT) :: B(:,:) !<On exit, the spherical part of A i.e., B=sphA
    REAL(SP), INTENT(OUT) :: C(:,:) !<On exit, the deviatoric part of A i.e., C=devA
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j
    REAL(SP) :: factor
    
    ENTERS("DecomposeSphericalDeviatoricMatrix2SP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(B,1).OR.SIZE(A,1)/=SIZE(C,1).OR.SIZE(A,2)/=SIZE(B,2).OR.SIZE(A,2)/=SIZE(C,2)) &
      & CALL FlagError("The sizes of the A, B & C matrices are not the same.",err,error,*999)
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("The A matrix is not square.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1))
    CASE(1)
      !Spherical part
      B(1,1)=A(1,1)
      !Deviatoric part
      C(1,1)=0.0_SP
    CASE(2)
      factor=(A(1,1)+A(2,2))/2.0_SP
      !Spherical part
      B(1,1)=factor
      B(1,2)=0.0_SP
      B(2,1)=0.0_SP
      B(2,2)=factor
      !Deviatoric part
      C(1,1)=A(1,1)-factor
      C(1,2)=A(1,2)
      C(2,1)=A(2,1)
      C(2,2)=A(2,2)-factor
    CASE(3)
      factor=(A(1,1)+A(2,2)+A(3,3))/3.0_SP
      !Spherical part
      B(1,1)=factor
      B(1,2)=0.0_SP
      B(1,3)=0.0_SP
      B(2,1)=0.0_SP
      B(2,2)=factor
      B(2,3)=0.0_SP
      B(3,1)=0.0_SP
      B(3,2)=0.0_SP
      B(3,3)=factor
      !Deviatoric part
      C(1,1)=A(1,1)-factor
      C(1,2)=A(1,2)
      C(1,3)=A(1,3)
      C(2,1)=A(2,1)
      C(2,2)=A(2,2)-factor
      C(2,3)=A(2,3)
      C(3,1)=A(3,1)
      C(3,2)=A(3,2)
      C(3,3)=A(3,3)-factor
    CASE DEFAULT
      CALL Trace(A,factor,err,error,*999)
      factor=factor/SIZE(A,1)
      DO i=1,SIZE(A,1)
        DO j=1,SIZE(A,1)
          B(i,j) = 0.0_SP
          C(i,j) = A(i,j)
        ENDDO !j
        B(i,i) = factor
        C(i,i) = C(i,i)-factor
      ENDDO !i
    END SELECT
 
    EXITS("DecomposeSphericalDeviatoricMatrix2SP")
    RETURN
999 ERRORSEXITS("DecomposeSphericalDeviatoricMatrix2SP",err,error)
    RETURN 1
    
  END SUBROUTINE DecomposeSphericalDeviatoricMatrix2SP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the spherical and deviatoric parts of a double precision matrix2 matrix A i.e., sphA in B and devA in C.
  SUBROUTINE DecomposeSphericalDeviatoricMatrix2DP(A,B,C,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The matrix A to decompose
    REAL(DP), INTENT(OUT) :: B(:,:) !<On exit, the spherical part of A i.e., B=sphA
    REAL(DP), INTENT(OUT) :: C(:,:) !<On exit, the deviatoric part of A i.e., C=devA
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j
    REAL(DP) :: factor
    
    ENTERS("DecomposeSphericalDeviatoricMatrix2DP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(B,1).OR.SIZE(A,1)/=SIZE(C,1).OR.SIZE(A,2)/=SIZE(B,2).OR.SIZE(A,2)/=SIZE(C,2)) &
      & CALL FlagError("The sizes of the A, B & C matrices are not the same.",err,error,*999)
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("The A matrix is not square.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1))
    CASE(1)
      !Spherical part
      B(1,1)=A(1,1)
      !Deviatoric part
      C(1,1)=0.0_DP
    CASE(2)
      factor=(A(1,1)+A(2,2))/2.0_DP
      !Spherical part
      B(1,1)=factor
      B(1,2)=0.0_DP
      B(2,1)=0.0_DP
      B(2,2)=factor
      !Deviatoric part
      C(1,1)=A(1,1)-factor
      C(1,2)=A(1,2)
      C(2,1)=A(2,1)
      C(2,2)=A(2,2)-factor
    CASE(3)
      factor=(A(1,1)+A(2,2)+A(3,3))/3.0_DP
      !Spherical part
      B(1,1)=factor
      B(1,2)=0.0_DP
      B(1,3)=0.0_DP
      B(2,1)=0.0_DP
      B(2,2)=factor
      B(2,3)=0.0_DP
      B(3,1)=0.0_DP
      B(3,2)=0.0_DP
      B(3,3)=factor
      !Deviatoric part
      C(1,1)=A(1,1)-factor
      C(1,2)=A(1,2)
      C(1,3)=A(1,3)
      C(2,1)=A(2,1)
      C(2,2)=A(2,2)-factor
      C(2,3)=A(2,3)
      C(3,1)=A(3,1)
      C(3,2)=A(3,2)
      C(3,3)=A(3,3)-factor
    CASE DEFAULT
      CALL Trace(A,factor,err,error,*999)
      factor=factor/SIZE(A,1)
      DO i=1,SIZE(A,1)
        DO j=1,SIZE(A,1)
          B(i,j) = 0.0_DP
          C(i,j) = A(i,j)
        ENDDO !j
        B(i,i) = factor
        C(i,i) = C(i,i)-factor
      ENDDO !i
    END SELECT
 
    EXITS("DecomposeSphericalDeviatoricMatrix2DP")
    RETURN
999 ERRORSEXITS("DecomposeSphericalDeviatoricMatrix2DP",err,error)
    RETURN 1
    
  END SUBROUTINE DecomposeSphericalDeviatoricMatrix2DP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the symmetric and skew-symmetric parts of a single precision matrix2 matrix A i.e., symA in B and skwA in C.
  SUBROUTINE DecomposeSymmetricSkewMatrix2SP(A,B,C,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The matrix A to decompose
    REAL(SP), INTENT(OUT) :: B(:,:) !<On exit, the symmetric part of A i.e., B=symA
    REAL(SP), INTENT(OUT) :: C(:,:) !<On exit, the skew-symmetric part of A i.e., C=skwA
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j
    
    ENTERS("DecomposeSymmetricSkewMatrix2SP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(B,1).OR.SIZE(A,1)/=SIZE(C,1).OR.SIZE(A,2)/=SIZE(B,2).OR.SIZE(A,2)/=SIZE(C,2)) &
      & CALL FlagError("The sizes of the A, B & C matrices are not the same.",err,error,*999)
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("The A matrix is not square.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1))
    CASE(1)
      !Symmetric part
      B(1,1)=A(1,1)
      !Skew-symmetric part
      C(1,1)=0.0_SP
    CASE(2)
      !Symmetric part
      B(1,1)=A(1,1)
      B(1,2)=(A(1,2)+A(2,1))/2.0_SP
      B(2,1)=B(1,2)
      B(2,2)=A(2,2)
      !Skew-symmetric part
      C(1,1)=0.0_SP
      C(1,2)=(A(1,2)-A(2,1))/2.0_SP
      C(2,1)=(A(2,1)-A(1,2))/2.0_SP
      C(2,2)=0.0_SP
    CASE(3)
      !Symmetric part
      B(1,1)=A(1,1)
      B(1,2)=(A(1,2)+A(2,1))/2.0_SP
      B(1,3)=(A(1,3)+A(3,1))/2.0_SP
      B(2,1)=B(1,2)
      B(2,2)=A(2,2)
      B(2,3)=(A(2,3)+A(3,2))/2.0_SP
      B(3,1)=B(1,3)
      B(3,2)=B(2,3)
      B(3,3)=A(3,3)
      !Skew-symmetric part
      C(1,1)=0.0_SP
      C(1,2)=(A(1,2)-A(2,1))/2.0_SP
      C(1,3)=(A(1,3)-A(3,1))/2.0_SP
      C(2,1)=(A(2,1)-A(2,1))/2.0_SP
      C(2,2)=0.0_SP
      C(2,3)=(A(2,3)-A(3,2))/2.0_SP
      C(3,1)=(A(3,1)-A(1,3))/2.0_SP
      C(3,2)=(A(3,2)-A(2,3))/2.0_SP
      C(3,3)=0.0_SP
    CASE DEFAULT
      DO i=1,SIZE(A,1)
        DO j=1,SIZE(A,2)
          !Symmetric part
          B(i,j)=(A(i,j)+A(j,i))/2.0_SP
          !Skew-symmetric part
          C(i,j)=(A(i,j)-A(j,i))/2.0_SP
        ENDDO !j
      ENDDO !i
    END SELECT
 
    EXITS("DecomposeSymmetricSkewMatrix2SP")
    RETURN
999 ERRORSEXITS("DecomposeSymmetricSkewMatrix2SP",err,error)
    RETURN 1
    
  END SUBROUTINE DecomposeSymmetricSkewMatrix2SP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the symmetric and skew-symmetric parts of a double precision matrix2 matrix A i.e., symA in B and skwA in C.
  SUBROUTINE DecomposeSymmetricSkewMatrix2DP(A,B,C,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The matrix A to decompose
    REAL(DP), INTENT(OUT) :: B(:,:) !<On exit, the symmetric part of A i.e., B=symA
    REAL(DP), INTENT(OUT) :: C(:,:) !<On exit, the skew-symmetric part of A i.e., C=skwA
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j
    
    ENTERS("DecomposeSymmetricSkewMatrix2DP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(B,1).OR.SIZE(A,1)/=SIZE(C,1).OR.SIZE(A,2)/=SIZE(B,2).OR.SIZE(A,2)/=SIZE(C,2)) &
      & CALL FlagError("The sizes of the A, B & C matrices are not the same.",err,error,*999)
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("The A matrix is not square.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1))
    CASE(1)
      !Symmetric part
      B(1,1)=A(1,1)
      !Skew-symmetric part
      C(1,1)=0.0_DP
    CASE(2)
      !Symmetric part
      B(1,1)=A(1,1)
      B(1,2)=(A(1,2)+A(2,1))/2.0_DP
      B(2,1)=B(1,2)
      B(2,2)=A(2,2)
      !Skew-symmetric part
      C(1,1)=0.0_DP
      C(1,2)=(A(1,2)-A(2,1))/2.0_DP
      C(2,1)=(A(2,1)-A(1,2))/2.0_DP
      C(2,2)=0.0_DP
    CASE(3)
      !Symmetric part
      B(1,1)=A(1,1)
      B(1,2)=(A(1,2)+A(2,1))/2.0_DP
      B(1,3)=(A(1,3)+A(3,1))/2.0_DP
      B(2,1)=B(1,2)
      B(2,2)=A(2,2)
      B(2,3)=(A(2,3)+A(3,2))/2.0_DP
      B(3,1)=B(1,3)
      B(3,2)=B(2,3)
      B(3,3)=A(3,3)
      !Skew-symmetric part
      C(1,1)=0.0_DP
      C(1,2)=(A(1,2)-A(2,1))/2.0_DP
      C(1,3)=(A(1,3)-A(3,1))/2.0_DP
      C(2,1)=(A(2,1)-A(2,1))/2.0_DP
      C(2,2)=0.0_DP
      C(2,3)=(A(2,3)-A(3,2))/2.0_DP
      C(3,1)=(A(3,1)-A(1,3))/2.0_DP
      C(3,2)=(A(3,2)-A(2,3))/2.0_DP
      C(3,3)=0.0_DP
    CASE DEFAULT
      DO i=1,SIZE(A,1)
        DO j=1,SIZE(A,2)
          !Symmetric part
          B(i,j)=(A(i,j)+A(j,i))/2.0_DP
          !Skew-symmetric part
          C(i,j)=(A(i,j)-A(j,i))/2.0_DP
        ENDDO !j
      ENDDO !i
    END SELECT
 
    EXITS("DecomposeSymmetricSkewMatrix2DP")
    RETURN
999 ERRORSEXITS("DecomposeSymmetricSkewMatrix2DP",err,error)
    RETURN 1
    
  END SUBROUTINE DecomposeSymmetricSkewMatrix2DP

  !
  !================================================================================================================================
  !

  !>Returns the determinant of a full integer matrix A.
  SUBROUTINE DeterminantFullIntg(A,detA,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: A(:,:) !<The matrix to find the determinant of
    INTEGER(INTG), INTENT(OUT) :: detA !<On exit, the determinant of A
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    
    ENTERS("DeterminantFullIntg",err,error,*999)

    detA=0_INTG

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix is not square.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1))
    CASE(1)
      detA=A(1,1)
    CASE(2)
      detA=A(1,1)*A(2,2)-A(2,1)*A(1,2)
    CASE(3)
      detA=A(1,1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+ &
        &  A(1,3)*A(2,1)*A(3,2)-A(1,1)*A(3,2)*A(2,3)- &
        &  A(2,1)*A(1,2)*A(3,3)-A(3,1)*A(2,2)*A(1,3)
    CASE DEFAULT
      CALL FlagError("Matrix size not implemented.",err,error,*999)
    END SELECT

    EXITS("DeterminantFullIntg")
    RETURN
999 ERRORSEXITS("DeterminantFullIntg",err,error)
    RETURN 1
    
  END SUBROUTINE DeterminantFullIntg
  
  !
  !================================================================================================================================
  !

  !>Returns the determinant of a full single precision matrix A.
  SUBROUTINE DeterminantFullSP(A,detA,err,error,*)
    
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The matrix to find the determinant of
    REAL(SP), INTENT(OUT) :: detA !<On exit, the determinant of A
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("DeterminantFullSP",err,error,*999)

    detA=0.0_SP

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix is not square.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1))
    CASE(1)
      detA=A(1,1)
    CASE(2)
      detA=A(1,1)*A(2,2)-A(2,1)*A(1,2)
    CASE(3)
      detA=A(1,1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+ &
        &  A(1,3)*A(2,1)*A(3,2)-A(1,1)*A(3,2)*A(2,3)- &
        &  A(2,1)*A(1,2)*A(3,3)-A(3,1)*A(2,2)*A(1,3)
    CASE DEFAULT
      CALL FlagError("Matrix size not implemented.",err,error,*999)
    END SELECT

    EXITS("DeterminantFullSP")
    RETURN
999 ERRORSEXITS("DeterminantFullSP",err,error)
    RETURN 1
    
  END SUBROUTINE DeterminantFullSP
  
  !
  !================================================================================================================================
  !

  !>Returns the determinant of a full double precision matrix A
  SUBROUTINE DeterminantFullDP(A,detA,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The matrix to find the determinant of
    REAL(DP), INTENT(OUT) :: detA !<On exit, the determinant of A
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    
    ENTERS("DeterminantFullDP",err,error,*999)

    detA=0.0_DP

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix is not square.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1))
    CASE(1)
      detA=A(1,1)
    CASE(2)
      detA=A(1,1)*A(2,2)-A(2,1)*A(1,2)
    CASE(3)
      detA=A(1,1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+ &
        &  A(1,3)*A(2,1)*A(3,2)-A(1,1)*A(3,2)*A(2,3)- &
        &  A(2,1)*A(1,2)*A(3,3)-A(3,1)*A(2,2)*A(1,3)
    CASE DEFAULT
      CALL FlagError("Matrix size not implemented.",err,error,*999)
    END SELECT

    EXITS("DeterminantFullDP")
    RETURN
999 ERRORSEXITS("DeterminantFullDP",err,error)
    RETURN 1
    
  END SUBROUTINE DeterminantFullDP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the matrix2-matrix2 dot (matrix) product of the single precision matrix A*B in C.
  SUBROUTINE DotProductMatrix2Matrix2SP(A,B,C,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The first matrix A
    REAL(SP), INTENT(IN) :: B(:,:) !<The second matrix B
    REAL(SP), INTENT(OUT) :: C(:,:) !<On exit, the dot product matrix C=A*B
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,k
    
    ENTERS("DotProductMatrix2Matrix2SP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,2)/=SIZE(B,1).OR.SIZE(A,1)/=SIZE(C,1).OR.SIZE(B,2)/=SIZE(C,2)) CALL FlagError("Invalid matrix sizes.",err,error,*999)
#endif    
    
    IF(SIZE(A,1)==SIZE(A,2).AND.SIZE(A,1)==SIZE(B,2)) THEN
      !A and B matrices are square and the same size
      SELECT CASE(SIZE(A,1))
      CASE(1)
        C(1,1)=A(1,1)*B(1,1)
      CASE(2)
        C(1,1)=A(1,1)*B(1,1)+A(1,2)*B(2,1)
        C(1,2)=A(1,1)*B(1,2)+A(1,2)*B(2,2)
        C(2,1)=A(2,1)*B(1,1)+A(2,2)*B(2,1)
        C(2,2)=A(2,1)*B(1,2)+A(2,2)*B(2,2)
      CASE(3)
        C(1,1)=A(1,1)*B(1,1)+A(1,2)*B(2,1)+A(1,3)*B(3,1)
        C(1,2)=A(1,1)*B(1,2)+A(1,2)*B(2,2)+A(1,3)*B(3,2)
        C(1,3)=A(1,1)*B(1,3)+A(1,2)*B(2,3)+A(1,3)*B(3,3)        
        C(2,1)=A(2,1)*B(1,1)+A(2,2)*B(2,1)+A(2,3)*B(3,1)
        C(2,2)=A(2,1)*B(1,2)+A(2,2)*B(2,2)+A(2,3)*B(3,2)
        C(2,3)=A(2,1)*B(1,3)+A(2,2)*B(2,3)+A(2,3)*B(3,3)
        C(3,1)=A(3,1)*B(1,1)+A(3,2)*B(2,1)+A(3,3)*B(3,1)
        C(3,2)=A(3,1)*B(1,2)+A(3,2)*B(2,2)+A(3,3)*B(3,2)
        C(3,3)=A(3,1)*B(1,3)+A(3,2)*B(2,3)+A(3,3)*B(3,3)
      CASE DEFAULT
        DO i=1,SIZE(A,1)
          DO j=1,SIZE(A,1)
            C(i,j) = 0.0_SP
            DO k=1,SIZE(A,1)
              C(i,j) = C(i,j) + A(i,k)*B(k,j)
            ENDDO !k
          ENDDO !j
        ENDDO !i
      END SELECT
    ELSE
      !A and B matrices are not square or not the same size
      DO i=1,SIZE(A,1)
        DO j=1,SIZE(B,2)
          C(i,j) = 0.0_SP
          DO k=1,SIZE(A,2)
            C(i,j) = C(i,j) + A(i,k)*B(k,j)
          ENDDO !k
        ENDDO !j
      ENDDO !i
    ENDIF
 
    EXITS("DotProductMatrix2Matrix2SP")
    RETURN
999 ERRORSEXITS("DotProductMatrix2Matrix2SP",err,error)
    RETURN 1
    
  END SUBROUTINE DotProductMatrix2Matrix2SP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the matrix2-matrix2 dot(matrix) product of the double precision matrix A*B in C.
  SUBROUTINE DotProductMatrix2Matrix2DP(A,B,C,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The A matrix
    REAL(DP), INTENT(IN) :: B(:,:) !<The B matrix
    REAL(DP), INTENT(OUT) :: C(:,:) !<On exit, the product matrix C=A*B
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,k
        
    ENTERS("DotProductMatrix2Matrix2DP",err,error,*999)

#ifdef WITH_PRECHEKCS    
    IF(SIZE(A,2)/=SIZE(B,1).OR.SIZE(A,1)/=SIZE(C,1).OR.SIZE(B,2)/=SIZE(C,2)) CALL FlagError("Invalid matrix sizes.",err,error,*999)
#endif    
    
    IF(SIZE(A,1)==SIZE(A,2).AND.SIZE(A,1)==SIZE(B,2)) THEN
      !A and B matrices are square and the same size
      SELECT CASE(SIZE(A,1))
      CASE(1)
        C(1,1)=A(1,1)*B(1,1)
      CASE(2)
        C(1,1)=A(1,1)*B(1,1)+A(1,2)*B(2,1)
        C(1,2)=A(1,1)*B(1,2)+A(1,2)*B(2,2)
        C(2,1)=A(2,1)*B(1,1)+A(2,2)*B(2,1)
        C(2,2)=A(2,1)*B(1,2)+A(2,2)*B(2,2)
      CASE(3)
        C(1,1)=A(1,1)*B(1,1)+A(1,2)*B(2,1)+A(1,3)*B(3,1)
        C(1,2)=A(1,1)*B(1,2)+A(1,2)*B(2,2)+A(1,3)*B(3,2)
        C(1,3)=A(1,1)*B(1,3)+A(1,2)*B(2,3)+A(1,3)*B(3,3)        
        C(2,1)=A(2,1)*B(1,1)+A(2,2)*B(2,1)+A(2,3)*B(3,1)
        C(2,2)=A(2,1)*B(1,2)+A(2,2)*B(2,2)+A(2,3)*B(3,2)
        C(2,3)=A(2,1)*B(1,3)+A(2,2)*B(2,3)+A(2,3)*B(3,3)
        C(3,1)=A(3,1)*B(1,1)+A(3,2)*B(2,1)+A(3,3)*B(3,1)
        C(3,2)=A(3,1)*B(1,2)+A(3,2)*B(2,2)+A(3,3)*B(3,2)
        C(3,3)=A(3,1)*B(1,3)+A(3,2)*B(2,3)+A(3,3)*B(3,3)
      CASE DEFAULT
        DO i=1,SIZE(A,1)
          DO j=1,SIZE(A,1)
            C(i,j) = 0.0_DP
            DO k=1,SIZE(A,1)
              C(i,j) = C(i,j) + A(i,k)*B(k,j)
            ENDDO !k
          ENDDO !j
        ENDDO !i
      END SELECT
    ELSE
      !A and B matrices are not square or not the same size
      DO i=1,SIZE(A,1)
        DO j=1,SIZE(B,2)
          C(i,j) = 0.0_DP
          DO k=1,SIZE(A,2)
            C(i,j) = C(i,j) + A(i,k)*B(k,j)
          ENDDO !k
        ENDDO !j
      ENDDO !i      
    ENDIF

    EXITS("DotProductMatrix2Matrix2DP")
    RETURN
999 ERRORSEXITS("DotProductMatrix2Matrix2DP",err,error)
    RETURN 1
    
  END SUBROUTINE DotProductMatrix2Matrix2DP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the matrix-vector dot product of the single precision matrice and vector, A.b in c.
  SUBROUTINE DotProductMatrixVectorSP(A,b,c,err,error,*)

    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The A matrix
    REAL(SP), INTENT(IN) :: b(:) !<The b vector
    REAL(SP), INTENT(OUT) :: c(:) !<On exit, the dot product vector c=A.b
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j

    ENTERS("DotProductMatrixVectorSP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,2)/=SIZE(b,1)) CALL FlagError("The number of columns in A does not match the number of rows in b.",err,error,*999)
    IF(SIZE(A,1)/=SIZE(c,1)) CALL FlagError("The number of rows of A does not match the number of rows in c.",err,error,*999)
#endif    
    
    IF(SIZE(A,1)==SIZE(A,2)) THEN
      !A matrix is square
      SELECT CASE(SIZE(A,1))
      CASE(1)
        c(1)=A(1,1)*b(1)
      CASE(2)
        c(1)=A(1,1)*b(1)+A(1,2)*b(2)
        c(2)=A(2,1)*b(1)+A(2,2)*b(2)
      CASE(3)
        c(1)=A(1,1)*b(1)+A(1,2)*b(2)+A(1,3)*b(3)
        c(2)=A(2,1)*b(1)+A(2,2)*b(2)+A(2,3)*b(3)
        c(3)=A(3,1)*b(1)+A(3,2)*b(2)+A(3,3)*b(3)
      CASE DEFAULT
        DO i=1,SIZE(A,1)
          c(i) = 0.0_SP
          DO j=1,SIZE(A,1)
            c(i) = c(i) + A(i,j)*b(j)
          ENDDO !j
        ENDDO !i
      END SELECT
    ELSE
      !A matrix is not square
      DO i=1,SIZE(A,1)
        c(i) = 0.0_SP
        DO j=1,SIZE(A,2)
          c(i) = c(i) + A(i,j)*b(j)
        ENDDO !j
      ENDDO !i
    ENDIF
    
    EXITS("DotProductMatrixVectorSP")
    RETURN
999 ERRORSEXITS("DotProductMatrixVectorSP",err,error)
    RETURN 1
    
  END SUBROUTINE DotProductMatrixVectorSP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the matrix-vector dot product of the double precision matrice and vector, A.b in c.
  SUBROUTINE DotProductMatrixVectorDP(A,b,c,err,error,*)

    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The A matrix
    REAL(DP), INTENT(IN) :: b(:)   !<The b vector
    REAL(DP), INTENT(OUT) :: c(:)  !<On exit, the inner product vector c=A*b
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j

    ENTERS("DotProductMatrixVectorDP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,2)/=SIZE(b,1)) CALL FlagError("The number of columns in A does not match the number of rows in b.",err,error,*999)
    IF(SIZE(A,1)/=SIZE(c,1)) CALL FlagError("The number of rows of A does not match the number of rows in c.",err,error,*999)
#endif    
    
    IF(SIZE(A,1)==SIZE(A,2)) THEN
      !A matrix is square
      SELECT CASE(SIZE(A,1))
      CASE(1)
        c(1)=A(1,1)*b(1)
      CASE(2)
        c(1)=A(1,1)*b(1)+A(1,2)*b(2)
        c(2)=A(2,1)*b(1)+A(2,2)*b(2)
      CASE(3)
        c(1)=A(1,1)*b(1)+A(1,2)*b(2)+A(1,3)*b(3)
        c(2)=A(2,1)*b(1)+A(2,2)*b(2)+A(2,3)*b(3)
        c(3)=A(3,1)*b(1)+A(3,2)*b(2)+A(3,3)*b(3)
      CASE DEFAULT
        DO i=1,SIZE(A,1)
          c(i) = 0.0_DP
          DO j=1,SIZE(A,1)
            c(i) = c(i) + A(i,j)*b(j)
          ENDDO !j
        ENDDO !i
      END SELECT
    ELSE
      !A matrix is not square
      DO i=1,SIZE(A,1)
        c(i) = 0.0_DP
        DO j=1,SIZE(A,2)
          c(i) = c(i) + A(i,j)*b(j)
        ENDDO !j
      ENDDO !i
    ENDIF
    
    EXITS("DotProductMatrixVectorDP")
    RETURN
999 ERRORSEXITS("DotProductMatrixVectorDP",err,error)
    RETURN 1
    
  END SUBROUTINE DotProductMatrixVectorDP

  !
  !
  !================================================================================================================================
  !

  !>Calculates and returns the vector-vector dot product of the single precision vectors a.b in c.
  SUBROUTINE DotProductVectorVectorSP(a,b,c,err,error,*)

    !Argument variables
    REAL(SP), INTENT(IN) :: a(:) !<The a vector
    REAL(SP), INTENT(IN) :: b(:) !<The b vector
    REAL(SP), INTENT(OUT) :: c !<On exit, the vector dot product c=a.b
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("DotProductVectorVectorSP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(a,1)/=SIZE(b,1)) CALL FlagError("The a and b vectors are of different lengths.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(a,1))
    CASE(1)
      c = a(1)*b(1)
    CASE(2)
      c = a(1)*b(1)+a(2)*b(2)
    CASE(3)
      c = a(1)*b(1)+a(2)*b(2)+a(3)*b(3) 
    CASE DEFAULT
      c = DOT_PRODUCT(a,b)
    END SELECT
    
    EXITS("DotProductVectorVectorSP")
    RETURN
999 ERRORSEXITS("DotProductVectorVectorSP",err,error)
    RETURN 1
    
  END SUBROUTINE DotProductVectorVectorSP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the vector-vector dot product of the double precision vectors a.b in c.
  SUBROUTINE DotProductVectorVectorDP(a,b,c,err,error,*)

    !Argument variables
    REAL(DP), INTENT(IN) :: a(:) !<The a vector
    REAL(DP), INTENT(IN) :: b(:) !<The b vector
    REAL(DP), INTENT(OUT) :: c !<On exit, the vector dot product c=a.b
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("DotProductVectorVectorDP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(a,1)/=SIZE(b,1)) CALL FlagError("The a and b vectors are of different lengths.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(a,1))
    CASE(1)
      c = a(1)*b(1)
    CASE(2)
      c = a(1)*b(1)+a(2)*b(2)
    CASE(3)
      c = a(1)*b(1)+a(2)*b(2)+a(3)*b(3) 
    CASE DEFAULT
      c = DOT_PRODUCT(a,b)
    END SELECT
    
    EXITS("DotProductVectorVectorDP")
    RETURN
999 ERRORSEXITS("DotProductVectorVectorDP",err,error)
    RETURN 1
    
  END SUBROUTINE DotProductVectorVectorDP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the matrix-transpose dot(matrix) product of the single precision matrix A^T*B in C.
  SUBROUTINE DotTransposeProductMatrix2Matrix2SP(A,B,C,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The first matrix A
    REAL(SP), INTENT(IN) :: B(:,:) !<The second matrix B
    REAL(SP), INTENT(OUT) :: C(:,:) !<On exit, the dot product matrix C=A^T*B
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,k
    
    ENTERS("DotTransposeProductMatrix2Matrix2SP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,2)/=SIZE(B,1).OR.SIZE(A,1)/=SIZE(C,1).OR.SIZE(B,2)/=SIZE(C,2)) CALL FlagError("Invalid matrix sizes.",err,error,*999)
#endif    

    IF(SIZE(A,1)==SIZE(A,2).AND.SIZE(A,1)==SIZE(B,2)) THEN
      !A and B matrices are square and the same size
      SELECT CASE(SIZE(A,1))
      CASE(1)
        C(1,1)=A(1,1)*B(1,1)
      CASE(2)
        C(1,1)=A(1,1)*B(1,1)+A(2,1)*B(2,1)
        C(1,2)=A(1,1)*B(1,2)+A(2,1)*B(2,2)
        C(2,1)=A(1,2)*B(1,1)+A(2,2)*B(2,1)
        C(2,2)=A(1,2)*B(1,2)+A(2,2)*B(2,2)
      CASE(3)
        C(1,1)=A(1,1)*B(1,1)+A(2,1)*B(2,1)+A(3,1)*B(3,1)
        C(1,2)=A(1,1)*B(1,2)+A(2,1)*B(2,2)+A(3,1)*B(3,2)
        C(1,3)=A(1,1)*B(1,3)+A(2,1)*B(2,3)+A(3,1)*B(3,3)        
        C(2,1)=A(1,2)*B(1,1)+A(2,2)*B(2,1)+A(3,2)*B(3,1)
        C(2,2)=A(1,2)*B(1,2)+A(2,2)*B(2,2)+A(3,2)*B(3,2)
        C(2,3)=A(1,2)*B(1,3)+A(2,2)*B(2,3)+A(3,2)*B(3,3)
        C(3,1)=A(1,3)*B(1,1)+A(2,3)*B(2,1)+A(3,3)*B(3,1)
        C(3,2)=A(1,3)*B(1,2)+A(2,3)*B(2,2)+A(3,3)*B(3,2)
        C(3,3)=A(1,3)*B(1,3)+A(2,3)*B(2,3)+A(3,3)*B(3,3)
      CASE DEFAULT
        DO i=1,SIZE(A,1)
          DO j=1,SIZE(A,1)
            C(i,j) = 0.0_SP
            DO k=1,SIZE(A,1)
              C(i,j) = C(i,j) + A(k,i)*B(k,j)
            ENDDO !k
          ENDDO !j
        ENDDO !i
      END SELECT
    ELSE
      !A and B matrices are not square or not the same size
      DO i=1,SIZE(A,2)
        DO j=1,SIZE(B,2)
          C(i,j) = 0.0_SP
          DO k=1,SIZE(A,1)
            C(i,j) = C(i,j) + A(k,i)*B(k,j)
          ENDDO !k
        ENDDO !j
      ENDDO !i      
    ENDIF

    EXITS("DotTransposeProductMatrix2Matrix2SP")
    RETURN
999 ERRORSEXITS("DotTransposeProductMatrix2Matrix2SP",err,error)
    RETURN 1
    
  END SUBROUTINE DotTransposeProductMatrix2Matrix2SP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the matrix-transpose dot(matrix) product of the double precision matrix A^T*B in C.
  SUBROUTINE DotTransposeProductMatrix2Matrix2DP(A,B,C,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The first matrix A
    REAL(DP), INTENT(IN) :: B(:,:) !<The second matrix B
    REAL(DP), INTENT(OUT) :: C(:,:) !<On exit, the product matrix C=A^T*B
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,k
    
    ENTERS("DotTransposeProductMatrix2Matrix2DP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,2)/=SIZE(B,1).OR.SIZE(A,1)/=SIZE(C,1).OR.SIZE(B,2)/=SIZE(C,2)) CALL FlagError("Invalid matrix sizes.",err,error,*999)
#endif    

    IF(SIZE(A,1)==SIZE(A,2).AND.SIZE(A,1)==SIZE(B,2)) THEN
      !A and B matrices are square and the same size
      SELECT CASE(SIZE(A,1))
      CASE(1)
        C(1,1)=A(1,1)*B(1,1)
      CASE(2)
        C(1,1)=A(1,1)*B(1,1)+A(2,1)*B(2,1)
        C(1,2)=A(1,1)*B(1,2)+A(2,1)*B(2,2)
        C(2,1)=A(1,2)*B(1,1)+A(2,2)*B(2,1)
        C(2,2)=A(1,2)*B(1,2)+A(2,2)*B(2,2)
      CASE(3)
        C(1,1)=A(1,1)*B(1,1)+A(2,1)*B(2,1)+A(3,1)*B(3,1)
        C(1,2)=A(1,1)*B(1,2)+A(2,1)*B(2,2)+A(3,1)*B(3,2)
        C(1,3)=A(1,1)*B(1,3)+A(2,1)*B(2,3)+A(3,1)*B(3,3)        
        C(2,1)=A(1,2)*B(1,1)+A(2,2)*B(2,1)+A(3,2)*B(3,1)
        C(2,2)=A(1,2)*B(1,2)+A(2,2)*B(2,2)+A(3,2)*B(3,2)
        C(2,3)=A(1,2)*B(1,3)+A(2,2)*B(2,3)+A(3,2)*B(3,3)
        C(3,1)=A(1,3)*B(1,1)+A(2,3)*B(2,1)+A(3,3)*B(3,1)
        C(3,2)=A(1,3)*B(1,2)+A(2,3)*B(2,2)+A(3,3)*B(3,2)
        C(3,3)=A(1,3)*B(1,3)+A(2,3)*B(2,3)+A(3,3)*B(3,3)
      CASE DEFAULT
        DO i=1,SIZE(A,1)
          DO j=1,SIZE(A,1)
            C(i,j) = 0.0_SP
            DO k=1,SIZE(A,1)
              C(i,j) = C(i,j) + A(k,i)*B(k,j)
            ENDDO !k
          ENDDO !j
        ENDDO !i
      END SELECT
    ELSE
      !A and B matrices are not square or not the same size
      DO i=1,SIZE(A,2)
        DO j=1,SIZE(B,2)
          C(i,j) = 0.0_SP
          DO k=1,SIZE(A,1)
            C(i,j) = C(i,j) + A(k,i)*B(k,j)
          ENDDO !k
        ENDDO !j
      ENDDO !i      
    ENDIF

    EXITS("DotTransposeProductMatrix2Matrix2DP")
    RETURN
999 ERRORSEXITS("DotTransposeProductMatrix2Matrix2DP",err,error)
    RETURN 1
    
  END SUBROUTINE DotTransposeProductMatrix2Matrix2DP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the matrix-transpose vector dot product of the single precision matrix and vector A^T*b in c.
  SUBROUTINE DotTransposeProductMatrixVectorSP(A,b,c,err,error,*)

    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The A matrix
    REAL(SP), INTENT(IN) :: b(:) !<The b vector
    REAL(SP), INTENT(OUT) :: c(:) !<On exit, the dot product vector c=A^T.b
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j

    ENTERS("DotTransposeProductMatrixVectorSP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(b,1)) CALL FlagError("The number of rows in A does not match the number of rows in b.",err,error,*999)
    IF(SIZE(A,2)/=SIZE(c,1)) CALL FlagError("The number of columns in A does not match the number of rows in c.",err,error,*999)
#endif    
  
    IF(SIZE(A,1)==SIZE(A,2)) THEN
      !A matrix is square
      SELECT CASE(SIZE(A,2))
      CASE(1)
        c(1)=A(1,1)*b(1)
      CASE(2)
        c(1)=A(1,1)*b(1)+A(2,1)*b(2)
        c(2)=A(1,2)*b(1)+A(2,2)*b(2)
      CASE(3)
        c(1)=A(1,1)*b(1)+A(2,1)*b(2)+A(3,1)*b(3)
        c(2)=A(1,2)*b(1)+A(2,2)*b(2)+A(3,2)*b(3)
        c(3)=A(1,3)*b(1)+A(2,3)*b(2)+A(3,3)*b(3)
      CASE DEFAULT
        DO i=1,SIZE(A,1)
          c(i) = 0.0_SP
          DO j=1,SIZE(A,1)
            c(i) = c(i) + A(j,i)*b(j)
          ENDDO !j
        ENDDO !i
      END SELECT
    ELSE
      !A matrix is not square
      DO i=1,SIZE(A,2)
        c(i) = 0.0_SP
        DO j=1,SIZE(A,1)
          c(i) = c(i) + A(j,i)*b(j)
        ENDDO !j
      ENDDO !i
    ENDIF

    EXITS("DotTransposeProductMatrixVectorSP")
    RETURN
999 ERRORSEXITS("DotTransposeProductMatrixVectorSP",err,error)
    RETURN 1
    
  END SUBROUTINE DotTransposeProductMatrixVectorSP
   
  !
  !================================================================================================================================
  !

  !>Calculates and returns the matrix-transpose vector dot product of the double precision matrix and vector A^T.b in c.
  SUBROUTINE DotTransposeProductMatrixVectorDP(A,b,c,err,error,*)

    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The A matrix
    REAL(DP), INTENT(IN) :: b(:) !<The b vector
    REAL(DP), INTENT(OUT) :: c(:) !<On exit, the inner product vector c=A^T.b
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j

    ENTERS("DotTransposeProductMatrixVectorDP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(b,1)) CALL FlagError("The number of rows in A does not match the number of rows in b.",err,error,*999)
    IF(SIZE(A,2)/=SIZE(c,1)) CALL FlagError("The number of columns in A does not match the number of rows in c.",err,error,*999)
#endif    
  
    IF(SIZE(A,1)==SIZE(A,2)) THEN
      !A matrix is square
      SELECT CASE(SIZE(A,2))
      CASE(1)
        c(1)=A(1,1)*b(1)
      CASE(2)
        c(1)=A(1,1)*b(1)+A(2,1)*b(2)
        c(2)=A(1,2)*b(1)+A(2,2)*b(2)
      CASE(3)
        c(1)=A(1,1)*b(1)+A(2,1)*b(2)+A(3,1)*b(3)
        c(2)=A(1,2)*b(1)+A(2,2)*b(2)+A(3,2)*b(3)
        c(3)=A(1,3)*b(1)+A(2,3)*b(2)+A(3,3)*b(3)
      CASE DEFAULT
        DO i=1,SIZE(A,1)
          c(i) = 0.0_DP
          DO j=1,SIZE(A,1)
            c(i) = c(i) + A(j,i)*b(j)
          ENDDO !j
        ENDDO !i
      END SELECT
    ELSE
      !A matrix is not square
      DO i=1,SIZE(A,2)
        c(i) = 0.0_DP
        DO j=1,SIZE(A,1)
          c(i) = c(i) + A(j,i)*b(j)
        ENDDO !j
      ENDDO !i
    ENDIF

    EXITS("DotTransposeProductMatrixVectorDP")
    RETURN
999 ERRORSEXITS("DotTransposeProductMatrixVectorDP",err,error)
    RETURN 1
    
  END SUBROUTINE DotTransposeProductMatrixVectorDP
   
  !
  !================================================================================================================================
  !
  
  !>Calculates and returns the matrix-product-transpose of the single precision matrix A*B^T in C.
  SUBROUTINE DotProductTransposeMatrix2Matrix2SP(A,B,C,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The first matrix A
    REAL(SP), INTENT(IN) :: B(:,:) !<The second matrix B
    REAL(SP), INTENT(OUT) :: C(:,:) !<On exit, the product matrix C=A*B^T
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,k
    
    ENTERS("DotProductTransposeMatrix2Matrix2SP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,2)/=SIZE(B,2).OR.SIZE(A,1)/=SIZE(C,1).OR.SIZE(B,1)/=SIZE(C,2)) CALL FlagError("Invalid matrix sizes.",err,error,*999)
#endif    

    IF(SIZE(A,1)==SIZE(A,2).AND.SIZE(A,1)==SIZE(B,1)) THEN
      !A and B matrices are square and the same size
      SELECT CASE(SIZE(A,1))
      CASE(1)
        C(1,1)=A(1,1)*B(1,1)
      CASE(2)
        C(1,1)=A(1,1)*B(1,1)+A(1,2)*B(1,2)
        C(1,2)=A(1,1)*B(2,1)+A(1,2)*B(2,2)
        C(2,1)=A(2,1)*B(1,1)+A(2,2)*B(1,2)
        C(2,2)=A(2,1)*B(2,1)+A(2,2)*B(2,2)
      CASE(3)
        C(1,1)=A(1,1)*B(1,1)+A(1,2)*B(1,2)+A(1,3)*B(1,3)
        C(1,2)=A(1,1)*B(2,1)+A(1,2)*B(2,2)+A(1,3)*B(2,3)
        C(1,3)=A(1,1)*B(3,1)+A(1,2)*B(3,2)+A(1,3)*B(3,3)        
        C(2,1)=A(2,1)*B(1,1)+A(2,2)*B(1,2)+A(2,3)*B(1,3)
        C(2,2)=A(2,1)*B(2,1)+A(2,2)*B(2,2)+A(2,3)*B(2,3)
        C(2,3)=A(2,1)*B(3,1)+A(2,2)*B(3,2)+A(2,3)*B(3,3)
        C(3,1)=A(3,1)*B(1,1)+A(3,2)*B(1,2)+A(3,3)*B(1,3)
        C(3,2)=A(3,1)*B(2,1)+A(3,2)*B(2,2)+A(3,3)*B(2,3)
        C(3,3)=A(3,1)*B(3,1)+A(3,2)*B(3,2)+A(3,3)*B(3,3)
      CASE DEFAULT
        DO i=1,SIZE(A,1)
          DO j=1,SIZE(A,1)
            C(i,j) = 0.0_SP
            DO k=1,SIZE(A,1)
              C(i,j) = C(i,j) + A(i,k)*B(j,k)
            ENDDO !k
          ENDDO !j
        ENDDO !i            
      END SELECT
    ELSE
      !A and B matrices are not square or not the same size
      DO i=1,SIZE(A,1)
        DO j=1,SIZE(B,1)
          C(i,j) = 0.0_SP
          DO k=1,SIZE(A,2)
            C(i,j) = C(i,j) + A(i,k)*B(j,k)
          ENDDO !k
        ENDDO !j
      ENDDO !i            
    ENDIF

    EXITS("DotProductTransposeMatrix2Matrix2SP")
    RETURN
999 ERRORSEXITS("DotProductTransposeMatrix2Matrix2SP",err,error)
    RETURN 1
    
  END SUBROUTINE DotProductTransposeMatrix2Matrix2SP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the matrix-product-transpose of the double precision matrix A*B^T in C.
  SUBROUTINE DotProductTransposeMatrix2Matrix2DP(A,B,C,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The A matrix
    REAL(DP), INTENT(IN) :: B(:,:) !<The B matrix
    REAL(DP), INTENT(OUT) :: C(:,:) !<On exit, the product matrix C=A*B^T
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,k
       
    ENTERS("DotProductTransposeMatrix2Matrix2DP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,2)/=SIZE(B,2).OR.SIZE(A,1)/=SIZE(C,1).OR.SIZE(B,1)/=SIZE(C,2)) CALL FlagError("Invalid matrix sizes.",err,error,*999)
#endif    

    IF(SIZE(A,1)==SIZE(A,2).AND.SIZE(A,1)==SIZE(B,1)) THEN
      !A and B matrices are square and the same size
      SELECT CASE(SIZE(A,1))
      CASE(1)
        C(1,1)=A(1,1)*B(1,1)
      CASE(2)
        C(1,1)=A(1,1)*B(1,1)+A(1,2)*B(1,2)
        C(1,2)=A(1,1)*B(2,1)+A(1,2)*B(2,2)
        C(2,1)=A(2,1)*B(1,1)+A(2,2)*B(1,2)
        C(2,2)=A(2,1)*B(2,1)+A(2,2)*B(2,2)
      CASE(3)
        C(1,1)=A(1,1)*B(1,1)+A(1,2)*B(1,2)+A(1,3)*B(1,3)
        C(1,2)=A(1,1)*B(2,1)+A(1,2)*B(2,2)+A(1,3)*B(2,3)
        C(1,3)=A(1,1)*B(3,1)+A(1,2)*B(3,2)+A(1,3)*B(3,3)        
        C(2,1)=A(2,1)*B(1,1)+A(2,2)*B(1,2)+A(2,3)*B(1,3)
        C(2,2)=A(2,1)*B(2,1)+A(2,2)*B(2,2)+A(2,3)*B(2,3)
        C(2,3)=A(2,1)*B(3,1)+A(2,2)*B(3,2)+A(2,3)*B(3,3)
        C(3,1)=A(3,1)*B(1,1)+A(3,2)*B(1,2)+A(3,3)*B(1,3)
        C(3,2)=A(3,1)*B(2,1)+A(3,2)*B(2,2)+A(3,3)*B(2,3)
        C(3,3)=A(3,1)*B(3,1)+A(3,2)*B(3,2)+A(3,3)*B(3,3)
      CASE DEFAULT
        DO i=1,SIZE(A,1)
          DO j=1,SIZE(A,1)
            C(i,j) = 0.0_DP
            DO k=1,SIZE(A,1)
              C(i,j) = C(i,j) + A(i,k)*B(j,k)
            ENDDO !k
          ENDDO !j
        ENDDO !i            
      END SELECT
    ELSE
      !A and B matrices are not square or not the same size
      DO i=1,SIZE(A,1)
        DO j=1,SIZE(B,1)
          C(i,j) = 0.0_DP
          DO k=1,SIZE(A,2)
            C(i,j) = C(i,j) + A(i,k)*B(j,k)
          ENDDO !k
        ENDDO !j
      ENDDO !i            
    ENDIF

    EXITS("DotProductTransposeMatrix2Matrix2DP")
    RETURN
999 ERRORSEXITS("DotProductTransposeMatrix2Matrix2DP",err,error)
    RETURN 1
    
  END SUBROUTINE DotProductTransposeMatrix2Matrix2DP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the matrix2-matrix2 double dot product of the single precision matrices A:B in c.
  SUBROUTINE DoubleDotProductMatrix2Matrix2SP(A,B,c,err,error,*)

    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The A matrix
    REAL(SP), INTENT(IN) :: B(:,:) !<The B matrix
    REAL(SP), INTENT(OUT) :: c !<On exit, the matrix double dot product c=A:B
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j
    
    ENTERS("DoubleDotProductMatrix2Matrix2SP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("The A matrix is not square.",err,error,*999)
    IF(SIZE(A,1)/=SIZE(B,1).OR.SIZE(A,2)/=SIZE(B,2)) &
      & CALL FlagError("The B matrix is not the same size as the A matrix.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1))
    CASE(1)
      c = A(1,1)*B(1,1)
    CASE(2)
      c = A(1,1)*B(1,1)+A(2,1)*B(2,1)+ &
        & A(1,2)*B(1,2)+A(2,2)*B(2,2)
    CASE(3)
      c = A(1,1)*B(1,1)+A(2,1)*B(2,1)+A(3,1)*B(3,1) + &
        & A(1,2)*B(1,2)+A(2,2)*B(2,2)+A(3,2)*B(3,2) + &
        & A(1,3)*B(1,3)+A(2,3)*B(2,3)+A(3,3)*B(3,3)
    CASE DEFAULT
      c = 0.0_SP
      DO i=1,SIZE(A,1)
        DO j=1,SIZE(A,2)
          c = c + A(i,j)*B(i,j)
        ENDDO !j
      ENDDO !i
    END SELECT
    
    EXITS("DoubleDotProductMatrix2Matrix2SP")
    RETURN
999 ERRORSEXITS("DoubleDotProductMatrix2Matrix2SP",err,error)
    RETURN 1
    
  END SUBROUTINE DoubleDotProductMatrix2Matrix2SP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the matrix2-matrix2 double dot product of the double precision matrices A:B in c.
  SUBROUTINE DoubleDotProductMatrix2Matrix2DP(A,B,c,err,error,*)

    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The A matrix
    REAL(DP), INTENT(IN) :: B(:,:) !<The B matrix
    REAL(DP), INTENT(OUT) :: c !<On exit, the inner product c=A:B
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j

    ENTERS("DoubleDotProductMatrix2Matrix2DP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("The A matrix is not square.",err,error,*999)
    IF(SIZE(A,1)/=SIZE(B,1).OR.SIZE(A,2)/=SIZE(B,2)) &
      & CALL FlagError("The B matrix is not the same size as the A matrix.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1))
    CASE(1)
      c = A(1,1)*B(1,1)
    CASE(2)
      c = A(1,1)*B(1,1)+A(2,1)*B(2,1)+ &
        & A(1,2)*B(1,2)+A(2,2)*B(2,2)
    CASE(3)
      c = A(1,1)*B(1,1)+A(2,1)*B(2,1)+A(3,1)*B(3,1) + &
        & A(1,2)*B(1,2)+A(2,2)*B(2,2)+A(3,2)*B(3,2) + &
        & A(1,3)*B(1,3)+A(2,3)*B(2,3)+A(3,3)*B(3,3)
    CASE DEFAULT
      c = 0.0_DP
      DO i=1,SIZE(A,1)
        DO j=1,SIZE(A,2)
          c = c + A(i,j)*B(i,j)
        ENDDO !j
      ENDDO !i
    END SELECT
    
    EXITS("DoubleDotProductMatrix2Matrix2DP")
    RETURN
999 ERRORSEXITS("DoubleDotProductMatrix2Matrix2DP",err,error)
    RETURN 1
    
  END SUBROUTINE DoubleDotProductMatrix2Matrix2DP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the matrix2-matrix4 double dot product of the single precision matrices A:B in C.
  SUBROUTINE DoubleDotProductMatrix2Matrix4SP(A,B,C,err,error,*)

    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The A matrix
    REAL(SP), INTENT(IN) :: B(:,:,:,:) !<The B matrix
    REAL(SP), INTENT(OUT) :: C(:,:) !<On exit, the matrix double dot product C=A:B
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,k,l
    
    ENTERS("DoubleDotProductMatrix2Matrix4SP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(B,1).OR.SIZE(A,2)/=SIZE(B,2)) &
      & CALL FlagError("The number of dimensions of the A matrix to no conform to the first two dimensions of the B matrix.", &
      & err,error,*999)
    IF(SIZE(C,1)/=SIZE(B,3).OR.SIZE(C,2)/=SIZE(B,4))  &
      & CALL FlagError("The C matrix is not the same size as the last two dimensions of the B matrix.",err,error,*999)
#endif    
    
    DO i=1,SIZE(B,3)
      DO j=1,SIZE(B,4)
        C(i,j) = 0.0_SP
        DO k=1,SIZE(B,1)
          DO l=1,SIZE(B,2)
            C(i,j) = C(i,j) + A(k,l)*B(k,l,i,j)
          ENDDO !l
        ENDDO !k
      ENDDO !j
    ENDDO !i
    
    EXITS("DoubleDotProductMatrix2Matrix4SP")
    RETURN
999 ERRORSEXITS("DoubleDotProductMatrix2Matrix4SP",err,error)
    RETURN 1
    
  END SUBROUTINE DoubleDotProductMatrix2Matrix4SP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the matrix2-matrix4 double dot product of the double precision matrices A:B in C.
  SUBROUTINE DoubleDotProductMatrix2Matrix4DP(A,B,C,err,error,*)

    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The A matrix
    REAL(DP), INTENT(IN) :: B(:,:,:,:) !<The B matrix
    REAL(DP), INTENT(OUT) :: C(:,:) !<On exit, the matrix double dot product C=A:B
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,k,l
    
    ENTERS("DoubleDotProductMatrix2Matrix4DP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(B,1).OR.SIZE(A,2)/=SIZE(B,2)) &
      & CALL FlagError("The number of dimensions of the A matrix to no conform to the first two dimensions of the B matrix.", &
      & err,error,*999)
    IF(SIZE(C,1)/=SIZE(B,3).OR.SIZE(C,2)/=SIZE(B,4))  &
      & CALL FlagError("The C matrix is not the same size as the last two dimensions of the B matrix.",err,error,*999)
#endif    
    
    DO i=1,SIZE(B,3)
      DO j=1,SIZE(B,4)
        C(i,j) = 0.0_DP
        DO k=1,SIZE(B,1)
          DO l=1,SIZE(B,2)
            C(i,j) = C(i,j) + A(k,l)*B(k,l,i,j)
          ENDDO !l
        ENDDO !k
      ENDDO !j
    ENDDO !i
    
    EXITS("DoubleDotProductMatrix2Matrix4DP")
    RETURN
999 ERRORSEXITS("DoubleDotProductMatrix2Matrix4DP",err,error)
    RETURN 1
    
  END SUBROUTINE DoubleDotProductMatrix2Matrix4DP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the matrix4-matrix2 double dot product of the single precision matrices A:B in C.
  SUBROUTINE DoubleDotProductMatrix4Matrix2SP(A,B,C,err,error,*)

    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:,:,:) !<The A matrix
    REAL(SP), INTENT(IN) :: B(:,:) !<The B matrix
    REAL(SP), INTENT(OUT) :: C(:,:) !<On exit, the matrix double dot product C=A:B
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,k,l
    
    ENTERS("DoubleDotProductMatrix4Matrix2SP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,3)/=SIZE(B,1).OR.SIZE(A,4)/=SIZE(B,2)) &
      & CALL FlagError("The last two dimensions of the A matrix does not conform to the number of dimensions of the B matrix.", &
      & err,error,*999)
    IF(SIZE(C,1)/=SIZE(A,1).OR.SIZE(C,2)/=SIZE(A,2))  &
      & CALL FlagError("The C matrix is not the same size as the first two dimensions of the A matrix.",err,error,*999)
#endif    
    
    DO i=1,SIZE(A,1)
      DO j=1,SIZE(A,2)
        C(i,j) = 0.0_SP
        DO k=1,SIZE(A,3)
          DO l=1,SIZE(A,4)
            C(i,j) = C(i,j) + A(i,j,k,l)*B(k,l)
          ENDDO !l
        ENDDO !k
      ENDDO !j
    ENDDO !i
    
    EXITS("DoubleDotProductMatrix4Matrix2SP")
    RETURN
999 ERRORSEXITS("DoubleDotProductMatrix4Matrix2SP",err,error)
    RETURN 1
    
  END SUBROUTINE DoubleDotProductMatrix4Matrix2SP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the matrix4-matrix2 double dot product of the double precision matrices A:B in C.
  SUBROUTINE DoubleDotProductMatrix4Matrix2DP(A,B,C,err,error,*)

    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:,:,:) !<The A matrix
    REAL(DP), INTENT(IN) :: B(:,:) !<The B matrix
    REAL(DP), INTENT(OUT) :: C(:,:) !<On exit, the matrix double dot product C=A:B
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,k,l
    
    ENTERS("DoubleDotProductMatrix4Matrix2DP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,3)/=SIZE(B,1).OR.SIZE(A,4)/=SIZE(B,2)) &
      & CALL FlagError("The last two dimensions of the A matrix does not conform to the number of dimensions of the B matrix.", &
      & err,error,*999)
    IF(SIZE(C,1)/=SIZE(A,1).OR.SIZE(C,2)/=SIZE(A,2))  &
      & CALL FlagError("The C matrix is not the same size as the first two dimensions of the A matrix.",err,error,*999)
#endif    
    
    DO i=1,SIZE(A,1)
      DO j=1,SIZE(A,2)
        C(i,j) = 0.0_DP
        DO k=1,SIZE(A,3)
          DO l=1,SIZE(A,4)
            C(i,j) = C(i,j) + A(i,j,k,l)*B(k,l)
          ENDDO !l
        ENDDO !k
      ENDDO !j
    ENDDO !i
    
    EXITS("DoubleDotProductMatrix4Matrix2DP")
    RETURN
999 ERRORSEXITS("DoubleDotProductMatrix4Matrix2DP",err,error)
    RETURN 1
    
  END SUBROUTINE DoubleDotProductMatrix4Matrix2DP
  
  !
  !================================================================================================================================
  !

  !>Calculates the elliptic integral of the second kind - E(m), for a double precision argument.
  PURE FUNCTION EdpDP(x)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: x !<The value to evaluate the function at
    !Function variable
    REAL(DP) :: EdpDP
    !Local variables
    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    REAL(DP), PARAMETER :: a1=0.44325141463_DP
    REAL(DP), PARAMETER :: a2=0.06260601220_DP
    REAL(DP), PARAMETER :: a3=0.04757383546_DP
    REAL(DP), PARAMETER :: a4=0.01736506451_DP
    REAL(DP), PARAMETER :: b1=0.24998368310_DP
    REAL(DP), PARAMETER :: b2=0.09200180037_DP
    REAL(DP), PARAMETER :: b3=0.04069697526_DP
    REAL(DP), PARAMETER :: b4=0.00526449639_DP
    REAL(DP) :: term1,term2,x1
    
    x1=1.0_DP-x
    term1=1.0_DP+(a1+(a2+(a3+a4*x1)*x1)*x1)*x1
    term2=(b1+(b2+(b3+b4*x1)*x1)*x1)*x1
    EdpDP=term1+term2*LOG(1.0_DP/x1)
    
    RETURN
    
  END FUNCTION EdpDP
  
  !
  !================================================================================================================================
  !

  !>Calculates the elliptic integral of the second kind - E(m), for a single precision argument.
  PURE FUNCTION EdpSP(x)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: x !<The value to evaluate the function at
    !Function variable
    REAL(SP) :: EdpSP
    !Local variables
    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    REAL(SP), PARAMETER :: a1=0.44325141463_SP
    REAL(SP), PARAMETER :: a2=0.06260601220_SP
    REAL(SP), PARAMETER :: a3=0.04757383546_SP
    REAL(SP), PARAMETER :: a4=0.01736506451_SP
    REAL(SP), PARAMETER :: b1=0.24998368310_SP
    REAL(SP), PARAMETER :: b2=0.09200180037_SP
    REAL(SP), PARAMETER :: b3=0.04069697526_SP
    REAL(SP), PARAMETER :: b4=0.00526449639_SP
    REAL(SP) :: term1,term2,x1
    
    x1=1.0_SP-x
    term1=1.0_SP+(a1+(a2+(a3+a4*x1)*x1)*x1)*x1
    term2=(b1+(b2+(b3+b4*x1)*x1)*x1)*x1
    EdpSP=term1+term2*LOG(1.0_SP/x1)
        
    RETURN
    
  END FUNCTION EdpSP
  
  !
  !================================================================================================================================
  !

  !>Returns the eigenvalues of a full single precision matrix A.
  SUBROUTINE EigenvalueFullSP(A,eValues,err,error,*)
    
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The matrix to find the eignenvalues for
    REAL(SP), INTENT(OUT) :: eValues(:) !<On exit, the eignevalues
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i
    REAL(SP) :: angle,b2,b3,c1,c2,d,q,q3,r,ri1,ri2,ri3,ri4,rq,temp,theta
    
    ENTERS("EigenvalueFullSP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix is not square.",err,error,*999)    
    IF(SIZE(A,1)>SIZE(eValues,1)) CALL FlagError("Evalues is too small.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1))
    CASE(1)
      eValues(1)=A(1,1)
    CASE(2)
      IF(ABS(A(1,2))>ZERO_TOLERANCE_SP) THEN
        ri1=A(1,1)+A(2,2)
        ri2=A(1,1)*A(2,2)-A(1,2)**2
        b2=ri1/2.0_SP
        c1=ri1*ri1
        c2=4.0_SP*ri2
        IF(c2>c1) CALL FlagError("Complex roots found in quadratic equation.",err,error,*999)
        b3=SQRT(c1-c2)/2.0_SP
        eValues(1)=b2+b3
        eValues(2)=b2-b3
      ELSE
        eValues(1)=A(1,1)
        eValues(2)=A(2,2)
      ENDIF
      IF(ABS(eValues(2))>ABS(eValues(1))) THEN
        temp=eValues(1)
        eValues(1)=eValues(2)
        eValues(2)=temp
      ENDIF
    CASE(3)
      ri1=A(1,1)+A(2,2)+A(3,3)
      ri2=A(1,1)*A(2,2)+A(2,2)*A(3,3)+A(3,3)*A(1,1)-(A(1,2)**2+A(2,3)**2+A(3,1)**2)
      CALL Determinant(A,ri3,err,error,*999)
      ri4=ri1/3.0_SP
      q=ri4*ri4-ri2/3.0_SP   
      r=ri4*(ri4*ri4-ri2/2.0_SP)+ri3/2.0_SP
      q3=q*q*q
      d=r*r-q3
      IF(ABS(d)>ZERO_TOLERANCE_SP) CALL FlagError("Complex roots found in solution of cubic equation.",err,error,*999)
      rq=SQRT(ABS(q))
      IF(ABS(q)<ZERO_TOLERANCE_SP) THEN
        theta=0.0_SP
      ELSE
        theta=ACOS(r/SQRT(ABS(q3)))/3.0_SP
      ENDIF
      angle=2.0_SP*REAL(PI,SP)/3.0_SP
      eValues(1)=2.0_SP*rq*COS(theta)+ri4
      eValues(2)=2.0_SP*rq*COS(theta+angle)+ri4
      eValues(3)=2.0_SP*rq*COS(theta+2.0_SP*angle)+ri4
      DO i=1,2
        IF(ABS(eValues(3))>ABS(eValues(i))) THEN
          temp=eValues(i)
          eValues(i)=eValues(3)
          eValues(3)=temp
        ENDIF
      ENDDO !i
    CASE DEFAULT
      CALL FlagError("Matrix size not implemented.",err,error,*999)
    END SELECT

    EXITS("EigenvalueFullSP")
    RETURN
999 ERRORSEXITS("EigenvalueFullSP",err,error)
    RETURN 1
    
  END SUBROUTINE EigenvalueFullSP

  !
  !================================================================================================================================
  !

  !>Returns the eigenvalues of a full double precision matrix A.
  SUBROUTINE EigenvalueFullDP(A,eValues,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The matrix to find the eigenvalues of
    REAL(DP), INTENT(OUT) :: eValues(:) !<On exit, the eigenvalues of the matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i
    REAL(DP) :: angle,b2,b3,c1,c2,d,q,q3,r,ri1,ri2,ri3,ri4,rq,temp,theta
    
    ENTERS("EigenvalueFullDP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix is not square.",err,error,*999)    
    IF(SIZE(A,1)>SIZE(eValues,1)) CALL FlagError("Evalues is too small.",err,error,*999)
#endif
    
    SELECT CASE(SIZE(A,1))
    CASE(1)
      eValues(1)=A(1,1)
    CASE(2)
      IF(ABS(A(1,2))>ZERO_TOLERANCE_DP) THEN
        ri1=A(1,1)+A(2,2)
        ri2=A(1,1)*A(2,2)-A(1,2)**2
        b2=ri1/2.0_DP
        c1=ri1*ri1
        c2=4.0_DP*ri2
        IF(c2>c1) CALL FlagError("Complex roots found in quadratic equation.",err,error,*999)
        b3=SQRT(c1-c2)/2.0_DP
        eValues(1)=b2+b3
        eValues(2)=b2-b3
      ELSE
        eValues(1)=A(1,1)
        eValues(2)=A(2,2)
      ENDIF
      IF(ABS(eValues(2))>ABS(eValues(1))) THEN
        temp=eValues(1)
        eValues(1)=eValues(2)
        eValues(2)=temp
      ENDIF
    CASE(3)
      ri1=A(1,1)+A(2,2)+A(3,3)
      ri2=A(1,1)*A(2,2)+A(2,2)*A(3,3)+A(3,3)*A(1,1)-(A(1,2)**2+A(2,3)**2+A(3,1)**2)
      CALL Determinant(A,ri3,err,error,*999)
      ri4=ri1/3.0_DP
      q=ri4*ri4-ri2/3.0_DP   
      r=ri4*(ri4*ri4-ri2/2.0_DP)+ri3/2.0_DP
      q3=q*q*q
      d=r*r-q3
      IF(ABS(d)>ZERO_TOLERANCE_DP) CALL FlagError("Complex roots found in solution of cubic equation.",err,error,*999)
      rq=SQRT(ABS(q))
      IF(ABS(q)<ZERO_TOLERANCE_DP) THEN
        theta=0.0_DP
      ELSE
        theta=ACOS(r/SQRT(ABS(q3)))/3.0_DP
      ENDIF
      angle=2.0_DP*PI/3.0_DP
      eValues(1)=2.0_DP*rq*COS(theta)+ri4
      eValues(2)=2.0_DP*rq*COS(theta+angle)+ri4
      eValues(3)=2.0_DP*rq*COS(theta+2.0_DP*angle)+ri4
      DO i=1,2
        IF(ABS(eValues(3))>ABS(eValues(i))) THEN
          temp=eValues(i)
          eValues(i)=eValues(3)
          eValues(3)=temp
        ENDIF
      ENDDO !i
    CASE DEFAULT
      CALL FlagError("Matrix size not implemented.",err,error,*999)
    END SELECT

    EXITS("EigenvalueFullDP")
    RETURN
999 ERRORSEXITS("EigenvalueFullDP",err,error)
    RETURN 1
    
  END SUBROUTINE EigenvalueFullDP

  !
  !================================================================================================================================
  !

  !>Returns the normalised eigenvector of a full single precision symmetric matrix A that corresponds to the eigenvalue eValue. 
  SUBROUTINE EigenvectorFullSP(A,eValue,eVector,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The matrix to find the eignevectors for
    REAL(SP), INTENT(IN) :: eValue !<The eigenvalue to find the eignevector for
    REAL(SP), INTENT(OUT) :: eVector(:) !<On exit, the eigenvector corresponding the the eigenvalue
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,i1,i2,i3,iCycle(3,3)
    REAL(SP) :: al,b(SIZE(A,1)),sum,U(SIZE(A,1),SIZE(A,2)),x(SIZE(A,1))

    DATA iCycle /2,1,1,3,1,2,1,2,1/
    
    ENTERS("EigenvectorFullSP",err,error,*999)

!!THIS NEEDS TO BE CHECKED

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix is not square.",err,error,*999)
    IF(SIZE(A,1)>SIZE(eVector,1)) CALL FlagError("Evector is too small.",err,error,*999)
#endif
    
    SELECT CASE(SIZE(A,1))
    CASE(1)
      eVector(1)=1.0_SP
    CASE(2)
      IF(ABS(A(1,2))>ZERO_TOLERANCE_SP) THEN
        IF(ABS(A(1,1)-eValue)>ABS(A(2,2)-eValue)) THEN
          al=SQRT(A(1,2)**2+(A(1,1)-eValue)**2)
          eVector(1)=A(1,2)/al
          eVector(2)=(eValue-A(1,1))/al
        ELSE
          al=SQRT(A(1,2)**2+(A(2,2)-eValue)**2)
          eVector(1)=(eValue-A(2,2))/al
          eVector(2)=A(1,2)/al
        ENDIF
      ELSE IF(ABS(eValue-A(1,1))<ZERO_TOLERANCE_SP) THEN
        eVector(1)=1.0_SP
        eVector(2)=0.0_SP
      ELSE IF(ABS(eValue-A(2,2))<ZERO_TOLERANCE_SP) THEN
        eVector(1)=0.0_SP
        eVector(2)=1.0_DP
      ENDIF
    CASE(3)
      IF(ABS(A(1,2))<ZERO_TOLERANCE_SP.AND.ABS(A(1,3))<ZERO_TOLERANCE_SP.AND.ABS(A(2,3))<ZERO_TOLERANCE_SP) THEN
        eVector=0.0_SP
        CALL FlagError("Zero matrix?? Eigenvectors undetermined.",err,error,*999)
      ELSE
        DO i=1,3
          U(i,:)=A(i,:)
          U(i,i)=U(i,i)-eValue
        ENDDO !i
        DO i=1,3
          x(i)=1.0_SP
          i1=iCycle(i,1)
          i2=iCycle(i,2)
          i3=iCycle(i,3)
          b(1)=-1.0_SP*U(i1,i)
          b(2)=-1.0_SP*U(i2,i)
          CALL SolveSmallLinearSystem(U(i1:i2:i3,i1:i2:i3),x,b,err,error,*999)
          sum=DOT_PRODUCT(U(i,:),x)
          IF(ABS(sum)<ZERO_TOLERANCE_SP) THEN
            CALL Normalise(x,eVector,err,error,*999)
          ENDIF
        ENDDO !i
      ENDIF
    CASE DEFAULT
      CALL FlagError("Matrix size not implemented.",err,error,*999)
    END SELECT

    EXITS("EigenvectorFullSP")
    RETURN
999 ERRORSEXITS("EigenvectorFullSP",err,error)
    RETURN 1
    
  END SUBROUTINE EigenvectorFullSP

  !
  !================================================================================================================================
  !

  !>Returns the normalised eigenvector of a full double precision symmetric matrix A that corresponds to the eigenvalue eValue.
  SUBROUTINE EigenvectorFullDP(A,eValue,eVector,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The matrix to find the eignevectors for
    REAL(DP), INTENT(IN) :: eValue !<The eigenvalue to find the eignevector for
    REAL(DP), INTENT(OUT) :: eVector(:) !<On exit, the eigenvector corresponding the the eigenvalue
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,i1,i2,i3,iCycle(3,3)
    REAL(DP) :: al,b(SIZE(A,1)),sum,U(SIZE(A,1),SIZE(A,2)),x(SIZE(A,1))

    DATA iCycle /2,1,1,3,1,2,1,2,1/
    
    ENTERS("EigenvectorFullDP",err,error,*999)

!!THIS NEEDS TO BE CHECKED

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix is not square.",err,error,*999)
    IF(SIZE(A,1)>SIZE(eVector,1)) CALL FlagError("Evector is too small.",err,error,*999)
#endif
    
    SELECT CASE(SIZE(A,1))
    CASE(1)
      eVector(1)=1.0_DP
    CASE(2)
      IF(ABS(A(1,2))>ZERO_TOLERANCE_DP) THEN
        IF(ABS(A(1,1)-eValue)>ABS(A(2,2)-eValue)) THEN
          al=SQRT(A(1,2)**2+(A(1,1)-eValue)**2)
          eVector(1)=A(1,2)/al
          eVector(2)=(eValue-A(1,1))/al
        ELSE
          al=SQRT(A(1,2)**2+(A(2,2)-eValue)**2)
          eVector(1)=(eValue-A(2,2))/al
          eVector(2)=A(1,2)/al
        ENDIF
      ELSE IF(ABS(eValue-A(1,1))<ZERO_TOLERANCE_DP) THEN
        eVector(1)=1.0_DP
        eVector(2)=0.0_DP
      ELSE IF(ABS(eValue-A(2,2))<ZERO_TOLERANCE_DP) THEN
        eVector(1)=0.0_DP
        eVector(2)=1.0_DP
      ENDIF
    CASE(3)
      IF(ABS(A(1,2))<ZERO_TOLERANCE_DP.AND.ABS(A(1,3))<ZERO_TOLERANCE_DP.AND.ABS(A(2,3))<ZERO_TOLERANCE_DP) THEN
        eVector=0.0_DP
        CALL FlagError("Zero matrix?? Eigenvectors undetermined.",err,error,*999)
      ELSE
        DO i=1,3
          U(i,:)=A(i,:)
          U(i,i)=U(i,i)-eValue
        ENDDO !i
        DO i=1,3
          x(i)=1.0_DP
          i1=iCycle(i,1)
          i2=iCycle(i,2)
          i3=iCycle(i,3)
          b(1)=-1.0_DP*U(i1,i)
          b(2)=-1.0_DP*U(i2,i)
          CALL SolveSmallLinearSystem(U(i1:i2:i3,i1:i2:i3),x,b,err,error,*999)
          sum=DOT_PRODUCT(U(i,:),x)
          IF(ABS(sum)<ZERO_TOLERANCE_DP) THEN
            CALL Normalise(x,eVector,err,error,*999)
            IF(err /= 0) GOTO 999
          ENDIF
        ENDDO !i
      ENDIF
    CASE DEFAULT
      CALL FlagError("Matrix size not implemented.",err,error,*999)
    END SELECT

    EXITS("EigenvectorFullDP")
    RETURN
999 ERRORSEXITS("EigenvectorFullDP",err,error)
    RETURN 1
    
  END SUBROUTINE EigenvectorFullDP
  
  !
  !================================================================================================================================
  !
  
  !>Calculates the modified Bessel function of the first kind of order 0 using the approximation of Abromowitz and Stegun,
  !>for a double precision argument.
  PURE FUNCTION I0DP(x)
      
    !Argument variables
    REAL(DP), INTENT(IN) :: x !<The value to evaluate the function at
    !Function variable
    REAL(DP) :: I0DP
    !Local variables
    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    REAL(DP), PARAMETER :: a1=1.0_DP
    REAL(DP), PARAMETER :: a2=3.5156229_DP
    REAL(DP), PARAMETER :: a3=3.0899424_DP
    REAL(DP), PARAMETER :: a4=1.2067492_DP
    REAL(DP), PARAMETER :: a5=0.2659732_DP
    REAL(DP), PARAMETER :: a6=0.0360768_DP
    REAL(DP), PARAMETER :: a7=0.0045813_DP
    REAL(DP) :: t

    !Calculate I0(x) for x < 3.75
    t=x*x/(3.75_DP*3.75_DP)
    I0DP=a1+(a2+(a3+(a4+(a5+(a6+a7*t)*t)*t)*t)*t)*t
    
    RETURN
    
  END FUNCTION I0DP
  
  !
  !================================================================================================================================
  !
  
  !>Calculates the modified Bessel function of the first kind of order 0 using the approximation of Abromowitz and Stegun,
  !>for a single precision argument.
  PURE FUNCTION I0SP(x)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: x !<The value to evaluate the function at
    !Function variable
    REAL(SP) :: I0SP
    !Local variables
    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    REAL(SP), PARAMETER :: a1=1.0_SP
    REAL(SP), PARAMETER :: a2=3.5156229_SP
    REAL(SP), PARAMETER :: a3=3.0899424_SP
    REAL(SP), PARAMETER :: a4=1.2067492_SP
    REAL(SP), PARAMETER :: a5=0.2659732_SP
    REAL(SP), PARAMETER :: a6=0.0360768_SP
    REAL(SP), PARAMETER :: a7=0.0045813_SP
    REAL(SP) :: t

    !Calculate I0(x) for x < 3.75
    t=x*x/(3.75_SP*3.75_SP)
    I0SP=a1+(a2+(a3+(a4+(a5+(a6+a7*t)*t)*t)*t)*t)*t
    
    RETURN
    
  END FUNCTION I0SP
  
  !
  !================================================================================================================================
  !
  
  !>Calculates the modified Bessel function of the first kind of order 1 using the approximation of Abromowitz and Stegun,
  !>for a double precision argument.
  PURE FUNCTION I1DP(x)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: x !<The value to evaluate the function at
    !Function variable
    REAL(DP) :: I1DP
    !Local variables
    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    REAL(DP), PARAMETER :: a1=0.5_DP
    REAL(DP), PARAMETER :: a2=0.87890594_DP
    REAL(DP), PARAMETER :: a3=0.51498869_DP
    REAL(DP), PARAMETER :: a4=0.15084934_DP
    REAL(DP), PARAMETER :: a5=0.02658733_DP
    REAL(DP), PARAMETER :: a6=0.00301532_DP
    REAL(DP), PARAMETER :: a7=0.00032411_DP
    REAL(DP) :: t

    !Calculate I1(x)
    t=(x/3.75_DP)*(x/3.75_DP)
    I1DP=(a1+(a2+(a3+(a4+(a5+(a6+a7*t)*t)*t)*t)*t)*t)*x
    
    RETURN
    
  END FUNCTION I1DP
  
  !
  !================================================================================================================================
  !
  
  !>Calculates the modified Bessel function of the first kind of order 1 using the approximation of Abromowitz and Stegun,
  !>for a single precision argument.
  PURE FUNCTION I1SP(x)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: x !<The value to evaluate the function at
    !Function variable
    REAL(SP) :: I1SP
    !Local variables
    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    REAL(SP), PARAMETER :: a1=0.5_SP
    REAL(SP), PARAMETER :: a2=0.87890594_SP
    REAL(SP), PARAMETER :: a3=0.51498869_SP
    REAL(SP), PARAMETER :: a4=0.15084934_SP
    REAL(SP), PARAMETER :: a5=0.02658733_SP
    REAL(SP), PARAMETER :: a6=0.00301532_SP
    REAL(SP), PARAMETER :: a7=0.00032411_SP
    REAL(SP) :: t

    !Calculate I1(x)
    t=(x/3.75_SP)*(x/3.75_SP)
    I1SP=(a1+(a2+(a3+(a4+(a5+(a6+7*t)*t)*t)*t)*t)*t)*x
    
    RETURN
    
  END FUNCTION I1SP
  
  !
  !================================================================================================================================
  !

  !>Returns an identity matrix
  SUBROUTINE IdentityMatrixSP(A,err,error,*)
    
    !Argument variables
    REAL(SP), INTENT(OUT) :: A(:,:) !<On exit, the identity matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i
    
    ENTERS("IdentityMatrixSP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix A is not square.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1)) 
    CASE(1)
      A(1,1)=1.0_SP
    CASE(2)
      A(1,1)=1.0_SP
      A(2,1)=0.0_SP
      A(1,2)=0.0_SP
      A(2,2)=1.0_SP
    CASE(3)
      A(1,1)=1.0_SP
      A(2,1)=0.0_SP
      A(3,1)=0.0_SP
      A(1,2)=0.0_SP
      A(2,2)=1.0_SP
      A(3,2)=0.0_SP
      A(1,3)=0.0_SP
      A(2,3)=0.0_SP
      A(3,3)=1.0_SP
    CASE DEFAULT
      A=0.0_DP
      DO i=1,SIZE(A,1)
        A(i,i)=1.0_SP
      ENDDO !i
    END SELECT

    EXITS("IdentityMatrixSP")
    RETURN
999 ERRORSEXITS("IdentityMatrixSP",err,error)
    RETURN 1
    
  END SUBROUTINE IdentityMatrixSP

  !
  !================================================================================================================================
  !

  !>Returns an identity matrix
  SUBROUTINE IdentityMatrixDP(A,err,error,*)
    
    !Argument variables
    REAL(DP), INTENT(OUT) :: A(:,:) !<On exit, the identity matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i
    
    ENTERS("IdentityMatrixDP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix A is not square.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1)) 
    CASE(1)
      A(1,1)=1.0_DP
    CASE(2)
      A(1,1)=1.0_DP
      A(2,1)=0.0_DP
      A(1,2)=0.0_DP
      A(2,2)=1.0_DP
    CASE(3)
      A(1,1)=1.0_DP
      A(2,1)=0.0_DP
      A(3,1)=0.0_DP
      A(1,2)=0.0_DP
      A(2,2)=1.0_DP
      A(3,2)=0.0_DP
      A(1,3)=0.0_DP
      A(2,3)=0.0_DP
      A(3,3)=1.0_DP
    CASE DEFAULT
      A=0.0_DP
      DO i=1,SIZE(A,1)
        A(i,i)=1.0_DP
      ENDDO !i
    END SELECT

    EXITS("IdentityMatrixDP")
    RETURN
999 ERRORSEXITS("IdentityMatrixDP",err,error)
    RETURN 1
    
  END SUBROUTINE IdentityMatrixDP

  !
  !================================================================================================================================
  !

  !Returns the first invariant of a single precision second order matrix A
  SUBROUTINE InvariantsOne2SP(A,I1,err,error,*)
    
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The second order A matrix to get the invariants for
    REAL(SP), INTENT(OUT) :: I1 !<On exit, the first invariant of A i.e., I1 = tr A
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("InvariantsOne2SP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix A is not square.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1)) 
    CASE(1)
      I1=A(1,1)
    CASE(2)
      I1=A(1,1)+A(2,2)
    CASE(3)
      I1=A(1,1)+A(2,2)+A(3,3)
    CASE DEFAULT
      CALL FlagError("Matrix size is not implemented.",err,error,*999)
    END SELECT

    EXITS("InvariantsOne2SP")
    RETURN
999 ERRORSEXITS("InvariantsOne2SP",err,error)
    RETURN 1
    
  END SUBROUTINE InvariantsOne2SP

  !
  !================================================================================================================================
  !

  !Returns the first invariant of a double precision second order matrix A
  SUBROUTINE InvariantsOne2DP(A,I1,err,error,*)
    
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The second order A matrix to get the invariants for
    REAL(DP), INTENT(OUT) :: I1 !<On exit, the first invariant of A i.e., I1 = tr A
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("InvariantsOne2DP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix A is not square.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1)) 
    CASE(1)
      I1=A(1,1)
    CASE(2)
      I1=A(1,1)+A(2,2)
    CASE(3)
      I1=A(1,1)+A(2,2)+A(3,3)
    CASE DEFAULT
      CALL FlagError("Matrix size is not implemented.",err,error,*999)
    END SELECT

    EXITS("InvariantsOne2DP")
    RETURN
999 ERRORSEXITS("InvariantsOne2DP",err,error)
    RETURN 1
    
  END SUBROUTINE InvariantsOne2DP

  !
  !================================================================================================================================
  !

  !Returns three invariants of a single precision second order matrix A
  SUBROUTINE InvariantsThree2SP(A,I1,I2,I3,err,error,*)
    
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The second order A matrix to get the invariants for
    REAL(SP), INTENT(OUT) :: I1 !<On exit, the first invariant of A i.e., I1 = tr A
    REAL(SP), INTENT(OUT) :: I2 !<On exit, the second invariant of A i.e., I2 = 1/2.((tr A)^2-tr(A.A))
    REAL(SP), INTENT(OUT) :: I3 !<On exit, the third invariant of A i.e., I3 = det A
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("InvariantsThree2SP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix A is not square.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1)) 
    CASE(1)
      I1=A(1,1)
      I2=0.0_SP
      I3=A(1,1)
    CASE(2)
      I1=A(1,1)+A(2,2)
      I2=A(1,1)*A(2,2)-A(1,2)*A(2,1)
      I3=A(1,1)*A(2,2)-A(1,2)*A(2,1)
    CASE(3)
      I1=A(1,1)+A(2,2)+A(3,3)
      I2=A(1,1)*A(2,2)+A(1,1)*A(3,3)+A(2,2)*A(3,3)-A(1,2)*A(2,1)-A(1,3)*A(3,1)-A(2,3)*A(3,2)
      I3=A(1,1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+A(1,3)*A(2,1)*A(3,2)-A(1,1)*A(3,2)*A(2,3)-A(2,1)*A(1,2)*A(3,3)- &
        & A(3,1)*A(2,2)*A(1,3)
    CASE DEFAULT
      CALL FlagError("Matrix size is not implemented.",err,error,*999)
    END SELECT

    EXITS("InvariantsThree2SP")
    RETURN
999 ERRORSEXITS("InvariantsThree2SP",err,error)
    RETURN 1
    
  END SUBROUTINE InvariantsThree2SP

  !
  !================================================================================================================================
  !

  !Returns three invariants of a double precision second order matrix A
  SUBROUTINE InvariantsThree2DP(A,I1,I2,I3,err,error,*)
    
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The second order A matrix to get the invariants for
    REAL(DP), INTENT(OUT) :: I1 !<On exit, the first invariant of A i.e., I1 = tr A
    REAL(DP), INTENT(OUT) :: I2 !<On exit, the second invariant of A i.e., I2 = 1/2.((tr A)^2-tr(A.A))
    REAL(DP), INTENT(OUT) :: I3 !<On exit, the third invariant of A i.e., I3 = det A
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("InvariantsThree2DP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix A is not square.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1)) 
    CASE(1)
      I1=A(1,1)
      I2=0.0_DP
      I3=A(1,1)
    CASE(2)
      I1=A(1,1)+A(2,2)
      I2=A(1,1)*A(2,2)-A(1,2)*A(2,1)
      I3=A(1,1)*A(2,2)-A(1,2)*A(2,1)
    CASE(3)
      I1=A(1,1)+A(2,2)+A(3,3)
      I2=A(1,1)*A(2,2)+A(1,1)*A(3,3)+A(2,2)*A(3,3)-A(1,2)*A(2,1)-A(1,3)*A(3,1)-A(2,3)*A(3,2)
      I3=A(1,1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+A(1,3)*A(2,1)*A(3,2)-A(1,1)*A(3,2)*A(2,3)-A(2,1)*A(1,2)*A(3,3)- &
        & A(3,1)*A(2,2)*A(1,3)
    CASE DEFAULT
      CALL FlagError("Matrix size is not implemented.",err,error,*999)
    END SELECT

    EXITS("InvariantsThree2DP")
    RETURN
999 ERRORSEXITS("InvariantsThree2DP",err,error)
    RETURN 1
    
  END SUBROUTINE InvariantsThree2DP

  !
  !================================================================================================================================
  !

  !Returns the first two invariants of a single precision second order matrix A
  SUBROUTINE InvariantsTwo2SP(A,I1,I2,err,error,*)
    
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The second order A matrix to get the invariants for
    REAL(SP), INTENT(OUT) :: I1 !<On exit, the first invariant of A i.e., I1 = tr A
    REAL(SP), INTENT(OUT) :: I2 !<On exit, the second invariant of A i.e., I2 = 1/2.((tr A)^2-tr(A.A))
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("InvariantsTwo2SP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix A is not square.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1)) 
    CASE(1)
      I1=A(1,1)
      I2=0.0_SP
    CASE(2)
      I1=A(1,1)+A(2,2)
      I2=A(1,1)*A(2,2)-A(1,2)*A(2,1)
    CASE(3)
      I1=A(1,1)+A(2,2)+A(3,3)
      I2=A(1,1)*A(2,2)+A(1,1)*A(3,3)+A(2,2)*A(3,3)-A(1,2)*A(2,1)-A(1,3)*A(3,1)-A(2,3)*A(3,2)
    CASE DEFAULT
      CALL FlagError("Matrix size is not implemented.",err,error,*999)
    END SELECT

    EXITS("InvariantsTwo2SP")
    RETURN
999 ERRORSEXITS("InvariantsTwo2SP",err,error)
    RETURN 1
    
  END SUBROUTINE InvariantsTwo2SP

  !
  !================================================================================================================================
  !

  !Returns the first two invariants of a double precision second order matrix A
  SUBROUTINE InvariantsTwo2DP(A,I1,I2,err,error,*)
    
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The second order A matrix to get the invariants for
    REAL(DP), INTENT(OUT) :: I1 !<On exit, the first invariant of A i.e., I1 = tr A
    REAL(DP), INTENT(OUT) :: I2 !<On exit, the second invariant of A i.e., I2 = 1/2.((tr A)^2-tr(A.A))
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("InvariantsTwo2DP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix A is not square.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1)) 
    CASE(1)
      I1=A(1,1)
      I2=0.0_DP
    CASE(2)
      I1=A(1,1)+A(2,2)
      I2=A(1,1)*A(2,2)-A(1,2)*A(2,1)
    CASE(3)
      I1=A(1,1)+A(2,2)+A(3,3)
      I2=A(1,1)*A(2,2)+A(1,1)*A(3,3)+A(2,2)*A(3,3)-A(1,2)*A(2,1)-A(1,3)*A(3,1)-A(2,3)*A(3,2)
    CASE DEFAULT
      CALL FlagError("Matrix size is not implemented.",err,error,*999)
    END SELECT

    EXITS("InvariantsTwo2DP")
    RETURN
999 ERRORSEXITS("InvariantsTwo2DP",err,error)
    RETURN 1
    
  END SUBROUTINE InvariantsTwo2DP

  !
  !================================================================================================================================
  !

  !>Inverts a full single precision matrix A to give matrix B and returns the determinant of A in det.
  SUBROUTINE InvertFullSP(A,B,det,err,error,*)
    
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The A matrix to invert
    REAL(SP), INTENT(OUT) :: B(:,:) !<On exit, the inverse of A
    REAL(SP), INTENT(OUT) :: det !<On exit, the determinant of A
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("InvertFullSP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix A is not square.",err,error,*999)
    IF(SIZE(B,1)/=SIZE(A,1).OR.SIZE(B,2)/=SIZE(A,2)) CALL FlagError("Matrix B is not the same size as matrix A.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1)) 
    CASE(1)
      det=A(1,1)
      IF(ABS(det)>ZERO_TOLERANCE_SP) THEN
        B(1,1)=1.0_SP/A(1,1)
      ELSE
        CALL FlagWarning("Matrix A is zero and cannot be inverted.",err,error,*999)
        B(1,1)=0.0_SP
      ENDIF
    CASE(2)
      det=A(1,1)*A(2,2)-A(1,2)*A(2,1)
      IF(ABS(det)>ZERO_TOLERANCE_SP) THEN
        B(1,1)=A(2,2)/det
        B(1,2)=-A(1,2)/det
        B(2,1)=-A(2,1)/det
        B(2,2)=A(1,1)/det
      ELSE
        CALL FlagWarning("Zero Determinant for matrix A.",err,error,*999)
        B=0.0_SP
      ENDIF
    CASE(3)
      det=A(1,1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+A(1,3)*A(2,1)*A(3,2)-A(1,1)*A(3,2)*A(2,3)-A(2,1)*A(1,2)*A(3,3)- &
        & A(3,1)*A(2,2)*A(1,3)
      IF(ABS(det)>ZERO_TOLERANCE_SP) THEN
        B(1,1)=(A(2,2)*A(3,3)-A(3,2)*A(2,3))/det
        B(2,1)=(A(2,3)*A(3,1)-A(3,3)*A(2,1))/det
        B(3,1)=(A(2,1)*A(3,2)-A(3,1)*A(2,2))/det
        B(1,2)=(A(3,2)*A(1,3)-A(1,2)*A(3,3))/det
        B(2,2)=(A(3,3)*A(1,1)-A(1,3)*A(3,1))/det
        B(3,2)=(A(3,1)*A(1,2)-A(1,1)*A(3,2))/det
        B(1,3)=(A(1,2)*A(2,3)-A(2,2)*A(1,3))/det
        B(2,3)=(A(1,3)*A(2,1)-A(2,3)*A(1,1))/det
        B(3,3)=(A(1,1)*A(2,2)-A(2,1)*A(1,2))/det
      ELSE
        CALL FlagWarning("Zero Determinant for matrix A.",err,error,*999)
        B=0.0_SP
      ENDIF
    CASE DEFAULT
      CALL FlagError("Matrix size is not implemented.",err,error,*999)
    END SELECT

    EXITS("InvertFullSP")
    RETURN
999 ERRORSEXITS("InvertFullSP",err,error)
    RETURN 1
    
  END SUBROUTINE InvertFullSP

  !
  !================================================================================================================================
  !

  !>Inverts a full double precision matrix A to give matrix B and returns the determinant of A in det.
  SUBROUTINE InvertFullDP(A,B,det,err,error,*)
    
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The matrix A to invert
    REAL(DP), INTENT(OUT) :: B(:,:) !<On exit, the inverse of A
    REAL(DP), INTENT(OUT) :: det !<On exit, the determinant of A
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("InvertFullDP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix A is not square.",err,error,*999)
    IF(SIZE(B,1)/=SIZE(A,1).OR.SIZE(B,2)/=SIZE(A,2)) CALL FlagError("Matrix B is not the same size as matrix A.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1)) 
    CASE(1)
      det=A(1,1)
      IF(ABS(det)>ZERO_TOLERANCE_DP) THEN
        B(1,1)=1.0_DP/A(1,1)
      ELSE
        CALL FlagWarning("Matrix A is zero and cannot be inverted.",err,error,*999)
        B(1,1)=0.0_DP
      ENDIF
    CASE(2)
      det=A(1,1)*A(2,2)-A(2,1)*A(1,2)
      IF(ABS(det)>ZERO_TOLERANCE_DP) THEN
        B(1,1)=A(2,2)/det
        B(1,2)=-A(1,2)/det
        B(2,1)=-A(2,1)/det
        B(2,2)=A(1,1)/det
      ELSE
        CALL FlagWarning("Zero Determinant for matrix A.",err,error,*999)
        B=0.0_DP
      ENDIF
    CASE(3)
      det=A(1,1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+A(1,3)*A(2,1)*A(3,2)-A(1,1)*A(3,2)*A(2,3)-A(2,1)*A(1,2)*A(3,3)- &
        & A(3,1)*A(2,2)*A(1,3)
      IF(ABS(det)>ZERO_TOLERANCE_DP) THEN
        B(1,1)=(A(2,2)*A(3,3)-A(3,2)*A(2,3))/det
        B(2,1)=(A(2,3)*A(3,1)-A(3,3)*A(2,1))/det
        B(3,1)=(A(2,1)*A(3,2)-A(3,1)*A(2,2))/det
        B(1,2)=(A(3,2)*A(1,3)-A(1,2)*A(3,3))/det
        B(2,2)=(A(3,3)*A(1,1)-A(1,3)*A(3,1))/det
        B(3,2)=(A(3,1)*A(1,2)-A(1,1)*A(3,2))/det
        B(1,3)=(A(1,2)*A(2,3)-A(2,2)*A(1,3))/det
        B(2,3)=(A(1,3)*A(2,1)-A(2,3)*A(1,1))/det
        B(3,3)=(A(1,1)*A(2,2)-A(2,1)*A(1,2))/det
      ELSE
        CALL FlagWarning("Zero Determinant for matrix A.",err,error,*999)
        B=0.0_DP
      ENDIF
    CASE DEFAULT
      CALL FlagError("Matrix size is not implemented.",err,error,*999)
    END SELECT

    EXITS("InvertFullDP")
    RETURN
999 ERRORSEXITS("InvertFullDP",err,error)
    RETURN 1
    
  END SUBROUTINE InvertFullDP
  
  !
  !================================================================================================================================
  !

  !>Inverts the transpose of a full single precision matrix A to give matrix B and returns the determinant of A in det.
  SUBROUTINE InvertTransposeFullSP(A,B,det,err,error,*)
    
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The A matrix to invert the tranpose of
    REAL(SP), INTENT(OUT) :: B(:,:) !<On exit, the inverse of A^T
    REAL(SP), INTENT(OUT) :: det !<On exit, the determinant of A^T
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("InvertTransposeFullSP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix A is not square.",err,error,*999)
    IF(SIZE(B,1)/=SIZE(A,2).OR.SIZE(B,2)/=SIZE(A,1)) CALL FlagError("Matrix B is not the same size as matrix A.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1)) 
    CASE(1)
      det=A(1,1)
      IF(ABS(det)>ZERO_TOLERANCE_SP) THEN
        B(1,1)=1.0_SP/A(1,1)
      ELSE
        CALL FlagWarning("Matrix A is zero and cannot be inverted.",err,error,*999)
        B(1,1)=0.0_SP
      ENDIF
    CASE(2)
      det=A(1,1)*A(2,2)-A(2,1)*A(1,2)
      IF(ABS(det)>ZERO_TOLERANCE_SP) THEN
        B(1,1)=A(2,2)/det
        B(1,2)=-A(2,1)/det
        B(2,1)=-A(1,2)/det
        B(2,2)=A(1,1)/det
      ELSE
        CALL FlagWarning("Zero Determinant for matrix A.",err,error,*999)
        B=0.0_SP
      ENDIF
    CASE(3)
      det=A(1,1)*A(2,2)*A(3,3)+A(2,1)*A(3,2)*A(1,3)+A(3,1)*A(1,2)*A(2,3)-A(1,1)*A(2,3)*A(3,2)-A(1,2)*A(2,1)*A(3,3)- &
        & A(1,3)*A(2,2)*A(3,1)
      IF(ABS(det)>ZERO_TOLERANCE_SP) THEN
        B(1,1)=(A(2,2)*A(3,3)-A(2,3)*A(3,2))/det
        B(2,1)=(A(3,2)*A(1,3)-A(3,3)*A(1,2))/det
        B(3,1)=(A(1,2)*A(2,3)-A(1,3)*A(2,2))/det
        B(1,2)=(A(2,3)*A(3,1)-A(2,1)*A(3,3))/det
        B(2,2)=(A(3,3)*A(1,1)-A(3,1)*A(1,3))/det
        B(3,2)=(A(1,3)*A(2,1)-A(1,1)*A(2,3))/det
        B(1,3)=(A(2,1)*A(3,2)-A(2,2)*A(3,1))/det
        B(2,3)=(A(3,1)*A(1,2)-A(3,2)*A(1,1))/det
        B(3,3)=(A(1,1)*A(2,2)-A(1,2)*A(2,1))/det
      ELSE
        CALL FlagWarning("Zero Determinant for matrix A.",err,error,*999)
        B=0.0_SP
      ENDIF
    CASE DEFAULT
      CALL FlagError("Matrix size is not implemented.",err,error,*999)
    END SELECT

    EXITS("InvertTransposeFullSP")
    RETURN
999 ERRORSEXITS("InvertTransposeFullSP",err,error)
    RETURN 1
    
  END SUBROUTINE InvertTransposeFullSP

  !
  !================================================================================================================================
  !

  !>Inverts the transpose of a full double precision matrix A to give matrix B and returns the determinant of A in det.
  SUBROUTINE InvertTransposeFullDP(A,B,det,err,error,*)
    
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The A matrix to invert the tranpose of
    REAL(DP), INTENT(OUT) :: B(:,:) !<On exit, the inverse of A^T
    REAL(DP), INTENT(OUT) :: det !<On exit, the determinant of A^T
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("InvertTransposeFullDP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix A is not square.",err,error,*999)
    IF(SIZE(B,1)/=SIZE(A,2).OR.SIZE(B,2)/=SIZE(A,1)) CALL FlagError("Matrix B is not the same size as matrix A.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1)) 
    CASE(1)
      det=A(1,1)
      IF(ABS(det)>ZERO_TOLERANCE_DP) THEN
        B(1,1)=1.0_DP/A(1,1)
      ELSE
        CALL FlagWarning("Matrix A is zero and cannot be inverted.",err,error,*999)
        B(1,1)=0.0_DP
      ENDIF
    CASE(2)
      det=A(1,1)*A(2,2)-A(2,1)*A(1,2)
      IF(ABS(det)>ZERO_TOLERANCE_DP) THEN
        B(1,1)=A(2,2)/det
        B(1,2)=-A(2,1)/det
        B(2,1)=-A(1,2)/det
        B(2,2)=A(1,1)/det
      ELSE
        CALL FlagWarning("Zero Determinant for matrix A.",err,error,*999)
        B=0.0_DP
      ENDIF
    CASE(3)
      det=A(1,1)*A(2,2)*A(3,3)+A(2,1)*A(3,2)*A(1,3)+A(3,1)*A(1,2)*A(2,3)-A(1,1)*A(2,3)*A(3,2)-A(1,2)*A(2,1)*A(3,3)- &
        & A(1,3)*A(2,2)*A(3,1)
      IF(ABS(det)>ZERO_TOLERANCE_DP) THEN
        B(1,1)=(A(2,2)*A(3,3)-A(2,3)*A(3,2))/det
        B(2,1)=(A(3,2)*A(1,3)-A(3,3)*A(1,2))/det
        B(3,1)=(A(1,2)*A(2,3)-A(1,3)*A(2,2))/det
        B(1,2)=(A(2,3)*A(3,1)-A(2,1)*A(3,3))/det
        B(2,2)=(A(3,3)*A(1,1)-A(3,1)*A(1,3))/det
        B(3,2)=(A(1,3)*A(2,1)-A(1,1)*A(2,3))/det
        B(1,3)=(A(2,1)*A(3,2)-A(2,2)*A(3,1))/det
        B(2,3)=(A(3,1)*A(1,2)-A(3,2)*A(1,1))/det
        B(3,3)=(A(1,1)*A(2,2)-A(1,2)*A(2,1))/det
      ELSE
        CALL FlagWarning("Zero Determinant for matrix A.",err,error,*999)
        B=0.0_DP
      ENDIF
    CASE DEFAULT
      CALL FlagError("Matrix size is not implemented.",err,error,*999)
    END SELECT

    EXITS("InvertTransposeFullDP")
    RETURN
999 ERRORSEXITS("InvertTransposeFullDP",err,error)
    RETURN 1
    
  END SUBROUTINE InvertTransposeFullDP

  !
  !================================================================================================================================
  !

  !>Calculates the modified Bessel function of the second kind of order 1 using the approximation of Abromowitz and Stegun,
  !>for a double precision argument.
  PURE FUNCTION K0DP(x)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: x !<The value to evaluate the function at
    !Function variable
    REAL(DP) :: K0DP
    !Local variables
    REAL(DP) :: a1,a2,a3,a4,a5,a6,a7,t

    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    !Calculate K0(x)
    IF(x<=2.0_DP) THEN
      t=x*x/4.0_DP
      a1=-0.57721566_DP
      a2=0.42278420_DP
      a3=0.23069756_DP
      a4=0.03488590_DP
      a5=0.00262698_DP
      a6=0.00010750_DP
      a7=0.00000740_DP
      K0DP=-LOG(x/2.0_DP)*I0(x)+(a1+(a2+(a3+(a4+(a5+(a6+a7*t)*t)*t)*t)*t)*t)
    ELSE
      IF(x>174.0_DP) THEN
        K0DP=0.0_DP
      ELSE
        t=2.0_DP/x
        a1=1.25331414_DP
        a2=-0.07832358_DP
        a3=0.02189568_DP
        a4=-0.01062446_DP
        a5=0.00587872_DP
        a6=-0.00251540_DP
        a7=0.00053208_DP
        K0DP=(a1+(a2+(a3+(a4+(a5+(a6+a7*t)*t)*t)*t)*t)*t)*EXP(-x)/SQRT(x)
      ENDIF
    ENDIF
    
    RETURN
    
  END FUNCTION K0DP
  
  !
  !================================================================================================================================
  !
  
  !>Calculates the modified Bessel function of the second kind of order 0 using the approximation of Abromowitz and Stegun,
  !>for a single precision argument.
  PURE FUNCTION K0SP(x)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: x !<The value to evaluate the function at
    !Function variable
    REAL(SP) :: K0SP
    !Local variables
    REAL(SP) :: a1,a2,a3,a4,a5,a6,a7,t

    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    !Calculate K0(x)
    IF(x<=2.0_SP) THEN
      t=x*x/4.0_SP
      a1=-0.57721566_SP
      a2=0.42278420_SP
      a3=0.23069756_SP
      a4=0.03488590_SP
      a5=0.00262698_SP
      a6=0.00010750_SP
      a7=0.00000740_SP
      K0SP=-LOG(x/2.0_SP)*I0(x)+(a1+(a2+(a3+(a4+(a5+(a6+a7*t)*t)*t)*t)*t)*t)
    ELSE
      IF(x>174.0_SP) THEN
        K0SP=0.0_SP
      ELSE
        t=2.0_SP/x
        a1=1.25331414_SP
        a2=-0.07832358_SP
        a3=0.02189568_SP
        a4=-0.01062446_SP
        a5=0.00587872_SP
        a6=-0.00251540_SP
        a7=0.00053208_SP
        K0SP=(a1+(a2+(a3+(a4+(a5+(a6+a7*t)*t)*t)*t)*t)*t)*EXP(-x)/SQRT(x)
      ENDIF
    ENDIF
    
    RETURN
    
  END FUNCTION K0SP
  
  !
  !================================================================================================================================
  !
  
  !>Calculates the modified Bessel function of the second kind of order 1 using the approximation of Abromowitz and Stegun,
  !>for a double precision argument.
  PURE FUNCTION K1DP(x)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: x !<The value to evaluate the function at
    !Function variable
    REAL(DP) :: K1DP
    !Local variables
    REAL(DP) :: a1,a2,a3,a4,a5,a6,a7,a8,t

    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    !Calculate K1(x)
    IF(x<=2.0_DP) THEN
      t=(x/2.0_DP)*(x/2.0_DP)
      a1=LOG(x/2.0_DP)*I1(x)
      a2=1.0_DP
      a3=0.15443144_DP
      a4=-0.67278579_DP
      a5=-0.18156897_DP
      a6=-0.01919402_DP
      a7=-0.00110404_DP
      a8=-0.00004686_DP
      K1DP=a1+a2/x+(a3+(a4+(a5+(a6+(a7+a8*t)*t)*t)*t)*t)*x/4.0_DP
    ELSE
      IF (x>174.0_DP) THEN
        K1DP=0.0_DP
      ELSE
        t=2.0_DP/x
        a1=1.25331414_DP
        a2=0.23498619_DP
        a3=-0.03655620_DP
        a4=0.01504268_DP
        a5=-0.00780355_DP
        a6=0.00325614_DP
        a7=-0.00068245_DP
        K1DP=(a1+(a2+(a3+(a4+(a5+(a6+a7*t)*t)*t)*t)*t)*t)*EXP(-x)/SQRT(x)
      ENDIF
    ENDIF
    
    RETURN
    
  END FUNCTION K1DP
  
  !
  !================================================================================================================================
  !
  
  !>Calculates the modified Bessel function of the second kind of order 1 using the approximation of Abromowitz and Stegun,
  !>for a single precision argument.
  PURE FUNCTION K1SP(x)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: x !<The value to evaluate the function at
    !Function variable
    REAL(SP) :: K1SP
    !Local variables
    REAL(SP) :: a1,a2,a3,a4,a5,a6,a7,a8,t

    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    !Calculate K1(x)
    IF(x<=2.0_SP) THEN
      t=(x/2.0_SP)*(x/2.0_SP)
      a1=LOG(x/2.0_SP)*I1(x)
      a2=1.0_SP
      a3=0.15443144_SP
      a4=-0.67278579_SP
      a5=-0.18156897_SP
      a6=-0.01919402_SP
      a7=-0.00110404_SP
      a8=-0.00004686_SP
      K1SP=a1+a2/x+(a3+(a4+(a5+(a6+(a7+a8*t)*t)*t)*t)*t)*x/4.0_SP
    ELSE
      IF (x>174.0_SP) THEN
        K1SP=0.0_SP
      ELSE
        t=2.0_SP/x
        a1=1.25331414_SP
        a2=0.23498619_SP
        a3=-0.03655620_SP
        a4=0.01504268_SP
        a5=-0.00780355_SP
        a6=0.00325614_SP
        a7=-0.00068245_SP
        K1SP=(a1+(a2+(a3+(a4+(a5+(a6+a7*t)*t)*t)*t)*t)*t)*EXP(-x)/SQRT(x)
      ENDIF
    ENDIF
    
    RETURN
    
  END FUNCTION K1SP
  
  !
  !================================================================================================================================
  !
  
  !>Calculates the elliptic integral of the first kind - K(m), for a double precision argument.
  PURE FUNCTION KdpDP(x)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: x !<The value to evaluate the function at
    !Function variable
    REAL(DP) :: KdpDP
    !Local variables
    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    REAL(DP), PARAMETER :: a0=1.38629436112_DP
    REAL(DP), PARAMETER :: a1=0.09666344259_DP
    REAL(DP), PARAMETER :: a2=0.03590092383_DP
    REAL(DP), PARAMETER :: a3=0.03742563713_DP
    REAL(DP), PARAMETER :: a4=0.01451196212_DP
    REAL(DP), PARAMETER :: b0=0.5_DP
    REAL(DP), PARAMETER :: b1=0.12498593597_DP
    REAL(DP), PARAMETER :: b2=0.06880248576_DP
    REAL(DP), PARAMETER :: b3=0.03328355346_DP
    REAL(DP), PARAMETER :: b4=0.00441787012_DP
    REAL(DP) :: term1,term2,x1

    x1=1.0_DP-x
    term1=a0+(a1+(a2+(a3+a4*x)*x)*x)*x
    term2=b0+(b1+(b2+(b3+b4*x)*x)*x)*x
    KdpDP=term1+term2*LOG(1.0_DP/x)    
    
    RETURN
    
  END FUNCTION KdpDP
  
  !
  !================================================================================================================================
  !
  
  !>Calculates the elliptic integral of the first kind - K(m), for a single precision argument.
  PURE FUNCTION KdpSP(x)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: x !<The value to evaluate the function at
    !Function variable
    REAL(SP) :: KdpSP
    !Local variables
    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    REAL(SP), PARAMETER :: a0=1.38629436112_SP
    REAL(SP), PARAMETER :: a1=0.09666344259_SP
    REAL(SP), PARAMETER :: a2=0.03590092383_SP
    REAL(SP), PARAMETER :: a3=0.03742563713_SP
    REAL(SP), PARAMETER :: a4=0.01451196212_SP
    REAL(SP), PARAMETER :: b0=0.5_SP
    REAL(SP), PARAMETER :: b1=0.12498593597_SP
    REAL(SP), PARAMETER :: b2=0.06880248576_SP
    REAL(SP), PARAMETER :: b3=0.03328355346_SP
    REAL(SP), PARAMETER :: b4=0.00441787012_SP
    REAL(SP) :: term1,term2,x1

    x1=1.0_SP-x
    term1=a0+(a1+(a2+(a3+a4*x)*x)*x)*x
    term2=B0+(b1+(b2+(b3+b4*x)*x)*x)*x
    KdpSP=term1+term2*LOG(1.0_SP/x)    
    
    RETURN
    
  END FUNCTION KdpSP
  
  !
  !================================================================================================================================
  !
  
  !>Returns the L1-norm of the single precision matrix A.
  SUBROUTINE L1NormMatrixSP(A,L1Norm,err,error,*)

    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The matrix to calculate the L1 norm of
    REAL(SP), INTENT(OUT) :: L1Norm !<On exit, the L1 norm of A
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j
    REAL(SP) :: columnSum

    ENTERS("L1NormMatrixSP",err,error,*999)

    L1Norm=0.0_SP
    DO j=1,SIZE(A,2)
      columnSum=0.0_SP
      DO i=1,SIZE(A,1)
        columnSum=columnSum+ABS(A(i,j))
      ENDDO !i
      IF(columnSum>L1Norm) L1Norm=columnSum
    ENDDO !j

    EXITS("L1NormMatrixSP")
    RETURN
999 ERRORSEXITS("L1NormMatrixSP",err,error)
    RETURN 1
    
  END SUBROUTINE L1NormMatrixSP

  !
  !================================================================================================================================
  !
  
  !>Returns the L1-norm of the double precision matrix A.
  SUBROUTINE L1NormMatrixDP(A,L1Norm,err,error,*)

    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The matrix to calculate the L1 norm of
    REAL(DP), INTENT(OUT) :: L1Norm !<On exit, the L1 norm of A
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j
    REAL(DP) :: columnSum

    ENTERS("L1NormMatrixDP",err,error,*999)

    L1Norm=0.0_DP
    DO j=1,SIZE(A,2)
      columnSum=0.0_DP
      DO i=1,SIZE(A,1)
        columnSum=columnSum+ABS(A(i,j))
      ENDDO !i
      IF(columnSum>L1Norm) L1Norm=columnSum
    ENDDO !j

    EXITS("L1NormMatrixDP")
    RETURN
999 ERRORSEXITS("L1NormMatrixDP",err,error)
    RETURN 1
    
  END SUBROUTINE L1NormMatrixDP

  !
  !================================================================================================================================
  !
  
  !>Returns the L1-norm of the single precision vector a.
  SUBROUTINE L1NormVectorSP(a,L1Norm,err,error,*)

    !Argument variables
    REAL(SP), INTENT(IN) :: a(:) !<The vector to calculate the L1 norm of
    REAL(SP), INTENT(OUT) :: L1Norm !<On exit, the L1 norm of a
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i

    ENTERS("L1NormVectorSP",err,error,*999)

    L1Norm=0.0_SP
    DO i=1,SIZE(a,1)
      L1Norm=L1Norm+ABS(a(i))
    ENDDO !i

    EXITS("L1NormVectorSP")
    RETURN
999 ERRORSEXITS("L1NormVectorSP",err,error)
    RETURN 1
    
  END SUBROUTINE L1NormVectorSP

  !
  !================================================================================================================================
  !

  !>Returns the L1-norm of the double precision vector a.
  SUBROUTINE L1NormVectorDP(a,L1Norm,err,error,*)

    !Argument variables
    REAL(DP), INTENT(IN) :: a(:) !<The vector to calculate the L1 norm of
    REAL(DP), INTENT(OUT) :: L1Norm !<On exit, the L1 norm of a
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i

    ENTERS("L1NormVectorDP",err,error,*999)

    L1Norm=0.0_DP
    DO i=1,SIZE(a,1)
      L1Norm=L1Norm+ABS(a(i))
    ENDDO !i

    EXITS("L1NormVectorDP")
    RETURN
999 ERRORSEXITS("L1NormVectorDP",err,error)
    RETURN 1
    
  END SUBROUTINE L1NormVectorDP
  
  !
  !================================================================================================================================
  !
  
  !>Returns the L2-norm of the single precision matrix A.
  SUBROUTINE L2NormMatrixSP(a,L2Norm,err,error,*)

    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The matrix to calculate the L2 norm of
    REAL(SP), INTENT(OUT) :: L2Norm !<On exit, the L2 norm of A
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i
    REAL(SP) :: ATA(SIZE(A,2),SIZE(A,2)),eValues(SIZE(A,2))

    ENTERS("L2NormMatrixSP",err,error,*999)

    CALL MatrixTransposeProduct(A,A,ATA,err,error,*999)
    CALL Eigenvalue(ATA,eValues,err,error,*999)
    L2Norm=0.0_SP
    DO i=1,SIZE(eValues,1)
      IF(eValues(i)>L2Norm) L2Norm=eValues(i)
    ENDDO !i
    
    EXITS("L2NormMatrixSP")
    RETURN
999 ERRORSEXITS("L2NormMatrixSP",err,error)
    RETURN 1
    
  END SUBROUTINE L2NormMatrixSP

  !
  !================================================================================================================================
  !
  
  !>Returns the L2-norm of the double precision matrix A.
  SUBROUTINE L2NormMatrixDP(A,L2Norm,err,error,*)

    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The matrix to calculate the L2 norm of
    REAL(DP), INTENT(OUT) :: L2Norm !<On exit, the L2 norm of A
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i
    REAL(DP) :: ATA(SIZE(A,2),SIZE(A,2)),eValues(SIZE(A,2))

    ENTERS("L2NormMatrixDP",err,error,*999)

    CALL MatrixTransposeProduct(A,A,ATA,err,error,*999)
    CALL Eigenvalue(ATA,eValues,err,error,*999)
    L2Norm=0.0_DP
    DO i=1,SIZE(eValues,1)
      IF(eValues(i)>L2Norm) L2Norm=eValues(i)
    ENDDO !i
    
    EXITS("L2NormMatrixDP")
    RETURN
999 ERRORSEXITS("L2NormMatrixDP",err,error)
    RETURN 1
    
  END SUBROUTINE L2NormMatrixDP

  !
  !================================================================================================================================
  !
  
  !>Returns the L2-norm of the single precision vector a.
  SUBROUTINE L2NormVectorSP(a,L2Norm,err,error,*)

    !Argument variables
    REAL(SP), INTENT(IN) :: a(:) !<The vector to calculate the L2 norm of
    REAL(SP), INTENT(OUT) :: L2Norm !<On exit, the L2 norm of a
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    REAL(SP) :: asum
    
    ENTERS("L2NormVectorSP",err,error,*999)

    asum=SUM(a*a,1)
    L2Norm=SQRT(asum)

    EXITS("L2NormVectorSP")
    RETURN
999 ERRORSEXITS("L2NormVectorSP",err,error)
    RETURN 1
    
  END SUBROUTINE L2NormVectorSP

  !
  !================================================================================================================================
  !

  !>Returns the L2-norm of the double precision vector a.
  SUBROUTINE L2NormVectorDP(a,L2Norm,err,error,*)

    !Argument variables
    REAL(DP), INTENT(IN) :: a(:) !<The vector to calculate the L2 norm of
    REAL(DP), INTENT(OUT) :: L2Norm !<On exit, the L2 norm of a
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    REAL(DP) :: asum

    ENTERS("L2NormVectorDP",err,error,*999)

    asum=SUM(a*a,1)
    L2Norm=SQRT(asum)

    EXITS("L2NormVectorDP")
    RETURN
999 ERRORSEXITS("L2NormVectorDP",err,error)
    RETURN 1
    
  END SUBROUTINE L2NormVectorDP

  !
  !================================================================================================================================
  !
  
  !>Returns the Linf-norm of the single precision matrix A.
  SUBROUTINE LInfNormMatrixSP(A,LInfNorm,err,error,*)

    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The matrix to calculate the Linf norm of
    REAL(SP), INTENT(OUT) :: LInfNorm !<On exit, the Linf norm of A
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j
    REAL(SP) :: rowSum

    ENTERS("LInfNormMatrixSP",err,error,*999)

    LInfNorm=0.0_SP
    DO i=1,SIZE(A,1)
      rowSum=0.0_SP
      DO j=1,SIZE(A,2)
        rowSum=rowSum+ABS(A(i,j))
      ENDDO !i
      IF(rowSum>LInfNorm) LInfNorm=rowSum
    ENDDO !j

    EXITS("LInfNormMatrixSP")
    RETURN
999 ERRORSEXITS("LInfNormMatrixSP",err,error)
    RETURN 1
    
  END SUBROUTINE LInfNormMatrixSP

  !
  !================================================================================================================================
  !
  
  !>Returns the Linf-norm of the double precision matrix A.
  SUBROUTINE LInfNormMatrixDP(A,LInfNorm,err,error,*)

    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The matrix to calculate the Linf norm of
    REAL(DP), INTENT(OUT) :: LInfNorm !<On exit, the Linf norm of A
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j
    REAL(DP) :: rowSum

    ENTERS("LInfNormMatrixDP",err,error,*999)

    LInfNorm=0.0_DP
    DO i=1,SIZE(A,1)
      rowSum=0.0_DP
      DO j=1,SIZE(A,2)
        rowSum=rowSum+ABS(A(i,j))
      ENDDO !i
      IF(rowSum>LInfNorm) LInfNorm=rowSum
    ENDDO !j

    EXITS("LInfNormMatrixDP")
    RETURN
999 ERRORSEXITS("LInfNormMatrixDP",err,error)
    RETURN 1
    
  END SUBROUTINE LInfNormMatrixDP

  !
  !================================================================================================================================
  !
  
  !>Returns the Linf-norm of the single precision vector a.
  SUBROUTINE LInfNormVectorSP(a,LInfNorm,err,error,*)

    !Argument variables
    REAL(SP), INTENT(IN) :: a(:) !<The vector to calculate the Linf norm of
    REAL(SP), INTENT(OUT) :: LInfNorm !<On exit, the Linf norm of a
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i

    ENTERS("LInfNormVectorSP",err,error,*999)

    LInfNorm=0.0_SP
    DO i=1,SIZE(a,1)
      IF(ABS(a(i))>LInfNorm) LInfNorm=ABS(a(i))
    ENDDO !i
    
    EXITS("LInfNormVectorSP")
    RETURN
999 ERRORSEXITS("LInfNormVectorSP",err,error)
    RETURN 1
    
  END SUBROUTINE LInfNormVectorSP
 
  !
  !================================================================================================================================
  !

  !>Returns the Linf-norm of the double precision vector a.
  SUBROUTINE LInfNormVectorDP(a,LInfNorm,err,error,*)

    !Argument variables
    REAL(DP), INTENT(IN) :: a(:) !<The vector to calculate the Linf norm of
    REAL(DP), INTENT(OUT) :: LInfNorm !<On exit, the Linf norm of a
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i

    ENTERS("LInfNormVectorDP",err,error,*999)

    LInfNorm=0.0_DP
    DO i=1,SIZE(a,1)
      IF(ABS(a(i))>LInfNorm) LInfNorm=ABS(a(i))
    ENDDO !i
    
    EXITS("LInfNormVectorDP")
    RETURN
999 ERRORSEXITS("LInfNormVectorDP",err,error)
    RETURN 1
    
  END SUBROUTINE LInfNormVectorDP

  !
  !================================================================================================================================
  !
  
  !>Calculates the Macaulay bracket <x> of a single precision number x
  PURE FUNCTION MacaulayBracketSP(x)

    !Argument variables
    REAL(SP), INTENT(IN) :: x !<argument to calculate <x> with
    !Function variable
    REAL(SP) :: MacaulayBracketSP
    
    MacaulayBracketSP=MAX(x,0.0_SP)

    RETURN
    
  END FUNCTION MacaulayBracketSP

  !
  !================================================================================================================================
  !
  
  !>Calculates the Macaulay bracket <x> of a double precision number x
  PURE FUNCTION MacaulayBracketDP(x)

    !Argument variables
    REAL(DP), INTENT(IN) :: x !<argument to calculate <x> with
    !Function variable
    REAL(DP) :: MacaulayBracketDP
    
    MacaulayBracketDP=MAX(x,0.0_DP)

    RETURN
    
  END FUNCTION MacaulayBracketDP

  !
  !================================================================================================================================
  !
  
  !>Normalises a real single precision matrix A.
  SUBROUTINE NormaliseMatrixSP(A,normA,err,error,*)

    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The matrix to normalise
    REAL(SP), INTENT(OUT) :: normA(:,:) !<On exit, the normalised matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    REAL(SP) :: detA
    
    ENTERS("NormaliseMatrixSP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(normA,1).OR.SIZE(A,2)/=SIZE(normA,2)) CALL FlagError("Invalid matrix sizes.",err,error,*999)
#endif    
    
    CALL Determinant(A,detA,err,error,*999)
    IF(ABS(detA)<ZERO_TOLERANCE_SP) THEN
      normA=A
      CALL FlagError("Can not normalise A as it has a zero determinant.",err,error,*999)
    ELSE
      normA=A/detA
    ENDIF

    EXITS("NormaliseMatrixSP")
    RETURN
999 ERRORSEXITS("NormaliseMatrixSP",err,error)
    RETURN 1
    
  END SUBROUTINE NormaliseMatrixSP

 !
  !================================================================================================================================
  !
  
  !>Normalises a real double precision matrix A.
  SUBROUTINE NormaliseMatrixDP(A,normA,err,error,*)

    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The matrix to normalise
    REAL(DP), INTENT(OUT) :: normA(:,:) !<On exit, the normalised matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    REAL(DP) :: detA
    
    ENTERS("NormaliseMatrixDP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(normA,1).OR.SIZE(A,2)/=SIZE(normA,2)) CALL FlagError("Invalid matrix sizes.",err,error,*999)
#endif    
    
    CALL Determinant(A,detA,err,error,*999)
    IF(ABS(detA)<ZERO_TOLERANCE) THEN
      normA=A
      CALL FlagError("Can not normalise A as it has a zero determinant.",err,error,*999)
    ELSE
      normA=A/detA
    ENDIF

    EXITS("NormaliseMatrixDP")
    RETURN
999 ERRORSEXITS("NormaliseMatrixDP",err,error)
    RETURN 1
    
  END SUBROUTINE NormaliseMatrixDP

  !
  !================================================================================================================================
  !
  
  !>Normalises a real single precision vector a.
  SUBROUTINE NormaliseVectorSP(a,norma,err,error,*)

    !Argument variables
    REAL(SP), INTENT(IN) :: a(:) !<The vector to normalise
    REAL(SP), INTENT(OUT) :: norma(:) !<On exit, the normalised vector
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    REAL(SP) :: length
    
    ENTERS("NormaliseVectorSP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(a,1)/=SIZE(normA,1)) CALL FlagError("Invalid vector sizes.",err,error,*999)
#endif    
    
    CALL L2Norm(a,length,err,error,*999)
    IF(ABS(length)<ZERO_TOLERANCE_SP) THEN
      norma=0.0_SP
    ELSE
      norma=a/length
    ENDIF

    EXITS("NormaliseVectorSP")
    RETURN
999 ERRORSEXITS("NormaliseVectorSP",err,error)
    RETURN 1
    
  END SUBROUTINE NormaliseVectorSP

  !
  !================================================================================================================================
  !

  !>Normalises a real double precision vector a.
  SUBROUTINE NormaliseVectorDP(a,norma,err,error,*)

    !Argument variables
    REAL(DP), INTENT(IN) :: a(:) !<The vector to normalise
    REAL(DP), INTENT(OUT) :: norma(:) !<On exit, the vector to normalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    REAL(DP) :: length
    
    ENTERS("NormaliseVectorDP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(a,1)/=SIZE(normA,1)) CALL FlagError("Invalid vector sizes.",err,error,*999)
#endif    
    
    CALL L2Norm(a,length,err,error,*999)
    IF(ABS(length)<ZERO_TOLERANCE_DP) THEN
      norma=0.0_DP
    ELSE
      norma=a/length
    ENDIF

    EXITS("NormaliseVectorDP")
    RETURN
999 ERRORSEXITS("NormaliseVectorDP",err,error)
    RETURN 1
    
  END SUBROUTINE NormaliseVectorDP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the normalised vector cross-prouct of the single precision vectors a x b in c.
  SUBROUTINE NormaliseCrossProductSP(a,b,c,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: a(:) !<The first vector in the cross product
    REAL(SP), INTENT(IN) :: b(:) !<The second vector in the cross product
    REAL(SP), INTENT(OUT) :: c(:) !<On exit, the normalised cross product of the first and second vectors
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    REAL(SP) :: temp(SIZE(c,1))
    
    ENTERS("NormaliseCrossProductSP",err,error,*999)

    CALL CrossProduct(a,b,temp,err,error,*999)
    CALL Normalise(temp,c,err,error,*999)

    EXITS("NormaliseCrossProductSP")
    RETURN
999 ERRORSEXITS("NormaliseCrossProductSP",err,error)
    RETURN 1
    
  END SUBROUTINE NormaliseCrossProductSP
  
  !
  !================================================================================================================================
  !

  !>Calculates and returns the normalised vector cross-prouct of the double precision vectors a x b in c.
  SUBROUTINE NormaliseCrossProductDP(a,b,c,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: a(:) !<The first vector in the cross product
    REAL(DP), INTENT(IN) :: b(:) !<The second vector in the cross product
    REAL(DP), INTENT(OUT) :: c(:) !<On exit, the normalised cross product of the first and second vectors
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    REAL(DP) :: temp(SIZE(c,1))
    
    ENTERS("NormaliseCrossProductDP",err,error,*999)

    CALL CrossProduct(a,b,temp,err,error,*999)
    CALL Normalise(temp,c,err,error,*999)

    EXITS("NormaliseCrossProductDP")
    RETURN
999 ERRORSEXITS("NormaliseCrossProductDP",err,error)
    RETURN 1
    
  END SUBROUTINE NormaliseCrossProductDP
  
  !
  !================================================================================================================================
  !
  
  !>Calculates single precision hyperbolic secant function
  PURE FUNCTION SechSP(a)

    !Argument variables
    REAL(SP), INTENT(IN) :: a !<argument to perform sech() on
    !Function variable
    REAL(SP) :: SechSP
    
    SechSP=1.0/COSH(a)

    RETURN
    
  END FUNCTION SechSP

  !
  !================================================================================================================================
  !

  !>Calculates double precision hyperbolic secant function
  PURE FUNCTION SechDP(a)

    !Argument variables
    REAL(DP), INTENT(IN) :: a !<argument to perform sech() on
    !Function variable
    REAL(DP) :: SechDP

    SechDP=1.0_DP/COSH(a)

    RETURN
    
  END FUNCTION SechDP

  !
  !================================================================================================================================
  !

  !>Finds the solution to a single precision eigenproblem 
  SUBROUTINE SolveEigenproblemSP(A,eigenvalues,eigenvectors,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The A matrix to solve the eigenproblem for
    REAL(SP), INTENT(OUT) :: eigenvalues(:) !<On exit, the eigenvalues of A
    REAL(SP), INTENT(OUT) :: eigenvectors(:,:) !<On exit, the eigenvectors of A
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: matrixSize,lWork
    INTEGER(INTG), PARAMETER :: lWMax=1000
    REAL(SP) :: B(SIZE(A,1),SIZE(A,2)),det,factor,trace,work(lWMax)
    
    ENTERS("SolveEigenproblemSP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix is not square.",err,error,*999)
    IF(SIZE(A,1)>SIZE(eigenvalues,1)) &
      & CALL FlagError("Size of eigenvalues is not the same as the number of rows in A.",err,error,*999)    
    IF(SIZE(A,1)>SIZE(eigenvectors,1)) CALL FlagError("Eigenvectors is too small.",err,error,*999)
    IF(SIZE(A,2)>SIZE(eigenvectors,2)) CALL FlagError("Eigenvectors is too small.",err,error,*999)
#endif    

    matrixSize=SIZE(A,1)
    SELECT CASE(matrixSize)
    CASE(1)
      eigenvalues(1)=A(1,1)
      eigenvectors(1,1)=0.0_SP
    CASE(2)
!!TODO: This needs to be checked      
      trace=A(1,1)+A(2,2)
      det=A(1,1)*A(2,2)-A(2,1)*A(1,2)
      factor=trace*trace/4.0_SP-det
      IF(factor<=ZERO_TOLERANCE_SP) CALL FlagError("Complex roots encountered.",err,error,*999)
      eigenvalues(1)=trace/2.0_SP+SQRT(factor)
      eigenvalues(2)=trace/2.0_SP-SQRT(factor)
      IF(ABS(A(2,1)>ZERO_TOLERANCE_SP)) THEN
        eigenvectors(1,1)=eigenvalues(1)-A(2,2)
        eigenvectors(2,1)=A(2,1)
        eigenvectors(1,2)=eigenvalues(2)-A(2,2)
        eigenvectors(2,2)=A(2,1)
      ELSE
        IF(ABS(A(1,2))>ZERO_TOLERANCE_SP) THEN
          eigenvectors(1,1)=A(1,2)
          eigenvectors(2,1)=eigenvalues(1)-A(1,1)
          eigenvectors(1,2)=A(1,2)
          eigenvectors(2,2)=eigenvalues(2)-A(1,1)
        ELSE
          eigenvectors(1,1)=1.0_SP
          eigenvectors(2,1)=0.0_SP
          eigenvectors(1,2)=0.0_SP
          eigenvectors(2,2)=1.0_SP
         ENDIF
      ENDIF      
    CASE DEFAULT
      B(1:matrixSize,1:matrixSize)=A(1:matrixSize,1:matrixSize)
      CALL SSYEV('V','U',matrixSize,B,matrixSize,eigenvalues,work,-1,err)
      IF(err/=0) CALL FlagError("Error in eigenvalue computation SSYEV.",err,error,*999)
      lWork=MIN(lWMax,INT(work(1)))
      CALL SSYEV('V','U',matrixSize,B,matrixSize,eigenvalues,work,lWork,err)
      IF(err/=0) CALL FlagError("Error in eigenvectror computation SSYEV.",err,error,*999)
      eigenvectors(1:matrixSize,1)=B(1:matrixSize,1)
      eigenvectors(1:matrixSize,2)=B(1:matrixSize,2)
      eigenvectors(1:matrixSize,3)=B(1:matrixSize,3)
    END SELECT

    EXITS("SolveEigenproblemSP")
    RETURN
999 ERRORSEXITS("SolveEigenproblemSP",err,error)
    RETURN 1
    
  END SUBROUTINE SolveEigenproblemSP

  !
  !================================================================================================================================
  !

  !>Finds the solution to a double precision eigenproblem 
  SUBROUTINE SolveEigenproblemDP(A,eigenvalues,eigenvectors,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The A matrix to solve the eigenproblem for
    REAL(DP), INTENT(OUT) :: eigenvalues(:) !<On exit, the eigenvalues of A
    REAL(DP), INTENT(OUT) :: eigenvectors(:,:) !<On exit, the eigenvectors of A
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: matrixSize,lWork
    INTEGER(INTG), PARAMETER :: lWMax=1000
    REAL(DP) :: B(SIZE(A,1),SIZE(A,2)),det,factor,trace,work(lWMax)
    
    ENTERS("SolveEigenproblemDP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix is not square.",err,error,*999)
    IF(SIZE(A,1)>SIZE(eigenvalues,1)) &
      & CALL FlagError("Size of eigenvalues is not the same as the number of rows in A.",err,error,*999)    
    IF(SIZE(A,1)>SIZE(eigenvectors,1)) CALL FlagError("Eigenvectors is too small.",err,error,*999)
    IF(SIZE(A,2)>SIZE(eigenvectors,2)) CALL FlagError("Eigenvectors is too small.",err,error,*999)
#endif    

    matrixSize=SIZE(A,1)
    SELECT CASE(matrixSize)
    CASE(1)
      eigenvalues(1)=A(1,1)
      eigenvectors(1,1)=0.0_DP
    CASE(2)
!!TODO: This needs to be checked      
      trace=A(1,1)+A(2,2)
      det=A(1,1)*A(2,2)-A(2,1)*A(1,2)
      factor=trace*trace/4.0_DP-det
      IF(factor<=ZERO_TOLERANCE) CALL FlagError("Complex roots encountered.",err,error,*999)
      eigenvalues(1)=trace/2.0_DP+SQRT(factor)
      eigenvalues(2)=trace/2.0_DP-SQRT(factor)
      IF(ABS(A(2,1)>ZERO_TOLERANCE)) THEN
        eigenvectors(1,1)=eigenvalues(1)-A(2,2)
        eigenvectors(2,1)=A(2,1)
        eigenvectors(1,2)=eigenvalues(2)-A(2,2)
        eigenvectors(2,2)=A(2,1)
      ELSE
        IF(ABS(A(1,2))>ZERO_TOLERANCE) THEN
          eigenvectors(1,1)=A(1,2)
          eigenvectors(2,1)=eigenvalues(1)-A(1,1)
          eigenvectors(1,2)=A(1,2)
          eigenvectors(2,2)=eigenvalues(2)-A(1,1)
        ELSE
          eigenvectors(1,1)=1.0_DP
          eigenvectors(2,1)=0.0_DP
          eigenvectors(1,2)=0.0_DP
          eigenvectors(2,2)=1.0_DP
         ENDIF
      ENDIF      
    CASE DEFAULT
      B(1:matrixSize,1:matrixSize)=A(1:matrixSize,1:matrixSize)
      CALL DSYEV('V','U',matrixSize,B,matrixSize,eigenvalues,work,-1,err)
      IF(err/=0) CALL FlagError("Error in eigenvalue computation DSYEV.",err,error,*999)
      lWork=MIN(lWMax,INT(work(1)))
      CALL DSYEV('V','U',matrixSize,B,matrixSize,eigenvalues,work,lWork,err)
      IF(err/=0) CALL FlagError("Error in eigenvector computation DSYEV.",err,error,*999)
      eigenvectors(1:matrixSize,1)=B(1:matrixSize,1)
      eigenvectors(1:matrixSize,2)=B(1:matrixSize,2)
      eigenvectors(1:matrixSize,3)=B(1:matrixSize,3)
    END SELECT

    EXITS("SolveEigenproblemDP")
    RETURN
999 ERRORSEXITS("SolveEigenproblemDP",err,error)
    RETURN 1
    
  END SUBROUTINE SolveEigenproblemDP

  !
  !================================================================================================================================
  !

  !>Finds the solution to a small single precision linear system Ax=b.
  SUBROUTINE SolveSmallLinearSystemSP(A,x,b,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(INOUT) :: A(:,:) !<The A matrix
    REAL(SP), INTENT(OUT) :: x(:) !<On exit, the solution vector x
    REAL(SP), INTENT(INOUT) :: b(:) !<The RHS vector b
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: info,iPivot(SIZE(A,1))
    REAL(SP) :: AInv(SIZE(A,1),SIZE(A,2)),det
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolveSmallLinearSystemSP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix is not square.",err,error,*999)
    IF(SIZE(A,1)/=SIZE(b,1)) CALL FlagError("Size of b is not the same as the number of rows in A.",err,error,*999)    
    IF(SIZE(A,1)>SIZE(x,1)) CALL FlagError("x is too small.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1))
    CASE(1)
      IF(ABS(A(1,1))<ZERO_TOLERANCE_SP) CALL FlagError("The A matrix is zero.",err,error,*999)
      x(1)=b(1)/A(1,1)
    CASE(2:3)
      CALL Invert(A,AInv,det,err,error,*999)
      x=MATMUL(AInv,b)
    CASE DEFAULT
      !Use LAPACK LU
      CALL DGETRF(SIZE(A,1),SIZE(A,2),A,SIZE(A,1),iPivot,info)
      IF(info/=0) THEN
        IF(info<0) THEN
          localError="Parameter number "//TRIM(NumberToVString(-info,"*",err,error))//" is illegal."
        ELSE
          localError="The A matrix is singular. The factorised diagonal entry for row "// &
            & TRIM(NumberToVString(info,"*",err,error))//" is zero."          
        ENDIF
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Solve the LU problem
      CALL DGETRS('N',1,A,SIZE(A,1),iPivot,b,SIZE(b,1),info)
      IF(info/=0) THEN
        localError="Parameter number "//TRIM(NumberToVString(-info,"*",err,error))//" is illegal."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      x=b
    END SELECT

    EXITS("SolveSmallLinearSystemSP")
    RETURN
999 ERRORSEXITS("SolveSmallLinearSystemSP",err,error)
    RETURN 1
    
  END SUBROUTINE SolveSmallLinearSystemSP

  !
  !================================================================================================================================
  !

  !>Finds the solution to a small double precision linear system Ax=b.
  SUBROUTINE SolveSmallLinearSystemDP(A,x,b,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(INOUT) :: A(:,:) !<The A matrix
    REAL(DP), INTENT(OUT) :: x(:) !<On exit, the solution vector x
    REAL(DP), INTENT(INOUT) :: b(:) !<The RHS vector b
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: info,iPivot(SIZE(A,1))
    REAL(DP) :: AInv(SIZE(A,1),SIZE(A,2)),det
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolveSmallLinearSystemDP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix is not square.",err,error,*999)
    IF(SIZE(A,1)/=SIZE(b,1)) CALL FlagError("Size of b is not the same as the number of rows in A.",err,error,*999)    
    IF(SIZE(A,1)>SIZE(x,1)) CALL FlagError("x is too small.",err,error,*999)
#endif    

    SELECT CASE(SIZE(A,1))
    CASE(1)
      IF(ABS(A(1,1))<ZERO_TOLERANCE_DP) CALL FlagError("The A matrix is zero.",err,error,*999)
      x(1)=b(1)/A(1,1)
    CASE(2:3)
      CALL Invert(A,AInv,det,err,error,*999)
      x=MATMUL(AInv,b)
    CASE DEFAULT
      !Use LAPACK LU
      CALL DGETRF(SIZE(A,1),SIZE(A,2),A,SIZE(A,1),iPivot,info)
      IF(info/=0) THEN
        IF(info<0) THEN
          localError="Parameter number "//TRIM(NumberToVString(-info,"*",err,error))//" is illegal."
        ELSE
          localError="The A matrix is singular. The factorised diagonal entry for row "// &
            & TRIM(NumberToVString(info,"*",err,error))//" is zero."          
        ENDIF
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Solve the LU problem
      CALL DGETRS('N',1,A,SIZE(A,1),iPivot,b,SIZE(b,1),info)
      IF(info/=0) THEN
        localError="Parameter number "//TRIM(NumberToVString(-info,"*",err,error))//" is illegal."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      x=b
    END SELECT

    EXITS("SolveSmallLinearSystemDP")
    RETURN
999 ERRORSEXITS("SolveSmallLinearSystemDP",err,error)
    RETURN 1
    
  END SUBROUTINE SolveSmallLinearSystemDP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the matrix-matrix tensor product of the single precision matrices AxoB in C.
  SUBROUTINE TensorProductMatrix2Matrix2SP(A,B,C,err,error,*)

    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The A matrix
    REAL(SP), INTENT(IN) :: B(:,:) !<The B matrix
    REAL(SP), INTENT(OUT) :: C(:,:,:,:) !<On exit, the matrix tensor product C=AxoB
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,k,l

    ENTERS("TensorProductMatrix2Matrix2SP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(B,1).OR.SIZE(A,2)/=SIZE(B,2)) &
      & CALL FlagError("The B matrix is not the same size as the A matrix.",err,error,*999)
    IF(SIZE(C,1)/=SIZE(A,1).OR.SIZE(C,2)/=SIZE(A,1).OR.SIZE(C,3)/=SIZE(A,1).OR.SIZE(C,4)/=SIZE(A,1)) &
      & CALL FlagError("The C matrix is not of the same dimension as the A and B matrices.",err,error,*999)
#endif    
    
    DO i=1,SIZE(A,1)
      DO j=1,SIZE(A,2)
        DO k=1,SIZE(B,1)
          DO l=1,SIZE(B,2)
            C(i,j,k,l)=A(i,j)*B(k,l)
          ENDDO !l
        ENDDO !k
      ENDDO !j
    ENDDO !i
    
    EXITS("TensorProductMatrix2Matrix2SP")
    RETURN
999 ERRORSEXITS("TensorProductMatrix2Matrix2SP",err,error)
    RETURN 1
    
  END SUBROUTINE TensorProductMatrix2Matrix2SP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the matrix-matrix tensor product of the double precision matrices AxoB in C.
  SUBROUTINE TensorProductMatrix2Matrix2DP(A,B,C,err,error,*)

    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The A matrix
    REAL(DP), INTENT(IN) :: B(:,:) !<The B matrix
    REAL(DP), INTENT(OUT) :: C(:,:,:,:) !<On exit, the matrix tensor product C=AxoB
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,k,l

    ENTERS("TensorProductMatrix2Matrix2DP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(B,1).OR.SIZE(A,2)/=SIZE(B,2)) &
      & CALL FlagError("The B matrix is not the same size as the A matrix.",err,error,*999)
    IF(SIZE(C,1)/=SIZE(A,1).OR.SIZE(C,2)/=SIZE(A,1).OR.SIZE(C,3)/=SIZE(A,1).OR.SIZE(C,4)/=SIZE(A,1)) &
      & CALL FlagError("The C matrix is not of the same dimension as the A and B matrices.",err,error,*999)
#endif    
    
    DO i=1,SIZE(A,1)
      DO j=1,SIZE(A,2)
        DO k=1,SIZE(B,1)
          DO l=1,SIZE(B,2)
            C(i,j,k,l)=A(i,j)*B(k,l)
          ENDDO !l
        ENDDO !k
      ENDDO !j
    ENDDO !i
    
    EXITS("TensorProductMatrix2Matrix2DP")
    RETURN
999 ERRORSEXITS("TensorProductMatrix2Matrix2DP",err,error)
    RETURN 1
    
  END SUBROUTINE TensorProductMatrix2Matrix2DP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the vector-vector tensor product of the single precision vectors axob in C.
  SUBROUTINE TensorProductVectorVectorSP(a,b,C,err,error,*)

    !Argument variables
    REAL(SP), INTENT(IN) :: a(:) !<The a vector
    REAL(SP), INTENT(IN) :: b(:) !<The b vector
    REAL(SP), INTENT(OUT) :: C(:,:) !<On exit, the vector-vector tensor product C=axob
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j

    ENTERS("TensorProductVectorVectorSP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(C,1)/=SIZE(a,1).OR.SIZE(C,2)/=SIZE(b,1)) &
      & CALL FlagError("The C matrix is not of the same dimension as the a and b vectors.",err,error,*999)
#endif    
    
    DO i=1,SIZE(a,1)
      DO j=1,SIZE(b,1)
        C(i,j)=a(i)*b(j)
      ENDDO !j
    ENDDO !i
   
    EXITS("TensorProductVectorVectorSP")
    RETURN
999 ERRORSEXITS("TensorProductVectorVectorSP",err,error)
    RETURN 1
    
  END SUBROUTINE TensorProductVectorVectorSP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the vector-vector tensor product of the double precision vectors axob in C.
  SUBROUTINE TensorProductVectorVectorDP(a,b,C,err,error,*)

    !Argument variables
    REAL(DP), INTENT(IN) :: a(:) !<The a vector
    REAL(DP), INTENT(IN) :: b(:) !<The b vector
    REAL(DP), INTENT(OUT) :: C(:,:) !<On exit, the vector-vector tensor product C=axob
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j

    ENTERS("TensorProductVectorVectorDP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(C,1)/=SIZE(a,1).OR.SIZE(C,2)/=SIZE(b,1)) &
      & CALL FlagError("The C matrix is not of the same dimension as the a and b vectors.",err,error,*999)
#endif    
    
    DO i=1,SIZE(a,1)
      DO j=1,SIZE(b,1)
        C(i,j)=a(i)*b(j)
      ENDDO !j
    ENDDO !i
   
    EXITS("TensorProductVectorVectorDP")
    RETURN
999 ERRORSEXITS("TensorProductVectorVectorDP",err,error)
    RETURN 1
    
  END SUBROUTINE TensorProductVectorVectorDP

  !
  !================================================================================================================================
  !

  !>Converts a single precision second order tensor into a Voigt matrix form of the tensor. 
  SUBROUTINE TensorToVoigt2SP(n,tensorIndexTypes,A,AV,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: n !<The number of dimensions
    INTEGER(INTG), INTENT(IN) :: tensorIndexTypes(:) !<tensorIndexTypes(rankIdx). The index type of the rankIdx'th rank of the second order tensor \see Constants_TensorIndexTypes
    REAL(SP), INTENT(IN) :: A(:,:) !<A(iIdx,jIdx). The second order tensor to convert to Voigt form.
    
    REAL(SP), INTENT(OUT) :: AV(:) !<AV(voigtIdx). On exit, the Voigt form of the second order tensor.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("TensorToVoigt2SP",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(n<1.OR.n>3) THEN
      localError="The specified number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))// &
        & " is invalid. The number of dimensions should be >= 1 and <= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(tensorIndexTypes,1)<2) THEN
      localError="The size of the tensor index types array of "// &
        & TRIM(NumberToVString(SIZE(tensorIndexTypes,1),"*",err,error))// &
        & " is too small. The size should be >= 2 for a second order tensor."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,1)<n) THEN
      localError="The size of the first index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,1),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,2)<n) THEN
      localError="The size of the second index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,2),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(AV,1)<NUMBER_OF_VOIGT(n)) THEN
      localError="The size of the Voigt array of "//TRIM(NumberToVString(SIZE(AV,1),"*",err,error))// &
        & " is too small for a second order Voigt array with "//TRIM(NumberToVString(n,"*",err,error))// &
        & " dimensions. The size of the Voigt array should be >= "//TRIM(NumberToVString(NUMBER_OF_VOIGT(n),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    SELECT CASE(n)
    CASE(1)
      AV(1)=A(VOIGT_TO_TENSOR1(1,1),VOIGT_TO_TENSOR1(2,1))
    CASE(2)
      AV(1)=A(VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(2,1))
      AV(2)=A(VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(2,2))
      SELECT CASE(tensorIndexTypes(1))
      CASE(TENSOR_CONTRAVARIANT_INDEX)
        SELECT CASE(tensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Tensor is (contravariant,contravariant)
          AV(3)=A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,3))
        CASE(TENSOR_COVARIANT_INDEX)
          !Tensor is (contravariant,covariant). Don't apply any correction factors to mixed second order tensors.
          AV(3)=A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,3))
        CASE DEFAULT
          localError="The second tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(TENSOR_COVARIANT_INDEX)
        SELECT CASE(tensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Tensor is (covariant,contravariant). Don't apply any correction factors to mixed second order tensors.
          AV(3)=A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,3))
        CASE(TENSOR_COVARIANT_INDEX)
          !Tensor is (covariant,covariant)
          AV(3)=2.0_SP*A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,3))
        CASE DEFAULT
          localError="The second tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The first tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(1),"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(3)
      AV(1)=A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,1))
      AV(2)=A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,2))
      AV(3)=A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,3))
      SELECT CASE(tensorIndexTypes(1))
      CASE(TENSOR_CONTRAVARIANT_INDEX)
        SELECT CASE(tensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Tensor is (contravariant,contravariant)
          AV(4)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,4))
          AV(5)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,5))
          AV(6)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,6))
        CASE(TENSOR_COVARIANT_INDEX)
          !Tensor is (contravariant,covariant). Don't apply any correction factors to mixed second order tensors.
          AV(4)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,4))
          AV(5)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,5))
          AV(6)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,6))
        CASE DEFAULT
          localError="The second tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(TENSOR_COVARIANT_INDEX)
        SELECT CASE(tensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Tensor is (covariant,contravariant). Don't apply any correction factors to mixed second order tensors.
          AV(4)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,4))
          AV(5)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,5))
          AV(6)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,6))
        CASE(TENSOR_COVARIANT_INDEX)
          !Tensor is (covariant,covariant)
          AV(4)=2.0_SP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,4))
          AV(5)=2.0_SP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,5))
          AV(6)=2.0_SP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,6))
        CASE DEFAULT
          localError="The second tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The first tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(1),"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT    

    EXITS("TensorToVoigt2SP")
    RETURN
999 ERRORSEXITS("TensorToVoigt2SP",err,error)
    RETURN 1
    
  END SUBROUTINE TensorToVoigt2SP

  !
  !================================================================================================================================
  !

  !>Converts a double precision second order tensor into a Voigt matrix form of the tensor. 
  SUBROUTINE TensorToVoigt2DP(n,tensorIndexTypes,A,AV,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: n !<The number of dimensions
    INTEGER(INTG), INTENT(IN) :: tensorIndexTypes(:) !<tensorIndexTypes(rankIdx). The index type of the rankIdx'th rank of the second order tensor \see Constants_TensorIndexTypes
    REAL(DP), INTENT(IN) :: A(:,:) !<A(iIdx,jIdx). The second order tensor to convert to Voigt form.
    
    REAL(DP), INTENT(OUT) :: AV(:) !<AV(voigtIdx). On exit, the Voigt form of the second order tensor.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("TensorToVoigt2DP",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(n<1.OR.n>3) THEN
      localError="The specified number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))// &
        & " is invalid. The number of dimensions should be >= 1 and <= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(tensorIndexTypes,1)<2) THEN
      localError="The size of the tensor index types array of "// &
        & TRIM(NumberToVString(SIZE(tensorIndexTypes,1),"*",err,error))// &
        & " is too small. The size should be >= 2 for a second order tensor."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,1)<n) THEN
      localError="The size of the first index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,1),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,2)<n) THEN
      localError="The size of the second index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,2),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(AV,1)<NUMBER_OF_VOIGT(n)) THEN
      localError="The size of the Voigt array of "//TRIM(NumberToVString(SIZE(AV,1),"*",err,error))// &
        & " is too small for a second order Voigt array with "//TRIM(NumberToVString(n,"*",err,error))// &
        & " dimensions. The size of the Voigt array should be >= "//TRIM(NumberToVString(NUMBER_OF_VOIGT(n),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    SELECT CASE(n)
    CASE(1)
      AV(1)=A(VOIGT_TO_TENSOR1(1,1),VOIGT_TO_TENSOR1(2,1))
    CASE(2)
      AV(1)=A(VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(2,1))
      AV(2)=A(VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(2,2))
      SELECT CASE(tensorIndexTypes(1))
      CASE(TENSOR_CONTRAVARIANT_INDEX)
        SELECT CASE(tensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Tensor is (contravariant,contravariant)
          AV(3)=A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,3))
        CASE(TENSOR_COVARIANT_INDEX)
          !Tensor is (contravariant,covariant). Don't apply any correction factors to mixed second order tensors.
          AV(3)=A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,3))
        CASE DEFAULT
          localError="The second tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(TENSOR_COVARIANT_INDEX)
        SELECT CASE(tensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Tensor is (covariant,contravariant). Don't apply any correction factors to mixed second order tensors.
          AV(3)=A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,3))
        CASE(TENSOR_COVARIANT_INDEX)
          !Tensor is (covariant,covariant)
          AV(3)=2.0_DP*A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,3))
        CASE DEFAULT
          localError="The second tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The first tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(1),"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(3)
      AV(1)=A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,1))
      AV(2)=A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,2))
      AV(3)=A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,3))
      SELECT CASE(tensorIndexTypes(1))
      CASE(TENSOR_CONTRAVARIANT_INDEX)
        SELECT CASE(tensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Tensor is (contravariant,contravariant)
          AV(4)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,4))
          AV(5)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,5))
          AV(6)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,6))
        CASE(TENSOR_COVARIANT_INDEX)
          !Tensor is (contravariant,covariant). Don't apply any correction factors to mixed second order tensors.
          AV(4)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,4))
          AV(5)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,5))
          AV(6)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,6))
        CASE DEFAULT
          localError="The second tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(TENSOR_COVARIANT_INDEX)
        SELECT CASE(tensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Tensor is (covariant,contravariant). Don't apply any correction factors to mixed second order tensors.
          AV(4)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,4))
          AV(5)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,5))
          AV(6)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,6))
        CASE(TENSOR_COVARIANT_INDEX)
          !Tensor is (covariant,covariant)
          AV(4)=2.0_DP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,4))
          AV(5)=2.0_DP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,5))
          AV(6)=2.0_DP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,6))
        CASE DEFAULT
          localError="The second tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The first tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(1),"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT    

    EXITS("TensorToVoigt2DP")
    RETURN
999 ERRORSEXITS("TensorToVoigt2DP",err,error)
    RETURN 1
    
  END SUBROUTINE TensorToVoigt2DP

  !
  !================================================================================================================================
  !

  !>Converts a single precision fourth order tensor into a Voigt matrix form of the tensor. 
  SUBROUTINE TensorToVoigt4SP(n,tensorIndexTypes,A,AV,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: n !<The number of dimensions
    INTEGER(INTG), INTENT(IN) :: tensorIndexTypes(:) !<tensorIndexTypes(rankIdx). The index type of the rankIdx'th rank of the fourth order tensor \see Constants_TensorIndexTypes
    REAL(SP), INTENT(IN) :: A(:,:,:,:) !<A(iIdx,jIdx,kIdx,lIdx). The fourth order tensor to convert to Voigt form.
    
    REAL(SP), INTENT(OUT) :: AV(:,:) !<AV(voigtIdx,voigtIdx). On exit, the Voigt form of the fourth order tensor.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("TensorToVoigt4SP",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(n<1.OR.n>3) THEN
      localError="The specified number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))// &
        & " is invalid. The number of dimensions should be >= 1 and <= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(tensorIndexTypes,1)<4) THEN
      localError="The size of the tensor index types array of "// &
        & TRIM(NumberToVString(SIZE(tensorIndexTypes,1),"*",err,error))// &
        & " is too small. The size should be >= 4 for a fourth order tensor."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,1)<n) THEN
      localError="The size of the first index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,1),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,2)<n) THEN
      localError="The size of the second index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,2),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,3)<n) THEN
      localError="The size of the third index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,3),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,4)<n) THEN
      localError="The size of the fourth index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,4),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(AV,1)<NUMBER_OF_VOIGT(n)) THEN
      localError="The size of the first index of the Voigt array of "//TRIM(NumberToVString(SIZE(AV,1),"*",err,error))// &
        & " is too small for a fourth order Voigt array with "//TRIM(NumberToVString(n,"*",err,error))// &
        & " dimensions. The size of the Voigt array should be >= "//TRIM(NumberToVString(NUMBER_OF_VOIGT(n),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(AV,2)<NUMBER_OF_VOIGT(n)) THEN
      localError="The size of the second index of the Voigt array of "//TRIM(NumberToVString(SIZE(AV,2),"*",err,error))// &
        & " is too small for a fourth order Voigt array with "//TRIM(NumberToVString(n,"*",err,error))// &
        & " dimensions. The size of the Voigt array should be >= "//TRIM(NumberToVString(NUMBER_OF_VOIGT(n),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    SELECT CASE(n)
    CASE(1)
      AV(1,1)=A(VOIGT_TO_TENSOR1(1,1),VOIGT_TO_TENSOR1(1,1),VOIGT_TO_TENSOR1(2,1),VOIGT_TO_TENSOR1(2,1))
    CASE(2)
      AV(1,1)=A(VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(2,1),VOIGT_TO_TENSOR2(2,1))
      AV(2,1)=A(VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(2,1),VOIGT_TO_TENSOR2(2,1))
      AV(1,2)=A(VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(2,2),VOIGT_TO_TENSOR2(2,2))
      AV(2,2)=A(VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(2,2),VOIGT_TO_TENSOR2(2,2))
      SELECT CASE(tensorIndexTypes(1))
      CASE(TENSOR_CONTRAVARIANT_INDEX)
        SELECT CASE(tensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Tensor is (contravariant,contravariant)
          AV(3,1)=A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,1),VOIGT_TO_TENSOR2(2,1))
          AV(3,2)=A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,2),VOIGT_TO_TENSOR2(2,2))
          AV(1,3)=A(VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
          AV(2,3)=A(VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
          AV(3,3)=A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
        CASE(TENSOR_COVARIANT_INDEX)
          !Tensor is (contravariant,covariant).
          AV(3,1)=A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,1),VOIGT_TO_TENSOR2(2,1))
          AV(3,2)=A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,2),VOIGT_TO_TENSOR2(2,2))
          AV(1,3)=2.0_SP*A(VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
          AV(2,3)=2.0_SP*A(VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
          AV(3,3)=2.0_SP*A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
        CASE DEFAULT
          localError="The second tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(TENSOR_COVARIANT_INDEX)
        SELECT CASE(tensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Tensor is (covariant,contravariant). Don't apply any correction factors to mixed second order tensors.
          AV(3,1)=2.0_SP*A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,1),VOIGT_TO_TENSOR2(2,1))
          AV(3,2)=2.0_SP*A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,2),VOIGT_TO_TENSOR2(2,2))
          AV(1,3)=A(VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
          AV(2,3)=A(VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
          AV(3,3)=2.0_SP*A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
        CASE(TENSOR_COVARIANT_INDEX)
          !Tensor is (covariant,covariant)
          AV(3,1)=2.0_SP*A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,1),VOIGT_TO_TENSOR2(2,1))
          AV(3,2)=2.0_SP*A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,2),VOIGT_TO_TENSOR2(2,2))
          AV(1,3)=2.0_SP*A(VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
          AV(2,3)=2.0_SP*A(VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
          AV(3,3)=4.0_SP*A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
        CASE DEFAULT
          localError="The second tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The first tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(1),"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(3)
      AV(1,1)=A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
      AV(2,1)=A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
      AV(3,1)=A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
      AV(1,2)=A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
      AV(2,2)=A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
      AV(3,2)=A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
      AV(1,3)=A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
      AV(2,3)=A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
      AV(3,3)=A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
      SELECT CASE(tensorIndexTypes(1))
      CASE(TENSOR_CONTRAVARIANT_INDEX)
        SELECT CASE(tensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Tensor is (contravariant,contravariant)
          AV(4,1)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(5,1)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(6,1)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(4,2)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(5,2)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(6,2)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(4,3)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(5,3)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(6,3)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(1,4)=A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(2,4)=A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(3,4)=A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(4,4)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(5,4)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(6,4)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(1,5)=A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(2,5)=A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(3,5)=A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(4,5)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(5,5)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(6,5)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(1,6)=A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(2,6)=A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(3,6)=A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(4,6)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(5,6)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(6,6)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
        CASE(TENSOR_COVARIANT_INDEX)
          !Tensor is (contravariant,covariant). Don't apply any correction factors to mixed second order tensors.
          AV(4,1)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(5,1)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(6,1)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(4,2)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(5,2)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(6,2)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(4,3)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(5,3)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(6,3)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(1,4)=2.0_SP*A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(2,4)=2.0_SP*A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(3,4)=2.0_SP*A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(4,4)=2.0_SP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(5,4)=2.0_SP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(6,4)=2.0_SP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(1,5)=2.0_SP*A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(2,5)=2.0_SP*A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(3,5)=2.0_SP*A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(4,5)=2.0_SP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(5,5)=2.0_SP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(6,5)=2.0_SP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(1,6)=2.0_SP*A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(2,6)=2.0_SP*A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(3,6)=2.0_SP*A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(4,6)=2.0_SP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(5,6)=2.0_SP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(6,6)=2.0_SP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
        CASE DEFAULT
          localError="The second tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(TENSOR_COVARIANT_INDEX)
        SELECT CASE(tensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Tensor is (covariant,contravariant). Don't apply any correction factors to mixed second order tensors.
          AV(4,1)=2.0_SP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(5,1)=2.0_SP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(6,1)=2.0_SP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(4,2)=2.0_SP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(5,2)=2.0_SP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(6,2)=2.0_SP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(4,3)=2.0_SP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(5,3)=2.0_SP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(6,3)=2.0_SP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(1,4)=A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(2,4)=A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(3,4)=A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(4,4)=2.0_SP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(5,4)=2.0_SP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(6,4)=2.0_SP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(1,5)=A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(2,5)=A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(3,5)=A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(4,5)=2.0_SP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(5,5)=2.0_SP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(6,5)=2.0_SP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(1,6)=A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(2,6)=A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(3,6)=A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(4,6)=2.0_SP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(5,6)=2.0_SP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(6,6)=2.0_SP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
        CASE(TENSOR_COVARIANT_INDEX)
          !Tensor is (covariant,covariant)
          AV(4,1)=2.0_SP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(5,1)=2.0_SP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(6,1)=2.0_SP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(4,2)=2.0_SP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(5,2)=2.0_SP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(6,2)=2.0_SP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(4,3)=2.0_SP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(5,3)=2.0_SP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(6,3)=2.0_SP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(1,4)=2.0_SP*A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(2,4)=2.0_SP*A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(3,4)=2.0_SP*A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(4,4)=4.0_SP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(5,4)=4.0_SP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(6,4)=4.0_SP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(1,5)=2.0_SP*A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(2,5)=2.0_SP*A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(3,5)=2.0_SP*A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(4,5)=4.0_SP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(5,5)=4.0_SP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(6,5)=4.0_SP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(1,6)=2.0_SP*A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(2,6)=2.0_SP*A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(3,6)=2.0_SP*A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(4,6)=4.0_SP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(5,6)=4.0_SP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(6,6)=4.0_SP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
        CASE DEFAULT
          localError="The second tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The first tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(1),"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT    

    EXITS("TensorToVoigt4SP")
    RETURN
999 ERRORSEXITS("TensorToVoigt4SP",err,error)
    RETURN 1
    
  END SUBROUTINE TensorToVoigt4SP

  !
  !================================================================================================================================
  !

  !>Converts a double precision fourth order tensor into a Voigt matrix form of the tensor. 
  SUBROUTINE TensorToVoigt4DP(n,tensorIndexTypes,A,AV,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: n !<The number of dimensions
    INTEGER(INTG), INTENT(IN) :: tensorIndexTypes(:) !<tensorIndexTypes(rankIdx). The index type of the rankIdx'th rank of the fourth order tensor \see Constants_TensorIndexTypes
    REAL(DP), INTENT(IN) :: A(:,:,:,:) !<A(iIdx,jIdx,kIdx,lIdx). The fourth order tensor to convert to Voigt form.
    
    REAL(DP), INTENT(OUT) :: AV(:,:) !<AV(voigtIdx,voigtIdx). On exit, the Voigt form of the fourth order tensor.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("TensorToVoigt4DP",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(n<1.OR.n>3) THEN
      localError="The specified number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))// &
        & " is invalid. The number of dimensions should be >= 1 and <= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(tensorIndexTypes,1)<4) THEN
      localError="The size of the tensor index types array of "// &
        & TRIM(NumberToVString(SIZE(tensorIndexTypes,1),"*",err,error))// &
        & " is too small. The size should be >= 4 for a fourth order tensor."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,1)<n) THEN
      localError="The size of the first index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,1),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,2)<n) THEN
      localError="The size of the second index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,2),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,3)<n) THEN
      localError="The size of the third index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,3),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,4)<n) THEN
      localError="The size of the fourth index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,4),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(AV,1)<NUMBER_OF_VOIGT(n)) THEN
      localError="The size of the first index of the Voigt array of "//TRIM(NumberToVString(SIZE(AV,1),"*",err,error))// &
        & " is too small for a fourth order Voigt array with "//TRIM(NumberToVString(n,"*",err,error))// &
        & " dimensions. The size of the Voigt array should be >= "//TRIM(NumberToVString(NUMBER_OF_VOIGT(n),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(AV,2)<NUMBER_OF_VOIGT(n)) THEN
      localError="The size of the second index of the Voigt array of "//TRIM(NumberToVString(SIZE(AV,2),"*",err,error))// &
        & " is too small for a fourth order Voigt array with "//TRIM(NumberToVString(n,"*",err,error))// &
        & " dimensions. The size of the Voigt array should be >= "//TRIM(NumberToVString(NUMBER_OF_VOIGT(n),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    SELECT CASE(n)
    CASE(1)
      AV(1,1)=A(VOIGT_TO_TENSOR1(1,1),VOIGT_TO_TENSOR1(1,1),VOIGT_TO_TENSOR1(2,1),VOIGT_TO_TENSOR1(2,1))
    CASE(2)
      AV(1,1)=A(VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(2,1),VOIGT_TO_TENSOR2(2,1))
      AV(2,1)=A(VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(2,1),VOIGT_TO_TENSOR2(2,1))
      AV(1,2)=A(VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(2,2),VOIGT_TO_TENSOR2(2,2))
      AV(2,2)=A(VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(2,2),VOIGT_TO_TENSOR2(2,2))
      SELECT CASE(tensorIndexTypes(1))
      CASE(TENSOR_CONTRAVARIANT_INDEX)
        SELECT CASE(tensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Tensor is (contravariant,contravariant)
          AV(3,1)=A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,1),VOIGT_TO_TENSOR2(2,1))
          AV(3,2)=A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,2),VOIGT_TO_TENSOR2(2,2))
          AV(1,3)=A(VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
          AV(2,3)=A(VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
          AV(3,3)=A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
        CASE(TENSOR_COVARIANT_INDEX)
          !Tensor is (contravariant,covariant).
          AV(3,1)=A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,1),VOIGT_TO_TENSOR2(2,1))
          AV(3,2)=A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,2),VOIGT_TO_TENSOR2(2,2))
          AV(1,3)=2.0_DP*A(VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
          AV(2,3)=2.0_DP*A(VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
          AV(3,3)=2.0_DP*A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
        CASE DEFAULT
          localError="The second tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(TENSOR_COVARIANT_INDEX)
        SELECT CASE(tensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Tensor is (covariant,contravariant). Don't apply any correction factors to mixed second order tensors.
          AV(3,1)=2.0_DP*A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,1),VOIGT_TO_TENSOR2(2,1))
          AV(3,2)=2.0_DP*A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,2),VOIGT_TO_TENSOR2(2,2))
          AV(1,3)=A(VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
          AV(2,3)=A(VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
          AV(3,3)=2.0_DP*A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
        CASE(TENSOR_COVARIANT_INDEX)
          !Tensor is (covariant,covariant)
          AV(3,1)=2.0_DP*A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,1),VOIGT_TO_TENSOR2(2,1))
          AV(3,2)=2.0_DP*A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,2),VOIGT_TO_TENSOR2(2,2))
          AV(1,3)=2.0_DP*A(VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(1,1),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
          AV(2,3)=2.0_DP*A(VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(1,2),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
          AV(3,3)=4.0_DP*A(VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(1,3),VOIGT_TO_TENSOR2(2,3),VOIGT_TO_TENSOR2(2,3))
        CASE DEFAULT
          localError="The second tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The first tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(1),"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(3)
      AV(1,1)=A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
      AV(2,1)=A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
      AV(3,1)=A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
      AV(1,2)=A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
      AV(2,2)=A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
      AV(3,2)=A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
      AV(1,3)=A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
      AV(2,3)=A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
      AV(3,3)=A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
      SELECT CASE(tensorIndexTypes(1))
      CASE(TENSOR_CONTRAVARIANT_INDEX)
        SELECT CASE(tensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Tensor is (contravariant,contravariant)
          AV(4,1)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(5,1)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(6,1)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(4,2)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(5,2)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(6,2)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(4,3)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(5,3)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(6,3)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(1,4)=A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(2,4)=A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(3,4)=A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(4,4)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(5,4)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(6,4)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(1,5)=A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(2,5)=A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(3,5)=A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(4,5)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(5,5)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(6,5)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(1,6)=A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(2,6)=A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(3,6)=A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(4,6)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(5,6)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(6,6)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
        CASE(TENSOR_COVARIANT_INDEX)
          !Tensor is (contravariant,covariant). Don't apply any correction factors to mixed second order tensors.
          AV(4,1)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(5,1)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(6,1)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(4,2)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(5,2)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(6,2)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(4,3)=A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(5,3)=A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(6,3)=A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(1,4)=2.0_DP*A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(2,4)=2.0_DP*A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(3,4)=2.0_DP*A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(4,4)=2.0_DP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(5,4)=2.0_DP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(6,4)=2.0_DP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(1,5)=2.0_DP*A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(2,5)=2.0_DP*A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(3,5)=2.0_DP*A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(4,5)=2.0_DP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(5,5)=2.0_DP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(6,5)=2.0_DP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(1,6)=2.0_DP*A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(2,6)=2.0_DP*A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(3,6)=2.0_DP*A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(4,6)=2.0_DP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(5,6)=2.0_DP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(6,6)=2.0_DP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
        CASE DEFAULT
          localError="The second tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(TENSOR_COVARIANT_INDEX)
        SELECT CASE(tensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Tensor is (covariant,contravariant). Don't apply any correction factors to mixed second order tensors.
          AV(4,1)=2.0_DP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(5,1)=2.0_DP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(6,1)=2.0_DP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(4,2)=2.0_DP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(5,2)=2.0_DP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(6,2)=2.0_DP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(4,3)=2.0_DP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(5,3)=2.0_DP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(6,3)=2.0_DP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(1,4)=A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(2,4)=A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(3,4)=A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(4,4)=2.0_DP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(5,4)=2.0_DP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(6,4)=2.0_DP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(1,5)=A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(2,5)=A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(3,5)=A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(4,5)=2.0_DP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(5,5)=2.0_DP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(6,5)=2.0_DP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(1,6)=A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(2,6)=A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(3,6)=A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(4,6)=2.0_DP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(5,6)=2.0_DP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(6,6)=2.0_DP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
        CASE(TENSOR_COVARIANT_INDEX)
          !Tensor is (covariant,covariant)
          AV(4,1)=2.0_DP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(5,1)=2.0_DP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(6,1)=2.0_DP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,1),VOIGT_TO_TENSOR3(2,1))
          AV(4,2)=2.0_DP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(5,2)=2.0_DP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(6,2)=2.0_DP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,2),VOIGT_TO_TENSOR3(2,2))
          AV(4,3)=2.0_DP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(5,3)=2.0_DP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(6,3)=2.0_DP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,3),VOIGT_TO_TENSOR3(2,3))
          AV(1,4)=2.0_DP*A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(2,4)=2.0_DP*A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(3,4)=2.0_DP*A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(4,4)=4.0_DP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(5,4)=4.0_DP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(6,4)=4.0_DP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,4),VOIGT_TO_TENSOR3(2,4))
          AV(1,5)=2.0_DP*A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(2,5)=2.0_DP*A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(3,5)=2.0_DP*A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(4,5)=4.0_DP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(5,5)=4.0_DP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(6,5)=4.0_DP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,5),VOIGT_TO_TENSOR3(2,5))
          AV(1,6)=2.0_DP*A(VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(1,1),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(2,6)=2.0_DP*A(VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(1,2),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(3,6)=2.0_DP*A(VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(1,3),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(4,6)=4.0_DP*A(VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(1,4),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(5,6)=4.0_DP*A(VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(1,5),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
          AV(6,6)=4.0_DP*A(VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(1,6),VOIGT_TO_TENSOR3(2,6),VOIGT_TO_TENSOR3(2,6))
        CASE DEFAULT
          localError="The second tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The first tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(1),"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT    

    EXITS("TensorToVoigt4DP")
    RETURN
999 ERRORSEXITS("TensorToVoigt4DP",err,error)
    RETURN 1
    
  END SUBROUTINE TensorToVoigt4DP

  !
  !================================================================================================================================
  !

  !>Transforms a vector, a(nu), in nu coordinates to a vector, b(xi), in xi coordinates. 
  SUBROUTINE TensorTransformVectorDP0(n,vectorIndexType,a,dXidNu,dNudXi,b,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: n !<The number of dimensions
    INTEGER(INTG), INTENT(IN) :: vectorIndexType !<The type (e.g., covariant or contravariant) of the vector index to transform \see Constants_TensorIndexTypes
    REAL(DP), INTENT(IN) :: a(:) !<a(nuIdx). The orginal vector a in nu coordinates to transform to xi coordinates.
    REAL(DP), INTENT(IN) :: dXidNu(:,:) !<dXidNu(xiIdx,nuIdx). The rate of change matrix of xi coordinates with nu coordinates.
    REAL(DP), INTENT(IN) :: dNudXi(:,:) !<dNudXi(nuIdx,xiIdx). The rate of change matrix of nu coordinates with xi coordinates.
    REAL(DP), INTENT(OUT) :: b(:) !<b(xiIdx). On exit, the transformed vector b in xi coordinates.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("TensorTransformVectorDP0",err,error,*999)

    CALL TensorTransformVectorDP1(n,[vectorIndexType],a,dXidNu,dNudXi,b,err,error,*999)
       
    EXITS("TensorTransformVectorDP0")
    RETURN
999 ERRORS("TensorTransformVectorDP0",err,error)
    EXITS("TensorTransformVectorDP0")
    RETURN 1
    
  END SUBROUTINE TensorTransformVectorDP0

  !
  !================================================================================================================================
  !

  !>Transforms a vector, a(nu), in nu coordinates to a vector, b(xi), in xi coordinates. 
  SUBROUTINE TensorTransformVectorDP1(n,vectorIndexType,a,dXidNu,dNudXi,b,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: n !<The number of dimensions
    INTEGER(INTG), INTENT(IN) :: vectorIndexType(1) !<The type (e.g., covariant or contravariant) of the vector index to transform \see Constants_TensorIndexTypes
    REAL(DP), INTENT(IN) :: a(:) !<a(nuIdx). The orginal vector a in nu coordinates to transform to xi coordinates.
    REAL(DP), INTENT(IN) :: dXidNu(:,:) !<dXidNu(xiIdx,nuIdx). The rate of change matrix of xi coordinates with nu coordinates.
    REAL(DP), INTENT(IN) :: dNudXi(:,:) !<dNudXi(nuIdx,xiIdx). The rate of change matrix of nu coordinates with xi coordinates.
    REAL(DP), INTENT(OUT) :: b(:) !<b(xiIdx). On exit, the transformed vector b in xi coordinates.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("TensorTransformVectorDP1",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(n<1.OR.n>3) THEN
      localError="The specified number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))// &
        & " is invalid. The number of dimensions should be >= 1 and <= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(a,1)<n) THEN
      localError="The size of the vector a to transform of "//TRIM(NumberToVString(SIZE(a,1),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(dXidNu,1)<n) THEN
      localError="The size of the first index of the dXidNu array of "//TRIM(NumberToVString(SIZE(dXidNu,1),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(dXidNu,2)<n) THEN
      localError="The size of the second index of the dXidNu array of "//TRIM(NumberToVString(SIZE(dXidNu,2),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(dNudXi,1)<n) THEN
      localError="The size of the first index of the dNudXi array of "//TRIM(NumberToVString(SIZE(dNudXi,1),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(dNudXi,2)<n) THEN
      localError="The size of the second index of the dNudXi array of "//TRIM(NumberToVString(SIZE(dNudXi,2),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(b,1)<n) THEN
      localError="The size of the transformed vector b of "//TRIM(NumberToVString(SIZE(b,1),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    !Transform the vector a in nu coordinates to the vector b in xi coordinates according to the appropriate transformation rule
    SELECT CASE(vectorIndexType(1))
    CASE(TENSOR_CONTRAVARIANT_INDEX)
      !Contravariant vector i.e., b^{r}=\delby{\xi^{r}}{\nu^{a}}.a^{a}
      CALL MatrixVectorProduct(dXidNu(1:n,1:n),a(1:n),b(1:n),err,error,*999)
    CASE(TENSOR_COVARIANT_INDEX)
      !Convariant vector i.e., b_{r}=\delby{\nu^{a}}{\xi^{r}}.a_{a}
      CALL MatrixTransposeVectorProduct(dNudXi(1:n,1:n),a(1:n),b(1:n),err,error,*999)
    CASE DEFAULT
      localError="The vector index type of "//TRIM(NumberToVString(vectorIndexType(1),"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("TensorTransformVectorDP1")
    RETURN
999 ERRORS("TensorTransformVectorDP1",err,error)
    EXITS("TensorTransformVectorDP1")
    RETURN 1
    
  END SUBROUTINE TensorTransformVectorDP1

  !
  !================================================================================================================================
  !

  !>Transforms a second order tenso, A(nu1,nu2), in nu coordinates to a second order tensor, B(xi1,xi2), in xi coordinates. 
  SUBROUTINE TensorTransformTensorTwoDP(n,tensorIndexTypes,A,dXidNu,dNudXi,B,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: n !<The number of dimensions
    INTEGER(INTG), INTENT(IN) :: tensorIndexTypes(:) !<tensorIndexTypes(rankIdx). The type (e.g., covariant or contravariant) of the rankIdx'th index to transform \see Constants_TensorIndexTypes
    REAL(DP), INTENT(IN) :: A(:,:) !<A(nuIdx1,nuIdx2). The orginal second order tensor A in nu coordinates to transform to xi coordinates.
    REAL(DP), INTENT(IN) :: dXidNu(:,:) !<dXidNu(xiIdx,nuIdx). The rate of change matrix of xi coordinates with nu coordinates.
    REAL(DP), INTENT(IN) :: dNudXi(:,:) !<dNudXi(nuIdx,xiIdx). The rate of change matrix of nu coordinates with xi coordinates.
    REAL(DP), INTENT(OUT) :: B(:,:) !<B(xiIdx1,xiIdx2). On exit, the transformed second order tensor B in xi coordinates.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: tempTensor(3,3)
    TYPE(VARYING_STRING) :: localError

    ENTERS("TensorTransformTensorTwoDP",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(n<1.OR.n>3) THEN
      localError="The specified number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))// &
        & " is invalid. The number of dimensions should be >= 1 and <= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(tensorIndexTypes,1)<2) THEN
      localError="The size of the specified tensor index types array is "// &
        & NumberToVString(SIZE(tensorIndexTypes,1),"*",err,error)//" and needs to be at least two for a second order tensor."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,1)<n) THEN
      localError="The size of the first index of the tensor A to transform of "//TRIM(NumberToVString(SIZE(A,1),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,2)<n) THEN
      localError="The size of the second index of the tensor A to transform of "//TRIM(NumberToVString(SIZE(A,2),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(dXidNu,1)<n) THEN
      localError="The size of the first index of the dXidNu array of "//TRIM(NumberToVString(SIZE(dXidNu,1),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(dXidNu,2)<n) THEN
      localError="The size of the second index of the dXidNu array of "//TRIM(NumberToVString(SIZE(dXidNu,2),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(dNudXi,1)<n) THEN
      localError="The size of the first index of the dNudXi array of "//TRIM(NumberToVString(SIZE(dNudXi,1),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(dNudXi,2)<n) THEN
      localError="The size of the second index of the dNudXi array of "//TRIM(NumberToVString(SIZE(dNudXi,2),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(B,1)<n) THEN
      localError="The size of the first index of the transformed tensor B of "//TRIM(NumberToVString(SIZE(B,1),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(B,2)<n) THEN
      localError="The size of the second index of the transformed tensor B of "//TRIM(NumberToVString(SIZE(B,2),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    !Transform the tensor A in nu coordinates to the tensor B in xi coordinates according to the appropriate transformation rule
    SELECT CASE(tensorIndexTypes(1))
    CASE(TENSOR_CONTRAVARIANT_INDEX)
      SELECT CASE(tensorIndexTypes(2))
      CASE(TENSOR_CONTRAVARIANT_INDEX)
        !(Contravariant,Contravariant) tensor i.e.,
        !B^{rs}=\delby{\xi^{r}}{\nu^{a}}\delby{\xi^{s}}{\nu^{b}}.A^{ab}
        CALL MatrixProductTranspose(A(1:n,1:n),dXidNu(1:n,1:n),tempTensor(1:n,1:n),err,error,*999)
        CALL MatrixProduct(dXidNu(1:n,1:n),tempTensor(1:n,1:n),B(1:n,1:n),err,error,*999)
      CASE(TENSOR_COVARIANT_INDEX)
        !(Contravariant,Covariant) tensor i.e.,
        !B^{r.}_{.s}=\delby{\xi^{r}}{\nu^{a}}\delby{\nu^{b}}{\xi^{s}}.A^{a.}_{.b}
        CALL MatrixProduct(A(1:n,1:n),dNudXi(1:n,1:n),tempTensor(1:n,1:n),err,error,*999)
        CALL MatrixProduct(dXidNu(1:n,1:n),tempTensor(1:n,1:n),B(1:n,1:n),err,error,*999)
      CASE(TENSOR_TWO_POINT_CONTRAVARIANT_INDEX)
        !Two-point tensor. Transform only as contravariant in the second index.
      CASE(TENSOR_TWO_POINT_COVARIANT_INDEX)
        !Two-point tensor. Transform only as covariant in the second index.
      CASE DEFAULT
        localError="The second tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(2),"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(TENSOR_COVARIANT_INDEX)
      SELECT CASE(tensorIndexTypes(2))
      CASE(TENSOR_CONTRAVARIANT_INDEX)
        !(Covariant,Contravariant) tensor i.e.,
        !B^{.s}_{r.}=\delby{\nu^{a}}{\xi^{r}}\delby{\xi^{s}}{\nu^{b}}.A^{a.}_{.b}
        CALL MatrixProduct(A(1:n,1:n),dNudXi(1:n,1:n),tempTensor(1:n,1:n),err,error,*999)
        CALL MatrixTransposeProduct(dNudXi(1:n,1:n),tempTensor(1:n,1:n),B(1:n,1:n),err,error,*999)
      CASE(TENSOR_COVARIANT_INDEX)
        !(Covariant,Covariant) tensor i.e., 
        !B_{rs}=\delby{\nu^{a}}{\xi^{r}}\delby{\nu^{b}}{\xi^{s}}.A_{ab}
        CALL MatrixProduct(A(1:n,1:n),dNudXi(1:n,1:n),tempTensor(1:n,1:n),err,error,*999)
        CALL MatrixTransposeProduct(dNudXi(1:n,1:n),tempTensor(1:n,1:n),B(1:n,1:n),err,error,*999)
      CASE(TENSOR_TWO_POINT_CONTRAVARIANT_INDEX)
        !Two-point tensor. Transform only as contravariant in the second index i.e.,
        !B^{as}=\delby{\xi^{s}}{\nu^{b}}.A^{ab}
        CALL MatrixProductTranspose(A(1:n,1:n),dXidNu(1:n,1:n),B(1:n,1:n),err,error,*999)
      CASE(TENSOR_TWO_POINT_COVARIANT_INDEX)
        !Two-point tensor. Transform only as covariant in the second index i.e.,
        !B_{as}=\delby{\nu^{b}}{\xi^{s}}.A_{ab}
        CALL MatrixProduct(A(1:n,1:n),dNudXi(1:n,1:n),B(1:n,1:n),err,error,*999)
     CASE DEFAULT       
        localError="The second tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(2),"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(TENSOR_TWO_POINT_CONTRAVARIANT_INDEX)
      !Two-point tensor. Transform only as contravariant in the first index i.e.,
      !B^{rb}=\delby{\xi^{r}}{\nu^{a}}.A^{ab}
      SELECT CASE(tensorIndexTypes(2))
      CASE(TENSOR_CONTRAVARIANT_INDEX, &
        & TENSOR_COVARIANT_INDEX)
        CALL MatrixProduct(dNudXi(1:n,1:n),A(1:n,1:n),B(1:n,1:n),err,error,*999)
      CASE DEFAULT
        localError="The second tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(2),"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT       
    CASE(TENSOR_TWO_POINT_COVARIANT_INDEX)
      !Two-point tensor. Transform only as covariant in the first index i.e.,
      !B_{rb}=\delby{\nu^{a}}{\xi^{r}}.A_{ab}
      SELECT CASE(tensorIndexTypes(2))
      CASE(TENSOR_CONTRAVARIANT_INDEX, &
        & TENSOR_COVARIANT_INDEX)
        CALL MatrixTransposeProduct(dNudXi(1:n,1:n),A(1:n,1:n),B(1:n,1:n),err,error,*999)
      CASE DEFAULT
        localError="The second tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(2),"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT       
   CASE DEFAULT
      localError="The first tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(1),"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("TensorTransformTensorTwoDP")
    RETURN
999 ERRORS("TensorTransformTensorTwoDP",err,error)
    EXITS("TensorTransformTensorTwoDP")
    RETURN 1
    
  END SUBROUTINE TensorTransformTensorTwoDP

  !
  !================================================================================================================================
  !

  !>Transforms a fourth order tenso, A(nu1,nu2,nu3,nu4), in nu coordinates to a fourth order tensor, B(xi1,xi2,xi3,xi4), in xi coordinates. 
  SUBROUTINE TensorTransformTensorFourDP(n,tensorIndexTypes,A,dXidNu,dNudXi,B,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: n !<The number of dimensions
    INTEGER(INTG), INTENT(IN) :: tensorIndexTypes(:) !<tensorIndexTypes(rankIdx). The type (e.g., covariant or contravariant) of the rankIdx'th index to transform \see Constants_TensorIndexTypes
    REAL(DP), INTENT(IN) :: A(:,:,:,:) !<A(nuIdx1,nuIdx2,nuIdx3,nuIdx4). The orginal fourth order tensor A in nu coordinates to transform to xi coordinates.
    REAL(DP), INTENT(IN) :: dXidNu(:,:) !<dXidNu(xiIdx,nuIdx). The rate of change matrix of xi coordinates with nu coordinates.
    REAL(DP), INTENT(IN) :: dNudXi(:,:) !<dNudXi(nuIdx,xiIdx). The rate of change matrix of nu coordinates with xi coordinates.
    REAL(DP), INTENT(OUT) :: B(:,:,:,:) !<B(xiIdx1,xiIdx2,xiIdx3,xiIdx4). On exit, the transformed fourth order tensor B in xi coordinates.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: aIdx,bIdx,cIdx,dIdx,rIdx,sIdx,tIdx,uIdx
    REAL(DP) :: sum
    TYPE(VARYING_STRING) :: localError

    ENTERS("TensorTransformTensorFourDP",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(n<1.OR.n>3) THEN
      localError="The specified number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))// &
        & " is invalid. The number of dimensions should be >= 1 and <= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(tensorIndexTypes,1)<4) THEN
      localError="The size of the specified tensor index types array is "// &
        & NumberToVString(SIZE(tensorIndexTypes,1),"*",err,error)//" and needs to be at least four for a fourth order tensor."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,1)<n) THEN
      localError="The size of the first index of the tensor A to transform of "//TRIM(NumberToVString(SIZE(A,1),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,2)<n) THEN
      localError="The size of the second index of the tensor A to transform of "//TRIM(NumberToVString(SIZE(A,2),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,3)<n) THEN
      localError="The size of the third index of the tensor A to transform of "//TRIM(NumberToVString(SIZE(A,3),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,4)<n) THEN
      localError="The size of the fourth index of the tensor A to transform of "//TRIM(NumberToVString(SIZE(A,4),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(dXidNu,1)<n) THEN
      localError="The size of the first index of the dXidNu array of "//TRIM(NumberToVString(SIZE(dXidNu,1),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(dXidNu,2)<n) THEN
      localError="The size of the second index of the dXidNu array of "//TRIM(NumberToVString(SIZE(dXidNu,2),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(dNudXi,1)<n) THEN
      localError="The size of the first index of the dNudXi array of "//TRIM(NumberToVString(SIZE(dNudXi,1),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(dNudXi,2)<n) THEN
      localError="The size of the second index of the dNudXi array of "//TRIM(NumberToVString(SIZE(dNudXi,2),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(B,1)<n) THEN
      localError="The size of the first index of the transformed tensor B of "//TRIM(NumberToVString(SIZE(B,1),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(B,2)<n) THEN
      localError="The size of the second index of the transformed tensor B of "//TRIM(NumberToVString(SIZE(B,2),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(B,3)<n) THEN
      localError="The size of the third index of the transformed tensor B of "//TRIM(NumberToVString(SIZE(B,3),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(B,4)<n) THEN
      localError="The size of the fourth index of the transformed tensor B of "//TRIM(NumberToVString(SIZE(B,4),"*",err,error))// &
        & " is too small. The size should be >= the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    !Transform the tensor A in nu coordinates to the tensor B in xi coordinates according to the appropriate transformation rule
    SELECT CASE(tensorIndexTypes(1))
    CASE(TENSOR_CONTRAVARIANT_INDEX)
      SELECT CASE(tensorIndexTypes(2))
      CASE(TENSOR_CONTRAVARIANT_INDEX)
        SELECT CASE(tensorIndexTypes(3))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          SELECT CASE(tensorIndexTypes(4))
          CASE(TENSOR_CONTRAVARIANT_INDEX)
            !(Contravariant,Contravariant,Contravariant,Contravariant) tensor
            DO rIdx=1,n
              DO sIdx=1,n
                DO tIdx=1,n
                  DO uIdx=1,n
                    sum=0.0_DP
                    DO aIdx=1,n
                      DO bIdx=1,n
                        DO cIdx=1,n
                          DO dIdx=1,n
                            sum=sum+dXidNu(rIdx,aIdx)* &
                              & dXidNu(sIdx,bIdx)* &
                              & dXidNu(tIdx,cIdx)* &
                              & dXidNu(uIdx,dIdx)* &
                              & A(aIdx,bIdx,cIdx,dIdx)
                          ENDDO !dIdx
                        ENDDO !cIdx
                      ENDDO !bIdx
                    ENDDO !aIdx
                    B(rIdx,sIdx,tIdx,uIdx)=sum
                 ENDDO !uIdx
                ENDDO !tIdx
              ENDDO !sIdx
            ENDDO !rIdx
          CASE(TENSOR_COVARIANT_INDEX)
            !(Contravariant,Contravariant,Contravariant,Covariant) tensor
            DO rIdx=1,n
              DO sIdx=1,n
                DO tIdx=1,n
                  DO uIdx=1,n
                    sum=0.0_DP
                    DO aIdx=1,n
                      DO bIdx=1,n
                        DO cIdx=1,n
                          DO dIdx=1,n
                            sum=sum+dXidNu(rIdx,aIdx)* &
                              & dXidNu(sIdx,bIdx)* &
                              & dXidNu(tIdx,cIdx)* &
                              & dNudXi(dIdx,uIdx)* &
                              & A(aIdx,bIdx,cIdx,dIdx)
                          ENDDO !dIdx
                        ENDDO !cIdx
                      ENDDO !bIdx
                    ENDDO !aIdx
                    B(rIdx,sIdx,tIdx,uIdx)=sum
                  ENDDO !uIdx
                ENDDO !tIdx
              ENDDO !sIdx
            ENDDO !rIdx
          CASE DEFAULT
            localError="The fourth tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(4),"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(TENSOR_COVARIANT_INDEX)
          SELECT CASE(tensorIndexTypes(4))
          CASE(TENSOR_CONTRAVARIANT_INDEX)
            !(Contravariant,Contravariant,Covariant,Contravariant) tensor
            DO rIdx=1,n
              DO sIdx=1,n
                DO tIdx=1,n
                  DO uIdx=1,n
                    sum=0.0_DP
                    DO aIdx=1,n
                      DO bIdx=1,n
                        DO cIdx=1,n
                          DO dIdx=1,n
                            sum=sum+dXidNu(rIdx,aIdx)* &
                              & dXidNu(sIdx,bIdx)* &
                              & dNudXi(cIdx,tIdx)* &
                              & dXidNu(uIdx,dIdx)* &
                              & A(aIdx,bIdx,cIdx,dIdx)
                          ENDDO !dIdx
                        ENDDO !cIdx
                      ENDDO !bIdx
                    ENDDO !aIdx
                    B(rIdx,sIdx,tIdx,uIdx)=sum
                  ENDDO !uIdx
                ENDDO !tIdx
              ENDDO !sIdx
            ENDDO !rIdx
          CASE(TENSOR_COVARIANT_INDEX)
            !(Contravariant,Contravariant,Covariant,Covariant) tensor
            DO rIdx=1,n
              DO sIdx=1,n
                DO tIdx=1,n
                  DO uIdx=1,n
                    sum=0.0_DP
                    DO aIdx=1,n
                      DO bIdx=1,n
                        DO cIdx=1,n
                          DO dIdx=1,n
                            sum=sum+dXidNu(rIdx,aIdx)* &
                              & dXidNu(sIdx,bIdx)* &
                              & dNudXi(cIdx,tIdx)* &
                              & dNudXi(dIdx,uIdx)* &
                              & A(aIdx,bIdx,cIdx,dIdx)
                          ENDDO !dIdx
                        ENDDO !cIdx
                      ENDDO !bIdx
                   ENDDO !aIdx
                   B(rIdx,sIdx,tIdx,uIdx)=sum
                  ENDDO !uIdx
                ENDDO !tIdx
              ENDDO !sIdx
            ENDDO !rIdx
          CASE DEFAULT
            localError="The fourth tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(4),"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The third tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(3),"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(TENSOR_COVARIANT_INDEX)
        SELECT CASE(tensorIndexTypes(3))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          SELECT CASE(tensorIndexTypes(4))
          CASE(TENSOR_CONTRAVARIANT_INDEX)
            !(Contravariant,Covariant,Contravariant,Contravariant) tensor
            DO rIdx=1,n
              DO sIdx=1,n
                DO tIdx=1,n
                  DO uIdx=1,n
                    sum=0.0_DP
                    DO aIdx=1,n
                      DO bIdx=1,n
                        DO cIdx=1,n
                          DO dIdx=1,n
                            sum=sum+dXidNu(rIdx,aIdx)* &
                              & dNudXi(bIdx,sIdx)* &
                              & dXidNu(tIdx,cIdx)* &
                              & dXidNu(uIdx,dIdx)* &
                              & A(aIdx,bIdx,cIdx,dIdx)
                          ENDDO !dIdx
                        ENDDO !cIdx
                      ENDDO !bIdx
                    ENDDO !aIdx
                    B(rIdx,sIdx,tIdx,uIdx)=sum
                  ENDDO !uIdx
                ENDDO !tIdx
              ENDDO !sIdx
            ENDDO !rIdx
           CASE(TENSOR_COVARIANT_INDEX)
            !(Contravariant,Covariant,Contravariant,Covariant) tensor
            DO rIdx=1,n
              DO sIdx=1,n
                DO tIdx=1,n
                  DO uIdx=1,n
                    sum=0.0_DP
                    DO aIdx=1,n
                      DO bIdx=1,n
                        DO cIdx=1,n
                          DO dIdx=1,n
                            sum=sum+dXidNu(rIdx,aIdx)* &
                              & dNudXi(bIdx,sIdx)* &
                              & dXidNu(tIdx,cIdx)* &
                              & dNudXi(dIdx,uIdx)* &
                              & A(aIdx,bIdx,cIdx,dIdx)
                          ENDDO !dIdx
                        ENDDO !cIdx
                      ENDDO !bIdx
                    ENDDO !aIdx
                    B(rIdx,sIdx,tIdx,uIdx)=sum
                 ENDDO !uIdx
                ENDDO !tIdx
              ENDDO !sIdx
            ENDDO !rIdx
          CASE DEFAULT
            localError="The fourth tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(4),"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(TENSOR_COVARIANT_INDEX)
          SELECT CASE(tensorIndexTypes(4))
          CASE(TENSOR_CONTRAVARIANT_INDEX)
            !(Contravariant,Covariant,Covariant,Contravariant) tensor
            DO rIdx=1,n
              DO sIdx=1,n
                DO tIdx=1,n
                  DO uIdx=1,n
                    sum=0.0_DP
                    DO aIdx=1,n
                      DO bIdx=1,n
                        DO cIdx=1,n
                          DO dIdx=1,n
                            sum=sum+dXidNu(rIdx,aIdx)* &
                              & dNudXi(bIdx,sIdx)* &
                              & dNudXi(cIdx,tIdx)* &
                              & dXidNu(uIdx,dIdx)* &
                              & A(aIdx,bIdx,cIdx,dIdx)
                          ENDDO !dIdx
                        ENDDO !cIdx
                      ENDDO !bIdx
                    ENDDO !aIdx
                    B(rIdx,sIdx,tIdx,uIdx)=sum
                  ENDDO !uIdx
                ENDDO !tIdx
              ENDDO !sIdx
            ENDDO !rIdx
          CASE(TENSOR_COVARIANT_INDEX)
            !(Contravariant,Covariant,Covariant,Covariant) tensor
            DO rIdx=1,n
              DO sIdx=1,n
                DO tIdx=1,n
                  DO uIdx=1,n
                    sum=0.0_DP
                    DO aIdx=1,n
                      DO bIdx=1,n
                        DO cIdx=1,n
                          DO dIdx=1,n
                            sum=sum+dXidNu(rIdx,aIdx)* &
                              & dNudXi(bIdx,sIdx)* &
                              & dNudXi(cIdx,tIdx)* &
                              & dNudXi(dIdx,uIdx)* &
                              & A(aIdx,bIdx,cIdx,dIdx)
                          ENDDO !dIdx
                        ENDDO !cIdx
                      ENDDO !bIdx
                    ENDDO !aIdx
                    B(rIdx,sIdx,tIdx,uIdx)=sum
                   ENDDO !uIdx
                ENDDO !tIdx
              ENDDO !sIdx
            ENDDO !rIdx
          CASE DEFAULT
             localError="The fourth tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(4),"*",err,error))//" is invalid."
             CALL FlagError(localError,err,error,*999)
           END SELECT
         CASE DEFAULT
           localError="The third tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(3),"*",err,error))//" is invalid."
           CALL FlagError(localError,err,error,*999)
         END SELECT
       CASE DEFAULT
         localError="The second tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(2),"*",err,error))//" is invalid."
         CALL FlagError(localError,err,error,*999)
       END SELECT
     CASE(TENSOR_COVARIANT_INDEX)
       SELECT CASE(tensorIndexTypes(2))
       CASE(TENSOR_CONTRAVARIANT_INDEX)
         SELECT CASE(tensorIndexTypes(3))
         CASE(TENSOR_CONTRAVARIANT_INDEX)
           SELECT CASE(tensorIndexTypes(4))
           CASE(TENSOR_CONTRAVARIANT_INDEX)
             !(Covariant,Contravariant,Contravariant,Contravariant) tensor i.e.,
             DO rIdx=1,n
               DO sIdx=1,n
                 DO tIdx=1,n
                   DO uIdx=1,n
                     sum=0.0_DP
                     DO aIdx=1,n
                       DO bIdx=1,n
                         DO cIdx=1,n
                           DO dIdx=1,n
                             sum=sum+dNudXi(aIdx,rIdx)* &
                               & dXidNu(sIdx,bIdx)* &
                               & dXidNu(tIdx,cIdx)* &
                               & dXidNu(uIdx,dIdx)* &
                               & A(aIdx,bIdx,cIdx,dIdx)
                           ENDDO !dIdx
                         ENDDO !cIdx
                       ENDDO !bIdx
                     ENDDO !aIdx
                     B(rIdx,sIdx,tIdx,uIdx)=sum
                   ENDDO !uIdx
                 ENDDO !tIdx
               ENDDO !sIdx
             ENDDO !rIdx
           CASE(TENSOR_COVARIANT_INDEX)
             !(Covariant,Contravariant,Contravariant,Covariant) tensor i.e.,
            DO rIdx=1,n
               DO sIdx=1,n
                 DO tIdx=1,n
                   DO uIdx=1,n
                     sum=0.0_DP
                     DO aIdx=1,n
                       DO bIdx=1,n
                         DO cIdx=1,n
                           DO dIdx=1,n
                             sum=sum+dNudXi(aIdx,rIdx)* &
                               & dXidNu(sIdx,bIdx)* &
                               & dXidNu(tIdx,cIdx)* &
                               & dNudXi(dIdx,uIdx)* &
                               & A(aIdx,bIdx,cIdx,dIdx)
                           ENDDO !dIdx
                         ENDDO !cIdx
                       ENDDO !bIdx
                     ENDDO !aIdx
                     B(rIdx,sIdx,tIdx,uIdx)=sum
                   ENDDO !uIdx
                 ENDDO !tIdx
               ENDDO !sIdx
             ENDDO !rIdx
         CASE DEFAULT
             localError="The fourth tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(4),"*",err,error))//" is invalid."
             CALL FlagError(localError,err,error,*999)
           END SELECT
         CASE(TENSOR_COVARIANT_INDEX)
           SELECT CASE(tensorIndexTypes(4))
           CASE(TENSOR_CONTRAVARIANT_INDEX)
             !(Covariant,Contravariant,Covariant,Contravariant) tensor
             DO rIdx=1,n
               DO sIdx=1,n
                 DO tIdx=1,n
                   DO uIdx=1,n
                     sum=0.0_DP
                     DO aIdx=1,n
                       DO bIdx=1,n
                         DO cIdx=1,n
                           DO dIdx=1,n
                             sum=sum+dNudXi(aIdx,rIdx)* &
                               & dXidNu(sIdx,bIdx)* &
                               & dNudXi(cIdx,tIdx)* &
                               & dXidNu(uIdx,dIdx)* &
                               & A(aIdx,bIdx,cIdx,dIdx)
                           ENDDO !dIdx
                         ENDDO !cIdx
                       ENDDO !bIdx
                     ENDDO !aIdx
                     B(rIdx,sIdx,tIdx,uIdx)=sum
                   ENDDO !uIdx
                 ENDDO !tIdx
               ENDDO !sIdx
             ENDDO !rIdx
           CASE(TENSOR_COVARIANT_INDEX)
             !(Covariant,Contravariant,Covariant,Covariant) tensor
             DO rIdx=1,n
               DO sIdx=1,n
                 DO tIdx=1,n
                   DO uIdx=1,n
                     sum=0.0_DP
                      DO aIdx=1,n
                       DO bIdx=1,n
                         DO cIdx=1,n
                           DO dIdx=1,n
                             sum=sum+dNudXi(aIdx,rIdx)* &
                               & dXidNu(sIdx,bIdx)* &
                               & dNudXi(cIdx,tIdx)* &
                               & dNudXi(dIdx,uIdx)* &
                               & A(aIdx,bIdx,cIdx,dIdx)
                           ENDDO !dIdx
                         ENDDO !cIdx
                       ENDDO !bIdx
                     ENDDO !aIdx
                     B(rIdx,sIdx,tIdx,uIdx)=sum
                   ENDDO !uIdx
                 ENDDO !tIdx
               ENDDO !sIdx
             ENDDO !rIdx
          CASE DEFAULT
             localError="The fourth tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(4),"*",err,error))//" is invalid."
             CALL FlagError(localError,err,error,*999)
           END SELECT
         CASE DEFAULT
           localError="The third tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(3),"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(TENSOR_COVARIANT_INDEX)
        SELECT CASE(tensorIndexTypes(3))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          SELECT CASE(tensorIndexTypes(4))
          CASE(TENSOR_CONTRAVARIANT_INDEX)
            !(Covariant,Covariant,Contravariant,Contravariant) tensor
            DO rIdx=1,n
              DO sIdx=1,n
                DO tIdx=1,n
                  DO uIdx=1,n
                    sum=0.0_DP
                    DO aIdx=1,n
                      DO bIdx=1,n
                        DO cIdx=1,n
                          DO dIdx=1,n
                            sum=sum+dNudXi(aIdx,rIdx)* &
                              & dNudXi(bIdx,sIdx)* &
                              & dXidNu(tIdx,cIdx)* &
                              & dXidNu(uIdx,dIdx)* &
                              & A(aIdx,bIdx,cIdx,dIdx)
                          ENDDO !dIdx
                        ENDDO !cIdx
                      ENDDO !bIdx
                   ENDDO !aIdx
                   B(rIdx,sIdx,tIdx,uIdx)=sum
                  ENDDO !uIdx
                ENDDO !tIdx
              ENDDO !sIdx
            ENDDO !rIdx
          CASE(TENSOR_COVARIANT_INDEX)
            !(Covariant,Covariant,Contravariant,Covariant) tensor
            DO rIdx=1,n
              DO sIdx=1,n
                DO tIdx=1,n
                  DO uIdx=1,n
                    sum=0.0_DP
                    DO aIdx=1,n
                      DO bIdx=1,n
                        DO cIdx=1,n
                          DO dIdx=1,n
                            sum=sum+dNudXi(aIdx,rIdx)* &
                              & dNudXi(bIdx,sIdx)* &
                              & dXidNu(tIdx,cIdx)* &
                              & dNudXi(dIdx,uIdx)* &
                              & A(aIdx,bIdx,cIdx,dIdx)
                          ENDDO !dIdx
                        ENDDO !cIdx
                      ENDDO !bIdx
                    ENDDO !aIdx
                    B(rIdx,sIdx,tIdx,uIdx)=sum
                 ENDDO !uIdx
                ENDDO !tIdx
              ENDDO !sIdx
            ENDDO !rIdx
         CASE DEFAULT
            localError="The fourth tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(4),"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(TENSOR_COVARIANT_INDEX)
          SELECT CASE(tensorIndexTypes(4))
          CASE(TENSOR_CONTRAVARIANT_INDEX)
            !(Covariant,Covariant,Covariant,Contravariant) tensor
            DO rIdx=1,n
              DO sIdx=1,n
                DO tIdx=1,n
                  DO uIdx=1,n
                    sum=0.0_DP
                    DO aIdx=1,n
                      DO bIdx=1,n
                        DO cIdx=1,n
                          DO dIdx=1,n
                            sum=sum+dNudXi(aIdx,rIdx)* &
                              & dNudXi(bIdx,sIdx)* &
                              & dNudXi(cIdx,tIdx)* &
                              & dXidNu(uIdx,dIdx)* &
                              & A(aIdx,bIdx,cIdx,dIdx)
                          ENDDO !dIdx
                        ENDDO !cIdx
                      ENDDO !bIdx
                    ENDDO !aIdx
                  ENDDO !uIdx
                ENDDO !tIdx
              ENDDO !sIdx
            ENDDO !rIdx
            B(rIdx,sIdx,tIdx,uIdx)=sum
          CASE(TENSOR_COVARIANT_INDEX)
            !(Covariant,Covariant,Covariant,Covariant) tensor
            DO rIdx=1,n
              DO sIdx=1,n
                DO tIdx=1,n
                  DO uIdx=1,n
                    sum=0.0_DP
                    DO aIdx=1,n
                      DO bIdx=1,n
                        DO cIdx=1,n
                          DO dIdx=1,n
                            sum=sum+dNudXi(aIdx,rIdx)* &
                              & dNudXi(bIdx,sIdx)* &
                              & dNudXi(cIdx,tIdx)* &
                              & dNudXi(dIdx,uIdx)* &
                              & A(aIdx,bIdx,cIdx,dIdx)
                          ENDDO !dIdx
                        ENDDO !cIdx
                      ENDDO !bIdx
                   ENDDO !aIdx
                   B(rIdx,sIdx,tIdx,uIdx)=sum
                   ENDDO !uIdx
                ENDDO !tIdx
              ENDDO !sIdx
            ENDDO !rIdx
          CASE DEFAULT
            localError="The fourth tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(4),"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The third tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(3),"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The second tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(2),"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The first tensor index type of "//TRIM(NumberToVString(tensorIndexTypes(1),"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("TensorTransformTensorFourDP")
    RETURN
999 ERRORS("TensorTransformTensorFourDP",err,error)
    EXITS("TensorTransformTensorFourDP")
    RETURN 1
    
  END SUBROUTINE TensorTransformTensorFourDP

  !
  !================================================================================================================================
  !

  !>Returns the trace of an integer matrix A.
  SUBROUTINE TraceIntg(A,traceA,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: A(:,:) !<The matrix to find the trace of
    INTEGER(INTG), INTENT(OUT) :: traceA !<On exit, the trace of the A matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i
    
    ENTERS("TraceIntg",err,error,*999)

    traceA=0_INTG

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix is not square.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1))
    CASE(1)
      traceA = A(1,1)
    CASE(2)
      traceA = A(1,1)+A(2,2)
    CASE(3)
      traceA = A(1,1)+A(2,2)+A(3,3)
    CASE DEFAULT
      DO i=1,SIZE(A,1)
        traceA = traceA + A(i,i)
      ENDDO !i
    END SELECT

    EXITS("TraceIntg")
    RETURN
999 ERRORSEXITS("TraceIntg",err,error)
    RETURN 1
    
  END SUBROUTINE TraceIntg
  
  !
  !================================================================================================================================
  !

  !>Returns the trace of a single precision matrix A.
  SUBROUTINE TraceSP(A,traceA,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The matrix to find the trace of
    REAL(SP), INTENT(OUT) :: traceA !<On exit, the trace of the A matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i
    
    ENTERS("TraceSP",err,error,*999)

    traceA = 0.0_SP

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix is not square.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1))
    CASE(1)
      traceA = A(1,1)
    CASE(2)
      traceA = A(1,1)+A(2,2)
    CASE(3)
      traceA = A(1,1)+A(2,2)+A(3,3)
    CASE DEFAULT
      DO i=1,SIZE(A,1)
        traceA = traceA + A(i,i)
      ENDDO !i
    END SELECT

    EXITS("TraceSP")
    RETURN
999 ERRORSEXITS("TraceSP",err,error)
    RETURN 1
    
  END SUBROUTINE TraceSP
  
  !
  !================================================================================================================================
  !

  !>Returns the trace of a double precision matrix A.
  SUBROUTINE TraceDP(A,traceA,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The matrix to find the trace of
    REAL(DP), INTENT(OUT) :: traceA !<On exit, the trace of the A matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i
    
    ENTERS("TraceDP",err,error,*999)

    traceA = 0.0_DP

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("Matrix is not square.",err,error,*999)
#endif    
    
    SELECT CASE(SIZE(A,1))
    CASE(1)
      traceA = A(1,1)
    CASE(2)
      traceA = A(1,1)+A(2,2)
    CASE(3)
      traceA = A(1,1)+A(2,2)+A(3,3)
    CASE DEFAULT
      DO i=1,SIZE(A,1)
        traceA = traceA + A(i,i)
      ENDDO !i
    END SELECT

    EXITS("TraceDP")
    RETURN
999 ERRORSEXITS("TraceDP",err,error)
    RETURN 1
    
  END SUBROUTINE TraceDP
  
  !
  !================================================================================================================================
  !

  !>Returns the transpose of a single precision matrix A in AT.
  SUBROUTINE TransposeMatrixSP(A,AT,err,error,*)
          
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The matrix to take the transpose of
    REAL(SP), INTENT(OUT) :: AT(:,:) !<On exit, the transpose of the matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j
    
    ENTERS("TransposeMatrixSP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(AT,2).OR.SIZE(A,2)/=SIZE(AT,1)) CALL FlagError("Invalid matrix sizes.",err,error,*999)
#endif    
    
    IF(SIZE(A,1)==SIZE(A,2)) THEN
      !A matrix is square
      SELECT CASE(SIZE(A,1))
      CASE(1)
        AT(1,1)=A(1,1)
      CASE(2)
        AT(1,1)=A(1,1)
        AT(1,2)=A(2,1)
        AT(2,1)=A(1,2)
        AT(2,2)=A(2,2)
      CASE(3)
        AT(1,1)=A(1,1)
        AT(1,2)=A(2,1)
        AT(1,3)=A(3,1)
        AT(2,1)=A(1,2)
        AT(2,2)=A(2,2)
        AT(2,3)=A(3,2)
        AT(3,1)=A(1,3)
        AT(3,2)=A(2,3)
        AT(3,3)=A(3,3)
      CASE DEFAULT
        DO i=1,SIZE(A,1)
          DO j=1,SIZE(A,1)
            AT(j,i)=A(i,j)
          ENDDO !j
        ENDDO !i
      END SELECT
    ELSE
      !A matrix is not square
      DO i=1,SIZE(A,1)
        DO j=1,SIZE(A,2)
          AT(j,i)=A(i,j)
        ENDDO !j
      ENDDO !i
    ENDIF
 
    EXITS("TransposeMatrixSP")
    RETURN
999 ERRORSEXITS("TransposeMatrixSP",err,error)
    RETURN 1
    
  END SUBROUTINE TransposeMatrixSP

  !
  !================================================================================================================================
  !

  !>Returns the transpose of a double precision matrix A in AT.
  SUBROUTINE TransposeMatrixDP(A,AT,err,error,*)
    
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The matrix to take the transpose of
    REAL(DP), INTENT(OUT) :: AT(:,:) !<On exit, the transpose of the matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j
            
    ENTERS("TransposeMatrixDP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(AT,2).OR.SIZE(A,2)/=SIZE(AT,1)) CALL FlagError("Invalid matrix sizes.",err,error,*999)
#endif    
    
    IF(SIZE(A,1)==SIZE(A,2)) THEN
      !A matrix is square
      SELECT CASE(SIZE(A,1))
      CASE(1)
        AT(1,1)=A(1,1)
      CASE(2)
        AT(1,1)=A(1,1)
        AT(1,2)=A(2,1)
        AT(2,1)=A(1,2)
        AT(2,2)=A(2,2)
      CASE(3)
        AT(1,1)=A(1,1)
        AT(1,2)=A(2,1)
        AT(1,3)=A(3,1)
        AT(2,1)=A(1,2)
        AT(2,2)=A(2,2)
        AT(2,3)=A(3,2)
        AT(3,1)=A(1,3)
        AT(3,2)=A(2,3)
        AT(3,3)=A(3,3)
      CASE DEFAULT
        DO i=1,SIZE(A,1)
          DO j=1,SIZE(A,1)
            AT(j,i)=A(i,j)
          ENDDO !j
        ENDDO !i
      END SELECT
    ELSE
      !A matrix is not square
      DO i=1,SIZE(A,1)
        DO j=1,SIZE(A,2)
          AT(j,i)=A(i,j)
        ENDDO !j
      ENDDO !i
    ENDIF
    
    EXITS("TransposeMatrixDP")
    RETURN
999 ERRORSEXITS("TransposeMatrixDP",err,error)
    RETURN 1
    
  END SUBROUTINE TransposeMatrixDP

  !
  !================================================================================================================================
  !

  !>Returns the unimodular part of a single precision matrix/tensor A.
  SUBROUTINE UnimodularMatrixSP(A,unimodA,err,error,*)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The matrix to find the unimodular part of
    REAL(SP), INTENT(OUT) :: unimodA(:,:) !<On exit, the unimodular part of the A matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    REAL(SP) :: detA,factor
    
    ENTERS("UnimodularMatrixSP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("The A matrix is not square.",err,error,*999)
    IF(SIZE(A,1)/=SIZE(unimodA,1).OR.SIZE(A,2)/=SIZE(unimodA,2)) &
      & CALL FlagError("The size of the unimodular matrix is not the same size as the A matrix.",err,error,*999)
#endif    

    SELECT CASE(SIZE(A,1))
    CASE(1)
      unimodA(1,1) = 1.0_SP
    CASE(2)
      CALL Determinant(A,detA,err,error,*999)
      factor = SIGN(SQRT(ABS(detA)),detA)
      unimodA = factor*A
    CASE(3)
      CALL Determinant(A,detA,err,error,*999)
      factor = SIGN(ABS(detA)**(-1.0_SP/3.0_SP),detA)
      unimodA = factor*A      
    CASE DEFAULT
      CALL FlagError("Matrix size not implemented.",err,error,*999)
    END SELECT

    EXITS("UnimodularMatrixSP")
    RETURN
999 ERRORSEXITS("UnimodularMatrixSP",err,error)
    RETURN 1
    
  END SUBROUTINE UnimodularMatrixSP
  
  !
  !================================================================================================================================
  !

  !>Returns the unimodular part of a double precision matrix/tensor A.
  SUBROUTINE UnimodularMatrixDP(A,unimodA,err,error,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The matrix to find the unimodular part of
    REAL(DP), INTENT(OUT) :: unimodA(:,:) !<On exit, the unimodular part of the A matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    REAL(DP) :: detA,factor
    
    ENTERS("UnimodularMatrixDP",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(SIZE(A,1)/=SIZE(A,2)) CALL FlagError("The A matrix is not square.",err,error,*999)
    IF(SIZE(A,1)/=SIZE(unimodA,1).OR.SIZE(A,2)/=SIZE(unimodA,2)) &
      & CALL FlagError("The size of the unimodular matrix is not the same size as the A matrix.",err,error,*999)
#endif    

    SELECT CASE(SIZE(A,1))
    CASE(1)
      unimodA(1,1) = 1.0_DP
    CASE(2)
      CALL Determinant(A,detA,err,error,*999)
      factor = SIGN(SQRT(ABS(detA)),detA)
      unimodA = factor*A
    CASE(3)
      CALL Determinant(A,detA,err,error,*999)
      factor = SIGN(ABS(detA)**(-1.0_DP/3.0_DP),detA)
      unimodA = factor*A      
    CASE DEFAULT
      CALL FlagError("Matrix size not implemented.",err,error,*999)
    END SELECT

    EXITS("UnimodularMatrixDP")
    RETURN
999 ERRORSEXITS("UnimodularMatrixDP",err,error)
    RETURN 1
    
  END SUBROUTINE UnimodularMatrixDP
  
  !
  !================================================================================================================================
  !

  !>Converts a Voigt representation of a single precision second order tensor into a matrix form of the tensor. \todo Think about Voigt scale factors.
  SUBROUTINE VoigtToTensor2SP0(n,voigtTensorIndexType,AV,A,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: n !<The number of dimensions
    INTEGER(INTG), INTENT(IN) :: voigtTensorIndexType !<The index type of the Voigt second order tensor \see Constants_TensorIndexTypes
    REAL(SP), INTENT(IN) :: AV(:) !<AV(voigtIdx). The second order tensor in Voigt form to convert
    REAL(SP), INTENT(OUT) :: A(:,:) !<A(iIdx,jIdx). On exit, the second order tensor converted from Voigt form.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("VoigtToTensor2SP0",err,error,*999)

    CALL VoigtToTensor2SP1(n,[voigtTensorIndexType],AV,A,err,error,*999)

    EXITS("VoigtToTensor2SP0")
    RETURN
999 ERRORSEXITS("VoigtToTensor2SP0",err,error)
    RETURN 1
    
  END SUBROUTINE VoigtToTensor2SP0

  !
  !================================================================================================================================
  !

  !>Converts a Voigt representation of a single precision second order tensor into a matrix form of the tensor.
  SUBROUTINE VoigtToTensor2SP1(n,voigtTensorIndexTypes,AV,A,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: n !<The number of dimensions
    INTEGER(INTG), INTENT(IN) :: voigtTensorIndexTypes(:) !<The index type of the Voigt second order tensor \see Constants_TensorIndexTypes
    REAL(SP), INTENT(IN) :: AV(:) !<AV(voigtIdx). The second order tensor in Voigt form to convert
    REAL(SP), INTENT(OUT) :: A(:,:) !<A(iIdx,jIdx). On exit, the second order tensor converted from Voigt form.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("VoigtToTensor2SP1",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(n<1.OR.n>3) THEN
      localError="The specified number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))// &
        & " is invalid. The number of dimensions should be >= 1 and <= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(voigtTensorIndexTypes,1)<1) THEN
      localError="The size of the Voigt tensor index types array of "// &
        & TRIM(NumberToVString(SIZE(voigtTensorIndexTypes,1),"*",err,error))// &
        & " is too small. The size should be >= 1 for a second order Voigt array."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(AV,1)<NUMBER_OF_VOIGT(n)) THEN
      localError="The size of the Voigt array of "//TRIM(NumberToVString(SIZE(AV,1),"*",err,error))// &
        & " is too small for a second order Voigt array with "//TRIM(NumberToVString(n,"*",err,error))// &
        & " dimensions. The size of the Voigt array should be >= "//TRIM(NumberToVString(NUMBER_OF_VOIGT(n),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,1)<n) THEN
      localError="The size of the first index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,1),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,2)<n) THEN
      localError="The size of the second index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,2),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    SELECT CASE(n)
    CASE(1)
      A(1,1)=AV(TENSOR_TO_VOIGT1(1,1))
    CASE(2)
      A(1,1)=AV(TENSOR_TO_VOIGT2(1,1))
      A(2,2)=AV(TENSOR_TO_VOIGT2(2,2))
      SELECT CASE(voigtTensorIndexTypes(1))
      CASE(TENSOR_CONTRAVARIANT_INDEX)
        A(2,1)=AV(TENSOR_TO_VOIGT2(2,1))
        A(1,2)=AV(TENSOR_TO_VOIGT2(1,2))
      CASE(TENSOR_COVARIANT_INDEX)
        A(2,1)=0.5_SP*AV(TENSOR_TO_VOIGT2(2,1))
        A(1,2)=0.5_SP*AV(TENSOR_TO_VOIGT2(1,2))
      CASE DEFAULT
        localError="The Voigt tensor index type of "//TRIM(NumberToVString(voigtTensorIndexTypes(1),"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(3)
      A(1,1)=AV(TENSOR_TO_VOIGT3(1,1))
      A(2,2)=AV(TENSOR_TO_VOIGT3(2,2))
      A(3,3)=AV(TENSOR_TO_VOIGT3(3,3))
      SELECT CASE(voigtTensorIndexTypes(1))
      CASE(TENSOR_CONTRAVARIANT_INDEX)
        A(2,1)=AV(TENSOR_TO_VOIGT3(2,1))
        A(3,1)=AV(TENSOR_TO_VOIGT3(3,1))
        A(1,2)=AV(TENSOR_TO_VOIGT3(1,2))
        A(3,2)=AV(TENSOR_TO_VOIGT3(3,2))
        A(1,3)=AV(TENSOR_TO_VOIGT3(1,3))
        A(2,3)=AV(TENSOR_TO_VOIGT3(2,3))
      CASE(TENSOR_COVARIANT_INDEX)
        A(2,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,1))
        A(3,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,1))
        A(1,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,2))
        A(3,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,2))
        A(1,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,3))
        A(2,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,3))
      CASE DEFAULT
        localError="The Voigt tensor index type of "//TRIM(NumberToVString(voigtTensorIndexTypes(1),"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT    

    EXITS("VoigtToTensor2SP1")
    RETURN
999 ERRORSEXITS("VoigtToTensor2SP1",err,error)
    RETURN 1
    
  END SUBROUTINE VoigtToTensor2SP1

  !
  !================================================================================================================================
  !

  !>Converts a Voigt representation of a double precision second order tensor into a matrix form of the tensor. \todo Think about Voigt scale factors.
  SUBROUTINE VoigtToTensor2DP0(n,voigtTensorIndexType,AV,A,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: n !<The number of dimensions
    INTEGER(INTG), INTENT(IN) :: voigtTensorIndexType !<The index type of the Voigt second order tensor \see Constants_TensorIndexTypes
    REAL(DP), INTENT(IN) :: AV(:) !<AV(voigtIdx). The second order tensor in Voigt form to convert
    REAL(DP), INTENT(OUT) :: A(:,:) !<A(iIdx,jIdx). On exit, the second order tensor converted from Voigt form.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("VoigtToTensor2DP0",err,error,*999)

    CALL VoigtToTensor2DP1(n,[voigtTensorIndexType],AV,A,err,error,*999)

    EXITS("VoigtToTensor2DP0")
    RETURN
999 ERRORSEXITS("VoigtToTensor2DP0",err,error)
    RETURN 1
    
  END SUBROUTINE VoigtToTensor2DP0

  !
  !================================================================================================================================
  !

  !>Converts a Voigt representation of a double precision second order tensor into a matrix form of the tensor.
  SUBROUTINE VoigtToTensor2DP1(n,voigtTensorIndexTypes,AV,A,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: n !<The number of dimensions
    INTEGER(INTG), INTENT(IN) :: voigtTensorIndexTypes(:) !<The index type of the Voigt second order tensor \see Constants_TensorIndexTypes
    REAL(DP), INTENT(IN) :: AV(:) !<AV(voigtIdx). The second order tensor in Voigt form to convert
    REAL(DP), INTENT(OUT) :: A(:,:) !<A(iIdx,jIdx). On exit, the second order tensor converted from Voigt form.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("VoigtToTensor2DP1",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(n<1.OR.n>3) THEN
      localError="The specified number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))// &
        & " is invalid. The number of dimensions should be >= 1 and <= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(voigtTensorIndexTypes,1)<1) THEN
      localError="The size of the Voigt tensor index types array of "// &
        & TRIM(NumberToVString(SIZE(voigtTensorIndexTypes,1),"*",err,error))// &
        & " is too small. The size should be >= 1 for a second order Voigt array."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(AV,1)<NUMBER_OF_VOIGT(n)) THEN
      localError="The size of the Voigt array of "//TRIM(NumberToVString(SIZE(AV,1),"*",err,error))// &
        & " is too small for a second order Voigt array with "//TRIM(NumberToVString(n,"*",err,error))// &
        & " dimensions. The size of the Voigt array should be >= "//TRIM(NumberToVString(NUMBER_OF_VOIGT(n),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,1)<n) THEN
      localError="The size of the first index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,1),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,2)<n) THEN
      localError="The size of the second index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,2),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    SELECT CASE(n)
    CASE(1)
      A(1,1)=AV(TENSOR_TO_VOIGT1(1,1))
    CASE(2)
      A(1,1)=AV(TENSOR_TO_VOIGT2(1,1))
      A(2,2)=AV(TENSOR_TO_VOIGT2(2,2))
      SELECT CASE(voigtTensorIndexTypes(1))
      CASE(TENSOR_CONTRAVARIANT_INDEX)
        A(2,1)=AV(TENSOR_TO_VOIGT2(2,1))
        A(1,2)=AV(TENSOR_TO_VOIGT2(1,2))
      CASE(TENSOR_COVARIANT_INDEX)
        A(2,1)=0.5_DP*AV(TENSOR_TO_VOIGT2(2,1))
        A(1,2)=0.5_DP*AV(TENSOR_TO_VOIGT2(1,2))
      CASE DEFAULT
        localError="The Voigt tensor index type of "//TRIM(NumberToVString(voigtTensorIndexTypes(1),"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(3)
      A(1,1)=AV(TENSOR_TO_VOIGT3(1,1))
      A(2,2)=AV(TENSOR_TO_VOIGT3(2,2))
      A(3,3)=AV(TENSOR_TO_VOIGT3(3,3))
      SELECT CASE(voigtTensorIndexTypes(1))
      CASE(TENSOR_CONTRAVARIANT_INDEX)
        A(2,1)=AV(TENSOR_TO_VOIGT3(2,1))
        A(3,1)=AV(TENSOR_TO_VOIGT3(3,1))
        A(1,2)=AV(TENSOR_TO_VOIGT3(1,2))
        A(3,2)=AV(TENSOR_TO_VOIGT3(3,2))
        A(1,3)=AV(TENSOR_TO_VOIGT3(1,3))
        A(2,3)=AV(TENSOR_TO_VOIGT3(2,3))
      CASE(TENSOR_COVARIANT_INDEX)
        A(2,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,1))
        A(3,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,1))
        A(1,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,2))
        A(3,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,2))
        A(1,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,3))
        A(2,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,3))
      CASE DEFAULT
        localError="The Voigt tensor index type of "//TRIM(NumberToVString(voigtTensorIndexTypes(1),"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT    

    EXITS("VoigtToTensor2DP1")
    RETURN
999 ERRORSEXITS("VoigtToTensor2DP1",err,error)
    RETURN 1
    
  END SUBROUTINE VoigtToTensor2DP1

  !
  !================================================================================================================================
  !

  !>Converts a Voigt representation of a single precision fourth order tensor into a matrix form of the tensor.
  SUBROUTINE VoigtToTensor4SP(n,voigtTensorIndexTypes,AV,A,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: n !<The number of dimensions
    INTEGER(INTG), INTENT(IN) :: voigtTensorIndexTypes(:) !<The index types of the Voigt fourth order tensor \see Constants_TensorIndexTypes
    REAL(SP), INTENT(IN) :: AV(:,:) !<AV(voigtIdx,voigtIdx). The fourth order tensor in Voigt form to convert
    REAL(SP), INTENT(OUT) :: A(:,:,:,:) !<A(iIdx,jIdx,kIdx,lIdx). On exit, the fourth order tensor converted from Voigt form.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("VoigtToTensor4SP",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(n<1.OR.n>3) THEN
      localError="The specified number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))// &
        & " is invalid. The number of dimensions should be >= 1 and <= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(voigtTensorIndexTypes,1)<2) THEN
      localError="The size of the Voigt tensor index types array of "// &
        & TRIM(NumberToVString(SIZE(voigtTensorIndexTypes,1),"*",err,error))// &
        & " is too small. The size should be >= 2 for a fourth order Voigt array."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(AV,1)<NUMBER_OF_VOIGT(n)) THEN
      localError="The size of the first index of the Voigt array of "//TRIM(NumberToVString(SIZE(AV,1),"*",err,error))// &
        & " is too small for a fourth order Voigt array with "//TRIM(NumberToVString(n,"*",err,error))// &
        & " dimensions. The size of the first index of the Voigt array should be >= "// &
        & TRIM(NumberToVString(NUMBER_OF_VOIGT(n),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(AV,2)<NUMBER_OF_VOIGT(n)) THEN
      localError="The size of the second index of the Voigt array of "//TRIM(NumberToVString(SIZE(AV,2),"*",err,error))// &
        & " is too small for a fourth order Voigt array with "//TRIM(NumberToVString(n,"*",err,error))// &
        & " dimensions. The size of the second index of the Voigt array should be >= "// &
        & TRIM(NumberToVString(NUMBER_OF_VOIGT(n),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,1)<n) THEN
      localError="The size of the first index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,1),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,2)<n) THEN
      localError="The size of the second index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,2),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,3)<n) THEN
      localError="The size of the third index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,3),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,4)<n) THEN
      localError="The size of the fourth index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,4),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    SELECT CASE(n)
    CASE(1)
      A(1,1,1,1)=AV(TENSOR_TO_VOIGT1(1,1),TENSOR_TO_VOIGT1(1,1))
    CASE(2)
      A(1,1,1,1)=AV(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(1,1))
      A(2,2,1,1)=AV(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(1,1))
      A(1,1,2,2)=AV(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(2,2))
      A(2,2,2,2)=AV(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(2,2))     
      SELECT CASE(voigtTensorIndexTypes(1))
      CASE(TENSOR_CONTRAVARIANT_INDEX)
        SELECT CASE(voigtTensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Voigt array is (contravariant,contravariant)
          A(2,1,1,1)=AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(1,1))
          A(1,2,1,1)=AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(1,1))
          A(1,1,2,1)=AV(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(2,1))
          A(2,1,2,1)=AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(2,1))
          A(1,2,2,1)=AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(2,1))
          A(2,2,2,1)=AV(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(2,1))
          A(1,1,1,2)=AV(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(1,2))
          A(2,1,1,2)=AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(1,2))
          A(1,2,1,2)=AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(1,2))
          A(2,2,1,2)=AV(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(1,2))
          A(2,1,2,2)=AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(2,2))
          A(1,2,2,2)=AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(2,2))          
        CASE(TENSOR_COVARIANT_INDEX)
          !Voigt array is (contravariant,covariant)
          A(2,1,1,1)=AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(1,1))
          A(1,2,1,1)=AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(1,1))
          A(1,1,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(2,1))
          A(2,1,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(2,1))
          A(1,2,2,1)=AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(2,1))
          A(2,2,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(2,1))
          A(1,1,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(1,2))
          A(2,1,1,2)=AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(1,2))
          A(1,2,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(1,2))
          A(2,2,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(1,2))
          A(2,1,2,2)=AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(2,2))
          A(1,2,2,2)=AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(2,2))          
        CASE DEFAULT
          localError="The second Voigt tensor index type of "//TRIM(NumberToVString(voigtTensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(TENSOR_COVARIANT_INDEX)
        SELECT CASE(voigtTensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Voigt array is (covariant,contravariant)
          A(2,1,1,1)=0.5_SP*AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(1,1))
          A(1,2,1,1)=0.5_SP*AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(1,1))
          A(1,1,2,1)=AV(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(2,1))
          A(2,1,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(2,1))
          A(1,2,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(2,1))
          A(2,2,2,1)=AV(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(2,1))
          A(1,1,1,2)=AV(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(1,2))
          A(2,1,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(1,2))
          A(1,2,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(1,2))
          A(2,2,1,2)=AV(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(1,2))
          A(2,1,2,2)=0.5_SP*AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(2,2))
          A(1,2,2,2)=0.5_SP*AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(2,2))          
        CASE(TENSOR_COVARIANT_INDEX)
          !Voigt array is (covariant,covariant)
          A(2,1,1,1)=0.5_SP*AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(1,1))
          A(1,2,1,1)=0.5_SP*AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(1,1))
          A(1,1,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(2,1))
          A(2,1,2,1)=0.25_SP*AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(2,1))
          A(1,2,2,1)=0.25_SP*AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(2,1))
          A(2,2,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(2,1))
          A(1,1,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(1,2))
          A(2,1,1,2)=0.25_SP*AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(1,2))
          A(1,2,1,2)=0.25_SP*AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(1,2))
          A(2,2,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(1,2))
          A(2,1,2,2)=0.5_SP*AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(2,2))
          A(1,2,2,2)=0.5_SP*AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(2,2))          
        CASE DEFAULT
          localError="The second Voigt tensor index type of "//TRIM(NumberToVString(voigtTensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The first Voigt tensor index type of "//TRIM(NumberToVString(voigtTensorIndexTypes(1),"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(3)
      A(1,1,1,1)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(1,1))
      A(2,2,1,1)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(1,1))
      A(3,3,1,1)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(1,1))
      A(1,1,2,2)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(2,2))
      A(2,2,2,2)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(2,2))
      A(3,3,2,2)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(2,2))
      A(1,1,3,3)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(3,3))
      A(2,2,3,3)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(3,3))
      A(3,3,3,3)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(3,3))
      SELECT CASE(voigtTensorIndexTypes(1))
      CASE(TENSOR_CONTRAVARIANT_INDEX)
        SELECT CASE(voigtTensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Voigt array is (contravariant,contravariant)
          A(2,1,1,1)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,1))
          A(3,1,1,1)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,1))
          A(1,2,1,1)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,1))
          A(3,2,1,1)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,1))
          A(1,3,1,1)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,1))
          A(2,3,1,1)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,1))
          A(1,1,2,1)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(2,1))
          A(2,1,2,1)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,1))
          A(3,1,2,1)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,1))
          A(1,2,2,1)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,1))
          A(2,2,2,1)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(2,1))
          A(3,2,2,1)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,1))
          A(1,3,2,1)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,1))
          A(2,3,2,1)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,1))
          A(3,3,2,1)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(2,1))
          A(1,1,3,1)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(3,1))
          A(2,1,3,1)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,1))
          A(3,1,3,1)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,1))
          A(1,2,3,1)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,1))
          A(2,2,3,1)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(3,1))
          A(3,2,3,1)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,1))
          A(1,3,3,1)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,1))
          A(2,3,3,1)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,1))
          A(3,3,3,1)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(3,1))
          A(1,1,1,2)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(1,2))
          A(2,1,1,2)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,2))
          A(3,1,1,2)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,2))
          A(1,2,1,2)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,2))
          A(2,2,1,2)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(1,2))
          A(3,2,1,2)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,2))
          A(1,3,1,2)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,2))
          A(2,3,1,2)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,2))
          A(3,3,1,2)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(1,2))
          A(2,1,2,2)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,2))
          A(3,1,2,2)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,2))
          A(1,2,2,2)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,2))
          A(3,2,2,2)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,2))
          A(1,3,2,2)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,2))
          A(2,3,2,2)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,2))
          A(1,1,3,2)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(3,2))
          A(2,1,3,2)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,2))
          A(3,1,3,2)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,2))
          A(1,2,3,2)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,2))
          A(2,2,3,2)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(3,2))
          A(3,2,3,2)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,2))
          A(1,3,3,2)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,2))
          A(2,3,3,2)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,2))
          A(3,3,3,2)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(3,2))
          A(1,1,1,3)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(1,3))
          A(2,1,1,3)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,3))
          A(3,1,1,3)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,3))
          A(1,2,1,3)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,3))
          A(2,2,1,3)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(1,3))
          A(3,2,1,3)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,3))
          A(1,3,1,3)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,3))
          A(2,3,1,3)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,3))
          A(3,3,1,3)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(1,3))
          A(1,1,2,3)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(2,3))
          A(2,1,2,3)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,3))
          A(3,1,2,3)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,3))
          A(1,2,2,3)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,3))
          A(2,2,2,3)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(2,3))
          A(3,2,2,3)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,3))
          A(1,3,2,3)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,3))
          A(2,3,2,3)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,3))
          A(3,3,2,3)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(2,3))
          A(2,1,3,3)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,3))
          A(3,1,3,3)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,3))
          A(1,2,3,3)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,3))
          A(3,2,3,3)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,3))
          A(1,3,3,3)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,3))
          A(2,3,3,3)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,3))
        CASE(TENSOR_COVARIANT_INDEX)
          !Voigt array is (contravariant,covariant)
          A(2,1,1,1)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,1))
          A(3,1,1,1)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,1))
          A(1,2,1,1)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,1))
          A(3,2,1,1)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,1))
          A(1,3,1,1)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,1))
          A(2,3,1,1)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,1))
          A(1,1,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(2,1))
          A(2,1,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,1))
          A(3,1,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,1))
          A(1,2,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,1))
          A(2,2,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(2,1))
          A(3,2,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,1))
          A(1,3,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,1))
          A(2,3,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,1))
          A(3,3,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(2,1))
          A(1,1,3,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(3,1))
          A(2,1,3,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,1))
          A(3,1,3,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,1))
          A(1,2,3,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,1))
          A(2,2,3,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(3,1))
          A(3,2,3,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,1))
          A(1,3,3,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,1))
          A(2,3,3,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,1))
          A(3,3,3,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(3,1))
          A(1,1,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(1,2))
          A(2,1,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,2))
          A(3,1,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,2))
          A(1,2,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,2))
          A(2,2,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(1,2))
          A(3,2,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,2))
          A(1,3,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,2))
          A(2,3,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,2))
          A(3,3,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(1,2))
          A(2,1,2,2)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,2))
          A(3,1,2,2)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,2))
          A(1,2,2,2)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,2))
          A(3,2,2,2)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,2))
          A(1,3,2,2)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,2))
          A(2,3,2,2)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,2))
          A(1,1,3,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(3,2))
          A(2,1,3,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,2))
          A(3,1,3,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,2))
          A(1,2,3,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,2))
          A(2,2,3,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(3,2))
          A(3,2,3,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,2))
          A(1,3,3,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,2))
          A(2,3,3,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,2))
          A(3,3,3,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(3,2))
          A(1,1,1,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(1,3))
          A(2,1,1,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,3))
          A(3,1,1,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,3))
          A(1,2,1,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,3))
          A(2,2,1,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(1,3))
          A(3,2,1,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,3))
          A(1,3,1,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,3))
          A(2,3,1,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,3))
          A(3,3,1,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(1,3))
          A(1,1,2,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(2,3))
          A(2,1,2,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,3))
          A(3,1,2,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,3))
          A(1,2,2,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,3))
          A(2,2,2,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(2,3))
          A(3,2,2,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,3))
          A(1,3,2,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,3))
          A(2,3,2,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,3))
          A(3,3,2,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(2,3))
          A(2,1,3,3)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,3))
          A(3,1,3,3)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,3))
          A(1,2,3,3)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,3))
          A(3,2,3,3)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,3))
          A(1,3,3,3)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,3))
          A(2,3,3,3)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,3))
        CASE DEFAULT
          localError="The second Voigt tensor index type of "//TRIM(NumberToVString(voigtTensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(TENSOR_COVARIANT_INDEX)
        SELECT CASE(voigtTensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Voigt array is (covariant,contravariant)
          A(2,1,1,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,1))
          A(3,1,1,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,1))
          A(1,2,1,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,1))
          A(3,2,1,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,1))
          A(1,3,1,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,1))
          A(2,3,1,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,1))
          A(1,1,2,1)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(2,1))
          A(2,1,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,1))
          A(3,1,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,1))
          A(1,2,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,1))
          A(2,2,2,1)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(2,1))
          A(3,2,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,1))
          A(1,3,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,1))
          A(2,3,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,1))
          A(3,3,2,1)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(2,1))
          A(1,1,3,1)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(3,1))
          A(2,1,3,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,1))
          A(3,1,3,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,1))
          A(1,2,3,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,1))
          A(2,2,3,1)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(3,1))
          A(3,2,3,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,1))
          A(1,3,3,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,1))
          A(2,3,3,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,1))
          A(3,3,3,1)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(3,1))
          A(1,1,1,2)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(1,2))
          A(2,1,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,2))
          A(3,1,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,2))
          A(1,2,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,2))
          A(2,2,1,2)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(1,2))
          A(3,2,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,2))
          A(1,3,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,2))
          A(2,3,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,2))
          A(3,3,1,2)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(1,2))
          A(2,1,2,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,2))
          A(3,1,2,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,2))
          A(1,2,2,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,2))
          A(3,2,2,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,2))
          A(1,3,2,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,2))
          A(2,3,2,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,2))
          A(1,1,3,2)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(3,2))
          A(2,1,3,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,2))
          A(3,1,3,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,2))
          A(1,2,3,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,2))
          A(2,2,3,2)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(3,2))
          A(3,2,3,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,2))
          A(1,3,3,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,2))
          A(2,3,3,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,2))
          A(3,3,3,2)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(3,2))
          A(1,1,1,3)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(1,3))
          A(2,1,1,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,3))
          A(3,1,1,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,3))
          A(1,2,1,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,3))
          A(2,2,1,3)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(1,3))
          A(3,2,1,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,3))
          A(1,3,1,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,3))
          A(2,3,1,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,3))
          A(3,3,1,3)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(1,3))
          A(1,1,2,3)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(2,3))
          A(2,1,2,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,3))
          A(3,1,2,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,3))
          A(1,2,2,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,3))
          A(2,2,2,3)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(2,3))
          A(3,2,2,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,3))
          A(1,3,2,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,3))
          A(2,3,2,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,3))
          A(3,3,2,3)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(2,3))
          A(2,1,3,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,3))
          A(3,1,3,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,3))
          A(1,2,3,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,3))
          A(3,2,3,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,3))
          A(1,3,3,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,3))
          A(2,3,3,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,3))
        CASE(TENSOR_COVARIANT_INDEX)
          !Voigt array is (covariant,covariant)
          A(2,1,1,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,1))
          A(3,1,1,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,1))
          A(1,2,1,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,1))
          A(3,2,1,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,1))
          A(1,3,1,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,1))
          A(2,3,1,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,1))
          A(1,1,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(2,1))
          A(2,1,2,1)=0.25_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,1))
          A(3,1,2,1)=0.25_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,1))
          A(1,2,2,1)=0.25_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,1))
          A(2,2,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(2,1))
          A(3,2,2,1)=0.25_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,1))
          A(1,3,2,1)=0.25_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,1))
          A(2,3,2,1)=0.25_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,1))
          A(3,3,2,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(2,1))
          A(1,1,3,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(3,1))
          A(2,1,3,1)=0.25_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,1))
          A(3,1,3,1)=0.25_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,1))
          A(1,2,3,1)=0.25_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,1))
          A(2,2,3,1)=0.25_SP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(3,1))
          A(3,2,3,1)=0.25_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,1))
          A(1,3,3,1)=0.25_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,1))
          A(2,3,3,1)=0.25_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,1))
          A(3,3,3,1)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(3,1))
          A(1,1,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(1,2))
          A(2,1,1,2)=0.25_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,2))
          A(3,1,1,2)=0.25_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,2))
          A(1,2,1,2)=0.25_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,2))
          A(2,2,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(1,2))
          A(3,2,1,2)=0.25_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,2))
          A(1,3,1,2)=0.25_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,2))
          A(2,3,1,2)=0.25_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,2))
          A(3,3,1,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(1,2))
          A(2,1,2,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,2))
          A(3,1,2,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,2))
          A(1,2,2,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,2))
          A(3,2,2,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,2))
          A(1,3,2,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,2))
          A(2,3,2,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,2))
          A(1,1,3,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(3,2))
          A(2,1,3,2)=0.25_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,2))
          A(3,1,3,2)=0.25_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,2))
          A(1,2,3,2)=0.25_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,2))
          A(2,2,3,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(3,2))
          A(3,2,3,2)=0.25_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,2))
          A(1,3,3,2)=0.25_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,2))
          A(2,3,3,2)=0.25_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,2))
          A(3,3,3,2)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(3,2))
          A(1,1,1,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(1,3))
          A(2,1,1,3)=0.25_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,3))
          A(3,1,1,3)=0.25_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,3))
          A(1,2,1,3)=0.25_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,3))
          A(2,2,1,3)=0.25_SP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(1,3))
          A(3,2,1,3)=0.25_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,3))
          A(1,3,1,3)=0.25_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,3))
          A(2,3,1,3)=0.25_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,3))
          A(3,3,1,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(1,3))
          A(1,1,2,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(2,3))
          A(2,1,2,3)=0.25_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,3))
          A(3,1,2,3)=0.25_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,3))
          A(1,2,2,3)=0.25_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,3))
          A(2,2,2,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(2,3))
          A(3,2,2,3)=0.25_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,3))
          A(1,3,2,3)=0.25_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,3))
          A(2,3,2,3)=0.25_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,3))
          A(3,3,2,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(2,3))
          A(2,1,3,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,3))
          A(3,1,3,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,3))
          A(1,2,3,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,3))
          A(3,2,3,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,3))
          A(1,3,3,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,3))
          A(2,3,3,3)=0.5_SP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,3))
        CASE DEFAULT
          localError="The second Voigt tensor index type of "//TRIM(NumberToVString(voigtTensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The first Voigt tensor index type of "//TRIM(NumberToVString(voigtTensorIndexTypes(1),"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT    

    EXITS("VoigtToTensor4SP")
    RETURN
999 ERRORSEXITS("VoigtToTensor4SP",err,error)
    RETURN 1
    
  END SUBROUTINE VoigtToTensor4SP

  !
  !================================================================================================================================
  !

  !>Converts a Voigt representation of a double precision fourth order tensor into a matrix form of the tensor.
  SUBROUTINE VoigtToTensor4DP(n,voigtTensorIndexTypes,AV,A,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: n !<The number of dimensions
    INTEGER(INTG), INTENT(IN) :: voigtTensorIndexTypes(:) !<The index types of the Voigt fourth order tensor \see Constants_TensorIndexTypes
    REAL(DP), INTENT(IN) :: AV(:,:) !<AV(voigtIdx,voigtIdx). The fourth order tensor in Voigt form to convert
    REAL(DP), INTENT(OUT) :: A(:,:,:,:) !<A(iIdx,jIdx,kIdx,lIdx). On exit, the fourth order tensor converted from Voigt form.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("VoigtToTensor4DP",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(n<1.OR.n>3) THEN
      localError="The specified number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))// &
        & " is invalid. The number of dimensions should be >= 1 and <= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(voigtTensorIndexTypes,1)<2) THEN
      localError="The size of the Voigt tensor index types array of "// &
        & TRIM(NumberToVString(SIZE(voigtTensorIndexTypes,1),"*",err,error))// &
        & " is too small. The size should be >= 2 for a fourth order Voigt array."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(AV,1)<NUMBER_OF_VOIGT(n)) THEN
      localError="The size of the first index of the Voigt array of "//TRIM(NumberToVString(SIZE(AV,1),"*",err,error))// &
        & " is too small for a fourth order Voigt array with "//TRIM(NumberToVString(n,"*",err,error))// &
        & " dimensions. The size of the first index of the Voigt array should be >= "// &
        & TRIM(NumberToVString(NUMBER_OF_VOIGT(n),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(AV,2)<NUMBER_OF_VOIGT(n)) THEN
      localError="The size of the second index of the Voigt array of "//TRIM(NumberToVString(SIZE(AV,2),"*",err,error))// &
        & " is too small for a fourth order Voigt array with "//TRIM(NumberToVString(n,"*",err,error))// &
        & " dimensions. The size of the second index of the Voigt array should be >= "// &
        & TRIM(NumberToVString(NUMBER_OF_VOIGT(n),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,1)<n) THEN
      localError="The size of the first index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,1),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,2)<n) THEN
      localError="The size of the second index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,2),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,3)<n) THEN
      localError="The size of the third index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,3),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(A,4)<n) THEN
      localError="The size of the fourth index of the tensor matrix of "//TRIM(NumberToVString(SIZE(A,4),"*",err,error))// &
        & " is too small. The size should be >= than the number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    SELECT CASE(n)
    CASE(1)
      A(1,1,1,1)=AV(TENSOR_TO_VOIGT1(1,1),TENSOR_TO_VOIGT1(1,1))
    CASE(2)
      A(1,1,1,1)=AV(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(1,1))
      A(2,2,1,1)=AV(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(1,1))
      A(1,1,2,2)=AV(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(2,2))
      A(2,2,2,2)=AV(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(2,2))     
      SELECT CASE(voigtTensorIndexTypes(1))
      CASE(TENSOR_CONTRAVARIANT_INDEX)
        SELECT CASE(voigtTensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Voigt array is (contravariant,contravariant)
          A(2,1,1,1)=AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(1,1))
          A(1,2,1,1)=AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(1,1))
          A(1,1,2,1)=AV(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(2,1))
          A(2,1,2,1)=AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(2,1))
          A(1,2,2,1)=AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(2,1))
          A(2,2,2,1)=AV(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(2,1))
          A(1,1,1,2)=AV(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(1,2))
          A(2,1,1,2)=AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(1,2))
          A(1,2,1,2)=AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(1,2))
          A(2,2,1,2)=AV(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(1,2))
          A(2,1,2,2)=AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(2,2))
          A(1,2,2,2)=AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(2,2))          
        CASE(TENSOR_COVARIANT_INDEX)
          !Voigt array is (contravariant,covariant)
          A(2,1,1,1)=AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(1,1))
          A(1,2,1,1)=AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(1,1))
          A(1,1,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(2,1))
          A(2,1,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(2,1))
          A(1,2,2,1)=AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(2,1))
          A(2,2,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(2,1))
          A(1,1,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(1,2))
          A(2,1,1,2)=AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(1,2))
          A(1,2,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(1,2))
          A(2,2,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(1,2))
          A(2,1,2,2)=AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(2,2))
          A(1,2,2,2)=AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(2,2))          
        CASE DEFAULT
          localError="The second Voigt tensor index type of "//TRIM(NumberToVString(voigtTensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(TENSOR_COVARIANT_INDEX)
        SELECT CASE(voigtTensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Voigt array is (covariant,contravariant)
          A(2,1,1,1)=0.5_DP*AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(1,1))
          A(1,2,1,1)=0.5_DP*AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(1,1))
          A(1,1,2,1)=AV(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(2,1))
          A(2,1,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(2,1))
          A(1,2,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(2,1))
          A(2,2,2,1)=AV(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(2,1))
          A(1,1,1,2)=AV(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(1,2))
          A(2,1,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(1,2))
          A(1,2,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(1,2))
          A(2,2,1,2)=AV(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(1,2))
          A(2,1,2,2)=0.5_DP*AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(2,2))
          A(1,2,2,2)=0.5_DP*AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(2,2))          
        CASE(TENSOR_COVARIANT_INDEX)
          !Voigt array is (covariant,covariant)
          A(2,1,1,1)=0.5_DP*AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(1,1))
          A(1,2,1,1)=0.5_DP*AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(1,1))
          A(1,1,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(2,1))
          A(2,1,2,1)=0.25_DP*AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(2,1))
          A(1,2,2,1)=0.25_DP*AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(2,1))
          A(2,2,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(2,1))
          A(1,1,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT2(1,1),TENSOR_TO_VOIGT2(1,2))
          A(2,1,1,2)=0.25_DP*AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(1,2))
          A(1,2,1,2)=0.25_DP*AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(1,2))
          A(2,2,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT2(2,2),TENSOR_TO_VOIGT2(1,2))
          A(2,1,2,2)=0.5_DP*AV(TENSOR_TO_VOIGT2(2,1),TENSOR_TO_VOIGT2(2,2))
          A(1,2,2,2)=0.5_DP*AV(TENSOR_TO_VOIGT2(1,2),TENSOR_TO_VOIGT2(2,2))          
        CASE DEFAULT
          localError="The second Voigt tensor index type of "//TRIM(NumberToVString(voigtTensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The first Voigt tensor index type of "//TRIM(NumberToVString(voigtTensorIndexTypes(1),"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(3)
      A(1,1,1,1)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(1,1))
      A(2,2,1,1)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(1,1))
      A(3,3,1,1)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(1,1))
      A(1,1,2,2)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(2,2))
      A(2,2,2,2)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(2,2))
      A(3,3,2,2)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(2,2))
      A(1,1,3,3)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(3,3))
      A(2,2,3,3)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(3,3))
      A(3,3,3,3)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(3,3))
      SELECT CASE(voigtTensorIndexTypes(1))
      CASE(TENSOR_CONTRAVARIANT_INDEX)
        SELECT CASE(voigtTensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Voigt array is (contravariant,contravariant)
          A(2,1,1,1)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,1))
          A(3,1,1,1)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,1))
          A(1,2,1,1)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,1))
          A(3,2,1,1)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,1))
          A(1,3,1,1)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,1))
          A(2,3,1,1)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,1))
          A(1,1,2,1)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(2,1))
          A(2,1,2,1)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,1))
          A(3,1,2,1)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,1))
          A(1,2,2,1)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,1))
          A(2,2,2,1)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(2,1))
          A(3,2,2,1)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,1))
          A(1,3,2,1)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,1))
          A(2,3,2,1)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,1))
          A(3,3,2,1)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(2,1))
          A(1,1,3,1)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(3,1))
          A(2,1,3,1)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,1))
          A(3,1,3,1)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,1))
          A(1,2,3,1)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,1))
          A(2,2,3,1)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(3,1))
          A(3,2,3,1)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,1))
          A(1,3,3,1)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,1))
          A(2,3,3,1)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,1))
          A(3,3,3,1)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(3,1))
          A(1,1,1,2)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(1,2))
          A(2,1,1,2)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,2))
          A(3,1,1,2)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,2))
          A(1,2,1,2)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,2))
          A(2,2,1,2)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(1,2))
          A(3,2,1,2)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,2))
          A(1,3,1,2)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,2))
          A(2,3,1,2)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,2))
          A(3,3,1,2)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(1,2))
          A(2,1,2,2)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,2))
          A(3,1,2,2)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,2))
          A(1,2,2,2)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,2))
          A(3,2,2,2)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,2))
          A(1,3,2,2)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,2))
          A(2,3,2,2)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,2))
          A(1,1,3,2)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(3,2))
          A(2,1,3,2)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,2))
          A(3,1,3,2)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,2))
          A(1,2,3,2)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,2))
          A(2,2,3,2)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(3,2))
          A(3,2,3,2)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,2))
          A(1,3,3,2)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,2))
          A(2,3,3,2)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,2))
          A(3,3,3,2)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(3,2))
          A(1,1,1,3)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(1,3))
          A(2,1,1,3)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,3))
          A(3,1,1,3)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,3))
          A(1,2,1,3)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,3))
          A(2,2,1,3)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(1,3))
          A(3,2,1,3)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,3))
          A(1,3,1,3)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,3))
          A(2,3,1,3)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,3))
          A(3,3,1,3)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(1,3))
          A(1,1,2,3)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(2,3))
          A(2,1,2,3)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,3))
          A(3,1,2,3)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,3))
          A(1,2,2,3)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,3))
          A(2,2,2,3)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(2,3))
          A(3,2,2,3)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,3))
          A(1,3,2,3)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,3))
          A(2,3,2,3)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,3))
          A(3,3,2,3)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(2,3))
          A(2,1,3,3)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,3))
          A(3,1,3,3)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,3))
          A(1,2,3,3)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,3))
          A(3,2,3,3)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,3))
          A(1,3,3,3)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,3))
          A(2,3,3,3)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,3))
        CASE(TENSOR_COVARIANT_INDEX)
          !Voigt array is (contravariant,covariant)
          A(2,1,1,1)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,1))
          A(3,1,1,1)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,1))
          A(1,2,1,1)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,1))
          A(3,2,1,1)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,1))
          A(1,3,1,1)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,1))
          A(2,3,1,1)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,1))
          A(1,1,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(2,1))
          A(2,1,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,1))
          A(3,1,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,1))
          A(1,2,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,1))
          A(2,2,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(2,1))
          A(3,2,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,1))
          A(1,3,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,1))
          A(2,3,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,1))
          A(3,3,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(2,1))
          A(1,1,3,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(3,1))
          A(2,1,3,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,1))
          A(3,1,3,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,1))
          A(1,2,3,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,1))
          A(2,2,3,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(3,1))
          A(3,2,3,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,1))
          A(1,3,3,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,1))
          A(2,3,3,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,1))
          A(3,3,3,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(3,1))
          A(1,1,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(1,2))
          A(2,1,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,2))
          A(3,1,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,2))
          A(1,2,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,2))
          A(2,2,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(1,2))
          A(3,2,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,2))
          A(1,3,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,2))
          A(2,3,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,2))
          A(3,3,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(1,2))
          A(2,1,2,2)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,2))
          A(3,1,2,2)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,2))
          A(1,2,2,2)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,2))
          A(3,2,2,2)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,2))
          A(1,3,2,2)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,2))
          A(2,3,2,2)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,2))
          A(1,1,3,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(3,2))
          A(2,1,3,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,2))
          A(3,1,3,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,2))
          A(1,2,3,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,2))
          A(2,2,3,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(3,2))
          A(3,2,3,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,2))
          A(1,3,3,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,2))
          A(2,3,3,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,2))
          A(3,3,3,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(3,2))
          A(1,1,1,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(1,3))
          A(2,1,1,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,3))
          A(3,1,1,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,3))
          A(1,2,1,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,3))
          A(2,2,1,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(1,3))
          A(3,2,1,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,3))
          A(1,3,1,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,3))
          A(2,3,1,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,3))
          A(3,3,1,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(1,3))
          A(1,1,2,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(2,3))
          A(2,1,2,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,3))
          A(3,1,2,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,3))
          A(1,2,2,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,3))
          A(2,2,2,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(2,3))
          A(3,2,2,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,3))
          A(1,3,2,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,3))
          A(2,3,2,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,3))
          A(3,3,2,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(2,3))
          A(2,1,3,3)=AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,3))
          A(3,1,3,3)=AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,3))
          A(1,2,3,3)=AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,3))
          A(3,2,3,3)=AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,3))
          A(1,3,3,3)=AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,3))
          A(2,3,3,3)=AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,3))
        CASE DEFAULT
          localError="The second Voigt tensor index type of "//TRIM(NumberToVString(voigtTensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(TENSOR_COVARIANT_INDEX)
        SELECT CASE(voigtTensorIndexTypes(2))
        CASE(TENSOR_CONTRAVARIANT_INDEX)
          !Voigt array is (covariant,contravariant)
          A(2,1,1,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,1))
          A(3,1,1,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,1))
          A(1,2,1,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,1))
          A(3,2,1,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,1))
          A(1,3,1,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,1))
          A(2,3,1,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,1))
          A(1,1,2,1)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(2,1))
          A(2,1,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,1))
          A(3,1,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,1))
          A(1,2,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,1))
          A(2,2,2,1)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(2,1))
          A(3,2,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,1))
          A(1,3,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,1))
          A(2,3,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,1))
          A(3,3,2,1)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(2,1))
          A(1,1,3,1)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(3,1))
          A(2,1,3,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,1))
          A(3,1,3,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,1))
          A(1,2,3,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,1))
          A(2,2,3,1)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(3,1))
          A(3,2,3,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,1))
          A(1,3,3,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,1))
          A(2,3,3,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,1))
          A(3,3,3,1)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(3,1))
          A(1,1,1,2)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(1,2))
          A(2,1,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,2))
          A(3,1,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,2))
          A(1,2,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,2))
          A(2,2,1,2)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(1,2))
          A(3,2,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,2))
          A(1,3,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,2))
          A(2,3,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,2))
          A(3,3,1,2)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(1,2))
          A(2,1,2,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,2))
          A(3,1,2,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,2))
          A(1,2,2,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,2))
          A(3,2,2,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,2))
          A(1,3,2,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,2))
          A(2,3,2,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,2))
          A(1,1,3,2)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(3,2))
          A(2,1,3,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,2))
          A(3,1,3,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,2))
          A(1,2,3,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,2))
          A(2,2,3,2)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(3,2))
          A(3,2,3,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,2))
          A(1,3,3,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,2))
          A(2,3,3,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,2))
          A(3,3,3,2)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(3,2))
          A(1,1,1,3)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(1,3))
          A(2,1,1,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,3))
          A(3,1,1,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,3))
          A(1,2,1,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,3))
          A(2,2,1,3)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(1,3))
          A(3,2,1,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,3))
          A(1,3,1,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,3))
          A(2,3,1,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,3))
          A(3,3,1,3)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(1,3))
          A(1,1,2,3)=AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(2,3))
          A(2,1,2,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,3))
          A(3,1,2,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,3))
          A(1,2,2,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,3))
          A(2,2,2,3)=AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(2,3))
          A(3,2,2,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,3))
          A(1,3,2,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,3))
          A(2,3,2,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,3))
          A(3,3,2,3)=AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(2,3))
          A(2,1,3,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,3))
          A(3,1,3,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,3))
          A(1,2,3,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,3))
          A(3,2,3,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,3))
          A(1,3,3,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,3))
          A(2,3,3,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,3))
        CASE(TENSOR_COVARIANT_INDEX)
          !Voigt array is (covariant,covariant)
          A(2,1,1,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,1))
          A(3,1,1,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,1))
          A(1,2,1,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,1))
          A(3,2,1,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,1))
          A(1,3,1,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,1))
          A(2,3,1,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,1))
          A(1,1,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(2,1))
          A(2,1,2,1)=0.25_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,1))
          A(3,1,2,1)=0.25_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,1))
          A(1,2,2,1)=0.25_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,1))
          A(2,2,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(2,1))
          A(3,2,2,1)=0.25_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,1))
          A(1,3,2,1)=0.25_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,1))
          A(2,3,2,1)=0.25_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,1))
          A(3,3,2,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(2,1))
          A(1,1,3,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(3,1))
          A(2,1,3,1)=0.25_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,1))
          A(3,1,3,1)=0.25_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,1))
          A(1,2,3,1)=0.25_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,1))
          A(2,2,3,1)=0.25_DP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(3,1))
          A(3,2,3,1)=0.25_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,1))
          A(1,3,3,1)=0.25_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,1))
          A(2,3,3,1)=0.25_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,1))
          A(3,3,3,1)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(3,1))
          A(1,1,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(1,2))
          A(2,1,1,2)=0.25_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,2))
          A(3,1,1,2)=0.25_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,2))
          A(1,2,1,2)=0.25_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,2))
          A(2,2,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(1,2))
          A(3,2,1,2)=0.25_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,2))
          A(1,3,1,2)=0.25_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,2))
          A(2,3,1,2)=0.25_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,2))
          A(3,3,1,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(1,2))
          A(2,1,2,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,2))
          A(3,1,2,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,2))
          A(1,2,2,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,2))
          A(3,2,2,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,2))
          A(1,3,2,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,2))
          A(2,3,2,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,2))
          A(1,1,3,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(3,2))
          A(2,1,3,2)=0.25_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,2))
          A(3,1,3,2)=0.25_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,2))
          A(1,2,3,2)=0.25_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,2))
          A(2,2,3,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(3,2))
          A(3,2,3,2)=0.25_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,2))
          A(1,3,3,2)=0.25_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,2))
          A(2,3,3,2)=0.25_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,2))
          A(3,3,3,2)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(3,2))
          A(1,1,1,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(1,3))
          A(2,1,1,3)=0.25_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(1,3))
          A(3,1,1,3)=0.25_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(1,3))
          A(1,2,1,3)=0.25_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(1,3))
          A(2,2,1,3)=0.25_DP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(1,3))
          A(3,2,1,3)=0.25_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(1,3))
          A(1,3,1,3)=0.25_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(1,3))
          A(2,3,1,3)=0.25_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(1,3))
          A(3,3,1,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(1,3))
          A(1,1,2,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,1),TENSOR_TO_VOIGT3(2,3))
          A(2,1,2,3)=0.25_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(2,3))
          A(3,1,2,3)=0.25_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(2,3))
          A(1,2,2,3)=0.25_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(2,3))
          A(2,2,2,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,2),TENSOR_TO_VOIGT3(2,3))
          A(3,2,2,3)=0.25_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(2,3))
          A(1,3,2,3)=0.25_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(2,3))
          A(2,3,2,3)=0.25_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(2,3))
          A(3,3,2,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,3),TENSOR_TO_VOIGT3(2,3))
          A(2,1,3,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,1),TENSOR_TO_VOIGT3(3,3))
          A(3,1,3,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,1),TENSOR_TO_VOIGT3(3,3))
          A(1,2,3,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,2),TENSOR_TO_VOIGT3(3,3))
          A(3,2,3,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(3,2),TENSOR_TO_VOIGT3(3,3))
          A(1,3,3,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(1,3),TENSOR_TO_VOIGT3(3,3))
          A(2,3,3,3)=0.5_DP*AV(TENSOR_TO_VOIGT3(2,3),TENSOR_TO_VOIGT3(3,3))
        CASE DEFAULT
          localError="The second Voigt tensor index type of "//TRIM(NumberToVString(voigtTensorIndexTypes(2),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The first Voigt tensor index type of "//TRIM(NumberToVString(voigtTensorIndexTypes(1),"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The number of dimensions of "//TRIM(NumberToVString(n,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT    

    EXITS("VoigtToTensor4DP")
    RETURN
999 ERRORSEXITS("VoigtToTensor4DP",err,error)
    RETURN 1
    
  END SUBROUTINE VoigtToTensor4DP

  !
  !================================================================================================================================
  !

  !>Returns a zero matrix
  SUBROUTINE ZeroMatrixSP(A,err,error,*)
    
    !Argument variables
    REAL(SP), INTENT(OUT) :: A(:,:) !<On exit, the zero matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("ZeroMatrixSP",err,error,*999)

    A=0.0_SP

    EXITS("ZeroMatrixSP")
    RETURN
999 ERRORSEXITS("ZeroMatrixSP",err,error)
    RETURN 1
    
  END SUBROUTINE ZeroMatrixSP

  !
  !================================================================================================================================
  !

  !>Returns a zero matrix
  SUBROUTINE ZeroMatrixDP(A,err,error,*)
    
    !Argument variables
    REAL(DP), INTENT(OUT) :: A(:,:) !<On exit, the zero matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
   
    ENTERS("ZeroMatrixDP",err,error,*999)

    A=0.0_DP

    EXITS("ZeroMatrixDP")
    RETURN
999 ERRORSEXITS("ZeroMatrixDP",err,error)
    RETURN 1
    
  END SUBROUTINE ZeroMatrixDP

  !
  !================================================================================================================================
  !

  !> Calculates second derivatives of a cubic spline function for a tabulated function y(x). Call  spline_cubic_val to evaluate at t values.
  !> algorithm adapted from John Burkhardt's spline_cubic_set routine from the SPLINE package (http://people.sc.fsu.edu/~jburkardt/f_src/spline/spline.html)
  SUBROUTINE spline_cubic_set (n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp, err, error, *)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: n !< size of x,y arrays to interpolate values from
    REAL(DP), INTENT(IN) :: t(n) !< t array: known values
    REAL(DP), INTENT(IN) :: y(n) !< y array: values to interpolate
    INTEGER(INTG), INTENT(IN) :: ibcbeg !< left boundary condition flag
    REAL(DP), INTENT(IN) :: ybcbeg !< 1st derivative interpolating function at point 1 (left boundary)
    INTEGER(INTG), INTENT(IN) :: ibcend !< right boundary condition flag
    REAL(DP), INTENT(IN) :: ybcend !< 1st derivative interpolating function at point n (right boundary)
    REAL(DP), INTENT(OUT) :: ypp(n) !< 2nd derivatives of interpolating function at x values
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    REAL(DP) :: diag(n)
    REAL(DP) :: sub(2:n)
    REAL(DP) :: sup(1:n-1)
    INTEGER(INTG) :: i
    TYPE(VARYING_STRING) :: localError

    ENTERS("spline_cubic_set",err,error,*999)

    ! Sanity checks
    IF ( n <= 1 ) then
      localError="spline interpolation requires at least 2 knots- user supplied "//TRIM(NumberToVString(n,"*",err,error))
      CALL FlagError(localError,err,error,*999)
    ENDIF
    DO i = 1, n-1
      IF ( t(i) >= t(i+1) ) then
        localError="Non-increasing knots supplied for cubic spline interpolation."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDDO

    !  Set the first equation.
    IF ( ibcbeg == 0 ) then
      ypp(1) = 0.0E+00_DP
      diag(1) = 1.0E+00_DP
      sup(1) = -1.0E+00_DP
    ELSE IF ( ibcbeg == 1 ) then
      ypp(1) = ( y(2) - y(1) ) / ( t(2) - t(1) ) - ybcbeg
      diag(1) = ( t(2) - t(1) ) / 3.0E+00_DP
      sup(1) = ( t(2) - t(1) ) / 6.0E+00_DP
    ELSE IF ( ibcbeg == 2 ) then
      ypp(1) = ybcbeg
      diag(1) = 1.0E+00_DP
      sup(1) = 0.0E+00_DP
    ELSE
      localError="The boundary flag IBCBEG must be 0, 1 or 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    !  Set the intermediate equations.
    DO i = 2, n-1
      ypp(i) = ( y(i+1) - y(i) ) / ( t(i+1) - t(i) ) &
       &     - ( y(i) - y(i-1) ) / ( t(i) - t(i-1) )
      sub(i) = ( t(i) - t(i-1) ) / 6.0E+00_DP
      diag(i) = ( t(i+1) - t(i-1) ) / 3.0E+00_DP
      sup(i) = ( t(i+1) - t(i) ) / 6.0E+00_DP
    ENDDO

    !  Set the last equation.
    IF ( ibcend == 0 ) then
      ypp(n) = 0.0E+00_DP
      sub(n) = -1.0E+00_DP
      diag(n) = 1.0E+00_DP
    ELSE IF ( ibcend == 1 ) then
      ypp(n) = ybcend - ( y(n) - y(n-1) ) / ( t(n) - t(n-1) )
      sub(n) = ( t(n) - t(n-1) ) / 6.0E+00_DP
      diag(n) = ( t(n) - t(n-1) ) / 3.0E+00_DP
    ELSE IF ( ibcend == 2 ) then
      ypp(n) = ybcend
      sub(n) = 0.0E+00_DP
      diag(n) = 1.0E+00_DP
    ELSE
      localError="The boundary flag IBCEND must be 0, 1 or 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    !  Special case:
    !    N = 2, IBCBEG = IBCEND = 0.
    IF ( n == 2 .and. ibcbeg == 0 .and. ibcend == 0 ) then
      ypp(1) = 0.0E+00_DP
      ypp(2) = 0.0E+00_DP

    !  Solve the linear system.
    ELSE
      CALL s3_fs ( sub, diag, sup, n, ypp, ypp, err, error, *999 )
    ENDIF

    EXITS("spline_cubic_set")
    RETURN
999 ERRORSEXITS("spline_cubic_set",err,error)
    RETURN 1
    
  END SUBROUTINE spline_cubic_set

  !
  !================================================================================================================================
  !

  !> S3_FS factors and solves a tridiagonal linear system.
  !> algorithm adapted from John Burkhardt's s3_fs routine from the SPLINE package (http://people.sc.fsu.edu/~jburkardt/f_src/spline/spline.html)
  SUBROUTINE s3_fs ( a1, a2, a3, n, b, x, err, error, *)

    !Argument variables
    REAL(DP), INTENT(INOUT) :: a1(2:n) !< IN: nonzero diagonal of linear system OUT: factorization info
    REAL(DP), INTENT(INOUT) :: a2(1:n) !< IN: nonzero diagonal of linear system OUT: factorization info
    REAL(DP), INTENT(INOUT) :: a3(1:n-1) !< IN: nonzero diagonal of linear system OUT: factorization info
    INTEGER(INTG), INTENT(IN) :: n !< size of x,y arrays to interpolate values from
    REAL(DP), INTENT(INOUT) :: b(n) !< IN: RHS of linear system OUT: factorization info
    REAL(DP), INTENT(OUT) :: x(n) !< solution of linear system
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i
    REAL(DP) :: xmult
    TYPE(VARYING_STRING) :: localError

    ENTERS("s3_fs",err,error,*999)

    !  The diagonal entries can't be zero.
    DO i = 1, n
      IF ( ABS(a2(i)) < ZERO_TOLERANCE ) then
        localError="Zero diagonal entry in tridiagonal linear system."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDDO

    DO i = 2, n-1
      xmult = a1(i) / a2(i-1)
      a2(i) = a2(i) - xmult * a3(i-1)
      b(i) = b(i) - xmult * b(i-1)
    ENDDO

    xmult = a1(n) / a2(n-1)
    a2(n) = a2(n) - xmult * a3(n-1)
    x(n) = ( b(n) - xmult * b(n-1) ) / a2(n)
    DO i = n-1, 1, -1
      x(i) = ( b(i) - a3(i) * x(i+1) ) / a2(i)
    ENDDO

    EXITS("s3_fs")
    RETURN
999 ERRORSEXITS("s3_fs",err,error)
    RETURN 1
    
  END SUBROUTINE s3_fs

  !
  !================================================================================================================================
  !
  
  !> Evaluates a cubic spline at a specified point. First call spline_cubic_set to calculate derivatives
  !> algorithm adapted from John Burkhardt's spline_cubic_val routine from the SPLINE package (http://people.sc.fsu.edu/~jburkardt/f_src/spline/spline.html)
  SUBROUTINE spline_cubic_val (n, t, y, ypp, tval, yval, ypval, yppval, err, error, *)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: n !< size of t,y arrays to interpolate values from
    REAL(DP), INTENT(IN) :: t(n) !< t array: known knot values
    REAL(DP), INTENT(IN) :: y(n) !< y array: data values to interpolate at the knots
    REAL(DP), INTENT(IN) :: ypp(n) !< 2nd derivatives of interpolating function at t values
    REAL(DP), INTENT(IN) :: tval !< point in t at which spline is to be evaluated
    REAL(DP), INTENT(OUT) :: yval !< spline interpolated y value at tval
    REAL(DP), INTENT(OUT) :: ypval !< first derivative of spline interpolated y value at tval
    REAL(DP), INTENT(OUT) :: yppval !< second derivative of spline interpolated y value at tval
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    REAL(DP) :: dt
    REAL(DP) :: h
    INTEGER(INTG) :: i
    INTEGER(INTG) :: left
    INTEGER(INTG) :: right
    LOGICAL :: foundInterval

    ENTERS("spline_cubic_val",err,error,*999)

    !  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
    !  Values below T(1) or above T(N) use extrapolation.
    foundInterval = .FALSE.
    DO i = 2, n - 1
      IF ( tval < t(i) ) THEN
        foundInterval=.TRUE.
        left = i - 1
        right = i
        EXIT
      ENDIF
    ENDDO
    IF (foundInterval .EQV. .FALSE.) THEN
      left = n - 1
      right = n
    ENDIF

    !  Evaluate the polynomial.
    dt = tval - t(left)
    h = t(right) - t(left)

    yval = y(left) &
     &   + dt * ( ( y(right) - y(left) ) / h &
     &          - ( ypp(right) / 6.0E+00_DP + ypp(left) / 3.0E+00_DP ) * h &
     &   + dt * ( 0.5E+00 * ypp(left) &
     &   + dt * ( ( ypp(right) - ypp(left) ) / ( 6.0E+00_DP * h ) ) ) )

    ypval = ( y(right) - y(left) ) / h &
     &   - ( ypp(right) / 6.0E+00_DP + ypp(left) / 3.0E+00_DP ) * h &
     &   + dt * ( ypp(left) &
     &   + dt * ( 0.5E+00_DP * ( ypp(right) - ypp(left) ) / h ) )

    yppval = ypp(left) + dt * ( ypp(right) - ypp(left) ) / h 

    EXITS("spline_cubic_val")
    RETURN
999 ERRORSEXITS("spline_cubic_val",err,error)
    RETURN 1
    
  END SUBROUTINE spline_cubic_val

  !
  !================================================================================================================================
  !

END MODULE Maths

