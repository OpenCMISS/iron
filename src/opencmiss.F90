!> \file
!> \author Chris Bradley
!> \brief The top level OpenCMISS module.
!>
!> \mainpage OpenCMISS Documentation
!>
!> An open source interactive computer program for Continuum Mechanics, Image analysis, Signal processing and System
!> Identification. Target usage: Bioengineering application of finite element analysis, boundary element and collocation
!> techniques.
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
!>

!> \defgroup OpenCMISS
!> The top level OpenCMISS module. This module is the buffer Fortran module between the OpenCMISS library and user code.
MODULE OpenCMISS
  
  USE KINDS

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !> \addtogroup OpenCMISS_KindConstants OpenCMISS::Constants::Kinds
  !> \brief Kind constants.
  !>@{
  !> \addtogroup OpenCMISS_IntegerKinds OpenCMISS::Constants::Kinds::Integers
  !> \brief Kind parameters for integer data types.
  !> \see OpenCMISS::Constants::Kinds,OpenCMISS
  !>@{
  INTEGER, PARAMETER :: CMISSIntg = INTG !<Standard integer kind. \see OpenCMISS_IntegerKinds,OpenCMISS
  INTEGER, PARAMETER :: CMISSSIntg = SINTG !<Short integer kind. \see OpenCMISS_IntegerKinds,OpenCMISS
  INTEGER, PARAMETER :: CMISSLIntg = LINTG !<Long integer kind. \see OpenCMISS_IntegerKinds,OpenCMISS
  INTEGER, PARAMETER :: CMISSPtr = PTR !<Pointer integer kind. \see OpenCMISS_IntegerKinds,OpenCMISS
  INTEGER, PARAMETER :: CMISSIdx = IDX !<Index integer kind. \see OpenCMISS_IntegerKinds,OpenCMISS
  INTEGER, PARAMETER :: CMISSLIdx = LIDX !<Long index integer kind. \see OpenCMISS_IntegerKinds,OpenCMISS
  !>@}
  !> \addtogroup OpenCMISS_RealKinds OpenCMISS::Constants::Kinds::Reals
  !> \brief Kind parameters for real data types.
  !> \see OpenCMISS::Constants::Kinds,OpenCMISS
  !>@{
  INTEGER, PARAMETER :: CMISSSP = SP !<Single precision real kind. \see OpenCMISS_RealKinds,OpenCMISS
  INTEGER, PARAMETER :: CMISSDP = DP !<Double precision real kind. \see OpenCMISS_RealKinds,OpenCMISS
  INTEGER, PARAMETER :: CMISSQP = QP !<Quadruple precision real kind. \see OpenCMISS_RealKinds,OpenCMISS
  INTEGER, PARAMETER :: CMISSRP = RP !<Working precision real kind. \see OpenCMISS_RealKinds,OpenCMISS
  !>@}
  !> \addtogroup OpenCMISS_ComplexKinds OpenCMISS::Constants::Kinds::Complex
  !> \brief Kind parameters for complex data types
  !> \see OpenCMISS::Constants::Kinds,OpenCMISS
  !>@{
  INTEGER, PARAMETER :: CMISSSPC = SPC !<Single precision complex kind. \see OpenCMISS_ComplexKinds,OpenCMISS
  INTEGER, PARAMETER :: CMISSDPC = DPC !<Double precision complex kind. \see OpenCMISS_ComplexKinds,OpenCMISS
  INTEGER, PARAMETER :: CMISSRPC = RPC !<Working precision complex kind. \see OpenCMISS_ComplexKinds,OpenCMISS
  !>@}
  !>@}

  !Module types

  !Module variables

  !Interfaces

  PUBLIC CMISSIntg,CMISSSIntg,CMISSLIntg,CMISSPtr,CMISSIdx,CMISSLIdx

  PUBLIC CMISSSP,CMISSDP,CMISSQP,CMISSRP

  PUBLIC CMISSSPC,CMISSDPC,CMISSRPC

CONTAINS

END MODULE OpenCMISS
