!> \file
!> \author 
!> \brief This module contains some hardcoded cell models and integration routines for cardiac electrophysiology.
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
!> Contributor(s): Sander Land, Chris Bradley
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

!>This module handles all hard coded electrophysiology cell routines.
MODULE ElectrophysiologyCellRoutines
  
  USE BaseRoutines
  USE Constants
  USE DecompositionAccessRoutines
  USE FieldAccessRoutines
  USE ISO_VARYING_STRING
  use Kinds
  use Strings
  use Types

#include "macros.h"  

  IMPLICIT NONE
  
  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC Electrophysiology_BuenoOrovioInitialise

  PUBLIC Electrophysiology_BuenoOrovioIntegrate
  
  PUBLIC Electrophysiology_TenTusscher06Initialise

  PUBLIC Electrophysiology_TenTusscher06Integrate

  INTERFACE Pow
    MODULE PROCEDURE PowDP
    MODULE PROCEDURE PowIntg
  END INTERFACE Pow

CONTAINS
  
  !
  !================================================================================================================================
  !

  ! for auto generated code

  REAL(DP) FUNCTION PowDP(a,b)
    REAL(DP), INTENT(IN) :: a,b
    PowDP = a**b
  END FUNCTION PowDP
  
  REAL(DP) FUNCTION PowIntg(a,b)
    REAL(DP), INTENT(IN) :: a
    INTEGER(INTG), INTENT(IN) :: b
    PowIntg = a**b
  END FUNCTION PowIntg

  !
  !================================================================================================================================
  !

  ! Bueno-Orovio (2008) minimal cell model. Membrane potential rescaled, partial evaluation applied.

  !>Initialize a field for the Bueno-Orovio 2008 cell model.
  SUBROUTINE Electrophysiology_BuenoOrovioInitialise(field,err,error,*)
    
    !Argument variables
    TYPE(FieldType), INTENT(INOUT), POINTER :: field !<The field to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables    
    INTEGER(INTG) :: componentIdx
    REAL(DP), PARAMETER, DIMENSION(1:4) :: y0 = [ -8.09242e+01,  9.99993e-01, 2.16366e-02, 9.84320e-01 ] ! paced initial condition

    ENTERS("Electrophysiology_BuenoOrovioInitialise",err,error,*999)
    
    DO componentIdx=1,4
      CALL Field_ComponentValuesInitialise(field,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,componentIdx,y0(componentIdx), &
        & err,error,*999)
    ENDDO !componentIdx

    EXITS("Electrophysiology_BuenoOrovioInitialise")
    RETURN
999 ERRORSEXITS("Electrophysiology_BuenoOrovioInitialise",err,error)
    RETURN 1

  END SUBROUTINE Electrophysiology_BuenoOrovioInitialise

  !
  !================================================================================================================================
  !

  !>Evaluate the rates for a Bueno Orovio cell model
  SUBROUTINE Electrophysiology_BuenoOroviodYdt(t,y,dydt,activationFactor)
    
    !Argument variables
    REAL(DP), INTENT(IN) :: y(:) !<Current state
    REAL(DP), INTENT(IN) :: t !<Current time 
    REAL(DP), INTENT(OUT) :: dydt(:) !<On exit, Derivatives
    REAL(DP), INTENT(IN) :: activationFactor !<Activation factor. 1 for default stimulus current, 0 for none
    !Local Variables
    REAL(DP), PARAMETER :: u_m=3.0e-01, u_p=1.3e-01, u_q=6.0e-03, u_r=6.0e-03, vscale=81.0, period = 1000.0
    REAL(DP) :: J_fi, J_si, J_so, Vm, bv, m, p, q, r, tau_o, tau_s, tau_so, tau_v_minus, tau_w_minus, w_inf, J_stim
  
    IF(MOD(t,period)>=0.AND.MOD(t,period)<=0.999999) THEN ! apply stimulus in first ms
      J_stim = -0.65 * activationFactor ! stimulus: ~52 mV for 1 ms
    ELSE
      J_stim  = 0
    ENDIF

    bv = (1+(Y(1)/vscale)); !< rescale, original Y(1) in model
    IF(bv < u_m) THEN
      m = 0
    ELSE
      m = 1
    ENDIF
    IF(bv < u_q) THEN
      q = 0
    ELSE
      q = 1
    ENDIF
    IF(bv < u_p) THEN
      p = 0
    ELSE
      p = 1
    ENDIF

    IF(bv < u_r) THEN
      r = 0
    ELSE
      r = 1
    ENDIF

    J_fi = (9.09090909090909e+00*(-m)*Y(2)*(-3.0e-01+bv)*(1.55e+00-bv))
    tau_v_minus = ((1150*q)+(60*(1-q)))
    dydt(2) = (((1-m)*(0-Y(2))/tau_v_minus)-(6.89655172413793e-01*m*Y(2)))
    tau_o = ((400*(1-r))+(6*r))
    tau_so = (3.002e+01+(-1.4512e+01*(1+TANH((2.046e+00*(-6.5e-01+bv))))))
    J_so = ((bv*(1-p)/tau_o)+(p/tau_so))
    J_si = (5.29801324503311e-01*(-p)*Y(4)*Y(3))
    dydt(1) = ((-J_fi-J_so-J_si-J_stim)*vscale)
    Vm = (-83+(8.57e+01*bv))
    tau_s = ((2.7342e+00*(1-p))+(16*p))
    dydt(3) = (((5.0e-01*(1+TANH((2.0994e+00*(-9.087e-01+bv)))))-Y(3))/tau_s)
    w_inf = (((1-r)*(1-(1.42857142857143e+01*bv)))+(9.4e-01*r))
    tau_w_minus = (60+(-2.25e+01*(1+TANH((65*(-3.0e-02+bv))))))
    dydt(4) = (((1-r)*(w_inf-Y(4))/tau_w_minus)-(5.0e-03*r*Y(4)))
    
  END SUBROUTINE Electrophysiology_BuenoOroviodYdt

  !
  !================================================================================================================================
  !

  !>Integrates a Bueno Orovio cell model
  SUBROUTINE Electrophysiology_BuenoOrovioIntegrate(cellsField,materialsField,t0,t1,err,error,*)
    
    !Argument variables    
    TYPE(FieldType), INTENT(INOUT), POINTER :: cellsField !<Independent field storing the cell data
    TYPE(FieldType), INTENT(INOUT), POINTER :: materialsField !<Material field component 1 the activation flags
    REAL(DP), INTENT(IN)    :: t0, t1 !<Integrate from time t0 to t1 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: cellIdx,componentIdx,localDOFIdx,numberOfCells
    INTEGER(INTG), PARAMETER :: cellDimension = 4
    REAL(DP) :: activationFactor,dt,dydt(1:cellDimension),t,y(1:cellDimension) 
    REAL(DP), POINTER :: cellData(:), activationData(:)
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(FieldType), POINTER :: activationField
    TYPE(FieldVariableType), POINTER :: activationVariable,cellsVariable

    ENTERS('Electrophysiology_BuenoOrovioIntegrate',err,error,*999)

    NULLIFY(decomposition)
    CALL Field_DecompositionGet(cellsField,decomposition,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,1,domain,err,error,*999)
    NULLIFY(domainTopology)
    CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(domainNodes)
    CALL DomainTopoloy_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
    CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfCells,err,error,*999)

    NULLIFY(cellsVariable)
    CALL Field_VariableGet(cellsField,FIELD_U_VARIABLE_TYPE,cellsVariable,err,error,*999)
    NULLIFY(activationVariable)
    CALL Field_VariableGet(activationField,FIELD_U_VARIABLE_TYPE,activationVariable,err,error,*999)
 
    CALL FieldVariable_ParameterSetDataGet(cellsVariable,FIELD_VALUES_SET_TYPE,cellData,err,error,*999)
    CALL FieldVariable_ParameterSetDataGet(activationVariable,FIELD_VALUES_SET_TYPE,activationData,err,error,*999)

    DO cellIdx=1,numberOfCells
      !   field ->   y
      DO componentIdx=1,cellDimension
        !Default to version 1 of each node derivative
        CALL FieldVariable_LocalNodeDOFGet(cellsVariable,1,1,cellIdx,componentIdx,localDOFIdx,err,error,*999)
        y(componentIdx) = cellData(localDOFIdx)
      ENDDO !componentIdx
      !Default to version 1 of each node derivative
      CALL FieldVariable_LocalNodeDOFGet(activationVariable,1,1,cellIdx,1,localDOFIdx,err,error,*999)
      activationFactor = activationData(localDOFIdx)
      
      t = t0
      DO WHILE (t < (t1 - ZERO_TOLERANCE))
        ! integrate one cell, one time step. 
        ! just use adaptive forward euler, tested in matlab
        call Electrophysiology_BuenoOroviodYdt(t,y,dydt,activationFactor)
        dt = MIN(t1-t,0.5)  ! default max
        dt = MIN(dt,1 / ABS(dydt(1))) ! maximum increase in Vm
        y = y + dt * dydt ! FWE
        t = t + dt
      ENDDO
      !   y -> field  
      DO componentIdx=1,cellDimension
        !Default to version 1 of each node derivative
        CALL FieldVariable_LocalNodeDOFGet(cellsVariable,1,1,cellIdx,componentIdx,localDOFIdx,err,error,*999)
        cellData(localDOFIdx) = y(componentIdx)
      ENDDO !componentIdx
      
    ENDDO !cellIdx
    
    CALL FieldVariable_ParameterSetDataRestore(cellsVariable,FIELD_VALUES_SET_TYPE,cellData,err,error,*999)
    CALL FieldVariable_ParameterSetDataRestore(activationVariable,FIELD_VALUES_SET_TYPE,activationData,err,error,*999)

    EXITS('Electrophysiology_BuenoOrovioIntegrate')
    RETURN
999 ERRORSEXITS('Electrophysiology_BuenoOrovioIntegrate',err,error)
    RETURN 1
    
  END SUBROUTINE Electrophysiology_BuenoOrovioIntegrate
  
  !
  !================================================================================================================================
  !

  ! Ten Tusscher and Panfilov (2006) human ventricular cell model

  !>Initialize a filed for the Ten Tusscher and Panfilov 2006 cell model
  SUBROUTINE Electrophysiology_TenTusscher06Initialise(field,err,error,*)
    
    !Argument variables
    TYPE(FieldType), INTENT(INOUT), POINTER :: field !<The field to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables    
    INTEGER(INTG) :: componentIdx
    REAL(DP), PARAMETER, DIMENSION(1:19) :: y0 = [ -85.23, 0.9755, 0.9953, 0.7888, 3.64, 0.000126, 0.00036, 0.9073, 0.7444,&
        & 0.7045, 0.00172,  3.373e-5, 136.89, 0.00621, 0.4712, 0.0095, 8.604, 2.42e-8, 0.999998 ] ! paced initial condition

    ENTERS("Electrophysiology_TenTusscher06Initialise",err,error,*999)
    
    DO componentIdx=1,19
      CALL Field_ComponentValuesInitialise(field,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,componentIdx,y0(componentIdx), &
        & err,error,*999)
    ENDDO !componentIdx

    EXITS("Electrophysiology_TenTusscher06Initialise")
    RETURN
999 ERRORSEXITS("Electrophysiology_TenTusscher06Initialise",err,error)
    RETURN 1

  END SUBROUTINE Electrophysiology_TenTusscher06Initialise
  
  !
  !================================================================================================================================
  !

  !>Evaluate the rates for a Ten Tusscher 2006 cell model
  SUBROUTINE Electrophysiology_TenTusscher06dYdt(t,y,dydt,activationFactor)

    !Argument variables    
    REAL(DP), INTENT(IN) :: t !<Current time 
    REAL(DP), dimension(:), INTENT(IN) :: y(:) !<Current state
    REAL(DP), DIMENSION(:), INTENT(OUT) :: dydt(:) !<On exit, the derivatives
    REAL(DP), INTENT(IN) :: activationFactor  !<Activation factor. 1 for default stimulus current, 0 for none
    !Local Variables
    REAL(DP), PARAMETER :: period=1000
    REAL(DP) :: Ca_i_bufc, Ca_sr_bufsr, Ca_ss_bufss, E_Ca, E_K, E_Ks, E_Na, O, alpha_K1, alpha_d, alpha_m, alpha_xr1, alpha_xr2,&
    & alpha_xs, beta_K1, beta_d, beta_m, beta_xr1, beta_xr2, beta_xs, d_inf, f2_inf, fCass_inf, f_inf, gamma_d, h_inf, i_CaL,&
    & i_K1, i_Kr, i_Ks, i_Na, i_NaCa, i_NaK, i_Stim, i_b_Ca, i_b_Na, i_leak, i_p_Ca, i_p_K, i_rel, i_to, i_up, i_xfer, j_inf, k1,&
    & k2, kcasr, m_inf, r_inf, s_inf, tau_d, tau_f, tau_f2, tau_fCass, tau_h, tau_j, tau_m, tau_r, tau_s, tau_xr1, tau_xr2,tau_xs,&
    & xK1_inf, xr1_inf, xr2_inf, xs_inf

    IF(MOD(t,period)>=0.AND.MOD(t,period)<=1.999999) THEN ! apply stimulus in 2 ms
      i_Stim = -3.57142857142857e+01 * activationFactor
    ELSE
      i_Stim = 0
    ENDIF

    i_CaL = (5.75002020961245e-01*Y(12)*Y(4)*Y(2)*Y(3)*(-15+Y(1))*(-2+(2.5e-01*Y(7)*EXP((7.48677816454909e-02*(-15+Y(1))))))/&
      & (-1+EXP((7.48677816454909e-02*(-15+Y(1))))))
    
    d_inf = (1/(1+EXP((1.33333333333333e-01*(-8-Y(1))))))
    alpha_d = (2.5e-01+(1.4e+00/(1+EXP((7.69230769230769e-02*(-35-Y(1)))))))
    beta_d = (1.4e+00/(1+EXP((2.0e-01*(5+Y(1))))))
    gamma_d = (1/(1+EXP((5.0e-02*(50-Y(1))))))
    tau_d = ((1*alpha_d*beta_d)+gamma_d)
    f2_inf = (3.3e-01+(6.7e-01/(1+EXP((1.42857142857143e-01*(35+Y(1)))))))
    tau_f2 = ((562*EXP((4.16666666666667e-03*(-pow((27+Y(1)),2)))))+(31/(1+EXP((1.0e-01*(25-Y(1))))))+&
      & (80/(1+EXP((1.0e-01*(30+Y(1)))))))
    dydt(2) = ((f2_inf-Y(2))/tau_f2)
    fCass_inf = (4.0e-01+(6.0e-01/(1+pow((20*Y(7)),2))))
    tau_fCass = (2+(80/(1+pow((20*Y(7)),2))))
    dydt(3) = ((fCass_inf-Y(3))/tau_fCass)
    f_inf = (1/(1+EXP((1.42857142857143e-01*(20+Y(1))))))
    tau_f = (20+(1.1025e+03*EXP((4.44444444444444e-03*(-pow((27+Y(1)),2)))))+(200/(1+EXP((1.0e-01*(13-Y(1))))))+&
      &(180/(1+EXP((1.0e-01*(30+Y(1)))))))
    dydt(4) = ((f_inf-Y(4))/tau_f)
    E_Ca = (1.33568803298478e+01*LOG((2/Y(6))))
    i_b_Ca = (5.92e-04*(Y(1)-E_Ca))
    kcasr = (2.5e+00-(1.5e+00/(1+pow((1.5e+00/Y(5)),2))))
    k1 = (1.5e-01/kcasr)
    O = (k1*pow(Y(7),2)*Y(8)/(6.0e-02+(k1*pow(Y(7),2))))
    i_rel = (1.02e-01*O*(Y(5)-Y(7)))
    i_up = (6.375e-03/(1+(6.25e-08/pow(Y(6),2))))
    i_leak = (3.6e-04*(Y(5)-Y(6)))
    i_xfer = (3.8e-03*(Y(7)-Y(6)))
    k2 = (4.5e-02*kcasr)
    dydt(8) = (((-k2)*Y(7)*Y(8))+(5.0e-03*(1-Y(8))))
    Ca_i_bufc = (1/(1+(2.0e-04/pow((1.0e-03+Y(6)),2))))
    Ca_sr_bufsr = (1/(1+(3/pow((3.0e-01+Y(5)),2))))
    Ca_ss_bufss = (1/(1+(1.0e-04/pow((2.5e-04+Y(7)),2))))
    i_p_Ca = (1.238e-01*Y(6)/(5.0e-04+Y(6)))
    i_NaCa = (8.66622022994244e-05*((2*EXP((1.31018617879609e-02*Y(1)))*pow(Y(17),3))-(6860000*&
      &exp((-2.43320290347846e-02*Y(1)))*Y(6)))/(1+(1.0e-01*EXP((-2.43320290347846e-02*Y(1))))))
    dydt(6) = (Ca_i_bufc*((6.66910509631797e-02*(i_leak-i_up))+i_xfer-(5.84427487220097e-05*(i_b_Ca+i_p_Ca-(2*i_NaCa)))))
    dydt(5) = (Ca_sr_bufsr*(i_up-i_rel-i_leak))
    dydt(7) = (Ca_ss_bufss*((-1.75328246166029e-02*i_CaL)+(2.00073152889539e+01*i_rel)-(300*i_xfer)))
    E_Na = (2.67137606596956e+01*LOG((140/Y(17))))
    i_Na = (1.4838e+01*pow(Y(11),3)*Y(9)*Y(10)*(Y(1)-E_Na))
    h_inf = (1/pow((1+EXP((1.34589502018843e-01*(7.155e+01+Y(1))))),2))
    IF(Y(1) < -40.0) THEN
      tau_h = (1/((5.7e-02*EXP((1.47058823529412e-01*(-80-Y(1)))))+(2.7e+00*EXP((7.9e-02*Y(1))))+(310000*EXP((3.485e-01*Y(1))))))
    ELSE
      tau_h = (1.68831168831169e-01*(1+EXP((-9.00900900900901e-02*(1.066e+01+Y(1))))))
    ENDIF
    dydt(9) = ((h_inf-Y(9))/tau_h)
    j_inf = (1/pow((1+EXP((1.34589502018843e-01*(7.155e+01+Y(1))))),2))
    IF(Y(1) < -40.0) THEN
      tau_j = (1/((1*((-25428*EXP((2.444e-01*Y(1))))-(6.948e-06*EXP((-4.391e-02*Y(1)))))*(3.778e+01+Y(1))/&
        &(1+EXP((3.11e-01*(7.923e+01+Y(1))))))+(2.424e-02*EXP((-1.052e-02*Y(1)))/(1+EXP((-1.378e-01*(4.014e+01+Y(1))))))))
    ELSE
      tau_j = (1.66666666666667e+00/EXP((5.7e-02*Y(1)))*(1+EXP((-1.0e-01*(32+Y(1))))))
    ENDIF
    dydt(10) = ((j_inf-Y(10))/tau_j)
    m_inf = (1/pow((1+EXP((1.10741971207087e-01*(-5.686e+01-Y(1))))),2))
    alpha_m = (1/(1+EXP((2.0e-01*(-60-Y(1))))))
    beta_m = ((1.0e-01/(1+EXP((2.0e-01*(35+Y(1))))))+(1.0e-01/(1+EXP((5.0e-03*(-50+Y(1)))))))
    tau_m = (1*alpha_m*beta_m)
    dydt(11) = ((m_inf-Y(11))/tau_m)
    E_K = (2.67137606596956e+01*LOG((5.4e+00/Y(13))))
    alpha_K1 = (1.0e-01/(1+EXP((6.0e-02*(-200+Y(1)-E_K)))))
    beta_K1 = (((3*EXP((2.0e-04*(100+Y(1)-E_K))))+EXP((1.0e-01*(-10+Y(1)-E_K))))/(1+EXP((-5.0e-01*(Y(1)-E_K)))))
    xK1_inf = (alpha_K1/(alpha_K1+beta_K1))
    i_K1 = (5.405e+00*xK1_inf*(Y(1)-E_K))
    i_to = (2.94e-01*Y(18)*Y(19)*(Y(1)-E_K))
    i_Kr = (1.53e-01*Y(14)*Y(15)*(Y(1)-E_K))
    E_Ks = (2.67137606596956e+01*LOG((9.6e+00/(Y(13)+(3.0e-02*Y(17))))))
    i_Ks = (3.92e-01*pow(Y(16),2)*(Y(1)-E_Ks))
    i_NaK = (2.298375e+00*Y(17)/(40+Y(17))/(1+(1.245e-01*EXP((-3.74338908227455e-03*Y(1))))+&
      &(3.53e-02*EXP((3.74338908227455e-02*(-Y(1)))))))
    i_b_Na = (2.9e-04*(Y(1)-E_Na))
    i_p_K = (1.46e-02*(Y(1)-E_K)/(1+EXP((1.67224080267559e-01*(25-Y(1))))))
    dydt(12) = ((d_inf-Y(12))/tau_d)
    dydt(1) = (-1*(i_K1+i_to+i_Kr+i_Ks+i_CaL+i_NaK+i_Na+i_b_Na+i_NaCa+i_b_Ca+i_p_K+i_p_Ca+i_Stim))
    dydt(13) = (-1.16885497444019e-04*(i_K1+i_to+i_Kr+i_Ks+i_p_K+i_Stim-(2*i_NaK)))
    xr1_inf = (1/(1+EXP((1.42857142857143e-01*(-26-Y(1))))))
    alpha_xr1 = (450/(1+EXP((1.0e-01*(-45-Y(1))))))
    beta_xr1 = (6/(1+EXP((8.69565217391304e-02*(30+Y(1))))))
    tau_xr1 = (1*alpha_xr1*beta_xr1)
    dydt(14) = ((xr1_inf-Y(14))/tau_xr1)
    xr2_inf = (1/(1+EXP((4.16666666666667e-02*(88+Y(1))))))
    alpha_xr2 = (3/(1+EXP((5.0e-02*(-60-Y(1))))))
    beta_xr2 = (1.12e+00/(1+EXP((5.0e-02*(-60+Y(1))))))
    tau_xr2 = (1*alpha_xr2*beta_xr2)
    dydt(15) = ((xr2_inf-Y(15))/tau_xr2)
    xs_inf = (1/(1+EXP((7.14285714285714e-02*(-5-Y(1))))))
    alpha_xs = (1400/SQRT((1+EXP((1.66666666666667e-01*(5-Y(1)))))))
    beta_xs = (1/(1+EXP((6.66666666666667e-02*(-35+Y(1))))))
    tau_xs = (80+(1*alpha_xs*beta_xs))
    dydt(16) = ((xs_inf-Y(16))/tau_xs)
    dydt(17) = (-1.16885497444019e-04*(i_Na+i_b_Na+(3*i_NaK)+(3*i_NaCa)))
    r_inf = (1/(1+EXP((1.66666666666667e-01*(20-Y(1))))))
    tau_r = (8.0e-01+(9.5e+00*EXP((5.55555555555556e-04*(-pow((40+Y(1)),2))))))
    dydt(18) = ((r_inf-Y(18))/tau_r)
    s_inf = (1/(1+EXP((2.0e-01*(20+Y(1))))))
    tau_s = (3+(85*EXP((3.125e-03*(-pow((45+Y(1)),2)))))+(5/(1+EXP((2.0e-01*(-20+Y(1)))))))
    dydt(19) = ((s_inf-Y(19))/tau_s)
    
  END SUBROUTINE Electrophysiology_TenTusscher06dYdt

  !
  !================================================================================================================================
  !

  !>Integrate a Ten Tusscher 06 cell model.
  SUBROUTINE Electrophysiology_TenTusscher06Integrate(cellsField,materialsField,t0,t1,err,error,*)
    
    !Argument variables
    TYPE(FieldType), INTENT(INOUT), POINTER :: cellsField !<Independent field storing the cell data
    TYPE(FieldType), INTENT(INOUT), POINTER :: materialsField !<Material field component 1 the activation flags
    REAL(DP), INTENT(IN)    :: t0, t1 !<Integrate from time t0 to t1 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: cellIdx,componentIdx,localDOFIdx,numberOfCells
    INTEGER, PARAMETER :: cellDimension = 19
    REAL(DP) :: dYdt(1:cellDimension)
    REAL(DP), POINTER :: y(:)
    REAL(DP) :: t, dt, activationFactor, m_inf, d_inf, m_inf0, d_inf0, prev_v
    REAL(DP), DIMENSION(:), POINTER :: cellData, activationData
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(FieldType), POINTER :: activationField
    TYPE(FieldVariableType), POINTER :: activationVariable,cellsVariable
    
    ENTERS('Electrophysiology_TenTusscher06Integrate',err,error,*999)
    
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(cellsField,decomposition,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,1,domain,err,error,*999)
    NULLIFY(domainTopology)
    CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(domainNodes)
    CALL DomainTopoloy_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
    CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfCells,err,error,*999)

    NULLIFY(cellsVariable)
    CALL Field_VariableGet(cellsField,FIELD_U_VARIABLE_TYPE,cellsVariable,err,error,*999)
    NULLIFY(activationVariable)
    CALL Field_VariableGet(activationField,FIELD_U_VARIABLE_TYPE,activationVariable,err,error,*999)
    

    CALL FieldVariable_ParameterSetDataGet(cellsVariable,FIELD_VALUES_SET_TYPE,cellData,err,error,*999)
    CALL FieldVariable_ParameterSetDataGet(activationVariable,FIELD_VALUES_SET_TYPE,activationData,err,error,*999)

    DO cellIdx=1,numberOfCells
      !Default to version 1 of each node derivative
      CALL FieldVariable_LocalNodeDOFGet(cellsVariable,1,1,cellIdx,1,localDOFIdx,err,error,*999)
      y => cellData(localDOFIdx:localDOFIdx+cellDimension-1)
      prev_v = y(1)
      !Default to version 1 of each node derivative
      CALL FieldVariable_LocalNodeDOFGet(activationVariable,1,1,cellIdx,1,localDOFIdx,err,error,*999)
      activationFactor = activationData(localDOFIdx)
      
  
      t = t0
      DO WHILE (t < (t1 - ZERO_TOLERANCE))
        ! integrate one cell, one time step. 
        ! just use adaptive forward euler, tested in matlab
        call Electrophysiology_TenTusscher06dYdt(t,y,dYdt,activationFactor)
        dt = MIN(t1-t,1.0/3)  ! default max
        dt = MIN(dt,1 / ABS(dYdt(1))) ! maximum increase in Vm

        m_inf0 = (1/pow((1+EXP((1.10741971207087e-01*(-5.686e+01-Y(1))))),2))
        d_inf0 = (1/(1+EXP((1.33333333333333e-01*(-8-Y(1))))))
        
        y = y + dt * dYdt ! FWE
        
        m_inf = (1/pow((1+EXP((1.10741971207087e-01*(-5.686e+01-Y(1))))),2))
        d_inf = (1/(1+EXP((1.33333333333333e-01*(-8-Y(1))))))
        
        ! FWEoo on m,d gates : prevent gates from crossing midpoint steady state
        
        m_inf = (m_inf0 + m_inf)/2
        d_inf = (d_inf0 + d_inf)/2
        IF(dYdt(11) * (y(11) - m_inf) > 0) THEN
          y(11) = m_inf
        ENDIF
        IF(dYdt(12) * (y(12) - d_inf) > 0) THEN
          y(12) = d_inf
        ENDIF
        
        t = t + dt
      ENDDO
      IF(prev_v < 0 .AND. y(1) > 0) THEN ! store activation times, where else?
        !Default to version 1 of each node derivative
        CALL FieldVariable_LocalNodeDOFGet(activationVariable,1,1,cellIdx,7,localDOFIdx,err,error,*999)
        activationData(localDOFIdx)=t1
      ENDIF
    ENDDO !cellIdx
    
    CALL FieldVariable_ParameterSetDataRestore(cellsVariable,FIELD_VALUES_SET_TYPE,cellData,err,error,*999)
    CALL FieldVariable_ParameterSetDataRestore(activationVariable,FIELD_VALUES_SET_TYPE,activationData,err,error,*999)

    EXITS('Electrophysiology_TenTusscher06Integrate')
    RETURN
999 ERRORSEXITS('Electrophysiology_TenTusscher06Integrate',err,error)
    RETURN 1
    
  END SUBROUTINE Electrophysiology_TenTusscher06Integrate
  
  !
  !================================================================================================================================
  !
  
END MODULE ElectrophysiologyCellRoutines


