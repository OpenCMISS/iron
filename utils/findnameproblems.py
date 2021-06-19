#!/usr/bin/env python

"""
Find problems with the diagnostics enters, exits and errors or DLLEXPORT statemtnts in OpenCMISS Fortran source.
This script will read all *.F90 source files in a given source path and check the following
for subroutines and functions:
 - That the names in DLLEXPORT statements match the routine name.
 - That the names in ENTERS, EXITS, ERRORS and ERRORSEXITS macros match the routine name.
 - That each function or subroutine has an ENTERS and an EXITS
 - If there is an EXITS, that the number of RETURN statements matches the number of EXITS statements
 - That RETURN statements do not have line numbers which could then be used to circumvent an EXITS

PURE FUNCTIONS and routines inside INTERFACE/END INTERFACE blocks are not checked. 

A source file can also be omitted from checking by adding the name to the ignore_sources list.

A routine can also be omitted from checking by adding the name to the ignore_routines list.

Debugging can be controlled by setting the debug_level variable [0-5].

Chris Bradley 16/2/2017.
"""

import sys
import os
import re

debug_level = 0

ignore_sources = ["fieldml_input_routines.F90","fieldml_output_routines.F90","fieldml_util_routines.F90",
		"field_IO_routines.F90","Hamilton_Jacobi_equations_routines.F90","electrophysiology_cell_routines.F90"]

ignore_routines = ["Enters","Exits","Errors","ExtractErrorMessageVS","ExtractErrorMessageC",
		"ExtractErrorStackVS","ExtractErrorStackC","FlagErrorC","FlagErrorVS",
		"FlagWarningVS","FlagWarningC","BaseRoutines_Finalise","BaseRoutines_Initialise",
		"LogicalToCharacter","LogicalToVString","NumberToCharacterIntg","NumberToCharacterLIntg",
		"NumberToCharacterSP","NumberToCharacterDP","NumberToVStringIntg","NumberToVStringLIntg",
		"NumberToVStringSP","NumberToVStringDP","StringToDoubleC","StringToDoubleVS",
		"StringToIntegerC","StringToIntegerVS","StringToLongIntegerC","StringToLongIntegerVS",
		"StringToLogicalC","StringToLogicalVS","StringToSingleC","StringToSingleVS",
		"WriteError","WriteStr","SNESSetJacobianBuffer","WriteStringC","WriteStringVS",
		"WriteStringValueC","WriteStringValueVS","WriteStringValueDP","WriteStringValueIntg",
		"WriteStringValueLIntg","WriteStringValueL","WriteStringValueSP","WriteStringValueVS",
		"WriteStringTwoValueCC","WriteStringTwoValueCDP","WriteStringTwoValueCIntg",
		"WriteStringTwoValueCL","WriteStringTwoValueCSP","WriteStringTwoValueCVS",
		"WriteStringTwoValueDPC","WriteStringTwoValueDPDP","WriteStringTwoValueDPIntg",
		"WriteStringTwoValueDPL","WriteStringTwoValueDPSP","WriteStringTwoValueDPVS",
		"WriteStringTwoValueIntgC","WriteStringTwoValueIntgDP","WriteStringTwoValueIntgIntg",
		"WriteStringTwoValueIntgL","WriteStringTwoValueIntgSP","WriteStringTwoValueIntgVS",
		"WriteStringTwoValueLC","WriteStringTwoValueLDP","WriteStringTwoValueLIntg",
		"WriteStringTwoValueLL","WriteStringTwoValueLSP","WriteStringTwoValueLVS",
		"WriteStringTwoValueSPC","WriteStringTwoValueSPDP","WriteStringTwoValueSPIntg",
		"WriteStringTwoValueSPL","WriteStringTwoValueSPSP","WriteStringTwoValueSPVS",
		"WriteStringTwoValueVSC","WriteStringTwoValueVSDP","WriteStringTwoValueVSIntg",
		"WriteStringTwoValueVSL","WriteStringTwoValueVSSP","WriteStringTwoValueVSVS",
		"WriteStringFmtValueC","WriteStringFmtValueDP","WriteStringFmtValueIntg",
		"WriteStringFmtValueLIntg","WriteStringFmtValueL","WriteStringFmtValueSP",
		"WriteStringFmtValueVS","WriteStringFmtTwoValueCC","WriteStringFmtTwoValueCDP",
		"WriteStringFmtTwoValueCIntg","WriteStringFmtTwoValueCL","WriteStringFmtTwoValueCSP",
		"WriteStringFmtTwoValueCVS","WriteStringFmtTwoValueDPC","WriteStringFmtTwoValueDPDP",
		"WriteStringFmtTwoValueDPIntg","WriteStringFmtTwoValueDPL","WriteStringFmtTwoValueDPSP",
		"WriteStringFmtTwoValueDPVS","WriteStringFmtTwoValueIntgC","WriteStringFmtTwoValueIntgDP",
		"WriteStringFmtTwoValueIntgIntg","WriteStringFmtTwoValueIntgL","WriteStringFmtTwoValueIntgSP",
		"WriteStringFmtTwoValueIntgVS","WriteStringFmtTwoValueLC","WriteStringFmtTwoValueLDP",
		"WriteStringFmtTwoValueLIntg","WriteStringFmtTwoValueLL","WriteStringFmtTwoValueLSP",
		"WriteStringFmtTwoValueLVS","WriteStringFmtTwoValueSPC","WriteStringFmtTwoValueSPDP",
		"WriteStringFmtTwoValueSPIntg","WriteStringFmtTwoValueSPL","WriteStringFmtTwoValueSPSP",
		"WriteStringFmtTwoValueSPVS","WriteStringFmtTwoValueVSC","WriteStringFmtTwoValueVSDP",
		"WriteStringFmtTwoValueVSIntg","WriteStringFmtTwoValueVSL","WriteStringFmtTwoValueVSSP",
		"WriteStringFmtTwoValueVSVS","WriteStringVectorDP","WriteStringVectorIntg",
		"WriteStringVectorLIntg","WriteStringVectorL","WriteStringVectorSP","WriteStringIdxVectorDP",
		"WriteStringIdxVectorIntg","WriteStringIdxVectorLIntg","WriteStringIdxVectorL",
		"WriteStringIdxVectorSP","WriteStringMatrixDP","WriteStringMatrixIntg","WriteStringMatrixLIntg",
		"WriteStringMatrixL","WriteStringMatrixSP","cmfe_FinaliseNumber","cmfe_FinaliseObj","cmfe_InitialiseNumber",
		"cmfe_InitialiseObj","cmfe_ExtractErrorMessageC","cmfe_ExtractErrorMessageVS","cmfe_ExtractErrorStackC",
		"cmfe_ExtractErrorStackVS","cmfe_Finalise_","cmfe_Initialise_","cmfe_HandleError","CMISSC2FString",
		"CMISSF2CString","CMISSC2FStrings",
		"clooping","ElasticDil","StateVarsIntegrator","deviatorBepr","UpdateFibersHomeo","nonlinearfuncs",
		"alpha1B","ElasticMeasures","Constitutive","umat","SmoothMultiPhase","ElasticTrialDist","EquivStrain",
		"EffDistStrain","Beprpr","OneThirdAlpha","Beprime","DeltatGamma","CauchyStress","hardening","Tangent",
		"JuxtaTensVec","JuxtaTensTens","TransTens","DetTens","TraceTens","InvTens","DotTensTens","IdentityTens",
		"ZeroTens","DevTens","MUnimodular","CrossVecVec","DotVecVec","TensProd","TensProd33","oplus","ominus",
	        "DotTens2Tens4","DotTens4Tens2","InvTens2by2","MacBrackets",
                "Problem_SolverJacobianEvaluatePetsc","Problem_SolverJacobianFDCalculatePetsc",
                "Problem_SolverObjectiveEvaluatePetsc","Problem_SolverResidualEvaluatePetsc",
                "Problem_SolverConvergenceTestPetsc","Problem_SolverDAECellMLRHSPetsc",
                "Problem_SolverNonlinearMonitorPETSC","Problem_SolverOptimiserMonitorPETSC"]

def _join_lines(source):
    """Remove Fortran line continuations"""

    return re.sub(r'[\t ]*&[\t ]*[\r\n]+[\t ]*&[\t ]*', ' ', source)

if len(sys.argv) > 1:
   source_path = sys.argv[1]
else:
   sys.stderr.write('Usage: ' + sys.argv[0] + ' source_path\n')
   exit(1)

num_fails = 0

iron_source_path = source_path + os.sep + 'src'
source_files = [
	     iron_source_path + os.sep + file_name
             for file_name in os.listdir(iron_source_path)
             if not (file_name in ignore_sources) and file_name.endswith('.F90')]

startinterface_re = re.compile(r'^\s*INTERFACE\s*!*', re.IGNORECASE)
endinterface_re   = re.compile(r'^\s*END\s*INTERFACE\s*!*', re.IGNORECASE)
startroutine_re   = re.compile(r'^\s*(PURE)*\s*(RECURSIVE)*\s*(SUBROUTINE|FUNCTION)+\s+([A-Za-z0-9_]+)\(')
endroutine_re     = re.compile(r'^\s*END\s*(SUBROUTINE|FUNCTION)+', re.IGNORECASE)

dllexport_re   = re.compile(r'^\s*!DLLEXPORT\(([A-Za-z0-9_]+)\)')
enters_re      = re.compile(r'^\s*[0-9]*\s*ENTERS\("([A-Za-z0-9_]+)",')
exits_re       = re.compile(r'^\s*[0-9]*\s*EXITS\("([A-Za-z0-9_]+)"\)')
errors_re      = re.compile(r'^\s*[0-9]*\s*ERRORS\("([A-Za-z0-9_]+)",')
errorsexits_re = re.compile(r'^\s*[0-9]*\s*ERRORSEXITS\("([A-Za-z0-9_]+)",')
return_re      = re.compile(r'^\s*([0-9]*)\s*RETURN\s*["A-Z0-9_,]*', re.IGNORECASE)

for source in source_files:
    print("Processing source file : ",source,"\n")
    source_lines = _join_lines(open(source, 'r').read()).splitlines()
    in_interface = 0
    in_routine = 0
    pure_routine = 0
    routinename = ""
    for (line_number, line) in enumerate(source_lines):
        if debug_level >= 3:
            print("  Line : ",line,"\n")
        if in_interface == 0:
            startinterface_match = startinterface_re.search(line)
            if debug_level >= 4:
                print("      Start interface match    : ",startinterface_match,"\n")
            if startinterface_match:
                in_interface = 1
        else:
            endinterface_match = endinterface_re.search(line)
            if debug_level >= 4:
                print("      End interface match      : ",endinterface_match,"\n")
            if endinterface_match:
                in_interface = 0
        if in_routine == 0:
            startroutine_match = startroutine_re.search(line)
            if debug_level >= 4:
                print("      Start routine match      : ",startroutine_match,"\n")
            if startroutine_match:
                if debug_level >= 5:
                    print("        Subroutine match group 1: ",startroutine_match.group(1),"\n")
                    print("        Subroutine match group 2: ",startroutine_match.group(2),"\n")
                    print("        Subroutine match group 3: ",startroutine_match.group(3),"\n")
                    print("        Subroutine match group 4: ",startroutine_match.group(4),"\n")
                routinename = startroutine_match.group(4)
                if startroutine_match.group(1) != None:
                    pure_routine = 1
                in_routine = 1
                dllexports = []
                enters = []
                exits = []
                errors= []
                errorsexits = []
                num_dllexports = 0
                num_returns = 0
                problem_return = 0
                if debug_level >= 2:
                    print("  Routine : ",routinename,"\n")
                    print("    Start routine ....\n") 
        else:
            endroutine_match = endroutine_re.search(line)
            if debug_level >= 4:
                print("      End routine match        : ",endroutine_match,"\n")
            if endroutine_match:
                if debug_level >=5:
                    print("        in_interface    = ",in_interface,"\n")
                    print("        pure_routine    = ",pure_routine,"\n")
                if routinename in ignore_routines:
                    if debug_level >=5:
                        print("        Routine name is in ignore list\n")
                else:
                    if in_interface == 0:
                        if pure_routine == 0:
                            num_dllexports = len(dllexports)
                            num_enters = len(enters)
                            num_exits = len(exits)
                            num_errors = len(errors)
                            num_errorsexits = len(errorsexits)
                            total_num_exits = num_exits + num_errorsexits
                            if debug_level >=5:
                                print("        num_dllexports  = ",num_dllexports,"\n")
                                print("        num_enters      = ",num_enters,"\n")
                                print("        num_exits       = ",num_exits,"\n")
                                print("        num_errors      = ",num_errors,"\n")
                                print("        num_errorsexits = ",num_errorsexits,"\n")
                                print("        total_num_exits = ",total_num_exits,"\n")
                                print("        num_returns     = ",num_returns,"\n")
                            if num_dllexports > 1:
                                print("ERROR: ",routinename," - ",num_dllexports," DLLEXPORT statements found\n")
                                num_fails = num_fails + 1
                            elif num_dllexports == 1:
                                for dllexportname in dllexports:
                                    if(routinename != dllexportname):
                                        print("ERROR: ",routinename," - DLLEXPORT name mismatch\n")
                                        num_fails = num_fails + 1
                            if num_enters == 0:
                                print("ERROR: ",routinename," - No ENTERS found\n")
                                num_fails = num_fails + 1
                            else:
                                if num_enters > 1:
                                    print("ERROR: ",routinename," - ",num_enters," ENTERS found\n")
                                    num_fails = num_fails + 1
                                    for entername in enters:
                                        if(routinename != entername):
                                            print("ERROR: ",routinename," - ENTERS name mismatch\n")
                                            num_fails = num_fails + 1
                            if total_num_exits == 0:
                                print("ERROR: ",routinename," - No EXITS found\n")
                                num_fails = num_fails + 1
                            else:
                                for exitname in exits:
                                    if(routinename != exitname):
                                        print("ERROR: ",routinename," - EXITS name mismatch\n")
                                        num_fails = num_fails + 1
                            for errorsname in errors:
                                if(routinename != errorsname):
                                    print("ERROR: ",routinename," - ERRORS name mismatch\n")
                                    num_fails = num_fails + 1
                            for errorsexitsname in errorsexits:
                                if(routinename != errorsexitsname):
                                    print("ERROR: ",routinename," - ERRORSEXITS name mismatch\n")
                                    num_fails = num_fails + 1
                            if total_num_exits > 0:
                                if total_num_exits != num_returns:
                                    print("ERROR: ",routinename," - Number of returns and exits do not match\n")
                                    num_fails = num_fails + 1
                            if problem_return != 0:
                                print("ERROR: ",routinename," - Possible EXITS circumvention\n")
                                num_fails = num_fails + 1
                if debug_level >= 2:
                    print("    Finish routine\n") 
                in_routine = 0
                pure_routine = 0
            else:
                dllexport_match = dllexport_re.search(line)
                enters_match = enters_re.search(line)
                exits_match = exits_re.search(line)
                errors_match = errors_re.search(line)
                errorsexits_match = errorsexits_re.search(line)
                return_match = return_re.search(line)
                if debug_level >= 4:
                    print("      DLLEXPORT routine match  : ",dllexport_match,"\n")
                    print("      Enters routine match     : ",enters_match,"\n")
                    print("      Exits routine match      : ",exits_match,"\n")
                    print("      Errors routine match     : ",errors_match,"\n")
                    print("      Errorsexits routine match: ",errorsexits_match,"\n")
                    print("      Return routine match     : ",return_match,"\n")
                if dllexport_match:
                    if debug_level >= 5:
                        print("        DLLEXPORT match group 1: ", dllexport_match.group(1),"\n")
                    dllexports.append(dllexport_match.group(1))
                    num_dllexports = num_dllexports + 1
                elif enters_match:
                    if debug_level >= 5:
                        print("        Enters match group 1: ", enters_match.group(1),"\n")
                    enters.append(enters_match.group(1))
                elif exits_match:
                    if debug_level >= 5:
                        print("        Exits match group 1: ",exits_match.group(1),"\n")
                    exits.append(exits_match.group(1))
                elif errors_match:
                    if debug_level >= 5:
                        print("        Errors match group 1: ",errors_match.group(1),"\n")
                    errors.append(errors_match.group(1))
                elif errorsexits_match:
                    if debug_level >= 5:
                        print("        Errorsexits match group 1: ",errorsexits_match.group(1),"\n")
                    errorsexits.append(errorsexits_match.group(1))
                elif return_match:
                    if debug_level >= 5:
                        print("        Return match group 1: ",return_match.group(1),"\n")
                    if len(return_match.group(1)) > 0:
                        problem_return = 1
                    num_returns = num_returns + 1
if num_fails == 0:
    print("NO ERRORS FOUND - PASS\n")
    exit(0)
else:
    print(num_fails," ERRORS FOUND - FAIL\n")
    exit(1)
 
