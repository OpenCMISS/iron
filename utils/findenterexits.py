#!/usr/bin/env python

"""
Find problems with the diagnostics enters, exits and errors in OpenCMISS Fortran source.
This script will read all *.f90 source files in a given source path and check the following
for subroutines and functions:
 - That the names in ENTERS, EXITS, ERRORS and ERRORS exits macros match the routine name.
 - That each function or subroutine has an ENTERS and an EXITS
 - If there is an EXITS, that the number of RETURN statements matches the number of EXITS statements
 - That RETURN statements do not have line numbers which could then be used to circumvent an EXITS

PURE FUNCTIONS and routines inside INTERFACE/END INTERFACE blocks are not checked. 

A routine can also be omitted from checking by adding the name to the ignore_routines list.

Debugging can be controlled by setting the debug_level variable [0-5].

Chris Bradley 16/2/17.
"""

import sys
import os
import re

debug_level = 0

ignore_routines = ["ENTERS","EXITS","ERRORS","EXTRACT_ERROR_MESSAGE_VS","EXTRACT_ERROR_MESSAGE_C",
		"ExtractErrorStackVS","ExtractErrorStackC","FLAG_ERROR_C","FLAG_ERROR_VS",
		"FLAG_WARNING_VS","FLAG_WARNING_C","BASE_ROUTINES_FINALISE","BASE_ROUTINES_INITIALISE",
		"WriteError","WRITE_STR","SNESSetJacobianBuffer","WRITE_STRING_C","WRITE_STRING_VS",
		"WRITE_STRING_VALUE_C","WRITE_STRING_VALUE_VS","WRITE_STRING_VALUE_DP","WRITE_STRING_VALUE_INTG",
		"WRITE_STRING_VALUE_LINTG","WRITE_STRING_VALUE_L","WRITE_STRING_VALUE_SP","WRITE_STRING_VALUE_VS",
		"WRITE_STRING_TWO_VALUE_C_C","WRITE_STRING_TWO_VALUE_C_DP","WRITE_STRING_TWO_VALUE_C_INTG",
		"WRITE_STRING_TWO_VALUE_C_L","WRITE_STRING_TWO_VALUE_C_SP","WRITE_STRING_TWO_VALUE_C_VS",
		"WRITE_STRING_TWO_VALUE_DP_C","WRITE_STRING_TWO_VALUE_DP_DP","WRITE_STRING_TWO_VALUE_DP_INTG",
		"WRITE_STRING_TWO_VALUE_DP_L","WRITE_STRING_TWO_VALUE_DP_SP","WRITE_STRING_TWO_VALUE_DP_VS",
		"WRITE_STRING_TWO_VALUE_INTG_C","WRITE_STRING_TWO_VALUE_INTG_DP","WRITE_STRING_TWO_VALUE_INTG_INTG",
		"WRITE_STRING_TWO_VALUE_INTG_L","WRITE_STRING_TWO_VALUE_INTG_SP","WRITE_STRING_TWO_VALUE_INTG_VS",
		"WRITE_STRING_TWO_VALUE_L_C","WRITE_STRING_TWO_VALUE_L_DP","WRITE_STRING_TWO_VALUE_L_INTG",
		"WRITE_STRING_TWO_VALUE_L_L","WRITE_STRING_TWO_VALUE_L_SP","WRITE_STRING_TWO_VALUE_L_VS",
		"WRITE_STRING_TWO_VALUE_SP_C","WRITE_STRING_TWO_VALUE_SP_DP","WRITE_STRING_TWO_VALUE_SP_INTG",
		"WRITE_STRING_TWO_VALUE_SP_L","WRITE_STRING_TWO_VALUE_SP_SP","WRITE_STRING_TWO_VALUE_SP_VS",
		"WRITE_STRING_TWO_VALUE_VS_C","WRITE_STRING_TWO_VALUE_VS_DP","WRITE_STRING_TWO_VALUE_VS_INTG",
		"WRITE_STRING_TWO_VALUE_VS_L","WRITE_STRING_TWO_VALUE_VS_SP","WRITE_STRING_TWO_VALUE_VS_VS",
		"WRITE_STRING_FMT_VALUE_C","WRITE_STRING_FMT_VALUE_DP","WRITE_STRING_FMT_VALUE_INTG",
		"WRITE_STRING_FMT_VALUE_LINTG","WRITE_STRING_FMT_VALUE_L","WRITE_STRING_FMT_VALUE_SP",
		"WRITE_STRING_FMT_VALUE_VS","WRITE_STRING_FMT_TWO_VALUE_C_C","WRITE_STRING_FMT_TWO_VALUE_C_DP",
		"WRITE_STRING_FMT_TWO_VALUE_C_INTG","WRITE_STRING_FMT_TWO_VALUE_C_L","WRITE_STRING_FMT_TWO_VALUE_C_SP",
		"WRITE_STRING_FMT_TWO_VALUE_C_VS","WRITE_STRING_FMT_TWO_VALUE_DP_C","WRITE_STRING_FMT_TWO_VALUE_DP_DP",
		"WRITE_STRING_FMT_TWO_VALUE_DP_INTG","WRITE_STRING_FMT_TWO_VALUE_DP_L","WRITE_STRING_FMT_TWO_VALUE_DP_SP",
		"WRITE_STRING_FMT_TWO_VALUE_DP_VS","WRITE_STRING_FMT_TWO_VALUE_INTG_C","WRITE_STRING_FMT_TWO_VALUE_INTG_DP",
		"WRITE_STRING_FMT_TWO_VALUE_INTG_INTG","WRITE_STRING_FMT_TWO_VALUE_INTG_L","WRITE_STRING_FMT_TWO_VALUE_INTG_SP",
		"WRITE_STRING_FMT_TWO_VALUE_INTG_VS","WRITE_STRING_FMT_TWO_VALUE_L_C","WRITE_STRING_FMT_TWO_VALUE_L_DP",
		"WRITE_STRING_FMT_TWO_VALUE_L_INTG","WRITE_STRING_FMT_TWO_VALUE_L_L","WRITE_STRING_FMT_TWO_VALUE_L_SP",
		"WRITE_STRING_FMT_TWO_VALUE_L_VS","WRITE_STRING_FMT_TWO_VALUE_SP_C","WRITE_STRING_FMT_TWO_VALUE_SP_DP",
		"WRITE_STRING_FMT_TWO_VALUE_SP_INTG","WRITE_STRING_FMT_TWO_VALUE_SP_L","WRITE_STRING_FMT_TWO_VALUE_SP_SP",
		"WRITE_STRING_FMT_TWO_VALUE_SP_VS","WRITE_STRING_FMT_TWO_VALUE_VS_C","WRITE_STRING_FMT_TWO_VALUE_VS_DP",
		"WRITE_STRING_FMT_TWO_VALUE_VS_INTG","WRITE_STRING_FMT_TWO_VALUE_VS_L","WRITE_STRING_FMT_TWO_VALUE_VS_SP",
		"WRITE_STRING_FMT_TWO_VALUE_VS_VS","WRITE_STRING_VECTOR_DP","WRITE_STRING_VECTOR_INTG",
		"WRITE_STRING_VECTOR_LINTG","WRITE_STRING_VECTOR_L","WRITE_STRING_VECTOR_SP","WRITE_STRING_IDX_VECTOR_DP",
		"WRITE_STRING_IDX_VECTOR_INTG","WRITE_STRING_IDX_VECTOR_LINTG","WRITE_STRING_IDX_VECTOR_L",
		"WRITE_STRING_IDX_VECTOR_SP","WRITE_STRING_MATRIX_DP","WRITE_STRING_MATRIX_INTG","WRITE_STRING_MATRIX_LINTG",
		"WRITE_STRING_MATRIX_L","WRITE_STRING_MATRIX_SP","cmfe_Finalise","cmfe_InitialiseNumber",
		"cmfe_InitialiseObj","cmfe_ExtractErrorMessageC","cmfe_ExtractErrorMessageVS","cmfe_ExtractErrorStackC",
		"cmfe_ExtractErrorStackVS","cmfe_Finalise_","cmfe_Initialise_","cmfe_HandleError","CMISSC2FString",
		"CMISSF2CString","CMISSC2FStrings"]

def _join_lines(source):
    """Remove Fortran line continuations"""

    return re.sub(r'[\t ]*&[\t ]*[\r\n]+[\t ]*&[\t ]*', ' ', source)

if len(sys.argv) >= 1:
   source_path = sys.argv[1]
else:
   sys.stderr.write('Usage: %s source_path\n' % sys.argv[0])
   exit(1)

iron_source_path = source_path + os.sep + 'src'
source_files = [
	     iron_source_path + os.sep + file_name
             for file_name in os.listdir(iron_source_path)
             if file_name.endswith('.f90')]

startinterface_re = re.compile(r'^\s*INTERFACE\s*!*', re.IGNORECASE)
endinterface_re = re.compile(r'^\s*END\s*INTERFACE\s*!*', re.IGNORECASE)
startroutine_re = re.compile(r'^\s*(PURE)*\s*(RECURSIVE)*\s*(SUBROUTINE|FUNCTION)+\s+([A-Za-z0-9_]+)\(')
endroutine_re = re.compile(r'^\s*END\s*(SUBROUTINE|FUNCTION)+', re.IGNORECASE)

enters_re      = re.compile(r'^\s*[0-9]*\s*ENTERS\("([A-Za-z0-9_]+)",')
exits_re       = re.compile(r'^\s*[0-9]*\s*EXITS\("([A-Za-z0-9_]+)"\)')
errors_re      = re.compile(r'^\s*[0-9]*\s*ERRORS\("([A-Za-z0-9_]+)",')
errorsexits_re = re.compile(r'^\s*[0-9]*\s*ERRORSEXITS\("([A-Za-z0-9_]+)",')
return_re = re.compile(r'^\s*([0-9]*)\s*RETURN\s*["A-Z0-9_,]*', re.IGNORECASE)

for source in source_files:
    print("Processing source file : %s\n") %source
    source_lines = _join_lines(open(source, 'r').read()).splitlines()
    in_interface = 0
    in_routine = 0
    pure_routine = 0
    routinename = ""
    for (line_number, line) in enumerate(source_lines):
        if debug_level >= 3:
	    print("  Line : %s\n") % line
	if in_interface == 0:
           startinterface_match = startinterface_re.search(line)
	   if debug_level >= 4:
	      print("      Start interface match    : %s\n") % startinterface_match
	   if startinterface_match:
	      in_interface = 1
	else:
           endinterface_match = endinterface_re.search(line)
	   if debug_level >= 4:
	      print("      End interface match      : %s\n") % endinterface_match
	   if endinterface_match:
	      in_interface = 0
    	if in_routine == 0:
           startroutine_match = startroutine_re.search(line)
	   if debug_level >= 4:
	      print("      Start routine match      : %s\n") % startroutine_match
	   if startroutine_match:
	      if debug_level >= 5:
	      	 print("        Subroutine match group 1: %s\n") % startroutine_match.group(1)
	      	 print("        Subroutine match group 2: %s\n") % startroutine_match.group(2)
	      	 print("        Subroutine match group 3: %s\n") % startroutine_match.group(3)
	      	 print("        Subroutine match group 4: %s\n") % startroutine_match.group(4)
	      routinename = startroutine_match.group(4)
	      if startroutine_match.group(1) != None:
	          pure_routine = 1
	      in_routine = 1
	      enters = []
	      exits = []
	      errors = []
	      errorsexits = []
	      num_returns = 0
	      problem_return = 0
              if debug_level >= 2:
	          print("  Routine : %s\n") % routinename
	          print("    Start routine ....\n") 
	else:
           endroutine_match = endroutine_re.search(line)
	   if debug_level >= 4:
	      print("      End routine match        : %s\n") % endroutine_match
	   if endroutine_match:
	      if debug_level >=5:
	           print("        in_interface    = %d\n") % in_interface
	           print("        pure_routine    = %d\n") % pure_routine
	      if routinename in ignore_routines:
	          if debug_level >=5:
		       print("        Routine name is in ignore list\n")
              else:
	          if in_interface == 0:
		      if pure_routine == 0:
	      	          num_enters = len(enters)
	      	          num_exits = len(exits)
	      	          num_errors = len(errors)
	      	          num_errorsexits = len(errorsexits)
	      	          total_num_exits = num_exits + num_errorsexits
	      		  if debug_level >=5:
	           	      print("        num_enters      = %d\n") % num_enters
	           	      print("        num_exits       = %d\n") % num_exits
	           	      print("        num_errors      = %d\n") % num_errors
	           	      print("        num_errorsexits = %d\n") % num_errorsexits
	           	      print("        total_num_exits = %d\n") % total_num_exits
	           	      print("        num_returns     = %d\n") % num_returns
	      	          if num_enters == 0:
   	      	              print("ERROR: %s - No ENTERS found\n") % routinename
	      	          else:
		              if num_enters > 1:
   	      	                  print("ERROR: %s - %d ENTERS found\n") % routinename,num_enters
		                  for entername in enters:
	      	       	              if(routinename != entername):
   	                                  print("ERROR: %s - ENTERS name mismatch\n") % routinename
	                  if total_num_exits == 0:
   	      	              print("ERROR: %s - No EXITS found\n") % routinename
		          else:
		              for exitname in exits:
	      	                  if(routinename != exitname):
   	                              print("ERROR: %s - EXITS name mismatch\n") % routinename
	                  for errorsname in errors:
	      	              if(routinename != errorsname):
   	                          print("ERROR: %s - ERRORS name mismatch\n") % routinename 
 	                  for errorsexitsname in errorsexits:
	      	              if(routinename != errorsexitsname):
   	                          print("ERROR: %s - ERRORSEXITS name mismatch\n") % routinename
	                  if total_num_exits > 0:
                              if total_num_exits != num_returns:
   	      	      	          print("ERROR: %s - Number of returns and exits do not match\n") % routinename	    		 
	      	          if problem_return != 0:
   	      	              print("ERROR: %s - Possible EXITS circumvention\n") % routinename
	      if debug_level >= 2:
                  print("    Finish routine\n") 
	      in_routine = 0
              pure_routine = 0
	   else:
	      enters_match = enters_re.search(line)
              exits_match = exits_re.search(line)
              errors_match = errors_re.search(line)
              errorsexits_match = errorsexits_re.search(line)
              return_match = return_re.search(line)
	      if debug_level >= 4:
	          print("      Enters routine match     : %s\n") % enters_match
	          print("      Exits routine match      : %s\n") % exits_match
	          print("      Errors routine match     : %s\n") % errors_match
	          print("      Errorsexits routine match: %s\n") % errorsexits_match
	          print("      Return routine match     : %s\n") % return_match
              if enters_match:
	          if debug_level >= 5:
	      	      print("        Enters match group 1: %s\n") % enters_match.group(1)
	      	  enters.append(enters_match.group(1))
              elif exits_match:
	          if debug_level >= 5:
	      	      print("        Exits match group 1: %s\n") % exits_match.group(1)
	       	  exits.append(exits_match.group(1))
              elif errors_match:
	          if debug_level >= 5:
	      	      print("        Errors match group 1: %s\n") % errors_match.group(1)
	       	  errors.append(errors_match.group(1))
              elif errorsexits_match:
	          if debug_level >= 5:
	      	      print("        Errorsexits match group 1: %s\n") % errorsexits_match.group(1)
	       	  errorsexits.append(errorsexits_match.group(1))
              elif return_match:
	          if debug_level >= 5:
	      	      print("        Return match group 1: %s\n") % return_match.group(1)
	      	  if len(return_match.group(1)) > 0:
		      problem_return = 1
	       	  num_returns = num_returns + 1
    #exit()
 
