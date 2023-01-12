/*
 * \file
 * \author Chris Bradley
 * \brief This is an example program which sets up a field which uses a more complex mesh using OpenCMISS calls from C.
 *
 * \section LICENSE
 *
 * Version: MPL 1.1/GPL 2.0/LGPL 2.1
 *
 * The contents of this file are subject to the Mozilla Public License
 * Version 1.1 (the "License"); you may not use this file except in
 * compliance with the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS"
 * basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
 * License for the specific language governing rights and limitations
 * under the License.
 *
 * The Original Code is OpenCMISS
 *
 * The Initial Developer of the Original Code is University of Auckland,
 * Auckland, New Zealand and University of Oxford, Oxford, United
 * Kingdom. Portions created by the University of Auckland and University
 * of Oxford are Copyright (C) 2007 by the University of Auckland and
 * the University of Oxford. All Rights Reserved.
 *
 * Contributor(s):
 *
 * Alternatively, the contents of this file may be used under the terms of
 * either the GNU General Public License Version 2 or later (the "GPL"), or
 * the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
 * in which case the provisions of the GPL or the LGPL are applicable instead
 * of those above. If you wish to allow use of your version of this file only
 * under the terms of either the GPL or the LGPL, and not to allow others to
 * use your version of this file under the terms of the MPL, indicate your
 * decision by deleting the provisions above and replace them with the notice
 * and other provisions required by the GPL or the LGPL. If you do not delete
 * the provisions above, a recipient may use your version of this file under
 * the terms of any one of the MPL, the GPL or the LGPL.
 *
 */
#include <stdlib.h>
#include <stdio.h>

#include "opencmiss/iron.h"

#define STRING_SIZE 20

#define CONTEXT_USER_NUMBER 1
#define REGION_USER_NUMBER 2

int main() 
{

  cmfe_ContextType context=(cmfe_ContextType)NULL;
  cmfe_RegionType region=(cmfe_RegionType)NULL,worldRegion=(cmfe_RegionType)NULL;
  char label[STRING_SIZE];
  int err;

  err = cmfe_Initialise();
  OPENCMISS_CHECK_ERROR(err,"Initialising OpenCMISS");
  err = cmfe_Context_Initialise(&context);
  OPENCMISS_CHECK_ERROR(err,"Initialising context");
  err = cmfe_Context_Create(CONTEXT_USER_NUMBER,context);
  OPENCMISS_CHECK_ERROR(err,"Creating context");
  err = cmfe_Region_Initialise(&worldRegion);
  OPENCMISS_CHECK_ERROR(err,"Initialising world region");
  err = cmfe_Context_WorldRegionGet(context,worldRegion);
  OPENCMISS_CHECK_ERROR(err,"Get world region");
  
  err = cmfe_Region_LabelGet(worldRegion,STRING_SIZE,label);
  printf("The world region label is '%s'.\n",label);
  
  err = cmfe_Region_Initialise(&region);
  err = cmfe_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region);
  err = cmfe_Region_LabelSet(region,8,"Testing");
  err = cmfe_Region_CreateFinish(region);
  
  err = cmfe_Region_LabelGet(region,STRING_SIZE,label);	       
  printf("The region label is '%s'.\n",label);
  
  /* Destroy the region */
  err = cmfe_Region_Finalise(&region);
  /* Destroy the context */
  err = cmfe_Context_Destroy(context);
  /* Finalise OpenCMISS */
  err = cmfe_Finalise();

  return err;
}
