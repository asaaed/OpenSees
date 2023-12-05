/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.11 $
// $Date: 2010-02-04 00:27:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/tcl/myCommands.cpp,v $

// Written: fmk 
// Created: 04/98
//
// Description: This file contains the function myCommands().
// myCommands() is called in g3AppInit() - all new user commands
// are to be placed in here.
//
// What: "@(#) myCommands.C, revA"

#include <Domain.h>
#include "TclModelBuilder.h"
#include "TclUniaxialMaterialTester.h"
#include "TclPlaneStressMaterialTester.h"
#include "TclSectionTester.h"
#include "stepEnlargement.h" /* stepEnlargement added by Aram */

#include <tcl.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern ModelBuilder *theBuilder;

#ifdef _PARALLEL_PROCESSING
#include <PartitionedDomain.h>
extern PartitionedDomain theDomain;
#else
extern Domain theDomain;
#endif

int
specifyModelBuilder(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

extern int
OPS_ResetInput(ClientData clientData, 
	       Tcl_Interp *interp,  
	       int cArg, 
	       int mArg, 
	       TCL_Char **argv, 
	       Domain *domain,
	       TclModelBuilder *builder);

// reduce command added here
#include "stepEnlargement.h"
static int StepEnlargement_Cmd(
	ClientData clientData,	/* Main window for application. */
	Tcl_Interp *interp,		/* Current interpreter. */
	int objc,			    /* Number of arguments. */
	Tcl_Obj *const objv[])	/* Argument values. */
{
	Tcl_Obj *resultList;

	int nFactor = 0;
	int enlargementType = 0;
	const char *fileNameInput;
	vector<double> datar;
	int length = 0;
	// check if input arg are sufficcient
	if (objc < 4) {
		Tcl_WrongNumArgs(interp, 1, objv, "$enlargementFactor $filePath");
		return TCL_ERROR;
	}
	else {
		if (Tcl_GetIntFromObj(interp, objv[2], &nFactor) != TCL_OK) {
			Tcl_WrongNumArgs(interp, 1, objv, ">> enlargementFactor should be \'Integer\'");
			return TCL_ERROR;
		}
		else if (Tcl_GetIntFromObj(interp, objv[1], &enlargementType) != TCL_OK) {
			Tcl_WrongNumArgs(interp, 1, objv, ">> enlargementType should be \'Integer\'");
			return TCL_ERROR;
		}
		else {
			Tcl_GetIntFromObj(interp, objv[1], &enlargementType);
			Tcl_GetIntFromObj(interp, objv[2], &nFactor);
		}
		fileNameInput = Tcl_GetStringFromObj(objv[3], &length);
	}
	datar = readPath(fileNameInput);
	if (enlargementType == 1) {
		datar = stepEnlargement_INT(nFactor, datar);
	}
	resultList = Tcl_NewListObj(0, NULL);
	for (double &points : datar) {
		Tcl_ListObjAppendElement(interp, resultList, Tcl_NewDoubleObj(points));
	}
	Tcl_SetObjResult(interp, resultList);
	// if writeToFile name is available:
	if (objc == 5) {
		const char *fileNameOutput = Tcl_GetStringFromObj(objv[4], &length);
		writetofile(fileNameOutput, datar);
	}
	
	return TCL_OK;
}
// implement a technique method
vector<double> stepEnlargement_INT(int nFactor, vector<double>& datav)
{
	int nFactorPrime{ 0 };
	int nPrime{ 0 };
	double sigmaf{ 0.0 };
	size_t nStepOrg = datav.size();
	//    check eq 35 of article
	if (nFactor % 2 == 0)
	{
		nFactorPrime = nFactor / 2;
	}
	else
	{
		nFactorPrime = (nFactor - 1) / 2;
	}
	// calculate number of steps of new series
	double dStepNew = ((double)nStepOrg - 1.0) / (double)nFactor + 1.0;
	size_t nStepNew = ceil(dStepNew);
	int nStepOrgNew = (nStepNew - 1)*nFactor + 1;
	while (nStepOrgNew > datav.size()) {
		datav.push_back(0.0);
	}
	// create vector
	vector<double> new_datav(nStepNew, 0);
	// assign f0 to f'0
	new_datav[0] = datav[0];
	// assign fend to f'end
	if (nStepOrgNew == nStepOrg)
	{
		new_datav[nStepNew - 1] = datav[nStepOrg - 1];
	}
	else
	{
		new_datav[nStepNew - 1] = 0.000;
	}

	int iLocNew{ 0 };
	int iLocOrg{ 0 };

	for (int i = 1; i<nStepNew - 1; i++)
	{
		iLocNew = i;
		// calculate Step_i of original series on new series
		iLocOrg = (iLocNew - 1)*nFactor + nFactor;
		// check step 2 and end-1
		if ((i == 1) | (i == nStepNew - 2))
		{
			nPrime = (nFactor - 1);
		}
		else
		{
			nPrime = nFactorPrime;
		}
		sigmaf = 0;
		for (int inneri = 1; inneri <= nPrime; inneri++)
		{
			sigmaf = sigmaf + (1.0 / (4.0 * nPrime))*(datav[iLocOrg + inneri] + datav[iLocOrg - inneri]);
		}
		new_datav[i] = (1.0 / 2.0)*datav[iLocOrg] + sigmaf;
	}
	return new_datav;
}

vector<double> stepEnlargement_MY(int nFactor, vector<double>& datav)
{
	size_t nStepOrg = datav.size();
	// calculate number of steps of new series
	double dStepNew = ((double)nStepOrg - 1.0) / (double)nFactor + 1.0;
	size_t nStepNew = ceil(dStepNew);
	int nStepOrgNew = (nStepNew - 1)*nFactor + 1;
	while (nStepOrgNew > datav.size()) {
		datav.push_back(0.0);
	}
	// create vector
	vector<double> new_datav(nStepNew, 0);

	double sigmaf{ 0.0 };
	double sigmaN{ 0.0 };

	// assign f0 to f'0
	new_datav[0] = datav[0];
	// assign fend to f'end
	// assign fend to f'end
	if (nStepOrgNew == nStepOrg)
	{
		new_datav[nStepNew - 1] = datav[nStepOrg - 1];
	}
	else
	{
		new_datav[nStepNew - 1] = 0.000;;
	}
	int iLocNew{ 0 };
	int iLocOrg{ 0 };

	for (int i = 1; i<nStepNew - 1; i++)
	{
		iLocNew = i;
		// calculate Step_i of original series on new series
		iLocOrg = (iLocNew - 1)*nFactor + nFactor;
		sigmaf = 0;
		sigmaN = 0;
		for (int inneri = 1; inneri < nFactor; inneri++)
		{
			sigmaf = sigmaf + pow(2.0, (nFactor - inneri - 1))*(datav[iLocOrg + inneri] + datav[iLocOrg - inneri]);
			sigmaN = sigmaN + 2 * (pow(2.0, (nFactor - inneri - 1)));
		}
		new_datav[i] = (pow(2.0, (nFactor - 1))*datav[iLocOrg] + sigmaf) / (sigmaN + pow(2.0, (nFactor - 1)));
	}

	return new_datav;
}

// read data by path
vector<double> readPath(string fileNameInput)
{
	double datapoint{ 0.0 };
	int counterData{ 0 };
	//    opens file for reding
	fstream inFile;
	inFile.open(fileNameInput, ios::in);
	//    count the data point
	while (inFile >> datapoint)
	{
		counterData++;
	}
	//    calculate number of steps of original series
	int nStepOrg = counterData;

	vector<double> datav(nStepOrg, 0.0);
	// rewind text file to start point
	inFile.clear();
	inFile.seekg(0);
	// read data
	for (int i = 0; i < nStepOrg; i++)
	{
		inFile >> datav[i];
	}
	inFile.close();
	return datav;
}

// write data by filename
void writetofile(string fileNameOutput, vector<double>& datav)
{

	fstream outFile;
	outFile.open(fileNameOutput, ios::out);
	for (int i = 0; i < datav.size(); i++)
	{
		outFile << fixed << scientific << setw(28) << setprecision(12) << datav[i] << endl;
	}
	outFile.close();
}
int myCommands(Tcl_Interp *interp) {
    Tcl_CreateCommand(interp, "model", specifyModelBuilder,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

	//opserr << " |-----------------------------------------------------------------------------|\n";
	//opserr << " | This version of OpenSees is a custom compiled version, which apply a method |\n";
	//opserr << " |      proposed by [1] to reduce the run time of time history analysis.       |\n";
	//opserr << " |                                                                             |\n";
	//opserr << " |                 Reference: [1] https://doi.org/10.1002/cnm.1097             |\n";
	//opserr << " |               Code Developed by: Aram Saaed (a.saaed@iiees.ac.ir)           |\n";
	//opserr << " |-----------------------------------------------------------------------------|\n";
//	opserr << " | Original command:                                                           |\n";
	//opserr << " | Command (https://opensees.berkeley.edu/wiki/index.php/Path_TimeSeries):     |\n";
//	opserr << " | https://opensees.berkeley.edu/wiki/index.php/Path_TimeSeries                |\n";
//	opserr << " |                                                                             |\n";
	//opserr << " | timeSeries Path $tsTag -dt $dt -values {list_of_values} ...                 |\n";
//	opserr << " |  <-useLast> <-prependZero> <-startTime $tStart>                             |\n";
	//opserr << " |                                                                             |\n";
	//opserr << " | New use of command:                                                         |\n";
//	opserr << " |                                                                             |\n";
	//opserr << " | timeSeries Path $tsTag -dt $dt                                              |\n";
//	opserr << " | <-factor $cFactor> <-useLast> <-prependZero> <-startTime $tStart>           |\n";
	//opserr << " | -values [stepEnlargement $Type $enlargementFactor $filePath <$outFilePath>] |\n";
	//opserr << " |                                                                             |\n";
	//opserr << " | Parameters definition:                                                      |\n";
	//opserr << " |    $Type             : type of reduction ( = 1)                             |\n";
	//opserr << " |    $enlargementFactor: ratio of step size to original step size (integer)   |\n";
	//opserr << " |    $filePath         : file containing the load factors values              |\n";
	//opserr << " |    $outFilePath      : file name, to save new generated record (optional)   |\n";
	//opserr << " |-----------------------------------------------------------------------------|\n";
	//opserr << "\n\n";
    // https://stackoverflow.com/questions/4053837/colorizing-text-in-the-console-with-c
	//HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
	opserr << " |-----------------------------------------------------------------------------|\n";
	//SetConsoleTextAttribute(hConsole, 12);
	opserr << "   This version of OpenSees is a custom compiled version, which apply a method  \n";
	opserr << "        proposed by [1] to reduce the run time of time history analysis.        \n";
	opserr << "                                                                                \n";
	opserr << "                  Reference: [1] https://doi.org/10.1002/cnm.1097               \n";
	opserr << "                 Code Developed by: Aram Saaed (a.saaed@iiees.ac.ir)            \n\n";
	//SetConsoleTextAttribute(hConsole, 7);
	opserr << "   Original Command: \n";
	opserr << "   (https://opensees.berkeley.edu/wiki/index.php/Path_TimeSeries): \n";
	opserr << "   timeSeries Path $tsTag -dt $dt -values ";
	//SetConsoleTextAttribute(hConsole, 14);
	opserr << "{list_of_values}";
	//SetConsoleTextAttribute(hConsole, 7);
	opserr << "<-useLast>  \n";
	opserr << "   <-prependZero> <-startTime $tStart>\n";
	opserr << "                                                                                \n";
	
	opserr << "   New use of command:                                                          \n";
	opserr << "   timeSeries Path $tsTag -dt $dt  \n";
	//opserr << "   -values [stepEnlargement $Type $enlargementFactor $filePath <$outFilePath>]  \n";
	opserr << "   -values ";
	//SetConsoleTextAttribute(hConsole, 14);
	opserr << "[stepEnlargement $Type $enlargementFactor $filePath <$outFilePath>]  \n";
	//SetConsoleTextAttribute(hConsole, 7);
	//opserr << "   -values [stepEnlargement $Type $enlargementFactor $filePath <$outFilePath>]  \n";
	opserr << "   <-factor $cFactor> <-useLast> <-prependZero> <-startTime $tStart> \n";
	opserr << "                                                                                \n";
	opserr << "   Parameters definition:                                                       \n";
	opserr << "      $Type             : type of reduction ( Type = 1)                         \n";
	opserr << "      $enlargementFactor: ratio of step size to original step size (integer)    \n";
	opserr << "      $filePath         : file containing the load factors values               \n";
	opserr << "      $outFilePath      : file name, to save new generated record (optional)    \n";
	
	opserr << " |-----------------------------------------------------------------------------|\n";
	opserr << "\n\n";
	Tcl_CreateObjCommand(interp, "stepEnlargement", &StepEnlargement_Cmd, NULL, NULL);
    return 0;
}


int
specifyModelBuilder(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  int cArg = 0;

  // make sure at least one other argument to contain model builder type given
  if (argc < 2) {
    opserr << "WARNING need to specify a model type, valid types:\n";
    opserr << "\tBasicBuilder\n";
    return TCL_ERROR;
  }    
  
  // invoke the descructor on the old builder
  if (theBuilder != 0) {
    delete theBuilder;
    theBuilder = 0;
  }
  
  // check argv[1] for type of ModelBuilder and create the object 
  if ((strcmp(argv[1],"basic") == 0) || (strcmp(argv[1],"BasicBuilder") == 0) ||
      (strcmp(argv[1],"Basic") == 0) || (strcmp(argv[1],"basicBuilder") == 0)) {
    int ndm =0;
    int ndf = 0;
    
    if (argc < 4) {
      opserr << "WARNING incorrect number of command arguments\n";
      opserr << "model modelBuilderType -ndm ndm? <-ndf ndf?> \n";
      return TCL_ERROR;
    }
    
    int argPos = 2;
    while (argPos < argc) {
      if (strcmp(argv[argPos],"-ndm") == 0 ||
	  strcmp(argv[argPos],"-NDM") == 0) {	
	argPos++;
	if (argPos < argc)
	  if (Tcl_GetInt(interp, argv[argPos], &ndm) != TCL_OK) {
	    opserr << "WARNING error reading ndm: " << argv[argPos];
	    opserr << "\nmodel modelBuilderType -ndm ndm? <-ndf ndf?>\n";
	    return TCL_ERROR;
	  }	  
	argPos++;
      }
      
      else if (strcmp(argv[argPos],"-ndf") == 0 ||
	       strcmp(argv[argPos],"-NDF") == 0) {	
	argPos++;
	if (argPos < argc)
	  if (Tcl_GetInt(interp, argv[argPos], &ndf) != TCL_OK) {
	    opserr << "WARNING error reading ndf: " << argv[argPos];
	    opserr << "\nmodel modelBuilderType -ndm ndm? <-ndf ndf?>\n";
	      return TCL_ERROR;
	  }	  
	argPos++;
      }
      
      else // Advance to next input argument if there are no matches -- MHS
	argPos++;
    }
    
    // check that ndm was specified
    if (ndm == 0) {
      opserr << "WARNING need to specify ndm\n";
      opserr << "model modelBuilderType -ndm ndm? <-ndf ndf?>\n";
      return TCL_ERROR;
    }
    
    // check for ndf, if not assume one
    if (ndf == 0) {
      if (ndm == 1) 
	ndf = 1;
      else if (ndm == 2)
	ndf = 3;
      else if (ndm == 3)
	ndf = 6;
      else {
	opserr << "WARNING specified ndm, " << ndm << ", will not work\n";
	opserr << "with any elements in BasicBuilder\n";
	return TCL_ERROR;
      }
    }
    
    // create the model builder
    TclModelBuilder *theTclBuilder = new TclModelBuilder(theDomain, interp, ndm, ndf);
    
    if (theTclBuilder != 0) {
      OPS_ResetInput(clientData, interp, 0, argc, argv, &theDomain, theTclBuilder);
      theBuilder = theTclBuilder;
    }

    if (theBuilder == 0) {
      opserr << "WARNING ran out of memory in creating BasicBuilder model\n";
      return TCL_ERROR;
    }
  }
  
  else if ((strcmp(argv[1],"test") == 0) || (strcmp(argv[1],"TestUniaxial") == 0) ||
	   (strcmp(argv[1],"testUniaxial") == 0) || (strcmp(argv[1],"UniaxialMaterialTest") == 0)) {
    int count = 1;
    if (argc == 3) {
      if (Tcl_GetInt(interp, argv[2], &count) != TCL_OK) {
	return TCL_ERROR;
      }	  
    }
    theBuilder = new TclUniaxialMaterialTester(theDomain, interp, count);
    if (theBuilder == 0) {
      opserr << "WARNING ran out of memory in creating TclUniaxialMAterialTester model\n";
      return TCL_ERROR;
    }
  }

  else if ((strcmp(argv[1],"testPlaneStress") == 0) || 
	   (strcmp(argv[1],"PlaneStressMaterialTest") == 0)) {
    int count = 1;
    if (argc == 3) {
      if (Tcl_GetInt(interp, argv[2], &count) != TCL_OK) {
	return TCL_ERROR;
      }	  
    }
    
    theBuilder = new TclPlaneStressMaterialTester(theDomain, interp, count);
    if (theBuilder == 0) {
      opserr << "WARNING ran out of memory in creating TclUniaxialMAterialTester model\n";
      return TCL_ERROR;
    }
  }
  
  
  else if ((strcmp(argv[1],"sectionTest") == 0) || (strcmp(argv[1],"TestSection") == 0) ||
	   (strcmp(argv[1],"testSection") == 0) || (strcmp(argv[1],"SectionForceDeformationTest") == 0)) {
    int count = 1;
    if (argc == 3) {
      if (Tcl_GetInt(interp, argv[2], &count) != TCL_OK) {
	return TCL_ERROR;
      }	  
    }
    theBuilder = new TclSectionTester(theDomain, interp, count);
    if (theBuilder == 0) {
      opserr << "WARNING ran out of memory in creating TclUniaxialMAterialTester model\n";
      return TCL_ERROR;
    }
  }
  
  else {
    opserr <<  "WARNING unknown model builder type\n";
    
    opserr << "WARNING model builder type " << argv[1]
	   << " not supported\n";
    return TCL_ERROR;
  }
  
  
  
  return TCL_OK;
}





