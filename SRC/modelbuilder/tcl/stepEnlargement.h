#ifndef __MAIN_H__
#define __MAIN_H__

// #include <windows.h>

#include <tcl.h>

#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>

#include <float.h>

using std::ifstream;
using namespace std;

// read data from file
vector<double> readPath(string fileNameInput);
// write data to file
void writetofile(string fileNameOutput, vector<double>& datav);
// method base on binary
vector<double> stepEnlargement_MY(int nFactor, vector<double>& datav);
// a technique method
vector<double> stepEnlargement_INT(int nFactor, vector<double>& datav);
/*  To use this exported function of dll, include this header
*  in your project.
*/

#define USE_TCL_STUBS






#endif // __MAIN_H__
#pragma once
