#include "stepEnlargement.h"

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