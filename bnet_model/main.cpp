#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include "Resources.h"
#include "ModelParameters.h"
#include "EdgeGraph.h"
#include "CachedEnergy.h"
#include "Agent.h"
#include "BNetModel.h"

int main(int argc, char *argv[])
{
	/*
	 * In order for the data to match the graph node labels, labelling must be
	 * done contiguously with NO missing nodes. So if the graph has 100 nodes
	 * nodes must take values from 0-99. If any numbers are missing the
	 * program will fail.
	 */
	using namespace std;

	string sParamFile("");
	if (argc == 2)
	{
		sParamFile = argv[1];
	}
	else
	{
		sParamFile = "test_params.dat";
	}

	ModelParameters modelparParameters(sParamFile);
	BNetModel bnetModel(modelparParameters);
	bnetModel.Engage();

//	cout << "Accept: " << accept << endl;
//	cout << "Flip: " << flip << endl;
//	cout << "Stay: " << stay << endl;
//	cout << "Total: " << accept + flip + stay << endl;

	return 0;
}
