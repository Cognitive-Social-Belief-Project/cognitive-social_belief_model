#include "ModelParameters.h"

using namespace std;

const int BUFFSIZE = 80;

// ModelParameters default constructor
ModelParameters::ModelParameters()
{
	m_sParamFile = "bnetparams.dat";
	GetParams();
}

// ModelParameters constructor
ModelParameters::ModelParameters(string sParamFile)
{
	m_sParamFile = sParamFile;
	GetParams();
}

ModelParameters::~ModelParameters()
{
	// Destroys allocated memory for Assignments
	for (vector<Assignment*>::iterator itIII = m_vpassNodeAssignments.begin(); itIII != m_vpassNodeAssignments.end(); ++itIII)
	{
		delete *itIII;
	}
}

// Read parameter file and write in data
void ModelParameters::GetParams()
{
	// Read file
	ifstream ifParamFile(m_sParamFile.c_str());

	if (!ifParamFile)
	{
		cerr << "Could not read parameter file:" << m_sParamFile << endl;
		exit(1);
	}

	// Buffer to retrieve lines
	string sBuff;
	stringstream ssLineOfFile(stringstream::in | stringstream::out);
	ssLineOfFile.clear();

	// Get output file name
	getline(ifParamFile, sBuff);
	ssLineOfFile << sBuff;
	ssLineOfFile >> m_sOutputFile;
	ssLineOfFile.clear();

	// Get sizes output file name
	getline(ifParamFile, sBuff);
	ssLineOfFile << sBuff;
	ssLineOfFile >> m_sSizesFile;
	ssLineOfFile.clear();

	// Get global energy output file name
	getline(ifParamFile, sBuff);
	ssLineOfFile << sBuff;
	ssLineOfFile >> m_sOutEnergyFile;
	ssLineOfFile.clear();

//	// Get individual energy output file name
//	getline(ifParamFile, sBuff);
//	ssLineOfFile << sBuff;
//	ssLineOfFile >> m_sIndivEnergyFile;
//	ssLineOfFile.clear();

	// Get I
	getline(ifParamFile, sBuff);
	ssLineOfFile << sBuff;
	ssLineOfFile >> m_dI;
	ssLineOfFile.clear();

	// Get T
	getline(ifParamFile, sBuff);
	ssLineOfFile << sBuff;
	ssLineOfFile >> m_dT;
	ssLineOfFile.clear();

	// Get J
	getline(ifParamFile, sBuff);
	ssLineOfFile << sBuff;
	ssLineOfFile >> m_dJ;
	ssLineOfFile.clear();

	// Get number of data points
	getline(ifParamFile, sBuff);
	ssLineOfFile << sBuff;
	ssLineOfFile >> m_ulNumPoints;
	ssLineOfFile.clear();

	// Get number of iterations
	getline(ifParamFile, sBuff);
	ssLineOfFile << sBuff;
	ssLineOfFile >> m_ulNumIterations;
	ssLineOfFile.clear();

	// Get runs
	getline(ifParamFile, sBuff);
	ssLineOfFile << sBuff;
	ssLineOfFile >> m_unRuns;
	ssLineOfFile.clear();

	// Get energy filename
	getline(ifParamFile, sBuff);
	ssLineOfFile << sBuff;
	ssLineOfFile >> m_sEnergyFile;
	ssLineOfFile.clear();

	// Get update method
	getline(ifParamFile, sBuff);
	ssLineOfFile << sBuff;
	ssLineOfFile >> m_sUpdateMethod;
	ssLineOfFile.clear();

	// Get graph filename
	getline(ifParamFile, sBuff);
	ssLineOfFile << sBuff;
	ssLineOfFile >> m_sGraphFile;
	ssLineOfFile.clear();

	// Get number of sets
	getline(ifParamFile, sBuff);
	ssLineOfFile << sBuff;
	ssLineOfFile >> m_unNumSets;
	ssLineOfFile.clear();

	// Get initial belief configuration;
	getline(ifParamFile, sBuff);
	ssLineOfFile << sBuff;
	int ntempHolder;
	ssLineOfFile >> ntempHolder;
	m_bRandomBeliefs = (ntempHolder != 0);
	ssLineOfFile.clear();
	ssLineOfFile.str("");

	// If there is an initial configuration read it, else move on
	if (m_bRandomBeliefs == false)
	{
		getline(ifParamFile, sBuff);
		ssLineOfFile << sBuff;
		ssLineOfFile >> m_nNumBeliefTypes;
		ssLineOfFile.clear();
		ssLineOfFile.str("");

		// Loop through each set
		for (unsigned int zzz=0U; zzz<m_unNumSets; zzz++)
		{
			// Loop through each Belief System
			vector<vector<bool> > vvbBeliefSet;
			for (unsigned int iii=0U; iii<m_nNumBeliefTypes; iii++)
			{

				getline(ifParamFile, sBuff);
				ssLineOfFile << sBuff;
				vector<bool> vbBeliefSystem;
				// Read in the belief system associations
				while (!ssLineOfFile.eof())
				{
					char chBuff[BUFFSIZE];
					ssLineOfFile.getline(chBuff, BUFFSIZE, ',');

					bool bAssociation = (atoi(chBuff) != -1);
					vbBeliefSystem.push_back(bAssociation);
				}

				// Add belief system to stack
				vvbBeliefSet.push_back(vbBeliefSystem);
				ssLineOfFile.str("");
				ssLineOfFile.clear();
			}
			m_vvvbBeliefSets.push_back(vvbBeliefSet);
		}


		// Read in ratios of the initial belief configurations
		getline(ifParamFile, sBuff);
		ssLineOfFile << sBuff;
		for (unsigned int iii=0U; iii<m_nNumBeliefTypes; iii++)
		{

			char chBuff[BUFFSIZE];
			ssLineOfFile.getline(chBuff, BUFFSIZE, ',');

			m_vdBeliefRatios.push_back(atof(chBuff));
		}
	}

	// Get initial zealot configuration
	ssLineOfFile.clear();
	ssLineOfFile.str("");
	getline(ifParamFile, sBuff);
	ssLineOfFile << sBuff;
	ntempHolder = 0;
	ssLineOfFile >> ntempHolder;
	m_bZealots = (ntempHolder != 0);
	ssLineOfFile.clear();

	// If there is an initial zealot configuration read it, else move on
	if (m_bZealots == true)
	{
		getline(ifParamFile, sBuff);
		ssLineOfFile << sBuff;
		ssLineOfFile >> m_nNumZealotTypes;
		ssLineOfFile.clear();
		ssLineOfFile.str("");

		// Loop through all sets
		for (unsigned int zzz=0U; zzz<m_unNumSets; zzz++)
		{
			// Loop through each Zealots belief system
			vector<vector<bool> > vvbZealotSet;
			for (unsigned int iii=0U; iii<m_nNumZealotTypes; iii++)
			{

				getline(ifParamFile, sBuff);
				ssLineOfFile.clear();
				ssLineOfFile.str("");
				ssLineOfFile << sBuff;
				vector<bool> vbBeliefSystem;
				// Read in the belief system associations
				while (!ssLineOfFile.eof())
				{
					char chBuff[BUFFSIZE];
					ssLineOfFile.getline(chBuff, BUFFSIZE, ',');

					bool bAssociation = (atoi(chBuff) != -1);
					vbBeliefSystem.push_back(bAssociation);
				}

				// Add belief system to stack
				vvbZealotSet.push_back(vbBeliefSystem);
			}
			m_vvvbZealotSets.push_back(vvbZealotSet);
		}


		// Read in number of initial zealots for each type of zealot
		getline(ifParamFile, sBuff);
		ssLineOfFile.clear();
		ssLineOfFile.str("");
		ssLineOfFile << sBuff;
		for (unsigned int iii=0U; iii<m_nNumZealotTypes; iii++)
		{
			char chBuff[BUFFSIZE];
			ssLineOfFile.getline(chBuff, BUFFSIZE, ',');

			m_vdZealotRatios.push_back(atof(chBuff));
		}
	}

	// Get size of graph
	ssLineOfFile.clear();
	ssLineOfFile.str("");
	getline(ifParamFile, sBuff);
	ssLineOfFile << sBuff;
	ssLineOfFile >> m_ulGraphSize;
	ssLineOfFile.clear();

	// Get whether there is a seed
	getline(ifParamFile, sBuff);
	ssLineOfFile << sBuff;
	ntempHolder = 0;
	ssLineOfFile >> ntempHolder;
	m_bSeed = (ntempHolder != 0);
	ssLineOfFile.clear();

	// If there is supposed to be a seed, get it
	if (m_bSeed)
	{
		// Get seed
		getline(ifParamFile, sBuff);
		ssLineOfFile << sBuff;
		ssLineOfFile >> m_unSeed;
		ssLineOfFile.clear();
	}

	// Get whether there is node assignments
	ssLineOfFile.clear();
	ssLineOfFile.str("");
	getline(ifParamFile, sBuff);
	ssLineOfFile << sBuff;
	ntempHolder = 0;
	ssLineOfFile >> ntempHolder;
	m_bIsAssignment = (ntempHolder != 0);
	ssLineOfFile.clear();

	// If there are node assignments, get them
	if (m_bIsAssignment == true)
	{
		// Loop through all nodes for assignment
		for (unsigned long iii = 0UL; iii < m_ulGraphSize; iii++)
		{
			getline(ifParamFile, sBuff);
			ssLineOfFile.clear();
			ssLineOfFile.str("");
			ssLineOfFile << sBuff;
			// Read in the belief system associations
			// Create NEW 3 element array to hold buffer
			std::string *pstrAssignments = new std::string [3];
			int nIndex = 0;
			while (!ssLineOfFile.eof())
			{
				char chBuff[BUFFSIZE];
				ssLineOfFile.getline(chBuff, BUFFSIZE, ',');

				pstrAssignments[nIndex] = std::string(chBuff);
				nIndex++;
			}
			// Dump 3 elements into the NEW struct
			ssLineOfFile.clear();
			ssLineOfFile.str("");
			ssLineOfFile << pstrAssignments[0];
			unsigned long ulTempNode = 0UL;
			ssLineOfFile >> ulTempNode;

			ssLineOfFile << pstrAssignments[2];
			unsigned int unTempSubset = 0U;
			ssLineOfFile >> unTempSubset;

			::Assignment *passAssignments = new ::Assignment();
			passAssignments->sType = pstrAssignments[1];
			passAssignments->ulNode = ulTempNode;
			passAssignments->unSubset = unTempSubset;
			// Erase the array
			delete[] pstrAssignments;
			// Add struct to vector
			m_vpassNodeAssignments.push_back(passAssignments);
		}
	}

	// Get which data to write (Time series)
	ssLineOfFile.clear();
	ssLineOfFile.str("");
	getline(ifParamFile, sBuff);
	ssLineOfFile << sBuff;
	ntempHolder = 0;
	ssLineOfFile >> ntempHolder;
	m_bTimeSeries = (ntempHolder != 0);
	ssLineOfFile.clear();

	// Get Final state
	ssLineOfFile.clear();
	ssLineOfFile.str("");
	getline(ifParamFile, sBuff);
	ssLineOfFile << sBuff;
	ntempHolder = 0;
	ssLineOfFile >> ntempHolder;
	m_bFinalState = (ntempHolder != 0);
	ssLineOfFile.clear();

	// TS Sizes
	ssLineOfFile.clear();
	ssLineOfFile.str("");
	getline(ifParamFile, sBuff);
	ssLineOfFile << sBuff;
	ntempHolder = 0;
	ssLineOfFile >> ntempHolder;
	m_bTSSizes = (ntempHolder != 0);
	ssLineOfFile.clear();

	// IsAppendOut()
	ssLineOfFile.clear();
	ssLineOfFile.str("");
	getline(ifParamFile, sBuff);
	ssLineOfFile << sBuff;
	ntempHolder = 0;
	ssLineOfFile >> ntempHolder;
	m_bAppendOut = (ntempHolder != 0);
	ssLineOfFile.clear();
//	// Preproccess Neighbors
//	ssLineOfFile.clear();
//	ssLineOfFile.str("");
//	getline(ifParamFile, sBuff);
//	ssLineOfFile << sBuff;
//	ntempHolder = 0;
//	ssLineOfFile >> ntempHolder;
//	m_bPreProccessNeighors = (ntempHolder != 0);
//	ssLineOfFile.clear();

}

// Set values
void ModelParameters::Set_I(const double dI)
{
	m_dI = dI;
}

void ModelParameters::Set_T(const double dT)
{
	m_dT = dT;
}

void ModelParameters::Set_J(const double dJ)
{
	m_dJ = dJ;
}

void ModelParameters::Set_EnergyFile(const string sFilename)
{
	m_sEnergyFile = sFilename;
}

void ModelParameters::Set_UpdateMethod(const string sMethod)
{
	m_sUpdateMethod = sMethod;
}

void ModelParameters::Set_GraphFile(const string sFilename)
{
	m_sGraphFile = sFilename;
}

void ModelParameters::Set_NonRandomBeliefs(bool bChoice)
{
	m_bRandomBeliefs = bChoice;
}

void ModelParameters::Set_Zealots(bool bChoice)
{
	m_bZealots = bChoice;
}

// Get stuff
const string &ModelParameters::GetOutFile() const
{
	return m_sOutputFile;
}

const string &ModelParameters::GetSizesFile() const
{
	return m_sSizesFile;
}

const string &ModelParameters::GetOutEnergyFile() const
{
	return m_sOutEnergyFile;
}

//bool ModelParameters::IsPreproccessNeighbors() const
//{
//	return m_bPreProccessNeighors;
//}

double &ModelParameters::GetI()
{
	return m_dI;
}

double &ModelParameters::GetT()
{
	return m_dT;
}

double &ModelParameters::GetJ()
{
	return m_dJ;
}

unsigned long ModelParameters::GetNumPoints() const
{
	return m_ulNumPoints;
}

unsigned long long ModelParameters::GetNumIter() const
{
	return m_ulNumIterations;
}

const unsigned int &ModelParameters::GetNumRuns() const
{
	return m_unRuns;
}

const string &ModelParameters::GetUpdateMethod() const
{
	return m_sUpdateMethod;
}

const string &ModelParameters::GetEnergyFile() const
{
	return m_sEnergyFile;
}

const string &ModelParameters::GetGraphFile() const
{
	return m_sGraphFile;
}

const unsigned int &ModelParameters::GetNumSets() const
{
	return m_unNumSets;
}

bool ModelParameters::IsRandom() const
{
	return m_bRandomBeliefs;
}

bool ModelParameters::IsZealots() const
{
	return m_bZealots;
}

vector<double>& ModelParameters::GetBeliefRatios()
{
	return m_vdBeliefRatios;
}

vector<vector<vector<bool> > >& ModelParameters::GetBeliefSets()
{
	return m_vvvbBeliefSets;
}

vector<bool> & ModelParameters::GetBeliefSetVector(unsigned int unSet, unsigned int unBeliefType)
{
	return m_vvvbBeliefSets[unSet][unBeliefType];
}

vector<double>& ModelParameters::GetZealotRatios()
{
	return m_vdZealotRatios;
}

vector<vector<vector<bool> > >& ModelParameters::GetZealotSets()
{
	return m_vvvbZealotSets;
}

vector<bool> & ModelParameters::GetZealotSetVector(unsigned int unSet, unsigned int unBeliefType)
{
	return m_vvvbZealotSets[unSet][unBeliefType];
}

bool ModelParameters::IsSeed() const
{
	return m_bSeed;
}

unsigned int ModelParameters::GetSeed() const
{
	return m_unSeed;
}

vector<Assignment*> & ModelParameters::GetAssignmentList()
{
	return m_vpassNodeAssignments;
}

Assignment* ModelParameters::GetNodeAssignment(unsigned long ulIndex)
{
	return m_vpassNodeAssignments[ulIndex];
}

bool ModelParameters::IsAssigned() const
{
	return m_bIsAssignment;
}

bool ModelParameters::IsFinalStates() const
{
	return m_bFinalState;
}

bool ModelParameters::IsTimeSeries() const
{
	return m_bTimeSeries;
}

bool ModelParameters::IsTSSizes() const
{
	return m_bTSSizes;
}

bool ModelParameters::IsAppendOut() const
{
	return m_bAppendOut;
}

//------------------------------------------------------------------------------
//
//	Friend functions
//
//------------------------------------------------------------------------------
ostream& operator<< (ostream &rout, ModelParameters &rcModelParameters)
{
    // Since operator<< is a friend of the Point class, we can access
    // Point's members directly.

	string sOutString = "";
	stringstream stream; // for converting types to strings
	string sConverted = "";

	sOutString += "Parameter file:\t\t" + rcModelParameters.m_sParamFile + "\n";

	// Write out Output filename
	stream << rcModelParameters.m_sOutputFile;
	stream >> sConverted;
	stream.clear();
	sOutString += "Output File:\t\t" + sConverted + "\n";

	// Write out sizes filename
	stream << rcModelParameters.m_sSizesFile;
	stream >> sConverted;
	stream.clear();
	sOutString += "Sizes File:\t\t" + sConverted + "\n";

	// Write out output energy filename
	stream << rcModelParameters.m_sOutEnergyFile;
	stream >> sConverted;
	stream.clear();
	sOutString += "Out-Energy File:\t\t" + sConverted + "\n";

	// Write out I
	stream << rcModelParameters.m_dI;
	stream >> sConverted;
	stream.clear();
	sOutString += "Peer-Inf (I):\t\t" + sConverted + "\n";

	// Write out Num Points
	stream << rcModelParameters.m_ulNumPoints;
	stream >> sConverted;
	stream.clear();
	sOutString += "Number of points:\t\t" + sConverted + "\n";

	// Write out number of iterations
	stream << rcModelParameters.m_ulNumIterations;
	stream >> sConverted;
	stream.clear();
	sOutString += "Number of Iterations:\t\t" + sConverted + "\n";

	// Write out number of runs
	stream << rcModelParameters.m_unRuns;
	stream >> sConverted;
	stream.clear();
	sOutString += "Number of runs:\t\t" + sConverted + "\n";

	// Write out energy file name
	sOutString += "Energy file:\t\t" + rcModelParameters.m_sEnergyFile + "\n";
	sOutString += "Update Method:\t\t" + rcModelParameters.m_sUpdateMethod + "\n";
	sOutString += "Graph file:\t\t" + rcModelParameters.m_sGraphFile + "\n";

	// Write out number of sets
	stream << rcModelParameters.m_unNumSets;
	stream >> sConverted;
	stream.clear();
	sOutString += "Number of Sets:\t\t" + sConverted + "\n";

	// Write out choice for initial belief configuration
	stream << rcModelParameters.m_bRandomBeliefs;
	stream >> sConverted;
	stream.clear();
	sOutString += "Random Beliefs:\t\t" + sConverted + "\n";

	// Write out choice for zealots
	stream << rcModelParameters.m_bZealots;
	stream >> sConverted;
	stream.clear();
	sOutString += "Zealots:\t\t" + sConverted + "\n";

	// Write out intial belief configurations
	if (!rcModelParameters.m_bRandomBeliefs)
	{
		// determine number of belief types
		stream.clear();
		stream.str("");
		stream << rcModelParameters.m_nNumBeliefTypes;
		stream >> sConverted;
		stream.clear();
		sOutString += "Number of Belief-sys:\t\t" + sConverted + "\n";

		for (unsigned int iii=0U; iii<rcModelParameters.m_nNumBeliefTypes; iii++)
		{
			stream << rcModelParameters.m_vdBeliefRatios[iii];
			stream >> sConverted;
			stream.clear();
			sOutString += "\tBelief Ratios:\t\t" + sConverted + "\n";
		}


		for (vector<vector<vector<bool> > >::iterator itIII = rcModelParameters.m_vvvbBeliefSets.begin(); itIII != rcModelParameters.m_vvvbBeliefSets.end(); ++itIII)
		{

			sOutString += "\n";
			for (vector<vector<bool> >::iterator itJJJ = itIII->begin(); itJJJ != itIII->end(); ++itJJJ)
			{
				stream.clear();
				stream.str("");
				sOutString += "\t[ ";
				for (vector<bool>::iterator itKKK = itJJJ->begin(); itKKK != itJJJ->end(); ++itKKK)
				{
					stream << *itKKK;
					stream >> sConverted;
					stream.clear();
					sOutString += sConverted + ' ';
				}
				sOutString += "]\n";
			}
			sOutString += "\n";
		}
	}

	// Write out zealot configurations
	if (rcModelParameters.m_bZealots)
	{
		stream.clear();
		stream.str("");
		stream << rcModelParameters.m_nNumZealotTypes;
		stream >> sConverted;
		stream.clear();
		sOutString += "Number of Zealot-Types:\t\t" + sConverted + "\n";
		for (unsigned int iii=0U; iii<rcModelParameters.m_nNumZealotTypes; iii++)
		{
			stream << rcModelParameters.m_vdZealotRatios[iii];
			stream >> sConverted;
			stream.clear();
			sOutString += "\tNumber of Zealots:\t\t" + sConverted + "\n";
		}

		for (vector<vector<vector<bool> > >::iterator itIII = rcModelParameters.m_vvvbZealotSets.begin(); itIII != rcModelParameters.m_vvvbZealotSets.end(); ++itIII)
		{

			sOutString += "\n";
			for (vector<vector<bool> >::iterator itJJJ = itIII->begin(); itJJJ != itIII->end(); ++itJJJ)
			{
				stream.clear();
				stream.str("");
				sOutString += "\t[ ";
				for (vector<bool>::iterator itKKK = itJJJ->begin(); itKKK != itJJJ->end(); ++itKKK)
				{
					stream << *itKKK;
					stream >> sConverted;
					stream.clear();
					sOutString += sConverted + ' ';
				}
				sOutString += "]\n";
			}
			sOutString += "\n";
		}
	}

    rout << sOutString;
    return rout;
}
