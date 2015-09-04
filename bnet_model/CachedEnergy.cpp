#include "CachedEnergy.h"

using namespace std;

CachedEnergy::CachedEnergy()
{
}

CachedEnergy::CachedEnergy(string sEnergyFile, unsigned int unSeed /*=static_cast<unsigned int>(std::time(0))*/)
{
	// Read file
	ifstream ifEnergyFile(sEnergyFile.c_str());

	if (!ifEnergyFile)
	{
		cerr << "Could not read parameter file:" << sEnergyFile << endl;
		exit(1);
	}

	// Buffer to retrieve lines
	string sBuff;
	stringstream ssLineOfFile;

	// Read number of associations (1st line)
	getline(ifEnergyFile, sBuff);
	ssLineOfFile << sBuff;
	ssLineOfFile >> m_unNumAssociations;
	ssLineOfFile.clear();
	ssLineOfFile.str("");

	// Read in energies and beliefs
	vector<string> vsLineOfFile;
	while (getline(ifEnergyFile, sBuff))
	{
		vsLineOfFile.push_back(sBuff);
	}

	// Split lines by ',' and add to member variables
	for (vector<string>::iterator itIII = vsLineOfFile.begin(); itIII != vsLineOfFile.end(); ++itIII)
	{

		vector<bool> vbListOfStates;
		vector<string> vsSplitLines = StringSplit(*itIII, ',');
		for (unsigned int jjj = 0U; jjj < m_unNumAssociations; jjj++)
		{
			vbListOfStates.push_back((vsSplitLines[jjj] != "-1"));
		}

		m_vdListOfEnergies.push_back(atof(vsSplitLines.back().c_str()));
		m_vvbListOfBeliefSys.push_back(vbListOfStates);
	}

	m_ulNumEnergies = m_vdListOfEnergies.size();

	// Create random belief selector from data
	unsigned long ulMin = 0UL;
	m_crngRandomBeliefSelector = new BoostRandomInt<unsigned long>(ulMin, m_ulNumEnergies-1, unSeed);

	// Create random state selector from data
	unsigned int unMin= 0U;
	m_crngRandomStateSelector = new BoostRandomInt<unsigned int>(unMin, m_unNumAssociations-1, unSeed);

	// Generate hash table look-up (unordered_
	m_pumapBeliefSysLookUp = new boost::unordered_map<vector<bool>, unsigned long >;

	// Fill table (vector<bool> -> unsigned long)
	for (unsigned long iii = 0UL; iii < m_vvbListOfBeliefSys.size(); iii++)
	{
		(*m_pumapBeliefSysLookUp).insert(make_pair(m_vvbListOfBeliefSys[iii], iii));
	}

	// Fill table (unsigned long -> vector<int>)
	m_pumapBeliefVect = new boost::unordered_map<unsigned long, std::vector<int>* >;
	for (unsigned long iii = 0UL; iii < m_vvbListOfBeliefSys.size(); iii++)
	{
		(*m_pumapBeliefVect)[iii] = GetIntBelief(m_vvbListOfBeliefSys[iii]);
	}

}

CachedEnergy::~CachedEnergy()
{
	delete m_crngRandomStateSelector;
	delete m_crngRandomBeliefSelector;
	delete m_pumapBeliefSysLookUp;
//	for_each( m_pumapBeliefVect->begin(),m_pumapBeliefVect.end(),
//          DeleteVector<vector<int>*>());
	for (unsigned long iii = 0L; iii < m_pumapBeliefVect->size(); iii++)
	{
		delete (*m_pumapBeliefVect)[iii];
	}
	m_pumapBeliefVect->clear();
	delete m_pumapBeliefVect;
}

unsigned long CachedEnergy::RandomBelief() const
{
	return (*m_crngRandomBeliefSelector)();
}

unsigned int CachedEnergy::RandomState() const
{
	return (*m_crngRandomStateSelector)();
}

vector<bool>& CachedEnergy::GetBeliefSys(unsigned long ulBeliefIndex)
{
	return m_vvbListOfBeliefSys[ulBeliefIndex];
}

vector<int>* CachedEnergy::GetBeliefVector(unsigned long ulBeliefIndex)
{
	return (*m_pumapBeliefVect)[ulBeliefIndex];
}


bool CachedEnergy::GetStateValue(unsigned long ulBeliefIndex, unsigned int unStateIndex) const
{
	return m_vvbListOfBeliefSys[ulBeliefIndex][unStateIndex];
}

double CachedEnergy::GetEnergy(unsigned long ulBeliefIndex) const
{
	return m_vdListOfEnergies[ulBeliefIndex];
}

unsigned int CachedEnergy::GetNumAssoc() const
{
	return m_unNumAssociations;
}

unsigned long CachedEnergy::Search(vector<bool> &key) const
{
	return (*m_pumapBeliefSysLookUp).at(key);
}

unsigned long CachedEnergy::GetNeighborState(unsigned long ulBeliefIndex, unsigned int unStateToFlip)
{
	// Temporarily flip state
	m_vvbListOfBeliefSys[ulBeliefIndex][unStateToFlip].flip();

	// Search for this flipped state
	unsigned long ulNeigborState = Search(m_vvbListOfBeliefSys[ulBeliefIndex]);

	// Flip state back
	m_vvbListOfBeliefSys[ulBeliefIndex][unStateToFlip].flip();

	// Return neighbor locations
	return ulNeigborState;
}

void CachedEnergy::KeyValueCheck()
{
	string sOutString = "";
	stringstream stream; // for converting types to strings
	string sConverted = "";

//	for (vector<vector<bool> >::iterator itIII = m_vvbListOfBeliefSys.begin(); itIII != m_vvbListOfBeliefSys.end(); ++itIII)
//	{
//		sOutString += "[ ";
//		stream.str("");
//		stream.clear();
//		sConverted = "";
//		for (vector<bool>::iterator itJJJ = itIII->begin(); itJJJ != itIII->end(); ++itJJJ)
//		{
//			stream << *itJJJ;
//			stream >> sConverted;
//			stream.clear();
//			sOutString += sConverted + " ";
//		}
//		sOutString += "] ";
//		stream << Search(*itIII);
//		stream >> sConverted;
//		sOutString += sConverted + "\n";
//	}

	sOutString += "\n\n\n";
	for (unsigned long iii = 0UL; iii < m_ulNumEnergies; iii++)
	{
		sOutString += "[ ";
		stream.str("");
		stream.clear();
		sConverted = "";
		for (unsigned int jjj = 0UL; jjj < m_unNumAssociations; jjj++)
		{
			stream << m_vvbListOfBeliefSys[iii][jjj];
			stream >> sConverted;
			stream.clear();
			sOutString += sConverted + " ";
		}
		sOutString += "] ";

		stream << Search(m_vvbListOfBeliefSys[iii]);
		stream >> sConverted;
		stream.clear();
		sOutString += sConverted;

		stream << iii;
		stream >> sConverted;
		stream.clear();
		sOutString += " " + sConverted;

		stream << m_vdListOfEnergies[iii];
		stream >> sConverted;
		stream.clear();
		sOutString += " " + sConverted + "\n";
	}

	cout << sOutString << endl;
}

vector<int>* CachedEnergy::GetIntBelief(vector<bool> &rvbBelief)
{
	vector<int> *pvnBelief = new vector<int>;

	for (vector<bool>::iterator itIII = rvbBelief.begin(); itIII != rvbBelief.end(); ++itIII)
	{
		if (*itIII)
		{
			pvnBelief->push_back(1);
		}
		else
		{
			pvnBelief->push_back(-1);
		}
	}

	return pvnBelief;
}

//------------------------------------------------------------------------------
//
//	Friend functions
//
//------------------------------------------------------------------------------
ostream& operator<< (ostream &rout, CachedEnergy &rcCachedEnergy)
{
	string sOutString = "";
	stringstream stream; // for converting types to strings
	string sConverted = "";

	stream << rcCachedEnergy.m_unNumAssociations;
	stream >> sConverted;
	stream.clear();
	sOutString += "Number of Assoc: " + sConverted + "\n";

	stream << rcCachedEnergy.m_ulNumEnergies;
	stream >> sConverted;
	stream.clear();
	sOutString += "Number of Belief-Sys: " + sConverted + "\n";

	// Print out belief systems and corresponding energies
	for (unsigned long iii = 0UL; iii < rcCachedEnergy.m_ulNumEnergies; iii++)
	{
		stream.clear();
		stream.str("");
		sOutString += "\t[ ";

		for (unsigned int jjj = 0; jjj < rcCachedEnergy.m_unNumAssociations; jjj++)
		{
			stream << rcCachedEnergy.m_vvbListOfBeliefSys[iii][jjj];
			stream >> sConverted;
			stream.clear();
			sOutString += sConverted + " ";
		}

		stream << rcCachedEnergy.m_vdListOfEnergies[iii];
		stream >> sConverted;
		stream.clear();
		stringstream counter;
		counter << iii;
		sOutString += " ] " + sConverted + "\t\t " + counter.str() + "\n";
	}

	rout << sOutString;
	return rout;
}
