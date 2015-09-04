#include "BNetModel.h"
#include <cassert>
#include <numeric>
#include <algorithm>

using namespace std;

BNetModel::BNetModel(ModelParameters &rcParameters)
	: m_rcBNetParameters(rcParameters)
{
	// Generate Energy states and graph
	if (m_rcBNetParameters.IsSeed())
	{
		m_pcEnergyStates = new CachedEnergy(m_rcBNetParameters.GetEnergyFile(), m_rcBNetParameters.GetSeed());
		m_pcGraph = new EdgeGraph(m_rcBNetParameters.GetGraphFile(), m_rcBNetParameters.GetSeed());

		// Create random number generators
		m_prngReal0to1 = new BoostRandomReal<double>(0.0,1.0, m_rcBNetParameters.GetSeed());
		m_prngUInt0or1 = new BoostRandomInt<unsigned int>(0,1, m_rcBNetParameters.GetSeed());
	}
	else
	{
		m_pcEnergyStates = new CachedEnergy(m_rcBNetParameters.GetEnergyFile());
		m_pcGraph = new EdgeGraph(m_rcBNetParameters.GetGraphFile());

		// Create random number generators
		m_prngReal0to1 = new BoostRandomReal<double>(0.0,1.0);
		m_prngUInt0or1 = new BoostRandomInt<unsigned int>(0,1);
	}

	// Set to null
	m_eEmit.state = 0U;
	m_eEmit.value = false;
	m_ullFirstBeliefIndex = 0ULL;

	// Populate with initial set of random agents
	InitialPopulate();

	// Create Agent neighbors hash
	CreateNeighborHash();
	AssignAgentNeighbors();

	// Determine initial energy
	m_dTotalEnergy = CalcTotalEnergy();
}

BNetModel::~BNetModel()
{
	for (vector<Agent*>::iterator itIII = m_vpagentAgentList.begin(); itIII != m_vpagentAgentList.end(); ++itIII)
	{
		delete *itIII;
	}
	delete m_pcEnergyStates;
	delete m_pcGraph;
	delete m_prngReal0to1;
	delete m_prngUInt0or1;
	for (unsigned long long iii = 0ULL; iii < m_pumapAgentNeighbors->size(); iii++)
	{
		delete (*m_pumapAgentNeighbors)[iii];
	}
	m_pumapAgentNeighbors->clear();
	delete m_pumapAgentNeighbors;
}

void BNetModel::CreateNeighborHash()
{
	// Create new hash
	m_pumapAgentNeighbors = new boost::unordered_map<unsigned long long, std::vector<Agent*>* >;

	// For each edge in network
	for (unsigned long long iii = 0ULL; iii < m_pcGraph->GetNumEdges(); iii++)
	{
		// ADD FOR FIRST NODE
		// if node in hash push neighbor into vector
		if (m_pumapAgentNeighbors->find(m_pcGraph->PickNode(iii, 0U)) != m_pumapAgentNeighbors->end())
		{
			Agent* pagentNeighbor = m_vpagentAgentList[m_pcGraph->PickNode(iii, 1U)];
			(*m_pumapAgentNeighbors)[m_pcGraph->PickNode(iii, 0U)]->push_back(pagentNeighbor);
		}
		// if node not in hash create vector and push neighbor into vector
		else
		{
			vector<Agent*> *pvagentNeighbors = new vector<Agent*>;
			Agent* pagentNeighbor = m_vpagentAgentList[m_pcGraph->PickNode(iii, 1U)];
			pvagentNeighbors->push_back(pagentNeighbor);
			m_pumapAgentNeighbors->insert(make_pair(m_pcGraph->PickNode(iii, 0U), pvagentNeighbors));
		}

		// ADD FOR SECOND NODE
		// if node in hash push neighbor into vector
		if (m_pumapAgentNeighbors->find(m_pcGraph->PickNode(iii, 1U)) != m_pumapAgentNeighbors->end())
		{
			Agent* pagentNeighbor = m_vpagentAgentList[m_pcGraph->PickNode(iii, 0U)];
			(*m_pumapAgentNeighbors)[m_pcGraph->PickNode(iii, 1U)]->push_back(pagentNeighbor);
		}
		// if node not in hash create vector and push neighbor into vector
		else
		{
			vector<Agent*> *pvagentNeighbors = new vector<Agent*>;
			Agent* pagentNeighbor = m_vpagentAgentList[m_pcGraph->PickNode(iii, 0U)];
			pvagentNeighbors->push_back(pagentNeighbor);
			m_pumapAgentNeighbors->insert(make_pair(m_pcGraph->PickNode(iii, 1U), pvagentNeighbors));
		}
	}
}

void BNetModel::AssignAgentNeighbors()
{
	// Loop through each agent and call agents to assign their appropriate neighbors
	for (unsigned long long iii = 0ULL; iii < m_vpagentAgentList.size(); iii++)
	{
		m_vpagentAgentList[iii]->SetNeighbors((*m_pumapAgentNeighbors)[iii]);
	}
}

void BNetModel::InitialPopulate()
{

	// Create random initial population (before configuring for specific types)
	for (unsigned long long iii = 0ULL; iii < m_pcGraph->GetNumNodes(); iii++)
	{
		Agent *pagentAddAgent = new Agent(m_pcEnergyStates, &(m_rcBNetParameters.GetI()), &(m_rcBNetParameters.GetT()), &(m_rcBNetParameters.GetJ()), m_prngReal0to1);
		m_vpagentAgentList.push_back(pagentAddAgent);
	}

	// If zealots are present, find how many there are and set the zealot index.
	if (m_rcBNetParameters.IsAssigned())
	{
		// Set the FirstBeliefIndex here:

		// Loop through each assignment
		m_ullFirstBeliefIndex = 0ULL;
		vector<Assignment*> &rvpassNodeAssignments = m_rcBNetParameters.GetAssignmentList();
		for (vector<Assignment*>::iterator itIII = rvpassNodeAssignments.begin(); itIII != rvpassNodeAssignments.end(); ++itIII)
		{
			if ((*itIII)->sType == "z")
			{
				m_ullFirstBeliefIndex++;
			}
		}
	}
	else if (m_rcBNetParameters.IsZealots())
	{
		// Important, to be consistent we must convert ratio->numZeal for each ratio BEFORE adding those numbers
		// this must be done so that no float addition errors carry over to the total number of zealots
		vector<double> &rvdZealotRatios = m_rcBNetParameters.GetZealotRatios();

		// Convert ratios to indexes for agent list
		vector<unsigned long long> vullRatioConversion;
		for (unsigned int iii = 0U; iii < rvdZealotRatios.size(); iii++)
		{
			vullRatioConversion.push_back(rvdZealotRatios[iii] * m_pcGraph->GetNumNodes());
		}

		// Find total number of zealots
		for (unsigned int iii = 0U; iii < vullRatioConversion.size(); iii++)
		{
			m_ullFirstBeliefIndex += vullRatioConversion[iii];
		}
	}
}

void BNetModel::RePopulate(unsigned int unSetNumber)
{

	// If random agents
	if (m_rcBNetParameters.IsRandom())
	{
		// have agents get new random value
		for (vector<Agent*>::iterator itIII = m_vpagentAgentList.begin(); itIII != m_vpagentAgentList.end(); ++itIII)
		{
			// Calls default version of set belief which assigns random belief
			(*itIII)->SetBelief();
		}
	}
	else if (m_rcBNetParameters.IsAssigned())
	{
		// Loop through each assignment
		vector<Assignment*> &rvpassNodeAssignments = m_rcBNetParameters.GetAssignmentList();
		for (vector<Assignment*>::iterator itIII = rvpassNodeAssignments.begin(); itIII != rvpassNodeAssignments.end(); ++itIII)
		{
			if ((*itIII)->sType == "z")
			{
				m_vpagentAgentList[(*itIII)->ulNode]->SetBelief(m_pcEnergyStates->Search(m_rcBNetParameters.GetZealotSetVector(unSetNumber, (*itIII)->unSubset)));
			}
			else
			{
				m_vpagentAgentList[(*itIII)->ulNode]->SetBelief(m_pcEnergyStates->Search(m_rcBNetParameters.GetBeliefSetVector(unSetNumber, (*itIII)->unSubset)));
			}
		}
	}
	else
	{
		// get ratios for current set
		assert(unSetNumber <= m_rcBNetParameters.GetNumSets());

		// loop through ratios and multiply them by size o network
		vector<double> &rvdRatios = m_rcBNetParameters.GetBeliefRatios();

		// Convert ratios to indexes for agent list as a cumulative sum
		vector<unsigned long long> vullRatioConversion;
		vullRatioConversion.push_back(m_ullFirstBeliefIndex); // Start at index following last zealot
		for (unsigned int iii = 0U; iii < rvdRatios.size(); iii++)
		{ // this will break if the user gives ratios that do not add to 1.0
			vullRatioConversion.push_back(static_cast<unsigned long long>(rvdRatios[iii] * m_pcGraph->GetNumNodes()) - vullRatioConversion[iii]);
		}
		// Set final index to max size of network
		vullRatioConversion[rvdRatios.size()] = m_pcGraph->GetNumNodes();

		// change agent beliefs by looping through indexes
		for (unsigned int iii = 0U; iii < vullRatioConversion.size() - 1; iii++)
		{
			for (unsigned long long jjj = vullRatioConversion[iii]; jjj < vullRatioConversion[iii+1]; jjj++)
			{
				// Get vector -> search for vector # -> set agent to #
				m_vpagentAgentList[jjj]->SetBelief(m_pcEnergyStates->Search(m_rcBNetParameters.GetBeliefSetVector(unSetNumber, iii)));
			}
		}
	}

	// If zealots included
	if (m_rcBNetParameters.IsZealots() && !m_rcBNetParameters.IsAssigned())
	{

		// Re-assign the beliefs of the zealots for this particular set
		vector<double> &rvdZealotRatios = m_rcBNetParameters.GetZealotRatios();

		// Convert ratios to indexes
		vector<unsigned long long> vullNumZealots;
		vullNumZealots.push_back(0U);
		for (unsigned int iii = 0U; iii < rvdZealotRatios.size(); iii++)
		{
			vullNumZealots.push_back(static_cast<unsigned long long>(rvdZealotRatios[iii] * m_pcGraph->GetNumNodes()) + vullNumZealots[iii]);
		}
		// Final zealot index to max
		vullNumZealots[rvdZealotRatios.size()] = m_ullFirstBeliefIndex;

		// Loop through each zealot type
		for (unsigned int iii = 0U; iii < vullNumZealots.size() - 1; iii++)
		{
			// This will loop through each set of zealots, up to, but not including the first belief index
			for (unsigned long long jjj = vullNumZealots[iii]; jjj < vullNumZealots[iii+1]; jjj++)
			{
				m_vpagentAgentList[jjj]->SetBelief(m_pcEnergyStates->Search(m_rcBNetParameters.GetZealotSetVector(unSetNumber, iii)));
			}
		}

	}
}

void BNetModel::Iterate()
{
	// Pick random edge
	unsigned long long ullChosenEdge = m_pcGraph->RandomEdge();

	// Pick Receiver
	unsigned int nReceiver = (*m_prngUInt0or1)();

	// Get chosen node
	unsigned long long ullChosenReceiver = m_pcGraph->PickNode(ullChosenEdge, nReceiver);

	// Check if reciever is a zealot (if not, then commence alteration, else continue)
	if (ullChosenReceiver >= m_ullFirstBeliefIndex)
	{

		// Get sender
		unsigned long long ullChosenSender(0ULL);
		if (nReceiver == 0U)
		{
			ullChosenSender = m_pcGraph->PickNode(ullChosenEdge, 1U );
		}
		else
		{
			ullChosenSender = m_pcGraph->PickNode(ullChosenEdge, 0U );
		}
		// m_eEmit is modified in place
		m_vpagentAgentList[ullChosenSender]->Emit(m_eEmit);
		m_vpagentAgentList[ullChosenReceiver]->Input(m_eEmit, m_dTotalEnergy);
	}

//
//	// If check if zealots are active
//	if (m_rcBNetParameters.IsZealots())
//	{
//
//	}
//	else
//	{
//		// Get chosen node
//		unsigned long long ullChosenReceiver = m_pcGraph->PickNode(ullChosenEdge, nReceiver);
//
//		// Get sender
//		unsigned long long ullChosenSender(0ULL);
//		if (nReceiver == 0U)
//		{
//			ullChosenSender = m_pcGraph->PickNode(ullChosenEdge, 1U );
//		}
//		else
//		{
//			ullChosenSender = m_pcGraph->PickNode(ullChosenEdge, 0U );
//		}
//
//		// m_eEmit is modified in place
//		m_vpagentAgentList[ullChosenSender]->Emit(m_eEmit);
//		m_vpagentAgentList[ullChosenReceiver]->Input(m_eEmit);
//	}
}

void BNetModel::RunSet(unsigned int unSetID)
{
	// Run model for chosen number of iterations
	unsigned long long ullMaxIter = m_rcBNetParameters.GetNumIter();
	unsigned long ulNumPoints = m_rcBNetParameters.GetNumPoints();
	unsigned long long ullUpdatePoint = ullMaxIter / static_cast<unsigned long long>(ulNumPoints);
	for (unsigned long long iii = 0ULL; iii < ullMaxIter; iii++)
	{
		if (((iii % ullUpdatePoint) == 0) && m_rcBNetParameters.IsTimeSeries())
		{
			// Append states to file
			if (m_rcBNetParameters.IsAppendOut())
			{
				AppendCurrentStates(m_rcBNetParameters.GetOutFile());
			}
			if (m_rcBNetParameters.IsTSSizes())
			{
				AppendSizesToFile(m_rcBNetParameters.GetSizesFile(), unSetID);
			}
//			if (!m_rcBNetParameters.IsAssigned())
//			{
			AppendEnergyToFile(m_rcBNetParameters.GetOutEnergyFile());
//			}

			Iterate();
		}
		else
		{
			Iterate();
		}
		if ((iii == ullMaxIter - 1) && m_rcBNetParameters.IsFinalStates())
		{
			AppendCurrentStates(m_rcBNetParameters.GetOutFile());
			AppendEnergyToFile(m_rcBNetParameters.GetOutEnergyFile());
			AppendSizesToFile(m_rcBNetParameters.GetSizesFile(), unSetID);
		}
	}
}

void BNetModel::Engage()
{
	for (unsigned int iii = 0U; iii < m_rcBNetParameters.GetNumSets(); iii++)
	{
		// optional callout
		cout << "Set: " << iii << "/" << m_rcBNetParameters.GetNumSets() << endl;
		//
		RePopulate(iii);
		RunSet(iii);
	}
}

void BNetModel::AppendCurrentStates(const string &rsFilename)
{
	// Open single file for appending
	ofstream ofDataFile;
	ofDataFile.open(rsFilename.c_str(), ios::out | ios::app | ios::binary);

	if (ofDataFile.is_open())
	{
		// Get string of belief and dump to file
		string sOutLine;
		for (std::vector<Agent*>::iterator itIII = m_vpagentAgentList.begin(); itIII != m_vpagentAgentList.end(); ++itIII)
		{
			// Convert to string (w/ \n)
			sOutLine += ConvertBeliefToString((*itIII)->GetBeliefIndex());
			// Add the energy
			stringstream ss;
			ss << (*itIII)->GetEnergy();
			sOutLine += ',' + ss.str() + "\n";
		}

		ofDataFile << sOutLine;
		ofDataFile.close();
	}
	else
	{
		cout << "error writting to file: " + rsFilename << endl;
		exit(1);
	}
}

string BNetModel::ConvertBeliefToString(unsigned long ulBelief)
{
	string sBelief("");
	vector<bool> &rvbBelief = m_pcEnergyStates->GetBeliefSys(ulBelief);
	for (unsigned int iii = 0U; iii < rvbBelief.size(); iii++)
	{
		if (iii != 0U)
		{
			sBelief += ',';
		}
		if (rvbBelief[iii])
		{
			sBelief += "1";
		}
		else
		{
			sBelief += "-1";
		}
	}

	return sBelief;
}

void BNetModel::UniqueAgents(boost::unordered_map<unsigned long, unsigned long long >* pbeliefHash)
{

	for (vector<Agent*>::iterator itIII = m_vpagentAgentList.begin(); itIII != m_vpagentAgentList.end(); ++itIII)
	{
		if (pbeliefHash->find((*itIII)->GetBeliefIndex()) != pbeliefHash->end())
		{
			(*pbeliefHash)[(*itIII)->GetBeliefIndex()] += 1ULL;
		}
		else
		{
			pbeliefHash->insert(make_pair((*itIII)->GetBeliefIndex(),1ULL));
		}
	}
}

void BNetModel::AppendSizesToFile(const string& rsSizesFile, unsigned int unSetID)
{
	// Create hash of agents with values that correspond to the number of agents with that belief
	boost::unordered_map<unsigned long, unsigned long long >* pbeliefHash = new boost::unordered_map<unsigned long, unsigned long long >;
	UniqueAgents(pbeliefHash);

	// Open single file for appending
	ofstream ofDataFile;
	ofDataFile.open(rsSizesFile.c_str(), ios::out | ios::app | ios::binary);

	if (ofDataFile.is_open())
	{
		// Get string of belief and dump to file
		string sOutLine("");
		stringstream ssConvert;
		ssConvert << unSetID;
		sOutLine += ssConvert.str() + ",NEWPOINT\n"; // Necessary to know when a new data point begins
		for (boost::unordered_map<unsigned long, unsigned long long >::iterator itIII = pbeliefHash->begin(); itIII != pbeliefHash->end(); ++itIII)
		{
			// Convert to string (w/ \n)
			sOutLine += ConvertBeliefToString(itIII->first);
			stringstream ss;
			ss << itIII->second;
			sOutLine += ',' + ss.str() + "\n";
		}

		ofDataFile << sOutLine;
		ofDataFile.close();
	}
	else
	{
		cout << "error writing to file: " + rsSizesFile << endl;
		exit(1);
	}

	delete pbeliefHash;
}

//void BNetModel::AppendTopNSizesToFile(const string& rsSizesFile, unsigned int unSetID, unsigned long ulTopN = 10UL)
//{
//	// Create hash of agents with values that correspond to the number of agents with that belief
//	boost::unordered_map<unsigned long, unsigned long long >* pbeliefHash = new boost::unordered_map<unsigned long, unsigned long long >;
//	UniqueAgents(pbeliefHash);
//
//	// Get sorted values
//	vector<unsigned long long> vullSizes;
//	for (boost::unordered_map<unsigned long, unsigned long long >::iterator itIII = pbeliefHash->begin(); itIII != pbeliefHash->end(); ++itIII)
//	{
//		vullSizes.push_back(itIII->second);
//	}
//	sort(vullSizes.begin(), vullSizes.end());
//
//	// Write largest groups to file
//	// Open single file for appending
//	ofstream ofDataFile;
//	ofDataFile.open(rsSizesFile.c_str(), ios::out | ios::app | ios::binary);
//
//	if (ofDataFile.is_open())
//	{
//		// Get string of belief and dump to file
//		string sOutLine("");
//		stringstream ssConvert;
//		ssConvert << unSetID;
//		sOutLine += ssConvert.str() + ",NEWPOINT\n"; // Necessary to know when a new data point begins
//
//		for (vector<unsigned long long>::iterator itIII = vullSizes.begin(); itIII != vullSizes.end(); ++itIII)
//		{
//			// Convert to string (w/ \n)
//			sOutLine += ConvertBeliefToString(itIII->first);
//			stringstream ss;
//			ss << itIII->second;
//			sOutLine += ',' + ss.str() + "\n";
//		}
//	}
//
//
//
//}

void BNetModel::AppendEnergyToFile(const std::string& rsEnergyFile)
{
	// Open single file for appending
	ofstream ofDataFile;
	ofDataFile.open(rsEnergyFile.c_str(), ios::out | ios::app | ios::binary);

	if (ofDataFile.is_open())
	{
		// Get string of belief and dump to file
		stringstream ss;

		ss << CalcGlobalEnergy();
		ss << ',';
		ss << CalcAvgIndivEnergy();
		ss << ',';
		ss << m_dTotalEnergy;
		ss << "\n";

		ofDataFile << ss.str();
		ofDataFile.close();
	}
	else
	{
		cout << "error writting to file: " + rsEnergyFile << endl;
		exit(1);
	}
}

double BNetModel::CalcGlobalEnergy()
{
	// Only works for edge graphs
	// sum variable long long
	long long llGlobalEnergy = 0LL;
	long long llMinGlobalEnergy = -1LL * static_cast<long long>(m_pcGraph->GetNumEdges()) * static_cast<long long>(m_pcEnergyStates->GetNumAssoc());
	unsigned long long* ullpEdgePointer = m_pcGraph->EdgeListPntr();

	// For each edge in network
	for (unsigned long long iii = 0ULL; iii < m_pcGraph->GetNumEdges(); iii++)
	{
		// take node pair, which is in vector list
		// get vector of int for beliefs of both
		vector<int> vnFirstAgent = m_vpagentAgentList[ullpEdgePointer[2ULL*iii]]->GetIntBelief();
		vector<int> vnSecondAgent = m_vpagentAgentList[ullpEdgePointer[2ULL*iii + 1ULL]]->GetIntBelief();

		// calculate the inner_product and add to sum
		llGlobalEnergy -= static_cast<long long>(inner_product(vnFirstAgent.begin(), vnFirstAgent.end(), vnSecondAgent.begin(), 0LL));
	}

	return llGlobalEnergy - llMinGlobalEnergy;
}

double BNetModel::CalcAvgIndivEnergy()
{
	// sum variable double
	double dEnergySum = 0.0;
	// For each agent in network
	for (vector<Agent*>::iterator itIII = m_vpagentAgentList.begin(); itIII != m_vpagentAgentList.end(); ++itIII)
	{
		// return energy of agent and add to total
		dEnergySum += (*itIII)->GetEnergy();
	}

	// return total / #agents
	return dEnergySum / m_vpagentAgentList.size();
}

double BNetModel::CalcTotalEnergy()
{
	// Determine individual energies
	double dTotalEnergy = 0.0;
	for (vector<Agent*>::iterator itIII = m_vpagentAgentList.begin(); itIII != m_vpagentAgentList.end(); ++itIII)
	{
		dTotalEnergy += static_cast<long double>((*itIII)->GetEnergy());
	}

	// Adjust for J
	dTotalEnergy *= m_rcBNetParameters.GetJ();

	// Determine global energy
	dTotalEnergy += static_cast<double>(CalcGlobalEnergy()) * m_rcBNetParameters.GetI();

	return dTotalEnergy;
}
