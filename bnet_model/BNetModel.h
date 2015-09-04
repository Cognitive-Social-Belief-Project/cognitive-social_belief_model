#ifndef BNETMODEL_H_
#define BNETMODEL_H_

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>
#include "ModelParameters.h"
#include "EdgeGraph.h"
#include "Agent.h"
#include "CachedEnergy.h"
#include "Resources.h"
#include <boost/unordered_map.hpp>

class BNetModel
{
private:
	ModelParameters &m_rcBNetParameters;
	CachedEnergy *m_pcEnergyStates;
	EdgeGraph *m_pcGraph;

	// vector of agents
	std::vector<Agent*> m_vpagentAgentList;
	// Vector of neighbor pointers
	boost::unordered_map<unsigned long long, std::vector<Agent*>* >* m_pumapAgentNeighbors;

	// Random real number generators
	BoostRandomReal<double> *m_prngReal0to1;
	BoostRandomInt<unsigned int> *m_prngUInt0or1;

	// emitor to be modified by agents
	Emitor m_eEmit;

	// Hamiltonian (total energy)
	double m_dTotalEnergy;

	// Immutable handler: Start at index 0 and end at M (excluding)
	// m_ullFirstBeliefIndex is the index after the last zealot (for looping purposes)
	unsigned long long m_ullFirstBeliefIndex;

	void CreateNeighborHash();
	void AssignAgentNeighbors();

public:

	BNetModel(ModelParameters &rcParameters);
	~BNetModel();

	void InitialPopulate();
	void RePopulate(unsigned int unSetNumber);
	void Iterate();
	void RunSet(unsigned int unSetID);
	void Engage();

	void UniqueAgents(boost::unordered_map<unsigned long, unsigned long long >* pagentHash);
	void AppendCurrentStates(const std::string &rsFilename);
	void AppendSizesToFile(const std::string& rsSizesFile, unsigned int unSetID);
//	void AppendTopNSizesToFile(const std::string& rsSizesFile, unsigned int unSetID, unsigned long ulTopN);
	void AppendEnergyToFile(const std::string& rsEnergyFile);

	double CalcGlobalEnergy();
	double CalcAvgIndivEnergy();
	double CalcTotalEnergy();

	std::string ConvertBeliefToString(unsigned long ulBelief);
};

#endif /* BNETMODEL_H_ */
