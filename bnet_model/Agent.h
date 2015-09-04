#ifndef AGENT_H_
#define AGENT_H_

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>
#include "CachedEnergy.h"

const double EPSILON = 0.00001;

struct Emitor
{
	unsigned int state;
	bool value;
};

class Agent
{
private:
	CachedEnergy *m_pcEnergyStates;
	double *m_pdI;
	double *m_pdT;
	double *m_pdJ;
	BoostRandomReal<double> *m_pcReal0to1;
	// Index that points to the state in the s_cEnergyStates class (vector)
	unsigned long m_ulBelief;
	// Pointer that points to the vector<int> value in EnergyStates map
	std::vector<int> *m_pvBeliefVector;
	// Pointer to a vector of pointers to neighboring Agents
	std::vector<Agent*> *m_pvagentNeighbors;

public:

	// Copy constructor
	Agent(const Agent &ragentOther);
	// Randomized constructor
	Agent(CachedEnergy *pcEnergyStates, double *pdI, double *pdT, double *pdJ, BoostRandomReal<double> *pcReal0to1);
	// Construct for specific initial state
	Agent(CachedEnergy *pcEnergyStates, double *pdI, double *pdT, double *pdJ, BoostRandomReal<double> *pcReal0to1, unsigned long ulBelief);

	~Agent();

	void Input(Emitor &eRecieved, double &rdTotalEnergy);
	void Emit(Emitor &reEmit);

	// get functions
	unsigned long GetBeliefIndex() const;
	std::vector<bool> &GetBoolBelief();
	double GetEnergy() const;
	std::vector<Agent*> * GetNeighborVector() const;
	CachedEnergy* GetPrivateEnergyStates() const;
	double* GetPrivateI() const;
	double* GetPrivateT() const;
	double* GetPrivateJ() const;
	BoostRandomReal<double>* GetPrivateRNG() const;
	unsigned long GetPrivateBelief() const;
	std::vector<int> GetIntBelief();
	std::vector<int>* GetBeliefVectorPointer() const;

	// set function
	void SetBelief();
	void SetBelief(unsigned long ulNewBelief);
	void SetNeighbors(std::vector<Agent*> *pvpagentNeighbors);
	
	// Operator overloads
	friend std::ostream& operator<< (std::ostream &rout, Agent &rcAgent);
	Agent& operator=(const Agent& ragentOther);

};

// == operator overload for class (for copying agents)
inline bool operator==(const Agent& lhs, const Agent& rhs) { return lhs.GetBeliefIndex() == rhs.GetBeliefIndex(); };

#endif /* AGENT_H_ */
