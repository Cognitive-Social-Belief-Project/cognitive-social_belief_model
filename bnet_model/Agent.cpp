#include <numeric>
#include "Agent.h"
#include "Resources.h"

using namespace std;

Agent::Agent(const Agent &ragentOther)
{
	*m_pcEnergyStates = *(ragentOther.GetPrivateEnergyStates());
	*m_pdI = *(ragentOther.GetPrivateI());
	*m_pdT = *(ragentOther.GetPrivateT());
	*m_pdJ = *(ragentOther.GetPrivateJ());
	*m_pcReal0to1 = *(ragentOther.GetPrivateRNG());
	m_ulBelief = ragentOther.GetPrivateBelief();
	*m_pvBeliefVector = *(ragentOther.GetBeliefVectorPointer());
	*m_pvagentNeighbors = *(ragentOther.GetNeighborVector());
}

Agent::Agent(CachedEnergy *pcEnergyStates, double *pdI, double *pdT, double *pdJ, BoostRandomReal<double> *pcReal0to1)
	: m_pcEnergyStates(pcEnergyStates), m_pdI(pdI), m_pdT(pdT), m_pdJ(pdJ), m_pcReal0to1(pcReal0to1), m_pvagentNeighbors(NULL)
{
	m_ulBelief = m_pcEnergyStates->RandomBelief();
	m_pvBeliefVector = m_pcEnergyStates->GetBeliefVector(m_ulBelief);
}

Agent::Agent(CachedEnergy *pcEnergyStates, double *pdI, double *pdT, double *pdJ, BoostRandomReal<double> *pcReal0to1, unsigned long ulBelief)
	: m_pcEnergyStates(pcEnergyStates), m_pdI(pdI), m_pdT(pdT), m_pdJ(pdJ), m_pcReal0to1(pcReal0to1), m_pvagentNeighbors(NULL)
{
	m_ulBelief = ulBelief;
	m_pvBeliefVector = m_pcEnergyStates->GetBeliefVector(m_ulBelief);
}

Agent::~Agent()
{
	// Nothing new
}

void Agent::Emit(Emitor &reEmit)
{
	// Change emitor state
	reEmit.state = m_pcEnergyStates->RandomState(); // Picks a location in the association vector
	reEmit.value = m_pcEnergyStates->GetStateValue(m_ulBelief, reEmit.state);
}

void Agent::Input(Emitor &reRecieved, double &rdTotalEnergy)
{

	// If the recieved state is different that the current state, apply acceptance rules
	if (m_pcEnergyStates->GetStateValue(m_ulBelief, reRecieved.state) != reRecieved.value)
	{
		// Set initial energy sum
		double dCurrentEnergy = 0.0;
		double dFutureEnergy = 0.0;

		// Generate neighbor belief system
		unsigned long ulNewBelief = m_pcEnergyStates->GetNeighborState(m_ulBelief, reRecieved.state);
		// Test -< get random belief and check
//		unsigned long ulNewBelief = m_pcEnergyStates->RandomBelief(); // test remove
		// test get flipped version (only for binary)
//		unsigned long ulNewBelief = 1;
//		if (m_ulBelief == 1UL)
//		{
//			ulNewBelief = 0;
//		}

		// Determine collective energy
		/*
		 * Can optimize by exchanging inner_product for a sum of vectors and then 1 inner_product
		 * Requires another iterator though to hold the temporary vector sums. Or maybe some other elaborate design that can avoid it
		 */
		vector<int>* pvnFutureBelief = m_pcEnergyStates->GetBeliefVector(ulNewBelief); // 1, -1 is the whole vector for ising
//		cout << endl;//test
		for (vector<Agent*>::iterator itIII = m_pvagentNeighbors->begin(); itIII != m_pvagentNeighbors->end(); ++itIII)
		{

			double dCurrentInner = inner_product(m_pvBeliefVector->begin(), m_pvBeliefVector->end(), (*itIII)->m_pvBeliefVector->begin(), 0.0);
			dCurrentEnergy += dCurrentInner;
			double dFutureInner = inner_product(pvnFutureBelief->begin(), pvnFutureBelief->end(), (*itIII)->m_pvBeliefVector->begin(), 0.0);
			dFutureEnergy += dFutureInner;

//			// test start
//			cout << "Current: (" << m_ulBelief << ", " << (*itIII)->m_ulBelief << ") " << " Prod: " << dCurrentInner << "  \t" <<
//					"<" << (*m_pvBeliefVector)[0] << "> " << "<" << (*(*itIII)->m_pvBeliefVector)[0] << ">" << endl;
//			cout << "Future: (" << ulNewBelief << ", " << (*itIII)->m_ulBelief << ") " << " Prod: " << dFutureInner << "  \t" <<
//					"<" << (*pvnFutureBelief)[0] << "> " << "<" << (*(*itIII)->m_pvBeliefVector)[0] << ">" << endl;
//			// test end
		}

		// Determine local energy
		/*
 	 	 * H = (- Sum_belief_triangles) + I * difference_inner_products
		 */
		dCurrentEnergy = *m_pdJ * m_pcEnergyStates->GetEnergy(m_ulBelief) - dCurrentEnergy * *m_pdI;
		dFutureEnergy = *m_pdJ * m_pcEnergyStates->GetEnergy(ulNewBelief) - dFutureEnergy * *m_pdI;

//		cout << m_ulBelief << "-Current: " << dCurrentEnergy << "\t" << ulNewBelief << "-Future: " << dFutureEnergy << "\t"; //test

		// start test
//		for (vector<Agent*>::iterator itIII = m_pvagentNeighbors->begin(); itIII != m_pvagentNeighbors->end(); ++itIII)
//		{
//			cout << (*itIII)->m_ulBelief << " ";
//		}
//		cout << endl; //test
		// test end

		// Accept state if it has lower energy
		if (dFutureEnergy < dCurrentEnergy)
		{
			SetBelief(ulNewBelief);
			rdTotalEnergy += dFutureEnergy - dCurrentEnergy;
//			accept++;//test
//			cout << " - lower (accept)\t" << endl; //test
		}
		else
		{
			double dAcceptanceProb = Boltzmann( dFutureEnergy - dCurrentEnergy, *m_pdT);
			double dRndVar = (*m_pcReal0to1)();
			if ( dRndVar < dAcceptanceProb)
			{
				SetBelief(ulNewBelief);
				rdTotalEnergy += dFutureEnergy - dCurrentEnergy;
//				flip++;//test
//				cout << " - flip " << endl; //test
			}
//			else
//			{
//				stay++;//test
//				cout << " - stay" << endl; //test
//			}
		}
//		SetBelief(1);//test
	}
}

unsigned long Agent::GetBeliefIndex() const
{
	return m_ulBelief;
}

vector<bool> &Agent::GetBoolBelief()
{
	return m_pcEnergyStates->GetBeliefSys(m_ulBelief);
}

double Agent::GetEnergy() const
{
	return m_pcEnergyStates->GetEnergy(m_ulBelief);
}

std::vector<Agent*> * Agent::GetNeighborVector() const
{
	return m_pvagentNeighbors;
}

void Agent::SetBelief()
{
	// Sets random belief
	m_ulBelief = m_pcEnergyStates->RandomBelief();
	// Changes belief vector pointer
	m_pvBeliefVector = m_pcEnergyStates->GetBeliefVector(m_ulBelief);
}

void Agent::SetNeighbors(vector<Agent*> *pvpagentNeighbors)
{
	m_pvagentNeighbors = pvpagentNeighbors;
}

void Agent::SetBelief(unsigned long ulNewBelief)
{
	m_ulBelief = ulNewBelief;
	// Change belief vector pointer
	m_pvBeliefVector = m_pcEnergyStates->GetBeliefVector(ulNewBelief);
}

CachedEnergy* Agent::GetPrivateEnergyStates() const
{
	return m_pcEnergyStates;
}

double* Agent::GetPrivateI() const
{
	return m_pdI;
}

double* Agent::GetPrivateT() const
{
	return m_pdT;
}

double* Agent::GetPrivateJ() const
{
	return m_pdJ;
}

BoostRandomReal<double>* Agent::GetPrivateRNG() const
{
	return m_pcReal0to1;
}

unsigned long Agent::GetPrivateBelief() const
{
	return m_ulBelief;
}

vector<int>* Agent::GetBeliefVectorPointer() const
{
	return m_pvBeliefVector;
}

vector<int> Agent::GetIntBelief()
{
	vector<int> vnBelief;
	vector<bool> &rvbBelief = GetBoolBelief();
	for (vector<bool>::iterator itIII = rvbBelief.begin(); itIII != rvbBelief.end(); ++itIII)
	{
		if (*itIII)
		{
			vnBelief.push_back(1);
		}
		else
		{
			vnBelief.push_back(-1);
		}
	}

	return vnBelief;
}

Agent& Agent::operator=(const Agent& ragentOther)
{
	if (this != &ragentOther)
	{
		*m_pcEnergyStates = *(ragentOther.GetPrivateEnergyStates());
		*m_pdI = *(ragentOther.GetPrivateI());
		*m_pdT = *(ragentOther.GetPrivateT());
		*m_pcReal0to1 = *(ragentOther.GetPrivateRNG());
		m_ulBelief = ragentOther.GetPrivateBelief();
	}

	return *this;
}

//------------------------------------------------------------------------------
//
//	Friend/Support functions
//
//------------------------------------------------------------------------------
ostream& operator<< (ostream &rout, Agent &rcAgent)
{
	string sOutString = "";
	stringstream stream; // for converting types to strings
	string sConverted = "";

	stream << rcAgent.m_ulBelief;
	stream >> sConverted;
	stream.clear();
	sOutString += "Agent Index: " + sConverted + "\n";

	rout << sOutString;
	return rout;
}
