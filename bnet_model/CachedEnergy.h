
#ifndef CACHEDENERGY_H_
#define CACHEDENERGY_H_

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <ctime>
#include <boost/unordered_map.hpp>
#include "Resources.h"

class CachedEnergy
{
private:

	unsigned int m_unNumAssociations;
	unsigned long m_ulNumEnergies;
	/*
	 * Note: Do not reference or point to the addresses of elements
	 * of the vector<bool> as this will cause errors.
	 *
	 *    bool* pb = &vb[1]; // conversion error - do not use
	 *    bool& refb = vb[1];   // conversion error - do not use
	 */
	std::vector<std::vector<bool> > m_vvbListOfBeliefSys;
	std::vector<double>	m_vdListOfEnergies;
	boost::unordered_map<std::vector<bool>, unsigned long > *m_pumapBeliefSysLookUp;
	boost::unordered_map<unsigned long, std::vector<int>* > *m_pumapBeliefVect;

	// Random belief selector
	BoostRandomInt<unsigned long> *m_crngRandomBeliefSelector;
	// Random state selector
	BoostRandomInt<unsigned int> *m_crngRandomStateSelector;

	CachedEnergy();
	std::vector<int>* GetIntBelief(std::vector<bool> &rvbBelief);

public:
	CachedEnergy(std::string sEnergyFile, unsigned int unSeed=static_cast<unsigned int>(std::time(0)));
	~CachedEnergy();

	unsigned long RandomBelief() const; // labeled as const because it doesn't change anything in the class
	unsigned int RandomState() const;
	std::vector<bool>& GetBeliefSys(unsigned long ulBeliefIndex);
	std::vector<int>* GetBeliefVector(unsigned long ulBeliefIndex);
	bool GetStateValue(unsigned long ulBeliefIndex, unsigned int unStateIndex) const;
	double GetEnergy(unsigned long ulBeliefIndex) const;
	unsigned long Search(std::vector<bool> &key) const;
	unsigned long GetNeighborState(unsigned long ulBeliefIndex, unsigned int unStateToFlip);
	unsigned int GetNumAssoc() const;

	// Will probably need a seed reset function for bnetmodel to call
	void KeyValueCheck();
	// Operator overloads
	friend std::ostream& operator<< (std::ostream &rout, CachedEnergy &rcCachedEnergy);

};


#endif /* CACHEDENERGY_H_ */
