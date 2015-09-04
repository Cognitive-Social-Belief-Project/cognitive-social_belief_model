
#ifndef MODELPARAMETERS_H_
#define MODELPARAMETERS_H_

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>

struct Assignment{
	unsigned long ulNode;
	std::string sType;
	unsigned int unSubset;
};

class ModelParameters
{
private:
	// General information for loading in parameters and cached data
	std::string m_sParamFile;
	std::string m_sSizesFile;
	std::string m_sOutEnergyFile;
	std::string m_sOutputFile;
	double m_dI;
	double m_dT;
	double m_dJ;
	unsigned long m_ulNumPoints;
	unsigned long long m_ulNumIterations;
	unsigned int m_unRuns;
	std::string m_sEnergyFile;
	std::string m_sUpdateMethod;
	std::string m_sGraphFile;
	unsigned int m_unNumSets;
	unsigned long m_ulGraphSize;
	bool m_bSeed;
	unsigned int m_unSeed;
	bool m_bFinalState;
	bool m_bTimeSeries;
	bool m_bTSSizes;
	bool m_bAppendOut;
//	bool m_bPreProccessNeighors;

	// Initial belief system configuration
	bool m_bRandomBeliefs;
	unsigned int m_nNumBeliefTypes;
	std::vector<double> m_vdBeliefRatios;
	std::vector<std::vector<std::vector<bool> > > m_vvvbBeliefSets;

	// Zealot configuration
	bool m_bZealots;
	unsigned int m_nNumZealotTypes;
	std::vector<double> m_vdZealotRatios;
	std::vector<std::vector<std::vector<bool> > > m_vvvbZealotSets;

	// Assignment info
	bool m_bIsAssignment;
	std::vector<Assignment*> m_vpassNodeAssignments;

public:
	// Constructors
	ModelParameters();
	ModelParameters(std::string sParamFile);
	~ModelParameters();

	// Read parameter file and set values
	void GetParams();

	// Get values
	const std::string &GetOutFile() const;
	const std::string &GetSizesFile() const;
	const std::string &GetOutEnergyFile() const;
	double &GetI();
	double &GetT();
	double &GetJ();
	unsigned long GetNumPoints() const;
	unsigned long long GetNumIter() const;
	const unsigned int &GetNumRuns() const;
	const std::string &GetUpdateMethod() const;
	const std::string &GetEnergyFile() const;
	const std::string &GetGraphFile() const;
	const unsigned int &GetNumSets() const;
	bool IsRandom() const;
	bool IsZealots() const;
	bool IsSeed() const;
	bool IsPreproccessNeighbors() const;
	unsigned int GetSeed() const;
	std::vector<double>& GetBeliefRatios();
	std::vector<std::vector<std::vector<bool> > >& GetBeliefSets();
	std::vector<bool> & GetBeliefSetVector(unsigned int unSet, unsigned int unBeliefType);
	std::vector<double>& GetZealotRatios();
	std::vector<std::vector<std::vector<bool> > >& GetZealotSets();
	std::vector<bool>& GetZealotSetVector(unsigned int unSet, unsigned int unBeliefType);
	Assignment* GetNodeAssignment(unsigned long ulIndex);
	std::vector<Assignment*>& GetAssignmentList();
	bool IsAssigned() const;
	bool IsTimeSeries() const;
	bool IsFinalStates() const;
	bool IsTSSizes() const;
	bool IsAppendOut() const;

	// Set values
	void Set_I(const double dI);
	void Set_T(const double dT);
	void Set_J(const double dJ);
	void Set_EnergyFile(const std::string sFilename);
	void Set_UpdateMethod(const std::string sMethod);
	void Set_GraphFile(const std::string sFilename);
	void Set_NonRandomBeliefs(bool bChoice);
	void Set_Zealots(bool bChoice);

	// Operator overloads
	friend std::ostream& operator<< (std::ostream &rout, ModelParameters &rcModelParameters);
};


#endif /* MODELPARAMETERS_H_ */
