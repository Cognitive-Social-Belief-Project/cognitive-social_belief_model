
#ifndef EDGEGRAPH_H_
#define EDGEGRAPH_H_

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include "Resources.h"
#include <boost/unordered_map.hpp>

class EdgeGraph
{
private:
	unsigned long long m_ullNumNodes;
	unsigned long long m_ullNumEdges;
	unsigned long long m_ullArrayLength;
	unsigned long long *m_aullEdgeList;
//	boost::unordered_map<unsigned long long, std::vector<unsigned long long>* > *m_phashNodeNeighbors;

	// Random number generator for selecting edges
	BoostRandomInt<unsigned long long> *m_prngEdgeSelector;
	//BoostRandomInt<unsigned long long> *m_prngNodeSelector;

public:

	EdgeGraph(std::string sGraphFile, unsigned int unSeed=static_cast<unsigned int>(std::time(0)));
	~EdgeGraph();

	// Random functions
	unsigned long long RandomEdge() const;
	//unsigned long long RandomNode() const;

	// Return node from edge and pair information
	unsigned long long PickNode(unsigned long long ullEdge, unsigned unPair) const;

	// Generate Node neighbor has
	void CreateNeighborHash();

	// Return reference to vector of node neighbors
	std::vector<unsigned long long> &NodeNeighbors(unsigned long long ullNode);

	// Get functions
	const unsigned long long &GetNumNodes() const;
	const unsigned long long &GetNumEdges() const;
	const unsigned long long &GetArrayLen() const;
	unsigned long long *EdgeListPntr();

	// Operator overloads
	friend std::ostream& operator<< (std::ostream &rout, EdgeGraph &rcEdgeGraph);

};


#endif /* EDGEGRAPH_H_ */
