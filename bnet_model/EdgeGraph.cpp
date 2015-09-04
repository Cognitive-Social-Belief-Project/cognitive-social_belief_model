#include "EdgeGraph.h"

using namespace std;

EdgeGraph::EdgeGraph(string sGraphFile, unsigned int unSeed /*=static_cast<unsigned int>(std::time(0))*/)
{
	// Read file
	ifstream ifGraphFile(sGraphFile.c_str());

	if (!ifGraphFile)
	{
		cerr << "Could not read parameter file:" << sGraphFile << endl;
		exit(1);
	}

	// Buffer to retrieve lines
	string sBuff;
	stringstream ssLineOfFile;
	// Create stream and buffer for comma delineated columns in file.
	string sColumn;
	stringstream ssColumnOfLine;

	// Read the graph dimensions
	getline(ifGraphFile, sBuff);
	ssLineOfFile << sBuff;

	getline(ssLineOfFile, sColumn, ',');
	ssColumnOfLine << sColumn;
	ssColumnOfLine >> m_ullNumNodes;
	ssColumnOfLine.clear();

	getline(ssLineOfFile, sColumn, ',');
	ssColumnOfLine << sColumn;
	ssColumnOfLine >> m_ullNumEdges;
	ssColumnOfLine.clear();

	ssLineOfFile.clear();
	ssLineOfFile.str("");

	// Read edges into graph
	m_aullEdgeList = new unsigned long long[m_ullNumEdges*2ULL];
	m_ullArrayLength = m_ullNumEdges*2ULL;
	unsigned long long ulNode;
	// This produces a list that can be accessed using:
	// Row iii, Col jjj:
	// m_aullEdgeList[iii*2 + jjj]
	for (unsigned long long iii=0ULL; iii<m_ullNumEdges; iii++)
	{

		getline(ifGraphFile, sBuff);
		ssLineOfFile << sBuff;

		// Get first connecting node
		getline(ssLineOfFile, sColumn, ',');
		ssColumnOfLine << sColumn;
		ssColumnOfLine >> ulNode;
		ssColumnOfLine.clear();
		m_aullEdgeList[iii*2] = ulNode;

		// Get second connecting node
		getline(ssLineOfFile, sColumn, ',');
		ssColumnOfLine << sColumn;
		ssColumnOfLine >> ulNode;
		ssColumnOfLine.clear();
		m_aullEdgeList[iii*2 + 1] = ulNode;

		ssLineOfFile.clear();
		ssLineOfFile.str("");
		ssColumnOfLine.str("");

	}

	/*
	 * This reads in the edges from file
	 * To allocate: aulEdgeList = new unsigned long[width*height]
	 * to call: aulEdgeList[row*width + column]
	 * to clean: delete [] aulEdgeList
	 */

	// Create random number generator
	unsigned long long ullMin = 0ULL;
	m_prngEdgeSelector = new BoostRandomInt<unsigned long long>(ullMin, m_ullNumEdges-1, unSeed);
	//m_prngNodeSelector = new BoostRandomInt<unsigned long long>(ullMin, m_ullNumEdges-1, unSeed);

}

EdgeGraph::~EdgeGraph()
{
	delete m_prngEdgeSelector;
	delete [] m_aullEdgeList;
}

unsigned long long EdgeGraph::RandomEdge() const
{
	return (*m_prngEdgeSelector)();
}

//unsigned long long EdgeGraph::RandomNode() const
//{
//	return (*m_prngNodeSelector)();
//}

//void EdgeGraph::CreateNeighborHash()
//{
//
//	// Need to check if vectors are created on heap or stack. (determine if we need new)
//
//	// Create new hash
//	m_phashNodeNeighbors = new boost::unordered_map<unsigned long long, vector<unsigned long long>* >;
//
//	// For edge in network -> loop
//	for (unsigned long long iii = 0ULL; iii < m_ullNumEdges; iii++)
//	{
//		// if node in hash push neighbor into vector
//		if (m_phashNodeNeighbors->find(iii) != m_phashNodeNeighbors->end())
//		{
//			(*m_phashNodeNeighbors)[iii]->push_back(PickNode(iii, 1U));
//		}
//		// if node not in hash create vector and push neighbor into vector
//		else
//		{
//			vector<unsigned long long> *pvullNeighors = new vector<unsigned long long>;
//			pvullNeighors->push_back(PickNode(iii, 1U));
//			m_phashNodeNeighbors->insert(make_pair(iii, pvullNeighors));
//		}
//	}
//
//}

//vector<unsigned long long> &EdgeGraph::NodeNeighbors(unsigned long long ullNode)
//{
//	return *(*m_phashNodeNeighbors)[ullNode];
//}

const unsigned long long &EdgeGraph::GetNumNodes() const
{
	return m_ullNumNodes;
}

const unsigned long long &EdgeGraph::GetNumEdges() const
{
	return m_ullNumEdges;
}

const unsigned long long &EdgeGraph::GetArrayLen() const
{
	return m_ullArrayLength;
}

unsigned long long *EdgeGraph::EdgeListPntr()
{
	return m_aullEdgeList;
}

unsigned long long EdgeGraph::PickNode(unsigned long long ullEdge, unsigned int unPair) const
{
	return m_aullEdgeList[2ULL * ullEdge + unPair];
}

//------------------------------------------------------------------------------
//
//	Friend functions
//
//------------------------------------------------------------------------------
ostream& operator<< (ostream &rout, EdgeGraph &rcEdgeGraph)
{

	string sOutString = "";
	stringstream stream; // for converting types to strings
	string sConverted = "";

	stream << rcEdgeGraph.m_ullNumNodes;
	stream >> sConverted;
	stream.clear();
	sOutString += "Number of Nodes: " + sConverted + "\n";

	stream << rcEdgeGraph.m_ullNumEdges;
	stream >> sConverted;
	stream.clear();
	sOutString += "Number of Edges: " + sConverted + "\n";

	// Print edge list
	for (unsigned long long iii=0ULL; iii < rcEdgeGraph.m_ullNumEdges; iii++)
	{
		stream << rcEdgeGraph.m_aullEdgeList[iii*2];
		stream >> sConverted;
		stream.clear();
		sOutString += sConverted + ", ";

		stream << rcEdgeGraph.m_aullEdgeList[iii*2 + 1];
		stream >> sConverted;
		stream.clear();
		sOutString += sConverted + "\n";
	}

	rout << sOutString;
	return rout;
}
