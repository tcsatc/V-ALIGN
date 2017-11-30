// Author: Shreyansh Chhajer 
// Tata Consultancy Services Ltd
using namespace std;


#include <bits/stdc++.h>
#include <stdio.h>
#include <string.h>

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <list>
#include <queue>
#include <set>
#include <stack>
#include <map>
#include <utility>


#define ff first
#define ss second

class graphClass{
public:
	// Negative Infinity
	float INF = -1.0 * FLT_MAX;
	// Adjacency List
	vector<int> *adj;
	// Incoming Edges
	vector<int> *inAdj;
	// Label of nodes, i.e., sequence of bases
	string *label;
	int E;
	// Number of edges in the virtual graph								
	int Ea;
	// Number of nodes in the gfa file
	int V = 0;															
	// Number of nodes in the genome graph
	int N;
	// Number of nodes in the Graph with each node having label size of 1
	int Va;
	// node Shortest Vertex Weighted Distance matrix
	int **pathWt;
	// Path Wt in reverse
	int **pathWtR;

	// Loop Cost
	int *loopCost;

	// Mapping string index to int index
	map<string, int> mapIx;
	// Reverse Map
	string *revIx;

	// Alignment Strings
	string resGraph, resSeq, exSeq;

	// Sequence Alignment Parameters
	float **M, **Q, **R;
	// Stores the flags, i.e., stores from where a certain value has arrived
	// 0 -> Q
	// 1 -> R
	// 2 -> M
	bool **fM, **fQ, **fR;
	pair<int, int> **afvs;
	// Stores the index to move to next
	int **iM, **iR;
	// For mapping a node index to it's corresponding index in topological sorted order
	vector<int> *ixTS;

	vector<int> **shortestPath;

	// Store the Aligned Path for the given sequence
	vector<pair<int, int> > path, tpPath;

	// Stores the minimum feedback vertex set
	vector<int> Vc;

	// Topological Sorting of V - V'
	vector<pair<int, int> > topOrder;

	// Boolean vector for checking the status of a vertex
	vector<int> vis;

	// Gap Penalty Parameters
	float gapOpen = 10.0;
	float gapExt = 0.5;

	// Score Parameters
	vector<vector<float> > sc;
	vector<int> chMap;

	// Sets up the parameters for computing topological order
	void computeTopologicalOrder();
	
	// DFS like implementation of topological sorting
	void topologicalSort(int v);

	graphClass(int type, string inputFile, string mfvsName, string scoreMatrixName);

	~graphClass();

	// Delta Function
	inline float delta(int k);

	// Score Function
	inline float score(char w, char xj);

	// Computes shortest vertex weighted distance between any two nodes 
	void computeVertexWeightedDistance();

	// Shortest Edge Distance from ith base of node u to jth base of node v
	int edgeDist(int u, int i, int v, int j);

	// Returns the minimum feedback vertex set as computed using it's external implementation
	vector<int> mfvsExternalImplementation();

	// Computing Matrices M, R, Q
	void computeMatricesMRQ(string seq, bool globalAlign, pair<float, float>);

	// Aligning a sequence to the directed graph
	void alignSequence(string seq, bool globalAlign, pair<float, float>, fstream&, fstream&);

	// For displaying shortest distance matrices, minimum feedback vertices, etc in a tabular format
	// void printGraph(fstream &debug);

	// For displaying M, R, Q Matrices in a tabular format
	void printMRQ(string x, fstream &debug);

	void printVirtualGraph(fstream &debug);

	// To create dot file for better visualization of the path aligned with the sequence
	void writeDotFile(fstream &debug);

	// Computes reverse complement of a string label
	string revComp(string);

	// Writes a dot file to visualize the gfa file
	void visualizeGFA(fstream &debug);

	// Writes a dot file to visualize mfvs node
	void mfvsVisualize(fstream &dotFile, int depth);

	// returns the complement nucleotide for the given nucleotide
	char revChar(char ch);
};
