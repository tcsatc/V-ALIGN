#include <bits/stdc++.h>
// For MFVS External Implementation
#include "graph.h"
#include "utils.h"

#include "graphClass.h"

// #include "graphClassMemberFunctions.h"

using namespace std;

#define ff first
#define ss second
timespec st,en;
double runtime;

int main(int argc, char** argv){

	int _;
	// clock_gettime(CLOCK_MONOTONIC,&st);
	
	if(argc == 1){
		printf("usage: ./genfvs graphfile\n");
		return 0;
	}

	string adjList;
	string mfvsFile = "";

	if(argc != 2){
		printf("Incorrect usage\n");
		printf("usage: ./genfvs graphfile\n");
		exit(1);
	}
	adjList = argv[1];

	mfvsFile = adjList;
	mfvsFile.append(".fvs");


	int type = 1;
	string suff = adjList.substr(adjList.length() - 3, 3);
	//cout << "Suff : " << suff << endl;
	if(suff == "gfa") {
		type = 0;
	}else if(suff == "dot"){
		type = 2;
	}

	graphClass G(type, adjList);

	vector<string> fvs = G.mfvsExternalImplementation(type);

	ofstream of;
	of.open(mfvsFile);
	for(int i=0; i < fvs.size(); ++i){
		of << fvs[i] << endl;
	}
	of.close();
	
	cout << "Feedback vertex set in file " << mfvsFile << endl;


	return 0;
}
