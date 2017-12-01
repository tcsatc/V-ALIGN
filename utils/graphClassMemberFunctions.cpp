#include <bits/stdc++.h>
#include "graphClass.h"
#include "graph.h"

// returns the complement nucleotide for the given nucleotide
char graphClass::revChar(char ch){
	if(ch == 'A')
		return 'T';
	else if(ch == 'T')
		return 'A';
	else if(ch == 'G')
		return 'C';
	else if(ch == 'C')
		return 'G';
	else
		return ch;
}

// Returns the reverse complement of the given label
string graphClass::revComp(string lbl){
	reverse(lbl.begin(), lbl.end());
	for(int i = 0; i < lbl.size(); i++)
		lbl[i] = revChar(lbl[i]);

	return lbl;
}


graphClass::graphClass(int type, string inputFile){
	// GFA Input
	if(type == 0){
		fstream gfaFile(inputFile, fstream::in);

		if(!gfaFile){
			perror("The following error occured");
			exit(1);
		}

		string line, _, ix, lbl;

		V = 0, Va = 0;

		while(getline(gfaFile, line)){
			if(line == "/*"){
				
				while(getline(gfaFile, line)){
					if(line == "*/")
						break;
				}

				continue;
			}
			istringstream iss(line);
			iss >> _;
			if(_ == "S")
				V += 1;
		}

		// Sets up the objects
		adj = new vector<int>[(V << 1) + 2];
		inAdj = new vector<int>[(V << 1) + 2];
		revIx = new string[(V << 1) + 2];

		gfaFile.clear();
		gfaFile.seekg(0, ios::beg);

		N = 1, Va = 0;

		// Dummy Node Label
		//label[N - 1] = "0";

		while(getline(gfaFile, line)){
			if(line == "/*"){
				
				while(getline(gfaFile, line)){
					if(line == "*/")
						break;
				}

				continue;
			}
			istringstream iss(line);
			iss >> _;
			if(_ != "S")
				continue;
			iss >> ix >> lbl;
			mapIx[ix] = N; 
			revIx[N] = ix;
			revIx[V + N] = ix;

			//cout << "Created two new vertices " << N << " " << (V + N) << endl;
			++N;
		}

		N = (V << 1);

		E = 0;
		gfaFile.clear();
		gfaFile.seekg(0, ios::beg);

		while(getline(gfaFile, line)){
			if(line == "/*"){
				
				while(getline(gfaFile, line)){
					if(line == "*/")
						break;
				}

				continue;
			}
			istringstream iss(line);
			iss >> _;
			if(_ != "L")
				continue;
			string u, v, t1, t2;
			iss >> u >> t1 >> v >> t2;
			// urc - reverse complement of u
			// vrc - reverse complement of v

			if(t1 == "+" && t2 == "+"){
				adj[mapIx[u]].push_back(mapIx[v]);		// u -> v 
				//cout << "Edge from " << mapIx[u] << " to " << mapIx[v] << endl;
				inAdj[mapIx[v]].push_back(mapIx[u]);
			}
			else if(t1 == "+" && t2 == "-"){
				adj[mapIx[u]].push_back(mapIx[v] + V);	// u -> vrc
				//cout << "Edge from " << mapIx[u] << " to " << (mapIx[v] + V) << endl;
				inAdj[mapIx[v] + V].push_back(mapIx[u]);
			}
			else if(t1 == "-" && t2 == "+"){
				adj[mapIx[u] + V].push_back(mapIx[v]);			// urc -> v
				//cout << "Edge from " << (mapIx[u] + V)  << " to " << mapIx[v] << endl;
				inAdj[mapIx[v]].push_back(mapIx[u] + V);   
			}
			else{
				adj[mapIx[u] + V].push_back(mapIx[v] + V);		// urc -> vrc
				//cout << "Edge from " << (mapIx[u] + V) << " to " << (mapIx[v] + V) << endl;
				inAdj[mapIx[v] + V].push_back(mapIx[u] + V);
			}

			E += 1;
		}

		for(int i = 1; i <= N; i++){
			if(inAdj[i].empty()){
				adj[0].push_back(i);
				inAdj[i].push_back(0);
			}
		}

		gfaFile.close();
	}
	else if(type == 1){		// Adjacency List Input
		fstream inputGraph;

		inputGraph.open(inputFile, fstream::in);

		if(!inputGraph){
			perror("The following error occured");
			return;
		}

		string line;

		N = 1, Va = 0;

		while(getline(inputGraph, line)){
			if(line == "/*"){
				
				while(getline(inputGraph, line)){
					if(line == "*/")
						break;
				}

				continue;
			}
			string ix, lbl;
			istringstream iss(line);
			iss >> ix >> lbl;
			if(lbl == "-")
				continue;
			Va += lbl.size();
			N += 1;
		} 

		N -= 1;

		// Sets up the objects
		adj = new vector<int>[N + 1];
		inAdj = new vector<int>[N + 1];
		
		ixTS = new vector<int>[N + 1];
		revIx = new string[N + 1];

		inputGraph.clear();
		inputGraph.seekg(0, ios::beg);

		N = 1;
		Va = 0;
		

		while(getline(inputGraph, line)){
			if(line == "/*"){
				
				while(getline(inputGraph, line)){
					if(line == "*/")
						break;
				}

				continue;
			}
			string ix, lbl;
			istringstream iss(line);
			iss >> ix >> lbl;
			if(lbl == "-")
				continue;
			mapIx[ix] = N++;
			revIx[N - 1] = ix;
		}

		N -= 1;
		
		inputGraph.clear();
		inputGraph.seekg(0, ios::beg);

		E = 0;

		while(getline(inputGraph, line)){
			if(line == "/*"){
				
				while(getline(inputGraph, line)){
					if(line == "*/")
						break;
				}

				continue;
			}
			string u, tmp, v;
			istringstream iss(line);
			iss >> u >> tmp;
			if(tmp != "-")
				continue;
			iss >> v;
			adj[mapIx[u]].push_back(mapIx[v]);
			inAdj[mapIx[v]].push_back(mapIx[u]);		
			E += 1;
		}

		for(int i = 1; i <= N; i++){
			if(inAdj[i].empty()){
				adj[0].push_back(i);
				inAdj[i].push_back(0);
			}
		}

		V = N;

		inputGraph.close();
	}
	else{		// Dot File Input
		fstream dotFile;

		dotFile.open(inputFile, fstream::in);

		if(!dotFile){
			perror("The following error occured");
			return;
		}

		string line;
		N = 0, Va = 0;

		// Assuming there are no repeating node id in the dot file

		while(getline(dotFile, line)){
			if(line == "/*"){
				
				while(getline(dotFile, line)){
					if(line == "*/")
						break;
				}

				continue;
			}
			string ix, lbl;
			istringstream iss(line);
			iss >> ix >> lbl;
			if(lbl.empty() || lbl[0] != '[')
				continue;
			N += 1;
		}

		// Sets up the objects
		adj = new vector<int>[N + 1];
		inAdj = new vector<int>[N + 1];

		
		revIx = new string[N + 1];

		dotFile.clear();
		dotFile.seekg(0, ios::beg);

		N = 1;
		Va = 0;
		
		// Dummy Node Label
		//label[N - 1] = "0";

		map<string, string> dotMap;

		while(getline(dotFile, line)){
			if(line == "/*"){
				
				while(getline(dotFile, line)){
					if(line == "*/")
						break;
				}

				continue;
			}
			string ix, lbl, id, _;
			istringstream iss(line);
			iss >> ix >> lbl;
			if(lbl.empty() || lbl[0] != '[')
				continue;
			stringstream ss(lbl);
			getline(ss, _, '"');
			getline(ss, lbl, '"');

			ss.str(lbl);
			getline(ss, _, ':');
			getline(ss, lbl, ':');
			if(_ != lbl)
				id = _;
			else
				id = ix;

			dotMap[ix] = id;
			mapIx[id] = N++;
			revIx[N - 1] = id;
		}

		N -= 1;
		
		dotFile.clear();
		dotFile.seekg(0, ios::beg);

		E = 0;

		while(getline(dotFile, line)){
			if(line == "/*"){
				
				while(getline(dotFile, line)){
					if(line == "*/")
						break;
				}

				continue;
			}
			string u, tmp, v;
			istringstream iss(line);
			iss >> u >> tmp;
			if(tmp != "->")
				continue;
			iss >> v;
			adj[mapIx[dotMap[u]]].push_back(mapIx[dotMap[v]]);
			inAdj[mapIx[dotMap[v]]].push_back(mapIx[dotMap[u]]);		
			E += 1;
		}

		for(int i = 1; i <= N; i++){
			if(inAdj[i].empty()){
				adj[0].push_back(i);
				inAdj[i].push_back(0);
			}
		}

		V = N;

		dotFile.close();
	}

}



// Computing shortest vertex weighted path considering each node as source using Dijkstra's Algorithm (Only for feedback Vertex Set)
void graphClass::computeVertexWeightedDistance(){
	for(int i = 1; i <= Vc.size(); i++)
		for(int j = 0; j <= N; j++)
			pathWt[i][j] = INT_MAX;

	for(int i = 1; i <= Vc.size(); i++)
		for(int j = 0; j <= N; j++)
			pathWtR[i][j] = INT_MAX;

	for(int i = 1; i <= Vc.size(); i++)
		loopCost[i] = INT_MAX;

	// O(N)
	vector<pair<int, int> > seq(N + 1);

	// O(N)
	for(int i = 0; i <= N; i++)
		seq[i] = make_pair(INT_MAX, i);
	
	// O(NElog(N) + N*Nlog(N))
	for(int t = 0; t < Vc.size(); t++){
		int i = Vc[t];

		// O(N)
		set<pair<int, int> > st(seq.begin(), seq.end());

		// O(log(N))
		st.erase(make_pair(INT_MAX, i));
		st.insert(make_pair(0, i));
		pathWt[t + 1][i] = (int)label[i].size();
		shortestPath[t + 1][i].push_back(i);

		// O(N)
		// Visited array to check if the shortest distance of a node has been calculated or not
		vector<bool> vis(N + 1, 0);

		// O(Elog(N) + Nlog(N))
		while(!st.empty()){

			pair<int, int> top = *(st.begin());
			st.erase(st.begin());
			vis[top.ss] = 1;
			if(top.ff == INT_MAX)
				break;

			for(int j = 0; j < inAdj[top.ss].size(); j++){

				if(vis[inAdj[top.ss][j]] == 0){

					st.erase({pathWt[t + 1][inAdj[top.ss][j]], inAdj[top.ss][j]});
					if(pathWt[t + 1][inAdj[top.ss][j]] > pathWt[t + 1][top.ss] + (int)label[inAdj[top.ss][j]].size()){
						pathWt[t + 1][inAdj[top.ss][j]] = pathWt[t + 1][top.ss] + (int)label[inAdj[top.ss][j]].size();
						shortestPath[t + 1][inAdj[top.ss][j]] = shortestPath[t + 1][top.ss];
						shortestPath[t + 1][inAdj[top.ss][j]].push_back(inAdj[top.ss][j]);
					}

					st.insert({pathWt[t + 1][inAdj[top.ss][j]], inAdj[top.ss][j]});
				}
			}
		}
	}

	// O(NElog(N) + N*Nlog(N))
	for(int t = 0; t < Vc.size(); t++){
		int i = Vc[t];

		// O(N)
		set<pair<int, int> > st(seq.begin(), seq.end());

		// O(log(N))
		st.erase(make_pair(INT_MAX, i));
		st.insert(make_pair(0, i));
		pathWtR[t + 1][i] = (int)label[i].size();

		// O(N)
		// Visited array to check if the shortest distance of a node has been calculated or not
		vector<bool> vis(N + 1, 0);

		// O(Elog(N) + Nlog(N))
		while(!st.empty()){

			pair<int, int> top = *(st.begin());
			st.erase(st.begin());
			vis[top.ss] = 1;
			if(top.ff == INT_MAX)
				break;

			for(int j = 0; j < adj[top.ss].size(); j++){

				if(vis[adj[top.ss][j]] == 0){

					st.erase({pathWtR[t + 1][adj[top.ss][j]], adj[top.ss][j]});
					pathWtR[t + 1][adj[top.ss][j]] = min(pathWtR[t + 1][adj[top.ss][j]], pathWtR[t + 1][top.ss] + (int)label[adj[top.ss][j]																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																														].size());
					st.insert({pathWtR[t + 1][adj[top.ss][j]], adj[top.ss][j]});
				}
			}
		}
	}

	for(int t = 0; t < Vc.size(); t++){
		for(int i = 0; i < inAdj[Vc[t]].size(); i++){
			if(pathWt[t + 1][inAdj[Vc[t]][i]] == INT_MAX || pathWtR[t + 1][inAdj[Vc[t]][i]] == INT_MAX)
				continue;
			if(inAdj[Vc[t]][i] == Vc[t]){
				loopCost[t + 1] = 1;
				break;
			}
			loopCost[t + 1] = min(loopCost[t + 1], pathWtR[t + 1][inAdj[Vc[t]][i]] + pathWt[t + 1][inAdj[Vc[t]][i]] - (int)(label[Vc[t]].size() << 1) - (int)label[inAdj[Vc[t]][i]].size() + 1);
		}	
	}
}

// Returns the minimum feedback vertex set as computed using it's external implementation
vector<string> graphClass::mfvsExternalImplementation(int type){
	Graph g(N + 1);
	for(int i = 0; i <= N; i++)
		g.addVertex(i);

	for(int i = 0; i <= N; i++){
		for(int j = 0; j < adj[i].size(); j++)
			g.addEdge(i, adj[i][j]);
	}

	vector<int> vec =  g.minimumFeedbackVertexSet(true);
	vector<string> *retvec = new vector<string>;
	if(type == 0){
		for(int i=0; i <  vec.size(); ++i){
			(*retvec).push_back(to_string(vec[i]));
		}
	}else{
		for(int i=0; i <  vec.size(); ++i){
			(*retvec).push_back(revIx[vec[i]]);
		}
	}
	return (*retvec);
}

// Shortest Edge Distance from ith base of node u to jth base of node v
int graphClass::edgeDist(int u, int i, int v, int j){
	if(u == v && j < i){
		if(loopCost[v] == INT_MAX)
			return INT_MAX;
		return loopCost[v] + label[u].size() + j - i - 1;
	}

	if(pathWt[v][u] == INT_MAX)
		return INT_MAX;

	return j + pathWt[v][u] - i - label[Vc[v - 1]].size();
}

// Sets up the parameters for computing topological order
void graphClass::computeTopologicalOrder(){
	vis = vector<int>(N + 1, 0);
	
	for(int i = 1; i <= N; i++)
		ixTS[i] = vector<int>(label[i].size());

	topOrder.clear();
	
	int p = 1;

	topOrder.push_back({0, 0});

	for(int i = 0; i < Vc.size(); i++){
		topOrder.push_back({Vc[i], 0});
		p++;
		vis[Vc[i]] = 1;
	}

	for(int i = 0; i < Vc.size(); i++)
		for(int j = 1; j < label[Vc[i]].size(); j++)
			topOrder.push_back({Vc[i], j}), p++;

	for(int i = 1; i <= N; i++){

		if(vis[i] == false){
			topologicalSort(i);
			reverse(topOrder.begin() + p, topOrder.end());
			p = topOrder.size();
		}
	}

	// reverse(topOrder.begin() + p, topOrder.end());

	ixTS[0].push_back(0);

	for(int i = 1; i <= Va; i++)
		ixTS[topOrder[i].ff][topOrder[i].ss] = i;
}

// DFS like mfvsExternalImplementationon of topological sorting
void graphClass::topologicalSort(int v){
	
	vis[v] = 1;

	for(int i = 0; i < adj[v].size(); i++)
		if(vis[adj[v][i]] == false)
			topologicalSort(adj[v][i]);

	for(int i = label[v].size() - 1; i >= 0; i--)
		topOrder.push_back({v, i});
}

inline float graphClass::delta(int k){
	return gapOpen + k * gapExt;
}	

inline float graphClass::score(char w, char xj){
	return sc[chMap[w - 'A']][chMap[xj - 'A']];
}

// Computing Matrices M, R, Q
void graphClass::computeMatricesMRQ(string x, bool globalAlignment, pair<float, float> penalty){
	x = '_' + x;
	gapOpen = penalty.ff;
	gapExt = penalty.ss;

	int m = x.size() - 1;

	M = new float*[Va + 1];
	Q = new float*[Va + 1];
	R = new float*[Va + 1];

	fM = new bool*[Va + 1];
	fQ = new bool*[Va + 1];
	fR = new bool*[Va + 1];

	iM = new int*[Va + 1];
	iR = new int*[Va + 1];

	afvs = new pair<int, int>*[Vc.size() + 1];

	for(int i = 1; i <= Vc.size(); i++)
		afvs[i] = new pair<int, int>[m + 1];

	for(int i = 0; i <= Va; i++){
		M[i] = new float[m + 1];
		Q[i] = new float[m + 1];
		R[i] = new float[m + 1];
		fM[i] = new bool[m + 1]();
		fQ[i] = new bool[m + 1]();
		fR[i] = new bool[m + 1]();
		iM[i] = new int[m + 1];
		iR[i] = new int[m + 1];
	}

	for(int i = 0; i <= Va; i++){
		for(int j = 0; j <= m; j++){
			M[i][j] = INF;
			Q[i][j] = INF;
			R[i][j] = INF;
			iM[i][j] = -1;
			iR[i][j] = -1;
		}
	}

	// Initializing Matrix M
	for(int i = 0; i <= Va; i++)
		M[i][0] = 0.0F;

	for(int j = 1; j <= m; j++)
		M[0][j] = (globalAlignment == true ? -1.0 * delta(j) : 0.0F);

	float tmp = 0.0;

	for(int j = 1; j <= m; j++){
		for(int w = 1; w <= Vc.size(); w++){

			int v = topOrder[w].ff, ix = topOrder[w].ss;
			// Q (Seq Discard)
			Q[w][j] = M[w][j - 1] - delta(1);
			fQ[w][j] = 1;

			tmp = (j == 1 ? INF : Q[w][j - 1] - gapExt);
			if(Q[w][j] < tmp){	
				Q[w][j] = tmp;
				fQ[w][j] = 0;
			}

			M[w][j] = Q[w][j];

			// M Substitutions
			float scr = score(label[v][ix], x[j]);
			for(int i = 0; i < inAdj[v].size(); i++){

				int u = ixTS[inAdj[v][i]][(int)label[inAdj[v][i]].size() - 1];
				tmp = (M[u][j - 1] == INF ? INF : M[u][j - 1] + scr);
				if(M[w][j] < tmp){
					M[w][j] = tmp;
					iM[w][j] = u;
				}
			}

			// R Calculation
			for(int u = 0; u <= N; u++){
				for(int i = 0; i < label[u].size() - 1; i++){
					if(u == v && i + 1 == ix)
						continue;
					tmp = (M[ixTS[u][i]][j - 1] == INF || edgeDist(u, i + 1, w, ix) == INT_MAX ? INF : M[ixTS[u][i]][j - 1] + score(label[u][i + 1], x[j]) - delta(edgeDist(u, i + 1, w, ix)));
					if(R[w][j] < tmp){
						R[w][j] = tmp;
						fR[w][j] = 0;
						iR[w][j] = ixTS[u][i];
						afvs[w][j] = {u, i + 1};
					}
				}
				for(int k = 0; k < adj[u].size(); k++){
					if(adj[u][k] == v && ix == 0)
						continue;
					tmp = (M[ixTS[u][(int)label[u].size() - 1]][j - 1] == INF || edgeDist(adj[u][k], 0, w, ix) == INT_MAX ? INF : M[ixTS[u][(int)label[u].size() - 1]][j - 1] + score(label[adj[u][k]][0], x[j]) - delta(edgeDist(adj[u][k], 0, w, topOrder[w].ss)));
					if(R[w][j] < tmp){
						R[w][j] = tmp;
						fR[w][j] = 0;
						iR[w][j] = ixTS[u][(int)label[u].size() - 1];
						afvs[w][j] = {adj[u][k], 0};
					}
				}
			}

			if(M[w][j] < R[w][j]){
				M[w][j] = R[w][j];
				iM[w][j] = -1;
				fM[w][j] = 1;
			}

			M[w][j] = max(M[w][j], (globalAlignment == true ? INF : 0.0F));
		}

		for(int w = Vc.size() + 1; w <= Va; w++){

			int v = topOrder[w].ff, ix = topOrder[w].ss;

			// Q (Seq Discard)
			Q[w][j] = M[w][j - 1] - delta(1);
			fQ[w][j] = 1;

			tmp = (j == 1 ? INF : Q[w][j - 1] - gapExt);
			if(Q[w][j] < tmp){	
				Q[w][j] = tmp;
				fQ[w][j] = 0;
			}

			M[w][j] = Q[w][j];

			float scr = score(label[v][ix], x[j]);

			if(ix > 0){
				// M (Substitution)
				int u = ixTS[v][ix - 1];
				tmp = (M[u][j - 1] == INF ? INF : M[u][j - 1] + scr);
				if(M[w][j] < tmp){
					M[w][j] = tmp;	
					iM[w][j] = u;
				}

				// R Calculation
				tmp = (M[u][j] == INF ? INF : M[u][j] - delta(1));
				if(R[w][j] < tmp){
					R[w][j] = tmp;
					fR[w][j] = 0;
					iR[w][j] = u;
				}

				tmp = (R[u][j] == INF ? INF : R[u][j] - gapExt);
				if(R[w][j] < tmp){
					R[w][j] = tmp;
					fR[w][j] = 1;
					iR[w][j] = u;
				}
			}
			else{
				for(int i = 0; i < inAdj[v].size(); i++){

					// M Calculation
					int u = ixTS[inAdj[v][i]][(int)label[inAdj[v][i]].size() - 1];
					tmp = (M[u][j - 1] == INF ? INF : M[u][j - 1] + scr);
					if(M[w][j] < tmp){
						M[w][j] = tmp;
						iM[w][j] = u;
					}

					// R Calculation
					tmp = (M[u][j] == INF ? INF : M[u][j] - delta(1));
					if(R[w][j] < tmp){
						R[w][j] = tmp;
						fR[w][j] = 0;
						iR[w][j] = u;
					}

					tmp = (R[u][j] == INF ? INF : R[u][j] - gapExt);
					if(R[w][j] < tmp){
						R[w][j] = tmp;
						fR[w][j] = 1;
						iR[w][j] = u;
					}
				}
			}
		
			if(M[w][j] < R[w][j]){
				M[w][j] = R[w][j];
				fM[w][j] = 1;
				iM[w][j] = -1;
			}
			M[w][j] = max(M[w][j], (globalAlignment == true ? INF : 0.0F));
		}
	}

	/*
	cout << "\nM Matrix :: \n";
	
	cout << setw(10) << left << " ";
	for(int i = 0; i <= m; i++)
		cout << setw(10) << left << '(' + to_string(i) + ", " + x[i] + ')';
	cout << "\n\n";
	
	for(int i = 0; i <= Va; i++){
	
		cout << setw(10) << left << '(' + to_string(topOrder[i].ff) + ", " + to_string(topOrder[i].ss) + ", " + label[topOrder[i].ff][topOrder[i].ss] + ')';
		for(int j = 0; j <= m; j++){
	
			if(M[i][j] == INF)
				cout << setw(10) << left << "INF";
			else
				cout << setw(10) << left << M[i][j];
		}
		cout << "\n";
	}
	cout << "\n";

	cout << "Q Matrix :: \n";
	cout << setw(10) << left << " ";
	
	for(int i = 0; i <= m; i++)
		cout << setw(10) << left << '(' + to_string(i) + ", " + x[i] + ')';
	cout << "\n\n";
	
	for(int i = 0; i <= Va; i++){
	
		cout << setw(10) << left << '(' + to_string(topOrder[i].ff) + ", " + to_string(topOrder[i].ss) + ", " + label[topOrder[i].ff][topOrder[i].ss] + ')';
		for(int j = 0; j <= m; j++){
	
			if(Q[i][j] == INF)
				cout << setw(10) << left << "INF";
			else
				cout << setw(10) << left << Q[i][j];
		}
		cout << "\n";
	}
	cout << "\n";

	cout << "R Matrix :: \n";
	cout << setw(10) << left << " ";
	
	for(int i = 0; i <= m; i++)
		cout << setw(10) << left << '(' + to_string(i) + ", " + x[i] + ')';
	cout << "\n\n";
	
	for(int i = 0; i <= Va; i++){
	
		cout << setw(10) << left << '(' + to_string(topOrder[i].ff) + ", " + to_string(topOrder[i].ss) + ", " + label[topOrder[i].ff][topOrder[i].ss] + ')';
		for(int j = 0; j <= m; j++){
	
			if(R[i][j] == INF)
				cout << setw(10) << left << "INF";
			else
				cout << setw(10) << left << R[i][j];
		}
		cout << "\n";
	}
	cout << "\n";
	*/
}


void graphClass::alignSequence(string x, bool globalAlignment, pair<float, float> penalty, fstream &out, fstream &debug){
	printf("Starting MRQ Computation.\n");
	computeMatricesMRQ(x, globalAlignment, penalty);

	printf("MRQ Computation Done.\n");

	out << "\ta) Input Sequence: " << x << "\n";
	
	int m = x.size();
	float alignScore;

	x = '_' + x;

	exSeq = "", resSeq = "", resGraph = "";
	string exSeqOut = "";

	path.clear();
	tpPath.clear();

	if(globalAlignment == true){
		int ix = 0;
		alignScore = INF;
		for(int i = 0; i <= Va; i++){
			if(alignScore < M[i][m])
				alignScore = M[i][m], ix = i;
		}

		for(int i = ix, j = m; i > 0 && j > 0; ){

			if(iM[i][j] != -1){
				resSeq += x[j];
				resGraph += label[topOrder[i].ff][topOrder[i].ss];
				exSeq += label[topOrder[i].ff][topOrder[i].ss];
				exSeqOut += label[topOrder[i].ff][topOrder[i].ss];
				path.push_back(topOrder[i]);
				tpPath.push_back(topOrder[i]);
				i = iM[i][j], j -= 1;
			}
			else if(fM[i][j] == 0){
		
				while(fQ[i][j] == 0){
					resGraph += '-';
					resSeq += x[j];
					j -= 1;
				}
				resGraph += '-';
				resSeq += x[j];
				j -= 1;
			}
			else if(fM[i][j] == 1){

				while(fR[i][j] == 1){
					resGraph += label[topOrder[i].ff][topOrder[i].ss];
					exSeq += label[topOrder[i].ff][topOrder[i].ss];
					tpPath.push_back(topOrder[i]);
					resSeq += '-';
					int u = iR[i][j], v = j;
					if(i <= Vc.size()){

						if(topOrder[i].ff == afvs[i][j].ff){

							for(int k = topOrder[i].ss - 1; k > afvs[i][j].ss; k--){
								exSeq += label[topOrder[i].ff][k];
								resGraph += label[topOrder[i].ff][k];
								resSeq += '-';
								tpPath.push_back({topOrder[i].ff, k});
							}
						}
						else{
							for(int k = topOrder[i].ss - 1; k >= 0; k--){
								exSeq += label[topOrder[i].ff][k];
								resGraph += label[topOrder[i].ff][k];
								resSeq += '-';
								tpPath.push_back({topOrder[i].ff, k});
							}

							for(int k = 1; k < shortestPath[i][afvs[i][j].ff].size() - 1; k++){

								for(int p = label[shortestPath[i][afvs[i][j].ff][k]].size() - 1; p >= 0; p--){
									exSeq += label[shortestPath[i][afvs[i][j].ff][k]][p];
									resGraph += label[shortestPath[i][afvs[i][j].ff][k]][p];
									resSeq += '-';
									tpPath.push_back({shortestPath[i][afvs[i][j].ff][k], p});
								}
							}

							for(int k = label[afvs[i][j].ff].size() - 1; k > afvs[i][j].ss; k--){
								exSeq += label[afvs[i][j].ff][k];
								resGraph += label[afvs[i][j].ff][k];
								resSeq += '-';
								tpPath.push_back({afvs[i][j].ff, k});
							}
						}	

						exSeq += label[afvs[i][j].ff][afvs[i][j].ss];
						resGraph += label[afvs[i][j].ff][afvs[i][j].ss];
						resSeq += x[j];
						// path.push_back(afvs[i][j]);
						tpPath.push_back(afvs[i][j]);
						v -= 1;
					}
					i = u, j = v;
				}

				resGraph += label[topOrder[i].ff][topOrder[i].ss];
				exSeq += label[topOrder[i].ff][topOrder[i].ss];
				tpPath.push_back(topOrder[i]);
				resSeq += '-';
				int u = iR[i][j], v = j;
				if(i <= Vc.size()){

					if(topOrder[i].ff == afvs[i][j].ff){

						for(int k = topOrder[i].ss - 1; k > afvs[i][j].ss; k--){
							exSeq += label[topOrder[i].ff][k];
							resGraph += label[topOrder[i].ff][k];
							resSeq += '-';
							tpPath.push_back({topOrder[i].ff, k});
						}
					}
					else{
						for(int k = topOrder[i].ss - 1; k >= 0; k--){
							exSeq += label[topOrder[i].ff][k];
							resGraph += label[topOrder[i].ff][k];
							resSeq += '-';
							tpPath.push_back({topOrder[i].ff, k});
						}

						for(int k = 1; k < shortestPath[i][afvs[i][j].ff].size() - 1; k++){

							for(int p = label[shortestPath[i][afvs[i][j].ff][k]].size() - 1; p >= 0; p--){
								exSeq += label[shortestPath[i][afvs[i][j].ff][k]][p];
								resGraph += label[shortestPath[i][afvs[i][j].ff][k]][p];
								resSeq += '-';
								tpPath.push_back({shortestPath[i][afvs[i][j].ff][k], p});
							}
						}

						for(int k = label[afvs[i][j].ff].size() - 1; k > afvs[i][j].ss; k--){
							exSeq += label[afvs[i][j].ff][k];
							resGraph += label[afvs[i][j].ff][k];
							resSeq += '-';
							tpPath.push_back({afvs[i][j].ff, k});
						}
					}	

					exSeq += label[afvs[i][j].ff][afvs[i][j].ss];
					resGraph += label[afvs[i][j].ff][afvs[i][j].ss];
					resSeq += x[j];
					// path.push_back(afvs[i][j]);
					tpPath.push_back(afvs[i][j]);
					v -= 1;
				}
				i = u, j = v;
			}
		}
	}	
	else{
		
		int ix, iy;
		alignScore = INF;
		for(int i = 0; i <= Va; i++){
			for(int j = 0; j <= m; j++)
				if(alignScore < M[i][j])
					alignScore = M[i][j], ix = i, iy = j;
		}

		for(int i = ix, j = iy; M[i][j] != 0 && i > 0 && j > 0; ){
			
			if(iM[i][j] != -1){
				resSeq += x[j];
				resGraph += label[topOrder[i].ff][topOrder[i].ss];
				exSeq += label[topOrder[i].ff][topOrder[i].ss];
				exSeqOut += label[topOrder[i].ff][topOrder[i].ss];
				path.push_back(topOrder[i]);
				tpPath.push_back(topOrder[i]);
				i = iM[i][j], j -= 1;
			}
			else if(fM[i][j] == 0){
		
				while(fQ[i][j] == 0){
					resGraph += '-';
					resSeq += x[j];
					j -= 1;
				}
				resGraph += '-';
				resSeq += x[j];
				j -= 1;
			}
			else if(fM[i][j] == 1){

				while(fR[i][j] == 1){
					resGraph += label[topOrder[i].ff][topOrder[i].ss];
					exSeq += label[topOrder[i].ff][topOrder[i].ss];
					tpPath.push_back(topOrder[i]);
					resSeq += '-';
					int u = iR[i][j], v = j;
					if(i <= Vc.size()){

						if(topOrder[i].ff == afvs[i][j].ff){

							for(int k = topOrder[i].ss - 1; k > afvs[i][j].ss; k--){
								exSeq += label[topOrder[i].ff][k];
								resGraph += label[topOrder[i].ff][k];
								resSeq += '-';
								tpPath.push_back({topOrder[i].ff, k});
							}
						}
						else{
							for(int k = topOrder[i].ss - 1; k >= 0; k--){
								exSeq += label[topOrder[i].ff][k];
								resGraph += label[topOrder[i].ff][k];
								resSeq += '-';
								tpPath.push_back({topOrder[i].ff, k});
							}

							for(int k = 1; k < shortestPath[i][afvs[i][j].ff].size() - 1; k++){

								for(int p = label[shortestPath[i][afvs[i][j].ff][k]].size() - 1; p >= 0; p--){
									exSeq += label[shortestPath[i][afvs[i][j].ff][k]][p];
									resGraph += label[shortestPath[i][afvs[i][j].ff][k]][p];
									resSeq += '-';
									tpPath.push_back({shortestPath[i][afvs[i][j].ff][k], p});
								}
							}

							for(int k = label[afvs[i][j].ff].size() - 1; k > afvs[i][j].ss; k--){
								exSeq += label[afvs[i][j].ff][k];
								resGraph += label[afvs[i][j].ff][k];
								resSeq += '-';
								tpPath.push_back({afvs[i][j].ff, k});
							}
						}	

						exSeq += label[afvs[i][j].ff][afvs[i][j].ss];
						resGraph += label[afvs[i][j].ff][afvs[i][j].ss];
						resSeq += x[j];
						// path.push_back(afvs[i][j]);
						tpPath.push_back(afvs[i][j]);
						v -= 1;
					}
					i = u, j = v;
				}

				resGraph += label[topOrder[i].ff][topOrder[i].ss];
				exSeq += label[topOrder[i].ff][topOrder[i].ss];
				tpPath.push_back(topOrder[i]);
				resSeq += '-';
				int u = iR[i][j], v = j;
				if(i <= Vc.size()){

					if(topOrder[i].ff == afvs[i][j].ff){

						for(int k = topOrder[i].ss - 1; k > afvs[i][j].ss; k--){
							exSeq += label[topOrder[i].ff][k];
							resGraph += label[topOrder[i].ff][k];
							resSeq += '-';
							tpPath.push_back({topOrder[i].ff, k});
						}
					}
					else{
						for(int k = topOrder[i].ss - 1; k >= 0; k--){
							exSeq += label[topOrder[i].ff][k];
							resGraph += label[topOrder[i].ff][k];
							resSeq += '-';
							tpPath.push_back({topOrder[i].ff, k});
						}

						for(int k = 1; k < shortestPath[i][afvs[i][j].ff].size() - 1; k++){

							for(int p = label[shortestPath[i][afvs[i][j].ff][k]].size() - 1; p >= 0; p--){
								exSeq += label[shortestPath[i][afvs[i][j].ff][k]][p];
								resGraph += label[shortestPath[i][afvs[i][j].ff][k]][p];
								resSeq += '-';
								tpPath.push_back({shortestPath[i][afvs[i][j].ff][k], p});
							}
						}

						for(int k = label[afvs[i][j].ff].size() - 1; k > afvs[i][j].ss; k--){
							exSeq += label[afvs[i][j].ff][k];
							resGraph += label[afvs[i][j].ff][k];
							resSeq += '-';
							tpPath.push_back({afvs[i][j].ff, k});
						}
					}	

					exSeq += label[afvs[i][j].ff][afvs[i][j].ss];
					resGraph += label[afvs[i][j].ff][afvs[i][j].ss];
					resSeq += x[j];
					// path.push_back(afvs[i][j]);
					tpPath.push_back(afvs[i][j]);
					v -= 1;
				}
				i = u, j = v;
			}
		}
	}	

	reverse(resGraph.begin(), resGraph.end());
	reverse(exSeq.begin(), exSeq.end());
	reverse(exSeqOut.begin(), exSeqOut.end());
	reverse(resSeq.begin(), resSeq.end());
	reverse(path.begin(), path.end());
	reverse(tpPath.begin(), tpPath.end());

	out << "\tb) Alignment Score: " << alignScore << "\n";
	out << "\tc) Exact Sequence from the graph: " << exSeq << "\n";
	out << "\td) Path Encoding the Sequence: ";

	if(!path.empty()){
		int i = 0, st = path[0].ss, ed = path[0].ss, j = i;
		string seq;
		seq += exSeqOut[i];
		for(j = i + 1; j < path.size() && path[j].ff == path[i].ff && path[j].ss == path[j - 1].ss + 1; j++)
			seq += exSeqOut[j], ed = path[j].ss;
		out << '<' << revIx[path[i].ff] << (path[i].ff <= V ? "" : "\'") << ':' << to_string(st) << ':' << seq << ':' << to_string(ed) << "> ";
		i = j;
		for(; i < path.size();){
			seq = "";
			seq += exSeqOut[i];
			st = path[i].ss, ed = path[i].ss;
			for(j = i + 1; j < path.size() && path[j].ff == path[i].ff && path[j].ss == path[j - 1].ss + 1; j++)
				seq += exSeqOut[j], ed = path[j].ss;
			out << "-> <" << revIx[path[i].ff] << (path[i].ff <= V ? "" : "\'") << ':' << to_string(st) << ':' << seq << ':' << to_string(ed) << "> ";
			i = j;
		} 
	}

	out << "\n\te) Alignment: \n";
	out << "\t\tGraph Sequence:   ";
	out << setw(5) << right << resGraph;
	out << "\n\t\tAligned Sequence: ";
	out << setw(5) << right << resSeq;
	out << "\n\n";

	if(!debug)	
		printMRQ(x, debug);

	for(int i = 0; i <= Va; i++){
		delete []M[i];
		delete []Q[i];
		delete []R[i];
		delete []fM[i];
		delete []fQ[i];
		delete []fR[i];
		delete []iM[i];
		delete []iR[i];
	}

	for(int i = 1; i <= Vc.size(); i++)
		delete []afvs[i];

	printf("Match - %d\n", (resSeq == resGraph));
}

void graphClass::visualizeGFA(fstream &dotFile){
	dotFile << "digraph graphname {\n";
	dotFile << "\trankdir=LR;\n";
	for(int i = 1; i <= N; i++){
		
		dotFile << "\t" << i <<" [label=<";

		dotFile << "<font color=";
		if(i <= V)
			dotFile << "\"black\"";
		else
			dotFile << "\"violet\"";	
		dotFile << ">" << revIx[i] << ":</font><font>" << label[i] << "</font>";
		dotFile << ">];\n";
	}

	for(int i = 1; i <= N; i++){
		for(int j = 0; j < adj[i].size(); j++){
			dotFile << "\t" << i << " -> " << adj[i][j] << " [label=";

			if(i <= V && adj[i][j] <= V) 		// u -> v
				dotFile << "\"(" << revIx[i] << "+" << revIx[adj[i][j]] << "+)\"";
			else if(i <= V && adj[i][j] > V)	// u -> vrc
				dotFile << "\"(" << revIx[i] << "+" << revIx[adj[i][j]] << "-)\"";
			else if(i > V && adj[i][j] <= V)	// urc -> v
				dotFile << "\"(" << revIx[i] << "-" << revIx[adj[i][j]] << "+)\"";
			else if(i > V && adj[i][j] > V)		// urc -> vrc
				dotFile << "\"(" << revIx[i] << "-" << revIx[adj[i][j]] << "-)\"";

			dotFile << "];\n";
		}
	}

	dotFile << "}\n";
}

void graphClass::mfvsVisualize(fstream &dotFile, int depth){
	dotFile << "digraph graphname {\n";
	// dotFile << "\trankdir=LR;\n";

	vis = vector<int>(N + 1, 0);
	vector<int> vertices;
	set<pair<int, int> > edges;

	for(int i = 0; i < Vc.size(); i++){
		queue<int> q, lvl;
		if(vis[Vc[i]] == 0)
			q.push(Vc[i]), lvl.push(0);

		while(!q.empty()){
			int node = q.front(), level = lvl.front();
			q.pop(), lvl.pop();

			vertices.push_back(node);

			if(level == depth)
				continue;

			for(int j = 0; j < adj[node].size(); j++){
				if(vis[adj[node][j]] == 0){
					vis[adj[node][j]] = 1;
					edges.insert({node, adj[node][j]});
					q.push(adj[node][j]);
					lvl.push(level + 1);
				}
				else
					edges.insert({node, adj[node][j]});
			}
		}

		q.push(Vc[i]), lvl.push(0);

		while(!q.empty()){
			int node = q.front(), level = lvl.front();
			q.pop(), lvl.pop();

			vertices.push_back(node);

			if(level == depth)
				continue;

			for(int j = 0; j < inAdj[node].size(); j++){
				if(vis[inAdj[node][j]] == 0){
					vis[inAdj[node][j]] = 1;
					edges.insert({inAdj[node][j], node});
					q.push(inAdj[node][j]);
					lvl.push(level + 1);
				}
				else
					edges.insert({inAdj[node][j], node});
			}
		}		
	}

	for(int i = 0; i < vertices.size(); i++){
		
		dotFile << "\t" << vertices[i] <<" [label=<";

		dotFile << "<font color=";
		if(vertices[i] <= V)
			dotFile << "\"black\"";
		else
			dotFile << "\"violet\"";	
		dotFile << ">" << revIx[vertices[i]] << ":</font><font>" << label[vertices[i]] << "</font>";
		dotFile << ">];\n";
	}

	for(set<pair<int, int> >::iterator it = edges.begin(); it != edges.end(); it++)
		dotFile << "\t" << it -> ff << " -> " << it -> ss << "\n";

	dotFile << "}\n";
}


void graphClass::writeDotFile(fstream &dotFile){
	dotFile << "digraph graphname {\n";
	int cnt = 0, p = 0;

	if(!tpPath.empty()){
		int i = 0, st = tpPath[0].ss, ed = tpPath[0].ss, j = i;
		string seq;
		seq += exSeq[i];
		for(j = i + 1; j < tpPath.size() && tpPath[j].ff == tpPath[i].ff && tpPath[j].ss == tpPath[j - 1].ss + 1; j++)
			seq += exSeq[j], ed = tpPath[j].ss;
		
		if(tpPath[i].ff <= V)	
			dotFile << "\t" << cnt + 1 << " [label=<<font>" << revIx[tpPath[i].ff] << ":</font>";
		else
			dotFile << "\t" << cnt + 1 << " [label=<<font color=\"violet\">" << revIx[tpPath[i].ff] << ":</font>";
		
		if(tpPath[i].ff <= (V << 1)){
			for(int k = 0; k < st; k++)
				dotFile << "<font color=\"grey\">" << label[tpPath[i].ff][k] << "</font>";

			for(int k = 0; k < seq.size(); k++){

				if(resSeq[p] == '-')
					dotFile << "<font color=\"red\">" << seq[k] << "</font>", p++;
				else if(resGraph[p] == resSeq[p])
					dotFile << "<font color=\"black\">" << seq[k] << "</font>", p++;
				else 
					dotFile << "<font color=\"blue\">" << seq[k] << "</font>", p++;

				if(p < resGraph.size() && resGraph[p] == '-'){
					// dotFile << "<font>-</font>";
					while(p < resGraph.size() && resGraph[p] == '-')
						dotFile << "<font>-</font>", p += 1;
				}
			}

			for(int k = ed + 1; k < label[tpPath[i].ff].size(); k++)
				dotFile << "<font color=\"grey\">" << label[tpPath[i].ff][k] << "</font>";

		}
		else{

			string col;

			for(int k = 0; k < st; k++)
				col += 'G';

			for(int k = 0; k < seq.size(); k++){

				if(resSeq[p] == '-')
					col += 'R', p++;
				else if(resGraph[p] == resSeq[p])
					col += 'B', p++;
				else 
					col += 'V', p++;

				if(p < resGraph.size() && resGraph[p] == '-'){
					// col += '-';
					while(p < resGraph.size() && resGraph[p] == '-')
						col += '-', p += 1;
				}
			}

			for(int k = ed + 1; k < label[tpPath[i].ff].size(); k++)
				col += 'G';

			reverse(col.begin(), col.end());

			for(int k = 0, l = 0; k < col.size(); k++){
				if(col[k] == '-')
					dotFile << "<font>-</font>";
				else if(col[k] == 'G')
					dotFile << "<font color=\"grey\">" << label[tpPath[i].ff - V][l++] << "</font>";
				else if(col[k] == 'R')
					dotFile << "<font color=\"red\">" << label[tpPath[i].ff - V][l++] << "</font>";
				else if(col[k] == 'B')
					dotFile << "<font color=\"black\">" << label[tpPath[i].ff - V][l++] << "</font>";
				else
					dotFile << "<font color=\"blue\">" << label[tpPath[i].ff - V][l++] << "</font>";
			}
			
		}


		dotFile << ">];\n";
		cnt += 1;

		i = j;
		for(; i < tpPath.size();){
			seq = "";
			seq += exSeq[i];
			st = tpPath[i].ss, ed = tpPath[i].ss;
			for(j = i + 1; j < tpPath.size() && tpPath[j].ff == tpPath[i].ff && tpPath[j].ss == tpPath[j - 1].ss + 1; j++)
				seq += exSeq[j], ed = tpPath[j].ss;

			if(tpPath[i].ff <= V)	
				dotFile << "\t" << cnt + 1 << " [label=<<font>" << revIx[tpPath[i].ff] << ":</font>";
			else
				dotFile << "\t" << cnt + 1 << " [label=<<font color=\"violet\">" << revIx[tpPath[i].ff] << ":</font>";
			
			if(tpPath[i].ff <= (V << 1)){
				for(int k = 0; k < st; k++)
					dotFile << "<font color=\"grey\">" << label[tpPath[i].ff][k] << "</font>";

				for(int k = 0; k < seq.size(); k++){

					if(resSeq[p] == '-')
						dotFile << "<font color=\"red\">" << seq[k] << "</font>", p++;
					else if(resGraph[p] == resSeq[p])
						dotFile << "<font color=\"black\">" << seq[k] << "</font>", p++;
					else 
						dotFile << "<font color=\"blue\">" << seq[k] << "</font>", p++;

					if(p < resGraph.size() && resGraph[p] == '-'){
						// dotFile << "<font>-</font>";
						while(p < resGraph.size() && resGraph[p] == '-')
							dotFile << "<font>-</font>", p += 1;
					}
				}

				for(int k = ed + 1; k < label[tpPath[i].ff].size(); k++)
					dotFile << "<font color=\"grey\">" << label[tpPath[i].ff][k] << "</font>";

			}
			else{

				string col;

				for(int k = 0; k < st; k++)
					col += 'G';

				for(int k = 0; k < seq.size(); k++){

					if(resSeq[p] == '-')
						col += 'R', p++;
					else if(resGraph[p] == resSeq[p])
						col += 'B', p++;
					else 
						col += 'V', p++;

					if(p < resGraph.size() && resGraph[p] == '-'){
						// col += '-';
						while(p < resGraph.size() && resGraph[p] == '-')
							col += '-', p += 1;
					}
				}

				for(int k = ed + 1; k < label[tpPath[i].ff].size(); k++)
					col += 'G';

				reverse(col.begin(), col.end());

				for(int k = 0, l = 0; k < col.size(); k++){
					if(col[k] == '-')
						dotFile << "<font>-</font>";
					else if(col[k] == 'G')
						dotFile << "<font color=\"grey\">" << label[tpPath[i].ff - V][l++] << "</font>";
					else if(col[k] == 'R')
						dotFile << "<font color=\"red\">" << label[tpPath[i].ff - V][l++] << "</font>";
					else if(col[k] == 'B')
						dotFile << "<font color=\"black\">" << label[tpPath[i].ff - V][l++] << "</font>";
					else
						dotFile << "<font color=\"blue\">" << label[tpPath[i].ff - V][l++] << "</font>";
				}
				
			}


			dotFile << ">];\n";
			cnt += 1;
			i = j;
		} 
	}

	for(int i = 1; i < cnt; i++)
		dotFile << "\t" << i << " -> " << i + 1 << "\n";

	dotFile << "}";
}	


void graphClass::printVirtualGraph(fstream &debug){
	debug << "Number of Nodes: " << Va << "\n";
	debug << "Number of Edges: " << E + Va - N << "\n";
	debug << "\nAdjacency List: \n\n";

	for(int i = 1; i <= N; i++){
		for(int j = 0; j < label[i].size() - 1; j++)
			debug << setw(50) << left << "(" + to_string(i) + ", " + to_string(j) + ") - (" + to_string(i) + ", " + to_string(j + 1) + ")" << "\n"; 
		
		for(int j = 0; j < adj[i].size(); j++)
			debug << setw(50) << left << "(" + to_string(i) + ", " + to_string((int)label[i].size() - 1) + ") - (" + to_string(adj[i][j]) + ", " + to_string(0) + ")" << "\n"; 	
	}		
}

void graphClass::printMRQ(string x, fstream &debug){
	x = '_' + x;

	int m = x.size() - 1;

	debug << "\nM Matrix :: \n";
	
	debug << setw(10) << left << " ";
	for(int i = 0; i <= m; i++)
		debug << setw(10) << left << '(' + to_string(i) + ", " + x[i] + ')';
	debug << "\n\n";
	
	for(int i = 0; i <= Va; i++){
	
		debug << setw(10) << left << '(' + to_string(topOrder[i].ff) + ", " + to_string(topOrder[i].ss) + ", " + label[topOrder[i].ff][topOrder[i].ss] + ')';
		for(int j = 0; j <= m; j++){
	
			if(M[i][j] == INF)
				debug << setw(10) << left << "INF";
			else
				debug << setw(10) << left << M[i][j];
		}
		debug << "\n";
	}
	debug << "\n";

	debug << "Q Matrix :: \n";
	debug << setw(10) << left << " ";
	
	for(int i = 0; i <= m; i++)
		debug << setw(10) << left << '(' + to_string(i) + ", " + x[i] + ')';
	debug << "\n\n";
	
	for(int i = 0; i <= Va; i++){
	
		debug << setw(10) << left << '(' + to_string(topOrder[i].ff) + ", " + to_string(topOrder[i].ss) + ", " + label[topOrder[i].ff][topOrder[i].ss] + ')';
		for(int j = 0; j <= m; j++){
	
			if(Q[i][j] == INF)
				debug << setw(10) << left << "INF";
			else
				debug << setw(10) << left << Q[i][j];
		}
		debug << "\n";
	}
	debug << "\n";

	debug << "R Matrix :: \n";
	debug << setw(10) << left << " ";
	
	for(int i = 0; i <= m; i++)
		debug << setw(10) << left << '(' + to_string(i) + ", " + x[i] + ')';
	debug << "\n\n";
	
	for(int i = 0; i <= Va; i++){
	
		debug << setw(10) << left << '(' + to_string(topOrder[i].ff) + ", " + to_string(topOrder[i].ss) + ", " + label[topOrder[i].ff][topOrder[i].ss] + ')';
		for(int j = 0; j <= m; j++){
	
			if(R[i][j] == INF)
				debug << setw(10) << left << "INF";
			else
				debug << setw(10) << left << R[i][j];
		}
		debug << "\n";
	}
	debug << "\n";
}

