// Author: Shreyansh Chhajer 
// Tata Consultancy Services Ltd

#include <vector>
#include <map>
#include <utility>
#include <time.h>

#include "graphClass.h"

using namespace std;

#define ff first
#define ss second
timespec st,en;
double runtime;

void printHelp(){
	printf("Usage: ./valign [OPTION]... [FILE]...");
	printf("\nArguments.\n");
	cout << endl << setw(25) << left << "\t-g filePath" << setw(25) << left << "tag g or G represents that the given path is either for a adjacency file or a" << endl;
	cout << setw(25) << left << "\t-G filePath" << setw(25) << left << "gfa file or a dot file." << endl;
	
	cout << endl << setw(25) << left << "\t-x filePath" << setw(25) << left << "tag x or X represents that the path is for a test sequence file." << endl;
	cout << setw(25) << left << "\t-X filePath" << setw(25) << left << "test sequence file can have newline separated sequences for alignment." << endl;

	cout << endl << setw(25) << left << "\t-go realNumber" << setw(25) << left << "tag go or GO represents that the given real number is the gap open cost." << endl;
	cout << setw(25) << left << "\t-GO realNumber" << setw(25) << left << "By default the gap open cost is 10.0" << endl;

	cout << endl << setw(25) << left << "\t-ge realNumber" << setw(25) << left << "tag ge or GE represents that the given real number is the gap extension cost." << endl;
	cout << setw(25) << left << "\t-GE realNumber" << setw(25) << left << "By default the gap extension cost is 0.5" << endl;

	cout << endl << setw(25) << left << "\t-global" << setw(25) << left << "tag global or GLOBAL represents that the alignment of the test sequences must be global." << endl;
	cout << setw(25) << left << "\t-GLOBAL" << setw(25) << left << "By default it is global alignment" << endl;

	cout << endl << setw(25) << left << "\t-local" << setw(25) << left << "tag local or LOCAL represents that the alignment of the test sequences must be local." << endl;
	cout << setw(25) << left << "\t-LOCAL" << setw(25) << left << "By default it is global alignment" << endl;

	cout << endl << setw(25) << left << "\t-o filePath" << setw(25) << left << "tag o or O represents that the path is for a output file." << endl;
	cout << setw(25) << left << "\t-O filePath" << setw(25) << left << "By default the output will be stored in a \"out.txt\" file in the current directory."<< endl;

	cout << endl << setw(25) << left << "\t-d filePath" << setw(25) << left << "tag d or D represents that the path is for a debug file." << endl;
	cout << setw(25) << left << "\t-D filePath" << setw(25) << left << "By default the output will be stored in a \"debug.txt\" file in the current directory."<< endl;

	cout << endl << setw(25) << left << "\t-v filePath" << setw(25) << left << "tag v or V represents that the path is for a Feedback Vertex Set file." << endl;
	cout << setw(25) << left << "\t-V filePath" << setw(25) << left << "This is not mandatory argument. If this is included then the function to calculate FVS will be bypassed."<< endl;

	cout << endl << setw(25) << left << "\t-dot directoryPath" << setw(25) << left << "tag dot or DOT represents that the path is for a directory to store the generated dot files." << endl;
	cout << setw(25) << left << "\t-DOT directoryPath" << setw(25) << left << "It will create a folder with a suffix \"DotVisuals\" in the directory given in the path and \n"; 
	cout << setw(25) << left << " " << setw(25) << left << "\tstores the dot files and shell script in the created directory."<< endl;	

}

int main(int argc, char** argv){

	int _; // = system("clear");
	// clock_gettime(CLOCK_MONOTONIC,&st);
	
	if(argc == 1){
		printf("Please provide required input arguments.\n");
		printf("For usage: ./valign --help\n");
		return 0;
	}

	if(string(argv[1]) == "--help"){
		printHelp();
		return 0;
	}

	string seqFile, adjList, outFile = "out.txt", debugFile = "";
	pair<float, float> penaltyParam = {10, 0.5};
	bool globalAlign = true;
	string mfvsFile = "", scoreFile = "";

	string::size_type sz;

	if(argc >= 3) cout << "Input graph : " << argv[2] << endl;
	string dotName = "";

	bool genadjlist = false;

	for(int i = 1; i < argc - 1; i += 2){
		string tp = argv[i];
		if(tp == "-x" || tp == "-X")
			seqFile = argv[i + 1];
		else if(tp == "-g" || tp == "-G")
			adjList = argv[i + 1];
		else if(tp == "-go" || tp == "-GO")
			penaltyParam.ff = stod(argv[i + 1], &sz);
		else if(tp == "-ge" || tp == "-GE")
			penaltyParam.ss = stod(argv[i + 1], &sz);
		else if(tp == "-global" || tp == "-GLOBAL")
			globalAlign = true, i -= 1;
		else if(tp == "-local" || tp == "-LOCAL")
			globalAlign = false, i -= 1;
		else if(tp == "-o" || tp == "-O")
			outFile = argv[i + 1];
		else if(tp == "-d" || tp == "-D")
			debugFile = argv[i + 1];
		else if(tp == "-v" || tp == "-V")
			mfvsFile = argv[i + 1];
		else if(tp == "-dot" || tp == "-DOT")
			dotName = argv[i + 1];
		else if(tp == "-s" || tp == "-S")
			scoreFile = argv[i + 1];
		else if(tp == "-genadjlist")
			genadjlist = true;
	}
	if(argc > 1){
		string tp = argv[argc - 1];
		if(tp == "-global" || tp == "-GLOBAL")
			globalAlign = true;
		else if(tp == "-local" || tp == "-LOCAL")
			globalAlign = false;
	}

	
	if(adjList == ""){
		cout << "Please provide graph file " << endl;
		exit(1);
	}

	if(genadjlist){
		string listfile = adjList;
		listfile.append(".list");
		exit(0);
	}

	if(seqFile == ""){
		printf("Please provide the query file.\n");
		return 0;
	}

	if(dotName.length() == 0) dotName = "dotVisuals";
	dotName.append("/");
	//dotName += name + "DotVisuals/";
	string systemCall = "mkdir " + dotName; 

	_ = system(systemCall.c_str());

	int type = 1; //adj format
	string suff = adjList.substr(adjList.length() - 3, 3);
	if(suff == "gfa") {
		type = 0;
	}else if(suff == "dot"){
		type = 2;
	}

	//cout << "File format : " << type << " " << suff << endl;

	graphClass G(type, adjList, mfvsFile, scoreFile);


	// fstream mfvsDot("mfvsDot.dot", fstream::out);
	// G.mfvsVisualize(mfvsDot, 5);

	// mfvsDot.close();

	// return 0;

	if(type == 0){
		adjList.pop_back();
		adjList.pop_back();
		adjList.pop_back();
		adjList += "visual.dot";

		fstream dotFile(adjList, fstream::out);

		G.visualizeGFA(dotFile);
	}

	fstream out, seq, debug;

	if(!debugFile.empty())
		debug.open(debugFile, fstream::out);

	out.open(outFile, fstream::out);
	seq.open(seqFile, fstream::in);

	string x;
	out << "E: " << G.E << "\nEa: " << G.Ea << "\nN: " << G.N << "\nVa: " << G.Va << "\n\n";

	int testCase = 1;

	//cout << "Align Time in Milliseconds =" << 
	long tot_time = 0;
	while(getline(seq, x)){
		out << "Qry Seq: " << to_string(testCase) << "\n";
		long start_time = chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now().time_since_epoch()).count();
		G.alignSequence(x, globalAlign, penaltyParam, out, debug);	
		long end_time = chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now().time_since_epoch()).count();
		tot_time += end_time - start_time;

		string dotFileName = dotName  + "qry." + to_string(testCase) + ".dot";
		
		fstream dotFile(dotFileName, fstream::out);
		G.writeDotFile(dotFile);	
		dotFile.close();

		printf("Qry seq - %d Done\n", testCase);
		
		testCase += 1;
	}
	cout << "Avg alignment time  = " << ((double)(tot_time)/((testCase -1) * 1000.0)) << " Sec." << endl; 

	out.close();
	seq.close();

	fstream dotFileRun(dotName + "dotRun" +  ".sh", fstream::out);

	for(int i = 1; i < testCase; i++)
		dotFileRun << "dot " << "qry." + to_string(i) + ".dot" << " -Tpdf -o " <<  "qry." + to_string(i) + ".pdf\n";

	dotFileRun.close();

	cout << "\nAlignment output is in " <<  outFile << endl;
	cout << "Dot files for visualization are present in the folder " <<  dotName << endl;
	


	// clock_gettime(CLOCK_MONOTONIC,&en);
	// runtime = en.tv_sec - st.tv_sec+(en.tv_nsec-st.tv_nsec)/(1e9);

	// fstream timeFile("execTime.txt", fstream::out | fstream::app);
	// timeFile << argv[2] << " " << runtime << "\n";

	// timeFile.close();

	return 0;
}
