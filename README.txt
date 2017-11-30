
- System Requirements: Linux, g++ compiler, graphviz tool support (if you want to visualize the alignment output)

- valign tool required the MVFS set to be input as a separate file consisting of a new line separated list of feedback vertex set.

[ For instance, the MFVS implementation from the MFVS project (https://github.com/ablondin/mfvs/) could be used to create the MFVS vertex list ]

- edna.mat is from ftp://ftp.ncbi.nih.gov/blast/matrices/NUC.4.4
        



Usage: ./valign [OPTION]... [FILE]...
Arguments.

	-g filePath             tag g or G represents that the given path is either for a adjacency file or a
	-G filePath             gfa file or a dot file.  

	-x filePath             tag x or X represents that the path is for a test sequence file.
	-X filePath             test sequence file can have newline separated sequences for alignment.

	-go realNumber          tag go or GO represents that the given real number is the gap open cost.
	-GO realNumber          By default the gap open cost is 10.0

	-ge realNumber          tag ge or GE represents that the given real number is the gap extension cost.
	-GE realNumber          By default the gap extension cost is 0.5

	-global                 tag global or GLOBAL represents that the alignment of the test sequences must be global.
	-GLOBAL                 By default it is global alignment

	-local                  tag local or LOCAL represents that the alignment of the test sequences must be local.
	-LOCAL                  By default it is global alignment

	-o filePath             tag o or O represents that the path is for a output file.
	-O filePath             By default the output will be stored in a "out.txt" file in the current directory.

	-d filePath             tag d or D represents that the path is for a debug file.
	-D filePath             By default the output will be stored in a "debug.txt" file in the current directory.

	-v filePath             tag v or V represents that the path is for a Feedback Vertex Set file.
	-V filePath             This is not mandatory argument. If this is included then the function to calculate FVS will be bypassed.

	-dot directoryPath      tag dot or DOT represents that the path is for a directory to store the generated dot files.
	-DOT directoryPath      It will create a folder with a suffix "DotVisuals" in the directory given in the path and 
                         	stores the dot files and shell script in the created directory.

Reference: 
Sequence Alignment On Directed Graphs
Kavya Vaddadi, Naveen Sivadasan, Kshitij Tayal, Rajgopal Srinivasan
http://www.biorxiv.org/content/early/2017/04/06/124941

