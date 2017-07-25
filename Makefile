CXX=g++
CPPFLAGS=-c -std=c++11 -O9
RM=rm -rf
all: valign
valign: sequenceAlignment.o graphClassMemberFunctions.o utils.o graph.o
	$(CXX) -o valign sequenceAlignment.o graphClassMemberFunctions.o utils.o graph.o
sequenceAlignment.o: sequenceAlignment.cpp
	$(CXX) $(CPPFLAGS) sequenceAlignment.cpp
graphClassMemberFunctions.o: graphClassMemberFunctions.cpp
	$(CXX) $(CPPFLAGS) graphClassMemberFunctions.cpp
utils.o: MFVSImp/utils.cpp
	$(CXX) $(CPPFLAGS) MFVSImp/utils.cpp
graph.o: MFVSImp/graph.cpp
	$(CXX) $(CPPFLAGS) MFVSImp/graph.cpp
clean:
	$(RM) *.o valign 
