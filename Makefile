CXX=g++
#CPPFLAGS=-c -std=c++11 -Ofast -march=native -flto -fwhole-program
CPPFLAGS=-c -std=c++11 -O3 -march=native -flto -fwhole-program
RM=rm -rf
all: valign
valign: sequenceAlignment.o graphClassMemberFunctions.o
	$(CXX) -o valign sequenceAlignment.o graphClassMemberFunctions.o
sequenceAlignment.o: sequenceAlignment.cpp
	$(CXX) $(CPPFLAGS) sequenceAlignment.cpp
graphClassMemberFunctions.o: graphClassMemberFunctions.cpp
	$(CXX) $(CPPFLAGS) graphClassMemberFunctions.cpp
clean:
	$(RM) *.o valign 
