CXX=$(shell root-config --cxx)
CXXFLAGS=-Wall -O2 -g

ROOTCXXFLAGS=$(shell root-config --cflags) # already has root-config --incdir appended
ROOTLIBS=$(shell root-config --libs)

all: createRootFiles createPtRootFiles calcPtSF calcSF ratio systCheck

createRootFiles: createRootFiles.C
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -o $@ $^ $(ROOTLIBS)

createPtRootFiles: createPtRootFiles.C
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -o $@ $^ $(ROOTLIBS)

calcPtSF: calcPtSF.C
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -o $@ $^ $(ROOTLIBS)

calcSF: calcSF.C
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -o $@ $^ $(ROOTLIBS)


ratio: ratio.C
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -o $@ $^ $(ROOTLIBS)

systCheck: systCheck.C
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -o $@ $^ $(ROOTLIBS)

clean:
	rm create*RootFiles calc*SF ratio systCheck
