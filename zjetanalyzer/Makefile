include Makefile.arch
#----------------------------------------------------------------------------------
#TESTS           := parton_analyzer.C
#TESTS           := zjetanalyzer.C
TESTS           := zjanalyzr.C
#TESTS           := npileup.C
#TESTS           := parton_analyzer.C
TESTO           := $(TESTS:.C=.o)
#TEST            := parton_analyzer
#TEST            := zjetanalyzer
TEST            := zjanalyzr
#TEST            := npileup
#TEST            := parton_analyzer

OBJ             = $(TESTO)
PROGRAMS        = $(TEST)

CXXFLAGS      += 
CXXFLAGS      += 
GLIBS         += 
#----------------------------------------------------------------------------------
.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(PROGRAMS)

$(TEST):        $(TESTO)
		$(LD) $(LDFLAGS) $^ $(GLIBS) $(OutPutOpt)$@
		@echo "$@ done"

clean:
		@rm -f $(OBJ) $(TESTO) core *~ src/*~ include/*~

.$(SrcSuf).$(ObjSuf):
		$(CXX) $(CXXFLAGS) -c $< -o $@
