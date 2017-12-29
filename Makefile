.SUFFIXES:
.SUFFIXES: .o .cpp
#============================================================

#Your path to Eigen
EIGEN=/usr/local/include/eigen3/Eigen

TARGET	    = adsrun
C_SOURCES   = ads_graph.cpp 
C_OBJS      = ads_graph.o  
MY_INCLUDES = graph.h util.h cg.h cg_multishift.h eigen.h graph.h update.h 

CCX = g++ -std=c++11
CXXFLAGS = -O2 -g -Wall -std=c++11 -I${EIGEN} -I. -Wall -Wno-sign-compare

#============================================================
all: $(TARGET)

#ads_graph.o: graph.h util.h cg.h cg_multishift.h eigen.h graph.h update.h

ads_graph.o: $(MY_INCLUDES)

.o:.cpp	$(MY_INCLUDES)
	$(CCX) -c $(CXXFLAGS) $<

$(TARGET) : $(C_OBJS)
	$(CCX) $(CXXFLAGS)  $^ $(LIBDIRS)  -o $@

# Implicit rules: $@ = target name, $< = first prerequisite name, $^ = name of all prerequisites
#============================================================

ALL_SOURCES = Makefile $(C_SOURCES) $(MY_INCLUDES) 

NOTES =
%= otherstuff.np 

clean:
	rm -f $(TARGET) $(C_OBJS) core *.~*~

tar: $(ALL_SOURCES) $(NOTES)
	tar cvf $(TARGET).tar $(ALL_SOURCES) $(NOTES)

$(TARGET).ps: $(ALL SOURCES)
	enscript -pcode.ps $(ALL_SOURCES)
