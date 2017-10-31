.SUFFIXES:
.SUFFIXES: .o .cpp
#============================================================
TARGET	=  adsrun_temporal

C_SOURCES =  ads_graph_temporal.cpp 
C_OBJS     =  ads_graph_temporal.o  
MY_INCLUDES = graph.h eigen.h

CCX = g++
CXXFLAGS = -O2 -g -Wall -std=c++11 -I./Eigen/ -I.


#============================================================
all: $(TARGET)


ads_graph_temporal.o: $(MY_INCLUDES)

.o:.cpp	$(MY_INCLUDES)
	$(CCX)  -c  $(CXXFLAGS) $<  

$(TARGET) :   $(C_OBJS)
	$(CCX) $(CXXFLAGS)  $^ $(LIBDIRS)  -o $@

# Implicit rules: $@ = target name, $< = first prerequisite name, $^ = name of all prerequisites
#============================================================

ALL_SOURCES = Makefile $(C_SOURCES) $(MY_INCLUDES) 

NOTES =
%= otherstuff.np 

clean:
	rm -f $(TARGET) $(C_OBJS) core *.~*~

tar: $(ALL_SOURCES) $(NOTES)
	tar cvf $(TARGET).tar $(ALL_SOURCES) test.p $(NOTES)

$(TARGET).ps: $(ALL SOURCES)
	enscript -pcode.ps $(ALL_SOURCES)


