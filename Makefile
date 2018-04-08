#Your path to Eigen
EIGEN=/usr/include/eigen3

#Your path to GSL
GSL=/scratch/CPviolator/gsl-2.4
GSL_LIBS= -lgsl -lgslcblas -lm

#Your path to CUDA
CUDA=/usr/local/cuda-9.0
CUDA_LIBS= -lcuda -lcurand -lcudart 


TARGET	 = adsrun

SOURCES  = main.cpp util.cpp cg.cpp hyp_util.cpp graph.cpp data_io.cpp data_proc.cpp mcPhiFourth2D.cpp phiFourth2D.cpp mcIsing2D.cpp gpuMcPhiFourth2D.cu
OBJS     = main.o util.o cg.o hyp_util.o graph.o data_io.o data_proc.o mcPhiFourth2D.o phiFourth2D.o mcIsing2D.o ising2D.o gpuMcPhiFourth2D.o
INCLUDES = util.h graph.h cg.h eigen.h mcPhiFourth2D.h phiFourth2D.h mcIsing2D.h ising2D.h monte_carlo_ads_local.h gpuMcPhiFourth2D.cuh

LIBS= -L${GSL} ${GSL_LIBS} -L${CUDA}/lib64 ${CUDA_LIBS}

ERRS=-Wall -Wno-sign-compare -Wno-int-in-bool-context -Wno-unused-but-set-variable -Wno-unknown-warning-option

CXX = g++
OMP_FLAG= -DUSE_OMP -fopenmp
GPU_FLAG= -DUSE_GPU
CXXFLAGS = -O3 -g -Wall -std=c++11  -I. -I${EIGEN} -I${GSL} -I${CUDA}/include ${ERRS} ${OMP_FLAG} ${GPU_FLAG} 

#-maxrregcount=38
NVCC = /usr/local/cuda-9.0/bin/nvcc 
NVCCFLAGS = -O3 -std=c++11  -I. -I${EIGEN} -I${GSL} -I${CUDA}/include ${GPU_FLAG} -arch compute_60 

#============================================================

all: $(TARGET)

${TARGET}: ${OBJS}
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

main.o: main.cpp $(INCLUDES) 
	${CXX} ${CXXFLAGS} -c main.cpp

graph.o: graph.cpp $(INCLUDES) 
	${CXX} ${CXXFLAGS} -c graph.cpp

util.o: util.cpp $(INCLUDES) 
	${CXX} ${CXXFLAGS} -c util.cpp

hyp_util.o: hyp_util.cpp $(INCLUDES) 
	${CXX} ${CXXFLAGS} -c hyp_util.cpp

data_io.o: data_io.cpp $(INCLUDES) 
	${CXX} ${CXXFLAGS} -c data_io.cpp

data_proc.o: data_proc.cpp $(INCLUDES) 
	${CXX} ${CXXFLAGS} -c data_proc.cpp

cg.o: cg.cpp $(INCLUDES) 
	${CXX} ${CXXFLAGS} -c cg.cpp

ising2D.o: ising2D.cpp $(INCLUDES) 
	${CXX} ${CXXFLAGS} -c ising2D.cpp

mcIsing2D.o: mcIsing2D.cpp $(INCLUDES) 
	${CXX} ${CXXFLAGS} -c mcIsing2D.cpp

mcPhiFourth2D.o: mcPhiFourth2D.cpp $(INCLUDES) 
	${CXX} ${CXXFLAGS} -c mcPhiFourth2D.cpp

phiFourth2D.o: phiFourth2D.cpp $(INCLUDES) 
	${CXX} ${CXXFLAGS} -c phiFourth2D.cpp

monte_carlo_ads_local.o: monte_carlo_ads_local.cpp $(INCLUDES) 
	${CXX} ${CXXFLAGS} -c monte_carlo_ads_local.cpp 

gpuMcPhiFourth2D.o: gpuMcPhiFourth2D.cu $(INCLUDES) 
	${NVCC} ${NVCCFLAGS} -c gpuMcPhiFourth2D.cu

ALL_SOURCES = Makefile $(SOURCES) $(INCLUDES) 

clean:
	rm -f $(TARGET) $(OBJS) core *.~*~

tar: $(ALL_SOURCES) $(NOTES)
	tar cvf $(TARGET).tar $(ALL_SOURCES) $(NOTES)
