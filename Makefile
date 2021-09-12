# If the GSL and FFTW3 are not in the standard location, specify their paths here.
GSL_INCLUDE_PATH=.
GSL_LIB_PATH=.
FFTW3_INCLUDE_PATH=.
FFTW3_LIB_PATH=.

# The rest usually does not need to be modified.
CC=g++
CFLAGS=-O3 -Wall -pedantic -std=c++11 -fopenmp

INCLUDE=-I $(GSL_INCLUDE_PATH) -I $(FFTW3_INCLUDE_PATH)
LINK=-L $(GSL_LIB_PATH) -lgsl -lgslcblas -L $(FFTW3_LIB_PATH) -lfftw3_omp -lfftw3 -ffast-math -flto -march=native

DOX_NAME=Doxyfile
MAKE_NAME=Makefile
DOC_NAME=documentation

SRC_PATH=src
EXE_PATH=exe
DOC_PATH=doc

SRC_HEADERS=$(wildcard $(SRC_PATH)/*.hpp)
SRC_SOURCES=$(wildcard $(SRC_PATH)/*.cpp)
SRC_TEMPLATES=$(wildcard $(SRC_PATH)/*.tpp)
SRC_OBJECTS=$(SRC_SOURCES:.cpp=.o)
SRC_DEPENDENCIES=$(SRC_OBJECTS:.o=.d)

EXE_SOURCES=$(wildcard $(EXE_PATH)/*.cpp)
EXECUTABLES=$(EXE_SOURCES:.cpp=.x)
EXE_DEPENDENCIES=$(EXECUTABLES:.x=.d)

CLEAN_FILES=$(SRC_OBJECTS) $(SRC_DEPENDENCIES) $(EXECUTABLES) $(EXE_DEPENDENCIES)

all: $(SRC_OBJECTS) $(EXECUTABLES)

-include $(SRC_DEPENDENCIES) $(EXE_DEPENDENCIES)

%.o: %.cpp $(MAKE_NAME)
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@ $(LINK)
	@$(CC) -MM $< > $*.d
	@\sed -i"" "1s|^|$(SRC_PATH)/|" $*.d

%.x: %.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $< -o $@ $(SRC_OBJECTS) $(LINK)
	@$(CC) -MM $< > $*.d
	@\sed -i"" "1s|^|$(EXE_PATH)/|" $*.d
	@\sed -i"" "s:\.o:\.x:g" $*.d

$(DOC_PATH)/$(DOC_NAME).html: $(SRC_HEADERS) $(SRC_SOURCES) $(SRC_TEMPLATES)
	\doxygen $(DOX_NAME)
	@\ln -sf html/index.html $(DOC_PATH)/$(DOC_NAME).html

doc: $(DOC_PATH)/$(DOC_NAME).html

clean:
	\rm -f $(CLEAN_FILES)