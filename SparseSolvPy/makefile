# MakeFile

CXX = icpx
ARX = llvm-ar
CFLAGS = -O3 -std=c++14 -xhost -qopenmp -I./BasicDefines -I../ -fPIC  -I/opt/intel/oneapi/intelpython/latest/include  -I/opt/intel/oneapi/intelpython/latest/include/python3.9
OPTS = -qmkl -ipo -shared
LIBS = -L`/opt/intel/oneapi/intelpython/latest/bin/python3-config --cflags --ldflags` -L/opt/intel/oneapi/intelpython/latest/lib/

#CXX = g++
#ARX = ar
#CFLAGS = -O3 -std=c++14 -mavx2 -fopenmp -I./BasicDefines -I../ -fPIC
#OPTS = -lgomp -shared
#LIBS = -L../

OBJDIR = ./obj

SOURCES  = $(wildcard *.cpp)
OBJECTS  = $(addprefix $(OBJDIR)/, $(SOURCES:.cpp=.o))

#target file
TARGET = ../SparseJMSPy.so

$(TARGET): $(OBJECTS) 
	$(ARX) rcs $@ $^
$(OBJDIR)/%.o: %.cpp
	@if [ ! -d $(OBJDIR) ]; \
		then echo "mkdir -p $(OBJDIR)"; mkdir -p $(OBJDIR); \
	fi
	@[ -d $(OBJDIR) ]
	$(CXX) -c $(CFLAGS) $(INDS) -o $@ $<

clean:
	rm $(OBJDIR)/*.o
	rm $(TARGET)
depend:
	makedepend -- -Y -- $(SOURCES)


# DO NOT DELETE
