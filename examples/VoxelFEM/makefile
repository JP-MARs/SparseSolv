# MakeFile

CXX = icpx
#CXX = g++
CFLAGS = -O3 -std=c++14 -xhost -qopenmp -I../../ -I../../SparseSolv
#CFLAGS = -O3 -std=c++14 -fopenmp -I../000_BaseLibraries
OPTS = -qmkl -ipo
#OPTS = -lgomp
LIBS = -L../../ -lSparseSolv
OBJDIR = ./obj

SOURCES  = $(wildcard *.cpp) 
OBJECTS  = $(addprefix $(OBJDIR)/, $(SOURCES:.cpp=.o))

#target file
TARGET = run_VoxelFEM

$(TARGET): $(OBJECTS) 
	$(CXX) -o $@ $^ $(LDFLAGS) $(LIBS) $(OPTS)
$(OBJDIR)/%.o: %.cpp
	@if [ ! -d $(OBJDIR) ]; \
		then echo "mkdir -p $(OBJDIR) dirs"; mkdir -p $(OBJDIR); \
	fi
	@[ -d $(OBJDIR) ]
	$(CXX) -c $(CFLAGS) $(INDS) -o $@ $<

clean:
	rm $(OBJDIR)/*.o
	rm $(TARGET)
depend:
	makedepend -- -Y -- $(SOURCES)


# DO NOT DELETE

