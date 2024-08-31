# MakeFile

CXX = icpx
ARX = llvm-ar
CFLAGS = -O3 -std=c++14 -xhost -qopenmp -I./BasicDefines -I../
OPTS = -qmkl -ipo
LIBS = -L../

#CXX = g++
#ARX = ar
#CFLAGS = -O3 -std=c++14 -mavx2 -fopenmp -I./BasicDefines -I../
#OPTS = -lgomp
#LIBS = -L../

OBJDIR = ./obj

SOURCES  = $(wildcard *.cpp)
OBJECTS  = $(addprefix $(OBJDIR)/, $(SOURCES:.cpp=.o))

#target file
TARGET = ../libSparseSolv.a

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
