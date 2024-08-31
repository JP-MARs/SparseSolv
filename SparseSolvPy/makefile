# MakeFile

CXX = icpx
#ARX = llvm-ar
CFLAGS = -O3 -std=c++14 -xhost -qopenmp -I./BasicDefines -I../ -fPIC  -I/opt/anaconda3/include  -I/opt/anaconda3/include/python3.12 
OPTS = -qmkl -ipo -shared 
LIBS = -L/opt/anaconda3/lib/python3.12/config-3.12-x86_64-linux-gnu -L/opt/anaconda3/lib  -lpthread -ldl  -lutil -lm


#CXX = g++
#ARX = ar
#CFLAGS = -O3 -std=c++14 -mavx2 -fopenmp -I./BasicDefines -I../ -fPIC  -I/opt/anaconda3/include  -I/opt/anaconda3/include/python3.12
#OPTS =  -lgomp -shared
#LIBS = `/opt/anaconda3/bin/python3-config --cflags --ldflags`

OBJDIR = ./obj

SOURCES  = $(wildcard *.cpp)
OBJECTS  = $(addprefix $(OBJDIR)/, $(SOURCES:.cpp=.o))

#target file
TARGET = ../SparseSolvPy.so

$(TARGET): $(OBJECTS) 
	$(CXX) -o $@ $^ $(LDFLAGS) $(LIBS) $(OPTS)
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
