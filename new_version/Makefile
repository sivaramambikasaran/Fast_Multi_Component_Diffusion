CC         =g++
CFLAGS	   =-c -Wall -Ofast -fopenmp -ffast-math -ffinite-math-only
LDFLAGS	   =-Ofast -fopenmp
SOURCES	   =./example/basic_example.cpp
OBJECTS	   =$(SOURCES:.cpp=.o)
EXECUTABLE =./exec/example

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) -I $(EIGEN_PATH) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) -I $(EIGEN_PATH) $(CFLAGS) $< -o $@

clean:
	rm ./exec/example ./example/*.o
