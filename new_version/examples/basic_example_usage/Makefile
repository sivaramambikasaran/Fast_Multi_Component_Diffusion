CC         =g++
CFLAGS	   =-c -Wall -Ofast -fopenmp -ffast-math -ffinite-math-only
LDFLAGS	   =-Ofast -fopenmp
SOURCES	   =./main.cpp
OBJECTS	   =$(SOURCES:.cpp=.o)
EXECUTABLE =./example

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) -I $(EIGEN_PATH) -I ./../../ $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) -I $(EIGEN_PATH) -I ./../../ $(CFLAGS) $< -o $@

clean:
	rm ./example ./*.o
