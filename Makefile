CC         =g++
CFLAGS	   =-c -Wall -Ofast -fopenmp -ffast-math -ffinite-math-only
LDFLAGS	   =-Ofast -fopenmp
SOURCES	   =./fast_diffusion.cpp
OBJECTS	   =$(SOURCES:.cpp=.o)
EXECUTABLE =./fast_diffusion

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) -I $(EIGEN_PATH) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) -I $(EIGEN_PATH) $(CFLAGS) $< -o $@

clean:
	rm fast_diffusion *.o
