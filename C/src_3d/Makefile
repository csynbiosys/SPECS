# For chemotaxis simulation
# F.Tostevin 10/07/2010
# modified by m. flores 12/07/2010 to add load_params function

# required system libs
SYS_LIBS=-lgsl -lgslcblas 
CFLAGS=-g0 -O2 -c
LDFLAGS=
CC= g++

EXECUTABLE=chemo
SOURCES=chemo.cpp cell.cpp load_params.cpp
OBJECTS=$(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	${CC} $(OBJECTS) ${LDFLAGS} ${SYS_LIBS} -o $@

.cpp.o:
	${CC} ${CFLAGS} $< -o $@

clean:
	rm -f *.o $(EXECUTABLE)
