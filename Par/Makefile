CC	= g++
CFLAGS	= -std=c++14 -Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors
OPFLAGS = -O3 -DNDEBUG
EXFLAGS = -fopenmp

OBJS = par.o

all: nasteroids

nasteroids: $(OBJS)

par.o: nasteroids-par.cpp
	$(CC) $(CFLAGS) $(OPFLAGS) $(EXFLAGS) -o par.o nasteroids-par.cpp

clean:
	rm -f *.o *.txt
