CC	= g++
CFLAGS	= -std=c++14 -Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors
OPFLAGS = -O3 -DNDEBUG

OBJS = seq.o 

all: nasteroids

nasteroids: $(OBJS)

seq.o: nasteroids-seq.cpp
	$(CC) $(CFLAGS) $(OPFLAGS) -o seq.o nasteroids-seq.cpp

clean:
	rm -f *.o *.txt


