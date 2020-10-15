CC = gcc -g -Wall

CFLAGS = -I$(FLINT_INCLUDE)
LDFLAGS = 
LDLIBS = -lmpir -lflint

genus4_search: points.o fields.o curves.o format_time.o

points.o: genus4.h

fields.o: genus4.h

curves.o: genus4.h

format_time.o: 

search2: genus4_search
	./genus4_search 2 genus4_GF2.data

search3: genus4_search
	./genus4_search 3 genus4_GF3.data

search4: genus4_search
	./genus4_search 4 genus4_GF4.data

clean:
	rm -rf *.o
	rm -rf *.dSYM
	rm -f genus4_search
