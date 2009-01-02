CC=gcc
CFLAGS=-I/Library/Frameworks/R.framework/Headers

default: shared_lib 

shared_lib: arlequin.o
	$(CC) arlequin.o -dynamiclib -o arlequin.so

test: arlequin.o test.o
	$(CC) arlequin.o test.o -o $@ 

test.o: arlequin.h test.c

arlequin.o: arlequin.h

clean:
	rm -f *.o *.so
	rm -f test
