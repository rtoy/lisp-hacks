# Simple makefile to build a shared librar of the reduction routines.

CFLAGS = -m32 -fPIC
LDFLAGS = -m32

libred.so : e_rem_pio2.o k_rem_pio2.o rem_pi.o
	$(CC) -shared $(LDFLAGS) -o $@ $^

clean :
	$(RM) *.o
	$(RM) libred.so
