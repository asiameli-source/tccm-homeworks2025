CC=gcc
SRC=MP2_energy.c main.c
CFLAGS = -I/usr/local/include
LDFLAGS = -L/usr/local/lib
LDLIBS = -ltrexio
OBJ=$(patsubst %.c, %.o, $(SRC))


MP2_energy: $(OBJ)
	$(CC) -o MP2_energy $(OBJ) $(LDFLAGS) $(LDLIBS)

%.o: %.c
	$(CC) -c $< $(CFLAGS)

MP2_energy.o: MP2_energy.h

clean:
	rm -f $(OBJ) MP2_energy
