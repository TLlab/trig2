# Makefile

CC     := g++
CFLAGS := -O2 -Wall -std=c99
OBJ    := ProcessAlignment.o delta.o fastx_read.o extractCDR3.o vdjreader.o
EXE    := ProcessAlignment

$(EXE):$(OBJ)

.PHONY:
	clean

clean:
	rm -f $(OBJ) $(EXE)
