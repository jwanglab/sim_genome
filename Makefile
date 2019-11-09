CC=gcc
CFLAGS=-std=c99

OBJECTS = sg

all: $(OBJECTS)

sg: main.c
	$(CC) $(CFLAGS) -lz main.c vcf.c -o sg -g

.PHONY: clean
clean:
	-rm $(OBJECTS)
