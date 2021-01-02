CC = gcc
CFLAGS = -O2 -c -g -Wall -DBINTRAJ
LDFLAGS = -lnetcdf
LIBFLAGS = -lm -fopenmp
TARGETNAME = "hydra_site.x"

SOURCES = input_proc.c grid_proc.c main.c utils.c AmberNetcdf.c xdrfile.c xdrfile_xtc.c

OBJECTS = $(SOURCES:.c=.o)

all: $(TARGETNAME)

$(TARGETNAME): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(LIBFLAGS) -o $@

.c.o:
	$(CC) $(LDFLAGS) $(CFLAGS) $< $(LIBFLAGS) -o $@

clean:
	/bin/rm -f $(OBJECTS) $(TARGETNAME)
