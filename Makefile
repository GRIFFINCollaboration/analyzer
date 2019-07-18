# currently set up for packages/old-versions/midas-unknown-grsmid00/linux-m32

DRV_DIR         = $(MIDASSYS)/drivers/bus
INC_DIR 	= $(MIDASSYS)/include
#LIB_DIR 	= $(MIDASSYS)/linux/lib
LIB_DIR 	= $(MIDASSYS)/linux-m32/lib

OSFLAGS = -DOS_LINUX -Dextname -DUSE_ROODY
CFLAGS = -O -m32 -Wall
#LIBS = -lm -lz -lutil -lnsl -lpthread # Now need -lrt for some reason ???
LIBS = -lm -lz -lutil -lnsl -lpthread -lrt

LIB = $(LIB_DIR)/libmidas.a

CC = gcc
CXX = g++
CFLAGS += -g -I$(INC_DIR) -I$(DRV_DIR) -I$(VME_DIR)/include
LDFLAGS += 
LDFEFLAGS += -L$(VME_DIR)/lib -lvme

MODULES 	= angrif.o andesc.o histogram.o web_server.o

all: analyzer 

analyzer: $(LIB) $(LIB_DIR)/mana.o analyzer.o $(MODULES) Makefile
	$(CC) $(CFLAGS) -o $@ $(LIB_DIR)/mana.o analyzer.o $(MODULES) \
	$(LIB) $(LDFLAGS) $(LIBS)

analyzer.o: analyzer.c
	$(CC) $(USERFLAGS) $(CFLAGS) $(OSFLAGS) -o $@ -c $<

angrif.o: angrif.c web_server.h histogram.h
	$(CC) $(USERFLAGS) $(CFLAGS) $(OSFLAGS) -o $@ -c $<

andesc.o: andesc.c web_server.h histogram.h
	$(CC) $(USERFLAGS) $(CFLAGS) $(OSFLAGS) -o $@ -c $<

web_server.o: web_server.c web_server.h histogram.h
	$(CC) $(USERFLAGS) $(CFLAGS) $(OSFLAGS) -o $@ -c $<

histogram.o: histogram.c histogram.h
	$(CC) $(USERFLAGS) $(CFLAGS) $(OSFLAGS) -o $@ -c $<

clean::
	rm -f *.o *~ \#* $(FEMODULES) $(MODULES)

#end file
