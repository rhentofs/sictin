CPPFLAGS = -O3
CC = g++

BIN_DIR = bin

TARGETS = $(BIN_DIR)/make_footprint $(BIN_DIR)/make_footprint_from_BED $(BIN_DIR)/access_signal $(BIN_DIR)/build_binaries $(BIN_DIR)/minmax


##
## use either of the below to change the representation of the signals, 
##
## INT_REP for unsigned intergers,
## DBL_REP for doubles,
## SHR_REP for unsigned short, default value if not set.  
##

##DFLAGS = -DINT_REP
##DFLAGS = -DDBL_REP
DFLAGS = -DSHR_REP

all: $(TARGETS)

clean: 
	-rm -rf $(TARGETS)
	-rm -rf $(BIN_DIR)/*.o


utils.o: src/utils.cpp
	$(CC) $(CPPFLAGS) $(DFLAGS) -c src/utils.cpp -o $(BIN_DIR)/utils.o

sam_utils.o: src/sam_utils.cpp
	$(CC) $(CPPFLAGS) $(DFLAGS) -c src/sam_utils.cpp -o $(BIN_DIR)/sam_utils.o

gff_utils.o: src/gff_utils.cpp
	$(CC) $(CPPFLAGS) $(DFLAGS) -c src/gff_utils.cpp -o $(BIN_DIR)/gff_utils.o

bed_utils.o: src/bed_utils.cpp
	$(CC) $(CPPFLAGS) $(DFLAGS) -c src/bed_utils.cpp -o $(BIN_DIR)/bed_utils.o

wig_utils.o: src/wig_utils.cpp
	$(CC) $(CPPFLAGS) $(DFLAGS) -c src/wig_utils.cpp -o $(BIN_DIR)/wig_utils.o


$(BIN_DIR)/make_footprint: src/make_footprint.cpp utils.o bed_utils.o
	$(CC) $(CPPFLAGS) $(DFLAGS) src/make_footprint.cpp $(BIN_DIR)/utils.o $(BIN_DIR)/bed_utils.o -o $(BIN_DIR)/make_footprint

$(BIN_DIR)/make_footprint_from_BED: src/make_footprint_from_bed.cpp utils.o bed_utils.o 
	$(CC) $(CPPFLAGS) $(DFLAGS) src/make_footprint_from_bed.cpp $(BIN_DIR)/utils.o $(BIN_DIR)/bed_utils.o -o $(BIN_DIR)/make_footprint_from_bed


$(BIN_DIR)/access_signal: src/access_signal.cpp utils.o bed_utils.o
	$(CC) $(CPPFLAGS) $(DFLAGS) src/access_signal.cpp $(BIN_DIR)/utils.o $(BIN_DIR)/bed_utils.o -o $(BIN_DIR)/access_signal

$(BIN_DIR)/build_binaries: src/build_binaries.cpp utils.o sam_utils.o gff_utils.o bed_utils.o wig_utils.o
	$(CC) $(CPPFLAGS) $(DFLAGS) src/build_binaries.cpp $(BIN_DIR)/utils.o $(BIN_DIR)/sam_utils.o $(BIN_DIR)/gff_utils.o $(BIN_DIR)/bed_utils.o $(BIN_DIR)/wig_utils.o -o $(BIN_DIR)/build_binaries

$(BIN_DIR)/minmax: src/minmax.cpp utils.o bed_utils.o
	$(CC) $(CPPFLAGS) $(DFLAGS) src/minmax.cpp $(BIN_DIR)/utils.o $(BIN_DIR)/bed_utils.o -o $(BIN_DIR)/minmax
