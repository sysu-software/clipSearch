CXXC=g++
LIBS=-lm -lz
CFLAGS = -O3 -g
HG_DEFS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE
HG_WARN=-Wformat -Wreturn-type
UTILITIES_DIR = ../../thirdUtils
BIN_DIR = ../../bin
BIO_DIR = ../../bioUtils
INCLUDES = -I$(UTILITIES_DIR)/RNAfoldLib \
           -I$(BIO_DIR)
BIO_LIBS   = -L$(UTILITIES_DIR)/RNAfoldLib/ -lRNAfold \
             -L$(BIO_DIR)/ -lbiotools \

clipSearch: clipSearch.o clipSearchMain.o
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -o ${BIN_DIR}/clipSearch clipSearchMain.o clipSearch.o \
	$(BIO_LIBS) $(LIBS) 

clipSearch.o: clipSearch.cpp clipSearch.h
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -c clipSearch.cpp
	
clipSearchMain.o: clipSearchMain.cpp clipSearch.h
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -c clipSearchMain.cpp
	
clean:
	rm -f *.o
