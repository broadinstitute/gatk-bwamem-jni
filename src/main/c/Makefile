ifeq ($(JAVA_HOME),)
$(error $$(JAVA_HOME) is not defined)
endif

JNI_INCLUDE_DIRS=$(addprefix -I,$(shell find $(JAVA_HOME)/include -type d))
CFLAGS=-ggdb -O2 -Wall -fPIC $(JNI_INCLUDE_DIRS)
CC=gcc

#OS-dependent extension lookup
UNAME := $(shell uname)
ifeq ($(UNAME),Darwin)
LIB_EXT=Darwin.dylib
else
LIB_EXT=Linux.so
endif

BWA_MEM_COMMIT=cb950614ce7217788780b9a8d445c64cd4d8f62e
JNI_BASE_NAME=org_broadinstitute_hellbender_utils_bwa_BwaMemIndex

all: libbwa.$(LIB_EXT)

libbwa.$(LIB_EXT): $(JNI_BASE_NAME).o jnibwa.o bwa/libbwa.a
	$(CC) -ggdb -dynamiclib -shared -o $@ $^ -lm -lz -lpthread

bwa:
	git clone https://github.com/lh3/bwa && cd bwa && git checkout $(BWA_MEM_COMMIT) && echo '#define BWA_COMMIT "'$(BWA_MEM_COMMIT)'"' > bwa_commit.h
	sed -i.bak -e's/\(LOBJS=.*\)/\1 bwtindex.o rle.o rope.o bwt.o is.o/g' bwa/Makefile

bwa/libbwa.a: bwa
	$(MAKE) CFLAGS="$(CFLAGS)" -C bwa libbwa.a

$(JNI_BASE_NAME).o: $(JNI_BASE_NAME).c jnibwa.h bwa

jnibwa.o: jnibwa.c jnibwa.h bwa

clean:
	rm -rf bwa *.o *.$(LIB_EXT)

.PHONY: all clean
