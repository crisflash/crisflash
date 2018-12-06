CC=gcc
RM=rm
MAKE=make
CFLAGS=-std=c99 -O3 -pthread
all:
	$(MAKE) --no-print-directory ./bin/crisflash
	$(MAKE) --no-print-directory crisflashVcf

./bin/crisflash: ./src/main_crisflash.c ./src/read.c ./src/nary_tree.c ./src/vcf.c ./src/test.c
	mkdir -p bin
	$(CC) $(CFLAGS) ./src/main_crisflash.c ./src/read.c ./src/nary_tree.c ./src/vcf.c -o ./bin/crisflash

crisflashVcf: ./src/main_vcf.c ./src/vcf.c ./src/nary_tree.c ./src/read.c
	$(CC) $(CFLAGS) ./src/main_vcf.c ./src/vcf.c ./src/read.c ./src/nary_tree.c -o ./bin/crisflashVcf

test:
	$(CC) $(CFLAGS) ./src/test.c ./src/read.c ./src/vcf.c ./src/nary_tree.c -o ./src/test

clean: ./bin/crisflash
	$(RM) ./bin/crisflash
	$(RM) ./bin/crisflashVcf
