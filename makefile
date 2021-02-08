EXEC=k-center
CFLAGS=-ansi -Wall -Wextra -pthread -Wconversion
OFLAGS=--std=c++17 -no-pie -g -O2 -DNDEBUG -ldl
LDFLAGS=-lm -L/usr/lib/x86_64-linux-gnu -pthread
LIBFLAGS=-D_REENTRANT 
CC=g++
BIN=bin/
SRC=src
HDR=headers
VPATH=$(SRC):$(HDR)


.PHONY: clean mrproper

default: $(EXEC)


$(BIN)main.o: main.c query.h point.h utils.h algo_fully_adv.h algo_lp_adv.h algo_disk_adv.h data_fully_adv.h set.h lookup.h

$(BIN)algo_fully_adv.o: algo_fully_adv.c query.h point.h utils.h algo_fully_adv.h data_fully_adv.h set.h

$(BIN)algo_lp_adv.o: algo_lp_adv.c query.h point.h utils.h algo_lp_adv.h data_fully_adv.h gurobi_c++.h

$(BIN)algo_disk_adv.o: algo_disk_adv.c query.h point.h utils.h algo_disk_adv.h data_fully_adv.h

$(BIN)point.o: point.c point.h

$(BIN)data_fully_adv.o:  data_fully_adv.c data_fully_adv.h point.h utils.h

$(BIN)query.o: query.c query.h utils.h point.h set.h lookup.h

$(BIN)set.o: set.c set.h utils.h

$(BIN)lookup.o: lookup.c lookup.h utils.h

$(EXEC): $(BIN)main.o $(BIN)algo_fully_adv.o $(BIN)algo_lp_adv.o $(BIN)algo_disk_adv.o $(BIN)query.o $(BIN)utils.o $(BIN)point.o $(BIN)data_fully_adv.o $(BIN)set.o $(BIN)lookup.o
	$(CC) -o $@ $^ lib/libgurobi_c++.a lib/libgurobi90.so  $(CFLAGS) $(LDFLAGS) $(OFLAGS)

$(BIN)utils.o: utils.c utils.h

clean: 
	rm -f $(BIN)*.o
	rm -f src/*~

mrproper: clean
	rm -f $(EXEC)



$(BIN)%.o : %.c 
	$(CC) -c $< $(CFLAGS) $(OFLAGS) $(LDFLAGS) -o $@
	
