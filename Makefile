include include.mk

INC_PARAMS=$(foreach d, $(INCLUDE_DIRS), -I$d)

EXE = bin/mpfm
EXE_DEBUG = mpfm_debug
TEST_EXE = test/run_tests

CC = ${CXX}
CFLAGS = -std=c++11 -Wall $(INC_PARAMS) -O3 -g -ligraph

DEBUG_FLAG = -D __DEBUG__


.PHONY: run test debug

all:  clean $(EXE)
	@echo  Compiled successfully to $(EXE)

debug: $(EXE_DEBUG)
	@echo Compiled successfully with __DEBUG__ flag to $(EXE)

$(EXE): $(OBJS)
	@echo Building mpfm...
	@echo
	@$(CC) $(CFLAGS) -o $(EXE) $(OBJS)

OBJS = $(SRCS:.c=.o) $(MAIN:.c=.o)

.c.o:
	$(CC) $(CFLAGS) -c $<  -o $@

clean:
	$(RM) *.o *~ $(EXE)

run: $(EXE)
	./${EXE}
