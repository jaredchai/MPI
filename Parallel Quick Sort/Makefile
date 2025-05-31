.PHONY: release clean

CC = mpicxx
CCFLAGS = -lm
SRC = pqsort.cpp
BIN = pqsort


release: ${SRC}
	@echo "Checking.."
	${CC} -o ${BIN} $< ${CCFLAGS}


clean:
	@echo "Cleaning up..."
	rm -rvf *.o ${BIN}
