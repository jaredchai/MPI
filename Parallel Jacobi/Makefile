.PHONY: release clean

CC = mpicxx
CCFLAGS = -lm
SRC = pjacobi.cpp utils.h
BIN = pjacobi


release: ${SRC}
	@echo "Checking.."
	${CC} -o ${BIN} $< ${CCFLAGS}


clean:
	@echo "Cleaning up..."
	rm -rvf *.o ${BIN}
