OMP_CXX = /usr/local/opt/llvm/bin/clang++
MPI_CXX = mpic++

MPI_SRC = rr.cpp

all: rr

rr: $(MPI_SRC)
	$(MPI_CXX) -std=c++17 -O2 $(MPI_SRC) -o rr

clean:
	rm -f rr
