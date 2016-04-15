# Builds all flavors and runs evaluation experiments.
all: default fast
	./default.o
	./fast.o log_fast_

# Builds executable program without any optimizations.
default:
	clang-omp++ main.cpp -fopenmp -o default.o
	
# Builds high-performance executable program.
fast:
	clang-omp++ main.cpp -fopenmp -o fast.o -O3

# Removes all unnecessary or temp files.
clean:
	rm *.o