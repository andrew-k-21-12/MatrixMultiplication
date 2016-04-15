# Builds executable program without any optimizations.
default:
	clang-omp++ main.cpp -fopenmp -o main.o
	
# Builds high-performance executable program.
fast:
	clang-omp++ main.cpp -fopenmp -o main.o -O3

# Removes all unnecessary or temp files.
clean:
	rm main.o