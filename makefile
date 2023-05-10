FLAGS= -DDEBUG
LIBS= -lm
ALWAYS_REBUILD=makefile

nbody: nbody.cu compute.cu
	nvcc $(FLAGS) $^ -o nbody
clean:
	rm -f nbody 
