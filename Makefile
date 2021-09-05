# Set these values to point to your installation of Google Test
GOOGLETEST=~/tmp/googletest/googletest
GOOGLETESTLIB=~/tmp/googletest/mybuild/lib

all: threads mpi

serial: float-serial double-serial

threads: float-threads double-threads

mpi: float-mpi double-mpi

float-serial:
	mkdir -p lib
	gcc -Wall -Wextra -pedantic -shared -fPIC -O2 src/agdeblend.c -Iinclude -lfftw3f -lm -o lib/libagdeblend.so

float-threads:
	mkdir -p lib
	gcc -Wall -Wextra -pedantic -shared -fPIC -O2 -DAGD_THREADS src/agdeblend.c -Iinclude -lfftw3f -fopenmp -lm -o lib/libagdeblend.so

double-serial:
	mkdir -p lib
	gcc -Wall -Wextra -pedantic -shared -fPIC -O2 -DAGD_DOUBLE src/agdeblend.c -Iinclude -lfftw3 -lm -o lib/libagdeblend_double.so

double-threads:
	mkdir -p lib
	gcc -Wall -Wextra -pedantic -shared -fPIC -O2 -DAGD_DOUBLE -DAGD_THREADS src/agdeblend.c -Iinclude -lfftw3 -fopenmp -lm -o lib/libagdeblend_double.so

float-mpi:
	mkdir -p lib
	mpicc -Wall -Wextra -pedantic -shared -fPIC -O2 -DAGD_MPI src/agdeblend.c -Iinclude -lfftw3f -lm -o lib/libagdeblend_mpi.so

double-mpi:
	mkdir -p lib
	mpicc -Wall -Wextra -pedantic -shared -fPIC -O2 -DAGD_DOUBLE -DAGD_MPI src/agdeblend.c -Iinclude -lfftw3 -lm -o lib/libagdeblend_mpi_double.so

clean:
	rm -rf lib

test:
	g++ -Wall -Wextra -pedantic src/test.cpp $(GOOGLETEST)/src/gtest_main.cc -DAGD_DOUBLE -Iinclude -I $(GOOGLETEST)/include/ -L $(GOOGLETESTLIB) -lgtest -lpthread -lfftw3 -lm -O2
	./a.out
	rm ./a.out

test-mpi:
	mpic++ -Wall -Wextra -pedantic src/mpi_test.cpp $(GOOGLETEST)/src/gtest_main.cc -DAGD_DOUBLE -DAGD_MPI -Iinclude -I $(GOOGLETEST)/include/ -L $(GOOGLETESTLIB) -lgtest -lpthread -lfftw3 -lm -O2
	mpirun -n 2 ./a.out
	rm ./a.out
