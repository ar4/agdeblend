MPICC = mpicc
MPIFC = mpif90

data = data_1_blended_data_2.bin data_2_blended_data.bin data_3_blended_data_0.bin data_3_trace_types_1.bin data_4_shottimes_0.bin \
       data_1_blended_data.bin data_2_channels.bin data_3_blended_data_1.bin data_3_true_data_0.bin data_4_shottimes_1.bin \
       data_1_channels.bin data_2_shottimes.bin data_3_channels_0.bin data_3_true_data_1.bin data_4_trace_types_0.bin \
       data_1_shottimes.bin data_2_trace_types.bin data_3_channels_1.bin data_4_blended_data_0.bin data_4_trace_types_1.bin \
       data_1_shottimes_out.bin data_2_true_data.bin data_3_shottimes_0.bin data_4_blended_data_1.bin data_4_true_data_0.bin \
       data_1_trace_types.bin data_2_wavelet.bin data_3_shottimes_1.bin data_4_channels_0.bin data_4_true_data_1.bin \
       data_1_true_data.bin data_2_wavelet_idxs.bin data_3_trace_types_0.bin data_4_channels_1.bin

all: make_data c_examples cpp_examples python_examples julia_examples fortran_examples

clean:
	rm -rf *.bin out example_*_* *.o *.mod check_example_? check_example_?_? make_data

c_examples: example_1_c example_2_c example_3_c example_4_c example_5_c example_6_c example_7_c example_8_c example_9_c

cpp_examples: example_1_cpp

python_examples: example_1_py example_2_py example_3_py example_4_py example_5_py example_6_py example_7_py example_8_py example_9_py

julia_examples: example_1_jl example_2_jl example_3_jl example_4_jl example_5_jl example_6_jl example_7_jl example_8_jl example_9_jl

fortran_examples: example_1_f example_2_f example_3_f example_3_f08 example_4_f example_5_f example_6_f example_6_f08 example_7_f example_8_f example_9_f

$(data): make_data

make_data: make_data.c
	$(CC) make_data.c -I ../include -L ../lib -lagdeblend -lfftw3f -lm -o make_data
	./make_data

check_example_1: check_example.c
	$(CC) check_example.c -DFILENAME="\"data_1_blended_data.bin\"" -DN=65536 -o check_example_1

check_example_2: check_example.c
	$(CC) check_example.c -DFILENAME="\"data_1_blended_data_2.bin\"" -DN=98304 -o check_example_2

check_example_4: check_example_rmse.c
	$(CC) check_example_rmse.c -DFILENAME="\"data_1_true_data.bin\"" -DOUT_FILENAME="\"data_1_deblended_data.bin\"" -DN=65536 -lm -o check_example_4

check_example_7: check_example_rmse.c
	$(CC) check_example_rmse.c -DFILENAME="\"data_2_true_data.bin\"" -DOUT_FILENAME="\"data_2_deblended_data.bin\"" -DN=65536 -lm -o check_example_7

check_example_8_0: check_example_rmse.c
	$(CC) check_example_rmse.c -DFILENAME="\"data_3_true_data_0.bin\"" -DOUT_FILENAME="\"data_3_deblended_data_0.bin\"" -DN=32768 -DMAX_ERROR_FRAC=0.05f -lm -o check_example_8_0

check_example_8_1: check_example_rmse.c
	$(CC) check_example_rmse.c -DFILENAME="\"data_3_true_data_1.bin\"" -DOUT_FILENAME="\"data_3_deblended_data_1.bin\"" -DN=49152 -DMAX_ERROR_FRAC=0.05f -lm -o check_example_8_1

check_example_9_0: check_example_rmse.c
	$(CC) check_example_rmse.c -DFILENAME="\"data_4_true_data_0.bin\"" -DOUT_FILENAME="\"data_4_deblended_data_0.bin\"" -DN=196608 -DMAX_ERROR_FRAC=0.05f -lm -o check_example_9_0

check_example_9_1: check_example_rmse.c
	$(CC) check_example_rmse.c -DFILENAME="\"data_4_true_data_1.bin\"" -DOUT_FILENAME="\"data_4_deblended_data_1.bin\"" -DN=196608 -DMAX_ERROR_FRAC=0.05f -lm -o check_example_9_1

example_1_c: $(data) check_example_1
	mkdir -p out
	rm -f out/*
	$(CC) example_1.c -I ../include -L ../lib -lagdeblend -o example_1_c
	./example_1_c
	./check_example_1

example_2_c: $(data) check_example_2
	mkdir -p out
	rm -f out/*
	$(CC) example_2.c -I ../include -L ../lib -lagdeblend -o example_2_c
	./example_2_c
	./check_example_2

example_3_c: $(data) check_example_1
	mkdir -p out
	rm -f out/*
	$(MPICC) example_3.c -DAGD_MPI -I ../include -L ../lib -lagdeblend_mpi -o example_3_c
	mpirun -np 2 ./example_3_c
	./check_example_1

example_4_c: $(data) check_example_4
	mkdir -p out
	rm -f out/*
	$(CC) example_4.c -I ../include -L ../lib -lagdeblend -o example_4_c
	./example_4_c
	./check_example_4

example_5_c: $(data) check_example_4
	mkdir -p out
	rm -f out/*
	$(CC) example_5.c -I ../include -L ../lib -lagdeblend -o example_5_c
	./example_5_c
	./check_example_4

example_6_c: $(data) check_example_4
	mkdir -p out
	rm -f out/*
	$(MPICC) example_6.c -DAGD_MPI -I ../include -L ../lib -lagdeblend_mpi -o example_6_c
	mpirun -np 2 ./example_6_c
	./check_example_4

example_7_c: $(data) check_example_7
	mkdir -p out
	rm -f out/*
	$(CC) example_7.c -I ../include -L ../lib -lagdeblend -o example_7_c
	./example_7_c
	./check_example_7

example_8_c: $(data) check_example_8_0 check_example_8_1
	mkdir -p out
	rm -f out/*
	$(CC) example_8.c -I ../include -L ../lib -lagdeblend -o example_8_c
	./example_8_c
	./check_example_8_0
	./check_example_8_1

example_9_c: $(data) check_example_9_0 check_example_9_1
	mkdir -p out
	rm -f out/*
	$(CC) example_9.c -I ../include -L ../lib -lagdeblend -o example_9_c
	./example_9_c
	./check_example_9_0
	./check_example_9_1

example_1_cpp: $(data) check_example_1
	mkdir -p out
	rm -f out/*
	$(CXX) example_1.cpp -I ../include -L ../lib -lagdeblend -o example_1_cpp
	./example_1_cpp
	./check_example_1

example_1_py: $(data) check_example_1
	mkdir -p out
	rm -f out/*
	python3 example_1.py
	./check_example_1

example_2_py: $(data) check_example_2
	mkdir -p out
	rm -f out/*
	python3 example_2.py
	./check_example_2

example_3_py: $(data) check_example_1
	mkdir -p out
	rm -f out/*
	mpirun -np 2 python3 example_3.py
	./check_example_1

example_4_py: $(data) check_example_4
	mkdir -p out
	rm -f out/*
	python3 example_4.py
	./check_example_4

example_5_py: $(data) check_example_4
	mkdir -p out
	rm -f out/*
	python3 example_5.py
	./check_example_4

example_6_py: $(data) check_example_4
	mkdir -p out
	rm -f out/*
	mpirun -np 2 python3 example_6.py
	./check_example_4

example_7_py: $(data) check_example_7
	mkdir -p out
	rm -f out/*
	python3 example_7.py
	./check_example_7

example_8_py: $(data) check_example_8_0 check_example_8_1
	mkdir -p out
	rm -f out/*
	python3 example_8.py
	./check_example_8_0
	./check_example_8_1

example_9_py: $(data) check_example_9_0 check_example_9_1
	mkdir -p out
	rm -f out/*
	python3 example_9.py
	./check_example_9_0
	./check_example_9_1

example_1_jl: $(data) check_example_1
	mkdir -p out
	rm -f out/*
	julia example_1.jl
	./check_example_1

example_2_jl: $(data) check_example_2
	mkdir -p out
	rm -f out/*
	julia example_2.jl
	./check_example_2

example_3_jl: $(data) check_example_1
	mkdir -p out
	rm -f out/*
	mpirun -np 2 julia example_3.jl
	./check_example_1

example_4_jl: $(data) check_example_4
	mkdir -p out
	rm -f out/*
	julia example_4.jl
	./check_example_4

example_5_jl: $(data) check_example_4
	mkdir -p out
	rm -f out/*
	julia example_5.jl
	./check_example_4

example_6_jl: $(data) check_example_4
	mkdir -p out
	rm -f out/*
	mpirun -np 2 julia example_6.jl
	./check_example_4

example_7_jl: $(data) check_example_7
	mkdir -p out
	rm -f out/*
	julia example_7.jl
	./check_example_7

example_8_jl: $(data) check_example_8_0 check_example_8_1
	mkdir -p out
	rm -f out/*
	julia example_8.jl
	./check_example_8_0
	./check_example_8_1

example_9_jl: $(data) check_example_9_0 check_example_9_1
	mkdir -p out
	rm -f out/*
	julia example_9.jl
	./check_example_9_0
	./check_example_9_1

example_1_f: $(data) check_example_1
	mkdir -p out
	rm -f out/*
	$(CC) ../wrappers/fortran/agdeblend_c.c -c -o agdeblend_c.o -I ../include/
	$(FC) agdeblend_c.o ../wrappers/fortran/agdeblend.F90 example_1.f90 -L ../lib/ -lagdeblend -o example_1_f
	./example_1_f
	./check_example_1

example_2_f: $(data) check_example_2
	mkdir -p out
	rm -f out/*
	$(CC) ../wrappers/fortran/agdeblend_c.c -c -o agdeblend_c.o -I ../include/
	$(FC) agdeblend_c.o ../wrappers/fortran/agdeblend.F90 example_2.f90 -L ../lib/ -lagdeblend -o example_2_f
	./example_2_f
	./check_example_2

example_3_f: $(data) check_example_1
	mkdir -p out
	rm -f out/*
	$(MPICC) ../wrappers/fortran/agdeblend_c.c -c -DAGD_MPI -o agdeblend_c.o -I ../include/
	$(MPIFC) agdeblend_c.o ../wrappers/fortran/agdeblend.F90 example_3.F90 -DAGD_MPI -L ../lib/ -lagdeblend_mpi -o example_3_f
	mpirun -np 2 ./example_3_f
	./check_example_1

example_3_f08: $(data) check_example_1
	mkdir -p out
	rm -f out/*
	$(MPICC) ../wrappers/fortran/agdeblend_c.c -c -DAGD_MPI -o agdeblend_c.o -I ../include/
	$(MPIFC) agdeblend_c.o ../wrappers/fortran/agdeblend.F90 example_3.F90 -DAGD_MPI -DAGD_F08 -L ../lib/ -lagdeblend_mpi -o example_3_f
	mpirun -np 2 ./example_3_f
	./check_example_1

example_4_f: $(data) check_example_4
	mkdir -p out
	rm -f out/*
	$(CC) ../wrappers/fortran/agdeblend_c.c -c -o agdeblend_c.o -I ../include/
	$(FC) agdeblend_c.o ../wrappers/fortran/agdeblend.F90 example_4.f90 -L ../lib/ -lagdeblend -o example_4_f
	./example_4_f
	./check_example_4

example_5_f: $(data) check_example_4
	mkdir -p out
	rm -f out/*
	$(CC) ../wrappers/fortran/agdeblend_c.c -c -o agdeblend_c.o -I ../include/
	$(FC) agdeblend_c.o ../wrappers/fortran/agdeblend.F90 example_5.f90 -L ../lib/ -lagdeblend -o example_5_f
	./example_5_f
	./check_example_4

example_6_f: $(data) check_example_4
	mkdir -p out
	rm -f out/*
	$(MPICC) ../wrappers/fortran/agdeblend_c.c -c -DAGD_MPI -o agdeblend_c.o -I ../include/
	$(MPIFC) agdeblend_c.o ../wrappers/fortran/agdeblend.F90 example_6.F90 -DAGD_MPI -L ../lib/ -lagdeblend_mpi -o example_6_f
	mpirun -np 2 ./example_6_f
	./check_example_4

example_6_f08: $(data) check_example_4
	mkdir -p out
	rm -f out/*
	$(MPICC) ../wrappers/fortran/agdeblend_c.c -c -DAGD_MPI -o agdeblend_c.o -I ../include/
	$(MPIFC) agdeblend_c.o ../wrappers/fortran/agdeblend.F90 example_6.F90 -DAGD_MPI -DAGD_F08 -L ../lib/ -lagdeblend_mpi -o example_6_f
	mpirun -np 2 ./example_6_f
	./check_example_4

example_7_f: $(data) check_example_7
	mkdir -p out
	rm -f out/*
	$(CC) ../wrappers/fortran/agdeblend_c.c -c -o agdeblend_c.o -I ../include/
	$(FC) agdeblend_c.o ../wrappers/fortran/agdeblend.F90 example_7.f90 -L ../lib/ -lagdeblend -o example_7_f
	./example_7_f
	./check_example_7

example_8_f: $(data) check_example_8_0 check_example_8_1
	mkdir -p out
	rm -f out/*
	$(CC) ../wrappers/fortran/agdeblend_c.c -c -o agdeblend_c.o -I ../include/
	$(FC) agdeblend_c.o ../wrappers/fortran/agdeblend.F90 example_8.f90 -L ../lib/ -lagdeblend -o example_8_f
	./example_8_f
	./check_example_8_0
	./check_example_8_1

example_9_f: $(data) check_example_9_0 check_example_9_1
	mkdir -p out
	rm -f out/*
	$(CC) ../wrappers/fortran/agdeblend_c.c -c -o agdeblend_c.o -I ../include/
	$(FC) agdeblend_c.o ../wrappers/fortran/agdeblend.F90 example_9.f90 -L ../lib/ -lagdeblend -o example_9_f
	./example_9_f
	./check_example_9_0
	./check_example_9_1

.PHONY: clean c_examples python_examples julia_examples fortran_examples make_data example_1_c example_2_c example_3_c example_4_c example_5_c example_6_c example_7_c example_8_c example_9_c example_1_py example_2_py example_3_py example_4_py example_5_py example_6_py example_7_py example_8_py example_9_py example_1_jl example_2_jl example_3_jl example_4_jl example_5_jl example_6_jl example_7_jl example_8_jl example_9_jl example_1_f example_2_f example_3_f example_3_f08 example_4_f example_5_f example_6_f example_6_f08 example_7_f example_8_f example_9_f
