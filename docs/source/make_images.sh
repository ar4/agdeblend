#! /bin/sh

cd ../../examples
python example_1.py
python plot_for_docs.py 1
python example_2.py
python plot_for_docs.py 2
mpirun -np 2 python example_3.py
python plot_for_docs.py 3
python example_4.py
python plot_for_docs.py 4
python example_5.py
python plot_for_docs.py 5
mpirun -np 2 python example_6.py
python plot_for_docs.py 6
python example_7.py
python plot_for_docs.py 7
python example_8.py
python plot_for_docs.py 8
#python example_9.py
#python plot_for_docs.py 9
mv *.jpg *.png ../docs/source/
