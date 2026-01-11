make clean
make
export OMP_NUM_THREADS=4
./ilmatar-cfd
python plot_tests.py
