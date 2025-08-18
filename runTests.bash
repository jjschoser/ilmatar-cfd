make clean
make
export OMP_NUM_THREADS=4
./simple-cfd
python plotTests.py
