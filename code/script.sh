echo "Problem size is set in source files. By default problem size = 10^7 integers."
echo "Sample sort is complilation:"
mpicc samplesort.c -o samplesort
echo "Sample sort :"
mpiexec -n 2 ./samplesort
mpiexec -n 4 ./samplesort
mpiexec -n 8 ./samplesort
mpiexec -n 12 ./samplesort
mpiexec -n 16 ./samplesort
echo "Quicksort:"
gcc qsort.c -o qsort
./qsort