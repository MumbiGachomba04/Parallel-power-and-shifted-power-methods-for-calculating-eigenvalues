#Parallel-power-and-shifted-power-methods-for-calculating-eigenvalues

This program is an implementation of the power and shifted power methods used to calculate the largest and second largest eigenvalues of a matrix respectively.
I have parallelised the program using OpenMP + OpenMPI

<b>How to run the program</b>
<br>
```export OMP_NUM_THREADS=<number of threads> ```
<br>
```mpic++ matrixPar.cc -o <binary file> -fopenmp ```
<br>
```mpirun -n <number of processes> ./<binary file> <input file> ```
