#!/bin/bash

# Testing setup
reruns=1
it=500
ra=1000.0
N1=200
N2=200
N3=100
filename='testresults.csv'

# Create necessary folders
echo "Creating folders..."
mkdir -p ./build
mkdir -p ./out
mkdir -p ./in
mkdir -p ./expected

echo "Command,Reruns,Average Runtime,Average Temperature,Max Error">$filename

# Run programs before main loop
echo "Compiling code..."
g++    --std=c++11 -O3 -Wall -o ./build/gs_secv    ./gs_optimized/main.cpp
g++    --std=c++11 -O3 -Wall -o ./build/gs_omp     ./gs_optimized_omp/main.cpp -fopenmp
mpic++ --std=c++11 -O3 -Wall -o ./build/gs_mpi     ./gs_optimized_mpi/main.cpp
mpic++ --std=c++11 -O3 -Wall -o ./build/gs_mpi_omp ./gs_optimized_mpi_omp/main.cpp -fopenmp

echo "Generating input matrices..."
python3 ./generatemat.py $N1 $N2 $N3 0.0 0.0 0.0 ./in/t.bin
python3 ./generatemat.py $N1 $N2 $N3 0.0 0.0 0.0 ./in/u.bin
python3 ./generatemat.py $N1 $N2 $N3 0.0 0.0 0.0 ./in/v.bin

echo "Computing the validation matrices..."
_=$(./build/gs_secv $ra 1.0 -1.0 $it -1 ./in/u.bin ./in/v.bin ./in/t.bin ./out/u ./out/v ./out/t)
cp ./out/t_iter_$it.bin ./expected/t_expected.bin
cp ./out/u_iter_$it.bin ./expected/u_expected.bin
cp ./out/v_iter_$it.bin ./expected/v_expected.bin

echo "Finished preparing the environment!"

# List of programs and their arguments
declare -a commands=(
    "./build/gs_secv $ra 1.0 -1.0 $it -1 ./in/u.bin ./in/v.bin ./in/t.bin ./out/u ./out/v ./out/t"
    "./build/gs_omp  $ra 1.0 -1.0 $it -1 ./in/u.bin ./in/v.bin ./in/t.bin ./out/u ./out/v ./out/t 2"
    "./build/gs_omp  $ra 1.0 -1.0 $it -1 ./in/u.bin ./in/v.bin ./in/t.bin ./out/u ./out/v ./out/t 4"
    "./build/gs_omp  $ra 1.0 -1.0 $it -1 ./in/u.bin ./in/v.bin ./in/t.bin ./out/u ./out/v ./out/t 8"
    # "./build/gs_omp  $ra 1.0 -1.0 $it -1 ./in/u.bin ./in/v.bin ./in/t.bin ./out/u ./out/v ./out/t 10"
    # "./build/gs_omp  $ra 1.0 -1.0 $it -1 ./in/u.bin ./in/v.bin ./in/t.bin ./out/u ./out/v ./out/t 16"
    # "./build/gs_omp  $ra 1.0 -1.0 $it -1 ./in/u.bin ./in/v.bin ./in/t.bin ./out/u ./out/v ./out/t 20"
    "mpirun -n 4  --hostfile ./hostfile.txt ./build/gs_mpi $ra 1.0 -1.0 $it -1 2 2 1 ./in/u.bin ./in/v.bin ./in/t.bin ./out/u ./out/v ./out/t"
    "mpirun -n 8  --hostfile ./hostfile.txt ./build/gs_mpi $ra 1.0 -1.0 $it -1 2 2 2 ./in/u.bin ./in/v.bin ./in/t.bin ./out/u ./out/v ./out/t"
    # "mpirun -n 18 --hostfile ./hostfile.txt ./build/gs_mpi $ra 1.0 -1.0 $it -1 3 3 2 ./in/u.bin ./in/v.bin ./in/t.bin ./out/u ./out/v ./out/t"
    # "mpirun -n 27 --hostfile ./hostfile.txt ./build/gs_mpi $ra 1.0 -1.0 $it -1 3 3 3 ./in/u.bin ./in/v.bin ./in/t.bin ./out/u ./out/v ./out/t"
    "mpirun -n 4  --hostfile ./hostfile.txt ./build/gs_mpi_omp $ra 1.0 -1.0 $it -1 2 2 1 ./in/u.bin ./in/v.bin ./in/t.bin ./out/u ./out/v ./out/t 4"
    # "mpirun -n 4  --hostfile ./hostfile.txt ./build/gs_mpi_omp $ra 1.0 -1.0 $it -1 2 2 1 ./in/u.bin ./in/v.bin ./in/t.bin ./out/u ./out/v ./out/t 10"
    "mpirun -n 8  --hostfile ./hostfile.txt ./build/gs_mpi_omp $ra 1.0 -1.0 $it -1 2 2 2 ./in/u.bin ./in/v.bin ./in/t.bin ./out/u ./out/v ./out/t 2"
    # "mpirun -n 8  --hostfile ./hostfile.txt ./build/gs_mpi_omp $ra 1.0 -1.0 $it -1 2 2 2 ./in/u.bin ./in/v.bin ./in/t.bin ./out/u ./out/v ./out/t 10"
)

run_prog() {
    local avgruntime=0
    local avgtemp=0
    local maxerr=-1
    local count=$reruns
    
    for (( i=0; i<$reruns; i++ ))
    do
        output=$($1)
        runtime=$(echo "$output" | head -n 1)
        maxtemp=$(echo "$output" | tail -n 1)
        
        err=$(python3 compmat.py ./out/t_iter_$it.bin ./expected/t_expected.bin)
        maxerr=$(awk "BEGIN {print ($err > $maxerr) ? $err : $maxerr}")

        avgruntime=$(awk "BEGIN {print $avgruntime + $runtime}")
        avgtemp=$(awk "BEGIN {print $avgtemp + $maxtemp}")
    done
    
    avgruntime=$(awk "BEGIN {print $avgruntime / $count}")
    avgtemp=$(awk "BEGIN {print $avgtemp / $count}")
    
    echo "$1,$count,$avgruntime,$avgtemp,$maxerr">>$filename
}

# Main loop to process each command
for cmd in "${commands[@]}"
do
    echo "Running: '$cmd'"
    run_prog "$cmd"
done
