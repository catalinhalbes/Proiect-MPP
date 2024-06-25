#!/bin/bash
./compile.sh
mkdir -p ./out/t
mkdir -p ./out/u
mkdir -p ./out/v
mpirun -n 27 --hostfile hosts.txt ./build/main 1000.0 1.0 -1.0 2000 50 3 3 3 ./in/u.bin ./in/v.bin ./in/t.bin ./out/u/u ./out/v/v ./out/t/t
