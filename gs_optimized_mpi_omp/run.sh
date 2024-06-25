#!/bin/bash
./compile.sh
mpirun -n 4 --hostfile hosts.txt ./build/main 1000.0 1.0 -1.0 2000 50 2 2 1 ../in/u.bin ../in/v.bin ../in/t.bin ../out/u/u ../out/v/v ../out/t/t 2
