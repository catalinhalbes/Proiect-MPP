#!/bin/bash
./compile.sh
mpirun -n 8 --hostfile hosts.txt ./build/main 1000.0 1.0 -1.0 2000 50 2 2 2 ../in/u.bin ../in/v.bin ../in/t.bin ../out/u/u ../out/v/v ../out/t/t
