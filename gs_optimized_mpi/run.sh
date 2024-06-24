#!/bin/bash
./compile.sh
mkdir -p ./out/t
mkdir -p ./out/u
mkdir -p ./out/v
mpirun -n 8 --hostfile hosts.txt ./build/main 500.0 1.0 -1.0 1000 50 2 2 2 ./in/test_u.bin ./in/test_v.bin ./in/test_t.bin ./out/u/u_test ./out/v/v_test ./out/t/t_test
