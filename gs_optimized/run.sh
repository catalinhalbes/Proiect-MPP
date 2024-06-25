#!/bin/bash
./compile.sh
mkdir -p ./out/t
mkdir -p ./out/u
mkdir -p ./out/v
./build/main 1000 1.0 -1.0 2000 50 ./in/u.bin ./in/v.bin ./in/t.bin ./out/u/u ./out/v/v ./out/t/t
# ./build/main2 ./in/u.bin ./in/v.bin ./in/t.bin ./out/u/u ./out/v/v ./out/t/t
