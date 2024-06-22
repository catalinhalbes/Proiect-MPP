#!/bin/bash
g++ --std=c++2a -O3 -Wall -march=native -mtune=native -flto -fuse-linker-plugin -o ./build/main main.cpp -fopenmp
