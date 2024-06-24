#!/bin/bash
mpic++ --std=c++11 -g -Wall -march=native -mtune=native -flto -fuse-linker-plugin -o ./build/main main.cpp
