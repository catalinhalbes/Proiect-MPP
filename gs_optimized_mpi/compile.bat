g++ --std=c++2a -O3 main.cpp -Wall -I"C:/Program Files (x86)/Microsoft SDKs/MPI/Include" "%MSMPI_LIB64%msmpi.lib" -march=native -mtune=native -flto -fuse-linker-plugin -o .\build\main
