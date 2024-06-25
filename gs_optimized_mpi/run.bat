mkdir -p out/t
mkdir -p out/u
mkdir -p out/v
mpiexec -n 8 build/main.exe 500.0 1.0 -1.0 1000 50 2 2 2 in/u.bin in/v.bin in/t.bin out/u/u out/v/v out/t/t
