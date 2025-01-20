#include "mpi.h"
#include <iostream>

constexpr int Cx = 3;
constexpr int Cy = 3;
constexpr int Cz = 3;
constexpr int P = Cx * Cy * Cz;

int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv);

	int rank, p;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	//std::cout << "Rank: " << rank << "/" << p << std::endl;

	if (p != P) {
		if (rank == 0) {
			std::cerr << "MPI program was not executed with expected number of processes!" << "\n"
				<< "Expected: " << P << "\n"
				<< "Received: " << p << "\n"
				<< "Verify mpiexec parameters or source code constants (Cx, Cy, Cz)!" << std::endl;
		}
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// create grid
	int dim_sizes[3] = { Cx, Cy, Cz };
	int wrap_around[3] = { 0, 0 , 0 };
	MPI_Comm MPI_COMM_GRID;
	MPI_Cart_create(MPI_COMM_WORLD, 3, dim_sizes, wrap_around, 0, &MPI_COMM_GRID);

	// get coords
	int coords[3];
	int grid_rank;
	MPI_Comm_rank(MPI_COMM_GRID, &grid_rank);
	MPI_Cart_coords(MPI_COMM_GRID, grid_rank, 3, coords);
	int cx = coords[0];
	int cy = coords[1];
	int cz = coords[2];
	// X, Y, Z
	int left_cell, right_cell, back_cell, front_cell, down_cell, up_cell;
	MPI_Cart_shift(MPI_COMM_GRID, 0, 1, &left_cell, &right_cell);
	MPI_Cart_shift(MPI_COMM_GRID, 1, 1, &back_cell, &front_cell);
	MPI_Cart_shift(MPI_COMM_GRID, 2, 1, &down_cell, &up_cell);

	MPI_Finalize();
	return 0;
}