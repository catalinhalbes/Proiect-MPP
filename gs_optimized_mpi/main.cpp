#include <iostream>
#include <utility>
#include <stdexcept>
#include <cstdint>
#include <string>
#include <chrono>
#include <mpi.h>

struct indices {
    size_t x, y, z;
};

class Matrix3D {
    public:
        double* elems;
        size_t original_dim_X, original_dim_Y, original_dim_Z;
        size_t original_stride_X, original_stride_Y;

        // TODO: Make sure that the values below are used when computing the iterations
        size_t begin_X, begin_Y, begin_Z;
        size_t dim_X, dim_Y, dim_Z;
        size_t stride_X, stride_Y;
        size_t size;

        size_t rank;
    
        static inline void check_size(const Matrix3D& mat1, const Matrix3D& mat2) {
            if (
                mat1.dim_X != mat2.dim_X ||
                mat1.dim_Y != mat2.dim_Y ||
                mat1.dim_Z != mat2.dim_Z
            ) {
                throw std::invalid_argument(
                    "Matrix of sizes (" + 
                    std::to_string(mat1.dim_X) + ", " + 
                    std::to_string(mat1.dim_Y) + ", " + 
                    std::to_string(mat1.dim_Z) + 
                    ") is incompatible with matrix of size (" +
                    std::to_string(mat2.dim_X) + ", " + 
                    std::to_string(mat2.dim_Y) + ", " + 
                    std::to_string(mat2.dim_Z) + ")"
                );
            }
        }

        Matrix3D(const std::string& filename, const indices& processes, const indices& process_idx) {
            std::FILE* f = std::fopen(filename.c_str(), "rb");
            if (f == nullptr) {
                throw std::runtime_error("Unable to open file: " + filename);
            }

            rank = process_idx.x * (processes.y + processes.z) + process_idx.y * processes.z + process_idx.z;

            std::fread(&original_dim_X, sizeof(size_t), 1, f);
            std::fread(&original_dim_Y, sizeof(size_t), 1, f);
            std::fread(&original_dim_Z, sizeof(size_t), 1, f);

            original_stride_X = original_dim_Y * original_dim_Z;
            original_stride_Y = original_dim_Z;

            size_t dim_div_proc_x = original_dim_X / processes.x;
            size_t dim_div_proc_y = original_dim_Y / processes.y;
            size_t dim_div_proc_z = original_dim_Z / processes.z;

            size_t dim_mod_proc_x = original_dim_X % processes.x;
            size_t dim_mod_proc_y = original_dim_Y % processes.y;
            size_t dim_mod_proc_z = original_dim_Z % processes.z;

            // to find the size divide with the number of processes on the given axis and add 1 if the remainder is greater than the index of the current process
            dim_X = dim_div_proc_x + (dim_mod_proc_x > process_idx.x); 
            dim_Y = dim_div_proc_y + (dim_mod_proc_y > process_idx.y); 
            dim_Z = dim_div_proc_z + (dim_mod_proc_z > process_idx.z); 

            stride_X = dim_Y * dim_Z;
            stride_Y = dim_Z;

            size = dim_X * dim_Y * dim_Z;

            elems = new double[size];

            // to find the begin offset fins the size of all block before
            begin_X = dim_div_proc_x * process_idx.x + (dim_mod_proc_x - 1 > process_idx.x? dim_mod_proc_x: process_idx.x);
            begin_Y = dim_div_proc_y * process_idx.y + (dim_mod_proc_y - 1 > process_idx.y? dim_mod_proc_y: process_idx.y);
            begin_Z = dim_div_proc_z * process_idx.z + (dim_mod_proc_z - 1 > process_idx.z? dim_mod_proc_z: process_idx.z);

            printf("[%lu] beginX=%lu x=%lu beginY=%lu y=%lu beginZ=%lu z=%lu\n", rank, begin_X, dim_X, begin_Y, dim_Y, begin_Z, dim_Z);

            for (size_t i = 0; i < dim_X; i++) {
                for (size_t j = 0; j < dim_Y; j++) {
                    if (fseek(f, 8 * (3 + (i + begin_X) * original_stride_X + (j + begin_Y) * original_stride_Y + begin_Z), SEEK_SET) != 0) {
                        fclose(f);
                        delete[] elems;
                        throw std::runtime_error(
                            "[" + std::to_string(rank) + "] Could not fseek for read in '" + filename + 
                            "' loc (" + std::to_string(i + begin_X) + ", " + std::to_string(j + begin_Y) + ", " + std::to_string(begin_Z) + ")"
                        );
                    }

                    if (fread(&(elems[i * stride_X + j * stride_Y]), sizeof(double), dim_Z, f) != dim_Z) {
                        fclose(f);
                        delete[] elems;
                        throw std::runtime_error(
                            "[" + std::to_string(rank) + "] Could not read all elements from '" + filename + 
                            "' loc (" + std::to_string(i + begin_X) + ", " + std::to_string(j + begin_Y) + ", " + std::to_string(begin_Z) + ")"
                        );
                    }
                }
            }

            std::fclose(f);
        }

        virtual ~Matrix3D() {
            if (elems == nullptr) {
                return;
            }

            dim_X = dim_Y = dim_Z = 0;
            delete[] elems;
            elems = nullptr;
        }

        void write_in_file(const std::string& filename) {
            std::FILE* f = std::fopen(filename.c_str(), "rb+");

            if (f == nullptr) {
                throw std::runtime_error("Unable to open file: " + filename);
            }

            for (size_t i = 0; i < dim_X; i++) {
                for (size_t j = 0; j < dim_Y; j++) {
                    if (fseek(f, 8 * (3 + (i + begin_X) * original_stride_X + (j + begin_Y) * original_stride_Y + begin_Z), SEEK_SET) != 0) {
                        fclose(f);
                        throw std::runtime_error(
                            "[" + std::to_string(rank) + "] Could not fseek for write in '" + filename + 
                            "' loc (" + std::to_string(i + begin_X) + ", " + std::to_string(j + begin_Y) + ", " + std::to_string(begin_Z) + ")"
                        );
                    }

                    if (fwrite(&(elems[i * stride_X + j * stride_Y]), sizeof(double), dim_Z, f) != dim_Z) {
                        fclose(f);
                        throw std::runtime_error(
                            "[" + std::to_string(rank) + "] Could not write all elements in '" + filename + 
                            "' loc (" + std::to_string(i + begin_X) + ", " + std::to_string(j + begin_Y) + ", " + std::to_string(begin_Z) + ")"
                        );
                    }
                }
            }
            std::fclose(f);
        }

        void create_matrix_file_preallocated(const std::string& filename) {
            if (rank != 0)
                fprintf(stderr, "File '%s' is being created by process [%lu] which might not be the intended one", filename.c_str(), rank);
            
            FILE* f = fopen(filename.c_str(), "wb");

            if (f == nullptr) {
                throw std::runtime_error("Unable to create file: " + filename);
            }

            size_t s = dim_X * dim_Y * dim_Z * 8 + 3 * 8;

            std::fwrite(&original_dim_X, sizeof(size_t), 1, f);
            std::fwrite(&original_dim_Y, sizeof(size_t), 1, f);
            std::fwrite(&original_dim_Z, sizeof(size_t), 1, f);

            if (fseek(f, s - 1, SEEK_SET) != 0) {
                fclose(f);
                throw std::runtime_error("Failed to seek in: " + filename);
            }

            if (fwrite("", 1, 1, f) != 1) {
                fclose(f);
                throw std::runtime_error("Failed to write at the end of: " + filename);
            }

            fclose(f);
        }
};

struct errs {
    double err_u;
    double err_v;
    double err_t;
};

inline errs updateCells(Matrix3D& u, Matrix3D& v, Matrix3D& t, size_t idx, const double RA, const double DELTA) {
    double* u_mat = u.elems;
    double* v_mat = v.elems;
    double* t_mat = t.elems;

    const double DELTA2 = DELTA * DELTA;

    // all strides should be equal
    const size_t x_stride = u.stride_X;
    const size_t y_stride = u.stride_Y;
    const size_t z_stride = 1;

    // precompute the neighbors
    const size_t up_idx     = idx + z_stride; // z + 1
    const size_t down_idx   = idx - z_stride; // z - 1
    const size_t right_idx  = idx + y_stride; // y + 1
    const size_t left_idx   = idx - y_stride; // y - 1
    const size_t back_idx   = idx + x_stride; // x + 1
    const size_t front_idx  = idx - x_stride; // x - 1

    double u_old = u_mat[idx];
    double v_old = v_mat[idx];
    double t_old = t_mat[idx];
    
    u_mat[idx] = 
        RA * DELTA / 12.0 * (t_mat[right_idx] - t_mat[left_idx]) + 
        (1.0 / 6.0) * (
            u_mat[front_idx] + u_mat[back_idx] + 
            u_mat[left_idx] + u_mat[right_idx] + 
            u_mat[down_idx] + u_mat[up_idx]
        );

    v_mat[idx] = 
        -RA * DELTA / 12.0 * (t_mat[back_idx] - t_mat[front_idx]) + 
        (1.0 / 6.0) * (
            v_mat[front_idx] + v_mat[back_idx] + 
            v_mat[left_idx] + v_mat[right_idx] + 
            v_mat[down_idx] + v_mat[up_idx]
        );

    t_mat[idx] =  
        (1.0 / 6.0) * (
            t_mat[front_idx] + t_mat[back_idx] + 
            t_mat[left_idx] + t_mat[right_idx] + 
            t_mat[down_idx] + t_mat[up_idx]
            + DELTA2 - (1.0 / 4.0) * (
                (u_mat[up_idx] - u_mat[down_idx]) * (t_mat[right_idx] - t_mat[left_idx]) -
                (u_mat[right_idx] - u_mat[left_idx]) * (t_mat[up_idx] - t_mat[down_idx]) +
                (v_mat[back_idx] - v_mat[front_idx]) * (t_mat[up_idx] - t_mat[down_idx]) -
                (v_mat[up_idx] - v_mat[down_idx]) * (t_mat[back_idx] - t_mat[front_idx])
            )
        );

    double err_u = std::abs(u_mat[idx] - u_old);
    double err_v = std::abs(v_mat[idx] - v_old);
    double err_t = std::abs(t_mat[idx] - t_old);
    return errs {err_u, err_v, err_t};
}

int main(int argc, char* argv[]) {
    auto start = std::chrono::high_resolution_clock::now();

    if (argc != 15) {
        printf("Usage: %s <Ra> <height> <min_diff> <max_it> <save_every> <proc_x> <proc_y> <proc_z> <u_mat_in> <v_mat_in> <t_mat_in> <u_out> <v_out> <t_out>\n", argv[0]);
        printf("Ra         - the Raylenght constant\n");
        printf("height     - the z height\n");
        printf("min_diff   - minimum difference between iterations to continue the simulation (<= 0 disabled)\n");
        printf("max_it     - the maximum number of iterations\n");
        printf("save_every - specify how often to output intermediary results (<= 0 disabled)\n");
        printf("proc_x     - the number of splits on the X axis\n");
        printf("proc_y     - the number of splits on the Y axis\n");
        printf("proc_z     - the number of splits on the Z axis\n");

        printf("\ninput matrices:  <u_mat_in> <v_mat_in> <t_mat_in>\n");
        printf("output matrices: <u_out> <v_out> <t_out>\n");
        exit(1);
    }

    const double RA       = std::stod(argv[1]);
    const double HEIGHT   = std::stod(argv[2]);
    const double MIN_DIFF = std::stod(argv[3]);
    const int64_t MAX_IT  = std::stoll(argv[4]);
    const int64_t STEP    = std::stoll(argv[5]);

    const int32_t Cx_aux  = std::stoi(argv[6]);
    const int32_t Cy_aux  = std::stoi(argv[7]);
    const int32_t Cz_aux  = std::stoi(argv[8]);

    if (MAX_IT <= 0) {
        fprintf(stderr, "Invalid maximum iterations! It should be strictly positive!\n");
        exit(1);
    }

    if (Cx_aux <= 0) {
        fprintf(stderr, "proc_x must be strictly positive\n");
        exit(1);
    }

    if (Cy_aux <= 0) {
        fprintf(stderr, "proc_y must be strictly positive\n");
        exit(1);
    }

    if (Cz_aux <= 0) {
        fprintf(stderr, "proc_z must be strictly positive\n");
        exit(1);
    }

    // avoid warnings
    const size_t Cx = (size_t) Cx_aux;
    const size_t Cy = (size_t) Cx_aux;
    const size_t Cz = (size_t) Cx_aux;
    const size_t P  = Cx * Cy * Cz;

    const std::string u_in    = argv[9];
    const std::string v_in    = argv[10];
    const std::string t_in    = argv[11];
    const std::string u_out   = argv[12];
    const std::string v_out   = argv[13];
    const std::string t_out   = argv[14];

    const bool USE_DIFF = MIN_DIFF > 0.0;
    const bool DO_STEP = STEP > 0;

    MPI_Init(&argc, &argv);

    int rank, p;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	if ((size_t)p != P) {
		if (rank == 0) {
            fprintf(stderr, "MPI program was executed a number of processes that does not match the given split!\n");
            fprintf(stderr, "MPI processes: %d\n", p);
            fprintf(stderr, "Split: (%ld, %ld, %ld) => %ld\n", Cx, Cy, Cz, P);
		}
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

    // create grid
	int dim_sizes[3] = { (int)Cx, (int)Cy, (int)Cz };
	int wrap_around[3] = { 0, 0, 0 };
	MPI_Comm MPI_COMM_GRID;
	MPI_Cart_create(MPI_COMM_WORLD, 3, dim_sizes, wrap_around, 0, &MPI_COMM_GRID);

	// get coords
	int coords[3];
	int grid_rank;
	MPI_Comm_rank(MPI_COMM_GRID, &grid_rank);
	MPI_Cart_coords(MPI_COMM_GRID, grid_rank, 3, coords);
	size_t cx = (size_t)coords[0];
	size_t cy = (size_t)coords[1];
	size_t cz = (size_t)coords[2];
	// X, Y, Z
    // left  = X - 1
    // right = X + 1
    // back  = Y - 1
    // front = Y + 1
    // down  = Z - 1
    // up    = Z + 1
	int left_cell, right_cell, back_cell, front_cell, down_cell, up_cell;
	MPI_Cart_shift(MPI_COMM_GRID, 0, 1, &left_cell, &right_cell);
	MPI_Cart_shift(MPI_COMM_GRID, 1, 1, &back_cell, &front_cell);
	MPI_Cart_shift(MPI_COMM_GRID, 2, 1, &down_cell, &up_cell);

    Matrix3D u(u_in, {Cx, Cy, Cz}, {cx, cy, cz});
    Matrix3D v(v_in, {Cx, Cy, Cz}, {cx, cy, cz});
    Matrix3D t(t_in, {Cx, Cy, Cz}, {cx, cy, cz});

    Matrix3D::check_size(u, v);
    Matrix3D::check_size(u, t);

    const size_t N1 = u.dim_X;
    const size_t N2 = u.dim_Y;
    const size_t N3 = u.dim_Z;

    const double DELTA = HEIGHT / (double) N3;

    const size_t STRIDE_X = u.stride_X;
    const size_t STRIDE_Y = u.stride_Y;

    int64_t nr_it = 0;

    // while (nr_it < MAX_IT) {
    //     double err_u = 0;
    //     double err_v = 0;
    //     double err_t = 0;

    //     if (DO_STEP && nr_it % STEP == 0) {
    //         const std::string t_file = t_out + "_iter_" + std::to_string(nr_it) + ".bin";
    //         const std::string u_file = u_out + "_iter_" + std::to_string(nr_it) + ".bin";
    //         const std::string v_file = v_out + "_iter_" + std::to_string(nr_it) + ".bin";

    //         if (rank == 0) {
    //             t.create_matrix_file_preallocated(t_file);
    //             u.create_matrix_file_preallocated(u_file);
    //             v.create_matrix_file_preallocated(v_file);
    //         }
    //         MPI_Barrier(MPI_COMM_GRID);
    //         t.write_in_file(t_file);
    //         u.write_in_file(u_file);
    //         v.write_in_file(v_file);
    //     }

    //     nr_it += 1;

    //     for (size_t i = 1; i < N1 - 1; i++) {
    //         const size_t i_idx = i * STRIDE_X;
    //         for (size_t j = 1; j < N2 - 1; j++) {
    //             const size_t j_idx = j * STRIDE_Y;
    //             for (size_t k = 1 + (i + j) % 2; k < N3 - 1; k += 2) {
    //                 errs errors = updateCells(u, v, t, i_idx + j_idx + k, RA, DELTA);
    //                 err_u = std::max(err_u, errors.err_u);
    //                 err_v = std::max(err_v, errors.err_v);
    //                 err_t = std::max(err_t, errors.err_t);
    //             }
    //         }
    //     }

    //     for (size_t i = 1; i < N1 - 1; i++) {
    //         const size_t i_idx = i * STRIDE_X;
    //         for (size_t j = 1; j < N2 - 1; j++) {
    //             const size_t j_idx = j * STRIDE_Y;
    //             for (size_t k = 1 + (i + j + 1) % 2; k < N3 - 1; k += 2) {
    //                 errs errors = updateCells(u, v, t, i_idx + j_idx + k, RA, DELTA);
    //                 err_u = std::max(err_u, errors.err_u);
    //                 err_v = std::max(err_v, errors.err_v);
    //                 err_t = std::max(err_t, errors.err_t);
    //             }
    //         }
    //     }

    //     // adiabatic conditions
    //     size_t yn1_idx = (N2 - 1) * STRIDE_Y;
    //     size_t yn2_idx = (N2 - 2) * STRIDE_Y;
    //     for (size_t i = 0; i < N1; i++) {
    //         size_t x_idx = i * STRIDE_X;
    //         for (size_t k = 0; k < N3; k++) {
    //             t.elems[x_idx + k] = t.elems[x_idx + STRIDE_Y + k];
    //             t.elems[x_idx + yn1_idx + k] = t.elems[x_idx + yn2_idx + k];
    //         }
    //     }

    //     size_t zn1_idx = N3 - 1;
    //     size_t zn2_idx = N3 - 2;
    //     for (size_t i = 0; i < N1; i++) {
    //         size_t x_idx = i * STRIDE_X;
    //         for (size_t j = 0; j < N2; j++) {
    //             size_t j_idx = j * STRIDE_Y;
    //             t.elems[x_idx + j_idx] = t.elems[x_idx + j_idx + 1];
    //             t.elems[x_idx + j_idx + zn1_idx] = t.elems[x_idx + j_idx + zn2_idx];
    //         }
    //     }

    //     if (USE_DIFF &&
    //         err_u < MIN_DIFF &&
    //         err_v < MIN_DIFF &&
    //         err_t < MIN_DIFF
    //     ) { break; }
    // }

    const std::string t_file = t_out + "_iter_" + std::to_string(nr_it) + ".bin";
    const std::string u_file = u_out + "_iter_" + std::to_string(nr_it) + ".bin";
    const std::string v_file = v_out + "_iter_" + std::to_string(nr_it) + ".bin";

    if (rank == 0) {
        t.create_matrix_file_preallocated(t_file);
        u.create_matrix_file_preallocated(u_file);
        v.create_matrix_file_preallocated(v_file);
    }
    MPI_Barrier(MPI_COMM_GRID);
    t.write_in_file(t_file);
    u.write_in_file(u_file);
    v.write_in_file(v_file);

    MPI_Finalize();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    printf("%.2f\n", duration.count() * 1000);
    return 0;
}