#include <iostream>
#include <stdexcept>
#include <cstdint>
#include <string>
#include <chrono>
#include <mpi.h>

constexpr int U_TAG = 1;
constexpr int V_TAG = 2;
constexpr int T_TAG = 3;

double abs(double val) {
    return val < 0.0? -val: val;
}

struct indices {
    size_t x, y, z;
};

class Matrix3D {
    private:
        void do_cleanup() {
            if (elems == nullptr) {
                return;
            }

            dim_X = dim_Y = dim_Z = 0;
            delete[] elems;
            delete[] up_down_recv_buf;
            delete[] left_right_recv_buf;
            delete[] back_front_recv_buf;
            delete[] up_down_send_buf;
            delete[] left_right_send_buf;
            delete[] back_front_send_buf;
            elems = nullptr;
            up_down_recv_buf = nullptr;
            left_right_recv_buf = nullptr;
            back_front_recv_buf = nullptr;
            up_down_send_buf = nullptr;
            left_right_send_buf = nullptr;
            back_front_send_buf = nullptr;
        }

    public:
        double* elems;

        size_t original_dim_X, original_dim_Y, original_dim_Z;
        size_t original_stride_X, original_stride_Y;

        size_t actual_dim_X, actual_dim_Y, actual_dim_Z;

        size_t begin_X, begin_Y, begin_Z;
        size_t dim_X, dim_Y, dim_Z;
        size_t stride_X, stride_Y;
        size_t size;

        size_t rank;

        int left_cell, right_cell, back_cell, front_cell, down_cell, up_cell;
        MPI_Comm comm;

        size_t up_down_size, left_right_size, back_front_size;
        double *up_down_recv_buf, *left_right_recv_buf, *back_front_recv_buf;
        double *up_down_send_buf, *left_right_send_buf, *back_front_send_buf;
    
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

        Matrix3D(
            const std::string& filename, 
            const indices& processes, const indices& process_idx, 
            int left_cell, int right_cell, 
            int back_cell, int front_cell, 
            int down_cell, int up_cell, 
            MPI_Comm comm
        ) {
            this->left_cell  = left_cell;
            this->right_cell = right_cell;
            this->back_cell  = back_cell;
            this->front_cell = front_cell;
            this->down_cell  = down_cell;
            this->up_cell    = up_cell;
            this->comm       = comm;
            
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

            size_t dim_div_proc_x = (original_dim_X - 2) / processes.x;
            size_t dim_div_proc_y = (original_dim_Y - 2) / processes.y;
            size_t dim_div_proc_z = (original_dim_Z - 2) / processes.z;

            size_t dim_mod_proc_x = (original_dim_X - 2) % processes.x;
            size_t dim_mod_proc_y = (original_dim_Y - 2) % processes.y;
            size_t dim_mod_proc_z = (original_dim_Z - 2) % processes.z;

            // to find the size divide with the number of processes on the given axis and add 1 if the remainder is greater than the index of the current process
            // dimension includes the size of the shared borders
            actual_dim_X = dim_div_proc_x + (process_idx.x < dim_mod_proc_x); 
            actual_dim_Y = dim_div_proc_y + (process_idx.y < dim_mod_proc_y); 
            actual_dim_Z = dim_div_proc_z + (process_idx.z < dim_mod_proc_z); 

            dim_X = actual_dim_X + 2;
            dim_Y = actual_dim_Y + 2;
            dim_Z = actual_dim_Z + 2;

            // to find the begin offset fins the size of all block before
            // begin starts from the shared borders
            begin_X = dim_div_proc_x * process_idx.x + std::min(dim_mod_proc_x, process_idx.x);
            begin_Y = dim_div_proc_y * process_idx.y + std::min(dim_mod_proc_y, process_idx.y);
            begin_Z = dim_div_proc_z * process_idx.z + std::min(dim_mod_proc_z, process_idx.z);

            up_down_size    = actual_dim_X * actual_dim_Y;
            left_right_size = actual_dim_Y * actual_dim_Z;
            back_front_size = actual_dim_X * actual_dim_Z;

            up_down_recv_buf    = new double[up_down_size];
            left_right_recv_buf = new double[left_right_size];
            back_front_recv_buf = new double[back_front_size];

            up_down_send_buf    = new double[up_down_size];
            left_right_send_buf = new double[left_right_size];
            back_front_send_buf = new double[back_front_size];

            stride_X = dim_Y * dim_Z;
            stride_Y = dim_Z;
            size = dim_X * dim_Y * dim_Z;

            elems = new double[size];

            // printf("[%lu] beginX=%lu x=%lu beginY=%lu y=%lu beginZ=%lu z=%lu ud=%lu lr=%lu bf=%lu\n", rank, begin_X, dim_X, begin_Y, dim_Y, begin_Z, dim_Z, up_down_size, left_right_size, back_front_size);

            for (size_t i = 0; i < dim_X; i++) {
                for (size_t j = 0; j < dim_Y; j++) {
                    if (fseek(f, 8 * (3 + (i + begin_X) * original_stride_X + (j + begin_Y) * original_stride_Y + begin_Z), SEEK_SET) != 0) {
                        fclose(f);
                        do_cleanup();
                        throw std::runtime_error(
                            "[" + std::to_string(rank) + "] Could not fseek for read in '" + filename + 
                            "' loc (" + std::to_string(i + begin_X) + ", " + std::to_string(j + begin_Y) + ", " + std::to_string(begin_Z) + ")"
                        );
                    }

                    if (fread(&(elems[i * stride_X + j * stride_Y]), sizeof(double), dim_Z, f) != dim_Z) {
                        fclose(f);
                        do_cleanup();
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
            do_cleanup();
        }

        void write_in_file(const std::string& filename) {
            std::FILE* f = std::fopen(filename.c_str(), "rb+");

            if (f == nullptr) {
                throw std::runtime_error("Unable to open file: " + filename);
            }

            size_t x_stride_original_min2 = (original_dim_Y - 2) * (original_dim_Z - 2);
            size_t y_stride_original_min2 = (original_dim_Z - 2);
            size_t begin_x_min1 = begin_X - 1;
            size_t begin_y_min1 = begin_Y - 1;

            for (size_t i = 1; i < dim_X - 1; i++) {
                for (size_t j = 1; j < dim_Y - 1; j++) {

                    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    // WARNING: THIS CODE ASSUMES THAT ONLY THE MIDDLE PART OF THE MATRIX IS OUTPUTTED AND IGNORES THE WALLS
                    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    //if (fseek(f, 8 * (3 + (i + begin_X) * original_stride_X + (j + begin_Y) * original_stride_Y + begin_Z + 1), SEEK_SET) != 0) {
                    if (fseek(f, 8 * (3 + (i + begin_x_min1) * x_stride_original_min2 + (j + begin_y_min1) * y_stride_original_min2 + begin_Z), SEEK_SET) != 0) {
                        fclose(f);
                        throw std::runtime_error(
                            "[" + std::to_string(rank) + "] Could not fseek for write in '" + filename + 
                            "' loc (" + std::to_string(i + begin_X) + ", " + std::to_string(j + begin_Y) + ", " + std::to_string(begin_Z) + ")"
                        );
                    }

                    if (fwrite(&(elems[i * stride_X + j * stride_Y + 1]), sizeof(double), actual_dim_Z, f) != actual_dim_Z) {
                        fclose(f);
                        throw std::runtime_error(
                            "[" + std::to_string(rank) + "] Could not write all elements in '" + filename + 
                            "' loc (" + std::to_string(i + begin_X) + ", " + std::to_string(j + begin_Y) + ", " + std::to_string(begin_Z) + ")"
                        );
                    }
                }
            }

            // to keep: this code tries to output the outer walls of the full matrix without overlaps

            // // TODO: update the throw message to better reflect the actual position and state of the program

            // if (down_cell < 0) { // write z = 0 wall
            //     for (size_t i = 0 + (left_cell >= 0); i < dim_X - 1 - (right_cell >= 0); i++) {
            //         for (size_t j = 0 + (back_cell >= 0); j < dim_Y - 1 - (front_cell >= 0); j++) {
            //             if (fseek(f, 8 * (3 + (i + begin_X) * original_stride_X + (j + begin_Y) * original_stride_Y + begin_Z), SEEK_SET) != 0) {
            //                 fclose(f);
            //                 throw std::runtime_error(
            //                     "[" + std::to_string(rank) + "] Could not fseek for write in '" + filename + 
            //                     "' loc (" + std::to_string(i + begin_X) + ", " + std::to_string(j + begin_Y) + ", " + std::to_string(begin_Z) + ")"
            //                 );
            //             }

            //             if (fwrite(&(elems[i * stride_X + j * stride_Y]), sizeof(double), 1, f) != 1) {
            //                 fclose(f);
            //                 throw std::runtime_error(
            //                     "[" + std::to_string(rank) + "] Could not write all elements in '" + filename + 
            //                     "' loc (" + std::to_string(i + begin_X) + ", " + std::to_string(j + begin_Y) + ", " + std::to_string(begin_Z) + ")"
            //                 );
            //             }
            //         }
            //     }
            // }
            // if (up_cell < 0) { // write z = dim - 1 wall
            //     for (size_t i = 0 + (left_cell >= 0); i < dim_X - 1 - (right_cell >= 0); i++) {
            //         for (size_t j = 0 + (back_cell >= 0); j < dim_Y - 1 - (front_cell >= 0); j++) {
            //             if (fseek(f, 8 * (3 + (i + begin_X) * original_stride_X + (j + begin_Y) * original_stride_Y + begin_Z + dim_Z - 1), SEEK_SET) != 0) {
            //                 fclose(f);
            //                 throw std::runtime_error(
            //                     "[" + std::to_string(rank) + "] Could not fseek for write in '" + filename + 
            //                     "' loc (" + std::to_string(i + begin_X) + ", " + std::to_string(j + begin_Y) + ", " + std::to_string(begin_Z) + ")"
            //                 );
            //             }

            //             if (fwrite(&(elems[i * stride_X + j * stride_Y + dim_Z - 1]), sizeof(double), 1, f) != 1) {
            //                 fclose(f);
            //                 throw std::runtime_error(
            //                     "[" + std::to_string(rank) + "] Could not write all elements in '" + filename + 
            //                     "' loc (" + std::to_string(i + begin_X) + ", " + std::to_string(j + begin_Y) + ", " + std::to_string(begin_Z) + ")"
            //                 );
            //             }
            //         }
            //     }
            // }
            // if (left_cell < 0) { // write x = 0 wall
            //     int k = (int)(down_cell >= 0);
            //     int current_dim_z = dim_Z - k - (up_cell >= 0);
            //     for (size_t j = 0 + (back_cell >= 0); j < dim_Y - 1 - (front_cell >= 0); j++) {
            //         if (fseek(f, 8 * (3 + begin_X * original_stride_X + (j + begin_Y) * original_stride_Y + k + begin_Z), SEEK_SET) != 0) {
            //             fclose(f);
            //             throw std::runtime_error(
            //                 "[" + std::to_string(rank) + "] Could not fseek for write in '" + filename + 
            //                 "' loc (" + std::to_string(begin_X) + ", " + std::to_string(j + begin_Y) + ", " + std::to_string(begin_Z) + ")"
            //             );
            //         }

            //         if (fwrite(&(elems[j * stride_Y + k]), sizeof(double), current_dim_z, f) != (size_t)current_dim_z) {
            //             fclose(f);
            //             throw std::runtime_error(
            //                 "[" + std::to_string(rank) + "] Could not write all elements in '" + filename + 
            //                 "' loc (" + std::to_string(begin_X) + ", " + std::to_string(j + begin_Y) + ", " + std::to_string(begin_Z) + ")"
            //             );
            //         }
            //     }
            // }
            // if (right_cell < 0) { // write x = dim - 1 wall
            //     int k = (int)(down_cell >= 0);
            //     int current_dim_z = dim_Z - k - (up_cell >= 0);
            //     for (size_t j = 0 + (back_cell >= 0); j < dim_Y - 1 - (front_cell >= 0); j++) {
            //         if (fseek(f, 8 * (3 + (dim_X - 1 + begin_X) * original_stride_X + (j + begin_Y) * original_stride_Y + k + begin_Z), SEEK_SET) != 0) {
            //             fclose(f);
            //             throw std::runtime_error(
            //                 "[" + std::to_string(rank) + "] Could not fseek for write in '" + filename + 
            //                 "' loc (" + std::to_string(dim_X - 1 + begin_X) + ", " + std::to_string(j + begin_Y) + ", " + std::to_string(begin_Z) + ")"
            //             );
            //         }

            //         if (fwrite(&(elems[(dim_X - 1) * stride_X + j * stride_Y + k]), sizeof(double), current_dim_z, f) != (size_t)current_dim_z) {
            //             fclose(f);
            //             throw std::runtime_error(
            //                 "[" + std::to_string(rank) + "] Could not write all elements in '" + filename + 
            //                 "' loc (" + std::to_string(dim_X - 1 + begin_X) + ", " + std::to_string(j + begin_Y) + ", " + std::to_string(begin_Z) + ")"
            //             );
            //         }
            //     }
            // }
            // if (back_cell < 0) { // write y = 0 wall
            //     int k = (int)(down_cell >= 0);
            //     int current_dim_z = dim_Z - k - (up_cell >= 0);
            //     for (size_t i = 0 + (left_cell >= 0); i < dim_X - 1 - (right_cell >= 0); i++) {
            //         if (fseek(f, 8 * (3 + (i + begin_X) * original_stride_X + begin_Y * original_stride_Y + k + begin_Z), SEEK_SET) != 0) {
            //             fclose(f);
            //             throw std::runtime_error(
            //                 "[" + std::to_string(rank) + "] Could not fseek for write in '" + filename + 
            //                 "' loc (" + std::to_string(begin_X) + ", " + std::to_string(begin_Y) + ", " + std::to_string(begin_Z) + ")"
            //             );
            //         }

            //         if (fwrite(&(elems[i * stride_X + k]), sizeof(double), current_dim_z, f) != (size_t)current_dim_z) {
            //             fclose(f);
            //             throw std::runtime_error(
            //                 "[" + std::to_string(rank) + "] Could not write all elements in '" + filename + 
            //                 "' loc (" + std::to_string(begin_X) + ", " + std::to_string(begin_Y) + ", " + std::to_string(begin_Z) + ")"
            //             );
            //         }
            //     }
            // }
            // if (front_cell < 0) { // write y = dim - 1 wall
            //     int k = (int)(down_cell >= 0);
            //     int current_dim_z = dim_Z - k - (up_cell >= 0);
            //     for (size_t i = 0 + (left_cell >= 0); i < dim_X - 1 - (right_cell >= 0); i++) {
            //         if (fseek(f, 8 * (3 + (i + begin_X) * original_stride_X + (dim_Y - 1 + begin_Y) * original_stride_Y + k + begin_Z), SEEK_SET) != 0) {
            //             fclose(f);
            //             throw std::runtime_error(
            //                 "[" + std::to_string(rank) + "] Could not fseek for write in '" + filename + 
            //                 "' loc (" + std::to_string(begin_X) + ", " + std::to_string(dim_Y - 1 + begin_Y) + ", " + std::to_string(begin_Z) + ")"
            //             );
            //         }

            //         if (fwrite(&(elems[i * stride_X + (dim_Y - 1) * stride_Y + k]), sizeof(double), current_dim_z, f) != (size_t)current_dim_z) {
            //             fclose(f);
            //             throw std::runtime_error(
            //                 "[" + std::to_string(rank) + "] Could not write all elements in '" + filename + 
            //                 "' loc (" + std::to_string(begin_X) + ", " + std::to_string(begin_Y) + ", " + std::to_string(begin_Z) + ")"
            //             );
            //         }
            //     }
            // }
            std::fclose(f);
        }

        void create_matrix_file_preallocated(const std::string& filename) {
            if (rank != 0)
                fprintf(stderr, "File '%s' is being created by process [%lu] which might not be the intended one", filename.c_str(), rank);
            
            FILE* f = fopen(filename.c_str(), "wb");

            if (f == nullptr) {
                throw std::runtime_error("Unable to create file: " + filename);
            }

             // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            // WARNING: THIS CODE ASSUMES THAT ONLY THE MIDDLE PART OF THE MATRIX IS OUTPUTTED AND IGNORES THE WALLS
            //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            size_t ogx = original_dim_X - 2;
            size_t ogy = original_dim_Y - 2;
            size_t ogz = original_dim_Z - 2;
            size_t s = ogx * ogy * ogz * 8 + 3 * 8;
            // size_t s = original_dim_X * original_dim_Y * original_dim_Z * 8 + 3 * 8;

            std::fwrite(&ogx, sizeof(size_t), 1, f);
            std::fwrite(&ogy, sizeof(size_t), 1, f);
            std::fwrite(&ogz, sizeof(size_t), 1, f);

            // std::fwrite(&original_dim_X, sizeof(size_t), 1, f);
            // std::fwrite(&original_dim_Y, sizeof(size_t), 1, f);
            // std::fwrite(&original_dim_Z, sizeof(size_t), 1, f);

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

        void communicate(int tag) {
            // TODO: don't do copy if the process is at an edge
            // TODO: precompute the offsets that transform the actual coord to local total coords

            // z = 0 face (down)
            if (down_cell >= 0) {
                for (size_t i = 0; i < actual_dim_X; i++) {
                    for (size_t j = 0; j < actual_dim_Y; j++) {
                        up_down_send_buf[i * actual_dim_Y + j] = elems[(i + 1) * stride_X + (j + 1) * stride_Y + 1];
                    }
                }
            }

            MPI_Sendrecv(
                up_down_send_buf, up_down_size, MPI_DOUBLE, down_cell, tag, // send
                up_down_recv_buf, up_down_size, MPI_DOUBLE, up_cell, tag,   // receive
                comm, MPI_STATUS_IGNORE
            );

            if (up_cell >= 0) {
                for (size_t i = 0; i < actual_dim_X; i++) {
                    for (size_t j = 0; j < actual_dim_Y; j++) {
                        elems[(i + 1) * stride_X + (j + 1) * stride_Y + dim_Z - 1] = up_down_recv_buf[i * actual_dim_Y + j];
                    }
                }
            }

            // z = dim - 1 face (up)
            if (up_cell >= 0) {
                for (size_t i = 0; i < actual_dim_X; i++) {
                    for (size_t j = 0; j < actual_dim_Y; j++) {
                        up_down_send_buf[i * actual_dim_Y + j] = elems[(i + 1) * stride_X + (j + 1) * stride_Y + dim_Z - 2];
                    }
                }
            }

            MPI_Sendrecv(
                up_down_send_buf, up_down_size, MPI_DOUBLE, up_cell, tag,   // send
                up_down_recv_buf, up_down_size, MPI_DOUBLE, down_cell, tag, // receive
                comm, MPI_STATUS_IGNORE
            );

            if (down_cell >= 0) {
                for (size_t i = 0; i < actual_dim_X; i++) {
                    for (size_t j = 0; j < actual_dim_Y; j++) {
                        elems[(i + 1) * stride_X + (j + 1) * stride_Y] = up_down_recv_buf[i * actual_dim_Y + j];
                    }
                }
            }

            // x = 0 face (left)
            if (left_cell >= 0) {
                for (size_t j = 0; j < actual_dim_Y; j++) {
                    for (size_t k = 0; k < actual_dim_Z; k++) {
                        left_right_send_buf[j * actual_dim_Z + k] = elems[stride_X + (j + 1) * stride_Y + k + 1];
                    }
                }
            }

            MPI_Sendrecv(
                left_right_send_buf, left_right_size, MPI_DOUBLE, left_cell, tag,  // send
                left_right_recv_buf, left_right_size, MPI_DOUBLE, right_cell, tag, // receive
                comm, MPI_STATUS_IGNORE
            );

            if (right_cell >= 0) {
                for (size_t j = 0; j < actual_dim_Y; j++) {
                    for (size_t k = 0; k < actual_dim_Z; k++) {
                        elems[(dim_X - 1) * stride_X + (j + 1) * stride_Y + k + 1] = left_right_recv_buf[j * actual_dim_Z + k];
                    }
                }
            }

            // x = dim - 1 face (right)
            if (right_cell >= 0) {
                for (size_t j = 0; j < actual_dim_Y; j++) {
                    for (size_t k = 0; k < actual_dim_Z; k++) {
                        left_right_send_buf[j * actual_dim_Z + k] = elems[(dim_X - 2) * stride_X + (j + 1) * stride_Y + k + 1];
                    }
                }
            }

            MPI_Sendrecv(
                left_right_send_buf, left_right_size, MPI_DOUBLE, right_cell, tag, // send
                left_right_recv_buf, left_right_size, MPI_DOUBLE, left_cell, tag,  // receive
                comm, MPI_STATUS_IGNORE
            );

            if (left_cell >= 0) {
                for (size_t j = 0; j < actual_dim_Y; j++) {
                    for (size_t k = 0; k < actual_dim_Z; k++) {
                        elems[(j + 1) * stride_Y + k + 1] = left_right_recv_buf[j * actual_dim_Z + k];
                    }
                }
            }

            // y = 0 face (back)
            if (back_cell >= 0) {
                for (size_t i = 0; i < actual_dim_X; i++) {
                    for (size_t k = 0; k < actual_dim_Z; k++) {
                        back_front_send_buf[i * actual_dim_Z + k] = elems[(i + 1) * stride_X + stride_Y + k + 1];
                    }
                }
            }

            MPI_Sendrecv(
                back_front_send_buf, back_front_size, MPI_DOUBLE, back_cell, tag,  // send
                back_front_recv_buf, back_front_size, MPI_DOUBLE, front_cell, tag, // receive
                comm, MPI_STATUS_IGNORE
            );

            if (front_cell >= 0) {
                for (size_t i = 0; i < actual_dim_X; i++) {
                    for (size_t k = 0; k < actual_dim_Z; k++) {
                        elems[(i + 1) * stride_X + (dim_Y - 1) * stride_Y + k + 1] = back_front_recv_buf[i * actual_dim_Z + k];
                    }
                }
            }

            // y = dim - 1 face (front)
            if (front_cell >= 0) {
                for (size_t i = 0; i < actual_dim_X; i++) {
                    for (size_t k = 0; k < actual_dim_Z; k++) {
                        back_front_send_buf[i * actual_dim_Z + k] = elems[(i + 1) * stride_X + (dim_Y - 2) * stride_Y + k + 1];
                    }
                }
            }

            MPI_Sendrecv(
                back_front_send_buf, back_front_size, MPI_DOUBLE, front_cell, tag, // send
                back_front_recv_buf, back_front_size, MPI_DOUBLE, back_cell, tag,  // receive
                comm, MPI_STATUS_IGNORE
            );

            if (back_cell >= 0) {
                for (size_t i = 0; i < actual_dim_X; i++) {
                    for (size_t k = 0; k < actual_dim_Z; k++) {
                        elems[(i + 1) * stride_X + k + 1] = back_front_recv_buf[i * actual_dim_Z + k];
                    }
                }
            }
        }
};

inline double updateCells(Matrix3D& u, Matrix3D& v, Matrix3D& t, size_t idx, const double RA, const double DELTA) {
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

    double err_u = abs(u_mat[idx] - u_old);
    double err_v = abs(v_mat[idx] - v_old);
    double err_t = abs(t_mat[idx] - t_old);

    double err = std::max(err_u, err_v);
    err = std::max(err, err_t);

    return err;
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
    const size_t Cy = (size_t) Cy_aux;
    const size_t Cz = (size_t) Cz_aux;
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
	MPI_Comm mpi_comm_grid;
	MPI_Cart_create(MPI_COMM_WORLD, 3, dim_sizes, wrap_around, 0, &mpi_comm_grid);

	// get coords
	int coords[3];
	int grid_rank;
	MPI_Comm_rank(mpi_comm_grid, &grid_rank);
	MPI_Cart_coords(mpi_comm_grid, grid_rank, 3, coords);
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
	MPI_Cart_shift(mpi_comm_grid, 0, 1, &left_cell, &right_cell);
	MPI_Cart_shift(mpi_comm_grid, 1, 1, &back_cell, &front_cell);
	MPI_Cart_shift(mpi_comm_grid, 2, 1, &down_cell, &up_cell);

    // // 3D Cartesian coordinates
	// std::cout << "Rank:" << rank << "\n"
	// 	<< "Coords: (" << cx << ", " << cy << ", " << cz << ")" << std::endl;

    // // neighbors
	// std::cout << "Rank: " << rank << "\n"
	// 	<< "Left: " << left_cell << "; Right: " << right_cell << "\n"
	// 	<< "Back: " << back_cell << "; Front: " << front_cell << "\n"
	// 	<< "Down: " << down_cell << "; Up: " << up_cell << "\n"
	// 	<< std::endl;

    Matrix3D u(u_in, {Cx, Cy, Cz}, {cx, cy, cz}, left_cell, right_cell, back_cell, front_cell, down_cell, up_cell, mpi_comm_grid);
    Matrix3D v(v_in, {Cx, Cy, Cz}, {cx, cy, cz}, left_cell, right_cell, back_cell, front_cell, down_cell, up_cell, mpi_comm_grid);
    Matrix3D t(t_in, {Cx, Cy, Cz}, {cx, cy, cz}, left_cell, right_cell, back_cell, front_cell, down_cell, up_cell, mpi_comm_grid);

    Matrix3D u_new(u_in, {Cx, Cy, Cz}, {cx, cy, cz}, left_cell, right_cell, back_cell, front_cell, down_cell, up_cell, mpi_comm_grid);
    Matrix3D v_new(v_in, {Cx, Cy, Cz}, {cx, cy, cz}, left_cell, right_cell, back_cell, front_cell, down_cell, up_cell, mpi_comm_grid);
    Matrix3D t_new(t_in, {Cx, Cy, Cz}, {cx, cy, cz}, left_cell, right_cell, back_cell, front_cell, down_cell, up_cell, mpi_comm_grid);

    Matrix3D::check_size(u, v);
    Matrix3D::check_size(u, t);

    const size_t N1 = u.dim_X;
    const size_t N2 = u.dim_Y;
    const size_t N3 = u.dim_Z;

    const double DELTA = HEIGHT / (double) u.original_dim_Z;

    const size_t STRIDE_X = u.stride_X;
    const size_t STRIDE_Y = u.stride_Y;

    int64_t nr_it = 0;

    while (nr_it < MAX_IT) {
        double err= 0;

        nr_it += 1;

        if (DO_STEP && nr_it % STEP == 0) {
            const std::string t_file = t_out + "_iter_" + std::to_string(nr_it) + ".bin";
            const std::string u_file = u_out + "_iter_" + std::to_string(nr_it) + ".bin";
            const std::string v_file = v_out + "_iter_" + std::to_string(nr_it) + ".bin";

            if (rank == 0) {
                t.create_matrix_file_preallocated(t_file);
                u.create_matrix_file_preallocated(u_file);
                v.create_matrix_file_preallocated(v_file);
            }
            MPI_Barrier(mpi_comm_grid);
            t.write_in_file(t_file);
            u.write_in_file(u_file);
            v.write_in_file(v_file);
        }

        // communicate
        if (nr_it > 1) {
            u.communicate(U_TAG);
            v.communicate(V_TAG);
            t.communicate(T_TAG);
        }

        for (size_t i = 1; i < N1 - 1; i++) {
            const size_t i_idx = i * STRIDE_X;
            for (size_t j = 1; j < N2 - 1; j++) {
                const size_t j_idx = j * STRIDE_Y;
                for (size_t k = 1 + (i + t.begin_X + j + t.begin_Y + t.begin_Z) % 2 ; k < N3 - 1; k += 2) {
                    double new_err = updateCells(u, v, t, i_idx + j_idx + k, RA, DELTA);
                    err = std::max(err, new_err);
                }
            }
        }

        // communicate
        u.communicate(U_TAG);
        v.communicate(V_TAG);
        t.communicate(T_TAG);

        for (size_t i = 1; i < N1 - 1; i++) {
            const size_t i_idx = i * STRIDE_X;
            for (size_t j = 1; j < N2 - 1; j++) {
                const size_t j_idx = j * STRIDE_Y;
                for (size_t k = 1 + (i + t.begin_X + j + t.begin_Y + t.begin_Z + 1) % 2; k < N3 - 1; k += 2) {
                    double new_err = updateCells(u, v, t, i_idx + j_idx + k, RA, DELTA);
                    err = std::max(err, new_err);
                }
            }
        }

        // adiabatic conditions
        size_t yn1_idx = (N2 - 1) * STRIDE_Y;
        size_t yn2_idx = (N2 - 2) * STRIDE_Y;

        if (back_cell < 0) { // no neighbor y = 0
            for (size_t i = 0; i < N1; i++) {
                size_t x_idx = i * STRIDE_X;
                for (size_t k = 0; k < N3; k++) {
                    t.elems[x_idx + k] = t.elems[x_idx + STRIDE_Y + k];
                }
            }
        }

        if (front_cell < 0) { // no neighbor y = dim - 1
            for (size_t i = 0; i < N1; i++) {
                size_t x_idx = i * STRIDE_X;
                for (size_t k = 0; k < N3; k++) {
                    t.elems[x_idx + yn1_idx + k] = t.elems[x_idx + yn2_idx + k];
                }
            }
        }

        size_t zn1_idx = N3 - 1;
        size_t zn2_idx = N3 - 2;

        if (down_cell < 0) { // no neighbor
            for (size_t i = 0; i < N1; i++) {
                size_t x_idx = i * STRIDE_X;
                for (size_t j = 0; j < N2; j++) {
                    size_t j_idx = j * STRIDE_Y;
                    t.elems[x_idx + j_idx] = t.elems[x_idx + j_idx + 1];
                }
            }
        }

        if (up_cell < 0) { // no neighbor
            for (size_t i = 0; i < N1; i++) {
                size_t x_idx = i * STRIDE_X;
                for (size_t j = 0; j < N2; j++) {
                    size_t j_idx = j * STRIDE_Y;
                    t.elems[x_idx + j_idx + zn1_idx] = t.elems[x_idx + j_idx + zn2_idx];
                }
            }
        }

        if (USE_DIFF) {
            double max_err = 0;
            MPI_Reduce(&err, &max_err, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm_grid);

            int should_stop = 0;
            if (rank == 0 && max_err < MIN_DIFF) {
                should_stop = 1;
            }

            MPI_Bcast(&should_stop, 1, MPI_INT, 0, mpi_comm_grid);

            if (should_stop == 1) {
                break;
            }
        }
    }

    const std::string t_file = t_out + "_iter_" + std::to_string(nr_it) + ".bin";
    const std::string u_file = u_out + "_iter_" + std::to_string(nr_it) + ".bin";
    const std::string v_file = v_out + "_iter_" + std::to_string(nr_it) + ".bin";

    if (rank == 0) {
        t.create_matrix_file_preallocated(t_file);
        u.create_matrix_file_preallocated(u_file);
        v.create_matrix_file_preallocated(v_file);
    }
    MPI_Barrier(mpi_comm_grid);
    t.write_in_file(t_file);
    u.write_in_file(u_file);
    v.write_in_file(v_file);

    // find max
    double local_max = t.elems[0];
    for (size_t i = 0; i < t.size; i++) {
        local_max = std::max(local_max, t.elems[i]);
    }

    double global_max = local_max;
    MPI_Reduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm_grid);

    MPI_Finalize();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    if (rank == 0)
        printf("%.2f\n%.16f\n", duration.count() * 1000, global_max);
    return 0;
}