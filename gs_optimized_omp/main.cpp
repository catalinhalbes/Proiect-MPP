#include <iostream>
#include <utility>
#include <stdexcept>
#include <cstdint>
#include <string>
#include <chrono>

#include "omp.h"

class Matrix3D {
    public:
        double* elems;
        size_t dim_X, dim_Y, dim_Z;
        size_t stride_X, stride_Y;
        size_t size;

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

        inline void check_bounds(size_t x, size_t y, size_t z) const {
            std::string err;

            if (x >= dim_X)
                err += "X index (" + std::to_string(x) + ") is out of range for dimension (" + std::to_string(dim_X) + ")\n";

            if (y >= dim_Y)
                err += "Y index (" + std::to_string(y) + ") is out of range for dimension (" + std::to_string(dim_Y) + ")\n";

            if (z >= dim_Z)
                err += "Z index (" + std::to_string(z) + ") is out of range for dimension (" + std::to_string(dim_Z) + ")\n";

            if (!err.empty()) {
                throw std::out_of_range(err);
            }
        }
    
    // public:
        Matrix3D(size_t dim_X, size_t dim_Y, size_t dim_Z, bool init_zero = false) {
            this->dim_X = dim_X;
            this->dim_Y = dim_Y;
            this->dim_Z = dim_Z;
            this->stride_X = dim_Y * dim_Z;
            this->stride_Y = dim_Z;
            this->size = dim_X * dim_Y * dim_Z;

            // printf("Allocating %llu elements\n", dim_X * dim_Y * dim_Z);

            if (init_zero) {
                elems = new double[size]();
            }
            else {
                elems = new double[size];
            }
        }

        Matrix3D(const Matrix3D& ot) {
            this->dim_X = ot.dim_X;
            this->dim_Y = ot.dim_Y;
            this->dim_Z = ot.dim_Z;
            this->stride_X = ot.stride_X;
            this->stride_Y = ot.stride_Y;
            this->size = ot.size;

            this->elems = new double[size];
            std::copy(ot.elems, ot.elems + ot.size, this->elems);
        }

        Matrix3D(Matrix3D&& ot) {
            this->dim_X = ot.dim_X;
            this->dim_Y = ot.dim_Y;
            this->dim_Z = ot.dim_Z;
            this->stride_X = ot.stride_X;
            this->stride_Y = ot.stride_Y;
            this->size = ot.size;

            this->elems = ot.elems;
            ot.elems = nullptr;
        }

        Matrix3D(const std::string& filename) {
            std::FILE* f = std::fopen(filename.c_str(), "rb");
            if (f == nullptr) {
                throw std::runtime_error("Unable to open file: " + filename);
            }

            std::fread(&dim_X, sizeof(size_t), 1, f);
            std::fread(&dim_Y, sizeof(size_t), 1, f);
            std::fread(&dim_Z, sizeof(size_t), 1, f);

            stride_X = dim_Y * dim_Z;
            stride_Y = dim_Z;
            size = dim_X * dim_Y * dim_Z;

            if (size == 0) {
                std::fclose(f);
                throw std::runtime_error("Invalid size for elems: (" + 
                    std::to_string(dim_X) + ", " + 
                    std::to_string(dim_Y) + ", " +
                    std::to_string(dim_Z) + ")"
                ); 
            }
            elems = new double[size];

            if (std::fread(elems, sizeof(double), size, f) != size) {
                std::fclose(f);
                delete[] elems;
                throw std::runtime_error("Unable to read all elements!");
            }
            std::fclose(f);
        }

        Matrix3D& operator= (const Matrix3D& ot) {
            // Guard self assignment
            if (this == &ot)
                return *this;
        
            if (this->size != ot.size) { // resource in *this cannot be reused
                this->size = ot.size;
                delete[]  elems;
                elems = new double[ot.size];
            } 

            this->dim_X = ot.dim_X;
            this->dim_Y = ot.dim_Y;
            this->dim_Z = ot.dim_Z;
            this->stride_X = ot.stride_X;
            this->stride_Y = ot.stride_Y;
        
            std::copy(ot.elems, ot.elems + ot.size, this->elems);
            return *this;
        }

        virtual ~Matrix3D() {
            if (elems == nullptr) {
                return;
            }

            dim_X = dim_Y = dim_Z = 0;
            delete[] elems;
            elems = nullptr;
        }

        void write_to_file(const std::string& filename) const {
            std::FILE* f = std::fopen(filename.c_str(), "wb");

            if (f == nullptr) {
                throw std::runtime_error("Unable to open file: " + filename);
            }

            size_t dx = dim_X - 2;
            size_t dy = dim_Y - 2;
            size_t dz = dim_Z - 2;

            std::fwrite(&dx, sizeof(size_t), 1, f);
            std::fwrite(&dy, sizeof(size_t), 1, f);
            std::fwrite(&dz, sizeof(size_t), 1, f);

            for (size_t i = 1; i < dim_X - 1; i++) {
                for (size_t j = 1; j < dim_Y - 1; j++) {
                    if (fseek(f, 8 * (3 + (i - 1) * dx * dz + (j - 1) * dz), SEEK_SET) != 0) {
                        fclose(f);
                        throw std::runtime_error(
                            "Unable seek in '" + filename + 
                            "' loc (" + std::to_string(i) + ", " + std::to_string(j) + ")"
                        );
                    }

                    if (fwrite(&(elems[i * stride_X + j * stride_Y + 1]), sizeof(double), dz, f) != dz) {
                        fclose(f);
                        throw std::runtime_error(
                            "Unable to write all element in '" + filename + 
                            "' loc (" + std::to_string(i) + ", " + std::to_string(j) + ")"
                        );
                    }
                }
            }

            // if (std::fwrite(elems, sizeof(double), size, f) != size) {
            //     std::fclose(f);
            //     throw std::runtime_error("Unable to write all elements!");
            // }
            std::fclose(f);
        }

        inline double get(size_t x, size_t y, size_t z) const {
            check_bounds(x, y, z);
            return elems[x * stride_X + y * stride_Y + z];
        }

        inline void set(size_t x, size_t y, size_t z, double val) {
            check_bounds(x, y, z);
            elems[x * stride_X + y * stride_Y + z] = val;
        }

        void swap(Matrix3D& ot) {
            check_size(*this, ot);
            double* aux;
            aux = elems;
            elems = ot.elems;
            ot.elems = aux;
        }

        double max() const {
            double max = elems[0];

            for (size_t i = 0; i < size; i++) {
                if (elems[i] > max) 
                    max = elems[i];
            }

            return max;
        }

        static double max_diff(const Matrix3D& mat1, const Matrix3D& mat2) {
            check_size(mat1, mat2);

            size_t s = mat1.dim_X * mat1.dim_Y * mat1.dim_Z;
            double max_diff = -1.0;

            for (size_t i = 0; i < s; i++) {
                double diff = abs(mat1.elems[i] - mat2.elems[i]);

                if (diff > max_diff) {
                    max_diff = diff;
                }
            }

            return max_diff;
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

    // printf("(%llu, %llu, %llu)\n", idx / x_stride, (idx - (idx / x_stride * x_stride)) / y_stride, (idx - (idx / x_stride * x_stride)) % y_stride);

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
    return errs {err_u, err_v, err_t};
}

int omp_thread_count() {
    int n = 0;
    #pragma omp parallel reduction(+:n)
    n += 1;
    return n;
}

int main(int argc, char* argv[]) {
    auto start = std::chrono::high_resolution_clock::now();

    if (argc != 13) {
        printf("Usage: %s <Ra> <height> <min_diff> <max_it> <save_every> <u_mat_in> <v_mat_in> <t_mat_in> <u_out> <v_out> <t_out> <num_threads>\n", argv[0]);
        printf("Ra          - the Raylenght constant\n");
        printf("height      - the z height\n");
        printf("min_diff    - minimum difference between iterations to continue the simulation (<= 0 disabled)\n");
        printf("max_it      - the maximum number of iterations\n");
        printf("save_every  - specify how often to output intermediary results (<= 0 disabled)\n");
        printf("num_threads - the number of threads to use\n");
        printf("\ninput matrices:  <u_mat_in> <v_mat_in> <t_mat_in>\n");
        printf("output matrices: <u_out> <v_out> <t_out>\n");
        exit(1);
    }

    const double RA       = std::stod(argv[1]);
    const double HEIGHT   = std::stod(argv[2]);
    const double MIN_DIFF = std::stod(argv[3]);
    const int64_t MAX_IT  = std::stoll(argv[4]);
    const int64_t STEP    = std::stoll(argv[5]);

    if (MAX_IT <= 0) {
        printf("Invalid maximum iterations! It should be strictly positive!\n");
        exit(1);
    }

    const std::string u_in    = argv[6];
    const std::string v_in    = argv[7];
    const std::string t_in    = argv[8];
    const std::string u_out   = argv[9];
    const std::string v_out   = argv[10];
    const std::string t_out   = argv[11];

    const int64_t NR_THR = std::stoll(argv[12]);

    if (NR_THR <= 0) {
        printf("Invalid number of threads! It should be strictly positive!\n");
        exit(1);
    }

    omp_set_num_threads(NR_THR);

    const bool USE_DIFF = MIN_DIFF > 0.0;
    const bool DO_STEP = STEP > 0;

    Matrix3D u(u_in);
    Matrix3D v(v_in);
    Matrix3D t(t_in);

    // printf("u: x: %llu y: %llu z: %llu size: %llu strideX: %llu strideY: %llu\n", u.dim_X, u.dim_Y, u.dim_Z, u.size, u.stride_X, u.stride_Y);
    // printf("v: x: %llu y: %llu z: %llu size: %llu strideX: %llu strideY: %llu\n", v.dim_X, v.dim_Y, v.dim_Z, v.size, v.stride_X, v.stride_Y);
    // printf("t: x: %llu y: %llu z: %llu size: %llu strideX: %llu strideY: %llu\n", t.dim_X, t.dim_Y, t.dim_Z, t.size, t.stride_X, t.stride_Y);

    Matrix3D::check_size(u, v);
    Matrix3D::check_size(u, t);

    const size_t N1 = u.dim_X;
    const size_t N2 = u.dim_Y;
    const size_t N3 = u.dim_Z;

    const double DELTA = HEIGHT / (double) N3;

    const size_t STRIDE_X = u.stride_X;
    const size_t STRIDE_Y = u.stride_Y;

    int64_t nr_it = 0;

    while (/*!stop &&*/ nr_it < MAX_IT) {
        double err_u = 0;
        double err_v = 0;
        double err_t = 0;

        nr_it += 1;

        if (DO_STEP && nr_it % STEP == 0) {
            t.write_to_file(t_out + "_iter_" + std::to_string(nr_it) + ".bin");
            u.write_to_file(u_out + "_iter_" + std::to_string(nr_it) + ".bin");
            v.write_to_file(v_out + "_iter_" + std::to_string(nr_it) + ".bin");
        }

        #pragma omp parallel for collapse(1) schedule(static, 8 * 32)
        for (size_t i = 1; i < N1 - 1; i++) {
            const size_t i_idx = i * STRIDE_X;
            for (size_t j = 1; j < N2 - 1; j++) {
                const size_t j_idx = j * STRIDE_Y;
                for (size_t k = 1 + (i + j) % 2; k < N3 - 1; k += 2) {
                    errs errors = updateCells(u, v, t, i_idx + j_idx + k, RA, DELTA);
                    err_u = std::max(err_u, errors.err_u);
                    err_v = std::max(err_v, errors.err_v);
                    err_t = std::max(err_t, errors.err_t);
                }
            }
        }

        #pragma omp parallel for collapse(2) schedule(static, 8 * 32)
        for (size_t i = 1; i < N1 - 1; i++) {
            const size_t i_idx = i * STRIDE_X;
            for (size_t j = 1; j < N2 - 1; j++) {
                const size_t j_idx = j * STRIDE_Y;
                for (size_t k = 1 + (i + j + 1) % 2; k < N3 - 1; k += 2) {
                    errs errors = updateCells(u, v, t, i_idx + j_idx + k, RA, DELTA);
                    err_u = std::max(err_u, errors.err_u);
                    err_v = std::max(err_v, errors.err_v);
                    err_t = std::max(err_t, errors.err_t);
                }
            }
        }

        // adiabatic conditions
        size_t yn1_idx = (N2 - 1) * STRIDE_Y;
        size_t yn2_idx = (N2 - 2) * STRIDE_Y;

        #pragma omp parallel for collapse(2) schedule(static, 8 * 32)
        for (size_t i = 0; i < N1; i++) {
            size_t x_idx = i * STRIDE_X;
            for (size_t k = 0; k < N3; k++) {
                t.elems[x_idx + k] = t.elems[x_idx + STRIDE_Y + k];
                t.elems[x_idx + yn1_idx + k] = t.elems[x_idx + yn2_idx + k];
            }
        }

        size_t zn1_idx = N3 - 1;
        size_t zn2_idx = N3 - 2;

        #pragma omp parallel for collapse(2) schedule(static, 8 * 32)
        for (size_t i = 0; i < N1; i++) {
            size_t x_idx = i * STRIDE_X;
            for (size_t j = 0; j < N2; j++) {
                size_t j_idx = j * STRIDE_Y;
                t.elems[x_idx + j_idx] = t.elems[x_idx + j_idx + 1];
                t.elems[x_idx + j_idx + zn1_idx] = t.elems[x_idx + j_idx + zn2_idx];
            }
        }

        if (USE_DIFF &&
            err_u < MIN_DIFF &&
            err_v < MIN_DIFF &&
            err_t < MIN_DIFF
        ) { break; }
    }

    t.write_to_file(t_out + "_iter_" + std::to_string(nr_it) + ".bin");
    u.write_to_file(u_out + "_iter_" + std::to_string(nr_it) + ".bin");
    v.write_to_file(v_out + "_iter_" + std::to_string(nr_it) + ".bin");

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    printf("%.2f\n", duration.count() * 1000);
}