#include <iostream>
#include <utility>
#include <stdexcept>
#include <cstdint>
#include <string>

#include <chrono>
#include "omp.h"

constexpr double RA = 500.0;
constexpr double HEIGHT = 1.0;
constexpr double MIN_DIFF = 1e-15;
constexpr double COLD_WALL_TEMP = -1.0;
constexpr double HOT_WALL_TEMP = 1.0;

constexpr size_t MAX_IT = 10000;
// how often to output results in iterations
constexpr size_t STEP = 20;

constexpr size_t X_ASPECT_RATIO = 1;
constexpr size_t Y_ASPECT_RATIO = 1;
constexpr size_t Z_ASPECT_RATIO = 1;

constexpr size_t N = 100;
constexpr size_t N1 = N * X_ASPECT_RATIO;
constexpr size_t N2 = N * Y_ASPECT_RATIO;
constexpr size_t N3 = N * Z_ASPECT_RATIO;

constexpr double DELTA = HEIGHT / (double) N;

constexpr double N_THREADS = 16;

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

            this->elems = std::exchange(ot.elems, nullptr);
        }

        Matrix3D(const std::string& filename) {
            std::FILE* f = std::fopen(filename.c_str(), "rb");
            if (f == nullptr) {
                throw std::runtime_error("Unable to open file: " + filename);
            }

            std::fread(&dim_X, sizeof(size_t), 1, f);
            std::fread(&dim_Y, sizeof(size_t), 1, f);
            std::fread(&dim_Z, sizeof(size_t), 1, f);

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

            std::fwrite(&dim_X, sizeof(size_t), 1, f);
            std::fwrite(&dim_Y, sizeof(size_t), 1, f);
            std::fwrite(&dim_Z, sizeof(size_t), 1, f);

            if (std::fwrite(elems, sizeof(double), size, f) != size) {
                std::fclose(f);
                throw std::runtime_error("Unable to write all elements!");
            }
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
            this->elems = std::exchange(ot.elems, this->elems);
            this->dim_X = std::exchange(ot.dim_X, this->dim_X);
            this->dim_Y = std::exchange(ot.dim_Y, this->dim_Y);
            this->dim_Z = std::exchange(ot.dim_Z, this->dim_Z);
            this->stride_X = std::exchange(ot.stride_X, this->stride_X);
            this->stride_Y = std::exchange(ot.stride_Y, this->stride_Y);
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

errs updateCells(Matrix3D& u, Matrix3D& v, Matrix3D& t, size_t i, size_t j, size_t k) {
    double u_old = u.get(i, j, k);
    double v_old = v.get(i, j, k);
    double t_old = t.get(i, j, k);
    
    u.set(i, j, k,
        RA * DELTA / 12.0 * (t.get(i, j + 1, k) - t.get(i, j - 1, k)) + 
        (1.0 / 6.0) * (
            u.get(i - 1, j, k) + u.get(i + 1, j, k) + 
            u.get(i, j - 1, k) + u.get(i, j + 1, k) + 
            u.get(i, j, k - 1) + u.get(i, j, k + 1)
        )
    );

    v.set(i, j, k,
        -RA * DELTA / 12.0 * (t.get(i + 1, j, k) - t.get(i - 1, j, k)) + 
        (1.0 / 6.0) * (
            v.get(i - 1, j, k) + v.get(i + 1, j, k) + 
            v.get(i, j - 1, k) + v.get(i, j + 1, k) + 
            v.get(i, j, k - 1) + v.get(i, j, k + 1)
        )
    );

    t.set(i, j, k, 
        (1.0 / 6.0) * (
            t.get(i - 1, j, k) + t.get(i + 1, j, k) + 
            t.get(i, j - 1, k) + t.get(i, j + 1, k) + 
            t.get(i, j, k - 1) + t.get(i, j, k + 1)
            + (DELTA * DELTA) - (1.0 / 4.0) * (
                (u.get(i, j, k + 1) - u.get(i, j, k - 1)) * (t.get(i, j + 1, k) - t.get(i, j - 1, k)) -
                (u.get(i, j + 1, k) - u.get(i, j - 1, k)) * (t.get(i, j, k + 1) - t.get(i, j, k - 1)) +
                (v.get(i + 1, j, k) - v.get(i - 1, j, k)) * (t.get(i, j, k + 1) - t.get(i, j, k - 1)) -
                (v.get(i, j, k + 1) - v.get(i, j, k - 1)) * (t.get(i + 1, j, k) - t.get(i - 1, j, k))
            )
        )
    );

    double err_u = std::abs(u.get(i, j, k) - u_old);
    double err_v = std::abs(v.get(i, j, k) - v_old);
    double err_t = std::abs(t.get(i, j, k) - t_old);
    return errs {err_u, err_v, err_t};
}

int main(int argc, char* argv[]) {
    if (argc < 7) {
        printf("Usage: %s [u_mat_in] [v_mat_in] [t_mat_in] [u_out] [v_out] [t_out]\n", argv[0]);
        return 1;
    }

    auto t1 = std::chrono::high_resolution_clock::now();

    std::string u_in = argv[1];
    std::string v_in = argv[2];
    std::string t_in = argv[3];
    std::string u_out = argv[4];
    std::string v_out = argv[5];
    std::string t_out = argv[6];

    omp_set_num_threads(N_THREADS);

    Matrix3D u(N1, N2, N3, true);
    Matrix3D v(N1, N2, N3, true);
    Matrix3D t(N1, N2, N3, true);

    // x = 0 is the HOT wall
//#pragma omp parallel private(i, j, k) shared(local_A, local_B, local_C, block_len)
    #pragma omp parallel for collapse(2)
    for (size_t j = 0; j < N2; j += 1) {
        for (size_t k = 0; k < N3; k++) {
            t.set(0, j, k, HOT_WALL_TEMP);
        }
    } 

    // x = (N1-1) is the COLD wall
    #pragma omp parallel for collapse(2)
    for (size_t j = 0; j < N2; j += 1) {
        for (size_t k = 0; k < N3; k++) {
            t.set(N1 - 1, j, k, COLD_WALL_TEMP);
        }
    }

    size_t nr_it = 0;
    bool stop = false;

    double err_u = 0;
    double err_v = 0;
    double err_t = 0;

    while (/*!stop &&*/ nr_it < MAX_IT) {
        if (nr_it % STEP == 0) {
            t.write_to_file(t_out + ".iter_" + std::to_string(nr_it) + ".bin");
            u.write_to_file(u_out + ".iter_" + std::to_string(nr_it) + ".bin");
            v.write_to_file(v_out + ".iter_" + std::to_string(nr_it) + ".bin");
        }

        nr_it += 1;
    
        #pragma omp parallel for collapse(2)
        for (size_t i = 1; i < N1 - 1; i++) {
            for (size_t j = 1; j < N2 - 1; j++) {
                for (size_t k = 1 + (i + j) % 2; k < N3 - 1; k += 2) {
                    errs errors = updateCells(u, v, t, i, j, k);
                    err_u = std::max(err_u, errors.err_u);
                    err_v = std::max(err_v, errors.err_v);
                    err_t = std::max(err_t, errors.err_t);
                }
            }
        }

        #pragma omp parallel for collapse(2)
        for (size_t i = 1; i < N1 - 1; i++) {
            for (size_t j = 1; j < N2 - 1; j++) {
                for (size_t k = 1 + (i + j + 1) % 2; k < N3 - 1; k += 2) {
                    errs errors = updateCells(u, v, t, i, j, k);
                    err_u = std::max(err_u, errors.err_u);
                    err_v = std::max(err_v, errors.err_v);
                    err_t = std::max(err_t, errors.err_t);
                }
            }
        }

        // adiabatic conditions

        // y = 0 wall == y = 1 wall
        #pragma omp parallel for collapse(2)
        for (size_t i = 0; i < N1; i++) {
            for (size_t k = 0; k < N3; k++) {
                t.set(i, 0, k, t.get(i, 1, k));
            }
        }

        // y = (N2 - 1) wall == y = (N2 - 2) wall
        #pragma omp parallel for collapse(2)
        for (size_t i = 0; i < N1; i++) {
            for (size_t k = 0; k < N3; k++) {
                t.set(i, N2 - 1, k, t.get(i, N2 - 2, k));
            }
        }

        // z = 0 wall == z = 1 wall
        #pragma omp parallel for collapse(2)
        for (size_t i = 0; i < N1; i++) {
            for (size_t j = 0; j < N2; j++) {
                t.set(i, j, 0, t.get(i, j, 1));
            }
        }

        // z = (N3 - 1) wall == z = (N3 - 2) wall
        #pragma omp parallel for collapse(2)
        for (size_t i = 0; i < N1; i++) {
            for (size_t j = 0; j < N2; j++) {
                t.set(i, j, N3 - 1, t.get(i, j, N3 - 2));
            }
        }

        if (
            err_u < MIN_DIFF &&
            err_v < MIN_DIFF &&
            err_t < MIN_DIFF
        ) {
            stop = true;
        }
    }

    printf("Iterations: %lu\n", nr_it);

    u.write_to_file(u_out);
    v.write_to_file(v_out);
    t.write_to_file(t_out);

    auto t2 = std::chrono::high_resolution_clock::now();
    auto total = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(t2 - t1);
    std::cout << "Total time taken: " << total.count() << " ms" << std::endl;
}