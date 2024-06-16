#include <iostream>
#include <stdexcept>
#include <cstdint>
#include <string>

constexpr double RA = 100.0;
constexpr double HEIGHT = 1.0;
constexpr double MIN_DIFF = 1e-9;

constexpr size_t MAX_IT = 1000;

constexpr size_t X_ASPECT_RATIO = 1;
constexpr size_t Y_ASPECT_RATIO = 1;
constexpr size_t Z_ASPECT_RATIO = 1;

constexpr size_t N = 25;
constexpr size_t N1 = N * X_ASPECT_RATIO;
constexpr size_t N2 = N * Y_ASPECT_RATIO;
constexpr size_t N3 = N * Z_ASPECT_RATIO;

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
            size_t aux;
            double* aux_elems;

            aux_elems = elems;
            elems = ot.elems;
            ot.elems = aux_elems;

            aux = dim_X;
            dim_X = ot.dim_X;
            ot.dim_X = aux;

            aux = dim_Y;
            dim_Y = ot.dim_Y;
            ot.dim_Y = aux;

            aux = dim_Z;
            dim_Z = ot.dim_Z;
            ot.dim_Z = aux;

            aux = stride_X;
            stride_X = ot.stride_X;
            ot.stride_X = aux;

            aux = stride_Y;
            stride_Y = ot.stride_Y;
            ot.stride_Y = aux;
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

int main(int argc, char* argv[]) {
    if (argc < 7) {
        printf("Usage: %s [u_mat_in] [v_mat_in] [t_mat_in] [u_out] [v_out] [t_out]\n", argv[0]);
        return 1;
    }

    std::string u_in = argv[1];
    std::string v_in = argv[2];
    std::string t_in = argv[3];
    std::string u_out = argv[4];
    std::string v_out = argv[5];
    std::string t_out = argv[6];

    Matrix3D u(N1, N2, N3, true);
    Matrix3D v(N1, N2, N3, true);
    Matrix3D t(N1, N2, N3, true);

    Matrix3D u_new(N1, N2, N3);
    Matrix3D v_new(N1, N2, N3);
    Matrix3D t_new(N1, N2, N3);

    size_t nr_it = 0;
    bool stop = false;

    while (!stop && nr_it < MAX_IT) {
        nr_it += 1;

        for (size_t i = 1; i < N1 - 1; i++) {
            for (size_t j = 1; j < N2 - 1; j++) {
                for (size_t k = 1; k < N3 - 1; k++) {
                    u_new.set(i, j, k,
                        RA * HEIGHT / 12.0 * 
                        (t.get(i, j + 1, k) - t.get(i, j - 1, k)) + 
                        (1.0 / 6.0) * (
                            u.get(i - 1, j, k) + u.get(i + 1, j, k) + 
                            u.get(i, j - 1, k) + u.get(i, j + 1, k) + 
                            u.get(i, j, k - 1) + u.get(i, j, k + 1)
                        )
                    );

                    v_new.set(i, j, k,
                        RA * HEIGHT / 12.0 * 
                        (t.get(i + 1, j, k) - t.get(i - 1, j, k)) + 
                        (1.0 / 6.0) * (
                            v.get(i - 1, j, k) + v.get(i + 1, j, k) + 
                            v.get(i, j - 1, k) + v.get(i, j + 1, k) + 
                            v.get(i, j, k - 1) + v.get(i, j, k + 1)
                        )
                    );

                    t_new.set(i, j, k, 
                        (1.0 / 6.0) * (
                            t.get(i - 1, j, k) + t.get(i + 1, j, k) + 
                            t.get(i, j - 1, k) + t.get(i, j + 1, k) + 
                            t.get(i, j, k - 1) + t.get(i, j, k + 1)
                        ) + (HEIGHT * HEIGHT) - (1.0 / 4.0) * (
                            (u.get(i, j, k + 1) - u.get(i, j, k - 1)) * (t.get(i, j + 1, k) - t.get(i, j - 1, k)) -
                            (u.get(i, j + 1, k) - u.get(i, j - 1, k)) * (t.get(i, j, k + 1) - t.get(i, j, k - 1)) +
                            (v.get(i + 1, j, k) - v.get(i - 1, j, k)) * (t.get(i, j, k + 1) - t.get(i, j, k - 1)) -
                            (v.get(i, j, k + 1) - v.get(i, j, k - 1)) * (t.get(i + 1, j, k) - t.get(i - 1, j, k))
                        )
                    );

                    /**
                     * the example was modifying the values both at the end and inside the for loops
                     * changing the values inside the loops will make the iterations depend on the order of execution 
                     * this is propably unintended, so this overwriting will be omitted
                    */

                    // u.set(i, j, k, u_new.get(i, j, k));
                    // v.set(i, j, k, v_new.get(i, j, k));
                    // t.set(i, j, k, t_new.get(i, j, k));
                }
            }
        }

        // adiabatic condition
        for (size_t j = 1; j < N2 - 1; j++) {
            for (size_t k = 1; k < N3 - 1; k++) {
                t_new.set(0, j, k, t_new.get(1, j, k));
            }
        }

        for (size_t j = 1; j < N2 - 1; j++) {
            for (size_t k = 1; k < N3 - 1; k++) {
                t_new.set(N1 - 1, j, k, t_new.get(N1 - 2, j, k));
            }
        }

        for (size_t i = 1; i < N1 - 1; i++) {
            for (size_t k = 1; k < N3 - 1; k++) {
                t_new.set(i, 0, k, t_new.get(i, 1, k));
            }
        }

        for (size_t i = 1; i < N1 - 1; i++) {
            for (size_t k = 1; k < N3 - 1; k++) {
                t_new.set(i, N2 - 1, k, t_new.get(i, N2 - 2, k));
            }
        }

        if (
            Matrix3D::max_diff(u, u_new) < MIN_DIFF &&
            Matrix3D::max_diff(v, v_new) < MIN_DIFF &&
            Matrix3D::max_diff(t, t_new) < MIN_DIFF
        ) {
            stop = true;
        }

        // swap buffers to avoid reallocating (the _new matrices will be overriden anyways)
        u.swap(u_new);
        v.swap(v_new);
        t.swap(t_new);
    }

    printf("Iterations: %lu\n", nr_it);

    u.write_to_file(u_out);
    v.write_to_file(v_out);
    t.write_to_file(t_out);
}