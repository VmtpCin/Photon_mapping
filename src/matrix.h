#pragma once

#include <bits/stdc++.h>

template <typename T>
struct Matrix {
    std::vector<std::vector<T>> data;

    // Constructors
    Matrix() = default;
    Matrix(size_t m, size_t n) { data.assign(m, std::vector<T>(n)); }
    Matrix(const std::vector<std::vector<T>> &v) : data(v) {}

    // Access operators
    std::vector<T>& operator[](size_t idx)             { return data[idx]; }
    const std::vector<T>& operator[](size_t idx) const { return data[idx]; }

    // Get dimensions
    inline size_t rows() const { return data.size(); }
    inline size_t cols() const { return data.empty() ? 0 : data[0].size(); }

    Matrix<T> operator+(const Matrix<T> &other) const {
        if (rows() != other.rows() || cols() != other.cols()) {
            std::cerr << "Matrix dimensions mismatch for sum: A(" << rows()
                      << ", " << cols() << "); B(" << other.rows() << ", " << other.cols() << ")\n";
            return Matrix<T>(0, 0);
        }

        Matrix<T> result(rows(), cols());

        for (size_t i = 0; i < rows(); ++i)
            for (size_t j = 0; j < cols(); ++j)
                result[i][j] = (*this)[i][j] + other[i][j];

        return result;
    }

    Matrix<T>& operator+=(const Matrix<T> &other) {
        if (rows() != other.rows() || cols() != other.cols()) {
            std::cerr << "Matrix dimensions mismatch for sum: A(" << rows()
                      << ", " << cols() << "); B(" << other.rows() << ", " << other.cols() << ")\n";
            return *this;
        }

        for (size_t i = 0; i < rows(); ++i)
            for (size_t j = 0; j < cols(); ++j)
                (*this)[i][j] += other[i][j];

        return *this;
    }

    Matrix<T> operator-(const Matrix<T> &other) const {
        if (rows() != other.rows() || cols() != other.cols()) {
            std::cerr << "Matrix dimensions mismatch for minus: A(" << rows()
                      << ", " << cols() << "); B(" << other.rows() << ", " << other.cols() << ")\n";
            return Matrix<T>(0, 0);
        }

        Matrix<T> result(rows(), cols());

        for (size_t i = 0; i < rows(); ++i)
            for (size_t j = 0; j < cols(); ++j)
                result[i][j] = (*this)[i][j] - other[i][j];

        return result;
    }

    Matrix<T>& operator-=(const Matrix<T> &other) {
        if (rows() != other.rows() || cols() != other.cols()) {
            std::cerr << "Matrix dimensions mismatch for minus: A(" << rows()
                      << ", " << cols() << "); B(" << other.rows() << ", " << other.cols() << ")\n";
            return *this;
        }

        for (size_t i = 0; i < rows(); ++i)
            for (size_t j = 0; j < cols(); ++j)
                (*this)[i][j] -= other[i][j];

        return *this;
    }
    
    template<typename U>
    std::vector<U> operator*(const std::vector<U> &v) const {
        if (cols() != v.size()) {
            std::cerr << "Vector size mismatch for multiplication: A(" << rows()
                      << ", " << cols() << "); B(" << v.rows() << ", " << v.cols() << ")\n";
            return std::vector<U>(0);
        }

        std::vector<U> result(rows(), U{}); // Default-construct each element

        for (size_t i = 0; i < rows(); ++i)
            for (size_t j = 0; j < cols(); ++j)
                result[i] += data[i][j] * v[j];

        return result;
    }

    // Matrix-matrix multiplication
    template<typename U>
    Matrix<U> operator*(const Matrix<U> &other) const {
        if (cols() != other.rows()) {
            std::cerr << "Vector size mismatch for multiplication: A(" << rows()
                      << ", " << cols() << "); B(" << other.rows() << ", " << other.cols() << ")\n";
            return Matrix<U>(0, 0);
        }

        Matrix<U> result(rows(), other.cols());

        for (size_t i = 0; i < rows(); ++i)
            for (size_t j = 0; j < other.cols(); ++j)
                for (size_t k = 0; k < cols(); ++k)
                    result[i][j] += data[i][k] * other[k][j];

        return result;
    }

    // Get the cofactor of the matrix excluding the given row and column
    Matrix<T> getCofactor(int p, int q) const {
        Matrix<T> result(rows() - 1, cols() - 1);
        int i = 0, j = 0;

        for (int row = 0; row < rows(); ++row) {
            for (int col = 0; col < cols(); ++col) {
                if (row != p && col != q) {
                    result[i][j++] = data[row][col];
                    if (j == cols() - 1) {
                        j = 0;
                        ++i;
                    }
                }
            }
        }
        
        return result;
    }

    // Compute the determinant of the matrix
    T determinant() const {
        if (rows() != cols()) {
            std::cerr << "Determinant can only be computed for square matrices (" << rows() << ", " << cols() << ")\n";
            return 0;
        }

        if (rows() == 1)
            return data[0][0];

        T det = 0;
        int sign = 1;

        for (int f = 0; f < cols(); ++f) {
            det += sign * data[0][f] * getCofactor(0, f).determinant();
            sign = -sign;
        }

        return det;
    }

    // Compute the inverse of the matrix
    Matrix<T> inverse() const {
        T det = determinant();
        if (det == 0) {
            std::cerr << "Inverse cannot be found as the matrix is singular.\n";
            return Matrix<T>(0, 0);
        }

        if (rows() != cols()) {
            std::cerr << "Inverse cannot be found as the matrix is not square (" << rows() << ", " << cols() << ")\n";
            return Matrix<T>(0, 0);
        }

        Matrix<T> inv(rows(), cols());

        for (int i = 0; i < rows(); ++i) {
            for (int j = 0; j < cols(); ++j) {
                int sign = ((i + j) % 2 == 0) ? 1 : -1;
                inv[j][i] = sign * (getCofactor(i, j).determinant()) / det;
            }
        }

        return inv;
    }

    Matrix<T> transpose() const {
        Matrix<T> result(cols(), rows());

        for (size_t i = 0; i < rows(); ++i)
            for (size_t j = 0; j < cols(); ++j)
                result[j][i] = (*this)[i][j];

        return result;
    }

    T sumValuesAbsolute() const {
        T result = 0;

        for (size_t i = 0; i < rows(); ++i)
            for (size_t j = 0; j < cols(); ++j)
                result += std::abs((*this)[i][j]);

        return result;
    }

    // Print the matrix to the console
    void print() const {
        for (const auto &row : data) {
            for (const auto &val : row)
                printf("%6.2f ", val);
            printf("\n");
        }
        printf("\n");
    }
};
