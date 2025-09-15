#pragma once

#include <cstddef>  // std::size_t
#include <vector>

class Matrix {
public:
    Matrix();
    Matrix(std::size_t nRow, std::size_t nCol, const double val = 0.0);
    Matrix(std::size_t nRow, std::size_t nCol, const std::vector<double> &data);

    Matrix(const Matrix& mat) noexcept;
    Matrix(Matrix&& mat) noexcept;
    
    ~Matrix();

    std::size_t rows() const noexcept;
    std::size_t cols() const noexcept;
    std::size_t size() const noexcept;

    double  operator()(std::size_t i, std::size_t j) const;
    double& operator()(std::size_t i, std::size_t j);
    double  operator()(std::size_t i) const;
    double& operator()(std::size_t i);
    
    double  at(std::size_t i, std::size_t j) const;
    double& at(std::size_t i, std::size_t j);
    double  at(std::size_t i) const;
    double& at(std::size_t i);

    const double* data() const;
    double* data();

    Matrix& operator=(const Matrix& mat) noexcept;
    Matrix& operator=(Matrix&& mat) noexcept;

    Matrix& operator+=(const double val);
    Matrix& operator-=(const double val);
    Matrix& operator*=(const double val);
    Matrix& operator/=(const double val);
    Matrix& operator+=(const Matrix &mat);
    Matrix& operator-=(const Matrix &mat);
    Matrix& operator*=(const Matrix &mat);
    Matrix& operator/=(const Matrix &mat);

    // set all elems the same value
    void fill(const double val);

    double det() const;
    double sum() const;

    // Frobenius norm, sqrt(\sum |a_{ij}|^2)
    double norm() const;

    // matrix manipulation
    Matrix trans() const;  // transpose
    Matrix inver() const;  // inverse, square matrix only

    Matrix diagg() const;  // diagonal matrix, N x 1
    double trace() const;  // trace, sum of diagonal elements

    Matrix reshape(std::size_t rRow, std::size_t rCol) const;

    Matrix sub_matrix(std::size_t rBeg, std::size_t rEnd,
                      std::size_t cBeg, std::size_t cEnd) const;

    void trans_self();
    void inver_self();
    void reshape_self(std::size_t rRow, std::size_t rCol);

    // Ax = b, assert A.nRow == A.nCol == b.nRow
    Matrix solve(const Matrix &b) const;
    
    friend Matrix operator+(const Matrix &mat);  // positive
    friend Matrix operator-(const Matrix &mat);  // negative
    friend Matrix operator+(const Matrix &mat, const double val);
    friend Matrix operator-(const Matrix &mat, const double val);
    friend Matrix operator*(const Matrix &mat, const double val);
    friend Matrix operator/(const Matrix &mat, const double val);
    friend Matrix operator+(const double val, const Matrix &mat);
    friend Matrix operator-(const double val, const Matrix &mat);
    friend Matrix operator*(const double val, const Matrix &mat);
    friend Matrix operator/(const double val, const Matrix &mat);
    friend Matrix operator+(const Matrix &m1, const Matrix &m2);
    friend Matrix operator-(const Matrix &m1, const Matrix &m2);
    friend Matrix operator*(const Matrix &m1, const Matrix &m2);
    friend Matrix operator/(const Matrix &m1, const Matrix &m2);

    // matrix product, assign: m1.col == m2.row
    friend Matrix operator%(const Matrix &m1, const Matrix &m2);

    // eigen solver
    friend class SelfAdjointEigenSolver;
    friend class GeneralSelfAdjointEigenSolver;

private:
    std::size_t nRow = 0;
    std::size_t nCol = 0;

    // use **mkl_malloc** to make sure memory is aligned
    double* _data = nullptr;
};

Matrix operator+(const Matrix &mat);
Matrix operator-(const Matrix &mat);
Matrix operator+(const Matrix &mat, const double val);
Matrix operator-(const Matrix &mat, const double val);
Matrix operator*(const Matrix &mat, const double val);
Matrix operator/(const Matrix &mat, const double val);
Matrix operator+(const double val, const Matrix &mat);
Matrix operator-(const double val, const Matrix &mat);
Matrix operator*(const double val, const Matrix &mat);
Matrix operator/(const double val, const Matrix &mat);
Matrix operator+(const Matrix &m1, const Matrix &m2);
Matrix operator-(const Matrix &m1, const Matrix &m2);
Matrix operator*(const Matrix &m1, const Matrix &m2);
Matrix operator/(const Matrix &m1, const Matrix &m2);

// matrix product, assign: m1.nCol == m2.nRow
Matrix operator%(const Matrix &m1, const Matrix &m2);


// HC = CE,
// where H is symmetric.
class SelfAdjointEigenSolver {
public:
    SelfAdjointEigenSolver(const Matrix &H);

    Matrix eigen_val() const { return eigenvalues; }
    Matrix eigen_vec() const { return eigenvectors; }
private:
    Matrix eigenvalues;
    Matrix eigenvectors;
};


// HC = SCE, 
// where H is symmetric, S is positive definite.
class GeneralSelfAdjointEigenSolver {
public:
    GeneralSelfAdjointEigenSolver(const Matrix &H, 
                                  const Matrix &S);

    Matrix eigen_val() const { return eigenvalues; }
    Matrix eigen_vec() const { return eigenvectors; }
private:
    Matrix eigenvalues;
    Matrix eigenvectors;
};