// TODO: size_t to MKL_INT      need checking
// TODO: size_t to lapack_int   need checking

#include "matrix_mkl.hpp"
#include <cstddef>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <cmath>
#include <utility>
#include <mkl.h>

Matrix::Matrix(): nRow(0), nCol(0), _data(nullptr) {}

Matrix::Matrix(std::size_t nRow, std::size_t nCol, const double val)
: nRow(nRow), nCol(nCol), _data(nullptr)
{
    std::size_t nSize = nRow * nCol;
    if (nSize == 0) return;

    _data = reinterpret_cast<double*>(mkl_malloc(nSize * sizeof(double), 64));

    if (_data == NULL) {
        std::cerr << "Matrix::Matrix(std::size_t nRow, std::size_t nCol, const double val)" << std::endl;
        std::cerr << "error: can't allocate memory." << std::endl;
        std::abort();
    }

    // TODO: performance improvement
    for (std::size_t i = 0; i < nSize; ++i) {
        _data[i] = val;
    }
}

Matrix::Matrix(std::size_t nRow, std::size_t nCol, const std::vector<double> &data)
: nRow(nRow), nCol(nCol), _data(nullptr)
{
    std::size_t nSize = nRow * nCol;

    if (data.size() != nSize) {
        std::cerr << "Matrix::Matrix(std::size_t nRow, std::size_t nCol, const std::vector<double> &data)" << std::endl;
        std::cerr << "error: data.size() != nRow * nCol" << std::endl;
        std::abort();
    }

    if (nSize == 0) return;

    _data = reinterpret_cast<double*>(mkl_malloc(nSize * sizeof(double), 64));

    if (_data == NULL) {
        std::cerr << "Matrix::Matrix(std::size_t nRow, std::size_t nCol, const std::vector<double> &data)" << std::endl;
        std::cerr << "error: can't allocate memory." << std::endl;
        std::abort();
    }

    // TODO: performance improvement
    std::memcpy(_data, data.data(), nSize * sizeof(double));
}

Matrix::Matrix(const Matrix& mat) noexcept
: nRow(mat.nRow), nCol(mat.nCol), _data(nullptr)
{
    std::size_t nSize = nRow * nCol;

    if (nSize == 0) return;

    _data = reinterpret_cast<double*>(mkl_malloc(nSize * sizeof(double), 64));

    if (_data == NULL) {
        std::cerr << "Matrix::Matrix(const Matrix& mat) noexcept" << std::endl;
        std::cerr << "error: can't allocate memory." << std::endl;
        std::abort();
    }

    // TODO: performance improvement
    std::memcpy(_data, mat._data, nSize * sizeof(double));
}

Matrix::Matrix(Matrix&& mat) noexcept
: nRow(mat.nRow), nCol(mat.nCol), _data(mat._data)
{
    mat.nRow = 0;
    mat.nCol = 0;
    mat._data = nullptr;
}

Matrix::~Matrix() {
    if (_data != nullptr) {
        mkl_free(_data);
        _data = nullptr;
    }
}

std::size_t Matrix::rows() const noexcept { return nRow; }

std::size_t Matrix::cols() const noexcept { return nCol; }

std::size_t Matrix::size() const noexcept { return nRow*nCol; }


double Matrix::operator()(std::size_t i, std::size_t j) const {
    assert(i < nRow);
    assert(j < nCol);
    return _data[i * nCol + j];
}

double& Matrix::operator()(std::size_t i, std::size_t j) {
    assert(i < nRow);
    assert(j < nCol);
    return _data[i * nCol + j];
}

double Matrix::operator()(std::size_t i) const {
    assert(i < nRow * nCol);
    return _data[i];
}

double& Matrix::operator()(std::size_t i) {
    assert(i < nRow * nCol);
    return _data[i];
}

double Matrix::at(std::size_t i, std::size_t j) const {
    if (i >= nRow || j >= nCol) {
        std::cerr << "double Matrix::at(std::size_t i, std::size_t j) const" << std::endl;
        std::cerr << "error: index out of range: i = " << i << ", j = " << j << std::endl;
        std::cerr << "where nRow = " << nRow << ", nCol = " << nCol << std::endl;
        std::abort();
    }
    return _data[i * nCol + j];
}

double& Matrix::at(std::size_t i, std::size_t j) {
    if (i >= nRow || j >= nCol) {
        std::cerr << "double& Matrix::at(std::size_t i, std::size_t j)" << std::endl;
        std::cerr << "error: index out of range: i = " << i << ", j = " << j << std::endl;
        std::cerr << "where nRow = " << nRow << ", nCol = " << nCol << std::endl;
        std::abort();
    }
    return _data[i * nCol + j];
}

double Matrix::at(std::size_t i) const {
    if (i >= nRow * nCol) {
        std::cerr << "Matrix::at(std::size_t i)" << std::endl;
        std::cerr << "error: index out of range: i = " << i << std::endl;
        std::cerr << "where size = " << nRow * nCol << std::endl;
        std::abort();
    }
    return _data[i];
}

double& Matrix::at(std::size_t i) {
    if (i >= nRow * nCol) {
        std::cerr << "MatrixZ::at(std::size_t i)" << std::endl;
        std::cerr << "error: index out of range: i = " << i << std::endl;
        std::cerr << "where size = " << nRow * nCol << std::endl;
        std::abort();
    }
    return _data[i];
}

const double* Matrix::data() const {
    return _data;
}

double* Matrix::data() {
    return _data;
}

Matrix& Matrix::operator=(const Matrix& mat) noexcept {
    if (this == &mat) return *this;

    // TODO: performance improvement if data is the size
    if (_data != nullptr) {
        mkl_free(_data);
        _data = nullptr;
    }

    nRow = mat.nRow;
    nCol = mat.nCol;

    std::size_t nSize = nRow * nCol;
    if (nSize == 0) return *this;

    _data = reinterpret_cast<double*>(mkl_malloc(nSize * sizeof(double), 64));

    if (_data == NULL) {
        std::cerr << "Matrix& Matrix::operator=(const Matrix& mat)" << std::endl;
        std::cerr << "error: can't allocate memory." << std::endl;
        std::abort();
    }

    // TODO: performance improvement
    std::memcpy(_data, mat._data, nSize * sizeof(double));

    return *this;
}

Matrix& Matrix::operator=(Matrix&& mat) noexcept {
    if (this == &mat) return *this;
    
    if (_data != nullptr) {
        mkl_free(_data);
        _data = nullptr;
    }

    nRow = mat.nRow;
    nCol = mat.nCol;
    _data = mat._data;

    mat.nRow = 0;
    mat.nCol = 0;
    mat._data = nullptr;


    // following implementation is not good,
    // because it not make sure mat is empty after move.

    // std::swap(nRow, mat.nRow);
    // std::swap(nCol, mat.nCol);
    // std::swap(_data, mat._data);

    return *this;
}


Matrix& Matrix::operator+=(const double val) {
    Matrix tmp(rows(), cols(), val);
    (*this) += tmp;
    return *this;
}

Matrix& Matrix::operator-=(const double val) {
    Matrix tmp(rows(), cols(), val);
    (*this) -= tmp;
    return *this;
}

Matrix& Matrix::operator*=(const double val) {
    cblas_dscal((const MKL_INT) nRow * nCol,
                (const double) val,
                (double *) _data,
                (const MKL_INT) 1);

    return *this;
}

Matrix& Matrix::operator/=(const double val) {
    const double valInv = 1.0 / val;
    (*this) *= valInv;
    return *this;
}


Matrix& Matrix::operator+=(const Matrix &mat) {
    assert(mat.nRow == nRow);
    assert(mat.nCol == nCol);

    if (mat.nRow != nRow) {
        std::cerr << "Matrix& Matrix::operator+=(const Matrix &mat)" << std::endl;
        std::cerr << "error: mat.nRow != nRow." << std::endl;
        std::abort();
    }

    if (mat.nCol != nCol) {
        std::cerr << "Matrix& Matrix::operator+=(const Matrix &mat)" << std::endl;
        std::cerr << "error: mat.nCol != nCol." << std::endl;
        std::abort();
    }

    vdAdd((const MKL_INT) nRow*nCol, 
          (const double*) _data, 
          (const double*) mat._data, 
          (double*) _data);

    return *this;
}

Matrix& Matrix::operator-=(const Matrix &mat) {
    assert(mat.nRow == nRow);
    assert(mat.nCol == nCol);

    if (mat.nRow != nRow) {
        std::cerr << "Matrix& Matrix::operator-=(const Matrix &mat)" << std::endl;
        std::cerr << "error: mat.nRow != nRow." << std::endl;
        std::abort();
    }

    if (mat.nCol != nCol) {
        std::cerr << "Matrix& Matrix::operator-=(const Matrix &mat)" << std::endl;
        std::cerr << "error: mat.nCol != nCol." << std::endl;
        std::abort();
    }

    vdSub((const MKL_INT) nRow*nCol, 
          (const double*) _data, 
          (const double*) mat._data, 
          (double*) _data);

    return *this;
}

Matrix& Matrix::operator*=(const Matrix &mat) {
    assert(mat.nRow == nRow);
    assert(mat.nCol == nCol);

    if (mat.nRow != nRow) {
        std::cerr << "Matrix& Matrix::operator*=(const Matrix &mat)" << std::endl;
        std::cerr << "error: mat.nRow != nRow." << std::endl;
        std::abort();
    }

    if (mat.nCol != nCol) {
        std::cerr << "Matrix& Matrix::operator*=(const Matrix &mat)" << std::endl;
        std::cerr << "error: mat.nCol != nCol." << std::endl;
        std::abort();
    }

    vdMul((const MKL_INT) nRow*nCol, 
          (const double*) _data, 
          (const double*) mat._data, 
          (double*) _data);

    return *this;
}

Matrix& Matrix::operator/=(const Matrix &mat) {
    assert(mat.nRow == nRow);
    assert(mat.nCol == nCol);

    if (mat.nRow != nRow) {
        std::cerr << "Matrix& Matrix::operator/=(const Matrix &mat)" << std::endl;
        std::cerr << "error: mat.nRow != nRow." << std::endl;
        std::abort();
    }

    if (mat.nCol != nCol) {
        std::cerr << "Matrix& Matrix::operator/=(const Matrix &mat)" << std::endl;
        std::cerr << "error: mat.nCol != nCol." << std::endl;
        std::abort();
    }

    vdDiv((const MKL_INT) nRow*nCol, 
          (const double*) _data, 
          (const double*) mat._data, 
          (double*) _data);

    return *this;
}

void Matrix::fill(const double val) {
    std::size_t nSize = nRow * nCol;
    
    // TODO: performance improvement
    for (std::size_t i = 0; i < nSize; ++i) {
        _data[i] = val;
    }
}

double Matrix::det() const {
    // TODO: some check
    assert(nRow == nCol);

    if (nRow*nCol == 0) {
        return 0.0;
    }

    std::size_t N = nRow;
    int* ipiv = (int *) std::malloc(N * sizeof(int));
    int info = -1;

    // 注意：MKL 的 dgetrf 认为输入是列优先的！
    // 所以我们要告知它我们的数据其实是转置后的（列主序视角）

    // 调用 MKL 的 LU 分解（列主序矩阵 A 的 LU 分解）

    Matrix tmp(*this);

    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 
                          (lapack_int) N, 
                          (lapack_int) N, 
                          (double *) tmp._data, 
                          (lapack_int) N, 
                          (lapack_int *) ipiv);

    if (info != 0) {
        std::cerr << "LU decomposition failed: info = " << info << std::endl;
        delete[] ipiv;
        return 0.0;
    }

    double det = 1.0;
    for (std::size_t i = 0; i < N; ++i) {
        det *= tmp(i, i);  // 行优先对角元素位置
    }

    // 计算置换的符号（奇数次交换为负）
    int sign = 1;
    for (std::size_t i = 0; i < N; ++i) {
        if (ipiv[i] != i + 1) {  // ipiv 从 1 开始计数
            sign = -sign;        // 每次交换改变符号
        }
    }
    det *= sign;

    free(ipiv);
    return det;
}

double Matrix::sum() const {
    double ret = 0.0;
    std::size_t nSize = nRow * nCol;

    // TODO: performance improvement
    for (std::size_t i = 0; i < nSize; ++i) {
        ret += _data[i];
    }

    return ret;
}

double Matrix::norm() const {
    double ret = 0.0;
    std::size_t nSize = nRow * nCol;

    // TODO: performance improvement
    for (std::size_t i = 0; i < nSize; ++i) {
        ret += _data[i] * _data[i];
    }

    return std::sqrt(ret);
}

Matrix Matrix::trans() const {
    Matrix ret;
    ret.nRow = nCol;
    ret.nCol = nRow;

    std::size_t nSize = ret.nRow * ret.nCol;
    if (nSize == 0) return ret;

    ret._data = reinterpret_cast<double*>(mkl_malloc(nSize * sizeof(double), 64));

    if (ret._data == NULL) {
        std::cerr << "Matrix Matrix::trans() const" << std::endl;
        std::cerr << "error: can't allocate memory." << std::endl;
        std::abort();
    }

    const double alpha = 1.0;
    
    // B := alpha*op(A)
    mkl_domatcopy('R',      // row major
                  'T',      // conjugate transpose
                  nRow,     // The number of rows in matrix A (the input matrix).
                  nCol,     // The number of columns in matrix A (the input matrix).
                  alpha,    // B := alpha*op(A)
                  _data, nCol, 
                  ret._data, nRow);

    return ret;
}


Matrix Matrix::inver() const {
    if (nRow != nCol) {
        std::cerr << "Matrix Matrix::inver() const" << std::endl;
        std::cerr << "error: matrix is not square." << std::endl;
        std::abort();
    }

    std::size_t N = nRow;
    Matrix ret(*this);

    if (nRow * nCol == 0) return ret;

    lapack_int info = -1;
    lapack_int *ipiv = reinterpret_cast<lapack_int*>(std::malloc(N * sizeof(lapack_int)));

    // step 1: Computes the LU factorization
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 
                          (lapack_int) N, 
                          (lapack_int) N, 
                          (double *) ret._data, 
                          (lapack_int) N, 
                          (lapack_int *) ipiv);

    if (info != 0) {
        std::cerr << "Matrix Matrix::inver() const" << std::endl;
        std::cerr << "error: LAPACKE_dgetrf failed with info = " << info << std::endl;
        std::abort();
    }

    info = -1;

    // step 2: Computes the inverse of an LU-factored general matrix
    info = LAPACKE_dgetri(LAPACK_ROW_MAJOR,
                          (lapack_int) N, 
                          (double *) ret._data, 
                          (lapack_int) N,
                          (const lapack_int *) ipiv);

    if (info != 0) {
        std::cerr << "Matrix Matrix::inver() const" << std::endl;
        std::cerr << "error: LAPACKE_dgetri failed with info = " << info << std::endl;
        std::abort();
    }

    std::free(ipiv);

    return ret;
}

// diagonal matrix, N x 1
Matrix Matrix::diagg() const {
    std::size_t N = std::min(nRow, nCol);

    Matrix ret;
    ret.nRow = N;
    ret.nCol = 1;

    // diag of a empty matrix
    if (N == 0) return ret;

    ret._data = reinterpret_cast<double*>(mkl_malloc(N * sizeof(double), 64));

    if (ret._data == NULL) {
        std::cerr << "Matrix Matrix::diagg() const" << std::endl;
        std::cerr << "error: can't allocate memory." << std::endl;
        std::abort();
    }

    // void cblas_dcopy(const MKL_INT n, const double *x, const MKL_INT incx, 
    //                                         double *y, const MKL_INT incy);
    
    cblas_dcopy(
        (const MKL_INT) N,
        (const double *) _data,
        (const MKL_INT) (nCol + 1),
        (double *) ret._data,
        (const MKL_INT) 1
    );

    return ret;
}


// double trace() const;  // trace, sum of diagonal elements


Matrix Matrix::reshape(std::size_t rRow, std::size_t rCol) const {
    if (rRow * rCol != nRow * nCol) {
        std::cerr << "Matrix Matrix::reshape(std::size_t rRow, std::size_t rCol) const" << std::endl;
        std::cerr << "error: reshape size not match." << std::endl;
        std::abort();
    }

    Matrix ret(*this);
    ret.nRow = rRow;
    ret.nCol = rCol;

    return ret;
}


Matrix Matrix::solve(const Matrix &b) const {
    if (nRow != nCol) {
        std::cerr << "Matrix Matrix::solve(const Matrix &b) const" << std::endl;
        std::cerr << "error: Ax = b, where A is not square." << std::endl;
        std::abort();
    }

    if (nRow != b.nRow) {
        std::cerr << "Matrix Matrix::solve(const Matrix &b) const" << std::endl;
        std::cerr << "error: Ax = b, where A.nRow != b.nRow." << std::endl;
        std::abort();
    }

    if (b.size() == 0) return b; 

    // Ax = B
    Matrix Acp(*this);
    Matrix Bcp(b);

    lapack_int info = -1;
    lapack_int *ipiv = reinterpret_cast<lapack_int*>(std::malloc(nRow * sizeof(lapack_int)));

    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR,
                         (lapack_int) nRow,
                         (lapack_int) Bcp.nCol,
                         (double *) Acp._data, (lapack_int) nCol,
                         (lapack_int *) ipiv,
                         (double *) Bcp._data, (lapack_int) Bcp.nCol);
    
    if (info != 0) {
        std::cerr << "Matrix Matrix::solve(const Matrix &b) const" << std::endl;
        std::cerr << "error: LAPACKE_zgesv failed with info = " << info << std::endl;
        std::abort();
    }

    std::free(ipiv);

    return Bcp;
}


/******************************************************************************/
/*                          positive  and  negative                           */
/******************************************************************************/

Matrix operator+(const Matrix &mat) {
    Matrix ret(mat);
    return ret;
}

Matrix operator-(const Matrix &mat) {
    Matrix ret(mat);

    // x = a*x
    cblas_dscal((const MKL_INT) ret.size(),
                (const double) -1.0,
                (double *) ret._data,
                (const MKL_INT) 1);

    return ret;
}


/******************************************************************************/
/*                           Matrix  and  Value                               */
/******************************************************************************/

// implementation by Matrix and Matrix
Matrix operator+(const Matrix &mat, const double val) {
    Matrix ret(mat.nRow, mat.nCol, val);
    ret += mat;
    return ret;
}

Matrix operator-(const Matrix &mat, const double val) {
    Matrix ret(mat.nRow, mat.nCol, -val);
    ret += mat;  // -val + mat
    return ret;
}

// implementation by Matrix and Value
Matrix operator*(const Matrix &mat, const double val) {
    Matrix ret(mat);
    ret *= val;
    return ret;
}

Matrix operator/(const Matrix &mat, const double val) {
    const double valInv = 1.0 / val;

    Matrix ret(mat);
    ret *= valInv;
    return ret;
}


/******************************************************************************/
/*                           Value  and  Matrix                               */
/******************************************************************************/

// implementation by Matrix and Matrix
Matrix operator+(const double val, const Matrix &mat) {
    Matrix ret(mat.nRow, mat.nCol, val);
    ret += mat;
    return ret;
}

Matrix operator-(const double val, const Matrix &mat) {
    Matrix ret(mat.nRow, mat.nCol, val);
    ret -= mat;
    return ret;
}

// implementation by Matrix and Value
Matrix operator*(const double val, const Matrix &mat) {
    Matrix ret(mat);
    ret *= val;
    return ret;
}

Matrix operator/(const double val, const Matrix &mat) {
    Matrix ret(mat.nRow, mat.nCol, val);
    ret /= mat;
    return ret;
}


/******************************************************************************/
/*                           Matrix  and  Matrix                              */
/******************************************************************************/

Matrix operator+(const Matrix &m1, const Matrix &m2) {
    assert(m1.nRow == m2.nRow);
    assert(m1.nCol == m2.nCol);

    if (m1.nRow != m2.nRow) {
        std::cerr << "Matrix operator+(const Matrix &m1, const Matrix &m2)" << std::endl;
        std::cerr << "error: m1.nRow != m2.nRow" << std::endl;
        std::abort();
    }

    if (m1.nCol != m2.nCol) {
        std::cerr << "Matrix operator+(const Matrix &m1, const Matrix &m2)" << std::endl;
        std::cerr << "error: m1.nCol != m2.nCol" << std::endl;
        std::abort();
    }

    const std::size_t M = m1.nRow;
    const std::size_t N = m1.nCol;
    
    Matrix ret;
    ret.nRow = M;
    ret.nCol = N;

    std::size_t nSize = M * N;
    if (nSize == 0) return ret;

    ret._data = reinterpret_cast<double*>(mkl_malloc(nSize*sizeof(double), 64));
    
    if (ret._data == NULL) {
        std::cerr << "Matrix operator+(const Matrix &m1, const Matrix &m2)" << std::endl;
        std::cerr << "error: can't allocate memory." << std::endl;
        std::abort();
    }

    vdAdd((const MKL_INT) nSize, 
          (const double*) m1._data, 
          (const double*) m2._data, 
          (double*) ret._data);

    return ret;
}

Matrix operator-(const Matrix &m1, const Matrix &m2) {
    assert(m1.nRow == m2.nRow);
    assert(m1.nCol == m2.nCol);

    if (m1.nRow != m2.nRow) {
        std::cerr << "Matrix operator-(const Matrix &m1, const Matrix &m2)" << std::endl;
        std::cerr << "error: m1.nRow != m2.nRow" << std::endl;
        std::abort();
    }

    if (m1.nCol != m2.nCol) {
        std::cerr << "Matrix operator-(const Matrix &m1, const Matrix &m2)" << std::endl;
        std::cerr << "error: m1.nCol != m2.nCol" << std::endl;
        std::abort();
    }

    const std::size_t M = m1.nRow;
    const std::size_t N = m1.nCol;
    
    Matrix ret;
    ret.nRow = M;
    ret.nCol = N;

    std::size_t nSize = M * N;
    if (nSize == 0) return ret;

    ret._data = reinterpret_cast<double*>(mkl_malloc(nSize*sizeof(double), 64));
    
    if (ret._data == NULL) {
        std::cerr << "Matrix operator-(const Matrix &m1, const Matrix &m2)" << std::endl;
        std::cerr << "error: can't allocate memory." << std::endl;
        std::abort();
    }

    vdSub((const MKL_INT) nSize, 
          (const double*) m1._data, 
          (const double*) m2._data, 
          (double*) ret._data);

    return ret;
}

Matrix operator*(const Matrix &m1, const Matrix &m2) {
    assert(m1.nRow == m2.nRow);
    assert(m1.nCol == m2.nCol);

    if (m1.nRow != m2.nRow) {
        std::cerr << "Matrix operator*(const Matrix &m1, const Matrix &m2)" << std::endl;
        std::cerr << "error: m1.nRow != m2.nRow" << std::endl;
        std::abort();
    }

    if (m1.nCol != m2.nCol) {
        std::cerr << "Matrix operator*(const Matrix &m1, const Matrix &m2)" << std::endl;
        std::cerr << "error: m1.nCol != m2.nCol" << std::endl;
        std::abort();
    }

    const std::size_t M = m1.nRow;
    const std::size_t N = m1.nCol;
    
    Matrix ret;
    ret.nRow = M;
    ret.nCol = N;

    std::size_t nSize = M * N;
    if (nSize == 0) return ret;

    ret._data = reinterpret_cast<double*>(mkl_malloc(nSize*sizeof(double), 64));
    
    if (ret._data == NULL) {
        std::cerr << "Matrix operator*(const Matrix &m1, const Matrix &m2)" << std::endl;
        std::cerr << "error: can't allocate memory." << std::endl;
        std::abort();
    }

    vdMul((const MKL_INT) nSize, 
          (const double*) m1._data, 
          (const double*) m2._data, 
          (double*) ret._data);

    return ret;
}

Matrix operator/(const Matrix &m1, const Matrix &m2) {
    assert(m1.nRow == m2.nRow);
    assert(m1.nCol == m2.nCol);

    if (m1.nRow != m2.nRow) {
        std::cerr << "Matrix operator/(const Matrix &m1, const Matrix &m2)" << std::endl;
        std::cerr << "error: m1.nRow != m2.nRow" << std::endl;
        std::abort();
    }

    if (m1.nCol != m2.nCol) {
        std::cerr << "Matrix operator/(const Matrix &m1, const Matrix &m2)" << std::endl;
        std::cerr << "error: m1.nCol != m2.nCol" << std::endl;
        std::abort();
    }

    const std::size_t M = m1.nRow;
    const std::size_t N = m1.nCol;
    
    Matrix ret;
    ret.nRow = M;
    ret.nCol = N;

    std::size_t nSize = M * N;
    if (nSize == 0) return ret;

    ret._data = reinterpret_cast<double*>(mkl_malloc(nSize*sizeof(double), 64));
    
    if (ret._data == NULL) {
        std::cerr << "Matrix operator/(const Matrix &m1, const Matrix &m2)" << std::endl;
        std::cerr << "error: can't allocate memory." << std::endl;
        std::abort();
    }

    vdDiv((const MKL_INT) nSize, 
          (const double*) m1._data, 
          (const double*) m2._data, 
          (double*) ret._data);

    return ret;
}

// matrix product, assign: m1.nCol == m2.nRow
Matrix operator%(const Matrix &m1, const Matrix &m2) {
    if (m1.nCol != m2.nRow) {
        std::cerr << "Matrix operator%(const Matrix &m1, const Matrix &m2)" << std::endl;
        std::cerr << "error: m1.nCol != m2.nRow" << std::endl;
        std::abort();
    }

    const std::size_t M = m1.nRow;
    const std::size_t N = m2.nCol;
    const std::size_t K = m1.nCol;

    if (M == 0 || N == 0 || K == 0) {
        return Matrix(M, N);
    }

    Matrix ret;
    ret.nRow = M;
    ret.nCol = N;

    std::size_t nSize = ret.nRow * ret.nCol;
    if (nSize == 0) return ret;

    ret._data = reinterpret_cast<double*>(mkl_malloc(nSize * sizeof(double), 64));

    if (ret._data == NULL) {
        std::cerr << "Matrix operator%(const Matrix &m1, const Matrix &m2)" << std::endl;
        std::cerr << "error: can't allocate memory." << std::endl;
        std::abort();
    }

    // C := alpha*op(A)*op(B) + beta*C
    cblas_dgemm(
        CblasRowMajor,                  // Layout
        CblasNoTrans,                   // transa
        CblasNoTrans,                   // transb
        (const MKL_INT) m1.nRow,        // m
        (const MKL_INT) m2.nCol,        // n
        (const MKL_INT) m1.nCol,        // k
        1.0,                            // alpha
        m1._data, (MKL_INT) m1.nCol,    // A, lda
        m2._data, (MKL_INT) m2.nCol,    // B, ldb
        0.0,                            // beta
        ret._data, (MKL_INT) ret.nCol   // C, ldc
    );

    return ret;
}


// HC = CE,
// where H is symmetric.
SelfAdjointEigenSolver::
SelfAdjointEigenSolver(const Matrix &H) 
: eigenvalues(H.nRow, 1), eigenvectors(H)
{
    if (H.nRow != H.nCol) {
        std::cerr << "SelfAdjointEigenSolver::"
                     "SelfAdjointEigenSolver(const Matrix &H)" << std::endl;
        std::cerr << "error: H is not square." << std::endl;
        std::abort();
    }

    if (H.nRow == 0) return;

    lapack_int info = -1;
    info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U',
                         (lapack_int) H.nRow,
                         (double *) eigenvectors._data, 
                         (lapack_int) eigenvectors.nCol,
                         (double *) eigenvalues._data);

    if (info != 0) {
        std::cerr << "SelfAdjointEigenSolver::"
                     "SelfAdjointEigenSolver(const Matrix &H)" << std::endl;
        std::cerr << "error: LAPACKE_dsyev failed with info = " << info << std::endl;
        std::abort();
    }
}


// HC = SCE, 
// where H is symmetric, S is positive definite.
GeneralSelfAdjointEigenSolver::
GeneralSelfAdjointEigenSolver(const Matrix &H, const Matrix &S)
: eigenvalues(H.nRow, 1), eigenvectors(H)
{
    if (H.nRow != H.nCol) {
        std::cerr << "GeneralSelfAdjointEigenSolver::"
                     "GeneralSelfAdjointEigenSolver(const Matrix &H, const Matrix &S)" << std::endl;
        std::cerr << "error: H is not square." << std::endl;
        std::abort();
    }

    if (S.nRow != S.nCol) {
        std::cerr << "GeneralSelfAdjointEigenSolver::"
                     "GeneralSelfAdjointEigenSolver(const Matrix &H, const Matrix &S)" << std::endl;
        std::cerr << "error: S is not square." << std::endl;
        std::abort();
    }

    if (H.nRow != S.nRow) {
        std::cerr << "GeneralSelfAdjointEigenSolver::"
                     "GeneralSelfAdjointEigenSolver(const Matrix &H, const Matrix &S)" << std::endl;
        std::cerr << "error: H.nRow != S.nRow." << std::endl;
        std::abort();
    }

    if (H.nRow == 0) return;

    std::size_t N = H.nRow;

    // LAPACKE_dsygv will change S, so make a backup
    Matrix Scp(S);

    // lapack_int LAPACKE_dsygv (int matrix_layout, lapack_int itype, 
    //                           char jobz, char uplo, lapack_int n, 
    //                           double* a, lapack_int lda, 
    //                           double* b, lapack_int ldb, double* w);

    lapack_int info = -1;
    info = LAPACKE_dsygv(LAPACK_ROW_MAJOR, 
                         (lapack_int) 1,    // itype=1: A*x = lambda*B*x
                         'V',               // calc eigenvectors also
                         'U',               // use upper triangle of H
                         (lapack_int) N,    // shape of H is (N, N) 
                         (double *) eigenvectors._data,  (lapack_int) N, 
                         (double *) Scp._data,           (lapack_int) N, 
                         (double *) eigenvalues._data );
    
    if (info != 0) {
        std::cerr << "GeneralSelfAdjointEigenSolver::"
                     "GeneralSelfAdjointEigenSolver(const Matrix &H, const Matrix &S)" << std::endl;
        std::cerr << "error: LAPACKE_dsygv failed with info = " << info << std::endl;
        std::abort();
    }
}