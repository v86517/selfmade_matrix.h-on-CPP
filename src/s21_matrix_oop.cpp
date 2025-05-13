#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() : rows_(0), cols_(0), matrix_(nullptr) {}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows < 0 || cols < 0) {
    throw std::out_of_range("cols and rows should be more than zero");
  }
  matrix_ = new double*[rows_]();
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_]();
  }
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : rows_(other.rows_), cols_(other.cols_) {
  matrix_ = new double*[rows_]();
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_]();
  }
  int row = rows_ < other.rows_ ? rows_ : other.rows_;
  int col = cols_ < other.cols_ ? cols_ : other.cols_;
  for (int i = 0; i < row; i++)
    for (int j = 0; j < col; j++) matrix_[i][j] = other(i, j);
}

S21Matrix::S21Matrix(S21Matrix&& other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
}

S21Matrix::~S21Matrix() {
  if (matrix_ != nullptr) {
    for (int i = 0; i < rows_; i++) {
      delete[] matrix_[i];
    }
  }
  delete[] matrix_;
  rows_ = 0;
  cols_ = 0;
  matrix_ = nullptr;
}

bool S21Matrix::EqMatrix(const S21Matrix& other) const {
  bool res = true;

  if ((this->rows_ != other.rows_) || (this->cols_ != other.cols_))
    res = false;
  else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        if (s21_fabs(this->matrix_[i][j] - other(i, j)) > EPS) {
          res = false;
          break;
        }
      }

      if (res == false) break;
    }
  }

  return res;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::out_of_range("different matrix dimensions");

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      this->matrix_[i][j] += other(i, j);
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::out_of_range("different matrix dimensions");

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      this->matrix_[i][j] -= other(i, j);
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      this->matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.rows_) {
    throw std::out_of_range(
        "The number of columns of the first matrix is not equal to the number "
        "of rows of the second matrix");
  }
  S21Matrix res = S21Matrix(rows_, other.cols_);
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < other.cols_; j++)
      for (int k = 0; k < cols_; k++)
        res.matrix_[i][j] += matrix_[i][k] * other(k, j);
  *this = res;
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix res = S21Matrix(cols_, rows_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      res.matrix_[j][i] = matrix_[i][j];
    }
  }
  return res;
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_) {
    throw std::out_of_range("The matrix is ​​not square");
  }
  double temp = 0;
  S21Matrix res = S21Matrix(cols_, rows_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      S21Matrix minor_matrix = GetMinorMatrix(j, i);
      temp = minor_matrix.Determinant();
      minor_matrix.~S21Matrix();
      res.matrix_[i][j] = ((i + j) % 2 == 1) ? -temp : temp;
    }
  }
  return res;
}

S21Matrix S21Matrix::GetMinorMatrix(int col, int row) {
  S21Matrix res = S21Matrix(cols_ - 1, rows_ - 1);
  for (int i = 0, k = 0; i < this->rows_; ++i) {
    if (row != i) {
      for (int j = 0, l = 0; j < this->cols_; ++j) {
        if (col != j) {
          res.matrix_[k][l] = matrix_[i][j];
          l++;
        }
      }
      k++;
    }
  }
  return res;
}

double S21Matrix::Determinant() {
  if (rows_ != cols_) {
    throw std::out_of_range("The matrix is not square");
  }
  double res = 0, temp = 0;
  if (cols_ == 1) {
    res = matrix_[0][0];
  } else if (cols_ == 2) {
    res = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  } else {
    for (int i = 0; i < cols_; i++) {
      S21Matrix minor_matrix = GetMinorMatrix(i, 0);
      temp = minor_matrix.Determinant();
      minor_matrix.~S21Matrix();
      res = (i % 2 == 1) ? res - matrix_[0][i] * temp
                         : res + matrix_[0][i] * temp;
    }
  }
  return res;
}

S21Matrix S21Matrix::InverseMatrix() {
  double determin = Determinant();
  if (s21_fabs(determin) < EPS) {
    throw std::out_of_range("Matrix determinant is 0");
  }
  S21Matrix res(*this);
  if (rows_ == 1 && cols_ == 1) {
    res.matrix_[0][0] = 1 / matrix_[0][0];
  } else {
    S21Matrix alg_dop = CalcComplements();
    res = alg_dop.Transpose();
    res.MulNumber(1 / determin);
  }
  return res;
}

double S21Matrix::operator()(int i, int j) const {
  // if ((i < 0 || j < 0) || ((i > (rows_ - 1)) || (j > (cols_ - 1))))
  //   throw std::out_of_range("Index is outside the matrix");
  return matrix_[i][j];
}

double& S21Matrix::operator()(int i, int j) {
  if ((i < 0 || j < 0) || ((i > (rows_ - 1)) || (j > (cols_ - 1))))
    throw std::out_of_range("Index is outside the matrix");
  return matrix_[i][j];
}

S21Matrix S21Matrix::operator=(const S21Matrix& other) {
  S21Matrix temp(other);

  for (int i = 0; i < rows_; i++) {
    delete[] matrix_[i];
  }
  delete[] matrix_;

  rows_ = temp.rows_;
  cols_ = temp.cols_;
  matrix_ = temp.matrix_;
  temp.matrix_ = nullptr;

  return *this;
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) const {
  S21Matrix res(*this);
  res.SumMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) const {
  S21Matrix res(*this);
  res.SubMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) const {
  S21Matrix res(*this);
  res.MulMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator*(const double num) const {
  S21Matrix res(*this);
  res.MulNumber(num);
  return res;
}

bool S21Matrix::operator==(const S21Matrix& other) const {
  return EqMatrix(other);
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  SumMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  SubMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  MulMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const double num) {
  MulNumber(num);
  return *this;
}

int S21Matrix::GetRows() const { return rows_; }

int S21Matrix::GetCols() const { return cols_; }

void S21Matrix::SetRows(int rows) {
  if (rows < 0) {
    throw std::out_of_range("Rows should be more than zero");
  }
  S21Matrix tmp(rows, cols_);
  int r = rows <= rows_ ? rows : rows_;
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < cols_; j++) {
      tmp.matrix_[i][j] = matrix_[i][j];
    }
  }
  *this = tmp;
}

void S21Matrix::SetCols(int cols) {
  if (cols < 0) {
    throw std::out_of_range("Cols should be more than zero");
  }
  S21Matrix tmp(rows_, cols);
  int c = cols <= cols_ ? cols : cols_;
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < c; j++) {
      tmp.matrix_[i][j] = matrix_[i][j];
    }
  }
  *this = tmp;
}

long double s21_fabs(double x) {
  x = x >= 0 ? x : -x;
  return (long double)x;
}
