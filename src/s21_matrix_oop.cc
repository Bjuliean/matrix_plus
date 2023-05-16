#include "s21_matrix_oop.h"

S21Matrix::S21Matrix(int rows) : rows_(rows), cols_(rows) {
  if (rows_ < 1)
    throw std::logic_error(
        "Wrong matrix size. Rows and columns should be >= 1");
  CreateMatrix();
}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows_ < 1 || cols_ < 1)
    throw std::logic_error(
        "Wrong matrix size. Rows and columns should be >= 1");
  CreateMatrix();
}

S21Matrix::S21Matrix(const S21Matrix &other)
    : rows_(other.rows_), cols_(other.cols_) {
  CreateMatrix();
  CopyMatrix(other);
}

S21Matrix::S21Matrix(S21Matrix &&other) noexcept
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

double &S21Matrix::operator()(int i, int j) {
  if (i < 0 || j < 0 || i >= rows_ || j >= cols_)
    throw std::logic_error("Invalid index");
  return matrix_[i][j];
}

S21Matrix S21Matrix::operator+(const S21Matrix &other) {
  S21Matrix nw(*this);
  nw.SumMatrix(other);
  return nw;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) {
  S21Matrix nw(*this);
  nw.SubMatrix(other);
  return nw;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) {
  S21Matrix nw(*this);
  nw.MulMatrix(other);
  return nw;
}

S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix nw(*this);
  nw.MulNumber(num);
  return nw;
}

S21Matrix &S21Matrix::operator+=(const S21Matrix &other) {
  SumMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator-=(const S21Matrix &other) {
  SubMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator*=(const S21Matrix &other) {
  MulMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator*=(const double num) {
  MulNumber(num);
  return *this;
}

S21Matrix S21Matrix::operator=(const S21Matrix &other) {
  if (!(this == &other)) {
    if (IsCorrectMatrix()) DeleteMatrix();
    rows_ = other.rows_;
    cols_ = other.cols_;
    CreateMatrix();
    CopyMatrix(other);
  }
  return *this;
}

void S21Matrix::SetRows(int rows) {
  if (rows < 1)
    throw std::logic_error("Attempt to set an incorrect matrix size");
  S21Matrix temp(*this);
  DeleteMatrix();
  rows_ = rows;
  cols_ = temp.cols_;
  CreateMatrix();
  MutatorCopy(temp);
}

void S21Matrix::SetCols(int cols) {
  if (cols < 1)
    throw std::logic_error("Attempt to set an incorrect matrix size");
  S21Matrix temp(*this);
  DeleteMatrix();
  cols_ = cols;
  rows_ = temp.rows_;
  CreateMatrix();
  MutatorCopy(temp);
}

bool S21Matrix::EqMatrix(const S21Matrix &other) {
  bool res = true;
  if (IsSuitableMatrix(other)) {
    for (int i = 0; i < rows_; i++) {
      for (int z = 0; z < cols_; z++) {
        if (fabs(this->matrix_[i][z] - other.matrix_[i][z]) > 1E-7) {
          res = false;
          break;
        }
      }
    }
  } else {
    res = false;
  }
  return res;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (!IsSuitableMatrix(other))
    throw std::logic_error("Addition of unsuitable matrices");
  for (int i = 0; i < rows_; i++) {
    for (int z = 0; z < cols_; z++) {
      matrix_[i][z] += other.matrix_[i][z];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (!IsSuitableMatrix(other))
    throw std::logic_error("Subtraction of unsuitable matrices");
  for (int i = 0; i < rows_; i++) {
    for (int z = 0; z < cols_; z++) {
      matrix_[i][z] -= other.matrix_[i][z];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++) {
    for (int z = 0; z < cols_; z++) {
      matrix_[i][z] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (cols_ != other.rows_)
    throw std::logic_error("Multiplication of unsuitable matrices");
  S21Matrix temp(rows_, other.cols_);
  for (int i = 0; i < rows_; i++) {
    for (int z = 0; z < other.cols_; z++) {
      for (int k = 0; k < cols_; k++) {
        temp.matrix_[i][z] += matrix_[i][k] * other.matrix_[k][z];
      }
    }
  }
  *this = temp;
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix temp(cols_, rows_);
  for (int i = 0; i < rows_; i++) {
    for (int z = 0; z < cols_; z++) {
      temp.matrix_[z][i] = matrix_[i][z];
    }
  }
  *this = temp;
  return *this;
}

S21Matrix S21Matrix::GetMinorMatrix(const S21Matrix &other, int oi, int oz) {
  S21Matrix nw(other.rows_ - 1, other.cols_ - 1);
  int x = 0, y = 0;
  for (int i = 0; i < other.rows_; i++) {
    x = 0;
    for (int z = 0; z < other.cols_; z++) {
      if (z != oz && i != oi) {
        nw.matrix_[y][x++] = other.matrix_[i][z];
      }
    }
    if (i != oi) y++;
  }
  return nw;
}

double S21Matrix::Determinant() {
  double res = 0;
  if (rows_ != cols_)
    throw std::logic_error(
        "Calculation of the determinant of a non-square matrix");
  if (this->rows_ == 1) {
    res = this->matrix_[0][0];
  } else if (this->rows_ == 2) {
    res = (this->matrix_[0][0] * this->matrix_[1][1]) -
          (this->matrix_[1][0] * this->matrix_[0][1]);
  } else {
    for (int i = 0; i < cols_; i++) {
      S21Matrix nw = GetMinorMatrix(*this, 0, i);
      res += pow(-1, i) * matrix_[0][i] * nw.Determinant();
    }
  }
  return res;
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_)
    throw std::logic_error(
        "Calculation of the matrix of algebraic additions in a non-square "
        "matrix");
  S21Matrix nw(rows_, cols_);
  if (rows_ == 1) {
    nw.matrix_[0][0] = matrix_[0][0];
  } else {
    for (int i = 0; i < nw.rows_; i++) {
      for (int z = 0; z < nw.cols_; z++) {
        S21Matrix temp = GetMinorMatrix(*this, i, z);
        nw.matrix_[i][z] = temp.Determinant() * pow(-1, i + z);
      }
    }
  }
  return nw;
}

S21Matrix S21Matrix::InverseMatrix() {
  if (rows_ != cols_ || this->Determinant() == 0)
    throw std::logic_error("Calculation of the inverse matrix is impossible");
  double det = Determinant();
  S21Matrix nw(rows_, cols_);
  if (rows_ == 1) {
    nw(0, 0) = 1 / det;
  } else {
    nw = this->CalcComplements();
    nw.Transpose();
    nw.MulNumber(1.0 / det);
  }
  return nw;
}

void S21Matrix::CreateMatrix() {
  if (rows_ > 0 && cols_ > 0) {
    matrix_ = new double *[rows_];
    for (int i = 0; i < rows_; i++) matrix_[i] = new double[cols_]();
  }
}

void S21Matrix::DeleteMatrix() {
  if (IsCorrectMatrix()) {
    for (int i = 0; i < rows_; i++) delete[] matrix_[i];
    delete[] matrix_;
    rows_ = 0;
    cols_ = 0;
    matrix_ = nullptr;
  }
}

bool S21Matrix::IsCorrectMatrix() const {
  return rows_ > 0 && cols_ > 0 && matrix_ != nullptr;
}

bool S21Matrix::IsSuitableMatrix(const S21Matrix &other) {
  bool res = true;
  if (!IsCorrectMatrix() || !other.IsCorrectMatrix()) res = false;
  if (rows_ != other.rows_ || cols_ != other.cols_) res = false;
  return res;
}

void S21Matrix::CopyMatrix(const S21Matrix &other) {
  if (IsSuitableMatrix(other)) {
    for (int i = 0; i < rows_; i++) {
      for (int z = 0; z < cols_; z++) {
        matrix_[i][z] = other.matrix_[i][z];
      }
    }
  }
}

void S21Matrix::MutatorCopy(const S21Matrix &other) {
  for (int i = 0; i < rows_; i++) {
    for (int z = 0; z < cols_; z++) {
      if (i >= other.rows_ || z >= other.cols_)
        matrix_[i][z] = 0;
      else
        matrix_[i][z] = other.matrix_[i][z];
    }
  }
}