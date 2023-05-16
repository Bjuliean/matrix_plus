#ifndef SRC_S21_MATRIX_OOP_H_
#define SRC_S21_MATRIX_OOP_H_

#include <cmath>
#include <iostream>

class S21Matrix {
 public:
  S21Matrix() : rows_(0), cols_(0), matrix_(nullptr){};
  S21Matrix(int rows);
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix &other);
  S21Matrix(S21Matrix &&other) noexcept;
  ~S21Matrix() { DeleteMatrix(); }

  double &operator()(int i, int j);
  S21Matrix operator+(const S21Matrix &other);
  S21Matrix operator-(const S21Matrix &other);
  S21Matrix operator*(const S21Matrix &other);
  S21Matrix operator*(const double num);
  friend S21Matrix operator*(double num, S21Matrix &matrix) {
    return matrix * num;
  };
  S21Matrix &operator+=(const S21Matrix &other);
  S21Matrix &operator-=(const S21Matrix &other);
  S21Matrix &operator*=(const S21Matrix &other);
  S21Matrix &operator*=(const double num);
  S21Matrix operator=(const S21Matrix &other);
  bool operator==(const S21Matrix &other) { return EqMatrix(other); }
  bool operator!=(const S21Matrix &other) { return !EqMatrix(other); }

  int GetRows() { return rows_; }
  int GetCols() { return cols_; }
  void SetRows(int rows);
  void SetCols(int cols);

  bool EqMatrix(const S21Matrix &other);
  void SumMatrix(const S21Matrix &other);
  void SubMatrix(const S21Matrix &other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix &other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();

 private:
  int rows_, cols_;
  double **matrix_;

  void CreateMatrix();
  void DeleteMatrix();
  void CopyMatrix(const S21Matrix &other);
  bool IsSuitableMatrix(const S21Matrix &other);
  void MutatorCopy(const S21Matrix &other);
  S21Matrix GetMinorMatrix(const S21Matrix &other, int oi, int oz);
  bool IsCorrectMatrix() const;
};

#endif  // SRC_S21_MATRIX_OOP_H_