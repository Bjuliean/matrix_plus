#ifndef SRC_S21_MATRIX_OOP_H_
#define SRC_S21_MATRIX_OOP_H_

#include <iostream>
#include <cmath>

class S21Matrix {
    private:
        int rows_, cols_;
        double **matrix_;

        void create_matrix();
        void delete_matrix();
        // bool is_correct_matrix(const S21Matrix &other);        

    public:
        S21Matrix();
        ~S21Matrix();
        //S21Matrix(int n);
        S21Matrix(int rows, int cols);
        S21Matrix(const S21Matrix& other);
        S21Matrix(S21Matrix&& other);

        void init();
        void debug();
        void copy_matrix(const S21Matrix &other);
        bool is_correct_matrix() const;
        bool is_suitable_matrix(const S21Matrix &other);
        void fill_matrix();
        S21Matrix get_minor_matrix(const S21Matrix &other, int oi, int oz);

        S21Matrix operator+(const S21Matrix &other);
        S21Matrix operator-(const S21Matrix &other);
        S21Matrix operator*(const S21Matrix &other);
        S21Matrix& operator+=(const S21Matrix &other);
        S21Matrix& operator-=(const S21Matrix &other);
        S21Matrix& operator*=(const S21Matrix &other); // get set
        S21Matrix operator=(const S21Matrix &other);
        bool operator ==(const S21Matrix &other);

        bool EqMatrix(const S21Matrix& other);
        void SumMatrix(const S21Matrix& other);
        void SubMatrix(const S21Matrix& other);
        void MulNumber(const double num);
        void MulMatrix(const S21Matrix& other);
        S21Matrix Transpose();
        S21Matrix CalcComplements();
        double Determinant();
        S21Matrix InverseMatrix();
};


#endif // SRC_S21_MATRIX_OOP_H_