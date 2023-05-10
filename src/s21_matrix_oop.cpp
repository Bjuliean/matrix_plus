#include "s21_matrix_oop.h"

// int main(void) {

//     return 0;
// }

// void S21Matrix::init() {
//     for(int i = 0; i < rows_; i++) {
//         for(int z = 0; z < cols_; z++) {
//             std::cin >> matrix_[i][z];
//         }
//     }
// }

// void S21Matrix::debug() {
//     std::cout << std::endl;
//     for(int i = 0; i < rows_; i++) {
//         for(int z = 0; z < cols_; z++) {
//             std::cout << matrix_[i][z] << ' ';
//         }
//         std::cout << std::endl;
//     }
//     std::cout << std::endl;
// }

void S21Matrix::create_matrix() {
    if(rows_ < 1 || cols_ < 1)
        throw std::logic_error("Trying to create wrong matrix");
    matrix_ = new double* [rows_];
    for(int i = 0; i < rows_; i++)
        matrix_[i] = new double [cols_];
    fill_matrix();
}

void S21Matrix::delete_matrix() {
    if(is_correct_matrix()) {
        for(int i = 0; i < rows_; i++)
            delete [] matrix_[i];
        delete [] matrix_;
        rows_ = 0;
        cols_ = 0;
        matrix_ = nullptr;
    }
}

void S21Matrix::fill_matrix() {
    for(int i = 0; i < rows_; i++) {
        for(int z = 0; z < cols_; z++) {
            matrix_[i][z] = 0;
        }
    }
}

bool S21Matrix::is_correct_matrix() const {
    return rows_ > 0 && cols_ > 0 && matrix_ != nullptr;
}

bool S21Matrix::is_suitable_matrix(const S21Matrix &other) {
    bool res = true;
    if(!is_correct_matrix() || !other.is_correct_matrix())
        res = false;
    if(rows_ != other.rows_ || cols_ != other.cols_)
        res = false;
    return res;
}

void S21Matrix::copy_matrix(const S21Matrix &other) {
    if(is_suitable_matrix(other)) {
        for(int i = 0; i < rows_; i++) {
            for(int z = 0; z < cols_; z++) {
                matrix_[i][z] = other.matrix_[i][z];
            }
        }
    }
}

void S21Matrix::mutator_copy(const S21Matrix &other) {
    for(int i = 0; i < rows_; i++) {
        for(int z = 0; z < cols_; z++) {
            if(i >= other.rows_ || z >= other.cols_)
                matrix_[i][z] = 0;
            else
                matrix_[i][z] = other.matrix_[i][z];
        }
    }
}

S21Matrix::S21Matrix() {
    rows_ = 2;
    cols_ = 2;
    create_matrix();
}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
    if(rows_ < 1 || cols_ < 1)
        throw std::logic_error("Wrong matrix size. Rows and columns should be >= 1");
    create_matrix();
}

S21Matrix::S21Matrix(const S21Matrix &other) : rows_(other.rows_), cols_(other.cols_) {
    if(rows_ < 1 || cols_ < 1)
        throw std::logic_error("Copying from the wrong matrix");
    create_matrix();
    copy_matrix(other);
}

S21Matrix::S21Matrix(S21Matrix &&other) : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
    other.rows_ = 0;
    other.cols_ = 0;
    other.matrix_ = nullptr;
}

S21Matrix::~S21Matrix() {
    delete_matrix();
}

bool S21Matrix::EqMatrix(const S21Matrix &other) {
    bool res = true;
    if(is_suitable_matrix(other)) {
        for(int i = 0; i < rows_; i++) {
            for(int z = 0; z < cols_; z++) {
                if(this->matrix_[i][z] - other.matrix_[i][z] > 1E-7) {
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
    if(!is_suitable_matrix(other))
        throw std::logic_error("Addition of unsuitable matrices");
    for(int i = 0; i < rows_; i++) {
        for(int z = 0; z < cols_; z++) {
            matrix_[i][z] += other.matrix_[i][z];
        }
    }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
    if(!is_suitable_matrix(other))
        throw std::logic_error("Subtraction of unsuitable matrices");
    for(int i = 0; i < rows_; i++) {
        for(int z = 0; z < cols_; z++) {
            matrix_[i][z] -= other.matrix_[i][z];
        }
    }
}

void S21Matrix::MulNumber(const double num) {
    for(int i = 0; i < rows_; i++) {
        for(int z = 0; z < cols_; z++) {
            matrix_[i][z] *= num;
        }
    }
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
    if(cols_ != other.rows_)
        throw std::logic_error("Multiplication of unsuitable matrices");
    S21Matrix other_copy(other);
    other_copy.Transpose();
    S21Matrix temp(rows_, rows_);
    for(int i = 0; i < rows_; i++) {
        for(int j = 0; j < rows_; j++) {
            for(int k = 0; k < cols_; k++) {
                temp.matrix_[i][j] += matrix_[i][k] * other_copy.matrix_[j][k];
            }
        }
    }
    *this = temp;
}

S21Matrix S21Matrix::Transpose() {
    S21Matrix temp(cols_, rows_);
    for(int i = 0; i < rows_; i++) {
        for(int z = 0; z < cols_; z++) {
            temp.matrix_[z][i] = matrix_[i][z];
        }
    }
    *this = temp;
    return *this;
}

S21Matrix S21Matrix::get_minor_matrix(const S21Matrix &other, int oi, int oz) {
    S21Matrix nw(other.rows_ - 1, other.cols_ - 1);
    int x = 0, y = 0;
    for(int i = 0; i < other.rows_; i++) {
        x = 0;
        for(int z = 0; z < other.cols_; z++) {
            if(z != oz && i != oi) {
                nw.matrix_[y][x++] = other.matrix_[i][z];
            }
        }
        if(i != oi)
            y++;
    }
    return nw;
}

double S21Matrix::Determinant() {
    double res = 0;
    if(rows_ != cols_)
        throw std::logic_error("Calculation of the determinant of a non-square matrix");
    if(this->rows_ == 1) {
        res = this->matrix_[0][0];
    }
    else if(this->rows_ == 2) {
        res = (this->matrix_[0][0] * this->matrix_[1][1]) - (this->matrix_[1][0] * this->matrix_[0][1]);
    }
    else {
        for(int i = 0; i < cols_; i++) {
            S21Matrix nw = get_minor_matrix(*this, 0, i);
            res += pow(-1, i) * matrix_[0][i] * nw.Determinant();
        }
    }
    return res;
}

S21Matrix S21Matrix::CalcComplements() {
    if(rows_ != cols_)
        throw std::logic_error("Calculation of the matrix of algebraic additions in a non-square matrix");
    S21Matrix nw(rows_, cols_);
    if(rows_ == 1) {
        nw.matrix_[0][0] = matrix_[0][0];
    } else {
        for(int i = 0; i < nw.rows_; i++) {
            for(int z = 0; z < nw.cols_; z++) {
                S21Matrix temp = get_minor_matrix(*this, i, z);
                nw.matrix_[i][z] = temp.Determinant() * pow(-1, i + z);
            }
        }
    }
    return nw;
}

S21Matrix S21Matrix::InverseMatrix() {
    if(rows_ != cols_ || this->Determinant() == 0)
        throw std::logic_error("Calculation of the inverse matrix is impossible");
    double det = Determinant();
    S21Matrix nw(rows_, cols_);
    nw = this->CalcComplements();
    nw.Transpose();
    nw.MulNumber(1.0 / det);
    return nw;
}

int S21Matrix::GetRows() {
    return rows_;
}

int S21Matrix::GetCols() {
    return cols_;
}

void S21Matrix::SetRows(int rows) {
    if(rows < 1)
        throw std::logic_error("Attempt to set an incorrect matrix size");
    S21Matrix temp(*this);
    delete_matrix();
    rows_ = rows;
    cols_ = temp.cols_;
    create_matrix();
    mutator_copy(temp);
}

void S21Matrix::SetCols(int cols) {
    if(cols < 1)
        throw std::logic_error("Attempt to set an incorrect matrix size");
    S21Matrix temp(*this);
    delete_matrix();
    cols_ = cols;
    rows_ = temp.rows_;
    create_matrix();
    mutator_copy(temp);
}

double& S21Matrix::operator ()(int i, int j) {
    if(i < 0 || j < 0 || i >= rows_ || j >= cols_)
        throw std::logic_error("Invalid index");
    return matrix_[i][j];
}

bool S21Matrix::operator ==(const S21Matrix &other) {
    return EqMatrix(other);
}

S21Matrix& S21Matrix::operator +=(const S21Matrix &other) {
    SumMatrix(other);
    return *this;
}

S21Matrix& S21Matrix::operator -=(const S21Matrix &other) {
    SubMatrix(other);
    return *this;
}

S21Matrix& S21Matrix::operator *=(const S21Matrix &other) {
    MulMatrix(other);
    return *this;
}

S21Matrix& S21Matrix::operator *=(const double num) {
    MulNumber(num);
    return *this;
}

S21Matrix S21Matrix::operator +(const S21Matrix &other) {
    S21Matrix nw(*this);
    nw.SumMatrix(other);
    return nw;
}

S21Matrix S21Matrix::operator -(const S21Matrix &other) {
    S21Matrix nw(*this);
    nw.SubMatrix(other);
    return nw;
}

S21Matrix S21Matrix::operator *(const S21Matrix &other) {
    S21Matrix nw(*this);
    nw.MulMatrix(other);
    return nw;
}

S21Matrix S21Matrix::operator *(const double num) {
    S21Matrix nw(*this);
    nw.MulNumber(num);
    return nw;
}

S21Matrix operator*(double num, S21Matrix& matrix) {
    return matrix * num;
}

S21Matrix S21Matrix::operator =(const S21Matrix &other) {
    if(!(this == &other)) {
        if(is_correct_matrix())
            delete_matrix();
        rows_ = other.rows_;
        cols_ = other.cols_;
        create_matrix();
        copy_matrix(other);
    }
    return *this;
}