#include "s21_matrix_oop.h"

int main() {

    int a, b;
    std::cin >> a;
    std::cin >> b;
    std::cout << std::endl;
    S21Matrix m(a, b);
    m.init();

    // m.Transpose();
    // m.debug();
    std::cout << std::endl;
    std::cin >> a;
    std::cin >> b;
    S21Matrix n(a, b);
    n.init();
    m.MulMatrix(n);

    m.debug();

    return 0;
}

void S21Matrix::init() {
    for(int i = 0; i < rows_; i++) {
        for(int z = 0; z < cols_; z++) {
            std::cin >> matrix_[i][z];
        }
    }
}

void S21Matrix::debug() {
    std::cout << std::endl;
    for(int i = 0; i < rows_; i++) {
        for(int z = 0; z < cols_; z++) {
            std::cout << matrix_[i][z] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void S21Matrix::create_matrix() {
    if(!(rows_ > 0 && cols_ > 0))
        throw std::logic_error("Trying to create wrong matrix");
    matrix_ = new double* [rows_];
    for(int i = 0; i < rows_; i++)
        matrix_[i] = new double [cols_];
    fill_matrix();
}

void S21Matrix::delete_matrix() {
    if(!is_correct_matrix())
        throw std::logic_error("Memory free error");
    for(int i = 0; i < rows_; i++)
        delete [] matrix_[i];
    delete [] matrix_;
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

// bool S21Matrix::is_correct_matrix(const S21Matrix &other) {
//     return is_correct_matrix(matrix) && is_correct_matrix(other);
// }

S21Matrix::S21Matrix() {
    rows_ = 2;
    cols_ = 2;
    create_matrix();
}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
    if(!(rows_ > 0 && cols_ > 0))
        throw std::logic_error("Wrong matrix size. Rows and columns should be >= 1");
    create_matrix();
}

S21Matrix::S21Matrix(const S21Matrix &other) : rows_(other.rows_), cols_(other.cols_) {
    if(!(rows_ > 0 && cols_ > 0))
        throw std::logic_error("Copying from the wrong matrix");
    rows_ = other.rows_;
    cols_ = other.cols_;
    create_matrix();
    copy_matrix(other);
}

S21Matrix::~S21Matrix() {
    delete_matrix();
}

bool S21Matrix::EqMatrix(const S21Matrix &other) {
    bool res = true;
    if(is_suitable_matrix(other)) {
        for(int i = 0; i < rows_; i++) {
            for(int z = 0; z < cols_; z++) {
                if(this->matrix_[i][z] != other.matrix_[i][z]) {
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
    // std::cout << "1DEBUG!!!" << std::endl;
    // other_copy.debug();
    other_copy.Transpose();
    S21Matrix temp(rows_, rows_);
    // std::cout << "DEBUG!!!" << std::endl;
    // other_copy.debug();
    for(int i = 0; i < rows_; i++) {
        for(int j = 0; j < rows_; j++) {
            for(int k = 0; k < cols_; k++) {
                //std::cout << matrix_[i][k] << "*" << other_copy.matrix_[j][k] << "+";
                //std::cout << "(" << temp.matrix_[i][j] << ")" << matrix_[i][k] << "*" << other_copy.matrix_[j][k] << "=" << matrix_[i][k] * other_copy.matrix_[j][k] << " + ";
                temp.matrix_[i][j] += matrix_[i][k] * other_copy.matrix_[j][k];
            }
            //std::cout << "\t";
        }
        //std::cout << std::endl;
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

bool S21Matrix::operator ==(const S21Matrix &other) {
    return EqMatrix(other);
}

S21Matrix S21Matrix::operator +(const S21Matrix &other) {
    SumMatrix(other);
    return *this;
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