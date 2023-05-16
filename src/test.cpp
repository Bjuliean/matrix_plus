#include <gtest/gtest.h>

#include "s21_matrix_oop.h"

TEST(EqMatrix, True) {
  S21Matrix matrix_a(3, 3);
  S21Matrix matrix_b(3, 3);
  ASSERT_TRUE(matrix_a == matrix_b);
}
TEST(EqMatrix, False) {
  S21Matrix matrix_a(3, 3);
  S21Matrix matrix_b(2, 2);
  ASSERT_FALSE(matrix_a == matrix_b);
}
TEST(SumMatrix, True) {
  S21Matrix matrix_a(2, 2);
  S21Matrix matrix_b(2, 2);
  S21Matrix result(2, 2);

  matrix_a(0, 0) = 3.14;
  matrix_a(0, 1) = 0.56;
  matrix_a(1, 0) = -69.3;
  matrix_a(1, 1) = 0;

  matrix_b(0, 0) = -78.14;
  matrix_b(0, 1) = 0;
  matrix_b(1, 0) = -0.3;
  matrix_b(1, 1) = 2;

  result(0, 0) = -75;
  result(0, 1) = 0.56;
  result(1, 0) = -69.6;
  result(1, 1) = 2;
  matrix_a.SumMatrix(matrix_b);
  ASSERT_TRUE(matrix_a == result);
}
TEST(SumMatrix, False) {
  S21Matrix matrix_a(1, 2);
  S21Matrix matrix_b(2, 2);

  matrix_a(0, 0) = 3.14;
  matrix_a(0, 1) = 0.56;

  matrix_b(0, 0) = -78.14;
  matrix_b(0, 1) = 0;
  matrix_b(1, 0) = -0.3;
  matrix_b(1, 1) = 2;

  EXPECT_THROW(matrix_a.SumMatrix(matrix_b), std::logic_error);
}
TEST(SubMatrix, True) {
  S21Matrix matrix_a(2, 2);
  S21Matrix matrix_b(2, 2);
  S21Matrix result(2, 2);

  matrix_a(0, 0) = 3;
  matrix_a(0, 1) = 2;
  matrix_a(1, 0) = -6;
  matrix_a(1, 1) = 0;

  matrix_b(0, 0) = -7;
  matrix_b(0, 1) = 0;
  matrix_b(1, 0) = -3;
  matrix_b(1, 1) = 2;

  result(0, 0) = 10;
  result(0, 1) = 2;
  result(1, 0) = -3;
  result(1, 1) = -2;

  matrix_a.SubMatrix(matrix_b);

  ASSERT_TRUE(matrix_a == result);
}
TEST(SubMatrix, False) {
  S21Matrix matrix_a(1, 2);
  S21Matrix matrix_b(2, 2);

  matrix_a(0, 0) = 3;
  matrix_a(0, 1) = 2;

  matrix_b(0, 0) = -7;
  matrix_b(0, 1) = 0;
  matrix_b(1, 0) = -3;
  matrix_b(1, 1) = 2;

  EXPECT_THROW(matrix_a.SubMatrix(matrix_b), std::logic_error);
}
TEST(MulNumber, True) {
  S21Matrix matrix_a(2, 2);
  S21Matrix result(2, 2);

  matrix_a(0, 0) = -78.14;
  matrix_a(0, 1) = 0;
  matrix_a(1, 0) = -0.3;
  matrix_a(1, 1) = 2;

  result(0, 0) = -781.4;
  result(0, 1) = 0;
  result(1, 0) = -3;
  result(1, 1) = 20;

  matrix_a.MulNumber(10);

  ASSERT_TRUE(matrix_a == result);
}
TEST(MulMatrix, True) {
  S21Matrix matrix_a(2, 2);
  S21Matrix matrix_b(2, 2);
  S21Matrix result(2, 2);

  matrix_a(0, 0) = 3;
  matrix_a(0, 1) = 2;
  matrix_a(1, 0) = -6.6;
  matrix_a(1, 1) = 0;

  matrix_b(0, 0) = -7;
  matrix_b(0, 1) = 0;
  matrix_b(1, 0) = -3.5;
  matrix_b(1, 1) = 2;

  result(0, 0) = -28;
  result(0, 1) = 4;
  result(1, 0) = 46.2;
  result(1, 1) = 0;

  matrix_a.MulMatrix(matrix_b);

  ASSERT_TRUE(matrix_a == result);
}
TEST(MulMatrix, False) {
  S21Matrix matrix_a(2, 1);
  S21Matrix matrix_b(2, 2);

  matrix_a(0, 0) = 3;
  matrix_a(1, 0) = -6.6;

  matrix_b(0, 0) = -7;
  matrix_b(0, 1) = 0;
  matrix_b(1, 0) = -3.5;
  matrix_b(1, 1) = 2;

  EXPECT_THROW(matrix_a.MulMatrix(matrix_b), std::logic_error);
}
TEST(OperatorParentheses, True) {
  S21Matrix matrix_a(2, 2);

  matrix_a(0, 0) = 3;
  matrix_a(0, 1) = 2;
  matrix_a(1, 0) = -6.6;
  matrix_a(1, 1) = 0;
  ASSERT_EQ(matrix_a(0, 1), 2);
}
TEST(OperatorParentheses, False) {
  S21Matrix matrix_a(2, 2);

  matrix_a(0, 0) = 3;
  matrix_a(0, 1) = 2;
  matrix_a(1, 0) = -6.6;
  matrix_a(1, 1) = 0;
  ASSERT_NE(matrix_a(0, 1), 10);
}
TEST(OperatorPlus, True) {
  S21Matrix matrix_a(2, 2);
  S21Matrix matrix_b(2, 2);
  S21Matrix result(2, 2);

  matrix_a(0, 0) = 3.14;
  matrix_a(0, 1) = 0.56;
  matrix_a(1, 0) = -69.3;
  matrix_a(1, 1) = 0;

  matrix_b(0, 0) = -78.14;
  matrix_b(0, 1) = 0;
  matrix_b(1, 0) = -0.3;
  matrix_b(1, 1) = 2;

  result(0, 0) = -75;
  result(0, 1) = 0.56;
  result(1, 0) = -69.6;
  result(1, 1) = 2;

  ASSERT_TRUE((matrix_a + matrix_b) == result);
}
TEST(OperatorMinus, True) {
  S21Matrix matrix_a(2, 2);
  S21Matrix matrix_b(2, 2);
  S21Matrix result(2, 2);

  matrix_a(0, 0) = 3;
  matrix_a(0, 1) = 2;
  matrix_a(1, 0) = -6;
  matrix_a(1, 1) = 0;

  matrix_b(0, 0) = -7;
  matrix_b(0, 1) = 0;
  matrix_b(1, 0) = -3;
  matrix_b(1, 1) = 2;

  result(0, 0) = 10;
  result(0, 1) = 2;
  result(1, 0) = -3;
  result(1, 1) = -2;

  ASSERT_TRUE((matrix_a - matrix_b) == result);
}
TEST(OperatorMultiplyMatrix, True) {
  S21Matrix matrix_a(2, 2);
  S21Matrix matrix_b(2, 2);
  S21Matrix result(2, 2);

  matrix_a(0, 0) = 3;
  matrix_a(0, 1) = 2;
  matrix_a(1, 0) = -6.6;
  matrix_a(1, 1) = 0;

  matrix_b(0, 0) = -7;
  matrix_b(0, 1) = 0;
  matrix_b(1, 0) = -3.5;
  matrix_b(1, 1) = 2;

  result(0, 0) = -28;
  result(0, 1) = 4;
  result(1, 0) = 46.2;
  result(1, 1) = 0;

  ASSERT_TRUE((matrix_a * matrix_b) == result);
}
TEST(OperatorMultiplyNumber, True) {
  S21Matrix matrix_a(2, 2);
  S21Matrix result(2, 2);

  matrix_a(0, 0) = -78.14;
  matrix_a(0, 1) = 0;
  matrix_a(1, 0) = -0.3;
  matrix_a(1, 1) = 2;

  result(0, 0) = -781.4;
  result(0, 1) = 0;
  result(1, 0) = -3;
  result(1, 1) = 20;

  ASSERT_TRUE((matrix_a * 10.0) == result);
}
TEST(OperatorMultiplyNum, True) {
  S21Matrix matrix_a(2, 2);
  S21Matrix result(2, 2);

  matrix_a(0, 0) = -78.14;
  matrix_a(0, 1) = 0;
  matrix_a(1, 0) = -0.3;
  matrix_a(1, 1) = 2;

  result(0, 0) = -781.4;
  result(0, 1) = 0;
  result(1, 0) = -3;
  result(1, 1) = 20;

  ASSERT_TRUE((matrix_a * 10.0) == result);
}
TEST(OperatorEqMatrix, True) {
  S21Matrix matrix_a(3, 3);
  S21Matrix matrix_b(1, 3);
  matrix_a = matrix_b;
  ASSERT_TRUE(matrix_a == matrix_b);
}

TEST(OperatorSumMatrix, True) {
  S21Matrix matrix_a(2, 2);
  S21Matrix matrix_b(2, 2);
  S21Matrix result(2, 2);

  matrix_a(0, 0) = 3.14;
  matrix_a(0, 1) = 0.56;
  matrix_a(1, 0) = -69.3;
  matrix_a(1, 1) = 0;

  matrix_b(0, 0) = -78.14;
  matrix_b(0, 1) = 0;
  matrix_b(1, 0) = -0.3;
  matrix_b(1, 1) = 2;

  result(0, 0) = -75;
  result(0, 1) = 0.56;
  result(1, 0) = -69.6;
  result(1, 1) = 2;

  ASSERT_TRUE((matrix_a += matrix_b) == result);
}
TEST(OperatorSumMatrix, False) {
  S21Matrix matrix_a(1, 2);
  S21Matrix matrix_b(2, 2);

  matrix_a(0, 0) = 3.14;
  matrix_a(0, 1) = 0.56;

  matrix_b(0, 0) = -78.14;
  matrix_b(0, 1) = 0;
  matrix_b(1, 0) = -0.3;
  matrix_b(1, 1) = 2;

  EXPECT_THROW((matrix_a += matrix_b), std::logic_error);
}

TEST(OperatorSubMatrix, True) {
  S21Matrix matrix_a(2, 2);
  S21Matrix matrix_b(2, 2);
  S21Matrix result(2, 2);

  matrix_a(0, 0) = 3;
  matrix_a(0, 1) = 2;
  matrix_a(1, 0) = -6;
  matrix_a(1, 1) = 0;

  matrix_b(0, 0) = -7;
  matrix_b(0, 1) = 0;
  matrix_b(1, 0) = -3;
  matrix_b(1, 1) = 2;

  result(0, 0) = 10;
  result(0, 1) = 2;
  result(1, 0) = -3;
  result(1, 1) = -2;

  ASSERT_TRUE((matrix_a -= matrix_b) == result);
}
TEST(OperatorSubMatrix, False) {
  S21Matrix matrix_a(1, 2);
  S21Matrix matrix_b(2, 2);

  matrix_a(0, 0) = 3;
  matrix_a(0, 1) = 2;

  matrix_b(0, 0) = -7;
  matrix_b(0, 1) = 0;
  matrix_b(1, 0) = -3;
  matrix_b(1, 1) = 2;

  EXPECT_THROW(matrix_a -= matrix_b, std::logic_error);
}
TEST(OperatorMulNumber, True) {
  S21Matrix matrix_a(2, 2);
  S21Matrix result(2, 2);

  matrix_a(0, 0) = -78.14;
  matrix_a(0, 1) = 0;
  matrix_a(1, 0) = -0.3;
  matrix_a(1, 1) = 2;

  result(0, 0) = -781.4;
  result(0, 1) = 0;
  result(1, 0) = -3;
  result(1, 1) = 20;

  ASSERT_TRUE((matrix_a *= 10) == result);
}
TEST(OperatorMulMatrix, True) {
  S21Matrix matrix_a(2, 2);
  S21Matrix matrix_b(2, 2);
  S21Matrix result(2, 2);

  matrix_a(0, 0) = 3;
  matrix_a(0, 1) = 2;
  matrix_a(1, 0) = -6.6;
  matrix_a(1, 1) = 0;

  matrix_b(0, 0) = -7;
  matrix_b(0, 1) = 0;
  matrix_b(1, 0) = -3.5;
  matrix_b(1, 1) = 2;

  result(0, 0) = -28.0;
  result(0, 1) = 4.0;
  result(1, 0) = 46.2;
  result(1, 1) = 0.0;

  ASSERT_TRUE((matrix_a *= matrix_b) == result);
}
TEST(OperatorMulMatrix, False) {
  S21Matrix matrix_a(2, 1);
  S21Matrix matrix_b(2, 2);

  matrix_a(0, 0) = 3;
  matrix_a(1, 0) = -6.6;

  matrix_b(0, 0) = -7;
  matrix_b(0, 1) = 0;
  matrix_b(1, 0) = -3.5;
  matrix_b(1, 1) = 2;

  EXPECT_THROW(matrix_a *= matrix_b, std::logic_error);
}

TEST(MulMatrixOperationTest1, True) {
  S21Matrix matrix_1 = S21Matrix(2, 2);
  matrix_1(0, 0) = 1.0;
  matrix_1(0, 1) = 2.0;
  matrix_1(1, 0) = 3.0;
  matrix_1(1, 1) = 4.0;

  S21Matrix matrix_2 = S21Matrix(2, 2);
  matrix_2(0, 0) = 4.0;
  matrix_2(0, 1) = 5.6;
  matrix_2(1, 0) = 6.0;
  matrix_2(1, 1) = 10.0;

  S21Matrix result = S21Matrix(2, 2);
  result(0, 0) = 16.0;
  result(0, 1) = 25.6;
  result(1, 0) = 36.0;
  result(1, 1) = 56.8;
  matrix_1.MulMatrix(matrix_2);

  EXPECT_TRUE(matrix_1 == result);
}

TEST(MulMatrixOperationTest2, True) {
  S21Matrix matrix_1 = S21Matrix(3, 2);
  matrix_1(0, 0) = 1.123456;
  matrix_1(0, 1) = 4.543564;
  matrix_1(1, 0) = 2.546356;
  matrix_1(1, 1) = 5.454325;
  matrix_1(2, 0) = 3.254562;
  matrix_1(2, 1) = 6.252452;

  S21Matrix matrix_2 = S21Matrix(2, 3);
  matrix_2(0, 0) = 1.254352;
  matrix_2(0, 1) = 6.245223;
  matrix_2(0, 2) = -1.14245;
  matrix_2(1, 0) = 2.153252;
  matrix_2(1, 1) = 3.132411;
  matrix_2(1, 2) = 4.413214;

  S21Matrix result = S21Matrix(3, 3);
  result(0, 0) = 11.19264755064;
  result(0, 1) = 21.248543103492;
  result(0, 2) = 18.768227947496;
  result(1, 0) = 14.938562956212;
  result(1, 1) = 32.987748684963;
  result(1, 2) = 21.16201903835;
  result(2, 0) = 17.545471127728;
  result(2, 1) = 39.910714879098;
  result(2, 2) = 23.875234343828;

  matrix_1.MulMatrix(matrix_2);
  EXPECT_TRUE(matrix_1 == result);
}

TEST(OperatorMulMatrix1, True) {
  S21Matrix matrix_1 = S21Matrix(2, 2);
  matrix_1(0, 0) = 1.0;
  matrix_1(0, 1) = 2.0;
  matrix_1(1, 0) = 3.0;
  matrix_1(1, 1) = 4.0;

  S21Matrix matrix_2 = S21Matrix(2, 2);
  matrix_2(0, 0) = 4.0;
  matrix_2(0, 1) = 5.6;
  matrix_2(1, 0) = 6.0;
  matrix_2(1, 1) = 10.0;

  S21Matrix matrix_3 = matrix_1 * matrix_2;

  S21Matrix result = S21Matrix(2, 2);
  result(0, 0) = 16.0;
  result(0, 1) = 25.6;
  result(1, 0) = 36.0;
  result(1, 1) = 56.8;

  EXPECT_TRUE(matrix_3 == result);
}

TEST(OperatorMulMatrix2, True) {
  S21Matrix matrix_1 = S21Matrix(2, 2);
  matrix_1(0, 0) = 1.0;
  matrix_1(0, 1) = 2.0;
  matrix_1(1, 0) = 3.0;
  matrix_1(1, 1) = 4.0;

  S21Matrix matrix_2 = S21Matrix(2, 2);
  matrix_2(0, 0) = 4.0;
  matrix_2(0, 1) = 5.6;
  matrix_2(1, 0) = 6.0;
  matrix_2(1, 1) = 10.0;

  matrix_1 *= matrix_2;

  S21Matrix result = S21Matrix(2, 2);
  result(0, 0) = 16.0;
  result(0, 1) = 25.6;
  result(1, 0) = 36.0;
  result(1, 1) = 56.8;

  EXPECT_TRUE(matrix_1 == result);
}

TEST(Transpose, True) {
  S21Matrix matrix_a(4, 3);
  S21Matrix result(3, 4);

  matrix_a(0, 0) = 7;
  matrix_a(0, 1) = -98;
  matrix_a(0, 2) = 0.5;
  matrix_a(1, 0) = 0;
  matrix_a(1, 1) = 5.4;
  matrix_a(1, 2) = 32;
  matrix_a(2, 0) = 3.12;
  matrix_a(2, 1) = 2323;
  matrix_a(2, 2) = 23;
  matrix_a(3, 0) = -78;
  matrix_a(3, 1) = 476.4;
  matrix_a(3, 2) = 21;

  result(0, 0) = 7;
  result(0, 1) = 0;
  result(0, 2) = 3.12;
  result(0, 3) = -78;
  result(1, 0) = -98;
  result(1, 1) = 5.4;
  result(1, 2) = 2323;
  result(1, 3) = 476.4;
  result(2, 0) = 0.5;
  result(2, 1) = 32;
  result(2, 2) = 23;
  result(2, 3) = 21;
  S21Matrix res = matrix_a.Transpose();
  ASSERT_TRUE(res == result);
}

TEST(Determinant1, True) {
  S21Matrix matrix(1, 1);
  matrix(0, 0) = 234.2312;
  EXPECT_DOUBLE_EQ(matrix.Determinant(), 234.2312);
}

TEST(Determinant2, True) {
  S21Matrix matrix(2, 2);
  matrix(0, 0) = 1.23;
  matrix(0, 1) = 2.83;
  matrix(1, 0) = 3.63;
  matrix(1, 1) = 4.33;
  EXPECT_DOUBLE_EQ(matrix.Determinant(), -4.947);
}

TEST(Determinant3, True) {
  S21Matrix matrix(3, 3);
  matrix(0, 0) = 1.0;
  matrix(0, 1) = 32.12;
  matrix(0, 2) = 432.12;
  matrix(1, 0) = 2.32;
  matrix(1, 1) = 123.12;
  matrix(1, 2) = 0.453;
  matrix(2, 0) = 3.34;
  matrix(2, 1) = 42.12;
  matrix(2, 2) = 345.21;
  EXPECT_DOUBLE_EQ(matrix.Determinant(), -118663.3809096);
}

TEST(Determinant4, True) {
  S21Matrix matrix(4, 4);
  matrix(0, 0) = 1.0;
  matrix(0, 1) = 5.2;
  matrix(0, 2) = 2.6;
  matrix(0, 3) = 7.2;
  matrix(1, 0) = 1.4;
  matrix(1, 1) = 5.2;
  matrix(1, 2) = 7.3;
  matrix(1, 3) = 3.1;
  matrix(2, 0) = 6.2;
  matrix(2, 1) = 5.2;
  matrix(2, 2) = 6.2;
  matrix(2, 3) = 2.2;
  matrix(3, 0) = 2.2;
  matrix(3, 1) = 6.2;
  matrix(3, 2) = 0.1;
  matrix(3, 3) = 5.1;
  EXPECT_DOUBLE_EQ(matrix.Determinant(), -672.3544);
}

TEST(Inverse1, False) {
  S21Matrix matrix_a(4, 3);
  S21Matrix result(3, 4);

  matrix_a(0, 0) = 7;
  matrix_a(0, 1) = -98;
  matrix_a(0, 2) = 0.5;
  matrix_a(1, 0) = 0;
  matrix_a(1, 1) = 5.4;
  matrix_a(1, 2) = 32;
  matrix_a(2, 0) = 3.12;
  matrix_a(2, 1) = 2323;
  matrix_a(2, 2) = 23;
  matrix_a(3, 0) = -78;
  matrix_a(3, 1) = 476.4;
  matrix_a(3, 2) = 21;

  EXPECT_THROW(matrix_a.InverseMatrix(), std::logic_error);
}
TEST(Inverse2, True) {
  S21Matrix matrix_a(3, 3);
  S21Matrix result(3, 3);

  matrix_a(0, 0) = 2;
  matrix_a(0, 1) = 5;
  matrix_a(0, 2) = 7;
  matrix_a(1, 0) = 6;
  matrix_a(1, 1) = 3;
  matrix_a(1, 2) = 4;
  matrix_a(2, 0) = 5;
  matrix_a(2, 1) = -2;
  matrix_a(2, 2) = -3;

  result(0, 0) = 1;
  result(0, 1) = -1;
  result(0, 2) = 1;
  result(1, 0) = -38;
  result(1, 1) = 41;
  result(1, 2) = -34;
  result(2, 0) = 27;
  result(2, 1) = -29;
  result(2, 2) = 24;

  ASSERT_TRUE(matrix_a.InverseMatrix().EqMatrix(result));

  S21Matrix matrix_b(3, 3);
  matrix_b(0, 0) = 1;
  matrix_b(0, 1) = 2;
  matrix_b(0, 2) = 3;
  matrix_b(1, 0) = 4;
  matrix_b(1, 1) = 5;
  matrix_b(1, 2) = 6;
  matrix_b(2, 0) = 7;
  matrix_b(2, 1) = 8;
  matrix_b(2, 2) = 9;

  EXPECT_THROW(matrix_b.InverseMatrix(), std::logic_error);
}

TEST(Inverse3, True) {
  S21Matrix matrix_1(1, 1);
  matrix_1(0, 0) = 432.1;

  S21Matrix matrix_2 = matrix_1.InverseMatrix();
  S21Matrix result(1, 1);
  result(0, 0) = 0.002314279102;

  std::cout << matrix_2(0, 0) << std::endl;

  EXPECT_TRUE(result == matrix_2);
}

TEST(Inverse4, True) {
  S21Matrix matrix_1(3, 3);
  matrix_1(0, 0) = 2.0;
  matrix_1(0, 1) = 5.0;
  matrix_1(0, 2) = 7.0;
  matrix_1(1, 0) = 6.0;
  matrix_1(1, 1) = 3.0;
  matrix_1(1, 2) = 4.0;
  matrix_1(2, 0) = 5.0;
  matrix_1(2, 1) = -2.0;
  matrix_1(2, 2) = -3.0;

  S21Matrix matrix_2 = matrix_1.InverseMatrix();
  S21Matrix result(3, 3);
  result(0, 0) = 1.0;
  result(0, 1) = -1.0;
  result(0, 2) = 1.0;
  result(1, 0) = -38.0;
  result(1, 1) = 41.0;
  result(1, 2) = -34.0;
  result(2, 0) = 27.0;
  result(2, 1) = -29.0;
  result(2, 2) = 24.0;

  EXPECT_TRUE(result == matrix_2);
}

TEST(Get, True) {
  S21Matrix matrix_a(3, 3);

  matrix_a(0, 0) = 2;
  matrix_a(0, 1) = 5;
  matrix_a(0, 2) = 7;
  matrix_a(1, 0) = 6;
  matrix_a(1, 1) = 3;
  matrix_a(1, 2) = 4;
  matrix_a(2, 0) = 5;
  matrix_a(2, 1) = -2;
  matrix_a(2, 2) = -3;

  ASSERT_EQ(matrix_a.GetRows(), 3);
  ASSERT_EQ(matrix_a.GetCols(), 3);
}
TEST(Set, True) {
  S21Matrix matrix_a(3, 3);
  S21Matrix result(3, 2);

  matrix_a(0, 0) = 2;
  matrix_a(0, 1) = 5;
  matrix_a(0, 2) = 7;
  matrix_a(1, 0) = 6;
  matrix_a(1, 1) = 3;
  matrix_a(1, 2) = 4;
  matrix_a(2, 0) = 5;
  matrix_a(2, 1) = -2;
  matrix_a(2, 2) = -3;

  result(0, 0) = 2;
  result(0, 1) = 5;

  result(1, 0) = 6;
  result(1, 1) = 3;

  result(2, 0) = 5;
  result(2, 1) = -2;
  matrix_a.SetCols(2);

  ASSERT_TRUE(matrix_a == result);

  S21Matrix matrix_b(3, 3);
  S21Matrix result_b(2, 3);

  matrix_b(0, 0) = 2;
  matrix_b(0, 1) = 5;
  matrix_b(0, 2) = 7;
  matrix_b(1, 0) = 6;
  matrix_b(1, 1) = 3;
  matrix_b(1, 2) = 4;
  matrix_b(2, 0) = 5;
  matrix_b(2, 1) = -2;
  matrix_b(2, 2) = -3;

  result_b(0, 0) = 2;
  result_b(0, 1) = 5;
  result_b(0, 2) = 7;
  result_b(1, 0) = 6;
  result_b(1, 1) = 3;
  result_b(1, 2) = 4;

  matrix_b.SetRows(2);
  ASSERT_TRUE(matrix_b == result_b);
}

TEST(CalcComplements, True) {
  S21Matrix matrix_1(3, 3);
  matrix_1(0, 0) = 1.0;
  matrix_1(0, 1) = 2.0;
  matrix_1(0, 2) = 3.0;
  matrix_1(1, 0) = 0.0;
  matrix_1(1, 1) = 4.0;
  matrix_1(1, 2) = 2.0;
  matrix_1(2, 0) = 5.0;
  matrix_1(2, 1) = 2.0;
  matrix_1(2, 2) = 1.0;
  S21Matrix matrix_2 = matrix_1.CalcComplements();

  S21Matrix result(3, 3);
  result(0, 0) = 0.0;
  result(0, 1) = 10.0;
  result(0, 2) = -20.0;
  result(1, 0) = 4.0;
  result(1, 1) = -14.0;
  result(1, 2) = 8.0;
  result(2, 0) = -8.0;
  result(2, 1) = -2.0;
  result(2, 2) = 4.0;

  EXPECT_TRUE(result == matrix_2);
}

TEST(MulNumberOperationTest, True) {
  S21Matrix matrix = S21Matrix(2, 2);
  matrix(0, 0) = 1.0;
  matrix(0, 1) = 2.0;
  matrix(1, 0) = 3.0;
  matrix(1, 1) = 4.0;
  matrix.MulNumber(12.4);

  S21Matrix result = S21Matrix(2, 2);
  result(0, 0) = 12.4;
  result(0, 1) = 24.8;
  result(1, 0) = 37.2;
  result(1, 1) = 49.6;

  EXPECT_TRUE(matrix == result);
}

TEST(MulMatrixNumberOperatorTest, True) {
  S21Matrix matrix_1 = S21Matrix(2, 2);
  matrix_1(0, 0) = 1.0;
  matrix_1(0, 1) = 2.0;
  matrix_1(1, 0) = 3.0;
  matrix_1(1, 1) = 4.0;
  S21Matrix matrix_2 = matrix_1 * 12.4;

  S21Matrix result = S21Matrix(2, 2);
  result(0, 0) = 12.4;
  result(0, 1) = 24.8;
  result(1, 0) = 37.2;
  result(1, 1) = 49.6;

  EXPECT_TRUE(matrix_2 == result);
}

TEST(MulMatrixNumberAssignOperatorTest, True) {
  S21Matrix matrix = S21Matrix(2, 2);
  matrix(0, 0) = 1.0;
  matrix(0, 1) = 2.0;
  matrix(1, 0) = 3.0;
  matrix(1, 1) = 4.0;
  matrix *= 12.4;

  S21Matrix result = S21Matrix(2, 2);
  result(0, 0) = 12.4;
  result(0, 1) = 24.8;
  result(1, 0) = 37.2;
  result(1, 1) = 49.6;

  EXPECT_TRUE(matrix == result);
}

TEST(MulNumberMatrixOperatorTest, True) {
  S21Matrix matrix_1 = S21Matrix(2, 2);
  matrix_1(0, 0) = 1.0;
  matrix_1(0, 1) = 2.0;
  matrix_1(1, 0) = 3.0;
  matrix_1(1, 1) = 4.0;
  S21Matrix matrix_2 = 12.4 * matrix_1;

  S21Matrix result = S21Matrix(2, 2);
  result(0, 0) = 12.4;
  result(0, 1) = 24.8;
  result(1, 0) = 37.2;
  result(1, 1) = 49.6;

  EXPECT_TRUE(matrix_2 == result);
}

TEST(Constructor, Default) {
  S21Matrix test = S21Matrix();
  EXPECT_TRUE(test.GetRows() == 0);
  EXPECT_TRUE(test.GetCols() == 0);
}

TEST(Constructor, By2Args) {
  S21Matrix test = S21Matrix(3, 3);
  EXPECT_TRUE(test.GetRows() == 3);
  EXPECT_TRUE(test.GetCols() == 3);
}

TEST(Constructor, Copy) {
  S21Matrix test = S21Matrix(3, 3);
  test(0, 0) = 1;
  S21Matrix temp = S21Matrix(test);
  EXPECT_TRUE(test == temp);
}

TEST(ConstructorTest, MoveConstructor) {
  S21Matrix matrix_1 = S21Matrix(2, 2);
  matrix_1(0, 0) = 1.0;
  matrix_1(0, 1) = 2.0;
  matrix_1(1, 0) = 3.0;
  matrix_1(1, 1) = 4.0;
  S21Matrix matrix_2 = matrix_1 * 12.4;

  matrix_1 = S21Matrix(std::move(matrix_2));
  S21Matrix result = S21Matrix(2, 2);
  result(0, 0) = 12.4;
  result(0, 1) = 24.8;
  result(1, 0) = 37.2;
  result(1, 1) = 49.6;

  EXPECT_TRUE(matrix_1 == result);
  EXPECT_EQ(matrix_2.GetCols(), 0);
  EXPECT_EQ(matrix_2.GetRows(), 0);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}