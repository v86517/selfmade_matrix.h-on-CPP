#include <gtest/gtest.h>

#include "s21_matrix_oop.h"

TEST(Inverse, T1) {
  S21Matrix mat(1, 1);
  S21Matrix exp(1, 1);
  S21Matrix res(1, 1);
  mat(0, 0) = 2;
  exp(0, 0) = 0.5;

  res = mat.InverseMatrix();

  ASSERT_TRUE(res == exp);
}

TEST(Inverse, T2) {
  S21Matrix mat(3, 3);
  S21Matrix exp(3, 3);
  S21Matrix res(3, 3);
  mat(0, 0) = 2;
  mat(0, 1) = 5;
  mat(0, 2) = 7;
  mat(1, 0) = 6;
  mat(1, 1) = 3;
  mat(1, 2) = 4;
  mat(2, 0) = 5;
  mat(2, 1) = -2;
  mat(2, 2) = -3;
  exp(0, 0) = 1;
  exp(0, 1) = -1;
  exp(0, 2) = 1;
  exp(1, 0) = -38;
  exp(1, 1) = 41;
  exp(1, 2) = -34;
  exp(2, 0) = 27;
  exp(2, 1) = -29;
  exp(2, 2) = 24;
  res = mat.InverseMatrix();

  ASSERT_TRUE(res == exp);
}

TEST(Inverse, Throw) {
  S21Matrix mat(3, 3);
  mat(0, 0) = -3;
  mat(0, 1) = -1;
  mat(0, 2) = -1;
  mat(1, 0) = 3;
  mat(1, 1) = 1;
  mat(1, 2) = 1;
  mat(2, 0) = 6;
  mat(2, 1) = 2;
  mat(2, 2) = 1;

  EXPECT_THROW(mat.InverseMatrix(), std::out_of_range);
}
TEST(Determinant, Throw) {
  S21Matrix mat(2, 3);
  EXPECT_THROW(mat.Determinant(), std::out_of_range);
}

TEST(CalcComplements, Throw) {
  S21Matrix mat(2, 3);
  EXPECT_THROW(mat.CalcComplements(), std::out_of_range);
}

TEST(MulMatrix, T1) {
  S21Matrix mat1(2, 2);
  S21Matrix mat2(2, 2);
  S21Matrix res(2, 2);
  S21Matrix exp(2, 2);
  mat1(0, 0) = 1;
  mat1(0, 1) = 2;
  mat1(1, 0) = 3;
  mat1(1, 1) = 4;
  mat2(0, 0) = 5;
  mat2(0, 1) = 6;
  mat2(1, 0) = 7;
  mat2(1, 1) = 8;
  exp(0, 0) = 19;
  exp(0, 1) = 22;
  exp(1, 0) = 43;
  exp(1, 1) = 50;

  res = mat1 * mat2;

  ASSERT_TRUE(res == exp);
}

TEST(MulMatrix, T2) {
  S21Matrix mat1(2, 2);
  S21Matrix mat2(2, 2);
  S21Matrix exp(2, 2);
  mat1(0, 0) = 1;
  mat1(0, 1) = 2;
  mat1(1, 0) = 3;
  mat1(1, 1) = 4;
  mat2(0, 0) = 5;
  mat2(0, 1) = 6;
  mat2(1, 0) = 7;
  mat2(1, 1) = 8;
  exp(0, 0) = 19;
  exp(0, 1) = 22;
  exp(1, 0) = 43;
  exp(1, 1) = 50;

  mat1 *= mat2;

  ASSERT_TRUE(mat1 == exp);
}

TEST(MulNumber, T1) {
  S21Matrix mat(2, 2);
  double num = 2;
  S21Matrix res(2, 2);
  S21Matrix exp(2, 2);
  mat(0, 0) = 14;
  mat(0, 1) = 1;
  mat(1, 0) = -12.2;
  mat(1, 1) = -9;
  exp(0, 0) = 28;
  exp(0, 1) = 2;
  exp(1, 0) = -24.4;
  exp(1, 1) = -18;

  mat *= num;

  ASSERT_TRUE(mat == exp);
}

TEST(MulNumber, T2) {
  S21Matrix mat(2, 2);
  double num = 2;
  S21Matrix res(2, 2);
  S21Matrix exp(2, 2);
  mat(0, 0) = 2;
  mat(0, 1) = 2.5;
  mat(1, 0) = 3;
  mat(1, 1) = -3.5;
  exp(0, 0) = 4;
  exp(0, 1) = 5;
  exp(1, 0) = 6;
  exp(1, 1) = -7;

  res = mat * num;

  ASSERT_TRUE(res == exp);
}

TEST(MulMatrix, Throw) {
  S21Matrix mat1(2, 3);
  S21Matrix mat2(4, 5);
  EXPECT_THROW(mat1 * mat2, std::out_of_range);
}

TEST(SumMatrix, T1) {
  S21Matrix mat1(2, 2);
  S21Matrix mat2(2, 2);
  S21Matrix res(2, 2);
  S21Matrix exp(2, 2);
  mat1(0, 0) = 4;
  mat1(0, 1) = 1;
  mat1(1, 0) = 8;
  mat1(1, 1) = 2;
  mat2(0, 0) = 3;
  mat2(0, 1) = 9;
  mat2(1, 0) = 2;
  mat2(1, 1) = 3;
  exp(0, 0) = 7;
  exp(0, 1) = 10;
  exp(1, 0) = 10;
  exp(1, 1) = 5;

  res = mat1 + mat2;
  ASSERT_TRUE(res == exp);
}

TEST(SumMatrix, T2) {
  S21Matrix mat1(2, 2);
  S21Matrix mat2(2, 2);
  S21Matrix exp(2, 2);
  mat1(0, 0) = 4;
  mat1(0, 1) = 1;
  mat1(1, 0) = 8;
  mat1(1, 1) = 2;
  mat2(0, 0) = 3;
  mat2(0, 1) = 9;
  mat2(1, 0) = 2;
  mat2(1, 1) = 3;
  exp(0, 0) = 7;
  exp(0, 1) = 10;
  exp(1, 0) = 10;
  exp(1, 1) = 5;

  mat1 += mat2;
  ASSERT_TRUE(mat1 == exp);
}

TEST(SumMatrix, Throw) {
  S21Matrix mat1(2, 3);
  S21Matrix mat2(4, 5);
  EXPECT_THROW(mat1 + mat2, std::out_of_range);
}

TEST(SubMatrix, T1) {
  S21Matrix mat1(2, 2);
  S21Matrix mat2(2, 2);
  S21Matrix res(2, 2);
  S21Matrix exp(2, 2);
  mat1(0, 0) = 14;
  mat1(0, 1) = 1;
  mat1(1, 0) = -12.2;
  mat1(1, 1) = -9;
  mat2(0, 0) = 14;
  mat2(0, 1) = 1;
  mat2(1, 0) = -12.2;
  mat2(1, 1) = -9;
  exp(0, 0) = 0;
  exp(0, 1) = 0;
  exp(1, 0) = 0;
  exp(1, 1) = 0;

  res = mat1 - mat2;
  ASSERT_TRUE(res == exp);
}

TEST(SubMatrix, T2) {
  S21Matrix mat1(2, 2);
  S21Matrix mat2(2, 2);
  S21Matrix exp(2, 2);
  mat1(0, 0) = 14;
  mat1(0, 1) = 1;
  mat1(1, 0) = -12.2;
  mat1(1, 1) = -9;
  mat2(0, 0) = 14;
  mat2(0, 1) = 1;
  mat2(1, 0) = -12.2;
  mat2(1, 1) = -9;
  exp(0, 0) = 0;
  exp(0, 1) = 0;
  exp(1, 0) = 0;
  exp(1, 1) = 0;

  mat1 -= mat2;
  ASSERT_TRUE(mat1 == exp);
}

TEST(SubMatrix, Throw) {
  S21Matrix mat1(2, 3);
  S21Matrix mat2(4, 5);
  EXPECT_THROW(mat1 - mat2, std::out_of_range);
}

TEST(EqMatrix, T1) {
  S21Matrix matrix_a(3, 3);
  S21Matrix matrix_b(2, 2);
  EXPECT_FALSE(matrix_a == matrix_b);
}

TEST(EqMatrix2, T2) {
  S21Matrix matrix_a(3, 3);
  S21Matrix matrix_b(3, 3);
  matrix_b(0, 0) = 0.01;
  EXPECT_FALSE(matrix_a == matrix_b);
}

TEST(constructor, T1) {
  S21Matrix mat1;
  S21Matrix mat2(2, 2);
  S21Matrix sec = std::move(mat2);
  S21Matrix mat3(sec);
  ASSERT_TRUE(mat3 == sec);
}

TEST(constructor, Throw) {
  EXPECT_THROW(S21Matrix mat(1, -1), std::out_of_range);
}
TEST(GetColsRows, T1) {
  S21Matrix matrix_a(3, 3);
  EXPECT_EQ(matrix_a.GetRows(), 3);
  EXPECT_EQ(matrix_a.GetCols(), 3);
  EXPECT_THROW(matrix_a.SetRows(-1), std::out_of_range);
  EXPECT_THROW(matrix_a.SetCols(-1), std::out_of_range);
}

TEST(GetColsRows, T2) {
  S21Matrix matrix_a(3, 3);
  matrix_a.SetRows(5);
  matrix_a.SetCols(2);
  EXPECT_EQ(matrix_a.GetRows(), 5);
  EXPECT_EQ(matrix_a.GetCols(), 2);
}

TEST(OperatorIndex, Throw1) {
  S21Matrix mat(2, 3);
  EXPECT_THROW(mat(4, 5), std::out_of_range);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}