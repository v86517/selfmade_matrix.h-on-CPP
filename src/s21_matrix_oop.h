
#include <iostream>
// #include <stdexcept>

#define EPS 1e-7

class S21Matrix {
 private:
  // attributes
  int rows_, cols_;  // rows and columns attributes
  double** matrix_;  // pointer to the memory where the matrix will be allocated

  // methods
  S21Matrix GetMinorMatrix(int col, int row);

 public:
  // constructors and destructor
  S21Matrix();                        // default constructor
  S21Matrix(int rows, int cols);      // parameterized constructor
  S21Matrix(const S21Matrix& other);  // copy cnstructor
  S21Matrix(S21Matrix&& other);       // move cnstructor
  ~S21Matrix();                       // destructor

  bool EqMatrix(const S21Matrix& other) const;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();

  // accesors and mutators
  void SetRows(int rows);
  int GetRows() const;
  void SetCols(int cols);
  int GetCols() const;

  S21Matrix operator=(const S21Matrix& other);  // assignment operator overload
  double& operator()(int i, int j);             // index operator overload
  double operator()(int i, int j) const;
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix operator+(const S21Matrix& other) const;
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other) const;
  S21Matrix& operator*=(const double num);
  S21Matrix operator*(const double num) const;
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix operator*(const S21Matrix& other) const;
  bool operator==(const S21Matrix& other) const;
};

long double s21_fabs(double x);
