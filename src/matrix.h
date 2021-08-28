#ifndef DEF_matrix
#define DEF_matrix

#include <vector>
#include <functional>

class Matrix
{

private:
    
    std::vector<double> matrix;
    static void f_random(double &x);
    unsigned int nCols;
    unsigned int nRows;

public:

    Matrix(const size_t cols_num, const size_t row_num);
    Matrix(Matrix const &copy);
    Matrix &operator=(Matrix const &copy);
    Matrix();
    //Matrix& Matrix::operator=(Matrix const& copy);
    ~Matrix();

    size_t getCols();
    size_t getRows();
    void Reset(const size_t cols_num, const size_t row_num);

    void print() const;
    void set(double val, unsigned int x, unsigned int y);
    double get(unsigned int x, unsigned int y);

    void init(double const &valueOfInit);

    void T();

    static Matrix add(Matrix const &m1, Matrix const &m2);

    static Matrix multiply(Matrix const &m1, Matrix const &m2);

    // static bool operator==(Matrix const& m1, Matrix const& m2);
    // static Matrix& operator+(Matrix const& m1, Matrix const& m2);
    void nRand(double local, double scale);
    void randomize();

    void apply(std::function<void(double &)> func);

    void add_by(double const &adding_Num);

    void multiply_by(double const &multiplying_Num);

    static Matrix normal_multiplication(Matrix const &m1, Matrix const &m2);

    double *to_ptr(bool arrangement_type);
    void fromArray(double arr[], size_t size);
    void fromVector(std::vector<double> arr);
    std::vector<double> toVector();
    friend bool operator==(Matrix const &m1, Matrix const &m2);
};

std::pair<double, double> Random_normal_disturbution(double local, double scale);
std::pair<double, double> normal_distrbution(double x, double y, double mean, double scale);

bool operator==(Matrix const &m1, Matrix const &m2);
bool operator!=(Matrix const &m1, Matrix const &m2);

Matrix operator+(Matrix const &m1, Matrix const &m2);
Matrix operator*(Matrix const &m1, Matrix const &m2);

#endif