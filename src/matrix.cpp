#include <iostream>
#include <random>
#include <algorithm>
#include <iterator>
#include "matrix.h"

using namespace std;

#define default_Row 1
#define default_Col 1

Matrix::Matrix(const size_t cols_num, const size_t row_num) : nCols(1), nRows(1)
{
    Reset(cols_num, row_num);
}
Matrix::Matrix()
{
}
void Matrix::Reset(const size_t cols_num, const size_t row_num)
{

    const size_t size = cols_num * row_num;
    if (size != 0)
    {
        if (size != (nCols * nRows))
        {
            matrix.resize(size);
        }
        nCols = cols_num;
        nRows = row_num;
    }
    else
    {
        matrix.resize(1);
        nCols = 1;
        nRows = 1;
    }
}

void Matrix::print() const
{

    for (int i = 0; i < nRows; i++)
    {
        cout << "[ ";
        for (int j = 0; j < nCols; j++)
        {
            cout << matrix[i * nCols + j] << " ";
        }
        cout << "]\n";
    }
}

void Matrix::init(double const &valueOfInit = 0)
{
    for (int i = 0; i < matrix.size(); i++)
        matrix[i] = valueOfInit;
}

Matrix::~Matrix()
{
}

void Matrix::T()
{
    if (nCols == nRows)
    {
        //for efficency
        const unsigned int n = nCols;
        double temp;
        for (int i = 0; i < n; i++)
        {
            for (int j = i; j < n; j++)
            {
                temp = matrix[i * n + j];
                matrix[i * n + j] = matrix[j * n + i];
                matrix[j * n + i] = temp;
                //cout << ":" << matrix[j * nCols + i] << " = " << matrix[i * nCols + j] << endl;
            }
        }
    }
    else
    {
        vector<double> newM(nCols * nRows);
        //this->print();
        for (int i = 0; i < nRows; i++)
        {
            for (int j = 0; j < nCols; j++)
            {
                newM[j * nRows + i] = matrix[i * nCols + j];
            }
        }

        swap(nCols, nRows);

        //matrix.clear();

        matrix = newM;
    }
}

Matrix::Matrix(Matrix const &copy) : nCols(1), nRows(1)
{
    Reset(copy.nCols, copy.nRows);
    matrix = copy.matrix;
}
Matrix &Matrix::operator=(Matrix const &copy)
{
    if (this != &copy)
    {
        //later
        nCols = copy.nCols;
        nRows = copy.nRows;
        matrix = copy.matrix;
    }
    return *this;
}

void Matrix::set(double val, unsigned int x, unsigned int y)
{
    if (x < nRows && y < nCols)
        matrix[x * nCols + y] = val;
    else
        cout << "error, indexes should be less than the actual size" << endl;
}
double Matrix::get(unsigned int x, unsigned int y)
{
    return matrix[x * nCols + y];
}

Matrix Matrix::add(Matrix const &m1, Matrix const &m2)
{
    unsigned int row_num = 0, col_num = 0;
    if (m1.nRows != m2.nRows)
    {
        cout << "Both matrix's rows number don't match";
        return m1;
    }
    else if (m1.nCols != m2.nCols)
    {
        cout << "Both matrix's columns number don't match";
        return m2;
    }
    else
    {
        row_num = m1.nRows;
        col_num = m1.nCols;
    }

    Matrix result(col_num, row_num);
    transform(m1.matrix.begin() , m1.matrix.end() , m2.matrix.begin() , result.matrix.begin() , plus<double>() );
    return result;
}

Matrix Matrix::multiply(Matrix const &m1, Matrix const &m2)
{
    unsigned int row_num = m1.nRows, col_num = m2.nCols;
    if (m1.nCols != m2.nRows)
    {
        cout << "Both matrix's don't match";
        Matrix result(col_num, row_num);
        return result;
    }
    Matrix result(col_num, row_num);
    for (int i = 0; i < row_num; i++)
    {
        for (int j = 0; j < col_num; j++)
        {
            result.matrix[i * col_num + j] = 0;
            for (int k = 0; k < m1.nCols; k++)
                result.matrix[i * col_num + j] += m1.matrix[i * m1.nCols + j] * m2.matrix[i * col_num + j];
        }
    }
    return result;
}

void Matrix ::f_random(double &x)
{
    x = 2.0f * ((double)rand() / RAND_MAX) - 1.0f;
}

void Matrix::randomize()
{ //generate numbers randomly between -1 and 1
    apply(f_random);
}

void Matrix::apply(std::function<void(double &)> func)
{
    for_each(matrix.begin(), matrix.end(), func);
}

void Matrix::add_by(double const &adding_Num)
{
    for (int i = 0; i < nRows; i++)
    {
        for (int j = 0; j < nCols; j++)
        {
            matrix[i * nCols + j] += adding_Num;
        }
    }
}

void Matrix::multiply_by(double const &multiplying_Num)
{
    for (int i = 0; i < nRows; i++)
    {
        for (int j = 0; j < nCols; j++)
        {
            matrix[i * nCols + j] *= multiplying_Num;
        }
    }
}

double *Matrix::to_ptr(bool arrangement_type = true)
{
    return &matrix[0];
}



void Matrix::fromArray(double arr[], size_t size = 1)
{
    double *p = arr;

    if (size < 1)
    {
        cout << "Error in size" << endl;
    }

    Reset(1, size);

    copy(arr, arr + size, matrix.begin());
}

void  Matrix::fromVector(vector<double> arr){
    matrix = arr;
}

vector<double>  Matrix::toVector(){
    return matrix;
}

bool operator==(Matrix const &m1, Matrix const &m2)
{
    if ((m1.nRows != m2.nRows) && (m1.nCols != m2.nCols))
        return false;
    return m1.matrix == m2.matrix;
}
bool operator!=(Matrix const &m1, Matrix const &m2)
{
    return !(m1 == m2);
}



Matrix operator+(Matrix const &m1, Matrix const &m2)
{
    return Matrix::add(m1, m2);
}

Matrix operator*(Matrix const &m1, Matrix const &m2)
{
    return Matrix::multiply(m1, m2);
}
