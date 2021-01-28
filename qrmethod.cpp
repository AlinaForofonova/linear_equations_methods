#include <iostream>
#include <cmath>
#include "matrixOperation.h"

#define myType float

using namespace std;


//���������� ������� � ������������������ ���� � ��������� ������� ��������//
//���������� true ���� ��� �������
bool upTriangleMatrix_QR(myType** matrix, myType** T, int n);
//����� �������� ��������
//���������� true ���� � ������� ��� ���� �� ���� �� ������� �������
bool selectMaxElement_QR(myType** matrix, myType** T, int n, int k);
//�������� ����� �������
void swapLine_QR(myType** matrix, int k, int index);
//��������� ���������, ������� B ����������������
void matrixMultMatrix_QR(myType** A, myType** B, int n);
//��������� ������� �������� �� ������ �����
void matrixMultB_QR(myType** T, myType** matrix, int n);
//�������� ��������� �������
void idMatrix_QR(myType** matrix, int n);
//����������������
void transpose_QR(myType** matrix, int n);

//�������� ���//
void backward_QR(myType** matrix, myType* x, int n);
//����� ��� �������� ����
myType sumBackward_QR(myType** matrix, myType* x, int n, int i);


//������� �������
bool qrFactorization(myType** matrix, myType** T, myType* x, int n)
{
    myType* b = new myType [n];
    for(int i = 0; i < n; i++)
        b[i] = matrix[i][n];

    if(upTriangleMatrix_QR(matrix, T, n))
    {
        for(int i = 0; i < n; i++)
            matrix[i][n] = b[i];
        matrixMultB_QR(T, matrix, n);
        backward_QR(matrix, x, n);
        //Q = T^-1 ��� Q = Transpose(T), ��� ��� T - �������������
        transpose_QR(T, n);
        delete [] b;
    }
    else
        return false;

    return true;
}

//���������� ������� � ������������������ ���� � ��������� ������� ��������
bool upTriangleMatrix_QR(myType** matrix, myType** T, int n)
{
    //c - cos, s - sin, hyp - ����������
    myType c, s, hyp;
    //veci/j - ���������� ������ ��� ������������ ������� �� ������� ��������
    myType veci, vecj;

    //�������� ������� �������� � ���������
    idMatrix_QR(T, n);

    for(int i = 0; i < n-1; i++)
        for(int j = i+1; j < n; j++)
        {
            if(selectMaxElement_QR(matrix, T, n, i+1))
                hyp = sqrt(matrix[i][i]*matrix[i][i] + matrix[j][i]*matrix[j][i]);
            else
                return false;

            c = matrix[i][i] / hyp;
            s = matrix[j][i] / hyp;

            for(int k = 0; k < n; k++)
            {
                veci = matrix[i][k];
                vecj = matrix[j][k];
                matrix[i][k] =  c*veci + s*vecj;
                matrix[j][k] = -s*veci + c*vecj;

                veci = T[i][k];
                vecj = T[j][k];
                T[i][k] =  c*veci + s*vecj;
                T[j][k] = -s*veci + c*vecj;
            }

            matrix[j][i] = 0;
        }

    if(fabs(matrix[n-1][n-1]) < 1e-15)
        return false;

    return true;
}

//����� �������� ��������
//n - �����������, k - �������
bool selectMaxElement_QR(myType** matrix, myType** T, int n, int k)
{
    myType max = matrix[k-1][k-1];
    int index = k-1;

    for(int i = k; i < n; i++)
        if(fabs(matrix[i][k-1]) > fabs(max))
        {
            max = matrix[i][k-1];
            index = i;
        }

    if(fabs(max) < 1e-15)
        return false;

    swapLine_QR(matrix, k, index);
    swapLine_QR(T, k, index);

    return true;
}

//�������� ����� �������
void swapLine_QR(myType** matrix, int k, int index)
{
    myType* adress = matrix[index];
    matrix[index]  = matrix[k-1];
    matrix[k-1]    = adress;
}

//��������� ���������, ������� B ����������������
void matrixMultMatrix_QR(myType** A, myType** B, int n)
{
    myType** C = new myType* [n];
        for(int i = 0; i < n; i++)
            C[i] = new myType [n];

    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
        {
            C[i][j] = 0;
            for(int k = 0; k < n; k++)
                C[i][j] += A[i][k] * B[k][j];
        }

    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            B[i][j] = C[i][j];

    for(int i = 0; i < n; i++)
        delete [] C[i];
}

//��������� ������ ����� �� ������� ��������, ������ ����� ����������������
void matrixMultB_QR(myType** T, myType** matrix, int n)
{
    myType* vec = new myType [n];

    for(int i = 0; i < n; i++)
    {
        vec[i] = 0;
        for(int j = 0; j < n; j++)
            vec[i] += T[i][j] * matrix[j][n];
    }

    for(int i = 0; i < n; i++)
        matrix[i][n] = vec[i];

    delete [] vec;
}

//�������� ��������� �������
void idMatrix_QR(myType** matrix, int n)
{
    for(int i = 0; i < n; i++)
        matrix[i][i] = 1;

    for(int i = 0; i < n; i++)
        for(int j = i+1; j < n; j++)
            matrix[i][j] = matrix[j][i] = 0;
}
//����������������
void transpose_QR(myType** matrix, int n)
{
    myType x;

    for(int i = 0; i < n; i++)
        for(int j = i + 1; j < n; j++)
        {
            x = matrix[i][j];
            matrix[i][j] = matrix[j][i];
            matrix[j][i] = x;
        }
}

//�������� ���
void backward_QR(myType** matrix, myType* x, int n)
{
    for(int i = n-1; i > -1; i--)
        x[i] = (matrix[i][n] - sumBackward_QR(matrix, x, n, i))/matrix[i][i];
}

//����� ��� �������� ����
myType sumBackward_QR(myType** matrix, myType* x, int n, int i)
{
    myType sum = 0;

    for(int j = i+1; j < n; j++)
        sum = sum + x[j] * matrix[i][j];

    return sum;
}
