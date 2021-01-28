#include <iostream>
#include <cmath>

#define myType float

using namespace std;


//���������� ������� � ������������������ ����//
//���������� true ���� ��� �������
bool upTriangleMatrix_Gauss(myType** matrix, int n);
//��������� ����� �������� ��������, n-k = ���������� ��������������� �����
//���������� true ���� � ������� ��� ���� �� ���� �� ������� �������
bool selectMaxElement_Gauss(myType** matrix, int n, int k);
//�������� ����� �������
void swapLine_Gauss(myType** matrix, int k, int index);

//�������� ���//
void backward_Gauss(myType** matrix, myType* x, int n);
//����� ��� �������� ����
myType sumBackward_Gauss(myType** matrix, myType* x, int n, int i);


//������� ������� ������
bool gauss(myType** matrix, myType* x, int n)
{
    if(upTriangleMatrix_Gauss(matrix, n))
        backward_Gauss(matrix, x, n);
    else
        return false;

    return true;
}

//���������� ������� � ������������������ ����
bool upTriangleMatrix_Gauss(myType** matrix, int n)
{
    //����������� c[i][k] = a[i][k]/a[k][k] (�� ���������, 7 ��������)
    myType c;

    //� ������ ����� k ��������� � ����������� n-k = ���������� ��������������� �����
    for(int k = 1; k < n+1; k++)
    {
        if(selectMaxElement_Gauss(matrix, n, k))
            for(int i = k; i < n; i++)
            {
                c = matrix[i][k-1]/matrix[k-1][k-1];
                matrix[i][k-1] = 0;
                for(int j = k; j < n+1; j++)
                    matrix[i][j] -= c * matrix[k-1][j];
            }
        else
            return false;
    }

    return true;
}

//����� �������� ��������
//n - �����������, k - �������
bool selectMaxElement_Gauss(myType** matrix, int n, int k)
{
    myType max = matrix[k-1][k-1];
    int index = k-1;

    for(int i = k; i < n; i++)
        if(fabs(matrix[i][k-1]) > fabs(max)){
            max = matrix[i][k-1];
            index = i;
        }

    if(fabs(max) < 1e-15)
        return false;

    swapLine_Gauss(matrix, k, index);

    return true;
}

//�������� ����� �������
void swapLine_Gauss(myType** matrix, int k, int index)
{
    myType* adress = matrix[index];
    matrix[index]  = matrix[k-1];
    matrix[k-1]    = adress;
}

//�������� ���
void backward_Gauss(myType** matrix, myType* x, int n)
{
    for(int i = n-1; i > -1; i--)
        x[i] = (matrix[i][n] - sumBackward_Gauss(matrix, x, n, i))/matrix[i][i];
}

//����� ��� �������� ����
myType sumBackward_Gauss(myType** matrix, myType* x, int n, int i)
{
    myType sum = 0;

    for(int j = i+1; j < n; j++)
        sum += x[j] * matrix[i][j];

    return sum;
}
