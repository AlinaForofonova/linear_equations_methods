#pragma once

#define myType float

bool qrFactorization(myType** matrix, myType** T, myType* x, int n);

//����������������
void transpose_QR(myType** matrix, int n);
//��������� ������ ����� �� ������� ��������, ������ ����� ����������������
void matrixMultB_QR(myType** T, myType** matrix, int n);

//�������� ���//
void backward_QR(myType** matrix, myType* x, int n);
