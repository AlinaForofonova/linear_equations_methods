#pragma once

#define myType float

using namespace std;

//������� ������//
//������ ������� c ������ ������
void printMatrixAB(myType** matrix, int n);
//������ ���������� �������
void printMatrix(myType** matrix, int n);
//������ �������
void printVector(myType* x, int n);
//����� ��������� ��� ����������� ������
void degenerateMatrix(myType** matrix, int n);


//����������� ������� � �������
void testSolution(myType** matrix, myType* x, myType* b1, int n);
//��������� ������ (�������� ���������, ��� ������� � �������)
void MatrMult(myType** A, myType** B, int n);
//���������� ������ �����
void perturbationRightPartMax(myType** matrix, int n);
void perturbationRightPartMin(myType** matrix, int n);


//����� �������� � ����������//
//�������������� �����
myType normVector1(myType* vec, int n);
myType normMatrix1(myType** matrix, int n);
//������� �����
myType normVector2(myType* vec, int n);
//���������� �����
myType normVectorInf(myType* vec, int n);
myType normMatrixInf(myType** matrix, int n);


//������ � �������//
//����� ����� �� ������
void choiceFile(string &testName);
//������ �� ����� � �������� ���������� �������
void readTestFile(string testName, myType** &matrix, int &n);
