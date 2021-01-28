#pragma once

#define myType float

bool qrFactorization(myType** matrix, myType** T, myType* x, int n);

//транспонирование
void transpose_QR(myType** matrix, int n);
//Умножение правой части на матрицу поворота, правая часть перезаписывается
void matrixMultB_QR(myType** T, myType** matrix, int n);

//напюрмши унд//
void backward_QR(myType** matrix, myType* x, int n);
