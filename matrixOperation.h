#pragma once

#define myType float

using namespace std;

//ФУНКЦИИ ПЕЧАТИ//
//печать матрицы c правой частью
void printMatrixAB(myType** matrix, int n);
//печать квадратной матрицы
void printMatrix(myType** matrix, int n);
//печать решения
void printVector(myType* x, int n);
//Вывод сообщения при вырожденном случае
void degenerateMatrix(myType** matrix, int n);


//подстановка решения в систему
void testSolution(myType** matrix, myType* x, myType* b1, int n);
//Умножение матриц (печатает результат, без запииси в матрицу)
void MatrMult(myType** A, myType** B, int n);
//Возмущение правой части
void perturbationRightPartMax(myType** matrix, int n);
void perturbationRightPartMin(myType** matrix, int n);


//НОРМЫ ВЕКТОРОВ И ОПЕРАТОРОВ//
//Октаэдрическая норма
myType normVector1(myType* vec, int n);
myType normMatrix1(myType** matrix, int n);
//Шаровая норма
myType normVector2(myType* vec, int n);
//Кубическая норма
myType normVectorInf(myType* vec, int n);
myType normMatrixInf(myType** matrix, int n);


//РАБОТА С ФАЙЛАМИ//
//выбор файла из списка
void choiceFile(string &testName);
//Чтение из файла и создание двумерного массива
void readTestFile(string testName, myType** &matrix, int &n);
