#include <iostream>
#include "matrixOperation.h"
#include "gaussmethod.h"
#include "qrmethod.h"

#define myType float

using namespace std;

void inverseMatrix(myType** matrix, myType** Q, myType** invMatrix, int n);
void printResult(string &testName, myType** matrix, myType* x, myType* b1, myType* discr, int n);

int main()
{
	setlocale(LC_ALL, "RUSSIAN");

	//имя файла с системой
	string testName;
        choiceFile(testName);
	//размерность системы, выбор метода
	int n = 0;
	//создание двумерного массива под систему (вместе с правой частью)
	myType **matrix, *x, *b1, *discr;
        readTestFile(testName, matrix, n);

    discr = new myType [n];
    x  = new myType [n];
    b1 = new myType [n];
    for(int i = 0; i < n; i++)
        b1[i] = 0;


    cout << "МЕТОД ГАУССА" << endl;
    cout << "Система в изначальном виде:" << endl;
    printMatrixAB(matrix, n);
    if(gauss(matrix, x, n))
        printResult(testName, matrix, x, b1, discr, n);
    else
        degenerateMatrix(matrix, n);


    cout << endl << endl << endl << endl << endl;


    readTestFile(testName, matrix, n);
    for(int i = 0; i < n; i++)
        b1[i] = 0;

    //Q - ортогональная, R - верхнетреугольная
    myType **Q = new myType* [n];
        for(int i = 0; i < n; i++)
            Q[i] = new myType [n];

    cout << "МЕТОД QR-РАЗЛОЖЕНИЯ" << endl;
    cout << "Система в изначальном виде:" << endl;
    printMatrixAB(matrix, n);
    if(qrFactorization(matrix, Q, x, n))
    {
        cout << "Ортогональная матрица Q:" << endl;
        printMatrix(Q, n);
        printResult(testName, matrix, x, b1, discr, n);
    }
    else
        degenerateMatrix(matrix, n);


    cout << endl << endl << endl << endl << endl;


    readTestFile(testName, matrix, n);

    cout << "ОЦЕНКА ЧИСЛА ОБУСЛОВЛЕННОСТИ" << endl;
    cout << "Матрица A:" << endl;
    printMatrix(matrix, n);

    myType** invMatrix = new myType* [n];
        for(int i = 0; i < n; i++)
            invMatrix[i] = new myType[n];


    readTestFile(testName, matrix, n);
    qrFactorization(matrix, Q, x, n);
    cout << "Обратная матрица" << endl;
    inverseMatrix(matrix, Q, invMatrix, n);
    printMatrix(invMatrix, n);

    readTestFile(testName, matrix, n);
    cout << "Произведение матриц A^-1*A" << endl;
    MatrMult(invMatrix, matrix, n);

    cout << "Cond_1   = " << normMatrix1(matrix, n) * normMatrix1(invMatrix, n) << endl;
    cout << "Cond_inf = " << normMatrixInf(matrix, n) * normMatrixInf(invMatrix, n) << endl;

    //ОЦЕНКА ЧЕРЕЗ РЕШЕНИЕ СИСТЕМЫ С ШЕВЕЛЕНИЕМ ПРАВОЙ ЧАСТИ//
    myType maxCond1 = 0, maxCond2 = 0, maxCondInf = 0, subMaxCond1, subMaxCond2, subMaxCondInf;
    myType *x1 = new myType [n];
    myType *subB = new myType [n];
        for(int i = 0; i < n; i++)\
            subB[i] = matrix[i][n];

    qrFactorization(matrix, Q, x, n);
    transpose_QR(Q, n);

    for(int i = 0; i < n; i++)
    {
        if(i != 0)
            subB[i-1] -= 0.001;
        subB[i] += 0.001;

        for(int j = 0; j < n; j++)
            matrix[j][n] = subB[j];

        matrixMultB_QR(Q, matrix, n);
        backward_QR(matrix, x1, n);

        for(int j = 0; j < n; j++)
            discr[j] = x[j] - x1[j];
        subMaxCond1   = 1000*normVector1(discr, n)   * normVector1(b1, n)   / normVector1(x, n);
        subMaxCondInf = 1000*normVectorInf(discr, n) * normVectorInf(b1, n) / normVectorInf(x, n);
		subMaxCond2 = 1000 * normVector2(discr, n)   * normVector2(b1, n) / normVector2(x, n);

        if(subMaxCond1 > maxCond1)
            maxCond1 = subMaxCond1;
        if(subMaxCondInf > maxCondInf)
            maxCondInf = subMaxCondInf;
		if (subMaxCond2 > maxCond2)
			maxCond2 = subMaxCond2;
    }
    subB[n-1] -= 0.001;

    for(int i = 0; i < n; i++)
    {
        if(i != 0)
            subB[i-1] += 0.001;
        subB[i] -= 0.001;

        for(int j = 0; j < n; j++)
            matrix[j][n] = subB[j];

        matrixMultB_QR(Q, matrix, n);
        backward_QR(matrix, x1, n);

        for(int j = 0; j < n; j++)
            discr[j] = x[j] - x1[j];
        subMaxCond1   = 1000*normVector1(discr, n)   * normVector1(b1, n)   / normVector1(x, n);
        subMaxCondInf = 1000*normVectorInf(discr, n) * normVectorInf(b1, n) / normVectorInf(x, n);
		subMaxCond2 = 1000 * normVector2(discr, n)   * normVector2(b1, n) / normVector2(x, n);

        if(subMaxCond1 > maxCond1)
            maxCond1 = subMaxCond1;
        if(subMaxCondInf > maxCondInf)
            maxCondInf = subMaxCondInf;
		if (subMaxCond2 > maxCond2)
			maxCond2 = subMaxCond2;
    }
    subB[n-1] += 0.001;

    cout << "dx/db_1 = " << maxCond1 << endl;
    cout << "dx/db_inf = " << maxCondInf << endl;
	cout << "dx/db_2 = " << maxCond2 << endl;

    delete [] subB;


/*
    myType maxCond1 = 0, maxCondInf = 0, subMaxCond1, subMaxCondInf;
    myType *x1 = new myType [n];

    for(int i = 0; i < n; i++)
    {
        readTestFile(testName, matrix, n);
        for(int j = 0; j < n; j++)
            b1[j] = matrix[j][n];
        matrix[i][n] += 0.001;
        gauss(matrix, x1, n);

        for(int j = 0; j < n; j++)
            discr[j] = x[j] - x1[j];
        subMaxCond1   = 1000*normVector1(discr, n)   * normVector1(b1, n)   / normVector1(x, n);
        subMaxCondInf = 1000*normVectorInf(discr, n) * normVectorInf(b1, n) / normVectorInf(x, n);

        if(subMaxCond1 > maxCond1)
            maxCond1 = subMaxCond1;
        if(subMaxCondInf > maxCondInf)
            maxCondInf = subMaxCondInf;
    }

    for(int i = 0; i < n; i++)
    {
        readTestFile(testName, matrix, n);
        for(int j = 0; j < n; j++)
            b1[j] = matrix[j][n];
        matrix[i][n] -= 0.001;
        gauss(matrix, x1, n);

        for(int j = 0; j < n; j++)
            discr[j] = x[j] - x1[j];
        subMaxCond1   = 1000*normVector1(discr, n)   * normVector1(b1, n)   / normVector1(x, n);
        subMaxCondInf = 1000*normVectorInf(discr, n) * normVectorInf(b1, n) / normVectorInf(x, n);

        if(subMaxCond1 > maxCond1)
            maxCond1 = subMaxCond1;
        if(subMaxCondInf > maxCondInf)
            maxCondInf = subMaxCondInf;
    }

    cout << "cond1 = " << maxCond1 << endl;
    cout << "condInf = " << maxCondInf << endl;

*/

	//очистка памяти
	for(int i = 0; i < n; i++)
    {
		delete[] matrix[i];
		delete[] invMatrix[i];
		delete[] Q[i];
    }
	delete[] x;
	delete[] x1;
	delete[] b1;
	delete[] discr;

	cin.get();
	cin.get();

	return 0;
}

void inverseMatrix(myType** matrix, myType** Q, myType** invMatrix, int n)
{
    myType* x = new myType [n];

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
            matrix[j][n] = Q[i][j];

        backward_QR(matrix, x, n);

        for(int j = 0; j < n; j++)
            invMatrix[j][i] = x[j];
    }

    delete[] x;
}

void printResult(string &testName, myType** matrix, myType* x, myType* b1, myType* discr, int n)
{
    cout << "Система приведённая к верхнетреугольному виду:" << endl;
    printMatrixAB(matrix, n);

    cout << "x = ";
    printVector(x, n);

    readTestFile(testName, matrix, n);

    cout << "b1 = ";
    testSolution(matrix, x, b1, n);
    printVector(b1, n);

    cout << "b - b1 = ";
    for(int i = 0; i < n; i++)
        discr[i] = matrix[i][n] - b1[i];
    printVector(discr, n);

    cout << "||b - b1||_1   = " << normVector1(discr, n)   << endl;
    cout << "||b - b1||_inf = " << normVectorInf(discr, n) << endl;
	cout << "||b - b1||_2   = " << normVector2(discr, n)   << endl;
}
