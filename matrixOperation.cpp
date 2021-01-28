#include <iostream>
#include <fstream>
#include <cmath>

#define myType float

using namespace std;

//������� ������//

//������ ������� c ������ ������
void printMatrixAB(myType** matrix, int n)
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n+1; j++)
        cout << matrix[i][j] << "       ";
        cout << endl;
    }
    cout << endl;
}

//������ ���������� �������
void printMatrix(myType** matrix, int n)
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        cout << matrix[i][j] << "       ";
        cout << endl;
    }
    cout << endl;
}

//������ �������
void printVector(myType* x, int n)
{
	for (int i = 0; i < n; i++)
		cout << x[i] << "   ";
	cout << endl << endl;
}

//����� ��������� ��� ����������� ������
void degenerateMatrix(myType** matrix, int n)
{
    cout << "����������� ������:" << endl;
    printMatrixAB(matrix, n);
    cout << "cond(A) = inf" << endl;
}



//������� ���������������� � ��������//

//����������� ������� � �������
void testSolution(myType** matrix, myType* x, myType* b1, int n)
{
	for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            b1[i] += x[j] * matrix[i][j];
}

//��������� ������ (�������� ���������, ��� ������� � �������)
void MatrMult(myType** A, myType** B, int n)
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
    {
        for(int j = 0; j < n; j++)
            cout << C[i][j] << "    ";
        cout << endl;
    }
    cout << endl;

    for(int i = 0; i < n; i++)
        delete [] C[i];
}

//���������� ������ ����� ��� ���������� ||b||
void perturbationRightPartMax(myType** matrix, int n)
{
    for(int i = 0; i < n; i++)
        if(matrix[i][n] > 0)
            matrix[i][n] += 0.001;
        else
            matrix[i][n] -= 0.001;
}

//���������� ������ ����� ��� ���������� ||b||
void perturbationRightPartMin(myType** matrix, int n)
{
    for(int i = 0; i < n; i++)
        if(matrix[i][n] < 0)
            matrix[i][n] += 0.001;
        else
            matrix[i][n] -= 0.001;
}





//����� �������� � ����������//

//�������������� ����� �������
myType normVector1(myType* vec, int n)
{
    myType sum = 0;

    for(int i = 0; i < n; i++)
        sum += fabs(vec[i]);

    return sum;
}

//�������������� ����� �������
myType normMatrix1(myType** matrix, int n)
{
    myType* vec = new myType [n];
    myType max = 0, norm;

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
            vec[j] = matrix[j][i];

        norm = normVector1(vec, n);

        if(norm > max)
            max = norm;
    }

    delete [] vec;

    return max;
}

//������� �����
myType normVector2(myType* vec, int n)
{
    myType sum = 0;

    for(int i = 0; i < n; i++)
        sum += vec[i]*vec[i];

    return sqrt(sum);
}

//���������� ����� �������
myType normVectorInf(myType* vec, int n)
{
    myType max = 0;

    for(int i = 0; i < n; i++)
        if(fabs(vec[i]) > max)
            max = fabs(vec[i]);

    return max;
}

//���������� ����� �������
myType normMatrixInf(myType** matrix, int n)
{
    myType* vec = new myType [n];
    myType max = 0, norm;

    for(int i = 0; i < n; i++)
    {
        norm = normVector1(matrix[i], n);

        if(norm > max)
            max = norm;
    }

    delete [] vec;

    return max;
}




//������� ���������� � �������//

//����� ����� �� ������
void choiceFile(string &testName)
{
	int num;

	cout << "�������� ����:" << endl
		<< "1. D1.txt" << endl
		<< "2. D2.txt" << endl
		<< "3. D3.txt" << endl
		<< "4. D4.txt" << endl
		<< "5. D5.txt" << endl
		<< "6. DATA11.txt" << endl
		<< "7. P_DAT11.txt" << endl
		<< "8. test.txt" << endl
		<< "9. test1.txt" << endl;

	cin >> num;
	if (num > 9 || num < 1)
		num = 1;

	switch (num)
	{
        case 1: testName = "tests\\D1.txt";
            break;
        case 2: testName = "tests\\D2.txt";
            break;
        case 3: testName = "tests\\D3.txt";
            break;
        case 4: testName = "tests\\D4.txt";
            break;
        case 5: testName = "tests\\D5.txt";
            break;
        case 6: testName = "tests\\DATA11.txt";
            break;
        case 7: testName = "tests\\P_DAT11.txt";
            break;
        case 8: testName = "tests\\test.txt";
            break;
        case 9: testName = "tests\\test1.txt";
            break;
	}
}

//������ �� ����� � �������� ���������� �������
void readTestFile(string testName, myType** &matrix, int &n)
{
	ifstream test(testName.c_str());
	int m;
	test >> m;

	if (n != m)
	for (int i = 0; i < n; i++)
		delete[] matrix[i];

    n = m;

	matrix = new myType*[n];
        for (int i = 0; i < n; i++)
            matrix[i] = new myType[n + 1];

	for (int i = 0; i < n; i++)
	for (int j = 0; j < n + 1; j++)
		test >> matrix[i][j];
}
