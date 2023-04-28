// метод нахождения собственных векторов симметричной матрицы - Метод вращений (Метод Якоби)

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iomanip>

using namespace std;

void init_matrix(double** matrix, int N)
{
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            cin >> matrix[i][j];
        }
    }
}

void is_matrix_symmetrical(int N, double** A)
{
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            if(A[i][j] != A[j][i])
            {
                cout << "Matrix is not symmetrical!!!" << endl;
                exit(0);
            }
        }
    }
}

// for symmetrical matrices
int rotation_method( double **A, int N, double **solution, double precision )
// A - ИЗНАЧАЛЬНАЯ МАТРИЦА, размерность матрицы - N, МАТРИЦА ИЗ СОБСТВЕННЫЗ ВЕКТОРОВ - solution, точность - precision
{
    int result = 1; // количество итерация
    int i, j, k;
    int maxI, maxJ;
    double max, fi;

    double** rotation_matrix;
    rotation_matrix = new double*[N];
    for ( i = 0; i < N; i ++ ) {
        rotation_matrix[i] = new double[N];}

    double** temp;
    temp = new double*[N];
    for ( i = 0; i < N; i ++ ) {
        temp[i] = new double[N];}

    // текущая точность - матричная мера
    double fault = 0.0;
    for ( i = 0; i < N; i ++ ) {
        for ( j = i + 1; j < N; j ++ ) {
            fault = fault + A[i][j]*A[i][j];}}
    fault = sqrt( 2*fault );

    // основной итерационный процесс-цикл
    while ( fault > precision ) {

        // находим максимальный наддиагональный элемент
        max = 0.0;
        for ( i = 0; i < N; i ++ )
        {
            for ( j = i + 1; j < N; j ++ )
            {
                if ( A[i][j] > 0 && A[i][j] > max )
                {
                    max = A[i][j];
                    maxI = i;
                    maxJ = j;
                }
                else if ( A[i][j] < 0 && - A[i][j] > max )
                {
                    max = - A[i][j];
                    maxI = i;
                    maxJ = j;
                }
            }
        }

        // создаём матрицу вращений

        // делаем матрицу rotation_matrix едничной
        for ( i = 0; i < N; i ++ )
        {
            for ( j = 0; j < N; j ++ )
            {
                rotation_matrix[i][j] = 0;
            }
            rotation_matrix[i][i] = 1;
        }

        // вставляем в соответствующие места под индексами максимального наддиагонального значения
        if ( A[maxI][maxI] == A[maxJ][maxJ] )
        {
            rotation_matrix[maxI][maxI] = rotation_matrix[maxJ][maxJ] =
            rotation_matrix[maxJ][maxI] = sqrt( 2.0 ) / 2.0;
            rotation_matrix[maxI][maxJ] = - sqrt( 2.0 ) / 2.0;
        }
        else
        {
            // http://mathhelpplanet.com/static.php?p=metody-resheniya-zadach-o-sobstvennykh-znacheniyakh-i-vektorakh-matritsy
            fi = 0.5 * atan((2.0 * A[maxI][maxJ]) / (A[maxI][maxI] - A[maxJ][maxJ])); // вычисляем угол fi
            // ставим на соответствующие места косинусы и синусы
            rotation_matrix[maxI][maxI] = rotation_matrix[maxJ][maxJ] = cos( fi );
            rotation_matrix[maxI][maxJ] = -sin( fi );
            rotation_matrix[maxJ][maxI] = sin( fi );
        }

        // заполняем матрицу temp нулями
        for ( i = 0; i < N; i ++ )
        {
            for ( j = 0; j < N; j ++ )
            {
                temp[i][j] = 0.0;
            }
        }

        // temp = rotation_matrix * A - умножение матрицы поворота слева
        for ( i = 0; i < N; i ++ )
        {
            for ( j = 0; j < N; j ++ )
            {
                for ( k = 0; k < N; k ++ )
                {
                    temp[i][j] = temp[i][j] + rotation_matrix[k][i] * A[k][j];
                }
            }
        }

        // ЗАНУЛЯЕМ МАТРИЦУ A
        for ( i = 0; i < N; i ++ )
        {
            for ( j = 0; j < N; j ++ )
            {
                A[i][j] = 0.0;
            }
        }

        // A = temp * rotation_matrix - умножение матрицы поворота справа
        for ( i = 0; i < N; i ++ )
        {
            for ( j = 0; j < N; j ++ )
            {
                for ( k = 0; k < N; k ++ )
                {
                    A[i][j] = A[i][j] + temp[i][k] * rotation_matrix[k][j];
                }
            }
        }

        // вычисляем текущую точность

        // корень из суммы квадратов всех элементов матрицы
        fault = 0.0;
        for ( i = 0; i < N; i ++ )
        {
            for ( j = i + 1; j < N; j ++ )
            {
                fault = fault + A[i][j]*A[i][j];
            }
        }
        fault = sqrt( 2*fault );

        // обнуляем элепенты матрицы temp
        for ( i = 0; i < N; i ++ )
        {
            for ( j = 0; j < N; j ++ )
            {
                temp[i][j] = 0.0;
            }
        }

        // temp = solution * rotation_matrix - вычисляем матрицу, в которой хранятся собственные вектор
        // Solution = H1 * H2 * ...
        for ( i = 0; i < N; i ++ )
        {
            for ( j = 0; j < N; j ++ )
            {
                for ( k = 0; k < N; k ++ )
                {
                    temp[i][j] = temp[i][j] + solution[i][k] * rotation_matrix[k][j];
                }
            }
        }

        // результат тройного цикла выше записываем в матрицу собственных векторов
        for ( i = 0; i < N; i ++ )
        {
            for ( j = 0; j < N; j ++ )
            {
                solution[i][j] = temp[i][j];
            }
        }

        result++; // количество итераций
    }

    return result;
}

int main()
{
    int i, j;
    int size;
    double **A, **solution, precision;
    cin >> size;
    A = new double*[size];
    solution = new double*[size]; // Матрица из собственных векторов
    for ( i = 0; i < size; i++ )
    {
        A[i] = new double[size];
        solution[i] = new double[size];
    }

    for ( i = 0; i < size; i ++ )
    {
        for ( j = 0; j < size; j ++ )
        {
            solution[i][j] = 0;
        }
        solution[i][i] = 1;
    }

    init_matrix(A, size);

    precision = 0.00000001;
    is_matrix_symmetrical(size, A);

        int steps = rotation_method( A, size, solution, precision );
        for ( i = 0; i < size; i++ )
        {
            cout << "Eigenvector v_" << i + 1 << ":\n";
            for ( j = 0; j < size; j ++ )
            {
                cout << fixed << setprecision(8) << solution[j][i] << "\n";
            }
        }
        cout << "Eigenvalues: \n";
        for ( i = 0; i < size; i++ )
        {
            cout << fixed << setprecision(8) << A[i][i] << "\n";
        }
        cout << "Iterations: " << steps;

        return 0;
}