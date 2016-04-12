#include <iostream>
#include <chrono>

using namespace std;



//! Describes side's size of used matrix.
const int CONFIG_MATRIX_SIZE = 1440; // 2880
const double CONFIG_MAX_M_VAL =  10;
const double CONFIG_MIN_M_VAL = -10;



//! Generates some random value in range.
double genRandom(double from, double to)
{
    double f = (double) rand() / RAND_MAX;
    return (from + f * (to - from));
}

//! Evaluates the number of elements in the triangular matrix with provided size.
int evalItemsCount(int size)
{
    int k = 0;
    
    for (int i = 0; i < size; i++)
        for (int j = 0; j <= i; j++)
            k++;
    
    return k;
}

//! Fills provided matrix as a vector.
void fillM(int elements, double * matrix)
{
    for (int i = 0; i < elements; i++)
        matrix[i] = genRandom(CONFIG_MIN_M_VAL, CONFIG_MAX_M_VAL);
}



//! Recursive func to convert square index to the vector's type. 100% works.
int getVectorIndexOfLTM(int i, int j)
{
    if (i == 0)
        return 0;
    
    return i + j + getVectorIndexOfLTM(i - 1, 0);
}

//! Converts matrix's indexes to vector's. Non-recursive.
int getFastVectorIndexOfLTM(int i, int j, int sizeArrFormat)
{
    if (i == 0)
        return 0;
    
    return i * (sizeArrFormat + 1)            // square of the previous rect
    // Skipped empties (rect + triangle).
    - ((sizeArrFormat - i) * i                // rect
       + (int) (0.5 * (i + 0.5) * (i + 0.5))) // triangle
    + j;                                      // current index shift
}

//! Reads a value from LTM vector with full index as if it was a symmetric matrix. 100% works.
double getValueOfSymmetricLTM(double * matrix, int i, int j, int sizeArrFormat)
{
    if (j > i)
        return matrix[getFastVectorIndexOfLTM(j, i, sizeArrFormat)];
    else
        return matrix[getFastVectorIndexOfLTM(i, j, sizeArrFormat)];
}



//! Recursive func to get vector index from matrix-type. 100% works.
int getVectorIndexOfUTM(int i, int j, int sizeArrFormat)
{
    if (i == 0)
        return j;
    
    return getVectorIndexOfUTM(i - 1, sizeArrFormat, sizeArrFormat) + j - (i - 1);
}

//! Fast implementation of matrix coords to vector convertion.
int getFastVectorIndexOfUTM(int i, int j, int sizeArrFormat)
{
    return (i + 1) * (sizeArrFormat + 1)   // square of the total rect
    - 1                                    // required shift
    - (sizeArrFormat - j)                  // right shift
    - (int) (0.5 * (i + 0.5) * (i + 0.5)); // full excluded points
    
}

//! Reads a value from UTM vector. 100% works.
double getValueOfUTM(double * matrix, int i, int j, int sizeArrFormat)
{
    if (i > j)
        return 0;
    
    return matrix[getFastVectorIndexOfUTM(i, j, sizeArrFormat)];
}



auto  start = chrono::high_resolution_clock::now();
auto finish = chrono::high_resolution_clock::now();
void setTimerStart()
{
    start = chrono::high_resolution_clock::now();
}
void getElapsedTimerTime()
{
    finish = chrono::high_resolution_clock::now();
    cout << chrono::duration_cast<chrono::nanoseconds>(finish - start).count() << " ns" << endl;
}



//! Implements multiplication of two matrix.
void simpleMultM(double * A, double * B, double * C, int mSize)
{
    int sizeArrFormat = mSize - 1;
    
    for (int i = 0; i < mSize; i++)
    {
        int lineIndex = i * mSize;
        
        for (int j = 0; j < mSize; j++)
        {
            cout << "Evaluating now " << i << " " << j << endl;
            
            for (int k = 0; k < mSize; k++)
                C[lineIndex + j] += (getValueOfSymmetricLTM(A, i, k, sizeArrFormat) * getValueOfUTM(B, k, j, sizeArrFormat));
        }
    }
}



int main(int argc, char * argv[])
{
    // Seeding random values.
    srand(time(NULL));
    
    // Saving number of elements in the triangular matrix with current size.
    int itemsCount = evalItemsCount(CONFIG_MATRIX_SIZE);
    // Size of the result matrix.
    int squareItemsCount = CONFIG_MATRIX_SIZE * CONFIG_MATRIX_SIZE;
    
    // Generating source and result template matrix.
    double * A = new double[itemsCount]; fillM(itemsCount, A); // symmetric, stored as lower triangular
    double * B = new double[itemsCount]; fillM(itemsCount, B); // upper triangular
    double * C = new double[squareItemsCount];                 // result matrix, rectangular
    
    // Everything above work 100%.
    
    
    
    // Clearing result matrix before evaluations.
    fill_n(C, squareItemsCount, 0);
    
    // Setting timer to start counting.
    setTimerStart();
    
    // Making simple multiplication of source matrix without any optimizations.
    simpleMultM(A, B, C, CONFIG_MATRIX_SIZE);
    
    // Getting op. time.
    getElapsedTimerTime();
    
    
    
    // Releasing used memory.
    delete[] A;
    delete[] B;
    delete[] C;
    
    
    
    // 2. Перемножить матрицы по блокам (n - размер матрицы, m - размер блока, n mod m = 0, брать с интервалом хотя бы 4)
    
    // 3. Распараллелить по OpenMP : перемножение блоков (каждое), параллелить операции над блоками. Многоядерность и векторизация. Также попробовать для этого опции компилятора.
    
    // 4. Определить опт. размер блока для каждого случая
    
    // 5. Брать несколько измерений времени и находить среднее, график зависимости времени от размера блока
    
    return 0;
}
