#include <iostream>
#include <chrono>
#include <omp.h>

using namespace std;



//! Describes side's size of used matrix.
const int CONFIG_MATRIX_SIZE = 240; // 1440; // 2880

const double CONFIG_MAX_M_VAL =  1000000;    // range of matrix values
const double CONFIG_MIN_M_VAL = -1000000;

const int CONFIG_INITIAL_BLOCK_SIZE = 2;     // inital block size and its change
const int CONFIG_STEP_TO_CHANGE_BLOCK = 1;

const int CONFIG_NUMBER_OF_EXPERIMENTS = 16; // number of evals to find mean



#define EXP_START totalTime = 0;                                         \
                  for (int i = 0; i < CONFIG_NUMBER_OF_EXPERIMENTS; i++) \
                  {                                                      \
                    fill_n(C, squareItemsCount, 0);                      \
                    setTimerStart();

#define EXP_END totalTime += getElapsedTimerTime(); \
                }

#define BLOCK_IT_BEGIN blocksSize = CONFIG_INITIAL_BLOCK_SIZE;     \
                       while (blocksSize < CONFIG_MATRIX_SIZE / 4) \
                       {                                           \
                         if (CONFIG_MATRIX_SIZE % blocksSize == 0) \
                         {

#define BLOCK_IT_END } blocksSize += CONFIG_STEP_TO_CHANGE_BLOCK; }



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



//! Recursive func to convert square index to the vector's type.
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

//! Reads a value from LTM vector with full index as if it was a symmetric matrix.
double getValueOfSymmetricLTM(double * matrix, int i, int j, int sizeArrFormat)
{
    if (j > i)
        return matrix[getFastVectorIndexOfLTM(j, i, sizeArrFormat)];
    else
        return matrix[getFastVectorIndexOfLTM(i, j, sizeArrFormat)];
}



//! Recursive func to get vector index from matrix-type.
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

//! Reads a value from UTM vector.
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
long getElapsedTimerTime()
{
    finish = chrono::high_resolution_clock::now();
    return chrono::duration_cast<chrono::nanoseconds>(finish - start).count();
}



//! Implements multiplication of two matrix. No any optimizations.
void simpleMultM(double * A, double * B, double * C, int n)
{
    int sizeArrFormat = n - 1;
    
    for (int i = 0; i < n; i++)
    {
        int lineIndex = i * n;
        
        for (int j = 0; j < n; j++)
        {
            int totalLineIndex = lineIndex + j;
            
            for (int k = 0; k < n; k++)
                C[totalLineIndex] += (getValueOfSymmetricLTM(A, i, k, sizeArrFormat) * getValueOfUTM(B, k, j, sizeArrFormat));
        }
    }
}

//! Implements multiplication of two matrix using nested blocks.
void blocksMultM(double * A, double * B, double * C, int n, int m)
{
    int sizeArrFormat = n - 1;
    int blocksCount = n / m;
    
    // Iterating throw blocks.
    for (int i = 0; i < blocksCount; i++)
    {
        int im = i * m;
        
        for (int j = 0; j < blocksCount; j++)
        {
            int jm = j * m;
            
            for (int k = 0; k < blocksCount; k++)
            {
                int km = k * m;
            
                // Iterating throw blocks' elements.
                for (int ii = 0; ii < m; ii++)
                {
                    int imii = im + ii;
                    int lineIndex = (imii) * n;
                    
                    for (int jj = 0; jj < m; jj++)
                    {
                        int jmjj = jm + jj;
                        int totalLineIndex = lineIndex + jmjj;
                        
                        for (int kk = 0; kk < m; kk++)
                        {
                            int kmkk = km + kk;
                            
                            C[totalLineIndex] +=
                            getValueOfSymmetricLTM(A, imii, kmkk, sizeArrFormat) *
                            getValueOfUTM(B, kmkk, jmjj, sizeArrFormat);
                        }
                    }
                }
            }
        }
    }
}

//! Implements multiplication of two matrix using nested blocks and also tries to make computations in parallel mode.
void parallelizedBlocksMultM(double * A, double * B, double * C, int n, int m)
{
    int sizeArrFormat = n - 1;
    int blocksCount = n / m;
    
    // Iterating throw blocks.
    for (int i = 0; i < blocksCount; i++)
    {
        int im = i * m;
        
        for (int j = 0; j < blocksCount; j++)
        {
            int jm = j * m;
            
            #pragma omp parallel for
            for (int k = 0; k < blocksCount; k++)
            {
                int km = k * m;
                
                // Iterating throw blocks' elements.
                for (int ii = 0; ii < m; ii++)
                {
                    int imii = im + ii;
                    int lineIndex = (imii) * n;
                    
                    for (int jj = 0; jj < m; jj++)
                    {
                        int jmjj = jm + jj;
                        int totalLineIndex = lineIndex + jmjj;
                        
                        for (int kk = 0; kk < m; kk++)
                        {
                            int kmkk = km + kk;
                            
                            C[totalLineIndex] +=
                            getValueOfSymmetricLTM(A, imii, kmkk, sizeArrFormat) *
                            getValueOfUTM(B, kmkk, jmjj, sizeArrFormat);
                        }
                    }
                }
            }
        }
    }
}



//! Adds some spaces before number, if needed.
char * outputFormatted(int number)
{
    if (number < 10)
        cout << "  " << number;
    else if (number < 100)
        cout << " " << number;
    else cout << number;
    
    return "";
}

//! Outputs provided usual matrix into the console.
void printM(double * matrix, int size)
{
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            cout << matrix[i * size + j] << " ";
        }
        cout << endl;
    }
}

//! Outputs provided symmetric LT matrix into the console.
void printSymmetricLTM(double * matrix, int size)
{
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            cout << getValueOfSymmetricLTM(matrix, i, j, size - 1) << " ";
        }
        cout << endl;
    }
}

//! Outputs provided UT matrix into the console.
void printUTM(double * matrix, int size)
{
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            cout << getValueOfUTM(matrix, i, j, size - 1) << " ";
        }
        cout << endl;
    }
}

//! Outputs all working matrix.
void printSourceAndResult(double * A, double * B, double * C, int size)
{
    cout << endl;
    printSymmetricLTM(A, size);
    cout << endl;
    printUTM(B, size);
    cout << endl;
    printM(C, size);
    cout << endl;
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
    
    
    
    // Disclaimer.
    cout << "Number of experiments for each method - " << CONFIG_NUMBER_OF_EXPERIMENTS << ", matrix size - " << CONFIG_MATRIX_SIZE << endl;
    
    // For further time estimations.
    long totalTime;
    
    // For further tiling.
    int blocksSize;
    
    
    
    // Simple multiplication.
    
    EXP_START
    
        // Making simple multiplication of source matrix without any optimizations.
        simpleMultM(A, B, C, CONFIG_MATRIX_SIZE);

    EXP_END
    
    // Printing average eval time into the console.
    cout << "Average time of the simple multiplication: " << (totalTime / CONFIG_NUMBER_OF_EXPERIMENTS) << endl;
    
    
    
    // Blocks multiplication experiments.
    
    cout << "No parallelization multiplications" << endl;
    
    BLOCK_IT_BEGIN
    
        EXP_START
    
            // Making nested blocks multiplication of source matrix without any optimizations.
            blocksMultM(A, B, C, CONFIG_MATRIX_SIZE, blocksSize);
    
        EXP_END
            
        // Printing average eval time.
        cout << outputFormatted(blocksSize) << " blocks, time: " << (totalTime / CONFIG_NUMBER_OF_EXPERIMENTS) << endl;
    
    BLOCK_IT_END
    
    
    
    // Simple parallelization experiments.
    
    cout << "Using OpenMP pragma for" << endl;
    
    BLOCK_IT_BEGIN
    
        EXP_START
    
            // Making nested blocks multiplication of source matrix with parallel computation optimization.
            parallelizedBlocksMultM(A, B, C, CONFIG_MATRIX_SIZE, blocksSize);
    
        EXP_END
            
        // Printing average eval time.
        cout << outputFormatted(blocksSize) << " blocks, time: " << (totalTime / CONFIG_NUMBER_OF_EXPERIMENTS) << endl;
    
    BLOCK_IT_END
    
    
    
    // Releasing used memory.
    delete[] A;
    delete[] B;
    delete[] C;
    
    
    
    // 2. Сделать нормальный вывод расчетов в файл и подбор размеров блоков. Поставить нужный размер матрицы.
    
    // 3. Распараллелить по OpenMP : перемножение блоков (каждое), параллелить операции над блоками. Многоядерность и векторизация. Также попробовать для этого опции компилятора.
    
    // 4. Определить опт. размер блока для каждого случая
    
    // 5. Брать несколько измерений времени и находить среднее, график зависимости времени от размера блока
    
    // 6. Добавить подробнейшую информацию о спецификации компьютера.
    
    return 0;
}
