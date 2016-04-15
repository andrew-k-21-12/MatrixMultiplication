#include <iostream>
#include <string>
#include <chrono>
#include <omp.h>
#include <fstream>

using namespace std;



//! Describes side's size of used matrix.
const int CONFIG_MATRIX_SIZE = 1440; // 2880

const double CONFIG_MAX_M_VAL =  1000000;    // range of matrix values
const double CONFIG_MIN_M_VAL = -1000000;

const int CONFIG_INITIAL_BLOCK_SIZE = 2;     // inital block size and its change
const int CONFIG_STEP_TO_CHANGE_BLOCK = 1;

const int CONFIG_NUMBER_OF_EXPERIMENTS = 8;  // number of evals to find mean

const string CONFIG_LOG_PATH = "logs/";      // where to store evaluations results



#define EXP_START totalTime = 0;                                         \
                  for (int i = 0; i < CONFIG_NUMBER_OF_EXPERIMENTS; i++) \
                  {                                                      \
                    fill_n(C, squareItemsCount, 0);                      \
                    setTimerStart();

#define EXP_END totalTime += getElapsedTimerTime(); \
                }

#define BLOCK_IT_BEGIN blocksSize = CONFIG_INITIAL_BLOCK_SIZE;      \
                       while (blocksSize < CONFIG_MATRIX_SIZE / 10) \
                       {                                            \
                         if (CONFIG_MATRIX_SIZE % blocksSize == 0)  \
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
void parallelizedMainBlocksMultM(double * A, double * B, double * C, int n, int m)
{
    int sizeArrFormat = n - 1;
    int blocksCount = n / m;
    
    // Iterating throw blocks.
    #pragma omp parallel for
    for (int i = 0; i < blocksCount; i++)
    {
        int im = i * m;
        
        #pragma omp parallel for
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
                        
                        // Very bad idea to make parallelizations here.
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
void parallelizedMainLinesBlocksMultM(double * A, double * B, double * C, int n, int m)
{
    int sizeArrFormat = n - 1;
    int blocksCount = n / m;
    
    // Iterating throw blocks.
    #pragma omp parallel for
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
                        
                        // Very bad idea to make parallelizations here.
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
void parallelizedMainColsBlocksMultM(double * A, double * B, double * C, int n, int m)
{
    int sizeArrFormat = n - 1;
    int blocksCount = n / m;
    
    // Iterating throw blocks.
    for (int i = 0; i < blocksCount; i++)
    {
        int im = i * m;
        
        #pragma omp parallel for
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
                        
                        // Very bad idea to make parallelizations here.
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
void parallelizedNestedLinesBlocksMultM(double * A, double * B, double * C, int n, int m)
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
                #pragma omp parallel for
                for (int ii = 0; ii < m; ii++)
                {
                    int imii = im + ii;
                    int lineIndex = (imii) * n;
                    
                    for (int jj = 0; jj < m; jj++)
                    {
                        int jmjj = jm + jj;
                        int totalLineIndex = lineIndex + jmjj;
                        
                        // Very bad idea to make parallelizations here.
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
void parallelizedMainLinesVectorBlocksMultM(double * A, double * B, double * C, int n, int m)
{
    int sizeArrFormat = n - 1;
    int blocksCount = n / m;
    
    // Iterating throw blocks.
    #pragma omp parallel for simd
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
                        
                        // Very bad idea to make parallelizations here.
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
const char * outputFormatted(int number)
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
    // To make different logs creating a dir.
    string logFileNamePrefix = "default";
    if (argc >= 2)
        logFileNamePrefix = argv[1];
    string pathToStoreLog = CONFIG_LOG_PATH + logFileNamePrefix;
    system(("mkdir -p " + pathToStoreLog).c_str());
    pathToStoreLog += "/";
    
    
    
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
    
    // Saving result to the log file.
    ofstream sf;
    sf.open(pathToStoreLog + "simple.csv");
    sf << "Experiments,Size(n),Time(ns)\n";
    sf << CONFIG_NUMBER_OF_EXPERIMENTS << "," << CONFIG_MATRIX_SIZE << "," << (totalTime / CONFIG_NUMBER_OF_EXPERIMENTS) << "\n";
    sf.close();
    
    // To check parallel computation results.
    double * checkC = new double[squareItemsCount];
    copy(C, C + squareItemsCount, checkC);
    
    
    
    // Blocks multiplication experiments.
    
    cout << "No parallelization multiplications" << endl;
    
    // Saving result to the log file.
    ofstream bf;
    bf.open(pathToStoreLog + "blocks.csv");
    bf << "Experiments,Size(n),Blocksize(m),Time(ns)\n";
    
    BLOCK_IT_BEGIN
    
        EXP_START
    
            // Making nested blocks multiplication of source matrix without any optimizations.
            blocksMultM(A, B, C, CONFIG_MATRIX_SIZE, blocksSize);
    
        EXP_END
            
        // Printing average eval time.
        cout << outputFormatted(blocksSize) << " blocks, time: " << (totalTime / CONFIG_NUMBER_OF_EXPERIMENTS) << endl;
    
        // Writing result to file.
        bf << CONFIG_NUMBER_OF_EXPERIMENTS << "," << CONFIG_MATRIX_SIZE << "," << blocksSize << "," << (totalTime / CONFIG_NUMBER_OF_EXPERIMENTS) << "\n";
    
    BLOCK_IT_END
    
    bf.close();
    
    
    
    // Main lines and columns parallelization experiments.
    
    cout << "Using OpenMP pragma for main lines and columns" << endl;
    
    // Saving result to the log file.
    ofstream bmlcf;
    bmlcf.open(pathToStoreLog + "blocks_p_main_lines_cols.csv");
    bmlcf << "Experiments,Size(n),Blocksize(m),Time(ns)\n";
    
    BLOCK_IT_BEGIN
    
        EXP_START
    
            // Making nested blocks multiplication of source matrix with parallel computation optimization.
            parallelizedMainBlocksMultM(A, B, C, CONFIG_MATRIX_SIZE, blocksSize);
    
        EXP_END
            
        // Printing average eval time.
        cout << outputFormatted(blocksSize) << " blocks, time: " << (totalTime / CONFIG_NUMBER_OF_EXPERIMENTS) << endl;
    
        // Writing result to file.
        bmlcf << CONFIG_NUMBER_OF_EXPERIMENTS << "," << CONFIG_MATRIX_SIZE << "," << blocksSize << "," << (totalTime / CONFIG_NUMBER_OF_EXPERIMENTS) << "\n";
    
    BLOCK_IT_END
    
    bmlcf.close();
    
    
    
    // Main lines parallelization experiments.
    
    cout << "Using OpenMP pragma for main lines" << endl;
    
    // Saving result to the log file.
    ofstream bmlf;
    bmlf.open(pathToStoreLog + "blocks_p_main_lines.csv");
    bmlf << "Experiments,Size(n),Blocksize(m),Time(ns)\n";
    
    BLOCK_IT_BEGIN
    
        EXP_START
    
            // Making nested blocks multiplication of source matrix with parallel computation optimization.
            parallelizedMainLinesBlocksMultM(A, B, C, CONFIG_MATRIX_SIZE, blocksSize);
    
        EXP_END
    
        // Printing average eval time.
        cout << outputFormatted(blocksSize) << " blocks, time: " << (totalTime / CONFIG_NUMBER_OF_EXPERIMENTS) << endl;
    
        // Writing result to file.
        bmlf << CONFIG_NUMBER_OF_EXPERIMENTS << "," << CONFIG_MATRIX_SIZE << "," << blocksSize << "," << (totalTime / CONFIG_NUMBER_OF_EXPERIMENTS) << "\n";
    
    BLOCK_IT_END
    
    bmlf.close();
    
    
    
    // Main columns parallelization experiments.
    
    cout << "Using OpenMP pragma for main columns" << endl;
    
    // Saving result to the log file.
    ofstream bmcf;
    bmcf.open(pathToStoreLog + "blocks_p_main_cols.csv");
    bmcf << "Experiments,Size(n),Blocksize(m),Time(ns)\n";
    
    BLOCK_IT_BEGIN
    
        EXP_START
    
            // Making nested blocks multiplication of source matrix with parallel computation optimization.
            parallelizedMainColsBlocksMultM(A, B, C, CONFIG_MATRIX_SIZE, blocksSize);
    
        EXP_END
    
        // Printing average eval time.
        cout << outputFormatted(blocksSize) << " blocks, time: " << (totalTime / CONFIG_NUMBER_OF_EXPERIMENTS) << endl;
    
        // Writing result to file.
        bmcf << CONFIG_NUMBER_OF_EXPERIMENTS << "," << CONFIG_MATRIX_SIZE << "," << blocksSize << "," << (totalTime / CONFIG_NUMBER_OF_EXPERIMENTS) << "\n";
    
    BLOCK_IT_END
    
    bmcf.close();
    
    
    
    // Nested lines parallelization experiments.
    
    cout << "Using OpenMP pragma for nested lines" << endl;
    
    // Saving result to the log file.
    ofstream bnlf;
    bnlf.open(pathToStoreLog + "blocks_p_nested_lines.csv");
    bnlf << "Experiments,Size(n),Blocksize(m),Time(ns)\n";
    
    BLOCK_IT_BEGIN
    
        EXP_START
    
            // Making nested blocks multiplication of source matrix with parallel computation optimization.
            parallelizedNestedLinesBlocksMultM(A, B, C, CONFIG_MATRIX_SIZE, blocksSize);
    
        EXP_END
    
        // Printing average eval time.
        cout << outputFormatted(blocksSize) << " blocks, time: " << (totalTime / CONFIG_NUMBER_OF_EXPERIMENTS) << endl;
    
        // Writing result to file.
        bnlf << CONFIG_NUMBER_OF_EXPERIMENTS << "," << CONFIG_MATRIX_SIZE << "," << blocksSize << "," << (totalTime / CONFIG_NUMBER_OF_EXPERIMENTS) << "\n";
    
    BLOCK_IT_END
    
    bnlf.close();
    
    
    
    // Main lines parallelization experiments with vectorization.
    
    cout << "Using OpenMP pragma for main lines with vectorization" << endl;
    
    // Saving result to the log file.
    ofstream bmlvf;
    bmlvf.open(pathToStoreLog + "blocks_p_main_lines_vector.csv");
    bmlvf << "Experiments,Size(n),Blocksize(m),Time(ns)\n";
    
    BLOCK_IT_BEGIN
    
        EXP_START
    
            // Making nested blocks multiplication of source matrix with parallel computation optimization.
            parallelizedMainLinesVectorBlocksMultM(A, B, C, CONFIG_MATRIX_SIZE, blocksSize);
    
        EXP_END
    
        // Printing average eval time.
        cout << outputFormatted(blocksSize) << " blocks, time: " << (totalTime / CONFIG_NUMBER_OF_EXPERIMENTS) << endl;
    
        // Writing result to file.
        bmlvf << CONFIG_NUMBER_OF_EXPERIMENTS << "," << CONFIG_MATRIX_SIZE << "," << blocksSize << "," << (totalTime / CONFIG_NUMBER_OF_EXPERIMENTS) << "\n";
    
    BLOCK_IT_END
    
    bmlvf.close();
    
    
    
    // Checking if parallel computation is correct.
    bool failure = false;
    for (int i = 0; i < squareItemsCount; i++)
        if (C[i] != checkC[i])
        {
            failure = true;
            break;
        }
    if (failure)
        cout << "Parallel computations are invalid!" << endl;
    
    
    
    // Releasing used memory.
    delete[] A;
    delete[] B;
    delete[] C;
    delete[] checkC;
    
    
    
    // 1. Определить опт. размер блока для каждого случая
    
    // 2. График зависимости времени от размера блока
    
    // 3. Добавить подробнейшую информацию о спецификации компьютера
    
    return 0;
}
