#include <iostream>

#include "libdynamix_outputs.h"

#define DIM 3
#define N 2
#define M 3

int main () {
 // file name variable for all outputs
 char * fileName;

 double scalar = 1.0;
 double vector [DIM] = {1.0, 2.0, 3.0};
 double times [DIM]= {0.0, 1.0, 2.0};
 double matrix [N*M]= {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};

 std::cout << "\nStarting tests...\n";

 fileName = "scalar.out";
 std::cout << "\nWriting a scalar to " << fileName << "...\n";
 outputDScalar(fileName, scalar);

 fileName = "vector.out";
 std::cout << "\nWriting a vector to " << fileName << "...\n";
 outputDVector(fileName, vector, DIM);

 fileName = "vectorT.out";
 std::cout << "\nWriting a vector transpose to " << fileName << "...\n";
 outputDVectorT(fileName, vector, DIM);

 fileName = "timeVectorT.out";
 std::cout << "\nWriting a time-dependent vector to " << fileName << "...\n";
 outputTimeDVector(fileName, times, vector, DIM);

 fileName = "matrix.out";
 std::cout << "\nWriting a matrix to " << fileName << "...\n";
 outputDMatrix(fileName, matrix, N, M);

 fileName = "matrixT.out";
 std::cout << "\nWriting a matrix transpose to " << fileName << "...\n";
 outputDMatrixT(fileName, matrix, N, M);

 std::cout << "\nTests finished!\n\n";

 return 0;
}
