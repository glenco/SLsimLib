#ifndef _UTILITIES_H
#define _UTILITIES_H

#include "standard.h"

//typedef double PosType;

template <class T>
void Matrix(T **matrix, long rows, long cols){
  matrix = new T*[rows];
  matrix[0] = new T[rows*cols];
  for (long i = 1; i < rows; ++i)
    matrix[i] = matrix[0] + i * cols;
}

template <class T>
void free_Matrix(T **matrix, long rows, long cols){
  if (rows) delete[] matrix[0];
  delete []matrix;
}

PosType **PosTypeMatrix(long rows, long cols);

void free_PosTypeMatrix(PosType **matrix, long rows, long cols);

#endif
