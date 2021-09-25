#ifndef SUBROUTINES_H
#define SUBROUTINES_H

#include <vector>
#include <iostream>
#include <math.h>
#include <float.h> 
#include <cmath>
#include <limits>
#include <iomanip>
#include <random>
#include <ctime>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// making life a bit simpler
typedef std::vector<double> Vec;
typedef std::vector<Vec> Mat;

struct generalResult{
    std::vector< std::string > errors;
    double determinant;
};

struct iterativeResult{
    std::vector< std::string > errors;
    int total_iterations = 0;
    std::vector<double> TOL_history;
};

struct basicResult{
    double estimate;
    std::vector< std::string > errors;
};

/* returns the magnetude(Euclidean norm) of a vector of doubles */
double absVec(Vec &v){
    double square_sum = 0;
    for(size_t i=0; i<v.size(); ++i){
        square_sum  += v[i]*v[i];
    }
    return sqrt(square_sum);
}

/* subtracts values from v by the respective u element and stores in a third vector result (result = v - u) */
void subtractVec(Vec &v, Vec &u, Vec &result){
    for(size_t i=0; i<v.size(); ++i){
        result[i] = v[i] - u[i];
    }
}

/* multiplies two matrices and stores in a third matrix */
void productMatrices(Mat &A, Mat &B, Mat &result){
    size_t A_rows = A.size();
    size_t A_columns = A[0].size();
    size_t B_rows = B[0].size();

    for(size_t j=0; j<B_rows; ++j){
        for(size_t i=0; i<A_rows; ++i){

            result[i][j] = 0;
            for(size_t k=0; k<A_columns; ++k){
                result[i][j] += A[i][k] * B[k][j];
            }

        }
    }

}

/* multiplies a matrix and a vector and stores the result in another vector */
void productMatrixVector(Mat &A, Vec &B, Vec &result){
    size_t A_rows = A.size();
    size_t A_columns = A[0].size();

    for(size_t i=0; i<A_rows; ++i){

        result[i] = 0;
        for(size_t k=0; k<A_columns; ++k){
            result[i] += A[i][k] * B[k];
        }

    }

}

/* inverts lower triangular matrix and stores in another matrix */
void invertLowerMat(Mat &L, Mat &L_inv){

    size_t L_rows = L.size();
    double sum;
    for (size_t i=0; i<L_rows; ++i){

        L_inv[i][i] = 1/L[i][i];
        for (size_t j=0; j<i; ++j){

            sum = 0;
            for (size_t k=j; k<i; ++k){
                sum += L[i][k] * L_inv[k][j];
            }

            L_inv[i][j] = -sum * L_inv[i][i];
        }

    }

}


/* used by determinantOfMatrix */
void subMatrix(Mat &A, Mat &temp, int p, int q, int order) {
   int i = 0, j = 0;
   // filling the sub matrix
   for (int row = 0; row < order; row++) {
      for (int col = 0; col < order; col++) {
         // skipping if the current row or column is not equal to the current
         // element row and column
         if (row != p && col != q) {
            temp[i][j++] = A[row][col];
            if (j == order - 1) {
               j = 0;
               i++;
            }
         }
      }
   }
}

/* calculates determinant of A by subdividing the matrix recursively */
double determinantOfMatrix(Mat &A, int order) {
   double determinant = 0;
   if (order == 1) {
      return A[0][0];
   }
   if (order == 2) {
      return (A[0][0] * A[1][1]) - (A[0][1] * A[1][0]);
   }
   Mat temp(order, Vec(order, 0.0));
   double sign = 1;
   for (int i = 0; i < order; i++) {
      subMatrix(A, temp, 0, i, order);
      determinant += sign * A[0][i] * determinantOfMatrix(temp, order - 1);
      sign = -sign;
   }
   return determinant;
}

/* Check if matrix A is diagonally dominant. */
bool diagonallyDominant(Mat &A, int &order){

    for(size_t i=0; i<order; ++i){

        double row_summation = 0;
        double column_summation = 0;
        for(size_t j=0; j<order; ++j){
            if( j!=i ){
                row_summation += std::abs(A[i][j]);
                column_summation += std::abs(A[j][i]);
            }
        }

        if(A[i][i] < row_summation || A[i][i] < column_summation){
            return false;
        }

    }

    // if the loop completes, A is diagonal dominant
    return true;

}

/* calculates determinant of matrix A(=LU) by multiplying the determinant of L by the determinant of U */
double calculateLUDeterminant(Mat &L, Mat &U, int &order){
    double L_determinant = 1;
    double U_determinant = 1;
    for(size_t i=0; i<order; ++i){
        L_determinant *= L[i][i];
        U_determinant *= U[i][i];
    }

    return L_determinant * U_determinant;
}

/* calculates determinant of A based on its eigenvalues */
double eigenDeterminant(Mat &A, Vec &values, int &order){
    double A_det = 1;
    for(size_t i=0; i<order; ++i){
        A_det *= values[i];
    }
    return A_det;
}

/* check if matrix is symmetric */
bool symmetricMatrix(Mat &A){
    int rows = A.size();
    int columns = A[0].size();

    for(size_t i=0; i<rows; ++i){
        for(size_t j=i+1; j<columns; ++j){
            if(A[i][j] != A[j][i]){
                return false;
            }
        }
    }
    return true;

}

/* returns the phi value based on indices of elements of matrix */
double phiAngle(Mat &A, int &i, int &j){

    if(A[i][i] == A[j][j]){
        return M_PI/4;
    }else{
        return atan(  (2 * A[i][j]) / (A[i][i] - A[j][j])  ) / 2; 
    }

}

/* returns a copy transposed of matrix A */
Mat transposeMatrix(Mat &A){
    int rows = A.size();
    int columns = A[0].size();
    Mat result(rows, Vec(columns, 0.0));

    for(size_t i=0; i<rows; ++i){
        result[i][i] = A[i][i];
        for(size_t j=i+1; j<columns; ++j){
            result[i][j] = A[j][i];
            result[j][i] = A[i][j];
        }
    }

    return result;

}

/* creates an identity matrix based on give order */
Mat basicIdentity(int order){
    Mat I(order, Vec(order, 0.0));

    for(size_t i=0; i<order; ++i){
        I[i][i] = 1;
    }

    return I;
}



/* solvig for Y in LY = B (Y being equal to UX) */
void forwardSubstitution(Mat &L, Vec &B, Vec &Y, int &order){
    std::cout << "forward substitution..." << std::endl;

    Y[0] = B[0] / L[0][0];
    for(int i=1; i<order; ++i){
        Y[i] = B[i];

        // only the lower triangle
        for(int j=0; j<=i-1; ++j){
            Y[i] -= L[i][j] * Y[j];
        }

        Y[i] /= L[i][i];
    }

    std::cout << "matrix Y: " << std::endl;
    for(size_t i=0; i<order; ++i){
        std::cout << ", " << std::setprecision(9) << Y[i];
    }
    std::cout << std::endl;

}

/* solving for X in UX = Y (assuming Y has already been calculated) */
void backwardSubstitution(Mat &U, Vec &Y, Vec &X, int &order){
    std::cout << "backward substitution..." << std::endl;

    X[order-1] = Y[order-1] / U[order-1][order-1];
    for(int i=order-2; i>=0; --i){
        X[i] = Y[i];

        // only the upper triangle
        for(int j=i+1; j<=order; ++j){
            X[i] -= U[i][j] * X[j];
        }

        X[i] /= U[i][i];   
    }

}


#endif
