#include "subRoutines.hpp"

generalResult result_cholesky;
bool cholesky_success = true;

/* reading data from matrix A and generating matrix L and L* (L_t) from said data */
void CholeskyDecomposition(Mat &A, Mat &L, Mat &L_t, int &order){
    std::cout << "Setting L and L*..." << std::endl;

    for(int i=0; i<order; ++i){
        double sum = 0;
        for(int k=0; k<=i-1; ++k){
            sum += L[i][k] * L[i][k];
        }
        
        double root_element = A[i][i] - sum;
        // check if A is not positive-definite
        if(root_element <= 0){
            // failure
            cholesky_success = false;
            break;
        }

        L[i][i] = sqrt(root_element);
        L_t[i][i] = sqrt(root_element);

        for(int j=i+1; j<order; ++j){

            double summation = 0;
            for(int k=0; k<=i-1; ++k){
                summation += L[i][k] * L[j][k];
            }

            L[j][i] = (A[i][j] - summation) / (double) L[i][i];
            L_t[i][j] = (A[i][j] - summation) / (double) L[i][i];
        }

    }


    double L_det = 1;
    double L_t_det = 1;
    for(int i=0; i<order; ++i){
        L_det *= L[i][i];
        L_t_det *= L_t[i][i];
    }

    
    std::cout << "det(L): " << L_det << std::endl;
    std::cout << "det(L_t): " << L_t_det << std::endl;


}


/* Solves for X in AX = B by Cholesky (LL*) Decomposition. Returns status (1=success, 0=failure) and determinant. */
generalResult solverCholesky( Mat &A, Vec &B, Vec &X, int &order){
    std::cout << "Cholesky Decomposition selected." << std::endl;
    // creating intermediate varibles for the LU decomposition
    Vec Y(order, 0.0);
    Mat L(order, Vec(order, 0.0)), L_t(order, Vec(order, 0.0));

    CholeskyDecomposition(A, L, L_t, order);

    result_cholesky.determinant = calculateLUDeterminant(L, L_t, order);
    // check if A was not positive-definite
    if(!cholesky_success){
        // failure
        std::cout << "Cholesky Decomposition failed. A is not positive-definite." << std::endl;
        result_cholesky.errors.push_back("A is not positive-definite.");
        return result_cholesky;
    }else{
        forwardSubstitution(L, B, Y, order);
        backwardSubstitution(L_t, Y, X, order);
        std::cout << "Cholesky Decomposition successfully done." << std::endl;
        return result_cholesky;
    }
}
