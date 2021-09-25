#include "subRoutines.hpp"

generalResult result_lu;

/* reading data from matrix A and generating triangle matrices L and U from said data */
void LUDecomposition(Mat &A, Mat &L, Mat &U, int &order){
    std::cout << "Setting L and U..." << std::endl;
    // initialization step  
    // line index
    for(int i=0; i<order; ++i){
        // column index
        for (int j=0; j<order; ++j){
            // U
            if( i==0 ){ 
                U[i][j] = A[i][j]; 
            }else { 
                U[i][j] = 0;
            }
            // L
            if( j==0 ){ 
                L[i][j] = A[i][j] / U[0][j];
            }else if( i==j ){ 
                L[i][j] = 1; 
            }else { 
                L[i][j] = 0;
            }
        }
    }

    double sum;
    // line index
    for(int i=1; i<order; ++i){
        // column index
        for(int j=1; j<order; ++j){

            sum = 0;
            for(int k=0; k<j; ++k){
                sum += L[i][k] * U[k][j]; 
            }

            if( i>j ){
                L[i][j] = (A[i][j] - sum) / U[j][j];
            }
            else{
                U[i][j] = A[i][j] - sum;
            }
        }   
    }

    double L_det = 1;
    double U_det = 1;
    for(int i=0; i<order; ++i){
        L_det *= L[i][i];
        U_det *= U[i][i];
    }

    std::cout << "det(L): " << L_det << std::endl;
    std::cout << "det(U): " << U_det << std::endl;

}

/* Solves for X in AX = B by LU Decomposition. Returns status (1=success, 0=failure) and determinant. */
generalResult solverLU( Mat &A, Vec &B, Vec &X, int &order){
    std::cout << "LU Decomposition selected." << std::endl;
    // creating intermediate varibles for the LU decomposition
    Vec Y(order, 0.0);
    Mat L(order, Vec(order, 0.0)), U(order, Vec(order, 0.0));

    LUDecomposition(A, L, U, order);
    
    result_lu.determinant = calculateLUDeterminant(L, U, order);
    // check for degenerate matrix
    if(result_lu.determinant == 0){
        // failure
        std::cout << "LU Decomposition failed. A is degenerate/singular." << std::endl;
        result_lu.errors.push_back("A is degenerate/singular.");
        return result_lu;
    }else{
        forwardSubstitution(L, B, Y, order);
        backwardSubstitution(U, Y, X, order);
        std::cout << "LU Decomposition successfully done." << std::endl;
        return result_lu;
    }

    
}
