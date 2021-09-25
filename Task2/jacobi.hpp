#include "subRoutines.hpp"

iterativeResult result_jacobi;

/* check if exists an element not in the diagonal with absolute value greater than or equal to tol */
bool checkNonDiagonal(Mat &A, int &order, double &tol){

    for(size_t i=0; i<order; ++i){
        for(size_t j=i+1; j<order; ++j){
            if(fabs(A[i][j]) >= tol || fabs(A[j][i]) >= tol){
                return true;
            }
        }
    }
    return false;
}

/* returns the index {row, column} of the abs max not in the diagonal */
std::vector<int> indexOfMaxNonDiagonal(Mat &A, int &order){

    double candidate = 0;
    std::vector<int> indices = {0, 0};

    for(size_t i=0; i<order; ++i){
        for(size_t j=i+1; j<order; ++j){
            
            if(fabs(A[i][j]) > candidate){
                candidate = fabs(A[i][j]);
                indices[0] = i;
                indices[1] = j;
            }
            if(fabs(A[j][i]) > candidate){
                candidate = fabs(A[j][i]);
                indices[0] = i;
                indices[1] = j;
            }

        }
    }

    return indices;

}

/* returns a Givens Matrix (P) given its special indices (i and j) and the reference matrix A */
Mat customGivensMatrix(int &order, int &i, int &j, Mat &A){

    Mat G(order, Vec(order, 0.0));
    double phi = phiAngle(A, i, j);
    double c = cos(phi);
    double s = sin(phi);

    // diagonal
    for(size_t k=0; k<order; ++k){
        G[k][k] = 1;
    }

    G[i][i] = c;
    G[j][j] = c;
    G[i][j] = -s;
    G[j][i] = s;
     
    return G;

}

void jacobiMethod(Mat &A, Vec &result_values, Mat &result_vectors, int &order, double &precision){
    //double residue = 1e9;

    Mat X = basicIdentity(order);
    Mat P(order, Vec(order));

    std::vector<int> indices = {0, 0};

    while( checkNonDiagonal(A, order, precision) ){

        // element to be zeroed
        indices = indexOfMaxNonDiagonal(A, order);

        // new P
        P = customGivensMatrix(order, indices[0], indices[1], A);

        // temp = Pt * A
        Mat temp(order, Vec(order));
        Mat P_t = transposeMatrix(P);
        productMatrices(P_t, A, temp);
        // A = temp * P
        productMatrices(temp, P, A);
        // which should be equivalent to A = Pt * A * P

        temp = X;
        productMatrices(temp, P, X);

        result_jacobi.total_iterations++;
    }

    // save diagonal
    for(size_t i=0; i<order; ++i){
        result_values.push_back(A[i][i]);
    }

    // save X
    result_vectors = X;
}

iterativeResult mainJacobiMethod(Mat &A, Vec &result_values, Mat &result_vectors, int &order, double &precision){
    std::cout << "Jacobi Method selected." << std::endl;

    if(symmetricMatrix(A)){
        jacobiMethod(A, result_values, result_vectors, order, precision);
    }else{
        result_jacobi.errors.push_back("MATRIX A IS NOT SYMMETRIC");
    }

    return result_jacobi;
}
