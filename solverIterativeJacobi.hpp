#include "subRoutines.hpp"

iterativeResult result_jacobi;

void iterativeJacobi(Mat &A, Vec &B, Vec &X, int &order, double &precision){
    double residue = 1e9;
    Vec V(order, 0.0), diff_vec(order, 0.0);
    
    // initialize X (initial solution)
    for(size_t i=0; i<order; ++i){
        X[i] = 1;
    }

    while(residue > precision){
        std::cout << "jacobi residue: " << residue << std::endl;
        for(size_t i=0; i<order; ++i){
            // defining new solution
            V[i] = B[i];
            for(size_t j=0; j<order; ++j){
                if( j!=i ){
                    V[i] -= A[i][j] * X[j];
                }
            }
            // new solution
            V[i] /= A[i][i];
        }

        subtractVec(V, X, diff_vec);
        // updating residue and solution
        residue = absVec(diff_vec) / absVec(V);
        X = V;
        result_jacobi.total_iterations++;
        result_jacobi.TOL_history.push_back( residue );
    }
}

iterativeResult solverIterativeJacobi(Mat &A, Vec &B, Vec &X, int &order, double &precision){
    std::cout << "Iterative Jacobi selected." << std::endl;

    // necessary condition for convergence of solution
    if(!diagonallyDominant(A, order)){
        // failure
        std::cout << "Iterative Jacobi failed. A is not diagonally dominant." << std::endl;
        return result_jacobi;
    }else{
        iterativeJacobi(A, B, X, order, precision);
        std::cout << "Iterative Jacobi successfully done." << std::endl;
        return result_jacobi;
    }
}
