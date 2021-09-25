#include "subRoutines.hpp"

iterativeResult result_gaussseidel;

void iterativeGaussSeidel(Mat &A, Vec &B, Vec &X, int &order, double &precision){
    double residue = 1e9;
    
    Mat L(order, Vec(order, 0.0)),  U(order, Vec(order, 0.0)),  L_inv(order, Vec(order,0));
    Vec V(order, 0.0), diff_vec(order, 0.0);

    // initialize X (initial solution)
    for(size_t i=0; i<order; ++i){
        X[i] = 1;
    }

    while(residue > precision){
        //std::cout << "gauss-seidel residue: " << residue << std::endl;
        for(int i=0; i<order; ++i){
            // defining new solution
            V[i] = B[i];
        
            for(int j=0; j<=i-1; ++j){
                V[i] -= A[i][j] * V[j];
            }

            for(int j=i+1; j<order; ++j){
                V[i] -= A[i][j] * X[j];
            }

            // new solution
            V[i] /= A[i][i];
        }

        subtractVec(V, X, diff_vec);
        // updating residue and solution
        residue = absVec(diff_vec) / absVec(V);
        X = V;
        result_gaussseidel.total_iterations++;
        result_gaussseidel.TOL_history.push_back( residue );
    }
}

/* solves for X in AX = B such that A = L* + U */
iterativeResult solverIterativeGaussSeidel(Mat &A, Vec &B, Vec &X, int &order, double &precision){
    std::cout << "Iterative Gauss-Seidel selected." << std::endl;

    // necessary condition for convergence of solution
    if(!diagonallyDominant(A, order)){
        // failure
        std::cout << "Iterative Gauss-Seidel failed. A is not diagonally dominant." << std::endl;
        return result_gaussseidel;
    }else{
        iterativeGaussSeidel(A, B, X, order, precision);
        std::cout << "Iterative Gauss-Seidel successfully done." << std::endl;
        return result_gaussseidel;
    }
}
