#include "subRoutines.hpp"

iterativeResult result_power;

void powerMethod(Mat &A, Vec &result_values, Mat &result_vectors, int &order, double &precision){
    double residue = 1e9;

    Vec eigenvector(order), Y(order), diff_vec(order);

    // initializing eigenvector
    for(size_t i=0; i<order; ++i){
        eigenvector[i] = 1;
    }
    double old_eigen = 1;
    
    result_values.push_back(1);
    result_vectors.push_back(eigenvector);

    while(residue > precision){
        
        // Y = AX
        productMatrixVector(A, eigenvector, Y);
        
        // normalize Y
        for(int i=1; i<order; ++i){
            Y[i] /= Y[0];
        }
        result_values.push_back(Y[0]);

        residue = fabs(Y[0] - old_eigen) /  fabs(Y[0]);
        old_eigen = Y[0];

        // divide by itself
        Y[0] = 1.0;

        result_vectors.push_back(Y);
        eigenvector = Y;

        result_power.total_iterations++;
        result_power.TOL_history.push_back( residue );

    }
    
}


iterativeResult mainPowerMethod(Mat &A, Vec &result_values, Mat &result_vectors, int &order, double &precision){
    std::cout << "Power Method selected." << std::endl;

    powerMethod(A, result_values, result_vectors, order, precision);

    return result_power;
}
