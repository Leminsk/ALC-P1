#include "subRoutines.hpp"

basicResult result_lagrange;

double lagrangeInterpolation(Vec &X, Vec &Y, double &p_x, int &amount){
    double solution = 0;
    double phi;
    for(size_t i=0; i<amount; ++i){
        phi = 1;
        for(size_t j=0; j<amount; ++j){
            if( i!=j ){
                if( (X[i] - X[j]) != 0 ){
                    phi *= (p_x - X[j]) / (X[i] - X[j]);
                }else{
                    result_lagrange.errors.push_back("CANNOT ESTIMATE: TWO PAIRS HAVE THE SAME X VALUE.");
                }
            }
        }
        solution = solution + Y[i]*phi;
    }

    return solution;
}

/* given a set of coordinates (pairs), estimate the function value of p_x by Lagrange Interpolation */
basicResult mainLagrange(Mat &pairs, double &p_x, int &amount){
    std::cout << "Lagrange Interpolation selected." << std::endl;
    Vec X(amount), Y(amount);

    for(size_t i=0; i<amount; ++i){
        X[i] = pairs[i][0];
        Y[i] = pairs[i][1];
    }

    result_lagrange.estimate = lagrangeInterpolation(X, Y, p_x, amount);
    return result_lagrange;
}
