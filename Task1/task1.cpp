#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>
#include "routines.hpp"

const std::string INPUT_FILE_NAME = "main_input.txt";
const std::string OUTPUT_FILE_NAME = "main_output.txt";
const std::string SEPARATOR = ",";

int N;     // matrix order
int ICOD;  // method code
bool IDET; // determinant of matrix A calculation flag

// AX = B
std::vector< std::vector<double> > matrix_A;
std::vector<double> matrix_X; // unknown solution
std::vector<double> matrix_B;
double TOLm;            // maximum tolerance for iterative solution

// initialize determinant as nan
double A_determinant = std::numeric_limits<double>::quiet_NaN();


bool output_determinant = false;
bool success = true;
std::vector<double> TOL_history;
int total_iteratios = 0;
std::vector<std::string> errors;

/* INPUT_FILE.txt description:
    line 0: N order if the system of equations
    line 1: ICOD -> 1 (LU Decomposition), 2 (Cholesky Decomposition), 3 (Iterative Jacobi), 4 (Iterative Gauss-Seidel)
    line 2: IDET -> 0 (no determinant), >0 (calculate determinant)
    line 3: file name for matrix A
    line 4: file name for matrix B
    line 5: TOLm (maximum tolerance for an iterative solution)

*/

/* matrix_A.txt description:
    Each line represents a line of the matrix. Each value (or column) must be separated by a comma ', '.
*/

/* matrix_B.txt description:
    Since B is a single column matrix, each line represents a value of the column.
*/

/* OUTPUT_FILE.txt description:
    line 0: solution of system X starts here
    line n: possible errors
    line n+k*: determinant (if prompted)
    line n+k*+j: number of iterations for convergence; history of variation of TOL
    
*/

/* separates an input_string using its delimiter and returns a new string vector containing the elements inbetween each occurance of the delimiter */
std::vector<std::string> stringSplit(std::string input_string, std::string delimiter){

    std::string input_copy = input_string;
    size_t pos = 0;
    std::string token;
    std::vector<std::string> result;

    while( (pos = input_copy.find(delimiter)) != std::string::npos ){
        result.push_back( input_copy.substr(0, pos) );
        input_copy.erase(0, pos + delimiter.length());
    }
    result.push_back(input_copy);

    return result;
}


/* populate matrix_A */
void readInputA(int order, std::string file_name){
    std::string line;
    std::ifstream file(file_name);
    std::vector<std::string> lines;
    bool correct_size = true;

    // storing lines in temporary vector
    while( std::getline(file, line) ){
        lines.push_back(line);
    }

    if(lines.size() == N){

        // iterate through lines first
        for(int line=0; line<order; ++line){
            std::vector<double> temp_vec;
            std::vector<std::string> split_string = stringSplit(lines[line], SEPARATOR);

            // build row
            for(int column=0; column<order; ++column){
                temp_vec.push_back( stod(split_string[column]) );
            }

            // add new row to the matrix
            matrix_A.push_back(temp_vec);
        }

    }else{
        errors.push_back("ORDER (N) INPUT DOES NOT MATCH MATRIX A'S SIZE");
        success = false;
    }

}

/* populate vec representing matrix_B */
void readInputB(int order, std::string file_name){
    std::string line;
    std::ifstream file(file_name);
    std::vector<std::string> lines;
    bool correct_size = true;

    // storing lines in temporary vector
    while( std::getline(file, line) ){
        lines.push_back(line);
    }

    if(lines.size() == N){
        // add by line
        for(int line=0; line<order; ++line){
            matrix_B.push_back( stod(lines[line]) );
        }
    }else{
        errors.push_back("ORDER (N) INPUT DOES NOT MATCH VECTOR B'S SIZE");
        success = false;
    }
    

}

/* read from txt, populate variables and set global flags */
void readMainInput(std::string input_file_name){
    std::string line;
    std::ifstream file(input_file_name);
    std::vector<std::string> file_lines;

    // save main input lines to string vector
    while( std::getline(file, line) ){
        file_lines.push_back(line);
    }

    // check for all arguments
    if(file_lines.size() == 6){

        N = stoi(file_lines[0]);

        ICOD = stoi(file_lines[1]);
        if( stoi(file_lines[2]) > 0 ){
            IDET = true;
        }else{
            IDET = false;
        }
        
        // populate variables matrix_A, matrix_B and matrix_X
        readInputA( N, file_lines[3] );
        readInputB( N, file_lines[4] );
        if(success){
            matrix_X.resize(N, 0);
            TOLm = stod(file_lines[5]);
        }else{
            // force failure
            ICOD = -1;
            IDET= false;
        }

    }else{
        errors.push_back("INCORRECT NUMBER OF ARGUMENTS IN "+INPUT_FILE_NAME);
        // force failure
        ICOD = -1;
        IDET = false;
    }

    
    

}

/* Write the necessary variables to a .txt. Will write only errors to output if any exist. */
void writeOutput(std::string output_file_name){

    std::ofstream file;
    file.open(output_file_name);

    if(success){

        file << "X: ";
        for(size_t i=0; i<N-1; ++i){
            file << matrix_X[i] << ", ";
        }
        file << matrix_X[N-1] << "\n";
      
        if(output_determinant){
            file << "Determinant of A: " << A_determinant << "\n";
        }

        if(total_iteratios > 0){
            file << "Total Iterations: " << total_iteratios << "\n";
            file << "TOL history: " << "\n";
            for(size_t i=0; i<TOL_history.size(); ++i){
                file << TOL_history[i] << "\n";
            }
        }


    }else{
        file << "Calculation failed. " << "ICOD: " << ICOD << "\n";
        for(size_t i=0; i<errors.size(); ++i){
            file << errors[i] << "\n";
        }
    }
    
    file.close();

}

void evaluateICOD(){
    generalResult gen_result;
    iterativeResult it_result;
    switch(ICOD){
        case 1:
            gen_result = solverLU(matrix_A, matrix_B, matrix_X, N);
            if(!gen_result.errors.empty()){
                success = false;
                errors = gen_result.errors;
                break;
            }
            A_determinant = gen_result.determinant;
            break;
        case 2:
            gen_result = solverCholesky(matrix_A, matrix_B, matrix_X, N);
            if(!gen_result.errors.empty()){
                success = false;
                errors = gen_result.errors;
                break;
            }
            A_determinant = gen_result.determinant;
            break;
        case 3:
            it_result = solverIterativeJacobi(matrix_A, matrix_B, matrix_X, N, TOLm);
            if(!it_result.errors.empty()){
                success = false;
                errors = it_result.errors;
                break;
            }
            TOL_history = it_result.TOL_history;
            total_iteratios = it_result.total_iterations;
            break;
        case 4:
            it_result = solverIterativeGaussSeidel(matrix_A, matrix_B, matrix_X, N, TOLm);
            if(!it_result.errors.empty()){
                success = false;
                errors = it_result.errors;
                break;
            }
            TOL_history = it_result.TOL_history;
            total_iteratios = it_result.total_iterations;
            break;
        default:
            std::cout << "INVALID ICOD FROM: " << INPUT_FILE_NAME << std::endl;
            errors.push_back("INVALID ICOD: "+std::to_string(ICOD));
            success = false;
    } 

}

void evaluateIDET(){
    if(IDET){
        if(std::isnan(A_determinant) == 1){
            A_determinant = determinantOfMatrix(matrix_A, N);
        }

        // write determinant to output file
        output_determinant = true;
        std::cout << "A_determinant: " << A_determinant << std::endl;
    }
}

int main(){
    
    std::cout << "============================================================================" << std::endl;

    readMainInput(INPUT_FILE_NAME);
    
    evaluateICOD();
    evaluateIDET();
    
    writeOutput(OUTPUT_FILE_NAME);
    
    std::cout << "============================================================================" << std::endl;
    std::cout << "Task1 main TERMINATED" << std::endl;
    std::cout << "============================================================================" << std::endl;

    // Windows only
    system("pause");
    
    return 0;
}
