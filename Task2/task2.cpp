#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>
#include "routines.hpp"

const std::string INPUT_FILE_NAME = "main_input.txt";
const std::string OUTPUT_FILE_NAME = "main_output.txt";
const std::string SEPARATOR = ","; // separates the values of rows of the input

int N;     // matrix order
int ICOD;  // method code
bool IDET; // determinant of matrix A calculation flag

// Square Matrix A
std::vector< std::vector<double> > matrix_A;
std::vector< double > eigenvalues;               // unkown at compilation time 
std::vector< std::vector<double> > eigenvectors; // unkown at compilation time
double TOLm;            // maximum tolerance for iterative solution

// initialize determinant as nan
double A_determinant = std::numeric_limits<double>::quiet_NaN();


bool output_determinant = false;
bool success = true;
std::vector<double> TOL_history;
int total_iterations = 0;
std::vector<std::string> errors;

/* INPUT_FILE.txt description:
    line 0: N order if the system of equations
    line 1: ICOD -> 1 (Power Method), 2 (Jacobi Method)
    line 2: IDET -> 0 (no determinant), >0 (calculate determinant)
    line 3: file name for matrix A
    line 4: TOLm (maximum tolerance for an iterative solution)

*/

/* matrix_A.txt description:
    Each line represents a line of the matrix. Each value (or column) must be separated by a ' '.
*/

/* OUTPUT_FILE.txt description:
    line 0: eigenvalues
    line 1: eigenvectors start here
    line 1+k*: determinant (if prompted)
    line 1+k*+j: number of iterations for convergence; history of variation of TOL
    
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
void readInputA(int &order, std::string file_name){
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
    if(file_lines.size() == 5){

        N = stoi(file_lines[0]);

        ICOD = stoi(file_lines[1]);
        if( stoi(file_lines[2]) > 0 ){
            IDET = true;
        }else{
            IDET = false;
        }
        
        // populate variables matrix_A, matrix_B and matrix_X
        readInputA( N, file_lines[3] );

        if(success){
            TOLm = stod(file_lines[4]);
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

/* Write the necessary variables to a .txt. Write only errors to output if any exist. */
void writeOutput(std::string output_file_name){

    std::ofstream file;
    file.open(output_file_name);

    if(success){

        file << "Eigenvalues:\n";
        for(size_t i=0; i<eigenvalues.size()-1; ++i){
            file << eigenvalues[i] << ", ";
        }
        file << eigenvalues[eigenvalues.size()-1] << "\n";
        
        file << "Eigenvectors:\n";
        for(size_t i=0; i<eigenvectors.size(); ++i){
            for(size_t j=0; j<eigenvectors[i].size()-1; ++j){
                file << eigenvectors[i][j] << ", ";
            }
            file << eigenvectors[i][eigenvectors[i].size()-1] <<"\n";
        }
        
        
        if(output_determinant){
            file << "Determinant of A: " << A_determinant << "\n";
        }

        if(total_iterations > 0){
            file << "Total Iterations: " << total_iterations << "\n";
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
    iterativeResult it_result;
    switch(ICOD){
        case 1:
            it_result = mainPowerMethod(matrix_A, eigenvalues, eigenvectors, N, TOLm);
            total_iterations = it_result.total_iterations;
            TOL_history = it_result.TOL_history;
            errors = it_result.errors;
            break;
        case 2:
            it_result = mainJacobiMethod(matrix_A, eigenvalues, eigenvectors, N, TOLm);
            total_iterations = it_result.total_iterations;
            TOL_history = it_result.TOL_history;
            errors = it_result.errors;
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
    }
}

int main(){
    
    std::cout << "============================================================================" << std::endl;

    readMainInput(INPUT_FILE_NAME);
    
    evaluateICOD();
    evaluateIDET();
    
    writeOutput(OUTPUT_FILE_NAME);
        
    std::cout << "============================================================================" << std::endl;
    std::cout << "Task2 main TERMINATED" << std::endl;
    std::cout << "============================================================================" << std::endl;

    // Windows only
    system("pause");

    
    return 0;
}