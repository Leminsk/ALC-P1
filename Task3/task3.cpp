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

int N;     // total amount of coordinate pairs
int ICOD;  // method code

double estimate_x;
std::vector< std::vector<double> > coordinates;
double final_estimation;

bool success = true;
std::vector<std::string> errors;

/* INPUT_FILE.txt description:
    line 0: N total amount of coordinate pairs
    line 1: ICOD -> 1 (Langrange Interpolation), 2 (Multilinear Regression)
    line 2: input file containing said coordinates (one pair per line)
    line 3: coordinate estimate x to calculate the value of the function (y)

*/

/* coordinates.txt description:
    Each line represents a line of the matrix. Each value (or column) must be separated by a ','.
*/

/* OUTPUT_FILE.txt description:
    line 0: estimated value of y
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


/* populate coordinates vector */
void readInputCoordinates(int &order, std::string file_name){
    std::string line;
    std::ifstream file(file_name);
    std::vector<std::string> lines;
    bool correct_size = true;

    // storing lines in temporary vector
    while( std::getline(file, line) ){
        lines.push_back(line);
    }

    if(lines.size() == order){
        // iterate through lines first
        for(int line=0; line<order; ++line){
            std::vector<double> temp_vec;
            std::vector<std::string> split_string = stringSplit(lines[line], SEPARATOR);

            if( split_string.size() == 2 ){
                // build coordinate
                temp_vec.push_back( stod(split_string[0]) );
                temp_vec.push_back( stod(split_string[1]) );
                coordinates.push_back(temp_vec);
            }else{
                errors.push_back("DETECTED ERRONEOUS COORDINATE. NOT A PAIR");
                success = false;
                break;
            }

        }

    }else{
        errors.push_back("AMOUNT (N) INPUT DOES NOT MATCH AMOUNT OF COORDINATES FROM INPUT FILE");
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
    if(file_lines.size() == 4){

        N = stoi(file_lines[0]);
        ICOD = stoi(file_lines[1]);
        
        // populate coordinates
        readInputCoordinates( N, file_lines[2] );

        if(success){
            estimate_x = stod(file_lines[3]);
        }else{
            // force failure
            ICOD = -1;
        }

    }else{
        errors.push_back("INCORRECT NUMBER OF ARGUMENTS IN "+INPUT_FILE_NAME);
        // force failure
        ICOD = -1;
    }


}

/* Write the necessary variables to a .txt. Write only errors to output if any exist. */
void writeOutput(std::string output_file_name){

    std::ofstream file;
    file.open(output_file_name);

    if(success){
        file << "Estimation: " << final_estimation;
    }else{
        file << "Calculation failed. " << "ICOD: " << ICOD << "\n";
        for(size_t i=0; i<errors.size(); ++i){
            file << errors[i] << "\n";
        }
    }
    
    file.close();

}

void evaluateICOD(){
    basicResult result;
    switch(ICOD){
        case 1:
            result = mainLagrange(coordinates, estimate_x, N);
            errors = result.errors;
            final_estimation = result.estimate;
            break;
        case 2:
            break;
        default:
            std::cout << "INVALID ICOD FROM: " << INPUT_FILE_NAME << std::endl;
            errors.push_back("INVALID ICOD: "+std::to_string(ICOD));
            success = false;
    } 

}


int main(){
    
    std::cout << "============================================================================" << std::endl;

    readMainInput(INPUT_FILE_NAME);
    
    evaluateICOD();
    
    writeOutput(OUTPUT_FILE_NAME);
        
    std::cout << "============================================================================" << std::endl;
    std::cout << "Task3 main TERMINATED" << std::endl;
    std::cout << "============================================================================" << std::endl;

    // Windows only
    system("pause");

    
    return 0;
}
