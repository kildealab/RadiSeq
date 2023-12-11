#include "support_functions.h"
#include "random_generator.h"
#include "summary_report.h"

#include <iostream>
#include <string>
#include <sstream>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <cstdio>
#include <algorithm>
#include <zlib.h>
#include <stdexcept>


//------------------------------------------------------------------------------------------------------
// This function is used to print the ASCII-art name of the program in the beginning of the output. (Font:Roman)
//------------------------------------------------------------------------------------------------------
void ascii_art(){
std::cout<<"|---------------------------------------------------------------------------------------------|\n"
         <<"|                                                                                             |\n"
         <<"|                                    .o8   o8o   .oooooo..o                                   |\n"
         <<"|                                   \"888   `\"'  d8P'    `Y8                                   |\n"
         <<"|           oooo d8b  .oooo.    .oooo888  oooo  Y88bo.       .ooooo.   .ooooo oo              |\n"
         <<"|           `888\"\"8P `P  )88b  d88' `888  `888   `\"Y8888o.  d88' `88b d88' `888               |\n"
         <<"|            888      .oP\"888  888   888   888       `\"Y88b 888ooo888 888   888               |\n"
         <<"|            888     d8(  888  888   888   888  oo     .d8P 888    .o 888   888               |\n"
         <<"|           d888b    `Y888\"\"8o `Y8bod88P\" o888o 8\"\"88888P'  `Y8bod8P' `V8bod888               |\n"
         <<"|                                                                           888.              |\n"
         <<"|                                                                           8P'               |\n"
         <<"|                                                                           \"                 |\n"
         <<"|                                                                                             |\n"
         <<"|   .oooooo..o  o8o                                oooo                .                      |\n"
         <<"|  d8P'    `Y8  `\"'                                `888              .o8                      |\n"
         <<"|  Y88bo.      oooo  ooo. .oo.  .oo.   oooo  oooo   888   .oooo.   .o888oo  .ooooo.  oooo d8b |\n"
         <<"|   `\"Y8888o.  `888  `888P\"Y88bP\"Y88b  `888  `888   888  `P  )88b    888   d88' `88b `888\"\"8P |\n"
         <<"|       `\"Y88b  888   888   888   888   888   888   888   .oP\"888    888   888   888  888     |\n"
         <<"|  oo     .d8P  888   888   888   888   888   888   888  d8(  888    888 . 888   888  888     |\n"
         <<"|  8\"\"88888P'  o888o o888o o888o o888o  `V88V\"V8P' o888o `Y888\"\"8o   \"888\" `Y8bod8P' d888b    |\n"
         <<"|                                                                                             |\n"
         <<"|---------------------------------------------------------------------------------------------|\n";
}
//------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------
// This function will check if the argument passed to the main function is appropriate. One and only 
// argument corresponding to the path to the user-defined parameter file is expected. If there are
// no arguments, more than one arguments, or an incorrect filepath: print error and exit. Else continue
//------------------------------------------------------------------------------------------------------
void checkArgument(int argc, char** argv){
    if (argc < 2){
        std::cerr<< "\n ERROR: A Parameter file is expected as an argument\n";
        exit(EXIT_FAILURE);
    }
    else if(argc > 2){
        std::cerr<< "\n ERROR: More than one argument is given\n"
                 << " (Expected: "<<argv[0]<<" path_to_ParameterFile.txt)\n";
        exit(EXIT_FAILURE);
    }
    else{
        std::cout<< "\n Parameter file specified is : "<<argv[1]<<'\n';
    } 
}
//------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------
// This function will convert any string to lowecase
//------------------------------------------------------------------------------------------------------
std::string lowercaseString(std::string* Pstr){
    std::string str{*Pstr};                                                 // Temporary variable to hold a string
    std::transform(str.begin(), str.end(), str.begin(),                                             
                    [](unsigned char c){return std::tolower(c);});
    return(str);
}
//------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------
// This function will convert any string referenced to uppercase
//------------------------------------------------------------------------------------------------------
void uppercaseString(std::string& str){
    std::transform(str.begin(), str.end(), str.begin(),                                             
                    [](unsigned char c){return std::toupper(c);});
}
//------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------
// This function can be used to convert any string that has dilimiter seperated tokens into a vector
// eg: "This,is,an,example" is a seperatedString with 4 tokens delimited with ','. When passed to this
// function, we get a vector ["This","is","an","example"]
// We are using template specialization to perform different tasks depending on the type of the vector
// so that we can convert "0.5, 0.6" --> [0.5,0.6] as well. 
// If the same vector is passed twice, the second set of string values will be overwritten to the vector
//------------------------------------------------------------------------------------------------------
template<>  // template to follow if type of the vector is string
void stringToVec(char delimiter, std::string* seperatedString, std::vector<std::string>& vec){
    std::istringstream ss(*seperatedString);
    std::string token;
    if(!vec.empty()){                                                       // if the same vector is passed twice (not empty)
        vec.clear();                                                        // empty the vector to hold new data 
    }                                                                       // i.e, (overriding default with user-given values)
    while (std::getline(ss, token, delimiter)){
        token.erase(0, token.find_first_not_of(" "));                       // removing leading whitespaces from token
        token.erase(token.find_last_not_of(" ")+1);                         // removing trailing whitespaces from token
        vec.push_back(token);                                               // a string is pushed back
    }
}
template<>  // template to follow if type of the vector is double
void stringToVec(char delimiter, std::string* seperatedString, std::vector<double>& vec){
    std::istringstream ss(*seperatedString);
    std::string token;
    if(!vec.empty()){                                                       // if the same vector is passed twice (not empty)
        vec.clear();                                                        // empty the vector to hold new data 
    }                                                                       // i.e, (overriding default with user-given values)
    while (std::getline(ss, token, delimiter)){
        token.erase(0, token.find_first_not_of(" "));                       // removing leading whitespaces from token
        token.erase(token.find_last_not_of(" ")+1);                         // removing trailing whitespaces from token
        vec.push_back(std::stod(token));                                    // a double is pushed back
    }
}
template<>  // template to follow if type of the vector is int
void stringToVec(char delimiter, std::string* seperatedString, std::vector<int>& vec){
    std::istringstream ss(*seperatedString);
    std::string token;
    if(!vec.empty()){                                                       // if the same vector is passed twice (not empty)
        vec.clear();                                                        // empty the vector to hold new data 
    }                                                                       // i.e, (overriding default with user-given values)
    while (std::getline(ss, token, delimiter)){
        token.erase(0, token.find_first_not_of(" "));                       // removing leading whitespaces from token
        token.erase(token.find_last_not_of(" ")+1);                         // removing trailing whitespaces from token
        vec.push_back(std::stoi(token));                                    // an int is pushed back
    }
}
//------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------
// This function will randomly remove an element from the vector passed. One an only one element will
// be removed from the vector
//------------------------------------------------------------------------------------------------------
template<> // Template for integer vector
void removeAnElement(std::vector<int>& vec){
    if(vec.empty()){
        std::cerr<<"\n ERROR: Vector passed to remove an element is empty\n";
        return;
    }
    int indexToRemove = rng::rand_int(0, (vec.size()-1));                   // Random index to remove. Index of a vector of size 10 goes from 0-9
    vec.erase(vec.begin() + indexToRemove);
}
template<> // Template for long vector 
void removeAnElement(std::vector<long>& vec){
    if(vec.empty()){
        std::cerr<<"\n ERROR: Vector passed to remove an element is empty\n";
        return;
    }
    int indexToRemove = rng::rand_int(0, (vec.size()-1));
    vec.erase(vec.begin() + indexToRemove);
}
template<> // Template for string vector 
void removeAnElement(std::vector<std::string>& vec){
    if(vec.empty()){return;}if(vec.empty()){
        std::cerr<<"\n ERROR: Vector passed to remove an element is empty\n";
        return;
    }
    int indexToRemove = rng::rand_int(0, (vec.size()-1));
    vec.erase(vec.begin() + indexToRemove);
}
//------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------
// This function will check if the two integer numbers passed are between two same consecutive elements
// of the vector provided. i.e, if vec=[10,20,30,40], and num1=12 and num2=15, since both these numbers
// are between 10 and 20 (two same consecutive elements), function will return True. 
//------------------------------------------------------------------------------------------------------
bool is_nums_in_same_interval(const std::vector<long>& vec, int num1, int num2){
    for (size_t i = 0; i<vec.size(); i++){
        if((i+1)<vec.size()){                                               // Proceed only if it is not the last element in the vector
            if ((num1>vec[i] && num1<vec[i+1]) && (num2>vec[i] && num2<vec[i+1])){
                return true;                                                // True if num1 and num2 are between same consecutive elements
            }
        }
    }
    return false;
}
//------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------
// This function will first sort the vector passed and then remove any duplicates if there is any
//------------------------------------------------------------------------------------------------------
void sortNremoveDuplicates_inVector(std::vector<long>& vec){
    std::sort(vec.begin(), vec.end());                                      // Sort the vector

    auto last = std::unique(vec.begin(), vec.end());                        // Remove adjacent duplicates
    vec.erase(last, vec.end());                                             // Erase the duplicates
}
//------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------
// This function will remove a non-empty directory and it's contents
//------------------------------------------------------------------------------------------------------
void remove_directory(const std::string& path){
    DIR* dir = opendir(path.c_str());
    if (dir == nullptr){                                                    // Directory does not exist or cannot be opened
        return;
    }

    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr){
        std::string filename = entry->d_name;
        if (filename != "." && filename != ".."){
            std::string filepath = path + "/" + filename;
            if (entry->d_type == DT_DIR){
                remove_directory(filepath);                                 // Recursively remove subdirectories
            } else{
                remove(filepath.c_str());                                   // Remove regular files
            }
        }
    }

    closedir(dir);
    rmdir(path.c_str());                                                    // Remove the directory itself
}
//------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------
// This function will check if the disk space available is enough to store all the fasta files and the
// fastq files produced during the run. If the disk space is not enough, program will exit with an error.
//------------------------------------------------------------------------------------------------------
void checkStorageSize(NGSParameters& parameter, NGSsdd& SDDdata, long one_fasta_size){
    long total_fasta_size = one_fasta_size * (1+SDDdata.get_num_of_exposures());                        // Total size for all the fasta files (damaged+undamaged cells)
    int size_readOfLength_one = 32;                                                                     // Estimated size needed for a fastq.gz file with one read of size one
    long one_read_size = size_readOfLength_one + (2*(parameter.get_read_length()-1));                   // For every new base, 2 bytes will get added to the size needed for one base
    long total_reads = (parameter.get_total_read_coverage()*report_ref_seq_length)/parameter.get_read_length();
    long total_fastq_size = one_read_size * total_reads;                                                // Total size for all the fastq files (all read data)
    long long total_storageSpace_needed = total_fasta_size + total_fastq_size;                          // Total space needed for storage (bytes)
    std::string formattedStorageNeeded = formatBytes(total_storageSpace_needed);                        // Convert bytes into human readable format
    
    const std::string& directoryString = *parameter.get_output_directory();                             // Output storage directory
    const char* directory = directoryString.c_str();
    char command[256];                                                                                  // Character array to store the command to be executed
    snprintf(command, sizeof(command), "df -k \"%s\" | awk 'NR==2 {print $4}'", directory);             // Command string; get the disk space available for the directory (in Kilobytes)

    FILE* pipe = popen(command, "r");                                                                   // Execute the command
    if (!pipe || !total_fasta_size){                                                                    // If we can't execute the command or if the fasta size was zero, then exit
        std::cerr<<"\n WARNING: Unable to determine the disk space available for storage \n"
                 <<" Estimated storage space needed to run the program with the given parameters is "<<formattedStorageNeeded<<"\n"
                 <<" If the available space is less than the estimate, the program might stop with an error at the end of the run\n";
        return;
    }

    char buffer[128];                                                                                   // Character array to store the command output read from the pipe
    std::string result = "";
    if (fgets(buffer, 128, pipe) != nullptr)
        result = buffer;

    pclose(pipe);                                                                                       // Close pipe 

    if (result.empty()){                                                                                // If failed to determine the disk space, return
        std::cerr<<"\n WARNING: Unable to determine the disk space available for storage \n"
                 <<" Estimated storage space needed to run the program with the given parameters is "<<formattedStorageNeeded<<"\n"
                 <<" If the available space is less than the estimate, the program might stop with an error at the end of the run\n";
        return;
    }

    long long spaceAvailable = std::stoll(result)*1000;                                                 // Disk space available (bytes)
    std::string formattedSpaceAvailable = formatBytes(spaceAvailable);                                  // Convert bytes into human readable format
    
    if (total_storageSpace_needed>spaceAvailable){                                                      // If available space is not enough 
        std::cerr<<"\n ERROR: Inadequate storage space ("<<formattedSpaceAvailable<<")\n"
                 <<" Estimated storage space needed to run the program with the given parameters is "<<formattedStorageNeeded<<"\n";
        exit(EXIT_FAILURE);
    }
    return ;
}
//------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------
// Function takes the number of bytes as input and returns a string representing the formatted size.
// The function calculates the appropriate suffix (e.g., B, KB, MB, etc.) based on powers of 1024
//------------------------------------------------------------------------------------------------------
std::string formatBytes(long long bytes){
    const int base = 1024;
    const char* suffixes[] = {"B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"};
    if (bytes == 0){
        return "0 B";                                                                                   // Handle special case of 0 bytes
    }
    int index = static_cast<int>(std::log(bytes) / std::log(base));
    double value = bytes / std::pow(base, index);

    std::ostringstream oss;
    oss << std::fixed << value;
    std::string formattedValue = oss.str();

    formattedValue.erase(formattedValue.find_last_not_of('0') + 1, std::string::npos);                  // Trim trailing zeros
    if (formattedValue.back() == '.'){
        formattedValue.pop_back();                                                                      // Remove decimal point if no fractional part
    }

    return formattedValue + " " + suffixes[index];
}
//------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------
// Function takes string data and return compressed string instead. This function is needed to create a 
// gzip text file
//------------------------------------------------------------------------------------------------------
void compressStringData(const std::string& input, std::string& output){
    const size_t CHUNK_SIZE = 16384;                                                                    // Chunk size for compressing data
    z_stream strm;                                                                                      // Temporary compression stream object to hold the state information for the compression    
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;

    if (deflateInit2(&strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 15 | 16, 8, Z_DEFAULT_STRATEGY) != Z_OK){// Initiate compression stream with necessary flags
        throw std::runtime_error("Failed to initialize zlib");
    }

    strm.next_in = (Bytef*)input.data();                                                                // Setting the input data
    strm.avail_in = input.size();                                                                       // Setting the size of the input data

    do{                                                                                                 // This is the loop that performs the compression
        char out[CHUNK_SIZE];
        strm.next_out = reinterpret_cast<Bytef*>(out);
        strm.avail_out = CHUNK_SIZE;

        if (deflate(&strm, Z_FINISH) == Z_STREAM_ERROR){
            deflateEnd(&strm);
            throw std::runtime_error("Failed to compress string");
        }

        size_t compressedBytes = CHUNK_SIZE - strm.avail_out;
        if (compressedBytes > 0) {
            output.append(out, compressedBytes);
        }
    } while (strm.avail_out == 0);                                                                      // The loop continues until strm.avail_out becomes non-zero, indicating that there is still space in the output buffer

    deflateEnd(&strm);                                                                                  // Release resources associated with the compression stream
}
//------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------
// Function takes the arguments: beta value, minimum and maximum bounds and the mode value of a beta distribution
// you want to create and generate a vector containing normalized probaility densities of each value between the bounds.
// The returned weight vector will have probabilities for each integer value between minimum and maximum bounds in order.
// eg: if the bound is [50,150], then generated vector will have 100 elements corresponding to the probabilities
// for 51,52,23...150 according to the beta distribution described by the beta value and mode.
//------------------------------------------------------------------------------------------------------
std::vector<double> beta_distribution_proabalities(double beta, int minSize, int maxSize, double realMode){
    int steps = maxSize - minSize;                                                                      // Determine how many integer steps are between max and min values
    double stepSize = 1.0/steps;                                                                        // The step size for these many steps in the interval [0,1] is range/steps
    double betaMode = (realMode-minSize)*stepSize;                                                      // Determine where the mode should be in the range [0,1] if actual mode in range [min,max] is known
    double alpha =  ((-betaMode*beta)+(2*betaMode)-1)/(betaMode-1);                                     // For a beta distribution Mode = (a-1)/(a+b-2). This equation can be used to find alpha (a)
    
    std::vector<double> weights(steps);                                                                 // Temporary vector to store the normalized probabilities for each integer value corresponding to each step
    double multiplier = tgamma(alpha+beta)/(tgamma(alpha)*tgamma(beta));                                // This is the multiplier constant in the beta probability density function
    double totalWeight{0};                                                                              // Temporary value to hold the sum of of probability densities to use for normalization later
    for(int i=0; i<steps; i++){
        double x = i*stepSize;                                                                          // Determine the value each x values (steps) in range [0,1] with the given step size
        weights[i] = multiplier * pow(x,(alpha-1)) * pow((1-x),(beta-1));                               // Beta PDF:  f(x) = Const * x^(q-1) * (1-x)^(b-1)
        totalWeight += weights[i];                                                                      // Get the total
    }
    
    for (double& weight : weights){                                                                     // Normalize the Beta PDF to have the integral under the curve to be 1
        weight/=totalWeight;                                        
    }
    return (weights);                                                                                   // Return the normalized probability vector 
}
//------------------------------------------------------------------------------------------------------