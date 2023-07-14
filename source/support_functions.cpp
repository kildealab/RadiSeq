#include "support_functions.h"
#include "random_generator.h"

#include <iostream>
#include <string>
#include <sstream>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

//------------------------------------------------------------------------------------------------------
// This function is used to print the ASCII-art name of the program in the beginning of the output. (Font:Roman)
//------------------------------------------------------------------------------------------------------
void ascii_art(){
std::cout<<"|---------------------------------------------------------------------------------------------|\n"
         <<"|                                                                                             |\n"
         <<"|   .oooooo..o oooooooooo.   oooooooooo.           ooooo      ooo   .oooooo.     .oooooo..o   |\n"
         <<"|  d8P'    `Y8 `888'   `Y8b  `888'   `Y8b          `888b.     `8'  d8P'  `Y8b   d8P'    `Y8   |\n"
         <<"|  Y88bo.       888      888  888      888          8 `88b.    8  888           Y88bo.        |\n"
         <<"|   `\"Y8888o.   888      888  888      888          8   `88b.  8  888            `\"Y8888o.    |\n"
         <<"|       `\"Y88b  888      888  888      888 8888888  8     `88b.8  888     ooooo      `\"Y88b   |\n"
         <<"|  oo     .d8P  888     d88'  888     d88'          8       `888  `88.    .88'  oo     .d8P   |\n"
         <<"|  8\"\"88888P'  o888bood8P'   o888bood8P'           o8o        `8   `Y8bood8P'   8\"\"88888P'    |\n"
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
bool is_nums_in_same_interval(const std::vector<long>& vec, int num1, int num2) {
    for (int i = 0; i<vec.size(); i++) {
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
void sortNremoveDuplicates_inVector(std::vector<long>& vec) {
    std::sort(vec.begin(), vec.end());                                      // Sort the vector

    auto last = std::unique(vec.begin(), vec.end());                        // Remove adjacent duplicates
    vec.erase(last, vec.end());                                             // Erase the duplicates
}
//------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------
// This function will remove a non-empty directory and it's contents
//------------------------------------------------------------------------------------------------------
void remove_directory(const std::string& path) {
    DIR* dir = opendir(path.c_str());
    if (dir == nullptr) {                                                   // Directory does not exist or cannot be opened
        return;
    }

    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        std::string filename = entry->d_name;
        if (filename != "." && filename != "..") {
            std::string filepath = path + "/" + filename;
            if (entry->d_type == DT_DIR) {
                remove_directory(filepath);                                 // Recursively remove subdirectories
            } else {
                remove(filepath.c_str());                                   // Remove regular files
            }
        }
    }

    closedir(dir);
    rmdir(path.c_str());                                                    // Remove the directory itself
}
//------------------------------------------------------------------------------------------------------