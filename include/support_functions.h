#ifndef SUPPORT_FUNCTIONS_H
#define SUPPORT_FUNCTIONS_H

#include <string>
#include <vector>

void ascii_art();                                                            // Function to print the ASCII art of the program name
void checkArgument(int, char**);                                             // Function to check the main() function arguments
std::string lowercaseString(std::string*);                                   // Function to convert any string to lowercase
void uppercaseString(std::string&);                                          // Function will convert any string referenced to uppercases 

// Using template specialization to make the funcition work differently depending on the type of the vector passed as an argument
template<typename T>                                                         
void stringToVec(char, std::string*, std::vector<T>&);                       // Function to convert any string elements seperated by a delimeter 
                                                                             // into a vector. eg: "this,is" --> ["this","is"]
template<typename T>
void removeAnElement(std::vector<T>&);                                       // Function to remove one element from a vector at random

bool is_nums_in_same_interval(const std::vector<long>&, int, int);           // Function is used to check if opposite SSBs are on same chromosome
void sortNremoveDuplicates_inVector(std::vector<long>&);                     // Function will sort() and remove duplicates from any long vector passed
void remove_directory(const std::string&);                                   // Function to remove a non-empty directory and its contents


#endif