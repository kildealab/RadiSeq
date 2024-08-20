#ifndef FILEIO_H
#define FILEIO_H

#include <iostream>
#include <string>
#include <map>

#include "parameter_handler.h"
#include "sddfile_handler.h"

class NGSParameters;                                                                                    // Forward declaration of the class before definition
class NGSsdd;                                                                                           // Forward declaration of the class before definition

bool checkFileExists(const std::string*);                                                               // function to check if a file exists and in a good condition. A pointer to a string (filename) is passed
bool checkFolderExists(const char*);                                                                    // function to check if a folder already exists
void readParameterFile(const std::string*, NGSParameters& );                                            // function to read a parameter file. A pointer to a string (filename) is passed
int countExposuresSDD(const std::string*, std::vector<std::streampos>&);                                // function to count the number of exposures in every SDD file and store the lines corresponding to new exposures 
std::vector<double> readActualDosefile(const std::string*, int);                                        // function to read a file containing actual dose in each exposure/irradiation
void readSDDfileHeader(const std::string*, NGSsdd&);                                                    // function to read an SDD file header and obtain values from fields 11 and 15 of the header
void readSDDfileData(const std::string*, NGSsdd&, int, int, std::vector<std::string>&, int);            // function to read an SDD file data
std::pair<int, std::vector<double>> readFragmentSizeDist(const std::string*);
char* createMemoryMappedFile(const std::string&, size_t);                                               // function to create a file with a memory-map and returns the pointer to the memory-map
char* generateInputFileMemoryMap(const std::string&, size_t&);                                          // function to generate a memory-map of an input file and returns the pointer to the memory-map
void resizeMemoryMappedFile(char*&, int&, size_t&, size_t);                                             // function to re-map an existing memory mapped file when data exceeds its allotted size 
int getNextChromSeq(std::ifstream&, std::string&, std::string&);                                        // function to read the reference sequence file, one chromosome sequence at a time
int getNextChromSeq_MM(const char*, size_t, size_t&, std::string&, std::string&);                       // function to read the reference sequence memory-map, one chromosome sequence at a time
void writeBatchToFile(std::vector<std::string>&, std::ofstream&, bool);                                 // function to write to files in batch from a buffer 
void writeBatchToMMFile(std::vector<std::string>&, char*&, char*&, size_t&, std::string filePath="NA"); // function to write to memory mapped Files in batch from a buffer. filePath is an optional argument 
int readFastaTemplate(std::ifstream&, std::string&, std::string&, std::string&, std::string&);          // function to read the fasta file template file, one chrom seq and chrom ID at a time 
int readFastaMemoryMap(const char*, size_t, size_t&, std::string&, std::string&, std::string&, std::string&);
void make_quality_distribution(std::ifstream&, std::vector<std::map<unsigned int, unsigned short>>&);   // function to make a read quaity distribution vector based on the quality profile file passed to the function
void generate_run_summaryReport(NGSParameters&, NGSsdd& );
long fileSize_bytes(const std::string&);                                                                // Function to determine the size of a file that is passed. Return in bytes

#endif