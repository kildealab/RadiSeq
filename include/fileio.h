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
int countExposuresSDD(const std::string*);                                                              // function to count the number of exposures in every SDD file. 
std::vector<double> readActualDosefile(const std::string*, int);                                        // function to read a file containing actual dose in each exposure/irradiation
void readSDDfileHeader(const std::string*, NGSsdd&);                                                    // function to read an SDD file header and obtain values from fields 11 and 15 of the header
void readSDDfileData(const std::string*, NGSsdd&, int);                                                 // function to read an SDD file data
int getNextChromSeq(std::ifstream&, std::string&, std::string&);                                        // function to read the reference sequence file, one chromosome sequence at a time
int readFastaTemplate(std::ifstream&, std::string&, std::string&, std::string&, std::string&);          // function to read the fasta file template file, one chrom seq and chrom ID at a time 
void make_quality_distribution(std::ifstream&, std::vector<std::map<unsigned int, unsigned short>>&);   // function to make a read quaity distribution vector based on the quality profile file passed to the function
void generate_run_summaryReport(NGSParameters&, NGSsdd& );

#endif