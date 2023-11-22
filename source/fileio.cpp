#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <sstream>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <vector>
#include <numeric>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstring>
#include <omp.h>

#include "fileio.h"
#include "summary_report.h"

//------------------------------------------------------------------------------------------------------------------------
// A function to check if a user-specified file can be opened successfully or not
//------------------------------------------------------------------------------------------------------------------------
bool checkFileExists(const std::string* filename){
    std::ifstream file(*filename);
    if (file.good()){
        file.close();
        return true;
    } else {
        file.close();
        return false;
    }
}
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------
// A function to check if a user-specified folder exists
//------------------------------------------------------------------------------------------------------------------------
bool checkFolderExists(const char* folderPath) {
    struct stat info;
    if (stat(folderPath, &info) != 0) {                                                                 // check if any info available on the path
        return false;                                                                                   // failed to get information about the path
    }
    return (info.st_mode & S_IFDIR) != 0;                                                               // if the path corresponds to a folder and not a file
}
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------
// A function to read any file that is formatted like a parameter file. 
// i.e, in the format:   ParameterName = Value  # comments
//------------------------------------------------------------------------------------------------------------------------
void readParameterFile(const std::string* paramfile, NGSParameters& parameter){
    std::ifstream Pfile(*paramfile);
    std::string line;

    // Loop through each line in the file and extract value for each parameter. Ignore comments starting with '#'
    while(std::getline(Pfile, line)){
        // Process only non-empty lines
        if (!line.empty()){

            // Find the first '#' character in a line and ignore everything after '#' 
            line = line.substr(0, line.find_first_of('#'));

            // Seperate ParameterName and Value from each line
            std::string data_before_equalSign{line.substr(0,line.find("="))};                            // line data upto '=' from the start of the line
            std::string data_after_equalSign{line.substr(line.find("=")+1)};                             // line data after '=' to the end of the line

            // Remove whitespaces before and after both ParameterName and Value
            std::string space{" \t\n\r"};                                                                // Whitespace characters to be removed
            size_t parameterName_begin{data_before_equalSign.find_first_not_of(space)};                  // Get the beginning of parameter name
            size_t parameterName_end{data_before_equalSign.find_last_not_of(space)};                     // Get the end of parameter name
            if (parameterName_begin != std::string::npos && parameterName_end !=std::string::npos){      // Process only lines with non-whitespace characters
                std::string parameterName{data_before_equalSign.substr(parameterName_begin, parameterName_end-parameterName_begin+1)};
                
                size_t paramValue_begin{data_after_equalSign.find_first_not_of(space)};                  // Get the beginning of parameter value
                size_t paramValue_end{data_after_equalSign.find_last_not_of(space)};                     // Get the end of parameter value
                if (paramValue_begin != std::string::npos && paramValue_end !=std::string::npos){        // Process if the value field is not empty
                    std::string parameterValue{data_after_equalSign.substr(paramValue_begin, paramValue_end-paramValue_begin+1)};
                    parameter.set_parameters(&parameterName, &parameterValue);
                }else{                                                                                   // If parameter specified but value is empty, print a error and exit
                    std::cerr<<"\n ERROR: Parameter \""<<parameterName<<"\" cannot be specified with an empty value\n";
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
    Pfile.close();
}
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------
// A function to count the number of exposure entries included in an SDD file and also to build a list of lines locations
// that corresponds to the starting of new exposures. This list will have (number of exposure +1) elements. The last 
// element would be the end position of the SDD file
//------------------------------------------------------------------------------------------------------------------------
int countExposuresSDD(const std::string* sddfile, std::vector<std::streampos>& exposureLines){
    if(!checkFileExists(sddfile)){
        std::cerr<<"\n ERROR: The SDD file provided: "<<*sddfile<<" is invalid\n";
        exit(EXIT_FAILURE);
    }
    std::ifstream Sfile(*sddfile);                                                                        // Starting an ifstream instance with SDD filename
    std::string line;                                                                                     // Temporary variable to hold each line read
    bool end_of_header_found(false);                                                                      // A flag to check if the end of the header reached
    int exp_count(0);                                                                                     // Counter for every new exposure found. Default is 0.
    std::streampos current_line;                                                                          // Temporary variable to store the starting position of the line that is currntly read

    while(std::getline(Sfile, line)){                                                                     // Loop through each line in the file
        if(line.empty()){
            current_line = Sfile.tellg();                                                                 // Get the end position of the current line that is being read. This would be the starting position of the next
            continue;
        }else{                                                                                            // Process only non-empty lines
            if(!end_of_header_found){
                if(line[0]=='*'){                                                                         // End of the header line starts with '*' character
                    end_of_header_found = true;
                }
                current_line = Sfile.tellg();                                                             // Get the end position of the current line that is being read. This would be the starting position of the next
                continue;
            }
            if(line[0] == '2'){                                                                           // New exposure entries starts with 2 in data field 1
                exp_count++;                                                                              // Increment the exposure counter
                exposureLines.push_back(current_line);                                                    // Add current line to the list of all new exposure lines
            }
            current_line = Sfile.tellg();                                                                 // Get the position following the current line that was just read. This would be the starting position of next line
        }
    }
    exposureLines.push_back(current_line);                                                                // Store the final position in the file after all lines are read
    return exp_count;
}
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------
// A function to read any file that is formatted like an actual dose file with 1D data. This function will only read
// number of lines specified by the value lineCount and append these lines into a vector and returns the vector not a copy. 
// If the data entries in the file are fewer than the lineCount, print ERROR and exit 
//------------------------------------------------------------------------------------------------------------------------
std::vector<double> readActualDosefile(const std::string* dosefile, int lineCount){
    if(!checkFileExists(dosefile)){
        std::cerr<<"\n ERROR: The dose file provided: "<<*dosefile<<" is invalid\n";
        exit(EXIT_FAILURE);
    }
    std::ifstream Dfile(*dosefile);                                                                       // Starting an ifstream instance with SDD filename
    std::string line;                                                                                     // Temporary variable to hold each line read
    std::vector<double> vec;                                                                              // Temporary vector to store the read lines (dose values)
    int missingEntries(lineCount);                                                                        // Temporary variable to hold the number of missing entries in file

    // Loop through each line in the file
    while(std::getline(Dfile, line)){
        
        if(!Dfile.eof() && !line.empty() && lineCount>=1){                                                // If the non-empty line is not the last line and we need dose values still
            vec.push_back(std::stod(line));                                                               // Append the dose value to the temporary vector
            lineCount--;                                                                                  // Decrement the number of lines to read
        }else if(Dfile.eof() && !line.empty() && lineCount==1){                                           // If it is the last non-empty line, but we only need one more data value (same line)
            vec.push_back(std::stod(line));                                                               // Append the dose value to the temporary vector
            lineCount--;                                                                                  // Decrement the number of lines to read
        }

        if(lineCount==0){break;}                                                                          // If enough dose values are already read, break the loop
        if(!line.empty()){missingEntries--;}                                                              // count the number of non-empty lines and subtract it from total required entries to get missing data count 
    }
    if(lineCount!=0){                                                                                     // If enough dose values are not present in the file, print error and exit   
        std::cerr<<"\n ERROR: The "<<*dosefile<<" is missing "<<missingEntries<<" more data entries than required \n";
        exit(EXIT_FAILURE); 
    }
    return (vec);                                                                                         // return the vector
}
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------
// A function to read any file that is formatted like an SDD file and get specific header fields 
//------------------------------------------------------------------------------------------------------------------------
void readSDDfileHeader(const std::string* sddfile, NGSsdd& SDDdata){
    if(!checkFileExists(sddfile)){
        std::cerr<<"\n ERROR: The SDD file provided: "<<*sddfile<<" is invalid\n";
        exit(EXIT_FAILURE);
    }
    std::ifstream Sfile(*sddfile);                                                                        // Starting an ifstream instance with SDD filename
    std::string line;                                                                                     // Temporary variable to hold each line read
    
    while(std::getline(Sfile, line)){                                                                     // Loop through each line in the file
        if (!line.empty()){                                                                               // Process only non-empty lines
            if(line[0]!='*'){                                                                             // End of the header line starts with '*' character
                if(line.find("Dose or fluence")==0){                                                      // Field 11 of the SDD header
                    SDDdata.set_expected_dose_gy(&line); 
                }else if(line.find("Chromosome sizes")==0){                                               // Field 15 of the SDD header
                    SDDdata.set_chrom_size_bp(&line);
                    SDDdata.set_cell_ploidy_and_chrom_mappping(&line);
                }
            }else{                                                                                        // when end of the header is reached
                break;
            }                                                                                 
        }
    }
}
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------
// This function will read one exposure data from one SDD file at a time completely independantly depending on the values
// of SDDfilenumber and exposureNum parameters passed. Function then sets all the essential SDD data fields 
//------------------------------------------------------------------------------------------------------------------------
void readSDDfileData(const std::string* sddfile, NGSsdd& SDDdata, int SDDfilenumber, int exposureNum, std::vector<std::string>& lineStack, int groupTID){
    #pragma omp single
    {
        std::ifstream Sfile(*sddfile);                                                                    // Starting an ifstream instance with SDD filename
        std::string line;                                                                                 // Temporary variable to hold each line read
        std::streampos start_pos = SDDdata.get_exposureLines(SDDfilenumber, exposureNum);                 // Get the line number corresponding to the start of exposure in this SDD file
        std::streampos end_pos = SDDdata.get_exposureLines(SDDfilenumber, exposureNum+1);                 // Get the line number corresponding to the end of this exposure data in this SDD file
        lineStack.clear();                                                                                // Empty the vector for each SDD file
        Sfile.seekg(start_pos);                                                                           // Set the position of the stream to the starting of new exposure data 
        std::streampos currentPos = Sfile.tellg();                                                             
        while(currentPos < end_pos){                                                                      // Process only the lines that correspond to this exposure data
            std::getline(Sfile, line);  
            if(!line.empty()){
                lineStack.push_back(line);                                                                // Populate the lineStack vector
            }
            currentPos = Sfile.tellg();
        }
    }
    #pragma omp barrier

    if(!lineStack.empty()){                                                                               // If the line stack is not empty
        #pragma omp single
        {
            int nWorkThreads = omp_get_num_threads();
            SDDdata.reset_workThread_data_holders(groupTID);                                              // Reset all sub-group thread data holding parameters
            SDDdata.init_set_workThread_data_holders(groupTID, nWorkThreads);                             // Resize all sub-group thread data holding vectors
        }
        #pragma omp barrier
        #pragma omp for
        for(size_t i=0; i<lineStack.size(); i++){                                                         // Iterate over the line stack and process each line independantly
            int workerTID = omp_get_thread_num();                                                         // Get the worker thread ID
            SDDdata.set_all_sdd_data_fields(&lineStack[i], groupTID, workerTID);
        }
    }

    return;
}
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------
// A function to create a file with an associated memory map. This function takes in the path to the file you want to create,
// and the file size you are expecting to reserve. A pointer to the memory map of the generated file will be returned.
//------------------------------------------------------------------------------------------------------------------------
char* createMemoryMappedFile(const std::string& filePath, size_t file_size) {
    // This is a system call to open the file. O_RDWR: File should be opened with read and write rights; O_CREAT: Create the file if not already existing. S_IRUSR | S_IWUSR: Give the user (ownder) the read and right previlages respectively
    int fileData = open(filePath.c_str(), O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);                                                                                                        
    if (fileData == -1) {                                                                                 // If the file cannot be opened
        perror(("\nERROR: Failed to create the file "+filePath+" \n").c_str());
        exit(EXIT_FAILURE);                                                                               // Exit with an error
    }

    if (ftruncate(fileData, file_size) == -1) {                                                           // Get the file in appropriate size
        perror(("\nERROR: Failed to successfully create the file "+filePath+". \n").c_str());;
        close(fileData);
        exit(EXIT_FAILURE);                                                                               // If failed, exit with error
    }

    // mmap function is used to map the file. PROT_READ | PROT_WRITE flags are for read and write permissions, MAP_SHARED sets the memory map accessible by other proceses and 0 is an offset to indicate mapping starts at the beginning of the file
    char* mapped_data = static_cast<char*>(mmap(nullptr, file_size, PROT_READ | PROT_WRITE, MAP_SHARED, fileData, 0));
    if (mapped_data == MAP_FAILED) {                                                                      // If the mapping was unsuccessful
        perror(("\nERROR: Failed to generate the memory map of the file "+filePath+". \n").c_str());;
        close(fileData);
        exit(EXIT_FAILURE);                                                                               // Exit with error
    }

    close(fileData);                                                                                      // Close the opened file

    return mapped_data;                                                                                   // Return the pointer to the mapped data for that file
}
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------
// A function to generate the memory-map of an input file that we want to read multiple times. This function takes in the 
// path to the input file and a dummy size_t variable to hold the file size. A pointer to the memory map is returned. 
//------------------------------------------------------------------------------------------------------------------------
char* generateInputFileMemoryMap(const std::string& inFilePath, size_t& inFileSize) {
    int inFile_fd = open(inFilePath.c_str(), O_RDONLY);                                                   // This is the system call to open the file for read-only
    if (inFile_fd == -1) {                                                                                // Checks if the file open call was unsuccessful. Exit with error if it was.
        perror(("\nERROR: Error opening input file: "+inFilePath+"\n").c_str());
        exit(EXIT_FAILURE);
    }

    struct stat inFile_sb;                                                                                // Struct object to hold the opened file information
    if (fstat(inFile_fd, &inFile_sb) == -1) {                                                             // Check if the system call 'fstat' can successfully access the file information
        perror(("\nERROR: Error getting input file information: "+inFilePath+"\n").c_str());
        exit(EXIT_FAILURE);                                                                               // Exit with error if the file information can't be retrieved
    }

    inFileSize = inFile_sb.st_size;                                                                       // Getting the size of the input file
    // mmap function is used to generate the memory map of the input file. PROT_READ for read-only permissions, MAP_PRIVATE for keep original file unaffected by the modifications to the mapped data, and 0 is the offset
    char* inFile_data = static_cast<char*>(mmap(0, inFileSize, PROT_READ, MAP_PRIVATE, inFile_fd, 0));
    if (inFile_data == MAP_FAILED) {                                                                      // If memory mapping fails
        perror(("\nERROR: Failed mapping input file to memory: "+inFilePath+"\n").c_str());
        close(inFile_fd);                                                                                 // Close the opened file
        exit(EXIT_FAILURE);                                                                               // Exit with error
    }
    close(inFile_fd);                                                                                     // Close the opened file

    return inFile_data;                                                                                   // Return the pointer to the memory map of the inputFile
}
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------
// A function to re-map the memory-map of an input file when the data exceeds allotted filesize. This function takes in the 
// initially mapped data, respective file descriptor, current file size and the new required filesize. It will update the 
// mapped_data and currentFileSize variables
//------------------------------------------------------------------------------------------------------------------------
void resizeMemoryMappedFile(char*& mapped_data, int& fileDescriptor, size_t& currentFileSize, size_t newFileSize) {
    if (ftruncate(fileDescriptor, newFileSize) == -1) {                                                   // Change the size of the physical file in storage to the new size
        perror("\nERROR: Failed to resize the memory-mapped file.\n");                                    // If this cannot be done, exit with error
        exit(EXIT_FAILURE);
    }

    char* new_mapped_data = static_cast<char*>(mmap(nullptr, newFileSize, PROT_READ | PROT_WRITE, MAP_SHARED, fileDescriptor, 0));// Map the re-sized file again to a new pointer
    if (new_mapped_data == MAP_FAILED) {                                                                  // Exit with an error if this mapping of the new file failed
        perror("\nERROR: Failed to remap the memory-mapped file with a new size.\n");
        exit(EXIT_FAILURE);
    }

    munmap(mapped_data, currentFileSize);                                                                 // Unmap the old mapping from RAM

    mapped_data = new_mapped_data;                                                                        // Replace the old memory-map pointer with the new one
    currentFileSize = newFileSize;                                                                        // Update the file size to correspond to the change
}
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------
// A function to read each chromosome sequence one-by-one from the reference sequence. It will get the chromosome IDs, 
// read the chromosome sequence, and merge multiple lines together. The string passed by reference will get the value of 
// the chromosome sequence. By calling this function in a loop, we can get all the chromosome seqs one-by-one
//------------------------------------------------------------------------------------------------------------------------
int getNextChromSeq(std::ifstream& ref_seq, std::string& chrom_seq, std::string& chromSeq_ID){
    chrom_seq.clear();                                                                                    // Clear the previously stored value if any
    chromSeq_ID.clear();                                                                                  // Clear the previously stored value of sequence ID if any
    std::string tmpString;                                                                                // Temporary vector to hold each line value
    int success_flag{0};                                                                                  // This will return 0 if no more sequences left to read
    char tmpChar = '-';
    std::ostringstream tmpStringBuffer;                                                                   // Temporary Buffer string to store sequence strings to be combined

    while(tmpChar!='>' && !ref_seq.eof() ){                                                               // Skip any empty fields before the chromosome ID
        ref_seq.get(tmpChar);
    }
    if(tmpChar!='>'){                                                                                     // Return if '>' char is nowhere to be found in the file
        return 0;
    }
    getline(ref_seq,chromSeq_ID);                                                                         // Read the chromosome ID; Important step to advance to the next line
    
    while((!ref_seq.eof()) && (ref_seq.peek()!='>')){                                                     // Read the sequences over multiple lines until the next ID or end
        tmpString.clear();
        getline(ref_seq,tmpString);
        if(tmpString[0]!='>' && !tmpString.empty()){                                                      // If it is really a sequence line and not empty, then
            if (tmpString[tmpString.size() - 1] == '\r'){ tmpString.resize(tmpString.size() - 1); }       // Remove the carriage return function from the end of each line
            tmpStringBuffer << tmpString;                                                                 // Merge multiple lines of chrom seq to one buffer string
            success_flag++;
        }
    }
    chrom_seq = tmpStringBuffer.str();                                                                    // Passing the combined string to the chrom_seq variable
    return(success_flag);                                                                                 // Return 0 if either reached end or if no sequence found
}
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------
// A function to read each chromosome sequence one-by-one from the reference sequence Memory Map. It will get the chromosome IDs, 
// read the chromosome sequence, and merge multiple lines together. The string passed by reference will get the value of 
// the chromosome sequence. By calling this function in a loop, we can get all the chromosome seqs one-by-one
//------------------------------------------------------------------------------------------------------------------------
int getNextChromSeq_MM(const char* refSeqData, size_t refSize, size_t& position, std::string& chrom_seq, std::string& chromSeq_ID) {
    chrom_seq.clear();                                                                                    // Clear the previously stored value if any
    chromSeq_ID.clear();                                                                                  // Clear the previously stored value of sequence ID if any
    std::ostringstream tmpStringBuffer;                                                                   // Temporary Buffer string to store sequence strings to be combined
    int success_flag = 0;                                                                                 // This will return 0 if no more sequences are left to read

    if (position >= refSize) {                                                                            // If we've reached the end of the memory-mapped data, return 0.
        return 0;  
    }

    while (refSeqData[position] != '>' && refSeqData[position] != '\0') {
        position ++;
    }

    if (refSeqData[position] == '>') {                                                                    // If the character is '>', it is the start of the chrom ID
        size_t end = position;                                                                            // A variable to hold the end position of the chrom ID
        while (end < refSize && refSeqData[end] != '\n') {                                                // Find the end of the chrom ID line
            end++;
        }
        chromSeq_ID = std::string(&refSeqData[position], end - position);                                 // Get a substring from 'position' to 'end' from genomeTempate_data. Extracting chrom ID
        chromSeq_ID.erase(chromSeq_ID.find_last_not_of(" \n\r\t") + 1);                                   // Trim whitespaces at the end
        position = end + 1;                                                                               // Move the position to the next line
        
        // Get the strand sequence
        while (position < refSize && refSeqData[position] != '>') {
            if (!(std::isspace(refSeqData[position]) || refSeqData[position] == '\n' || refSeqData[position] == '\r' || refSeqData[position] == '\t')) {
                chrom_seq += refSeqData[position];
                position++;                                                                               // Move to the next character
            } else {
                 position++;                                                                              // Skip whitespace and escape characters
            }
        }
        success_flag++;
    }

    if (!success_flag) {                                                                                  // Return false if the forward strand cannot be read
        return 0;
    }

    return success_flag;  // Return 0 if either reached the end or if no sequence found
}
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------
// A function to write to file in batch. It takes a batch_buffer, which is a buffer vector of strings that needs to be
// written to the outputFile. This function writes the buffer to the file and clears the buffer for new data
//------------------------------------------------------------------------------------------------------------------------
void writeBatchToFile(std::vector<std::string>& batch_buffer, std::ofstream& outputFile){
    for (const std::string& line : batch_buffer) {                                                        // Parsing through the batch_buffer 
        if (!line.empty()) {                                                                              // Skip empty lines if any
            outputFile << line;                                                                           // Writing each line to the outputfile
        }
    }
    batch_buffer.clear();                                                                                 // Clear the buffer after writing it to file
}
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------
// A function to write to a memory-mapped file in batch. It takes a batch_buffer, which is a buffer vector of strings that
// needs to be written to the outputFile. This function writes the buffer to the MMfile and clears the buffer for new data.
// If writing will exceed the allocated file size, and if the path to the physical file is also provided, then this funciton
// will call the resizeMemoryMappedFile function which will re-size and re-map the file, and return the new mapped file. 
//------------------------------------------------------------------------------------------------------------------------
void writeBatchToMMFile(std::vector<std::string>& batch_buffer, char*& position_in_File, char*& fileMapping, size_t& fileSize, std::string filePath){
    for (const std::string& data : batch_buffer){                                                         // Parsing through the batch_buffer 
        if (!data.empty()) {                                                                              // Skip empty lines if any
            long spaceLeft_in_MM = (reinterpret_cast<char*>(fileMapping+fileSize)-(position_in_File+data.size())); // Remaining space in MM = end pointer - current pointer. -ve if the allotted storage limit exceeded. It's a pointer operation.
            if (spaceLeft_in_MM < 0){                                                                     // If there is not enough space in the memory-mapped file to write this data
                if (filePath != "NA"){                                                                    // If a valid file path is given to indicate memory-map resizing should be done, then
                    int fileDescriptor = open(filePath.c_str(), O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);     // Reopen the file that we want to re-size to get the file descriptor
                    if (fileDescriptor == -1) {
                        perror("ERROR: Failed to reopen the file");
                        exit(EXIT_FAILURE);
                    }
                    size_t newFileSize = (fileSize-spaceLeft_in_MM) + 10000000;                           // What is immediately needed (initial fize size + more required for the current data) + 10 MB extra (in bytes)
                    ptrdiff_t  offset = position_in_File - fileMapping;                                   // Offset determines how far was the current position in the initially mapped file. ptddiff_t is a pointer difference type
                    resizeMemoryMappedFile(fileMapping,fileDescriptor,fileSize,newFileSize);              // Re-map the file after changing it's size
                    position_in_File = fileMapping + offset;                                              // Modify the position in file where we want to write next to the new map + offset (end of the data in map)
                }else{
                    std::cerr << "\nError: Insufficient space in memory-mapped file\n";                   // If there is not enough space to write the next set of data and memory map is not asked to be resized, print error and break
                    exit(EXIT_FAILURE);
                }
            }
            memcpy(position_in_File, data.c_str(), data.size());                                          // Copy the data from the vector to the memory-mapped file
            position_in_File += data.size();                                                              // Increment the memory map file pointer so that the next line to write can start from there
        }
    }
    batch_buffer.clear();                                                                                 // Clear the buffer after writing it to file
}
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------
// A function to read the undamaged genome template. Every odd and even sequences will be read simultaneously. Odd
// sequence corresponds to the forward strand sequence and the even ones corresponds to the reverse complementary strand sequence. 
// They will be stored to A and B tagged strings respectively (Not optimized version)
//------------------------------------------------------------------------------------------------------------------------
/* int readFastaTemplate(std::ifstream& fastaFile, std::string& chrom_ID_A, std::string& chrom_seq_A, std::string& chrom_ID_B, std::string& chrom_seq_B) {
    chrom_seq_A.clear();                                                                                  // Clear the previously stored value if any
    chrom_seq_B.clear();
    std::string tmpString;                                                                                // Temporary vector to hold each line value
    int success_flag{0};                                                                                  // This will return 0 if no more sequences left to read
    char tmpChar = '-';
   
    while(tmpChar!='>' && !fastaFile.eof() ){                                                             // Skip any empty fields before the chromosome ID
        fastaFile.get(tmpChar);
    }
    
    if(tmpChar!='>'){                                                                                     // Return if '>' char is nowhere to be found in the file
        return 0;
    }
    
    // Get the forward strand ID and sequence
    getline(fastaFile,chrom_ID_A);                                                                        // Read the chromosome ID and store it in the string passed
    chrom_ID_A.insert(0,1,'>');                                                                           // Reintroduce '>' in the beginning of the chrom A ID, since the get() function removed it
    while((!fastaFile.eof()) && (fastaFile.peek()!='>')){                                                 // Read the sequences over multiple lines until the next ID or end
        tmpString.clear();
        getline(fastaFile,tmpString);
        if(tmpString[0]!='>' && !tmpString.empty()){                                                      // If it is really a sequence line and not empty, then
            if (tmpString[tmpString.size() - 1] == '\r'){ tmpString.resize(tmpString.size() - 1); }       // Remove the carriage return function from the end of each line
            chrom_seq_A.append(tmpString);                                                                // Merge multiple lines of chrom seq to one
            success_flag++;
        }
    }
    
    if (!success_flag){return 0;}                                                                         // Return false if the forward strand cannot be read
    
    // Get the reverse complementary strand ID and sequence
    getline(fastaFile,chrom_ID_B);                                                                        // Read the chromosome ID and store it in the string passed
    while((!fastaFile.eof()) && (fastaFile.peek()!='>')){                                                 // Read the sequences over multiple lines until the next ID or end
        tmpString.clear();
        getline(fastaFile,tmpString);
        if(tmpString[0]!='>' && !tmpString.empty()){                                                      // If it is really a sequence line and not empty, then
            if (tmpString[tmpString.size() - 1] == '\r'){ tmpString.resize(tmpString.size() - 1); }       // Remove the carriage return function from the end of each line
            chrom_seq_B.append(tmpString);                                                                // Merge multiple lines of chrom seq to one
            success_flag++;
        }
    }
    return(success_flag);                                                                                 // Return 0 if either reached end or if no sequence found
} */
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------
// A function to read the undamaged genome template. Every odd and even sequences will be read simultaneously. Odd
// sequence corresponds to the forward strand sequence and the even ones corresponds to the reverse complementary strand sequence. 
// They will be stored to A and B tagged strings respectively (This is an optimized version)
//------------------------------------------------------------------------------------------------------------------------
/* int readFastaTemplate(std::ifstream& fastaFile, std::string& chrom_ID_A, std::string& chrom_seq_A, std::string& chrom_ID_B, std::string& chrom_seq_B) {
    chrom_seq_A.clear();                                                                                  // Clear the previously stored value if any
    chrom_seq_B.clear();
    std::string tmpString;                                                                                // Temporary vector to hold each line value
    int success_flag{0};                                                                                  // This will return 0 if no more sequences left to read
    char tmpChar;

    // Get the forward strand ID
    while (fastaFile >> tmpChar) {                                                                        // Assign each character from file to tmpChr
        if (tmpChar == '>') {                                                                             // If the '>' is found in the file
            getline(fastaFile, chrom_ID_A);                                                               // Read the chromosome ID and store it in the string passed
            chrom_ID_A = '>' + chrom_ID_A;                                                                // Reintroduce '>' in the beginning of the chrom A ID, since the get() function removed it
            break;                                                                                        // Break the loop if chrom_ID_A is found
        }
    }

    if (chrom_ID_A.empty()) {                                                                             // If no '>' character is found in the entire file, then exit
        return 0;
    }

    // Get the forward strand sequence
    while (fastaFile >> tmpString) {                                                                      // Assign the next whole string to tmpString variable
        tmpString.erase(std::remove(tmpString.begin(), tmpString.end(), '\r'), tmpString.end());          // Remove any carriage return ('\r') characters from tmpString
        chrom_seq_A += tmpString;                                                                         // Forward sequence
        success_flag++;
        break;                                                                                            // Break the loop before proceeding to the next line in the file
    }

    if (!success_flag) {                                                                                  // Return false if the forward strand cannot be read
        return 0;
    }

    // Read the reverse complementary strand ID
    while (fastaFile >> tmpChar) {
        if (tmpChar == '>') {
            getline(fastaFile, chrom_ID_B);
            chrom_ID_B = '>' + chrom_ID_B;                                                                // Reintroduce '>' in the beginning of the chrom A ID, since the get() function removed it
            break;                                                                                        // Break the loop if the second '>' character corresponding to the reverse strand ID is found
        }
    }

    // Read the reverse complementary strand sequence
    while (fastaFile >> tmpString) {                                                                      // Assign the next whole string to tmpString variable
        tmpString.erase(std::remove(tmpString.begin(), tmpString.end(), '\r'), tmpString.end());          // Remove any carriage return ('\r') characters from tmpString
        chrom_seq_B += tmpString;                                                                         // Reverse sequence
        success_flag++;
        break;                                                                                            // Break the loop before proceeding to the next line in the file
    }

    return success_flag;                                                                                  // Return 0 if either reached end or if no sequence found
} */
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------
// A function to read the undamaged genome from a memory map instead of reading from file. Every odd and even sequences will 
// be read simultaneously. Odd sequence corresponds to the forward strand sequence and the even ones corresponds to the 
// reverse complementary strand sequence. They will be stored to A and B tagged strings respectively 
//------------------------------------------------------------------------------------------------------------------------
int readFastaMemoryMap(const char* genomeTemplate_data, size_t templateSize, size_t& position, std::string& chrom_ID_A, std::string& chrom_seq_A, std::string& chrom_ID_B, std::string& chrom_seq_B) {
    chrom_seq_A.clear();                                                                                  // Clear the previously stored value if any
    chrom_seq_B.clear();
    chrom_ID_A.clear();
    chrom_ID_B.clear();
    int success_flag{0};                                                                                  // This will return 0 if no more sequences left to read
    
    if (position >= templateSize) {                                                                       // If we've reached the end of the memory-mapped data, return 0.
        return 0;  
    }

    // Process the forward strand ID and sequence
    if (genomeTemplate_data[position] == '>') {                                                           // If the character is '>', it is the start of the chrom ID
        size_t end = position;                                                                            // A variable to hold the end position of the chrom ID
        while (end < templateSize && genomeTemplate_data[end] != '\n') {                                  // Find the end of the chrom ID line
            end++;
        }
        chrom_ID_A = std::string(&genomeTemplate_data[position], end - position);                         // Get a substring from 'position' to 'end' from genomeTempate_data. Extracting chrom ID
        chrom_ID_A.erase(chrom_ID_A.find_last_not_of(" \n\r\t") + 1);                                     // Trim whitespaces at the end
        position = end + 1;                                                                               // Move the position to the next line
        
        // Get the forward strand sequence
        size_t start = position;
        while (position < templateSize && genomeTemplate_data[position] != '>') {                         // Find the end of the chrom sequence by looking for the next '>'
            position++;
        }
        chrom_seq_A = std::string(&genomeTemplate_data[start], position - start);                         // Get the substring corresponding to the chrom_seq
        chrom_seq_A.erase(chrom_seq_A.find_last_not_of(" \n\r\t") + 1);                                   // Trim whitespaces at the end
        
        success_flag++;
    }

    if (!success_flag) {                                                                                  // Return false if the forward strand cannot be read
        return 0;
    }

    // Process the reverse complementary strand ID and sequence
    if (position < templateSize && genomeTemplate_data[position] == '>') {                                // If the character is '>', it is the start of the chrom ID
        // Get the reverse complementary strand ID
        size_t end = position;
        while (end < templateSize && genomeTemplate_data[end] != '\n') {                                  // Find the end of the chrom ID line
            end++;
        }
        chrom_ID_B = std::string(&genomeTemplate_data[position], end - position);                         // Get the substring corresponding to the chrom ID
        chrom_ID_B.erase(chrom_ID_B.find_last_not_of(" \n\r\t") + 1);                                     // Trim whitespace at the end
        position = end + 1;                                                                               // Move the position to the next line
        
        // Get the reverse complementary strand sequence
        size_t start = position;
        while (position < templateSize && genomeTemplate_data[position] != '>') {                         // Find the end of the chrom sequence by looking for the next '>'
            position++;
        }
        chrom_seq_B = std::string(&genomeTemplate_data[start], position - start);                         // Get the substring corresponding to the chrom_seq  
        chrom_seq_B.erase(chrom_seq_B.find_last_not_of(" \n\r\t") + 1);                                   // Trim whitespaces at the end
        
        success_flag++;
    }
    return success_flag;
}
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------
// This function makes read quality distribution vector based on the read quality profile provided in the 'read_quality_file'.
// These read quality files are expected to be formatted according to the ART read simulation program profiler. 
// The data must be arrayed in pairs such that each line in the pair has the same identifier and position number.  
// The first line in a pair is a list of quality scores in ascending order and the second line are the corresponding cumulative
// frequencies of the quality scores. Identifier can be 'A','T','G','C' and '.' . The identifier indicate if the quality score
// is specific to a nitrogen base or if it is specified for the combination of all base with '.' 
// This function will only consider quality scores with identifier '.' in this program. Need modifications if base-specific
// score to be considered.
//------------------------------------------------------------------------------------------------------------------------
void make_quality_distribution(std::ifstream& read_quality_file, std::vector<std::map<unsigned int, unsigned short>>& quality_distribution_vec){
    char quality_identifier;                                                                              // Variable to hold the quality identifier character
    int read_pos;                                                                                         // Variable to store the read position for which the score is accociated
    
    while(!read_quality_file.eof()){
        std::string line;
        
        // Process the first line in the line pair corresponding to the quality score 
        getline(read_quality_file, line);                                                                 // Get the first line from the read_quality_profile file provided
        if(line.length() == 0) continue;                                                                  // Continue if the line is empty
		if(line[0] != '.') continue;                                                                      // Process the quality score lines starting with '.' only. These lines are the quality scores that are not specific to a nitrogen base
        std::istringstream lineString(line);                                                              // Convert the line into a string stream
        lineString>>quality_identifier;                                                                   // Assign quality identifier from the lineString. Important step to remove the character from the string
        lineString>>read_pos;                                                                             // Assign read position from the lineString. This is needed to remove the value from the string

        unsigned short score;                                                                             // Individual quality score values from the line string
        std::vector<unsigned short> quality_score;                                                        // Vector to hold all the quality scores
        while (lineString >> score){                                                                      // Get each score value from the string stream and assign it to score variable
            quality_score.push_back(score);                                                               // Populate the quality score vector
        }

        // Process the second line in the line pair corresponding to the cumulative frequency of the score
        getline(read_quality_file, line);                                                                 // Get the second line from the line pair
        lineString.clear();                                                                             
        lineString.str(line);                                                                             // Assign the new line as lineString
        lineString>>quality_identifier;
        lineString>>read_pos;

        unsigned long frequency;                                                                          // Individual cumulative frequency value
        std::vector<unsigned long> cumulative_frequency;                                                  // Vector to hold all the cumulative frequencies
        while (lineString >> frequency){  
            cumulative_frequency.push_back(frequency); 
        }

        double denominator = cumulative_frequency[cumulative_frequency.size()-1] / 1000000.0;             // denominator = Total cumulative frequency / 1000000
        std::map<unsigned int, unsigned short> distribution;                                                

        for(size_t i=0; i<cumulative_frequency.size(); i++){
            unsigned int cc = static_cast<unsigned int>(ceil(cumulative_frequency[i]/denominator));           
            distribution[cc] = quality_score[i];
        }
        if(distribution.size()>0){
            quality_distribution_vec.push_back(distribution);
        }
    }
}
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------
// This function will generate a summary report at the end of the run. It will create a file called "Run_summary_report.txt"
// in the output directory
//------------------------------------------------------------------------------------------------------------------------
void generate_run_summaryReport(NGSParameters& parameter, NGSsdd& SDDdata){
    std::string report_filename = *parameter.get_output_directory()+"/Run_summary_report.txt";
    std::ofstream report_file(report_filename.c_str());
    std::time_t currentTime = std::time(nullptr);                                                         // Get the current system time
    std::string timeString = std::ctime(&currentTime);                                                    // Convert the current time to a string representation
    report_file<<"------------------------------------------------------------------------------------------------------------\n"
               <<"                                   RADISEQ SIMULATOR RUN SUMMARY REPORT             \n"
               <<"------------------------------------------------------------------------------------------------------------\n"
               <<" This report was prepared on "<<timeString<<"\n\n"
               <<" Time duration of this run                                            : "<<(report_cpu_time_used/60)<<" min\n";
    if(parameter.get_random_seed() == 0){
        report_file<<" The seed used for this run                                           : Random seed\n";
    }else{
        report_file<<" The random seed used for this run                                    : "<<parameter.get_random_seed()<<"\n";
    }
    if(parameter.get_number_of_threads() == 1){
        report_file<<" Multi-threading option                                               : Disabled\n";
    }else{
        report_file<<" Multi-threading option                                               : Enabled\n"
                   <<" Number of threads used for the run                                   : "<<parameter.get_number_of_threads()<<"\n";
    }
    report_file<<" The parameter file specified                                         : "<<report_parameterFileName<<"\n\n"
               <<" ------ SDD file information -------\n"
               <<" Number of cells irradiated in the Monte Carlo Simulation             : "<<SDDdata.get_num_of_exposures()<<"\n";
                   
    if(parameter.get_merge_damages_from_particles()){
        report_file<<" Number of SDD files given to merge damages from                      : "<<parameter.get_num_of_particles_to_merge()<<"\n"
                   <<" The names of these SDD files are                                     : ";
        for(int i=0;i<parameter.get_num_of_particles_to_merge();i++){
            if(i!=0) report_file<<", ";
            report_file<<parameter.get_sddfile_path()[i];
        }
        report_file<<"\n These SDD files are generated with the primary particles             : ";
        for(int i=0;i<parameter.get_num_of_particles_to_merge();i++){
            if(i!=0) report_file<<", ";
            if(i>static_cast<int>(parameter.get_names_of_particles_to_merge()->size())-1) report_file<<"Unspecified";
            else {
                const std::vector<std::string>* particle_names = parameter.get_names_of_particles_to_merge();
                report_file<<(*particle_names)[i];
            }
        }
        report_file<<"\n The relative dose contributions of these particles                   : ";
        for(int i=0;i<parameter.get_num_of_particles_to_merge();i++){
            if(i!=0) report_file<<", ";
            if(i>static_cast<int>(parameter.get_relative_dose_contributions().size())-1) report_file<<"0";
            else {
                report_file<<parameter.get_relative_dose_contributions()[i];
            }
        }    
    }else{
        report_file<<" Name of the SDD file provided                                        : "<<parameter.get_sddfile_path()[0]<<"\n";
        const std::vector<std::string>* particle_names = parameter.get_names_of_particles_to_merge();
        report_file<<" Name of the primary particle used for the irradiation                : "<<(*particle_names)[0];              
    }
    report_file<<"\n The dose delivered to a single cell during the irradiation           : "<<SDDdata.get_expected_dose_gy()<<" Gy\n";
    if(parameter.get_adjust_damages_with_actual_dose()){
        report_file<<" Are the radiation-induced damages adjusted for actual delivered dose : Yes \n"
                   <<" Number of damages are adjusted for delivered dose using the files    : ";
        for(int i=0;i<parameter.get_num_of_particles_to_merge();i++){
            if(i!=0) report_file<<", ";
            const std::vector<std::string>* actual_dose_files = parameter.get_actual_dosefile_path();
            report_file<<(*actual_dose_files)[i];
        }
        report_file<<"\n";     
    }else report_file<<" Are the radiation-induced damages adjusted for actual delivered dose : No \n";

    report_file<<" The genome length of the Monte Carlo cell model                      : "<<SDDdata.get_sdd_genome_length()<<" bp\n"
               <<" Name of the reference genome file used                               : "<<*parameter.get_reference_genome()<<"\n"
               <<" The length of the reference genome                                   : "<<report_ref_seq_length<<" bp\n\n"
               <<" ------ Sequencing information -------\n"
               <<" Name of the Illumina sequencer used for NGS simulation               : "<<*parameter.get_sequencer()<<"\n";
    if(*parameter.get_sequencing_mode() == "single"){
        report_file<<" The sequencing protocol used                                         : Single-cell Whole Genome Sequencing\n";
    }else report_file<<" The sequencing protocol used                                         : Bulk-cell Whole Genome Sequencing\n";
    
    if(parameter.get_paired_end_sequencing()){
        report_file<<" The read generation modality used                                    : Paired-end read generation\n"
                   <<" The mean length of DNA fragments created in the sample               : "<<parameter.get_mean_DNA_fragment_length()
                   <<" bp with std. deviation of "<<parameter.get_std_dev_DNA_fragment_length()<<" bp\n";
    }else report_file<<" The read generation modality used                                    : Single-end read generation\n";
    
    int total_reads = std::round(parameter.get_total_read_coverage()*report_ref_seq_length)/parameter.get_read_length();
    report_file<<" Total number of cells in the sample                                  : "<<parameter.get_num_of_cells_in_sample()<<"\n"
               <<" Number of cells randomly sampled from the sample pool for sequencing : "<<parameter.get_num_of_cells_to_sequence()<<"\n"
               <<" Length of the reads simulated                                        : "<<parameter.get_read_length()<<" bp\n"
               <<" Total number reads generated in the NGS simulation (expected)        : "<<total_reads<<"\n";
    if(*parameter.get_sequencing_mode() == "bulk"){
        report_file<<" The number of reads actually generated                               : "<<std::accumulate(report_readsGenerated_perCell.begin(),report_readsGenerated_perCell.end(),0)<<"\n";
    }
    report_file<<" Average number of reads generated per cell                           : "<<total_reads/parameter.get_num_of_cells_to_sequence()<<"\n"
               <<" Total read coverage obtained                                         : "<<parameter.get_total_read_coverage()<<"\n"
               <<" The mean read coverage per cell                                      : "<<parameter.get_total_read_coverage()/parameter.get_num_of_cells_to_sequence()<<"\n";
    
    report_file<<" Maximum number of errors allowed per read                            : ";
    if(parameter.get_max_errors_in_read()<0) report_file<<"Unlimited\n";
    else report_file<<parameter.get_max_errors_in_read()<<"\n";

    if(parameter.get_paired_end_sequencing()){
        report_file<<" The error-rate used for adding insertions in read 1                  : "<<parameter.get_insertion_error_rate_read1()<<"\n"
                   <<" The error-rate used for adding insertions in read 2                  : "<<parameter.get_insertion_error_rate_read2()<<"\n"
                   <<" The error-rate used for deletions in read 1                          : "<<parameter.get_deletion_error_rate_read1()<<"\n"
                   <<" The error-rate used for deletions in read 2                          : "<<parameter.get_deletion_error_rate_read2()<<"\n";
    }else{
        report_file<<" The error-rate used for adding insertions in the read                : "<<parameter.get_insertion_error_rate_read1()<<"\n"
                   <<" The error-rate used for deletions in the read                        : "<<parameter.get_deletion_error_rate_read1()<<"\n";
    } 
    
    if(*parameter.get_sequencing_mode() == "single"){
        report_file<<" The directory where the sequenced FASTQ files are stored             : "<<*parameter.get_output_directory()<<"\n\n"
                   <<" ------ Output data information -------\n\n"
                   <<"------------------------------------------------------------------------------------------------------------\n"
                   <<"|    Name of the cell sequenced                   |    Output FASTQ filename                               |\n"
                   <<"------------------------------------------------------------------------------------------------------------\n";
        for(size_t i=0; i<report_cells_sequenced.size(); i++){
            size_t cell_name_size = report_cells_sequenced[i].size();
            report_file<<std::left<<std::setw(50)<<"|    "+report_cells_sequenced[i].erase(cell_name_size-3)<<"|    "
                       <<std::left<<std::setw(52)<<report_fastq_output[i]+"_R1.fastq.gz"<<"|\n";
            if(parameter.get_paired_end_sequencing()){
                report_file<<std::left<<std::setw(50)<<"|    "<<"|    "<<std::left<<std::setw(52)<<report_fastq_output[i]+"_R2.fastq.gz"<<"|\n";
            }
            report_file<<"|                                                                                                          |\n";
        }
        report_file<<"------------------------------------------------------------------------------------------------------------\n";
    }else{ 
        report_file<<" The directory where the sequenced FASTQ files are stored             : "<<*parameter.get_output_directory()<<"\n\n"
                   <<" ------ Output data information -------\n\n"
                   <<" Name of the output FASTQ file                                        : "<<report_fastq_output[0]<<"_R1.fastq.gz";
        if(parameter.get_paired_end_sequencing()){
            report_file<<", "<<report_fastq_output[0]<<"_R2.fastq.gz";
        }
        report_file<<"\n------------------------------------------------------------------------------------------------------------\n"
                   <<"|    Name of the cell sequenced                   |    Percentage read contribution of the cell            |\n"
                   <<"------------------------------------------------------------------------------------------------------------\n";
        std::map<std::string, int> read_contribution;
        for (size_t i = 0; i < report_cells_sequenced.size(); i++) {
            read_contribution[report_cells_sequenced[i]] += report_readsGenerated_perCell[i];
        }
        for (const auto& entry : read_contribution) {
            std::string cell_name = entry.first;
            size_t cell_name_size = cell_name.size();
            double percent_read = (entry.second/static_cast<double>(total_reads))*100;
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << percent_read;
            std::string truncated_percentRead = ss.str();
            report_file<<std::left<<std::setw(50)<<"|    "+cell_name.erase(cell_name_size-3)<<"|    "
                       <<std::left<<std::setw(52)<<truncated_percentRead+"%"<<"|\n";
        }
        report_file<<"------------------------------------------------------------------------------------------------------------\n";
    }
}
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------
// A function to determine the size of a file that is passed. Return value in bytes
//------------------------------------------------------------------------------------------------------
long fileSize_bytes(const std::string& filename){
    std::ifstream file(filename, std::ifstream::ate | std::ifstream::binary);
    if (!file) {                                                                                          // Return 0 if file cannot be opened
        return 0;
    }
    std::streampos fileSize = file.tellg();
    if (fileSize == -1) {                                                                                 // Return 0 if file size cannot be determined
        return 0;
    }
    file.close();
    return fileSize;
}
//------------------------------------------------------------------------------------------------------