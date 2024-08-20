#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <zlib.h>

// This function will read one line at a time from a .gz file and passes it through the line variable. 
bool gzReadLine(gzFile file, std::string &line) {
    char buffer[8192];                                                                                  // Temporary buffer to store the line data
    char *result = gzgets(file, buffer, sizeof(buffer));                                                // Reads the line to the buffer, without null terminator
    if (!result) return false;                                                                          // Return false if the line cannot be read successfully
    line = buffer;
    if (!line.empty() && line.back() == '\n'){                                                          // If the last character is '\n', then remove it
        line.pop_back();
    }
    return true;
}

int main(int argc, char* argv[]) {
    if (argc != 5){                                                                                     // Check for the correct number of arguments
        std::cerr << "Usage: " << argv[0] << " -f <fastq.gz> -o <outputFile.txt>" << std::endl;
        return 1;
    }

    std::string inputFile;                                                                              // Variable to hold the name of the fatsq.gz file
    std::string outputFile;                                                                             // Variable to hold the name of the output file to generate

    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-f" && i+1 < argc){                                                                 // If the flag -f is seen, read the next argument as the fastq.gz file name
            inputFile = argv[++i];
        } else if (arg == "-o" && i+1 < argc){                                                          // If the flag -o is seen, read the next argument as the output file name
            outputFile = argv[++i];
        } else {                                                                                        // If appropriate arguments are not seen, then print the error
            std::cerr << "Usage: " << argv[0] << " -f <fastq.gz> -o <outputFile.txt>" << std::endl;
            return 1;
        }
    }

    // Check if both files are provided
    if (inputFile.empty() || outputFile.empty()) {
        std::cerr << "Error: Both input and output files must be specified." << std::endl;
        std::cerr << "Usage: " << argv[0] << " -f <fastq.gz> -o <outputFile.txt>" << std::endl;
        return 1;
    }

    // Check to open the fastq.gz file
    gzFile file = gzopen(inputFile.c_str(), "r");
    if (!file) {
        std::cerr << "Error opening input file: " << inputFile << std::endl;
        return 1;
    }

    // Check if the output file is already present. Else, make one. 
    std::ofstream outFile(outputFile);
    if (!outFile) {
        std::cerr << "Error opening output file: " << outputFile << std::endl;
        return 1;
    }

    std::string line;                                                                                   // Temporary string to hold each line from the fastq file
    std::vector<std::vector<int>> quality_distribution;                                                 // Vector to hold the quality distribution vectors

    while (true){
        // Read the identifier line
        if (!gzReadLine(file, line)) break;                                                             // Break the loop if the function cannot read the file or if the end is reached
  
        // Read the sequence line
        if (!gzReadLine(file, line)) break;
        std::string sequence = line;

        // Read the plus line
        if (!gzReadLine(file, line)) break;
  
        // Read the quality line
        if (!gzReadLine(file, line)) break;
        std::string quality = line;

        // Ensure the quality distribution vector is large enough
        if (quality_distribution.size() < quality.size()){                                              // There should be an inside vector per the number of bases in the read
            quality_distribution.resize(quality.size(), std::vector<int>(256, 0));                      // If the current number of inside vectors is smaller, add more inside vectors
        }

        // Accumulate the quality scores
        for (size_t i = 0; i < quality.size(); i++){                                                    // Iterate over each base position in a read
            int q = (static_cast<int>(quality[i]))-33;                                                  // -33 to get quality score from Phred score
            quality_distribution[i][q]++;                                                               // Increment the count of a quality score in the corresponding inside vector
        }
    }

    gzclose(file);

    // Output the quality distribution to the output file
    for (size_t i = 0; i < quality_distribution.size(); i++){
        outFile << ".   "<<i;                                                                           // Start the first line with a period followed by positions
        for (size_t q = 0; q < quality_distribution[i].size(); q++){
            if (quality_distribution[i][q] > 0){
                outFile << "    " << q;                                                                 // Output the positions of non-zero values
            }
        }
        outFile << std::endl;                                                                           // End the first line

        outFile << ".   "<<i;                                                                           // Start the second line with a period followed by values
        for (size_t q = 0; q < quality_distribution[i].size(); q++){
            if (quality_distribution[i][q] > 0){
                outFile << "    " << quality_distribution[i][q];                                        // Output the non-zero values
            }
        }
        outFile << std::endl;                                                                           // End the second line
    }

    outFile.close();

    return 0;
}
