#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <sstream>
#include <cmath>
#include <ctime>
#include <iomanip>

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
// A function to count the number of exposure entries included in an SDD file
//------------------------------------------------------------------------------------------------------------------------
int countExposuresSDD(const std::string* sddfile){
    std::ifstream Sfile(*sddfile);                                                                        // Starting an ifstream instance with SDD filename
    std::string line;                                                                                     // Temporary variable to hold each line read
    bool end_of_header_found(false);                                                                      // A flag to check if the end of the header reached
    int exp_count(1);                                                                                     // Counter for every new exposure found. Default is 1.

    while(std::getline(Sfile, line)){                                                                     // Loop through each line in the file
        if (!line.empty()){                                                                               // Process only non-empty lines
            if(!end_of_header_found){
                if(line[0]=='*'){                                                                         // End of the header line starts with '*' character
                    end_of_header_found = true;
                }
                continue;
            }
            if(line[0] == '2'){                                                                           // New exposure entries starts with 2 in data field 1
                exp_count++;                                                                              // Increment the exposure counter
            }
        }
    }
    return exp_count;
}
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------
// A function to read any file that is formatted like an actual dose file with 1D data. This function will only read
// number of lines specified by the value lineCount and append these lines into a vector and returns the vector not a copy. 
// If the data entries in the file are fewer than the lineCount, print ERROR and exit 
//------------------------------------------------------------------------------------------------------------------------
std::vector<double> readActualDosefile(const std::string* dosefile, int lineCount){
    std::ifstream Dfile(*dosefile);                                                                       // Starting an ifstream instance with SDD filename
    std::string line;                                                                                     // Temporary variable to hold each line read
    std::vector<double> vec;                                                                              // Temporary vector to store the read lines (dose values)
    int missingEntries(lineCount);                                                                        // Temporary variable to hold the number of missing entries in file

    // Loop through each line in the file
    while(std::getline(Dfile, line)){
        
        if(!Dfile.eof() && !line.empty() && lineCount>=1){                                               // If the non-empty line is not the last line and we need dose values still
            vec.push_back(std::stod(line));                                                              // Append the dose value to the temporary vector
            lineCount--;                                                                                 // Decrement the number of lines to read
        }else if(Dfile.eof() && !line.empty() && lineCount==1){                                          // If it is the last non-empty line, but we only need one more data value (same line)
            vec.push_back(std::stod(line));                                                              // Append the dose value to the temporary vector
            lineCount--;                                                                                 // Decrement the number of lines to read
        }

        if(lineCount==0){break;}                                                                         // If enough dose values are already read, break the loop
        if(!line.empty()){missingEntries--;}                                                             // count the number of non-empty lines and subtract it from total required entries to get missing data count 
    }
    if(lineCount!=0){                                                                                    // If enough dose values are not present in the file, print error and exit   
        std::cerr<<"\n ERROR: The "<<*dosefile<<" is missing "<<missingEntries<<" more data entries than required \n";
        exit(EXIT_FAILURE); 
    }
    return std::move(vec);                                                                                // move() transfers the ownership of the memory to the caller of the function. This avoids copying the entire vector
}
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------
// A function to read any file that is formatted like an SDD file and get specific header fields 
//------------------------------------------------------------------------------------------------------------------------
void readSDDfileHeader(const std::string* sddfile, NGSsdd& SDDdata){
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
// A function to read any file that is formatted like an SDD file and get the data fields. This function will read
// one exposure data at a time. The second time the function is called, it will read the second exposure data and so on. 
//------------------------------------------------------------------------------------------------------------------------
void readSDDfileData(const std::string* sddfile, NGSsdd& SDDdata, int SDDfilenumber){
    std::ifstream Sfile(*sddfile);                                                                        // Starting an ifstream instance with SDD filename
    std::string line;                                                                                     // Temporary variable to hold each line read
    SDDdata.reset_original_num_damages();                                                                 // Reset the damage counter before the new exposure data or new SDD is read
    std::vector<std::streampos>& pos = *SDDdata.get_lines_read_sdd();
    bool end_of_header_found(true);                                                                       // A flag to check if the end of the header reached
    
    if (pos[SDDfilenumber]== 0){                                                                          // If the position is set to be the beggining of the file, then
        end_of_header_found = false;                                                                      // Set the flag to say that the header is not found yet
    }
    Sfile.seekg(pos[SDDfilenumber]);                                                                      // Get the position of the last line that was read for the SDD
    while(std::getline(Sfile, line)){
        if(!line.empty()){
            if(!end_of_header_found){                                                                     // Ignore everything before the end of header. Find header in the first pass
                if(line[0]=='*'){                                                                         // End of the header line starts with '*' character
                    end_of_header_found = true;
                }
                continue;
            }

            // If the line is after header and before the new exposure entries
            SDDdata.set_all_sdd_data_fields(&line);

            if(Sfile.peek() == '2'){                                                                      // Process only upto the beginning of a new exposure data
                SDDdata.set_lines_read_sdd(Sfile.tellg(), SDDfilenumber);                                 // Store the position of the new exposure data and return from function
                return;
            }
        }
    }
}
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------
// A function to read each chromosome sequence one-by-one from the reference sequence. It will get the chromosome IDs, 
// read the chromosome sequence, and merge multiple lines together. The string passed by reference will get the value of 
// the chromosome sequence. By calling this function in a loop, we can get all the chromosome seqs one-by-one
//------------------------------------------------------------------------------------------------------------------------
int getNextChromSeq(std::ifstream& ref_seq, std::string& chrom_seq, std::string& chromSeq_ID) {
    chrom_seq.clear();                                                                                    // Clear the previously stored value if any
    chromSeq_ID.clear();                                                                                  // Clear the previously stored value of sequence ID if any
    std::string tmpString;                                                                                // Temporary vector to hold each line value
    int success_flag{0};                                                                                  // This will return 0 if no more sequences left to read
    char tmpChar = '-';
   
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
            chrom_seq.append(tmpString);                                                                  // Merge multiple lines of chrom seq to one
            success_flag++;
        }
    }
    return(success_flag);                                                                                 // Return 0 if either reached end or if no sequence found
}
//------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------
// A function to read the undamaged genome template. Every odd and even sequences will be read simultaneously. Odd
// sequence corresponds to the forward strand sequence and the even ones corresponds to the reverse complementary strand sequence. 
// They will be stored to A and B tagged strings respectively
//------------------------------------------------------------------------------------------------------------------------
int readFastaTemplate(std::ifstream& fastaFile, std::string& chrom_ID_A, std::string& chrom_seq_A, std::string& chrom_ID_B, std::string& chrom_seq_B) {
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

        for(int i=0; i<cumulative_frequency.size(); i++){
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
               <<"                                   SDD-NGS SIMULATOR RUN SUMMARY REPORT             \n"
               <<"------------------------------------------------------------------------------------------------------------\n"
               <<" This report was prepared on "<<timeString<<"\n\n"
               <<" Time duration of this run                                            : "<<report_cpu_time_used<<" sec\n";
    if(parameter.get_random_seed() == 0){
        report_file<<" The seed used for this run                                           : Random seed\n";
    }else{
        report_file<<" The random seed used for this run                                    : "<<parameter.get_random_seed()<<"\n";
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
            if(i>parameter.get_names_of_particles_to_merge()->size()-1) report_file<<"Unspecified";
            else {
                const std::vector<std::string>* particle_names = parameter.get_names_of_particles_to_merge();
                report_file<<(*particle_names)[i];
            }
        }
        report_file<<"\n The relative dose contributions of these particles                   : ";
        for(int i=0;i<parameter.get_num_of_particles_to_merge();i++){
            if(i!=0) report_file<<", ";
            if(i>parameter.get_relative_dose_contributions().size()-1) report_file<<"0";
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
                   <<" Number of damages are adjusted for delievered dose using the files   : ";
        for(int i=0;i<parameter.get_num_of_particles_to_merge();i++){
            if(i!=0) report_file<<", ";
            const std::vector<std::string>* actual_dose_files = parameter.get_actual_dosefile_path();
            report_file<<(*actual_dose_files)[i];
        }
        report_file<<"\n";     
    }else report_file<<" Are the radiation-induced damages adjusted for actual delivered dose : No \n";

    report_file<<" The genome length of the Monte Carlo cell model                      : "<<SDDdata.get_sdd_genome_length()<<" bp\n"
               <<" Name of the reference genome file used                               : "<<*parameter.get_reference_genome()<<"\n"
               <<" The length of the reference gonome                                   : "<<report_ref_seq_length<<" bp\n\n"
               <<" ------ Sequencing information -------\n"
               <<" Name of the Illumina sequencer used for NGS simulation               : "<<*parameter.get_sequencer()<<"\n";
    if(*parameter.get_sequencing_mode() == "single"){
        report_file<<" The sequencing protocol used                                         : Single-cell Whole Genome Sequencing\n";
    }else report_file<<" The sequencing protocol used                                         : Bulk-cell Whole Genome Sequencing\n";
    
    report_file<<" Length of the reads simulated                                        : "<<parameter.get_read_length()<<" bp\n";

    if(parameter.get_paired_end_sequencing()){
        report_file<<" The read generation modality used                                    : Paired-end read generation\n"
                   <<" The mean length of DNA fragments created in the sample               : "<<parameter.get_mean_DNA_fragment_length()
                   <<" bp with std. deviation of "<<parameter.get_std_dev_DNA_fragment_length()<<" bp\n";
    }else report_file<<" The read generation modality used                                    : Single-end read generation\n";
    
    int total_reads = (parameter.get_total_read_coverage()*report_ref_seq_length)/parameter.get_read_length();
    report_file<<" Total number of cells in the sample                                  : "<<parameter.get_num_of_cells_in_sample()<<"\n"
               <<" Number of cells randomly sampled from the sample pool for sequencing : "<<parameter.get_num_of_cells_to_sequence()<<"\n"
               <<" Total number reads generated in the NGS simulation (expected)        : "<<total_reads<<"\n"
               <<" Average number of reads generated per cell                           : "<<total_reads/parameter.get_num_of_cells_to_sequence()<<"\n"
               <<" Maximum number of errors allowed per read                            : ";
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
        report_file<<" The mean read coverage per cell                                      : "<<parameter.get_total_read_coverage()/parameter.get_num_of_cells_to_sequence()<<"\n"
                   <<" The directory where the sequenced FASTQ files are stored             : "<<*parameter.get_output_directory()<<"\n\n"
                   <<" ------ Output data information -------\n\n"
                   <<"------------------------------------------------------------------------------------------------------------\n"
                   <<"|    Name of the cell sequenced                   |    Output FASTQ filename                               |\n"
                   <<"------------------------------------------------------------------------------------------------------------\n";
        for(int i=0; i<report_cells_sequenced.size(); i++){
            int cell_name_size = report_cells_sequenced[i].size();
            report_file<<std::left<<std::setw(50)<<"|    "+report_cells_sequenced[i].erase(cell_name_size-3)<<"|    "
                       <<std::left<<std::setw(52)<<report_fastq_output[i]+"_R1.fastq.gz"<<"|\n";
            if(parameter.get_paired_end_sequencing()){
                report_file<<std::left<<std::setw(50)<<"|    "<<"|    "<<std::left<<std::setw(52)<<report_fastq_output[i]+"_R2.fastq.gz"<<"|\n";
            }
            report_file<<"|                                                                                                          |\n";
        }
        report_file<<"------------------------------------------------------------------------------------------------------------\n";
    }else{ 
        report_file<<" Total read coverage obtained                                         : "<<parameter.get_total_read_coverage()<<"\n"
                   <<" The directory where the sequenced FASTQ files are stored             : "<<*parameter.get_output_directory()<<"\n\n"
                   <<" ------ Output data information -------\n\n"
                   <<" Name of the output FASTQ file                                        : "<<report_fastq_output[0]<<"_R1.fastq.gz";
        if(parameter.get_paired_end_sequencing()){
            report_file<<", "<<report_fastq_output[0]<<"_R2.fastq.gz";
        }
        report_file<<"\n------------------------------------------------------------------------------------------------------------\n"
                   <<"|    Name of the cell sequenced                   |    Percentage read contribution of the cell            |\n"
                   <<"------------------------------------------------------------------------------------------------------------\n";
        std::map<std::string, int> read_contribution;
        for (int i = 0; i < report_cells_sequenced.size(); i++) {
            read_contribution[report_cells_sequenced[i]] += report_readsGenerated_perCell[i];
        }
        for (const auto& entry : read_contribution) {
            std::string cell_name = entry.first;
            int cell_name_size = cell_name.size();
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