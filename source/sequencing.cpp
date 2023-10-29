#include "sequencing.h"
#include "fileio.h"
#include "art_framework.h"
#include "random_generator.h"
#include "summary_report.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <sys/mman.h>

//--------------------------------------------------------------------------------------------
// This is the function that will perform the single-cell sequencing of all the cells that are
// sampled from the pool of cells one-by-one. This function generates reads according to the
// coverage per cell for every cell. The generated reads will have indel errors as well as
// base call errors based on the sequencer's error quality profile. The read data will be 
// stored to fastq.gz files that correspond to each cell.
//--------------------------------------------------------------------------------------------
void single_cell_sequencing(NGSParameters& parameter, const std::vector<std::string>& cellGenomes_to_be_sequenced, long ref_seq_length){
    double coverage_per_cell = parameter.get_total_read_coverage()/parameter.get_num_of_cells_to_sequence();
    ART::set_read_quality_distribution(*parameter.get_r1_quality_profile(), *parameter.get_r2_quality_profile());
    ART::set_baseCall_error_probability();                                                              // Generate base call error probabilty distribution (same for all cells)

    for (int i=0; i<cellGenomes_to_be_sequenced.size(); i++){                                           // Iterate though each cell genome that is to be sequenced to generate reads
        
        std::string cell_fasta_filename = (*parameter.get_output_directory())+"/temp/"+cellGenomes_to_be_sequenced[i];
        if (checkFileExists(&cell_fasta_filename)){;                                                    // Check if the fasta file of the cell to be sequenced can be read
            std::cout<<"\n Sequencing of cell "<<i+1<<" is in progress \n";
        }else{
            std::cerr<<"\n WARNING: Unable to sequence cell "<<i+1<<". This cell will be ignored.\n";
            continue;
        }

        size_t position{0};                                                                             // Temporary variable to hold the last-read position in the memory map of the current fasta file
        size_t fastaFileSize;                                                                           // A variable to hold the size of the fasta file of the current cell, during memory-mapping
        void* fastaFileMM = generateInputFileMemoryMap(cell_fasta_filename, fastaFileSize);             // Create the memory-map of the fasta file
        const char* cell_fastaFileData = static_cast<char*>(fastaFileMM);                               // Casting the memory-map void pointer to a const char pointer for further processing

        std::string chromSegmentSeq;                                                                    // Temporary variable to hold each chromosome segment sequence from the fasta file one at a time
        std::string chromSegmentSeq_ID;                                                                 // Temporary variable to hold IDs of each hromosome segment sequence

        ART read1;                                                                                      // Creating an ART class object and setting the insertion and deletion probability vectors for that read object
        read1.set_read_error_probability(parameter.get_read_length(), parameter.get_insertion_error_rate_read1(), read1.insertion_probability_vec, parameter.get_max_errors_in_read());
        read1.set_read_error_probability(parameter.get_read_length(), parameter.get_deletion_error_rate_read1(), read1.deletion_probability_vec, parameter.get_max_errors_in_read());

        std::vector<short> read1_quality_score_vec;                                                     // Vector to hold the quality scores for read 1

        std::string output_fastq_R1_filename = (*parameter.get_output_directory())+"/"+(*parameter.get_output_fastq_filename_prefix())+"_"+std::to_string(i)+"_R1.fastq.gz";
        std::ofstream fastq_R1_file(output_fastq_R1_filename.c_str(),std::ios::binary);                 // ofstream object of the output fastq file for read 1
        
        const int batchSize{4000};                                                                      // Define a batch size for writing reads to the output file. These much data will be stored in cache before writing it on the file

        // Single-end sequencing
        if(!parameter.get_paired_end_sequencing()){                                                     // If user asked for single-end (not paired-end) sequencing
            std::vector<std::string> batch_buffer;                                                      // Create a buffer for storing read data.
            while(getNextChromSeq_MM(cell_fastaFileData, fastaFileSize, position, chromSegmentSeq, chromSegmentSeq_ID)){// Iterate through each chromsome segment sequence in the fasta file memory map of the cell that is being sequenced
                long num_reads_per_segment = static_cast<long>((coverage_per_cell*chromSegmentSeq.size())/parameter.get_read_length());
                if(!read1.set_chromSegmentSeq(parameter.get_read_length(), chromSegmentSeq)){           // If chromSegmentSeq is smaller than the read length, then continue with the next segment
                    continue;
                };
                
                while(num_reads_per_segment>0){
                    read1.generate_read_with_indel();                                                   // Make a read with random indel errors
                    read1_quality_score_vec.clear();
                    read1.get_read_quality(read1_quality_score_vec, 1);                                 // Get the read quality scores for the read positions
                    read1.add_baseCall_error(read1_quality_score_vec);                                  // Add base call errors to the read based on the quality scores
                    
                    std::string read_data = "@"+chromSegmentSeq_ID+"_read"+std::to_string(num_reads_per_segment)+"\n"; // @readID
                    read_data += (*read1.get_final_read_sequence())+ "\n+\n";                                          // read sequence and +
                    for(int k=0; k<read1_quality_score_vec.size(); k++){                                               // read quality scores
					    read_data += static_cast<char>(read1_quality_score_vec[k]+33);
				    }
                    read_data += "\n";
                    
                    num_reads_per_segment--;
                    
                    batch_buffer.push_back(read_data);                                                  // Add the read data to the batch buffer.                    
                    if (batch_buffer.size() >= batchSize) {                                             // Check if the batch buffer is full, and write it to the file if needed.
                        writeBatchToFile(batch_buffer, fastq_R1_file);
                    }
                }
            }
            writeBatchToFile(batch_buffer, fastq_R1_file);                                              // If there are unwritten data in batch buffer, write that too when the loop ends
        }
        // Paired-end sequencing
        else{                                                                                           // If user asked to perform paired-end sequencing
            ART read2;                                                                                  // Creating an ART class object and setting the insertion and deletion probability vectors for that read object
            read2.set_read_error_probability(parameter.get_read_length(), parameter.get_insertion_error_rate_read1(), read2.insertion_probability_vec, parameter.get_max_errors_in_read());
            read2.set_read_error_probability(parameter.get_read_length(), parameter.get_deletion_error_rate_read1(), read2.deletion_probability_vec, parameter.get_max_errors_in_read());

            std::vector<short> read2_quality_score_vec;                                                 // Vector to hold the quality scores for read 2

            std::string output_fastq_R2_filename = (*parameter.get_output_directory())+"/"+(*parameter.get_output_fastq_filename_prefix())+"_"+std::to_string(i)+"_R2.fastq.gz";
            std::ofstream fastq_R2_file(output_fastq_R2_filename.c_str(),std::ios::binary);             // ofstream object of the output fastq file for read 2

            std::vector<std::string> batch_buffer_r1;                                                   // Create a buffer for storing read 1 data.
            std::vector<std::string> batch_buffer_r2;                                                   // Create a buffer for storing read 2 data.
    
            while(getNextChromSeq_MM(cell_fastaFileData, fastaFileSize, position, chromSegmentSeq, chromSegmentSeq_ID)){// Iterate through each chromsome segment sequence in the fasta file memory map of the cell that is being sequenced
                long num_reads_per_segment = static_cast<long>((coverage_per_cell*chromSegmentSeq.size())/parameter.get_read_length());
                if(!read1.set_chromSegmentSeq(parameter.get_read_length(), chromSegmentSeq)){           // If chromSegmentSeq is smaller than the read length, then continue with the next segment
                    continue;
                };
                read2.set_chromSegmentSeq(parameter.get_read_length(), chromSegmentSeq);                // Giving a copy of these for read 2 as well
                
                while(num_reads_per_segment>0){
                    
                    ART::generate_paired_reads_with_indel(read1, read2, parameter.get_mean_DNA_fragment_length(), parameter.get_std_dev_DNA_fragment_length());
                                                                                                        // Make two paired-reads from the same DNA fragment with indel errors
                    // Process read 1 first
                    read1_quality_score_vec.clear();
                    read1.get_read_quality(read1_quality_score_vec, 1);                                 // Get the read quality scores for the read positions for read 1
                    read1.add_baseCall_error(read1_quality_score_vec);                                  // Add base call errors to the read based on the quality scores on read 2
                    
                    std::string read1_data = "@"+chromSegmentSeq_ID+"_read"+std::to_string(num_reads_per_segment)+"\n"; // @readID
                    read1_data += (*read1.get_final_read_sequence()) + "\n+\n";                                         // read sequence and +
                    for(int k=0; k<read1_quality_score_vec.size(); k++){                                                // read quality scores
					    read1_data += static_cast<char>(read1_quality_score_vec[k]+33);
				    }
                    read1_data += "\n";

                    // Process read 2
                    read2_quality_score_vec.clear();
                    read2.get_read_quality(read2_quality_score_vec, 2);                                 // Get the read quality scores for the read positions for read 2
                    read2.add_baseCall_error(read2_quality_score_vec);                                  // Add base call errors to the read based on the quality scores on read 2
                    
                    std::string read2_data = "@"+chromSegmentSeq_ID+"_read"+std::to_string(num_reads_per_segment)+"\n"; // @readID
                    read2_data += (*read2.get_final_read_sequence()) + "\n+\n";                                         // read sequence and +
                    for(int k=0; k<read2_quality_score_vec.size(); k++){                                                // read quality scores
					    read2_data += static_cast<char>(read2_quality_score_vec[k]+33);
				    }
                    read2_data += "\n";
                    num_reads_per_segment -= 2;

                    batch_buffer_r1.push_back(read1_data);                                              // Add the read 1 data to the batch buffer 1.   
                    batch_buffer_r2.push_back(read2_data);                                              // Add the read 2 data to the batch buffer 2.                 
                    if (batch_buffer_r1.size() >= batchSize) {                                          // Check if the batch buffer is full, and write it to the file if needed.
                        writeBatchToFile(batch_buffer_r1, fastq_R1_file);
                        writeBatchToFile(batch_buffer_r2, fastq_R2_file);
                    }
                }
            }
            writeBatchToFile(batch_buffer_r1, fastq_R1_file);                                           // If there are unwritten data in batch buffer 1, write that too when the loop ends
            writeBatchToFile(batch_buffer_r2, fastq_R2_file);                                           // If there are unwritten data in batch buffer 1, write that too when the loop ends

            fastq_R2_file.close();
        }
        fastq_R1_file.close();

        report_cells_sequenced.push_back(cellGenomes_to_be_sequenced[i]);                               // To make the report, making a vector of all the cells that were sequenced
        report_fastq_output.push_back((*parameter.get_output_fastq_filename_prefix())+"_"+std::to_string(i));

        munmap(fastaFileMM, fastaFileSize);                                                             // Unmap the fasta file of the current cell to avoid memory-leaks
    }  
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This is the function that will perform the Bulk-cell sequencing of all the cells that are
// sampled from the pool of cells. This function will randomly pick one cell from the
// sample, then it will randomly pick a chromosome fragment and randomly generate one read. 
// The process will be repeated until enough reads are generated randomly.
//--------------------------------------------------------------------------------------------
void bulk_cell_sequencing(NGSParameters& parameter, const std::vector<std::string>& cellGenomes_to_be_sequenced, const std::vector<int>& lines_in_cells_to_sequence, long ref_seq_length ){
    int total_num_reads = std::round((parameter.get_total_read_coverage()*ref_seq_length)/parameter.get_read_length());
    int max_num_reads_per_cell = std::round(total_num_reads/parameter.get_num_of_cells_to_sequence());  // Temporary variable to hold the maximum number of reads per cell
    
    ART::set_read_quality_distribution(*parameter.get_r1_quality_profile(), *parameter.get_r2_quality_profile());
    ART::set_baseCall_error_probability();                                                              // Generate base call error probabilty distribution (same for all cells)
    
    ART read1;                                                                                          // Creating an ART class object and setting the insertion and deletion probability vectors for that read object
    read1.set_read_error_probability(parameter.get_read_length(), parameter.get_insertion_error_rate_read1(), read1.insertion_probability_vec, parameter.get_max_errors_in_read());
    read1.set_read_error_probability(parameter.get_read_length(), parameter.get_deletion_error_rate_read1(), read1.deletion_probability_vec, parameter.get_max_errors_in_read());

    std::vector<short> read1_quality_score_vec;                                                         // Vector to hold the quality scores for read 1

    std::string output_fastq_R1_filename = (*parameter.get_output_directory())+"/"+(*parameter.get_output_fastq_filename_prefix())+"_R1.fastq.gz";
    std::ofstream fastq_R1_file(output_fastq_R1_filename.c_str(),std::ios::binary);                     // ofstream object of the output fastq file for read 1

    // Generate a list of lines that we want to read randomly from all the cell file prior to reading each file. This is done to avoid processing the same fasta files multiple times to generate reads from it one by one. 
    // Here we generate a map, that indicates which all lines from one fasta files should be read so that we call read all these lines in one go. 
    int readPerFrag{1};                                                                                 // A variable to indicate how many reads per fragments we need to create. It will get 1 if single-end sequencing
    if(parameter.get_paired_end_sequencing()){
        readPerFrag = 2;                                                                                // This variable gets 2 if we want paired-end sequencing
    }
    std::unordered_map<std::string, std::vector<int>> rand_lines_map;                                   // A map to hold the randoms lines in each cell files that we want to read to generate the reads
    while(total_num_reads>0){
        int cell_index = rng::rand_int(1,cellGenomes_to_be_sequenced.size()) - 1;                       // Randomly choose a cell to generate reads. 'cell_index' can take values from 0 to cellGenomes_to_be_sequenced.size minus 1
        int reads_to_generate_in_cell = rng::rand_int(0,max_num_reads_per_cell);                        // Randomly determine how many reads should be generated from this particular cell that is selected randomly
        if(max_num_reads_per_cell < readPerFrag ){                                                      // If max_num_reads_per_cell is zero, this will prevent program getting stuck
            reads_to_generate_in_cell = rng::rand_int(0,total_num_reads);
        } 
        while(reads_to_generate_in_cell>total_num_reads){                                               // If the number of reads to generate exceeds the number of reads remaining to make in the current loop, then regenerate the number
           reads_to_generate_in_cell = rng::rand_int(0,total_num_reads);
        }               
        
        int reads_will_be_generated_in_cell{0};                                                         // A temporary variable to hold the number of reads that will actually be created from each cell. For single-end seq this number will be equal to reads_to_generate_in_cell
        int totalLines_in_fasta = lines_in_cells_to_sequence[cell_index];                               // Variable to hold the number of lines in the fasra file of the cell selected
        while (reads_to_generate_in_cell > (readPerFrag-1)){                                            // Randomly generate a list of line numbers that we want to use pick the chrom fragments from the fasta file
            int line_number = rng::rand_int(1,totalLines_in_fasta);                                     // line_number can take any value from 0 to length of the fasta file
            if (line_number %2 == 0){                                                                   // If the generated line number is even, it corresponds to the sequence line, not sequence ID line
                line_number -= 1;                                                                       // Store the line number corresponding to the sequence ID instead. If the line number is odd, it is already corresponding to the sequence ID postiontion in fasta file
            }
            auto it = rand_lines_map.find(cellGenomes_to_be_sequenced[cell_index]);                     // Check in the map if the file corresponding to the cell was previously randomly selected and if in the map already
            if (it != rand_lines_map.end()) {                                                           // If cell_index exists in the map, append the line number corresponding to it's file
                it->second.push_back(line_number);
            } else {                                                                                    // If cell_index doesn't exist, create a new entry in the map for the corresponding cell file
                rand_lines_map[cellGenomes_to_be_sequenced[cell_index]] = {line_number};
            }
            reads_to_generate_in_cell -= readPerFrag;
            reads_will_be_generated_in_cell += readPerFrag;
        }

        total_num_reads -= reads_will_be_generated_in_cell;
    }  
 
    const int batchSize{4000};                                                                          // Define a batch size for writing reads to the output file. These much data will be stored in cache before writing it on the file

    // Single-end sequencing
    if(!parameter.get_paired_end_sequencing()){
        std::vector<std::string> batch_buffer;                                                          // Create a buffer for storing read data.
        
        for (const auto& mapEntry : rand_lines_map) {                                                   // Iterate through each entries in the map. entry.first will be the cell index and entry.second is the list of the lines we need to read from that cell
            std::string cell_fasta_filename = (*parameter.get_output_directory())+"/temp/"+mapEntry.first;
            size_t position{0};                                                                         // Temporary variable to hold the last-read position in the memory map of the current fasta file
            size_t fastaFileSize;                                                                       // A variable to hold the size of the fasta file of the current cell, during memory-mapping
            void* fastaFileMM = generateInputFileMemoryMap(cell_fasta_filename, fastaFileSize);         // Create the memory-map of the fasta file
            const char* cell_fastaFileData = static_cast<char*>(fastaFileMM);                           // Casting the memory-map void pointer to a const char pointer for further processing

            std::string chromSegmentSeq;                                                                // Temporary variable to hold each chromosome segment sequence from the fasta file one at a time
            std::string chromSegmentSeq_ID;                                                             // Temporary variable to hold each chromsome segment ID
            
            int line_count{1};                                                                          // Temporary variable to count the number of lines read from the fasta files using the getNextChromSeq function
            int reads_actually_generated{0};                                                            // Temporary counter variable to count the number of reads actually generated for each cell; for the report
            
            std::vector<int> rand_lines_list = mapEntry.second;                                         // Copy the lines list from the map into a vector for processing
            std::sort(rand_lines_list.begin(), rand_lines_list.end());                                  // Sort the line numbers in the vector in ascending order

            while(getNextChromSeq_MM(cell_fastaFileData, fastaFileSize, position, chromSegmentSeq, chromSegmentSeq_ID)){// Iterate through each chromsome segment sequence in the fasta file memory map of the cell that is being sequenced
                if(rand_lines_list.empty()) break;                                                      // If no more reads to be generated from the current cell, then break the loop
                
                if(!read1.set_chromSegmentSeq(parameter.get_read_length(), chromSegmentSeq)){           // If chromSegmentSeq is smaller than the read length, then continue with the next segment
                    continue;
                };
                
                std::vector<int>::iterator it;                                                          // Temporary iterator object for the rand_lines_list vector
                while ((it = std::find(rand_lines_list.begin(), rand_lines_list.end(), line_count)) != rand_lines_list.end()) {
                    read1.generate_read_with_indel();                                                   // Make a read with random indel (insertions and deletions) errors
                    read1_quality_score_vec.clear();
                    read1.get_read_quality(read1_quality_score_vec, 1);                                 // Get the read quality scores for the read positions
                    read1.add_baseCall_error(read1_quality_score_vec);                                  // Add base call errors to the read based on the quality scores
                    
                    std::string read_data = "@"+mapEntry.first+"_"+chromSegmentSeq_ID+"_read"+std::to_string(reads_actually_generated)+"\n";  // @readID
                    read_data += (*read1.get_final_read_sequence())+ "\n+\n";                                                                 // read sequence and +
                    for(int k=0; k<read1_quality_score_vec.size(); k++){                                                                      // read quality scores
					    read_data +=static_cast<char>(read1_quality_score_vec[k]+33);
				    }
                    read_data +="\n";
                
                    rand_lines_list.erase(it);                                                          // Remove the line from the list after making a read
                    total_num_reads--;                                                                  // When a read is generated, decrement the count
                    reads_actually_generated++;
                    
                    batch_buffer.push_back(read_data);                                                  // Add the read data to the batch buffer.                    
                    if (batch_buffer.size() >= batchSize) {                                             // Check if the batch buffer is full, and write it to the file if needed.
                        for (const std::string& read : batch_buffer) {
                            fastq_R1_file << read;
                        }
                        batch_buffer.clear();
                    }
                }
                line_count += 2;                                                                        // Since two lines (ID and sequence) were read from the fasta file, increment line count by 2
            }
            munmap(fastaFileMM, fastaFileSize);                                                         // Unmap the fasta file of the current cell to avoid memory-leaks
    
            report_cells_sequenced.push_back(mapEntry.first);                                           // Storing the names of the cells that were actually sequenced, for the report
            report_readsGenerated_perCell.push_back(reads_actually_generated);                          // Storing the number of reads generated per cell in a vector for the final summary report
        }
        for (const std::string& read : batch_buffer) {                                                  // If there are unwritten data in batch buffer, write that too when the loop ends
            fastq_R1_file << read;
        } 
    }
    // Paired-end sequencing
    else{
        ART read2;                                                                                      // Creating an ART class object and setting the insertion and deletion probability vectors for that read object
        read2.set_read_error_probability(parameter.get_read_length(), parameter.get_insertion_error_rate_read1(), read2.insertion_probability_vec, parameter.get_max_errors_in_read());
        read2.set_read_error_probability(parameter.get_read_length(), parameter.get_deletion_error_rate_read1(), read2.deletion_probability_vec, parameter.get_max_errors_in_read());

        std::vector<short> read2_quality_score_vec;                                                     // Vector to hold the quality scores for read 2

        std::string output_fastq_R2_filename = (*parameter.get_output_directory())+"/"+(*parameter.get_output_fastq_filename_prefix())+"_R2.fastq.gz";
        std::ofstream fastq_R2_file(output_fastq_R2_filename.c_str(),std::ios::binary);                 // ofstream object of the output fastq file for read 2

        std::vector<std::string> batch_buffer_r1;                                                       // Create a buffer for storing read 1 data.
        std::vector<std::string> batch_buffer_r2;                                                       // Create a buffer for storing read 2 data.
    
        for (const auto& mapEntry : rand_lines_map) {                                                   // Iterate through each entries in the map. entry.first will be the cell index and entry.second is the list of the lines we need to read from that cell
            std::string cell_fasta_filename = (*parameter.get_output_directory())+"/temp/"+mapEntry.first;
            size_t position{0};                                                                         // Temporary variable to hold the last-read position in the memory map of the current fasta file
            size_t fastaFileSize;                                                                       // A variable to hold the size of the fasta file of the current cell, during memory-mapping
            void* fastaFileMM = generateInputFileMemoryMap(cell_fasta_filename, fastaFileSize);         // Create the memory-map of the fasta file
            const char* cell_fastaFileData = static_cast<char*>(fastaFileMM);                           // Casting the memory-map void pointer to a const char pointer for further processing

            std::string chromSegmentSeq;                                                                // Temporary variable to hold each chromosome segment sequence from the fasta file one at a time
            std::string chromSegmentSeq_ID;                                                             // Temporary variable to hold each chromsome segment ID

            int line_count{1};                                                                          // Temporary vector to count the number of lines read from the fasta files using the getNextChromSeq function
            int reads_actually_generated{0};                                                            // Temporary counter variable to count the number of reads actually generated for each cell; for the report
            
            std::vector<int> rand_lines_list = mapEntry.second;                                         // Copy the lines list from the map into a vector for processing
            std::sort(rand_lines_list.begin(), rand_lines_list.end());                                  // Sort the line numbers in the vector in ascending order

            while(getNextChromSeq_MM(cell_fastaFileData, fastaFileSize, position, chromSegmentSeq, chromSegmentSeq_ID)){// Iterate through each chromsome segment sequence in the fasta file memory map of the cell that is being sequenced
                if(rand_lines_list.empty()) break;                                                      // If no more reads to be generated from the current cell, then break the loop
                
                if(!read1.set_chromSegmentSeq(parameter.get_read_length(), chromSegmentSeq)){           // If chromSegmentSeq is smaller than the read length, then continue with the next segment
                    continue;
                };
                read2.set_chromSegmentSeq(parameter.get_read_length(), chromSegmentSeq);                // Giving a copy of these for read 2 as well

                std::vector<int>::iterator it;                                                          // Temporary iterator object for the rand_lines_list vector
                while ((it = std::find(rand_lines_list.begin(), rand_lines_list.end(), line_count)) != rand_lines_list.end()) {
                    ART::generate_paired_reads_with_indel(read1, read2, parameter.get_mean_DNA_fragment_length(), parameter.get_std_dev_DNA_fragment_length());
                                                                                                        // Make two paired-reads from the same DNA fragment with indel errors
                    // Process read 1 first
                    read1_quality_score_vec.clear();
                    read1.get_read_quality(read1_quality_score_vec, 1);                                 // Get the read quality scores for the read positions for read 1
                    read1.add_baseCall_error(read1_quality_score_vec);                                  // Add base call errors to the read based on the quality scores on read 2
                    
                    std::string read1_data ="@"+mapEntry.first+"_"+chromSegmentSeq_ID+"_read"+std::to_string(reads_actually_generated)+"\n";   // @readID
                    read1_data +=(*read1.get_final_read_sequence())+"\n+\n";                                                                   // read sequence and +
                    for(int k=0; k<read1_quality_score_vec.size(); k++){                                                                       // read quality scores
					    read1_data +=static_cast<char>(read1_quality_score_vec[k]+33);
				    }
                    read1_data +="\n";

                    // Process read 2
                    read2_quality_score_vec.clear();
                    read2.get_read_quality(read2_quality_score_vec, 2);                                 // Get the read quality scores for the read positions for read 2
                    read2.add_baseCall_error(read2_quality_score_vec);                                  // Add base call errors to the read based on the quality scores on read 2
                    
                    std::string read2_data ="@"+mapEntry.first+"_"+chromSegmentSeq_ID+"_read"+std::to_string(reads_actually_generated+1)+"\n";   // @readID
                    read2_data +=(*read2.get_final_read_sequence())+"\n+\n";                                                                     // read sequence and +
                    for(int k=0; k<read2_quality_score_vec.size(); k++){                                                                         // read quality scores
					    read2_data +=static_cast<char>(read2_quality_score_vec[k]+33);
				    }
                    read2_data +="\n";
                    
                    rand_lines_list.erase(it);                                                          // Remove the line from the list after making a read
                    total_num_reads -= 2;                                                               // When paired reads (2) are generated, decrement the count
                    reads_actually_generated += 2;
                
                    batch_buffer_r1.push_back(read1_data);                                              // Add the read 1 data to the batch buffer 1.   
                    batch_buffer_r2.push_back(read2_data);                                              // Add the read 2 data to the batch buffer 2.                 
                    if (batch_buffer_r1.size() >= batchSize) {                                          // Check if the batch buffer is full, and write it to the file if needed.
                        for (const std::string& read1 : batch_buffer_r1) {
                            fastq_R1_file << read1;
                        }
                        batch_buffer_r1.clear();

                        for (const std::string& read2 : batch_buffer_r2) {
                            fastq_R2_file << read2;
                        }
                        batch_buffer_r2.clear();
                    }
                }
                line_count += 2;                                                                        // Since two lines (ID and sequence) were read from the fasta file, increment line count by 2
            } 
            munmap(fastaFileMM, fastaFileSize);                                                         // Unmap the fasta file of the current cell to avoid memory-leaks
    
            report_cells_sequenced.push_back(mapEntry.first);                                           // Storing the names of the cells that were actually sequenced, for the report
            report_readsGenerated_perCell.push_back(reads_actually_generated);                          // Storing the number of reads generated per cell in a vector for the final summary report
        }
     
        for (const std::string& read1 : batch_buffer_r1) {                                              // If there are unwritten data in batch buffer 1, write that too when the loop ends
            fastq_R1_file << read1;
        }  
        for (const std::string& read2 : batch_buffer_r2) {                                              // If there are unwritten data in batch buffer 1, write that too when the loop ends
            fastq_R2_file << read2;
        }

        fastq_R2_file.close();
    }
    fastq_R1_file.close();

    report_fastq_output.push_back(*parameter.get_output_fastq_filename_prefix());                       // Output filename without the extension for the report
}
//--------------------------------------------------------------------------------------------