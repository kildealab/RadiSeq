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

        std::ifstream cell_fasta_file(cell_fasta_filename);                                             // ifstream object of the fasta file to be sequenced
        std::string chromSegmentSeq;                                                                    // Temporary variable to hold each chromosome segment sequence from the fasta file one at a time
        std::string chromSegmentSeq_ID;                                                                 // Temporary variable to hold IDs of each hromosome segment sequence

        ART read1;                                                                                      // Creating an ART class object and setting the insertion and deletion probability vectors for that read object
        read1.set_read_error_probability(parameter.get_read_length(), parameter.get_insertion_error_rate_read1(), read1.insertion_probability_vec, parameter.get_max_errors_in_read());
        read1.set_read_error_probability(parameter.get_read_length(), parameter.get_deletion_error_rate_read1(), read1.deletion_probability_vec, parameter.get_max_errors_in_read());

        std::vector<short> read1_quality_score_vec;                                                     // Vector to hold the quality scores for read 1

        std::string output_fastq_R1_filename = (*parameter.get_output_directory())+"/"+(*parameter.get_output_fastq_filename_prefix())+"_"+std::to_string(i)+"_R1.fastq.gz";
        std::ofstream fastq_R1_file(output_fastq_R1_filename.c_str(),std::ios::binary);                 // ofstream object of the output fastq file for read 1
        

        // Single-end sequencing
        if(!parameter.get_paired_end_sequencing()){                                                     // If user asked for single-end (not paired-end) sequencing
            while(getNextChromSeq(cell_fasta_file, chromSegmentSeq, chromSegmentSeq_ID)){               // Iterate through each chromsome segment sequence in the fasta file of the cell that is being sequenced
                long num_reads_per_segment = static_cast<long>((coverage_per_cell*chromSegmentSeq.size())/parameter.get_read_length());
                if(!ART::init_set(parameter.get_read_length(), chromSegmentSeq)){                       // If chromSegmentSeq is smaller than the read length, then continue with the next segment
                    continue;
                };
            
                while(num_reads_per_segment>0){
                    read1.generate_read_with_indel();                                                   // Make a read with random indel errors
                    read1_quality_score_vec.clear();
                    read1.get_read_quality(read1_quality_score_vec, 1);                                 // Get the read quality scores for the read positions
                    read1.add_baseCall_error(read1_quality_score_vec);                                  // Add base call errors to the read based on the quality scores
                    
                    fastq_R1_file<<"@"<<chromSegmentSeq_ID<<"_read"<<std::to_string(num_reads_per_segment)<<"\n" // @readID
                                 <<(*read1.get_final_read_sequence())                                            // read sequence
                                 <<"\n+\n";                                                                      // +
                    for(int k=0; k<read1_quality_score_vec.size(); k++){                                         // read quality scores
					    fastq_R1_file<<static_cast<char>(read1_quality_score_vec[k]+33);
				    }
                    fastq_R1_file<<"\n";
                    
                    num_reads_per_segment--;
                }
            }
        }
        // Paired-end sequencing
        else{                                                                                           // If user asked to perform paired-end sequencing
            ART read2;                                                                                  // Creating an ART class object and setting the insertion and deletion probability vectors for that read object
            read2.set_read_error_probability(parameter.get_read_length(), parameter.get_insertion_error_rate_read1(), read2.insertion_probability_vec, parameter.get_max_errors_in_read());
            read2.set_read_error_probability(parameter.get_read_length(), parameter.get_deletion_error_rate_read1(), read2.deletion_probability_vec, parameter.get_max_errors_in_read());

            std::vector<short> read2_quality_score_vec;                                                 // Vector to hold the quality scores for read 2

            std::string output_fastq_R2_filename = (*parameter.get_output_directory())+"/"+(*parameter.get_output_fastq_filename_prefix())+"_"+std::to_string(i)+"_R2.fastq.gz";
            std::ofstream fastq_R2_file(output_fastq_R2_filename.c_str(),std::ios::binary);             // ofstream object of the output fastq file for read 2

            while(getNextChromSeq(cell_fasta_file, chromSegmentSeq, chromSegmentSeq_ID)){               // Iterate through each chromsome segment sequence in the fasta file of the cell that is being sequenced
                long num_reads_per_segment = static_cast<long>((coverage_per_cell*chromSegmentSeq.size())/parameter.get_read_length());
                if(!ART::init_set(parameter.get_read_length(), chromSegmentSeq)){                       // If chromSegmentSeq is smaller than the read length, then continue with the next segment
                    continue;
                };
            
                while(num_reads_per_segment>0){
                    
                    ART::generate_paired_reads_with_indel(read1, read2, parameter.get_mean_DNA_fragment_length(), parameter.get_std_dev_DNA_fragment_length());
                                                                                                        // Make two paired-reads from the same DNA fragment with indel errors
                    // Process read 1 first
                    read1_quality_score_vec.clear();
                    read1.get_read_quality(read1_quality_score_vec, 1);                                 // Get the read quality scores for the read positions for read 1
                    read1.add_baseCall_error(read1_quality_score_vec);                                  // Add base call errors to the read based on the quality scores on read 2
                    
                    fastq_R1_file<<"@"<<chromSegmentSeq_ID<<"_read"<<std::to_string(num_reads_per_segment)<<"\n" // @readID
                                 <<(*read1.get_final_read_sequence())                                            // read sequence
                                 <<"\n+\n";                                                                      // +
                    for(int k=0; k<read1_quality_score_vec.size(); k++){                                         // read quality scores
					    fastq_R1_file<<static_cast<char>(read1_quality_score_vec[k]+33);
				    }
                    fastq_R1_file<<"\n";

                    // Process read 2
                    read2_quality_score_vec.clear();
                    read2.get_read_quality(read2_quality_score_vec, 2);                                 // Get the read quality scores for the read positions for read 2
                    read2.add_baseCall_error(read2_quality_score_vec);                                  // Add base call errors to the read based on the quality scores on read 2
                    
                    fastq_R2_file<<"@"<<chromSegmentSeq_ID<<"_read"<<std::to_string(num_reads_per_segment)<<"\n" // @readID
                                 <<(*read2.get_final_read_sequence())                                            // read sequence
                                 <<"\n+\n";                                                                      // +
                    for(int k=0; k<read2_quality_score_vec.size(); k++){                                         // read quality scores
					    fastq_R2_file<<static_cast<char>(read2_quality_score_vec[k]+33);
				    }
                    fastq_R2_file<<"\n";
                    num_reads_per_segment -= 2;
                }
            }
            fastq_R2_file.close();
        }
        fastq_R1_file.close();

        report_cells_sequenced.push_back(cellGenomes_to_be_sequenced[i]);                              // To make the report, making a vector of all the cells that were sequenced
        report_fastq_output.push_back((*parameter.get_output_fastq_filename_prefix())+"_"+std::to_string(i));

        cell_fasta_file.close();
    }  
}
//--------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
// This is the function that will perform the Bulk-cell sequencing of all the cells that are
// sampled from the pool of cells. This function will randomly pick one cell from the
// sample, then it will randomly pick a chromosome fragment and randomly generate one read. 
// The process will be repeated until enough reads are generated randomly.
//--------------------------------------------------------------------------------------------
void bulk_cell_sequencing(NGSParameters& parameter, const std::vector<std::string>& cellGenomes_to_be_sequenced, long ref_seq_length ){
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
     
    // Single-end sequencing
    if(!parameter.get_paired_end_sequencing()){
        while(total_num_reads>0){
            int cell_index = rng::rand_int(1,cellGenomes_to_be_sequenced.size()) - 1;                   // Randomly choose a cell to generate reads. 'cell' can take values from 0 to cellGenomes_to_be_sequenced.size minus 1
            std::string cell_fasta_filename = (*parameter.get_output_directory())+"/temp/"+cellGenomes_to_be_sequenced[cell_index];

            int reads_to_generate = rng::rand_int(0,max_num_reads_per_cell);                            // Randomly determine how many reads should be generated from this particular cell that is selected randomly
            if(max_num_reads_per_cell == 0 ){                                                           // If max_num_reads_per_cell is zero, this will prevent program getting stuck
                reads_to_generate = rng::rand_int(0,total_num_reads);
            }
            while(reads_to_generate>total_num_reads){                                                   // If the number of reads to generate exceeds the number of reads remaining to make in the current loop, then regenerate the number
                reads_to_generate = rng::rand_int(0,total_num_reads);
            }
            
            std::ifstream cell_fasta_file(cell_fasta_filename);                                         // ifstream object of the fasta file to be sequenced
            std::string chromSegmentSeq;                                                                // Temporary variable to hold each chromosome segment sequence from the fasta file one at a time
            std::string chromSegmentSeq_ID;                                                             // Temporary variable to hold each chromsome segment ID
            
            int totalLines_in_fasta = 0;                                                                // Variable to hold the number of lines in the file file of the cell selected
            while (std::getline(cell_fasta_file, chromSegmentSeq)) {                                    // Counting the number of lines in the fasta file for the cell
                totalLines_in_fasta++;
            }
            cell_fasta_file.clear();                                                                    // Reset the error flags of the stream so that subsequent operations can be performed without any issues
            cell_fasta_file.seekg(0);                                                                   // Seeking back to the beginning of the file, ensuring that the stream is in a valid state

            std::vector<int> rand_lines_list;                                                           // Vector to hold the fasta file lines that we want to use to generate reads at random
            while (reads_to_generate > 0){                                                              // Randomly generate a list of line numbers that we want to use pick the chrom fragments from the fasta file
                int line_number = rng::rand_int(1,totalLines_in_fasta);                                 // line_number can take any value from 0 to length of the fasta file
                if (line_number %2 == 0){                                                               // If the generated line number is even, it corresponds to the sequence line, not sequence ID line
                    rand_lines_list.push_back(line_number-1);                                           // Store the line number corresponding to the sequence ID instead
                }else{
                    rand_lines_list.push_back(line_number);                                             // If the line number is odd, it is already corresponding to the sequence ID postiontion in fasta file
                }
                reads_to_generate--;
            }
            std::sort(rand_lines_list.begin(), rand_lines_list.end());                                  // This is a sorted list of randomly obtained line numbers for the cell fasta file. To use to generate reads
            
            int line_count = 1;                                                                         // Temporary vector to count the number of lines read from the fasta files using the getNextChromSeq function
            int reads_actually_generated = 0;                                                           // Temporary counter variable to count the number of reads actually generated for each cell; for the report
            while(getNextChromSeq(cell_fasta_file, chromSegmentSeq, chromSegmentSeq_ID)){               // Iterate through each chromsome segment sequence in the fasta file of the cell that is being sequenced
                if(rand_lines_list.empty()) break;                                                      // If no more reads to be generated from the current cell, then break the loop
                
                if(!ART::init_set(parameter.get_read_length(), chromSegmentSeq)){                       // If chromSegmentSeq is smaller than the read length, then continue with the next segment
                    continue;
                };
                
                std::vector<int>::iterator it;                                                          // Temporary iterator object for the rand_lines_list vector
                while ((it = std::find(rand_lines_list.begin(), rand_lines_list.end(), line_count)) != rand_lines_list.end()) {
                    read1.generate_read_with_indel();                                                   // Make a read with random indel errors
                    read1_quality_score_vec.clear();
                    read1.get_read_quality(read1_quality_score_vec, 1);                                 // Get the read quality scores for the read positions
                    read1.add_baseCall_error(read1_quality_score_vec);                                  // Add base call errors to the read based on the quality scores
                    
                    fastq_R1_file<<"@"<<chromSegmentSeq_ID<<"_read"<<std::to_string(line_count)<<"\n"   // @readID
                                 <<(*read1.get_final_read_sequence())                                   // read sequence
                                 <<"\n+\n";                                                             // +
                    for(int k=0; k<read1_quality_score_vec.size(); k++){                                // read quality scores
					    fastq_R1_file<<static_cast<char>(read1_quality_score_vec[k]+33);
				    }
                    fastq_R1_file<<"\n";
                    
                    rand_lines_list.erase(it);                                                          // Remove the line from the list after making a read
                    total_num_reads--;                                                                  // When a read is generated, decrement the count
                    reads_actually_generated++;
                }
                line_count += 2;                                                                        // Since two lines (ID and sequence) were read from the fasta file, increment line count by 2
            } 
            cell_fasta_file.close();

            report_cells_sequenced.push_back(cellGenomes_to_be_sequenced[cell_index]);                  // Storing the names of the cells that were actually sequenced, for the report
            report_readsGenerated_perCell.push_back(reads_actually_generated);                          // Storing the number of reads generated per cell in a vector for the final summary report

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

        while(total_num_reads>0){
            int cell_index = rng::rand_int(1,cellGenomes_to_be_sequenced.size()) - 1;                   // Randomly choose a cell to generate reads. 'cell' can take values from 0 to cellGenomes_to_be_sequenced.size minus 1
            std::string cell_fasta_filename = (*parameter.get_output_directory())+"/temp/"+cellGenomes_to_be_sequenced[cell_index];

            int reads_to_generate = rng::rand_int(0,max_num_reads_per_cell);                            // Randomly determine how many reads should be generated from this particular cell that is selected randomly. This is to ensure uniformity
            if(max_num_reads_per_cell < 2){                                                             // If max_num_reads_per_cell is zero, this will prevent program getting stuck
                reads_to_generate = rng::rand_int(0,total_num_reads);
            }
            while(reads_to_generate>total_num_reads){                                                   // If the number of reads to generate exceeds the number of reads remaining to make in the current loop, then regenerate the number
                reads_to_generate = rng::rand_int(0,total_num_reads);
            }
            
            std::ifstream cell_fasta_file(cell_fasta_filename);                                         // ifstream object of the fasta file to be sequenced
            std::string chromSegmentSeq;                                                                // Temporary variable to hold each chromosome segment sequence from the fasta file one at a time
            std::string chromSegmentSeq_ID;                                                             // Temporary variable to hold each chromsome segment ID
            
            int totalLines_in_fasta = 0;                                                                // Variable to hold the number of lines in the file file of the cell selected
            while (std::getline(cell_fasta_file, chromSegmentSeq)) {                                    // Counting the number of lines in the fasta file for the cell
                totalLines_in_fasta++;
            }
            cell_fasta_file.clear();                                                                    // Reset the error flags of the stream so that subsequent operations can be performed without any issues
            cell_fasta_file.seekg(0);                                                                   // Seeking back to the beginning of the file, ensuring that the stream is in a valid state

            std::vector<int> rand_lines_list;                                                           // Vector to hold the fasta file lines that we want to use to generate reads at random
            while (reads_to_generate > 1){                                                              // reads_to_generate > 1 is used to avoid over generating reads per cell, since it's paired-end
                int line_number = rng::rand_int(1,totalLines_in_fasta);                                 // line_number can take any value from 0 to length of the fasta file
                if (line_number %2 == 0){                                                               // If the generated line number is even, it corresponds to the sequence line, not sequence ID line
                    rand_lines_list.push_back(line_number-1);                                           // Store the line number corresponding to the sequence ID instead
                }else{
                    rand_lines_list.push_back(line_number);                                             // If the line number is odd, it is already corresponding to the sequence ID postiontion in fasta file
                }
                reads_to_generate -= 2;                                                                 // Paired-end will create two reads from one fasta file line (chrom segment sequence)
            }
            std::sort(rand_lines_list.begin(), rand_lines_list.end());                                  // This is a sorted list of randomly obtained line numbers for the cell fasta file. To use to generate reads
            
            int line_count = 1;                                                                         // Temporary vector to count the number of lines read from the fasta files using the getNextChromSeq function
            int reads_actually_generated = 0;                                                           // Temporary counter variable to count the number of reads actually generated for each cell; for the report
            while(getNextChromSeq(cell_fasta_file, chromSegmentSeq, chromSegmentSeq_ID)){               // Iterate through each chromsome segment sequence in the fasta file of the cell that is being sequenced
                if(rand_lines_list.empty()) break;                                                      // If no more reads to be generated from the current cell, then break the loop
                
                if(!ART::init_set(parameter.get_read_length(), chromSegmentSeq)){                       // If chromSegmentSeq is smaller than the read length, then continue with the next segment
                    continue;
                };
                
                std::vector<int>::iterator it;                                                          // Temporary iterator object for the rand_lines_list vector
                while ((it = std::find(rand_lines_list.begin(), rand_lines_list.end(), line_count)) != rand_lines_list.end()) {
                    ART::generate_paired_reads_with_indel(read1, read2, parameter.get_mean_DNA_fragment_length(), parameter.get_std_dev_DNA_fragment_length());
                                                                                                        // Make two paired-reads from the same DNA fragment with indel errors
                    // Process read 1 first
                    read1_quality_score_vec.clear();
                    read1.get_read_quality(read1_quality_score_vec, 1);                                 // Get the read quality scores for the read positions for read 1
                    read1.add_baseCall_error(read1_quality_score_vec);                                  // Add base call errors to the read based on the quality scores on read 2
                    
                    fastq_R1_file<<"@"<<chromSegmentSeq_ID<<"_read"<<std::to_string(line_count)<<"\n"   // @readID
                                 <<(*read1.get_final_read_sequence())                                   // read sequence
                                 <<"\n+\n";                                                             // +
                    for(int k=0; k<read1_quality_score_vec.size(); k++){                                // read quality scores
					    fastq_R1_file<<static_cast<char>(read1_quality_score_vec[k]+33);
				    }
                    fastq_R1_file<<"\n";

                    // Process read 2
                    read2_quality_score_vec.clear();
                    read2.get_read_quality(read2_quality_score_vec, 2);                                 // Get the read quality scores for the read positions for read 2
                    read2.add_baseCall_error(read2_quality_score_vec);                                  // Add base call errors to the read based on the quality scores on read 2
                    
                    fastq_R2_file<<"@"<<chromSegmentSeq_ID<<"_read"<<std::to_string(line_count)<<"\n"   // @readID
                                 <<(*read2.get_final_read_sequence())                                   // read sequence
                                 <<"\n+\n";                                                             // +
                    for(int k=0; k<read2_quality_score_vec.size(); k++){                                // read quality scores
					    fastq_R2_file<<static_cast<char>(read2_quality_score_vec[k]+33);
				    }
                    fastq_R2_file<<"\n";
                    
                    rand_lines_list.erase(it);                                                          // Remove the line from the list after making a read
                    total_num_reads -= 2;                                                               // When paired reads (2) are generated, decrement the count
                    reads_actually_generated += 2;
                }
                line_count += 2;                                                                        // Since two lines (ID and sequence) were read from the fasta file, increment line count by 2
            } 
            cell_fasta_file.close();

            report_cells_sequenced.push_back(cellGenomes_to_be_sequenced[cell_index]);                  // Storing the names of the cells that were actually sequenced, for the report
            report_readsGenerated_perCell.push_back(reads_actually_generated);                          // Storing the number of reads generated per cell in a vector for the final summary report
        }
        fastq_R2_file.close();
    }
    fastq_R1_file.close();

    report_fastq_output.push_back(*parameter.get_output_fastq_filename_prefix());                       // Output filename without the extension for the report
}
//--------------------------------------------------------------------------------------------