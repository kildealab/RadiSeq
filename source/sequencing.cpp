#include "sequencing.h"
#include "fileio.h"
#include "art_framework.h"
#include "random_generator.h"
#include "summary_report.h"
#include "support_functions.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <sys/mman.h>
#include <cmath>
#include <omp.h>

//--------------------------------------------------------------------------------------------
// This is the function that will perform the single-cell sequencing of all the cells that are
// sampled from the pool of cells one-by-one. This function generates reads according to the
// coverage per cell for every cell. The generated reads will have indel errors as well as
// base call errors based on the sequencer's error quality profile. The read data will be 
// stored to fastq.gz files that correspond to each cell.
//--------------------------------------------------------------------------------------------
void single_cell_sequencing(NGSParameters& parameter, const std::vector<std::string>& cellGenomes_to_be_sequenced){
    double coverage_per_cell = parameter.get_total_read_coverage()*0.5/parameter.get_num_of_cells_to_sequence(); // Halved because reads can come from both strands
    ART::initiate_read_generation(parameter.get_read_length(), parameter.get_GC_binSize(), parameter.get_fraction_nonFR_read_pairs(), parameter.get_read_artifacts_rate());
    ART::set_read_quality_distribution(*parameter.get_r1_quality_profile(), *parameter.get_r2_quality_profile());
    //ART::set_baseCall_error_probability();                                                              // Generate base call error probabilty distribution (same for all cells)
    std::vector<double> fragment_weights;                                                                 // Vector to hold the normalized fragment weights
    int min_DNA_fragment_length{1};                                                                       // Variable to hold the minimum DNA fragment length to be generated
    if(parameter.get_is_fragment_distribution_from_file()){                                               // If the fragment size distribution needs to be created from the file provided for paired-end sequencing
        auto fragmentData = readFragmentSizeDist(parameter.get_fragment_size_distribution_path());        // Read the fragment size distribution file from the path
        min_DNA_fragment_length = fragmentData.first;                                                     // First in the returned value pair is the minimum DNA fragment length
        fragment_weights = fragmentData.second;                                                           // Second is the fragment weights
    }else{                                                                                                // Generate the fragmentation probability vector for paired-end sequencing using a beta function if that is needed
        fragment_weights = beta_distribution_proabalities(parameter.get_beta_of_beta_distribution(), parameter.get_min_DNA_fragment_length(), parameter.get_max_DNA_fragment_length(),parameter.get_mode_DNA_fragment_length());
        min_DNA_fragment_length = parameter.get_min_DNA_fragment_length();
    }
    
    for (size_t i=0; i<cellGenomes_to_be_sequenced.size(); i++){                                        // Iterate though each cell genome that is to be sequenced to generate reads
        std::string cell_fasta_filename = (*parameter.get_output_directory())+"/temp/"+cellGenomes_to_be_sequenced[i];
        if (checkFileExists(&cell_fasta_filename)){;                                                    // Check if the fasta file of the cell to be sequenced can be read
            std::cout<<"\n Sequencing of cell "<<i+1<<" is in progress \n";
        }else{
            std::cerr<<"\n WARNING: Unable to sequence cell "<<i+1<<". This cell will be ignored.\n";
            continue;
        }

        size_t position{0};                                                                             // Temporary variable to hold the last-read position in the memory map of the current fasta file
        size_t fastaFileSize{0};                                                                        // A variable to hold the size of the fasta file of the current cell, during memory-mapping
        void* fastaFileMM = generateInputFileMemoryMap(cell_fasta_filename, fastaFileSize);             // Create the memory-map of the fasta file
        const char* cell_fastaFileData = static_cast<char*>(fastaFileMM);                               // Casting the memory-map void pointer to a const char pointer for further processing

        std::string output_fastq_R1_filename = (*parameter.get_output_directory())+"/"+(*parameter.get_output_fastq_filename_prefix())+"_"+std::to_string(i)+"_R1.fastq.gz";
        std::ofstream fastq_R1_file(output_fastq_R1_filename.c_str(),std::ios::binary);                 // ofstream object of the output fastq file for read 1
        
        ART read1;                                                                                      // Creating an ART class object and setting the insertion and deletion probability vectors for that read object
        read1.set_read_error_rates(parameter.get_insertion_error_rate_read1(), parameter.get_deletion_error_rate_read1());
        //read1.set_read_error_probability(parameter.get_read_length(), parameter.get_insertion_error_rate_read1(), read1.insertion_probability_vec, parameter.get_max_errors_in_read());
        //read1.set_read_error_probability(parameter.get_read_length(), parameter.get_deletion_error_rate_read1(), read1.deletion_probability_vec, parameter.get_max_errors_in_read());

        std::string chromSegSeq;                                                                        // Temporary variable to hold each chromosome segment sequence from the fasta file one at a time
        std::string chromSegSeq_ID;                                                                     // Temporary variable to hold IDs of each hromosome segment sequence
        int end_flag;                                                                                   // Flag to indicate the end of memory map sequences
        bool goodSegment{true};                                                                         // Flag to indicate if the chromSegSeq is good for further processing

        const int batchSize{2000};                                                                      // Define a batch size for writing reads to the output file. These much data will be stored in cache before writing it on the file
        
        int nThreads_User = parameter.get_number_of_threads();                                          // Variable holding the number of threads the user requesting for parallel processing
        int batchSize_thread = std::round(batchSize/nThreads_User);                                     // Devide the total cache size for the buffer equally for all the threads

        // Single-end sequencing
        if(!parameter.get_paired_end_sequencing()){                                                     // If user asked for single-end (not paired-end) sequencing
            std::vector<std::vector<std::string>> batch_buffer(nThreads_User);                          // Create a buffer vector with vectors for each thread to store read data
            #pragma omp parallel shared(chromSegSeq,chromSegSeq_ID,end_flag,goodSegment)                // Start the parallel region. Threads will be generated and they will get the specific shared variables
            {
                int local_end_flag = 1;                                                                 // A flag for each thread to check if the end of memory map is reached. 
                while(local_end_flag){
                    #pragma omp master                                                                  // Only one thread should execute this block; other threads will skip this section
                    {
                        end_flag = getNextChromSeq_MM(cell_fastaFileData, fastaFileSize, position, chromSegSeq, chromSegSeq_ID); // 0 if end of the memory map is reached
                    }
                    #pragma omp barrier                                                                 // Threads should wait here till all the threads reach this point
                    #pragma omp flush(end_flag)                                                         // Make sure the end_flag variable gets the latest updated value and not get cached
                    local_end_flag = end_flag;                                                          // Update the local flag so that all the threads get the updated value
                    long num_reads_per_segment = static_cast<long>((coverage_per_cell*chromSegSeq.size())/parameter.get_read_length());
                
                    if((chromSegSeq.size()-parameter.get_read_length()) > 0){                           // Proceed only if chromSegmentSeq is bigger than the read length, otherwise continue with the next segment
                        #pragma omp master                                                              // This section needs to be done by only one thread
                        {
                            int nThreads_omp = omp_get_num_threads();                                   // Get the actual number of threads available for OpenMP
                            goodSegment = true;                                                         // Reset flag for each segment
                            goodSegment = goodSegment && read1.int_set(chromSegSeq, nThreads_omp);      // Boolean AND operation ensures that TRUE only when both conditons are TRUE
                            //goodSegment = goodSegment && ART::get_N_mask_and_GCbias(parameter.get_N_threshold_in_reads());// Mask regions in the chromSegSeq where the number of N's exceed the threshold to avoid making reads
                            goodSegment = goodSegment && ART::GCbias_maker();                           // Calculate the GC bias for the chromSegSeq 
                            
                            if(*parameter.get_coverage_distribution() == "mda"){                        // Prepare the read distribution if the selected coverage distributon is MDA
                                goodSegment = goodSegment && ART::MDA_distribution_maker(num_reads_per_segment, coverage_per_cell);
                            }
                        }
                        #pragma omp barrier                                                             // Make sure all threads reach this point, before proceeding with the rest
                        #pragma omp flush(goodSegment)                                                  // Make sure the goodSegment variable gets the latest updated value and not get cached
                        bool local_goodSegment = goodSegment;                                           // Update the local flag so that all the threads get the updated value
                        
                        if(local_goodSegment){                                                              // Proceed only if the chromoSegSeq is good
                            #pragma omp for schedule(dynamic)                                               // This is where we are splitting the iterations of the for loop to each thread
                            for(int j=0; j<num_reads_per_segment; j++){
                                int threadID = omp_get_thread_num();                                        // Get the ID of each thread being tracked
                                read1.generate_read_with_indel(threadID);                                   // Make a read with random indel errors
                                std::vector<short> read1_quality_score_vec;                                 // Vector to hold the quality scores for read 1
                                read1.get_read_quality(read1_quality_score_vec, 1, threadID);               // Get the read quality scores for the read positions
                                read1.add_baseCall_error(read1_quality_score_vec, threadID);                // Add base call errors to the read based on the quality scores
                                
                                std::string chromID = chromSegSeq_ID; chromID.erase(0,1);                   // Remove the '>' symbol from the chrom ID
                                std::string read_data = "@"+chromID+"_read"+std::to_string(j)+"\n";         // @readID
                                read_data += (*read1.get_final_read_sequence(threadID))+ "\n+\n";           // read sequence and +
                                for(size_t k=0; k<(*read1.get_final_read_sequence(threadID)).size(); k++){  // read quality scores; insert only as many quality values as with the length of sequence
                                    read_data += static_cast<char>(read1_quality_score_vec[k]+32);          // +33 to get the phred score
                                }
                                read_data += "\n";
                                
                                batch_buffer[threadID].push_back(read_data);                                // Add the read data to the buffer vector of the respective thread                    
                                if (batch_buffer[threadID].size() >= static_cast<size_t>(batchSize_thread)){// Check if the thread's batch buffer is full, and write it to the file if needed.
                                    #pragma omp critical(section1)
                                    {
                                        writeBatchToFile(batch_buffer[threadID], fastq_R1_file, true);
                                    }
                                }
                            }
                        }
                    }
                    #pragma omp barrier
                }
            }
            for (size_t l=0;l<batch_buffer.size();l++){
                writeBatchToFile(batch_buffer[l], fastq_R1_file, true);                                 // If there are unwritten data in batch buffer, write that too when the loop ends
            }
        }
        // Paired-end sequencing
        else{                                                                                           // If user asked to perform paired-end sequencing
            std::string output_fastq_R2_filename = (*parameter.get_output_directory())+"/"+(*parameter.get_output_fastq_filename_prefix())+"_"+std::to_string(i)+"_R2.fastq.gz";
            std::ofstream fastq_R2_file(output_fastq_R2_filename.c_str(),std::ios::binary);             // ofstream object of the output fastq file for read 2
            
            ART read2;                                                                                  // Creating an ART class object and setting the insertion and deletion probability vectors for that read object
            read2.set_read_error_rates(parameter.get_insertion_error_rate_read2(),parameter.get_deletion_error_rate_read2());
            //read2.set_read_error_probability(parameter.get_read_length(), parameter.get_insertion_error_rate_read2(), read2.insertion_probability_vec, parameter.get_max_errors_in_read());
            //read2.set_read_error_probability(parameter.get_read_length(), parameter.get_deletion_error_rate_read2(), read2.deletion_probability_vec, parameter.get_max_errors_in_read());

            int good_iterations{0};                                                                     // Temporary counter variable to count the number of iterations where a read was generated. Read will not be generated if some condition is not met
    
            std::vector<std::vector<std::string>> batch_buffer_r1(nThreads_User);                       // Create a buffer for each thread for storing read 1 data.
            std::vector<std::vector<std::string>> batch_buffer_r2(nThreads_User);                       // Create a buffer for each thread for storing read 2 data.
            
            #pragma omp parallel shared(chromSegSeq,chromSegSeq_ID,end_flag,goodSegment,good_iterations)// Start the parallel region. Threads will be generated and they will get the specific shared variables
            {
                int local_end_flag = 1;                                                                 // A flag for each thread to check if the end of memory map is reached. 
                while(local_end_flag){
                    #pragma omp master                                                                  // Only the master thread should execute this block; other threads will skip this section
                    {   
                        end_flag = getNextChromSeq_MM(cell_fastaFileData, fastaFileSize, position, chromSegSeq, chromSegSeq_ID); // 0 if end of the memory map is reached
                    }
                    #pragma omp barrier                                                                 // Threads should wait here till all the threads reach this point
                    #pragma omp flush(end_flag)                                                         // Make sure the end_flag variable gets the latest updated value and not get cached
                    local_end_flag = end_flag;                                                          // Update the local flag so that all the threads get the updated value
                    long num_reads_per_segment = static_cast<long>((coverage_per_cell*chromSegSeq.size())/parameter.get_read_length());
                    
                    if(num_reads_per_segment>=1 && (static_cast<int>(chromSegSeq.size())>parameter.get_read_length())){   // Proceed only if at least one read is needed from chromSegmentSeq and segment is bigger than the read length, otherwise continue with the next segment
                        #pragma omp master                                                              // This section needs to be done by only one thread
                        {   
                            int nThreads_omp = omp_get_num_threads();                                   // Get the actual number of threads available for OpenMP
                            goodSegment = true;                                                         // Reset flag for each segment
                            goodSegment = goodSegment && read1.int_set(chromSegSeq, nThreads_omp);
                            goodSegment = goodSegment && read2.int_set(chromSegSeq, nThreads_omp);
                            //goodSegment = goodSegment && ART::get_N_mask_and_GCbias(parameter.get_N_threshold_in_reads(), parameter.get_min_DNA_fragment_length()); // Mask regions in the chromSegSeq where the number of N's exceed the threshold to avoid making reads
                            goodSegment = goodSegment && ART::GCbias_maker(min_DNA_fragment_length);
                            if(*parameter.get_coverage_distribution() == "mda" && goodSegment){         // Divided the number of reads per segment by 2 coz each site generate two reads in paired-end seq
                                goodSegment = goodSegment && ART::MDA_distribution_maker(num_reads_per_segment/2, coverage_per_cell);
                            }
                            good_iterations = 0;                                                        // Reset the counter for each chrom segment
                        }
                        #pragma omp barrier                                                             // Make sure all threads reach this point, before proceeding with the rest
                        #pragma omp flush(goodSegment)                                                  // Make sure the goodSegment variable gets the latest updated value and not get cached
                        bool local_goodSegment = goodSegment;                                           // Update the local flag so that all the threads get the updated value
                        
                        if(local_goodSegment){                                                          // Proceed only if the chromoSegSeq is good
                            while(good_iterations<num_reads_per_segment){
                                int target_reads = num_reads_per_segment - good_iterations;             // Remaining reads to be generated in each while loop
                                int local_good_iterations{0};                                           // Thread-local counter for iterations that generated reads
                                #pragma omp for schedule(dynamic) 
                                for(int j=0; j<target_reads; j+=2){                                     // Create as many reads per segment in a loop. J is incremented by 2 since two reads are created in one loop (paired end)
                                    int threadID = omp_get_thread_num();                                               // Get the ID of each thread being tracked
                                    if(ART::generate_paired_reads_with_indel(read1, read2, min_DNA_fragment_length, fragment_weights, threadID)){
                                                                                                                       // Make two paired-reads from the same DNA fragment with indel errors
                                        // Process read 1 first
                                        std::vector<short> read1_quality_score_vec;                                    // Vector to hold the quality scores for read 1
                                        read1.get_read_quality(read1_quality_score_vec, 1,threadID);                   // Get the read quality scores for the read positions for read 1
                                        read1.add_baseCall_error(read1_quality_score_vec,threadID);                    // Add base call errors to the read based on the quality scores on read 2
                                        
                                        std::string chromID = chromSegSeq_ID; chromID.erase(0,1);                      // Remove the '>' symbol from the chrom ID
                                        std::string read1_data = "@"+chromID+"_read"+std::to_string(j)+"/1\n";         // @readID
                                        read1_data += (*read1.get_final_read_sequence(threadID)) + "\n+\n";            // read sequence and +
                                        for(size_t k=0; k<(*read1.get_final_read_sequence(threadID)).size(); k++){     // read quality scores
                                            read1_data += static_cast<char>(read1_quality_score_vec[k]+32);            // +33 to get the phred score
                                        }
                                        read1_data += "\n";

                                        // Process read 2
                                        std::vector<short> read2_quality_score_vec;                                    // Vector to hold the quality scores for read 2
                                        read2.get_read_quality(read2_quality_score_vec, 2,threadID);                   // Get the read quality scores for the read positions for read 2
                                        read2.add_baseCall_error(read2_quality_score_vec,threadID);                    // Add base call errors to the read based on the quality scores on read 2
                                        
                                        std::string read2_data = "@"+chromID+"_read"+std::to_string(j)+"/2\n";         // @readID
                                        read2_data += (*read2.get_final_read_sequence(threadID)) + "\n+\n";            // read sequence and +
                                        for(size_t k=0; k<(*read2.get_final_read_sequence(threadID)).size(); k++){     // read quality scores
                                            read2_data += static_cast<char>(read2_quality_score_vec[k]+32);            // +33 to get the phred score
                                        }
                                        read2_data += "\n";
                                        
                                        local_good_iterations+=2;                                                      // Increment the thread-local counter

                                        batch_buffer_r1[threadID].push_back(read1_data);                               // Add the read 1 data to the thread's buffer 1.   
                                        batch_buffer_r2[threadID].push_back(read2_data);                               // Add the read 2 data to the thread's buffer 2.                 
                                        if (batch_buffer_r1[threadID].size() >= static_cast<size_t>(batchSize_thread)){// Check if the batch buffer is full, and write it to the file if needed.
                                            #pragma omp critical(section1)
                                            {
                                                writeBatchToFile(batch_buffer_r1[threadID], fastq_R1_file, true);
                                                writeBatchToFile(batch_buffer_r2[threadID], fastq_R2_file, true);
                                            }
                                        }
                                    }
                                }
                                #pragma omp atomic                                                                     // Ensure only one thread modifies the shared variable at a time
                                good_iterations += local_good_iterations;
                                #pragma omp barrier   
                            }
                            #pragma omp barrier
                        }
                    }
                    #pragma omp barrier
                }
            }
            for (size_t l=0;l<batch_buffer_r1.size();l++){
                writeBatchToFile(batch_buffer_r1[l], fastq_R1_file, true);                              // If there are unwritten data in batch buffer 1, write that too when the loop ends
                writeBatchToFile(batch_buffer_r2[l], fastq_R2_file, true);                              // If there are unwritten data in batch buffer 1, write that too when the loop ends
            }
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
void bulk_cell_sequencing(NGSParameters& parameter, const std::vector<std::string>& cellGenomes_to_be_sequenced, const std::vector<std::vector<double>>& line_weight_in_cells_to_sequence, long ref_seq_length ){ 
    int total_num_reads = std::round((parameter.get_total_read_coverage()*(ref_seq_length))/parameter.get_read_length());
    int max_num_reads_per_cell = std::round(total_num_reads/parameter.get_num_of_cells_to_sequence());  // Temporary variable to hold the maximum number of reads per cell
    
    ART::initiate_read_generation(parameter.get_read_length(), parameter.get_GC_binSize(), parameter.get_fraction_nonFR_read_pairs(), parameter.get_read_artifacts_rate());
    ART::set_read_quality_distribution(*parameter.get_r1_quality_profile(), *parameter.get_r2_quality_profile());
    //ART::set_baseCall_error_probability();                                                              // Generate base call error probabilty distribution (same for all cells)
    
    ART read1;                                                                                          // Creating an ART class object and setting the insertion and deletion probability vectors for that read object
    read1.set_read_error_rates(parameter.get_insertion_error_rate_read1(), parameter.get_deletion_error_rate_read1());
    //read1.set_read_error_probability(parameter.get_read_length(), parameter.get_insertion_error_rate_read1(), read1.insertion_probability_vec, parameter.get_max_errors_in_read());
    //read1.set_read_error_probability(parameter.get_read_length(), parameter.get_deletion_error_rate_read1(), read1.deletion_probability_vec, parameter.get_max_errors_in_read());
    
    std::string output_fastq_R1_filename = (*parameter.get_output_directory())+"/"+(*parameter.get_output_fastq_filename_prefix())+"_R1.fastq.gz";
    std::ofstream fastq_R1_file(output_fastq_R1_filename.c_str(),std::ios::binary);                     // ofstream object of the output fastq file for read 1

    // Generate a list of lines that we want to read randomly from all the cell file prior to reading each file. This is done to avoid processing the same fasta files multiple times to generate reads from it one by one. 
    // Here we generate a map, that indicates which all lines from one fasta files should be read so that we call read all these lines in one go. 
    int readPerFrag{1};                                                                                 // A variable to indicate how many reads per fragments we need to create. It will get 1 if single-end sequencing
    if(parameter.get_paired_end_sequencing()){
        readPerFrag = 2;                                                                                // This variable gets 2 if we want paired-end sequencing
    }
    std::vector<std::string> cellGenomes;                                                               // A vector to hold the names of cells that are randomly selected to sequence
    std::vector<std::unordered_map<int, int>> line_frequency_map_vec;                                   // A vector of maps to hold the randoms lines and its frequencies (ie, <linenumber, frequency>) in each cell files that we want to read to generate the reads, in order of cell in cellGenome vector
    while (total_num_reads > (readPerFrag-1)){
        int cell_index = rng::rand_int(1,static_cast<int>(cellGenomes_to_be_sequenced.size())) - 1;     // Randomly choose a cell to generate reads. 'cell_index' can take values from 0 to cellGenomes_to_be_sequenced.size minus 1
        int reads_to_generate_in_cell = rng::rand_int(0,max_num_reads_per_cell);                        // Randomly determine how many reads should be generated from this particular cell that is selected randomly
        
        if(max_num_reads_per_cell < readPerFrag ){                                                      // If max_num_reads_per_cell is zero, this will prevent program getting stuck
            reads_to_generate_in_cell = rng::rand_int(0,total_num_reads);
        } 
        while(reads_to_generate_in_cell>total_num_reads){                                               // If the number of reads to generate exceeds the number of reads remaining to make in the current loop, then regenerate the number
           reads_to_generate_in_cell = rng::rand_int(0,total_num_reads);
        }               
       
        int reads_will_be_generated_in_cell{0};                                                         // A temporary variable to hold the number of reads that will actually be created from each cell. For single-end seq this number will be equal to reads_to_generate_in_cell
        while (reads_to_generate_in_cell > (readPerFrag-1)){                                            // Randomly generate a list of line numbers that we want to use pick the chrom fragments from the fasta file
            int line_number = rng::weighted_rand_int(line_weight_in_cells_to_sequence[cell_index]);     // line_number can take any value from 1 to length of the fasta file
            if (line_number==0){                                                                        // weighted_rand_int returns 0 if the weight vector of the cell is inappropriate
                break;                                                                                  // Skip the cell if its weight vector is bad
            }
            if(line_number %2 == 0){                                                                    // If the generated line number is even, it corresponds to the sequence line, not sequence ID line
                line_number -= 1;                                                                       // Store the line number corresponding to the sequence ID instead. If the line number is odd, it is already corresponding to the sequence ID postiontion in fasta file
            }
            auto it = std::find(cellGenomes.begin(), cellGenomes.end(), cellGenomes_to_be_sequenced[cell_index]); 
            if(it != cellGenomes.end()){                                                                // Check vector if this cell was previously selected randomly and in the list already
                int index = std::distance(cellGenomes.begin(), it);                                     // If the cell already exists in the vector, get the index in the vector
                auto& line_freq_map = line_frequency_map_vec[index];                                    // Use the index to access the corresponding map
                if (line_freq_map.find(line_number) != line_freq_map.end()){                            // Check if the line is already an entry in the map
                    line_freq_map[line_number] += 1;                                                    // If the line already exists in the frequency map, then increment the frequency
                }else{                                                                                  // If the line is not in the map already, initialize with a frequency of 1
                    line_freq_map[line_number] = 1;   
                }
            }else{                                                                                      // If the cell doesn't exist in the list already, create a new entry in both vectors
                cellGenomes.push_back(cellGenomes_to_be_sequenced[cell_index]);                         // Add cell name to the list
                std::unordered_map<int, int> newMap;                                                    // Make a new map for the cell
                newMap[line_number] = 1;
                line_frequency_map_vec.push_back(newMap);
            }
            reads_to_generate_in_cell -= readPerFrag;
            reads_will_be_generated_in_cell += readPerFrag;
        }
        total_num_reads -= reads_will_be_generated_in_cell;
    } 
    const int batchSize{2000};                                                                          // Define a batch size for writing reads to the output file. These much data will be stored in cache before writing it on the file

    int nThreads_User = parameter.get_number_of_threads();                                              // Variable holding the number of threads the user requesting for parallel processing
    int batchSize_thread = std::round(batchSize/nThreads_User);                                         // Devide the total cache size for the buffer equally for all the threads

    // Single-end sequencing
    if(!parameter.get_paired_end_sequencing()){
        std::vector<std::vector<std::string>> batch_buffer(nThreads_User);                              // Create a vector for each thread buffer for storing read data.
        
        for (size_t cell=0; cell<cellGenomes.size(); cell++){                                           // Iterate through each cell that we want to sequence to do the sequencing
            std::string cell_fasta_filename = (*parameter.get_output_directory())+"/temp/"+cellGenomes[cell];
            size_t position{0};                                                                         // Temporary variable to hold the last-read position in the memory map of the current fasta file
            size_t fastaFileSize{0};                                                                    // A variable to hold the size of the fasta file of the current cell, during memory-mapping
            void* fastaFileMM = generateInputFileMemoryMap(cell_fasta_filename, fastaFileSize);         // Create the memory-map of the fasta file
            const char* cell_fastaFileData = static_cast<char*>(fastaFileMM);                           // Casting the memory-map void pointer to a const char pointer for further processing

            int line_inFile{1};                                                                         // Temporary variable to count the number of lines read from the fasta memory map
            int reads_actually_generated{0};                                                            // Temporary counter variable to count the number of reads actually generated for each cell; for the report
            
            std::unordered_map<int, int> line_frequency_map = line_frequency_map_vec[cell];             // Get the line frequency map of the cell being processed
            size_t entries_in_freqMap = line_frequency_map.size();                                      // Get the number of entries in the line frequency map for later use
            
            std::string chromSegSeq;                                                                    // Temporary variable to hold each chromosome segment sequence from the fasta file one at a time
            std::string chromSegSeq_ID;                                                                 // Temporary variable to hold IDs of each hromosome segment sequence
            int end_flag;                                                                               // Flag to indicate the end of memory map sequences
            bool goodSegment{true};                                                                     // Flag to indicate if the chromSegSeq is good for further processing

            #pragma omp parallel shared(chromSegSeq,chromSegSeq_ID,end_flag,line_inFile,goodSegment)    // Start the parallel region. Threads will be generated and they will get the specific shared variables
            {
                int line_num = line_inFile;                                                             // A counter for each thread to check which chromosome segment is being processed from the memory map
                int local_end_flag = 1;                                                                 // A flag for each thread to check if the end of memory map is reached. 
                while(local_end_flag){
                    #pragma omp master                                                                  // Only one thread should execute this block; other threads will skip this section
                    {
                        end_flag = getNextChromSeq_MM(cell_fastaFileData, fastaFileSize, position, chromSegSeq, chromSegSeq_ID); // 0 if end of the memory map is reached
                        if(entries_in_freqMap == 1){end_flag = 0;}                                      // If there is only 1 entry left in the map, set the flag so that loop exits after processing it (there will be atleast one entry per cell). 0 --> no more reads to be generated from this sequence
                    }
                    #pragma omp barrier                                                                 // Threads should wait here till all the threads reach this point
                    #pragma omp flush(end_flag)                                                         // Make sure the end_flag variable gets the latest updated value and not get cached
                    local_end_flag = end_flag;                                                          // Update the local flag so that all the threads get the updated value
                    
                    auto iter = line_frequency_map.find(line_num);                                      // Check if the current sequence seqment line is amoung the list of lines that we want to sequence
                    if (iter != line_frequency_map.end()){                                              // If the line is in the list to be sequenced, perform the sequencing
                        if(static_cast<int>(chromSegSeq.size())>parameter.get_read_length()){           // Proceed only if chromSegmentSeq is bigger than the read length, otherwise continue with the next segment
                            #pragma omp master                                                          // This section needs to be done by only one thread
                            {
                                int nThreads_omp = omp_get_num_threads();                               // Get the actual number of threads available for OpenMP
                                goodSegment = true;                                                     // Reset flag for each segment
                                goodSegment = goodSegment && read1.int_set(chromSegSeq, nThreads_omp);  // TRUE only when both conditions are TRUE
                                //goodSegment = goodSegment && ART::get_N_mask_and_GCbias(parameter.get_N_threshold_in_reads());// Mask regions in the chromSegSeq where the number of N's exceed the threshold to avoid making reads
                                goodSegment = goodSegment && ART::GCbias_maker();                       // Determine the GC bias distribution in the chromSegSeq
                            }
                            #pragma omp barrier                                                         // Make sure all threads reach this point, before proceeding with the rest
                            #pragma omp flush(goodSegment)                                              // Make sure the goodSegment variable gets the latest updated value and not get cached
                            if(goodSegment){                                                            // Proceed only if the chromSegSeq is good
                                #pragma omp for schedule(dynamic) reduction (+:reads_actually_generated)// This is where we are splitting the iterations of the for loop to each thread with a reduction variable to specify that it is shared
                                for(int i=0; i<line_frequency_map[line_num]; i++){                      // line_frequency_map[line_num] gives the frequency of that line in the list, which is equal to the number of reads we want to generate from that line
                                    int threadID = omp_get_thread_num();                                // Get the ID of each thread being tracked
                                    read1.generate_read_with_indel(threadID);                           // Make a read with random indel (insertions and deletions) errors
                                    std::vector<short> read1_quality_score_vec;                         // Vector to hold the quality scores for read 1
                                    read1.get_read_quality(read1_quality_score_vec, 1, threadID);       // Get the read quality scores for the read positions
                                    read1.add_baseCall_error(read1_quality_score_vec, threadID);        // Add base call errors to the read based on the quality scores
                                    
                                    std::string cellName = cellGenomes[cell]; cellName.erase(cellName.length()-3);
                                    std::string chromID = chromSegSeq_ID; chromID.erase(0,1);           // Remove the '>' symbol from the chrom ID
                                    std::string read_data = "@"+cellName+"_"+chromID+"_read"+std::to_string(reads_actually_generated)+"\n";  // @readID
                                    read_data += (*read1.get_final_read_sequence(threadID))+ "\n+\n";                                        // read sequence and +
                                    for(size_t k=0; k<(*read1.get_final_read_sequence(threadID)).size(); k++){                               // read quality scores
                                        read_data +=static_cast<char>(read1_quality_score_vec[k]+32);                                        // +33 to get the phred score
                                    }
                                    read_data +="\n";
                                
                                    reads_actually_generated +=1;
                                    
                                    batch_buffer[threadID].push_back(read_data);                                // Add the read data to the batch buffer.                    
                                    if (batch_buffer[threadID].size() >= static_cast<size_t>(batchSize_thread)){// Check if the batch buffer is full, and write it to the file if needed.
                                        #pragma omp critical(section1)
                                        {
                                            writeBatchToFile(batch_buffer[threadID], fastq_R1_file, true);
                                        }
                                    }
                                }
                            }
                            
                        }
                        #pragma omp master                                                              // Only one thread should increment the counters
                        {
                            entries_in_freqMap--;                                                       // One entry in the map is processed, so decrement the number of entries remaining in the map
                            line_inFile +=2;                                                            // Increment the number of lines read from the fasta memory map                            
                        }    
                    }else{                                                                              // If the line currently got from the memory map need not be sequenced
                        #pragma omp single
                        {
                            line_inFile +=2;                                                            // Increment the number of lines read from the fasta memory map                            
                        }
                    }
                    #pragma omp barrier                                                                 // Wait till all the threads reach here
                    line_num = line_inFile;                                                             // Update the local counter so that all threads get the updated value
                }
            }
            munmap(fastaFileMM, fastaFileSize);                                                         // Unmap the fasta file of the current cell to avoid memory-leaks
    
            report_cells_sequenced.push_back(cellGenomes[cell]);                                        // Storing the names of the cells that were actually sequenced, for the report
            report_readsGenerated_perCell.push_back(reads_actually_generated);                          // Storing the number of reads generated per cell in a vector for the final summary report
        }
        for(size_t l=0; l<batch_buffer.size(); l++){
            writeBatchToFile(batch_buffer[l], fastq_R1_file, true);
        } 
    }
    // Paired-end sequencing
    else{
        ART read2;                                                                                      // Creating an ART class object and setting the insertion and deletion probability vectors for that read object
        read2.set_read_error_rates(parameter.get_insertion_error_rate_read2(), parameter.get_deletion_error_rate_read2());
        //read2.set_read_error_probability(parameter.get_read_length(), parameter.get_insertion_error_rate_read2(), read2.insertion_probability_vec, parameter.get_max_errors_in_read());
        //read2.set_read_error_probability(parameter.get_read_length(), parameter.get_deletion_error_rate_read2(), read2.deletion_probability_vec, parameter.get_max_errors_in_read());
        
        std::string output_fastq_R2_filename = (*parameter.get_output_directory())+"/"+(*parameter.get_output_fastq_filename_prefix())+"_R2.fastq.gz";
        std::ofstream fastq_R2_file(output_fastq_R2_filename.c_str(),std::ios::binary);                 // ofstream object of the output fastq file for read 2
        
        std::vector<double> fragment_weights;                                                           // Vector to hold the normalized fragment weights
        int min_DNA_fragment_length{1};                                                                 // Variable to hold the minimum DNA fragment length to be generated
        if(parameter.get_is_fragment_distribution_from_file()){                                         // If the fragment size distribution needs to be created from the file provided for paired-end sequencing
            auto fragmentData = readFragmentSizeDist(parameter.get_fragment_size_distribution_path());  // Read the fragment size distribution file from the path
            min_DNA_fragment_length = fragmentData.first;                                               // First in the returned value pair is the minimum DNA fragment length
            fragment_weights = fragmentData.second;                                                     // Second is the fragment weights
        }else{                                                                                          // Generate the fragmentation probability vector for paired-end sequencing using a beta function if that is needed
            min_DNA_fragment_length = parameter.get_min_DNA_fragment_length();
            fragment_weights = beta_distribution_proabalities(parameter.get_beta_of_beta_distribution(), min_DNA_fragment_length, parameter.get_max_DNA_fragment_length(),parameter.get_mode_DNA_fragment_length());
        }
        
        std::vector<std::vector<std::string>> batch_buffer_r1(nThreads_User);                           // Create a buffer for each thread for storing read 1 data.
        std::vector<std::vector<std::string>> batch_buffer_r2(nThreads_User);                           // Create a buffer for each thread for storing read 2 data.
        for (size_t cell=0; cell<cellGenomes.size(); cell++){                                           // Iterate through each cell that we want to sequence to do the sequencing
            std::string cell_fasta_filename = (*parameter.get_output_directory())+"/temp/"+cellGenomes[cell];
            size_t position{0};                                                                         // Temporary variable to hold the last-read position in the memory map of the current fasta file
            size_t fastaFileSize{0};                                                                    // A variable to hold the size of the fasta file of the current cell, during memory-mapping
            void* fastaFileMM = generateInputFileMemoryMap(cell_fasta_filename, fastaFileSize);         // Create the memory-map of the fasta file
            const char* cell_fastaFileData = static_cast<char*>(fastaFileMM);                           // Casting the memory-map void pointer to a const char pointer for further processing

            int line_inFile{1};                                                                         // Temporary variable to count the number of lines read from the fasta memory map
            int reads_actually_generated{0};                                                            // Temporary counter variable to count the number of reads actually generated for each cell; for the report
            int good_iterations{0};                                                                     // Temporary counter variable to count the number of iterations where a read was generated. Read will not be generated if some condition is not met

            std::unordered_map<int, int> line_frequency_map = line_frequency_map_vec[cell];             // Get the line frequency map of the cell being processed
            size_t entries_in_freqMap = line_frequency_map.size();                                      // Get the number of entries in the line frequency map for later use
            
            std::string chromSegSeq;                                                                    // Temporary variable to hold each chromosome segment sequence from the fasta file one at a time
            std::string chromSegSeq_ID;                                                                 // Temporary variable to hold IDs of each hromosome segment sequence
            int end_flag;                                                                               // Flag to indicate the end of memory map sequences
            bool goodSegment{true};                                                                     // Flag to indicate if the chromSegSeq is good for further processing
            #pragma omp parallel shared(chromSegSeq,chromSegSeq_ID,end_flag,line_inFile,goodSegment,good_iterations)// Start the parallel region. Threads will be generated and they will get the specific shared variables
            {
                int line_num = line_inFile;                                                             // A counter for each thread to check which chromosome segment is being processed from the memory map
                int local_end_flag = 1;                                                                 // A flag for each thread to check if the end of memory map is reached. 
                while(local_end_flag){
                    #pragma omp master                                                                  // Only one thread should execute this block; other threads will skip this section
                    {
                        end_flag = getNextChromSeq_MM(cell_fastaFileData, fastaFileSize, position, chromSegSeq, chromSegSeq_ID); // 0 if end of the memory map is reached
                        if(entries_in_freqMap == 1){end_flag = 0;}                                      // If there is only 1 entry left in the map, set the flag so that loop exits after processing it (there will be atleast one entry per cell). 0 --> no more reads to be generated from this sequence
                    }
                    #pragma omp barrier                                                                 // Threads should wait here till all the threads reach this point
                    #pragma omp flush(end_flag)                                                         // Make sure the end_flag variable gets the latest updated value and not get cached
                    local_end_flag = end_flag;                                                          // Update the local flag so that all the threads get the updated value
                    
                    auto iter = line_frequency_map.find(line_num);                                      // Check if the current sequence seqment line is amoung the list of lines that we want to sequence
                    if (iter != line_frequency_map.end()){                                              // If the line is in the list to be sequenced, perform the sequencing
                        if(static_cast<int>(chromSegSeq.size())>parameter.get_read_length()){           // Proceed only if chromSegmentSeq is bigger than the read length, otherwise continue with the next segment
                            #pragma omp master                                                          // This section needs to be done by only one thread
                            {
                                int nThreads_omp = omp_get_num_threads();                               // Get the actual number of threads available for OpenMP
                                goodSegment = true;                                                     // Reset flag for each segment
                                goodSegment = goodSegment && read1.int_set(chromSegSeq, nThreads_omp);
                                goodSegment = goodSegment && read2.int_set(chromSegSeq, nThreads_omp);
                                //goodSegment = goodSegment && ART::get_N_mask_and_GCbias(parameter.get_N_threshold_in_reads(), parameter.get_min_DNA_fragment_length());   // Mask regions in the chromSegSeq where the number of N's exceed the threshold to avoid making reads
                                goodSegment = goodSegment && ART::GCbias_maker(min_DNA_fragment_length);// Determine the GC bias distribution in the chromSegSeq
                                good_iterations = 0;                                                    // Reset the counter for each chrom segment
                            }
                            #pragma omp barrier                                                         // Make sure all threads reach this point, before proceeding with the rest
                            #pragma omp flush(goodSegment)                                              // Make sure the goodSegment variable gets the latest updated value and not get cached
                            if(goodSegment){
                                while(good_iterations<line_frequency_map[line_num]){                    // line_frequency_map[line_inFile] gives the frequency of that line in the list, which is equal to the number of reads we want to generate from that line
                                    int target_reads = line_frequency_map[line_num] - good_iterations;  // Remaining reads to be generated in each while loop
                                    int local_good_iterations{0};                                       // Thread-local counter for iterations that generated reads
                                    #pragma omp for schedule(dynamic) reduction(+:reads_actually_generated)           // This is where we are splitting the iterations of the for loop to each thread with a reduction variable to specify that it is shared
                                    for(int i=0; i<target_reads; i++){                   
                                        int threadID = omp_get_thread_num();                                          // Get the ID of each thread being tracked
                                        if(ART::generate_paired_reads_with_indel(read1, read2, min_DNA_fragment_length, fragment_weights, threadID)){
                                                                                                                      // Make two paired-reads from the same DNA fragment with indel errors
                                            // Process read 1 first
                                            std::vector<short> read1_quality_score_vec;                               // Vector to hold the quality scores for read 1
                                            read1.get_read_quality(read1_quality_score_vec, 1, threadID);             // Get the read quality scores for the read positions for read 1
                                            read1.add_baseCall_error(read1_quality_score_vec, threadID);              // Add base call errors to the read based on the quality scores on read 2
                                            
                                            std::string cellName = cellGenomes[cell]; cellName.erase(cellName.length()-3);
                                            std::string chromID = chromSegSeq_ID; chromID.erase(0,1);                 // Remove the '>' symbol from the chrom ID
                                            std::string read1_data ="@"+cellName+"_"+chromID+"_read"+std::to_string(reads_actually_generated)+"/1\n";   // @readID
                                            read1_data +=(*read1.get_final_read_sequence(threadID))+"\n+\n";                                            // read sequence and +
                                            for(size_t k=0; k<(*read1.get_final_read_sequence(threadID)).size(); k++){                                  // read quality scores
                                                read1_data +=static_cast<char>(read1_quality_score_vec[k]+32);                                          // +33 to get the phred score
                                            }
                                            read1_data +="\n";
                                            
                                            // Process read 2
                                            std::vector<short> read2_quality_score_vec;                               // Vector to hold the quality scores for read 2
                                            read2.get_read_quality(read2_quality_score_vec, 2, threadID);             // Get the read quality scores for the read positions for read 2
                                            read2.add_baseCall_error(read2_quality_score_vec, threadID);              // Add base call errors to the read based on the quality scores on read 2
                                            
                                            std::string read2_data ="@"+cellName+"_"+chromID+"_read"+std::to_string(reads_actually_generated)+"/2\n";   // @readID
                                            read2_data +=(*read2.get_final_read_sequence(threadID))+"\n+\n";                                            // read sequence and +
                                            for(size_t k=0; k<(*read2.get_final_read_sequence(threadID)).size(); k++){                                  // read quality scores
                                                read2_data +=static_cast<char>(read2_quality_score_vec[k]+32);                                          // +33 to get the phred score
                                            }
                                            read2_data +="\n";
                                            
                                            reads_actually_generated += 2;
                                            local_good_iterations++;

                                            batch_buffer_r1[threadID].push_back(read1_data);                          // Add the read 1 data to the batch buffer 1.   
                                            batch_buffer_r2[threadID].push_back(read2_data);                          // Add the read 2 data to the batch buffer 2.                 
                                            if (batch_buffer_r1[threadID].size() >= static_cast<size_t>(batchSize_thread)){// Check if the batch buffer is full, and write it to the file if needed.
                                                #pragma omp critical(section1)
                                                {
                                                    writeBatchToFile(batch_buffer_r1[threadID], fastq_R1_file, true);
                                                    writeBatchToFile(batch_buffer_r2[threadID], fastq_R2_file, true);
                                                }
                                            }
                                        }
                                    }
                                    #pragma omp atomic                                                                     // Ensure only one thread modifies the shared variable at a time
                                    good_iterations += local_good_iterations;
                                    #pragma omp barrier
                                }
                            } 
                        }
                        #pragma omp single                                                              // Only one thread should increment the counters
                        {
                            entries_in_freqMap--;                                                       // One entry in the map is processed, so decrement the number of entries remaining in the map
                            line_inFile +=2;                                                            // Increment the number of lines read from the fasta memory map                            
                        }    
                    }else{                                                                              // If the line currently got from the memory map need not be sequenced
                        #pragma omp single
                        {
                            line_inFile +=2;                                                            // Increment the number of lines read from the fasta memory map                            
                        }
                    }
                    #pragma omp barrier                                                                 // Wait till all the threads reach here
                    line_num = line_inFile;                                                             // Update the local counter so that all threads get the updated value
                }
            }
            munmap(fastaFileMM, fastaFileSize);                                                         // Unmap the fasta file of the current cell to avoid memory-leaks
    
            report_cells_sequenced.push_back(cellGenomes[cell]);                                        // Storing the names of the cells that were actually sequenced, for the report
            report_readsGenerated_perCell.push_back(reads_actually_generated);                          // Storing the number of reads generated per cell in a vector for the final summary report
        }
     
        for (size_t l=0;l<batch_buffer_r1.size();l++){
            writeBatchToFile(batch_buffer_r1[l], fastq_R1_file, true);                                  // If there are unwritten data in batch buffer 1, write that too when the loop ends
            writeBatchToFile(batch_buffer_r2[l], fastq_R2_file, true);                                  // If there are unwritten data in batch buffer 2, write that too when the loop ends
        }

        fastq_R2_file.close();
    }
    fastq_R1_file.close();

    report_fastq_output.push_back(*parameter.get_output_fastq_filename_prefix());                       // Output filename without the extension for the report
}
//--------------------------------------------------------------------------------------------