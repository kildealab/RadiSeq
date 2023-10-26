#include <iostream>
#include <string>
#include <sys/stat.h>
#include <fstream>
#include <ctime>
#include <unistd.h>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <sys/mman.h>

#include "support_functions.h"
#include "fileio.h"
#include "parameter_handler.h"
#include "sddfile_handler.h"
#include "random_generator.h"
#include "fastafile_handler.h"
#include "sequencing.h"
#include "summary_report.h"


int main(int argc, char* argv[]){
    
    clock_t start_time, end_time;                                                                               // Variables to hold start and end times of this program
    start_time = clock();                                                                                       // Getting the starting time of the program
    
    const char* dataFolder = std::getenv("RADISEQ_DATA_DIR");                                                   // Getting the environment variable "RADISEQ_DATA" which holds the position of the radiSeqData folder path
    std::string dataFolderPath;                                                                                 // Variable to hold the string value of the environment variable
    if (dataFolder != nullptr){                                                                                 // Check if the environment variable is set. If not, exit with error message
        dataFolderPath = dataFolder;                                                                            // Convert the environment variable to a C++ string
    }else{                                                                                
        std::cerr<<"\n ERROR: The environment variable \"RADISEQ_DATA_DIR\" is not set correctly \n"
                 <<" Refer to the installation instructions in README \n";
        exit(EXIT_FAILURE);
    }

    ascii_art();                                                                                                // Prints the program name in output       
    checkArgument(argc, argv);                                                                                  // Check if parameter file is provided as an argument. If not, print error and exit. 
    const std::string user_parameter_file{argv[1]};                                                             // Defining the path to the user-specified parameter file
    
    //-------------- Stage 1: Reading the UserParameter file -----------------//
    NGSParameters parameters;                                                                                   // Initializing parameters of the class NGSParameters
    std::cout<<"\n ----- Initiating parameter file processing ----- \n";
    parameters.process_parameterFile(&user_parameter_file, parameters, &dataFolderPath);                        // Set sequencing parameters using user-specified parameter file. Undefined parameters will default
    std::cout<<"\n Successfully completed the parameter file processing \n";

    //-------------- Stage 2: Reading all the SDD files -----------------//
    NGSsdd SDDdata;                                                                                             // Initializing SDDdata of the class NGSsdd
    std::cout<<"\n ----- Initiating SDD file processing ----- \n";
    SDDdata.process_sddfile(SDDdata, parameters);                                                               // Process the data from all the SDD files provided

    std::cout<<"\n Found "<<SDDdata.get_num_of_exposures()<<" damaged cells in the SDD file(s)\n";
    std::string output_directory = *parameters.get_output_directory();                                          // Output directory path for all output files
    if(!checkFolderExists(output_directory.c_str())){                                                           // If output directory does not exist already
        mkdir(output_directory.c_str(), 0700);                                                                  // Create one; 0700:only owner will have read, write, and execute permissions on the created directory.
    }
    // Creating a temporary directory to store damaged cell genome FASTA files that will get generated
    std::string tempFolderPath = output_directory + "/temp";                                                    // Create a path for the temporary folder
    mkdir(tempFolderPath.c_str(), 0700);                                                                        // Create a temp directory for individual fasta files for damaged cells
    
    std::cout<<"\n Constructing an un-damaged cell model\n";
    // Generating an undamaged cell genome template. This is used later to create damaged cell genomes. 
    long ref_genomeFile_size = fileSize_bytes(*parameters.get_reference_genome());                              // Find the size of the reference sequence file
    std::string genomeTemplatePath = tempFolderPath+"/Undamaged_cell.fa";                                       // Name of the undamaged fasta file template
    size_t templateSize = static_cast<size_t>(ref_genomeFile_size*4);                                           // The size of an Undamaged file is estimated to be 4 times the size of the reference sequence file
    char* genomeTemplate_data = createMemoryMappedFile(genomeTemplatePath,templateSize);                        // Generate a memory-map placeholder to store the memory map of the undamaged fasta file as it gets created later
    // Build the UndamagedGenomeTemplate file and the memory map
    long ref_seq_length = buildUndamagedGenomeTemplate_MM(genomeTemplate_data, templateSize, SDDdata.get_num_chrom(), SDDdata.get_chrom_mapping(), parameters.get_reference_genome());
    long one_fasta_size = fileSize_bytes(tempFolderPath+"/Undamaged_cell.fa");                                  // Calculate the size of the undamaged cell fasta file. This will be the size of every fasta file
    checkStorageSize(parameters, SDDdata, one_fasta_size);                                                      // Check if there is enough storage space to run this program with the given parameters
    // Make sure the difference between the reference genome length and the MC model length is within the required limit
    double percent_diff_seq_length = ((std::abs(SDDdata.get_sdd_genome_length()-ref_seq_length))/ref_seq_length)*100;
    if(*parameters.get_sequencer() != "test"){                                                                  // For all scenarios other than the test run,
        if(percent_diff_seq_length>parameters.get_max_acceptable_seq_length_difference()){                      // If the reference seq length and the monte carlo model seq length are different more than the value specified
            std::cerr<<"\n ERROR: The reference sequence length ("<<ref_seq_length<<" bp) and "
                     <<"the Monte Carlo model genome length("<<SDDdata.get_sdd_genome_length()<<" bp) \n"
                     <<" are significantly different (>"<<parameters.get_max_acceptable_seq_length_difference() <<"%) \n";
            exit(EXIT_FAILURE);
        }else{
            if (SDDdata.get_sdd_genome_length()>ref_seq_length){                                                // If the reference sequence is smaller than the MC model's genome
                std::cerr<<"\n WARNING: Reference sequence length ("<<ref_seq_length<<" bp) is smaller "
                         <<"than the Monte Carlo model genome("<<SDDdata.get_sdd_genome_length()<<" bp). \n"
                         <<" So damages beyond the reference sequence length will be ignored from sequencing. \n";
            }else{                                                                                              // If the reference sequence is smaller than the MC model's genome
                std::cerr<<"\n WARNING: Reference sequence length ("<<ref_seq_length<<" bp) is bigger "
                         <<"than the Monte Carlo model genome("<<SDDdata.get_sdd_genome_length()<<" bp). \n";
            }
        }
    }
    std::cout<<"\n Successfully completed the SDD file processing \n";

    // Read each cell (exposure) damage data from the SDD file, adjust the damages according to the relative dose and actual dose delivered if necessary,
    // then combine multiple radiation damages on the same cell if needed, find DSB locations and then build a damaged genome FASTA file for each cell   
    std::vector<double>& rel_dose = parameters.get_relative_dose_contributions();
    std::vector<std::string>& sdd_paths = parameters.get_sddfile_path();                          
    std::vector<int> lines_in_cell_files{};                                                                     // Vector to store the total number of lines written in each damage cell's fasta file
    std::cout<<"\n ----- Building damaged genomes of the irradiated cells ----- \n";
    
    
    for(int i=0; i<SDDdata.get_num_of_exposures(); i++){                                                        // Iterate over each exposure (cell) data
        for(int j=0; j<SDDdata.get_num_of_SDDs(); j++){                                                         // Iterate through every SDD file (damage files for a cell) given for the same exposure (cell)
            readSDDfileData(&sdd_paths[j], SDDdata, j);                                                         // Read SDD data fields and get damages in each exposure
            SDDdata.adjust_damages_data(rel_dose[j],i,j,parameters.get_adjust_damages_with_actual_dose());      // Adjust the number of damages if needed and store damage to a permanent vector
            SDDdata.reset_temporary_damage_vecs();                                                              // Reset the temporary nested vectors used to hold the damage values before the next SDD file of the same cell
        }
        //SDDdata.find_DNA_breakPoints(parameters.get_dsb_threshold());
        //-------------- Stage 3: Generating damaged cell genomes -----------------//
        std::string fastaFileName = "/Damaged_cell_" + std::to_string(i+1) + ".fa";
        lines_in_cell_files.push_back(buildDamagedCellGenome_from_MM(SDDdata, tempFolderPath, fastaFileName, genomeTemplate_data, templateSize));
        std::cout<<"\n Built damage genomes of "<<std::to_string(i+1)<<" cells \n"; 
        SDDdata.reset_permanent_damage_vecs();                                                                  // Reset all the bigger permanent damage vectors including DNAbreakpoints before processing the next cell
    }
    munmap(genomeTemplate_data, templateSize);                                                                  // Unmap the memory-map to avoid memory leaks after use
    std::cout<<"\n Building of all the damaged cell genomes is now complete \n";

    //-------------- Stage 4: Integrating ART pipeline -----------------//
    
    // Randomly sample 'num_of_cells_to_sequence' from 'num_of_cells_in_sample'. First generate a random number in the range of number of cells in sample. If that number is less than the number of damaged cells, then add the 
    // damaged cell fasta filename to the list. If that cell is already added, then repeat. If the random number is greater than the num of damaged cells, then add the undamaged fasta filename to the list. 
    std::vector<std::string> cellGenomes_to_be_sequenced{};                                                     // Vector to store the fasta filenames of the cells to be sequenced
    std::vector<int> lines_in_cells_to_sequence{};                                                              // Vector to store the total number of lines in the fasta files of cells to be sequenced in order
    while (cellGenomes_to_be_sequenced.size() < parameters.get_num_of_cells_to_sequence()){                     // Iterate till enough number of cells are randomly sampled from the pool of cells for sequencing    
        int cell_id = rng::rand_int(1, parameters.get_num_of_cells_in_sample());                                // Randomly pick an ID for a cell to be sequenced from the sample     
        if (cell_id <= SDDdata.get_num_of_exposures()){                                                         // If the ID corresponds to a damaged cell, then
            std::string damaged_fasta_filename = "Damaged_cell_" + std::to_string(cell_id) + ".fa";
            if (std::find(cellGenomes_to_be_sequenced.begin(), cellGenomes_to_be_sequenced.end(), damaged_fasta_filename) != cellGenomes_to_be_sequenced.end()){
                continue;                                                                                       // Check if the filename is already in the vector. If yes, continue and pick another cell ID
            }else{
                cellGenomes_to_be_sequenced.push_back(damaged_fasta_filename);                                  // Else, add the corresponding damaged cell name to the list
                lines_in_cells_to_sequence.push_back(lines_in_cell_files[cell_id-1]);                           // Add the number of lines that damaged cell file has to the list
            }
        }else{                                                                                                  // If the picked cell ID is beyond the number of damaged cells, add undamaged genome fasta to the list
            cellGenomes_to_be_sequenced.push_back("Undamaged_cell.fa");
            lines_in_cells_to_sequence.push_back(SDDdata.get_num_chrom()*4);                                    // Undamaged cell file will have a total number of lines = 4 times the number of chromosomes
        }
    }
    
    // Perfrom single-cell or bulk-cell sequencing as required
    if (*parameters.get_sequencing_mode() == "single"){
        std::cout<<"\n ----- Initiating Single-cell sequencing of "<<parameters.get_num_of_cells_to_sequence()<<" cells -----\n";
        single_cell_sequencing(parameters, cellGenomes_to_be_sequenced, ref_seq_length);
        std::cout<<"\n Single-cell sequencing of cells are now complete. You can find the sequenced FASTQ files in the output folder: "<<output_directory<<"\n";
    }else{                                                                                                      // If not single, then the other option is only bulk
        std::cout<<"\n ----- Initiating Bulk-cell sequencing of "<<parameters.get_num_of_cells_to_sequence()<<" cells -----\n";
        bulk_cell_sequencing(parameters, cellGenomes_to_be_sequenced, lines_in_cells_to_sequence, ref_seq_length);
        std::cout<<"\n Bulk-cell sequencing of cells are now complete. You can find the sequenced FASTQ files in the output folder: \""<<output_directory<<"\"\n\n";
    }
    
    end_time = clock();                                                                                         // Finding the end time of the run

    // Make the summary report of the run if specified
    if (parameters.get_summary_report()){
        report_cpu_time_used = (static_cast<double>(end_time - start_time)) / CLOCKS_PER_SEC;
        report_parameterFileName = user_parameter_file;
        report_ref_seq_length = ref_seq_length;
        generate_run_summaryReport(parameters, SDDdata);
    }

    // Remove the temporary directory (temp) that stores the fasta file once processing is done
    //remove_directory(tempFolderPath); 
}