#ifndef SDDFILE_HANDLER_H
#define SDDFILE_HANDLER_H

#include "parameter_handler.h"

class NGSsdd{
    int num_of_SDDs;                                                                 // Variable to hold the number of SDD files/ particles to merge
    int num_of_exposures;                                                            // Variable to hold the number of exposures in SDD. This will be the lowest value among all SDDs
    int num_of_damagedCells_toBuild;                                                 // Variable to hold the number of damaged cells we want to build from all the exposures in the SDD file. This value will be <= num of exposures
    std::vector<std::vector<double>> actual_dose_delivered;                          // Vector that holds vectors than contain actual dose delivered in each exposure for every SDD file
    std::vector<double> expected_dose_gy;                                            // Vector to hold the dose in Gy that was expected to be delivered in each SDD file
    std::vector<long> chrom_size_bp;                                                 // Vector that holds the chromosome length in bp from the MC model
    std::vector<std::streampos> lines_read_sdd;                                      // This vector will hold the location in the SDD file where we last left off; for all SDDs
    std::vector<std::vector<std::streampos>> exposureLine_vec;                       // This vector will hold vectors that has the positions of new exposure data lines in each SDD files. 
    std::vector<std::vector<long>> damage_start_loc_vec;                             // Vector holding the location of the start of the damage block within a chromosome for the line that is processing
    std::vector<std::vector<int>> chrom_ID_vec;                                      // Vector to hold the chromosome ID of the damage block corrently processing. Starts from zero   
    std::vector<std::vector<std::vector<long>>> temp_backbone1_breaks_vec;           // This temporary vector will hold all the locations of breaks in backbone 1 for one exposure in 1 SDD
    std::vector<std::vector<std::vector<long>>> temp_basestrand1_damages_vec;        // This temporary vector will hold all the locations of base damages in the 5'-3' strand for one exposure in 1 SDD
    std::vector<std::vector<std::vector<long>>> temp_basestrand2_damages_vec;        // This temporary vector will hold all the locations of base damages in the 3'-5' strand for one exposure in 1 SDD
    std::vector<std::vector<std::vector<long>>> temp_backbone2_breaks_vec;           // This temporary vector will hold all the locations of breaks in backbone 1 for one exposure in 1 SDD
    std::vector<int> original_num_damages;                                           // Vector to hold the original number of damages in one exposure from one SDD before any adjustments. One value per thread
    std::vector<std::vector<long>> backbone1_break_loc_vec;                          // Vector will hold all the locations of breaks in backbone 1. Each element will be indicating break location in the genome. From all SDDs in 1 exposure
    std::vector<std::vector<long>> basestrand1_damage_loc_vec;                       // Vector will hold all the locations of base damages in the 5'-3' strand. From all SDDs in 1 exposure 
    std::vector<std::vector<long>> basestrand2_damage_loc_vec;                       // Vector will hold all the locations of base damages in the 3'-5' strand. From all SDDs in 1 exposure
    std::vector<std::vector<long>> backbone2_break_loc_vec;                          // Vector will hold all the locations of breaks in backbone 1. Each element will be indicating break location in the genome. From all SDDs in 1 exposure
    std::vector<long> chrom_end_loc{0};                                              // Vector stores the location of the ends of each chromosome in genome (in unit of bp). First element is zero. 
    int num_chroms;                                                                  // Variable sets the number of chromosomes in the cell model
    int chrom_mapping_type{0};                                                       // 0 if haploid, 1 if chrom mapping is 1,1,2,2......22,22,X,Y; 2 if chrom mapping is 1,2,......22,1,2,......22,X,Y
    //std::vector<long> DNA_breakPoints;                                               // Vector to hold all the DSB locations of each cell at a time


public:
    NGSsdd();                                                                        // Default constructor
    void process_sddfile(NGSsdd&, NGSParameters&);                                   // Function to process all the SDD files

    void set_num_of_SDDs(int);                                                       // function to set 'num_of_SDDs'
    int get_num_of_SDDs();                                                           // function to get 'num_of_SDDs'

    void set_num_of_exposures(int);                                                  // function to set 'num_of_exposures'
    int get_num_of_exposures();                                                      // function to get 'num_of_exposures'

    void set_num_of_damagedCells_toBuild(int);                                       // function to set 'num_of_damagedCells_toBuild'
    int get_num_of_damagedCells_toBuild();                                           // function to get 'num_of_damagedCells_toBuild'

    const std::vector<std::vector<double>>* get_actual_dose_delivered();             // function to get 'actual_dose_delivered' vector

    void set_expected_dose_gy(std::string*);                                         // function to set expected dose in Gray value
    double get_expected_dose_gy(int);                                                // function to get expected dose in Gray

    void set_chrom_size_bp(std::string*);                                            // function to set the chromosone sizes in bp from SDD
    const std::vector<long>* get_chrom_size_bp();                                    // function to get 'chrom_size_bp' vector

    void set_cell_ploidy_and_chrom_mappping(std::string*);                           // function to calculate the ploidy of the cell (if there are paired chrosomosomes) and how chrom is arranged in SDD file if paired
    int get_chrom_mapping();                                                         // function to get the chromosome mapping used to generate the SDD file if diploid
    int get_num_chrom();                                                             // functio to get 'num_chroms'

    void set_chrom_end_loc();                                                        // function to set the chromosome end locations
    const std::vector<long>* get_chrom_end_loc();                                    // funtion to get the vector with chromosome end locations
    long get_sdd_genome_length();                                                    // funtion to get the total legth of the genome of the cell model from SDD

    void set_lines_read_sdd(std::streampos, int);                                    // function to set the lines_read_sdd vector
    std::vector<std::streampos>* get_lines_read_sdd();                               // function to get 'ines_read_sdd' vector
    std::streampos get_exposureLines(int, int);                                      // function will return one line position from the vector that stores new expsure lines of each SDD file

    void init_set_data_holders(int);                                                 // function to set/initialize all the SDD data holders for processing
    void init_set_workThread_data_holders(int, int);                                 // function to set/initialize all the temporary data holders for the worker threads
    void reset_workThread_data_holders(int);                                         // function to empty all the temporary worker thread data holders
    void merge_workerThread_vectors(int);                                            // function to concatenate vectors that store worker-thread wise information into a single vector for the group
    
    void set_all_sdd_data_fields(std::string*, int, int);                            // this is a mother set function that will call individual set functions for the SDD data fields 
    void set_sddData_field3(std::string*, int, int);                                 // function to process field3 sddData
    void set_sddData_field4(std::string*, int, int);                                 // function to process field4 sddData
    void set_sddData_field7(std::string*, int, int);                                 // function to process field7 sddData

    int get_original_num_damages(int);                                               // function to get the original number of damages
    
    void reset_temporary_damage_vecs(int);                                           // function to empty the temporary vectors after each SDD
    void reset_permanent_damage_vecs(int);                                           // function to empty the permanent vectors after each exposure
    void adjust_damages_data(int, int, bool, int);                                   // function to remove damages according to the actual dose if needed
    
    //void find_DNA_breakPoints(int);                                                  // function to find DSBs and Chromosome ends to divide the genome into fragments
    //std::vector<long>& get_DNA_breakPoints();                                        // function to get the DNA break points for each cell in process at the moment

    std::vector<long>& get_backbone1_break_loc(int);                                 // function to get the backbone 1 break locations
    std::vector<long>& get_backbone2_break_loc(int);                                 // function to get the backbone 2 break locations
    std::vector<long>& get_basestrand1_damage_loc(int);                              // function to get the base strand 1 damage locations
    std::vector<long>& get_basestrand2_damage_loc(int);                              // function to get the base strand 2 damage locations

};

#endif