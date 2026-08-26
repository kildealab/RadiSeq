#ifndef INDUCE_SEQ_H
#define INDUCE_SEQ_H

#include <vector>
#include <map>
#include "parameter_handler.h"
#include "sddfile_handler.h"

// Object that stores data and functions used for the INDUCE-seq simulation. During a simulation, data for the simulation is stored in 
// vectors such as dsb_locations, which are indexed by groupTID, the id for the group of threads used for that radiation exposure/set of DNA damage. 
// One InduceSeq object simultaneously holds data for multiple simulations and exposures, indexed by different groupTID values. 
// This matches the format of the NGSsdd object sddData. During the program, DNA damage data is to be stored in the sddData attribute,
// and initialized / created for each exposure before the INDUCE-seq simulation is run. The groupTID indices of the sddData data arrays should match
// the groupTID indices of the InduceSeq object. 
// In this object, and in NGSsdd objects, base positions start at 1, and do not reset at chromosome boundaries. 
// Break site positions start at 1 being the site between the first and second base. In general, the break site directly after the nth base has position n. 
class InduceSeq {

public:
    InduceSeq(NGSsdd& sddData, NGSParameters parameters, std::string tempFolderPath); // Sets the parameter member and calls set_genome_data

    void init_set_data_holders(int nGroupThreads);                                   // function to set/initialize all the DSB data holders for processing
    void reset_permanent_damage_vecs(int groupTID);                                  // function to empty the permanent vectors after each exposure

    void run_simulation(int cell_number, int groupTID, int NworkerThreads, int threadIDOffset);
    std::vector<std::vector<long>>& get_dsb_locations(int groupTID);                 // function to get the DSB locations
    void close();                                                                    // function to release resources held by this object (e.g. unmap genome_fasta)

private:
    void set_genome_data(std::string& tempFolderPath);                               // function to set genome data by memory-mapping the reference genome
    void find_DSBs(int DSBthreshold, int groupTID);                                  // function to find DSBs from backbone breaks on opposite strands within a threshold and on the same chromosome
    void save_dsb_locations(int cell_number, int groupTID);                          // function to save every element of dsb_locations[groupTID] to a csv file
    void load_blunted_ends(int groupTID);                                            // function to determine blunted ends from the DSB locations
    void save_dsb_blunted_ends(int cell_number, int groupTID);                       // function to save every element of dsb_blunted_ends[groupTID] to a csv file
    void load_dsb_fragments(int groupTID, int threadID);                             // function to calculate and save the DSB fragments from the blunted ends
    void filter_dsb_fragments(int groupTID, int threadID);                           // function to filter dsb_fragments_left/right by fragment size, using probability_of_sequencing_function, called after load_dsb_fragments
    void filter_dsb_strands_ssd(int groupTID);                                       // function to filter DSB strands by single-strand damage
    void find_base_pair_damages(int groupTID);                                       // function to find base pair damages on the DSB strands
    void get_dna_sequence(std::string& dna_seq, std::vector<long>& bp_damages, std::vector<long>& dsb_strand, bool is_left);  // function to extract the DNA sequence for a DSB strand from the genome
    void generate_simulation_output(int cell_number, int groupTID, int num_available_threads, int threadIDOffset);  // function to induce sequencing from DSB fragments
    int get_random_fragment_length(int threadID);                                    // function to sample a fragment length from fragment_size_distribution
    void set_fragment_size_distribution_from_file();                                 // function to read the induce_seq fragment size distribution file and populate fragment_size_distribution
    void set_probability_of_keeping_from_file();                                     // function to read the induce_seq probability of sequencing file and populate probability_of_sequencing_function

    // Data used and created during an INDUCE-seq simulation. All of these are indexed by groupTID, as explained in the InduceSeq description. 
    std::vector<std::vector<std::vector<long>>> dsb_locations;                       // Stores double strand breaks. Each entry is {strand1_break_location, strand2_break_location, chromosome_index, is_previous_step_dsb}, where chromosome_index starts from 0, and is_previous_step_dsb marks whether the previous entry is part of the same interconnected DSB cluster (see find_DSBs)
    std::vector<std::vector<std::vector<long>>> dsb_blunted_ends;                    // Stores blunted DSB ends, one entry per set of 2 DSB ends (which can be the result of 1 DSB or an interconnected cluster of DSBs). Each entry is {left_edge_location, right_edge_location, chromosome_index, left_dsb_strand1_location, left_dsb_strand2_location, right_dsb_strand1_location, right_dsb_strand2_location}, 
                                                                                     // where the left/right dsb_strand fields are the strand1/strand2 break locations of the (possibly different) DSBs that caused the left and right edges (see load_blunted_ends)
    std::vector<std::vector<std::vector<long>>> dsb_fragments_left;                  // Stores left DSB fragments, one entry per fragment. Each entry is {blunted_end_location, other_end_location, chromosome_index, causing_dsb_strand1_location, causing_dsb_strand2_location}; for left fragments element 0 > element 1 (see load_dsb_fragments). After filter_dsb_strands_ssd runs, entries that don't survive single-strand-damage filtering are removed, leaving only the surviving denatured DNA strands (same per-entry format)
    std::vector<std::vector<std::vector<long>>> dsb_fragments_right;                 // Stores right DSB fragments, in the same per-entry format as dsb_fragments_left, except element 1 > element 0 (see load_dsb_fragments). Filtered in place by filter_dsb_strands_ssd, as with dsb_fragments_left
    std::vector<std::vector<std::vector<long>>> base_pair_damages_left;              // Stores base pair damages on the left DSB strands. Indexed the same as dsb_fragments_left (post single-strand-damage filtering): entry i is the (possibly empty) list of base-pair damage locations found on dsb_fragments_left's i-th strand
    std::vector<std::vector<std::vector<long>>> base_pair_damages_right;             // Stores base pair damages on the right DSB strands. Indexed the same as dsb_fragments_right (post single-strand-damage filtering): entry i is the (possibly empty) list of base-pair damage locations found on dsb_fragments_right's i-th strand
    
    
    std::map<float, int> fragment_size_distribution;
    std::map<int, double> probability_of_sequencing_function;                        // Length -> probability of keeping (retaining) a fragment of that length, read from the induce_seq probability of sequencing file
    NGSParameters parameter;                                                         // Holds the simulation parameters
    NGSsdd& sdd_data;                                                                // Reference to the shared NGSsdd instance holding genome/backbone break data
    std::vector<long> chrom_end_loc;                                                 // Copy of *sdd_data.get_chrom_end_loc(), set at construction time
    char* genome_fasta{nullptr};                                                     // Path to the reference genome FASTA file
    size_t genome_fasta_size{0};                                                     // Size in bytes of the genome_fasta memory map
    std::vector<int> cum_chrom_header_sizes;                                         // Cumulative chromosome header sizes
    std::vector<std::string> chrom_headers;                                          // Chromosome headers, excluding the trailing '\n'
    std::string P7_adapter_sequence;                                                 // P7 adapter sequence, appended to a read (with any 'I' placeholders replaced by random bases) when its DSB fragment is shorter than the read length
};
#endif
