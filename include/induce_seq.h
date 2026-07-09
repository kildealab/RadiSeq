#ifndef INDUCE_SEQ_H
#define INDUCE_SEQ_H

#include <vector>
#include <map>
#include "parameter_handler.h"
#include "sddfile_handler.h"

class InduceSeq {

public:
    InduceSeq(NGSsdd& sddData);                                                       // Lightweight constructor: binds sdd_data only. Use when induce_seq is not requested, to avoid the cost of set_genome_data
    InduceSeq(NGSsdd& sddData, NGSParameters parameters, std::string tempFolderPath); // Also sets the parameter member and calls set_genome_data, which is computationally intensive

    void init_set_data_holders(int nGroupThreads);                                   // function to set/initialize all the DSB data holders for processing
    void reset_permanent_damage_vecs(int groupTID);                                   // function to empty the permanent vectors after each exposure

    void set_parameter(NGSParameters param);                                         // function to set the parameter member
    void set_genome_data(std::string& tempFolderPath);                               // function to set genome data by memory-mapping the reference genome
    void find_DSBs(int DSBthreshold, int groupTID);                                  // function to find DSBs from backbone breaks on opposite strands within a threshold and on the same chromosome
    void get_blunted_ends(int groupTID);                                             // function to determine blunted ends from the DSB locations
    void get_dsb_fragments(int groupTID, int threadID);                              // function to get the DSB fragments from the blunted ends
    void filter_dsb_strands_ssd(int groupTID);                                       // function to filter DSB strands by single-strand damage
    void find_base_pair_damages(int groupTID);                                       // function to find base pair damages on the DSB strands
    void get_dna_sequence(std::string& dna_seq, std::vector<long>& bp_damages, std::vector<long>& dsb_strand, bool is_left);  // function to extract the DNA sequence for a DSB strand from the genome
    void generate_simulation_output(int cell_number, int groupTID, int num_available_threads, int threadIDOffset);  // function to induce sequencing from DSB fragments
    void run_simulation(int cell_number, int groupTID, int threadID, int NworkerThreads, int threadIDOffset);
    std::vector<std::vector<long>>& get_dsb_locations(int groupTID);                 // function to get the DSB locations
    int get_random_fragment_length(int threadID);                                    // function to sample a fragment length from fragment_size_distribution
    void close();                                                                    // function to release resources held by this object (e.g. unmap genome_fasta)

private:
    void set_fragment_size_distribution_from_file();                                 // function to read the induce_seq fragment size distribution file and populate fragment_size_distribution

    std::vector<std::vector<std::vector<long>>> dsb_locations;                       // Stores double strand breaks. First index is groupTID, each element is a list of [backbone1_site, backbone2_site, chrom indx]. chrom indx starts from 0, 1, 2, ...
    std::vector<std::vector<std::vector<long>>> dsb_blunted_ends;                    // Stores blunted DSB ends. First index is groupTID
    std::vector<std::vector<std::vector<long>>> dsb_fragments_left;                  // Stores left DSB fragments. First index is groupTID
    std::vector<std::vector<std::vector<long>>> dsb_fragments_right;                 // Stores right DSB fragments. First index is groupTID
    std::vector<std::vector<std::vector<long>>> dsb_strands_left;                    // Stores left DSB strands. First index is groupTID
    std::vector<std::vector<std::vector<long>>> dsb_strands_right;                   // Stores right DSB strands. First index is groupTID
    std::vector<std::vector<std::vector<long>>> base_pair_damages_left;              // Stores base pair damages on the left DSB strands. First index is groupTID
    std::vector<std::vector<std::vector<long>>> base_pair_damages_right;             // Stores base pair damages on the right DSB strands. First index is groupTID
    std::map<float, int> fragment_size_distribution{{0.2f, 5}, {0.6f, 6}, {1.0f, 7}};
    NGSParameters parameter;                                                         // Holds the simulation parameters
    NGSsdd& sdd_data;                                                                // Reference to the shared NGSsdd instance holding genome/backbone break data
    char* genome_fasta{nullptr};                                                     // Path to the reference genome FASTA file
    size_t genome_fasta_size{0};                                                     // Size in bytes of the genome_fasta memory map
    std::vector<int> cum_chrom_header_sizes;                                         // Cumulative chromosome header sizes
    std::vector<std::string> chrom_headers;                                          // Chromosome headers, excluding the trailing '\n'


    int second_size_filter{3};                                                       // Second size filter threshold
    

};
#endif
