#include "induce_seq.h"
#include "art_framework.h"
#include "fastafile_handler.h"
#include "fileio.h"
#include "random_generator.h"
#include "summary_report.h"

#include <fstream>
#include <iostream>
#include <map>
#include <algorithm>
#include <cmath>
#include <omp.h>
#include <sys/mman.h>


InduceSeq::InduceSeq(NGSsdd& sddData, NGSParameters parameters, std::string tempFolderPath) : parameter(parameters), sdd_data(sddData) {
    //genome data is only needed for generating reads
    if (parameter.get_generate_reads()) {
        set_genome_data(tempFolderPath);
    }
    set_fragment_size_distribution_from_file();
    set_probability_of_keeping_from_file();
    P7_adapter_sequence = *parameter.get_P7_adapter_sequence();

    if (parameter.get_generate_reads()) {
        // Initialize parameters for read generation. GC_binSize and fraction_nonFR_read_pairs are not used for Induce Seq, but are neaded as parameters for the function
        ART::initiate_read_generation(parameter.get_read_length(), parameter.get_GC_binSize(), parameter.get_fraction_nonFR_read_pairs(), parameter.get_read_artifacts_rate());
        // r2_quality_profile is not used in InduceSeq, it is needed as an argument for the function
        ART::set_read_quality_distribution(*parameter.get_r1_quality_profile(), *parameter.get_r2_quality_profile());
    } else {
        // If reads are generated, then chrom_end_loc are determined from the reference genome.
        // In this case, they must be read from the SDD file
        chrom_end_loc = *sddData.get_chrom_end_loc();
        // No genome FASTA is loaded in this mode (set_genome_data is never called above), so there is no reference genome file to report
        report_reference_genome_used = "Not applicable (read generation is disabled)";
        report_ref_seq_length = chrom_end_loc.back();
    }
}



// Reads the DNA fragment size distribution file (same format as fragment_size_distribution_path)
// from the path set by the dsb_end_fragment_size_distribution_path parameter, and converts the returned
// (min_fragment_length, normalized_counts) pair into a cumulative-probability -> length lookup map, replacing
// the hard-coded default fragment_size_distribution. get_random_fragment_length samples from this map.
void InduceSeq::set_fragment_size_distribution_from_file() {
    auto fragmentData = readFragmentSizeDist(parameter.get_dsb_end_fragment_size_distribution_path());
    int min_fragment_length = fragmentData.first;
    const std::vector<double>& normalized_counts = fragmentData.second;

    fragment_size_distribution.clear();
    double cumulative_probability = 0.0;
    for (size_t i = 0; i < normalized_counts.size(); i++) {
        cumulative_probability += normalized_counts[i];                            // Running sum of the normalized counts turns the PMF into a CDF
        fragment_size_distribution[static_cast<float>(cumulative_probability)] = min_fragment_length + static_cast<int>(i);
    }
}

// Reads the induce_seq probability-of-sequencing file (see readProbabilityOfSequencing) from the path
// set by the induce_seq_probability_of_sequencing_path parameter, and stores it as a length -> probability
// lookup, probability_of_sequencing_function.
void InduceSeq::set_probability_of_keeping_from_file() {
    probability_of_sequencing_function = readProbabilityOfSequencing(parameter.get_induce_seq_probability_of_sequencing_path());
}


std::vector<std::vector<long>>& InduceSeq::get_dsb_locations(int groupTID) {
    return(dsb_locations[groupTID]);
}

// Sets the shape for data vectors, making each vector contain nGroupThreads empty vectors. 
// This is because the data held in InduceSeq class is indexed by thread group ID, in the same way as NGSsdd. 
// Each thread group processes one damaged genome that is to be sequenced, and accesses only the data corresponding to it. 
void InduceSeq::init_set_data_holders(int nGroupThreads){
    dsb_locations.clear();
    dsb_locations.resize(nGroupThreads);

    dsb_blunted_ends.clear();
    dsb_blunted_ends.resize(nGroupThreads);

    dsb_fragments_left.clear();
    dsb_fragments_left.resize(nGroupThreads);

    dsb_fragments_right.clear();
    dsb_fragments_right.resize(nGroupThreads);

    base_pair_damages_left.clear();
    base_pair_damages_left.resize(nGroupThreads);

    base_pair_damages_right.clear();
    base_pair_damages_right.resize(nGroupThreads);
}

//clears all data for a given groupTID. 
void InduceSeq::reset_permanent_damage_vecs(int groupTID){                                                 // empty all the permanent vectors before the next exposure
    dsb_locations[groupTID].clear();
    dsb_blunted_ends[groupTID].clear();
    dsb_fragments_left[groupTID].clear();
    dsb_fragments_right[groupTID].clear();
    base_pair_damages_left[groupTID].clear();
    base_pair_damages_right[groupTID].clear();
}

//runs the full INDUCE-seq simulation pipeline, including saving output, for a given groupTID.
//cell_number is the used for naming output files results
//NumWorkerThreads is the number of threads available to the given thread group (specified by groupTID). 
//threadIDOffset is the globally-unique thread base index reserved for this group. When multiple threads are spawned, 
//their IDs are determined as threadIDOffset + local_thread_index, where local_thread_index will go from 0 to NumWorkerThreads-1.
//These IDs will be unique across across all concurrent calls of run_simulations, and can be used to access thread-specific data without causing thread-races. 
void InduceSeq::run_simulation(int cell_number, int groupTID, int NumWorkerThreads, int threadIDOffset) {
    find_DSBs(parameter.get_dsb_threshold(), groupTID);
    // save_dsb_locations(cell_number, groupTID);
    load_blunted_ends(groupTID);
    if (parameter.get_output_dsbs()) {
        save_dsb_blunted_ends(cell_number, groupTID);
    }
    //threadIDOffset is used as the threadID for these 2 single-threaded functions, since this will be the global ID of the master (0th) thread in the thread group
    load_dsb_fragments(groupTID, threadIDOffset); 
    filter_dsb_fragments(groupTID, threadIDOffset);
    if (parameter.get_remove_strands_with_SSBs()) {
        filter_dsb_strands_ssd(groupTID);
    }
    find_base_pair_damages(groupTID);
    generate_simulation_output(cell_number, groupTID, NumWorkerThreads, threadIDOffset);
}

// Creates a fasta file containing the genome data in a particular format, as described in buildUndamagedGenomeTemplate_ForwardOnly_MM
// If induce_seq_genome_fasta_path is set and points to an existing file, that previously-built file is loaded instead of
// rebuilding it. If it is set but the file does not exist yet, the file is built and saved to that path (instead of
// tempFolderPath) so that it can be re-used by later runs. If it is unset (default), the file is built fresh in
// tempFolderPath. Sets genome_fasta to point to the first character of the memory map of the fasta file.
// Initializes cum_chrom_header_sizes and chrom_headers (described in the header file). If the chromosome sizes
// listed in the SDD file don't match the constructed genome fasta file, chrom_end_loc (along with
// cum_chrom_header_sizes and chrom_headers) is recomputed from the fasta file's actual chromosome sizes instead,
// after printing a warning, rather than aborting the run (see calculateChromEndLoc).
void InduceSeq::set_genome_data(std::string& tempFolderPath) {
    std::string empty_path = "\"\"";                                                                          // Sentinel value that marks an unset path-type parameter, as read from the default parameter file
    const std::string& savedGenomeFastaPath = *parameter.get_induce_seq_genome_fasta_path();
    bool pathIsSet = (savedGenomeFastaPath != empty_path);

    if (pathIsSet && checkFileExists(&savedGenomeFastaPath)) {
        genome_fasta = generateInputFileMemoryMap(savedGenomeFastaPath, genome_fasta_size);                    // Load the previously saved genome fasta memory-map instead of rebuilding it
        report_reference_genome_used = savedGenomeFastaPath;                                                   // An INDUCE-seq formatted genome already existed, so report that as the genome file used
    } else {
        long ref_genomeFile_size = fileSize_bytes(*parameter.get_reference_genome());
        std::string genomeTemplatePath = pathIsSet ? savedGenomeFastaPath : tempFolderPath+"/genome_spaceless.fa";
        genome_fasta_size = static_cast<size_t>(ref_genomeFile_size*2 + 1000);                                        // The size of an Undamaged file is estimated to be 2 times the size of the reference sequence file, plus potentiall 1000 for chromosome headers 
        std::cout<<"\n Building an induce_seq-formatted reference genome. "<<std::endl;
        if (!pathIsSet) {                                                                                      // Only when building to the temp folder (which is deleted at the end of the run), rather than a persistent, user-specified path
            std::cout<<" This process can be avoided in future runs by saving the formatted genome using the induce_seq_genome_fasta_path parameter"<<std::endl;
        }
        genome_fasta = createMemoryMappedFile(genomeTemplatePath, genome_fasta_size);                          // Generate a memory-map placeholder to store the memory map of the undamaged fasta file as it gets created later
        buildUndamagedGenomeTemplate_ForwardOnly_MM(genome_fasta, genome_fasta_size, sdd_data.get_num_chrom(), sdd_data.get_chrom_mapping(), parameter.get_reference_genome(), cum_chrom_header_sizes);
        report_reference_genome_used = *parameter.get_reference_genome();                                      // A new INDUCE-seq formatted genome was constructed from the unformatted reference genome, so report that instead
        if (pathIsSet) {                                                                                       // Only when saved to a persistent, user-specified path rather than the temp folder
            std::cout<<"\n Saved the induce_seq-formatted reference genome FASTA file to: "<<savedGenomeFastaPath<<" for reuse in future runs"<<std::endl;
        }
    }

    calculateChromEndLoc(cum_chrom_header_sizes, chrom_headers, chrom_end_loc, genome_fasta, genome_fasta_size);
    report_ref_seq_length = chrom_end_loc.back();                                                              // Length of the reference genome, for the summary report
    
    // Make sure the difference between the reference genome length and the MC model length is within the required limit
    long ref_seq_length = chrom_end_loc.back();
    double percent_diff_seq_length = (std::abs(sdd_data.get_sdd_genome_length()-ref_seq_length)/static_cast<double>(ref_seq_length))*100;
                                                               // For all scenarios other than the test run,
    if(percent_diff_seq_length>parameter.get_max_acceptable_seq_length_difference()){                      // If the reference seq length and the monte carlo model seq length are different more than the value specified
        std::cerr<<"\n ERROR: The reference sequence length ("<<ref_seq_length<<" bp) and "
                <<"the Monte Carlo model genome length("<<sdd_data.get_sdd_genome_length()<<" bp) \n"
                <<" are significantly different (>"<<parameter.get_max_acceptable_seq_length_difference() <<"%) \n";
        exit(EXIT_FAILURE);
    }
    

    std::cout<<"\n Successfully completed the SDD file processing "<<std::endl;
}

// finds DSBs from the ssd data, for a given groupTID. Populates the dsb_locations[groupTID] vector (format described in header file).
// DSBthreshold is the maximum distance between a pair of strand breaks for them to qualify as a dsb.
// ssd_data should have its backbone break vectors initialized and sorted least position to greatest position before this function is called.
// Damages beyond the end of the last chromosome (chrom_end_loc.back(), which can be smaller than what the SDD file's
// own damage positions assume if its declared chromosome sizes didn't match the constructed genome fasta; see
// set_genome_data) are not considered: since both break vectors are sorted ascending, the loop simply stops as soon
// as either site reaches such a position, rather than processing the (invalid) remainder of either vector.
void InduceSeq::find_DSBs(int DSBthreshold, int groupTID){
    std::vector<long>& backbone1_breaks = sdd_data.get_backbone1_break_loc(groupTID);
    std::vector<long>& backbone2_breaks = sdd_data.get_backbone2_break_loc(groupTID);
    std::vector<long>::iterator site1 = backbone1_breaks.begin();
	std::vector<long>::iterator site2 = backbone2_breaks.begin();
    const std::vector<long>& chromEnds = chrom_end_loc;                                                         // chrom_end_loc is sorted in ascending order
    long genome_end = chromEnds.back();                                                                         // End position of the last chromosome; no valid damage position can be larger than this

    // variable to keep track of whether the previous step in the while loop was a dsb
    // used for keeping track of dsbs that are part of interconnected dsb, as explained below
    bool prevStepDSB = 0;
    
    //go through backbone breaks, checking for dsbs.
    while (site1 != backbone1_breaks.end() && site2 != backbone2_breaks.end() && *site1 <= genome_end && *site2 <= genome_end){
        int siteDiff = *site1 - *site2;                                                                 // separation in number of bp
		bool isDSB{0};                                                                                  // initiating with zero
        long chromIdx1;
        if(abs(siteDiff) <= DSBthreshold){
            // finds the index of 
            chromIdx1 = std::upper_bound(chromEnds.begin(), chromEnds.end(), *site1) - chromEnds.begin() - 1;  // index of the largest chromEnds element <= *site1 (binary search since chromEnds is sorted)
            long chromIdx2 = std::upper_bound(chromEnds.begin(), chromEnds.end(), *site2) - chromEnds.begin() - 1;
            isDSB = (chromIdx1 == chromIdx2);                                                          // same chromosome if both sites fall in the same interval
        }

        if(isDSB){
            // save dsb to data vector
            dsb_locations[groupTID].push_back({*site1, *site2, chromIdx1, prevStepDSB});

            
            // One site is incremented, to move onto the next strand break. The choice of which site to increment has the form below for the following reason. 
            // Consider the strand breaks as nodes in a graph, with an edge between any two breaks that satifsy the distance criterion for a dsb. 
            // If dsbs are all isolated, then the graph is a set of 2-node connected components, and no special consideration is needed.
            // However, if many breaks are close together, then the graph can include a connected component with many nodes. 
            // In such a component, we are only interested in the first and the last pair of connected nodes, in terms of position in the genome:
            // these two pairs will form the edges of two dna segments, and all the DNA in between them will get all broken up, and will be irrelevent to induce-seq
            // The simplest way to find and record these pairs is with the incrementing method below.
            // This method is gauranteed to detect the first and last connected pair in any connected component, and on every iteration step in between, it will detect a dsb
            // The fourth data field in recorded dsbs is whether the previous step was a dsb. 
            // Thus, when parsing through dsbs, if consecutive dsbs have 1 in that field, they are part of a set of interconnected dsbs.
            // In that case, the last dsb in that sequence and dsb directly before that, which has a 0 in the fourth field, are the dsbs to consider for InduceSeq
            bool site1_has_next = (site1 + 1) != backbone1_breaks.end() && *(site1 + 1) <= genome_end;
            bool site2_has_next = (site2 + 1) != backbone2_breaks.end() && *(site2 + 1) <= genome_end;
            if (site1_has_next && site2_has_next) {
                if ( (*(site1+1) - *site2) < (*(site2+1) - *site1)) {
                    site1++;
                } else {
                    site2++;
                }
            } else if (site1_has_next) {
                // site2 has no more breaks to compare against; advance site1 so its next value can still be checked against the current site2
                site1++;
            } else {
                // Either site2 has a next break to check against the current site1, or neither side has one left (in which case
                // this just advances site1 to end the loop; site2's value has already been fully considered by this point)
                site2++;
            }
        }else{
            // move to the next strand break
            if (siteDiff>0){site2++;}                                                                   // if site1 is after site2; update site2
            else if (siteDiff<0){site1++;}                                                              // if site1 is before site2; update site1
        }
        prevStepDSB = isDSB;
    }
}

// Saves every element of dsb_locations[groupTID] to a csv file, one row per DSB (see find_DSBs for
// how these fields are populated, including what the 'is_previous_step_dsb' field means for
// interconnected DSB clusters).
void InduceSeq::save_dsb_locations(int cell_number, int groupTID) {
    std::string filename = (*parameter.get_output_directory())+"/"+(*parameter.get_output_fastq_filename_prefix())+"_"+std::to_string(cell_number)+"_dsbs.csv";
    std::ofstream dsbs_file(filename.c_str());
    dsbs_file << "strand1_break_location,strand2_break_location,chromosome_index,is_previous_step_dsb\n";
    for (const std::vector<long>& dsb : dsb_locations[groupTID]) {
        dsbs_file << dsb[0] << "," << dsb[1] << "," << dsb[2] << "," << dsb[3] << "\n";
    }
    dsbs_file.close();
}

// Saves every element of dsb_blunted_ends[groupTID] to a csv file, one row per blunted end (see
// load_blunted_ends for how these fields are populated, and for the field format). The
// left_edge/right_edge locations are also reported relative to their chromosome's start (chrom_end_loc[chrom_idx]
// is the cumulative length of all preceding chromosomes, so subtracting it converts a genome-wide location
// to a 1-based location within the chromosome; see load_dsb_fragments for the same conversion).
void InduceSeq::save_dsb_blunted_ends(int cell_number, int groupTID) {
    std::string filename = (*parameter.get_output_directory())+"/"+(*parameter.get_output_fastq_filename_prefix())+"_"+std::to_string(cell_number)+"_dsb_blunted_ends.csv";
    std::ofstream blunted_ends_file(filename.c_str());
    blunted_ends_file << "left edge location,right edge location,left edge location relative to chromosome start,right edge location relative to chromosome start,chromosome index (starting at 0 in the order listed in SDD file)\n";
    for (const std::vector<long>& blunted_end : dsb_blunted_ends[groupTID]) {
        long chrom_start = chrom_end_loc[blunted_end[2]];
        blunted_ends_file << blunted_end[0] << "," << blunted_end[1] << "," << (blunted_end[0] - chrom_start) << "," << (blunted_end[1] - chrom_start) << "," << blunted_end[2] << "\n";
    }
    blunted_ends_file.close();
}

void InduceSeq::close() {
    if (genome_fasta != nullptr) {                                                                             // genome_fasta stays null when generate_reads is False, since set_genome_data is never called
        munmap(genome_fasta, genome_fasta_size);
    }
}

// For a given groupTID, populated the dsb_blunted_ends array to contain all the DNA blunted ends. 
// an entry of dsb_blunted_ends is in the form {location of base on the left edge, location of base on the right edge, chromosome index,
// strand1 (backbone1) location of the dsb that caused the left edge, strand2 (backbone2) location of the dsb that caused the left edge,
// strand1 (backbone1) location of the dsb that caused the right edge, strand2 (backbone2) location of the dsb that caused the right edge}
// A cluster of chained/connected dsbs can have the left edge caused by a different dsb than the right edge, so both are tracked separately.
void InduceSeq::load_blunted_ends(int groupTID) {
    dsb_blunted_ends[groupTID].clear();
    std::vector<std::vector<long>> dsb_locs = get_dsb_locations(groupTID);
    // base locations start at 1, and chromosome indices at 0
    std::vector<long> new_dsb_blunted_ends = {0, 0, 0, 0, 0, 0, 0};
    for (size_t i = 0; i < dsb_locs.size(); i++) {
        std::vector<long> dsb = dsb_locs[i];
        // dsb[3] contains information about compound dsbs (dsbs that form one larger dsb featuring more than 2 single strand breaks)
        // Any new dsb/compound dsb will have the first dsb have dsb[3] = 0. All dsbs belonging to the same compound dsb after that one will have dsb[3] = 1. 
        if (dsb[3] == 0) {
            new_dsb_blunted_ends[0] = dsb[1];
            new_dsb_blunted_ends[3] = dsb[0];                                     // strand1 location of the dsb that starts this cluster (causes the left edge)
            new_dsb_blunted_ends[4] = dsb[1];                                     // strand2 location of the dsb that starts this cluster (causes the left edge)
        }
        // when the next dsb has dsb[3] = 0, it is known that the current dsb is the last dsb in a compound dsb
        if (i + 1 == dsb_locs.size() || dsb_locs[i+1][3] == 0) { 
            new_dsb_blunted_ends[1] = dsb_locs[i][0] + 1;
            new_dsb_blunted_ends[2] = dsb_locs[i][2];
            new_dsb_blunted_ends[5] = dsb_locs[i][0];                             // strand1 location of the dsb that closes this cluster (causes the right edge)
            new_dsb_blunted_ends[6] = dsb_locs[i][1];                             // strand2 location of the dsb that closes this cluster (causes the right edge)
            dsb_blunted_ends[groupTID].push_back(new_dsb_blunted_ends);
        }
    }
}

// Populates the dsb_fragments_left and dsb_fragments_right arrays with all the dsb fragments in the genome after DNA fragmentation.
// Fragments are stored as vectors in the following format: {location of the fragment's end at the DSB-caused blunted end, location of the
// fragment's other end, chromosome index, strand1 (backbone1) location of the dsb that caused this blunted end,
// strand2 (backbone2) location of the dsb that caused this blunted end}. For left fragments, element
// 0 > element 1 ; for right fragments, element 1 > element 0 
void InduceSeq::load_dsb_fragments(int groupTID, int threadID) {
    // left fragments are those that extend from a dsb end towards lesser genome positions. 
    // they extend leftward when seen on a DNA diagram drawn with the usual convention where the top strand is 5' to 3'
    dsb_fragments_left[groupTID].clear(); 
    dsb_fragments_right[groupTID].clear();

    // State carried over from the previous dsb's right fragment, so the current dsb's left fragment
    // can be checked against it for overlap. previous_chrom_idx starts at -1 (never a valid chrom index)
    // so the very first dsb is never mistakenly treated as overlapping with a "previous" fragment.
    long previous_right_start = 0;
    long previous_right_end = 0;
    int previous_chrom_idx = -1;
    long previous_right_dsb_strand1 = 0;                                            // The right-causing dsb's strand1/strand2 locations for the previous blunted end,
    long previous_right_dsb_strand2 = 0;                                            // carried over in case the previous dsb's right fragment needs to be re-pushed after an overlap merge
    bool previous_right_appended = false;                                          // Whether the previous dsb's right fragment is currently the last element of dsb_fragments_right[groupTID]

    // Set by the right-fragment check below when dsb i's right fragment runs into dsb i+1's left
    // blunted end, so that dsb i+1's left fragment is skipped entirely on the next loop iteration.
    bool skip_left_fragment = false;

    size_t i = 0;
    while (i < dsb_blunted_ends[groupTID].size()) {
        std::vector<long> dsb_blunted_end = dsb_blunted_ends[groupTID][i];
        int chrom_idx = dsb_blunted_end[2];
        long left_dsb_strand1 = dsb_blunted_end[3];                                // strand1 location of the dsb that caused the left edge
        long left_dsb_strand2 = dsb_blunted_end[4];                                // strand2 location of the dsb that caused the left edge
        long right_dsb_strand1 = dsb_blunted_end[5];                               // strand1 location of the dsb that caused the right edge
        long right_dsb_strand2 = dsb_blunted_end[6];                               // strand2 location of the dsb that caused the right edge

        bool skip_this_left_fragment = skip_left_fragment;                         // Capture the flag set by the previous iteration before resetting it for this one
        skip_left_fragment = false;

        int fragment_length = get_random_fragment_length(threadID) - parameter.get_P5_adapter_length();
        if (fragment_length > 1 && !skip_this_left_fragment) {
            long left_start = dsb_blunted_end[0];
            long left_end = left_start - fragment_length + 1;
            long chrom_start = chrom_end_loc[chrom_idx] + 1;
            bool drop_current_left = false;
            if (left_end < chrom_start) {
                left_end = chrom_start;
            } else if (previous_right_appended && chrom_idx == previous_chrom_idx && left_end < previous_right_end) {
                // The current left fragment overlaps with the previous dsb's right fragment.
                long overlap = previous_right_end - left_end;                      // How far the two fragments overlap
                long blunted_end_gap = left_start - previous_right_start;          // Raw distance between the two dsb break points (not the fragment ends)
                if (overlap > parameter.get_maximum_overlap_fragment_generation() * (blunted_end_gap + 2*parameter.get_P5_adapter_length())) {
                    // Overlap is large; the section in between these two DSBs is taken to not have been fragmented. No reads produced from this fragment, since it contains 2 P5 adapters. 
                    // right fragment that was already pushed, and skip pushing the current left fragment.
                    dsb_fragments_right[groupTID].pop_back();
                    drop_current_left = true;
                } else {
                    // Overlap is small enough to salvage: both fragments are shortened to meet at the midpoint of their two ends.
                    long merged_end = static_cast<long>(0.5 * (left_end + previous_right_end));
                    left_end = merged_end;

                    // Re-push the previous right fragment with its new (shortened) end
                    dsb_fragments_right[groupTID].pop_back();
                    dsb_fragments_right[groupTID].push_back({previous_right_start, merged_end, previous_chrom_idx, previous_right_dsb_strand1, previous_right_dsb_strand2});
                }
            }

            if (!drop_current_left) {
                dsb_fragments_left[groupTID].push_back({left_start, left_end, chrom_idx, left_dsb_strand1, left_dsb_strand2});
            }
        }

        fragment_length = get_random_fragment_length(threadID) - parameter.get_P5_adapter_length();
        long right_start = dsb_blunted_end[1];
        long right_end = right_start + fragment_length - 1;
        bool right_appended = false;
        if (fragment_length > 1) {
            long chrom_end = chrom_end_loc[chrom_idx + 1];
            if (right_end > chrom_end) right_end = chrom_end;

            // Check whether this right fragment runs into the next dsb's left blunted end (plus the adapter length).
            // If so, drop this right fragment and remember to also skip the next dsb's left fragment.
            bool overlaps_next_dsb = false;
            if (i + 1 < dsb_blunted_ends[groupTID].size()) {
                const std::vector<long>& next_dsb_blunted_end = dsb_blunted_ends[groupTID][i+1];
                if (chrom_idx == next_dsb_blunted_end[2] && right_end > next_dsb_blunted_end[0] + parameter.get_P5_adapter_length()) {
                    overlaps_next_dsb = true;
                    skip_left_fragment = true;
                }
            }

            if (!overlaps_next_dsb) {
                dsb_fragments_right[groupTID].push_back({right_start, right_end, chrom_idx, right_dsb_strand1, right_dsb_strand2});
                right_appended = true;
            }
        }

        previous_right_start = right_start;
        previous_right_end = right_end;
        previous_chrom_idx = chrom_idx;
        previous_right_dsb_strand1 = right_dsb_strand1;
        previous_right_dsb_strand2 = right_dsb_strand2;
        previous_right_appended = right_appended;

        i++;
    }
}

// Filters dsb_fragments_left/right[groupTID], first by a flat probability_of_sequencing_multiplier
// chance (independent of fragment size) that a general fragment doesn't produce a read, then by fragment
// size: for each remaining fragment, looks up its keep probability in probability_of_sequencing_function
// (a length -> probability map, indexed by fragment length alone, not including the P5 adapter). If
// the fragment size is larger than the largest size in the map, the largest size's probability is used instead.
// The fragment is then removed with probability (1 - keep probability). 
void InduceSeq::filter_dsb_fragments(int groupTID, int threadID) {
    std::vector<std::vector<long>>& left_fragments = dsb_fragments_left[groupTID];
    std::vector<std::vector<long>>& right_fragments = dsb_fragments_right[groupTID];

    // Flat (size-independent) chance that a fragment doesn't survive sequencing, applied before the
    // size-dependent filtering below.
    double keep_probability_flat = parameter.get_probability_of_sequencing_multiplier();
    auto should_remove_flat = [threadID, keep_probability_flat](const std::vector<long>&) {
        return rng::rand_double(0.0, 1.0, threadID) >= keep_probability_flat;
    };
    left_fragments.erase(std::remove_if(left_fragments.begin(), left_fragments.end(), should_remove_flat), left_fragments.end());
    right_fragments.erase(std::remove_if(right_fragments.begin(), right_fragments.end(), should_remove_flat), right_fragments.end());

    int largest_tabulated_size = probability_of_sequencing_function.rbegin()->first;

    auto should_remove = [this, threadID, largest_tabulated_size](long size) {
        double keep_probability = (size > largest_tabulated_size)
            ? probability_of_sequencing_function.rbegin()->second
            : probability_of_sequencing_function.lower_bound(size)->second;
        return rng::rand_double(0.0, 1.0, threadID) >= keep_probability;
    };

    left_fragments.erase(
        std::remove_if(left_fragments.begin(), left_fragments.end(), [this, &should_remove](const std::vector<long>& frag) {
            return should_remove(frag[0] - frag[1] + 1);
        }),
        left_fragments.end()
    );

    right_fragments.erase(
        std::remove_if(right_fragments.begin(), right_fragments.end(), [this, &should_remove](const std::vector<long>& frag) {
            return should_remove(frag[1] - frag[0] + 1);
        }),
        right_fragments.end()
    );
}

// Filters dsb_fragments_left/right[groupTID] in place, removing fragments based on whether a single strand break is
// present on the fragment, which would cause the fragment (a denatured DNA strand) to not be sequenced. Surviving
// entries are left unchanged, in the same per-entry format as before filtering (see load_dsb_fragments).
void InduceSeq::filter_dsb_strands_ssd(int groupTID) {
    std::vector<std::vector<long>>& left_fragments = dsb_fragments_left[groupTID];
    std::vector<std::vector<long>>& right_fragments = dsb_fragments_right[groupTID];
    std::vector<long> strand_1_breaks = sdd_data.get_backbone1_break_loc(groupTID);
    std::vector<long> strand_2_breaks = sdd_data.get_backbone2_break_loc(groupTID);

    size_t ssd_i = 0;
    left_fragments.erase(
        std::remove_if(left_fragments.begin(), left_fragments.end(), [&](const std::vector<long>& dsb_frag) {
            long frag_end = dsb_frag[1];
            long causing_break = dsb_frag[3];                                      // backbone1 break that formed this dsb

            bool is_good = true;

            // checks for ssbs
            // Because of overhang fill-in during blunting, the edge of a fragment can be at a different location than the original dsb edge.
            // damages with positions in between these two positions would not be on the fragment, since that part of the fragment is remade.
            // That is why causing_break is used as a bound, instead of the end of the break.
            while (ssd_i < strand_1_breaks.size() && strand_1_breaks[ssd_i] < causing_break) {
                if (strand_1_breaks[ssd_i] >= frag_end) {
                    is_good = false;
                }
                ssd_i++;
            }

            return !is_good;
        }),
        left_fragments.end()
    );

    ssd_i = 0;
    right_fragments.erase(
        std::remove_if(right_fragments.begin(), right_fragments.end(), [&](const std::vector<long>& dsb_frag) {
            long frag_start = dsb_frag[0];
            long frag_end = dsb_frag[1];
            long causing_break = dsb_frag[4];                                      // strand2 location of the dsb that caused the right edge

            bool is_good = true;

            // Because of overhang fill-in during blunting, the edge of a fragment can be at a different location than the original dsb edge.
            // damages with positions in between these two positions would not be on the fragment, since that part of the fragment is remade. This loop filters them out.
            while (ssd_i < strand_2_breaks.size() && strand_2_breaks[ssd_i] <= causing_break) {
                ssd_i++;
            }
            // Breaks strictly after the causing break are genuinely
            while (ssd_i < strand_2_breaks.size() && strand_2_breaks[ssd_i] < frag_end) {
                if (strand_2_breaks[ssd_i] >= frag_start) {
                   is_good = false;
                }
                ssd_i++;
            }

            return !is_good;
        }),
        right_fragments.end()
    );
}


// Populates base_pair_damages_left/right[groupTID], one entry per surviving fragment in dsb_fragments_left/right[groupTID]
// (same index, in the same order), with the list of base-pair damage locations that are on that fragment, on the strand that
// produces an INDUCEseq read: basestrand1 for left fragments, basestrand2 for right fragments. The list can be empty for a fragment if
// the fragment contains no damage. Assumes dsb_fragments_left/right[groupTID] and the basestrand1/2 damage location
// vectors are all sorted ascending: bp_i is advanced monotonically across all fragments of a side (not reset per
// fragment), so this runs in a single linear pass over each damage-location vector rather than re-scanning it per fragment.
void InduceSeq::find_base_pair_damages(int groupTID) {
    base_pair_damages_left[groupTID].clear();
    base_pair_damages_right[groupTID].clear();

    size_t bp_i = 0;
    std::vector<long> bp_damages = sdd_data.get_basestrand1_damage_loc(groupTID);
    for (std::vector<long> dsb_strand : dsb_fragments_left[groupTID]) {
        std::vector<long> bp_damages_in_strand;
        while (bp_i < bp_damages.size() && bp_damages[bp_i] <= dsb_strand[0]) {
            if (bp_damages[bp_i] >= dsb_strand[1]) {
                bp_damages_in_strand.push_back(bp_damages[bp_i]);
            }
            bp_i++;
        }
        base_pair_damages_left[groupTID].push_back(bp_damages_in_strand);
    }

    bp_i = 0;
    bp_damages = sdd_data.get_basestrand2_damage_loc(groupTID);
    for (std::vector<long> dsb_strand : dsb_fragments_right[groupTID]) {
        std::vector<long> bp_damages_in_strand;
        while (bp_i < bp_damages.size() && bp_damages[bp_i] <= dsb_strand[1]) {
            if (bp_damages[bp_i] >= dsb_strand[0]) {
                bp_damages_in_strand.push_back(bp_damages[bp_i]);
            }
            bp_i++;
        }
        base_pair_damages_right[groupTID].push_back(bp_damages_in_strand);
    }
}

// Generates output files containing generated reads and other output, for a given groupTID.
// cell_number is used for file naming purposes. 
// Function is multithreaded, num_available_threads describes how many threads are available for it. 
// threadIDOffset is used to obtain global threadIDs for all worker threads; these IDS are threadIDOffset + localThreadID,
// where localThreadID goes from 0, 1, 2 to num_available_threads.  
// a global threadID is needed to generate reads from an ART object, since other calls of generate_simulation_output might be running in parallel
void InduceSeq::generate_simulation_output(int cell_number, int groupTID, int num_available_threads, int threadIDOffset) {

    // ofstream object of the output fastq file for read 1. Not opened/created at all if generate_reads is False.
    std::string output_file_ending;
    std::string output_fastq_R1_filename;
    std::ofstream fastq_R1_file;
    if (parameter.get_generate_reads()) {
        if (parameter.get_compress_output()) {
            output_file_ending = ".fastq.gz";
            output_fastq_R1_filename = (*parameter.get_output_directory())+"/"+(*parameter.get_output_fastq_filename_prefix())+"_"+std::to_string(cell_number)+"_R1" + output_file_ending;
            fastq_R1_file.open(output_fastq_R1_filename.c_str(),std::ios::binary);
        } else {
            output_file_ending = ".fastq";
            output_fastq_R1_filename = (*parameter.get_output_directory())+"/"+(*parameter.get_output_fastq_filename_prefix())+"_"+std::to_string(cell_number)+"_R1" + output_file_ending;
            fastq_R1_file.open(output_fastq_R1_filename.c_str());
        }
    }

    // ofstream object of the output file listing the DSBs that were sequenced, if requested. Never compressed.
    std::string output_sequenced_dsbs_filename;
    std::ofstream sequenced_dsbs_file;
    if (parameter.get_output_sequenced_dsbs()) {
        output_sequenced_dsbs_filename = (*parameter.get_output_directory())+"/"+(*parameter.get_output_fastq_filename_prefix())+"_"+std::to_string(cell_number)+"_sequenced_dsbs.csv";
        sequenced_dsbs_file.open(output_sequenced_dsbs_filename.c_str());
        sequenced_dsbs_file << "fragment start location,fragment end location,fragment start location relative to chromosome start,fragment end location relative to chromosome start,chromosome index (starting at 0 in order of chromosome sizes listed in sdd file),is forward read,read number\n";
    }

    ART read1;                                                                                      // Creating an ART class object and setting the insertion and deletion probability vectors for that read object
    if (parameter.get_generate_reads()) {
        read1.set_read_error_rates(parameter.get_insertion_error_rate_read1(), parameter.get_deletion_error_rate_read1());
        read1.set_read_error_probability(parameter.get_read_length(), parameter.get_insertion_error_rate_read1(), read1.insertion_probability_vec, parameter.get_max_errors_in_read());
        read1.set_read_error_probability(parameter.get_read_length(), parameter.get_deletion_error_rate_read1(), read1.deletion_probability_vec, parameter.get_max_errors_in_read());
        read1.resize_vectors(parameter.get_number_of_threads());                                     // Sized to the full global thread count since threadID below spans that range, not just this call's num_available_threads
    }

    std::string chromSegSeq;                                                                        // Temporary variable to hold each chromosome segment sequence from the fasta file one at a time
    std::string chromSegSeq_ID;                                                                     // Temporary variable to hold IDs of each hromosome segment sequence


    const int batchSize{2000};                                                                      // Define a batch size for writing reads to the output file. These much data will be stored in cache before writing it on the file

    int batchSize_thread = std::round(batchSize/num_available_threads);                                     // Devide the total cache size for the buffer equally for all the threads

    std::vector<std::vector<std::string>> batch_buffer(num_available_threads);
    std::vector<std::vector<std::string>> dsb_batch_buffer(num_available_threads);                   // Per-thread buffer for the sequenced-DSBs CSV lines, flushed the same way as batch_buffer

    #pragma omp parallel for num_threads(num_available_threads)
    for (size_t i=0; i<dsb_fragments_left[groupTID].size() + dsb_fragments_right[groupTID].size(); i++) {
        int localTID = omp_get_thread_num();                                                       // ID local to this call's team (0..num_available_threads-1); safe to index batch_buffer, which is private to this call
        int threadID = threadIDOffset + localTID;                                                  // Globally-unique ID (0..nThreads_User-1) required by rng::local_mt and ART's per-thread buffers, which are shared across all concurrently-running groups
        std::vector<long> dsb_strand;
        std::vector<long> bp_damages;
        bool is_left;
        if (i < dsb_fragments_left[groupTID].size()) {
            dsb_strand = dsb_fragments_left[groupTID][i];
            bp_damages = base_pair_damages_left[groupTID][i];
            is_left = true;
        } else {
            dsb_strand = dsb_fragments_right[groupTID][i - dsb_fragments_left[groupTID].size()];
            bp_damages = base_pair_damages_right[groupTID][i - dsb_fragments_left[groupTID].size()];
            is_left = false;
        }

        if (parameter.get_output_sequenced_dsbs()) {
            long chrom_start = chrom_end_loc[dsb_strand[2]];
            std::string dsb_data = std::to_string(dsb_strand[0])+","+std::to_string(dsb_strand[1])+","+std::to_string(dsb_strand[0]-chrom_start)+","+std::to_string(dsb_strand[1]-chrom_start)+","+std::to_string(dsb_strand[2])+","+std::to_string(!is_left)+","+std::to_string(i)+"\n";
            dsb_batch_buffer[localTID].push_back(dsb_data);                         // Add the DSB data to the buffer vector of the respective thread

            if (dsb_batch_buffer[localTID].size() >= static_cast<size_t>(batchSize_thread)) {// Check if the batch buffer is full, and write it to the file if needed.
                #pragma omp critical(section2)
                {
                    writeBatchToFile(dsb_batch_buffer[localTID], sequenced_dsbs_file, false);
                }
            }
        }

        if (parameter.get_generate_reads()) {
            std::string dna_seq;
            get_dna_sequence(dna_seq, bp_damages, dsb_strand, is_left);
            read1.generate_read_with_indel_from_frag(dna_seq, P7_adapter_sequence, threadID);               // Make a read with random indel errors
            std::vector<short> read1_quality_score_vec;                                 // Vector to hold the quality scores for read 1
            read1.get_read_quality(read1_quality_score_vec, 1, threadID);               // Get the read quality scores for the read positions
            read1.add_baseCall_error(read1_quality_score_vec, threadID);                // Add base call errors to the read based on the quality scores

            std::string chromID = chrom_headers[dsb_strand[2]];
            std::string read_data = "@"+chromID+"_read"+std::to_string(i)+"\n";                           // @readID
            read_data += (*read1.get_final_read_sequence(threadID))+ "\n+\n";           // read sequence and +
            for(size_t k=0; k<(*read1.get_final_read_sequence(threadID)).size(); k++){  // read quality scores; insert only as many quality values as with the length of sequence
                read_data += static_cast<char>(read1_quality_score_vec[k]+32);          // +33 to get the phred score
            }
            read_data += "\n";

            batch_buffer[localTID].push_back(read_data);                                // Add the read data to the buffer vector of the respective thread

            if (batch_buffer[localTID].size() >= static_cast<size_t>(batchSize_thread)) {// Check if the batch buffer is full, and write it to the file if needed.
                if (parameter.get_compress_output()) {
                    // compression can be done in paralell, since there is no shared memory.
                    std::string compressed_batch = getCompressedBatch(batch_buffer[localTID]);
                    // writing can only be done by one thread at a time.
                    // writeBatchToFile with compression = true is not used so that compression and writing can be done in separate blocks
                    #pragma omp critical(section1)
                    {
                        fastq_R1_file.write(compressed_batch.c_str(), compressed_batch.size());
                    }
                } else {
                    #pragma omp critical(section1)
                    {
                        writeBatchToFile(batch_buffer[localTID], fastq_R1_file, false);
                    }
                }
            }
        }
    }
    if (parameter.get_generate_reads()) {
        for (size_t l=0;l<batch_buffer.size();l++){
            writeBatchToFile(batch_buffer[l], fastq_R1_file, parameter.get_compress_output());      // If there are unwritten data in batch buffer, write that too when the loop ends
        }
        fastq_R1_file.close();

        // Record this cell's data for the summary report. generate_simulation_output can run concurrently for different
        // cells (each handled by a different thread group), so the shared report_* vectors are guarded with a critical section.
        long reads_generated_this_cell = static_cast<long>(dsb_fragments_left[groupTID].size() + dsb_fragments_right[groupTID].size());
        #pragma omp critical(summary_report_data)
        {
            report_cells_sequenced.push_back("Damaged_cell_"+std::to_string(cell_number));
            report_fastq_output.push_back((*parameter.get_output_fastq_filename_prefix())+"_"+std::to_string(cell_number));
            report_readsGenerated_perCell.push_back(static_cast<int>(reads_generated_this_cell));
        }
    }
    if (parameter.get_output_sequenced_dsbs()) {
        for (size_t l=0;l<dsb_batch_buffer.size();l++){
            writeBatchToFile(dsb_batch_buffer[l], sequenced_dsbs_file, false);      // If there is unwritten data in the dsb batch buffer, write that too when the loop ends
        }
        sequenced_dsbs_file.close();
    }
}

// Sets dna_seq to be the DNA sequence of a dsb_strand, in the order that it will generate a read. bp_damages is an array specifying the 
// locations of damages on the strand; these bases are replaced by N. is_left specifies whether the strand came from dsb_fragments_left (generates a backward read) or dsb_fragments_right (forward read)
void InduceSeq::get_dna_sequence(std::string& dna_seq, std::vector<long>& bp_damages, std::vector<long>& dsb_strand, bool is_left) {
    int chrom_idx = dsb_strand[2];
    long start_char_i = dsb_strand[0] + cum_chrom_header_sizes[chrom_idx] - 1;
    long end_char_i = dsb_strand[1] + cum_chrom_header_sizes[chrom_idx] - 1;

    if (is_left) {
        dna_seq = std::string(genome_fasta + end_char_i, genome_fasta + start_char_i + 1);
        for (long bp_damage : bp_damages) {
            int dam_i = bp_damage - dsb_strand[1];
            dna_seq[dam_i] = 'N';
        }
        std::reverse(dna_seq.begin(), dna_seq.end());
        for (char& base : dna_seq) {
            if      (base == 'A') base = 'T';
            else if (base == 'T') base = 'A';
            else if (base == 'C') base = 'G';
            else if (base == 'G') base = 'C';
        }
    } else {
        dna_seq = std::string(genome_fasta + start_char_i, genome_fasta + end_char_i + 1);
        for (long bp_damage : bp_damages) {
            int dam_i = bp_damage - dsb_strand[0];
            dna_seq[dam_i] = 'N';
        }
    }
}

// returns a random fragment length according to the fragment size distribution parameter
// thread ID is a global thread index (from 0 to total number of threads used in program - 1)
int InduceSeq::get_random_fragment_length(int threadID) {
    float rand_val = rng::rand_float(0.0f, 1.0f, threadID);
    int fragment_length = fragment_size_distribution.upper_bound(rand_val)->second;
    return fragment_length;
}
