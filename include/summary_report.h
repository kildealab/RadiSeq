#ifndef SUMMARY_REPORT_H
#define SUMMARY_REPORT_H

#include <string>
#include <vector>

// Having these variables defined in a header file lets us use them in all cpp files as a global variable. Much needed to generate the summary report
inline double report_cpu_time_used;                                                                     // Variable to hold the duration of the run time of this program
inline std::string report_parameterFileName;                                                            // Parameter file name for the report
inline long report_ref_seq_length;                                                                      // Reference sequence length for the report
inline std::vector<std::string> report_cells_sequenced;                                                 // Names of cells that were actually sequenced
inline std::vector<std::string> report_fastq_output;                                                    // Output file names without "_R1.fastq.gz" extension
inline std::vector<int> report_readsGenerated_perCell;                                                  // Number of reads generated from each cell in bulk-sequencing

#endif