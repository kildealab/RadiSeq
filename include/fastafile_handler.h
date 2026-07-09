#ifndef FASTAFILE_HANDLER_H
#define FASTAFILE_HANDLER_H

#include <string>

#include "sddfile_handler.h"

//long buildUndamagedGenomeTemplate(const std::string&, int, int, const std::string*);
long buildUndamagedGenomeTemplate_MM(char*, std::size_t, int, int, const std::string*, std::vector<double>&, double*, int);
int buildUndamagedGenomeTemplate_ForwardOnly_MM(char*, std::size_t, int, int, const std::string*, std::vector<int>&, const std::vector<long>&);
double getReverseComplementarySeq(const std::string&, std::string&, int GC_binSize=0);                               // GC_binSize is optional
int calculateCumChromHeaderSizes(std::vector<int>&, std::vector<std::string>&, char*, const std::vector<long>&);
//int buildDamagedCellGenome(NGSsdd&, const std::string&, const std::string&);
std::vector<double> buildDamagedCellGenome_from_MM(NGSsdd&, const std::string&, const std::string&, char*, size_t, long, int);


#endif