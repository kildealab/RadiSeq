#ifndef FASTAFILE_HANDLER_H
#define FASTAFILE_HANDLER_H

#include <string>

#include "sddfile_handler.h"

//long buildUndamagedGenomeTemplate(const std::string&, int, int, const std::string*);
long buildUndamagedGenomeTemplate_MM(char*, std::size_t, int, int, const std::string*);
void getReverseComplementarySeq(std::string&, std::string&);
//int buildDamagedCellGenome(NGSsdd&, const std::string&, const std::string&);
int buildDamagedCellGenome_from_MM(NGSsdd&, const std::string&, const std::string&, char*, size_t);


#endif