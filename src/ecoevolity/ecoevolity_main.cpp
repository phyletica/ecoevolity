#include <ncl/nxsmultiformat.h>

#include "util.hpp"

int main(int argc, char *argv[]) {
    std::cout << "Hello World!" << std::endl;
    MultiFormatReader nexusReader(-1, NxsReader::WARNINGS_TO_STDERR);
    nexusReader.ReadFilepath(argv[1], MultiFormatReader::NEXUS_FORMAT);
    
    int numTaxaBlocks = nexusReader.GetNumTaxaBlocks();
    std::cout << numTaxaBlocks << " TAXA block(s) read.\n";
    for (int i = 0; i < numTaxaBlocks; ++i) {
        NxsTaxaBlock * taxaBlock = nexusReader.GetTaxaBlock(i);
        std::string taxaBlockTitle = taxaBlock->GetTitle();
        std::cout << "Taxa block index " << i << " has the Title \"" << taxaBlockTitle << "\"\n";
    
        const unsigned nCharBlocks = nexusReader.GetNumCharactersBlocks(taxaBlock);
        std::cout  <<  nCharBlocks << " CHARACTERS/DATA block(s) refer to this TAXA block\n";
        for (unsigned j = 0; j < nCharBlocks; ++j) {
            const NxsCharactersBlock * charBlock = nexusReader.GetCharactersBlock(taxaBlock, j);
            std::string charBlockTitle = taxaBlock->GetTitle();
            std::cout << "Taxa block index " << j << " has the Title \"" << charBlockTitle << "\"\n";
        }
    }
    return 0;
}
