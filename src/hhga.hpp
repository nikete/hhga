#include <map>
#include <vector>
#include <string>
#include <getopt.h>
#include <stdlib.h>
#include <iostream>
#include "bamtools/api/BamMultiReader.h"
#include "bamtools/api/BamWriter.h"
#include "bamtools/api/SamReadGroup.h"
#include "Fasta.h"
#include "Variant.h"
#include <cmath>

namespace hhga {

using namespace BamTools;
using namespace std;

#define PHRED_MAX 1000

short qualityChar2ShortInt(char c);
long double qualityChar2LongDouble(char c);
long double lnqualityChar2ShortInt(char c);
char qualityInt2Char(short i);
long double ln2log10(long double prob);
long double log102ln(long double prob);
long double phred2ln(int qual);
long double ln2phred(long double prob);
long double phred2float(int qual);
long double float2phred(long double prob);
void parse_region(const string& region,
                  string& startSeq,
                  int& startPos,
                  int& stopPos);
void set_region(BamMultiReader& reader, const string& region_str);
void set_region(vcflib::VariantCallFile& vcffile, const string& region_str);

class HHGA {
public:
    string chrom_name;
    int start_pos;
    int stop_pos;

    typedef vcflib::VariantAllele allele_type;
    set<allele_type> alleles;
    vector<BamAlignment> alignments;
    map<BamAlignment*, vector<allele_type> > alignment_alleles;
    HHGA(const string& region_string, BamMultiReader& bam_reader, FastaReference& fasta_ref);
    const string str(void);
};

}
