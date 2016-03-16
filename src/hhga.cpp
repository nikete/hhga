#include "hhga.hpp"

namespace hhga {

short qualityChar2ShortInt(char c) {
    return static_cast<short>(c) - 33;
}

long double qualityChar2LongDouble(char c) {
    return static_cast<long double>(c) - 33;
}

long double lnqualityChar2ShortInt(char c) {
    return log(static_cast<short>(c) - 33);
}

char qualityInt2Char(short i) {
    return static_cast<char>(i + 33);
}

long double ln2log10(long double prob) {
    return M_LOG10E * prob;
}

long double log102ln(long double prob) {
    return M_LN10 * prob;
}

long double phred2ln(int qual) {
    return M_LN10 * qual * -.1;
}

long double ln2phred(long double prob) {
    return -10 * M_LOG10E * prob;
}

long double phred2float(int qual) {
    return pow(10, qual * -.1);
}

long double float2phred(long double prob) {
    if (prob == 1)
        return PHRED_MAX;  // guards against "-0"
    long double p = -10 * (long double) log10(prob);
    if (p < 0 || p > PHRED_MAX) // int overflow guard
        return PHRED_MAX;
    else
        return p;
}

void parse_region(
    const string& region,
    string& startSeq,
    int& startPos,
    int& stopPos) {

    size_t foundFirstColon = region.find(":");

    // we only have a single string, use the whole sequence as the target
    if (foundFirstColon == string::npos) {
        startSeq = region;
        startPos = 0;
        stopPos = -1;
    } else {
        startSeq = region.substr(0, foundFirstColon);
        string sep = "..";
        size_t foundRangeSep = region.find(sep, foundFirstColon);
        if (foundRangeSep == string::npos) {
            sep = "-";
            foundRangeSep = region.find("-", foundFirstColon);
        }
        if (foundRangeSep == string::npos) {
            startPos = atoi(region.substr(foundFirstColon + 1).c_str());
            // differ from bamtools in this regard, in that we process only
            // the specified position if a range isn't given
            stopPos = startPos + 1;
        } else {
            startPos = atoi(region.substr(foundFirstColon + 1, foundRangeSep - foundFirstColon).c_str());
            // if we have range sep specified, but no second number, read to the end of sequence
            if (foundRangeSep + sep.size() != region.size()) {
                stopPos = atoi(region.substr(foundRangeSep + sep.size()).c_str()); // end-exclusive, bed-format
            } else {
                //stopPos = reference.sequenceLength(startSeq);
                stopPos = -1;
            }
        }
    }
}

void set_region(BamMultiReader& reader,
                const string& region_str) {

    // parse the region string
    if (!region_str.empty()) {

        map<string, int> refLength;
        map<string, int> refID;

        int id = 0;
        RefVector references = reader.GetReferenceData();
        for (RefVector::iterator r = references.begin(); r != references.end(); ++r) {
            refLength[r->RefName] = r->RefLength;
            refID[r->RefName] = id++;
        }

        // parse the region string
        string startSeq;
        int startPos;
        int stopPos;

        parse_region(region_str, startSeq, startPos, stopPos);

        if (stopPos == -1) {
            stopPos = refLength[startSeq];
        }

        int startSeqRefID = refID[startSeq];

        if (!reader.LocateIndexes()) {
            cerr << "[hhga] could not load BAM index" << endl;
            exit(1);
        } else {
            reader.SetRegion(startSeqRefID, startPos, startSeqRefID, stopPos);
        }
    }
}

void set_region(vcflib::VariantCallFile& vcffile, const string& region_str) {
    if (vcffile.is_open()) {
        vcffile.setRegion(region_str);
    }
}


HHGA::HHGA(const string& region_str, BamMultiReader& bam_reader, FastaReference& fasta_ref) {

    // store the names of all the reference sequences in the BAM file
    map<int, string> referenceIDToName;
    vector<RefData> referenceSequences = bam_reader.GetReferenceData();
    int i = 0;
    for (RefVector::iterator r = referenceSequences.begin(); r != referenceSequences.end(); ++r) {
        referenceIDToName[i] = r->RefName;
        ++i;
    }

    // parse the region string
    parse_region(region_str, chrom_name, start_pos, stop_pos);

    // set up our readers
    set_region(bam_reader, region_str);

    long int lowestReferenceBase = 0;
    long unsigned int referenceBases = 0;
    unsigned int currentRefSeqID = 0;

    // get the alignments at the locus
    BamAlignment aln;
    while (bam_reader.GetNextAlignment(aln)) {
        if (aln.IsMapped()) {
            alignments.push_back(aln);
        }
    }

    for (auto& aln : alignments) {
        // iterate through the alignment
        // converting it into a series of alleles

        // record the qualities
        //for (string::iterator c = aln.Qualities.begin(); c != al.Qualities.end(); ++c) {
        //    ++qual_hist[qualityChar2ShortInt(*c)];
        //

        long unsigned int endpos = aln.GetEndPosition();
        string refseq = fasta_ref.getSubSequence(referenceIDToName[aln.RefID],
                                                 aln.Position,
                                                 aln.GetEndPosition() - (aln.Position - 1));
        const string& readseq = aln.QueryBases;
        int rel_pos = aln.Position - this->start_pos;

        int rp = 0; int sp = 0;

        vector<allele_type>& aln_alleles = alignment_alleles[&aln];

        vector<CigarOp>::const_iterator cigarIter = aln.CigarData.begin();
        vector<CigarOp>::const_iterator cigarEnd  = aln.CigarData.end();
        for ( ; cigarIter != cigarEnd; ++cigarIter ) {
            unsigned int len = cigarIter->Length;
            char t = cigarIter->Type;        
            switch (t) {
            case 'I':
                aln_alleles.push_back(allele_type("", readseq.substr(sp, len), rp + rel_pos));
                sp += len;
                break;
            case 'D':
                aln_alleles.push_back(allele_type(refseq.substr(rp, len), "", rp + rel_pos));
                rp += len;
                break;
            case 'M':
                {
                    for (int i = 0; i < len; ++i) {
                        allele_type a = allele_type(refseq.substr(rp + i, 1),
                                                    readseq.substr(sp + i, 1),
                                                    rp + i + rel_pos);
                        aln_alleles.push_back(a);
                    }
                }
                rp += len;
                sp += len;
                break;
            case 'S':
                rp += len;
                sp += len;
                break;
            default:
                cerr << "do not recognize cigar element " << t <<":"<< len << endl;
                break;
            }

        }
    }
}

const string HHGA::str(void) {
    //return std::to_string(alleles.size());
    stringstream out;
    for (auto& aln : alignments) {
        out << "<<<<<----------------------------------------" << endl;
        out << aln.Name << endl;
        out << aln.QueryBases << endl;
        for (auto& allele : alignment_alleles[&aln]) {
           out << allele << ", ";
        }
        out << endl;
        out << "----------------------------------------->>>>" << endl;
    }
    return out.str();
}

}
