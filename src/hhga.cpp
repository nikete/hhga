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

void setRegion(
    const string& regionStr,
    BamMultiReader& reader,
    vcflib::VariantCallFile& vcffile) {

    // parse the region string
    if (!regionStr.empty()) {

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

        parse_region(regionStr, startSeq, startPos, stopPos);

        if (stopPos == -1) {
            stopPos = refLength[startSeq];
        }

        int startSeqRefID = refID[startSeq];

        if (!reader.LocateIndexes()) {
            cerr << "region specified, but could not open load BAM index" << endl;
            exit(1);
        } else {
            reader.SetRegion(startSeqRefID, startPos, startSeqRefID, stopPos);
        }

        if (vcffile.is_open()) {
            vcffile.setRegion(regionStr);
        }
    }
}

}
