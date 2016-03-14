#include "hhga.hpp"

using namespace std;
using namespace hhga;

void printUsage(int argc, char** argv) {

    cerr << "usage: " << argv[0] << " [-b FILE]" << endl
         << endl
         << "options:" << endl
         << "    -h, --help            this dialog" << endl
         << "    -f, --fasta-reference FILE  the reference sequence" << endl
         << "    -b, --bam FILE        use this BAM as input (multiple allowed)" << endl
         << "    -r, --region REGION   limit output to those in this region (chr:start-end)" << endl
         << endl
         << "Generates reports on the rate of putative mutations or errors in the input alignment data." << endl
         << "Alignments are read from the specified files, or stdin if none are specified" << endl
         << endl
         << "authors: Erik Garrison <erik.garrison@gmail.com> and Nicolás Della Penna <nikete@gmail.com>" << endl;

}

int main(int argc, char** argv) {

    vector<string> inputFilenames;
    string vcf_file_name;

    string regionStr;

    string fastaFile;

    // parse command-line options
    int c;

    while (true) {

        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"bam",  required_argument, 0, 'b'},
            {"vcf", required_argument, 0, 'v'},
            {"region", required_argument, 0, 'r'},
            {"fasta-reference", required_argument, 0, 'f'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hb:r:f:v:",
                         long_options, &option_index);

        if (c == -1)
            break;
 
        switch (c) {

        case '?':
            printUsage(argc, argv);
            return 0;
            break;

        case 'h':
            printUsage(argc, argv);
            return 0;
            break;

        case 'b':
            inputFilenames.push_back(optarg);
            break;

        case 'v':
            vcf_file_name = optarg;
            break;

        case 'r':
            regionStr = optarg;
            break;

        case 'f':
            fastaFile = optarg;
            break;

        default:
            return 1;
            break;
        }
    }

    if (fastaFile.empty()) {
        cerr << "no FASTA reference specified" << endl;
        return 1;
    }

    if (inputFilenames.empty()) {
        cerr << "no input files specified" << endl;
        return 1;
    }

    BamMultiReader reader;
    if (!reader.Open(inputFilenames)) {
        cerr << "could not open input BAM files" << endl;
        return 1;
    }

    vcflib::VariantCallFile vcf_file;
    if (!vcf_file_name.empty()) {
        vcf_file.open(vcf_file_name);
        if (!vcf_file.is_open()) {
            cerr << "could not open " << vcf_file_name << endl;
            return 1;
        }
    }

    setRegion(regionStr, reader, vcf_file);

    // store the names of all the reference sequences in the BAM file
    map<int, string> referenceIDToName;
    vector<RefData> referenceSequences = reader.GetReferenceData();
    int i = 0;
    for (RefVector::iterator r = referenceSequences.begin(); r != referenceSequences.end(); ++r) {
        referenceIDToName[i] = r->RefName;
        ++i;
    }

    FastaReference fr;
    fr.open(fastaFile);

    long unsigned int alignedBases = 0;
    long unsigned int mismatchCount = 0;
    long unsigned int gapCount = 0;
    map<int, long unsigned int> mismatches;
    map<int, long unsigned int> gaps;

    long int lowestReferenceBase = 0;
    long unsigned int referenceBases = 0;
    unsigned int currentRefSeqID = 0;

    map<short, uint64_t> qual_hist;

    BamAlignment al;
    while (reader.GetNextAlignment(al)) {
        if (al.IsMapped()) {

            // record the qualities
            for (string::iterator c = al.Qualities.begin(); c != al.Qualities.end(); ++c) {
                ++qual_hist[qualityChar2ShortInt(*c)];
            }
            
            long unsigned int endpos = al.GetEndPosition();
            // this happens when we switch reference sequences
            if (currentRefSeqID != al.RefID) {
                //cout << "al.GetEndPosition() = " << endpos << "  lowestReferenceBase = " << lowestReferenceBase << "  reset" << endl;
                currentRefSeqID = al.RefID;
                referenceBases += lowestReferenceBase;
                lowestReferenceBase = 0;
            } else if (endpos > lowestReferenceBase) {
                //cout << "al.GetEndPosition() = " << endpos << "  lowestReferenceBase = " << lowestReferenceBase << "  adding " << endpos - lowestReferenceBase << endl;
                referenceBases += endpos - lowestReferenceBase;
                lowestReferenceBase = endpos;
            }

            //cout << al.Position << endl;
            //cout << al.AlignedBases << endl;
            string refsequence = fr.getSubSequence(referenceIDToName[al.RefID], al.Position, al.GetEndPosition() - (al.Position - 1));
            //cout << refsequence << endl;

            alignedBases += al.QueryBases.size();

            int rp = 0; int sp = 0;
            vector<CigarOp>::const_iterator cigarIter = al.CigarData.begin();
            vector<CigarOp>::const_iterator cigarEnd  = al.CigarData.end();
            for ( ; cigarIter != cigarEnd; ++cigarIter ) {
                unsigned int l = cigarIter->Length;
                char t = cigarIter->Type;

                if (t == 'M') { // match or mismatch

                    int firstMismatch = -1;

                    for (int i=0; i<l; i++) {

                        // extract aligned base
                        char b = al.QueryBases.at(rp);

                        // get reference allele
                        char sb = refsequence.at(sp);

                        // record mismatch if we have a mismatch here
                        if (firstMismatch >= 0) {
                            if (b == sb) {
                                // mismatch termination
                                // register multi-base mismatch
                                int length = rp - firstMismatch;
                                //string qualstr = rQual.substr(rp - length, length);
                                ++mismatches[length];
                                mismatchCount += length;
                                firstMismatch = -1;
                            } else {
                                // mismatch extension
                            }
                        } else {
                            if (b != sb) {
                                // mismatch state
                                firstMismatch = rp;
                            } else {
                                // match state
                            }
                        }

                        // update positions
                        ++sp;
                        ++rp;
                    }

                } else if (t == 'D') {
                    ++gaps[-l];
                    ++gapCount;
                    sp += l;
                } else if (t == 'I') {
                    ++gaps[l];
                    ++gapCount;
                    rp += l;
                } else if (t == 'S') {
                    rp += l;
                } else if (t == 'H') {
                } else if (t == 'N') {
                    sp += l;
                    rp += l;
                }
            }
        }
    }

    reader.Close();

    cout << "reference bases:\t" << referenceBases << endl;
    cout << "total aligned bases:\t" << alignedBases << endl;
    cout << "mean alignment depth:\t" << (long double) alignedBases / (long double) referenceBases << endl;
    cout << "total mismatched bases:\t" << mismatchCount << endl;
    cout << "total gap bases:\t" << gapCount << endl;
    cout << "mismatch rate per aligned bp:\t" << (long double) mismatchCount / (long double) alignedBases << endl;
    cout << "gap rate per aligned bp:\t" << (long double) gapCount / (long double) alignedBases << endl;
    cout << "mismatch + gap per aligned bp:\t" << (long double) ( gapCount + mismatchCount ) / (long double) alignedBases << endl;
    cout << endl;

    cout << "mismatch length distribution" << endl;
    cout << "length\tcount\trate (per aligned base)" << endl;
    for (map<int, long unsigned int>::iterator p = mismatches.begin(); p != mismatches.end(); ++p) {
        cout << p->first << "\t" << p->second << "\t" << (long double) p->second / (long double) alignedBases << endl;
    }
    cout << endl;

    cout << "gap length distribution:" << endl;
    cout << "length\tcount\trate (per aligned base)" << endl;
    for (map<int, long unsigned int>::iterator p = gaps.begin(); p != gaps.end(); ++p) {
        cout << p->first << "\t" << p->second << "\t" << (long double) p->second / (long double) alignedBases << endl;
    }

    cout << endl;
    cout << "quality histogram of alignmed reads" << endl;
    cout << "#qual\tcount" << endl;
    for (map<short, uint64_t>::const_iterator q = qual_hist.begin(); q != qual_hist.end(); ++q) {
        cout << q->first << "\t" << q->second << endl;
    }
    //cout << endl;

    return 0;

}
