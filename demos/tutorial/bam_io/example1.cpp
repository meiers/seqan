#include <seqan/bam_io.h>

using namespace seqan;

int main()
{
    // Open input file, BamFileIn can read SAM and BAM files.
    BamFileIn bamFileIn("example.sam");

    // Open output file, BamFileOut accepts also an ostream and a format tag.
    BamFileOut bamFileOut(std::cout, Sam());

    // Copy header.
    BamHeader header;
    readRecord(header, bamFileIn);
    writeRecord(bamFileOut, header);

    // Copy records.
    BamAlignmentRecord record;
    while (!atEnd(bamFileIn))
    {
        readRecord(record, bamFileIn);
        writeRecord(bamFileOut, record);
    }

    return 0;
}
