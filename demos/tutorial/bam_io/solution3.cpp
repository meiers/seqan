#include <seqan/bam_io.h>

int main()
{
    // Open input file.
    seqan::BamFileIn bamFileIn;
    if (!open(bamFileIn, "example.sam"))
    {
        std::cerr << "ERROR: Could not open example.sam!" << std::endl;
        return 1;
    }

    try
    {
        // Read header.
        seqan::BamHeader header;
        readRecord(header, bamFileIn);

        unsigned tagIdx = 0;
        unsigned numXXtags = 0;

        // Rear records.
        seqan::BamAlignmentRecord record;
        while (!atEnd(bamFileIn))
        {
            readRecord(record, bamFileIn);
            seqan::BamTagsDict tagsDict(record.tags);

            if (findTagKey(tagIdx, tagsDict, "XX"))
                numXXtags += 1;
        }
    }
    catch (seqan::IOError const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "Number of records with the XX tag: " << numXXtags << "\n";

    return 0;
}
